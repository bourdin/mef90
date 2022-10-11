#include "../MEF90/mef90.inc"
Program HeatXfer
#include "petsc/finclude/petsc.h"
   Use m_MEF90
   Use m_MEF90_HeatXfer
   Use m_MEF90_HeatXferCtx
   Use m_MEF90_HeatXferDefault
   Use petsc
   Implicit NONE   

   PetscErrorCode                                     :: ierr
   Type(MEF90Ctx_Type),target                         :: MEF90Ctx

   Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions
   Type(MEF90HeatXferCtx_Type)                        :: MEF90HeatXferCtx
   Type(MEF90HeatXferGlobalOptions_Type),Pointer      :: MEF90HeatXferGlobalOptions
                                                         
   Type(tDM)                                          :: dm,temperatureDM
   Type(tIS)                                          :: setIS
   PetscInt,Dimension(:),Pointer                      :: setID
   PetscInt                                           :: set
   PetscReal,Dimension(:),Pointer                     :: time,energy,cellWork,faceWork

   PetscBool                                          :: flg
   Character(len=MEF90MXSTRLEN)                       :: IOBuffer
   Type(tPetscViewer)                                 :: logViewer
   
   Type(tSNES)                                        :: temperatureSNES
   Type(tTS)                                          :: temperatureTS
   !Type(tTSAdapt)                                     :: temperatureTSAdapt
   Type(tVec)                                         :: temperature,temperatureResidual

   PetscReal                                          :: temperatureTSInitialStep
   !PetscInt                                           :: tsTempmaxIter
   !PetscReal                                          :: t
   
   PetscInt                                           :: step
   PetscInt                                           :: dim
      
   !!! Initialize MEF90
   PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
   PetscCallA(MEF90Initialize(ierr))

   !!! Get all MEF90-wide options
   PetscCallA(MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90DefaultGlobalOptions,ierr))
   PetscCallA(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr))

   PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm,MEF90Ctx%geometryFile,PETSC_NULL_CHARACTER,PETSC_TRUE,dm,ierr))
   PetscCallA(DMPlexDistributeSetDefault(dm,PETSC_FALSE,ierr))
   PetscCallA(DMSetFromOptions(dm,ierr))
   PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-heatXfer_dm_view",ierr))

   Inquire(file=MEF90Ctx%resultFile,exist=flg)
   If (flg) Then
      ! we assume that the output file exists and is formatted
      PetscCallA(MEF90CtxOpenEXO(MEF90Ctx,MEF90Ctx%resultViewer,FILE_MODE_APPEND,ierr))
   Else
      ! we need to create the output file
      EXOFormat: block
         PetscInt                                            :: numNodalVar = 1, numCellVar = 1, numFaceVar = 2, numGVar = 0
         Character(len=MEF90MXSTRLEN),Dimension(:),Pointer   :: faceVarName, nodalVarName, cellVarName, gVarName

         Allocate(nodalVarName(numNodalVar))
         Allocate(cellVarName(numCellVar))
         Allocate(faceVarName(numFaceVar))
         Allocate(gVarName(numGVar))
         nodalVarName = ["Temperature        "]
         cellVarName  = ["Flux               "]
         faceVarName  = ["boundaryFlux       ","externalTemperature"]
         PetscCallA(MEF90CtxOpenEXO(MEF90Ctx,MEF90Ctx%resultViewer,FILE_MODE_WRITE,ierr))
         PetscCallA(MEF90EXODMView(dm,MEF90Ctx%resultViewer,MEF90GlobalOptions%elementOrder,ierr))
         PetscCallA(MEF90EXOFormat(MEF90Ctx%resultViewer,gVarName,cellVarName,nodalVarName,faceVarName,time,ierr))
         DeAllocate(nodalVarName)
         DeAllocate(cellVarName)
         DeAllocate(faceVarName)
         DeAllocate(gVarName)
      End block EXOFormat
   End If
   distribute: Block 
       Type(tDM),target                    :: dmDist
       PetscInt                            :: ovlp = 0
       Type(tPetscSF)                      :: naturalPointSF

       If (MEF90Ctx%NumProcs > 1) Then
           PetscCallA(DMPlexDistribute(dm,ovlp,naturalPointSF,dmDist,ierr))
           PetscCallA(DMPlexSetMigrationSF(dmDist,naturalPointSF, ierr))
           PetscCallA(PetscSFDestroy(naturalPointSF,ierr))
           PetscCallA(DMDestroy(dm,ierr))
           dm = dmDist
       End If
   End Block distribute
   PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-heatXfer_dm_view",ierr))

   !!! Create HeatXfer context, get all HeatXfer options
   PetscCallA(MEF90HeatXferCtxCreate(MEF90HeatXferCtx,dm,MEF90Ctx,ierr))
   PetscCallA(MEF90HeatXferCtxSetFromOptions(MEF90HeatXferCtx,PETSC_NULL_CHARACTER,MEF90HeatXferDefaultGlobalOptions,MEF90HeatXferDefaultCellSetOptions,MEF90HeatXferDefaultFaceSetOptions,MEF90HeatXferDefaultVertexSetOptions,ierr))
   !!! We no longer need the DM. We have the megaDM in MEF90HeatXferCtx
   PetscCallA(DMDestroy(dm,ierr))
   PetscCallA(PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr))

   !!! Get parse all materials data from the command line
   PetscCallA(DMGetDimension(MEF90HeatXferCtx%megaDM,dim,ierr))
   If (dim == 2) Then
      PetscCallA(MEF90MatPropBagSetFromOptions(MEF90HeatXferCtx%MaterialPropertiesBag,MEF90HeatXferCtx%megaDM,MEF90Mathium2D,MEF90Ctx,ierr))
   Else
      PetscCallA(MEF90MatPropBagSetFromOptions(MEF90HeatXferCtx%MaterialPropertiesBag,MEF90HeatXferCtx%megaDM,MEF90Mathium3D,MEF90Ctx,ierr))
   End If   

   !!! Create GLOBAL vectors for the unknown (temperature), residuals, etc
   PetscCallA(VecGetDM(MEF90HeatXferCtx%temperatureLocal,temperatureDM,ierr)) 
   !!! This only borrows a reference so we do not need to delete it
   PetscCallA(DMCreateGlobalVector(temperatureDM,temperature,ierr))
   PetscCallA(PetscObjectSetName(temperature,"Temperature",ierr))
   PetscCallA(VecDuplicate(temperature,temperatureResidual,ierr))
   PetscCallA(PetscObjectSetName(temperatureResidual,"temperatureResidual",ierr))

   PetscCallA(MEF90CtxGetTime(MEF90Ctx,time,ierr))

   !!! 
   !!! Create SNES or TS, Mat and set KSP default options
   !!!
   If (MEF90HeatXferGlobalOptions%timeSteppingType == MEF90HeatXFer_timeSteppingTypeSteadyState) Then
      PetscCallA(MEF90HeatXferCreateSNES(MEF90HeatXferCtx,temperatureSNES,temperatureResidual,ierr))
   Else
      temperatureTSInitialStep = (time(size(time))-time(1)) / (size(time) - 1.0_Kr) / 10.0_Kr
      temperatureTSInitialStep = time(1)
      PetscCallA(MEF90HeatXferCreateTS(MEF90HeatXferCtx,temperatureTS,temperatureResidual,temperatureTSInitialStep,temperatureTSInitialStep,ierr))
      !PetscCallA(TSGetAdapt(temperatureTS,temperatureTSAdapt,ierr))
      !PetscCallA(TSAdaptSetFromOptions(temperatureTSAdapt,ierr))
   End If
   
   !!! 
   !!! Allocate array of works and energies
   !!!
   Allocate(energy(size(MEF90HeatXferCtx%CellSetOptionsBag)))
   Allocate(cellWork(size(MEF90HeatXferCtx%CellSetOptionsBag)))
   Allocate(faceWork(size(MEF90HeatXferCtx%FaceSetOptionsBag)))

   !!!
   !!! Actual computations / time stepping
   !!!
   If (MEF90GlobalOptions%timeSkip > 0) Then
      ! PetscCallA(DMGetLocalVector(MEF90HeatXferCtx%DMScal,localVec,ierr))
      ! PetscCallA(VecLoadExodusVertex(MEF90HeatXferCtx%DMScal,localVec,MEF90HeatXferCtx%MEF90Ctx%IOcomm, &
      !                          MEF90HeatXferCtx%MEF90Ctx%fileExoUnit,MEF90GlobalOptions%timeSkip,MEF90HeatXferGlobalOptions%TempOffset,ierr))
      ! PetscCallA(DMLocalToGlobalBegin(MEF90HeatXferCtx%DMScal,localVec,INSERT_VALUES,MEF90HeatXferCtx%Temperature,ierr))
      ! PetscCallA(DMLocalToGlobalEnd(MEF90HeatXferCtx%DMScal,localVec,INSERT_VALUES,MEF90HeatXferCtx%Temperature,ierr))
      ! PetscCallA(DMRestoreLocalVector(MEF90HeatXferCtx%DMScal,localVec,ierr))
      If (MEF90HeatXferGlobalOptions%timeSteppingType == MEF90HeatXFer_timeSteppingTypeTransient) Then
         PetscCallA(TSSetTime(temperatureTS,time(MEF90GlobalOptions%timeSkip),ierr))
      End If
   End If

   Do step = MEF90GlobalOptions%timeSkip+1,MEF90GlobalOptions%timeNumStep
      Select Case (MEF90HeatXferGlobalOptions%timeSteppingType)
      Case (MEF90HeatXFer_timeSteppingTypeSteadyState) 
         Write(IOBuffer,100) step,time(step)
         PetscCallA(PetscPrintf(MEF90Ctx%comm,IOBuffer,ierr))
         !!! Update fields
         PetscCallA(MEF90HeatXferUpdateTransients(MEF90HeatXferCtx,step,time(step),ierr))
         PetscCallA(DMLocalToGlobal(temperatureDM,MEF90HeatXferCtx%temperatureLocal,INSERT_VALUES,temperature,ierr))
         !!! Solve SNES
         PetscCallA(SNESSolve(temperatureSNES,PETSC_NULL_VEC,temperature,ierr))
      Case (MEF90HeatXFer_timeSteppingTypeTransient)
         Write(IOBuffer,200) step,time(step)
         PetscCallA(PetscPrintf(MEF90Ctx%comm,IOBuffer,ierr))
      !    If (step > 1) Then
      !       !!! Update fields
      !       PetscCallA(MEF90HeatXferSetTransients(MEF90HeatXferCtx,step,time(step),ierr))
      !       PetscCallA(MEF90HeatXferUpdateboundaryTemperature(MEF90HeatXferCtx%temperature,MEF90HeatXferCtx,ierr))
      !       !!! Make sure TS does not overstep
      !       PetscCallA(TSGetTime(temperatureTS,t,ierr))
      !       If (t < time(step)) Then
      !          !PetscCallA(TSAdaptSetStepLimits(tsAdaptTemp,PETSC_DECIDE,(time(step)-time)/2.0_Kr,ierr))
      !          !!! Something is up here. 
      !          !!! replacing the constant 10000 with a variable leads to divergence of TSAdapt
      !          !!! when using gcc
      !          PetscCallA(TSSetDuration(temperatureTS,10000,time(step),ierr))
      !          PetscCallA(TSSolve(temperatureTS,MEF90HeatXferCtx%temperature,time(step),ierr))
      !          PetscCallA(TSGetTime(temperatureTS,t,ierr))
      !          time(step) = t
      !       Else
      !          Write(IOBuffer,*) 'TS exceeded analysis time. Skipping step\n'
      !          PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
      !       End If
      !    End If
      End Select

      !!! Compute energies
      PetscCallA(MEF90HeatXFerEnergy(MEF90HeatXferCtx,energy,cellWork,faceWork,ierr))
      PetscCall(DMGetLabelIdIS(temperatureDM,MEF90CellSetLabelName,setIS,ierr))
      Call MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm,setIS,ierr);CHKERRQ(ierr)
      PetscCallA(ISGetIndicesF90(setIS,setID,ierr))
      Do set = 1, size(setID)
         Write(IOBuffer,101) setID(set),energy(set),cellWork(set),energy(set)-cellWork(set)
         PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))
      End Do
      PetscCallA(ISRestoreIndicesF90(setIS,setID,ierr))
      PetscCallA(ISDestroy(setIS,ierr))

      PetscCall(DMGetLabelIdIS(temperatureDM,MEF90FaceSetLabelName,setIS,ierr))
      Call MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm,setIS,ierr);CHKERRQ(ierr)
      PetscCallA(ISGetIndicesF90(setIS,setID,ierr))
      Do set = 1, size(setID)
         Write(IOBuffer,103) setID(set),faceWork(set)
         PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))
      End Do
      PetscCallA(ISRestoreIndicesF90(setIS,setID,ierr))
      PetscCallA(ISDestroy(setIS,ierr))

      Write(IOBuffer,102) sum(energy),sum(cellWork)+sum(faceWork),sum(energy)-sum(cellWork)-sum(faceWork)
      PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))
      !!! Save results
      PetscCallA(MEF90HeatXferViewEXO(MEF90HeatXferCtx,step,ierr))
   End Do
100 Format("Solving steady state step ",I4,", t=",ES12.5,"\n")
200 Format("Solving transient step ",I4,", t=",ES12.5,"\n")
101 Format("cell set ",I4," thermal energy: ",ES12.5," flux work: ",ES12.5," total: ",ES12.5,"\n")
103 Format("face set ",I4,"                              flux work: ",ES12.5,"\n")
102 Format("======= Total thermal energy: ",ES12.5," flux work: ",ES12.5," total: ",ES12.5,"\n")
   
   DeAllocate(time)
   DeAllocate(energy)
   DeAllocate(cellWork)
   DeAllocate(faceWork)
   PetscCallA(PetscViewerASCIIOpen(MEF90Ctx%comm,trim(MEF90FilePrefix(MEF90Ctx%resultFile))//'.log',logViewer, ierr))
   PetscCallA(PetscLogView(logViewer,ierr))
   PetscCallA(PetscViewerFlush(logViewer,ierr))
   PetscCallA(PetscViewerDestroy(logViewer,ierr))

   !!! Clean up and exit nicely
   If (MEF90HeatXferGlobalOptions%timeSteppingType == MEF90HeatXFer_timeSteppingTypeSteadyState) Then
      PetscCallA(SNESDestroy(temperatureSNES,ierr))
   Else
      PetscCallA(TSDestroy(temperatureTS,ierr))
   End If

   PetscCallA(VecDestroy(temperatureResidual,ierr))
   PetscCallA(VecDestroy(temperature,ierr))
   PetscCallA(DMDestroy(temperatureDM,ierr))
   PetscCallA(MEF90HeatXferCtxDestroy(MEF90HeatXferCtx,ierr))
   PetscCallA(MEF90CtxDestroy(MEF90Ctx,ierr))
   PetscCallA(MEF90Finalize(ierr))
   PetscCallA(PetscFinalize(ierr))
End Program HeatXfer
