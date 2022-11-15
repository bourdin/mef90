#include "../MEF90/mef90.inc"
Program ThermoElasticity
#include "petsc/finclude/petsc.h"
   Use m_MEF90
   Use m_MEF90_DefMechCtx
   Use m_MEF90_DefMech
   Use m_MEF90_HeatXfer
   Use m_MEF90_HeatXferCtx
   Use m_vDefDefault
   Use petsc
   Implicit NONE

   PetscErrorCode                                     :: ierr
   Type(MEF90Ctx_Type),target                         :: MEF90Ctx
   Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions

   !!! Defect mechanics contexts
   Type(MEF90DefMechCtx_Type)                         :: MEF90DefMechCtx
   Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
   !!! HeatXfer contexts
   Type(MEF90HeatXferCtx_Type)                        :: MEF90HeatXferCtx
   Type(MEF90HeatXferGlobalOptions_Type),Pointer      :: MEF90HeatXferGlobalOptions

   Type(tDM),target                                   :: dm,temperatureDM,displacementDM
   Type(tIS)                                          :: setIS
   PetscInt,Dimension(:),Pointer                      :: setID
   PetscInt                                           :: set
   PetscReal,Dimension(:),Pointer                     :: time,energy,cellWork,faceWork

   Type(tSNES)                                        :: SNESDisplacement
   SNESConvergedReason                                :: SNESDisplacementConvergedReason
   Type(tVec)                                         :: displacement,displacementResidual

   Type(tSNES)                                        :: SNESTemperature
   Type(tTS)                                          :: TSTemperature
   Type(tTSAdapt)                                     :: TSAdaptTemperature
   Type(tVec)                                         :: temperature,temperatureResidual

   PetscReal                                          :: tsTemperatureInitialTimeStep,tsTemperatureInitialTime
   PetscInt                                           :: tsTemperatureMaxIter
   PetscReal                                          :: t

   PetscBool                                          :: flg
   Character(len=MEF90MXSTRLEN)                       :: IOBuffer
   Type(tPetscViewer)                                 :: logViewer
   PetscInt                                           :: numfield

   PetscInt                                           :: step
   PetscInt                                           :: dim

   !!! Initialize MEF90
   PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
   PetscCallA(MEF90Initialize(ierr))

   !!! Get all MEF90-wide options
   PetscCallA(MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90CtxDefaultGlobalOptions,ierr))
   PetscCallA(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr))

   PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm,MEF90Ctx%geometryFile,PETSC_NULL_CHARACTER,PETSC_TRUE,dm,ierr))
   PetscCallA(DMPlexDistributeSetDefault(dm,PETSC_FALSE,ierr))
   PetscCallA(DMSetUseNatural(dm,PETSC_TRUE,ierr))
   PetscCallA(DMSetFromOptions(dm,ierr))
   PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-mef90_dm_view",ierr))

   PetscCallA(MEF90CtxGetTime(MEF90Ctx,time,ierr))

   Inquire(file=MEF90Ctx%resultFile,exist=flg)
   If (flg) Then
      ! we assume that the output file exists and is formatted
      PetscCallA(MEF90CtxOpenEXO(MEF90Ctx,MEF90Ctx%resultViewer,FILE_MODE_APPEND,ierr))
   Else
      ! we need to create the output file
      EXOFormat: block
         PetscInt                                            :: numNodalVar = 2, numCellVar = 3, numFaceVar = 3, numGVar = 0
         Character(len=MEF90MXSTRLEN),Dimension(:),Pointer   :: faceVarName, nodalVarName, cellVarName, gVarName

         Allocate(nodalVarName(numNodalVar))
         Allocate(cellVarName(numCellVar))
         Allocate(faceVarName(numFaceVar))
         Allocate(gVarName(numGVar))
         nodalVarName = ["Displacement       ","Temperature        "]
         cellVarName  = ["bodyForce          ","Stress             ","HeatFlux           "]
         faceVarName  = ["boundaryForce      ","boundaryFlux       ","externalTemperature"]
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
       PetscInt                            :: ovlp = 0_Ki
       Type(tPetscSF)                      :: naturalPointSF

       If (MEF90Ctx%NumProcs > 1) Then
           PetscCallA(DMPlexDistribute(dm,ovlp,naturalPointSF,dmDist,ierr))
           PetscCallA(DMPlexSetMigrationSF(dmDist,naturalPointSF, ierr))
           PetscCallA(PetscSFDestroy(naturalPointSF,ierr))
           PetscCallA(DMDestroy(dm,ierr))
           dm = dmDist
       End If
   End Block distribute
   PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-mef90_dm_view",ierr))

   !!! Create HeatXfer context, get all HeatXfer options
   ! PetscCallA(MEF90HeatXferCtxCreate(MEF90HeatXferCtx,dm,MEF90Ctx,ierr))
   ! PetscCallA(MEF90HeatXferCtxSetFromOptions(MEF90HeatXferCtx,PETSC_NULL_CHARACTER,HeatXferDefaultGlobalOptions,HeatXferDefaultCellSetOptions,HeatXferDefaultFaceSetOptions,HeatXferDefaultVertexSetOptions,ierr))

   !!! Create DefMechCtx, get anll defMech options
   PetscCallA(MEF90DefMechCtxCreate(MEF90DefMechCtx,dm,MEF90Ctx,ierr))
   PetscCallA(MEF90DefMechCtxSetFromOptions(MEF90DefMechCtx,PETSC_NULL_CHARACTER,DefMechDefaultGlobalOptions,DefMechDefaultCellSetOptions,DefMechDefaultFaceSetOptions,DefMechDefaultVertexSetOptions,ierr))

   !!! We no longer need the DM. We have the megaDM in MEF90HeatXferCtx and MEF90DefMechCtx
   PetscCallA(DMDestroy(dm,ierr))
   PetscCallA(PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr))

   !!! Get parse all materials data from the command line
   PetscCallA(DMGetDimension(MEF90HeatXferCtx%megaDM,dim,ierr))
   If (dim == 2) Then
      PetscCallA(MEF90MatPropBagSetFromOptions(MEF90HeatXferCtx%MaterialPropertiesBag,MEF90HeatXferCtx%megaDM,MEF90Mathium2D,MEF90Ctx,ierr))
   Else
      PetscCallA(MEF90MatPropBagSetFromOptions(MEF90HeatXferCtx%MaterialPropertiesBag,MEF90HeatXferCtx%megaDM,MEF90Mathium3D,MEF90Ctx,ierr))
   End If   

   !!! Create GLOBAL vectors for the unknowns (temperature,displacements), residuals, etc
   PetscCallA(VecGetDM(MEF90HeatXferCtx%temperatureLocal,temperatureDM,ierr)) 
   !!! This only borrows a reference so we do not need to delete it
   PetscCallA(DMCreateGlobalVector(temperatureDM,temperature,ierr))
   PetscCallA(PetscObjectSetName(temperature,"Temperature",ierr))
   PetscCallA(VecDuplicate(temperature,temperatureResidual,ierr))
   PetscCallA(PetscObjectSetName(temperatureResidual,"temperatureResidual",ierr))

   PetscCallA(VecGetDM(MEF90DefMechCtx%displacementLocal,displacementDM,ierr)) 
   !!! This only borrows a reference so we do not need to delete it
   PetscCallA(DMCreateGlobalVector(displacementDM,displacement,ierr))
   PetscCallA(PetscObjectSetName(displacement,"displacement",ierr))
   PetscCallA(VecDuplicate(displacement,displacementResidual,ierr))
   PetscCallA(PetscObjectSetName(displacementResidual,"displacementResidual",ierr))

   !!! 
   !!! Create SNES or TS, Mat and set KSP default options
   !!!
   If (MEF90HeatXferGlobalOptions%timeSteppingType == MEF90HeatXFer_timeSteppingTypeSteadyState) Then
      PetscCallA(MEF90HeatXferCreateSNES(MEF90HeatXferCtx,SNESTemperature,temperatureResidual,ierr))
   Else
      TSTemperatureInitialTimeStep = (time(size(time))-time(1)) / (size(time) - 1.0_Kr) / 10.0_Kr
      TSTemperatureInitialTime = time(1)
      PetscCallA(MEF90HeatXferCreateTS(MEF90HeatXferCtx,TSTemperature,temperatureResidual,TSTemperatureInitialTime,TSTemperatureInitialTimeStep,ierr))
      !PetscCallA(TSGetAdapt(temperatureTS,temperatureTSAdapt,ierr))
      !PetscCallA(TSAdaptSetFromOptions(temperatureTSAdapt,ierr))
   End If
   
   !!! 
   !!! Allocate array of works and energies
   !!!
   Allocate(energy(size(MEF90HeatXferCtx%CellSetOptionsBag)))
   Allocate(cellWork(size(MEF90HeatXferCtx%CellSetOptionsBag)))
   Allocate(faceWork(size(MEF90HeatXferCtx%FaceSetOptionsBag)))


   !!! -------------------------------------



   !!!
   !!! Actual computations / time stepping
   !!!
   If (MEF90DefMechGlobalOptions%timeSteppingType == MEF90DefMech_timeSteppingTypeQuasiStatic) Then
      Do step = 1,MEF90GlobalOptions%timeNumStep
         Write(IOBuffer,100) step,time(step)
         PetscCallA(PetscPrintf(MEF90Ctx%comm,IOBuffer,ierr))

         !!! Solve for temperature
         Select Case (MEF90HeatXferGlobalOptions%timeSteppingType)
         Case (MEF90HeatXfer_timeSteppingTypeSteadyState)
            !!! Update fields
            ! PetscCallA(MEF90HeatXferSetTransients(MEF90HeatXferCtx,step,time(step),ierr))
            ! !!! Solve SNES
            ! PetscCallA(MEF90HeatXferUpdateboundaryTemperature(MEF90HeatXferCtx%temperature,MEF90HeatXferCtx,ierr))
            ! PetscCallA(SNESSolve(snesTemp,PETSC_NULL_OBJECT,MEF90HeatXferCtx%temperature,ierr))

            ! !!! Compute thermal energy
            ! PetscCallA(MEF90HeatXFerEnergy(MEF90HeatXferCtx%temperature,time(step),MEF90HeatXferCtx,energy,work,ierr))
            ! PetscCallA(DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Cell Sets',CellSetGlobalIS,ierr))
            ! PetscCallA(MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr))
            ! PetscCallA(ISGetIndicesF90(CellSetGlobalIS,setID,ierr))
            ! Do set = 1, size(setID)
            !    Write(IOBuffer,101) setID(set),energy(set),work(set),energy(set)-work(set)
            !    PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))
            ! End Do
            ! PetscCallA(ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr))
            ! PetscCallA(ISDestroy(CellSetGlobalIS,ierr))
            ! Write(IOBuffer,102) sum(energy),sum(work),sum(energy)-sum(work)
            ! PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))
            ! !!! Save results
            ! PetscCallA(MEF90HeatXferViewEXO(MEF90HeatXferCtx,step,ierr))
         Case (MEF90HeatXfer_timeSteppingTypeTransient)
            If (step > 1) Then
               ! !!! Update fields
               ! PetscCallA(MEF90HeatXferSetTransients(MEF90HeatXferCtx,step,time(step),ierr))
               ! PetscCallA(MEF90HeatXferUpdateboundaryTemperature(MEF90HeatXferCtx%temperature,MEF90HeatXferCtx,ierr))
               ! !!! Make sure TS does not overstep
               ! PetscCallA(TSGetTime(tsTemp,t,ierr))
               ! If (t < time(step)) Then
               !    PetscCallA(TSAdaptSetStepLimits(tsAdaptTemp,PETSC_DECIDE,(time(step)-time)/2.0_Kr,ierr))
               !    !!! Something is up here.
               !    !!! replacing the constant 10000 with a variable leads to divergence of TSAdapt
               !    !!! when using gcc
               !    PetscCallA(TSSetDuration(tsTemp,10000,time(step),ierr))
               !    PetscCallA(TSSolve(tsTemp,MEF90HeatXferCtx%temperature,time(step),ierr))
               !    PetscCallA(TSGetTime(tsTemp,t,ierr))
               !    time(step) = t
               ! Else
               !    Write(IOBuffer,*) 'TS exceeded analysis time. Skipping step\n'
               !    PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
               ! End If
            End If

            !!! Compute thermal energy
            ! PetscCallA(MEF90HeatXFerEnergy(MEF90HeatXferCtx%temperature,time(step),MEF90HeatXferCtx,energy,work,ierr))
            ! PetscCallA(DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Cell Sets',CellSetGlobalIS,ierr))
            ! PetscCallA(MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr))
            ! PetscCallA(ISGetIndicesF90(CellSetGlobalIS,setID,ierr))
            ! Do set = 1, size(setID)
            !    Write(IOBuffer,101) setID(set),energy(set),work(set),energy(set)-work(set)
            !    PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))
            ! End Do
            ! PetscCallA(ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr))
            ! PetscCallA(ISDestroy(CellSetGlobalIS,ierr))
            ! Write(IOBuffer,102) sum(energy),sum(work),sum(energy)-sum(work)
            ! PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))
            ! !!! Save results
            ! PetscCallA(MEF90HeatXferViewEXO(MEF90HeatXferCtx,step,ierr))
         End Select

         !!! Solve for displacement
         Select case(MEF90DefMechGlobalOptions%timeSteppingType)
         Case (MEF90DefMech_timeSteppingTypeQuasiStatic)
            ! !!! Update fields
            ! PetscCallA(MEF90DefMechSetTransients(MEF90DefMechCtx,step,time(step),ierr))
            ! PetscCallA(MEF90DefMechUpdateboundaryDisplacement(MEF90DefMechCtx%displacement,MEF90DefMechCtx,ierr))

            ! !!! Solve SNES
            ! PetscCallA(SNESSolve(snesDisp,PETSC_NULL_OBJECT,MEF90DefMechCtx%displacement,ierr))
            ! PetscCallA(SNESGetConvergedReason(snesDisp,snesDispConvergedReason,ierr))
            ! Write(IOBuffer,*) "SNESConvergedReason returned ",snesDispConvergedReason,"\n"
            ! PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))

            ! !!! Compute energies
            ! energy = 0.0_Kr
            ! work = 0.0_Kr
            ! PetscCallA(MEF90DefMechWork(MEF90DefMechCtx%displacement,MEF90DefMechCtx,work,ierr))
            ! PetscCallA(MEF90DefMechElasticEnergy(MEF90DefMechCtx%displacement,MEF90DefMechCtx,energy,ierr))
            ! PetscCallA(DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Cell Sets',CellSetGlobalIS,ierr))
            ! PetscCallA(MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr))
            ! PetscCallA(ISGetIndicesF90(CellSetGlobalIS,setID,ierr))
            ! Do set = 1, size(setID)
            !    Write(IOBuffer,201) setID(set),energy(set),work(set),energy(set)-work(set)
            !    PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))
            ! End Do
            ! PetscCallA(ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr))
            ! PetscCallA(ISDestroy(CellSetGlobalIS,ierr))
            ! Write(IOBuffer,202) sum(energy),sum(work),sum(energy)-sum(work)
            ! PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))

            ! !!! Save results and boundary Values
            ! PetscCallA(MEF90DefMechViewEXO(MEF90DefMechCtx,step,ierr))
         End Select
      End Do
   End If
100 Format("Solving steady state step ",I4,", t=",ES12.5,"\n")
101 Format("cell set ",I4," thermal energy: ",ES12.5," fluxes work: ",ES12.5," total: ",ES12.5,"\n")
102 Format("======= Total thermal energy: ",ES12.5," fluxes work: ",ES12.5," total: ",ES12.5,"\n")
201 Format("cell set ",I4," elastic energy: ",ES12.5," work: ",ES12.5," total: ",ES12.5,"\n")
202 Format("======= Total elastic energy: ",ES12.5," work: ",ES12.5," total: ",ES12.5,"\n")

   !!! Clean up and exit nicely
   Select case(MEF90DefMechGlobalOptions%timeSTeppingType)
   Case (MEF90DefMech_timeSTeppingTypeQuasiStatic)
      PetscCallA(SNESDestroy(SNESDisplacement,ierr))
      PetscCallA(VecDestroy(displacementResidual,ierr))
      PetscCallA(VecDestroy(displacement,ierr))
   End Select

   !!! Clean up and exit nicely
   If (MEF90HeatXferGlobalOptions%timeSteppingType == MEF90HeatXFer_timeSteppingTypeSteadyState) Then
      PetscCallA(SNESDestroy(SNESTemperature,ierr))
   Else
      PetscCallA(TSDestroy(TSTemperature,ierr))
   End If
   PetscCallA(VecDestroy(temperatureResidual,ierr))
   PetscCallA(VecDestroy(temperature,ierr))

   DeAllocate(time)
   DeAllocate(energy)
   DeAllocate(cellWork)
   DeAllocate(faceWork)
   PetscCallA(MEF90DefMechCtxDestroy(MEF90DefMechCtx,ierr))
   PetscCallA(MEF90HeatXferCtxDestroy(MEF90HeatXferCtx,ierr))

   
   PetscCallA(PetscViewerASCIIOpen(MEF90Ctx%comm,trim(MEF90FilePrefix(MEF90Ctx%resultFile))//'.log',logViewer, ierr))
   PetscCallA(PetscLogView(logViewer,ierr))
   PetscCallA(PetscViewerDestroy(logViewer,ierr))


   PetscCallA(MEF90CtxDestroy(MEF90Ctx,ierr))
   PetscCallA(MEF90Finalize(ierr))
   PetscCallA(PetscFinalize(ierr))
End Program ThermoElasticity
