#include "../MEF90/mef90.inc"
Program HeatXfer
   #include "petsc/finclude/petsc.h"
   Use m_MEF90
   Use m_MEF90_HeatXfer
   Use m_MEF90_HeatXferCtx
   Use petsc
   Implicit NONE   

   PetscErrorCode                                     :: ierr
   Type(MEF90Ctx_Type),target                         :: MEF90Ctx
   Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions
   Type(MEF90CtxGlobalOptions_Type),Parameter         :: MEF90DefaultGlobalOptions = MEF90CtxGlobalOptions_Type( &
                                                         1,                               & ! verbose
                                                         PETSC_FALSE,                     & ! dryrun
                                                         MEF90TimeInterpolation_linear,   & ! timeInterpolation
                                                         0.0_Kr,                          & ! timeMin
                                                         1.0_Kr,                          & ! timeMax
                                                         11,                              & ! timeNumStep
                                                         0,                               & ! timeSkip
                                                         1.0_Kr,                          & ! timeNumCycle
                                                         MEF90ElementFamilyLagrange,      & ! elementFamily
                                                         1)                                 ! elementOrder

   Type(MEF90HeatXferCtx_Type)                        :: MEF90HeatXferCtx
   Type(MEF90HeatXferGlobalOptions_Type),Pointer      :: MEF90HeatXferGlobalOptions
   Type(MEF90HeatXferGlobalOptions_Type),Parameter    :: MEF90HeatXferDefaultGlobalOptions = MEF90HeatXferGlobalOptions_Type( &
                                                         MEF90HeatXFer_timeSteppingTypeSteadyState, & ! timeSteppingType
                                                         PETSC_FALSE,         & ! addNullSpace
                                                         0.,                  & ! initialTemperature
                                                         MEF90Scaling_Linear, & ! boundaryTempScaling
                                                         MEF90Scaling_Linear, & ! externalTempScaling
                                                         MEF90Scaling_Linear, & ! fluxScaling
                                                         MEF90Scaling_Linear)   ! boundaryFluxScaling

   Type(MEF90HeatXferCellSetOptions_Type),Parameter   :: MEF90HeatXferDefaultCellSetOptions = MEF90HeatXferCellSetOptions_Type( &
                                                         0.0_Kr,        & ! flux
                                                         0.0_Kr,        & ! surfaceThermalConductivity
                                                         0.0_Kr,        & ! externalTemp
                                                         PETSC_FALSE,   & ! Has BC
                                                         0.0_Kr,        & ! boundaryTemperature
                                                         [0.0_Kr,0.0_Kr,0.0_Kr]) ! AdvectionVector
                                                         
   Type(MEF90HeatXferVertexSetOptions_Type),Parameter :: MEF90HeatXferDefaultVertexSetOptions = MEF90HeatXferVertexSetOptions_Type( &
                                                         PETSC_FALSE,   & ! Has BC
                                                         0.0_Kr)          ! boundaryTemp
                                                         
   Type(tDM),target                                   :: dm
   Type(tIS)                                          :: setIS,cellIS,CellSetGlobalIS
   PetscInt,Dimension(:),Pointer                      :: setID
   PetscInt                                           :: numset,set
   PetscReal,Dimension(:),Pointer                     :: time,energy,work

          
   PetscBool                                          :: flg
   Character(len=MEF90MXSTRLEN)                       :: IOBuffer
   Type(tPetscViewer)                                 :: logViewer
   Integer                                            :: numfield
   
   Type(tSNES)                                        :: snesTemp
   Type(tTS)                                          :: tsTemp
   Type(tTSAdapt)                                     :: tsAdaptTemp
   Type(tVec)                                         :: residualTemp,localVec

   PetscReal                                          :: tsTempInitialStep,tsTempInitialTime
   PetscInt                                           :: tsTempmaxIter
   PetscReal                                          :: t
   
   Integer                                            :: step
   PetscInt                                           :: dim
      
   !!! Initialize MEF90
   PetscCall(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
   PetscCall(MEF90Initialize(ierr))

   !!! Get all MEF90-wide options
   PetscCall(MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90DefaultGlobalOptions,ierr))
   PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr))

   PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm,MEF90Ctx%geometryfile,PETSC_NULL_CHARACTER,interpolate,dm,ierr))
   PetscCallA(DMPlexDistributeSetDefault(dm,PETSC_FALSE,ierr))
   PetscCallA(DMSetFromOptions(dm,ierr))
   distribute: Block 
       Type(tDM),target                    :: dmDist
       PetscInt                            :: ovlp = 0
       If (MEF90Ctx%NumProcs > 1) Then
           PetscCallA(DMPlexDistribute(dm,ovlp,PETSC_NULL_SF,dmDist,ierr))
           PetscCallA(DMDestroy(dm,ierr))
           dm = dmDist
       End If
   End Block distribute

   !!! Open output file
   !PetscCall(MEF90CtxOpenEXO(MEF90Ctx,Mesh,ierr))

   
   !!! Create HeatXfer context, get all HeatXfer options
   PetscCall(MEF90HeatXferCtxCreate(MEF90HeatXferCtx,Mesh,MEF90Ctx,ierr))
   PetscCall(MEF90HeatXferCtxSetFromOptions(MEF90HeatXferCtx,PETSC_NULL_CHARACTER,MEF90HeatXferDefaultGlobalOptions, &
                                       MEF90HeatXferDefaultCellSetOptions,MEF90HeatXferDefaultVertexSetOptions,ierr))
   PetscCall(PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr))

   !!! Get material properties bags
   PetscCall(DMMeshGetDimension(Mesh,dim,ierr))
   If (dim == 2) Then
      PetscCall(MEF90MatPropBagSetFromOptions(MEF90HeatXferCtx%MaterialPropertiesBag,MEF90HeatXferCtx%DM,MEF90Mathium2D,MEF90Ctx,ierr))
   Else
      PetscCall(MEF90MatPropBagSetFromOptions(MEF90HeatXferCtx%MaterialPropertiesBag,MEF90HeatXferCtx%DM,MEF90Mathium3D,MEF90Ctx,ierr))
   End If   

   PetscCall(MEF90CtxGetTime(MEF90Ctx,time,ierr))

   !!! 
   !!! Create SNES or TS, Mat and set KSP default options
   !!!

   PetscCall(VecDuplicate(MEF90HeatXferCtx%temperature,residualTemp,ierr))
   PetscCall(PetscObjectSetName(residualTemp,"residualTemperature",ierr))
   If (MEF90HeatXferGlobalOptions%timeSteppingType == MEF90HeatXFer_timeSteppingTypeSteadyState) Then
      PetscCall(MEF90HeatXferCreateSNES(MEF90HeatXferCtx,snesTemp,residualTemp,ierr)
   Else
      tsTempInitialStep = (time(size(time))-time(1)) / (size(time) - 1.0_Kr) / 10.0_Kr
      tsTempInitialTime = time(1)
      PetscCall(MEF90HeatXferCreateTS(MEF90HeatXferCtx,tsTemp,residualTemp,tsTempInitialTime,tsTempInitialStep,ierr)
      PetscCall(TSGetAdapt(tsTemp,tsAdaptTemp,ierr))
      PetscCall(TSAdaptSetFromOptions(tsAdaptTemp,ierr))
   End If
   
   !!! 
   !!! Allocate array of works and energies
   !!!
   Allocate(energy(size(MEF90HeatXferCtx%CellSetOptionsBag)))
   Allocate(work(size(MEF90HeatXferCtx%CellSetOptionsBag)))

   !!!
   !!! Try to figure out if the file was formatted
   !!!
   If (MEF90Ctx%rank == 0) Then
      PetscCall(EXGVP(MEF90Ctx%fileExoUnit,"N",numfield,ierr)
   End If
   PetscCall(MPI_Bcast(numfield,1,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)

   If (numfield == 0) Then
      PetscCall(MEF90HeatXferFormatEXO(MEF90HeatXferCtx,ierr)
   End If
   
   !!!
   !!! Actual computations / time stepping
   !!!
   If (MEF90GlobalOptions%timeSkip > 0) Then
      PetscCall(DMGetLocalVector(MEF90HeatXferCtx%DMScal,localVec,ierr))
      PetscCall(VecLoadExodusVertex(MEF90HeatXferCtx%DMScal,localVec,MEF90HeatXferCtx%MEF90Ctx%IOcomm, &
                               MEF90HeatXferCtx%MEF90Ctx%fileExoUnit,MEF90GlobalOptions%timeSkip,MEF90HeatXferGlobalOptions%TempOffset,ierr))
      PetscCall(DMLocalToGlobalBegin(MEF90HeatXferCtx%DMScal,localVec,INSERT_VALUES,MEF90HeatXferCtx%Temperature,ierr))
      PetscCall(DMLocalToGlobalEnd(MEF90HeatXferCtx%DMScal,localVec,INSERT_VALUES,MEF90HeatXferCtx%Temperature,ierr))
      PetscCall(DMRestoreLocalVector(MEF90HeatXferCtx%DMScal,localVec,ierr))
      If (MEF90HeatXferGlobalOptions%timeSteppingType == MEF90HeatXFer_timeSteppingTypeTransient) Then
         PetscCall(TSSetTime(tsTemp,time(MEF90GlobalOptions%timeSkip),ierr))
      End If
   End If

   Do step = MEF90GlobalOptions%timeSkip+1,MEF90GlobalOptions%timeNumStep

      Select Case (MEF90HeatXferGlobalOptions%timeSteppingType)
      Case (MEF90HeatXFer_timeSteppingTypeSteadyState) 
         Write(IOBuffer,100) step,time(step)
         PetscCall(PetscPrintf(MEF90Ctx%comm,IOBuffer,ierr))
         !!! Update fields
         PetscCall(MEF90HeatXferSetTransients(MEF90HeatXferCtx,step,time(step),ierr))
         !!! Solve SNES
         PetscCall(MEF90HeatXferUpdateboundaryTemperature(MEF90HeatXferCtx%temperature,MEF90HeatXferCtx,ierr))
         PetscCall(SNESSolve(snesTemp,PETSC_NULL_OBJECT,MEF90HeatXferCtx%temperature,ierr))
      Case (MEF90HeatXFer_timeSteppingTypeTransient)
         Write(IOBuffer,200) step,time(step)
         PetscCall(PetscPrintf(MEF90Ctx%comm,IOBuffer,ierr))
         If (step > 1) Then
            !!! Update fields
            PetscCall(MEF90HeatXferSetTransients(MEF90HeatXferCtx,step,time(step),ierr)
            PetscCall(MEF90HeatXferUpdateboundaryTemperature(MEF90HeatXferCtx%temperature,MEF90HeatXferCtx,ierr);
            !!! Make sure TS does not overstep
            PetscCall(TSGetTime(tsTemp,t,ierr))
            If (t < time(step)) Then
               !PetscCall(TSAdaptSetStepLimits(tsAdaptTemp,PETSC_DECIDE,(time(step)-time)/2.0_Kr,ierr))
               !!! Something is up here. 
               !!! replacing the constant 10000 with a variable leads to divergence of TSAdapt
               !!! when using gcc
               PetscCall(TSSetDuration(tsTemp,10000,time(step),ierr))
               PetscCall(TSSolve(tsTemp,MEF90HeatXferCtx%temperature,time(step),ierr))
               PetscCall(TSGetTime(tsTemp,t,ierr))
               time(step) = t
            Else
               Write(IOBuffer,*) 'TS exceeded analysis time. Skipping step\n'
               PetscCall(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
            End If
         End If
      End Select
      !!! Compute energies
      PetscCall(MEF90HeatXFerEnergy(MEF90HeatXferCtx%temperature,time(step),MEF90HeatXferCtx,energy,work,ierr))
      PetscCall(DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Cell Sets',CellSetGlobalIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr)) 
      PetscCall(ISGetIndicesF90(CellSetGlobalIS,setID,ierr))
      Do set = 1, size(setID)
         Write(IOBuffer,101) setID(set),energy(set),work(set),energy(set)-work(set)
         PetscCall(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))
      End Do
      PetscCall(ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr))
      PetscCall(ISDestroy(CellSetGlobalIS,ierr))
      Write(IOBuffer,102) sum(energy),sum(work),sum(energy)-sum(work)
      PetscCall(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))
      !!! Save results
      PetscCall(MEF90HeatXferViewEXO(MEF90HeatXferCtx,step,ierr)
   End Do
100 Format("Solving steady state step ",I4,", t=",ES12.5,"\n")
200 Format("Solving transient step ",I4,", t=",ES12.5,"\n")
101 Format("cell set ",I4," thermal energy: ",ES12.5," fluxes work: ",ES12.5," total: ",ES12.5,"\n")
102 Format("======= Total thermal energy: ",ES12.5," fluxes work: ",ES12.5," total: ",ES12.5,"\n")
   !!! Clean up and exit nicely
   If (MEF90HeatXferGlobalOptions%timeSteppingType == MEF90HeatXFer_timeSteppingTypeSteadyState) Then
      PetscCall(SNESDestroy(snesTemp,ierr))
   Else
      PetscCall(TSDestroy(tsTemp,ierr))
   End If
   PetscCall(VecDestroy(residualTemp,ierr))
   PetscCall(MEF90HeatXferCtxDestroyVectors(MEF90HeatXferCtx,ierr))
   
   DeAllocate(time)
   DeAllocate(energy)
   DeAllocate(work)
   PetscCall(PetscViewerASCIIOpen(MEF90Ctx%comm,trim(MEF90FilePrefix(MEF90Ctx%resultFile))//'.log',logViewer, ierr))
   PetscCall(PetscLogView(logViewer,ierr))
   PetscCall(PetscViewerFlush(logViewer,ierr))
   PetscCall(PetscViewerDestroy(logViewer,ierr))
   PetscCall(MEF90HeatXferCtxDestroy(MEF90HeatXferCtx,ierr))
   PetscCall(MEF90CtxCloseEXO(MEF90Ctx,ierr))
   PetscCall(MEF90CtxDestroy(MEF90Ctx,ierr))
   PetscCall(MEF90Finalize(ierr))
   PetscCall(PetscFinalize(ierr))
End Program HeatXfer
