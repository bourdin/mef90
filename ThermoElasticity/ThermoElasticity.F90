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

   Type(tSNES)                                        :: displacementSNES
   !SNESConvergedReason                                :: SNESDisplacementConvergedReason
   Type(tVec)                                         :: displacement,displacementResidual

   Type(tSNES)                                        :: temperatureSNES
   !SNESConvergedReason                                :: SNESTemperatureConvergedReason
   Type(tTS)                                          :: temperatureTS
   Type(tTSAdapt)                                     :: temperatureTSAdapt
   Type(tVec)                                         :: temperature,temperatureResidual

   PetscReal                                          :: temperatureInitialTimeStep,temperatureInitialTime
   !PetscInt                                           :: tsTemperatureMaxIter

   PetscBool                                          :: flg
   Character(len=MEF90MXSTRLEN)                       :: IOBuffer
   Type(tPetscViewer)                                 :: logViewer

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
   PetscCallA(MEF90HeatXferCtxCreate(MEF90HeatXferCtx,dm,MEF90Ctx,ierr))
   PetscCallA(MEF90HeatXferCtxSetFromOptions(MEF90HeatXferCtx,PETSC_NULL_CHARACTER,HeatXferDefaultGlobalOptions,HeatXferDefaultCellSetOptions,HeatXferDefaultFaceSetOptions,HeatXferDefaultVertexSetOptions,ierr))
   PetscCallA(PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr))

   !!! Create DefMechCtx, get all defMech options
   PetscCallA(MEF90DefMechCtxCreate(MEF90DefMechCtx,dm,MEF90Ctx,ierr))
   PetscCallA(MEF90DefMechCtxSetFromOptions(MEF90DefMechCtx,PETSC_NULL_CHARACTER,DefMechDefaultGlobalOptions,DefMechDefaultCellSetOptions,DefMechDefaultFaceSetOptions,DefMechDefaultVertexSetOptions,ierr))
   PetscCallA(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr))

   !!! We no longer need the DM. We have the megaDM in MEF90HeatXferCtx and MEF90DefMechCtx
   PetscCallA(DMDestroy(dm,ierr))

   ! !!! DeAllocate MEF90DefMechCtx%temperatureLocal and make it point to MEF90HeatXferCtx%temperatureLocal
   ! PetscCall(VecDestroy(MEF90DefMechCtx%temperatureLocal,ierr))
   ! DeAllocate(MEF90DefMechCtx%temperatureLocal)
   ! Nullify(MEF90DefMechCtx%temperatureLocal)
   ! MEF90DefMechCtx%temperatureLocal => MEF90HeatXferCtx%temperatureLocal

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
      PetscCallA(MEF90HeatXferCreateSNES(MEF90HeatXferCtx,temperatureSNES,temperatureResidual,ierr))
   Else
      temperatureInitialTimeStep = (time(size(time))-time(1)) / (size(time) - 1.0_Kr) / 10.0_Kr
      temperatureInitialTime = time(1)
      PetscCallA(MEF90HeatXferCreateTS(MEF90HeatXferCtx,temperatureTS,temperatureResidual,temperatureInitialTime,temperatureInitialTimeStep,ierr))
      PetscCallA(TSGetAdapt(temperatureTS,temperatureTSAdapt,ierr))
   End If

   Select case(MEF90DefMechGlobalOptions%timeSteppingType)
   Case (MEF90DefMech_timeSTeppingTypeQuasiStatic)
      PetscCallA(MEF90DefMechCreateSNESDisplacement(MEF90DefMechCtx,displacementSNES,displacementResidual,ierr))
   End Select
   
   !!! 
   !!! Allocate array of works and energies
   !!!
   Allocate(energy(size(MEF90HeatXferCtx%CellSetOptionsBag)))
   Allocate(cellWork(size(MEF90HeatXferCtx%CellSetOptionsBag)))
   Allocate(faceWork(size(MEF90HeatXferCtx%FaceSetOptionsBag)))

!    !!!
!    !!! Actual computations / time stepping
!    !!!
!    If (MEF90DefMechGlobalOptions%timeSteppingType == MEF90DefMech_timeSteppingTypeQuasiStatic) Then
!       Do step = 1,MEF90GlobalOptions%timeNumStep
!          Write(IOBuffer,100) step,time(step)
!          PetscCallA(PetscPrintf(MEF90Ctx%comm,IOBuffer,ierr))

!          !!! Solve for temperature
!          Select Case (MEF90HeatXferGlobalOptions%timeSteppingType)
!          Case (MEF90HeatXfer_timeSteppingTypeSteadyState)
!             PetscCallA(MEF90HeatXferUpdateTransients(MEF90HeatXferCtx,step,time(step),ierr))
!             PetscCallA(DMLocalToGlobal(temperatureDM,MEF90HeatXferCtx%temperatureLocal,INSERT_VALUES,temperature,ierr))
!             !!! Solve SNES
!             PetscCallA(SNESSolve(temperatureSNES,PETSC_NULL_VEC,temperature,ierr))
!             PetscCallA(DMGlobalToLocal(temperatureDM,temperature,INSERT_VALUES,MEF90HeatXferCtx%temperatureLocal,ierr))
!          Case (MEF90HeatXfer_timeSteppingTypeTransient)
!             If (step > 1) Then
!                Write(IOBuffer,200) step,time(step)
!                PetscCallA(PetscPrintf(MEF90Ctx%comm,IOBuffer,ierr))
!                If (step > 1) Then
!                   !!! Update fields
!                   PetscCallA(MEF90HeatXferUpdateTransients(MEF90HeatXferCtx,step,time(step),ierr))
!                   PetscCallA(TSSetMaxTime(temperatureTS,time(step),ierr))
!                   PetscCallA(DMLocalToGlobal(temperatureDM,MEF90HeatXferCtx%temperatureLocal,INSERT_VALUES,temperature,ierr))
!                   PetscCallA(TSSolve(temperatureTS,temperature,ierr))
!                   PetscCallA(DMGlobalToLocal(temperatureDM,temperature,INSERT_VALUES,MEF90HeatXferCtx%temperatureLocal,ierr))
!                   End If
!                End If
!          End Select

!          !!! Compute energies
!          PetscCallA(MEF90HeatXFerEnergy(MEF90HeatXferCtx,energy,cellWork,faceWork,ierr))
!          PetscCall(DMGetLabelIdIS(temperatureDM,MEF90CellSetLabelName,setIS,ierr))
!          Call MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm,setIS,ierr);CHKERRQ(ierr)
!          PetscCallA(ISGetIndicesF90(setIS,setID,ierr))
!          Do set = 1, size(setID)
!             Write(IOBuffer,101) setID(set),energy(set),cellWork(set),energy(set)-cellWork(set)
!             PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))
!          End Do
!          PetscCallA(ISRestoreIndicesF90(setIS,setID,ierr))
!          PetscCallA(ISDestroy(setIS,ierr))
   
!          PetscCall(DMGetLabelIdIS(temperatureDM,MEF90FaceSetLabelName,setIS,ierr))
!          Call MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm,setIS,ierr);CHKERRQ(ierr)
!          PetscCallA(ISGetIndicesF90(setIS,setID,ierr))
!          Do set = 1, size(setID)
!             Write(IOBuffer,103) setID(set),faceWork(set)
!             PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))
!          End Do
!          PetscCallA(ISRestoreIndicesF90(setIS,setID,ierr))
!          PetscCallA(ISDestroy(setIS,ierr))
   
!          Write(IOBuffer,102) sum(energy),sum(cellWork)+sum(faceWork),sum(energy)-sum(cellWork)-sum(faceWork)
!          PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))
!          !!! Save results
!          PetscCallA(MEF90HeatXferViewEXO(MEF90HeatXferCtx,step,ierr))
   
!          !!! Solve for displacement
!          Select case(MEF90DefMechGlobalOptions%timeSteppingType)
!          Case (MEF90DefMech_timeSteppingTypeQuasiStatic)
!             !PetscCallA(MEF90DefMechUpdateTransients(MEF90DefMechCtx,step,time(step),ierr))
!             PetscCallA(DMLocalToGlobal(displacementDM,MEF90DefMechCtx%displacementLocal,INSERT_VALUES,displacement,ierr))
!             !!! Solve SNES
!             ! PetscCallA(SNESSolve(SNESDßsplacement,PETSC_NULL_VEC,displacement,ierr))
!             ! PetscCallA(DMGlobalToLocal(displacementDM,displacement,INSERT_VALUES,MEF90DefMechCtx%displacementLocal,ierr))

!             !!! Compute energies
!             energy   = 0.0_Kr
!             cellWork = 0.0_Kr
!             faceWork = 0.0_Kr
!             ! PetscCallA(MEF90DefMechWork(MEF90DefMechCtx%displacement,MEF90DefMechCtx,work,ierr))
!             ! PetscCallA(MEF90DefMechElasticEnergy(MEF90DefMechCtx%displacement,MEF90DefMechCtx,energy,ierr))
!             PetscCall(DMGetLabelIdIS(displacementDM,MEF90CellSetLabelName,setIS,ierr))
!             PetscCallA(MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr))
!             PetscCallA(ISGetIndicesF90(setIS,setID,ierr))
!             Do set = 1, size(setID)
!                Write(IOBuffer,201) setID(set),energy(set),cellwork(set),energy(set)-cellwork(set)
!                PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))
!             End Do
!             PetscCallA(ISRestoreIndicesF90(setIS,setID,ierr))
!             PetscCallA(ISDestroy(setIS,ierr))

!             PetscCall(DMGetLabelIdIS(displacementDM,MEF90FaceSetLabelName,setIS,ierr))
!             PetscCallA(MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr))
!             PetscCallA(ISGetIndicesF90(setIS,setID,ierr))
!             Do set = 1, size(setID)
!                Write(IOBuffer,203) setID(set),faceWork(set)
!                PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))
!             End Do
!             PetscCallA(ISRestoreIndicesF90(setIS,setID,ierr))
!             PetscCallA(ISDestroy(setIS,ierr))
!             Write(IOBuffer,202) sum(energy),sum(cellWork)+sum(faceWork),sum(energy)-sum(cellWork)-sum(faceWork)
!             PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))

!             !!! Save results and boundary Values
!             !PetscCallA(MEF90DefMechViewEXO(MEF90DefMechCtx,step,ierr))
!          End Select
!       End Do
!    End If
! 100 Format("\nSolving steady state step ",I4,", t=",ES12.5,"\n")
! 200 Format("\nSolving transient step ",I4,", t=",ES12.5,"\n")
! 101 Format("cell set ",I4," thermal energy: ",ES12.5," flux: ",ES12.5," total: ",ES12.5,"\n")
! 102 Format("======= Total thermal energy: ",ES12.5," flux: ",ES12.5," total: ",ES12.5,"\n")
! 103 Format("face set ",I4,"                              flux: ",ES12.5,"\n")
! 201 Format("cell set ",I4," elastic energy: ",ES12.5," work: ",ES12.5," total: ",ES12.5,"\n")
! 203 Format("face set ",I4,"                              work: ",ES12.5," \n")
! 202 Format("======= Total elastic energy: ",ES12.5," work: ",ES12.5," total: ",ES12.5,"\n")

   !!! Clean up and exit nicely
   Select case(MEF90DefMechGlobalOptions%timeSteppingType)
   Case (MEF90DefMech_timeSTeppingTypeQuasiStatic)
      PetscCallA(SNESDestroy(displacementSNES,ierr))
   End Select

   PetscCallA(VecDestroy(displacementResidual,ierr))
   PetscCallA(VecDestroy(displacement,ierr))

   !!! Clean up and exit nicely
   If (MEF90HeatXferGlobalOptions%timeSteppingType == MEF90HeatXFer_timeSteppingTypeSteadyState) Then
      PetscCallA(SNESDestroy(temperatureSNES,ierr))
   Else
      PetscCallA(TSDestroy(temperatureTS,ierr))
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
