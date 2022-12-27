#include "../MEF90/mef90.inc"
Program vDef
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

   Type(tDM),target                                   :: dm,temperatureDM,displacementDM,damageDM
   Type(tIS)                                          :: setIS
   PetscInt,Dimension(:),Pointer                      :: setID
   PetscInt                                           :: set
   PetscReal,Dimension(:),Pointer                     :: time,elasticEnergy,bodyForceWork,boundaryForceWork,cohesiveEnergy,surfaceEnergy

   Type(tSNES)                                        :: displacementSNES,damageSNES
   SNESConvergedReason                                :: displacementSNESConvergedReason,damageSNESConvergedReason
   Type(tVec)                                         :: displacement,displacementResidual,damage,damageResidual
   Type(tVec)                                         :: damageAltMinOld

   Type(tSNES)                                        :: temperatureSNES
   SNESConvergedReason                                :: temperatureSNESConvergedReason
   Type(tTS)                                          :: temperatureTS
   Type(tTSAdapt)                                     :: temperatureTSAdapt
   Type(tVec)                                         :: temperature,temperatureResidual

   PetscReal                                          :: temperatureInitialTimeStep,temperatureInitialTime
   !PetscInt                                           :: tsTemperatureMaxIter
   PetscLogStage                                      :: logStageHeatXfer,logStageDamage,logStageDisplacement,logStageEnergy,logStageIO


   PetscBool                                          :: flg
   Character(len=MEF90MXSTRLEN)                       :: IOBuffer
   Type(tPetscViewer)                                 :: logViewer

   PetscInt                                           :: step
   PetscInt                                           :: AltMinIter,AltMinStep=0_Ki
   PetscReal                                          :: damageMaxChange,damageMin,damageMax
   
   
   PetscInt                                           :: dim

   !!! Initialize MEF90
   PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
   PetscCallA(MEF90Initialize(ierr))
   PetscCallA(PetscLogStageRegister('HeatXfer    ',logStageHeatXfer,ierr))
   PetscCallA(PetscLogStageRegister('Damage      ',logStageDamage,ierr))
   PetscCallA(PetscLogStageRegister('Displacement',logStageDisplacement,ierr))
   PetscCallA(PetscLogStageRegister('Energy      ',logStageEnergy,ierr))
   PetscCallA(PetscLogStageRegister('IO          ',logStageIO,ierr))


   !!! Get all MEF90-wide options
   PetscCallA(MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90CtxDefaultGlobalOptions,ierr))
   PetscCallA(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr))

   PetscCallA(PetscLogStagePush(logStageIO,ierr))
   If (MEF90GlobalOptions%verbose > 1) Then
      PetscCallA(PetscPrintf(PETSC_COMM_WORLD,"Reading geometry\n",ierr))
   End If
   PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm,MEF90Ctx%geometryFile,PETSC_NULL_CHARACTER,PETSC_TRUE,dm,ierr))
   PetscCallA(DMPlexDistributeSetDefault(dm,PETSC_FALSE,ierr))
   PetscCallA(DMSetUseNatural(dm,PETSC_TRUE,ierr))
   PetscCallA(DMSetFromOptions(dm,ierr))
   PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-mef90_dm_view",ierr))

   PetscCallA(DMGetDimension(dm,dim,ierr))
   PetscCallA(MEF90CtxGetTime(MEF90Ctx,time,ierr))

   !!! Create HeatXfer context, get all HeatXfer options
   PetscCallA(MEF90HeatXferCtxCreate(MEF90HeatXferCtx,dm,MEF90Ctx,ierr))
   PetscCallA(MEF90HeatXferCtxSetFromOptions(MEF90HeatXferCtx,PETSC_NULL_CHARACTER,HeatXferDefaultGlobalOptions,HeatXferDefaultCellSetOptions,HeatXferDefaultFaceSetOptions,HeatXferDefaultVertexSetOptions,ierr))
   PetscCallA(PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr))

   !!! Create DefMechCtx, get all defMech options
   PetscCallA(MEF90DefMechCtxCreate(MEF90DefMechCtx,dm,MEF90Ctx,ierr))
   PetscCallA(MEF90DefMechCtxSetFromOptions(MEF90DefMechCtx,PETSC_NULL_CHARACTER,DefMechDefaultGlobalOptions,DefMechDefaultCellSetOptions,DefMechDefaultFaceSetOptions,DefMechDefaultVertexSetOptions,ierr))
   PetscCallA(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr))

   !!! Calling Inquire on all MPI ranks followed by exopen_par (MEF90CtxOpenEXO) can lead to a strange race condition
   !!! Strangely enough, adding an MPI_Barrier does not help.
   !!! There is no real good reason to call Inquire on all ranks anyway.
   If (MEF90Ctx%rank == 0) Then
      Inquire(file=MEF90Ctx%resultFile,exist=flg)
   End If
   PetscCallMPIA(MPI_Bcast(flg,1,MPI_LOGICAL,0,MEF90Ctx%Comm,ierr))
   If (flg) Then
      ! we assume that the output file is formatted
      If (MEF90GlobalOptions%verbose > 1) Then
         PetscCallA(PetscPrintf(PETSC_COMM_WORLD,"Opening result file\n",ierr))
      End If
      PetscCallA(MEF90CtxOpenEXO(MEF90Ctx,MEF90Ctx%resultViewer,FILE_MODE_APPEND,ierr))
   Else
      ! we need to create the output file
      If (MEF90GlobalOptions%verbose > 1) Then
         PetscCallA(PetscPrintf(PETSC_COMM_WORLD,"Creating result file\n",ierr))
      End If
      PetscCallA(PetscViewerDestroy(MEF90Ctx%resultViewer,ierr))
      PetscCallA(MEF90CtxOpenEXO(MEF90Ctx,MEF90Ctx%resultViewer,FILE_MODE_WRITE,ierr))
      PetscCallA(MEF90EXODMView(dm,MEF90Ctx%resultViewer,MEF90GlobalOptions%elementOrder,ierr))

      If (MEF90GlobalOptions%verbose > 1) Then
         PetscCallA(PetscPrintf(PETSC_COMM_WORLD,"Formatting result file\n",ierr))
      End If
      PetscCallA(MEF90DefMechFormatEXO(MEF90DefMechCtx,time,ierr))
      If (MEF90GlobalOptions%verbose > 1) Then
         PetscCallA(PetscPrintf(PETSC_COMM_WORLD,"Done Formatting result file\n",ierr))
      End If
   End If
   PetscCallA(PetscLogStagePop(ierr))

   distribute: Block 
      Type(tDM),target                    :: dmDist
      PetscInt                            :: ovlp = 0_Ki
      Type(tPetscSF)                      :: naturalPointSF

      If (MEF90Ctx%NumProcs > 1) Then
         If (MEF90GlobalOptions%verbose > 1) Then
            PetscCallA(PetscPrintf(PETSC_COMM_WORLD,"Distributing mesh\n",ierr))
         End If
         PetscCallA(DMSetUseNatural(dm,PETSC_TRUE,ierr))
         PetscCallA(DMPlexDistribute(dm,ovlp,naturalPointSF,dmDist,ierr))
         PetscCallA(DMPlexSetMigrationSF(dmDist,naturalPointSF, ierr))
         PetscCallA(PetscSFDestroy(naturalPointSF,ierr))
         PetscCallA(DMDestroy(dm,ierr))
         dm = dmDist
      End If
   End Block distribute
   PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-mef90_dm_view",ierr))

   !!! We no longer need the DM. We have the megaDM in MEF90HeatXferCtx and MEF90DefMechCtx
   PetscCallA(DMDestroy(dm,ierr))

   !!! Get parse all materials data from the command line
   If (dim == 2) Then
      PetscCallA(MEF90MatPropBagSetFromOptions(MEF90DefMechCtx%MaterialPropertiesBag,MEF90DefMechCtx%megaDM,MEF90Mathium2D,MEF90Ctx,ierr))
   Else
      PetscCallA(MEF90MatPropBagSetFromOptions(MEF90DefMechCtx%MaterialPropertiesBag,MEF90DefMechCtx%megaDM,MEF90Mathium3D,MEF90Ctx,ierr))
   End If   
   MEF90HeatXferCtx%MaterialPropertiesBag => MEF90DefMechCtx%MaterialPropertiesBag

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

   PetscCallA(VecGetDM(MEF90DefMechCtx%damageLocal,damageDM,ierr)) 
   !!! This only borrows a reference so we do not need to delete it
   PetscCallA(DMCreateGlobalVector(damageDM,damage,ierr))
   PetscCallA(PetscObjectSetName(damage,"damage",ierr))
   PetscCallA(VecDuplicate(damage,damageResidual,ierr))
   PetscCallA(VecDuplicate(damage,damageAltMinOld,ierr))

   !!! 
   !!! Create SNES or TS, Mat and set KSP default options
   !!!
   Select Case (MEF90HeatXferGlobalOptions%timeSteppingType)
   Case (MEF90HeatXFer_timeSteppingTypeSteadyState) 
      PetscCallA(MEF90HeatXferCreateSNES(MEF90HeatXferCtx,temperatureSNES,temperatureResidual,ierr))
   Case (MEF90HeatXFer_timeSteppingTypeTransient)
      temperatureInitialTimeStep = (time(size(time))-time(1)) / (size(time) - 1.0_Kr) / 10.0_Kr
      temperatureInitialTime = time(1)
      PetscCallA(MEF90HeatXferCreateTS(MEF90HeatXferCtx,temperatureTS,temperatureResidual,temperatureInitialTime,temperatureInitialTimeStep,ierr))
      PetscCallA(TSGetAdapt(temperatureTS,temperatureTSAdapt,ierr))
   Case(MEF90HeatXfer_timeSteppingTypeNULL)
      Continue
   End Select

   Select Case(MEF90DefMechGlobalOptions%timeSteppingType)
   Case (MEF90DefMech_TimeSteppingTypeQuasiStatic)
      PetscCallA(MEF90DefMechCreateSNESDisplacement(MEF90DefMechCtx,displacementSNES,displacementResidual,ierr))
      PetscCallA(MEF90DefMechCreateSNESDamage(MEF90DefMechCtx,damageSNES,damageResidual,ierr))
   Case (MEF90DefMech_TimeSteppingTypeNULL)
      Continue
   End Select
   !!! 
   !!! Allocate array of works and energies
   !!!
   Allocate(elasticEnergy(size(MEF90HeatXferCtx%CellSetOptionsBag)))
   Allocate(bodyForceWork(size(MEF90HeatXferCtx%CellSetOptionsBag)))
   Allocate(cohesiveEnergy(size(MEF90HeatXferCtx%CellSetOptionsBag)))
   Allocate(surfaceEnergy(size(MEF90HeatXferCtx%CellSetOptionsBag)))
   Allocate(boundaryForceWork(size(MEF90HeatXferCtx%FaceSetOptionsBag)))

   !!!
   !!! Actual computations / time stepping
   !!!
   If ((MEF90DefMechGlobalOptions%timeSteppingType /= MEF90DefMech_TimeSteppingTypeNULL) .OR. (MEF90HeatXferGlobalOptions%timeSteppingType /= MEF90DefMech_TimeSteppingTypeNULL))  Then
      Do step = 1,MEF90GlobalOptions%timeNumStep
         Write(IOBuffer,100) step,time(step)
         PetscCallA(PetscPrintf(MEF90Ctx%comm,IOBuffer,ierr))

         !!! Solve for temperature
         PetscCallA(PetscLogStagePush(logStageHeatXfer,ierr))
         Select Case (MEF90HeatXferGlobalOptions%timeSteppingType)
         Case (MEF90HeatXfer_timeSteppingTypeSteadyState)
            PetscCallA(MEF90HeatXferSetTransients(MEF90HeatXferCtx,step,time(step),ierr))
            PetscCallA(DMLocalToGlobal(temperatureDM,MEF90HeatXferCtx%temperatureLocal,INSERT_VALUES,temperature,ierr))
            !!! Solve SNES
            PetscCallA(SNESSolve(temperatureSNES,PETSC_NULL_VEC,temperature,ierr))
            PetscCallA(SNESGetConvergedReason(temperatureSNES,temperatureSNESConvergedReason,ierr))
            If (temperatureSNESConvergedReason < 0) Then  
               Write(IOBuffer,400) "temperature",temperatureSNESConvergedReason
               PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))
            End If

            PetscCallA(DMGlobalToLocal(temperatureDM,temperature,INSERT_VALUES,MEF90HeatXferCtx%temperatureLocal,ierr))
            PetscCallA(VecCopy(MEF90HeatXferCtx%temperatureLocal,MEF90DefMechCtx%temperatureLocal,ierr))
         Case (MEF90HeatXfer_timeSteppingTypeTransient)
            If (step > 1) Then
               Write(IOBuffer,200) step,time(step)
               PetscCallA(PetscPrintf(MEF90Ctx%comm,IOBuffer,ierr))
               If (step > 1) Then
                  !!! Update fields
                  PetscCallA(MEF90HeatXferSetTransients(MEF90HeatXferCtx,step,time(step),ierr))
                  PetscCallA(TSSetMaxTime(temperatureTS,time(step),ierr))
                  PetscCallA(DMLocalToGlobal(temperatureDM,MEF90HeatXferCtx%temperatureLocal,INSERT_VALUES,temperature,ierr))
                  PetscCallA(TSSolve(temperatureTS,temperature,ierr))
                  PetscCallA(DMGlobalToLocal(temperatureDM,temperature,INSERT_VALUES,MEF90HeatXferCtx%temperatureLocal,ierr))
                  End If
               End If
         End Select

         If (MEF90HeatXferGlobalOptions%timeSteppingType /= MEF90DefMech_TimeSteppingTypeNULL) Then
            !!! Compute energies
            PetscCallA(MEF90HeatXFerEnergy(MEF90HeatXferCtx,elasticEnergy,bodyForceWork,boundaryForceWork,ierr))
            PetscCallA(DMGetLabelIdIS(temperatureDM,MEF90CellSetLabelName,setIS,ierr))
            PetscCallA(MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm,setIS,ierr))
            PetscCallA(ISGetIndicesF90(setIS,setID,ierr))
            Do set = 1, size(setID)
               Write(IOBuffer,101) setID(set),elasticEnergy(set),bodyForceWork(set),elasticEnergy(set)-bodyForceWork(set)
               PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))
            End Do
            PetscCallA(ISRestoreIndicesF90(setIS,setID,ierr))
            PetscCallA(ISDestroy(setIS,ierr))
      
            PetscCallA(DMGetLabelIdIS(temperatureDM,MEF90FaceSetLabelName,setIS,ierr))
            PetscCallA(MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm,setIS,ierr))
            PetscCallA(ISGetIndicesF90(setIS,setID,ierr))
            Do set = 1, size(setID)
               Write(IOBuffer,103) setID(set),boundaryForceWork(set)
               PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))
            End Do
            PetscCallA(ISRestoreIndicesF90(setIS,setID,ierr))
            PetscCallA(ISDestroy(setIS,ierr))
      
            Write(IOBuffer,102) sum(elasticEnergy),sum(bodyForceWork)+sum(boundaryForceWork),sum(elasticEnergy)-sum(bodyForceWork)-sum(boundaryForceWork)
            PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))
            !!! Save results
            PetscCallA(MEF90HeatXferViewEXO(MEF90HeatXferCtx,step,ierr))
         End If
         PetscCallA(PetscLogStagePop(ierr))

         !!! Solve for displacement and damage
         PetscCallA(MEF90DefMechSetTransients(MEF90DefMechCtx,step,time(step),ierr))
         PetscCallA(MEF90DefMechUpdateDamageBounds(MEF90DefMechCtx,damageSNES,damage,ierr))
         PetscCallA(DMLocalToGlobal(displacementDM,MEF90DefMechCtx%displacementLocal,INSERT_VALUES,displacement,ierr))
         !PetscCallA(VecCopy(MEF90DefMechCtx%displacement,MEF90DefMechCtx%displacementPreviousStep,ierr))
         !PetscCallA(VecCopy(damage,damagePreviousStep,ierr))

         Select case(MEF90DefMechGlobalOptions%timeSteppingType)
         Case (MEF90DefMech_timeSteppingTypeQuasiStatic)
            Select case(MEF90DefMechGlobalOptions%SolverType)
            Case(MEF90DefMech_SolverTypeAltMin)
               !PetscCallA(SNESSetLagPreconditioner(displacementSNES,1,ierr))
               PetscCallA(SNESSetLagPreconditioner(damageSNES,1_Ki,ierr))

               AltMin: Do AltMinIter = 1, MEF90DefMechGlobalOptions%maxit
                  AltMinStep = AltMinStep + 1
                  Write(IObuffer,208) AltMinIter
                  PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))

                  If (mod(AltMinIter-1,MEF90DefMechGlobalOptions%PCLag) == 0) Then
                     !PetscCallA(SNESSetLagPreconditioner(displacementSNES,-2_Ki,ierr))
                     PetscCallA(SNESSetLagPreconditioner(damageSNES,-2_Ki,ierr))
                  End If 

                  !!! Solve SNES displacement
                  PetscCallA(PetscLogStagePush(logStageDisplacement,ierr))
                  PetscCallA(SNESSolve(displacementSNES,PETSC_NULL_VEC,displacement,ierr))
                  PetscCallA(SNESGetConvergedReason(displacementSNES,displacementSNESConvergedReason,ierr))
                  If (displacementSNESConvergedReason < 0) Then  
                     Write(IOBuffer,400) "displacement",displacementSNESConvergedReason
                     PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))
                  End If

                  PetscCallA(DMGlobalToLocal(displacementDM,displacement,INSERT_VALUES,MEF90DefMechCtx%displacementLocal,ierr))
                  PetscCallA(DMLocalToGlobal(damageDM,MEF90DefMechCtx%damageLocal,INSERT_VALUES,damage,ierr))
                  PetscCallA(PetscLogStagePop(ierr))

                  PetscCallA(VecCopy(damage,damageAltMinOld,ierr))
                  !!! Solve SNES damage
                  PetscCallA(PetscLogStagePush(logStageDamage,ierr))
                  PetscCallA(SNESSolve(damageSNES,PETSC_NULL_VEC,damage,ierr))
                  PetscCallA(SNESGetConvergedReason(damageSNES,damageSNESConvergedReason,ierr))
                  If (damageSNESConvergedReason < 0) Then  
                     Write(IOBuffer,400) "damage",damageSNESConvergedReason
                     PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))
                  End If
                  PetscCallA(DMGlobalToLocal(damageDM,damage,INSERT_VALUES,MEF90DefMechCtx%damageLocal,ierr))
                  PetscCallA(PetscLogStagePop(ierr))

                  PetscCallA(VecMin(damage,PETSC_NULL_INTEGER,damageMin,ierr))
                  PetscCallA(VecMax(damage,PETSC_NULL_INTEGER,damageMax,ierr))
                  
                  PetscCallA(VecAxPy(damageAltMinOld,-1.0_Kr,damage,ierr))
                  ! damageAltMinOld = damageAltMinOld - damage
                  PetscCallA(VecNorm(damageAltMinOld,NORM_INFINITY,damageMaxChange,ierr))
                  Write(IOBuffer,209) damageMin,damageMax,damageMaxChange                  
                  PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))

                  If (damageMaxChange <= MEF90DefMechGlobalOptions%damageATol) Then
                     EXIT
                  End If

                  If (mod(AltMinIter,25) == 0) Then
                     !!! Save results and boundary Values
                     PetscCallA(PetscLogStagePush(logStageIO,ierr))
                     PetscCallA(MEF90DefMechViewEXO(MEF90DefMechCtx,step,ierr))
                     PetscCallA(PetscLogStagePop(ierr))
                  End If
               End Do AltMin
            Case default
               Write(IOBuffer,*) "Unimplemented DefMech solver type: ", MEF90DefMechGlobalOptions%SolverType, "\n"
               PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))
               STOP
            End Select ! solverType
            !!! Compute energies
            PetscCallA(PetscLogStagePush(logStageEnergy,ierr))
            elasticEnergy     = 0.0_Kr
            bodyForceWork     = 0.0_Kr
            boundaryForceWork = 0.0_Kr
            surfaceEnergy     = 0.0_Kr
            cohesiveEnergy    = 0.0_Kr
            PetscCallA(MEF90DefMechWork(MEF90DefMechCtx,bodyForceWork,boundaryForceWork,ierr))
            PetscCallA(MEF90DefMechElasticEnergy(MEF90DefMechCtx,elasticEnergy,ierr))
            PetscCallA(MEF90DefMechSurfaceEnergy(MEF90DefMechCtx,surfaceEnergy,ierr))

            PetscCallA(MEF90DefMechStress(MEF90DefMechCtx,MEF90DefMechCtx%stress,ierr))
            PetscCallA(DMGetLabelIdIS(displacementDM,MEF90CellSetLabelName,setIS,ierr))
            PetscCallA(MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr))
            PetscCallA(ISGetIndicesF90(setIS,setID,ierr))
            Do set = 1, size(setID)
               Write(IOBuffer,201) setID(set),elasticEnergy(set),bodyForceWork(set),cohesiveEnergy(set),surfaceEnergy(set),elasticEnergy(set)-bodyForceWork(set)+cohesiveEnergy(set)+surfaceEnergy(set)
               PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))
            End Do
            PetscCallA(ISRestoreIndicesF90(setIS,setID,ierr))
            PetscCallA(ISDestroy(setIS,ierr))

            PetscCallA(DMGetLabelIdIS(displacementDM,MEF90FaceSetLabelName,setIS,ierr))
            PetscCallA(MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr))
            PetscCallA(ISGetIndicesF90(setIS,setID,ierr))
            Do set = 1, size(setID)
               Write(IOBuffer,203) setID(set),boundaryForceWork(set)
               PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))
            End Do
            PetscCallA(ISRestoreIndicesF90(setIS,setID,ierr))
            PetscCallA(ISDestroy(setIS,ierr))
            Write(IOBuffer,202) sum(elasticEnergy),sum(bodyForceWork)+sum(boundaryForceWork),sum(cohesiveEnergy),sum(surfaceEnergy),sum(elasticEnergy)-sum(bodyForceWork)-sum(boundaryForceWork)+sum(cohesiveEnergy)+sum(surfaceEnergy)
            PetscCallA(PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr))
            PetscCallA(PetscLogStagePop(ierr))

            !!! Save results and boundary Values
            PetscCallA(PetscLogStagePush(logStageIO,ierr))
            PetscCallA(MEF90DefMechViewEXO(MEF90DefMechCtx,step,ierr))
            PetscCallA(PetscLogStagePop(ierr))

            PetscCallA(PetscViewerASCIIOpen(MEF90Ctx%comm,trim(MEF90FilePrefix(MEF90Ctx%resultFile))//'.log',logViewer, ierr))
            PetscCallA(PetscLogView(logViewer,ierr))
            PetscCallA(PetscViewerDestroy(logViewer,ierr)) 
         End Select ! timeStepingType
      End Do ! step
   End If ! timeSteppingType
   Write(IOBuffer,*) 'Total number of alternate minimizations:',AltMinStep,'\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)         
100 Format("\nSolving steady state step ",I4,", t=",ES12.5,"\n")
200 Format("\nSolving transient step ",I4,", t=",ES12.5,"\n")
101 Format("cell set ",I4," thermal energy: ",ES12.5," flux: ",ES12.5," total: ",ES12.5,"\n")
102 Format("======= Total thermal energy: ",ES12.5," flux: ",ES12.5," total: ",ES12.5,"\n")
103 Format("face set ",I4,"                              flux: ",ES12.5,"\n")
201 Format("cell set ",I4,"  elastic energy: ",ES12.5," work: ",ES12.5," cohesive: ",ES12.5," surface: ",ES12.5," total: ",ES12.5,"\n")
203 Format("face set ",I4,"                      boundary work: ",ES12.5,"\n")
202 Format("======= Total: elastic energy: ",ES12.5," work: ",ES12.5," cohesive: ",ES12.5," surface: ",ES12.5," total: ",ES12.5,"\n")

208 Format("   Alt. Min. step ",I5," ")
209 Format(" alpha min / max", ES12.5, " / ", ES12.5, ", max change ", ES12.5,"\n")

400 Format(" [ERROR]: ",A," SNESSolve failed with SNESConvergedReason ",I2,". \n Check https://petsc.org/release/docs/manualpages/SNES/SNESConvergedReason/ for error code meaning.\n")

   !!! Clean up and exit nicely
   Select case(MEF90DefMechGlobalOptions%timeSteppingType)
   Case (MEF90DefMech_timeSTeppingTypeQuasiStatic)
      PetscCallA(SNESDestroy(displacementSNES,ierr))
      PetscCallA(VecDestroy(displacementResidual,ierr))
      PetscCallA(VecDestroy(displacement,ierr))
   End Select

   !!! Clean up and exit nicely
   Select Case (MEF90HeatXferGlobalOptions%timeSteppingType)
   Case (MEF90HeatXFer_timeSteppingTypeSteadyState) 
      PetscCallA(SNESDestroy(temperatureSNES,ierr))
   Case (MEF90HeatXFer_timeSteppingTypeTransient)
      PetscCallA(TSDestroy(temperatureTS,ierr))
   End Select
   PetscCallA(VecDestroy(temperatureResidual,ierr))
   PetscCallA(VecDestroy(temperature,ierr))

   Select Case(MEF90DefMechGlobalOptions%timeSteppingType)
   Case (MEF90DefMech_TimeSteppingTypeQuasiStatic)
      PetscCallA(SNESDestroy(displacementSNES,ierr))
      PetscCallA(SNESDestroy(damageSNES,ierr))
   End Select
   PetscCallA(VecDestroy(displacementResidual,ierr))
   PetscCallA(VecDestroy(displacement,ierr))
   PetscCallA(VecDestroy(damageResidual,ierr))
   PetscCallA(VecDestroy(damageAltMinOld,ierr))
   PetscCallA(VecDestroy(damage,ierr))

   DeAllocate(time)
   DeAllocate(elasticEnergy)
   DeAllocate(bodyForceWork)
   DeAllocate(boundaryForceWork)
   DeAllocate(cohesiveEnergy)
   DeAllocate(surfaceEnergy)
   PetscCallA(MEF90DefMechCtxDestroy(MEF90DefMechCtx,ierr))
   PetscCallA(MEF90HeatXferCtxDestroy(MEF90HeatXferCtx,ierr))
   
   PetscCallA(PetscViewerDestroy(MEF90Ctx%resultViewer,ierr))
   If (.NOT. MEF90GlobalOptions%dryrun) Then
      Call PetscViewerASCIIOpen(MEF90Ctx%comm,trim(MEF90FilePrefix(MEF90Ctx%resultFile))//'.log',logViewer, ierr);CHKERRQ(ierr)
      Call PetscLogView(logViewer,ierr);CHKERRQ(ierr)
      Call PetscViewerDestroy(logViewer,ierr);CHKERRQ(ierr)
   End If
   PetscCallA(MEF90CtxDestroy(MEF90Ctx,ierr))
   PetscCallA(MEF90Finalize(ierr))
   PetscCallA(PetscFinalize(ierr))
End Program vDef
