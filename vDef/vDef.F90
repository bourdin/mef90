#include "../MEF90/mef90.inc"
program vDef
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petsctao.h"
   use m_MEF90
   use m_MEF90_DefMechCtx
   use m_MEF90_DefMech
   use m_MEF90_HeatXferCtx
   use m_MEF90_HeatXfer
   use m_vDefDefault
   use petsc
   use petsctao
   implicit none(type)

   PetscErrorCode                                     :: ierr
   type(MEF90Ctx_Type), target                        :: MEF90Ctx
   type(MEF90CtxGlobalOptions_Type), pointer          :: MEF90GlobalOptions

   !!! Defect mechanics contexts
   type(MEF90DefMechCtx_Type)                         :: MEF90DefMechCtx
   type(MEF90DefMechGlobalOptions_Type), pointer      :: MEF90DefMechGlobalOptions
   !!! HeatXfer contexts
   type(MEF90HeatXferCtx_Type)                        :: MEF90HeatXferCtx
   type(MEF90HeatXferGlobalOptions_Type), pointer     :: MEF90HeatXferGlobalOptions

   type(tDM), target                                  :: dm, temperatureDM, displacementDM, damageDM
   type(tIS)                                          :: setIS
   PetscInt, dimension(:), pointer                    :: setID
   PetscInt                                           :: set
   PetscReal, dimension(:), pointer                   :: time, elasticEnergy, bodyForceWork, boundaryForceWork, cohesiveEnergy, surfaceEnergy

   type(tSNES)                                        :: displacementSNES, damageSNES
   type(tTao)                                         :: damageTAO
   type(eSNESConvergedReason)                         :: displacementSNESConvergedReason, damageSNESConvergedReason
   type(eTaoConvergedReason)                          :: damageTAOConvergedReason
   type(tVec)                                         :: displacement, displacementResidual, damage, damageResidual
   type(tVec)                                         :: damageAltMinOld
   type(tVec)                                         :: damageLB, damageUB
   PetscReal, dimension(:), pointer                   :: damageArray, damageAltMinOldArray, damageLBArray, damageUBArray
   PetscInt                                           :: iDof
   PetscReal                                          :: SOROmega, mySOROmega

   type(tSNES)                                        :: temperatureSNES
   type(eSNESConvergedReason)                         :: temperatureSNESConvergedReason
   type(tTS)                                          :: temperatureTS
   type(tTSAdapt)                                     :: temperatureTSAdapt
   type(tVec)                                         :: temperature, temperatureResidual

   PetscReal                                          :: temperatureInitialTimeStep, temperatureInitialTime
   !PetscInt                                           :: tsTemperatureMaxIter
   PetscLogStage                                      :: logStageHeatXfer, logStageDamage, logStageDisplacement, logStageEnergy, logStageIO

   PetscBool                                          :: flg, EXONeedsFormatting = PETSC_FALSE
   character(len=MEF90MXSTRLEN)                       :: IOBuffer
   Type(tPetscViewer)                                 :: logViewer

   PetscInt                                           :: step
   PetscExodusIIInt                                   :: EXOstep
   PetscInt                                           :: AltMinIter, AltMinStep = 0_ki
   PetscReal                                          :: damageMaxChange, damageMin, damageMax

   !!! Initialize MEF90
   PetscCallA(PetscInitialize(ierr))
   PetscCallA(MEF90Initialize(PETSC_COMM_WORLD, ierr))
   PetscCallA(PetscLogStageRegister('HeatXfer    ', logStageHeatXfer, ierr))
   PetscCallA(PetscLogStageRegister('Damage      ', logStageDamage, ierr))
   PetscCallA(PetscLogStageRegister('Displacement', logStageDisplacement, ierr))
   PetscCallA(PetscLogStageRegister('Energy      ', logStageEnergy, ierr))
   PetscCallA(PetscLogStageRegister('IO          ', logStageIO, ierr))

   !!! Get all MEF90-wide options
   PetscCallA(MEF90CtxCreate(PETSC_COMM_WORLD, MEF90Ctx, MEF90CtxDefaultGlobalOptions, ierr))
   PetscCallA(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag, MEF90GlobalOptions, ierr))

   if (MEF90GlobalOptions%verbose > 0) then
      PetscCallA(PetscPrintf(MEF90Ctx%comm, "Reading geometry\n", ierr))
   end if
   PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm, MEF90Ctx%geometryFile, PETSC_NULL_CHARACTER, PETSC_TRUE, dm, ierr))
   PetscCallA(DMPlexDistributeSetDefault(dm, PETSC_FALSE, ierr))
   PetscCallA(DMSetUseNatural(dm, PETSC_TRUE, ierr))
   PetscCallA(DMSetFromOptions(dm, ierr))
   PetscCallA(DMViewFromOptions(dm, PETSC_NULL_OBJECT, "-mef90_dm_view", ierr))

   !!! Calling Inquire on all MPI ranks followed by exopen_par (MEF90CtxOpenEXO) can lead to a strange race condition
   !!! Strangely enough, adding an MPI_Barrier does not help.
   !!! There is no real good reason to call Inquire on all ranks anyway.
   if (MEF90Ctx%rank == 0) then
      inquire (file=MEF90Ctx%resultFile, exist=flg)
   end if
   PetscCallMPIA(MPI_Bcast(flg, 1, MPI_LOGICAL, 0, MEF90Ctx%Comm, ierr))
   if (flg) then
      ! we assume that the output file is formatted
      if (MEF90GlobalOptions%verbose > 0) then
         PetscCallA(PetscPrintf(MEF90Ctx%comm, "Opening result file\n", ierr))
      end if
      PetscCallA(MEF90CtxOpenEXO(MEF90Ctx, MEF90Ctx%resultViewer, FILE_MODE_APPEND, ierr))
   else
      ! we need to create the output file
      if (MEF90GlobalOptions%verbose > 0) then
         PetscCallA(PetscPrintf(MEF90Ctx%comm, "Creating result file\n", ierr))
      end if
      ! PetscCallA(PetscViewerDestroy(MEF90Ctx%resultViewer,ierr))
      PetscCallA(MEF90CtxOpenEXO(MEF90Ctx, MEF90Ctx%resultViewer, FILE_MODE_WRITE, ierr))
      PetscCallA(MEF90EXODMView(dm, MEF90Ctx%resultViewer, MEF90GlobalOptions%elementOrder, ierr))
      EXONeedsFormatting = PETSC_TRUE
   end if

   distribute: block
      type(tDM), target                   :: dmDist
      PetscInt                            :: ovlp = 0_ki
      type(tPetscSF)                      :: naturalPointSF

      if (MEF90Ctx%NumProcs > 1) then
         if (MEF90GlobalOptions%verbose > 0) then
            PetscCallA(PetscPrintf(MEF90Ctx%comm, "Distributing mesh\n", ierr))
         end if
         PetscCallA(DMSetUseNatural(dm, PETSC_TRUE, ierr))
         PetscCallA(DMPlexDistribute(dm, ovlp, naturalPointSF, dmDist, ierr))
         PetscCallA(DMPlexSetMigrationSF(dmDist, naturalPointSF, ierr))
         PetscCallA(PetscSFDestroy(naturalPointSF, ierr))
         PetscCallA(DMDestroy(dm, ierr))
         dm = dmDist
      end if
   end block distribute
   PetscCallA(DMViewFromOptions(dm, PETSC_NULL_OBJECT, "-mef90_dm_view", ierr))

   !!! Create HeatXfer context, get all HeatXfer options
   PetscCallA(MEF90HeatXferCtxCreate(MEF90HeatXferCtx, dm, MEF90Ctx, ierr))
   PetscCallA(MEF90HeatXferCtxSetFromOptions(MEF90HeatXferCtx, PETSC_NULL_CHARACTER, HeatXferDefaultGlobalOptions, HeatXferDefaultCellSetOptions, HeatXferDefaultFaceSetOptions, HeatXferDefaultVertexSetOptions, ierr))
   PetscCallA(PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag, MEF90HeatXferGlobalOptions, ierr))

   !!! Create DefMechCtx, get all defMech options
   PetscCallA(MEF90DefMechCtxCreate(MEF90DefMechCtx, dm, MEF90Ctx, ierr))
   PetscCallA(MEF90DefMechCtxSetFromOptions(MEF90DefMechCtx, PETSC_NULL_CHARACTER, DefMechDefaultGlobalOptions, DefMechDefaultCellSetOptions, DefMechDefaultFaceSetOptions, DefMechDefaultVertexSetOptions, ierr))
   PetscCallA(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag, MEF90DefMechGlobalOptions, ierr))
   PetscCallA(VecDestroy(MEF90DefMechCtx%temperatureLocal, ierr))
   deallocate (MEF90DefMechCtx%temperatureLocal)
   MEF90DefMechCtx%temperatureLocal => MEF90HeatXferCtx%temperatureLocal

   !!! We no longer need the DM. We have the megaDM in MEF90HeatXferCtx and MEF90DefMechCtx
   PetscCallA(DMDestroy(dm, ierr))

   !!! Get parse all materials data from the command line
   if (MEF90DefMechCtx%dim == 2) then
      PetscCallA(MEF90MatPropBagSetFromOptions(MEF90DefMechCtx%MaterialPropertiesBag, MEF90DefMechCtx%megaDM, MEF90Mathium2D, MEF90Ctx, ierr))
   else
      PetscCallA(MEF90MatPropBagSetFromOptions(MEF90DefMechCtx%MaterialPropertiesBag, MEF90DefMechCtx%megaDM, MEF90Mathium3D, MEF90Ctx, ierr))
   end if
   MEF90HeatXferCtx%MaterialPropertiesBag => MEF90DefMechCtx%MaterialPropertiesBag

   !!! Create GLOBAL vectors for the unknowns (temperature,displacements), residuals, etc
   PetscCallA(VecGetDM(MEF90HeatXferCtx%temperatureLocal, temperatureDM, ierr))
   !!! This only borrows a reference so we do not need to delete it
   PetscCallA(DMCreateGlobalVector(temperatureDM, temperature, ierr))
   PetscCallA(PetscObjectSetName(temperature, "Temperature", ierr))
   PetscCallA(VecDuplicate(temperature, temperatureResidual, ierr))
   PetscCallA(PetscObjectSetName(temperatureResidual, "temperatureResidual", ierr))

   PetscCallA(VecGetDM(MEF90DefMechCtx%displacementLocal, displacementDM, ierr))
   !!! This only borrows a reference so we do not need to delete it
   PetscCallA(DMCreateGlobalVector(displacementDM, displacement, ierr))
   PetscCallA(PetscObjectSetName(displacement, "displacement", ierr))
   PetscCallA(VecDuplicate(displacement, displacementResidual, ierr))
   PetscCallA(PetscObjectSetName(displacementResidual, "displacementResidual", ierr))

   PetscCallA(VecGetDM(MEF90DefMechCtx%damageLocal, damageDM, ierr))
   !!! This only borrows a reference so we do not need to delete it
   PetscCallA(DMCreateGlobalVector(damageDM, damage, ierr))
   PetscCallA(PetscObjectSetName(damage, "damage", ierr))
   PetscCallA(VecDuplicate(damage, damageResidual, ierr))
   PetscCallA(VecDuplicate(damage, damageAltMinOld, ierr))

   !!!
   !!! Create SNES or TS, Mat and set KSP default options
   !!!
   select case (MEF90HeatXferGlobalOptions%timeSteppingType)
   case (MEF90HeatXFer_timeSteppingTypeSteadyState)
      PetscCallA(MEF90HeatXferCreateSNES(MEF90HeatXferCtx, temperatureSNES, temperatureResidual, ierr))
   case (MEF90HeatXFer_timeSteppingTypeTransient)
      temperatureInitialTimeStep = (time(size(time)) - time(1)) / (size(time) - 1.0_kr) / 10.0_kr
      temperatureInitialTime = time(1)
      PetscCallA(MEF90HeatXferCreateTS(MEF90HeatXferCtx, temperatureTS, temperatureResidual, temperatureInitialTime, temperatureInitialTimeStep, ierr))
      PetscCallA(TSGetAdapt(temperatureTS, temperatureTSAdapt, ierr))
   case (MEF90HeatXfer_timeSteppingTypeNULL)
      continue
   end select

   select case (MEF90DefMechGlobalOptions%timeSteppingType)
   case (MEF90DefMech_TimeSteppingTypeQuasiStatic)
      PetscCallA(MEF90DefMechCreateSNESDisplacement(MEF90DefMechCtx, displacementSNES, displacementResidual, ierr))
      select case (MEF90DefMechGlobalOptions%damageSolverType)
      case (MEF90DefMech_DamageSolverTypeSNES)
         PetscCallA(MEF90DefMechCreateSNESDamage(MEF90DefMechCtx, damageSNES, damageResidual, ierr))
      case (MEF90DefMech_DamageSolverTypeTao)
         PetscCallA(MEF90DefMechCreateTAODamage(MEF90DefMechCtx, damageTAO, damageResidual, ierr))
         PetscCallA(TAOSetSolution(damageTAO, damage, ierr))
      end select ! MEF90DefMechGlobalOptions%damageSolverType
   case (MEF90DefMech_TimeSteppingTypeNULL)
      continue
   end select

   !!!
   !!! Allocate array of works and energies
   !!!
   allocate (elasticEnergy(size(MEF90HeatXferCtx%CellSetOptionsBag)))
   allocate (bodyForceWork(size(MEF90HeatXferCtx%CellSetOptionsBag)))
   allocate (cohesiveEnergy(size(MEF90HeatXferCtx%CellSetOptionsBag)))
   allocate (surfaceEnergy(size(MEF90HeatXferCtx%CellSetOptionsBag)))
   allocate (boundaryForceWork(size(MEF90HeatXferCtx%FaceSetOptionsBag)))

   !!!
   !!! Format Exodus file if needed
   !!!
   PetscCallA(DMGetDimension(MEF90DefMechCtx%megaDM, MEF90DefMechCtx%dim, ierr))
   PetscCallA(MEF90CtxGetTime(MEF90Ctx, time, ierr))
   if (EXONeedsFormatting) then
      if (MEF90GlobalOptions%verbose > 0) then
         PetscCallA(PetscPrintf(MEF90Ctx%comm, "Formatting result file\n", ierr))
      end if
      PetscCallA(MEF90DefMechFormatEXO(MEF90DefMechCtx, time, ierr))
      if (MEF90GlobalOptions%verbose > 0) then
         PetscCallA(PetscPrintf(MEF90Ctx%comm, "Done Formatting result file\n", ierr))
      end if
   end if
   if (MEF90GlobalOptions%verbose > 0) then
      PetscCallA(PetscViewerView(MEF90Ctx%resultViewer, PETSC_VIEWER_STDOUT_WORLD, ierr))
   end if

   !!!
   !!! Actual computations / time stepping
   !!!
   if (((MEF90DefMechGlobalOptions%timeSteppingType /= MEF90DefMech_TimeSteppingTypeNULL) .or. (MEF90HeatXferGlobalOptions%timeSteppingType /= MEF90DefMech_TimeSteppingTypeNULL)) .and. &
       (.not. MEF90GlobalOptions%dryrun)) then

         !!! Reload current state if necessary
      if (MEF90GlobalOptions%timeSkip > 0) then
         EXOstep = MEF90GlobalOptions%timeSkip
         select case (MEF90HeatXferGlobalOptions%timeSteppingType)
         case (MEF90HeatXfer_timeSteppingTypeSteadyState)
            PetscCallA(MEF90EXOVecLoad(MEF90HeatXferCtx%temperatureLocal, MEF90HeatXferCtx%temperatureToIOSF, MEF90HeatXferCtx%IOToTemperatureSF, MEF90Ctx%resultViewer, EXOstep, 1_ki, ierr))
         case (MEF90HeatXfer_timeSteppingTypeTransient)
            PetscCallA(TSSetTime(temperatureTS, time(EXOstep), ierr))
         end select

         select case (MEF90DefMechGlobalOptions%timeSteppingType)
         case (MEF90DefMech_timeSteppingTypeQuasiStatic)
            PetscCallA(MEF90EXOVecLoad(MEF90DefMechCtx%displacementLocal, MEF90DefMechCtx%displacementToIOSF, MEF90DefMechCtx%IOToDisplacementSF, MEF90Ctx%resultViewer, EXOstep, MEF90DefMechCtx%dim, ierr))
            PetscCallA(MEF90EXOVecLoad(MEF90DefMechCtx%damageLocal, MEF90DefMechCtx%damageToIOSF, MEF90DefMechCtx%IOToDamageSF, MEF90Ctx%resultViewer, EXOstep, 1_ki, ierr))
         end select
      end if

      step = MEF90GlobalOptions%timeSkip + 1
      mainloopQS: do
         MEF90DefMechCtx%analysisTime = time(step)
         if (step > 1) then
            MEF90DefMechCtx%timeStep = time(step) - time(step - 1)
         else
            MEF90DefMechCtx%timeStep = 0.0_kr
         end if

         write (IOBuffer, 100) step, time(step)
         PetscCallA(PetscPrintf(MEF90Ctx%comm, IOBuffer, ierr))

         !!! Solve for temperature
         PetscCallA(PetscLogStagePush(logStageHeatXfer, ierr))
         select case (MEF90HeatXferGlobalOptions%timeSteppingType)
         case (MEF90HeatXfer_timeSteppingTypeSteadyState)
            PetscCallA(MEF90HeatXferSetTransients(MEF90HeatXferCtx, step, time(step), ierr))
            PetscCallA(DMLocalToGlobal(temperatureDM, MEF90HeatXferCtx%temperatureLocal, INSERT_VALUES, temperature, ierr))
            !!! Solve SNES
            PetscCallA(SNESSolve(temperatureSNES, PETSC_NULL_VEC, temperature, ierr))
            PetscCallA(SNESGetConvergedReason(temperatureSNES, temperatureSNESConvergedReason, ierr))
            if (temperatureSNESConvergedReason%v < 0) then
               write (IOBuffer, 400) "temperature", temperatureSNESConvergedReason
               PetscCallA(PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr))
            end if

            PetscCallA(DMGlobalToLocal(temperatureDM, temperature, INSERT_VALUES, MEF90HeatXferCtx%temperatureLocal, ierr))
            PetscCallA(VecCopy(MEF90HeatXferCtx%temperatureLocal, MEF90DefMechCtx%temperatureLocal, ierr))
         case (MEF90HeatXfer_timeSteppingTypeTransient)
            if (step > 1) then
               write (IOBuffer, 200) step, time(step)
               PetscCallA(PetscPrintf(MEF90Ctx%comm, IOBuffer, ierr))
               if (step > 1) then
                  !!! Update fields
                  PetscCallA(MEF90HeatXferSetTransients(MEF90HeatXferCtx, step, time(step), ierr))
                  PetscCallA(TSSetMaxTime(temperatureTS, time(step), ierr))
                  PetscCallA(DMLocalToGlobal(temperatureDM, MEF90HeatXferCtx%temperatureLocal, INSERT_VALUES, temperature, ierr))
                  PetscCallA(TSSolve(temperatureTS, temperature, ierr))
                  PetscCallA(DMGlobalToLocal(temperatureDM, temperature, INSERT_VALUES, MEF90HeatXferCtx%temperatureLocal, ierr))
               end if
            end if
         end select
         PetscCallA(PetscLogStagePop(ierr))

         if (MEF90HeatXferGlobalOptions%timeSteppingType /= MEF90DefMech_TimeSteppingTypeNULL) then
            !!! Compute energies
            PetscCallA(PetscLogStagePush(logStageEnergy, ierr))
            PetscCallA(MEF90HeatXFerEnergy(MEF90HeatXferCtx, elasticEnergy, bodyForceWork, boundaryForceWork, ierr))
            PetscCallA(DMGetLabelIdIS(temperatureDM, MEF90CellSetLabelName, setIS, ierr))
            PetscCallA(MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm, setIS, ierr))
            PetscCallA(ISGetIndices(setIS, setID, ierr))
            do set = 1, size(setID)
               write (IOBuffer, 101) setID(set), elasticEnergy(set), bodyForceWork(set), elasticEnergy(set) - bodyForceWork(set)
               PetscCallA(PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr))
            end do
            PetscCallA(ISRestoreIndices(setIS, setID, ierr))
            PetscCallA(ISDestroy(setIS, ierr))

            PetscCallA(DMGetLabelIdIS(temperatureDM, MEF90FaceSetLabelName, setIS, ierr))
            PetscCallA(MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm, setIS, ierr))
            PetscCallA(ISGetIndices(setIS, setID, ierr))
            do set = 1, size(setID)
               write (IOBuffer, 103) setID(set), boundaryForceWork(set)
               PetscCallA(PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr))
            end do
            PetscCallA(ISRestoreIndices(setIS, setID, ierr))
            PetscCallA(ISDestroy(setIS, ierr))

            write (IOBuffer, 102) sum(elasticEnergy), sum(bodyForceWork) + sum(boundaryForceWork), sum(elasticEnergy) - sum(bodyForceWork) - sum(boundaryForceWork)
            PetscCallA(PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr))
            PetscCallA(PetscLogStagePop(ierr))

            !!! Save results
            EXOstep = step
            PetscCallA(PetscLogStagePush(logStageIO, ierr))
            PetscCallA(MEF90HeatXferViewEXO(MEF90HeatXferCtx, EXOstep, ierr))
            PetscCallA(PetscLogStagePop(ierr))
         end if

         !!! Solve for displacement and damage
         PetscCallA(MEF90DefMechSetTransients(MEF90DefMechCtx, step, time(step), ierr))
         select case (MEF90DefMechGlobalOptions%damageSolverType)
         case (MEF90DefMech_DamageSolverTypeSNES)
            PetscCallA(MEF90DefMechUpdateDamageBounds(MEF90DefMechCtx, damageSNES, damage, ierr))
         case (MEF90DefMech_DamageSolverTypeTao)
            PetscCallA(MEF90DefMechTAOUpdateDamageBounds(MEF90DefMechCtx, damageTAO, damage, ierr))
         end select ! MEF90DefMechGlobalOptions%damageSolverType
         PetscCallA(DMLocalToGlobal(displacementDM, MEF90DefMechCtx%displacementLocal, INSERT_VALUES, displacement, ierr))

         select case (MEF90DefMechGlobalOptions%timeSteppingType)
         case (MEF90DefMech_timeSteppingTypeQuasiStatic)
            select case (MEF90DefMechGlobalOptions%SolverType)
            case (MEF90DefMech_SolverTypeAltMin)
               PetscCallA(SNESSetLagPreconditioner(displacementSNES, 1_ki, ierr))
               if (MEF90DefMechGlobalOptions%damageSolverType == MEF90DefMech_DamageSolverTypeSNES) then
                  PetscCallA(SNESSetLagPreconditioner(damageSNES, 1_ki, ierr))
               end if

               AltMin: do AltMinIter = 1, MEF90DefMechGlobalOptions%damageMaxIt
                  AltMinStep = AltMinStep + 1
                  write (IObuffer, 208) AltMinIter
                  PetscCallA(PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr))

                  if (mod(AltMinIter - 1, MEF90DefMechGlobalOptions%PCLag) == 0) then
                     PetscCallA(SNESSetLagPreconditioner(displacementSNES, -2_ki, ierr))
                     if (MEF90DefMechGlobalOptions%damageSolverType == MEF90DefMech_DamageSolverTypeSNES) then
                        PetscCallA(SNESSetLagPreconditioner(damageSNES, 2_ki, ierr))
                     end if
                  end if

                  !!! Solve SNES displacement
                  PetscCallA(PetscLogStagePush(logStageDisplacement, ierr))
                  PetscCallA(SNESSolve(displacementSNES, PETSC_NULL_VEC, displacement, ierr))
                  PetscCallA(SNESGetConvergedReason(displacementSNES, displacementSNESConvergedReason, ierr))
                  if (displacementSNESConvergedReason%v < 0) then
                     write (IOBuffer, 400) "displacement", displacementSNESConvergedReason
                     PetscCallA(PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr))
                  end if

                  PetscCallA(DMGlobalToLocal(displacementDM, displacement, INSERT_VALUES, MEF90DefMechCtx%displacementLocal, ierr))
                  PetscCallA(DMLocalToGlobal(damageDM, MEF90DefMechCtx%damageLocal, INSERT_VALUES, damage, ierr))
                  PetscCallA(PetscLogStagePop(ierr))

                  PetscCallA(PetscLogStagePush(logStageDamage, ierr))
                  PetscCallA(VecCopy(damage, damageAltMinOld, ierr))

                  !!! Solve for damage field
                  PetscCallA(DMLocalToGlobal(damageDM, MEF90DefMechCtx%damageLocal, INSERT_VALUES, damage, ierr))
                  PetscCallA(VecCopy(damage, damageAltMinOld, ierr))
                  select case (MEF90DefMechGlobalOptions%damageSolverType)
                  case (MEF90DefMech_DamageSolverTypeSNES)
                     PetscCallA(SNESSolve(damageSNES, PETSC_NULL_VEC, damage, ierr))
                     PetscCallA(SNESGetConvergedReason(damageSNES, damageSNESConvergedReason, ierr))
                     if (damageSNESConvergedReason%v < 0) then
                        write (IOBuffer, 400) "damage", damageSNESConvergedReason
                        PetscCallA(PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr))
                     end if
                  case (MEF90DefMech_DamageSolverTypeTao)
                     PetscCallA(TAOSolve(damageTAO, ierr))
                     PetscCallA(TAOGetConvergedReason(damageTAO, damageTAOConvergedReason, ierr))
                     if (damageTAOConvergedReason%v < 0) then
                        write (IOBuffer, 401) "damage", damageTAOConvergedReason
                        PetscCallA(PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr))
                     end if
                     PetscCallA(TAOGetSolution(damageTAO, damage, ierr))
                  end select ! MEF90DefMechGlobalOptions%damageSolverType
                  PetscCallA(DMGlobalToLocal(damageDM, damage, INSERT_VALUES, MEF90DefMechCtx%damageLocal, ierr))

                  !!! Over relaxation of the damage variable
                  if (AltMinIter > 1) then
                     if ((MEF90DefMechGlobalOptions%SOROmega > 0.0_kr) .and. (MEF90DefMechGlobalOptions%SOROmega /= 1.0)) then
                        mySOROmega = MEF90DefMechGlobalOptions%SOROmega
                        !!! LIMITED SOR
                        PetscCallA(SNESVIGetVariableBounds(damageSNES, damageLB, damageUB, ierr))
                        PetscCallA(VecGetArrayRead(damageLB, damageLBArray, ierr))
                        PetscCallA(VecGetArrayRead(damageUB, damageUBArray, ierr))
                        PetscCallA(VecGetArrayRead(damageAltMinOld, damageAltMinOldArray, ierr))
                        PetscCallA(VecGetArrayRead(damage, damageArray, ierr))
                        do iDof = 1, size(damageArray)
                           if (damageArray(iDof) > damageAltMinOldArray(iDof)) then
                              mySOROmega = min(mySOROmega, (damageUBArray(iDof) - damageAltMinOldArray(iDof)) / (damageArray(iDof) - damageAltMinOldArray(iDof)))
                           else if (damageArray(iDof) < damageAltMinOldArray(iDof)) then
                              mySOROmega = min(mySOROmega, (damageLBArray(iDof) - damageAltMinOldArray(iDof)) / (damageArray(iDof) - damageAltMinOldArray(iDof)))
                           end if
                        end do
                        PetscCallA(VecRestoreArrayRead(damage, damageArray, ierr))
                        PetscCallA(VecRestoreArrayRead(damageAltMinOld, damageAltMinOldArray, ierr))
                        PetscCallA(VecRestoreArrayRead(damageUB, damageUBArray, ierr))
                        PetscCallA(VecRestoreArrayRead(damageLB, damageLBArray, ierr))
                        PetscCallA(MPI_AllReduce(mySOROmega, SOROmega, 1, MPIU_SCALAR, MPI_MIN, MEF90Ctx%comm, ierr))
                        PetscCallA(VecAXPBY(damage, 1.0_kr - SOROmega, SOROmega, damageAltMinOld, ierr))
                     else if (MEF90DefMechGlobalOptions%SOROmega < 0.0_kr) then
                        !!! PROJECTED SOR
                        SOROmega = -MEF90DefMechGlobalOptions%SOROmega
                        PetscCallA(VecAXPBY(damage, 1.0_kr - SOROmega, SOROmega, damageAltMinOld, ierr))
                        PetscCallA(SNESVIGetVariableBounds(damageSNES, damageLB, damageUB, ierr))
                        PetscCallA(VecPointwiseMax(damage, damage, damageLB, ierr))
                        PetscCallA(VecPointwiseMin(damage, damage, damageUB, ierr))
                     end if
                  end if

                  !!! Monitor the progress of the Alt Min algorithm
                  PetscCallA(VecMin(damage, PETSC_NULL_INTEGER, damageMin, ierr))
                  PetscCallA(VecMax(damage, PETSC_NULL_INTEGER, damageMax, ierr))
                  PetscCallA(VecAxPy(damageAltMinOld, -1.0_kr, damage, ierr))
                  PetscCallA(VecNorm(damageAltMinOld, NORM_INFINITY, damageMaxChange, ierr))
                  write (IOBuffer, 209) damageMin, damageMax, damageMaxChange
                  PetscCallA(PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr))
                  PetscCallA(PetscLogStagePop(ierr))

                  !!! Test for convergence based on the L^\infty norm of the increment
                  if (damageMaxChange <= MEF90DefMechGlobalOptions%damageATol) then
                     exit altMin
                  end if

                  if (mod(AltMinIter, 25_Ki) == 0) then
                     EXOstep = step
                     !!! Save results and boundary Values
                     PetscCallA(PetscLogStagePush(logStageIO, ierr))
                     PetscCallA(MEF90DefMechViewEXO(MEF90DefMechCtx, EXOstep, ierr))
                     PetscCallA(PetscLogStagePop(ierr))
                  end if
               end do AltMin
            case default
               write (IOBuffer, *) "Unimplemented DefMech solver type: ", MEF90DefMechGlobalOptions%SolverType, "\n"
               PetscCallA(PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr))
               stop
            end select ! solverType

            !!! Compute energies
            PetscCallA(PetscLogStagePush(logStageEnergy, ierr))
            elasticEnergy = 0.0_kr
            bodyForceWork = 0.0_kr
            boundaryForceWork = 0.0_kr
            surfaceEnergy = 0.0_kr
            cohesiveEnergy = 0.0_kr
            PetscCallA(MEF90DefMechWork(MEF90DefMechCtx, bodyForceWork, boundaryForceWork, ierr))
            PetscCallA(MEF90DefMechElasticEnergy(MEF90DefMechCtx, elasticEnergy, ierr))
            PetscCallA(MEF90DefMechSurfaceEnergy(MEF90DefMechCtx, surfaceEnergy, ierr))

            PetscCallA(DMGetLabelIdIS(displacementDM, MEF90CellSetLabelName, setIS, ierr))
            PetscCallA(MEF90ISAllGatherMerge(MEF90Ctx%comm, setIS, ierr))
            PetscCallA(ISGetIndices(setIS, setID, ierr))
            do set = 1, size(setID)
               write (IOBuffer, 201) setID(set), elasticEnergy(set), bodyForceWork(set), cohesiveEnergy(set), surfaceEnergy(set), elasticEnergy(set) - bodyForceWork(set) + cohesiveEnergy(set) + surfaceEnergy(set)
               PetscCallA(PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr))
               write (IOBuffer, 500) step, time(step), elasticEnergy(set), bodyForceWork(set), cohesiveEnergy(set), surfaceEnergy(set), elasticEnergy(set) - bodyForceWork(set) + cohesiveEnergy(set) + surfaceEnergy(set)
               PetscCallA(PetscViewerASCIIPrintf(MEF90DefMechCtx%setEnergyViewer(set), IOBuffer, ierr))
               PetscCallA(PetscViewerFlush(MEF90DefMechCtx%setEnergyViewer(set), ierr))
            end do
            PetscCallA(ISRestoreIndices(setIS, setID, ierr))
            PetscCallA(ISDestroy(setIS, ierr))

            PetscCallA(DMGetLabelIdIS(displacementDM, MEF90FaceSetLabelName, setIS, ierr))
            PetscCallA(MEF90ISAllGatherMerge(MEF90Ctx%comm, setIS, ierr))
            PetscCallA(ISGetIndices(setIS, setID, ierr))
            do set = 1, size(setID)
               write (IOBuffer, 203) setID(set), boundaryForceWork(set)
               PetscCallA(PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr))
            end do
            PetscCallA(ISRestoreIndices(setIS, setID, ierr))
            PetscCallA(ISDestroy(setIS, ierr))

            write (IOBuffer, 202) sum(elasticEnergy), sum(bodyForceWork) + sum(boundaryForceWork), sum(cohesiveEnergy), sum(surfaceEnergy), sum(elasticEnergy) - sum(bodyForceWork) - sum(boundaryForceWork) + sum(cohesiveEnergy) + sum(surfaceEnergy)
            PetscCallA(PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr))
            write (IOBuffer, 500) step, time(step), sum(elasticEnergy), sum(bodyForceWork) + sum(boundaryForceWork), sum(cohesiveEnergy), sum(surfaceEnergy), sum(elasticEnergy) - sum(bodyForceWork) - sum(boundaryForceWork) + sum(cohesiveEnergy) + sum(surfaceEnergy)
            PetscCallA(PetscViewerASCIIPrintf(MEF90DefMechCtx%globalEnergyViewer, IOBuffer, ierr))
            PetscCallA(PetscViewerFlush(MEF90DefMechCtx%globalEnergyViewer, ierr))
            PetscCallA(PetscLogStagePop(ierr))

            !!! Save results and boundary Values
            if (MEF90DefMechGlobalOptions%stressExport) then
               PetscCallA(MEF90DefMechStress(MEF90DefMechCtx, MEF90DefMechCtx%stress, ierr))
            end if
            PetscCallA(PetscLogStagePush(logStageIO, ierr))
            EXOstep = step
            PetscCallA(MEF90DefMechViewEXO(MEF90DefMechCtx, EXOstep, ierr))
            PetscCallA(PetscLogStagePop(ierr))

            PetscCallA(PetscViewerASCIIOpen(MEF90Ctx%comm,trim(MEF90FilePrefix(MEF90Ctx%resultFile))//'.log',logViewer, ierr))
            PetscCallA(PetscLogView(logViewer,ierr))
            PetscCallA(PetscViewerDestroy(logViewer,ierr))
         end select ! timeStepingType

         if (step == MEF90GlobalOptions%timeNumStep) then
            exit mainloopQS
         else
            step = step + 1
         end if
      end do MainloopQS
   end if ! timeSteppingType
   write (IOBuffer, *) 'Total number of alternate minimizations:', AltMinStep, '\n'
   PetscCallA(PetscPrintf(MEF90Ctx%comm, IOBuffer, ierr))

100 format("\nSolving steady state step ", I4, ", t=", ES12.5, "\n")
200 format("\nSolving transient step ", I4, ", t=", ES12.5, "\n")
101 format("cell set ", I4, " thermal energy: ", ES12.5, " flux: ", ES12.5, " total: ", ES12.5, "\n")
102 format("======= Total thermal energy: ", ES12.5, " flux: ", ES12.5, " total: ", ES12.5, "\n")
103 format("face set ", I4, "                              flux: ", ES12.5, "\n")
201 format("cell set ", I4, "  elastic energy: ", ES12.5, " work: ", ES12.5, " cohesive: ", ES12.5, " surface: ", ES12.5, " total: ", ES12.5, "\n")
203 format("face set ", I4, "                      boundary work: ", ES12.5, "\n")
202 format("======= Total: elastic energy: ", ES12.5, " work: ", ES12.5, " cohesive: ", ES12.5, " surface: ", ES12.5, " total: ", ES12.5, "\n")
500 format(I6, 6(ES16.5), "\n")

208 format("   Alt. Min. step ", I5, " ")
209 format(" alpha min / max", ES12.5, " / ", ES12.5, ", max change ", ES12.5, "\n")

400 format(" [ERROR]: ", A, " SNESSolve failed with SNESConvergedReason ", I2, ". \n Check https://petsc.org/release/docs/manualpages/SNES/SNESConvergedReason/ for error code meaning.\n")
401 format(" [ERROR]: ", A, " TAOSolve failed with TAOConvergedReason ", I2, ". \n Check https://petsc.org/release/docs/manualpages/TAO/TAOConvergedReason/ for error code meaning.\n")

   !!! Clean up and exit nicely
   select case (MEF90DefMechGlobalOptions%timeSteppingType)
   case (MEF90DefMech_timeSTeppingTypeQuasiStatic)
      PetscCallA(SNESDestroy(displacementSNES, ierr))
      PetscCallA(VecDestroy(displacementResidual, ierr))
      PetscCallA(VecDestroy(displacement, ierr))
   end select

   !!! Clean up and exit nicely
   select case (MEF90HeatXferGlobalOptions%timeSteppingType)
   case (MEF90HeatXFer_timeSteppingTypeSteadyState)
      PetscCallA(SNESDestroy(temperatureSNES, ierr))
   case (MEF90HeatXFer_timeSteppingTypeTransient)
      PetscCallA(TSDestroy(temperatureTS, ierr))
   end select
   PetscCallA(VecDestroy(temperatureResidual, ierr))
   PetscCallA(VecDestroy(temperature, ierr))

   select case (MEF90DefMechGlobalOptions%timeSteppingType)
   case (MEF90DefMech_TimeSteppingTypeQuasiStatic)
      PetscCallA(SNESDestroy(displacementSNES, ierr))
      select case (MEF90DefMechGlobalOptions%damageSolverType)
      case (MEF90DefMech_DamageSolverTypeSNES)
         PetscCallA(SNESDestroy(damageSNES, ierr))
      case (MEF90DefMech_DamageSolverTypeTao)
         PetscCallA(TAODestroy(damageTAO, ierr))
      end select ! MEF90DefMechGlobalOptions%damageSolverType
   end select
   PetscCallA(VecDestroy(displacementResidual, ierr))
   PetscCallA(VecDestroy(displacement, ierr))
   PetscCallA(VecDestroy(damageResidual, ierr))
   PetscCallA(VecDestroy(damageAltMinOld, ierr))
   PetscCallA(VecDestroy(damage, ierr))

   deallocate (time)
   deallocate (elasticEnergy)
   deallocate (bodyForceWork)
   deallocate (boundaryForceWork)
   deallocate (cohesiveEnergy)
   deallocate (surfaceEnergy)
   PetscCallA(MEF90DefMechCtxDestroy(MEF90DefMechCtx, ierr))
   nullify (MEF90HeatXferCtx%temperatureLocal)
   PetscCallA(MEF90HeatXferCtxDestroy(MEF90HeatXferCtx, ierr))

   PetscCallA(PetscViewerDestroy(MEF90Ctx%resultViewer, ierr))
   If (.NOT. MEF90GlobalOptions%dryrun) Then
      PetscCallA(PetscViewerASCIIOpen(MEF90Ctx%comm,trim(MEF90FilePrefix(MEF90Ctx%resultFile))//'.log',logViewer, ierr))
      PetscCallA(PetscLogView(logViewer,ierr))
      PetscCallA(PetscViewerDestroy(logViewer,ierr))
   End If
   PetscCallA(MEF90CtxDestroy(MEF90Ctx, ierr))
   PetscCallA(MEF90Finalize(ierr))
   PetscCallA(PetscFinalize(ierr))
end program vDef
