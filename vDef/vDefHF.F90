#include "../MEF90/mef90.inc"
program vDefHF
#include <finclude/petscdef.h>

   use petsc
   use m_MEF90
   use m_MEF90_DefMechCtx
   use m_MEF90_DefMech
   use m_MEF90_HeatXferCtx
   use m_MEF90_HeatXfer
   use m_vDefDefault
   implicit none(type, external)

   PetscErrorCode                                     :: ierr
   type(MEF90Ctx_Type), target                         :: MEF90Ctx
   type(MEF90CtxGlobalOptions_Type), pointer           :: MEF90GlobalOptions

   !!! Defect mechanics contexts
   type(MEF90DefMechCtx_Type)                         :: MEF90DefMechCtx
   type(MEF90DefMechGlobalOptions_Type), pointer       :: MEF90DefMechGlobalOptions

   !!! HeatXfer contexts
   type(MEF90HeatXferCtx_Type)                        :: MEF90HeatXferCtx
   type(MEF90HeatXferGlobalOptions_Type), pointer      :: MEF90HeatXferGlobalOptions

   type(DM), target                                    :: Mesh
   type(IS)                                           :: setIS, CellSetGlobalIS
   PetscInt, dimension(:), pointer                      :: setID
   PetscInt                                           :: numset, set
   PetscReal, dimension(:), pointer                     :: time
   PetscReal, dimension(:), pointer                     :: thermalEnergySet, heatFluxWorkSet
   PetscReal, dimension(:), pointer                     :: elasticEnergySet, surfaceEnergySet, forceWorkSet, cohesiveEnergySet
   PetscReal, dimension(:), pointer                     :: elasticEnergy, surfaceEnergy, forceWork, cohesiveEnergy, totalMechanicalEnergy

   SNESConvergedReason                                :: snesDispConvergedReason
   type(Vec)                                          :: residualDisp
   SNESConvergedReason                                :: snesDamageConvergedReason
   type(Vec)                                          :: residualDamage, damageAltMinOld, localVec
   PetscInt                                           :: AltMinIter
   PetscReal                                          :: damageMaxChange, CrackPressureMaxChange

   type(SNES)                                         :: snesTemp
   SNESConvergedReason                                :: snesTempConvergedReason
   type(TS)                                           :: tsTemp
   TSConvergedReason                                  :: tsTempConvergedReason
   type(TSAdapt)                                      :: tsAdaptTemp
   type(Vec)                                          :: residualTemp

   PetscReal                                          :: tsTempInitialStep, tsTempInitialTime
   PetscInt                                           :: tsTempmaxIter
   PetscReal                                          :: t

   PetscBool                                          :: flg
   character(len=MEF90MXSTRLEN)                      :: IOBuffer
   type(PetscViewer)                                  :: logViewer, pressureViewer
   integer                                            :: numfield

   integer                                            :: step
   PetscInt                                           :: dim

   PetscReal                                          :: alphaMaxChange, alphaMin, alphaMax
   PetscLogStage                                      :: logStageHeat, logStageDamage, logStageDisplacement, logStageEnergy, logStageIO

   !!! WorkControlled and CrackPressure variables
   PetscReal                                          :: WorkControlledRescaling, ErrorEstimationWorkControlled
   PetscBool, dimension(:), pointer                     :: ActivatedWorkControlledBlocksList, ActivatedCrackPressureBlocksList
   PetscReal, dimension(:), pointer                     :: CrackVolumeSet, WorkControlled

   type(MEF90DefMechCellSetOptions_Type), pointer      :: cellSetOptions

   !!! Secant method for sneddon
   PetscReal, dimension(3)                              :: CrackPressureSave, CrackVolumeSave
   PetscInt                                            :: CrackVolumeIter, I1, I2, I3
   PetscReal                                           :: InjectedVolumeRelativeError = 1
   type(Vec)                                           :: CrackPressureMask

   !!! Initialize MEF90
   call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
   call MEF90Initialize(PETSC_COMM_WORLD, ierr)
   call PetscPrintf(PETSC_COMM_WORLD, " # vDefHF: numerical implementation of variational models of fracture for Hydraulic Fracturing\n", ierr); CHKERRQ(ierr)
   call PetscLogStageRegister('HeatXfer    ', logStageHeat, ierr); CHKERRQ(ierr)
   call PetscLogStageRegister('Damage      ', logStageDamage, ierr); CHKERRQ(ierr)
   call PetscLogStageRegister('Displacement', logStageDisplacement, ierr); CHKERRQ(ierr)
   call PetscLogStageRegister('Energy      ', logStageEnergy, ierr); CHKERRQ(ierr)
   call PetscLogStageRegister('IO          ', logStageIO, ierr); CHKERRQ(ierr)

   !!! Get all MEF90-wide options
   call MEF90CtxCreate(PETSC_COMM_WORLD, MEF90Ctx, MEF90CtxDefaultGlobalOptions, ierr); CHKERRQ(ierr)
   call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx % GlobalOptionsBag, MEF90GlobalOptions, ierr); CHKERRQ(ierr)

   !!! Get DM from mesh
   call MEF90CtxGetDMMeshEXO(MEF90Ctx, Mesh, ierr); CHKERRQ(ierr)
   call DMMeshGetDimension(Mesh, dim, ierr); CHKERRQ(ierr)
   call DMMeshSetMaxDof(Mesh, dim, ierr); CHKERRQ(ierr)
   call DMSetBlockSize(Mesh, dim, ierr); CHKERRQ(ierr)

   !!! Open output file
   call MEF90CtxOpenEXO(MEF90Ctx, Mesh, ierr)

   !!! Create DefMech context, get all DefMech options
   call MEF90DefMechCtxCreate(MEF90DefMechCtx, Mesh, MEF90Ctx, ierr); CHKERRQ(ierr)
   if (dim == 2) then
      call MEF90DefMechCtxSetFromOptions(MEF90DefMechCtx, PETSC_NULL_CHARACTER, vDefDefMechDefaultGlobalOptions2D, &
                                         DefMechDefaultCellSetOptions, DefMechDefaultVertexSetOptions, ierr)
   else
      call MEF90DefMechCtxSetFromOptions(MEF90DefMechCtx, PETSC_NULL_CHARACTER, vDefDefMechDefaultGlobalOptions3D, &
                                         DefMechDefaultCellSetOptions, DefMechDefaultVertexSetOptions, ierr)
   end if
   call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx % GlobalOptionsBag, MEF90DefMechGlobalOptions, ierr); CHKERRQ(ierr)

   !!! Create HeatXfer context, get all HeatXfer options
   call MEF90HeatXferCtxCreate(MEF90HeatXferCtx, Mesh, MEF90Ctx, ierr); CHKERRQ(ierr)
   call MEF90HeatXferCtxSetFromOptions(MEF90HeatXferCtx, PETSC_NULL_CHARACTER, HeatXferDefaultGlobalOptions, &
                                       HeatXferDefaultCellSetOptions, HeatXferDefaultVertexSetOptions, ierr)
   call PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx % GlobalOptionsBag, MEF90HeatXferGlobalOptions, ierr); CHKERRQ(ierr)

   !!! Get material properties bags
   if (dim == 2) then
      call MEF90MatPropBagSetFromOptions(MEF90DefMechCtx % MaterialPropertiesBag, MEF90DefMechCtx % DMVect, MEF90Mathium2D, MEF90Ctx, ierr)
   else
      call MEF90MatPropBagSetFromOptions(MEF90DefMechCtx % MaterialPropertiesBag, MEF90DefMechCtx % DMVect, MEF90Mathium3D, MEF90Ctx, ierr)
   end if
   MEF90HeatXferCtx % MaterialPropertiesBag => MEF90DefMechCtx % MaterialPropertiesBag

   !!! Create time array from global options
   call MEF90CtxGetTime(MEF90Ctx, time, ierr)

   !!! Create sections, vectors, and solvers for DefMech Context
   call MEF90DefMechCtxSetSections(MEF90DefMechCtx, ierr)
   call MEF90DefMechCtxCreateVectors(MEF90DefMechCtx, ierr)
   call VecDuplicate(MEF90DefMechCtx % damage, damageAltMinOld, ierr); CHKERRQ(ierr)
   call VecDuplicate(MEF90DefMechCtx % displacement, residualDisp, ierr); CHKERRQ(ierr)
   call PetscObjectSetName(residualDisp, "residualDisp", ierr); CHKERRQ(ierr)
   call MEF90DefMechCreateSNESDisplacement(MEF90DefMechCtx, MEF90DefMechCtx % snesDisp, residualDisp, ierr)
   call VecDuplicate(MEF90DefMechCtx % damage, residualDamage, ierr); CHKERRQ(ierr)
   call PetscObjectSetName(residualDamage, "residualDamage", ierr); CHKERRQ(ierr)
   call MEF90DefMechCreateSNESDamage(MEF90DefMechCtx, MEF90DefMechCtx % snesDamage, residualDamage, ierr)
   deallocate (MEF90DefMechCtx % temperature)

   !!! Create sections, vectors, and solvers for HeatXfer Context
   if (MEF90HeatXferGlobalOptions % timeSteppingType /= MEF90HeatXfer_timeSteppingTypeNULL) then
      call MEF90HeatXferCtxSetSections(MEF90HeatXferCtx, ierr)
      call MEF90HeatXferCtxCreateVectors(MEF90HeatXferCtx, ierr)
      call VecDuplicate(MEF90HeatXferCtx % temperature, residualTemp, ierr); CHKERRQ(ierr)
      call PetscObjectSetName(residualTemp, "residualTemp", ierr); CHKERRQ(ierr)
      select case (MEF90HeatXferGlobalOptions % timeSteppingType)
      case (MEF90HeatXfer_timeSteppingTypeSteadyState)
         call MEF90HeatXferCreateSNES(MEF90HeatXferCtx, snesTemp, residualTemp, ierr)
      case (MEF90HeatXfer_timeSteppingTypeTransient)
         tsTempInitialStep = (time(size(time)) - time(1)) / (size(time) - 1.0_kr) / 10.0_kr
         tsTempInitialTime = time(1)
         call MEF90HeatXferCreateTS(MEF90HeatXferCtx, tsTemp, residualTemp, tsTempInitialTime, tsTempInitialStep, ierr)
         call TSGetAdapt(tsTemp, tsAdaptTemp, ierr); CHKERRQ(ierr)
         call TSAdaptSetFromOptions(tsAdaptTemp, ierr); CHKERRQ(ierr)
      end select

      !!! Link the temperature field in the DefMechContext with that of the HeatXfer
      MEF90DefMechCtx % temperature => MEF90HeatXferCtx % temperature
   end if
   !!!
   !!! Allocate array of works and energies
   !!!
   allocate (elasticEnergySet(size(MEF90DefMechCtx % CellSetOptionsBag)))
   elasticEnergySet = 0.0_kr
   allocate (surfaceEnergySet(size(MEF90DefMechCtx % CellSetOptionsBag)))
   surfaceEnergySet = 0.0_kr
   allocate (forceWorkSet(size(MEF90DefMechCtx % CellSetOptionsBag)))
   forceWorkSet = 0.0_kr
   allocate (cohesiveEnergySet(size(MEF90DefMechCtx % CellSetOptionsBag)))
   cohesiveEnergySet = 0.0_kr
   allocate (thermalEnergySet(size(MEF90DefMechCtx % CellSetOptionsBag)))
   thermalEnergySet = 0.0_kr
   allocate (heatFluxWorkSet(size(MEF90DefMechCtx % CellSetOptionsBag)))
   heatFluxWorkSet = 0.0_kr
   allocate (CrackVolumeSet(size(MEF90DefMechCtx % CellSetOptionsBag)))
   CrackVolumeSet = 0.0_kr

   allocate (elasticEnergy(MEF90GlobalOptions % timeNumStep))
   elasticEnergy = 0.0_kr
   allocate (surfaceEnergy(MEF90GlobalOptions % timeNumStep))
   surfaceEnergy = 0.0_kr
   allocate (forceWork(MEF90GlobalOptions % timeNumStep))
   forceWork = 0.0_kr
   allocate (cohesiveEnergy(MEF90GlobalOptions % timeNumStep))
   cohesiveEnergy = 0.0_kr
   allocate (totalMechanicalEnergy(MEF90GlobalOptions % timeNumStep))
   totalMechanicalEnergy = 0.0_kr
   allocate (WorkControlled(MEF90GlobalOptions % timeNumStep))
   WorkControlled = 0.0_kr

   !!!
   !!! Try to figure out if the file was formatted
   !!!
   if (MEF90Ctx % rank == 0) then
      call EXGVP(MEF90Ctx % fileExoUnit, "N", numfield, ierr)
   end if
   call MPI_Bcast(numfield, 1, MPIU_INTEGER, 0, MEF90Ctx % comm, ierr)
   if (numfield == 0) then
      call MEF90DefMechFormatEXO(MEF90DefMechCtx, time, ierr)
   end if

   !!! Create logical list of blocks where crack pressure or work control is activated
   allocate (ActivatedCrackPressureBlocksList(size(MEF90DefMechCtx % CellSetOptionsBag)))
   allocate (ActivatedWorkControlledBlocksList(size(MEF90DefMechCtx % CellSetOptionsBag)))
   call DMmeshGetLabelIdIS(MEF90DefMechCtx % CellDMVect, 'Cell Sets', CellSetGlobalIS, ierr); CHKERRQ(ierr)
   call MEF90ISAllGatherMerge(PETSC_COMM_WORLD, CellSetGlobalIS, ierr); CHKERRQ(ierr)
   call ISGetIndices(CellSetGlobalIS, setID, ierr); CHKERRQ(ierr)
   do set = 1, size(setID)
      call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx % CellSetOptionsBag(set), cellSetOptions, ierr); CHKERRQ(ierr)
      ActivatedCrackPressureBlocksList(set) = cellSetOptions % CrackVolumeControlled
      ActivatedWorkControlledBlocksList(set) = cellSetOptions % WorkControlled
   end do
   call ISRestoreIndices(CellSetGlobalIS, setID, ierr); CHKERRQ(ierr)
   call ISDestroy(CellSetGlobalIS, ierr); CHKERRQ(ierr)

   call VecDuplicate(MEF90DefMechCtx % CrackPressure, CrackPressureMask, ierr); CHKERRQ(ierr)

   !!! Create a viewer for the volume-pressure curve
   call PetscViewerASCIIOpen(MEF90Ctx % comm, trim(MEF90FilePrefix(MEF90Ctx % resultFile))//'.pres', pressureViewer, ierr); CHKERRQ(ierr)
   call PetscViewerASCIIPrintf(pressureViewer, "# Iteration   volume   pressure\n", ierr); CHKERRQ(ierr)
   call PetscViewerFlush(pressureViewer, ierr); CHKERRQ(ierr)

   !!! Actual computations / time stepping
   !!!
   if (.not. MEF90GlobalOptions % dryrun) then
      if (MEF90GlobalOptions % timeSkip > 0) then
         !!! Restore state from file.
         call DMGetLocalVector(MEF90DefMechCtx % DMScal, localVec, ierr); CHKERRQ(ierr)
         call VecLoadExodusVertex(MEF90DefMechCtx % DMScal, localVec, MEF90DefMechCtx % MEF90Ctx % IOcomm, &
                                  MEF90DefMechCtx % MEF90Ctx % fileExoUnit, MEF90GlobalOptions % timeSkip, MEF90DefMechGlobalOptions % damageOffset, ierr); CHKERRQ(ierr)
         call DMLocalToGlobalBegin(MEF90DefMechCtx % DMScal, localVec, INSERT_VALUES, MEF90DefMechCtx % damage, ierr); CHKERRQ(ierr)
         call DMLocalToGlobalEnd(MEF90DefMechCtx % DMScal, localVec, INSERT_VALUES, MEF90DefMechCtx % damage, ierr); CHKERRQ(ierr)
         call VecCopy(MEF90DefMechCtx % damage, damageAltMinOld, ierr); CHKERRQ(ierr)
         call DMRestoreLocalVector(MEF90DefMechCtx % DMScal, localVec, ierr); CHKERRQ(ierr)

         if (MEF90DefMechGlobalOptions % crackPressureOffset > 0) then
            call DMGetLocalVector(MEF90DefMechCtx % cellDMScal, localVec, ierr); CHKERRQ(ierr)
            call VecLoadExodusCell(MEF90DefMechCtx % cellDMScal, localVec, MEF90DefMechCtx % MEF90Ctx % IOcomm, &
                                   MEF90DefMechCtx % MEF90Ctx % fileExoUnit, MEF90GlobalOptions % timeSkip, MEF90DefMechGlobalOptions % crackPressureOffset, ierr); CHKERRQ(ierr)
            call DMLocalToGlobalBegin(MEF90DefMechCtx % cellDMScal, localVec, INSERT_VALUES, MEF90DefMechCtx % crackPressure, ierr); CHKERRQ(ierr)
            call DMLocalToGlobalEnd(MEF90DefMechCtx % cellDMScal, localVec, INSERT_VALUES, MEF90DefMechCtx % crackPressure, ierr); CHKERRQ(ierr)
            call DMRestoreLocalVector(MEF90DefMechCtx % cellDMScal, localVec, ierr); CHKERRQ(ierr)
         end if

         call DMGetLocalVector(MEF90DefMechCtx % DMVect, localVec, ierr); CHKERRQ(ierr)
         call VecLoadExodusVertex(MEF90DefMechCtx % DMVect, localVec, MEF90DefMechCtx % MEF90Ctx % IOcomm, &
                                  MEF90DefMechCtx % MEF90Ctx % fileExoUnit, MEF90GlobalOptions % timeSkip, MEF90DefMechGlobalOptions % displacementOffset, ierr); CHKERRQ(ierr)
         call DMLocalToGlobalBegin(MEF90DefMechCtx % DMVect, localVec, INSERT_VALUES, MEF90DefMechCtx % Displacement, ierr); CHKERRQ(ierr)
         call DMLocalToGlobalEnd(MEF90DefMechCtx % DMVect, localVec, INSERT_VALUES, MEF90DefMechCtx % Displacement, ierr); CHKERRQ(ierr)
         call DMRestoreLocalVector(MEF90DefMechCtx % DMVect, localVec, ierr); CHKERRQ(ierr)
      end if
      step = MEF90GlobalOptions % timeSkip + 1
      mainloopQS: do
         if (step > 1) then
            MEF90DefMechCtx % timeStep = time(step) - time(step - 1)
         else
            MEF90DefMechCtx % timeStep = 0.0_kr
         end if

         call PetscLogStagePush(logStageHeat, ierr); CHKERRQ(ierr)
         !!! Solve for temperature
         select case (MEF90HeatXferGlobalOptions % timeSteppingType)
         case (MEF90HeatXfer_timeSteppingTypeSteadyState)
            write (IOBuffer, 100) step, time(step)
            call PetscPrintf(MEF90Ctx % comm, IOBuffer, ierr); CHKERRQ(ierr)

            !!! Update fields
            call MEF90HeatXferSetTransients(MEF90HeatXferCtx, step, time(step), ierr)
            !!! Solve SNES
            call MEF90HeatXferUpdateboundaryTemperature(MEF90HeatXferCtx % temperature, MEF90HeatXferCtx, ierr); 
            call SNESSolve(snesTemp, PETSC_NULL_OBJECT, MEF90HeatXferCtx % temperature, ierr); CHKERRQ(ierr)
            call SNESGetConvergedReason(snesTemp, snesTempConvergedReason, ierr); CHKERRQ(ierr)
            if (snesTempConvergedReason < 0) then
               write (IOBuffer, 400) "temperature", snesTempConvergedReason
               call PetscPrintf(MEF90Ctx % Comm, IOBuffer, ierr); CHKERRQ(ierr)
            end if

            !!! Compute thermal energy
            call MEF90HeatXFerEnergy(MEF90HeatXferCtx % temperature, time(step), MEF90HeatXferCtx, thermalEnergySet, heatFluxWorkSet, ierr); CHKERRQ(ierr)
            call DMmeshGetLabelIdIS(MEF90HeatXferCtx % DM, 'Cell Sets', CellSetGlobalIS, ierr); CHKERRQ(ierr)
            call MEF90ISAllGatherMerge(MEF90Ctx % Comm, CellSetGlobalIS, ierr); CHKERRQ(ierr)
            call ISGetIndices(CellSetGlobalIS, setID, ierr); CHKERRQ(ierr)
            call PetscPrintf(MEF90Ctx % Comm, "\nThermal energies: \n", ierr); CHKERRQ(ierr)
            do set = 1, size(setID)
               write (IOBuffer, 101) setID(set), thermalEnergySet(set), heatFluxWorkSet(set), thermalEnergySet(set) - heatFluxWorkSet(set)
               call PetscPrintf(MEF90Ctx % Comm, IOBuffer, ierr); CHKERRQ(ierr)
            end do
            call ISRestoreIndices(CellSetGlobalIS, setID, ierr); CHKERRQ(ierr)
            call ISDestroy(CellSetGlobalIS, ierr); CHKERRQ(ierr)
            write (IOBuffer, 102) sum(thermalEnergySet), sum(heatFluxWorkSet), sum(thermalEnergySet) - sum(heatFluxWorkSet)
            call PetscPrintf(MEF90Ctx % Comm, IOBuffer, ierr); CHKERRQ(ierr)

            !!! Save results
            call MEF90HeatXferViewEXO(MEF90HeatXferCtx, step, ierr)
         case (MEF90HeatXfer_timeSteppingTypeTransient)
            if (step > 1) then
               write (IOBuffer, 110) step, time(step)
               call PetscPrintf(MEF90Ctx % comm, IOBuffer, ierr); CHKERRQ(ierr)
               !!! Update fields
               call MEF90HeatXferSetTransients(MEF90HeatXferCtx, step, time(step), ierr)
               call MEF90HeatXferUpdateboundaryTemperature(MEF90HeatXferCtx % temperature, MEF90HeatXferCtx, ierr); 
               !!! Make sure TS does not overstep
               call TSGetTime(tsTemp, t, ierr); CHKERRQ(ierr)
               if (t < time(step)) then
                  call TSAdaptSetStepLimits(tsAdaptTemp, PETSC_DECIDE, (time(step) - time) / 2.0_kr, ierr); CHKERRQ(ierr)
                  !!! Something is up here.
                  !!! replacing the constant 10000 with a variable leads to divergence of TSAdapt
                  !!! when using gcc
                  call TSSetDuration(tsTemp, 10000, time(step), ierr); CHKERRQ(ierr)
                  call TSSolve(tsTemp, MEF90HeatXferCtx % temperature, time(step), ierr); CHKERRQ(ierr)
                  call TSGetConvergedReason(tsTemp, tsTempConvergedReason, ierr); CHKERRQ(ierr)
                  if (tsTempConvergedReason < 0) then
                     write (IOBuffer, 410) "temperature", tsTempConvergedReason
                     call PetscPrintf(MEF90Ctx % Comm, IOBuffer, ierr); CHKERRQ(ierr)
                  end if
                  call TSGetTime(tsTemp, t, ierr); CHKERRQ(ierr)
                  time(step) = t
               else
                  write (IOBuffer, *) 'TS exceeded analysis time. Skipping step\n'
                  call PetscPrintf(MEF90Ctx % Comm, IOBuffer, ierr); CHKERRQ(ierr)
               end if
            end if

            !!! Compute thermal energy
            call MEF90HeatXFerEnergy(MEF90HeatXferCtx % temperature, time(step), MEF90HeatXferCtx, thermalEnergySet, heatFluxWorkSet, ierr); CHKERRQ(ierr)
            call DMmeshGetLabelIdIS(MEF90HeatXferCtx % DM, 'Cell Sets', CellSetGlobalIS, ierr); CHKERRQ(ierr)
            call MEF90ISAllGatherMerge(MEF90Ctx % Comm, CellSetGlobalIS, ierr); CHKERRQ(ierr)
            call ISGetIndices(CellSetGlobalIS, setID, ierr); CHKERRQ(ierr)
            call PetscPrintf(MEF90Ctx % Comm, "\nThermal energies: \n", ierr); CHKERRQ(ierr)
            do set = 1, size(setID)
               write (IOBuffer, 101) setID(set), thermalEnergySet(set), heatFluxWorkSet(set), thermalEnergySet(set) - heatFluxWorkSet(set)
               call PetscPrintf(MEF90Ctx % Comm, IOBuffer, ierr); CHKERRQ(ierr)
            end do
            call ISRestoreIndices(CellSetGlobalIS, setID, ierr); CHKERRQ(ierr)
            call ISDestroy(CellSetGlobalIS, ierr); CHKERRQ(ierr)
            write (IOBuffer, 102) sum(thermalEnergySet), sum(heatFluxWorkSet), sum(thermalEnergySet) - sum(heatFluxWorkSet)
            call PetscPrintf(MEF90Ctx % Comm, IOBuffer, ierr); CHKERRQ(ierr)
            !!! Save results
            call MEF90HeatXferViewEXO(MEF90HeatXferCtx, step, ierr)
         case (MEF90HeatXfer_timeSteppingTypeNULL)
            continue
         case default
            write (IOBuffer, *) "Implemented HeatXfer mode: ", MEF90HeatXferGlobalOptions % timeSteppingType, "\n"
            call PetscPrintf(MEF90Ctx % Comm, IOBuffer, ierr); CHKERRQ(ierr)
            stop
         end select
         call PetscLogStagePop(ierr); CHKERRQ(ierr)

         !!! Solve for displacement and damage
         write (IOBuffer, 200) step, time(step)
         call PetscPrintf(MEF90Ctx % Comm, IOBuffer, ierr); CHKERRQ(ierr)
         damageMaxChange = 1.0d+20
         call MEF90DefMechUpdateDamageBounds(MEF90DefMechCtx, MEF90DefMechCtx % snesDamage, MEF90DefMechCtx % damage, ierr); CHKERRQ(ierr)

         !!! Update fields
         call PetscLogStagePush(logStageIO, ierr); CHKERRQ(ierr)
         call MEF90DefMechSetTransients(MEF90DefMechCtx, step, time(step), ierr)
         call MEF90DefMechUpdateboundaryDisplacement(MEF90DefMechCtx % displacement, MEF90DefMechCtx, ierr)
         call MEF90DefMechUpdateboundaryDamage(MEF90DefMechCtx % damage, MEF90DefMechCtx, ierr)
         call PetscLogStagePop(ierr); CHKERRQ(ierr)

         select case (MEF90DefMechGlobalOptions % timeSteppingType)
         case (MEF90DefMech_timeSteppingTypeQuasiStatic)
            !Call SNESSetLagPreconditioner(MEF90DefMechCtx%snesDamage,1,ierr);CHKERRQ(ierr)
            call SNESSetLagPreconditioner(MEF90DefMechCtx % snesDisp, 1, ierr); CHKERRQ(ierr)

            call VecCopy(MEF90DefMechCtx % CrackPressure, CrackPressureMask, ierr); CHKERRQ(ierr)

            AltMin: do AltMinIter = 1, MEF90DefMechGlobalOptions % maxit
               write (IObuffer, 208) AltMinIter
               call PetscPrintf(MEF90Ctx % Comm, IOBuffer, ierr); CHKERRQ(ierr)

               if (mod(AltMinIter - 1, MEF90DefMechGlobalOptions % PCLag) == 0) then
                  !Call SNESSetLagPreconditioner(MEF90DefMechCtx%snesDisp,-2,ierr);CHKERRQ(ierr)
                  call SNESSetLagPreconditioner(MEF90DefMechCtx % snesDamage, -2, ierr); CHKERRQ(ierr)
               end if

               call PetscLogStagePush(logStageDisplacement, ierr); CHKERRQ(ierr)
               !!! Solve for displacement, possibly using a Volume pressure equilibrium loop based on secant method

               !!! CrackPressure Block independent of alternate proj because exit if convergence
               if (any(ActivatedCrackPressureBlocksList)) then
                  !! Initialization Variable Secant Method
                  CrackPressureSave = [0.0_kr, 1.0_kr, 2.0_kr]
                  CrackVolumeSave = [0.0_kr, 1.0_kr, 2.0_kr]
                  if (step > 1) then
                     CrackVolumeSave = time(step - 1) * CrackVolumeSave
                  end if

                  SecantMthd: do CrackVolumeIter = 1, 5
                     I1 = mod(CrackVolumeIter + 2, 3) + 1
                     I2 = mod(CrackVolumeIter, 3) + 1
                     I3 = mod(CrackVolumeIter + 1, 3) + 1

                     if (CrackVolumeSave(I2) == CrackVolumeSave(I1)) then
                        write (*, *) "[WARNING] dividing by zero in the secant method, premature exit"
                        exit
                     end if

                     CrackPressureSave(I3) = CrackPressureSave(I2) - (CrackVolumeSave(I2) - time(step)) * (CrackPressureSave(I2) - CrackPressureSave(I1)) / (CrackVolumeSave(I2) - CrackVolumeSave(I1))

                     call VecCopy(CrackPressureMask, MEF90DefMechCtx % CrackPressure, ierr); CHKERRQ(ierr)
                     call VecScale(MEF90DefMechCtx % CrackPressure, CrackPressureSave(I3), ierr); CHKERRQ(ierr)

                     !!! Solve displacement SNES
                     call SNESSolve(MEF90DefMechCtx % snesDisp, PETSC_NULL_OBJECT, MEF90DefMechCtx % displacement, ierr); CHKERRQ(ierr)
                     call SNESGetConvergedReason(MEF90DefMechCtx % snesDisp, snesDispConvergedReason, ierr); CHKERRQ(ierr)
                     if (snesDispConvergedReason < 0) then
                        write (IOBuffer, 400) "displacement", snesDispConvergedReason
                        call PetscPrintf(MEF90Ctx % Comm, IOBuffer, ierr); CHKERRQ(ierr)
                     end if
                     !!!! Calculate the injected volume
                     CrackVolumeSet = 0.0_kr
                     call MEF90DefMechCrackVolume(MEF90DefMechCtx % displacement, MEF90DefMechCtx, CrackVolumeSet, ierr); CHKERRQ(ierr)
                     CrackVolumeSave(I3) = sum(CrackVolumeSet, MASK=ActivatedCrackPressureBlocksList)
                     InjectedVolumeRelativeError = abs(time(step) - CrackVolumeSave(I3)) / (1.0_kr + time(step))

                     !!!! Condition to exit loop for the secant method
                     if (InjectedVolumeRelativeError <= MEF90DefMechGlobalOptions % InjectedVolumeAtol) then
                        exit
                     end if
                  end do SecantMthd
                  write (IOBuffer, 302) CrackPressureSave(I3)
                  call PetscPrintf(MEF90Ctx % Comm, IOBuffer, ierr); CHKERRQ(ierr)
               else
                  call SNESSolve(MEF90DefMechCtx % snesDisp, PETSC_NULL_OBJECT, MEF90DefMechCtx % displacement, ierr); CHKERRQ(ierr)
                  call SNESGetConvergedReason(MEF90DefMechCtx % snesDisp, snesDispConvergedReason, ierr); CHKERRQ(ierr)
                  if (snesDispConvergedReason < 0) then
                     write (IOBuffer, 400) "displacement", snesDispConvergedReason
                     call PetscPrintf(MEF90Ctx % Comm, IOBuffer, ierr); CHKERRQ(ierr)
                  end if
                  !!! WorkControlled Block
                  if (any(ActivatedWorkControlledBlocksList)) then
                     forceWorkSet = 0.0_kr
                     call MEF90DefMechWork(MEF90DefMechCtx % displacement, MEF90DefMechCtx, forceWorkSet, ierr); CHKERRQ(ierr)
                     WorkControlled(step) = sum(forceWorkSet, MASK=ActivatedWorkControlledBlocksList)
                     WorkControlledRescaling = sqrt(time(step) / WorkControlled(step))
                     call VecScale(MEF90DefMechCtx % pressureForce, WorkControlledRescaling, ierr); CHKERRQ(ierr)
                     ErrorEstimationWorkControlled = ((abs(time(step) - WorkControlled(step))) / (1.0_kr + time(step)))
                  end if
               end if
               call PetscLogStagePop(ierr); CHKERRQ(ierr)

               call PetscLogStagePush(logStageDamage, ierr); CHKERRQ(ierr)
               call VecCopy(MEF90DefMechCtx % damage, damageAltMinOld, ierr); CHKERRQ(ierr)
               call SNESSolve(MEF90DefMechCtx % snesDamage, PETSC_NULL_OBJECT, MEF90DefMechCtx % damage, ierr); CHKERRQ(ierr)
               call SNESGetConvergedReason(MEF90DefMechCtx % snesDamage, snesDamageConvergedReason, ierr); CHKERRQ(ierr)
               if (snesDamageConvergedReason < 0) then
                  write (IOBuffer, 400) "damage field", snesDamageConvergedReason
                  call PetscPrintf(MEF90Ctx % Comm, IOBuffer, ierr); CHKERRQ(ierr)
               end if

               call VecMin(MEF90DefMechCtx % damage, PETSC_NULL_INTEGER, alphaMin, ierr); CHKERRQ(ierr)
               call VecMax(MEF90DefMechCtx % damage, PETSC_NULL_INTEGER, alphaMax, ierr); CHKERRQ(ierr)
               call VecAxPy(damageAltMinOld, -1.0_kr, MEF90DefMechCtx % damage, ierr); CHKERRQ(ierr)
               call VecNorm(damageAltMinOld, NORM_INFINITY, damageMaxChange, ierr); CHKERRQ(ierr)
               write (IOBuffer, 209) alphamin, alphamax, damageMaxChange
               call PetscPrintf(MEF90Ctx % Comm, IOBuffer, ierr); CHKERRQ(ierr)
               call PetscLogStagePop(ierr); CHKERRQ(ierr)

               !!! Conditions to exit the loop in alpha
               if (damageMaxChange <= MEF90DefMechGlobalOptions % damageATol) then
                  exit
               end if
               if (mod(AltMinIter, 25) == 0) then
                  call PetscLogStagePush(logStageIO, ierr); CHKERRQ(ierr)
                  call MEF90DefMechViewEXO(MEF90DefMechCtx, step, ierr)
                  call PetscLogStagePop(ierr); CHKERRQ(ierr)
               end if
            end do AltMin
            write (IOBuffer, 303) step, CrackVolumeSave(I3), CrackPressureSave(I3)
            call PetscViewerASCIIPrintf(pressureViewer, IOBuffer, ierr); CHKERRQ(ierr)
            call PetscViewerFlush(pressureViewer, ierr); CHKERRQ(ierr)

            if (AltMinIter == MEF90DefMechGlobalOptions % maxit) then
               write (IOBuffer, 412) MEF90DefMechGlobalOptions % maxit
               call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, ierr); CHKERRQ(ierr)
            end if

            !!! Compute energies
            elasticEnergySet = 0.0_kr
            forceWorkSet = 0.0_kr
            surfaceEnergySet = 0.0_kr
            cohesiveEnergySet = 0.0_kr
            CrackVolumeSet = 0.0_kr

            call PetscLogStagePush(logStageEnergy, ierr); CHKERRQ(ierr)
            call MEF90DefMechElasticEnergy(MEF90DefMechCtx % displacement, MEF90DefMechCtx, elasticEnergySet, ierr); CHKERRQ(ierr)
            call MEF90DefMechWork(MEF90DefMechCtx % displacement, MEF90DefMechCtx, forceWorkSet, ierr); CHKERRQ(ierr)
            call MEF90DefMechSurfaceEnergy(MEF90DefMechCtx % damage, MEF90DefMechCtx, surfaceEnergySet, ierr); CHKERRQ(ierr)
            call MEF90DefMechCohesiveEnergy(MEF90DefMechCtx % displacement, MEF90DefMechCtx, cohesiveEnergySet, ierr); CHKERRQ(ierr)

            elasticEnergy(step) = sum(elasticEnergySet)
            forceWork(step) = sum(forceWorkSet)
            surfaceEnergy(step) = sum(surfaceEnergySet)
            cohesiveEnergy(step) = sum(cohesiveEnergySet)
            totalMechanicalEnergy(step) = elasticEnergy(step) - forceWork(step) + cohesiveEnergy(step) + surfaceEnergy(step)
            !!!
            !!! Print and save energies
            !!!
            call DMmeshGetLabelIdIS(MEF90DefMechCtx % DMVect, 'Cell Sets', CellSetGlobalIS, ierr); CHKERRQ(ierr)
            call MEF90ISAllGatherMerge(MEF90Ctx % Comm, CellSetGlobalIS, ierr); CHKERRQ(ierr)
            call ISGetIndices(CellSetGlobalIS, setID, ierr); CHKERRQ(ierr)
            call PetscPrintf(MEF90Ctx % Comm, "\nMechanical energies: \n", ierr); CHKERRQ(ierr)
            do set = 1, size(setID)
               write (IOBuffer, 201) setID(set), elasticEnergySet(set), forceWorkSet(set), cohesiveEnergySet(set), surfaceEnergySet(set), elasticEnergySet(set) - forceWorkSet(set) + cohesiveEnergySet(set) + surfaceEnergySet(set)
               call PetscPrintf(MEF90Ctx % Comm, IOBuffer, ierr); CHKERRQ(ierr)

               write (IOBuffer, 500) step, time(step), elasticEnergySet(set), forceWorkSet(set), cohesiveEnergySet(set), surfaceEnergySet(set), elasticEnergySet(set) - forceWorkSet(set) + cohesiveEnergySet(set) + surfaceEnergySet(set)
               call PetscViewerASCIIPrintf(MEF90DefMechCtx % setEnergyViewer(set), IOBuffer, ierr); CHKERRQ(ierr)
               call PetscViewerFlush(MEF90DefMechCtx % setEnergyViewer(set), ierr); CHKERRQ(ierr)
            end do
            call ISRestoreIndices(CellSetGlobalIS, setID, ierr); CHKERRQ(ierr)
            call ISDestroy(CellSetGlobalIS, ierr); CHKERRQ(ierr)
            write (IOBuffer, 202) elasticEnergy(step), forceWork(step), cohesiveEnergy(step), surfaceEnergy(step), totalMechanicalEnergy(step)
            call PetscPrintf(MEF90Ctx % Comm, IOBuffer, ierr); CHKERRQ(ierr)
            write (IOBuffer, 500) step, time(step), elasticEnergy(step), cohesiveEnergy(step), forceWork(step), surfaceEnergy(step), totalMechanicalEnergy(step)
            call PetscViewerASCIIPrintf(MEF90DefMechCtx % globalEnergyViewer, IOBuffer, ierr); CHKERRQ(ierr)
            call PetscViewerFlush(MEF90DefMechCtx % globalEnergyViewer, ierr); CHKERRQ(ierr)
            call PetscLogStagePop(ierr); CHKERRQ(ierr)

         case (MEF90DefMech_timeSteppingTypeNULL)
            continue
         case default
            write (IOBuffer, *) "Implemented DefMech time stepping: ", MEF90DefMechGlobalOptions % timeSteppingType, "\n"
            call PetscPrintf(MEF90Ctx % Comm, IOBuffer, ierr); CHKERRQ(ierr)
            stop
         end select

         !!!
         !!! Save results and boundary Values
         !!!
         call PetscLogStagePush(logStageIO, ierr); CHKERRQ(ierr)
         if (MEF90DefMechGlobalOptions % stressOffset > 0) then
            call MEF90DefMechStress(MEF90DefMechCtx % displacement, MEF90DefMechCtx, MEF90DefMechCtx % stress, ierr)
         end if

         call MEF90DefMechViewEXO(MEF90DefMechCtx, step, ierr)
         call PetscLogStagePop(ierr); CHKERRQ(ierr)

         !!!
         !!! Save performance log file
         !!!
         call PetscViewerASCIIOpen(MEF90Ctx % comm, trim(MEF90FilePrefix(MEF90Ctx % resultFile))//'.log', logViewer, ierr); CHKERRQ(ierr)
         call PetscLogView(logViewer, ierr); CHKERRQ(ierr)
         call PetscViewerFlush(logViewer, ierr); CHKERRQ(ierr)

         if (step == MEF90GlobalOptions % timeNumStep) then
            exit
         else
            step = step + 1
         end if
      end do MainloopQS
   end if
   write (IOBuffer, *) 'Total number of alternate minimizations:', AltMinIter, '\n'
   call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, ierr); CHKERRQ(ierr)

   !!! Clean up and exit nicely
   select case (MEF90DefMechGlobalOptions % timeSteppingType)
   case (MEF90DefMech_timeSteppingTypeQuasiStatic)
      call SNESDestroy(MEF90DefMechCtx % snesDisp, ierr); CHKERRQ(ierr)
      call SNESDestroy(MEF90DefMechCtx % snesDamage, ierr); CHKERRQ(ierr)
      call VecDestroy(residualDisp, ierr); CHKERRQ(ierr)
      call VecDestroy(residualDamage, ierr); CHKERRQ(ierr)
   end select

   select case (MEF90HeatXferGlobalOptions % timeSteppingType)
   case (MEF90HeatXfer_timeSteppingTypeSteadyState)
      call SNESDestroy(snesTemp, ierr); CHKERRQ(ierr)
   case (MEF90HeatXfer_timeSteppingTypeTransient)
      call TSDestroy(tsTemp, ierr); CHKERRQ(ierr)
   end select

   call MEF90DefMechCtxDestroyVectors(MEF90DefMechCtx, ierr)
   nullify (MEF90HeatXferCtx % temperature)
   call MEF90HeatXferCtxDestroyVectors(MEF90HeatXferCtx, ierr)
   call VecDestroy(damageAltMinOld, ierr); CHKERRQ(ierr)
   call VecDestroy(CrackPressureMask, ierr); CHKERRQ(ierr)

   call DMDestroy(Mesh, ierr); CHKERRQ(ierr)

   deallocate (elasticEnergySet)
   deallocate (surfaceEnergySet)
   deallocate (forceWorkSet)
   deallocate (cohesiveEnergySet)
   deallocate (thermalEnergySet)
   deallocate (heatFluxWorkSet)
   deallocate (ActivatedCrackPressureBlocksList)

   deallocate (elasticEnergy)
   deallocate (forceWork)
   deallocate (cohesiveEnergy)
   deallocate (surfaceEnergy)
   deallocate (totalMechanicalEnergy)

   call MEF90DefMechCtxDestroy(MEF90DefMechCtx, ierr); CHKERRQ(ierr)
   call MEF90HeatXferCtxDestroy(MEF90HeatXferCtx, ierr); CHKERRQ(ierr)
   call MEF90CtxCloseEXO(MEF90Ctx, ierr)

   call PetscLogView(logViewer, ierr); CHKERRQ(ierr)
   call PetscViewerDestroy(logViewer, ierr); CHKERRQ(ierr)
   call PetscViewerDestroy(pressureViewer, ierr); CHKERRQ(ierr)
   call MEF90CtxDestroy(MEF90Ctx, ierr)
   call MEF90Finalize(ierr)
   call PetscFinalize(ierr)
100 format("\nHeat transfer: solving steady state step ", I4, ", t=", ES12.5, "\n")
101 format("cell set ", I4, " thermal energy: ", ES12.5, " fluxes work: ", ES12.5, " total: ", ES12.5, "\n")
102 format("======= Total thermal energy: ", ES12.5, " fluxes work: ", ES12.5, " total: ", ES12.5, "\n")
110 format("\nHeat transfer: step ", I4, ", until t=", ES12.5, "\n")
200 format("\nMechanics: step ", I4, ", t=", ES12.5, "\n")
201 format("cell set ", I4, "  elastic energy: ", ES12.5, " work: ", ES12.5, " cohesive: ", ES12.5, " surface: ", ES12.5, " total: ", ES12.5, "\n")
202 format("======= Total: elastic energy: ", ES12.5, " work: ", ES12.5, " cohesive: ", ES12.5, " surface: ", ES12.5, " total: ", ES12.5, "\n")
208 format("   Alt. Min. step ", I5)
209 format(" alpha min / max", ES12.5, " / ", ES12.5, ", max change ", ES12.5, "\n")
301 format("   CrackPressure Iter", I5, " error ", ES12.5, "\n")
302 format(" crack pressure ", ES12.5)
303 format(I6, 2(ES16.5), "\n")
308 format("   Alt. Proj. step ", I5, " ")
400 format(" [ERROR]: ", A, " SNESSolve failed with SNESConvergedReason ", I2, ". \n Check http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESConvergedReason.html for error code meaning.\n")
410 format(" [ERROR]: ", A, " TSSolve failed with TSConvergedReason ", I2, ". \n Check http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESConvergedReason.html for error code meaning.\n")
412 format(" [ERROR]: Alternate minimizations algorithm did NOT converge in ", I5, "iterations.\n")
500 format(I6, 6(ES16.5), "\n")
end program vDefHF
