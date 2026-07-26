#include "../MEF90/mef90.inc"
program CoupledPlasticityDamage
#include <finclude/petscdef.h>

   use petsc
   use m_MEF90
   use m_MEF90_DefMech_class
   use m_MEF90_DefMech
   use m_MEF90_HeatXfer_class
   use m_MEF90_HeatXfer
   implicit none(type, external)

   PetscErrorCode                                     :: ierr
   type(MEF90Ctx_Type), target                         :: MEF90Ctx
   type(MEF90CtxGlobalOptions_Type), pointer           :: MEF90GlobalOptions

   !!! Defect mechanics contexts
   type(MEF90DefMech_Type), target                     :: MEF90DefMechCtx
   type(MEF90DefMechGlobalOptions_Type)                :: MEF90DefMechGlobalOptions

   !!! HeatXfer contexts
   type(MEF90HeatXfer_Type), target                    :: MEF90HeatXferCtx
   type(MEF90HeatXferGlobalOptions_Type)               :: MEF90HeatXferGlobalOptions

   type(DM), target                                    :: Mesh
   type(IS)                                           :: setIS, CellSetGlobalIS
   PetscInt, dimension(:), pointer                      :: setID
   PetscInt                                           :: numCellSet
   PetscInt                                           :: numset, set
   PetscReal, dimension(:), pointer                     :: time
   PetscReal, dimension(:), pointer                     :: thermalEnergySet, heatFluxWorkSet
   PetscReal, dimension(:), pointer                     :: elasticEnergySet, surfaceEnergySet, forceWorkSet, cohesiveEnergySet
   PetscReal, dimension(:), pointer                     :: elasticEnergy, surfaceEnergy, forceWork, cohesiveEnergy, totalMechanicalEnergy

   type(SNES)                                         :: snesDisp
   SNESConvergedReason                                :: snesDispConvergedReason
   type(Vec)                                          :: residualDisp, displacementOld
   PetscReal                                          :: displacementMaxChange
   type(SNES)                                         :: snesDamage
   SNESConvergedReason                                :: snesDamageConvergedReason
   type(Vec)                                          :: residualDamage, damageOld
   PetscInt                                           :: AltMinIter
   PetscReal                                          :: damageMaxChange

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
   type(PetscViewer)                                  :: logViewer
   integer                                            :: numfield

   integer                                            :: step
   PetscInt                                           :: dim

   PetscReal                                          :: alphaMaxChange, alphaMin, alphaMax

   !!! plastic variables
   PetscReal, dimension(:), pointer                     :: plasticDissipationSet
   PetscReal, dimension(:), pointer                     :: plasticDissipation
   type(Vec)                                          :: plasticStrainOld
   PetscInt                                           :: InnerLoopIter
   PetscReal                                          :: PlasticStrainMaxChange, RelativeAbsoluteplasticStrainATol
   type(Vec)                                          :: plasticstrainerror
   type(Vec)                                          :: plasticStrainPrevious
   type(Vec)                                          :: localVec

   !!! cumulatedDissipatedPlasticEnergy
   type(Vec)                                          :: cumulatedDissipatedPlasticEnergyOld
   type(Vec)                                          :: cumulatedDissipatedPlasticEnergyVariation

   type(MEF90DefMechCellSetOptions_Type)              :: cellSetOptions

   !!! Initialize MEF90
   call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
   call MEF90Initialize(PETSC_COMM_WORLD, ierr)
   call PetscPrintf(PETSC_COMM_WORLD, " # vDefUpa: numerical implementation of variational models of Ductile Defect Mechanics\n", ierr); CHKERRQ(ierr)

   !!! Get all MEF90-wide options
   call MEF90CtxCreate(PETSC_COMM_WORLD, MEF90Ctx, "", ierr); CHKERRQ(ierr)
   call MEF90Ctx%setFromOptions(ierr); CHKERRQ(ierr)
   MEF90GlobalOptions => MEF90Ctx%globalOptions

   !!! Get DM from mesh
   call MEF90CtxGetDMMeshEXO(MEF90Ctx, Mesh, ierr); CHKERRQ(ierr)
   call DMMeshGetDimension(Mesh, dim, ierr); CHKERRQ(ierr)
   call DMMeshSetMaxDof(Mesh, dim, ierr); CHKERRQ(ierr)
   call DMSetBlockSize(Mesh, dim, ierr); CHKERRQ(ierr)

   !!! Open output file
   call MEF90CtxOpenEXO(MEF90Ctx, Mesh, ierr)

   !!! Create DefMech context, get all DefMech options
   call MEF90DefMechCreate(MEF90DefMechCtx, Mesh, MEF90Ctx, "", ierr); CHKERRQ(ierr)
   call MEF90DefMechCtx%setFromOptions(ierr); CHKERRQ(ierr)
   MEF90DefMechGlobalOptions = MEF90DefMechCtx%globalOptions

   !!! Create HeatXfer context, get all HeatXfer options
   call MEF90HeatXferCreate(MEF90HeatXferCtx, Mesh, MEF90Ctx, "", ierr); CHKERRQ(ierr)
   !!! vDef does not export the temperature: the DefMech context owns that field
   MEF90HeatXferCtx%globalOptions%temperatureExport = PETSC_FALSE
   call MEF90HeatXferCtx%setFromOptions(ierr); CHKERRQ(ierr)
   MEF90HeatXferGlobalOptions = MEF90HeatXferCtx%globalOptions

   !!! Get material properties bags
   if (dim == 2) then

   !!! Create time array from global options
   call MEF90CtxGetTime(MEF90Ctx, time, ierr)

   !!! Create sections, vectors, and solvers for DefMech Context
   call MEF90DefMechCtxSetSections(MEF90DefMechCtx, ierr)
   call MEF90DefMechCreateVectors(MEF90DefMechCtx, ierr)
   call VecDuplicate(MEF90DefMechCtx%damage, damageOld, ierr); CHKERRQ(ierr)
   call VecDuplicate(MEF90DefMechCtx%displacement, residualDisp, ierr); CHKERRQ(ierr)
   call PetscObjectSetName(residualDisp, "residualDisp", ierr); CHKERRQ(ierr)
   call VecDuplicate(MEF90DefMechCtx%displacement, displacementOld, ierr); CHKERRQ(ierr)
   call PetscObjectSetName(displacementOld, "displacementOld", ierr); CHKERRQ(ierr)
   call MEF90DefMechCreateSNESDisplacement(MEF90DefMechCtx, snesDisp, residualDisp, ierr)
   call VecDuplicate(MEF90DefMechCtx%damage, residualDamage, ierr); CHKERRQ(ierr)
   call PetscObjectSetName(residualDamage, "residualDamage", ierr); CHKERRQ(ierr)
   call MEF90DefMechCreateSNESDamage(MEF90DefMechCtx, snesDamage, residualDamage, ierr)

   !!!cumulatedDissipatedPlasticEnergy Vectors
   call VecDuplicate(MEF90DefMechCtx%cumulatedPlasticDissipation, cumulatedDissipatedPlasticEnergyOld, ierr); CHKERRQ(ierr)
   call VecDuplicate(MEF90DefMechCtx%cumulatedPlasticDissipation, cumulatedDissipatedPlasticEnergyVariation, ierr); CHKERRQ(ierr)
   call VecCopy(MEF90DefMechCtx%cumulatedPlasticDissipation, cumulatedDissipatedPlasticEnergyOld, ierr); CHKERRQ(ierr)
   deallocate (MEF90DefMechCtx%temperature)

   call VecDuplicate(MEF90DefMechCtx%plasticStrain, plasticStrainOld, ierr); CHKERRQ(ierr)
   call VecDuplicate(MEF90DefMechCtx%plasticStrain, plasticStrainPrevious, ierr); CHKERRQ(ierr)

   !!! Create sections, vectors, and solvers for HeatXfer Context
   if (MEF90HeatXferGlobalOptions%timeSteppingType /= MEF90HeatXfer_timeSteppingTypeNULL) then
      call MEF90HeatXferCtxSetSections(MEF90HeatXferCtx, ierr)
      call MEF90HeatXferCreateVectors(MEF90HeatXferCtx, ierr)
      call VecDuplicate(MEF90HeatXferCtx%temperature, residualTemp, ierr); CHKERRQ(ierr)
      call PetscObjectSetName(residualTemp, "residualTemp", ierr); CHKERRQ(ierr)
      select case (MEF90HeatXferGlobalOptions%timeSteppingType)
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
      MEF90DefMechCtx%temperature => MEF90HeatXferCtx%temperature
   end if
   !!!
   !!! Allocate array of works and energies
   !!!
   call MEF90DMGetNumSets(MEF90DefMechCtx%megaDM, MEF90CellSetLabelName, numCellSet, ierr); CHKERRQ(ierr)
   allocate (elasticEnergySet(numCellSet))
   elasticEnergySet = 0.0_kr
   allocate (surfaceEnergySet(numCellSet))
   surfaceEnergySet = 0.0_kr
   allocate (forceWorkSet(numCellSet))
   forceWorkSet = 0.0_kr
   allocate (cohesiveEnergySet(numCellSet))
   cohesiveEnergySet = 0.0_kr
   allocate (thermalEnergySet(numCellSet))
   thermalEnergySet = 0.0_kr
   allocate (heatFluxWorkSet(numCellSet))
   heatFluxWorkSet = 0.0_kr
   allocate (plasticDissipationSet(numCellSet))
   plasticDissipationSet = 0.0_kr

   allocate (elasticEnergy(MEF90GlobalOptions%timeNumStep))
   elasticEnergy = 0.0_kr
   allocate (surfaceEnergy(MEF90GlobalOptions%timeNumStep))
   surfaceEnergy = 0.0_kr
   allocate (forceWork(MEF90GlobalOptions%timeNumStep))
   forceWork = 0.0_kr
   allocate (cohesiveEnergy(MEF90GlobalOptions%timeNumStep))
   cohesiveEnergy = 0.0_kr
   allocate (totalMechanicalEnergy(MEF90GlobalOptions%timeNumStep))
   totalMechanicalEnergy = 0.0_kr
   allocate (plasticDissipation(MEF90GlobalOptions%timeNumStep))
   plasticDissipation = 0.0_kr

   !!!
   !!! Try to figure out if the file was formatted
   !!!
   if (MEF90Ctx%rank == 0) then
      call EXGVP(MEF90Ctx%fileExoUnit, "N", numfield, ierr)
   end if
   call MPI_Bcast(numfield, 1, MPIU_INTEGER, 0, MEF90Ctx%comm, ierr)
   if (numfield == 0) then
      call MEF90DefMechFormatEXO(MEF90DefMechCtx, time, ierr)
   end if

   !!! Actual computations / time stepping
   !!!
   if (.not. MEF90GlobalOptions%dryrun) then
      if (MEF90GlobalOptions%timeSkip > 0) then
         !!! Restore state from file.
         call DMGetLocalVector(MEF90DefMechCtx%DMScal, localVec, ierr); CHKERRQ(ierr)
         call VecLoadExodusVertex(MEF90DefMechCtx%DMScal, localVec, MEF90DefMechCtx%MEF90Ctx%IOcomm, &
                                  MEF90DefMechCtx%MEF90Ctx%fileExoUnit, MEF90GlobalOptions%timeSkip, MEF90DefMechGlobalOptions%damageOffset, ierr); CHKERRQ(ierr)
         call DMLocalToGlobalBegin(MEF90DefMechCtx%DMScal, localVec, INSERT_VALUES, MEF90DefMechCtx%damage, ierr); CHKERRQ(ierr)
         call DMLocalToGlobalEnd(MEF90DefMechCtx%DMScal, localVec, INSERT_VALUES, MEF90DefMechCtx%damage, ierr); CHKERRQ(ierr)
         call VecCopy(MEF90DefMechCtx%damage, damageOld, ierr); CHKERRQ(ierr)
         call DMRestoreLocalVector(MEF90DefMechCtx%DMScal, localVec, ierr); CHKERRQ(ierr)

         if (MEF90DefMechGlobalOptions%cumulatedPlasticDissipationOffset > 0) then
            call DMGetLocalVector(MEF90DefMechCtx%cellDMScal, localVec, ierr); CHKERRQ(ierr)
            call VecLoadExodusCell(MEF90DefMechCtx%cellDMScal, localVec, MEF90DefMechCtx%MEF90Ctx%IOcomm, &
                                   MEF90DefMechCtx%MEF90Ctx%fileExoUnit, MEF90GlobalOptions%timeSkip, MEF90DefMechGlobalOptions%cumulatedPlasticDissipationOffset, ierr); CHKERRQ(ierr)
            call DMLocalToGlobalBegin(MEF90DefMechCtx%cellDMScal, localVec, INSERT_VALUES, MEF90DefMechCtx%cumulatedPlasticDissipation, ierr); CHKERRQ(ierr)
            call DMLocalToGlobalEnd(MEF90DefMechCtx%cellDMScal, localVec, INSERT_VALUES, MEF90DefMechCtx%cumulatedPlasticDissipation, ierr); CHKERRQ(ierr)
            call VecCopy(MEF90DefMechCtx%cumulatedPlasticDissipation, cumulatedDissipatedPlasticEnergyOld, ierr); CHKERRQ(ierr)
            call DMRestoreLocalVector(MEF90DefMechCtx%cellDMScal, localVec, ierr); CHKERRQ(ierr)
         end if

         if (MEF90DefMechGlobalOptions%plasticStrainOffset > 0) then
            call DMGetLocalVector(MEF90DefMechCtx%cellDMMatS, localVec, ierr); CHKERRQ(ierr)
            call VecLoadExodusCell(MEF90DefMechCtx%cellDMMatS, localVec, MEF90DefMechCtx%MEF90Ctx%IOcomm, &
                                   MEF90DefMechCtx%MEF90Ctx%fileExoUnit, MEF90GlobalOptions%timeSkip, MEF90DefMechGlobalOptions%plasticStrainOffset, ierr); CHKERRQ(ierr)
            call DMLocalToGlobalBegin(MEF90DefMechCtx%cellDMMatS, localVec, INSERT_VALUES, MEF90DefMechCtx%plasticStrain, ierr); CHKERRQ(ierr)
            call DMLocalToGlobalEnd(MEF90DefMechCtx%cellDMMatS, localVec, INSERT_VALUES, MEF90DefMechCtx%plasticStrain, ierr); CHKERRQ(ierr)
            call VecCopy(MEF90DefMechCtx%plasticStrain, plasticStrainOld, ierr); CHKERRQ(ierr)
            call DMRestoreLocalVector(MEF90DefMechCtx%cellDMMatS, localVec, ierr); CHKERRQ(ierr)
         end if

         call DMGetLocalVector(MEF90DefMechCtx%DMVect, localVec, ierr); CHKERRQ(ierr)
         call VecLoadExodusVertex(MEF90DefMechCtx%DMVect, localVec, MEF90DefMechCtx%MEF90Ctx%IOcomm, &
                                  MEF90DefMechCtx%MEF90Ctx%fileExoUnit, MEF90GlobalOptions%timeSkip, MEF90DefMechGlobalOptions%displacementOffset, ierr); CHKERRQ(ierr)
         call DMLocalToGlobalBegin(MEF90DefMechCtx%DMVect, localVec, INSERT_VALUES, MEF90DefMechCtx%Displacement, ierr); CHKERRQ(ierr)
         call DMLocalToGlobalEnd(MEF90DefMechCtx%DMVect, localVec, INSERT_VALUES, MEF90DefMechCtx%Displacement, ierr); CHKERRQ(ierr)
         call DMRestoreLocalVector(MEF90DefMechCtx%DMVect, localVec, ierr); CHKERRQ(ierr)
      end if
      step = MEF90GlobalOptions%timeSkip + 1
      mainloopQS: do
         !!! Solve for temperature
         select case (MEF90HeatXferGlobalOptions%timeSteppingType)
         case (MEF90HeatXfer_timeSteppingTypeSteadyState)
            write (IOBuffer, 100) step, time(step)
            call PetscPrintf(MEF90Ctx%comm, IOBuffer, ierr); CHKERRQ(ierr)

            !!! Update fields
            call MEF90HeatXferSetTransients(MEF90HeatXferCtx, step, time(step), ierr)
            !!! Solve SNES
            call MEF90HeatXferUpdateboundaryTemperature(MEF90HeatXferCtx%temperature, MEF90HeatXferCtx, ierr); 
            call SNESSolve(snesTemp, PETSC_NULL_OBJECT, MEF90HeatXferCtx%temperature, ierr); CHKERRQ(ierr)
            call SNESGetConvergedReason(snesTemp, snesTempConvergedReason, ierr); CHKERRQ(ierr)
            if (snesTempConvergedReason < 0) then
               write (IOBuffer, 400) "temperature", snesTempConvergedReason
               call PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr); CHKERRQ(ierr)
            end if

            !!! Compute thermal energy
            call MEF90HeatXFerEnergy(MEF90HeatXferCtx%temperature, time(step), MEF90HeatXferCtx, thermalEnergySet, heatFluxWorkSet, ierr); CHKERRQ(ierr)
            call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM, 'Cell Sets', CellSetGlobalIS, ierr); CHKERRQ(ierr)
            call MEF90ISAllGatherMerge(MEF90Ctx%Comm, CellSetGlobalIS, ierr); CHKERRQ(ierr)
            call ISGetIndices(CellSetGlobalIS, setID, ierr); CHKERRQ(ierr)
            call PetscPrintf(MEF90Ctx%Comm, "\nThermal energies: \n", ierr); CHKERRQ(ierr)
            do set = 1, size(setID)
               write (IOBuffer, 101) setID(set), thermalEnergySet(set), heatFluxWorkSet(set), thermalEnergySet(set) - heatFluxWorkSet(set)
               call PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr); CHKERRQ(ierr)
            end do
            call ISRestoreIndices(CellSetGlobalIS, setID, ierr); CHKERRQ(ierr)
            call ISDestroy(CellSetGlobalIS, ierr); CHKERRQ(ierr)
            write (IOBuffer, 102) sum(thermalEnergySet), sum(heatFluxWorkSet), sum(thermalEnergySet) - sum(heatFluxWorkSet)
            call PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr); CHKERRQ(ierr)

            !!! Save results
            call MEF90HeatXferViewEXO(MEF90HeatXferCtx, step, ierr)
         case (MEF90HeatXfer_timeSteppingTypeTransient)
            if (step > 1) then
               write (IOBuffer, 110) step, time(step)
               call PetscPrintf(MEF90Ctx%comm, IOBuffer, ierr); CHKERRQ(ierr)
               !!! Update fields
               call MEF90HeatXferSetTransients(MEF90HeatXferCtx, step, time(step), ierr)
               call MEF90HeatXferUpdateboundaryTemperature(MEF90HeatXferCtx%temperature, MEF90HeatXferCtx, ierr); 
               !!! Make sure TS does not overstep
               call TSGetTime(tsTemp, t, ierr); CHKERRQ(ierr)
               if (t < time(step)) then
                  call TSAdaptSetStepLimits(tsAdaptTemp, PETSC_DECIDE, (time(step) - time) / 2.0_kr, ierr); CHKERRQ(ierr)
                  !!! Something is up here.
                  !!! replacing the constant 10000 with a variable leads to divergence of TSAdapt
                  !!! when using gcc
                  call TSSetDuration(tsTemp, 10000, time(step), ierr); CHKERRQ(ierr)
                  call TSSolve(tsTemp, MEF90HeatXferCtx%temperature, time(step), ierr); CHKERRQ(ierr)
                  call TSGetConvergedReason(tsTemp, tsTempConvergedReason, ierr); CHKERRQ(ierr)
                  if (tsTempConvergedReason < 0) then
                     write (IOBuffer, 410) "temperature", tsTempConvergedReason
                     call PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr); CHKERRQ(ierr)
                  end if
                  call TSGetTime(tsTemp, t, ierr); CHKERRQ(ierr)
                  time(step) = t
               else
                  write (IOBuffer, *) 'TS exceeded analysis time. Skipping step\n'
                  call PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr); CHKERRQ(ierr)
               end if
            end if

            !!! Compute thermal energy
            call MEF90HeatXFerEnergy(MEF90HeatXferCtx%temperature, time(step), MEF90HeatXferCtx, thermalEnergySet, heatFluxWorkSet, ierr); CHKERRQ(ierr)
            call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM, 'Cell Sets', CellSetGlobalIS, ierr); CHKERRQ(ierr)
            call MEF90ISAllGatherMerge(MEF90Ctx%Comm, CellSetGlobalIS, ierr); CHKERRQ(ierr)
            call ISGetIndices(CellSetGlobalIS, setID, ierr); CHKERRQ(ierr)
            call PetscPrintf(MEF90Ctx%Comm, "\nThermal energies: \n", ierr); CHKERRQ(ierr)
            do set = 1, size(setID)
               write (IOBuffer, 101) setID(set), thermalEnergySet(set), heatFluxWorkSet(set), thermalEnergySet(set) - heatFluxWorkSet(set)
               call PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr); CHKERRQ(ierr)
            end do
            call ISRestoreIndices(CellSetGlobalIS, setID, ierr); CHKERRQ(ierr)
            call ISDestroy(CellSetGlobalIS, ierr); CHKERRQ(ierr)
            write (IOBuffer, 102) sum(thermalEnergySet), sum(heatFluxWorkSet), sum(thermalEnergySet) - sum(heatFluxWorkSet)
            call PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr); CHKERRQ(ierr)
            !!! Save results
            call MEF90HeatXferViewEXO(MEF90HeatXferCtx, step, ierr)
         case (MEF90HeatXfer_timeSteppingTypeNULL)
            continue
         case default
            write (IOBuffer, *) "Implemented HeatXfer mode: ", MEF90HeatXferGlobalOptions%timeSteppingType, "\n"
            call PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr); CHKERRQ(ierr)
            stop
         end select

         !!! Solve for displacement and damage
         select case (MEF90DefMechGlobalOptions%timeSteppingType)
         case (MEF90DefMech_timeSTeppingTypeQuasiStatic)
            write (IOBuffer, 200) step, time(step)
            call PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr); CHKERRQ(ierr)
            damageMaxChange = 1.0d+20

            call MEF90DefMechUpdateDamageBounds(MEF90DefMechCtx, snesDamage, MEF90DefMechCtx%damage, ierr); CHKERRQ(ierr)

            !!! Update fields
            call MEF90DefMechSetTransients(MEF90DefMechCtx, step, time(step), ierr)
            call MEF90DefMechUpdateboundaryDisplacement(MEF90DefMechCtx%displacement, MEF90DefMechCtx, ierr)
            call MEF90DefMechUpdateboundaryDamage(MEF90DefMechCtx%damage, MEF90DefMechCtx, ierr)

            call SNESSetLagPreconditioner(snesDisp, 1, ierr); CHKERRQ(ierr)
            AltMin: do AltMinIter = 1, MEF90DefMechGlobalOptions%maxit
               write (IObuffer, 208) AltMinIter
               call PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr); CHKERRQ(ierr)

               if (mod(AltMinIter - 1, MEF90DefMechGlobalOptions%PCLag) == 0) then
                  call SNESSetLagPreconditioner(snesDamage, -2, ierr); CHKERRQ(ierr)
               end if
               call VecCopy(MEF90DefMechCtx%displacement, displacementOld, ierr); CHKERRQ(ierr)
               call SNESSolve(snesDisp, PETSC_NULL_OBJECT, MEF90DefMechCtx%displacement, ierr); CHKERRQ(ierr)
               call SNESGetConvergedReason(snesDisp, snesDispConvergedReason, ierr); CHKERRQ(ierr)
               if (snesDispConvergedReason < 0) then
                  write (IOBuffer, 400) "displacement", snesDispConvergedReason
                  call PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr); CHKERRQ(ierr)
               end if

               call SNESSetLagPreconditioner(snesDamage, 1, ierr); CHKERRQ(ierr)
               call SNESSetLagJacobian(snesDamage, 1, ierr); CHKERRQ(ierr)
               InnerLoop: do InnerLoopIter = 1, MEF90DefMechGlobalOptions%maxit
                  write (IObuffer, 308) InnerLoopIter
                  !!! Since u does not change in this loop, the Jacobian for alpha is constant
                  !!! so we only evaluate at the first iteration
                  if (InnerLoopIter > 1) then
                     call SNESSetLagJacobian(snesDamage, -1, ierr); CHKERRQ(ierr)
                  end if
                  if (mod(InnerLoopIter - 1, MEF90DefMechGlobalOptions%PCLag) == 0) then
                     call SNESSetLagPreconditioner(snesDisp, -2, ierr); CHKERRQ(ierr)
                  end if

                  !!! Solve PlasticProjection
                  !!! Save the PlasticStrainPrevious, add damage in DefMechPlasticStrainUpdate, error on the norm PlasticStrain
                  !!! calculate cumulatedDissipatedPlasticEnergy =  cumulatedDissipatedPlasticEnergyPrevious + cumulatedDissipatedPlasticEnergyVariation
                  !!! Absolute/Relative error
                  !!! || p_i - p_{i-1} ||_L_inifnity / (1+ || p_i ||_L_inifnity) < tolerance
                  call VecCopy(MEF90DefMechCtx%plasticStrain, plasticStrainPrevious, ierr); CHKERRQ(ierr)
                  call MEF90DefMechPlasticStrainUpdate(MEF90DefMechCtx, MEF90DefMechCtx%plasticStrain, MEF90DefMechCtx%displacement, plasticStrainOld, plasticStrainPrevious, cumulatedDissipatedPlasticEnergyVariation, cumulatedDissipatedPlasticEnergyOld, ierr); CHKERRQ(ierr)
                  call VecAxPy(plasticStrainPrevious, -1.0_kr, MEF90DefMechCtx%plasticStrain, ierr); CHKERRQ(ierr)
                  call VecNorm(plasticStrainPrevious, NORM_INFINITY, PlasticStrainMaxChange, ierr); CHKERRQ(ierr)
                  write (IOBuffer, 210) InnerLoopIter, PlasticStrainMaxChange
                  call PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr); CHKERRQ(ierr)
                  !!! not sure about this...
                  call VecNorm(MEF90DefMechCtx%plasticStrain, NORM_INFINITY, RelativeAbsoluteplasticStrainATol, ierr); CHKERRQ(ierr)
                  RelativeAbsoluteplasticStrainATol = (1.0_kr + RelativeAbsoluteplasticStrainATol) * MEF90DefMechGlobalOptions%plasticStrainATol

                  call VecWAXPY(MEF90DefMechCtx%cumulatedPlasticDissipation, 1.0_kr, cumulatedDissipatedPlasticEnergyOld, cumulatedDissipatedPlasticEnergyVariation, ierr); CHKERRQ(ierr)

                  !!! Minimization with respect to alpha
                  call VecCopy(MEF90DefMechCtx%damage, damageOld, ierr); CHKERRQ(ierr)
                  call SNESSolve(snesDamage, PETSC_NULL_OBJECT, MEF90DefMechCtx%damage, ierr); CHKERRQ(ierr)
                  call SNESGetConvergedReason(snesDamage, snesDamageConvergedReason, ierr); CHKERRQ(ierr)
                  if (snesDamageConvergedReason < 0) then
                     write (IOBuffer, 400) "damage", snesDamageConvergedReason
                     call PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr); CHKERRQ(ierr)
                  end if
                  call VecMin(MEF90DefMechCtx%damage, PETSC_NULL_INTEGER, alphaMin, ierr); CHKERRQ(ierr)
                  call VecMax(MEF90DefMechCtx%damage, PETSC_NULL_INTEGER, alphaMax, ierr); CHKERRQ(ierr)
                  call VecAxPy(damageOld, -1.0_kr, MEF90DefMechCtx%damage, ierr); CHKERRQ(ierr)
                  call VecNorm(damageOld, NORM_INFINITY, damageMaxChange, ierr); CHKERRQ(ierr)
                  write (IOBuffer, 209) alphamin, alphamax, damageMaxChange
                  call PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr); CHKERRQ(ierr)
                  call VecNorm(damageOld, NORM_INFINITY, damageMaxChange, ierr); CHKERRQ(ierr)

                  if ((PlasticStrainMaxChange <= RelativeAbsoluteplasticStrainATol) .and. (damageMaxChange <= MEF90DefMechGlobalOptions%damageATol)) then
                     exit
                  end if
               end do InnerLoop

               !!! Test on displacement change. It would probably make more sense to check the change on constrain.
               call VecAxPy(displacementOld, -1.0_kr, MEF90DefMechCtx%displacement, ierr); CHKERRQ(ierr)
               call VecNorm(displacementOld, NORM_INFINITY, displacementMaxChange, ierr); CHKERRQ(ierr)
               write (IOBuffer, 211) displacementMaxChange
               call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, ierr); CHKERRQ(ierr)

               !!! Conditions to exit the loop in alpha
               if (displacementMaxChange <= MEF90DefMechGlobalOptions%damageATol) then
                  exit
               end if

               if (mod(AltMinIter, 10) == 0) then
                  call MEF90DefMechViewEXO(MEF90DefMechCtx, step, ierr)
               end if
            end do AltMin
            if (AltMinIter == MEF90DefMechGlobalOptions%maxit) then
               write (IOBuffer, 412) MEF90DefMechGlobalOptions%maxit
               call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, ierr); CHKERRQ(ierr)
            end if

            !!! Compute energies
            elasticEnergySet = 0.0_kr
            forceWorkSet = 0.0_kr
            surfaceEnergySet = 0.0_kr
            cohesiveEnergySet = 0.0_kr
            plasticDissipationSet = 0.0_kr

            call MEF90DefMechElasticEnergy(MEF90DefMechCtx%displacement, MEF90DefMechCtx, elasticEnergySet, ierr); CHKERRQ(ierr)
            call MEF90DefMechWork(MEF90DefMechCtx%displacement, MEF90DefMechCtx, forceWorkSet, ierr); CHKERRQ(ierr)
            call MEF90DefMechSurfaceEnergy(MEF90DefMechCtx%damage, MEF90DefMechCtx, surfaceEnergySet, ierr); CHKERRQ(ierr)
            call MEF90DefMechCohesiveEnergy(MEF90DefMechCtx%displacement, MEF90DefMechCtx, cohesiveEnergySet, ierr); CHKERRQ(ierr)
            call MEF90DefMechPlasticDissipation(MEF90DefMechCtx%displacement, MEF90DefMechCtx, plasticStrainOld, plasticDissipationSet, ierr); CHKERRQ(ierr)

            plasticDissipation(step) = sum(plasticDissipationSet)
            elasticEnergy(step) = sum(elasticEnergySet)
            forceWork(step) = sum(forceWorkSet)
            surfaceEnergy(step) = sum(surfaceEnergySet)
            cohesiveEnergy(step) = sum(cohesiveEnergySet)
            totalMechanicalEnergy(step) = elasticEnergy(step) - forceWork(step) + cohesiveEnergy(step) + surfaceEnergy(step) + plasticDissipation(step)
            !!!
            !!! Print and save energies
            !!!
            call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect, 'Cell Sets', CellSetGlobalIS, ierr); CHKERRQ(ierr)
            call MEF90ISAllGatherMerge(MEF90Ctx%Comm, CellSetGlobalIS, ierr); CHKERRQ(ierr)
            call ISGetIndices(CellSetGlobalIS, setID, ierr); CHKERRQ(ierr)
            call PetscPrintf(MEF90Ctx%Comm, "\nMechanical energies: \n", ierr); CHKERRQ(ierr)
            do set = 1, size(setID)
               plasticDissipationSet(set) = plasticDissipationSet(set)
               write (IOBuffer, 201) setID(set), elasticEnergySet(set), forceWorkSet(set), cohesiveEnergySet(set), surfaceEnergySet(set), elasticEnergySet(set) - forceWorkSet(set) + cohesiveEnergySet(set) + surfaceEnergySet(set) + plasticDissipationSet(set), plasticDissipationSet(set)
               call PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr); CHKERRQ(ierr)

               write (IOBuffer, 500) step, time(step), elasticEnergySet(set), forceWorkSet(set), cohesiveEnergySet(set), surfaceEnergySet(set), elasticEnergySet(set) - forceWorkSet(set) + cohesiveEnergySet(set) + surfaceEnergySet(set) + plasticDissipationSet(set), plasticDissipationSet(set)
               call PetscViewerASCIIPrintf(MEF90DefMechCtx%setEnergyViewer(set), IOBuffer, ierr); CHKERRQ(ierr)
               call PetscViewerFlush(MEF90DefMechCtx%setEnergyViewer(set), ierr); CHKERRQ(ierr)
            end do
            call ISRestoreIndices(CellSetGlobalIS, setID, ierr); CHKERRQ(ierr)
            call ISDestroy(CellSetGlobalIS, ierr); CHKERRQ(ierr)
            write (IOBuffer, 202) elasticEnergy(step), forceWork(step), cohesiveEnergy(step), surfaceEnergy(step), totalMechanicalEnergy(step), plasticDissipation(step)
            call PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr); CHKERRQ(ierr)
            write (IOBuffer, 500) step, time(step), elasticEnergy(step), cohesiveEnergy(step), forceWork(step), surfaceEnergy(step), totalMechanicalEnergy(step), plasticDissipation(step)
            call PetscViewerASCIIPrintf(MEF90DefMechCtx%globalEnergyViewer, IOBuffer, ierr); CHKERRQ(ierr)
            call PetscViewerFlush(MEF90DefMechCtx%globalEnergyViewer, ierr); CHKERRQ(ierr)

         case (MEF90DefMech_timeSteppingTypeNULL)
            continue
         case default
            write (IOBuffer, *) "Implemented DefMech time stepping type: ", MEF90DefMechGlobalOptions%timeSteppingType, "\n"
            call PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr); CHKERRQ(ierr)
            stop
         end select
         !!!
         !!! Save results and boundary Values
         !!!
         if (MEF90DefMechGlobalOptions%stressOffset > 0) then
            call MEF90DefMechStress(MEF90DefMechCtx%displacement, MEF90DefMechCtx, MEF90DefMechCtx%stress, ierr)
         end if

         !!! Update plasticstrainold & cumulatedDissipatedPlasticEnergy
         call VecCopy(MEF90DefMechCtx%plasticStrain, plasticStrainOld, ierr); CHKERRQ(ierr)
         call VecCopy(MEF90DefMechCtx%cumulatedPlasticDissipation, cumulatedDissipatedPlasticEnergyOld, ierr); CHKERRQ(ierr)

         call MEF90DefMechViewEXO(MEF90DefMechCtx, step, ierr)

         !!!
         !!! Save performance log file
         !!!
         call PetscViewerASCIIOpen(MEF90Ctx%comm, trim(MEF90FilePrefix(MEF90Ctx%resultFile))//'.log', logViewer, ierr); CHKERRQ(ierr)
         call PetscLogView(logViewer, ierr); CHKERRQ(ierr)
         call PetscViewerFlush(logViewer, ierr); CHKERRQ(ierr)

         if (step == MEF90GlobalOptions%timeNumStep) then
            exit
         else
            step = step + 1
         end if
      end do MainloopQS
   end if
   !!! Clean up and exit nicely
   select case (MEF90DefMechGlobalOptions%timeSteppingType)
   case (MEF90DefMech_timeSteppingTypeQuasiStatic)
      call SNESDestroy(snesDisp, ierr); CHKERRQ(ierr)
      call VecDestroy(residualDisp, ierr); CHKERRQ(ierr)
      call VecDestroy(displacementOld, ierr); CHKERRQ(ierr)
   end select

   select case (MEF90HeatXferGlobalOptions%timeSteppingType)
   case (MEF90HeatXfer_timeSteppingTypeSteadyState)
      call SNESDestroy(snesTemp, ierr); CHKERRQ(ierr)
   case (MEF90HeatXfer_timeSteppingTypeTransient)
      call TSDestroy(tsTemp, ierr); CHKERRQ(ierr)
   end select

   call MEF90DefMechDestroyVectors(MEF90DefMechCtx, ierr)
   nullify (MEF90HeatXferCtx%temperature)
   call MEF90HeatXferDestroyVectors(MEF90HeatXferCtx, ierr)
   call VecDestroy(damageOld, ierr); CHKERRQ(ierr)
   call VecDestroy(plasticStrainPrevious, ierr); CHKERRQ(ierr)
   call VecDestroy(cumulatedDissipatedPlasticEnergyOld, ierr); CHKERRQ(ierr)
   call VecDestroy(cumulatedDissipatedPlasticEnergyVariation, ierr); CHKERRQ(ierr)

   call DMDestroy(Mesh, ierr); CHKERRQ(ierr)

   deallocate (elasticEnergySet)
   deallocate (surfaceEnergySet)
   deallocate (forceWorkSet)
   deallocate (cohesiveEnergySet)
   deallocate (thermalEnergySet)
   deallocate (heatFluxWorkSet)
   deallocate (plasticDissipationSet)

   deallocate (elasticEnergy)
   deallocate (forceWork)
   deallocate (cohesiveEnergy)
   deallocate (surfaceEnergy)
   deallocate (totalMechanicalEnergy)
   deallocate (plasticDissipation)

   call MEF90DefMechDestroy(MEF90DefMechCtx, ierr); CHKERRQ(ierr)
   call MEF90HeatXferDestroy(MEF90HeatXferCtx, ierr); CHKERRQ(ierr)
   call MEF90CtxCloseEXO(MEF90Ctx, ierr)

   call PetscLogView(logViewer, ierr); CHKERRQ(ierr)
   call PetscViewerDestroy(logViewer, ierr); CHKERRQ(ierr)
   call MEF90CtxDestroy(MEF90Ctx, ierr)
   call MEF90Finalize(ierr)
   call PetscFinalize(ierr)
100 format("\nHeat transfer: solving steady state step ", I4, ", t=", ES12.5, "\n")
101 format("cell set ", I4, " thermal energy: ", ES12.5, " fluxes work: ", ES12.5, " total: ", ES12.5, "\n")
102 format("======= Total thermal energy: ", ES12.5, " fluxes work: ", ES12.5, " total: ", ES12.5, "\n")
110 format("\nHeat transfer: step ", I4, ", until t=", ES12.5, "\n")
200 format("\nMechanics: step ", I4, ", t=", ES12.5, "\n")
201 format("cell set ", I4, "  elastic energy: ", ES12.5, " work: ", ES12.5, " cohesive: ", ES12.5, " surface: ", ES12.5, " total: ", ES12.5, " plastic dissipation: ", ES12.5, "\n")
202 format("======= Total: elastic energy: ", ES12.5, " work: ", ES12.5, " cohesive: ", ES12.5, " surface: ", ES12.5, " total: ", ES12.5, " plastic dissipation: ", ES12.5, "\n")
208 format("   Alt. Min. step  ", I5, " \n")
209 format(".  Damage min / max", ES12.5, " / ", ES12.5, ", max change ", ES12.5, "\n")
210 format("      Inner loop ", I4": plastic strain max change ", ES12.5)
211 format("      Displacement max change ", ES12.5, "\n")
308 format("   p-alpha loop  ", I5, " ")
400 format(" [ERROR]: ", A, " SNESSolve failed with SNESConvergedReason ", I2, ". \n Check http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESConvergedReason.html for error code meaning.\n")
410 format(" [ERROR]: ", A, " TSSolve failed with TSConvergedReason ", I2, ". \n Check http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESConvergedReason.html for error code meaning.\n")
412 format(" [ERROR]: Alternate minimizations algorithm did NOT converge in ", I5, "iterations.\n")
500 format(I6, 7(ES16.5), "\n")
end program CoupledPlasticityDamage
