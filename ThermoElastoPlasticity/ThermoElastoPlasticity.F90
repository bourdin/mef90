#include "../MEF90/mef90.inc"
program ThermoElastoPlasticity
#include <finclude/petscdef.h>

   use m_MEF90
   use m_MEF90_DefMech_class
   use m_MEF90_DefMech
   use m_MEF90_HeatXfer
   use m_MEF90_HeatXfer_class
   use petsc
   implicit none(type, external)

   PetscErrorCode                                     :: ierr
   type(MEF90Ctx_Type), target                         :: MEF90Ctx
   type(MEF90CtxGlobalOptions_Type), pointer           :: MEF90GlobalOptions
   type(MEF90CtxGlobalOptions_Type), parameter         :: MEF90DefaultGlobalOptions = MEF90CtxGlobalOptions_Type( &
                                                          1, & ! verbose
                                                          PETSC_FALSE, & ! validate
                                                          MEF90TimeInterpolation_linear, & ! timeInterpolation
                                                          0.0_kr, & ! timeMin
                                                          1.0_kr, & ! timeMax
                                                          11, & ! timeNumStep
                                                          0, & ! timeSkip
                                                          1.0_kr, & ! frequency
                                                          MEF90FileFormat_EXOSingle)       ! fileFormat

   !!! Defect mechanics contexts
   type(MEF90DefMech_Type), target                    :: MEF90DefMechCtx
   type(MEF90DefMechGlobalOptions_Type)               :: MEF90DefMechGlobalOptions
   type(MEF90DefMechGlobalOptions_Type), parameter    :: MEF90DefMechDefaultGlobalOptions2D = MEF90DefMechGlobalOptions_Type( &
                                                          MEF90DefMech_TimeSteppingTypeQuasiStatic, & ! solverType
                                                          MEF90DefMech_SolverTypeAltMin, & ! timeSteppingType
                                                          PETSC_TRUE, & ! disp_addNullSpace
                                                          3, & ! DisplacementOffset
                                                          2, & ! DamageOffset
                                                          3, & ! boundaryDisplacementOffset
                                                          0, & ! boundaryDamageOffset
                                                          1, & ! temperatureOffset
                                                          4, & ! ForceOffset
                                                          3, & ! pressureForceOffset
                                                          0, & ! CrackPressureOffset
                                                          0, & ! plasticStrainOffset
                                                          6, & ! StressOffset
                                                          MEF90Scaling_Linear, & ! boundaryDisplacementScaling
                                                          MEF90Scaling_CST, & ! boundaryDamageScaling
                                                          MEF90Scaling_Linear, & ! ForceScaling
                                                          MEF90Scaling_Linear, & ! pressureForceScaling
                                                          MEF90Scaling_Linear, & ! CrackPressureScaling
                                                          1e-3, & ! damage_atol
                                                          1000, & ! maxit
                                                          10, & ! PCLag
                                                          1.0_kr, & ! SOROmega
                                                          0., & ! irrevThres
                                                          MEF90DefMech_BTTypeNULL, & ! BTType
                                                          -1, & ! BTInt
                                                          -1, & ! BTScope
                                                          1.0e-2, & ! BTTol
                                                          1.0e-4, & ! plasticStrainAtol
                                                          1, & ! cumulatedDissipatedPlasticEnergyOffset
                                                          1.0e-3, & ! InjectedVolumeAtol
                                                          0.0_kr, & ! dampingCoefficientDisplacement
                                                          0.0_kr)                    ! dampingCoefficientDamage

   type(MEF90DefMechGlobalOptions_Type), parameter    :: MEF90DefMechDefaultGlobalOptions3D = MEF90DefMechGlobalOptions_Type( &
                                                          MEF90DefMech_TimeSteppingTypeQuasiStatic, & ! solverType
                                                          MEF90DefMech_SolverTypeAltMin, & ! timeSteppingType
                                                          PETSC_TRUE, & ! disp_addNullSpace
                                                          3, & ! DisplacementOffset
                                                          2, & ! DamageOffset
                                                          3, & ! boundaryDisplacementOffset
                                                          0, & ! boundaryDamageOffset
                                                          1, & ! temperatureOffset
                                                          4, & ! ForceOffset
                                                          3, & ! pressureForceOffset
                                                          0, & ! CrackPressureOffset
                                                          0, & ! plasticStrainOffset
                                                          6, & ! StressOffset
                                                          MEF90Scaling_Linear, & ! boundaryDisplacementScaling
                                                          MEF90Scaling_CST, & ! boundaryDamageScaling
                                                          MEF90Scaling_Linear, & ! ForceScaling
                                                          MEF90Scaling_Linear, & ! pressureForceScaling
                                                          MEF90Scaling_Linear, & ! CrackPressureScaling
                                                          1e-3, & ! damage_atol
                                                          1000, & ! maxit
                                                          10, & ! PCLag
                                                          1.0_kr, & ! SOROmega
                                                          0., & ! irrevThres
                                                          MEF90DefMech_BTTypeNULL, & ! BTType
                                                          -1, & ! BTInt
                                                          -1, & ! BTScope
                                                          1.0e-2, & ! BTTol
                                                          1.0e-4, & ! plasticStrainAtol
                                                          1, & ! cumulatedDissipatedPlasticEnergyOffset
                                                          1.0e-3, & ! InjectedVolumeAtol
                                                          0.0_kr, & ! dampingCoefficientDisplacement
                                                          0.0_kr)                    ! dampingCoefficientDamage

   type(MEF90DefMechCellSetOptions_Type), parameter   :: MEF90DefMechDefaultCellSetOptions = MEF90DefMechCellSetOptions_Type( &
                                                          -1, & ! elemTypeShortIDDispl will be overriden
                                                          -1, & ! elemTypeShortIDDamage will be overriden
                                                          [0.0_kr, 0.0_kr, 0.0_kr], & ! force
                                                          0.0_kr, & ! pressureForce
                                                          0.0_kr, & ! CrackPressure
                                                          MEF90DefMech_damageTypeAT1Elastic, & ! damageType
                                                          MEF90DefMech_plasticityTypeNone, & ! plasticityType
                                                          MEF90DefMech_unilateralContactTypeNone, & ! unilateralContactType
                                                          [PETSC_FALSE, PETSC_FALSE, PETSC_FALSE], & ! Has Displacement BC
                                                          [0.0_kr, 0.0_kr, 0.0_kr], & ! boundary Displacement
                                                          PETSC_FALSE, & ! Has Damage BC
                                                          PETSC_FALSE, & ! IsCrackPressureActivated
                                                          PETSC_FALSE, & ! IsWorkControlledActivated
                                                          0._kr)                                      ! Boundary Damage
   type(MEF90DefMechVertexSetOptions_Type), parameter :: MEF90DefMechDefaultVertexSetOptions = MEF90DefMechVertexSetOptions_Type( &
                                                          [PETSC_FALSE, PETSC_FALSE, PETSC_FALSE], & ! Has Displacement BC
                                                          [0.0_kr, 0.0_kr, 0.0_kr], & ! boundary Displacement
                                                          PETSC_FALSE, & ! Has Damage BC
                                                          0.0_kr)                                    ! boundary Damage

   !!! HeatXfer contexts
   type(MEF90HeatXfer_Type), target                   :: MEF90HeatXferCtx
   type(MEF90HeatXferGlobalOptions_Type)              :: MEF90HeatXferGlobalOptions
   type(MEF90HeatXferGlobalOptions_Type), parameter   :: MEF90HeatXferDefaultGlobalOptions = MEF90HeatXferGlobalOptions_Type( &
                                                          MEF90HeatXFer_timeSteppingTypeSteadyState, & ! mode
                                                          PETSC_FALSE, & ! addNullSpace
                                                          1, & ! tempOffset
                                                          0., & ! initialTemperature
                                                          MEF90Scaling_Linear, & ! boundaryTempScaling
                                                          0, & ! boundaryTempOffset
                                                          MEF90Scaling_Linear, & ! externalTempScaling
                                                          2, & ! externalTempOffset
                                                          MEF90Scaling_Linear, & ! fluxScaling
                                                          1)                     ! fluxOffset
   type(MEF90HeatXferCellSetOptions_Type), parameter  :: MEF90HeatXferDefaultCellSetOptions = MEF90HeatXferCellSetOptions_Type( &
                                                          -1, & ! elemTypeShortID will be overriden
                                                          0.0_kr, & ! flux
                                                          0.0_kr, & ! surfaceThermalConductivity
                                                          0.0_kr, & ! externalTemp
                                                          PETSC_FALSE, & ! Has BC
                                                          0.0_kr)          ! boundaryTemp

   type(MEF90HeatXferVertexSetOptions_Type), parameter :: MEF90HeatXferDefaultVertexSetOptions = MEF90HeatXferVertexSetOptions_Type( &
                                                          PETSC_FALSE, & ! Has BC
                                                          0.0_kr)          ! boundaryTemp

   type(DM), target                                    :: Mesh
   type(IS)                                           :: setIS, CellSetGlobalIS
   PetscInt, dimension(:), pointer                      :: setID
   PetscInt                                           :: numCellSet
   PetscInt                                           :: numset, set
   PetscReal, dimension(:), pointer                     :: time, energy, work, plasticDissipation, plasticDissipationvariation
   type(SNES)                                         :: snesDisp
   SNESConvergedReason                                :: snesDispConvergedReason
   type(Vec)                                          :: residualDisp

   !! Plastic
   type(Vec)                                          :: plasticStrainOld
   type(SNES)                                         :: snesTemp
   type(TS)                                           :: tsTemp
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

   !!! Alternate minimization
   PetscInt                                           :: AltMinIter
   PetscReal                                          :: PlasticStrainMaxChange
   type(Vec)                                          :: plasticstrainerror
   type(Vec)                                          :: plasticStrainPrevious

   !!! cumulatedDissipatedPlasticEnergy
   type(Vec)                                          :: cumulatedDissipatedPlasticEnergyOld
   type(Vec)                                          :: cumulatedDissipatedPlasticEnergyVariation

   !!! Initialize MEF90
   call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
   call MEF90Initialize(PETSC_COMM_WORLD, ierr)

   !!! Get all MEF90-wide options
   call MEF90CtxCreate(PETSC_COMM_WORLD, MEF90Ctx, MEF90DefaultGlobalOptions, ierr); CHKERRQ(ierr)
   call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag, MEF90GlobalOptions, ierr); CHKERRQ(ierr)

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
   call MEF90HeatXferCtx%setFromOptions(ierr); CHKERRQ(ierr)
   MEF90HeatXferGlobalOptions = MEF90HeatXferCtx%globalOptions

   !!! Get material properties bags
   if (dim == 2) then
      call MEF90MatPropBagSetFromOptions(MEF90DefMechCtx%MaterialPropertiesBag, MEF90DefMechCtx%DMVect, MEF90Mathium2D, MEF90Ctx, ierr)
   else
      call MEF90MatPropBagSetFromOptions(MEF90DefMechCtx%MaterialPropertiesBag, MEF90DefMechCtx%DMVect, MEF90Mathium3D, MEF90Ctx, ierr)
   end if
   MEF90HeatXferCtx%MaterialPropertiesBag => MEF90DefMechCtx%MaterialPropertiesBag

   !!! Create time array from global options
   call MEF90CtxGetTime(MEF90Ctx, time, ierr)

   !!! Set the data layout
   call MEF90DefMechCtxSetSections(MEF90DefMechCtx, ierr)
   call MEF90HeatXferCtxSetSections(MEF90HeatXferCtx, ierr)

   !!! Create vectors
   call MEF90DefMechCreateVectors(MEF90DefMechCtx, ierr)
   call MEF90HeatXferCreateVectors(MEF90HeatXferCtx, ierr)

   !!! Link the temperature field from the DefMech context to that of the HeatXfer context
   deallocate (MEF90DefMechCtx%temperature)
   !MEF90DefMechCtx%temperature => MEF90HeatXferCtx%temperature

   !!!
   !!! Create sections, vectors, and solvers for DefMech Context
   !!!
   call VecDuplicate(MEF90DefMechCtx%displacement, residualDisp, ierr); CHKERRQ(ierr)
   call PetscObjectSetName(residualDisp, "residualDisp", ierr); CHKERRQ(ierr)
   call MEF90DefMechCreateSNESDisplacement(MEF90DefMechCtx, snesDisp, residualDisp, ierr)
   call VecDuplicate(MEF90HeatXferCtx%temperature, residualTemp, ierr); CHKERRQ(ierr)
   call PetscObjectSetName(residualTemp, "residualTemp", ierr); CHKERRQ(ierr)
      !!!cumulatedDissipatedPlasticEnergy Vectors
   call VecDuplicate(MEF90DefMechCtx%cumulatedPlasticDissipation, cumulatedDissipatedPlasticEnergyOld, ierr); CHKERRQ(ierr)
   call VecDuplicate(MEF90DefMechCtx%cumulatedPlasticDissipation, cumulatedDissipatedPlasticEnergyVariation, ierr); CHKERRQ(ierr)
   call VecCopy(MEF90DefMechCtx%cumulatedPlasticDissipation, cumulatedDissipatedPlasticEnergyOld, ierr); CHKERRQ(ierr)
   call VecDuplicate(MEF90DefMechCtx%plasticStrain, plasticStrainOld, ierr); CHKERRQ(ierr)
   call VecDuplicate(MEF90DefMechCtx%plasticStrain, plasticStrainPrevious, ierr); CHKERRQ(ierr)

   if (MEF90HeatXferGlobalOptions%timeSteppingType /= MEF90HeatXfer_timeSteppingTypeNULL) then
      call MEF90HeatXferCtxSetSections(MEF90HeatXferCtx, ierr)
      call MEF90HeatXferCreateVectors(MEF90HeatXferCtx, ierr)
      call VecDuplicate(MEF90HeatXferCtx%temperature, residualTemp, ierr); CHKERRQ(ierr)
      call PetscObjectSetName(residualTemp, "residualTemp", ierr); CHKERRQ(ierr)
      select case (MEF90HeatXferGlobalOptions%timeSteppingType)
      case (MEF90HeatXFer_timeSteppingTypeSteadyState)
         call MEF90HeatXferCreateSNES(MEF90HeatXferCtx, snesTemp, residualTemp, ierr)
      case (MEF90HeatXFer_timeSteppingTypeTransient)
         call MEF90HeatXferCreateTS(MEF90HeatXferCtx, tsTemp, residualTemp, ierr)
         tsTempInitialStep = (time(size(time)) - time(1)) / (size(time) + 0.0_kr) / 10.0_kr
         tsTempInitialTime = time(1)
         call TSSetInitialTimeStep(tsTemp, tsTempInitialTime, tsTempInitialStep, ierr); CHKERRQ(ierr)
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
   allocate (energy(numCellSet))
   energy = 0.0_kr

   allocate (plasticDissipationvariation(numCellSet))
   plasticDissipationvariation = 0.0_kr

   allocate (plasticDissipation(numCellSet))
   plasticDissipation = 0.0_kr

   allocate (work(numCellSet))
   work = 0.0_kr

   !!!
   !!! Try to figure out if the file was formatted
   !!!
   if (MEF90Ctx%rank == 0) then
      call EXGVP(MEF90Ctx%fileExoUnit, "N", numfield, ierr)
   end if
   call MPI_Bcast(numfield, 1, MPIU_INTEGER, 0, MEF90Ctx%comm, ierr)
   if (numfield == 0) then
      call MEF90DefMechFormatEXO(MEF90DefMechCtx, time, ierr)
      !!! Will have to figure out this one
   end if

   !!!
   !!! Actual computations / time stepping
   !!!
   if (MEF90DefMechGlobalOptions%timeSteppingType == MEF90DefMech_timeSteppingTypeQuasiStatic) then
      do step = 1, MEF90GlobalOptions%timeNumStep
         write (IOBuffer, 100) step, time(step)
         call PetscPrintf(MEF90Ctx%comm, IOBuffer, ierr); CHKERRQ(ierr)

         !!! Solve for temperature
         select case (MEF90HeatXferGlobalOptions%timeSteppingType)
         case (MEF90HeatXFer_timeSteppingTypeSteadyState)
            !!! Update fields
            call MEF90HeatXferSetTransients(MEF90HeatXferCtx, step, time(step), ierr)
            !!! Solve SNES
            call MEF90HeatXferUpdateboundaryTemperature(MEF90HeatXferCtx%temperature, MEF90HeatXferCtx, ierr); 
            call SNESSolve(snesTemp, PETSC_NULL_OBJECT, MEF90HeatXferCtx%temperature, ierr); CHKERRQ(ierr)

            !!! Compute thermal energy
            call MEF90HeatXFerEnergy(MEF90HeatXferCtx%temperature, time(step), MEF90HeatXferCtx, energy, work, ierr); CHKERRQ(ierr)
            call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM, 'Cell Sets', CellSetGlobalIS, ierr); CHKERRQ(ierr)
            call MEF90ISAllGatherMerge(PETSC_COMM_WORLD, CellSetGlobalIS, ierr); CHKERRQ(ierr)
            call ISGetIndices(CellSetGlobalIS, setID, ierr); CHKERRQ(ierr)
            do set = 1, size(setID)
               write (IOBuffer, 101) setID(set), energy(set), work(set), energy(set) - work(set)
               call PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr); CHKERRQ(ierr)
            end do
            call ISRestoreIndices(CellSetGlobalIS, setID, ierr); CHKERRQ(ierr)
            call ISDestroy(CellSetGlobalIS, ierr); CHKERRQ(ierr)
            write (IOBuffer, 102) sum(energy), sum(work), sum(energy) - sum(work)
            call PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr); CHKERRQ(ierr)
            !!! Save results
            call MEF90HeatXferViewEXO(MEF90HeatXferCtx, step, ierr)
         case (MEF90HeatXFer_timeSteppingTypeTransient)
            if (step > 1) then
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
                  call TSGetTime(tsTemp, t, ierr); CHKERRQ(ierr)
                  time(step) = t
               else
                  write (IOBuffer, *) 'TS exceeded analysis time. Skipping step\n'
                  call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, ierr); CHKERRQ(ierr)
               end if
            end if

            !!! Compute thermal energy
            call MEF90HeatXFerEnergy(MEF90HeatXferCtx%temperature, time(step), MEF90HeatXferCtx, energy, work, ierr); CHKERRQ(ierr)
            call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM, 'Cell Sets', CellSetGlobalIS, ierr); CHKERRQ(ierr)
            call MEF90ISAllGatherMerge(PETSC_COMM_WORLD, CellSetGlobalIS, ierr); CHKERRQ(ierr)
            call ISGetIndices(CellSetGlobalIS, setID, ierr); CHKERRQ(ierr)
            do set = 1, size(setID)
               write (IOBuffer, 101) setID(set), energy(set), work(set), energy(set) - work(set)
               call PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr); CHKERRQ(ierr)
            end do
            call ISRestoreIndices(CellSetGlobalIS, setID, ierr); CHKERRQ(ierr)
            call ISDestroy(CellSetGlobalIS, ierr); CHKERRQ(ierr)
            write (IOBuffer, 102) sum(energy), sum(work), sum(energy) - sum(work)
            call PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr); CHKERRQ(ierr)
            !!! Save results
            call MEF90HeatXferViewEXO(MEF90HeatXferCtx, step, ierr)
         end select

         !!! Solve for displacement
         select case (MEF90DefMechGlobalOptions%timeSteppingType)
         case (MEF90DefMech_timeSteppingTypeQuasiStatic)
            !!! Update fields
            call MEF90DefMechSetTransients(MEF90DefMechCtx, step, time(step), ierr)
            call MEF90DefMechUpdateboundaryDisplacement(MEF90DefMechCtx%displacement, MEF90DefMechCtx, ierr)

            !!! beginning of the alternate minimization

            AltMin: do AltMinIter = 1, MEF90DefMechGlobalOptions%maxit
               call PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr); CHKERRQ(ierr)
               !! Duplicated PlasticStraincumulated
               !Call VecDuplicate(MEF90DefMechCtx%plasticStrain,plasticStrainCumulatedDuringAltMin,ierr);CHKERRQ(ierr)

               !!! Solve SNES
               call SNESSolve(snesDisp, PETSC_NULL_OBJECT, MEF90DefMechCtx%displacement, ierr); CHKERRQ(ierr)
               call SNESGetConvergedReason(snesDisp, snesDispConvergedReason, ierr); CHKERRQ(ierr)
               if (snesDispConvergedReason < 0) then
                  write (IOBuffer, *) "SNESConvergedReason returned ", snesDispConvergedReason, "\n"
                  call PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr); CHKERRQ(ierr)
               end if

               !!! Save the PlasticStrainPrevious
               call VecCopy(MEF90DefMechCtx%plasticStrain, plasticStrainPrevious, ierr); CHKERRQ(ierr)
               !!! Solve PlasticProjection
               call MEF90DefMechPlasticStrainUpdate(MEF90DefMechCtx, MEF90DefMechCtx%PlasticStrain, MEF90DefMechCtx%displacement, PlasticStrainOld, plasticStrainPrevious, cumulatedDissipatedPlasticEnergyVariation, cumulatedDissipatedPlasticEnergyOld, ierr); CHKERRQ(ierr)
               !!! add damage in DefMechPlasticStrainUpdate

               call VecAxPy(plasticStrainPrevious, -1.0_kr, MEF90DefMechCtx%plasticStrain, ierr); CHKERRQ(ierr)
               !!! Calculate the Infinity norm in error on PlasticStrain
               call VecNorm(plasticStrainPrevious, NORM_INFINITY, PlasticStrainMaxChange, ierr); CHKERRQ(ierr)
               call VecWAXPY(MEF90DefMechCtx%cumulatedPlasticDissipation, 1.0_kr, cumulatedDissipatedPlasticEnergyOld, cumulatedDissipatedPlasticEnergyVariation, ierr); CHKERRQ(ierr)
               !!! Sum The PlasticStrainPrevious with PlasticStrain
               !Call VecAxPy(MEF90DefMechCtx%plasticStrain,1.0_Kr,plasticStrainPrevious,ierr);CHKERRQ(ierr)
               !Call VecView(MEF90DefMechCtx%plasticStrain,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)

               write (*, *) 'PlasticStrainMaxChange:          ', PlasticStrainMaxChange

               if (PlasticStrainMaxChange <= MEF90DefMechGlobalOptions%plasticStrainATol) then
                  exit
               end if

               if (mod(AltMinIter, 25) == 0) then
                  call MEF90DefMechViewEXO(MEF90DefMechCtx, step, ierr)
               end if

            end do AltMin

            !Call VecView(plasticStrainOld,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)

            !!! --> end of the alternate minimization

            !!! Compute energies
            energy = 0.0_kr
            plasticDissipationvariation = 0.0_kr
            work = 0.0_kr
            call MEF90DefMechWork(MEF90DefMechCtx%displacement, MEF90DefMechCtx, work, ierr); CHKERRQ(ierr)
            call MEF90DefMechElasticEnergy(MEF90DefMechCtx%displacement, MEF90DefMechCtx, energy, ierr); CHKERRQ(ierr)

            call MEF90DefMechPlasticDissipation(MEF90DefMechCtx%displacement, MEF90DefMechCtx, plasticStrainOld, plasticDissipationvariation, ierr); CHKERRQ(ierr)

            call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect, 'Cell Sets', CellSetGlobalIS, ierr); CHKERRQ(ierr)
            call MEF90ISAllGatherMerge(PETSC_COMM_WORLD, CellSetGlobalIS, ierr); CHKERRQ(ierr)
            call ISGetIndices(CellSetGlobalIS, setID, ierr); CHKERRQ(ierr)
            do set = 1, size(setID)
               plasticDissipation(set) = plasticDissipation(set) + plasticDissipationvariation(set)
               write (IOBuffer, 201) setID(set), energy(set), work(set), energy(set) - work(set), plasticDissipation(set)
               call PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr); CHKERRQ(ierr)
            end do

            call ISRestoreIndices(CellSetGlobalIS, setID, ierr); CHKERRQ(ierr)
            call ISDestroy(CellSetGlobalIS, ierr); CHKERRQ(ierr)
            write (IOBuffer, 202) sum(energy), sum(work), sum(energy) - sum(work)
            call PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr); CHKERRQ(ierr)

            !!! Save results and boundary Values
            if (MEF90DefMechGlobalOptions%stressOffset > 0) then
               call MEF90DefMechStress(MEF90DefMechCtx%displacement, MEF90DefMechCtx, MEF90DefMechCtx%stress, ierr)
            end if

            call MEF90DefMechViewEXO(MEF90DefMechCtx, step, ierr)

            !!! Update plasticstrainold & cumulatedDissipatedPlasticEnergy
            call VecCopy(MEF90DefMechCtx%plasticStrain, plasticStrainOld, ierr); CHKERRQ(ierr)
            call VecCopy(MEF90DefMechCtx%cumulatedPlasticDissipation, cumulatedDissipatedPlasticEnergyOld, ierr); CHKERRQ(ierr)

         end select

         !!! update plasticstrainold
         call VecCopy(MEF90DefMechCtx%plasticStrain, plasticStrainOld, ierr); CHKERRQ(ierr)

      end do !!step
   end if
100 format("Solving steady state step ", I4, ", t=", ES12.5, "\n")
101 format("cell set ", I4, " thermal energy: ", ES12.5, " fluxes work: ", ES12.5, " total: ", ES12.5, "\n")
102 format("======= Total thermal energy: ", ES12.5, " fluxes work: ", ES12.5, " total: ", ES12.5, "\n")
201 format("cell set ", I4, " elastic energy: ", ES12.5, " work: ", ES12.5, " total: ", ES12.5, " plastic dissipation: ", ES12.5, "\n")
202 format("======= Total elastic energy: ", ES12.5, " work: ", ES12.5, " total: ", ES12.5, "\n")

   !!! Clean up and exit nicely
   select case (MEF90DefMechGlobalOptions%timeSteppingType)
   case (MEF90DefMech_timeSteppingTypeQuasiStatic)
      call SNESDestroy(snesDisp, ierr); CHKERRQ(ierr)
      call VecDestroy(residualDisp, ierr); CHKERRQ(ierr)
   end select

   select case (MEF90HeatXferGlobalOptions%timeSteppingType)
   case (MEF90HeatXFer_timeSteppingTypeSteadyState)
      call SNESDestroy(snesTemp, ierr); CHKERRQ(ierr)
   case (MEF90HeatXFer_timeSteppingTypeTransient)
      call TSDestroy(tsTemp, ierr); CHKERRQ(ierr)
   end select

   call MEF90DefMechDestroyVectors(MEF90DefMechCtx, ierr)
   nullify (MEF90HeatXferCtx%temperature)
   call MEF90HeatXferDestroyVectors(MEF90HeatXferCtx, ierr)
   call VecDestroy(residualDisp, ierr); CHKERRQ(ierr)
   call VecDestroy(residualTemp, ierr); CHKERRQ(ierr)
   call VecDestroy(plasticStrainPrevious, ierr); CHKERRQ(ierr)
   call VecDestroy(cumulatedDissipatedPlasticEnergyOld, ierr); CHKERRQ(ierr)
   call VecDestroy(cumulatedDissipatedPlasticEnergyVariation, ierr); CHKERRQ(ierr)

   call DMDestroy(Mesh, ierr); CHKERRQ(ierr)

   deallocate (time)
   deallocate (energy)
   deallocate (plasticDissipation)
   deallocate (work)
   deallocate (plasticDissipationvariation)

   call MEF90DefMechDestroy(MEF90DefMechCtx, ierr); CHKERRQ(ierr)
   call MEF90HeatXferDestroy(MEF90HeatXferCtx, ierr); CHKERRQ(ierr)
   call MEF90CtxCloseEXO(MEF90Ctx, ierr)

   call PetscViewerASCIIOpen(MEF90Ctx%comm, trim(MEF90FilePrefix(MEF90Ctx%resultFile))//'.log', logViewer, ierr); CHKERRQ(ierr)
   call PetscLogView(logViewer, ierr); CHKERRQ(ierr)
   call PetscViewerDestroy(logViewer, ierr); CHKERRQ(ierr)
   call MEF90CtxDestroy(MEF90Ctx, ierr)
   call MEF90Finalize(ierr)
   call PetscFinalize(ierr)
end program ThermoElastoPlasticity
