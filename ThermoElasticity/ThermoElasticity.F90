#include "../MEF90/mef90.inc"
program ThermoElasticity
#include "petsc/finclude/petsc.h"
   use m_MEF90
   use m_MEF90_DefMech_class
   use m_MEF90_DefMech
   use m_MEF90_HeatXfer
   use m_MEF90_HeatXfer_class
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

   type(tDM), target                                   :: dm, temperatureDM, displacementDM
   type(tIS)                                          :: setIS
   PetscInt, dimension(:), pointer                      :: setID
   PetscInt                                           :: numCellSet
   PetscInt                                           :: numFaceSet
   PetscInt                                           :: set
   PetscReal, dimension(:), pointer                     :: time, energy, bodyForceWork, boundaryForceWork

   type(tSNES)                                        :: displacementSNES
   type(eSNESConvergedReason)                         :: displacementSNESConvergedReason
   type(tVec)                                         :: displacement, displacementResidual

   type(tSNES)                                        :: temperatureSNES
   type(eSNESConvergedReason)                         :: temperatureSNESConvergedReason
   type(tTS)                                          :: temperatureTS
   type(tTSAdapt)                                     :: temperatureTSAdapt
   type(tVec)                                         :: temperature, temperatureResidual

   PetscReal                                          :: temperatureInitialTimeStep, temperatureInitialTime
   !PetscInt                                           :: tsTemperatureMaxIter

   PetscBool                                          :: flg, EXONeedsFormatting = PETSC_FALSE
   character(len=MEF90MXSTRLEN)                       :: IOBuffer
   type(tPetscViewer)                                 :: logViewer

   PetscInt                                           :: step

   !!! Initialize MEF90
   PetscCallA(PetscInitialize(ierr))
   PetscCallA(MEF90Initialize(PETSC_COMM_WORLD, ierr))

   !!! Get all MEF90-wide options
   PetscCallA(MEF90CtxCreate(PETSC_COMM_WORLD, MEF90Ctx, "", ierr))
   PetscCallA(MEF90Ctx%setFromOptions(ierr))
   MEF90GlobalOptions => MEF90Ctx%globalOptions

   if (MEF90GlobalOptions%verbose > 1) then
      PetscCallA(PetscPrintf(PETSC_COMM_WORLD, "Reading geometry\n", ierr))
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
      if (MEF90GlobalOptions%verbose > 1) then
         PetscCallA(PetscPrintf(PETSC_COMM_WORLD, "Opening result file\n", ierr))
      end if
      PetscCallA(MEF90CtxOpenEXO(MEF90Ctx, MEF90Ctx%resultViewer, FILE_MODE_APPEND, ierr))
   else
      ! we need to create the output file
      if (MEF90GlobalOptions%verbose > 1) then
         PetscCallA(PetscPrintf(PETSC_COMM_WORLD, "Creating result file\n", ierr))
      end if
      PetscCallA(PetscViewerDestroy(MEF90Ctx%resultViewer, ierr))
      PetscCallA(MEF90CtxOpenEXO(MEF90Ctx, MEF90Ctx%resultViewer, FILE_MODE_WRITE, ierr))
      PetscCallA(MEF90EXODMView(dm, MEF90Ctx%resultViewer, MEF90GlobalOptions%elementOrder, ierr))
      EXONeedsFormatting = PETSC_TRUE
   end if

   distribute: block
      type(tDM), target                    :: dmDist
      PetscInt                            :: ovlp = 0_ki
      type(tPetscSF)                      :: naturalPointSF

      if (MEF90Ctx%NumProcs > 1) then
         if (MEF90GlobalOptions%verbose > 1) then
            PetscCallA(PetscPrintf(PETSC_COMM_WORLD, "Distributing mesh\n", ierr))
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
   PetscCallA(MEF90HeatXferCreate(MEF90HeatXferCtx, dm, MEF90Ctx, "", ierr))
   PetscCallA(MEF90HeatXferCtx%setFromOptions(ierr))
   MEF90HeatXferGlobalOptions = MEF90HeatXferCtx%globalOptions

   !!! Create DefMechCtx, get all defMech options
   PetscCallA(MEF90DefMechCreate(MEF90DefMechCtx, dm, MEF90Ctx, "", ierr))
   PetscCallA(MEF90DefMechCtx%setFromOptions(ierr))
   PetscCallA(MEF90DefMechGlobalOptionsSetFromOptions(MEF90DefMechCtx%comm, trim(MEF90DefMechCtx%prefix), MEF90DefMechGlobalOptions, ierr))
   PetscCallA(VecDestroy(MEF90DefMechCtx%temperatureLocal, ierr))
   deallocate (MEF90DefMechCtx%temperatureLocal)
   MEF90DefMechCtx%temperatureLocal => MEF90HeatXferCtx%temperatureLocal

   PetscCallA(DMGetDimension(dm, MEF90DefMechCtx%dim, ierr))
   PetscCallA(MEF90CtxGetTime(MEF90Ctx, time, ierr))
   if (EXONeedsFormatting) then
      if (MEF90GlobalOptions%verbose > 1) then
         PetscCallA(PetscPrintf(PETSC_COMM_WORLD, "Formatting result file\n", ierr))
      end if
      PetscCallA(MEF90DefMechFormatEXO(MEF90DefMechCtx, time, ierr))
      if (MEF90GlobalOptions%verbose > 1) then
         PetscCallA(PetscPrintf(PETSC_COMM_WORLD, "Done Formatting result file\n", ierr))
      end if
   end if

   !!! We no longer need the DM. We have the megaDM in MEF90HeatXferCtx and MEF90DefMechCtx
   PetscCallA(DMDestroy(dm, ierr))


   !!! Create GLOBAL vectors for the unknowns (temperature,displacements), residuals, etc
   PetscCallA(VecGetDM(MEF90HeatXferCtx%temperatureLocal, temperatureDM, ierr))
   !!! This only borrows a reference so we do not need to delete it
   PetscCallA(DMCreateGlobalVector(temperatureDM, temperature, ierr))
   PetscCallA(PetscObjectSetName(temperature, "Temperature", ierr))
   PetscCallA(VecDuplicate(temperature, temperatureResidual, ierr))

   PetscCallA(VecGetDM(MEF90DefMechCtx%displacementLocal, displacementDM, ierr))
   !!! This only borrows a reference so we do not need to delete it
   PetscCallA(DMCreateGlobalVector(displacementDM, displacement, ierr))
   PetscCallA(PetscObjectSetName(displacement, "displacement", ierr))
   PetscCallA(VecDuplicate(displacement, displacementResidual, ierr))

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
   case (MEF90DefMech_TimeSteppingTypeNULL)
      continue
   end select
   !!!
   !!! Allocate array of works and energies
   !!!
   PetscCallA(MEF90DMGetNumSets(MEF90HeatXferCtx%megaDM, MEF90CellSetLabelName, numCellSet, ierr))
   PetscCallA(MEF90DMGetNumSets(MEF90HeatXferCtx%megaDM, MEF90FaceSetLabelName, numFaceSet, ierr))
   allocate (energy(numCellSet))
   allocate (bodyForceWork(numCellSet))
   allocate (boundaryForceWork(numFaceSet))

   !!!
   !!! Actual computations / time stepping
   !!!
   if ((MEF90DefMechGlobalOptions%timeSteppingType /= MEF90DefMech_TimeSteppingTypeNULL) .or. (MEF90HeatXferGlobalOptions%timeSteppingType /= MEF90DefMech_TimeSteppingTypeNULL)) then
      do step = 1, MEF90GlobalOptions%timeNumStep
         write (IOBuffer, 100) step, time(step)
         PetscCallA(PetscPrintf(MEF90Ctx%comm, IOBuffer, ierr))

         !!! Solve for temperature
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

         if (MEF90HeatXferGlobalOptions%timeSteppingType /= MEF90DefMech_TimeSteppingTypeNULL) then
            !!! Compute energies
            PetscCallA(MEF90HeatXFerEnergy(MEF90HeatXferCtx, energy, bodyForceWork, boundaryForceWork, ierr))
            PetscCallA(DMGetLabelIdIS(temperatureDM, MEF90CellSetLabelName, setIS, ierr))
            PetscCallA(MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm, setIS, ierr))
            PetscCallA(ISGetIndices(setIS, setID, ierr))
            do set = 1, size(setID)
               write (IOBuffer, 101) setID(set), energy(set), bodyForceWork(set), energy(set) - bodyForceWork(set)
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

            write (IOBuffer, 102) sum(energy), sum(bodyForceWork) + sum(boundaryForceWork), sum(energy) - sum(bodyForceWork) - sum(boundaryForceWork)
            PetscCallA(PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr))
            !!! Save results
            PetscCallA(MEF90HeatXferViewEXO(MEF90HeatXferCtx, step, ierr))
         end if

         !!! Solve for displacement
         select case (MEF90DefMechGlobalOptions%timeSteppingType)
         case (MEF90DefMech_timeSteppingTypeQuasiStatic)
            PetscCallA(MEF90DefMechSetTransients(MEF90DefMechCtx, step, time(step), ierr))
            PetscCallA(DMLocalToGlobal(displacementDM, MEF90DefMechCtx%displacementLocal, INSERT_VALUES, displacement, ierr))
            !!! Solve SNES
            PetscCallA(SNESSolve(displacementSNES, PETSC_NULL_VEC, displacement, ierr))
            PetscCallA(SNESGetConvergedReason(displacementSNES, displacementSNESConvergedReason, ierr))
            if (displacementSNESConvergedReason%v < 0) then
               write (IOBuffer, 400) "displacement", displacementSNESConvergedReason
               PetscCallA(PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr))
            end if
            PetscCallA(DMGlobalToLocal(displacementDM, displacement, INSERT_VALUES, MEF90DefMechCtx%displacementLocal, ierr))

            !!! Compute energies
            energy = 0.0_kr
            bodyForceWork = 0.0_kr
            boundaryForceWork = 0.0_kr
            PetscCallA(MEF90DefMechWork(MEF90DefMechCtx, bodyForceWork, boundaryForceWork, ierr))
            PetscCallA(MEF90DefMechElasticEnergy(MEF90DefMechCtx, energy, ierr))
            PetscCallA(MEF90DefMechStress(MEF90DefMechCtx, MEF90DefMechCtx%stress, ierr))
            PetscCallA(DMGetLabelIdIS(displacementDM, MEF90CellSetLabelName, setIS, ierr))
            PetscCallA(MEF90ISAllGatherMerge(PETSC_COMM_WORLD, setIS, ierr))
            PetscCallA(ISGetIndices(setIS, setID, ierr))
            do set = 1, size(setID)
               write (IOBuffer, 201) setID(set), energy(set), bodyForceWork(set), energy(set) - bodyForceWork(set)
               PetscCallA(PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr))
            end do
            PetscCallA(ISRestoreIndices(setIS, setID, ierr))
            PetscCallA(ISDestroy(setIS, ierr))

            PetscCallA(DMGetLabelIdIS(displacementDM, MEF90FaceSetLabelName, setIS, ierr))
            PetscCallA(MEF90ISAllGatherMerge(PETSC_COMM_WORLD, setIS, ierr))
            PetscCallA(ISGetIndices(setIS, setID, ierr))
            do set = 1, size(setID)
               write (IOBuffer, 203) setID(set), boundaryForceWork(set)
               PetscCallA(PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr))
            end do
            PetscCallA(ISRestoreIndices(setIS, setID, ierr))
            PetscCallA(ISDestroy(setIS, ierr))
            write (IOBuffer, 202) sum(energy), sum(bodyForceWork) + sum(boundaryForceWork), sum(energy) - sum(bodyForceWork) - sum(boundaryForceWork)
            PetscCallA(PetscPrintf(MEF90Ctx%Comm, IOBuffer, ierr))

            !!! Save results and boundary Values
            PetscCallA(MEF90DefMechViewEXO(MEF90DefMechCtx, step, ierr))
         end select
      end do
   end if
100 format("\nSolving steady state step ", I4, ", t=", ES12.5, "\n")
200 format("\nSolving transient step ", I4, ", t=", ES12.5, "\n")
101 format("cell set ", I4, " thermal energy: ", ES12.5, " flux: ", ES12.5, " total: ", ES12.5, "\n")
102 format("======= Total thermal energy: ", ES12.5, " flux: ", ES12.5, " total: ", ES12.5, "\n")
103 format("face set ", I4, "                              flux: ", ES12.5, "\n")
201 format("cell set ", I4, " elastic energy: ", ES12.5, " work: ", ES12.5, " total: ", ES12.5, "\n")
203 format("face set ", I4, "                              work: ", ES12.5, " \n")
202 format("======= Total elastic energy: ", ES12.5, " work: ", ES12.5, " total: ", ES12.5, "\n")
400 format(" [ERROR]: ", A, " SNESSolve failed with SNESConvergedReason ", I2, ". \n Check https://petsc.org/release/docs/manualpages/SNES/SNESConvergedReason/ for error code meaning.\n")

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
   end select
   PetscCallA(VecDestroy(displacementResidual, ierr))
   PetscCallA(VecDestroy(displacement, ierr))

   deallocate (time)
   deallocate (energy)
   deallocate (bodyForceWork)
   deallocate (boundaryForceWork)
   PetscCallA(MEF90DefMechDestroy(MEF90DefMechCtx, ierr))
   nullify (MEF90HeatXferCtx%temperatureLocal)
   PetscCallA(MEF90HeatXferDestroy(MEF90HeatXferCtx, ierr))

   PetscCallA(PetscViewerASCIIOpen(MEF90Ctx%comm, trim(MEF90FilePrefix(MEF90Ctx%resultFile))//'.log', logViewer, ierr))
   PetscCallA(PetscLogView(logViewer, ierr))
   PetscCallA(PetscViewerDestroy(logViewer, ierr))
   PetscCallA(PetscViewerDestroy(MEF90Ctx%resultViewer, ierr))
   PetscCallA(MEF90CtxDestroy(MEF90Ctx, ierr))
   PetscCallA(MEF90Finalize(ierr))
   PetscCallA(PetscFinalize(ierr))
end program ThermoElasticity
