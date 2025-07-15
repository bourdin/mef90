module localFunctions
#include <petsc/finclude/petsc.h>
   use petsc
   use m_MEF90
   implicit none(type, external)

contains
   subroutine MyVecView(v, ierr)
      type(tVec), intent(IN)               :: v
      PetscErrorCode, intent(INOUT)        :: ierr

      type(tDM)                           :: dm
      PetscInt                            :: p, pStart, pEnd, numDofClosure
      character(len=MEF90MXSTRLEN)        :: IOBuffer
      PetscScalar, dimension(:), pointer    :: vArray
      PetscInt                            :: height = 0

      PetscCall(VecGetDM(v, dm, ierr))

      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Cell closure\n", ierr))
      PetscCall(DMPlexGetHeightStratum(dm, height, pStart, pEnd, ierr))
      write (*, *) 'cells: ', pStart, pEnd
      do p = pStart, pEnd - 1
         PetscCall(MEF90VecGetClosureSize(v, p, numDofClosure, ierr))
         if (numDofClosure > 0) then
            PetscCall(DMPlexVecGetClosure(dm, PETSC_NULL_SECTION, v, p, vArray, ierr))
            write (IOBuffer, *) p, vArray, "\n"
            PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
            PetscCall(DMPlexVecRestoreClosure(dm, PETSC_NULL_SECTION, v, p, vArray, ierr))
         end if
      end do

      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Point Values\n", ierr))
      PetscCall(DMPlexGetChart(dm, pStart, pEnd, ierr))
      do p = pStart, pEnd - 1
         PetscCall(MEF90VecGetClosureSize(v, p, numDofClosure, ierr))
         if (numDofClosure > 0) then
            PetscCall(DMPlexVecGetClosure(dm, PETSC_NULL_SECTION, v, p, vArray, ierr))
            write (IOBuffer, *) p, vArray, "\n"
            PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
            PetscCall(DMPlexVecRestoreClosure(dm, PETSC_NULL_SECTION, v, p, vArray, ierr))
         end if
      end do
   end subroutine MyVecView
end module localFunctions

program TestDMPlexVecGetClosure
#include <petsc/finclude/petsc.h>
   use m_MEF90
   use petsc
   use localFunctions
   implicit none(type, external)

   PetscErrorCode                      :: ierr
   type(MEF90Ctx_Type), target          :: MEF90Ctx
   type(MEF90CtxGlobalOptions_Type)    :: MEF90GlobalOptions_default

   type(tDM)                           :: dm
   PetscBool                           :: interpolate = PETSC_TRUE
   type(tVec)                          :: v
   type(tPetscSection)                 :: section
   PetscInt                            :: pStart, pEnd, p, numDofClosure
   PetscReal, dimension(:), pointer      :: time

   PetscCallA(PetscInitialize(ierr))
   PetscCallA(MEF90Initialize(PETSC_COMM_WORLD, ierr))

   MEF90GlobalOptions_default % verbose = 1_ki
   MEF90GlobalOptions_default % dryrun = PETSC_FALSE
   MEF90GlobalOptions_default % timeMin = 0.0_kr
   MEF90GlobalOptions_default % timeMax = 1.0_kr
   MEF90GlobalOptions_default % timeNumStep = 11_ki
   MEF90GlobalOptions_default % timeInterpolation = MEF90TimeInterpolation_linear
   MEF90GlobalOptions_default % timeSkip = 0_ki
   MEF90GlobalOptions_default % timeNumCycle = 1_ki
   MEF90GlobalOptions_default % elementFamily = MEF90ElementFamilyLagrange
   MEF90GlobalOptions_default % elementOrder = 1_ki

   PetscCallA(MEF90CtxCreate(PETSC_COMM_WORLD, MEF90Ctx, MEF90GlobalOptions_default, ierr))
   PetscCallA(MEF90CtxGetTime(MEF90Ctx, time, ierr))

   PetscCallA(DMPlexCreateFromFile(MEF90Ctx % Comm, MEF90Ctx % geometryfile, PETSC_NULL_CHARACTER, interpolate, dm, ierr))
   PetscCallA(DMPlexDistributeSetDefault(dm, PETSC_FALSE, ierr))
   PetscCallA(DMSetFromOptions(dm, ierr))
   PetscCallA(DMViewFromOptions(dm, PETSC_NULL_OBJECT, "-dm_view", ierr))

   distribute: block
      type(tDM), target                    :: dmDist
      PetscInt                            :: ovlp = 0
      if (MEF90Ctx % NumProcs > 1) then
         PetscCallA(DMPlexDistribute(dm, ovlp, PETSC_NULL_SF, dmDist, ierr))
         PetscCallA(DMDestroy(dm, ierr))
         dm = dmDist
      end if
   end block distribute
   PetscCallA(DMViewFromOptions(dm, PETSC_NULL_OBJECT, "-mef90dm_view", ierr))

   PetscCallA(PetscSectionCreate(PETSC_COMM_WORLD, section, ierr))
   PetscCallA(DMPlexGetChart(dm, pStart, pEnd, ierr))
   PetscCallA(PetscSectionSetChart(section, pStart, pEnd, ierr))

   PetscCallA(PetscSectionSetDof(section, pStart, 1_ki, ierr))

   PetscCallA(PetscSectionSetUp(section, ierr))

   PetscCallA(DMSetLocalSection(dm, section, ierr))
   PetscCallA(PetscObjectViewFromOptions(section, PETSC_NULL_OBJECT, "-mef90section_view", ierr))

   PetscCallA(DMGetLocalVector(dm, v, ierr))

   PetscCallA(VecSet(v, -1.0_kr, ierr))
   PetscCallA(VecViewFromOptions(v, PETSC_NULL_OBJECT, "-mef90vec_view", ierr))

   do p = pStart, pEnd - 1
      PetscCallA(MEF90VecGetClosureSize(v, p, numDofClosure, ierr))
   end do
   PetscCallA(MyVecView(v, ierr))

   PetscCallA(PetscSectionDestroy(section, ierr))
   PetscCallA(DMRestoreLocalVector(dm, v, ierr))
   PetscCallA(DMDestroy(dm, ierr))

   call MEF90CtxDestroy(MEF90Ctx, ierr)
   call MEF90Finalize(ierr)
   call PetscFinalize(ierr)
end program TestDMPlexVecGetClosure

