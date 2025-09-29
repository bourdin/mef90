program TestHookesLaw_Type
#include "../MEF90/mef90.inc"
#include "petsc/finclude/petsc.h"
   use m_MEF90

   implicit none(type, external)

   character(len = MEF90MXSTRLEN) :: prefix1 = 'cs0001_'
   character(len = MEF90MXSTRLEN) :: prefix2 = 'cs0002_'
   class(MEF90HookesLaw_Type), allocatable :: HookesLaw1, HookesLaw2, HookesLaw3
   PetscErrorCode :: ierr

   PetscCallA(PetscInitialize(ierr))
   PetscCallA(MEF90Initialize(PETSC_COMM_WORLD, ierr))

   call MEF90GetHookesLaw(PETSC_COMM_WORLD, prefix1, 2, HookesLaw1, ierr)
   call HookesLaw1%SetFromOptions(ierr)
   call HookesLaw1%view(PETSC_VIEWER_STDOUT_WORLD, ierr)

   call MEF90GetHookesLaw(PETSC_COMM_WORLD, prefix2, 2, HookesLaw2, ierr)
   call HookesLaw2%SetFromOptions(ierr)
   call HookesLaw2%view(PETSC_VIEWER_STDOUT_WORLD, ierr)

   HookesLaw3 = MEF90HookesLawSum(HookesLaw1, HookesLaw2)
   call HookesLaw3%view(PETSC_VIEWER_STDOUT_WORLD, ierr)

   PetscCallA(MEF90Finalize(ierr))
   PetscCallA(PetscFinalize(ierr))
end program TestHookesLaw_Type
