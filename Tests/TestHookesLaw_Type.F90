program TestHookesLaw_Type
#include "../MEF90/mef90.inc"
#include "petsc/finclude/petsc.h"
   use m_MEF90

   implicit none(type, external)

   character(len = MEF90MXSTRLEN) :: prefix = 'cs0001_'
   class(MEF90HookesLaw_Type), allocatable :: HookesLaw
   PetscErrorCode :: ierr

   PetscCallA(PetscInitialize(ierr))
   PetscCallA(MEF90Initialize(PETSC_COMM_WORLD, ierr))

   call MEF90GetHookesLaw(PETSC_COMM_WORLD, prefix, 2, HookesLaw, ierr)
   call HookesLaw%SetFromOptions(ierr)
   call HookesLaw%view(PETSC_VIEWER_STDOUT_WORLD, ierr)

   PetscCallA(MEF90Finalize(ierr))
   PetscCallA(PetscFinalize(ierr))
end program TestHookesLaw_Type
