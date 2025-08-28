module m_MEF90_HookesLaw_class
#include "petsc/finclude/petsc.h"
   use m_MEF90_Parameters
   use m_MEF90_Utils
   ! use m_MEF90_LinAlg
   use m_MEF90_BaseClass
   implicit none(type, external)
   private
   public :: MEF90HookesLaw_Type
   public :: MEF90HookesLaw_setFromOptions

!!!
!!!
!!!  MEF90HookesLaw_Type: The abstract class used to define a stress-strain relation
!!!
!!!  (c) 2025 Blaise Bourdin bourdin@mcmaster.ca
!!!

   type, abstract, extends(MEF90Object_Type) :: MEF90HookesLaw_Type
   end type MEF90HookesLaw_Type

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLaw_setFromOptions"
!!!
!!!
!!!  MEF90HookesLaw_setFromOptions: initializes a MEF90_DefMechAT_Type from options
!!!  (c) 2025 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine MEF90HookesLaw_setFromOptions(self,ierr)
      class(MEF90HookesLaw_Type), intent(inout) :: self
      PetscErrorCode,intent(inout) :: ierr

      PetscBool :: printHelp
      PetscCall(PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-verbose", printHelp, PETSC_NULL_BOOL, ierr))
      if (printHelp) then
         call self%view(PETSC_VIEWER_STDOUT_WORLD,ierr)
      end if
   end subroutine MEF90HookesLaw_setFromOptions
end module m_MEF90_HookesLaw_class
