module m_MEF90_baseClass
#include "petsc/finclude/petsc.h"
   use petscsys
   use m_MEF90_Parameters
   use, intrinsic :: iso_c_binding
   implicit none(type, external)

   private
   public :: MEF90Object

!!!
!!!
!!!  MEF90_baseClass: The base class for all MEF90 classes
!!!
!!!  (c) 2025 Blaise Bourdin bourdin@mcmaster.ca
!!!

   type, abstract :: MEF90Object
      MPIU_Comm :: comm
      character(len=MEF90MXSTRLEN, kind=c_char) :: prefix
   contains
      procedure(MEF90ObjectSetFromOptionsInterface), pass(self), deferred :: setFromOptions
      procedure(MEF90ObjectViewInterface), pass(self), deferred           :: view
   end type MEF90Object

   abstract interface
      subroutine MEF90ObjectSetFromOptionsInterface(self,ierr)
         use petscsys
         import :: MEF90Object
         class(MEF90Object), intent(inout) :: self
         PetscErrorCode, intent(inout)     :: ierr
      end subroutine MEF90ObjectSetFromOptionsInterface

      subroutine MEF90ObjectViewInterface(self, viewer,ierr)
         use petscsys
         import :: MEF90Object
         class(MEF90Object), intent(in) :: self
         type(tPetscViewer), intent(in) :: viewer
         PetscErrorCode, intent(inout)  :: ierr
      end subroutine MEF90ObjectViewInterface
   end interface
end module m_MEF90_baseClass
