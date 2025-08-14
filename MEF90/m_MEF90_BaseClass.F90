module m_MEF90_baseClass
#include "petsc/finclude/petsc.h"
   use petscsys
   use m_MEF90_Parameters
   use, intrinsic :: iso_c_binding
   implicit none(type, external)

   private
   public :: MEF90Object_Type

!!!
!!!
!!!  MEF90_baseClass: The base class for all MEF90 classes
!!!
!!!  (c) 2025 Blaise Bourdin bourdin@mcmaster.ca
!!!

   type, abstract :: MEF90Object_Type
      MPI_Comm :: comm
      character(len=MEF90MXSTRLEN, kind=c_char) :: prefix
   contains
      procedure(MEF90Object_SetFromOptionsInterface), pass(self), deferred :: setFromOptions
      procedure(MEF90Object_ViewInterface), pass(self), deferred :: view
   end type MEF90Object_Type

   abstract interface
      subroutine MEF90Object_setFromOptionsInterface(self,ierr)
         use petscsys
         import :: MEF90Object_Type
         class(MEF90Object_Type), intent(inout) :: self
         PetscErrorCode, intent(inout) :: ierr
      end subroutine MEF90Object_setFromOptionsInterface

      subroutine MEF90Object_ViewInterface(self, viewer,ierr)
         use petscsys
         import :: MEF90Object_Type
         class(MEF90Object_Type), intent(in) :: self
         type(tPetscViewer), intent(in) :: viewer
         PetscErrorCode, intent(inout) :: ierr
      end subroutine MEF90Object_ViewInterface
   end interface
end module m_MEF90_baseClass
