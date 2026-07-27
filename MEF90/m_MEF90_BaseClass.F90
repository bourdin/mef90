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
      character(len=MEF90MXSTRLEN, kind=c_char) :: name = "NULL"
   contains
      procedure, pass(self)                                               :: setFromOptions => MEF90Object_setFromOptions
      procedure(MEF90ObjectViewInterface), pass(self), deferred           :: view_internal
      procedure, pass(self)                                               :: view => MEF90Object_view
   end type MEF90Object

   abstract interface
      subroutine MEF90ObjectViewInterface(self, viewer, ierr)
         use petscsys
         import :: MEF90Object
         class(MEF90Object), intent(in) :: self
         type(tPetscViewer), intent(in) :: viewer
         PetscErrorCode, intent(inout)  :: ierr
      end subroutine MEF90ObjectViewInterface
   end interface
contains

#undef __FUNCT__
#define __FUNCT__ "MEF90Object_view"
!!!
!!!
!!!  MEF90Object_view: Views a MEF90Object and flagged it as viewed so that it is not viewed multiple times when -verbose is set to a high value
!!!  (c) 2026 Blaise Bourdin bourdin@mcmaster.ca
!!!

subroutine MEF90Object_view(self, viewer, ierr)
      class(MEF90Object), intent(in) :: self
      type(tPetscViewer), intent(in) :: viewer
      PetscErrorCode, intent(inout)  :: ierr

      PetscBool                      :: helpPrinted, flg

      PetscCall(PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-" // trim(self%name) // "_viewed", helpPrinted, PETSC_NULL_BOOL, ierr))
      if (.not. helpPrinted) then
         call self%view_internal(viewer, ierr)
         PetscCall(PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-" //trim(self%name) // "_viewed", "true", ierr))
         ! This is just to silence a warning about the option being set but not used. The option is used to prevent printing the same object multiple times when -verbose is set to a high value.
         PetscCall(PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-" // trim(self%name) // "_viewed", flg, PETSC_NULL_BOOL, ierr))
      end if
   end subroutine MEF90Object_view

! contains
#undef __FUNCT__
#define __FUNCT__ "MEF90Object_setFromOptions"
!!!
!!!
!!!  MEF90Object_setFromOptions: generic version of the setFromOptions method for all MEF90Objects, which just prints the object if -verbose is set to a high value
!!!  Can be overridden in classes if they have options to set.
!!!  (c) 2026 Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90Object_setFromOptions(self, ierr)
      class(MEF90Object), intent(inout)    :: self
      PetscErrorCode, intent(inout)        :: ierr

      PetscViewer                          :: stdoutViewer
      PetscInt                             :: verbose = 0

      PetscCall(PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-verbose", verbose, PETSC_NULL_BOOL, ierr))
      if (verbose > 0) then
         PetscCall(PetscViewerASCIIGetStdout(self%comm, stdoutViewer, ierr))
         call self%view(stdoutViewer, ierr)
      end if
   end subroutine MEF90Object_setFromOptions
end module m_MEF90_baseClass
