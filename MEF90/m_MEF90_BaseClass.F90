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
      procedure(MEF90ObjectSetFromOptionsInterface), pass(self), deferred :: setFromOptions
      procedure(MEF90ObjectViewInterface), pass(self), deferred           :: view_internal
      procedure, pass(self)                                               :: view => MEF90Object_view
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
contains
   subroutine MEF90Object_view(self, viewer, ierr)
      class(MEF90Object), intent(in) :: self
      type(tPetscViewer), intent(in) :: viewer
      PetscErrorCode, intent(inout)  :: ierr

      PetscBool                      :: helpPrinted, flg
      ! PetscViewer                    :: MEF90ObjectViewer

      PetscCall(PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-" // trim(self%name) // "_viewed", helpPrinted, PETSC_NULL_BOOL, ierr))
      if (.not. helpPrinted) then   
         ! if (.not. present(viewer)) then
         !    PetscCall(PetscViewerASCIIGetStdout(self%comm, MEF90ObjectViewer, ierr))
         ! else 
         !    MEF90ObjectViewer = viewer
         ! end if
         call self%view_internal(viewer, ierr)
         PetscCall(PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-" //trim(self%name) // "_viewed", "true", ierr))
         ! This is just to silence a warning about the option being set but not used. The option is used to prevent printing the same object multiple times when -verbose is set to a high value.
         PetscCall(PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-" // trim(self%name) // "_viewed", flg, PETSC_NULL_BOOL, ierr))
      end if
   end subroutine MEF90Object_view  
end module m_MEF90_baseClass
