module m_MEF90_HookesLaw_class
#include "petsc/finclude/petsc.h"
   use m_MEF90_Parameters
   use m_MEF90_BaseClass, only: MEF90Object
   use petscsys

   implicit none(type, external)
   private
   public :: MEF90HookesLaw
   ! public :: MEF90HookesLaw_setFromOptions

!!!
!!!
!!!  MEF90HookesLaw: The abstract class used to define a stress-strain relation
!!!
!!!  (c) 2025 Blaise Bourdin bourdin@mcmaster.ca
!!!

   type, abstract, extends(MEF90Object) :: MEF90HookesLaw
   contains
      procedure(HookesLawMultInterface), deferred     :: mult
      procedure(HookesLawMultMultInterface), deferred :: multmult
   end type MEF90HookesLaw

   abstract interface
      subroutine HookesLawMultInterface(A, phi, Aphi, ierr)
         use m_MEF90_LinAlg_class, only: mef90Mat
         use petscsys
         import :: MEF90HookesLaw

         class(MEF90HookesLaw), intent(in)         :: A
         class(mef90Mat), intent(in)               :: phi
         class(mef90Mat), allocatable, intent(out) :: Aphi
         PetscErrorCode, intent(inout)             :: ierr
      end subroutine HookesLawMultInterface

      subroutine HookesLawMultMultInterface(A, phi, psi, Aphipsi, ierr)
         use m_MEF90_LinAlg_class, only: mef90Mat
         use petscsys
         import :: MEF90HookesLaw

         class(MEF90HookesLaw), intent(in)        :: A
         class(mef90Mat), intent(in)              :: phi, psi
         PetscReal, intent(out)                   :: Aphipsi
         PetscErrorCode, intent(inout)            :: ierr
      end subroutine HookesLawMultMultInterface
   end interface

! contains
! #undef __FUNCT__
! #define __FUNCT__ "MEF90HookesLaw_setFromOptions"
! !!!
! !!!
! !!!  MEF90HookesLaw_setFromOptions: initializes a MEF90Hookes from options
! !!!  (c) 2025 Blaise Bourdin bourdin@mcmaster.ca
! !!!

!    subroutine MEF90HookesLaw_setFromOptions(self, ierr)
!       class(MEF90HookesLaw), intent(inout) :: self
!       PetscErrorCode,intent(inout)         :: ierr

!       PetscViewer                          :: stdoutViewer
!       PetscInt                             :: verbose

!       PetscCall(PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-verbose", verbose, PETSC_NULL_BOOL, ierr))
!       if (verbose > 0) then
!          PetscCall(PetscViewerASCIIGetStdout(self%comm, stdoutViewer, ierr))
!          call self%view(stdoutViewer, ierr)
!       end if
!    end subroutine MEF90HookesLaw_setFromOptions
end module m_MEF90_HookesLaw_class
