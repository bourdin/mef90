module m_MEF90_HookesLawZero
#include "petsc/finclude/petsc.h"
   use m_MEF90_Parameters
   use m_MEF90_Utils
   use m_MEF90_LinAlg
   use m_MEF90_HookesLaw_Class
   use iso_c_binding

   implicit none(type)
   private
   public :: MEF90HookesLawZero

   type, extends(MEF90HookesLaw) :: MEF90HookesLawZero
   contains
         procedure :: setFromOptions => MEF90HookesLawZero_setFromOptions
         procedure :: view => MEF90HookesLawZero_view
         procedure :: mult => MEF90HookesLawZero_mult
         procedure :: multmult => MEF90HookesLawZero_multmult
   end type MEF90HookesLawZero

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLawZero_setFromOptions"
!!!
!!!
!!!  MEF90HookesLawZero_setFromOptions: initializes a MEF90HookesLawZero from options
!!!  (c) 2025 Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90HookesLawZero_setFromOptions(self, ierr)
      class(MEF90HookesLawZero), intent(inout)        :: self
      PetscErrorCode,intent(inout)                    :: ierr

      PetscInt                                        :: printHelp

      ! no options
      PetscCall(PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-verbose", printHelp, PETSC_NULL_BOOL, ierr))
      if (printHelp > 1) then
         call self%view(PETSC_VIEWER_STDOUT_WORLD,ierr)
      end if
   end subroutine MEF90HookesLawZero_setFromOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLawZero_view"
!!!
!!!
!!!  MEF90HookesLawZero_view: the default viewer for a MEF90_DefMechAT_Type
!!!  (c) 2025 Blaise Bourdin bourdin@mcmaster.ca
!!!

subroutine MEF90HookesLawZero_View(self, viewer, ierr)
      class(MEF90HookesLawZero), intent(in)        :: self
      type(tPetscViewer), intent(in)               :: viewer
      PetscErrorCode, intent(inout)                :: ierr

      character(len=MEF90MXSTRLEN, kind=c_char)    :: IOBuffer
      character(len=MEF90MXSTRLEN, kind=c_char)    :: viewerType

      PetscCall(PetscViewerGetType(viewer, viewerType, ierr))
      if (viewerType == 'ascii') then
         write(IOBuffer, "(A,': Options for MEF90HookesLaw\n')") trim(self%prefix)//"HookesLaw"
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         write(IOBuffer, "('         Type: Zero\n')")
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
      end if
   end subroutine MEF90HookesLawZero_View

   subroutine MEF90HookesLawZero_mult(A, phi, Aphi, ierr)
      class(MEF90HookesLawZero), intent(in)        :: A
      class(mef90Mat), intent(in)                  :: phi
      class(mef90Mat), allocatable, intent(out)    :: Aphi
      character(len=MEF90MXSTRLEN, kind=c_char)    :: IOBuffer
      PetscErrorCode, intent(inout)                :: ierr

      select type (phinD => phi)
         type is (MatS2D)
            Aphi = MatS2D()
         type is (MatS3D)
            Aphi = MatS3D()
         class default
            write (IOBuffer, *) "Incompatible arguments in "//__FUNCT__//'\n'
            PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
      end select ! phi
   end subroutine MEF90HookesLawZero_mult

   subroutine MEF90HookesLawZero_multmult(A, phi, psi, Aphipsi, ierr)
      class(MEF90HookesLawZero), intent(in)        :: A
      class(mef90Mat), intent(in)                  :: phi, psi
      PetscReal, intent(out)                       :: Aphipsi
      PetscErrorCode, intent(inout)                :: ierr

      Aphipsi = 0.0_Kr

   end subroutine MEF90HookesLawZero_multmult
end module m_MEF90_HookesLawZero
