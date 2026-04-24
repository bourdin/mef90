#include "../MEF90/mef90.inc"
module m_MEF90_DefMechSplitNone
#include "petsc/finclude/petsc.h"
use m_MEF90_DefMechSplit_class
use m_MEF90_Materials
use m_MEF90_HookesLaw
implicit none(type)
private
public :: MEF90DefMechSplitNone

type, extends(MEF90DefMechSplit) :: MEF90DefMechSplitNone
contains
   procedure              :: setFromOptions => MEF90DefMechSplitNone_setFromOptions
   procedure              :: view => MEF90DefMechSplitNone_view
   procedure, pass(self)  :: setup => setupNONE
   procedure, pass(self)  :: EED => EEDNone
   procedure, pass(self)  :: DEED => DEEDNone
   procedure, pass(self)  :: D2EED => D2EEDNone
end type MEF90DefMechSplitNone


contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSplitNone_setFromOptions"
!!!
!!!
!!!  MEF90DefMechSplitNone_setFromOptions: the default constructor for a MEF90DefMechSplitNone
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!      2026 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine MEF90DefMechSplitNone_setFromOptions(self, ierr)
      class(MEF90DefMechSplitNone), intent(inout) :: self
      PetscErrorCode, intent(inout)               :: ierr
      PetscBool                                   :: printHelp

      ! self%damageOrder = 0
      self%quadratureOrder = 2
      self%type = 'MEF90DefMechSplitNone'

      !!! MEF90DefMechSplitNone has no options
      PetscCall(PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-verbose", printHelp, PETSC_NULL_BOOL, ierr))
      if (printHelp) then
         call self%view(PETSC_VIEWER_STDOUT_WORLD,ierr)
      end if

   end subroutine MEF90DefMechSplitNone_setFromOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSplitNone_view"
!!!
!!!
!!!  MEF90DefMechSplitNone_view: view a MEF90DefMechSplitNone
!!!  (c) 2026 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine MEF90DefMechSplitNone_view(self, viewer, ierr)
      class(MEF90DefMechSplitNone), intent(in)    :: self
      type(tPetscViewer), intent(in)              :: viewer
      PetscErrorCode,intent(inout)                :: ierr

      character(len=MEF90MXSTRLEN, kind=c_char)   :: IOBuffer
      character(len=MEF90MXSTRLEN, kind=c_char)   :: viewerType

      PetscCall(PetscViewerGetType(viewer, viewerType, ierr))
      if (viewerType == 'ascii') then
         write(IOBuffer, "(A,': Options for MEF90DefMechSplit\n')") trim(self%prefix) // "split"
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         write(IOBuffer, "('         type: none\n')")
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         write(IOBuffer, "('         No options\n')")
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
      end if
   end subroutine MEF90DefMechSplitNone_view

#undef __FUNCT__
#define __FUNCT__ "setupNONE"
!!!
!!!
!!!  setupNONE: the setup routine for a MEF90DefMechSplitNone, which does nothing since there is no split
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!      2026 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine setupNONE(self, Strain, ierr)
      use m_MEF90
      implicit none(type, external)

      class(MEF90DefMechSplitNone), intent(inout) :: self
      class(mef90Mat), intent(IN)                 :: Strain
      PetscErrorCode, intent(inout)               :: ierr

      self%strain = Strain
   end subroutine setupNONE


#undef __FUNCT__
#define __FUNCT__ "EEDNone"
!!!
!!!
!!!  EEDNone: Compute the positive and negative part of the elastic energy density associated with a strain tensor
!!!           without a split, we have EEDPlus  = 1/2 HookesLaw Strain \cdot Strain
!!!                                    EEDMinus = 0
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!      2026 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine EEDNone(self, HookesLaw, phi, EEDPlus, EEDMinus, ierr)
      class(MEF90DefMechSplitNone), intent(IN) :: self
      class(MEF90HookesLaw), intent(IN)        :: HookesLaw
      class(mef90Mat), intent(IN)              :: phi
      PetscReal, intent(OUT)                   :: EEDPlus, EEDMinus
      PetscErrorCode, intent(inout)            :: ierr

      Call HookesLaw%multmult(phi, phi, EEDPlus, ierr)
      EEDPlus = EEDPlus * 0.5_kr
      EEDMinus = 0.0_kr
   end subroutine EEDNone

#undef __FUNCT__
#define __FUNCT__ "DEEDNone"
!!!
!!!
!!!  DEEDNone: Compute the derivative of the positive and negative part of the elastic energy density (positive and negative stress)
!!!               evaluated at the strain tensor Strain.
!!!           without a split, we have DEEDPlus  = HookesLaw Strain
!!!                                    DEEDMinus = 0
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!      2026 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine DEEDNone(self, HookesLaw, phi, DEEDPlus, DEEDMinus, ierr)
      class(MEF90DefMechSplitNone), intent(IN) :: self
      class(MEF90HookesLaw), intent(IN)        :: HookesLaw
      class(mef90Mat), intent(IN)              :: phi
      PetscReal, intent(OUT)                   :: DEEDPlus, DEEDMinus
      PetscErrorCode, intent(inout)            :: ierr

      call HookesLaw%multmult(self%strain, phi, DEEDPlus, ierr)
      DEEDminus = 0.0_kr
   end subroutine DEEDNone

#undef __FUNCT__
#define __FUNCT__ "D2EEDNone"
!!!
!!!
!!!  D2EEDNone: Compute the second derivative of the positive and negative part of the elastic energy density (positive and negative stress)
!!!               evaluated at the strain tensor Strain.
!!!               without a split, D2EEDPlus = HookesLaw, D2EEDMinus = 0
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   subroutine D2EEDNone(self, HookesLaw, phi, psi, D2EEDPlus, D2EEDMinus, ierr)
      class(MEF90DefMechSplitNone), intent(IN) :: self
      class(MEF90HookesLaw), intent(IN)        :: HookesLaw
      class(mef90Mat), intent(IN)              :: phi, psi
      PetscReal, intent(OUT)                   :: D2EEDPlus, D2EEDMinus
      PetscErrorCode, intent(inout)            :: ierr

      call HookesLaw%multmult(phi, psi, D2EEDPlus, ierr)
      D2EEDMinus = 0.0_kr
   end subroutine D2EEDNone
end module m_MEF90_DefMechSplitNone
