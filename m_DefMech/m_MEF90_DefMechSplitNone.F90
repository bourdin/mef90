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

      self%damageOrder = 0
      self%strainOrder = 2
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
#define __FUNCT__ "EEDNone"
!!!
!!!
!!!  EEDNone: Compute the positive and negative part of the elastic energy density associated with a strain tensor
!!!           without a split, we have EEDPlus  = 1/2 HookesLaw Strain \cdot Strain
!!!                                    EEDMinus = 0
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!      2026 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine EEDNone(self, Strain, HookesLaw, EEDPlus, EEDMinus)
      class(MEF90DefMechSplitNone), intent(IN)      :: self
      class(mef90Mat), intent(IN)                   :: Strain
      class(MEF90HookesLaw), allocatable            :: HookesLaw
      PetscReal, intent(OUT)                        :: EEDPlus, EEDMinus

      PetscErrorCode                                :: ierr

      Call HookesLaw%multmult(Strain, Strain, EEDPlus, ierr)
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
   subroutine DEEDNone(self, Strain, HookesLaw, DEEDPlus, DEEDMinus)
      class(MEF90DefMechSplitNone), intent(IN)      :: self
      class(mef90Mat), intent(IN)                   :: Strain
      class(MEF90HookesLaw), allocatable            :: HookesLaw
      class(mef90Mat), allocatable, intent(OUT)     :: DEEDPlus, DEEDMinus

      PetscErrorCode                                :: ierr

      select type (Strain)
         type is (MatS2D)
            DEEDMinus = MatS2D()
            DEEDPlus = MatS2D()
         type is (MatS3D)
            DEEDMinus = MatS3D()
            DEEDPlus = MatS3D()
      end select
      call HookesLaw%mult(Strain, DEEDPlus, ierr)
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
   subroutine D2EEDNone(self, Strain, HookesLaw, D2EEDPlus, D2EEDMinus)
      class(MEF90DefMechSplitNone), intent(IN)        :: self
      class(mef90Mat), intent(IN)                     :: Strain
      class(MEF90HookesLaw), allocatable, intent(IN)  :: HookesLaw
      class(MEF90HookesLaw), allocatable, intent(OUT) :: D2EEDPlus, D2EEDMinus

      D2EEDPlus = HookesLaw
      D2EEDMinus = MEF90HookesLawZero(comm = HookesLaw%comm,  prefix = HookesLaw%prefix)
   end subroutine D2EEDNone
end module m_MEF90_DefMechSplitNone
