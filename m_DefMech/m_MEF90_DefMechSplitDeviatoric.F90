#include "../MEF90/mef90.inc"
module m_MEF90_DefMechSplitDeviatoric
#include "petsc/finclude/petsc.h"
use m_MEF90_DefMechSplit_class
use m_MEF90_Materials
use m_MEF90_HookesLaw
implicit none(type)
private
public :: MEF90DefMechSplitDeviatoric

type, extends(MEF90DefMechSplit) :: MEF90DefMechSplitDeviatoric
   PetscReal :: gamma = 1.0e-1
contains
   procedure, pass(self)  :: setFromOptions => setFromOptionsDeviatoric
   procedure, pass(self)  :: view_internal => viewDeviatoric
   procedure, pass(self)  :: setup => setupDeviatoric
   procedure, pass(self)  :: EED => EEDDeviatoric
   procedure, pass(self)  :: DEED => DEEDDeviatoric
   procedure, pass(self)  :: D2EED => D2EEDDeviatoric
end type MEF90DefMechSplitDeviatoric


contains
#undef __FUNCT__
#define __FUNCT__ "setFromOptionsDeviatoric"
!!!
!!!
!!!  setFromOptionsDeviatoric: the default constructor for a MEF90DefMechSplitDeviatoric
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!      2026 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine setFromOptionsDeviatoric(self, ierr)
      class(MEF90DefMechSplitDeviatoric), intent(inout) :: self
      PetscErrorCode, intent(inout)             :: ierr
      PetscInt                                  :: verbose = 0

      ! self%damageOrder = 0
      self%quadratureOrder = 2
      self%type = 'MEF90DefMechSplitDeviatoric'

      PetscCall(PetscOptionsBegin(self%comm, trim(self%prefix) // "split_Deviatoric_", "Options for MEF90DefMechSplitDeviatoric_type", "mef90DefMech", ierr))
         PetscCall(PetscOptionsBool("-hybrid", "Use a hybrid split", "MEF90", PETSC_FALSE, self%isHybrid, PETSC_NULL_BOOL, ierr))
      PetscCall(PetscOptionsEnd(ierr))


      PetscCall(PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-verbose", verbose, PETSC_NULL_BOOL, ierr))
      if (verbose > 0) then
         call self%view(PETSC_VIEWER_STDOUT_WORLD,ierr)
      end if
   end subroutine setFromOptionsDeviatoric

#undef __FUNCT__
#define __FUNCT__ "viewDeviatoric"
!!!
!!!
!!!  viewDeviatoric: view a MEF90DefMechSplitDeviatoric
!!!  (c) 2026 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine viewDeviatoric(self, viewer, ierr)
      class(MEF90DefMechSplitDeviatoric), intent(in)    :: self
      type(tPetscViewer), intent(in)              :: viewer
      PetscErrorCode,intent(inout)                :: ierr

      character(len=MEF90MXSTRLEN, kind=c_char)   :: IOBuffer
      character(len=MEF90MXSTRLEN, kind=c_char)   :: viewerType

      PetscCall(PetscViewerGetType(viewer, viewerType, ierr))
      if (viewerType == 'ascii') then
         write(IOBuffer, "(A,': Options for MEF90DefMechSplit\n')") trim(self%prefix) // "split"
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
         write(IOBuffer, "('         type: Deviatoric\n')")
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
         write(IOBuffer, "('         hybrid: ', L1, '\n')") self%isHybrid
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      end if
   end subroutine viewDeviatoric

#undef __FUNCT__
#define __FUNCT__ "setupDeviatoric"
!!!
!!!
!!!  setupDeviatoric: the setup routine for a MEF90DefMechSplitDeviatoric, which does nothing since there is no split
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!      2026 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine setupDeviatoric(self, Strain, ierr)
      use m_MEF90
      implicit none(type, external)

      class(MEF90DefMechSplitDeviatoric), intent(inout) :: self
      class(mef90Mat), intent(IN)               :: Strain
      PetscErrorCode, intent(inout)             :: ierr

      self%strain = Strain
   end subroutine setupDeviatoric


#undef __FUNCT__
#define __FUNCT__ "EEDDeviatoric"
!!!
!!!
!!!  EEDDeviatoric: Compute the positive and negative part of the elastic energy density associated with a strain tensor
!!!           EEDMinus = SmoothPositiveSquare(-trace(e)) AI.I/2/N^2
!!!           EEDPlus  = Ae.e/2 - EEDMinus
!!!  with N = 3 for 3D and N = 2 for 2D, where e is the strain tensor, A is the Hookes law tensor
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!      2026 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine EEDDeviatoric(self, HookesLaw, phi, EEDPlus, EEDMinus, ierr)
      class(MEF90DefMechSplitDeviatoric), intent(IN) :: self
      class(MEF90HookesLaw), intent(IN)      :: HookesLaw
      class(mef90Mat), intent(IN)            :: phi
      PetscReal, intent(OUT)                 :: EEDPlus, EEDMinus
      PetscErrorCode, intent(inout)          :: ierr

      select type(e => phi)
      type is (MatS2D)
         call HookesLaw%multmult(DeviatoricPart(e), e, EEDPlus, ierr)
         call HookesLaw%multmult(HydrostaticPart(e), e, EEDMinus, ierr)
      type is (MatS3D)
         call HookesLaw%multmult(DeviatoricPart(e), e, EEDPlus, ierr)
         call HookesLaw%multmult(HydrostaticPart(e), e, EEDMinus, ierr)
      end select ! self%strain

      EEDMinus = EEDMinus / 2.0_Kr
      EEDPlus  = EEDPlus  / 2.0_Kr
   end subroutine EEDDeviatoric

#undef __FUNCT__
#define __FUNCT__ "DEEDDeviatoric"
!!!
!!!
!!!  DEEDDeviatoric: Compute the directional derivative of the positive and negative part of the elastic energy density in the direction phi at the point defined by the strain tensor
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!      2026 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine DEEDDeviatoric(self, HookesLaw, phi, DEEDPlus, DEEDMinus, ierr)
      class(MEF90DefMechSplitDeviatoric), intent(IN) :: self
      class(MEF90HookesLaw), intent(IN)      :: HookesLaw
      class(mef90Mat), intent(IN)            :: phi
      PetscReal, intent(OUT)                 :: DEEDPlus, DEEDMinus
      PetscErrorCode, intent(inout)          :: ierr

      select type(e => self%strain)
      type is (MatS2D)
         call HookesLaw%multmult(DeviatoricPart(e), phi, DEEDPlus, ierr)
         call HookesLaw%multmult(HydrostaticPart(e), phi, DEEDMinus, ierr)
         !!! Note that we use the orthogonality property A phi^s.psi^D = 0
      type is (MatS3D)
         call HookesLaw%multmult(DeviatoricPart(e), phi, DEEDPlus, ierr)
         call HookesLaw%multmult(HydrostaticPart(e), phi, DEEDMinus, ierr)
      end select ! self%strain
   end subroutine DEEDDeviatoric

#undef __FUNCT__
#define __FUNCT__ "D2EEDDeviatoric"
!!!
!!!
!!!  D2EEDDeviatoric: Compute the second derivative of the positive and negative part of the elastic energy density (positive and negative stress)
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   subroutine D2EEDDeviatoric(self, HookesLaw, phi, psi, D2EEDPlus, D2EEDMinus, ierr)
      class(MEF90DefMechSplitDeviatoric), intent(IN) :: self
      class(MEF90HookesLaw), intent(IN)      :: HookesLaw
      class(mef90Mat), intent(IN)            :: phi, psi
      PetscReal, intent(OUT)                 :: D2EEDPlus, D2EEDMinus
      PetscErrorCode, intent(inout)          :: ierr

      select type(e => psi)
      type is (MatS2D)
         call HookesLaw%multmult(DeviatoricPart(e), phi, D2EEDPlus, ierr)
         call HookesLaw%multmult(HydrostaticPart(e), phi, D2EEDMinus, ierr)
         !!! Note that we use the orthogonality property A phi^s.psi^D = 0
      type is (MatS3D)
         call HookesLaw%multmult(DeviatoricPart(e), phi, D2EEDPlus, ierr)
         call HookesLaw%multmult(HydrostaticPart(e), phi, D2EEDMinus, ierr)
      end select ! self%strain
   end subroutine D2EEDDeviatoric
end module m_MEF90_DefMechSplitDeviatoric
