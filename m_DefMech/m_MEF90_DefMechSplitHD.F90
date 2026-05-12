#include "../MEF90/mef90.inc"
module m_MEF90_DefMechSplitHD
#include "petsc/finclude/petsc.h"
use m_MEF90_DefMechSplit_class
use m_MEF90_Materials
use m_MEF90_HookesLaw
implicit none(type)
private
public :: MEF90DefMechSplitHD

type, extends(MEF90DefMechSplit) :: MEF90DefMechSplitHD
   PetscReal :: gamma = 1.0e-1
contains
   procedure, pass(self)  :: setFromOptions => setFromOptionsHD
   procedure, pass(self)  :: view => viewHD
   procedure, pass(self)  :: setup => setupHD
   procedure, pass(self)  :: EED => EEDHD
   procedure, pass(self)  :: DEED => DEEDHD
   procedure, pass(self)  :: D2EED => D2EEDHD
end type MEF90DefMechSplitHD


contains
#undef __FUNCT__
#define __FUNCT__ "setFromOptionsHD"
!!!
!!!
!!!  setFromOptionsHD: the default constructor for a MEF90DefMechSplitHD
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!      2026 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine setFromOptionsHD(self, ierr)
      class(MEF90DefMechSplitHD), intent(inout) :: self
      PetscErrorCode, intent(inout)             :: ierr
      PetscBool                                 :: printHelp

      ! self%damageOrder = 0
      self%quadratureOrder = 2
      self%type = 'MEF90DefMechSplitHD'

      PetscCall(PetscOptionsBegin(self%comm, trim(self%prefix) // "split_HydrostaticDeviatoric_", "Options for MEF90DefMechSplitHD_type", "mef90DefMech", ierr))
         PetscCall(PetscOptionsReal('-gamma', 'gamma parameter', '[]', self%gamma, self%gamma, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsBool("-hybrid", "Use a hybrid split", "MEF90", PETSC_FALSE, self%isHybrid, PETSC_NULL_BOOL, ierr))
      PetscCall(PetscOptionsEnd(ierr))

      PetscCall(PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-verbose", printHelp, PETSC_NULL_BOOL, ierr))
      if (printHelp) then
         call self%view(PETSC_VIEWER_STDOUT_WORLD,ierr)
      end if
   end subroutine setFromOptionsHD

#undef __FUNCT__
#define __FUNCT__ "viewHD"
!!!
!!!
!!!  viewHD: view a MEF90DefMechSplitHD
!!!  (c) 2026 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine viewHD(self, viewer, ierr)
      class(MEF90DefMechSplitHD), intent(in)    :: self
      type(tPetscViewer), intent(in)              :: viewer
      PetscErrorCode,intent(inout)                :: ierr

      character(len=MEF90MXSTRLEN, kind=c_char)   :: IOBuffer
      character(len=MEF90MXSTRLEN, kind=c_char)   :: viewerType

      PetscCall(PetscViewerGetType(viewer, viewerType, ierr))
      if (viewerType == 'ascii') then
         write(IOBuffer, "(A,': Options for MEF90DefMechSplit\n')") trim(self%prefix) // "split"
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         write(IOBuffer, "('         type: HydrostaticDeviatoric\n')")
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         write(IOBuffer, "('         hybrid: ', L1, '\n')") self%isHybrid
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         write(IOBuffer, "('         gamma: ',ES12.5,' []\n')") self%gamma
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
      end if
   end subroutine viewHD

#undef __FUNCT__
#define __FUNCT__ "setupHD"
!!!
!!!
!!!  setupHD: the setup routine for a MEF90DefMechSplitHD, which does nothing since there is no split
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!      2026 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine setupHD(self, Strain, ierr)
      use m_MEF90
      implicit none(type, external)

      class(MEF90DefMechSplitHD), intent(inout) :: self
      class(mef90Mat), intent(IN)               :: Strain
      PetscErrorCode, intent(inout)             :: ierr

      self%strain = Strain
   end subroutine setupHD


#undef __FUNCT__
#define __FUNCT__ "EEDHD"
!!!
!!!
!!!  EEDHD: Compute the positive and negative part of the elastic energy density associated with a strain tensor
!!!           EEDMinus = SmoothPositiveSquare(-trace(e)) AI.I/2/N^2 
!!!           EEDPlus  = Ae.e/2 - EEDMinus
!!!  with N = 3 for 3D and N = 2 for 2D, where e is the strain tensor, A is the Hookes law tensor
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!      2026 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine EEDHD(self, HookesLaw, phi, EEDPlus, EEDMinus, ierr)
      class(MEF90DefMechSplitHD), intent(IN) :: self
      class(MEF90HookesLaw), intent(IN)      :: HookesLaw
      class(mef90Mat), intent(IN)            :: phi
      PetscReal, intent(OUT)                 :: EEDPlus, EEDMinus
      PetscErrorCode, intent(inout)          :: ierr
      
      PetscReal                              :: AIIN2 ! AI.I/N^2
      PetscReal                              :: tre ! tr(e)

      AIIN2 = 0.0_Kr
      tre = 0.0_Kr
      select type(HookesLaw)
      type is (MEF90HookesLawIsotropic2D)
         AIIN2 = HookesLaw%BulkModulus
      type is (MEF90HookesLawIsotropic3D)
         AIIN2 = HookesLaw%BulkModulus
      class default
         select type(e => phi)
         type is (MatS2D)
            call HookesLaw%multmult(MEF90MatS2DIdentity, MEF90MatS2DIdentity, AIIN2, ierr)
            AIIN2 = AIIN2 / 4.0_Kr
         type is (MatS3D)
            call HookesLaw%multmult(MEF90MatS3DIdentity, MEF90MatS3DIdentity, AIIN2, ierr)
            AIIN2 = AIIN2 / 9.0_Kr
         end select ! self%strain
      end select ! HookesLaw

      select type(e => phi)
      type is (MatS2D)
         tre = trace(e)
      type is (MatS3D)
         tre = trace(e)
      end select ! phi

      EEDMinus = AIIN2 * MEF90DefMechSplit_SmoothPositiveSquare(-tre, self%gamma) / 2.0_Kr

      call HookesLaw%multmult(phi, phi, EEDPlus, ierr)
      EEDPlus = EEDPlus / 2.0_Kr - EEDMinus
   end subroutine EEDHD

#undef __FUNCT__
#define __FUNCT__ "DEEDHD"
!!!
!!!
!!!  DEEDHD: Compute the directional derivative of the positive and negative part of the elastic energy density in the direction phi at the point defined by the strain tensor
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!      2026 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine DEEDHD(self, HookesLaw, phi, DEEDPlus, DEEDMinus, ierr)
      class(MEF90DefMechSplitHD), intent(IN) :: self
      class(MEF90HookesLaw), intent(IN)      :: HookesLaw
      class(mef90Mat), intent(IN)            :: phi
      PetscReal, intent(OUT)                 :: DEEDPlus, DEEDMinus
      PetscErrorCode, intent(inout)          :: ierr

      PetscReal                              :: AIIN2 ! AI.I/N^2
      PetscReal                              :: trStrain ! tr(self%Strain)
      PetscReal                              :: trPhi ! tr(phi)

      AIIN2 = 0.0_Kr
      trStrain = 0.0_Kr
      trPhi = 0.0_Kr
      select type(HookesLaw)
      type is (MEF90HookesLawIsotropic2D)
         AIIN2 = HookesLaw%BulkModulus
      type is (MEF90HookesLawIsotropic3D)
         AIIN2 = HookesLaw%BulkModulus
      class default
         select type(e => self%strain)
         type is (MatS2D)
            call HookesLaw%multmult(MEF90MatS2DIdentity, MEF90MatS2DIdentity, AIIN2, ierr)
            AIIN2 = AIIN2 / 4.0_Kr
         type is (MatS3D)
            call HookesLaw%multmult(MEF90MatS3DIdentity, MEF90MatS3DIdentity, AIIN2, ierr)
            AIIN2 = AIIN2 / 9.0_Kr
         end select ! self%strain
      end select ! HookesLaw

      select type(e => self%strain)
      type is (MatS2D)
         trStrain = trace(e)
      type is (MatS3D)
         trStrain = trace(e)
      end select ! self%strain

      select type(p =>phi)
      type is (MatS2D)
         trPhi = trace(p)
      type is (MatS3D)
         trPhi = trace(p)
      end select

      DEEDMinus = -AIIN2 * MEF90DefMechSplit_DSmoothPositiveSquare(-trStrain, self%gamma) * trPhi
      call HookesLaw%multmult(phi, self%strain, DEEDPlus, ierr)
      DEEDPlus = DEEDPlus - DEEDMinus 
   end subroutine DEEDHD

#undef __FUNCT__
#define __FUNCT__ "D2EEDHD"
!!!
!!!
!!!  D2EEDHD: Compute the second derivative of the positive and negative part of the elastic energy density (positive and negative stress)
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   subroutine D2EEDHD(self, HookesLaw, phi, psi, D2EEDPlus, D2EEDMinus, ierr)
      class(MEF90DefMechSplitHD), intent(IN) :: self
      class(MEF90HookesLaw), intent(IN)      :: HookesLaw
      class(mef90Mat), intent(IN)            :: phi, psi
      PetscReal, intent(OUT)                 :: D2EEDPlus, D2EEDMinus
      PetscErrorCode, intent(inout)          :: ierr

      PetscReal                              :: AIIN2 ! AI.I/N^2
      PetscReal                              :: trStrain ! tr(self%Strain)
      PetscReal                              :: trPhi, trPsi ! tr(phi), tr(psi)

      AIIN2 = 0.0_Kr
      trStrain = 0.0_Kr
      trPhi = 0.0_Kr
      trPsi = 0.0_Kr
      select type(HookesLaw)
      type is (MEF90HookesLawIsotropic2D)
         AIIN2 = HookesLaw%BulkModulus
      type is (MEF90HookesLawIsotropic3D)
         AIIN2 = HookesLaw%BulkModulus
      class default
         select type(e => self%strain)
         type is (MatS2D)
            call HookesLaw%multmult(MEF90MatS2DIdentity, MEF90MatS2DIdentity, AIIN2, ierr)
            AIIN2 = AIIN2 / 4.0_Kr
            trStrain = trace(e)
         type is (MatS3D)
            call HookesLaw%multmult(MEF90MatS3DIdentity, MEF90MatS3DIdentity, AIIN2, ierr)
            AIIN2 = AIIN2 / 9.0_Kr
            trStrain = trace(e)
         end select ! self%strain
      end select ! HookesLaw

      select type(e => self%strain)
      type is (MatS2D)
         trStrain = trace(e)
      type is (MatS3D)
         trStrain = trace(e)
      end select ! self%strain

      select type(p =>phi)
      type is (MatS2D)
         trPhi = trace(p)
      type is (MatS3D)
         trPhi = trace(p)
      end select

      select type(p =>psi)
      type is (MatS2D)
         trPsi = trace(p)
      type is (MatS3D)
         trPsi = trace(p)
      end select

      D2EEDMinus = AIIN2 * MEF90DefMechSplit_D2SmoothPositiveSquare(-trStrain, self%gamma) * trPhi * trPsi 
      call HookesLaw%multmult(phi, psi, D2EEDPlus, ierr)
      D2EEDPlus = D2EEDPlus - D2EEDMinus
   end subroutine D2EEDHD
end module m_MEF90_DefMechSplitHD
