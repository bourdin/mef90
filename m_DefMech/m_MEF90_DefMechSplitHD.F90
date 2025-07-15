#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
module MEF90_APPEND(m_MEF90_DefMechSplitHD,MEF90_DIM)D
#define MEF90_DEFMECHSPLITHD_CONSTRUCTOR MEF90_APPEND(m_MEF90_DefMechSplitHD_Constructor,MEF90_DIM)D
#include "petsc/finclude/petsc.h"

use MEF90_APPEND(m_MEF90_DefMechSplit_class,MEF90_DIM)D
use m_MEF90_Materials
implicit none(type)
private
public :: MEF90_DEFMECHSPLITHD

type, extends(MEF90_DEFMECHSPLIT)                   :: MEF90_DEFMECHSPLITHD
   PetscReal                                        :: gamma
contains
   procedure, pass(self)                            :: EED => EEDHD
   procedure, pass(self)                            :: DEED => DEEDHD
   procedure, pass(self)                            :: D2EED => D2EEDHD
end type MEF90_DEFMECHSPLITHD

interface MEF90_DEFMECHSPLITHD
   module procedure MEF90_DEFMECHSPLITHD_CONSTRUCTOR
end interface

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_DEFMECHSPLITHD_CONSTRUCTOR"
!!!
!!!
!!!  MEF90_DEFMECHSPLITHD_CONSTRUCTOR: the default constructor for a MEF90_DEFMECHSPLITHD
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
type(MEF90_DEFMECHSPLITHD) function MEF90_DEFMECHSPLITHD_CONSTRUCTOR(gamma)
   PetscReal, intent(IN)                             :: gamma

   MEF90_DEFMECHSPLITHD_CONSTRUCTOR % gamma = gamma
   MEF90_DEFMECHSPLITHD_CONSTRUCTOR % damageOrder = 3
   MEF90_DEFMECHSPLITHD_CONSTRUCTOR % strainOrder = 2
   MEF90_DEFMECHSPLITHD_CONSTRUCTOR % type = 'MEF90DefMech_unilateralContactTypeHD'
end function MEF90_DEFMECHSPLITHD_CONSTRUCTOR

#undef __FUNCT__
#define __FUNCT__ "EEDHD"
!!!
!!!
!!!  EEDHD: Compute the positive and negative part of the elastic energy density associated with a strain tensor Epsilon
!!!              following the expression from [Amor et al., 2009] W^+ = [sigma(Epsilon)]^+ . Epsilon^+ = sigma(Epsilon) . Epsilon^+
!!!              by orthogality of the masonry projection
!!! [Amor et al., 2009] Amor, H., Marigo, J.-J., and Maurini, C. (2009). Regularized formulation of the variational brittle fracture with unilateral contact:
!!!                     Numerical experiments. J. Mech. Phys. Solids, 57(8):1209 – 1229.
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
subroutine EEDHD(self, Strain, HookesLaw, EEDPlus, EEDMinus)
   class(MEF90_DEFMECHSPLITHD), intent(IN)           :: self
   type(MEF90_MATS), intent(IN)                      :: Strain
   type(MEF90_HOOKESLAW), intent(IN)                 :: HookesLaw
   PetscReal, intent(OUT)                           :: EEDPlus, EEDMinus

   PetscErrorCode                                   :: ierr
   character(len=MEF90MXSTRLEN)                     :: IOBuffer

   if (HookesLaw % type /= MEF90HookesLawTypeIsotropic) then
      write (IOBuffer, *) "Hydrostatic-Deviatoric projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
      PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
   end if

   EEDMinus = MEF90_DefMechSplit_SmoothPositiveSquare(-trace(Strain), self % gamma) * (HookesLaw % lambda + 2.0_kr * HookesLaw % mu / MEF90_DIM) * 0.5_kr
   EEDPlus = ((HookesLaw * Strain) .dotP.Strain) * 0.5_kr - EEDMinus
end subroutine

#undef __FUNCT__
#define __FUNCT__ "DEEDHD"
!!!
!!!
!!!  DEEDHD: Compute the derivative of the positive and negative part of the elastic energy density (positive and negative stress)
!!!               evaluated at the strain tensor Strain.
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
subroutine DEEDHD(self, Strain, HookesLaw, DEEDPlus, DEEDMinus)
   class(MEF90_DEFMECHSPLITHD), intent(IN)           :: self
   type(MEF90_MATS), intent(IN)                      :: Strain
   type(MEF90_HOOKESLAW), intent(IN)                 :: HookesLaw
   type(MEF90_MATS), intent(OUT)                     :: DEEDPlus, DEEDMinus

   PetscErrorCode                                   :: ierr
   character(len=MEF90MXSTRLEN)                     :: IOBuffer

   if (HookesLaw % type /= MEF90HookesLawTypeIsotropic) then
      write (IOBuffer, *) "Hydrostatic-Deviatoric projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
      PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
   end if

   DEEDMinus = (-MEF90_DefMechSplit_DSmoothPositiveSquare(-trace(Strain), self % gamma) * (HookesLaw % lambda + 2.0_kr * HookesLaw % mu / MEF90_DIM) * 0.5_kr) * MEF90_MATS_IDENTITY
   DEEDPlus = (HookesLaw * Strain) - DEEDMinus
end subroutine

#undef __FUNCT__
#define __FUNCT__ "D2EEDHD"
!!!
!!!
!!!  D2EEDHD: Compute the derivative of the positive and negative part of the elastic energy density (positive and negative stress)
!!!               evaluated at the strain tensor Strain.
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
subroutine D2EEDHD(self, Strain, HookesLaw, D2EEDPlus, D2EEDMinus)
   class(MEF90_DEFMECHSPLITHD), intent(IN)           :: self
   type(MEF90_MATS), intent(IN)                      :: Strain
   type(MEF90_HOOKESLAW), intent(IN)                 :: HookesLaw
   type(MEF90_HOOKESLAW), intent(OUT)                :: D2EEDPlus, D2EEDMinus

   PetscErrorCode                                   :: ierr
   character(len=MEF90MXSTRLEN)                     :: IOBuffer

   if (HookesLaw % type /= MEF90HookesLawTypeIsotropic) then
      write (IOBuffer, *) "Hydrostatic-Deviatoric projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
      PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
   end if

   D2EEDPlus % fullTensor = 0.0_kr
   D2EEDMinus % fullTensor = 0.0_kr
   D2EEDPlus % type = MEF90HookesLawTypeIsotropic
   D2EEDMinus % type = MEF90HookesLawTypeIsotropic
   D2EEDMinus % lambda = MEF90_DefMechSplit_D2SmoothPositiveSquare(-trace(Strain), self % gamma) * (HookesLaw % lambda + 2.0_kr * HookesLaw % mu / MEF90_DIM) * 0.5_kr
   D2EEDMinus % mu = 0.0_kr
#if MEF90_DIM == 2
   D2EEDMinus % isPlaneStress = HookesLaw % isPlaneStress
   D2EEDMinus % YoungsModulus = 2.0_kr * D2EEDMinus % mu * (1.0_kr + D2EEDMinus % PoissonRatio)
   D2EEDMinus % BulkModulus = D2EEDMinus % lambda + D2EEDMinus % mu
   if (HookesLaw % isPlaneStress) then
      D2EEDMinus % PoissonRatio = D2EEDMinus % lambda / (D2EEDMinus % lambda + D2EEDMinus % mu) * 0.5_kr
   else
      D2EEDMinus % PoissonRatio = D2EEDMinus % lambda / (D2EEDMinus % lambda + 2.0_kr * D2EEDMinus % mu) * 0.5_kr
   end if
#else
   D2EEDMinus % PoissonRatio = D2EEDMinus % lambda / (D2EEDMinus % lambda + D2EEDMinus % mu) * 0.5_kr
   D2EEDMinus % YoungsModulus = D2EEDMinus % mu * (3.0_kr * D2EEDMinus % lambda + 2.0_kr * D2EEDMinus % mu) / (D2EEDMinus % lambda + D2EEDMinus % mu)
   D2EEDMinus % BulkModulus = D2EEDMinus % lambda + D2EEDMinus % mu * 2.0_kr / 3.0_kr
#endif
   D2EEDPlus = HookesLaw - D2EEDMinus
end subroutine
end module MEF90_APPEND(m_MEF90_DefMechSplitHD,MEF90_DIM)D
