#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
module MEF90_APPEND(m_MEF90_DefMechSplitDeviatoric,MEF90_DIM)D
#define MEF90_DEFMECHSPLITDEVIATORIC_CONSTRUCTOR MEF90_APPEND(m_MEF90_DefMechSplitDeviatoric_Constructor,MEF90_DIM)D
#include "petsc/finclude/petsc.h"

use m_MEF90_Materials
use MEF90_APPEND(m_MEF90_DefMechSplit_class,MEF90_DIM)D
implicit none(type)
private
public :: MEF90_DEFMECHSPLITDEVIATORIC

type, extends(MEF90_DEFMECHSPLIT)                   :: MEF90_DEFMECHSPLITDEVIATORIC
contains
   procedure, pass(self)                            :: EED => EEDDeviatoric
   procedure, pass(self)                            :: DEED => DEEDDeviatoric
   procedure, pass(self)                            :: D2EED => D2EEDDeviatoric
end type MEF90_DEFMECHSPLITDEVIATORIC

interface MEF90_DEFMECHSPLITDEVIATORIC
   module procedure MEF90_DEFMECHSPLITDEVIATORIC_CONSTRUCTOR
end interface

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_DEFMECHSPLITDEVIATORIC_CONSTRUCTOR"
!!!
!!!
!!!  MEF90_DEFMECHSPLITDEVIATORIC_CONSTRUCTOR: the default constructor for a MEF90_DEFMECHSPLITDeviatoric
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
type(MEF90_DEFMECHSPLITDEVIATORIC) function MEF90_DEFMECHSPLITDEVIATORIC_CONSTRUCTOR()

   MEF90_DEFMECHSPLITDEVIATORIC_CONSTRUCTOR % damageOrder = 0
   MEF90_DEFMECHSPLITDEVIATORIC_CONSTRUCTOR % strainOrder = 2
   MEF90_DEFMECHSPLITDEVIATORIC_CONSTRUCTOR % type = 'MEF90DefMech_unilateralContactTypeDeviatoric'
end function MEF90_DEFMECHSPLITDEVIATORIC_CONSTRUCTOR

#undef __FUNCT__
#define __FUNCT__ "EEDDeviatoric"
!!!
!!!
!!!  EEDDeviatoric: Compute the positive and negative part of the elastic energy density associated with a strain tensor Epsilon
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
subroutine EEDDeviatoric(self, Strain, HookesLaw, EEDPlus, EEDMinus)
   class(MEF90_DEFMECHSPLITDEVIATORIC), intent(IN)   :: self
   type(MEF90_MATS), intent(IN)                      :: Strain
   type(MEF90_HOOKESLAW), intent(IN)                 :: HookesLaw
   PetscReal, intent(OUT)                           :: EEDPlus, EEDMinus

   PetscErrorCode                                   :: ierr
   character(len=MEF90MXSTRLEN)                     :: IOBuffer

   if (HookesLaw % type /= MEF90HookesLawTypeIsotropic) then
      write (IOBuffer, *) "Deviatoric projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
      PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
   end if

   EEDMinus = trace(Strain)**2 * (HookesLaw % lambda + 2.0_kr * HookesLaw % mu / MEF90_DIM) * 0.5_kr ! Ae^s.e^s /2
   EEDPlus = ((HookesLaw * Strain) .dotP.Strain) * 0.5_kr - EEDMinus
end subroutine

#undef __FUNCT__
#define __FUNCT__ "DEEDDeviatoric"
!!!
!!!
!!!  DEEDDeviatoric: Compute the derivative of the positive and negative part of the elastic energy density (positive and negative stress)
!!!               evaluated at the strain tensor Strain.
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
subroutine DEEDDeviatoric(self, Strain, HookesLaw, DEEDPlus, DEEDMinus)
   class(MEF90_DEFMECHSPLITDEVIATORIC), intent(IN)   :: self
   type(MEF90_MATS), intent(IN)                      :: Strain
   type(MEF90_HOOKESLAW), intent(IN)                 :: HookesLaw
   type(MEF90_MATS), intent(OUT)                     :: DEEDPlus, DEEDMinus

   PetscErrorCode                                   :: ierr
   character(len=MEF90MXSTRLEN)                     :: IOBuffer

   if (HookesLaw % type /= MEF90HookesLawTypeIsotropic) then
      write (IOBuffer, *) "Deviatoric projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
      PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
   end if

   DEEDMinus = (trace(Strain) * (HookesLaw % lambda + 2.0_kr * HookesLaw % mu / MEF90_DIM)) * MEF90_MATS_IDENTITY
   DEEDPlus = (HookesLaw * Strain) - DEEDMinus
end subroutine

#undef __FUNCT__
#define __FUNCT__ "D2EEDDeviatoric"
!!!
!!!
!!!  D2EEDDeviatoric: Compute the derivative of the positive and negative part of the elastic energy density (positive and negative stress)
!!!               evaluated at the strain tensor Strain.
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
subroutine D2EEDDeviatoric(self, Strain, HookesLaw, D2EEDPlus, D2EEDMinus)
   class(MEF90_DEFMECHSPLITDEVIATORIC), intent(IN)   :: self
   type(MEF90_MATS), intent(IN)                      :: Strain
   type(MEF90_HOOKESLAW), intent(IN)                 :: HookesLaw
   type(MEF90_HOOKESLAW), intent(OUT)                :: D2EEDPlus, D2EEDMinus

   PetscErrorCode                                   :: ierr
   character(len=MEF90MXSTRLEN)                     :: IOBuffer

   if (HookesLaw % type /= MEF90HookesLawTypeIsotropic) then
      write (IOBuffer, *) "Deviatoric projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
      PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
   end if

   D2EEDPlus % fullTensor = 0.0_kr
   D2EEDMinus % fullTensor = 0.0_kr
   D2EEDPlus % type = MEF90HookesLawTypeIsotropic
   D2EEDMinus % type = MEF90HookesLawTypeIsotropic
   D2EEDMinus % lambda = (HookesLaw % lambda + 2.0_kr * HookesLaw % mu / MEF90_DIM) * 0.5_kr
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
end module MEF90_APPEND(m_MEF90_DefMechSplitDeviatoric,MEF90_DIM)D
