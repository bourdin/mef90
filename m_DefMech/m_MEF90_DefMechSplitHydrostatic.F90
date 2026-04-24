#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
module MEF90_APPEND(m_MEF90_DefMechSplitHydrostatic,MEF90_DIM)D
#define MEF90_DEFMECHSPLITHYDROSTATIC_CONSTRUCTOR MEF90_APPEND(m_MEF90_DefMechSplitHydrostatic_Constructor,MEF90_DIM)D
#include "petsc/finclude/petsc.h"

use MEF90_APPEND(m_MEF90_DefMechSplit_class,MEF90_DIM)D
use m_MEF90_Materials
implicit none(type)
private
public :: MEF90_DEFMECHSPLITHYDROSTATIC

type, extends(MEF90_DEFMECHSPLIT)                   :: MEF90_DEFMECHSPLITHYDROSTATIC
   PetscReal                                        :: gamma
contains
   procedure, pass(self)                            :: EED => EEDHydrostatic
   procedure, pass(self)                            :: DEED => DEEDHydrostatic
   procedure, pass(self)                            :: D2EED => D2EEDHydrostatic
end type MEF90_DEFMECHSPLITHYDROSTATIC

interface MEF90_DEFMECHSPLITHYDROSTATIC
   module procedure MEF90_DEFMECHSPLITHYDROSTATIC_CONSTRUCTOR
end interface

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_DEFMECHSPLITHYDROSTATIC_CONSTRUCTOR"
!!!
!!!
!!!  MEF90_DEFMECHSPLITHYDROSTATIC_CONSTRUCTOR: the default constructor for a MEF90_DEFMECHSPLITHydrostatic
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
type(MEF90_DEFMECHSPLITHYDROSTATIC) function MEF90_DEFMECHSPLITHYDROSTATIC_CONSTRUCTOR(gamma)
   PetscReal, intent(IN)                             :: gamma

   MEF90_DEFMECHSPLITHYDROSTATIC_CONSTRUCTOR%gamma = gamma
   ! MEF90_DEFMECHSPLITHYDROSTATIC_CONSTRUCTOR%damageOrder = 3
   MEF90_DEFMECHSPLITHYDROSTATIC_CONSTRUCTOR%quadratureOrder = 2
   MEF90_DEFMECHSPLITHYDROSTATIC_CONSTRUCTOR%type = 'MEF90DefMech_unilateralContactTypeHydrostatic'
end function MEF90_DEFMECHSPLITHYDROSTATIC_CONSTRUCTOR

#undef __FUNCT__
#define __FUNCT__ "EEDHydrostatic"
!!!
!!!
!!!  EEDHydrostatic: Compute the positive and negative part of the elastic energy density associated with a strain tensor Epsilon
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
subroutine EEDHydrostatic(self, Strain, HookesLaw, EEDPlus, EEDMinus)
   class(MEF90_DEFMECHSPLITHYDROSTATIC), intent(IN) :: self
   type(MEF90_MATS), intent(IN)                     :: Strain
   type(MEF90_HOOKESLAW), intent(IN)                :: HookesLaw
   PetscReal, intent(OUT)                           :: EEDPlus, EEDMinus

   PetscErrorCode                                   :: ierr
   character(len=MEF90MXSTRLEN)                     :: IOBuffer

   if (HookesLaw%type /= MEF90HookesLawTypeIsotropic) then
      write (IOBuffer, *) "Hydrostatic projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
      PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
   end if

   EEDPlus = MEF90_DefMechSplit_SmoothPositiveSquare(trace(Strain), self%gamma) * (HookesLaw%lambda + 2.0_kr * HookesLaw%mu / MEF90_DIM) * 0.5_kr
   ! Ae^s.e^s /2
   EEDMinus = ((HookesLaw * Strain) .dotP.Strain) * 0.5_kr - EEDPlus
end subroutine

#undef __FUNCT__
#define __FUNCT__ "DEEDHydrostatic"
!!!
!!!
!!!  DEEDHydrostatic: Compute the derivative of the positive and negative part of the elastic energy density (positive and negative stress)
!!!               evaluated at the strain tensor Strain.
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
subroutine DEEDHydrostatic(self, Strain, HookesLaw, DEEDPlus, DEEDMinus)
   class(MEF90_DEFMECHSPLITHYDROSTATIC), intent(IN) :: self
   type(MEF90_MATS), intent(IN)                     :: Strain
   type(MEF90_HOOKESLAW), intent(IN)                :: HookesLaw
   type(MEF90_MATS), intent(OUT)                    :: DEEDPlus, DEEDMinus

   PetscErrorCode                                   :: ierr
   character(len=MEF90MXSTRLEN)                     :: IOBuffer

   if (HookesLaw%type /= MEF90HookesLawTypeIsotropic) then
      write (IOBuffer, *) "Hydrostatic projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
      PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
   end if

   DEEDPlus = MEF90_DefMechSplit_DSmoothPositiveSquare(trace(Strain), self%gamma) * ((HookesLaw%lambda + 2.0_kr * HookesLaw%mu / MEF90_DIM) * 0.5_kr) * MEF90_MATS_IDENTITY
   DEEDMinus = (HookesLaw * Strain) - DEEDPlus
end subroutine

#undef __FUNCT__
#define __FUNCT__ "D2EEDHydrostatic"
!!!
!!!
!!!  D2EEDHydrostatic: Compute the derivative of the positive and negative part of the elastic energy density (positive and negative stress)
!!!               evaluated at the strain tensor Strain.
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
subroutine D2EEDHydrostatic(self, Strain, HookesLaw, D2EEDPlus, D2EEDMinus)
   class(MEF90_DEFMECHSPLITHYDROSTATIC), intent(IN) :: self
   type(MEF90_MATS), intent(IN)                     :: Strain
   type(MEF90_HOOKESLAW), intent(IN)                :: HookesLaw
   type(MEF90_HOOKESLAW), intent(OUT)               :: D2EEDPlus, D2EEDMinus

   PetscErrorCode                                   :: ierr
   character(len=MEF90MXSTRLEN)                     :: IOBuffer

   if (HookesLaw%type /= MEF90HookesLawTypeIsotropic) then
      write (IOBuffer, *) "Hydrostatic projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
      PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
   end if

   D2EEDMinus%fullTensor = 0.0_kr
   D2EEDPlus%fullTensor = 0.0_kr
   D2EEDPlus%type = MEF90HookesLawTypeIsotropic
   D2EEDPlus%type = MEF90HookesLawTypeIsotropic
   D2EEDPlus%lambda = MEF90_DefMechSplit_D2SmoothPositiveSquare(trace(Strain), self%gamma) * (HookesLaw%lambda + 2.0_kr * HookesLaw%mu / MEF90_DIM) * 0.5_kr
   D2EEDPlus%mu = 0.0_kr
#if MEF90_DIM == 2
   D2EEDPlus%isPlaneStress = HookesLaw%isPlaneStress
   D2EEDPlus%YoungsModulus = 2.0_kr * D2EEDMinus%mu * (1.0_kr + D2EEDMinus%PoissonRatio)
   D2EEDPlus%BulkModulus = D2EEDMinus%lambda + D2EEDMinus%mu
   if (HookesLaw%isPlaneStress) then
      D2EEDPlus%PoissonRatio = D2EEDPlus%lambda / (D2EEDPlus%lambda + D2EEDPlus%mu) * 0.5_kr
   else
      D2EEDPlus%PoissonRatio = D2EEDPlus%lambda / (D2EEDPlus%lambda + 2.0_kr * D2EEDPlus%mu) * 0.5_kr
   end if
#else
   D2EEDPlus%PoissonRatio = D2EEDPlus%lambda / (D2EEDPlus%lambda + D2EEDPlus%mu) * 0.5_kr
   D2EEDPlus%YoungsModulus = D2EEDPlus%mu * (3.0_kr * D2EEDPlus%lambda + 2.0_kr * D2EEDPlus%mu) / (D2EEDPlus%lambda + D2EEDPlus%mu)
   D2EEDPlus%BulkModulus = D2EEDPlus%lambda + D2EEDPlus%mu * 2.0_kr / 3.0_kr
#endif
   D2EEDMinus = HookesLaw - D2EEDPlus
end subroutine
end module MEF90_APPEND(m_MEF90_DefMechSplitHydrostatic,MEF90_DIM)D
