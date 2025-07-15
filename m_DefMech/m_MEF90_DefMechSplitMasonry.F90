#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
module MEF90_APPEND(m_MEF90_DefMechSplitMasonry,MEF90_DIM)D
#include "petsc/finclude/petsc.h"
use MEF90_APPEND(m_MEF90_DefMechSplit_class,MEF90_DIM)D
use m_MEF90_Materials
#define MEF90_DEFMECHSPLITMASONBY_CONSTRUCTOR MEF90_APPEND(m_MEF90_DefMechSplitMasonry_Constructor,MEF90_DIM)D
implicit none(type)
private
public :: MEF90_DEFMECHSPLITMASONRY

type, extends(MEF90_DEFMECHSPLIT)                   :: MEF90_DEFMECHSPLITMASONRY
contains
   procedure, pass(self)                            :: EED => EEDMasonry
   procedure, pass(self)                            :: DEED => DEEDMasonry
   procedure, pass(self)                            :: D2EED => D2EEDMasonry
end type MEF90_DEFMECHSPLITMASONRY

interface MEF90_DEFMECHSPLITMASONRY
   module procedure MEF90_DEFMECHSPLITMASONRY_CONSTRUCTOR
end interface

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_DEFMECHSPLITMASONRY_CONSTRUCTOR"
!!!
!!!
!!!  MEF90_DEFMECHSPLITMASONRY_CONSTRUCTOR: the default constructor for a MEF90_DEFMECHSPLITMASONRY
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
type(MEF90_DEFMECHSPLITMASONRY) function MEF90_DEFMECHSPLITMASONRY_CONSTRUCTOR()
   MEF90_DEFMECHSPLITMASONRY_CONSTRUCTOR % damageOrder = 0
   MEF90_DEFMECHSPLITMASONRY_CONSTRUCTOR % strainOrder = 2
   MEF90_DEFMECHSPLITMASONRY_CONSTRUCTOR % type = 'MEF90DefMech_unilateralContactTypeMasonry'
end function MEF90_DEFMECHSPLITMASONRY_CONSTRUCTOR

#undef __FUNCT__
#define __FUNCT__ "EEDMasonry"
!!!
!!!
!!!  EEDMasonry: Compute the positive and negative part of the elastic energy density associated with a strain tensor Epsilon
!!!              following the expression from [Li, 2016] W^+ = [sigma(Epsilon)]^+ . Epsilon^+ = sigma(Epsilon) . Epsilon^+
!!!              by orthogality of the masonry projection
!!!
!!!  [Li, 2016] Li, T. (2016). Gradient Damage Modeling of Dynamic Brittle Fracture. PhD thesis, Universite Paris-Saclay – Ecole Polytechnique.
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
subroutine EEDMasonry(self, Strain, HookesLaw, EEDPlus, EEDMinus)
   class(MEF90_DEFMECHSPLITMASONRY), intent(IN)        :: self
   type(MEF90_MATS), intent(IN)                        :: Strain
   type(MEF90_HOOKESLAW), intent(IN)                   :: HookesLaw
   PetscReal, intent(OUT)                             :: EEDPlus, EEDMinus

   type(MEF90_MATS)                                   :: D, DPlus
   type(MEF90_MAT)                                    :: Pinv
   PetscReal                                          :: nu
#if MEF90_DIM == 2
   PetscReal                                          :: alpha
#endif
   PetscErrorCode                                     :: ierr
   character(len=MEF90MXSTRLEN)                       :: IOBuffer

   if (HookesLaw % type /= MEF90HookesLawTypeIsotropic) then
      write (IOBuffer, *) "Masonry projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
      PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
   end if
   call Diagonalize(Strain, Pinv, D)
      !!! D is the strain tensor in the principal basis

   nu = HookesLaw % PoissonRatio
   DPlus = 0.0_kr

#if MEF90_DIM == 2
   if (HookesLaw % isPlaneStress) then
      alpha = nu / (1.0_kr - nu)
   else
      alpha = nu / (1.0_kr - 2.0_kr * nu)
   end if
   if (D % XX >= 0.0_kr) then
      DPlus = D
   else if ((1.0_kr + alpha) * D % YY + alpha * D % XX >= 0.0_kr) then
      Dplus % YY = alpha / (1.0_kr + alpha) * D % XX + D % YY
   end if
#else
   if (D % XX >= 0.0_kr) then
      DPlus = D
   else if (nu * D % XX + D % YY >= 0.0_kr) then
      DPlus % YY = nu * D % XX + D % YY
      DPlus % ZZ = nu * D % XX + D % ZZ
   else if (nu * (D % XX + D % YY) + (1.0_kr - nu) * D % ZZ >= 0.0_kr) then
      DPlus % ZZ = nu / (1.0_kr - nu) * (D % YY + D % XX) + D % ZZ
   end if
#endif
      !!! We compute the Elastic energy density in the principal basis
   EEDPlus = (HookesLaw * D.dotP.DPlus) * 0.5_kr
   EEDMinus = (HookesLaw * D.dotP. (D - DPlus)) * 0.5_kr
end subroutine EEDMasonry

#undef __FUNCT__
#define __FUNCT__ "DEEDMasonry"
!!!
!!!
!!!  DEEDMasonry: Compute the derivative of the positive and negative part of the elastic energy density (positive and negative stress)
!!!               evaluated at the strain tensor Strain.
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
subroutine DEEDMasonry(self, Strain, HookesLaw, DEEDPlus, DEEDMinus)
   class(MEF90_DEFMECHSPLITMASONRY), intent(IN)        :: self
   type(MEF90_MATS), intent(IN)                        :: Strain
   type(MEF90_HOOKESLAW), intent(IN)                   :: HookesLaw
   type(MEF90_MATS), intent(OUT)                       :: DEEDPlus, DEEDMinus

   type(MEF90_MATS)                                   :: D, StrainPlus
   type(MEF90_MAT)                                    :: Pinv
   PetscReal                                          :: E, nu
#if MEF90_DIM == 2
   PetscReal                                          :: alpha
#endif
   PetscErrorCode                                     :: ierr
   character(len=MEF90MXSTRLEN)                       :: IOBuffer

   if (HookesLaw % type /= MEF90HookesLawTypeIsotropic) then
      write (IOBuffer, *) "Masonry projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
      PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
   end if

   call Diagonalize(Strain, Pinv, D)
      !!! D is the strain tensor in the principal basis

   E = HookesLaw % YoungsModulus
   nu = HookesLaw % PoissonRatio
   StrainPlus = 0.0_kr
   DEEDPlus = 0.0_kr
   DEEDMinus = 0.0_kr
#if MEF90_DIM == 2
   if (HookesLaw % isPlaneStress) then
      alpha = nu / (1.0_kr - nu)
   else
      alpha = nu / (1.0_kr - 2.0_kr * nu)
   end if
   if (D % XX >= 0.0_kr) then
      DEEDPlus = HookesLaw * Strain
      DEEDMinus = 0.0_kr
   else if ((1.0_kr + alpha) * D % YY + alpha * D % XX >= 0.0_kr) then
         !!! Compute the projection of the strain in the principal basis
      StrainPlus % XX = 0.0_kr
      StrainPlus % YY = alpha / (1.0_kr + alpha) * D % XX + D % YY
         !!! Compute Sigma(Strain^+) then pull back to the canonical basis
      DEEDPlus = MEF90MatRaRt(HookesLaw * StrainPlus, Pinv)
      DEEDMinus = HookesLaw * Strain - DEEDPlus
   else
      DEEDPlus = 0.0_kr
      DEEDMinus = HookesLaw * Strain
   end if
#else
   if (D % XX >= 0.0_kr) then
      DEEDPlus = HookesLaw * Strain
      DEEDMinus = 0.0_kr
   else if (nu * D % XX + D % YY >= 0.0_kr) then
      StrainPlus % YY = nu * D % XX + D % YY
      StrainPlus % ZZ = nu * D % XX + D % ZZ
         !!! Compute Sigma(Strain^+) then pull back to the canonical basis
      DEEDPlus = MEF90MatRaRt(HookesLaw * StrainPlus, Pinv)
      DEEDMinus = HookesLaw * Strain - DEEDPlus
   else if (nu * (D % XX + D % YY) + (1.0_kr - nu) * D % ZZ >= 0.0_kr) then
      StrainPlus % ZZ = nu / (1.0_kr - nu) * (D % XX + D % YY) + D % ZZ
         !!! Compute Sigma(Strain^+) then pull back to the canonical basis
      DEEDPlus = MEF90MatRaRt(HookesLaw * StrainPlus, Pinv)
      DEEDMinus = HookesLaw * Strain - DEEDPlus
   else
      DEEDPlus = 0.0_kr
      DEEDMinus = HookesLaw * Strain
   end if
#endif
end subroutine DEEDMasonry

#undef __FUNCT__
#define __FUNCT__ "D2EEDMasonry"
!!!
!!!
!!!  D2EEDMasonry: Compute the second derivative of the positive and negative part of the elastic energy density (positive and negative stress)
!!!               evaluated at the strain tensor Strain.
!!!               writing Sigma^{+,-}(E) = (g_1^{+,-}(e1,e2,e3), g_2^{+,-}(21,e2,e3), g_3^{+,-}(e1,e2,e3)) in the basis of principal strains
!!!               we get that D2EED = \partial Sigma^+/partial E i.e.
!!!               D2EED^{+,-}_{iijj} = g^{+,-}_{i,j}
!!!               D2EED^{+,-}_{ijij} = H^{+,-}_{ij}
!!!               with H_{ij}(e) = (g_i(a) - g_j(a)) / (a_i - a_j) if a_i /= a_j
!!!                                 g_{i,i}(a) - g_{i,j}(a)        if a_i  = a_j
!!!                                 0                              if i = j (which is consistant with the fact that g_{i,j} = g_{j,i} here)
!!! [Silhavy, 1997] Silhavy, M. (1997). The Mechanics and Thermodynamics of Continuous Media. Springer Berlin Heidelberg.
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
subroutine D2EEDMasonry(self, Strain, HookesLaw, D2EEDPlus, D2EEDMinus)
   class(MEF90_DEFMECHSPLITMASONRY), intent(IN)        :: self
   type(MEF90_MATS), intent(IN)                        :: Strain
   type(MEF90_HOOKESLAW), intent(IN)                   :: HookesLaw
   type(MEF90_HOOKESLAW), intent(OUT)                  :: D2EEDPlus, D2EEDMinus

   type(MEF90_TENS4OS)                                :: A
   type(MEF90_MATS)                                   :: D
   type(MEF90_MAT)                                    :: Pinv
   PetscReal                                          :: E, nu
#if MEF90_DIM == 2
   PetscReal                                          :: alpha
#endif
   PetscErrorCode                                     :: ierr
   character(len=MEF90MXSTRLEN)                       :: IOBuffer

   if (HookesLaw % type /= MEF90HookesLawTypeIsotropic) then
      write (IOBuffer, *) "Masonry projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
      PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
   end if

   call Diagonalize(Strain, Pinv, D)
      !!! D is the strain tensor in the principal basis

   E = HookesLaw % YoungsModulus
   nu = HookesLaw % PoissonRatio
   A = 0.0_kr
   D2EEDPlus % fullTensor = 0.0_kr
   D2EEDMinus % fullTensor = 0.0_kr
#if MEF90_DIM == 2
   if (HookesLaw % isPlaneStress) then
      alpha = nu / (1.0_kr - nu)
   else
      alpha = nu / (1.0_kr - 2.0_kr * nu)
   end if
   if (D % XX >= 0.0_kr) then
      D2EEDPlus % type = MEF90HookesLawTypeIsotropic
      D2EEDMinus % type = MEF90HookesLawTypeIsotropic
      D2EEDPlus % YoungsModulus = E
      D2EEDPlus % PoissonRatio = nu
      D2EEDPlus % lambda = HookesLaw % lambda
      D2EEDPlus % mu = HookesLaw % mu

      D2EEDMinus % YoungsModulus = 0.0_kr
      D2EEDMinus % PoissonRatio = 0.0_kr
      D2EEDMinus % lambda = 0.0_kr
      D2EEDMinus % mu = 0.0_kr
   else if (alpha * D % XX + (1.0_kr + alpha) * D % YY >= 0.0_kr) then
         !!! Note that since D%XX < 0, we cannot have D%XX = D%YY or (1.0_Kr + alpha) * D%YY + D%XX would be < 0
         !!! so we do not need to worry about the case a_i  = a_j in the computation of H_{ij}
      A % XXXX = E * alpha**2 / (1.0_kr + nu) / (1.0_kr + alpha)
      A % XXYY = E * alpha / (1.0_kr + nu)
      A % YYYY = E * (1.0_kr + alpha) / (1.0_kr + nu)
      A % XYXY = E / (1.0_kr + nu) / (1.0_kr + alpha) * ((1.0_kr + alpha) * D % YY + alpha * D % XX) / (D % YY - D % XX)

      D2EEDPlus % type = MEF90HookesLawTypeFull
      D2EEDMinus % type = MEF90HookesLawTypeFull
      D2EEDPlus % fullTensor = Tens4OSTransform(A, Pinv)
      call MEF90HookeLawIsoLambdaMu2D(D2EEDMinus % fullTensor, HookesLaw % lambda, HookesLaw % mu)
      D2EEDMinus % fullTensor = D2EEDMinus % fullTensor - D2EEDPlus % fullTensor

   else
      D2EEDPlus % type = MEF90HookesLawTypeIsotropic
      D2EEDMinus % type = MEF90HookesLawTypeIsotropic

      D2EEDPlus % YoungsModulus = 0.0_kr
      D2EEDPlus % PoissonRatio = 0.0_kr
      D2EEDPlus % lambda = 0.0_kr
      D2EEDPlus % mu = 0.0_kr

      D2EEDMinus % YoungsModulus = E
      D2EEDMinus % PoissonRatio = nu
      D2EEDMinus % lambda = HookesLaw % lambda
      D2EEDMinus % mu = HookesLaw % mu
   end if
#else
   if (D % XX >= 0.0_kr) then
      D2EEDPlus % type = MEF90HookesLawTypeIsotropic
      D2EEDMinus % type = MEF90HookesLawTypeIsotropic

      D2EEDPlus % YoungsModulus = E
      D2EEDPlus % PoissonRatio = nu
      D2EEDPlus % lambda = HookesLaw % lambda
      D2EEDPlus % mu = HookesLaw % mu

      D2EEDMinus % YoungsModulus = 0.0_kr
      D2EEDMinus % PoissonRatio = 0.0_kr
      D2EEDMinus % lambda = 0.0_kr
      D2EEDMinus % mu = 0.0_kr
   else if (nu * D % XX + D % YY >= 0.0_kr) then
      A % XXXX = 2.0_kr * E * nu**2 / (1.0_kr + nu) / (1.0_kr - 2.0_kr * nu)
      A % XXYY = E * nu / (1.0_kr + nu) / (1.0_kr - 2.0_kr * nu)
      A % XXZZ = A % XXYY
      A % YYYY = E * (1.0_kr - nu) / (1.0_kr + nu) / (1.0_kr - 2.0_kr * nu)
      A % YYZZ = A % XXYY
      A % ZZZZ = A % YYYY
      A % YZYZ = E / (1.0_kr + nu)
      A % XZXZ = E / (1.0_kr + nu) * (D % ZZ + nu * D % XX) / (D % ZZ - D % XX)
      A % XYXY = E / (1.0_kr + nu) * (D % YY + nu * D % XX) / (D % YY - D % XX)

      D2EEDPlus % type = MEF90HookesLawTypeFull
      D2EEDMinus % type = MEF90HookesLawTypeFull
      D2EEDPlus % fullTensor = Tens4OSTransform(A, Pinv)
      call MEF90HookeLawIsoLambdaMu3D(D2EEDMinus % fullTensor, HookesLaw % lambda, HookesLaw % mu)
      D2EEDMinus % fullTensor = D2EEDMinus % fullTensor - D2EEDPlus % fullTensor
   else if (nu * (D % XX + D % YY) + (1.0_kr - nu) * D % ZZ >= 0.0_kr) then
      A % XXXX = E * nu**2 / (1.0_kr - nu**2) / (1.0_kr - 2.0_kr * nu)
      A % XXYY = A % XXXX
      A % XXZZ = E * nu / (1.0_kr + nu) / (1.0_kr - 2.0_kr * nu)
      A % YYYY = A % XXXX
      A % YYZZ = A % XXZZ
      A % ZZZZ = E * (1.0_kr - nu) / (1.0_kr + nu) / (1.0_kr - 2.0_kr * nu)
      A % YZYZ = E * nu / (1.0_kr - nu**2) * (nu * (D % XX + D % YY) + (1.0_kr - nu) * D % XX) / (D % ZZ - D % YY)
      A % XZXZ = E * nu / (1.0_kr - nu**2) * (nu * (D % XX + D % YY) + (1.0_kr - nu) * D % XX) / (D % ZZ - D % XX)
      A % XYXY = 0.0_kr

      D2EEDPlus % type = MEF90HookesLawTypeFull
      D2EEDMinus % type = MEF90HookesLawTypeFull
      D2EEDPlus % fullTensor = Tens4OSTransform(A, Pinv)
      call MEF90HookeLawIsoLambdaMu3D(D2EEDMinus % fullTensor, HookesLaw % lambda, HookesLaw % mu)
      D2EEDMinus % fullTensor = D2EEDMinus % fullTensor - D2EEDPlus % fullTensor
   else
      D2EEDPlus % type = MEF90HookesLawTypeIsotropic
      D2EEDMinus % type = MEF90HookesLawTypeIsotropic
      D2EEDPlus % YoungsModulus = 0.0_kr
      D2EEDPlus % PoissonRatio = 0.0_kr
      D2EEDPlus % lambda = 0.0_kr
      D2EEDPlus % mu = 0.0_kr

      D2EEDMinus % YoungsModulus = E
      D2EEDMinus % PoissonRatio = nu
      D2EEDMinus % lambda = HookesLaw % lambda
      D2EEDMinus % mu = HookesLaw % mu
   end if
#endif
end subroutine D2EEDMasonry

end module MEF90_APPEND(m_MEF90_DefMechSplitMasonry,MEF90_DIM)D
