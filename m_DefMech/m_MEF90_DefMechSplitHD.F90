#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
Module MEF90_APPEND(m_MEF90_DefMechSplitHD,MEF90_DIM)D
#define MEF90_DEFMECHSPLITHD_CONSTRUCTOR MEF90_APPEND(m_MEF90_DefMechSplitHD_Constructor,MEF90_DIM)D
#include "petsc/finclude/petsc.h"

   Use m_MEF90
   Use MEF90_APPEND(m_MEF90_DefMechSplit_class,MEF90_DIM)D
   implicit none

   Type, extends(MEF90_DEFMECHSPLIT)                   :: MEF90_DEFMECHSPLITHD
      PetscReal                                        :: gamma
   Contains
      Procedure, pass(self)                            :: EED   => EEDHD
      Procedure, pass(self)                            :: DEED  => DEEDHD
      Procedure, pass(self)                            :: D2EED => D2EEDHD
   End Type

   interface MEF90_DEFMECHSPLITHD
      module procedure MEF90_DEFMECHSPLITHD_CONSTRUCTOR
   end interface

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_DEFMECHSPLITHD_CONSTRUCTOR"
!!!
!!!  
!!!  MEF90_DEFMECHSPLITHD_CONSTRUCTOR: the default constructor for a MEF90_DEFMECHSPLITHD
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Type(MEF90_DEFMECHSPLITHD) Function MEF90_DEFMECHSPLITHD_CONSTRUCTOR(gamma)
      PetscReal,Intent(IN)                             :: gamma

      MEF90_DEFMECHSPLITHD_CONSTRUCTOR%gamma       = gamma
      MEF90_DEFMECHSPLITHD_CONSTRUCTOR%damageOrder = 3
      MEF90_DEFMECHSPLITHD_CONSTRUCTOR%strainOrder = 2
      MEF90_DEFMECHSPLITHD_CONSTRUCTOR%type        = 'MEF90DefMech_unilateralContactTypeHD'
   End Function MEF90_DEFMECHSPLITHD_CONSTRUCTOR



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
   Subroutine EEDHD(self,Strain,HookesLaw,EEDPlus,EEDMinus)
      Class(MEF90_DEFMECHSPLITHD),Intent(IN)           :: self
      Type(MEF90_MATS),Intent(IN)                      :: Strain
      Type(MEF90_HOOKESLAW),Intent(IN)                 :: HookesLaw
      PetscReal, Intent(OUT)                           :: EEDPlus,EEDMinus

      PetscErrorCode                                   :: ierr
      Character(len=MEF90MXSTRLEN)                     :: IOBuffer

      If (HookesLaw%type /= MEF90HookesLawTypeIsotropic) Then
         Write(IOBuffer,*) "Hydrostatic-Deviatoric projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
         PetscCall(PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr))
         SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,IOBuffer)
      End If

      EEDMinus  = MEF90_DefMechSplit_SmoothPositiveSquare(-trace(Strain),self%gamma) * (HookesLaw%lambda + 2.0_Kr * HookesLaw%mu / MEF90_DIM)  * 0.5_Kr
      EEDPlus   = ((HookesLaw * Strain) .dotP. Strain) * 0.5_Kr - EEDMinus
   End Subroutine

#undef __FUNCT__
#define __FUNCT__ "DEEDHD"
!!!
!!!  
!!!  DEEDHD: Compute the derivative of the positive and negative part of the elastic energy density (positive and negative stress)
!!!               evaluated at the strain tensor Strain.
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine DEEDHD(self,Strain,HookesLaw,DEEDPlus,DEEDMinus)
      Class(MEF90_DEFMECHSPLITHD),Intent(IN)           :: self
      Type(MEF90_MATS),Intent(IN)                      :: Strain
      Type(MEF90_HOOKESLAW),Intent(IN)                 :: HookesLaw
      Type(MEF90_MATS),Intent(OUT)                     :: DEEDPlus,DEEDMinus

      PetscErrorCode                                   :: ierr
      Character(len=MEF90MXSTRLEN)                     :: IOBuffer

      If (HookesLaw%type /= MEF90HookesLawTypeIsotropic) Then
         Write(IOBuffer,*) "Hydrostatic-Deviatoric projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
         PetscCall(PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr))
         SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,IOBuffer)
      End If

      DEEDMinus = (-MEF90_DefMechSplit_DSmoothPositiveSquare(-trace(Strain),self%gamma) * (HookesLaw%lambda + 2.0_Kr * HookesLaw%mu / MEF90_DIM) * 0.5_Kr) * MEF90_MATS_IDENTITY
      DEEDPlus  =  (HookesLaw * Strain) - DEEDMinus
   End Subroutine

#undef __FUNCT__
#define __FUNCT__ "D2EEDHD"
!!!
!!!  
!!!  D2EEDHD: Compute the derivative of the positive and negative part of the elastic energy density (positive and negative stress)
!!!               evaluated at the strain tensor Strain.
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine D2EEDHD(self,Strain,HookesLaw,D2EEDPlus,D2EEDMinus)
      Class(MEF90_DEFMECHSPLITHD),Intent(IN)           :: self
      Type(MEF90_MATS),Intent(IN)                      :: Strain
      Type(MEF90_HOOKESLAW),Intent(IN)                 :: HookesLaw
      Type(MEF90_HOOKESLAW),Intent(OUT)                :: D2EEDPlus,D2EEDMinus

      PetscErrorCode                                   :: ierr
      Character(len=MEF90MXSTRLEN)                     :: IOBuffer

      If (HookesLaw%type /= MEF90HookesLawTypeIsotropic) Then
         Write(IOBuffer,*) "Hydrostatic-Deviatoric projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
         PetscCall(PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr))
         SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,IOBuffer)
      End If

      D2EEDPlus%fullTensor  = 0.0_Kr
      D2EEDMinus%fullTensor = 0.0_Kr
      D2EEDPlus%Type    = MEF90HookesLawTypeIsotropic
      D2EEDMinus%Type   = MEF90HookesLawTypeIsotropic
      D2EEDMinus%lambda = MEF90_DefMechSplit_D2SmoothPositiveSquare(-trace(Strain),self%gamma) * (HookesLaw%lambda + 2.0_Kr * HookesLaw%mu / MEF90_DIM) * 0.5_Kr
      D2EEDMinus%mu     = 0.0_Kr
#if MEF90_DIM == 2
      D2EEDMinus%isPlaneStress = HookesLaw%isPlaneStress
      D2EEDMinus%YoungsModulus = 2.0_Kr * D2EEDMinus%mu * (1.0_Kr + D2EEDMinus%PoissonRatio)
      D2EEDMinus%BulkModulus   = D2EEDMinus%lambda + D2EEDMinus%mu
      If (HookesLaw%isPlaneStress) Then
         D2EEDMinus%PoissonRatio  = D2EEDMinus%lambda / (D2EEDMinus%lambda + D2EEDMinus%mu) * 0.5_Kr
      Else
         D2EEDMinus%PoissonRatio  = D2EEDMinus%lambda / (D2EEDMinus%lambda + 2.0_Kr * D2EEDMinus%mu) * 0.5_Kr
      End If
#else
      D2EEDMinus%PoissonRatio  = D2EEDMinus%lambda / (D2EEDMinus%lambda + D2EEDMinus%mu) * 0.5_Kr
      D2EEDMinus%YoungsModulus = D2EEDMinus%mu * (3.0_Kr * D2EEDMinus%lambda + 2.0_Kr * D2EEDMinus%mu) / (D2EEDMinus%lambda + D2EEDMinus%mu)
      D2EEDMinus%BulkModulus   = D2EEDMinus%lambda + D2EEDMinus%mu * 2.0_Kr / 3.0_Kr
#endif
      D2EEDPlus = HookesLaw - D2EEDMinus
   End Subroutine
End Module MEF90_APPEND(m_MEF90_DefMechSplitHD,MEF90_DIM)D
