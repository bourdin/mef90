#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
Module MEF90_APPEND(m_MEF90_DefMechSplitHydrostatic,MEF90_DIM)D
#define MEF90_DEFMECHSPLITHYDROSTATIC_CONSTRUCTOR MEF90_APPEND(m_MEF90_DefMechSplitHydrostatic_Constructor,MEF90_DIM)D
#include "petsc/finclude/petsc.h"

   Use m_MEF90
   Use MEF90_APPEND(m_MEF90_DefMechSplit_class,MEF90_DIM)D
   implicit none

   Type, extends(MEF90_DEFMECHSPLIT)                   :: MEF90_DEFMECHSPLITHYDROSTATIC
      PetscReal                                        :: gamma
   Contains
      Procedure, pass(self)                            :: EED   => EEDHydrostatic
      Procedure, pass(self)                            :: DEED  => DEEDHydrostatic
      Procedure, pass(self)                            :: D2EED => D2EEDHydrostatic
   End Type

   interface MEF90_DEFMECHSPLITHYDROSTATIC
      module procedure MEF90_DEFMECHSPLITHYDROSTATIC_CONSTRUCTOR
   end interface

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_DEFMECHSPLITHYDROSTATIC_CONSTRUCTOR"
!!!
!!!  
!!!  MEF90_DEFMECHSPLITHYDROSTATIC_CONSTRUCTOR: the default constructor for a MEF90_DEFMECHSPLITHydrostatic
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Type(MEF90_DEFMECHSPLITHYDROSTATIC) Function MEF90_DEFMECHSPLITHYDROSTATIC_CONSTRUCTOR(gamma)
      PetscReal,Intent(IN)                             :: gamma

      MEF90_DEFMECHSPLITHYDROSTATIC_CONSTRUCTOR%gamma       = gamma
      MEF90_DEFMECHSPLITHYDROSTATIC_CONSTRUCTOR%damageOrder = 3
      MEF90_DEFMECHSPLITHYDROSTATIC_CONSTRUCTOR%strainOrder = 2
      MEF90_DEFMECHSPLITHYDROSTATIC_CONSTRUCTOR%type        = 'MEF90DefMech_unilateralContactTypeHydrostatic'
   End Function MEF90_DEFMECHSPLITHYDROSTATIC_CONSTRUCTOR



#undef __FUNCT__
#define __FUNCT__ "EEDHydrostatic"
!!!
!!!  
!!!  EEDHydrostatic: Compute the positive and negative part of the elastic energy density associated with a strain tensor Epsilon
!!!    
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!    
   Subroutine EEDHydrostatic(self,Strain,HookesLaw,EEDPlus,EEDMinus)
      Class(MEF90_DEFMECHSPLITHYDROSTATIC),Intent(IN)  :: self
      Type(MEF90_MATS),Intent(IN)                      :: Strain
      Type(MEF90_HOOKESLAW),Intent(IN)                 :: HookesLaw
      PetscReal, Intent(OUT)                           :: EEDPlus,EEDMinus

      PetscErrorCode                                   :: ierr
      Character(len=MEF90MXSTRLEN)                     :: IOBuffer

      If (HookesLaw%type /= MEF90HookesLawTypeIsotropic) Then
         Write(IOBuffer,*) "Hydrostatic projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
         PetscCall(PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr))
         SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,IOBuffer)
      End If

      EEDPlus   = MEF90_DefMechSplit_SmoothPositiveSquare(trace(Strain),self%gamma) * (HookesLaw%lambda + 2.0_Kr * HookesLaw%mu / MEF90_DIM)  * 0.5_Kr
      ! Ae^s.e^s /2
      EEDMinus  = ((HookesLaw * Strain) .dotP. Strain) * 0.5_Kr - EEDPlus
   End Subroutine

#undef __FUNCT__
#define __FUNCT__ "DEEDHydrostatic"
!!!
!!!  
!!!  DEEDHydrostatic: Compute the derivative of the positive and negative part of the elastic energy density (positive and negative stress)
!!!               evaluated at the strain tensor Strain.
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine DEEDHydrostatic(self,Strain,HookesLaw,DEEDPlus,DEEDMinus)
      Class(MEF90_DEFMECHSPLITHYDROSTATIC),Intent(IN)  :: self
      Type(MEF90_MATS),Intent(IN)                      :: Strain
      Type(MEF90_HOOKESLAW),Intent(IN)                 :: HookesLaw
      Type(MEF90_MATS),Intent(OUT)                     :: DEEDPlus,DEEDMinus

      PetscErrorCode                                   :: ierr
      Character(len=MEF90MXSTRLEN)                     :: IOBuffer

      If (HookesLaw%type /= MEF90HookesLawTypeIsotropic) Then
         Write(IOBuffer,*) "Hydrostatic projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
         PetscCall(PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr))
         SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,IOBuffer)
      End If

      DEEDPlus  = MEF90_DefMechSplit_DSmoothPositiveSquare(trace(Strain),self%gamma) * ((HookesLaw%lambda + 2.0_Kr * HookesLaw%mu / MEF90_DIM) * 0.5_Kr) * MEF90_MATS_IDENTITY
      DEEDMinus =  (HookesLaw * Strain) - DEEDPlus
      PetscCall(PetscLogFlops(6._pflop,ierr))
   End Subroutine

#undef __FUNCT__
#define __FUNCT__ "D2EEDHydrostatic"
!!!
!!!  
!!!  D2EEDHydrostatic: Compute the derivative of the positive and negative part of the elastic energy density (positive and negative stress)
!!!               evaluated at the strain tensor Strain.
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine D2EEDHydrostatic(self,Strain,HookesLaw,D2EEDPlus,D2EEDMinus)
      Class(MEF90_DEFMECHSPLITHYDROSTATIC),Intent(IN)  :: self
      Type(MEF90_MATS),Intent(IN)                      :: Strain
      Type(MEF90_HOOKESLAW),Intent(IN)                 :: HookesLaw
      Type(MEF90_HOOKESLAW),Intent(OUT)                :: D2EEDPlus,D2EEDMinus

      PetscErrorCode                                   :: ierr
      Character(len=MEF90MXSTRLEN)                     :: IOBuffer

      If (HookesLaw%type /= MEF90HookesLawTypeIsotropic) Then
         Write(IOBuffer,*) "Hydrostatic projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
         PetscCall(PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr))
         SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,IOBuffer)
      End If

      D2EEDMinus%fullTensor = 0.0_Kr
      D2EEDPlus%fullTensor  = 0.0_Kr
      D2EEDPlus%Type        = MEF90HookesLawTypeIsotropic
      D2EEDPlus%Type        = MEF90HookesLawTypeIsotropic
      D2EEDPlus%lambda      = MEF90_DefMechSplit_D2SmoothPositiveSquare(trace(Strain),self%gamma) * (HookesLaw%lambda + 2.0_Kr * HookesLaw%mu / MEF90_DIM) * 0.5_Kr
      D2EEDPlus%mu          = 0.0_Kr
#if MEF90_DIM == 2
      D2EEDPlus%isPlaneStress = HookesLaw%isPlaneStress
      D2EEDPlus%YoungsModulus = 2.0_Kr * D2EEDMinus%mu * (1.0_Kr + D2EEDMinus%PoissonRatio)
      D2EEDPlus%BulkModulus   = D2EEDMinus%lambda + D2EEDMinus%mu
      If (HookesLaw%isPlaneStress) Then
         D2EEDPlus%PoissonRatio  = D2EEDPlus%lambda / (D2EEDPlus%lambda + D2EEDPlus%mu) * 0.5_Kr
         PetscCall(PetscLogFlops(14._pflop,ierr))
      Else
         D2EEDPlus%PoissonRatio  = D2EEDPlus%lambda / (D2EEDPlus%lambda + 2.0_Kr * D2EEDPlus%mu) * 0.5_Kr
         PetscCall(PetscLogFlops(17._pflop,ierr))
      End If
#else
      D2EEDPlus%PoissonRatio  = D2EEDPlus%lambda / (D2EEDPlus%lambda + D2EEDPlus%mu) * 0.5_Kr
      D2EEDPlus%YoungsModulus = D2EEDPlus%mu * (3.0_Kr * D2EEDPlus%lambda + 2.0_Kr * D2EEDPlus%mu) / (D2EEDPlus%lambda + D2EEDPlus%mu)
      D2EEDPlus%BulkModulus   = D2EEDPlus%lambda + D2EEDPlus%mu * 2.0_Kr / 3.0_Kr
      PetscCall(PetscLogFlops(19._pflop,ierr))
#endif
      D2EEDMinus = HookesLaw - D2EEDPlus
   End Subroutine
End Module MEF90_APPEND(m_MEF90_DefMechSplitHydrostatic,MEF90_DIM)D
