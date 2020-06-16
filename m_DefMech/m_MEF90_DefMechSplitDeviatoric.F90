#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
Module MEF90_APPEND(m_MEF90_DefMechSplitDeviatoric,MEF90_DIM)D
#define MEF90_DEFMECHSPLITDEVIATORIC_CONSTRUCTOR MEF90_APPEND(m_MEF90_DefMechSplitDeviatoric_Constructor,MEF90_DIM)D
#include "finclude/petscdef.h"

   Use m_MEF90
   Use MEF90_APPEND(m_MEF90_DefMechSplit_class,MEF90_DIM)D
   implicit none

   Type, extends(MEF90_DEFMECHSPLIT)                   :: MEF90_DEFMECHSPLITDEVIATORIC
   Contains
      Procedure, pass(self)                            :: EED   => EEDDeviatoric
      Procedure, pass(self)                            :: DEED  => DEEDDeviatoric
      Procedure, pass(self)                            :: D2EED => D2EEDDeviatoric
   End Type

   interface MEF90_DEFMECHSPLITDEVIATORIC
      module procedure MEF90_DEFMECHSPLITDEVIATORIC_CONSTRUCTOR
   end interface

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_DEFMECHSPLITDEVIATORIC_CONSTRUCTOR"
!!!
!!!  
!!!  MEF90_DEFMECHSPLITDEVIATORIC_CONSTRUCTOR: the default constructor for a MEF90_DEFMECHSPLITDeviatoric
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Type(MEF90_DEFMECHSPLITDEVIATORIC) Function MEF90_DEFMECHSPLITDEVIATORIC_CONSTRUCTOR()

      MEF90_DEFMECHSPLITDEVIATORIC_CONSTRUCTOR%damageOrder = 3
      MEF90_DEFMECHSPLITDEVIATORIC_CONSTRUCTOR%strainOrder = 2
      MEF90_DEFMECHSPLITDEVIATORIC_CONSTRUCTOR%type        = 'MEF90DefMech_unilateralContactTypeDeviatoric'
   End Function MEF90_DEFMECHSPLITDEVIATORIC_CONSTRUCTOR



#undef __FUNCT__
#define __FUNCT__ "EEDDeviatoric"
!!!
!!!  
!!!  EEDDeviatoric: Compute the positive and negative part of the elastic energy density associated with a strain tensor Epsilon
!!!    
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!    
   Subroutine EEDDeviatoric(self,Strain,HookesLaw,EEDPlus,EEDMinus)
      Class(MEF90_DEFMECHSPLITDEVIATORIC),Intent(IN)   :: self
      Type(MEF90_MATS),Intent(IN)                      :: Strain
      Type(MEF90_HOOKESLAW),Intent(IN)                 :: HookesLaw
      PetscReal, Intent(OUT)                           :: EEDPlus,EEDMinus

      PetscErrorCode                                   :: ierr
      Character(len=MEF90_MXSTRLEN)                    :: IOBuffer

      If (HookesLaw%type /= MEF90HookesLawTypeIsotropic) Then
         Write(IOBuffer,*) "Deviatoric projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
         Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr)
         SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,IOBuffer,ierr)
      End If

      EEDMinus  = (HookesLaw%lambda + 2.0_Kr * HookesLaw%mu / MEF90_DIM)  * 0.5_Kr ! Ae^s.e^s /2
      EEDPlus   = ((HookesLaw * Strain) .dotP. Strain) * 0.5_Kr - EEDMinus
   End Subroutine

#undef __FUNCT__
#define __FUNCT__ "DEEDDeviatoric"
!!!
!!!  
!!!  DEEDDeviatoric: Compute the derivative of the positive and negative part of the elastic energy density (positive and negative stress)
!!!               evaluated at the strain tensor Strain.
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine DEEDDeviatoric(self,Strain,HookesLaw,DEEDPlus,DEEDMinus)
      Class(MEF90_DEFMECHSPLITDEVIATORIC),Intent(IN)   :: self
      Type(MEF90_MATS),Intent(IN)                      :: Strain
      Type(MEF90_HOOKESLAW),Intent(IN)                 :: HookesLaw
      Type(MEF90_MATS),Intent(OUT)                     :: DEEDPlus,DEEDMinus

      PetscErrorCode                                   :: ierr
      Character(len=MEF90_MXSTRLEN)                    :: IOBuffer

      If (HookesLaw%type /= MEF90HookesLawTypeIsotropic) Then
         Write(IOBuffer,*) "Deviatoric projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
         Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr)
         SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,IOBuffer,ierr)
      End If

      DEEDMinus = ((HookesLaw%lambda + 2.0_Kr * HookesLaw%mu / MEF90_DIM) * 0.5_Kr) * MEF90_MATS_IDENTITY
      DEEDPlus  =  (HookesLaw * Strain) - DEEDMinus
      Call PetscLogFlops(6._pflop,ierr);CHKERRQ(ierr)
   End Subroutine

#undef __FUNCT__
#define __FUNCT__ "D2EEDDeviatoric"
!!!
!!!  
!!!  D2EEDDeviatoric: Compute the derivative of the positive and negative part of the elastic energy density (positive and negative stress)
!!!               evaluated at the strain tensor Strain.
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine D2EEDDeviatoric(self,Strain,HookesLaw,D2EEDPlus,D2EEDMinus)
      Class(MEF90_DEFMECHSPLITDEVIATORIC),Intent(IN)   :: self
      Type(MEF90_MATS),Intent(IN)                      :: Strain
      Type(MEF90_HOOKESLAW),Intent(IN)                 :: HookesLaw
      Type(MEF90_HOOKESLAW),Intent(OUT)                :: D2EEDPlus,D2EEDMinus

      PetscErrorCode                                   :: ierr
      Character(len=MEF90_MXSTRLEN)                    :: IOBuffer

      If (HookesLaw%type /= MEF90HookesLawTypeIsotropic) Then
         Write(IOBuffer,*) "Deviatoric projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
         Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr)
         SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,IOBuffer,ierr)
      End If

      D2EEDPlus%fullTensor  = 0.0_Kr
      D2EEDMinus%fullTensor = 0.0_Kr
      D2EEDPlus%Type    = MEF90HookesLawTypeIsotropic
      D2EEDMinus%Type   = MEF90HookesLawTypeIsotropic
      D2EEDMinus%lambda = (HookesLaw%lambda + 2.0_Kr * HookesLaw%mu / MEF90_DIM) * 0.5_Kr
      D2EEDMinus%mu     = 0.0_Kr
#if MEF90_DIM == 2
      If (HookesLaw%isPlaneStress) Then
         D2EEDMinus%PoissonRatio  = D2EEDMinus%lambda / (D2EEDMinus%lambda + D2EEDMinus%mu) * 0.5_Kr
         D2EEDMinus%YoungsModulus = 2.0_Kr * D2EEDMinus%mu * (1.0_Kr + D2EEDMinus%PoissonRatio)
         D2EEDMinus%BulkModulus   = D2EEDMinus%lambda + D2EEDMinus%mu
         Call PetscLogFlops(14._pflop,ierr);CHKERRQ(ierr)
      Else
         D2EEDMinus%PoissonRatio  = D2EEDMinus%lambda / (D2EEDMinus%lambda + 2.0_Kr * D2EEDMinus%mu) * 0.5_Kr
         D2EEDMinus%YoungsModulus = 2.0_Kr * D2EEDMinus%mu * (1.0_Kr + D2EEDMinus%PoissonRatio)
         D2EEDMinus%BulkModulus   = D2EEDMinus%lambda + D2EEDMinus%mu
         Call PetscLogFlops(17._pflop,ierr);CHKERRQ(ierr)
      End If
#else
      D2EEDMinus%PoissonRatio  = D2EEDMinus%lambda / (D2EEDMinus%lambda + D2EEDMinus%mu) * 0.5_Kr
      D2EEDMinus%YoungsModulus = D2EEDMinus%mu * (3.0_Kr * D2EEDMinus%lambda + 2.0_Kr * D2EEDMinus%mu) / (D2EEDMinus%lambda + D2EEDMinus%mu)
      D2EEDMinus%BulkModulus   = D2EEDMinus%lambda + D2EEDMinus%mu * 2.0_Kr / 3.0_Kr
      Call PetscLogFlops(19._pflop,ierr);CHKERRQ(ierr)
#endif
      D2EEDPlus = HookesLaw - D2EEDMinus
   End Subroutine
End Module MEF90_APPEND(m_MEF90_DefMechSplitDeviatoric,MEF90_DIM)D
