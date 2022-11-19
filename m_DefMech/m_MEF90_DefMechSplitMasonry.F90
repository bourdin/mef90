#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
Module MEF90_APPEND(m_MEF90_DefMechSplitMasonry,MEF90_DIM)D
#include "petsc/finclude/petsc.h"
   Use m_MEF90
   Use MEF90_APPEND(m_MEF90_DefMechSplit_class,MEF90_DIM)D
#define MEF90_DEFMECHSPLITMASONBY_CONSTRUCTOR MEF90_APPEND(m_MEF90_DefMechSplitMasonry_Constructor,MEF90_DIM)D
   implicit none

   Type, extends(MEF90_DEFMECHSPLIT)                   :: MEF90_DEFMECHSPLITMASONRY
   Contains
      Procedure, pass(self)                            :: EED   => EEDMasonry
      Procedure, pass(self)                            :: DEED  => DEEDMasonry
      Procedure, pass(self)                            :: D2EED => D2EEDMasonry
   End Type

   interface MEF90_DEFMECHSPLITMASONRY
      module procedure MEF90_DEFMECHSPLITMASONRY_CONSTRUCTOR
   end interface

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_DEFMECHSPLITMASONRY_CONSTRUCTOR"
!!!
!!!  
!!!  MEF90_DEFMECHSPLITMASONRY_CONSTRUCTOR: the default constructor for a MEF90_DEFMECHSPLITMASONRY
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Type(MEF90_DEFMECHSPLITMASONRY) Function MEF90_DEFMECHSPLITMASONRY_CONSTRUCTOR()
      MEF90_DEFMECHSPLITMASONRY_CONSTRUCTOR%damageOrder  = 0
      MEF90_DEFMECHSPLITMASONRY_CONSTRUCTOR%strainOrder  = 2
      MEF90_DEFMECHSPLITMASONRY_CONSTRUCTOR%type         = 'MEF90DefMech_unilateralContactTypeMasonry'
   End Function MEF90_DEFMECHSPLITMASONRY_CONSTRUCTOR

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
   Subroutine EEDMasonry(self,Strain,HookesLaw,EEDPlus,EEDMinus)
      Class(MEF90_DEFMECHSPLITMASONRY),Intent(IN)        :: self
      Type(MEF90_MATS),Intent(IN)                        :: Strain
      Type(MEF90_HOOKESLAW),Intent(IN)                   :: HookesLaw
      PetscReal, Intent(OUT)                             :: EEDPlus,EEDMinus

      Type(MEF90_MATS)                                   :: D,DPlus
      Type(MEF90_MAT)                                    :: Pinv
      PetscReal                                          :: nu
#if MEF90_DIM == 2
      PetscReal                                          :: alpha
#endif
      PetscErrorCode                                     :: ierr
      Character(len=MEF90MXSTRLEN)                       :: IOBuffer

      If (HookesLaw%type /= MEF90HookesLawTypeIsotropic) Then
         Write(IOBuffer,*) "Masonry projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
         PetscCall(PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr))
         SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,IOBuffer)
      End If
      Call Diagonalize(Strain,Pinv,D)
      !!! D is the strain tensor in the principal basis

      nu = HookesLaw%PoissonRatio
      DPlus = 0.0_Kr

#if MEF90_DIM == 2
      If (HookesLaw%isPlaneStress) Then
         alpha = nu / (1.0_Kr - nu)
      Else
         alpha = nu / (1.0_Kr - 2.0_Kr * nu)
      End If
      If (D%XX >= 0.0_Kr) Then
         DPlus  = D
      Else If ((1.0_Kr + alpha) * D%YY + alpha * D%XX >= 0.0_Kr ) Then
         Dplus%YY = alpha / (1.0_Kr + alpha) * D%XX + D%YY
      End If
#else
      If (D%XX >= 0.0_Kr) Then
         DPlus  = D
      Else If (nu * D%XX + D%YY >= 0.0_Kr ) Then
         DPlus%YY = nu * D%XX + D%YY
         DPlus%ZZ = nu * D%XX + D%ZZ
      Else If (nu * (D%XX + D%YY) + (1.0_Kr - nu) * D%ZZ >= 0.0_Kr ) Then
         DPlus%ZZ = nu / (1.0_Kr - nu) * (D%YY + D%XX) + D%ZZ
      End If
#endif
      !!! We compute the Elastic energy density in the principal basis
      EEDPlus  = (HookesLaw * D .dotP. DPlus) * 0.5_Kr
      EEDMinus = (HookesLaw * D .dotP. (D-DPlus)) * 0.5_Kr
   End Subroutine EEDMasonry


#undef __FUNCT__
#define __FUNCT__ "DEEDMasonry"
!!!
!!!  
!!!  DEEDMasonry: Compute the derivative of the positive and negative part of the elastic energy density (positive and negative stress)
!!!               evaluated at the strain tensor Strain.
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine DEEDMasonry(self,Strain,HookesLaw,DEEDPlus,DEEDMinus)
      Class(MEF90_DEFMECHSPLITMASONRY),Intent(IN)        :: self
      Type(MEF90_MATS),Intent(IN)                        :: Strain
      Type(MEF90_HOOKESLAW),Intent(IN)                   :: HookesLaw
      Type(MEF90_MATS),Intent(OUT)                       :: DEEDPlus,DEEDMinus

      Type(MEF90_MATS)                                   :: D,StrainPlus
      Type(MEF90_MAT)                                    :: Pinv
      PetscReal                                          :: E, nu
#if MEF90_DIM == 2
      PetscReal                                          :: alpha
#endif
      PetscErrorCode                                     :: ierr
      Character(len=MEF90MXSTRLEN)                       :: IOBuffer

      If (HookesLaw%type /= MEF90HookesLawTypeIsotropic) Then
         Write(IOBuffer,*) "Masonry projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
         PetscCall(PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr))
         SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,IOBuffer)
      End If

      Call Diagonalize(Strain,Pinv,D)
      !!! D is the strain tensor in the principal basis

      E  = HookesLaw%YoungsModulus 
      nu = HookesLaw%PoissonRatio
      StrainPlus = 0.0_Kr
      DEEDPlus   = 0.0_Kr
      DEEDMinus  = 0.0_Kr
#if MEF90_DIM == 2
      If (HookesLaw%isPlaneStress) Then
         alpha = nu / (1.0_Kr - nu)
      Else
         alpha = nu / (1.0_Kr - 2.0_Kr * nu)
      End If
      If (D%XX >= 0.0_Kr) Then
         DEEDPlus  = HookesLaw * Strain
         DEEDMinus = 0.0_Kr
      Else If ((1.0_Kr + alpha) * D%YY + alpha * D%XX >= 0.0_Kr) Then
         !!! Compute the projection of the strain in the principal basis
         StrainPlus%XX = 0.0_Kr
         StrainPlus%YY = alpha / (1.0_Kr + alpha) * D%XX + D%YY
         !!! Compute Sigma(Strain^+) then pull back to the canonical basis
         DEEDPlus  = MEF90MatRaRt(HookesLaw * StrainPlus,Pinv)
         DEEDMinus = HookesLaw * Strain - DEEDPlus
      Else
         DEEDPlus  = 0.0_Kr
         DEEDMinus = HookesLaw * Strain
      End If
#else
      If (D%XX >= 0.0_Kr) Then
         DEEDPlus  = HookesLaw * Strain
         DEEDMinus = 0.0_Kr
      Else If (nu * D%XX + D%YY >= 0.0_Kr ) Then
         StrainPlus%YY = nu * D%XX + D%YY
         StrainPlus%ZZ = nu * D%XX + D%ZZ
         !!! Compute Sigma(Strain^+) then pull back to the canonical basis
         DEEDPlus  = MEF90MatRaRt(HookesLaw * StrainPlus,Pinv)
         DEEDMinus = HookesLaw * Strain - DEEDPlus
      Else If (nu * (D%XX + D%YY) + (1.0_Kr - nu) * D%ZZ >= 0.0_Kr ) Then
         StrainPlus%ZZ = nu / (1.0_Kr - nu) * (D%XX + D%YY) + D%ZZ
         !!! Compute Sigma(Strain^+) then pull back to the canonical basis
         DEEDPlus  = MEF90MatRaRt(HookesLaw * StrainPlus,Pinv)
         DEEDMinus = HookesLaw * Strain - DEEDPlus
      Else 
         DEEDPlus  = 0.0_Kr
         DEEDMinus = HookesLaw * Strain
      End If
#endif
   End Subroutine DEEDMasonry

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
   Subroutine D2EEDMasonry(self,Strain,HookesLaw,D2EEDPlus,D2EEDMinus)
      Class(MEF90_DEFMECHSPLITMASONRY),Intent(IN)        :: self
      Type(MEF90_MATS),Intent(IN)                        :: Strain
      Type(MEF90_HOOKESLAW),Intent(IN)                   :: HookesLaw
      Type(MEF90_HOOKESLAW),Intent(OUT)                  :: D2EEDPlus,D2EEDMinus

      Type(MEF90_TENS4OS)                                :: A
      Type(MEF90_MATS)                                   :: D
      Type(MEF90_MAT)                                    :: Pinv
      PetscReal                                          :: E, nu
#if MEF90_DIM == 2
      PetscReal                                          :: alpha
#endif
      PetscErrorCode                                     :: ierr
      Character(len=MEF90MXSTRLEN)                       :: IOBuffer

      If (HookesLaw%type /= MEF90HookesLawTypeIsotropic) Then
         Write(IOBuffer,*) "Masonry projection not implemented for non isotropic Hooke laws: "//__FUNCT__//"\n"
         PetscCall(PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr))
         SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,IOBuffer)
      End If

      Call Diagonalize(Strain,Pinv,D)
      !!! D is the strain tensor in the principal basis

      E  = HookesLaw%YoungsModulus
      nu = HookesLaw%PoissonRatio
      A  = 0.0_Kr
      D2EEDPlus%fullTensor  = 0.0_Kr
      D2EEDMinus%fullTensor = 0.0_Kr
#if MEF90_DIM == 2
      If (HookesLaw%isPlaneStress) Then
         alpha = nu / (1.0_Kr - nu)
      Else
         alpha = nu / (1.0_Kr - 2.0_Kr * nu)
      End If
      If (D%XX >= 0.0_Kr) Then
         D2EEDPlus%Type  = MEF90HookesLawTypeIsotropic
         D2EEDMinus%Type = MEF90HookesLawTypeIsotropic
         D2EEDPlus%YoungsModulus  = E
         D2EEDPlus%PoissonRatio   = nu
         D2EEDPlus%lambda         = HookesLaw%lambda 
         D2EEDPlus%mu             = HookesLaw%mu

         D2EEDMinus%YoungsModulus = 0.0_Kr
         D2EEDMinus%PoissonRatio  = 0.0_Kr
         D2EEDMinus%lambda        = 0.0_Kr 
         D2EEDMinus%mu            = 0.0_Kr
      Else If (alpha * D%XX + (1.0_Kr + alpha) * D%YY >= 0.0_Kr ) Then
         !!! Note that since D%XX < 0, we cannot have D%XX = D%YY or (1.0_Kr + alpha) * D%YY + D%XX would be < 0
         !!! so we do not need to worry about the case a_i  = a_j in the computation of H_{ij}
         A%XXXX = E * alpha**2 / (1.0_Kr + nu) / (1.0_Kr + alpha)
         A%XXYY = E * alpha / (1.0_Kr + nu)
         A%YYYY = E * (1.0_Kr + alpha) / (1.0_Kr + nu)
         A%XYXY = E / (1.0_Kr + nu) / (1.0_Kr + alpha) * ((1.0_Kr + alpha) * D%YY + alpha * D%XX) / (D%YY - D%XX)

         D2EEDPlus%Type  = MEF90HookesLawTypeFull
         D2EEDMinus%Type = MEF90HookesLawTypeFull
         D2EEDPlus%fullTensor  = Tens4OSTransform(A,Pinv)
         Call MEF90HookeLawIsoLambdaMu2D(D2EEDMinus%fullTensor,HookesLaw%lambda,HookesLaw%mu)
         D2EEDMinus%fullTensor = D2EEDMinus%fullTensor - D2EEDPlus%fullTensor

      Else
         D2EEDPlus%Type  = MEF90HookesLawTypeIsotropic
         D2EEDMinus%Type = MEF90HookesLawTypeIsotropic

         D2EEDPlus%YoungsModulus = 0.0_Kr
         D2EEDPlus%PoissonRatio  = 0.0_Kr
         D2EEDPlus%lambda        = 0.0_Kr 
         D2EEDPlus%mu            = 0.0_Kr

         D2EEDMinus%YoungsModulus = E
         D2EEDMinus%PoissonRatio  = nu
         D2EEDMinus%lambda        = HookesLaw%lambda 
         D2EEDMinus%mu            = HookesLaw%mu
      End If
#else
      If (D%XX >= 0.0_Kr) Then
         D2EEDPlus%Type  = MEF90HookesLawTypeIsotropic
         D2EEDMinus%Type = MEF90HookesLawTypeIsotropic

         D2EEDPlus%YoungsModulus  = E
         D2EEDPlus%PoissonRatio   = nu
         D2EEDPlus%lambda         = HookesLaw%lambda 
         D2EEDPlus%mu             = HookesLaw%mu

         D2EEDMinus%YoungsModulus = 0.0_Kr
         D2EEDMinus%PoissonRatio  = 0.0_Kr
         D2EEDMinus%lambda        = 0.0_Kr 
         D2EEDMinus%mu            = 0.0_Kr
      Else If (nu * D%XX + D%YY >= 0.0_Kr ) Then
         A%XXXX = 2.0_Kr * E * nu**2 / (1.0_Kr + nu) / (1.0_Kr - 2.0_Kr * nu)
         A%XXYY = E * nu / (1.0_Kr + nu) / (1.0_Kr - 2.0_Kr * nu)
         A%XXZZ = A%XXYY
         A%YYYY = E * (1.0_Kr - nu) / (1.0_Kr + nu) / (1.0_Kr - 2.0_Kr * nu)
         A%YYZZ = A%XXYY
         A%ZZZZ = A%YYYY
         A%YZYZ = E / (1.0_Kr + nu)
         A%XZXZ = E / (1.0_Kr + nu) * (D%ZZ + nu * D%XX) / (D%ZZ - D%XX)
         A%XYXY = E / (1.0_Kr + nu) * (D%YY + nu * D%XX) / (D%YY - D%XX)

         D2EEDPlus%Type  = MEF90HookesLawTypeFull
         D2EEDMinus%Type = MEF90HookesLawTypeFull
         D2EEDPlus%fullTensor  = Tens4OSTransform(A,Pinv)
         Call MEF90HookeLawIsoLambdaMu3D(D2EEDMinus%fullTensor,HookesLaw%lambda,HookesLaw%mu)
         D2EEDMinus%fullTensor = D2EEDMinus%fullTensor - D2EEDPlus%fullTensor
      Else If (nu * (D%XX + D%YY) + (1.0_Kr - nu) * D%ZZ >= 0.0_Kr ) Then
         A%XXXX = E * nu**2 / (1.0_Kr - nu**2) / (1.0_Kr - 2.0_Kr * nu)
         A%XXYY = A%XXXX
         A%XXZZ = E * nu / (1.0_Kr + nu) / (1.0_Kr - 2.0_Kr * nu)
         A%YYYY = A%XXXX
         A%YYZZ = A%XXZZ
         A%ZZZZ = E * (1.0_Kr - nu) / (1.0_Kr + nu) / (1.0_Kr - 2.0_Kr * nu)
         A%YZYZ = E * nu / (1.0_Kr - nu**2) * (nu * (D%XX + D%YY) + (1.0_Kr - nu) * D%XX) / (D%ZZ - D%YY)
         A%XZXZ = E * nu / (1.0_Kr - nu**2) * (nu * (D%XX + D%YY) + (1.0_Kr - nu) * D%XX) / (D%ZZ - D%XX)
         A%XYXY = 0.0_Kr

         D2EEDPlus%Type  = MEF90HookesLawTypeFull
         D2EEDMinus%Type = MEF90HookesLawTypeFull
         D2EEDPlus%fullTensor  = Tens4OSTransform(A,Pinv)
         Call MEF90HookeLawIsoLambdaMu3D(D2EEDMinus%fullTensor,HookesLaw%lambda,HookesLaw%mu)
         D2EEDMinus%fullTensor = D2EEDMinus%fullTensor - D2EEDPlus%fullTensor
      Else
         D2EEDPlus%Type  = MEF90HookesLawTypeIsotropic
         D2EEDMinus%Type = MEF90HookesLawTypeIsotropic
         D2EEDPlus%YoungsModulus = 0.0_Kr
         D2EEDPlus%PoissonRatio  = 0.0_Kr
         D2EEDPlus%lambda        = 0.0_Kr 
         D2EEDPlus%mu            = 0.0_Kr

         D2EEDMinus%YoungsModulus = E
         D2EEDMinus%PoissonRatio  = nu
         D2EEDMinus%lambda        = HookesLaw%lambda 
         D2EEDMinus%mu            = HookesLaw%mu
      End If
#endif
   End Subroutine D2EEDMasonry

End Module MEF90_APPEND(m_MEF90_DefMechSplitMasonry,MEF90_DIM)D
