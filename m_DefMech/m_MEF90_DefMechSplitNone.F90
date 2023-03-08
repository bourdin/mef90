#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
Module MEF90_APPEND(m_MEF90_DefMechSplitNone,MEF90_DIM)D
#include "petsc/finclude/petsc.h"
   Use m_MEF90
   Use MEF90_APPEND(m_MEF90_DefMechSplit_class,MEF90_DIM)D
#define MEF90_DEFMECHSPLITNONE_CONSTRUCTOR MEF90_APPEND(m_MEF90_DefMechSplitNone_Constructor,MEF90_DIM)D
   implicit none

   Type, extends(MEF90_DEFMECHSPLIT)                   :: MEF90_DEFMECHSPLITNONE
   Contains
      Procedure, pass(self)                            :: EED   => EEDNone
      Procedure, pass(self)                            :: DEED  => DEEDNone
      Procedure, pass(self)                            :: D2EED => D2EEDNone
   end Type

   interface MEF90_DEFMECHSPLITNONE
      module procedure MEF90_DEFMECHSPLITNONE_CONSTRUCTOR
   end interface

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_DEFMECHSPLITNONE_CONSTRUCTOR"
!!!
!!!  
!!!  MEF90_DEFMECHSPLITNONE_CONSTRUCTOR: the default constructor for a MEF90_DEFMECHSPLITNONE
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Type(MEF90_DEFMECHSPLITNONE) Function MEF90_DEFMECHSPLITNONE_CONSTRUCTOR()
      MEF90_DEFMECHSPLITNONE_CONSTRUCTOR%damageOrder  = 0
      MEF90_DEFMECHSPLITNONE_CONSTRUCTOR%strainOrder  = 2
      MEF90_DEFMECHSPLITNONE_CONSTRUCTOR%type         = 'MEF90DefMech_unilateralContactTypeNone'
   End Function MEF90_DEFMECHSPLITNONE_CONSTRUCTOR

#undef __FUNCT__
#define __FUNCT__ "EEDNone"
!!!
!!!  
!!!  EEDNone: Compute the positive and negative part of the elastic energy density associated with a strain tensor 
!!!           without a split, we have EEDPlus  = 1/2 HookesLaw Strain \cdot Strain 
!!!                                    EEDMinus = 0
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine EEDNone(self,Strain,HookesLaw,EEDPlus,EEDMinus)
      Class(MEF90_DEFMECHSPLITNONE),Intent(IN)           :: self
      Type(MEF90_MATS),Intent(IN)                        :: Strain
      Type(MEF90_HOOKESLAW),Intent(IN)                   :: HookesLaw
      PetscReal, Intent(OUT)                             :: EEDPlus,EEDMinus

      EEDPlus  = ((HookesLaw * Strain) .dotP. Strain) * 0.5_Kr
      EEDMinus = 0.0_Kr
   End Subroutine EEDNone


#undef __FUNCT__
#define __FUNCT__ "DEEDNone"
!!!
!!!  
!!!  DEEDNone: Compute the derivative of the positive and negative part of the elastic energy density (positive and negative stress)
!!!               evaluated at the strain tensor Strain.
!!!           without a split, we have DEEDPlus  = HookesLaw Strain 
!!!                                    DEEDMinus = 0
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine DEEDNone(self,Strain,HookesLaw,DEEDPlus,DEEDMinus)
      Class(MEF90_DEFMECHSPLITNONE),Intent(IN)           :: self
      Type(MEF90_MATS),Intent(IN)                        :: Strain
      Type(MEF90_HOOKESLAW),Intent(IN)                   :: HookesLaw
      Type(MEF90_MATS),Intent(OUT)                       :: DEEDPlus,DEEDMinus

      DEEDPlus  = HookesLaw * Strain
      DEEDMinus = 0.0_Kr
   End Subroutine DEEDNone

#undef __FUNCT__
#define __FUNCT__ "D2EEDNone"
!!!
!!!  
!!!  D2EEDNone: Compute the second derivative of the positive and negative part of the elastic energy density (positive and negative stress)
!!!               evaluated at the strain tensor Strain.
!!!               without a split, D2EEDPlus = HookesLaw, D2EEDMinus = 0
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine D2EEDNone(self,Strain,HookesLaw,D2EEDPlus,D2EEDMinus)
      Class(MEF90_DEFMECHSPLITNONE),Intent(IN)           :: self
      Type(MEF90_MATS),Intent(IN)                        :: Strain
      Type(MEF90_HOOKESLAW),Intent(IN)                   :: HookesLaw
      Type(MEF90_HOOKESLAW),Intent(OUT)                  :: D2EEDPlus,D2EEDMinus

      D2EEDPlus = HookesLaw

      
#if MEF90_DIM==2
      D2EEDMinus%isPlaneStress = HookesLaw%isPlaneStress
#endif
      D2EEDMinus%type          = HookesLaw%type
      Select Case(HookesLaw%type)
      Case(MEF90HookesLawTypeIsotropic)
         D2EEDMinus%YoungsModulus = 0.0_Kr * HookesLaw%YoungsModulus 
         D2EEDMinus%PoissonRatio  = 0.0_Kr * HookesLaw%PoissonRatio
         D2EEDMinus%lambda        = 0.0_Kr * HookesLaw%lambda
         D2EEDMinus%mu            = 0.0_Kr * HookesLaw%mu
         D2EEDMinus%BulkModulus   = 0.0_Kr * HookesLaw%BulkModulus
      Case(MEF90HookesLawTypeFull)
         D2EEDMinus%FullTensor   = 0.0_Kr * HookesLaw%FullTensor
         D2EEDMinus%fullTensorLocal   = 0.0_Kr * HookesLaw%fullTensorLocal
      End Select !HookesLaw%type
   End Subroutine D2EEDNone

End Module MEF90_APPEND(m_MEF90_DefMechSplitNone,MEF90_DIM)D
