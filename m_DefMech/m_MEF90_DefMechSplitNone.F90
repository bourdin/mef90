#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
Module MEF90_APPEND(m_MEF90_DefMechNone,MEF90_DIM)D
#include "finclude/petscdef.h"
   Use m_MEF90
   Use MEF90_APPEND(m_MEF90_DefMechSplit_class,MEF90_DIM)D
   implicit none

   Type, extends(MEF90_DEFMECHSPLIT)                   :: MEF90_DEFMECHNONESPLIT
   Contains
      Procedure, pass(self)                            :: EED   => EEDNone
      Procedure, pass(self)                            :: DEED  => DEEDNone
      Procedure, pass(self)                            :: D2EED => D2EEDNone
   end Type

Contains
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
      Class(MEF90_DEFMECHNONESPLIT),Intent(IN)           :: self
      Type(MEF90_MATS),Intent(IN)                        :: Strain
      Type(MEF90_HOOKESLAW),Intent(IN)                   :: HookesLaw
      PetscReal, Intent(OUT)                             :: EEDPlus,EEDMinus

      EEDPlus  = (HookesLaw * Strain .dotP. Strain) * 0.5_Kr
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
      Class(MEF90_DEFMECHNONESPLIT),Intent(IN)           :: self
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
      Class(MEF90_DEFMECHNONESPLIT),Intent(IN)           :: self
      Type(MEF90_MATS),Intent(IN)                        :: Strain
      Type(MEF90_HOOKESLAW),Intent(IN)                   :: HookesLaw
      Type(MEF90_HOOKESLAW),Intent(OUT)                  :: D2EEDPlus,D2EEDMinus

      Type(MEF90_TENS4OS)                                :: A
      Type(MEF90_MATS)                                   :: D
      Type(MEF90_MAT)                                    :: Pinv
      PetscReal                                          :: E, nu,alpha
      PetscErrorCode                                     :: ierr

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
   End Subroutine D2EEDNone

End Module MEF90_APPEND(m_MEF90_DefMechNone,MEF90_DIM)D
