#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
Module MEF90_APPEND(m_MEF90_DefMechSplit_class,MEF90_DIM)D
#include "finclude/petscdef.h"
#define EEDINTERFACE   MEF90_APPEND(EED,MEF90_DIM)D
#define DEEDINTERFACE  MEF90_APPEND(DEED,MEF90_DIM)D
#define D2EEDINTERFACE MEF90_APPEND(D2EED,MEF90_DIM)D

   Use m_MEF90
   Implicit none

!!!
!!!  
!!!  MEF90_DefMechSplit_Type: The abstract class used to define an energy split for
!!!                           handling unilateral contact
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!

   Type, abstract :: MEF90_DEFMECHSPLIT
   Contains
      Procedure(EEDINTERFACE), pass(self), deferred     :: EED
      Procedure(DEEDINTERFACE), pass(self), deferred    :: DEED
      Procedure(D2EEDINTERFACE), pass(self), deferred   :: D2EED
   End Type

   Abstract Interface
      Subroutine EEDINTERFACE(self,Strain,HookesLaw,EEDPlus,EEDMinus)
         Use m_MEF90
         import :: MEF90_DEFMECHSPLIT
         Class(MEF90_DEFMECHSPLIT),Intent(IN)           :: self
         Type(MEF90_MATS),Intent(IN)                    :: Strain
         Type(MEF90_HOOKESLAW),Intent(IN)               :: HookesLaw
         PetscReal, Intent(OUT)                         :: EEDPlus,EEDMinus
      End Subroutine EEDINTERFACE

      Subroutine DEEDINTERFACE(self,Strain,HookesLaw,DEEDPlus,DEEDMinus)
         Use m_MEF90
         import :: MEF90_DEFMECHSPLIT
         Class(MEF90_DEFMECHSPLIT),Intent(IN)           :: self
         Type(MEF90_MATS),Intent(IN)                    :: Strain
         Type(MEF90_HOOKESLAW),Intent(IN)               :: HookesLaw
         Type(MEF90_MATS),Intent(OUT)                   :: DEEDPlus,DEEDMinus
      End Subroutine DEEDINTERFACE

      Subroutine D2EEDINTERFACE(self,Strain,HookesLaw,D2EEDPlus,D2EEDMinus)
         Use m_MEF90
         import :: MEF90_DEFMECHSPLIT
         Class(MEF90_DEFMECHSPLIT),Intent(IN)           :: self
         Type(MEF90_MATS),Intent(IN)                    :: Strain
         Type(MEF90_HOOKESLAW),Intent(IN)               :: HookesLaw
         Type(MEF90_HOOKESLAW),Intent(OUT)              :: D2EEDPlus,D2EEDMinus
      End Subroutine D2EEDINTERFACE
   End Interface
End Module MEF90_APPEND(m_MEF90_DefMechSplit_class,MEF90_DIM)D

