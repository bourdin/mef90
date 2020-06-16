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
      Character(len=MEF90_MXSTRLEN)                     :: type
      Integer                                           :: damageOrder
      Integer                                           :: strainOrder
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

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_DefMechSplit_SmoothPositiveSquare"
!!!
!!!  
!!!  MEF90_DefMechSplit_SmoothPositiveSquare: a low order polynomial C^2 regularization of (max(x,0))**2, defined by
!!!       0                         if x \le -gamma/2
!!!       (x+\gamma/2)^3/3/\gamma   if -\gamma/2 < x \le \gamma/2
!!!       x^2+gamma^2/12            otherwise
!!!    
!!!  (c) 2018-2020 Blaise Bourdin bourdin@lsu.edu
!!!

   Function MEF90_DefMechSplit_SmoothPositiveSquare(x,gamma)
      PetscReal,Intent(IN)                             :: x
      PetscReal,Intent(IN)                             :: gamma
      PetscReal                                        :: MEF90_DefMechSplit_SmoothPositiveSquare

      PetscErrorCode                                   :: ierr
      PetscReal                                        :: gammaOver2 

      gammaOver2 = gamma * 0.5_Kr
      If (x <= -gammaOver2) Then
         MEF90_DefMechSplit_SmoothPositiveSquare = 0.0_Kr
         Call PetscLogFlops(2._pflop,ierr);CHKERRQ(ierr)
      Else if (x <= gammaOver2) Then
         MEF90_DefMechSplit_SmoothPositiveSquare = (x+gammaOver2)**3 / 3.0_Kr / gamma
         Call PetscLogFlops(6._pflop,ierr);CHKERRQ(ierr)
      Else
         MEF90_DefMechSplit_SmoothPositiveSquare = x**2 + gammaOver2**2/3.0_Kr
         Call PetscLogFlops(6._pflop,ierr);CHKERRQ(ierr)
      End If
   End Function MEF90_DefMechSplit_SmoothPositiveSquare

#undef __FUNCT__
#define __FUNCT__ "MEF90_DefMechSplit_DSmoothPositiveSquare"
!!!
!!!  
!!!  MEF90_DefMechSplit_DSmoothPositiveSquare: the first derivative of MEF90_DefMechSplitHD_PenaltyFunction
!!!    
!!!  (c) 2018-2020 Blaise Bourdin bourdin@lsu.edu
!!!

   Function MEF90_DefMechSplit_DSmoothPositiveSquare(x,gamma)
      PetscReal,Intent(IN)                             :: x
      PetscReal,Intent(IN)                             :: gamma
      PetscReal                                        :: MEF90_DefMechSplit_DSmoothPositiveSquare

      PetscErrorCode                                   :: ierr
      PetscReal                                        :: gammaOver2 

      gammaOver2 = gamma * 0.5_Kr
      If (x <= -gammaOver2) Then
         MEF90_DefMechSplit_DSmoothPositiveSquare = 0.0_Kr
         Call PetscLogFlops(2._pflop,ierr);CHKERRQ(ierr)
      Else if (x <= gammaOver2) Then
         MEF90_DefMechSplit_DSmoothPositiveSquare = (x+gammaOver2)**2 / gamma
         Call PetscLogFlops(6._pflop,ierr);CHKERRQ(ierr)
      Else
         MEF90_DefMechSplit_DSmoothPositiveSquare = 2.0_Kr * x
         Call PetscLogFlops(3._pflop,ierr);CHKERRQ(ierr)
      End If
   End Function MEF90_DefMechSplit_DSmoothPositiveSquare

#undef __FUNCT__
#define __FUNCT__ "MEF90_DefMechSplit_D2SmoothPositiveSquare"
!!!
!!!  
!!!  MEF90_DefMechSplit_D2SmoothPositiveSquare: the second derivative of MEF90_DefMechSplitHD_PenaltyFunction
!!!
!!!  (c) 2018-2020 Blaise Bourdin bourdin@lsu.edu
!!!

   Function MEF90_DefMechSplit_D2SmoothPositiveSquare(x,gamma)
      PetscReal,Intent(IN)                             :: x
      PetscReal,Intent(IN)                             :: gamma
      PetscReal                                        :: MEF90_DefMechSplit_D2SmoothPositiveSquare

      PetscErrorCode                                   :: ierr
      PetscReal                                        :: gammaOver2 

      gammaOver2 = gamma * 0.5_Kr
      If (x <= -gammaOver2) Then
         MEF90_DefMechSplit_D2SmoothPositiveSquare = 0.0_Kr
         Call PetscLogFlops(2._pflop,ierr);CHKERRQ(ierr)
      Else if (x <= gammaOver2) Then
         MEF90_DefMechSplit_D2SmoothPositiveSquare = 1.0_Kr + x / gammaOver2
         Call PetscLogFlops(4._pflop,ierr);CHKERRQ(ierr)
      Else
         MEF90_DefMechSplit_D2SmoothPositiveSquare = 2.0_Kr 
         Call PetscLogFlops(2._pflop,ierr);CHKERRQ(ierr)
      End If
   End Function MEF90_DefMechSplit_D2SmoothPositiveSquare

End Module MEF90_APPEND(m_MEF90_DefMechSplit_class,MEF90_DIM)D

