#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
Module MEF90_APPEND(m_MEF90_DefMechHydrostaticDeviatoric,MEF90_DIM)D
#include "finclude/petscdef.h"
   Use m_MEF90
   Use MEF90_APPEND(m_MEF90_DefMechSplit_class,MEF90_DIM)D
   implicit none

   Type, extends(MEF90_DEFMECHSPLIT)                   :: MEF90_DEFMECHHYDROSTATICDEVIATORICSPLIT
   PetscReal                                           :: gamma
   Contains
      Procedure, pass(self)                            :: EED   => EEDHydrostaticDeviatoric
      Procedure, pass(self)                            :: DEED  => DEEDHydrostaticDeviatoric
      Procedure, pass(self)                            :: D2EED => D2EEDHydrostaticDeviatoric
   End Type

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90HydrostaticDeviatoricPenaltyFunction"
!!!
!!!  
!!!  MEF90HydrostaticDeviatoricPenaltyFunction: a low order polynomial C^2 regularization of (max(x,0))**2, defined by
!!!       0                         if x \le 0
!!!       x^3/3/\gamma             if 0 < x \le \gamma
!!!       (x-\gamma/2)^2+gamma^2/12 otherwise
!!!    
!!!  (c) 2018 Blaise Bourdin bourdin@lsu.edu
!!!

   Pure Function MEF90HydrostaticDeviatoricPenaltyFunction(x)
      PetscReal,Intent(IN) :: x
      PetscReal            :: MEF90HydrostaticDeviatoricPenaltyFunction

      If (x <= 0.0_Kr) Then
         MEF90HydrostaticDeviatoricPenaltyFunction = 0.0_Kr
      Else if (x <= MEF90_HDRegularization) Then
         MEF90HydrostaticDeviatoricPenaltyFunction = x**3/3.0_Kr / MEF90_HDRegularization
      Else
         MEF90HydrostaticDeviatoricPenaltyFunction = (x-MEF90_HDRegularization/2.0_Kr)**2 + MEF90_HDRegularization**2/12.0_Kr
      End If
   End Function MEF90HydrostaticDeviatoricPenaltyFunction

#undef __FUNCT__
#define __FUNCT__ "MEF90HydrostaticDeviatoricDPenaltyFunction"
!!!
!!!  
!!!  MEF90HydrostaticDeviatoricDPenaltyFunction: the first derivative of MEF90HydrostaticDeviatoricPenaltyFunction
!!!    
!!!  (c) 2018 Blaise Bourdin bourdin@lsu.edu
!!!

   Pure Function MEF90HydrostaticDeviatoricDPenaltyFunction(x)
      PetscReal,Intent(IN) :: x
      PetscReal            :: MEF90HydrostaticDeviatoricDPenaltyFunction

      If (x <= 0.0_Kr) Then
         MEF90HydrostaticDeviatoricDPenaltyFunction = 0.0_Kr
      Else if (x <= MEF90_HDRegularization) Then
         MEF90HydrostaticDeviatoricDPenaltyFunction = x**2 / MEF90_HDRegularization
      Else
         MEF90HydrostaticDeviatoricDPenaltyFunction = 2.0_Kr * x - MEF90_HDRegularization
      End If
   End Function MEF90HydrostaticDeviatoricDPenaltyFunction

#undef __FUNCT__
#define __FUNCT__ "MEF90HydrostaticDeviatoricD2PenaltyFunction"
!!!
!!!  
!!!  MEF90HydrostaticDeviatoricD2PenaltyFunction: the second derivative of MEF90HydrostaticDeviatoricPenaltyFunction
!!!    
!!!  (c) 2018 Blaise Bourdin bourdin@lsu.edu
!!!

   Pure Function MEF90HydrostaticDeviatoricD2PenaltyFunction(x)
      PetscReal,Intent(IN) :: x
      PetscReal            :: MEF90HydrostaticDeviatoricD2PenaltyFunction

      If (x <= 0.0_Kr) Then
         MEF90HydrostaticDeviatoricD2PenaltyFunction = 0.0_Kr
      Else if (x <= MEF90_HDRegularization) Then
         MEF90HydrostaticDeviatoricD2PenaltyFunction = 2.0_Kr * x / MEF90_HDRegularization
      Else
         MEF90HydrostaticDeviatoricD2PenaltyFunction = 2.0_Kr 
      End If
   End Function MEF90HydrostaticDeviatoricD2PenaltyFunction

End Module MEF90_APPEND(m_MEF90_DefMechHydrostaticDeviatoric,MEF90_DIM)D
