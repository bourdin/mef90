#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
module m_MEF90_DefMechAT1exp
#include "finclude/petscdef.h"
   Use m_MEF90
   Use m_MEF90_DefMechAT_class
   implicit none


!!! AT1exp, a variant of AT1 model with an exponential stiffness interpolation 
!!! function:
!!!
!!! a_b(s)  = 1 + (e^{-bs} - 1) / (1 - e^-b) if b /= 0
!!! a_0(s)  = 1-s
!!!
!!! a is convex if b > 1 and
!!! a_b'(0) = -b / (1-e^{-b}) < -2 if b < 1.5
!!!
   Type, extends(MEF90_DefMechAT_Type)                 :: MEF90_DefMechAT1exp_Type
      PetscReal                                        :: b
   Contains
      Procedure, pass(self)                            :: a   => aAT1exp
      Procedure, pass(self)                            :: Da  => DaAT1exp
      Procedure, pass(self)                            :: D2a => D2aAT1exp

      Procedure, pass(self)                            :: w   => wAT1exp
      Procedure, pass(self)                            :: Dw  => DwAT1exp
      Procedure, pass(self)                            :: D2w => D2wAT1exp
   end Type

   interface MEF90_DefMechAT1exp_Type
      module procedure MEF90_DefMechAT1exp_Constructor
   end interface

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_DefMechAT1exp_Constructor"
!!!
!!!  
!!!  MEF90_DefMechAT1exp_Constructor: the default constructor for a MEF90_DefMechAT1exp_Type
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Type(MEF90_DefMechAT1exp_Type) Function MEF90_DefMechAT1exp_Constructor(b)
      PetscReal,Intent(IN)                             :: b
      
      MEF90_DefMechAT1exp_Constructor%b                 = b
      MEF90_DefMechAT1exp_Constructor%cw                = 2.0_Kr / 3.0_Kr
      MEF90_DefMechAT1exp_Constructor%aorder            = 2
      MEF90_DefMechAT1exp_Constructor%worder            = 1
      MEF90_DefMechAT1exp_Constructor%type              = 'MEF90_DefMechAT1exp'
   End Function MEF90_DefMechAT1exp_Constructor

#undef __FUNCT__
#define __FUNCT__ "aAT1exp"
!!!
!!!  
!!!  aAT1exp: the "a" function of the standard AT1exp model, i.e. a(\alpha) = (1-\alpha)^2
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function aAT1exp(self,alpha)
      Class(MEF90_DefMechAT1exp_Type),Intent(IN)       :: self
      PetscReal                                        :: alpha

      PetscLogDouble                                   :: flops
      PetscErrorCode                                   :: ierr

      If (self%b == 0.0_Kr) Then
         aAT1exp = 1.0_Kr - alpha
         flops = 1.0
      Else
         aAT1exp = 1.0_Kr + (exp(-self%b * alpha) - 1.0_Kr) / (1.0_kr - exp(-self%b))
         flops = 9.0
      EndIf
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End function aAT1exp

#undef __FUNCT__
#define __FUNCT__ "DaAT1exp"
!!!
!!!  
!!!  DaAT1exp: the derivative of the "a" function of the standard AT1exp model, i.e. a(\alpha) = (1-\alpha)^2
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function DaAT1exp(self,alpha)
      Class(MEF90_DefMechAT1exp_Type),Intent(IN)       :: self
      PetscReal                                        :: alpha

      PetscLogDouble                                   :: flops
      PetscErrorCode                                   :: ierr

      If (self%b == 0.0_Kr) Then
         DaAT1exp = -1.0_Kr
         flops = 0.0
      Else
         DaAT1exp = -self%b * exp(-self%b * alpha) / (1.0_kr - exp(-self%b))
         flops = 9.0
      End If
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End function DaAT1exp

#undef __FUNCT__
#define __FUNCT__ "D2aAT1exp"
!!!
!!!  
!!!  D2aAT1exp: the second derivative of the "a" function of the standard AT1exp model, i.e. a(\alpha) = (1-\alpha)^2
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function D2aAT1exp(self,alpha)
      Class(MEF90_DefMechAT1exp_Type),Intent(IN)       :: self
      PetscReal                                        :: alpha

      PetscLogDouble                                   :: flops
      PetscErrorCode                                   :: ierr
      If (self%b == 0.0_Kr) Then
         D2aAT1exp = 0.0_Kr
         flops = 0.0
      Else
         D2aAT1exp = self%b**2 * exp(-self%b * alpha) / (1.0_kr - exp(-self%b))
         flops = 9.0
      End If
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End function D2aAT1exp

#undef __FUNCT__
#define __FUNCT__ "wAT1exp"
!!!
!!!  
!!!  wAT1exp: the "w" function of the standard AT1exp model, i.e. w(\alpha) = \alpha
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function wAT1exp(self,alpha)
      Class(MEF90_DefMechAT1exp_Type),Intent(IN)       :: self
      PetscReal                                        :: alpha

      wAT1exp = alpha
   End function wAT1exp

#undef __FUNCT__
#define __FUNCT__ "DwAT1exp"
!!!
!!!  
!!!  DwAT1exp: the derivative of the "w" function of the standard AT1exp model, i.e. w(\alpha) = \alpha
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function DwAT1exp(self,alpha)
      Class(MEF90_DefMechAT1exp_Type),Intent(IN)       :: self
      PetscReal                                        :: alpha

      DwAT1exp = 1.0_Kr
   End function DwAT1exp

#undef __FUNCT__
#define __FUNCT__ "D2wAT1exp"
!!!
!!!  
!!!  D2wAT1exp: the second derivative of the "w" function of the standard AT1exp model, i.e. w(\alpha) = \alpha
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function D2wAT1exp(self,alpha)
      Class(MEF90_DefMechAT1exp_Type),Intent(IN)       :: self
      PetscReal                                        :: alpha

      D2wAT1exp = 0.0_Kr
   End function D2wAT1exp
End module m_MEF90_DefMechAT1exp
