#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
module m_MEF90_DefMechAT1exp
#include "petsc/finclude/petsc.h"
   ! Use m_MEF90
   use m_MEF90_DefMechAT_class
   implicit none(type, external)
   private
   public :: MEF90DefMechAT1exp_Type

!!! AT1exp, a variant of AT1 model with an exponential stiffness interpolation
!!! function:
!!!
!!! a_b(s)  = 1 + (e^{-bs} - 1) / (1 - e^-b) if b /= 0
!!! a_0(s)  = 1-s
!!!
!!! a is convex if b > 1 and
!!! a_b'(0) = -b / (1-e^{-b}) < -2 if b < 1.5
!!!
   type, extends(MEF90DefMechAT_Type)                  :: MEF90DefMechAT1exp_Type
      PetscReal                                        :: b
   contains
      procedure, pass(self)                            :: a => aAT1exp
      procedure, pass(self)                            :: Da => DaAT1exp
      procedure, pass(self)                            :: D2a => D2aAT1exp

      procedure, pass(self)                            :: w => wAT1exp
      procedure, pass(self)                            :: Dw => DwAT1exp
      procedure, pass(self)                            :: D2w => D2wAT1exp
   end type MEF90DefMechAT1exp_Type

   interface MEF90DefMechAT1exp_Type
      module procedure MEF90DefMechAT1exp_Constructor
   end interface

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechAT1exp_Constructor"
!!!
!!!
!!!  MEF90DefMechAT1exp_Constructor: the default constructor for a MEF90_DefMechAT1exp_Type
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   type(MEF90DefMechAT1exp_Type) function MEF90DefMechAT1exp_Constructor(b)
      PetscReal, intent(IN)                             :: b

      MEF90DefMechAT1exp_Constructor%b = b
      MEF90DefMechAT1exp_Constructor%cw = 2.0_kr / 3.0_kr
      MEF90DefMechAT1exp_Constructor%aorder = 2
      MEF90DefMechAT1exp_Constructor%worder = 1
      MEF90DefMechAT1exp_Constructor%type = 'MEF90_DefMechAT1exp'
   end function MEF90DefMechAT1exp_Constructor

#undef __FUNCT__
#define __FUNCT__ "aAT1exp"
!!!
!!!
!!!  aAT1exp: the "a" function of the standard AT1exp model, i.e. a(\alpha) = (1-\alpha)^2
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function aAT1exp(self, alpha)
      class(MEF90DefMechAT1exp_Type), intent(IN)        :: self
      PetscReal                                        :: alpha

      if (self%b == 0.0_kr) then
         aAT1exp = 1.0_kr - alpha
      else
         aAT1exp = 1.0_kr + (exp(-self%b * alpha) - 1.0_kr) / (1.0_kr - exp(-self%b))
      end if
   end function aAT1exp

#undef __FUNCT__
#define __FUNCT__ "DaAT1exp"
!!!
!!!
!!!  DaAT1exp: the derivative of the "a" function of the standard AT1exp model, i.e. a(\alpha) = (1-\alpha)^2
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function DaAT1exp(self, alpha)
      class(MEF90DefMechAT1exp_Type), intent(IN)        :: self
      PetscReal                                        :: alpha

      if (self%b == 0.0_kr) then
         DaAT1exp = -1.0_kr
      else
         DaAT1exp = -self%b * exp(-self%b * alpha) / (1.0_kr - exp(-self%b))
      end if
   end function DaAT1exp

#undef __FUNCT__
#define __FUNCT__ "D2aAT1exp"
!!!
!!!
!!!  D2aAT1exp: the second derivative of the "a" function of the standard AT1exp model, i.e. a(\alpha) = (1-\alpha)^2
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function D2aAT1exp(self, alpha)
      class(MEF90DefMechAT1exp_Type), intent(IN)        :: self
      PetscReal                                        :: alpha

      if (self%b == 0.0_kr) then
         D2aAT1exp = 0.0_kr
      else
         D2aAT1exp = self%b**2 * exp(-self%b * alpha) / (1.0_kr - exp(-self%b))
      end if
   end function D2aAT1exp

#undef __FUNCT__
#define __FUNCT__ "wAT1exp"
!!!
!!!
!!!  wAT1exp: the "w" function of the standard AT1exp model, i.e. w(\alpha) = \alpha
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function wAT1exp(self, alpha)
      class(MEF90DefMechAT1exp_Type), intent(IN)         :: self
      PetscReal                                        :: alpha

      wAT1exp = alpha
   end function wAT1exp

#undef __FUNCT__
#define __FUNCT__ "DwAT1exp"
!!!
!!!
!!!  DwAT1exp: the derivative of the "w" function of the standard AT1exp model, i.e. w(\alpha) = \alpha
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function DwAT1exp(self, alpha)
      class(MEF90DefMechAT1exp_Type), intent(IN)        :: self
      PetscReal                                        :: alpha

      DwAT1exp = 1.0_kr
   end function DwAT1exp

#undef __FUNCT__
#define __FUNCT__ "D2wAT1exp"
!!!
!!!
!!!  D2wAT1exp: the second derivative of the "w" function of the standard AT1exp model, i.e. w(\alpha) = \alpha
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function D2wAT1exp(self, alpha)
      class(MEF90DefMechAT1exp_Type), intent(IN)        :: self
      PetscReal                                        :: alpha

      D2wAT1exp = 0.0_kr
   end function D2wAT1exp
end module m_MEF90_DefMechAT1exp
