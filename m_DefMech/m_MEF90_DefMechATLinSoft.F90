#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
module m_MEF90_DefMechATLinSoft
#include "petsc/finclude/petsc.h"
   ! Use m_MEF90
   use m_MEF90_DefMechAT_class
   implicit none(type, external)
   private
   public :: MEF90DefMechATLinSoft_Type

!!! LinSoft, a model with linear softening phase
!!! function:
!!!
!!! a_k(s) = (1-s)**2/(k+(1-k)(1-s)**2)
!!! w(s)   = 1-(1-s)^2
!!! cw = pi/4
   type, extends(MEF90DefMechAT_Type)                  :: MEF90DefMechATLinSoft_Type
      PetscReal                                        :: k
   contains
      procedure, pass(self)                            :: a => aLinSoft
      procedure, pass(self)                            :: Da => DaLinSoft
      procedure, pass(self)                            :: D2a => D2aLinSoft

      procedure, pass(self)                            :: w => wLinSoft
      procedure, pass(self)                            :: Dw => DwLinSoft
      procedure, pass(self)                            :: D2w => D2wLinSoft
   end type MEF90DefMechATLinSoft_Type

   interface MEF90DefMechATLinSoft_Type
      module procedure MEF90DefMechATLinSoft_Constructor
   end interface

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechATLinSoft_Constructor"
!!!
!!!
!!!  MEF90DefMechATLinSoft_Constructor: the default constructor for a MEF90_DefMechATLinSoft_Type
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   type(MEF90DefMechATLinSoft_Type) function MEF90DefMechATLinSoft_Constructor(k)
      PetscReal, intent(IN)                                :: k

      MEF90DefMechATLinSoft_Constructor % k = k
      MEF90DefMechATLinSoft_Constructor % cw = PETSC_PI / 4.0_kr
      MEF90DefMechATLinSoft_Constructor % aorder = 2
      MEF90DefMechATLinSoft_Constructor % worder = 2
      MEF90DefMechATLinSoft_Constructor % type = 'MEF90DefMechATLinSoft'
   end function MEF90DefMechATLinSoft_Constructor

#undef __FUNCT__
#define __FUNCT__ "aLinSoft"
!!!
!!!
!!!  aLinSoft: the "a" function of the LinSoft model
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function aLinSoft(self, alpha)
      class(MEF90DefMechATLinSoft_Type), intent(IN)     :: self
      PetscReal                                        :: alpha

      aLinSoft = (1.0_kr - alpha)**2 / (self % k + (1.0_kr - self % k) * (1.0_kr * alpha)**2)
   end function aLinSoft

#undef __FUNCT__
#define __FUNCT__ "DaLinSoft"
!!!
!!!
!!!  DaLinSoft: the derivative of the "a" function of the LinSoft model
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function DaLinSoft(self, alpha)
      class(MEF90DefMechATLinSoft_Type), intent(IN)     :: self
      PetscReal                                        :: alpha

      DaLinSoft = -2.0_kr * alpha * (1.0_kr - alpha)**2 * (self % k - 1.0_kr) / (alpha**2 * (self % k - 1.0_kr) + 1.0_kr)**2 &
                  + (2.0_kr * alpha - 2.0_kr) / (alpha**2 * (self % k - 1.0_kr) + 1.0_kr)
   end function DaLinSoft

#undef __FUNCT__
#define __FUNCT__ "D2aLinSoft"
!!!
!!!
!!!  D2aLinSoft: the second derivative of the "a" function of the LinSoft model
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function D2aLinSoft(self, alpha)
      class(MEF90DefMechATLinSoft_Type), intent(IN)     :: self
      PetscReal                                        :: alpha

      D2aLinSoft = 8.0_kr * alpha**2 * (1.0_kr - alpha)**2 * (self % k - 1.0_kr)**2 / (alpha**2 * (self % k - 1.0_kr) + 1.0_kr)**3 &
                   - 4.0_kr * alpha * (2.0_kr * alpha - 2.0_kr) * (self % k - 1.0_kr) / (alpha**2 * (self % k - 1.0_kr) + 1.0_kr)**2 &
                   - 2.0_kr * (1.0_kr - alpha)**2 * (self % k - 1.0_kr) / (alpha**2 * (self % k - 1.0_kr) + 1.0_kr)**2 &
                   + 2.0_kr / (alpha**2 * (self % k - 1.0_kr) + 1.0_kr)
   end function D2aLinSoft

#undef __FUNCT__
#define __FUNCT__ "wLinSoft"
!!!
!!!
!!!  wLinSoft: the "w" function of the LinSoft model
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function wLinSoft(self, alpha)
      class(MEF90DefMechATLinSoft_Type), intent(IN)     :: self
      PetscReal                                        :: alpha

      wLinSoft = alpha * (2.0_kr - alpha)

   end function wLinSoft

#undef __FUNCT__
#define __FUNCT__ "DwLinSoft"
!!!
!!!
!!!  DwLinSoft: the derivative of the "w" function of the LinSoft model
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function DwLinSoft(self, alpha)
      class(MEF90DefMechATLinSoft_Type), intent(IN)     :: self
      PetscReal                                        :: alpha

      DwLinSoft = 2.0_kr * (1.0_kr - alpha)

   end function DwLinSoft

#undef __FUNCT__
#define __FUNCT__ "D2wLinSoft"
!!!
!!!
!!!  D2wLinSoft: the second derivative of the "w" function of the LinSoft model
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function D2wLinSoft(self, alpha)
      class(MEF90DefMechATLinSoft_Type), intent(IN)     :: self
      PetscReal                                        :: alpha

      D2wLinSoft = -2.0_kr
   end function D2wLinSoft
end module m_MEF90_DefMechATLinSoft
