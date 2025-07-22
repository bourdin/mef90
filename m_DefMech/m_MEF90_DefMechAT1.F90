#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
module m_MEF90_DefMechAT1
#include "petsc/finclude/petsc.h"
   ! Use m_MEF90
   use m_MEF90_DefMechAT_class
   implicit none(type, external)
   private
   public :: MEF90DefMechAT1_Type

   type, extends(MEF90DefMechAT_Type)                  :: MEF90DefMechAT1_Type
   contains
      procedure, pass(self)                            :: a => aAT1
      procedure, pass(self)                            :: Da => DaAT1
      procedure, pass(self)                            :: D2a => D2aAT1

      procedure, pass(self)                            :: w => wAT1
      procedure, pass(self)                            :: Dw => DwAT1
      procedure, pass(self)                            :: D2w => D2wAT1
   end type MEF90DefMechAT1_Type

   interface MEF90DefMechAT1_Type
      module procedure MEF90DefMechAT1_Constructor
   end interface

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechAT1_Constructor"
!!!
!!!
!!!  MEF90_DefMechAT1_Constructor: the default constructor for a MEF90_DefMechAT1_Type
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   type(MEF90DefMechAT1_Type) function MEF90DefMechAT1_Constructor()
      MEF90DefMechAT1_Constructor%cw = 2.0_kr / 3.0_kr
      MEF90DefMechAT1_Constructor%aorder = 2
      MEF90DefMechAT1_Constructor%worder = 1
      MEF90DefMechAT1_Constructor%type = 'MEF90DefMechAT1'
   end function MEF90DefMechAT1_Constructor

#undef __FUNCT__
#define __FUNCT__ "aAT1"
!!!
!!!
!!!  aAT1: the "a" function of the standard AT1 model, i.e. a(\alpha) = (1-\alpha)^2
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function aAT1(self, alpha)
      class(MEF90DefMechAT1_Type), intent(IN)           :: self
      PetscReal                                        :: alpha

      aAT1 = (1.0_kr - alpha)**2
   end function aAT1

#undef __FUNCT__
#define __FUNCT__ "DaAT1"
!!!
!!!
!!!  DaAT1: the derivative of the "a" function of the standard AT1 model, i.e. a(\alpha) = (1-\alpha)^2
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function DaAT1(self, alpha)
      class(MEF90DefMechAT1_Type), intent(IN)           :: self
      PetscReal                                        :: alpha

      DaAT1 = -2.0_kr * (1.0_kr - alpha)
   end function DaAT1

#undef __FUNCT__
#define __FUNCT__ "D2aAT1"
!!!
!!!
!!!  D2aAT1: the second derivative of the "a" function of the standard AT1 model, i.e. a(\alpha) = (1-\alpha)^2
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function D2aAT1(self, alpha)
      class(MEF90DefMechAT1_Type), intent(IN)           :: self
      PetscReal                                        :: alpha

      D2aAT1 = 2.0_kr
   end function D2aAT1

#undef __FUNCT__
#define __FUNCT__ "wAT1"
!!!
!!!
!!!  wAT1: the "w" function of the standard AT1 model, i.e. w(\alpha) = \alpha
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function wAT1(self, alpha)
      class(MEF90DefMechAT1_Type), intent(IN)           :: self
      PetscReal                                        :: alpha

      wAT1 = alpha
   end function wAT1

#undef __FUNCT__
#define __FUNCT__ "DwAT1"
!!!
!!!
!!!  DwAT1: the derivative of the "w" function of the standard AT1 model, i.e. w(\alpha) = \alpha
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function DwAT1(self, alpha)
      class(MEF90DefMechAT1_Type), intent(IN)           :: self
      PetscReal                                        :: alpha

      DwAT1 = 1.0_kr
   end function DwAT1

#undef __FUNCT__
#define __FUNCT__ "D2wAT1"
!!!
!!!
!!!  D2wAT1: the second derivative of the "w" function of the standard AT1 model, i.e. w(\alpha) = \alpha
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function D2wAT1(self, alpha)
      class(MEF90DefMechAT1_Type), intent(IN)           :: self
      PetscReal                                        :: alpha

      D2wAT1 = 0.0_kr
   end function D2wAT1
end module m_MEF90_DefMechAT1
