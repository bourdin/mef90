#include "../MEF90/mef90.inc"
module m_MEF90_DefMechAT2
#include "petsc/finclude/petsc.h"
   use m_MEF90_DefMechAT_class
   implicit none(type, external)
   private
   public :: MEF90DefMechAT2_Type

   type, extends(MEF90DefMechAT_Type)                  :: MEF90DefMechAT2_Type
   contains
      procedure, pass(self)                            :: a => aAT2
      procedure, pass(self)                            :: Da => DaAT2
      procedure, pass(self)                            :: D2a => D2aAT2

      procedure, pass(self)                            :: w => wAT2
      procedure, pass(self)                            :: Dw => DwAT2
      procedure, pass(self)                            :: D2w => D2wAT2
   end type MEF90DefMechAT2_Type

   interface MEF90DefMechAT2_Type
      module procedure MEF90DefMechAT2_Constructor
   end interface

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechAT2_Constructor"
!!!
!!!
!!!  MEF90DefMechAT2_Constructor: the default constructor for a MEF90_DefMechAT2_Type
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   type(MEF90DefMechAT2_Type) function MEF90DefMechAT2_Constructor()
      MEF90DefMechAT2_Constructor%cw = 0.5_kr
      MEF90DefMechAT2_Constructor%aorder = 2
      MEF90DefMechAT2_Constructor%worder = 2
      MEF90DefMechAT2_Constructor%type = 'MEF90DefMechAT2'
   end function MEF90DefMechAT2_Constructor

#undef __FUNCT__
#define __FUNCT__ "aAT2"
!!!
!!!
!!!  aAT2: the "a" function of the standard AT2 model, i.e. a(\alpha) = (1-\alpha)^2
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function aAT2(self, alpha)
      class(MEF90DefMechAT2_Type), intent(IN)           :: self
      PetscReal                                        :: alpha

      aAT2 = (1.0_kr - alpha)**2
   end function aAT2

#undef __FUNCT__
#define __FUNCT__ "DaAT2"
!!!
!!!
!!!  DaAT2: the derivative of the "a" function of the standard AT2 model, i.e. a(\alpha) = (1-\alpha)^2
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function DaAT2(self, alpha)
      class(MEF90DefMechAT2_Type), intent(IN)           :: self
      PetscReal                                        :: alpha

      DaAT2 = -2.0_kr * (1.0_kr - alpha)
   end function DaAT2

#undef __FUNCT__
#define __FUNCT__ "D2aAT2"
!!!
!!!
!!!  D2aAT2: the second derivative of the "a" function of the standard AT2 model, i.e. a(\alpha) = (1-\alpha)^2
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function D2aAT2(self, alpha)
      class(MEF90DefMechAT2_Type), intent(IN)           :: self
      PetscReal                                        :: alpha

      D2aAT2 = 2.0_kr
   end function D2aAT2

#undef __FUNCT__
#define __FUNCT__ "wAT2"
!!!
!!!
!!!  wAT2: the "w" function of the standard AT2 model, i.e. w(\alpha) = \alpha
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function wAT2(self, alpha)
      class(MEF90DefMechAT2_Type), intent(IN)           :: self
      PetscReal                                        :: alpha

      wAT2 = alpha**2
   end function wAT2

#undef __FUNCT__
#define __FUNCT__ "DwAT2"
!!!
!!!
!!!  DwAT2: the derivative of the "w" function of the standard AT2 model, i.e. w(\alpha) = \alpha^2
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function DwAT2(self, alpha)
      class(MEF90DefMechAT2_Type), intent(IN)           :: self
      PetscReal                                        :: alpha

      DwAT2 = 2.0_kr * alpha
   end function DwAT2

#undef __FUNCT__
#define __FUNCT__ "D2wAT2"
!!!
!!!
!!!  D2wAT2: the second derivative of the "w" function of the standard AT2 model, i.e. w(\alpha) = \alpha^2
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function D2wAT2(self, alpha)
      class(MEF90DefMechAT2_Type), intent(IN)           :: self
      PetscReal                                        :: alpha

      D2wAT2 = 2.0_kr
   end function D2wAT2
end module m_MEF90_DefMechAT2
