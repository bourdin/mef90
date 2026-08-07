#include "../MEF90/mef90.inc"
module m_MEF90_DefMechAT1
#include "petsc/finclude/petsc.h"
   use petscsys
   use m_MEF90_Parameters
   use m_MEF90_LinAlg

   use m_MEF90_DefMechAT_class
   implicit none(type)
   private
   public :: MEF90DefMechAT1_Type

   type, extends(MEF90DefMechAT_Type) :: MEF90DefMechAT1_Type

   contains
      procedure, pass(self) :: a => aAT1
      procedure, pass(self) :: Da => DaAT1
      procedure, pass(self) :: D2a => D2aAT1

      procedure, pass(self) :: w => wAT1
      procedure, pass(self) :: Dw => DwAT1
      procedure, pass(self) :: D2w => D2wAT1

      procedure, pass(self) :: setFromOptions => MEF90DefMechAT1_setFromOptions
   end type MEF90DefMechAT1_Type

contains

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechAT1_setFromOptions"
!!!
!!!
!!!  MEF90DefMechAT1_setFromOptions: initializes a MEF90_DefMechAT1_Type from options
!!!  (c) 2025 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine MEF90DefMechAT1_setFromOptions(self, ierr)
      class(MEF90DefMechAT1_Type), intent(inout) :: self
      PetscErrorCode, intent(inout) :: ierr

      self%cw = 2.0_kr / 3.0_kr
      self%aorder = 2
      self%worder = 1
      self%type = 'MEF90DefMechAT1'

      PetscCall(MEF90DefMechAT_setFromOptions(self, ierr))
   end subroutine MEF90DefMechAT1_setFromOptions

#undef __FUNCT__
#define __FUNCT__ "aAT1"
!!! author: Blaise Bourdin (bourdin@lsu.edu, bourdin@mcmaster.ca)
!!! date: 2020
!!! author: somebody else (bourdin@lsu.edu, bourdin@mcmaster.ca)
!!! date: 2026
!!!
!!!  aAT1: the "a" function of the standard AT1 model, i.e. $a(\alpha) = (1-\alpha)^2$
!!!
   PetscReal function aAT1(self, alpha)
      class(MEF90DefMechAT1_Type), intent(IN) :: self
      PetscReal :: alpha

      aAT1 = self%residualStiffness + (1.0_Kr - self%residualStiffness) * (1.0_kr - alpha)**2
   end function aAT1

#undef __FUNCT__
#define __FUNCT__ "DaAT1"
!!!
!!!
!!!  DaAT1: the derivative of the "a" function of the standard AT1 model, i.e. a(\alpha) = (1-\alpha)^2
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!      2026 Blaise Bourdin bourdin@mcmaster.ca
!!!
   PetscReal function DaAT1(self, alpha)
      class(MEF90DefMechAT1_Type), intent(IN) :: self
      PetscReal :: alpha

      DaAT1 = -2.0_kr * (1.0_Kr - self%residualStiffness) * (1.0_kr - alpha)
   end function DaAT1

#undef __FUNCT__
#define __FUNCT__ "D2aAT1"
!!!
!!!
!!!  D2aAT1: the second derivative of the "a" function of the standard AT1 model, i.e. a(\alpha) = (1-\alpha)^2
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!      2026 Blaise Bourdin bourdin@mcmaster.ca
!!!
   PetscReal function D2aAT1(self, alpha)
      class(MEF90DefMechAT1_Type), intent(IN) :: self
      PetscReal :: alpha

      D2aAT1 = 2.0_kr * (1.0_Kr - self%residualStiffness)
   end function D2aAT1

#undef __FUNCT__
#define __FUNCT__ "wAT1"
!!!
!!!
!!!  wAT1: the "w" function of the standard AT1 model, i.e. w(\alpha) = \alpha
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!      2026 Blaise Bourdin bourdin@mcmaster.ca
!!!
   PetscReal function wAT1(self, alpha)
      class(MEF90DefMechAT1_Type), intent(IN) :: self
      PetscReal :: alpha

      wAT1 = alpha
   end function wAT1

#undef __FUNCT__
#define __FUNCT__ "DwAT1"
!!!
!!!
!!!  DwAT1: the derivative of the "w" function of the standard AT1 model, i.e. w(\alpha) = \alpha
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!      2026 Blaise Bourdin bourdin@mcmaster.ca
!!!
   PetscReal function DwAT1(self, alpha)
      class(MEF90DefMechAT1_Type), intent(IN) :: self
      PetscReal :: alpha

      DwAT1 = 1.0_kr
   end function DwAT1

#undef __FUNCT__
#define __FUNCT__ "D2wAT1"
!!!
!!!
!!!  D2wAT1: the second derivative of the "w" function of the standard AT1 model, i.e. w(\alpha) = \alpha
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!      2026 Blaise Bourdin bourdin@mcmaster.ca
!!!
   PetscReal function D2wAT1(self, alpha)
      class(MEF90DefMechAT1_Type), intent(IN) :: self
      PetscReal :: alpha

      D2wAT1 = 0.0_kr
   end function D2wAT1
end module m_MEF90_DefMechAT1
