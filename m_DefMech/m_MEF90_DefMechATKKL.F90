#include "../MEF90/mef90.inc"
!!!
!!! KKL (v2) model, which can be rewritten as a regular GD model with
!!!  $a(\alpha) = g(\alpha)$ and $w(\alpha) = 1-g(\alpha)$
!!!  where $g(\alpha) = 4(1-\alpha)^3 - 3 (1-\alpha)^4$
!!! Karma, A., Kessler, D. A., and Levine, H. (2001).
!!! Phase-field model of mode III dynamic fracture. Phys. Rev. Lett., 87(4):045501. doi:10.1103/PhysRevLett.87.045501
!!! 
module m_MEF90_DefMechATKKL
#include "petsc/finclude/petsc.h"
   use petscsys
   use m_MEF90_Parameters
   use m_MEF90_LinAlg

   use m_MEF90_DefMechAT_class
   implicit none(type)
   private
   public :: MEF90DefMechATKKL_Type
   type, extends(MEF90DefMechAT_Type) :: MEF90DefMechATKKL_Type
   contains
      procedure, pass(self) :: a => aKKL
      procedure, pass(self) :: Da => DaKKL
      procedure, pass(self) :: D2a => D2aKKL

      procedure, pass(self) :: w => wKKL
      procedure, pass(self) :: Dw => DwKKL
      procedure, pass(self) :: D2w => D2wKKL

      procedure, pass(self) :: setFromOptions => MEF90DefMechATKKL_setFromOptions
   end type MEF90DefMechATKKL_Type

contains

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechATKKL_setFromOptions"
!!! author: Blaise Bourdin (2025-26, bourdin@mcmaster.ca)
!!!
!!!  MEF90DefMechATKKL_setFromOptions: initializes a MEF90_DefMechATKKL_Type from options
!!!
   subroutine MEF90DefMechATKKL_setFromOptions(self,ierr)
      class(MEF90DefMechATKKL_Type), intent(inout) :: self
      PetscErrorCode,intent(inout) :: ierr

      self%cw = 0.7165753016381484_kr
      self%aorder = 3
      self%worder = 3
      self%type = 'MEF90DefMechATKKL'

      PetscCall(MEF90DefMechAT_setFromOptions(self, ierr))
   end subroutine MEF90DefMechATKKL_setFromOptions

#undef __FUNCT__
#define __FUNCT__ "aKKL"
!!! author: Blaise Bourdin (2020, bourdin@lsu.edu)
!!!
!!!  aKKL: the "a" function of the KKL (v2) model,
!!!
   PetscReal function aKKL(self, alpha)
      class(MEF90DefMechATKKL_Type), intent(IN) :: self
      PetscReal :: alpha

      aKKL = 4.0_kr * (1.0_kr - alpha)**3 - 3.0_kr * (1.0_kr - alpha)**4
   end function aKKL

#undef __FUNCT__
#define __FUNCT__ "DaKKL"
!!! author: Blaise Bourdin (2020, bourdin@lsu.edu)
!!!
!!!  DaKKL: the derivative of the "a" function of the KKL (v2) model
!!!
   PetscReal function DaKKL(self, alpha)
      class(MEF90DefMechATKKL_Type), intent(IN) :: self
      PetscReal :: alpha

      DaKKL = -12.0_kr * ((1.0_kr - alpha)**2 - (1.0_kr - alpha)**3)
   end function DaKKL

#undef __FUNCT__
#define __FUNCT__ "D2aKKL"
!!! author: Blaise Bourdin (2020, bourdin@lsu.edu)
!!!
!!!  D2aKKL: the second derivative of the "a" function of the KKL (v2) model
!!!
   PetscReal function D2aKKL(self, alpha)
      class(MEF90DefMechATKKL_Type), intent(IN) :: self
      PetscReal :: alpha

      D2aKKL = 24.0_kr * (1.0_kr - alpha) - 36.0_kr * (1.0_kr - alpha)**2
   end function D2aKKL

#undef __FUNCT__
#define __FUNCT__ "wKKL"
!!! author: Blaise Bourdin (2020, bourdin@lsu.edu)
!!!
!!!  wKKL: the "w" function of the of the KKL (v2) model
!!!
   PetscReal function wKKL(self, alpha)
      class(MEF90DefMechATKKL_Type), intent(IN) :: self
      PetscReal :: alpha

      wKKL = 1.0_kr - 4.0_kr * (1.0_kr - alpha)**3 + 3.0_kr * (1.0_kr - alpha)**4
   end function wKKL

#undef __FUNCT__
#define __FUNCT__ "DwKKL"
!!! author: Blaise Bourdin (2020, bourdin@lsu.edu)
!!!
!!!  DwKKL: the derivative of the "w" function of the KKL (v2) model
!!!
   PetscReal function DwKKL(self, alpha)
      class(MEF90DefMechATKKL_Type), intent(IN) :: self
      PetscReal :: alpha

      DwKKL = 12.0_kr * ((1.0_kr - alpha)**2 - (1.0_kr - alpha)**3)
   end function DwKKL

#undef __FUNCT__
#define __FUNCT__ "D2wKKL"
!!! author: Blaise Bourdin (2020, bourdin@lsu.edu)
!!!
!!!  D2wKKL: the second derivative of the "w" function of the KKL (v2) model
!!!
   PetscReal function D2wKKL(self, alpha)
      class(MEF90DefMechATKKL_Type), intent(IN) :: self
      PetscReal :: alpha

      D2wKKL = -24.0_kr * (1.0_kr - alpha) + 36.0_kr * (1.0_kr - alpha)**2
   end function D2wKKL
end module m_MEF90_DefMechATKKL
