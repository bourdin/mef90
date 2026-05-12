#include "../MEF90/mef90.inc"
module m_MEF90_DefMechSplit_class
#include "petsc/finclude/petsc.h"

use m_MEF90_LinAlg_class
use m_MEF90_Parameters
use m_MEF90_BaseClass
use m_MEF90_HookesLaw_Class
implicit none(type, external)
private
public :: MEF90DefMechSplit
public :: MEF90DefMechSplit_SmoothPositiveSquare
public :: MEF90DefMechSplit_DSmoothPositiveSquare
public :: MEF90DefMechSplit_D2SmoothPositiveSquare
!!!
!!!
!!!  MEF90DefMechSplit: The abstract class used to define an energy split for
!!!                           handling unilateral contact
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!      2026 Blaise Bourdin bourdin@mcmaster.ca
!!!

type, extends(MEF90Object), abstract :: MEF90DefMechSplit
   character(len=MEF90MXSTRLEN)                          :: type = ''
   ! integer                                               :: damageOrder = 0
   integer                                               :: quadratureOrder = 0
   PetscBool                                             :: isHybrid = PETSC_FALSE
   class(mef90Mat), allocatable                          :: strain
contains
   procedure(setupInterface), pass(self), deferred       :: setup
   procedure(EEDInterface), pass(self), deferred         :: EED
   procedure(DEEDInterface), pass(self), deferred        :: DEED
   procedure(D2EEDInterface), pass(self), deferred       :: D2EED
end type MEF90DefMechSplit

abstract interface

   subroutine setupInterface(self, Strain, ierr)
      use m_MEF90
      import :: MEF90DefMechSplit
      implicit none(type, external)

      class(MEF90DefMechSplit), intent(inout) :: self
      class(mef90Mat), intent(IN)             :: Strain
      PetscErrorCode, intent(inout)           :: ierr
   end subroutine setupInterface
   
   subroutine EEDInterface(self, HookesLaw, phi, EEDPlus, EEDMinus, ierr)
      use m_MEF90
      import :: MEF90DefMechSplit
      implicit none(type, external)

      class(MEF90DefMechSplit), intent(IN) :: self
      class(MEF90HookesLaw), intent(IN)    :: HookesLaw
      class(mef90Mat), intent(IN)          :: phi
      PetscReal, intent(OUT)               :: EEDPlus, EEDMinus
      PetscErrorCode, intent(inout)        :: ierr
   end subroutine EEDInterface

   subroutine DEEDInterface(self, HookesLaw, phi, DEEDPlus, DEEDMinus, ierr)
      use m_MEF90
      import :: MEF90DefMechSplit
      implicit none(type, external)

      class(MEF90DefMechSplit), intent(IN) :: self
      class(MEF90HookesLaw), intent(IN)    :: HookesLaw
      class(mef90Mat), intent(IN)          :: phi
      PetscReal, intent(OUT)               :: DEEDPlus, DEEDMinus
      PetscErrorCode, intent(inout)        :: ierr
   end subroutine DEEDInterface

   subroutine D2EEDInterface(self, HookesLaw, phi, psi, D2EEDPlus, D2EEDMinus, ierr)
      use m_MEF90
      import :: MEF90DefMechSplit
      implicit none(type, external)

      class(MEF90DefMechSplit), intent(IN) :: self
      class(MEF90HookesLaw), intent(IN)    :: HookesLaw
      class(mef90Mat), intent(IN)          :: phi, psi
      PetscReal, intent(OUT)               :: D2EEDPlus, D2EEDMinus
      PetscErrorCode, intent(inout)        :: ierr
   end subroutine D2EEDInterface
end interface

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSplit_SmoothPositiveSquare"
!!!
!!!
!!!  MEF90DefMechSplit_SmoothPositiveSquare: a low order polynomial C^2 regularization of (max(x,0))**2, defined by
!!!       0                         if x \le -gamma/2
!!!       (x+\gamma/2)^3/3/\gamma   if -\gamma/2 < x \le \gamma/2
!!!       x^2+gamma^2/12            otherwise
!!!
!!!  (c) 2018-2020 Blaise Bourdin bourdin@lsu.edu
!!!      2026 Blaise Bourdin bourdin@mcmaster.ca
!!!

function MEF90DefMechSplit_SmoothPositiveSquare(x, gamma)
   PetscReal, intent(IN) :: x
   PetscReal, intent(IN) :: gamma
   PetscReal             :: MEF90DefMechSplit_SmoothPositiveSquare

   PetscReal             :: gammaOver2

   gammaOver2 = gamma * 0.5_kr
   if (x <= -gammaOver2) then
      MEF90DefMechSplit_SmoothPositiveSquare = 0.0_kr
   else if (x <= gammaOver2) then
      MEF90DefMechSplit_SmoothPositiveSquare = (x + gammaOver2)**3 / 3.0_kr / gamma
   else
      MEF90DefMechSplit_SmoothPositiveSquare = x**2 + gammaOver2**2 / 3.0_kr
   end if
end function MEF90DefMechSplit_SmoothPositiveSquare

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSplit_DSmoothPositiveSquare"
!!!
!!!
!!!  MEF90DefMechSplit_DSmoothPositiveSquare: the first derivative of MEF90_DefMechSplitHD_PenaltyFunction
!!!
!!!  (c) 2018-2020 Blaise Bourdin bourdin@lsu.edu
!!!      2026 Blaise Bourdin bourdin@mcmaster.ca
!!!

function MEF90DefMechSplit_DSmoothPositiveSquare(x, gamma)
   PetscReal, intent(IN) :: x
   PetscReal, intent(IN) :: gamma
   PetscReal             :: MEF90DefMechSplit_DSmoothPositiveSquare

   PetscReal             :: gammaOver2

   gammaOver2 = gamma * 0.5_kr
   if (x <= -gammaOver2) then
      MEF90DefMechSplit_DSmoothPositiveSquare = 0.0_kr
   else if (x <= gammaOver2) then
      MEF90DefMechSplit_DSmoothPositiveSquare = (x + gammaOver2)**2 / gamma
   else
      MEF90DefMechSplit_DSmoothPositiveSquare = 2.0_kr * x
   end if
end function MEF90DefMechSplit_DSmoothPositiveSquare

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSplit_D2SmoothPositiveSquare"
!!!
!!!
!!!  MEF90DefMechSplit_D2SmoothPositiveSquare: the second derivative of MEF90_DefMechSplitHD_PenaltyFunction
!!!
!!!  (c) 2018-2020 Blaise Bourdin bourdin@lsu.edu
!!!      2026 Blaise Bourdin bourdin@mcmaster.ca
!!!

function MEF90DefMechSplit_D2SmoothPositiveSquare(x, gamma)
   PetscReal, intent(IN) :: x
   PetscReal, intent(IN) :: gamma
   PetscReal             :: MEF90DefMechSplit_D2SmoothPositiveSquare

   PetscReal             :: gammaOver2

   gammaOver2 = gamma * 0.5_kr
   if (x <= -gammaOver2) then
      MEF90DefMechSplit_D2SmoothPositiveSquare = 0.0_kr
   else if (x <= gammaOver2) then
      MEF90DefMechSplit_D2SmoothPositiveSquare = 1.0_kr + x / gammaOver2
   else
      MEF90DefMechSplit_D2SmoothPositiveSquare = 2.0_kr
   end if
end function MEF90DefMechSplit_D2SmoothPositiveSquare

end module m_MEF90_DefMechSplit_class

