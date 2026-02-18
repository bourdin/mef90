#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
module MEF90_APPEND(m_MEF90_DefMechSplit_class,MEF90_DIM)D
#include "petsc/finclude/petsc.h"
#define EEDINTERFACE   MEF90_APPEND(EED,MEF90_DIM)D
#define DEEDINTERFACE  MEF90_APPEND(DEED,MEF90_DIM)D
#define D2EEDINTERFACE MEF90_APPEND(D2EED,MEF90_DIM)D

use m_MEF90_Materials
implicit none(type, external)
private
public :: MEF90_DEFMECHSPLIT
public :: MEF90_DefMechSplit_SmoothPositiveSquare
public :: MEF90_DefMechSplit_DSmoothPositiveSquare
public :: MEF90_DefMechSplit_D2SmoothPositiveSquare
!!!
!!!
!!!  MEF90_DefMechSplit_Type: The abstract class used to define an energy split for
!!!                           handling unilateral contact
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!

type, abstract :: MEF90_DEFMECHSPLIT
   character(len=MEF90MXSTRLEN)                      :: type
   integer                                           :: damageOrder
   integer                                           :: strainOrder
contains

   procedure(EEDINTERFACE), pass(self), deferred     :: EED
   procedure(DEEDINTERFACE), pass(self), deferred    :: DEED
   procedure(D2EEDINTERFACE), pass(self), deferred   :: D2EED
end type MEF90_DEFMECHSPLIT

abstract interface
   subroutine EEDINTERFACE(self, Strain, HookesLaw, EEDPlus, EEDMinus)
      use m_MEF90
      import :: MEF90_DEFMECHSPLIT
      implicit none(type, external)

      class(MEF90_DEFMECHSPLIT), intent(IN)          :: self
      type(MEF90_MATS), intent(IN)                   :: Strain
      ! type(MEF90_HOOKESLAW), intent(IN)              :: HookesLaw
      class(MEF90HookesLaw_Type), allocatable        :: HookesLaw

      PetscReal, intent(OUT)                         :: EEDPlus, EEDMinus
   end subroutine EEDINTERFACE

   subroutine DEEDINTERFACE(self, Strain, HookesLaw, DEEDPlus, DEEDMinus)
      use m_MEF90
      import :: MEF90_DEFMECHSPLIT
      implicit none(type, external)

      class(MEF90_DEFMECHSPLIT), intent(IN)           :: self
      type(MEF90_MATS), intent(IN)                    :: Strain
      ! type(MEF90_HOOKESLAW), intent(IN)               :: HookesLaw
      type(MEF90_MATS), intent(OUT)                   :: DEEDPlus, DEEDMinus
      class(MEF90HookesLaw_Type), allocatable         :: HookesLaw

   end subroutine DEEDINTERFACE

   subroutine D2EEDINTERFACE(self, Strain, HookesLaw, D2EEDPlus, D2EEDMinus)
      use m_MEF90
      import :: MEF90_DEFMECHSPLIT
      implicit none(type, external)

      class(MEF90_DEFMECHSPLIT), intent(IN)           :: self
      type(MEF90_MATS), intent(IN)                    :: Strain
      class(MEF90HookesLaw_Type), allocatable         :: HookesLaw, D2EEDPlus, D2EEDMinus
   end subroutine D2EEDINTERFACE
end interface

contains
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

function MEF90_DefMechSplit_SmoothPositiveSquare(x, gamma)
   PetscReal, intent(IN)                             :: x
   PetscReal, intent(IN)                             :: gamma
   PetscReal                                        :: MEF90_DefMechSplit_SmoothPositiveSquare

   PetscReal                                        :: gammaOver2

   gammaOver2 = gamma * 0.5_kr
   if (x <= -gammaOver2) then
      MEF90_DefMechSplit_SmoothPositiveSquare = 0.0_kr
   else if (x <= gammaOver2) then
      MEF90_DefMechSplit_SmoothPositiveSquare = (x + gammaOver2)**3 / 3.0_kr / gamma
   else
      MEF90_DefMechSplit_SmoothPositiveSquare = x**2 + gammaOver2**2 / 3.0_kr
   end if
end function MEF90_DefMechSplit_SmoothPositiveSquare

#undef __FUNCT__
#define __FUNCT__ "MEF90_DefMechSplit_DSmoothPositiveSquare"
!!!
!!!
!!!  MEF90_DefMechSplit_DSmoothPositiveSquare: the first derivative of MEF90_DefMechSplitHD_PenaltyFunction
!!!
!!!  (c) 2018-2020 Blaise Bourdin bourdin@lsu.edu
!!!

function MEF90_DefMechSplit_DSmoothPositiveSquare(x, gamma)
   PetscReal, intent(IN)                             :: x
   PetscReal, intent(IN)                             :: gamma
   PetscReal                                        :: MEF90_DefMechSplit_DSmoothPositiveSquare

   PetscReal                                        :: gammaOver2

   gammaOver2 = gamma * 0.5_kr
   if (x <= -gammaOver2) then
      MEF90_DefMechSplit_DSmoothPositiveSquare = 0.0_kr
   else if (x <= gammaOver2) then
      MEF90_DefMechSplit_DSmoothPositiveSquare = (x + gammaOver2)**2 / gamma
   else
      MEF90_DefMechSplit_DSmoothPositiveSquare = 2.0_kr * x
   end if
end function MEF90_DefMechSplit_DSmoothPositiveSquare

#undef __FUNCT__
#define __FUNCT__ "MEF90_DefMechSplit_D2SmoothPositiveSquare"
!!!
!!!
!!!  MEF90_DefMechSplit_D2SmoothPositiveSquare: the second derivative of MEF90_DefMechSplitHD_PenaltyFunction
!!!
!!!  (c) 2018-2020 Blaise Bourdin bourdin@lsu.edu
!!!

function MEF90_DefMechSplit_D2SmoothPositiveSquare(x, gamma)
   PetscReal, intent(IN)                             :: x
   PetscReal, intent(IN)                             :: gamma
   PetscReal                                        :: MEF90_DefMechSplit_D2SmoothPositiveSquare

   PetscReal                                        :: gammaOver2

   gammaOver2 = gamma * 0.5_kr
   if (x <= -gammaOver2) then
      MEF90_DefMechSplit_D2SmoothPositiveSquare = 0.0_kr
   else if (x <= gammaOver2) then
      MEF90_DefMechSplit_D2SmoothPositiveSquare = 1.0_kr + x / gammaOver2
   else
      MEF90_DefMechSplit_D2SmoothPositiveSquare = 2.0_kr
   end if
end function MEF90_DefMechSplit_D2SmoothPositiveSquare

end module MEF90_APPEND(m_MEF90_DefMechSplit_class,MEF90_DIM)D

