#include "../MEF90/mef90.inc"
module m_MEF90_DefMechAT_class
#include "petsc/finclude/petsc.h"

   use petscsys
   use m_MEF90_Parameters
   use iso_c_binding
   implicit none(type, external)

!!!
!!!
!!!  MEF90_DefMechAT_Type: The abstract class used to define a generalized Ambrosio-Tortorelli phase field model
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!

   type, abstract :: MEF90DefMechAT_Type
      PetscReal                                        :: cw
      PetscInt                                         :: aOrder
      PetscInt                                         :: wOrder
      character(len=MEF90MXSTRLEN)                     :: type
   contains
      procedure(ATInterface), pass(self), deferred     :: a
      procedure(ATInterface), pass(self), deferred     :: Da
      procedure(ATInterface), pass(self), deferred     :: D2a
      procedure(ATInterface), pass(self), deferred     :: w
      procedure(ATInterface), pass(self), deferred     :: Dw
      procedure(ATInterface), pass(self), deferred     :: D2w
   end type MEF90DefMechAT_Type

   abstract interface
      PetscReal function ATInterface(self, alpha)
         use petscsys
         import :: MEF90DefMechAT_Type
         class(MEF90DefMechAT_Type), intent(IN)         :: self
         PetscReal                                     :: alpha
      end function ATInterface
   end interface
end module m_MEF90_DefMechAT_class

