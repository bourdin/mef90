#include "../MEF90/mef90.inc"
module m_MEF90_DefMechDrivingForce_class
#include "petsc/finclude/petsc.h"

   use m_MEF90_Parameters
   use petscsys
   implicit none(type, external)

!!!
!!!
!!!  MEF90_DefMechDrivingForce_Type: The abstract class used to define a generalized Ambrosio-Tortorelli phase field model
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!

   type, abstract :: MEF90_DefMechDrivingForce_Type
      character(len=MEF90MXSTRLEN)                          :: type
   contains
      procedure(DrivingForceInterface), pass(self), deferred :: df
      procedure(DrivingForceInterface), pass(self), deferred :: Ddf
   end type

   abstract interface
      PetscReal function DrivingForceInterface(self, alpha)
         import :: MEF90_DefMechDrivingForce_Type
         class(MEF90_DefMechDrivingForce_Type), intent(IN)    :: self
         PetscReal                                           :: alpha
      end function DrivingForceInterface
   end interface
end module m_MEF90_DefMechDrivingForce_class

