#include "../MEF90/mef90.inc"
Module m_MEF90_DefMechDrivingForce_class
#include "petsc/finclude/petsc.h"

   Use m_MEF90
   Implicit none

!!!
!!!  
!!!  MEF90_DefMechDrivingForce_Type: The abstract class used to define a generalized Ambrosio-Tortorelli phase field model
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!

   Type, abstract :: MEF90_DefMechDrivingForce_Type
      Character(len=MEF90MXSTRLEN)                          :: type
   Contains
      Procedure(DrivingForceInterface), pass(self), deferred :: df
      Procedure(DrivingForceInterface), pass(self), deferred :: Ddf
   End Type

   Abstract Interface
      PetscReal function DrivingForceInterface(self,alpha)
         import :: MEF90_DefMechDrivingForce_Type
         Class(MEF90_DefMechDrivingForce_Type),Intent(IN)    :: self
         PetscReal                                           :: alpha
      End function DrivingForceInterface
   End Interface
End Module m_MEF90_DefMechDrivingForce_class

