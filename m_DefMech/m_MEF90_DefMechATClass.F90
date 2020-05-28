#include "../MEF90/mef90.inc"
Module m_MEF90_DefMechAT_class
#include "finclude/petscdef.h"

   Use m_MEF90
   Implicit none

!!!
!!!  
!!!  MEF90_DefMechAT_Type: The abstract class used to define a generalized Ambrosio-Tortorelli phase field model
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!

   Type, abstract :: MEF90_DefMechAT_Type
      PetscInt                                          :: aOrder
      PetscInt                                          :: wOrder      
   Contains
      Procedure(ATInterface), pass(self), deferred      :: ATa
      Procedure(ATInterface), pass(self), deferred      :: ATDa
      Procedure(ATInterface), pass(self), deferred      :: ATD2a
      Procedure(ATInterface), pass(self), deferred      :: ATw
      Procedure(ATInterface), pass(self), deferred      :: ATDw
      Procedure(ATInterface), pass(self), deferred      :: ATD2w
   End Type

   Abstract Interface
      PetscReal function ATInterface(self,alpha)
         import :: MEF90_DefMechAT_Type
         Class(MEF90_DefMechAT_Type),Intent(IN)         :: self
         PetscReal                                      :: alpha
      End function ATInterface
   End Interface
End Module m_MEF90_DefMechAT_class

