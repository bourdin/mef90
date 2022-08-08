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
      PetscReal                                        :: cw
      PetscInt                                         :: aOrder
      PetscInt                                         :: wOrder  
      Character(len=MEF90MXSTRLEN)                    :: type
   Contains
      Procedure(ATInterface), pass(self), deferred     :: a
      Procedure(ATInterface), pass(self), deferred     :: Da
      Procedure(ATInterface), pass(self), deferred     :: D2a
      Procedure(ATInterface), pass(self), deferred     :: w
      Procedure(ATInterface), pass(self), deferred     :: Dw
      Procedure(ATInterface), pass(self), deferred     :: D2w
   End Type

   Abstract Interface
      PetscReal function ATInterface(self,alpha)
         import :: MEF90_DefMechAT_Type
         Class(MEF90_DefMechAT_Type),Intent(IN)        :: self
         PetscReal                                     :: alpha
      End function ATInterface
   End Interface
End Module m_MEF90_DefMechAT_class

