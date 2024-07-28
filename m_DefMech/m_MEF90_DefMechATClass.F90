#include "../MEF90/mef90.inc"
Module m_MEF90_DefMechAT_class
#include "petsc/finclude/petsc.h"

   Use m_MEF90
   Implicit none

!!!
!!!  
!!!  MEF90_DefMechAT_Type: The abstract class used to define a generalized Ambrosio-Tortorelli phase field model
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!  (c) 2024 Bliase Bourdin bourdin@mcmaster.ca
!!!

   Type, abstract :: MEF90DefMechAT_Type
      PetscReal                                        :: cw
      PetscReal                                        :: internalLength
      PetscInt                                         :: aOrder
      PetscInt                                         :: wOrder  
      Type(MPI_comm)                                   :: comm
      Character(len=MEF90MXSTRLEN)                     :: prefix
      Character(len=MEF90MXSTRLEN)                     :: type
   Contains
      Procedure(ATInterface), pass(self), deferred     :: a
      Procedure(ATInterface), pass(self), deferred     :: Da
      Procedure(ATInterface), pass(self), deferred     :: D2a
      Procedure(ATInterface), pass(self), deferred     :: w
      Procedure(ATInterface), pass(self), deferred     :: Dw
      Procedure(ATInterface), pass(self), deferred     :: D2w
      ! Procedure(UnaryInterface), pass(self), deferred  :: view
   End Type

   Abstract Interface
      PetscReal function ATInterface(self,alpha)
         import :: MEF90DefMechAT_Type
         Class(MEF90DefMechAT_Type),Intent(IN)         :: self
         PetscReal,Intent(IN)                          :: alpha
      End function ATInterface
   End Interface

   ! Abstract Interface
   !    PetscReal function UnaryInterface(self)
   !       import :: MEF90DefMechAT_Type
   !       Class(MEF90DefMechAT_Type),Intent(IN)         :: self
   !    End function UnaryInterface
   ! End Interface
End Module m_MEF90_DefMechAT_class

