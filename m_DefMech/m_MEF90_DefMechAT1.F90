#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
module m_MEF90_DefMechAT1
#include "finclude/petscdef.h"
   Use m_MEF90
   Use m_MEF90_DefMechAT_class
   implicit none

   Type, extends(MEF90_DefMechAT_Type)                 :: MEF90_DefMechAT1_Type
   Contains
      Procedure, pass(self)                            :: a   => aAT1
      Procedure, pass(self)                            :: Da  => DaAT1
      Procedure, pass(self)                            :: D2a => D2aAT1

      Procedure, pass(self)                            :: w   => wAT1
      Procedure, pass(self)                            :: Dw  => DwAT1
      Procedure, pass(self)                            :: D2w => D2wAT1
   end Type

   interface MEF90_DefMechAT1_Type
      module procedure MEF90_DefMechAT1_Constructor
   end interface

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_DefMechAT1_Constructor"
!!!
!!!  
!!!  MEF90_DefMechAT1_Constructor: the default constructor for a MEF90_DefMechAT1_Type
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Type(MEF90_DefMechAT1_Type) Function MEF90_DefMechAT1_Constructor()
      MEF90_DefMechAT1_Constructor%cw                = 2.0_Kr / 3.0_Kr
      MEF90_DefMechAT1_Constructor%aorder            = 2
      MEF90_DefMechAT1_Constructor%worder            = 1
   End Function MEF90_DefMechAT1_Constructor

#undef __FUNCT__
#define __FUNCT__ "aAT1"
!!!
!!!  
!!!  aAT1: the "a" function of the standard AT1 model, i.e. a(\alpha) = (1-\alpha)^2
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function aAT1(self,alpha)
      Class(MEF90_DefMechAT1_Type),Intent(IN)          :: self
      PetscReal                                        :: alpha

      PetscLogDouble                                   :: flops
      PetscErrorCode                                   :: ierr

      aAT1 = (1.0_kr - alpha)**2
      flops = 2.0
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End function aAT1

#undef __FUNCT__
#define __FUNCT__ "DaAT1"
!!!
!!!  
!!!  DaAT1: the derivative of the "a" function of the standard AT1 model, i.e. a(\alpha) = (1-\alpha)^2
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function DaAT1(self,alpha)
      Class(MEF90_DefMechAT1_Type),Intent(IN)          :: self
      PetscReal                                        :: alpha

      PetscLogDouble                                   :: flops
      PetscErrorCode                                   :: ierr

      DaAT1 = -2.0_Kr * (1.0_Kr - alpha)
      flops = 2.0
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End function DaAT1

#undef __FUNCT__
#define __FUNCT__ "D2aAT1"
!!!
!!!  
!!!  D2aAT1: the second derivative of the "a" function of the standard AT1 model, i.e. a(\alpha) = (1-\alpha)^2
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function D2aAT1(self,alpha)
      Class(MEF90_DefMechAT1_Type),Intent(IN)          :: self
      PetscReal                                        :: alpha

      D2aAT1 = 2.0_Kr
   End function D2aAT1

#undef __FUNCT__
#define __FUNCT__ "wAT1"
!!!
!!!  
!!!  wAT1: the "w" function of the standard AT1 model, i.e. w(\alpha) = \alpha
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function wAT1(self,alpha)
      Class(MEF90_DefMechAT1_Type),Intent(IN)          :: self
      PetscReal                                        :: alpha

      wAT1 = alpha
   End function wAT1

#undef __FUNCT__
#define __FUNCT__ "DwAT1"
!!!
!!!  
!!!  DwAT1: the derivative of the "w" function of the standard AT1 model, i.e. w(\alpha) = \alpha
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function DwAT1(self,alpha)
      Class(MEF90_DefMechAT1_Type),Intent(IN)          :: self
      PetscReal                                        :: alpha

      DwAT1 = 1.0_Kr
   End function DwAT1

#undef __FUNCT__
#define __FUNCT__ "D2wAT1"
!!!
!!!  
!!!  D2wAT1: the second derivative of the "w" function of the standard AT1 model, i.e. w(\alpha) = \alpha
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function D2wAT1(self,alpha)
      Class(MEF90_DefMechAT1_Type),Intent(IN)          :: self
      PetscReal                                        :: alpha

      D2wAT1 = 0.0_Kr
   End function D2wAT1
End module m_MEF90_DefMechAT1
