#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
!!!
!!! KKL (v2) model, which can be rewritten as a regular GD model with 
!!!  a(\alpha) = g(\alpha) and w(\alpha) = 1-g(\alpha)
!!!  where g(\alpha) = 4(1-\alpha)^3 - 3 (1-\alpha)^4
!!! [Karma et al., 2001] Karma, A., Kessler, D. A., and Levine, H. (2001). 
!!! Phase-field model of mode III dynamic fracture. Phys. Rev. Lett., 87(4):045501.

module m_MEF90_DefMechATKKL
#include "petsc/finclude/petsc.h"
   Use m_MEF90
   Use m_MEF90_DefMechAT_class
   implicit none

   Type, extends(MEF90DefMechAT_Type)                 :: MEF90DefMechATKKL_Type
   Contains
      Procedure, pass(self)                            :: a   => aKKL
      Procedure, pass(self)                            :: Da  => DaKKL
      Procedure, pass(self)                            :: D2a => D2aKKL

      Procedure, pass(self)                            :: w   => wKKL
      Procedure, pass(self)                            :: Dw  => DwKKL
      Procedure, pass(self)                            :: D2w => D2wKKL
   end Type

   interface MEF90DefMechATKKL_Type
      module procedure MEF90DefMechATKKL_Constructor
   end interface

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechATKKL_Constructor"
!!!
!!!  MEF90DefMechATKKL_Constructor: the default constructor for a MEF90_DefMechATKKL_Type
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Type(MEF90DefMechATKKL_Type) Function MEF90DefMechATKKL_Constructor()
      MEF90DefMechATKKL_Constructor%cw                = 0.7165753016381484_Kr
      MEF90DefMechATKKL_Constructor%aorder            = 3
      MEF90DefMechATKKL_Constructor%worder            = 1
      MEF90DefMechATKKL_Constructor%type              = 'MEF90DefMechKKL'
   End Function MEF90DefMechATKKL_Constructor

#undef __FUNCT__
#define __FUNCT__ "aKKL"
!!!
!!!  aKKL: the "a" function of the KKL (v2) model, 
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function aKKL(self,alpha)
      Class(MEF90DefMechATKKL_Type),Intent(IN)         :: self
      PetscReal                                        :: alpha

      aKKL = 4.0_Kr * (1.0_Kr - alpha)**3 - 3.0_Kr * (1.0_Kr - alpha)**4
   End function aKKL

#undef __FUNCT__
#define __FUNCT__ "DaKKL"
!!!
!!!  DaKKL: the derivative of the "a" function of the KKL (v2) model
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function DaKKL(self,alpha)
      Class(MEF90DefMechATKKL_Type),Intent(IN)         :: self
      PetscReal                                        :: alpha

      DaKKL = -12.0_Kr * ((1.0_Kr - alpha)**2 - (1.0_Kr - alpha)**3)
   End function DaKKL

#undef __FUNCT__
#define __FUNCT__ "D2aKKL"
!!!
!!!  D2aKKL: the second derivative of the "a" function of the KKL (v2) model
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function D2aKKL(self,alpha)
      Class(MEF90DefMechATKKL_Type),Intent(IN)         :: self
      PetscReal                                        :: alpha

      D2aKKL = 24.0_Kr * (1.0_Kr - alpha) - 36.0_Kr * (1.0_Kr - alpha)**2
   End function D2aKKL

#undef __FUNCT__
#define __FUNCT__ "wKKL"
!!!
!!!  
!!!  wKKL: the "w" function of the of the KKL (v2) model
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function wKKL(self,alpha)
      Class(MEF90DefMechATKKL_Type),Intent(IN)         :: self
      PetscReal                                        :: alpha

      wKKL = 1.0_Kr - 4.0_Kr * (1.0_Kr - alpha)**3 + 3.0_Kr * (1.0_Kr - alpha)**4
   End function wKKL

#undef __FUNCT__
#define __FUNCT__ "DwKKL"
!!!
!!!  
!!!  DwKKL: the derivative of the "w" function of the KKL (v2) model
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function DwKKL(self,alpha)
      Class(MEF90DefMechATKKL_Type),Intent(IN)         :: self
      PetscReal                                        :: alpha

      DwKKL = 12.0_Kr * ((1.0_Kr - alpha)**2 - (1.0_Kr - alpha)**3)
   End function DwKKL

#undef __FUNCT__
#define __FUNCT__ "D2wKKL"
!!!
!!!  
!!!  D2wKKL: the second derivative of the "w" function of the KKL (v2) model
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function D2wKKL(self,alpha)
      Class(MEF90DefMechATKKL_Type),Intent(IN)         :: self
      PetscReal                                        :: alpha

      D2wKKL = 24.0_Kr * (1.0_Kr - alpha) - 36.0_Kr * (1.0_Kr - alpha)**2
   End function D2wKKL
End module m_MEF90_DefMechATKKL
