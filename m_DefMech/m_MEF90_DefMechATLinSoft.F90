#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
module m_MEF90_DefMechATLinSoft
#include "petsc/finclude/petsc.h"
   Use m_MEF90
   Use m_MEF90_DefMechAT_class
   implicit none


!!! LinSoft, a model with linear softening phase
!!! function:
!!!
!!! a_k(s) = (1-s)**2/(k+(1-k)(1-s)**2)
!!! w(s)   = 1-(1-s)^2
!!! cw = pi/4
   Type, extends(MEF90DefMechAT_Type)                 :: MEF90DefMechATLinSoft_Type
      PetscReal                                        :: k
   Contains
      Procedure, pass(self)                            :: a   => aLinSoft
      Procedure, pass(self)                            :: Da  => DaLinSoft
      Procedure, pass(self)                            :: D2a => D2aLinSoft

      Procedure, pass(self)                            :: w   => wLinSoft
      Procedure, pass(self)                            :: Dw  => DwLinSoft
      Procedure, pass(self)                            :: D2w => D2wLinSoft
   end Type

   interface MEF90DefMechATLinSoft_Type
      module procedure MEF90DefMechATLinSoft_Constructor
   end interface

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechATLinSoft_Constructor"
!!!
!!!  
!!!  MEF90DefMechATLinSoft_Constructor: the default constructor for a MEF90_DefMechATLinSoft_Type
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Type(MEF90DefMechATLinSoft_Type) Function MEF90DefMechATLinSoft_Constructor(k)
      PetscReal,Intent(IN)                                :: k
      
      MEF90DefMechATLinSoft_Constructor%k                 = k
      MEF90DefMechATLinSoft_Constructor%cw                = PETSC_PI/4.0_Kr
      MEF90DefMechATLinSoft_Constructor%aorder            = 2
      MEF90DefMechATLinSoft_Constructor%worder            = 2
      MEF90DefMechATLinSoft_Constructor%type              = 'MEF90DefMechATLinSoft'
   End Function MEF90DefMechATLinSoft_Constructor

#undef __FUNCT__
#define __FUNCT__ "aLinSoft"
!!!
!!!  
!!!  aLinSoft: the "a" function of the LinSoft model
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function aLinSoft(self,alpha)
      Class(MEF90DefMechATLinSoft_Type),Intent(IN)     :: self
      PetscReal                                        :: alpha

      aLinSoft = (1.0_Kr - alpha)**2 / (self%k + (1.0_Kr - self%k) * (1.0_Kr * alpha)**2)
   End function aLinSoft

#undef __FUNCT__
#define __FUNCT__ "DaLinSoft"
!!!
!!!  
!!!  DaLinSoft: the derivative of the "a" function of the LinSoft model
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function DaLinSoft(self,alpha)
      Class(MEF90DefMechATLinSoft_Type),Intent(IN)     :: self
      PetscReal                                        :: alpha

      DaLinSoft = -2.0_Kr*alpha*(1.0_Kr - alpha)**2 * (self%k - 1.0_Kr) / (alpha**2 * (self%k - 1.0_Kr) + 1.0_Kr)**2 &
                + (2.0_Kr*alpha - 2.0_Kr) / (alpha**2 * (self%k - 1.0_Kr) + 1.0_Kr)
   End function DaLinSoft

#undef __FUNCT__
#define __FUNCT__ "D2aLinSoft"
!!!
!!!  
!!!  D2aLinSoft: the second derivative of the "a" function of the LinSoft model
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function D2aLinSoft(self,alpha)
      Class(MEF90DefMechATLinSoft_Type),Intent(IN)     :: self
      PetscReal                                        :: alpha

      D2aLinSoft = 8.0_Kr*alpha**2 * (1.0_Kr - alpha)**2 * (self%k - 1.0_Kr)**2 / (alpha**2 * (self%k - 1.0_Kr) + 1.0_Kr)**3 &
                  - 4.0_Kr*alpha * (2.0_Kr*alpha - 2.0_Kr) * (self%k - 1.0_Kr) / (alpha**2 * (self%k - 1.0_Kr) + 1.0_Kr)**2   &
                  - 2.0_Kr*(1.0_Kr - alpha)**2 * (self%k - 1.0_Kr) / (alpha**2 * (self%k - 1.0_Kr) + 1.0_Kr)**2               &
                  + 2.0_Kr / (alpha**2 * (self%k - 1.0_Kr) + 1.0_Kr)
   End function D2aLinSoft

#undef __FUNCT__
#define __FUNCT__ "wLinSoft"
!!!
!!!  
!!!  wLinSoft: the "w" function of the LinSoft model
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function wLinSoft(self,alpha)
      Class(MEF90DefMechATLinSoft_Type),Intent(IN)     :: self
      PetscReal                                        :: alpha

      wLinSoft = alpha * (2.0_Kr - alpha)

   End function wLinSoft

#undef __FUNCT__
#define __FUNCT__ "DwLinSoft"
!!!
!!!  
!!!  DwLinSoft: the derivative of the "w" function of the LinSoft model
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function DwLinSoft(self,alpha)
      Class(MEF90DefMechATLinSoft_Type),Intent(IN)     :: self
      PetscReal                                        :: alpha

      DwLinSoft = 2.0_Kr * (1.0_Kr - alpha)

   End function DwLinSoft

#undef __FUNCT__
#define __FUNCT__ "D2wLinSoft"
!!!
!!!  
!!!  D2wLinSoft: the second derivative of the "w" function of the LinSoft model
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function D2wLinSoft(self,alpha)
      Class(MEF90DefMechATLinSoft_Type),Intent(IN)     :: self
      PetscReal                                        :: alpha

      D2wLinSoft = -2.0_Kr
   End function D2wLinSoft
End module m_MEF90_DefMechATLinSoft
