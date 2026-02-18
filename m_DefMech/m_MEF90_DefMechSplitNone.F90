#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
module MEF90_APPEND(m_MEF90_DefMechSplitNone,MEF90_DIM)D
#include "petsc/finclude/petsc.h"
use MEF90_APPEND(m_MEF90_DefMechSplit_class,MEF90_DIM)D
use m_MEF90_Materials
#define MEF90_DEFMECHSPLITNONE_CONSTRUCTOR MEF90_APPEND(m_MEF90_DefMechSplitNone_Constructor,MEF90_DIM)D
implicit none(type)
private
public :: MEF90_DEFMECHSPLITNONE

type, extends(MEF90_DEFMECHSPLIT)                   :: MEF90_DEFMECHSPLITNONE
contains
   procedure, pass(self)                            :: EED => EEDNone
   procedure, pass(self)                            :: DEED => DEEDNone
   procedure, pass(self)                            :: D2EED => D2EEDNone
end type MEF90_DEFMECHSPLITNONE

interface MEF90_DEFMECHSPLITNONE
   module procedure MEF90_DEFMECHSPLITNONE_CONSTRUCTOR
end interface

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_DEFMECHSPLITNONE_CONSTRUCTOR"
!!!
!!!
!!!  MEF90_DEFMECHSPLITNONE_CONSTRUCTOR: the default constructor for a MEF90_DEFMECHSPLITNONE
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
type(MEF90_DEFMECHSPLITNONE) function MEF90_DEFMECHSPLITNONE_CONSTRUCTOR()
   MEF90_DEFMECHSPLITNONE_CONSTRUCTOR%damageOrder = 0
   MEF90_DEFMECHSPLITNONE_CONSTRUCTOR%strainOrder = 2
   MEF90_DEFMECHSPLITNONE_CONSTRUCTOR%type = 'MEF90DefMech_unilateralContactTypeNone'
end function MEF90_DEFMECHSPLITNONE_CONSTRUCTOR

#undef __FUNCT__
#define __FUNCT__ "EEDNone"
!!!
!!!
!!!  EEDNone: Compute the positive and negative part of the elastic energy density associated with a strain tensor
!!!           without a split, we have EEDPlus  = 1/2 HookesLaw Strain \cdot Strain
!!!                                    EEDMinus = 0
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
subroutine EEDNone(self, Strain, HookesLaw, EEDPlus, EEDMinus)
   class(MEF90_DEFMECHSPLITNONE), intent(IN)           :: self
   type(MEF90_MATS), intent(IN)                        :: Strain
   type(MEF90_HOOKESLAW), intent(IN)                   :: HookesLaw
   PetscReal, intent(OUT)                              :: EEDPlus, EEDMinus

   EEDPlus = ((HookesLaw * Strain) .dotP.Strain) * 0.5_kr
   EEDMinus = 0.0_kr
end subroutine EEDNone

#undef __FUNCT__
#define __FUNCT__ "DEEDNone"
!!!
!!!
!!!  DEEDNone: Compute the derivative of the positive and negative part of the elastic energy density (positive and negative stress)
!!!               evaluated at the strain tensor Strain.
!!!           without a split, we have DEEDPlus  = HookesLaw Strain
!!!                                    DEEDMinus = 0
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
subroutine DEEDNone(self, Strain, HookesLaw, DEEDPlus, DEEDMinus)
   class(MEF90_DEFMECHSPLITNONE), intent(IN)           :: self
   type(MEF90_MATS), intent(IN)                        :: Strain
   type(MEF90_HOOKESLAW), intent(IN)                   :: HookesLaw
   type(MEF90_MATS), intent(OUT)                       :: DEEDPlus, DEEDMinus

   DEEDPlus = HookesLaw * Strain
   DEEDMinus = 0.0_kr
end subroutine DEEDNone

#undef __FUNCT__
#define __FUNCT__ "D2EEDNone"
!!!
!!!
!!!  D2EEDNone: Compute the second derivative of the positive and negative part of the elastic energy density (positive and negative stress)
!!!               evaluated at the strain tensor Strain.
!!!               without a split, D2EEDPlus = HookesLaw, D2EEDMinus = 0
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
subroutine D2EEDNone(self, Strain, HookesLaw, D2EEDPlus, D2EEDMinus)
   class(MEF90_DEFMECHSPLITNONE), intent(IN)           :: self
   type(MEF90_MATS), intent(IN)                        :: Strain
   type(MEF90_HOOKESLAW), intent(IN)                   :: HookesLaw
   type(MEF90_HOOKESLAW), intent(OUT)                  :: D2EEDPlus, D2EEDMinus

   D2EEDPlus = HookesLaw

#if MEF90_DIM==2
   D2EEDMinus%isPlaneStress = HookesLaw%isPlaneStress
#endif
   D2EEDMinus%type = HookesLaw%type
   select case (HookesLaw%type)
   case (MEF90HookesLawTypeIsotropic)
      D2EEDMinus%YoungsModulus = 0.0_kr * HookesLaw%YoungsModulus
      D2EEDMinus%PoissonRatio = 0.0_kr * HookesLaw%PoissonRatio
      D2EEDMinus%lambda = 0.0_kr * HookesLaw%lambda
      D2EEDMinus%mu = 0.0_kr * HookesLaw%mu
      D2EEDMinus%BulkModulus = 0.0_kr * HookesLaw%BulkModulus
   case (MEF90HookesLawTypeFull)
      D2EEDMinus%FullTensor = 0.0_kr * HookesLaw%FullTensor
      D2EEDMinus%fullTensorLocal = 0.0_kr * HookesLaw%fullTensorLocal
   end select !HookesLaw%type
end subroutine D2EEDNone

end module MEF90_APPEND(m_MEF90_DefMechSplitNone,MEF90_DIM)D
