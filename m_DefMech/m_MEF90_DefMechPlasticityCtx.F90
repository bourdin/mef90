#include "../MEF90/mef90.inc"
module MEF90_APPEND(m_MEF90_DefMechPlasticityCtx,MEF90_DIM)D
#include "petsc/finclude/petsc.h"
use m_MEF90_Materials

type :: MEF90DefMechPlasticityCtx
   type(MEF90_HOOKESLAW)       :: HookesLaw
   real(Kind=Kr)               :: YieldStress
   real(Kind=Kr)               :: DuctileCouplingPower
   type(MEF90_MATS)            :: totalStrain
   type(MEF90_MATS)            :: PlasticStrainOld
   type(MEF90_MATS)            :: plasticStrainPrevious
   real(Kind=Kr)               :: Damage
   real(Kind=Kr)               :: residualStiffness
   real(Kind=Kr)               :: residualYieldStress
   real(Kind=Kr)               :: CoefficientLinSoft
   real(Kind=Kr)               :: CoefficientDruckerPrager
   real(Kind=Kr)               :: CoefficientCapModel0
   real(Kind=Kr)               :: CoefficientCapModel1
   real(Kind=Kr)               :: CoefficientCapModel2
   real(Kind=Kr)               :: CoefficientCapModelD
   PetscBool                   :: isPlaneStress
   real(Kind=Kr)               :: cumulatedDissipatedPlasticEnergy
   PetscBool                   :: isLinearIsotropicHardening
   real(Kind=Kr)               :: CoeffF
   real(Kind=Kr)               :: CoeffG
   real(Kind=Kr)               :: CoeffH
   real(Kind=Kr)               :: CoeffM
   real(Kind=Kr)               :: CoeffN
   real(Kind=Kr)               :: CoeffL
   real(Kind=Kr)               :: YieldTau0
   real(Kind=Kr)               :: residualYieldTau0
   real(Kind=Kr)               :: phi1
   real(Kind=Kr)               :: phi2
   real(Kind=Kr)               :: Phi
   real(Kind=Kr)               :: delta
   PetscBool                   :: isNoPlCoupling
   type(MEF90RotationMatrix3D) :: RotationMatrix3D
   PetscBool                   :: isViscousPlasticity
   real(Kind=Kr)               :: ViscosityGamma0
   real(Kind=Kr)               :: ViscosityN
   real(Kind=Kr)               :: Viscositydt
   real(Kind=Kr)               :: viscousCumulatedDissipatedPlasticEnergyVariation
   real(Kind=Kr)               :: m
end type MEF90DefMechPlasticityCtx
end module MEF90_APPEND(m_MEF90_DefMechPlasticityCtx,MEF90_DIM)D
