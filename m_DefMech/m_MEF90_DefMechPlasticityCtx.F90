#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
Module MEF90_APPEND(m_MEF90_DefMechPlasticityCtx,MEF90_DIM)D
#include "finclude/petscdef.h"
   use m_MEF90

type :: MEF90DefMechPlasticityCtx
   Type(MEF90_HOOKESLAW)       :: HookesLaw
   real(Kind = Kr)             :: YieldStress
   real(Kind = Kr)             :: DuctileCouplingPower
   Type(MEF90_MATS)            :: InelasticStrain
   Type(MEF90_MATS)            :: PlasticStrainOld
   Type(MEF90_MATS)            :: plasticStrainPrevious
   real(Kind = Kr)             :: Damage
   real(Kind = Kr)             :: residualStiffness
   real(Kind = Kr)             :: residualYieldStress
   real(Kind = Kr)             :: CoefficientLinSoft
   real(Kind = Kr)             :: CoefficientDruckerPrager
   real(Kind = Kr)             :: CoefficientCapModel0
   real(Kind = Kr)             :: CoefficientCapModel1
   real(Kind = Kr)             :: CoefficientCapModel2
   real(Kind = Kr)             :: CoefficientCapModelD
   logical(Kind = Kr)          :: isPlaneStress
   real(Kind = Kr)             :: cumulatedDissipatedPlasticEnergy
   logical(Kind = Kr)          :: isLinearIsotropicHardening
   real(Kind = Kr)             :: CoeffF
   real(Kind = Kr)             :: CoeffG
   real(Kind = Kr)             :: CoeffH
   real(Kind = Kr)             :: CoeffM
   real(Kind = Kr)             :: CoeffN
   real(Kind = Kr)             :: CoeffL
   real(Kind = Kr)             :: YieldTau0
   real(Kind = Kr)             :: residualYieldTau0
   real(Kind = Kr)             :: phi1
   real(Kind = Kr)             :: phi2
   real(Kind = Kr)             :: Phi
   real(Kind = Kr)             :: delta
   logical(Kind = Kr)          :: isNoPlCoupling
   Type(MEF90RotationMatrix3D) :: RotationMatrix3D
   logical(Kind = Kr)          :: isViscousPlasticity
   real(Kind = Kr)             :: ViscosityGamma0
   real(Kind = Kr)             :: ViscosityN
   real(Kind = Kr)             :: Viscositydt
   real(Kind = Kr)             :: viscousCumulatedDissipatedPlasticEnergyVariation
   real(Kind = Kr)             :: m
end type MEF90DefMechPlasticityCtx

End Module MEF90_APPEND(m_MEF90_DefMechPlasticityCtx,MEF90_DIM)D
