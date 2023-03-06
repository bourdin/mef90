#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
Module MEF90_APPEND(m_MEF90_DefMechPlasticityCtx,MEF90_DIM)D
#include "petsc/finclude/petsc.h"
   use m_MEF90

   Type :: MEF90DefMechPlasticityCtx
      Type(MEF90_HOOKESLAW)       :: HookesLaw
      Real(Kind = Kr)             :: YieldStress
      Real(Kind = Kr)             :: DuctileCouplingPower
      Type(MEF90_MATS)            :: totalStrain
      Type(MEF90_MATS)            :: PlasticStrainOld
      Type(MEF90_MATS)            :: plasticStrainPrevious
      Real(Kind = Kr)             :: Damage
      Real(Kind = Kr)             :: residualStiffness
      Real(Kind = Kr)             :: residualYieldStress
      Real(Kind = Kr)             :: CoefficientLinSoft
      Real(Kind = Kr)             :: CoefficientDruckerPrager
      Real(Kind = Kr)             :: CoefficientCapModel0
      Real(Kind = Kr)             :: CoefficientCapModel1
      Real(Kind = Kr)             :: CoefficientCapModel2
      Real(Kind = Kr)             :: CoefficientCapModelD
      PetscBool                   :: isPlaneStress
      Real(Kind = Kr)             :: cumulatedPlasticDensity
      PetscBool                   :: isLinearIsotropicHardening
      Real(Kind = Kr)             :: CoeffF
      Real(Kind = Kr)             :: CoeffG
      Real(Kind = Kr)             :: CoeffH
      Real(Kind = Kr)             :: CoeffM
      Real(Kind = Kr)             :: CoeffN
      Real(Kind = Kr)             :: CoeffL
      Real(Kind = Kr)             :: YieldTau0
      Real(Kind = Kr)             :: residualYieldTau0
      Real(Kind = Kr)             :: phi1
      Real(Kind = Kr)             :: phi2
      Real(Kind = Kr)             :: Phi
      Real(Kind = Kr)             :: delta
      PetscBool                   :: isNoPlCoupling
      Type(MEF90RotationMatrix3D) :: RotationMatrix3D
      PetscBool                   :: isViscousPlasticity
      Real(Kind = Kr)             :: ViscosityGamma0
      Real(Kind = Kr)             :: ViscosityN
      Real(Kind = Kr)             :: Viscositydt
      Real(Kind = Kr)             :: viscousCumulatedPlasticDissipationIncrement
      Real(Kind = Kr)             :: m
   End Type MEF90DefMechPlasticityCtx
End Module MEF90_APPEND(m_MEF90_DefMechPlasticityCtx,MEF90_DIM)D
