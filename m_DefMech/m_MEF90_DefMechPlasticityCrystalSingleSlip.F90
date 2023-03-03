#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
Module MEF90_APPEND(m_MEF90_DefMechPlasticityCrystalSingleSlip,MEF90_DIM)D
#include "petsc/finclude/petsc.h"

use m_MEF90
   use m_MEF90_DefMechCtx
   use MEF90_APPEND(m_MEF90_DefMechPlasticityCtx,MEF90_DIM)D
   implicit NONE

Contains
#undef __FUNCT__
#define __FUNCT__ "FHG_CRYSTALSINGLESLIP"
!!!
!!!
!!!  fhg: Single slip Crystal Plasticity: m=[100] and n=[010]
!!!
!!!  (c) 2022 Jean-Michel Scherer, Caltech <scherer@caltech.edu>,
!!!           Blaise Bourdin, LSU, <bourdin@lsu.edu>
!!!
!!!

   subroutine FHG_CRYSTALSINGLESLIP(x,f,h,g,myctx) bind(c)
      use,intrinsic :: iso_c_binding
      use m_MEF90

      real(kind=c_double)                       :: x(*)
      real(kind=c_double)                       :: f(*)
      real(kind=c_double)                       :: h(*)
      real(kind=c_double)                       :: g(*)
      real(kind = Kr)                           :: StiffnessA,StiffnessB
      
      type(c_ptr),intent(in),value              :: myctx
      type(MEF90DefMechPlasticityCtx),pointer   :: myctx_ptr
      type(MEF90_MATS)                          :: xMatS
      type(MEF90_MATS)                          :: PlasticStrainFlow, Stress
      type(MatS3D)                              :: Strain3D, PlasticStrainFlow3D, Stress3D, Stress3DCrystal, MatrixMu, TotalPlasticIncrement, TotalPlasticIncrementCrystal
      real(Kind = Kr),dimension(1)              :: ResolvedShearStress,PlasticSlipIncrement 
      integer(Kind = Kr)                        :: s
      type(Vect3D)                              :: n
      type(Vect3D)                              :: m
      real(Kind = Kr)                           :: normS, dt
      real(Kind = Kr),dimension(3,1)            :: n_s = reshape((/ 0., 1., 0./), (/3,1/)) !!! slip planes normals
      real(Kind = Kr),dimension(3,1)            :: m_s = reshape((/ 1., 0., 0./), (/3,1/)) !!! slip directions
#if MEF90_DIM==2
      real(Kind = Kr)                           :: E,nu,lambda,mu
#endif
      
      !!! Casting x into a MEF90_MATS
      xMatS = x(1:SIZEOFMEF90_MATS)
      !!! This is the fortran equivalent of casting ctx into a c_ptr
      call c_f_pointer(myctx,myctx_ptr)

      !!! Select which softening young model
      if (myctx_ptr%CoefficientLinSoft==0) then
         StiffnessA = (1.0_Kr - myctx_ptr%Damage)**2 + myctx_ptr%residualStiffness
         StiffnessB = (1.0_Kr - myctx_ptr%Damage)**myctx_ptr%DuctileCouplingPower + myctx_ptr%residualStiffness
      else
         StiffnessA = ( (1.0_Kr - myctx_ptr%Damage)**2 /( 1.0_Kr + ( myctx_ptr%CoefficientLinSoft - 1.0_Kr )*(1.0_Kr - (1.0_Kr - myctx_ptr%Damage)**2 ) ) ) + myctx_ptr%residualStiffness
         StiffnessB = (1.0_Kr - myctx_ptr%Damage)**myctx_ptr%DuctileCouplingPower + myctx_ptr%residualStiffness
      endif

      if (myctx_ptr%isNoPlCoupling .eqv. .true.) then
         StiffnessB = 1.0_Kr
      else
         StiffnessB = ( (1.0_Kr-myctx_ptr%residualYieldStress)*StiffnessB + myctx_ptr%residualYieldStress )
      endif

      PlasticStrainFlow = xMatS-myctx_ptr%PlasticStrainOld
      Stress = (myctx_ptr%HookesLaw*(myctx_ptr%totalStrain-xMatS)) !*StiffnessA

#if MEF90_DIM==2
      !!! If plane strain
      E      = myctx_ptr%HookesLaw%YoungsModulus
      nu     = myctx_ptr%HookesLaw%PoissonRatio
      mu     = E / (1.0_Kr + nu) * .5_Kr
      lambda = E * nu / (1.0_Kr + nu) / (1 - 2.0_Kr * nu)

      Strain3D      = 0.0_Kr
      Strain3D%XX   = myctx_ptr%totalStrain%XX
      Strain3D%YY   = myctx_ptr%totalStrain%YY
      Strain3D%XY   = myctx_ptr%totalStrain%XY

      PlasticStrainFlow3D    = 0.0_Kr
      PlasticStrainFlow3D%XX = xMatS%XX - myctx_ptr%PlasticStrainOld%XX
      PlasticStrainFlow3D%YY = xMatS%YY - myctx_ptr%PlasticStrainOld%YY
      PlasticStrainFlow3D%XY = xMatS%XY - myctx_ptr%PlasticStrainOld%XY
      PlasticStrainFlow3D%ZZ = -( PlasticStrainFlow3D%XX + PlasticStrainFlow3D%YY )

      Stress3D     = 0.0_Kr
      Stress3D%XX  = Stress%XX
      Stress3D%XY  = Stress%XY
      Stress3D%YY  = Stress%YY
      Stress3D%ZZ  = lambda*(Trace(Strain3D)) + 2*mu*(Strain3D%ZZ+Trace(xMatS))
#elif MEF90_DIM==3
      Stress3D             = Stress
      Strain3D             = myctx_ptr%totalStrain
      PlasticStrainFlow3D  = xMatS - myctx_ptr%PlasticStrainOld
#endif

      Stress3DCrystal = MEF90MatRaRt(Stress3D,myctx_ptr%RotationMatrix3D%fullTensor)

      TotalPlasticIncrementCrystal = 0.0_Kr
      TotalPlasticIncrement        = 0.0_Kr
      myctx_ptr%viscouscumulatedDissipatedPlasticEnergyVariation = 0.0_Kr
      
      normS = (2.0_Kr*SQRT(1.0_Kr))
      
      dt = myctx_ptr%Viscositydt 

      f(1) = 0.5_Kr * StiffnessA * Stress .DotP. (myctx_ptr%totalStrain-xMatS)
      Do s = 1,1
         m = m_s(:,s) 
         n = n_s(:,s) 
         MatrixMu = ((m .TensP. n) + (n .TensP. m)) / normS 
         ResolvedShearStress(s) =  Stress3DCrystal .DotP. MatrixMu
         if (myctx_ptr%isViscousPlasticity) then 
            PlasticSlipIncrement(s) = dt * myctx_ptr%ViscosityGamma0 * SIGN(1.0_Kr, ResolvedShearStress(s)) *&
                                    & MAX( (ABS(StiffnessA*ResolvedShearStress(s)) -  StiffnessB*myctx_ptr%YieldTau0) / myctx_ptr%YieldTau0 , 0. )**myctx_ptr%ViscosityN
            TotalPlasticIncrementCrystal = TotalPlasticIncrementCrystal + (PlasticSlipIncrement(s) * MatrixMu)
            !myctx_ptr%plasticSlipsVariation(s) = PlasticSlipIncrement(s)   
            myctx_ptr%viscouscumulatedDissipatedPlasticEnergyVariation = myctx_ptr%viscouscumulatedDissipatedPlasticEnergyVariation + ResolvedShearStress(s)*PlasticSlipIncrement(s)
         else
            print *, "Rate-independent crystal plasticity is not implemented."
            !PlasticSlipIncrement(s) = dt * myctx_ptr%eta * SIGN(1.0_Kr, ResolvedShearStress(s)) *&
            !                        & MAX( (ABS(StiffnessA*ResolvedShearStress(s)) -  StiffnessB*myctx_ptr%YieldTau0) / myctx_ptr%YieldTau0 , 0. )**myctx_ptr%ViscosityN
            !TotalPlasticIncrementCrystal = TotalPlasticIncrementCrystal + (PlasticSlipIncrement(s) * MatrixMu)
            !myctx_ptr%viscouscumulatedDissipatedPlasticEnergyVariation = myctx_ptr%viscouscumulatedDissipatedPlasticEnergyVariation + ResolvedShearStress(s)*PlasticSlipIncrement(s)
         endif 
      End Do
      f(1) = f(1) + StiffnessB*myctx_ptr%viscouscumulatedDissipatedPlasticEnergyVariation
      TotalPlasticIncrement = MEF90MatRtaR(TotalPlasticIncrementCrystal,myctx_ptr%RotationMatrix3D%fullTensor)
#if MEF90_DIM==2
      h(1) = PlasticStrainFlow3D%XX - TotalPlasticIncrement%XX
      h(2) = PlasticStrainFlow3D%XY - TotalPlasticIncrement%XY
      h(3) = PlasticStrainFlow3D%YY - TotalPlasticIncrement%YY
#elif MEF90_DIM==3
      h(1) = NORM(PlasticStrainFlow3D - TotalPlasticIncrement)
#endif
   end subroutine FHG_CRYSTALSINGLESLIP

End Module MEF90_APPEND(m_MEF90_DefMechPlasticityCrystalSingleSlip,MEF90_DIM)D
