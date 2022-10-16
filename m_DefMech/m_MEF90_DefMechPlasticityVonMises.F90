#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
Module MEF90_APPEND(m_MEF90_DefMechPlasticityVonMises,MEF90_DIM)D
#include "petsc/finclude/petsc.h"

use m_MEF90
   use m_MEF90_DefMechCtx
   use MEF90_APPEND(m_MEF90_DefMechPlasticityCtx,MEF90_DIM)D
   implicit NONE

Contains
#undef __FUNCT__
#define __FUNCT__ "FHG_VONMISES"
!!!
!!!
!!!  fhg: VonMises
!!!
!!!  (c) 2015 Erwan Tanne : erwan.tanne@gmail.com
!!!
!!!

   subroutine FHG_VONMISES(x,f,h,g,myctx) bind(c)
      use,intrinsic :: iso_c_binding
      use m_MEF90

      real(kind=c_double)                       :: x(*)
      real(kind=c_double)                       :: f(*)
      real(kind=c_double)                       :: h(*)
      real(kind=c_double)                       :: g(*)
      real(Kind = Kr)                           :: StiffnessA,PlasticStrainCumulated
      real(Kind = Kr)                           :: StiffnessB,Stiffness
      type(c_ptr),intent(in),value              :: myctx

      type(MEF90DefMechPlasticityCtx),pointer   :: myctx_ptr
      type(MEF90_MATS)                          :: xMatS
      type(MEF90_MATS)                          :: Stress


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

      Stress=myctx_ptr%HookesLaw*(myctx_ptr%InelasticStrain-xMatS)

      If (myctx_ptr%isLinearIsotropicHardening .eqv. .true. ) then
         PlasticStrainCumulated = (myctx_ptr%cumulatedDissipatedPlasticEnergy+(Stress .DotP. (xMatS-myctx_ptr%PlasticStrainOld)))/myctx_ptr%YieldStress
      else
         PlasticStrainCumulated=0.0_Kr
      endif


      Stiffness=StiffnessB*(1.0_Kr+ PlasticStrainCumulated )

      f(1) = ( (myctx_ptr%HookesLaw *(xMatS-myctx_ptr%PlasticStrainOld)) .DotP. (xMatS-myctx_ptr%PlasticStrainOld) )
      g(1) = StiffnessA * sqrt( MEF90_DIM / (MEF90_DIM - 1.0_kr)  * ( deviatoricPart(Stress)  .DotP.  deviatoricPart(Stress) ))  - ( (1.0_Kr-myctx_ptr%residualYieldStress)*Stiffness + myctx_ptr%residualYieldStress )*myctx_ptr%YieldStress
      g(2) = -PlasticStrainCumulated
      h(1) = Trace(xMatS)
   end subroutine FHG_VONMISES



#undef __FUNCT__
#define __FUNCT__ "FHG_VONMISESPLANETHEORY"
!!!
!!!
!!!  fhg: VonMises
!!!
!!!  (c) 2015 Erwan Tanne : erwan.tanne@gmail.com
!!!
!!!

   subroutine FHG_VONMISESPLANETHEORY(x,f,h,g,myctx) bind(c)
      use,intrinsic :: iso_c_binding
      use m_MEF90

      real(kind=c_double)                       :: x(*)
      real(kind=c_double)                       :: f(*)
      real(kind=c_double)                       :: h(*)
      real(kind=c_double)                       :: g(*)
      real(Kind = Kr)                           :: StiffnessA         !! Stiffness = a(alpha)/b(alpha)
      real(Kind = Kr)                           :: StiffnessB
      type(c_ptr),intent(in),value              :: myctx

      type(MEF90DefMechPlasticityCtx),pointer   :: myctx_ptr
      type(MEF90_MATS)                          :: xMatS
      real(Kind = Kr)                           :: lambda,mu,E,nu
      type(MatS3D)                              :: Strain,PlasticStrainFlow,Stress


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

      E     = myctx_ptr%HookesLaw%YoungsModulus
      nu    = myctx_ptr%HookesLaw%PoissonRatio
      mu    = E / (1.0_Kr + nu) * .5_Kr
      lambda  = E * nu / (1.0_Kr + nu) / (1 - 2.0_Kr * nu)

      if ( myctx_ptr%isPlaneStress .eqv. .true.) then
         !!! If plane stress
         Strain      = 0.0_Kr
         Strain%XX   = myctx_ptr%InelasticStrain%XX
         Strain%YY   = myctx_ptr%InelasticStrain%YY
         Strain%XY   = myctx_ptr%InelasticStrain%XY
         Strain%ZZ   = (-lambda*Trace(myctx_ptr%InelasticStrain) - 2*mu*Trace(myctx_ptr%plasticStrainPrevious))/(lambda + 2*mu)

         PlasticStrainFlow    = 0.0_Kr
         PlasticStrainFlow%XX = xMatS%XX-myctx_ptr%PlasticStrainOld%XX
         PlasticStrainFlow%YY = xMatS%YY-myctx_ptr%PlasticStrainOld%YY
         PlasticStrainFlow%XY = xMatS%XY-myctx_ptr%PlasticStrainOld%XY
         PlasticStrainFlow%ZZ = -( PlasticStrainFlow%XX + PlasticStrainFlow%YY )

         Stress    = 0.0_Kr
         Stress%XX = lambda*(Trace(Strain)) + 2*mu*(Strain%XX-xMatS%XX)
         Stress%YY = lambda*(Trace(Strain)) + 2*mu*(Strain%YY-xMatS%YY)
         Stress%XY = 2*mu*(Strain%XY-xMatS%XY)
         Stress%ZZ = lambda*(Trace(Strain)) + 2*mu*(Strain%ZZ + Trace(xMatS) )
      else
         !!! If plane strain
         Strain      = 0.0_Kr
         Strain%XX   = myctx_ptr%InelasticStrain%XX
         Strain%YY   = myctx_ptr%InelasticStrain%YY
         Strain%XY   = myctx_ptr%InelasticStrain%XY

         PlasticStrainFlow    = 0.0_Kr
         PlasticStrainFlow%XX = xMatS%XX-myctx_ptr%PlasticStrainOld%XX
         PlasticStrainFlow%YY = xMatS%YY-myctx_ptr%PlasticStrainOld%YY
         PlasticStrainFlow%XY = xMatS%XY-myctx_ptr%PlasticStrainOld%XY
         PlasticStrainFlow%ZZ = -( PlasticStrainFlow%XX + PlasticStrainFlow%YY )

         Stress    = 0.0_Kr
         Stress%XX = lambda*(Trace(Strain)) + 2*mu*(Strain%XX-xMatS%XX)
         Stress%YY = lambda*(Trace(Strain)) + 2*mu*(Strain%YY-xMatS%YY)
         Stress%XY = 2*mu*(Strain%XY-xMatS%XY)
         Stress%ZZ = lambda*(Trace(Strain)) + 2*mu*(Strain%ZZ+ Trace(xMatS))
      endif

      if ( myctx_ptr%isNoPlCoupling .eqv. .true.) then
         f(1) = ( PlasticStrainFlow .DotP. PlasticStrainFlow )
         g(1) = StiffnessA * sqrt( (3.0/2.0)*( deviatoricPart(Stress) .dotP. deviatoricPart(Stress) ) ) - myctx_ptr%YieldStress
      else
         f(1) = ( PlasticStrainFlow .DotP. PlasticStrainFlow )
         g(1) = StiffnessA * sqrt( (3.0/2.0)*( deviatoricPart(Stress) .dotP. deviatoricPart(Stress) ) ) - ( (1.0_Kr-myctx_ptr%residualYieldStress)*StiffnessB + myctx_ptr%residualYieldStress )*myctx_ptr%YieldStress
      endif

   end subroutine FHG_VONMISESPLANETHEORY



#undef __FUNCT__
#define __FUNCT__ "FHG_VONMISES1D"
!!!
!!!
!!!  fhg: VonMises
!!!
!!!  (c) 2015 Erwan Tanne : erwan.tanne@gmail.com
!!!
!!!

   subroutine FHG_VONMISES1D(x,f,h,g,myctx) bind(c)
      use,intrinsic :: iso_c_binding
      use m_MEF90

      real(kind=c_double)                       :: x(*)
      real(kind=c_double)                       :: f(*)
      real(kind=c_double)                       :: h(*)
      real(kind=c_double)                       :: g(*)
      real(Kind = Kr)                           :: StiffnessA
      real(Kind = Kr)                           :: StiffnessB
      real(Kind = Kr)                           :: Stiffness,PlasticStrainCumulated
      type(c_ptr),intent(in),value              :: myctx

      type(MEF90DefMechPlasticityCtx),pointer   :: myctx_ptr
      type(MEF90_MATS)                          :: xMatS,Stress

      xMatS = x(1:SIZEOFMEF90_MATS)
      !!! This is the fortran equivalent of casting ctx into a c_ptr
      call c_f_pointer(myctx,myctx_ptr)

       !!! Select which softening young model
      if (myctx_ptr%CoefficientLinSoft==0) then
         StiffnessA = (1.0_Kr - myctx_ptr%Damage)**2
         StiffnessB = (1.0_Kr - myctx_ptr%Damage)**myctx_ptr%DuctileCouplingPower
      else
         StiffnessA = ( (1.0_Kr - myctx_ptr%Damage)**2 /( 1.0_Kr + ( myctx_ptr%CoefficientLinSoft - 1.0_Kr )*(1.0_Kr - (1.0_Kr - myctx_ptr%Damage)**2 ) ) )
         StiffnessB = (1.0_Kr - myctx_ptr%Damage)**myctx_ptr%DuctileCouplingPower
      endif

      Stress = myctx_ptr%HookesLaw*(myctx_ptr%InelasticStrain-xMatS)

      If (myctx_ptr%isLinearIsotropicHardening .eqv. .true. ) then
         PlasticStrainCumulated = (myctx_ptr%cumulatedDissipatedPlasticEnergy+(Stress .DotP. (xMatS-myctx_ptr%PlasticStrainOld)))/myctx_ptr%YieldStress
      else
         PlasticStrainCumulated=0.0_Kr
      endif

      Stiffness=StiffnessB*( 1.0_Kr+ PlasticStrainCumulated )


      f(1) = ( (myctx_ptr%HookesLaw *(xMatS-myctx_ptr%PlasticStrainOld)) .DotP. (xMatS-myctx_ptr%PlasticStrainOld) )
      g(1) = StiffnessA * sqrt( MEF90_DIM / (MEF90_DIM - 1.0_kr)  * ( deviatoricPart(Stress)  .DotP.  deviatoricPart(Stress) ))  - myctx_ptr%YieldStress*Stiffness
      g(2) = -PlasticStrainCumulated
      h(1) = xMatS%YY
      h(2) = xMatS%XY
   end subroutine FHG_VONMISES1D

End Module MEF90_APPEND(m_MEF90_DefMechPlasticityVonMises,MEF90_DIM)D