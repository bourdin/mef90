#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
Module MEF90_APPEND(m_MEF90_DefMechPlasticityCrystalBCC,MEF90_DIM)D
#include "petsc/finclude/petsc.h"

use m_MEF90
   use m_MEF90_DefMechCtx
   use MEF90_APPEND(m_MEF90_DefMechPlasticityCtx,MEF90_DIM)D
   implicit NONE

Contains
#undef __FUNCT__
#define __FUNCT__ "FHG_CRYSTALBCC"
!!!
!!!
!!!  fhg: BCC Crystal Plasticity (octahedral): equivalent to FCC in small strains
!!!
!!!  (c) 2022 Jean-Michel Scherer, Caltech <scherer@caltech.edu>,
!!!           Blaise Bourdin, LSU, <bourdin@lsu.edu>
!!!
!!!

   subroutine FHG_CRYSTALBCC(x,f,h,g,myctx) bind(c)
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
      type(MEF90_MATS)                          :: PlasticStrainFlow,Stress
      type(MatS3D)                              :: Strain3D,PlasticStrainFlow3D,Stress3D,Stress3DCrystal,TotalPlasticIncrement,TotalPlasticIncrementCrystal
      type(MatS3D)                              :: PlasticStrain3D,PlasticStrainFlow3DCrystal
      type(MatS3D),dimension(12)                :: MatrixMu
      real(Kind = Kr),dimension(12)             :: ResolvedShearStress,PlasticSlipIncrement
      integer(Kind = Kr)                        :: s,active
      type(Vect3D)                              :: n
      type(Vect3D)                              :: m
      real(Kind = Kr)                           :: normS,dt,CRSS
                                                !!! slip systems nomenclature:
                                                !!!                   Bd,B4         Ba,B2         Bc,B5         Db,D4         Dc,D1         Da,D6
                                                !!!                   Ab,A2         Ad,A6         Ac,A3         Cb,C5         Ca,C3         Cd,C1
      real(Kind = Kr),dimension(3,12)           :: n_s = reshape((/-1., 0., 1.,   0.,-1., 1.,  -1., 1., 0.,  -1., 0., 1.,   0., 1., 1.,   1., 1., 0., &
                                                                  & 0.,-1., 1.,   1., 1., 0.,   1., 0., 1.,  -1., 1., 0.,   1., 0., 1.,   0., 1., 1.  /), (/3,12/)) !!! slip planes normals
      real(Kind = Kr),dimension(3,12)           :: m_s = reshape((/ 1., 1., 1.,   1., 1., 1.,   1., 1., 1.,   1.,-1., 1.,   1.,-1., 1.,   1.,-1., 1., &
                                                                  &-1., 1., 1.,  -1., 1., 1.,  -1., 1., 1.,  -1.,-1., 1.,  -1.,-1., 1.,  -1.,-1., 1.  /), (/3,12/)) !!! slip directions
#if MEF90_DIM==2
      real(Kind = Kr)                           :: E,nu,lambda,mu
#endif
      real(Kind = Kr)                           :: taueq
      
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

      if (myctx_ptr%isNoPlCoupling) then
         StiffnessB = 1.0_Kr
      else
         StiffnessB = ( (1.0_Kr-myctx_ptr%residualYieldStress)*StiffnessB + myctx_ptr%residualYieldStress )
      endif
      
      PlasticStrainFlow = xMatS-myctx_ptr%PlasticStrainOld
      Stress = (myctx_ptr%HookesLaw*(myctx_ptr%totalStrain-xMatS))

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
      
      PlasticStrain3D = 0.0_Kr
      PlasticStrain3D%XX = xMatS%XX
      PlasticStrain3D%YY = xMatS%YY
      PlasticStrain3D%XY = xMatS%XY
      PlasticStrain3D%ZZ = -(PlasticStrain3D%XX + PlasticStrain3D%YY)
#elif MEF90_DIM==3
      Stress3D             = Stress
      Strain3D             = myctx_ptr%totalStrain
      PlasticStrainFlow3D  = xMatS - myctx_ptr%PlasticStrainOld
      PlasticStrain3D      = xMatS
#endif
      
      Stress3DCrystal = MEF90MatRaRt(Stress3D,myctx_ptr%RotationMatrix3D%fullTensor)
      Stress3DCrystal = MEF90MatRaRt(myctx_ptr%RotationMatrix3D%fullTensor,myctx_ptr%RotationMatrix3D%fullTensor)
!!! Not sure why the following works fine but the line above doesn't
!Stress3DCrystal = myctx_ptr%RotationMatrix3D%fullTensor*Stress3D*transpose(myctx_ptr%RotationMatrix3D%fullTensor)
      PlasticStrainFlow3DCrystal = MEF90MatRaRt(PlasticStrainFlow3D,myctx_ptr%RotationMatrix3D%fullTensor)

      normS = (2.0_Kr*SQRT(6.0_Kr))  
      Do s=1,12
         m = m_s(:,s) 
         n = n_s(:,s)
         MatrixMu(s) = ((m .TensP. n) + (n .TensP. m)) / normS
      end do
       
      !if ( .NOT. (myctx_ptr%YieldQ==0.0_Kr) ) then
      !   PlasticStrain3DCrystal = MEF90MatRaRt(PlasticStrain3D,myctx_ptr%RotationMatrix3D%fullTensor)         
      !   PlasticSlips(s) = PlasticStrain3DCrystal .DotP. MatrixMu(s)
      !end if

      TotalPlasticIncrementCrystal = 0.0_Kr
      TotalPlasticIncrement        = 0.0_Kr
      myctx_ptr%viscouscumulatedDissipatedPlasticEnergyVariation = 0.0_Kr

      dt = myctx_ptr%Viscositydt
      active = 0
      taueq = 0.0_Kr
      Do s = 1,12
         ResolvedShearStress(s) =  Stress3DCrystal .DotP. MatrixMu(s)
         CRSS = myctx_ptr%YieldTau0
         !if ( .NOT. (myctx_ptr%YieldQ==0.0_Kr) ) then
         !   Do k=1,12
         !      CRSS = CRSS + myctx_ptr%YieldQ * myctx_ptr%InteractionMatrix%him(s,k) * (1.0_Kr - EXP(-myctx_ptr%Yieldb * ABS(PlasticSlips(k)) ))
         !   end do
         !end if
         if (myctx_ptr%isViscousPlasticity) then    
            PlasticSlipIncrement(s) = dt * myctx_ptr%ViscosityGamma0 * SIGN(1.0_Kr, ResolvedShearStress(s)) *&
                                    & MAX( (ABS(StiffnessA*ResolvedShearStress(s)) -  StiffnessB*CRSS) / myctx_ptr%YieldTau0 , 0. )**myctx_ptr%ViscosityN
            TotalPlasticIncrementCrystal = TotalPlasticIncrementCrystal + (PlasticSlipIncrement(s) * MatrixMu(s))
            !myctx_ptr%plasticSlipsVariation(s) = PlasticSlipIncrement(s)
            myctx_ptr%viscouscumulatedDissipatedPlasticEnergyVariation = myctx_ptr%viscouscumulatedDissipatedPlasticEnergyVariation + ResolvedShearStress(s)*PlasticSlipIncrement(s)
         else
            taueq = taueq + ABS( (StiffnessA * ResolvedShearStress(s)) /  CRSS ) ** myctx_ptr%m
            !myctx_ptr%PlasticSlipsVariation(s) = SIGN(1.0_KR, ResolvedShearStress(s)) * (Stress3DCrystal .DotP. PlasticStrainFlow3DCrystal) * (ABS( (StiffnessA * ResolvedShearStress(s)) / CRSS ) ** (myctx_ptr%m - 1.0_Kr)) / CRSS
         endif
         if ((ABS(StiffnessA*ResolvedShearStress(s)) - StiffnessB*CRSS)>0) then
            active = active + 1
         end if
      End Do
      if (myctx_ptr%isViscousPlasticity) then
         f(1) = (0.5_Kr * StiffnessA * Stress .DotP. (myctx_ptr%totalStrain-xMatS)) + (StiffnessB * myctx_ptr%viscouscumulatedDissipatedPlasticEnergyVariation)
         TotalPlasticIncrement = MEF90MatRtaR(TotalPlasticIncrementCrystal,myctx_ptr%RotationMatrix3D%fullTensor)
#if MEF90_DIM==2
         h(1) = PlasticStrainFlow3D%XX - TotalPlasticIncrement%XX
         h(2) = PlasticStrainFlow3D%XY - TotalPlasticIncrement%XY
         h(3) = PlasticStrainFlow3D%YY - TotalPlasticIncrement%YY
#elif MEF90_DIM==3
         h(1) = NORM(PlasticStrainFlow3D - TotalPlasticIncrement)
#endif
      else
         f(1) = ( (myctx_ptr%HookesLaw *(xMatS-myctx_ptr%PlasticStrainOld)) .DotP. (xMatS-myctx_ptr%PlasticStrainOld) )
         g(1) = ((taueq) ** (1.0_Kr/myctx_ptr%m)) - StiffnessB
         h(1) = Trace(xMatS)
      end if
   end subroutine FHG_CRYSTALBCC
End Module MEF90_APPEND(m_MEF90_DefMechPlasticityCrystalBCC,MEF90_DIM)D
