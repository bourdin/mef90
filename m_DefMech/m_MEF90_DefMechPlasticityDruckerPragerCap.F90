#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
Module MEF90_APPEND(m_MEF90_DefMechPlasticityDruckerPragerCap,MEF90_DIM)D
#include "petsc/finclude/petsc.h"

   use m_MEF90
   use m_MEF90_DefMechCtx
   use MEF90_APPEND(m_MEF90_DefMechPlasticityCtx,MEF90_DIM)D
   implicit NONE

Contains
#undef __FUNCT__
#define __FUNCT__ "FHG_DRUCKERPRAGERCAPMODEL"
!!!
!!!
!!!  fhg: VonMises
!!!
!!!  (c) 2015 Erwan Tanne : erwan.tanne@gmail.com
!!!
!!!

   subroutine FHG_DRUCKERPRAGERCAPMODEL(x,f,h,g,myctx) bind(c)
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

      Stress=(myctx_ptr%HookesLaw*(myctx_ptr%totalStrain-xMatS))*StiffnessA
      f(1) = ( (myctx_ptr%HookesLaw *(xMatS-myctx_ptr%PlasticStrainOld)) .DotP. (xMatS-myctx_ptr%PlasticStrainOld) ) * StiffnessA / 2.0
      g(1) =  sqrt( MEF90_DIM / (MEF90_DIM - 1.0_kr) * ( deviatoricPart(Stress)  .DotP.  deviatoricPart(Stress) ))  - myctx_ptr%YieldStress*StiffnessB -  myctx_ptr%CoefficientDruckerPrager*Trace(Stress)
      g(2) = myctx_ptr%CoefficientCapModelD * sqrt( MEF90_DIM / (MEF90_DIM - 1.0_kr) * ( deviatoricPart(Stress)  .DotP.  deviatoricPart(Stress) ))  - myctx_ptr%CoefficientCapModel0*StiffnessB -  myctx_ptr%CoefficientCapModel1*Trace(Stress) - myctx_ptr%CoefficientCapModel2*Trace(Stress)**2.0
   end subroutine FHG_DRUCKERPRAGERCAPMODEL
End Module MEF90_APPEND(m_MEF90_DefMechPlasticityDruckerPragerCap,MEF90_DIM)D
