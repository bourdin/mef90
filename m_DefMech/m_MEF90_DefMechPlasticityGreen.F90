#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
Module MEF90_APPEND(m_MEF90_DefMechPlasticityGreen,MEF90_DIM)D
#include "petsc/finclude/petsc.h"

use m_MEF90
   use m_MEF90_DefMechCtx
   use MEF90_APPEND(m_MEF90_DefMechPlasticityCtx,MEF90_DIM)D
   implicit NONE

Contains
#undef __FUNCT__
#define __FUNCT__ "FHG_GREEN"

!!!
!!!
!!!  fhg: Green
!!!
!!!  (c) 2017 Stella Brach, Caltech <brach@caltech.edu>,
!!!           Blaise Bourdin, LSU, <bourdin@lsu.edu>
!!!
!!!

   subroutine FHG_GREEN(x,f,h,g,myctx) bind(c)
      use,intrinsic :: iso_c_binding
      use m_MEF90

      real(kind=c_double)                       :: x(*)
      real(kind=c_double)                       :: f(*)
      real(kind=c_double)                       :: h(*)
      real(kind=c_double)                       :: g(*)
      real(Kind = Kr)                           :: StiffnessA

      type(c_ptr),intent(in),value              :: myctx
      type(MEF90DefMechPlasticityCtx),pointer   :: myctx_ptr
      type(MEF90_MATS)                          :: xMatS
      type(MEF90_MATS)                          :: Stress

      !!! Casting x into a MEF90_MATS
      xMatS = x(1:SIZEOFMEF90_MATS)
      !!! This is the fortran equivalent of casting ctx into a c_ptr
      call c_f_pointer(myctx,myctx_ptr)

      !!! Select which softening young model
      if (myctx_ptr%CoefficientLinSoft==0) then
         StiffnessA = (1.0_Kr - myctx_ptr%Damage)**2 + myctx_ptr%residualStiffness
      else
         StiffnessA = ( (1.0_Kr - myctx_ptr%Damage)**2 /( 1.0_Kr + ( myctx_ptr%CoefficientLinSoft - 1.0_Kr )*(1.0_Kr - (1.0_Kr - myctx_ptr%Damage)**2 ) ) ) + myctx_ptr%residualStiffness
      endif

      Stress=(myctx_ptr%HookesLaw*(myctx_ptr%InelasticStrain-xMatS))*StiffnessA
      f(1) = ( (myctx_ptr%HookesLaw *(xMatS-myctx_ptr%PlasticStrainOld)) .DotP. (xMatS-myctx_ptr%PlasticStrainOld) )

      if ( myctx_ptr%Damage == 0.0_Kr) then
         !!! If porosity equals zero, Green reduces to von Mises. A residual damage term is added in order to avoid instability on plastic admissibility
         g(1) = sqrt((1.0_Kr/2.0_Kr)*(deviatoricPart(Stress)  .DotP.  deviatoricPart(Stress))+(myctx_ptr%delta*Trace(Stress))**2)-myctx_ptr%YieldStress
      else
         !!! If porosity is greater than zero, the Green criterion holds and no plastic admissibility has to be considered
         g(1) = sqrt((1.0_Kr/2.0_Kr)*(deviatoricPart(Stress)  .DotP.  deviatoricPart(Stress))+(myctx_ptr%Damage*Trace(Stress))**2)-myctx_ptr%YieldStress
      endif


   end subroutine FHG_GREEN


End Module MEF90_APPEND(m_MEF90_DefMechPlasticityGreen,MEF90_DIM)D
