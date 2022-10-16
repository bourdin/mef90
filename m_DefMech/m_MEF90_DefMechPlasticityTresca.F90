#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
Module MEF90_APPEND(m_MEF90_DefMechPlasticityTresca,MEF90_DIM)D
#include "petsc/finclude/petsc.h"

use m_MEF90
   use m_MEF90_DefMechCtx
   use MEF90_APPEND(m_MEF90_DefMechPlasticityCtx,MEF90_DIM)D
   implicit NONE

Contains
#undef __FUNCT__
#define __FUNCT__ "FHG_TRESCA"
!!!
!!!
!!!  fhg: Tresca
!!!
!!!  (c) 2015 Erwan Tanne : erwan.tanne@gmail.com
!!!
!!!

   subroutine FHG_TRESCA(x,f,h,g,myctx) bind(c)
      use,intrinsic :: iso_c_binding
      use m_MEF90

      real(kind=c_double)                       :: x(*)
      real(kind=c_double)                       :: f(*)
      real(kind=c_double)                       :: h(*)
      real(kind=c_double)                       :: g(*)
      type(c_ptr),intent(in),value              :: myctx

      type(MEF90DefMechPlasticityCtx),pointer   :: myctx_ptr
      type(MEF90_MATS)                          :: xMatS
      type(MEF90_MAT)                           :: MatProjLocBasisToPrincipalBasis
      type(MEF90_MATS)                          :: MatDiagPrincipalBasis

      xMatS = x(1:SIZEOFMEF90_MATS)
      !!! This is the fortran equivalent of casting ctx into a c_ptr
      call c_f_pointer(myctx,myctx_ptr)


      !write(*,*) 'A.e(u):         ', myctx_ptr%HookesLaw*myctx_ptr%Strain
      ! D=P^(-1).A.P
      call Diagonalize(deviatoricPart(myctx_ptr%HookesLaw*myctx_ptr%InelasticStrain),MatProjLocBasisToPrincipalBasis,MatDiagPrincipalBasis)

      f(1) = ( (myctx_ptr%HookesLaw * (xMatS-myctx_ptr%PlasticStrainOld)) .DotP. (xMatS-myctx_ptr%PlasticStrainOld) ) /2.
      h(1) = Trace(xMatS)

#if MEF90_DIM == 2
      g(1) = +(MatDiagPrincipalBasis%XX-MatDiagPrincipalBasis%YY) - myctx_ptr%YieldStress
      g(2) = -(MatDiagPrincipalBasis%XX-MatDiagPrincipalBasis%YY) - myctx_ptr%YieldStress
      g(3) = +(MatDiagPrincipalBasis%YY)                 - myctx_ptr%YieldStress
      g(4) = -(MatDiagPrincipalBasis%YY)                 - myctx_ptr%YieldStress
      g(5) = +(MatDiagPrincipalBasis%XX)                 - myctx_ptr%YieldStress
      g(6) = -(MatDiagPrincipalBasis%XX)                 - myctx_ptr%YieldStress
#else
      Write(*,*) 'Tresca3D is NOT implemented'
#endif
   end subroutine FHG_TRESCA


End Module MEF90_APPEND(m_MEF90_DefMechPlasticityTresca,MEF90_DIM)D
