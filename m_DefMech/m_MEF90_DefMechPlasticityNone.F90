#include "../MEF90/mef90.inc"
module MEF90_APPEND(m_MEF90_DefMechPlasticityNone,MEF90_DIM)D
#include "petsc/finclude/petsc.h"

implicit none(type, external)

contains
#undef __FUNCT__
#define __FUNCT__ "FHG_NONE"
!!! author: Erwan Tanne (2015, erwan.tanne@gmail.com)
!!!
!!!  fhg: VonMises
!!!

subroutine FHG_NONE(x, f, h, g, myctx) bind(c)
   use, intrinsic :: iso_c_binding

   real(kind=c_double)                       :: x(*)
   real(kind=c_double)                       :: f(*)
   real(kind=c_double)                       :: h(*)
   real(kind=c_double)                       :: g(*)
   type(c_ptr), intent(in), value              :: myctx

   continue
end subroutine FHG_NONE

end module MEF90_APPEND(m_MEF90_DefMechPlasticityNone,MEF90_DIM)D
