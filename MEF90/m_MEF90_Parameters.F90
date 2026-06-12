module m_MEF90_Parameters
#include "petsc/finclude/petsc.h"
#include "../mef90version.h"

   use petscsys
   implicit none(type, external)
! #include "exodusII.inc"
#include "../mef90version.h"

   private
   public :: Ki, Kr, MEF90MXSTRLEN, MEF90INFINITY, MEF90NINFINITY

   ! The following ensures that mef90 and PETSC real types are compatible:
   ! thanks to Michael Metcalf in comp.lang.fortran
   ! PetscReal, parameter          :: PReal = 1.0
   ! integer, parameter            :: Kr = selected_real_kind(precision(PReal))

   ! PetscInt, parameter           :: PInt = 1
   ! integer, parameter            :: Ki = kind(PInt)

   ! PetscLogDouble, parameter     :: flop = 1.0
   ! integer, parameter            :: PFlop = selected_real_kind(precision(flop))

   integer, parameter            :: Kr = PETSC_REAL_KIND
   integer, parameter            :: Ki = PETSC_INT_KIND

   PetscInt, parameter           :: MEF90MXSTRLEN = 256

   PetscReal, parameter          :: MEF90INFINITY = huge(1.0_kr)
   PetscReal, parameter          :: MEF90NINFINITY = -huge(1.0_kr)
end module m_MEF90_Parameters
