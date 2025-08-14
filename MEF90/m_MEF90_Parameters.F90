module m_MEF90_Parameters
#include "petsc/finclude/petsc.h"
#include "../mef90version.h"

   use petscsys
   implicit none(type, external)
#include "exodusII.inc"
#include "../mef90version.h"

   ! The following ensures that mef90 and PETSC real types are compatible:
   ! thanks to Michael Metcalf in comp.lang.fortran
   PetscReal, parameter                  :: PReal = 1.0
   integer, parameter, public            :: Kr = selected_real_kind(precision(PReal))

   PetscInt, parameter                   :: PInt = 1
   integer, parameter, public            :: Ki = kind(PInt)

   PetscLogDouble, parameter             :: flop = 1.0
   integer, parameter, public            :: PFlop = selected_real_kind(precision(flop))

   PetscInt, parameter, public           :: MEF90MXSTRLEN = 256

   PetscReal, parameter, public          :: MEF90INFINITY = huge(1.0_kr)
   PetscReal, parameter, public          :: MEF90NINFINITY = -huge(1.0_kr)
end module m_MEF90_Parameters
