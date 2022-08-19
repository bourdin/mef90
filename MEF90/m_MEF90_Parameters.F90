Module m_MEF90_Parameters
#include "petsc/finclude/petsc.h"
#include "../mef90version.h"

   Use petsc
   IMPLICIT NONE
#include "exodusII.inc"  
#include "../mef90version.h"
   
   ! The following ensures that mef90 and PETSC real types are compatible:
   ! thanks to Michael Metcalf in comp.lang.fortran
   PetscReal,Parameter                 :: PReal = 1.0
   Integer,Parameter,Public            :: Kr = Selected_Real_Kind(Precision(PReal))
                                         
   PetscInt,Parameter                  :: PInt = 1
   Integer,Parameter,Public            :: Ki = kind(PInt)
   
   PetscLogDouble,Parameter            :: flop = 1.0
   Integer,Parameter,Public            :: PFlop = Selected_Real_Kind(Precision(flop))
                                         
   PetscInt,Parameter,Public           :: MEF90MXSTRLEN = 256

   PetscReal,Parameter,Public          :: MEF90INFINITY  = huge(1.0_Kr)
   PetscReal,Parameter,Public          :: MEF90NINFINITY = -huge(1.0_Kr)
End Module m_MEF90_Parameters
