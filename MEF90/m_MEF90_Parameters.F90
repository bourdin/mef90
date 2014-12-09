Module m_MEF90_Parameters
#include "finclude/petscdef.h"
#include "../mef90version.h"

   Use petsc
   IMPLICIT NONE
#include "exodusII.inc"  
   
   ! The following ensures that mef90 and PETSC real types are compatible:
   ! thanks to Michael Metcalf in comp.lang.fortran
   PetscReal,Parameter                 :: PReal = 1.0
   Integer,Parameter,Public            :: Kr = Selected_Real_Kind(Precision(PReal))
                                         
   PetscInt,Parameter                  :: PInt = 1
   Integer,Parameter,Public            :: Ki = Selected_Int_Kind(PInt)
   
   PetscLogDouble,Parameter            :: flop = 1.0
   Integer,Parameter,Public            :: PFlop = Selected_Real_Kind(Precision(flop))
                                         
   PetscInt,Parameter,Public           :: MEF90_MXSTRLEN = 256
End Module m_MEF90_Parameters
