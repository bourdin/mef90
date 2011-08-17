Module m_MEF_Parameters
#include "finclude/petscdef.h"

   Use petsc
   IMPLICIT NONE

   include "exodusII.inc"  
   
   ! The following ensures that mef90 and PETSC real types are compatible:
   ! thank's to Michael Metcalf in comp.lang.fortran
   PetscReal, Parameter                  :: PReal = 1.0
   Integer, Parameter, Public            :: Kr = Selected_Real_Kind(Precision(PReal))
                                         
   PetscInt, Parameter                   :: PInt = 1
   Integer, Parameter, Public            :: Ki = Selected_Int_Kind(PInt)
   
   Integer, Parameter, Public           :: F_In   = 60
   ! Default unit used to read in files
   
   Integer, Parameter, Public           :: F_Out  = 61
   ! Default unit to write files
   
   PetscInt, Public                     :: MEF90_GaussOrder
   PetscInt, Public                     :: MEF90_MyRank
   PetscInt, Public                     :: MEF90_NumProcs
   PetscInt, Parameter, Public          :: MEF90_MXSTRLEN = 256

   !!!   
   !!! Element names
   !!!
   PetscInt, Parameter, Public                   :: MEF90_P1_Lagrange = 1
   PetscInt, Parameter, Public                   :: MEF90_P2_Lagrange = 2
   
!   PetscInt, Parameter, Public                   :: MEF90_Q1_Lagrange = 3
!   PetscInt, Parameter, Public                   :: MEF90_Q2_Lagrange = 4
End Module m_MEF_Parameters
