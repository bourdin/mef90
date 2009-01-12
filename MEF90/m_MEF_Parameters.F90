Module m_MEF_Parameters
#include "finclude/petscdef.h"
   Use petsc
   IMPLICIT NONE
!   PRIVATE

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
   
   Integer, Public                      :: MEF90_GaussOrder
   Integer, Public                      :: MEF90_MyRank
   Integer, Public                      :: MEF90_NumProcs
End Module m_MEF_Parameters
