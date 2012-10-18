Program HeatXfer
#include <finclude/petscdef.h>
   Use m_MEF90_HeatXfer
   Use petsc
   Implicit NONE   

   PetscErrorCode    :: ierr
   Type(Vect2D)      :: V2
   Type(Vect2D)      :: V3
   
   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   Call MEF90_Initialize(ierr)
   V2 = -1.0_Kr
   V3 = 1.2345_Kr
   
   Write(*,*) "V2 in main ", V2
   Call m_MEF90_HeatXferOperatorAssembly(V2,ierr)
   Write(*,*) "V3 in main ", V3
   Call m_MEF90_HeatXferOperatorAssembly(V3,ierr)

   Call MEF90_Finalize(ierr)
   Call PetscFinalize()
End Program HeatXfer