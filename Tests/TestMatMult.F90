program testMatMult
#include "../MEF90/mef90.inc"
#include "petsc/finclude/petsc.h"
   use m_MEF90
   implicit none

   class(mef90Mat), allocatable  :: M2, M3
   type(Mat2D)                   :: K2
   class(mef90Vect), allocatable :: u2, v2, u3, v3

   PetscErrorCode :: ierr
    
   PetscCallA(PetscInitialize(ierr))
   PetscCallA(MEF90Initialize(PETSC_COMM_WORLD, ierr))

   M2 = Mat2D(1., 2., 3., 4.)
   K2 = Mat2D(1,2,3,4)
   u2 = Vect2D(1,1)
   write(*,*) trace(M2)
   write(*,*) trace(K2)
   
   PetscCallA(MEF90Finalize(ierr))
   PetscCallA(PetscFinalize(ierr))
end program testMatMult
