program testLinAlg
#include "../MEF90/mef90.inc"
#include "petsc/finclude/petsc.h"
   use m_MEF90
   implicit none

   class(mef90Mat), allocatable :: M
   type(MatS2D) :: M2D
   type(MatS3D) :: M3D
   PetscReal, dimension(:), allocatable :: Marray
   PetscErrorCode :: ierr
    
   PetscCallA(PetscInitialize(ierr))
   PetscCallA(MEF90Initialize(PETSC_COMM_WORLD, ierr))

   ! allocate(Marray(3), source = [1.0_Kr, 1.0_Kr, 0.0_Kr])
   allocate(Marray(3))
   Marray = MEF90MatS2DIdentity
   write(*,*) 'Marray: ', Marray
   M2D = Marray
   write(*,*) 'M2D:    ', M2D
   M = MatS2D(1.0_Kr,1.0_Kr,0.0_Kr)
   ! M = Marray

   PetscCallA(MEF90Finalize(ierr))
   PetscCallA(PetscFinalize(ierr))
end program testLinAlg