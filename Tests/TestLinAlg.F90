program testLinAlg
#include "../MEF90/mef90.inc"
#include "petsc/finclude/petsc.h"
   use m_MEF90
   implicit none

   class(mef90Mat), allocatable :: M
   type(MatS2D) :: M2D
   type(MatS3D) :: M3D
   type(Vect2D) :: v2D 
   type(Vect3D) :: v3D 
   PetscReal,dimension(3) :: M2
   PetscReal,dimension(6) :: M3
   
   PetscReal, dimension(:), allocatable :: Marray
   integer :: i, j
   Type(MEF90_MATS), dimension(3) :: eij2D = [MEF90_MATS(1.0_Kr, 0.0_Kr, 0.0_Kr), MEF90_MATS(0.0_Kr, 1.0_Kr, 0.0_Kr), MEF90_MATS(0.0_Kr, 0.0_Kr, 1.0_Kr)]
   PetscErrorCode :: ierr
    
   PetscCallA(PetscInitialize(ierr))
   PetscCallA(MEF90Initialize(PETSC_COMM_WORLD, ierr))

   ! allocate(Marray(3), source = [1.0_Kr, 1.0_Kr, 0.0_Kr])
   allocate(Marray(3))
   Marray = MEF90MatS2DIdentity
   write(*,*) 'Marray: ', Marray
   M2D = Marray
   write(*,*) 'M2D:    ', M2D
   ! M = MatS2D(1.0_Kr,1.0_Kr,0.0_Kr)
   M = MatS3D()

   V2D%X = 1.0_Kr
   V2D%Y = 1.0_Kr

   V3D%X = 1.0_Kr
   V3D%Y = 1.0_Kr
   V3D%Z = 1.0_Kr

   select type(k => M)
      type is (MatS2D)
         write(*,*) 'Mv2D:   ', k*V2D
      type is (MatS3D)
         write(*,*) 'Mv3D:   ', k*V3D
   end select  

   write(*,*) eij2D
   do i = 1, 3
      write(*,*) merge(1.0_Kr, 0.0_Kr, i==1)
   end do

   write([(merge(1.0_Kr, 0.0_Kr, i==1), i=1,3)])
   
   
   PetscCallA(MEF90Finalize(ierr))
   PetscCallA(PetscFinalize(ierr))
end program testLinAlg