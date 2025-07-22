program TestSpectralDecomposition
#include "petsc/finclude/petsc.h"
   use m_MEF90
   implicit none(type, external)

   type(MatS2D)                     :: M2D, D2D
   type(Mat2D)                      :: P2D
   type(MatS3D)                     :: M3D, D3D
   type(Mat3D)                      :: P3D
   type(MatS2D), dimension(2)        :: ppleDirections2D
   type(MatS3D), dimension(3)        :: ppleDirections3D
   PetscReal, dimension(2)           :: ppleValues2D
   PetscReal, dimension(3)           :: ppleValues3D
   PetscBool                        :: flg
   PetscErrorCode                   :: ierr
   PetscInt                         :: i, j, n = 10
   PetscRandom                      :: RdmCtx

   PetscCallA(PetscInitialize(ierr))
   PetscCallA(MEF90Initialize(PETSC_COMM_WORLD, ierr))
   PetscCallA(PetscRandomCreate(PETSC_COMM_WORLD, RdmCtx, ierr))
   PetscCallA(PetscRandomSetFromOptions(RdmCtx, ierr))
   PetscCallA(PetscRandomSetInterval(RdmCtx, -1.0_kr, 1.0_kr, ierr))

   PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS, '', '-n', n, flg, ierr))
   write (*, *) 'Testing SpectralDecomposition'
   do i = 1, n
      PetscCallA(PetscRandomGetValue(RdmCtx, M2D%XX, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, M2D%YY, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, M2D%XY, ierr))
      call SpectralDecomposition(M2D, ppleValues2D, ppleDirections2D)
      do j = 1, 2
         M2D = M2D - ppleValues2D(j) * ppleDirections2D(j)
      end do

      PetscCallA(PetscRandomGetValue(RdmCtx, M3D%XX, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, M3D%YY, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, M3D%ZZ, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, M3D%YZ, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, M3D%XZ, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, M3D%XY, ierr))
      call SpectralDecomposition(M3D, ppleValues3D, ppleDirections3D)
      do j = 1, 3
         M3D = M3D - ppleValues3D(j) * ppleDirections3D(j)
      end do
      write (*, *) i, norm(M2D), norm(M3D)
   end do

   write (*, *) 'Testing EigenVectorValues and MatRaRt'
   do i = 1, n
      PetscCallA(PetscRandomGetValue(RdmCtx, M2D%XX, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, M2D%YY, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, M2D%XY, ierr))
      call Diagonalize(M2D, P2D, D2D)
      M2D = M2D - MEF90MatRaRt(D2D, P2D)
      PetscCallA(PetscRandomGetValue(RdmCtx, M3D%XX, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, M3D%YY, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, M3D%ZZ, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, M3D%YZ, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, M3D%XZ, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, M3D%XY, ierr))
      call Diagonalize(M3D, P3D, D3D)
      M3D = M3D - MEF90MatRaRt(D3D, P3D)
      write (*, *) i, norm(M2D), norm(M3D)
   end do

   PetscCallA(PetscRandomDestroy(RdmCtx, ierr))
   PetscCallA(MEF90Finalize(ierr))
   PetscCallA(PetscFinalize(ierr))
end program TestSpectralDecomposition
