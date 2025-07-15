program HookeLaw
#include <petsc/finclude/petsc.h>
   use m_MEF90
   use petsc
   implicit none(type, external)

   PetscReal                        :: E, nu
   type(Tens4OS2D)                  :: HookeLaw2D, HookeLawInv2D
   type(Tens4OS3D)                  :: HookeLaw3D, HookeLawInv3D
   type(MatS2D)                     :: sigma2D, Epsilon2D
   type(MatS3D)                     :: sigma3D, Epsilon3D
   PetscBool                        :: flg, mef90
   PetscErrorCode                   :: ierr
   PetscInt                         :: i, n
   PetscRandom                      :: RdmCtx

   PetscCallA(PetscInitialize(ierr))
   PetscCallA(MEF90Initialize(PETSC_COMM_WORLD, ierr))
   PetscCallA(PetscRandomCreate(PETSC_COMM_WORLD, RdmCtx, ierr))
   PetscCallA(PetscRandomSetFromOptions(RdmCtx, ierr))
   PetscCallA(PetscRandomSetInterval(RdmCtx, -1.0_kr, 1.0_kr, ierr))

   E = 1.0_kr
   !PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-n',numMat,flg,ierr))
   PetscCallA(PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-E', E, flg, ierr))
   nu = .3_kr
   PetscCallA(PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-nu', nu, flg, ierr))
   mef90 = PETSC_FALSE
   PetscCallA(PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-mef90', mef90, flg, ierr))
   n = 10
   PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-n', n, flg, ierr))

   PetscCallA(MEF90HookeLawIsoENu3D(HookeLaw3D, E, nu))
   HookeLawInv3D = Invert(HookeLaw3D)
   PetscCallA(MEF90HookeLawIsoEnu2DPlaneStrain(HookeLaw2D, E, nu))
   HookeLawInv2D = Invert(HookeLaw2D)
   do i = 1, n
      PetscCallA(PetscRandomGetValue(RdmCtx, epsilon2D % XX, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, epsilon2D % YY, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, epsilon2D % XY, ierr))
      sigma2D = HookeLaw2D * epsilon2D
      epsilon2D = epsilon2D - HookeLawInv2D * sigma2D

      PetscCallA(PetscRandomGetValue(RdmCtx, epsilon3D % XX, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, epsilon3D % YY, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, epsilon3D % ZZ, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, epsilon3D % YZ, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, epsilon3D % XZ, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, epsilon3D % XY, ierr))
      sigma3D = HookeLaw3D * epsilon3D
      epsilon3D = epsilon3D - HookeLawInv3D * sigma3D
      write (*, *) i, norm(epsilon2D), norm(epsilon3D)
   end do

   HookeLawInv2D = sqrt(HookeLaw2D)
   HookeLawInv3D = sqrt(HookeLaw3D)
   do i = 1, n
      PetscCallA(PetscRandomGetValue(RdmCtx, epsilon2D % XX, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, epsilon2D % YY, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, epsilon2D % XY, ierr))
      sigma2D = HookeLaw2D * epsilon2D
      epsilon2D = sigma2D - HookeLawInv2D * (HookeLawInv2D * epsilon2D)

      PetscCallA(PetscRandomGetValue(RdmCtx, epsilon3D % XX, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, epsilon3D % YY, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, epsilon3D % ZZ, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, epsilon3D % YZ, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, epsilon3D % XZ, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, epsilon3D % XY, ierr))
      sigma3D = HookeLaw3D * epsilon3D
      epsilon3D = sigma3D - HookeLawInv3D * (HookeLawInv3D * epsilon3D)
      write (*, *) i, norm(epsilon2D), norm(epsilon3D)
   end do

   epsilon3D % XX = 1.0_kr
   epsilon3D % YY = 1.0_kr
   epsilon3D % ZZ = 1.0_kr
   epsilon3D % YZ = 0.0_kr
   epsilon3D % XZ = 0.0_kr
   epsilon3D % XY = 0.0_kr
   write (*, *) HookeLaw3D
   sigma3D = HookeLaw3D * epsilon3D
   write (*, *) 'epsilon', epsilon3D
   write (*, *) 'sigma  ', sigma3D

   PetscCallA(PetscRandomDestroy(RdmCtx, ierr))
   PetscCallA(MEF90Finalize(ierr))
   PetscCallA(PetscFinalize(ierr))
end program HookeLaw
