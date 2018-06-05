Program HookeLaw
#include "finclude/petscdef.h"
   Use m_MEF90
   IMPLICIT NONE

   PetscReal                        :: E, nu
   Type(Tens4OS2D)                  :: HookeLaw2D,HookeLawInv2D
   Type(Tens4OS3D)                  :: HookeLaw3D,HookeLawInv3D
   Type(MatS2D)                     :: sigma2D,Epsilon2D
   Type(MatS3D)                     :: sigma3D,Epsilon3D
   PetscBool                        :: flg,mef90
   Character(len=1024)              :: IOBuffer
   PetscInt                         :: ierr
   PetscInt                         :: i,n
   PetscRandom                      :: RdmCtx

   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   Call MEF90Initialize(ierr)
   Call PetscRandomCreate(PETSC_COMM_WORLD,RdmCtx,ierr);CHKERRQ(ierr)
   Call PetscRandomSetFromOptions(RdmCtx,ierr);CHKERRQ(ierr)
   Call PetscRandomSetInterval(RdmCtx,-1.0_Kr,1.0_Kr,ierr);CHKERRQ(ierr)
   
   E = 1.0_Kr
   Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-E',E,flg,ierr);CHKERRQ(ierr);
   nu = .3_Kr
   Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-nu',nu,flg,ierr);CHKERRQ(ierr);
   mef90 = PETSC_FALSE
   Call PetscOptionsGetBool(PETSC_NULL_CHARACTER,'-mef90',mef90,flg,ierr);CHKERRQ(ierr);
   n = 10
   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-n',n,flg,ierr);CHKERRQ(ierr)
   
   Call MEF90HookeLawIsoENu3D(HookeLaw3D,E,nu)
   HookeLawInv3D = Invert(HookeLaw3D)
   Call MEF90HookeLawIsoEnu2DPlaneStrain(HookeLaw2D,E,nu)
   HookeLawInv2D = Invert(HookeLaw2D)
   Do i = 1,n
      Call PetscRandomGetValue(RdmCtx,epsilon2D%XX,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,epsilon2D%YY,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,epsilon2D%XY,ierr);CHKERRQ(ierr)
      sigma2D = HookeLaw2D * epsilon2D
      epsilon2D = epsilon2D - HookeLawInv2D * sigma2D

      Call PetscRandomGetValue(RdmCtx,epsilon3D%XX,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,epsilon3D%YY,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,epsilon3D%ZZ,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,epsilon3D%YZ,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,epsilon3D%XZ,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,epsilon3D%XY,ierr);CHKERRQ(ierr)
      sigma3D = HookeLaw3D * epsilon3D
      epsilon3D = epsilon3D - HookeLawInv3D * sigma3D
      Write(*,*) i, norm(epsilon2D),norm(epsilon3D)
   End Do

   HookeLawInv2D = sqrt(HookeLaw2D)
   HookeLawInv3D = sqrt(HookeLaw3D)
   Do i = 1,n
      Call PetscRandomGetValue(RdmCtx,epsilon2D%XX,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,epsilon2D%YY,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,epsilon2D%XY,ierr);CHKERRQ(ierr)
      sigma2D = HookeLaw2D * epsilon2D
      epsilon2D = sigma2D - HookeLawInv2D * (HookeLawInv2D * epsilon2D)

      Call PetscRandomGetValue(RdmCtx,epsilon3D%XX,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,epsilon3D%YY,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,epsilon3D%ZZ,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,epsilon3D%YZ,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,epsilon3D%XZ,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,epsilon3D%XY,ierr);CHKERRQ(ierr)
      sigma3D = HookeLaw3D * epsilon3D
      epsilon3D = sigma3D - HookeLawInv3D * (HookeLawInv3D * epsilon3D)
      Write(*,*) i, norm(epsilon2D),norm(epsilon3D)
   End Do

   epsilon3D%XX = 1.0_Kr
   epsilon3D%YY = 1.0_Kr
   epsilon3D%ZZ = 1.0_Kr
   epsilon3D%YZ = 0.0_Kr
   epsilon3D%XZ = 0.0_Kr
   epsilon3D%XY = 0.0_Kr
   write(*,*) HookeLaw3D
   sigma3D = HookeLaw3D * epsilon3D
   write(*,*) 'epsilon', epsilon3D
   write(*,*) 'sigma  ', sigma3D


   Call PetscRandomDestroy(RdmCtx,ierr);CHKERRQ(ierr)
   Call MEF90Finalize(ierr)
   Call PetscFinalize()
End Program HookeLaw