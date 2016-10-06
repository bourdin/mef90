Program TestMasonryProjection
#include "finclude/petscdef.h"
   Use m_MEF90
   IMPLICIT NONE

   PetscReal                        :: E,nu
   Type(MEF90HookesLaw2D)           :: A2D
   Type(MEF90HookesLaw3D)           :: A3D

   Type(MatS2D)                     :: M2D,M2DPlus,M2DMinus
   Type(MatS3D)                     :: M3D,M3DPlus,M3DMinus

   PetscBool                        :: flg,mef90
   Character(len=1024)              :: IOBuffer
   PetscInt                         :: ierr
   PetscInt                         :: i,j,n=10
   PetscRandom                      :: RdmCtx

   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   Call MEF90Initialize(ierr)
   Call PetscRandomCreate(PETSC_COMM_WORLD,RdmCtx,ierr);CHKERRQ(ierr)
   Call PetscRandomSetFromOptions(RdmCtx,ierr);CHKERRQ(ierr)
   Call PetscRandomSetInterval(RdmCtx,-1.0_Kr,1.0_Kr,ierr);CHKERRQ(ierr)
   
   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-n',n,flg,ierr);CHKERRQ(ierr)
   E = 1.0_Kr
   Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-E',E,flg,ierr);CHKERRQ(ierr);
   nu = .3_Kr
   Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-nu',nu,flg,ierr);CHKERRQ(ierr);

   A2D%type = MEF90HookesLawTypeIsotropic
   !! plane stress
   A2D%lambda = E * nu / (1.0_Kr - nu**2) 
   A2D%mu     = E / (1.0_Kr + nu) * .5_Kr
   A2D%YoungsModulus = E
   A2D%PoissonRatio  = nu

   A3D%type = MEF90HookesLawTypeIsotropic
   A3D%lambda = E * nu / (1.0_Kr + nu) / (1.0_Kr - 2.0_Kr * nu)
   A3D%mu     = E / (1.0_Kr + nu) * .5_Kr    
   A3D%YoungsModulus = E
   A3D%PoissonRatio  = nu

   Do i = 1,n
      Call PetscRandomGetValue(RdmCtx,M2D%XX,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M2D%YY,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M2D%XY,ierr);CHKERRQ(ierr)
      Call MasonryProjection(M2D,A2D,M2DPlus,M2DMinus)
      Write(*,*)i, (A2D * M2DPlus) .DotP. M2DMinus

      Call PetscRandomGetValue(RdmCtx,M3D%XX,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M3D%YY,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M3D%ZZ,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M3D%YZ,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M3D%XZ,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M3D%XY,ierr);CHKERRQ(ierr)
      Call MasonryProjection(M3D,A3D,M3DPlus,M3DMinus)
      Write(*,*)i, (A3D * M3DPlus) .DotP. M3DMinus
   End Do


   Call PetscRandomDestroy(RdmCtx,ierr);CHKERRQ(ierr)
   Call MEF90Finalize(ierr)
   Call PetscFinalize()
End Program TestMasonryProjection