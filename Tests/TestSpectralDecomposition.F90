Program TestSpectralDecomposition
#include "petsc/finclude/petsc.h"
   Use m_MEF90
   IMPLICIT NONE

   PetscReal                        :: E, nu
   Type(MatS2D)                     :: M2D,D2D
   Type(Mat2D)                      :: P2D
   Type(MatS3D)                     :: M3D,D3D
   Type(Mat3D)                      :: P3D
   Type(MatS2D),Dimension(2)        :: ppleDirections2D
   Type(MatS3D),Dimension(3)        :: ppleDirections3D
   PetscReal,Dimension(2)           :: ppleValues2D
   PetscReal,Dimension(3)           :: ppleValues3D
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
   Write(*,*) 'Testing SpectralDecomposition'
   Do i = 1,n
      Call PetscRandomGetValue(RdmCtx,M2D%XX,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M2D%YY,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M2D%XY,ierr);CHKERRQ(ierr)
      Call SpectralDecomposition(M2D,ppleValues2D,ppleDirections2D)
      Do j = 1,2
         M2D = M2D - ppleValues2D(j) * ppleDirections2D(j)   
      End Do

      Call PetscRandomGetValue(RdmCtx,M3D%XX,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M3D%YY,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M3D%ZZ,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M3D%YZ,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M3D%XZ,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M3D%XY,ierr);CHKERRQ(ierr)
      Call SpectralDecomposition(M3D,ppleValues3D,ppleDirections3D)
      Do j = 1,3
         M3D = M3D - ppleValues3D(j) * ppleDirections3D(j)   
      End Do
      Write(*,*) i, norm(M2D),norm(M3D)
   End Do

   Write(*,*) 'Testing EigenVectorValues and MatRaRt'
   Do i = 1,n
      Call PetscRandomGetValue(RdmCtx,M2D%XX,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M2D%YY,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M2D%XY,ierr);CHKERRQ(ierr)
      Call Diagonalize(M2D,P2D,D2D)
      M2D = M2D - MatRaRt(D2D,P2D)
      Call PetscRandomGetValue(RdmCtx,M3D%XX,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M3D%YY,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M3D%ZZ,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M3D%YZ,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M3D%XZ,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M3D%XY,ierr);CHKERRQ(ierr)
      Call Diagonalize(M3D,P3D,D3D)
      M3D = M3D - MatRaRt(D3D,P3D)
      Write(*,*) i, norm(M2D),norm(M3D)
   End Do

   Call PetscRandomDestroy(RdmCtx,ierr);CHKERRQ(ierr)
   Call MEF90Finalize(ierr)
   Call PetscFinalize()
End Program TestSpectralDecomposition