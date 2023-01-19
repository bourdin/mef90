Program TestSpectralDecomposition
#include "petsc/finclude/petsc.h"
   Use m_MEF90
   IMPLICIT NONE

   Type(MatS2D)                     :: M2D,D2D
   Type(Mat2D)                      :: P2D
   Type(MatS3D)                     :: M3D,D3D
   Type(Mat3D)                      :: P3D
   Type(MatS2D),Dimension(2)        :: ppleDirections2D
   Type(MatS3D),Dimension(3)        :: ppleDirections3D
   PetscReal,Dimension(2)           :: ppleValues2D
   PetscReal,Dimension(3)           :: ppleValues3D
   PetscBool                        :: flg
   PetscErrorCode                   :: ierr
   PetscInt                         :: i,j,n=10
   PetscRandom                      :: RdmCtx

   PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
   PetscCallA(MEF90Initialize(ierr))
   PetscCallA(PetscRandomCreate(PETSC_COMM_WORLD,RdmCtx,ierr))
   PetscCallA(PetscRandomSetFromOptions(RdmCtx,ierr))
   PetscCallA(PetscRandomSetInterval(RdmCtx,-1.0_Kr,1.0_Kr,ierr))
   
   PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,'','-n',n,flg,ierr))
   Write(*,*) 'Testing SpectralDecomposition'
   Do i = 1,n
      PetscCallA(PetscRandomGetValue(RdmCtx,M2D%XX,ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx,M2D%YY,ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx,M2D%XY,ierr))
      Call SpectralDecomposition(M2D,ppleValues2D,ppleDirections2D)
      Do j = 1,2
         M2D = M2D - ppleValues2D(j) * ppleDirections2D(j)   
      End Do

      PetscCallA(PetscRandomGetValue(RdmCtx,M3D%XX,ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx,M3D%YY,ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx,M3D%ZZ,ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx,M3D%YZ,ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx,M3D%XZ,ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx,M3D%XY,ierr))
      Call SpectralDecomposition(M3D,ppleValues3D,ppleDirections3D)
      Do j = 1,3
         M3D = M3D - ppleValues3D(j) * ppleDirections3D(j)   
      End Do
      Write(*,*) i, norm(M2D),norm(M3D)
   End Do

   Write(*,*) 'Testing EigenVectorValues and MatRaRt'
   Do i = 1,n
      PetscCallA(PetscRandomGetValue(RdmCtx,M2D%XX,ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx,M2D%YY,ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx,M2D%XY,ierr))
      Call Diagonalize(M2D,P2D,D2D)
      M2D = M2D - MEF90MatRaRt(D2D,P2D)
      PetscCallA(PetscRandomGetValue(RdmCtx,M3D%XX,ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx,M3D%YY,ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx,M3D%ZZ,ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx,M3D%YZ,ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx,M3D%XZ,ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx,M3D%XY,ierr))
      Call Diagonalize(M3D,P3D,D3D)
      M3D = M3D - MEF90MatRaRt(D3D,P3D)
      Write(*,*) i, norm(M2D),norm(M3D)
   End Do

   PetscCallA(PetscRandomDestroy(RdmCtx,ierr))
   PetscCallA(MEF90Finalize(ierr))
   PetscCallA(PetscFinalize(ierr))
End Program TestSpectralDecomposition