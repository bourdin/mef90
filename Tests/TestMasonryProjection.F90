Program TestMasonryProjection
#include "../MEF90/mef90.inc"
#include "finclude/petscdef.h"
   Use m_MEF90
   Use MEF90_APPEND(m_MEF90_DefMechMasonry,MEF90_DIM)D

   IMPLICIT NONE

   PetscReal                        :: E,nu
   Type(MEF90_HOOKESLAW)            :: A,APlus,AMinus

   Type(MEF90_MATS)                 :: M,Sigma,SigmaPlus,SigmaMinus,D,DPlus,Dminus
   Type(MEF90_MAT)                  :: Pinv,PPinv
   PetscReal                        :: EED,EEDPlus,EEDMinus

   PetscBool                        :: flg,mef90
   Character(len=1024)              :: IOBuffer
   PetscInt                         :: ierr
   PetscInt                         :: i,j,n=10
   PetscRandom                      :: RdmCtx
   PetscReal                        :: eps = 2.0_Kr * epsilon(1.0_Kr)
   PetscInt,dimension(8)            :: iSeed
   Integer,Parameter                :: sizeOfMatS = SIZEOFMEF90_MATS
   Integer,Parameter                :: dim = MEF90_DIM

   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   Call MEF90Initialize(ierr)
   Call PetscRandomCreate(PETSC_COMM_WORLD,RdmCtx,ierr);CHKERRQ(ierr)
   Call PetscRandomSetFromOptions(RdmCtx,ierr);CHKERRQ(ierr)
   Call PetscRandomSetInterval(RdmCtx,-1.0_Kr,1.0_Kr,ierr);CHKERRQ(ierr)
   Call date_and_time(VALUES=iSeed)
   Call PetscRandomSetSeed(RdmCtx,iSeed(8) * 1000 * iseed(7),ierr);CHKERRQ(ierr)
   Call PetscRandomSeed(RdmCtx,ierr);CHKERRQ(ierr)

   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-n',n,flg,ierr);CHKERRQ(ierr)
   E = 1.0_Kr
   Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-E',E,flg,ierr);CHKERRQ(ierr);
   nu = .3_Kr
   Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-nu',nu,flg,ierr);CHKERRQ(ierr);

   A%type = MEF90HookesLawTypeIsotropic
   A%YoungsModulus = E
   A%PoissonRatio  = nu
#if MEF90_DIM == 2
   A%IsPlaneStress = .false.
   If (A%IsPlaneStress) Then
      A%lambda = E * nu / (1.0_Kr - nu**2)
      A%mu     = E / (1.0_Kr + nu) * .5_Kr
   Else
      A%lambda = E * nu / (1.0_Kr + nu) / (1.0_Kr - 2.0_Kr * nu)
      A%mu     = E / (1.0_Kr + nu) * .5_Kr
   End If
#else
   A%lambda = E * nu / (1.0_Kr + nu) / (1 - 2.0_Kr * nu)
   A%mu     = E / (1.0_Kr + nu) * .5_Kr
#endif

   Do i = 0, n-1
      Call PetscRandomGetValue(RdmCtx,M%XX,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M%YY,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M%XY,ierr);CHKERRQ(ierr)
      ! M%XX = cos(2.* i * PETSC_PI / (n-1.))
      ! M%YY = sin(2.* i * PETSC_PI / (n-1.))
      ! Write(*,'(''theta: '',F12.5)') 360. * i / (n-1.)

#if MEF90_DIM == 3
      Call PetscRandomGetValue(RdmCtx,M%ZZ,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M%YZ,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M%XZ,ierr);CHKERRQ(ierr)
#endif
      Write(*,'(A,<sizeOfMatS>(E12.5,2x))') 'M: ', M
      Call Diagonalize(A*M,Pinv,D)

      EED = 0.5_Kr * (A * M) .dotP. M
      Call EEDMasonry(M,A,EEDPlus,EEDMinus)
      Call DEEDMasonry(M,A,Sigmaplus,Sigmaminus)
      Call D2EEDMasonry(M,A,APlus,AMinus)
      Write(*,'(A,<sizeOfMatS>(E12.5,2x))') '                   A^+ M: ', APlus * M
      Write(*,'(A,<sizeOfMatS>(E12.5,2x))') '                 Sigma^+: ', SigmaPlus
      Call Diagonalize(SigmaPlus,Pinv,DPlus)

      Write(*,'(A,<sizeOfMatS>(E12.5,2x))') '                   A^- M: ', AMinus * M
      Write(*,'(A,<sizeOfMatS>(E12.5,2x))') '              SigmaMinus: ', SigmaMinus

      Call Diagonalize(SigmaMinus,Pinv,DMinus)
#if MEF90_DIM == 2
      Write(*,'(A,2(E12.5,2x))') '      Principal stresses: ', D%XX,D%YY
      Write(*,'(A,2(E12.5,2x))') '    Sigma^+ (pple basis): ', DPlus%XX,DPlus%YY
      Write(*,'(A,2(E12.5,2x))') '    Sigma^- (pple basis): ', DMinus%XX,DMinus%YY
      Write(*,'(A,2(E12.5,2x))') '        Sum (pple basis): ', DMinus%XX + DPlus%XX,DMinus%YY + DPlus%YY
#else
      Write(*,'(A,3(E12.5,2x))') '      Principal stresses: ', D%XX,D%YY,D%ZZ
      Write(*,'(A,3(E12.5,2x))') '    Sigma^+ (pple basis): ', DPlus%XX,DPlus%YY,DPlus%ZZ
      Write(*,'(A,3(E12.5,2x))') '    Sigma^- (pple basis): ', DMinus%XX,DMinus%YY,DMinus%ZZ
      Write(*,'(A,3(E12.5,2x))') '        Sum (pple basis): ', DMinus%XX + DPlus%XX,DMinus%YY + DPlus%YY,DMinus%ZZ + DPlus%ZZ
#endif
      Write(*,'(A,<sizeOfMatS>(E12.5,2x))') '           A^+ M + A^- M: ', APlus * M + AMinus * M
      Write(*,'(A,<sizeOfMatS>(E12.5,2x))') '                      AM: ', A*M
      Write(*,'(A,<sizeOfMatS>(E12.5,2x))') '       Sigma^+ + Sigma^-: ', SigmaPlus + SigmaMinus

      Write(*,'(A,3(E12.5,2x))') '         EED^+ EED^- EED: ', EEDPlus, EEDMinus, EED
      Write(*,'(A,3(E12.5,2x))') '1/2 Sigma^+ M, Sigma^- M: ', SigmaPlus .dotP. M * 0.5_Kr, SigmaMinus .dotP. M * 0.5_Kr, (SigmaPlus + SigmaMinus) .dotP. M * 0.5_Kr 

      Write(*,'(A,3(E12.5,2x))') '    1/2 A^+ M M, A^- M M: ', APlus * M .dotP. M * 0.5_Kr, AMinus * M .dotP. M * 0.5_Kr, (APlus * M .dotP. M) * 0.5_Kr + (AMinus * M .dotP. M) * 0.5_Kr

      Write(*,*)
   End Do

   Call PetscRandomDestroy(RdmCtx,ierr);CHKERRQ(ierr)
   Call MEF90Finalize(ierr)
   Call PetscFinalize()
End Program TestMasonryProjection