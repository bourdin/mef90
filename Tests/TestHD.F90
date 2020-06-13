Program TestSplit
#include "../MEF90/mef90.inc"
#include "../m_DefMech/mef90DefMech.inc"
#include "finclude/petscdef.h"
   Use m_MEF90
   Use MEF90_APPEND(m_MEF90_DefMechSplit,MEF90_DIM)D

   IMPLICIT NONE

   PetscReal                        :: E,nu
   Type(MEF90_HOOKESLAW)            :: A,ASph,ADev,APlus,AMinus

   Type(MEF90_MATS)                 :: M,MSph,MDev,DEEDPlus,DEEDMinus
   Type(MEF90_MAT)                  :: Pinv,PPinv
   PetscReal                        :: EED,EEDSph,EEDDev,EEDPlus,EEDMinus
   PetscReal                        :: Ndim = MEF90_DIM
   PetscReal                        :: lambda,mu

   PetscBool                        :: flg
   Character(len=1024)              :: IOBuffer
   PetscInt                         :: ierr
   PetscInt                         :: i,j,n=10
   PetscRandom                      :: RdmCtx
   PetscReal                        :: eps = 2.0_Kr * epsilon(1.0_Kr)
   PetscReal                        :: gamma
   PetscInt,dimension(8)            :: iSeed
   Integer,Parameter                :: sizeOfMatS = SIZEOFMEF90_MATS
   Integer,Parameter                :: dim = MEF90_DIM

   Class(MEF90_DEFMECHSPLIT),Allocatable :: Split

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


   gamma = 1.0e-10
   split = MEF90_DEFMECHSPLITHD(gamma)

   !A%type = MEF90HookesLawTypeFull
   A%type = MEF90HookesLawTypeIsotropic
   A%YoungsModulus = E
   A%PoissonRatio  = nu

   ASph%type = MEF90HookesLawTypeFull
   ASph%type = MEF90HookesLawTypeIsotropic
   ADev%type = MEF90HookesLawTypeFull
   ADev%type = MEF90HookesLawTypeIsotropic
   ASph%fullTensor = 0.0_Kr
   ADev%fullTensor = 0.0_Kr

#if MEF90_DIM == 2
   A%IsPlaneStress = .false.
   If (A%IsPlaneStress) Then
      lambda = E * nu / (1.0_Kr - nu**2)
      mu     = E / (1.0_Kr + nu) * .5_Kr
   Else
      lambda = E * nu / (1.0_Kr + nu) / (1.0_Kr - 2.0_Kr * nu)
      mu     = E / (1.0_Kr + nu) * .5_Kr
   End If
#else
   lambda = E * nu / (1.0_Kr + nu) / (1 - 2.0_Kr * nu)
   mu     = E / (1.0_Kr + nu) * .5_Kr
#endif
   A%lambda    = lambda
   A%mu        = mu
   ASph%lambda = lambda + 2.0_Kr * mu / Ndim
   ASph%mu     = 0.0_Kr
   ADev%lambda = (1.0_Kr - Ndim) * mu / Ndim
   ADev%mu     = mu

   Write(*,*) 'lambda', A%lambda, ASph%lambda, ADev%lambda
   Write(*,*) 'mu    ', A%mu,     ASph%mu,     ADev%mu

   Do i = 0, n-1
      Call PetscRandomGetValue(RdmCtx,M%XX,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M%YY,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M%XY,ierr);CHKERRQ(ierr)
#if MEF90_DIM == 3
      Call PetscRandomGetValue(RdmCtx,M%ZZ,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M%YZ,ierr);CHKERRQ(ierr)
      Call PetscRandomGetValue(RdmCtx,M%XZ,ierr);CHKERRQ(ierr)
#endif
      Write(*,'(A,<sizeOfMatS>(E12.5,2x))') 'M: ', M


      MSph   = HydrostaticPart(M)
      MDev   = DeviatoricPart(M)
      EED    = ((A * M) .dotP. M) * 0.5_Kr
      EEDSph = ((A * MSph) .dotP. MSph) * 0.5_Kr
      EEDDev = ((A * MDev) .dotP. MDev) * 0.5_Kr
     
      Write(*,'(A,4(E12.5,2x))') 'W^S, W^D, W^S + W^D , W ', EEDSph, EEDDev, EEDSph + EEDDev, EED
      EED    = (A%lambda * trace(M)**2 + 2.0_Kr * A%mu * Trace(M*M)) * 0.5_Kr
      EEDSph = (A%lambda + 2.0_Kr * A%mu / Ndim) * trace(M)**2 * 0.5_Kr
      EEDDev = A%mu * (trace(M*M) - trace(M)**2/Ndim)
      Write(*,'(A,4(E12.5,2x))') 'W^S, W^D, W^S + W^D , W ', EEDSph, EEDDev, EEDSph + EEDDev, EED
      Write(*,'(A,3(E12.5,2x))') 'W^S, A^SM.M, A^S M^S.M^S', EEDSph, ((ASph * M) .dotP. M) * 0.5_Kr, ((ASph * MSph) .dotP. MSph) * 0.5_Kr
      Write(*,'(A,3(E12.5,2x))') 'W^D, A^DM.M, A^D M^D.M^D', EEDDev, ((ADev * M) .dotP. M) * 0.5_Kr, ((ADev * MDev) .dotP. MDev) * 0.5_Kr
      If (trace(M) > 0.0_Kr) Then
         EEDPlus  = EEDSph + EEDDev
         EEDMinus = 0.0_Kr
      Else
         EEDPlus  = EEDDev
         EEDMinus = EEDSph
      End If
      Write(*,'(A,(E12.5,2x))') '  tr(M)      ', trace(M)
      Write(*,'(A,4(E12.5,2x))') 'Direct: EED, EED^+, EED^-, Sum', EED, EEDPlus, EEDMinus, EEDPlus + EEDMinus

      Call Split%EED(M,A,EEDPlus,EEDMinus)
      Write(*,'(A,4(E12.5,2x))') '   EED: EED, EED^+, EED^-, Sum', EED, EEDPlus, EEDMinus, EEDPlus + EEDMinus

      Call Split%DEED(M,A,DEEDPlus,DEEDMinus)
      EEDPlus  = (DEEDPlus  .dotP. M) * 0.5_Kr
      EEDMinus = (DEEDMinus .dotP. M) * 0.5_Kr
      Write(*,'(A,4(E12.5,2x))') '  DEED: EED, EED^+, EED^-, Sum', EED, EEDPlus, EEDMinus, EEDPlus + EEDMinus

      Call Split%D2EED(M,A,APlus,AMinus)
      EEDPlus  = ((APlus  * M) .dotP. M) * 0.5_Kr
      EEDMinus = ((AMinus * M) .dotP. M) * 0.5_Kr
      Write(*,'(A,4(E12.5,2x))') 'D2EED:  EED, EED^+, EED^-, Sum', EED, EEDPlus, EEDMinus, EEDPlus + EEDMinus
      Write(*,'(A,2(E12.5,2x))') 'A   lambda, mu                ', A%lambda,A%mu
      Write(*,'(A,2(E12.5,2x))') 'A^+ lambda, mu                ', APlus%lambda,APlus%mu
      Write(*,'(A,2(E12.5,2x))') 'A^- lambda, mu                ', AMinus%lambda,AMinus%mu

      Write(*,'(A,<sizeOfMatS>(E12.5,2x))') ' Direct: A^+ M: ', (APlus  * M)
      Write(*,'(A,<sizeOfMatS>(E12.5,2x))') '   DEED: A^+ M: ', DEEDPlus
      Write(*,'(A,<sizeOfMatS>(E12.5,2x))') ' Direct: A^- M: ', (AMinus  * M)
      Write(*,'(A,<sizeOfMatS>(E12.5,2x))') '   DEED: A^- M: ', DEEDMinus
      Write(*,*) 
      !Write(*,*) A%fullTensor - ADev%FullTensor - ASph%fullTensor

      !Write(*,'(A,3(E12.5,2x))') '    1/2 A^+ M M, A^- M M: ', APlus * M .dotP. M * 0.5_Kr, AMinus * M .dotP. M * 0.5_Kr, (APlus * M .dotP. M) * 0.5_Kr + (AMinus * M .dotP. M) * 0.5_Kr

      Write(*,*)
   End Do

   Call PetscRandomDestroy(RdmCtx,ierr);CHKERRQ(ierr)
   Call MEF90Finalize(ierr)
   Call PetscFinalize()
End Program TestSplit
