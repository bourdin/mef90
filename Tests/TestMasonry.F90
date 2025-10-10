program TestSplit
#include "../MEF90/mef90.inc"
#include "../m_DefMech/mef90DefMech.inc"
#include "petsc/finclude/petsc.h"
   use m_MEF90
   use MEF90_APPEND(m_MEF90_DefMechSplit,MEF90_DIM)D

   implicit none(type, external)

   PetscReal                        :: E, nu
   type(MEF90_HOOKESLAW)            :: A, APlus, AMinus

   type(MEF90_MATS)                 :: M, Sigma, SigmaPlus, SigmaMinus, D, DPlus, Dminus
   type(MEF90_MAT)                  :: Pinv, PPinv
   PetscReal                        :: EED, EEDPlus, EEDMinus

   PetscBool                        :: flg, mef90
   character(len=1024)              :: IOBuffer
   PetscErrorCode                   :: ierr
   PetscInt                         :: i, j, n = 10
   PetscRandom                      :: RdmCtx
   PetscReal                        :: eps = 2.0_kr * epsilon(1.0_kr)
   PetscInt, dimension(8)            :: iSeed
   integer, parameter                :: sizeOfMatS = SIZEOFMEF90_MATS
   integer, parameter                :: dim = MEF90_DIM
   PetscReal                        :: gamma

   class(MEF90_DEFMECHSPLIT), allocatable :: Split

   call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
   call MEF90Initialize(PETSC_COMM_WORLD, ierr)
   call PetscRandomCreate(PETSC_COMM_WORLD, RdmCtx, ierr); CHKERRQ(ierr)
   call PetscRandomSetFromOptions(RdmCtx, ierr); CHKERRQ(ierr)
   call PetscRandomSetInterval(RdmCtx, -1.0_kr, 1.0_kr, ierr); CHKERRQ(ierr)
   call date_and_time(VALUES=iSeed)
   call PetscRandomSetSeed(RdmCtx, iSeed(8) * 1000 * iseed(7), ierr); CHKERRQ(ierr)
   call PetscRandomSeed(RdmCtx, ierr); CHKERRQ(ierr)

   call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-n', n, flg, ierr); CHKERRQ(ierr)
   E = 1.0_kr
   call PetscOptionsGetReal(PETSC_NULL_CHARACTER, '-E', E, flg, ierr); CHKERRQ(ierr); 
   nu = .3_kr
   call PetscOptionsGetReal(PETSC_NULL_CHARACTER, '-nu', nu, flg, ierr); CHKERRQ(ierr); 
   write (*, *) 'Masonry split:'
   !split = MEF90_DEFMECHMASONRYSPLIT()
   gamma = 1.0e-5
   split = MEF90_DEFMECHSPLITHD(gamma)

   A%type = MEF90HookesLawTypeIsotropic
   A%YoungsModulus = E
   A%PoissonRatio = nu
#if MEF90_DIM == 2
   A%IsPlaneStress = .false.
   if (A%IsPlaneStress) then
      A%lambda = E * nu / (1.0_kr - nu**2)
      A%mu = E / (1.0_kr + nu)*.5_kr
   else
      A%lambda = E * nu / (1.0_kr + nu) / (1.0_kr - 2.0_kr * nu)
      A%mu = E / (1.0_kr + nu)*.5_kr
   end if
#else
   A%lambda = E * nu / (1.0_kr + nu) / (1 - 2.0_kr * nu)
   A%mu = E / (1.0_kr + nu)*.5_kr
#endif

   do i = 0, n - 1
      call PetscRandomGetValue(RdmCtx, M%XX, ierr); CHKERRQ(ierr)
      call PetscRandomGetValue(RdmCtx, M%YY, ierr); CHKERRQ(ierr)
      call PetscRandomGetValue(RdmCtx, M%XY, ierr); CHKERRQ(ierr)
#if MEF90_DIM == 3
      call PetscRandomGetValue(RdmCtx, M%ZZ, ierr); CHKERRQ(ierr)
      call PetscRandomGetValue(RdmCtx, M%YZ, ierr); CHKERRQ(ierr)
      call PetscRandomGetValue(RdmCtx, M%XZ, ierr); CHKERRQ(ierr)
#endif
      write (*, '(A,<sizeOfMatS>(E12.5,2x))') 'M: ', M
      call Diagonalize(A * M, Pinv, D)

      EED = 0.5_kr * ((A * M) .dotP.M)
      call Split%EED(M, A, EEDPlus, EEDMinus)
      call Split%DEED(M, A, Sigmaplus, Sigmaminus)
      call Split%D2EED(M, A, APlus, AMinus)
      write (*, '(A,<sizeOfMatS>(E12.5,2x))') '                   A^+ M: ', APlus * M
      write (*, '(A,<sizeOfMatS>(E12.5,2x))') '                 Sigma^+: ', SigmaPlus
      call Diagonalize(SigmaPlus, Pinv, DPlus)

      write (*, '(A,<sizeOfMatS>(E12.5,2x))') '                   A^- M: ', AMinus * M
      write (*, '(A,<sizeOfMatS>(E12.5,2x))') '              SigmaMinus: ', SigmaMinus

      call Diagonalize(SigmaMinus, Pinv, DMinus)
      ! If (Trace(M) > 0.0_Kr) Then
      !    write(*,'(A,3(ES12.5,2x))') '1/2 AM.M, EED+, EED-', EED, 0.5_Kr * ((A * M) .dotP. M)
      ! End If
#if MEF90_DIM == 2
      write (*, '(A,2(E12.5,2x))') '      Principal stresses: ', D%XX, D%YY
      write (*, '(A,2(E12.5,2x))') '    Sigma^+ (pple basis): ', DPlus%XX, DPlus%YY
      write (*, '(A,2(E12.5,2x))') '    Sigma^- (pple basis): ', DMinus%XX, DMinus%YY
      write (*, '(A,2(E12.5,2x))') '        Sum (pple basis): ', DMinus%XX + DPlus%XX, DMinus%YY + DPlus%YY
#else
      write (*, '(A,3(E12.5,2x))') '      Principal stresses: ', D%XX, D%YY, D%ZZ
      !Write(*,'(A,3(E12.5,2x))') '    Sigma^+ (pple basis): ', DPlus%XX,DPlus%YY,DPlus%ZZ
      !Write(*,'(A,3(E12.5,2x))') '    Sigma^- (pple basis): ', DMinus%XX,DMinus%YY,DMinus%ZZ
      !Write(*,'(A,3(E12.5,2x))') '        Sum (pple basis): ', DMinus%XX + DPlus%XX,DMinus%YY + DPlus%YY,DMinus%ZZ + DPlus%ZZ
#endif
      !Write(*,'(A,<sizeOfMatS>(E12.5,2x))') '           A^+ M + A^- M: ', APlus * M + AMinus * M
      !Write(*,'(A,<sizeOfMatS>(E12.5,2x))') '                      AM: ', A*M
      !Write(*,'(A,<sizeOfMatS>(E12.5,2x))') '       Sigma^+ + Sigma^-: ', SigmaPlus + SigmaMinus

      write (*, '(A,3(E12.5,2x))') '         EED^+ EED^- EED: ', EEDPlus, EEDMinus, EED
      !Write(*,'(A,3(E12.5,2x))') '1/2 Sigma^+ M, Sigma^- M: ', SigmaPlus .dotP. M * 0.5_Kr, SigmaMinus .dotP. M * 0.5_Kr, (SigmaPlus + SigmaMinus) .dotP. M * 0.5_Kr

      !Write(*,'(A,3(E12.5,2x))') '    1/2 A^+ M M, A^- M M: ', APlus * M .dotP. M * 0.5_Kr, AMinus * M .dotP. M * 0.5_Kr, (APlus * M .dotP. M) * 0.5_Kr + (AMinus * M .dotP. M) * 0.5_Kr

      write (*, *)
   end do

   call PetscRandomDestroy(RdmCtx, ierr); CHKERRQ(ierr)
   call MEF90Finalize(ierr)
   call PetscFinalize()
end program TestSplit
