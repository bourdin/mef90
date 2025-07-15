program TestSplit
#include "../MEF90/mef90.inc"
#include "../m_DefMech/mef90DefMech.inc"
#include "petsc/finclude/petsc.h"
   use m_MEF90
   use MEF90_APPEND(m_MEF90_DefMechSplit,MEF90_DIM)D

   implicit none(type, external)

   PetscReal                        :: E, nu
   type(MEF90_HOOKESLAW)            :: A, ASph, ADev, APlus, AMinus

   type(MEF90_MATS)                 :: M, MSph, MDev, DEEDPlus, DEEDMinus
   type(MEF90_MAT)                  :: Pinv, PPinv
   PetscReal                        :: EED, EEDSph, EEDDev, EEDPlus, EEDMinus
   PetscReal                        :: Ndim = MEF90_DIM
   PetscReal                        :: lambda, mu

   PetscBool                        :: flg
   character(len=1024)              :: IOBuffer
   PetscInt                         :: ierr
   PetscInt                         :: i, j, n = 10
   PetscRandom                      :: RdmCtx
   PetscReal                        :: eps = 2.0_kr * epsilon(1.0_kr)
   PetscReal                        :: gamma
   PetscInt, dimension(8)            :: iSeed
   integer, parameter                :: sizeOfMatS = SIZEOFMEF90_MATS
   integer, parameter                :: dim = MEF90_DIM

   class(MEF90_DEFMECHSPLIT), allocatable :: Split

   PetscCallA(PetscInitialize(ierr))
   PetscCallA(MEF90Initialize(PETSC_COMM_WORLD, ierr))
   PetscCallA(PetscRandomCreate(PETSC_COMM_WORLD, RdmCtx, ierr))
   PetscCallA(PetscRandomSetFromOptions(RdmCtx, ierr))
   PetscCallA(PetscRandomSetInterval(RdmCtx, -1.0_kr, 1.0_kr, ierr))
   call date_and_time(VALUES=iSeed)
   PetscCallA(PetscRandomSetSeed(RdmCtx, iSeed(8) * 1000 * iseed(7), ierr))
   PetscCallA(PetscRandomSeed(RdmCtx, ierr))

   PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-n', n, flg, ierr))
   E = 1.0_kr
   PetscCallA(PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-E', E, flg, ierr))
   nu = .3_kr
   PetscCallA(PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-nu', nu, flg, ierr))

   gamma = 1.0e-10
   split = MEF90_DEFMECHSPLITHD(gamma)

   !A%type = MEF90HookesLawTypeFull
   A % type = MEF90HookesLawTypeIsotropic
   A % YoungsModulus = E
   A % PoissonRatio = nu

   ASph % type = MEF90HookesLawTypeFull
   ASph % type = MEF90HookesLawTypeIsotropic
   ADev % type = MEF90HookesLawTypeFull
   ADev % type = MEF90HookesLawTypeIsotropic
   ASph % fullTensor = 0.0_kr
   ADev % fullTensor = 0.0_kr

#if MEF90_DIM == 2
   A % IsPlaneStress = .false.
   if (A % IsPlaneStress) then
      lambda = E * nu / (1.0_kr - nu**2)
      mu = E / (1.0_kr + nu)*.5_kr
   else
      lambda = E * nu / (1.0_kr + nu) / (1.0_kr - 2.0_kr * nu)
      mu = E / (1.0_kr + nu)*.5_kr
   end if
#else
   lambda = E * nu / (1.0_kr + nu) / (1 - 2.0_kr * nu)
   mu = E / (1.0_kr + nu)*.5_kr
#endif
   A % lambda = lambda
   A % mu = mu
   ASph % lambda = lambda + 2.0_kr * mu / Ndim
   ASph % mu = 0.0_kr
   ADev % lambda = (1.0_kr - Ndim) * mu / Ndim
   ADev % mu = mu

   write (*, *) 'lambda', A % lambda, ASph % lambda, ADev % lambda
   write (*, *) 'mu    ', A % mu, ASph % mu, ADev % mu

   do i = 0, n - 1
      PetscCallA(PetscRandomGetValue(RdmCtx, M % XX, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, M % YY, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, M % XY, ierr))
#if MEF90_DIM == 3
      PetscCallA(PetscRandomGetValue(RdmCtx, M % ZZ, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, M % YZ, ierr))
      PetscCallA(PetscRandomGetValue(RdmCtx, M % XZ, ierr))
#endif
      write (*, '(A,6(E12.5,2x))') 'M: ', M

      MSph = HydrostaticPart(M)
      MDev = DeviatoricPart(M)
      EED = ((A * M) .dotP.M) * 0.5_kr
      EEDSph = ((A * MSph) .dotP.MSph) * 0.5_kr
      EEDDev = ((A * MDev) .dotP.MDev) * 0.5_kr

      write (*, '(A,4(E12.5,2x))') 'W^S, W^D, W^S + W^D , W ', EEDSph, EEDDev, EEDSph + EEDDev, EED
      EED = (A % lambda * trace(M)**2 + 2.0_kr * A % mu * Trace(M * M)) * 0.5_kr
      EEDSph = (A % lambda + 2.0_kr * A % mu / Ndim) * trace(M)**2 * 0.5_kr
      EEDDev = A % mu * (trace(M * M) - trace(M)**2 / Ndim)
      write (*, '(A,4(E12.5,2x))') 'W^S, W^D, W^S + W^D , W ', EEDSph, EEDDev, EEDSph + EEDDev, EED
      write (*, '(A,3(E12.5,2x))') 'W^S, A^SM.M, A^S M^S.M^S', EEDSph, ((ASph * M) .dotP.M) * 0.5_kr, ((ASph * MSph) .dotP.MSph) * 0.5_kr
      write (*, '(A,3(E12.5,2x))') 'W^D, A^DM.M, A^D M^D.M^D', EEDDev, ((ADev * M) .dotP.M) * 0.5_kr, ((ADev * MDev) .dotP.MDev) * 0.5_kr
      if (trace(M) > 0.0_kr) then
         EEDPlus = EEDSph + EEDDev
         EEDMinus = 0.0_kr
      else
         EEDPlus = EEDDev
         EEDMinus = EEDSph
      end if
      write (*, '(A,(E12.5,2x))') '  tr(M)      ', trace(M)
      write (*, '(A,4(E12.5,2x))') 'Direct: EED, EED^+, EED^-, Sum', EED, EEDPlus, EEDMinus, EEDPlus + EEDMinus

      call Split % EED(M, A, EEDPlus, EEDMinus)
      write (*, '(A,4(E12.5,2x))') '   EED: EED, EED^+, EED^-, Sum', EED, EEDPlus, EEDMinus, EEDPlus + EEDMinus

      call Split % DEED(M, A, DEEDPlus, DEEDMinus)
      EEDPlus = (DEEDPlus.dotP.M) * 0.5_kr
      EEDMinus = (DEEDMinus.dotP.M) * 0.5_kr
      write (*, '(A,4(E12.5,2x))') '  DEED: EED, EED^+, EED^-, Sum', EED, EEDPlus, EEDMinus, EEDPlus + EEDMinus

      call Split % D2EED(M, A, APlus, AMinus)
      EEDPlus = ((APlus * M) .dotP.M) * 0.5_kr
      EEDMinus = ((AMinus * M) .dotP.M) * 0.5_kr
      write (*, '(A,4(E12.5,2x))') 'D2EED:  EED, EED^+, EED^-, Sum', EED, EEDPlus, EEDMinus, EEDPlus + EEDMinus
      write (*, '(A,2(E12.5,2x))') 'A   lambda, mu                ', A % lambda, A % mu
      write (*, '(A,2(E12.5,2x))') 'A^+ lambda, mu                ', APlus % lambda, APlus % mu
      write (*, '(A,2(E12.5,2x))') 'A^- lambda, mu                ', AMinus % lambda, AMinus % mu

      write (*, '(A,6(E12.5,2x))') ' Direct: A^+ M: ', (APlus * M)
      write (*, '(A,6(E12.5,2x))') '   DEED: A^+ M: ', DEEDPlus
      write (*, '(A,6(E12.5,2x))') ' Direct: A^- M: ', (AMinus * M)
      write (*, '(A,6(E12.5,2x))') '   DEED: A^- M: ', DEEDMinus
      write (*, *)
      !Write(*,*) A%fullTensor - ADev%FullTensor - ASph%fullTensor

      !Write(*,'(A,3(E12.5,2x))') '    1/2 A^+ M M, A^- M M: ', APlus * M .dotP. M * 0.5_Kr, AMinus * M .dotP. M * 0.5_Kr, (APlus * M .dotP. M) * 0.5_Kr + (AMinus * M .dotP. M) * 0.5_Kr

      write (*, *)
   end do

   PetscCallA(PetscRandomDestroy(RdmCtx, ierr))
   PetscCallA(MEF90Finalize(ierr))
   PetscCallA(PetscFinalize(ierr))
end program TestSplit
