#include "../MEF90/mef90.inc"
module m_MEF90_HeatXfer
#include "petsc/finclude/petsc.h"
   use petscsnes
   use petsctao
   use m_MEF90_EXO
   use m_MEF90_HeatXfer_class
   use, intrinsic :: iso_c_binding
   use m_MEF90_HeatXferAssembly2D, &
      MEF90HeatXferEnergy2D => MEF90HeatXferEnergy, &
      MEF90HeatXferOperator2D => MEF90HeatXferOperator, &
      MEF90HeatXferBilinearForm2D => MEF90HeatXferBilinearForm, &
      MEF90HeatXferIFunction2D => MEF90HeatXferIFunction, &
      MEF90HeatXferIJacobian2D => MEF90HeatXferIJacobian
   use m_MEF90_HeatXferAssembly3D, &
      MEF90HeatXferEnergy3D => MEF90HeatXFerEnergy, &
      MEF90HeatXferOperator3D => MEF90HeatXferOperator, &
      MEF90HeatXferBilinearForm3D => MEF90HeatXferBilinearForm, &
      MEF90HeatXferIFunction3D => MEF90HeatXferIFunction, &
      MEF90HeatXferIJacobian3D => MEF90HeatXferIJacobian

   implicit none(type)

   !Private
   public MEF90HeatXferOperator
   public MEF90HeatXferBilinearForm
   public MEF90HeatXferEnergy
   public MEF90HeatXferSetTransients
   public MEF90HeatXferIFunction
   public MEF90HeatXferIJacobian
   public MEF90HeatXferViewEXO
   public MEF90HeatXferCreateSNES
   public MEF90HeatXferCreateTS
contains

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferSetTransients"
!!!
!!!
!!!  MEF90HeatXferSetTransients: Update all transient data (boundary / external temperature and fluxes)
!!!                              using the proper scaling law
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90HeatXferSetTransients(MEF90HeatXferCtx, step, time, ierr)
      type(MEF90HeatXfer_Type), intent(INOUT)          :: MEF90HeatXferCtx
      PetscInt, intent(IN)                             :: step
      PetscReal, intent(IN)                            :: time
      PetscErrorCode, intent(INOUT)                      :: ierr

      type(MEF90HeatXferGlobalOptions_Type)            :: MEF90HeatXferGlobalOptions
      type(MEF90CtxGlobalOptions_Type)                 :: MEF90GlobalOptions
      type(tDM)                                        :: dm
      type(tVec)                                       :: tmpVec
      PetscExodusIIInt                                 :: exoStep

      EXOStep = step
      PetscCall(MEF90CtxGlobalOptionsSetFromOptions(MEF90HeatXferCtx%MEF90Ctx%comm, trim(MEF90HeatXferCtx%MEF90Ctx%prefix), MEF90GlobalOptions, ierr))
      MEF90HeatXferGlobalOptions = MEF90HeatXferCtx%globalOptions

      PetscCall(VecGetDM(MEF90HeatXferCtx%temperatureLocal, dm, ierr))

      select case (MEF90HeatXferGlobalOptions%boundaryTemperatureScaling)
      case (MEF90Scaling_File)
         PetscCall(DMGetLocalVector(dm, tmpVec, ierr))
         PetscCall(PetscObjectSetName(tmpVec, "Temperature", ierr))
         PetscCall(MEF90EXOVecLoad(tmpVec, MEF90HeatXferCtx%temperatureToIOSF, MEF90HeatXferCtx%IOToTemperatureSF, MEF90HeatXferCtx%MEF90Ctx%resultViewer, EXOstep, 1_ki, ierr))
         PetscCall(MEF90VecCopySF(tmpVec, MEF90HeatXferCtx%temperatureLocal, MEF90HeatXferCtx%boundaryToTemperatureSF, ierr))
         PetscCall(DMRestoreLocalVector(dm, tmpVec, ierr))
      case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90HeatXferCtx%temperatureLocal, time, ierr))
      case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90HeatXferCtx%temperatureLocal, 1.0_kr, ierr))
      case (MEF90Scaling_Expr)
         PetscCall(MEF90VecSetBCValuesFromOptionsExpr(MEF90HeatXferCtx%temperatureLocal, time, ierr))
      end select

      select case (MEF90HeatXferGlobalOptions%externalTemperatureScaling)
      case (MEF90Scaling_File)
         SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB, "External temperature from file not implemented yet "//__FUNCT__)
      case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90HeatXferCtx%externalTemperatureLocal, time, ierr))
      case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90HeatXferCtx%externalTemperatureLocal, 1.0_kr, ierr))
      case (MEF90Scaling_Expr)
         PetscCall(MEF90VecSetValuesFromOptionsExpr(MEF90HeatXferCtx%externalTemperatureLocal, time, ierr))
      end select

      select case (MEF90HeatXferGlobalOptions%fluxScaling)
      case (MEF90Scaling_File)
         SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB, "Flux from file not implemented yet "//__FUNCT__)
      case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90HeatXferCtx%fluxLocal, time, ierr))
      case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90HeatXferCtx%fluxLocal, 1.0_kr, ierr))
      case (MEF90Scaling_Expr)
         PetscCall(MEF90VecSetValuesFromOptionsExpr(MEF90HeatXferCtx%fluxLocal, time, ierr))
      end select

      select case (MEF90HeatXferGlobalOptions%boundaryFluxScaling)
      case (MEF90Scaling_File)
         PetscCall(MEF90EXOVecLoad(MEF90HeatXferCtx%boundaryFluxLocal, MEF90HeatXferCtx%boundaryFluxToIOSF, MEF90HeatXferCtx%IOToBoundaryFluxSF, MEF90HeatXferCtx%MEF90Ctx%resultViewer, EXOstep, 1_ki, ierr))
      case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90HeatXferCtx%boundaryFluxLocal, time, ierr))
      case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90HeatXferCtx%boundaryFluxLocal, 1.0_kr, ierr))
      case (MEF90Scaling_Expr)
         PetscCall(MEF90VecSetValuesFromOptionsExpr(MEF90HeatXferCtx%boundaryFluxLocal, time, ierr))
      end select
   end subroutine MEF90HeatXferSetTransients

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferOperator"
!!!
!!!
!!!  MEF90HeatXferOperator: wraps calls to MEF90HeatXferOperator from m_MEF90_HeatXferAssembly
!!!                         since overloading cannot be used here
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90HeatXferOperator(snesTemp, x, residual, PETScCtx, ierr)
      type(tSNES), intent(IN)                             :: snesTemp
      type(tVec), intent(IN)                              :: x
      type(tVec), intent(INOUT)                           :: residual
      type(c_ptr), intent(IN)                             :: PETScCtx
      type(MEF90HeatXfer_Type), pointer                   :: MEF90HeatXferCtx
      PetscErrorCode, intent(INOUT)                       :: ierr

      call c_f_pointer(PETScCtx, MEF90HeatXferCtx)
      if (MEF90HeatXferCtx%dim == 2) then
         PetscCall(MEF90HeatXferOperator2D(snesTemp, x, residual, MEF90HeatXferCtx, ierr))
      else if (MEF90HeatXferCtx%dim == 3) then
         PetscCall(MEF90HeatXferOperator3D(snesTemp, x, residual, MEF90HeatXferCtx, ierr))
      end if
   end subroutine MEF90HeatXferOperator

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferBilinearForm"
!!!
!!!
!!!  MEF90HeatXferBilinearForm: wraps calls to MEF90HeatXferBilinearForm from m_MEF90_HeatXferAssembly
!!!                             since overloading cannot be used here
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90HeatXferBilinearForm(snesTemp, x, A, M, PETScCtx, ierr)
      type(tSNES), intent(IN)                             :: snesTemp
      type(tVec), intent(IN)                              :: x
      type(tMat), intent(INOUT)                           :: A, M
      type(c_ptr), intent(IN)                             :: PETScCtx
      type(MEF90HeatXfer_Type), pointer                   :: MEF90HeatXferCtx
      PetscErrorCode, intent(INOUT)                       :: ierr

      call c_f_pointer(PETScCtx, MEF90HeatXferCtx)
      if (MEF90HeatXferCtx%dim == 2) then
         PetscCall(MEF90HeatXferBilinearForm2D(snesTemp, x, A, M, MEF90HeatXferCtx, ierr))
      else if (MEF90HeatXferCtx%dim == 3) then
         PetscCall(MEF90HeatXferBilinearForm3D(snesTemp, x, A, M, MEF90HeatXferCtx, ierr))
      end if
   end subroutine MEF90HeatXferBilinearForm

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXFerEnergy"
!!!
!!!
!!!  MEF90HeatXFerEnergy: wraps calls to MEF90HeatXferEnergy from m_MEF90_HeatXferAssembly
!!!                       since overloading cannot be used here
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90HeatXFerEnergy(MEF90HeatXferCtx, energy, bodyWork, surfaceWork, ierr)
      type(MEF90HeatXfer_Type), intent(IN)                 :: MEF90HeatXferCtx
      PetscReal, dimension(:), pointer                     :: energy, bodyWork, surfaceWork
      PetscErrorCode, intent(INOUT)                       :: ierr

      if (MEF90HeatXferCtx%dim == 2) then
         PetscCall(MEF90HeatXFerEnergy2D(MEF90HeatXferCtx, energy, bodyWork, surfaceWork, ierr))
      else if (MEF90HeatXferCtx%dim == 3) then
         PetscCall(MEF90HeatXFerEnergy3D(MEF90HeatXferCtx, energy, bodyWork, surfaceWork, ierr))
      end if
   end subroutine MEF90HeatXFerEnergy

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXFerIFunction"
!!!
!!!
!!!  MEF90HeatXFerIFunction: wraps calls to MEF90HeatXFerIFunction from m_MEF90_HeatXferAssembly
!!!                         since overloading cannot be used here
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90HeatXFerIFunction(tempTS, time, x, xdot, F, PETScCtx, ierr)
      type(tTS), intent(IN)                               :: tempTS
      PetscReal, intent(IN)                               :: time
      type(tVec), intent(IN)                              :: x, xdot
      type(tVec), intent(INOUT)                           :: F
      type(c_ptr), intent(IN)                             :: PETScCtx
      type(MEF90HeatXfer_Type), pointer                   :: MEF90HeatXferCtx
      PetscErrorCode, intent(INOUT)                       :: ierr

      call c_f_pointer(PETScCtx, MEF90HeatXferCtx)
      if (MEF90HeatXferCtx%dim == 2) then
         PetscCall(MEF90HeatXFerIFunction2D(tempTS, time, x, xdot, F, MEF90HeatXferCtx, ierr))
      else if (MEF90HeatXferCtx%dim == 3) then
         PetscCall(MEF90HeatXFerIFunction3D(tempTS, time, x, xdot, F, MEF90HeatXferCtx, ierr))
      end if
   end subroutine MEF90HeatXFerIFunction

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferIJacobian"
!!!
!!!
!!!  MEF90HeatXferIJacobian: wraps calls to MEF90HeatXferIJacobian from m_MEF90_HeatXferAssembly
!!!                         since overloading cannot be used here
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90HeatXferIJacobian(tempTS, t, x, xdot, shift, A, M, PETScCtx, ierr)
      type(tTS), intent(IN)                               :: tempTS
      PetscReal, intent(IN)                               :: t
      type(tVec), intent(IN)                              :: x, xdot
      PetscReal, intent(IN)                               :: shift
      type(tMat), intent(INOUT)                           :: A, M
      type(c_ptr), intent(IN)                             :: PETScCtx
      type(MEF90HeatXfer_Type), pointer                   :: MEF90HeatXferCtx
      PetscErrorCode, intent(INOUT)                       :: ierr

      call c_f_pointer(PETScCtx, MEF90HeatXferCtx)
      if (MEF90HeatXferCtx%dim == 2) then
         PetscCall(MEF90HeatXferIJacobian2D(tempTS, t, x, xdot, shift, A, M, MEF90HeatXferCtx, ierr))
      else if (MEF90HeatXferCtx%dim == 3) then
         PetscCall(MEF90HeatXferIJacobian3D(tempTS, t, x, xdot, shift, A, M, MEF90HeatXferCtx, ierr))
      end if
   end subroutine MEF90HeatXferIJacobian

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferViewEXO"
!!!
!!!
!!!  MEF90HeatXferViewEXO:
!!!
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022 Blaise Bourdin bourdin@mcmaster.ca
!!!      2022 Alexis Marboeuf marboeua@mcmaster.ca
!!!

   subroutine MEF90HeatXferViewEXO(MEF90HeatXferCtx, step, ierr)
      type(MEF90HeatXfer_Type), intent(IN)                :: MEF90HeatXferCtx
      PetscExodusIIInt, intent(IN)                        :: step
      PetscErrorCode, intent(INOUT)                         :: ierr

      type(MEF90HeatXferGlobalOptions_Type)               :: MEF90HeatXferGlobalOptions

      MEF90HeatXferGlobalOptions = MEF90HeatXferCtx%globalOptions

      if (MEF90HeatXferGlobalOptions%temperatureExport) then
         PetscCall(MEF90EXOVecView(MEF90HeatXferCtx%temperatureLocal, MEF90HeatXferCtx%temperatureToIOSF, MEF90HeatXferCtx%IOToTemperatureSF, MEF90HeatXferCtx%MEF90Ctx%resultViewer, step, 1_Ki, ierr))
      end if
   end subroutine MEF90HeatXferViewEXO

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferCreateSNES"
!!!
!!!
!!!  MEF90HeatXferCreateSNES:
!!!
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022 Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90HeatXferCreateSNES(MEF90HeatXferCtx, snesTemp, residual, ierr)
      type(MEF90HeatXfer_Type), target, intent(IN)        :: MEF90HeatXferCtx
      type(tSNES), intent(OUT)                            :: snesTemp
      type(tVec), intent(IN)                              :: residual
      PetscErrorCode, intent(INOUT)                         :: ierr

      type(MEF90HeatXferGlobalOptions_Type)              :: MEF90HeatXferGlobalOptions
      type(tDM)                                          :: dm
      type(tMat)                                         :: matTemp
      type(tMatNullSpace)                                :: nspTemp
      type(tKSP)                                         :: kspTemp
      PetscReal                                          :: rtol, dtol

      MEF90HeatXferGlobalOptions = MEF90HeatXferCtx%globalOptions
      PetscCall(VecGetDM(MEF90HeatXferCtx%temperatureLocal, dm, ierr))
      PetscCall(DMCreateMatrix(dm, matTemp, ierr))
      PetscCall(MatSetOptionsPrefix(matTemp, "Temperature_", ierr))
      !!! The matrix is not symmetric if the advection vector is /= 0
      PetscCall(MatSetOption(matTemp, MAT_SPD, PETSC_TRUE, ierr))
      PetscCall(MatSetOption(matTemp, MAT_SYMMETRY_ETERNAL, PETSC_TRUE, ierr))
      PetscCall(MatSetOption(matTemp, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE, ierr))
      if (MEF90HeatXferGlobalOptions%addNullSpace) then
         PetscCall(MatNullSpaceCreate(MEF90HeatXferCtx%MEF90Ctx%Comm, PETSC_TRUE, 0_ki, PETSC_NULL_VEC_ARRAY, nspTemp, ierr))
         PetscCall(MatSetNullSpace(matTemp, nspTemp, ierr))
         PetscCall(MatNullSpaceDestroy(nspTemp, ierr))
      end if
      PetscCall(MatSetFromOptions(matTemp, ierr))

      PetscCall(SNESCreate(MEF90HeatXferCtx%MEF90Ctx%Comm, snesTemp, ierr))
      PetscCall(SNESSetApplicationContext(snesTemp, MEF90HeatXferCtx%PETScCtx, ierr))
      PetscCall(SNESSetDM(snesTemp, dm, ierr))
      PetscCall(SNESSetType(snesTemp, SNESKSPONLY, ierr))
      PetscCall(SNESSetOptionsPrefix(snesTemp, 'Temperature_', ierr))

      PetscCall(SNESSetFunction(snesTemp, residual, MEF90HeatXferOperator, MEF90HeatXferCtx%PETScCtx, ierr))
      PetscCall(SNESSetJacobian(snesTemp, matTemp, matTemp, MEF90HeatXferBilinearForm, MEF90HeatXferCtx%PETScCtx, ierr))
      PetscCall(SNESSetFromOptions(snesTemp, ierr))
      !!!
      !!! Set some KSP options
      !!!
      PetscCall(SNESGetKSP(snesTemp, kspTemp, ierr))
      PetscCall(KSPSetType(kspTemp, KSPCG, ierr))
      PetscCall(KSPSetInitialGuessNonzero(kspTemp, PETSC_TRUE, ierr))
      rtol = 1.0d-8
      dtol = 1.0d+10
      PetscCall(KSPSetTolerances(kspTemp, rtol, PETSC_DEFAULT_REAL, dtol, PETSC_DEFAULT_INTEGER, ierr))
      PetscCall(KSPSetFromOptions(kspTemp, ierr))
      PetscCall(MatDestroy(matTemp, ierr))
   end subroutine MEF90HeatXferCreateSNES

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferCreateTS"
!!!
!!!
!!!  MEF90HeatXferCreateTS:
!!!
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90HeatXferCreateTS(MEF90HeatXferCtx, tsTemp, residual, initialTime, initialStep, ierr)

      type(MEF90HeatXfer_Type), target, intent(IN)        :: MEF90HeatXferCtx
      type(tTS), intent(OUT)                              :: tsTemp
      type(tVec), intent(IN)                              :: residual
      PetscReal, intent(IN)                               :: initialTime, initialStep
      PetscErrorCode, intent(INOUT)                         :: ierr

      type(MEF90HeatXferGlobalOptions_Type)              :: MEF90HeatXferGlobalOptions
      type(tDM)                                          :: dm
      type(tMat)                                         :: matTemp
      type(tMatNullSpace)                                :: nspTemp
      type(tSNES)                                        :: snesTemp
      type(tKSP)                                         :: kspTemp
      PetscReal                                          :: rtol, dtol

      MEF90HeatXferGlobalOptions = MEF90HeatXferCtx%globalOptions
      PetscCall(VecGetDM(MEF90HeatXferCtx%temperatureLocal, dm, ierr))
      PetscCall(DMCreateMatrix(dm, matTemp, ierr))
      PetscCall(MatSetOptionsPrefix(matTemp, "Temperature_", ierr))
      PetscCall(MatSetOption(matTemp, MAT_SPD, PETSC_TRUE, ierr))
      PetscCall(MatSetOption(matTemp, MAT_SYMMETRY_ETERNAL, PETSC_TRUE, ierr))
      PetscCall(MatSetOption(matTemp, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE, ierr))
      if (MEF90HeatXferGlobalOptions%addNullSpace) then
         PetscCall(MatNullSpaceCreate(MEF90HeatXferCtx%MEF90Ctx%Comm, PETSC_TRUE, 0_ki, PETSC_NULL_VEC_ARRAY, nspTemp, ierr))
         PetscCall(MatSetNullSpace(matTemp, nspTemp, ierr))
         PetscCall(MatNullSpaceDestroy(nspTemp, ierr))
      end if
      PetscCall(MatSetFromOptions(matTemp, ierr))

      PetscCall(TSCreate(MEF90HeatXferCtx%MEF90Ctx%Comm, tsTemp, ierr))
      PetscCall(TSSetDM(tsTemp, dm, ierr))
      PetscCall(TSSetOptionsPrefix(tsTemp, 'Temperature_', ierr))
      PetscCall(TSGetSNES(tsTemp, snesTemp, ierr))

      PetscCall(TSSetIFunction(tsTemp, residual, MEF90HeatXFerIFunction, MEF90HeatXferCtx%PETScCtx, ierr))
      PetscCall(TSSetIJacobian(tsTemp, matTemp, matTemp, MEF90HeatXFerIJacobian, MEF90HeatXferCtx%PETScCtx, ierr))

      PetscCall(TSSetType(tsTemp, 'rosw', ierr))
      PetscCall(TSRosWSetType(tsTemp, 'ra3pw', ierr))
      PetscCall(TSSetProblemType(tsTemp, TS_LINEAR, ierr))
      PetscCall(VecSet(MEF90HeatXferCtx%temperatureLocal, MEF90HeatXferGlobalOptions%initialTemperature, ierr))
      PetscCall(TSSetSolution(tsTemp, MEF90HeatXferCtx%temperatureLocal, ierr))
      PetscCall(TSSetTime(tsTemp, initialTime, ierr))
      PetscCall(TSSetTimeStep(tsTemp, initialStep, ierr))

      PetscCall(TSSetExactFinalTime(tsTemp, TS_EXACTFINALTIME_MATCHSTEP, ierr))
      PetscCall(TSSetFromOptions(tsTemp, ierr))
      !!!
      !!! Set some KSP options
      !!!
      PetscCall(SNESGetKSP(snesTemp, kspTemp, ierr))
      PetscCall(KSPSetType(kspTemp, KSPCG, ierr))
      PetscCall(KSPSetInitialGuessNonzero(kspTemp, PETSC_TRUE, ierr))
      rtol = 1.0d-8
      dtol = 1.0d+10
      PetscCall(KSPSetTolerances(kspTemp, rtol, PETSC_DEFAULT_REAL, dtol, PETSC_DEFAULT_INTEGER, ierr))
      PetscCall(KSPSetFromOptions(kspTemp, ierr))
      PetscCall(MatDestroy(matTemp, ierr))
   end subroutine MEF90HeatXferCreateTS
end module m_MEF90_HeatXfer
