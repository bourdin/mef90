#include "../MEF90/mef90.inc"
module m_MEF90_DefMech
#include "petsc/finclude/petsc.h"
   use petscsnes
   use petsctao

   use m_MEF90_EXO
   use m_MEF90_DefMechCtx

   use m_MEF90_DefMechAssembly2D, &
      MEF90DefMechOperatorDisplacement2D => MEF90DefMechOperatorDisplacement, &
      MEF90DefMechBilinearFormDisplacement2D => MEF90DefMechBilinearFormDisplacement, &
      MEF90DefMechWork2D => MEF90DefMechWork, &
      MEF90DefMechCohesiveEnergy2D => MEF90DefMechCohesiveEnergy, &
      MEF90DefMechPlasticdissipation2D => MEF90DefMechPlasticDissipation, &
      MEF90DefMechElasticEnergy2D => MEF90DefMechElasticEnergy, &
      MEF90DefMechOperatorDamage2D => MEF90DefMechOperatorDamage, &
      MEF90DefMechTAOGradientDamage2D => MEF90DefMechTAOGradientDamage, &
      MEF90DefMechBilinearFormDamage2D => MEF90DefMechBilinearFormDamage, &
      MEF90DefMechTAOHessianDamage2D => MEF90DefMechTAOHessianDamage, &
      MEF90DefMechSurfaceEnergy2D => MEF90DefMechSurfaceEnergy, &
      MEF90DefMechTAOObjectiveDamage2D => MEF90DefMechTAOObjectiveDamage, &
      MEF90DefMechCrackVolume2D => MEF90DefMechCrackVolume, &
      MEF90DefMechStress2D => MEF90DefMechStress
   use m_MEF90_DefMechPlasticity2D, &
      MEF90DefMechPlasticStrainUpdate2D => MEF90DefMechPlasticStrainUpdate
   use m_MEF90_DefMechAssembly3D, &
      MEF90DefMechOperatorDisplacement3D => MEF90DefMechOperatorDisplacement, &
      MEF90DefMechBilinearFormDisplacement3D => MEF90DefMechBilinearFormDisplacement, &
      MEF90DefMechWork3D => MEF90DefMechWork, &
      MEF90DefMechCohesiveEnergy3D => MEF90DefMechCohesiveEnergy, &
      MEF90DefMechPlasticdissipation3D => MEF90DefMechPlasticDissipation, &
      MEF90DefMechElasticEnergy3D => MEF90DefMechElasticEnergy, &
      MEF90DefMechOperatorDamage3D => MEF90DefMechOperatorDamage, &
      MEF90DefMechTAOGradientDamage3D => MEF90DefMechTAOGradientDamage, &
      MEF90DefMechBilinearFormDamage3D => MEF90DefMechBilinearFormDamage, &
      MEF90DefMechTAOHessianDamage3D => MEF90DefMechTAOHessianDamage, &
      MEF90DefMechSurfaceEnergy3D => MEF90DefMechSurfaceEnergy, &
      MEF90DefMechTAOObjectiveDamage3D => MEF90DefMechTAOObjectiveDamage, &
      MEF90DefMechCrackVolume3D => MEF90DefMechCrackVolume, &
      MEF90DefMechStress3D => MEF90DefMechStress
   use m_MEF90_DefMechPlasticity3D, &
      MEF90DefMechPlasticStrainUpdate3D => MEF90DefMechPlasticStrainUpdate

   implicit none(type)
   ! private
   public :: MEF90DefMechSetTransients
   public :: MEF90DefMechOperatorDisplacement
   public :: MEF90DefMechBilinearFormDisplacement
   public :: MEF90DefMechCreateSNESDisplacement

   public :: MEF90DefMechOperatorDamage
   public :: MEF90DefMechBilinearFormDamage
   public :: MEF90DefMechCreateSNESDamage
   public :: MEF90DefMechUpdateDamageBounds

   public :: MEF90DefMechTAOObjectiveDamage
   public :: MEF90DefMechTAOGradientDamage
   public :: MEF90DefMechTAOHessianDamage
   public :: MEF90DefMechCreateTAODamage
   public :: MEF90DefMechTAOUpdateDamageBounds

   public :: MEF90DefMechViewEXO
   public :: MEF90DefMechSurfaceEnergy
   public :: MEF90DefMechElasticEnergy
   public :: MEF90DefMechWork
   public :: MEF90DefMechCohesiveEnergy
   public :: MEF90DefMechPlasticDissipation
   public :: MEF90DefMechCrackVolume
   public :: MEF90DefMechStress
   public :: MEF90DefMechPlasticStrainUpdate

   public :: MEF90DefMechFormatEXO

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSetTransients"
!!!
!!!
!!!  MEF90DefMechSetTransients:
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!  (c)    2022 Alexis Marboeuf marboeua@mcmaster.ca
!!!

   subroutine MEF90DefMechSetTransients(MEF90DefMechCtx, step, time, ierr)
      type(MEF90DefMechCtx_Type), intent(INOUT)        :: MEF90DefMechCtx
      PetscInt, intent(IN)                             :: step
      PetscReal, intent(IN)                            :: time
      PetscErrorCode, intent(INOUT)                    :: ierr

      type(MEF90DefMechGlobalOptions_Type), pointer    :: MEF90DefMechGlobalOptions
      type(MEF90CtxGlobalOptions_Type), pointer        :: MEF90GlobalOptions
      type(tDM)                                       :: dmDisplacement, dmDamage, dmCohesiveDisplacement
      type(tVec)                                      :: tmpVec
      character(len=MEF90MXSTRLEN)                    :: IOBuffer

      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag, MEF90GlobalOptions, ierr))
      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag, MEF90DefMechGlobalOptions, ierr))

      PetscCall(VecGetDM(MEF90DefMechCtx%damageLocal, dmDamage, ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%displacementLocal, dmDisplacement, ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%cohesiveDisplacement, dmCohesiveDisplacement, ierr))

      select case (MEF90DefMechGlobalOptions%boundaryDisplacementScaling)
      case (MEF90Scaling_File)
         PetscCall(DMGetLocalVector(dmDisplacement, tmpVec, ierr))
         PetscCall(PetscObjectSetName(tmpVec, "Displacement", ierr))
         PetscCall(MEF90EXOVecLoad(tmpVec, MEF90DefMechCtx%displacementToIOSF, MEF90DefMechCtx%IOToDisplacementSF, MEF90DefMechCtx%MEF90Ctx%resultViewer, step, MEF90DefMechCtx%dim, ierr))
         PetscCall(MEF90VecCopySF(tmpVec, MEF90DefMechCtx%displacementLocal, MEF90DefMechCtx%displacementConstraintsSF, ierr))
         PetscCall(DMRestoreLocalVector(dmDisplacement, tmpVec, ierr))
      case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%displacementLocal, time, ierr))
      case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%displacementLocal, 1.0_kr, ierr))
      case (MEF90Scaling_Expr)
         PetscCall(MEF90VecSetBCValuesFromOptionsExpr(MEF90DefMechCtx%displacementLocal, time, ierr))
      end select

      select case (MEF90DefMechGlobalOptions%boundaryDamageScaling)
      case (MEF90Scaling_File)
         PetscCall(DMGetLocalVector(dmDamage, tmpVec, ierr))
         PetscCall(PetscObjectSetName(tmpVec, "Damage", ierr))
         PetscCall(MEF90EXOVecLoad(tmpVec, MEF90DefMechCtx%damageToIOSF, MEF90DefMechCtx%IOToDamageSF, MEF90DefMechCtx%MEF90Ctx%resultViewer, step, 1_ki, ierr))
         PetscCall(MEF90VecCopySF(tmpVec, MEF90DefMechCtx%damageLocal, MEF90DefMechCtx%damageConstraintsSF, ierr))
         PetscCall(DMRestoreLocalVector(dmDamage, tmpVec, ierr))
      case (MEF90Scaling_Linear)
         write (IOBuffer, '((A),": linear scaling of damage does not make any sense.\n")') __FUNCT__
         SETERRQ(MEF90DefMechCtx%MEF90Ctx%Comm, PETSC_ERR_ARG_WRONG, IOBuffer)
      case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%damageLocal, 1.0_kr, ierr))
      case (MEF90Scaling_Expr)
         PetscCall(MEF90VecSetBCValuesFromOptionsExpr(MEF90DefMechCtx%damageLocal, time, ierr))
      end select

      select case (MEF90DefMechGlobalOptions%cohesiveDisplacementScaling)
      case (MEF90Scaling_File)
         PetscCall(MEF90EXOVecLoad(MEF90DefMechCtx%cohesiveDisplacement, MEF90DefMechCtx%cohesiveDisplacementToIOSF, MEF90DefMechCtx%IOToCohesiveDisplacementSF, MEF90DefMechCtx%MEF90Ctx%resultViewer, step, MEF90DefMechCtx%dim, ierr))
      case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90DefMechCtx%cohesiveDisplacement, time, ierr))
      case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90DefMechCtx%cohesiveDisplacement, 1.0_kr, ierr))
      case (MEF90Scaling_Expr)
         PetscCall(MEF90VecSetValuesFromOptionsExpr(MEF90DefMechCtx%cohesiveDisplacement, time, ierr))
      end select

      select case (MEF90DefMechGlobalOptions%displacementLowerBoundScaling)
      case (MEF90Scaling_File)
         write (*, *) __FUNCT__, ": file scaling of displacement lower bound does not make any sense."
         stop
      case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%displacementLowerBoundLocal, time, ierr))
      case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%displacementLowerBoundLocal, 1.0_kr, ierr))
      case (MEF90Scaling_Expr)
         PetscCall(MEF90VecSetBCValuesFromOptionsExpr(MEF90DefMechCtx%displacementLowerBoundLocal, time, ierr))
      end select

      select case (MEF90DefMechGlobalOptions%displacementUpperBoundScaling)
      case (MEF90Scaling_File)
         write (*, *) __FUNCT__, ": file scaling of displacement upper bound does not make any sense."
         stop
      case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%displacementUpperBoundLocal, time, ierr))
      case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%displacementUpperBoundLocal, 1.0_kr, ierr))
      case (MEF90Scaling_Expr)
         PetscCall(MEF90VecSetBCValuesFromOptionsExpr(MEF90DefMechCtx%displacementUpperBoundLocal, time, ierr))
      end select

      select case (MEF90DefMechGlobalOptions%bodyForceScaling)
      case (MEF90Scaling_File)
         PetscCall(MEF90EXOVecLoad(MEF90DefMechCtx%bodyForce, MEF90DefMechCtx%bodyForceToIOSF, MEF90DefMechCtx%IOToBodyForceSF, MEF90DefMechCtx%MEF90Ctx%resultViewer, step, MEF90DefMechCtx%dim, ierr))
      case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90DefMechCtx%bodyForce, time, ierr))
      case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90DefMechCtx%bodyForce, 1.0_kr, ierr))
      case (MEF90Scaling_Expr)
         PetscCall(MEF90VecSetValuesFromOptionsExpr(MEF90DefMechCtx%bodyForce, time, ierr))
      end select

      select case (MEF90DefMechGlobalOptions%boundaryForceScaling)
      case (MEF90Scaling_File)
         SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB, "Boundary force from file not implemented yet "//__FUNCT__)
         ! PetscCall(MEF90EXOVecLoad(MEF90DefMechCtx%boundaryForce,MEF90DefMechCtx%boundaryForceToIOSF,MEF90DefMechCtx%IOToBoundaryForceSF,MEF90DefMechCtx%MEF90Ctx%resultViewer,step,MEF90DefMechCtx%dim,ierr))
      case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90DefMechCtx%boundaryForce, time, ierr))
      case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90DefMechCtx%boundaryForce, 1.0_kr, ierr))
      case (MEF90Scaling_Expr)
         PetscCall(MEF90VecSetValuesFromOptionsExpr(MEF90DefMechCtx%boundaryForce, time, ierr))
      end select

      select case (MEF90DefMechGlobalOptions%pressureForceScaling)
      case (MEF90Scaling_File)
         SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB, "Pressure force from file not implemented yet "//__FUNCT__)
         ! PetscCall(MEF90EXOVecLoad(MEF90DefMechCtx%pressureForce,MEF90DefMechCtx%pressureForceToIOSF,MEF90DefMechCtx%IOToPressureForceSF,MEF90DefMechCtx%MEF90Ctx%resultViewer,step,1_Ki,ierr))
      case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90DefMechCtx%pressureForce, time, ierr))
      case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90DefMechCtx%pressureForce, 1.0_kr, ierr))
      case (MEF90Scaling_Expr)
         PetscCall(MEF90VecSetValuesFromOptionsExpr(MEF90DefMechCtx%pressureForce, time, ierr))
      end select
   end subroutine MEF90DefMechSetTransients

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDisplacement"
!!!
!!!
!!!  MEF90DefMechOperatorDisplacement: wraps calls to MEF90DefMechOperatorDisplacement from m_MEF90_DefMechAssembly
!!!                        since overloading cannot be used here
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   subroutine MEF90DefMechOperatorDisplacement(snesTemp, x, residual, MEF90DefMechCtx, ierr)
      type(tSNES), intent(IN)                             :: snesTemp
      type(tVec), intent(IN)                              :: x
      type(tVec), intent(INOUT)                           :: residual
      type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode, intent(INOUT)                       :: ierr

      if (MEF90DefMechCtx%dim == 2) then
         PetscCall(MEF90DefMechOperatorDisplacement2D(snesTemp, x, residual, MEF90DefMechCtx, ierr))
      else if (MEF90DefMechCtx%dim == 3) then
         PetscCall(MEF90DefMechOperatorDisplacement3D(snesTemp, x, residual, MEF90DefMechCtx, ierr))
      end if
   end subroutine MEF90DefMechOperatorDisplacement

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDisplacement"
!!!
!!!
!!!  MEF90DefMechBilinearFormDisplacement: wraps calls to MEF90DefMechBilinearFormDisplacement from m_MEF90_DefMechAssembly
!!!                            since overloading cannot be used here
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   subroutine MEF90DefMechBilinearFormDisplacement(snesDispl, x, A, M, MEF90DefMechCtx, ierr)
      type(tSNES), intent(IN)                             :: snesDispl
      type(tVec), intent(IN)                              :: x
      type(tMat), intent(INOUT)                           :: A, M
      type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode, intent(INOUT)                       :: ierr

      if (MEF90DefMechCtx%dim == 2) then
         PetscCall(MEF90DefMechBilinearFormDisplacement2D(snesDispl, x, A, M, MEF90DefMechCtx, ierr))
      else if (MEF90DefMechCtx%dim == 3) then
         PetscCall(MEF90DefMechBilinearFormDisplacement3D(snesDispl, x, A, M, MEF90DefMechCtx, ierr))
      end if
   end subroutine MEF90DefMechBilinearFormDisplacement

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechWork"
!!!
!!!
!!!  MEF90DefMechWork: wraps calls to MEF90DefMechWork from m_MEF90_DefMechAssembly
!!                       since overloading cannot be used here
!!!
!!!  (c) 2012-22 Blaise Bourdin bourdin@lsu.edu
!!!

   subroutine MEF90DefMechWork(MEF90DefMechCtx, bodyForceWork, boundaryForceWork, ierr)
      type(MEF90DefMechCtx_Type), intent(IN)           :: MEF90DefMechCtx
      PetscReal, dimension(:), pointer                  :: bodyForceWork, boundaryForceWork
      PetscErrorCode, intent(INOUT)                    :: ierr

      if (MEF90DefMechCtx%dim == 2) then
         PetscCall(MEF90DefMechWork2D(MEF90DefMechCtx, bodyForceWork, boundaryForceWork, ierr))
      else if (MEF90DefMechCtx%dim == 3) then
         PetscCall(MEF90DefMechWork3D(MEF90DefMechCtx, bodyForceWork, boundaryForceWork, ierr))
      end if
   end subroutine MEF90DefMechWork

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCohesiveEnergy"
!!!
!!!
!!!  MEF90DefMechCohesiveEnergy: wraps calls to MEF90DefMechCohesiveEnergy from m_MEF90_DefMechAssembly
!!                       since overloading cannot be used here
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   subroutine MEF90DefMechCohesiveEnergy(MEF90DefMechCtx, cohesiveEnergy, ierr)
      type(MEF90DefMechCtx_Type), intent(IN)           :: MEF90DefMechCtx
      PetscReal, dimension(:), pointer                  :: cohesiveEnergy
      PetscErrorCode, intent(INOUT)                    :: ierr

      if (MEF90DefMechCtx%dim == 2) then
         PetscCall(MEF90DefMechCohesiveEnergy2D(MEF90DefMechCtx, cohesiveEnergy, ierr))
      else if (MEF90DefMechCtx%dim == 3) then
         PetscCall(MEF90DefMechCohesiveEnergy3D(MEF90DefMechCtx, cohesiveEnergy, ierr))
      end if
   end subroutine MEF90DefMechCohesiveEnergy

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechElasticEnergy"
!!!
!!!
!!!  MEF90DefMechElasticEnergy: wraps calls to MEF90DefMechElasticEnergy from m_MEF90_DefMechAssembly
!!                       since overloading cannot be used here
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   subroutine MEF90DefMechElasticEnergy(MEF90DefMechCtx, energy, ierr)
      type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
      PetscReal, dimension(:), pointer                     :: energy
      PetscErrorCode, intent(INOUT)                       :: ierr

      if (MEF90DefMechCtx%dim == 2) then
         PetscCall(MEF90DefMechElasticEnergy2D(MEF90DefMechCtx, energy, ierr))
      else if (MEF90DefMechCtx%dim == 3) then
         PetscCall(MEF90DefMechElasticEnergy3D(MEF90DefMechCtx, energy, ierr))
      end if
   end subroutine MEF90DefMechElasticEnergy

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechPlasticDissipation"
!!!
!!!
!!!  MEF90DefMechPlasticDissipation: wraps calls to MEF90DefMechPlasticDissipation from m_MEF90_DefMechAssembly
!!                       since overloading cannot be used here
!!!
!!!  (c) 2012-14 Erwan TANNE erwan.tanne@gmail.com
!!!

   subroutine MEF90DefMechPlasticDissipation(x, MEF90DefMechCtx, plasticStrainOld, energy, ierr)
      type(tVec), intent(IN)                              :: x
      type(tVec), intent(IN)                              :: plasticStrainOld
      type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
      PetscReal, dimension(:), pointer                     :: energy
      PetscErrorCode, intent(INOUT)                       :: ierr

      if (MEF90DefMechCtx%dim == 2) then
         ! Call MEF90DefMechPlasticDissipation2D(x,MEF90DefMechCtx,plasticStrainOld,energy,ierr)
      else if (MEF90DefMechCtx%dim == 3) then
         ! Call MEF90DefMechPlasticDissipation3D(x,MEF90DefMechCtx,plasticStrainOld,energy,ierr)
      end if
   end subroutine MEF90DefMechPlasticDissipation

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechStress"
!!!
!!!
!!!  MEF90DefMechStress: wraps calls to MEF90DefMechElasticEnergy from m_MEF90_DefMechAssembly
!!                       since overloading cannot be used here
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   subroutine MEF90DefMechStress(MEF90DefMechCtx, stress, ierr)
      type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
      type(tVec), intent(IN)                              :: stress
      PetscErrorCode, intent(INOUT)                       :: ierr

      if (MEF90DefMechCtx%dim == 2) then
         PetscCall(MEF90DefMechStress2D(MEF90DefMechCtx, stress, ierr))
      else if (MEF90DefMechCtx%dim == 3) then
         PetscCall(MEF90DefMechStress3D(MEF90DefMechCtx, stress, ierr))
      end if
   end subroutine MEF90DefMechStress

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCrackVolume"
!!!
!!!
!!!  MEF90DefMechCrackVolume:
!!!
!!!  (c) 2016 erwan
!!!

   subroutine MEF90DefMechCrackVolume(MEF90DefMechCtx, CrackVolume, ierr)
      type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
      PetscReal, dimension(:), pointer                     :: CrackVolume
      PetscErrorCode, intent(INOUT)                       :: ierr

      if (MEF90DefMechCtx%dim == 2) then
         PetscCall(MEF90DefMechCrackVolume2D(MEF90DefMechCtx, CrackVolume, ierr))
      else if (MEF90DefMechCtx%dim == 3) then
         PetscCall(MEF90DefMechCrackVolume3D(MEF90DefMechCtx, CrackVolume, ierr))
      end if
   end subroutine MEF90DefMechCrackVolume

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamage"
!!!
!!!
!!!  MEF90DefMechOperatorDamage: wraps calls to MEF90DefMechOperatorDamage from m_MEF90_DefMechAssembly
!!!                        since overloading cannot be used here
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   subroutine MEF90DefMechOperatorDamage(snesDamage, damage, residual, MEF90DefMechCtx, ierr)
      type(tSNES), intent(IN)                             :: snesDamage
      type(tVec), intent(IN)                              :: damage
      type(tVec), intent(INOUT)                           :: residual
      type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode, intent(INOUT)                       :: ierr

      if (MEF90DefMechCtx%dim == 2) then
         PetscCall(MEF90DefMechOperatorDamage2D(snesDamage, damage, residual, MEF90DefMechCtx, ierr))
      else if (MEF90DefMechCtx%dim == 3) then
         PetscCall(MEF90DefMechOperatorDamage3D(snesDamage, damage, residual, MEF90DefMechCtx, ierr))
      end if
   end subroutine MEF90DefMechOperatorDamage

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechTAOGradientDamage"
!!!
!!!
!!!  MEF90DefMechTAOGradientDamage:
!!!
!!!  (c) 2022 Alexis Marboeuf marboeua@mcmaster.ca
!!!

   subroutine MEF90DefMechTAOGradientDamage(taoDamage, damage, residual, MEF90DefMechCtx, ierr)
      type(tTao), intent(IN)                              :: taoDamage
      type(tVec), intent(IN)                              :: damage
      type(tVec), intent(INOUT)                           :: residual
      type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode, intent(INOUT)                       :: ierr

      if (MEF90DefMechCtx%dim == 2) then
         PetscCall(MEF90DefMechTAOGradientDamage2D(taoDamage, damage, residual, MEF90DefMechCtx, ierr))
      else if (MEF90DefMechCtx%dim == 3) then
         PetscCall(MEF90DefMechTAOGradientDamage3D(taoDamage, damage, residual, MEF90DefMechCtx, ierr))
      end if
   end subroutine MEF90DefMechTAOGradientDamage

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamage"
!!!
!!!
!!!  MEF90DefMechBilinearFormDamage: wraps calls to MEF90DefMechBilinearFormDamage from m_MEF90_DefMechAssembly
!!!                            since overloading cannot be used here
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   subroutine MEF90DefMechBilinearFormDamage(snesDamage, damage, A, M, MEF90DefMechCtx, ierr)
      type(tSNES), intent(IN)                             :: snesDamage
      type(tVec), intent(IN)                              :: damage
      type(tMat), intent(INOUT)                           :: A, M
      type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode, intent(INOUT)                       :: ierr

      if (MEF90DefMechCtx%dim == 2) then
         PetscCall(MEF90DefMechBilinearFormDamage2D(snesDamage, damage, A, M, MEF90DefMechCtx, ierr))
      else if (MEF90DefMechCtx%dim == 3) then
         PetscCall(MEF90DefMechBilinearFormDamage3D(snesDamage, damage, A, M, MEF90DefMechCtx, ierr))
      end if
   end subroutine MEF90DefMechBilinearFormDamage

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechTAOHessianDamage"
!!!
!!!
!!!  MEF90DefMechTAOHessianDamage:
!!!
!!!  (c) 2022 Alexis Marboeuf marboeua@mcmaster.ca
!!!

   subroutine MEF90DefMechTAOHessianDamage(taoDamage, damage, A, M, MEF90DefMechCtx, ierr)
      type(tTao), intent(IN)                              :: taoDamage
      type(tVec), intent(IN)                              :: damage
      type(tMat), intent(INOUT)                           :: A, M
      type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode, intent(INOUT)                       :: ierr

      if (MEF90DefMechCtx%dim == 2) then
         PetscCall(MEF90DefMechTAOHessianDamage2D(taoDamage, damage, A, M, MEF90DefMechCtx, ierr))
      else if (MEF90DefMechCtx%dim == 3) then
         PetscCall(MEF90DefMechTAOHessianDamage3D(taoDamage, damage, A, M, MEF90DefMechCtx, ierr))
      end if
   end subroutine MEF90DefMechTAOHessianDamage

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSurfaceEnergy"
!!!
!!!
!!!  MEF90DefMechSurfaceEnergy: wraps calls to MEF90DefMechSurfaceEnergy from m_MEF90_DefMechAssembly
!!                       since overloading cannot be used here
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   subroutine MEF90DefMechSurfaceEnergy(MEF90DefMechCtx, energy, ierr)
      type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
      PetscReal, dimension(:), pointer                     :: energy
      PetscErrorCode, intent(INOUT)                       :: ierr

      if (MEF90DefMechCtx%dim == 2) then
         PetscCall(MEF90DefMechSurfaceEnergy2D(MEF90DefMechCtx, energy, ierr))
      else if (MEF90DefMechCtx%dim == 3) then
         PetscCall(MEF90DefMechSurfaceEnergy3D(MEF90DefMechCtx, energy, ierr))
      end if
   end subroutine MEF90DefMechSurfaceEnergy

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechTAOObjectiveDamage"
!!!
!!!
!!!  MEF90DefMechTAOObjectiveDamage:
!!!
!!!  (c) 2022 Alexis Marboeuf marboeua@mcmaster.ca
!!!

   subroutine MEF90DefMechTAOObjectiveDamage(taoDamage, damage, energy, MEF90DefMechCtx, ierr)
      type(tTao), intent(IN)                              :: taoDamage
      type(tVec), intent(IN)                              :: damage
      PetscReal, intent(INOUT)                            :: energy
      type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode, intent(INOUT)                       :: ierr

      if (MEF90DefMechCtx%dim == 2) then
         PetscCall(MEF90DefMechTAOObjectiveDamage2D(taoDamage, damage, energy, MEF90DefMechCtx, ierr))
      else if (MEF90DefMechCtx%dim == 3) then
         PetscCall(MEF90DefMechTAOObjectiveDamage3D(taoDamage, damage, energy, MEF90DefMechCtx, ierr))
      end if
   end subroutine MEF90DefMechTAOObjectiveDamage

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechFormatEXO"
!!!
!!!
!!!  MEF90DefMechFormatEXO:
!!!
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022 Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90DefMechFormatEXO(MEF90DefMechCtx, time, ierr)
      type(MEF90DefMechCtx_Type), intent(INOUT)           :: MEF90DefMechCtx
      PetscReal, dimension(:), pointer                     :: time
      PetscErrorCode, intent(OUT)                         :: ierr

      character(len=MXSTLN), dimension(:), pointer         :: nameG, nameN, nameC
      type(MEF90DefMechGlobalOptions_Type), pointer       :: MEF90DefMechGlobalOptions
      PetscInt                                           :: numFields, offset

      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag, MEF90DefMechGlobalOptions, ierr))
      allocate (nameG(0))

      numFields = 0
      if (MEF90DefMechGlobalOptions%displacementExport) then
         numFields = numFields + MEF90DefMechCtx%dim
      end if
      if (MEF90DefMechGlobalOptions%damageExport) then
         numFields = numFields + 1
      end if
      if (MEF90DefMechGlobalOptions%temperatureExport) then
         numFields = numFields + 1
      end if

      allocate (nameN(numFields))
      offset = 1
      if (MEF90DefMechGlobalOptions%displacementExport) then
         nameN(offset + 0) = "Displacement_X"
         nameN(offset + 1) = "Displacement_Y"
         if (MEF90DefMechCtx%dim == 3) then
            nameN(offset + 2) = "Displacement_Z"
         end if
         offset = offset + MEF90DefMechCtx%dim
      end if
      if (MEF90DefMechGlobalOptions%damageExport) then
         nameN(offset) = "Damage"
         offset = offset + 1
      end if
      if (MEF90DefMechGlobalOptions%temperatureExport) then
         nameN(offset) = "Temperature"
         offset = offset + 1
      end if

      numFields = 0
      if (MEF90DefMechGlobalOptions%stressExport) then
         numFields = numFields + MEF90DefMechCtx%dim * (MEF90DefMechCtx%dim + 1) / 2
      end if
      if (MEF90DefMechGlobalOptions%plasticStrainExport) then
         numFields = numFields + MEF90DefMechCtx%dim * (MEF90DefMechCtx%dim + 1) / 2
      end if
      if (MEF90DefMechGlobalOptions%cumulatedPlasticDissipationExport) then
         numFields = numFields + MEF90DefMechCtx%dim * (MEF90DefMechCtx%dim + 1) / 2
      end if

      allocate (nameC(numFields))
      offset = 1
      if (MEF90DefMechGlobalOptions%stressExport) then
         if (MEF90DefMechCtx%dim == 2) then
            nameC(offset + 0) = "Stress_XX"
            nameC(offset + 1) = "Stress_YY"
            nameC(offset + 2) = "Stress_XY"
         else
            nameC(offset + 0) = "Stress_XX"
            nameC(offset + 1) = "Stress_YY"
            nameC(offset + 2) = "Stress_ZZ"
            nameC(offset + 3) = "Stress_YZ"
            nameC(offset + 4) = "Stress_XZ"
            nameC(offset + 5) = "Stress_XY"
         end if
         offset = offset + MEF90DefMechCtx%dim * (MEF90DefMechCtx%dim + 1) / 2
      end if

      if (MEF90DefMechGlobalOptions%plasticStrainExport) then
         if (MEF90DefMechCtx%dim == 2) then
            nameC(offset + 0) = "PlasticStrain_XX"
            nameC(offset + 1) = "PlasticStrain_YY"
            nameC(offset + 2) = "PlasticStrain_XY"
         else
            nameC(offset + 0) = "PlasticStrain_XX"
            nameC(offset + 1) = "PlasticStrain_YY"
            nameC(offset + 2) = "PlasticStrain_ZZ"
            nameC(offset + 3) = "PlasticStrain_YZ"
            nameC(offset + 4) = "PlasticStrain_XZ"
            nameC(offset + 5) = "PlasticStrain_XY"
         end if
         offset = offset + MEF90DefMechCtx%dim * (MEF90DefMechCtx%dim + 1) / 2
      end if

      if (MEF90DefMechGlobalOptions%cumulatedPlasticDissipationExport) then
         if (MEF90DefMechCtx%dim == 2) then
            nameC(offset + 0) = "CumulatedPlasticDissipation_XX"
            nameC(offset + 1) = "CumulatedPlasticDissipation_YY"
            nameC(offset + 2) = "CumulatedPlasticDissipation_XY"
         else
            nameC(offset + 0) = "CumulatedPlasticDissipation_XX"
            nameC(offset + 1) = "CumulatedPlasticDissipation_YY"
            nameC(offset + 2) = "CumulatedPlasticDissipation_ZZ"
            nameC(offset + 3) = "CumulatedPlasticDissipation_YZ"
            nameC(offset + 4) = "CumulatedPlasticDissipation_XZ"
            nameC(offset + 5) = "CumulatedPlasticDissipation_XY"
         end if
      end if
      PetscCall(MEF90EXOFormat(MEF90DefMechCtx%MEF90Ctx%resultViewer, nameG, nameC, nameN, time, ierr))
      deallocate (nameG)
      deallocate (nameN)
      deallocate (nameC)
   end subroutine MEF90DefMechFormatEXO

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechViewEXO"
!!!
!!!
!!!  MEF90DefMechViewEXO: Save all fields in a MEF90DefMechCtx_Type in an exodus file
!!!
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022 Alexis Marboeuf, marboeua@mcmaster.ca
!!!

   subroutine MEF90DefMechViewEXO(MEF90DefMechCtx, step, ierr)
      type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
      PetscInt, intent(IN)                                :: step
      PetscErrorCode, intent(INOUT)                       :: ierr

      type(MEF90DefMechGlobalOptions_Type), pointer       :: MEF90DefMechGlobalOptions

      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag, MEF90DefMechGlobalOptions, ierr))

      if (MEF90DefMechGlobalOptions%displacementExport) then
         PetscCall(MEF90EXOVecView(MEF90DefMechCtx%displacementLocal, MEF90DefMechCtx%displacementToIOSF, MEF90DefMechCtx%IOToDisplacementSF, MEF90DefMechCtx%MEF90Ctx%resultViewer, step, MEF90DefMechCtx%dim, ierr))
      end if
      if (MEF90DefMechGlobalOptions%damageExport) then
         PetscCall(MEF90EXOVecView(MEF90DefMechCtx%damageLocal, MEF90DefMechCtx%damageToIOSF, MEF90DefMechCtx%IOToDamageSF, MEF90DefMechCtx%MEF90Ctx%resultViewer, step, 1_ki, ierr))
      end if
      if (MEF90DefMechGlobalOptions%stressExport) then
         PetscCall(MEF90EXOVecView(MEF90DefMechCtx%stress, MEF90DefMechCtx%stressToIOSF, MEF90DefMechCtx%IOToStressSF, MEF90DefMechCtx%MEF90Ctx%resultViewer, step, MEF90DefMechCtx%dim * (MEF90DefMechCtx%dim + 1) / 2, ierr))
      end if
      if (MEF90DefMechGlobalOptions%plasticStrainExport) then
         PetscCall(MEF90EXOVecView(MEF90DefMechCtx%plasticStrain, MEF90DefMechCtx%plasticStrainToIOSF, MEF90DefMechCtx%IOToPlasticStrainSF, MEF90DefMechCtx%MEF90Ctx%resultViewer, step, MEF90DefMechCtx%dim * (MEF90DefMechCtx%dim + 1) / 2, ierr))
      end if
      if (MEF90DefMechGlobalOptions%cumulatedPlasticDissipationExport) then
         PetscCall(MEF90EXOVecView(MEF90DefMechCtx%cumulatedPlasticDissipation, MEF90DefMechCtx%cumulatedPlasticDissToIOSF, MEF90DefMechCtx%IOToCumulatedPlasticDissSF, MEF90DefMechCtx%MEF90Ctx%resultViewer, step, MEF90DefMechCtx%dim * (MEF90DefMechCtx%dim + 1) / 2, ierr))
      end if
   end subroutine MEF90DefMechViewEXO

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCreateSNESDisplacement"
!!!
!!!
!!!  MEF90DefMechCreateSNESDisplacement:
!!!
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022 Alexis Marboeuf, marboeua@mcmaster.ca
!!!

   subroutine MEF90DefMechCreateSNESDisplacement(MEF90DefMechCtx, snesDisplacement, residual, ierr)
      type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
      type(tSNES), intent(OUT)                            :: snesDisplacement
      type(tVec), intent(IN)                              :: residual
      PetscErrorCode, intent(INOUT)                       :: ierr

      type(MEF90DefMechGlobalOptions_Type), pointer       :: MEF90DefMechGlobalOptions
      type(tDM)                                          :: dm
      type(tMat)                                         :: matDisplacement
      type(tMatNullSpace)                                :: nspDisplacement
      type(tKSP)                                         :: kspDisplacement
      type(tVec)                                         :: gCoord
      PetscReal                                          :: rtol, dtol, atol, stol

      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag, MEF90DefMechGlobalOptions, ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%displacementLocal, dm, ierr))
      PetscCall(DMCreateMatrix(dm, matDisplacement, ierr))
      PetscCall(MatSetOptionsPrefix(matDisplacement, "Displacement_", ierr))
      PetscCall(MatSetOption(matDisplacement, MAT_SPD, PETSC_TRUE, ierr))
      PetscCall(MatSetOption(matDisplacement, MAT_SYMMETRY_ETERNAL, PETSC_TRUE, ierr))
      PetscCall(MatSetOption(matDisplacement, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE, ierr))
      !!!
      !!! Set the matrix near null-space consisting of all rigid motions.
      !!!
      PetscCall(DMGetGlobalVector(dm, gCoord, ierr))
      PetscCall(MEF90DefMechProjectCoordinates_Private(gCoord, ierr))
      PetscCall(MatNullSpaceCreateRigidBody(gCoord, nspDisplacement, ierr))
      PetscCall(MatSetNearNullSpace(matDisplacement, nspDisplacement, ierr))
      PetscCall(MatNullSpaceDestroy(nspDisplacement, ierr))
      PetscCall(DMRestoreGlobalVector(dm, gCoord, ierr))
      PetscCall(MatSetFromOptions(matDisplacement, ierr))

      PetscCall(SNESCreate(MEF90DefMechCtx%MEF90Ctx%Comm, snesDisplacement, ierr))
      PetscCall(SNESSetApplicationContext(snesDisplacement, MEF90DefMechCtx, ierr))
      PetscCall(SNESSetDM(snesDisplacement, dm, ierr))
      PetscCall(SNESSetType(snesDisplacement, SNESKSPONLY, ierr))
      PetscCall(SNESSetOptionsPrefix(snesDisplacement, 'Displacement_', ierr))

      PetscCall(SNESSetFunction(snesDisplacement, residual, MEF90DefMechOperatorDisplacement, MEF90DefMechCtx, ierr))
      PetscCall(SNESSetJacobian(snesDisplacement, matDisplacement, matDisplacement, MEF90DefMechBilinearFormDisplacement, MEF90DefMechCtx, ierr))
      atol = 1.0d-7
      rtol = 1.0d-5
      stol = 1.0d-7
      PetscCall(SNESSetTolerances(snesDisplacement, atol, rtol, stol, PETSC_DEFAULT_INTEGER, PETSC_DEFAULT_INTEGER, ierr))
      PetscCall(SNESSetFromOptions(snesDisplacement, ierr))
      !!!
      !!! Set some KSP options
      !!!
      PetscCall(SNESGetKSP(snesDisplacement, kspDisplacement, ierr))
      PetscCall(KSPSetType(kspDisplacement, KSPCG, ierr))
      PetscCall(KSPSetInitialGuessNonzero(kspDisplacement, PETSC_TRUE, ierr))
      atol = 1.0d-8
      rtol = 1.0d-8
      dtol = 1.0d+10
      PetscCall(KSPSetTolerances(kspDisplacement, rtol, atol, dtol, PETSC_DEFAULT_INTEGER, ierr))
      PetscCall(KSPSetFromOptions(kspDisplacement, ierr))
      PetscCall(MatDestroy(matDisplacement, ierr))
   end subroutine MEF90DefMechCreateSNESDisplacement

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCreateSNESDamage"
!!!
!!!
!!!  MEF90DefMechCreateSNESDamage:
!!!
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022 Alexis Marboeuf, marboeua@mcmaster.ca
!!!

   subroutine MEF90DefMechCreateSNESDamage(MEF90DefMechCtx, snesDamage, residual, ierr)
      type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
      type(tSNES), intent(OUT)                            :: snesDamage
      type(tVec), intent(IN)                              :: residual
      PetscErrorCode, intent(INOUT)                       :: ierr

      type(MEF90DefMechGlobalOptions_Type), pointer       :: MEF90DefMechGlobalOptions
      type(tDM)                                          :: dm
      type(tMat)                                         :: matDamage
      type(tKSP)                                         :: kspDamage
      type(tVec)                                         :: UB, LB
      PetscReal                                          :: rtol, dtol, atol, stol

      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag, MEF90DefMechGlobalOptions, ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%damageLocal, dm, ierr))
      PetscCall(DMCreateMatrix(dm, matDamage, ierr))
      PetscCall(MatSetOptionsPrefix(matDamage, "Damage_", ierr))
      !!! The matrix is not symmetric if the advection vector is /= 0
      PetscCall(MatSetOption(matDamage, MAT_SPD, PETSC_TRUE, ierr))
      PetscCall(MatSetOption(matDamage, MAT_SYMMETRY_ETERNAL, PETSC_TRUE, ierr))
      PetscCall(MatSetOption(matDamage, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE, ierr))
      PetscCall(MatSetFromOptions(matDamage, ierr))

      PetscCall(SNESCreate(MEF90DefMechCtx%MEF90Ctx%Comm, snesDamage, ierr))
      PetscCall(SNESSetApplicationContext(snesDamage, MEF90DefMechCtx, ierr))
      PetscCall(SNESSetDM(snesDamage, dm, ierr))
      PetscCall(SNESSetType(snesDamage, SNESVINEWTONRSLS, ierr))
      PetscCall(SNESSetOptionsPrefix(snesDamage, 'Damage_', ierr))

      PetscCall(DMCreateGlobalVector(dm, LB, ierr))
      PetscCall(VecDuplicate(LB, UB, ierr))
      PetscCall(VecSet(LB, 0.0_kr, ierr))
      PetscCall(VecSet(UB, 1.0_kr, ierr))
      PetscCall(SNESVISetVariableBounds(snesDamage, LB, UB, ierr))
      PetscCall(DMRestoreGlobalVector(dm, LB, ierr))
      PetscCall(VecDestroy(UB, ierr))

      PetscCall(SNESSetFunction(snesDamage, residual, MEF90DefMechOperatorDamage, MEF90DefMechCtx, ierr))
      PetscCall(SNESSetJacobian(snesDamage, matDamage, matDamage, MEF90DefMechBilinearFormDamage, MEF90DefMechCtx, ierr))
      atol = 1.0d-7
      rtol = 1.0d-5
      stol = 1.0d-7
      PetscCall(SNESSetTolerances(snesDamage, atol, rtol, stol, PETSC_DEFAULT_INTEGER, PETSC_DEFAULT_INTEGER, ierr))
      PetscCall(SNESSetFromOptions(snesDamage, ierr))
      !!!
      !!! Set some KSP options
      !!!
      PetscCall(SNESGetKSP(snesDamage, kspDamage, ierr))
      PetscCall(KSPSetType(kspDamage, KSPCG, ierr))
      PetscCall(KSPSetInitialGuessNonzero(kspDamage, PETSC_TRUE, ierr))
      atol = 1.0d-8
      rtol = 1.0d-8
      dtol = 1.0d+10
      PetscCall(KSPSetTolerances(kspDamage, rtol, atol, dtol, PETSC_DEFAULT_INTEGER, ierr))
      PetscCall(KSPSetFromOptions(kspDamage, ierr))
      PetscCall(MatDestroy(matDamage, ierr))
   end subroutine MEF90DefMechCreateSNESDamage

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCreateTAODamage"
!!!
!!!
!!!  MEF90DefMechCreateTAODamage:
!!!
!!!  (c) 2022 Alexis Marboeuf, marboeua@mcmaster.ca
!!!

   subroutine MEF90DefMechCreateTAODamage(MEF90DefMechCtx, taoDamage, residual, ierr)
      type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
      type(tTao), intent(OUT)                             :: taoDamage
      type(tVec), intent(IN)                              :: residual
      PetscErrorCode, intent(INOUT)                       :: ierr

      type(MEF90DefMechGlobalOptions_Type), pointer       :: MEF90DefMechGlobalOptions
      type(tDM)                                          :: dm
      type(tMat)                                         :: matDamage
      type(tKSP)                                         :: kspDamage
      type(tVec)                                         :: UB, LB
      PetscReal                                          :: rtol, dtol, atol, stol

      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag, MEF90DefMechGlobalOptions, ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%damageLocal, dm, ierr))
      PetscCall(DMCreateMatrix(dm, matDamage, ierr))
      PetscCall(MatSetOptionsPrefix(matDamage, "Damage_", ierr))
      !!! The matrix is not symmetric if the advection vector is /= 0
      PetscCall(MatSetOption(matDamage, MAT_SPD, PETSC_TRUE, ierr))
      PetscCall(MatSetOption(matDamage, MAT_SYMMETRY_ETERNAL, PETSC_TRUE, ierr))
      PetscCall(MatSetOption(matDamage, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE, ierr))
      PetscCall(MatSetFromOptions(matDamage, ierr))

      PetscCall(TAOCreate(MEF90DefMechCtx%MEF90Ctx%Comm, taoDamage, ierr))
      PetscCall(TAOSetApplicationContext(taoDamage, MEF90DefMechCtx, ierr))
      ! PetscCall(TAOSetDM(taoDamage,dm,ierr))
      PetscCall(TAOSetType(taoDamage, TAOBNTR, ierr))
      PetscCall(TAOSetOptionsPrefix(taoDamage, 'Damage_', ierr))

      PetscCall(DMCreateGlobalVector(dm, LB, ierr))
      PetscCall(VecDuplicate(LB, UB, ierr))
      PetscCall(VecSet(LB, 0.0_kr, ierr))
      PetscCall(VecSet(UB, 1.0_kr, ierr))
      PetscCall(TAOSetVariableBounds(taoDamage, LB, UB, ierr))
      PetscCall(DMRestoreGlobalVector(dm, LB, ierr))
      PetscCall(VecDestroy(UB, ierr))

      PetscCall(TAOSetObjective(taoDamage, MEF90DefMechTAOObjectiveDamage, MEF90DefMechCtx, ierr))
      PetscCall(TAOSetGradient(taoDamage, residual, MEF90DefMechTAOGradientDamage, MEF90DefMechCtx, ierr))
      PetscCall(TAOSetHessian(taoDamage, matDamage, matDamage, MEF90DefMechTAOHessianDamage, MEF90DefMechCtx, ierr))
      atol = 1.0d-5
      rtol = 1.0d-4
      stol = 1.0d-5
      PetscCall(TAOSetTolerances(taoDamage, atol, rtol, stol, ierr))
      PetscCall(TAOSetFromOptions(taoDamage, ierr))
      ! !!!
      ! !!! Set some KSP options
      ! !!!
      PetscCall(TAOGetKSP(taoDamage, kspDamage, ierr))
      PetscCall(KSPSetType(kspDamage, KSPSTCG, ierr))
      PetscCall(KSPSetInitialGuessNonzero(kspDamage, PETSC_TRUE, ierr))
      atol = 1.0d-8
      rtol = 1.0d-8
      dtol = 1.0d+10
      PetscCall(KSPSetTolerances(kspDamage, rtol, atol, dtol, PETSC_DEFAULT_INTEGER, ierr))
      PetscCall(KSPSetFromOptions(kspDamage, ierr))
      PetscCall(MatDestroy(matDamage, ierr))
   end subroutine MEF90DefMechCreateTAODamage

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechUpdateDamageBounds"
!!!
!!!
!!!  MEF90DefMechUpdateDamageBounds:
!!!
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022 Alexis Marboeuf marboeua@mcmaster.ca
!!!

   subroutine MEF90DefMechUpdateDamageBounds(MEF90DefMechCtx, snesDamage, alpha, ierr)
      type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
      type(tSNES), intent(INOUT)                          :: snesDamage
      type(tVec), intent(IN)                              :: alpha
      PetscErrorCode, intent(INOUT)                       :: ierr

      type(tDM)                                          :: dm
      type(tVec)                                         :: LB, UB
      PetscReal, dimension(:), pointer                     :: LBPtr
      PetscInt                                           :: i
      type(MEF90DefMechGlobalOptions_Type), pointer       :: MEF90DefMechGlobalOptions

      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag, MEF90DefMechGlobalOptions, ierr))
      PetscCall(VecGetDM(alpha, dm, ierr))
      PetscCall(DMGetGlobalVector(dm, LB, ierr))
      PetscCall(DMGetGlobalVector(dm, UB, ierr))

      PetscCall(VecSet(UB, 1.0_kr, ierr))
      PetscCall(VecCopy(alpha, LB, ierr))
      if (MEF90DefMechGlobalOptions%irrevthres > 0.0_kr) then
         PetscCall(VecGetArray(LB, LBPtr, ierr))
         do i = 1, size(LBPtr)
            if (LBPtr(i) <= MEF90DefMechGlobalOptions%irrevthres) then
               LBPtr(i) = 0.0_kr
            end if
         end do
         PetscCall(VecRestoreArray(LB, LBPtr, ierr))
      end if
      PetscCall(SNESVISetVariableBounds(snesDamage, LB, UB, ierr))
      PetscCall(DMRestoreGlobalVector(dm, LB, ierr))
      PetscCall(DMRestoreGlobalVector(dm, UB, ierr))
   end subroutine MEF90DefMechUpdateDamageBounds

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechTAOUpdateDamageBounds"
!!!
!!!
!!!  MEF90DefMechTAOUpdateDamageBounds:
!!!
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022 Alexis Marboeuf marboeua@mcmaster.ca
!!!

   subroutine MEF90DefMechTAOUpdateDamageBounds(MEF90DefMechCtx, taoDamage, alpha, ierr)
      type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
      type(tTao), intent(INOUT)                           :: taoDamage
      type(tVec), intent(IN)                              :: alpha
      PetscErrorCode, intent(INOUT)                       :: ierr

      type(tDM)                                          :: dm
      type(tVec)                                         :: LB, UB
      PetscReal, dimension(:), pointer                     :: LBPtr
      PetscInt                                           :: i
      type(MEF90DefMechGlobalOptions_Type), pointer       :: MEF90DefMechGlobalOptions

      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag, MEF90DefMechGlobalOptions, ierr))
      PetscCall(VecGetDM(alpha, dm, ierr))
      PetscCall(DMGetGlobalVector(dm, LB, ierr))
      PetscCall(DMGetGlobalVector(dm, UB, ierr))

      PetscCall(VecSet(UB, 1.0_kr, ierr))
      PetscCall(VecCopy(alpha, LB, ierr))
      if (MEF90DefMechGlobalOptions%irrevthres > 0.0_kr) then
         PetscCall(VecGetArray(LB, LBPtr, ierr))
         do i = 1, size(LBPtr)
            if (LBPtr(i) <= MEF90DefMechGlobalOptions%irrevthres) then
               LBPtr(i) = 0.0_kr
            end if
         end do
         PetscCall(VecRestoreArray(LB, LBPtr, ierr))
      end if
      PetscCall(TAOSetVariableBounds(taoDamage, LB, UB, ierr))
      PetscCall(DMRestoreGlobalVector(dm, LB, ierr))
      PetscCall(DMRestoreGlobalVector(dm, UB, ierr))
   end subroutine MEF90DefMechTAOUpdateDamageBounds

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechPlasticStrainUpdate"
!!!
!!!
!!!  MEF90DefMechPlasticStrainUpdate: wraps calls to MEF90DefMechPlasticStrainUpdate from m_MEF90_DefMechPlasticity
!!!                        since overloading cannot be used here
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   subroutine MEF90DefMechPlasticStrainUpdate(MEF90DefMechCtx, plasticStrain, x, PlasticStrainOld, plasticStrainPrevious, cumulatedDissipatedPlasticEnergyVariation, cumulatedDissipatedPlasticEnergyOld, ierr)
      type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
      type(tVec), intent(INOUT)                           :: plasticStrain
      type(tVec), intent(IN)                              :: x, PlasticStrainOld, plasticStrainPrevious, cumulatedDissipatedPlasticEnergyVariation, cumulatedDissipatedPlasticEnergyOld
      PetscErrorCode, intent(INOUT)                       :: ierr

      PetscInt                                           :: dim
      PetscCall(DMGetDimension(MEF90DefMechCtx%megaDM, dim, ierr))
      if (dim == 2) then
         ! Call MEF90DefMechPlasticStrainUpdate2D(MEF90DefMechCtx,plasticStrain,x,PlasticStrainOld,plasticStrainPrevious,cumulatedDissipatedPlasticEnergyVariation,cumulatedDissipatedPlasticEnergyOld,ierr)
      else if (dim == 3) then
         ! Call MEF90DefMechPlasticStrainUpdate3D(MEF90DefMechCtx,plasticStrain,x,PlasticStrainOld,plasticStrainPrevious,cumulatedDissipatedPlasticEnergyVariation,cumulatedDissipatedPlasticEnergyOld,ierr)
      end if
   end subroutine MEF90DefMechPlasticStrainUpdate

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechProjectCoordinates_Private"

   subroutine MEF90DefMechProjectCoordinates_Private(v, ierr)
      type(tVec), intent(INOUT)           :: v
      PetscErrorCode, intent(INOUT)       :: ierr

      PetscInt                           :: pStart, pEnd, p, numDof, cNumDof, i
      type(tDM)                          :: dm
      type(tPetscSection)                :: coordSection, s
      type(tVec)                         :: coordVec, locV
      PetscScalar, dimension(:), pointer   :: coordArray, vArray
      PetscScalar, dimension(3)           :: xyz
      PetscInt                           :: dim, pOffset

      PetscCall(VecGetDM(v, dm, ierr))
      PetscCall(DMGetLocalSection(dm, s, ierr))
      PetscCall(PetscSectionGetChart(s, pStart, pEnd, ierr))
      PetscCall(DMGetCoordinateSection(dm, coordSection, ierr))
      PetscCall(DMGetCoordinatesLocal(dm, coordVec, ierr))
      PetscCall(DMGetDimension(dm, dim, ierr))

      PetscCall(DMGetLocalVector(dm, locV, ierr))
      PetscCall(VecGetArray(locV, vArray, ierr))

      do p = pStart, pEnd - 1
         PetscCall(PetscSectionGetDof(s, p, numDof, ierr))
         PetscCall(PetscSectionGetConstraintDof(s, p, cNumDof, ierr))
         if ((numDof > 0) .and. (cNumDof == 0)) then
            !!! trick: the coordinate of a point is the average of the coordinates of the points in its closure
            PetscCall(DMPlexVecGetClosure(dm, coordSection, coordVec, p, PETSC_NULL_INTEGER, coordArray, ierr))
            do i = 1, dim
               xyz(i) = sum(coordArray(i:size(coordArray):dim)) * dim / size(coordArray)
            end do
            PetscCall(DMPlexVecRestoreClosure(dm, coordSection, coordVec, p, PETSC_NULL_INTEGER, coordArray, ierr))

            PetscCall(PetscSectionGetOffset(s, p, pOffset, ierr))
            do i = 1, numDof
               vArray(pOffset + i) = xyz(i)
            end do
         end if
      end do
      PetscCall(VecRestoreArray(locV, vArray, ierr))
      PetscCall(DMLocalToGlobalBegin(dm, locV, INSERT_VALUES, v, ierr))
      PetscCall(DMLocalToGlobalEnd(dm, locV, INSERT_VALUES, v, ierr))
      PetscCall(DMRestoreLocalVector(dm, locV, ierr))
   end subroutine MEF90DefMechProjectCoordinates_Private
end module m_MEF90_DefMech
