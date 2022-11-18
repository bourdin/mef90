#include "../MEF90/mef90.inc"
Module m_MEF90_DefMech
#include "petsc/finclude/petsc.h"
   Use m_MEF90
   Use m_MEF90_DefMechCtx
   Use m_MEF90_DefMechAT

   Use m_MEF90_DefMechAssembly2D, &
      MEF90DefMechOperatorDisplacement2D     => MEF90DefMechOperatorDisplacement,      &
      MEF90DefMechBilinearFormDisplacement2D => MEF90DefMechBilinearFormDisplacement,  &     
      MEF90DefMechWork2D                     => MEF90DefMechWork,                      &
      MEF90DefMechCohesiveEnergy2D           => MEF90DefMechCohesiveEnergy,            &
      MEF90DefMechPlasticdissipation2D       => MEF90DefMechPlasticDissipation,        &
      MEF90DefMechElasticEnergy2D            => MEF90DefMechElasticEnergy,             &
      MEF90DefMechOperatorDamage2D           => MEF90DefMechOperatorDamage,            &
      MEF90DefMechBilinearFormDamage2D       => MEF90DefMechBilinearFormDamage,        &
      MEF90DefMechSurfaceEnergy2D            => MEF90DefMechSurfaceEnergy,             &
      MEF90DefMechCrackVolume2D              => MEF90DefMechCrackVolume,               &
      MEF90DefMechStress2D                   => MEF90DefMechStress
   Use m_MEF90_DefMechPlasticity2D, &
      MEF90DefMechPlasticStrainUpdate2D      => MEF90DefMechPlasticStrainUpdate
   Use m_MEF90_DefMechAssembly3D, &
      MEF90DefMechOperatorDisplacement3D     => MEF90DefMechOperatorDisplacement,      &
      MEF90DefMechBilinearFormDisplacement3D => MEF90DefMechBilinearFormDisplacement,  &     
      MEF90DefMechWork3D                     => MEF90DefMechWork,                      &
      MEF90DefMechCohesiveEnergy3D           => MEF90DefMechCohesiveEnergy,            &
      MEF90DefMechPlasticdissipation3D       => MEF90DefMechPlasticDissipation,        &
      MEF90DefMechElasticEnergy3D            => MEF90DefMechElasticEnergy,             &
      MEF90DefMechOperatorDamage3D           => MEF90DefMechOperatorDamage,            &
      MEF90DefMechBilinearFormDamage3D       => MEF90DefMechBilinearFormDamage,        &
      MEF90DefMechSurfaceEnergy3D            => MEF90DefMechSurfaceEnergy,             &
      MEF90DefMechCrackVolume3D              => MEF90DefMechCrackVolume,               &
      MEF90DefMechStress3D                   => MEF90DefMechStress
   Use m_MEF90_DefMechPlasticity3D, &
      MEF90DefMechPlasticStrainUpdate3D      => MEF90DefMechPlasticStrainUpdate

   Implicit none
   !Private
   Public :: MEF90DefMechSetTransients
   Public :: MEF90DefMechOperatorDisplacement
   Public :: MEF90DefMechBilinearFormDisplacement
   Public :: MEF90DefMechCreateSNESDisplacement

   Public :: MEF90DefMechOperatorDamage
   Public :: MEF90DefMechBilinearFormDamage
   Public :: MEF90DefMechCreateSNESDamage
   Public :: MEF90DefMechUpdateDamageBounds
   
   Public :: MEF90DefMechFormatEXO
   Public :: MEF90DefMechViewEXO
   Public :: MEF90DefMechSurfaceEnergy
   Public :: MEF90DefMechElasticEnergy
   Public :: MEF90DefMechWork
   Public :: MEF90DefMechCohesiveEnergy
   Public :: MEF90DefMechPlasticDissipation
   Public :: MEF90DefMechCrackVolume
   Public :: MEF90DefMechStress
   Public :: MEF90DefMechPlasticStrainUpdate

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSetTransients"
!!!
!!!  
!!!  MEF90DefMechSetTransients:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!  (c)    2022 Alexis Marboeuf marboeua@mcmaster.ca
!!!

   Subroutine MEF90DefMechSetTransients(MEF90DefMechCtx,step,time,ierr)
      Type(MEF90DefMechCtx_Type),Intent(INOUT)        :: MEF90DefMechCtx
      PetscInt,Intent(IN)                             :: step
      PetscReal,Intent(IN)                            :: time
      PetscErrorCode,Intent(INOUT)                    :: ierr
   
      Type(MEF90DefMechGlobalOptions_Type),pointer    :: MEF90DefMechGlobalOptions
      Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90GlobalOptions
      Type(tDM)                                       :: dmDisplacement,dmDamage,dmCohesiveDisplacement
      Type(tVec)                                      :: tmpVec

      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr))
      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr))

      PetscCall(VecGetDM(MEF90DefMechCtx%damageLocal,dmDamage,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%displacementLocal,dmDisplacement,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%cohesiveDisplacement,dmCohesiveDisplacement,ierr))

      Select case (MEF90DefMechGlobalOptions%boundaryDisplacementScaling)
      Case (MEF90Scaling_File)
         PetscCall(DMGetLocalVector(dmDisplacement,tmpVec,ierr))
         PetscCall(PetscObjectSetName(tmpVec,"Displacement",ierr))
         PetscCall(MEF90EXOVecLoad(tmpVec,MEF90DefMechCtx%displacementToIOSF,MEF90DefMechCtx%IOToDisplacementSF,MEF90DefMechCtx%MEF90Ctx%resultViewer,step,ierr))
         PetscCall(MEF90VecCopySF(tmpVec,MEF90DefMechCtx%displacementLocal,MEF90DefMechCtx%boundaryToDisplacementSF,ierr))
         PetscCall(DMRestoreLocalVector(dmDisplacement,tmpVec,ierr))
      Case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%displacementLocal,time,ierr))
      Case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%displacementLocal,1.0_Kr,ierr))
      End Select

      Select case (MEF90DefMechGlobalOptions%boundaryDamageScaling)
      Case (MEF90Scaling_File)
         PetscCall(DMGetLocalVector(dmDamage,tmpVec,ierr))
         PetscCall(PetscObjectSetName(tmpVec,"Damage",ierr))
         PetscCall(MEF90EXOVecLoad(tmpVec,MEF90DefMechCtx%damageToIOSF,MEF90DefMechCtx%IOToDamageSF,MEF90DefMechCtx%MEF90Ctx%resultViewer,step,ierr))
         PetscCall(MEF90VecCopySF(tmpVec,MEF90DefMechCtx%damageLocal,MEF90DefMechCtx%boundaryToDamageSF,ierr))
         PetscCall(DMRestoreLocalVector(dmDamage,tmpVec,ierr))
      Case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%damageLocal,time,ierr))
      Case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%damageLocal,1.0_Kr,ierr))
      End Select

      Select case (MEF90DefMechGlobalOptions%cohesiveDisplacementScaling)
      Case (MEF90Scaling_File)
         PetscCall(MEF90EXOVecLoad(MEF90DefMechCtx%cohesiveDisplacement,MEF90DefMechCtx%cohesiveDisplacementToIOSF,MEF90DefMechCtx%IOToCohesiveDisplacementSF,MEF90DefMechCtx%MEF90Ctx%resultViewer,step,ierr))
      Case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%cohesiveDisplacement,time,ierr))
      Case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%cohesiveDisplacement,1.0_Kr,ierr))
      End Select

      Select case (MEF90DefMechGlobalOptions%displacementLowerBoundScaling)
      Case (MEF90Scaling_File)
         Write(*,*) __FUNCT__,": file scaling of displacement lower bound does not make any sense."
         STOP
      Case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%displacementLowerBoundLocal,time,ierr))
      Case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%displacementLowerBoundLocal,1.0_Kr,ierr))
      End Select

      Select case (MEF90DefMechGlobalOptions%displacementUpperBoundScaling)
      Case (MEF90Scaling_File)
         Write(*,*) __FUNCT__,": file scaling of displacement upper bound does not make any sense."
         STOP
      Case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%displacementUpperBoundLocal,time,ierr))
      Case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%displacementUpperBoundLocal,1.0_Kr,ierr))
      End Select

      Select case (MEF90DefMechGlobalOptions%bodyForceScaling)
      Case (MEF90Scaling_File)
         PetscCall(MEF90EXOVecLoad(MEF90DefMechCtx%bodyForce,MEF90DefMechCtx%bodyForceToIOSF,MEF90DefMechCtx%IOToBodyForceSF,MEF90DefMechCtx%MEF90Ctx%resultViewer,step,ierr))
      Case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%bodyForce,time,ierr))
      Case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%bodyForce,1.0_Kr,ierr))
      End Select

      Select case (MEF90DefMechGlobalOptions%boundaryForceScaling)
      Case (MEF90Scaling_File)
         PetscCall(MEF90EXOVecLoad(MEF90DefMechCtx%boundaryForce,MEF90DefMechCtx%boundaryForceToIOSF,MEF90DefMechCtx%IOToBoundaryForceSF,MEF90DefMechCtx%MEF90Ctx%resultViewer,step,ierr))
      Case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%boundaryForce,time,ierr))
      Case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%boundaryForce,1.0_Kr,ierr))
      End Select

      Select case (MEF90DefMechGlobalOptions%pressureForceScaling)
      Case (MEF90Scaling_File)
         PetscCall(MEF90EXOVecLoad(MEF90DefMechCtx%pressureForce,MEF90DefMechCtx%pressureForceToIOSF,MEF90DefMechCtx%IOToPressureForceSF,MEF90DefMechCtx%MEF90Ctx%resultViewer,step,ierr))
      Case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%pressureForce,time,ierr))
      Case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%pressureForce,1.0_Kr,ierr))
      End Select

   End Subroutine MEF90DefMechSetTransients

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDisplacement"
!!!
!!!  
!!!  MEF90DefMechOperatorDisplacement: wraps calls to MEF90DefMechOperatorDisplacement from m_MEF90_DefMechAssembly
!!!                        since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDisplacement(snesTemp,x,residual,MEF90DefMechCtx,ierr)
      Type(tSNES),Intent(IN)                             :: snesTemp
      Type(tVec),Intent(IN)                              :: x
      Type(tVec),Intent(INOUT)                           :: residual
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(INOUT)                       :: ierr
      
      PetscInt                                           :: dim      
      PetscCall(DMGetDimension(MEF90DefMechCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         PetscCall(MEF90DefMechOperatorDisplacement2D(snesTemp,x,residual,MEF90DefMechCtx,ierr))
      Else If (dim == 3) Then
         PetscCall(MEF90DefMechOperatorDisplacement3D(snesTemp,x,residual,MEF90DefMechCtx,ierr))
      End If      
   End Subroutine MEF90DefMechOperatorDisplacement
   
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDisplacement"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDisplacement: wraps calls to MEF90DefMechBilinearFormDisplacement from m_MEF90_DefMechAssembly
!!!                            since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDisplacement(snesDispl,x,A,M,flg,MEF90DefMechCtx,ierr)
      Type(tSNES),Intent(IN)                             :: snesDispl
      Type(tVec),Intent(IN)                              :: x
      Type(tMat),Intent(INOUT)                           :: A,M
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(INOUT)                       :: ierr

      PetscInt                                           :: dim      
      PetscCall(DMGetDimension(MEF90DefMechCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         PetscCall(MEF90DefMechBilinearFormDisplacement2D(snesDispl,x,A,M,flg,MEF90DefMechCtx,ierr))
      Else If (dim == 3) Then
         PetscCall(MEF90DefMechBilinearFormDisplacement3D(snesDispl,x,A,M,flg,MEF90DefMechCtx,ierr))
      End If      
   End Subroutine MEF90DefMechBilinearFormDisplacement

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechWork"
!!!
!!!  
!!!  MEF90DefMechWork: wraps calls to MEF90DefMechWork from m_MEF90_DefMechAssembly
!!                       since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechWork(MEF90DefMechCtx,work,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)           :: MEF90DefMechCtx
      PetscReal,Dimension(:),Pointer                  :: work
      PetscErrorCode,Intent(INOUT)                    :: ierr

      PetscInt                                        :: dim      
      PetscCall(DMGetDimension(MEF90DefMechCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         PetscCall(MEF90DefMechWork2D(MEF90DefMechCtx,work,ierr))
      Else If (dim == 3) Then
         PetscCall(MEF90DefMechWork3D(MEF90DefMechCtx,work,ierr))
      End If      
   End Subroutine MEF90DefMechWork

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCohesiveEnergy"
!!!
!!!  
!!!  MEF90DefMechCohesiveEnergy: wraps calls to MEF90DefMechCohesiveEnergy from m_MEF90_DefMechAssembly
!!                       since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechCohesiveEnergy(MEF90DefMechCtx,cohesiveEnergy,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)           :: MEF90DefMechCtx
      PetscReal,Dimension(:),Pointer                  :: cohesiveEnergy
      PetscErrorCode,Intent(INOUT)                    :: ierr

      PetscInt                                        :: dim      
      PetscCall(DMGetDimension(MEF90DefMechCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         PetscCall(MEF90DefMechCohesiveEnergy2D(MEF90DefMechCtx,cohesiveEnergy,ierr))
      Else If (dim == 3) Then
         PetscCall(MEF90DefMechCohesiveEnergy3D(MEF90DefMechCtx,cohesiveEnergy,ierr))
      End If      
   End Subroutine MEF90DefMechCohesiveEnergy

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechElasticEnergy"
!!!
!!!  
!!!  MEF90DefMechElasticEnergy: wraps calls to MEF90DefMechElasticEnergy from m_MEF90_DefMechAssembly
!!                       since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechElasticEnergy(MEF90DefMechCtx,energy,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,dimension(:),Pointer                     :: energy
      PetscErrorCode,Intent(INOUT)                       :: ierr

      PetscInt                                           :: dim      
      PetscCall(DMGetDimension(MEF90DefMechCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         PetscCall(MEF90DefMechElasticEnergy2D(MEF90DefMechCtx,energy,ierr))
      Else If (dim == 3) Then
         PetscCall(MEF90DefMechElasticEnergy3D(MEF90DefMechCtx,energy,ierr))
      End If      
   End Subroutine MEF90DefMechElasticEnergy


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechPlasticDissipation"
!!!
!!!  
!!!  MEF90DefMechPlasticDissipation: wraps calls to MEF90DefMechPlasticDissipation from m_MEF90_DefMechAssembly
!!                       since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Erwan TANNE erwan.tanne@gmail.com
!!!

   Subroutine MEF90DefMechPlasticDissipation(x,MEF90DefMechCtx,plasticStrainOld,energy,ierr)
      Type(tVec),Intent(IN)                              :: x
      Type(tVec),Intent(IN)                              :: plasticStrainOld
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,dimension(:),Pointer                     :: energy
      PetscErrorCode,Intent(INOUT)                       :: ierr

      PetscInt                                           :: dim      
      PetscCall(DMGetDimension(MEF90DefMechCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         ! Call MEF90DefMechPlasticDissipation2D(x,MEF90DefMechCtx,plasticStrainOld,energy,ierr)
      Else If (dim == 3) Then
         ! Call MEF90DefMechPlasticDissipation3D(x,MEF90DefMechCtx,plasticStrainOld,energy,ierr)
      End If      
   End Subroutine MEF90DefMechPlasticDissipation

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechStress"
!!!
!!!  
!!!  MEF90DefMechStress: wraps calls to MEF90DefMechElasticEnergy from m_MEF90_DefMechAssembly
!!                       since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechStress(MEF90DefMechCtx,stress,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      Type(tVec),Intent(IN)                              :: stress
      PetscErrorCode,Intent(INOUT)                       :: ierr

      PetscInt                                           :: dim      
      PetscCall(DMGetDimension(MEF90DefMechCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         PetscCall(MEF90DefMechStress2D(MEF90DefMechCtx,stress,ierr))
      Else If (dim == 3) Then
         PetscCall(MEF90DefMechStress3D(MEF90DefMechCtx,stress,ierr))
      End If      
   End Subroutine MEF90DefMechStress


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCrackVolume"
!!!
!!!  
!!!  MEF90DefMechCrackVolume: 
!!!  
!!!  (c) 2016 erwan
!!!

   Subroutine MEF90DefMechCrackVolume(MEF90DefMechCtx,CrackVolume,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,dimension(:),Pointer                     :: CrackVolume
      PetscErrorCode,Intent(INOUT)                       :: ierr

      PetscInt                                           :: dim      
      PetscCall(DMGetDimension(MEF90DefMechCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         PetscCall(MEF90DefMechCrackVolume2D(MEF90DefMechCtx,CrackVolume,ierr))
      Else If (dim == 3) Then
         PetscCall(MEF90DefMechCrackVolume3D(MEF90DefMechCtx,CrackVolume,ierr))
      End If      
   End Subroutine MEF90DefMechCrackVolume

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamage"
!!!
!!!  
!!!  MEF90DefMechOperatorDamage: wraps calls to MEF90DefMechOperatorDamage from m_MEF90_DefMechAssembly
!!!                        since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDamage(snesTemp,x,residual,MEF90DefMechCtx,ierr)
      Type(tSNES),Intent(IN)                             :: snesTemp
      Type(tVec),Intent(IN)                              :: x
      Type(tVec),Intent(INOUT)                           :: residual
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(INOUT)                       :: ierr
      
      PetscInt                                           :: dim      
      PetscCall(DMGetDimension(MEF90DefMechCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         PetscCall(MEF90DefMechOperatorDamage2D(snesTemp,x,residual,MEF90DefMechCtx,ierr))
      Else If (dim == 3) Then
         PetscCall(MEF90DefMechOperatorDamage3D(snesTemp,x,residual,MEF90DefMechCtx,ierr))
      End If      
   End Subroutine MEF90DefMechOperatorDamage
   
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamage"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamage: wraps calls to MEF90DefMechBilinearFormDamage from m_MEF90_DefMechAssembly
!!!                            since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDamage(snesDispl,x,A,M,flg,MEF90DefMechCtx,ierr)
      Type(tSNES),Intent(IN)                             :: snesDispl
      Type(tVec),Intent(IN)                              :: x
      Type(tMat),Intent(INOUT)                           :: A,M
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(INOUT)                       :: ierr

      PetscInt                                           :: dim      
      PetscCall(DMGetDimension(MEF90DefMechCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         PetscCall(MEF90DefMechBilinearFormDamage2D(snesDispl,x,A,M,flg,MEF90DefMechCtx,ierr))
      Else If (dim == 3) Then
         PetscCall(MEF90DefMechBilinearFormDamage3D(snesDispl,x,A,M,flg,MEF90DefMechCtx,ierr))
      End If      
   End Subroutine MEF90DefMechBilinearFormDamage

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSurfaceEnergy"
!!!
!!!  
!!!  MEF90DefMechSurfaceEnergy: wraps calls to MEF90DefMechSurfaceEnergy from m_MEF90_DefMechAssembly
!!                       since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechSurfaceEnergy(MEF90DefMechCtx,energy,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,dimension(:),Pointer                     :: energy
      PetscErrorCode,Intent(INOUT)                       :: ierr

      PetscInt                                           :: dim      
      PetscCall(DMGetDimension(MEF90DefMechCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         PetscCall(MEF90DefMechSurfaceEnergy2D(MEF90DefMechCtx,energy,ierr))
      Else If (dim == 3) Then
         PetscCall(MEF90DefMechSurfaceEnergy3D(MEF90DefMechCtx,energy,ierr))
      End If      
   End Subroutine MEF90DefMechSurfaceEnergy

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechViewEXO"
!!!
!!!  
!!!  MEF90DefMechViewEXO: Save all fields in a MEF90DefMechCtx_Type in an exodus file
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022 Alexis Marboeuf, marboeua@mcmaster.ca
!!!

   Subroutine MEF90DefMechViewEXO(MEF90DefMechCtx,step,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscInt,Intent(IN)                                :: step
      PetscErrorCode,Intent(INOUT)                       :: ierr

      PetscCall(MEF90EXOVecView(MEF90DefMechCtx%displacementLocal,MEF90DefMechCtx%displacementToIOSF,MEF90DefMechCtx%IOToDisplacementSF,MEF90DefMechCtx%MEF90Ctx%resultViewer,step,ierr))
      PetscCall(MEF90EXOVecView(MEF90DefMechCtx%damageLocal,MEF90DefMechCtx%damageToIOSF,MEF90DefMechCtx%IOToDamageSF,MEF90DefMechCtx%MEF90Ctx%resultViewer,step,ierr))
   End Subroutine MEF90DefMechViewEXO
   
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechFormatEXO"
!!!
!!!  
!!!  MEF90DefMechFormatEXO:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechFormatEXO(MEF90DefMechCtx,time,ierr)
      Type(MEF90DefMechCtx_Type),Intent(INOUT)           :: MEF90DefMechCtx
      PetscReal,Dimension(:),Pointer                     :: time
      PetscErrorCode,Intent(INOUT)                       :: ierr

      ! Character(len=MXSTLN),Dimension(:),Pointer         :: nameG,nameV,nameC
      ! Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
      ! Integer                                            :: dim,numfield,step

      ! Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)
      ! PetscCall(DMGetDimension(MEF90DefMechCtx%megaDM,dim,ierr))
      ! Allocate(nameG(0))
      ! !nameG(1) = "Elastic Energy"
      ! !nameG(2) = "Work"
      ! !nameG(3) = "Surface Energy"
      ! !nameG(4) = "Total Energy"
   
      ! numfield = max(MEF90DefMechGlobalOptions%displacementOffset+dim-1, &
      !                MEF90DefMechGlobalOptions%damageOffset,&
      !                MEF90DefMechGlobalOptions%boundaryDisplacementOffset+dim-1,&
      !                MEF90DefMechGlobalOptions%boundaryDamageOffset,&
      !                MEF90DefMechGlobalOptions%temperatureOffset)
      ! Allocate(nameV(numfield))
      ! nameV = "empty"
      ! If ((MEF90DefMechGlobalOptions%boundaryDisplacementOffset > 0) .AND. &
      !     (MEF90DefMechGlobalOptions%boundaryDisplacementOffset /= MEF90DefMechGlobalOptions%displacementOffset)) Then
      !    nameV(MEF90DefMechGlobalOptions%boundaryDisplacementOffset+0)    = "Boundary_Displacement_X"
      !    nameV(MEF90DefMechGlobalOptions%boundaryDisplacementOffset+1)    = "Boundary_Displacement_Y"
      !    If (dim == 3) Then
      !       nameV(MEF90DefMechGlobalOptions%boundaryDisplacementOffset+2) = "Boundary_Displacement_Z"
      !    End If
      ! End If
      ! If (MEF90DefMechGlobalOptions%displacementOffset > 0) Then
      !    nameV(MEF90DefMechGlobalOptions%displacementOffset+0)            = "Displacement_X"
      !    nameV(MEF90DefMechGlobalOptions%displacementOffset+1)            = "Displacement_Y"
      !    If (dim == 3) Then
      !       nameV(MEF90DefMechGlobalOptions%displacementOffset+2)         = "Displacement_Z"
      !    End If
      ! End If
      ! If (MEF90DefMechGlobalOptions%damageOffset > 0) Then
      !    nameV(MEF90DefMechGlobalOptions%damageOffset)                    = "Damage"
      ! End If
      ! If ((MEF90DefMechGlobalOptions%boundaryDamageOffset > 0) .AND. &
      !     (MEF90DefMechGlobalOptions%boundaryDamageOffset /= MEF90DefMechGlobalOptions%damageOffset))Then
      !    nameV(MEF90DefMechGlobalOptions%boundaryDamageOffset)            = "Boundary_Damage"
      ! End If
      ! If (MEF90DefMechGlobalOptions%temperatureOffset > 0) Then
      !    nameV(MEF90DefMechGlobalOptions%temperatureOffset)               = "Temperature"
      ! End If
      
      ! numfield = 0
      ! If (MEF90DefMechGlobalOptions%forceOffset > 0) Then
      !    numfield = max(numfield,MEF90DefMechGlobalOptions%forceOffset+dim-1)
      ! End If
      
      ! If (MEF90DefMechGlobalOptions%pressureForceOffset > 0) Then
      !    numfield = max(numfield,MEF90DefMechGlobalOptions%pressureForceOffset)
      ! End If
      
      ! If (MEF90DefMechGlobalOptions%CrackPressureOffset > 0) Then
      !    numfield = max(numfield,MEF90DefMechGlobalOptions%CrackPressureOffset)
      ! End If

      ! If (MEF90DefMechGlobalOptions%stressOffset > 0) Then
      !    numfield = max(numfield,MEF90DefMechGlobalOptions%stressOffset+(dim*(dim+1))/2-1)
      ! End If

      ! If (MEF90DefMechGlobalOptions%plasticStrainOffset > 0) Then
      !    numfield = max(numfield,MEF90DefMechGlobalOptions%plasticStrainOffset+(dim*(dim+1))/2-1)
      ! End If
      ! If (MEF90DefMechGlobalOptions%cumulatedPlasticDissipationOffset > 0) Then
      !    numfield = max(numfield,MEF90DefMechGlobalOptions%cumulatedPlasticDissipationOffset)
      ! End If
      ! Allocate(nameC(numfield))

      ! If (numfield > 0) Then
      !    nameC = "empty"
      ! End If
      ! If (MEF90DefMechGlobalOptions%forceOffset > 0) Then
      !    nameC(MEF90DefMechGlobalOptions%forceOffset+0)                                = "Force_X"
      !    nameC(MEF90DefMechGlobalOptions%forceOffset+1)                                = "Force_Y"
      !    If (dim == 3) Then
      !       nameC(MEF90DefMechGlobalOptions%forceOffset+2)                             = "Force_Z"
      !    End If
      ! End If
      
      ! If (MEF90DefMechGlobalOptions%pressureForceOffset > 0) Then
      !    nameC(MEF90DefMechGlobalOptions%pressureForceOffset)                          = "Pressure_Force"
      ! End If
      
      ! If (MEF90DefMechGlobalOptions%CrackPressureOffset > 0) Then
      !    nameC(MEF90DefMechGlobalOptions%CrackPressureOffset)                          = "Crack_Pressure"
      ! End If

      ! If (MEF90DefMechGlobalOptions%stressOffset > 0) Then
      !    If (dim == 2) Then
      !       nameC(MEF90DefMechGlobalOptions%stressOffset+0)                            = "Stress_XX"
      !       nameC(MEF90DefMechGlobalOptions%stressOffset+1)                            = "Stress_YY"
      !       nameC(MEF90DefMechGlobalOptions%stressOffset+2)                            = "Stress_XY"
      !    Else
      !       nameC(MEF90DefMechGlobalOptions%stressOffset+0)                            = "Stress_XX"
      !       nameC(MEF90DefMechGlobalOptions%stressOffset+1)                            = "Stress_YY"
      !       nameC(MEF90DefMechGlobalOptions%stressOffset+2)                            = "Stress_ZZ"
      !       nameC(MEF90DefMechGlobalOptions%stressOffset+3)                            = "Stress_YZ"
      !       nameC(MEF90DefMechGlobalOptions%stressOffset+4)                            = "Stress_XZ"
      !       nameC(MEF90DefMechGlobalOptions%stressOffset+5)                            = "Stress_XY"
      !    End If
      ! End If
      ! If (MEF90DefMechGlobalOptions%plasticStrainOffset > 0) Then
      !    If (dim == 2) Then
      !       nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+0)                     = "plasticStrain_XX"
      !       nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+1)                     = "plasticStrain_YY"
      !       nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+2)                     = "plasticStrain_XY"
      !    Else
      !       nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+0)                     = "plasticStrain_XX"
      !       nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+1)                     = "plasticStrain_YY"
      !       nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+2)                     = "plasticStrain_ZZ"
      !       nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+3)                     = "plasticStrain_YZ"
      !       nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+4)                     = "plasticStrain_XZ"
      !       nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+5)                     = "plasticStrain_XY"
      !    End If
      ! End If

      ! If (MEF90DefMechGlobalOptions%cumulatedPlasticDissipationOffset > 0) Then
      !    nameC(MEF90DefMechGlobalOptions%cumulatedPlasticDissipationOffset)       = "Cumulated_Plastic_Dissipation"
      ! End If
      
      ! Call MEF90EXOFormat(MEF90DefMechCtx%MEF90Ctx%fileEXOUNIT,nameG,nameC,nameV,ierr)
      ! !!! This makes no sense, but there seems to be a bug in exodus / OSX where
      ! !!! formatting is not flushed to the drive
      ! Call MEF90CtxCloseEXO(MEF90DefMechCtx%MEF90Ctx,ierr)
      ! Call MEF90CtxOpenEXO(MEF90DefMechCtx%MEF90Ctx,MEF90DefMechCtx%DM,ierr)
      ! If (MEF90DefMechCtx%MEF90Ctx%rank == 0) Then
      !    !Call EXUPDA(MEF90DefMechCtx%MEF90Ctx%fileExoUnit,ierr)
      !    If (associated(time)) then
      !       Do step = 1, size(time)
      !          Call EXPTIM(MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,time(step),ierr)
      !       End Do
      !    End If
      ! End If
      ! DeAllocate(nameG)
      ! DeAllocate(nameV)
      ! DeAllocate(nameC)
   End Subroutine MEF90DefMechFormatEXO
   
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCreateSNESDisplacement"
!!!
!!!  
!!!  MEF90DefMechCreateSNESDisplacement:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022 Alexis Marboeuf, marboeua@mcmaster.ca
!!!

   Subroutine MEF90DefMechCreateSNESDisplacement(MEF90DefMechCtx,snesDisp,residual,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      Type(tSNES),Intent(OUT)                            :: snesDisp
      Type(tVec),Intent(IN)                              :: residual
      PetscErrorCode,Intent(INOUT)                       :: ierr
      
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
      Type(tDM)                                          :: dm
      Type(tMat)                                         :: matTemp
      Type(tMatNullSpace)                                :: nspTemp
      Type(tKSP)                                         :: kspTemp
      PetscReal                                          :: rtol,dtol
      
      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%displacementLocal,dm,ierr))
      PetscCall(DMCreateMatrix(dm,matTemp,iErr))
      PetscCall(MatSetOptionsPrefix(matTemp,"Displacement_",ierr))
      !!! The matrix is not symmetric if the advection vector is /= 0
      PetscCall(MatSetOption(matTemp,MAT_SPD,PETSC_TRUE,ierr))
      PetscCall(MatSetOption(matTemp,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr))
      PetscCall(MatSetOption(matTemp,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr))
      If (MEF90DefMechGlobalOptions%addDisplacementNullSpace) Then
         PetscCall(MatNullSpaceCreate(MEF90DefMechCtx%MEF90Ctx%Comm,PETSC_TRUE,0_Ki,PETSC_NULL_VEC,nspTemp,ierr))
         PetscCall(MatSetNullSpace(matTemp,nspTemp,ierr))
      End If
      PetscCall(MatSetFromOptions(matTemp,ierr))

      PetscCall(SNESCreate(MEF90DefMechCtx%MEF90Ctx%Comm,snesDisp,ierr))
      PetscCall(SNESSetApplicationContext(snesDisp,MEF90DefMechCtx,ierr))
      PetscCall(SNESSetDM(snesDisp,dm,ierr))
      PetscCall(SNESSetType(snesDisp,SNESKSPONLY,ierr))
      PetscCall(SNESSetOptionsPrefix(snesDisp,'Displacement_',ierr))

      PetscCall(SNESSetFunction(snesDisp,residual,MEF90DefMechOperatorDisplacement,MEF90DefMechCtx,ierr))
      PetscCall(SNESSetJacobian(snesDisp,matTemp,matTemp,MEF90DefMechBilinearFormDisplacement,MEF90DefMechCtx,ierr))
      PetscCall(SNESSetFromOptions(snesDisp,ierr))
      !!! 
      !!! Set some KSP options
      !!!
      PetscCall(SNESGetKSP(snesDisp,kspTemp,ierr))
      PetscCall(KSPSetType(kspTemp,KSPCG,ierr))
      PetscCall(KSPSetInitialGuessNonzero(kspTemp,PETSC_TRUE,ierr))
      rtol = 1.0D-8
      dtol = 1.0D+10
      PetscCall(KSPSetTolerances(kspTemp,rtol,PETSC_DEFAULT_REAL,dtol,PETSC_DEFAULT_INTEGER,ierr))
      PetscCall(KSPSetFromOptions(kspTemp,ierr))
      PetscCall(MatDestroy(matTemp,ierr))
   End Subroutine MEF90DefMechCreateSNESDisplacement

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCreateSNESDamage"
!!!
!!!  
!!!  MEF90DefMechCreateSNESDamage:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022 Alexis Marboeuf, marboeua@mcmaster.ca
!!!

   Subroutine MEF90DefMechCreateSNESDamage(MEF90DefMechCtx,snesDamage,residual,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      Type(tSNES),Intent(OUT)                            :: snesDamage
      Type(tVec),Intent(IN)                              :: residual
      PetscErrorCode,Intent(INOUT)                       :: ierr
      
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
      Type(tDM)                                          :: dm
      Type(tMat)                                         :: matTemp
      Type(tKSP)                                         :: kspTemp
      PetscReal                                          :: rtol,dtol
      
      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%damageLocal,dm,ierr))
      PetscCall(DMCreateMatrix(dm,matTemp,iErr))
      PetscCall(MatSetOptionsPrefix(matTemp,"Damage_",ierr))
      !!! The matrix is not symmetric if the advection vector is /= 0
      PetscCall(MatSetOption(matTemp,MAT_SPD,PETSC_TRUE,ierr))
      PetscCall(MatSetOption(matTemp,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr))
      PetscCall(MatSetOption(matTemp,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr))
      ! If (MEF90DefMechGlobalOptions%addDamageNullSpace) Then
      !    PetscCall(MatNullSpaceCreate(MEF90DefMechCtx%MEF90Ctx%Comm,PETSC_TRUE,0_Ki,PETSC_NULL_VEC,nspTemp,ierr))
      !    PetscCall(MatSetNullSpace(matTemp,nspTemp,ierr))
      ! End If
      PetscCall(MatSetFromOptions(matTemp,ierr))

      PetscCall(SNESCreate(MEF90DefMechCtx%MEF90Ctx%Comm,snesDamage,ierr))
      PetscCall(SNESSetApplicationContext(snesDamage,MEF90DefMechCtx,ierr))
      PetscCall(SNESSetDM(snesDamage,dm,ierr))
      PetscCall(SNESSetType(snesDamage,SNESKSPONLY,ierr))
      PetscCall(SNESSetOptionsPrefix(snesDamage,'Damage_',ierr))

      PetscCall(SNESSetFunction(snesDamage,residual,MEF90DefMechOperatorDamage,MEF90DefMechCtx,ierr))
      PetscCall(SNESSetJacobian(snesDamage,matTemp,matTemp,MEF90DefMechBilinearFormDamage,MEF90DefMechCtx,ierr))
      PetscCall(SNESSetFromOptions(snesDamage,ierr))
      !!! 
      !!! Set some KSP options
      !!!
      PetscCall(SNESGetKSP(snesDamage,kspTemp,ierr))
      PetscCall(KSPSetType(kspTemp,KSPCG,ierr))
      PetscCall(KSPSetInitialGuessNonzero(kspTemp,PETSC_TRUE,ierr))
      rtol = 1.0D-8
      dtol = 1.0D+10
      PetscCall(KSPSetTolerances(kspTemp,rtol,PETSC_DEFAULT_REAL,dtol,PETSC_DEFAULT_INTEGER,ierr))
      PetscCall(KSPSetFromOptions(kspTemp,ierr))
      PetscCall(MatDestroy(matTemp,ierr))
   End Subroutine MEF90DefMechCreateSNESDamage

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechUpdateDamageBounds"
!!!
!!!  
!!!  MEF90DefMechUpdateDamageBounds:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechUpdateDamageBounds(MEF90DefMechCtx,snesDamage,alpha,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      Type(tSNES),Intent(OUT)                            :: snesDamage
      Type(tVec),Intent(IN)                              :: alpha
      PetscErrorCode,Intent(INOUT)                       :: ierr
      
      ! Type(Vec)                                          :: LB,UB
      ! PetscReal,Dimension(:),Pointer                     :: LBPtr
      ! PetscInt                                           :: i
      ! Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions

      ! Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)      
      ! Call DMGetGlobalVector(MEF90DefMechCtx%DMScal,LB,ierr);CHKERRQ(ierr)
      ! Call DMGetGlobalVector(MEF90DefMechCtx%DMScal,UB,ierr);CHKERRQ(ierr)

      ! Call VecSet(UB,1.0_Kr,ierr);CHKERRQ(ierr)
      ! Call VecCopy(alpha,LB,ierr);CHKERRQ(ierr)
      ! If (MEF90DefMechGlobalOptions%irrevthres > 0.0_Kr) Then
      !    Call VecGetArrayF90(LB,LBPtr,ierr);CHKERRQ(ierr)
      !    Do i = 1, size(LBPtr)
      !       If (LBPtr(i) <= MEF90DefMechGlobalOptions%irrevthres) Then
      !          LBPtr(i) = 0.0_Kr
      !       End If
      !    End Do
      !    Call VecRestoreArrayF90(LB,LBPtr,ierr);CHKERRQ(ierr)
      ! End If
      ! Call SNESVISetVariableBounds(snesDamage,LB,UB,ierr);CHKERRQ(ierr)
      ! Call DMRestoreGlobalVector(MEF90DefMechCtx%DMScal,LB,ierr);CHKERRQ(ierr)
      ! Call DMRestoreGlobalVector(MEF90DefMechCtx%DMScal,UB,ierr);CHKERRQ(ierr)      
   End Subroutine MEF90DefMechUpdateDamageBounds

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechPlasticStrainUpdate"
!!!
!!!  
!!!  MEF90DefMechPlasticStrainUpdate: wraps calls to MEF90DefMechPlasticStrainUpdate from m_MEF90_DefMechPlasticity
!!!                        since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechPlasticStrainUpdate(MEF90DefMechCtx,plasticStrain,x,PlasticStrainOld,plasticStrainPrevious,cumulatedDissipatedPlasticEnergyVariation,cumulatedDissipatedPlasticEnergyOld,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      Type(tVec),Intent(INOUT)                           :: plasticStrain
      Type(tVec),Intent(IN)                              :: x,PlasticStrainOld,plasticStrainPrevious,cumulatedDissipatedPlasticEnergyVariation,cumulatedDissipatedPlasticEnergyOld
      PetscErrorCode,Intent(INOUT)                       :: ierr
      
      PetscInt                                           :: dim      
      PetscCall(DMGetDimension(MEF90DefMechCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         ! Call MEF90DefMechPlasticStrainUpdate2D(MEF90DefMechCtx,plasticStrain,x,PlasticStrainOld,plasticStrainPrevious,cumulatedDissipatedPlasticEnergyVariation,cumulatedDissipatedPlasticEnergyOld,ierr)
      Else If (dim == 3) Then
         ! Call MEF90DefMechPlasticStrainUpdate3D(MEF90DefMechCtx,plasticStrain,x,PlasticStrainOld,plasticStrainPrevious,cumulatedDissipatedPlasticEnergyVariation,cumulatedDissipatedPlasticEnergyOld,ierr)
      End If      
   End Subroutine MEF90DefMechPlasticStrainUpdate
End Module m_MEF90_DefMech
