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
   Private
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
         ! Call MEF90DefMechOperatorDisplacement2D(snesTemp,x,residual,MEF90DefMechCtx,ierr)
      Else If (dim == 3) Then
         ! Call MEF90DefMechOperatorDisplacement3D(snesTemp,x,residual,MEF90DefMechCtx,ierr)
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
         ! Call MEF90DefMechBilinearFormDisplacement2D(snesDispl,x,A,M,flg,MEF90DefMechCtx,ierr)
      Else If (dim == 3) Then
         ! Call MEF90DefMechBilinearFormDisplacement3D(snesDispl,x,A,M,flg,MEF90DefMechCtx,ierr)
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

   Subroutine MEF90DefMechWork(DisplacementVec,MEF90DefMechCtx,work,ierr)
      Type(tVec),Intent(IN)                           :: DisplacementVec
      Type(MEF90DefMechCtx_Type),Intent(IN)           :: MEF90DefMechCtx
      PetscReal,Dimension(:),Pointer                  :: work
      PetscErrorCode,Intent(INOUT)                    :: ierr

      PetscInt                                        :: dim      
      PetscCall(DMGetDimension(MEF90DefMechCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         ! Call MEF90DefMechWork2D(DisplacementVec,MEF90DefMechCtx,work,ierr)
      Else If (dim == 3) Then
         ! Call MEF90DefMechWork3D(DisplacementVec,MEF90DefMechCtx,work,ierr)
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

   Subroutine MEF90DefMechCohesiveEnergy(DisplacementVec,MEF90DefMechCtx,cohesiveEnergy,ierr)
      Type(tVec),Intent(IN)                           :: DisplacementVec
      Type(MEF90DefMechCtx_Type),Intent(IN)           :: MEF90DefMechCtx
      PetscReal,Dimension(:),Pointer                  :: cohesiveEnergy
      PetscErrorCode,Intent(INOUT)                    :: ierr

      PetscInt                                        :: dim      
      PetscCall(DMGetDimension(MEF90DefMechCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         ! Call MEF90DefMechCohesiveEnergy2D(DisplacementVec,MEF90DefMechCtx,cohesiveEnergy,ierr)
      Else If (dim == 3) Then
         ! Call MEF90DefMechCohesiveEnergy3D(DisplacementVec,MEF90DefMechCtx,cohesiveEnergy,ierr)
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

   Subroutine MEF90DefMechElasticEnergy(x,MEF90DefMechCtx,energy,ierr)
      Type(tVec),Intent(IN)                              :: x
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,dimension(:),Pointer                     :: energy
      PetscErrorCode,Intent(INOUT)                       :: ierr

      PetscInt                                           :: dim      
      PetscCall(DMGetDimension(MEF90DefMechCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         ! Call MEF90DefMechElasticEnergy2D(x,MEF90DefMechCtx,energy,ierr)
      Else If (dim == 3) Then
         ! Call MEF90DefMechElasticEnergy3D(x,MEF90DefMechCtx,energy,ierr)
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

   Subroutine MEF90DefMechStress(x,MEF90DefMechCtx,stress,ierr)
      Type(tVec),Intent(IN)                              :: x
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      Type(tVec),Intent(IN)                              :: stress
      PetscErrorCode,Intent(INOUT)                       :: ierr

      PetscInt                                           :: dim      
      PetscCall(DMGetDimension(MEF90DefMechCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         ! Call MEF90DefMechStress2D(x,MEF90DefMechCtx,stress,ierr)
      Else If (dim == 3) Then
         ! Call MEF90DefMechStress3D(x,MEF90DefMechCtx,stress,ierr)
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

   Subroutine MEF90DefMechCrackVolume(x,MEF90DefMechCtx,CrackVolume,ierr)
      Type(tVec),Intent(IN)                              :: x
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,dimension(:),Pointer                     :: CrackVolume
      PetscErrorCode,Intent(INOUT)                       :: ierr

      PetscInt                                           :: dim      
      PetscCall(DMGetDimension(MEF90DefMechCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         ! Call MEF90DefMechCrackVolume2D(x,MEF90DefMechCtx,CrackVolume,ierr)
      Else If (dim == 3) Then
         ! Call MEF90DefMechCrackVolume3D(x,MEF90DefMechCtx,CrackVolume,ierr)
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
         ! Call MEF90DefMechOperatorDamage2D(snesTemp,x,residual,MEF90DefMechCtx,ierr)
      Else If (dim == 3) Then
         ! Call MEF90DefMechOperatorDamage3D(snesTemp,x,residual,MEF90DefMechCtx,ierr)
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
         ! Call MEF90DefMechBilinearFormDamage2D(snesDispl,x,A,M,flg,MEF90DefMechCtx,ierr)
      Else If (dim == 3) Then
         ! Call MEF90DefMechBilinearFormDamage3D(snesDispl,x,A,M,flg,MEF90DefMechCtx,ierr)
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

   Subroutine MEF90DefMechSurfaceEnergy(x,MEF90DefMechCtx,energy,ierr)
      Type(tVec),Intent(IN)                              :: x
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,dimension(:),Pointer                     :: energy
      PetscErrorCode,Intent(INOUT)                       :: ierr

      PetscInt                                           :: dim      
      PetscCall(DMGetDimension(MEF90DefMechCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         ! Call MEF90DefMechSurfaceEnergy2D(x,MEF90DefMechCtx,energy,ierr)
      Else If (dim == 3) Then
         ! Call MEF90DefMechSurfaceEnergy3D(x,MEF90DefMechCtx,energy,ierr)
      End If      
   End Subroutine MEF90DefMechSurfaceEnergy

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechViewEXO"
!!!
!!!  
!!!  MEF90DefMechViewEXO: Save all fields in a MEF90DefMechCtx_Type in an exodus file
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechViewEXO(MEF90DefMechCtx,step,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscInt,Intent(IN)                                :: step
      PetscErrorCode,Intent(INOUT)                       :: ierr

      ! Type(Vec)                                          :: localVec
      ! Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions

      ! Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)

      ! !!! cell-based fields are located before nodal ones in exodus files, so saving them first.
      ! !!! they have no ghost points, so there is no real difference between their local and global vectors
      ! !!!
      ! If (MEF90DefMechGlobalOptions%ForceOffset > 0) Then
      !    If (Associated(MEF90DefMechCtx%Force)) Then
      !       Call VecViewExodusCell(MEF90DefMechCtx%cellDMVect,MEF90DefMechCtx%Force,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
      !                              MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%forceOffset,ierr);CHKERRQ(ierr)
      !    Else
      !       Call PetscPrintf(PETSC_COMM_WORLD,"[WARNING] force field not associated, not saving. Use -force_offset 0 \n",ierr);CHKERRQ(ierr)
      !    End If
      ! End If

      ! If (MEF90DefMechGlobalOptions%pressureForceOffset > 0) Then
      !    If (Associated(MEF90DefMechCtx%pressureForce)) Then
      !       Call VecViewExodusCell(MEF90DefMechCtx%cellDMScal,MEF90DefMechCtx%pressureForce,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
      !                              MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%pressureForceOffset,ierr);CHKERRQ(ierr)
      !    Else
      !       Call PetscPrintf(PETSC_COMM_WORLD,"[WARNING] pressureForce field not associated, not saving. Use -pressureForce_offset 0 \n",ierr);CHKERRQ(ierr)
      !    End If
      ! End If

      ! If (MEF90DefMechGlobalOptions%CrackPressureOffset > 0) Then
      !    If (Associated(MEF90DefMechCtx%CrackPressure)) Then
      !       Call VecViewExodusCell(MEF90DefMechCtx%cellDMScal,MEF90DefMechCtx%CrackPressure,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
      !                              MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%CrackPressureOffset,ierr);CHKERRQ(ierr)
      !    Else
      !       Call PetscPrintf(PETSC_COMM_WORLD,"[WARNING] CrackPressure field not associated, not saving. Use -CrackPressure_offset 0 \n",ierr);CHKERRQ(ierr)
      !    End If
      ! End If

      ! If (MEF90DefMechGlobalOptions%plasticStrainOffset > 0) Then
      !    If (Associated(MEF90DefMechCtx%plasticStrain)) Then
      !       Call VecViewExodusCell(MEF90DefMechCtx%cellDMMatS,MEF90DefMechCtx%plasticStrain,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
      !                              MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%plasticStrainOffset,ierr);CHKERRQ(ierr)
      !    Else
      !       Call PetscPrintf(PETSC_COMM_WORLD,"[WARNING] plasticStrain field not associated, not saving. Use -plasticStrain_offset 0 \n",ierr);CHKERRQ(ierr)
      !    End If
      ! End If


      ! If (MEF90DefMechGlobalOptions%cumulatedPlasticDissipationOffset > 0) Then
      !    If (Associated(MEF90DefMechCtx%cumulatedPlasticDissipation)) Then
      !       Call VecViewExodusCell(MEF90DefMechCtx%cellDMScal,MEF90DefMechCtx%cumulatedPlasticDissipation,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
      !                              MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%cumulatedPlasticDissipationOffset,ierr);CHKERRQ(ierr)
      !    Else
      !       Call PetscPrintf(PETSC_COMM_WORLD,"[WARNING] cumulatedDissipatedPlasticEnergy field not associated, not saving. Use -cumulatedPlasticDissipationOffset 0 \n",ierr);CHKERRQ(ierr)
      !    End If
      ! End If

      ! If (MEF90DefMechGlobalOptions%stressOffset > 0) Then
      !    If (Associated(MEF90DefMechCtx%stress)) Then
      !       Call VecViewExodusCell(MEF90DefMechCtx%cellDMMatS,MEF90DefMechCtx%stress,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
      !                              MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%stressOffset,ierr);CHKERRQ(ierr)
      !    Else
      !       Call PetscPrintf(PETSC_COMM_WORLD,"[WARNING] stress field not associated, not saving. Use -stress_offset 0 \n",ierr);CHKERRQ(ierr)
      !    End If
      ! End If


      ! If ((MEF90DefMechGlobalOptions%boundaryDisplacementOffset > 0) .AND. &
      !     (MEF90DefMechGlobalOptions%boundaryDisplacementOffset /= MEF90DefMechGlobalOptions%displacementOffset)) Then
      !    If (Associated(MEF90DefMechCtx%boundaryDisplacement)) Then
      !       Call DMGetLocalVector(MEF90DefMechCtx%DMVect,localVec,ierr);CHKERRQ(ierr)
      !       Call DMGlobalToLocalBegin(MEF90DefMechCtx%DMVect,MEF90DefMechCtx%boundaryDisplacement,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
      !       Call DMGlobalToLocalEnd(MEF90DefMechCtx%DMVect,MEF90DefMechCtx%boundaryDisplacement,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
      !       Call VecViewExodusVertex(MEF90DefMechCtx%DMVect,localVec,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
      !                                MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%boundaryDisplacementOffset,ierr);CHKERRQ(ierr)
      !       Call DMRestoreLocalVector(MEF90DefMechCtx%DMVect,localVec,ierr);CHKERRQ(ierr)
      !    Else
      !       Call PetscPrintf(PETSC_COMM_WORLD,"[WARNING] boundaryDisplacement field not associated, not saving. Use -boundaryDisplacement_offset 0 \n",ierr);CHKERRQ(ierr)
      !    End If
      ! End If

      ! If (MEF90DefMechGlobalOptions%displacementOffset > 0) Then
      !    If (Associated(MEF90DefMechCtx%displacement)) Then
      !       Call DMGetLocalVector(MEF90DefMechCtx%DMVect,localVec,ierr);CHKERRQ(ierr)
      !       Call DMGlobalToLocalBegin(MEF90DefMechCtx%DMVect,MEF90DefMechCtx%Displacement,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
      !       Call DMGlobalToLocalEnd(MEF90DefMechCtx%DMVect,MEF90DefMechCtx%Displacement,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
      !       Call VecViewExodusVertex(MEF90DefMechCtx%DMVect,localVec,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
      !                                MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%displacementOffset,ierr);CHKERRQ(ierr)
      !       Call DMRestoreLocalVector(MEF90DefMechCtx%DMVect,localVec,ierr);CHKERRQ(ierr)
      !    Else
      !       Call PetscPrintf(PETSC_COMM_WORLD,"[WARNING] displacement field not associated, not saving. Use -displacement_offset 0 \n",ierr);CHKERRQ(ierr)
      !    End If
      ! End If

      ! If (MEF90DefMechGlobalOptions%damageOffset > 0) Then
      !    If (Associated(MEF90DefMechCtx%damage)) Then
      !       Call DMGetLocalVector(MEF90DefMechCtx%DMScal,localVec,ierr);CHKERRQ(ierr)
      !       Call DMGlobalToLocalBegin(MEF90DefMechCtx%DMScal,MEF90DefMechCtx%damage,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
      !       Call DMGlobalToLocalEnd(MEF90DefMechCtx%DMScal,MEF90DefMechCtx%damage,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
      !       Call VecViewExodusVertex(MEF90DefMechCtx%DMScal,localVec,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
      !                                MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%damageOffset,ierr);CHKERRQ(ierr)
      !       Call DMRestoreLocalVector(MEF90DefMechCtx%DMScal,localVec,ierr);CHKERRQ(ierr)
      !    Else
      !       Call PetscPrintf(PETSC_COMM_WORLD,"[WARNING] damage field not associated, not saving. Use -damage_offset 0 \n",ierr);CHKERRQ(ierr)
      !    End If
      ! End If

      ! If ((MEF90DefMechGlobalOptions%boundaryDamageOffset > 0) .AND. &
      !     (MEF90DefMechGlobalOptions%boundaryDamageOffset /= MEF90DefMechGlobalOptions%damageOffset))Then
      !    If (Associated(MEF90DefMechCtx%boundaryDamage)) Then
      !       Call DMGetLocalVector(MEF90DefMechCtx%DMScal,localVec,ierr);CHKERRQ(ierr)
      !       Call DMGlobalToLocalBegin(MEF90DefMechCtx%DMScal,MEF90DefMechCtx%boundaryDamage,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
      !       Call DMGlobalToLocalEnd(MEF90DefMechCtx%DMScal,MEF90DefMechCtx%boundaryDamage,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
      !       Call VecViewExodusVertex(MEF90DefMechCtx%DMScal,localVec,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
      !                                MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%boundaryDamageOffset,ierr);CHKERRQ(ierr)
      !       Call DMRestoreLocalVector(MEF90DefMechCtx%DMScal,localVec,ierr);CHKERRQ(ierr)
      !    Else
      !       Call PetscPrintf(PETSC_COMM_WORLD,"[WARNING] boundaryDamage field not associated, not saving. Use -boundaryDamage_offset 0 \n",ierr);CHKERRQ(ierr)
      !    End If
      ! End If

      ! If (MEF90DefMechGlobalOptions%temperatureOffset > 0) Then
      !    If (Associated(MEF90DefMechCtx%temperature)) Then
      !       Call DMGetLocalVector(MEF90DefMechCtx%DMScal,localVec,ierr);CHKERRQ(ierr)
      !       Call DMGlobalToLocalBegin(MEF90DefMechCtx%DMScal,MEF90DefMechCtx%temperature,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
      !       Call DMGlobalToLocalEnd(MEF90DefMechCtx%DMScal,MEF90DefMechCtx%temperature,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
      !       Call VecViewExodusVertex(MEF90DefMechCtx%DMScal,localVec,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
      !                                MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%temperatureOffset,ierr);CHKERRQ(ierr)
      !       Call DMRestoreLocalVector(MEF90DefMechCtx%DMScal,localVec,ierr);CHKERRQ(ierr)
      !    Else
      !       Call PetscPrintf(PETSC_COMM_WORLD,"[WARNING] Temperature field not associated, not saving. Use -temperature_offset 0 \n",ierr);CHKERRQ(ierr)
      !    End If
      ! End If
      ! If (MEF90DefMechCtx%MEF90Ctx%rank == 0) Then
      !    Call EXUPDA(MEF90DefMechCtx%MEF90Ctx%fileExoUnit,ierr)
      ! End If
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
!!!

   Subroutine MEF90DefMechCreateSNESDisplacement(MEF90DefMechCtx,snesDisp,residual,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      Type(tSNES),Intent(OUT)                            :: snesDisp
      Type(tVec),Intent(IN)                              :: residual
      PetscErrorCode,Intent(INOUT)                       :: ierr
      
      ! Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
      ! Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions
      ! Type(Mat)                                          :: matDisp
      ! Type(SectionReal)                                  :: coordSec
      ! Type(Vec)                                          :: CoordVec
      ! PetscReal,Dimension(:,:),Pointer                   :: CoordPtr
      ! Type(VecScatter)                                   :: ScatterSecToVec
      ! Type(MatNullSpace)                                 :: nspDisp
      ! Type(Vec)                                          :: residualDisp
      ! Type(KSP)                                          :: kspDisp
      ! Type(PC)                                           :: pcDisp
      ! SNESLineSearch                                     :: lsDisp
      ! PetscReal                                          :: atol,rtol,dtol
      ! PetscReal,Dimension(:),Pointer                     :: CoordPCPtr
      ! PetscInt                                           :: dim
      
      ! Call DMMeshGetDimension(MEF90DefMechCtx%DMVect,dim,ierr);CHKERRQ(ierr)
      ! Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)
      ! Call DMCreateMatrix(MEF90DefMechCtx%DMVect,MATAIJ,matDisp,iErr);CHKERRQ(iErr)
      ! Call MatSetOptionsPrefix(matDisp,"Disp_",ierr);CHKERRQ(ierr)
      ! Call MatSetOption(matDisp,MAT_SPD,PETSC_TRUE,ierr);CHKERRQ(ierr)
      ! Call MatSetOption(matDisp,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr);CHKERRQ(ierr)
      ! Call MatSetOption(matDisp,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr);CHKERRQ(ierr)
      ! If (MEF90DefMechGlobalOptions%addDisplacementNullSpace) Then
      !    Call DMMeshGetSectionReal(MEF90DefMechCtx%DMVect,'coordinates',coordSec,ierr);CHKERRQ(ierr)
      !    Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMVect,coordSec,ScatterSecToVec,ierr);CHKERRQ(ierr)
      !    Call DMCreateGlobalVector(MEF90DefMechCtx%DMVect,coordVec,ierr)
      !    Call SectionRealToVec(coordSec,ScatterSecToVec,SCATTER_FORWARD,coordVec,ierr);CHKERRQ(ierr)
      !    Call MatNullSpaceCreateRigidBody(coordVec,nspDisp,ierr);CHKERRQ(ierr)
      !    Call MatSetNearNullSpace(matDisp,nspDisp,ierr);CHKERRQ(ierr)
      !    !!!Call MatSetNullSpace(matDisp,nspDisp,ierr);CHKERRQ(ierr)
      !    Call MatNullSpaceDestroy(nspDisp,ierr);CHKERRQ(ierr)
      !    Call SectionRealDestroy(coordSec,ierr);CHKERRQ(ierr)
      !    Call VecDestroy(coordVec,ierr);CHKERRQ(ierr)
      !    Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
      ! End If

      ! Call MatSetFromOptions(matDisp,ierr);CHKERRQ(ierr)

      ! If (MEF90DefMechGlobalOptions%timeSteppingType == MEF90DefMech_TimeSteppingTypeQuasiStatic) Then
      !    Call SNESCreate(PETSC_COMM_WORLD,snesDisp,ierr);CHKERRQ(ierr)
      !    Call SNESSetApplicationContext(snesDisp,MEF90DefMechCtx,ierr);CHKERRQ(ierr)
      !    Call SNESSetDM(snesDisp,MEF90DefMechCtx%DMVect,ierr);CHKERRQ(ierr)
      !    Call SNESSetOptionsPrefix(snesDisp,'Disp_',ierr);CHKERRQ(ierr)
      !    !Call SNESSetType(snesDisp,SNESKSPONLY,ierr);CHKERRQ(ierr)
      !    Call SNESSetType(snesDisp,SNESLS,ierr);CHKERRQ(ierr)
      !    Call SNESGetSNESLineSearch(snesDisp,lsDisp,ierr);CHKERRQ(ierr)
      !    Call SNESLineSearchSetType(lsDisp,SNESLINESEARCHL2,ierr);CHKERRQ(ierr)
         
      !    Call SNESSetFunction(snesDisp,residual,MEF90DefMechOperatorDisplacement,MEF90DefMechCtx,ierr);CHKERRQ(ierr)
      !    Call SNESSetJacobian(snesDisp,matDisp,matDisp,MEF90DefMechBilinearFormDisplacement,MEF90DefMechCtx,ierr);CHKERRQ(ierr)
      !    atol = 1.0D-7
      !    rtol = 1.0D-5
      !    Call SNESSetTolerances(snesDisp,atol,rtol,PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_INTEGER,PETSC_DEFAULT_INTEGER,ierr);CHKERRQ(ierr)
      !    Call SNESSetFromOptions(snesDisp,ierr);CHKERRQ(ierr)

      !    !!! 
      !    !!! Set some KSP options
      !    !!!
      !    Call SNESGetKSP(snesDisp,kspDisp,ierr);CHKERRQ(ierr)
      !    Call KSPSetType(kspDisp,KSPCG,ierr);CHKERRQ(ierr)
      !    Call KSPSetInitialGuessNonzero(kspDisp,PETSC_TRUE,ierr);CHKERRQ(ierr)
      !    atol = 1.0D-8
      !    rtol = 1.0D-8
      !    dtol = 1.0D+10
      !    Call KSPSetTolerances(kspDisp,rtol,atol,dtol,PETSC_DEFAULT_INTEGER,ierr);CHKERRQ(ierr)
      !    Call KSPSetFromOptions(kspDisp,ierr);CHKERRQ(ierr)
      ! End If

      ! !!! set coordinates in PC for GAMG
      ! !!! For some reason, this makes gamg convergence worse, when the null space is specified.
      ! !!! Will investigate later
      ! Call KSPGetPC(kspDisp,pcDisp,ierr);CHKERRQ(ierr)
      ! !Call DMMeshGetCoordinatesF90(MEF90DefMechCtx%DMVect,coordPtr,ierr);CHKERRQ(ierr)
      ! !Allocate(coordPCPtr(size(CoordPtr)))
      ! !coordPCPtr = reshape(transpose(coordPtr),[size(CoordPtr)])
      ! !coordPCPtr = reshape((coordPtr),[size(CoordPtr)])
      ! !Call PCSetCoordinates(pcDisp,dim,size(coordPtr),coordPCPtr,ierr);CHKERRQ(ierr)
      ! !DeAllocate(coordPCPtr)
      ! !Call DMMeshRestoreCoordinatesF90(MEF90DefMechCtx%DMVect,coordPtr,ierr);CHKERRQ(ierr)
      ! Call PCSetFromOptions(pcDisp,ierr);CHKERRQ(ierr)

      ! Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)
      ! !!! SNESView seems to segfault when PC ml is used, so we setup an additional debug level so we can bypass SNESView 
      ! If (MEF90GlobalOptions%verbose > 1) Then
      !    Call SNESView(snesDisp,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
      ! End If
   End Subroutine MEF90DefMechCreateSNESDisplacement

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCreateSNESDamage"
!!!
!!!  
!!!  MEF90DefMechCreateSNESDamage:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechCreateSNESDamage(MEF90DefMechCtx,snesDamage,residual,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      Type(tSNES),Intent(OUT)                            :: snesDamage
      Type(tVec),Intent(IN)                              :: residual
      PetscErrorCode,Intent(INOUT)                       :: ierr
      
      ! Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
      ! Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions
      ! Type(Mat)                                          :: matDamage
      ! PetscReal,Dimension(:,:),Pointer                   :: CoordPtr
      ! Type(Vec)                                          :: residualDamage
      ! Type(KSP)                                          :: kspDamage
      ! Type(PC)                                           :: pcDamage
      ! SNESLineSearch                                     :: lsDamage
      ! PetscReal                                          :: atol,rtol,dtol
      ! PetscInt                                           :: dim
      ! Type(Vec)                                          :: LB,UB
      
      ! Call DMMeshGetDimension(MEF90DefMechCtx%DMScal,dim,ierr);CHKERRQ(ierr)
      ! Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)
      ! Call DMCreateMatrix(MEF90DefMechCtx%DMScal,MATAIJ,matDamage,iErr);CHKERRQ(iErr)
      ! Call MatSetOptionsPrefix(matDamage,"damage_",ierr);CHKERRQ(ierr)
      ! Call MatSetOption(matDamage,MAT_SPD,PETSC_TRUE,ierr);CHKERRQ(ierr)
      ! Call MatSetOption(matDamage,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr);CHKERRQ(ierr)
      ! Call MatSetOption(matDamage,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr);CHKERRQ(ierr)
      ! Call MatSetFromOptions(matDamage,ierr);CHKERRQ(ierr)

      ! If (MEF90DefMechGlobalOptions%timeSteppingType == MEF90DefMech_TimeSteppingTypeQuasiStatic) Then
      !    Call SNESCreate(PETSC_COMM_WORLD,snesDamage,ierr);CHKERRQ(ierr)
      !    Call SNESSetApplicationContext(snesDamage,MEF90DefMechCtx,ierr);CHKERRQ(ierr)
      !    Call SNESSetDM(snesDamage,MEF90DefMechCtx%DMScal,ierr);CHKERRQ(ierr)
      !    Call SNESSetOptionsPrefix(snesDamage,'damage_',ierr);CHKERRQ(ierr)
      !    Call SNESGetSNESLineSearch(snesDamage,lsDamage,ierr);CHKERRQ(ierr)
      !    Call SNESLineSearchSetType(lsDamage,SNESLINESEARCHL2,ierr);CHKERRQ(ierr)
         
      !    !!! Set default bounds for the damage field
      !    Call DMCreateGlobalVector(MEF90DefMechCtx%DMScal,LB,ierr);CHKERRQ(ierr)
      !    Call VecDuplicate(LB,UB,ierr);CHKERRQ(ierr)
      !    Call VecSet(LB,0.0_Kr,ierr);CHKERRQ(ierr)
      !    Call VecSet(UB,1.0_Kr,ierr);CHKERRQ(ierr)
      !    Call SNESVISetVariableBounds(snesDamage,LB,UB,ierr);CHKERRQ(ierr)


      !    Call SNESSetFunction(snesDamage,residual,MEF90DefMechOperatorDamage,MEF90DefMechCtx,ierr);CHKERRQ(ierr)
      !    Call SNESSetJacobian(snesDamage,matDamage,matDamage,MEF90DefMechBilinearFormDamage,MEF90DefMechCtx,ierr);CHKERRQ(ierr)
      !    atol = 1.0D-7
      !    rtol = 1.0D-5
      !    Call SNESSetTolerances(snesDamage,atol,rtol,PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_INTEGER,PETSC_DEFAULT_INTEGER,ierr);CHKERRQ(ierr)
      !    Call SNESSetFromOptions(snesDamage,ierr);CHKERRQ(ierr)

      !    !!! 
      !    !!! Set some KSP options
      !    !!!
      !    Call SNESGetKSP(snesDamage,kspDamage,ierr);CHKERRQ(ierr)
      !    Call KSPSetType(kspDamage,KSPCG,ierr);CHKERRQ(ierr)
      !    Call KSPSetInitialGuessNonzero(kspDamage,PETSC_TRUE,ierr);CHKERRQ(ierr)
      !    rtol = 1.0D-8
      !    atol = 1.0D-8
      !    dtol = 1.0D+10
      !    Call KSPSetTolerances(kspDamage,rtol,atol,dtol,PETSC_DEFAULT_INTEGER,ierr);CHKERRQ(ierr)
      !    Call KSPSetFromOptions(kspDamage,ierr);CHKERRQ(ierr)
      ! End If
      

      ! !!! set coordinates in PC for GAMG
      ! !!! For some reason, this makes gamg convergence worse, when the null space is specified.
      ! !!! Will investigate later
      ! Call KSPGetPC(kspDamage,pcDamage,ierr);CHKERRQ(ierr)
      ! !Call DMMeshGetCoordinatesF90(MEF90DefMechCtx%DMScal,coordPtr,ierr);CHKERRQ(ierr)
      ! !Allocate(coordPCPtr(size(CoordPtr)))
      ! !coordPCPtr = reshape(transpose(coordPtr),[size(CoordPtr)])
      ! !coordPCPtr = reshape((coordPtr),[size(CoordPtr)])
      ! !Call PCSetCoordinates(pcDamage,dim,size(coordPtr),coordPCPtr,ierr);CHKERRQ(ierr)
      ! !DeAllocate(coordPCPtr)
      ! !Call DMMeshRestoreCoordinatesF90(MEF90DefMechCtx%DMScal,coordPtr,ierr);CHKERRQ(ierr)
      ! Call PCSetFromOptions(pcDamage,ierr);CHKERRQ(ierr)

      ! Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)
      ! !!! SNESView seems to segfault when PC ml is used, so we setup an additional debug level so we can bypass SNESView 
      ! If (MEF90GlobalOptions%verbose > 1) Then
      !    Call SNESView(snesDamage,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
      ! End If
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
