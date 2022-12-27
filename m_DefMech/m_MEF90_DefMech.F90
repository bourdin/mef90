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
      Character(len=MEF90MXSTRLEN)                    :: IOBuffer

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
         Write(IOBuffer,'((A),": linear scaling of damage does not make any sense.\n")') __FUNCT__
         SETERRQ(MEF90DefMechCtx%MEF90Ctx%Comm,PETSC_ERR_ARG_WRONG,IOBuffer)
      Case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%damageLocal,1.0_Kr,ierr))
      End Select

      Select case (MEF90DefMechGlobalOptions%cohesiveDisplacementScaling)
      Case (MEF90Scaling_File)
         PetscCall(MEF90EXOVecLoad(MEF90DefMechCtx%cohesiveDisplacement,MEF90DefMechCtx%cohesiveDisplacementToIOSF,MEF90DefMechCtx%IOToCohesiveDisplacementSF,MEF90DefMechCtx%MEF90Ctx%resultViewer,step,ierr))
      Case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90DefMechCtx%cohesiveDisplacement,time,ierr))
      Case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90DefMechCtx%cohesiveDisplacement,1.0_Kr,ierr))
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
         PetscCall(MEF90VecSetValuesFromOptions(MEF90DefMechCtx%bodyForce,time,ierr))
      Case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90DefMechCtx%bodyForce,1.0_Kr,ierr))
      End Select

      Select case (MEF90DefMechGlobalOptions%boundaryForceScaling)
      Case (MEF90Scaling_File)
         PetscCall(MEF90EXOVecLoad(MEF90DefMechCtx%boundaryForce,MEF90DefMechCtx%boundaryForceToIOSF,MEF90DefMechCtx%IOToBoundaryForceSF,MEF90DefMechCtx%MEF90Ctx%resultViewer,step,ierr))
      Case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90DefMechCtx%boundaryForce,time,ierr))
      Case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90DefMechCtx%boundaryForce,1.0_Kr,ierr))
      End Select

      Select case (MEF90DefMechGlobalOptions%pressureForceScaling)
      Case (MEF90Scaling_File)
         PetscCall(MEF90EXOVecLoad(MEF90DefMechCtx%pressureForce,MEF90DefMechCtx%pressureForceToIOSF,MEF90DefMechCtx%IOToPressureForceSF,MEF90DefMechCtx%MEF90Ctx%resultViewer,step,ierr))
      Case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90DefMechCtx%pressureForce,time,ierr))
      Case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90DefMechCtx%pressureForce,1.0_Kr,ierr))
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

   Subroutine MEF90DefMechBilinearFormDisplacement(snesDispl,x,A,M,MEF90DefMechCtx,ierr)
      Type(tSNES),Intent(IN)                             :: snesDispl
      Type(tVec),Intent(IN)                              :: x
      Type(tMat),Intent(INOUT)                           :: A,M
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(INOUT)                       :: ierr

      PetscInt                                           :: dim      
      PetscCall(DMGetDimension(MEF90DefMechCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         PetscCall(MEF90DefMechBilinearFormDisplacement2D(snesDispl,x,A,M,MEF90DefMechCtx,ierr))
      Else If (dim == 3) Then
         PetscCall(MEF90DefMechBilinearFormDisplacement3D(snesDispl,x,A,M,MEF90DefMechCtx,ierr))
      End If      
   End Subroutine MEF90DefMechBilinearFormDisplacement

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechWork"
!!!
!!!  
!!!  MEF90DefMechWork: wraps calls to MEF90DefMechWork from m_MEF90_DefMechAssembly
!!                       since overloading cannot be used here
!!!  
!!!  (c) 2012-22 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechWork(MEF90DefMechCtx,bodyForceWork,boundaryForceWork,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)           :: MEF90DefMechCtx
      PetscReal,Dimension(:),Pointer                  :: bodyForceWork,boundaryForceWork
      PetscErrorCode,Intent(INOUT)                    :: ierr

      PetscInt                                        :: dim      
      PetscCall(DMGetDimension(MEF90DefMechCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         PetscCall(MEF90DefMechWork2D(MEF90DefMechCtx,bodyForceWork,boundaryForceWork,ierr))
      Else If (dim == 3) Then
         PetscCall(MEF90DefMechWork3D(MEF90DefMechCtx,bodyForceWork,boundaryForceWork,ierr))
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

   Subroutine MEF90DefMechBilinearFormDamage(snesDispl,x,A,M,MEF90DefMechCtx,ierr)
      Type(tSNES),Intent(IN)                             :: snesDispl
      Type(tVec),Intent(IN)                              :: x
      Type(tMat),Intent(INOUT)                           :: A,M
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(INOUT)                       :: ierr

      PetscInt                                           :: dim      
      PetscCall(DMGetDimension(MEF90DefMechCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         PetscCall(MEF90DefMechBilinearFormDamage2D(snesDispl,x,A,M,MEF90DefMechCtx,ierr))
      Else If (dim == 3) Then
         PetscCall(MEF90DefMechBilinearFormDamage3D(snesDispl,x,A,M,MEF90DefMechCtx,ierr))
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
#define __FUNCT__ "MEF90DefMechFormatEXO"
!!!
!!!  
!!!  MEF90DefMechFormatEXO:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022 Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90DefMechFormatEXO(MEF90DefMechCtx,time,ierr)
      Type(MEF90DefMechCtx_Type),Intent(INOUT)           :: MEF90DefMechCtx
      PetscReal,Dimension(:),Pointer                     :: time
      PetscErrorCode,Intent(OUT)                         :: ierr

      Character(len=MXSTLN),Dimension(:),Pointer         :: nameG,nameN,nameC,nameF
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
      PetscInt                                           :: dim,numFields,offset

      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr))
      PetscCall(DMGetDimension(MEF90DefMechCtx%megaDM,dim,ierr))
      Allocate(nameG(0))
      Allocate(nameF(0))

      numFields = 0
      If (MEF90DefMechGlobalOptions%displacementExport) Then
         numFields = numFields+dim
      End If
      If (MEF90DefMechGlobalOptions%damageExport) Then
         numFields = numFields+1
      End If
      If (MEF90DefMechGlobalOptions%temperatureExport) Then
         numFields = numFields+1
      End If

      Allocate(nameN(numFields))
      offset = 1
      If (MEF90DefMechGlobalOptions%displacementExport) Then
         nameN(offset+0) = "Displacement_X"
         nameN(offset+1) = "Displacement_Y"
         If (dim == 3) Then
            nameN(offset+2) = "Displacement_Z"
         End If
         offset = offset + dim
      End If
      If (MEF90DefMechGlobalOptions%damageExport) Then
         nameN(offset) = "Damage"
         offset = offset + 1
      End If
      If (MEF90DefMechGlobalOptions%temperatureExport) Then
         nameN(offset) = "Temperature"
      End If
      
      numFields = 0
      If (MEF90DefMechGlobalOptions%stressExport) Then
         numFields = numFields + dim * (dim+1) / 2
      End If
      If (MEF90DefMechGlobalOptions%plasticStrainExport) Then
         numFields = numFields + dim * (dim+1) / 2
      End If
      If (MEF90DefMechGlobalOptions%cumulatedPlasticDissipationExport) Then
         numFields = numFields + dim * (dim+1) / 2
      End If

      If (MEF90DefMechGlobalOptions%stressExport) Then
         numFields = numFields + dim * (dim+1) / 2
      End If
      If (MEF90DefMechGlobalOptions%plasticStrainExport) Then
         numFields = numFields + dim * (dim+1) / 2
      End If
      If (MEF90DefMechGlobalOptions%cumulatedPlasticDissipationExport) Then
         numFields = 1
      End If
      Allocate(nameC(numFields))

      offset = 1
      If (MEF90DefMechGlobalOptions%stressExport) Then
         If (dim == 2) Then
            nameC(offset+0) = "Stress_XX"
            nameC(offset+1) = "Stress_YY"
            nameC(offset+2) = "Stress_XY"
         Else
            nameC(offset+0) = "Stress_XX"
            nameC(offset+1) = "Stress_YY"
            nameC(offset+2) = "Stress_ZZ"
            nameC(offset+3) = "Stress_YZ"
            nameC(offset+4) = "Stress_XZ"
            nameC(offset+5) = "Stress_XY"
         End If
         offset = offset + dim * (dim+1) / 2
      End If

      If (MEF90DefMechGlobalOptions%plasticStrainExport) Then
         If (dim == 2) Then
            nameC(offset+0) = "PlasticStrain_XX"
            nameC(offset+1) = "PlasticStrain_YY"
            nameC(offset+2) = "PlasticStrain_XY"
         Else
            nameC(offset+0) = "PlasticStrain_XX"
            nameC(offset+1) = "PlasticStrain_YY"
            nameC(offset+2) = "PlasticStrain_ZZ"
            nameC(offset+3) = "PlasticStrain_YZ"
            nameC(offset+4) = "PlasticStrain_XZ"
            nameC(offset+5) = "PlasticStrain_XY"
         End If
         offset = offset + dim * (dim+1) / 2
      End If

      If (MEF90DefMechGlobalOptions%cumulatedPlasticDissipationExport) Then
         nameC(offset) = "CumulatedPlasticDissipation"
      End If

      PetscCallA(MEF90EXOFormat(MEF90DefMechCtx%MEF90Ctx%resultViewer,nameG,nameC,nameN,nameF,time,ierr))
      DeAllocate(nameG)
      DeAllocate(nameN)
      DeAllocate(nameC)
      DeAllocate(nameF)
   End Subroutine MEF90DefMechFormatEXO
   
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

      Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions

      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr))

      If (MEF90DefMechGlobalOptions%displacementExport) Then
         PetscCall(MEF90EXOVecView(MEF90DefMechCtx%displacementLocal,MEF90DefMechCtx%displacementToIOSF,MEF90DefMechCtx%IOToDisplacementSF,MEF90DefMechCtx%MEF90Ctx%resultViewer,step,ierr))
      End If
      If (MEF90DefMechGlobalOptions%damageExport) Then
         PetscCall(MEF90EXOVecView(MEF90DefMechCtx%damageLocal,MEF90DefMechCtx%damageToIOSF,MEF90DefMechCtx%IOToDamageSF,MEF90DefMechCtx%MEF90Ctx%resultViewer,step,ierr))
      End If
      If (MEF90DefMechGlobalOptions%stressExport) Then
         PetscCall(MEF90EXOVecView(MEF90DefMechCtx%stress,MEF90DefMechCtx%stressToIOSF,MEF90DefMechCtx%IOToStressSF,MEF90DefMechCtx%MEF90Ctx%resultViewer,step,ierr))
      End If
      If (MEF90DefMechGlobalOptions%plasticStrainExport) Then
         PetscCall(MEF90EXOVecView(MEF90DefMechCtx%plasticStrain,MEF90DefMechCtx%plasticStrainToIOSF,MEF90DefMechCtx%IOToPlasticStrainSF,MEF90DefMechCtx%MEF90Ctx%resultViewer,step,ierr))
      End If
      If (MEF90DefMechGlobalOptions%cumulatedPlasticDissipationExport) Then
         PetscCall(MEF90EXOVecView(MEF90DefMechCtx%cumulatedPlasticDissipation,MEF90DefMechCtx%cumulatedPlasticDissToIOSF,MEF90DefMechCtx%IOToCumulatedPlasticDissSF,MEF90DefMechCtx%MEF90Ctx%resultViewer,step,ierr))
      End If
   End Subroutine MEF90DefMechViewEXO
   
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCreateSNESDisplacement"
!!!
!!!  
!!!  MEF90DefMechCreateSNESDisplacement:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022 Alexis Marboeuf, marboeua@mcmaster.ca
!!!

   Subroutine MEF90DefMechCreateSNESDisplacement(MEF90DefMechCtx,snesDisplacement,residual,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      Type(tSNES),Intent(OUT)                            :: snesDisplacement
      Type(tVec),Intent(IN)                              :: residual
      PetscErrorCode,Intent(INOUT)                       :: ierr
      
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
      Type(tDM)                                          :: dm
      Type(tMat)                                         :: matDisplacement
      Type(tMatNullSpace)                                :: nspDisplacement
      Type(tKSP)                                         :: kspDisplacement
      Type(tVec)                                         :: gCoord
      PetscReal                                          :: rtol,dtol,atol,stol
      
      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%displacementLocal,dm,ierr))
      PetscCall(DMCreateMatrix(dm,matDisplacement,ierr))
      PetscCall(MatSetOptionsPrefix(matDisplacement,"Displacement_",ierr))
      PetscCall(MatSetOption(matDisplacement,MAT_SPD,PETSC_TRUE,ierr))
      PetscCall(MatSetOption(matDisplacement,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr))
      PetscCall(MatSetOption(matDisplacement,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr))
      !!!
      !!! Set the matrix near null-space consisting of all rigid motions.
      !!!
      PetscCall(DMGetGlobalVector(dm,gCoord,ierr))
      PetscCall(MEF90DefMechProjectCoordinates_Private(gCoord,ierr))
      PetscCall(MatNullSpaceCreateRigidBody(gCoord,nspDisplacement,ierr))
      PetscCall(MatSetNearNullSpace(matDisplacement,nspDisplacement,ierr))
      PetscCall(MatNullSpaceDestroy(nspDisplacement,ierr))
      PetscCall(DMRestoreGlobalVector(dm,gCoord,ierr))
      PetscCall(MatSetFromOptions(matDisplacement,ierr))

      PetscCall(SNESCreate(MEF90DefMechCtx%MEF90Ctx%Comm,snesDisplacement,ierr))
      PetscCall(SNESSetApplicationContext(snesDisplacement,MEF90DefMechCtx,ierr))
      PetscCall(SNESSetDM(snesDisplacement,dm,ierr))
      PetscCall(SNESSetType(snesDisplacement,SNESKSPONLY,ierr))
      PetscCall(SNESSetOptionsPrefix(snesDisplacement,'Displacement_',ierr))

      PetscCall(SNESSetFunction(snesDisplacement,residual,MEF90DefMechOperatorDisplacement,MEF90DefMechCtx,ierr))
      PetscCall(SNESSetJacobian(snesDisplacement,matDisplacement,matDisplacement,MEF90DefMechBilinearFormDisplacement,MEF90DefMechCtx,ierr))
      atol = 1.0D-7
      rtol = 1.0D-5
      stol = 1.0D-7
      PetscCall(SNESSetTolerances(snesDisplacement,atol,rtol,stol,PETSC_DEFAULT_INTEGER,PETSC_DEFAULT_INTEGER,ierr))
      PetscCall(SNESSetFromOptions(snesDisplacement,ierr))
      !!! 
      !!! Set some KSP options
      !!!
      PetscCall(SNESGetKSP(snesDisplacement,kspDisplacement,ierr))
      PetscCall(KSPSetType(kspDisplacement,KSPCG,ierr))
      PetscCall(KSPSetInitialGuessNonzero(kspDisplacement,PETSC_TRUE,ierr))
      atol = 1.0D-8
      rtol = 1.0D-8
      dtol = 1.0D+10
      PetscCall(KSPSetTolerances(kspDisplacement,rtol,atol,dtol,PETSC_DEFAULT_INTEGER,ierr))
      PetscCall(KSPSetFromOptions(kspDisplacement,ierr))
      PetscCall(MatDestroy(matDisplacement,ierr))
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
      Type(tMat)                                         :: matDamage
      Type(tKSP)                                         :: kspDamage
      Type(tVec)                                         :: UB,LB
      PetscReal                                          :: rtol,dtol,atol,stol
      
      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%damageLocal,dm,ierr))
      PetscCall(DMCreateMatrix(dm,matDamage,iErr))
      PetscCall(MatSetOptionsPrefix(matDamage,"Damage_",ierr))
      !!! The matrix is not symmetric if the advection vector is /= 0
      PetscCall(MatSetOption(matDamage,MAT_SPD,PETSC_TRUE,ierr))
      PetscCall(MatSetOption(matDamage,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr))
      PetscCall(MatSetOption(matDamage,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr))
      PetscCall(MatSetFromOptions(matDamage,ierr))

      PetscCall(SNESCreate(MEF90DefMechCtx%MEF90Ctx%Comm,snesDamage,ierr))
      PetscCall(SNESSetApplicationContext(snesDamage,MEF90DefMechCtx,ierr))
      PetscCall(SNESSetDM(snesDamage,dm,ierr))
      PetscCall(SNESSetType(snesDamage,SNESVINEWTONRSLS,ierr))
      PetscCall(SNESSetOptionsPrefix(snesDamage,'Damage_',ierr))

      PetscCall(DMCreateGlobalVector(dm,LB,ierr))
      PetscCall(VecDuplicate(LB,UB,ierr))
      PetscCall(VecSet(LB,0.0_Kr,ierr))
      PetscCall(VecSet(UB,1.0_Kr,ierr))
      PetscCall(SNESVISetVariableBounds(snesDamage,LB,UB,ierr))
      PetscCall(DMRestoreGlobalVector(dm,LB,ierr))
      PetscCall(VecDestroy(UB,ierr))

      PetscCall(SNESSetFunction(snesDamage,residual,MEF90DefMechOperatorDamage,MEF90DefMechCtx,ierr))
      PetscCall(SNESSetJacobian(snesDamage,matDamage,matDamage,MEF90DefMechBilinearFormDamage,MEF90DefMechCtx,ierr))
      atol = 1.0D-7
      rtol = 1.0D-5
      stol = 1.0D-7
      PetscCall(SNESSetTolerances(snesDamage,atol,rtol,stol,PETSC_DEFAULT_INTEGER,PETSC_DEFAULT_INTEGER,ierr))
      PetscCall(SNESSetFromOptions(snesDamage,ierr))
      !!! 
      !!! Set some KSP options
      !!!
      PetscCall(SNESGetKSP(snesDamage,kspDamage,ierr))
      PetscCall(KSPSetType(kspDamage,KSPCG,ierr))
      PetscCall(KSPSetInitialGuessNonzero(kspDamage,PETSC_TRUE,ierr))
      atol = 1.0D-8
      rtol = 1.0D-8
      dtol = 1.0D+10
      PetscCall(KSPSetTolerances(kspDamage,rtol,atol,dtol,PETSC_DEFAULT_INTEGER,ierr))
      PetscCall(KSPSetFromOptions(kspDamage,ierr))
      PetscCall(MatDestroy(matDamage,ierr))
   End Subroutine MEF90DefMechCreateSNESDamage

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechUpdateDamageBounds"
!!!
!!!  
!!!  MEF90DefMechUpdateDamageBounds:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022 Alexis Marboeuf marboeua@mcmaster.ca
!!!

   Subroutine MEF90DefMechUpdateDamageBounds(MEF90DefMechCtx,snesDamage,alpha,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      Type(tSNES),Intent(INOUT)                          :: snesDamage
      Type(tVec),Intent(IN)                              :: alpha
      PetscErrorCode,Intent(INOUT)                       :: ierr
      
      Type(tDM)                                          :: dm
      Type(tVec)                                         :: LB,UB
      PetscReal,Dimension(:),Pointer                     :: LBPtr
      PetscInt                                           :: i
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions

      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr))
      PetscCall(VecGetDM(alpha,dm,ierr))
      PetscCall(DMGetGlobalVector(dm,LB,ierr))
      PetscCall(DMGetGlobalVector(dm,UB,ierr))

      PetscCall(VecSet(UB,1.0_Kr,ierr))
      PetscCall(VecCopy(alpha,LB,ierr))
      If (MEF90DefMechGlobalOptions%irrevthres > 0.0_Kr) Then
         PetscCall(VecGetArrayF90(LB,LBPtr,ierr))
         Do i = 1, size(LBPtr)
            If (LBPtr(i) <= MEF90DefMechGlobalOptions%irrevthres) Then
               LBPtr(i) = 0.0_Kr
            End If
         End Do
         PetscCall(VecRestoreArrayF90(LB,LBPtr,ierr))
      End If
      PetscCall(SNESVISetVariableBounds(snesDamage,LB,UB,ierr))
      PetscCall(DMRestoreGlobalVector(dm,LB,ierr))
      PetscCall(DMRestoreGlobalVector(dm,UB,ierr)) 
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

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechProjectCoordinates_Private"

   Subroutine MEF90DefMechProjectCoordinates_Private(v,ierr)
      Type(tVec),intent(INOUT)           :: v
      PetscErrorCode,intent(INOUT)       :: ierr

      PetscInt                           :: pStart,pEnd,p,numDof,cNumDof,i
      Type(tDM)                          :: dm
      Type(tPetscSection)                :: coordSection,s
      Type(tVec)                         :: coordVec,locV
      PetscScalar,dimension(:),Pointer   :: coordArray,vArray
      PetscScalar,dimension(3)           :: xyz
      PetscInt                           :: dim,pOffset

      PetscCall(VecGetDM(v,dm,ierr))
      PetscCall(DMGetLocalSection(dm,s,ierr))
      PetscCall(PetscSectionGetChart(s,pStart,pEnd,ierr))
      PetscCall(DMGetCoordinateSection(dm,coordSection,ierr))
      PetscCall(DMGetCoordinatesLocal(dm,coordVec,ierr))
      PetscCall(DMGetDimension(dm,dim,ierr))

      PetscCall(DMGetLocalVector(dm,locV,ierr))
      PetscCall(VecGetArrayF90(locV,vArray,ierr))

      Do p = pStart,pEnd-1
         PetscCall(PetscSectionGetDof(s,p,numDof,ierr))
         PetscCall(PetscSectionGetConstraintDof(s,p,cNumDof,ierr))
         If ((numDof > 0) .AND. (cNumDof == 0)) Then
            !!! trick: the coordinate of a point is the average of the coordinates of the points in its closure
            PetscCall(DMPlexVecGetClosure(dm,coordSection,coordVec,p,coordArray,ierr))
            Do i = 1,dim
               xyz(i) = sum(coordArray(i:size(coordArray):dim)) * dim / size(coordArray)
            End Do
            PetscCall(DMPlexVecRestoreClosure(dm,coordSection,coordVec,p,coordArray,ierr))

            PetscCall(PetscSectionGetOffset(s,p,pOffset,ierr))
            Do i = 1,numDof
               vArray(pOffset+i) = xyz(i)
            End Do
         End If
      End Do
      PetscCall(VecRestoreArrayF90(locV,vArray,ierr))
      PetscCall(DMLocalToGlobalBegin(dm,locV,INSERT_VALUES,v,ierr))
      PetscCall(DMLocalToGlobalEnd(dm,locV,INSERT_VALUES,v,ierr))
      PetscCall(DMRestoreLocalVector(dm,locV,ierr))
      !!! Of course, this does not use informations from the section, so it does over-write constrained values
   End Subroutine MEF90DefMechProjectCoordinates_Private
End Module m_MEF90_DefMech
