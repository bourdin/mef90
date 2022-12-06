#include "../MEF90/mef90.inc"
#include "../m_DefMech/mef90DefMech.inc"
Module MEF90_APPEND(m_MEF90_DefMechAssembly,MEF90_DIM)D
#include "petsc/finclude/petsc.h"
#define MEF90_HDRegularization 0.01_Kr

   Use m_MEF90
   Use m_MEF90_DefMechCtx
   Use MEF90_APPEND(m_MEF90_DefMechSplit,MEF90_DIM)D
   Use m_MEF90_DefMechAT
   
   Implicit none
   Private
   Public MEF90DefMechOperatorDisplacement,     &
          MEF90DefMechBilinearFormDisplacement, &
          MEF90DefMechWork,                     &
          MEF90DefMechCohesiveEnergy,           &
          MEF90DefMechPlasticDissipation,       &
          MEF90DefMechElasticEnergy,            &
          MEF90DefMechOperatorDamage,           &
          MEF90DefMechBilinearFormDamage,       &
          MEF90DefMechSurfaceEnergy,            &
          MEF90DefMechCrackVolume,              &
          MEF90DefMechStress

Contains

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDisplacement"
!!!
!!!  
!!!  MEF90DefMechOperatorDisplacement: Build the operator. When called in SNES, the solution time should always match the target time, 
!!!                                    so there is no need for interpolation of the forcees, external, and boundary values
!!!  
!!!  (c) 2012-20 Blaise Bourdin bourdin@lsu.edu, Erwan Tanne erwan.tanne@gmail.com
!!!         2022 Blaise Bourdin bourdin@mcmaster.ca
!!!         2022 Alexis Marboeuf, marboeua@mcmaster.ca
!!!

   Subroutine MEF90DefMechOperatorDisplacement(snesDisplacement,displacement,residual,MEF90DefMechCtx,ierr)
      Type(tSNES),Intent(IN)                             :: snesDisplacement
      Type(tVec),Intent(IN)                              :: displacement
      Type(tVec),Intent(INOUT)                           :: residual
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(INOUT)                       :: ierr

      Type(tDM)                                          :: dmDisplacement,dmDamage,dmTemperature,dmCohesiveDisplacement,dmBodyForce,dmBoundaryForce,dmPressureForce,dmPlasticStrain
      Type(tPetscSection)                                :: sectionBodyForce,sectionBoundaryForce,sectionPressureForce,sectionPlasticStrain
      PetscReal,Dimension(:),Pointer                     :: cohesiveDisplacementDof,displacementDof,damageDof,temperatureDof
      PetscReal,Dimension(:),Pointer                     :: bodyForceArray,boundaryForceArray,pressureForceArray,plasticStrainArray
      Type(MEF90_MATS)                                   :: plasticStrainCell
      Type(tIS)                                          :: setIS,setPointIS
      PetscInt,dimension(:),Pointer                      :: setID,setPointID
      PetscInt                                           :: set,QuadratureOrder,cell,dim,vecOffset
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90DefMechFaceSetOptions_Type),Pointer      :: faceSetOptions
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemVect
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemScal
      DMPolytopeType                                     :: cellGeometry
      Type(MEF90ElementType)                             :: elemVectType,elemScalType
      PetscReal,Dimension(:),Pointer                     :: residualDof
      
      Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90CtxGlobalOptions
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
      Type(tVec)                                         :: locResidual
      Class(MEF90DefMechAT_Type),Allocatable             :: ATModel
      Class(MEF90_DEFMECHSPLIT),Allocatable              :: Split
      PetscBool                                          :: cellIsElastic
      Type(MEF90_MATS)                                   :: inelasticStrainGauss,stressGaussPlus,stressGaussMinus,stressGauss
      Type(MEF90_VECT)                                   :: U0Gauss,bodyForce,boundaryForce,pressureForce
      PetscReal                                          :: damageGauss,temperatureGauss
      PetscInt                                           :: iDof,iGauss,numDofDisplacement,numDofDamage,numGauss
    
      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr))
      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%displacementLocal,dmDisplacement,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%cohesiveDisplacement,dmCohesiveDisplacement,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%damageLocal,dmDamage,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%TemperatureLocal,dmTemperature,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%bodyForce,dmBodyForce,ierr))
      PetscCall(DMGetLocalSection(dmBodyForce,sectionBodyForce,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%plasticStrain,dmPlasticStrain,ierr))
      PetscCall(DMGetLocalSection(dmPlasticStrain,sectionPlasticStrain,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%boundaryForce,dmBoundaryForce,ierr))
      PetscCall(DMGetLocalSection(dmBoundaryForce,sectionBoundaryForce,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%pressureForce,dmPressureForce,ierr))
      PetscCall(DMGetLocalSection(dmPressureForce,sectionPressureForce,ierr))

      PetscCall(DMGetDimension(dmDisplacement,dim,ierr))

      PetscCall(DMGetLocalVector(dmDisplacement,locResidual,ierr))
      PetscCall(DMGlobalToLocal(dmDisplacement,displacement,INSERT_VALUES,MEF90DefMechCtx%displacementLocal,ierr))

      PetscCall(VecSet(residual,0.0_Kr,ierr))
      PetscCall(VecSet(locResidual,0.0_Kr,ierr))

      !!! get IS for cell sets
      PetscCall(DMGetLabelIdIS(dmDisplacement,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm,setIS,ierr))

      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(VecGetArrayF90(MEF90DefMechCtx%plasticStrain,plasticStrainArray,ierr))
         PetscCall(VecGetArrayF90(MEF90DefMechCtx%bodyForce,bodyForceArray,ierr))
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         Do set = 1,size(setID)
            PetscCall(DMGetStratumIS(dmDisplacement,MEF90CellSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
               PetscCall(PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr))
               PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr))
         
               PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(DMPlexGetCellType(dmDisplacement,setPointID(1),cellGeometry,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellGeometry,elemVectType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellGeometry,elemScalType,ierr))

               !!! get the ATModel and split objects
               PetscCall(MEF90DefMechGetATModel(cellSetOptions,ATModel,cellIsElastic,ierr))
               If (cellSetOptions%unilateralContactHybrid) Then
                  Split = MEF90_DEFMECHSPLITNONE()
               Else
                  PetscCall(MEF90DefMechGetSplit(cellSetOptions,Split,ierr))
               End If
               !!! Allocate elements 
               QuadratureOrder = max(2*elemVectType%order,split%damageOrder + split%strainOrder)
               PetscCall(MEF90ElementCreate(dmDisplacement,setPointIS,elemVect,QuadratureOrder,elemVectType,ierr))
               PetscCall(MEF90ElementCreate(dmDamage,setPointIS,elemScal,QuadratureOrder,elemScalType,ierr))

               numDofDisplacement = size(elemVect(1)%BF(:,1))
               numDofDamage = size(elemScal(1)%BF(:,1))
               numGauss = size(elemVect(1)%Gauss_C)

               Allocate(residualDof(numDofDisplacement))

               Do cell = 1,size(setPointID)
                  residualDof = 0.0_Kr
                  Do iGauss = 1,numGauss
                     !!! Main term: [a(\alpha) sigma^+(u) + sigma^-(u)] . e(v)
                     damageGauss = 0.0_Kr
                     If (.NOT. cellIsElastic) Then
                        PetscCall(DMPlexVecGetClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                           Do iDof = 1,numDofDamage
                              damageGauss = damageGauss + damageDof(iDof) * elemScal(cell)%BF(iDof,iGauss)
                           End Do ! iDof numDofDamage
                        PetscCall(DMPlexVecRestoreClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                     End If

                     inelasticStrainGauss = 0.0_Kr
                     PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                     Do iDof = 1,numDofDisplacement
                        inelasticStrainGauss = inelasticStrainGauss + displacementDof(iDof) * elemVect(cell)%GradS_BF(iDof,iGauss)
                     End Do ! iDof numDofDisplacement

                     PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))

                     temperatureGauss = 0.0_Kr
                     PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,MEF90DefMechCtx%TemperatureLocal,setPointID(cell),temperatureDof,ierr))
                     Do iDof = 1,numDofDamage
                        temperatureGauss = temperatureGauss + temperatureDof(iDof) * elemScal(cell)%BF(iDof,iGauss)
                     End Do ! iDof numDofDamage

                     PetscCall(DMPlexVecRestoreClosure(dmTemperature,PETSC_NULL_SECTION,MEF90DefMechCtx%TemperatureLocal,setPointID(cell),temperatureDof,ierr))
                     inelasticStrainGauss = inelasticStrainGauss - (temperatureGauss * matpropSet%LinearThermalExpansion)

   ! #if MEF90_DIM == 2
   !!! We need something along these lines
   !!! Adding terms in planestrain for plasticity with tr(p) = 0
   ! If (.NOT. matProp%HookesLaw%isPlaneStress) Then
   !    stressGauss = stressGauss +  stiffness * ( matProp%HookesLaw%lambda*trace(plasticStrainCell)*MEF90MatS2DIdentity )
   ! End If
   ! #endif

                     PetscCall(PetscSectionGetOffset(sectionPlasticStrain,setPointID(cell),vecOffset,ierr))
                     plasticStrainCell = plasticStrainArray(vecOffset:vecOffset+SIZEOFMEF90_MATS)

                     Call Split%DEED(inelasticStrainGauss - plasticStrainCell,matpropSet%HookesLaw,stressGaussPlus,stressGaussMinus)
                     If (cellIsElastic) Then
                        stressGauss = stressGaussPlus + stressGaussMinus
                     Else
                        stressGauss = (ATModel%a(damageGauss) + matpropSet%residualStiffness) * stressGaussPlus + stressGaussMinus
                     End If

                     Do iDof = 1,numDofDisplacement
                        residualDof(iDof) = residualDof(iDof) + elemVect(cell)%Gauss_C(iGauss) * (stressGauss .DotP. elemVect(cell)%GradS_BF(iDof,iGauss))
                     End Do ! iDof numDofDisplacement
                  End Do ! iGauss
                  PetscCall(DMPlexVecSetClosure(dmDisplacement,PETSC_NULL_SECTION,locResidual,setPointID(cell),residualDof,ADD_VALUES,ierr))
               End Do ! cell

               If (norm2(cellSetOptions%bodyForce) /= 0.0_Kr) Then
                  Do cell = 1,size(setPointID)
                     residualDof = 0.0_Kr
                     !!! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !!! If this happens, we will need to protect this loop
                     PetscCall(PetscSectionGetOffset(sectionBodyForce,setPointID(cell),vecOffset,ierr))
                     bodyForce = bodyForceArray(vecOffset:vecOffset+SIZEOFMEF90_VECT)
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDisplacement
                           residualDof(iDof) = residualDof(iDof) - elemVect(cell)%Gauss_C(iGauss) * (bodyForce .DotP. elemVect(cell)%BF(iDof,iGauss))
                        End Do ! iDof numDofDisplacement
                     End Do ! iGauss
                     PetscCall(DMPlexVecSetClosure(dmDisplacement,PETSC_NULL_SECTION,locResidual,setPointID(cell),residualDof,ADD_VALUES,ierr))
                  End Do ! cell
               End If ! cellSetOptions%bodyForce

               ! !!! Damping
               ! If (MEF90DefMechGlobalOptions%dampingCoefficientDisplacement * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
               !    Do cell = 1,size(setPointID)
               !       residualDof = 0.0_Kr
               !       Do iGauss = 1,numGauss
               !          displacementDampingGauss = 0.0_Kr
               !          PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,locDisplacement,setPointID(cell),displacementDof,ierr))
               !          PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementPreviousStepLocal,setPointID(cell),displacementPreviousStepDof,ierr))
               !          Do iDof = 1,numDofDisplacement
               !             displacementDampingGauss = displacementDampingGauss + ((displacementDof(iDof) - displacementPreviousStepDof(iDof)) * elemVect(cell)%BF(iDof,iGauss))
               !          End Do ! iDof numDofDisplacement
               !          PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,locDisplacement,setPointID(cell),displacementDof,ierr))
               !          PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementPreviousStepLocal,setPointID(cell),displacementPreviousStepDof,ierr))
               !          displacementDampingGauss = displacementDampingGauss * MEF90DefMechGlobalOptions%dampingCoefficientDisplacement / MEF90DefMechCtx%timeStep
               !          Do iDof = 1,numDofDisplacement
               !             residualDof(iDof) = residualDof(iDof) + elemVect(cell)%Gauss_C(iGauss) * &
               !                                  ( displacementDampingGauss .DotP. elemVect(cell)%BF(iDof,iGauss) )
               !          End Do ! iDof numDofDisplacement
               !       End Do ! iGauss
               !       PetscCall(DMPlexVecSetClosure(dmDisplacement,PETSC_NULL_SECTION,locResidual,setPointID(cell),residualDof,ADD_VALUES,ierr))
               !    End Do ! cell
               ! End If ! damping

               !!! Cohesive force
               If ((norm2(cellSetOptions%cohesiveDisplacement) /= 0.0_Kr) .AND. (matpropSet%cohesiveStiffness /= 0.0_Kr)) Then
                  Do cell = 1,size(setPointID)
                     residualDof = 0.0_Kr
                     U0Gauss = 0.0_Kr
                     PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                     PetscCall(DMPlexVecGetClosure(dmCohesiveDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%cohesiveDisplacement,setPointID(cell),cohesiveDisplacementDof,ierr))
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDisplacement
                           U0Gauss = U0Gauss + (displacementDof(iDof) - cohesiveDisplacementDof(iDof)) * elemVect(cell)%BF(iDof,iGauss)
                        End Do ! iDof numDofDisplacement
                     End Do ! iGauss

                     PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                     PetscCall(DMPlexVecRestoreClosure(dmCohesiveDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%cohesiveDisplacement,setPointID(cell),cohesiveDisplacementDof,ierr))
                     U0Gauss = U0Gauss * matpropSet%cohesiveStiffness
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDisplacement
                           residualDof(iDof) = residualDof(iDof) + elemVect(cell)%Gauss_C(iGauss) * (U0Gauss .DotP. elemVect(cell)%BF(iDof,iGauss))
                        End Do ! iDof numDofDisplacement
                     End Do ! iGauss
                     PetscCall(DMPlexVecSetClosure(dmDisplacement,PETSC_NULL_SECTION,locResidual,setPointID(cell),residualDof,ADD_VALUES,ierr))
                  End Do ! cell
               End If ! cohesiveDisplacement

               ! !!! crack pressure
               ! If (cellSetOptions%crackPressure /= 0.0_Kr)) Then
               !    Do cell = 1,size(setPointID)
               !       residualDof = 0.0_Kr
               !       Do iGauss = 1,numGauss
               !          gradientDamageGauss = 0.0_Kr
               !          PetscCall(DMPlexVecGetClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
               !          Do iDof = 1,numDofDamage
               !             GradientDamageGauss = GradientDamageGauss + damageDof(iDof) * elemScal(cell)%Grad_BF(iDof,iGauss) 
               !          End Do ! iDof numDofDamage
               !          PetscCall(DMPlexVecRestoreClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
               !          CrackPressureCell = 0.0_Kr
               !          PetscCall(DMPlexVecGetClosure(dmCrackPressure,PETSC_NULL_SECTION,MEF90DefMechCtx%crackPressure,setPointID(cell),crackPressureDof,ierr))
               !          Do iDof = 1,numDofDamage
               !             CrackPressureCell = CrackPressureCell + crackPressureDof(iDof) * elemScal(cell)%Grad_BF(iDof,iGauss) 
               !          End Do ! iDof numDofDamage
               !          PetscCall(DMPlexVecRestoreClosure(dmCrackPressure,PETSC_NULL_SECTION,MEF90DefMechCtx%crackPressure,setPointID(cell),crackPressureDof,ierr))
               !          Do iDof = 1,numDofDisplacement
               !             residualDof(iDof) = residualDof(iDof) + elemVect(cell)%Gauss_C(iGauss) * &
               !                                  CrackPressureCell * (GradientDamageGauss  .DotP. elemVect(cell)%BF(iDof,iGauss))
               !          End Do ! iDof numDofDisplacement
               !       End Do ! iGauss
               !       PetscCall(DMPlexVecSetClosure(dmDisplacement,PETSC_NULL_SECTION,locResidual,setPointID(cell),residualDof,ADD_VALUES,ierr))
               !    End Do ! cell
               ! End If ! crack Pressure

               DeAllocate(residualDof)

               PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(MEF90ElementDestroy(elemVect,ierr))
               PetscCall(MEF90ElementDestroy(elemScal,ierr))
            End If ! setPointIS
            PetscCall(ISDestroy(setPointIS,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
         PetscCall(VecRestoreArrayF90(MEF90DefMechCtx%plasticStrain,plasticStrainArray,ierr))
         PetscCall(VecRestoreArrayF90(MEF90DefMechCtx%bodyForce,bodyForceArray,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))

      !!! face-based contributions
      PetscCall(DMGetLabelIdIS(dmDisplacement,MEF90FaceSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm,setIS,ierr))
      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         PetscCall(VecGetArrayF90(MEF90DefMechCtx%boundaryForce,boundaryForceArray,ierr))
         PetscCall(VecGetArrayF90(MEF90DefMechCtx%pressureForce,pressureForceArray,ierr))
         Do set = 1,size(setID)
            PetscCall(DMGetStratumIS(dmDisplacement,MEF90FaceSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
              PetscCall(PetscBagGetDataMEF90DefMechCtxFaceSetOptions(MEF90DefMechCtx%FaceSetOptionsBag(set),faceSetOptions,ierr))
              PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
              PetscCall(DMPlexGetCellType(dmDisplacement,setPointID(1),cellGeometry,ierr))
              PetscCall(MEF90ElementGetTypeBoundary(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellGeometry,elemVectType,ierr))
              PetscCall(MEF90ElementGetTypeBoundary(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellGeometry,elemScalType,ierr))

              !!! Allocate elements 
              QuadratureOrder = 2*elemVectType%order
              PetscCall(MEF90ElementCreate(dmDisplacement,setPointIS,elemVect,QuadratureOrder,elemVectType,ierr))
              PetscCall(MEF90ElementCreate(dmDamage,setPointIS,elemScal,QuadratureOrder,elemScalType,ierr))

              numDofDisplacement = size(elemVect(1)%BF(:,1))
              numDofDamage = size(elemScal(1)%BF(:,1))
              numGauss = size(elemVect(1)%Gauss_C)

              Allocate(residualDof(numDofDisplacement))

               If (norm2(faceSetOptions%boundaryForce) /= 0.0_Kr) Then
                  Do cell = 1,size(setPointID)
                     residualDof = 0.0_Kr
                     PetscCall(PetscSectionGetOffset(sectionBoundaryForce,setPointID(cell),vecOffset,ierr))
                     boundaryForce = boundaryForceArray(vecOffset:vecOffset+SIZEOFMEF90_VECT)
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDisplacement
                           residualDof(iDof) = residualDof(iDof) - elemVect(cell)%Gauss_C(iGauss) * (boundaryForce .DotP. elemVect(cell)%BF(iDof,iGauss))
                        End Do ! iDof numDofDisplacement
                     End Do ! iGauss
                     PetscCall(DMPlexVecSetClosure(dmDisplacement,PETSC_NULL_SECTION,locResidual,setPointID(cell),residualDof,ADD_VALUES,ierr))
                  End Do ! cell
               End If ! faceSetOptions%boundaryForce

               If (faceSetOptions%pressureForce /= 0.0_Kr) Then
                  Do cell = 1,size(setPointID)
                     residualDof = 0.0_Kr
                     PetscCall(PetscSectionGetOffset(sectionPressureForce,setPointID(cell),vecOffset,ierr))
                     pressureForce = pressureForceArray(vecOffset)*elemVect(cell)%outerNormal
                     Do iGauss = 1,numGauss  
                        Do iDof = 1,numDofDisplacement
                           residualDof(iDof) = residualDof(iDof) - elemVect(cell)%Gauss_C(iGauss) * (pressureForce .DotP. elemVect(cell)%BF(iDof,iGauss))
                        End Do ! iDof numDofDisplacement
                     End Do ! iGauss
                     PetscCall(DMPlexVecSetClosure(dmDisplacement,PETSC_NULL_SECTION,locResidual,setPointID(cell),residualDof,ADD_VALUES,ierr))
                  End Do ! cell
               End If ! faceSetOptions%boundaryForce

               DeAllocate(residualDof)
               PetscCall(MEF90ElementDestroy(elemVect,ierr))
               PetscCall(MEF90ElementDestroy(elemScal,ierr))
               PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
            End If ! setPointIS
            PetscCall(ISDestroy(setPointIS,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
         PetscCall(VecRestoreArrayF90(MEF90DefMechCtx%boundaryForce,boundaryForceArray,ierr))
         PetscCall(VecRestoreArrayF90(MEF90DefMechCtx%pressureForce,pressureForceArray,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))
      PetscCall(DMLocalToGlobalBegin(dmDisplacement,locResidual,ADD_VALUES,residual,ierr))
      PetscCall(DMLocalToGlobalEnd(dmDisplacement,locResidual,ADD_VALUES,residual,ierr))
      PetscCall(DMRestoreLocalVector(dmDisplacement,locResidual,ierr))
      !PetscCall(DMRestoreLocalVector(dmDisplacement,locDisplacement,ierr))
   End Subroutine MEF90DefMechOperatorDisplacement

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDisplacement"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDisplacement:
!!!  
!!!  (c) 2012-20 Blaise Bourdin bourdin@lsu.edu,Erwan Tanne erwan.tanne@gmail.com
!!!         2022 Blaise Bourdin bourdin@mcmaster.ca
!!!         2022 Alexis Marboeuf, marboeua@mcmaster.ca
!!!

   Subroutine MEF90DefMechBilinearFormDisplacement(snesDisplacement,displacement,A,M,MEF90DefMechCtx,ierr)
      Type(tSNES),Intent(IN)                             :: snesDisplacement
      Type(tVec),Intent(IN)                              :: displacement
      Type(tMat),Intent(INOUT)                           :: A,M
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(INOUT)                       :: ierr

      Type(tDM)                                          :: dmDisplacement,dmDamage,dmTemperature,dmPlasticStrain
      Type(tPetscSection)                                :: sectionPlasticStrain
      Type(tIS)                                          :: setIS,setPointIS
      PetscInt                                           :: numDofDisplacement,numDofDamage,numGauss,set,cell,iGauss,iDof,jDof,QuadratureOrder,dim,vecOffset
      PetscReal                                          :: damageGauss,temperatureGauss
      PetscBool                                          :: cellIsElastic
      PetscInt,dimension(:),Pointer                      :: setID,setPointID
      PetscReal,Dimension(:),Pointer                     :: matDof,displacementDof,damageDof,temperatureDof,plasticStrainArray
      Type(MEF90_HOOKESLAW)                              :: AGaussPlus,AGaussMinus,AGauss
      Type(MEF90_MATS)                                   :: inelasticStrainGauss,plasticStrainCell,AGradS_BF
      Type(MEF90_VECT)                                   :: U0Gauss
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemVect
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemScal
      DMPolytopeType                                     :: cellGeometry
      Type(MEF90ElementType)                             :: elemVectType,elemScalType
      Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90CtxGlobalOptions
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
      Class(MEF90DefMechAT_Type),Allocatable             :: ATModel
      Class(MEF90_DEFMECHSPLIT),Allocatable              :: Split

      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr))
      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%displacementLocal,dmDisplacement,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%damageLocal,dmDamage,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%temperatureLocal,dmTemperature,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%plasticStrain,dmPlasticStrain,ierr))
      PetscCall(DMGetLocalSection(dmPlasticStrain,sectionPlasticStrain,ierr))

      PetscCall(DMGetDimension(dmDisplacement,dim,ierr))

      PetscCall(DMGlobalToLocal(dmDisplacement,displacement,INSERT_VALUES,MEF90DefMechCtx%displacementLocal,ierr))

      PetscCall(MatZeroEntries(A,ierr))

      !!! get IS for cell sets
      PetscCall(DMGetLabelIdIS(dmDisplacement,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm,setIS,ierr))

      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         PetscCall(VecGetArrayF90(MEF90DefMechCtx%plasticStrain,plasticStrainArray,ierr))
         Do set = 1,size(setID)
            PetscCall(DMGetStratumIS(dmDisplacement,MEF90CellSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
               PetscCall(PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr))
               PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr))
         
               PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(DMPlexGetCellType(dmDisplacement,setPointID(1),cellGeometry,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellGeometry,elemVectType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellGeometry,elemScalType,ierr))

               !!! get the ATModel and split objects
               PetscCall(MEF90DefMechGetATModel(cellSetOptions,ATModel,cellIsElastic,ierr))
               If (cellSetOptions%unilateralContactHybrid) Then
                  Split = MEF90_DEFMECHSPLITNONE()
               Else
                  PetscCall(MEF90DefMechGetSplit(cellSetOptions,Split,ierr))
               End If
               !!! Allocate elements 
               QuadratureOrder = max(2*elemVectType%order,split%damageOrder + split%strainOrder)
               PetscCall(MEF90ElementCreate(dmDisplacement,setPointIS,elemVect,QuadratureOrder,elemVectType,ierr))
               PetscCall(MEF90ElementCreate(dmDamage,setPointIS,elemScal,QuadratureOrder,elemScalType,ierr))

               numDofDisplacement = size(elemVect(1)%BF(:,1))
               numDofDamage = size(elemScal(1)%BF(:,1))
               numGauss = size(elemVect(1)%Gauss_C)

               Allocate(matDof(numDofDisplacement*numDofDisplacement))

               Do cell = 1,size(setPointID)
                  matDof = 0.0_Kr
                  Do iGauss = 1,numGauss
                     !!! Main term: [a(\alpha) sigma^+(u) + sigma^-(u)] . e(v)
                     damageGauss = 0.0_Kr
                     If (.NOT. cellIsElastic) Then
                        PetscCall(DMPlexVecGetClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                        Do iDof = 1,numDofDamage
                           damageGauss = damageGauss + damageDof(iDof) * elemScal(cell)%BF(iDof,iGauss)
                        End Do ! iDof numDofDamage
                        PetscCall(DMPlexVecRestoreClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                     End If

                     inelasticStrainGauss = 0.0_Kr
                     PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                     Do iDof = 1,numDofDisplacement
                        inelasticStrainGauss = inelasticStrainGauss + displacementDof(iDof) * elemVect(cell)%GradS_BF(iDof,iGauss)
                     End Do ! iDof numDofDisplacement
                     PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))

                     temperatureGauss = 0.0_Kr
                     PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,MEF90DefMechCtx%temperatureLocal,setPointID(cell),temperatureDof,ierr))
                     Do iDof = 1,numDofDamage
                        temperatureGauss = temperatureGauss + temperatureDof(iDof) * elemScal(cell)%BF(iDof,iGauss)
                     End Do ! iDof numDofDamage
                     PetscCall(DMPlexVecRestoreClosure(dmTemperature,PETSC_NULL_SECTION,MEF90DefMechCtx%temperatureLocal,setPointID(cell),temperatureDof,ierr))
                     inelasticStrainGauss = inelasticStrainGauss - (temperatureGauss * matpropSet%LinearThermalExpansion)

! #if MEF90_DIM == 2
!!! We need something along these lines
!!! Adding terms in planestrain for plasticity with tr(p) = 0
! If (.NOT. matProp%HookesLaw%isPlaneStress) Then
!    stressGauss = stressGauss +  stiffness * ( matProp%HookesLaw%lambda*trace(plasticStrainCell)*MEF90MatS2DIdentity )
! End If
! #endif

                     PetscCall(PetscSectionGetOffset(sectionPlasticStrain,setPointID(cell),vecOffset,ierr))
                     plasticStrainCell = plasticStrainArray(vecOffset:vecOffset+SIZEOFMEF90_MATS)

                     If (cellIsElastic) Then
                        AGauss = matpropSet%HookesLaw
                     Else
                        Call Split%D2EED(inelasticStrainGauss - plasticStrainCell,matpropSet%HookesLaw,AGaussPlus,AGaussMinus)
                        AGauss = (ATModel%a(damageGauss) + matpropSet%residualStiffness) * AGaussPlus + AGaussMinus
                     End If
                     Do jDof = 0,numDofDisplacement-1
                        Do iDof = 1,numDofDisplacement
                           AGradS_BF = AGauss * elemVect(cell)%GradS_BF(iDof,iGauss)
                           matDof(jDof*numDofDisplacement+iDof) = matDof(jDof*numDofDisplacement+iDof) + elemVect(cell)%Gauss_C(iGauss) * (AGradS_BF .DotP. elemVect(cell)%GradS_BF(jDof+1,iGauss))
                        End Do ! jDof numDofDisplacement
                     End Do ! iDof numDofDisplacement
                  End Do ! iGauss
                  PetscCall(DMPlexMatSetClosure(dmDisplacement,PETSC_NULL_SECTION,PETSC_NULL_SECTION,A,setPointID(cell),matDof,ADD_VALUES,ierr))
               End Do ! cell

               ! !!! Damping
               ! If (MEF90DefMechGlobalOptions%dampingCoefficientDisplacement * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
               !    Do cell = 1,size(setPointID)
               !       matDof = 0.0_Kr
               !       Do iGauss = 1,numGauss
               !          Do jDof = 0,numDofDisplacement-1
               !             Do iDof = 1,numDofDisplacement
               !                U0Gauss = elemVect(cell)%BF(iDof,iGauss) * MEF90DefMechGlobalOptions%dampingCoefficientDisplacement / MEF90DefMechCtx%timeStep
               !                matDof(jDof*numDofDisplacement+iDof) = matDof(jDof*numDofDisplacement+iDof) + elemVect(cell)%Gauss_C(iGauss) * (U0Gauss .DotP. elemVect(cell)%BF(jDof+1,iGauss))
               !             End Do ! iDof numDofDisplacement
               !          End Do ! jDof numDofDisplacement
               !       End Do ! iGauss
               !       PetscCall(DMPlexMatSetClosure(dmDisplacement,PETSC_NULL_SECTION,PETSC_NULL_SECTION,A,setPointID(cell),matDof,ADD_VALUES,ierr))
               !    End Do ! cell
               ! End If ! damping

               !!! Cohesive force
               If ((norm2(cellSetOptions%cohesiveDisplacement) /= 0.0_Kr) .AND. (matpropSet%cohesiveStiffness /= 0.0_Kr)) Then
                  Do cell = 1,size(setPointID)
                     matDof = 0.0_Kr
                     Do iGauss = 1,numGauss
                        Do jDof = 0,numDofDisplacement-1
                           Do iDof = 1,numDofDisplacement
                              U0Gauss = matpropSet%cohesiveStiffness * elemVect(cell)%BF(iDof,iGauss)
                              matDof(jDof*numDofDisplacement+iDof) = matDof(jDof*numDofDisplacement+iDof) + elemVect(cell)%Gauss_C(iGauss) * (U0Gauss .DotP. elemVect(cell)%BF(jDof+1,iGauss))
                           End Do ! iDof numDofDisplacement
                        End Do ! jDof numDofDisplacement
                     End Do ! iGauss
                     PetscCall(DMPlexMatSetClosure(dmDisplacement,PETSC_NULL_SECTION,PETSC_NULL_SECTION,A,setPointID(cell),matDof,ADD_VALUES,ierr))
                  End Do ! cell
               End If ! cohesiveDisplacement

               DeAllocate(matDof)

               PetscCall(MEF90ElementDestroy(elemVect,ierr))
               PetscCall(MEF90ElementDestroy(elemScal,ierr))
               PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
            End If ! setPointIS
            PetscCall(ISDestroy(setPointIS,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
         PetscCall(VecRestoreArrayF90(MEF90DefMechCtx%plasticStrain,plasticStrainArray,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))
      PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr))
      PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr))
      PetscCall(MatCopy(A,M,SAME_NONZERO_PATTERN,ierr))
   End Subroutine MEF90DefMechBilinearFormDisplacement

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechWork"
!!!
!!!  
!!!  MEF90DefMechWork:
!!!  
!!!  (c) 2014-2022 Blaise Bourdin bourdin@lsu.edu
!!!         2022 Blaise Bourdin bourdin@mcmaster.ca
!!!      2022 Alexis Marboeuf, marboeua@mcmaster.ca
!!!

   Subroutine MEF90DefMechWork(MEF90DefMechCtx,bodyForceWork,boundaryForceWork,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,Dimension(:),Pointer                     :: bodyForceWork,boundaryForceWork
      PetscErrorCode,Intent(INOUT)                       :: ierr

      Type(tDM)                                          :: dmDisplacement,dmBodyForce,dmBoundaryForce,dmPressureForce
      Type(tPetscSection)                                :: sectionBodyForce,sectionBoundaryForce,sectionPressureForce
      PetscReal,Dimension(:),Pointer                     :: displacementDof,bodyForceArray,boundaryForceArray,pressureForceArray
      Type(tIS)                                          :: setIS,setPointIS
      PetscInt,dimension(:),Pointer                      :: setID,setPointID
      PetscInt                                           :: set,QuadratureOrder,cell,dim,vecOffset
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90DefMechFaceSetOptions_Type),Pointer      :: faceSetOptions
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemVect
      DMPolytopeType                                     :: cellGeometry
      Type(MEF90ElementType)                             :: elemVectType
      
      Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90CtxGlobalOptions
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
      Type(MEF90_VECT)                                   :: bodyForce,boundaryForce,pressureForce
      PetscReal                                          :: myWork
      PetscInt                                           :: iDof,iGauss,numDofDisplacement,numGauss
    
      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr))
      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%displacementLocal,dmDisplacement,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%bodyForce,dmBodyForce,ierr))
      PetscCall(DMGetLocalSection(dmBodyForce,sectionBodyForce,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%boundaryForce,dmBoundaryForce,ierr))
      PetscCall(DMGetLocalSection(dmBoundaryForce,sectionBoundaryForce,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%pressureForce,dmPressureForce,ierr))
      PetscCall(DMGetLocalSection(dmPressureForce,sectionPressureForce,ierr))

      PetscCall(DMGetDimension(dmDisplacement,dim,ierr))

      !!! get IS for cell sets
      PetscCall(DMGetLabelIdIS(dmDisplacement,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm,setIS,ierr))

      bodyForceWork = 0.0_Kr
      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         PetscCall(VecGetArrayF90(MEF90DefMechCtx%bodyForce,bodyForceArray,ierr))
         Do set = 1,size(setID)
            myWork = 0.0_Kr
            PetscCall(DMGetStratumIS(dmDisplacement,MEF90CellSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
               PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr))
         
               PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(DMPlexGetCellType(dmDisplacement,setPointID(1),cellGeometry,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellGeometry,elemVectType,ierr))

               !!! Allocate elements 
               QuadratureOrder = 2 * (elemVectType%order - 1)
               PetscCall(MEF90ElementCreate(dmDisplacement,setPointIS,elemVect,QuadratureOrder,elemVectType,ierr))

               numDofDisplacement = size(elemVect(1)%BF(:,1))
               numGauss = size(elemVect(1)%Gauss_C)

               If (norm2(cellSetOptions%bodyForce) /= 0.0_Kr) Then
                  Do cell = 1,size(setPointID)
                     PetscCall(PetscSectionGetOffset(sectionBodyForce,setPointID(cell),vecOffset,ierr))
                     bodyForce = bodyForceArray(vecOffset:vecOffset+SIZEOFMEF90_VECT)
                     PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDisplacement
                           myWork = myWork + elemVect(cell)%Gauss_C(iGauss) * displacementDof(iDof) * (bodyForce .DotP. elemVect(cell)%BF(iDof,iGauss))
                        End Do ! iDof numDofDisplacement
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                  End Do ! cell
               End If ! cellSetOptions%bodyForce

               PetscCall(MEF90ElementDestroy(elemVect,ierr))
            End If ! setPointIS
            PetscCall(ISDestroy(setPointIS,ierr))
            PetscCallMPI(MPI_AllReduce(myWork,bodyForceWork(set),1,MPIU_SCALAR,MPI_SUM,MEF90DefMechCtx%MEF90Ctx%comm,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
         PetscCall(VecRestoreArrayF90(MEF90DefMechCtx%bodyForce,bodyForceArray,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))

      !!! face-based contributions
      PetscCall(DMGetLabelIdIS(dmDisplacement,MEF90FaceSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm,setIS,ierr))

      boundaryForceWork = 0.0_Kr
      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         PetscCall(VecGetArrayF90(MEF90DefMechCtx%boundaryForce,boundaryForceArray,ierr))
         PetscCall(VecGetArrayF90(MEF90DefMechCtx%pressureForce,pressureForceArray,ierr))
         Do set = 1,size(setID)
            myWork = 0.0_Kr
            PetscCall(DMGetStratumIS(dmDisplacement,MEF90FaceSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
               PetscCall(PetscBagGetDataMEF90DefMechCtxFaceSetOptions(MEF90DefMechCtx%FaceSetOptionsBag(set),faceSetOptions,ierr))
               PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(DMPlexGetCellType(dmDisplacement,setPointID(1),cellGeometry,ierr))
               PetscCall(MEF90ElementGetTypeBoundary(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellGeometry,elemVectType,ierr))

               !!! Allocate elements 
               QuadratureOrder = 2*elemVectType%order
               PetscCall(MEF90ElementCreate(dmDisplacement,setPointIS,elemVect,QuadratureOrder,elemVectType,ierr))

               numDofDisplacement = size(elemVect(1)%BF(:,1))
               numGauss = size(elemVect(1)%Gauss_C)

               If (norm2(faceSetOptions%boundaryForce) /= 0.0_Kr) Then
                  Do cell = 1,size(setPointID)
                     PetscCall(PetscSectionGetOffset(sectionBoundaryForce,setPointID(cell),vecOffset,ierr))
                     boundaryForce = boundaryForceArray(vecOffset:vecOffset+SIZEOFMEF90_VECT)
                     PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDisplacement
                           myWork = myWork + elemVect(cell)%Gauss_C(iGauss) * displacementDof(iDof) * (boundaryForce .DotP. elemVect(cell)%BF(iDof,iGauss))
                        End Do ! iDof numDofDisplacement
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                  End Do ! cell
               End If ! faceSetOptions%boundaryForce

               If (faceSetOptions%pressureForce /= 0.0_Kr) Then
                  Do cell = 1,size(setPointID)
                     PetscCall(PetscSectionGetOffset(sectionPressureForce,setPointID(cell),vecOffset,ierr))
                     PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                     pressureForce = boundaryForceArray(vecOffset)*elemVect(cell)%outerNormal
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDisplacement
                           myWork = myWork + elemVect(cell)%Gauss_C(iGauss) * displacementDof(iDof) * (pressureForce .DotP. elemVect(cell)%BF(iDof,iGauss))
                        End Do ! iDof numDofDisplacement
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                  End Do ! cell
               End If ! faceSetOptions%pressureForce

               PetscCall(MEF90ElementDestroy(elemVect,ierr))
               PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
            End If ! setPointIS
            PetscCallMPI(MPI_AllReduce(myWork,boundaryForceWork(set),1,MPIU_SCALAR,MPI_SUM,MEF90DefMechCtx%MEF90Ctx%comm,ierr))
         End Do ! set
         PetscCall(VecRestoreArrayF90(MEF90DefMechCtx%boundaryForce,boundaryForceArray,ierr))
         PetscCall(VecRestoreArrayF90(MEF90DefMechCtx%pressureForce,pressureForceArray,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))

   End Subroutine MEF90DefMechWork

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCohesiveEnergy"
!!!
!!!  
!!!  MEF90DefMechCohesiveEnergy:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022 Blaise Bourdin bourdin@mcmaster.ca
!!!      2022 Alexis Marboeuf, marboeua@mcmaster.ca
!!!

   Subroutine MEF90DefMechCohesiveEnergy(MEF90DefMechCtx,cohesiveEnergy,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,Dimension(:),Pointer                     :: cohesiveEnergy
      PetscErrorCode,Intent(INOUT)                       :: ierr

      Type(tDM)                                          :: dmDisplacement,dmCohesiveDisplacement
      PetscReal,Dimension(:),Pointer                     :: displacementDof,cohesiveDisplacementDof
      Type(tIS)                                          :: setIS,setPointIS
      PetscInt,dimension(:),Pointer                      :: setID,setPointID
      PetscInt                                           :: set,QuadratureOrder,cell
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemVect
      DMPolytopeType                                     :: cellGeometry
      Type(MEF90ElementType)                             :: elemVectType
      
      Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90CtxGlobalOptions
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
      Type(MEF90_VECT)                                   :: U0Gauss
      PetscReal                                          :: myCohesiveEnergy
      PetscInt                                           :: iDof,iGauss,numDofDisplacement,numGauss
    
      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr))
      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%displacementLocal,dmDisplacement,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%cohesiveDisplacement,dmCohesiveDisplacement,ierr))

      !!! get IS for cell sets
      PetscCall(DMGetLabelIdIS(dmDisplacement,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm,setIS,ierr))

      cohesiveEnergy = 0.0_Kr
      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         Do set = 1,size(setID)
            myCohesiveEnergy= 0.0_Kr
            PetscCall(DMGetStratumIS(dmDisplacement,MEF90CellSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
               PetscCall(PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr))
               PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr))
         
               PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(DMPlexGetCellType(dmDisplacement,setPointID(1),cellGeometry,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellGeometry,elemVectType,ierr))

               !!! Allocate elements 
               QuadratureOrder = 2 * (elemVectType%order - 1)
               PetscCall(MEF90ElementCreate(dmDisplacement,setPointIS,elemVect,QuadratureOrder,elemVectType,ierr))

               numDofDisplacement = size(elemVect(1)%BF(:,1))
               numGauss = size(elemVect(1)%Gauss_C)

               If ((norm2(cellSetOptions%cohesiveDisplacement) /= 0.0_Kr) .AND. (matpropSet%cohesiveStiffness /= 0.0_Kr)) Then
                  Do cell = 1,size(setPointID)
                     U0Gauss = 0.0_Kr
                     PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                     PetscCall(DMPlexVecGetClosure(dmCohesiveDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%cohesiveDisplacement,setPointID(cell),cohesiveDisplacementDof,ierr))
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDisplacement
                           U0Gauss = U0Gauss + (displacementDof(iDof) - cohesiveDisplacementDof(iDof)) * elemVect(cell)%BF(iDof,iGauss)
                        End Do ! iDof numDofDisplacement
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmCohesiveDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%cohesiveDisplacement,setPointID(cell),cohesiveDisplacementDof,ierr))
                     U0Gauss = U0Gauss * matpropSet%cohesiveStiffness
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDisplacement
                           myCohesiveEnergy = myCohesiveEnergy + elemVect(cell)%Gauss_C(iGauss) * displacementDof(iDof) * (U0Gauss .DotP. elemVect(cell)%BF(iDof,iGauss))
                        End Do ! iDof numDofDisplacement
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                  End Do ! cell
               End If ! cohesiveDisplacement

               PetscCall(MEF90ElementDestroy(elemVect,ierr))
            End If ! setPointIS
            PetscCall(ISDestroy(setPointIS,ierr))
            PetscCallMPI(MPI_AllReduce(myCohesiveEnergy,cohesiveEnergy(set),1,MPIU_SCALAR,MPI_SUM,MEF90DefMechCtx%MEF90Ctx%comm,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))
   
   End Subroutine MEF90DefMechCohesiveEnergy


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechPlasticDissipation"
!!!
!!!  
!!!  MEF90DefMechPlasticDissipation: 
!!!  
!!!  (c) 2015 Erwan TANNE erwan.tanne@gmail.com
!!!

   Subroutine MEF90DefMechPlasticDissipation(MEF90DefMechCtx,plasticStrainOld,energy,ierr)
      Type(tVec),Intent(IN)                              :: plasticStrainOld
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,dimension(:),Pointer                     :: energy
      PetscErrorCode,Intent(INOUT)                       :: ierr
   
   End Subroutine MEF90DefMechPlasticDissipation



#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechElasticEnergy"
!!!
!!!  
!!!  MEF90DefMechElasticEnergy: 
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!         2022 Blaise Bourdin bourdin@mcmaster.ca
!!!         2022 Alexis Marboeuf, marboeua@mcmaster.ca
!!!

   Subroutine MEF90DefMechElasticEnergy(MEF90DefMechCtx,energy,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,dimension(:),Pointer                     :: energy
      PetscErrorCode,Intent(INOUT)                       :: ierr

      Type(tDM)                                          :: dmDisplacement,dmDamage,dmTemperature,dmPlasticStrain
      Type(tPetscSection)                                :: sectionPlasticStrain
      PetscReal,Dimension(:),Pointer                     :: displacementDof,damageDof,temperatureDof,plasticStrainArray
      PetscInt                                           :: vecOffset
      Type(MEF90_MATS)                                   :: plasticStrainCell
      Type(tIS)                                          :: setIS,setPointIS
      PetscInt,dimension(:),Pointer                      :: setID,setPointID
      PetscInt                                           :: set,QuadratureOrder,cell,dim
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemVect
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemScal
      DMPolytopeType                                     :: cellGeometry
      Type(MEF90ElementType)                             :: elemVectType,elemScalType
      
      Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90CtxGlobalOptions
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
      Class(MEF90DefMechAT_Type),Allocatable             :: ATModel
      Class(MEF90_DEFMECHSPLIT),Allocatable              :: Split
      PetscBool                                          :: cellIsElastic
      Type(MEF90_MATS)                                   :: inelasticStrainGauss
      PetscReal                                          :: damageGauss,temperatureGauss,myEnergy,EEDPlus,EEDMinus,elasticEnergyDensityGauss
      PetscInt                                           :: iDof,iGauss,numDofDisplacement,numDofDamage,numGauss
    
      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr))
      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%displacementLocal,dmDisplacement,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%damageLocal,dmDamage,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%TemperatureLocal,dmTemperature,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%plasticStrain,dmPlasticStrain,ierr))
      PetscCall(DMGetLocalSection(dmPlasticStrain,sectionPlasticStrain,ierr))

      PetscCall(DMGetDimension(dmDisplacement,dim,ierr))

      !!! get IS for cell sets
      PetscCall(DMGetLabelIdIS(dmDisplacement,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm,setIS,ierr))

      energy = 0.0_Kr
      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(VecGetArrayF90(MEF90DefMechCtx%plasticStrain,plasticStrainArray,ierr))
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         Do set = 1,size(setID)
            myEnergy = 0.0_Kr
            PetscCall(DMGetStratumIS(dmDisplacement,MEF90CellSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
               PetscCall(PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr))
               PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr))
         
               PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(DMPlexGetCellType(dmDisplacement,setPointID(1),cellGeometry,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellGeometry,elemVectType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellGeometry,elemScalType,ierr))

               !!! get the ATModel and split objects
               PetscCall(MEF90DefMechGetATModel(cellSetOptions,ATModel,cellIsElastic,ierr))
               If (cellSetOptions%unilateralContactHybrid) Then
                  Split = MEF90_DEFMECHSPLITNONE()
               Else
                  PetscCall(MEF90DefMechGetSplit(cellSetOptions,Split,ierr))
               End If
               !!! Allocate elements 
               QuadratureOrder = max(2*elemVectType%order,split%damageOrder + split%strainOrder)
               PetscCall(MEF90ElementCreate(dmDisplacement,setPointIS,elemVect,QuadratureOrder,elemVectType,ierr))
               PetscCall(MEF90ElementCreate(dmDamage,setPointIS,elemScal,QuadratureOrder,elemScalType,ierr))

               numDofDisplacement = size(elemVect(1)%BF(:,1))
               numDofDamage = size(elemScal(1)%BF(:,1))
               numGauss = size(elemVect(1)%Gauss_C)
               Do cell = 1,size(setPointID)
                  PetscCall(DMPlexVecGetClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                  PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                  PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,MEF90DefMechCtx%TemperatureLocal,setPointID(cell),temperatureDof,ierr))

                  PetscCall(PetscSectionGetOffset(sectionPlasticStrain,setPointID(cell),vecOffset,ierr))
                  plasticStrainCell = plasticStrainArray(vecOffset:vecOffset+SIZEOFMEF90_MATS)
                  Do iGauss = 1, numGauss
                     damageGauss = 0.0_Kr
                        If (.NOT. cellIsElastic) Then
                        Do iDof = 1,numDofDamage
                           damageGauss = damageGauss + damageDof(iDof) * elemScal(cell)%BF(iDof,iGauss)
                        End Do ! iDof numDofDamage
                     End If

                     inelasticStrainGauss = 0.0_Kr
                     Do iDof = 1,numDofDisplacement
                        inelasticStrainGauss = inelasticStrainGauss + displacementDof(iDof) * elemVect(cell)%GradS_BF(iDof,iGauss)
                     End Do ! iDof numDofDisplacement

                     temperatureGauss = 0.0_Kr
                     Do iDof = 1,numDofDamage
                        temperatureGauss = temperatureGauss + temperatureDof(iDof) * elemScal(cell)%BF(iDof,iGauss)
                     End Do ! iDof numDofDamage
                     inelasticStrainGauss = inelasticStrainGauss - (temperatureGauss * matpropSet%LinearThermalExpansion)

! #if MEF90_DIM == 2
!!! We need something along these lines
!!! Adding terms in planestrain for plasticity with tr(p) = 0
! If (.NOT. matProp%HookesLaw%isPlaneStress) Then
!    stressGauss = stressGauss +  stiffness * ( matProp%HookesLaw%lambda*trace(plasticStrainCell)*MEF90MatS2DIdentity )
! End If
! #endif

                     Call Split%EED(inelasticStrainGauss - plasticStrainCell,matpropSet%HookesLaw,EEDPlus,EEDMinus)
                     If (cellIsElastic) Then
                        elasticEnergyDensityGauss = EEDPlus + EEDMinus
                     Else
                        elasticEnergyDensityGauss = (ATmodel%a(damageGauss) + matpropSet%residualStiffness) * EEDPlus + EEDMinus
                     End If

                     myEnergy = myEnergy + elemVect(cell)%Gauss_C(iGauss) * elasticEnergyDensityGauss
                  End Do ! iGauss

                  PetscCall(DMPlexVecRestoreClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                  PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                  PetscCall(DMPlexVecRestoreClosure(dmTemperature,PETSC_NULL_SECTION,MEF90DefMechCtx%TemperatureLocal,setPointID(cell),temperatureDof,ierr))
               End Do ! cell
               PetscCall(MEF90ElementDestroy(elemVect,ierr))
               PetscCall(MEF90ElementDestroy(elemScal,ierr))
            End If ! setPointIS
            PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
            PetscCall(ISDestroy(setPointIS,ierr))
            PetscCallMPI(MPI_AllReduce(myEnergy,energy(set),1,MPIU_SCALAR,MPI_SUM,MEF90DefMechCtx%MEF90Ctx%comm,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
         PetscCall(VecRestoreArrayF90(MEF90DefMechCtx%plasticStrain,plasticStrainArray,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))   
   End Subroutine MEF90DefMechElasticEnergy

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechStress"
!!!
!!!  
!!!  MEF90DefMechStress: 
!!!  
!!!  (c) 2012-20 Blaise Bourdin bourdin@lsu.edu
!!!         2022 Blaise Bourdin bourdin@mcmaster.ca
!!!         2022 Alexis Marboeuf, marboeua@mcmaster.ca
!!!

   Subroutine MEF90DefMechStress(MEF90DefMechCtx,stress,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(INOUT)                       :: ierr
      Type(tVec),Intent(IN)                              :: stress

      Type(tDM)                                          :: dmDisplacement,dmDamage,dmTemperature,dmPlasticStrain,dmStress
      Type(tPetscSection)                                :: sectionPlasticStrain
      PetscReal,Dimension(:),Pointer                     :: displacementDof,damageDof,temperatureDof,plasticStrainArray
      Type(MEF90_MATS)                                   :: plasticStrainCell
      Type(tIS)                                          :: setIS,setPointIS
      PetscInt,dimension(:),Pointer                      :: setID,setPointID
      PetscInt                                           :: set,QuadratureOrder,cell,dim,vecOffset
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemVect
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemScal
      DMPolytopeType                                     :: cellGeometry,cellDamageType
      Type(MEF90ElementType)                             :: elemVectType,elemScalType
      PetscReal,Dimension(:),Pointer                     :: stressDof
      
      Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90CtxGlobalOptions
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
      Class(MEF90DefMechAT_Type),Allocatable             :: ATModel
      Class(MEF90_DEFMECHSPLIT),Allocatable              :: Split
      PetscBool                                          :: cellIsElastic
      Type(MEF90_MATS)                                   :: inelasticStrainGauss,stressGaussPlus,stressGaussMinus,stressGauss
      PetscReal                                          :: damageGauss,temperatureGauss
      PetscInt                                           :: iDof,iGauss,numDofDisplacement,numDofDamage,numGauss
    
      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr))
      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr))
      PetscCall(VecGetDM(stress,dmStress,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%displacementLocal,dmDisplacement,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%damageLocal,dmDamage,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%TemperatureLocal,dmTemperature,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%plasticStrain,dmPlasticStrain,ierr))
      PetscCall(DMGetLocalSection(dmPlasticStrain,sectionPlasticStrain,ierr))

      PetscCall(DMGetDimension(dmDisplacement,dim,ierr))

      PetscCall(VecSet(stress,0.0_Kr,ierr))

      !!! get IS for cell sets
      PetscCall(DMGetLabelIdIS(dmDisplacement,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm,setIS,ierr))

      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         PetscCall(VecGetArrayF90(MEF90DefMechCtx%plasticStrain,plasticStrainArray,ierr))
         Do set = 1,size(setID)
            PetscCall(DMGetStratumIS(dmDisplacement,MEF90CellSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
               PetscCall(PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr))
               PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr))
         
               PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(DMPlexGetCellType(dmDisplacement,setPointID(1),cellGeometry,ierr))
               PetscCall(DMPlexGetCellType(dmDamage,setPointID(1),cellDamageType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellGeometry,elemVectType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellDamageType,elemScalType,ierr))

               !!! get the ATModel and split objects
               PetscCall(MEF90DefMechGetATModel(cellSetOptions,ATModel,cellIsElastic,ierr))
               If (cellSetOptions%unilateralContactHybrid) Then
                  Split = MEF90_DEFMECHSPLITNONE()
               Else
                  PetscCall(MEF90DefMechGetSplit(cellSetOptions,Split,ierr))
               End If
               !!! Allocate elements 
               QuadratureOrder = max(2*elemVectType%order,split%damageOrder + split%strainOrder)
               PetscCall(MEF90ElementCreate(dmDisplacement,setPointIS,elemVect,QuadratureOrder,elemVectType,ierr))
               PetscCall(MEF90ElementCreate(dmDamage,setPointIS,elemScal,QuadratureOrder,elemScalType,ierr))

               numDofDisplacement = size(elemVect(1)%BF(:,1))
               numDofDamage = size(elemScal(1)%BF(:,1))
               numGauss = size(elemVect(1)%Gauss_C)

               Allocate(stressDof((dim*(dim+1_Ki))/2_Ki))

               Do cell = 1,size(setPointID)
                  damageGauss = 0.0_Kr
                  If (.NOT. cellIsElastic) Then
                     PetscCall(DMPlexVecGetClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDamage
                           damageGauss = damageGauss + damageDof(iDof) * elemScal(cell)%BF(iDof,iGauss)
                        End Do ! iDof numDofDamage
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                  End If

                  inelasticStrainGauss = 0.0_Kr
                  PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDisplacement
                        inelasticStrainGauss = inelasticStrainGauss + displacementDof(iDof) * elemVect(cell)%GradS_BF(iDof,iGauss)
                     End Do ! iDof numDofDisplacement
                  End Do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))

                  temperatureGauss = 0.0_Kr
                  PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,MEF90DefMechCtx%TemperatureLocal,setPointID(cell),temperatureDof,ierr))
                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDamage
                        temperatureGauss = temperatureGauss + temperatureDof(iDof) * elemScal(cell)%BF(iDof,iGauss)
                     End Do ! iDof numDofDamage
                  End Do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmTemperature,PETSC_NULL_SECTION,MEF90DefMechCtx%TemperatureLocal,setPointID(cell),temperatureDof,ierr))
                  inelasticStrainGauss = inelasticStrainGauss - (temperatureGauss * matpropSet%LinearThermalExpansion)

! #if MEF90_DIM == 2
!!! We need something along these lines
!!! Adding terms in planestrain for plasticity with tr(p) = 0
! If (.NOT. matProp%HookesLaw%isPlaneStress) Then
!    stressGauss = stressGauss +  stiffness * ( matProp%HookesLaw%lambda*trace(plasticStrainCell)*MEF90MatS2DIdentity )
! End If
! #endif

                  PetscCall(PetscSectionGetOffset(sectionPlasticStrain,setPointID(cell),vecOffset,ierr))
                  plasticStrainCell = plasticStrainArray(vecOffset:vecOffset+SIZEOFMEF90_MATS)

                  Call Split%DEED(inelasticStrainGauss - plasticStrainCell,matpropSet%HookesLaw,stressGaussPlus,stressGaussMinus)
                  If (cellIsElastic) Then
                     stressGauss = stressGaussPlus + stressGaussMinus
                  Else
                     stressGauss = (ATModel%a(damageGauss) + matpropSet%residualStiffness) * stressGaussPlus + stressGaussMinus
                  End If

                  ! Cast MATS into array
                  stressDof = stressGauss
                  PetscCall(DMPlexVecSetClosure(dmStress,PETSC_NULL_SECTION,stress,setPointID(cell),stressDof,INSERT_VALUES,ierr))
               End Do ! cell

               DeAllocate(stressDof)

               PetscCall(MEF90ElementDestroy(elemVect,ierr))
               PetscCall(MEF90ElementDestroy(elemScal,ierr))
            End If ! setPointIS
            PetscCall(ISDestroy(setPointIS,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
         PetscCall(VecRestoreArrayF90(MEF90DefMechCtx%plasticStrain,plasticStrainArray,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))
   End Subroutine MEF90DefMechStress

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamage"
!!!
!!!  
!!!  MEF90DefMechOperatorDamage: Build the operator. When called in SNES, the solution time should always match the target time, 
!!!                                    so there is no need for interpolation of the forcees, external, and boundary values
!!!  
!!!  (c) 2012-20 Blaise Bourdin bourdin@lsu.edu, Erwan Tanne erwan.tanne@gmail.com
!!!         2022 Blaise Bourdin bourdin@mcmaster.ca
!!!         2022 Alexis Marboeuf marboeua@mcmaster.ca
!!!

   Subroutine  MEF90DefMechOperatorDamage(snesDamage,damage,residual,MEF90DefMechCtx,ierr)
      Type(tSNES),Intent(IN)                             :: snesDamage
      Type(tVec),Intent(IN)                              :: damage
      Type(tVec),Intent(INOUT)                           :: residual
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(tDM)                                          :: dmDisplacement,dmDamage,dmTemperature,dmCohesiveDisplacement,dmPlasticStrain
      Type(tPetscSection)                                :: sectionPlasticStrain
      PetscReal,Dimension(:),Pointer                     :: displacementDof,damageDof,temperatureDof,plasticStrainArray
      Type(MEF90_MATS)                                   :: plasticStrainCell
      Type(tIS)                                          :: setIS,setPointIS
      PetscInt,dimension(:),Pointer                      :: setID,setPointID
      PetscInt                                           :: set,QuadratureOrder,cell,dim,vecOffset
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemVect
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemScal
      DMPolytopeType                                     :: cellGeometry,cellDamageType
      Type(MEF90ElementType)                             :: elemVectType,elemScalType
      PetscReal,Dimension(:),Pointer                     :: residualDof
      
      Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90CtxGlobalOptions
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
      Type(tVec)                                         :: locResidual
      Class(MEF90DefMechAT_Type),Allocatable             :: ATModel
      Class(MEF90_DEFMECHSPLIT),Allocatable              :: Split
      PetscBool                                          :: cellIsElastic
      Type(MEF90_MATS)                                   :: inelasticStrainGauss,C2
      Type(MEF90_VECT)                                   :: gradDamageGauss
      PetscReal                                          :: damageGauss,temperatureGauss,EEDGaussMinus,EEDGaussPlus,C1,C3
      PetscInt                                           :: iDof,iGauss,numDofDisplacement,numDofDamage,numGauss
    
      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr))
      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%displacementLocal,dmDisplacement,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%cohesiveDisplacement,dmCohesiveDisplacement,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%damageLocal,dmDamage,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%TemperatureLocal,dmTemperature,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%plasticStrain,dmPlasticStrain,ierr))
      PetscCall(DMGetLocalSection(dmPlasticStrain,sectionPlasticStrain,ierr))

      PetscCall(DMGetDimension(dmDisplacement,dim,ierr))

      PetscCall(DMGetLocalVector(dmDamage,locResidual,ierr))
      PetscCall(DMGlobalToLocal(dmDamage,damage,INSERT_VALUES,MEF90DefMechCtx%damageLocal,ierr))
      
      PetscCall(VecSet(residual,0.0_Kr,ierr))
      PetscCall(VecSet(locResidual,0.0_Kr,ierr))

      !!! get IS for cell sets
      PetscCall(DMGetLabelIdIS(dmDamage,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm,setIS,ierr))

      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         PetscCall(VecGetArrayF90(MEF90DefMechCtx%plasticStrain,plasticStrainArray,ierr))
         Do set = 1,size(setID)
            PetscCall(DMGetStratumIS(dmDamage,MEF90CellSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
               PetscCall(PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr))
               PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr))
         
               PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(DMPlexGetCellType(dmDisplacement,setPointID(1),cellGeometry,ierr))
               PetscCall(DMPlexGetCellType(dmDamage,setPointID(1),cellDamageType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellGeometry,elemVectType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellDamageType,elemScalType,ierr))

               !!! get the ATModel and split objects
               PetscCall(MEF90DefMechGetATModel(cellSetOptions,ATModel,cellIsElastic,ierr))
               If (cellSetOptions%unilateralContactHybrid) Then
                  Split = MEF90_DEFMECHSPLITNONE()
               Else
                  PetscCall(MEF90DefMechGetSplit(cellSetOptions,Split,ierr))
               End If
               !!! Allocate elements 
               QuadratureOrder = max(2*elemVectType%order,split%damageOrder + split%strainOrder)
               PetscCall(MEF90ElementCreate(dmDisplacement,setPointIS,elemVect,QuadratureOrder,elemVectType,ierr))
               PetscCall(MEF90ElementCreate(dmDamage,setPointIS,elemScal,QuadratureOrder,elemScalType,ierr))

               numDofDisplacement = size(elemVect(1)%BF(:,1))
               numDofDamage = size(elemScal(1)%BF(:,1))
               numGauss = size(elemScal(1)%Gauss_C)

               Allocate(residualDof(numDofDamage))

               Do cell = 1,size(setPointID)
                  residualDof = 0.0_Kr
                  damageGauss = 0.0_Kr
                  gradDamageGauss = 0.0_Kr
                  If (.NOT. cellIsElastic) Then
                     PetscCall(DMPlexVecGetClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDamage
                           damageGauss = damageGauss + damageDof(iDof) * elemScal(cell)%BF(iDof,iGauss)
                           gradDamageGauss = gradDamageGauss + damageDof(iDof) * elemScal(cell)%Grad_BF(iDof,iGauss)
                        End Do ! iDof numDofDamage
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                  End If

                  inelasticStrainGauss = 0.0_Kr
                  PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDisplacement
                        inelasticStrainGauss = inelasticStrainGauss + displacementDof(iDof) * elemVect(cell)%GradS_BF(iDof,iGauss)
                     End Do ! iDof numDofDisplacement
                  End Do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))

                  temperatureGauss = 0.0_Kr
                  PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,MEF90DefMechCtx%TemperatureLocal,setPointID(cell),temperatureDof,ierr))
                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDamage
                        temperatureGauss = temperatureGauss + temperatureDof(iDof) * elemScal(cell)%BF(iDof,iGauss)
                     End Do ! iDof numDofDamage
                  End Do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmTemperature,PETSC_NULL_SECTION,MEF90DefMechCtx%TemperatureLocal,setPointID(cell),temperatureDof,ierr))
                  inelasticStrainGauss = inelasticStrainGauss - (temperatureGauss * matpropSet%LinearThermalExpansion)

! #if MEF90_DIM == 2
!!! We need something along these lines
!!! Adding terms in planestrain for plasticity with tr(p) = 0
! If (.NOT. matProp%HookesLaw%isPlaneStress) Then
!    stressGauss = stressGauss +  stiffness * ( matProp%HookesLaw%lambda*trace(plasticStrainCell)*MEF90MatS2DIdentity )
! End If
! #endif

                  PetscCall(PetscSectionGetOffset(sectionPlasticStrain,setPointID(cell),vecOffset,ierr))
                  plasticStrainCell = plasticStrainArray(vecOffset:vecOffset+SIZEOFMEF90_MATS)

                  If (.NOT. cellIsElastic) Then
                     Call Split%EED(inelasticStrainGauss - plasticStrainCell,matpropSet%HookesLaw,EEDGaussPlus,EEDGaussMinus)
                  Else
                     EEDGaussPlus = 0.0_Kr
                  End If
                  C1 = matpropSet%fractureToughness / ATModel%cw * 0.25_Kr / matpropSet%internalLength
                  C2 = matpropSet%fractureToughness / ATModel%cw * 0.5_Kr * matpropSet%internalLength * matpropSet%toughnessAnisotropyMatrix
                  C3 = ATModel%Da(damageGauss) * EEDGaussPlus + C1 * ATModel%Dw(damageGauss)
                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDamage
                        residualDof(iDof) = residualDof(iDof) + elemScal(cell)%Gauss_C(iGauss) * ( &
                                             C3 * elemScal(cell)%BF(iDof,iGauss) + (C2 * gradDamageGauss .DotP. elemScal(cell)%Grad_BF(iDof,iGauss)))
                     End Do ! iDof numDofDamage
                  End Do ! iGauss
                  PetscCall(DMPlexVecSetClosure(dmDamage,PETSC_NULL_SECTION,locResidual,setPointID(cell),residualDof,ADD_VALUES,ierr))
               End Do ! cell

               ! !!! crack pressure
               ! If (cellSetOptions%crackPressure /= 0.0_Kr)) Then
               !    Do cell = 1,size(setPointID)
               !       residualDof = 0.0_Kr
               !       Do iGauss = 1,numGauss
               !          CrackPressureCell = 0.0_Kr
               !          PetscCall(DMPlexVecGetClosure(dmCrackPressure,PETSC_NULL_SECTION,MEF90DefMechCtx%crackPressure,setPointID(cell),crackPressureDof,ierr))
               !          Do iDof = 1,numDofDamage
               !             CrackPressureCell = CrackPressureCell + crackPressureDof(iDof) * elemScal(cell)%BF(iDof,iGauss) 
               !          End Do ! iDof numDofDamage
               !          PetscCall(DMPlexVecRestoreClosure(dmCrackPressure,PETSC_NULL_SECTION,MEF90DefMechCtx%crackPressure,setPointID(cell),crackPressureDof,ierr))
               !          displacementCell = 0.0_Kr
               !          PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
               !          Do iDof = 1,numDofDisplacement
               !             displacementCell = displacementCell + displacementDof(iDof) * elemVect(cell)%BF(iDof,iGauss) 
               !          End Do ! iDof numDofDisplacement
               !          PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
               !          Do iDof = 1,numDofDamage
               !             residualDof(iDof) = residualDof(iDof) + elemScal(cell)%Gauss_C(iGauss) * &
               !                                  CrackPressureCell * (displacementCell  .DotP. elemScal(cell)%Grad_BF(iDof,iGauss))
               !          End Do ! iDof numDofDamage
               !       End Do ! iGauss
               !       PetscCall(DMPlexVecSetClosure(dmDamage,PETSC_NULL_SECTION,locResidual,setPointID(cell),residualDof,ADD_VALUES,ierr))
               !    End Do ! cell
               ! End If ! crack Pressure

               DeAllocate(residualDof)

               PetscCall(MEF90ElementDestroy(elemVect,ierr))
               PetscCall(MEF90ElementDestroy(elemScal,ierr))
            End If ! setPointIS
            PetscCall(ISDestroy(setPointIS,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
         PetscCall(VecRestoreArrayF90(MEF90DefMechCtx%plasticStrain,plasticStrainArray,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))
      PetscCall(DMLocalToGlobalBegin(dmDamage,locResidual,ADD_VALUES,residual,ierr))
      PetscCall(DMLocalToGlobalEnd(dmDamage,locResidual,ADD_VALUES,residual,ierr))
      PetscCall(DMRestoreLocalVector(dmDamage,locResidual,ierr))
   End Subroutine MEF90DefMechOperatorDamage

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamage"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamage:
!!!  
!!!  (c) 2012-19 Blaise Bourdin bourdin@lsu.edu, Erwan Tanne erwan.tanne@gmail.com
!!!         2022 Blaise Bourdin bourdin@mcmaster.ca
!!!         2022 Alexis Marboeuf marboeua@mcmaster.ca
!!!

   Subroutine MEF90DefMechBilinearFormDamage(snesDamage,damage,A,M,flg,MEF90DefMechCtx,ierr)
      Type(tSNES),Intent(IN)                             :: snesDamage
      Type(tVec),Intent(IN)                              :: damage
      Type(tMat),Intent(INOUT)                           :: A,M
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(tDM)                                          :: dmDisplacement,dmDamage,dmTemperature,dmCohesiveDisplacement,dmPlasticStrain
      Type(tPetscSection)                                :: sectionPlasticStrain
      PetscReal,Dimension(:),Pointer                     :: displacementDof,damageDof,temperatureDof,plasticStrainArray
      Type(MEF90_MATS)                                   :: plasticStrainCell
      Type(tIS)                                          :: setIS,setPointIS
      PetscInt,dimension(:),Pointer                      :: setID,setPointID
      PetscInt                                           :: set,QuadratureOrder,cell,dim,vecOffset
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemVect
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemScal
      DMPolytopeType                                     :: cellGeometry,cellDamageType
      Type(MEF90ElementType)                             :: elemVectType,elemScalType
      PetscReal,Dimension(:),Pointer                     :: matDof
      
      Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90CtxGlobalOptions
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
      Type(tVec)                                         :: locResidual
      Class(MEF90DefMechAT_Type),Allocatable             :: ATModel
      Class(MEF90_DEFMECHSPLIT),Allocatable              :: Split
      PetscBool                                          :: cellIsElastic
      Type(MEF90_MATS)                                   :: inelasticStrainGauss,C2
      Type(MEF90_VECT)                                   :: gradDamageGauss
      PetscReal                                          :: damageGauss,temperatureGauss,EEDGaussMinus,EEDGaussPlus,C1,C3
      PetscInt                                           :: iDof,jDof,iGauss,numDofDisplacement,numDofDamage,numGauss
    
      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr))
      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%displacementLocal,dmDisplacement,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%cohesiveDisplacement,dmCohesiveDisplacement,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%damageLocal,dmDamage,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%TemperatureLocal,dmTemperature,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%plasticStrain,dmPlasticStrain,ierr))
      PetscCall(DMGetLocalSection(dmPlasticStrain,sectionPlasticStrain,ierr))

      PetscCall(DMGetDimension(dmDisplacement,dim,ierr))

      PetscCall(DMGetLocalVector(dmDamage,locResidual,ierr))
      PetscCall(DMGlobalToLocal(dmDamage,damage,INSERT_VALUES,MEF90DefMechCtx%damageLocal,ierr))

      PetscCall(MatZeroEntries(A,ierr))

      !!! get IS for cell sets
      PetscCall(DMGetLabelIdIS(dmDamage,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm,setIS,ierr))

      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         PetscCall(VecGetArrayF90(MEF90DefMechCtx%plasticStrain,plasticStrainArray,ierr))
         Do set = 1,size(setID)
            PetscCall(DMGetStratumIS(dmDamage,MEF90CellSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
               PetscCall(PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr))
               PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr))
         
               PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(DMPlexGetCellType(dmDisplacement,setPointID(1),cellGeometry,ierr))
               PetscCall(DMPlexGetCellType(dmDamage,setPointID(1),cellDamageType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellGeometry,elemVectType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellDamageType,elemScalType,ierr))

               !!! get the ATModel and split objects
               PetscCall(MEF90DefMechGetATModel(cellSetOptions,ATModel,cellIsElastic,ierr))
               If (cellSetOptions%unilateralContactHybrid) Then
                  Split = MEF90_DEFMECHSPLITNONE()
               Else
                  PetscCall(MEF90DefMechGetSplit(cellSetOptions,Split,ierr))
               End If
               !!! Allocate elements 
               QuadratureOrder = max(2*elemVectType%order,split%damageOrder + split%strainOrder)
               PetscCall(MEF90ElementCreate(dmDisplacement,setPointIS,elemVect,QuadratureOrder,elemVectType,ierr))
               PetscCall(MEF90ElementCreate(dmDamage,setPointIS,elemScal,QuadratureOrder,elemScalType,ierr))

               numDofDisplacement = size(elemVect(1)%BF(:,1))
               numDofDamage = size(elemScal(1)%BF(:,1))
               numGauss = size(elemScal(1)%Gauss_C)

               Allocate(matDof(numDofDamage*numDofDamage))
               
               Do cell = 1,size(setPointID)
                  matDof = 0.0_Kr
                  damageGauss = 0.0_Kr
                  gradDamageGauss = 0.0_Kr
                  If (.NOT. cellIsElastic) Then
                     PetscCall(DMPlexVecGetClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDamage
                           damageGauss = damageGauss + damageDof(iDof) * elemScal(cell)%BF(iDof,iGauss)
                           gradDamageGauss = gradDamageGauss + damageDof(iDof) * elemScal(cell)%Grad_BF(iDof,iGauss)
                        End Do ! iDof numDofDamage
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                  End If

                  inelasticStrainGauss = 0.0_Kr
                  PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDisplacement
                        inelasticStrainGauss = inelasticStrainGauss + displacementDof(iDof) * elemVect(cell)%GradS_BF(iDof,iGauss)
                     End Do ! iDof numDofDisplacement
                  End Do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))

                  temperatureGauss = 0.0_Kr
                  PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,MEF90DefMechCtx%TemperatureLocal,setPointID(cell),temperatureDof,ierr))
                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDamage
                        temperatureGauss = temperatureGauss + temperatureDof(iDof) * elemScal(cell)%BF(iDof,iGauss)
                     End Do ! iDof numDofDamage
                  End Do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmTemperature,PETSC_NULL_SECTION,MEF90DefMechCtx%TemperatureLocal,setPointID(cell),temperatureDof,ierr))
                  inelasticStrainGauss = inelasticStrainGauss - (temperatureGauss * matpropSet%LinearThermalExpansion)

! #if MEF90_DIM == 2
!!! We need something along these lines
!!! Adding terms in planestrain for plasticity with tr(p) = 0
! If (.NOT. matProp%HookesLaw%isPlaneStress) Then
!    stressGauss = stressGauss +  stiffness * ( matProp%HookesLaw%lambda*trace(plasticStrainCell)*MEF90MatS2DIdentity )
! End If
! #endif

                  PetscCall(PetscSectionGetOffset(sectionPlasticStrain,setPointID(cell),vecOffset,ierr))
                  plasticStrainCell = plasticStrainArray(vecOffset:vecOffset+SIZEOFMEF90_MATS)

                  If (.NOT. cellIsElastic) Then
                     Call Split%EED(inelasticStrainGauss - plasticStrainCell,matpropSet%HookesLaw,EEDGaussPlus,EEDGaussMinus)
                  Else
                     EEDGaussPlus = 0.0_Kr
                  End If
                  C1 = matpropSet%fractureToughness / ATModel%cw * 0.25_Kr / matpropSet%internalLength
                  C2 = matpropSet%fractureToughness / ATModel%cw * 0.5_Kr * matpropSet%internalLength * matpropSet%toughnessAnisotropyMatrix
                  C3 = ATModel%D2a(damageGauss) * EEDGaussPlus + C1 * ATModel%D2w(damageGauss)
                  Do iGauss = 1,numGauss
                     Do jDof = 0,numDofDamage-1
                        Do iDof = 1,numDofDamage
                           matDof(jDof*numDofDamage+iDof)  = matDof(jDof*numDofDamage+iDof)  + elemScal(cell)%Gauss_C(iGauss) * ( &
                                             C3 * elemScal(cell)%BF(iDof,iGauss) * elemScal(cell)%BF(jDof+1,iGauss) + (C2 * elemScal(cell)%Grad_BF(iDof,iGauss) .DotP. elemScal(cell)%Grad_BF(jDof+1,iGauss)))
                        End Do ! iDof numDofDamage
                     End Do ! jDof numDofDamage
                  End Do ! iGauss
                  PetscCall(DMPlexMatSetClosure(dmDamage,PETSC_NULL_SECTION,PETSC_NULL_SECTION,A,setPointID(cell),matDof,ADD_VALUES,ierr))
               End Do ! cell
               
               DeAllocate(matDof)

               PetscCall(MEF90ElementDestroy(elemVect,ierr))
               PetscCall(MEF90ElementDestroy(elemScal,ierr))
            End If ! setPointIS
            PetscCall(ISDestroy(setPointIS,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
         PetscCall(VecRestoreArrayF90(MEF90DefMechCtx%plasticStrain,plasticStrainArray,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))
      PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr))
      PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr))
      PetscCall(MatCopy(A,M,SAME_NONZERO_PATTERN,ierr))
   End Subroutine MEF90DefMechBilinearFormDamage

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSurfaceEnergy"
!!!
!!!  
!!!  MEF90DefMechSurfaceEnergy:
!!!  
!!!  (c) 2014-2020 Blaise Bourdin bourdin@lsu.edu
!!!           2022 Blaise Bourdin bourdin@mcmaster.ca
!!!           2022 Alexis Marboeuf marboeua@mcmaster.ca
!!!

   Subroutine MEF90DefMechSurfaceEnergy(MEF90DefMechCtx,energy,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,dimension(:),Pointer                     :: energy
      PetscErrorCode,Intent(INOUT)                       :: ierr

      Type(tDM)                                          :: dmDamage
      PetscReal,Dimension(:),Pointer                     :: damageDof
      Type(tIS)                                          :: setIS,setPointIS
      PetscInt,dimension(:),Pointer                      :: setID,setPointID
      PetscInt                                           :: set,QuadratureOrder,cell,dim
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemScal
      DMPolytopeType                                     :: cellDamageType
      Type(MEF90ElementType)                             :: elemScalType
      
      Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90CtxGlobalOptions
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
      Class(MEF90DefMechAT_Type),Allocatable             :: ATModel
      PetscBool                                          :: cellIsElastic
      Type(MEF90_MATS)                                   :: C2
      Type(MEF90_VECT)                                   :: gradDamageGauss
      PetscReal                                          :: damageGauss,C1,C3,myEnergy
      PetscInt                                           :: iDof,iGauss,numDofDamage,numGauss
    
      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr))
      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%damageLocal,dmDamage,ierr))

      PetscCall(DMGetDimension(dmDamage,dim,ierr))

      !!! get IS for cell sets
      PetscCall(DMGetLabelIdIS(dmDamage,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm,setIS,ierr))

      energy = 0.0_Kr
      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         Do set = 1,size(setID)
            myEnergy = 0.0_Kr
            PetscCall(DMGetStratumIS(dmDamage,MEF90CellSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
               PetscCall(PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr))
               PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr))
         
               PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(DMPlexGetCellType(dmDamage,setPointID(1),cellDamageType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellDamageType,elemScalType,ierr))

               !!! get the ATModel and split objects
               PetscCall(MEF90DefMechGetATModel(cellSetOptions,ATModel,cellIsElastic,ierr))
               !!! Allocate elements 
               QuadratureOrder = max(ATmodel%wOrder,2*(elemScalType%order-1))
               PetscCall(MEF90ElementCreate(dmDamage,setPointIS,elemScal,QuadratureOrder,elemScalType,ierr))

               numDofDamage = size(elemScal(1)%BF(:,1))
               numGauss = size(elemScal(1)%Gauss_C)

               Do cell = 1,size(setPointID)
                  damageGauss = 0.0_Kr
                  gradDamageGauss = 0.0_Kr
                  If (.NOT. cellIsElastic) Then
                     PetscCall(DMPlexVecGetClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDamage
                           damageGauss = damageGauss + damageDof(iDof) * elemScal(cell)%BF(iDof,iGauss)
                           gradDamageGauss = gradDamageGauss + damageDof(iDof) * elemScal(cell)%Grad_BF(iDof,iGauss)
                        End Do ! iDof numDofDamage
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                  End If

                  C1 = matpropSet%fractureToughness / ATModel%cw * 0.25_Kr / matpropSet%internalLength
                  C2 = matpropSet%fractureToughness / ATModel%cw * 0.25_Kr * matpropSet%internalLength * matpropSet%toughnessAnisotropyMatrix
                  C3 = C1 * ATModel%w(damageGauss)
                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDamage
                        myEnergy = myEnergy + elemScal(cell)%Gauss_C(iGauss) * ( &
                                          C3 * elemScal(cell)%BF(iDof,iGauss) + (C2 * gradDamageGauss .DotP. gradDamageGauss))
                     End Do ! iDof numDofDamage
                  End Do ! iGauss
               End Do ! cell

               PetscCall(MEF90ElementDestroy(elemScal,ierr))
            End If ! setPointIS
            PetscCall(ISDestroy(setPointIS,ierr))
            PetscCallMPI(MPI_AllReduce(myEnergy,energy(set),1,MPIU_SCALAR,MPI_SUM,MEF90DefMechCtx%MEF90Ctx%comm,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))
   End Subroutine MEF90DefMechSurfaceEnergy


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCrackVolume"
!!!
!!!  
!!!  MEF90DefMechCrackVolume:
!!!  
!!!  (c) 2016-2021 Erwan erwan.tanne@gmail.com, Blaise Bourdin bourdin@lsu.edu  
!!!         2022 Blaise Bourdin bourdin@mcmaster.ca
!!!         2022 Alexis Marboeuf marboeua@mcmaster.ca
!!!

   Subroutine MEF90DefMechCrackVolume(MEF90DefMechCtx,CrackVolume,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,Dimension(:),Pointer                     :: CrackVolume
      PetscErrorCode,Intent(INOUT)                       :: ierr

      Type(tDM)                                          :: dmDamage,dmDisplacement
      PetscReal,Dimension(:),Pointer                     :: damageDof,displacementDof
      Type(tIS)                                          :: setIS,setPointIS
      PetscInt,dimension(:),Pointer                      :: setID,setPointID
      PetscInt                                           :: set,QuadratureOrder,cell,dim
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemVect
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemScal
      DMPolytopeType                                     :: cellGeometry,cellDamageType
      Type(MEF90ElementType)                             :: elemVectType,elemScalType
      
      Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90CtxGlobalOptions
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
      Class(MEF90DefMechAT_Type),Allocatable             :: ATModel
      PetscBool                                          :: cellIsElastic
      Type(MEF90_VECT)                                   :: gradDamageGauss,displacementCell
      PetscReal                                          :: myCrackVolume
      PetscInt                                           :: iDof,iGauss,numDofDamage,numDofDisplacement,numGauss
    
      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr))
      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%damageLocal,dmDamage,ierr))
      PetscCall(VecGetDM(MEF90DefMechCtx%displacementLocal,dmDisplacement,ierr))

      PetscCall(DMGetDimension(dmDamage,dim,ierr))

      !!! get IS for cell sets
      PetscCall(DMGetLabelIdIS(dmDamage,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm,setIS,ierr))

      CrackVolume = 0.0_Kr
      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         Do set = 1,size(setID)
            myCrackVolume = 0.0_Kr
            PetscCall(DMGetStratumIS(dmDamage,MEF90CellSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
               PetscCall(PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr))
               PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr))
         
               PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(DMPlexGetCellType(dmDisplacement,setPointID(1),cellGeometry,ierr))
               PetscCall(DMPlexGetCellType(dmDamage,setPointID(1),cellDamageType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellGeometry,elemVectType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellDamageType,elemScalType,ierr))

               !!! get the ATModel and split objects
               PetscCall(MEF90DefMechGetATModel(cellSetOptions,ATModel,cellIsElastic,ierr))
               !!! Allocate elements 
               QuadratureOrder = 2 * (elemScalType%order - 1)
               PetscCall(MEF90ElementCreate(dmDisplacement,setPointIS,elemVect,QuadratureOrder,elemVectType,ierr))
               PetscCall(MEF90ElementCreate(dmDamage,setPointIS,elemScal,QuadratureOrder,elemScalType,ierr))

               numDofDisplacement = size(elemVect(1)%BF(:,1))
               numDofDamage = size(elemScal(1)%BF(:,1))
               numGauss = size(elemScal(1)%Gauss_C)

               Do cell = 1,size(setPointID)
                  gradDamageGauss = 0.0_Kr
                  If (.NOT. cellIsElastic) Then
                     PetscCall(DMPlexVecGetClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDamage
                           gradDamageGauss = gradDamageGauss + damageDof(iDof) * elemScal(cell)%Grad_BF(iDof,iGauss)
                        End Do ! iDof numDofDamage
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                  End If

                  displacementCell = 0.0_Kr
                  PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDisplacement
                        displacementCell = displacementCell + displacementDof(iDof) * elemVect(cell)%BF(iDof,iGauss)
                     End Do ! iDof numDofDisplacement
                  End Do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))

                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDamage
                        myCrackVolume = myCrackVolume + elemScal(cell)%Gauss_C(iGauss) * ( &
                                             gradDamageGauss .DotP. displacementCell )
                     End Do ! iDof numDofDamage
                  End Do ! iGauss
               End Do ! cell

               PetscCall(MEF90ElementDestroy(elemScal,ierr))
               PetscCall(MEF90ElementDestroy(elemVect,ierr))
            End If ! setPointIS
            PetscCall(ISDestroy(setPointIS,ierr))
            PetscCallMPI(MPI_AllReduce(myCrackVolume,CrackVolume(set),1,MPIU_SCALAR,MPI_SUM,MEF90DefMechCtx%MEF90Ctx%comm,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))
   End Subroutine MEF90DefMechCrackVolume
End Module MEF90_APPEND(m_MEF90_DefMechAssembly,MEF90_DIM)D
