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

   PetscReal,Parameter,Private                         :: cwLinSoft = 0.7853981634_Kr
Contains


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDisplacement"
!!!
!!!  
!!!  MEF90DefMechOperatorDisplacement: Build the operator. When called in SNES, the solution time should always match the target time, 
!!!                                    so there is no need for interpolation of the forcees, external, and boundary values
!!!  
!!!  (c) 2012-20 Blaise Bourdin bourdin@lsu.edu, Erwan Tanne erwan.tanne@gmail.com
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
      PetscReal,Dimension(:),Pointer                     :: cohesiveDisplacementDof,displacementDof,damageDof,temperatureDof,bodyForceDof,boundaryForceDof,pressureForceDof,plasticStrainDof
      Type(MEF90_MATS)                                   :: plasticStrainCell
      Type(tIS)                                          :: setIS,setPointIS
      PetscInt,dimension(:),Pointer                      :: setID,setPointID
      PetscInt,Dimension(:),Pointer                      :: vecOffset
      PetscInt                                           :: set,QuadratureOrder,cell,dim
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90DefMechFaceSetOptions_Type),Pointer      :: faceSetOptions
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemDamage
      DMPolytopeType                                     :: cellDisplacementType,cellDamageType
      Type(MEF90ElementType)                             :: elemDisplacementType,elemDamageType
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

      PetscCallA(DMGetDimension(dmDisplacement,dim,ierr))

      !PetscCall(DMGetLocalVector(dmDisplacement,locDisplacement,ierr))
      PetscCall(DMGetLocalVector(dmDisplacement,locResidual,ierr))
      !PetscCall(MEF90VecGlobalToLocalConstraint(displacement,MEF90DefMechCtx%displacementLocal,locDisplacement,ierr))
      PetscCallA(DMGlobalToLocal(dmDisplacement,displacement,INSERT_VALUES,MEF90DefMechCtx%displacementLocal,ierr))

      PetscCall(VecSet(residual,0.0_Kr,ierr))
      PetscCall(VecSet(locResidual,0.0_Kr,ierr))

      !!! get IS for cell sets
      PetscCall(DMGetLabelIdIS(dmDisplacement,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm,setIS,ierr))

      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         Do set = 1,size(setID)
            PetscCall(DMGetStratumIS(dmDisplacement,MEF90CellSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
               PetscCall(PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr))
               PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr))
         
               PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(DMPlexGetCellType(dmDisplacement,setPointID(1),cellDisplacementType,ierr))
               PetscCall(DMPlexGetCellType(dmDamage,setPointID(1),cellDamageType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellDisplacementType,elemDisplacementType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellDamageType,elemDamageType,ierr))

               !!! get the ATModel and split objects
               PetscCall(MEF90DefMechGetATModel(cellSetOptions,ATModel,cellIsElastic,ierr))
               If (cellSetOptions%unilateralContactHybrid) Then
                  Split = MEF90_DEFMECHSPLITNONE()
               Else
                  PetscCall(MEF90DefMechGetSplit(cellSetOptions,Split,ierr))
               End If
               !!! Allocate elements 
               QuadratureOrder = max(2*elemDisplacementType%order,split%damageOrder + split%strainOrder)
               PetscCall(MEF90ElementCreate(dmDisplacement,setPointIS,elemDisplacement,QuadratureOrder,elemDisplacementType,ierr))
               PetscCall(MEF90ElementCreate(dmDamage,setPointIS,elemDamage,QuadratureOrder,elemDamageType,ierr))

               numDofDisplacement = size(elemDisplacement(1)%BF(:,1))
               numDofDamage = size(elemDamage(1)%BF(:,1))
               numGauss = size(elemDisplacement(1)%Gauss_C)

               Allocate(residualDof(numDofDisplacement))
               Allocate(bodyForceDof(dim))
               Allocate(plasticStrainDof((dim*(dim+1_Ki))/2_Ki))
               Allocate(vecOffset(1_Ki))

               Do cell = 1,size(setPointID)
                  residualDof = 0.0_Kr

                  !!! Main term: [a(\alpha) sigma^+(u) + sigma^-(u)] . e(v)
                  damageGauss = 0.0_Kr
                  If (.NOT. cellIsElastic) Then
                     PetscCall(DMPlexVecGetClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDamage
                           damageGauss = damageGauss + damageDof(iDof) * elemDamage(cell)%BF(iDof,iGauss)
                        End Do ! iDof numDofDamage
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                  End If

                  inelasticStrainGauss = 0.0_Kr
                  PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDisplacement
                        inelasticStrainGauss = inelasticStrainGauss + displacementDof(iDof) * elemDisplacement(cell)%GradS_BF(iDof,iGauss)
                     End Do ! iDof numDofDisplacement
                  End Do ! iGauss

                  PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))

                  temperatureGauss = 0.0_Kr
                  PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,MEF90DefMechCtx%TemperatureLocal,setPointID(cell),temperatureDof,ierr))
                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDamage
                        temperatureGauss = temperatureGauss + temperatureDof(iDof) * elemDamage(cell)%BF(iDof,iGauss)
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

                  PetscCall(PetscSectionGetOffset(sectionPlasticStrain,setPointID(cell),vecOffset(1),ierr))
                  PetscCall(VecGetValues(MEF90DefMechCtx%plasticStrain,(dim*(dim+1_Ki))/2_Ki,vecOffset,plasticStrainDof,ierr))
                  plasticStrainCell = plasticStrainDof

                  Call Split%DEED(inelasticStrainGauss - plasticStrainCell,matpropSet%HookesLaw,stressGaussPlus,stressGaussMinus)
                  If (cellIsElastic) Then
                     stressGauss = stressGaussPlus + stressGaussMinus
                  Else
                     stressGauss = (ATModel%a(damageGauss) + matpropSet%residualStiffness) * stressGaussPlus + stressGaussMinus
                  End If

                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDisplacement
                        residualDof(iDof) = residualDof(iDof) + elemDisplacement(cell)%Gauss_C(iGauss) * (stressGauss .DotP. elemDisplacement(cell)%GradS_BF(iDof,iGauss))
                     End Do ! iDof numDofDisplacement
                  End Do ! iGauss

                  PetscCall(DMPlexVecSetClosure(dmDisplacement,PETSC_NULL_SECTION,locResidual,setPointID(cell),residualDof,ADD_VALUES,ierr))
               End Do ! cell

               If (norm2(cellSetOptions%bodyForce) /= 0.0_Kr) Then
                  Do cell = 1,size(setPointID)
                     residualDof = 0.0_Kr
                     !!! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !!! If this happens, we will need to protect this loop
                     PetscCall(PetscSectionGetOffset(sectionBodyForce,setPointID(cell),vecOffset(1),ierr))
                     PetscCall(VecGetValues(MEF90DefMechCtx%bodyForce,dim,vecOffset,bodyForceDof,ierr))
                     bodyForce = bodyForceDof
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDisplacement
                           residualDof(iDof) = residualDof(iDof) - elemDisplacement(cell)%Gauss_C(iGauss) * (bodyForce .DotP. elemDisplacement(cell)%BF(iDof,iGauss))
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
               !             displacementDampingGauss = displacementDampingGauss + ((displacementDof(iDof) - displacementPreviousStepDof(iDof)) * elemDisplacement(cell)%BF(iDof,iGauss))
               !          End Do ! iDof numDofDisplacement
               !          PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,locDisplacement,setPointID(cell),displacementDof,ierr))
               !          PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementPreviousStepLocal,setPointID(cell),displacementPreviousStepDof,ierr))
               !          displacementDampingGauss = displacementDampingGauss * MEF90DefMechGlobalOptions%dampingCoefficientDisplacement / MEF90DefMechCtx%timeStep
               !          Do iDof = 1,numDofDisplacement
               !             residualDof(iDof) = residualDof(iDof) + elemDisplacement(cell)%Gauss_C(iGauss) * &
               !                                  ( displacementDampingGauss .DotP. elemDisplacement(cell)%BF(iDof,iGauss) )
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
                           U0Gauss = U0Gauss + (displacementDof(iDof) - cohesiveDisplacementDof(iDof)) * elemDisplacement(cell)%BF(iDof,iGauss)
                        End Do ! iDof numDofDisplacement
                     End Do ! iGauss

                     PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                     PetscCall(DMPlexVecRestoreClosure(dmCohesiveDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%cohesiveDisplacement,setPointID(cell),cohesiveDisplacementDof,ierr))
                     U0Gauss = U0Gauss * matpropSet%cohesiveStiffness
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDisplacement
                           residualDof(iDof) = residualDof(iDof) + elemDisplacement(cell)%Gauss_C(iGauss) * (U0Gauss .DotP. elemDisplacement(cell)%BF(iDof,iGauss))
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
               !             GradientDamageGauss = GradientDamageGauss + damageDof(iDof) * elemDamage(cell)%Grad_BF(iDof,iGauss) 
               !          End Do ! iDof numDofDamage
               !          PetscCall(DMPlexVecRestoreClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
               !          CrackPressureCell = 0.0_Kr
               !          PetscCall(DMPlexVecGetClosure(dmCrackPressure,PETSC_NULL_SECTION,MEF90DefMechCtx%crackPressure,setPointID(cell),crackPressureDof,ierr))
               !          Do iDof = 1,numDofDamage
               !             CrackPressureCell = CrackPressureCell + crackPressureDof(iDof) * elemDamage(cell)%Grad_BF(iDof,iGauss) 
               !          End Do ! iDof numDofDamage
               !          PetscCall(DMPlexVecRestoreClosure(dmCrackPressure,PETSC_NULL_SECTION,MEF90DefMechCtx%crackPressure,setPointID(cell),crackPressureDof,ierr))
               !          Do iDof = 1,numDofDisplacement
               !             residualDof(iDof) = residualDof(iDof) + elemDisplacement(cell)%Gauss_C(iGauss) * &
               !                                  CrackPressureCell * (GradientDamageGauss  .DotP. elemDisplacement(cell)%BF(iDof,iGauss))
               !          End Do ! iDof numDofDisplacement
               !       End Do ! iGauss
               !       PetscCall(DMPlexVecSetClosure(dmDisplacement,PETSC_NULL_SECTION,locResidual,setPointID(cell),residualDof,ADD_VALUES,ierr))
               !    End Do ! cell
               ! End If ! crack Pressure

               DeAllocate(residualDof)
               DeAllocate(bodyForceDof)
               DeAllocate(plasticStrainDof)
               DeAllocate(vecOffset)

               PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(MEF90ElementDestroy(elemDisplacement,ierr))
               PetscCall(MEF90ElementDestroy(elemDamage,ierr))
            End If ! setPointIS
            PetscCall(ISDestroy(setPointIS,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))

      !!! face-based contributions
      PetscCall(DMGetLabelIdIS(dmDisplacement,MEF90FaceSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm,setIS,ierr))
      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         Do set = 1,size(setID)
            PetscCall(DMGetStratumIS(dmDisplacement,MEF90FaceSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
              PetscCall(PetscBagGetDataMEF90DefMechCtxFaceSetOptions(MEF90DefMechCtx%FaceSetOptionsBag(set),faceSetOptions,ierr))
              PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
              PetscCall(DMPlexGetCellType(dmDisplacement,setPointID(1),cellDisplacementType,ierr))
              PetscCall(DMPlexGetCellType(dmDamage,setPointID(1),cellDamageType,ierr))
              PetscCall(MEF90ElementGetTypeBoundary(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellDisplacementType,elemDisplacementType,ierr))
              PetscCall(MEF90ElementGetTypeBoundary(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellDamageType,elemDamageType,ierr))

              !!! Allocate elements 
              QuadratureOrder = 2*elemDisplacementType%order
              PetscCall(MEF90ElementCreate(dmDisplacement,setPointIS,elemDisplacement,QuadratureOrder,elemDisplacementType,ierr))
              PetscCall(MEF90ElementCreate(dmDamage,setPointIS,elemDamage,QuadratureOrder,elemDamageType,ierr))

              numDofDisplacement = size(elemDisplacement(1)%BF(:,1))
              numDofDamage = size(elemDamage(1)%BF(:,1))
              numGauss = size(elemDisplacement(1)%Gauss_C)

              Allocate(residualDof(numDofDisplacement))
              Allocate(boundaryForceDof(dim))
              Allocate(pressureForceDof(1_Ki))
              Allocate(vecOffset(1_Ki))

               If (norm2(faceSetOptions%boundaryForce) /= 0.0_Kr) Then
                  Do cell = 1,size(setPointID)
                     residualDof = 0.0_Kr
                     PetscCall(PetscSectionGetOffset(sectionBoundaryForce,setPointID(cell),vecOffset(1),ierr))
                     PetscCall(VecGetValues(MEF90DefMechCtx%boundaryForce,dim,vecOffset,boundaryForceDof,ierr))
                     boundaryForce = boundaryForceDof
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDisplacement
                           residualDof(iDof) = residualDof(iDof) - elemDisplacement(cell)%Gauss_C(iGauss) * (boundaryForce .DotP. elemDisplacement(cell)%BF(iDof,iGauss))
                        End Do ! iDof numDofDisplacement
                     End Do ! iGauss
                     PetscCall(DMPlexVecSetClosure(dmDisplacement,PETSC_NULL_SECTION,locResidual,setPointID(cell),residualDof,ADD_VALUES,ierr))
                  End Do ! cell
               End If ! faceSetOptions%boundaryForce

               If (faceSetOptions%pressureForce /= 0.0_Kr) Then
                  Do cell = 1,size(setPointID)
                     residualDof = 0.0_Kr
                     PetscCall(PetscSectionGetOffset(sectionPressureForce,setPointID(cell),vecOffset(1),ierr))
                     PetscCall(VecGetValues(MEF90DefMechCtx%pressureForce,1_Ki,vecOffset,pressureForceDof,ierr))
                     pressureForce = pressureForceDof(1)*elemDisplacement(cell)%outerNormal
                     Do iGauss = 1,numGauss  
                        Do iDof = 1,numDofDisplacement
                           residualDof(iDof) = residualDof(iDof) - elemDisplacement(cell)%Gauss_C(iGauss) * (pressureForce .DotP. elemDisplacement(cell)%BF(iDof,iGauss))
                        End Do ! iDof numDofDisplacement
                     End Do ! iGauss
                     PetscCall(DMPlexVecSetClosure(dmDisplacement,PETSC_NULL_SECTION,locResidual,setPointID(cell),residualDof,ADD_VALUES,ierr))
                  End Do ! cell
               End If ! faceSetOptions%boundaryForce

               DeAllocate(residualDof)
               DeAllocate(boundaryForceDof)
               DeAllocate(pressureForceDof)
               DeAllocate(vecOffset)
               PetscCall(MEF90ElementDestroy(elemDisplacement,ierr))
               PetscCall(MEF90ElementDestroy(elemDamage,ierr))
               PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
            End If ! setPointIS
            PetscCall(ISDestroy(setPointIS,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
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
      PetscInt                                           :: numDofDisplacement,numDofDamage,numGauss,set,cell,iGauss,iDof,jDof,QuadratureOrder,dim
      PetscReal                                          :: damageGauss,temperatureGauss
      PetscBool                                          :: cellIsElastic
      PetscInt,dimension(:),Pointer                      :: setID,setPointID,vecOffset
      PetscReal,Dimension(:),Pointer                     :: matDof,displacementDof,damageDof,temperatureDof,plasticStrainDof
      Type(MEF90_HOOKESLAW)                              :: AGaussPlus,AGaussMinus,AGauss
      Type(MEF90_MATS)                                   :: inelasticStrainGauss,plasticStrainCell,AGradS_BF
      Type(MEF90_VECT)                                   :: U0Gauss
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemDamage
      DMPolytopeType                                     :: cellDisplacementType,cellDamageType
      Type(MEF90ElementType)                             :: elemDisplacementType,elemDamageType
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

      PetscCallA(DMGetDimension(dmDisplacement,dim,ierr))

      !PetscCall(DMGetLocalVector(dmDisplacement,locDisplacement,ierr))
      !PetscCall(MEF90VecGlobalToLocalConstraint(displacement,MEF90DefMechCtx%displacementLocal,locDisplacement,ierr))
      PetscCallA(DMGlobalToLocal(dmDisplacement,displacement,INSERT_VALUES,MEF90DefMechCtx%displacementLocal,ierr))

      PetscCall(MatZeroEntries(A,ierr))

      !!! get IS for cell sets
      PetscCall(DMGetLabelIdIS(dmDisplacement,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm,setIS,ierr))

      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         Do set = 1,size(setID)
            PetscCall(DMGetStratumIS(dmDisplacement,MEF90CellSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
               PetscCall(PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr))
               PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr))
         
               PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(DMPlexGetCellType(dmDisplacement,setPointID(1),cellDisplacementType,ierr))
               PetscCall(DMPlexGetCellType(dmDamage,setPointID(1),cellDamageType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellDisplacementType,elemDisplacementType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellDamageType,elemDamageType,ierr))

               !!! get the ATModel and split objects
               PetscCall(MEF90DefMechGetATModel(cellSetOptions,ATModel,cellIsElastic,ierr))
               If (cellSetOptions%unilateralContactHybrid) Then
                  Split = MEF90_DEFMECHSPLITNONE()
               Else
                  PetscCall(MEF90DefMechGetSplit(cellSetOptions,Split,ierr))
               End If
               !!! Allocate elements 
               QuadratureOrder = max(2*elemDisplacementType%order,split%damageOrder + split%strainOrder)
               PetscCall(MEF90ElementCreate(dmDisplacement,setPointIS,elemDisplacement,QuadratureOrder,elemDisplacementType,ierr))
               PetscCall(MEF90ElementCreate(dmDamage,setPointIS,elemDamage,QuadratureOrder,elemDamageType,ierr))

               numDofDisplacement = size(elemDisplacement(1)%BF(:,1))
               numDofDamage = size(elemDamage(1)%BF(:,1))
               numGauss = size(elemDisplacement(1)%Gauss_C)

               Allocate(matDof(numDofDisplacement*numDofDisplacement))
               Allocate(plasticStrainDof((dim*(dim+1_Ki))/2_Ki))
               Allocate(vecOffset(1_Ki))

               Do cell = 1,size(setPointID)
                  matDof = 0.0_Kr
                  !!! Main term: [a(\alpha) sigma^+(u) + sigma^-(u)] . e(v)
                  damageGauss = 0.0_Kr
                  If (.NOT. cellIsElastic) Then
                     PetscCall(DMPlexVecGetClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDamage
                           damageGauss = damageGauss + damageDof(iDof) * elemDamage(cell)%BF(iDof,iGauss)
                        End Do ! iDof numDofDamage
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                  End If

                  inelasticStrainGauss = 0.0_Kr
                  PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDisplacement
                        inelasticStrainGauss = inelasticStrainGauss + displacementDof(iDof) * elemDisplacement(cell)%GradS_BF(iDof,iGauss)
                     End Do ! iDof numDofDisplacement
                  End Do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))

                  temperatureGauss = 0.0_Kr
                  PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,MEF90DefMechCtx%temperatureLocal,setPointID(cell),temperatureDof,ierr))
                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDamage
                        temperatureGauss = temperatureGauss + temperatureDof(iDof) * elemDamage(cell)%BF(iDof,iGauss)
                     End Do ! iDof numDofDamage
                  End Do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmTemperature,PETSC_NULL_SECTION,MEF90DefMechCtx%temperatureLocal,setPointID(cell),temperatureDof,ierr))
                  inelasticStrainGauss = inelasticStrainGauss - (temperatureGauss * matpropSet%LinearThermalExpansion)

! #if MEF90_DIM == 2
!!! We need something along these lines
!!! Adding terms in planestrain for plasticity with tr(p) = 0
! If (.NOT. matProp%HookesLaw%isPlaneStress) Then
!    stressGauss = stressGauss +  stiffness * ( matProp%HookesLaw%lambda*trace(plasticStrainCell)*MEF90MatS2DIdentity )
! End If
! #endif

                  PetscCall(PetscSectionGetOffset(sectionPlasticStrain,setPointID(cell),vecOffset(1),ierr))
                  PetscCall(VecGetValues(MEF90DefMechCtx%plasticStrain,(dim*(dim+1_Ki))/2_Ki,vecOffset,plasticStrainDof,ierr))
                  plasticStrainCell = plasticStrainDof

                  If (cellIsElastic) Then
                     AGauss = matpropSet%HookesLaw
                  Else
                     Call Split%D2EED(inelasticStrainGauss - plasticStrainCell,matpropSet%HookesLaw,AGaussPlus,AGaussMinus)
                     AGauss = (ATModel%a(damageGauss) + matpropSet%residualStiffness) * AGaussPlus + AGaussMinus
                  End If
                  Do iGauss = 1,numGauss
                     Do jDof = 0,numDofDisplacement-1
                        Do iDof = 1,numDofDisplacement
                           AGradS_BF = AGauss * elemDisplacement(cell)%GradS_BF(iDof,iGauss)
                           matDof(jDof*numDofDisplacement+iDof) = matDof(jDof*numDofDisplacement+iDof) + elemDisplacement(cell)%Gauss_C(iGauss) * (AGradS_BF .DotP. elemDisplacement(cell)%GradS_BF(jDof+1,iGauss))
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
               !                U0Gauss = elemDisplacement(cell)%BF(iDof,iGauss) * MEF90DefMechGlobalOptions%dampingCoefficientDisplacement / MEF90DefMechCtx%timeStep
               !                matDof(jDof*numDofDisplacement+iDof) = matDof(jDof*numDofDisplacement+iDof) + elemDisplacement(cell)%Gauss_C(iGauss) * (U0Gauss .DotP. elemDisplacement(cell)%BF(jDof+1,iGauss))
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
                              U0Gauss = matpropSet%cohesiveStiffness * elemDisplacement(cell)%BF(iDof,iGauss)
                              matDof(jDof*numDofDisplacement+iDof) = matDof(jDof*numDofDisplacement+iDof) + elemDisplacement(cell)%Gauss_C(iGauss) * (U0Gauss .DotP. elemDisplacement(cell)%BF(jDof+1,iGauss))
                           End Do ! iDof numDofDisplacement
                        End Do ! jDof numDofDisplacement
                     End Do ! iGauss
                     PetscCall(DMPlexMatSetClosure(dmDisplacement,PETSC_NULL_SECTION,PETSC_NULL_SECTION,A,setPointID(cell),matDof,ADD_VALUES,ierr))
                  End Do ! cell
               End If ! cohesiveDisplacement

               DeAllocate(matDof)
               DeAllocate(plasticStrainDof)
               DeAllocate(vecOffset)

               PetscCall(MEF90ElementDestroy(elemDisplacement,ierr))
               PetscCall(MEF90ElementDestroy(elemDamage,ierr))
               PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
            End If ! setPointIS
            PetscCall(ISDestroy(setPointIS,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))
      !PetscCall(DMRestoreLocalVector(dmDisplacement,locDisplacement,ierr))
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
!!!      2022 Alexis Marboeuf, marboeua@mcmaster.ca
!!!

   Subroutine MEF90DefMechWork(MEF90DefMechCtx,bodyForceWork,boundaryForceWork,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,Dimension(:),Pointer                     :: bodyForceWork,boundaryForceWork
      PetscErrorCode,Intent(INOUT)                       :: ierr

      Type(tDM)                                          :: dmDisplacement,dmBodyForce,dmBoundaryForce,dmPressureForce
      Type(tPetscSection)                                :: sectionBodyForce,sectionBoundaryForce,sectionPressureForce
      PetscReal,Dimension(:),Pointer                     :: displacementDof,bodyForceDof,boundaryForceDof,pressureForceDof
      Type(tIS)                                          :: setIS,setPointIS
      PetscInt,dimension(:),Pointer                      :: setID,setPointID,vecOffset
      PetscInt                                           :: set,QuadratureOrder,cell,dim
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90DefMechFaceSetOptions_Type),Pointer      :: faceSetOptions
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      DMPolytopeType                                     :: cellDisplacementType
      Type(MEF90ElementType)                             :: elemDisplacementType
      
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

      PetscCallA(DMGetDimension(dmDisplacement,dim,ierr))

      !!! get IS for cell sets
      PetscCall(DMGetLabelIdIS(dmDisplacement,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm,setIS,ierr))

      bodyForceWork = 0.0_Kr
      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         Do set = 1,size(setID)
            myWork = 0.0_Kr
            PetscCall(DMGetStratumIS(dmDisplacement,MEF90CellSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
               PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr))
         
               PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(DMPlexGetCellType(dmDisplacement,setPointID(1),cellDisplacementType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellDisplacementType,elemDisplacementType,ierr))

               !!! Allocate elements 
               QuadratureOrder = 2 * (elemDisplacementType%order - 1)
               PetscCall(MEF90ElementCreate(dmDisplacement,setPointIS,elemDisplacement,QuadratureOrder,elemDisplacementType,ierr))

               numDofDisplacement = size(elemDisplacement(1)%BF(:,1))
               numGauss = size(elemDisplacement(1)%Gauss_C)

               Allocate(bodyForceDof(dim))
               Allocate(vecOffset(1_Ki))

               If (norm2(cellSetOptions%bodyForce) /= 0.0_Kr) Then
                  Do cell = 1,size(setPointID)
                     PetscCall(PetscSectionGetOffset(sectionBodyForce,setPointID(cell),vecOffset(1),ierr))
                     PetscCall(VecGetValues(MEF90DefMechCtx%bodyForce,dim,vecOffset,bodyForceDof,ierr))
                     PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                     bodyForce = bodyForceDof
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDisplacement
                           myWork = myWork + elemDisplacement(cell)%Gauss_C(iGauss) * displacementDof(iDof) * (bodyForce .DotP. elemDisplacement(cell)%BF(iDof,iGauss))
                        End Do ! iDof numDofDisplacement
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                  End Do ! cell
               End If ! cellSetOptions%bodyForce

               DeAllocate(bodyForceDof)
               DeAllocate(vecOffset)

               PetscCall(MEF90ElementDestroy(elemDisplacement,ierr))
            End If ! setPointIS
            PetscCall(ISDestroy(setPointIS,ierr))
            PetscCallMPI(MPI_AllReduce(myWork,bodyForceWork(set),1,MPIU_SCALAR,MPI_SUM,MEF90DefMechCtx%MEF90Ctx%comm,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))

      !!! face-based contributions
      PetscCall(DMGetLabelIdIS(dmDisplacement,MEF90FaceSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm,setIS,ierr))

      boundaryForceWork = 0.0_Kr
      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         Do set = 1,size(setID)
            myWork = 0.0_Kr
            PetscCall(DMGetStratumIS(dmDisplacement,MEF90FaceSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
              PetscCall(PetscBagGetDataMEF90DefMechCtxFaceSetOptions(MEF90DefMechCtx%FaceSetOptionsBag(set),faceSetOptions,ierr))
              PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
              PetscCall(DMPlexGetCellType(dmDisplacement,setPointID(1),cellDisplacementType,ierr))
              PetscCall(MEF90ElementGetTypeBoundary(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellDisplacementType,elemDisplacementType,ierr))

              !!! Allocate elements 
              QuadratureOrder = 2*elemDisplacementType%order
              PetscCall(MEF90ElementCreate(dmDisplacement,setPointIS,elemDisplacement,QuadratureOrder,elemDisplacementType,ierr))

              numDofDisplacement = size(elemDisplacement(1)%BF(:,1))
              numGauss = size(elemDisplacement(1)%Gauss_C)

              Allocate(boundaryForceDof(dim))
              Allocate(pressureForceDof(1_Ki))
              Allocate(vecOffset(1_Ki))

               If (norm2(faceSetOptions%boundaryForce) /= 0.0_Kr) Then
                  Do cell = 1,size(setPointID)
                     PetscCall(PetscSectionGetOffset(sectionBoundaryForce,setPointID(cell),vecOffset(1),ierr))
                     PetscCall(VecGetValues(MEF90DefMechCtx%boundaryForce,dim,vecOffset,boundaryForceDof,ierr))
                     PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                     boundaryForce = boundaryForceDof
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDisplacement
                           myWork = myWork + elemDisplacement(cell)%Gauss_C(iGauss) * displacementDof(iDof) * (boundaryForce .DotP. elemDisplacement(cell)%BF(iDof,iGauss))
                        End Do ! iDof numDofDisplacement
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                  End Do ! cell
               End If ! faceSetOptions%boundaryForce

               If (faceSetOptions%pressureForce /= 0.0_Kr) Then
                  Do cell = 1,size(setPointID)
                     PetscCall(PetscSectionGetOffset(sectionPressureForce,setPointID(cell),vecOffset(1),ierr))
                     PetscCall(VecGetValues(MEF90DefMechCtx%pressureForce,1_Ki,vecOffset,pressureForceDof,ierr))
                     PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                     pressureForce = pressureForceDof(1)*elemDisplacement(cell)%outerNormal
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDisplacement
                           myWork = myWork + elemDisplacement(cell)%Gauss_C(iGauss) * displacementDof(iDof) * (pressureForce .DotP. elemDisplacement(cell)%BF(iDof,iGauss))
                        End Do ! iDof numDofDisplacement
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                  End Do ! cell
               End If ! faceSetOptions%pressureForce

               DeAllocate(boundaryForceDof)
               DeAllocate(pressureForceDof)
               DeAllocate(vecOffset)
               PetscCall(MEF90ElementDestroy(elemDisplacement,ierr))
               PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
            End If ! setPointIS
            PetscCallMPI(MPI_AllReduce(myWork,boundaryForceWork(set),1,MPIU_SCALAR,MPI_SUM,MEF90DefMechCtx%MEF90Ctx%comm,ierr))
         End Do ! set
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
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      DMPolytopeType                                     :: cellDisplacementType
      Type(MEF90ElementType)                             :: elemDisplacementType
      
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
               PetscCall(DMPlexGetCellType(dmDisplacement,setPointID(1),cellDisplacementType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellDisplacementType,elemDisplacementType,ierr))

               !!! Allocate elements 
               QuadratureOrder = 2 * (elemDisplacementType%order - 1)
               PetscCall(MEF90ElementCreate(dmDisplacement,setPointIS,elemDisplacement,QuadratureOrder,elemDisplacementType,ierr))

               numDofDisplacement = size(elemDisplacement(1)%BF(:,1))
               numGauss = size(elemDisplacement(1)%Gauss_C)

               If ((norm2(cellSetOptions%cohesiveDisplacement) /= 0.0_Kr) .AND. (matpropSet%cohesiveStiffness /= 0.0_Kr)) Then
                  Do cell = 1,size(setPointID)
                     U0Gauss = 0.0_Kr
                     PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                     PetscCall(DMPlexVecGetClosure(dmCohesiveDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%cohesiveDisplacement,setPointID(cell),cohesiveDisplacementDof,ierr))
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDisplacement
                           U0Gauss = U0Gauss + (displacementDof(iDof) - cohesiveDisplacementDof(iDof)) * elemDisplacement(cell)%BF(iDof,iGauss)
                        End Do ! iDof numDofDisplacement
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmCohesiveDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%cohesiveDisplacement,setPointID(cell),cohesiveDisplacementDof,ierr))
                     U0Gauss = U0Gauss * matpropSet%cohesiveStiffness
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDisplacement
                           myCohesiveEnergy = myCohesiveEnergy + elemDisplacement(cell)%Gauss_C(iGauss) * displacementDof(iDof) * (U0Gauss .DotP. elemDisplacement(cell)%BF(iDof,iGauss))
                        End Do ! iDof numDofDisplacement
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                  End Do ! cell
               End If ! cohesiveDisplacement

               PetscCall(MEF90ElementDestroy(elemDisplacement,ierr))
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
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemDamage
      DMPolytopeType                                     :: cellDisplacementType,cellDamageType
      Type(MEF90ElementType)                             :: elemDisplacementType,elemDamageType
      
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

      PetscCallA(DMGetDimension(dmDisplacement,dim,ierr))

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
               PetscCall(DMPlexGetCellType(dmDisplacement,setPointID(1),cellDisplacementType,ierr))
               PetscCall(DMPlexGetCellType(dmDamage,setPointID(1),cellDamageType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellDisplacementType,elemDisplacementType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellDamageType,elemDamageType,ierr))

               !!! get the ATModel and split objects
               PetscCall(MEF90DefMechGetATModel(cellSetOptions,ATModel,cellIsElastic,ierr))
               If (cellSetOptions%unilateralContactHybrid) Then
                  Split = MEF90_DEFMECHSPLITNONE()
               Else
                  PetscCall(MEF90DefMechGetSplit(cellSetOptions,Split,ierr))
               End If
               !!! Allocate elements 
               QuadratureOrder = max(2*elemDisplacementType%order,split%damageOrder + split%strainOrder)
               PetscCall(MEF90ElementCreate(dmDisplacement,setPointIS,elemDisplacement,QuadratureOrder,elemDisplacementType,ierr))
               PetscCall(MEF90ElementCreate(dmDamage,setPointIS,elemDamage,QuadratureOrder,elemDamageType,ierr))

               numDofDisplacement = size(elemDisplacement(1)%BF(:,1))
               numDofDamage = size(elemDamage(1)%BF(:,1))
               numGauss = size(elemDisplacement(1)%Gauss_C)
               Do cell = 1,size(setPointID)
                  damageGauss = 0.0_Kr
                  If (.NOT. cellIsElastic) Then
                     PetscCall(DMPlexVecGetClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDamage
                           damageGauss = damageGauss + damageDof(iDof) * elemDamage(cell)%BF(iDof,iGauss)
                        End Do ! iDof numDofDamage
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                  End If

                  inelasticStrainGauss = 0.0_Kr
                  PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDisplacement
                        inelasticStrainGauss = inelasticStrainGauss + displacementDof(iDof) * elemDisplacement(cell)%GradS_BF(iDof,iGauss)
                     End Do ! iDof numDofDisplacement
                  End Do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))

                  temperatureGauss = 0.0_Kr
                  PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,MEF90DefMechCtx%TemperatureLocal,setPointID(cell),temperatureDof,ierr))
                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDamage
                        temperatureGauss = temperatureGauss + temperatureDof(iDof) * elemDamage(cell)%BF(iDof,iGauss)
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

                  Call Split%EED(inelasticStrainGauss - plasticStrainCell,matpropSet%HookesLaw,EEDPlus,EEDMinus)
                  If (cellIsElastic) Then
                     elasticEnergyDensityGauss = EEDPlus + EEDMinus
                  Else
                     elasticEnergyDensityGauss = (ATmodel%a(damageGauss) + matpropSet%residualStiffness) * EEDPlus + EEDMinus
                  End If

                  Do iGauss = 1,numGauss
                     myEnergy = myEnergy + elemDisplacement(cell)%Gauss_C(iGauss) * elasticEnergyDensityGauss
                  End Do ! iGauss
               End Do ! cell
               PetscCall(MEF90ElementDestroy(elemDisplacement,ierr))
               PetscCall(MEF90ElementDestroy(elemDamage,ierr))
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
!!!         2022 Alexis Marboeuf, marboeua@mcmaster.ca
!!!

   Subroutine MEF90DefMechStress(MEF90DefMechCtx,stress,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(INOUT)                       :: ierr
      Type(tVec),Intent(IN)                              :: stress

      Type(tDM)                                          :: dmDisplacement,dmDamage,dmTemperature,dmPlasticStrain,dmStress
      Type(tPetscSection)                                :: sectionPlasticStrain
      PetscReal,Dimension(:),Pointer                     :: displacementDof,damageDof,temperatureDof,plasticStrainDof
      Type(MEF90_MATS)                                   :: plasticStrainCell
      Type(tIS)                                          :: setIS,setPointIS
      PetscInt,dimension(:),Pointer                      :: setID,setPointID
      PetscInt,Dimension(:),Pointer                      :: vecOffset
      PetscInt                                           :: set,QuadratureOrder,cell,dim
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemDamage
      DMPolytopeType                                     :: cellDisplacementType,cellDamageType
      Type(MEF90ElementType)                             :: elemDisplacementType,elemDamageType
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

      PetscCallA(DMGetDimension(dmDisplacement,dim,ierr))

      PetscCall(VecSet(stress,0.0_Kr,ierr))

      !!! get IS for cell sets
      PetscCall(DMGetLabelIdIS(dmDisplacement,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm,setIS,ierr))

      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         Do set = 1,size(setID)
            PetscCall(DMGetStratumIS(dmDisplacement,MEF90CellSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
               PetscCall(PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr))
               PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr))
         
               PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(DMPlexGetCellType(dmDisplacement,setPointID(1),cellDisplacementType,ierr))
               PetscCall(DMPlexGetCellType(dmDamage,setPointID(1),cellDamageType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellDisplacementType,elemDisplacementType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellDamageType,elemDamageType,ierr))

               !!! get the ATModel and split objects
               PetscCall(MEF90DefMechGetATModel(cellSetOptions,ATModel,cellIsElastic,ierr))
               If (cellSetOptions%unilateralContactHybrid) Then
                  Split = MEF90_DEFMECHSPLITNONE()
               Else
                  PetscCall(MEF90DefMechGetSplit(cellSetOptions,Split,ierr))
               End If
               !!! Allocate elements 
               QuadratureOrder = max(2*elemDisplacementType%order,split%damageOrder + split%strainOrder)
               PetscCall(MEF90ElementCreate(dmDisplacement,setPointIS,elemDisplacement,QuadratureOrder,elemDisplacementType,ierr))
               PetscCall(MEF90ElementCreate(dmDamage,setPointIS,elemDamage,QuadratureOrder,elemDamageType,ierr))

               numDofDisplacement = size(elemDisplacement(1)%BF(:,1))
               numDofDamage = size(elemDamage(1)%BF(:,1))
               numGauss = size(elemDisplacement(1)%Gauss_C)

               Allocate(plasticStrainDof((dim*(dim+1_Ki))/2_Ki))
               Allocate(stressDof((dim*(dim+1_Ki))/2_Ki))
               Allocate(vecOffset(1_Ki))

               Do cell = 1,size(setPointID)
                  damageGauss = 0.0_Kr
                  If (.NOT. cellIsElastic) Then
                     PetscCall(DMPlexVecGetClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDamage
                           damageGauss = damageGauss + damageDof(iDof) * elemDamage(cell)%BF(iDof,iGauss)
                        End Do ! iDof numDofDamage
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                  End If

                  inelasticStrainGauss = 0.0_Kr
                  PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDisplacement
                        inelasticStrainGauss = inelasticStrainGauss + displacementDof(iDof) * elemDisplacement(cell)%GradS_BF(iDof,iGauss)
                     End Do ! iDof numDofDisplacement
                  End Do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))

                  temperatureGauss = 0.0_Kr
                  PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,MEF90DefMechCtx%TemperatureLocal,setPointID(cell),temperatureDof,ierr))
                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDamage
                        temperatureGauss = temperatureGauss + temperatureDof(iDof) * elemDamage(cell)%BF(iDof,iGauss)
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

                  PetscCall(PetscSectionGetOffset(sectionPlasticStrain,setPointID(cell),vecOffset(1),ierr))
                  PetscCall(VecGetValues(MEF90DefMechCtx%plasticStrain,(dim*(dim+1_Ki))/2_Ki,vecOffset,plasticStrainDof,ierr))
                  plasticStrainCell = plasticStrainDof

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

               DeAllocate(plasticStrainDof)
               DeAllocate(vecOffset)
               DeAllocate(stressDof)

               PetscCall(MEF90ElementDestroy(elemDisplacement,ierr))
               PetscCall(MEF90ElementDestroy(elemDamage,ierr))
            End If ! setPointIS
            PetscCall(ISDestroy(setPointIS,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))
         
   End Subroutine MEF90DefMechStress


! #undef __FUNCT__
! #define __FUNCT__ "MEF90DefMechBilinearFormDamageLinSoftLoc"
! !!!
! !!!  
! !!!  MEF90DefMechBilinearFormDamageLinSoftLoc: Assembles the bilinear form for LinSoft, that is:
! !!!    <K(alpha) beta , delta> = \int a''(\alpha) W(e(u)) beta \delta 
! !!!                            + Gc/pi/\ell w''(alpha) beta \delta
! !!!                            + 2 Gc/pi*\ell \nabla beta \nabla \delta
! !!! with a(\alpha) = (1-w(\alpha)) / ( 1 + (k-1) w(\alpha)) and w(\alpha) = 1 - (1-\alpha)^2
! !!! i.e. a''(\alpha) = 2 k * ( k + 3(k-1) v^2) / [k + (1-k) v^2]^3, w''(\alpha) = -2
! !!!  
! !!!  (c) 2014-2015 Erwan Tanne erwan.tanne@gmail.com, Blaise Bourdin bourdin@lsu.edu
! !!!

!    Subroutine MEF90DefMechBilinearFormDamageLinSoftLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
!       PetscReal,Dimension(:,:),Pointer                   :: Aloc
!       PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
!       Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
!       PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
!       Type(MEF90_MATPROP),Intent(IN)                     :: matprop
!       Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
!       Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

!       PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
!       PetscReal                                          :: elasticEnergyDensityGauss,temperatureGauss,damageGauss
!       Type(MEF90_MATS)                                   :: inelasticStrainGauss
!       PetscReal                                          :: C2,C1,N,k,C0
!       PetscLogDouble                                     :: flops
!       PetscErrorCode                                     :: ierr

!       numDofDisplacement = size(elemDisplacement%BF,1)
!       numDofDamage = size(elemDamage%BF,1)
!       numGauss = size(elemDamage%BF,2)
      
!       k = matprop%CoefficientLinSoft
!       C1 = 2.0_Kr * matprop%fractureToughness / matprop%internalLength / PETSC_PI
!       C2 = 2.0_Kr * matprop%fractureToughness * matprop%internalLength / PETSC_PI
!       N  = matprop%DuctileCouplingPower

!       Do iGauss = 1,numGauss
!          temperatureGauss = 0.0_Kr
!          If (Associated(temperatureDof)) Then
!             Do iDoF1 = 1,numDofDamage
!                temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
!             End Do
!          End If

!          inelasticStrainGauss = 0.0_Kr
!          Do iDoF1 = 1,numDofDisplacement
!             inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
!          End Do
!          inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion - plasticStrainCell
!          elasticEnergyDensityGauss = (matprop%HookesLaw * inelasticStrainGauss) .DotP. inelasticStrainGauss

!          damageGauss = 0.0_Kr
!          Do iDof1 = 1, numDofDamage
!             damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
!          End Do
!          C0 = k * (k + 3.0_Kr * (k - 1.0_Kr) * (1.0_Kr - damageGauss)**2) / &
!               (k + (1.0_Kr - k) * (1.0_Kr - damageGauss)**2)**3
!          !!! This is really twice the elastic energy density

!          Do iDoF1 = 1,numDofDamage
!             Do iDoF2 = 1,numDofDamage
!                Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) * ( & 
!                                   (elasticEnergyDensityGauss * C0 - C1 + &
!                                   (N * (N - 1.0) *( (1.0_Kr-damageGauss)**(N - 2.0) ) ) * cumulatedDissipatedPlasticEnergyCell ) * elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss) &
!                                   + C2 * ((matprop%toughnessAnisotropyMatrix * elemDamage%Grad_BF(iDoF1,iGauss)) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
!             End Do
!          End Do
!       End Do
!       !flops = 2 * numGauss * numDofDisplacement**2
!       !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
!    End Subroutine MEF90DefMechBilinearFormDamageLinSoftLoc

! #undef __FUNCT__
! #define __FUNCT__ "MEF90DefMechOperatorDamageLinSoftLoc"
! !!!
! !!!  
! !!!  MEF90DefMechOperatorDamageLinSoftLoc: Assembles the operator for LinSoft:
! !!!    F(\alpha)beta = \int [k(\alpha-1)W(e(u))] / [k-(k-1)(1-\alpha)^2]^2 beta 
! !!!                   + 2Gc/pi (1-\alpha) beta / \ell 
! !!!                   + 2Gc/pi \nabla \alpha \nabla beta * \ell
! !!!  
! !!!  (c) 2014-2015 Erwan Tanne erwan.tanne@gmail.com, Blaise Bourdin bourdin@lsu.edu
! !!!

!    Subroutine MEF90DefMechOperatorDamageLinSoftLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
!       PetscReal,Dimension(:),Pointer                     :: residualLoc
!       PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
!       Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
!       PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
!       Type(MEF90_MATPROP),Intent(IN)                     :: matprop
!       Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
!       Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
!       PetscReal,Intent(IN)                               :: CrackPressureCell

!       PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
!       PetscReal                                          :: elasticEnergyDensityGauss,temperatureGauss,damageGauss
!       Type(MEF90_MATS)                                   :: inelasticStrainGauss
!       Type(MEF90_VECT)                                   :: gradientDamageGauss
!       PetscReal                                          :: C1,C2,k,C0Operator,N
!       PetscLogDouble                                     :: flops
!       PetscErrorCode                                     :: ierr


!       numDofDisplacement = size(elemDisplacement%BF,1)
!       numDofDamage = size(elemDamage%BF,1)
!       numGauss = size(elemDamage%BF,2)
      
!       C1 = 2.0_Kr * matprop%fractureToughness / matprop%internalLength / PETSC_PI
!       C2 = 2.0_Kr * matprop%fractureToughness * matprop%internalLength / PETSC_PI
!       k = matprop%CoefficientLinSoft
!       N = matprop%DuctileCouplingPower

!       Do iGauss = 1,numGauss
!          temperatureGauss = 0.0_Kr
!          If (Associated(temperatureDof)) Then
!             Do iDoF1 = 1,numDofDamage
!                temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
!             End Do
!          End If

!          inelasticStrainGauss = 0.0_Kr
!          Do iDoF1 = 1,numDofDisplacement
!             inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
!          End Do
!          inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion - plasticStrainCell
!          elasticEnergyDensityGauss = (matprop%HookesLaw * inelasticStrainGauss) .DotP. inelasticStrainGauss
!          !!! This is really twice the elastic energy density

!          damageGauss = 0.0_Kr
!          gradientDamageGauss = 0.0_Kr
!          Do iDoF1 = 1,numDofDamage
!             damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
!             gradientDamageGauss = gradientDamageGauss + elemDamage%Grad_BF(iDoF1,iGauss) * xDof(iDoF1)
!          End Do
!          C0Operator = k * (damageGauss - 1.0_Kr) * elasticEnergyDensityGauss / &
!                      (k + (1.0_Kr - k) * (1.0_kr - damageGauss)**2)**2

!          Do iDoF2 = 1,numDofDamage
!             residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * ( &
!                                  (C0Operator + C1 * (1.0_Kr - damageGauss) &
!                                  - ( N * cumulatedDissipatedPlasticEnergyCell * (1.0_Kr - damageGauss)**(DBLE(N - 1.0))  ) ) * elemDamage%BF(iDoF2,iGauss) & 
!                                   + C2 * ((matprop%toughnessAnisotropyMatrix * gradientDamageGauss) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
!          End Do
!       End Do
!       !flops = 2 * numGauss * numDofDisplacement**2
!       !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
!    End Subroutine MEF90DefMechOperatorDamageLinSoftLoc



#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamage"
!!!
!!!  
!!!  MEF90DefMechOperatorDamage: Build the operator. When called in SNES, the solution time should always match the target time, 
!!!                                    so there is no need for interpolation of the forcees, external, and boundary values
!!!  
!!!  (c) 2012-20 Blaise Bourdin bourdin@lsu.edu, Erwan Tanne erwan.tanne@gmail.com
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
      PetscReal,Dimension(:),Pointer                     :: displacementDof,damageDof,temperatureDof,plasticStrainDof
      Type(MEF90_MATS)                                   :: plasticStrainCell
      Type(tIS)                                          :: setIS,setPointIS
      PetscInt,dimension(:),Pointer                      :: setID,setPointID,vecOffset
      PetscInt                                           :: set,QuadratureOrder,cell,dim
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemDamage
      DMPolytopeType                                     :: cellDisplacementType,cellDamageType
      Type(MEF90ElementType)                             :: elemDisplacementType,elemDamageType
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

      PetscCallA(DMGetDimension(dmDisplacement,dim,ierr))

      !PetscCall(DMGetLocalVector(dmDamage,locDamage,ierr))
      PetscCall(DMGetLocalVector(dmDamage,locResidual,ierr))
      !PetscCall(MEF90VecGlobalToLocalConstraint(damage,MEF90DefMechCtx%damageLocal,locDamage,ierr))
      PetscCall(DMGlobalToLocal(dmDamage,damage,INSERT_VALUES,MEF90DefMechCtx%damageLocal,ierr))
      
      PetscCall(VecSet(residual,0.0_Kr,ierr))
      PetscCall(VecSet(locResidual,0.0_Kr,ierr))

      !!! get IS for cell sets
      PetscCall(DMGetLabelIdIS(dmDamage,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm,setIS,ierr))

      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         Do set = 1,size(setID)
            PetscCall(DMGetStratumIS(dmDamage,MEF90CellSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
               PetscCall(PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr))
               PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr))
         
               PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(DMPlexGetCellType(dmDisplacement,setPointID(1),cellDisplacementType,ierr))
               PetscCall(DMPlexGetCellType(dmDamage,setPointID(1),cellDamageType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellDisplacementType,elemDisplacementType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellDamageType,elemDamageType,ierr))

               !!! get the ATModel and split objects
               PetscCall(MEF90DefMechGetATModel(cellSetOptions,ATModel,cellIsElastic,ierr))
               If (cellSetOptions%unilateralContactHybrid) Then
                  Split = MEF90_DEFMECHSPLITNONE()
               Else
                  PetscCall(MEF90DefMechGetSplit(cellSetOptions,Split,ierr))
               End If
               !!! Allocate elements 
               QuadratureOrder = max(2*elemDisplacementType%order,split%damageOrder + split%strainOrder)
               PetscCall(MEF90ElementCreate(dmDisplacement,setPointIS,elemDisplacement,QuadratureOrder,elemDisplacementType,ierr))
               PetscCall(MEF90ElementCreate(dmDamage,setPointIS,elemDamage,QuadratureOrder,elemDamageType,ierr))

               numDofDisplacement = size(elemDisplacement(1)%BF(:,1))
               numDofDamage = size(elemDamage(1)%BF(:,1))
               numGauss = size(elemDamage(1)%Gauss_C)

               Allocate(residualDof(numDofDamage))
               Allocate(plasticStrainDof((dim*(dim+1_Ki))/2_Ki))
               Allocate(vecOffset(1_Ki))

               Do cell = 1,size(setPointID)
                  residualDof = 0.0_Kr
                  damageGauss = 0.0_Kr
                  gradDamageGauss = 0.0_Kr
                  If (.NOT. cellIsElastic) Then
                     PetscCall(DMPlexVecGetClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDamage
                           damageGauss = damageGauss + damageDof(iDof) * elemDamage(cell)%BF(iDof,iGauss)
                           gradDamageGauss = gradDamageGauss + damageDof(iDof) * elemDamage(cell)%Grad_BF(iDof,iGauss)
                        End Do ! iDof numDofDamage
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                  End If

                  inelasticStrainGauss = 0.0_Kr
                  PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDisplacement
                        inelasticStrainGauss = inelasticStrainGauss + displacementDof(iDof) * elemDisplacement(cell)%GradS_BF(iDof,iGauss)
                     End Do ! iDof numDofDisplacement
                  End Do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))

                  temperatureGauss = 0.0_Kr
                  PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,MEF90DefMechCtx%TemperatureLocal,setPointID(cell),temperatureDof,ierr))
                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDamage
                        temperatureGauss = temperatureGauss + temperatureDof(iDof) * elemDamage(cell)%BF(iDof,iGauss)
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

                  PetscCall(PetscSectionGetOffset(sectionPlasticStrain,setPointID(cell),vecOffset(1),ierr))
                  PetscCall(VecGetValues(MEF90DefMechCtx%plasticStrain,(dim*(dim+1_Ki))/2_Ki,vecOffset,plasticStrainDof,ierr))
                  plasticStrainCell = plasticStrainDof

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
                        residualDof(iDof) = residualDof(iDof) + elemDamage(cell)%Gauss_C(iGauss) * ( &
                                             C3 * elemDamage(cell)%BF(iDof,iGauss) + (C2 * gradDamageGauss .DotP. elemDamage(cell)%Grad_BF(iDof,iGauss)))
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
               !             CrackPressureCell = CrackPressureCell + crackPressureDof(iDof) * elemDamage(cell)%BF(iDof,iGauss) 
               !          End Do ! iDof numDofDamage
               !          PetscCall(DMPlexVecRestoreClosure(dmCrackPressure,PETSC_NULL_SECTION,MEF90DefMechCtx%crackPressure,setPointID(cell),crackPressureDof,ierr))
               !          displacementCell = 0.0_Kr
               !          PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
               !          Do iDof = 1,numDofDisplacement
               !             displacementCell = displacementCell + displacementDof(iDof) * elemDisplacement(cell)%BF(iDof,iGauss) 
               !          End Do ! iDof numDofDisplacement
               !          PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
               !          Do iDof = 1,numDofDamage
               !             residualDof(iDof) = residualDof(iDof) + elemDamage(cell)%Gauss_C(iGauss) * &
               !                                  CrackPressureCell * (displacementCell  .DotP. elemDamage(cell)%Grad_BF(iDof,iGauss))
               !          End Do ! iDof numDofDamage
               !       End Do ! iGauss
               !       PetscCall(DMPlexVecSetClosure(dmDamage,PETSC_NULL_SECTION,locResidual,setPointID(cell),residualDof,ADD_VALUES,ierr))
               !    End Do ! cell
               ! End If ! crack Pressure

               DeAllocate(residualDof)
               DeAllocate(plasticStrainDof)
               DeAllocate(vecOffset)

               PetscCall(MEF90ElementDestroy(elemDisplacement,ierr))
               PetscCall(MEF90ElementDestroy(elemDamage,ierr))
            End If ! setPointIS
            PetscCall(ISDestroy(setPointIS,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
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
      PetscReal,Dimension(:),Pointer                     :: displacementDof,damageDof,temperatureDof,plasticStrainDof
      Type(MEF90_MATS)                                   :: plasticStrainCell
      Type(tIS)                                          :: setIS,setPointIS
      PetscInt,dimension(:),Pointer                      :: setID,setPointID,vecOffset
      PetscInt                                           :: set,QuadratureOrder,cell,dim
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemDamage
      DMPolytopeType                                     :: cellDisplacementType,cellDamageType
      Type(MEF90ElementType)                             :: elemDisplacementType,elemDamageType
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

      PetscCallA(DMGetDimension(dmDisplacement,dim,ierr))

      !PetscCall(DMGetLocalVector(dmDamage,locDamage,ierr))
      PetscCall(DMGetLocalVector(dmDamage,locResidual,ierr))
      !PetscCall(MEF90VecGlobalToLocalConstraint(damage,MEF90DefMechCtx%damageLocal,locDamage,ierr))
      PetscCall(DMGlobalToLocal(dmDamage,damage,INSERT_VALUES,MEF90DefMechCtx%damageLocal,ierr))

      PetscCall(MatZeroEntries(A,ierr))

      !!! get IS for cell sets
      PetscCall(DMGetLabelIdIS(dmDamage,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm,setIS,ierr))

      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         Do set = 1,size(setID)
            PetscCall(DMGetStratumIS(dmDamage,MEF90CellSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
               PetscCall(PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr))
               PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr))
         
               PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(DMPlexGetCellType(dmDisplacement,setPointID(1),cellDisplacementType,ierr))
               PetscCall(DMPlexGetCellType(dmDamage,setPointID(1),cellDamageType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellDisplacementType,elemDisplacementType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellDamageType,elemDamageType,ierr))

               !!! get the ATModel and split objects
               PetscCall(MEF90DefMechGetATModel(cellSetOptions,ATModel,cellIsElastic,ierr))
               If (cellSetOptions%unilateralContactHybrid) Then
                  Split = MEF90_DEFMECHSPLITNONE()
               Else
                  PetscCall(MEF90DefMechGetSplit(cellSetOptions,Split,ierr))
               End If
               !!! Allocate elements 
               QuadratureOrder = max(2*elemDisplacementType%order,split%damageOrder + split%strainOrder)
               PetscCall(MEF90ElementCreate(dmDisplacement,setPointIS,elemDisplacement,QuadratureOrder,elemDisplacementType,ierr))
               PetscCall(MEF90ElementCreate(dmDamage,setPointIS,elemDamage,QuadratureOrder,elemDamageType,ierr))

               numDofDisplacement = size(elemDisplacement(1)%BF(:,1))
               numDofDamage = size(elemDamage(1)%BF(:,1))
               numGauss = size(elemDamage(1)%Gauss_C)

               Allocate(matDof(numDofDamage*numDofDamage))
               Allocate(plasticStrainDof((dim*(dim+1_Ki))/2_Ki))
               Allocate(vecOffset(1_Ki))

               Do cell = 1,size(setPointID)
                  matDof = 0.0_Kr
                  damageGauss = 0.0_Kr
                  gradDamageGauss = 0.0_Kr
                  If (.NOT. cellIsElastic) Then
                     PetscCall(DMPlexVecGetClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDamage
                           damageGauss = damageGauss + damageDof(iDof) * elemDamage(cell)%BF(iDof,iGauss)
                           gradDamageGauss = gradDamageGauss + damageDof(iDof) * elemDamage(cell)%Grad_BF(iDof,iGauss)
                        End Do ! iDof numDofDamage
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                  End If

                  inelasticStrainGauss = 0.0_Kr
                  PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDisplacement
                        inelasticStrainGauss = inelasticStrainGauss + displacementDof(iDof) * elemDisplacement(cell)%GradS_BF(iDof,iGauss)
                     End Do ! iDof numDofDisplacement
                  End Do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))

                  temperatureGauss = 0.0_Kr
                  PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,MEF90DefMechCtx%TemperatureLocal,setPointID(cell),temperatureDof,ierr))
                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDamage
                        temperatureGauss = temperatureGauss + temperatureDof(iDof) * elemDamage(cell)%BF(iDof,iGauss)
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

                  PetscCall(PetscSectionGetOffset(sectionPlasticStrain,setPointID(cell),vecOffset(1),ierr))
                  PetscCall(VecGetValues(MEF90DefMechCtx%plasticStrain,(dim*(dim+1_Ki))/2_Ki,vecOffset,plasticStrainDof,ierr))
                  plasticStrainCell = plasticStrainDof

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
                           matDof(jDof*numDofDamage+iDof)  = matDof(jDof*numDofDamage+iDof)  + elemDamage(cell)%Gauss_C(iGauss) * ( &
                                             C3 * elemDamage(cell)%BF(iDof,iGauss) * elemDamage(cell)%BF(jDof+1,iGauss) + (C2 * elemDamage(cell)%Grad_BF(iDof,iGauss) .DotP. elemDamage(cell)%Grad_BF(jDof+1,iGauss)))
                        End Do ! iDof numDofDamage
                     End Do ! jDof numDofDamage
                  End Do ! iGauss
                  PetscCall(DMPlexMatSetClosure(dmDamage,PETSC_NULL_SECTION,PETSC_NULL_SECTION,A,setPointID(cell),matDof,ADD_VALUES,ierr))
               End Do ! cell
               
               DeAllocate(matDof)
               DeAllocate(plasticStrainDof)
               DeAllocate(vecOffset)

               PetscCall(MEF90ElementDestroy(elemDisplacement,ierr))
               PetscCall(MEF90ElementDestroy(elemDamage,ierr))
            End If ! setPointIS
            PetscCall(ISDestroy(setPointIS,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
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
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemDamage
      DMPolytopeType                                     :: cellDamageType
      Type(MEF90ElementType)                             :: elemDamageType
      
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

      PetscCallA(DMGetDimension(dmDamage,dim,ierr))

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
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellDamageType,elemDamageType,ierr))

               !!! get the ATModel and split objects
               PetscCall(MEF90DefMechGetATModel(cellSetOptions,ATModel,cellIsElastic,ierr))
               !!! Allocate elements 
               QuadratureOrder = max(ATmodel%wOrder,2*(elemDamageType%order-1))
               PetscCall(MEF90ElementCreate(dmDamage,setPointIS,elemDamage,QuadratureOrder,elemDamageType,ierr))

               numDofDamage = size(elemDamage(1)%BF(:,1))
               numGauss = size(elemDamage(1)%Gauss_C)

               Do cell = 1,size(setPointID)
                  damageGauss = 0.0_Kr
                  gradDamageGauss = 0.0_Kr
                  If (.NOT. cellIsElastic) Then
                     PetscCall(DMPlexVecGetClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDamage
                           damageGauss = damageGauss + damageDof(iDof) * elemDamage(cell)%BF(iDof,iGauss)
                           gradDamageGauss = gradDamageGauss + damageDof(iDof) * elemDamage(cell)%Grad_BF(iDof,iGauss)
                        End Do ! iDof numDofDamage
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                  End If

                  C1 = matpropSet%fractureToughness / ATModel%cw * 0.25_Kr / matpropSet%internalLength
                  C2 = matpropSet%fractureToughness / ATModel%cw * 0.25_Kr * matpropSet%internalLength * matpropSet%toughnessAnisotropyMatrix
                  C3 = C1 * ATModel%w(damageGauss)
                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDamage
                        myEnergy = myEnergy + elemDamage(cell)%Gauss_C(iGauss) * ( &
                                          C3 * elemDamage(cell)%BF(iDof,iGauss) + (C2 * gradDamageGauss .DotP. gradDamageGauss))
                     End Do ! iDof numDofDamage
                  End Do ! iGauss
               End Do ! cell

               PetscCall(MEF90ElementDestroy(elemDamage,ierr))
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
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemDamage
      DMPolytopeType                                     :: cellDisplacementType,cellDamageType
      Type(MEF90ElementType)                             :: elemDisplacementType,elemDamageType
      
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

      PetscCallA(DMGetDimension(dmDamage,dim,ierr))

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
               PetscCall(DMPlexGetCellType(dmDisplacement,setPointID(1),cellDisplacementType,ierr))
               PetscCall(DMPlexGetCellType(dmDamage,setPointID(1),cellDamageType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellDisplacementType,elemDisplacementType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellDamageType,elemDamageType,ierr))

               !!! get the ATModel and split objects
               PetscCall(MEF90DefMechGetATModel(cellSetOptions,ATModel,cellIsElastic,ierr))
               !!! Allocate elements 
               QuadratureOrder = 2 * (elemDamageType%order - 1)
               PetscCall(MEF90ElementCreate(dmDisplacement,setPointIS,elemDisplacement,QuadratureOrder,elemDisplacementType,ierr))
               PetscCall(MEF90ElementCreate(dmDamage,setPointIS,elemDamage,QuadratureOrder,elemDamageType,ierr))

               numDofDisplacement = size(elemDisplacement(1)%BF(:,1))
               numDofDamage = size(elemDamage(1)%BF(:,1))
               numGauss = size(elemDamage(1)%Gauss_C)

               Do cell = 1,size(setPointID)
                  gradDamageGauss = 0.0_Kr
                  If (.NOT. cellIsElastic) Then
                     PetscCall(DMPlexVecGetClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                     Do iGauss = 1,numGauss
                        Do iDof = 1,numDofDamage
                           gradDamageGauss = gradDamageGauss + damageDof(iDof) * elemDamage(cell)%Grad_BF(iDof,iGauss)
                        End Do ! iDof numDofDamage
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),damageDof,ierr))
                  End If

                  displacementCell = 0.0_Kr
                  PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))
                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDisplacement
                        displacementCell = displacementCell + displacementDof(iDof) * elemDisplacement(cell)%BF(iDof,iGauss)
                     End Do ! iDof numDofDisplacement
                  End Do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),displacementDof,ierr))

                  Do iGauss = 1,numGauss
                     Do iDof = 1,numDofDamage
                        myCrackVolume = myCrackVolume + elemDamage(cell)%Gauss_C(iGauss) * ( &
                                             gradDamageGauss .DotP. displacementCell )
                     End Do ! iDof numDofDamage
                  End Do ! iGauss
               End Do ! cell

               PetscCall(MEF90ElementDestroy(elemDamage,ierr))
               PetscCall(MEF90ElementDestroy(elemDisplacement,ierr))
            End If ! setPointIS
            PetscCall(ISDestroy(setPointIS,ierr))
            PetscCallMPI(MPI_AllReduce(myCrackVolume,CrackVolume(set),1,MPIU_SCALAR,MPI_SUM,MEF90DefMechCtx%MEF90Ctx%comm,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))
   
   End Subroutine MEF90DefMechCrackVolume
End Module MEF90_APPEND(m_MEF90_DefMechAssembly,MEF90_DIM)D
