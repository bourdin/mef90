#include "../MEF90/mef90.inc"
#include "../m_DefMech/mef90DefMech.inc"
module MEF90_APPEND(m_MEF90_DefMechAssembly,MEF90_DIM)D
#include "petsc/finclude/petsc.h"
#define MEF90_HDRegularization 0.01_Kr

use petscsnes
use petsctao
use m_MEF90_DefMechCtx
use m_MEF90_Materials
use MEF90_APPEND(m_MEF90_DefMechSplit,MEF90_DIM)D
use m_MEF90_DefMechAT

implicit none(type, external)
! Private
public MEF90DefMechOperatorDisplacement, &
   MEF90DefMechBilinearFormDisplacement, &
   MEF90DefMechWork, &
   MEF90DefMechCohesiveEnergy, &
   MEF90DefMechPlasticDissipation, &
   MEF90DefMechElasticEnergy, &
   MEF90DefMechOperatorDamage, &
   MEF90DefMechTAOGradientDamage, &
   MEF90DefMechBilinearFormDamage, &
   MEF90DefMechTAOHessianDamage, &
   MEF90DefMechSurfaceEnergy, &
   MEF90DefMechTAOObjectiveDamage, &
   MEF90DefMechCrackVolume, &
   MEF90DefMechStress

contains

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

subroutine MEF90DefMechOperatorDisplacement(snesDisplacement, displacement, residual, MEF90DefMechCtx, ierr)
   type(tSNES), intent(IN)                             :: snesDisplacement
   type(tVec), intent(IN)                              :: displacement
   type(tVec), intent(INOUT)                           :: residual
   type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
   PetscErrorCode, intent(INOUT)                       :: ierr

   type(tDM)                                           :: dmDisplacement, dmDamage, dmTemperature, dmCohesiveDisplacement, dmBodyForce, dmBoundaryForce, dmPressureForce, dmPlasticStrain
   type(tPetscSection)                                 :: sectionBodyForce, sectionBoundaryForce, sectionPressureForce, sectionPlasticStrain
   PetscReal, dimension(:), pointer                    :: cohesiveDisplacementDof, displacementDof, damageDof, temperatureDof
   PetscReal, dimension(:), pointer                    :: bodyForceArray, boundaryForceArray, pressureForceArray, plasticStrainArray
   type(MEF90_MATS)                                    :: plasticStrainCell
   type(tIS)                                           :: setIS, setPointIS
   PetscInt, dimension(:), pointer                     :: setID, setPointID
   PetscInt                                            :: set, QuadratureOrder, cell, dim, vecOffset
   type(MEF90_MATPROP), pointer                        :: matpropSet
   type(MEF90DefMechCellSetOptions_Type), pointer      :: cellSetOptions
   type(MEF90DefMechFaceSetOptions_Type), pointer      :: faceSetOptions
   type(MEF90_ELEMENT_ELAST), dimension(:), pointer    :: elemVect
   type(MEF90_ELEMENT_SCAL), dimension(:), pointer     :: elemScal
   type(eDMPolytopeType)                               :: cellGeometry
   type(MEF90ElementType)                              :: elemVectType, elemScalType
   PetscReal, dimension(:), pointer                    :: residualDof

   type(MEF90CtxGlobalOptions_Type), pointer           :: MEF90CtxGlobalOptions
   type(MEF90DefMechGlobalOptions_Type), pointer       :: MEF90DefMechGlobalOptions
   type(tVec)                                          :: locResidual
   class(MEF90DefMechAT_Type), allocatable             :: ATModel
   class(MEF90_DEFMECHSPLIT), allocatable              :: Split
   type(MEF90_MATS)                                    :: totalStrainGauss, stressGaussPlus, stressGaussMinus, stressGauss
   type(MEF90_VECT)                                    :: U0Gauss, bodyForce, boundaryForce, pressureForce
   PetscReal                                           :: damageGauss, temperatureGauss
   PetscInt                                            :: iDof, iGauss, numDofDisplacement, numDofDamage, numGauss
   character(len=MEF90MXSTRLEN)                        :: prefix

   PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag, MEF90CtxGlobalOptions, ierr))
   PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag, MEF90DefMechGlobalOptions, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%displacementLocal, dmDisplacement, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%cohesiveDisplacement, dmCohesiveDisplacement, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%damageLocal, dmDamage, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%TemperatureLocal, dmTemperature, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%bodyForce, dmBodyForce, ierr))
   PetscCall(DMGetLocalSection(dmBodyForce, sectionBodyForce, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%plasticStrain, dmPlasticStrain, ierr))
   PetscCall(DMGetLocalSection(dmPlasticStrain, sectionPlasticStrain, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%boundaryForce, dmBoundaryForce, ierr))
   PetscCall(DMGetLocalSection(dmBoundaryForce, sectionBoundaryForce, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%pressureForce, dmPressureForce, ierr))
   PetscCall(DMGetLocalSection(dmPressureForce, sectionPressureForce, ierr))

   PetscCall(DMGetDimension(dmDisplacement, dim, ierr))

   PetscCall(DMGetLocalVector(dmDisplacement, locResidual, ierr))
   PetscCall(DMGlobalToLocal(dmDisplacement, displacement, INSERT_VALUES, MEF90DefMechCtx%displacementLocal, ierr))

   PetscCall(VecSet(residual, 0.0_kr, ierr))
   PetscCall(VecSet(locResidual, 0.0_kr, ierr))

   !!! get IS for cell sets
   PetscCall(DMGetLabelIdIS(dmDisplacement, MEF90CellSetLabelName, setIS, ierr))
   PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm, setIS, ierr))

   if (.not. PetscObjectIsNull(setIS)) then
      PetscCall(VecGetArray(MEF90DefMechCtx%plasticStrain, plasticStrainArray, ierr))
      PetscCall(VecGetArray(MEF90DefMechCtx%bodyForce, bodyForceArray, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         PetscCall(DMGetStratumIS(dmDisplacement, MEF90CellSetLabelName, setID(set), setPointIS, ierr))
         if (.not. PetscObjectIsNull(setPointIS)) then
            PetscCall(PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set), matpropSet, ierr))
            PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set), cellSetOptions, ierr))

            PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
            PetscCall(DMPlexGetCellType(dmDisplacement, setPointID(1), cellGeometry, ierr))
            PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellGeometry, elemVectType, ierr))
            PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellGeometry, elemScalType, ierr))

            !!! get the ATModel and split objects
            write(prefix,'("cs",I4.4,"_")') setID(set)
            PetscCall(MEF90DefMechGetATModel(MEF90DefMechCtx%MEF90Ctx%comm, prefix, ATModel, ierr))
            PetscCall(ATModel%setFromOptions(ierr))

            if (cellSetOptions%unilateralContactHybrid) then
               Split = MEF90_DEFMECHSPLITNONE()
            else
               PetscCall(MEF90DefMechGetSplit(cellSetOptions, Split, ierr))
            end if

            !!! Allocate elements
            QuadratureOrder = max(2 * elemVectType%order, split%damageOrder + split%strainOrder)
            PetscCall(MEF90ElementCreate(dmDisplacement, setPointIS, elemVect, QuadratureOrder, elemVectType, ierr))
            PetscCall(MEF90ElementCreate(dmDamage, setPointIS, elemScal, QuadratureOrder, elemScalType, ierr))

            numDofDisplacement = size(elemVect(1)%BF(:, 1))
            numDofDamage = size(elemScal(1)%BF(:, 1))
            numGauss = size(elemVect(1)%Gauss_C)

            allocate (residualDof(numDofDisplacement))

            do cell = 1, size(setPointID)
               residualDof = 0.0_kr
               do iGauss = 1, numGauss
                  !!! Main term: [a(\alpha) sigma^+(u) + sigma^-(u)] . e(v)
                  damageGauss = 0.0_kr
                  if (.not. ATModel%isElastic) then
                     PetscCall(DMPlexVecGetClosure(dmDamage, PETSC_NULL_SECTION, MEF90DefMechCtx%damageLocal, setPointID(cell), PETSC_NULL_INTEGER, damageDof, ierr))
                     do iDof = 1, numDofDamage
                        damageGauss = damageGauss + damageDof(iDof) * elemScal(cell)%BF(iDof, iGauss)
                     end do ! iDof numDofDamage
                     PetscCall(DMPlexVecRestoreClosure(dmDamage, PETSC_NULL_SECTION, MEF90DefMechCtx%damageLocal, setPointID(cell), PETSC_NULL_INTEGER, damageDof, ierr))
                  end if

                  totalStrainGauss = 0.0_kr
                  PetscCall(DMPlexVecGetClosure(dmDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%displacementLocal, setPointID(cell), PETSC_NULL_INTEGER, displacementDof, ierr))
                  do iDof = 1, numDofDisplacement
                     totalStrainGauss = totalStrainGauss + displacementDof(iDof) * elemVect(cell)%GradS_BF(iDof, iGauss)
                  end do ! iDof numDofDisplacement

                  PetscCall(DMPlexVecRestoreClosure(dmDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%displacementLocal, setPointID(cell), PETSC_NULL_INTEGER, displacementDof, ierr))

                  temperatureGauss = 0.0_kr
                  PetscCall(DMPlexVecGetClosure(dmTemperature, PETSC_NULL_SECTION, MEF90DefMechCtx%TemperatureLocal, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
                  do iDof = 1, numDofDamage
                     temperatureGauss = temperatureGauss + temperatureDof(iDof) * elemScal(cell)%BF(iDof, iGauss)
                  end do ! iDof numDofDamage

                  PetscCall(DMPlexVecRestoreClosure(dmTemperature, PETSC_NULL_SECTION, MEF90DefMechCtx%TemperatureLocal, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
                  totalStrainGauss = totalStrainGauss - (temperatureGauss * matpropSet%LinearThermalExpansion)

                  ! #if MEF90_DIM == 2
   !!! We need something along these lines
   !!! Adding terms in planestrain for plasticity with tr(p) = 0
                  ! If (.NOT. matProp%HookesLaw%isPlaneStress) Then
                  !    stressGauss = stressGauss +  stiffness * ( matProp%HookesLaw%lambda*trace(plasticStrainCell)*MEF90MatS2DIdentity )
                  ! End If
                  ! #endif

                  PetscCall(PetscSectionGetOffset(sectionPlasticStrain, setPointID(cell), vecOffset, ierr))
                  plasticStrainCell = plasticStrainArray(vecOffset + 1:vecOffset + 1 + SIZEOFMEF90_MATS)

                  call Split%DEED(totalStrainGauss - plasticStrainCell, matpropSet%HookesLaw, stressGaussPlus, stressGaussMinus)
                  if (ATModel%isElastic) then
                     stressGauss = stressGaussPlus + stressGaussMinus
                  else
                     stressGauss = (ATModel%a(damageGauss) + matpropSet%residualStiffness) * stressGaussPlus + stressGaussMinus
                  end if

                  do iDof = 1, numDofDisplacement
                     residualDof(iDof) = residualDof(iDof) + elemVect(cell)%Gauss_C(iGauss) * (stressGauss .DotP. elemVect(cell)%GradS_BF(iDof, iGauss))
                  end do ! iDof numDofDisplacement
               end do ! iGauss
               PetscCall(DMPlexVecSetClosure(dmDisplacement, PETSC_NULL_SECTION, locResidual, setPointID(cell), residualDof, ADD_VALUES, ierr))
            end do ! cell

            if (norm2(cellSetOptions%bodyForce) /= 0.0_kr) then
               do cell = 1, size(setPointID)
                  residualDof = 0.0_kr
                     !!! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !!! If this happens, we will need to protect this loop
                  PetscCall(PetscSectionGetOffset(sectionBodyForce, setPointID(cell), vecOffset, ierr))
                  bodyForce = bodyForceArray(vecOffset + 1:vecOffset + 1 + SIZEOFMEF90_VECT)
                  do iGauss = 1, numGauss
                     do iDof = 1, numDofDisplacement
                        residualDof(iDof) = residualDof(iDof) - elemVect(cell)%Gauss_C(iGauss) * (bodyForce .DotP. elemVect(cell)%BF(iDof, iGauss))
                     end do ! iDof numDofDisplacement
                  end do ! iGauss
                  PetscCall(DMPlexVecSetClosure(dmDisplacement, PETSC_NULL_SECTION, locResidual, setPointID(cell), residualDof, ADD_VALUES, ierr))
               end do ! cell
            end if ! cellSetOptions%bodyForce

            ! !!! Damping
            ! If (MEF90DefMechGlobalOptions%dampingCoefficientDisplacement * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
            !    Do cell = 1,size(setPointID)
            !       residualDof = 0.0_Kr
            !       Do iGauss = 1,numGauss
            !          displacementDampingGauss = 0.0_Kr
            !          PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,locDisplacement,setPointID(cell),PETSC_NULL_INTEGER,displacementDof,ierr))
            !          PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementPreviousStepLocal,setPointID(cell),PETSC_NULL_INTEGER,displacementPreviousStepDof,ierr))
            !          Do iDof = 1,numDofDisplacement
            !             displacementDampingGauss = displacementDampingGauss + ((displacementDof(iDof) - displacementPreviousStepDof(iDof)) * elemVect(cell)%BF(iDof,iGauss))
            !          End Do ! iDof numDofDisplacement
            !          PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,locDisplacement,setPointID(cell),PETSC_NULL_INTEGER,displacementDof,ierr))
            !          PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementPreviousStepLocal,setPointID(cell),PETSC_NULL_INTEGER,displacementPreviousStepDof,ierr))
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
            if ((norm2(cellSetOptions%cohesiveDisplacement) /= 0.0_kr) .and. (matpropSet%cohesiveStiffness /= 0.0_kr)) then
               do cell = 1, size(setPointID)
                  residualDof = 0.0_kr
                  U0Gauss = 0.0_kr
                  PetscCall(DMPlexVecGetClosure(dmDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%displacementLocal, setPointID(cell), PETSC_NULL_INTEGER, displacementDof, ierr))
                  PetscCall(DMPlexVecGetClosure(dmCohesiveDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%cohesiveDisplacement, setPointID(cell), PETSC_NULL_INTEGER, cohesiveDisplacementDof, ierr))
                  do iGauss = 1, numGauss
                     do iDof = 1, numDofDisplacement
                        U0Gauss = U0Gauss + (displacementDof(iDof) - cohesiveDisplacementDof(iDof)) * elemVect(cell)%BF(iDof, iGauss)
                     end do ! iDof numDofDisplacement
                     U0Gauss = U0Gauss * matpropSet%cohesiveStiffness
                     do iDof = 1, numDofDisplacement
                        residualDof(iDof) = residualDof(iDof) + elemVect(cell)%Gauss_C(iGauss) * (U0Gauss .DotP. elemVect(cell)%BF(iDof, iGauss))
                     end do ! iDof numDofDisplacement
                  end do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%displacementLocal, setPointID(cell), PETSC_NULL_INTEGER, displacementDof, ierr))
                  PetscCall(DMPlexVecRestoreClosure(dmCohesiveDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%cohesiveDisplacement, setPointID(cell), PETSC_NULL_INTEGER, cohesiveDisplacementDof, ierr))
                  PetscCall(DMPlexVecSetClosure(dmDisplacement, PETSC_NULL_SECTION, locResidual, setPointID(cell), residualDof, ADD_VALUES, ierr))
               end do ! cell
            end if ! cohesiveDisplacement

            ! !!! crack pressure
            ! If (cellSetOptions%crackPressure /= 0.0_Kr)) Then
            !    Do cell = 1,size(setPointID)
            !       residualDof = 0.0_Kr
            !       Do iGauss = 1,numGauss
            !          gradientDamageGauss = 0.0_Kr
            !          PetscCall(DMPlexVecGetClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),PETSC_NULL_INTEGER,damageDof,ierr))
            !          Do iDof = 1,numDofDamage
            !             GradientDamageGauss = GradientDamageGauss + damageDof(iDof) * elemScal(cell)%Grad_BF(iDof,iGauss)
            !          End Do ! iDof numDofDamage
            !          PetscCall(DMPlexVecRestoreClosure(dmDamage,PETSC_NULL_SECTION,MEF90DefMechCtx%damageLocal,setPointID(cell),PETSC_NULL_INTEGER,damageDof,ierr))
            !          CrackPressureCell = 0.0_Kr
            !          PetscCall(DMPlexVecGetClosure(dmCrackPressure,PETSC_NULL_SECTION,MEF90DefMechCtx%crackPressure,setPointID(cell),PETSC_NULL_INTEGER,crackPressureDof,ierr))
            !          Do iDof = 1,numDofDamage
            !             CrackPressureCell = CrackPressureCell + crackPressureDof(iDof) * elemScal(cell)%Grad_BF(iDof,iGauss)
            !          End Do ! iDof numDofDamage
            !          PetscCall(DMPlexVecRestoreClosure(dmCrackPressure,PETSC_NULL_SECTION,MEF90DefMechCtx%crackPressure,setPointID(cell),PETSC_NULL_INTEGER,crackPressureDof,ierr))
            !          Do iDof = 1,numDofDisplacement
            !             residualDof(iDof) = residualDof(iDof) + elemVect(cell)%Gauss_C(iGauss) * &
            !                                  CrackPressureCell * (GradientDamageGauss  .DotP. elemVect(cell)%BF(iDof,iGauss))
            !          End Do ! iDof numDofDisplacement
            !       End Do ! iGauss
            !       PetscCall(DMPlexVecSetClosure(dmDisplacement,PETSC_NULL_SECTION,locResidual,setPointID(cell),residualDof,ADD_VALUES,ierr))
            !    End Do ! cell
            ! End If ! crack Pressure

            deallocate (residualDof)

            PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
            PetscCall(MEF90ElementDestroy(elemVect, ierr))
            PetscCall(MEF90ElementDestroy(elemScal, ierr))
            PetscCall(ISDestroy(setPointIS, ierr))
         end if ! setPointIS
      end do ! set
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(VecRestoreArray(MEF90DefMechCtx%plasticStrain, plasticStrainArray, ierr))
      PetscCall(VecRestoreArray(MEF90DefMechCtx%bodyForce, bodyForceArray, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end if ! setIS

      !!! face-based contributions
   PetscCall(DMGetLabelIdIS(dmDisplacement, MEF90FaceSetLabelName, setIS, ierr))
   PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm, setIS, ierr))
   if (.not. PetscObjectIsNull(setIS)) then
      PetscCall(ISGetIndices(setIS, setID, ierr))
      PetscCall(VecGetArray(MEF90DefMechCtx%boundaryForce, boundaryForceArray, ierr))
      PetscCall(VecGetArray(MEF90DefMechCtx%pressureForce, pressureForceArray, ierr))
      do set = 1, size(setID)
         PetscCall(DMGetStratumIS(dmDisplacement, MEF90FaceSetLabelName, setID(set), setPointIS, ierr))
         if (.not. PetscObjectIsNull(setPointIS)) then
            PetscCall(PetscBagGetDataMEF90DefMechCtxFaceSetOptions(MEF90DefMechCtx%FaceSetOptionsBag(set), faceSetOptions, ierr))
            PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
            PetscCall(DMPlexGetCellType(dmDisplacement, setPointID(1), cellGeometry, ierr))
            PetscCall(MEF90ElementGetTypeBoundary(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellGeometry, elemVectType, ierr))
            PetscCall(MEF90ElementGetTypeBoundary(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellGeometry, elemScalType, ierr))

              !!! Allocate elements
            QuadratureOrder = 2 * elemVectType%order
            PetscCall(MEF90ElementCreate(dmDisplacement, setPointIS, elemVect, QuadratureOrder, elemVectType, ierr))
            PetscCall(MEF90ElementCreate(dmDamage, setPointIS, elemScal, QuadratureOrder, elemScalType, ierr))

            numDofDisplacement = size(elemVect(1)%BF(:, 1))
            numDofDamage = size(elemScal(1)%BF(:, 1))
            numGauss = size(elemVect(1)%Gauss_C)

            allocate (residualDof(numDofDisplacement))

            if (norm2(faceSetOptions%boundaryForce) /= 0.0_kr) then
               do cell = 1, size(setPointID)
                  residualDof = 0.0_kr
                  PetscCall(PetscSectionGetOffset(sectionBoundaryForce, setPointID(cell), vecOffset, ierr))
                  boundaryForce = boundaryForceArray(vecOffset + 1:vecOffset + 1 + SIZEOFMEF90_VECT)
                  do iGauss = 1, numGauss
                     do iDof = 1, numDofDisplacement
                        residualDof(iDof) = residualDof(iDof) - elemVect(cell)%Gauss_C(iGauss) * (boundaryForce .DotP. elemVect(cell)%BF(iDof, iGauss))
                     end do ! iDof numDofDisplacement
                  end do ! iGauss
                  PetscCall(DMPlexVecSetClosure(dmDisplacement, PETSC_NULL_SECTION, locResidual, setPointID(cell), residualDof, ADD_VALUES, ierr))
               end do ! cell
            end if ! faceSetOptions%boundaryForce

            if (faceSetOptions%pressureForce /= 0.0_kr) then
               do cell = 1, size(setPointID)
                  residualDof = 0.0_kr
                  PetscCall(PetscSectionGetOffset(sectionPressureForce, setPointID(cell), vecOffset, ierr))
                  pressureForce = pressureForceArray(vecOffset + 1) * elemVect(cell)%outerNormal
                  do iGauss = 1, numGauss
                     do iDof = 1, numDofDisplacement
                        residualDof(iDof) = residualDof(iDof) - elemVect(cell)%Gauss_C(iGauss) * (pressureForce .DotP. elemVect(cell)%BF(iDof, iGauss))
                     end do ! iDof numDofDisplacement
                  end do ! iGauss
                  PetscCall(DMPlexVecSetClosure(dmDisplacement, PETSC_NULL_SECTION, locResidual, setPointID(cell), residualDof, ADD_VALUES, ierr))
               end do ! cell
            end if ! faceSetOptions%boundaryForce

            deallocate (residualDof)
            PetscCall(MEF90ElementDestroy(elemVect, ierr))
            PetscCall(MEF90ElementDestroy(elemScal, ierr))
            PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
            PetscCall(ISDestroy(setPointIS, ierr))
         end if ! setPointIS
      end do ! set
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(VecRestoreArray(MEF90DefMechCtx%boundaryForce, boundaryForceArray, ierr))
      PetscCall(VecRestoreArray(MEF90DefMechCtx%pressureForce, pressureForceArray, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end if ! setIS
   PetscCall(DMLocalToGlobalBegin(dmDisplacement, locResidual, ADD_VALUES, residual, ierr))
   PetscCall(DMLocalToGlobalEnd(dmDisplacement, locResidual, ADD_VALUES, residual, ierr))
   PetscCall(DMRestoreLocalVector(dmDisplacement, locResidual, ierr))
end subroutine MEF90DefMechOperatorDisplacement

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

subroutine MEF90DefMechBilinearFormDisplacement(snesDisplacement, displacement, A, M, MEF90DefMechCtx, ierr)
   type(tSNES), intent(IN)                             :: snesDisplacement
   type(tVec), intent(IN)                              :: displacement
   type(tMat), intent(INOUT)                           :: A, M
   type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
   PetscErrorCode, intent(INOUT)                       :: ierr

   type(tDM)                                           :: dmDisplacement, dmDamage, dmTemperature, dmPlasticStrain
   type(tPetscSection)                                 :: sectionPlasticStrain
   type(tIS)                                           :: setIS, setPointIS
   PetscInt                                            :: numDofDisplacement, numDofDamage, numGauss, set, cell, iGauss, iDof, jDof, QuadratureOrder, dim, vecOffset
   PetscReal                                           :: damageGauss, temperatureGauss
   PetscInt, dimension(:), pointer                     :: setID, setPointID
   PetscReal, dimension(:), pointer                    :: matDof, displacementDof, damageDof, temperatureDof, plasticStrainArray
   type(MEF90_HOOKESLAW)                               :: AGaussPlus, AGaussMinus, AGauss
   type(MEF90_MATS)                                    :: totalStrainGauss, plasticStrainCell, AGradS_BF
   type(MEF90_VECT)                                    :: U0Gauss
   type(MEF90_MATPROP), pointer                        :: matpropSet
   type(MEF90DefMechCellSetOptions_Type), pointer      :: cellSetOptions
   type(MEF90_ELEMENT_ELAST), dimension(:), pointer    :: elemVect
   type(MEF90_ELEMENT_SCAL), dimension(:), pointer     :: elemScal
   type(eDMPolytopeType)                               :: cellGeometry
   type(MEF90ElementType)                              :: elemVectType, elemScalType
   type(MEF90CtxGlobalOptions_Type), pointer           :: MEF90CtxGlobalOptions
   type(MEF90DefMechGlobalOptions_Type), pointer       :: MEF90DefMechGlobalOptions
   class(MEF90DefMechAT_Type), allocatable             :: ATModel
   class(MEF90_DEFMECHSPLIT), allocatable              :: Split
   character(len=MEF90MXSTRLEN)                        :: prefix

   PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag, MEF90CtxGlobalOptions, ierr))
   PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag, MEF90DefMechGlobalOptions, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%displacementLocal, dmDisplacement, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%damageLocal, dmDamage, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%temperatureLocal, dmTemperature, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%plasticStrain, dmPlasticStrain, ierr))
   PetscCall(DMGetLocalSection(dmPlasticStrain, sectionPlasticStrain, ierr))

   PetscCall(DMGetDimension(dmDisplacement, dim, ierr))

   PetscCall(DMGlobalToLocal(dmDisplacement, displacement, INSERT_VALUES, MEF90DefMechCtx%displacementLocal, ierr))

   PetscCall(MatZeroEntries(A, ierr))

   !!! get IS for cell sets
   PetscCall(DMGetLabelIdIS(dmDisplacement, MEF90CellSetLabelName, setIS, ierr))
   PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm, setIS, ierr))

   if (.not. PetscObjectIsNull(setIS)) then
      PetscCall(ISGetIndices(setIS, setID, ierr))
      PetscCall(VecGetArray(MEF90DefMechCtx%plasticStrain, plasticStrainArray, ierr))
      do set = 1, size(setID)
         PetscCall(DMGetStratumIS(dmDisplacement, MEF90CellSetLabelName, setID(set), setPointIS, ierr))
         if (.not. PetscObjectIsNull(setPointIS)) then
            PetscCall(PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set), matpropSet, ierr))
            PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set), cellSetOptions, ierr))

            PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
            PetscCall(DMPlexGetCellType(dmDisplacement, setPointID(1), cellGeometry, ierr))
            PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellGeometry, elemVectType, ierr))
            PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellGeometry, elemScalType, ierr))

            !!! get the ATModel and split objects
            write(prefix,'("cs",I4.4,"_")') setID(set)
            PetscCall(MEF90DefMechGetATModel(MEF90DefMechCtx%MEF90Ctx%comm, prefix, ATModel, ierr))
            PetscCall(ATModel%setFromOptions(ierr))
            if (cellSetOptions%unilateralContactHybrid) then
               Split = MEF90_DEFMECHSPLITNONE()
            else
               PetscCall(MEF90DefMechGetSplit(cellSetOptions, Split, ierr))
            end if

            !!! Allocate elements
            QuadratureOrder = max(2 * elemVectType%order, split%damageOrder + split%strainOrder)
            PetscCall(MEF90ElementCreate(dmDisplacement, setPointIS, elemVect, QuadratureOrder, elemVectType, ierr))
            PetscCall(MEF90ElementCreate(dmDamage, setPointIS, elemScal, QuadratureOrder, elemScalType, ierr))

            numDofDisplacement = size(elemVect(1)%BF(:, 1))
            numDofDamage = size(elemScal(1)%BF(:, 1))
            numGauss = size(elemVect(1)%Gauss_C)

            allocate (matDof(numDofDisplacement * numDofDisplacement))

            do cell = 1, size(setPointID)
               matDof = 0.0_kr
               do iGauss = 1, numGauss
                     !!! Main term: [a(\alpha) sigma^+(u) + sigma^-(u)] . e(v)
                  damageGauss = 0.0_kr
                  if (.not. ATModel%isElastic) then
                     PetscCall(DMPlexVecGetClosure(dmDamage, PETSC_NULL_SECTION, MEF90DefMechCtx%damageLocal, setPointID(cell), PETSC_NULL_INTEGER, damageDof, ierr))
                     do iDof = 1, numDofDamage
                        damageGauss = damageGauss + damageDof(iDof) * elemScal(cell)%BF(iDof, iGauss)
                     end do ! iDof numDofDamage
                     PetscCall(DMPlexVecRestoreClosure(dmDamage, PETSC_NULL_SECTION, MEF90DefMechCtx%damageLocal, setPointID(cell), PETSC_NULL_INTEGER, damageDof, ierr))
                  end if

                  totalStrainGauss = 0.0_kr
                  PetscCall(DMPlexVecGetClosure(dmDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%displacementLocal, setPointID(cell), PETSC_NULL_INTEGER, displacementDof, ierr))
                  do iDof = 1, numDofDisplacement
                     totalStrainGauss = totalStrainGauss + displacementDof(iDof) * elemVect(cell)%GradS_BF(iDof, iGauss)
                  end do ! iDof numDofDisplacement
                  PetscCall(DMPlexVecRestoreClosure(dmDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%displacementLocal, setPointID(cell), PETSC_NULL_INTEGER, displacementDof, ierr))

                  temperatureGauss = 0.0_kr
                  PetscCall(DMPlexVecGetClosure(dmTemperature, PETSC_NULL_SECTION, MEF90DefMechCtx%temperatureLocal, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
                  do iDof = 1, numDofDamage
                     temperatureGauss = temperatureGauss + temperatureDof(iDof) * elemScal(cell)%BF(iDof, iGauss)
                  end do ! iDof numDofDamage
                  PetscCall(DMPlexVecRestoreClosure(dmTemperature, PETSC_NULL_SECTION, MEF90DefMechCtx%temperatureLocal, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
                  totalStrainGauss = totalStrainGauss - (temperatureGauss * matpropSet%LinearThermalExpansion)

! #if MEF90_DIM == 2
!!! We need something along these lines
!!! Adding terms in planestrain for plasticity with tr(p) = 0
! If (.NOT. matProp%HookesLaw%isPlaneStress) Then
!    stressGauss = stressGauss +  stiffness * ( matProp%HookesLaw%lambda*trace(plasticStrainCell)*MEF90MatS2DIdentity )
! End If
! #endif

                  PetscCall(PetscSectionGetOffset(sectionPlasticStrain, setPointID(cell), vecOffset, ierr))
                  plasticStrainCell = plasticStrainArray(vecOffset + 1:vecOffset + 1 + SIZEOFMEF90_MATS)

                  if (ATModel%isElastic) then
                     AGauss = matpropSet%HookesLaw
                  else
                     call Split%D2EED(totalStrainGauss - plasticStrainCell, matpropSet%HookesLaw, AGaussPlus, AGaussMinus)
                     AGauss = (ATModel%a(damageGauss) + matpropSet%residualStiffness) * AGaussPlus + AGaussMinus
                  end if
                  do jDof = 0, numDofDisplacement - 1
                     do iDof = 1, numDofDisplacement
                        AGradS_BF = AGauss * elemVect(cell)%GradS_BF(iDof, iGauss)
                        matDof(jDof * numDofDisplacement + iDof) = matDof(jDof * numDofDisplacement + iDof) + elemVect(cell)%Gauss_C(iGauss) * (AGradS_BF .DotP. elemVect(cell)%GradS_BF(jDof + 1, iGauss))
                     end do ! jDof numDofDisplacement
                  end do ! iDof numDofDisplacement
               end do ! iGauss
               PetscCall(DMPlexMatSetClosure(dmDisplacement, PETSC_NULL_SECTION, PETSC_NULL_SECTION, A, setPointID(cell), matDof, ADD_VALUES, ierr))
            end do ! cell

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
            if ((norm2(cellSetOptions%cohesiveDisplacement) /= 0.0_kr) .and. (matpropSet%cohesiveStiffness /= 0.0_kr)) then
               do cell = 1, size(setPointID)
                  matDof = 0.0_kr
                  do iGauss = 1, numGauss
                     do jDof = 0, numDofDisplacement - 1
                        do iDof = 1, numDofDisplacement
                           U0Gauss = matpropSet%cohesiveStiffness * elemVect(cell)%BF(iDof, iGauss)
                           matDof(jDof * numDofDisplacement + iDof) = matDof(jDof * numDofDisplacement + iDof) + elemVect(cell)%Gauss_C(iGauss) * (U0Gauss .DotP. elemVect(cell)%BF(jDof + 1, iGauss))
                        end do ! iDof numDofDisplacement
                     end do ! jDof numDofDisplacement
                  end do ! iGauss
                  PetscCall(DMPlexMatSetClosure(dmDisplacement, PETSC_NULL_SECTION, PETSC_NULL_SECTION, A, setPointID(cell), matDof, ADD_VALUES, ierr))
               end do ! cell
            end if ! cohesiveDisplacement

            deallocate (matDof)

            PetscCall(MEF90ElementDestroy(elemVect, ierr))
            PetscCall(MEF90ElementDestroy(elemScal, ierr))
            PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
            PetscCall(ISDestroy(setPointIS, ierr))
         end if ! setPointIS
      end do ! set
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(VecRestoreArray(MEF90DefMechCtx%plasticStrain, plasticStrainArray, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end if ! setIS
   PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr))
   PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr))
   PetscCall(MatCopy(A, M, SAME_NONZERO_PATTERN, ierr))
end subroutine MEF90DefMechBilinearFormDisplacement

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

subroutine MEF90DefMechWork(MEF90DefMechCtx, bodyForceWork, boundaryForceWork, ierr)
   type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
   PetscReal, dimension(:), pointer                    :: bodyForceWork, boundaryForceWork
   PetscErrorCode, intent(INOUT)                       :: ierr

   type(tDM)                                           :: dmDisplacement, dmBodyForce, dmBoundaryForce, dmPressureForce
   type(tPetscSection)                                 :: sectionBodyForce, sectionBoundaryForce, sectionPressureForce
   PetscReal, dimension(:), pointer                    :: displacementDof, bodyForceArray, boundaryForceArray, pressureForceArray
   type(tIS)                                           :: setIS, setPointIS
   PetscInt, dimension(:), pointer                     :: setID, setPointID
   PetscInt                                            :: set, QuadratureOrder, cell, dim, vecOffset
   type(MEF90DefMechCellSetOptions_Type), pointer      :: cellSetOptions
   type(MEF90DefMechFaceSetOptions_Type), pointer      :: faceSetOptions
   type(MEF90_ELEMENT_ELAST), dimension(:), pointer    :: elemVect
   type(eDMPolytopeType)                               :: cellGeometry
   type(MEF90ElementType)                              :: elemVectType

   type(MEF90CtxGlobalOptions_Type), pointer           :: MEF90CtxGlobalOptions
   type(MEF90DefMechGlobalOptions_Type), pointer       :: MEF90DefMechGlobalOptions
   type(MEF90_VECT)                                    :: bodyForce, boundaryForce, pressureForce
   PetscReal                                           :: myWork
   PetscInt                                            :: iDof, iGauss, numDofDisplacement, numGauss

   PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag, MEF90CtxGlobalOptions, ierr))
   PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag, MEF90DefMechGlobalOptions, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%displacementLocal, dmDisplacement, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%bodyForce, dmBodyForce, ierr))
   PetscCall(DMGetLocalSection(dmBodyForce, sectionBodyForce, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%boundaryForce, dmBoundaryForce, ierr))
   PetscCall(DMGetLocalSection(dmBoundaryForce, sectionBoundaryForce, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%pressureForce, dmPressureForce, ierr))
   PetscCall(DMGetLocalSection(dmPressureForce, sectionPressureForce, ierr))

   PetscCall(DMGetDimension(dmDisplacement, dim, ierr))

   !!! get IS for cell sets
   PetscCall(DMGetLabelIdIS(dmDisplacement, MEF90CellSetLabelName, setIS, ierr))
   PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm, setIS, ierr))

   bodyForceWork = 0.0_kr
   if (.not. PetscObjectIsNull(setIS)) then
      PetscCall(ISGetIndices(setIS, setID, ierr))
      PetscCall(VecGetArray(MEF90DefMechCtx%bodyForce, bodyForceArray, ierr))
      do set = 1, size(setID)
         myWork = 0.0_kr
         PetscCall(DMGetStratumIS(dmDisplacement, MEF90CellSetLabelName, setID(set), setPointIS, ierr))
         if (.not. PetscObjectIsNull(setPointIS)) then
            PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set), cellSetOptions, ierr))

            PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
            PetscCall(DMPlexGetCellType(dmDisplacement, setPointID(1), cellGeometry, ierr))
            PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellGeometry, elemVectType, ierr))

            !!! Allocate elements
            QuadratureOrder = 2 * (elemVectType%order - 1)
            PetscCall(MEF90ElementCreate(dmDisplacement, setPointIS, elemVect, QuadratureOrder, elemVectType, ierr))

            numDofDisplacement = size(elemVect(1)%BF(:, 1))
            numGauss = size(elemVect(1)%Gauss_C)

            if (norm2(cellSetOptions%bodyForce) /= 0.0_kr) then
               do cell = 1, size(setPointID)
                  PetscCall(PetscSectionGetOffset(sectionBodyForce, setPointID(cell), vecOffset, ierr))
                  bodyForce = bodyForceArray(vecOffset + 1:vecOffset + 1 + SIZEOFMEF90_VECT)
                  PetscCall(DMPlexVecGetClosure(dmDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%displacementLocal, setPointID(cell), PETSC_NULL_INTEGER, displacementDof, ierr))
                  do iGauss = 1, numGauss
                     do iDof = 1, numDofDisplacement
                        myWork = myWork + elemVect(cell)%Gauss_C(iGauss) * displacementDof(iDof) * (bodyForce .DotP. elemVect(cell)%BF(iDof, iGauss))
                     end do ! iDof numDofDisplacement
                  end do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%displacementLocal, setPointID(cell), PETSC_NULL_INTEGER, displacementDof, ierr))
               end do ! cell
            end if ! cellSetOptions%bodyForce

            PetscCall(MEF90ElementDestroy(elemVect, ierr))
            PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
            PetscCall(ISDestroy(setPointIS, ierr))
         end if ! setPointIS
         PetscCallMPI(MPI_AllReduce(myWork, bodyForceWork(set), 1, MPIU_SCALAR, MPI_SUM, MEF90DefMechCtx%MEF90Ctx%comm, ierr))
      end do ! set
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(VecRestoreArray(MEF90DefMechCtx%bodyForce, bodyForceArray, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end if ! setIS

   !!! face-based contributions
   PetscCall(DMGetLabelIdIS(dmDisplacement, MEF90FaceSetLabelName, setIS, ierr))
   PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm, setIS, ierr))

   boundaryForceWork = 0.0_kr
   if (.not. PetscObjectIsNull(setIS)) then
      PetscCall(ISGetIndices(setIS, setID, ierr))
      PetscCall(VecGetArray(MEF90DefMechCtx%boundaryForce, boundaryForceArray, ierr))
      PetscCall(VecGetArray(MEF90DefMechCtx%pressureForce, pressureForceArray, ierr))
      do set = 1, size(setID)
         myWork = 0.0_kr
         PetscCall(DMGetStratumIS(dmDisplacement, MEF90FaceSetLabelName, setID(set), setPointIS, ierr))
         if (.not. PetscObjectIsNull(setPointIS)) then
            PetscCall(PetscBagGetDataMEF90DefMechCtxFaceSetOptions(MEF90DefMechCtx%FaceSetOptionsBag(set), faceSetOptions, ierr))
            PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
            PetscCall(DMPlexGetCellType(dmDisplacement, setPointID(1), cellGeometry, ierr))
            PetscCall(MEF90ElementGetTypeBoundary(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellGeometry, elemVectType, ierr))

            !!! Allocate elements
            QuadratureOrder = 2 * elemVectType%order
            PetscCall(MEF90ElementCreate(dmDisplacement, setPointIS, elemVect, QuadratureOrder, elemVectType, ierr))

            numDofDisplacement = size(elemVect(1)%BF(:, 1))
            numGauss = size(elemVect(1)%Gauss_C)

            if (norm2(faceSetOptions%boundaryForce) /= 0.0_kr) then
               do cell = 1, size(setPointID)
                  PetscCall(PetscSectionGetOffset(sectionBoundaryForce, setPointID(cell), vecOffset, ierr))
                  boundaryForce = boundaryForceArray(vecOffset + 1:vecOffset + 1 + SIZEOFMEF90_VECT)
                  PetscCall(DMPlexVecGetClosure(dmDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%displacementLocal, setPointID(cell), PETSC_NULL_INTEGER, displacementDof, ierr))
                  do iGauss = 1, numGauss
                     do iDof = 1, numDofDisplacement
                        myWork = myWork + elemVect(cell)%Gauss_C(iGauss) * displacementDof(iDof) * (boundaryForce .DotP. elemVect(cell)%BF(iDof, iGauss))
                     end do ! iDof numDofDisplacement
                  end do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%displacementLocal, setPointID(cell), PETSC_NULL_INTEGER, displacementDof, ierr))
               end do ! cell
            end if ! faceSetOptions%boundaryForce

            if (faceSetOptions%pressureForce /= 0.0_kr) then
               do cell = 1, size(setPointID)
                  PetscCall(PetscSectionGetOffset(sectionPressureForce, setPointID(cell), vecOffset, ierr))
                  PetscCall(DMPlexVecGetClosure(dmDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%displacementLocal, setPointID(cell), PETSC_NULL_INTEGER, displacementDof, ierr))
                  pressureForce = boundaryForceArray(vecOffset + 1) * elemVect(cell)%outerNormal
                  do iGauss = 1, numGauss
                     do iDof = 1, numDofDisplacement
                        myWork = myWork + elemVect(cell)%Gauss_C(iGauss) * displacementDof(iDof) * (pressureForce .DotP. elemVect(cell)%BF(iDof, iGauss))
                     end do ! iDof numDofDisplacement
                  end do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%displacementLocal, setPointID(cell), PETSC_NULL_INTEGER, displacementDof, ierr))
               end do ! cell
            end if ! faceSetOptions%pressureForce

            PetscCall(MEF90ElementDestroy(elemVect, ierr))
            PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
            PetscCall(ISDestroy(setPointIS, ierr))
         end if ! setPointIS
         PetscCallMPI(MPI_AllReduce(myWork, boundaryForceWork(set), 1, MPIU_SCALAR, MPI_SUM, MEF90DefMechCtx%MEF90Ctx%comm, ierr))
      end do ! set
      PetscCall(VecRestoreArray(MEF90DefMechCtx%boundaryForce, boundaryForceArray, ierr))
      PetscCall(VecRestoreArray(MEF90DefMechCtx%pressureForce, pressureForceArray, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end if ! setIS

end subroutine MEF90DefMechWork

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

subroutine MEF90DefMechCohesiveEnergy(MEF90DefMechCtx, cohesiveEnergy, ierr)
   type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
   PetscReal, dimension(:), pointer                    :: cohesiveEnergy
   PetscErrorCode, intent(INOUT)                       :: ierr

   type(tDM)                                           :: dmDisplacement, dmCohesiveDisplacement
   PetscReal, dimension(:), pointer                    :: displacementDof, cohesiveDisplacementDof
   type(tIS)                                           :: setIS, setPointIS
   PetscInt, dimension(:), pointer                     :: setID, setPointID
   PetscInt                                            :: set, QuadratureOrder, cell
   type(MEF90_MATPROP), pointer                        :: matpropSet
   type(MEF90DefMechCellSetOptions_Type), pointer      :: cellSetOptions
   type(MEF90_ELEMENT_ELAST), dimension(:), pointer    :: elemVect
   type(eDMPolytopeType)                               :: cellGeometry
   type(MEF90ElementType)                              :: elemVectType

   type(MEF90CtxGlobalOptions_Type), pointer           :: MEF90CtxGlobalOptions
   type(MEF90DefMechGlobalOptions_Type), pointer       :: MEF90DefMechGlobalOptions
   type(MEF90_VECT)                                    :: U0Gauss
   PetscReal                                           :: myCohesiveEnergy
   PetscInt                                            :: iDof, iGauss, numDofDisplacement, numGauss

   PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag, MEF90CtxGlobalOptions, ierr))
   PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag, MEF90DefMechGlobalOptions, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%displacementLocal, dmDisplacement, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%cohesiveDisplacement, dmCohesiveDisplacement, ierr))

   !!! get IS for cell sets
   PetscCall(DMGetLabelIdIS(dmDisplacement, MEF90CellSetLabelName, setIS, ierr))
   PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm, setIS, ierr))

   cohesiveEnergy = 0.0_kr
   if (.not. PetscObjectIsNull(setIS)) then
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         myCohesiveEnergy = 0.0_kr
         PetscCall(DMGetStratumIS(dmDisplacement, MEF90CellSetLabelName, setID(set), setPointIS, ierr))
         if (.not. PetscObjectIsNull(setPointIS)) then
            PetscCall(PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set), matpropSet, ierr))
            PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set), cellSetOptions, ierr))

            PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
            PetscCall(DMPlexGetCellType(dmDisplacement, setPointID(1), cellGeometry, ierr))
            PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellGeometry, elemVectType, ierr))

            !!! Allocate elements
            QuadratureOrder = 2 * (elemVectType%order - 1)
            PetscCall(MEF90ElementCreate(dmDisplacement, setPointIS, elemVect, QuadratureOrder, elemVectType, ierr))

            numDofDisplacement = size(elemVect(1)%BF(:, 1))
            numGauss = size(elemVect(1)%Gauss_C)

            if ((norm2(cellSetOptions%cohesiveDisplacement) /= 0.0_kr) .and. (matpropSet%cohesiveStiffness /= 0.0_kr)) then
               do cell = 1, size(setPointID)
                  U0Gauss = 0.0_kr
                  PetscCall(DMPlexVecGetClosure(dmDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%displacementLocal, setPointID(cell), PETSC_NULL_INTEGER, displacementDof, ierr))
                  PetscCall(DMPlexVecGetClosure(dmCohesiveDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%cohesiveDisplacement, setPointID(cell), PETSC_NULL_INTEGER, cohesiveDisplacementDof, ierr))
                  do iGauss = 1, numGauss
                     do iDof = 1, numDofDisplacement
                        U0Gauss = U0Gauss + (displacementDof(iDof) - cohesiveDisplacementDof(iDof)) * elemVect(cell)%BF(iDof, iGauss)
                     end do ! iDof numDofDisplacement
                     do iDof = 1, numDofDisplacement
                        myCohesiveEnergy = myCohesiveEnergy + elemVect(cell)%Gauss_C(iGauss) * matpropSet%cohesiveStiffness * (U0Gauss .DotP. U0Gauss)
                     end do ! iDof numDofDisplacement
                  end do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmCohesiveDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%cohesiveDisplacement, setPointID(cell), PETSC_NULL_INTEGER, cohesiveDisplacementDof, ierr))
                  PetscCall(DMPlexVecRestoreClosure(dmDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%displacementLocal, setPointID(cell), PETSC_NULL_INTEGER, displacementDof, ierr))
               end do ! cell
            end if ! cohesiveDisplacement

            PetscCall(MEF90ElementDestroy(elemVect, ierr))
            PetscCall(ISDestroy(setPointIS, ierr))
         end if ! setPointIS
         PetscCallMPI(MPI_AllReduce(myCohesiveEnergy, cohesiveEnergy(set), 1, MPIU_SCALAR, MPI_SUM, MEF90DefMechCtx%MEF90Ctx%comm, ierr))
      end do ! set
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end if ! setIS

end subroutine MEF90DefMechCohesiveEnergy

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechPlasticDissipation"
!!!
!!!
!!!  MEF90DefMechPlasticDissipation:
!!!
!!!  (c) 2015 Erwan TANNE erwan.tanne@gmail.com
!!!

subroutine MEF90DefMechPlasticDissipation(MEF90DefMechCtx, plasticStrainOld, energy, ierr)
   type(tVec), intent(IN)                              :: plasticStrainOld
   type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
   PetscReal, dimension(:), pointer                    :: energy
   PetscErrorCode, intent(INOUT)                       :: ierr

end subroutine MEF90DefMechPlasticDissipation

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

subroutine MEF90DefMechElasticEnergy(MEF90DefMechCtx, energy, ierr)
   type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
   PetscReal, dimension(:), pointer                    :: energy
   PetscErrorCode, intent(INOUT)                       :: ierr

   type(tDM)                                           :: dmDisplacement, dmDamage, dmTemperature, dmPlasticStrain
   type(tPetscSection)                                 :: sectionPlasticStrain
   PetscReal, dimension(:), pointer                    :: displacementDof, damageDof, temperatureDof, plasticStrainArray
   PetscInt                                            :: vecOffset
   type(MEF90_MATS)                                    :: plasticStrainCell
   type(tIS)                                           :: setIS, setPointIS
   PetscInt, dimension(:), pointer                     :: setID, setPointID
   PetscInt                                            :: set, QuadratureOrder, cell, dim
   type(MEF90_MATPROP), pointer                        :: matpropSet
   type(MEF90DefMechCellSetOptions_Type), pointer      :: cellSetOptions
   type(MEF90_ELEMENT_ELAST), dimension(:), pointer    :: elemVect
   type(MEF90_ELEMENT_SCAL), dimension(:), pointer     :: elemScal
   type(eDMPolytopeType)                               :: cellGeometryVect, cellGeometryScal
   type(MEF90ElementType)                              :: elemVectType, elemScalType

   type(MEF90CtxGlobalOptions_Type), pointer           :: MEF90CtxGlobalOptions
   type(MEF90DefMechGlobalOptions_Type), pointer       :: MEF90DefMechGlobalOptions
   class(MEF90DefMechAT_Type), allocatable             :: ATModel
   class(MEF90_DEFMECHSPLIT), allocatable              :: Split
   type(MEF90_MATS)                                    :: totalStrainGauss
   PetscReal                                           :: damageGauss, temperatureGauss, myEnergy, EEDPlus, EEDMinus, elasticEnergyDensityGauss
   PetscInt                                            :: iDof, iGauss, numDofDisplacement, numDofDamage, numGauss
   character(len=MEF90MXSTRLEN)                        :: prefix

   PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag, MEF90CtxGlobalOptions, ierr))
   PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag, MEF90DefMechGlobalOptions, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%displacementLocal, dmDisplacement, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%damageLocal, dmDamage, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%TemperatureLocal, dmTemperature, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%plasticStrain, dmPlasticStrain, ierr))
   PetscCall(DMGetLocalSection(dmPlasticStrain, sectionPlasticStrain, ierr))

   PetscCall(DMGetDimension(dmDisplacement, dim, ierr))

   !!! get IS for cell sets
   PetscCall(DMGetLabelIdIS(dmDisplacement, MEF90CellSetLabelName, setIS, ierr))
   PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm, setIS, ierr))

   energy = 0.0_kr
   if (.not. PetscObjectIsNull(setIS)) then
      PetscCall(VecGetArray(MEF90DefMechCtx%plasticStrain, plasticStrainArray, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         myEnergy = 0.0_kr
         PetscCall(DMGetStratumIS(dmDisplacement, MEF90CellSetLabelName, setID(set), setPointIS, ierr))
         if (.not. PetscObjectIsNull(setPointIS)) then
            PetscCall(PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set), matpropSet, ierr))
            PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set), cellSetOptions, ierr))

            PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
            PetscCall(DMPlexGetCellType(dmDisplacement, setPointID(1), cellGeometryVect, ierr))
            PetscCall(DMPlexGetCellType(dmDamage, setPointID(1), cellGeometryScal, ierr))
            PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellGeometryVect, elemVectType, ierr))
            PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellGeometryScal, elemScalType, ierr))

            !!! get the ATModel and split objects
            write(prefix,'("cs",I4.4,"_")') setID(set)
            PetscCall(MEF90DefMechGetATModel(MEF90DefMechCtx%MEF90Ctx%comm, prefix, ATModel, ierr))
            PetscCall(ATModel%setFromOptions(ierr))

            if (cellSetOptions%unilateralContactHybrid) then
               Split = MEF90_DEFMECHSPLITNONE()
            else
               PetscCall(MEF90DefMechGetSplit(cellSetOptions, Split, ierr))
            end if

            !!! Allocate elements
            QuadratureOrder = max(2 * elemVectType%order, split%damageOrder + split%strainOrder)
            PetscCall(MEF90ElementCreate(dmDisplacement, setPointIS, elemVect, QuadratureOrder, elemVectType, ierr))
            PetscCall(MEF90ElementCreate(dmDamage, setPointIS, elemScal, QuadratureOrder, elemScalType, ierr))

            numDofDisplacement = size(elemVect(1)%BF(:, 1))
            numDofDamage = size(elemScal(1)%BF(:, 1))
            numGauss = size(elemVect(1)%Gauss_C)
            do cell = 1, size(setPointID)
               PetscCall(DMPlexVecGetClosure(dmDamage, PETSC_NULL_SECTION, MEF90DefMechCtx%damageLocal, setPointID(cell), PETSC_NULL_INTEGER, damageDof, ierr))
               PetscCall(DMPlexVecGetClosure(dmDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%displacementLocal, setPointID(cell), PETSC_NULL_INTEGER, displacementDof, ierr))
               PetscCall(DMPlexVecGetClosure(dmTemperature, PETSC_NULL_SECTION, MEF90DefMechCtx%TemperatureLocal, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))

               PetscCall(PetscSectionGetOffset(sectionPlasticStrain, setPointID(cell), vecOffset, ierr))
               plasticStrainCell = plasticStrainArray(vecOffset + 1:vecOffset + 1 + SIZEOFMEF90_MATS)
               do iGauss = 1, numGauss
                  damageGauss = 0.0_kr
                  if (.not. ATModel%isElastic) then
                     do iDof = 1, numDofDamage
                        damageGauss = damageGauss + damageDof(iDof) * elemScal(cell)%BF(iDof, iGauss)
                     end do ! iDof numDofDamage
                  end if

                  totalStrainGauss = 0.0_kr
                  do iDof = 1, numDofDisplacement
                     totalStrainGauss = totalStrainGauss + displacementDof(iDof) * elemVect(cell)%GradS_BF(iDof, iGauss)
                  end do ! iDof numDofDisplacement

                  temperatureGauss = 0.0_kr
                  do iDof = 1, numDofDamage
                     temperatureGauss = temperatureGauss + temperatureDof(iDof) * elemScal(cell)%BF(iDof, iGauss)
                  end do ! iDof numDofDamage
                  totalStrainGauss = totalStrainGauss - (temperatureGauss * matpropSet%LinearThermalExpansion)

! #if MEF90_DIM == 2
!!! We need something along these lines
!!! Adding terms in planestrain for plasticity with tr(p) = 0
! If (.NOT. matProp%HookesLaw%isPlaneStress) Then
!    stressGauss = stressGauss +  stiffness * ( matProp%HookesLaw%lambda*trace(plasticStrainCell)*MEF90MatS2DIdentity )
! End If
! #endif

                  call Split%EED(totalStrainGauss - plasticStrainCell, matpropSet%HookesLaw, EEDPlus, EEDMinus)
                  if (ATModel%isElastic) then
                     elasticEnergyDensityGauss = EEDPlus + EEDMinus
                  else
                     elasticEnergyDensityGauss = (ATmodel%a(damageGauss) + matpropSet%residualStiffness) * EEDPlus + EEDMinus
                  end if

                  myEnergy = myEnergy + elemVect(cell)%Gauss_C(iGauss) * elasticEnergyDensityGauss
               end do ! iGauss

               PetscCall(DMPlexVecRestoreClosure(dmDamage, PETSC_NULL_SECTION, MEF90DefMechCtx%damageLocal, setPointID(cell), PETSC_NULL_INTEGER, damageDof, ierr))
               PetscCall(DMPlexVecRestoreClosure(dmDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%displacementLocal, setPointID(cell), PETSC_NULL_INTEGER, displacementDof, ierr))
               PetscCall(DMPlexVecRestoreClosure(dmTemperature, PETSC_NULL_SECTION, MEF90DefMechCtx%TemperatureLocal, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
            end do ! cell
            PetscCall(MEF90ElementDestroy(elemVect, ierr))
            PetscCall(MEF90ElementDestroy(elemScal, ierr))
            PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
            PetscCall(ISDestroy(setPointIS, ierr))
         end if ! setPointIS
         PetscCallMPI(MPI_AllReduce(myEnergy, energy(set), 1, MPIU_SCALAR, MPI_SUM, MEF90DefMechCtx%MEF90Ctx%comm, ierr))
      end do ! set
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(VecRestoreArray(MEF90DefMechCtx%plasticStrain, plasticStrainArray, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end if ! setIS
end subroutine MEF90DefMechElasticEnergy

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

subroutine MEF90DefMechStress(MEF90DefMechCtx, stress, ierr)
   type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
   PetscErrorCode, intent(INOUT)                       :: ierr
   type(tVec), intent(IN)                              :: stress

   type(tDM)                                           :: dmDisplacement, dmDamage, dmTemperature, dmPlasticStrain, dmStress
   type(tPetscSection)                                 :: sectionPlasticStrain
   PetscReal, dimension(:), pointer                    :: displacementDof, damageDof, temperatureDof, plasticStrainArray
   type(MEF90_MATS)                                    :: plasticStrainCell
   type(tIS)                                           :: setIS, setPointIS
   PetscInt, dimension(:), pointer                     :: setID, setPointID
   PetscInt                                            :: set, QuadratureOrder, cell, dim, vecOffset
   type(MEF90_MATPROP), pointer                        :: matpropSet
   type(MEF90DefMechCellSetOptions_Type), pointer      :: cellSetOptions
   type(MEF90_ELEMENT_ELAST), dimension(:), pointer    :: elemVect
   type(MEF90_ELEMENT_SCAL), dimension(:), pointer     :: elemScal
   type(eDMPolytopeType)                               :: cellGeometryVect, cellGeometryScal
   type(MEF90ElementType)                              :: elemVectType, elemScalType
   PetscReal, dimension(:), pointer                    :: stressDof

   type(MEF90CtxGlobalOptions_Type), pointer           :: MEF90CtxGlobalOptions
   type(MEF90DefMechGlobalOptions_Type), pointer       :: MEF90DefMechGlobalOptions
   class(MEF90DefMechAT_Type), allocatable             :: ATModel
   class(MEF90_DEFMECHSPLIT), allocatable              :: Split
   type(MEF90_MATS)                                    :: totalStrainGauss, stressGaussPlus, stressGaussMinus, stressCell
   PetscReal                                           :: damageGauss, temperatureGauss, cellSize
   PetscInt                                            :: iDof, iGauss, numDofDisplacement, numDofDamage, numGauss
   character(len=MEF90MXSTRLEN)                        :: prefix

   PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag, MEF90CtxGlobalOptions, ierr))
   PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag, MEF90DefMechGlobalOptions, ierr))
   PetscCall(VecGetDM(stress, dmStress, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%displacementLocal, dmDisplacement, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%damageLocal, dmDamage, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%TemperatureLocal, dmTemperature, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%plasticStrain, dmPlasticStrain, ierr))
   PetscCall(DMGetLocalSection(dmPlasticStrain, sectionPlasticStrain, ierr))

   PetscCall(DMGetDimension(dmDisplacement, dim, ierr))

   PetscCall(VecSet(stress, 0.0_kr, ierr))

   !!! get IS for cell sets
   PetscCall(DMGetLabelIdIS(dmDisplacement, MEF90CellSetLabelName, setIS, ierr))
   PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm, setIS, ierr))

   if (.not. PetscObjectIsNull(setIS)) then
      PetscCall(ISGetIndices(setIS, setID, ierr))
      PetscCall(VecGetArray(MEF90DefMechCtx%plasticStrain, plasticStrainArray, ierr))
      do set = 1, size(setID)
         PetscCall(DMGetStratumIS(dmDisplacement, MEF90CellSetLabelName, setID(set), setPointIS, ierr))
         if (.not. PetscObjectIsNull(setPointIS)) then
            PetscCall(PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set), matpropSet, ierr))
            PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set), cellSetOptions, ierr))

            PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
            PetscCall(DMPlexGetCellType(dmDisplacement, setPointID(1), cellGeometryVect, ierr))
            PetscCall(DMPlexGetCellType(dmDamage, setPointID(1), cellGeometryScal, ierr))
            PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellGeometryVect, elemVectType, ierr))
            PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellGeometryScal, elemScalType, ierr))

            !!! get the ATModel and split objects
            write(prefix,'("cs",I4.4,"_")') setID(set)
            PetscCall(MEF90DefMechGetATModel(MEF90DefMechCtx%MEF90Ctx%comm, prefix, ATModel, ierr))
            PetscCall(ATModel%setFromOptions(ierr))

            if (cellSetOptions%unilateralContactHybrid) then
               Split = MEF90_DEFMECHSPLITNONE()
            else
               PetscCall(MEF90DefMechGetSplit(cellSetOptions, Split, ierr))
            end if

            !!! Allocate elements
            QuadratureOrder = max(2 * elemVectType%order, split%damageOrder + split%strainOrder)
            PetscCall(MEF90ElementCreate(dmDisplacement, setPointIS, elemVect, QuadratureOrder, elemVectType, ierr))
            PetscCall(MEF90ElementCreate(dmDamage, setPointIS, elemScal, QuadratureOrder, elemScalType, ierr))

            numDofDisplacement = size(elemVect(1)%BF(:, 1))
            numDofDamage = size(elemScal(1)%BF(:, 1))
            numGauss = size(elemVect(1)%Gauss_C)

            allocate (stressDof(SIZEOFMEF90_MATS))

            do cell = 1, size(setPointID)
               PetscCall(DMPlexVecGetClosure(dmDamage, PETSC_NULL_SECTION, MEF90DefMechCtx%damageLocal, setPointID(cell), PETSC_NULL_INTEGER, damageDof, ierr))
               PetscCall(DMPlexVecGetClosure(dmDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%displacementLocal, setPointID(cell), PETSC_NULL_INTEGER, displacementDof, ierr))
               PetscCall(DMPlexVecGetClosure(dmTemperature, PETSC_NULL_SECTION, MEF90DefMechCtx%TemperatureLocal, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
               stressCell = 0.0_kr
               cellSize = 0.0_kr
               do iGauss = 1, numGauss
                  damageGauss = 0.0_kr
                  if (.not. ATModel%isElastic) then
                     do iDof = 1, numDofDamage
                        damageGauss = damageGauss + damageDof(iDof) * elemScal(cell)%BF(iDof, iGauss)
                     end do ! iDof numDofDamage
                  end if

                  totalStrainGauss = 0.0_kr
                  do iDof = 1, numDofDisplacement
                     totalStrainGauss = totalStrainGauss + displacementDof(iDof) * elemVect(cell)%GradS_BF(iDof, iGauss)
                  end do ! iDof numDofDisplacement

                  temperatureGauss = 0.0_kr
                  do iDof = 1, numDofDamage
                     temperatureGauss = temperatureGauss + temperatureDof(iDof) * elemScal(cell)%BF(iDof, iGauss)
                  end do ! iDof numDofDamage
                  totalStrainGauss = totalStrainGauss - (temperatureGauss * matpropSet%LinearThermalExpansion)

! #if MEF90_DIM == 2
!!! We need something along these lines
!!! Adding terms in planestrain for plasticity with tr(p) = 0
! If (.NOT. matProp%HookesLaw%isPlaneStress) Then
!    stressGauss = stressGauss +  stiffness * ( matProp%HookesLaw%lambda*trace(plasticStrainCell)*MEF90MatS2DIdentity )
! End If
! #endif

                  PetscCall(PetscSectionGetOffset(sectionPlasticStrain, setPointID(cell), vecOffset, ierr))
                  plasticStrainCell = plasticStrainArray(vecOffset + 1:vecOffset + 1 + SIZEOFMEF90_MATS)

                  call Split%DEED(totalStrainGauss - plasticStrainCell, matpropSet%HookesLaw, stressGaussPlus, stressGaussMinus)
                  if (ATModel%isElastic) then
                     stressCell = stressCell + elemVect(cell)%Gauss_C(iGauss) * (stressGaussPlus + stressGaussMinus)
                  else
                     stressCell = stressCell + elemVect(cell)%Gauss_C(iGauss) * ((ATModel%a(damageGauss) + matpropSet%residualStiffness) * stressGaussPlus + stressGaussMinus)
                  end if
                  cellSize = cellSize + elemVect(cell)%Gauss_C(iGauss)
               end do ! iGauss
               stressCell = stressCell / cellSize
               stressDof = stressCell
               PetscCall(DMPlexVecSetClosure(dmStress, PETSC_NULL_SECTION, stress, setPointID(cell), stressDof, INSERT_VALUES, ierr))
               PetscCall(DMPlexVecRestoreClosure(dmDamage, PETSC_NULL_SECTION, MEF90DefMechCtx%damageLocal, setPointID(cell), PETSC_NULL_INTEGER, damageDof, ierr))
               PetscCall(DMPlexVecRestoreClosure(dmDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%displacementLocal, setPointID(cell), PETSC_NULL_INTEGER, displacementDof, ierr))
               PetscCall(DMPlexVecRestoreClosure(dmTemperature, PETSC_NULL_SECTION, MEF90DefMechCtx%TemperatureLocal, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
            end do ! cell
            deallocate (stressDof)

            PetscCall(MEF90ElementDestroy(elemVect, ierr))
            PetscCall(MEF90ElementDestroy(elemScal, ierr))
            PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
            PetscCall(ISDestroy(setPointIS, ierr))
         end if ! setPointIS
      end do ! set
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(VecRestoreArray(MEF90DefMechCtx%plasticStrain, plasticStrainArray, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end if ! setIS
end subroutine MEF90DefMechStress

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

subroutine MEF90DefMechOperatorDamage(snesDamage, damage, residual, MEF90DefMechCtx, ierr)
   type(tSNES), intent(IN)                             :: snesDamage
   type(tVec), intent(IN)                              :: damage
   type(tVec), intent(INOUT)                           :: residual
   type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
   PetscErrorCode, intent(OUT)                         :: ierr

   type(tDM)                                           :: dmDisplacement, dmDamage, dmTemperature, dmCohesiveDisplacement, dmPlasticStrain
   type(tPetscSection)                                 :: sectionPlasticStrain
   PetscReal, dimension(:), pointer                    :: displacementDof, damageDof, temperatureDof, plasticStrainArray
   type(MEF90_MATS)                                    :: plasticStrainCell
   type(tIS)                                           :: setIS, setPointIS
   PetscInt, dimension(:), pointer                     :: setID, setPointID
   PetscInt                                            :: set, QuadratureOrder, cell, dim, vecOffset
   type(MEF90_MATPROP), pointer                        :: matpropSet
   type(MEF90DefMechCellSetOptions_Type), pointer      :: cellSetOptions
   type(MEF90_ELEMENT_ELAST), dimension(:), pointer    :: elemVect
   type(MEF90_ELEMENT_SCAL), dimension(:), pointer     :: elemScal
   type(eDMPolytopeType)                               :: cellGeometryVect, cellGeometryScal
   type(MEF90ElementType)                              :: elemVectType, elemScalType
   PetscReal, dimension(:), pointer                    :: residualDof

   type(MEF90CtxGlobalOptions_Type), pointer           :: MEF90CtxGlobalOptions
   type(MEF90DefMechGlobalOptions_Type), pointer       :: MEF90DefMechGlobalOptions
   type(tVec)                                          :: locResidual
   class(MEF90DefMechAT_Type), allocatable             :: ATModel
   class(MEF90_DEFMECHSPLIT), allocatable              :: Split
   type(MEF90_MATS)                                    :: totalStrainGauss, C2
   type(MEF90_VECT)                                    :: gradDamageGauss
   PetscReal                                           :: damageGauss, temperatureGauss, EEDGaussMinus, EEDGaussPlus, C1, C3
   PetscInt                                            :: iDof, iGauss, numDofDisplacement, numDofDamage, numGauss
   character(len=MEF90MXSTRLEN)                        :: prefix

   PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag, MEF90CtxGlobalOptions, ierr))
   PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag, MEF90DefMechGlobalOptions, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%displacementLocal, dmDisplacement, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%cohesiveDisplacement, dmCohesiveDisplacement, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%damageLocal, dmDamage, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%TemperatureLocal, dmTemperature, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%plasticStrain, dmPlasticStrain, ierr))
   PetscCall(DMGetLocalSection(dmPlasticStrain, sectionPlasticStrain, ierr))

   PetscCall(DMGetDimension(dmDisplacement, dim, ierr))

   PetscCall(DMGetLocalVector(dmDamage, locResidual, ierr))
   PetscCall(DMGlobalToLocal(dmDamage, damage, INSERT_VALUES, MEF90DefMechCtx%damageLocal, ierr))

   PetscCall(VecSet(residual, 0.0_kr, ierr))
   PetscCall(VecSet(locResidual, 0.0_kr, ierr))

   !!! get IS for cell sets
   PetscCall(DMGetLabelIdIS(dmDamage, MEF90CellSetLabelName, setIS, ierr))
   PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm, setIS, ierr))

   if (.not. PetscObjectIsNull(setIS)) then
      PetscCall(ISGetIndices(setIS, setID, ierr))
      PetscCall(VecGetArray(MEF90DefMechCtx%plasticStrain, plasticStrainArray, ierr))
      do set = 1, size(setID)
         PetscCall(DMGetStratumIS(dmDamage, MEF90CellSetLabelName, setID(set), setPointIS, ierr))
         if (.not. PetscObjectIsNull(setPointIS)) then
            PetscCall(PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set), matpropSet, ierr))
            PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set), cellSetOptions, ierr))

            PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
            PetscCall(DMPlexGetCellType(dmDisplacement, setPointID(1), cellGeometryVect, ierr))
            PetscCall(DMPlexGetCellType(dmDamage, setPointID(1), cellGeometryScal, ierr))
            PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellGeometryVect, elemVectType, ierr))
            PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellGeometryScal, elemScalType, ierr))

            !!! get the ATModel and split objects
            write(prefix,'("cs",I4.4,"_")') setID(set)
            PetscCall(MEF90DefMechGetATModel(MEF90DefMechCtx%MEF90Ctx%comm, prefix, ATModel, ierr))
            PetscCall(ATModel%setFromOptions(ierr))

            if (cellSetOptions%unilateralContactHybrid) then
               Split = MEF90_DEFMECHSPLITNONE()
            else
               PetscCall(MEF90DefMechGetSplit(cellSetOptions, Split, ierr))
            end if
            
            !!! Allocate elements
            QuadratureOrder = max(2 * elemVectType%order, split%damageOrder + split%strainOrder)
            PetscCall(MEF90ElementCreate(dmDisplacement, setPointIS, elemVect, QuadratureOrder, elemVectType, ierr))
            PetscCall(MEF90ElementCreate(dmDamage, setPointIS, elemScal, QuadratureOrder, elemScalType, ierr))

            numDofDisplacement = size(elemVect(1)%BF(:, 1))
            numDofDamage = size(elemScal(1)%BF(:, 1))
            numGauss = size(elemScal(1)%Gauss_C)

            allocate (residualDof(numDofDamage))

            do cell = 1, size(setPointID)
               residualDof = 0.0_kr
               PetscCall(DMPlexVecGetClosure(dmDamage, PETSC_NULL_SECTION, MEF90DefMechCtx%damageLocal, setPointID(cell), PETSC_NULL_INTEGER, damageDof, ierr))
               PetscCall(DMPlexVecGetClosure(dmDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%displacementLocal, setPointID(cell), PETSC_NULL_INTEGER, displacementDof, ierr))
               PetscCall(DMPlexVecGetClosure(dmTemperature, PETSC_NULL_SECTION, MEF90DefMechCtx%TemperatureLocal, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
               do iGauss = 1, numGauss
                  damageGauss = 0.0_kr
                  gradDamageGauss = 0.0_kr
                  if (.not. ATModel%isElastic) then
                     do iDof = 1, numDofDamage
                        damageGauss = damageGauss + damageDof(iDof) * elemScal(cell)%BF(iDof, iGauss)
                        gradDamageGauss = gradDamageGauss + damageDof(iDof) * elemScal(cell)%Grad_BF(iDof, iGauss)
                     end do ! iDof numDofDamage
                  end if

                  totalStrainGauss = 0.0_kr
                  do iDof = 1, numDofDisplacement
                     totalStrainGauss = totalStrainGauss + displacementDof(iDof) * elemVect(cell)%GradS_BF(iDof, iGauss)
                  end do ! iDof numDofDisplacement

                  temperatureGauss = 0.0_kr
                  do iDof = 1, numDofDamage
                     temperatureGauss = temperatureGauss + temperatureDof(iDof) * elemScal(cell)%BF(iDof, iGauss)
                  end do ! iDof numDofDamage
                  totalStrainGauss = totalStrainGauss - (temperatureGauss * matpropSet%LinearThermalExpansion)

                  ! #if MEF90_DIM == 2
   !!! We need something along these lines
   !!! Adding terms in planestrain for plasticity with tr(p) = 0
                  ! If (.NOT. matProp%HookesLaw%isPlaneStress) Then
                  !    stressGauss = stressGauss +  stiffness * ( matProp%HookesLaw%lambda*trace(plasticStrainCell)*MEF90MatS2DIdentity )
                  ! End If
                  ! #endif

                  PetscCall(PetscSectionGetOffset(sectionPlasticStrain, setPointID(cell), vecOffset, ierr))
                  plasticStrainCell = plasticStrainArray(vecOffset + 1:vecOffset + 1 + SIZEOFMEF90_MATS)

                  if (.not. ATModel%isElastic) then
                     call Split%EED(totalStrainGauss - plasticStrainCell, matpropSet%HookesLaw, EEDGaussPlus, EEDGaussMinus)
                  else
                     EEDGaussPlus = 0.0_kr
                  end if
                  C1 = matpropSet%fractureToughness / ATModel%cw * 0.25_kr / ATModel%internalLength
                  C2 = matpropSet%fractureToughness / ATModel%cw * 0.5_kr * ATModel%internalLength * matpropSet%toughnessAnisotropyMatrix
                  C3 = ATModel%Da(damageGauss) * EEDGaussPlus + C1 * ATModel%Dw(damageGauss)
                  do iDof = 1, numDofDamage
                     residualDof(iDof) = residualDof(iDof) + elemScal(cell)%Gauss_C(iGauss) * ( &
                                         C3 * elemScal(cell)%BF(iDof, iGauss) + (C2 * gradDamageGauss .DotP. elemScal(cell)%Grad_BF(iDof, iGauss)))
                  end do ! iDof numDofDamage
               end do ! iGauss
               PetscCall(DMPlexVecRestoreClosure(dmTemperature, PETSC_NULL_SECTION, MEF90DefMechCtx%TemperatureLocal, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
               PetscCall(DMPlexVecRestoreClosure(dmDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%displacementLocal, setPointID(cell), PETSC_NULL_INTEGER, displacementDof, ierr))
               PetscCall(DMPlexVecRestoreClosure(dmDamage, PETSC_NULL_SECTION, MEF90DefMechCtx%damageLocal, setPointID(cell), PETSC_NULL_INTEGER, damageDof, ierr))
               PetscCall(DMPlexVecSetClosure(dmDamage, PETSC_NULL_SECTION, locResidual, setPointID(cell), residualDof, ADD_VALUES, ierr))
            end do ! cell

            ! !!! crack pressure
            ! If (cellSetOptions%crackPressure /= 0.0_Kr)) Then
            !    Do cell = 1,size(setPointID)
            !       residualDof = 0.0_Kr
            !       Do iGauss = 1,numGauss
            !          CrackPressureCell = 0.0_Kr
            !          PetscCall(DMPlexVecGetClosure(dmCrackPressure,PETSC_NULL_SECTION,MEF90DefMechCtx%crackPressure,setPointID(cell),PETSC_NULL_INTEGER,crackPressureDof,ierr))
            !          Do iDof = 1,numDofDamage
            !             CrackPressureCell = CrackPressureCell + crackPressureDof(iDof) * elemScal(cell)%BF(iDof,iGauss)
            !          End Do ! iDof numDofDamage
            !          PetscCall(DMPlexVecRestoreClosure(dmCrackPressure,PETSC_NULL_SECTION,MEF90DefMechCtx%crackPressure,setPointID(cell),PETSC_NULL_INTEGER,crackPressureDof,ierr))
            !          displacementCell = 0.0_Kr
            !          PetscCall(DMPlexVecGetClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),PETSC_NULL_INTEGER,displacementDof,ierr))
            !          Do iDof = 1,numDofDisplacement
            !             displacementCell = displacementCell + displacementDof(iDof) * elemVect(cell)%BF(iDof,iGauss)
            !          End Do ! iDof numDofDisplacement
            !          PetscCall(DMPlexVecRestoreClosure(dmDisplacement,PETSC_NULL_SECTION,MEF90DefMechCtx%displacementLocal,setPointID(cell),PETSC_NULL_INTEGER,displacementDof,ierr))
            !          Do iDof = 1,numDofDamage
            !             residualDof(iDof) = residualDof(iDof) + elemScal(cell)%Gauss_C(iGauss) * &
            !                                  CrackPressureCell * (displacementCell  .DotP. elemScal(cell)%Grad_BF(iDof,iGauss))
            !          End Do ! iDof numDofDamage
            !       End Do ! iGauss
            !       PetscCall(DMPlexVecSetClosure(dmDamage,PETSC_NULL_SECTION,locResidual,setPointID(cell),residualDof,ADD_VALUES,ierr))
            !    End Do ! cell
            ! End If ! crack Pressure
            deallocate (residualDof)

            PetscCall(MEF90ElementDestroy(elemVect, ierr))
            PetscCall(MEF90ElementDestroy(elemScal, ierr))
            PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
            PetscCall(ISDestroy(setPointIS, ierr))
         end if ! setPointIS
      end do ! set
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(VecRestoreArray(MEF90DefMechCtx%plasticStrain, plasticStrainArray, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end if ! setIS
   PetscCall(DMLocalToGlobalBegin(dmDamage, locResidual, ADD_VALUES, residual, ierr))
   PetscCall(DMLocalToGlobalEnd(dmDamage, locResidual, ADD_VALUES, residual, ierr))
   PetscCall(DMRestoreLocalVector(dmDamage, locResidual, ierr))
end subroutine MEF90DefMechOperatorDamage

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechTAOGradientDamage"
!!!
!!!
!!!  MEF90DefMechTAOGradientDamage:
!!!
!!!  (c) 2012-20 2022 Alexis Marboeuf marboeua@mcmaster.ca
!!!      2023 Blaise Bourdin  bourdin@mcmaster.ca
!!!

subroutine MEF90DefMechTAOGradientDamage(taoDamage, damage, residual, MEF90DefMechCtx, ierr)
   type(tTao), intent(IN)                              :: taoDamage
   type(tVec), intent(IN)                              :: damage
   type(tVec), intent(INOUT)                           :: residual
   type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
   PetscErrorCode, intent(OUT)                         :: ierr

   type(tSNES)                                         :: dummySNES
   type(tDM)                                           :: dmDamage

   PetscCall(VecGetDM(MEF90DefMechCtx%damageLocal, dmDamage, ierr))
   PetscCall(DMGlobalToLocal(dmDamage, damage, INSERT_VALUES, MEF90DefMechCtx%damageLocal, ierr))
   PetscCall(MEF90DefMechOperatorDamage(dummySNES, damage, residual, MEF90DefMechCtx, ierr))
end subroutine MEF90DefMechTAOGradientDamage

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

subroutine MEF90DefMechBilinearFormDamage(snesDamage, damage, A, M, MEF90DefMechCtx, ierr)
   type(tSNES), intent(IN)                             :: snesDamage
   type(tVec), intent(IN)                              :: damage
   type(tMat), intent(INOUT)                           :: A, M
   type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
   PetscErrorCode, intent(OUT)                         :: ierr

   type(tDM)                                           :: dmDisplacement, dmDamage, dmTemperature, dmCohesiveDisplacement, dmPlasticStrain
   type(tPetscSection)                                 :: sectionPlasticStrain
   PetscReal, dimension(:), pointer                    :: displacementDof, damageDof, temperatureDof, plasticStrainArray
   type(MEF90_MATS)                                    :: plasticStrainCell
   type(tIS)                                           :: setIS, setPointIS
   PetscInt, dimension(:), pointer                     :: setID, setPointID
   PetscInt                                            :: set, QuadratureOrder, cell, dim, vecOffset
   type(MEF90_MATPROP), pointer                        :: matpropSet
   type(MEF90DefMechCellSetOptions_Type), pointer      :: cellSetOptions
   type(MEF90_ELEMENT_ELAST), dimension(:), pointer    :: elemVect
   type(MEF90_ELEMENT_SCAL), dimension(:), pointer     :: elemScal
   type(eDMPolytopeType)                               :: cellGeometryVect, cellGeometryScal
   type(MEF90ElementType)                              :: elemVectType, elemScalType
   PetscReal, dimension(:), pointer                    :: matDof

   type(MEF90CtxGlobalOptions_Type), pointer           :: MEF90CtxGlobalOptions
   type(MEF90DefMechGlobalOptions_Type), pointer       :: MEF90DefMechGlobalOptions
   class(MEF90DefMechAT_Type), allocatable             :: ATModel
   class(MEF90_DEFMECHSPLIT), allocatable              :: Split
   type(MEF90_MATS)                                    :: totalStrainGauss, C2
   type(MEF90_VECT)                                    :: gradDamageGauss
   PetscReal                                           :: damageGauss, temperatureGauss, EEDGaussMinus, EEDGaussPlus, C1, C3
   PetscInt                                            :: iDof, jDof, iGauss, numDofDisplacement, numDofDamage, numGauss
   character(len=MEF90MXSTRLEN)                        :: prefix

   PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag, MEF90CtxGlobalOptions, ierr))
   PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag, MEF90DefMechGlobalOptions, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%displacementLocal, dmDisplacement, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%cohesiveDisplacement, dmCohesiveDisplacement, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%damageLocal, dmDamage, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%TemperatureLocal, dmTemperature, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%plasticStrain, dmPlasticStrain, ierr))
   PetscCall(DMGetLocalSection(dmPlasticStrain, sectionPlasticStrain, ierr))

   PetscCall(DMGetDimension(dmDisplacement, dim, ierr))

   PetscCall(DMGlobalToLocal(dmDamage, damage, INSERT_VALUES, MEF90DefMechCtx%damageLocal, ierr))

   PetscCall(MatZeroEntries(A, ierr))

   !!! get IS for cell sets
   PetscCall(DMGetLabelIdIS(dmDamage, MEF90CellSetLabelName, setIS, ierr))
   PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm, setIS, ierr))

   if (.not. PetscObjectIsNull(setIS)) then
      PetscCall(ISGetIndices(setIS, setID, ierr))
      PetscCall(VecGetArray(MEF90DefMechCtx%plasticStrain, plasticStrainArray, ierr))
      do set = 1, size(setID)
         PetscCall(DMGetStratumIS(dmDamage, MEF90CellSetLabelName, setID(set), setPointIS, ierr))
         if (.not. PetscObjectIsNull(setPointIS)) then
            PetscCall(PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set), matpropSet, ierr))
            PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set), cellSetOptions, ierr))

            PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
            PetscCall(DMPlexGetCellType(dmDisplacement, setPointID(1), cellGeometryVect, ierr))
            PetscCall(DMPlexGetCellType(dmDamage, setPointID(1), cellGeometryScal, ierr))
            PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellGeometryVect, elemVectType, ierr))
            PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellGeometryScal, elemScalType, ierr))

            !!! get the ATModel and split objects
            write(prefix,'("cs",I4.4,"_")') setID(set)
            PetscCall(MEF90DefMechGetATModel(MEF90DefMechCtx%MEF90Ctx%comm, prefix, ATModel, ierr))
            PetscCall(ATModel%setFromOptions(ierr))

            if (cellSetOptions%unilateralContactHybrid) then
               Split = MEF90_DEFMECHSPLITNONE()
            else
               PetscCall(MEF90DefMechGetSplit(cellSetOptions, Split, ierr))
            end if
               !!! Allocate elements
            QuadratureOrder = max(2 * elemVectType%order, split%damageOrder + split%strainOrder)
            PetscCall(MEF90ElementCreate(dmDisplacement, setPointIS, elemVect, QuadratureOrder, elemVectType, ierr))
            PetscCall(MEF90ElementCreate(dmDamage, setPointIS, elemScal, QuadratureOrder, elemScalType, ierr))

            numDofDisplacement = size(elemVect(1)%BF(:, 1))
            numDofDamage = size(elemScal(1)%BF(:, 1))
            numGauss = size(elemScal(1)%Gauss_C)

            allocate (matDof(numDofDamage * numDofDamage))

            do cell = 1, size(setPointID)
               matDof = 0.0_kr
               PetscCall(DMPlexVecGetClosure(dmDamage, PETSC_NULL_SECTION, MEF90DefMechCtx%damageLocal, setPointID(cell), PETSC_NULL_INTEGER, damageDof, ierr))
               PetscCall(DMPlexVecGetClosure(dmDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%displacementLocal, setPointID(cell), PETSC_NULL_INTEGER, displacementDof, ierr))
               PetscCall(DMPlexVecGetClosure(dmTemperature, PETSC_NULL_SECTION, MEF90DefMechCtx%TemperatureLocal, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
               do iGauss = 1, numGauss
                  damageGauss = 0.0_kr
                  gradDamageGauss = 0.0_kr
                  if (.not. ATModel%isElastic) then
                     do iDof = 1, numDofDamage
                        damageGauss = damageGauss + damageDof(iDof) * elemScal(cell)%BF(iDof, iGauss)
                        gradDamageGauss = gradDamageGauss + damageDof(iDof) * elemScal(cell)%Grad_BF(iDof, iGauss)
                     end do ! iDof numDofDamage
                  end if

                  totalStrainGauss = 0.0_kr
                  do iDof = 1, numDofDisplacement
                     totalStrainGauss = totalStrainGauss + displacementDof(iDof) * elemVect(cell)%GradS_BF(iDof, iGauss)
                  end do ! iDof numDofDisplacement

                  temperatureGauss = 0.0_kr
                  do iDof = 1, numDofDamage
                     temperatureGauss = temperatureGauss + temperatureDof(iDof) * elemScal(cell)%BF(iDof, iGauss)
                  end do ! iDof numDofDamage
                  totalStrainGauss = totalStrainGauss - (temperatureGauss * matpropSet%LinearThermalExpansion)

! #if MEF90_DIM == 2
!!! We need something along these lines
!!! Adding terms in planestrain for plasticity with tr(p) = 0
! If (.NOT. matProp%HookesLaw%isPlaneStress) Then
!    stressGauss = stressGauss +  stiffness * ( matProp%HookesLaw%lambda*trace(plasticStrainCell)*MEF90MatS2DIdentity )
! End If
! #endif

                  PetscCall(PetscSectionGetOffset(sectionPlasticStrain, setPointID(cell), vecOffset, ierr))
                  plasticStrainCell = plasticStrainArray(vecOffset + 1:vecOffset + 1 + SIZEOFMEF90_MATS)

                  if (.not. ATModel%isElastic) then
                     call Split%EED(totalStrainGauss - plasticStrainCell, matpropSet%HookesLaw, EEDGaussPlus, EEDGaussMinus)
                  else
                     EEDGaussPlus = 0.0_kr
                  end if
                  C1 = matpropSet%fractureToughness / ATModel%cw * 0.25_kr / ATModel%internalLength
                  C2 = matpropSet%fractureToughness / ATModel%cw * 0.5_kr * ATModel%internalLength * matpropSet%toughnessAnisotropyMatrix
                  C3 = ATModel%D2a(damageGauss) * EEDGaussPlus + C1 * ATModel%D2w(damageGauss)
                  do jDof = 0, numDofDamage - 1
                     do iDof = 1, numDofDamage
                        matDof(jDof * numDofDamage + iDof) = matDof(jDof * numDofDamage + iDof) + elemScal(cell)%Gauss_C(iGauss) * ( &
                                                             C3 * elemScal(cell)%BF(iDof, iGauss) * elemScal(cell)%BF(jDof + 1, iGauss) + (C2 * elemScal(cell)%Grad_BF(iDof, iGauss) .DotP. elemScal(cell)%Grad_BF(jDof + 1, iGauss)))
                     end do ! iDof numDofDamage
                  end do ! jDof numDofDamage
               end do ! iGauss
               PetscCall(DMPlexVecRestoreClosure(dmDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%displacementLocal, setPointID(cell), PETSC_NULL_INTEGER, displacementDof, ierr))
               PetscCall(DMPlexVecRestoreClosure(dmDamage, PETSC_NULL_SECTION, MEF90DefMechCtx%damageLocal, setPointID(cell), PETSC_NULL_INTEGER, damageDof, ierr))
               PetscCall(DMPlexVecRestoreClosure(dmTemperature, PETSC_NULL_SECTION, MEF90DefMechCtx%TemperatureLocal, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
               PetscCall(DMPlexMatSetClosure(dmDamage, PETSC_NULL_SECTION, PETSC_NULL_SECTION, A, setPointID(cell), matDof, ADD_VALUES, ierr))
            end do ! cell

            deallocate (matDof)

            PetscCall(MEF90ElementDestroy(elemVect, ierr))
            PetscCall(MEF90ElementDestroy(elemScal, ierr))
            PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
            PetscCall(ISDestroy(setPointIS, ierr))
         end if ! setPointIS
      end do ! set
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(VecRestoreArray(MEF90DefMechCtx%plasticStrain, plasticStrainArray, ierr))
   end if ! setIS
   PetscCall(ISDestroy(setIS, ierr))
   PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr))
   PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr))
   PetscCall(MatCopy(A, M, SAME_NONZERO_PATTERN, ierr))
end subroutine MEF90DefMechBilinearFormDamage

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechTAOHessianDamage"
!!!
!!!
!!!  MEF90DefMechTAOHessianDamage:
!!!
!!!  (c) 2022 Alexis Marboeuf marboeua@mcmaster.ca
!!!      2023 Blaise Bourdin  bourdin@mcmaster.ca
!!!

subroutine MEF90DefMechTAOHessianDamage(taoDamage, damage, A, M, MEF90DefMechCtx, ierr)
   type(tTao), intent(IN)                              :: taoDamage
   type(tVec), intent(IN)                              :: damage
   type(tMat), intent(INOUT)                           :: A, M
   type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
   PetscErrorCode, intent(OUT)                         :: ierr

   type(tSNES)                                         :: dummySNES
   type(tDM)                                           :: dmDamage

   PetscCall(VecGetDM(MEF90DefMechCtx%damageLocal, dmDamage, ierr))
   PetscCall(DMGlobalToLocal(dmDamage, damage, INSERT_VALUES, MEF90DefMechCtx%damageLocal, ierr))
   PetscCall(MEF90DefMechBilinearFormDamage(dummySNES, damage, A, M, MEF90DefMechCtx, ierr))
end subroutine MEF90DefMechTAOHessianDamage

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

subroutine MEF90DefMechSurfaceEnergy(MEF90DefMechCtx, energy, ierr)
   type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
   PetscReal, dimension(:), pointer                    :: energy
   PetscErrorCode, intent(INOUT)                       :: ierr

   type(tDM)                                           :: dmDamage
   PetscReal, dimension(:), pointer                    :: damageDof
   type(tIS)                                           :: setIS, setPointIS
   PetscInt, dimension(:), pointer                     :: setID, setPointID
   PetscInt                                            :: set, QuadratureOrder, cell, dim
   type(MEF90_MATPROP), pointer                        :: matpropSet
   type(MEF90DefMechCellSetOptions_Type), pointer      :: cellSetOptions
   type(MEF90_ELEMENT_SCAL), dimension(:), pointer     :: elemScal
   type(eDMPolytopeType)                               :: cellGeometryScal
   type(MEF90ElementType)                              :: elemScalType

   type(MEF90CtxGlobalOptions_Type), pointer           :: MEF90CtxGlobalOptions
   type(MEF90DefMechGlobalOptions_Type), pointer       :: MEF90DefMechGlobalOptions
   class(MEF90DefMechAT_Type), allocatable             :: ATModel
   type(MEF90_MATS)                                    :: C2
   type(MEF90_VECT)                                    :: gradDamageGauss
   PetscReal                                           :: damageGauss, C1, myEnergy
   PetscInt                                            :: iDof, iGauss, numDofDamage, numGauss
   character(len=MEF90MXSTRLEN)                        :: prefix

   PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag, MEF90CtxGlobalOptions, ierr))
   PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag, MEF90DefMechGlobalOptions, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%damageLocal, dmDamage, ierr))

   PetscCall(DMGetDimension(dmDamage, dim, ierr))

      !!! get IS for cell sets
   PetscCall(DMGetLabelIdIS(dmDamage, MEF90CellSetLabelName, setIS, ierr))
   PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm, setIS, ierr))

   energy = 0.0_kr
   if (.not. PetscObjectIsNull(setIS)) then
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         myEnergy = 0.0_kr
         PetscCall(DMGetStratumIS(dmDamage, MEF90CellSetLabelName, setID(set), setPointIS, ierr))
         if (.not. PetscObjectIsNull(setPointIS)) then
            PetscCall(PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set), matpropSet, ierr))
            PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set), cellSetOptions, ierr))

            PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
            PetscCall(DMPlexGetCellType(dmDamage, setPointID(1), cellGeometryScal, ierr))
            PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellGeometryScal, elemScalType, ierr))

            !!! get the ATModel object
            write(prefix,'("cs",I4.4,"_")') setID(set)
            PetscCall(MEF90DefMechGetATModel(MEF90DefMechCtx%MEF90Ctx%comm, prefix, ATModel, ierr))
            PetscCall(ATModel%setFromOptions(ierr))
            !!! Allocate elements
            QuadratureOrder = max(ATmodel%wOrder, 2 * (elemScalType%order - 1))
            PetscCall(MEF90ElementCreate(dmDamage, setPointIS, elemScal, QuadratureOrder, elemScalType, ierr))

            numDofDamage = size(elemScal(1)%BF(:, 1))
            numGauss = size(elemScal(1)%Gauss_C)

            C1 = matpropSet%fractureToughness / ATModel%cw * 0.25_kr / ATModel%internalLength
            C2 = matpropSet%fractureToughness / ATModel%cw * 0.25_kr * ATModel%internalLength * matpropSet%toughnessAnisotropyMatrix
            do cell = 1, size(setPointID)
               PetscCall(DMPlexVecGetClosure(dmDamage, PETSC_NULL_SECTION, MEF90DefMechCtx%damageLocal, setPointID(cell), PETSC_NULL_INTEGER, damageDof, ierr))
               do iGauss = 1, numGauss
                  damageGauss = 0.0_kr
                  gradDamageGauss = 0.0_kr
                  if (.not. ATModel%isElastic) then
                     do iDof = 1, numDofDamage
                        damageGauss = damageGauss + damageDof(iDof) * elemScal(cell)%BF(iDof, iGauss)
                        gradDamageGauss = gradDamageGauss + damageDof(iDof) * elemScal(cell)%Grad_BF(iDof, iGauss)
                     end do ! iDof numDofDamage
                  end if
                  myEnergy = myEnergy + elemScal(cell)%Gauss_C(iGauss) * &
                             (C1 * ATModel%w(damageGauss) + (C2 * gradDamageGauss .DotP. gradDamageGauss))
               end do ! iGauss
               PetscCall(DMPlexVecRestoreClosure(dmDamage, PETSC_NULL_SECTION, MEF90DefMechCtx%damageLocal, setPointID(cell), PETSC_NULL_INTEGER, damageDof, ierr))
            end do ! cell
            PetscCall(MEF90ElementDestroy(elemScal, ierr))
         end if ! setPointIS
         PetscCall(ISDestroy(setPointIS, ierr))
         PetscCallMPI(MPI_AllReduce(myEnergy, energy(set), 1, MPIU_SCALAR, MPI_SUM, MEF90DefMechCtx%MEF90Ctx%comm, ierr))
      end do ! set
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
   end if ! setIS
   PetscCall(ISDestroy(setIS, ierr))
end subroutine MEF90DefMechSurfaceEnergy

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechTAObjectiveDamage"
!!!
!!!
!!!  MEF90DefMechTAOObjectiveDamage:
!!!
!!!  (c) 2022 Alexis Marboeuf marboeua@mcmaster.ca
!!!      2023 Blaise Bourdin  bourdin@mcmaster.ca
!!!

subroutine MEF90DefMechTAOObjectiveDamage(taoDamage, damage, energy, MEF90DefMechCtx, ierr)
   type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
   type(tTao), intent(IN)                              :: taoDamage
   type(tVec), intent(IN)                              :: damage
   PetscReal, intent(INOUT)                            :: energy
   PetscErrorCode, intent(INOUT)                       :: ierr

   PetscReal, dimension(:), pointer                    :: surfaceEnergy, elasticEnergy
   type(tDM)                                           :: dmDamage
   type(tIS)                                           :: setIS
   PetscInt                                            :: numSet

   PetscCall(VecGetDM(MEF90DefMechCtx%damageLocal, dmDamage, ierr))
   PetscCall(DMGlobalToLocal(dmDamage, damage, INSERT_VALUES, MEF90DefMechCtx%damageLocal, ierr))

      !!! get IS for cell sets
   PetscCall(DMGetLabelIdIS(dmDamage, MEF90CellSetLabelName, setIS, ierr))
   PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm, setIS, ierr))

   numSet = 0
   if (.not. PetscObjectIsNull(setIS)) then
      PetscCall(ISGetSize(setIS, numSet, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end if

   allocate (surfaceEnergy(numSet))
   allocate (elasticEnergy(numSet))
   surfaceEnergy = 0.0_kr
   elasticEnergy = 0.0_kr
   PetscCall(MEF90DefMechSurfaceEnergy(MEF90DefMechCtx, surfaceEnergy, ierr))
   PetscCall(MEF90DefMechElasticEnergy(MEF90DefMechCtx, elasticEnergy, ierr))

   energy = sum(surfaceEnergy) + sum(elasticEnergy)
end subroutine MEF90DefMechTAOObjectiveDamage

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

subroutine MEF90DefMechCrackVolume(MEF90DefMechCtx, CrackVolume, ierr)
   type(MEF90DefMechCtx_Type), intent(IN)              :: MEF90DefMechCtx
   PetscReal, dimension(:), pointer                    :: CrackVolume
   PetscErrorCode, intent(INOUT)                       :: ierr

   type(tDM)                                           :: dmDamage, dmDisplacement
   PetscReal, dimension(:), pointer                    :: damageDof, displacementDof
   type(tIS)                                           :: setIS, setPointIS
   PetscInt, dimension(:), pointer                     :: setID, setPointID
   PetscInt                                            :: set, QuadratureOrder, cell, dim
   type(MEF90_MATPROP), pointer                        :: matpropSet
   type(MEF90DefMechCellSetOptions_Type), pointer      :: cellSetOptions
   type(MEF90_ELEMENT_ELAST), dimension(:), pointer    :: elemVect
   type(MEF90_ELEMENT_SCAL), dimension(:), pointer     :: elemScal
   type(eDMPolytopeType)                               :: cellGeometryVect
   type(eDMPolytopeType)                               :: cellGeometryScal
   type(MEF90ElementType)                              :: elemVectType, elemScalType

   type(MEF90CtxGlobalOptions_Type), pointer           :: MEF90CtxGlobalOptions
   type(MEF90DefMechGlobalOptions_Type), pointer       :: MEF90DefMechGlobalOptions
   class(MEF90DefMechAT_Type), allocatable             :: ATModel
   type(MEF90_VECT)                                    :: gradDamageGauss, displacementCell
   PetscReal                                           :: myCrackVolume
   PetscInt                                            :: iDof, iGauss, numDofDamage, numDofDisplacement, numGauss
   character(len=MEF90MXSTRLEN)                        :: prefix

   PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag, MEF90CtxGlobalOptions, ierr))
   PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag, MEF90DefMechGlobalOptions, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%damageLocal, dmDamage, ierr))
   PetscCall(VecGetDM(MEF90DefMechCtx%displacementLocal, dmDisplacement, ierr))

   PetscCall(DMGetDimension(dmDamage, dim, ierr))

      !!! get IS for cell sets
   PetscCall(DMGetLabelIdIS(dmDamage, MEF90CellSetLabelName, setIS, ierr))
   PetscCall(MEF90ISAllGatherMerge(MEF90DefMechCtx%MEF90Ctx%comm, setIS, ierr))

   CrackVolume = 0.0_kr
   if (.not. PetscObjectIsNull(setIS)) then
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         myCrackVolume = 0.0_kr
         PetscCall(DMGetStratumIS(dmDamage, MEF90CellSetLabelName, setID(set), setPointIS, ierr))
         if (.not. PetscObjectIsNull(setPointIS)) then
            PetscCall(PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set), matpropSet, ierr))
            PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set), cellSetOptions, ierr))

            PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
            PetscCall(DMPlexGetCellType(dmDisplacement, setPointID(1), cellGeometryVect, ierr))
            PetscCall(DMPlexGetCellType(dmDamage, setPointID(1), cellGeometryScal, ierr))
            PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellGeometryVect, elemVectType, ierr))
            PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellGeometryScal, elemScalType, ierr))

            !!! get the ATModel and split objects
            write(prefix,'("cs",I4.4,"_")') setID(set)
            PetscCall(MEF90DefMechGetATModel(MEF90DefMechCtx%MEF90Ctx%comm, prefix, ATModel, ierr))
            PetscCall(ATModel%setFromOptions(ierr))

            !!! Allocate elements
            QuadratureOrder = 2 * (elemScalType%order - 1)
            PetscCall(MEF90ElementCreate(dmDisplacement, setPointIS, elemVect, QuadratureOrder, elemVectType, ierr))
            PetscCall(MEF90ElementCreate(dmDamage, setPointIS, elemScal, QuadratureOrder, elemScalType, ierr))

            numDofDisplacement = size(elemVect(1)%BF(:, 1))
            numDofDamage = size(elemScal(1)%BF(:, 1))
            numGauss = size(elemScal(1)%Gauss_C)

            do cell = 1, size(setPointID)
               PetscCall(DMPlexVecGetClosure(dmDamage, PETSC_NULL_SECTION, MEF90DefMechCtx%damageLocal, setPointID(cell), PETSC_NULL_INTEGER, damageDof, ierr))
               PetscCall(DMPlexVecGetClosure(dmDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%displacementLocal, setPointID(cell), PETSC_NULL_INTEGER, displacementDof, ierr))
               do iGauss = 1, numGauss
                  gradDamageGauss = 0.0_kr
                  if (.not. ATModel%isElastic) then
                     do iDof = 1, numDofDamage
                        gradDamageGauss = gradDamageGauss + damageDof(iDof) * elemScal(cell)%Grad_BF(iDof, iGauss)
                     end do ! iDof numDofDamage
                  end if

                  displacementCell = 0.0_kr
                  do iDof = 1, numDofDisplacement
                     displacementCell = displacementCell + displacementDof(iDof) * elemVect(cell)%BF(iDof, iGauss)
                  end do ! iDof numDofDisplacement

                  do iDof = 1, numDofDamage
                     myCrackVolume = myCrackVolume + elemScal(cell)%Gauss_C(iGauss) * ( &
                                     gradDamageGauss .DotP. displacementCell)
                  end do ! iDof numDofDamage
               end do ! iGauss
               PetscCall(DMPlexVecRestoreClosure(dmDamage, PETSC_NULL_SECTION, MEF90DefMechCtx%damageLocal, setPointID(cell), PETSC_NULL_INTEGER, damageDof, ierr))
               PetscCall(DMPlexVecRestoreClosure(dmDisplacement, PETSC_NULL_SECTION, MEF90DefMechCtx%displacementLocal, setPointID(cell), PETSC_NULL_INTEGER, displacementDof, ierr))
            end do ! cell

            PetscCall(MEF90ElementDestroy(elemScal, ierr))
            PetscCall(MEF90ElementDestroy(elemVect, ierr))
            PetscCall(ISDestroy(setPointIS, ierr))
         end if ! setPointIS
         PetscCallMPI(MPI_AllReduce(myCrackVolume, CrackVolume(set), 1, MPIU_SCALAR, MPI_SUM, MEF90DefMechCtx%MEF90Ctx%comm, ierr))
      end do ! set
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end if ! setIS
end subroutine MEF90DefMechCrackVolume
end module MEF90_APPEND(m_MEF90_DefMechAssembly,MEF90_DIM)D
