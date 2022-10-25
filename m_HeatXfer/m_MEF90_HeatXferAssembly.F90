#include "../MEF90/mef90.inc"
Module MEF90_APPEND(m_MEF90_HeatXferAssembly,MEF90_DIM)D
#include "petsc/finclude/petsc.h"
   Use m_MEF90
   Use m_MEF90_HeatXferCtx
   Implicit none
   
   Private
   Public MEF90HeatXferOperator, &
          MEF90HeatXferBilinearForm, &
          MEF90HeatXFerEnergy, &
          MEF90HeatXFerIFunction, &
          MEF90HeatXferIJacobian

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferOperator"
!!!
!!!  
!!!  MEF90HeatXferOperator: Build the operator. When called in SNES, the solution time should always match the target time, 
!!!                         so there is no need for interpolation of the fluxes, external, and boundary values
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!      2022    Alexis Marboeuf marboeua@mcmaster.ca
!!!

   Subroutine MEF90HeatXferOperator(snesTemp,x,residual,MEF90HeatXferCtx,ierr)
      Type(tSNES),Intent(IN)                          :: snesTemp
      Type(tVec),Intent(IN)                           :: x
      Type(tVec),Intent(INOUT)                        :: residual
      Type(MEF90HeatXferCtx_Type),Intent(IN)          :: MEF90HeatXferCtx
      PetscErrorCode,Intent(INOUT)                    :: ierr
   
      Type(tDM)                                       :: dmTemperature
      Type(tVec)                                      :: locTemperature,locResidual
      Type(tIS)                                       :: setIS,setPointIS
      PetscInt,dimension(:),Pointer                   :: setID,setPointID
      PetscInt                                        :: set,QuadratureOrder
      PetscInt                                        :: cell, iDof,jDof,iGauss
      Type(MEF90HeatXferCellSetOptions_Type),pointer  :: cellSetOptions
      Type(MEF90HeatXferFaceSetOptions_Type),pointer  :: faceSetOptions
      Type(MEF90_MATPROP),pointer                     :: matpropSet
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer   :: elem
      Type(MEF90ElementType)                          :: elementType
      DMPolytopeType                                  :: cellType
      Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90CtxGlobalOptions
      Type(MEF90HeatXferGlobalOptions_Type),pointer   :: MEF90HeatXferGlobalOptions
      PetscReal,Dimension(:),Pointer                  :: temperatureDof,fluxDof,boundaryFluxDof,externalTemperatureDof,residualDof
      Type(MEF90_VECT)                                :: advectionVec
      
      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90HeatXferCtx%MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr))
      PetscCall(PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr))
      PetscCall(VecGetDM(MEF90HeatXferCtx%temperatureLocal,dmTemperature,ierr))

      PetscCall(DMGetLocalVector(dmTemperature,locTemperature,ierr))
      PetscCall(DMGetLocalVector(dmTemperature,locResidual,ierr))
      PetscCall(MEF90VecGlobalToLocalConstraint(x,MEF90HeatXferCtx%temperatureLocal,locTemperature,ierr))

      PetscCall(VecSet(residual,0.0_Kr,ierr))
      PetscCall(VecSet(locResidual,0.0_Kr,ierr))
         
      !!! cell-based contributions
      PetscCall(DMGetLabelIdIS(dmTemperature,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm,setIS,ierr))
      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         Do set = 1,size(setID)
            PetscCall(DMGetStratumIS(dmTemperature,MEF90CellSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
               PetscCall(PetscBagGetDataMEF90MatProp(MEF90HeatXferCtx%MaterialPropertiesBag(set),matpropSet,ierr))
               PetscCall(PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr))
         
               PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(DMPlexGetCellType(dmTemperature,setPointID(1),cellType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellType,elementType,ierr))
               QuadratureOrder = elementType%order * 2
               PetscCall(MEF90ElementCreate(dmTemperature,setPointIS,elem,QuadratureOrder,elementType,ierr))

               Allocate(residualDof(size(elem(1)%BF(:,1))))

               Do cell = 1, size(setPointID)
                  residualDof = 0.0_Kr
                  !!! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                  !!! If this happens, we will need to protect this loop
                  PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,locTemperature,setPointID(cell),temperatureDof,ierr))
                  Do iGauss = 1,size(elem(cell)%Gauss_C)
                     Do jDof = 1, size(elem(cell)%BF(:,1))
                        Do iDof = 1, size(elem(cell)%BF(:,1))
                           residualDof(jDof) = residualDof(jDof) + (matpropSet%ThermalConductivity*temperatureDof(iDof)*elem(cell)%Grad_BF(iDof,iGauss) .dotP. elem(cell)%Grad_BF(jDof,iGauss))*elem(cell)%Gauss_C(iGauss)
                        End Do ! iDof
                     End Do ! jDof
                  End Do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmTemperature,PETSC_NULL_SECTION,locTemperature,setPointID(cell),temperatureDof,ierr))
                  PetscCall(DMPlexVecSetClosure(dmTemperature,PETSC_NULL_SECTION,locResidual,setPointID(cell),residualDof,ADD_VALUES,ierr))
               End Do ! cell

               If (norm2(cellSetOptions%advectionVector) /= 0.0_Kr) Then
                  advectionVec = cellSetOptions%advectionVector
                  Do cell = 1, size(setPointID)
                     residualDof = 0.0_Kr
                     !!! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !!! If this happens, we will need to protect this loop
                     PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,locTemperature,setPointID(cell),temperatureDof,ierr))
                     Do iGauss = 1,size(elem(cell)%Gauss_C)
                        Do jDof = 1, size(elem(cell)%BF(:,1))
                           Do iDof = 1, size(elem(cell)%BF(:,1))
                              residualDof(jDof) = residualDof(jDof) - (matpropSet%Density*matpropSet%SpecificHeat*advectionVec .dotP. temperatureDof(iDof)*elem(cell)%Grad_BF(iDof,iGauss))*elem(cell)%BF(jDof,iGauss)*elem(cell)%Gauss_C(iGauss)
                           End Do ! iDof
                        End Do ! jDof
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmTemperature,PETSC_NULL_SECTION,locTemperature,setPointID(cell),temperatureDof,ierr))
                     PetscCall(DMPlexVecSetClosure(dmTemperature,PETSC_NULL_SECTION,locResidual,setPointID(cell),residualDof,ADD_VALUES,ierr))
                  End Do ! cell
               End If ! cellSetOptions%advectionVector

               If (cellSetOptions%flux /= 0.0_Kr) Then
                  Do cell = 1, size(setPointID)
                     residualDof = 0.0_Kr
                     !!! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !!! If this happens, we will need to protect this loop
                     PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,MEF90HeatXferCtx%fluxLocal,setPointID(cell),fluxDof,ierr))
                     Do iGauss = 1,size(elem(cell)%Gauss_C)
                        Do iDof = 1, size(elem(cell)%BF(:,1))
                           residualDof(iDof) = residualDof(iDof) - fluxDof(iDof)*elem(cell)%BF(iDof,iGauss)*elem(cell)%Gauss_C(iGauss)
                        End Do ! iDof
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmTemperature,PETSC_NULL_SECTION,MEF90HeatXferCtx%fluxLocal,setPointID(cell),fluxDof,ierr))
                     PetscCall(DMPlexVecSetClosure(dmTemperature,PETSC_NULL_SECTION,locResidual,setPointID(cell),residualDof,ADD_VALUES,ierr))
                  End Do ! cell
               End If ! cellSetOptions%flux

               DeAllocate(residualDof)
               PetscCall(MEF90ElementDestroy(elem,ierr))
               PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
            End If ! pointIS
            PetscCall(ISDestroy(setPointIS,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))

      !!! face-based contributions
      PetscCall(DMGetLabelIdIS(dmTemperature,MEF90FaceSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm,setIS,ierr))
      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         Do set = 1,size(setID)
            PetscCall(DMGetStratumIS(dmTemperature,MEF90FaceSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
              PetscCall(PetscBagGetDataMEF90HeatXferCtxFaceSetOptions(MEF90HeatXferCtx%FaceSetOptionsBag(set),faceSetOptions,ierr))

               If (faceSetOptions%boundaryFlux /= 0.0_Kr) Then
                  PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
                  PetscCall(DMPlexGetCellType(dmTemperature,setPointID(1),cellType,ierr))
                  PetscCall(MEF90ElementGetTypeBoundary(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellType,elementType,ierr))
                  QuadratureOrder = elementType%order * 2
                  PetscCall(MEF90ElementCreate(dmTemperature,setPointIS,elem,QuadratureOrder,elementType,ierr))
                  Allocate(residualDof(size(elem(1)%BF(:,1))))
                  Do cell = 1, size(setPointID)
                     residualDof = 0.0_Kr
                     !!! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !!! If this happens, we will need to protect this loop
                     PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,MEF90HeatXFerCtx%boundaryFluxLocal,setPointID(cell),boundaryFluxDof,ierr))
                     Do iGauss = 1,size(elem(cell)%Gauss_C)
                        Do iDof = 1, size(elem(cell)%BF(:,1))
                           residualDof(iDof) = residualDof(iDof) - boundaryFluxDof(iDof)*elem(cell)%BF(iDof,iGauss)*elem(cell)%Gauss_C(iGauss)
                        End Do ! iDof
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmTemperature,PETSC_NULL_SECTION,MEF90HeatXFerCtx%boundaryFluxLocal,setPointID(cell),boundaryFluxDof,ierr))
                     PetscCall(DMPlexVecSetClosure(dmTemperature,PETSC_NULL_SECTION,locResidual,setPointID(cell),residualDof,ADD_VALUES,ierr))
                  End Do ! cell
                  DeAllocate(residualDof)
                  PetscCall(MEF90ElementDestroy(elem,ierr))
                  PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
               End If ! faceSetOptions%boundaryFlux

               If (faceSetOptions%surfaceThermalConductivity /= 0.0_Kr) Then
                  PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
                  PetscCall(DMPlexGetCellType(dmTemperature,setPointID(1),cellType,ierr))
                  PetscCall(MEF90ElementGetTypeBoundary(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellType,elementType,ierr))
                  QuadratureOrder = elementType%order * 2
                  PetscCall(MEF90ElementCreate(dmTemperature,setPointIS,elem,QuadratureOrder,elementType,ierr))
                  Allocate(residualDof(size(elem(1)%BF(:,1))))
                  Do cell = 1, size(setPointID)
                     residualDof = 0.0_Kr
                     !!! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !!! If this happens, we will need to protect this loop
                     PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,locTemperature,setPointID(cell),temperatureDof,ierr))
                     PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,MEF90HeatXFerCtx%externalTemperatureLocal,setPointID(cell),externalTemperatureDof,ierr))
                     Do iGauss = 1,size(elem(cell)%Gauss_C)
                        Do jDof = 1, size(elem(cell)%BF(:,1))
                           Do iDof = 1, size(elem(cell)%BF(:,1))
                              residualDof(jDof) = residualDof(jDof) + faceSetOptions%surfaceThermalConductivity*temperatureDof(iDof)*elem(cell)%BF(iDof,iGauss)*elem(cell)%BF(jDof,iGauss)*elem(cell)%Gauss_C(iGauss)
                           End Do ! iDof
                           residualDof(jDof) = residualDof(jDof) - faceSetOptions%surfaceThermalConductivity*externalTemperatureDof(jDof)*elem(cell)%BF(jDof,iGauss)*elem(cell)%Gauss_C(iGauss)
                        End Do ! jDof
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmTemperature,PETSC_NULL_SECTION,MEF90HeatXFerCtx%externalTemperatureLocal,setPointID(cell),externalTemperatureDof,ierr))
                     PetscCall(DMPlexVecRestoreClosure(dmTemperature,PETSC_NULL_SECTION,locTemperature,setPointID(cell),temperatureDof,ierr))
                     PetscCall(DMPlexVecSetClosure(dmTemperature,PETSC_NULL_SECTION,locResidual,setPointID(cell),residualDof,ADD_VALUES,ierr))
                  End Do ! cell
                  DeAllocate(residualDof)
                  PetscCall(MEF90ElementDestroy(elem,ierr))
                  PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
               End If ! faceSetOptions%surfaceThermalConductivity
            End If ! pointIS
            PetscCall(ISDestroy(setPointIS,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))
      PetscCall(DMLocalToGlobalBegin(dmTemperature,locResidual,ADD_VALUES,residual,ierr))
      PetscCall(DMLocalToGlobalEnd(dmTemperature,locResidual,ADD_VALUES,residual,ierr))
      PetscCall(DMRestoreLocalVector(dmTemperature,locTemperature,ierr))
      PetscCall(DMRestoreLocalVector(dmTemperature,locResidual,ierr))
   End Subroutine MEF90HeatXferOperator
   
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferBilinearForm"
!!!
!!!  
!!!  MEF90HeatXferBilinearForm:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!      2022    Alexis Marboeuf marboeua@mcmaster.ca
!!!

   Subroutine MEF90HeatXferBilinearForm(snesTemp,x,A,M,MEF90HeatXferCtx,ierr)
      Type(tSNES),Intent(IN)                          :: snesTemp
      Type(tVec),Intent(IN)                           :: x
      Type(tMat),Intent(INOUT)                        :: A,M
      Type(MEF90HeatXferCtx_Type),Intent(IN)          :: MEF90HeatXferCtx
      PetscErrorCode,Intent(INOUT)                    :: ierr  
         
      Type(tDM)                                       :: dmTemperature
      Type(tIS)                                       :: setIS,setPointIS
      PetscInt,dimension(:),Pointer                   :: setID,setPointID
      PetscInt                                        :: set,QuadratureOrder
      PetscInt                                        :: cell,iDof,jDof,iGauss,nbDof
      Type(MEF90HeatXferCellSetOptions_Type),pointer  :: cellSetOptions
      Type(MEF90HeatXferFaceSetOptions_Type),pointer  :: faceSetOptions
      Type(MEF90_MATPROP),pointer                     :: matpropSet
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer   :: elem
      Type(MEF90ElementType)                          :: elementType
      DMPolytopeType                                  :: cellType
      Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90CtxGlobalOptions
      Type(MEF90HeatXferGlobalOptions_Type),pointer   :: MEF90HeatXferGlobalOptions
      PetscReal,Dimension(:),Pointer                  :: matDof
      Type(MEF90_VECT)                                :: advectionVec
      
      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90HeatXferCtx%MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr))
      PetscCall(PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr))
      PetscCall(VecGetDM(MEF90HeatXferCtx%temperatureLocal,dmTemperature,ierr))

      PetscCall(MatZeroEntries(A,ierr))
         
      !!! cell-based gradient contributions
      PetscCall(DMGetLabelIdIS(dmTemperature,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm,setIS,ierr))
      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         Do set = 1,size(setID)
            PetscCall(DMGetStratumIS(dmTemperature,MEF90CellSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
               PetscCall(PetscBagGetDataMEF90MatProp(MEF90HeatXferCtx%MaterialPropertiesBag(set),matpropSet,ierr))
               PetscCall(PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr))
         
               PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(DMPlexGetCellType(dmTemperature,setPointID(1),cellType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellType,elementType,ierr))
               QuadratureOrder = elementType%order * 2
               PetscCall(MEF90ElementCreate(dmTemperature,setPointIS,elem,QuadratureOrder,elementType,ierr))

               nbDof = size(elem(1)%BF(:,1))
               Allocate(matDof(nbDof*nbDof))

               Do cell = 1, size(setPointID)
                  matDof = 0.0_Kr
                  !!! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                  !!! If this happens, we will need to protect this loop
                  Do iGauss = 1,size(elem(cell)%Gauss_C)
                     Do jDof = 0, nbDof-1
                        Do iDof = 1, nbDof
                           matDof(jDof*nbDof+iDof) = matDof(jDof*nbDof+iDof) + (matpropSet%ThermalConductivity*elem(cell)%Grad_BF(iDof,iGauss) .dotP. elem(cell)%Grad_BF(jDof+1,iGauss))*elem(cell)%Gauss_C(iGauss)
                        End Do ! iDof
                     End Do ! jDof
                  End Do ! iGauss
                  PetscCall(DMPlexMatSetClosure(dmTemperature,PETSC_NULL_SECTION,PETSC_NULL_SECTION,A,setPointID(cell),matDof,ADD_VALUES,ierr))
               End Do ! cell

               If (norm2(cellSetOptions%advectionVector) /= 0.0_Kr) Then
                  advectionVec = cellSetOptions%advectionVector
                  Do cell = 1, size(setPointID)
                     matDof = 0.0_Kr
                     !!! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !!! If this happens, we will need to protect this loop
                     Do iGauss = 1,size(elem(cell)%Gauss_C)
                        Do jDof = 0, nbDof-1
                           Do iDof = 1, nbDof
                              matDof(jDof*nbDof+iDof) = matDof(jDof*nbDof+iDof) - (matpropSet%Density*matpropSet%SpecificHeat*advectionVec .dotP. elem(cell)%Grad_BF(iDof,iGauss))*elem(cell)%BF(jDof+1,iGauss)*elem(cell)%Gauss_C(iGauss)
                           End Do ! iDof
                        End Do ! jDof
                     End Do ! iGauss
                     PetscCall(DMPlexMatSetClosure(dmTemperature,PETSC_NULL_SECTION,PETSC_NULL_SECTION,A,setPointID(cell),matDof,ADD_VALUES,ierr))
                  End Do ! cell
               End If ! cellSetOptions%advectionVector

               DeAllocate(matDof)
               PetscCall(MEF90ElementDestroy(elem,ierr))
               PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
            End If ! pointIS
            PetscCall(ISDestroy(setPointIS,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))

      !!! face-based energies
      PetscCall(DMGetLabelIdIS(dmTemperature,MEF90FaceSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm,setIS,ierr))
      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         Do set = 1,size(setID)
            PetscCall(DMGetStratumIS(dmTemperature,MEF90FaceSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
              PetscCall(PetscBagGetDataMEF90HeatXferCtxFaceSetOptions(MEF90HeatXferCtx%FaceSetOptionsBag(set),faceSetOptions,ierr))

               If (faceSetOptions%surfaceThermalConductivity /= 0.0_Kr) Then
                  PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
                  PetscCall(DMPlexGetCellType(dmTemperature,setPointID(1),cellType,ierr))
                  PetscCall(MEF90ElementGetTypeBoundary(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellType,elementType,ierr))
                  QuadratureOrder = elementType%order * 2
                  PetscCall(MEF90ElementCreate(dmTemperature,setPointIS,elem,QuadratureOrder,elementType,ierr))

                  nbDof = size(elem(1)%BF(:,1))
                  Allocate(matDof(nbDof*nbDof))

                  Do cell = 1, size(setPointID)
                     matDof = 0.0_Kr
                     !!! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !!! If this happens, we will need to protect this loop
                     Do iGauss = 1,size(elem(cell)%Gauss_C)
                        Do jDof = 0, nbDof-1
                           Do iDof = 1, nbDof
                              matDof(jDof*nbDof+iDof) = matDof(jDof*nbDof+iDof) + faceSetOptions%surfaceThermalConductivity*elem(cell)%BF(iDof,iGauss)*elem(cell)%BF(jDof+1,iGauss)*elem(cell)%Gauss_C(iGauss)
                           End Do ! iDof
                        End Do ! jDof
                     End Do ! iGauss
                     PetscCall(DMPlexMatSetClosure(dmTemperature,PETSC_NULL_SECTION,PETSC_NULL_SECTION,A,setPointID(cell),matDof,ADD_VALUES,ierr))
                  End Do ! cell

                  DeAllocate(matDof)
                  PetscCall(MEF90ElementDestroy(elem,ierr))
                  PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
               End If ! faceSetOptions%surfaceThermalConductivity
            End If ! pointIS
            PetscCall(ISDestroy(setPointIS,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))   
      PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr))
      PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr))
      PetscCall(MatCopy(A,M,SAME_NONZERO_PATTERN,ierr))
   End Subroutine MEF90HeatXferBilinearForm

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXFerEnergy"
!!!
!!!  
!!!  MEF90HeatXFerEnergy:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90HeatXFerEnergy(MEF90HeatXferCtx,energy,bodyWork,surfaceWork,ierr)
      Type(MEF90HeatXferCtx_Type),Intent(IN)          :: MEF90HeatXferCtx
      PetscReal,Dimension(:),Pointer                  :: energy,bodyWork,surfaceWork
      PetscErrorCode,Intent(INOUT)                    :: ierr
   
      Type(tDM)                                       :: dmTemperature
      Type(tIS)                                       :: setIS,setPointIS
      PetscInt,dimension(:),Pointer                   :: setID,setPointID
      PetscInt                                        :: set,QuadratureOrder
      PetscInt                                        :: cell, iDof,iGauss
      Type(MEF90HeatXferCellSetOptions_Type),pointer  :: cellSetOptions
      Type(MEF90HeatXferFaceSetOptions_Type),pointer  :: faceSetOptions
      Type(MEF90_MATPROP),pointer                     :: matpropSet
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer   :: elem
      Type(MEF90ElementType)                          :: elementType
      DMPolytopeType                                  :: cellType
      Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90CtxGlobalOptions
      Type(MEF90HeatXferGlobalOptions_Type),pointer   :: MEF90HeatXferGlobalOptions
      PetscReal,Dimension(:),Pointer                  :: TemperatureDof,fluxDof,boundaryFluxDof
      Type(MEF90_VECT)                                :: gradTemperatureCell
      PetscReal                                       :: bodyWorkCell,surfaceWorkCell
      PetscReal                                       :: myEnergy,myBodyWork,mySurfaceWork
      
      
      energy      = 0.0_Kr
      bodyWork    = 0.0_Kr
      surfaceWork = 0.0_Kr
      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90HeatXferCtx%MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr))
      PetscCall(PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr))
      PetscCall(VecGetDM(MEF90HeatXferCtx%temperatureLocal,dmTemperature,ierr))
         
      !!! cell-based energies
      PetscCall(DMGetLabelIdIS(dmTemperature,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm,setIS,ierr))
      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         Do set = 1,size(setID)
            PetscCall(DMGetStratumIS(dmTemperature,MEF90CellSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
               PetscCall(PetscBagGetDataMEF90MatProp(MEF90HeatXferCtx%MaterialPropertiesBag(set),matpropSet,ierr))
               PetscCall(PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr))
         
               PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(DMPlexGetCellType(dmTemperature,setPointID(1),cellType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellType,elementType,ierr))
               QuadratureOrder = elementType%order * 2
               PetscCall(MEF90ElementCreate(dmTemperature,setPointIS,elem,QuadratureOrder,elementType,ierr))

               myEnergy      = 0.0_Kr
               Do cell = 1, size(setPointID)
                  !!! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                  !!! If this happens, we will need to protect this loop
                  PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,MEF90HeatXferCtx%temperatureLocal,setPointID(cell),temperatureDof,ierr))
                  Do iGauss = 1,size(elem(cell)%Gauss_C)
                     gradTemperatureCell = 0.0_Kr
                     Do iDof = 1, size(elem(cell)%BF(:,1))
                        gradTemperatureCell = gradTemperatureCell + temperatureDof(iDof) * elem(cell)%Grad_BF(iDof,iGauss)
                     End Do ! iDof
                     myEnergy = myEnergy + ((matpropSet%ThermalConductivity * gradTemperatureCell) .dotP. gradTemperatureCell) * elem(cell)%Gauss_C(iGauss)
                  End Do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmTemperature,PETSC_NULL_SECTION,MEF90HeatXferCtx%temperatureLocal,setPointID(cell),temperatureDof,ierr))
               End Do ! cell

               myBodyWork    = 0.0_Kr
               If (cellSetOptions%flux /= 0.0_Kr) Then
                  Do cell = 1, size(setPointID)
                     !!! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !!! If this happens, we will need to protect this loop
                     PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,MEF90HeatXferCtx%temperatureLocal,setPointID(cell),temperatureDof,ierr))
                     PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,MEF90HeatXferCtx%fluxLocal,setPointID(cell),fluxDof,ierr))
                     Do iGauss = 1,size(elem(cell)%Gauss_C)
                        bodyWorkCell = 0.0_Kr
                        Do iDof = 1, size(elem(cell)%BF(:,1))
                           bodyWorkCell = bodyWorkCell + fluxDof(iDof) * temperatureDof(iDof) * elem(cell)%BF(iDof,iGauss)
                        End Do ! iDof
                        myBodyWork = myBodyWork + bodyWorkCell * elem(cell)%Gauss_C(iGauss)
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmTemperature,PETSC_NULL_SECTION,MEF90HeatXferCtx%temperatureLocal,setPointID(cell),temperatureDof,ierr))
                     PetscCall(DMPlexVecRestoreClosure(dmTemperature,PETSC_NULL_SECTION,MEF90HeatXferCtx%fluxLocal,setPointID(cell),fluxDof,ierr))
                  End Do ! cell
               End If ! cellSetOptions%flux

               PetscCall(MEF90ElementDestroy(elem,ierr))
               PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
            End If ! pointIS
            PetscCall(ISDestroy(setPointIS,ierr))
            PetscCallMPI(MPI_AllReduce(myEnergy,energy(set),1,MPIU_SCALAR,MPI_SUM,MEF90HeatXferCtx%MEF90Ctx%comm,ierr))
            PetscCallMPI(MPI_AllReduce(myBodyWork,bodyWork(set),1,MPIU_SCALAR,MPI_SUM,MEF90HeatXferCtx%MEF90Ctx%comm,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))

      !!! face-based energies
      PetscCall(DMGetLabelIdIS(dmTemperature,MEF90FaceSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm,setIS,ierr))
      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         Do set = 1,size(setID)
            PetscCall(DMGetStratumIS(dmTemperature,MEF90FaceSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
              PetscCall(PetscBagGetDataMEF90HeatXferCtxFaceSetOptions(MEF90HeatXferCtx%FaceSetOptionsBag(set),faceSetOptions,ierr))
               mySurfaceWork    = 0.0_Kr
               If (faceSetOptions%boundaryFlux /= 0.0_Kr) Then
                  PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
                  PetscCall(DMPlexGetCellType(dmTemperature,setPointID(1),cellType,ierr))
                  PetscCall(MEF90ElementGetTypeBoundary(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellType,elementType,ierr))
                  QuadratureOrder = elementType%order * 2
                  PetscCall(MEF90ElementCreate(dmTemperature,setPointIS,elem,QuadratureOrder,elementType,ierr))

                  Do cell = 1, size(setPointID)
                     !!! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !!! If this happens, we will need to protect this loop
                     PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,MEF90HeatXferCtx%temperatureLocal,setPointID(cell),temperatureDof,ierr))
                     PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,MEF90HeatXferCtx%boundaryFluxLocal,setPointID(cell),boundaryFluxDof,ierr))
                     Do iGauss = 1,size(elem(cell)%Gauss_C)
                        surfaceWorkCell = 0.0_Kr
                        Do iDof = 1, size(elem(cell)%BF(:,1))
                           surfaceWorkCell = surfaceWorkCell + boundaryFluxDof(iDof) * temperatureDof(iDof) * elem(cell)%BF(iDof,iGauss)
                        End Do ! iDof
                        mySurfaceWork = mySurfaceWork + surfaceWorkCell * elem(cell)%Gauss_C(iGauss)
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmTemperature,PETSC_NULL_SECTION,MEF90HeatXferCtx%temperatureLocal,setPointID(cell),temperatureDof,ierr))
                     PetscCall(DMPlexVecRestoreClosure(dmTemperature,PETSC_NULL_SECTION,MEF90HeatXferCtx%boundaryFluxLocal,setPointID(cell),boundaryFluxDof,ierr))
                  End Do ! cell

                  PetscCall(MEF90ElementDestroy(elem,ierr))
                  PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
               End If ! faceSetOptions%boundaryFlux
            End If ! pointIS
            PetscCall(ISDestroy(setPointIS,ierr))
            PetscCallMPI(MPI_AllReduce(mySurfaceWork,surfaceWork(set),1,MPIU_SCALAR,MPI_SUM,MEF90HeatXferCtx%MEF90Ctx%comm,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))
   End Subroutine MEF90HeatXFerEnergy


#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXFerIFunction"
!!!
!!!  
!!!  MEF90HeatXFerIFunction:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022 Blaise Bourdin bourdin@mcmaster.ca
!!!      2022 Alexis Marboeuf marboeua@mcmaster.ca
!!!

   Subroutine MEF90HeatXFerIFunction(tempTS,time,x,xdot,F,MEF90HeatXferCtx,ierr)
      Type(tTS),Intent(IN)                            :: tempTS
      PetscReal,Intent(IN)                            :: time
      Type(tVec),Intent(IN)                           :: x,xdot
      Type(tVec),Intent(INOUT)                        :: F
      Type(MEF90HeatXferCtx_Type),Intent(IN)          :: MEF90HeatXferCtx
      PetscErrorCode,Intent(INOUT)                    :: ierr
         
      Type(tDM)                                       :: dmTemperature
      Type(tVec)                                      :: locTemperature,locTemperatureDot,locF
      Type(tIS)                                       :: setIS,setPointIS
      PetscInt,dimension(:),Pointer                   :: setID,setPointID
      PetscInt                                        :: set,QuadratureOrder
      PetscInt                                        :: cell, iDof,jDof,iGauss
      Type(MEF90HeatXferCellSetOptions_Type),pointer  :: cellSetOptions
      Type(MEF90HeatXferFaceSetOptions_Type),pointer  :: faceSetOptions
      Type(MEF90_MATPROP),pointer                     :: matpropSet
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer   :: elem
      Type(MEF90ElementType)                          :: elementType
      DMPolytopeType                                  :: cellType
      Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90CtxGlobalOptions
      Type(MEF90HeatXferGlobalOptions_Type),pointer   :: MEF90HeatXferGlobalOptions
      PetscReal,Dimension(:),Pointer                  :: temperatureDof,temperatureDotDof,fluxDof,boundaryFluxDof,externalTemperatureDof,residualDof
      Type(MEF90_VECT)                                :: advectionVec
      
      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90HeatXferCtx%MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr))
      PetscCall(PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr))
      PetscCall(VecGetDM(MEF90HeatXferCtx%temperatureLocal,dmTemperature,ierr))

      PetscCall(DMGetLocalVector(dmTemperature,locTemperature,ierr))
      PetscCall(DMGetLocalVector(dmTemperature,locTemperatureDot,ierr))
      PetscCall(DMGetLocalVector(dmTemperature,locF,ierr))
      PetscCall(MEF90VecGlobalToLocalConstraint(x,MEF90HeatXferCtx%temperatureLocal,locTemperature,ierr))
      PetscCall(VecSet(locTemperatureDot,0.0_Kr,ierr))
      PetscCall(DMGlobalToLocalBegin(dmTemperature,xdot,INSERT_VALUES,locTemperatureDot,ierr))
      PetscCall(DMGlobalToLocalEnd(dmTemperature,xdot,INSERT_VALUES,locTemperatureDot,ierr))

      PetscCall(VecSet(F,0.0_Kr,ierr))
      PetscCall(VecSet(locF,0.0_Kr,ierr))
         
      !!! cell-based contributions
      PetscCall(DMGetLabelIdIS(dmTemperature,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm,setIS,ierr))
      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         Do set = 1,size(setID)
            PetscCall(DMGetStratumIS(dmTemperature,MEF90CellSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
               PetscCall(PetscBagGetDataMEF90MatProp(MEF90HeatXferCtx%MaterialPropertiesBag(set),matpropSet,ierr))
               PetscCall(PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr))
         
               PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(DMPlexGetCellType(dmTemperature,setPointID(1),cellType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellType,elementType,ierr))
               QuadratureOrder = elementType%order * 2
               PetscCall(MEF90ElementCreate(dmTemperature,setPointIS,elem,QuadratureOrder,elementType,ierr))

               Allocate(residualDof(size(elem(1)%BF(:,1))))

               Do cell = 1, size(setPointID)
                  residualDof = 0.0_Kr
                  !!! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                  !!! If this happens, we will need to protect this loop
                  PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,locTemperatureDot,setPointID(cell),temperatureDotDof,ierr))
                  Do iGauss = 1,size(elem(cell)%Gauss_C)
                     Do jDof = 1, size(elem(cell)%BF(:,1))
                        Do iDof = 1, size(elem(cell)%BF(:,1))
                           residualDof(jDof) = residualDof(jDof) + matpropSet%Density*matpropSet%SpecificHeat*temperatureDotDof(iDof)*elem(cell)%BF(iDof,iGauss)*elem(cell)%BF(jDof,iGauss)*elem(cell)%Gauss_C(iGauss)
                        End Do ! iDof
                     End Do ! jDof
                  End Do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmTemperature,PETSC_NULL_SECTION,locTemperatureDot,setPointID(cell),temperatureDotDof,ierr))
                  PetscCall(DMPlexVecSetClosure(dmTemperature,PETSC_NULL_SECTION,locF,setPointID(cell),residualDof,ADD_VALUES,ierr))
               End Do ! cell

               Do cell = 1, size(setPointID)
                  residualDof = 0.0_Kr
                  !!! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                  !!! If this happens, we will need to protect this loop
                  PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,locTemperature,setPointID(cell),temperatureDof,ierr))
                  Do iGauss = 1,size(elem(cell)%Gauss_C)
                     Do jDof = 1, size(elem(cell)%BF(:,1))
                        Do iDof = 1, size(elem(cell)%BF(:,1))
                           residualDof(jDof) = residualDof(jDof) + (matpropSet%ThermalConductivity*temperatureDof(iDof)*elem(cell)%Grad_BF(iDof,iGauss) .dotP. elem(cell)%Grad_BF(jDof,iGauss))*elem(cell)%Gauss_C(iGauss)
                        End Do ! iDof
                     End Do ! jDof
                  End Do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmTemperature,PETSC_NULL_SECTION,locTemperature,setPointID(cell),temperatureDof,ierr))
                  PetscCall(DMPlexVecSetClosure(dmTemperature,PETSC_NULL_SECTION,locF,setPointID(cell),residualDof,ADD_VALUES,ierr))
               End Do ! cell

               If (norm2(cellSetOptions%advectionVector) /= 0.0_Kr) Then
                  advectionVec = cellSetOptions%advectionVector
                  Do cell = 1, size(setPointID)
                     residualDof = 0.0_Kr
                     !!! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !!! If this happens, we will need to protect this loop
                     PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,locTemperature,setPointID(cell),temperatureDof,ierr))
                     Do iGauss = 1,size(elem(cell)%Gauss_C)
                        Do jDof = 1, size(elem(cell)%BF(:,1))
                           Do iDof = 1, size(elem(cell)%BF(:,1))
                              residualDof(jDof) = residualDof(jDof) - (matpropSet%Density*matpropSet%SpecificHeat*advectionVec .dotP. temperatureDof(iDof)*elem(cell)%Grad_BF(iDof,iGauss))*elem(cell)%BF(jDof,iGauss)*elem(cell)%Gauss_C(iGauss)
                           End Do ! iDof
                        End Do ! jDof
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmTemperature,PETSC_NULL_SECTION,locTemperature,setPointID(cell),temperatureDof,ierr))
                     PetscCall(DMPlexVecSetClosure(dmTemperature,PETSC_NULL_SECTION,locF,setPointID(cell),residualDof,ADD_VALUES,ierr))
                  End Do ! cell
               End If ! cellSetOptions%advectionVector

               If (cellSetOptions%flux /= 0.0_Kr) Then
                  Do cell = 1, size(setPointID)
                     residualDof = 0.0_Kr
                     !!! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !!! If this happens, we will need to protect this loop
                     PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,MEF90HeatXferCtx%fluxLocal,setPointID(cell),fluxDof,ierr))
                     Do iGauss = 1,size(elem(cell)%Gauss_C)
                        Do iDof = 1, size(elem(cell)%BF(:,1))
                           residualDof(iDof) = residualDof(iDof) - fluxDof(iDof)*elem(cell)%BF(iDof,iGauss)*elem(cell)%Gauss_C(iGauss)
                        End Do ! iDof
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmTemperature,PETSC_NULL_SECTION,MEF90HeatXferCtx%fluxLocal,setPointID(cell),fluxDof,ierr))
                     PetscCall(DMPlexVecSetClosure(dmTemperature,PETSC_NULL_SECTION,locF,setPointID(cell),residualDof,ADD_VALUES,ierr))
                  End Do ! cell
               End If ! cellSetOptions%flux

               DeAllocate(residualDof)
               PetscCall(MEF90ElementDestroy(elem,ierr))
               PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
            End If ! pointIS
            PetscCall(ISDestroy(setPointIS,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))

      !!! face-based contributions
      PetscCall(DMGetLabelIdIS(dmTemperature,MEF90FaceSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm,setIS,ierr))
      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         Do set = 1,size(setID)
            PetscCall(DMGetStratumIS(dmTemperature,MEF90FaceSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
              PetscCall(PetscBagGetDataMEF90HeatXferCtxFaceSetOptions(MEF90HeatXferCtx%FaceSetOptionsBag(set),faceSetOptions,ierr))

               If (faceSetOptions%boundaryFlux /= 0.0_Kr) Then
                  PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
                  PetscCall(DMPlexGetCellType(dmTemperature,setPointID(1),cellType,ierr))
                  PetscCall(MEF90ElementGetTypeBoundary(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellType,elementType,ierr))
                  QuadratureOrder = elementType%order * 2
                  PetscCall(MEF90ElementCreate(dmTemperature,setPointIS,elem,QuadratureOrder,elementType,ierr))
                  Allocate(residualDof(size(elem(1)%BF(:,1))))
                  Do cell = 1, size(setPointID)
                     residualDof = 0.0_Kr
                     !!! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !!! If this happens, we will need to protect this loop
                     PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,MEF90HeatXFerCtx%boundaryFluxLocal,setPointID(cell),boundaryFluxDof,ierr))
                     Do iGauss = 1,size(elem(cell)%Gauss_C)
                        Do iDof = 1, size(elem(cell)%BF(:,1))
                           residualDof(iDof) = residualDof(iDof) - boundaryFluxDof(iDof)*elem(cell)%BF(iDof,iGauss)*elem(cell)%Gauss_C(iGauss)
                        End Do ! iDof
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmTemperature,PETSC_NULL_SECTION,MEF90HeatXFerCtx%boundaryFluxLocal,setPointID(cell),boundaryFluxDof,ierr))
                     PetscCall(DMPlexVecSetClosure(dmTemperature,PETSC_NULL_SECTION,locF,setPointID(cell),residualDof,ADD_VALUES,ierr))
                  End Do ! cell
                  DeAllocate(residualDof)
                  PetscCall(MEF90ElementDestroy(elem,ierr))
                  PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
               End If ! faceSetOptions%boundaryFlux

               If (faceSetOptions%surfaceThermalConductivity /= 0.0_Kr) Then
                  PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
                  PetscCall(DMPlexGetCellType(dmTemperature,setPointID(1),cellType,ierr))
                  PetscCall(MEF90ElementGetTypeBoundary(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellType,elementType,ierr))
                  QuadratureOrder = elementType%order * 2
                  PetscCall(MEF90ElementCreate(dmTemperature,setPointIS,elem,QuadratureOrder,elementType,ierr))
                  Allocate(residualDof(size(elem(1)%BF(:,1))))
                  Do cell = 1, size(setPointID)
                     residualDof = 0.0_Kr
                     !!! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !!! If this happens, we will need to protect this loop
                     PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,locTemperature,setPointID(cell),temperatureDof,ierr))
                     PetscCall(DMPlexVecGetClosure(dmTemperature,PETSC_NULL_SECTION,MEF90HeatXFerCtx%externalTemperatureLocal,setPointID(cell),externalTemperatureDof,ierr))
                     Do iGauss = 1,size(elem(cell)%Gauss_C)
                        Do jDof = 1, size(elem(cell)%BF(:,1))
                           Do iDof = 1, size(elem(cell)%BF(:,1))
                              residualDof(jDof) = residualDof(jDof) + faceSetOptions%surfaceThermalConductivity*temperatureDof(iDof)*elem(cell)%BF(iDof,iGauss)*elem(cell)%BF(jDof,iGauss)*elem(cell)%Gauss_C(iGauss)
                           End Do ! iDof
                           residualDof(jDof) = residualDof(jDof) - faceSetOptions%surfaceThermalConductivity*externalTemperatureDof(jDof)*elem(cell)%BF(jDof,iGauss)*elem(cell)%Gauss_C(iGauss)
                        End Do ! jDof
                     End Do ! iGauss
                     PetscCall(DMPlexVecRestoreClosure(dmTemperature,PETSC_NULL_SECTION,MEF90HeatXFerCtx%externalTemperatureLocal,setPointID(cell),externalTemperatureDof,ierr))
                     PetscCall(DMPlexVecRestoreClosure(dmTemperature,PETSC_NULL_SECTION,locTemperature,setPointID(cell),temperatureDof,ierr))
                     PetscCall(DMPlexVecSetClosure(dmTemperature,PETSC_NULL_SECTION,locF,setPointID(cell),residualDof,ADD_VALUES,ierr))
                  End Do ! cell
                  DeAllocate(residualDof)
                  PetscCall(MEF90ElementDestroy(elem,ierr))
                  PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
               End If ! faceSetOptions%surfaceThermalConductivity
            End If ! pointIS
            PetscCall(ISDestroy(setPointIS,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))
      PetscCall(DMLocalToGlobalBegin(dmTemperature,locF,ADD_VALUES,F,ierr))
      PetscCall(DMLocalToGlobalEnd(dmTemperature,locF,ADD_VALUES,F,ierr))
      PetscCall(DMRestoreLocalVector(dmTemperature,locTemperature,ierr))
      PetscCall(DMRestoreLocalVector(dmTemperature,locTemperatureDot,ierr))
      PetscCall(DMRestoreLocalVector(dmTemperature,locF,ierr))
   End Subroutine MEF90HeatXFerIFunction
   
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferIJacobian"
!!!
!!!  
!!!  MEF90HeatXferIJacobian:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022 Blaise Bourdin bourdin@mcmaster.ca
!!!      2022 Alexis Marboeuf marboeua@mcmaster.ca
!!!

   Subroutine MEF90HeatXferIJacobian(tempTS,t,x,xdot,shift,A,M,MEF90HeatXferCtx,ierr)
      Type(tTS),Intent(IN)                            :: tempTS
      PetscReal,Intent(IN)                            :: t
      Type(tVec),Intent(IN)                           :: x,xdot
      PetscReal,Intent(IN)                            :: shift
      Type(tMat),Intent(INOUT)                        :: A,M
      Type(MEF90HeatXferCtx_Type),Intent(IN)          :: MEF90HeatXferCtx
      PetscErrorCode,Intent(INOUT)                    :: ierr  

      Type(tDM)                                       :: dmTemperature
      Type(tIS)                                       :: setIS,setPointIS
      PetscInt,dimension(:),Pointer                   :: setID,setPointID
      PetscInt                                        :: set,QuadratureOrder
      PetscInt                                        :: cell,iDof,jDof,iGauss,nbDof
      Type(MEF90HeatXferCellSetOptions_Type),pointer  :: cellSetOptions
      Type(MEF90HeatXferFaceSetOptions_Type),pointer  :: faceSetOptions
      Type(MEF90_MATPROP),pointer                     :: matpropSet
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer   :: elem
      Type(MEF90ElementType)                          :: elementType
      DMPolytopeType                                  :: cellType
      Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90CtxGlobalOptions
      Type(MEF90HeatXferGlobalOptions_Type),pointer   :: MEF90HeatXferGlobalOptions
      PetscReal,Dimension(:),Pointer                  :: matDof
      Type(MEF90_VECT)                                :: advectionVec
      
      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90HeatXferCtx%MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr))
      PetscCall(PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr))
      PetscCall(VecGetDM(MEF90HeatXferCtx%temperatureLocal,dmTemperature,ierr))

      PetscCall(MatZeroEntries(A,ierr))
         
      !!! cell-based gradient contributions
      PetscCall(DMGetLabelIdIS(dmTemperature,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm,setIS,ierr))
      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         Do set = 1,size(setID)
            PetscCall(DMGetStratumIS(dmTemperature,MEF90CellSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
               PetscCall(PetscBagGetDataMEF90MatProp(MEF90HeatXferCtx%MaterialPropertiesBag(set),matpropSet,ierr))
               PetscCall(PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr))
         
               PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(DMPlexGetCellType(dmTemperature,setPointID(1),cellType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellType,elementType,ierr))
               QuadratureOrder = elementType%order * 2
               PetscCall(MEF90ElementCreate(dmTemperature,setPointIS,elem,QuadratureOrder,elementType,ierr))

               nbDof = size(elem(1)%BF(:,1))
               Allocate(matDof(nbDof*nbDof))

               Do cell = 1, size(setPointID)
                  matDof = 0.0_Kr
                  !!! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                  !!! If this happens, we will need to protect this loop
                  Do iGauss = 1,size(elem(cell)%Gauss_C)
                     Do jDof = 0, nbDof-1
                        Do iDof = 1, nbDof
                           matDof(jDof*nbDof+iDof) = matDof(jDof*nbDof+iDof) + shift*matpropSet%Density*matpropSet%SpecificHeat*elem(cell)%BF(iDof,iGauss)*elem(cell)%BF(jDof+1,iGauss)*elem(cell)%Gauss_C(iGauss)
                        End Do ! iDof
                     End Do ! jDof
                  End Do ! iGauss
                  PetscCall(DMPlexMatSetClosure(dmTemperature,PETSC_NULL_SECTION,PETSC_NULL_SECTION,A,setPointID(cell),matDof,ADD_VALUES,ierr))
               End Do ! cell

               Do cell = 1, size(setPointID)
                  matDof = 0.0_Kr
                  !!! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                  !!! If this happens, we will need to protect this loop
                  Do iGauss = 1,size(elem(cell)%Gauss_C)
                     Do jDof = 0, nbDof-1
                        Do iDof = 1, nbDof
                           matDof(jDof*nbDof+iDof) = matDof(jDof*nbDof+iDof) + (matpropSet%ThermalConductivity*elem(cell)%Grad_BF(iDof,iGauss) .dotP. elem(cell)%Grad_BF(jDof+1,iGauss))*elem(cell)%Gauss_C(iGauss)
                        End Do ! iDof
                     End Do ! jDof
                  End Do ! iGauss
                  PetscCall(DMPlexMatSetClosure(dmTemperature,PETSC_NULL_SECTION,PETSC_NULL_SECTION,A,setPointID(cell),matDof,ADD_VALUES,ierr))
               End Do ! cell

               If (norm2(cellSetOptions%advectionVector) /= 0.0_Kr) Then
                  advectionVec = cellSetOptions%advectionVector
                  Do cell = 1, size(setPointID)
                     matDof = 0.0_Kr
                     !!! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !!! If this happens, we will need to protect this loop
                     Do iGauss = 1,size(elem(cell)%Gauss_C)
                        Do jDof = 0, nbDof-1
                           Do iDof = 1, nbDof
                              matDof(jDof*nbDof+iDof) = matDof(jDof*nbDof+iDof) - (matpropSet%Density*matpropSet%SpecificHeat*advectionVec .dotP. elem(cell)%Grad_BF(iDof,iGauss))*elem(cell)%BF(jDof+1,iGauss)*elem(cell)%Gauss_C(iGauss)
                           End Do ! iDof
                        End Do ! jDof
                     End Do ! iGauss
                     PetscCall(DMPlexMatSetClosure(dmTemperature,PETSC_NULL_SECTION,PETSC_NULL_SECTION,A,setPointID(cell),matDof,ADD_VALUES,ierr))
                  End Do ! cell
               End If ! cellSetOptions%advectionVector

               DeAllocate(matDof)
               PetscCall(MEF90ElementDestroy(elem,ierr))
               PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
            End If ! pointIS
            PetscCall(ISDestroy(setPointIS,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))

      !!! face-based energies
      PetscCall(DMGetLabelIdIS(dmTemperature,MEF90FaceSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm,setIS,ierr))
      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         Do set = 1,size(setID)
            PetscCall(DMGetStratumIS(dmTemperature,MEF90FaceSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then
              PetscCall(PetscBagGetDataMEF90HeatXferCtxFaceSetOptions(MEF90HeatXferCtx%FaceSetOptionsBag(set),faceSetOptions,ierr))

               If (faceSetOptions%surfaceThermalConductivity /= 0.0_Kr) Then
                  PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
                  PetscCall(DMPlexGetCellType(dmTemperature,setPointID(1),cellType,ierr))
                  PetscCall(MEF90ElementGetTypeBoundary(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellType,elementType,ierr))
                  QuadratureOrder = elementType%order * 2
                  PetscCall(MEF90ElementCreate(dmTemperature,setPointIS,elem,QuadratureOrder,elementType,ierr))

                  nbDof = size(elem(1)%BF(:,1))
                  Allocate(matDof(nbDof*nbDof))

                  Do cell = 1, size(setPointID)
                     matDof = 0.0_Kr
                     !!! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !!! If this happens, we will need to protect this loop
                     Do iGauss = 1,size(elem(cell)%Gauss_C)
                        Do jDof = 0, nbDof-1
                           Do iDof = 1, nbDof
                              matDof(jDof*nbDof+iDof) = matDof(jDof*nbDof+iDof) + faceSetOptions%surfaceThermalConductivity*elem(cell)%BF(iDof,iGauss)*elem(cell)%BF(jDof+1,iGauss)*elem(cell)%Gauss_C(iGauss)
                           End Do ! iDof
                        End Do ! jDof
                     End Do ! iGauss
                     PetscCall(DMPlexMatSetClosure(dmTemperature,PETSC_NULL_SECTION,PETSC_NULL_SECTION,A,setPointID(cell),matDof,ADD_VALUES,ierr))
                  End Do ! cell

                  DeAllocate(matDof)
                  PetscCall(MEF90ElementDestroy(elem,ierr))
                  PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
               End If ! faceSetOptions%surfaceThermalConductivity
            End If ! pointIS
            PetscCall(ISDestroy(setPointIS,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))   
      PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr))
      PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr))
      PetscCall(MatCopy(A,M,SAME_NONZERO_PATTERN,ierr))
   End Subroutine MEF90HeatXferIJacobian
End Module MEF90_APPEND(m_MEF90_HeatXferAssembly,MEF90_DIM)D
