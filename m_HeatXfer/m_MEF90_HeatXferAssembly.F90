#include "../MEF90/mef90.inc"
Module MEF90_APPEND(m_MEF90_HeatXferAssembly,MEF90_DIM)D
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_HeatXferCtx
   Implicit none

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferOperator"
!!!
!!!  
!!!  MEF90HeatXferOperator: Build the operator. When called in SNES, the solution time should always match the target time, 
!!!                         so there is no need for interpolation of the fluxes, external, and boundary values
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90_APPEND(MEF90HeatXferOperator,MEF90_DIM)D(snesTemp,x,residual,MEF90HeatXferCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesTemp
      Type(Vec),Intent(IN)                               :: x
      Type(Vec),Intent(INOUT)                            :: residual
      Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(SectionReal)                                  :: xSec,residualSec
      Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),pointer                        :: matpropSet
      Type(MEF90HeatXferCellSetOptions_Type),pointer     :: cellSetOptions
      Type(MEF90HeatXferVertexSetOptions_Type),pointer   :: vertexSetOptions
      Type(VecScatter)                                   :: ScatterSecToVec
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elem
      Type(MEF90Element_Type)                            :: elemType
      PetscInt                                           :: dof
      PetscReal,Dimension(:),Pointer                     :: residualPtr,xPtr
      
      Call DMMeshGetSectionReal(MEF90HeatXferCtx%DM,'default',xSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(xSec,residualSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90HeatXferCtx%DM,xSec,ScatterSecToVec,ierr);CHKERRQ(ierr)
   
      Call SectionRealSet(residualSec,0.0_Kr,ierr);CHKERRQ(ierr)
      Call VecSet(residual,0.0_Kr,ierr);CHKERRQ(ierr)
   
      Call SectionRealToVec(xSec,ScatterSecToVec,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call DMMeshGetStratumIS(MEF90HeatXferCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
   
         Call PetscBagGetDataMEF90MatProp(MEF90HeatXferCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
   
         elemType = MEF90_knownElements(cellSetOptions%ElemTypeShortID)
         QuadratureOrder = elemType%order * 2
         Call MEF90Element_Create(MEF90HeatXferCtx%DM,setIS,elem,QuadratureOrder,CellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)
         Call MEF90DiffusionOperatorSet(residualSec,MEF90HeatXferCtx%DM,xSec,setIS,matpropSet%ThermalConductivity,cellSetOptions%SurfaceThermalConductivity,elem,elemType,ierr);CHKERRQ(ierr)
         !!! Since the flux vector is defined at the cells, it must be assembled in the RHS
         
         Call MEF90Element_Destroy(elem,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do     
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
      !!! "Ghost update" for the residual Section
      Call SectionRealComplete(residualSec,ierr);CHKERRQ(ierr)
      !!! Scatter back from SectionReal to Vec
      Call SectionRealToVec(residualSec,ScatterSecToVec,SCATTER_FORWARD,residual,ierr);CHKERRQ(ierr)
      
      
      !!!
      !!! Account for BC entries in the residual
      !!!
      Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90HeatXferCtxVertexSetOptions(MEF90HeatXferCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
         If (vertexSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(MEF90HeatXferCtx%DM,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
            Call ISGetIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
            Call DMMeshISCreateISglobaldof(MEF90HeatXferCtx%DM,setIS,0,setISdof,ierr);CHKERRQ(ierr)
            Call ISGetIndicesF90(setISdof,setdofIdx,ierr);CHKERRQ(ierr)
            Allocate(residualPtr(size(setIdx)))
            Do dof = 1, size(setIdx)
               Call SectionRealRestrict(xSec,setIdx(dof),xPtr,ierr);CHKERRQ(ierr)
               !!! At this point, I have two choices:
               !!! I can use SectionRealUpdate (but NOT complete, since this would add values)
               !!! or I can do a VecSetValues, provided that I have build  setISdof using DMMeshISCreateISglobaldof
               !!! going for the later               
               residualPtr(dof) = xPtr(1)
               Call SectionRealRestore(xSec,setIdx(dof),xPtr,ierr);CHKERRQ(ierr)
            End Do
            Call VecSetValues(residual,size(setIdx),setdofIdx,residualPtr,INSERT_VALUES,ierr);CHKERRQ(ierr)
            DeAllocate(residualPtr)
            Call ISRestoreIndicesF90(setISdof,setdofIdx,ierr);CHKERRQ(ierr)
            Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
            Call ISRestoreIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If
      End Do

      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
   
      Call VecAssemblyBegin(residual,ierr)
      Call VecAssemblyEnd(residual,ierr)
      Call SectionRealDestroy(residualSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90_APPEND(MEF90HeatXferOperator,MEF90_DIM)D
   
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferBilinearForm"
!!!
!!!  
!!!  MEF90HeatXferBilinearForm:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90_APPEND(MEF90HeatXferBilinearForm,MEF90_DIM)D(snesTemp,x,A,M,flg,MEF90HeatXferCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesTemp
      Type(Vec),Intent(IN)                               :: x
      Type(Mat),Intent(INOUT)                            :: A,M
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
      PetscErrorCode,Intent(OUT)                         :: ierr  
         
      Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),pointer                        :: matpropSet
      Type(MEF90HeatXferCellSetOptions_Type),pointer     :: cellSetOptions
      Type(MEF90HeatXferVertexSetOptions_Type),pointer   :: vertexSetOptions
      Type(DM)                                           :: mesh
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elem
      Type(MEF90Element_Type)                            :: elemType
      
      Call MatZeroEntries(A,ierr);CHKERRQ(ierr)
      Call SNESGetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
   
         Call PetscBagGetDataMEF90MatProp(MEF90HeatXferCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         elemType = MEF90_knownElements(cellSetOptions%ElemTypeShortID)
   
         QuadratureOrder = 2 * elemType%order
         Call MEF90Element_Create(mesh,setIS,elem,QuadratureOrder,CellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)
         Call MEF90DiffusionBilinearFormSet(A,mesh,setIS,matpropSet%ThermalConductivity,cellSetOptions%SurfaceThermalConductivity,elem,elemType,ierr);CHKERRQ(ierr)
         Call MEF90Element_Destroy(elem,ierr)
   
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do     
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
      Call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
   
      !!!
      !!! Boundary conditions
      !!!
      Call DMmeshGetLabelIdIS(mesh,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90HeatXferCtxVertexSetOptions(MEF90HeatXferCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
         If (vertexSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(mesh,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
            Call DMMeshISCreateISglobaldof(mesh,setIS,0,setISdof,ierr);CHKERRQ(ierr)
            Call MatZeroRowsColumnsIS(A,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
         End If
      End Do
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      
      flg = SAME_NONZERO_PATTERN
   End Subroutine MEF90_APPEND(MEF90HeatXferBilinearForm,MEF90_DIM)D

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferRHS"
!!!
!!!  
!!!  MEF90HeatXferRHS: Build the time dependent RHS for the SNES
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90_APPEND(MEF90HeatXferRHS,MEF90_DIM)D(rhs,t,MEF90HeatXferCtx,ierr)
      Type(Vec),Intent(IN)                               :: rhs
      PetscReal,Intent(IN)                               :: t
      Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof
      PetscInt,dimension(:),Pointer                      :: setID,SetIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder,numCell
      Type(MEF90HeatXferCellSetOptions_Type),pointer     :: cellSetOptions
      Type(MEF90HeatXferVertexSetOptions_Type),pointer   :: vertexSetOptions
      Type(MEF90_MATPROP),pointer                        :: matpropSet
      PetscReal,Dimension(:),Pointer                     :: tmpPtr1,tmpPtr2
      !PetscReal                                          :: scaling
      Type(VecScatter)                                   :: ScatterSecToVec,ScatterSecToVecCell
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elem
      Type(MEF90Element_Type)                            :: elemType
      Type(Vec)                                          :: flux
      Type(Vec)                                          :: modifiedFluxVec
      Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90CtxGlobalOptions
      Type(MEF90HeatXferGlobalOptions_Type),pointer      :: MEF90HeatXferGlobalOptions
      Type(SectionReal)                                  :: fluxPreviousSec,fluxTargetSec
      Type(SectionReal)                                  :: externalTemperaturePreviousSec,externalTemperatureTargetSec
      Type(SectionReal)                                  :: boundaryTemperaturePreviousSec,boundaryTemperatureTargetSec
      Type(SectionReal)                                  :: rhsSec,modifiedFluxSec
      PetscReal,Dimension(:),Pointer                     :: fluxPreviousPtr,fluxTargetPtr
      PetscReal,Dimension(:),Pointer                     :: externalTemperaturePreviousPtr,externalTemperatureTargetPtr
      PetscReal,Dimension(:),Pointer                     :: modifiedFluxPtr,rhsPtr
      PetscInt                                           :: dof
   
   
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90HeatXferCtx%MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr);CHKERRQ(ierr)

      !!! Create SectionReal to hold fields   
      Call DMMeshGetSectionReal(MEF90HeatXferCtx%DM,'default',rhsSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90HeatXferCtx%DM,rhsSec,ScatterSecToVec,ierr);CHKERRQ(ierr)
      Call DMMeshGetSectionReal(MEF90HeatXferCtx%CellDM,'default',modifiedFluxSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90HeatXferCtx%CellDM,modifiedFluxSec,ScatterSecToVecCell,ierr);CHKERRQ(ierr)

      !!! Scatter necessary fields from the Section to the matching Vec
      Call SectionRealDuplicate(modifiedFluxSec,fluxTargetSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(fluxTargetSec,ScatterSecToVecCell,SCATTER_REVERSE,MEF90HeatXferCtx%fluxTarget,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(modifiedFluxSec,externalTemperatureTargetSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(externalTemperatureTargetSec,ScatterSecToVecCell,SCATTER_REVERSE,MEF90HeatXferCtx%externalTemperatureTarget,ierr);CHKERRQ(ierr)

      Call SectionRealDuplicate(modifiedFluxSec,fluxPreviousSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(fluxPreviousSec,ScatterSecToVecCell,SCATTER_REVERSE,MEF90HeatXferCtx%fluxPrevious,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(modifiedFluxSec,externalTemperaturePreviousSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(externalTemperaturePreviousSec,ScatterSecToVecCell,SCATTER_REVERSE,MEF90HeatXferCtx%externalTemperaturePrevious,ierr);CHKERRQ(ierr)

      
      Call SectionRealSet(rhsSec,0.0_Kr,ierr);CHKERRQ(ierr)
      Call VecSet(rhs,0.0_Kr,ierr);CHKERRQ(ierr)
   
      Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

      Allocate(modifiedFluxPtr(1))
      Do set = 1,size(setID)
         Call DMMeshGetStratumIS(MEF90HeatXferCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call ISGetIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)

         Call PetscBagGetDataMEF90MatProp(MEF90HeatXferCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
   
         !!! Modified flux is flux + surfaceThermalConductivity * refTemp      
         If (t == MEF90HeatXferCtx%timePrevious) Then
            Do dof = 1,size(setIdx)
               Call SectionRealRestrict(fluxPreviousSec,setIdx(dof),tmpPtr1,ierr);CHKERRQ(ierr)
               Call SectionRealRestrict(externalTemperaturePreviousSec,setIdx(dof),tmpPtr2,ierr);CHKERRQ(ierr)
               modifiedFluxPtr =  tmpPtr1 + cellSetOptions%SurfaceThermalConductivity * tmpPtr2
               Call SectionRealUpdate(modifiedFluxSec,setIdx(dof),modifiedFluxPtr,INSERT_VALUES,ierr);CHKERRQ(ierr)
               Call SectionRealRestore(externalTemperaturePreviousSec,setIdx(dof),tmpPtr2,ierr);CHKERRQ(ierr)
               Call SectionRealRestore(fluxPreviousSec,setIdx(dof),tmpPtr1,ierr);CHKERRQ(ierr)
            End Do
         Else If (t == MEF90HeatXferCtx%timeTarget) Then
            Do dof = 1,size(setIdx)
               Call SectionRealRestrict(fluxTargetSec,setIdx(dof),tmpPtr1,ierr);CHKERRQ(ierr)
               Call SectionRealRestrict(externalTemperatureTargetSec,setIdx(dof),tmpPtr2,ierr);CHKERRQ(ierr)
               modifiedFluxPtr =  tmpPtr1 + cellSetOptions%SurfaceThermalConductivity * tmpPtr2
               Call SectionRealUpdate(modifiedFluxSec,setIdx(dof),modifiedFluxPtr,INSERT_VALUES,ierr);CHKERRQ(ierr)
               Call SectionRealRestore(externalTemperatureTargetSec,setIdx(dof),tmpPtr2,ierr);CHKERRQ(ierr)
               Call SectionRealRestore(fluxTargetSec,setIdx(dof),tmpPtr1,ierr);CHKERRQ(ierr)
            End Do
         Else  !!! tmpPtr1 between the steps
            Do dof = 1,size(setIdx)
               Call SectionRealRestrict(fluxPreviousSec,setIdx(dof),tmpPtr1,ierr);CHKERRQ(ierr)
               Call SectionRealRestrict(fluxTargetSec,setIdx(dof),tmpPtr2,ierr);CHKERRQ(ierr)
               modifiedFluxPtr =  (MEF90HeatXferCtx%timeTarget - t) * tmpPtr1 + (t - MEF90HeatXferCtx%timePrevious) * tmpPtr2
               Call SectionRealRestore(fluxTargetSec,setIdx(dof),tmpPtr2,ierr);CHKERRQ(ierr)
               Call SectionRealRestore(fluxPreviousSec,setIdx(dof),tmpPtr1,ierr);CHKERRQ(ierr)

               Call SectionRealRestrict(externalTemperaturePreviousSec,setIdx(dof),tmpPtr1,ierr);CHKERRQ(ierr)
               Call SectionRealRestrict(externalTemperatureTargetSec,setIdx(dof),tmpPtr2,ierr);CHKERRQ(ierr)
               modifiedFluxPtr =  modifiedFluxPtr + cellSetOptions%SurfaceThermalConductivity &
                                  * ((MEF90HeatXferCtx%timeTarget - t) * tmpPtr1 + (t - MEF90HeatXferCtx%timePrevious) * tmpPtr2)
               modifiedFluxPtr = modifiedFluxPtr / (MEF90HeatXferCtx%timeTarget - MEF90HeatXferCtx%timePrevious)
               Call SectionRealUpdate(modifiedFluxSec,setIdx(dof),modifiedFluxPtr,INSERT_VALUES,ierr);CHKERRQ(ierr)
               Call SectionRealRestore(externalTemperatureTargetSec,setIdx(dof),tmpPtr2,ierr);CHKERRQ(ierr)
               Call SectionRealRestore(externalTemperaturePreviousSec,setIdx(dof),tmpPtr1,ierr);CHKERRQ(ierr)
            End Do
         End If
          
         elemType = MEF90_knownElements(cellSetOptions%ElemTypeShortID)
         QuadratureOrder = elemType%order * 2
         Call MEF90Element_Create(MEF90HeatXferCtx%DM,setIS,elem,QuadratureOrder,CellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)   
         Call MEF90DiffusionRHSSetCell(rhssec,MEF90HeatXferCtx%DM,modifiedFluxSec,setIS,elem,elemType,ierr);CHKERRQ(ierr)
         Call MEF90Element_Destroy(elem,ierr)

         Call ISRestoreIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do     
      DeAllocate(modifiedFluxPtr)
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
      !!! "Ghost update" for the residual Section
      Call SectionRealComplete(rhsSec,ierr);CHKERRQ(ierr)
      !!! Scatter back from SectionReal to Vec
      Call SectionRealToVec(rhsSec,ScatterSecToVec,SCATTER_FORWARD,rhs,ierr);CHKERRQ(ierr)

      Call SectionRealDestroy(modifiedFluxSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(fluxPreviousSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(externalTemperaturePreviousSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(fluxTargetSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(externalTemperatureTargetSec,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVecCell,ierr);CHKERRQ(ierr)
      
      !!!
      !!! Account for BC entries
      !!!
      Call SectionRealDuplicate(rhsSec,boundaryTemperaturePreviousSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(boundaryTemperaturePreviousSec,ScatterSecToVec,SCATTER_REVERSE,MEF90HeatXferCtx%boundaryTemperaturePrevious,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(rhsSec,boundaryTemperatureTargetSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(boundaryTemperatureTargetSec,ScatterSecToVec,SCATTER_REVERSE,MEF90HeatXferCtx%boundaryTemperatureTarget,ierr);CHKERRQ(ierr)
      Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90HeatXferCtxVertexSetOptions(MEF90HeatXferCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
         If (vertexSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(MEF90HeatXferCtx%DM,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
            Call ISGetIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
            Call DMMeshISCreateISglobaldof(MEF90HeatXferCtx%DM,setIS,0,setISdof,ierr);CHKERRQ(ierr)
            Call ISGetIndicesF90(setISdof,setdofIdx,ierr);CHKERRQ(ierr)
            Allocate(rhsPtr(size(setIdx)))
            If (t == MEF90HeatXferCtx%timePrevious) Then
               Do dof = 1,size(setIdx)
                  Call SectionRealRestrict(boundaryTemperaturePreviousSec,setIdx(dof),tmpPtr1,ierr);CHKERRQ(ierr)
                  rhsPtr(dof) = tmpPtr1(1)
                  Call SectionRealRestore(boundaryTemperaturePreviousSec,setIdx(dof),tmpPtr1,ierr);CHKERRQ(ierr)
               End Do
            Else If (t == MEF90HeatXferCtx%timeTarget) Then
               Do dof = 1,size(setIdx)
                  Call SectionRealRestrict(boundaryTemperatureTargetSec,setIdx(dof),tmpPtr1,ierr);CHKERRQ(ierr)
                  rhsPtr(dof) = tmpPtr1(1)
                  Call SectionRealRestore(boundaryTemperatureTargetSec,setIdx(dof),tmpPtr1,ierr);CHKERRQ(ierr)
               End Do
            Else  !!! tmpPtr1 between the steps
               Do dof = 1,size(setIdx)
                  Call SectionRealRestrict(boundaryTemperaturePreviousSec,setIdx(dof),tmpPtr1,ierr);CHKERRQ(ierr)
                  Call SectionRealRestrict(boundaryTemperatureTargetSec,setIdx(dof),tmpPtr2,ierr);CHKERRQ(ierr)
                  rhsPtr(dof) = ((MEF90HeatXferCtx%timeTarget - t) * tmpPtr1(1) + (t - MEF90HeatXferCtx%timePrevious) * tmpPtr2(2)) &
                              /  (MEF90HeatXferCtx%timeTarget - MEF90HeatXferCtx%timePrevious)
                  Call SectionRealRestore(boundaryTemperatureTargetSec,setIdx(dof),tmpPtr2,ierr);CHKERRQ(ierr)
                  Call SectionRealRestore(boundaryTemperaturePreviousSec,setIdx(dof),tmpPtr1,ierr);CHKERRQ(ierr)
               End Do
            End If
            Call VecSetValues(rhs,size(setIdx),setdofIdx,rhsPtr,INSERT_VALUES,ierr);CHKERRQ(ierr)
            DeAllocate(rhsPtr)
            Call ISRestoreIndicesF90(setISdof,setdofIdx,ierr);CHKERRQ(ierr)
            Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
            Call ISRestoreIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If
      End Do

      !!!
      !!! Cleanup
      !!!
      Call VecAssemblyBegin(rhs,ierr)
      Call VecAssemblyEnd(rhs,ierr)

      Call SectionRealDestroy(rhsSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(boundaryTemperaturePreviousSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(boundaryTemperatureTargetSec,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90_APPEND(MEF90HeatXferRHS,MEF90_DIM)D

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXFerEnergy"
!!!
!!!  
!!!  MEF90HeatXFerEnergy:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90_APPEND(MEF90HeatXFerEnergy,MEF90_DIM)D(temperatureVec,t,MEF90HeatXferCtx,energy,work,ierr)
      Type(Vec),Intent(IN)                            :: temperatureVec
      PetscReal,Intent(IN)                            :: t
      Type(MEF90HeatXferCtx_Type),Intent(IN)          :: MEF90HeatXferCtx
      PetscReal,Dimension(:),Pointer                  :: energy,work
      PetscErrorCode,Intent(OUT)                      :: ierr
   
      Type(Vec)                                       :: fluxVec,externalTemperatureVec
      Type(SectionReal)                               :: temperatureSec,fluxSec,externalTemperatureSec
      Type(IS)                                        :: CellSetGlobalIS,setIS,setISdof
      PetscInt,dimension(:),Pointer                   :: setID,SetIdx
      PetscInt                                        :: set,QuadratureOrder,numCell
      Type(MEF90HeatXferCellSetOptions_Type),pointer  :: cellSetOptions
      Type(MEF90_MATPROP),pointer                     :: matpropSet
      Type(VecScatter)                                :: ScatterSecToVec,ScatterSecToVecCell
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer   :: elem
      Type(MEF90Element_Type)                         :: elemType
      Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90CtxGlobalOptions
      Type(MEF90HeatXferGlobalOptions_Type),pointer   :: MEF90HeatXferGlobalOptions
      PetscReal                                       :: myEnergy,myWork
      
      
      energy = 0.0_Kr
      work   = 0.0_Kr
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90HeatXferCtx%MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr);CHKERRQ(ierr)
   
      Call DMMeshGetSectionReal(MEF90HeatXferCtx%DM,'default',temperatureSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90HeatXferCtx%DM,temperatureSec,ScatterSecToVec,ierr);CHKERRQ(ierr)
      Call DMMeshGetSectionReal(MEF90HeatXferCtx%cellDM,'default',fluxSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(fluxSec,externalTemperatureSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90HeatXferCtx%CellDM,fluxSec,ScatterSecToVecCell,ierr);CHKERRQ(ierr)

      !!! Convert the fields to Section, interpolating the external temperature and fluxes if the current time
      !!! is neither the starting nor the target time
      Call SectionRealToVec(temperatureSec,ScatterSecToVec,SCATTER_REVERSE,temperatureVec,ierr);CHKERRQ(ierr)
      If (t == MEF90HeatXferCtx%timePrevious) Then
         Call SectionRealToVec(externalTemperatureSec,ScatterSecToVecCell,SCATTER_REVERSE,MEF90HeatXferCtx%externalTemperaturePrevious,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(fluxSec,ScatterSecToVecCell,SCATTER_REVERSE,MEF90HeatXferCtx%fluxPrevious,ierr);CHKERRQ(ierr)   
      Else If (t == MEF90HeatXferCtx%timeTarget) Then
         Call SectionRealToVec(externalTemperatureSec,ScatterSecToVecCell,SCATTER_REVERSE,MEF90HeatXferCtx%externalTemperatureTarget,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(fluxSec,ScatterSecToVecCell,SCATTER_REVERSE,MEF90HeatXferCtx%fluxTarget,ierr);CHKERRQ(ierr)   
      Else
         Call VecCopy(MEF90HeatXferCtx%externalTemperaturePrevious,externalTemperatureVec,ierr);CHKERRQ(ierr)
         Call VecScale(externalTemperatureVec,(MEF90HeatXferCtx%timeTarget - t) / (MEF90HeatXferCtx%timeTarget - MEF90HeatXferCtx%timePrevious),ierr);CHKERRQ(ierr)
         Call VecAXPY(externalTemperatureVec,(t - MEF90HeatXferCtx%timePrevious) / (MEF90HeatXferCtx%timeTarget - MEF90HeatXferCtx%timePrevious),MEF90HeatXferCtx%externalTemperatureTarget,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(externalTemperatureSec,ScatterSecToVecCell,SCATTER_REVERSE,externalTemperatureVec,ierr);CHKERRQ(ierr)

         Call VecCopy(MEF90HeatXferCtx%fluxPrevious,fluxVec,ierr);CHKERRQ(ierr)
         Call VecScale(fluxVec,(MEF90HeatXferCtx%timeTarget - t) / (MEF90HeatXferCtx%timeTarget - MEF90HeatXferCtx%timePrevious),ierr);CHKERRQ(ierr)
         Call VecAXPY(fluxVec,(t - MEF90HeatXferCtx%timePrevious) / (MEF90HeatXferCtx%timeTarget - MEF90HeatXferCtx%timePrevious),MEF90HeatXferCtx%fluxTarget,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(fluxSec,ScatterSecToVecCell,SCATTER_REVERSE,fluxVec,ierr);CHKERRQ(ierr)
      End If
      
      Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         myenergy = 0.0_Kr
         mywork = 0.0_Kr
         Call DMMeshGetStratumIS(MEF90HeatXferCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call PetscBagGetDataMEF90MatProp(MEF90HeatXferCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         
         elemType = MEF90_knownElements(cellSetOptions%ElemTypeShortID)
         QuadratureOrder = elemType%order * 2
         Call MEF90Element_Create(MEF90HeatXferCtx%DM,setIS,elem,QuadratureOrder,CellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)   
         
         Call MEF90DiffusionEnergySet(myenergy,temperatureSec,externalTemperatureSec,MEF90HeatXferCtx%DM,matpropSet%ThermalConductivity,cellSetOptions%SurfaceThermalConductivity,setIS,elem,elemType,ierr);CHKERRQ(ierr)
         Call MPI_AllReduce(myEnergy,energy(set),1,MPIU_SCALAR,MPI_SUM,MEF90HeatXferCtx%MEF90Ctx%comm,ierr);CHKERRQ(ierr)
         
         Call MEF90DiffusionWorkSetCell(mywork,temperatureSec,MEF90HeatXferCtx%DM,fluxSec,setIS,elem,elemType,ierr);CHKERRQ(ierr)
         Call MPI_AllReduce(myWork,work(set),1,MPIU_SCALAR,MPI_SUM,MEF90HeatXferCtx%MEF90Ctx%comm,ierr);CHKERRQ(ierr)
         Call MEF90Element_Destroy(elem,ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
      
      Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(externalTemperatureSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(fluxSec,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVecCell,ierr);CHKERRQ(ierr)
   End Subroutine MEF90_APPEND(MEF90HeatXFerEnergy,MEF90_DIM)D
End Module MEF90_APPEND(m_MEF90_HeatXferAssembly,MEF90_DIM)D
