#include "../MEF90/mef90.inc"
Module MEF90_APPEND(m_MEF90_ElasticityAssembly,MEF90_DIM)D
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_ElasticityCtx
   Implicit none
   Private
   Public MEF90ElasticityOperator,     &
          MEF90ElasticityBilinearForm
!          MEF90ElasticityEnergy

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90ElasticityOperator"
!!!
!!!  
!!!  MEF90ElasticityOperator: Build the operator. When called in SNES, the solution time should always match the target time, 
!!!                         so there is no need for interpolation of the fluxes, external, and boundary values
!!!  
!!!  (c) 2012-13 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90ElasticityOperator(snesTemp,x,residual,MEF90ElasticityCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesTemp
      Type(Vec),Intent(IN)                               :: x
      Type(Vec),Intent(INOUT)                            :: residual
      Type(MEF90ElasticityCtx_Type),Intent(IN)           :: MEF90ElasticityCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(DM)                                           :: mesh
      Type(SectionReal)                                  :: xSec,residualSec
      Type(SectionReal)                                  :: inealsticStrain,inelasticStrainCell
      Type(SectionReal)                                  :: boundaryDisplacementSec,forceSec,pressureForceSec
      PetscReal,Dimension(:),Pointer                     :: xPtr,residualPtr
      PetscReal,Dimension(:),Pointer                     :: inelasticStrainPtr,inelasticStrainCellPtr
      PetscReal,Dimension(:),Pointer                     :: boundaryDisplacementPtr,forcePtr,boundaryForcePtr
      Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),pointer                        :: matpropSet
      Type(MEF90ElasticityCellSetOptions_Type),pointer   :: cellSetOptions
      Type(MEF90ElasticityVertexSetOptions_Type),pointer :: vertexSetOptions
      Type(VecScatter)                                   :: ScatterSecToVec,ScatterSecToVecScal,ScatterSecToVecMatS
      Type(VecScatter)                                   :: ScatterSecToVecCell,ScatterSecToVecCellScal,ScatterSecToVecCellMatS
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elem
      Type(MEF90Element_Type)                            :: elemType
      PetscInt                                           :: cell,dof,c,dim

      Call DMMeshGetDimension(mesh,dim,ierr);CHKERRQ(ierr)
      !!! Create dof-based sections
      Call DMMeshGetSectionReal(MEF90ElasticityCtx%DMVect,'default',xSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90ElasticityCtx%DM,xSec,ScatterSecToVec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(xSec,residualSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(xSec,boundaryDisplacementSec,ierr);CHKERRQ(ierr)
      !Call DMMeshGetSectionReal(MEF90ElasticityCtx%DMMatS,'default',inelasticStrainSec,ierr);CHKERRQ(ierr)
      !Call DMMeshCreateGlobalScatter(MEF90ElasticityCtx%DMMatS,inelasticStrainSec,ScatterSecToVecMatS,ierr);CHKERRQ(ierr)

      !!! Create cell based sections, and allocate required pointers
      Call DMMeshGetSectionReal(MEF90ElasticityCtx%CellDMVect,'default',forceSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90ElasticityCtx%CellDMVect,forceSec,ScatterSecToVecCell,ierr);CHKERRQ(ierr)
      Call DMMeshGetSectionReal(MEF90ElasticityCtx%CellDMScal,'default',pressureForceSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90ElasticityCtx%CellDMScal,pressureForceSec,ScatterSecToVecCellScal,ierr);CHKERRQ(ierr)
      !Call DMMeshGetSectionReal(MEF90ElasticityCtx%CellDMMatS,'default',inelasticStrainCellSec,ierr);CHKERRQ(ierr)
      !Call DMMeshCreateGlobalScatter(MEF90ElasticityCtx%CellDMMatS,inelasticStrainCellSec,ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)
      Allocate(pressureForcePtr(1))
      Allocate(forcePtr(dim))
      Allocate(inelasticStrainCell(dim*(dim+1)/2)
   
      !!! Scatter data from Vec to Sec, or initialize
      Call SectionRealToVec(xSec,ScatterSecToVec,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr)
        
      Call SectionRealToVec(forceSec,ScatterSecToVecCell,SCATTER_REVERSE,MEF90ElasticityCtx%temperature,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(pressureForceSec,ScatterSecToVecCellScal,SCATTER_REVERSE,MEF90ElasticityCtx%externalTemperature,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(boundaryDisplacementSec,ScatterSecToVec,SCATTER_REVERSE,MEF90ElasticityCtx%boundaryTemperature,ierr);CHKERRQ(ierr)

      Call SectionRealSet(residualSec,0.0_Kr,ierr);CHKERRQ(ierr)
      Call VecSet(residual,0.0_Kr,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(MEF90ElasticityCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

      !!!
      !!! We loop over all element twice. The first time in order to assembly all not BC cell sets
      !!! In the second pass, we only update the BC where necessary
      !!! vertex cet BC are updated last, so that they override cell set BC
      !!!
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90ElasticityCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90ElasticityCtxCellSetOptions(MEF90ElasticityCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         
         If (.not. cellSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(MEF90ElasticityCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
            Call ISGetIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
         
            elemType = MEF90_knownElements(cellSetOptions%ElemTypeShortID)
            QuadratureOrder = elemType%order * 2
            Call MEF90Element_Create(MEF90ElasticityCtx%DM,setIS,elem,QuadratureOrder,CellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)
            Call MEF90DiffusionOperatorSet(residualSec,MEF90ElasticityCtx%DM,xSec,setIS,matpropSet%ThermalConductivity,cellSetOptions%SurfaceThermalConductivity,elem,elemType,ierr);CHKERRQ(ierr)

            !!! Modified flux is flux + surfaceThermalConductivity * refTemp      
            !!! I _could_ use a SecAXPY, but this would summ all values at all cells for each block
            !!! I _could_ also create Sections restricted to the cell set only
            !!! But this is going away with sections anyway...
            Do cell = 1,size(setIdx)
               Call SectionRealRestrict(fluxSec,setIdx(cell),fluxPtr,ierr);CHKERRQ(ierr)
               modifiedFluxPtr = - fluxPtr
               Call SectionRealUpdate(modifiedFluxSec,setIdx(cell),modifiedFluxPtr,ADD_VALUES,ierr);CHKERRQ(ierr)
               Call SectionRealRestore(fluxSec,setIdx(cell),fluxPtr,ierr);CHKERRQ(ierr)
            End Do
            Do cell = 1,size(setIdx)
               Call SectionRealRestrict(externalTemperatureSec,setIdx(cell),externalTemperaturePtr,ierr);CHKERRQ(ierr)
               modifiedFluxPtr = - cellSetOptions%SurfaceThermalConductivity * externalTemperaturePtr
               Call SectionRealUpdate(modifiedFluxSec,setIdx(cell),modifiedFluxPtr,ADD_VALUES,ierr);CHKERRQ(ierr)
               Call SectionRealRestore(externalTemperatureSec,setIdx(cell),externalTemperaturePtr,ierr);CHKERRQ(ierr)
            End Do
            Call MEF90DiffusionRHSSetCell(residualsec,MEF90ElasticityCtx%DM,modifiedFluxSec,setIS,elem,elemType,ierr);CHKERRQ(ierr)
          
            Call MEF90Element_Destroy(elem,ierr);CHKERRQ(ierr)
            Call ISRestoreIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If ! cellSetOptions%Has_BC
      End Do ! set
      
      Call SNESGetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90ElasticityCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90ElasticityCtxCellSetOptions(MEF90ElasticityCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         
         If (cellSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(MEF90ElasticityCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
            Call ISGetIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
            elemType = MEF90_knownElements(cellSetOptions%ElemTypeShortID)
            Allocate(boundaryTemperaturePtr(elemType%numDof))
            Allocate(xPtr(elemType%numDof))
            Allocate(residualPtr(elemType%numDof))
            Do cell = 1,size(setIdx)
               Call SectionRealRestrictClosure(boundaryTemperatureSec,mesh,setIdx(cell),elemType%numDof,boundaryTemperaturePtr,ierr);CHKERRQ(ierr)
               Call SectionRealRestrictClosure(xSec,mesh,setIdx(cell),elemType%numDof,xPtr,ierr);CHKERRQ(ierr)
               residualPtr = xPtr - BoundaryTemperaturePtr
               Call SectionRealUpdateClosure(residualSec,mesh,setIdx(cell),boundaryTemperaturePtr,INSERT_VALUES,ierr);CHKERRQ(iErr)
            End Do
            Call ISRestoreIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
            DeAllocate(boundaryTemperaturePtr)
            DeAllocate(xPtr)
            DeAllocate(residualPtr)
         End If ! cellSetOptions%Has_BC
      End Do ! set
      
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
      !!! "Ghost update" for the residual Section
      Call SectionRealComplete(residualSec,ierr);CHKERRQ(ierr)
      !!! Scatter back from SectionReal to Vec
      Call SectionRealToVec(residualSec,ScatterSecToVec,SCATTER_FORWARD,residual,ierr);CHKERRQ(ierr)
      
      
      !!!
      !!! Account for BC entries in the residual
      !!!
      Call DMmeshGetLabelIdIS(MEF90ElasticityCtx%DM,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90ElasticityCtxVertexSetOptions(MEF90ElasticityCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
         If (vertexSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(MEF90ElasticityCtx%DM,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
            Call ISGetIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
            Call DMMeshISCreateISglobaldof(MEF90ElasticityCtx%DM,setIS,0,setISdof,ierr);CHKERRQ(ierr)
            Call ISGetIndicesF90(setISdof,setdofIdx,ierr);CHKERRQ(ierr)
            Allocate(residualPtr(size(setIdx)))
            Do dof = 1, size(setIdx)
               Call SectionRealRestrict(xSec,setIdx(dof),xPtr,ierr);CHKERRQ(ierr)
               Call SectionRealRestrict(boundaryTemperatureSec,setIdx(dof),boundaryTemperaturePtr,ierr);CHKERRQ(ierr)
               residualPtr(dof) = xPtr(1) - boundaryTemperaturePtr(1)
               Call SectionRealRestore(xSec,setIdx(dof),xPtr,ierr);CHKERRQ(ierr)
               Call SectionRealRestore(boundaryTemperatureSec,setIdx(dof),boundaryTemperaturePtr,ierr);CHKERRQ(ierr)
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
      Call SectionRealDestroy(modifiedfluxSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(fluxSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(externalTemperatureSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(boundaryTemperatureSec,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVecCell,ierr);CHKERRQ(ierr)      
   End Subroutine MEF90ElasticityOperator
   

#undef __FUNCT__
#define __FUNCT__ "MEF90ElasticityBilinearForm"
!!!
!!!  
!!!  MEF90ElasticityBilinearForm:
!!!  
!!!  (c) 2012-13 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90ElasticityBilinearForm(snesTemp,x,A,M,flg,MEF90ElasticityCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesTemp
      Type(Vec),Intent(IN)                               :: x
      Type(Mat),Intent(INOUT)                            :: A,M
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90ElasticityCtx_Type),Intent(IN)           :: MEF90ElasticityCtx
      PetscErrorCode,Intent(OUT)                         :: ierr  
         
      Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),pointer                        :: matpropSet
      Type(MEF90ElasticityCellSetOptions_Type),pointer   :: cellSetOptions
      Type(MEF90ElasticityVertexSetOptions_Type),pointer :: vertexSetOptions
      Type(DM)                                           :: mesh
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elem
      Type(MEF90Element_Type)                            :: elemType
      PetscInt,Dimension(:),Pointer                      :: setIdx,bcIdx,Cone
      Type(IS)                                           :: bcIS
      PetscInt                                           :: cell,v,numBC,numDoF,numCell,c,dim
      
      Call MatZeroEntries(A,ierr);CHKERRQ(ierr)
      Call SNESGetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)

      Call DMMeshGetDimension(mesh,dim,ierr);CHKERRQ(ierr)
      Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90ElasticityCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90ElasticityCtxCellSetOptions(MEF90ElasticityCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         elemType = MEF90_knownElements(cellSetOptions%ElemTypeShortID)
         QuadratureOrder = 2 * elemType%order
         Call MEF90Element_Create(mesh,setIS,elem,QuadratureOrder,CellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)
         Call MEF90ElasticityBilinearFormSet(A,mesh,setIS,matPropSet%HookesLaw,elem,elemType,ierr);CHKERRQ(ierr)
         Call MEF90Element_Destroy(elem,ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do     
      Call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
      Call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
      !!!
      !!! Boundary conditions at cell sets
      !!!
      Call DMmeshGetStratumSize(mesh,"height",0,numCell,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90ElasticityCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90ElasticityCtxCellSetOptions(MEF90ElasticityCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call MEF90_ISCreateCelltoVertex(mesh,PETSC_COMM_WORLD,setIS,bcIS,ierr)
         Do c = 1, dim
            If (cellSetOptions%Has_displacementBC(c)) Then
               Call DMMeshISCreateISglobaldof(mesh,setIS,c,setISdof,ierr);CHKERRQ(ierr)
               Call MatZeroRowsColumnsIS(A,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
               Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
            End If ! cellSetOptions%Has_displacementBC
         End Do ! c
         Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do     
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
   
      !!!
      !!! Boundary conditions at vertex sets
      !!!
      Call DMmeshGetLabelIdIS(mesh,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90ElasticityCtxVertexSetOptions(MEF90ElasticityCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
!!! Adapt for BC / dof
         Do c = 1, dim
            If (vertexSetOptions%Has_displacementBC(c)) Then
               Call DMMeshGetStratumIS(mesh,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
               Call DMMeshISCreateISglobaldof(mesh,setIS,c,setISdof,ierr);CHKERRQ(ierr)
               Call MatZeroRowsColumnsIS(A,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
            End If
         End Do
      End Do
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      
      flg = SAME_NONZERO_PATTERN
   End Subroutine MEF90ElasticityBilinearForm
End Module MEF90_APPEND(m_MEF90_ElasticityAssembly,MEF90_DIM)D
