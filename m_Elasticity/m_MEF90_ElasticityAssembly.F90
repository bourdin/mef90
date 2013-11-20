#include "../MEF90/mef90.inc"
Module MEF90_APPEND(m_MEF90_ElasticityAssembly,MEF90_DIM)D
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_ElasticityCtx
   Implicit none
   Private
   Public MEF90ElasticityBilinearForm!MEF90ElasticityOperator!, &
!          MEF90ElasticityBilinearForm, &
!          MEF90ElasticityEnergy

Contains
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
      PetscInt                                           :: cell,v,numBC,numDoF,numCell
      
      Call MatZeroEntries(A,ierr);CHKERRQ(ierr)
      Call SNESGetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90ElasticityCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90ElasticityCtxCellSetOptions(MEF90ElasticityCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
!!! Adapt for BC / dof
         !If (.not. cellSetOptions%Has_displacementBC) Then
         !   Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         !   elemType = MEF90_knownElements(cellSetOptions%ElemTypeShortID)
         !   QuadratureOrder = 2 * elemType%order
         !   Call MEF90Element_Create(mesh,setIS,elem,QuadratureOrder,CellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)
         !   Call MEF90DiffusionBilinearFormSet(A,mesh,setIS,matpropSet%ThermalConductivity,cellSetOptions%SurfaceThermalConductivity,elem,elemType,ierr);CHKERRQ(ierr)
         !   Call MEF90Element_Destroy(elem,ierr)
         !   Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         !End If ! cellSetOptions%Has_displacementBC
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
!!! Adapt for BC / dof
         !If (cellSetOptions%Has_displacementBC) Then
         !   Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         !   Call MEF90_ISCreateCelltoVertex(mesh,PETSC_COMM_WORLD,setIS,bcIS,ierr)
         !   Call DMMeshISCreateISglobaldof(MEF90ElasticityCtx%DM,setIS,0,setISdof,ierr);CHKERRQ(ierr)
         !   Call MatZeroRowsColumnsIS(A,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
         !   Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
         !   Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         !   Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
         !End If ! cellSetOptions%Has_displacementBC
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
         !If (vertexSetOptions%Has_displacementBC(:)) Then
         !   Call DMMeshGetStratumIS(mesh,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         !   Call DMMeshISCreateISglobaldof(mesh,setIS,0,setISdof,ierr);CHKERRQ(ierr)
         !   Call MatZeroRowsColumnsIS(A,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
         !End If
      End Do
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      
      flg = SAME_NONZERO_PATTERN
   End Subroutine MEF90ElasticityBilinearForm
End Module MEF90_APPEND(m_MEF90_ElasticityAssembly,MEF90_DIM)D
