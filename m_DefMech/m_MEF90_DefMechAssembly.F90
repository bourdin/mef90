#include "../MEF90/mef90.inc"
Module MEF90_APPEND(m_MEF90_DefMechAssembly,MEF90_DIM)D
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_DefMechCtx
   Implicit none
   Private
   Public MEF90DefMechOperator,     &
          MEF90DefMechBilinearForm
!          MEF90DefMechEnergy

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperator"
!!!
!!!  
!!!  MEF90DefMechOperator: Build the operator. When called in SNES, the solution time should always match the target time, 
!!!                         so there is no need for interpolation of the fluxes, external, and boundary values
!!!  
!!!  (c) 2012-13 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechOperator(snesDisplacement,x,residual,MEF90DefMechCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesDisplacement
      Type(Vec),Intent(IN)                               :: x
      Type(Vec),Intent(INOUT)                            :: residual
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(DM)                                           :: mesh
      Type(SectionReal)                                  :: xSec,residualSec
      Type(SectionReal)                                  :: inelasticStrainSec,inelasticStrainCellSec
      Type(SectionReal)                                  :: boundaryDisplacementSec,forceSec,pressureForceSec
      PetscReal,Dimension(:),Pointer                     :: xPtr,residualPtr
      PetscReal,Dimension(:),Pointer                     :: inelasticStrainPtr,inelasticStrainCellPtr
      PetscReal,Dimension(:),Pointer                     :: boundaryDisplacementPtr,forcePtr,pressureForcePtr
      Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof,bcIS
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),pointer      :: cellSetOptions
      Type(MEF90DefMechVertexSetOptions_Type),pointer    :: vertexSetOptions
      Type(VecScatter)                                   :: ScatterSecToVec,ScatterSecToVecScal,ScatterSecToVecMatS
      Type(VecScatter)                                   :: ScatterSecToVecCell,ScatterSecToVecCellScal,ScatterSecToVecCellMatS
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemDamage
      Type(MEF90Element_Type)                            :: elemDisplacementType,elemDamageType
      PetscInt                                           :: p,dof,c,dim,numDof
      PetscInt                                           :: nVal

      Call SNESGetDM(snesDisplacement,mesh,ierr);CHKERRQ(ierr)
      Call DMMeshGetDimension(mesh,dim,ierr);CHKERRQ(ierr)
      !!! Create dof-based sections
      Call DMMeshGetSectionReal(MEF90DefMechCtx%DMVect,'default',xSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DM,xSec,ScatterSecToVec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(xSec,residualSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(xSec,boundaryDisplacementSec,ierr);CHKERRQ(ierr)
      !Call DMMeshGetSectionReal(MEF90DefMechCtx%DMMatS,'default',inelasticStrainSec,ierr);CHKERRQ(ierr)
      !Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMMatS,inelasticStrainSec,ScatterSecToVecMatS,ierr);CHKERRQ(ierr)

      !!! Create cell based sections, and allocate required pointers
      !!! I could check if any of these vectors are null 
      Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMVect,'default',forceSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMVect,forceSec,ScatterSecToVecCell,ierr);CHKERRQ(ierr)
      Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMScal,'default',pressureForceSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMScal,pressureForceSec,ScatterSecToVecCellScal,ierr);CHKERRQ(ierr)
      Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMMatS,'default',inelasticStrainCellSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMMatS,inelasticStrainCellSec,ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)
      Allocate(pressureForcePtr(1))
      Allocate(forcePtr(dim))
      Allocate(inelasticStrainCellPtr(dim*(dim+1)/2))
   
      !!! Scatter data from Vec to Sec, or initialize
      Call SectionRealToVec(xSec,ScatterSecToVec,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr) 
        
      Call SectionRealToVec(forceSec,ScatterSecToVecCell,SCATTER_REVERSE,MEF90DefMechCtx%force,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(pressureForceSec,ScatterSecToVecCellScal,SCATTER_REVERSE,MEF90DefMechCtx%pressureForce,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(boundaryDisplacementSec,ScatterSecToVec,SCATTER_REVERSE,MEF90DefMechCtx%boundaryDisplacement,ierr);CHKERRQ(ierr)

      Call SectionRealSet(residualSec,0.0_Kr,ierr);CHKERRQ(ierr)
      Call VecSet(residual,0.0_Kr,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

      !!!
      !!! We loop over all element twice. The first time in order to assembly all non BC cell sets
      !!! In the second pass, we only update the BC where necessary
      !!! vertex set BC are updated last, so that they override cell set BC
      !!!
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         elemDisplacementType = MEF90_knownElements(cellSetOptions%elemTypeShortIDDisplacement)
         elemDamageType = MEF90_knownElements(cellSetOptions%elemTypeShortIDDamage)

         !!! Call proper local assembly depending on the type of damage law
         Select Case (cellSetOptions%defectLaw)
         Case (MEF90DefMech_defectLawElasticity)
            QuadratureOrder = 2 * elemDisplacementType%order
            Call MEF90Element_Create(mesh,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
            !!! Contribution of the bilinear form
            Call MEF90ElasticityOperatorSet(residualSec,mesh,xSec,setIS,matPropSet%HookesLaw,elemDisplacement,elemDisplacementType,ierr);CHKERRQ(ierr)
            !!! Inelastic strains

            !!! Force
            Call MEF90ElasticityForceRHSSetCell(residualSec,mesh,forceSec,setIS,elemDisplacement,elemDisplacementType,ierr);CHKERRQ(ierr)

            !!! pressure Force
            
            Call MEF90Element_Destroy(elemDisplacement,ierr)
         Case (MEF90DefMech_defectLawBrittleFracture,MEF90DefMech_defectLawPlasticity)
            Print*,__FUNCT__,': Unimplemented Damage law',cellSetOptions%defectLaw
            STOP      
         End Select
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do   
      DeAllocate(pressureForcePtr)
      DeAllocate(forcePtr)  
      DeAllocate(inelasticStrainCellPtr)

      !!! "Ghost update" for the residual Section
      Call SectionRealComplete(residualSec,ierr);CHKERRQ(ierr)
      !!! Scatter back from SectionReal to Vec
      Call SectionRealToVec(residualSec,ScatterSecToVec,SCATTER_FORWARD,residual,ierr);CHKERRQ(ierr)

      !!!
      !!! Cell set BC
      !!!
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Do c = 1, dim
            If (cellSetOptions%Has_DisplacementBC(c)) Then
               Call DMMeshGetStratumIS(MEF90DefMechCtx%DMVect,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
               Call MEF90_ISCreateCelltoVertex(MEF90DefMechCtx%DMVect,PETSC_COMM_WORLD,setIS,bcIS,ierr)
               Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
               Allocate(xPtr(nval),stat=ierr)
               Allocate(residualPtr(nval),stat=ierr)
               Allocate(boundaryDisplacementPtr(nval),stat=ierr)
               Call MEF90_VecGetValuesISdof(MEF90DefMechCtx%DMVect,x,xPtr,bcIS,c,ierr)
               Call MEF90_VecGetValuesISdof(MEF90DefMechCtx%DMVect,MEF90DefMechCtx%boundaryDisplacement,boundaryDisplacementPtr,bcIS,c,ierr)
               residualPtr = xPtr - boundaryDisplacementPtr
               Call MEF90_VecSetValuesISdof(MEF90DefMechCtx%DMVect,residual,residualPtr,bcIS,c,INSERT_VALUES,ierr)
               DeAllocate(boundaryDisplacementPtr)
               DeAllocate(residualPtr)
               DeAllocate(xPtr)
               Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
               Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
            End If ! cellSetOptions%Has_BC
         End Do
      End Do ! set
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)

      
      !!!
      !!! Vertex set BC
      !!!
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxVertexSetOptions(MEF90DefMechCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
         Do c = 1, dim
            If (vertexSetOptions%Has_DisplacementBC(c)) Then
               Call DMMeshGetStratumIS(MEF90DefMechCtx%DMVect,'Vertex Sets',setID(set),bcIS,ierr);CHKERRQ(iErr)
               Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
               Allocate(xPtr(nval),stat=ierr)
               Allocate(residualPtr(nval),stat=ierr)
               Allocate(boundaryDisplacementPtr(nval),stat=ierr)
               Call MEF90_VecGetValuesISdof(MEF90DefMechCtx%DMVect,x,xPtr,bcIS,c,ierr)
               Call MEF90_VecGetValuesISdof(MEF90DefMechCtx%DMVect,MEF90DefMechCtx%boundaryDisplacement,boundaryDisplacementPtr,bcIS,c,ierr)
               residualPtr = xPtr - boundaryDisplacementPtr
               Call MEF90_VecSetValuesISdof(MEF90DefMechCtx%DMVect,residual,residualPtr,bcIS,c,INSERT_VALUES,ierr)
               DeAllocate(boundaryDisplacementPtr)
               DeAllocate(residualPtr)
               DeAllocate(xPtr)
               Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
            End If ! vertexSetOptions%Has_BC
         End Do
      End Do ! set
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
   
      Call VecAssemblyBegin(residual,ierr)
      Call VecAssemblyEnd(residual,ierr)

      Call SectionRealDestroy(inelasticStrainCellSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(pressureForceSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(forceSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(boundaryDisplacementSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(residualSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)

      Call VecScatterDestroy(ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)      
      Call VecScatterDestroy(ScatterSecToVecCellScal,ierr);CHKERRQ(ierr)      
      Call VecScatterDestroy(ScatterSecToVecCell,ierr);CHKERRQ(ierr)      
      Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
      
      !!! Check destructions
   End Subroutine MEF90DefMechOperator

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearForm"
!!!
!!!  
!!!  MEF90DefMechBilinearForm:
!!!  
!!!  (c) 2012-13 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechBilinearForm(snesDispl,x,A,M,flg,MEF90DefMechCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesDispl
      Type(Vec),Intent(IN)                               :: x
      Type(Mat),Intent(INOUT)                            :: A,M
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr  
         
      Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),pointer      :: cellSetOptions
      Type(MEF90DefMechVertexSetOptions_Type),pointer    :: vertexSetOptions
      Type(DM)                                           :: mesh
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemDamage
      Type(MEF90Element_Type)                            :: elemTypeDisplacement,elemTypeDamage
      PetscInt,Dimension(:),Pointer                      :: setIdx,bcIdx,Cone
      Type(IS)                                           :: bcIS
      PetscInt                                           :: cell,v,numBC,numDoF,numCell,c,dim
      
      Call MatZeroEntries(A,ierr);CHKERRQ(ierr)
      Call SNESGetDM(snesDispl,mesh,ierr);CHKERRQ(ierr)

      Call DMMeshGetDimension(mesh,dim,ierr);CHKERRQ(ierr)
      Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         elemTypeDisplacement = MEF90_knownElements(cellSetOptions%elemTypeShortIDDisplacement)
         elemTypeDamage = MEF90_knownElements(cellSetOptions%elemTypeShortIDDamage)

         !!! Call proper local assembly depending on the type of damage law
         Select Case (cellSetOptions%defectLaw)
         Case (MEF90DefMech_defectLawElasticity)
            QuadratureOrder = 2 * elemTypeDisplacement%order
            Call MEF90Element_Create(mesh,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
            Call MEF90ElasticityBilinearFormSet(A,mesh,setIS,matPropSet%HookesLaw,elemDisplacement,elemTypeDisplacement,ierr);CHKERRQ(ierr)
            Call MEF90Element_Destroy(elemDisplacement,ierr)
         Case (MEF90DefMech_defectLawBrittleFracture,MEF90DefMech_defectLawPlasticity)
            Print*,__FUNCT__,': Unimplemented Damage law',cellSetOptions%defectLaw
            STOP      
         End Select
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do     
      Call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
      Call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
      !!!
      !!! Boundary conditions at cell sets
      !!!
      Call DMmeshGetStratumSize(mesh,"height",0,numCell,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call MEF90_ISCreateCelltoVertex(mesh,PETSC_COMM_WORLD,setIS,bcIS,ierr)
         Do c = 1, dim
            If (cellSetOptions%Has_displacementBC(c)) Then
               Call DMMeshISCreateISglobaldof(mesh,bcIS,c-1,setISdof,ierr);CHKERRQ(ierr)
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
         Call PetscBagGetDataMEF90DefMechCtxVertexSetOptions(MEF90DefMechCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
         Do c = 1, dim
            If (vertexSetOptions%Has_displacementBC(c)) Then
               Call DMMeshGetStratumIS(mesh,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
               Call DMMeshISCreateISglobaldof(mesh,setIS,c-1,setISdof,ierr);CHKERRQ(ierr)
               Call MatZeroRowsColumnsIS(A,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
            End If
         End Do
      End Do
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      
      flg = SAME_NONZERO_PATTERN
   End Subroutine MEF90DefMechBilinearForm
End Module MEF90_APPEND(m_MEF90_DefMechAssembly,MEF90_DIM)D
