#include "../MEF90/mef90.inc"
Module MEF90_APPEND(m_MEF90_DefMechAssembly,MEF90_DIM)D
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_DefMechCtx
   Implicit none
   Private
   Public MEF90DefMechOperatorDisplacement,     &
          MEF90DefMechBilinearFormDisplacement, &
          MEF90DefMechWork,                     &
          MEF90DefMechElasticEnergy,            &
          MEF90DefMechOperatorDamage,           &
          MEF90DefMechBilinearFormDamage,       &
          MEF90DefMechSurfaceEnergy

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDisplacement"
!!!
!!!  
!!!  MEF90DefMechOperatorDisplacement: Build the operator. When called in SNES, the solution time should always match the target time, 
!!!                         so there is no need for interpolation of the forcees, external, and boundary values
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechOperatorDisplacement(snesDisplacement,x,residual,MEF90DefMechCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesDisplacement
      Type(Vec),Intent(IN)                               :: x
      Type(Vec),Intent(INOUT)                            :: residual
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(DM)                                           :: mesh
      Type(SectionReal)                                  :: xSec,residualSec
      Type(SectionReal)                                  :: plasticStrainSec,temperatureSec,damageSec
      Type(SectionReal)                                  :: boundaryDisplacementSec,forceSec,pressureForceSec
      PetscReal,Dimension(:),Pointer                     :: xPtr,residualPtr
      PetscReal,Dimension(:),Pointer                     :: boundaryDisplacementPtr
      Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof,bcIS
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),pointer      :: cellSetOptions
      Type(MEF90DefMechVertexSetOptions_Type),pointer    :: vertexSetOptions
      Type(VecScatter)                                   :: ScatterSecToVec,ScatterSecToVecScal
      Type(VecScatter)                                   :: ScatterSecToVecCell,ScatterSecToVecCellScal,ScatterSecToVecCellMatS
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemDamage
      Type(MEF90Element_Type)                            :: elemDisplacementType,elemDamageType
      PetscInt                                           :: p,dof,c,dim,numDof
      PetscInt                                           :: nVal

      Call SNESGetDM(snesDisplacement,mesh,ierr);CHKERRQ(ierr)
      Call DMMeshGetDimension(mesh,dim,ierr);CHKERRQ(ierr)
      !!! Create dof-based sections and Scatter from the Vecs
      Call DMMeshGetSectionReal(MEF90DefMechCtx%DMVect,'default',xSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMVect,xSec,ScatterSecToVec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(xSec,residualSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(xSec,boundaryDisplacementSec,ierr);CHKERRQ(ierr)
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call DMMeshGetSectionReal(MEF90DefMechCtx%DMScal,'default',temperatureSec,ierr);CHKERRQ(ierr)
         Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMScal,temperatureSec,ScatterSecToVecScal,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(temperatureSec,ScatterSecToVecScal,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr)   
      Else
         temperatureSec%v = 0
      End If

      !!! Create cell based sections, and allocate required pointers
      !!! I could check if any of these vectors are null 
      Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMVect,'default',forceSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMVect,forceSec,ScatterSecToVecCell,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(forceSec,ScatterSecToVecCell,SCATTER_REVERSE,MEF90DefMechCtx%force,ierr);CHKERRQ(ierr)
      
      If (Associated(MEF90DefMechCtx%pressureForce)) Then
         Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMScal,'default',pressureForceSec,ierr);CHKERRQ(ierr)
         Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMScal,pressureForceSec,ScatterSecToVecCellScal,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(pressureForceSec,ScatterSecToVecCellScal,SCATTER_REVERSE,MEF90DefMechCtx%pressureForce,ierr);CHKERRQ(ierr)
      Else
         pressureForceSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMMatS,'default',plasticStrainSec,ierr);CHKERRQ(ierr)
         Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMMatS,plasticStrainSec,ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticStrainSec,ScatterSecToVecCellMatS,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)          
      Else
         PlasticStrainSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%damage)) Then
         If (Associated(MEF90DefMechCtx%temperature)) Then
            Call SectionRealDuplicate(temperatureSec,damageSec,ierr);CHKERRQ(ierr)
         Else
            Call DMMeshGetSectionReal(MEF90DefMechCtx%DMScal,'default',damageSec,ierr);CHKERRQ(ierr)
            Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMScal,damageSec,ScatterSecToVecScal,ierr);CHKERRQ(ierr)
         End If
         Call SectionRealToVec(damageSec,ScatterSecToVecScal,SCATTER_REVERSE,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)
      Else
         damageSec%v = 0
      End If

      !!! Scatter data from Vec to Sec, or initialize
      Call SectionRealToVec(xSec,ScatterSecToVec,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr) 
      Call SectionRealToVec(boundaryDisplacementSec,ScatterSecToVec,SCATTER_REVERSE,MEF90DefMechCtx%boundaryDisplacement,ierr);CHKERRQ(ierr)

      Call SectionRealSet(residualSec,0.0_Kr,ierr);CHKERRQ(ierr)
      Call VecSet(residual,0.0_Kr,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
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
         elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
         elemDamageType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)

         !!! Call proper local assembly depending on the type of damage law
         Select Case (cellSetOptions%defectLaw)
         Case (MEF90DefMech_defectLawElasticity)
            If (Associated(MEF90DefMechCtx%temperature)) Then
               QuadratureOrder = 2 * max(elemDisplacementType%order - 1, elemDamageType%order)
            Else
               QuadratureOrder = 2 * elemDisplacementType%order - 2
            End If
            Call MEF90Element_Create(mesh,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
            Call MEF90Element_Create(mesh,setIS,elemDamage,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
            !!! Contribution of the bilinear form
            Call MEF90ElasticityOperatorSet(residualSec,mesh,xSec,setIS,matPropSet%HookesLaw,elemDisplacement,elemDisplacementType,ierr);CHKERRQ(ierr)
            !!! temperature
            If (Associated(MEF90DefMechCtx%temperature)) Then
               Call MEF90ElasticityInelasticStrainRHSSetVertex(residualSec,MEF90DefMechCtx%DM,MEF90DefMechCtx%DMScal,temperatureSec,matPropSet%HookesLaw*matPropSet%LinearThermalExpansion,setIS,elemDisplacement,elemDisplacementType,elemDamage,elemDamageType,ierr)
            End If
            !!! plastic Strain
            If (Associated(MEF90DefMechCtx%plasticStrain)) Then
               Call MEF90ElasticityInelasticStrainRHSSetCell(residualSec,MEF90DefMechCtx%DM,MEF90DefMechCtx%CellDMMatS,plasticStrainSec,matPropSet%HookesLaw,setIS,elemDisplacement,elemDisplacementType,ierr)
            End If
            !!! Force
            Call MEF90ElasticityForceRHSSetCell(residualSec,mesh,forceSec,setIS,elemDisplacement,elemDisplacementType,ierr);CHKERRQ(ierr)

            !!! pressure Force
            If (Associated(MEF90DefMechCtx%pressureForce) .AND. (elemDisplacementType%coDim > 0)) Then
               Call MEF90ElasticityPressureForceRHSSetCell(residualSec,mesh,pressureForceSec,setIS,elemDisplacement,elemDisplacementType,ierr);CHKERRQ(ierr)
            End If
            Call MEF90Element_Destroy(elemDisplacement,ierr)
            Call MEF90Element_Destroy(elemDamage,ierr)
         Case (MEF90DefMech_defectLawGradientDamage)
            If (Associated(MEF90DefMechCtx%temperature)) Then
               QuadratureOrder = 2 * max(elemDisplacementType%order - 1, elemDamageType%order) + 2 * elemDamageType%order
            Else
               QuadratureOrder = 2 * elemDisplacementType%order - 2 + 2 * elemDamageType%order
            End If
            Call MEF90Element_Create(mesh,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
            Call MEF90Element_Create(mesh,setIS,elemDamage,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
            !!! Contribution of the bilinear form
            Call MEF90GradDamageDispOperatorSet(residualSec,mesh,MEF90DefMechCtx%DMScal,xSec,damageSec,cellSetOptions%residualStiffness,setIS,matPropSet%HookesLaw,elemDisplacement,elemDisplacementType,elemDamage,elemDamageType,ierr)
            !!! temperature
            If (Associated(MEF90DefMechCtx%temperature)) Then
               Call MEF90GradDamageDispInelasticStrainRHSSetVertex(residualSec,mesh,MEF90DefMechCtx%DMScal,temperatureSec,matPropSet%HookesLaw*matPropSet%LinearThermalExpansion,damageSec,cellSetOptions%residualStiffness,setIS,elemDisplacement,elemDisplacementType,elemDamage,elemDamageType,ierr)
            End If
            !!! plastic Strain
            If (Associated(MEF90DefMechCtx%plasticStrain)) Then
               Call MEF90GradDamageDispInelasticStrainRHSSetCell(residualSec,MEF90DefMechCtx%DM,MEF90DefMechCtx%CellDMMatS,MEF90DefMechCtx%DMScal,plasticStrainSec,matPropSet%HookesLaw,damageSec,cellSetOptions%residualStiffness,setIS,elemDisplacement,elemDisplacementType,elemDamage,elemDamageType,ierr)
            End If
            !!! Force
            Call MEF90ElasticityForceRHSSetCell(residualSec,mesh,forceSec,setIS,elemDisplacement,elemDisplacementType,ierr);CHKERRQ(ierr)

            !!! pressure Force
            If (Associated(MEF90DefMechCtx%pressureForce) .AND. (elemDisplacementType%coDim > 0)) Then
               Call MEF90ElasticityPressureForceRHSSetCell(residualSec,mesh,pressureForceSec,setIS,elemDisplacement,elemDisplacementType,ierr);CHKERRQ(ierr)
            End If
            Call MEF90Element_Destroy(elemDisplacement,ierr)
            Call MEF90Element_Destroy(elemDamage,ierr)
         
         Case (MEF90DefMech_defectLawPlasticity)
            Print*,__FUNCT__,': Unimplemented Damage law',cellSetOptions%defectLaw
            STOP      
         End Select
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do   

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
               Call MEF90ISCreateCelltoVertex(MEF90DefMechCtx%DMVect,PETSC_COMM_WORLD,setIS,bcIS,ierr)
               Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
               Allocate(xPtr(nval),stat=ierr)
               Allocate(residualPtr(nval),stat=ierr)
               Allocate(boundaryDisplacementPtr(nval),stat=ierr)
               Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMVect,x,xPtr,bcIS,c,ierr)
               Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMVect,MEF90DefMechCtx%boundaryDisplacement,boundaryDisplacementPtr,bcIS,c,ierr)
               residualPtr = xPtr - boundaryDisplacementPtr
               Call MEF90VecSetValuesISdof(MEF90DefMechCtx%DMVect,residual,residualPtr,bcIS,c,INSERT_VALUES,ierr)
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
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
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
               Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMVect,x,xPtr,bcIS,c,ierr)
               Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMVect,MEF90DefMechCtx%boundaryDisplacement,boundaryDisplacementPtr,bcIS,c,ierr)
               residualPtr = xPtr - boundaryDisplacementPtr
               Call MEF90VecSetValuesISdof(MEF90DefMechCtx%DMVect,residual,residualPtr,bcIS,c,INSERT_VALUES,ierr)
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

      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
         Call VecScatterDestroy(ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)      
      End If
      
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVecScal,ierr);CHKERRQ(ierr)
      End If

      If (Associated(MEF90DefMechCtx%pressureForce)) Then
         Call SectionRealDestroy(pressureForceSec,ierr);CHKERRQ(ierr)
         Call VecScatterDestroy(ScatterSecToVecCellScal,ierr);CHKERRQ(ierr)      
      End If
      
      If (Associated(MEF90DefMechCtx%damage)) Then
         Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)
         Call VecScatterDestroy(ScatterSecToVecScal,ierr);CHKERRQ(ierr)
      End If
      
      Call SectionRealDestroy(forceSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(boundaryDisplacementSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(residualSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)

      Call VecScatterDestroy(ScatterSecToVecCell,ierr);CHKERRQ(ierr)      
      Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDisplacement

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDisplacement"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDisplacement:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechBilinearFormDisplacement(snesDispl,x,A,M,flg,MEF90DefMechCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesDispl
      Type(Vec),Intent(IN)                               :: x
      Type(Mat),Intent(INOUT)                            :: A,M
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr  
         
      Type(SectionReal)                                  :: damageSec
      Type(VecScatter)                                   :: ScatterSecToVecScal
      Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),pointer      :: cellSetOptions
      Type(MEF90DefMechVertexSetOptions_Type),pointer    :: vertexSetOptions
      Type(DM)                                           :: mesh
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemDamage
      Type(MEF90Element_Type)                            :: ElemDisplacementType,ElemDamageType
      PetscInt,Dimension(:),Pointer                      :: setIdx,bcIdx,Cone
      Type(IS)                                           :: bcIS
      PetscInt                                           :: cell,v,numBC,numDoF,numCell,c,dim
      
      Call MatZeroEntries(A,ierr);CHKERRQ(ierr)
      Call SNESGetDM(snesDispl,mesh,ierr);CHKERRQ(ierr)
      
      Call DMMeshGetSectionReal(MEF90DefMechCtx%DMScal,'default',damageSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMScal,damageSec,ScatterSecToVecScal,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(damageSec,ScatterSecToVecScal,SCATTER_REVERSE,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)


      Call DMMeshGetDimension(mesh,dim,ierr);CHKERRQ(ierr)
      Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         ElemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
         ElemDamageType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)
         !!! Call proper local assembly depending on the type of damage law
         Select Case (cellSetOptions%defectLaw)
         Case (MEF90DefMech_defectLawElasticity,MEF90DefMech_defectLawPlasticity)
            !!! Elements of codim > 0 have no contribution to the bilinear form, so we may as well skip them
            If (ElemDisplacementType%coDim == 0) Then
               QuadratureOrder = 2 * elemDisplacementType%order
               Call MEF90Element_Create(mesh,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
               Call MEF90ElasticityBilinearFormSet(A,mesh,setIS,matPropSet%HookesLaw,elemDisplacement,ElemDisplacementType,ierr);CHKERRQ(ierr)
               Call MEF90Element_Destroy(elemDisplacement,ierr)
            End If
         Case (MEF90DefMech_defectLawGradientDamage)
         !!! switch base don ATnum
            !!! Elements of codim > 0 have no contribution to the bilinear form, so we may as well skip them
            If (ElemDisplacementType%coDim == 0) Then
               QuadratureOrder = 2 * elemDisplacementType%order + 2 * ElemDamageType%order
               Call MEF90Element_Create(mesh,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
               Call MEF90Element_Create(mesh,setIS,elemDamage,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
               Call MEF90GradDamageDispBilinearFormSet(A,mesh,MEF90DefMechCtx%DMScal,setIS,matPropSet%HookesLaw,damageSec, &
                    cellSetOptions%residualStiffness,elemDisplacement,elemDisplacementType,elemDamage,elemDamageType,ierr)
               Call MEF90Element_Destroy(elemDisplacement,ierr)
               Call MEF90Element_Destroy(elemDamage,ierr)
            End If
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
         Call MEF90ISCreateCelltoVertex(mesh,PETSC_COMM_WORLD,setIS,bcIS,ierr)
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
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
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
      
      Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVecScal,ierr);CHKERRQ(ierr)
      
      flg = SAME_NONZERO_PATTERN
   End Subroutine MEF90DefMechBilinearFormDisplacement
   
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechWork"
!!!
!!!  
!!!  MEF90DefMechWork:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechWork(xVec,MEF90DefMechCtx,work,ierr)
      Type(Vec),Intent(IN)                               :: xVec
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,Dimension(:),Pointer                     :: work
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(SectionReal)                                  :: xSec
      Type(SectionReal)                                  :: forceSec,pressureForceSec
      Type(IS)                                           :: CellSetGlobalIS,setIS
      PetscInt,Dimension(:),Pointer                      :: setID
      PetscInt                                           :: set,QuadratureOrder,numCell
      Type(MEF90_MATPROP),pointer                        :: matpropSet
      Type(VecScatter)                                   :: ScatterSecToVec
      Type(VecScatter)                                   :: ScatterSecToVecCell,ScatterSecToVecCellScal
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90Element_Type)                            :: elemDisplacementType
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: globalOptions
      Type(MEF90DefMechCellSetOptions_Type),pointer      :: cellSetOptions
      PetscReal                                          :: myWork
      
      
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%globalOptionsBag,globalOptions,ierr);CHKERRQ(ierr)
      !!! Create dof-based sections
      Call DMMeshGetSectionReal(MEF90DefMechCtx%DMVect,'default',xSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMVect,xSec,ScatterSecToVec,ierr);CHKERRQ(ierr)

      !!! Create cell based sections, and allocate required pointers
      !!! I could check if any of these vectors are null 
      Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMVect,'default',ForceSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMVect,forceSec,ScatterSecToVecCell,ierr);CHKERRQ(ierr)
      Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMScal,'default',pressureForceSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMScal,pressureForceSec,ScatterSecToVecCellScal,ierr);CHKERRQ(ierr)
   
      !!! Scatter data from Vec to Sec, or initialize
      Call SectionRealToVec(xSec,ScatterSecToVec,SCATTER_REVERSE,xVec,ierr);CHKERRQ(ierr)         
      Call SectionRealToVec(forceSec,ScatterSecToVecCell,SCATTER_REVERSE,MEF90DefMechCtx%force,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(pressureForceSec,ScatterSecToVecCellScal,SCATTER_REVERSE,MEF90DefMechCtx%pressureForce,ierr);CHKERRQ(ierr)

      Work   = 0.0_Kr
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         mywork = 0.0_Kr
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)

         !!! Call proper local assembly depending on the type of damage law
         QuadratureOrder = elemDisplacementType%order
         Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)

         !!! Force work
         If (globalOptions%forceScaling /= MEF90Scaling_Null) Then
            Call MEF90ElasticityWorkSetCell(mywork,xSec,MEF90DefMechCtx%CellDMVect,ForceSec,setIS,elemDisplacement,elemDisplacementType,ierr)
         End If

         !!! pressure force work
         If ((globalOptions%pressureForceScaling /= MEF90Scaling_Null) .AND. (elemDisplacementType%coDim > 0)) Then
            Call MEF90ElasticityPressureWorkSetCell(mywork,xSec,MEF90DefMechCtx%CellDMVect,pressureForceSec,setIS,elemDisplacement,elemDisplacementType,ierr)
         End If
         Call MEF90Element_Destroy(elemDisplacement,ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)

         Call MPI_AllReduce(MPI_IN_PLACE,myWork,1,MPIU_SCALAR,MPI_SUM,MEF90DefMechCtx%MEF90Ctx%comm,ierr);CHKERRQ(ierr)
         work(set) = work(set) + mywork
      End Do ! set
      Call VecScatterDestroy(ScatterSecToVecCellScal,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVecCell,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(forceSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(pressureForceSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechWork   

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechElasticEnergy"
!!!
!!!  
!!!  MEF90DefMechElasticEnergy: 
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechElasticEnergy(x,MEF90DefMechCtx,energy,ierr)
      Type(Vec),Intent(IN)                               :: x
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,dimension(:),Pointer                     :: energy
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(SectionReal)                                  :: xSec
      Type(SectionReal)                                  :: damageSec,plasticStrainSec,temperatureSec
      Type(IS)                                           :: CellSetGlobalIS,setIS,setISdof,bcIS
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),pointer      :: cellSetOptions
      Type(VecScatter)                                   :: ScatterSecToVec,ScatterSecToVecScal,ScatterSecToVecCellMatS
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemScal
      Type(MEF90Element_Type)                            :: elemDisplacementType,elemScalType
      PetscReal                                          :: myenergy

      !!! Create dof-based sections
      Call DMMeshGetSectionReal(MEF90DefMechCtx%DMVect,'default',xSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMVect,xSec,ScatterSecToVec,ierr);CHKERRQ(ierr)

      !!! Scatter data from Vec to Sec, or initialize
      Call SectionRealToVec(xSec,ScatterSecToVec,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr) 

      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      
      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMMatS,'default',plasticStrainSec,ierr);CHKERRQ(ierr)
         Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMMatS,plasticStrainSec,ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticStrainSec,ScatterSecToVecCellMatS,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)          
      Else
         PlasticStrainSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call DMMeshGetSectionReal(MEF90DefMechCtx%DMScal,'default',temperatureSec,ierr);CHKERRQ(ierr)
         Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMScal,temperatureSec,ScatterSecToVecScal,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(temperatureSec,ScatterSecToVecScal,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr)          
      Else
         temperatureSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%damage)) Then
         If (Associated(MEF90DefMechCtx%temperature)) Then
            Call SectionRealDuplicate(temperatureSec,damageSec,ierr);CHKERRQ(ierr)
         Else
            Call DMMeshGetSectionReal(MEF90DefMechCtx%DMScal,'default',damageSec,ierr);CHKERRQ(ierr)
            Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMScal,damageSec,ScatterSecToVecScal,ierr);CHKERRQ(ierr)
         End If
         Call SectionRealToVec(damageSec,ScatterSecToVecScal,SCATTER_REVERSE,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)          
      Else
         damageSec%v = 0
      End If

      Do set = 1,size(setID)
         myenergy = 0.0_Kr
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DMVect,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
         elemScalType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)

         !!! Call proper local assembly depending on the type of damage law
         Select Case (cellSetOptions%defectLaw)
         Case (MEF90DefMech_defectLawElasticity,MEF90DefMech_defectLawPlasticity)
            If (elemDisplacementType%coDim == 0) Then
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  QuadratureOrder = 2 * max(elemDisplacementType%order - 1, elemScalType%order)
               Else
                  QuadratureOrder = 2 * elemDisplacementType%order - 2
               End If
               Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
               Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
               Call MEF90ElasticityEnergySet(myenergy,xSec,plasticStrainSec,temperatureSec,MEF90DefMechCtx%DMVect,MEF90DefMechCtx%DMScal,setIS, &
                                             matpropSet%HookesLaw,matpropSet%LinearThermalExpansion,elemDisplacement,elemDisplacementType,elemScal,elemScalType,ierr)
               Call MEF90Element_Destroy(elemDisplacement,ierr)
               Call MEF90Element_Destroy(elemScal,ierr)
            End If
         Case (MEF90DefMech_defectLawGradientDamage)
            If (elemDisplacementType%coDim == 0) Then
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  QuadratureOrder = 2 * max(elemDisplacementType%order - 1, elemScalType%order) + 2 * elemScalType%order
               Else
                  QuadratureOrder = 2 * elemDisplacementType%order - 2 + 2 * elemScalType%order
               End If
               Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
               Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
               Call MEF90GradDamageElasticEnergySet(myenergy,xSec,damageSec,plasticStrainSec,temperatureSec,MEF90DefMechCtx%DMVect,MEF90DefMechCtx%DMScal,setIS, &
                                             CellSetOptions%residualStiffness,matpropSet%HookesLaw,matpropSet%LinearThermalExpansion,elemDisplacement,elemDisplacementType,elemScal,elemScalType,ierr)
               Call MEF90Element_Destroy(elemDisplacement,ierr)
               Call MEF90Element_Destroy(elemScal,ierr)
            End If
         End Select
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         Call MPI_AllReduce(MPI_IN_PLACE,myenergy,1,MPIU_SCALAR,MPI_SUM,MEF90DefMechCtx%MEF90Ctx%comm,ierr);CHKERRQ(ierr)
         energy(set) = energy(set) + myenergy
      End Do ! set

      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
         Call VecScatterDestroy(ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)                   
      End If

      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
         Call VecScatterDestroy(ScatterSecToVecScal,ierr);CHKERRQ(ierr)
      End If
   End Subroutine MEF90DefMechElasticEnergy

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamage"
!!!
!!!  
!!!  MEF90DefMechOperatorDamage:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechOperatorDamage(snesDamage,alpha,residual,MEF90DefMechCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesDamage
      Type(Vec),Intent(IN)                               :: alpha
      Type(Vec),Intent(INOUT)                            :: residual
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(DM)                                           :: mesh
      Type(SectionReal)                                  :: alphaSec,residualSec
      Type(SectionReal)                                  :: plasticStrainSec,temperatureSec,displacementSec
      Type(SectionReal)                                  :: boundaryDamageSec
      PetscReal,Dimension(:),Pointer                     :: alphaPtr,residualPtr
      PetscReal,Dimension(:),Pointer                     :: boundaryDamagePtr
      Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof,bcIS
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),pointer      :: cellSetOptions
      Type(MEF90DefMechVertexSetOptions_Type),pointer    :: vertexSetOptions
      Type(VecScatter)                                   :: ScatterSecToVecScal,ScatterSecToVecVect
      Type(VecScatter)                                   :: ScatterSecToVecCell,ScatterSecToVecCellScal,ScatterSecToVecCellMatS
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemDamage
      Type(MEF90Element_Type)                            :: elemDisplacementType,elemDamageType
      PetscInt                                           :: p,dof,numDof
      PetscInt                                           :: nVal
      PetscReal                                          :: negOne = -1.0_Kr

      Call SNESGetDM(snesDamage,mesh,ierr);CHKERRQ(ierr)
      !!! Create dof-based sections and Scatter from the Vecs
      Call DMMeshGetSectionReal(MEF90DefMechCtx%DMScal,'default',alphaSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMScal,alphaSec,ScatterSecToVecScal,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(alphaSec,residualSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(alphaSec,boundaryDamageSec,ierr);CHKERRQ(ierr)
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDuplicate(alphaSec,temperatureSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(temperatureSec,ScatterSecToVecScal,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr)          
      Else
         temperatureSec%v = 0
      End If
      Call DMMeshGetSectionReal(MEF90DefMechCtx%DMVect,'default',displacementSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMVect,displacementSec,ScatterSecToVecVect,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(displacementSec,ScatterSecToVecVect,SCATTER_REVERSE,MEF90DefMechCtx%displacement,ierr);CHKERRQ(ierr)          
      
      !!! Create cell based sections, and allocate required pointers
      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMMatS,'default',plasticStrainSec,ierr);CHKERRQ(ierr)
         Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMMatS,plasticStrainSec,ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticStrainSec,ScatterSecToVecCellMatS,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)          
      Else
         PlasticStrainSec%v = 0
      End If

      !!! Scatter data from Vec to Sec, or initialize
      Call SectionRealToVec(alphaSec,ScatterSecToVecScal,SCATTER_REVERSE,alpha,ierr);CHKERRQ(ierr) 
      Call SectionRealToVec(boundaryDamageSec,ScatterSecToVecScal,SCATTER_REVERSE,MEF90DefMechCtx%boundaryDamage,ierr);CHKERRQ(ierr)

      Call SectionRealSet(residualSec,0.0_Kr,ierr);CHKERRQ(ierr)
      Call VecSet(residual,0.0_Kr,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
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
         elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
         elemDamageType       = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)

         !!! Call proper local assembly depending on the type of damage law
         Select Case (cellSetOptions%defectLaw)
         Case (MEF90DefMech_defectLawElasticity,MEF90DefMech_defectLawPlasticity)
            If (elemDisplacementType%coDim == 0) Then
               Select Case(cellSetOptions%gradientDamageLaw)
               Case(MEF90DefMech_defectLawGradientDamageAT1)
                  QuadratureOrder = max(elemDisplacementType%order, 2 * elemDamageType%order - 2)
                  Call MEF90Element_Create(mesh,setIS,elemDamage,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
                  Call MEF90GradDamageDamageOperatorSetAT1Elastic(residualSec,mesh,alphaSec,setIS,matPropSet%internalLength,matPropSet%FractureToughness,elemDamage,elemDamageType,ierr)
                  Call MEF90GradDamageDamageRHSSetAT1Elastic(residualSec,negone,mesh,setIS,matPropSet%internalLength,matPropSet%FractureToughness,elemDamage,elemDamageType,ierr)
                  Call MEF90Element_Destroy(elemDamage,ierr)
               Case(MEF90DefMech_defectLawGradientDamageAT2)
                  QuadratureOrder = 2 * elemDamageType%order
                  Call MEF90Element_Create(mesh,setIS,elemDamage,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
                  Call MEF90GradDamageDamageOperatorSetAT2Elastic(residualSec,mesh,alphaSec,setIS,matPropSet%internalLength,matPropSet%FractureToughness,elemDamage,elemDamageType,ierr)
                  Call MEF90Element_Destroy(elemDamage,ierr)
               Case default
                  Print*,__FUNCT__,': Unimplemented gradient damage law',cellSetOptions%gradientDamageLaw
                  STOP  
               End Select                
            End If
         Case (MEF90DefMech_defectLawGradientDamage)
            If (elemDisplacementType%coDim == 0) Then
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  QuadratureOrder = 2 * max(elemDisplacementType%order - 1, elemDamageType%order) + 2 * elemDamageType%order
               Else
                  QuadratureOrder = 2 * elemDisplacementType%order - 2 + 2 * elemDamageType%order
               End If
               Call MEF90Element_Create(mesh,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
               Call MEF90Element_Create(mesh,setIS,elemDamage,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
               Select Case(cellSetOptions%gradientDamageLaw)
               Case(MEF90DefMech_defectLawGradientDamageAT1)
                  Call MEF90GradDamageDamageOperatorSetAT1(residualSec,mesh,MEF90DefMechCtx%DMVect,alphaSec,setIS,displacementSec,temperatureSec,plasticStrainSec, & 
                                                           matPropSet%internalLength,matPropSet%HookesLaw,matPropSet%LinearThermalExpansion,matPropSet%FractureToughness, &
                                                           elemDamage,elemDamageType,elemDisplacement,elemDisplacementType,ierr)
                  Call MEF90GradDamageDamageRHSSetAT1(residualSec,negone,mesh,MEF90DefMechCtx%DMVect,setIS,displacementSec,temperatureSec,plasticStrainSec, &
                                                      matPropSet%internalLength,matPropSet%HookesLaw,matPropSet%LinearThermalExpansion,matPropSet%FractureToughness, &
                                                      elemDamage,elemDamageType,elemDisplacement,elemDisplacementType,ierr)
               Case(MEF90DefMech_defectLawGradientDamageAT2)
                  Call MEF90GradDamageDamageOperatorSetAT2(residualSec,mesh,MEF90DefMechCtx%DMVect,alphaSec,setIS,displacementSec,temperatureSec,plasticStrainSec, &
                                                           matPropSet%internalLength,matPropSet%HookesLaw,matPropSet%LinearThermalExpansion,matPropSet%FractureToughness, &
                                                           elemDamage,elemDamageType,elemDisplacement,elemDisplacementType,ierr)
                  Call MEF90GradDamageDamageRHSSetAT2(residualSec,negone,mesh,MEF90DefMechCtx%DMVect,setIS,displacementSec,temperatureSec,plasticStrainSec, &
                                                      matPropSet%internalLength,matPropSet%HookesLaw,matPropSet%LinearThermalExpansion,matPropSet%FractureToughness, &
                                                      elemDamage,elemDamageType,elemDisplacement,elemDisplacementType,ierr)
               Case default
                  Print*,__FUNCT__,': Unimplemented gradient damage law',cellSetOptions%gradientDamageLaw
                  STOP  
               End Select                
               Call MEF90Element_Destroy(elemDisplacement,ierr)
               Call MEF90Element_Destroy(elemDamage,ierr)
            End If
         End Select
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do   

      !!! "Ghost update" for the residual Section
      Call SectionRealComplete(residualSec,ierr);CHKERRQ(ierr)
      !!! Scatter back from SectionReal to Vec
      Call SectionRealToVec(residualSec,ScatterSecToVecScal,SCATTER_FORWARD,residual,ierr);CHKERRQ(ierr)

      !!!
      !!! Cell set BC
      !!!
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         If (cellSetOptions%Has_DamageBC) Then
            Call DMMeshGetStratumIS(MEF90DefMechCtx%DMScal,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
            Call MEF90ISCreateCelltoVertex(MEF90DefMechCtx%DMScal,PETSC_COMM_WORLD,setIS,bcIS,ierr)
            Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
            Allocate(alphaPtr(nval),stat=ierr)
            Allocate(residualPtr(nval),stat=ierr)
            Allocate(boundaryDamagePtr(nval),stat=ierr)
            Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMScal,alpha,alphaPtr,bcIS,1,ierr)
            Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMScal,MEF90DefMechCtx%boundaryDamage,boundaryDamagePtr,bcIS,1,ierr)
            residualPtr = alphaPtr - boundaryDamagePtr
            Call MEF90VecSetValuesISdof(MEF90DefMechCtx%DMScal,residual,residualPtr,bcIS,1,INSERT_VALUES,ierr)
            DeAllocate(boundaryDamagePtr)
            DeAllocate(residualPtr)
            DeAllocate(alphaPtr)
            Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If ! cellSetOptions%Has_BC
      End Do ! set
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)

      
      !!!
      !!! Vertex set BC
      !!!
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxVertexSetOptions(MEF90DefMechCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
         If (vertexSetOptions%Has_DamageBC) Then
            Call DMMeshGetStratumIS(MEF90DefMechCtx%DMScal,'Vertex Sets',setID(set),bcIS,ierr);CHKERRQ(iErr)
            Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
            Allocate(alphaPtr(nval),stat=ierr)
            Allocate(residualPtr(nval),stat=ierr)
            Allocate(boundaryDamagePtr(nval),stat=ierr)
            Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMScal,alpha,alphaPtr,bcIS,1,ierr)
            Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMScal,MEF90DefMechCtx%boundaryDamage,boundaryDamagePtr,bcIS,1,ierr)
            residualPtr = alphaPtr - boundaryDamagePtr
            Call MEF90VecSetValuesISdof(MEF90DefMechCtx%DMScal,residual,residualPtr,bcIS,1,INSERT_VALUES,ierr)
            DeAllocate(boundaryDamagePtr)
            DeAllocate(residualPtr)
            DeAllocate(alphaPtr)
            Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
         End If ! vertexSetOptions%Has_BC
      End Do ! set
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
   
      Call VecAssemblyBegin(residual,ierr)
      Call VecAssemblyEnd(residual,ierr)

      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
         Call VecScatterDestroy(ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)      
      End If
      
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVecScal,ierr);CHKERRQ(ierr)
      End If

      Call SectionRealDestroy(displacementSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(boundaryDamageSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(residualSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(alphaSec,ierr);CHKERRQ(ierr)

      Call VecScatterDestroy(ScatterSecToVecScal,ierr);CHKERRQ(ierr)      
      Call VecScatterDestroy(ScatterSecToVecVect,ierr);CHKERRQ(ierr)   
      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call VecScatterDestroy(ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)
      End If
   End Subroutine MEF90DefMechOperatorDamage


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamage"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamage:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechBilinearFormDamage(snesDamage,alpha,A,M,flg,MEF90DefMechCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesDamage
      Type(Vec),Intent(IN)                               :: alpha
      Type(Mat),Intent(INOUT)                            :: A,M
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr  
         
      Type(SectionReal)                                  :: displacementSec,temperatureSec,plasticStrainSec   
      Type(VecScatter)                                   :: ScatterSecToVecVect,ScatterSecToVecScal,ScatterSecToVecCellMatS
      Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),pointer      :: cellSetOptions
      Type(MEF90DefMechVertexSetOptions_Type),pointer    :: vertexSetOptions
      Type(DM)                                           :: mesh
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemDamage
      Type(MEF90Element_Type)                            :: elemDisplacementType,elemDamageType
      PetscInt,Dimension(:),Pointer                      :: setIdx,bcIdx,Cone
      Type(IS)                                           :: bcIS
      PetscInt                                           :: cell,v,numBC,numDoF
      
      Call MatZeroEntries(A,ierr);CHKERRQ(ierr)
      Call SNESGetDM(snesDamage,mesh,ierr);CHKERRQ(ierr)

      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call DMMeshGetSectionReal(MEF90DefMechCtx%DMScal,'default',temperatureSec,ierr);CHKERRQ(ierr)
         Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMScal,temperatureSec,ScatterSecToVecScal,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(temperatureSec,ScatterSecToVecScal,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr)          
      Else
         temperatureSec%v = 0
      End If
      Call DMMeshGetSectionReal(MEF90DefMechCtx%DMVect,'default',displacementSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMVect,displacementSec,ScatterSecToVecVect,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(displacementSec,ScatterSecToVecVect,SCATTER_REVERSE,MEF90DefMechCtx%displacement,ierr);CHKERRQ(ierr)          
      
      !!! Create cell based sections, and allocate required pointers
      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMMatS,'default',plasticStrainSec,ierr);CHKERRQ(ierr)
         Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMMatS,plasticStrainSec,ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticStrainSec,ScatterSecToVecCellMatS,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)          
      Else
         PlasticStrainSec%v = 0
      End If

      Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
         elemDamageType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)
         !!! Call proper local assembly depending on the type of damage law
         Select Case (cellSetOptions%defectLaw)
         Case (MEF90DefMech_defectLawElasticity,MEF90DefMech_defectLawPlasticity)
            If (elemDamageType%coDim == 0) Then
               Select Case(cellSetOptions%gradientDamageLaw)
               Case(MEF90DefMech_defectLawGradientDamageAT1)
                  QuadratureOrder = max(elemDisplacementType%order, 2 * elemDamageType%order - 2)
                  Call MEF90Element_Create(mesh,setIS,elemDamage,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
                  Call MEF90GradDamageDamageBilinearFormSetAT1Elastic(A,mesh,setIS,matPropSet%internalLength,matPropSet%FractureToughness,elemDamage,elemDamageType,ierr)
                  Call MEF90Element_Destroy(elemDamage,ierr)
               Case(MEF90DefMech_defectLawGradientDamageAT2)
                  QuadratureOrder = 2 * elemDamageType%order
                  Call MEF90Element_Create(mesh,setIS,elemDamage,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
                  Call MEF90GradDamageDamageBilinearFormSetAT2Elastic(A,mesh,setIS,matPropSet%internalLength,matPropSet%FractureToughness,elemDamage,elemDamageType,ierr)
                  Call MEF90Element_Destroy(elemDamage,ierr)
               End Select
            End If
         Case (MEF90DefMech_defectLawGradientDamage)
            !!! Elements of codim > 0 have no contribution to the bilinear form, so we need to skip them
            If (elemDamageType%coDim == 0) Then
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  QuadratureOrder = 2 * max(elemDisplacementType%order - 1, elemDamageType%order) + 2 * elemDamageType%order
               Else
                  QuadratureOrder = 2 * elemDisplacementType%order - 2 + 2 * elemDamageType%order
               End If
               Call MEF90Element_Create(mesh,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
               Call MEF90Element_Create(mesh,setIS,elemDamage,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
               Select Case(cellSetOptions%gradientDamageLaw)
               Case(MEF90DefMech_defectLawGradientDamageAT1)
                  Call MEF90GradDamageDamageBilinearFormSetAT1(A,mesh,MEF90DefMechCtx%DMVect,setIS,displacementSec,temperatureSec,plasticStrainSec, & 
                                                           matPropSet%internalLength,matPropSet%HookesLaw,matPropSet%LinearThermalExpansion,matPropSet%FractureToughness, &
                                                           elemDamage,elemDamageType,elemDisplacement,elemDisplacementType,ierr)
               Case(MEF90DefMech_defectLawGradientDamageAT2)
                  Call MEF90GradDamageDamageBilinearFormSetAT2(A,mesh,MEF90DefMechCtx%DMVect,setIS,displacementSec,temperatureSec,plasticStrainSec, & 
                                                           matPropSet%internalLength,matPropSet%HookesLaw,matPropSet%LinearThermalExpansion,matPropSet%FractureToughness, &
                                                           elemDamage,elemDamageType,elemDisplacement,elemDisplacementType,ierr)
               Case default
                  Print*,__FUNCT__,': Unimplemented gradient damage law',cellSetOptions%gradientDamageLaw
                  STOP  
               End Select    
               Call MEF90Element_Destroy(elemDisplacement,ierr)
               Call MEF90Element_Destroy(elemDamage,ierr)
            End If
         End Select
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do     
      Call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
      Call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
      !!!
      !!! Boundary conditions at cell sets
      !!!
      !Call DMmeshGetStratumSize(mesh,"height",0,numCell,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call MEF90ISCreateCelltoVertex(mesh,PETSC_COMM_WORLD,setIS,bcIS,ierr)
         If (cellSetOptions%Has_damageBC) Then
            Call DMMeshISCreateISglobaldof(mesh,bcIS,0,setISdof,ierr);CHKERRQ(ierr)
            Call MatZeroRowsColumnsIS(A,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
            Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
         End If ! cellSetOptions%Has_damageBC
         Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do     
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
   
      !!!
      !!! Boundary conditions at vertex sets
      !!!
      Call DMmeshGetLabelIdIS(mesh,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxVertexSetOptions(MEF90DefMechCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
         If (vertexSetOptions%Has_damageBC) Then
            Call DMMeshGetStratumIS(mesh,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
            Call DMMeshISCreateISglobaldof(mesh,setIS,0,setISdof,ierr);CHKERRQ(ierr)
            Call MatZeroRowsColumnsIS(A,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
         End If
      End Do
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      
      Call SectionRealDestroy(displacementSec,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVecVect,ierr);CHKERRQ(ierr)
      If (temperatureSec%v /= 0) Then
         Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
         Call VecScatterDestroy(ScatterSecToVecScal,ierr);CHKERRQ(ierr)
      End If
      If (plasticStrainSec%v /= 0) Then
         Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
         Call VecScatterDestroy(ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)
      End If
      
      flg = SAME_NONZERO_PATTERN
   End Subroutine MEF90DefMechBilinearFormDamage
   

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSurfaceEnergy"
!!!
!!!  
!!!  MEF90DefMechSurfaceEnergy:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechSurfaceEnergy(alpha,MEF90DefMechCtx,energy,ierr)
      Type(Vec),Intent(IN)                               :: alpha
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,dimension(:),Pointer                     :: energy
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(SectionReal)                                  :: alphaSec
      Type(IS)                                           :: CellSetGlobalIS,setIS,setISdof,bcIS
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),pointer      :: cellSetOptions
      Type(VecScatter)                                   :: ScatterSecToVecScal
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemScal
      Type(MEF90Element_Type)                            :: elemScalType
      PetscReal                                          :: myenergy

      !!! Create dof-based sections
      Call DMMeshGetSectionReal(MEF90DefMechCtx%DMScal,'default',alphaSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMScal,alphaSec,ScatterSecToVecScal,ierr);CHKERRQ(ierr)

      !!! Scatter data from Vec to Sec, or initialize
      Call SectionRealToVec(alphaSec,ScatterSecToVecScal,SCATTER_REVERSE,alpha,ierr);CHKERRQ(ierr) 

      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      
      Do set = 1,size(setID)
         myenergy = 0.0_Kr
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         elemScalType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)

         !!! Call proper local assembly depending on the type of damage law
         Select Case (cellSetOptions%defectLaw)
         Case (MEF90DefMech_defectLawGradientDamage)
            If (elemScalType%coDim == 0) Then
               Select Case(cellSetOptions%gradientDamageLaw)
               Case(MEF90DefMech_defectLawGradientDamageAT1)
                  QuadratureOrder = max(2 * elemScalType%order - 2, elemScalType%order)
                  Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
                  Call MEF90GradDamageSurfaceEnergySetAT1(myenergy,alphaSec,MEF90DefMechCtx%DMScal,setIS,matpropSet%internalLength,matpropSet%fractureToughness,elemScal,elemScalType,ierr)
                  Call MEF90Element_Destroy(elemScal,ierr)
               Case(MEF90DefMech_defectLawGradientDamageAT2)
                  QuadratureOrder = 2*elemScalType%order
                  Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
                  Call MEF90GradDamageSurfaceEnergySetAT1(myenergy,alphaSec,MEF90DefMechCtx%DMScal,setIS,matpropSet%internalLength,matpropSet%fractureToughness,elemScal,elemScalType,ierr)
                  Call MEF90Element_Destroy(elemScal,ierr)
               End Select
            End If
         End Select
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         Call MPI_AllReduce(MPI_IN_PLACE,myenergy,1,MPIU_SCALAR,MPI_SUM,MEF90DefMechCtx%MEF90Ctx%comm,ierr);CHKERRQ(ierr)
         energy(set) = energy(set) + myenergy
      End Do ! set

      Call SectionRealDestroy(alphaSec,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVecScal,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechSurfaceEnergy
End Module MEF90_APPEND(m_MEF90_DefMechAssembly,MEF90_DIM)D
