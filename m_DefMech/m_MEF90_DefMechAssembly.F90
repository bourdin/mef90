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
Subroutine MEF90DefMechOperator(elem)
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elem
   Continue
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
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elem
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
            Call MEF90Element_Create(mesh,setIS,elem,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
            Call MEF90ElasticityBilinearFormSet(A,mesh,setIS,matPropSet%HookesLaw,elem,elemTypeDisplacement,ierr);CHKERRQ(ierr)
            Call MEF90Element_Destroy(elem,ierr)
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
         Call PetscBagGetDataMEF90DefMechCtxVertexSetOptions(MEF90DefMechCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
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
   End Subroutine MEF90DefMechBilinearForm
End Module MEF90_APPEND(m_MEF90_DefMechAssembly,MEF90_DIM)D
