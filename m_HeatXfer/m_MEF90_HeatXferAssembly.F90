#include "../MEF90/mef90.inc"
Module MEF90_APPEND(m_MEF90_HeatXferAssembly,MEF90_DIM)D
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_HeatXferCtx

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferOperator"
!!!
!!!  
!!!  MEF90HeatXferOperator:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine MEF90HeatXferOperator(snesTemp,x,residual,MEF90HeatXferCtx,ierr)
   Type(SNES),Intent(IN)                              :: snesTemp
   Type(Vec),Intent(IN)                               :: x
   Type(Vec),Intent(INOUT)                            :: residual
   Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
   PetscErrorCode,Intent(OUT)                         :: ierr

   Type(SectionReal)                                  :: xSec,residualSec
   Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof
   PetscInt,dimension(:),Pointer                      :: setID
   PetscInt,Dimension(:),Pointer                      :: setIdx
   PetscInt                                           :: set,QuadratureOrder
   Type(MEF90_MATPROP),pointer                        :: matpropSet
   Type(MEF90HeatXferCellSetOptions_Type),pointer     :: cellSetOptions
   Type(MEF90HeatXferVertexSetOptions_Type),pointer   :: vertexSetOptions
   PetscReal,Dimension(:),Pointer                     :: BC
   Type(DM)                                           :: mesh
   Type(VecScatter)                                   :: ScatterSecToVec
   Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elem
   Type(MEF90Element_Type)                            :: elemType
   
   Call SNESGetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)

   Call DMMeshGetSectionReal(mesh,'default',xSec,ierr);CHKERRQ(ierr)
   Call DMMeshCreateGlobalScatter(mesh,xSec,ScatterSecToVec,ierr);CHKERRQ(ierr)

   Call SectionRealDuplicate(xSec,residualSec,ierr);CHKERRQ(ierr)
   Call SectionRealSet(residualSec,0.0_Kr,ierr);CHKERRQ(ierr)
   Call VecSet(residual,0.0_Kr,ierr);CHKERRQ(ierr)

   Call SectionRealToVec(xSec,ScatterSecToVec,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr)

   Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
   Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)

      Call PetscBagGetDataMEF90MatProp(MEF90HeatXferCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)

      elemType = MEF90_knownElements(cellSetOptions%ElemTypeShortID)
      QuadratureOrder = elemType%order * 2
      Call MEF90Element_Create(mesh,setIS,elem,QuadratureOrder,CellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)

      Call MEF90Diffusion_OperatorSet(residualSec,mesh,xSec,setIS,matpropSet%ThermalConductivity,cellSetOptions%SurfaceThermalConductivity,elem,elemType,ierr);CHKERRQ(ierr)

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
   !!! Zero-out BC entries in the residual
   !!!
   Call DMmeshGetLabelIdIS(mesh,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
   Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      Call PetscBagGetDataMEF90HeatXferCtxVertexSetOptions(MEF90HeatXferCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
      If (vertexSetOptions%Has_BC) Then
         Call DMMeshGetStratumIS(mesh,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call DMMeshISCreateISglobaldof(mesh,xSec,setIS,0,setISdof,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
         Allocate(BC(size(setIdx)))
         BC = 0.0_Kr
         Call VecSetValues(residual,size(setIdx),setIdx,BC,INSERT_VALUES,ierr);CHKERRQ(ierr)
         DeAllocate(BC)
         Call ISRestoreIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
      End If
   End Do
   Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)

   Call VecAssemblyBegin(residual,ierr)
   Call VecAssemblyEnd(residual,ierr)
   Call SectionRealDestroy(residualSec,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
   Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
End Subroutine MEF90HeatXferOperator

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferBilinearForm"
!!!
!!!  
!!!  MEF90HeatXferBilinearForm:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine MEF90HeatXferBilinearForm(snesTemp,x,A,M,flg,MEF90HeatXferCtx,ierr)
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
   Type(SectionReal)                                  :: xSec
   Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elem
   Type(MEF90Element_Type)                            :: elemType
   
   Call MatZeroEntries(A,ierr);CHKERRQ(ierr)
   Call SNESGetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)
   Call DMMeshGetSectionReal(mesh,'default',xSec,ierr);CHKERRQ(ierr)
   !!! 
   !!! No need to zero out Sec because it is only used as a PetscSection when assembling Mat
   !!!
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
      Call MEF90Diffusion_BilinearFormSet(A,mesh,xSec,setIS,matpropSet%ThermalConductivity,cellSetOptions%SurfaceThermalConductivity,elem,elemType,ierr);CHKERRQ(ierr)
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
         Call DMMeshISCreateISglobaldof(mesh,xSec,setIS,0,setISdof,ierr);CHKERRQ(ierr)
         Call MatZeroRowsColumnsIS(A,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
      End If
   End Do
   Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
   
   flg = SAME_NONZERO_PATTERN
End Subroutine MEF90HeatXferBilinearForm

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferRHS_Cst"
!!!
!!!  
!!!  MEF90HeatXferRHS_Cst: Build the time dependent RHS for the SNES with CST fluxes scaled by time
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine MEF90HeatXferRHS_Cst(snesTemp,rhs,t,MEF90HeatXferCtx,ierr)
   Type(SNES),Intent(IN)                              :: snesTemp
   Type(Vec),Intent(IN)                               :: rhs
   PetscReal,Intent(IN)                               :: t
   Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
   PetscErrorCode,Intent(OUT)                         :: ierr
   
   Type(SectionReal)                                  :: rhsSec
   Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof
   PetscInt,dimension(:),Pointer                      :: setID,setIdx
   PetscInt                                           :: set,QuadratureOrder
   Type(MEF90HeatXferCellSetOptions_Type),pointer     :: cellSetOptions
   Type(MEF90HeatXferVertexSetOptions_Type),pointer   :: vertexSetOptions
   PetscReal,Dimension(:),Pointer                     :: BC
   PetscReal                                          :: F
   Type(DM)                                           :: mesh
   Type(VecScatter)                                   :: ScatterSecToVec
   Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elem
   Type(MEF90Element_Type)                            :: elemType

   Call SNESGetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)

   Call DMMeshGetSectionReal(mesh,'default',rhsSec,ierr);CHKERRQ(ierr)
   Call DMMeshCreateGlobalScatter(mesh,rhsSec,ScatterSecToVec,ierr);CHKERRQ(ierr)

   Call SectionRealSet(rhsSec,0.0_Kr,ierr);CHKERRQ(ierr)
   Call VecSet(rhs,0.0_Kr,ierr);CHKERRQ(ierr)

   Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
   Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
      Call PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)

      elemType = MEF90_knownElements(cellSetOptions%ElemTypeShortID)
      QuadratureOrder = elemType%order * 2
      Call MEF90Element_Create(mesh,setIS,elem,QuadratureOrder,CellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)

      F = t * (cellSetOptions%Flux + cellSetOptions%SurfaceThermalConductivity * cellSetOptions%externalTemp) 
      Call MEF90Diffusion_RHSSet(rhsSec,mesh,F,setIS,elem,elemType,ierr);CHKERRQ(ierr)

      Call MEF90Element_Destroy(elem,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   End Do     
   Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
   !!! "Ghost update" for the residual Section
   Call SectionRealComplete(rhsSec,ierr);CHKERRQ(ierr)
   !!! Scatter back from SectionReal to Vec
   Call SectionRealToVec(rhsSec,ScatterSecToVec,SCATTER_FORWARD,rhs,ierr);CHKERRQ(ierr)
   
   
   !!!
   !!! Zero-out BC entries in the residual
   !!!
   Call DMmeshGetLabelIdIS(mesh,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
   Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      Call PetscBagGetDataMEF90HeatXferCtxVertexSetOptions(MEF90HeatXferCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
      If (vertexSetOptions%Has_BC) Then
         Call DMMeshGetStratumIS(mesh,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call DMMeshISCreateISglobaldof(mesh,rhsSec,setIS,0,setISdof,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
         Allocate(BC(size(setIdx)))
         BC = 0.0_Kr
         Call VecSetValues(rhs,size(setIdx),setIdx,BC,INSERT_VALUES,ierr);CHKERRQ(ierr)
         DeAllocate(BC)
         Call ISRestoreIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
      End If
   End Do
   Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
   
   !!!
   !!! Cleanup
   !!!
   Call VecAssemblyBegin(rhs,ierr)
   Call VecAssemblyEnd(rhs,ierr)
   Call SectionRealDestroy(rhsSec,ierr);CHKERRQ(ierr)
   Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
End Subroutine MEF90HeatXferRHS_Cst

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferRHS"
!!!
!!!  
!!!  MEF90HeatXferRHS_Cst: Build the time dependent RHS for the SNES with CST fluxes scaled by time
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine MEF90HeatXferRHS(snesTemp,rhs,flux,reftemp,MEF90HeatXferCtx,ierr)
   Type(SNES),Intent(IN)                              :: snesTemp
   Type(Vec),Intent(IN)                               :: rhs
   Type(Vec),Intent(IN)                               :: flux
   Type(Vec),Intent(IN)                               :: refTemp
   Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
   PetscErrorCode,Intent(OUT)                         :: ierr
   
   Type(SectionReal)                                  :: rhsSec,fSec
   Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof
   PetscInt,dimension(:),Pointer                      :: setID,SetIdx
   PetscInt                                           :: set,QuadratureOrder,numCell
   Type(MEF90HeatXferCellSetOptions_Type),pointer     :: cellSetOptions
   Type(MEF90HeatXferVertexSetOptions_Type),pointer   :: vertexSetOptions
   Type(MEF90_MATPROP),pointer                        :: matpropSet
   PetscReal,Dimension(:),Pointer                     :: BC
   Type(DM)                                           :: mesh
   Type(VecScatter)                                   :: ScatterSecToVec
   Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elem
   Type(MEF90Element_Type)                            :: elemType
   Type(Vec)                                          :: modifiedFlux

   Call SNESGetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)

   Call DMMeshGetSectionReal(mesh,'default',rhsSec,ierr);CHKERRQ(ierr)
   Call DMMeshGetSectionReal(mesh,'default',fSec,ierr);CHKERRQ(ierr)
   Call DMMeshCreateGlobalScatter(mesh,rhsSec,ScatterSecToVec,ierr);CHKERRQ(ierr)

   Call SectionRealSet(rhsSec,0.0_Kr,ierr);CHKERRQ(ierr)
   Call VecSet(rhs,0.0_Kr,ierr);CHKERRQ(ierr)
   Call VecDuplicate(flux,modifiedFLux,ierr);CHKERRQ(ierr)

   Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
   Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
      Call PetscBagGetDataMEF90MatProp(MEF90HeatXferCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)

      !!! Modified flux is flux + surfaceThermalConductivity * refTemp
      Call VecCopy(flux,modifiedFlux,ierr);CHKERRQ(ierr)
      Call VecAXPY(modifiedFlux,cellSetOptions%SurfaceThermalConductivity,refTemp,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(fSec,ScatterSecToVec,SCATTER_REVERSE,modifiedFlux,ierr);CHKERRQ(ierr)
      
      elemType = MEF90_knownElements(cellSetOptions%ElemTypeShortID)
      QuadratureOrder = elemType%order * 2
      Call MEF90Element_Create(mesh,setIS,elem,QuadratureOrder,CellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)

      Call MEF90Diffusion_RHSSet(rhssec,mesh,fsec,setIS,elem,elemType,ierr);CHKERRQ(ierr)

      Call MEF90Element_Destroy(elem,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   End Do     
   Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
   !!! "Ghost update" for the residual Section
   Call SectionRealComplete(rhsSec,ierr);CHKERRQ(ierr)
   !!! Scatter back from SectionReal to Vec
   Call SectionRealToVec(rhsSec,ScatterSecToVec,SCATTER_FORWARD,rhs,ierr);CHKERRQ(ierr)
   Call VecDestroy(modifiedFlux,ierr);CHKERRQ(ierr)
   
   !!!
   !!! Zero-out BC entries in the residual
   !!!
   Call DMmeshGetLabelIdIS(mesh,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
   Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      Call PetscBagGetDataMEF90HeatXferCtxVertexSetOptions(MEF90HeatXferCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
      If (vertexSetOptions%Has_BC) Then
         Call DMMeshGetStratumIS(mesh,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call DMMeshISCreateISglobaldof(mesh,rhsSec,setIS,0,setISdof,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
         Allocate(BC(size(setIdx)))
         BC = 0.0_Kr
         Call VecSetValues(rhs,size(setIdx),setIdx,BC,INSERT_VALUES,ierr);CHKERRQ(ierr)
         DeAllocate(BC)
         Call ISRestoreIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
      End If
   End Do
   Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
   
   !!!
   !!! Cleanup
   !!!
   Call VecAssemblyBegin(rhs,ierr)
   Call VecAssemblyEnd(rhs,ierr)
   Call SectionRealDestroy(rhsSec,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(fSec,ierr);CHKERRQ(ierr)
   Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
End Subroutine MEF90HeatXferRHS
End Module MEF90_APPEND(m_MEF90_HeatXferAssembly,MEF90_DIM)D
