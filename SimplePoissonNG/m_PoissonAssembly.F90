#include "SimplePoisson.inc"
Module M_POISSONASSEMBLY
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_PoissonGlobalProperties
   Use M_POISSONCELLSETPROPERTIES
   Use m_PoissonVertexSetProperties
   Implicit NONE
Contains

#undef __FUNCT__
#define __FUNCT__ "SimplePoissonBilinearForm"
!!!
!!!  
!!!  SimplePoissonBilinearForm:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine SimplePoissonBilinearForm(snesTemp,x,A,M,flg,PoissonCtx,ierr)
   Type(SNES),Intent(IN)                           :: snesTemp
   Type(Vec),Intent(IN)                            :: x
   Type(Mat),Intent(INOUT)                         :: A,M
   MatStructure,Intent(INOUT)                      :: flg
   Type(MEF90Ctx_Type),intent(IN)                  :: PoissonCtx
   PetscErrorCode,Intent(OUT)                      :: ierr  
   
   Type(IS)                                        :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof
   PetscInt,dimension(:),Pointer                   :: setID
   PetscInt                                        :: set,QuadratureOrder
   Type(MEF90_MATPROP),pointer                     :: matpropSet
   Type(PoissonCellSetProperties_Type),pointer     :: cellSetProperties
   Type(PoissonVertexSetProperties_Type),pointer   :: vertexSetProperties
   Type(DM)                                        :: mesh
   Type(SectionReal)                               :: xSec
   Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer   :: elem
   Type(Element_Type)                              :: elemType
   
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

      Call PetscBagGetDataMEF90_MatProp(PoissonCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataPoissonCellSetProperties(PoissonCtx%CellSetPropertiesBag(set),cellSetProperties,ierr);CHKERRQ(ierr)

      elemType = MEF90_knownElements(cellSetProperties%ElemTypeShortID)

      QuadratureOrder = 2 * elemType%order
      Call MEF90_ElementCreate(mesh,setIS,elem,QuadratureOrder,CellSetProperties%ElemTypeShortID,ierr);CHKERRQ(ierr)
      Call MEF90_DiffusionBilinearFormSet(A,mesh,xSec,setIS,matpropSet%ThermalConductivity,cellSetProperties%SurfaceThermalConductivity,elem,elemType,ierr);CHKERRQ(ierr)
      Call MEF90_ElementDestroy(elem,ierr)

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
      Call PetscBagGetDataPoissonVertexSetProperties(PoissonCtx%VertexSetPropertiesBag(set),vertexSetProperties,ierr);CHKERRQ(ierr)
      If (vertexSetProperties%Has_BC) Then
         Call DMMeshGetStratumIS(mesh,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call DMMeshISCreateISglobaldof(mesh,xSec,setIS,0,setISdof,ierr);CHKERRQ(ierr)
         Call MatZeroRowsColumnsIS(A,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
      End If
   End Do
   Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
   
   flg = SAME_NONZERO_PATTERN
End Subroutine SimplePoissonBilinearForm


#undef __FUNCT__
#define __FUNCT__ "SimplePoissonOperator"
!!!
!!!  
!!!  SimplePoissonOperator:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine SimplePoissonOperator(snesTemp,x,residual,PoissonCtx,ierr)
   Type(SNES),Intent(IN)                           :: snesTemp
   Type(Vec),Intent(IN)                            :: x
   Type(Vec),Intent(INOUT)                         :: residual
   Type(MEF90Ctx_Type),intent(IN)                  :: PoissonCtx
   PetscErrorCode,Intent(OUT)                      :: ierr

   Type(SectionReal)                               :: xSec,residualSec
   Type(IS)                                        :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof
   PetscInt,dimension(:),Pointer                   :: setID
   PetscInt,Dimension(:),Pointer                   :: setIdx
   PetscInt                                        :: set,QuadratureOrder
   Type(MEF90_MATPROP),pointer                     :: matpropSet
   Type(PoissonCellSetProperties_Type),pointer     :: cellSetProperties
   Type(PoissonVertexSetProperties_Type),pointer   :: vertexSetProperties
   PetscReal,Dimension(:),Pointer                  :: BC
   PetscReal                                       :: F
   Type(DM)                                        :: mesh
   Type(VecScatter)                                :: ScatterSecToVec
   Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer   :: elem
   Type(Element_Type)                              :: elemType
   
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

      Call PetscBagGetDataMEF90_MatProp(PoissonCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataPoissonCellSetProperties(PoissonCtx%CellSetPropertiesBag(set),cellSetProperties,ierr);CHKERRQ(ierr)

      elemType = MEF90_knownElements(cellSetProperties%ElemTypeShortID)
      QuadratureOrder = elemType%order * 2
      Call MEF90_ElementCreate(mesh,setIS,elem,QuadratureOrder,CellSetProperties%ElemTypeShortID,ierr);CHKERRQ(ierr)

      Call MEF90_DiffusionOperatorSet(residualSec,mesh,xSec,setIS,matpropSet%ThermalConductivity,cellSetProperties%SurfaceThermalConductivity,elem,elemType,ierr);CHKERRQ(ierr)

      Call MEF90_ElementDestroy(elem,ierr);CHKERRQ(ierr)
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
      Call PetscBagGetDataPoissonVertexSetProperties(PoissonCtx%VertexSetPropertiesBag(set),vertexSetProperties,ierr);CHKERRQ(ierr)
      If (vertexSetProperties%Has_BC) Then
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
End Subroutine SimplePoissonOperator

#undef __FUNCT__
#define __FUNCT__ "SimplePoissonRHS_Cst"
!!!
!!!  
!!!  SimplePoissonRHS_Cst: Build the time dependent RHS for the SNES with CST fluxes scaled by time
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine SimplePoissonRHS_Cst(snesTemp,rhs,t,PoissonCtx,ierr)
   Type(SNES),Intent(IN)                           :: snesTemp
   Type(Vec),Intent(IN)                            :: rhs
   PetscReal,Intent(IN)                            :: t
   Type(MEF90Ctx_Type),intent(IN)                  :: PoissonCtx
   PetscErrorCode,Intent(OUT)                      :: ierr

   Type(SectionReal)                               :: rhsSec
   Type(IS)                                        :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof
   PetscInt,dimension(:),Pointer                   :: setID,setIdx
   PetscInt                                        :: set,QuadratureOrder
   Type(PoissonCellSetProperties_Type),pointer     :: cellSetProperties
   Type(PoissonVertexSetProperties_Type),pointer   :: vertexSetProperties
   PetscReal,Dimension(:),Pointer                  :: BC
   PetscReal                                       :: F
   Type(DM)                                        :: mesh
   Type(VecScatter)                                :: ScatterSecToVec
   Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer   :: elem
   Type(Element_Type)                              :: elemType

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
      Call PetscBagGetDataPoissonCellSetProperties(PoissonCtx%CellSetPropertiesBag(set),cellSetProperties,ierr);CHKERRQ(ierr)

      elemType = MEF90_knownElements(cellSetProperties%ElemTypeShortID)
      QuadratureOrder = elemType%order * 2
      Call MEF90_ElementCreate(mesh,setIS,elem,QuadratureOrder,CellSetProperties%ElemTypeShortID,ierr);CHKERRQ(ierr)

      F = t * (cellSetProperties%Flux + cellSetProperties%SurfaceThermalConductivity * cellSetProperties%referenceTemp) 
      Call MEF90_DiffusionRHSSet(rhsSec,mesh,F,setIS,elem,elemType,ierr);CHKERRQ(ierr)

      Call MEF90_ElementDestroy(elem,ierr);CHKERRQ(ierr)
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
      Call PetscBagGetDataPoissonVertexSetProperties(PoissonCtx%VertexSetPropertiesBag(set),vertexSetProperties,ierr);CHKERRQ(ierr)
      If (vertexSetProperties%Has_BC) Then
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
End Subroutine SimplePoissonRHS_Cst

#undef __FUNCT__
#define __FUNCT__ "SimplePoissonRHS"
!!!
!!!  
!!!  SimplePoissonRHS_Cst: Build the time dependent RHS for the SNES with CST fluxes scaled by time
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine SimplePoissonRHS(snesTemp,rhs,flux,reftemp,PoissonCtx,ierr)
   Type(SNES),Intent(IN)                           :: snesTemp
   Type(Vec),Intent(IN)                            :: rhs
   Type(Vec),Intent(IN)                            :: flux
   Type(Vec),Intent(IN)                            :: refTemp
   Type(MEF90Ctx_Type),intent(IN)                  :: PoissonCtx
   PetscErrorCode,Intent(OUT)                      :: ierr

   Type(SectionReal)                               :: rhsSec,fSec
   Type(IS)                                        :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof
   PetscInt,dimension(:),Pointer                   :: setID,SetIdx
   PetscInt                                        :: set,QuadratureOrder,numCell
   Type(PoissonCellSetProperties_Type),pointer     :: cellSetProperties
   Type(PoissonVertexSetProperties_Type),pointer   :: vertexSetProperties
   PetscReal,Dimension(:),Pointer                  :: BC
   Type(DM)                                        :: mesh
   Type(VecScatter)                                :: ScatterSecToVec
   Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer   :: elem
   Type(Element_Type)                              :: elemType
   Type(Vec)                                       :: modifiedFlux

   Call SNESGetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)

   Call DMMeshGetSectionReal(mesh,'default',rhsSec,ierr);CHKERRQ(ierr)
   Call DMMeshGetSectionReal(mesh,'default',fSec,ierr);CHKERRQ(ierr)
   Call DMMeshCreateGlobalScatter(mesh,rhsSec,ScatterSecToVec,ierr);CHKERRQ(ierr)

   Call SectionRealSet(rhsSec,0.0_Kr,ierr);CHKERRQ(ierr)
   Call VecSet(rhs,0.0_Kr,ierr);CHKERRQ(ierr)
   Call VecDuplicate(flux,modifiedFLux,ierr);CHKERRQ(ierr)

   Allocate(elem(numCell))
   Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
   Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
      Call PetscBagGetDataPoissonCellSetProperties(PoissonCtx%CellSetPropertiesBag(set),cellSetProperties,ierr);CHKERRQ(ierr)

      !!! Modified flux is flux + surfaceThermalConductivity * refTemp
      Call VecCopy(flux,modifiedFlux,ierr);CHKERRQ(ierr)
      Call VecAXPY(modifiedFlux,cellSetProperties%SurfaceThermalConductivity,refTemp,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(fSec,ScatterSecToVec,SCATTER_REVERSE,modifiedFlux,ierr);CHKERRQ(ierr)
      
      elemType = MEF90_knownElements(cellSetProperties%ElemTypeShortID)
      QuadratureOrder = elemType%order * 2
      Call MEF90_ElementCreate(mesh,setIS,elem,QuadratureOrder,CellSetProperties%ElemTypeShortID,ierr);CHKERRQ(ierr)

      Call MEF90_DiffusionRHSSet(rhssec,mesh,fsec,setIS,elem,elemType,ierr);CHKERRQ(ierr)

      Call MEF90_ElementDestroy(elem,ierr);CHKERRQ(ierr)
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
      Call PetscBagGetDataPoissonVertexSetProperties(PoissonCtx%VertexSetPropertiesBag(set),vertexSetProperties,ierr);CHKERRQ(ierr)
      If (vertexSetProperties%Has_BC) Then
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
End Subroutine SimplePoissonRHS

#undef __FUNCT__
#define __FUNCT__ "SimplePoissonEnergies_Cst"
!!!
!!!  
!!!  SimplePoissonEnergies_Cst:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine SimplePoissonEnergies_Cst(snesTemp,x,fluxScaling,PoissonCtx,energy,work,ierr)
   Type(SNES),Intent(IN)                           :: snesTemp
   Type(Vec),Intent(IN)                            :: x
   PetscReal,Intent(IN)                            :: fluxScaling
   Type(MEF90Ctx_Type),intent(IN)                  :: PoissonCtx
   PetscReal,Dimension(:),Pointer                  :: energy,work
   PetscErrorCode,Intent(OUT)                      :: ierr

   Type(IS)                                        :: CellSetGlobalIS,setIS,setISdof
   PetscInt,dimension(:),Pointer                   :: setID
   PetscInt                                        :: set,QuadratureOrder
   Type(MEF90_MATPROP),pointer                     :: matpropSet
   Type(PoissonCellSetProperties_Type),pointer     :: cellSetProperties
   PetscReal                                       :: myenergy,mywork
   Type(DM)                                        :: mesh
   Type(SectionReal)                               :: xSec
   Type(VecScatter)                                :: ScatterSecToVec
   Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer   :: elem
   Type(Element_Type)                              :: elemType

   energy = 0.0_Kr
   work = 0.0_Kr   
   Call SNESGetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)
   Call DMMeshGetSectionReal(mesh,'default',xSec,ierr);CHKERRQ(ierr)
   Call DMMeshCreateGlobalScatter(mesh,xSec,ScatterSecToVec,ierr);CHKERRQ(ierr)
   Call SectionRealToVec(xSec,ScatterSecToVec,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr)

   Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)   
   Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      myenergy = 0.0_Kr
      mywork   = 0.0_Kr
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)

      Call PetscBagGetDataMEF90_MatProp(PoissonCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataPoissonCellSetProperties(PoissonCtx%CellSetPropertiesBag(set),cellSetProperties,ierr);CHKERRQ(ierr)

      elemType = MEF90_knownElements(cellSetProperties%ElemTypeShortID)
      QuadratureOrder = elemType%order * 2
      Call MEF90_ElementCreate(mesh,setIS,elem,QuadratureOrder,CellSetProperties%ElemTypeShortID,ierr);CHKERRQ(ierr)
      !!! Assembly part of the residual coming from the bilinear form on the blocks of 
      !!! codimension 0
      If (MEF90_knownElements(cellSetProperties%ElemTypeShortID)%coDim == 0) Then
         Call MEF90_DiffusionEnergySet(myenergy,xSec,mesh,matpropSet%ThermalConductivity,0.0_Kr,setIS,elem,elemType,ierr)
      End If
      Call MPI_AllReduce(myenergy,energy(set),1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)

      If (cellSetProperties%Flux /= 0.0_Kr) Then
         Call MEF90_DiffusionWorkSet(mywork,xSec,mesh,fluxScaling*cellSetProperties%Flux,setIS,elem,elemType,ierr)
      End If
      Call MPI_AllReduce(mywork,work(set),1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)
      Call MEF90_ElementDestroy(elem,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   End Do     
   Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)

   Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
   Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
End Subroutine SimplePoissonEnergies_Cst

#undef __FUNCT__
#define __FUNCT__ "SimplePoissonEnergies"
!!!
!!!  
!!!  SimplePoissonEnergies:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine SimplePoissonEnergies(snesTemp,x,flux,PoissonCtx,energy,work,ierr)
   Type(SNES),Intent(IN)                           :: snesTemp
   Type(Vec),Intent(IN)                            :: x
   Type(Vec),Intent(IN)                            :: flux
   Type(MEF90Ctx_Type),intent(IN)                  :: PoissonCtx
   PetscReal,Dimension(:),Pointer                  :: energy,work
   PetscErrorCode,Intent(OUT)                      :: ierr

   Type(IS)                                        :: CellSetGlobalIS,setIS,setISdof
   PetscInt,dimension(:),Pointer                   :: setID
   PetscInt                                        :: set,QuadratureOrder
   Type(MEF90_MATPROP),pointer                     :: matpropSet
   Type(PoissonCellSetProperties_Type),pointer     :: cellSetProperties
   PetscReal                                       :: myenergy,mywork
   Type(DM)                                        :: mesh
   Type(SectionReal)                               :: xSec,fluxSec
   Type(VecScatter)                                :: ScatterSecToVec
   Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer   :: elem
   Type(Element_Type)                              :: elemType

   energy = 0.0_Kr
   work = 0.0_Kr   
   Call SNESGetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)
   Call DMMeshGetSectionReal(mesh,'default',xSec,ierr);CHKERRQ(ierr)
   Call DMMeshGetSectionReal(mesh,'default',fluxSec,ierr);CHKERRQ(ierr)
   Call DMMeshCreateGlobalScatter(mesh,xSec,ScatterSecToVec,ierr);CHKERRQ(ierr)
   Call SectionRealToVec(xSec,ScatterSecToVec,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr)
   Call SectionRealToVec(fluxSec,ScatterSecToVec,SCATTER_REVERSE,flux,ierr);CHKERRQ(ierr)

   Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)   
   Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      myenergy = 0.0_Kr
      mywork   = 0.0_Kr
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)

      Call PetscBagGetDataMEF90_MatProp(PoissonCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataPoissonCellSetProperties(PoissonCtx%CellSetPropertiesBag(set),cellSetProperties,ierr);CHKERRQ(ierr)

      elemType = MEF90_knownElements(cellSetProperties%ElemTypeShortID)
      QuadratureOrder = elemType%order * 2
      Call MEF90_ElementCreate(mesh,setIS,elem,QuadratureOrder,CellSetProperties%ElemTypeShortID,ierr);CHKERRQ(ierr)
      !!! Assembly part of the residual coming from the bilinear form on the blocks of 
      !!! codimension 0
      If (MEF90_knownElements(cellSetProperties%ElemTypeShortID)%coDim == 0) Then
         Call MEF90_DiffusionEnergySet(myenergy,xSec,mesh,matpropSet%ThermalConductivity,0.0_Kr,setIS,elem,elemType,ierr)
      End If
      Call MPI_AllReduce(myenergy,energy(set),1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)

      Call MEF90_DiffusionWorkSet(mywork,xSec,mesh,fluxSec,setIS,elem,elemType,ierr)

      Call MPI_AllReduce(mywork,work(set),1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)
      Call MEF90_ElementDestroy(elem,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   End Do     
   Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)

   Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
   Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
End Subroutine SimplePoissonEnergies
End Module M_POISSONASSEMBLY