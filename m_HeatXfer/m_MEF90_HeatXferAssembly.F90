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
   
      Type(SectionReal)                                  :: xSec,residualSec,modifiedFluxSec
      Type(SectionReal)                                  :: fluxSec,externalTemperatureSec,boundaryTemperatureSec
      PetscReal,Dimension(:),Pointer                     :: xPtr,residualPtr,modifiedFluxPtr
      PetscReal,Dimension(:),Pointer                     :: fluxPtr,externalTemperaturePtr,boundaryTemperaturePtr
      Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),pointer                        :: matpropSet
      Type(MEF90HeatXferCellSetOptions_Type),pointer     :: cellSetOptions
      Type(MEF90HeatXferVertexSetOptions_Type),pointer   :: vertexSetOptions
      Type(VecScatter)                                   :: ScatterSecToVec,ScatterSecToVecCell
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elem
      Type(MEF90Element_Type)                            :: elemType
      PetscInt                                           :: dof

      !!! Create dof-based sections
      Call DMMeshGetSectionReal(MEF90HeatXferCtx%DM,'default',xSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(xSec,residualSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(xSec,boundaryTemperatureSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90HeatXferCtx%DM,xSec,ScatterSecToVec,ierr);CHKERRQ(ierr)

      !!! Create cell based sections, and allocate required pointers
      Allocate(modifiedFluxPtr(1))
      Call DMMeshGetSectionReal(MEF90HeatXferCtx%CellDM,'default',modifiedFluxSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90HeatXferCtx%CellDM,modifiedFluxSec,ScatterSecToVecCell,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(modifiedFluxSec,fluxSec,ierr);CHKERRQ(ierr)
   
      Call SectionRealDuplicate(modifiedFluxSec,externalTemperatureSec,ierr);CHKERRQ(ierr)
      Call SectionRealSet(modifiedFluxSec,0.0_Kr,ierr);CHKERRQ(ierr)
   
      !!! Scatter data from Vec to Sec, or initialize
      Call SectionRealToVec(xSec,ScatterSecToVec,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr)
        
      Call SectionRealSet(modifiedFluxSec,0.0_Kr,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(fluxSec,ScatterSecToVecCell,SCATTER_REVERSE,MEF90HeatXferCtx%flux,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(externalTemperatureSec,ScatterSecToVecCell,SCATTER_REVERSE,MEF90HeatXferCtx%externalTemperature,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(boundaryTemperatureSec,ScatterSecToVec,SCATTER_REVERSE,MEF90HeatXferCtx%boundaryTemperature,ierr);CHKERRQ(ierr)

      Call SectionRealSet(residualSec,0.0_Kr,ierr);CHKERRQ(ierr)
      Call VecSet(residual,0.0_Kr,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call DMMeshGetStratumIS(MEF90HeatXferCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call ISGetIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
         
         Call PetscBagGetDataMEF90MatProp(MEF90HeatXferCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
   
         elemType = MEF90_knownElements(cellSetOptions%ElemTypeShortID)
         QuadratureOrder = elemType%order * 2
         Call MEF90Element_Create(MEF90HeatXferCtx%DM,setIS,elem,QuadratureOrder,CellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)
         Call MEF90DiffusionOperatorSet(residualSec,MEF90HeatXferCtx%DM,xSec,setIS,matpropSet%ThermalConductivity,cellSetOptions%SurfaceThermalConductivity,elem,elemType,ierr);CHKERRQ(ierr)

         !!! Modified flux is flux + surfaceThermalConductivity * refTemp      
         !!! I _could_ use a SecAXPY, but this would summ all values at all cells for each block
         !!! I _could_ also create Sections restricted to the cell set only
         !!! But this is going away with sections anyway...
         Do dof = 1,size(setIdx)
            Call SectionRealRestrict(fluxSec,setIdx(dof),fluxPtr,ierr);CHKERRQ(ierr)
            Call SectionRealRestrict(externalTemperatureSec,setIdx(dof),externalTemperaturePtr,ierr);CHKERRQ(ierr)
            modifiedFluxPtr =  -fluxPtr - cellSetOptions%SurfaceThermalConductivity * externalTemperaturePtr
            Call SectionRealUpdate(modifiedFluxSec,setIdx(dof),modifiedFluxPtr,INSERT_VALUES,ierr);CHKERRQ(ierr)
            Call SectionRealRestore(externalTemperatureSec,setIdx(dof),externalTemperaturePtr,ierr);CHKERRQ(ierr)
            Call SectionRealRestore(fluxSec,setIdx(dof),fluxPtr,ierr);CHKERRQ(ierr)
         End Do
         Call MEF90DiffusionRHSSetCell(residualsec,MEF90HeatXferCtx%DM,modifiedFluxSec,setIS,elem,elemType,ierr);CHKERRQ(ierr)
          
         Call MEF90Element_Destroy(elem,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
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

      Call SectionRealToVec(temperatureSec,ScatterSecToVec,SCATTER_REVERSE,temperatureVec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(externalTemperatureSec,ScatterSecToVecCell,SCATTER_REVERSE,MEF90HeatXferCtx%externalTemperature,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(fluxSec,ScatterSecToVecCell,SCATTER_REVERSE,MEF90HeatXferCtx%flux,ierr);CHKERRQ(ierr)   
      
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
