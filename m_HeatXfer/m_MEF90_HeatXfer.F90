#include "../MEF90/mef90.inc"
Module m_MEF90_HeatXfer
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_HeatXferCtx
   Implicit none

   !Interface  m_MEF90_HeatXferOperatorAssembly
   !   Module Procedure m_MEF90_HeatXferOperatorAssembly2D,m_MEF90_HeatXferOperatorAssembly3D
   !End Interface  m_MEF90_HeatXferOperatorAssembly
   
Contains

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferInterpolateFields"
!!!
!!!  
!!!  MEF90HeatXferInterpolateFields:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferInterpolateFields(time,MEF90HeatXferCtx,ierr)
      PetscErrorCode,Intent(OUT)                            :: ierr
      PetscReal,Intent(IN)                                  :: time
      Type(MEF90HeatXferCtx_Type),intent(IN)                :: MEF90HeatXferCtx
   
      Type(MEF90HeatXferGlobalOptions_Type),pointer         :: MEF90HeatXferGlobalOptions
      Type(MEF90HeatXferVertexSetOptions_Type),pointer      :: vertexSetOptions
      Type(MEF90HeatXferCellSetOptions_Type),pointer        :: vertexSetOptions
      Type(IS)                                              :: VertexSetGlobalIS,cellSetGlobalIS,setIS,setISdof
      PetscInt,dimension(:),Pointer                         :: setID
      PetscInt,Dimension(:),Pointer                         :: setIdx
      PetscInt                                              :: set
      PetscReal,Dimension(:),Pointer                        :: val
      Type(SectionReal)                                     :: xSec
      
      Call PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr);CHKERRQ(ierr)
      Call DMMeshGetSectionReal(MEF90HeatXferCtx%DM,'default',xSec,ierr);CHKERRQ(ierr)

      !!! BoundaryTemperature is define at vertex sets
      Call VecSet(MEF90HeatXferCtx%BoundaryTemperature,0.0_Kr,ierr);CHKERRQ(ierr)
      Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90HeatXferCtxVertexSetOptions(MEF90HeatXferCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
         If (vertexSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(MEF90HeatXferCtx%DM,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
            Call DMMeshISCreateISglobaldof(MEF90HeatXferCtx%DM,xSec,setIS,0,setISdof,ierr);CHKERRQ(ierr)
            Call ISGetIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
            Allocate(val(size(setIdx)),stat=ierr)
            Select Case (MEF90HeatXferGlobalOptions%boundaryTempScaling) 
               Case(MEF90Scaling_CST)
                  val = vertexSetOptions%boundaryTemp
               Case(MEF90Scaling_Linear)
                  val = vertexSetOptions%boundaryTemp * time
               Case(MEF90Scaling_File)
                  !!! Interpolate boundary displacement
                  Call VecGetValues(MEF90HeatXferCtx%boundaryTemperature,size(setIdx),setIdx,val,ierr);CHKERRQ(ierr)
            End Select
            Call VecSetValues(MEF90HeatXferCtx%BoundaryTemperature,size(setIdx),setIdx,val,INSERT_VALUES,ierr);CHKERRQ(ierr)
            DeAllocate(val)
            Call ISRestoreIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
         End If !vertexSetOptions%Has_BC
      End Do
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call VecAssemblyBegin(MEF90HeatXferCtx%BoundaryTemperature,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(MEF90HeatXferCtx%BoundaryTemperature,ierr);CHKERRQ(ierr)
      
      !!! Flux and external temperature are constant on each cell set
      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90HeatXferInterpolateFields

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferSetBoundaryTemperature"
!!!
!!!  
!!!  MEF90HeatXferSetBoundaryTemperature:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferSetBoundaryTemperature(x,time,MEF90HeatXferCtx,ierr)
      Type(Vec),Intent(IN)                                  :: x
      PetscReal,Intent(IN)                                  :: time
      Type(MEF90HeatXferCtx_Type),intent(IN)                :: MEF90HeatXferCtx
      PetscErrorCode,Intent(OUT)                            :: ierr
   
      Type(MEF90HeatXferGlobalOptions_Type),pointer         :: MEF90HeatXferGlobalOptions
      Type(MEF90HeatXferVertexSetOptions_Type),pointer      :: vertexSetOptions
      Type(IS)                                              :: VertexSetGlobalIS,setIS,setISdof
      PetscInt,dimension(:),Pointer                         :: setID
      PetscInt,Dimension(:),Pointer                         :: setIdx
      PetscInt                                              :: set
      PetscReal,Dimension(:),Pointer                        :: val
      Type(SectionReal)                                     :: xSec
      
      Call PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr);CHKERRQ(ierr)
      Call DMMeshGetSectionReal(MEF90HeatXferCtx%DM,'default',xSec,ierr);CHKERRQ(ierr)
   
      Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90HeatXferCtxVertexSetOptions(MEF90HeatXferCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
         If (vertexSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(MEF90HeatXferCtx%DM,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
            Call DMMeshISCreateISglobaldof(MEF90HeatXferCtx%DM,xSec,setIS,0,setISdof,ierr);CHKERRQ(ierr)
            Call ISGetIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
            Allocate(val(size(setIdx)),stat=ierr)
            Call VecGetValues(MEF90HeatXferCtx%boundaryTemperature,size(setIdx),setIdx,val,ierr);CHKERRQ(ierr)
            Call VecSetValues(x,size(setIdx),setIdx,val,INSERT_VALUES,ierr);CHKERRQ(ierr)
            DeAllocate(val)
            Call ISRestoreIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
         End If !vertexSetOptions%Has_BC
      End Do
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
   
      Call VecAssemblyBegin(x,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(x,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90HeatXferSetBoundaryTemperature
End Module m_MEF90_HeatXfer