#include "../MEF90/mef90.inc"
Module m_MEF90_HeatXfer
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_HeatXferCtx
   Implicit none
   Private
   Public MEF90HeatXferGetTransients
   Public MEF90HeatXferSetBoundaryTemperature
   
Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferGetTransients"
!!!
!!!  
!!!  MEF90HeatXferGetTransients:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferGetTransients(MEF90HeatXferCtx,step,time,ierr)
      Type(MEF90HeatXferCtx_Type),Intent(INOUT)       :: MEF90HeatXferCtx
      PetscInt,Intent(IN)                             :: step
      PetscReal,Intent(IN)                            :: time
      PetscErrorCode,Intent(OUT)                      :: ierr
   
      Type(MEF90HeatXferGlobalOptions_Type),pointer   :: MEF90HeatXferGlobalOptions
      Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90GlobalOptions

      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90HeatXferCtx%MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr);CHKERRQ(ierr)

      If (step > 1) Then
         !!! This copy is not necessary if the scaling of each field is CST... 
         Call VecCopy(MEF90HeatXferCtx%fluxTarget,MEF90HeatXferCtx%fluxPrevious,ierr);CHKERRQ(ierr)
         Call VecCopy(MEF90HeatXferCtx%externalTemperatureTarget,MEF90HeatXferCtx%externalTemperaturePrevious,ierr);CHKERRQ(ierr)
         Call VecCopy(MEF90HeatXferCtx%boundaryTemperatureTarget,MEF90HeatXferCtx%boundaryTemperaturePrevious,ierr);CHKERRQ(ierr)
         MEF90HeatXferCtx%timePrevious = MEF90HeatXferCtx%timeTarget
      End If
      MEF90HeatXferCtx%timeTarget = time
      
      Select case (MEF90HeatXferGlobalOptions%fluxScaling)
         Case (MEF90Scaling_File)
            Call VecLoadExodusCell(MEF90HeatXferCtx%cellDM,MEF90HeatXferCtx%fluxTarget,MEF90HeatXferCtx%MEF90Ctx%IOcomm, &
                                   MEF90HeatXferCtx%MEF90Ctx%fileExoUnit,step,MEF90HeatXferGlobalOptions%fluxOffset,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_Linear)
            Call MEF90HeatXferGetfluxCst(MEF90HeatXferCtx%fluxTarget,MEF90HeatXferCtx,ierr)
            Call VecScale(MEF90HeatXferCtx%fluxTarget,time,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_CST)
            Call MEF90HeatXferGetfluxCst(MEF90HeatXferCtx%fluxTarget,MEF90HeatXferCtx,ierr)
      End Select
      Select case (MEF90HeatXferGlobalOptions%externalTempScaling)
         Case (MEF90Scaling_File)
            Call VecLoadExodusCell(MEF90HeatXferCtx%cellDM,MEF90HeatXferCtx%externalTemperatureTarget,MEF90HeatXferCtx%MEF90Ctx%IOcomm, &
                                   MEF90HeatXferCtx%MEF90Ctx%fileExoUnit,step,MEF90HeatXferGlobalOptions%externalTempOffset,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_Linear)
            Call MEF90HeatXferGetexternalTemperatureCst(MEF90HeatXferCtx%externalTemperatureTarget,MEF90HeatXferCtx,ierr)
            Call VecScale(MEF90HeatXferCtx%externalTemperatureTarget,time,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_CST)
            Call MEF90HeatXferGetexternalTemperatureCst(MEF90HeatXferCtx%externalTemperatureTarget,MEF90HeatXferCtx,ierr)
      End Select
      Select case (MEF90HeatXferGlobalOptions%boundaryTempScaling)
         Case (MEF90Scaling_File)
            Call VecLoadExodusVertex(MEF90HeatXferCtx%cellDM,MEF90HeatXferCtx%boundaryTemperatureTarget,MEF90HeatXferCtx%MEF90Ctx%IOcomm, &
                                     MEF90HeatXferCtx%MEF90Ctx%fileExoUnit,step,MEF90HeatXferGlobalOptions%boundaryTempOffset,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_Linear)
            Call MEF90HeatXferGetboundaryTemperatureCst(MEF90HeatXferCtx%boundaryTemperatureTarget,MEF90HeatXferCtx,ierr)
            Call VecScale(MEF90HeatXferCtx%boundaryTemperatureTarget,time,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_CST)
            Call MEF90HeatXferGetboundaryTemperatureCst(MEF90HeatXferCtx%boundaryTemperatureTarget,MEF90HeatXferCtx,ierr)
      End Select
      If (step == 1) Then
         Call VecCopy(MEF90HeatXferCtx%fluxTarget,MEF90HeatXferCtx%fluxPrevious,ierr);CHKERRQ(ierr)
         Call VecCopy(MEF90HeatXferCtx%externalTemperatureTarget,MEF90HeatXferCtx%externalTemperaturePrevious,ierr);CHKERRQ(ierr)
         Call VecCopy(MEF90HeatXferCtx%boundaryTemperatureTarget,MEF90HeatXferCtx%boundaryTemperaturePrevious,ierr);CHKERRQ(ierr)
         MEF90HeatXferCtx%timePrevious = MEF90HeatXferCtx%timeTarget
      End If
   End Subroutine MEF90HeatXferGetTransients
   
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferGetFluxCst"
!!!
!!!  
!!!  MEF90HeatXferGetFluxCst:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferGetFluxCst(x,MEF90HeatXferCtx,ierr)
      Type(Vec),Intent(IN)                            :: x
      Type(MEF90HeatXferCtx_Type),Intent(IN)          :: MEF90HeatXferCtx
      PetscErrorCode,Intent(OUT)                      :: ierr
   
      Type(MEF90HeatXferGlobalOptions_Type),pointer   :: MEF90HeatXferGlobalOptions
      Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90GlobalOptions
      Type(MEF90HeatXferCellSetOptions_Type),pointer  :: cellSetOptions
      Type(IS)                                        :: cellSetGlobalIS,setIS,setISdof
      PetscInt,dimension(:),Pointer                   :: setID
      PetscInt,Dimension(:),Pointer                   :: setIdx
      PetscInt                                        :: set
      PetscReal,Dimension(:),Pointer                  :: val
      
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90HeatXferCtx%MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr);CHKERRQ(ierr)
      
      !!! flux is cell-centered
      Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90HeatXferCtx%cellDM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call DMMeshISCreateISglobaldof(MEF90HeatXferCtx%cellDM,setIS,0,setISdof,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
         Allocate(val(size(setIdx)),stat=ierr)
         val = cellSetOptions%flux
         Call VecSetValues(x,size(setIdx),setIdx,val,INSERT_VALUES,ierr);CHKERRQ(ierr)
         DeAllocate(val)
         Call ISRestoreIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(cellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(cellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call VecAssemblyBegin(x,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(x,ierr);CHKERRQ(ierr)
   End Subroutine MEF90HeatXferGetFluxCst

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferGetexternalTemperatureCst"
!!!
!!!  
!!!  MEF90HeatXferGetexternalTemperatureCst:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferGetexternalTemperatureCst(x,MEF90HeatXferCtx,ierr)
      Type(Vec),Intent(IN)                            :: x
      Type(MEF90HeatXferCtx_Type),Intent(IN)          :: MEF90HeatXferCtx
      PetscErrorCode,Intent(OUT)                      :: ierr
   
      Type(MEF90HeatXferGlobalOptions_Type),pointer   :: MEF90HeatXferGlobalOptions
      Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90GlobalOptions
      Type(MEF90HeatXferCellSetOptions_Type),pointer  :: cellSetOptions
      Type(IS)                                        :: cellSetGlobalIS,setIS,setISdof
      PetscInt,dimension(:),Pointer                   :: setID
      PetscInt,Dimension(:),Pointer                   :: setIdx
      PetscInt                                        :: set
      PetscReal,Dimension(:),Pointer                  :: val
      !Type(SectionReal)                               :: xSec
      
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90HeatXferCtx%MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr);CHKERRQ(ierr)
      
      !!! externalTemperature is cell-centered
      Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90HeatXferCtx%cellDM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call DMMeshISCreateISglobaldof(MEF90HeatXferCtx%cellDM,setIS,0,setISdof,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
         Allocate(val(size(setIdx)),stat=ierr)
         val = cellSetOptions%externalTemp
         Call VecSetValues(x,size(setIdx),setIdx,val,INSERT_VALUES,ierr);CHKERRQ(ierr)
         DeAllocate(val)
         Call ISRestoreIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(cellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(cellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call VecAssemblyBegin(x,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(x,ierr);CHKERRQ(ierr)
   End Subroutine MEF90HeatXferGetexternalTemperatureCst

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferGetboundaryTemperatureCst"
!!!
!!!  
!!!  MEF90HeatXferGetboundaryTemperatureCst:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferGetboundaryTemperatureCst(x,MEF90HeatXferCtx,ierr)
      Type(Vec),Intent(IN)                               :: x
      Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(MEF90HeatXferGlobalOptions_Type),pointer      :: MEF90HeatXferGlobalOptions
      Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions
      Type(MEF90HeatXferVertexSetOptions_Type),pointer   :: vertexSetOptions
      Type(IS)                                           :: VertexSetGlobalIS,setIS,setISdof
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt,Dimension(:),Pointer                      :: setIdx
      PetscInt                                           :: set
      PetscReal,Dimension(:),Pointer                     :: val
      
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90HeatXferCtx%MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr);CHKERRQ(ierr)
      
      !!! boundaryTemperature is Vertex-centered
      Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90HeatXferCtxVertexSetOptions(MEF90HeatXferCtx%VertexSetOptionsBag(set),VertexSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90HeatXferCtx%DM,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call DMMeshISCreateISglobaldof(MEF90HeatXferCtx%DM,setIS,0,setISdof,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
         Allocate(val(size(setIdx)),stat=ierr)
         val = VertexSetOptions%boundaryTemp
         Call VecSetValues(x,size(setIdx),setIdx,val,INSERT_VALUES,ierr);CHKERRQ(ierr)
         DeAllocate(val)
         Call ISRestoreIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call VecAssemblyBegin(x,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(x,ierr);CHKERRQ(ierr)
   End Subroutine MEF90HeatXferGetboundaryTemperatureCst

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferSetBoundaryTemperature"
!!!
!!!  
!!!  MEF90HeatXferSetBoundaryTemperature:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
!!! FOLLOW THE LINES OF MEF90HeatXferOperator and use VecSetValues instead of SectionRealUpdate + Scatter
   Subroutine MEF90HeatXferSetBoundaryTemperature(x,t,MEF90HeatXferCtx,ierr)
      Type(Vec),Intent(IN)                               :: x
      PetscReal,Intent(IN)                               :: t
      Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(MEF90HeatXferGlobalOptions_Type),pointer      :: MEF90HeatXferGlobalOptions
      Type(MEF90HeatXferVertexSetOptions_Type),pointer   :: vertexSetOptions
      Type(IS)                                           :: VertexSetGlobalIS,setIS
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt,Dimension(:),Pointer                      :: setIdx
      PetscInt                                           :: set
      PetscReal                                          :: scaling
      Type(SectionReal)                                  :: boundaryTemperaturePreviousSec,boundaryTemperatureTargetSec,xSec
      PetscReal,Dimension(:),Pointer                     :: boundaryTemperaturePreviousPtr,boundaryTemperatureTargetPtr,xPtr
      Type(VecScatter)                                   :: ScatterSecToVec
      PetscInt                                           :: dof
      
      Call PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr);CHKERRQ(ierr)
      
      !!! Get Sections and VecScatter
      Call DMMeshGetSectionReal(MEF90HeatXferCtx%DM,'default',xSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(xSec,boundaryTemperaturePreviousSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(xSec,boundaryTemperatureTargetSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90HeatXferCtx%DM,xSec,ScatterSecToVec,ierr);CHKERRQ(ierr)

      !!! Scatter from the Vecs to the sections
      Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(xSec,ScatterSecToVec,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(boundaryTemperaturePreviousSec,ScatterSecToVec,SCATTER_REVERSE,MEF90HeatXferCtx%boundaryTemperaturePrevious,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(boundaryTemperatureTargetSec,ScatterSecToVec,SCATTER_REVERSE,MEF90HeatXferCtx%boundaryTemperatureTarget,ierr);CHKERRQ(ierr)

      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Allocate(xPtr(1))
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90HeatXferCtxVertexSetOptions(MEF90HeatXferCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
         If (vertexSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(MEF90HeatXferCtx%DM,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
            Call ISGetIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
            If (t == MEF90HeatXferCtx%timePrevious) Then
               Do dof = 1, size(setIdx)
                  Call SectionRealRestrict(boundaryTemperaturePreviousSec,setIdx(dof),boundaryTemperaturePreviousPtr,ierr);CHKERRQ(ierr)
                  xPtr = boundaryTemperaturePreviousPtr
                  Call SectionRealUpdate(xSec,setIdx(dof),xPtr,INSERT_VALUES,ierr);CHKERRQ(ierr)
                  Call SectionRealRestore(boundaryTemperaturePreviousSec,setIdx(dof),boundaryTemperaturePreviousPtr,ierr);CHKERRQ(ierr)                  
               End Do
            Else If (t == MEF90HeatXferCtx%timeTarget) Then
               Do dof = 1, size(setIdx)
                  Call SectionRealRestrict(boundaryTemperatureTargetSec,setIdx(dof),boundaryTemperatureTargetPtr,ierr);CHKERRQ(ierr)
                  xPtr = boundaryTemperatureTargetPtr
                  Call SectionRealUpdate(xSec,setIdx(dof),xPtr,INSERT_VALUES,ierr);CHKERRQ(ierr)
                  Call SectionRealRestore(boundaryTemperatureTargetSec,setIdx(dof),boundaryTemperatureTargetPtr,ierr);CHKERRQ(ierr)                  
               End Do
            Else
               Do dof = 1, size(setIdx)
                  Call SectionRealRestrict(boundaryTemperaturePreviousSec,setIdx(dof),boundaryTemperaturePreviousPtr,ierr);CHKERRQ(ierr)
                  Call SectionRealRestrict(boundaryTemperatureTargetSec,setIdx(dof),boundaryTemperatureTargetPtr,ierr);CHKERRQ(ierr)
                  xPtr = (boundaryTemperaturePreviousPtr * (MEF90HeatXferCtx%timeTarget - t) +                                      &
                          boundaryTemperatureTargetPtr * (MEF90HeatXferCtx%timePrevious - t)) /                                     &
                         (MEF90HeatXferCtx%timeTarget - MEF90HeatXferCtx%timePrevious)
                  Call SectionRealUpdate(xSec,setIdx(dof),xPtr,INSERT_VALUES,ierr);CHKERRQ(ierr)
                  Call SectionRealRestore(boundaryTemperatureTargetSec,setIdx(dof),boundaryTemperatureTargetPtr,ierr);CHKERRQ(ierr)                  
                  Call SectionRealRestore(boundaryTemperaturePreviousSec,setIdx(dof),boundaryTemperaturePreviousPtr,ierr);CHKERRQ(ierr)
               End Do
            End If
            Call ISRestoreIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If !vertexSetOptions%Has_BC
      End Do
      DeAllocate(xPtr)
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(xSec,ScatterSecToVec,SCATTER_FORWARD,x,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(boundaryTemperatureTargetSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(boundaryTemperaturePreviousSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90HeatXferSetBoundaryTemperature
End Module m_MEF90_HeatXfer