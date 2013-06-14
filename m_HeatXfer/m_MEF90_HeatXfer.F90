#include "../MEF90/mef90.inc"
Module m_MEF90_HeatXfer
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_HeatXferCtx
   Use m_MEF90_HeatXferAssembly2D!, &
      !MEF90HeatXferEnergy2D => MEF90HeatXferEnergy, &
      !MEF90HeatXferOperator2D => MEF90HeatXferOperator, &
      !MEF90HeatXferBilinearForm2D => MEF90HeatXferBilinearForm_private
   Use m_MEF90_HeatXferAssembly3D!, &
      !MEF90HeatXferEnergy3D => MEF90HeatXFerEnergy, &
      !MEF90HeatXferOperator3D => MEF90HeatXferOperator, &
      !MEF90HeatXferBilinearForm3D => MEF90HeatXferBilinearForm_private
      !!! For some reason, this trick seems to confuse intel 13.0...
   !Private
   !Public MEF90HeatXferSetTransients
   !Public MEF90HeatXferOperator
   !Public MEF90HeatXferBilinearForm
   !Public MEF90HeatXferEnergy
   !Public :: MEF90HeatXferCtx_Type
   !Public :: MEF90HeatXferGlobalOptions_Type
   !Public :: MEF90HeatXferCellSetOptions_Type
   !Public :: MEF90HeatXferVertexSetOptions_Type
   Implicit none
Contains

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferSetTransients"
!!!
!!!  
!!!  MEF90HeatXferSetTransients:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferSetTransients(MEF90HeatXferCtx,step,time,ierr)
      Type(MEF90HeatXferCtx_Type),Intent(INOUT)       :: MEF90HeatXferCtx
      PetscInt,Intent(IN)                             :: step
      PetscReal,Intent(IN)                            :: time
      PetscErrorCode,Intent(OUT)                      :: ierr
   
      Type(MEF90HeatXferGlobalOptions_Type),pointer   :: MEF90HeatXferGlobalOptions
      Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90GlobalOptions

      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90HeatXferCtx%MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr);CHKERRQ(ierr)

      Select case (MEF90HeatXferGlobalOptions%fluxScaling)
         Case (MEF90Scaling_File)
            Call VecLoadExodusCell(MEF90HeatXferCtx%cellDM,MEF90HeatXferCtx%flux,MEF90HeatXferCtx%MEF90Ctx%IOcomm, &
                                   MEF90HeatXferCtx%MEF90Ctx%fileExoUnit,step,MEF90HeatXferGlobalOptions%fluxOffset,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_Linear)
            Call MEF90HeatXferSetfluxCst(MEF90HeatXferCtx%flux,MEF90HeatXferCtx,ierr)
            Call VecScale(MEF90HeatXferCtx%flux,time,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_CST)
            Call MEF90HeatXferSetfluxCst(MEF90HeatXferCtx%flux,MEF90HeatXferCtx,ierr)
      End Select
      Select case (MEF90HeatXferGlobalOptions%externalTempScaling)
         Case (MEF90Scaling_File)
            Call VecLoadExodusCell(MEF90HeatXferCtx%cellDM,MEF90HeatXferCtx%externalTemperature,MEF90HeatXferCtx%MEF90Ctx%IOcomm, &
                                   MEF90HeatXferCtx%MEF90Ctx%fileExoUnit,step,MEF90HeatXferGlobalOptions%externalTempOffset,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_Linear)
            Call MEF90HeatXferSetexternalTemperatureCst(MEF90HeatXferCtx%externalTemperature,MEF90HeatXferCtx,ierr)
            Call VecScale(MEF90HeatXferCtx%externalTemperature,time,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_CST)
            Call MEF90HeatXferSetexternalTemperatureCst(MEF90HeatXferCtx%externalTemperature,MEF90HeatXferCtx,ierr)
      End Select
      Select case (MEF90HeatXferGlobalOptions%boundaryTempScaling)
         Case (MEF90Scaling_File)
            Call VecLoadExodusVertex(MEF90HeatXferCtx%cellDM,MEF90HeatXferCtx%boundaryTemperature,MEF90HeatXferCtx%MEF90Ctx%IOcomm, &
                                     MEF90HeatXferCtx%MEF90Ctx%fileExoUnit,step,MEF90HeatXferGlobalOptions%boundaryTempOffset,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_Linear)
            Call MEF90HeatXferSetboundaryTemperatureCst(MEF90HeatXferCtx%boundaryTemperature,MEF90HeatXferCtx,ierr)
            Call VecScale(MEF90HeatXferCtx%boundaryTemperature,time,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_CST)
            Call MEF90HeatXferSetboundaryTemperatureCst(MEF90HeatXferCtx%boundaryTemperature,MEF90HeatXferCtx,ierr)
      End Select
   End Subroutine MEF90HeatXferSetTransients
   
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferSetFluxCst"
!!!
!!!  
!!!  MEF90HeatXferSetFluxCst:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferSetFluxCst(x,MEF90HeatXferCtx,ierr)
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
   End Subroutine MEF90HeatXferSetFluxCst

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferSetexternalTemperatureCst"
!!!
!!!  
!!!  MEF90HeatXferSetexternalTemperatureCst:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferSetexternalTemperatureCst(x,MEF90HeatXferCtx,ierr)
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
   End Subroutine MEF90HeatXferSetexternalTemperatureCst

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferSetboundaryTemperatureCst"
!!!
!!!  
!!!  MEF90HeatXferSetboundaryTemperatureCst:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferSetboundaryTemperatureCst(x,MEF90HeatXferCtx,ierr)
      Type(Vec),Intent(IN)                               :: x
      Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(MEF90HeatXferGlobalOptions_Type),pointer      :: MEF90HeatXferGlobalOptions
      Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions
      Type(MEF90HeatXferVertexSetOptions_Type),pointer   :: vertexSetOptions
      Type(IS)                                           :: VertexSetGlobalIS,setIS,setISdof,bcIS
      Type(MEF90HeatXferCellSetOptions_Type),pointer     :: cellSetOptions
      Type(IS)                                           :: CellSetGlobalIS
      PetscInt,Dimension(:),Pointer                      :: setID
      PetscInt,Dimension(:),Pointer                      :: setIdx
      PetscInt                                           :: set
      PetscReal,Dimension(:),Pointer                     :: val
      PetscInt,Dimension(:),Pointer                      :: cone
      Type(MEF90Element_Type)                            :: elemType
      PetscInt                                           :: cell,dof
      
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90HeatXferCtx%MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr);CHKERRQ(ierr)
      
      !!! boundaryTemperature is Vertex-centered
      !!! We first set the boundary values inherited from cell sets, then that of vertex sets
      !!!
      Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

      !!!
      !!! We loop over all element twice. The first tim in order to assembly all not BC cell sets
      !!! In the second pass, we only update the BC where necessary
      !!! vertex cet BC are updated last, so that they override cell set BC
      !!!

      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)         
         If (cellSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(MEF90HeatXferCtx%dm,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
            Call MEF90_ISCreateCelltoVertex(MEF90HeatXferCtx%dm,PETSC_COMM_WORLD,setIS,bcIS,ierr)
            Call DMMeshISCreateISglobaldof(MEF90HeatXferCtx%DM,bcIS,0,setISdof,ierr);CHKERRQ(ierr)
            Call ISGetIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
            Allocate(val(size(setIdx)),stat=ierr)
            val = cellSetOptions%boundaryTemp
            Call VecSetValues(x,size(setIdx),setIdx,val,INSERT_VALUES,ierr);CHKERRQ(ierr)
            DeAllocate(val)
            Call ISRestoreIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
            Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
            Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If ! cellSetOptions%Has_BC
      End Do ! set
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
      
      !!! Vertex Sets
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
         Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call VecAssemblyBegin(x,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(x,ierr);CHKERRQ(ierr)
   End Subroutine MEF90HeatXferSetboundaryTemperatureCst

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferUpdateboundaryTemperature"
!!!
!!!  
!!!  MEF90HeatXferUpdateboundaryTemperature:
!!!  
!!!  (c) 2013 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferUpdateboundaryTemperature(x,MEF90HeatXferCtx,ierr)
      Type(Vec),Intent(IN)                               :: x
      Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
      PetscErrorCode,Intent(OUT)                         :: ierr

   
      Type(MEF90HeatXferGlobalOptions_Type),pointer      :: MEF90HeatXferGlobalOptions
      Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions
      Type(MEF90HeatXferVertexSetOptions_Type),pointer   :: vertexSetOptions
      Type(IS)                                           :: VertexSetGlobalIS,setIS,setISdof
      Type(MEF90HeatXferVertexSetOptions_Type),pointer   :: cellSetOptions
      Type(IS)                                           :: CellSetGlobalIS
      PetscInt,Dimension(:),Pointer                      :: setID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,dof,cell
      PetscReal,Dimension(:),Pointer                     :: boundaryTemperaturePtr,xPtr
      Type(SectionReal)                                  :: boundaryTemperatureSec
      Type(VecScatter)                                   :: ScatterSecToVec
      Type(MEF90Element_Type)                            :: elemType
      
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90HeatXferCtx%MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr);CHKERRQ(ierr)
      
      !!! boundaryTemperature is Vertex-centered
      Call DMMeshGetSectionReal(MEF90HeatXferCtx%DM,'default',boundaryTemperatureSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90HeatXferCtx%DM,boundaryTemperatureSec,ScatterSecToVec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(boundaryTemperatureSec,ScatterSecToVec,SCATTER_REVERSE,MEF90HeatXferCtx%boundaryTemperature,ierr);CHKERRQ(ierr)
      
      !!!
      !!! cell set temperature first, followed vertex sets
      !!!
      Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90HeatXferCtxVertexSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         If (cellSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(MEF90HeatXferCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
            Call ISGetIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
            Call DMMeshISCreateISglobaldof(MEF90HeatXferCtx%DM,setIS,0,setISdof,ierr);CHKERRQ(ierr)
            Call ISGetIndicesF90(setISdof,setdofIdx,ierr);CHKERRQ(ierr)
            Allocate(xPtr(size(setIdx)))
            Do dof = 1,size(setIdx)
               Call SectionRealRestrict(boundaryTemperatureSec,setIdx(dof),boundaryTemperaturePtr,ierr);CHKERRQ(ierr)
               xPtr(dof) = boundaryTemperaturePtr(1)
               Call SectionRealRestore(boundaryTemperatureSec,setIdx(dof),boundaryTemperaturePtr,ierr);CHKERRQ(ierr)
            End Do
            Call VecSetValues(x,size(setIdx),setdofIdx,xPtr,INSERT_VALUES,ierr);CHKERRQ(ierr)
            DeAllocate(xPtr)
            Call ISRestoreIndicesF90(setISdof,setdofIdx,ierr);CHKERRQ(ierr)
            Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
            Call ISRestoreIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If
      End Do
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)

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
            Allocate(xPtr(size(setIdx)))
            Do dof = 1, size(setIdx)
               Call SectionRealRestrict(boundaryTemperatureSec,setIdx(dof),boundaryTemperaturePtr,ierr);CHKERRQ(ierr)
               xPtr(dof) = boundaryTemperaturePtr(1)
               Call SectionRealRestore(boundaryTemperatureSec,setIdx(dof),boundaryTemperaturePtr,ierr);CHKERRQ(ierr)
            End Do
            Call VecSetValues(x,size(setIdx),setdofIdx,xPtr,INSERT_VALUES,ierr);CHKERRQ(ierr)
            DeAllocate(xPtr)
            Call ISRestoreIndicesF90(setISdof,setdofIdx,ierr);CHKERRQ(ierr)
            Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
            Call ISRestoreIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If
      End Do
      Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(boundaryTemperatureSec,ierr);CHKERRQ(ierr)      
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call VecAssemblyBegin(x,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(x,ierr);CHKERRQ(ierr)
End Subroutine MEF90HeatXferUpdateboundaryTemperature

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
      
      PetscInt                                           :: dim      
      Call DMMeshGetDimension(MEF90HeatXferCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90HeatXferOperator2D(snesTemp,x,residual,MEF90HeatXferCtx,ierr)
      Else If (dim == 3) Then
         Call MEF90HeatXferOperator3D(snesTemp,x,residual,MEF90HeatXferCtx,ierr)
      End If      
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

      PetscInt                                           :: dim      
      Call DMMeshGetDimension(MEF90HeatXferCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90HeatXferBilinearForm2D(snesTemp,x,A,M,flg,MEF90HeatXferCtx,ierr)
      Else If (dim == 3) Then
         Call MEF90HeatXferBilinearForm3D(snesTemp,x,A,M,flg,MEF90HeatXferCtx,ierr)
      End If      
   End Subroutine MEF90HeatXferBilinearForm

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXFerEnergy"
!!!
!!!  
!!!  MEF90HeatXFerEnergy:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXFerEnergy(temperatureVec,t,MEF90HeatXferCtx,energy,work,ierr)
      Type(Vec),Intent(IN)                            :: temperatureVec
      PetscReal,Intent(IN)                            :: t
      Type(MEF90HeatXferCtx_Type),Intent(IN)          :: MEF90HeatXferCtx
      PetscReal,Dimension(:),Pointer                  :: energy,work
      PetscErrorCode,Intent(OUT)                      :: ierr

      PetscInt                                           :: dim      
      Call DMMeshGetDimension(MEF90HeatXferCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90HeatXFerEnergy2D(temperatureVec,t,MEF90HeatXferCtx,energy,work,ierr)
      Else If (dim == 3) Then
         Call MEF90HeatXFerEnergy3D(temperatureVec,t,MEF90HeatXferCtx,energy,work,ierr)
      End If      
   End Subroutine MEF90HeatXFerEnergy
End Module m_MEF90_HeatXfer