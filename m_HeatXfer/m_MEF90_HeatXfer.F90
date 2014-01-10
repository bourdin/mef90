#include "../MEF90/mef90.inc"
Module m_MEF90_HeatXfer
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_HeatXferCtx
   Use m_MEF90_HeatXferAssembly2D, &
      MEF90HeatXferEnergy2D      => MEF90HeatXferEnergy, &
      MEF90HeatXferOperator2D    => MEF90HeatXferOperator, &
      MEF90HeatXferBilinearForm2D => MEF90HeatXferBilinearForm, &
      MEF90HeatXferIFunction2D    => MEF90HeatXferIFunction, &
      MEF90HeatXferIJacobian2D    => MEF90HeatXferIJacobian
   Use m_MEF90_HeatXferAssembly3D, &
      MEF90HeatXferEnergy3D      => MEF90HeatXFerEnergy, &
      MEF90HeatXferOperator3D    => MEF90HeatXferOperator, &
      MEF90HeatXferBilinearForm3D => MEF90HeatXferBilinearForm, &
      MEF90HeatXferIFunction3D    => MEF90HeatXferIFunction, &
      MEF90HeatXferIJacobian3D    => MEF90HeatXferIJacobian

   Implicit none

   Private
   Public MEF90HeatXferOperator
   Public MEF90HeatXferBilinearForm
   Public MEF90HeatXferEnergy
   Public MEF90HeatXferSetTransients
   Public MEF90HeatXferUpdateboundaryTemperature
   Public MEF90HeatXferIFunction
   Public MEF90HeatXferIJacobian
   Public MEF90HeatXferViewEXO
   Public MEF90HeatXferFormatEXO
   Public MEF90HeatXferCreateSNES
   Public MEF90HeatXferCreateTS
Contains

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferSetTransients"
!!!
!!!  
!!!  MEF90HeatXferSetTransients: Update all transient data (boundary / external temperature and fluxes)
!!!                              using the proper scaling law
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
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
!!!  MEF90HeatXferSetFluxCst: low level function called by MEF90HeatXferSetTransients
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
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
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(cellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(cellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call VecAssemblyBegin(x,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(x,ierr);CHKERRQ(ierr)
   End Subroutine MEF90HeatXferSetFluxCst

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferSetexternalTemperatureCst"
!!!
!!!  MEF90HeatXferSetexternalTemperatureCst: low level function called by MEF90HeatXferSetTransients
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
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
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(cellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(cellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call VecAssemblyBegin(x,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(x,ierr);CHKERRQ(ierr)
   End Subroutine MEF90HeatXferSetexternalTemperatureCst

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferSetboundaryTemperatureCst"
!!!
!!!  MEF90HeatXferSetboundaryTemperatureCst: low level function called by MEF90HeatXferSetTransients
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
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
      PetscInt                                           :: set,nval
      PetscReal,Dimension(:),Pointer                     :: val
      PetscInt,Dimension(:),Pointer                      :: cone
      Type(MEF90Element_Type)                            :: elemType
      PetscInt                                           :: cell,dof
      
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90HeatXferCtx%MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr);CHKERRQ(ierr)
      
      !!! boundaryTemperature is Vertex-centered
      !!! We first set the boundary values inherited from cell sets, then that of vertex sets
      !!!
      !!! Cell Sets
      Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)         
         If (cellSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(MEF90HeatXferCtx%dm,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
            Call MEF90_ISCreateCelltoVertex(MEF90HeatXferCtx%dm,PETSC_COMM_WORLD,setIS,bcIS,ierr)
            Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
            Allocate(val(nval),stat=ierr)
            val = cellSetOptions%boundaryTemp
            Call MEF90_VecSetValuesISdof(MEF90HeatXferCtx%DM,x,val,bcIS,1,INSERT_VALUES,ierr)
            DeAllocate(val)
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
         If (vertexSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(MEF90HeatXferCtx%dm,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
            Call ISGetSize(setIS,nval,ierr);CHKERRQ(ierr)
            Allocate(val(nval),stat=ierr)
            val = vertexSetOptions%boundaryTemp
            Call MEF90_VecSetValuesISdof(MEF90HeatXferCtx%DM,x,val,setIS,1,INSERT_VALUES,ierr)
            DeAllocate(val)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         EndIf
      End Do
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call VecAssemblyBegin(x,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(x,ierr);CHKERRQ(ierr)
   End Subroutine MEF90HeatXferSetboundaryTemperatureCst

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferUpdateboundaryTemperature"
!!!
!!!  
!!!  MEF90HeatXferUpdateboundaryTemperature: Update the solution vector x with the boundary temperature from 
!!!                                          MEF90HeatXferCtx%BoundaryTemperature
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
      Type(IS)                                           :: VertexSetGlobalIS,setIS
      Type(MEF90HeatXferCellSetOptions_Type),pointer     :: cellSetOptions
      Type(IS)                                           :: CellSetGlobalIS,bCIS
      PetscInt,Dimension(:),Pointer                      :: setID
      PetscInt                                           :: set,nval
      PetscReal,Dimension(:),Pointer                     :: boundaryTemperaturePtr,xPtr
      Type(MEF90Element_Type)                            :: elemType
      
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90HeatXferCtx%MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr);CHKERRQ(ierr)
      
      !!! Cell Sets
      Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)         
         If (cellSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(MEF90HeatXferCtx%dm,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
            Call MEF90_ISCreateCelltoVertex(MEF90HeatXferCtx%dm,PETSC_COMM_WORLD,setIS,bcIS,ierr)
            Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
            Allocate(xPtr(nval),stat=ierr)
            Call MEF90_VecGetValuesISdof(MEF90HeatXferCtx%DM,MEF90HeatXferCtx%boundaryTemperature,xPtr,bcIS,1,ierr)
            Call MEF90_VecSetValuesISdof(MEF90HeatXferCtx%DM,x,xPtr,bcIS,1,INSERT_VALUES,ierr)
            DeAllocate(xPtr)
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
         If (vertexSetOptions%Has_BC) Then
            Call DMMeshGetStratumIS(MEF90HeatXferCtx%dm,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
            Call ISGetSize(setIS,nval,ierr);CHKERRQ(ierr)
            Allocate(xPtr(nval),stat=ierr)
            Call MEF90_VecGetValuesISdof(MEF90HeatXferCtx%DM,MEF90HeatXferCtx%boundaryTemperature,xPtr,setIS,1,ierr)
            Call MEF90_VecSetValuesISdof(MEF90HeatXferCtx%DM,x,xPtr,setIS,1,INSERT_VALUES,ierr)
            DeAllocate(xPtr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         EndIf
      End Do
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)

      Call VecAssemblyBegin(x,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(x,ierr);CHKERRQ(ierr)
End Subroutine MEF90HeatXferUpdateboundaryTemperature

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferOperator"
!!!
!!!  
!!!  MEF90HeatXferOperator: wraps calls to MEF90HeatXferOperator from m_MEF90_HeatXferAssembly
!!!                         since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
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
!!!  MEF90HeatXferBilinearForm: wraps calls to MEF90HeatXferBilinearForm from m_MEF90_HeatXferAssembly
!!!                             since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
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
!!!  MEF90HeatXFerEnergy: wraps calls to MEF90HeatXferEnergy from m_MEF90_HeatXferAssembly
!!!                       since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXFerEnergy(temperatureVec,t,MEF90HeatXferCtx,energy,work,ierr)
      Type(Vec),Intent(IN)                               :: temperatureVec
      PetscReal,Intent(IN)                               :: t
      Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
      PetscReal,Dimension(:),Pointer                     :: energy,work
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      PetscInt                                           :: dim      

      Call DMMeshGetDimension(MEF90HeatXferCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90HeatXFerEnergy2D(temperatureVec,t,MEF90HeatXferCtx,energy,work,ierr)
      Else If (dim == 3) Then
         Call MEF90HeatXFerEnergy3D(temperatureVec,t,MEF90HeatXferCtx,energy,work,ierr)
      End If      
   End Subroutine MEF90HeatXFerEnergy

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXFerIFunction"
!!!
!!!  
!!!  MEF90HeatXFerIFunction: wraps calls to MEF90HeatXFerIFunction from m_MEF90_HeatXferAssembly
!!!                         since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXFerIFunction(tempTS,time,x,xdot,F,MEF90HeatXferCtx,ierr)
      Type(TS),Intent(IN)                                :: tempTS
      PetscReal,Intent(IN)                               :: time
      Type(Vec),Intent(IN)                               :: x,xdot
      Type(Vec),Intent(INOUT)                            :: F
      Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      PetscInt                                           :: dim      

      Call DMMeshGetDimension(MEF90HeatXferCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90HeatXFerIFunction2D(tempTS,time,x,xdot,F,MEF90HeatXferCtx,ierr)
      Else If (dim == 3) Then
         Call MEF90HeatXFerIFunction3D(tempTS,time,x,xdot,F,MEF90HeatXferCtx,ierr)
      End If      
   End Subroutine MEF90HeatXFerIFunction
   
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferIJacobian"
!!!
!!!  
!!!  MEF90HeatXferIJacobian: wraps calls to MEF90HeatXferIJacobian from m_MEF90_HeatXferAssembly
!!!                         since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferIJacobian(tempTS,t,x,xdot,shift,A,M,flg,MEF90HeatXferCtx,ierr)
      Type(TS),Intent(IN)                                :: tempTS
      PetscReal,Intent(IN)                               :: t
      Type(Vec),Intent(IN)                               :: x,xdot
      PetscReal,Intent(IN)                               :: shift
      Type(Mat),Intent(INOUT)                            :: A,M
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
      PetscErrorCode,Intent(OUT)                         :: ierr  
      
      PetscInt                                           :: dim      

      Call DMMeshGetDimension(MEF90HeatXferCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90HeatXferIJacobian2D(tempTS,t,x,xdot,shift,A,M,flg,MEF90HeatXferCtx,ierr)
      Else If (dim == 3) Then
         Call MEF90HeatXferIJacobian3D(tempTS,t,x,xdot,shift,A,M,flg,MEF90HeatXferCtx,ierr)
      End If      
   End Subroutine MEF90HeatXferIJacobian
   
   
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferViewEXO"
!!!
!!!  
!!!  MEF90HeatXferViewEXO:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferViewEXO(MEF90HeatXferCtx,step,ierr)
      Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
      PetscInt,Intent(IN)                                :: step
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(Vec)                                          :: localVec
      Type(MEF90HeatXferGlobalOptions_Type),pointer      :: MEF90HeatXferGlobalOptions

      Call PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr);CHKERRQ(ierr)
      Call DMGetLocalVector(MEF90HeatXferCtx%cellDM,localVec,ierr);CHKERRQ(ierr)
      If (MEF90HeatXferGlobalOptions%fluxOffset > 0) Then
         Call DMGlobalToLocalBegin(MEF90HeatXferCtx%cellDM,MEF90HeatXferCtx%flux,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalEnd(MEF90HeatXferCtx%cellDM,MEF90HeatXferCtx%flux,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call VecViewExodusCell(MEF90HeatXferCtx%cellDM,localVec,MEF90HeatXferCtx%MEF90Ctx%IOcomm, &
                                MEF90HeatXferCtx%MEF90Ctx%fileExoUnit,step,MEF90HeatXferGlobalOptions%fluxOffset,ierr);CHKERRQ(ierr)
      End If
      
      If (MEF90HeatXferGlobalOptions%externalTempOffset > 0) Then
         Call DMGlobalToLocalBegin(MEF90HeatXferCtx%cellDM,MEF90HeatXferCtx%externalTemperature,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalEnd(MEF90HeatXferCtx%cellDM,MEF90HeatXferCtx%externalTemperature,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call VecViewExodusCell(MEF90HeatXferCtx%cellDM,localVec,MEF90HeatXferCtx%MEF90Ctx%IOcomm, &
                                MEF90HeatXferCtx%MEF90Ctx%fileExoUnit,step,MEF90HeatXferGlobalOptions%externalTempOffset,ierr);CHKERRQ(ierr)
      End If
      Call DMRestoreLocalVector(MEF90HeatXferCtx%cellDM,localVec,ierr);CHKERRQ(ierr)
      
      Call DMGetLocalVector(MEF90HeatXferCtx%DM,localVec,ierr);CHKERRQ(ierr)
      If (MEF90HeatXferGlobalOptions%tempOffset > 0) Then
         Call DMGlobalToLocalBegin(MEF90HeatXferCtx%DM,MEF90HeatXferCtx%temperature,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalEnd(MEF90HeatXferCtx%DM,MEF90HeatXferCtx%temperature,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call VecViewExodusVertex(MEF90HeatXferCtx%DM,localVec,MEF90HeatXferCtx%MEF90Ctx%IOcomm, &
                                  MEF90HeatXferCtx%MEF90Ctx%fileExoUnit,step,MEF90HeatXferGlobalOptions%tempOffset,ierr);CHKERRQ(ierr)
      End If

      If (MEF90HeatXferGlobalOptions%boundaryTempOffset > 0) Then
         Call DMGlobalToLocalBegin(MEF90HeatXferCtx%DM,MEF90HeatXferCtx%boundaryTemperature,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalEnd(MEF90HeatXferCtx%DM,MEF90HeatXferCtx%boundaryTemperature,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call VecViewExodusVertex(MEF90HeatXferCtx%DM,localVec,MEF90HeatXferCtx%MEF90Ctx%IOcomm, &
                                  MEF90HeatXferCtx%MEF90Ctx%fileExoUnit,step,MEF90HeatXferGlobalOptions%boundaryTempOffset,ierr);CHKERRQ(ierr)
         End If
      Call DMRestoreLocalVector(MEF90HeatXferCtx%DM,localVec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90HeatXferViewEXO

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferFormatEXO"
!!!
!!!  
!!!  MEF90HeatXferFormatEXO:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferFormatEXO(MEF90HeatXferCtx,ierr)
      Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
      PetscErrorCode,Intent(OUT)                         :: ierr

      Character(len=MXSTLN),Dimension(:),Pointer         :: nameG,nameV,nameC
      Type(MEF90HeatXferGlobalOptions_Type),pointer      :: MEF90HeatXferGlobalOptions
      Integer                                            :: numfield

      Call PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr);CHKERRQ(ierr)
      Allocate(nameG(0))
      !nameG(1) = "Energy"
      !nameG(2) = "work"
   
      numfield = max(MEF90HeatXferGlobalOptions%tempOffset, &
                     MEF90HeatXferGlobalOptions%boundaryTempOffset)
      Allocate(nameV(numfield))
      nameV = "empty"
      If (MEF90HeatXferGlobalOptions%tempOffset > 0) Then
         nameV(MEF90HeatXferGlobalOptions%tempOffset) = "temperature"
      End If
      If (MEF90HeatXferGlobalOptions%boundaryTempOffset > 0) Then
         nameV(MEF90HeatXferGlobalOptions%boundaryTempOffset) = "boundary temperature"
      End If
                     
      numfield = max(MEF90HeatXferGlobalOptions%externalTempOffset, &
                     MEF90HeatXferGlobalOptions%fluxOffset)
      Allocate(nameC(numfield))
      nameC = "empty"
      If (MEF90HeatXferGlobalOptions%externalTempOffset > 0) Then
         nameC(MEF90HeatXferGlobalOptions%externalTempOffset) = "external temperature"
      End If
      If (MEF90HeatXferGlobalOptions%fluxOffset > 0) Then
         nameC(MEF90HeatXferGlobalOptions%fluxOffset) = "heat flux"
      End If
      
      Call MEF90EXOFormat(MEF90HeatXferCtx%MEF90Ctx%fileEXOUNIT,nameG,nameC,nameV,ierr)
      DeAllocate(nameG)
      DeAllocate(nameV)
      DeAllocate(nameC)
   End Subroutine MEF90HeatXferFormatEXO

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferCreateSNES"
!!!
!!!  
!!!  MEF90HeatXferCreateSNES:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferCreateSNES(MEF90HeatXferCtx,snesTemp,residual,ierr)
      Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
      Type(SNES),Intent(OUT)                             :: snesTemp
      Type(Vec),Intent(IN)                               :: residual
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(MEF90HeatXferGlobalOptions_Type),pointer      :: MEF90HeatXferGlobalOptions
      Type(Mat)                                          :: matTemp
      Type(MatNullSpace)                                 :: nspTemp
      Type(KSP)                                          :: kspTemp
      Type(PC)                                           :: pcTemp
      PetscReal                                          :: rtol,dtol

      Call PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr);CHKERRQ(ierr)
      Call DMCreateMatrix(MEF90HeatXferCtx%DM,MATAIJ,matTemp,iErr);CHKERRQ(iErr)
      Call MatSetOptionsPrefix(matTemp,"temp_",ierr);CHKERRQ(ierr)
      Call MatSetOption(matTemp,MAT_SPD,PETSC_TRUE,ierr);CHKERRQ(ierr)
      Call MatSetOption(matTemp,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr);CHKERRQ(ierr)
      Call MatSetOption(matTemp,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr);CHKERRQ(ierr)
      If (MEF90HeatXferGlobalOptions%addNullSpace) Then
         Call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,PETSC_NULL_OBJECT,nspTemp,ierr);CHKERRQ(ierr)
         Call MatSetNullSpace(matTemp,nspTemp,ierr);CHKERRQ(ierr)
      End If
      Call MatSetFromOptions(matTemp,ierr);CHKERRQ(ierr)

      Call SNESCreate(PETSC_COMM_WORLD,snesTemp,ierr);CHKERRQ(ierr)
      Call SNESSetApplicationContext(snesTemp,MEF90HeatXferCtx,ierr);CHKERRQ(ierr)
      Call SNESSetDM(snesTemp,MEF90HeatXferCtx%DM,ierr);CHKERRQ(ierr)
      Call SNESSetType(snesTemp,SNESKSPONLY,ierr);CHKERRQ(ierr)
      Call SNESSetOptionsPrefix(snesTemp,'temp_',ierr);CHKERRQ(ierr)

      Call SNESSetFunction(snesTemp,residual,MEF90HeatXferOperator,MEF90HeatXferCtx,ierr);CHKERRQ(ierr)
      Call SNESSetJacobian(snesTemp,matTemp,matTemp,MEF90HeatXferBilinearForm,MEF90HeatXferCtx,ierr);CHKERRQ(ierr)
      Call SNESSetFromOptions(snesTemp,ierr);CHKERRQ(ierr)
      !!! 
      !!! Set some KSP options
      !!!
      Call SNESGetKSP(snesTemp,kspTemp,ierr);CHKERRQ(ierr)
      Call KSPSetType(kspTemp,KSPCG,ierr);CHKERRQ(ierr)
      Call KSPSetInitialGuessNonzero(kspTemp,PETSC_TRUE,ierr);CHKERRQ(ierr)
      If (MEF90HeatXferGlobalOptions%addNullSpace) Then
         Call KSPSetNullSpace(kspTemp,nspTemp,ierr);CHKERRQ(ierr)
      End If
      rtol = 1.0D-8
      dtol = 1.0D+10
      Call KSPSetTolerances(kspTemp,rtol,PETSC_DEFAULT_DOUBLE_PRECISION,dtol,PETSC_DEFAULT_INTEGER,ierr);CHKERRQ(ierr)
      Call KSPSetFromOptions(kspTemp,ierr);CHKERRQ(ierr)
   End Subroutine MEF90HeatXferCreateSNES

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferCreateTS"
!!!
!!!  
!!!  MEF90HeatXferCreateTS:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferCreateTS(MEF90HeatXferCtx,tsTemp,residual,ierr)
      Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
      Type(TS),Intent(OUT)                               :: tsTemp
      Type(Vec),Intent(IN)                               :: residual
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(MEF90HeatXferGlobalOptions_Type),pointer      :: MEF90HeatXferGlobalOptions
      Type(Mat)                                          :: matTemp
      Type(MatNullSpace)                                 :: nspTemp
      Type(SNES)                                         :: snesTemp
      Type(KSP)                                          :: kspTemp
      Type(PC)                                           :: pcTemp
      PetscReal                                          :: rtol,dtol

      Call PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr);CHKERRQ(ierr)
      Call DMCreateMatrix(MEF90HeatXferCtx%DM,MATAIJ,matTemp,iErr);CHKERRQ(iErr)
      Call MatSetOptionsPrefix(matTemp,"temp_",ierr);CHKERRQ(ierr)
      Call MatSetOption(matTemp,MAT_SPD,PETSC_TRUE,ierr);CHKERRQ(ierr)
      Call MatSetOption(matTemp,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr);CHKERRQ(ierr)
      Call MatSetOption(matTemp,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr);CHKERRQ(ierr)
      If (MEF90HeatXferGlobalOptions%addNullSpace) Then
         Call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,PETSC_NULL_OBJECT,nspTemp,ierr);CHKERRQ(ierr)
         Call MatSetNullSpace(matTemp,nspTemp,ierr);CHKERRQ(ierr)
      End If
      Call MatSetFromOptions(matTemp,ierr);CHKERRQ(ierr)

      Call TSCreate(PETSC_COMM_WORLD,tsTemp,ierr);CHKERRQ(ierr)
      Call TSSetDM(tsTemp,MEF90HeatXferCtx%DM,ierr);CHKERRQ(ierr)
      Call TSSetOptionsPrefix(tsTemp,'temp_',ierr);CHKERRQ(ierr)
      Call TSGetSNES(tsTemp,snesTemp,ierr);CHKERRQ(ierr)

      Call TSSetIFunction(tsTemp,residual,MEF90HeatXFerIFunction,MEF90HeatXferCtx,ierr);CHKERRQ(ierr)
      Call TSSetIJacobian(tsTemp,matTemp,matTemp,MEF90HeatXFerIJacobian,MEF90HeatXferCtx,ierr);CHKERRQ(ierr)

      Call TSSetType(tsTemp,'rosw',ierr);CHKERRQ(ierr)
      Call TSRosWSetType(tsTemp,'ra3pw',ierr);CHKERRQ(ierr)
      Call TSSetProblemType(tsTemp,TS_LINEAR,ierr);CHKERRQ(ierr)
      Call VecSet(MEF90HeatXferCtx%temperature,MEF90HeatXferGlobalOptions%initialTemperature,ierr);CHKERRQ(ierr)
      Call TSSetSolution(tsTemp,MEF90HeatXferCtx%temperature,ierr);CHKERRQ(ierr)
      Call TSSetExactFinalTime(tsTemp,PETSC_TRUE,ierr);CHKERRQ(ierr)
      Call TSSetFromOptions(tsTemp,ierr);CHKERRQ(ierr)
      !!! 
      !!! Set some KSP options
      !!!
      Call SNESGetKSP(snesTemp,kspTemp,ierr);CHKERRQ(ierr)
      Call KSPSetType(kspTemp,KSPCG,ierr);CHKERRQ(ierr)
      Call KSPSetInitialGuessNonzero(kspTemp,PETSC_TRUE,ierr);CHKERRQ(ierr)
      If (MEF90HeatXferGlobalOptions%addNullSpace) Then
         Call KSPSetNullSpace(kspTemp,nspTemp,ierr);CHKERRQ(ierr)
      End If
      rtol = 1.0D-8
      dtol = 1.0D+10
      Call KSPSetTolerances(kspTemp,rtol,PETSC_DEFAULT_DOUBLE_PRECISION,dtol,PETSC_DEFAULT_INTEGER,ierr);CHKERRQ(ierr)
      Call KSPSetFromOptions(kspTemp,ierr);CHKERRQ(ierr)
   End Subroutine MEF90HeatXferCreateTS
End Module m_MEF90_HeatXfer