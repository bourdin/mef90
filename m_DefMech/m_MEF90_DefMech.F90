#include "../MEF90/mef90.inc"
Module m_MEF90_DefMech
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_DefMechCtx

   Use m_MEF90_DefMechAssembly2D, &
      MEF90DefMechOperator2D        => MEF90DefMechOperator,     &
      MEF90DefMechBilinearForm2D    => MEF90DefMechBilinearForm, &     
      MEF90DefMechWork2D            => MEF90DefMechWork,         &
      MEF90DefMechElasticEnergy2D   => MEF90DefMechElasticEnergy
   Use m_MEF90_DefMechAssembly3D, &
      MEF90DefMechOperator3D        => MEF90DefMechOperator,     &
      MEF90DefMechBilinearForm3D    => MEF90DefMechBilinearForm, &     
      MEF90DefMechWork3D            => MEF90DefMechWork,         &
      MEF90DefMechElasticEnergy3D   => MEF90DefMechElasticEnergy

   Implicit none
   
Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSetTransients"
!!!
!!!  
!!!  MEF90DefMechSetTransients:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechSetTransients(MEF90DefMechCtx,step,time,ierr)
      Type(MEF90DefMechCtx_Type),Intent(INOUT)        :: MEF90DefMechCtx
      PetscInt,Intent(IN)                             :: step
      PetscReal,Intent(IN)                            :: time
      PetscErrorCode,Intent(OUT)                      :: ierr
   
      Type(MEF90DefMechGlobalOptions_Type),pointer    :: MEF90DefMechGlobalOptions
      Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90GlobalOptions

      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)

      Select case (MEF90DefMechGlobalOptions%forceScaling)
         Case (MEF90Scaling_File)
            Call VecLoadExodusCell(MEF90DefMechCtx%cellDMVect,MEF90DefMechCtx%force,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
                                   MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%forceOffset,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_Linear)
            Call MEF90DefMechSetForceCst(MEF90DefMechCtx%force,MEF90DefMechCtx,ierr)
            Call VecScale(MEF90DefMechCtx%force,time,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_CST)
            Call MEF90DefMechSetforceCst(MEF90DefMechCtx%force,MEF90DefMechCtx,ierr)
      End Select
      Select case (MEF90DefMechGlobalOptions%pressureForceScaling)
         Case (MEF90Scaling_File)
            Call VecLoadExodusCell(MEF90DefMechCtx%cellDMScal,MEF90DefMechCtx%pressureForce,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
                                   MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%pressureForceOffset,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_Linear)
            Call MEF90DefMechSetPressureForceCst(MEF90DefMechCtx%pressureForce,MEF90DefMechCtx,ierr)
            Call VecScale(MEF90DefMechCtx%pressureForce,time,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_CST)
            Call MEF90DefMechSetforceCst(MEF90DefMechCtx%pressureForce,MEF90DefMechCtx,ierr)
      End Select
      Select case (MEF90DefMechGlobalOptions%boundaryDisplacementScaling)
         Case (MEF90Scaling_File)
            Call VecLoadExodusVertex(MEF90DefMechCtx%cellDMVect,MEF90DefMechCtx%boundaryDisplacement,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
                                     MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%boundaryDisplacementOffset,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_Linear)
            Call MEF90DefMechSetboundaryDisplacementCst(MEF90DefMechCtx%boundaryDisplacement,MEF90DefMechCtx,ierr)
            Call VecScale(MEF90DefMechCtx%boundaryDisplacement,time,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_CST)
            Call MEF90DefMechSetboundaryDisplacementCst(MEF90DefMechCtx%boundaryDisplacement,MEF90DefMechCtx,ierr)
      End Select
   End Subroutine MEF90DefMechSetTransients

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSetForceCst"
!!!
!!!  
!!!  MEF90DefMechSetForceCst:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechSetForceCst(x,MEF90DefMechCtx,ierr)
      Type(Vec),Intent(IN)                            :: x
      Type(MEF90DefMechCtx_Type),Intent(IN)           :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                      :: ierr
   
      Type(SectionReal)                               :: xSec
      Type(VecScatter)                                :: scatterSecToVec
      Type(MEF90DefMechGlobalOptions_Type),pointer    :: MEF90DefMechGlobalOptions
      Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90GlobalOptions
      Type(MEF90DefMechCellSetOptions_Type),pointer   :: cellSetOptions
      Type(IS)                                        :: cellSetGlobalIS,setIS
      PetscInt,Dimension(:),Pointer                   :: setID
      PetscInt,Dimension(:),Pointer                   :: setIdx
      PetscInt                                        :: set,c,dim
      PetscReal,Dimension(:),Pointer                  :: val
      
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)
      Call DMMeshGetDimension(MEF90DefMechCtx%cellDMVect,dim,ierr);CHKERRQ(ierr)

      Call DMMeshGetSectionReal(MEF90DefMechCtx%cellDMVect,'default',xSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%cellDMVect,xSec,ScatterSecToVec,ierr);CHKERRQ(ierr)

      !!! force is cell-centered
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%CellDMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%cellDMVect,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
        Call ISGetIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
         Do c = 1, size(setIdx)
            Call SectionRealRestrict(xSec,setIDx(c),val,ierr);CHKERRQ(ierr)
            val = cellSetOptions%force
            Call SectionRealRestore(xSec,setIDx(c),val,ierr);CHKERRQ(ierr)
         End Do
         Call ISRestoreIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(cellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(cellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(xSec,ScatterSecToVec,SCATTER_FORWARD,x,ierr);CHKERRQ(ierr) 
      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechSetForceCst

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSetPressureForceCst"
!!!
!!!  
!!!  MEF90DefMechSetPressureForceCst:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechSetPressureForceCst(x,MEF90DefMechCtx,ierr)
      Type(Vec),Intent(IN)                            :: x
      Type(MEF90DefMechCtx_Type),Intent(IN)           :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                      :: ierr
   
      Type(SectionReal)                               :: xSec
      Type(VecScatter)                                :: scatterSecToVec
      Type(MEF90DefMechGlobalOptions_Type),pointer    :: MEF90DefMechGlobalOptions
      Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90GlobalOptions
      Type(MEF90DefMechCellSetOptions_Type),pointer   :: cellSetOptions
      Type(IS)                                        :: cellSetGlobalIS,setIS
      PetscInt,Dimension(:),Pointer                   :: setID
      PetscInt,Dimension(:),Pointer                   :: setIdx
      PetscInt                                        :: set,c,dim
      PetscReal,Dimension(:),Pointer                  :: val
      
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)
      Call DMMeshGetDimension(MEF90DefMechCtx%cellDMVect,dim,ierr);CHKERRQ(ierr)

      Call DMMeshGetSectionReal(MEF90DefMechCtx%cellDMScal,'default',xSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%cellDMScal,xSec,ScatterSecToVec,ierr);CHKERRQ(ierr)

      !!! pressure force is cell-centered
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%CellDMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%cellDMVect,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
        Call ISGetIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
         Do c = 1, size(setIdx)
            Call SectionRealRestrict(xSec,setIDx(c),val,ierr);CHKERRQ(ierr)
            val = cellSetOptions%pressureForce
            Call SectionRealRestore(xSec,setIDx(c),val,ierr);CHKERRQ(ierr)
         End Do
         Call ISRestoreIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(cellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(cellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(xSec,ScatterSecToVec,SCATTER_FORWARD,x,ierr);CHKERRQ(ierr) 
      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechSetPressureForceCst

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSetBoundaryDisplacementCst"
!!!
!!!  
!!!  MEF90DefMechSetBoundaryDisplacementCst:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechSetBoundaryDisplacementCst(x,MEF90DefMechCtx,ierr)
      Type(Vec),Intent(IN)                               :: x
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
      Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions
      Type(MEF90DefMechVertexSetOptions_Type),pointer    :: vertexSetOptions
      Type(IS)                                           :: VertexSetGlobalIS,setIS,setISdof,bcIS
      Type(MEF90DefMechCellSetOptions_Type),pointer      :: cellSetOptions
      Type(IS)                                           :: CellSetGlobalIS
      PetscInt,Dimension(:),Pointer                      :: setID
      PetscInt,Dimension(:),Pointer                      :: setIdx
      PetscInt                                           :: set,nval
      PetscReal,Dimension(:),Pointer                     :: val
      PetscInt,Dimension(:),Pointer                      :: cone
      Type(MEF90Element_Type)                            :: elemType
      PetscInt                                           :: cell,dof,c,dim
      
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)
      Call DMMeshGetDimension(MEF90DefMechCtx%DMVect,dim,ierr);CHKERRQ(ierr)
      
      !!! boundaryDisplacement is Vertex-centered
      !!! We first set the boundary values inherited from cell sets, then that of vertex sets
      !!!
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)   
         Do c = 1, dim
            If (cellSetOptions%Has_displacementBC(c)) Then
               Call DMMeshGetStratumIS(MEF90DefMechCtx%DMVect,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
               Call MEF90_ISCreateCelltoVertex(MEF90DefMechCtx%DMVect,PETSC_COMM_WORLD,setIS,bcIS,ierr)
               Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
               Allocate(val(nval),stat=ierr)
               val =  cellSetOptions%boundaryDisplacement(c)
               Call MEF90_VecSetValuesISdof(MEF90DefMechCtx%DMVect,x,val,bcIS,c,INSERT_VALUES,ierr)
               DeAllocate(Val)
               Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
               Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
            End If! cellSetOptions%Has_BC
         End Do ! c
      End Do ! set
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
      
      !!! Vertex Sets
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxVertexSetOptions(MEF90DefMechCtx%vertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)   
         Do c = 1, dim
            If (vertexSetOptions%Has_displacementBC(c)) Then
               Call DMMeshGetStratumIS(MEF90DefMechCtx%DMVect,'Vertex Sets',setID(set),bcIS,ierr);CHKERRQ(iErr)
               Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
               Allocate(val(nval),stat=ierr)
               val =  vertexSetOptions%boundaryDisplacement(c)
               Call MEF90_VecSetValuesISdof(MEF90DefMechCtx%DMVect,x,val,bcIS,c,INSERT_VALUES,ierr)
               DeAllocate(Val)
               Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
            End If! cellSetOptions%Has_BC
         End Do ! c
      End Do ! set
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call VecAssemblyBegin(x,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(x,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechSetBoundaryDisplacementCst

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechUpdateboundaryDisplacement"
!!!
!!!  
!!!  MEF90DefMechUpdateboundaryDisplacement: Update the solution vector x with the boundary Displacement from 
!!!                                          MEF90DefMechCtx%BoundaryDisplacement
!!!  
!!!  (c) 2013 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechUpdateboundaryDisplacement(x,MEF90DefMechCtx,ierr)
      Type(Vec),Intent(IN)                               :: x
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr

   
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
      Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions
      Type(MEF90DefMechVertexSetOptions_Type),pointer    :: vertexSetOptions
      Type(IS)                                           :: VertexSetGlobalIS,setIS,setISdof
      Type(MEF90DefMechCellSetOptions_Type),pointer      :: cellSetOptions
      Type(IS)                                           :: CellSetGlobalIS,bcIS
      PetscInt,Dimension(:),Pointer                      :: setID
      PetscInt                                           :: set,c,dim,nval
      PetscReal,Dimension(:),Pointer                     :: boundaryDisplacementPtr
      Type(MEF90Element_Type)                            :: elemType
      
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)
      Call DMMeshGetDimension(MEF90DefMechCtx%DMVect,dim,ierr);CHKERRQ(ierr)
      
      !!!
      !!! cell set Displacement first, followed vertex sets
      !!!
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)   
         Do c = 1, dim
            If (cellSetOptions%Has_displacementBC(c)) Then
               Call DMMeshGetStratumIS(MEF90DefMechCtx%DMVect,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
               Call MEF90_ISCreateCelltoVertex(MEF90DefMechCtx%DMVect,PETSC_COMM_WORLD,setIS,bcIS,ierr)
               Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
               Allocate(boundaryDisplacementPtr(nval),stat=ierr)
               Call MEF90_VecGetValuesISdof(MEF90DefMechCtx%DMVect,MEF90DefMechCtx%boundaryDisplacement,boundaryDisplacementPtr,bcIS,c,ierr)
               Call MEF90_VecSetValuesISdof(MEF90DefMechCtx%DMVect,x,boundaryDisplacementPtr,bcIS,c,INSERT_VALUES,ierr)
               DeAllocate(boundaryDisplacementPtr)
               Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
               Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
            End If! cellSetOptions%Has_BC
         End Do ! c
      End Do ! set
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)

      !!! Vertex sets
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxVertexSetOptions(MEF90DefMechCtx%vertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)   
         Do c = 1, dim
            If (vertexSetOptions%Has_displacementBC(c)) Then
               Call DMMeshGetStratumIS(MEF90DefMechCtx%DMVect,'Vertex Sets',setID(set),bcIS,ierr);CHKERRQ(iErr)
               Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
               Allocate(boundaryDisplacementPtr(nval),stat=ierr)
               Call MEF90_VecGetValuesISdof(MEF90DefMechCtx%DMVect,MEF90DefMechCtx%boundaryDisplacement,boundaryDisplacementPtr,bcIS,c,ierr)
               Call MEF90_VecSetValuesISdof(MEF90DefMechCtx%DMVect,x,boundaryDisplacementPtr,bcIS,c,INSERT_VALUES,ierr)
               DeAllocate(boundaryDisplacementPtr)
               Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
            End If! cellSetOptions%Has_BC
         End Do ! c
      End Do ! set
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)

      Call VecAssemblyBegin(x,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(x,ierr);CHKERRQ(ierr)
End Subroutine MEF90DefMechUpdateboundaryDisplacement

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperator"
!!!
!!!  
!!!  MEF90DefMechOperator: wraps calls to MEF90DefMechOperator from m_MEF90_DefMechAssembly
!!!                        since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechOperator(snesTemp,x,residual,MEF90DefMechCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesTemp
      Type(Vec),Intent(IN)                               :: x
      Type(Vec),Intent(INOUT)                            :: residual
      Type(MEF90DefMechCtx_Type),Intent(IN)             :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      PetscInt                                           :: dim      
      Call DMMeshGetDimension(MEF90DefMechCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90DefMechOperator2D(snesTemp,x,residual,MEF90DefMechCtx,ierr)
      Else If (dim == 3) Then
         Call MEF90DefMechOperator3D(snesTemp,x,residual,MEF90DefMechCtx,ierr)
      End If      
   End Subroutine MEF90DefMechOperator
   
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearForm"
!!!
!!!  
!!!  MEF90DefMechBilinearForm: wraps calls to MEF90DefMechBilinearForm from m_MEF90_DefMechAssembly
!!!                            since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechBilinearForm(snesDispl,x,A,M,flg,MEF90DefMechCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesDispl
      Type(Vec),Intent(IN)                               :: x
      Type(Mat),Intent(INOUT)                            :: A,M
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr  

      PetscInt                                           :: dim      
      Call DMMeshGetDimension(MEF90DefMechCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90DefMechBilinearForm2D(snesDispl,x,A,M,flg,MEF90DefMechCtx,ierr)
      Else If (dim == 3) Then
         Call MEF90DefMechBilinearForm3D(snesDispl,x,A,M,flg,MEF90DefMechCtx,ierr)
      End If      
   End Subroutine MEF90DefMechBilinearForm

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechWork"
!!!
!!!  
!!!  MEF90DefMechWork: wraps calls to MEF90DefMechWork from m_MEF90_DefMechAssembly
!!                       since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechWork(DisplacementVec,MEF90DefMechCtx,work,ierr)
      Type(Vec),Intent(IN)                            :: DisplacementVec
      Type(MEF90DefMechCtx_Type),Intent(IN)           :: MEF90DefMechCtx
      PetscReal,Dimension(:),Pointer                  :: work
      PetscErrorCode,Intent(OUT)                      :: ierr

      PetscInt                                        :: dim      
      Call DMMeshGetDimension(MEF90DefMechCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90DefMechWork2D(DisplacementVec,MEF90DefMechCtx,work,ierr)
      Else If (dim == 3) Then
         Call MEF90DefMechWork3D(DisplacementVec,MEF90DefMechCtx,work,ierr)
      End If      
   End Subroutine MEF90DefMechWork

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechElasticEnergy"
!!!
!!!  
!!!  MEF90DefMechElasticEnergy: wraps calls to MEF90DefMechElasticEnergy from m_MEF90_DefMechAssembly
!!                       since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechElasticEnergy(x,MEF90DefMechCtx,energy,ierr)
      Type(Vec),Intent(IN)                               :: x
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,dimension(:),Pointer                     :: energy
      PetscErrorCode,Intent(OUT)                         :: ierr

      PetscInt                                        :: dim      
      Call DMMeshGetDimension(MEF90DefMechCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90DefMechElasticEnergy2D(x,MEF90DefMechCtx,energy,ierr)
      Else If (dim == 3) Then
         Call MEF90DefMechElasticEnergy3D(x,MEF90DefMechCtx,energy,ierr)
      End If      
   End Subroutine MEF90DefMechElasticEnergy
End Module m_MEF90_DefMech
