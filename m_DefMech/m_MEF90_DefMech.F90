#include "../MEF90/mef90.inc"
Module m_MEF90_DefMech
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_DefMechCtx

   Use m_MEF90_DefMechAssembly2D, &
      MEF90DefMechOperator2D     => MEF90DefMechOperator, &
      MEF90DefMechBilinearForm2D => MEF90DefMechBilinearForm      
   Use m_MEF90_DefMechAssembly3D, &
      MEF90DefMechOperator3D     => MEF90DefMechOperator, &
      MEF90DefMechBilinearForm3D => MEF90DefMechBilinearForm

   Implicit none
   
Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSetTransients"
!!!
!!!  
!!!  MEF90DefMechSetTransients:
!!!  
!!!  (c) 2012-13 Blaise Bourdin bourdin@lsu.edu
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
!!!  (c) 2012-13 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechSetForceCst(x,MEF90DefMechCtx,ierr)
      Type(Vec),Intent(IN)                            :: x
      Type(MEF90DefMechCtx_Type),Intent(IN)           :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                      :: ierr
   
      Type(MEF90DefMechGlobalOptions_Type),pointer    :: MEF90DefMechGlobalOptions
      Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90GlobalOptions
      Type(MEF90DefMechCellSetOptions_Type),pointer   :: cellSetOptions
      Type(IS)                                        :: cellSetGlobalIS,setIS,setISdof
      PetscInt,Dimension(:),Pointer                   :: setID
      PetscInt,Dimension(:),Pointer                   :: setIdx
      PetscInt                                        :: set,c,dim
      PetscReal,Dimension(:),Pointer                  :: val
      
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)
      Call DMMeshGetDimension(MEF90DefMechCtx%cellDMVect,dim,ierr);CHKERRQ(ierr)

      !!! force is cell-centered
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%CellDMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%cellDMVect,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Do c = 1,dim
            Call DMMeshISCreateISglobaldof(MEF90DefMechCtx%cellDMVect,setIS,c-1,setISdof,ierr);CHKERRQ(ierr)
            Call ISGetIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
            Allocate(val(size(setIdx)),stat=ierr)
            val = cellSetOptions%force(c)
            Call VecSetValues(x,size(setIdx),setIdx,val,INSERT_VALUES,ierr);CHKERRQ(ierr)
            DeAllocate(val)
            Call ISRestoreIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
            Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
         End Do
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(cellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(cellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call VecAssemblyBegin(x,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(x,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechSetForceCst

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSetPressureForceCst"
!!!
!!!  
!!!  MEF90DefMechSetPressureForceCst:
!!!  
!!!  (c) 2012-13 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechSetPressureForceCst(x,MEF90DefMechCtx,ierr)
      Type(Vec),Intent(IN)                            :: x
      Type(MEF90DefMechCtx_Type),Intent(IN)           :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                      :: ierr
   
      Type(MEF90DefMechGlobalOptions_Type),pointer    :: MEF90DefMechGlobalOptions
      Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90GlobalOptions
      Type(MEF90DefMechCellSetOptions_Type),pointer   :: cellSetOptions
      Type(IS)                                        :: cellSetGlobalIS,setIS,setISdof
      PetscInt,Dimension(:),Pointer                   :: setID
      PetscInt,Dimension(:),Pointer                   :: setIdx
      PetscInt                                        :: set
      PetscReal,Dimension(:),Pointer                  :: val
      
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)
      
      !!! pressureForce is cell-centered
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%cellDMScal,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%cellDMScal,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call DMMeshISCreateISglobaldof(MEF90DefMechCtx%cellDMScal,setIS,0,setISdof,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
         Allocate(val(size(setIdx)),stat=ierr)
         val = cellSetOptions%PressureForce
         Call VecSetValues(x,size(setIdx),setIdx,val,INSERT_VALUES,ierr);CHKERRQ(ierr)
         DeAllocate(val)
         Call ISRestoreIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
         Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(cellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(cellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call VecAssemblyBegin(x,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(x,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechSetPressureForceCst

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSetBoundaryDisplacementCst"
!!!
!!!  
!!!  MEF90DefMechSetBoundaryDisplacementCst:
!!!  
!!!  (c) 2012-13 Blaise Bourdin bourdin@lsu.edu
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
      PetscInt                                           :: set
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
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DMVect,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call MEF90_ISCreateCelltoVertex(MEF90DefMechCtx%DMVect,PETSC_COMM_WORLD,setIS,bcIS,ierr)
         Do c = 1, dim
            If (cellSetOptions%Has_displacementBC(c)) Then
               Call DMMeshISCreateISglobaldof(MEF90DefMechCtx%DMVect,bcIS,c-1,setISdof,ierr);CHKERRQ(ierr)
               Call ISGetIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
               Allocate(val(size(setIdx)),stat=ierr)
               val = cellSetOptions%boundaryDisplacement(c)
               Call VecSetValues(x,size(setIdx),setIdx,val,INSERT_VALUES,ierr);CHKERRQ(ierr)
               DeAllocate(val)
               Call ISRestoreIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
               Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
            End If ! cellSetOptions%Has_displacementBC
            End Do
         Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do ! set
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
      
      !!! Vertex Sets
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxVertexSetOptions(MEF90DefMechCtx%VertexSetOptionsBag(set),VertexSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DMVect,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Do c = 1, dim
            If (vertexSetOptions%Has_displacementBC(c)) Then
               Call DMMeshISCreateISglobaldof(MEF90DefMechCtx%DMVect,setIS,c-1,setISdof,ierr);CHKERRQ(ierr)
               Call ISGetIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
               Allocate(val(size(setIdx)),stat=ierr)
               val = VertexSetOptions%boundaryDisplacement(c)
               Call VecSetValues(x,size(setIdx),setIdx,val,INSERT_VALUES,ierr);CHKERRQ(ierr)
               DeAllocate(val)
               Call ISRestoreIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
               Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
            End If
         End Do
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call VecAssemblyBegin(x,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(x,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechSetBoundaryDisplacementCst

!#undef __FUNCT__
!#define __FUNCT__ "MEF90DefMechOperator"
!!!!
!!!!  
!!!!  MEF90DefMechOperator:
!!!!  
!!!!  (c) 2012-13 Blaise Bourdin bourdin@lsu.edu
!!!!
!   Subroutine MEF90DefMechOperator(snesTemp,x,residual,MEF90DefMechCtx,ierr)
!      Type(SNES),Intent(IN)                              :: snesTemp
!      Type(Vec),Intent(IN)                               :: x
!      Type(Vec),Intent(INOUT)                            :: residual
!      Type(MEF90DefMechCtx_Type),Intent(IN)             :: MEF90DefMechCtx
!      PetscErrorCode,Intent(OUT)                         :: ierr
!      
!      PetscInt                                           :: dim      
!      Call DMMeshGetDimension(MEF90DefMechCtx%DM,dim,ierr);CHKERRQ(ierr)
!      If (dim == 2) Then
!         Call MEF90DefMechOperator2D(snesTemp,x,residual,MEF90DefMechCtx,ierr)
!      Else If (dim == 3) Then
!         Call MEF90DefMechOperator3D(snesTemp,x,residual,MEF90DefMechCtx,ierr)
!      End If      
!   End Subroutine MEF90DefMechOperator
!   
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

      PetscInt                                           :: dim      
      Call DMMeshGetDimension(MEF90DefMechCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90DefMechBilinearForm2D(snesDispl,x,A,M,flg,MEF90DefMechCtx,ierr)
      Else If (dim == 3) Then
         Call MEF90DefMechBilinearForm3D(snesDispl,x,A,M,flg,MEF90DefMechCtx,ierr)
      End If      
   End Subroutine MEF90DefMechBilinearForm

!#undef __FUNCT__
!#define __FUNCT__ "MEF90DefMechEnergy"
!!!!
!!!!  
!!!!  MEF90DefMechEnergy:
!!!!  
!!!!  (c) 2012-13 Blaise Bourdin bourdin@lsu.edu
!!!!
!   Subroutine MEF90DefMechEnergy(temperatureVec,t,MEF90DefMechCtx,energy,work,ierr)
!      Type(Vec),Intent(IN)                            :: temperatureVec
!      PetscReal,Intent(IN)                            :: t
!      Type(MEF90DefMechCtx_Type),Intent(IN)          :: MEF90DefMechCtx
!      PetscReal,Dimension(:),Pointer                  :: energy,work
!      PetscErrorCode,Intent(OUT)                      :: ierr
!
!      PetscInt                                           :: dim      
!      Call DMMeshGetDimension(MEF90DefMechCtx%DM,dim,ierr);CHKERRQ(ierr)
!      If (dim == 2) Then
!         Call MEF90DefMechEnergy2D(temperatureVec,t,MEF90DefMechCtx,energy,work,ierr)
!      Else If (dim == 3) Then
!         Call MEF90DefMechEnergy3D(temperatureVec,t,MEF90DefMechCtx,energy,work,ierr)
!      End If      
!   End Subroutine MEF90DefMechEnergy


End Module m_MEF90_DefMech
