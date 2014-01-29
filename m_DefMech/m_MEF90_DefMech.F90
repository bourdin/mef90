#include "../MEF90/mef90.inc"
Module m_MEF90_DefMech
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_DefMechCtx

   Use m_MEF90_DefMechAssembly2D, &
      MEF90DefMechOperatorDisplacement2D     => MEF90DefMechOperatorDisplacement,      &
      MEF90DefMechBilinearFormDisplacement2D => MEF90DefMechBilinearFormDisplacement,  &     
      MEF90DefMechWork2D                     => MEF90DefMechWork,                      &
      MEF90DefMechElasticEnergy2D            => MEF90DefMechElasticEnergy,             &
      MEF90DefMechOperatorDamage2D           => MEF90DefMechOperatorDamage,            &
      MEF90DefMechBilinearFormDamage2D       => MEF90DefMechBilinearFormDamage,        &
      MEF90DefMechSurfaceEnergy2D            => MEF90DefMechSurfaceEnergy     
   Use m_MEF90_DefMechAssembly3D, &
      MEF90DefMechOperatorDisplacement3D     => MEF90DefMechOperatorDisplacement,      &
      MEF90DefMechBilinearFormDisplacement3D => MEF90DefMechBilinearFormDisplacement,  &     
      MEF90DefMechWork3D                     => MEF90DefMechWork,                      &
      MEF90DefMechElasticEnergy3D            => MEF90DefMechElasticEnergy,             &
      MEF90DefMechOperatorDamage3D           => MEF90DefMechOperatorDamage,            &
      MEF90DefMechBilinearFormDamage3D       => MEF90DefMechBilinearFormDamage,        &
      MEF90DefMechSurfaceEnergy3D            => MEF90DefMechSurfaceEnergy     

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
      Type(Vec)                                       :: localVec

      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)

      If (Associated(MEF90DefMechCtx%force)) Then
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
      End If
      If (Associated(MEF90DefMechCtx%pressureForce)) Then
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
      End If
      Select case (MEF90DefMechGlobalOptions%boundaryDisplacementScaling)
         Case (MEF90Scaling_File)
            Call DMGetLocalVector(MEF90DefMechCtx%DMVect,localVec,ierr);CHKERRQ(ierr)
            Call VecLoadExodusVertex(MEF90DefMechCtx%DMVect,localVec,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
                                     MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%boundaryDisplacementOffset,ierr);CHKERRQ(ierr)
            Call DMLocalToGlobalBegin(MEF90DefMechCtx%DMVect,localVec,INSERT_VALUES,MEF90DefMechCtx%boundaryDisplacement,ierr);CHKERRQ(ierr)
            Call DMLocalToGlobalEnd(MEF90DefMechCtx%DMVect,localVec,INSERT_VALUES,MEF90DefMechCtx%boundaryDisplacement,ierr);CHKERRQ(ierr)
            Call DMRestoreLocalVector(MEF90DefMechCtx%DMVect,localVec,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_Linear)
            Call MEF90DefMechSetboundaryDisplacementCst(MEF90DefMechCtx%boundaryDisplacement,MEF90DefMechCtx,ierr)
            Call VecScale(MEF90DefMechCtx%boundaryDisplacement,time,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_CST)
            Call MEF90DefMechSetboundaryDisplacementCst(MEF90DefMechCtx%boundaryDisplacement,MEF90DefMechCtx,ierr)
      End Select
      Select case (MEF90DefMechGlobalOptions%boundaryDamageScaling)
         Case (MEF90Scaling_File)
            Call DMGetLocalVector(MEF90DefMechCtx%DMScal,localVec,ierr);CHKERRQ(ierr)
            Call VecLoadExodusVertex(MEF90DefMechCtx%DMScal,localVec,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
                                     MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%boundaryDamageOffset,ierr);CHKERRQ(ierr)
            Call DMLocalToGlobalBegin(MEF90DefMechCtx%DMScal,localVec,INSERT_VALUES,MEF90DefMechCtx%boundaryDamage,ierr);CHKERRQ(ierr)
            Call DMLocalToGlobalEnd(MEF90DefMechCtx%DMScal,localVec,INSERT_VALUES,MEF90DefMechCtx%boundaryDamage,ierr);CHKERRQ(ierr)
            Call DMRestoreLocalVector(MEF90DefMechCtx%DMScal,localVec,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_CST)
            Call MEF90DefMechSetboundaryDamageCst(MEF90DefMechCtx%boundaryDamage,MEF90DefMechCtx,ierr)
         Case (MEF90Scaling_Linear)
            Write(*,*) __FUNCT__,": linear scaling of damage variable does not make any sense."
            STOP
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
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
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
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
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
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)   
         Do c = 1, dim
            If (cellSetOptions%Has_displacementBC(c)) Then
               Call DMMeshGetStratumIS(MEF90DefMechCtx%DMVect,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
               Call MEF90ISCreateCelltoVertex(MEF90DefMechCtx%DMVect,PETSC_COMM_WORLD,setIS,bcIS,ierr)
               Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
               Allocate(val(nval),stat=ierr)
               val =  cellSetOptions%boundaryDisplacement(c)
               Call MEF90VecSetValuesISdof(MEF90DefMechCtx%DMVect,x,val,bcIS,c,INSERT_VALUES,ierr)
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
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxVertexSetOptions(MEF90DefMechCtx%vertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)   
         Do c = 1, dim
            If (vertexSetOptions%Has_displacementBC(c)) Then
               Call DMMeshGetStratumIS(MEF90DefMechCtx%DMVect,'Vertex Sets',setID(set),bcIS,ierr);CHKERRQ(iErr)
               Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
               Allocate(val(nval),stat=ierr)
               val =  vertexSetOptions%boundaryDisplacement(c)
               Call MEF90VecSetValuesISdof(MEF90DefMechCtx%DMVect,x,val,bcIS,c,INSERT_VALUES,ierr)
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
#define __FUNCT__ "MEF90DefMechSetBoundaryDamageCst"
!!!
!!!  
!!!  MEF90DefMechSetBoundaryDamageCst:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechSetBoundaryDamageCst(x,MEF90DefMechCtx,ierr)
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
      PetscInt                                           :: cell
      
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)
      
      !!! boundaryDamage is Vertex-centered
      !!! We first set the boundary values inherited from cell sets, then that of vertex sets
      !!!
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)   
         If (cellSetOptions%Has_damageBC) Then
            Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
            Call MEF90ISCreateCelltoVertex(MEF90DefMechCtx%DMScal,PETSC_COMM_WORLD,setIS,bcIS,ierr)
            Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
            Allocate(val(nval),stat=ierr)
            val =  cellSetOptions%boundaryDamage
            Call MEF90VecSetValuesISdof(MEF90DefMechCtx%DMScal,x,val,bcIS,1,INSERT_VALUES,ierr)
            DeAllocate(Val)
            Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If! cellSetOptions%Has_BC
      End Do ! set
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
      
      !!! Vertex Sets
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DM,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxVertexSetOptions(MEF90DefMechCtx%vertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)   
         If (vertexSetOptions%Has_damageBC) Then
            Call DMMeshGetStratumIS(MEF90DefMechCtx%DMVect,'Vertex Sets',setID(set),bcIS,ierr);CHKERRQ(iErr)
            Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
            Allocate(val(nval),stat=ierr)
            val =  vertexSetOptions%boundaryDamage
            Call MEF90VecSetValuesISdof(MEF90DefMechCtx%DMScal,x,val,bcIS,1,INSERT_VALUES,ierr)
            DeAllocate(Val)
            Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
         End If! cellSetOptions%Has_BC
      End Do ! set
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call VecAssemblyBegin(x,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(x,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechSetBoundaryDamageCst

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
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)   
         Do c = 1, dim
            If (cellSetOptions%Has_displacementBC(c)) Then
               Call DMMeshGetStratumIS(MEF90DefMechCtx%DMVect,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
               Call MEF90ISCreateCelltoVertex(MEF90DefMechCtx%DMVect,PETSC_COMM_WORLD,setIS,bcIS,ierr)
               Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
               Allocate(boundaryDisplacementPtr(nval),stat=ierr)
               Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMVect,MEF90DefMechCtx%boundaryDisplacement,boundaryDisplacementPtr,bcIS,c,ierr)
               Call MEF90VecSetValuesISdof(MEF90DefMechCtx%DMVect,x,boundaryDisplacementPtr,bcIS,c,INSERT_VALUES,ierr)
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
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxVertexSetOptions(MEF90DefMechCtx%vertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)   
         Do c = 1, dim
            If (vertexSetOptions%Has_displacementBC(c)) Then
               Call DMMeshGetStratumIS(MEF90DefMechCtx%DMVect,'Vertex Sets',setID(set),bcIS,ierr);CHKERRQ(iErr)
               Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
               Allocate(boundaryDisplacementPtr(nval),stat=ierr)
               Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMVect,MEF90DefMechCtx%boundaryDisplacement,boundaryDisplacementPtr,bcIS,c,ierr)
               Call MEF90VecSetValuesISdof(MEF90DefMechCtx%DMVect,x,boundaryDisplacementPtr,bcIS,c,INSERT_VALUES,ierr)
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
#define __FUNCT__ "MEF90DefMechUpdateboundaryDamage"
!!!
!!!  
!!!  MEF90DefMechUpdateboundaryDamage: Update the solution vector x with the boundary Displacement from 
!!!                                          MEF90DefMechCtx%BoundaryDisplacement
!!!  
!!!  (c) 2013 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechUpdateboundaryDamage(x,MEF90DefMechCtx,ierr)
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
      PetscInt                                           :: set,nval
      PetscReal,Dimension(:),Pointer                     :: boundaryDamagePtr
      Type(MEF90Element_Type)                            :: elemType
      
   
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)
      
      !!!
      !!! cell set Displacement first, followed vertex sets
      !!!
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)   
         If (cellSetOptions%Has_damageBC) Then
            Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
            Call MEF90ISCreateCelltoVertex(MEF90DefMechCtx%DMScal,PETSC_COMM_WORLD,setIS,bcIS,ierr)
            Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
            Allocate(boundaryDamagePtr(nval),stat=ierr)
            Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMScal,MEF90DefMechCtx%boundaryDamage,boundaryDamagePtr,bcIS,1,ierr)
            Call MEF90VecSetValuesISdof(MEF90DefMechCtx%DMScal,x,boundaryDamagePtr,bcIS,1,INSERT_VALUES,ierr)
            DeAllocate(boundaryDamagePtr)
            Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If! cellSetOptions%Has_BC
      End Do ! set
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)

      !!! Vertex sets
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DM,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxVertexSetOptions(MEF90DefMechCtx%vertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)   
         If (vertexSetOptions%Has_damageBC) Then
            Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Vertex Sets',setID(set),bcIS,ierr);CHKERRQ(iErr)
            Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
            Allocate(boundaryDamagePtr(nval),stat=ierr)
            Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMScal,MEF90DefMechCtx%boundaryDamage,boundaryDamagePtr,bcIS,1,ierr)
            Call MEF90VecSetValuesISdof(MEF90DefMechCtx%DM,x,boundaryDamagePtr,bcIS,1,INSERT_VALUES,ierr)
            DeAllocate(boundaryDamagePtr)
            Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
         End If! cellSetOptions%Has_BC
      End Do ! set
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)

      Call VecAssemblyBegin(x,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(x,ierr);CHKERRQ(ierr)
End Subroutine MEF90DefMechUpdateboundaryDamage

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDisplacement"
!!!
!!!  
!!!  MEF90DefMechOperatorDisplacement: wraps calls to MEF90DefMechOperatorDisplacement from m_MEF90_DefMechAssembly
!!!                        since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechOperatorDisplacement(snesTemp,x,residual,MEF90DefMechCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesTemp
      Type(Vec),Intent(IN)                               :: x
      Type(Vec),Intent(INOUT)                            :: residual
      Type(MEF90DefMechCtx_Type),Intent(IN)             :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      PetscInt                                           :: dim      
      Call DMMeshGetDimension(MEF90DefMechCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90DefMechOperatorDisplacement2D(snesTemp,x,residual,MEF90DefMechCtx,ierr)
      Else If (dim == 3) Then
         Call MEF90DefMechOperatorDisplacement3D(snesTemp,x,residual,MEF90DefMechCtx,ierr)
      End If      
   End Subroutine MEF90DefMechOperatorDisplacement
   
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDisplacement"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDisplacement: wraps calls to MEF90DefMechBilinearFormDisplacement from m_MEF90_DefMechAssembly
!!!                            since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechBilinearFormDisplacement(snesDispl,x,A,M,flg,MEF90DefMechCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesDispl
      Type(Vec),Intent(IN)                               :: x
      Type(Mat),Intent(INOUT)                            :: A,M
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr  

      PetscInt                                           :: dim      
      Call DMMeshGetDimension(MEF90DefMechCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90DefMechBilinearFormDisplacement2D(snesDispl,x,A,M,flg,MEF90DefMechCtx,ierr)
      Else If (dim == 3) Then
         Call MEF90DefMechBilinearFormDisplacement3D(snesDispl,x,A,M,flg,MEF90DefMechCtx,ierr)
      End If      
   End Subroutine MEF90DefMechBilinearFormDisplacement

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

      PetscInt                                           :: dim      
      Call DMMeshGetDimension(MEF90DefMechCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90DefMechElasticEnergy2D(x,MEF90DefMechCtx,energy,ierr)
      Else If (dim == 3) Then
         Call MEF90DefMechElasticEnergy3D(x,MEF90DefMechCtx,energy,ierr)
      End If      
   End Subroutine MEF90DefMechElasticEnergy

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamage"
!!!
!!!  
!!!  MEF90DefMechOperatorDamage: wraps calls to MEF90DefMechOperatorDamage from m_MEF90_DefMechAssembly
!!!                        since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechOperatorDamage(snesTemp,x,residual,MEF90DefMechCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesTemp
      Type(Vec),Intent(IN)                               :: x
      Type(Vec),Intent(INOUT)                            :: residual
      Type(MEF90DefMechCtx_Type),Intent(IN)             :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      PetscInt                                           :: dim      
      Call DMMeshGetDimension(MEF90DefMechCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90DefMechOperatorDamage2D(snesTemp,x,residual,MEF90DefMechCtx,ierr)
      Else If (dim == 3) Then
         Call MEF90DefMechOperatorDamage3D(snesTemp,x,residual,MEF90DefMechCtx,ierr)
      End If      
   End Subroutine MEF90DefMechOperatorDamage
   
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamage"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamage: wraps calls to MEF90DefMechBilinearFormDamage from m_MEF90_DefMechAssembly
!!!                            since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechBilinearFormDamage(snesDispl,x,A,M,flg,MEF90DefMechCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesDispl
      Type(Vec),Intent(IN)                               :: x
      Type(Mat),Intent(INOUT)                            :: A,M
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr  

      PetscInt                                           :: dim      
      Call DMMeshGetDimension(MEF90DefMechCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90DefMechBilinearFormDamage2D(snesDispl,x,A,M,flg,MEF90DefMechCtx,ierr)
      Else If (dim == 3) Then
         Call MEF90DefMechBilinearFormDamage3D(snesDispl,x,A,M,flg,MEF90DefMechCtx,ierr)
      End If      
   End Subroutine MEF90DefMechBilinearFormDamage

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSurfaceEnergy"
!!!
!!!  
!!!  MEF90DefMechSurfaceEnergy: wraps calls to MEF90DefMechSurfaceEnergy from m_MEF90_DefMechAssembly
!!                       since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechSurfaceEnergy(x,MEF90DefMechCtx,energy,ierr)
      Type(Vec),Intent(IN)                               :: x
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,dimension(:),Pointer                     :: energy
      PetscErrorCode,Intent(OUT)                         :: ierr

      PetscInt                                           :: dim      
      Call DMMeshGetDimension(MEF90DefMechCtx%DM,dim,ierr);CHKERRQ(ierr)
      If (dim == 2) Then
         Call MEF90DefMechSurfaceEnergy2D(x,MEF90DefMechCtx,energy,ierr)
      Else If (dim == 3) Then
         Call MEF90DefMechSurfaceEnergy3D(x,MEF90DefMechCtx,energy,ierr)
      End If      
   End Subroutine MEF90DefMechSurfaceEnergy

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechViewEXO"
!!!
!!!  
!!!  MEF90DefMechViewEXO: Save all fields in a MEF90DefMechCtx_Type in an exodus file
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechViewEXO(MEF90DefMechCtx,step,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscInt,Intent(IN)                                :: step
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(Vec)                                          :: localVec
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions

      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)
      If ((MEF90DefMechGlobalOptions%boundaryDisplacementOffset > 0) .AND. &
          (MEF90DefMechGlobalOptions%boundaryDisplacementOffset /= MEF90DefMechGlobalOptions%displacementOffset)) Then
         Call DMGetLocalVector(MEF90DefMechCtx%DMVect,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalBegin(MEF90DefMechCtx%DMVect,MEF90DefMechCtx%boundaryDisplacement,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalEnd(MEF90DefMechCtx%DMVect,MEF90DefMechCtx%boundaryDisplacement,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call VecViewExodusVertex(MEF90DefMechCtx%DMVect,localVec,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
                                  MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%boundaryDisplacementOffset,ierr);CHKERRQ(ierr)
         Call DMRestoreLocalVector(MEF90DefMechCtx%DMVect,localVec,ierr);CHKERRQ(ierr)
      End If

      If (MEF90DefMechGlobalOptions%displacementOffset > 0) Then
         Call DMGetLocalVector(MEF90DefMechCtx%DMVect,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalBegin(MEF90DefMechCtx%DMVect,MEF90DefMechCtx%Displacement,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalEnd(MEF90DefMechCtx%DMVect,MEF90DefMechCtx%Displacement,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call VecViewExodusVertex(MEF90DefMechCtx%DMVect,localVec,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
                                  MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%displacementOffset,ierr);CHKERRQ(ierr)
         Call DMRestoreLocalVector(MEF90DefMechCtx%DMVect,localVec,ierr);CHKERRQ(ierr)
      End If

      If (MEF90DefMechGlobalOptions%damageOffset > 0) Then
         Call DMGetLocalVector(MEF90DefMechCtx%DMScal,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalBegin(MEF90DefMechCtx%DMScal,MEF90DefMechCtx%damage,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalEnd(MEF90DefMechCtx%DMScal,MEF90DefMechCtx%damage,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call VecViewExodusVertex(MEF90DefMechCtx%DMScal,localVec,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
                                  MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%damageOffset,ierr);CHKERRQ(ierr)
         Call DMRestoreLocalVector(MEF90DefMechCtx%DMScal,localVec,ierr);CHKERRQ(ierr)
      End If

      If ((MEF90DefMechGlobalOptions%boundaryDamageOffset > 0) .AND. &
          (MEF90DefMechGlobalOptions%boundaryDamageOffset /= MEF90DefMechGlobalOptions%damageOffset))Then
         Call DMGetLocalVector(MEF90DefMechCtx%DMScal,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalBegin(MEF90DefMechCtx%DMScal,MEF90DefMechCtx%boundaryDamage,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalEnd(MEF90DefMechCtx%DMScal,MEF90DefMechCtx%boundaryDamage,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call VecViewExodusVertex(MEF90DefMechCtx%DMScal,localVec,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
                                  MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%boundaryDamageOffset,ierr);CHKERRQ(ierr)
         Call DMRestoreLocalVector(MEF90DefMechCtx%DMScal,localVec,ierr);CHKERRQ(ierr)
      End If

      If (MEF90DefMechGlobalOptions%temperatureOffset > 0) Then
         Call DMGetLocalVector(MEF90DefMechCtx%DMScal,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalBegin(MEF90DefMechCtx%DMScal,MEF90DefMechCtx%temperature,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalEnd(MEF90DefMechCtx%DMScal,MEF90DefMechCtx%temperature,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call VecViewExodusVertex(MEF90DefMechCtx%DMScal,localVec,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
                                  MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%temperatureOffset,ierr);CHKERRQ(ierr)
         Call DMRestoreLocalVector(MEF90DefMechCtx%DMScal,localVec,ierr);CHKERRQ(ierr)
      End If

      If (MEF90DefMechGlobalOptions%forceOffset > 0) Then
         Call DMGetLocalVector(MEF90DefMechCtx%cellDMVect,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalBegin(MEF90DefMechCtx%cellDMVect,MEF90DefMechCtx%Force,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalEnd(MEF90DefMechCtx%cellDMVect,MEF90DefMechCtx%Force,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call VecViewExodusCell(MEF90DefMechCtx%cellDMVect,localVec,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
                                MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%forceOffset,ierr);CHKERRQ(ierr)
         Call DMRestoreLocalVector(MEF90DefMechCtx%cellDMVect,localVec,ierr);CHKERRQ(ierr)
      End If

      If (MEF90DefMechGlobalOptions%pressureForceOffset > 0) Then
         Call DMGetLocalVector(MEF90DefMechCtx%cellDMScal,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalBegin(MEF90DefMechCtx%cellDMScal,MEF90DefMechCtx%pressureForce,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalEnd(MEF90DefMechCtx%cellDMScal,MEF90DefMechCtx%pressureForce,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call VecViewExodusCell(MEF90DefMechCtx%cellDMScal,localVec,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
                                MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%pressureForceOffset,ierr);CHKERRQ(ierr)
         Call DMRestoreLocalVector(MEF90DefMechCtx%cellDMScal,localVec,ierr);CHKERRQ(ierr)
      End If

      If (MEF90DefMechGlobalOptions%plasticStrainOffset > 0) Then
         Call DMGetLocalVector(MEF90DefMechCtx%cellDMMatS,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalBegin(MEF90DefMechCtx%cellDMMatS,MEF90DefMechCtx%plasticStrain,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalEnd(MEF90DefMechCtx%cellDMMatS,MEF90DefMechCtx%plasticStrain,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call VecViewExodusCell(MEF90DefMechCtx%cellDMMatS,localVec,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
                                MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%plasticStrainOffset,ierr);CHKERRQ(ierr)
         Call DMRestoreLocalVector(MEF90DefMechCtx%cellDMMatS,localVec,ierr);CHKERRQ(ierr)
      End If
      If (MEF90DefMechCtx%MEF90Ctx%rank == 0) Then
         Call EXUPDA(MEF90DefMechCtx%MEF90Ctx%fileExoUnit,ierr)
      End If
   End Subroutine MEF90DefMechViewEXO
   
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechFormatEXO"
!!!
!!!  
!!!  MEF90DefMechFormatEXO:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechFormatEXO(MEF90DefMechCtx,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr

      Character(len=MXSTLN),Dimension(:),Pointer         :: nameG,nameV,nameC
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
      Integer                                            :: dim,numfield

      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)
      Call DMMeshGetDimension(MEF90DefMechCtx%DM,dim,ierr);CHKERRQ(ierr)
      Allocate(nameG(0))
      !nameG(1) = "Elastic Energy"
      !nameG(2) = "Work"
      !nameG(3) = "Surface Energy"
      !nameG(4) = "Total Energy"
   
      numfield = max(MEF90DefMechGlobalOptions%displacementOffset+dim-1, &
                     MEF90DefMechGlobalOptions%damageOffset,&
                     MEF90DefMechGlobalOptions%boundaryDisplacementOffset+dim-1,&
                     MEF90DefMechGlobalOptions%boundaryDamageOffset,&
                     MEF90DefMechGlobalOptions%temperatureOffset)
      Allocate(nameV(numfield))
      nameV = "empty"
      If ((MEF90DefMechGlobalOptions%boundaryDisplacementOffset > 0) .AND. &
          (MEF90DefMechGlobalOptions%boundaryDisplacementOffset /= MEF90DefMechGlobalOptions%displacementOffset)) Then
         nameV(MEF90DefMechGlobalOptions%boundaryDisplacementOffset+0)    = "Boundary_Displacement_X"
         nameV(MEF90DefMechGlobalOptions%boundaryDisplacementOffset+1)    = "Boundary_Displacement_Y"
         If (dim == 3) Then
            nameV(MEF90DefMechGlobalOptions%boundaryDisplacementOffset+2) = "Boundary_Displacement_Z"
         End If
      End If
      If (MEF90DefMechGlobalOptions%displacementOffset > 0) Then
         nameV(MEF90DefMechGlobalOptions%displacementOffset+0)            = "Displacement_X"
         nameV(MEF90DefMechGlobalOptions%displacementOffset+1)            = "Displacement_Y"
         If (dim == 3) Then
            nameV(MEF90DefMechGlobalOptions%displacementOffset+2)         = "Displacement_Z"
         End If
      End If
      If (MEF90DefMechGlobalOptions%damageOffset > 0) Then
         nameV(MEF90DefMechGlobalOptions%damageOffset)                    = "Damage"
      End If
      If ((MEF90DefMechGlobalOptions%boundaryDamageOffset > 0) .AND. &
          (MEF90DefMechGlobalOptions%boundaryDamageOffset /= MEF90DefMechGlobalOptions%damageOffset))Then
         nameV(MEF90DefMechGlobalOptions%boundaryDamageOffset)            = "Boundary_Damage"
      End If
      If (MEF90DefMechGlobalOptions%temperatureOffset > 0) Then
         nameV(MEF90DefMechGlobalOptions%temperatureOffset)               = "Temperature"
      End If
                     
      numfield = max(MEF90DefMechGlobalOptions%forceOffset+dim-1,&
                     MEF90DefMechGlobalOptions%pressureForceOffset,&
                     MEF90DefMechGlobalOptions%StressOffset+(dim*(dim+1))/2-1,&
                     MEF90DefMechGlobalOptions%plasticStrainOffset+(dim*(dim+1))/2-1)
      Allocate(nameC(numfield))
      nameC = "empty"
      If (MEF90DefMechGlobalOptions%forceOffset > 0) Then
         nameC(MEF90DefMechGlobalOptions%forceOffset+0)                 = "Force_X"
         nameC(MEF90DefMechGlobalOptions%forceOffset+1)                 = "Force_Y"
         If (dim == 3) Then
            nameC(MEF90DefMechGlobalOptions%forceOffset+2)              = "Force_Z"
         End If
      End If
      
      If (MEF90DefMechGlobalOptions%pressureForceOffset > 0) Then
         nameC(MEF90DefMechGlobalOptions%pressureForceOffset)           = "Pressure_Force"
      End If
      If (MEF90DefMechGlobalOptions%stressOffset > 0) Then
         If (dim == 2) Then
            nameC(MEF90DefMechGlobalOptions%stressOffset+0)             = "Stress_XX"
            nameC(MEF90DefMechGlobalOptions%stressOffset+1)             = "Stress_YY"
            nameC(MEF90DefMechGlobalOptions%stressOffset+2)             = "Stress_XY"
         Else
            nameC(MEF90DefMechGlobalOptions%stressOffset+0)             = "Stress_XX"
            nameC(MEF90DefMechGlobalOptions%stressOffset+1)             = "Stress_YY"
            nameC(MEF90DefMechGlobalOptions%stressOffset+2)             = "Stress_ZZ"
            nameC(MEF90DefMechGlobalOptions%stressOffset+3)             = "Stress_YZ"
            nameC(MEF90DefMechGlobalOptions%stressOffset+4)             = "Stress_XZ"
            nameC(MEF90DefMechGlobalOptions%stressOffset+5)             = "Stress_XY"
         End If
      End If
      If (MEF90DefMechGlobalOptions%plasticStrainOffset > 0) Then
         If (dim == 2) Then
            nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+0)      = "plasticStrain_XX"
            nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+1)      = "plasticStrain_YY"
            nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+2)      = "plasticStrain_XY"
         Else
            nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+0)      = "plasticStrain_XX"
            nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+1)      = "plasticStrain_YY"
            nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+2)      = "plasticStrain_ZZ"
            nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+3)      = "plasticStrain_YZ"
            nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+4)      = "plasticStrain_XZ"
            nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+5)      = "plasticStrain_XY"
         End If
      End If
      
      Call MEF90EXOFormat(MEF90DefMechCtx%MEF90Ctx%fileEXOUNIT,nameG,nameC,nameV,ierr)
      DeAllocate(nameG)
      DeAllocate(nameV)
      DeAllocate(nameC)
   End Subroutine MEF90DefMechFormatEXO
   
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCreateSolversDisp"
!!!
!!!  
!!!  MEF90DefMechCreateSolversDisp:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechCreateSolversDisp(MEF90DefMechCtx,snesDisp,residual,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      Type(SNES),Intent(OUT)                             :: snesDisp
      Type(Vec),Intent(IN)                               :: residual
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
      Type(Mat)                                          :: matDisp
      Type(SectionReal)                                  :: coordSec
      Type(Vec)                                          :: CoordVec
      PetscReal,Dimension(:,:),Pointer                   :: CoordPtr
      Type(VecScatter)                                   :: ScatterSecToVec
      Type(MatNullSpace)                                 :: nspDisp
      Type(Vec)                                          :: residualDisp
      Type(KSP)                                          :: kspDisp
      Type(PC)                                           :: pcDisp
      PetscReal                                          :: atol,rtol,dtol
      PetscReal,Dimension(:),Pointer                     :: CoordPCPtr
      PetscInt                                           :: dim
      
      Call DMMeshGetDimension(MEF90DefMechCtx%DMVect,dim,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)
      Call DMCreateMatrix(MEF90DefMechCtx%DMVect,MATAIJ,matDisp,iErr);CHKERRQ(iErr)
      Call MatSetOptionsPrefix(matDisp,"Disp_",ierr);CHKERRQ(ierr)
      Call MatSetOption(matDisp,MAT_SPD,PETSC_TRUE,ierr);CHKERRQ(ierr)
      Call MatSetOption(matDisp,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr);CHKERRQ(ierr)
      Call MatSetOption(matDisp,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr);CHKERRQ(ierr)
      If (MEF90DefMechGlobalOptions%addDisplacementNullSpace) Then
         Call DMMeshGetSectionReal(MEF90DefMechCtx%DMVect,'coordinates',coordSec,ierr);CHKERRQ(ierr)
         Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMVect,coordSec,ScatterSecToVec,ierr);CHKERRQ(ierr)
         Call DMCreateGlobalVector(MEF90DefMechCtx%DMVect,coordVec,ierr)
         Call SectionRealToVec(coordSec,ScatterSecToVec,SCATTER_FORWARD,coordVec,ierr);CHKERRQ(ierr)
         Call MatNullSpaceCreateRigidBody(coordVec,nspDisp,ierr);CHKERRQ(ierr)
         Call MatSetNearNullSpace(matDisp,nspDisp,ierr);CHKERRQ(ierr)
         !!!Call MatSetNullSpace(matDisp,nspDisp,ierr);CHKERRQ(ierr)
         Call MatNullSpaceDestroy(nspDisp,ierr);CHKERRQ(ierr)
         Call SectionRealDestroy(coordSec,ierr);CHKERRQ(ierr)
         Call VecDestroy(coordVec,ierr);CHKERRQ(ierr)
         Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
      End If

      Call MatSetFromOptions(matDisp,ierr);CHKERRQ(ierr)

      If (MEF90DefMechGlobalOptions%mode == MEF90DefMech_ModeQuasiStatic) Then
         Call SNESCreate(PETSC_COMM_WORLD,snesDisp,ierr);CHKERRQ(ierr)
         Call SNESSetApplicationContext(snesDisp,MEF90DefMechCtx,ierr);CHKERRQ(ierr)
         Call SNESSetDM(snesDisp,MEF90DefMechCtx%DMVect,ierr);CHKERRQ(ierr)
         Call SNESSetOptionsPrefix(snesDisp,'Disp_',ierr);CHKERRQ(ierr)
         Call SNESSetType(snesDisp,SNESKSPONLY,ierr);CHKERRQ(ierr)

         Call SNESSetFunction(snesDisp,residual,MEF90DefMechOperatorDisplacement,MEF90DefMechCtx,ierr);CHKERRQ(ierr)
         Call SNESSetJacobian(snesDisp,matDisp,matDisp,MEF90DefMechBilinearFormDisplacement,MEF90DefMechCtx,ierr);CHKERRQ(ierr)
         !atol = 1.0D-10
         !rtol = 1.0D-10
         !Call SNESSetTolerances(snesDisp,atol,PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_INTEGER,PETSC_DEFAULT_INTEGER,ierr);CHKERRQ(ierr)
         Call SNESSetFromOptions(snesDisp,ierr);CHKERRQ(ierr)

         !!! 
         !!! Set some KSP options
         !!!
         Call SNESGetKSP(snesDisp,kspDisp,ierr);CHKERRQ(ierr)
         Call KSPSetType(kspDisp,KSPCG,ierr);CHKERRQ(ierr)
         Call KSPSetInitialGuessNonzero(kspDisp,PETSC_TRUE,ierr);CHKERRQ(ierr)
         rtol = 1.0D-8
         atol = 1.0D-8
         dtol = 1.0D+10
         Call KSPSetTolerances(kspDisp,rtol,atol,dtol,PETSC_DEFAULT_INTEGER,ierr);CHKERRQ(ierr)
         Call KSPSetFromOptions(kspDisp,ierr);CHKERRQ(ierr)
      End If

      !!! set coordinates in PC for GAMG
      !!! For some reason, this makes gamg convergence worse, when the null space is specified.
      !!! Will investigate later
      !Call KSPGetPC(kspDisp,pcDisp,ierr);CHKERRQ(ierr)
      !Call DMMeshGetCoordinatesF90(MEF90DefMechCtx%DMVect,coordPtr,ierr);CHKERRQ(ierr)
      !Allocate(coordPCPtr(size(CoordPtr)))
      !coordPCPtr = reshape(transpose(coordPtr),[size(CoordPtr)])
      !coordPCPtr = reshape((coordPtr),[size(CoordPtr)])
      !Call PCSetCoordinates(pcDisp,dim,size(coordPtr),coordPCPtr,ierr);CHKERRQ(ierr)
      !DeAllocate(coordPCPtr)
      !Call DMMeshRestoreCoordinatesF90(MEF90DefMechCtx%DMVect,coordPtr,ierr);CHKERRQ(ierr)
      !Call PCSetFromOptions(pcDisp,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechCreateSolversDisp

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCreateSolversDamage"
!!!
!!!  
!!!  MEF90DefMechCreateSolversDamage:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechCreateSolversDamage(MEF90DefMechCtx,snesDamage,residual,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      Type(SNES),Intent(OUT)                             :: snesDamage
      Type(Vec),Intent(IN)                               :: residual
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
      Type(Mat)                                          :: matDamage
      Type(SectionReal)                                  :: CoordSec
      Type(Vec)                                          :: CoordVec
      PetscReal,Dimension(:,:),Pointer                   :: CoordPtr
      Type(VecScatter)                                   :: ScatterSecToVec
      Type(Vec)                                          :: residualDamage
      Type(KSP)                                          :: kspDamage
      Type(PC)                                           :: pcDamage
      PetscReal                                          :: atol,rtol,dtol
      PetscReal,Dimension(:),Pointer                     :: CoordPCPtr
      PetscInt                                           :: dim
      Type(Vec)                                          :: LB,UB
      
      Call DMMeshGetDimension(MEF90DefMechCtx%DMScal,dim,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)
      Call DMCreateMatrix(MEF90DefMechCtx%DMScal,MATAIJ,matDamage,iErr);CHKERRQ(iErr)
      Call MatSetOptionsPrefix(matDamage,"damage_",ierr);CHKERRQ(ierr)
      Call MatSetOption(matDamage,MAT_SPD,PETSC_TRUE,ierr);CHKERRQ(ierr)
      Call MatSetOption(matDamage,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr);CHKERRQ(ierr)
      Call MatSetOption(matDamage,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr);CHKERRQ(ierr)
      Call MatSetFromOptions(matDamage,ierr);CHKERRQ(ierr)

      If (MEF90DefMechGlobalOptions%mode == MEF90DefMech_ModeQuasiStatic) Then
         Call SNESCreate(PETSC_COMM_WORLD,snesDamage,ierr);CHKERRQ(ierr)
         Call SNESSetApplicationContext(snesDamage,MEF90DefMechCtx,ierr);CHKERRQ(ierr)
         Call SNESSetDM(snesDamage,MEF90DefMechCtx%DMScal,ierr);CHKERRQ(ierr)
         Call SNESSetOptionsPrefix(snesDamage,'damage_',ierr);CHKERRQ(ierr)
         !Call SNESSetType(snesDamage,SNESKSPONLY,ierr);CHKERRQ(ierr)
         
         !!! Set default bounds for the damage field
         Call DMCreateGlobalVector(MEF90DefMechCtx%DMScal,LB,ierr);CHKERRQ(ierr)
         Call VecDuplicate(LB,UB,ierr);CHKERRQ(ierr)
         Call VecSet(LB,0.0_Kr,ierr);CHKERRQ(ierr)
         Call VecSet(UB,1.0_Kr,ierr);CHKERRQ(ierr)
         Call SNESVISetVariableBounds(snesDamage,LB,UB,ierr);CHKERRQ(ierr)


         Call SNESSetFunction(snesDamage,residual,MEF90DefMechOperatorDamage,MEF90DefMechCtx,ierr);CHKERRQ(ierr)
         Call SNESSetJacobian(snesDamage,matDamage,matDamage,MEF90DefMechBilinearFormDamage,MEF90DefMechCtx,ierr);CHKERRQ(ierr)
         atol = 1.0D-8
         !rtol = 1.0D-10
         Call SNESSetTolerances(snesDamage,atol,PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_INTEGER,PETSC_DEFAULT_INTEGER,ierr);CHKERRQ(ierr)
         Call SNESSetFromOptions(snesDamage,ierr);CHKERRQ(ierr)

         !!! 
         !!! Set some KSP options
         !!!
         Call SNESGetKSP(snesDamage,kspDamage,ierr);CHKERRQ(ierr)
         Call KSPSetType(kspDamage,KSPCG,ierr);CHKERRQ(ierr)
         Call KSPSetInitialGuessNonzero(kspDamage,PETSC_TRUE,ierr);CHKERRQ(ierr)
         rtol = 1.0D-8
         atol = 1.0D-8
         dtol = 1.0D+10
         Call KSPSetTolerances(kspDamage,rtol,atol,dtol,PETSC_DEFAULT_INTEGER,ierr);CHKERRQ(ierr)
         Call KSPSetFromOptions(kspDamage,ierr);CHKERRQ(ierr)
      End If
      

      !!! set coordinates in PC for GAMG
      !!! For some reason, this makes gamg convergence worse, when the null space is specified.
      !!! Will investigate later
      !Call KSPGetPC(kspDamage,pcDamage,ierr);CHKERRQ(ierr)
      !Call DMMeshGetCoordinatesF90(MEF90DefMechCtx%DMScal,coordPtr,ierr);CHKERRQ(ierr)
      !Allocate(coordPCPtr(size(CoordPtr)))
      !coordPCPtr = reshape(transpose(coordPtr),[size(CoordPtr)])
      !coordPCPtr = reshape((coordPtr),[size(CoordPtr)])
      !Call PCSetCoordinates(pcDamage,dim,size(coordPtr),coordPCPtr,ierr);CHKERRQ(ierr)
      !DeAllocate(coordPCPtr)
      !Call DMMeshRestoreCoordinatesF90(MEF90DefMechCtx%DMScal,coordPtr,ierr);CHKERRQ(ierr)
      !Call PCSetFromOptions(pcDamage,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechCreateSolversDamage
   
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechUpdateDamageBounds"
!!!
!!!  
!!!  MEF90DefMechUpdateDamageBounds:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechUpdateDamageBounds(MEF90DefMechCtx,snesDamage,alpha,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      Type(SNES),Intent(OUT)                             :: snesDamage
      Type(Vec),Intent(IN)                               :: alpha
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      Type(Vec)                                          :: LB,UB
      PetscReal,Dimension(:),Pointer                     :: LBPtr
      PetscInt                                           :: i
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions

      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)      
      Call DMGetGlobalVector(MEF90DefMechCtx%DMScal,LB,ierr);CHKERRQ(ierr)
      Call DMGetGlobalVector(MEF90DefMechCtx%DMScal,UB,ierr);CHKERRQ(ierr)

      Call VecSet(UB,1.0_Kr,ierr);CHKERRQ(ierr)
      Call VecCopy(alpha,LB,ierr);CHKERRQ(ierr)
      If (MEF90DefMechGlobalOptions%irrevthres > 0.0_Kr) Then
         Call VecGetArrayF90(LB,LBPtr,ierr);CHKERRQ(ierr)
         Do i = 1, size(LBPtr)
            If (LBPtr(i) <= MEF90DefMechGlobalOptions%irrevthres) Then
               LBPtr(i) = 0.0_Kr
            End If
         End Do
         Call VecRestoreArrayF90(LB,LBPtr,ierr);CHKERRQ(ierr)
      End If
      Call SNESVISetVariableBounds(snesDamage,LB,UB,ierr);CHKERRQ(ierr)
      Call DMRestoreGlobalVector(MEF90DefMechCtx%DMScal,LB,ierr);CHKERRQ(ierr)
      Call DMRestoreGlobalVector(MEF90DefMechCtx%DMScal,UB,ierr);CHKERRQ(ierr)      
   End Subroutine MEF90DefMechUpdateDamageBounds
End Module m_MEF90_DefMech
