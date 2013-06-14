#include "../MEF90/mef90.inc"
Module m_MEF90_DefMEch
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_DefMechCtx

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSetTransients"
!!!
!!!  
!!!  MEF90DefMechSetTransients:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
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
            Call MEF90DefMechSetforceCst(MEF90DefMechCtx%force,MEF90DefMechCtx,ierr)
            Call VecScale(MEF90DefMechCtx%force,time,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_CST)
            Call MEF90DefMechSetforceCst(MEF90DefMechCtx%force,MEF90DefMechCtx,ierr)
      End Select
      Select case (MEF90DefMechGlobalOptions%pressureForceScaling)
         Case (MEF90Scaling_File)
            Call VecLoadExodusCell(MEF90DefMechCtx%cellDMScal,MEF90DefMechCtx%pressureForce,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
                                   MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%pressureForceOffset,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_Linear)
            !Call MEF90DefMechSetforceCst(MEF90DefMechCtx%pressureForce,MEF90DefMechCtx,ierr)
            Call VecScale(MEF90DefMechCtx%pressureForce,time,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_CST)
            Call MEF90DefMechSetforceCst(MEF90DefMechCtx%pressureForce,MEF90DefMechCtx,ierr)
      End Select
      Select case (MEF90DefMechGlobalOptions%boundaryDisplacementScaling)
         Case (MEF90Scaling_File)
            Call VecLoadExodusVertex(MEF90DefMechCtx%cellDMVect,MEF90DefMechCtx%boundaryDisplacement,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
                                     MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%boundaryDisplacementOffset,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_Linear)
            !Call MEF90DefMechSetboundaryDisplacementCst(MEF90DefMechCtx%boundaryDisplacement,MEF90DefMechCtx,ierr)
            Call VecScale(MEF90DefMechCtx%boundaryDisplacement,time,ierr);CHKERRQ(ierr)
         Case (MEF90Scaling_CST)
            !Call MEF90DefMechSetboundaryDisplacementCst(MEF90DefMechCtx%boundaryDisplacement,MEF90DefMechCtx,ierr)
      End Select
   End Subroutine MEF90DefMechSetTransients

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSetForceCst"
!!!
!!!  
!!!  MEF90DefMechSetForceCst:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechSetForceCst(x,MEF90DefMechCtx,ierr)
      Type(Vec),Intent(IN)                            :: x
      Type(MEF90DefMechCtx_Type),Intent(IN)           :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                      :: ierr
   
      Type(MEF90DefMechGlobalOptions_Type),pointer    :: MEF90DefMechGlobalOptions
      Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90GlobalOptions
      Type(MEF90DefMechCellSetOptions_Type),pointer   :: cellSetOptions
      Type(IS)                                        :: cellSetGlobalIS,setIS,setISdof
      PetscInt,dimension(:),Pointer                   :: setID
      PetscInt,Dimension(:),Pointer                   :: setIdx
      PetscInt                                        :: set
      PetscReal,Dimension(:),Pointer                  :: val
      
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)
      
      !!! force is cell-centered
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMScal,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%cellDMScal,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call DMMeshISCreateISglobaldof(MEF90DefMechCtx%cellDMScal,setIS,0,setISdof,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
         Allocate(val(size(setIdx)),stat=ierr)
         val = cellSetOptions%force
         Call VecSetValues(x,size(setIdx),setIdx,val,INSERT_VALUES,ierr);CHKERRQ(ierr)
         DeAllocate(val)
         Call ISRestoreIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
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
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechSetPressureForceCst(x,MEF90DefMechCtx,ierr)
      Type(Vec),Intent(IN)                            :: x
      Type(MEF90DefMechCtx_Type),Intent(IN)           :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                      :: ierr
   
      Type(MEF90DefMechGlobalOptions_Type),pointer    :: MEF90DefMechGlobalOptions
      Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90GlobalOptions
      Type(MEF90DefMechCellSetOptions_Type),pointer   :: cellSetOptions
      Type(IS)                                        :: cellSetGlobalIS,setIS,setISdof
      PetscInt,dimension(:),Pointer                   :: setID
      PetscInt,Dimension(:),Pointer                   :: setIdx
      PetscInt                                        :: set
      PetscReal,Dimension(:),Pointer                  :: val
      
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)
      
      !!! pressureForce is cell-centered
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMScal,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
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
      End Do
      Call ISRestoreIndicesF90(cellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(cellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call VecAssemblyBegin(x,ierr);CHKERRQ(ierr)
      Call VecAssemblyEnd(x,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechSetPressureForceCst


End Module m_MEF90_DefMech
