#include "../MEF90/mef90.inc"
Module m_MEF90_ElasticityCtx_Type
#include "finclude/petscdef.h"
   Use m_MEF90
   Use,Intrinsic :: iso_c_binding
   Implicit none
   !Private  
   !Public :: MEF90ElasticityCtx_Type
   !Public :: MEF90ElasticityGlobalOptions_Type
   !Public :: MEF90ElasticityCellSetOptions_Type
   !Public :: MEF90ElasticityVertexSetOptions_Type
   
   Type MEF90ElasticityCtx_Type
      PetscReal                        :: time

      Type(Vec),pointer                :: force
      Type(Vec),pointer                :: boundaryDisplacement
      Type(Vec),pointer                :: pressureForce
      Type(Vec),pointer                :: temperature      
      
      PetscBag                         :: GlobalOptionsBag
      PetscBag,Dimension(:),Pointer    :: CellSetOptionsBag
      PetscBag,Dimension(:),Pointer    :: VertexSetOptionsBag
      PetscBag,Dimension(:),Pointer    :: MaterialPropertiesBag
      Type(MEF90Ctx_Type),pointer      :: MEF90Ctx
      Type(DM),pointer                 :: DM
      Type(DM)                         :: cellDM
   End Type MEF90ElasticityCtx_Type
   
   Type MEF90ElasticityGlobalOptions_Type
      PetscInt                         :: mode
      PetscBool                        :: addNullSpace
      PetscInt                         :: DisplacementOffset
      PetscInt                         :: boundaryDisplacementScaling
      PetscInt                         :: boundaryDisplacementOffset
      PetscInt                         :: ForceScaling
      PetscInt                         :: ForceOffset
   End Type MEF90ElasticityGlobalOptions_Type

   Type MEF90ElasticityCellSetOptions_Type
      PetscInt                         :: elemTypeShortID
      PetscBool,Dimension(3)           :: Has_displacementBC
      PetscReal,Dimension(3)           :: boundaryDisplacement
      PetscReal,Dimension(3)           :: Force
      PetscReal                        :: pressureForce
   End Type MEF90ElasticityCellSetOptions_Type

   Type MEF90ElasticityVertexSetOptions_Type
      PetscBool,Dimension(3)           :: Has_displacementBC
      PetscReal,Dimension(3)           :: boundaryDisplacement
   End Type MEF90ElasticityVertexSetOptions_Type 
End Module m_MEF90_ElasticityCtx_Type

Module m_MEF90ElasticityGlobalOptions_Private
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_ElasticityCtx_Type
   Implicit None

   Private
   Public :: PetscBagGetDataMEF90ElasticityCtxGlobalOptions
   
   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF90_ElasticityCtx_Type
         PetscBag                                           :: bag
         Type(MEF90ElasticityGlobalOptions_Type),pointer    :: data
         PetscErrorCode                                     :: ierr
      End subroutine PetscBagGetData
   End interface
Contains

#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90ElasticityCtxGlobalOptions"
!!!
!!!  PetscBagGetDataMEF90ElasticityCtxGlobalOptions - Custom interface to PetscGetData
!!!
   Subroutine PetscBagGetDataMEF90ElasticityCtxGlobalOptions(bag,data,ierr)
      PetscBag                                        :: bag
      Type(MEF90ElasticityGlobalOptions_Type),pointer :: data
      PetscErrorCode                                  :: ierr
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataMEF90ElasticityCtxGlobalOptions
End Module m_MEF90ElasticityGlobalOptions_Private

Module m_MEF90ElasticityCellSetOptions_Private
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_ElasticityCtx_Type
   Implicit None

   Private
   Public :: PetscBagGetDataMEF90ElasticityCtxCellSetOptions
   
   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF90_ElasticityCtx_Type
         PetscBag                                           :: bag
         Type(MEF90ElasticityCellSetOptions_Type),pointer   :: data
         PetscErrorCode                                     :: ierr
      End subroutine PetscBagGetData
   End interface
Contains

#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90ElasticityCtxCellSetOptions"
!!!
!!!  PetscBagGetDataMEF90ElasticityCtxCellSetOptions - Custom interface to PetscGetData
!!!
   Subroutine PetscBagGetDataMEF90ElasticityCtxCellSetOptions(bag,data,ierr)
      PetscBag                                          :: bag
      Type(MEF90ElasticityCellSetOptions_Type),pointer  :: data
      PetscErrorCode                                    :: ierr
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataMEF90ElasticityCtxCellSetOptions
End Module m_MEF90ElasticityCellSetOptions_Private

Module m_MEF90ElasticityVertexSetOptions_Private
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_ElasticityCtx_Type
   Implicit None

   Private
   Public :: PetscBagGetDataMEF90ElasticityCtxVertexSetOptions
   
   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF90_ElasticityCtx_Type
         PetscBag                                           :: bag
         Type(MEF90ElasticityVertexSetOptions_Type),pointer :: data
         PetscErrorCode                                     :: ierr
      End subroutine PetscBagGetData
   End interface

Contains
#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90ElasticityCtxVertexSetOptions"
!!!
!!!  PetscBagGetDataMEF90ElasticityCtxVertexSetOptions - Custom interface to PetscGetData
!!!
   Subroutine PetscBagGetDataMEF90ElasticityCtxVertexSetOptions(bag,data,ierr)
      PetscBag                                           :: bag
      Type(MEF90ElasticityVertexSetOptions_Type),pointer :: data
      PetscErrorCode                                     :: ierr
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataMEF90ElasticityCtxVertexSetOptions
End Module m_MEF90ElasticityVertexSetOptions_Private

Module m_MEF90_ElasticityCtx
#include "finclude/petscdef.h"
   Use m_MEF90
   Use m_MEF90_ElasticityCtx_Type
   Use m_MEF90ElasticityGlobalOptions_Private
   Use m_MEF90ElasticityCellSetOptions_Private
   Use m_MEF90ElasticityVertexSetOptions_Private
   Implicit none

   PetscSizeT,protected   :: sizeofMEF90ElasticityGlobalOptions
   PetscSizeT,protected   :: sizeofMEF90ElasticityCellSetOptions
   PetscSizeT,protected   :: sizeofMEF90ElasticityVertexSetOptions
   
   Enum,bind(c)
      enumerator  :: MEF90Elasticity_ModeSteadyState = 0, &
                     MEF90Elasticity_ModeSteadyTransient
   End Enum
   Character(len = MEF90_MXSTRLEN),dimension(5),protected   :: MEF90Elasticity_ModeList
   
Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90ElasticityCtx_InitializePrivate"
!!!
!!!  
!!!  MEF90ElasticityCtx_InitializePrivate:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90ElasticityCtx_InitializePrivate(ierr)
      PetscErrorCode,Intent(OUT)                      :: ierr
   
      Type(MEF90ElasticityGlobalOptions_Type)         :: ElasticityGlobalOptions
      Type(MEF90ElasticityCellSetOptions_Type)        :: ElasticityCellSetOptions
      Type(MEF90ElasticityVertexSetOptions_Type)      :: ElasticityVertexSetOptions
      character(len=1),pointer                        :: dummychar(:)
      PetscSizeT                                      :: sizeofchar
      
      Call PetscDataTypeGetSize(PETSC_CHAR,sizeofchar,ierr)
      sizeofMEF90ElasticityGlobalOptions = size(transfer(ElasticityGlobalOptions,dummychar))*sizeofchar
      sizeofMEF90ElasticityCellSetOptions = size(transfer(ElasticityCellSetOptions,dummychar))*sizeofchar
      sizeofMEF90ElasticityVertexSetOptions = size(transfer(ElasticityVertexSetOptions,dummychar))*sizeofchar

      MEF90Elasticity_ModeList(1) = 'SteadyState'
      MEF90Elasticity_ModeList(2) = 'Transient'
      MEF90Elasticity_ModeList(3) = 'MEF90_Elasticity_Mode'
      MEF90Elasticity_ModeList(4) = '_MEF90_Elasticity_Mode'
      MEF90Elasticity_ModeList(5) = ''
   End Subroutine MEF90ElasticityCtx_InitializePrivate
   
#undef __FUNCT__
#define __FUNCT__ "MEF90ElasticityCtx_Create"
!!!
!!!  
!!!  MEF90ElasticityCtx_Create:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90ElasticityCtx_Create(ElasticityCtx,Mesh,MEF90Ctx,ierr)
      Type(MEF90ElasticityCtx_Type),Intent(OUT)          :: ElasticityCtx
      Type(DM),target,Intent(IN)                         :: Mesh
      Type(MEF90Ctx_Type),target,Intent(IN)              :: MEF90Ctx
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90CtxGlobalOptions
      Type(MEF90ElasticityGlobalOptions_Type),pointer    :: MEF90ElasticityGlobalOptions
      Type(MEF90ElasticityCellSetOptions_Type),pointer   :: MEF90ElasticityCellSetOptions
      Type(MEF90ElasticityVertexSetOptions_Type),pointer :: MEF90ElasticityVertexSetOptions
      Type(IS)                                           :: setIS
      PetscInt                                           :: set,numSet

      Call MEF90ElasticityCtx_InitializePrivate(ierr)
      ElasticityCtx%DM => Mesh
      ElasticityCtx%MEF90Ctx => MEF90Ctx
      Call DMMeshClone(ElasticityCtx%DM,ElasticityCtx%cellDM,ierr);CHKERRQ(ierr)
      Call DMMeshSetMaxDof(ElasticityCtx%cellDM,1,ierr);CHKERRQ(ierr) 

      Call PetscBagCreate(MEF90Ctx%comm,sizeofMEF90ElasticityGlobalOptions,ElasticityCtx%GlobalOptionsBag,ierr);CHKERRQ(ierr)
      
      !!! Call DMmeshGetLabelSize(Mesh,'Cell Sets',numSet,ierr);CHKERRQ(ierr)
      !!! NO: I need to allocate for the overall number of sets, not the local one

      Call DMmeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr) 
      Call ISGetSize(setIS,numSet,ierr);CHKERRQ(ierr)

      Allocate(ElasticityCtx%CellSetOptionsBag(numSet),stat=ierr)
      Do set = 1, numSet
         Call PetscBagCreate(MEF90Ctx%comm,sizeofMEF90ElasticityCellSetOptions,ElasticityCtx%CellSetOptionsBag(set),ierr);CHKERRQ(ierr)
      End Do
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(mesh,'Vertex Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr) 
      Call ISGetSize(setIS,numSet,ierr);CHKERRQ(ierr)
      Allocate(ElasticityCtx%VertexSetOptionsBag(numSet),stat=ierr)
      Do set = 1, numSet
         Call PetscBagCreate(MEF90Ctx%comm,sizeofMEF90ElasticityVertexSetOptions,ElasticityCtx%VertexSetOptionsBag(set),ierr);CHKERRQ(ierr)
      End Do
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      
      Nullify(ElasticityCtx%ForcePrevious)
      Nullify(ElasticityCtx%ForceTarget)
      Nullify(ElasticityCtx%boundaryDisplacementPrevious)
      Nullify(ElasticityCtx%boundaryDisplacementTarget)
   End Subroutine MEF90ElasticityCtx_Create
   
#undef __FUNCT__
#define __FUNCT__ "MEF90ElasticityCtx_Destroy"
!!!
!!!  
!!!  MEF90ElasticityCtx_Destroy:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90ElasticityCtx_Destroy(ElasticityCtx,ierr)
      Type(MEF90ElasticityCtx_Type),Intent(OUT)       :: ElasticityCtx
      PetscErrorCode,Intent(OUT)                      :: ierr
      
      PetscInt                                        :: set
   
      Call PetscBagDestroy(ElasticityCtx%GlobalOptionsBag,ierr);CHKERRQ(ierr)
      Do set = 1, size(ElasticityCtx%CellSetOptionsBag)
         Call PetscBagDestroy(ElasticityCtx%CellSetOptionsBag(set),ierr);CHKERRQ(ierr)
      End Do
      DeAllocate(ElasticityCtx%CellSetOptionsBag)
      Do set = 1, size(ElasticityCtx%VertexSetOptionsBag)
         Call PetscBagDestroy(ElasticityCtx%VertexSetOptionsBag(set),ierr);CHKERRQ(ierr)
      End Do
      DeAllocate(ElasticityCtx%VertexSetOptionsBag)
      
      Nullify(ElasticityCtx%DM)
      Nullify(ElasticityCtx%MEF90Ctx)
      Nullify(ElasticityCtx%ForcePrevious)
      Nullify(ElasticityCtx%ForceTarget)
      Nullify(ElasticityCtx%boundaryDisplacementPrevious)
      Nullify(ElasticityCtx%boundaryDisplacementTarget)
      Nullify(ElasticityCtx%DM)
      Call DMDestroy(ElasticityCtx%cellDM,ierr);CHKERRQ(ierr)
   End Subroutine MEF90ElasticityCtx_Destroy

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90ElasticityCtxGlobalOptions"
!!!
!!!  
!!!  PetscBagRegisterMEF90ElasticityCtxGlobalOptions:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine PetscBagRegisterMEF90ElasticityCtxGlobalOptions(bag,name,prefix,default,ierr)
      PetscBag                                           :: bag
      Character(len=*),Intent(IN)                        :: prefix,name
      Type(MEF90ElasticityGlobalOptions_Type),Intent(IN) :: default
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(MEF90ElasticityGlobalOptions_Type),pointer      :: ElasticityGlobalOptions
      Call PetscBagGetDataMEF90ElasticityCtxGlobalOptions(bag,ElasticityGlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagSetName(bag,trim(name),"ElasticityGlobalOptions MEF90 Heat transfer global options",ierr);CHKERRQ(ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix),ierr);CHKERRQ(ierr)

      Call PetscBagRegisterEnum(bag,ElasticityGlobalOptions%mode,MEF90Elasticity_ModeList,default%mode,'Elasticity_mode','Type of heat transfer computation',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterBool(bag,ElasticityGlobalOptions%addNullSpace,default%addNullSpace,'addNullSpace','Add null space to SNES',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,ElasticityGlobalOptions%DisplacementOffset,default%DisplacementOffset,'Displacement_Offset','Position of Displacement field in EXO file',ierr);CHKERRQ(ierr)

      Call PetscBagRegisterEnum(bag,ElasticityGlobalOptions%boundaryDisplacementScaling,MEF90ScalingList,default%boundaryDisplacementScaling,'boundaryDisplacement_scaling','Boundary Displacement scaling',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,ElasticityGlobalOptions%boundaryDisplacementOffset,default%boundaryDisplacementOffset,'boundaryDisplacement_Offset','Position of boundary Displacement field in EXO file',ierr);CHKERRQ(ierr)

      Call PetscBagRegisterEnum(bag,ElasticityGlobalOptions%ForceScaling,MEF90ScalingList,default%ForceScaling,'Force_scaling','Heat Force scaling',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt(bag,ElasticityGlobalOptions%ForceOffset,default%ForceOffset,'Force_Offset','Position of heat Force field in EXO file',ierr);CHKERRQ(ierr)      
   End Subroutine PetscBagRegisterMEF90ElasticityCtxGlobalOptions

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90ElasticityCtxCellSetOptions"
!!!
!!!  
!!!  PetscBagRegisterMEF90ElasticityCtxCellSetOptions:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine PetscBagRegisterMEF90ElasticityCtxCellSetOptions(bag,name,prefix,default,ierr)
      PetscBag                                             :: bag
      Character(len=*),Intent(IN)                          :: prefix,name
      Type(MEF90ElasticityCellSetOptions_Type),Intent(IN)  :: default
      PetscErrorCode,Intent(OUT)                           :: ierr

      Type(MEF90ElasticityCellSetOptions_Type),pointer      :: ElasticityCellSetOptions
      Call PetscBagGetDataMEF90ElasticityCtxCellSetOptions(bag,ElasticityCellSetOptions,ierr);CHKERRQ(ierr)
      Call PetscBagSetName(bag,trim(name),"ElasticityCellSetOptions MEF90 Heat transfer Cell Set options",ierr);CHKERRQ(ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix),ierr);CHKERRQ(ierr)

      Call PetscBagRegisterInt(bag,ElasticityCellSetOptions%ElemTypeShortID,default%ElemTypeShortID,'ShortID','Element type ShortID',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterRealArray(bag,ElasticityCellSetOptions%force,3,'Force','[N.m^(-3) / N.m^(-2) / N.m^(-1)] (f): body / boundary force',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,ElasticityCellSetOptions%pressureForce,default%pressureForce,'pressureForce','[N.m^(-2) / N.m^(-1)] (p): boundary pressureforce',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterBoolArray(bag,ElasticityCellSetOptions%Has_displacementBC,3,'DisplacementBC','Displacement has Dirichlet boundary Condition (Y/N)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterRealArray(bag,ElasticityCellSetOptions%boundaryDisplacement,3,'boundaryDisplacement','Displacement boundary value',ierr);CHKERRQ(ierr)
   End Subroutine PetscBagRegisterMEF90ElasticityCtxCellSetOptions

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90ElasticityCtxVertexSetOptions"
!!!
!!!  
!!!  PetscBagRegisterMEF90ElasticityCtxVertexSetOptions:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine PetscBagRegisterMEF90ElasticityCtxVertexSetOptions(bag,name,prefix,default,ierr)
      PetscBag                                              :: bag
      Character(len=*),Intent(IN)                           :: prefix,name
      Type(MEF90ElasticityVertexSetOptions_Type),Intent(IN) :: default
      PetscErrorCode,Intent(OUT)                            :: ierr

      Type(MEF90ElasticityVertexSetOptions_Type),pointer      :: ElasticityVertexSetOptions
      Call PetscBagGetDataMEF90ElasticityCtxVertexSetOptions(bag,ElasticityVertexSetOptions,ierr);CHKERRQ(ierr)
      Call PetscBagSetName(bag,trim(name),"ElasticityVertexSetOptions MEF90 Heat transfer Vertex Set options",ierr);CHKERRQ(ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix),ierr);CHKERRQ(ierr)

      Call PetscBagRegisterBoolArray(bag,ElasticityVertexSetOptions%Has_displacementBC,3,'DisplacementBC','Displacement has Dirichlet boundary Condition (Y/N)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterRealArray(bag,ElasticityVertexSetOptions%boundaryDisplacement,3,'boundaryDisplacement','Displacement boundary value',ierr);CHKERRQ(ierr)
   End Subroutine PetscBagRegisterMEF90ElasticityCtxVertexSetOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90ElasticityCtx_SetFromOptions"
!!!
!!!  
!!!  MEF90ElasticityCtx_SetFromOptions:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90ElasticityCtx_SetFromOptions(ElasticityCtx,prefix,defaultGlobalOptions, &
                                              defaultCellSetOptions,defaultVertexSetOptions,ierr)
      Type(MEF90ElasticityCtx_Type),Intent(OUT)               :: ElasticityCtx
      Character(len=*),Intent(IN)                             :: prefix
      Type(MEF90ElasticityGlobalOptions_Type),Intent(IN)      :: defaultGlobalOptions
      Type(MEF90ElasticityCellSetOptions_Type),Intent(IN)     :: defaultCellSetOptions
      Type(MEF90ElasticityVertexSetOptions_Type),Intent(IN)   :: defaultVertexSetOptions
      PetscErrorCode,Intent(OUT)                              :: ierr
   
      Type(MEF90CtxGlobalOptions_Type),pointer                :: MEF90CtxGlobalOptions
      Type(MEF90ElasticityCellSetOptions_Type)                :: myDefaultCellSetOptions
      Type(MEF90Element_Type),Dimension(:),Pointer            :: ElemType
      Type(IS)                                                :: setIS
      PetscInt,Dimension(:),Pointer                           :: setID
      PetscInt                                                :: set
      Character(len=MEF90_MXSTRLEN)                           :: IOBuffer,setName,setprefix

      Call PetscBagGetDataMEF90CtxGlobalOptions(ElasticityCtx%MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr);CHKERRQ(ierr)
      !!!
      !!! Registering Global Context
      !!!
      Call PetscBagRegisterMEF90ElasticityCtxGlobalOptions(ElasticityCtx%GlobalOptionsBag,"MEF90Elasticity Global Ctx",prefix,defaultGlobalOptions,ierr);CHKERRQ(ierr)

      If (MEF90CtxGlobalOptions%verbose > 0) Then
         Call PetscBagView(ElasticityCtx%GlobalOptionsBag,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
      End If

      !!!
      !!! Registering Cell Set Context
      !!! We override the default element type with that detected from the exodus file
      !!!
      Call DMmeshGetLabelIdIS(ElasticityCtx%DM,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(ElasticityCtx%MEF90Ctx%comm,setIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      
      Call EXOGetCellSetElementType_Scal(ElasticityCtx%MEF90Ctx,ElemType,ierr)
      Do set = 1, size(setID)
         Write(setName,100) setID(set)
         Write(setprefix,101) setID(set)
         mydefaultCellSetOptions = defaultCellSetOptions
         mydefaultCellSetOptions%ElemTypeShortID = ElemType(set)%ShortID
         Call PetscBagRegisterMEF90ElasticityCtxCellSetOptions(ElasticityCtx%CellSetOptionsBag(set),setName,setPrefix,mydefaultCellSetOptions,ierr)
         If (MEF90CtxGlobalOptions%verbose > 0) Then
            Write(IOBuffer,103) setID(set),trim(setprefix)
            Call PetscPrintf(ElasticityCtx%MEF90Ctx%comm,IOBuffer,ierr);CHKERRQ(ierr)
            Call PetscBagView(ElasticityCtx%CellSetOptionsBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
            Call PetscPrintf(ElasticityCtx%MEF90Ctx%comm,"\n",ierr);CHKERRQ(ierr)
         End if
      End Do
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      DeAllocate(ElemType)
      
      !!!
      !!! Registering Vertex Set Context
      !!! We override the default element type with that detected from the exodus file
      !!!
      Call DMmeshGetLabelIdIS(ElasticityCtx%DM,'Vertex Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      
      Do set = 1, size(setID)
         Write(setName,200) setID(set)
         Write(setprefix,201) setID(set)
         Call PetscBagRegisterMEF90ElasticityCtxVertexSetOptions(ElasticityCtx%VertexSetOptionsBag(set),setName,setPrefix,defaultVertexSetOptions,ierr)
         If (MEF90CtxGlobalOptions%verbose > 0) Then
            Write(IOBuffer,203) setID(set),trim(setprefix)
            Call PetscPrintf(ElasticityCtx%MEF90Ctx%comm,IOBuffer,ierr);CHKERRQ(ierr)
            Call PetscBagView(ElasticityCtx%VertexSetOptionsBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
            Call PetscPrintf(ElasticityCtx%MEF90Ctx%comm,"\n",ierr);CHKERRQ(ierr)
         End if
      End Do
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)

100 Format('Cell set ',I4)
101 Format('cs',I4.4,'_')
103 Format('Registering cell set ',I4,' prefix: ',A,'\n')
200 Format('Vertex set ',I4)
201 Format('vs',I4.4,'_')
203 Format('Registering vertex set ',I4,' prefix: ',A,'\n')
   End Subroutine MEF90ElasticityCtx_SetFromOptions

End Module m_MEF90_ElasticityCtx
