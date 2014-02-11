#include "../MEF90/mef90.inc"
Module m_MEF90_HeatXferCtx_Type
#include "finclude/petscdef.h"
   Use m_MEF90
   Use,Intrinsic :: iso_c_binding
   Implicit none
   Private  
   Public :: MEF90HeatXferCtx_Type
   Public :: MEF90HeatXferGlobalOptions_Type
   Public :: MEF90HeatXferCellSetOptions_Type
   Public :: MEF90HeatXferVertexSetOptions_Type
   
   Type MEF90HeatXferCtx_Type
      Type(Vec),pointer                :: temperature
      Type(Vec),pointer                :: flux
      Type(Vec),pointer                :: boundaryTemperature
      Type(Vec),pointer                :: externalTemperature

      PetscBag                         :: GlobalOptionsBag
      PetscBag,Dimension(:),Pointer    :: CellSetOptionsBag
      PetscBag,Dimension(:),Pointer    :: VertexSetOptionsBag
      PetscBag,Dimension(:),Pointer    :: MaterialPropertiesBag
      Type(MEF90Ctx_Type),pointer      :: MEF90Ctx
      Type(DM),pointer                 :: DM
      Type(DM)                         :: cellDM
   End Type MEF90HeatXferCtx_Type
   
   Type MEF90HeatXferGlobalOptions_Type
      PetscEnum                        :: mode
      PetscBool                        :: addNullSpace
      PetscInt                         :: tempOffset
      PetscReal                        :: initialTemperature
      PetscInt                         :: boundaryTempScaling
      PetscInt                         :: boundaryTempOffset
      PetscInt                         :: externalTempScaling
      PetscInt                         :: externalTempOffset
      PetscInt                         :: fluxScaling
      PetscInt                         :: fluxOffset
      !!! offset  = position in data file (required for exodus)
      !!! scaling = time (step) scaling law currently CST, Linear, Null (not present) File
   End Type MEF90HeatXferGlobalOptions_Type

   Type MEF90HeatXferCellSetOptions_Type
      PetscInt                         :: elemTypeShortID
      PetscReal                        :: flux
      PetscReal                        :: surfaceThermalConductivity
      PetscReal                        :: externalTemp
      PetscBool                        :: Has_BC
      PetscReal                        :: boundaryTemp
   End Type MEF90HeatXferCellSetOptions_Type

   Type MEF90HeatXferVertexSetOptions_Type
      PetscBool                        :: Has_BC
      PetscReal                        :: boundaryTemp
   End Type MEF90HeatXferVertexSetOptions_Type 
End Module m_MEF90_HeatXferCtx_Type

Module m_MEF90HeatXferGlobalOptions_Private
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_HeatXferCtx_Type
   Implicit None

   Private
   Public :: PetscBagGetDataMEF90HeatXferCtxGlobalOptions
   
   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF90_HeatXferCtx_Type
         PetscBag                                           :: bag
         Type(MEF90HeatXferGlobalOptions_Type),pointer      :: data
         PetscErrorCode                                     :: ierr
      End subroutine PetscBagGetData
   End interface
Contains

#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90HeatXferCtxGlobalOptions"
!!!
!!!  PetscBagGetDataMEF90HeatXferCtxGlobalOptions - Custom interface to PetscGetData
!!!
   Subroutine PetscBagGetDataMEF90HeatXferCtxGlobalOptions(bag,data,ierr)
      PetscBag                                        :: bag
      Type(MEF90HeatXferGlobalOptions_Type),pointer   :: data
      PetscErrorCode                                  :: ierr
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataMEF90HeatXferCtxGlobalOptions
End Module m_MEF90HeatXferGlobalOptions_Private

Module m_MEF90HeatXferCellSetOptions_Private
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_HeatXferCtx_Type
   Implicit None

   Private
   Public :: PetscBagGetDataMEF90HeatXferCtxCellSetOptions
   
   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF90_HeatXferCtx_Type
         PetscBag                                           :: bag
         Type(MEF90HeatXferCellSetOptions_Type),pointer     :: data
         PetscErrorCode                                     :: ierr
      End subroutine PetscBagGetData
   End interface
Contains

#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90HeatXferCtxCellSetOptions"
!!!
!!!  PetscBagGetDataMEF90HeatXferCtxCellSetOptions - Custom interface to PetscGetData
!!!
   Subroutine PetscBagGetDataMEF90HeatXferCtxCellSetOptions(bag,data,ierr)
      PetscBag                                        :: bag
      Type(MEF90HeatXferCellSetOptions_Type),pointer  :: data
      PetscErrorCode                                  :: ierr
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataMEF90HeatXferCtxCellSetOptions
End Module m_MEF90HeatXferCellSetOptions_Private

Module m_MEF90HeatXferVertexSetOptions_Private
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_HeatXferCtx_Type
   Implicit None

   Private
   Public :: PetscBagGetDataMEF90HeatXferCtxVertexSetOptions
   
   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF90_HeatXferCtx_Type
         PetscBag                                           :: bag
         Type(MEF90HeatXferVertexSetOptions_Type),pointer   :: data
         PetscErrorCode                                     :: ierr
      End subroutine PetscBagGetData
   End interface

Contains
#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90HeatXferCtxVertexSetOptions"
!!!
!!!  PetscBagGetDataMEF90HeatXferCtxVertexSetOptions - Custom interface to PetscGetData
!!!
   Subroutine PetscBagGetDataMEF90HeatXferCtxVertexSetOptions(bag,data,ierr)
      PetscBag                                           :: bag
      Type(MEF90HeatXferVertexSetOptions_Type),pointer   :: data
      PetscErrorCode                                     :: ierr
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataMEF90HeatXferCtxVertexSetOptions
End Module m_MEF90HeatXferVertexSetOptions_Private

Module m_MEF90_HeatXferCtx
#include "finclude/petscdef.h"
   Use m_MEF90
   Use m_MEF90_HeatXferCtx_Type
   Use m_MEF90HeatXferGlobalOptions_Private
   Use m_MEF90HeatXferCellSetOptions_Private
   Use m_MEF90HeatXferVertexSetOptions_Private
   Implicit none

   PetscSizeT,protected   :: sizeofMEF90HeatXferGlobalOptions
   PetscSizeT,protected   :: sizeofMEF90HeatXferCellSetOptions
   PetscSizeT,protected   :: sizeofMEF90HeatXferVertexSetOptions
   
   Enum,bind(c)
      enumerator  :: MEF90HeatXfer_ModeNULL = 0,    &
                     MEF90HeatXFer_ModeSteadyState, &
                     MEF90HeatXFer_ModeTransient
   End Enum
   Character(len = MEF90_MXSTRLEN),dimension(6),protected   :: MEF90HeatXFer_ModeList
   
Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferCtxInitialize_Private"
!!!
!!!  
!!!  MEF90HeatXferCtxInitialize_Private:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferCtxInitialize_Private(ierr)
      PetscErrorCode,Intent(OUT)                      :: ierr
   
      Type(MEF90HeatXferGlobalOptions_Type)           :: HeatXferGlobalOptions
      Type(MEF90HeatXferCellSetOptions_Type)          :: HeatXferCellSetOptions
      Type(MEF90HeatXferVertexSetOptions_Type)        :: HeatXferVertexSetOptions
      character(len=1),pointer                        :: dummychar(:)
      PetscSizeT                                      :: sizeofchar
      
      Call PetscDataTypeGetSize(PETSC_CHAR,sizeofchar,ierr)
      sizeofMEF90HeatXferGlobalOptions = size(transfer(HeatXferGlobalOptions,dummychar))*sizeofchar
      sizeofMEF90HeatXferCellSetOptions = size(transfer(HeatXferCellSetOptions,dummychar))*sizeofchar
      sizeofMEF90HeatXferVertexSetOptions = size(transfer(HeatXferVertexSetOptions,dummychar))*sizeofchar

      MEF90HeatXFer_ModeList(1) = 'null'
      MEF90HeatXFer_ModeList(2) = 'SteadyState'
      MEF90HeatXFer_ModeList(3) = 'Transient'
      MEF90HeatXFer_ModeList(4) = 'MEF90_HeatXFer_Mode'
      MEF90HeatXFer_ModeList(5) = '_MEF90_HeatXFer_Mode'
      MEF90HeatXFer_ModeList(6) = ''
   End Subroutine MEF90HeatXferCtxInitialize_Private
   
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferCtxCreate"
!!!
!!!  
!!!  MEF90HeatXferCtxCreate:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferCtxCreate(HeatXferCtx,Mesh,MEF90Ctx,ierr)
      Type(MEF90HeatXferCtx_Type),Intent(OUT)            :: HeatXferCtx
      Type(DM),target,Intent(IN)                         :: Mesh
      Type(MEF90Ctx_Type),target,Intent(IN)              :: MEF90Ctx
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90CtxGlobalOptions
      Type(MEF90HeatXferGlobalOptions_Type),pointer      :: MEF90HeatXferGlobalOptions
      Type(MEF90HeatXferCellSetOptions_Type),pointer     :: MEF90HeatXferCellSetOptions
      Type(MEF90HeatXferVertexSetOptions_Type),pointer   :: MEF90HeatXferVertexSetOptions
      Type(IS)                                           :: setIS
      PetscInt                                           :: set,numSet

      Call MEF90HeatXferCtxInitialize_Private(ierr)
      HeatXferCtx%DM => Mesh
      HeatXferCtx%MEF90Ctx => MEF90Ctx
      Call DMMeshClone(HeatXferCtx%DM,HeatXferCtx%cellDM,ierr);CHKERRQ(ierr)
      Call DMMeshSetMaxDof(HeatXferCtx%cellDM,1,ierr);CHKERRQ(ierr) 

      Call PetscBagCreate(MEF90Ctx%comm,sizeofMEF90HeatXferGlobalOptions,HeatXferCtx%GlobalOptionsBag,ierr);CHKERRQ(ierr)
      
      !!! Call DMmeshGetLabelSize(Mesh,'Cell Sets',numSet,ierr);CHKERRQ(ierr)
      !!! NO: I need to allocate for the overall number of sets, not the local one

      Call DMmeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr) 
      Call ISGetLocalSize(setIS,numSet,ierr);CHKERRQ(ierr)
      Allocate(HeatXferCtx%CellSetOptionsBag(numSet),stat=ierr)
      Do set = 1, numSet
         Call PetscBagCreate(MEF90Ctx%comm,sizeofMEF90HeatXferCellSetOptions,HeatXferCtx%CellSetOptionsBag(set),ierr);CHKERRQ(ierr)
      End Do
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(mesh,'Vertex Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr) 
      Call ISGetLocalSize(setIS,numSet,ierr);CHKERRQ(ierr)
      Allocate(HeatXferCtx%VertexSetOptionsBag(numSet),stat=ierr)
      Do set = 1, numSet
         Call PetscBagCreate(MEF90Ctx%comm,sizeofMEF90HeatXferVertexSetOptions,HeatXferCtx%VertexSetOptionsBag(set),ierr);CHKERRQ(ierr)
      End Do
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      
      Nullify(HeatXferCtx%flux)
      Nullify(HeatXferCtx%boundaryTemperature)
      Nullify(HeatXferCtx%externalTemperature)
   End Subroutine MEF90HeatXferCtxCreate
   
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferCtxDestroy"
!!!
!!!  
!!!  MEF90HeatXferCtxDestroy:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferCtxDestroy(HeatXferCtx,ierr)
      Type(MEF90HeatXferCtx_Type),Intent(OUT)         :: HeatXferCtx
      PetscErrorCode,Intent(OUT)                      :: ierr
      
      PetscInt                                        :: set
   
      Call PetscBagDestroy(HeatXferCtx%GlobalOptionsBag,ierr);CHKERRQ(ierr)
      Do set = 1, size(HeatXferCtx%CellSetOptionsBag)
         Call PetscBagDestroy(HeatXferCtx%CellSetOptionsBag(set),ierr);CHKERRQ(ierr)
      End Do
      DeAllocate(HeatXferCtx%CellSetOptionsBag)
      Do set = 1, size(HeatXferCtx%VertexSetOptionsBag)
         Call PetscBagDestroy(HeatXferCtx%VertexSetOptionsBag(set),ierr);CHKERRQ(ierr)
      End Do
      DeAllocate(HeatXferCtx%VertexSetOptionsBag)
      
      Nullify(HeatXferCtx%DM)
      Nullify(HeatXferCtx%MEF90Ctx)
      Nullify(HeatXferCtx%flux)
      Nullify(HeatXferCtx%boundaryTemperature)
      Nullify(HeatXferCtx%externalTemperature)
      Nullify(HeatXferCtx%DM)
      Call DMDestroy(HeatXferCtx%cellDM,ierr);CHKERRQ(ierr)
   End Subroutine MEF90HeatXferCtxDestroy

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90HeatXferCtxGlobalOptions"
!!!
!!!  
!!!  PetscBagRegisterMEF90HeatXferCtxGlobalOptions:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine PetscBagRegisterMEF90HeatXferCtxGlobalOptions(bag,name,prefix,default,ierr)
      PetscBag                                           :: bag
      Character(len=*),Intent(IN)                        :: prefix,name
      Type(MEF90HeatXferGlobalOptions_Type),Intent(IN)   :: default
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(MEF90HeatXferGlobalOptions_Type),pointer      :: HeatXferGlobalOptions
      Call PetscBagGetDataMEF90HeatXferCtxGlobalOptions(bag,HeatXferGlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagSetName(bag,trim(name),"HeatXferGlobalOptions MEF90 Heat transfer global options",ierr);CHKERRQ(ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix),ierr);CHKERRQ(ierr)

      Call PetscBagRegisterEnum(bag,HeatXferGlobalOptions%mode,MEF90HeatXFer_ModeList,default%mode,'heatxfer_mode','Type of heat transfer computation',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterBool(bag,HeatXferGlobalOptions%addNullSpace,default%addNullSpace,'heatxfer_addNullSpace','Add null space to SNES',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,HeatXferGlobalOptions%tempOffset,default%tempOffset,'temp_Offset','Position of temperature field in EXO file',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,HeatXferGlobalOptions%initialTemperature,default%initialTemperature,'heatxfer_initialTemp','[K] (T): Initial Temperature' ,ierr);CHKERRQ(ierr)

      Call PetscBagRegisterEnum(bag,HeatXferGlobalOptions%boundaryTempScaling,MEF90ScalingList,default%boundaryTempScaling,'boundaryTemp_scaling','Boundary temperature scaling',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,HeatXferGlobalOptions%boundaryTempOffset,default%boundaryTempOffset,'boundaryTemp_Offset','Position of boundary temperature field in EXO file',ierr);CHKERRQ(ierr)

      Call PetscBagRegisterEnum(bag,HeatXferGlobalOptions%externalTempScaling,MEF90ScalingList,default%externalTempScaling,'externalTemp_scaling','External Temperature scaling',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,HeatXferGlobalOptions%externalTempOffset,default%externalTempOffset,'externalTemp_Offset','Position of external temperature field in EXO file',ierr);CHKERRQ(ierr)

      Call PetscBagRegisterEnum(bag,HeatXferGlobalOptions%fluxScaling,MEF90ScalingList,default%fluxScaling,'flux_scaling','Heat flux scaling',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt(bag,HeatXferGlobalOptions%fluxOffset,default%fluxOffset,'flux_Offset','Position of heat flux field in EXO file',ierr);CHKERRQ(ierr)      

   End Subroutine PetscBagRegisterMEF90HeatXferCtxGlobalOptions

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90HeatXferCtxCellSetOptions"
!!!
!!!  
!!!  PetscBagRegisterMEF90HeatXferCtxCellSetOptions:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine PetscBagRegisterMEF90HeatXferCtxCellSetOptions(bag,name,prefix,default,ierr)
      PetscBag                                           :: bag
      Character(len=*),Intent(IN)                        :: prefix,name
      Type(MEF90HeatXferCellSetOptions_Type),Intent(IN)  :: default
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(MEF90HeatXferCellSetOptions_Type),pointer      :: HeatXferCellSetOptions
      Call PetscBagGetDataMEF90HeatXferCtxCellSetOptions(bag,HeatXferCellSetOptions,ierr);CHKERRQ(ierr)
      Call PetscBagSetName(bag,trim(name),"HeatXferCellSetOptions MEF90 Heat transfer Cell Set options",ierr);CHKERRQ(ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix),ierr);CHKERRQ(ierr)

      Call PetscBagRegisterInt(bag,HeatXferCellSetOptions%ElemTypeShortID,default%ElemTypeShortID,'ShortID','Element type ShortID',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,HeatXferCellSetOptions%Flux,default%Flux,'Flux','[J.s^(-1).m^(-3) / J.s^(-1).m^(-2)] (f): Internal / applied heat flux',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,HeatXferCellSetOptions%SurfaceThermalConductivity,default%SurfaceThermalConductivity,'SurfaceThermalConductivity','[J.s^(-2).m^(-1).K^(-1)] (H) Surface Thermal Conductivity',ierr)
      Call PetscBagRegisterReal(bag,HeatXferCellSetOptions%externalTemp,default%externalTemp,'externalTemp','Reference temperature T [K]',ierr)
      Call PetscBagRegisterBool(bag,HeatXferCellSetOptions%Has_BC,default%Has_BC,'TempBC','Temperature has Dirichlet boundary Condition (Y/N)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,HeatXferCellSetOptions%boundaryTemp,default%boundaryTemp,'boundaryTemp','Temperature boundary value',ierr);CHKERRQ(ierr)
   End Subroutine PetscBagRegisterMEF90HeatXferCtxCellSetOptions

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90HeatXferCtxVertexSetOptions"
!!!
!!!  
!!!  PetscBagRegisterMEF90HeatXferCtxVertexSetOptions:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine PetscBagRegisterMEF90HeatXferCtxVertexSetOptions(bag,name,prefix,default,ierr)
      PetscBag                                              :: bag
      Character(len=*),Intent(IN)                           :: prefix,name
      Type(MEF90HeatXferVertexSetOptions_Type),Intent(IN)   :: default
      PetscErrorCode,Intent(OUT)                            :: ierr

      Type(MEF90HeatXferVertexSetOptions_Type),pointer      :: HeatXferVertexSetOptions
      Call PetscBagGetDataMEF90HeatXferCtxVertexSetOptions(bag,HeatXferVertexSetOptions,ierr);CHKERRQ(ierr)
      Call PetscBagSetName(bag,trim(name),"HeatXferVertexSetOptions MEF90 Heat transfer Vertex Set options",ierr);CHKERRQ(ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix),ierr);CHKERRQ(ierr)

      Call PetscBagRegisterBool(bag,HeatXferVertexSetOptions%Has_BC,default%Has_BC,'TempBC','Temperature has Dirichlet boundary Condition (Y/N)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,HeatXferVertexSetOptions%boundaryTemp,default%boundaryTemp,'boundaryTemp','Temperature boundary value',ierr);CHKERRQ(ierr)
   End Subroutine PetscBagRegisterMEF90HeatXferCtxVertexSetOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferCtxSetFromOptions"
!!!
!!!  
!!!  MEF90HeatXferCtxSetFromOptions:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferCtxSetFromOptions(heatXferCtx,prefix,defaultGlobalOptions, &
                                              defaultCellSetOptions,defaultVertexSetOptions,ierr)
      Type(MEF90HeatXferCtx_Type),Intent(OUT)               :: heatXferCtx
      Character(len=*),Intent(IN)                           :: prefix
      Type(MEF90HeatXferGlobalOptions_Type),Intent(IN)      :: defaultGlobalOptions
      Type(MEF90HeatXferCellSetOptions_Type),Intent(IN)     :: defaultCellSetOptions
      Type(MEF90HeatXferVertexSetOptions_Type),Intent(IN)   :: defaultVertexSetOptions
      PetscErrorCode,Intent(OUT)                            :: ierr
   
      Type(MEF90CtxGlobalOptions_Type),pointer              :: MEF90CtxGlobalOptions
      Type(MEF90HeatXferCellSetOptions_Type)                :: myDefaultCellSetOptions
      Type(MEF90Element_Type),Dimension(:),Pointer          :: ElemType
      Type(IS)                                              :: setIS
      PetscInt,Dimension(:),Pointer                         :: setID
      PetscInt                                              :: set
      Character(len=MEF90_MXSTRLEN)                         :: IOBuffer,setName,setprefix

      Call PetscBagGetDataMEF90CtxGlobalOptions(heatXferCtx%MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr);CHKERRQ(ierr)
      !!!
      !!! Registering Global Context
      !!!
      Call PetscBagRegisterMEF90HeatXferCtxGlobalOptions(heatXferCtx%GlobalOptionsBag,"MEF90HeatXfer Global Ctx",prefix,defaultGlobalOptions,ierr);CHKERRQ(ierr)

      If (MEF90CtxGlobalOptions%verbose > 0) Then
         Call PetscBagView(heatXferCtx%GlobalOptionsBag,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
      End If

      !!!
      !!! Registering Cell Set Context
      !!! We override the default element type with that detected from the exodus file
      !!!
      Call DMmeshGetLabelIdIS(heatXferCtx%DM,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(heatXferCtx%MEF90Ctx%comm,setIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      
      Call EXOGetCellSetElementType_Scal(heatXferCtx%MEF90Ctx,ElemType,ierr)
      Do set = 1, size(setID)
         Write(setName,100) setID(set)
         Write(setprefix,101) setID(set)
         mydefaultCellSetOptions = defaultCellSetOptions
         mydefaultCellSetOptions%ElemTypeShortID = ElemType(set)%ShortID
         Call PetscBagRegisterMEF90HeatXferCtxCellSetOptions(heatXferCtx%CellSetOptionsBag(set),setName,setPrefix,mydefaultCellSetOptions,ierr)
         If (MEF90CtxGlobalOptions%verbose > 0) Then
            Write(IOBuffer,103) setID(set),trim(setprefix)
            Call PetscPrintf(heatXferCtx%MEF90Ctx%comm,IOBuffer,ierr);CHKERRQ(ierr)
            Call PetscBagView(heatXferCtx%CellSetOptionsBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
            Call PetscPrintf(heatXferCtx%MEF90Ctx%comm,"\n",ierr);CHKERRQ(ierr)
         End if
      End Do
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      DeAllocate(ElemType)
      
      !!!
      !!! Registering Vertex Set Context
      !!! We override the default element type with that detected from the exodus file
      !!!
      Call DMmeshGetLabelIdIS(heatXferCtx%DM,'Vertex Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      
      Do set = 1, size(setID)
         Write(setName,200) setID(set)
         Write(setprefix,201) setID(set)
         Call PetscBagRegisterMEF90HeatXferCtxVertexSetOptions(heatXferCtx%VertexSetOptionsBag(set),setName,setPrefix,defaultVertexSetOptions,ierr)
         If (MEF90CtxGlobalOptions%verbose > 0) Then
            Write(IOBuffer,203) setID(set),trim(setprefix)
            Call PetscPrintf(heatXferCtx%MEF90Ctx%comm,IOBuffer,ierr);CHKERRQ(ierr)
            Call PetscBagView(heatXferCtx%VertexSetOptionsBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
            Call PetscPrintf(heatXferCtx%MEF90Ctx%comm,"\n",ierr);CHKERRQ(ierr)
         End if
      End Do
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)

100 Format('Cell set ',I4)
101 Format('cs',I4.4,'_')
103 Format('\nRegistering cell set ',I4,' prefix: ',A,'\n')
200 Format('Vertex set ',I4)
201 Format('vs',I4.4,'_')
203 Format('\nRegistering vertex set ',I4,' prefix: ',A,'\n')
   End Subroutine MEF90HeatXferCtxSetFromOptions
   
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferCtxSetSections"
!!!
!!!  
!!!  MEF90HeatXferCtxSetSections:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferCtxSetSections(MEF90HeatXferCtx,ierr)
      Type(MEF90HeatXferCtx_Type),Intent(OUT)               :: MEF90HeatXferCtx
      PetscErrorCode,Intent(OUT)                            :: ierr
      
      Type(SectionReal)                                     :: defaultSection

      Call DMMeshGetVertexSectionReal(MEF90HeatXferCtx%DM,"default",1,defaultSection,ierr);CHKERRQ(ierr)
      Call DMMeshSetSectionReal(MEF90HeatXferCtx%DM,"default",defaultSection,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(defaultSection,ierr);CHKERRQ(ierr)
      Call DMSetBlockSize(MEF90HeatXferCtx%DM,1,ierr);CHKERRQ(ierr)

      Call DMMeshGetCellSectionReal(MEF90HeatXferCtx%cellDM,"default",1,defaultSection,ierr);CHKERRQ(ierr)
      Call DMMeshSetSectionReal(MEF90HeatXferCtx%cellDM,"default",defaultSection,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(defaultSection,ierr);CHKERRQ(ierr)
      Call DMSetBlockSize(MEF90HeatXferCtx%CellDM,1,ierr);CHKERRQ(ierr)
   End Subroutine MEF90HeatXferCtxSetSections


#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferCtxCreateVectors"
!!!
!!!  
!!!  MEF90HeatXferCtxCreateVectors:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferCtxCreateVectors(MEF90HeatXferCtx,ierr)
      Type(MEF90HeatXferCtx_Type),Intent(OUT)               :: MEF90HeatXferCtx
      PetscErrorCode,Intent(OUT)                            :: ierr

      Allocate(MEF90HeatXferCtx%temperature)
      Call DMCreateGlobalVector(MEF90HeatXferCtx%DM,MEF90HeatXferCtx%temperature,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(MEF90HeatXferCtx%temperature,"Temperature",ierr);CHKERRQ(ierr)
      Call VecSet(MEF90HeatXferCtx%temperature,0.0_Kr,ierr);CHKERRQ(ierr)
   
      Allocate(MEF90HeatXferCtx%boundaryTemperature)
      Call DMCreateGlobalVector(MEF90HeatXferCtx%DM,MEF90HeatXferCtx%boundaryTemperature,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(MEF90HeatXferCtx%boundaryTemperature,"Boundary_Temperature",ierr);CHKERRQ(ierr)
      Call VecSet(MEF90HeatXferCtx%boundaryTemperature,0.0_Kr,ierr);CHKERRQ(ierr)
   
      Allocate(MEF90HeatXferCtx%externalTemperature)
      Call DMCreateGlobalVector(MEF90HeatXferCtx%cellDM,MEF90HeatXferCtx%externalTemperature,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(MEF90HeatXferCtx%externalTemperature,"External_Temperature",ierr);CHKERRQ(ierr)
      Call VecSet(MEF90HeatXferCtx%externalTemperature,0.0_Kr,ierr);CHKERRQ(ierr)

      Allocate(MEF90HeatXferCtx%flux)
      Call DMCreateGlobalVector(MEF90HeatXferCtx%cellDM,MEF90HeatXferCtx%flux,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(MEF90HeatXferCtx%flux,"Flux",ierr);CHKERRQ(ierr)
      Call VecSet(MEF90HeatXferCtx%flux,0.0_Kr,ierr);CHKERRQ(ierr)
   End Subroutine MEF90HeatXferCtxCreateVectors   

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferCtxDestroyVectors"
!!!
!!!  
!!!  MEF90HeatXferCtxDestroyVectors:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferCtxDestroyVectors(MEF90HeatXferCtx,ierr)
      Type(MEF90HeatXferCtx_Type),Intent(OUT)               :: MEF90HeatXferCtx
      PetscErrorCode,Intent(OUT)                            :: ierr
   
      If (Associated(MEF90HeatXferCtx%temperature)) Then 
         Call VecDestroy(MEF90HeatXferCtx%temperature,ierr);CHKERRQ(ierr)
         DeAllocate(MEF90HeatXferCtx%temperature)
         Nullify(MEF90HeatXferCtx%temperature)
      End If   

      If (Associated(MEF90HeatXferCtx%boundaryTemperature)) Then 
         Call VecDestroy(MEF90HeatXferCtx%boundaryTemperature,ierr);CHKERRQ(ierr)
         DeAllocate(MEF90HeatXferCtx%boundaryTemperature)
         Nullify(MEF90HeatXferCtx%boundaryTemperature)
      End If   

      If (Associated(MEF90HeatXferCtx%externalTemperature)) Then 
         Call VecDestroy(MEF90HeatXferCtx%externalTemperature,ierr);CHKERRQ(ierr)
         DeAllocate(MEF90HeatXferCtx%externalTemperature)
         Nullify(MEF90HeatXferCtx%externalTemperature)
      End If   

      If (Associated(MEF90HeatXferCtx%flux)) Then 
         Call VecDestroy(MEF90HeatXferCtx%flux,ierr);CHKERRQ(ierr)
         DeAllocate(MEF90HeatXferCtx%flux)
         Nullify(MEF90HeatXferCtx%flux)
      End If   
   End Subroutine MEF90HeatXferCtxDestroyVectors
End Module m_MEF90_HeatXferCtx
