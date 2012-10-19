!!!
!!! What do I do with initial temperature?
!!!
#include "../MEF90/mef90.inc"
Module m_MEF90_HeatXferCtx_Type
#include "finclude/petscdef.h"
   Use m_MEF90
   Use,Intrinsic :: iso_c_binding
   Implicit none
   !Private  
   !Public :: MEF90HeatXferCtx_Type
   !Public :: MEF90HeatXferGlobalOptions_Type
   !Public :: MEF90HeatXferCellSetOptions_Type
   !Public :: MEF90HeatXferVertexSetOptions_Type
   
   Type MEF90HeatXferCtx_Type
      PetscReal                        :: timePrevious,timeTarget

      Type(Vec),pointer                :: fluxPrevious,fluxTarget,Flux
      Type(Vec),pointer                :: boundaryValuePrevious,boundaryValueTarget,boundaryValue
      Type(Vec),pointer                :: externalValuePrevious,externalValueTarget,externalValue
      !!! XXXPrevious     represents a field at the previous time step
      !!! XXXTarget  represents a field at the time step currently being computed
      !!! XXX        represents a field at current time, interpolated from XXXPrevious and XXXTarget
      
      PetscBag                         :: GlobalOptionsBag
      PetscBag,Dimension(:),Pointer    :: CellSetOptionsBag
      PetscBag,Dimension(:),Pointer    :: VertexSetOptionsBag
      Type(MEF90Ctx_Type),pointer      :: MEF90Ctx
      Type(DM),pointer                 :: DM
   End Type MEF90HeatXferCtx_Type
   
   Type MEF90HeatXferGlobalOptions_Type
      PetscInt                         :: mode
      PetscBool                        :: addNullSpace
      PetscInt                         :: tempOffset
      PetscInt                         :: boundaryTempOffset
      PetscInt                         :: externalTempOffset
      PetscInt                         :: fluxOffset
   End Type MEF90HeatXferGlobalOptions_Type

   Type MEF90HeatXferCellSetOptions_Type
      PetscInt                         :: elemTypeShortID
      PetscReal                        :: flux
      PetscReal                        :: surfaceThermalConductivity
      PetscReal                        :: externalTemp
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
   Public :: PetscBagGetDataMEF90HeatXferGlobalOptions
   
   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF90_HeatXferCtx_Type
         PetscBag                                           :: bag
         Type(MEF90HeatXferGlobalOptions_Type),pointer   :: data
         PetscErrorCode                                     :: ierr
      End subroutine PetscBagGetData
   End interface
Contains

#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90HeatXferGlobalOptions"
!!!
!!!  PetscBagGetDataMEF90HeatXferGlobalOptions - Custom interface to PetscGetData
!!!
   Subroutine PetscBagGetDataMEF90HeatXferGlobalOptions(bag,data,ierr)
      PetscBag                                        :: bag
      Type(MEF90HeatXferGlobalOptions_Type),pointer     :: data
      PetscErrorCode                                  :: ierr
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataMEF90HeatXferGlobalOptions
End Module m_MEF90HeatXferGlobalOptions_Private

Module m_MEF90HeatXferCellSetOptions_Private
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_HeatXferCtx_Type
   Implicit None

   Private
   Public :: PetscBagGetDataMEF90HeatXferCellSetOptions
   
   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF90_HeatXferCtx_Type
         PetscBag                                           :: bag
         Type(MEF90HeatXferCellSetOptions_Type),pointer   :: data
         PetscErrorCode                                     :: ierr
      End subroutine PetscBagGetData
   End interface
Contains

#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90HeatXferCellSetOptions"
!!!
!!!  PetscBagGetDataMEF90HeatXferCellSetOptions - Custom interface to PetscGetData
!!!
   Subroutine PetscBagGetDataMEF90HeatXferCellSetOptions(bag,data,ierr)
      PetscBag                                        :: bag
      Type(MEF90HeatXferCellSetOptions_Type),pointer     :: data
      PetscErrorCode                                  :: ierr
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataMEF90HeatXferCellSetOptions
End Module m_MEF90HeatXferCellSetOptions_Private

Module m_MEF90HeatXferVertexSetOptions_Private
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_HeatXferCtx_Type
   Implicit None

   Private
   Public :: PetscBagGetDataMEF90HeatXferVertexSetOptions
   
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
#define __FUNCT__ "PetscBagGetDataMEF90HeatXferVertexSetOptions"
!!!
!!!  PetscBagGetDataMEF90HeatXferVertexSetOptions - Custom interface to PetscGetData
!!!
   Subroutine PetscBagGetDataMEF90HeatXferVertexSetOptions(bag,data,ierr)
      PetscBag                                        :: bag
      Type(MEF90HeatXferVertexSetOptions_Type),pointer     :: data
      PetscErrorCode                                  :: ierr
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataMEF90HeatXferVertexSetOptions
End Module m_MEF90HeatXferVertexSetOptions_Private

Module m_MEF90_HeatXferCtx
#include "finclude/petscdef.h"
   Use m_MEF90
   Use m_MEF90_HeatXferCtx_Type
   Use m_MEF90HeatXferGlobalOptions_Private
   Use m_MEF90HeatXferCellSetOptions_Private
   Use m_MEF90HeatXferVertexSetOptions_Private
   Implicit none

   PetscInt,protected   :: sizeofMEF90HeatXferGlobalOptions
   PetscInt,protected   :: sizeofMEF90HeatXferCellSetOptions
   PetscInt,protected   :: sizeofMEF90HeatXferVertexSetOptions
   
   Enum,bind(c)
      enumerator  :: MEF90HeatXFer_ModeSteadyState = 0, &
                     MEF90HeatXFer_ModeSteadyTransient
   End Enum
   Character(len = MEF90_MXSTRLEN),dimension(5),protected   :: MEF90HeatXFer_ModeList
   
Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferCtx_InitializePrivate"
!!!
!!!  
!!!  MEF90HeatXferCtx_InitializePrivate:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferCtx_InitializePrivate(ierr)
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

      MEF90HeatXFer_ModeList(1) = 'SteadyState'
      MEF90HeatXFer_ModeList(2) = 'Transient'
      MEF90HeatXFer_ModeList(3) = 'MEF90_HeatXFer_Mode'
      MEF90HeatXFer_ModeList(4) = '_MEF90_HeatXFer_Mode'
      MEF90HeatXFer_ModeList(5) = ''
   End Subroutine MEF90HeatXferCtx_InitializePrivate
   
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferCtx_Create"
!!!
!!!  
!!!  MEF90HeatXferCtx_Create:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferCtx_Create(HeatXferCtx,Mesh,MEF90Ctx,ierr)
      Type(MEF90HeatXferCtx_Type),Intent(OUT)            :: HeatXferCtx
      Type(DM),target,Intent(IN)                         :: Mesh
      Type(MEF90Ctx_Type),target,Intent(IN)              :: MEF90Ctx
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90CtxGlobalOptions
      Type(MEF90HeatXferGlobalOptions_Type),pointer      :: MEF90HeatXferGlobalOptions
      Type(MEF90HeatXferCellSetOptions_Type),pointer     :: MEF90HeatXferCellSetOptions
      Type(MEF90HeatXferVertexSetOptions_Type),pointer   :: MEF90HeatXferVertexSetOptions
      PetscInt                                           :: set,numSet

      Call MEF90HeatXferCtx_InitializePrivate(ierr)
      HeatXferCtx%DM => Mesh
      HeatXferCtx%MEF90Ctx => MEF90Ctx
      
      Call PetscBagCreate(MEF90Ctx%comm,sizeofMEF90HeatXferGlobalOptions,HeatXferCtx%GlobalOptionsBag,ierr);CHKERRQ(ierr)
      Call DMmeshGetLabelSize(Mesh,'Cell Sets',numSet,ierr);CHKERRQ(ierr)
      Allocate(HeatXferCtx%CellSetOptionsBag(numSet),stat=ierr)
      Do set = 1, numSet
         Call PetscBagCreate(MEF90Ctx%comm,sizeofMEF90HeatXferCellSetOptions,HeatXferCtx%CellSetOptionsBag(set),ierr);CHKERRQ(ierr)
      End Do

      Call DMmeshGetLabelSize(Mesh,'Vertex Sets',numSet,ierr);CHKERRQ(ierr)
      Allocate(HeatXferCtx%VertexSetOptionsBag(numSet),stat=ierr)
      Do set = 1, numSet
         Call PetscBagCreate(MEF90Ctx%comm,sizeofMEF90HeatXferVertexSetOptions,HeatXferCtx%VertexSetOptionsBag(set),ierr);CHKERRQ(ierr)
      End Do
      
      Nullify(HeatXferCtx%fluxPrevious)
      Nullify(HeatXferCtx%fluxTarget)
      Nullify(HeatXferCtx%flux)
      Nullify(HeatXferCtx%boundaryValuePrevious)
      Nullify(HeatXferCtx%boundaryValueTarget)
      Nullify(HeatXferCtx%boundaryValue)
      Nullify(HeatXferCtx%externalValuePrevious)
      Nullify(HeatXferCtx%externalValueTarget)
      Nullify(HeatXferCtx%externalValue)
      Nullify(HeatXferCtx%DM)
      Nullify(HeatXferCtx%DM)
   End Subroutine MEF90HeatXferCtx_Create
   
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferCtx_Destroy"
!!!
!!!  
!!!  MEF90HeatXferCtx_Destroy:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90HeatXferCtx_Destroy(HeatXferCtx,ierr)
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
      Nullify(HeatXferCtx%fluxPrevious)
      Nullify(HeatXferCtx%fluxTarget)
      Nullify(HeatXferCtx%flux)
      Nullify(HeatXferCtx%boundaryValuePrevious)
      Nullify(HeatXferCtx%boundaryValueTarget)
      Nullify(HeatXferCtx%boundaryValue)
      Nullify(HeatXferCtx%externalValuePrevious)
      Nullify(HeatXferCtx%externalValueTarget)
      Nullify(HeatXferCtx%externalValue)
      Nullify(HeatXferCtx%DM)
      Nullify(HeatXferCtx%DM)
   End Subroutine MEF90HeatXferCtx_Destroy

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90HeatXferGlobalOptions"
!!!
!!!  
!!!  PetscBagRegisterMEF90HeatXferGlobalOptions:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine PetscBagRegisterMEF90HeatXferGlobalOptions(bag,name,prefix,default,ierr)
      PetscBag                                           :: bag
      Character(len=*),Intent(IN)                        :: prefix,name
      Type(MEF90HeatXferGlobalOptions_Type),Intent(IN)   :: default
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(MEF90HeatXferGlobalOptions_Type),pointer      :: HeatXferGlobalOptions
      Call PetscBagGetDataMEF90HeatXferGlobalOptions(bag,HeatXferGlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagSetName(bag,trim(name),"HeatXferGlobalOptions MEF90 Heat transfer global options",ierr);CHKERRQ(ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix),ierr);CHKERRQ(ierr)

      Call PetscBagRegisterEnum(bag,HeatXferGlobalOptions%mode,MEF90HeatXFer_ModeList,default%mode,'heatxfer_mode','Type of heat transfer computation',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterBool(bag,HeatXferGlobalOptions%addNullSpace,default%addNullSpace,'addNullSpace','Add null space to SNES',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt(bag,HeatXferGlobalOptions%tempOffset,default%tempOffset,'temp_Offset','Position of temperature field in EXO file',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt(bag,HeatXferGlobalOptions%boundaryTempOffset,default%boundaryTempOffset,'boundaryTemp_Offset','Position of boundary temperature field in EXO file',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt(bag,HeatXferGlobalOptions%externalTempOffset,default%externalTempOffset,'externalTemp_Offset','Position of external temperature field in EXO file',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt(bag,HeatXferGlobalOptions%fluxOffset,default%fluxOffset,'flux_Offset','Position of heat flux field in EXO file',ierr);CHKERRQ(ierr)      
   End Subroutine PetscBagRegisterMEF90HeatXferGlobalOptions

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90HeatXferCellSetOptions"
!!!
!!!  
!!!  PetscBagRegisterMEF90HeatXferCellSetOptions:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine PetscBagRegisterMEF90HeatXferCellSetOptions(bag,name,prefix,default,ierr)
      PetscBag                                           :: bag
      Character(len=*),Intent(IN)                        :: prefix,name
      Type(MEF90HeatXferCellSetOptions_Type),Intent(IN)  :: default
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(MEF90HeatXferCellSetOptions_Type),pointer      :: HeatXferCellSetOptions
      Call PetscBagGetDataMEF90HeatXferCellSetOptions(bag,HeatXferCellSetOptions,ierr);CHKERRQ(ierr)
      Call PetscBagSetName(bag,trim(name),"HeatXferCellSetOptions MEF90 Heat transfer Cell Set options",ierr);CHKERRQ(ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix),ierr);CHKERRQ(ierr)

      Call PetscBagRegisterInt(bag,HeatXferCellSetOptions%ElemTypeShortID,default%ElemTypeShortID,'ShortID','Element type ShortID',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,HeatXferCellSetOptions%Flux,default%Flux,'Flux','[J.s^(-1).m^(-3) / J.s^(-1).m^(-2)] (f): Internal / applied heat flux',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,HeatXferCellSetOptions%SurfaceThermalConductivity,default%SurfaceThermalConductivity,'SurfaceThermalConductivity','[J.s^(-1).m^(-2).K^(-1)] (H) Surface Thermal Conductivity',ierr)
      Call PetscBagRegisterReal(bag,HeatXferCellSetOptions%externalTemp,default%externalTemp,'externalTemp','Reference temperature T [K]',ierr)
   End Subroutine PetscBagRegisterMEF90HeatXferCellSetOptions

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90HeatXferVertexSetOptions"
!!!
!!!  
!!!  PetscBagRegisterMEF90HeatXferVertexSetOptions:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine PetscBagRegisterMEF90HeatXferVertexSetOptions(bag,name,prefix,default,ierr)
      PetscBag                                              :: bag
      Character(len=*),Intent(IN)                           :: prefix,name
      Type(MEF90HeatXferVertexSetOptions_Type),Intent(IN)   :: default
      PetscErrorCode,Intent(OUT)                            :: ierr

      Type(MEF90HeatXferVertexSetOptions_Type),pointer      :: HeatXferVertexSetOptions
      Call PetscBagGetDataMEF90HeatXferVertexSetOptions(bag,HeatXferVertexSetOptions,ierr);CHKERRQ(ierr)
      Call PetscBagSetName(bag,trim(name),"HeatXferVertexSetOptions MEF90 Heat transfer Vertex Set options",ierr);CHKERRQ(ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix),ierr);CHKERRQ(ierr)

      Call PetscBagRegisterBool(bag,HeatXferVertexSetOptions%Has_BC,default%Has_BC,'TempBC','Temperature has Dirichlet boundary Condition (Y/N)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,HeatXferVertexSetOptions%boundaryTemp,default%boundaryTemp,'boundaryTemp','Temperature boundary value',ierr);CHKERRQ(ierr)
   End Subroutine PetscBagRegisterMEF90HeatXferVertexSetOptions

End Module m_MEF90_HeatXferCtx
