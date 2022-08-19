Module m_MEF90_HeatXferCtx_Type
#include "petsc/finclude/petsc.h"
   Use m_MEF90
   Use,Intrinsic :: iso_c_binding
   Implicit none
   Private  
   Public :: MEF90HeatXferCtx_Type
   Public :: MEF90HeatXferGlobalOptions_Type
   Public :: MEF90HeatXferCellSetOptions_Type
   Public :: MEF90HeatXferVertexSetOptions_Type
   
   Type MEF90HeatXferCtx_Type
      Type(tVec),pointer               :: temperatureLocal
      Type(tVec),pointer               :: boundaryTemperatureLocal
      Type(tVec),pointer               :: fluxLocal
      Type(tVec),pointer               :: boundaryFluxLocal
      Type(tVec),pointer               :: externalTemperatureLocal

      PetscBag                         :: GlobalOptionsBag
      PetscBag,Dimension(:),Pointer    :: CellSetOptionsBag
      PetscBag,Dimension(:),Pointer    :: FaceSetOptionsBag
      PetscBag,Dimension(:),Pointer    :: VertexSetOptionsBag
      PetscBag,Dimension(:),Pointer    :: MaterialPropertiesBag
      Type(MEF90Ctx_Type),pointer      :: MEF90Ctx
      Type(tDM),pointer                :: megaDM
   End Type MEF90HeatXferCtx_Type
   
   Type MEF90HeatXferGlobalOptions_Type
      PetscEnum                        :: timeSteppingType
      PetscBool                        :: addNullSpace
      PetscReal                        :: initialTemperature
      PetscEnum                        :: boundaryTempScaling
      PetscEnum                        :: externalTempScaling
      PetscEnum                        :: fluxScaling
      PetscEnum                        :: boundaryFluxScaling
      !!! scaling = time (step) scaling law currently CST, Linear, Null (not present), File
   End Type MEF90HeatXferGlobalOptions_Type

   Type MEF90HeatXferCellSetOptions_Type
      PetscReal                        :: flux
      PetscReal                        :: surfaceThermalConductivity
      PetscReal                        :: externalTemp
      PetscBool                        :: Has_BC
      PetscReal                        :: boundaryTemp
      PetscReal,dimension(3)           :: advectionVector
   End Type MEF90HeatXferCellSetOptions_Type

   Type MEF90HeatXferVertexSetOptions_Type
      PetscBool                        :: Has_BC
      PetscReal                        :: boundaryTemp
   End Type MEF90HeatXferVertexSetOptions_Type 
End Module m_MEF90_HeatXferCtx_Type

Module m_MEF90HeatXferGlobalOptions_Private
#include "petsc/finclude/petsc.h"
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
      
      PetscCall(PetscBagGetData(bag,data,ierr))
   End Subroutine PetscBagGetDataMEF90HeatXferCtxGlobalOptions
End Module m_MEF90HeatXferGlobalOptions_Private

Module m_MEF90HeatXferCellSetOptions_Private
#include "petsc/finclude/petsc.h"
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
      
      PetscCall(PetscBagGetData(bag,data,ierr))
   End Subroutine PetscBagGetDataMEF90HeatXferCtxCellSetOptions
End Module m_MEF90HeatXferCellSetOptions_Private

Module m_MEF90HeatXferVertexSetOptions_Private
#include "petsc/finclude/petsc.h"
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
      
      PetscCall(PetscBagGetData(bag,data,ierr))
   End Subroutine PetscBagGetDataMEF90HeatXferCtxVertexSetOptions
End Module m_MEF90HeatXferVertexSetOptions_Private

Module m_MEF90_HeatXferCtx
#include "petsc/finclude/petsc.h"
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
      enumerator  :: MEF90HeatXfer_timeSteppingTypeNULL = 0,    &
                     MEF90HeatXFer_timeSteppingTypeSteadyState, &
                     MEF90HeatXFer_timeSteppingTypeTransient
   End Enum
   Character(len = MEF90MXSTRLEN),dimension(6),protected   :: MEF90HeatXFer_timeSteppingTypeList
   
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
      Character(len=1),pointer                        :: dummychar(:)
      PetscSizeT                                      :: sizeofchar
      
      PetscCall(PetscDataTypeGetSize(PETSC_CHAR,sizeofchar,ierr))
      sizeofMEF90HeatXferGlobalOptions = size(transfer(HeatXferGlobalOptions,dummychar))*sizeofchar
      sizeofMEF90HeatXferCellSetOptions = size(transfer(HeatXferCellSetOptions,dummychar))*sizeofchar
      sizeofMEF90HeatXferVertexSetOptions = size(transfer(HeatXferVertexSetOptions,dummychar))*sizeofchar

      MEF90HeatXFer_timeSteppingTypeList(1) = 'null'
      MEF90HeatXFer_timeSteppingTypeList(2) = 'SteadyState'
      MEF90HeatXFer_timeSteppingTypeList(3) = 'Transient'
      MEF90HeatXFer_timeSteppingTypeList(4) = 'MEF90_HeatXFer_timeSteppingType'
      MEF90HeatXFer_timeSteppingTypeList(5) = '_MEF90_HeatXFer_timeSteppingType'
      MEF90HeatXFer_timeSteppingTypeList(6) = ''
   End Subroutine MEF90HeatXferCtxInitialize_Private
   
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferCtxCreate"
!!!
!!!  
!!!  MEF90HeatXferCtxCreate:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90HeatXferCtxCreate(HeatXferCtx,dm,MEF90Ctx,ierr)
      Type(MEF90HeatXferCtx_Type),Intent(OUT)            :: HeatXferCtx
      Type(tDM),target,Intent(IN)                        :: dm
      Type(MEF90Ctx_Type),target,Intent(IN)              :: MEF90Ctx
      PetscErrorCode,Intent(INOUT)                       :: ierr

      Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions   
      Type(tIS)                                          :: setIS
      PetscInt                                           :: set,numSet
      Character(len=MEF90MXSTRLEN)                       :: vecName

      Call MEF90HeatXferCtxInitialize_Private(ierr)
      HeatXferCtx%MEF90Ctx => MEF90Ctx
      PetscCall(PetscBagCreate(MEF90Ctx%comm,sizeofMEF90HeatXferGlobalOptions,HeatXferCtx%GlobalOptionsBag,ierr))
      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr))
      
      !!! I need to allocate for the overall number of sets, not the local one
      PetscCall(DMGetLabelIdIS(dm,MEF90CellSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr)) 
      PetscCall(ISGetLocalSize(setIS,numSet,ierr))
      Allocate(HeatXferCtx%CellSetOptionsBag(numSet),stat=ierr)
      Do set = 1, numSet
         PetscCall(PetscBagCreate(MEF90Ctx%comm,sizeofMEF90HeatXferCellSetOptions,HeatXferCtx%CellSetOptionsBag(set),ierr))
      End Do
      PetscCall(ISDestroy(setIS,ierr))

      PetscCall(DMGetLabelIdIS(dm,MEF90FaceSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr)) 
      PetscCall(ISGetLocalSize(setIS,numSet,ierr))
      Allocate(HeatXferCtx%FaceSetOptionsBag(numSet),stat=ierr)
      Do set = 1, numSet
         PetscCall(PetscBagCreate(MEF90Ctx%comm,sizeofMEF90HeatXferCellSetOptions,HeatXferCtx%FaceSetOptionsBag(set),ierr))
      End Do
      PetscCall(ISDestroy(setIS,ierr))

      PetscCall(DMGetLabelIdIS(dm,MEF90VertexSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr)) 
      PetscCall(ISGetLocalSize(setIS,numSet,ierr))
      Allocate(HeatXferCtx%VertexSetOptionsBag(numSet),stat=ierr)
      Do set = 1, numSet
         PetscCall(PetscBagCreate(MEF90Ctx%comm,sizeofMEF90HeatXferVertexSetOptions,HeatXferCtx%VertexSetOptionsBag(set),ierr))
      End Do
      PetscCall(ISDestroy(setIS,ierr))
      
      Allocate(HeatXferCtx%temperatureLocal)
      vecName = "Temperature"
      PetscCall(MEF90CreateLocalVector(dm,MEF90GlobalOptions%elementFamily,MEF90GlobalOptions%elementOrder,1_Ki,vecName,HeatXferCtx%temperatureLocal,ierr))

      Allocate(HeatXferCtx%fluxLocal)
      vecName = "Flux"
      PetscCall(MEF90CreateCellVector(dm,1_Ki,vecName,HeatXferCtx%fluxLocal,ierr))

      Allocate(HeatXferCtx%externalTemperatureLocal)
      vecName = "External Temperature"
      PetscCall(MEF90CreateCellVector(dm,1_Ki,vecName,HeatXferCtx%externalTemperatureLocal,ierr))

      !!! Skipping boundary fields for now
      !!! Allocate(HeatXferCtx%boundaryFluxLocal)
      !!! vecName = "Boundary Flux"
      !!!Allocate(HeatXferCtx%boundaryTemperatureLocal)
      !!!vecName = "Boundary Temperature"
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
   
      !Nullify(HeatXferCtx%dm)
      If (Associated(HeatXferCtx%temperatureLocal)) Then
         PetscCall(VecDestroy(HeatXferCtx%temperatureLocal,ierr))
         Nullify(HeatXferCtx%temperatureLocal)
      End If
      If (Associated(HeatXferCtx%boundaryTemperatureLocal)) Then
         PetscCall(VecDestroy(HeatXferCtx%boundaryTemperatureLocal,ierr))
         Nullify(HeatXferCtx%boundaryTemperatureLocal)
      End If
      If (Associated(HeatXferCtx%ExternalTemperatureLocal)) Then
         PetscCall(VecDestroy(HeatXferCtx%ExternalTemperatureLocal,ierr))
         Nullify(HeatXferCtx%ExternalTemperatureLocal)
      End If
      If (Associated(HeatXferCtx%fluxLocal)) Then
         PetscCall(VecDestroy(HeatXferCtx%fluxLocal,ierr))
         Nullify(HeatXferCtx%fluxLocal)
      End If
      If (Associated(HeatXferCtx%boundaryFluxLocal)) Then
         PetscCall(VecDestroy(HeatXferCtx%boundaryFluxLocal,ierr))
         Nullify(HeatXferCtx%boundaryFluxLocal)
      End If

      Nullify(HeatXferCtx%MEF90Ctx)

      PetscCall(PetscBagDestroy(HeatXferCtx%GlobalOptionsBag,ierr))
      Do set = 1, size(HeatXferCtx%CellSetOptionsBag)
         PetscCall(PetscBagDestroy(HeatXferCtx%CellSetOptionsBag(set),ierr))
      End Do
      DeAllocate(HeatXferCtx%CellSetOptionsBag)
      Do set = 1, size(HeatXferCtx%FaceSetOptionsBag)
         PetscCall(PetscBagDestroy(HeatXferCtx%FaceSetOptionsBag(set),ierr))
      End Do
      DeAllocate(HeatXferCtx%FaceSetOptionsBag)
      Do set = 1, size(HeatXferCtx%VertexSetOptionsBag)
         PetscCall(PetscBagDestroy(HeatXferCtx%VertexSetOptionsBag(set),ierr))
      End Do
      DeAllocate(HeatXferCtx%VertexSetOptionsBag)
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
      PetscCall(PetscBagGetDataMEF90HeatXferCtxGlobalOptions(bag,HeatXferGlobalOptions,ierr))
      PetscCall(PetscBagSetName(bag,trim(name),"HeatXferGlobalOptions MEF90 Heat transfer global options",ierr))
      PetscCall(PetscBagSetOptionsPrefix(bag,trim(prefix),ierr))

      PetscCall(PetscBagRegisterEnum(bag,HeatXferGlobalOptions%timeSteppingType,MEF90HeatXFer_timeSteppingTypeList,default%timeSteppingType,'heatxfer_timeStepping_type','Type of heat transfer computation',ierr))
      PetscCall(PetscBagRegisterBool(bag,HeatXferGlobalOptions%addNullSpace,default%addNullSpace,'heatxfer_addNullSpace','Add null space to SNES',ierr))
      PetscCall(PetscBagRegisterReal(bag,HeatXferGlobalOptions%initialTemperature,default%initialTemperature,'heatxfer_initialTemp','[K] (T): Initial Temperature' ,ierr))

      PetscCall(PetscBagRegisterEnum(bag,HeatXferGlobalOptions%boundaryTempScaling,MEF90ScalingList,default%boundaryTempScaling,'boundaryTemp_scaling','Boundary temperature scaling',ierr))

      PetscCall(PetscBagRegisterEnum(bag,HeatXferGlobalOptions%externalTempScaling,MEF90ScalingList,default%externalTempScaling,'externalTemp_scaling','External Temperature scaling',ierr))

      PetscCall(PetscBagRegisterEnum(bag,HeatXferGlobalOptions%fluxScaling,MEF90ScalingList,default%fluxScaling,'flux_scaling','Heat flux scaling',ierr))

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
      PetscCall(PetscBagGetDataMEF90HeatXferCtxCellSetOptions(bag,HeatXferCellSetOptions,ierr))
      PetscCall(PetscBagSetName(bag,trim(name),"HeatXferCellSetOptions MEF90 Heat transfer Cell Set options",ierr))
      PetscCall(PetscBagSetOptionsPrefix(bag,trim(prefix),ierr))

      HeatXferCellSetOptions%advectionVector = default%advectionVector
      PetscCall(PetscBagRegisterReal(bag,HeatXferCellSetOptions%Flux,default%Flux,'Flux','[J.s^(-1).m^(-3) / J.s^(-1).m^(-2) / J.s^(-1).m^(-1)] (f): Internal / applied heat flux',ierr))
      PetscCall(PetscBagRegisterReal(bag,HeatXferCellSetOptions%SurfaceThermalConductivity,default%SurfaceThermalConductivity,'SurfaceThermalConductivity','[J.s^(-1).m^(-2).K^(-1) / J.s^(-1).m^(-1).K^(-1) ] (H) Surface Thermal Conductivity',ierr))
      PetscCall(PetscBagRegisterReal(bag,HeatXferCellSetOptions%externalTemp,default%externalTemp,'externalTemp','Reference temperature T [K]',ierr))
      PetscCall(PetscBagRegisterBool(bag,HeatXferCellSetOptions%Has_BC,default%Has_BC,'TempBC','Temperature has Dirichlet boundary Condition (Y/N)',ierr))
      PetscCall(PetscBagRegisterReal(bag,HeatXferCellSetOptions%boundaryTemp,default%boundaryTemp,'boundaryTemp','Temperature boundary value',ierr))
      PetscCall(PetscBagRegisterRealArray(bag,HeatXferCellSetOptions%advectionVector,3,'advectionVector','[m.s^(-1)] (V): advection vector',ierr))
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
      PetscCall(PetscBagGetDataMEF90HeatXferCtxVertexSetOptions(bag,HeatXferVertexSetOptions,ierr))
      PetscCall(PetscBagSetName(bag,trim(name),"HeatXferVertexSetOptions MEF90 Heat transfer Vertex Set options",ierr))
      PetscCall(PetscBagSetOptionsPrefix(bag,trim(prefix),ierr))

      PetscCall(PetscBagRegisterBool(bag,HeatXferVertexSetOptions%Has_BC,default%Has_BC,'TempBC','Temperature has Dirichlet boundary Condition (Y/N)',ierr))
      PetscCall(PetscBagRegisterReal(bag,HeatXferVertexSetOptions%boundaryTemp,default%boundaryTemp,'boundaryTemp','Temperature boundary value',ierr))
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
      Type(tIS)                                             :: setIS
      PetscInt,Dimension(:),Pointer                         :: setID
      PetscInt                                              :: set
      Character(len=MEF90MXSTRLEN)                         :: IOBuffer,setName,setprefix

      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(heatXferCtx%MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr))
      !!!
      !!! Registering Global Context
      !!!
      PetscCall(PetscBagRegisterMEF90HeatXferCtxGlobalOptions(heatXferCtx%GlobalOptionsBag,"MEF90HeatXfer Global Ctx",prefix,defaultGlobalOptions,ierr))

      If (MEF90CtxGlobalOptions%verbose > 0) Then
         PetscCall(PetscBagView(heatXferCtx%GlobalOptionsBag,PETSC_VIEWER_STDOUT_WORLD,ierr))
      End If

      !!!
      !!! Registering Cell Set Context
      !!! We override the default element type with that detected from the exodus file
      !!!
      PetscCall(DMGetLabelIdIS(heatXferCtx%megaDM,MEF90CellSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(heatXferCtx%MEF90Ctx%comm,setIS,ierr)) 
      PetscCall(ISGetIndicesF90(setIS,setID,ierr))
      
      Do set = 1, size(setID)
         Write(setName,100) setID(set)
         Write(setprefix,101) setID(set)
         mydefaultCellSetOptions = defaultCellSetOptions
         PetscCall(PetscBagRegisterMEF90HeatXferCtxCellSetOptions(heatXferCtx%CellSetOptionsBag(set),setName,setPrefix,mydefaultCellSetOptions,ierr))
         If (MEF90CtxGlobalOptions%verbose > 0) Then
            Write(IOBuffer,103) setID(set),trim(setprefix)
            PetscCall(PetscPrintf(heatXferCtx%MEF90Ctx%comm,IOBuffer,ierr))
            PetscCall(PetscBagView(heatXferCtx%CellSetOptionsBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr))
            PetscCall(PetscPrintf(heatXferCtx%MEF90Ctx%comm,"\n",ierr))
         End if
      End Do
      PetscCall(ISDestroy(setIS,ierr))
      
      !!!
      !!! Registering Face Set Context
      !!! We override the default element type with that detected from the exodus file
      !!!
      PetscCall(DMGetLabelIdIS(heatXferCtx%megaDM,MEF90FaceSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(heatXferCtx%MEF90Ctx%comm,setIS,ierr)) 
      PetscCall(ISGetIndicesF90(setIS,setID,ierr))
      
      Do set = 1, size(setID)
         Write(setName,300) setID(set)
         Write(setprefix,301) setID(set)
         mydefaultCellSetOptions = defaultCellSetOptions
         PetscCall(PetscBagRegisterMEF90HeatXferCtxCellSetOptions(heatXferCtx%FaceSetOptionsBag(set),setName,setPrefix,mydefaultCellSetOptions,ierr))
         If (MEF90CtxGlobalOptions%verbose > 0) Then
            Write(IOBuffer,303) setID(set),trim(setprefix)
            PetscCall(PetscPrintf(heatXferCtx%MEF90Ctx%comm,IOBuffer,ierr))
            PetscCall(PetscBagView(heatXferCtx%FaceSetOptionsBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr))
            PetscCall(PetscPrintf(heatXferCtx%MEF90Ctx%comm,"\n",ierr))
         End if
      End Do
      PetscCall(ISDestroy(setIS,ierr))
      
      !!!
      !!! Registering Vertex Set Context
      !!! We override the default element type with that detected from the exodus file
      !!!
      PetscCall(DMGetLabelIdIS(heatXferCtx%megaDM,MEF90VertexSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr)) 
      PetscCall(ISGetIndicesF90(setIS,setID,ierr))
      
      Do set = 1, size(setID)
         Write(setName,200) setID(set)
         Write(setprefix,201) setID(set)
         PetscCall(PetscBagRegisterMEF90HeatXferCtxVertexSetOptions(heatXferCtx%VertexSetOptionsBag(set),setName,setPrefix,defaultVertexSetOptions,ierr))
         If (MEF90CtxGlobalOptions%verbose > 0) Then
            Write(IOBuffer,203) setID(set),trim(setprefix)
            PetscCall(PetscPrintf(heatXferCtx%MEF90Ctx%comm,IOBuffer,ierr))
            PetscCall(PetscBagView(heatXferCtx%VertexSetOptionsBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr))
            PetscCall(PetscPrintf(heatXferCtx%MEF90Ctx%comm,"\n",ierr))
         End if
      End Do
      PetscCall(ISDestroy(setIS,ierr))

100 Format('Cell set ',I4)
101 Format('cs',I4.4,'_')
103 Format('\nRegistering cell set ',I4,' prefix: ',A,'\n')
200 Format('Vertex set ',I4)
201 Format('vs',I4.4,'_')
203 Format('\nRegistering vertex set ',I4,' prefix: ',A,'\n')
300 Format('Face set ',I4)
301 Format('fs',I4.4,'_')
303 Format('\nRegistering face set ',I4,' prefix: ',A,'\n')
   End Subroutine MEF90HeatXferCtxSetFromOptions
   
! #undef __FUNCT__
! #define __FUNCT__ "MEF90HeatXferCtxCreateVectors"
! !!!
! !!!  
! !!!  MEF90HeatXferCtxCreateVectors:
! !!!  
! !!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
! !!!

!    Subroutine MEF90HeatXferCtxCreateVectors(MEF90HeatXferCtx,ierr)
!       Type(MEF90HeatXferCtx_Type),Intent(OUT)               :: MEF90HeatXferCtx
!       PetscErrorCode,Intent(OUT)                            :: ierr

!       Allocate(MEF90HeatXferCtx%temperature)
!       PetscCall(DMCreateGlobalVector(MEF90HeatXferCtx%dmScal,MEF90HeatXferCtx%temperature,ierr))
!       PetscCall(PetscObjectSetName(MEF90HeatXferCtx%temperature,"Temperature",ierr))
!       PetscCall(VecSet(MEF90HeatXferCtx%temperature,0.0_Kr,ierr))
!       Allocate(MEF90HeatXferCtx%boundaryTemperature)
!       PetscCall(DMCreateGlobalVector(MEF90HeatXferCtx%dmScal,MEF90HeatXferCtx%boundaryTemperature,ierr))
!       PetscCall(PetscObjectSetName(MEF90HeatXferCtx%boundaryTemperature,"Boundary_Temperature",ierr))
!       PetscCall(VecSet(MEF90HeatXferCtx%boundaryTemperature,0.0_Kr,ierr))
   
!       Allocate(MEF90HeatXferCtx%externalTemperature)
!       PetscCall(DMCreateGlobalVector(MEF90HeatXferCtx%cellDMScal,MEF90HeatXferCtx%externalTemperature,ierr))
!       PetscCall(PetscObjectSetName(MEF90HeatXferCtx%externalTemperature,"External_Temperature",ierr))
!       PetscCall(VecSet(MEF90HeatXferCtx%externalTemperature,0.0_Kr,ierr))

!       Allocate(MEF90HeatXferCtx%flux)
!       PetscCall(DMCreateGlobalVector(MEF90HeatXferCtx%cellDMScal,MEF90HeatXferCtx%flux,ierr))
!       PetscCall(PetscObjectSetName(MEF90HeatXferCtx%flux,"Flux",ierr))
!       PetscCall(VecSet(MEF90HeatXferCtx%flux,0.0_Kr,ierr))
!    End Subroutine MEF90HeatXferCtxCreateVectors   

End Module m_MEF90_HeatXferCtx
