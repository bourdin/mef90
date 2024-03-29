Module m_MEF90_HeatXferCtx_Type
#include "petsc/finclude/petsc.h"
   Use m_MEF90
   Use,Intrinsic :: iso_c_binding
   Implicit none
   Private  
   Public :: MEF90HeatXferCtx_Type
   Public :: MEF90HeatXferGlobalOptions_Type
   Public :: MEF90HeatXferCellSetOptions_Type
   Public :: MEF90HeatXferFaceSetOptions_Type
   Public :: MEF90HeatXferVertexSetOptions_Type
   
   Type MEF90HeatXferCtx_Type
      Type(MEF90Ctx_Type),pointer      :: MEF90Ctx
      Type(tDM)                        :: megaDM
      PetscInt                         :: dim

      Type(tVec),pointer               :: temperatureLocal
      Type(tVec),pointer               :: externalTemperatureLocal
      Type(tVec),pointer               :: fluxLocal
      Type(tVec),pointer               :: boundaryFluxLocal

      Type(tPetscViewer)               :: viewer
      Type(tPetscSF)                   :: temperatureToIOSF,IOToTemperatureSF
      Type(tPetscSF)                   :: boundaryToTemperatureSF
      Type(tPetscSF)                   :: externalTemperatureToIOSF,IOToexternalTemperatureSF
      Type(tPetscSF)                   :: fluxToIOSF,IOTofluxSF
      Type(tPetscSF)                   :: boundaryFluxToIOSF,IOToboundaryFluxSF

      PetscBag                         :: GlobalOptionsBag
      PetscBag,Dimension(:),Pointer    :: CellSetOptionsBag
      PetscBag,Dimension(:),Pointer    :: FaceSetOptionsBag
      PetscBag,Dimension(:),Pointer    :: VertexSetOptionsBag
      PetscBag,Dimension(:),Pointer    :: MaterialPropertiesBag
   End Type MEF90HeatXferCtx_Type
   
   Type MEF90HeatXferGlobalOptions_Type
      PetscEnum                        :: timeSteppingType
      PetscBool                        :: addNullSpace
      PetscReal                        :: initialTemperature
      PetscEnum                        :: boundaryTemperatureScaling
      PetscEnum                        :: externalTemperatureScaling
      PetscEnum                        :: fluxScaling
      PetscEnum                        :: boundaryFluxScaling
      PetscBool                        :: temperatureExport
      !!! scaling = time (step) scaling law currently CST, Linear, Null (not present), File
   End Type MEF90HeatXferGlobalOptions_Type

   Type MEF90HeatXferCellSetOptions_Type
      PetscReal                        :: flux
      PetscBool                        :: Has_BC
      PetscReal                        :: boundaryTemperature
      PetscReal,dimension(3)           :: advectionVector
   End Type MEF90HeatXferCellSetOptions_Type

   Type MEF90HeatXferFaceSetOptions_Type
      PetscReal                        :: boundaryFlux
      PetscReal                        :: surfaceThermalConductivity
      PetscReal                        :: externalTemperature
      PetscBool                        :: Has_BC
      PetscReal                        :: boundaryTemperature
   End Type MEF90HeatXferFaceSetOptions_Type

   Type MEF90HeatXferVertexSetOptions_Type
      PetscBool                        :: Has_BC
      PetscReal                        :: boundaryTemperature
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

Module m_MEF90HeatXferFaceSetOptions_Private
#include "petsc/finclude/petsc.h"
   Use m_MEF90
   Use m_MEF90_HeatXferCtx_Type
   Implicit None

   Private
   Public :: PetscBagGetDataMEF90HeatXferCtxFaceSetOptions
   
   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF90_HeatXferCtx_Type
         PetscBag                                           :: bag
         Type(MEF90HeatXferFaceSetOptions_Type),pointer     :: data
         PetscErrorCode                                     :: ierr
      End subroutine PetscBagGetData
   End interface
Contains

#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90HeatXferCtxFaceSetOptions"
!!!
!!!  PetscBagGetDataMEF90HeatXferCtxFaceSetOptions - Custom interface to PetscGetData
!!!

   Subroutine PetscBagGetDataMEF90HeatXferCtxFaceSetOptions(bag,data,ierr)
      PetscBag                                        :: bag
      Type(MEF90HeatXferFaceSetOptions_Type),pointer  :: data
      PetscErrorCode                                  :: ierr
      
      PetscCall(PetscBagGetData(bag,data,ierr))
   End Subroutine PetscBagGetDataMEF90HeatXferCtxFaceSetOptions
End Module m_MEF90HeatXferFaceSetOptions_Private
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
   Use m_MEF90HeatXferFaceSetOptions_Private
   Use m_MEF90HeatXferVertexSetOptions_Private
   Implicit none

   PetscSizeT,protected   :: sizeofMEF90HeatXferGlobalOptions
   PetscSizeT,protected   :: sizeofMEF90HeatXferCellSetOptions
   PetscSizeT,protected   :: sizeofMEF90HeatXferFaceSetOptions
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
      Type(MEF90HeatXferFaceSetOptions_Type)          :: HeatXferFaceSetOptions
      Type(MEF90HeatXferVertexSetOptions_Type)        :: HeatXferVertexSetOptions
      Character(len=1),pointer                        :: dummychar(:)
      PetscSizeT                                      :: sizeofchar
      
      PetscCall(PetscDataTypeGetSize(PETSC_CHAR,sizeofchar,ierr))
      sizeofMEF90HeatXferGlobalOptions = size(transfer(HeatXferGlobalOptions,dummychar))*sizeofchar
      sizeofMEF90HeatXferCellSetOptions = size(transfer(HeatXferCellSetOptions,dummychar))*sizeofchar
      sizeofMEF90HeatXferFaceSetOptions = size(transfer(HeatXferFaceSetOptions,dummychar))*sizeofchar
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
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
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
      Type(tDM),Dimension(:),Pointer                     :: dmList
      Type(tPetscSF)                                     :: dummySF

      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr))

      PetscCall(MEF90HeatXferCtxInitialize_Private(ierr))
      HeatXferCtx%MEF90Ctx => MEF90Ctx
      PetscCall(PetscBagCreate(MEF90Ctx%comm,sizeofMEF90HeatXferGlobalOptions,HeatXferCtx%GlobalOptionsBag,ierr))
      
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

      PetscCall(DMGetDimension(dm,HeatXferCtx%dim,ierr))
      
      vecName = "Temperature"
      Allocate(HeatXferCtx%temperatureLocal)
      PetscCall(MEF90CreateLocalVector(dm,MEF90GlobalOptions%elementFamily,MEF90GlobalOptions%elementOrder,1_Ki,vecName,HeatXferCtx%temperatureLocal,ierr)) 
      vecName = "ExternalTemperature"
      Allocate(HeatXferCtx%externalTemperatureLocal)
      PetscCall(MEF90CreateBoundaryCellVector(dm,1_Ki,vecName,HeatXferCtx%externalTemperatureLocal,ierr))
      vecName = "Flux"
      Allocate(HeatXferCtx%fluxLocal)
      PetscCall(MEF90CreateCellVector(dm,1_Ki,vecName,HeatXferCtx%fluxLocal,ierr))
      vecName = "BoundaryFlux"
      Allocate(HeatXferCtx%boundaryFluxLocal)
      PetscCall(MEF90CreateBoundaryCellVector(dm,1_Ki,vecName,HeatXferCtx%boundaryFluxLocal,ierr))

      !!! Create the  unknowns and parameters superDM
      Allocate(dmList(4))
      PetscCall(VecGetDM(HeatXferCtx%temperatureLocal,dmList(1),ierr))
      PetscCall(VecGetDM(HeatXferCtx%externalTemperatureLocal,dmList(2),ierr))
      PetscCall(VecGetDM(HeatXferCtx%fluxLocal,dmList(3),ierr))
      PetscCall(VecGetDM(HeatXferCtx%boundaryFluxLocal,dmList(4),ierr))
      PetscCall(DMCreateSuperDM(dmList,4_kI,PETSC_NULL_IS,HeatXferCtx%megaDM,ierr))
      DeAllocate(dmList)

      ! !!! Create the IO SF for all fields
      PetscCall(MEF90IOSFCreate(MEF90Ctx,HeatXferCtx%TemperatureLocal,HeatXferCtx%temperatureToIOSF,HeatXferCtx%IOToTemperatureSF,ierr))
      ! PetscCall(MEF90IOSFCreate(MEF90Ctx,HeatXferCtx%externalTemperatureLocal,HeatXferCtx%externalTemperatureToIOSF,HeatXferCtx%IOToExternalTemperatureSF,ierr))
      PetscCall(MEF90IOSFCreate(MEF90Ctx,HeatXferCtx%fluxLocal,HeatXferCtx%fluxToIOSF,HeatXferCtx%IOToFluxSF,ierr))
      ! PetscCall(MEF90IOSFCreate(MEF90Ctx,HeatXferCtx%boundaryFluxLocal,HeatXferCtx%boundaryFluxToIOSF,HeatXferCtx%IOToBoundaryFluxSF,ierr))

      !!! Create the SF to exchange boundary values of the temperature. 
      PetscCall(MEF90ConstraintSFCreate(HeatXferCtx%MEF90Ctx,HeatXferCtx%TemperatureLocal,HeatXferCtx%temperatureLocal,HeatXferCtx%boundaryToTemperatureSF,dummySF,ierr))
      PetscCall(PetscSFDestroy(dummySF,ierr))
   End Subroutine MEF90HeatXferCtxCreate
   
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferCtxDestroy"
!!!
!!!  
!!!  MEF90HeatXferCtxDestroy:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90HeatXferCtxDestroy(HeatXferCtx,ierr)
      Type(MEF90HeatXferCtx_Type),Intent(INOUT)       :: HeatXferCtx
      PetscErrorCode,Intent(INOUT)                    :: ierr
      
      PetscInt                                        :: set
   
      If (Associated(HeatXferCtx%temperatureLocal)) Then
         PetscCall(VecDestroy(HeatXferCtx%temperatureLocal,ierr))
         DeAllocate(HeatXferCtx%temperatureLocal)
         Nullify(HeatXferCtx%temperatureLocal)
      End If
      If (Associated(HeatXferCtx%ExternalTemperatureLocal)) Then
         PetscCall(VecDestroy(HeatXferCtx%ExternalTemperatureLocal,ierr))
         DeAllocate(HeatXferCtx%ExternalTemperatureLocal)
         Nullify(HeatXferCtx%ExternalTemperatureLocal)
      End If
      If (Associated(HeatXferCtx%fluxLocal)) Then
         PetscCall(VecDestroy(HeatXferCtx%fluxLocal,ierr))
         DeAllocate(HeatXferCtx%fluxLocal)
         Nullify(HeatXferCtx%fluxLocal)
      End If
      If (Associated(HeatXferCtx%boundaryFluxLocal)) Then
         PetscCall(VecDestroy(HeatXferCtx%boundaryFluxLocal,ierr))
         DeAllocate(HeatXferCtx%boundaryFluxLocal)
         Nullify(HeatXferCtx%boundaryFluxLocal)
      End If

      !!! Destroy SFs
      PetscCall(PetscSFDestroy(HeatXferCtx%temperatureToIOSF,ierr))
      PetscCall(PetscSFDestroy(HeatXferCtx%IOToTemperatureSF,ierr))
      PetscCall(PetscSFDestroy(HeatXferCtx%externalTemperatureToIOSF,ierr))
      PetscCall(PetscSFDestroy(HeatXferCtx%IOToExternalTemperatureSF,ierr))
      PetscCall(PetscSFDestroy(HeatXferCtx%fluxToIOSF,ierr))
      PetscCall(PetscSFDestroy(HeatXferCtx%IOToFluxSF,ierr))
      PetscCall(PetscSFDestroy(HeatXferCtx%boundaryFluxToIOSF,ierr))
      PetscCall(PetscSFDestroy(HeatXferCtx%IOToBoundaryFluxSF,ierr))

      PetscCall(PetscSFDestroy(HeatXferCtx%boundaryToTemperatureSF,ierr))

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

      Do set = 1, size(HeatXferCtx%MaterialPropertiesBag)
         PetscCall(PetscBagDestroy(HeatXferCtx%MaterialPropertiesBag(set),ierr))
      End Do
      DeAllocate(HeatXferCtx%MaterialPropertiesBag)

      PetscCall(DMDestroy(HeatXferCtx%megaDM,ierr))
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
      PetscErrorCode,Intent(INOUT)                       :: ierr

      Type(MEF90HeatXferGlobalOptions_Type),pointer      :: HeatXferGlobalOptions
      PetscCall(PetscBagGetDataMEF90HeatXferCtxGlobalOptions(bag,HeatXferGlobalOptions,ierr))
      PetscCall(PetscBagSetName(bag,trim(name),"HeatXferGlobalOptions MEF90 Heat transfer global options",ierr))
      PetscCall(PetscBagSetOptionsPrefix(bag,trim(prefix),ierr))

      PetscCall(PetscBagRegisterBool(bag,HeatXferGlobalOptions%temperatureExport,default%temperatureExport,'temperature_export','Export temperature',ierr))

      PetscCall(PetscBagRegisterEnum(bag,HeatXferGlobalOptions%timeSteppingType,MEF90HeatXFer_timeSteppingTypeList,default%timeSteppingType,'heatxfer_timeStepping_type','Type of heat transfer computation',ierr))
      PetscCall(PetscBagRegisterBool(bag,HeatXferGlobalOptions%addNullSpace,default%addNullSpace,'heatxfer_addNullSpace','Add null space to SNES',ierr))
      PetscCall(PetscBagRegisterReal(bag,HeatXferGlobalOptions%initialTemperature,default%initialTemperature,'heatxfer_initialTemperature','[K] (T): Initial Temperature' ,ierr))

      PetscCall(PetscBagRegisterEnum(bag,HeatXferGlobalOptions%boundaryTemperatureScaling,MEF90ScalingList,default%boundaryTemperatureScaling,'boundaryTemperature_scaling','Boundary temperature scaling',ierr))

      PetscCall(PetscBagRegisterEnum(bag,HeatXferGlobalOptions%externalTemperatureScaling,MEF90ScalingList,default%externalTemperatureScaling,'externalTemperature_scaling','External Temperature scaling',ierr))

      PetscCall(PetscBagRegisterEnum(bag,HeatXferGlobalOptions%fluxScaling,MEF90ScalingList,default%fluxScaling,'flux_scaling','Heat flux scaling',ierr))
      PetscCall(PetscBagRegisterEnum(bag,HeatXferGlobalOptions%boundaryFluxScaling,MEF90ScalingList,default%boundaryFluxScaling,'boundaryFlux_scaling','Boundary heat flux scaling',ierr))
   End Subroutine PetscBagRegisterMEF90HeatXferCtxGlobalOptions

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90HeatXferCtxCellSetOptions"
!!!
!!!  
!!!  PetscBagRegisterMEF90HeatXferCtxCellSetOptions:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine PetscBagRegisterMEF90HeatXferCtxCellSetOptions(bag,name,prefix,default,ierr)
      PetscBag                                           :: bag
      Character(len=*),Intent(IN)                        :: prefix,name
      Type(MEF90HeatXferCellSetOptions_Type),Intent(IN)  :: default
      PetscErrorCode,Intent(INOUT)                       :: ierr

      Type(MEF90HeatXferCellSetOptions_Type),pointer      :: HeatXferCellSetOptions
      PetscCall(PetscBagGetDataMEF90HeatXferCtxCellSetOptions(bag,HeatXferCellSetOptions,ierr))
      PetscCall(PetscBagSetName(bag,trim(name),"HeatXferCellSetOptions MEF90 Heat transfer Cell Set options",ierr))
      PetscCall(PetscBagSetOptionsPrefix(bag,trim(prefix),ierr))

      HeatXferCellSetOptions%advectionVector = default%advectionVector
      PetscCall(PetscBagRegisterReal(bag,HeatXferCellSetOptions%Flux,default%Flux,'Flux','[J.s^(-1).m^(-3) / J.s^(-1).m^(-2)] (f): Internal / applied heat flux',ierr))
      PetscCall(PetscBagRegisterBool(bag,HeatXferCellSetOptions%Has_BC,default%Has_BC,'TemperatureBC','Temperature has Dirichlet boundary Condition (Y/N)',ierr))
      PetscCall(PetscBagRegisterReal(bag,HeatXferCellSetOptions%boundaryTemperature,default%boundaryTemperature,'boundaryTemperature','Temperature boundary value',ierr))
      PetscCall(PetscBagRegisterRealArray(bag,HeatXferCellSetOptions%advectionVector,3_Ki,'advectionVector','[m.s^(-1)] (V): advection vector',ierr))
   End Subroutine PetscBagRegisterMEF90HeatXferCtxCellSetOptions

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90HeatXferCtxFaceSetOptions"
!!!
!!!  
!!!  PetscBagRegisterMEF90HeatXferCtxFaceSetOptions:
!!!  
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine PetscBagRegisterMEF90HeatXferCtxFaceSetOptions(bag,name,prefix,default,ierr)
      PetscBag                                           :: bag
      Character(len=*),Intent(IN)                        :: prefix,name
      Type(MEF90HeatXferFaceSetOptions_Type),Intent(IN)  :: default
      PetscErrorCode,Intent(INOUT)                       :: ierr

      Type(MEF90HeatXferFaceSetOptions_Type),pointer      :: HeatXferFaceSetOptions
      PetscCall(PetscBagGetDataMEF90HeatXferCtxFaceSetOptions(bag,HeatXferFaceSetOptions,ierr))
      PetscCall(PetscBagSetName(bag,trim(name),"HeatXferFaceSetOptions MEF90 Heat transfer Face Set options",ierr))
      PetscCall(PetscBagSetOptionsPrefix(bag,trim(prefix),ierr))

      PetscCall(PetscBagRegisterReal(bag,HeatXferFaceSetOptions%boundaryFlux,default%boundaryFlux,'boundaryFlux','[J.s^(-1).m^(-2) / J.s^(-1).m^(-1)] (f): Internal / applied heat flux',ierr))
      PetscCall(PetscBagRegisterReal(bag,HeatXferFaceSetOptions%SurfaceThermalConductivity,default%SurfaceThermalConductivity,'SurfaceThermalConductivity','[J.s^(-1).m^(-2).K^(-1) / J.s^(-1).m^(-1).K^(-1) ] (H) Surface Thermal Conductivity',ierr))
      PetscCall(PetscBagRegisterReal(bag,HeatXferFaceSetOptions%externalTemperature,default%externalTemperature,'externalTemperature','Reference temperature T [K]',ierr))
      PetscCall(PetscBagRegisterBool(bag,HeatXferFaceSetOptions%Has_BC,default%Has_BC,'TemperatureBC','Temperature has Dirichlet boundary Condition (Y/N)',ierr))
      PetscCall(PetscBagRegisterReal(bag,HeatXferFaceSetOptions%boundaryTemperature,default%boundaryTemperature,'boundaryTemperature','Temperature boundary value',ierr))
   End Subroutine PetscBagRegisterMEF90HeatXferCtxFaceSetOptions
   
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
      PetscErrorCode,Intent(INOUT)                          :: ierr

      Type(MEF90HeatXferVertexSetOptions_Type),pointer      :: HeatXferVertexSetOptions
      PetscCall(PetscBagGetDataMEF90HeatXferCtxVertexSetOptions(bag,HeatXferVertexSetOptions,ierr))
      PetscCall(PetscBagSetName(bag,trim(name),"HeatXferVertexSetOptions MEF90 Heat transfer Vertex Set options",ierr))
      PetscCall(PetscBagSetOptionsPrefix(bag,trim(prefix),ierr))

      PetscCall(PetscBagRegisterBool(bag,HeatXferVertexSetOptions%Has_BC,default%Has_BC,'TemperatureBC','Temperature has Dirichlet boundary Condition (Y/N)',ierr))
      PetscCall(PetscBagRegisterReal(bag,HeatXferVertexSetOptions%boundaryTemperature,default%boundaryTemperature,'boundaryTemperature','Temperature boundary value',ierr))
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
                                              defaultCellSetOptions,defaultFaceSetOptions,defaultVertexSetOptions,ierr)
      Type(MEF90HeatXferCtx_Type),Intent(INOUT)             :: heatXferCtx
      Character(len=*),Intent(IN)                           :: prefix
      Type(MEF90HeatXferGlobalOptions_Type),Intent(IN)      :: defaultGlobalOptions
      Type(MEF90HeatXferCellSetOptions_Type),Intent(IN)     :: defaultCellSetOptions
      Type(MEF90HeatXferFaceSetOptions_Type),Intent(IN)     :: defaultFaceSetOptions
      Type(MEF90HeatXferVertexSetOptions_Type),Intent(IN)   :: defaultVertexSetOptions
      PetscErrorCode,Intent(INOUT)                          :: ierr
   
      Type(MEF90CtxGlobalOptions_Type),pointer              :: MEF90CtxGlobalOptions
      Type(tIS)                                             :: setIS
      PetscInt,Dimension(:),Pointer                         :: setID
      PetscInt                                              :: set
      Character(len=MEF90MXSTRLEN)                          :: IOBuffer,setName,setprefix

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
      !!!
      PetscCall(DMGetLabelIdIS(heatXferCtx%megaDM,MEF90CellSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(heatXferCtx%MEF90Ctx%comm,setIS,ierr)) 
      PetscCall(ISGetIndicesF90(setIS,setID,ierr))
      Do set = 1, size(setID)
         Write(setName,"('Cell set ',I4)") setID(set)
         Write(setprefix,"('cs',I4.4,'_')") setID(set)
         PetscCall(PetscBagRegisterMEF90HeatXferCtxCellSetOptions(heatXferCtx%CellSetOptionsBag(set),setName,setPrefix,defaultCellSetOptions,ierr))
         If (MEF90CtxGlobalOptions%verbose > 0) Then
            Write(IOBuffer,"('\nRegistering cell set ',I4,' prefix: ',A,'\n')") setID(set),trim(setprefix)
            PetscCall(PetscPrintf(heatXferCtx%MEF90Ctx%comm,IOBuffer,ierr))
            PetscCall(PetscBagView(heatXferCtx%CellSetOptionsBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr))
            PetscCall(PetscPrintf(heatXferCtx%MEF90Ctx%comm,"\n",ierr))
         End if
      End Do
      PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      PetscCall(ISDestroy(setIS,ierr))
      
      !!!
      !!! Registering Face Set Context
      !!!
      PetscCall(DMGetLabelIdIS(heatXferCtx%megaDM,MEF90FaceSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(heatXferCtx%MEF90Ctx%comm,setIS,ierr)) 
      PetscCall(ISGetIndicesF90(setIS,setID,ierr))
      Do set = 1, size(setID)
         Write(setName,"('Face set ',I4)") setID(set)
         Write(setprefix,"('fs',I4.4,'_')") setID(set)
         PetscCall(PetscBagRegisterMEF90HeatXferCtxFaceSetOptions(heatXferCtx%FaceSetOptionsBag(set),setName,setPrefix,defaultFaceSetOptions,ierr))
         If (MEF90CtxGlobalOptions%verbose > 0) Then
            Write(IOBuffer,"('\nRegistering face set ',I4,' prefix: ',A,'\n')") setID(set),trim(setprefix)
            PetscCall(PetscPrintf(heatXferCtx%MEF90Ctx%comm,IOBuffer,ierr))
            PetscCall(PetscBagView(heatXferCtx%FaceSetOptionsBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr))
            PetscCall(PetscPrintf(heatXferCtx%MEF90Ctx%comm,"\n",ierr))
         End if
      End Do
      PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      PetscCall(ISDestroy(setIS,ierr))
      
      !!!
      !!! Registering Vertex Set Context
      !!!
      PetscCall(DMGetLabelIdIS(heatXferCtx%megaDM,MEF90VertexSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr)) 
      PetscCall(ISGetIndicesF90(setIS,setID,ierr))
      Do set = 1, size(setID)
         Write(setName,"('Vertex set ',I4)") setID(set)
         Write(setprefix,"('vs',I4.4,'_')") setID(set)
         PetscCall(PetscBagRegisterMEF90HeatXferCtxVertexSetOptions(heatXferCtx%VertexSetOptionsBag(set),setName,setPrefix,defaultVertexSetOptions,ierr))
         If (MEF90CtxGlobalOptions%verbose > 0) Then
            Write(IOBuffer,"('\nRegistering vertex set ',I4,' prefix: ',A,'\n')") setID(set),trim(setprefix)
            PetscCall(PetscPrintf(heatXferCtx%MEF90Ctx%comm,IOBuffer,ierr))
            PetscCall(PetscBagView(heatXferCtx%VertexSetOptionsBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr))
            PetscCall(PetscPrintf(heatXferCtx%MEF90Ctx%comm,"\n",ierr))
         End if
      End Do
      PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      PetscCall(ISDestroy(setIS,ierr))
   End Subroutine MEF90HeatXferCtxSetFromOptions
End Module m_MEF90_HeatXferCtx
