module m_MEF90_HeatXferCtx_Type
#include "petsc/finclude/petsc.h"
   use m_MEF90_Ctx
   implicit none(type, external)
   private
   public :: MEF90HeatXferCtx_Type
   public :: MEF90HeatXferGlobalOptions_Type
   public :: MEF90HeatXferCellSetOptions_Type
   public :: MEF90HeatXferFaceSetOptions_Type
   public :: MEF90HeatXferVertexSetOptions_Type

   type MEF90HeatXferCtx_Type
      type(MEF90Ctx_Type), pointer          :: MEF90Ctx
      type(tDM)                            :: megaDM
      PetscInt                             :: dim

      type(tVec), pointer                   :: temperatureLocal
      type(tVec), pointer                   :: externalTemperatureLocal
      type(tVec), pointer                   :: fluxLocal
      type(tVec), pointer                   :: boundaryFluxLocal

      type(tPetscViewer)                   :: viewer
      type(tPetscSF)                       :: temperatureToIOSF, IOToTemperatureSF
      type(tPetscSF)                       :: boundaryToTemperatureSF
      type(tPetscSF)                       :: externalTemperatureToIOSF, IOToexternalTemperatureSF
      type(tPetscSF)                       :: fluxToIOSF, IOTofluxSF
      type(tPetscSF)                       :: boundaryFluxToIOSF, IOToboundaryFluxSF

      type(tPetscBag)                      :: GlobalOptionsBag
      type(tPetscBag), dimension(:), pointer :: CellSetOptionsBag
      type(tPetscBag), dimension(:), pointer :: FaceSetOptionsBag
      type(tPetscBag), dimension(:), pointer :: VertexSetOptionsBag
      type(tPetscBag), dimension(:), pointer :: MaterialPropertiesBag
   end type MEF90HeatXferCtx_Type

   type MEF90HeatXferGlobalOptions_Type
      PetscEnum                        :: timeSteppingType
      PetscBool                        :: addNullSpace
      PetscReal                        :: initialTemperature
      PetscEnum                        :: boundaryTemperatureScaling
      PetscEnum                        :: externalTemperatureScaling
      PetscEnum                        :: fluxScaling
      PetscEnum                        :: boundaryFluxScaling
      PetscBool                        :: temperatureExport
      !!! scaling = time (step) scaling law currently CST, Linear, Null (not present), File
   end type MEF90HeatXferGlobalOptions_Type

   type MEF90HeatXferCellSetOptions_Type
      PetscReal                        :: flux
      PetscBool                        :: Has_BC
      PetscReal                        :: boundaryTemperature
      PetscReal, dimension(3)           :: advectionVector
   end type MEF90HeatXferCellSetOptions_Type

   type MEF90HeatXferFaceSetOptions_Type
      PetscReal                        :: boundaryFlux
      PetscReal                        :: surfaceThermalConductivity
      PetscReal                        :: externalTemperature
      PetscBool                        :: Has_BC
      PetscReal                        :: boundaryTemperature
   end type MEF90HeatXferFaceSetOptions_Type

   type MEF90HeatXferVertexSetOptions_Type
      PetscBool                        :: Has_BC
      PetscReal                        :: boundaryTemperature
   end type MEF90HeatXferVertexSetOptions_Type
end module m_MEF90_HeatXferCtx_Type

module m_MEF90HeatXferGlobalOptions_Private
#include "petsc/finclude/petsc.h"
   use petscbag
   use m_MEF90_HeatXferCtx_Type
   implicit none(type, external)

   private
   public :: PetscBagGetDataMEF90HeatXferCtxGlobalOptions

contains

#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90HeatXferCtxGlobalOptions"
!!!
!!!  PetscBagGetDataMEF90HeatXferCtxGlobalOptions - Custom interface to PetscGetData
!!!

   subroutine PetscBagGetDataMEF90HeatXferCtxGlobalOptions(bag, data, ierr)
      type(tPetscBag)                                    :: bag
      type(MEF90HeatXferGlobalOptions_Type), pointer   :: data
      PetscErrorCode                                  :: ierr

      PetscCall(PetscBagGetData(bag, data, ierr))
   end subroutine PetscBagGetDataMEF90HeatXferCtxGlobalOptions
end module m_MEF90HeatXferGlobalOptions_Private

module m_MEF90HeatXferCellSetOptions_Private
#include "petsc/finclude/petsc.h"
   use petscbag
   use m_MEF90_HeatXferCtx_Type
   implicit none(type, external)

   private
   public :: PetscBagGetDataMEF90HeatXferCtxCellSetOptions

contains

#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90HeatXferCtxCellSetOptions"
!!!
!!!  PetscBagGetDataMEF90HeatXferCtxCellSetOptions - Custom interface to PetscGetData
!!!

   subroutine PetscBagGetDataMEF90HeatXferCtxCellSetOptions(bag, data, ierr)
      type(tPetscBag)                                    :: bag
      type(MEF90HeatXferCellSetOptions_Type), pointer  :: data
      PetscErrorCode                                  :: ierr

      PetscCall(PetscBagGetData(bag, data, ierr))
   end subroutine PetscBagGetDataMEF90HeatXferCtxCellSetOptions
end module m_MEF90HeatXferCellSetOptions_Private

module m_MEF90HeatXferFaceSetOptions_Private
#include "petsc/finclude/petsc.h"
   use petscbag
   use m_MEF90_HeatXferCtx_Type
   implicit none(type, external)

   private
   public :: PetscBagGetDataMEF90HeatXferCtxFaceSetOptions

contains

#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90HeatXferCtxFaceSetOptions"
!!!
!!!  PetscBagGetDataMEF90HeatXferCtxFaceSetOptions - Custom interface to PetscGetData
!!!

   subroutine PetscBagGetDataMEF90HeatXferCtxFaceSetOptions(bag, data, ierr)
      type(tPetscBag)                                 :: bag
      type(MEF90HeatXferFaceSetOptions_Type), pointer  :: data
      PetscErrorCode                                  :: ierr

      PetscCall(PetscBagGetData(bag, data, ierr))
   end subroutine PetscBagGetDataMEF90HeatXferCtxFaceSetOptions
end module m_MEF90HeatXferFaceSetOptions_Private

module m_MEF90HeatXferVertexSetOptions_Private
#include "petsc/finclude/petsc.h"
   use petscbag
   use m_MEF90_HeatXferCtx_Type
   implicit none(type, external)

   private
   public :: PetscBagGetDataMEF90HeatXferCtxVertexSetOptions

contains
#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90HeatXferCtxVertexSetOptions"
!!!
!!!  PetscBagGetDataMEF90HeatXferCtxVertexSetOptions - Custom interface to PetscGetData
!!!

   subroutine PetscBagGetDataMEF90HeatXferCtxVertexSetOptions(bag, data, ierr)
      type(tPetscBag)                                    :: bag
      type(MEF90HeatXferVertexSetOptions_Type), pointer   :: data
      PetscErrorCode                                     :: ierr

      PetscCall(PetscBagGetData(bag, data, ierr))
   end subroutine PetscBagGetDataMEF90HeatXferCtxVertexSetOptions
end module m_MEF90HeatXferVertexSetOptions_Private

module m_MEF90_HeatXferCtx
#include "petsc/finclude/petsc.h"
   use m_MEF90_Ctx
   use m_MEF90_DMPlex
   use m_MEF90_Materials
   use m_MEF90_HeatXferCtx_Type
   use m_MEF90HeatXferGlobalOptions_Private
   use m_MEF90HeatXferCellSetOptions_Private
   use m_MEF90HeatXferFaceSetOptions_Private
   use m_MEF90HeatXferVertexSetOptions_Private
   implicit none(type)

   PetscSizeT, protected   :: sizeofMEF90HeatXferGlobalOptions
   PetscSizeT, protected   :: sizeofMEF90HeatXferCellSetOptions
   PetscSizeT, protected   :: sizeofMEF90HeatXferFaceSetOptions
   PetscSizeT, protected   :: sizeofMEF90HeatXferVertexSetOptions

   enum, bind(c)
      enumerator  :: MEF90HeatXfer_timeSteppingTypeNULL = 0, &
         MEF90HeatXFer_timeSteppingTypeSteadyState, &
         MEF90HeatXFer_timeSteppingTypeTransient
   end enum
   character(len=MEF90MXSTRLEN), dimension(6), protected   :: MEF90HeatXFer_timeSteppingTypeList

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferCtxInitialize_Private"
!!!
!!!
!!!  MEF90HeatXferCtxInitialize_Private:
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   subroutine MEF90HeatXferCtxInitialize_Private(ierr)
      PetscErrorCode, intent(OUT)                      :: ierr

      type(MEF90HeatXferGlobalOptions_Type)           :: HeatXferGlobalOptions
      type(MEF90HeatXferCellSetOptions_Type)          :: HeatXferCellSetOptions
      type(MEF90HeatXferFaceSetOptions_Type)          :: HeatXferFaceSetOptions
      type(MEF90HeatXferVertexSetOptions_Type)        :: HeatXferVertexSetOptions
      character(len=1), pointer                        :: dummychar(:)
      PetscSizeT                                      :: sizeofchar

      PetscCall(PetscDataTypeGetSize(PETSC_CHAR, sizeofchar, ierr))
      sizeofMEF90HeatXferGlobalOptions = size(transfer(HeatXferGlobalOptions, dummychar)) * sizeofchar
      sizeofMEF90HeatXferCellSetOptions = size(transfer(HeatXferCellSetOptions, dummychar)) * sizeofchar
      sizeofMEF90HeatXferFaceSetOptions = size(transfer(HeatXferFaceSetOptions, dummychar)) * sizeofchar
      sizeofMEF90HeatXferVertexSetOptions = size(transfer(HeatXferVertexSetOptions, dummychar)) * sizeofchar

      MEF90HeatXFer_timeSteppingTypeList(1) = 'null'
      MEF90HeatXFer_timeSteppingTypeList(2) = 'SteadyState'
      MEF90HeatXFer_timeSteppingTypeList(3) = 'Transient'
      MEF90HeatXFer_timeSteppingTypeList(4) = 'MEF90_HeatXFer_timeSteppingType'
      MEF90HeatXFer_timeSteppingTypeList(5) = '_MEF90_HeatXFer_timeSteppingType'
      MEF90HeatXFer_timeSteppingTypeList(6) = ''
   end subroutine MEF90HeatXferCtxInitialize_Private

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferCtxCreate"
!!!
!!!
!!!  MEF90HeatXferCtxCreate:
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90HeatXferCtxCreate(HeatXferCtx, dm, MEF90Ctx, ierr)
      type(MEF90HeatXferCtx_Type), intent(OUT)            :: HeatXferCtx
      type(tDM), target, intent(IN)                        :: dm
      type(MEF90Ctx_Type), target, intent(IN)              :: MEF90Ctx
      PetscErrorCode, intent(INOUT)                       :: ierr

      type(MEF90CtxGlobalOptions_Type), pointer           :: MEF90GlobalOptions
      type(tIS)                                          :: setIS
      PetscInt                                           :: set, numSet
      character(len=MEF90MXSTRLEN)                       :: vecName
      type(tDM), dimension(:), pointer                     :: dmList
      type(tPetscSF)                                     :: dummySF

      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag, MEF90GlobalOptions, ierr))

      PetscCall(MEF90HeatXferCtxInitialize_Private(ierr))
      HeatXferCtx%MEF90Ctx => MEF90Ctx
      PetscCall(PetscBagCreate(MEF90Ctx%comm, sizeofMEF90HeatXferGlobalOptions, HeatXferCtx%GlobalOptionsBag, ierr))

      !!! I need to allocate for the overall number of sets, not the local one
      PetscCall(DMGetLabelIdIS(dm, MEF90CellSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(PETSC_COMM_WORLD, setIS, ierr))
      PetscCall(ISGetLocalSize(setIS, numSet, ierr))
      allocate (HeatXferCtx%CellSetOptionsBag(numSet), stat=ierr)
      do set = 1, numSet
         PetscCall(PetscBagCreate(MEF90Ctx%comm, sizeofMEF90HeatXferCellSetOptions, HeatXferCtx%CellSetOptionsBag(set), ierr))
      end do
      PetscCall(ISDestroy(setIS, ierr))

      PetscCall(DMGetLabelIdIS(dm, MEF90FaceSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(PETSC_COMM_WORLD, setIS, ierr))
      PetscCall(ISGetLocalSize(setIS, numSet, ierr))
      allocate (HeatXferCtx%FaceSetOptionsBag(numSet), stat=ierr)
      do set = 1, numSet
         PetscCall(PetscBagCreate(MEF90Ctx%comm, sizeofMEF90HeatXferCellSetOptions, HeatXferCtx%FaceSetOptionsBag(set), ierr))
      end do
      PetscCall(ISDestroy(setIS, ierr))

      PetscCall(DMGetLabelIdIS(dm, MEF90VertexSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(PETSC_COMM_WORLD, setIS, ierr))
      PetscCall(ISGetLocalSize(setIS, numSet, ierr))
      allocate (HeatXferCtx%VertexSetOptionsBag(numSet), stat=ierr)
      do set = 1, numSet
         PetscCall(PetscBagCreate(MEF90Ctx%comm, sizeofMEF90HeatXferVertexSetOptions, HeatXferCtx%VertexSetOptionsBag(set), ierr))
      end do
      PetscCall(ISDestroy(setIS, ierr))

      PetscCall(DMGetDimension(dm, HeatXferCtx%dim, ierr))

      vecName = "Temperature"
      allocate (HeatXferCtx%temperatureLocal)
      PetscCall(MEF90CreateLocalVector(dm, MEF90GlobalOptions%elementFamily, MEF90GlobalOptions%elementOrder, 1_ki, vecName, HeatXferCtx%temperatureLocal, ierr))
      vecName = "ExternalTemperature"
      allocate (HeatXferCtx%externalTemperatureLocal)
      PetscCall(MEF90CreateBoundaryCellVector(dm, 1_ki, vecName, HeatXferCtx%externalTemperatureLocal, ierr))
      vecName = "Flux"
      allocate (HeatXferCtx%fluxLocal)
      PetscCall(MEF90CreateCellVector(dm, 1_ki, vecName, HeatXferCtx%fluxLocal, ierr))
      vecName = "BoundaryFlux"
      allocate (HeatXferCtx%boundaryFluxLocal)
      PetscCall(MEF90CreateBoundaryCellVector(dm, 1_ki, vecName, HeatXferCtx%boundaryFluxLocal, ierr))

      !!! Create the  unknowns and parameters superDM
      allocate (dmList(4))
      PetscCall(VecGetDM(HeatXferCtx%temperatureLocal, dmList(1), ierr))
      PetscCall(VecGetDM(HeatXferCtx%externalTemperatureLocal, dmList(2), ierr))
      PetscCall(VecGetDM(HeatXferCtx%fluxLocal, dmList(3), ierr))
      PetscCall(VecGetDM(HeatXferCtx%boundaryFluxLocal, dmList(4), ierr))
      PetscCall(DMCreateSuperDM(dmList, 4_ki, PETSC_NULL_IS_POINTER, HeatXferCtx%megaDM, ierr))
      deallocate (dmList)

      ! !!! Create the IO SF for all fields
      PetscCall(MEF90IOSFCreate(MEF90Ctx, HeatXferCtx%TemperatureLocal, HeatXferCtx%temperatureToIOSF, HeatXferCtx%IOToTemperatureSF, ierr))
      ! PetscCall(MEF90IOSFCreate(MEF90Ctx,HeatXferCtx%externalTemperatureLocal,HeatXferCtx%externalTemperatureToIOSF,HeatXferCtx%IOToExternalTemperatureSF,ierr))
      PetscCall(MEF90IOSFCreate(MEF90Ctx, HeatXferCtx%fluxLocal, HeatXferCtx%fluxToIOSF, HeatXferCtx%IOToFluxSF, ierr))
      ! PetscCall(MEF90IOSFCreate(MEF90Ctx,HeatXferCtx%boundaryFluxLocal,HeatXferCtx%boundaryFluxToIOSF,HeatXferCtx%IOToBoundaryFluxSF,ierr))

      !!! Create the SF to exchange boundary values of the temperature.
      PetscCall(MEF90ConstraintSFCreate(HeatXferCtx%MEF90Ctx, HeatXferCtx%TemperatureLocal, HeatXferCtx%temperatureLocal, HeatXferCtx%boundaryToTemperatureSF, dummySF, ierr))
      PetscCall(PetscSFDestroy(dummySF, ierr))
   end subroutine MEF90HeatXferCtxCreate

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferCtxDestroy"
!!!
!!!
!!!  MEF90HeatXferCtxDestroy:
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90HeatXferCtxDestroy(HeatXferCtx, ierr)
      type(MEF90HeatXferCtx_Type), intent(INOUT)       :: HeatXferCtx
      PetscErrorCode, intent(INOUT)                    :: ierr

      PetscInt                                        :: set

      if (associated(HeatXferCtx%temperatureLocal)) then
         PetscCall(VecDestroy(HeatXferCtx%temperatureLocal, ierr))
         deallocate (HeatXferCtx%temperatureLocal)
         nullify (HeatXferCtx%temperatureLocal)
      end if
      if (associated(HeatXferCtx%ExternalTemperatureLocal)) then
         PetscCall(VecDestroy(HeatXferCtx%ExternalTemperatureLocal, ierr))
         deallocate (HeatXferCtx%ExternalTemperatureLocal)
         nullify (HeatXferCtx%ExternalTemperatureLocal)
      end if
      if (associated(HeatXferCtx%fluxLocal)) then
         PetscCall(VecDestroy(HeatXferCtx%fluxLocal, ierr))
         deallocate (HeatXferCtx%fluxLocal)
         nullify (HeatXferCtx%fluxLocal)
      end if
      if (associated(HeatXferCtx%boundaryFluxLocal)) then
         PetscCall(VecDestroy(HeatXferCtx%boundaryFluxLocal, ierr))
         deallocate (HeatXferCtx%boundaryFluxLocal)
         nullify (HeatXferCtx%boundaryFluxLocal)
      end if

      !!! Destroy SFs
      PetscCall(PetscSFDestroy(HeatXferCtx%temperatureToIOSF, ierr))
      PetscCall(PetscSFDestroy(HeatXferCtx%IOToTemperatureSF, ierr))
      PetscCall(PetscSFDestroy(HeatXferCtx%externalTemperatureToIOSF, ierr))
      PetscCall(PetscSFDestroy(HeatXferCtx%IOToExternalTemperatureSF, ierr))
      PetscCall(PetscSFDestroy(HeatXferCtx%fluxToIOSF, ierr))
      PetscCall(PetscSFDestroy(HeatXferCtx%IOToFluxSF, ierr))
      PetscCall(PetscSFDestroy(HeatXferCtx%boundaryFluxToIOSF, ierr))
      PetscCall(PetscSFDestroy(HeatXferCtx%IOToBoundaryFluxSF, ierr))

      PetscCall(PetscSFDestroy(HeatXferCtx%boundaryToTemperatureSF, ierr))

      nullify (HeatXferCtx%MEF90Ctx)

      PetscCall(PetscBagDestroy(HeatXferCtx%GlobalOptionsBag, ierr))
      do set = 1, size(HeatXferCtx%CellSetOptionsBag)
         PetscCall(PetscBagDestroy(HeatXferCtx%CellSetOptionsBag(set), ierr))
      end do
      deallocate (HeatXferCtx%CellSetOptionsBag)
      do set = 1, size(HeatXferCtx%FaceSetOptionsBag)
         PetscCall(PetscBagDestroy(HeatXferCtx%FaceSetOptionsBag(set), ierr))
      end do
      deallocate (HeatXferCtx%FaceSetOptionsBag)

      do set = 1, size(HeatXferCtx%VertexSetOptionsBag)
         PetscCall(PetscBagDestroy(HeatXferCtx%VertexSetOptionsBag(set), ierr))
      end do
      deallocate (HeatXferCtx%VertexSetOptionsBag)

      do set = 1, size(HeatXferCtx%MaterialPropertiesBag)
         PetscCall(PetscBagDestroy(HeatXferCtx%MaterialPropertiesBag(set), ierr))
      end do
      deallocate (HeatXferCtx%MaterialPropertiesBag)

      PetscCall(DMDestroy(HeatXferCtx%megaDM, ierr))
   end subroutine MEF90HeatXferCtxDestroy

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90HeatXferCtxGlobalOptions"
!!!
!!!
!!!  PetscBagRegisterMEF90HeatXferCtxGlobalOptions:
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   subroutine PetscBagRegisterMEF90HeatXferCtxGlobalOptions(bag, name, prefix, default, ierr)
      type(tPetscBag)                                    :: bag
      character(len=*), intent(IN)                        :: prefix, name
      type(MEF90HeatXferGlobalOptions_Type), intent(IN)   :: default
      PetscErrorCode, intent(INOUT)                       :: ierr

      type(MEF90HeatXferGlobalOptions_Type), pointer      :: HeatXferGlobalOptions
      PetscCall(PetscBagGetDataMEF90HeatXferCtxGlobalOptions(bag, HeatXferGlobalOptions, ierr))
      PetscCall(PetscBagSetName(bag, trim(name), "HeatXferGlobalOptions MEF90 Heat transfer global options", ierr))
      PetscCall(PetscBagSetOptionsPrefix(bag, trim(prefix), ierr))

      PetscCall(PetscBagRegisterBool(bag, HeatXferGlobalOptions%temperatureExport, default%temperatureExport, 'temperature_export', 'Export temperature', ierr))

      PetscCall(PetscBagRegisterEnum(bag, HeatXferGlobalOptions%timeSteppingType, MEF90HeatXFer_timeSteppingTypeList, default%timeSteppingType, 'heatxfer_timeStepping_type', 'Type of heat transfer computation', ierr))
      PetscCall(PetscBagRegisterBool(bag, HeatXferGlobalOptions%addNullSpace, default%addNullSpace, 'heatxfer_addNullSpace', 'Add null space to SNES', ierr))
      PetscCall(PetscBagRegisterReal(bag, HeatXferGlobalOptions%initialTemperature, default%initialTemperature, 'heatxfer_initialTemperature', '[K] (T): Initial Temperature', ierr))

      PetscCall(PetscBagRegisterEnum(bag, HeatXferGlobalOptions%boundaryTemperatureScaling, MEF90ScalingList, default%boundaryTemperatureScaling, 'boundaryTemperature_scaling', 'Boundary temperature scaling', ierr))

      PetscCall(PetscBagRegisterEnum(bag, HeatXferGlobalOptions%externalTemperatureScaling, MEF90ScalingList, default%externalTemperatureScaling, 'externalTemperature_scaling', 'External Temperature scaling', ierr))

      PetscCall(PetscBagRegisterEnum(bag, HeatXferGlobalOptions%fluxScaling, MEF90ScalingList, default%fluxScaling, 'flux_scaling', 'Heat flux scaling', ierr))
      PetscCall(PetscBagRegisterEnum(bag, HeatXferGlobalOptions%boundaryFluxScaling, MEF90ScalingList, default%boundaryFluxScaling, 'boundaryFlux_scaling', 'Boundary heat flux scaling', ierr))
   end subroutine PetscBagRegisterMEF90HeatXferCtxGlobalOptions

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90HeatXferCtxCellSetOptions"
!!!
!!!
!!!  PetscBagRegisterMEF90HeatXferCtxCellSetOptions:
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine PetscBagRegisterMEF90HeatXferCtxCellSetOptions(bag, name, prefix, default, ierr)
      type(tPetscBag)                                    :: bag
      character(len=*), intent(IN)                        :: prefix, name
      type(MEF90HeatXferCellSetOptions_Type), intent(IN)  :: default
      PetscErrorCode, intent(INOUT)                       :: ierr

      type(MEF90HeatXferCellSetOptions_Type), pointer      :: HeatXferCellSetOptions
      PetscCall(PetscBagGetDataMEF90HeatXferCtxCellSetOptions(bag, HeatXferCellSetOptions, ierr))
      PetscCall(PetscBagSetName(bag, trim(name), "HeatXferCellSetOptions MEF90 Heat transfer Cell Set options", ierr))
      PetscCall(PetscBagSetOptionsPrefix(bag, trim(prefix), ierr))

      HeatXferCellSetOptions%advectionVector = default%advectionVector
      PetscCall(PetscBagRegisterReal(bag, HeatXferCellSetOptions%Flux, default%Flux, 'Flux', '[J.s^(-1).m^(-3) / J.s^(-1).m^(-2)] (f): Internal / applied heat flux', ierr))
      PetscCall(PetscBagRegisterBool(bag, HeatXferCellSetOptions%Has_BC, default%Has_BC, 'TemperatureBC', 'Temperature has Dirichlet boundary Condition (Y/N)', ierr))
      PetscCall(PetscBagRegisterReal(bag, HeatXferCellSetOptions%boundaryTemperature, default%boundaryTemperature, 'boundaryTemperature', 'Temperature boundary value', ierr))
      PetscCall(PetscBagRegisterRealArray(bag, HeatXferCellSetOptions%advectionVector, 3_ki, 'advectionVector', '[m.s^(-1)] (V): advection vector', ierr))
   end subroutine PetscBagRegisterMEF90HeatXferCtxCellSetOptions

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90HeatXferCtxFaceSetOptions"
!!!
!!!
!!!  PetscBagRegisterMEF90HeatXferCtxFaceSetOptions:
!!!
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine PetscBagRegisterMEF90HeatXferCtxFaceSetOptions(bag, name, prefix, default, ierr)
      type(tPetscBag)                                    :: bag
      character(len=*), intent(IN)                        :: prefix, name
      type(MEF90HeatXferFaceSetOptions_Type), intent(IN)  :: default
      PetscErrorCode, intent(INOUT)                       :: ierr

      type(MEF90HeatXferFaceSetOptions_Type), pointer      :: HeatXferFaceSetOptions
      PetscCall(PetscBagGetDataMEF90HeatXferCtxFaceSetOptions(bag, HeatXferFaceSetOptions, ierr))
      PetscCall(PetscBagSetName(bag, trim(name), "HeatXferFaceSetOptions MEF90 Heat transfer Face Set options", ierr))
      PetscCall(PetscBagSetOptionsPrefix(bag, trim(prefix), ierr))

      PetscCall(PetscBagRegisterReal(bag, HeatXferFaceSetOptions%boundaryFlux, default%boundaryFlux, 'boundaryFlux', '[J.s^(-1).m^(-2) / J.s^(-1).m^(-1)] (f): Internal / applied heat flux', ierr))
      PetscCall(PetscBagRegisterReal(bag, HeatXferFaceSetOptions%SurfaceThermalConductivity, default%SurfaceThermalConductivity, 'SurfaceThermalConductivity', '[J.s^(-1).m^(-2).K^(-1) / J.s^(-1).m^(-1).K^(-1) ] (H) Surface Thermal Conductivity', ierr))
      PetscCall(PetscBagRegisterReal(bag, HeatXferFaceSetOptions%externalTemperature, default%externalTemperature, 'externalTemperature', 'Reference temperature T [K]', ierr))
      PetscCall(PetscBagRegisterBool(bag, HeatXferFaceSetOptions%Has_BC, default%Has_BC, 'TemperatureBC', 'Temperature has Dirichlet boundary Condition (Y/N)', ierr))
      PetscCall(PetscBagRegisterReal(bag, HeatXferFaceSetOptions%boundaryTemperature, default%boundaryTemperature, 'boundaryTemperature', 'Temperature boundary value', ierr))
   end subroutine PetscBagRegisterMEF90HeatXferCtxFaceSetOptions

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90HeatXferCtxVertexSetOptions"
!!!
!!!
!!!  PetscBagRegisterMEF90HeatXferCtxVertexSetOptions:
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   subroutine PetscBagRegisterMEF90HeatXferCtxVertexSetOptions(bag, name, prefix, default, ierr)
      type(tPetscBag)                                       :: bag
      character(len=*), intent(IN)                           :: prefix, name
      type(MEF90HeatXferVertexSetOptions_Type), intent(IN)   :: default
      PetscErrorCode, intent(INOUT)                          :: ierr

      type(MEF90HeatXferVertexSetOptions_Type), pointer      :: HeatXferVertexSetOptions
      PetscCall(PetscBagGetDataMEF90HeatXferCtxVertexSetOptions(bag, HeatXferVertexSetOptions, ierr))
      PetscCall(PetscBagSetName(bag, trim(name), "HeatXferVertexSetOptions MEF90 Heat transfer Vertex Set options", ierr))
      PetscCall(PetscBagSetOptionsPrefix(bag, trim(prefix), ierr))

      PetscCall(PetscBagRegisterBool(bag, HeatXferVertexSetOptions%Has_BC, default%Has_BC, 'TemperatureBC', 'Temperature has Dirichlet boundary Condition (Y/N)', ierr))
      PetscCall(PetscBagRegisterReal(bag, HeatXferVertexSetOptions%boundaryTemperature, default%boundaryTemperature, 'boundaryTemperature', 'Temperature boundary value', ierr))
   end subroutine PetscBagRegisterMEF90HeatXferCtxVertexSetOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferCtxSetFromOptions"
!!!
!!!
!!!  MEF90HeatXferCtxSetFromOptions:
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   subroutine MEF90HeatXferCtxSetFromOptions(heatXferCtx, prefix, defaultGlobalOptions, &
                                             defaultCellSetOptions, defaultFaceSetOptions, defaultVertexSetOptions, ierr)
      type(MEF90HeatXferCtx_Type), intent(INOUT)             :: heatXferCtx
      character(len=*), intent(IN)                           :: prefix
      type(MEF90HeatXferGlobalOptions_Type), intent(IN)      :: defaultGlobalOptions
      type(MEF90HeatXferCellSetOptions_Type), intent(IN)     :: defaultCellSetOptions
      type(MEF90HeatXferFaceSetOptions_Type), intent(IN)     :: defaultFaceSetOptions
      type(MEF90HeatXferVertexSetOptions_Type), intent(IN)   :: defaultVertexSetOptions
      PetscErrorCode, intent(INOUT)                          :: ierr

      type(MEF90CtxGlobalOptions_Type), pointer              :: MEF90CtxGlobalOptions
      type(tIS)                                             :: setIS
      PetscInt, dimension(:), pointer                         :: setID
      PetscInt                                              :: set
      character(len=MEF90MXSTRLEN)                          :: IOBuffer, setName, setprefix

      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(heatXferCtx%MEF90Ctx%GlobalOptionsBag, MEF90CtxGlobalOptions, ierr))
      !!!
      !!! Registering Global Context
      !!!
      PetscCall(PetscBagRegisterMEF90HeatXferCtxGlobalOptions(heatXferCtx%GlobalOptionsBag, "MEF90HeatXfer Global Ctx", prefix, defaultGlobalOptions, ierr))

      if (MEF90CtxGlobalOptions%verbose > 0) then
         PetscCall(PetscBagView(heatXferCtx%GlobalOptionsBag, PETSC_VIEWER_STDOUT_WORLD, ierr))
      end if

      !!!
      !!! Registering Cell Set Context
      !!!
      PetscCall(DMGetLabelIdIS(heatXferCtx%megaDM, MEF90CellSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(heatXferCtx%MEF90Ctx%comm, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         write (setName, "('Cell set ',I4)") setID(set)
         write (setprefix, "('cs',I4.4,'_')") setID(set)
         PetscCall(PetscBagRegisterMEF90HeatXferCtxCellSetOptions(heatXferCtx%CellSetOptionsBag(set), setName, setPrefix, defaultCellSetOptions, ierr))
         if (MEF90CtxGlobalOptions%verbose > 0) then
            write (IOBuffer, "('\nRegistering cell set ',I4,' prefix: ',A,'\n')") setID(set), trim(setprefix)
            PetscCall(PetscPrintf(heatXferCtx%MEF90Ctx%comm, IOBuffer, ierr))
            PetscCall(PetscBagView(heatXferCtx%CellSetOptionsBag(set), PETSC_VIEWER_STDOUT_WORLD, ierr))
            PetscCall(PetscPrintf(heatXferCtx%MEF90Ctx%comm, "\n", ierr))
         end if
      end do
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))

      !!!
      !!! Registering Face Set Context
      !!!
      PetscCall(DMGetLabelIdIS(heatXferCtx%megaDM, MEF90FaceSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(heatXferCtx%MEF90Ctx%comm, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         write (setName, "('Face set ',I4)") setID(set)
         write (setprefix, "('fs',I4.4,'_')") setID(set)
         PetscCall(PetscBagRegisterMEF90HeatXferCtxFaceSetOptions(heatXferCtx%FaceSetOptionsBag(set), setName, setPrefix, defaultFaceSetOptions, ierr))
         if (MEF90CtxGlobalOptions%verbose > 0) then
            write (IOBuffer, "('\nRegistering face set ',I4,' prefix: ',A,'\n')") setID(set), trim(setprefix)
            PetscCall(PetscPrintf(heatXferCtx%MEF90Ctx%comm, IOBuffer, ierr))
            PetscCall(PetscBagView(heatXferCtx%FaceSetOptionsBag(set), PETSC_VIEWER_STDOUT_WORLD, ierr))
            PetscCall(PetscPrintf(heatXferCtx%MEF90Ctx%comm, "\n", ierr))
         end if
      end do
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))

      !!!
      !!! Registering Vertex Set Context
      !!!
      PetscCall(DMGetLabelIdIS(heatXferCtx%megaDM, MEF90VertexSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(PETSC_COMM_WORLD, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         write (setName, "('Vertex set ',I4)") setID(set)
         write (setprefix, "('vs',I4.4,'_')") setID(set)
         PetscCall(PetscBagRegisterMEF90HeatXferCtxVertexSetOptions(heatXferCtx%VertexSetOptionsBag(set), setName, setPrefix, defaultVertexSetOptions, ierr))
         if (MEF90CtxGlobalOptions%verbose > 0) then
            write (IOBuffer, "('\nRegistering vertex set ',I4,' prefix: ',A,'\n')") setID(set), trim(setprefix)
            PetscCall(PetscPrintf(heatXferCtx%MEF90Ctx%comm, IOBuffer, ierr))
            PetscCall(PetscBagView(heatXferCtx%VertexSetOptionsBag(set), PETSC_VIEWER_STDOUT_WORLD, ierr))
            PetscCall(PetscPrintf(heatXferCtx%MEF90Ctx%comm, "\n", ierr))
         end if
      end do
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end subroutine MEF90HeatXferCtxSetFromOptions
end module m_MEF90_HeatXferCtx
