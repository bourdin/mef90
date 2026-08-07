#include "../MEF90/mef90.inc"
module m_MEF90_HeatXfer_class
#include "petsc/finclude/petsc.h"
   use petscsys
   use m_MEF90_BaseClass
   use m_MEF90_Ctx
   use m_MEF90_DMPlex
   use m_MEF90_LinAlg
   use, intrinsic :: iso_c_binding
   implicit none(type)

!!!
!!!
!!!  MEF90HeatXfer_Type: The class holding the state and the options of a heat transfer problem
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2026    Blaise Bourdin bourdin@mcmaster.ca
!!!

   enum, bind(c)
      enumerator  :: MEF90HeatXfer_timeSteppingTypeNULL = 0, &
         MEF90HeatXFer_timeSteppingTypeSteadyState, &
         MEF90HeatXFer_timeSteppingTypeTransient
   end enum
   character(len=MEF90MXSTRLEN), dimension(6), protected   :: MEF90HeatXFer_timeSteppingTypeList = [character(len=MEF90MXSTRLEN) :: &
      'null', &
      'SteadyState', &
      'Transient', &
      'MEF90_HeatXFer_timeSteppingType', &
      '_MEF90_HeatXFer_timeSteppingType', &
      '']

!!!
!!!  MEF90HeatXferGlobalOptions_Type: the problem-wide options of a MEF90HeatXfer_Type.
!!!                                   The values given here are the defaults used by setFromOptions
!!!
   type MEF90HeatXferGlobalOptions_Type
      PetscEnum                        :: timeSteppingType = MEF90HeatXFer_timeSteppingTypeSteadyState
      PetscBool                        :: addNullSpace = PETSC_FALSE
      PetscReal                        :: initialTemperature = 0.0_kr
      !!! scaling = time (step) scaling law currently CST, Linear, Null (not present), File
      PetscEnum                        :: boundaryTemperatureScaling = MEF90Scaling_Linear
      PetscEnum                        :: externalTemperatureScaling = MEF90Scaling_Linear
      PetscEnum                        :: fluxScaling = MEF90Scaling_Linear
      PetscEnum                        :: boundaryFluxScaling = MEF90Scaling_Linear
      PetscBool                        :: temperatureExport = PETSC_TRUE
   end type MEF90HeatXferGlobalOptions_Type

   type MEF90HeatXferCellSetOptions_Type
      PetscReal                        :: flux = 0.0_kr
      PetscBool                        :: Has_BC = PETSC_FALSE
      PetscReal                        :: boundaryTemperature = 0.0_kr
      PetscReal, dimension(3)          :: advectionVector = [0.0_kr, 0.0_kr, 0.0_kr]
      !!! Material properties. thermalConductivity is a symmetric matrix, so it is dimension
      !!! dependent and allocated by MEF90HeatXferCellSetOptionsSetFromOptions
      PetscReal                        :: density = 1.0_kr
      PetscReal                        :: specificHeat = 1.0_kr
      class(mef90Mat), allocatable     :: thermalConductivity
   end type MEF90HeatXferCellSetOptions_Type

   type MEF90HeatXferFaceSetOptions_Type
      PetscReal                        :: boundaryFlux = 0.0_kr
      PetscReal                        :: surfaceThermalConductivity = 0.0_kr
      PetscReal                        :: externalTemperature = 0.0_kr
      PetscBool                        :: Has_BC = PETSC_FALSE
      PetscReal                        :: boundaryTemperature = 0.0_kr
   end type MEF90HeatXferFaceSetOptions_Type

   type MEF90HeatXferVertexSetOptions_Type
      PetscBool                        :: Has_BC = PETSC_FALSE
      PetscReal                        :: boundaryTemperature = 0.0_kr
   end type MEF90HeatXferVertexSetOptions_Type

   type, extends(MEF90Object) :: MEF90HeatXfer_Type
      type(MEF90Ctx_Type), pointer            :: MEF90Ctx => null()
      type(tDM)                               :: megaDM
      PetscInt                                :: dim

      type(tVec), pointer                     :: temperatureLocal => null()
      type(tVec), pointer                     :: externalTemperatureLocal => null()
      type(tVec), pointer                     :: fluxLocal => null()
      type(tVec), pointer                     :: boundaryFluxLocal => null()

      type(tPetscViewer)                      :: viewer
      type(tPetscSF)                          :: temperatureToIOSF, IOToTemperatureSF
      type(tPetscSF)                          :: boundaryToTemperatureSF
      type(tPetscSF)                          :: externalTemperatureToIOSF, IOToexternalTemperatureSF
      type(tPetscSF)                          :: fluxToIOSF, IOTofluxSF
      type(tPetscSF)                          :: boundaryFluxToIOSF, IOToboundaryFluxSF

      type(MEF90HeatXferGlobalOptions_Type)   :: globalOptions
      !! Per-set options are not stored: they are read from the options database where they are needed,
      !! with MEF90HeatXfer[Cell,Face,Vertex]SetOptionsSetFromOptions. The number of sets is obtained
      !! from the megaDM with MEF90DMGetNumSets.

      !!! Handle on self, set once in MEF90HeatXferCreate with PETScCtx = c_loc(HeatXfer), and handed to
      !!! PETSc wherever an application context is expected: SNESSetFunction, SNESSetJacobian,
      !!! TSSetIFunction, TSSetIJacobian, ... The callbacks in m_MEF90_HeatXfer recover the context with
      !!!    call c_f_pointer(PETScCtx, MEF90HeatXferCtx)
      !!!
      !!! MEF90HeatXfer_Type cannot be passed to PETSc directly, the way MEF90HeatXferCtx_Type was up to
      !!! mef90 0.5.2, because extending MEF90Object gives it type-bound procedures and PETSc declares
      !!! its context arguments as assumed-type. This has to be a component rather than a local or an
      !!! inline c_loc(), for lifetime reasons that are easy to get wrong: see the note at the top of
      !!! MEF90/m_MEF90_BaseClass.F90 before changing any of this.
      type(c_ptr)                             :: PETScCtx = C_NULL_PTR
   contains
      procedure, pass(self) :: setFromOptions => MEF90HeatXferSetFromOptions
      procedure, pass(self) :: view_internal => MEF90HeatXferView
   end type MEF90HeatXfer_Type

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferCreate"
!!!
!!!
!!!  MEF90HeatXferCreate: allocates all the fields of a MEF90HeatXfer_Type
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2026    Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90HeatXferCreate(HeatXfer, dm, MEF90Ctx, prefix, ierr)
      !!! HeatXfer has the target attribute so that HeatXfer%PETScCtx can be handed over to PETSc as an
      !!! application context. The actual argument MUST therefore also have the target attribute.
      type(MEF90HeatXfer_Type), target, intent(OUT) :: HeatXfer
      type(tDM), target, intent(IN)                 :: dm
      type(MEF90Ctx_Type), target, intent(IN)       :: MEF90Ctx
      character(len=*), intent(IN)                  :: prefix
      PetscErrorCode, intent(INOUT)                 :: ierr

      type(MEF90CtxGlobalOptions_Type)              :: MEF90GlobalOptions
      character(len=MEF90MXSTRLEN)                  :: vecName
      type(tDM), dimension(:), pointer              :: dmList
      type(tPetscSF)                                :: dummySF

      PetscCallA(MEF90CtxGlobalOptionsSetFromOptions(MEF90Ctx%comm, trim(MEF90Ctx%prefix), MEF90GlobalOptions, ierr))
      HeatXfer%MEF90Ctx => MEF90Ctx
      HeatXfer%comm = MEF90Ctx%comm
      HeatXfer%prefix = prefix
      HeatXfer%name = trim(prefix)//"HeatXfer"
      HeatXfer%PETScCtx = c_loc(HeatXfer)

      PetscCall(DMGetDimension(dm, HeatXfer%dim, ierr))

      vecName = "Temperature"
      allocate (HeatXfer%temperatureLocal)
      PetscCall(MEF90CreateLocalVector(dm, MEF90GlobalOptions%elementFamily, MEF90GlobalOptions%elementOrder, 1_ki, vecName, HeatXfer%temperatureLocal, ierr))
      vecName = "ExternalTemperature"
      allocate (HeatXfer%externalTemperatureLocal)
      PetscCall(MEF90CreateBoundaryCellVector(dm, 1_ki, vecName, HeatXfer%externalTemperatureLocal, ierr))
      vecName = "Flux"
      allocate (HeatXfer%fluxLocal)
      PetscCall(MEF90CreateCellVector(dm, 1_ki, vecName, HeatXfer%fluxLocal, ierr))
      vecName = "BoundaryFlux"
      allocate (HeatXfer%boundaryFluxLocal)
      PetscCall(MEF90CreateBoundaryCellVector(dm, 1_ki, vecName, HeatXfer%boundaryFluxLocal, ierr))

      !! Create the  unknowns and parameters superDM
      allocate (dmList(4))
      PetscCall(VecGetDM(HeatXfer%temperatureLocal, dmList(1), ierr))
      PetscCall(VecGetDM(HeatXfer%externalTemperatureLocal, dmList(2), ierr))
      PetscCall(VecGetDM(HeatXfer%fluxLocal, dmList(3), ierr))
      PetscCall(VecGetDM(HeatXfer%boundaryFluxLocal, dmList(4), ierr))
      PetscCall(DMCreateSuperDM(dmList, 4_ki, PETSC_NULL_IS_POINTER, HeatXfer%megaDM, ierr))
      deallocate (dmList)

      ! !!! Create the IO SF for all fields
      PetscCall(MEF90IOSFCreate(MEF90Ctx, HeatXfer%TemperatureLocal, HeatXfer%temperatureToIOSF, HeatXfer%IOToTemperatureSF, ierr))
      ! PetscCall(MEF90IOSFCreate(MEF90Ctx,HeatXfer%externalTemperatureLocal,HeatXfer%externalTemperatureToIOSF,HeatXfer%IOToExternalTemperatureSF,ierr))
      PetscCall(MEF90IOSFCreate(MEF90Ctx, HeatXfer%fluxLocal, HeatXfer%fluxToIOSF, HeatXfer%IOToFluxSF, ierr))
      ! PetscCall(MEF90IOSFCreate(MEF90Ctx,HeatXfer%boundaryFluxLocal,HeatXfer%boundaryFluxToIOSF,HeatXfer%IOToBoundaryFluxSF,ierr))

      !! Create the SF to exchange boundary values of the temperature.
      PetscCall(MEF90ConstraintSFCreate(HeatXfer%MEF90Ctx, HeatXfer%TemperatureLocal, HeatXfer%temperatureLocal, HeatXfer%boundaryToTemperatureSF, dummySF, ierr))
      PetscCall(PetscSFDestroy(dummySF, ierr))
   end subroutine MEF90HeatXferCreate

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferDestroy"
!!!
!!!
!!!  MEF90HeatXferDestroy: destroys a MEF90HeatXfer_Type
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2026    Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90HeatXferDestroy(HeatXfer, ierr)
      type(MEF90HeatXfer_Type), intent(INOUT)          :: HeatXfer
      PetscErrorCode, intent(INOUT)                    :: ierr

      HeatXfer%PETScCtx = C_NULL_PTR
      if (associated(HeatXfer%temperatureLocal)) then
         PetscCall(VecDestroy(HeatXfer%temperatureLocal, ierr))
         deallocate (HeatXfer%temperatureLocal)
         nullify (HeatXfer%temperatureLocal)
      end if
      if (associated(HeatXfer%ExternalTemperatureLocal)) then
         PetscCall(VecDestroy(HeatXfer%ExternalTemperatureLocal, ierr))
         deallocate (HeatXfer%ExternalTemperatureLocal)
         nullify (HeatXfer%ExternalTemperatureLocal)
      end if
      if (associated(HeatXfer%fluxLocal)) then
         PetscCall(VecDestroy(HeatXfer%fluxLocal, ierr))
         deallocate (HeatXfer%fluxLocal)
         nullify (HeatXfer%fluxLocal)
      end if
      if (associated(HeatXfer%boundaryFluxLocal)) then
         PetscCall(VecDestroy(HeatXfer%boundaryFluxLocal, ierr))
         deallocate (HeatXfer%boundaryFluxLocal)
         nullify (HeatXfer%boundaryFluxLocal)
      end if

      !! Destroy SFs
      ! PetscCall(PetscSFDestroy(HeatXfer%temperatureToIOSF, ierr))
      ! PetscCall(PetscSFDestroy(HeatXfer%IOToTemperatureSF, ierr))
      ! PetscCall(PetscSFDestroy(HeatXfer%externalTemperatureToIOSF, ierr))
      ! PetscCall(PetscSFDestroy(HeatXfer%IOToExternalTemperatureSF, ierr))
      ! PetscCall(PetscSFDestroy(HeatXfer%fluxToIOSF, ierr))
      ! PetscCall(PetscSFDestroy(HeatXfer%IOToFluxSF, ierr))
      ! PetscCall(PetscSFDestroy(HeatXfer%boundaryFluxToIOSF, ierr))
      ! PetscCall(PetscSFDestroy(HeatXfer%IOToBoundaryFluxSF, ierr))

      PetscCall(PetscSFDestroy(HeatXfer%boundaryToTemperatureSF, ierr))

      nullify (HeatXfer%MEF90Ctx)

      PetscCall(DMDestroy(HeatXfer%megaDM, ierr))
   end subroutine MEF90HeatXferDestroy

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferGlobalOptionsSetFromOptions_Private"
!!!
!!!
!!!  MEF90HeatXferGlobalOptionsSetFromOptions_Private: reads the problem-wide options of a MEF90HeatXfer_Type
!!!
!!!  (c) 2026 Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90HeatXferGlobalOptionsSetFromOptions_Private(comm, prefix, options, ierr)
      MPIU_Comm, intent(IN)                                   :: comm
      character(len=*), intent(IN)                            :: prefix
      type(MEF90HeatXferGlobalOptions_Type), intent(INOUT)    :: options
      PetscErrorCode, intent(INOUT)                           :: ierr

      PetscCall(PetscOptionsBegin(comm, prefix, "Options for MEF90HeatXfer_Type", "mef90HeatXfer", ierr))
         PetscCall(PetscOptionsBool('-temperature_export', 'Export temperature', 'mef90HeatXfer', options%temperatureExport, options%temperatureExport, PETSC_NULL_BOOL, ierr))

         PetscCall(PetscOptionsEnum('-heatxfer_timeStepping_type', 'Type of heat transfer computation', 'mef90HeatXfer', MEF90HeatXFer_timeSteppingTypeList, options%timeSteppingType, options%timeSteppingType, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsBool('-heatxfer_addNullSpace', 'Add null space to SNES', 'mef90HeatXfer', options%addNullSpace, options%addNullSpace, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-heatxfer_initialTemperature', '[K] (T): Initial Temperature', 'mef90HeatXfer', options%initialTemperature, options%initialTemperature, PETSC_NULL_BOOL, ierr))

         PetscCall(PetscOptionsEnum('-boundaryTemperature_scaling', 'Boundary temperature scaling', 'mef90HeatXfer', MEF90ScalingList, options%boundaryTemperatureScaling, options%boundaryTemperatureScaling, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsEnum('-externalTemperature_scaling', 'External Temperature scaling', 'mef90HeatXfer', MEF90ScalingList, options%externalTemperatureScaling, options%externalTemperatureScaling, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsEnum('-flux_scaling', 'Heat flux scaling', 'mef90HeatXfer', MEF90ScalingList, options%fluxScaling, options%fluxScaling, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsEnum('-boundaryFlux_scaling', 'Boundary heat flux scaling', 'mef90HeatXfer', MEF90ScalingList, options%boundaryFluxScaling, options%boundaryFluxScaling, PETSC_NULL_BOOL, ierr))
      PetscCall(PetscOptionsEnd(ierr))
   end subroutine MEF90HeatXferGlobalOptionsSetFromOptions_Private

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferCellSetOptionsSetFromOptions"
!!!
!!!
!!!  MEF90HeatXferCellSetOptionsSetFromOptions: reads the options of a single cell set
!!!
!!!  (c) 2026 Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90HeatXferCellSetOptionsSetFromOptions(comm, prefix, dim, options, ierr)
      MPIU_Comm, intent(IN)                                   :: comm
      character(len=*), intent(IN)                            :: prefix
      PetscInt, intent(IN)                                    :: dim
      type(MEF90HeatXferCellSetOptions_Type), intent(OUT)     :: options
      PetscErrorCode, intent(INOUT)                           :: ierr

      PetscInt                                                :: nOpt
      PetscReal, dimension(:), allocatable                    :: tmpArray

      select case (dim)
      case (2)
         options%thermalConductivity = MEF90MatS2DIdentity
      case (3)
         options%thermalConductivity = MEF90MatS3DIdentity
      case default
         SETERRQ(comm, PETSC_ERR_ARG_OUTOFRANGE, "dim must be 2 or 3 in " // __FUNCT__)
      end select

      PetscCall(PetscOptionsBegin(comm, prefix, "Options for a MEF90HeatXfer_Type cell set", "mef90HeatXfer", ierr))
         PetscCall(PetscOptionsReal('-Flux', '[J.s^(-1).m^(-3) / J.s^(-1).m^(-2)] (f): Internal / applied heat flux', 'mef90HeatXfer', options%Flux, options%Flux, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsBool('-TemperatureBC', 'Temperature has Dirichlet boundary Condition (Y/N)', 'mef90HeatXfer', options%Has_BC, options%Has_BC, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-boundaryTemperature', 'Temperature boundary value', 'mef90HeatXfer', options%boundaryTemperature, options%boundaryTemperature, PETSC_NULL_BOOL, ierr))
         nOpt = 3
         PetscCall(PetscOptionsRealArray('-advectionVector', '[m.s^(-1)] (V): advection vector', 'mef90HeatXfer', options%advectionVector, nOpt, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-Density', '[kg.m^(-dim)] (rho) Density', 'mef90HeatXfer', options%density, options%density, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-SpecificHeat', '[J.kg^(-1).K^(-1)] (Cp) Specific heat', 'mef90HeatXfer', options%specificHeat, options%specificHeat, PETSC_NULL_BOOL, ierr))
         select type (K => options%thermalConductivity)
         type is (MatS2D)
            nOpt = 3
            allocate (tmpArray(nOpt))
            tmpArray = K
            PetscCall(PetscOptionsRealArray('-ThermalConductivity', '[J.m^(-1).s^(-1).K^(-1)] (K) Thermal conductivity', 'mef90HeatXfer', tmpArray, nOpt, PETSC_NULL_BOOL, ierr))
            K = tmpArray
            deallocate (tmpArray)
         type is (MatS3D)
            nOpt = 6
            allocate (tmpArray(nOpt))
            tmpArray = K
            PetscCall(PetscOptionsRealArray('-ThermalConductivity', '[J.m^(-1).s^(-1).K^(-1)] (K) Thermal conductivity', 'mef90HeatXfer', tmpArray, nOpt, PETSC_NULL_BOOL, ierr))
            K = tmpArray
            deallocate (tmpArray)
         end select
      PetscCall(PetscOptionsEnd(ierr))
   end subroutine MEF90HeatXferCellSetOptionsSetFromOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferFaceSetOptionsSetFromOptions"
!!!
!!!
!!!  MEF90HeatXferFaceSetOptionsSetFromOptions: reads the options of a single face set
!!!
!!!  (c) 2026 Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90HeatXferFaceSetOptionsSetFromOptions(comm, prefix, options, ierr)
      MPIU_Comm, intent(IN)                                   :: comm
      character(len=*), intent(IN)                            :: prefix
      type(MEF90HeatXferFaceSetOptions_Type), intent(OUT)     :: options
      PetscErrorCode, intent(INOUT)                           :: ierr

      PetscCall(PetscOptionsBegin(comm, prefix, "Options for a MEF90HeatXfer_Type face set", "mef90HeatXfer", ierr))
         PetscCall(PetscOptionsReal('-boundaryFlux', '[J.s^(-1).m^(-2) / J.s^(-1).m^(-1)] (f): Internal / applied heat flux', 'mef90HeatXfer', options%boundaryFlux, options%boundaryFlux, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-SurfaceThermalConductivity', '[J.s^(-1).m^(-2).K^(-1) / J.s^(-1).m^(-1).K^(-1) ] (H) Surface Thermal Conductivity', 'mef90HeatXfer', options%SurfaceThermalConductivity, options%SurfaceThermalConductivity, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-externalTemperature', 'Reference temperature T [K]', 'mef90HeatXfer', options%externalTemperature, options%externalTemperature, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsBool('-TemperatureBC', 'Temperature has Dirichlet boundary Condition (Y/N)', 'mef90HeatXfer', options%Has_BC, options%Has_BC, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-boundaryTemperature', 'Temperature boundary value', 'mef90HeatXfer', options%boundaryTemperature, options%boundaryTemperature, PETSC_NULL_BOOL, ierr))
      PetscCall(PetscOptionsEnd(ierr))
   end subroutine MEF90HeatXferFaceSetOptionsSetFromOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferVertexSetOptionsSetFromOptions"
!!!
!!!
!!!  MEF90HeatXferVertexSetOptionsSetFromOptions: reads the options of a single vertex set
!!!
!!!  (c) 2026 Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90HeatXferVertexSetOptionsSetFromOptions(comm, prefix, options, ierr)
      MPIU_Comm, intent(IN)                                   :: comm
      character(len=*), intent(IN)                            :: prefix
      type(MEF90HeatXferVertexSetOptions_Type), intent(OUT)   :: options
      PetscErrorCode, intent(INOUT)                           :: ierr

      PetscCall(PetscOptionsBegin(comm, prefix, "Options for a MEF90HeatXfer_Type vertex set", "mef90HeatXfer", ierr))
         PetscCall(PetscOptionsBool('-TemperatureBC', 'Temperature has Dirichlet boundary Condition (Y/N)', 'mef90HeatXfer', options%Has_BC, options%Has_BC, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-boundaryTemperature', 'Temperature boundary value', 'mef90HeatXfer', options%boundaryTemperature, options%boundaryTemperature, PETSC_NULL_BOOL, ierr))
      PetscCall(PetscOptionsEnd(ierr))
   end subroutine MEF90HeatXferVertexSetOptionsSetFromOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferSetFromOptions"
!!!
!!!
!!!  MEF90HeatXferSetFromOptions: initializes a MEF90HeatXfer_Type from options
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2026    Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90HeatXferSetFromOptions(self, ierr)
      class(MEF90HeatXfer_Type), intent(INOUT)                :: self
      PetscErrorCode, intent(INOUT)                           :: ierr

      type(tIS)                                               :: setIS
      PetscInt, dimension(:), pointer                         :: setID
      PetscInt                                                :: set
      character(len=MEF90MXSTRLEN)                            :: setPrefix
      type(MEF90HeatXferCellSetOptions_Type)                  :: cellSetOptions
      type(MEF90HeatXferFaceSetOptions_Type)                  :: faceSetOptions
      type(MEF90HeatXferVertexSetOptions_Type)                :: vertexSetOptions
      PetscBool                                               :: printHelp

      !!
      !! Problem-wide options
      !!
      PetscCall(MEF90HeatXferGlobalOptionsSetFromOptions_Private(self%comm, trim(self%prefix), self%globalOptions, ierr))

      !!
      !! The per-set options are read where they are needed, but they have to be visited once here so
      !! that they are registered with the options database
      !!
      PetscCall(DMGetLabelIdIS(self%megaDM, MEF90CellSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(self%comm, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         write (setPrefix, "(A,'cs',I4.4,'_')") trim(self%prefix), setID(set)
         PetscCall(MEF90HeatXferCellSetOptionsSetFromOptions(self%comm, trim(setPrefix), self%dim, cellSetOptions, ierr))
      end do
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))

      PetscCall(DMGetLabelIdIS(self%megaDM, MEF90FaceSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(self%comm, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         write (setPrefix, "(A,'fs',I4.4,'_')") trim(self%prefix), setID(set)
         PetscCall(MEF90HeatXferFaceSetOptionsSetFromOptions(self%comm, trim(setPrefix), faceSetOptions, ierr))
      end do
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))

      PetscCall(DMGetLabelIdIS(self%megaDM, MEF90VertexSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(self%comm, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         write (setPrefix, "(A,'vs',I4.4,'_')") trim(self%prefix), setID(set)
         PetscCall(MEF90HeatXferVertexSetOptionsSetFromOptions(self%comm, trim(setPrefix), vertexSetOptions, ierr))
      end do
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))

      PetscCall(PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-verbose", printHelp, PETSC_NULL_BOOL, ierr))
      if (printHelp) then
         call self%view(PETSC_VIEWER_STDOUT_WORLD, ierr)
      end if
   end subroutine MEF90HeatXferSetFromOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferView"
!!!
!!!
!!!  MEF90HeatXferView: the default viewer for a MEF90HeatXfer_Type
!!!
!!!  (c) 2026 Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90HeatXferView(self, viewer, ierr)
      class(MEF90HeatXfer_Type), intent(IN)                   :: self
      type(tPetscViewer), intent(IN)                          :: viewer
      PetscErrorCode, intent(INOUT)                           :: ierr

      character(len=MEF90MXSTRLEN, kind=c_char)               :: IOBuffer
      character(len=MEF90MXSTRLEN, kind=c_char)               :: viewerType
      type(tIS)                                               :: setIS
      PetscInt, dimension(:), pointer                         :: setID
      PetscInt                                                :: set
      character(len=MEF90MXSTRLEN)                            :: setPrefix
      type(MEF90HeatXferCellSetOptions_Type)                  :: cellSetOptions
      type(MEF90HeatXferFaceSetOptions_Type)                  :: faceSetOptions
      type(MEF90HeatXferVertexSetOptions_Type)                :: vertexSetOptions

      PetscCall(PetscViewerGetType(viewer, viewerType, ierr))
      if (viewerType /= 'ascii') return

      write (IOBuffer, "(A,': Options for MEF90HeatXfer_Type\n')") trim(self%prefix)//"heatxfer"
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         time stepping type: ',A,'\n')") trim(MEF90HeatXFer_timeSteppingTypeList(self%globalOptions%timeSteppingType + 1))
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         initial temperature: ',ES12.5,' [K]\n')") self%globalOptions%initialTemperature
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         add null space: ',L1,' export temperature: ',L1,'\n')") self%globalOptions%addNullSpace, self%globalOptions%temperatureExport
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         scaling boundary temperature / external temperature: ',A,' / ',A,'\n')") &
         trim(MEF90ScalingList(self%globalOptions%boundaryTemperatureScaling + 1)), trim(MEF90ScalingList(self%globalOptions%externalTemperatureScaling + 1))
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         scaling flux / boundary flux: ',A,' / ',A,'\n')") &
         trim(MEF90ScalingList(self%globalOptions%fluxScaling + 1)), trim(MEF90ScalingList(self%globalOptions%boundaryFluxScaling + 1))
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))

      PetscCall(DMGetLabelIdIS(self%megaDM, MEF90CellSetLabelName, setIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(self%comm, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         write (setPrefix, "(A,'cs',I4.4,'_')") trim(self%prefix), setID(set)
         PetscCall(MEF90HeatXferCellSetOptionsSetFromOptions(self%comm, trim(setPrefix), self%dim, cellSetOptions, ierr))
         write (IOBuffer, "(A,'cs',I4.4,': flux: ',ES12.5,' advection vector: ',3(ES12.5,' '),'\n')") &
            trim(self%prefix), setID(set), cellSetOptions%flux, cellSetOptions%advectionVector
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
         write (IOBuffer, "('         temperature BC: ',L1,' value: ',ES12.5,'\n')") &
            cellSetOptions%Has_BC, cellSetOptions%boundaryTemperature
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      end do
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))

      PetscCall(DMGetLabelIdIS(self%megaDM, MEF90FaceSetLabelName, setIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(self%comm, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         write (setPrefix, "(A,'fs',I4.4,'_')") trim(self%prefix), setID(set)
         PetscCall(MEF90HeatXferFaceSetOptionsSetFromOptions(self%comm, trim(setPrefix), faceSetOptions, ierr))
         write (IOBuffer, "(A,'fs',I4.4,': boundary flux: ',ES12.5,' surface thermal conductivity: ',ES12.5,'\n')") &
            trim(self%prefix), setID(set), faceSetOptions%boundaryFlux, faceSetOptions%surfaceThermalConductivity
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
         write (IOBuffer, "('         external temperature: ',ES12.5,' temperature BC: ',L1,' value: ',ES12.5,'\n')") &
            faceSetOptions%externalTemperature, faceSetOptions%Has_BC, faceSetOptions%boundaryTemperature
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      end do
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))

      PetscCall(DMGetLabelIdIS(self%megaDM, MEF90VertexSetLabelName, setIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(self%comm, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         write (setPrefix, "(A,'vs',I4.4,'_')") trim(self%prefix), setID(set)
         PetscCall(MEF90HeatXferVertexSetOptionsSetFromOptions(self%comm, trim(setPrefix), vertexSetOptions, ierr))
         write (IOBuffer, "(A,'vs',I4.4,': temperature BC: ',L1,' value: ',ES12.5,'\n')") &
            trim(self%prefix), setID(set), vertexSetOptions%Has_BC, vertexSetOptions%boundaryTemperature
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      end do
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end subroutine MEF90HeatXferView
end module m_MEF90_HeatXfer_class
