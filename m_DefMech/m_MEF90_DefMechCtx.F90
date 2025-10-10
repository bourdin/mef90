module m_MEF90_DefMechCtx_Type
#include "petsc/finclude/petsc.h"
   use m_MEF90_Ctx
   implicit none(type, external)

   type MEF90DefMechCtx_Type
      type(MEF90Ctx_Type), pointer             :: MEF90Ctx => null()
      type(tDM)                                :: megaDM
      PetscInt                                 :: dim
      PetscReal                                :: analysisTime, timeStep

      !!!  vertex based vec
      type(tVec), pointer                       :: displacementLocal => null()
      type(tVec), pointer                       :: displacementPreviousStepLocal => null()
      type(tVec), pointer                       :: damageLocal => null()
      type(tVec), pointer                       :: damagePreviousStepLocal => null()
      type(tVec), pointer                       :: displacementLowerBoundLocal => null()
      type(tVec), pointer                       :: displacementUpperBoundLocal => null()
      type(tVec), pointer                       :: temperatureLocal => null()

      !!! cell based vec
      type(tVec), pointer                       :: bodyForce => null()
      type(tVec), pointer                       :: boundaryForce => null()
      type(tVec), pointer                       :: pressureForce => null()
      !Type(tVec),pointer                       :: crackPressure => null()
      type(tVec), pointer                       :: cohesiveDisplacement => null()
      type(tVec), pointer                       :: plasticStrain => null()
      type(tVec), pointer                       :: cumulatedPlasticDissipation => null()
      type(tVec), pointer                       :: stress => null()

      type(tPetscSF)                            :: displacementToIOSF, IOToDisplacementSF
      type(tPetscSF)                            :: cohesiveDisplacementToIOSF, IOToCohesiveDisplacementSF
      type(tPetscSF)                            :: displacementConstraintsSF
      type(tPetscSF)                            :: damageToIOSF, IOToDamageSF
      type(tPetscSF)                            :: damageConstraintsSF
      type(tPetscSF)                            :: temperatureToIOSF, IOToTemperatureSF
      type(tPetscSF)                            :: bodyForceToIOSF, IOToBodyForceSF
      type(tPetscSF)                            :: boundaryForceToIOSF, IOToBoundaryForceSF
      type(tPetscSF)                            :: pressureForceToIOSF, IOToPressureForceSF
      type(tPetscSF)                            :: stressToIOSF, IOToStressSF
      type(tPetscSF)                            :: plasticStrainToIOSF, IOToPlasticStrainSF
      type(tPetscSF)                            :: cumulatedPlasticDissToIOSF, IOToCumulatedPlasticDissSF

      type(tPetscBag)                           :: GlobalOptionsBag
      type(tPetscBag), dimension(:), pointer    :: CellSetOptionsBag => null()
      type(tPetscBag), dimension(:), pointer    :: FaceSetOptionsBag => null()
      type(tPetscBag), dimension(:), pointer    :: VertexSetOptionsBag => null()
      type(tPetscBag), dimension(:), pointer    :: MaterialPropertiesBag => null()

      type(tPetscViewer)                        :: globalEnergyViewer
      type(tPetscViewer), dimension(:), pointer :: setEnergyViewer => null()

      PetscBool                                 :: hasDisplacementBounds
      PetscBool                                 :: hasUnilateralContact
   end type MEF90DefMechCtx_Type

   type MEF90DefMechGlobalOptions_Type
      PetscEnum                              :: timeSteppingType
      PetscEnum                              :: solverType
      PetscEnum                              :: damageSolverType

      !!! scaling = time (step) scaling law currently CST, Linear, or File
      PetscEnum                              :: boundaryDisplacementScaling
      PetscEnum                              :: displacementLowerBoundScaling
      PetscEnum                              :: displacementUpperBoundScaling
      PetscEnum                              :: cohesiveDisplacementScaling
      PetscEnum                              :: boundaryDamageScaling
      PetscEnum                              :: bodyForceScaling
      PetscEnum                              :: boundaryForceScaling
      PetscEnum                              :: pressureForceScaling
      PetscEnum                              :: CrackPressureScaling

      PetscReal                              :: damageATol
      PetscInt                               :: damageMaxIt
      PetscInt                               :: PCLag
      PetscReal                              :: SOROmega
      PetscReal                              :: irrevthres
      PetscEnum                              :: BTType
      PetscInt                               :: BTInterval
      PetscInt                               :: BTScope
      PetscReal                              :: BTTol
      PetscReal                              :: plasticStrainATol
      PetscReal                              :: InjectedVolumeATol
      PetscReal                              :: dampingCoefficientDisplacement
      PetscReal                              :: dampingCoefficientDamage
      PetscBool                              :: temperatureExport
      PetscBool                              :: displacementExport
      PetscBool                              :: damageExport
      PetscBool                              :: stressExport
      PetscBool                              :: plasticStrainExport
      PetscBool                              :: cumulatedPlasticDissipationExport
   end type MEF90DefMechGlobalOptions_Type

   type MEF90DefMechCellSetOptions_Type
      PetscReal, dimension(3)                 :: bodyforce
      PetscReal                               :: crackPressure
      PetscEnum                               :: damageType
      PetscEnum                               :: plasticityType
      PetscEnum                               :: unilateralContactType
      PetscReal                               :: unilateralContactHydrostaticDeviatoricGamma
      PetscBool                               :: unilateralContactHybrid
      PetscReal                               :: DamageATLinSoftk
      PetscReal                               :: DamageAT1expb
      PetscEnum                               :: drivingForceType
      PetscReal, dimension(3)                 :: cohesiveDisplacement
      PetscBool, dimension(3)                 :: Has_displacementBC
      PetscReal, dimension(3)                 :: boundaryDisplacement
      PetscReal, dimension(3)                 :: displacementLowerBound
      PetscReal, dimension(3)                 :: displacementUpperBound
      PetscBool                               :: Has_damageBC
      PetscReal                               :: boundaryDamage
      PetscBool                               :: CrackVolumeControlled
      PetscBool                               :: WorkControlled
   end type MEF90DefMechCellSetOptions_Type

   type MEF90DefMechFaceSetOptions_Type
      PetscReal, dimension(3)                 :: boundaryforce
      PetscReal                               :: pressureForce
      PetscBool, dimension(3)                 :: Has_displacementBC
      PetscReal, dimension(3)                 :: boundaryDisplacement
      PetscReal, dimension(3)                 :: displacementLowerBound
      PetscReal, dimension(3)                 :: displacementUpperBound
      PetscBool                               :: Has_damageBC
      PetscReal                               :: boundaryDamage
   end type MEF90DefMechFaceSetOptions_Type

   type MEF90DefMechVertexSetOptions_Type
      PetscBool, dimension(3)                 :: Has_displacementBC
      PetscReal, dimension(3)                 :: boundaryDisplacement
      PetscReal, dimension(3)                 :: displacementLowerBound
      PetscReal, dimension(3)                 :: displacementUpperBound
      PetscBool                               :: Has_damageBC
      PetscReal                               :: boundaryDamage
   end type MEF90DefMechVertexSetOptions_Type
end module m_MEF90_DefMechCtx_Type

module m_MEF90DefMechGlobalOptions_Private
#include "petsc/finclude/petsc.h"
   use petscbag
   use m_MEF90_DefMechCtx_Type
   implicit none(type, external)

   private
   public :: PetscBagGetDataMEF90DefMechCtxGlobalOptions

   interface PetscBagGetData
      subroutine PetscBagGetData(bag, data, ierr)
         use petscbag
         use m_MEF90_DefMechCtx_Type
         implicit none(type)
         type(tPetscBag)                                       :: bag
         type(MEF90DefMechGlobalOptions_Type), pointer         :: data
         PetscErrorCode                                        :: ierr
      end subroutine PetscBagGetData
   end interface
contains

#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90DefMechCtxGlobalOptions"
!!!
!!!  PetscBagGetDataMEF90DefMechCtxGlobalOptions - Custom interface to PetscGetData
!!!

   subroutine PetscBagGetDataMEF90DefMechCtxGlobalOptions(bag, data, ierr)
      type(tPetscBag), intent(IN)                            :: bag
      type(MEF90DefMechGlobalOptions_Type), pointer          :: data
      PetscErrorCode, intent(INOUT)                          :: ierr

      PetscCall(PetscBagGetData(bag, data, ierr))
   end subroutine PetscBagGetDataMEF90DefMechCtxGlobalOptions
end module m_MEF90DefMechGlobalOptions_Private

module m_MEF90DefMechCellSetOptions_Private
#include "petsc/finclude/petsc.h"
   use petscbag
   use m_MEF90_DefMechCtx_Type
   implicit none(type, external)

   private
   public :: PetscBagGetDataMEF90DefMechCtxCellSetOptions

   interface PetscBagGetData
      subroutine PetscBagGetData(bag, data, ierr)
         use petscbag
         use m_MEF90_DefMechCtx_Type
         implicit none(type, external)

         type(tPetscBag)                                       :: bag
         type(MEF90DefMechCellSetOptions_Type), pointer         :: data
         PetscErrorCode                                        :: ierr
      end subroutine PetscBagGetData
   end interface
contains

#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90DefMechCtxCellSetOptions"
!!!
!!!  PetscBagGetDataMEF90DefMechCtxCellSetOptions - Custom interface to PetscGetData
!!!

   subroutine PetscBagGetDataMEF90DefMechCtxCellSetOptions(bag, data, ierr)
      type(tPetscBag), intent(IN)                            :: bag
      type(MEF90DefMechCellSetOptions_Type), pointer         :: data
      PetscErrorCode, intent(INOUT)                          :: ierr

      PetscCall(PetscBagGetData(bag, data, ierr))
   end subroutine PetscBagGetDataMEF90DefMechCtxCellSetOptions
end module m_MEF90DefMechCellSetOptions_Private

module m_MEF90DefMechFaceSetOptions_Private
#include "petsc/finclude/petsc.h"
   use petscbag
   use m_MEF90_DefMechCtx_Type
   implicit none(type, external)

   private
   public :: PetscBagGetDataMEF90DefMechCtxFaceSetOptions

   interface PetscBagGetData
      subroutine PetscBagGetData(bag, data, ierr)
         use petscbag
         use m_MEF90_DefMechCtx_Type
         implicit none(type, external)

         type(tPetscBag)                                       :: bag
         type(MEF90DefMechFaceSetOptions_Type), pointer         :: data
         PetscErrorCode                                        :: ierr
      end subroutine PetscBagGetData
   end interface
contains

#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90DefMechCtxFaceSetOptions"
!!!
!!!  PetscBagGetDataMEF90DefMechCtxFaceSetOptions - Custom interface to PetscGetData
!!!

   subroutine PetscBagGetDataMEF90DefMechCtxFaceSetOptions(bag, data, ierr)
      type(tPetscBag), intent(IN)                            :: bag
      type(MEF90DefMechFaceSetOptions_Type), pointer         :: data
      PetscErrorCode, intent(INOUT)                          :: ierr

      PetscCall(PetscBagGetData(bag, data, ierr))
   end subroutine PetscBagGetDataMEF90DefMechCtxFaceSetOptions
end module m_MEF90DefMechFaceSetOptions_Private

module m_MEF90DefMechVertexSetOptions_Private
#include "petsc/finclude/petsc.h"
   use petscbag
   use m_MEF90_DefMechCtx_Type
   implicit none(type, external)

   private
   public :: PetscBagGetDataMEF90DefMechCtxVertexSetOptions

   interface PetscBagGetData
      subroutine PetscBagGetData(bag, data, ierr)
         use petscbag
         use m_MEF90_DefMechCtx_Type
         implicit none(type, external)

         type(tPetscBag)                                       :: bag
         type(MEF90DefMechVertexSetOptions_Type), pointer       :: data
         PetscErrorCode, intent(INOUT)                          :: ierr
      end subroutine PetscBagGetData
   end interface

contains
#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90DefMechCtxVertexSetOptions"
!!!
!!!  PetscBagGetDataMEF90DefMechCtxVertexSetOptions - Custom interface to PetscGetData
!!!

   subroutine PetscBagGetDataMEF90DefMechCtxVertexSetOptions(bag, data, ierr)
      type(tPetscBag), intent(IN)                            :: bag
      type(MEF90DefMechVertexSetOptions_Type), pointer       :: data
      PetscErrorCode, intent(INOUT)                          :: ierr

      PetscCall(PetscBagGetData(bag, data, ierr))
   end subroutine PetscBagGetDataMEF90DefMechCtxVertexSetOptions
end module m_MEF90DefMechVertexSetOptions_Private

module m_MEF90_DefMechCtx
#include "petsc/finclude/petsc.h"
   use m_MEF90_DMPlex
   use m_MEF90_DefMechCtx_Type
   use m_MEF90DefMechGlobalOptions_Private
   use m_MEF90DefMechCellSetOptions_Private
   use m_MEF90DefMechFaceSetOptions_Private
   use m_MEF90DefMechVertexSetOptions_Private
   implicit none(type)

   PetscSizeT, protected   :: sizeofMEF90DefMechGlobalOptions
   PetscSizeT, protected   :: sizeofMEF90DefMechCellSetOptions
   PetscSizeT, protected   :: sizeofMEF90DefMechFaceSetOptions
   PetscSizeT, protected   :: sizeofMEF90DefMechVertexSetOptions

   enum, bind(c)
      enumerator :: MEF90DefMech_SolverTypeAltMin = 0, &
         MEF90DefMech_SolverTypeQuasiNewton1, &
         MEF90DefMech_SolverTypeQuasiNewton2
   end enum
   character(len=MEF90MXSTRLEN), dimension(6), protected   :: MEF90DefMech_SolverTypeList

   enum, bind(c)
      enumerator :: MEF90DefMech_TimeSteppingTypeNULL = 0, &
         MEF90DefMech_TimeSteppingTypeQuasiStatic
   end enum
   character(len=MEF90MXSTRLEN), dimension(5), protected   :: MEF90DefMech_TimeSteppingTypeList

   enum, bind(c)
      enumerator :: MEF90DefMech_DamageSolverTypeSNES = 0, &
         MEF90DefMech_DamageSolverTypeTAO
   end enum
   character(len=MEF90MXSTRLEN), dimension(5), protected   :: MEF90DefMech_DamageSolverTypeList

   enum, bind(c)
      enumerator :: MEF90DefMech_BTTypeNULL = 0, &
         MEF90DefMech_BTTypeBackward, &
         MEF90DefMech_BTTypeForward
   end enum
   character(len=MEF90MXSTRLEN), dimension(6), protected   :: MEF90DefMech_BTTypeList

   enum, bind(c)
      enumerator :: MEF90DefMech_damageTypeAT1 = 0, &
         MEF90DefMech_damageTypeAT1exp, &
         MEF90DefMech_damageTypeAT2, &
         MEF90DefMech_damageTypeLinSoft, &
         MEF90DefMech_damageTypeKKL, &
         MEf90DefMech_damageTypeAT1Elastic, &
         MEf90DefMech_damageTypeAT1expElastic, &
         MEf90DefMech_damageTypeAT2Elastic, &
         MEF90DefMech_damageTypeLinSoftElastic, &
         MEF90DefMech_damageTypeKKLElastic

   end enum
   character(len=MEF90MXSTRLEN), dimension(13), protected   :: MEF90DefMech_damageTypeList

   enum, bind(c)
      enumerator :: MEF90DefMech_plasticityTypeNone = 0, &
         MEF90DefMech_plasticityTypeTresca, &
         MEF90DefMech_plasticityTypeVonMises, &
         MEF90DefMech_plasticityTypeVonMisesPlaneTheory, &
         MEF90DefMech_plasticityTypeCapModel, &
         MEF90DefMech_plasticityTypeDruckerPragerCapModel, &
         MEF90DefMech_plasticityTypeVonMises1D, &
         MEF90DefMech_plasticityTypeHillPlaneTheory, &
         MEF90DefMech_PlasticityTypeGreen, &
         MEF90DefMech_PlasticityTypeGurson
   end enum
   character(len=MEF90MXSTRLEN), dimension(13), protected   :: MEF90DefMech_plasticityTypeList

   enum, bind(c)
      enumerator :: MEF90DefMech_unilateralContactTypeNone = 0, &
         MEF90DefMech_unilateralContactTypeHydrostaticDeviatoric, &
         MEF90DefMech_unilateralContactTypeHydrostatic, &
         MEF90DefMech_unilateralContactTypeDeviatoric, &
         MEF90DefMech_unilateralContactTypePrincipalStrains, &
         MEF90DefMech_unilateralContactTypeMasonry
   end enum
   character(len=MEF90MXSTRLEN), dimension(9), protected   :: MEF90DefMech_unilateralContactTypeList

   enum, bind(c)
      enumerator :: MEF90DefMech_drivingForceTypeNone = 0, &
         MEF90_DefMechDrivingForceTypeDruckerPrager, &
         MEF90DefMech_drivingForceTypeDruckerPrager2
   end enum
   character(len=MEF90MXSTRLEN), dimension(6), protected   :: MEF90DefMech_drivingForceTypeList
contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCtxInitialize_Private"
!!!
!!!
!!!  MEF90DefMechCtxInitialize_Private:
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90DefMechCtxInitialize_Private(ierr)
      PetscErrorCode, intent(OUT)                         :: ierr

      type(MEF90DefMechGlobalOptions_Type)               :: DefMechGlobalOptions
      type(MEF90DefMechCellSetOptions_Type)              :: DefMechCellSetOptions
      type(MEF90DefMechFaceSetOptions_Type)              :: DefMechFaceSetOptions
      type(MEF90DefMechVertexSetOptions_Type)            :: DefMechVertexSetOptions
      character(len=1), pointer                          :: dummychar(:)
      PetscSizeT                                         :: sizeofchar

      PetscCall(PetscDataTypeGetSize(PETSC_CHAR, sizeofchar, ierr))
      sizeofMEF90DefMechGlobalOptions = size(transfer(DefMechGlobalOptions, dummychar)) * sizeofchar
      sizeofMEF90DefMechCellSetOptions = size(transfer(DefMechCellSetOptions, dummychar)) * sizeofchar
      sizeofMEF90DefMechFaceSetOptions = size(transfer(DefMechFaceSetOptions, dummychar)) * sizeofchar
      sizeofMEF90DefMechVertexSetOptions = size(transfer(DefMechVertexSetOptions, dummychar)) * sizeofchar

      MEF90DefMech_SolverTypeList(1) = 'AltMin'
      MEF90DefMech_SolverTypeList(2) = 'QuasiNewton1'
      MEF90DefMech_SolverTypeList(3) = 'QuasiNewton2'
      MEF90DefMech_SolverTypeList(4) = 'MEF90DefMech_SolverType'
      MEF90DefMech_SolverTypeList(5) = '_MEF90DefMech_SolverType'
      MEF90DefMech_SolverTypeList(6) = ''

      MEF90DefMech_TimeSteppingTypeList(1) = 'Null'
      MEF90DefMech_TimeSteppingTypeList(2) = 'QuasiStatic'
      MEF90DefMech_TimeSteppingTypeList(3) = 'MEF90DefMech_TimeSteppingType'
      MEF90DefMech_TimeSteppingTypeList(4) = '_MEF90DefMech_TimeSteppingType'
      MEF90DefMech_TimeSteppingTypeList(5) = ''

      MEF90DefMech_DamageSolverTypeList(1) = 'SNES'
      MEF90DefMech_DamageSolverTypeList(2) = 'Tao'
      MEF90DefMech_DamageSolverTypeList(3) = 'MEF90DefMech_DamageSolverType'
      MEF90DefMech_DamageSolverTypeList(4) = '_MEF90DefMech_DamageSolverType'
      MEF90DefMech_DamageSolverTypeList(5) = ''

      MEF90DefMech_BTTypeList(1) = 'Null'
      MEF90DefMech_BTTypeList(2) = 'Backward'
      MEF90DefMech_BTTypeList(3) = 'Forward'
      MEF90DefMech_BTTypeList(4) = 'MEF90DefMech_BTType'
      MEF90DefMech_BTTypeList(5) = '_MEF90DefMech_BTType'
      MEF90DefMech_BTTypeList(6) = ''

      MEF90DefMech_damageTypeList(1) = 'AT1'
      MEF90DefMech_damageTypeList(2) = 'AT1exp'
      MEF90DefMech_damageTypeList(3) = 'AT2'
      MEF90DefMech_damageTypeList(4) = 'LinSoft'
      MEF90DefMech_damageTypeList(5) = 'KKL'
      MEF90DefMech_damageTypeList(6) = 'AT1Elastic'
      MEF90DefMech_damageTypeList(7) = 'AT1expElastic'
      MEF90DefMech_damageTypeList(8) = 'AT2Elastic'
      MEF90DefMech_damageTypeList(9) = 'LinSoftElastic'
      MEF90DefMech_damageTypeList(10) = 'KKLElastic'
      MEF90DefMech_damageTypeList(11) = 'MEF90DefMech_damageType'
      MEF90DefMech_damageTypeList(12) = '_MEF90DefMech_damageType'
      MEF90DefMech_damageTypeList(13) = ''

      MEF90DefMech_plasticityTypeList(1) = 'None'
      MEF90DefMech_plasticityTypeList(2) = 'Tresca'
      MEF90DefMech_plasticityTypeList(3) = 'VonMises'
      MEF90DefMech_plasticityTypeList(4) = 'VonMisesPlaneTheory'
      MEF90DefMech_plasticityTypeList(5) = 'CapModel'
      MEF90DefMech_plasticityTypeList(6) = 'DruckerPragerCapModel'
      MEF90DefMech_plasticityTypeList(7) = 'VonMises1D'
      MEF90DefMech_plasticityTypeList(8) = 'HillPlaneTheory'
      MEF90DefMech_plasticityTypeList(9) = 'Green'
      MEF90DefMech_plasticityTypeList(10) = 'Gurson'
      MEF90DefMech_plasticityTypeList(11) = 'MEF90DefMech_plasticityType'
      MEF90DefMech_plasticityTypeList(12) = '_MEF90DefMech_plasticityType'
      MEF90DefMech_plasticityTypeList(13) = ''

      MEF90DefMech_unilateralContactTypeList(1) = 'None'
      MEF90DefMech_unilateralContactTypeList(2) = 'HydrostaticDeviatoric'
      MEF90DefMech_unilateralContactTypeList(3) = 'Hydrostatic'
      MEF90DefMech_unilateralContactTypeList(4) = 'Deviatoric'
      MEF90DefMech_unilateralContactTypeList(5) = 'PrincipalStrains'
      MEF90DefMech_unilateralContactTypeList(6) = 'Masonry'
      MEF90DefMech_unilateralContactTypeList(7) = 'MEF90DefMech_unilateralContactTypeList'
      MEF90DefMech_unilateralContactTypeList(8) = '_MEF90DefMech_unilateralContactTypeList'
      MEF90DefMech_unilateralContactTypeList(9) = ''

      MEF90DefMech_drivingForceTypeList(1) = 'None'
      MEF90DefMech_drivingForceTypeList(2) = 'DruckerPrager'
      MEF90DefMech_drivingForceTypeList(3) = 'DruckerPrager2'
      MEF90DefMech_drivingForceTypeList(4) = 'MEF90DefMech_drivingForceTypeList'
      MEF90DefMech_drivingForceTypeList(5) = '_MEF90DefMech_drivingForceTypeList'
      MEF90DefMech_drivingForceTypeList(6) = ''
   end subroutine MEF90DefMechCtxInitialize_Private

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCtxCreate"
!!!
!!!
!!!  MEF90DefMechCtxCreate:
!!!
!!!  (c) 2012-20 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90DefMechCtxCreate(DefMechCtx, dm, MEF90Ctx, ierr)
      type(MEF90DefMechCtx_Type), intent(OUT)                   :: DefMechCtx
      type(tDM), target, intent(IN)                             :: dm
      type(MEF90Ctx_Type), target, intent(IN)                   :: MEF90Ctx
      PetscErrorCode, intent(INOUT)                             :: ierr

      type(MEF90CtxGlobalOptions_Type), pointer                 :: MEF90CtxGlobalOptions
      type(tIS)                                                 :: setIS
      PetscInt                                                  :: set, numSet
      PetscInt, dimension(:), pointer                           :: setID
      character(len=MEF90MXSTRLEN)                              :: filename, IOBuffer
      character(len=MEF90MXSTRLEN)                              :: vecName
      type(tDM), dimension(:), pointer                          :: dmList
      type(tPetscSF)                                            :: dummySF

      PetscCall(MEF90DefMechCtxInitialize_Private(ierr))
      DefMechCtx%MEF90Ctx => MEF90Ctx
      PetscCall(PetscBagCreate(MEF90Ctx%comm, sizeofMEF90DefMechGlobalOptions, DefMechCtx%GlobalOptionsBag, ierr))
      PetscCall(DMGetLabelIdIS(dm, MEF90CellSetLabelName, setIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90Ctx%comm, setIS, ierr))
      PetscCall(ISGetLocalSize(setIS, numSet, ierr))
      allocate (DefMechCtx%CellSetOptionsBag(numSet), stat=ierr)
      do set = 1, numSet
         PetscCall(PetscBagCreate(MEF90Ctx%comm, sizeofMEF90DefMechCellSetOptions, DefMechCtx%CellSetOptionsBag(set), ierr))
      end do
      PetscCall(ISDestroy(setIS, ierr))

      PetscCall(DMGetLabelIdIS(dm, MEF90FaceSetLabelName, setIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90Ctx%comm, setIS, ierr))
      PetscCall(ISGetLocalSize(setIS, numSet, ierr))
      allocate (DefMechCtx%FaceSetOptionsBag(numSet), stat=ierr)
      do set = 1, numSet
         PetscCall(PetscBagCreate(MEF90Ctx%comm, sizeofMEF90DefMechFaceSetOptions, DefMechCtx%FaceSetOptionsBag(set), ierr))
      end do
      PetscCall(ISDestroy(setIS, ierr))

      PetscCall(DMGetLabelIdIS(dm, MEF90VertexSetLabelName, setIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90Ctx%comm, setIS, ierr))
      PetscCall(ISGetLocalSize(setIS, numSet, ierr))
      allocate (DefMechCtx%VertexSetOptionsBag(numSet), stat=ierr)
      do set = 1, numSet
         PetscCall(PetscBagCreate(MEF90Ctx%comm, sizeofMEF90DefMechVertexSetOptions, DefMechCtx%VertexSetOptionsBag(set), ierr))
      end do
      PetscCall(ISDestroy(setIS, ierr))

      !!!
      !!! Create energy viewers
      !!!
      filename = trim(MEF90FilePrefix(MEF90Ctx%resultFile))//'.ener'
      PetscCall(PetscViewerASCIIOpen(MEF90Ctx%comm, filename, DefMechCtx%globalEnergyViewer, ierr))
      PetscCall(PetscViewerASCIIPrintf(DefMechCtx%globalEnergyViewer, "# step     load            elastic energy  work            cohesive energy surface energy  total energy   plastic dissipation \n", ierr))
      PetscCall(PetscViewerFlush(DefMechCtx%globalEnergyViewer, ierr))

      PetscCall(DMGetLabelIdIS(dm, MEF90CellSetLabelName, setIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90Ctx%comm, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))

      ! PetscCall(ISGetLocalSize(setIS,numSet,ierr))
      allocate (DefMechCtx%setEnergyViewer(size(setID)), stat=ierr)
      do set = 1, size(setID)
         write (filename, 101) trim(MEF90FilePrefix(MEF90Ctx%resultFile)), setID(set)
         PetscCall(PetscViewerASCIIOpen(MEF90Ctx%comm, filename, DefMechCtx%setEnergyViewer(set), ierr))
         write (IOBuffer, 102) setID(set)
         PetscCall(PetscViewerASCIIPrintf(DefMechCtx%setEnergyViewer(set), IOBuffer, ierr))
         PetscCall(PetscViewerASCIIPrintf(DefMechCtx%setEnergyViewer(set), "# step     load            elastic energy  work            cohesive energy surface energy  total energy   plastic dissipation\n", ierr))
         PetscCall(PetscViewerFlush(DefMechCtx%setEnergyViewer(set), ierr))
      end do
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))
101   format(A, '-', I4.4, '.enerblk')
102   format("# cell set ", I4, "\n")
      DefMechCtx%analysisTime = 0.0_kr
      DefMechCtx%timeStep = 0.0_kr

      !!! Create Vecs and SF
      PetscCall(DMGetDimension(dm, DefMechCtx%dim, ierr))

      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(DefMechCtx%MEF90Ctx%GlobalOptionsBag, MEF90CtxGlobalOptions, ierr))

      vecName = "Displacement"
      allocate (DefMechCtx%displacementLocal, stat=ierr)
      PetscCall(MEF90CreateLocalVector(dm, MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, DefMechCtx%dim, vecName, DefMechCtx%displacementLocal, ierr))
      allocate (DefMechCtx%displacementPreviousStepLocal, stat=ierr)
      PetscCall(VecDuplicate(DefMechCtx%displacementLocal, DefMechCtx%displacementPreviousStepLocal, ierr))
      PetscCall(PetscObjectSetName(DefMechCtx%displacementPreviousStepLocal, "DisplacementPreviousStep", ierr))
      allocate (DefMechCtx%displacementLowerBoundLocal, stat=ierr)
      PetscCall(VecDuplicate(DefMechCtx%displacementLocal, DefMechCtx%displacementLowerBoundLocal, ierr))
      PetscCall(PetscObjectSetName(DefMechCtx%displacementLowerBoundLocal, "DisplacementLowerBound", ierr))
      allocate (DefMechCtx%displacementUpperBoundLocal, stat=ierr)
      PetscCall(VecDuplicate(DefMechCtx%displacementLocal, DefMechCtx%displacementUpperBoundLocal, ierr))
      PetscCall(PetscObjectSetName(DefMechCtx%displacementUpperBoundLocal, "DisplacementUpperBound", ierr))

      vecName = "Damage"
      allocate (DefMechCtx%damageLocal, stat=ierr)
      PetscCall(MEF90CreateLocalVector(dm, MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, 1_ki, vecName, DefMechCtx%damageLocal, ierr))
      allocate (DefMechCtx%damagePreviousStepLocal, stat=ierr)
      PetscCall(VecDuplicate(DefMechCtx%damageLocal, DefMechCtx%damagePreviousStepLocal, ierr))
      PetscCall(PetscObjectSetName(DefMechCtx%damagePreviousStepLocal, "damagePreviousStep", ierr))

      allocate (DefMechCtx%TemperatureLocal, stat=ierr)
      PetscCall(VecDuplicate(DefMechCtx%damageLocal, DefMechCtx%TemperatureLocal, ierr))
      PetscCall(PetscObjectSetName(DefMechCtx%TemperatureLocal, "Temperature", ierr))

      vecName = "cohesiveDisplacement"
      allocate (DefMechCtx%cohesiveDisplacement, stat=ierr)
      PetscCall(MEF90CreateLocalVector(dm, MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, DefMechCtx%dim, vecName, DefMechCtx%cohesiveDisplacement, ierr))
      vecName = "bodyForce"
      allocate (DefMechCtx%bodyForce, stat=ierr)
      PetscCall(MEF90CreateCellVector(dm, DefMechCtx%dim, vecName, DefMechCtx%bodyForce, ierr))
      vecName = "boundaryForce"
      allocate (DefMechCtx%boundaryForce, stat=ierr)
      PetscCall(MEF90CreateBoundaryCellVector(dm, DefMechCtx%dim, vecName, DefMechCtx%boundaryForce, ierr))
      vecName = "pressureForce"
      allocate (DefMechCtx%pressureForce, stat=ierr)
      PetscCall(MEF90CreateBoundaryCellVector(dm, 1_ki, vecName, DefMechCtx%pressureForce, ierr))

      vecName = "plasticStrain"
      allocate (DefMechCtx%plasticStrain, stat=ierr)
      PetscCall(MEF90CreateCellVector(dm, (DefMechCtx%dim * (DefMechCtx%dim + 1_ki)) / 2_ki, vecName, DefMechCtx%plasticStrain, ierr))
      allocate (DefMechCtx%cumulatedPlasticDissipation, stat=ierr)
      PetscCall(VecDuplicate(DefMechCtx%plasticStrain, DefMechCtx%cumulatedPlasticDissipation, ierr))
      PetscCall(PetscObjectSetName(DefMechCtx%cumulatedPlasticDissipation, "cumulatedPlasticDissipation", ierr))
      allocate (DefMechCtx%stress, stat=ierr)
      PetscCall(VecDuplicate(DefMechCtx%plasticStrain, DefMechCtx%stress, ierr))
      PetscCall(PetscObjectSetName(DefMechCtx%stress, "Stress", ierr))

      !!! Create megaDM
      allocate (dmList(7))
      PetscCall(VecGetDM(DefMechCtx%displacementLocal, dmList(1), ierr))
      PetscCall(VecGetDM(DefMechCtx%damageLocal, dmList(2), ierr))
      PetscCall(VecGetDM(DefMechCtx%temperatureLocal, dmList(3), ierr))
      PetscCall(VecGetDM(DefMechCtx%bodyForce, dmList(4), ierr))
      PetscCall(VecGetDM(DefMechCtx%boundaryForce, dmList(5), ierr))
      PetscCall(VecGetDM(DefMechCtx%pressureForce, dmList(6), ierr))
      PetscCall(VecGetDM(DefMechCtx%plasticStrain, dmList(7), ierr))
      PetscCall(DMCreateSuperDM(dmList, 7_ki, PETSC_NULL_IS_POINTER, DefMechCtx%megaDM, ierr))
      deallocate (dmList)

      !!! Create the IO SF for all fields
      PetscCall(MEF90IOSFCreate(MEF90Ctx, DefMechCtx%displacementLocal, DefMechCtx%displacementToIOSF, DefMechCtx%IOTodisplacementSF, ierr))
      PetscCall(MEF90IOSFCreate(MEF90Ctx, DefMechCtx%damageLocal, DefMechCtx%damageToIOSF, DefMechCtx%IOTodamageSF, ierr))

      PetscCall(MEF90IOSFCreate(MEF90Ctx, DefMechCtx%cohesiveDisplacement, DefMechCtx%cohesiveDisplacementToIOSF, DefMechCtx%IOToCohesiveDisplacementSF, ierr))

      PetscCall(MEF90IOSFCreate(MEF90Ctx, DefMechCtx%temperatureLocal, DefMechCtx%temperatureToIOSF, DefMechCtx%IOTotemperatureSF, ierr))
      PetscCall(MEF90IOSFCreate(MEF90Ctx, DefMechCtx%bodyForce, DefMechCtx%bodyForceToIOSF, DefMechCtx%IOTobodyForceSF, ierr))
      ! PetscCall(MEF90IOSFCreate(MEF90Ctx,DefMechCtx%boundaryForce,DefMechCtx%boundaryForceToIOSF,DefMechCtx%IOToboundaryForceSF,ierr))
      ! PetscCall(MEF90IOSFCreate(MEF90Ctx,DefMechCtx%pressureForce,DefMechCtx%pressureForceToIOSF,DefMechCtx%IOTopressureForceSF,ierr))
      PetscCall(MEF90IOSFCreate(MEF90Ctx, DefMechCtx%stress, DefMechCtx%stressToIOSF, DefMechCtx%IOToStressSF, ierr))
      PetscCall(MEF90IOSFCreate(MEF90Ctx, DefMechCtx%plasticStrain, DefMechCtx%plasticStrainToIOSF, DefMechCtx%IOToplasticStrainSF, ierr))
      PetscCall(MEF90IOSFCreate(MEF90Ctx, DefMechCtx%cumulatedPlasticDissipation, DefMechCtx%cumulatedPlasticDissToIOSF, DefMechCtx%IOToCumulatedPlasticDissSF, ierr))

      !!! Create the SF to exchange boundary values of the displacement and damage.
      PetscCall(MEF90ConstraintSFCreate(DefMechCtx%MEF90Ctx, DefMechCtx%displacementLocal, DefMechCtx%displacementLocal, DefMechCtx%displacementConstraintsSF, dummySF, ierr))
      PetscCall(PetscSFDestroy(dummySF, ierr))
      PetscCall(MEF90ConstraintSFCreate(DefMechCtx%MEF90Ctx, DefMechCtx%damageLocal, DefMechCtx%damageLocal, DefMechCtx%damageConstraintsSF, dummySF, ierr))
      PetscCall(PetscSFDestroy(dummySF, ierr))
   end subroutine MEF90DefMechCtxCreate

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCtxDestroy"
!!!
!!!
!!!  MEF90DefMechCtxDestroy: destroys a MEF90DefMechCtx_Type
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90DefMechCtxDestroy(DefMechCtx, ierr)
      type(MEF90DefMechCtx_Type), intent(INOUT)        :: DefMechCtx
      PetscErrorCode, intent(INOUT)                    :: ierr

      PetscInt                                        :: set

      PetscCall(PetscBagDestroy(DefMechCtx%GlobalOptionsBag, ierr))
      do set = 1, size(DefMechCtx%CellSetOptionsBag)
         PetscCall(PetscBagDestroy(DefMechCtx%CellSetOptionsBag(set), ierr))
      end do
      deallocate (DefMechCtx%CellSetOptionsBag)
      do set = 1, size(DefMechCtx%FaceSetOptionsBag)
         PetscCall(PetscBagDestroy(DefMechCtx%FaceSetOptionsBag(set), ierr))
      end do
      deallocate (DefMechCtx%FaceSetOptionsBag)
      do set = 1, size(DefMechCtx%VertexSetOptionsBag)
         PetscCall(PetscBagDestroy(DefMechCtx%VertexSetOptionsBag(set), ierr))
      end do
      deallocate (DefMechCtx%VertexSetOptionsBag)

      !!!
      !!! Close energy viewers
      !!!
      PetscCall(PetscViewerDestroy(DefMechCtx%globalEnergyViewer, ierr))
      do set = 1, size(DefMechCtx%setEnergyViewer)
         PetscCall(PetscViewerDestroy(DefMechCtx%setEnergyViewer(set), ierr))
      end do
      deallocate (DefMechCtx%setEnergyViewer)

      !!! Destroy Vecs and SF and deAllocate them
      if (associated(DefMechCtx%displacementLocal)) then
         PetscCall(VecDestroy(DefMechCtx%displacementLocal, ierr))
         deallocate (DefMechCtx%displacementLocal)
         nullify (DefMechCtx%displacementLocal)
      end if
      if (associated(DefMechCtx%displacementPreviousStepLocal)) then
         PetscCall(VecDestroy(DefMechCtx%displacementPreviousStepLocal, ierr))
         deallocate (DefMechCtx%displacementPreviousStepLocal)
         nullify (DefMechCtx%displacementPreviousStepLocal)
      end if
      if (associated(DefMechCtx%displacementLowerBoundLocal)) then
         PetscCall(VecDestroy(DefMechCtx%displacementLowerBoundLocal, ierr))
         deallocate (DefMechCtx%displacementLowerBoundLocal)
         nullify (DefMechCtx%displacementLowerBoundLocal)
      end if
      if (associated(DefMechCtx%displacementUpperBoundLocal)) then
         PetscCall(VecDestroy(DefMechCtx%displacementUpperBoundLocal, ierr))
         deallocate (DefMechCtx%displacementUpperBoundLocal)
         nullify (DefMechCtx%displacementUpperBoundLocal)
      end if

      if (associated(DefMechCtx%damageLocal)) then
         PetscCall(VecDestroy(DefMechCtx%damageLocal, ierr))
         deallocate (DefMechCtx%damageLocal)
         nullify (DefMechCtx%damageLocal)
      end if
      if (associated(DefMechCtx%damagePreviousStepLocal)) then
         PetscCall(VecDestroy(DefMechCtx%damagePreviousStepLocal, ierr))
         deallocate (DefMechCtx%damagePreviousStepLocal)
         nullify (DefMechCtx%damagePreviousStepLocal)
      end if

      if (associated(DefMechCtx%temperatureLocal)) then
         PetscCall(VecDestroy(DefMechCtx%temperatureLocal, ierr))
         deallocate (DefMechCtx%temperatureLocal)
         nullify (DefMechCtx%temperatureLocal)
      end if

      if (associated(DefMechCtx%cohesiveDisplacement)) then
         PetscCall(VecDestroy(DefMechCtx%cohesiveDisplacement, ierr))
         deallocate (DefMechCtx%cohesiveDisplacement)
         nullify (DefMechCtx%cohesiveDisplacement)
      end if

      if (associated(DefMechCtx%bodyForce)) then
         PetscCall(VecDestroy(DefMechCtx%bodyForce, ierr))
         deallocate (DefMechCtx%bodyForce)
         nullify (DefMechCtx%bodyForce)
      end if
      if (associated(DefMechCtx%boundaryForce)) then
         PetscCall(VecDestroy(DefMechCtx%boundaryForce, ierr))
         deallocate (DefMechCtx%boundaryForce)
         nullify (DefMechCtx%boundaryForce)
      end if
      if (associated(DefMechCtx%pressureForce)) then
         PetscCall(VecDestroy(DefMechCtx%pressureForce, ierr))
         deallocate (DefMechCtx%pressureForce)
         nullify (DefMechCtx%pressureForce)
      end if

      if (associated(DefMechCtx%plasticStrain)) then
         PetscCall(VecDestroy(DefMechCtx%plasticStrain, ierr))
         deallocate (DefMechCtx%plasticStrain)
         nullify (DefMechCtx%plasticStrain)
      end if
      if (associated(DefMechCtx%cumulatedPlasticDissipation)) then
         PetscCall(VecDestroy(DefMechCtx%cumulatedPlasticDissipation, ierr))
         deallocate (DefMechCtx%cumulatedPlasticDissipation)
         nullify (DefMechCtx%cumulatedPlasticDissipation)
      end if
      if (associated(DefMechCtx%Stress)) then
         PetscCall(VecDestroy(DefMechCtx%Stress, ierr))
         deallocate (DefMechCtx%Stress)
         nullify (DefMechCtx%Stress)
      end if

      !!! Destroy all PetscSF
      PetscCall(PetscSFDestroy(DefMechCtx%displacementToIOSF, ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%IOTodisplacementSF, ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%cohesiveDisplacementToIOSF, ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%IOToCohesiveDisplacementSF, ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%damageToIOSF, ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%IOTodamageSF, ierr))

      PetscCall(PetscSFDestroy(DefMechCtx%temperatureToIOSF, ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%IOTotemperatureSF, ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%bodyForceToIOSF, ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%IOTobodyForceSF, ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%boundaryForceToIOSF, ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%IOToboundaryForceSF, ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%pressureForceToIOSF, ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%IOTopressureForceSF, ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%stressToIOSF, ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%IOToStressSF, ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%plasticStrainToIOSF, ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%IOToplasticStrainSF, ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%cumulatedPlasticDissToIOSF, ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%IOToCumulatedPlasticDissSF, ierr))

      !!! Destroy the SF to exchange boundary values of the displacement and damage.
      PetscCall(PetscSFDestroy(DefMechCtx%displacementConstraintsSF, ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%damageConstraintsSF, ierr))

      !!! Destroy the megaDM
      PetscCall(DMDestroy(DefMechCtx%megaDM, ierr))
   end subroutine MEF90DefMechCtxDestroy

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90DefMechCtxGlobalOptions"
!!!
!!!
!!!  PetscBagRegisterMEF90DefMechCtxGlobalOptions:
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine PetscBagRegisterMEF90DefMechCtxGlobalOptions(bag, name, prefix, default, ierr)
      type(tPetscBag), intent(IN)                               :: bag
      character(len=*), intent(IN)                              :: prefix, name
      type(MEF90DefMechGlobalOptions_Type), intent(IN)          :: default
      PetscErrorCode, intent(INOUT)                             :: ierr

      type(MEF90DefMechGlobalOptions_Type), pointer             :: DefMechGlobalOptions

      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(bag, DefMechGlobalOptions, ierr))
      PetscCall(PetscBagSetName(bag, trim(name), "DefMechGlobalOptions MEF90 Defect Mechanics global options", ierr))
      PetscCall(PetscBagSetOptionsPrefix(bag, trim(prefix), ierr))

      PetscCall(PetscBagRegisterEnum(bag, DefMechGlobalOptions%timeSteppingType, MEF90DefMech_TimeSteppingTypeList, default%timeSteppingType, 'DefMech_TimeStepping_Type', 'Type of defect mechanics Time steping', ierr))
      PetscCall(PetscBagRegisterEnum(bag, DefMechGlobalOptions%solverType, MEF90DefMech_SolverTypeList, default%solverType, 'DefMech_solver_Type', 'Type of defect mechanics solver', ierr))
      PetscCall(PetscBagRegisterEnum(bag, DefMechGlobalOptions%damageSolverType, MEF90DefMech_DamageSolverTypeList, default%damageSolverType, 'DefMech_damageSolver_Type', 'Type of defect mechanics damage solver', ierr))

      PetscCall(PetscBagRegisterBool(bag, DefMechGlobalOptions%temperatureExport, default%temperatureExport, 'temperature_export', 'Export temperature in result file', ierr))
      PetscCall(PetscBagRegisterBool(bag, DefMechGlobalOptions%displacementExport, default%displacementExport, 'displacement_export', 'Export displacement in result file', ierr))
      PetscCall(PetscBagRegisterBool(bag, DefMechGlobalOptions%damageExport, default%damageExport, 'damage_export', 'Export damage in result file', ierr))
      PetscCall(PetscBagRegisterBool(bag, DefMechGlobalOptions%stressExport, default%stressExport, 'stress_export', 'Export stress in result file', ierr))
      PetscCall(PetscBagRegisterBool(bag, DefMechGlobalOptions%plasticStrainExport, default%plasticStrainExport, 'plasticstrain_export', 'Export plastic strain in result file', ierr))
      PetscCall(PetscBagRegisterBool(bag, DefMechGlobalOptions%cumulatedPlasticDissipationExport, default%cumulatedPlasticDissipationExport, 'cumulatedplasticdissipation_export', 'Export cumulated plastic dissipation in result file', ierr))

      PetscCall(PetscBagRegisterEnum(bag, DefMechGlobalOptions%boundaryDisplacementScaling, MEF90ScalingList, default%boundaryDisplacementScaling, 'boundaryDisplacement_scaling', 'Boundary displacement scaling', ierr))
      PetscCall(PetscBagRegisterEnum(bag, DefMechGlobalOptions%displacementLowerBoundScaling, MEF90ScalingList, default%displacementLowerBoundScaling, 'displacementlowerbound_scaling', 'Displacement lower bound scaling', ierr))
      PetscCall(PetscBagRegisterEnum(bag, DefMechGlobalOptions%displacementUpperBoundScaling, MEF90ScalingList, default%displacementUpperBoundScaling, 'displacementupperbound_scaling', 'Displacement upper bound scaling', ierr))
      PetscCall(PetscBagRegisterEnum(bag, DefMechGlobalOptions%cohesiveDisplacementScaling, MEF90ScalingList, default%cohesiveDisplacementScaling, 'cohesiveDisplacement_scaling', 'Cohesive displacement scaling', ierr))
      PetscCall(PetscBagRegisterEnum(bag, DefMechGlobalOptions%boundaryDamageScaling, MEF90ScalingList, default%boundaryDamageScaling, 'boundaryDamage_scaling', 'Boundary damage scaling', ierr))
      PetscCall(PetscBagRegisterEnum(bag, DefMechGlobalOptions%bodyForceScaling, MEF90ScalingList, default%bodyForceScaling, 'bodyforce_scaling', 'Body force scaling', ierr))
      PetscCall(PetscBagRegisterEnum(bag, DefMechGlobalOptions%boundaryForceScaling, MEF90ScalingList, default%boundaryForceScaling, 'boundaryforce_scaling', 'Boundary force scaling', ierr))
      PetscCall(PetscBagRegisterEnum(bag, DefMechGlobalOptions%pressureForceScaling, MEF90ScalingList, default%pressureforceScaling, 'pressureForce_scaling', 'Pressure force scaling', ierr))
      PetscCall(PetscBagRegisterEnum(bag, DefMechGlobalOptions%CrackPressureScaling, MEF90ScalingList, default%CrackPressureScaling, 'crackPressure_scaling', 'Crack Pressure scaling', ierr))

      PetscCall(PetscBagRegisterReal(bag, DefMechGlobalOptions%damageATol, default%damageATol, 'defmech_damage_atol', 'Absolute tolerance on damage error', ierr))
      PetscCall(PetscBagRegisterInt(bag, DefMechGlobalOptions%damageMaxIt, default%damageMaxIt, 'defmech_damage_maxit', 'Maximum number of alternate minimizations for damage', ierr))
      PetscCall(PetscBagRegisterInt(bag, DefMechGlobalOptions%PCLag, default%PCLag, 'defmech_pclag', 'Interval at which the PC is recomputed during alternate minimization', ierr))
      PetscCall(PetscBagRegisterReal(bag, DefMechGlobalOptions%SOROmega, default%SOROmega, 'defmech_SOR_Omega', 'Alterate Minimization over relaxation factor (>0 for limited, <0 for projected) ', ierr))
      PetscCall(PetscBagRegisterReal(bag, DefMechGlobalOptions%irrevthres, default%irrevthres, 'defmech_irrevThres', 'Threshold above which irreversibility is enforced (0 for monotonicity, .99 for equality)', ierr))

      PetscCall(PetscBagRegisterEnum(bag, DefMechGlobalOptions%BTType, MEF90DefMech_BTTypeList, default%BTType, 'BT_Type', 'Backtracking type', ierr))
      PetscCall(PetscBagRegisterInt(bag, DefMechGlobalOptions%BTInterval, default%BTInterval, 'BT_Interval', 'Interval at which Backtracking is run in inner loop (0 for outer loop)', ierr))
      PetscCall(PetscBagRegisterInt(bag, DefMechGlobalOptions%BTScope, default%BTScope, 'BT_Scope', 'Backtracking scope (0 for unlimited)', ierr))
      PetscCall(PetscBagRegisterReal(bag, DefMechGlobalOptions%BTTol, default%BTTol, 'BT_Tol', 'Backtracking relative tolerance', ierr))

      PetscCall(PetscBagRegisterReal(bag, DefMechGlobalOptions%plasticStrainATol, default%plasticStrainATol, 'defmech_plasticstrain_atol', 'Absolute tolerance on plastic error', ierr))

      PetscCall(PetscBagRegisterReal(bag, DefMechGlobalOptions%InjectedVolumeATol, default%InjectedVolumeATol, 'defmech_InjectedVolume_atol', 'Absolute tolerance on injected volume error', ierr))
      PetscCall(PetscBagRegisterReal(bag, DefMechGlobalOptions%dampingCoefficientDisplacement, default%dampingCoefficientDisplacement, 'defmech_dampingCoefficient_displacement', 'Damping coefficient on displacement field (0 for minimization, 1 for semi-implicit gradient flow)', ierr))
      PetscCall(PetscBagRegisterReal(bag, DefMechGlobalOptions%dampingCoefficientDamage, default%dampingCoefficientDamage, 'defmech_dampingCoefficient_damage', 'Damping coefficient on damage field (0 for minimization, 1 for semi-implicit gradient flow)', ierr))
   end subroutine PetscBagRegisterMEF90DefMechCtxGlobalOptions

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90DefMechCtxCellSetOptions"
!!!
!!!
!!!  PetscBagRegisterMEF90DefMechCtxCellSetOptions:
!!!
!!!  (c) 2012-18 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine PetscBagRegisterMEF90DefMechCtxCellSetOptions(bag, name, prefix, default, ierr)
      type(tPetscBag), intent(IN)                         :: bag
      character(len=*), intent(IN)                        :: prefix, name
      type(MEF90DefMechCellSetOptions_Type), intent(IN)   :: default
      PetscErrorCode, intent(INOUT)                       :: ierr

      type(MEF90DefMechCellSetOptions_Type), pointer      :: DefMechCellSetOptions

      PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(bag, DefMechCellSetOptions, ierr))
      PetscCall(PetscBagSetName(bag, trim(name), "DefMechCellSetOptions MEF90 Defect Mechanics Cell Set options", ierr))
      PetscCall(PetscBagSetOptionsPrefix(bag, trim(prefix), ierr))

      DefMechCellSetOptions%bodyForce = default%bodyForce
      DefMechCellSetOptions%boundaryDisplacement = default%boundaryDisplacement
      DefMechCellSetOptions%displacementLowerBound = default%displacementLowerBound
      DefMechCellSetOptions%displacementUpperBound = default%displacementUpperBound
      DefMechCellSetOptions%Has_displacementBC = default%Has_displacementBC

      PetscCall(PetscBagRegisterRealArray(bag, DefMechCellSetOptions%bodyForce, 3_ki, 'bodyForce', '[N.m^(-3) / N.m^(-2)] (f): body force', ierr))
      PetscCall(PetscBagRegisterReal(bag, DefMechCellSetOptions%CrackPressure, default%CrackPressure, 'CrackPressure', 'without unit: internal crack pressure', ierr))
      PetscCall(PetscBagRegisterReal(bag, DefMechCellSetOptions%DamageATLinSoftk, default%DamageATLinSoftk, 'damage_LinSoft_k', '[unit-less] (k): k parameter in the Linear Softening damage model', ierr))
      PetscCall(PetscBagRegisterReal(bag, DefMechCellSetOptions%DamageAT1expb, default%DamageAT1expb, 'damage_AT1exp_b', '[unit-less] (b): b parameter in tha AT1 model with exponential stiffness interpolation', ierr))
      PetscCall(PetscBagRegisterEnum(bag, DefMechCellSetOptions%damageType, MEF90DefMech_damageTypeList, default%damageType, 'damage_type', 'Type of damage law', ierr))
      PetscCall(PetscBagRegisterEnum(bag, DefMechCellSetOptions%plasticityType, MEF90DefMech_plasticityTypeList, default%plasticityType, 'plasticity_type', 'Type of plasticity law', ierr))
      PetscCall(PetscBagRegisterEnum(bag, DefMechCellSetOptions%unilateralContactType, MEF90DefMech_unilateralContactTypeList, default%unilateralContactType, 'unilateralContact_type', 'Type of handling of unilateral contact', ierr))
      PetscCall(PetscBagRegisterReal(bag, DefMechCellSetOptions%unilateralContactHydrostaticDeviatoricGamma, default%unilateralContactHydrostaticDeviatoricGamma, 'unilateralContact_hydrostaticDeviatoric_gamma', '[unit-less] (gamma): Hydrostatic Deviatoric regularization parameter', ierr))
      PetscCall(PetscBagRegisterBool(bag, DefMechCellSetOptions%unilateralContactHybrid, default%unilateralContactHybrid, 'unilateralContact_hybrid', 'Use hybrid unilateral contact formulation (Y/N)', ierr))
      PetscCall(PetscBagRegisterEnum(bag, DefMechCellSetOptions%drivingForceType, MEF90DefMech_drivingForceTypeList, default%drivingForceType, 'drivingForce_type', 'Type of nucleation driving force', ierr))
      PetscCall(PetscBagRegisterRealArray(bag, DefMechCellSetOptions%cohesiveDisplacement, 3_ki, 'cohesiveDisplacement', '[m] (U): Cohesive displacement value', ierr))
      PetscCall(PetscBagRegisterBoolArray(bag, DefMechCellSetOptions%Has_displacementBC, 3_ki, 'DisplacementBC', 'Displacement has Dirichlet boundary Condition (Y/N)', ierr))
      PetscCall(PetscBagRegisterRealArray(bag, DefMechCellSetOptions%boundaryDisplacement, 3_ki, 'boundaryDisplacement', '[m] (U): Displacement boundary value', ierr))
      PetscCall(PetscBagRegisterRealArray(bag, DefMechCellSetOptions%displacementLowerBound, 3_ki, 'displacementLowerBound', '[m] (U): Displacement lower bound', ierr))
      PetscCall(PetscBagRegisterRealArray(bag, DefMechCellSetOptions%displacementUpperBound, 3_ki, 'displacementUpperBound', '[m] (U): Displacement upper bound', ierr))
      PetscCall(PetscBagRegisterBool(bag, DefMechCellSetOptions%Has_DamageBC, default%Has_DamageBC, 'DamageBC', 'Damage has Dirichlet boundary Condition (Y/N)', ierr))
      PetscCall(PetscBagRegisterBool(bag, DefMechCellSetOptions%CrackVolumeControlled, default%CrackVolumeControlled, 'CrackVolumeControlled', 'Crack Pressure controlled by the crack volume in this block (Y/N)', ierr))
      PetscCall(PetscBagRegisterBool(bag, DefMechCellSetOptions%WorkControlled, default%WorkControlled, 'WorkControlled', 'Force magnitude controlled by its work in this block (Y/N)', ierr))
      PetscCall(PetscBagRegisterReal(bag, DefMechCellSetOptions%boundaryDamage, default%boundaryDamage, 'boundaryDamage', '[unit-less] (alpha): Damage boundary value', ierr))
   end subroutine PetscBagRegisterMEF90DefMechCtxCellSetOptions

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90DefMechCtxFaceSetOptions"
!!!
!!!
!!!  PetscBagRegisterMEF90DefMechCtxFaceSetOptions:
!!!
!!!  (c) 2012-18 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine PetscBagRegisterMEF90DefMechCtxFaceSetOptions(bag, name, prefix, default, ierr)
      type(tPetscBag), intent(IN)                         :: bag
      character(len=*), intent(IN)                        :: prefix, name
      type(MEF90DefMechFaceSetOptions_Type), intent(IN)   :: default
      PetscErrorCode, intent(INOUT)                       :: ierr

      type(MEF90DefMechFaceSetOptions_Type), pointer      :: DefMechFaceSetOptions

      PetscCall(PetscBagGetDataMEF90DefMechCtxFaceSetOptions(bag, DefMechFaceSetOptions, ierr))
      PetscCall(PetscBagSetName(bag, trim(name), "DefMechFaceSetOptions MEF90 Defect Mechanics Face Set options", ierr))
      PetscCall(PetscBagSetOptionsPrefix(bag, trim(prefix), ierr))

      DefMechFaceSetOptions%boundaryForce = default%boundaryForce
      DefMechFaceSetOptions%boundaryDisplacement = default%boundaryDisplacement
      DefMechFaceSetOptions%displacementLowerBound = default%displacementLowerBound
      DefMechFaceSetOptions%displacementUpperBound = default%displacementUpperBound
      DefMechFaceSetOptions%Has_displacementBC = default%Has_displacementBC

      PetscCall(PetscBagRegisterRealArray(bag, DefMechFaceSetOptions%boundaryForce, 3_ki, 'boundaryForce', '[N.m^(-2) / N.m^(-1)] (f): boundary force', ierr))
      PetscCall(PetscBagRegisterReal(bag, DefMechFaceSetOptions%pressureForce, default%pressureForce, 'pressureForce', 'without unit: internal crack pressure', ierr))
      PetscCall(PetscBagRegisterBoolArray(bag, DefMechFaceSetOptions%Has_displacementBC, 3_ki, 'DisplacementBC', 'Displacement has Dirichlet boundary Condition (Y/N)', ierr))
      PetscCall(PetscBagRegisterRealArray(bag, DefMechFaceSetOptions%boundaryDisplacement, 3_ki, 'boundaryDisplacement', '[m] (U): Displacement boundary value', ierr))
      PetscCall(PetscBagRegisterRealArray(bag, DefMechFaceSetOptions%displacementLowerBound, 3_ki, 'displacementLowerBound', '[m] (U): Displacement lower bound', ierr))
      PetscCall(PetscBagRegisterRealArray(bag, DefMechFaceSetOptions%displacementUpperBound, 3_ki, 'displacementUpperBound', '[m] (U): Displacement upper bound', ierr))
      PetscCall(PetscBagRegisterBool(bag, DefMechFaceSetOptions%Has_DamageBC, default%Has_DamageBC, 'DamageBC', 'Damage has Dirichlet boundary Condition (Y/N)', ierr))
      PetscCall(PetscBagRegisterReal(bag, DefMechFaceSetOptions%boundaryDamage, default%boundaryDamage, 'boundaryDamage', '[unit-less] (alpha): Damage boundary value', ierr))
   end subroutine PetscBagRegisterMEF90DefMechCtxFaceSetOptions

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90DefMechCtxVertexSetOptions"
!!!
!!!
!!!  PetscBagRegisterMEF90DefMechCtxVertexSetOptions:
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine PetscBagRegisterMEF90DefMechCtxVertexSetOptions(bag, name, prefix, default, ierr)
      type(tPetscBag), intent(IN)                            :: bag
      character(len=*), intent(IN)                           :: prefix, name
      type(MEF90DefMechVertexSetOptions_Type), intent(IN)    :: default
      PetscErrorCode, intent(INOUT)                          :: ierr

      type(MEF90DefMechVertexSetOptions_Type), pointer       :: DefMechVertexSetOptions
      PetscCall(PetscBagGetDataMEF90DefMechCtxVertexSetOptions(bag, DefMechVertexSetOptions, ierr))
      PetscCall(PetscBagSetName(bag, trim(name), "DefMechVertexSetOptions MEF90 Defect Mechanics Vertex Set options", ierr))
      PetscCall(PetscBagSetOptionsPrefix(bag, trim(prefix), ierr))

      DefMechVertexSetOptions%Has_displacementBC = default%Has_displacementBC
      DefMechVertexSetOptions%boundaryDisplacement = default%boundaryDisplacement
      DefMechVertexSetOptions%displacementLowerBound = default%displacementLowerBound
      DefMechVertexSetOptions%displacementUpperBound = default%displacementUpperBound
      PetscCall(PetscBagRegisterBoolArray(bag, DefMechVertexSetOptions%Has_displacementBC, 3_ki, 'DisplacementBC', 'Displacement has Dirichlet boundary Condition (Y/N)', ierr))
      PetscCall(PetscBagRegisterRealArray(bag, DefMechVertexSetOptions%boundaryDisplacement, 3_ki, 'boundaryDisplacement', '[m] (U): Displacement boundary value', ierr))
      PetscCall(PetscBagRegisterRealArray(bag, DefMechVertexSetOptions%displacementLowerBound, 3_ki, 'displacementLowerBound', '[m] (U): Displacement lower bound', ierr))
      PetscCall(PetscBagRegisterRealArray(bag, DefMechVertexSetOptions%displacementUpperBound, 3_ki, 'displacementUpperBound', '[m] (U): Displacement upper bound', ierr))
      PetscCall(PetscBagRegisterBool(bag, DefMechVertexSetOptions%Has_DamageBC, default%Has_DamageBC, 'DamageBC', 'Damage has Dirichlet boundary Condition (Y/N)', ierr))
      PetscCall(PetscBagRegisterReal(bag, DefMechVertexSetOptions%boundaryDamage, default%boundaryDamage, 'boundaryDamage', '[unit-less] (alpha): boundaryDamage', ierr))
   end subroutine PetscBagRegisterMEF90DefMechCtxVertexSetOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCtxSetFromOptions"
!!!
!!!
!!!  MEF90DefMechCtxSetFromOptions:
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90DefMechCtxSetFromOptions(DefMechCtx, prefix, defaultGlobalOptions, &
                                            defaultCellSetOptions, &
                                            defaultFaceSetOptions, &
                                            defaultVertexSetOptions, ierr)
      type(MEF90DefMechCtx_Type), intent(INOUT)              :: DefMechCtx
      character(len=*), intent(IN)                           :: prefix
      type(MEF90DefMechGlobalOptions_Type), intent(IN)       :: defaultGlobalOptions
      type(MEF90DefMechCellSetOptions_Type), intent(IN)      :: defaultCellSetOptions
      type(MEF90DefMechFaceSetOptions_Type), intent(IN)      :: defaultFaceSetOptions
      type(MEF90DefMechVertexSetOptions_Type), intent(IN)    :: defaultVertexSetOptions
      PetscErrorCode, intent(INOUT)                          :: ierr

      type(MEF90CtxGlobalOptions_Type), pointer              :: MEF90CtxGlobalOptions
      type(tIS)                                             :: setIS
      PetscInt, dimension(:), pointer                         :: setID
      PetscInt                                              :: set
      character(len=MEF90MXSTRLEN)                          :: IOBuffer, setName, setprefix

      type(MEF90DefMechCellSetOptions_Type), pointer         :: cellSetOptions
      type(MEF90DefMechFaceSetOptions_Type), pointer         :: faceSetOptions

      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(DefMechCtx%MEF90Ctx%GlobalOptionsBag, MEF90CtxGlobalOptions, ierr))
      !!!
      !!! Registering Global Context
      !!!
      PetscCall(PetscBagRegisterMEF90DefMechCtxGlobalOptions(DefMechCtx%GlobalOptionsBag, "MEF90DefMech Global Ctx", prefix, defaultGlobalOptions, ierr))

      if (MEF90CtxGlobalOptions%verbose > 0) then
         PetscCall(PetscBagView(DefMechCtx%GlobalOptionsBag, PETSC_VIEWER_STDOUT_WORLD, ierr))
      end if

      !!!
      !!! Registering Cell Set Context
      !!!
      PetscCall(DMGetLabelIdIS(DefMechCtx%megaDM, MEF90CellSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(DefMechCtx%MEF90Ctx%comm, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         write (setName, "('Cell set ',I4)") setID(set)
         write (setprefix, "('cs',I4.4,'_')") setID(set)
         PetscCall(PetscBagRegisterMEF90DefMechCtxCellSetOptions(DefMechCtx%CellSetOptionsBag(set), setName, setPrefix, defaultCellSetOptions, ierr))
         if (MEF90CtxGlobalOptions%verbose > 0) then
            write (IOBuffer, "('\nRegistering cell set ',I4,' prefix: ',A,'\n')") setID(set), trim(setprefix)
            PetscCall(PetscPrintf(DefMechCtx%MEF90Ctx%comm, IOBuffer, ierr))
            PetscCall(PetscBagView(DefMechCtx%CellSetOptionsBag(set), PETSC_VIEWER_STDOUT_WORLD, ierr))
            PetscCall(PetscPrintf(DefMechCtx%MEF90Ctx%comm, "\n", ierr))
         end if
      end do
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))

      !!!
      !!! Registering Face Set Context
      !!!
      PetscCall(DMGetLabelIdIS(DefMechCtx%megaDM, MEF90FaceSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(DefMechCtx%MEF90Ctx%comm, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         write (setName, "('Face set ',I4)") setID(set)
         write (setprefix, "('fs',I4.4,'_')") setID(set)
         PetscCall(PetscBagRegisterMEF90DefMechCtxFaceSetOptions(DefMechCtx%FaceSetOptionsBag(set), setName, setPrefix, defaultFaceSetOptions, ierr))
         if (MEF90CtxGlobalOptions%verbose > 0) then
            write (IOBuffer, "('\nRegistering face set ',I4,' prefix: ',A,'\n')") setID(set), trim(setprefix)
            PetscCall(PetscPrintf(DefMechCtx%MEF90Ctx%comm, IOBuffer, ierr))
            PetscCall(PetscBagView(DefMechCtx%FaceSetOptionsBag(set), PETSC_VIEWER_STDOUT_WORLD, ierr))
            PetscCall(PetscPrintf(DefMechCtx%MEF90Ctx%comm, "\n", ierr))
         end if
      end do
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))

      !!!
      !!! Registering Vertex Set Context
      !!!
      PetscCall(DMGetLabelIdIS(DefMechCtx%megaDM, MEF90VertexSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(DefMechCtx%MEF90Ctx%comm, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         write (setName, "('Vertex set ',I4)") setID(set)
         write (setprefix, "('vs',I4.4,'_')") setID(set)
         PetscCall(PetscBagRegisterMEF90DefMechCtxVertexSetOptions(DefMechCtx%VertexSetOptionsBag(set), setName, setPrefix, defaultVertexSetOptions, ierr))
         if (MEF90CtxGlobalOptions%verbose > 0) then
            write (IOBuffer, "('\nRegistering vertex set ',I4,' prefix: ',A,'\n')") setID(set), trim(setprefix)
            PetscCall(PetscPrintf(DefMechCtx%MEF90Ctx%comm, IOBuffer, ierr))
            PetscCall(PetscBagView(DefMechCtx%VertexSetOptionsBag(set), PETSC_VIEWER_STDOUT_WORLD, ierr))
            PetscCall(PetscPrintf(DefMechCtx%MEF90Ctx%comm, "\n", ierr))
         end if
      end do
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))

      !!! Identify if any of the displacement bounds are set
      DefMechCtx%hasDisplacementBounds = PETSC_FALSE

      PetscCall(DMGetLabelIdIS(DefMechCtx%megaDM, MEF90CellSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(DefMechCtx%MEF90Ctx%comm, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(DefMechCtx%CellSetOptionsBag(set), cellSetOptions, ierr))
         DefMechCtx%hasDisplacementBounds = any(cellSetOptions%displacementLowerBound /= MEF90NINFINITY) .or. DefMechCtx%hasDisplacementBounds
         DefMechCtx%hasDisplacementBounds = any(cellSetOptions%displacementUpperBound /= MEF90INFINITY) .or. DefMechCtx%hasDisplacementBounds
      end do
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))
      PetscCall(DMGetLabelIdIS(DefMechCtx%megaDM, MEF90FaceSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(DefMechCtx%MEF90Ctx%comm, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         PetscCall(PetscBagGetDataMEF90DefMechCtxFaceSetOptions(DefMechCtx%FaceSetOptionsBag(set), faceSetOptions, ierr))
         DefMechCtx%hasDisplacementBounds = any(faceSetOptions%displacementLowerBound /= MEF90NINFINITY) .or. DefMechCtx%hasDisplacementBounds
         DefMechCtx%hasDisplacementBounds = any(faceSetOptions%displacementUpperBound /= MEF90INFINITY) .or. DefMechCtx%hasDisplacementBounds
      end do
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end subroutine MEF90DefMechCtxSetFromOptions
end module m_MEF90_DefMechCtx
