Module m_MEF90_DefMechCtx_Type
#include "petsc/finclude/petsc.h"
   Use m_MEF90
   Implicit none
   
   Type MEF90DefMechCtx_Type
      Type(MEF90Ctx_Type),pointer             :: MEF90Ctx
      Type(tDM)                               :: megaDM
      PetscInt                                :: dim
      PetscReal                               :: analysisTime,timeStep

      !!!  vertex based vec
      Type(tVec),pointer                      :: displacementLocal,displacementPreviousStepLocal
      Type(tVec),pointer                      :: damageLocal,damagePreviousStepLocal
      Type(tVec),pointer                      :: displacementLowerBoundLocal,displacementUpperBoundLocal
      Type(tVec),Pointer                      :: temperatureLocal
      
      !!! cell based vec      
      Type(tVec),pointer                      :: bodyForce
      Type(tVec),pointer                      :: boundaryForce
      Type(tVec),pointer                      :: pressureForce
      !Type(tVec),pointer                      :: crackPressure
      Type(tVec),pointer                      :: cohesiveDisplacement
      Type(tVec),Pointer                      :: plasticStrain
      Type(tVec),Pointer                      :: cumulatedPlasticDissipation
      Type(tVec),Pointer                      :: stress

      Type(tPetscSF)                          :: displacementToIOSF,IOToDisplacementSF
      Type(tPetscSF)                          :: cohesiveDisplacementToIOSF,IOToCohesiveDisplacementSF
      Type(tPetscSF)                          :: displacementConstraintsSF
      Type(tPetscSF)                          :: damageToIOSF,IOToDamageSF
      Type(tPetscSF)                          :: damageConstraintsSF
      Type(tPetscSF)                          :: temperatureToIOSF,IOToTemperatureSF
      Type(tPetscSF)                          :: bodyForceToIOSF,IOToBodyForceSF
      Type(tPetscSF)                          :: boundaryForceToIOSF,IOToBoundaryForceSF
      Type(tPetscSF)                          :: pressureForceToIOSF,IOToPressureForceSF
      Type(tPetscSF)                          :: stressToIOSF,IOToStressSF
      Type(tPetscSF)                          :: plasticStrainToIOSF,IOToPlasticStrainSF
      Type(tPetscSF)                          :: cumulatedPlasticDissToIOSF,IOToCumulatedPlasticDissSF
      
      PetscBag                                :: GlobalOptionsBag
      PetscBag,Dimension(:),Pointer           :: CellSetOptionsBag
      PetscBag,Dimension(:),Pointer           :: FaceSetOptionsBag
      PetscBag,Dimension(:),Pointer           :: VertexSetOptionsBag
      PetscBag,Dimension(:),Pointer           :: MaterialPropertiesBag

      Type(tPetscViewer)                      :: globalEnergyViewer
      Type(tPetscViewer),Dimension(:),Pointer :: setEnergyViewer

      PetscBool                               :: hasDisplacementBounds
      PetscBool                               :: hasUnilateralContact
   End Type MEF90DefMechCtx_Type
   
   Type MEF90DefMechGlobalOptions_Type
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
   End Type MEF90DefMechGlobalOptions_Type

   Type MEF90DefMechCellSetOptions_Type
      PetscReal,Dimension(3)                 :: bodyforce
      PetscReal                              :: crackPressure
      PetscEnum                              :: damageType
      PetscEnum                              :: plasticityType
      PetscEnum                              :: unilateralContactType
      PetscReal                              :: unilateralContactHydrostaticDeviatoricGamma
      PetscBool                              :: unilateralContactHybrid
      PetscReal                              :: DamageATLinSoftk
      PetscReal                              :: DamageAT1expb
      PetscEnum                              :: drivingForceType
      PetscReal,Dimension(3)                 :: cohesiveDisplacement
      PetscBool,Dimension(3)                 :: Has_displacementBC
      PetscReal,Dimension(3)                 :: boundaryDisplacement
      PetscReal,Dimension(3)                 :: displacementLowerBound
      PetscReal,Dimension(3)                 :: displacementUpperBound
      PetscBool                              :: Has_damageBC
      PetscReal                              :: boundaryDamage
      PetscBool                              :: CrackVolumeControlled
      PetscBool                              :: WorkControlled
   End Type MEF90DefMechCellSetOptions_Type
   Type MEF90DefMechFaceSetOptions_Type
      PetscReal,Dimension(3)                 :: boundaryforce
      PetscReal                              :: pressureForce
      PetscBool,Dimension(3)                 :: Has_displacementBC
      PetscReal,Dimension(3)                 :: boundaryDisplacement
      PetscReal,Dimension(3)                 :: displacementLowerBound
      PetscReal,Dimension(3)                 :: displacementUpperBound
      PetscBool                              :: Has_damageBC
      PetscReal                              :: boundaryDamage
   End Type MEF90DefMechFaceSetOptions_Type

   Type MEF90DefMechVertexSetOptions_Type
      PetscBool,Dimension(3)                 :: Has_displacementBC
      PetscReal,Dimension(3)                 :: boundaryDisplacement
      PetscReal,Dimension(3)                 :: displacementLowerBound
      PetscReal,Dimension(3)                 :: displacementUpperBound
      PetscBool                              :: Has_damageBC
      PetscReal                              :: boundaryDamage
   End Type MEF90DefMechVertexSetOptions_Type 
End Module m_MEF90_DefMechCtx_Type

Module m_MEF90DefMechGlobalOptions_Private
#include "petsc/finclude/petsc.h"
   Use m_MEF90
   Use m_MEF90_DefMechCtx_Type
   Implicit None

   Private
   Public :: PetscBagGetDataMEF90DefMechCtxGlobalOptions
   
   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF90_DefMechCtx_Type
         PetscBag                                              :: bag
         Type(MEF90DefMechGlobalOptions_Type),pointer          :: data
         PetscErrorCode                                        :: ierr
      End subroutine PetscBagGetData
   End interface
Contains

#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90DefMechCtxGlobalOptions"
!!!
!!!  PetscBagGetDataMEF90DefMechCtxGlobalOptions - Custom interface to PetscGetData
!!!

   Subroutine PetscBagGetDataMEF90DefMechCtxGlobalOptions(bag,data,ierr)
      PetscBag                                              :: bag
      Type(MEF90DefMechGlobalOptions_Type),pointer          :: data
      PetscErrorCode                                        :: ierr
      
      PetscCall(PetscBagGetData(bag,data,ierr))
   End Subroutine PetscBagGetDataMEF90DefMechCtxGlobalOptions
End Module m_MEF90DefMechGlobalOptions_Private

Module m_MEF90DefMechCellSetOptions_Private
#include "petsc/finclude/petsc.h"
   Use m_MEF90
   Use m_MEF90_DefMechCtx_Type
   Implicit None

   Private
   Public :: PetscBagGetDataMEF90DefMechCtxCellSetOptions
   
   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF90_DefMechCtx_Type
         PetscBag                                              :: bag
         Type(MEF90DefMechCellSetOptions_Type),pointer         :: data
         PetscErrorCode                                        :: ierr
      End subroutine PetscBagGetData
   End interface
Contains

#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90DefMechCtxCellSetOptions"
!!!
!!!  PetscBagGetDataMEF90DefMechCtxCellSetOptions - Custom interface to PetscGetData
!!!

   Subroutine PetscBagGetDataMEF90DefMechCtxCellSetOptions(bag,data,ierr)
      PetscBag                                              :: bag
      Type(MEF90DefMechCellSetOptions_Type),pointer         :: data
      PetscErrorCode                                        :: ierr
      
      PetscCall(PetscBagGetData(bag,data,ierr))
   End Subroutine PetscBagGetDataMEF90DefMechCtxCellSetOptions
End Module m_MEF90DefMechCellSetOptions_Private

Module m_MEF90DefMechFaceSetOptions_Private
#include "petsc/finclude/petsc.h"
   Use m_MEF90
   Use m_MEF90_DefMechCtx_Type
   Implicit None

   Private
   Public :: PetscBagGetDataMEF90DefMechCtxFaceSetOptions
   
   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF90_DefMechCtx_Type
         PetscBag                                              :: bag
         Type(MEF90DefMechFaceSetOptions_Type),pointer         :: data
         PetscErrorCode                                        :: ierr
      End subroutine PetscBagGetData
   End interface
Contains
   
#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90DefMechCtxFaceSetOptions"
!!!
!!!  PetscBagGetDataMEF90DefMechCtxFaceSetOptions - Custom interface to PetscGetData
!!!

   Subroutine PetscBagGetDataMEF90DefMechCtxFaceSetOptions(bag,data,ierr)
      PetscBag                                              :: bag
      Type(MEF90DefMechFaceSetOptions_Type),pointer         :: data
      PetscErrorCode                                        :: ierr
      
      PetscCall(PetscBagGetData(bag,data,ierr))
   End Subroutine PetscBagGetDataMEF90DefMechCtxFaceSetOptions
End Module m_MEF90DefMechFaceSetOptions_Private
   
Module m_MEF90DefMechVertexSetOptions_Private
#include "petsc/finclude/petsc.h"
   Use m_MEF90
   Use m_MEF90_DefMechCtx_Type
   Implicit None

   Private
   Public :: PetscBagGetDataMEF90DefMechCtxVertexSetOptions
   
   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF90_DefMechCtx_Type
         PetscBag                                                 :: bag
         Type(MEF90DefMechVertexSetOptions_Type),pointer          :: data
         PetscErrorCode                                           :: ierr
      End subroutine PetscBagGetData
   End interface

Contains
#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90DefMechCtxVertexSetOptions"
!!!
!!!  PetscBagGetDataMEF90DefMechCtxVertexSetOptions - Custom interface to PetscGetData
!!!

   Subroutine PetscBagGetDataMEF90DefMechCtxVertexSetOptions(bag,data,ierr)
      PetscBag                                                 :: bag
      Type(MEF90DefMechVertexSetOptions_Type),pointer          :: data
      PetscErrorCode                                           :: ierr
      
      PetscCall(PetscBagGetData(bag,data,ierr))
   End Subroutine PetscBagGetDataMEF90DefMechCtxVertexSetOptions
End Module m_MEF90DefMechVertexSetOptions_Private

Module m_MEF90_DefMechCtx
#include "petsc/finclude/petsc.h"
   Use m_MEF90
   Use m_MEF90_DefMechCtx_Type
   Use m_MEF90DefMechGlobalOptions_Private
   Use m_MEF90DefMechCellSetOptions_Private
   Use m_MEF90DefMechFaceSetOptions_Private
   Use m_MEF90DefMechVertexSetOptions_Private
   Implicit none

   PetscSizeT,protected   :: sizeofMEF90DefMechGlobalOptions
   PetscSizeT,protected   :: sizeofMEF90DefMechCellSetOptions
   PetscSizeT,protected   :: sizeofMEF90DefMechFaceSetOptions
   PetscSizeT,protected   :: sizeofMEF90DefMechVertexSetOptions
   
   Enum,bind(c)
      enumerator  :: MEF90DefMech_SolverTypeAltMin = 0, &
                     MEF90DefMech_SolverTypeQuasiNewton1, &
                     MEF90DefMech_SolverTypeQuasiNewton2
   End Enum
   Character(len = MEF90MXSTRLEN),Dimension(6),protected   :: MEF90DefMech_SolverTypeList

   Enum,bind(c)
      enumerator  :: MEF90DefMech_TimeSteppingTypeNULL = 0,     &
                     MEF90DefMech_TimeSteppingTypeQuasiStatic
   End Enum
   Character(len = MEF90MXSTRLEN),Dimension(5),protected   :: MEF90DefMech_TimeSteppingTypeList
   
   Enum,bind(c)
      enumerator  :: MEF90DefMech_DamageSolverTypeSNES = 0, &
                     MEF90DefMech_DamageSolverTypeTAO
   End Enum
   Character(len = MEF90MXSTRLEN),Dimension(5),protected   :: MEF90DefMech_DamageSolverTypeList

   Enum,bind(c)
      enumerator  :: MEF90DefMech_BTTypeNULL = 0,    &
                     MEF90DefMech_BTTypeBackward,    &
                     MEF90DefMech_BTTypeForward
   End Enum
   Character(len = MEF90MXSTRLEN),Dimension(6),protected   :: MEF90DefMech_BTTypeList
   
   Enum,bind(c)
      enumerator  :: MEF90DefMech_damageTypeAT1 = 0,        &
                     MEF90DefMech_damageTypeAT1exp,         &
                     MEF90DefMech_damageTypeAT2,            &
                     MEF90DefMech_damageTypeLinSoft,        &
                     MEF90DefMech_damageTypeKKL,            &
                     MEf90DefMech_damageTypeAT1Elastic,     &
                     MEf90DefMech_damageTypeAT1expElastic,  &
                     MEf90DefMech_damageTypeAT2Elastic,     &
                     MEF90DefMech_damageTypeLinSoftElastic, &
                     MEF90DefMech_damageTypeKKLElastic

   End Enum
   Character(len = MEF90MXSTRLEN),Dimension(13),protected   :: MEF90DefMech_damageTypeList
   
   Enum,bind(c)
      enumerator  :: MEF90DefMech_plasticityTypeNone = 0,               &
                     MEF90DefMech_plasticityTypeTresca,                 &
                     MEF90DefMech_plasticityTypeVonMises,               &
                     MEF90DefMech_plasticityTypeVonMisesPlaneTheory,    &
                     MEF90DefMech_plasticityTypeCapModel,               &
                     MEF90DefMech_plasticityTypeDruckerPragerCapModel,  &
                     MEF90DefMech_plasticityTypeVonMises1D,             &
                     MEF90DefMech_plasticityTypeHillPlaneTheory,        &
                     MEF90DefMech_PlasticityTypeGreen,                  &
                     MEF90DefMech_PlasticityTypeGurson
   End Enum
   Character(len = MEF90MXSTRLEN),Dimension(13),protected   :: MEF90DefMech_plasticityTypeList
   
   Enum,bind(c)
      enumerator  :: MEF90DefMech_unilateralContactTypeNone = 0,                     &
                     MEF90DefMech_unilateralContactTypeHydrostaticDeviatoric,        &
                     MEF90DefMech_unilateralContactTypeHydrostatic,                  &
                     MEF90DefMech_unilateralContactTypeDeviatoric,                   &
                     MEF90DefMech_unilateralContactTypePrincipalStrains,             &
                     MEF90DefMech_unilateralContactTypeMasonry
   End Enum
   Character(len = MEF90MXSTRLEN),Dimension(9),protected   :: MEF90DefMech_unilateralContactTypeList

   Enum,bind(c)
      enumerator  :: MEF90DefMech_drivingForceTypeNone = 0,       &
                     MEF90_DefMechDrivingForceTypeDruckerPrager,  &
                     MEF90DefMech_drivingForceTypeDruckerPrager2
   End Enum
   Character(len = MEF90MXSTRLEN),Dimension(6),protected   :: MEF90DefMech_drivingForceTypeList
Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCtxInitialize_Private"
!!!
!!!  
!!!  MEF90DefMechCtxInitialize_Private:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90DefMechCtxInitialize_Private(ierr)
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(MEF90DefMechGlobalOptions_Type)               :: DefMechGlobalOptions
      Type(MEF90DefMechCellSetOptions_Type)              :: DefMechCellSetOptions
      Type(MEF90DefMechFaceSetOptions_Type)              :: DefMechFaceSetOptions
      Type(MEF90DefMechVertexSetOptions_Type)            :: DefMechVertexSetOptions
      character(len=1),pointer                           :: dummychar(:)
      PetscSizeT                                         :: sizeofchar
      
      PetscCall(PetscDataTypeGetSize(PETSC_CHAR,sizeofchar,ierr))
      sizeofMEF90DefMechGlobalOptions = size(transfer(DefMechGlobalOptions,dummychar))*sizeofchar
      sizeofMEF90DefMechCellSetOptions = size(transfer(DefMechCellSetOptions,dummychar))*sizeofchar
      sizeofMEF90DefMechFaceSetOptions = size(transfer(DefMechFaceSetOptions,dummychar))*sizeofchar
      sizeofMEF90DefMechVertexSetOptions = size(transfer(DefMechVertexSetOptions,dummychar))*sizeofchar

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
      
      MEF90DefMech_damageTypeList(1)  = 'AT1'
      MEF90DefMech_damageTypeList(2)  = 'AT1exp'
      MEF90DefMech_damageTypeList(3)  = 'AT2'
      MEF90DefMech_damageTypeList(4)  = 'LinSoft'
      MEF90DefMech_damageTypeList(5)  = 'KKL'
      MEF90DefMech_damageTypeList(6)  = 'AT1Elastic'
      MEF90DefMech_damageTypeList(7)  = 'AT1expElastic'
      MEF90DefMech_damageTypeList(8)  = 'AT2Elastic'
      MEF90DefMech_damageTypeList(9)  = 'LinSoftElastic'
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
      MEF90DefMech_plasticityTypeList(9)  = 'Green'
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
   End Subroutine MEF90DefMechCtxInitialize_Private
   
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCtxCreate"
!!!
!!!  
!!!  MEF90DefMechCtxCreate:
!!!  
!!!  (c) 2012-20 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90DefMechCtxCreate(DefMechCtx,dm,MEF90Ctx,ierr)
      Type(MEF90DefMechCtx_Type),Intent(OUT)                   :: DefMechCtx
      Type(tDM),target,Intent(IN)                              :: dm
      Type(MEF90Ctx_Type),target,Intent(IN)                    :: MEF90Ctx
      PetscErrorCode,Intent(INOUT)                             :: ierr
   
      Type(MEF90CtxGlobalOptions_Type),pointer                 :: MEF90CtxGlobalOptions
      Type(tIS)                                                :: setIS
      PetscInt                                                 :: set,numSet
      Character(len=MEF90MXSTRLEN)                             :: filename,IOBuffer
      Character(len=MEF90MXSTRLEN)                             :: vecName
      Type(tDM),DImension(:),Pointer                           :: dmList
      Type(tPetscSF)                                           :: dummySF

      PetscCall(MEF90DefMechCtxInitialize_Private(ierr))
      DefMechCtx%MEF90Ctx => MEF90Ctx
      PetscCall(PetscBagCreate(MEF90Ctx%comm,sizeofMEF90DefMechGlobalOptions,DefMechCtx%GlobalOptionsBag,ierr))      
      PetscCall(DMGetLabelIdIS(dm,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr)) 
      PetscCall(ISGetLocalSize(setIS,numSet,ierr))
      Allocate(DefMechCtx%CellSetOptionsBag(numSet),stat=ierr)
      Do set = 1, numSet
         PetscCall(PetscBagCreate(MEF90Ctx%comm,sizeofMEF90DefMechCellSetOptions,DefMechCtx%CellSetOptionsBag(set),ierr))
      End Do
      PetscCall(ISDestroy(setIS,ierr))

      PetscCall(DMGetLabelIdIS(dm,MEF90FaceSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr)) 
      PetscCall(ISGetLocalSize(setIS,numSet,ierr))
      Allocate(DefMechCtx%FaceSetOptionsBag(numSet),stat=ierr)
      Do set = 1, numSet
         PetscCall(PetscBagCreate(MEF90Ctx%comm,sizeofMEF90DefMechFaceSetOptions,DefMechCtx%FaceSetOptionsBag(set),ierr))
      End Do
      PetscCall(ISDestroy(setIS,ierr))

      PetscCall(DMGetLabelIdIS(dm,MEF90VertexSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr)) 
      PetscCall(ISGetLocalSize(setIS,numSet,ierr))
      Allocate(DefMechCtx%VertexSetOptionsBag(numSet),stat=ierr)
      Do set = 1, numSet
         PetscCall(PetscBagCreate(MEF90Ctx%comm,sizeofMEF90DefMechVertexSetOptions,DefMechCtx%VertexSetOptionsBag(set),ierr))
      End Do
      PetscCall(ISDestroy(setIS,ierr))

      !!! 
      !!! Create energy viewers
      !!!      
      filename = trim(MEF90FilePrefix(MEF90Ctx%resultFile))//'.ener'
      PetscCall(PetscViewerASCIIOpen(MEF90Ctx%comm,filename,DefMechCtx%globalEnergyViewer,ierr))
      PetscCall(PetscViewerASCIIPrintf(DefMechCtx%globalEnergyViewer,"# step     load            elastic energy  work            cohesive energy surface energy  total energy   plastic dissipation \n",ierr))
      PetscCall(PetscViewerFlush(DefMechCtx%globalEnergyViewer,ierr))
      
      PetscCall(DMGetLabelIdIS(dm,MEF90FaceSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr)) 
      PetscCall(ISGetLocalSize(setIS,numSet,ierr))
      Allocate(DefMechCtx%setEnergyViewer(numSet),stat=ierr)
      Do set = 1, numSet
         Write(filename,101) trim(MEF90FilePrefix(MEF90Ctx%resultFile)),set
         PetscCall(PetscViewerASCIIOpen(MEF90Ctx%comm,filename,DefMechCtx%setEnergyViewer(set),ierr))
         Write(IOBuffer,102) set
         PetscCall(PetscViewerASCIIPrintf(DefMechCtx%setEnergyViewer(set),IOBuffer,ierr))
         PetscCall(PetscViewerASCIIPrintf(DefMechCtx%setEnergyViewer(set),"# step     load            elastic energy  work            cohesive energy surface energy  total energy   plastic dissipation\n",ierr))
         PetscCall(PetscViewerFlush(DefMechCtx%setEnergyViewer(set),ierr))
      End Do
      PetscCall(ISDestroy(setIS,ierr))
101 Format(A,'-',I4.4,'.enerblk')
102 Format("# cell set ",I4,"\n")
      DefMechCtx%analysisTime = 0.0_Kr
      DefMechCtx%timeStep = 0.0_Kr
   
      !!! Create Vecs and SF   
      PetscCall(DMGetDimension(dm,DefMechCtx%dim,ierr))

      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr))

      vecName = "Displacement"
      Allocate(DefMechCtx%displacementLocal,stat=ierr)
      PetscCall(MEF90CreateLocalVector(dm,MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,DefMechCtx%dim,vecName,DefMechCtx%displacementLocal,ierr)) 
      Allocate(DefMechCtx%displacementPreviousStepLocal,stat=ierr)
      PetscCall(VecDuplicate(DefMechCtx%displacementLocal,DefMechCtx%displacementPreviousStepLocal,ierr))
      PetscCall(PetscObjectSetName(DefMechCtx%displacementPreviousStepLocal,"DisplacementPreviousStep",ierr))
      Allocate(DefMechCtx%displacementLowerBoundLocal,stat=ierr)
      PetscCall(VecDuplicate(DefMechCtx%displacementLocal,DefMechCtx%displacementLowerBoundLocal,ierr))
      PetscCall(PetscObjectSetName(DefMechCtx%displacementLowerBoundLocal,"DisplacementLowerBound",ierr))
      Allocate(DefMechCtx%displacementUpperBoundLocal,stat=ierr)
      PetscCall(VecDuplicate(DefMechCtx%displacementLocal,DefMechCtx%displacementUpperBoundLocal,ierr))
      PetscCall(PetscObjectSetName(DefMechCtx%displacementUpperBoundLocal,"DisplacementUpperBound",ierr))

      vecName = "Damage"
      Allocate(DefMechCtx%damageLocal,stat=ierr)
      PetscCall(MEF90CreateLocalVector(dm,MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,1_Ki,vecName,DefMechCtx%damageLocal,ierr)) 
      Allocate(DefMechCtx%damagePreviousStepLocal,stat=ierr)
      PetscCall(VecDuplicate(DefMechCtx%damageLocal,DefMechCtx%damagePreviousStepLocal,ierr))
      PetscCall(PetscObjectSetName(DefMechCtx%damagePreviousStepLocal,"damagePreviousStep",ierr))

      Allocate(DefMechCtx%TemperatureLocal,stat=ierr)
      PetscCall(VecDuplicate(DefMechCtx%damageLocal,DefMechCtx%TemperatureLocal,ierr))
      PetscCall(PetscObjectSetName(DefMechCtx%TemperatureLocal,"Temperature",ierr))

      vecName = "cohesiveDisplacement"
      Allocate(DefMechCtx%cohesiveDisplacement,stat=ierr)
      PetscCall(MEF90CreateLocalVector(dm,MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,DefMechCtx%dim,vecName,DefMechCtx%cohesiveDisplacement,ierr)) 
      vecName = "bodyForce"
      Allocate(DefMechCtx%bodyForce,stat=ierr)
      PetscCall(MEF90CreateCellVector(dm,DefMechCtx%dim,vecName,DefMechCtx%bodyForce,ierr))
      vecName = "boundaryForce"
      Allocate(DefMechCtx%boundaryForce,stat=ierr)
      PetscCall(MEF90CreateBoundaryCellVector(dm,DefMechCtx%dim,vecName,DefMechCtx%boundaryForce,ierr))
      vecName = "pressureForce"
      Allocate(DefMechCtx%pressureForce,stat=ierr)
      PetscCall(MEF90CreateBoundaryCellVector(dm,1_Ki,vecName,DefMechCtx%pressureForce,ierr))

      vecName = "plasticStrain"
      Allocate(DefMechCtx%plasticStrain,stat=ierr)
      PetscCall(MEF90CreateCellVector(dm,(DefMechCtx%dim*(DefMechCtx%dim+1_Ki))/2_Ki,vecName,DefMechCtx%plasticStrain,ierr))
      Allocate(DefMechCtx%cumulatedPlasticDissipation,stat=ierr)
      PetscCall(VecDuplicate(DefMechCtx%plasticStrain,DefMechCtx%cumulatedPlasticDissipation,ierr))
      PetscCall(PetscObjectSetName(DefMechCtx%cumulatedPlasticDissipation,"cumulatedPlasticDissipation",ierr))
      Allocate(DefMechCtx%stress,stat=ierr)
      PetscCall(VecDuplicate(DefMechCtx%plasticStrain,DefMechCtx%stress,ierr))
      PetscCall(PetscObjectSetName(DefMechCtx%stress,"Stress",ierr))

      !!! Create megaDM
      Allocate(dmList(7))
      PetscCall(VecGetDM(DefMechCtx%displacementLocal,dmList(1),ierr))
      PetscCall(VecGetDM(DefMechCtx%damageLocal,dmList(2),ierr))
      PetscCall(VecGetDM(DefMechCtx%temperatureLocal,dmList(3),ierr))
      PetscCall(VecGetDM(DefMechCtx%bodyForce,dmList(4),ierr))
      PetscCall(VecGetDM(DefMechCtx%boundaryForce,dmList(5),ierr))
      PetscCall(VecGetDM(DefMechCtx%pressureForce,dmList(6),ierr))
      PetscCall(VecGetDM(DefMechCtx%plasticStrain,dmList(7),ierr))
      PetscCall(DMCreateSuperDM(dmList,7_kI,PETSC_NULL_IS,DefMechCtx%megaDM,ierr))
      DeAllocate(dmList)

      !!! Create the IO SF for all fields
      PetscCall(MEF90IOSFCreate(MEF90Ctx,DefMechCtx%displacementLocal,DefMechCtx%displacementToIOSF,DefMechCtx%IOTodisplacementSF,ierr))
      PetscCall(MEF90IOSFCreate(MEF90Ctx,DefMechCtx%damageLocal,DefMechCtx%damageToIOSF,DefMechCtx%IOTodamageSF,ierr))

      PetscCall(MEF90IOSFCreate(MEF90Ctx,DefMechCtx%cohesiveDisplacement,DefMechCtx%cohesiveDisplacementToIOSF,DefMechCtx%IOToCohesiveDisplacementSF,ierr))

      PetscCall(MEF90IOSFCreate(MEF90Ctx,DefMechCtx%temperatureLocal,DefMechCtx%temperatureToIOSF,DefMechCtx%IOTotemperatureSF,ierr))
      PetscCall(MEF90IOSFCreate(MEF90Ctx,DefMechCtx%bodyForce,DefMechCtx%bodyForceToIOSF,DefMechCtx%IOTobodyForceSF,ierr))
      PetscCall(MEF90IOSFCreate(MEF90Ctx,DefMechCtx%boundaryForce,DefMechCtx%boundaryForceToIOSF,DefMechCtx%IOToboundaryForceSF,ierr))
      PetscCall(MEF90IOSFCreate(MEF90Ctx,DefMechCtx%pressureForce,DefMechCtx%pressureForceToIOSF,DefMechCtx%IOTopressureForceSF,ierr))
      PetscCall(MEF90IOSFCreate(MEF90Ctx,DefMechCtx%stress,DefMechCtx%stressToIOSF,DefMechCtx%IOToStressSF,ierr))
      PetscCall(MEF90IOSFCreate(MEF90Ctx,DefMechCtx%plasticStrain,DefMechCtx%plasticStrainToIOSF,DefMechCtx%IOToplasticStrainSF,ierr))
      PetscCall(MEF90IOSFCreate(MEF90Ctx,DefMechCtx%cumulatedPlasticDissipation,DefMechCtx%cumulatedPlasticDissToIOSF,DefMechCtx%IOToCumulatedPlasticDissSF,ierr))

      !!! Create the SF to exchange boundary values of the displacement and damage. 
      PetscCall(MEF90ConstraintSFCreate(DefMechCtx%MEF90Ctx,DefMechCtx%displacementLocal,DefMechCtx%displacementLocal,DefMechCtx%displacementConstraintsSF,dummySF,ierr))
      PetscCall(PetscSFDestroy(dummySF,ierr))
      PetscCall(MEF90ConstraintSFCreate(DefMechCtx%MEF90Ctx,DefMechCtx%damageLocal,DefMechCtx%damageLocal,DefMechCtx%damageConstraintsSF,dummySF,ierr))
      PetscCall(PetscSFDestroy(dummySF,ierr))
   End Subroutine MEF90DefMechCtxCreate

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCtxDestroy"
!!!
!!!  
!!!  MEF90DefMechCtxDestroy: destroys a MEF90DefMechCtx_Type
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90DefMechCtxDestroy(DefMechCtx,ierr)
      Type(MEF90DefMechCtx_Type),Intent(INOUT)        :: DefMechCtx
      PetscErrorCode,Intent(INOUT)                    :: ierr
      
      PetscInt                                        :: set
   
      PetscCall(PetscBagDestroy(DefMechCtx%GlobalOptionsBag,ierr))
      Do set = 1, size(DefMechCtx%CellSetOptionsBag)
         PetscCall(PetscBagDestroy(DefMechCtx%CellSetOptionsBag(set),ierr))
      End Do
      DeAllocate(DefMechCtx%CellSetOptionsBag)
      Do set = 1, size(DefMechCtx%FaceSetOptionsBag)
         PetscCall(PetscBagDestroy(DefMechCtx%FaceSetOptionsBag(set),ierr))
      End Do
      DeAllocate(DefMechCtx%FaceSetOptionsBag)
      Do set = 1, size(DefMechCtx%VertexSetOptionsBag)
         PetscCall(PetscBagDestroy(DefMechCtx%VertexSetOptionsBag(set),ierr))
      End Do
      DeAllocate(DefMechCtx%VertexSetOptionsBag)
      
      !!! 
      !!! Close energy viewers
      !!!      
      PetscCall(PetscViewerDestroy(DefMechCtx%globalEnergyViewer,ierr))
      Do set = 1, size(DefMechCtx%setEnergyViewer)
         PetscCall(PetscViewerDestroy(DefMechCtx%setEnergyViewer(set),ierr))
      End Do
      DeAllocate(DefMechCtx%setEnergyViewer)

      !!! Destroy Vecs and SF and deAllocate them      
      If (Associated(DefMechCtx%displacementLocal)) Then
         PetscCall(VecDestroy(DefMechCtx%displacementLocal,ierr))
         DeAllocate(DefMechCtx%displacementLocal)
         Nullify(DefMechCtx%displacementLocal) 
      End If
      If (Associated(DefMechCtx%displacementPreviousStepLocal)) Then
         PetscCall(VecDestroy(DefMechCtx%displacementPreviousStepLocal,ierr))
         DeAllocate(DefMechCtx%displacementPreviousStepLocal)
         Nullify(DefMechCtx%displacementPreviousStepLocal) 
      End If
      If (Associated(DefMechCtx%displacementLowerBoundLocal)) Then
         PetscCall(VecDestroy(DefMechCtx%displacementLowerBoundLocal,ierr))
         DeAllocate(DefMechCtx%displacementLowerBoundLocal)
         Nullify(DefMechCtx%displacementLowerBoundLocal) 
      End If
      If (Associated(DefMechCtx%displacementUpperBoundLocal)) Then
         PetscCall(VecDestroy(DefMechCtx%displacementUpperBoundLocal,ierr))
         DeAllocate(DefMechCtx%displacementUpperBoundLocal)
         Nullify(DefMechCtx%displacementUpperBoundLocal) 
      End If

      If (Associated(DefMechCtx%damageLocal)) Then
         PetscCall(VecDestroy(DefMechCtx%damageLocal,ierr))
         DeAllocate(DefMechCtx%damageLocal)
         Nullify(DefMechCtx%damageLocal) 
      End If
      If (Associated(DefMechCtx%damagePreviousStepLocal)) Then
         PetscCall(VecDestroy(DefMechCtx%damagePreviousStepLocal,ierr))
         DeAllocate(DefMechCtx%damagePreviousStepLocal)
         Nullify(DefMechCtx%damagePreviousStepLocal) 
      End If

      If (Associated(DefMechCtx%temperatureLocal)) Then
         PetscCall(VecDestroy(DefMechCtx%temperatureLocal,ierr))
         DeAllocate(DefMechCtx%temperatureLocal)
         Nullify(DefMechCtx%temperatureLocal) 
      End If

      If (Associated(DefMechCtx%cohesiveDisplacement)) Then
         PetscCall(VecDestroy(DefMechCtx%cohesiveDisplacement,ierr))
         DeAllocate(DefMechCtx%cohesiveDisplacement)
         Nullify(DefMechCtx%cohesiveDisplacement) 
      End If

      If (Associated(DefMechCtx%bodyForce)) Then
         PetscCall(VecDestroy(DefMechCtx%bodyForce,ierr))
         DeAllocate(DefMechCtx%bodyForce)
         Nullify(DefMechCtx%bodyForce) 
      End If
      If (Associated(DefMechCtx%boundaryForce)) Then
         PetscCall(VecDestroy(DefMechCtx%boundaryForce,ierr))
         DeAllocate(DefMechCtx%boundaryForce)
         Nullify(DefMechCtx%boundaryForce) 
      End If
      If (Associated(DefMechCtx%pressureForce)) Then
         PetscCall(VecDestroy(DefMechCtx%pressureForce,ierr))
         DeAllocate(DefMechCtx%pressureForce)
         Nullify(DefMechCtx%pressureForce) 
      End If

      If (Associated(DefMechCtx%plasticStrain)) Then
         PetscCall(VecDestroy(DefMechCtx%plasticStrain,ierr))
         DeAllocate(DefMechCtx%plasticStrain)
         Nullify(DefMechCtx%plasticStrain) 
      End If
      If (Associated(DefMechCtx%cumulatedPlasticDissipation)) Then
         PetscCall(VecDestroy(DefMechCtx%cumulatedPlasticDissipation,ierr))
         DeAllocate(DefMechCtx%cumulatedPlasticDissipation)
         Nullify(DefMechCtx%cumulatedPlasticDissipation) 
      End If
      If (Associated(DefMechCtx%Stress)) Then
         PetscCall(VecDestroy(DefMechCtx%Stress,ierr))
         DeAllocate(DefMechCtx%Stress)
         Nullify(DefMechCtx%Stress) 
      End If

      !!! Destroy all PetscSF
      PetscCall(PetscSFDestroy(DefMechCtx%displacementToIOSF,ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%IOTodisplacementSF,ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%cohesiveDisplacementToIOSF,ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%IOToCohesiveDisplacementSF,ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%damageToIOSF,ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%IOTodamageSF,ierr))

      PetscCall(PetscSFDestroy(DefMechCtx%temperatureToIOSF,ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%IOTotemperatureSF,ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%bodyForceToIOSF,ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%IOTobodyForceSF,ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%boundaryForceToIOSF,ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%IOToboundaryForceSF,ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%pressureForceToIOSF,ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%IOTopressureForceSF,ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%stressToIOSF,ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%IOToStressSF,ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%plasticStrainToIOSF,ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%IOToplasticStrainSF,ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%cumulatedPlasticDissToIOSF,ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%IOToCumulatedPlasticDissSF,ierr))

      !!! Destroy the SF to exchange boundary values of the displacement and damage. 
      PetscCall(PetscSFDestroy(DefMechCtx%displacementConstraintsSF,ierr))
      PetscCall(PetscSFDestroy(DefMechCtx%damageConstraintsSF,ierr))

      !!! Destroy the megaDM
      PetscCall(DMDestroy(DefMechCtx%megaDM,ierr))
   End Subroutine MEF90DefMechCtxDestroy
   


#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90DefMechCtxGlobalOptions"
!!!
!!!  
!!!  PetscBagRegisterMEF90DefMechCtxGlobalOptions:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine PetscBagRegisterMEF90DefMechCtxGlobalOptions(bag,name,prefix,default,ierr)
      PetscBag                                                 :: bag
      Character(len=*),Intent(IN)                              :: prefix,name
      Type(MEF90DefMechGlobalOptions_Type),Intent(IN)          :: default
      PetscErrorCode,Intent(INOUT)                             :: ierr

      Type(MEF90DefMechGlobalOptions_Type),pointer             :: DefMechGlobalOptions

      PetscCall(PetscBagGetDataMEF90DefMechCtxGlobalOptions(bag,DefMechGlobalOptions,ierr))
      PetscCall(PetscBagSetName(bag,trim(name),"DefMechGlobalOptions MEF90 Defect Mechanics global options",ierr))
      PetscCall(PetscBagSetOptionsPrefix(bag,trim(prefix),ierr))

      PetscCall(PetscBagRegisterEnum(bag,DefMechGlobalOptions%timeSteppingType,MEF90DefMech_TimeSteppingTypeList,default%timeSteppingType,'DefMech_TimeStepping_Type','Type of defect mechanics Time steping',ierr))
      PetscCall(PetscBagRegisterEnum(bag,DefMechGlobalOptions%solverType,MEF90DefMech_SolverTypeList,default%solverType,'DefMech_solver_Type','Type of defect mechanics solver',ierr))
      PetscCall(PetscBagRegisterEnum(bag,DefMechGlobalOptions%damageSolverType,MEF90DefMech_DamageSolverTypeList,default%damageSolverType,'DefMech_damageSolver_Type','Type of defect mechanics damage solver',ierr))

      PetscCall(PetscBagRegisterBool(bag,DefMechGlobalOptions%temperatureExport,default%temperatureExport,'temperature_export','Export temperature in result file',ierr))
      PetscCall(PetscBagRegisterBool(bag,DefMechGlobalOptions%displacementExport,default%displacementExport,'displacement_export','Export displacement in result file',ierr))
      PetscCall(PetscBagRegisterBool(bag,DefMechGlobalOptions%damageExport,default%damageExport,'damage_export','Export damage in result file',ierr))
      PetscCall(PetscBagRegisterBool(bag,DefMechGlobalOptions%stressExport,default%stressExport,'stress_export','Export stress in result file',ierr))
      PetscCall(PetscBagRegisterBool(bag,DefMechGlobalOptions%plasticStrainExport,default%plasticStrainExport,'plasticstrain_export','Export plastic strain in result file',ierr))
      PetscCall(PetscBagRegisterBool(bag,DefMechGlobalOptions%cumulatedPlasticDissipationExport,default%cumulatedPlasticDissipationExport,'cumulatedplasticdissipation_export','Export cumulated plastic dissipation in result file',ierr))

      PetscCall(PetscBagRegisterEnum(bag,DefMechGlobalOptions%boundaryDisplacementScaling,MEF90ScalingList,default%boundaryDisplacementScaling,'boundaryDisplacement_scaling','Boundary displacement scaling',ierr))
      PetscCall(PetscBagRegisterEnum(bag,DefMechGlobalOptions%displacementLowerBoundScaling,MEF90ScalingList,default%displacementLowerBoundScaling,'displacementlowerbound_scaling','Displacement lower bound scaling',ierr))
      PetscCall(PetscBagRegisterEnum(bag,DefMechGlobalOptions%displacementUpperBoundScaling,MEF90ScalingList,default%displacementUpperBoundScaling,'displacementupperbound_scaling','Displacement upper bound scaling',ierr))
      PetscCall(PetscBagRegisterEnum(bag,DefMechGlobalOptions%cohesiveDisplacementScaling,MEF90ScalingList,default%cohesiveDisplacementScaling,'cohesiveDisplacement_scaling','Cohesive displacement scaling',ierr))
      PetscCall(PetscBagRegisterEnum(bag,DefMechGlobalOptions%boundaryDamageScaling,MEF90ScalingList,default%boundaryDamageScaling,'boundaryDamage_scaling','Boundary damage scaling',ierr))
      PetscCall(PetscBagRegisterEnum(bag,DefMechGlobalOptions%bodyForceScaling,MEF90ScalingList,default%bodyForceScaling,'bodyforce_scaling','Body force scaling',ierr))
      PetscCall(PetscBagRegisterEnum(bag,DefMechGlobalOptions%boundaryForceScaling,MEF90ScalingList,default%boundaryForceScaling,'boundaryforce_scaling','Boundary force scaling',ierr))
      PetscCall(PetscBagRegisterEnum(bag,DefMechGlobalOptions%pressureForceScaling,MEF90ScalingList,default%pressureforceScaling,'pressureForce_scaling','Pressure force scaling',ierr))
      PetscCall(PetscBagRegisterEnum(bag,DefMechGlobalOptions%CrackPressureScaling,MEF90ScalingList,default%CrackPressureScaling,'crackPressure_scaling','Crack Pressure scaling',ierr))

      PetscCall(PetscBagRegisterReal(bag,DefMechGlobalOptions%damageATol,default%damageATol,'defmech_damage_atol','Absolute tolerance on damage error',ierr))
      PetscCall(PetscBagRegisterInt (bag,DefMechGlobalOptions%damageMaxIt,default%damageMaxIt,'defmech_damage_maxit','Maximum number of alternate minimizations for damage',ierr))
      PetscCall(PetscBagRegisterInt (bag,DefMechGlobalOptions%PCLag,default%PCLag,'defmech_pclag','Interval at which the PC is recomputed during alternate minimization',ierr))
      PetscCall(PetscBagRegisterReal(bag,DefMechGlobalOptions%SOROmega,default%SOROmega,'defmech_SOR_Omega','Alterate Minimization over relaxation factor (>0 for limited, <0 for projected) ',ierr))
      PetscCall(PetscBagRegisterReal(bag,DefMechGlobalOptions%irrevthres,default%irrevthres,'defmech_irrevThres','Threshold above which irreversibility is enforced (0 for monotonicity, .99 for equality)',ierr))

      PetscCall(PetscBagRegisterEnum(bag,DefMechGlobalOptions%BTType,MEF90DefMech_BTTypeList,default%BTType,'BT_Type','Backtracking type',ierr))
      PetscCall(PetscBagRegisterInt (bag,DefMechGlobalOptions%BTInterval,default%BTInterval,'BT_Interval','Interval at which Backtracking is run in inner loop (0 for outer loop)',ierr))
      PetscCall(PetscBagRegisterInt (bag,DefMechGlobalOptions%BTScope,default%BTScope,'BT_Scope','Backtracking scope (0 for unlimited)',ierr))
      PetscCall(PetscBagRegisterReal(bag,DefMechGlobalOptions%BTTol,default%BTTol,'BT_Tol','Backtracking relative tolerance',ierr))

      PetscCall(PetscBagRegisterReal(bag,DefMechGlobalOptions%plasticStrainATol,default%plasticStrainATol,'defmech_plasticstrain_atol','Absolute tolerance on plastic error',ierr))
      
      PetscCall(PetscBagRegisterReal(bag,DefMechGlobalOptions%InjectedVolumeATol,default%InjectedVolumeATol,'defmech_InjectedVolume_atol','Absolute tolerance on injected volume error',ierr))
      PetscCall(PetscBagRegisterReal(bag,DefMechGlobalOptions%dampingCoefficientDisplacement,default%dampingCoefficientDisplacement,'defmech_dampingCoefficient_displacement','Damping coefficient on displacement field (0 for minimization, 1 for semi-implicit gradient flow)',ierr))
      PetscCall(PetscBagRegisterReal(bag,DefMechGlobalOptions%dampingCoefficientDamage,default%dampingCoefficientDamage,'defmech_dampingCoefficient_damage','Damping coefficient on damage field (0 for minimization, 1 for semi-implicit gradient flow)',ierr))
   End Subroutine PetscBagRegisterMEF90DefMechCtxGlobalOptions

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90DefMechCtxCellSetOptions"
!!!
!!!  
!!!  PetscBagRegisterMEF90DefMechCtxCellSetOptions:
!!!  
!!!  (c) 2012-18 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine PetscBagRegisterMEF90DefMechCtxCellSetOptions(bag,name,prefix,default,ierr)
      PetscBag                                           :: bag
      Character(len=*),Intent(IN)                        :: prefix,name
      Type(MEF90DefMechCellSetOptions_Type),Intent(IN)   :: default
      PetscErrorCode,Intent(INOUT)                       :: ierr

      Type(MEF90DefMechCellSetOptions_Type),pointer      :: DefMechCellSetOptions
      
      PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(bag,DefMechCellSetOptions,ierr))
      PetscCall(PetscBagSetName(bag,trim(name),"DefMechCellSetOptions MEF90 Defect Mechanics Cell Set options",ierr))
      PetscCall(PetscBagSetOptionsPrefix(bag,trim(prefix),ierr))
      
      DefMechCellSetOptions%bodyForce              = default%bodyForce
      DefMechCellSetOptions%boundaryDisplacement   = default%boundaryDisplacement
      DefMechCellSetOptions%displacementLowerBound = default%displacementLowerBound
      DefMechCellSetOptions%displacementUpperBound = default%displacementUpperBound
      DefMechCellSetOptions%Has_displacementBC     = default%Has_displacementBC

      PetscCall(PetscBagRegisterRealArray(bag,DefMechCellSetOptions%bodyForce,3_Ki,'bodyForce','[N.m^(-3) / N.m^(-2)] (f): body force',ierr))
      PetscCall(PetscBagRegisterReal(bag,DefMechCellSetOptions%CrackPressure,default%CrackPressure,'CrackPressure','without unit: internal crack pressure',ierr))
      PetscCall(PetscBagRegisterReal(bag,DefMechCellSetOptions%DamageATLinSoftk,default%DamageATLinSoftk,'damage_LinSoft_k','[unit-less] (k): k parameter in the Linear Softening damage model',ierr))
      PetscCall(PetscBagRegisterReal(bag,DefMechCellSetOptions%DamageAT1expb,default%DamageAT1expb,'damage_AT1exp_b','[unit-less] (b): b parameter in tha AT1 model with exponential stiffness interpolation',ierr))
      PetscCall(PetscBagRegisterEnum(bag,DefMechCellSetOptions%damageType,MEF90DefMech_damageTypeList,default%damageType,'damage_type','Type of damage law',ierr))
      PetscCall(PetscBagRegisterEnum(bag,DefMechCellSetOptions%plasticityType,MEF90DefMech_plasticityTypeList,default%plasticityType,'plasticity_type','Type of plasticity law',ierr))
      PetscCall(PetscBagRegisterEnum(bag,DefMechCellSetOptions%unilateralContactType,MEF90DefMech_unilateralContactTypeList,default%unilateralContactType,'unilateralContact_type','Type of handling of unilateral contact',ierr))
      PetscCall(PetscBagRegisterReal(bag,DefMechCellSetOptions%unilateralContactHydrostaticDeviatoricGamma,default%unilateralContactHydrostaticDeviatoricGamma,'unilateralContact_hydrostaticDeviatoric_gamma','[unit-less] (gamma): Hydrostatic Deviatoric regularization parameter',ierr))
      PetscCall(PetscBagRegisterBool(bag,DefMechCellSetOptions%unilateralContactHybrid,default%unilateralContactHybrid,'unilateralContact_hybrid','Use hybrid unilateral contact formulation (Y/N)',ierr))
      PetscCall(PetscBagRegisterEnum(bag,DefMechCellSetOptions%drivingForceType,MEF90DefMech_drivingForceTypeList,default%drivingForceType,'drivingForce_type','Type of nucleation driving force',ierr))
      PetscCall(PetscBagRegisterRealArray(bag,DefMechCellSetOptions%cohesiveDisplacement,3_Ki,'cohesiveDisplacement','[m] (U): Cohesive displacement value',ierr))
      PetscCall(PetscBagRegisterBoolArray(bag,DefMechCellSetOptions%Has_displacementBC,3_Ki,'DisplacementBC','Displacement has Dirichlet boundary Condition (Y/N)',ierr))
      PetscCall(PetscBagRegisterRealArray(bag,DefMechCellSetOptions%boundaryDisplacement,3_Ki,'boundaryDisplacement','[m] (U): Displacement boundary value',ierr))
      PetscCall(PetscBagRegisterRealArray(bag,DefMechCellSetOptions%displacementLowerBound,3_Ki,'displacementLowerBound','[m] (U): Displacement lower bound',ierr))
      PetscCall(PetscBagRegisterRealArray(bag,DefMechCellSetOptions%displacementUpperBound,3_Ki,'displacementUpperBound','[m] (U): Displacement upper bound',ierr))
      PetscCall(PetscBagRegisterBool(bag,DefMechCellSetOptions%Has_DamageBC,default%Has_DamageBC,'DamageBC','Damage has Dirichlet boundary Condition (Y/N)',ierr))
      PetscCall(PetscBagRegisterBool(bag,DefMechCellSetOptions%CrackVolumeControlled,default%CrackVolumeControlled,'CrackVolumeControlled','Crack Pressure controlled by the crack volume in this block (Y/N)',ierr))
      PetscCall(PetscBagRegisterBool(bag,DefMechCellSetOptions%WorkControlled,default%WorkControlled,'WorkControlled','Force magnitude controlled by its work in this block (Y/N)',ierr))
      PetscCall(PetscBagRegisterReal(bag,DefMechCellSetOptions%boundaryDamage,default%boundaryDamage,'boundaryDamage','[unit-less] (alpha): Damage boundary value',ierr))
   End Subroutine PetscBagRegisterMEF90DefMechCtxCellSetOptions

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90DefMechCtxFaceSetOptions"
!!!
!!!  
!!!  PetscBagRegisterMEF90DefMechCtxFaceSetOptions:
!!!  
!!!  (c) 2012-18 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine PetscBagRegisterMEF90DefMechCtxFaceSetOptions(bag,name,prefix,default,ierr)
      PetscBag                                           :: bag
      Character(len=*),Intent(IN)                        :: prefix,name
      Type(MEF90DefMechFaceSetOptions_Type),Intent(IN)   :: default
      PetscErrorCode,Intent(INOUT)                       :: ierr

      Type(MEF90DefMechFaceSetOptions_Type),pointer      :: DefMechFaceSetOptions
      
      PetscCall(PetscBagGetDataMEF90DefMechCtxFaceSetOptions(bag,DefMechFaceSetOptions,ierr))
      PetscCall(PetscBagSetName(bag,trim(name),"DefMechFaceSetOptions MEF90 Defect Mechanics Face Set options",ierr))
      PetscCall(PetscBagSetOptionsPrefix(bag,trim(prefix),ierr))
      
      DefMechFaceSetOptions%boundaryForce          = default%boundaryForce
      DefMechFaceSetOptions%boundaryDisplacement   = default%boundaryDisplacement
      DefMechFaceSetOptions%displacementLowerBound = default%displacementLowerBound
      DefMechFaceSetOptions%displacementUpperBound = default%displacementUpperBound
      DefMechFaceSetOptions%Has_displacementBC     = default%Has_displacementBC

      PetscCall(PetscBagRegisterRealArray(bag,DefMechFaceSetOptions%boundaryForce,3_Ki,'boundaryForce','[N.m^(-2) / N.m^(-1)] (f): boundary force',ierr))
      PetscCall(PetscBagRegisterReal(bag,DefMechFaceSetOptions%pressureForce,default%pressureForce,'pressureForce','without unit: internal crack pressure',ierr))
      PetscCall(PetscBagRegisterBoolArray(bag,DefMechFaceSetOptions%Has_displacementBC,3_Ki,'DisplacementBC','Displacement has Dirichlet boundary Condition (Y/N)',ierr))
      PetscCall(PetscBagRegisterRealArray(bag,DefMechFaceSetOptions%boundaryDisplacement,3_Ki,'boundaryDisplacement','[m] (U): Displacement boundary value',ierr))
      PetscCall(PetscBagRegisterRealArray(bag,DefMechFaceSetOptions%displacementLowerBound,3_Ki,'displacementLowerBound','[m] (U): Displacement lower bound',ierr))
      PetscCall(PetscBagRegisterRealArray(bag,DefMechFaceSetOptions%displacementUpperBound,3_Ki,'displacementUpperBound','[m] (U): Displacement upper bound',ierr))
      PetscCall(PetscBagRegisterBool(bag,DefMechFaceSetOptions%Has_DamageBC,default%Has_DamageBC,'DamageBC','Damage has Dirichlet boundary Condition (Y/N)',ierr))
      PetscCall(PetscBagRegisterReal(bag,DefMechFaceSetOptions%boundaryDamage,default%boundaryDamage,'boundaryDamage','[unit-less] (alpha): Damage boundary value',ierr))
   End Subroutine PetscBagRegisterMEF90DefMechCtxFaceSetOptions

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90DefMechCtxVertexSetOptions"
!!!
!!!  
!!!  PetscBagRegisterMEF90DefMechCtxVertexSetOptions:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!
   Subroutine PetscBagRegisterMEF90DefMechCtxVertexSetOptions(bag,name,prefix,default,ierr)
      PetscBag                                              :: bag
      Character(len=*),Intent(IN)                           :: prefix,name
      Type(MEF90DefMechVertexSetOptions_Type),Intent(IN)    :: default
      PetscErrorCode,Intent(INOUT)                          :: ierr

      Type(MEF90DefMechVertexSetOptions_Type),pointer       :: DefMechVertexSetOptions
      PetscCall(PetscBagGetDataMEF90DefMechCtxVertexSetOptions(bag,DefMechVertexSetOptions,ierr))
      PetscCall(PetscBagSetName(bag,trim(name),"DefMechVertexSetOptions MEF90 Defect Mechanics Vertex Set options",ierr))
      PetscCall(PetscBagSetOptionsPrefix(bag,trim(prefix),ierr))

      DefMechVertexSetOptions%Has_displacementBC     = default%Has_displacementBC
      DefMechVertexSetOptions%boundaryDisplacement   = default%boundaryDisplacement
      DefMechVertexSetOptions%displacementLowerBound = default%displacementLowerBound
      DefMechVertexSetOptions%displacementUpperBound = default%displacementUpperBound
      PetscCall(PetscBagRegisterBoolArray(bag,DefMechVertexSetOptions%Has_displacementBC,3_Ki,'DisplacementBC','Displacement has Dirichlet boundary Condition (Y/N)',ierr))
      PetscCall(PetscBagRegisterRealArray(bag,DefMechVertexSetOptions%boundaryDisplacement,3_Ki,'boundaryDisplacement','[m] (U): Displacement boundary value',ierr))
      PetscCall(PetscBagRegisterRealArray(bag,DefMechVertexSetOptions%displacementLowerBound,3_Ki,'displacementLowerBound','[m] (U): Displacement lower bound',ierr))
      PetscCall(PetscBagRegisterRealArray(bag,DefMechVertexSetOptions%displacementUpperBound,3_Ki,'displacementUpperBound','[m] (U): Displacement upper bound',ierr))
      PetscCall(PetscBagRegisterBool(bag,DefMechVertexSetOptions%Has_DamageBC,default%Has_DamageBC,'DamageBC','Damage has Dirichlet boundary Condition (Y/N)',ierr))
      PetscCall(PetscBagRegisterReal(bag,DefMechVertexSetOptions%boundaryDamage,default%boundaryDamage,'boundaryDamage','[unit-less] (alpha): boundaryDamage',ierr))
   End Subroutine PetscBagRegisterMEF90DefMechCtxVertexSetOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCtxSetFromOptions"
!!!
!!!  
!!!  MEF90DefMechCtxSetFromOptions:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90DefMechCtxSetFromOptions(DefMechCtx,prefix,defaultGlobalOptions, &
                                              defaultCellSetOptions,    &
                                              defaultFaceSetOptions,    &
                                              defaultVertexSetOptions,ierr)
      Type(MEF90DefMechCtx_Type),Intent(INOUT)              :: DefMechCtx
      Character(len=*),Intent(IN)                           :: prefix
      Type(MEF90DefMechGlobalOptions_Type),Intent(IN)       :: defaultGlobalOptions
      Type(MEF90DefMechCellSetOptions_Type),Intent(IN)      :: defaultCellSetOptions
      Type(MEF90DefMechFaceSetOptions_Type),Intent(IN)      :: defaultFaceSetOptions
      Type(MEF90DefMechVertexSetOptions_Type),Intent(IN)    :: defaultVertexSetOptions
      PetscErrorCode,Intent(INOUT)                          :: ierr
   
      Type(MEF90CtxGlobalOptions_Type),pointer              :: MEF90CtxGlobalOptions
      Type(tIS)                                             :: setIS
      PetscInt,Dimension(:),Pointer                         :: setID
      PetscInt                                              :: set
      Character(len=MEF90MXSTRLEN)                          :: IOBuffer,setName,setprefix

      Type(MEF90DefMechCellSetOptions_Type),pointer         :: cellSetOptions    
      Type(MEF90DefMechFaceSetOptions_Type),pointer         :: faceSetOptions    

      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr))
      !!!
      !!! Registering Global Context
      !!!
      PetscCall(PetscBagRegisterMEF90DefMechCtxGlobalOptions(DefMechCtx%GlobalOptionsBag,"MEF90DefMech Global Ctx",prefix,defaultGlobalOptions,ierr))

      If (MEF90CtxGlobalOptions%verbose > 0) Then
         PetscCall(PetscBagView(DefMechCtx%GlobalOptionsBag,PETSC_VIEWER_STDOUT_WORLD,ierr))
      End If

      !!!
      !!! Registering Cell Set Context
      !!!
      PetscCall(DMGetLabelIdIS(DefMechCtx%megaDM,MEF90CellSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(DefMechCtx%MEF90Ctx%comm,setIS,ierr)) 
      PetscCall(ISGetIndicesF90(setIS,setID,ierr))
      Do set = 1, size(setID)
         Write(setName,"('Cell set ',I4)") setID(set)
         Write(setprefix,"('cs',I4.4,'_')") setID(set)
         PetscCall(PetscBagRegisterMEF90DefMechCtxCellSetOptions(DefMechCtx%CellSetOptionsBag(set),setName,setPrefix,defaultCellSetOptions,ierr))
         If (MEF90CtxGlobalOptions%verbose > 0) Then
            Write(IOBuffer,"('\nRegistering cell set ',I4,' prefix: ',A,'\n')") setID(set),trim(setprefix)
            PetscCall(PetscPrintf(DefMechCtx%MEF90Ctx%comm,IOBuffer,ierr))
            PetscCall(PetscBagView(DefMechCtx%CellSetOptionsBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr))
            PetscCall(PetscPrintf(DefMechCtx%MEF90Ctx%comm,"\n",ierr))
         End if
      End Do
      PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      PetscCall(ISDestroy(setIS,ierr))

      !!!
      !!! Registering Face Set Context
      !!!
      PetscCall(DMGetLabelIdIS(DefMechCtx%megaDM,MEF90FaceSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(DefMechCtx%MEF90Ctx%comm,setIS,ierr)) 
      PetscCall(ISGetIndicesF90(setIS,setID,ierr))
      Do set = 1, size(setID)
         Write(setName,"('Face set ',I4)") setID(set)
         Write(setprefix,"('fs',I4.4,'_')") setID(set)
         PetscCall(PetscBagRegisterMEF90DefMechCtxFaceSetOptions(DefMechCtx%FaceSetOptionsBag(set),setName,setPrefix,defaultFaceSetOptions,ierr))
         If (MEF90CtxGlobalOptions%verbose > 0) Then
            Write(IOBuffer,"('\nRegistering face set ',I4,' prefix: ',A,'\n')") setID(set),trim(setprefix)
            PetscCall(PetscPrintf(DefMechCtx%MEF90Ctx%comm,IOBuffer,ierr))
            PetscCall(PetscBagView(DefMechCtx%FaceSetOptionsBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr))
            PetscCall(PetscPrintf(DefMechCtx%MEF90Ctx%comm,"\n",ierr))
         End if
      End Do
      PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      PetscCall(ISDestroy(setIS,ierr))
      
      !!!
      !!! Registering Vertex Set Context
      !!!
      PetscCall(DMGetLabelIdIS(DefMechCtx%megaDM,MEF90VertexSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr)) 
      PetscCall(ISGetIndicesF90(setIS,setID,ierr))
      Do set = 1, size(setID)
         Write(setName,"('Vertex set ',I4)") setID(set)
         Write(setprefix,"('vs',I4.4,'_')") setID(set)
         PetscCall(PetscBagRegisterMEF90DefMechCtxVertexSetOptions(DefMechCtx%VertexSetOptionsBag(set),setName,setPrefix,defaultVertexSetOptions,ierr))
         If (MEF90CtxGlobalOptions%verbose > 0) Then
            Write(IOBuffer,"('\nRegistering vertex set ',I4,' prefix: ',A,'\n')") setID(set),trim(setprefix)
            PetscCall(PetscPrintf(DefMechCtx%MEF90Ctx%comm,IOBuffer,ierr))
            PetscCall(PetscBagView(DefMechCtx%VertexSetOptionsBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr))
            PetscCall(PetscPrintf(DefMechCtx%MEF90Ctx%comm,"\n",ierr))
         End if
      End Do
      PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      PetscCall(ISDestroy(setIS,ierr))
      
      !!! Identify if any of the displacement bounds are set
      DefMechCtx%hasDisplacementBounds = PETSC_FALSE

      PetscCall(DMGetLabelIdIS(DefMechCtx%megaDM,MEF90CellSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(DefMechCtx%MEF90Ctx%comm,setIS,ierr)) 
      PetscCall(ISGetIndicesF90(setIS,setID,ierr))
      Do set = 1, size(setID)
         PetscCall(PetscBagGetDataMEF90DefMechCtxCellSetOptions(DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr))
         DefMechCtx%hasDisplacementBounds = any(cellSetOptions%displacementLowerBound /= MEF90NINFINITY) .OR. DefMechCtx%hasDisplacementBounds
         DefMechCtx%hasDisplacementBounds = any(cellSetOptions%displacementUpperBound /= MEF90INFINITY)  .OR. DefMechCtx%hasDisplacementBounds
      End Do
      PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      PetscCall(ISDestroy(setIS,ierr))
      PetscCall(DMGetLabelIdIS(DefMechCtx%megaDM,MEF90FaceSetLabelName, SetIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(DefMechCtx%MEF90Ctx%comm,setIS,ierr)) 
      PetscCall(ISGetIndicesF90(setIS,setID,ierr))
      Do set = 1, size(setID)
         PetscCall(PetscBagGetDataMEF90DefMechCtxFaceSetOptions(DefMechCtx%FaceSetOptionsBag(set),faceSetOptions,ierr))
         DefMechCtx%hasDisplacementBounds = any(faceSetOptions%displacementLowerBound /= MEF90NINFINITY) .OR. DefMechCtx%hasDisplacementBounds
         DefMechCtx%hasDisplacementBounds = any(faceSetOptions%displacementUpperBound /= MEF90INFINITY)  .OR. DefMechCtx%hasDisplacementBounds
      End Do
      PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      PetscCall(ISDestroy(setIS,ierr))
   End Subroutine MEF90DefMechCtxSetFromOptions
End Module m_MEF90_DefMechCtx
