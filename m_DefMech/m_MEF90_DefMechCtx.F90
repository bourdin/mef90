#include "../MEF90/mef90.inc"
Module m_MEF90_DefMechCtx_Type
#include "finclude/petscdef.h"
   Use m_MEF90
   Implicit none
   
   Type MEF90DefMechCtx_Type
      PetscReal                              :: analysisTime,timeStep
      !!!  vertex based vec
      Type(Vec),pointer                      :: displacement,displacementPreviousStep
      Type(Vec),pointer                      :: damage,damagePreviousStep
      Type(Vec),pointer                      :: boundaryDisplacement
      Type(Vec),pointer                      :: displacementLowerBound,displacementUpperBound
      Type(Vec),pointer                      :: boundaryDamage
      Type(Vec),Pointer                      :: temperature
      
      !!! cell based vec      
      Type(Vec),pointer                      :: force
      Type(Vec),pointer                      :: pressureForce
      Type(Vec),pointer                      :: CrackPressure
      Type(Vec),Pointer                      :: plasticStrain
      Type(Vec),Pointer                      :: cumulatedDissipatedPlasticEnergy
      Type(Vec),Pointer                      :: stress
      
      PetscBag                               :: GlobalOptionsBag
      PetscBag,Dimension(:),Pointer          :: CellSetOptionsBag
      PetscBag,Dimension(:),Pointer          :: VertexSetOptionsBag
      PetscBag,Dimension(:),Pointer          :: MaterialPropertiesBag
      Type(MEF90Ctx_Type),pointer            :: MEF90Ctx
      Type(DM),pointer                       :: DM
      Type(DM)                               :: DMScal,DMVect           ! Remove all these
      Type(DM)                               :: cellDMScal,cellDMVect   ! after switching to 
      Type(DM)                               :: DMMatS,cellDMMatS       ! DMplex
      
      Type(VecScatter)                       :: DMScalScatter
      Type(VecScatter)                       :: DMVectScatter
      Type(VecScatter)                       :: DMMatSScatter
      Type(VecScatter)                       :: cellDMScalScatter
      Type(VecScatter)                       :: cellDMVectScatter
      Type(VecScatter)                       :: cellDMMatSScatter
      
      Type(SectionReal)                      :: DMSec
      Type(SectionReal)                      :: DMScalSec,DMVectSec           ! Remove all these
      Type(SectionReal)                      :: cellDMScalSec,cellDMVectSec   ! after switching to 
      Type(SectionReal)                      :: DMMatSSec,cellDMMatSSec       ! DMplex

      Type(PetscViewer)                      :: globalEnergyViewer
      Type(PetscViewer),Dimension(:),Pointer :: setEnergyViewer

      Type(SNES)                             :: SNESDisp,SNESDamage
      PetscBool                              :: hasDisplacementBounds
      PetscBool                              :: hasUnilateralContact
   End Type MEF90DefMechCtx_Type
   
   Type MEF90DefMechGlobalOptions_Type
      PetscInt                               :: timeSteppingType
      PetscInt                               :: solverType
      PetscBool                              :: addDisplacementNullSpace
      !!! Position of vertex-based vecs in exo file
      PetscInt                               :: displacementOffset
      PetscInt                               :: damageOffset
      PetscInt                               :: boundaryDisplacementOffset
      PetscInt                               :: boundaryDamageOffset
      PetscInt                               :: temperatureOffset
      !!! Position of cell-based vecs in exo files
      PetscInt                               :: forceOffset
      PetscInt                               :: pressureForceOffset
      PetscInt                               :: CrackPressureOffset
      PetscInt                               :: plasticStrainOffset
      PetscInt                               :: stressOffset
      PetscInt                               :: cumulatedPlasticDissipationOffset
      !!! scaling = time (step) scaling law currently CST, Linear, or File
      PetscInt                               :: boundaryDisplacementScaling
      PetscInt                               :: displacementLowerBoundScaling
      PetscInt                               :: displacementUpperBoundScaling
      PetscInt                               :: boundaryDamageScaling
      PetscInt                               :: forceScaling
      PetscInt                               :: pressureForceScaling
      PetscInt                               :: CrackPressureScaling
      PetscReal                              :: damageATol
      PetscInt                               :: maxit
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
   End Type MEF90DefMechGlobalOptions_Type

   Type MEF90DefMechCellSetOptions_Type
      PetscInt                               :: elemTypeShortIDDisplacement
      PetscInt                               :: elemTypeShortIDDamage
      PetscReal,Dimension(3)                 :: force
      PetscReal                              :: pressureForce
      PetscReal                              :: CrackPressure
      PetscEnum                              :: damageType
      PetscEnum                              :: plasticityType
      PetscEnum                              :: unilateralContactType
      PetscInt                               :: drivingForceType
      PetscBool,Dimension(3)                 :: Has_displacementBC
      PetscReal,Dimension(3)                 :: boundaryDisplacement
      PetscReal,Dimension(3)                 :: displacementLowerBound
      PetscReal,Dimension(3)                 :: displacementUpperBound
      PetscBool                              :: Has_damageBC
      PetscBool                              :: CrackVolumeControlled
      PetscBool                              :: WorkControlled
      PetscReal                              :: boundaryDamage
   End Type MEF90DefMechCellSetOptions_Type

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
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
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
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataMEF90DefMechCtxGlobalOptions
End Module m_MEF90DefMechGlobalOptions_Private

Module m_MEF90DefMechCellSetOptions_Private
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
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
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataMEF90DefMechCtxCellSetOptions
End Module m_MEF90DefMechCellSetOptions_Private

Module m_MEF90DefMechVertexSetOptions_Private
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
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
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataMEF90DefMechCtxVertexSetOptions
End Module m_MEF90DefMechVertexSetOptions_Private

Module m_MEF90_DefMechCtx
#include "finclude/petscdef.h"
   Use m_MEF90
   Use m_MEF90_DefMechCtx_Type
   Use m_MEF90DefMechGlobalOptions_Private
   Use m_MEF90DefMechCellSetOptions_Private
   Use m_MEF90DefMechVertexSetOptions_Private
   Implicit none

   PetscSizeT,protected   :: sizeofMEF90DefMechGlobalOptions
   PetscSizeT,protected   :: sizeofMEF90DefMechCellSetOptions
   PetscSizeT,protected   :: sizeofMEF90DefMechVertexSetOptions
   
   Enum,bind(c)
      enumerator  :: MEF90DefMech_SolverTypeAltMin = 0, &
                     MEF90DefMech_SolverTypeQuasiNewton1, &
                     MEF90DefMech_SolverTypeQuasiNewton2
   End Enum
   Character(len = MEF90_MXSTRLEN),Dimension(6),protected   :: MEF90DefMech_SolverTypeList

   Enum,bind(c)
      enumerator  :: MEF90DefMech_TimeSteppingTypeNULL = 0,     &
                     MEF90DefMech_TimeSteppingTypeQuasiStatic
   End Enum
   Character(len = MEF90_MXSTRLEN),Dimension(5),protected   :: MEF90DefMech_TimeSteppingTypeList
   
   Enum,bind(c)
      enumerator  :: MEF90DefMech_BTTypeNULL = 0,    &
                     MEF90DefMech_BTTypeBackward,    &
                     MEF90DefMech_BTTypeForward
   End Enum
   Character(len = MEF90_MXSTRLEN),Dimension(6),protected   :: MEF90DefMech_BTTypeList
   
   Enum,bind(c)
      enumerator  :: MEF90DefMech_damageTypeAT1 = 0,        &
                     MEF90DefMech_damageTypeAT2,            &
                     MEF90DefMech_damageTypeLinSoft,        &
                     MEF90DefMech_damageTypeKKL,            &
                     MEf90DefMech_damageTypeAT1Elastic,     &
                     MEf90DefMech_damageTypeAT2Elastic,     &
                     MEF90DefMech_damageTypeLinSoftElastic, &
                     MEF90DefMech_damageTypeKKLElastic

   End Enum
   Character(len = MEF90_MXSTRLEN),Dimension(11),protected   :: MEF90DefMech_damageTypeList
   
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
   Character(len = MEF90_MXSTRLEN),Dimension(13),protected   :: MEF90DefMech_plasticityTypeList
   
   Enum,bind(c)
      enumerator  :: MEF90DefMech_unilateralContactTypeNone = 0,                     &
                     MEF90DefMech_unilateralContactTypeHydrostaticDeviatoric,        &
                     MEF90DefMech_unilateralContactTypeBrittleDuctile,               &
                     MEF90DefMech_unilateralContactTypePositiveHydrostatic,          &
                     MEF90DefMech_unilateralContactTypeHybridHydrostaticDeviatoric,  &
                     MEF90DefMech_unilateralContactTypeDeviatoric,                   &
                     MEF90DefMech_unilateralContactTypePrincipalStrains,             &
                     MEF90DefMech_unilateralContactTypeMasonry
   End Enum
   Character(len = MEF90_MXSTRLEN),Dimension(11),protected   :: MEF90DefMech_unilateralContactTypeList

   Enum,bind(c)
      enumerator  :: MEF90DefMech_drivingForceTypeNone = 0,       &
                     MEF90DefMech_drivingForceTypeDruckerPrager,  &
                     MEF90DefMech_drivingForceTypeDruckerPrager2
   End Enum
   Character(len = MEF90_MXSTRLEN),Dimension(6),protected   :: MEF90DefMech_drivingForceTypeList
Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCtxInitialize_Private"
!!!
!!!  
!!!  MEF90DefMechCtxInitialize_Private:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechCtxInitialize_Private(ierr)
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(MEF90DefMechGlobalOptions_Type)               :: DefMechGlobalOptions
      Type(MEF90DefMechCellSetOptions_Type)              :: DefMechCellSetOptions
      Type(MEF90DefMechVertexSetOptions_Type)            :: DefMechVertexSetOptions
      character(len=1),pointer                           :: dummychar(:)
      PetscSizeT                                         :: sizeofchar
      
      Call PetscDataTypeGetSize(PETSC_CHAR,sizeofchar,ierr)
      sizeofMEF90DefMechGlobalOptions = size(transfer(DefMechGlobalOptions,dummychar))*sizeofchar
      sizeofMEF90DefMechCellSetOptions = size(transfer(DefMechCellSetOptions,dummychar))*sizeofchar
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
      
      MEF90DefMech_BTTypeList(1) = 'Null'
      MEF90DefMech_BTTypeList(2) = 'Backward'
      MEF90DefMech_BTTypeList(3) = 'Forward'
      MEF90DefMech_BTTypeList(4) = 'MEF90DefMech_BTType'
      MEF90DefMech_BTTypeList(5) = '_MEF90DefMech_BTType'
      MEF90DefMech_BTTypeList(6) = ''
      
      MEF90DefMech_damageTypeList(1)  = 'AT1'
      MEF90DefMech_damageTypeList(2)  = 'AT2'
      MEF90DefMech_damageTypeList(3)  = 'LinSoft'
      MEF90DefMech_damageTypeList(4)  = 'KKL'
      MEF90DefMech_damageTypeList(5)  = 'AT1Elastic'
      MEF90DefMech_damageTypeList(6)  = 'AT2Elastic'
      MEF90DefMech_damageTypeList(7)  = 'LinSoftElastic'
      MEF90DefMech_damageTypeList(8)  = 'KKLElastic'
      MEF90DefMech_damageTypeList(9)  = 'MEF90DefMech_damageType'
      MEF90DefMech_damageTypeList(10) = '_MEF90DefMech_damageType'
      MEF90DefMech_damageTypeList(11) = ''

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
      MEF90DefMech_unilateralContactTypeList(3) = 'BrittleDuctile'
      MEF90DefMech_unilateralContactTypeList(4) = 'PositiveHydrostatic'
      MEF90DefMech_unilateralContactTypeList(5) = 'HybridHydrostaticDeviatoric'
      MEF90DefMech_unilateralContactTypeList(6) = 'Deviatoric'
      MEF90DefMech_unilateralContactTypeList(7) = 'PrincipalStrains'
      MEF90DefMech_unilateralContactTypeList(8) = 'Masonry'
      MEF90DefMech_unilateralContactTypeList(9) = 'MEF90DefMech_unilateralContactTypeList'
      MEF90DefMech_unilateralContactTypeList(10) = '_MEF90DefMech_unilateralContactTypeList'
      MEF90DefMech_unilateralContactTypeList(11) = ''

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
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechCtxCreate(DefMechCtx,Mesh,MEF90Ctx,ierr)
      Type(MEF90DefMechCtx_Type),Intent(OUT)                   :: DefMechCtx
      Type(DM),target,Intent(IN)                               :: Mesh
      Type(MEF90Ctx_Type),target,Intent(IN)                    :: MEF90Ctx
      PetscErrorCode,Intent(OUT)                               :: ierr
   
      Type(MEF90CtxGlobalOptions_Type),pointer                 :: MEF90CtxGlobalOptions
      Type(MEF90DefMechGlobalOptions_Type),pointer             :: MEF90DefMechGlobalOptions
      Type(MEF90DefMechCellSetOptions_Type),pointer            :: MEF90DefMechCellSetOptions
      Type(MEF90DefMechVertexSetOptions_Type),pointer          :: MEF90DefMechVertexSetOptions
      Type(SectionReal)                                        :: defaultSection
      Type(IS)                                                 :: setIS
      PetscInt                                                 :: set,numSet
      PetscInt                                                 :: dim
      Character(len=MEF90_MXSTRLEN)                            :: filename,IOBuffer

      Call MEF90DefMechCtxInitialize_Private(ierr)
      Call DMMeshGetDimension(Mesh,dim,ierr);CHKERRQ(ierr)
      DefMechCtx%DM => Mesh
      DefMechCtx%MEF90Ctx => MEF90Ctx
      
      !!! Clone DM for each of the other layouts
      Call DMMeshClone(Mesh,DefMechCtx%cellDMVect,ierr);CHKERRQ(ierr)
      Call DMMeshSetMaxDof(DefMechCtx%cellDMVect,dim,ierr);CHKERRQ(ierr) 
      Call DMSetBlockSize(DefMechCtx%cellDMVect,dim,ierr);CHKERRQ(ierr)
      
      Call DMMeshClone(Mesh,DefMechCtx%cellDMScal,ierr);CHKERRQ(ierr)
      Call DMMeshSetMaxDof(DefMechCtx%cellDMScal,1,ierr);CHKERRQ(ierr) 
      Call DMSetBlockSize(DefMechCtx%cellDMScal,1,ierr);CHKERRQ(ierr)

      Call DMMeshClone(Mesh,DefMechCtx%DMVect,ierr);CHKERRQ(ierr)
      Call DMMeshSetMaxDof(DefMechCtx%DMVect,dim,ierr);CHKERRQ(ierr) 
      Call DMSetBlockSize(DefMechCtx%DMVect,dim,ierr);CHKERRQ(ierr)

      Call DMMeshClone(Mesh,DefMechCtx%DMScal,ierr);CHKERRQ(ierr)
      Call DMMeshSetMaxDof(DefMechCtx%DMScal,1,ierr);CHKERRQ(ierr) 
      Call DMSetBlockSize(DefMechCtx%DMScal,1,ierr);CHKERRQ(ierr)

      Call DMMeshClone(Mesh,DefMechCtx%DMMatS,ierr);CHKERRQ(ierr)
      Call DMMeshSetMaxDof(DefMechCtx%DMMatS,(dim*(dim+1))/2,ierr);CHKERRQ(ierr) 
      Call DMSetBlockSize(DefMechCtx%DMMatS,(dim*(dim+1))/2,ierr);CHKERRQ(ierr)

      Call DMMeshClone(Mesh,DefMechCtx%cellDMMatS,ierr);CHKERRQ(ierr)
      Call DMMeshSetMaxDof(DefMechCtx%cellDMMatS,(dim*(dim+1))/2,ierr);CHKERRQ(ierr) 
      Call DMSetBlockSize(DefMechCtx%cellDMMatS,(dim*(dim+1))/2,ierr);CHKERRQ(ierr)

      Call PetscBagCreate(MEF90Ctx%comm,sizeofMEF90DefMechGlobalOptions,DefMechCtx%GlobalOptionsBag,ierr);CHKERRQ(ierr)
      
      !!! Call DMmeshGetLabelSize(Mesh,'Cell Sets',numSet,ierr);CHKERRQ(ierr)
      !!! NO: I need to allocate for the overall number of sets, not the local one

      Call DMmeshGetLabelIdIS(Mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr) 
      Call ISGetLocalSize(setIS,numSet,ierr);CHKERRQ(ierr)
      Allocate(DefMechCtx%CellSetOptionsBag(numSet),stat=ierr)
      Do set = 1, numSet
         Call PetscBagCreate(MEF90Ctx%comm,sizeofMEF90DefMechCellSetOptions,DefMechCtx%CellSetOptionsBag(set),ierr);CHKERRQ(ierr)
      End Do
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(Mesh,'Vertex Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr) 
      Call ISGetLocalSize(setIS,numSet,ierr);CHKERRQ(ierr)
      Allocate(DefMechCtx%VertexSetOptionsBag(numSet),stat=ierr)
      Do set = 1, numSet
         Call PetscBagCreate(MEF90Ctx%comm,sizeofMEF90DefMechVertexSetOptions,DefMechCtx%VertexSetOptionsBag(set),ierr);CHKERRQ(ierr)
      End Do
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)


      !!! 
      !!! Create energy viewers
      !!!      
      filename = trim(MEF90FilePrefix(MEF90Ctx%resultFile))//'.ener'
      Call PetscViewerASCIIOpen(MEF90Ctx%comm,filename,DefMechCtx%globalEnergyViewer,ierr);CHKERRQ(ierr)
      Call PetscViewerASCIIPrintf(DefMechCtx%globalEnergyViewer,"# step     load            elastic energy  work            cohesive energy surface energy  total energy   dissipation plastic \n",ierr);CHKERRQ(ierr)
      Call PetscViewerFlush(DefMechCtx%globalEnergyViewer,ierr);CHKERRQ(ierr)
      
      Call DMmeshGetLabelIdIS(Mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr) 
      Call ISGetLocalSize(setIS,numSet,ierr);CHKERRQ(ierr)
      Allocate(DefMechCtx%setEnergyViewer(numSet),stat=ierr)
      Do set = 1, numSet
         Write(filename,101) trim(MEF90FilePrefix(MEF90Ctx%resultFile)),set
         Call PetscViewerASCIIOpen(MEF90Ctx%comm,filename,DefMechCtx%setEnergyViewer(set),ierr);CHKERRQ(ierr)
         Write(IOBuffer,102) set
         Call PetscViewerASCIIPrintf(DefMechCtx%setEnergyViewer(set),IOBuffer,ierr);CHKERRQ(ierr)
         Call PetscViewerASCIIPrintf(DefMechCtx%setEnergyViewer(set),"# step     load            elastic energy  work            cohesive energy surface energy  total energy   dissipation plastic \n",ierr);CHKERRQ(ierr)
         Call PetscViewerFlush(DefMechCtx%setEnergyViewer(set),ierr);CHKERRQ(ierr)
      End Do
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
101 Format(A,'-',I4.4,'.enerblk')
102 Format("# cell set ",I4,"\n")
      DefMechCtx%analysisTime = 0.0_Kr
      DefMechCtx%timeStep = 0.0_Kr
      Nullify(DefMechCtx%force)
      Nullify(DefMechCtx%pressureforce)
      Nullify(DefMechCtx%CrackPressure)
      Nullify(DefMechCtx%boundaryDisplacement)
      Nullify(DefMechCtx%displacementLowerBound)
      Nullify(DefMechCtx%displacementUpperBound)
      Nullify(DefMechCtx%boundaryDamage)
      Nullify(DefMechCtx%Displacement)
      Nullify(DefMechCtx%displacementPreviousStep)
      Nullify(DefMechCtx%Damage)
      Nullify(DefMechCtx%damagePreviousStep)
      Nullify(DefMechCtx%temperature)
      Nullify(DefMechCtx%plasticStrain)
      Nullify(DefMechCtx%cumulatedDissipatedPlasticEnergy)
   End Subroutine MEF90DefMechCtxCreate
   
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCtxSetSections"
!!!
!!!  
!!!  MEF90DefMechCtxSetSections: Set the data layout for each of the fields involved in a MEF90DefMechCtx_Type
!!!                               Uses Sieve convenience functions for now, but will be pulling layout informations
!!!                               from the element types when switching to DMComplex
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechCtxSetSections(DefMechCtx,ierr)
      Type(MEF90DefMechCtx_Type),Intent(INOUT)        :: DefMechCtx
      PetscErrorCode,Intent(OUT)                      :: ierr
   
      PetscInt                                        :: dim

      Call DMMeshGetDimension(DefMechCtx%DM,dim,ierr);CHKERRQ(ierr)

      Call DMMeshGetVertexSectionReal(DefMechCtx%DMVect,"default",dim,DefMechCtx%DMVectSec,ierr);CHKERRQ(ierr)
      Call DMMeshSetSectionReal(DefMechCtx%DMVect,"default",DefMechCtx%DMVectSec,ierr);CHKERRQ(ierr)
      Call DMSetBlockSize(DefMechCtx%DMVect,dim,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(DefMechCtx%DMVect,DefMechCtx%DMVectSec,DefMechCtx%DMVectScatter,ierr);CHKERRQ(ierr)
   
      Call DMMeshGetCellSectionReal(DefMechCtx%cellDMVect,"default",dim,DefMechCtx%cellDMVectSec,ierr);CHKERRQ(ierr)
      Call DMMeshSetSectionReal(DefMechCtx%cellDMVect,"default",DefMechCtx%cellDMVectSec,ierr);CHKERRQ(ierr)
      Call DMSetBlockSize(DefMechCtx%CellDMVect,dim,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(DefMechCtx%cellDMVect,DefMechCtx%cellDMVectSec,DefMechCtx%cellDMVectScatter,ierr);CHKERRQ(ierr)
      
      Call DMMeshGetVertexSectionReal(DefMechCtx%DMScal,"default",1,DefMechCtx%DMScalSec,ierr);CHKERRQ(ierr)
      Call DMMeshSetSectionReal(DefMechCtx%DMScal,"default",DefMechCtx%DMScalSec,ierr);CHKERRQ(ierr)
      Call DMSetBlockSize(DefMechCtx%DMScal,1,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(DefMechCtx%DMScal,DefMechCtx%DMScalSec,DefMechCtx%DMScalScatter,ierr);CHKERRQ(ierr)

      Call DMMeshGetCellSectionReal(DefMechCtx%cellDMScal,"default",1,DefMechCtx%cellDMScalSec,ierr);CHKERRQ(ierr)
      Call DMMeshSetSectionReal(DefMechCtx%cellDMScal,"default",DefMechCtx%cellDMScalSec,ierr);CHKERRQ(ierr)
      Call DMSetBlockSize(DefMechCtx%cellDMScal,1,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(DefMechCtx%cellDMScal,DefMechCtx%cellDMScalSec,DefMechCtx%cellDMScalScatter,ierr);CHKERRQ(ierr)

      Call DMMeshGetVertexSectionReal(DefMechCtx%DMMatS,"default",(dim*(dim+1))/2,DefMechCtx%DMMatSSec,ierr);CHKERRQ(ierr)
      Call DMMeshSetSectionReal(DefMechCtx%DMMatS,"default",DefMechCtx%DMMatSSec,ierr);CHKERRQ(ierr)
      Call DMSetBlockSize(DefMechCtx%DMMatS,(dim*(dim+1))/2,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(DefMechCtx%DMMatS,DefMechCtx%DMMatSSec,DefMechCtx%DMMatSScatter,ierr);CHKERRQ(ierr)

      Call DMMeshGetCellSectionReal(DefMechCtx%cellDMMatS,"default",(dim*(dim+1))/2,DefMechCtx%cellDMMatSSec,ierr);CHKERRQ(ierr)
      Call DMMeshSetSectionReal(DefMechCtx%cellDMMatS,"default",DefMechCtx%cellDMMatSSec,ierr);CHKERRQ(ierr)
      Call DMSetBlockSize(DefMechCtx%cellDMMatS,(dim*(dim+1))/2,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(DefMechCtx%cellDMMatS,DefMechCtx%cellDMMatSSec,DefMechCtx%cellDMMatSScatter,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechCtxSetSections

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCtxCreateVectors"
!!!
!!!  
!!!  MEF90DefMechCtxCreateVectors: Create a default set of vectors in a MEF90DefMechCtx_Type
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechCtxCreateVectors(DefMechCtx,ierr)
      Type(MEF90DefMechCtx_Type),Intent(INOUT)        :: DefMechCtx
      PetscErrorCode,Intent(OUT)                      :: ierr
      
      Allocate(DefMechCtx%Displacement,stat=ierr)
      Call DMCreateGlobalVector(DefMechCtx%DMVect,DefMechCtx%Displacement,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(DefMechCtx%Displacement,"Displacement",ierr);CHKERRQ(ierr)
      Call VecSet(DefMechCtx%Displacement,0.0_Kr,ierr);CHKERRQ(ierr)

      Allocate(DefMechCtx%displacementPreviousStep,stat=ierr)
      Call DMCreateGlobalVector(DefMechCtx%DMVect,DefMechCtx%displacementPreviousStep,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(DefMechCtx%displacementPreviousStep,"displacementPreviousStep",ierr);CHKERRQ(ierr)
      Call VecSet(DefMechCtx%displacementPreviousStep,0.0_Kr,ierr);CHKERRQ(ierr)

      Allocate(DefMechCtx%boundaryDisplacement,stat=ierr)
      Call DMCreateGlobalVector(DefMechCtx%DMVect,DefMechCtx%boundaryDisplacement,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(DefMechCtx%boundaryDisplacement,"boundary Displacement",ierr);CHKERRQ(ierr)
      Call VecSet(DefMechCtx%boundaryDisplacement,0.0_Kr,ierr);CHKERRQ(ierr)
   
      Allocate(DefMechCtx%displacementLowerBound,stat=ierr)
      Call DMCreateGlobalVector(DefMechCtx%DMVect,DefMechCtx%displacementLowerBound,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(DefMechCtx%displacementLowerBound,"displacement Lower Bound",ierr);CHKERRQ(ierr)
      Call VecSet(DefMechCtx%displacementLowerBound,0.0_Kr,ierr);CHKERRQ(ierr)
   
      Allocate(DefMechCtx%displacementUpperBound,stat=ierr)
      Call DMCreateGlobalVector(DefMechCtx%DMVect,DefMechCtx%displacementUpperBound,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(DefMechCtx%displacementUpperBound,"displacement Upper Bound",ierr);CHKERRQ(ierr)
      Call VecSet(DefMechCtx%displacementUpperBound,0.0_Kr,ierr);CHKERRQ(ierr)
   
      Allocate(DefMechCtx%force,stat=ierr)
      Call DMCreateGlobalVector(DefMechCtx%cellDMVect,DefMechCtx%force,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(DefMechCtx%force,"Force",ierr);CHKERRQ(ierr)
      Call VecSet(DefMechCtx%force,0.0_Kr,ierr);CHKERRQ(ierr)

      Allocate(DefMechCtx%pressureForce,stat=ierr)
      Call DMCreateGlobalVector(DefMechCtx%cellDMScal,DefMechCtx%pressureForce,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(DefMechCtx%pressureForce,"PressureForce",ierr);CHKERRQ(ierr)
      Call VecSet(DefMechCtx%pressureForce,0.0_Kr,ierr);CHKERRQ(ierr)

      Allocate(DefMechCtx%CrackPressure,stat=ierr)
      Call DMCreateGlobalVector(DefMechCtx%cellDMScal,DefMechCtx%CrackPressure,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(DefMechCtx%CrackPressure,"CrackPressure",ierr);CHKERRQ(ierr)
      Call VecSet(DefMechCtx%CrackPressure,0.0_Kr,ierr);CHKERRQ(ierr)

      Allocate(DefMechCtx%temperature,stat=ierr)
      Call DMCreateGlobalVector(DefMechCtx%DMScal,DefMechCtx%temperature,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(DefMechCtx%temperature,"temperature",ierr);CHKERRQ(ierr)
      Call VecSet(DefMechCtx%temperature,0.0_Kr,ierr);CHKERRQ(ierr)

      Allocate(DefMechCtx%plasticStrain,stat=ierr)
      Call DMCreateGlobalVector(DefMechCtx%CellDMMatS,DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(DefMechCtx%plasticStrain,"plasticStrain",ierr);CHKERRQ(ierr)
      Call VecSet(DefMechCtx%plasticStrain,0.0_Kr,ierr);CHKERRQ(ierr)

      Allocate(DefMechCtx%cumulatedDissipatedPlasticEnergy,stat=ierr)
      Call DMCreateGlobalVector(DefMechCtx%cellDMScal,DefMechCtx%cumulatedDissipatedPlasticEnergy,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(DefMechCtx%cumulatedDissipatedPlasticEnergy,"cumulatedDissipatedPlasticEnergy",ierr);CHKERRQ(ierr)
      Call VecSet(DefMechCtx%cumulatedDissipatedPlasticEnergy,0.0_Kr,ierr);CHKERRQ(ierr)
   
      Allocate(DefMechCtx%damage,stat=ierr)
      Call DMCreateGlobalVector(DefMechCtx%DMScal,DefMechCtx%damage,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(DefMechCtx%damage,"damage",ierr);CHKERRQ(ierr)
      Call VecSet(DefMechCtx%damage,0.0_Kr,ierr);CHKERRQ(ierr)

      Allocate(DefMechCtx%damagePreviousStep,stat=ierr)
      Call DMCreateGlobalVector(DefMechCtx%DMScal,DefMechCtx%damagePreviousStep,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(DefMechCtx%damagePreviousStep,"damagePreviousStep",ierr);CHKERRQ(ierr)
      Call VecSet(DefMechCtx%damagePreviousStep,0.0_Kr,ierr);CHKERRQ(ierr)

      Allocate(DefMechCtx%boundaryDamage,stat=ierr)
      Call DMCreateGlobalVector(DefMechCtx%DMScal,DefMechCtx%boundaryDamage,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(DefMechCtx%boundaryDamage,"boundaryDamage",ierr);CHKERRQ(ierr)
      Call VecSet(DefMechCtx%boundaryDamage,0.0_Kr,ierr);CHKERRQ(ierr)

      Allocate(DefMechCtx%stress,stat=ierr)
      Call DMCreateGlobalVector(DefMechCtx%CellDMMatS,DefMechCtx%stress,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(DefMechCtx%stress,"stress",ierr);CHKERRQ(ierr)
      Call VecSet(DefMechCtx%stress,0.0_Kr,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechCtxCreateVectors

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCtxDestroyVectors"
!!!
!!!  
!!!  MEF90DefMechCtxDestroyVectors: destroys the Vecs in a MEF90DefMechCtx_Type
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechCtxDestroyVectors(DefMechCtx,ierr)
      Type(MEF90DefMechCtx_Type),Intent(INOUT)        :: DefMechCtx
      PetscErrorCode,Intent(OUT)                      :: ierr

      If (Associated(DefMechCtx%Displacement)) Then 
         Call VecDestroy(DefMechCtx%Displacement,ierr);CHKERRQ(ierr)
         DeAllocate(DefMechCtx%Displacement)
         Nullify(DefMechCtx%Displacement)
      End If   

      If (Associated(DefMechCtx%displacementPreviousStep)) Then 
         Call VecDestroy(DefMechCtx%displacementPreviousStep,ierr);CHKERRQ(ierr)
         DeAllocate(DefMechCtx%displacementPreviousStep)
         Nullify(DefMechCtx%displacementPreviousStep)
      End If   

      If (Associated(DefMechCtx%boundaryDisplacement)) Then 
         Call VecDestroy(DefMechCtx%boundaryDisplacement,ierr);CHKERRQ(ierr)
         DeAllocate(DefMechCtx%boundaryDisplacement)
         Nullify(DefMechCtx%boundaryDisplacement)
      End If   

      If (Associated(DefMechCtx%displacementLowerBound)) Then 
         Call VecDestroy(DefMechCtx%displacementLowerBound,ierr);CHKERRQ(ierr)
         DeAllocate(DefMechCtx%displacementLowerBound)
         Nullify(DefMechCtx%displacementLowerBound)
      End If   

      If (Associated(DefMechCtx%displacementUpperBound)) Then 
         Call VecDestroy(DefMechCtx%displacementUpperBound,ierr);CHKERRQ(ierr)
         DeAllocate(DefMechCtx%displacementUpperBound)
         Nullify(DefMechCtx%displacementUpperBound)
      End If   

      If (Associated(DefMechCtx%boundaryDamage)) Then 
         Call VecDestroy(DefMechCtx%boundaryDamage,ierr);CHKERRQ(ierr)
         DeAllocate(DefMechCtx%boundaryDamage)
         Nullify(DefMechCtx%boundaryDamage)
      End If
   
      If (Associated(DefMechCtx%damage)) Then 
         Call VecDestroy(DefMechCtx%damage,ierr);CHKERRQ(ierr)
         DeAllocate(DefMechCtx%damage)
         Nullify(DefMechCtx%damage)
      End If
   
      If (Associated(DefMechCtx%damagePreviousStep)) Then 
         Call VecDestroy(DefMechCtx%damagePreviousStep,ierr);CHKERRQ(ierr)
         DeAllocate(DefMechCtx%damagePreviousStep)
         Nullify(DefMechCtx%damagePreviousStep)
      End If
   
      If (Associated(DefMechCtx%force)) Then 
         Call VecDestroy(DefMechCtx%force,ierr);CHKERRQ(ierr)   
         DeAllocate(DefMechCtx%force)
         Nullify(DefMechCtx%force)
      End If

      If (Associated(DefMechCtx%pressureForce)) Then 
         Call VecDestroy(DefMechCtx%pressureForce,ierr);CHKERRQ(ierr)      
         DeAllocate(DefMechCtx%pressureForce)
         Nullify(DefMechCtx%pressureForce)
      End If

      If (Associated(DefMechCtx%CrackPressure)) Then 
         Call VecDestroy(DefMechCtx%CrackPressure,ierr);CHKERRQ(ierr)      
         DeAllocate(DefMechCtx%CrackPressure)
         Nullify(DefMechCtx%CrackPressure)
      End If

      If (Associated(DefMechCtx%temperature)) Then 
         Call VecDestroy(DefMechCtx%temperature,ierr);CHKERRQ(ierr)   
         DeAllocate(DefMechCtx%temperature)
         Nullify(DefMechCtx%temperature)
      End If

      If (Associated(DefMechCtx%plasticStrain)) Then 
         Call VecDestroy(DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)   
         DeAllocate(DefMechCtx%plasticStrain)
         Nullify(DefMechCtx%plasticStrain)
      End If

      If (Associated(DefMechCtx%cumulatedDissipatedPlasticEnergy)) Then 
         Call VecDestroy(DefMechCtx%cumulatedDissipatedPlasticEnergy,ierr);CHKERRQ(ierr)   
         DeAllocate(DefMechCtx%cumulatedDissipatedPlasticEnergy)
         Nullify(DefMechCtx%cumulatedDissipatedPlasticEnergy)
      End If

      If (Associated(DefMechCtx%stress)) Then 
         Call VecDestroy(DefMechCtx%stress,ierr);CHKERRQ(ierr)   
         DeAllocate(DefMechCtx%stress)
         Nullify(DefMechCtx%stress)
      End If
   End Subroutine MEF90DefMechCtxDestroyVectors

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCtxDestroy"
!!!
!!!  
!!!  MEF90DefMechCtxDestroy: destroys a MEF90DefMechCtx_Type
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechCtxDestroy(DefMechCtx,ierr)
      Type(MEF90DefMechCtx_Type),Intent(OUT)          :: DefMechCtx
      PetscErrorCode,Intent(OUT)                      :: ierr
      
      PetscInt                                        :: set
   
      Call PetscBagDestroy(DefMechCtx%GlobalOptionsBag,ierr);CHKERRQ(ierr)
      Do set = 1, size(DefMechCtx%CellSetOptionsBag)
         Call PetscBagDestroy(DefMechCtx%CellSetOptionsBag(set),ierr);CHKERRQ(ierr)
      End Do
      DeAllocate(DefMechCtx%CellSetOptionsBag)
      Do set = 1, size(DefMechCtx%VertexSetOptionsBag)
         Call PetscBagDestroy(DefMechCtx%VertexSetOptionsBag(set),ierr);CHKERRQ(ierr)
      End Do
      DeAllocate(DefMechCtx%VertexSetOptionsBag)
      
      !!! 
      !!! Close energy viewers
      !!!      
      Call PetscViewerDestroy(DefMechCtx%globalEnergyViewer,ierr);CHKERRQ(ierr)
      Do set = 1, size(DefMechCtx%setEnergyViewer)
         Call PetscViewerDestroy(DefMechCtx%setEnergyViewer(set),ierr);CHKERRQ(ierr)
      End Do
      DeAllocate(DefMechCtx%setEnergyViewer)

      !!! 
      !!! We only nullify the Vecs since they may be shared with other contexts. 
      !!!
      Nullify(DefMechCtx%force) 
      Nullify(DefMechCtx%pressureforce)
      Nullify(DefMechCtx%CrackPressure)
      Nullify(DefMechCtx%boundaryDisplacement)
      Nullify(DefMechCtx%displacementLowerBound)
      Nullify(DefMechCtx%displacementUpperBound)
      Nullify(DefMechCtx%boundaryDamage)
      Nullify(DefMechCtx%Displacement)
      Nullify(DefMechCtx%displacementPreviousStep)
      Nullify(DefMechCtx%Damage)
      Nullify(DefMechCtx%damagePreviousStep)
      Nullify(DefMechCtx%temperature)
      Nullify(DefMechCtx%plasticStrain)
      Nullify(DefMechCtx%cumulatedDissipatedPlasticEnergy)

      Call VecScatterDestroy(DefMechCtx%DMScalScatter,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(DefMechCtx%cellDMScalScatter,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(DefMechCtx%DMVectScatter,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(DefMechCtx%cellDMVectScatter,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(DefMechCtx%DMMatSScatter,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(DefMechCtx%cellDMMatSScatter,ierr);CHKERRQ(ierr)

      !Call SectionRealDestroy(DefMechCtx%DMScalSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(DefMechCtx%cellDMScalSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(DefMechCtx%DMVectSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(DefMechCtx%cellDMVectSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(DefMechCtx%DMMatSSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(DefMechCtx%cellDMMatSSec,ierr);CHKERRQ(ierr)

      Call DMDestroy(DefMechCtx%DMScal,ierr);CHKERRQ(ierr)
      Call DMDestroy(DefMechCtx%cellDMScal,ierr);CHKERRQ(ierr)
      Call DMDestroy(DefMechCtx%DMVect,ierr);CHKERRQ(ierr)
      !Call DMDestroy(DefMechCtx%cellDMVect,ierr);CHKERRQ(ierr)
      Call DMDestroy(DefMechCtx%DMMatS,ierr);CHKERRQ(ierr)
      Call DMDestroy(DefMechCtx%cellDMMatS,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechCtxDestroy

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90DefMechCtxGlobalOptions"
!!!
!!!  
!!!  PetscBagRegisterMEF90DefMechCtxGlobalOptions:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine PetscBagRegisterMEF90DefMechCtxGlobalOptions(bag,name,prefix,default,ierr)
      PetscBag                                                 :: bag
      Character(len=*),Intent(IN)                              :: prefix,name
      Type(MEF90DefMechGlobalOptions_Type),Intent(IN)          :: default
      PetscErrorCode,Intent(OUT)                               :: ierr

      Type(MEF90DefMechGlobalOptions_Type),pointer      :: DefMechGlobalOptions

      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(bag,DefMechGlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagSetName(bag,trim(name),"DefMechGlobalOptions MEF90 Defect Mechanics global options",ierr);CHKERRQ(ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix),ierr);CHKERRQ(ierr)

      Call PetscBagRegisterEnum(bag,DefMechGlobalOptions%timeSteppingType,MEF90DefMech_TimeSteppingTypeList,default%timeSteppingType,'DefMech_TimeStepping_Type','Type of defect mechanics Time steping',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterEnum(bag,DefMechGlobalOptions%solverType,MEF90DefMech_SolverTypeList,default%solverType,'DefMech_solver_Type','Type of defect mechanics solver',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterBool(bag,DefMechGlobalOptions%addDisplacementNullSpace,default%addDisplacementNullSpace,'disp_addNullSpace','Add null space to SNES',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,DefMechGlobalOptions%displacementOffset,default%displacementOffset,'displacement_Offset','Position of displacement field in EXO file',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,DefMechGlobalOptions%damageOffset,default%damageOffset,'damage_Offset','Position of damage field in EXO file',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,DefMechGlobalOptions%boundaryDamageOffset,default%boundaryDamageOffset,'boundaryDamage_Offset','Position of boundary damage field in EXO file',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,DefMechGlobalOptions%stressOffset,default%stressOffset,'stress_Offset','Position of stress field in EXO file',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,DefMechGlobalOptions%temperatureOffset,default%temperatureOffset,'temperature_Offset','Position of temperature field in EXO file',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,DefMechGlobalOptions%plasticStrainOffset,default%plasticStrainOffset,'plasticStrain_Offset','Position of the plastic strain field in EXO file',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,DefMechGlobalOptions%cumulatedPlasticDissipationOffset,default%cumulatedPlasticDissipationOffset,'cumulatedPlasticDissipation_Offset','Position of the Cumulated Plastic Plastic Dissipation field in EXO file',ierr);CHKERRQ(ierr)

      Call PetscBagRegisterEnum(bag,DefMechGlobalOptions%boundaryDisplacementScaling,MEF90ScalingList,default%boundaryDisplacementScaling,'boundaryDisplacement_scaling','Boundary displacement scaling',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterEnum(bag,DefMechGlobalOptions%displacementLowerBoundScaling,MEF90ScalingList,default%displacementLowerBoundScaling,'displacementlowerbound_scaling','Displacement lower bound scaling',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterEnum(bag,DefMechGlobalOptions%displacementUpperBoundScaling,MEF90ScalingList,default%displacementUpperBoundScaling,'displacementupperbound_scaling','Displacement upper bound scaling',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,DefMechGlobalOptions%boundaryDisplacementOffset,default%boundaryDisplacementOffset,'boundaryDisplacement_Offset','Position of boundary displacement field in EXO file',ierr);CHKERRQ(ierr)

      Call PetscBagRegisterEnum(bag,DefMechGlobalOptions%boundaryDamageScaling,MEF90ScalingList,default%boundaryDamageScaling,'boundaryDamage_scaling','Boundary damage scaling',ierr);CHKERRQ(ierr)

      Call PetscBagRegisterEnum(bag,DefMechGlobalOptions%forceScaling,MEF90ScalingList,default%forceScaling,'force_scaling','Force scaling',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,DefMechGlobalOptions%forceOffset,default%forceOffset,'force_Offset','Position of force field in EXO file',ierr);CHKERRQ(ierr)

      Call PetscBagRegisterEnum(bag,DefMechGlobalOptions%pressureForceScaling,MEF90ScalingList,default%pressureforceScaling,'pressureForce_scaling','Pressure force scaling',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,DefMechGlobalOptions%pressureForceOffset,default%pressureForceOffset,'pressureForce_Offset','Position of pressure force field in EXO file',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterEnum(bag,DefMechGlobalOptions%CrackPressureScaling,MEF90ScalingList,default%CrackPressureScaling,'CrackPressure_scaling','Crack Pressure scaling',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,DefMechGlobalOptions%CrackPressureOffset,default%CrackPressureOffset,'CrackPressure_Offset','Position of Crack Pressure field in EXO file',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,DefMechGlobalOptions%damageATol,default%damageATol,'defmech_damage_atol','Absolute tolerance on damage error',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,DefMechGlobalOptions%maxit,default%maxit,'defmech_maxit','Maximum number of alternate minimizations for damage',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,DefMechGlobalOptions%PCLag,default%PCLag,'defmech_pclag','Interval at which the PC is recomputed during alternate minimization',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,DefMechGlobalOptions%irrevthres,default%irrevthres,'defmech_irrevThres','Threshold above which irreversibility is enforced (0 for monotonicity, .99 for equality)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,DefMechGlobalOptions%SOROmega,default%SOROmega,'defmech_SOR_Omega','Alterate Minimization over relaxation factor (>0 for limited, <0 for projected) ',ierr);CHKERRQ(ierr)

      Call PetscBagRegisterEnum(bag,DefMechGlobalOptions%BTType,MEF90DefMech_BTTypeList,default%BTType,'BT_Type','Backtracking type',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,DefMechGlobalOptions%BTInterval,default%BTInterval,'BT_Interval','Interval at which Backtracking is run in inner loop (0 for outer loop)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,DefMechGlobalOptions%BTScope,default%BTScope,'BT_Scope','Backtracking scope (0 for unlimited)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,DefMechGlobalOptions%BTTol,default%BTTol,'BT_Tol','Backtracking relative tolerance',ierr);CHKERRQ(ierr)

      Call PetscBagRegisterReal(bag,DefMechGlobalOptions%plasticStrainATol,default%plasticStrainATol,'defmech_plasticstrain_atol','Absolute tolerance on plastic error',ierr);CHKERRQ(ierr)
      
      Call PetscBagRegisterReal(bag,DefMechGlobalOptions%InjectedVolumeATol,default%InjectedVolumeATol,'defmech_InjectedVolume_atol','Absolute tolerance on injected volume error',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,DefMechGlobalOptions%dampingCoefficientDisplacement,default%dampingCoefficientDisplacement,'defmech_dampingCoefficient_displacement','Damping coefficient on displacement field (0 for minimization, 1 for semi-implicit gradient flow)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,DefMechGlobalOptions%dampingCoefficientDamage,default%dampingCoefficientDamage,'defmech_dampingCoefficient_damage','Damping coefficient on damage field (0 for minimization, 1 for semi-implicit gradient flow)',ierr);CHKERRQ(ierr)
   End Subroutine PetscBagRegisterMEF90DefMechCtxGlobalOptions

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90DefMechCtxCellSetOptions"
!!!
!!!  
!!!  PetscBagRegisterMEF90DefMechCtxCellSetOptions:
!!!  
!!!  (c) 2012-18 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine PetscBagRegisterMEF90DefMechCtxCellSetOptions(bag,name,prefix,default,ierr)
      PetscBag                                           :: bag
      Character(len=*),Intent(IN)                        :: prefix,name
      Type(MEF90DefMechCellSetOptions_Type),Intent(IN)   :: default
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(MEF90DefMechCellSetOptions_Type),pointer      :: DefMechCellSetOptions
      
      Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(bag,DefMechCellSetOptions,ierr);CHKERRQ(ierr)
      Call PetscBagSetName(bag,trim(name),"DefMechCellSetOptions MEF90 Defect Mechanics Cell Set options",ierr);CHKERRQ(ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix),ierr);CHKERRQ(ierr)
      
      DefMechCellSetOptions%force = default%Force
      DefMechCellSetOptions%boundaryDisplacement   = default%boundaryDisplacement
      DefMechCellSetOptions%displacementLowerBound = default%displacementLowerBound
      DefMechCellSetOptions%displacementUpperBound = default%displacementUpperBound
      DefMechCellSetOptions%Has_displacementBC     = default%Has_displacementBC

      Call PetscBagRegisterInt(bag,DefMechCellSetOptions%ElemTypeShortIDDisplacement,default%ElemTypeShortIDDisplacement,'ShortIDDisplacement','Displacement element type ShortID',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt(bag,DefMechCellSetOptions%ElemTypeShortIDDamage,default%ElemTypeShortIDDamage,'ShortIDDamage','Damage field element type ShortID',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterRealArray(bag,DefMechCellSetOptions%force,3,'Force','[N.m^(-3) / N.m^(-2) / N.m^(-1)] (f): body / boundary force',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,DefMechCellSetOptions%pressureForce,default%pressureForce,'pressureForce','[N.m^(-2) / N.m^(-1)] (p): boundary pressureforce',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,DefMechCellSetOptions%CrackPressure,default%CrackPressure,'CrackPressure','without unit: internal crack pressure',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterEnum(bag,DefMechCellSetOptions%damageType,MEF90DefMech_damageTypeList,default%damageType,'damage_type','Type of damage law',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterEnum(bag,DefMechCellSetOptions%plasticityType,MEF90DefMech_plasticityTypeList,default%plasticityType,'plasticity_type','Type of plasticity law',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterEnum(bag,DefMechCellSetOptions%unilateralContactType,MEF90DefMech_unilateralContactTypeList,default%unilateralContactType,'unilateralContact_type','Type of handling of unilateral contact',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterEnum(bag,DefMechCellSetOptions%drivingForceType,MEF90DefMech_drivingForceTypeList,default%drivingForceType,'drivingForce_type','Type of nucleation driving force',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterBoolArray(bag,DefMechCellSetOptions%Has_displacementBC,3,'DisplacementBC','Displacement has Dirichlet boundary Condition (Y/N)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterRealArray(bag,DefMechCellSetOptions%boundaryDisplacement,3,'boundaryDisplacement','[m] (U): Displacement boundary value',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterRealArray(bag,DefMechCellSetOptions%displacementLowerBound,3,'displacementLowerBound','[m] (U): Displacement lower bound',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterRealArray(bag,DefMechCellSetOptions%displacementUpperBound,3,'displacementUpperBound','[m] (U): Displacement upper bound',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterBool(bag,DefMechCellSetOptions%Has_DamageBC,default%Has_DamageBC,'DamageBC','Damage has Dirichlet boundary Condition (Y/N)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterBool(bag,DefMechCellSetOptions%CrackVolumeControlled,default%CrackVolumeControlled,'CrackVolumeControlled','Crack Pressure controlled by the crack volume in this block (Y/N)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterBool(bag,DefMechCellSetOptions%WorkControlled,default%WorkControlled,'WorkControlled','Force magnitude controlled by its work in this block (Y/N)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,DefMechCellSetOptions%boundaryDamage,default%boundaryDamage,'boundaryDamage','[unit-less] (alpha): Damage boundary value',ierr);CHKERRQ(ierr)
   End Subroutine PetscBagRegisterMEF90DefMechCtxCellSetOptions

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90DefMechCtxVertexSetOptions"
!!!
!!!  
!!!  PetscBagRegisterMEF90DefMechCtxVertexSetOptions:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine PetscBagRegisterMEF90DefMechCtxVertexSetOptions(bag,name,prefix,default,ierr)
      PetscBag                                              :: bag
      Character(len=*),Intent(IN)                           :: prefix,name
      Type(MEF90DefMechVertexSetOptions_Type),Intent(IN)    :: default
      PetscErrorCode,Intent(OUT)                            :: ierr

      Type(MEF90DefMechVertexSetOptions_Type),pointer       :: DefMechVertexSetOptions
      Call PetscBagGetDataMEF90DefMechCtxVertexSetOptions(bag,DefMechVertexSetOptions,ierr);CHKERRQ(ierr)
      Call PetscBagSetName(bag,trim(name),"DefMechVertexSetOptions MEF90 Defect Mechanics Vertex Set options",ierr);CHKERRQ(ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix),ierr);CHKERRQ(ierr)

      DefMechVertexSetOptions%Has_displacementBC     = default%Has_displacementBC
      DefMechVertexSetOptions%boundaryDisplacement   = default%boundaryDisplacement
      DefMechVertexSetOptions%displacementLowerBound = default%displacementLowerBound
      DefMechVertexSetOptions%displacementUpperBound = default%displacementUpperBound
      Call PetscBagRegisterBoolArray(bag,DefMechVertexSetOptions%Has_displacementBC,3,'DisplacementBC','Displacement has Dirichlet boundary Condition (Y/N)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterRealArray(bag,DefMechVertexSetOptions%boundaryDisplacement,3,'boundaryDisplacement','[m] (U): Displacement boundary value',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterRealArray(bag,DefMechVertexSetOptions%displacementLowerBound,3,'displacementLowerBound','[m] (U): Displacement lower bound',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterRealArray(bag,DefMechVertexSetOptions%displacementUpperBound,3,'displacementUpperBound','[m] (U): Displacement upper bound',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterBool(bag,DefMechVertexSetOptions%Has_DamageBC,default%Has_DamageBC,'DamageBC','Damage has Dirichlet boundary Condition (Y/N)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,DefMechVertexSetOptions%boundaryDamage,default%boundaryDamage,'boundaryDamage','[unit-less] (alpha): boundaryDamage',ierr);CHKERRQ(ierr)
   End Subroutine PetscBagRegisterMEF90DefMechCtxVertexSetOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCtxSetFromOptions"
!!!
!!!  
!!!  MEF90DefMechCtxSetFromOptions:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechCtxSetFromOptions(DefMechCtx,prefix,defaultGlobalOptions, &
                                              defaultCellSetOptions,    &
                                              defaultVertexSetOptions,ierr)
      Type(MEF90DefMechCtx_Type),Intent(OUT)                :: DefMechCtx
      Character(len=*),Intent(IN)                           :: prefix
      Type(MEF90DefMechGlobalOptions_Type),Intent(IN)       :: defaultGlobalOptions
      Type(MEF90DefMechCellSetOptions_Type),Intent(IN)      :: defaultCellSetOptions
      Character(len=MXSTLN),Dimension(:),Pointer            :: cellSetNames
      Type(MEF90DefMechVertexSetOptions_Type),Intent(IN)    :: defaultVertexSetOptions
      Character(len=MXSTLN),Dimension(:),Pointer            :: vertexSetNames
      PetscErrorCode,Intent(OUT)                            :: ierr
   
      Type(MEF90CtxGlobalOptions_Type),pointer              :: MEF90CtxGlobalOptions
      Type(MEF90DefMechCellSetOptions_Type)                 :: myDefaultCellSetOptions
      Type(MEF90Element_Type),Dimension(:),Pointer          :: ElemTypeScal,ElemTypeElast
      Type(IS)                                              :: setIS
      PetscInt,Dimension(:),Pointer                         :: setID
      PetscInt                                              :: set
      Character(len=MEF90_MXSTRLEN)                         :: IOBuffer,setName,setprefix
      Type(MEF90DefMechCellSetOptions_Type),pointer         :: cellSetOptions    
      Type(IS)                                              :: cellSetGlobalIS  

      Call PetscBagGetDataMEF90CtxGlobalOptions(DefMechCtx%MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr);CHKERRQ(ierr)
      !!!
      !!! Registering Global Context
      !!!
      Call PetscBagRegisterMEF90DefMechCtxGlobalOptions(DefMechCtx%GlobalOptionsBag,"MEF90DefMech Global Ctx",prefix,defaultGlobalOptions,ierr);CHKERRQ(ierr)

      If (MEF90CtxGlobalOptions%verbose > 0) Then
         Call PetscBagView(DefMechCtx%GlobalOptionsBag,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
      End If

      !!!
      !!! Registering Cell Set Context
      !!! We override the default element type with that detected from the exodus file
      !!!
      Call DMmeshGetLabelIdIS(DefMechCtx%DMVect,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(DefMechCtx%MEF90Ctx%comm,setIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      
      Call EXOGetCellSetElementType_Scal(DefMechCtx%MEF90Ctx,ElemTypeScal,ierr)
      Call EXOGetCellSetElementType_Elast(DefMechCtx%MEF90Ctx,ElemTypeElast,ierr)
      Do set = 1, size(setID)
         Write(setName,100) setID(set)
         Write(setprefix,101) setID(set)
         mydefaultCellSetOptions = defaultCellSetOptions
         mydefaultCellSetOptions%ElemTypeShortIDDisplacement = ElemTypeElast(set)%ShortID
         mydefaultCellSetOptions%ElemTypeShortIDDamage = ElemTypeScal(set)%ShortID
         Call PetscBagRegisterMEF90DefMechCtxCellSetOptions(DefMechCtx%CellSetOptionsBag(set),setName,setPrefix,mydefaultCellSetOptions,ierr)
         If (MEF90CtxGlobalOptions%verbose > 0) Then
            Write(IOBuffer,103) setID(set),trim(setprefix)
            Call PetscPrintf(DefMechCtx%MEF90Ctx%comm,IOBuffer,ierr);CHKERRQ(ierr)
            Call PetscBagView(DefMechCtx%CellSetOptionsBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
            Call PetscPrintf(DefMechCtx%MEF90Ctx%comm,"\n",ierr);CHKERRQ(ierr)
         End if
      End Do
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      DeAllocate(ElemTypeScal)
      DeAllocate(ElemTypeElast)
      
      !!!
      !!! Registering Vertex Set Context
      !!! We override the default element type with that detected from the exodus file
      !!!
      Call DMmeshGetLabelIdIS(DefMechCtx%DMVect,'Vertex Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      
      Do set = 1, size(setID)
         Write(setName,200) setID(set)
         Write(setprefix,201) setID(set)
         Call PetscBagRegisterMEF90DefMechCtxVertexSetOptions(DefMechCtx%VertexSetOptionsBag(set),setName,setPrefix,defaultVertexSetOptions,ierr)
         If (MEF90CtxGlobalOptions%verbose > 0) Then
            Write(IOBuffer,203) setID(set),trim(setprefix)
            Call PetscPrintf(DefMechCtx%MEF90Ctx%comm,IOBuffer,ierr);CHKERRQ(ierr)
            Call PetscBagView(DefMechCtx%VertexSetOptionsBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
            Call PetscPrintf(DefMechCtx%MEF90Ctx%comm,"\n",ierr);CHKERRQ(ierr)
         End if
      End Do
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)

      !!! Identify if any of the displacement bounds are set
      DefMechCtx%hasDisplacementBounds = PETSC_FALSE
      Call DMmeshGetLabelIdIS(DefMechCtx%CellDMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         DefMechCtx%hasDisplacementBounds = any(cellSetOptions%displacementLowerBound /= MEF90_NINFINITY) .OR. DefMechCtx%hasDisplacementBounds
         DefMechCtx%hasDisplacementBounds = any(cellSetOptions%displacementUpperBound /= MEF90_INFINITY)  .OR. DefMechCtx%hasDisplacementBounds
      EndDo
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)

100 Format('Cell set ',I4)
101 Format('cs',I4.4,'_')
103 Format('\nRegistering cell set ',I4,' prefix: ',A,'\n')
200 Format('Vertex set ',I4)
201 Format('vs',I4.4,'_')
203 Format('\nRegistering vertex set ',I4,' prefix: ',A,'\n')
   End Subroutine MEF90DefMechCtxSetFromOptions

End Module m_MEF90_DefMechCtx
