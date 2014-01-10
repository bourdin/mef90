#include "../MEF90/mef90.inc"
Module m_MEF90_DefMechCtx_Type
#include "finclude/petscdef.h"
   Use m_MEF90
   Use,Intrinsic :: iso_c_binding
   Implicit none
   !Private  
   !Public :: MEF90DefMechCtx_Type
   !Public :: MEF90DefMechGlobalOptions_Type
   !Public :: MEF90DefMechCellSetOptions_Type
   !Public :: MEF90DefMechVertexSetOptions_Type
   
   Type MEF90DefMechCtx_Type
      !!!  vertex based vec
      Type(Vec),pointer                :: displacement
      Type(Vec),pointer                :: damage
      Type(Vec),pointer                :: boundaryDisplacement
      Type(Vec),pointer                :: boundaryDamage
      Type(Vec),Pointer                :: temperature

      !!! cell based vec
      Type(Vec),pointer                :: force
      Type(Vec),pointer                :: pressureForce
      Type(Vec),Pointer                :: plasticStrain

      PetscBag                         :: GlobalOptionsBag
      PetscBag,Dimension(:),Pointer    :: CellSetOptionsBag
      PetscBag,Dimension(:),Pointer    :: VertexSetOptionsBag
      PetscBag,Dimension(:),Pointer    :: MaterialPropertiesBag
      Type(MEF90Ctx_Type),pointer      :: MEF90Ctx
      Type(DM),pointer                 :: DM
      Type(DM)                         :: DMScal,DMVect           ! Remove all these
      Type(DM)                         :: cellDMScal,cellDMVect   ! after switching to 
      Type(DM)                         :: DMMatS,cellDMMatS       ! DMcomplex
   End Type MEF90DefMechCtx_Type
   
   Type MEF90DefMechGlobalOptions_Type
      PetscInt                         :: mode
      PetscBool                        :: addDisplacementNullSpace
      !!! Position of vertex-based vecs in exo file
      PetscInt                         :: displacementOffset
      PetscInt                         :: damageOffset
      PetscInt                         :: boundaryDisplacementOffset
      PetscInt                         :: boundaryDamageOffset
      PetscInt                         :: temperatureOffset
      !!! Position of cell-based vecs in exo files
      PetscInt                         :: forceOffset
      PetscInt                         :: pressureForceOffset
      PetscInt                         :: plasticStrainOffset
      PetscInt                         :: stressOffset
      !!! scaling = time (step) scaling law currently CST, Linear, or File
      PetscInt                         :: boundaryDisplacementScaling
      PetscInt                         :: forceScaling
      PetscInt                         :: pressureForceScaling
   End Type MEF90DefMechGlobalOptions_Type

   Type MEF90DefMechCellSetOptions_Type
      PetscInt                         :: elemTypeShortIDDisplacement
      PetscInt                         :: elemTypeShortIDDamage
      PetscReal,Dimension(3)           :: force
      PetscReal                        :: pressureForce
      PetscEnum                        :: defectLaw
      PetscBool,Dimension(3)           :: Has_displacementBC
      PetscReal,Dimension(3)           :: boundaryDisplacement
      PetscBool                        :: Has_damageBC
      PetscReal                        :: boundaryDamage
   End Type MEF90DefMechCellSetOptions_Type

   Type MEF90DefMechVertexSetOptions_Type
      PetscBool,Dimension(3)           :: Has_displacementBC
      PetscReal,Dimension(3)           :: boundaryDisplacement
      PetscBool                        :: Has_damageBC
      PetscReal                        :: boundaryDamage
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
      enumerator  :: MEF90DefMech_ModeNULL = 0,    &
                     MEF90DefMech_ModeQuasiStatic, &
                     MEF90DefMech_ModeGradientFlow
   End Enum
   Character(len = MEF90_MXSTRLEN),Dimension(6),protected   :: MEF90DefMech_ModeList
   
   Enum,bind(c)
      enumerator  :: MEF90DefMech_defectLawElasticity = 0, &
                     MEF90DefMech_defectLawBrittleFracture, &
                     MEF90DefMech_defectLawPlasticity
   End Enum
   Character(len = MEF90_MXSTRLEN),Dimension(6),protected   :: MEF90DefMech_defectLawList
   
   Enum,bind(c)
      enumerator  :: MEF90DefMech_defectLawBrittleFractureAT1 = 0, &
                     MEF90DefMech_defectLawBrittleFractureAT2
   End Enum
   Character(len = MEF90_MXSTRLEN),Dimension(5),protected   :: MEF90DefMech_defectLawBrittleFractureList

   Enum,bind(c)
      enumerator  :: MEF90DefMech_defectLawPlasticityTresca = 0, &
                     MEF90DefMech_defectLawPlasticityVonMises
   End Enum
   Character(len = MEF90_MXSTRLEN),Dimension(5),protected   :: MEF90DefMech_defectLawPlasticityList
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

      MEF90DefMech_ModeList(1) = 'Null'
      MEF90DefMech_ModeList(2) = 'QuasiStatic'
      MEF90DefMech_ModeList(3) = 'GradientFlow'
      MEF90DefMech_ModeList(4) = 'MEF90_DefMech_Mode'
      MEF90DefMech_ModeList(5) = '_MEF90_DefMech_Mode'
      MEF90DefMech_ModeList(6) = ''
      
      MEF90DefMech_defectLawList(1) = 'Elasticity'
      MEF90DefMech_defectLawList(2) = 'BrittleFracture'
      MEF90DefMech_defectLawList(3) = 'Plasticity'
      MEF90DefMech_defectLawList(4) = 'MEF90DefMech_defectLaw'
      MEF90DefMech_defectLawList(5) = '_MEF90DefMech_defectLaw'
      MEF90DefMech_defectLawList(6) = ''

      MEF90DefMech_defectLawBrittleFractureList(1) = 'AT1'
      MEF90DefMech_defectLawBrittleFractureList(2) = 'AT2'
      MEF90DefMech_defectLawBrittleFractureList(3) = 'MEF90DefMech_defectLawBrittleFracture'
      MEF90DefMech_defectLawBrittleFractureList(4) = '_MEF90DefMech_defectLawBrittleFracture'
      MEF90DefMech_defectLawBrittleFractureList(5) = ''
      
      MEF90DefMech_defectLawPlasticityList(1) = 'Tresca'
      MEF90DefMech_defectLawPlasticityList(2) = 'VonMises'
      MEF90DefMech_defectLawPlasticityList(3) = 'MEF90DefMech_defectLawPlasticity'
      MEF90DefMech_defectLawPlasticityList(4) = '_MEF90DefMech_defectLawPlasticity'
      MEF90DefMech_defectLawPlasticityList(5) = ''
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
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr) 
      Call ISGetLocalSize(setIS,numSet,ierr);CHKERRQ(ierr)
      Allocate(DefMechCtx%CellSetOptionsBag(numSet),stat=ierr)
      Do set = 1, numSet
         Call PetscBagCreate(MEF90Ctx%comm,sizeofMEF90DefMechCellSetOptions,DefMechCtx%CellSetOptionsBag(set),ierr);CHKERRQ(ierr)
      End Do
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(Mesh,'Vertex Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr) 
      Call ISGetLocalSize(setIS,numSet,ierr);CHKERRQ(ierr)
      Allocate(DefMechCtx%VertexSetOptionsBag(numSet),stat=ierr)
      Do set = 1, numSet
         Call PetscBagCreate(MEF90Ctx%comm,sizeofMEF90DefMechVertexSetOptions,DefMechCtx%VertexSetOptionsBag(set),ierr);CHKERRQ(ierr)
      End Do
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      
      Nullify(DefMechCtx%force)
      Nullify(DefMechCtx%pressureforce)
      Nullify(DefMechCtx%boundaryDisplacement)
      Nullify(DefMechCtx%boundaryDamage)
      Nullify(DefMechCtx%Displacement)
      Nullify(DefMechCtx%Damage)
      Nullify(DefMechCtx%temperature)
      Nullify(DefMechCtx%plasticStrain)
   End Subroutine MEF90DefMechCtxCreate
   
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCtxSetSections"
!!!
!!!  
!!!  MEF90DefMechCtxSetSections: Set the data layout for each of the fields involved in a MEF90DefMechCtx_Type
!!!                               Uses Sieve convenience functions for now, but will be pulling layout informations
!!!                               from teh element types when switching to DMComplex
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechCtxSetSections(DefMechCtx,ierr)
      Type(MEF90DefMechCtx_Type),Intent(INOUT)        :: DefMechCtx
      PetscErrorCode,Intent(OUT)                      :: ierr
   
      Type(SectionReal)                               :: defaultSection
      PetscInt                                        :: dim

      Call DMMeshGetDimension(DefMechCtx%DM,dim,ierr);CHKERRQ(ierr)

      Call DMMeshGetVertexSectionReal(DefMechCtx%DMVect,"default",dim,defaultSection,ierr);CHKERRQ(ierr)
      Call DMMeshSetSectionReal(DefMechCtx%DMVect,"default",defaultSection,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(defaultSection,ierr);CHKERRQ(ierr)
      Call DMSetBlockSize(DefMechCtx%DMVect,dim,ierr);CHKERRQ(ierr)
   
      Call DMMeshGetCellSectionReal(DefMechCtx%cellDMVect,"default",dim,defaultSection,ierr);CHKERRQ(ierr)
      Call DMMeshSetSectionReal(DefMechCtx%cellDMVect,"default",defaultSection,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(defaultSection,ierr);CHKERRQ(ierr)
      Call DMSetBlockSize(DefMechCtx%CellDMVect,dim,ierr);CHKERRQ(ierr)
      
      Call DMMeshGetVertexSectionReal(DefMechCtx%DMScal,"default",1,defaultSection,ierr);CHKERRQ(ierr)
      Call DMMeshSetSectionReal(DefMechCtx%DMScal,"default",defaultSection,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(defaultSection,ierr);CHKERRQ(ierr)
      Call DMSetBlockSize(DefMechCtx%DMScal,1,ierr);CHKERRQ(ierr)

      Call DMMeshGetCellSectionReal(DefMechCtx%cellDMScal,"default",1,defaultSection,ierr);CHKERRQ(ierr)
      Call DMMeshSetSectionReal(DefMechCtx%cellDMScal,"default",defaultSection,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(defaultSection,ierr);CHKERRQ(ierr)
      Call DMSetBlockSize(DefMechCtx%cellDMScal,1,ierr);CHKERRQ(ierr)

      Call DMMeshGetCellSectionReal(DefMechCtx%cellDMMatS,"default",(dim*(dim+1))/2,defaultSection,ierr);CHKERRQ(ierr)
      Call DMMeshSetSectionReal(DefMechCtx%cellDMMatS,"default",defaultSection,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(defaultSection,ierr);CHKERRQ(ierr)   
      Call DMSetBlockSize(DefMechCtx%cellDMMatS,(dim*(dim+1))/2,ierr);CHKERRQ(ierr)
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

      Allocate(DefMechCtx%boundaryDisplacement,stat=ierr)
      Call DMCreateGlobalVector(DefMechCtx%DMVect,DefMechCtx%boundaryDisplacement,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(DefMechCtx%boundaryDisplacement,"boundary Displacement",ierr);CHKERRQ(ierr)
      Call VecSet(DefMechCtx%boundaryDisplacement,0.0_Kr,ierr);CHKERRQ(ierr)
   
      Allocate(DefMechCtx%force,stat=ierr)
      Call DMCreateGlobalVector(DefMechCtx%cellDMVect,DefMechCtx%force,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(DefMechCtx%force,"Force",ierr);CHKERRQ(ierr)
      Call VecSet(DefMechCtx%force,0.0_Kr,ierr);CHKERRQ(ierr)

      Allocate(DefMechCtx%pressureForce,stat=ierr)
      Call DMCreateGlobalVector(DefMechCtx%cellDMScal,DefMechCtx%pressureForce,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(DefMechCtx%pressureForce,"PressureForce",ierr);CHKERRQ(ierr)
      Call VecSet(DefMechCtx%pressureForce,0.0_Kr,ierr);CHKERRQ(ierr)

      Allocate(DefMechCtx%temperature,stat=ierr)
      Call DMCreateGlobalVector(DefMechCtx%DMScal,DefMechCtx%temperature,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(DefMechCtx%temperature,"temperature",ierr);CHKERRQ(ierr)
      Call VecSet(DefMechCtx%temperature,0.0_Kr,ierr);CHKERRQ(ierr)

      Allocate(DefMechCtx%plasticStrain,stat=ierr)
      Call DMCreateGlobalVector(DefMechCtx%CellDMMatS,DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(DefMechCtx%plasticStrain,"plasticStrain",ierr);CHKERRQ(ierr)
      Call VecSet(DefMechCtx%plasticStrain,0.0_Kr,ierr);CHKERRQ(ierr)
   
      Allocate(DefMechCtx%damage,stat=ierr)
      Call DMCreateGlobalVector(DefMechCtx%DMScal,DefMechCtx%damage,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(DefMechCtx%damage,"damage",ierr);CHKERRQ(ierr)
      Call VecSet(DefMechCtx%damage,0.0_Kr,ierr);CHKERRQ(ierr)

      Allocate(DefMechCtx%boundaryDamage,stat=ierr)
      Call DMCreateGlobalVector(DefMechCtx%DMScal,DefMechCtx%boundaryDamage,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(DefMechCtx%boundaryDamage,"boundaryDamage",ierr);CHKERRQ(ierr)
      Call VecSet(DefMechCtx%boundaryDamage,0.0_Kr,ierr);CHKERRQ(ierr)
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

      If (Associated(DefMechCtx%boundaryDisplacement)) Then 
         Call VecDestroy(DefMechCtx%boundaryDisplacement,ierr);CHKERRQ(ierr)
         DeAllocate(DefMechCtx%boundaryDisplacement)
         Nullify(DefMechCtx%boundaryDisplacement)
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
      
      Nullify(DefMechCtx%force)
      Nullify(DefMechCtx%pressureforce)
      Nullify(DefMechCtx%boundaryDisplacement)
      Nullify(DefMechCtx%boundaryDamage)
      Nullify(DefMechCtx%Displacement)
      Nullify(DefMechCtx%Damage)
      Nullify(DefMechCtx%temperature)
      Nullify(DefMechCtx%plasticStrain)
      Call DMDestroy(DefMechCtx%DMScal,ierr);CHKERRQ(ierr)
      Call DMDestroy(DefMechCtx%cellDMScal,ierr);CHKERRQ(ierr)
      Call DMDestroy(DefMechCtx%DMVect,ierr);CHKERRQ(ierr)
      Call DMDestroy(DefMechCtx%cellDMVect,ierr);CHKERRQ(ierr)
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
      Call PetscBagSetName(bag,trim(name),"DefMechGlobalOptions MEF90 Heat transfer global options",ierr);CHKERRQ(ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix),ierr);CHKERRQ(ierr)

      Call PetscBagRegisterEnum(bag,DefMechGlobalOptions%mode,MEF90DefMech_ModeList,default%mode,'DefMech_mode','Type of heat transfer computation',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterBool(bag,DefMechGlobalOptions%addDisplacementNullSpace,default%addDisplacementNullSpace,'disp_addNullSpace','Add null space to SNES',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,DefMechGlobalOptions%displacementOffset,default%displacementOffset,'displacement_Offset','Position of displacement field in EXO file',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,DefMechGlobalOptions%damageOffset,default%damageOffset,'damage_Offset','Position of damage field in EXO file',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,DefMechGlobalOptions%boundaryDamageOffset,default%damageOffset,'boundaryDamage_Offset','Position of boundary damage field in EXO file',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,DefMechGlobalOptions%stressOffset,default%stressOffset,'stress_Offset','Position of stress field in EXO file',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,DefMechGlobalOptions%temperatureOffset,default%temperatureOffset,'temperature_Offset','Position of temperature field in EXO file',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,DefMechGlobalOptions%plasticStrainOffset,default%plasticStrainOffset,'plasticStrain_Offset','Position of the plastic Strain field in EXO file',ierr);CHKERRQ(ierr)

      Call PetscBagRegisterEnum(bag,DefMechGlobalOptions%boundaryDisplacementScaling,MEF90ScalingList,default%boundaryDisplacementScaling,'boundaryDisplacement_scaling','Boundary displacement scaling',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,DefMechGlobalOptions%boundaryDisplacementOffset,default%boundaryDisplacementOffset,'boundaryDisplacement_Offset','Position of boundary displacement field in EXO file',ierr);CHKERRQ(ierr)

      Call PetscBagRegisterEnum(bag,DefMechGlobalOptions%forceScaling,MEF90ScalingList,default%forceScaling,'force_scaling','Force scaling',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,DefMechGlobalOptions%forceOffset,default%forceOffset,'force_Offset','Position of force field in EXO file',ierr);CHKERRQ(ierr)

      Call PetscBagRegisterEnum(bag,DefMechGlobalOptions%pressureForceScaling,MEF90ScalingList,default%pressureforceScaling,'pressureForce_scaling','Pressure force scaling',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,DefMechGlobalOptions%pressureForceOffset,default%pressureForceOffset,'pressureForce_Offset','Position of pressure force field in EXO file',ierr);CHKERRQ(ierr)
   End Subroutine PetscBagRegisterMEF90DefMechCtxGlobalOptions

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90DefMechCtxCellSetOptions"
!!!
!!!  
!!!  PetscBagRegisterMEF90DefMechCtxCellSetOptions:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine PetscBagRegisterMEF90DefMechCtxCellSetOptions(bag,name,prefix,default,ierr)
      PetscBag                                           :: bag
      Character(len=*),Intent(IN)                        :: prefix,name
      Type(MEF90DefMechCellSetOptions_Type),Intent(IN)   :: default
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(MEF90DefMechCellSetOptions_Type),pointer      :: DefMechCellSetOptions
      
      Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(bag,DefMechCellSetOptions,ierr);CHKERRQ(ierr)
      Call PetscBagSetName(bag,trim(name),"DefMechCellSetOptions MEF90 Heat transfer Cell Set options",ierr);CHKERRQ(ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix),ierr);CHKERRQ(ierr)
      
      DefMechCellSetOptions%force = default%Force
      DefMechCellSetOptions%Has_displacementBC   = default%Has_displacementBC
      DefMechCellSetOptions%boundaryDisplacement = default%boundaryDisplacement

      Call PetscBagRegisterInt(bag,DefMechCellSetOptions%ElemTypeShortIDDisplacement,default%ElemTypeShortIDDisplacement,'ShortIDDisplacement','Displacement element type ShortID',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt(bag,DefMechCellSetOptions%ElemTypeShortIDDamage,default%ElemTypeShortIDDamage,'ShortIDDamage','Damage field element type ShortID',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterRealArray(bag,DefMechCellSetOptions%force,3,'Force','[N.m^(-3) / N.m^(-2) / N.m^(-1)] (f): body / boundary force',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,DefMechCellSetOptions%pressureForce,default%pressureForce,'pressureForce','[N.m^(-2) / N.m^(-1)] (p): boundary pressureforce',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterEnum(bag,DefMechCellSetOptions%defectLaw,MEF90DefMech_defectLawList,default%defectLaw,'damageLaw','damage law',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterBoolArray(bag,DefMechCellSetOptions%Has_displacementBC,3,'DisplacementBC','Displacement has Dirichlet boundary Condition (Y/N)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterRealArray(bag,DefMechCellSetOptions%boundaryDisplacement,3,'boundaryDisplacement','[m] (U): Displacement boundary value',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterBool(bag,DefMechCellSetOptions%Has_DamageBC,default%Has_DamageBC,'DamageBC','Damage has Dirichlet boundary Condition (Y/N)',ierr);CHKERRQ(ierr)
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

      Type(MEF90DefMechVertexSetOptions_Type),pointer      :: DefMechVertexSetOptions
      Call PetscBagGetDataMEF90DefMechCtxVertexSetOptions(bag,DefMechVertexSetOptions,ierr);CHKERRQ(ierr)
      Call PetscBagSetName(bag,trim(name),"DefMechVertexSetOptions MEF90 Heat transfer Vertex Set options",ierr);CHKERRQ(ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix),ierr);CHKERRQ(ierr)

      DefMechVertexSetOptions%Has_displacementBC   = default%Has_displacementBC
      DefMechVertexSetOptions%boundaryDisplacement = default%boundaryDisplacement
      Call PetscBagRegisterBoolArray(bag,DefMechVertexSetOptions%Has_displacementBC,3,'DisplacementBC','Displacement has Dirichlet boundary Condition (Y/N)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterRealArray(bag,DefMechVertexSetOptions%boundaryDisplacement,3,'boundaryDisplacement','[m] (U): Displacement boundary value',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterBool(bag,DefMechVertexSetOptions%Has_DamageBC,default%Has_DamageBC,'DamageBC','Damage has Dirichlet boundary Condition (Y/N)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,DefMechVertexSetOptions%boundaryDamage,default%boundaryDamage,'[unit-less] (alpha): boundaryDamage','Damage boundary value',ierr);CHKERRQ(ierr)
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
      Character(len=MXSTLN),Dimension(:),POinter            :: cellSetNames
      Type(MEF90DefMechVertexSetOptions_Type),Intent(IN)    :: defaultVertexSetOptions
      Character(len=MXSTLN),Dimension(:),POinter            :: vertexSetNames
      PetscErrorCode,Intent(OUT)                            :: ierr
   
      Type(MEF90CtxGlobalOptions_Type),pointer              :: MEF90CtxGlobalOptions
      Type(MEF90DefMechCellSetOptions_Type)                 :: myDefaultCellSetOptions
      Type(MEF90Element_Type),Dimension(:),Pointer          :: ElemTypeScal,ElemTypeElast
      Type(IS)                                              :: setIS
      PetscInt,Dimension(:),Pointer                         :: setID
      PetscInt                                              :: set
      Character(len=MEF90_MXSTRLEN)                         :: IOBuffer,setName,setprefix

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
      Call MEF90_ISAllGatherMerge(DefMechCtx%MEF90Ctx%comm,setIS,ierr);CHKERRQ(ierr) 
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
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr) 
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

100 Format('Cell set ',I4)
101 Format('cs',I4.4,'_')
103 Format('Registering cell set ',I4,' prefix: ',A,'\n')
200 Format('Vertex set ',I4)
201 Format('vs',I4.4,'_')
203 Format('Registering vertex set ',I4,' prefix: ',A,'\n')
   End Subroutine MEF90DefMechCtxSetFromOptions
End Module m_MEF90_DefMechCtx
