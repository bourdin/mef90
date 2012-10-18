#include "SimplePoisson.inc"
Module m_PoissonCellSetProperties_Type
#include "finclude/petscdef.h"
   Use m_MEF90
   Implicit none
   Private  
   Public :: PoissonCellSetProperties_Type

   Type PoissonCellSetProperties_Type   
      PetscInt                                        :: ElemTypeShortID
      PetscReal                                       :: Flux
      PetscReal                                       :: SurfaceThermalConductivity
      PetscReal                                       :: referenceTemp
   End Type PoissonCellSetProperties_Type   
End Module m_PoissonCellSetProperties_Type

Module m_PoissonCellSetProperties   
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_PoissonGlobalProperties
   Use m_PoissonCellSetProperties_Type
   Implicit None
   Private

   PetscSizeT,protected,public :: sizeofPoissonCellSetProperties

   Public :: PoissonCellSetProperties_Type
   Public :: PoissonCellSetPropertiesInitialize
   Public :: PetscBagGetDataPoissonCellSetProperties
   Public :: PetscBagRegisterPoissonCellSetProperties
   Public :: MEF90CtxPoissonCellSetPropertiesCreate

   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_PoissonCellSetProperties_Type
         PetscBag                                     :: bag
         Type(PoissonCellSetProperties_Type),pointer  :: data
         PetscErrorCode                               :: ierr
      End subroutine PetscBagGetData
   End interface
Contains

#undef __FUNCT__
#define __FUNCT__ "PoissonCellSetPropertiesInitialize"
!!!
!!!  PoissonCellSetPropertiesInitialize - setup all protected variables related to m_PoissonCellSetProperties:
!!!                                      namely: sizeofPoissonCellSetProperties
!!!
   Subroutine PoissonCellSetPropertiesInitialize(ierr)
      PetscErrorCode,intent(OUT)                   :: ierr

      Type(PoissonCellSetProperties_Type),Target   :: CellSetProperties
      character(len=1),pointer                     :: dummychar(:)
      PetscSizeT                                   :: sizeofchar
   
      Call PetscDataTypeGetSize(PETSC_CHAR,sizeofchar,ierr)
      sizeofPoissonCellSetProperties = size(transfer(CellSetProperties,dummychar))*sizeofchar
   End Subroutine PoissonCellSetPropertiesInitialize

#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataPoissonCellSetProperties"
!!!
!!!  PetscBagGetDataPoissonCellSetProperties - Custom interface to PetscGetData
!!!
   Subroutine PetscBagGetDataPoissonCellSetProperties(bag,data,ierr)
      PetscBag                                        :: bag
      Type(PoissonCellSetProperties_Type),pointer     :: data
      PetscErrorCode                                  :: ierr
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataPoissonCellSetProperties

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterPoissonCellSetProperties"
!!!
!!!  PetscBagRegisterPoissonCellSetProperties - Register a PetscBag associated with a PoissonCellSetProperties_Type
!!!                                            Typically, the bag is in a MEF90Ctx_Type
!!!
#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterPoissonCellSetProperties"
   Subroutine PetscBagRegisterPoissonCellSetProperties(bag,name,prefix,default,ierr)
      PetscBag                                        :: bag
      Character(len=*),intent(IN)                     :: prefix,name
      type(PoissonCellSetProperties_Type),intent(IN)  :: default
      PetscErrorCode,intent(OUT)                      :: ierr

      Type(PoissonCellSetProperties_Type),pointer     :: CellSetProperties

      Call PetscBagGetDataPoissonCellSetProperties(bag,CellSetProperties,ierr);CHKERRQ(ierr)
      Call PetscBagSetName(bag,trim(name),"CellSetProperties object: Cell Set properties",ierr);CHKERRQ(ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix),ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt(bag,CellSetProperties%ElemTypeShortID,default%ElemTypeShortID,'ShortID','Element type ShortID',ierr);CHKERRQ(ierr)
      !Call PetscBagRegisterEnum(bag,CellSetProperties%ElemTypeShortIDEnum,MEF90_knownElementNames,default%ElemTypeShortID,'shortidenum','Element type name',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,CellSetProperties%Flux,default%Flux,'Flux','[J.s^(-1).m^(-3) / J.s^(-1).m^(-2)] (f): Internal / applied heat flux',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,CellSetProperties%SurfaceThermalConductivity,default%SurfaceThermalConductivity,'SurfaceThermalConductivity','[J.s^(-1).m^(-2).K^(-1)] (H) Surface Thermal Conductivity',ierr)
      Call PetscBagRegisterReal(bag,CellSetProperties%ReferenceTemp,default%ReferenceTemp,'ReferenceTemp','Reference temperature T [K]',ierr)
      !Call PetscBagSetFromOptions(bag,ierr);CHKERRQ(ierr)
   End Subroutine PetscBagRegisterPoissonCellSetProperties

#undef __FUNCT__
#define __FUNCT__ "MEF90CtxPoissonCellSetPropertiesCreate"
!!!
!!!  
!!!  MEF90CtxPoissonCellSetPropertiesCreate - Call PetscBagRegisterPoissonCellSetProperties for the CellSet properties
!!!                                           of a MEF90Ctx_Type and PetscBagRegisterMEF90_MatProp for its material properties
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90CtxPoissonCellSetPropertiesCreate(MEF90Ctx,snesTemp,defaultCellSetProperties,defaultElemType,ierr)
      Type(MEF90Ctx_Type),intent(OUT)                    :: MEF90Ctx
      Type(SNES),Intent(IN)                              :: snesTemp
      Type(PoissonCellSetProperties_Type),intent(IN)     :: defaultCellSetProperties
      Type(MEF90Element_Type),Dimension(:),Pointer       :: defaultElemType
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(PoissonCellSetProperties_Type)                :: mydefaultCellSetProperties
      Type(DM)                                           :: mesh
      Type(PoissonGlobalProperties_Type),pointer         :: GlobalProperties
      Type(IS)                                           :: CellSetCellSetIS
      PetscInt,Dimension(:),Pointer                      :: setID
      PetscInt                                           :: set
      Character(len=MEF90_MXSTRLEN)                      :: IOBuffer,setName,setprefix

      Call SNESGetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)

      Call PetscBagGetDataPoissonGlobalProperties(MEF90Ctx%GlobalPropertiesBag,GlobalProperties,ierr);CHKERRQ(ierr)
      Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetCellSetIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetCellSetIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetCellSetIS,setID,ierr);CHKERRQ(ierr)
      Allocate(MEF90Ctx%CellSetPropertiesBag(size(setID)),stat=ierr)
      Allocate(MEF90Ctx%MaterialPropertiesBag(size(setID)),stat=ierr)
            
      Do set = 1, size(setID)
         Write(setName,100) setID(set)
         Write(setprefix,101) setID(set)
         If (GlobalProperties%verbose > 0) Then
            Write(IOBuffer,103) setID(set),trim(setprefix)
            Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
         End if
         
         mydefaultCellSetProperties = defaultCellSetProperties
         mydefaultCellSetProperties%ElemTypeShortID = defaultElemType(set)%ShortID
         
         Call PetscBagCreate(PETSC_COMM_WORLD,sizeofPoissonCellSetProperties,MEF90Ctx%CellSetPropertiesBag(set),ierr)
         Call PetscBagRegisterPoissonCellSetProperties(MEF90Ctx%CellSetPropertiesBag(set),setName,setprefix,mydefaultCellSetProperties,ierr);CHKERRQ(ierr)

         Call PetscBagCreate(PETSC_COMM_WORLD,SIZEOFMATPROP,MEF90Ctx%MaterialPropertiesBag(set),ierr)
         Call PetscBagRegisterMEF90_MatProp(MEF90Ctx%MaterialPropertiesBag(set),setName,setprefix,DEFAULT_MATERIAL,ierr);CHKERRQ(ierr)
         if (GlobalProperties%verbose > 0) Then
            Call PetscBagView(MEF90Ctx%CellSetPropertiesBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
            Call PetscBagView(MEF90Ctx%MaterialPropertiesBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
            Call PetscPrintf(PETSC_COMM_WORLD,'\n',ierr);CHKERRQ(ierr)
         End If
      End Do
      Call ISRestoreIndicesF90(CellSetCellSetIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetCellSetIS,ierr);CHKERRQ(ierr)
100 Format('Cell set ',I4)
101 Format('cs',I4.4,'_')
103 Format('Registering cell set ',I4,' prefix: ',A,'\n')
   End Subroutine MEF90CtxPoissonCellSetPropertiesCreate
End Module m_PoissonCellSetProperties   

