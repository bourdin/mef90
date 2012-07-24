#include "SimplePoisson.inc"
Module m_PoissonGlobalProperties_Type
#include "finclude/petscdef.h"
   Use m_MEF90
   Implicit none
   Private  
   Public :: PoissonGlobalProperties_Type

   Type PoissonGlobalProperties_Type   
      PetscInt                                        :: verbose
      PetscBool                                       :: splitIO ! make enum?
   End Type PoissonGlobalProperties_Type   
End Module m_PoissonGlobalProperties_Type

Module m_PoissonGlobalProperties   
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_PoissonGlobalProperties_Type
   Implicit None
   Private

   PetscSizeT,protected,public :: sizeofPoissonGlobalProperties

   Public :: PoissonGlobalProperties_Type
   Public :: PoissonGlobalPropertiesInitialize
   Public :: PetscBagGetDataPoissonGlobalProperties
   Public :: PetscBagRegisterPoissonGlobalProperties
   Public :: MEF90CtxPoissonGlobalPropertiesCreate

   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_PoissonGlobalProperties_Type
         PetscBag                                     :: bag
         Type(PoissonGlobalProperties_Type),pointer   :: data
         PetscErrorCode                               :: ierr
      End subroutine PetscBagGetData
   End interface
Contains

#undef __FUNCT__
#define __FUNCT__ "PoissonGlobalPropertiesInitialize"
!!!
!!!  PoissonGlobalPropertiesInitialize - setup all protected variables related to m_PoissonGlobalProperties:
!!!                                      namely: sizeofPoissonGlobalProperties
!!!
   Subroutine PoissonGlobalPropertiesInitialize(ierr)
      PetscErrorCode,intent(OUT)                   :: ierr

      Type(PoissonGlobalProperties_Type),Target    :: GlobalProperties
      character(len=1),pointer                     :: dummychar(:)
      PetscSizeT                                   :: sizeofchar
   
      Call PetscDataTypeGetSize(PETSC_CHAR,sizeofchar,ierr)
      sizeofPoissonGlobalProperties = size(transfer(GlobalProperties,dummychar))*sizeofchar
   End Subroutine PoissonGlobalPropertiesInitialize

#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataPoissonGlobalProperties"
!!!
!!!  PetscBagGetDataPoissonGlobalProperties - Custom interface to PetscGetData
!!!
   Subroutine PetscBagGetDataPoissonGlobalProperties(bag,data,ierr)
      PetscBag                                        :: bag
      Type(PoissonGlobalProperties_Type),pointer      :: data
      PetscErrorCode                                  :: ierr
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataPoissonGlobalProperties

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterPoissonGlobalProperties"
!!!
!!!  PetscBagRegisterPoissonGlobalProperties - Register a PetscBag associated with a PoissonGlobalProperties_Type
!!!                                            Typically, the bag is in a MEF90Ctx_Type
!!!
   Subroutine PetscBagRegisterPoissonGlobalProperties(bag,name,prefix,default,ierr)
      PetscBag                                        :: bag
      Character(len=*),intent(IN)                     :: prefix,name
      type(PoissonGlobalProperties_Type),intent(IN)   :: default
      PetscErrorCode,intent(OUT)                      :: ierr

      Type(PoissonGlobalProperties_Type),pointer      :: GlobalProperties
      PetscBool                                       :: flg

      Call PetscBagGetDataPoissonGlobalProperties(bag,GlobalProperties,ierr);CHKERRQ(ierr)
      Call PetscBagSetName(bag,trim(name),"Global properties object",ierr);CHKERRQ(ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix),ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt(bag,GlobalProperties%verbose,default%verbose,'verbose','verbosity level',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterBool(bag,GlobalProperties%splitIO,default%splitio,'splitio','set to true for distributed IO split in multiple exo files',ierr);CHKERRQ(ierr)
   End Subroutine PetscBagRegisterPoissonGlobalProperties

#undef __FUNCT__
#define __FUNCT__ "MEF90CtxPoissonGlobalPropertiesCreate"
!!!
!!!  
!!!  MEF90CtxPoissonGlobalPropertiesCreate - Call PetscBagRegisterPoissonGlobalProperties for the global properties
!!!                                          of a MEF90Ctx_Type. Mostly for symmetry with CellSet and VertexSet functions
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90CtxPoissonGlobalPropertiesCreate(MEF90Ctx,defaultGlobalProperties,ierr)
      Type(MEF90Ctx_Type),intent(OUT)                 :: MEF90Ctx
      Type(PoissonGlobalProperties_Type),Intent(IN)   :: defaultGlobalProperties
      PetscErrorCode,Intent(OUT)                      :: ierr
   
      Call PetscBagCreate(PETSC_COMM_WORLD,sizeofPoissonGlobalProperties,MEF90Ctx%GlobalPropertiesBag,ierr)
      Call PetscBagRegisterPoissonGlobalProperties(MEF90Ctx%GlobalPropertiesBag,"SimplePoissonNG",PETSC_NULL_CHARACTER,defaultGlobalProperties,ierr);CHKERRQ(ierr)
   End Subroutine MEF90CtxPoissonGlobalPropertiesCreate
End Module m_PoissonGlobalProperties