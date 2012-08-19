#include "SimplePoisson.inc"
Module m_PoissonVertexSetProperties_Type
#include "finclude/petscdef.h"
   Use m_MEF90
   Implicit none
   Private  
   Public :: PoissonVertexSetProperties_Type

   Type PoissonVertexSetProperties_Type   
      PetscBool                                       :: Has_BC
      PetscReal                                       :: BC
   End Type PoissonVertexSetProperties_Type   
End Module m_PoissonVertexSetProperties_Type

Module m_PoissonVertexSetProperties   
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_PoissonGlobalProperties
   Use m_PoissonVertexSetProperties_Type
   Implicit None
   Private

   PetscSizeT,protected,public :: sizeofPoissonVertexSetProperties

   Public :: PoissonVertexSetProperties_Type
   Public :: PoissonVertexSetPropertiesInitialize
   Public :: PetscBagGetDataPoissonVertexSetProperties
   Public :: PetscBagRegisterPoissonVertexSetProperties
   Public :: MEF90CtxPoissonVertexSetPropertiesCreate

   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_PoissonVertexSetProperties_Type
         PetscBag                                        :: bag
         Type(PoissonVertexSetProperties_Type),pointer   :: data
         PetscErrorCode                                  :: ierr
      End subroutine PetscBagGetData
   End interface
Contains

#undef __FUNCT__
#define __FUNCT__ "PoissonVertexSetPropertiesInitialize"
!!!
!!!  PoissonVertexSetPropertiesInitialize - setup all protected variables related to m_PoissonVertexSetProperties:
!!!                                         namely: sizeofPoissonVertexSetProperties
!!!
   Subroutine PoissonVertexSetPropertiesInitialize(ierr)
      PetscErrorCode,intent(OUT)                   :: ierr

      Type(PoissonVertexSetProperties_Type),Target :: VertexSetProperties
      character(len=1),pointer                     :: dummychar(:)
      PetscSizeT                                   :: sizeofchar
   
      Call PetscDataTypeGetSize(PETSC_CHAR,sizeofchar,ierr)
      sizeofPoissonVertexSetProperties = size(transfer(VertexSetProperties,dummychar))*sizeofchar
   End Subroutine PoissonVertexSetPropertiesInitialize

#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataPoissonVertexSetProperties"
!!!
!!!  PetscBagGetDataPoissonVertexSetProperties - Custom interface to PetscGetData
!!!
   Subroutine PetscBagGetDataPoissonVertexSetProperties(bag,data,ierr)
      PetscBag                                           :: bag
      type(PoissonVertexSetProperties_Type),pointer      :: data
      PetscErrorCode                                     :: ierr
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataPoissonVertexSetProperties

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterPoissonVertexSetProperties"
!!!
!!!  PetscBagRegisterPoissonVertexSetProperties - Register a PetscBag associated with a PoissonVertexSetProperties_Type
!!!                                               Typically, the bag is in a MEF90Ctx_Type
!!!
   Subroutine PetscBagRegisterPoissonVertexSetProperties(bag,name,prefix,default,ierr)
      PetscBag                                           :: bag
      Character(len=*),intent(IN)                        :: prefix,name
      type(PoissonVertexSetProperties_Type),intent(IN)   :: default
      PetscErrorCode,intent(OUT)                         :: ierr

      Type(PoissonVertexSetProperties_Type),pointer      :: VertexSetProperties

      Call PetscBagGetDataPoissonVertexSetProperties(bag,VertexSetProperties,ierr);CHKERRQ(ierr)
      Call PetscBagSetName(bag,trim(name),"VertexSetProperties object: Vertex Set properties",ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix), ierr)
      Call PetscBagRegisterBool(bag,VertexSetProperties%Has_BC,default%Has_BC,'TempBC','Temperature has Dirichlet boundary Condition (Y/N)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,VertexSetProperties%BC,default%BC,'Temp','Temperature boundary value',ierr);CHKERRQ(ierr)
      !Call PetscBagSetFromOptions(bag,ierr);CHKERRQ(ierr)
   End Subroutine PetscBagRegisterPoissonVertexSetProperties

#undef __FUNCT__
#define __FUNCT__ "MEF90CtxPoissonVertexSetPropertiesCreate"
!!!
!!!  
!!!  MEF90CtxPoissonVertexSetPropertiesCreate - Call PetscBagRegisterPoissonVertexSetProperties for the VertexSet properties
!!!                                             of a MEF90Ctx_Type
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90CtxPoissonVertexSetPropertiesCreate(MEF90Ctx,snesTemp,defaultVertexSetProperties,ierr)
      Type(MEF90Ctx_Type),intent(OUT)                    :: MEF90Ctx
      Type(SNES),Intent(IN)                              :: snesTemp
      Type(PoissonVertexSetProperties_Type),intent(IN)   :: defaultVertexSetProperties
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(DM)                                           :: mesh
      Type(PoissonGlobalProperties_Type),pointer         :: GlobalProperties
      Type(IS)                                           :: VertexSetGlobalIS
      PetscInt,Dimension(:),Pointer                      :: setID
      PetscInt                                           :: set
      Character(len=MEF90_MXSTRLEN)                      :: IOBuffer,setName,setprefix

      Call SNESGetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)

      Call PetscBagGetDataPoissonGlobalProperties(MEF90Ctx%GlobalPropertiesBag,GlobalProperties,ierr);CHKERRQ(ierr)
      
      Call DMmeshGetLabelIdIS(mesh,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Allocate(MEF90Ctx%VertexSetPropertiesBag(size(setID)),stat=ierr)

      Do set = 1, size(setID)
         Write(setName,200) setID(set)
         Write(setprefix,201) setID(set)
         If (GlobalProperties%verbose > 0) Then
            Write(IOBuffer,203) setID(set),trim(setprefix)
            Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
         End if
         Call PetscBagCreate(PETSC_COMM_WORLD,sizeofPoissonVertexSetProperties,MEF90Ctx%VertexSetPropertiesBag(set),ierr)
         Call PetscBagRegisterPoissonVertexSetProperties(MEF90Ctx%VertexSetPropertiesBag(set),setname,setprefix,defaultVertexSetProperties,ierr);CHKERRQ(ierr)
         If (GlobalProperties%verbose > 0) Then
            Call PetscBagView(MEF90Ctx%VertexSetPropertiesBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
            Call PetscPrintf(PETSC_COMM_WORLD,'\n',ierr);CHKERRQ(ierr)
         End If
      End Do
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)

200 Format('Vertex set ',I4)
201 Format('vs',I4.4,'_')
203 Format('Registering vertex set ',I4,' prefix: ',A,'\n')
   End Subroutine MEF90CtxPoissonVertexSetPropertiesCreate

End Module m_PoissonVertexSetProperties   
