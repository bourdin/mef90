#include "SimplePoisson.inc"

Module M_POISSON_TYPES
#include "finclude/petscdef.h"
   Use m_MEF90
   Implicit none
   !Private  
   !Public :: PoissonCtx_Type
   !Public :: PoissonCellSetProperties_Type
   !Public :: PoissonVertexSetProperties_Type
   
   Type PoissonCtx_Type
      Type(DM)                                        :: mesh
      Type(IS)                                        :: CellSetGlobalIS,VertexSetGlobalIS
      Type(IS)                                        :: CellSetLocalIS,VertexSetLocalIS
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer   :: Elem
      Type(Field)                                     :: BCU
      !Type(Field)                                     :: F !Make forces cell centered in assembly routines
      PetscInt                                        :: verbose
      PetscBag,Dimension(:),Pointer                   :: CellSetPropertiesBag
      PetscBag,Dimension(:),Pointer                   :: VertexSetPropertiesBag
      PetscBag,Dimension(:),Pointer                   :: MaterialPropertiesBag
   End Type PoissonCtx_Type

   Type PoissonCellSetProperties_Type   
      PetscInt                                        :: ElemTypeShortID
      !PetscInt                                        :: ElemTypeName
      PetscReal                                       :: Force
   End Type PoissonCellSetProperties_Type   
   
   Type PoissonVertexSetProperties_Type   
      PetscBool                                       :: Has_BC
      PetscReal                                       :: BC
   End Type PoissonVertexSetProperties_Type   

   Type(PoissonCellSetProperties_Type),Parameter      :: defaultCellSetProperties = PoissonCellSetProperties_Type( DEFAULT_ELEMENT_SHORTID,0.0_Kr )
   Type(PoissonVertexSetProperties_Type),Parameter    :: defaultVertexSetProperties = PoissonVertexSetProperties_Type( PETSC_TRUE,0.0_Kr )
End Module M_POISSON_TYPES

Module M_POISSONCELLSETPROPERTY_INTERFACE   
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use M_POISSON_TYPES
   Implicit None
   Private
   Public :: PetscBagGetDataPoissonCellSetProperties

   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use M_POISSON_TYPES
         PetscBag                                     :: bag
         type(PoissonCellSetProperties_Type),pointer  :: data
         PetscErrorCode                               :: ierr
      End subroutine PetscBagGetData
   End interface
Contains
   Subroutine PetscBagGetDataPoissonCellSetProperties(bag,data,ierr)
      PetscBag                                        :: bag
      type(PoissonCellSetProperties_Type),pointer     :: data
      PetscErrorCode                                  :: ierr
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataPoissonCellSetProperties
End Module M_POISSONCELLSETPROPERTY_INTERFACE   


Module M_POISSONVERTEXSETPROPERTY_INTERFACE   
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use M_POISSON_TYPES
   Implicit None
   Private
   Public :: PetscBagGetDataPoissonVertexSetProperties

   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use M_POISSON_TYPES
         PetscBag                                        :: bag
         type(PoissonVertexSetProperties_Type),pointer   :: data
         PetscErrorCode                                  :: ierr
      End subroutine PetscBagGetData
   End interface
Contains
   Subroutine PetscBagGetDataPoissonVertexSetProperties(bag,data,ierr)
      PetscBag                                           :: bag
      type(PoissonVertexSetProperties_Type),pointer      :: data
      PetscErrorCode                                     :: ierr
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataPoissonVertexSetProperties
End Module M_POISSONVERTEXSETPROPERTY_INTERFACE   
!!!
!!!   
   
Module M_POISSON
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use M_POISSON_TYPES
   Use M_POISSONCELLSETPROPERTY_INTERFACE
   Use M_POISSONVERTEXSETPROPERTY_INTERFACE
   Use m_MEF90
   Implicit NONE
   !Private
   
   !Public :: m_Poisson_Initialize
   !Public :: PetscBagRegisterPoissonCellSetProperties
   !Public :: PetscBagGetDataPoissonCellSetProperties
   !Public :: PetscBagRegisterPoissonVertexSetProperties
   !Public :: PetscBagGetDataPoissonVertexSetProperties
   
   PetscSizeT,protected    :: sizeofPoissonCellSetProperties
   PetscSizeT,protected    :: sizeofPoissonVertexSetProperties
   PetscSizeT,protected    :: sizeofPoissonCtx
   
Contains
#undef __FUNCT__
#define __FUNCT__ "m_Poisson_Initialize"
   Subroutine m_Poisson_Initialize(ierr)
      PetscErrorCode,intent(OUT)          :: ierr

      Type(PoissonCellSetProperties_Type),Target   :: CellSetProperties
      Type(PoissonVertexSetProperties_Type),Target :: VertexSetProperties
      Type(PoissonCtx_Type),Target                 :: PoissonCtx
      character(len=1),pointer                     :: dummychar(:)
      PetscSizeT                                   :: sizeofchar
   
      Call PetscDataTypeGetSize(PETSC_CHAR,sizeofchar,ierr)
      sizeofPoissonCellSetProperties = size(transfer(CellSetProperties,dummychar))*sizeofchar
      sizeofPoissonVertexSetProperties = size(transfer(VertexSetProperties,dummychar))*sizeofchar
      sizeofPoissonCtx = size(transfer(PoissonCtx,dummychar))*sizeofchar
   End Subroutine m_Poisson_Initialize
   
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
      !Call PetscBagRegisterEnum(bag,CellSetProperties%ElemTypeName,MEF90_knownElementNames,default%ElemTypeShortID,'Element type name',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt(bag,CellSetProperties%ElemTypeShortID,default%ElemTypeShortID,'Element type ShortID',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,CellSetProperties%Force,default%force,'Force','Magnitude of the applied force',ierr);CHKERRQ(ierr)
   End Subroutine PetscBagRegisterPoissonCellSetProperties

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterPoissonVertexSetProperties"
   Subroutine PetscBagRegisterPoissonVertexSetProperties(bag,name,prefix,default,ierr)
      PetscBag                                           :: bag
      Character(len=*),intent(IN)                        :: prefix,name
      type(PoissonVertexSetProperties_Type),intent(IN)   :: default
      PetscErrorCode,intent(OUT)                         :: ierr

      Type(PoissonVertexSetProperties_Type),pointer      :: VertexSetProperties

      Call PetscBagGetDataPoissonVertexSetProperties(bag,VertexSetProperties,ierr);CHKERRQ(ierr)
      Call PetscBagSetName(bag,trim(name),"VertexSetProperties object: Vertex Set properties",ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix), ierr)
      Call PetscBagRegisterBool(bag,VertexSetProperties%Has_BC,PETSC_FALSE,'Temp_HasBC','Temperature has Dirichlet boundary Condition (Y/N)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,VertexSetProperties%BC,0.0_Kr,'Temp_BC','Temperature boundary value',ierr);CHKERRQ(ierr)
   End Subroutine PetscBagRegisterPoissonVertexSetProperties
   


#undef __FUNCT__
#define __FUNCT__ "PoissonCtxInit"
   !!!
   !!!  
   !!!  PoissonCtxInit:
   !!!  
   !!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
   !!!
   Subroutine PoissonCtxInit(PoissonCtx,prefix,ierr)
      Type(PoissonCtx_Type),intent(OUT)         :: PoissonCtx
      Character(len=MEF90_MXSTRLEN),intent(IN)  :: prefix
      PetscErrorCode,Intent(OUT)                :: ierr
     
      Integer                                   :: cpu_ws,io_ws,exoidIN=0
      Real                                      :: exo_version
      Character(len=MEF90_MXSTRLEN)             :: filename,name   
      Type(DM)                                  :: tmp_mesh
      PetscInt                                  :: numCell
      PetscInt,Dimension(:),Pointer             :: setID
      PetscInt                                  :: set
      PetscInt                                  :: numCellSetGlobal,numVertexSetGlobal

      !!! 1. Read and distribute the DM from the exo file <prefix>.gen
      If (MEF90_MyRank == 0) Then
         cpu_ws = 8
         io_ws = 8
         filename = Trim(prefix)//'.gen'
         exoidIN = EXOPEN(filename,EXREAD,cpu_ws,io_ws,exo_version,ierr)
      End If
      
      !!! Read the mesh from <prefix.gen> and partition it
      If (MEF90_NumProcs == 1) Then
         Call DMmeshCreateExodusNG(PETSC_COMM_WORLD,exoidIN,PoissonCtx%mesh,ierr);CHKERRQ(ierr)
      Else
         Call DMmeshCreateExodusNG(PETSC_COMM_WORLD,exoidIN,Tmp_mesh,ierr);CHKERRQ(ierr)
      
         Call DMmeshDistribute(Tmp_mesh,PETSC_NULL_CHARACTER,PoissonCtx%mesh,ierr);CHKERRQ(ierr)
         Call DMDestroy(Tmp_mesh,ierr);CHKERRQ(ierr)
      End If
      
      !!!
      !!! Get local and compute global IS for cell and vertex sets
      !!!
      Call DMmeshGetLabelIdIS(PoissonCtx%mesh,'Cell Sets',PoissonCtx%CellSetLocalIS,ierr);CHKERRQ(ierr)
      Call ISDuplicate(PoissonCtx%CellSetLocalIS,PoissonCtx%CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,PoissonCtx%CellSetGlobalIS,ierr);CHKERRQ(ierr)
      
      Call DMmeshGetLabelIdIS(PoissonCtx%mesh,'Vertex Sets',PoissonCtx%VertexSetLocalIS,ierr);CHKERRQ(ierr)
      Call ISDuplicate(PoissonCtx%VertexSetLocalIS,PoissonCtx%VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,PoissonCtx%VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      
      !!! Get number of cell, vertices and dimension
      !Call DMmeshGetStratumSize(PoissonCtx%mesh,"depth",0,numVertex,ierr);CHKERRQ(ierr)
      Call DMmeshGetStratumSize(PoissonCtx%mesh,"height",0,numCell,ierr);CHKERRQ(ierr)
      !Call DMMeshGetDimension(PoissonCtx%mesh,numDim,ierr);CHKERRQ(ierr)

      !!! Allocate elements array and bags
      Allocate(PoissonCtx%Elem(numCell))

      !!! Allocate and register bags:
      Call ISGetIndicesF90(PoissonCtx%CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Allocate(PoissonCtx%CellSetPropertiesBag(size(setID)))
      Allocate(PoissonCtx%MaterialPropertiesBag(size(setID)))
      Do set = 1, size(setID)
         !!! get set name from EXO file
         !!! Write a function for this in m_EXO
         !!! Set a reasonable default value for the default element here
         Call PetscBagRegisterPoissonCellSetProperties(PoissonCtx%CellSetPropertiesBag(set),name,prefix,defaultCellSetProperties,ierr);CHKERRQ(ierr)
         Call PetscBagRegisterMatProp(PoissonCtx%MaterialPropertiesBag(set),name,prefix,DEFAULT_MATERIAL,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(PoissonCtx%CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      
      !!! Allocate and register bags:
      Call ISGetIndicesF90(PoissonCtx%VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Allocate(PoissonCtx%VertexSetPropertiesBag(size(setID)))
      Do set = 1, size(setID)
         Call PetscBagRegisterPoissonVertexSetProperties(PoissonCtx%VertexSetPropertiesBag(set),name,prefix,defaultVertexSetProperties,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(PoissonCtx%VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)

   End Subroutine PoissonCtxInit
End Module M_POISSON
