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
      PetscBag,Dimension(:),Pointer                   :: CellSetPropertiesBag
      PetscBag,Dimension(:),Pointer                   :: VertexSetPropertiesBag
      PetscBag,Dimension(:),Pointer                   :: MaterialPropertiesBag
   End Type PoissonCtx_Type

   Type PoissonCellSetProperties_Type   
      PetscInt                                        :: ElemTypeShortID
      !!! Maybe this plus handling of short name + gauss quadrature order 
      !!! should be part of a ElementBag and handled in m_MEF_Elements?
      !PetscInt                                        :: ElemTypeShortIDEnum
      PetscReal                                       :: Force
   End Type PoissonCellSetProperties_Type   
   
   Type PoissonVertexSetProperties_Type   
      PetscBool                                       :: Has_BC
      PetscReal                                       :: BC
   End Type PoissonVertexSetProperties_Type   

   !Type(PoissonCellSetProperties_Type),Parameter      :: defaultCellSetProperties = PoissonCellSetProperties_Type( DEFAULT_ELEMENT_SHORTID,0.0_Kr )
   !Type(PoissonVertexSetProperties_Type),Parameter    :: defaultVertexSetProperties = PoissonVertexSetProperties_Type( PETSC_TRUE,0.0_Kr )
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
      Type(PoissonCellSetProperties_Type),pointer     :: data
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
      Call PetscBagRegisterInt(bag,CellSetProperties%ElemTypeShortID,default%ElemTypeShortID,'ShortID','Element type ShortID',ierr);CHKERRQ(ierr)
      !Call PetscBagRegisterEnum(bag,CellSetProperties%ElemTypeShortIDEnum,MEF90_knownElementNames,default%ElemTypeShortIDEnum,'shortidenum','Element type name',ierr);CHKERRQ(ierr)
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
      Call PetscBagRegisterBool(bag,VertexSetProperties%Has_BC,PETSC_FALSE,'TempBC','Temperature has Dirichlet boundary Condition (Y/N)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,VertexSetProperties%BC,0.0_Kr,'Temp','Temperature boundary value',ierr);CHKERRQ(ierr)
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
      Type(PoissonCtx_Type),intent(OUT)            :: PoissonCtx
      Character(len=*),intent(IN)                  :: prefix
      PetscErrorCode,Intent(OUT)                   :: ierr
     
      PetscInt                                     :: verbose
      PetscBool                                    :: flg
      Integer                                      :: cpu_ws,io_ws,exoidIN=0
      Real                                         :: exo_version
      Character(len=MXSTLN)                        :: filename,EXOElemType
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer,setName,setprefix
      Type(DM)                                     :: tmp_mesh
      PetscInt                                     :: numCell
      PetscInt,Dimension(:),Pointer                :: setID
      PetscInt                                     :: set
      PetscInt                                     :: numCellSetGlobal,numVertexSetGlobal
      Type(Element_Type)                           :: defaultElemType
      Type(PoissonCellSetProperties_Type)          :: defaultCellSetProperties
      Type(PoissonCellSetProperties_Type),pointer  :: CellSetProperties
      Type(PoissonVertexSetProperties_Type)        :: defaultVertexSetProperties
      integer                                      :: exoerr
      PetscInt,Parameter                           :: QuadratureOrder = 2
      Type(IS)                                     :: setIS
      PetscInt,Dimension(:),Pointer                :: cellID
      
      verbose = 0
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'--verbose',verbose,flg,ierr);CHKERRQ(ierr)

      !!! 1. Read and distribute the DM from the exo file <prefix>.gen
      If (MEF90_MyRank == 0) Then
         cpu_ws = 8
         io_ws = 8
         filename = Trim(prefix)//'.gen'
         exoidIN = EXOPEN(filename,EXREAD,cpu_ws,io_ws,exo_version,exoerr)
         If (verbose >0) Then
            Write(IOBuffer,99) exoidIN, exoerr
            Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr)
         End If
99 Format('EXOPEN status: ',I4,' exoidIN: ',I4,'\n')         
         If (exoerr < 0) Then
            Write(IOBuffer,*) '\n\nError opening EXO file ', trim(filename),'\n\n'
            Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr);
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,IOBuffer,ierr);
         EndIf
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
      
      
      !!! Allocate and register cell sets bags: CellSetProperties and MatProp
      !!! A loop on CellSetGlobalIS is required since PetscBarRegister is collective on PETSC_COMM_WORLD
      !!!
      Call ISGetIndicesF90(PoissonCtx%CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Allocate(PoissonCtx%CellSetPropertiesBag(size(setID)))
      Allocate(PoissonCtx%MaterialPropertiesBag(size(setID)))
      
      Do set = 1, size(setID)
         Write(setName,100) setID(set)
         Write(setprefix,101) setID(set)
         If (verbose > 0) Then
            Write(IOBuffer,103) setID(set),trim(setprefix)
            Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
         End if
         !!! Use a reasonable guess from the EXO file if nothing is given
         If (MEF90_MyRank == 0) Then
            Call EXOSetElementType_Load(exoidIN,setID(set),EXOElemType)
            Call EXO2Element_TypeScal(EXOElemType,MEF90_DIM,defaultElemType)
         End If
         Call MPI_BCast(defaultElemType%ShortID,MXSTLN,MPI_CHARACTER,0,PETSC_COMM_WORLD,ierr)

         defaultCellSetProperties = PoissonCellSetProperties_Type(defaultElemType%ShortID,0.0_Kr)
         Call PetscBagCreate(PETSC_COMM_WORLD,sizeofPoissonCellSetProperties,PoissonCtx%CellSetPropertiesBag(set),ierr)
         Call PetscBagRegisterPoissonCellSetProperties(PoissonCtx%CellSetPropertiesBag(set),setName,setprefix,defaultCellSetProperties,ierr);CHKERRQ(ierr)

         Call PetscBagCreate(PETSC_COMM_WORLD,sizeofMEF90_MatProp2D,PoissonCtx%MaterialPropertiesBag(set),ierr)
         Call PetscBagRegisterMEF90_MatProp(PoissonCtx%MaterialPropertiesBag(set),setName,setprefix,DEFAULT_MATERIAL,ierr);CHKERRQ(ierr)
         if (verbose > 0) Then
            Call PetscBagView(PoissonCtx%CellSetPropertiesBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
            Call PetscPrintf(PETSC_COMM_WORLD,'\n',ierr);CHKERRQ(ierr)
            Call PetscBagView(PoissonCtx%MaterialPropertiesBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
            Call PetscPrintf(PETSC_COMM_WORLD,'\n',ierr);CHKERRQ(ierr)
         End If
      End Do
      Call ISRestoreIndicesF90(PoissonCtx%CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

      !!! Allocate and register vertex properties bags:
      Call ISGetIndicesF90(PoissonCtx%VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Allocate(PoissonCtx%VertexSetPropertiesBag(size(setID)))
      defaultVertexSetProperties = PoissonVertexSetProperties_Type( PETSC_TRUE,0.0_Kr )
      Do set = 1, size(setID)
         Write(setName,200) setID(set)
         Write(setprefix,201) setID(set)
         If (verbose > 0) Then
            Write(IOBuffer,203) setID(set),trim(setprefix)
            Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
         End if
         Call PetscBagCreate(PETSC_COMM_WORLD,sizeofPoissonVertexSetProperties,PoissonCtx%VertexSetPropertiesBag(set),ierr)
         Call PetscBagRegisterPoissonVertexSetProperties(PoissonCtx%VertexSetPropertiesBag(set),setname,setprefix,defaultVertexSetProperties,ierr);CHKERRQ(ierr)
         if (verbose > 0) Then
            Call PetscBagView(PoissonCtx%VertexSetPropertiesBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
         End If
      End Do
      Call ISRestoreIndicesF90(PoissonCtx%VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
100 Format('Cell set ',I4)
101 Format('cs',I4.4,'_')
103 Format('=== Registering cell set ',I4,' prefix: ',A,'\n')
200 Format('Vertex set ',I4)
201 Format('vs',I4.4,'_')
203 Format('=== Registering vertex set ',I4,' prefix: ',A,'\n')

      !!! Geometry and properties is done
      !!! Now passing to FE part

      !!! Allocate elements array then initialize the elements
      Allocate(PoissonCtx%Elem(numCell))
      Call ISGetIndicesF90(PoissonCtx%CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      !!! This again needs to use the GLOBAL IS since the elemType is a global property
      Do set = 1, size(setID)
         Call PetscBagGetDataPoissonCellSetProperties(PoissonCtx%CellSetPropertiesBag(set),CellSetProperties,ierr)

         Call DMmeshGetStratumIS(PoissonCtx%mesh,'Cell Sets',setID(set),setIS,ierr); CHKERRQ(ierr)
         Call ElementInit(PoissonCtx%mesh,setIS,PoissonCtx%Elem,QuadratureOrder,CellSetProperties%ElemTypeShortID)
         !!! I could use a different quadrature order on different blocks if I'd add it to the CellSetProperties!

         !!! I need a function that allocates a SectionReal using the mesh and the dof layout from the elemType

         !!! Get IS for block setID(set)
         !Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
   
         !Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do


      !!! 
      !!! Create a Section consistent with the element choice
      !!!

   End Subroutine PoissonCtxInit
End Module M_POISSON
