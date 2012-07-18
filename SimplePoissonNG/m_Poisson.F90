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
      Type(SectionReal)                               :: Section
      Type(VecScatter)                                :: ScatterSecToVec
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
         Type(PoissonCellSetProperties_Type),pointer  :: data
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
   
   !external SimplePoissonNGBilinearFormAssembly
   
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
      !Call PetscBagSetFromOptions(bag,ierr);CHKERRQ(ierr)
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
      !Call PetscBagSetFromOptions(bag,ierr);CHKERRQ(ierr)
   End Subroutine PetscBagRegisterPoissonVertexSetProperties
   


#undef __FUNCT__
#define __FUNCT__ "PoissonCtxCreate"
   !!!
   !!!  
   !!!  PoissonCtxInit:
   !!!  
   !!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
   !!!
   Subroutine PoissonCtxCreate(PoissonCtx,exoid,ierr)
      Type(PoissonCtx_Type),intent(OUT)            :: PoissonCtx
      Integer,Intent(IN)                           :: exoid
      PetscErrorCode,Intent(OUT)                   :: ierr
     
      PetscInt                                     :: verbose
      PetscBool                                    :: flg
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
      PetscInt,Parameter                           :: QuadratureOrder = 4
      Type(IS)                                     :: setIS
      PetscInt,Dimension(:),Pointer                :: cellID
      
      verbose = 0
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'--verbose',verbose,flg,ierr);CHKERRQ(ierr)

      !!! Read the mesh from <prefix.gen> and partition it
      If (MEF90_NumProcs == 1) Then
         Call DMmeshCreateExodusNG(PETSC_COMM_WORLD,exoid,PoissonCtx%mesh,ierr);CHKERRQ(ierr)
      Else
         Call DMmeshCreateExodusNG(PETSC_COMM_WORLD,exoid,Tmp_mesh,ierr);CHKERRQ(ierr)
      
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
            Call EXOSetElementType_Load(exoid,setID(set),EXOElemType)
            Call EXO2Element_TypeScal(EXOElemType,MEF90_DIM,defaultElemType)
         End If
         Call MPI_BCast(defaultElemType%ShortID,MXSTLN,MPI_CHARACTER,0,PETSC_COMM_WORLD,ierr)

         defaultCellSetProperties = PoissonCellSetProperties_Type(defaultElemType%ShortID,0.0_Kr)
         Call PetscBagCreate(PETSC_COMM_WORLD,sizeofPoissonCellSetProperties,PoissonCtx%CellSetPropertiesBag(set),ierr)
         Call PetscBagRegisterPoissonCellSetProperties(PoissonCtx%CellSetPropertiesBag(set),setName,setprefix,defaultCellSetProperties,ierr);CHKERRQ(ierr)

         Call PetscBagCreate(PETSC_COMM_WORLD,SIZEOFMATPROP,PoissonCtx%MaterialPropertiesBag(set),ierr)
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

         Call DMmeshGetStratumIS(PoissonCtx%mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(ierr)
         Call ElementInit(PoissonCtx%mesh,setIS,PoissonCtx%Elem,QuadratureOrder,CellSetProperties%ElemTypeShortID)
         
         !!! I could use a different quadrature order on different blocks if I'd add it to the CellSetProperties!

         !!! I need a function that allocates a SectionReal using the mesh and the dof layout from the elemType
         !!! But with old style sieve, it is a bit of a pain in the ass, so I'll use th shortcuts for now.
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do


      !!! 
      !!! Create a Section consistent with the element choice
      !!!
      Call DMMeshGetVertexSectionReal(PoissonCtx%mesh,"default",1,PoissonCtx%Section,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(PoissonCtx%mesh,PoissonCtx%Section,PoissonCtx%ScatterSecToVec,ierr);CHKERRQ(ierr)
      !Call SectionRealView(PoissonCtx%Section,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr);
   End Subroutine PoissonCtxCreate
   
#undef __FUNCT__
#define __FUNCT__ "PoissonCtxDestroy"
!!!
!!!  
!!!  PoissonCtxDestroy:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine PoissonCtxDestroy(PoissonCtx,ierr)
   Type(PoissonCtx_Type)                                    :: PoissonCtx
   PetscErrorCode,Intent(OUT)                               :: ierr
   
   PetscInt                                                 :: e,set,nset
   
   Call ISGetLocalSize(PoissonCtx%CellSetGlobalIS,nset,ierr);CHKERRQ(ierr)
   Do set = 1, nset
      Call PetscBagDestroy(PoissonCtx%CellSetPropertiesBag(set),ierr);CHKERRQ(ierr)
      Call PetscBagDestroy(PoissonCtx%MaterialPropertiesBag(set),ierr);CHKERRQ(ierr)
   End Do
   DeAllocate(PoissonCtx%CellSetPropertiesBag)
   DeAllocate(PoissonCtx%MaterialPropertiesBag)
   Call ISDestroy(PoissonCtx%CellSetGlobalIS,ierr);CHKERRQ(ierr)
   Call ISDestroy(PoissonCtx%CellSetLocalIS,ierr);CHKERRQ(ierr)

   Call ISGetLocalSize(PoissonCtx%VertexSetGlobalIS,nset,ierr);CHKERRQ(ierr)
   Do set = 1, nset
      Call PetscBagDestroy(PoissonCtx%VertexSetPropertiesBag(set),ierr);CHKERRQ(ierr)
   End Do
   DeAllocate(PoissonCtx%VertexSetPropertiesBag)
   Call ISDestroy(PoissonCtx%VertexSetGlobalIS,ierr);CHKERRQ(ierr)
   
   Call SectionRealDestroy(PoissonCtx%Section,ierr);CHKERRQ(ierr)
   Call VecScatterDestroy(PoissonCtx%ScatterSecToVec,ierr);CHKERRQ(ierr)
   Do e = 1, size(PoissonCtx%Elem)
      Call ElementDestroy(PoissonCtx%Elem(e))
   End Do
   DeAllocate(PoissonCtx%Elem)
   
   Call DMDestroy(PoissonCtx%mesh,ierr);CHKERRQ(ierr);
End Subroutine PoissonCtxDestroy 

#undef __FUNCT__
#define __FUNCT__ "SimplePoissonNGBilinearFormAssembly"
!!!
!!!  
!!!  SimplePoissonNGBilinearFormAssembly:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine SimplePoissonNGBilinearFormAssembly(snesTemp,x,A,M,flg,PoissonCtx,ierr)
   Type(SNES),Intent(IN)                           :: snesTemp
   Type(Vec),Intent(IN)                            :: x
   Type(Mat),Intent(INOUT)                         :: A,M
   MatStructure,Intent(INOUT)                      :: flg
   Type(PoissonCtx_Type),intent(IN)                :: PoissonCtx
   PetscErrorCode,Intent(OUT)                      :: ierr  
   
   Type(IS)                                        :: setIS,setISdof
   PetscInt,dimension(:),Pointer                   :: setID
   PetscInt                                        :: set
   Type(MEF90_MATPROP),pointer                     :: matpropSet
   Type(PoissonCellSetProperties_Type),pointer     :: cellSetProperties
   Type(PoissonVertexSetProperties_Type),pointer   :: vertexSetProperties
   Type(Element_Type)                              :: elemType
   
   Call ISGetIndicesF90(PoissonCtx%CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      Call DMMeshGetStratumIS(PoissonCtx%mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)

      Call PetscBagGetDataMEF90_MatProp(PoissonCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataPoissonCellSetProperties(PoissonCtx%CellSetPropertiesBag(set),cellSetProperties,ierr);CHKERRQ(ierr)

      Call Element_TypeFindByID(cellSetProperties%ElemTypeShortID,elemType) 

      If (elemType%coDim == 0) Then
         Call MEF90_DiffusionBilinearFormSet(A,PoissonCtx%mesh,PoissonCtx%Section,setIS,matpropSet%ThermalDiffusivity,0.0_Kr,PoissonCtx%elem,elemType,ierr);CHKERRQ(ierr)
      End If
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   End Do     
   Call ISRestoreIndicesF90(PoissonCtx%CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
   Call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)

   !!!
   !!! Should this really be called here?
   !!!
   Call ISGetIndicesF90(PoissonCtx%VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      Call PetscBagGetDataPoissonVertexSetProperties(PoissonCtx%VertexSetPropertiesBag(set),vertexSetProperties,ierr);CHKERRQ(ierr)
      If (vertexSetProperties%Has_BC) Then
         Call DMMeshGetStratumIS(PoissonCtx%mesh,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call DMMeshISCreateISglobaldof(PoissonCtx%mesh,PoissonCtx%Section,setIS,0,setISdof,ierr);CHKERRQ(ierr)
         Call MatZeroRowsColumnsIS(A,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
      End If
   End Do
   Call ISRestoreIndicesF90(PoissonCtx%VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   
   flg = SAME_NONZERO_PATTERN
End Subroutine SimplePoissonNGBilinearFormAssembly

#undef __FUNCT__
#define __FUNCT__ "SimplePoissonNGOperatorAssembly"
!!!
!!!  
!!!  SimplePoissonNGOperatorAssembly:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine SimplePoissonNGOperatorAssembly(snesTemp,x,residual,PoissonCtx,ierr)
!!! Check that BC handling is correct (I am pretty sure it is not, BC should not be needed)
   Type(SNES),Intent(IN)                           :: snesTemp
   Type(Vec),Intent(IN)                            :: x
   Type(Vec),Intent(INOUT)                         :: residual
   Type(PoissonCtx_Type),intent(IN)                :: PoissonCtx
   PetscErrorCode,Intent(OUT)                      :: ierr

   Type(SectionReal)                               :: residualSec
   Type(IS)                                        :: setIS,setISdof
   PetscInt,dimension(:),Pointer                   :: setID
   PetscInt,Dimension(:),Pointer                   :: setIdx
   PetscInt                                        :: set
   Type(MEF90_MATPROP),pointer                     :: matpropSet
   Type(PoissonCellSetProperties_Type),pointer     :: cellSetProperties
   Type(PoissonVertexSetProperties_Type),pointer   :: vertexSetProperties
   Type(Element_Type)                              :: elemType
   PetscReal,Dimension(:),Pointer                  :: BC
   PetscReal                                       :: F
   
   Call SectionRealDuplicate(PoissonCtx%Section,residualSec,ierr);CHKERRQ(ierr)
   Call SectionRealSet(ResidualSec,0.0_Kr,ierr);CHKERRQ(ierr)
   Call VecSet(residual,0.0_Kr,ierr);CHKERRQ(ierr)
   ! probably not needed
   Call SectionRealToVec(PoissonCtx%Section,PoissonCtx%ScatterSecToVec,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr)

   Call ISGetIndicesF90(PoissonCtx%CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      Call DMMeshGetStratumIS(PoissonCtx%mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)

      Call PetscBagGetDataMEF90_MatProp(PoissonCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataPoissonCellSetProperties(PoissonCtx%CellSetPropertiesBag(set),cellSetProperties,ierr);CHKERRQ(ierr)
      Call Element_TypeFindByID(cellSetProperties%ElemTypeShortID,elemType) 

      !!! Assembly part of the residual coming from the bilinear form on the blocks of 
      !!! codimension 0
      If (elemType%coDim == 0) Then
         Call MEF90_DiffusionOperatorSet(residualSec,PoissonCtx%mesh,PoissonCtx%Section,setIS,matpropSet%ThermalDiffusivity,0.0_Kr,PoissonCtx%elem,elemType,ierr);CHKERRQ(ierr)
      End If
      !!! Assembly the force part
      If (cellSetProperties%Force /= 0.0_Kr) Then
         F = -cellSetProperties%Force
         Call MEF90_DiffusionRHSSet(residualSec,PoissonCtx%mesh,F,setIS,PoissonCtx%elem,elemType,ierr);CHKERRQ(ierr)
      End If
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   End Do     
   Call ISRestoreIndicesF90(PoissonCtx%CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   !!! "Ghost update" for the residual Section
   Call SectionRealComplete(residualSec,ierr);CHKERRQ(ierr)
   !!! Scatter back from SectionReal to Vec
   Call SectionRealToVec(residualSec,PoissonCtx%ScatterSecToVec,SCATTER_FORWARD,residual,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(residualSec,ierr);CHKERRQ(ierr)
   
   
   !!!
   !!! Zero-out BC entries in the residual
   !!!
   Call ISGetIndicesF90(PoissonCtx%VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      Call PetscBagGetDataPoissonVertexSetProperties(PoissonCtx%VertexSetPropertiesBag(set),vertexSetProperties,ierr);CHKERRQ(ierr)
      If (vertexSetProperties%Has_BC) Then
         Call DMMeshGetStratumIS(PoissonCtx%mesh,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call DMMeshISCreateISglobaldof(PoissonCtx%mesh,PoissonCtx%Section,setIS,0,setISdof,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
         Allocate(BC(size(setIdx)))
         BC = 0.0_Kr
         Call VecSetValues(residual,size(setIdx),setIdx,BC,INSERT_VALUES,ierr);CHKERRQ(ierr)
         DeAllocate(BC)
         Call ISRestoreIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
      End If
   End Do
   Call ISRestoreIndicesF90(PoissonCtx%VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Call VecAssemblyBegin(residual,ierr)
   Call VecAssemblyEnd(residual,ierr)
End Subroutine SimplePoissonNGOperatorAssembly


#undef __FUNCT__
#define __FUNCT__ "SimplePoissonFormInitialGuess"
!!!
!!!  
!!!  SimplePoissonFormInitialGuess:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine SimplePoissonFormInitialGuess(x,PoissonCtx,ierr)
   Type(Vec),Intent(IN)                            :: x
   Type(PoissonCtx_Type),intent(IN)                :: PoissonCtx
   PetscErrorCode,Intent(OUT)                      :: ierr

   Type(PoissonVertexSetProperties_Type),pointer   :: vertexSetProperties
   Type(IS)                                        :: setIS,setISdof
   PetscInt,dimension(:),Pointer                   :: setID
   PetscInt,Dimension(:),Pointer                   :: setIdx
   PetscInt                                        :: set
   PetscReal,Dimension(:),Pointer                  :: BC

   Call ISGetIndicesF90(PoissonCtx%VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      Call PetscBagGetDataPoissonVertexSetProperties(PoissonCtx%VertexSetPropertiesBag(set),vertexSetProperties,ierr);CHKERRQ(ierr)
      If (vertexSetProperties%Has_BC) Then
         Call DMMeshGetStratumIS(PoissonCtx%mesh,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call DMMeshISCreateISglobaldof(PoissonCtx%mesh,PoissonCtx%Section,setIS,0,setISdof,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
         Allocate(BC(size(setIdx)))
         BC = vertexSetProperties%BC
         Call VecSetValues(x,size(setIdx),setIdx,BC,INSERT_VALUES,ierr);CHKERRQ(ierr)
         DeAllocate(BC)
         Call ISRestoreIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
      End If
   End Do
   Call ISRestoreIndicesF90(PoissonCtx%VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Call VecAssemblyBegin(x,ierr);CHKERRQ(ierr)
   Call VecAssemblyEnd(x,ierr);CHKERRQ(ierr)
End Subroutine SimplePoissonFormInitialGuess

#undef __FUNCT__
#define __FUNCT__ "SimplePoissonComputeEnergies"
!!!
!!!  
!!!  SimplePoissonComputeEnergies:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine SimplePoissonComputeEnergies(x,PoissonCtx,energy,work,ierr)
   Type(Vec),Intent(IN)                            :: x
   Type(PoissonCtx_Type),intent(IN)                :: PoissonCtx
   PetscReal,Dimension(:),Pointer                  :: energy,work
   PetscErrorCode,Intent(OUT)                      :: ierr

   Type(IS)                                        :: setIS,setISdof
   PetscInt,dimension(:),Pointer                   :: setID
   PetscInt,Dimension(:),Pointer                   :: setIdx
   PetscInt                                        :: set
   Type(MEF90_MATPROP),pointer                     :: matpropSet
   Type(PoissonCellSetProperties_Type),pointer     :: cellSetProperties
   PetscReal                                       :: myenergy,mywork
   Type(Element_Type)                              :: elemType

   energy = 0.0_Kr
   work = 0.0_Kr
   
   Call SectionRealToVec(PoissonCtx%Section,PoissonCtx%ScatterSecToVec,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr)
   Call ISGetIndicesF90(PoissonCtx%CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      myenergy = 0.0_Kr
      mywork   = 0.0_Kr
      Call DMMeshGetStratumIS(PoissonCtx%mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)

      Call PetscBagGetDataMEF90_MatProp(PoissonCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataPoissonCellSetProperties(PoissonCtx%CellSetPropertiesBag(set),cellSetProperties,ierr);CHKERRQ(ierr)
      Call Element_TypeFindByID(cellSetProperties%ElemTypeShortID,elemType) 

      !!! Assembly part of the residual coming from the bilinear form on the blocks of 
      !!! codimension 0
      If (elemType%coDim == 0) Then
         Call MEF90_DiffusionEnergySet(myenergy,PoissonCtx%Section,PoissonCtx%mesh,matpropSet%ThermalDiffusivity,0.0_Kr,setIS,PoissonCtx%elem,elemType,ierr)
      End If
      Call MPI_AllReduce(myenergy,energy(set),1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)
      If (cellSetProperties%Force /= 0.0_Kr) Then
         Call MEF90_DiffusionWorkSet(mywork,PoissonCtx%Section,PoissonCtx%mesh,cellSetProperties%Force,setIS,PoissonCtx%elem,elemType,ierr)
      End If
      Call MPI_AllReduce(mywork,work(set),1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   End Do     
   Call ISRestoreIndicesF90(PoissonCtx%CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

   
End Subroutine SimplePoissonComputeEnergies
End Module M_POISSON
