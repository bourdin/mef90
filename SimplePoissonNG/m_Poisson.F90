#include "SimplePoisson.inc"

Module M_POISSON_TYPES
#include "finclude/petscdef.h"
   Use m_MEF90
   Implicit none
   Private  
   Public :: MEF90Ctx_Type
   Public :: PoissonCellSetProperties_Type
   Public :: PoissonVertexSetProperties_Type
   
   Type MEF90Ctx_Type
      !Type(IS)                                        :: CellSetGlobalIS,VertexSetGlobalIS
      !Type(IS)                                        :: CellSetLocalIS,VertexSetLocalIS
      PetscBag                                        :: GlobalPropertiesBag
      PetscBag,Dimension(:),Pointer                   :: CellSetPropertiesBag
      PetscBag,Dimension(:),Pointer                   :: VertexSetPropertiesBag
      PetscBag,Dimension(:),Pointer                   :: MaterialPropertiesBag
   End Type MEF90Ctx_Type

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
   Private
  
   Public :: MEF90Ctx_Type
    
   Public :: m_Poisson_Initialize
   Public :: PetscBagRegisterPoissonCellSetProperties
   Public :: PetscBagGetDataPoissonCellSetProperties
   Public :: PetscBagRegisterPoissonVertexSetProperties
   Public :: PetscBagGetDataPoissonVertexSetProperties
   Public :: PoissonCtxCreate
   Public :: PoissonCtxDestroy
   Public :: SimplePoissonBilinearForm
   Public :: SimplePoissonOperator
   Public :: SimplePoissonEnergies
   Public :: SimplePoissonFormInitialGuess
   
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
      Type(MEF90Ctx_Type),Target                   :: PoissonCtx
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
      !Call PetscBagRegisterEnum(bag,CellSetProperties%ElemTypeShortIDEnum,MEF90_knownElementNames,default%ElemTypeShortID,'shortidenum','Element type name',ierr);CHKERRQ(ierr)
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
   Subroutine PoissonCtxCreate(PoissonCtx,snesTemp,exoid,ierr)
      Type(MEF90Ctx_Type),intent(OUT)              :: PoissonCtx
      Type(SNES),Intent(IN)                        :: snesTemp
      Integer,Intent(IN)                           :: exoid
      PetscErrorCode,Intent(OUT)                   :: ierr
     
      Type(DM)                                     :: mesh
      PetscInt                                     :: verbose
      PetscBool                                    :: flg
      Character(len=MXSTLN)                        :: filename,EXOElemType
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer,setName,setprefix
      PetscInt,Dimension(:),Pointer                :: setID
      PetscInt                                     :: set
      Type(Element_Type)                           :: defaultElemType
      Type(PoissonCellSetProperties_Type)          :: defaultCellSetProperties
      Type(PoissonCellSetProperties_Type),pointer  :: CellSetProperties
      Type(PoissonVertexSetProperties_Type)        :: defaultVertexSetProperties
      Type(IS)                                     :: CellSetGlobalIS,VertexSetGlobalIS
      Type(SectionReal)                            :: Sec
      
      verbose = 0
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'--verbose',verbose,flg,ierr);CHKERRQ(ierr)

      Call SNESGetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)
      
      !!! Get number of cell, vertices and dimension
      !Call DMmeshGetStratumSize(mesh,"depth",0,numVertex,ierr);CHKERRQ(ierr)
      !Call DMmeshGetStratumSize(mesh,"height",0,numCell,ierr);CHKERRQ(ierr)
      !Call DMMeshGetDimension(mesh,numDim,ierr);CHKERRQ(ierr)
      
      
      !!!
      !!! Allocate and register cell sets bags: CellSetProperties and MatProp
      !!!
      Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Allocate(PoissonCtx%CellSetPropertiesBag(size(setID)),stat=ierr)
      Allocate(PoissonCtx%MaterialPropertiesBag(size(setID)),stat=ierr)
      
      Do set = 1, size(setID)
         Write(setName,100) setID(set)
         Write(setprefix,101) setID(set)
         If (verbose > 0) Then
            Write(IOBuffer,103) MEF90_MyRank,setID(set),trim(setprefix)
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
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)

      !!!
      !!! Allocate and register vertex properties bags:
      !!!
      Call DMmeshGetLabelIdIS(mesh,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Allocate(PoissonCtx%VertexSetPropertiesBag(size(setID)),stat=ierr)
      defaultVertexSetProperties = PoissonVertexSetProperties_Type( PETSC_TRUE,0.0_Kr )
      Do set = 1, size(setID)
         Write(setName,200) setID(set)
         Write(setprefix,201) setID(set)
         If (verbose > 0) Then
            Write(IOBuffer,203) MEF90_MyRank,setID(set),trim(setprefix)
            Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
         End if
         Call PetscBagCreate(PETSC_COMM_WORLD,sizeofPoissonVertexSetProperties,PoissonCtx%VertexSetPropertiesBag(set),ierr)
         Call PetscBagRegisterPoissonVertexSetProperties(PoissonCtx%VertexSetPropertiesBag(set),setname,setprefix,defaultVertexSetProperties,ierr);CHKERRQ(ierr)
         If (verbose > 0) Then
            Call PetscBagView(PoissonCtx%VertexSetPropertiesBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
         End If
      End Do
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      
100 Format('Cell set ',I4)
101 Format('cs',I4.4,'_')
103 Format('[',I4.4,'] Registering cell set ',I4,' prefix: ',A,'\n')
200 Format('Vertex set ',I4)
201 Format('vs',I4.4,'_')
203 Format('[',I4.4,'] Registering vertex set ',I4,' prefix: ',A,'\n')
   End Subroutine PoissonCtxCreate
   
#undef __FUNCT__
#define __FUNCT__ "PoissonCtxDestroy"
!!!
!!!  
!!!  PoissonCtxDestroy:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine PoissonCtxDestroy(PoissonCtx,snesTemp,ierr)
   Type(MEF90Ctx_Type)                                      :: PoissonCtx
   Type(SNES),Intent(IN)                                    :: snesTemp
   PetscErrorCode,Intent(OUT)                               :: ierr

   Type(IS)                                                 :: setIS   
   Type(DM)                                                 :: mesh
   PetscInt                                                 :: e,set,nset
   
   Call SNESGetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)

   Call DMmeshGetLabelIdIS(mesh,'Cell Sets',SetIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr) 
   Call ISGetLocalSize(setIS,nset,ierr);CHKERRQ(ierr)
   Do set = 1, nset
      Call PetscBagDestroy(PoissonCtx%CellSetPropertiesBag(set),ierr);CHKERRQ(ierr)
      Call PetscBagDestroy(PoissonCtx%MaterialPropertiesBag(set),ierr);CHKERRQ(ierr)
   End Do
   DeAllocate(PoissonCtx%CellSetPropertiesBag,stat=ierr)
   DeAllocate(PoissonCtx%MaterialPropertiesBag,stat=ierr)
   Call ISDestroy(setIS,ierr);CHKERRQ(ierr)

   Call DMmeshGetLabelIdIS(mesh,'Vertex Sets',SetIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr) 
   Call ISGetLocalSize(setIS,nset,ierr);CHKERRQ(ierr)
   Do set = 1, nset
      Call PetscBagDestroy(PoissonCtx%VertexSetPropertiesBag(set),ierr);CHKERRQ(ierr)
   End Do
   DeAllocate(PoissonCtx%VertexSetPropertiesBag)
   Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   
   !Call PetscBagDestroy(PoissonCtx%GlobalPropertiesBag,ierr);CHKERRQ(ierr)
End Subroutine PoissonCtxDestroy 

#undef __FUNCT__
#define __FUNCT__ "SimplePoissonBilinearForm"
!!!
!!!  
!!!  SimplePoissonBilinearForm:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine SimplePoissonBilinearForm(snesTemp,x,A,M,flg,PoissonCtx,ierr)
   Type(SNES),Intent(IN)                           :: snesTemp
   Type(Vec),Intent(IN)                            :: x
   Type(Mat),Intent(INOUT)                         :: A,M
   MatStructure,Intent(INOUT)                      :: flg
   Type(MEF90Ctx_Type),intent(IN)                  :: PoissonCtx
   PetscErrorCode,Intent(OUT)                      :: ierr  
   
   Type(IS)                                        :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof
   PetscInt,dimension(:),Pointer                   :: setID
   PetscInt                                        :: set,QuadratureOrder
   Type(MEF90_MATPROP),pointer                     :: matpropSet
   Type(PoissonCellSetProperties_Type),pointer     :: cellSetProperties
   Type(PoissonVertexSetProperties_Type),pointer   :: vertexSetProperties
   Type(DM)                                        :: mesh
   Type(SectionReal)                               :: xSec
   Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer   :: elem
   Type(Element_Type)                              :: elemType
   
   Call SNESGetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)
   Call DMMeshGetSectionReal(mesh,'default',xSec,ierr);CHKERRQ(ierr)
   !!! 
   !!! No need to zero out Sec because it is only used as a PetscSection when assembling Mat
   !!!
   Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
   Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)

      Call PetscBagGetDataMEF90_MatProp(PoissonCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataPoissonCellSetProperties(PoissonCtx%CellSetPropertiesBag(set),cellSetProperties,ierr);CHKERRQ(ierr)

      elemType = MEF90_knownElements(cellSetProperties%ElemTypeShortID)
      If (elemType%coDim == 0) Then
         !lambda = 0 in MEF90_DiffusionBilinearFormSet so the polynomial order is indeed 2(p-1) where p is the 
         ! polynomial order of the element
         QuadratureOrder = (elemType%order-1) * 2
         Call MEF90_ElementCreate(mesh,setIS,elem,QuadratureOrder,CellSetProperties%ElemTypeShortID,ierr);CHKERRQ(ierr)
         Call MEF90_DiffusionBilinearFormSet(A,mesh,xSec,setIS,matpropSet%ThermalDiffusivity,0.0_Kr,elem,elemType,ierr);CHKERRQ(ierr)
         Call MEF90_ElementDestroy(elem,ierr)
      End If
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   End Do     
   Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
   Call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)

   !!!
   !!! Boundary conditions
   !!!
   Call DMmeshGetLabelIdIS(mesh,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
   Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      Call PetscBagGetDataPoissonVertexSetProperties(PoissonCtx%VertexSetPropertiesBag(set),vertexSetProperties,ierr);CHKERRQ(ierr)
      If (vertexSetProperties%Has_BC) Then
         Call DMMeshGetStratumIS(mesh,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call DMMeshISCreateISglobaldof(mesh,xSec,setIS,0,setISdof,ierr);CHKERRQ(ierr)
         Call MatZeroRowsColumnsIS(A,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
      End If
   End Do
   Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
   
   flg = SAME_NONZERO_PATTERN
End Subroutine SimplePoissonBilinearForm

#undef __FUNCT__
#define __FUNCT__ "SimplePoissonOperator"
!!!
!!!  
!!!  SimplePoissonOperator:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine SimplePoissonOperator(snesTemp,x,residual,PoissonCtx,ierr)
   Type(SNES),Intent(IN)                           :: snesTemp
   Type(Vec),Intent(IN)                            :: x
   Type(Vec),Intent(INOUT)                         :: residual
   Type(MEF90Ctx_Type),intent(IN)                  :: PoissonCtx
   PetscErrorCode,Intent(OUT)                      :: ierr

   Type(SectionReal)                               :: xSec,residualSec
   Type(IS)                                        :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof
   PetscInt,dimension(:),Pointer                   :: setID
   PetscInt,Dimension(:),Pointer                   :: setIdx
   PetscInt                                        :: set,QuadratureOrder
   Type(MEF90_MATPROP),pointer                     :: matpropSet
   Type(PoissonCellSetProperties_Type),pointer     :: cellSetProperties
   Type(PoissonVertexSetProperties_Type),pointer   :: vertexSetProperties
   PetscReal,Dimension(:),Pointer                  :: BC
   PetscReal                                       :: F
   Type(DM)                                        :: mesh
   Type(VecScatter)                                :: ScatterSecToVec
   Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer   :: elem
   Type(Element_Type)                              :: elemType
   
   Call SNESGetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)

   Call DMMeshGetSectionReal(mesh,'default',xSec,ierr);CHKERRQ(ierr)
   Call DMMeshCreateGlobalScatter(mesh,xSec,ScatterSecToVec,ierr);CHKERRQ(ierr)

   Call SectionRealDuplicate(xSec,residualSec,ierr);CHKERRQ(ierr)
   Call SectionRealSet(residualSec,0.0_Kr,ierr);CHKERRQ(ierr)
   Call VecSet(residual,0.0_Kr,ierr);CHKERRQ(ierr)

   Call SectionRealToVec(xSec,ScatterSecToVec,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr)

   Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
   Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)

      Call PetscBagGetDataMEF90_MatProp(PoissonCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataPoissonCellSetProperties(PoissonCtx%CellSetPropertiesBag(set),cellSetProperties,ierr);CHKERRQ(ierr)

      elemType = MEF90_knownElements(cellSetProperties%ElemTypeShortID)
      QuadratureOrder = elemType%order * 2
      Call MEF90_ElementCreate(mesh,setIS,elem,QuadratureOrder,CellSetProperties%ElemTypeShortID,ierr);CHKERRQ(ierr)
      !!! Assembly part of the residual coming from the bilinear form on the blocks of 
      !!! codimension 0
      If (MEF90_knownElements(cellSetProperties%ElemTypeShortID)%coDim == 0) Then
         Call MEF90_DiffusionOperatorSet(residualSec,mesh,xSec,setIS,matpropSet%ThermalDiffusivity,0.0_Kr,elem,elemType,ierr);CHKERRQ(ierr)
      End If
      !!! Assembly the force part
      If (cellSetProperties%Force /= 0.0_Kr) Then
         F = -cellSetProperties%Force
         Call MEF90_DiffusionRHSSet(residualSec,mesh,F,setIS,elem,elemType,ierr);CHKERRQ(ierr)
      End If
      Call MEF90_ElementDestroy(elem,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   End Do     
   Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
   !!! "Ghost update" for the residual Section
   Call SectionRealComplete(residualSec,ierr);CHKERRQ(ierr)
   !!! Scatter back from SectionReal to Vec
   Call SectionRealToVec(residualSec,ScatterSecToVec,SCATTER_FORWARD,residual,ierr);CHKERRQ(ierr)
   
   
   !!!
   !!! Zero-out BC entries in the residual
   !!!
   Call DMmeshGetLabelIdIS(mesh,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
   Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      Call PetscBagGetDataPoissonVertexSetProperties(PoissonCtx%VertexSetPropertiesBag(set),vertexSetProperties,ierr);CHKERRQ(ierr)
      If (vertexSetProperties%Has_BC) Then
         Call DMMeshGetStratumIS(mesh,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call DMMeshISCreateISglobaldof(mesh,xSec,setIS,0,setISdof,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
         Allocate(BC(size(setIdx)))
         BC = 0.0_Kr
         Call VecSetValues(residual,size(setIdx),setIdx,BC,INSERT_VALUES,ierr);CHKERRQ(ierr)
         DeAllocate(BC)
         Call ISRestoreIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
      End If
   End Do
   Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)

   Call VecAssemblyBegin(residual,ierr)
   Call VecAssemblyEnd(residual,ierr)
   Call SectionRealDestroy(residualSec,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
   Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
End Subroutine SimplePoissonOperator


#undef __FUNCT__
#define __FUNCT__ "SimplePoissonFormInitialGuess"
!!!
!!!  
!!!  SimplePoissonFormInitialGuess:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine SimplePoissonFormInitialGuess(snesTemp,x,PoissonCtx,ierr)
   Type(SNES),Intent(IN)                           :: snesTemp
   Type(Vec),Intent(IN)                            :: x
   Type(MEF90Ctx_Type),intent(IN)                  :: PoissonCtx
   PetscErrorCode,Intent(OUT)                      :: ierr

   Type(PoissonVertexSetProperties_Type),pointer   :: vertexSetProperties
   Type(IS)                                        :: VertexSetGlobalIS,setIS,setISdof
   PetscInt,dimension(:),Pointer                   :: setID
   PetscInt,Dimension(:),Pointer                   :: setIdx
   PetscInt                                        :: set
   PetscReal,Dimension(:),Pointer                  :: BC
   Type(DM)                                        :: mesh
   Type(SectionReal)                               :: xSec
   
   Call SNESGetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)
   Call DMMeshGetSectionReal(mesh,'default',xSec,ierr);CHKERRQ(ierr)

   Call DMmeshGetLabelIdIS(mesh,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
   Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      Call PetscBagGetDataPoissonVertexSetProperties(PoissonCtx%VertexSetPropertiesBag(set),vertexSetProperties,ierr);CHKERRQ(ierr)
      If (vertexSetProperties%Has_BC) Then
         Call DMMeshGetStratumIS(mesh,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call DMMeshISCreateISglobaldof(mesh,xSec,setIS,0,setISdof,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
         Allocate(BC(size(setIdx)))
         BC = vertexSetProperties%BC
         Call VecSetValues(x,size(setIdx),setIdx,BC,INSERT_VALUES,ierr);CHKERRQ(ierr)
         DeAllocate(BC)
         Call ISRestoreIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
      End If
   End Do
   Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)

   Call VecAssemblyBegin(x,ierr);CHKERRQ(ierr)
   Call VecAssemblyEnd(x,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
End Subroutine SimplePoissonFormInitialGuess

#undef __FUNCT__
#define __FUNCT__ "SimplePoissonEnergies"
!!!
!!!  
!!!  SimplePoissonEnergies:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine SimplePoissonEnergies(snesTemp,x,PoissonCtx,energy,work,ierr)
   Type(SNES),Intent(IN)                           :: snesTemp
   Type(Vec),Intent(IN)                            :: x
   Type(MEF90Ctx_Type),intent(IN)                  :: PoissonCtx
   PetscReal,Dimension(:),Pointer                  :: energy,work
   PetscErrorCode,Intent(OUT)                      :: ierr

   Type(IS)                                        :: CellSetGlobalIS,setIS,setISdof
   PetscInt,dimension(:),Pointer                   :: setID
   PetscInt,Dimension(:),Pointer                   :: setIdx
   PetscInt                                        :: set,QuadratureOrder
   Type(MEF90_MATPROP),pointer                     :: matpropSet
   Type(PoissonCellSetProperties_Type),pointer     :: cellSetProperties
   PetscReal                                       :: myenergy,mywork
   Type(DM)                                        :: mesh
   Type(SectionReal)                               :: xSec
   Type(VecScatter)                                :: ScatterSecToVec
   Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer   :: elem
   Type(Element_Type)                              :: elemType

   energy = 0.0_Kr
   work = 0.0_Kr   
   Call SNESGetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)
   Call DMMeshGetSectionReal(mesh,'default',xSec,ierr);CHKERRQ(ierr)
   Call DMMeshCreateGlobalScatter(mesh,xSec,ScatterSecToVec,ierr);CHKERRQ(ierr)
   Call SectionRealToVec(xSec,ScatterSecToVec,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr)

   Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)   
   Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      myenergy = 0.0_Kr
      mywork   = 0.0_Kr
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)

      Call PetscBagGetDataMEF90_MatProp(PoissonCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataPoissonCellSetProperties(PoissonCtx%CellSetPropertiesBag(set),cellSetProperties,ierr);CHKERRQ(ierr)

      elemType = MEF90_knownElements(cellSetProperties%ElemTypeShortID)
      QuadratureOrder = elemType%order * 2
      Call MEF90_ElementCreate(mesh,setIS,elem,QuadratureOrder,CellSetProperties%ElemTypeShortID,ierr);CHKERRQ(ierr)
      !!! Assembly part of the residual coming from the bilinear form on the blocks of 
      !!! codimension 0
      If (MEF90_knownElements(cellSetProperties%ElemTypeShortID)%coDim == 0) Then
         Call MEF90_DiffusionEnergySet(myenergy,xSec,mesh,matpropSet%ThermalDiffusivity,0.0_Kr,setIS,elem,elemType,ierr)
      End If
      Call MPI_AllReduce(myenergy,energy(set),1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)
      If (cellSetProperties%Force /= 0.0_Kr) Then
         Call MEF90_DiffusionWorkSet(mywork,xSec,mesh,cellSetProperties%Force,setIS,elem,elemType,ierr)
      End If
      Call MPI_AllReduce(mywork,work(set),1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)
      Call MEF90_ElementDestroy(elem,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   End Do     
   Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)

   Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
   Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
End Subroutine SimplePoissonEnergies
End Module M_POISSON
