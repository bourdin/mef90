#if MEF90_DIM == 3
#define M_POISSON_TYPES m_Poisson3D_Types
#else
#define M_POISSON_TYPES m_Poisson2D_Types
#endif

Module M_POISSON_TYPES
#include "../MEF90/mef90.inc"
#include "finclude/petscdef.h"
   Use m_MEF90
   Implicit none
   Private  
   Public :: PoissonCtx_Type
   Public :: PoissonCellSetProperties_Type
   Public :: PoissonVertexSetProperties_Type
   
   Type PoissonCtx_Type
      Type(DM)                                        :: mesh
      Type(IS)                                        :: CellSetGlobalIS,VertexSetGlobalIS
      Type(IS)                                        :: CellSetLocalIS,VertexSetLocalIS
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer   :: Elem
      !Type(Element_Type),Dimension(:),Pointer         :: ElementType ! move to CellSetProperty_Type
      Type(Field)                                     :: BCU
      Type(Field)                                     :: F !Make forces cell centered in assembly routines
   End Type PoissonCtx_Type

   Type PoissonCellSetProperties_Type   
      Character(len=MEF90_MXSTRLEN)                   :: ElemTypeName
      PetscReal                                       :: Force
   End Type PoissonCellSetProperties_Type   
   
   Type PoissonVertexSetProperties_Type   
      PetscBool                                       :: Has_BC
      PetscReal                                       :: BC
   End Type PoissonVertexSetProperties_Type   
End Module M_POISSON_TYPES


#if MEF90_DIM == 3
#define M_POISSONCELLSETPROPERTY_INTERFACE m_PoissonCellSetProperty3D_Interface
#else
#define M_POISSONCELLSETPROPERTY_INTERFACE m_PoissonCellSetProperty2D_Interface
#endif
Module M_POISSONCELLSETPROPERTY_INTERFACE   
#include "../MEF90/mef90.inc"
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
      PetscBag                                     :: bag
      type(PoissonCellSetProperties_Type),pointer  :: data
      PetscErrorCode                               :: ierr
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataPoissonCellSetProperties
End Module M_POISSONCELLSETPROPERTY_INTERFACE   
!!!
!!!   
   
#if MEF90_DIM == 3
#define M_POISSON m_Poisson3D
#else
#define M_POISSON m_Poisson2D
#endif

Module M_POISSON
#include "../MEF90/mef90.inc"
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use M_POISSON_TYPES
   Use M_POISSONCELLSETPROPERTY_INTERFACE
   Use m_MEF90
   Implicit NONE
   Private
   
   !Public :: m_Poisson_Initialize
   Public :: PetscBagRegisterPoissonCellSetProperties
   Public :: PetscBagGetDataPoissonCellSetProperties
   !Public :: PetscBagRegisterPoissonVertexSetProperties
   !Public :: PetscBagGetDataPoissonVertexSetProperties
   !Public :: PetscBagRegisterPoissonCtx
   !Public :: PetscBagGetDataPoissonCtx
   
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
      Call PetscBagSetName(bag,trim(name),"CellSetProperties object: Cell Set properties",ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix), ierr)
      Call PetscBagRegisterEnum(bag,CellSetProperties%ElemTypeName,MEF90_knownElementNames, &
                                MEF90_P1_Lagrange_2D_Scal%shortID,'Element type')
      Call PetscBagRegisterReal(bag,CellSetProperties%Force,0.0_Kr,'Force','Magnitude of the applied force',ierr);CHKERRQ(ierr)
   End Subroutine PetscBagRegisterPoissonCellSetProperties
End Module M_POISSON
