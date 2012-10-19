#include "../MEF90/mef90.inc"
Module m_MEF90_HeatXferCtx_Type
#include "finclude/petscdef.h"
   Use m_MEF90
   Use,Intrinsic :: iso_c_binding
   Implicit none
   Private  
   Public :: MEF90_HeatXferCtx_Type
   
   PetscInt,protected   :: sizeofHeatXferGlobalOptions
   PetscInt,protected   :: sizeofHeatXferCellSetOptions
   PetscInt,protected   :: sizeofHeatXferVertexOptions
   
   Enum,bind(c)
      enumerator  :: MEF90_HeatXFer_ModeSteadyState = 0, &
                     MEF90_HeatXFer_ModeSteadyTransient
   End Enum
   Character(len = MEF90_MXSTRLEN),dimension(5),protected   :: MEF90_HeatXFer_ModeList
   
   Type MEF90_HeatXferCtx_Type
      PetscReal                        :: timeOld,timeTarget

      Type(Vec),pointer                :: fluxOld,fluxTarget,Flux
      Type(Vec),pointer                :: BoundaryValueOld,BoundaryValueTarget,BoundaryValue
      Type(Vec),pointer                :: ExternalValue,interpolationExternalValue
      !!! XXXOld     represents a field at the previous time step
      !!! XXXTarget  represents a field at the time step currently being computed
      !!! XXX        represents a field at current time, interpolated from XXXOld and XXXTarget
      
      PetscBag                         :: GlobalOptions
      PetscBag,Dimension(:),Pointer    :: CellSetOptions
      PetscBag,Dimension(:),Pointer    :: VertexSetOptions
   End Type MEF90_HeatXferCtx_Type
   
   Type MEF90_HeatXferGlobalOptions_Type
      PetscInt                         :: mode
      PetscBool                        :: addNullSpace
      PetscInt                         :: tempOffset
      PetscInt                         :: boundaryTempOffset
      PetscInt                         :: externalTempOffset
      PetscInt                         :: fluxOffset
   End Type MEF90_HeatXferGlobalOptions_Type

   Type MEF90_HeatXferCellSetOptions_Type
   End Type MEF90_HeatXferCellSetOptions_Type

   Type MEF90_HeatXferVertexSetOptions_Type
   End Type MEF90_HeatXferVertexSetOptions_Type
   
Contains 
#undef __FUNCT__
#define __FUNCT__ "MEF90_HeatXferCtxInitialize"
!!!
!!!  
!!!  MEF90_HeatXferCtxInitialize:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90_HeatXferCtxInitialize(ierr)
      PetscErrorCode,Intent(OUT)                      :: ierr
   
      Type(MEF90_HeatXferGlobalOptions_Type)          :: HeatXferGlobalOptions
      Type(MEF90_HeatXferCellSetOptions_Type)         :: HeatXferCellSetOptions
      Type(MEF90_HeatXferVertexSetOptions_Type)       :: HeatXferVertexSetOptions
      character(len=1),pointer                        :: dummychar(:)
      PetscSizeT                                      :: sizeofchar
      
      Call PetscDataTypeGetSize(PETSC_CHAR,sizeofchar,ierr)
      sizeofHeatXferGlobalOptions = size(transfer(HeatXferGlobalOptions,dummychar))*sizeofchar
      sizeofHeatXferCellSetOptions = size(transfer(HeatXferCellSetOptions,dummychar))*sizeofchar
      sizeofHeatXferVertexOptions = size(transfer(HeatXferVertexSetOptions,dummychar))*sizeofchar

      MEF90_HeatXFer_ModeList(1) = 'SteadyState'
      MEF90_HeatXFer_ModeList(2) = 'Transient'
      MEF90_HeatXFer_ModeList(3) = 'MEF90_HeatXFer_Mode'
      MEF90_HeatXFer_ModeList(4) = '_MEF90_HeatXFer_Mode'
      MEF90_HeatXFer_ModeList(5) = ''
   End Subroutine MEF90_HeatXferCtxInitialize

End Module m_MEF90_HeatXferCtx_Type
Module m_MEF90_HeatXferCtx
End Module m_MEF90_HeatXferCtx
