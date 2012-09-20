Module m_MEF_Ctx
#include "finclude/petscdef.h"
   Use, Intrinsic :: iso_c_binding
   Use m_MEF_Parameters
   Implicit none

   Private  
   Public :: MEF90Ctx_Type
   Public :: MEF90_CtxInitialize
      
   Enum,bind(c)
      Enumerator ::  MEF90Scaling_CST=0,        &
                     MEF90Scaling_Linear,       &  
                     MEF90Scaling_Quadratic,    &
                     MEF90Scaling_Logarithmic,  &
                     MEF90Scaling_txtFile,      &
                     MEF90Scaling_dataFile
   End Enum
   Character(len=MEF90_MXSTRLEN),dimension(9),protected  :: MEF90ScalingList
   
   Enum,bind(c)
      Enumerator ::  MEF90FileMode_Replace = 0, &
                     MEF90FileMode_Append
   End Enum
   Character(len=MEF90_MXSTRLEN),dimension(5),protected  :: MEF90FileModeList
   
   Enum,bind(c)
      Enumerator  :: MEF90EXOMode_Single = 0,    &
                     MEF90EXOMode_Split 
   End Enum
   Character(len=MEF90_MXSTRLEN),dimension(5),protected  :: MEF90EXOModeList

   Type MEF90Ctx_Type
      !PetscBag                                        :: MEF90Bag
      PetscBag                                        :: GlobalPropertiesBag
      PetscBag,Dimension(:),Pointer                   :: CellSetPropertiesBag
      PetscBag,Dimension(:),Pointer                   :: VertexSetPropertiesBag
      PetscBag,Dimension(:),Pointer                   :: MaterialPropertiesBag
      Type(C_Ptr)                                     :: AppCtx
   End Type MEF90Ctx_Type
   
Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_CtxInitialize"
!!!
!!!  
!!!  MEF90_CtxInitialize:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90_CtxInitialize(ierr)
      PetscErrorCode,Intent(OUT)                      :: ierr
   
      MEF90ScalingList(1) = 'constant'
      MEF90ScalingList(2) = 'linear'
      MEF90ScalingList(3) = 'quadratic'
      MEF90ScalingList(4) = 'logarithmic'
      MEF90ScalingList(5) = 'txtfile'
      MEF90ScalingList(6) = 'datafile'
      MEF90ScalingList(7) = 'scaling'
      MEF90ScalingList(8) = '_MEF90Scaling'
      MEF90ScalingList(9) = ''
      
      MEF90FileModeList(1) = 'Replace'
      MEF90FileModeList(2) = 'Append'
      MEF90FileModeList(3) = 'FileMode'
      MEF90FileModeList(4) = '_FileMode'
      MEF90FileModeList(5) = ''

      MEF90EXOModeList(1) = 'Single'
      MEF90EXOModeList(2) = 'Split'
      MEF90EXOModeList(3) = 'EXOMode'
      MEF90EXOModeList(4) = '_EXOMode'
      MEF90EXOModeList(5) = '' 
      ierr = 0
   End Subroutine MEF90_CtxInitialize
   
End Module m_MEF_Ctx

!  Shall I move verbose, fileformat, filemode and all other time informations to a MEF90
!  Bag? This is very likely to be shared between applications
!  Global, CellSet, VertexSet Properties should be moved into the appctx