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
                     MEF90Scaling_File
   End Enum
   Character(len=MEF90_MXSTRLEN),dimension(6),protected  :: MEF90ScalingList
   
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
   
   Enum,bind(c)
      Enumerator  :: MEF90TimeInterpolation_linear = 0,  &
                     MEF90TimeInterpolation_quadratic,   &
                     MEF90TimeInterpolation_logarithmic, &
                     MEF90TimeInterpolation_file
      !!! is file command line, exo data file 
      !!! -time 0,1 => two steps
      !!! -time 0,1 -time_dt .1 / -time_numsteps 11
      !!! -time 0,.1,.2,.3,.4,.5,.6,.7,.8,.9 is easy since it will go in an YAML file
      !!! -time_exo => read in exodus file
   End Enum
   Character(len=MEF90_MXSTRLEN),dimension(7),protected  :: MEF90TimeInterpolationList

   Type MEF90Ctx_Type
      PetscInt                                        :: verbose
      Character(len=MEF90_MXSTRLEN)                   :: prefix
      PetscReal,Dimension(:),Pointer                  :: time
      PetscReal,Target                                :: timeMin,timeMax
      PetscInt,Target                                 :: timeNumStep
      PetscEnum,Target                                :: timeInterpolation
      PetscEnum                                       :: fileFormat
      PetscEnum                                       :: fileMode
      Integer                                         :: fileExoUnitIn
      Integer                                         :: fileExoUnitOut

!!!   Keep for compatibility reasons until MEF90HeatXfer is working
      PetscBag                                        :: GlobalPropertiesBag
      PetscBag,Dimension(:),Pointer                   :: CellSetPropertiesBag
      PetscBag,Dimension(:),Pointer                   :: VertexSetPropertiesBag
      PetscBag,Dimension(:),Pointer                   :: MaterialPropertiesBag
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
      MEF90ScalingList(3) = 'file'
      MEF90ScalingList(4) = 'MEF90scaling'
      MEF90ScalingList(5) = '_MEF90Scaling'
      MEF90ScalingList(6) = ''
      
      MEF90FileModeList(1) = 'Replace'
      MEF90FileModeList(2) = 'Append'
      MEF90FileModeList(3) = 'MEF90FileMode'
      MEF90FileModeList(4) = '_MEF90FileMode'
      MEF90FileModeList(5) = ''

      MEF90EXOModeList(1) = 'Single'
      MEF90EXOModeList(2) = 'Split'
      MEF90EXOModeList(3) = 'MEF90EXOMode'
      MEF90EXOModeList(4) = '_MEF90EXOMode'
      MEF90EXOModeList(5) = '' 
      
      MEF90TimeInterpolationList(1) = 'linear'
      MEF90TimeInterpolationList(2) = 'quadratic'
      MEF90TimeInterpolationList(3) = 'logarithmic'
      MEF90TimeInterpolationList(4) = 'file'
      MEF90TimeInterpolationList(5) = 'MEF90TimeInterpolation'
      MEF90TimeInterpolationList(6) = '_MEF90TimeInterpolation'
      MEF90TimeInterpolationList(7) = ''
      
      ierr = 0
   End Subroutine MEF90_CtxInitialize
   
End Module m_MEF_Ctx

!  Shall I move verbose, fileformat, filemode and all other time informations to a MEF90
!  Bag? This is very likely to be shared between applications
!  Global, CellSet, VertexSet Properties should be moved into the appctx
