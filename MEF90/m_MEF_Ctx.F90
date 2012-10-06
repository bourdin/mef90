Module m_MEF_Ctx
#include "finclude/petscdef.h"
   Use, Intrinsic :: iso_c_binding
   Use m_MEF_Parameters
   Implicit none

   Private  
   Public :: MEF90Ctx_Type
   Public :: MEF90Ctx_InitializePrivate
   Public :: MEF90Ctx_GetTime
      
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
                     MEF90TimeInterpolation_exo
                     
                     
      !!! 
      !!! is file command line, exo data file 
      !!! -time 0,1 => two steps
      !!! -time 0 -time_dt .1 -time_numsteps 11
      !!! -time 0,.1,.2,.3,.4,.5,.6,.7,.8,.9 is easy since it will go in an YAML file
      !!! -time_exo => read in exodus file
      !!! -time_txt => read in external file
   End Enum
   Character(len=MEF90_MXSTRLEN),dimension(6),protected  :: MEF90TimeInterpolationList

   Type MEF90Ctx_Type
      PetscInt                                        :: verbose
      Character(len=MEF90_MXSTRLEN)                   :: prefix
      PetscEnum                                       :: timeInterpolation
      PetscReal                                       :: timeMin,timeMax
      PetscInt                                        :: timeNumStep
      PetscEnum                                       :: fileFormat
      PetscEnum                                       :: fileMode
      Integer                                         :: fileExoUnitIn
      Integer                                         :: fileExoUnitOut

      PetscReal,Dimension(:),Pointer                  :: time !!! This can't be in a bag since the size of the bag needs to be known a priori!
                                                              !!! It will have to be a global variable in all applications
                                                              !!! And I will need to implement MEF90GetTimeArray(t,MEF90Ctx)
      
!!!   Keep for compatibility reasons until MEF90HeatXfer is working
      PetscBag                                        :: GlobalPropertiesBag
      PetscBag,Dimension(:),Pointer                   :: CellSetPropertiesBag
      PetscBag,Dimension(:),Pointer                   :: VertexSetPropertiesBag
      PetscBag,Dimension(:),Pointer                   :: MaterialPropertiesBag
   End Type MEF90Ctx_Type
   
Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90Ctx_InitializePrivate"
!!!
!!!  
!!!  MEF90Ctx_InitializePrivate:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90Ctx_InitializePrivate(ierr)
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
      MEF90TimeInterpolationList(3) = 'exo'
      MEF90TimeInterpolationList(4) = 'MEF90TimeInterpolation'
      MEF90TimeInterpolationList(5) = '_MEF90TimeInterpolation'
      MEF90TimeInterpolationList(6) = ''
      
      ierr = 0
   End Subroutine MEF90Ctx_InitializePrivate
   
#undef __FUNCT__
#define __FUNCT__ "MEF90Ctx_GetTime"
!!!
!!!  
!!!  MEF90Ctx_GetTimeArray:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90Ctx_GetTime(MEF90Ctx,t,ierr)
      Type(MEF90Ctx_Type),Intent(INOUT)               :: MEF90Ctx
      PetscReal,Dimension(:),Pointer                  :: t
      PetscErrorCode,Intent(OUT)                      :: ierr
      
      PetscReal                                       :: dt
      Integer                                         :: i,n
      Real                                            :: dummyR
      Character(len=1)                                :: dummyS
      Integer                                         :: exoerr
      
      
      SelectCase(MEF90Ctx%timeInterpolation)
      Case (MEF90TimeInterpolation_linear)
         Allocate(t(MEF90Ctx%timeNumStep))
         dt = (MEF90Ctx%timeMax - MEF90Ctx%timeMin) / Real(MEF90Ctx%timeNumStep)
         t = (/ (MEF90Ctx%timeMin + Real(i) * dt, i = 0,MEF90Ctx%timeNumStep-1) /)
         t(MEF90Ctx%timeNumStep) = MEF90Ctx%timeMax
      Case (MEF90TimeInterpolation_quadratic)
         !!! Natural time scale for the heat equation
         Allocate(t(MEF90Ctx%timeNumStep))
         dt = (MEF90Ctx%timeMax**2 - MEF90Ctx%timeMin**2) / Real(MEF90Ctx%timeNumStep)
         t = (/ (sqrt(MEF90Ctx%timeMin**2 + Real(i) * dt), i = 0,MEF90Ctx%timeNumStep-1) /)
         t(MEF90Ctx%timeNumStep) = MEF90Ctx%timeMax
      Case (MEF90TimeInterpolation_exo)
         Call EXINQ(MEF90Ctx%fileExoUnitIn,EXTIMS,MEF90Ctx%timeNumStep,dummyR,dummyS,exoerr)
         Allocate(t(MEF90Ctx%timeNumStep))
         Call EXGATIM(MEF90Ctx%fileExoUnitIn,t,exoerr)
         MEF90Ctx%timeMin = t(1)
         MEF90Ctx%timeMax = t(MEF90Ctx%timeNumStep)
      End Select
      ierr = 0
   End Subroutine MEF90Ctx_GetTime
End Module m_MEF_Ctx
