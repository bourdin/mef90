Module m_MEF_Ctx_Type
#include "finclude/petscdef.h"
   Use m_MEF_Parameters
   Use,Intrinsic :: iso_c_binding
   Implicit none
   Private  
   Public :: MEF90Ctx_Type

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

      !PetscReal,Dimension(:),Pointer                  :: time !!! This can't be in a bag since the size of the bag needs to be known a priori!
                                                              !!! It will have to be a global variable in all applications
                                                              !!! And I will need to implement MEF90GetTimeArray(t,MEF90Ctx)
      
!!!   Keep for compatibility reasons until MEF90HeatXfer is working
      PetscBag                                        :: GlobalPropertiesBag
      PetscBag,Dimension(:),Pointer                   :: CellSetPropertiesBag
      PetscBag,Dimension(:),Pointer                   :: VertexSetPropertiesBag
      PetscBag,Dimension(:),Pointer                   :: MaterialPropertiesBag
   End Type MEF90Ctx_Type
End Module m_MEF_Ctx_Type

Module m_MEF_Ctx
#include "finclude/petscdef.h"
   Use, Intrinsic :: iso_c_binding
   Use m_MEF_Parameters
   Use m_MEF_Ctx_Type
   Implicit none

   Private  
   Public :: MEF90Ctx_Type
   Public :: MEF90Ctx_InitializePrivate
   Public :: PetscBagGetDataMEF90Ctx
   Public :: MEF90Ctx_GetTime
   Public :: sizeofMEF90Ctx
      
   PetscSizeT,protected    :: sizeofMEF90Ctx

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
      Enumerator  :: MEF90FileFormat_EXOSingle = 0,    &
                     MEF90FileFormat_EXOSplit 
   End Enum
   Character(len=MEF90_MXSTRLEN),dimension(5),protected  :: MEF90FileFormatList
   
   Enum,bind(c)
      Enumerator  :: MEF90TimeInterpolation_linear = 0,  &
                     MEF90TimeInterpolation_quadratic,   &
                     MEF90TimeInterpolation_exo
   End Enum
   Character(len=MEF90_MXSTRLEN),dimension(6),protected  :: MEF90TimeInterpolationList

   
   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF_Ctx_Type
         PetscBag                                     :: bag
         Type(MEF90Ctx_Type),pointer                  :: data
         PetscErrorCode                               :: ierr
      End subroutine PetscBagGetData
   End interface

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
   
      Type(MEF90Ctx_Type),Target                   :: MEF90Ctx
      character(len=1),pointer                     :: dummychar(:)
      PetscSizeT                                   :: sizeofchar
   
      Call PetscDataTypeGetSize(PETSC_CHAR,sizeofchar,ierr)
      sizeofMEF90Ctx = size(transfer(MEF90Ctx,dummychar))*sizeofchar

      MEF90ScalingList(1) = 'constant'
      MEF90ScalingList(2) = 'linear'
      MEF90ScalingList(3) = 'file'
      MEF90ScalingList(4) = 'MEF90scaling'
      MEF90ScalingList(5) = '_MEF90Scaling'
      MEF90ScalingList(6) = ''
      
      MEF90FileFormatList(1) = 'Replace'
      MEF90FileFormatList(2) = 'MEF90FileFormat'
      MEF90FileFormatList(3) = '_MEF90FileFormat'
      MEF90FileFormatList(4) = ''

      MEF90FileModeList(1) = 'Replace'
      MEF90FileModeList(2) = 'Append'
      MEF90FileModeList(3) = 'MEF90FileMode'
      MEF90FileModeList(4) = '_MEF90FileMode'
      MEF90FileModeList(5) = ''

      MEF90FileFormatList(1) = 'EXOSingle'
      MEF90FileFormatList(2) = 'EXOSplit'
      MEF90FileFormatList(3) = 'MEF90FileFormat'
      MEF90FileFormatList(4) = '_MEF90FileFormatList'
      MEF90FileFormatList(5) = '' 
      
      MEF90TimeInterpolationList(1) = 'linear'
      MEF90TimeInterpolationList(2) = 'quadratic'
      MEF90TimeInterpolationList(3) = 'exo'
      MEF90TimeInterpolationList(4) = 'MEF90TimeInterpolation'
      MEF90TimeInterpolationList(5) = '_MEF90TimeInterpolation'
      MEF90TimeInterpolationList(6) = ''
      
      ierr = 0
   End Subroutine MEF90Ctx_InitializePrivate
   
#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90Ctx"
!!!
!!!  
!!!  PetscBagGetDataMEF90Ctx:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine PetscBagGetDataMEF90Ctx(bag,data,ierr)
      PetscBag                                        :: bag
      Type(MEF90Ctx_Type),pointer                     :: data
      PetscErrorCode                                  :: ierr
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataMEF90Ctx
   
#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90Ctx"
!!!
!!!  
!!!  PetscBagRegisterMEF90Ctx:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine PetscBagRegisterMEF90Ctx(bag,name,prefix,default,ierr)
      PetscBag                                        :: bag
      Character(len=*),intent(IN)                     :: prefix,name
      Type(MEF90Ctx_Type),intent(IN)                  :: default
      PetscErrorCode,intent(OUT)                      :: ierr

      Type(MEF90Ctx_Type),pointer                     :: MEF90Ctx
      
      Call PetscBagGetDataMEF90Ctx(bag,MEF90Ctx,ierr);CHKERRQ(ierr)
      Call PetscBagSetName(bag,trim(name),"MEF90 Global properties object",ierr);CHKERRQ(ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix),ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,MEF90Ctx%verbose,default%verbose,'verbose','Verbosity: level',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterString(bag,MEF90Ctx%prefix,MEF90_MXSTRLEN,default%prefix,'prefix','prefix',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterEnum(bag,MEF90Ctx%timeInterpolation,MEF90TimeInterpolationList,default%timeInterpolation,'time_interpolation','Time: interpolation type',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,MEF90Ctx%timeMin,default%timeMin,'time_min','Time: min',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,MEF90Ctx%timeMax,default%timeMax,'time_max','Time: max',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,MEF90Ctx%timeNumStep,default%timeNumStep,'time_numstep','Time: number of time steps',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterEnum(bag,MEF90Ctx%fileFormat,MEF90FileFormatList,default%fileFormat,'file_format','I/O: file format.',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterEnum(bag,MEF90Ctx%fileMode,MEF90FileModeList,default%fileMode,'file_mode','I/O: file mode.',ierr);CHKERRQ(ierr)
   End Subroutine PetscBagRegisterMEF90Ctx

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
         Call EXGATM(MEF90Ctx%fileExoUnitIn,t,exoerr)
         MEF90Ctx%timeMin = t(1)
         MEF90Ctx%timeMax = t(MEF90Ctx%timeNumStep)
      End Select
      ierr = 0
   End Subroutine MEF90Ctx_GetTime
End Module m_MEF_Ctx
