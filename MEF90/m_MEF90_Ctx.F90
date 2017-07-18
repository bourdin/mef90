Module m_MEF90_Ctx_Type
#include "petsc/finclude/petsc.h"
   Use m_MEF90_Parameters
   Use,Intrinsic :: iso_c_binding
   Implicit none
   Private  
   Public :: MEF90Ctx_Type
   Public :: MEF90CtxGlobalOptions_Type
   
   Type MEF90Ctx_Type
      MPI_Comm                                        :: comm
      MPI_Comm                                        :: IOcomm
      Integer                                         :: rank
      Character(len=MEF90_MXSTRLEN,kind=C_char)       :: inputmesh      
      Character(len=MEF90_MXSTRLEN,kind=C_char)       :: prefix
      Integer                                         :: fileExoUnit
      PetscBag                                        :: GlobalOptionsBag      
   End Type MEF90Ctx_Type
   
   Type MEF90CtxGlobalOptions_Type
      PetscInt                                        :: verbose
      PetscBool                                       :: dryrun
      !Character(len=MEF90_MXSTRLEN,kind=C_char)       :: prefix
      !!! There seems to be an incompatibility between PetscBagRegisterString and intel fortran 13.0
      !!! Moving prefix to the context itself instead of the bag
      PetscEnum                                       :: timeInterpolation
      PetscReal                                       :: timeMin
      PetscReal                                       :: timeMax
      PetscInt                                        :: timeNumStep
      PetscInt                                        :: timeSkip
      PetscReal                                       :: timeFrequency
      PetscEnum                                       :: fileFormat
   End Type MEF90CtxGlobalOptions_Type
End Module m_MEF90_Ctx_Type

Module m_MEF90_Ctx
#include "petsc/finclude/petsc.h"
   Use, Intrinsic :: iso_c_binding
   Use m_MEF90_Parameters
   Use m_MEF90_Ctx_Type
   Use m_MEF90_Utils
   Implicit none

   !Private  
   Public :: MEF90Ctx_Type
   Public :: MEF90CtxGlobalOptions_Type
   Public :: MEF90CtxInitialize_Private
   Public :: PetscBagGetDataMEF90CtxGlobalOptions
   Public :: MEF90CtxGetTime
   Public :: sizeofMEF90CtxGlobalOptions

   Private :: PetscBagGetData
   !!! Very important. PetscGetData must remain private to this module or others will not be able to declare their own interface 
   !!! for other derived type
      
   PetscSizeT,protected    :: sizeofMEF90CtxGlobalOptions

   Enum,bind(c)
      Enumerator ::  MEF90Scaling_CST=0,        &
                     MEF90Scaling_Linear,       &  
                     MEF90Scaling_File,         &  
                     MEF90Scaling_Null
   End Enum
   Character(len=MEF90_MXSTRLEN),dimension(7),protected  :: MEF90ScalingList
   
   Enum,bind(c)
      Enumerator  :: MEF90FileFormat_EXOSingle = 0,    &
                     MEF90FileFormat_EXOSplit 
   End Enum
   Character(len=MEF90_MXSTRLEN),dimension(5),protected  :: MEF90FileFormatList
   
   Enum,bind(c)
      Enumerator  :: MEF90TimeInterpolation_linear = 0,  &
                     MEF90TimeInterpolation_cycling,     &
                     MEF90TimeInterpolation_quadratic,   &
                     MEF90TimeInterpolation_exo
   End Enum
   Character(len=MEF90_MXSTRLEN),dimension(7),protected  :: MEF90TimeInterpolationList
   
   

   
   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF90_Ctx_Type
         PetscBag                                     :: bag
         Type(MEF90CtxGlobalOptions_Type),pointer     :: data
         PetscErrorCode,Intent(OUT)                   :: ierr
      End subroutine PetscBagGetData
   End interface

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90Ctx_InitializePrivate"
!!!
!!!  
!!!  MEF90Ctx_InitializePrivate:
!!!  
!!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90CtxInitialize_Private(ierr)
      PetscErrorCode,Intent(OUT)                   :: ierr
   
      Type(MEF90CtxGlobalOptions_Type)             :: MEF90CtxGlobalOptions
      character(len=1),pointer                     :: dummychar(:)
      PetscSizeT                                   :: sizeofchar
   
      Call PetscDataTypeGetSize(PETSC_CHAR,sizeofchar,ierr)
      sizeofMEF90CtxGlobalOptions = size(transfer(MEF90CtxGlobalOptions,dummychar))*sizeofchar

      MEF90ScalingList(1) = 'constant'
      MEF90ScalingList(2) = 'linear'
      MEF90ScalingList(3) = 'file'
      MEF90ScalingList(4) = 'null'
      MEF90ScalingList(5) = 'MEF90scaling'
      MEF90ScalingList(6) = '_MEF90Scaling'
      MEF90ScalingList(7) = ''
      
      MEF90FileFormatList(1) = 'Replace'
      MEF90FileFormatList(2) = 'MEF90FileFormat'
      MEF90FileFormatList(3) = '_MEF90FileFormat'
      MEF90FileFormatList(4) = ''

      MEF90FileFormatList(1) = 'EXOSingle'
      MEF90FileFormatList(2) = 'EXOSplit'
      MEF90FileFormatList(3) = 'MEF90FileFormat'
      MEF90FileFormatList(4) = '_MEF90FileFormatList'
      MEF90FileFormatList(5) = '' 
      
      MEF90TimeInterpolationList(1) = 'linear'
      MEF90TimeInterpolationList(2) = 'cycling'
      MEF90TimeInterpolationList(3) = 'quadratic'
      MEF90TimeInterpolationList(4) = 'exo'
      MEF90TimeInterpolationList(5) = 'MEF90TimeInterpolation'
      MEF90TimeInterpolationList(6) = '_MEF90TimeInterpolation'
      MEF90TimeInterpolationList(7) = ''
   End Subroutine MEF90CtxInitialize_Private
   
#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90CtxGlobalOptions"
!!!
!!!  
!!!  PetscBagGetDataMEF90CtxGlobalOptions:
!!!  
!!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine PetscBagGetDataMEF90CtxGlobalOptions(bag,data,ierr)
      PetscBag                                        :: bag
      Type(MEF90CtxGlobalOptions_Type),pointer        :: data
      PetscErrorCode,Intent(OUT)                      :: ierr
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataMEF90CtxGlobalOptions
   
#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90CtxGlobalOptions"
!!!
!!!  
!!!  PetscBagRegisterMEF90CtxGlobalOptions:
!!!  
!!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine PetscBagRegisterMEF90CtxGlobalOptions(bag,name,prefix,default,ierr)
      PetscBag                                        :: bag
      Character(len=*),intent(IN)                     :: prefix,name
      Type(MEF90CtxGlobalOptions_Type),intent(IN)     :: default
      PetscErrorCode,intent(OUT)                      :: ierr

      Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90CtxGlobalOptions
      
      Call PetscBagGetDataMEF90CtxGlobalOptions(bag,MEF90CtxGlobalOptions,ierr);CHKERRQ(ierr)
      Call PetscBagSetName(bag,trim(name),"MEF90 Global options object",ierr);CHKERRQ(ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix),ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,MEF90CtxGlobalOptions%verbose,default%verbose,'verbose','Verbosity: level',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterBool(bag,MEF90CtxGlobalOptions%dryrun,default%dryrun,'dryrun','Dry run in order to validate the options file. Use in combination with -h to print help or -verbose 1 to check input deck',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterEnum(bag,MEF90CtxGlobalOptions%timeInterpolation,MEF90TimeInterpolationList,default%timeInterpolation,'time_interpolation','Time: interpolation type',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,MEF90CtxGlobalOptions%timeMin,default%timeMin,'time_min','Time: min',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,MEF90CtxGlobalOptions%timeMax,default%timeMax,'time_max','Time: max',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,MEF90CtxGlobalOptions%timeNumStep,default%timeNumStep,'time_numstep','Time: number of time steps',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt (bag,MEF90CtxGlobalOptions%timeSkip,   default%timeSkip,'time_skip','Time: number of time steps',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterReal(bag,MEF90CtxGlobalOptions%timeFrequency,default%timeFrequency,'time_frequency','Time: frequency   (1+t) cos(2 k Pi t)',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterEnum(bag,MEF90CtxGlobalOptions%fileFormat,MEF90FileFormatList,default%fileFormat,'file_format','I/O: file format.',ierr);CHKERRQ(ierr)
   End Subroutine PetscBagRegisterMEF90CtxGlobalOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90CtxCreate"
!!!
!!!  
!!!  MEF90CtxCreate:
!!!  
!!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90CtxCreate(comm,MEF90Ctx,default,ierr)
         MPI_Comm,Intent(IN)                          :: comm
         Type(MEF90Ctx_type),Intent(OUT)              :: MEF90Ctx
         Type(MEF90CtxGlobalOptions_Type),Intent(IN)  :: default
         PetscInt,Intent(OUT)                         :: ierr
       
         Type(MEF90CtxGlobalOptions_Type),pointer     :: GlobalOptions
         Character(len=MEF90_MXSTRLEN)                :: IOBuffer
         PetscBool                                    :: flg
      
         MEF90Ctx%comm = comm
         Call MPI_COMM_RANK(MEF90Ctx%comm,MEF90Ctx%rank,ierr)
         Call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-inputmesh',MEF90Ctx%inputmesh,flg,ierr);CHKERRQ(ierr)
         If (.NOT. flg) Then
            Call PetscPrintf(PETSC_COMM_WORLD,"no input mesh given (-inputmesh)\n",ierr);CHKERRQ(ierr)
            Call PetscFinalize(ierr)
            STOP
            !SETERRQ(comm,PETSC_ERR_FILE_OPEN,"no file prefix given\n")
         End If
         Call MEF90GetFilePrefix(MEF90Ctx%inputmesh,MEF90Ctx%prefix)

         Call PetscBagCreate(comm,sizeofMEF90CtxGlobalOptions,MEF90Ctx%GlobalOptionsBag,ierr);CHKERRQ(ierr)
         Call PetscBagRegisterMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,'MEF90Ctx',PETSC_NULL_CHARACTER,default,ierr);CHKERRQ(ierr)

         Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)
         If (GlobalOptions%verbose > 0) Then
            Write(IOBuffer,*) 'MEF90 Global Context: \n'
            Call PetscPrintf(comm,IOBuffer,ierr);CHKERRQ(ierr)
            Write(IOBuffer,100) trim(MEF90Ctx%prefix)
            Call PetscPrintf(comm,IOBuffer,ierr);CHKERRQ(ierr)
            Write(IOBuffer,102) trim(MEF90Ctx%prefix)//'.log'
            Call PetscPrintf(comm,IOBuffer,ierr);CHKERRQ(ierr)
            Call PetscBagView(MEF90Ctx%GlobalOptionsBag,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
            Call PetscPrintf(comm,"\n",ierr);CHKERRQ(ierr)
         End If
         If (GlobalOptions%fileformat == MEF90FileFormat_EXOSingle) Then
            MEF90Ctx%IOComm = MEF90Ctx%comm
         Else
            MEF90Ctx%IOComm = PETSC_COMM_SELF
         End If   
         MEF90Ctx%fileexounit = 0
      100 Format('  prefix:       ',(A),'\n')
      102 Format('  log file:     ',(A),'\n')
   End Subroutine MEF90CtxCreate

#undef __FUNCT__
#define __FUNCT__ "MEF90CtxDestroy"
!!!
!!!  
!!!  MEF90CtxDestroy:
!!!  
!!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90CtxDestroy(MEF90Ctx,ierr)
      Type(MEF90Ctx_Type),Intent(INOUT)               :: MEF90Ctx
      PetscErrorCode,Intent(OUT)                      :: ierr

      Type(MEF90CtxGlobalOptions_Type),pointer        :: GlobalOptions

      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)
      If (GlobalOptions%dryrun) Then
         Call PetscOptionsLeft(PETSC_NULL_OPTIONS,ierr);CHKERRQ(ierr)
      End If
      Call PetscBagDestroy(MEF90Ctx%GlobalOptionsBag,ierr);CHKERRQ(ierr)
   End Subroutine MEF90CtxDestroy

#undef __FUNCT__
#define __FUNCT__ "MEF90CtxGetTime"
!!!
!!!  
!!!  MEF90CtxGetTime:
!!!  
!!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90CtxGetTime(MEF90Ctx,t,ierr)
      Type(MEF90Ctx_Type),Intent(INOUT)               :: MEF90Ctx
      PetscReal,Dimension(:),Pointer                  :: t
      PetscErrorCode,Intent(OUT)                      :: ierr
      
      PetscReal                                       :: dt
      Integer                                         :: i
      Real                                            :: dummyR
      Character(len=1)                                :: dummyS
      Integer                                         :: exoerr
      Type(MEF90CtxGlobalOptions_Type),pointer        :: GlobalOptions

      i = 0 ! silence gfortran silly warning
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)
      Select Case(GlobalOptions%timeInterpolation)
      Case (MEF90TimeInterpolation_linear)
         Allocate(t(GlobalOptions%timeNumStep))
         dt = 0.0_Kr
         If (GlobalOptions%timeNumStep > 1) Then
            dt = (GlobalOptions%timeMax - GlobalOptions%timeMin) / Real(GlobalOptions%timeNumStep-1.0_Kr)
         End If
         t = [ ( GlobalOptions%timeMin + Real(i) * dt , i = 0,GlobalOptions%timeNumStep-1 ) ]
         t(GlobalOptions%timeNumStep) = GlobalOptions%timeMax

      Case (MEF90TimeInterpolation_cycling)
         Allocate(t(GlobalOptions%timeNumStep))
         dt = 0.0_Kr
         If (GlobalOptions%timeNumStep > 1) Then
            dt = (GlobalOptions%timeMax - GlobalOptions%timeMin) / Real(GlobalOptions%timeNumStep-1.0_Kr)
         End If
         t = [ ( (1.0_Kr + (GlobalOptions%timeMin + Real(i) * dt )) * (1.0_Kr - cos( 2.0_Kr*PETSC_PI*GlobalOptions%timeFrequency * (GlobalOptions%timeMin + Real(i) * dt ))), i = 0,GlobalOptions%timeNumStep-1) ]
         t(GlobalOptions%timeNumStep) = ( (1.0_Kr + (GlobalOptions%timeMax ))*( 1.0_Kr - cos(  2.0_Kr*PETSC_PI*GlobalOptions%timeFrequency * GlobalOptions%timeMax ) ))

      Case (MEF90TimeInterpolation_quadratic)
         !!! Natural time scale for the heat equation
         Allocate(t(GlobalOptions%timeNumStep))
         dt = 0.0_Kr
         If (GlobalOptions%timeNumStep > 1) Then
            dt = (sqrt(GlobalOptions%timeMax) - sqrt(GlobalOptions%timeMin)) / Real(GlobalOptions%timeNumStep-1.0_Kr)
         End If
         t = [ ((sqrt(GlobalOptions%timeMin) + Real(i) * dt)**2, i = 0,GlobalOptions%timeNumStep-1) ]
         t(GlobalOptions%timeNumStep) = GlobalOptions%timeMax

      Case (MEF90TimeInterpolation_exo)
         Select case(GlobalOptions%FileFormat)
         Case (MEF90FileFormat_EXOSingle)
            If (MEF90Ctx%rank == 0) Then
               If (MEF90Ctx%fileExoUnit /= 0) Then
                  Call EXINQ(MEF90Ctx%fileExoUnit,EXTIMS,GlobalOptions%timeNumStep,dummyR,dummyS,exoerr)
                  Allocate(t(GlobalOptions%timeNumStep))
                  Call EXGATM(MEF90Ctx%fileExoUnit,t,exoerr)
                  Call MPI_Bcast(GlobalOptions%timeNumStep,1,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)
                  Call MPI_Bcast(t,GlobalOptions%timeNumStep,MPIU_SCALAR,0,MEF90Ctx%comm,ierr)
               Else
                  Call PetscPrintf(PETSC_COMM_SELF,"EXO input file must be open prior to calling MEF90Ctx_GetTime\n",ierr);
                  SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"EXO input file must be open prior to calling MEF90Ctx_GetTime\n")
               End If
            Else
               Call MPI_Bcast(GlobalOptions%timeNumStep,1,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)
               Allocate(t(GlobalOptions%timeNumStep))
               Call MPI_Bcast(t,GlobalOptions%timeNumStep,MPIU_SCALAR,0,MEF90Ctx%comm,ierr)
            End If            
            GlobalOptions%timeMin = t(1)
            GlobalOptions%timeMax = t(GlobalOptions%timeNumStep)
         Case (MEF90FileFormat_EXOSplit)
            If (MEF90Ctx%fileExoUnit /= 0) Then
               Call EXINQ(MEF90Ctx%fileExoUnit,EXTIMS,GlobalOptions%timeNumStep,dummyR,dummyS,exoerr)
               Allocate(t(GlobalOptions%timeNumStep))
               Call EXGATM(MEF90Ctx%fileExoUnit,t,exoerr)
            Else
               Call PetscPrintf(PETSC_COMM_SELF,"EXO input file must be open prior to calling MEF90Ctx_GetTime\n",ierr);
               SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"EXO input file must be open prior to calling MEF90Ctx_GetTime\n")
            End If
         End Select
      End Select
      If ((GlobalOptions%verbose > 0) .AND. (MEF90Ctx%rank == 0)) Then
         Call PetscPrintf(PETSC_COMM_SELF,"Time values array:\n",ierr);CHKERRQ(ierr)  
         Call PetscRealView(GlobalOptions%timeNumStep,t,PETSC_VIEWER_STDOUT_SELF,ierr);CHKERRQ(ierr)
         Call PetscPrintf(PETSC_COMM_SELF,"===\n",ierr);CHKERRQ(ierr)  
      End If
   End Subroutine MEF90CtxGetTime
End Module m_MEF90_Ctx
