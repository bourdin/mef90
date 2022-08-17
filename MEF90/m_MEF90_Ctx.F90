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
      PetscMPIInt                                     :: rank,numProcs
      Character(len=MEF90MXSTRLEN,kind=C_char)        :: geometryfile,resultfile
      PetscBag                                        :: GlobalOptionsBag      
      End Type MEF90Ctx_Type
   
   Type MEF90CtxGlobalOptions_Type
      PetscInt                                        :: verbose
      PetscBool                                       :: dryrun
      PetscEnum                                       :: timeInterpolation
      PetscReal                                       :: timeMin
      PetscReal                                       :: timeMax
      PetscInt                                        :: timeNumStep
      PetscInt                                        :: timeSkip
      PetscInt                                        :: timeNumCycle
      PetscEnum                                       :: elementFamily
      PetscInt                                        :: elementOrder
   End Type MEF90CtxGlobalOptions_Type
End Module m_MEF90_Ctx_Type

Module m_MEF90_Ctx
#include "petsc/finclude/petsc.h"
   Use, Intrinsic :: iso_c_binding
   Use m_MEF90_Parameters
   Use m_MEF90_Ctx_Type
   Use m_MEF90_Utils
   Use m_MEF90_Elements
   Implicit none

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
   Character(len=MEF90MXSTRLEN),dimension(7),protected  :: MEF90ScalingList
      
   Enum,bind(c)
      Enumerator  :: MEF90TimeInterpolation_linear = 0,  &
                     MEF90TimeInterpolation_Vcycle,     &
                     MEF90TimeInterpolation_quadratic,   &
                     MEF90TimeInterpolation_exo
   End Enum
   Character(len=MEF90MXSTRLEN),dimension(7),protected  :: MEF90TimeInterpolationList
   
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
!!!      2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90CtxInitialize_Private(ierr)
      PetscErrorCode,Intent(OUT)                   :: ierr
   
      Type(MEF90CtxGlobalOptions_Type)             :: MEF90CtxGlobalOptions
      character(len=1),pointer                     :: dummychar(:)
      PetscSizeT                                   :: sizeofchar
   
      PetscCall(PetscDataTypeGetSize(PETSC_CHAR,sizeofchar,ierr))
      sizeofMEF90CtxGlobalOptions = size(transfer(MEF90CtxGlobalOptions,dummychar))*sizeofchar

      MEF90ScalingList(1) = 'constant'
      MEF90ScalingList(2) = 'linear'
      MEF90ScalingList(3) = 'file'
      MEF90ScalingList(4) = 'null'
      MEF90ScalingList(5) = 'MEF90scaling'
      MEF90ScalingList(6) = '_MEF90Scaling'
      MEF90ScalingList(7) = ''
      
      MEF90TimeInterpolationList(1) = 'linear'
      MEF90TimeInterpolationList(2) = 'Vcycle'
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
!!!      2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine PetscBagGetDataMEF90CtxGlobalOptions(bag,data,ierr)
      PetscBag                                        :: bag
      Type(MEF90CtxGlobalOptions_Type),pointer        :: data
      PetscErrorCode,Intent(OUT)                      :: ierr
      
      PetscCall(PetscBagGetData(bag,data,ierr))
   End Subroutine PetscBagGetDataMEF90CtxGlobalOptions
   
#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90CtxGlobalOptions"
!!!
!!!  
!!!  PetscBagRegisterMEF90CtxGlobalOptions:
!!!  
!!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine PetscBagRegisterMEF90CtxGlobalOptions(bag,name,prefix,default,ierr)
      PetscBag                                        :: bag
      Character(len=*),intent(IN)                     :: prefix,name
      Type(MEF90CtxGlobalOptions_Type),intent(IN)     :: default
      PetscErrorCode,intent(OUT)                      :: ierr

      Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90CtxGlobalOptions
      
      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(bag,MEF90CtxGlobalOptions,ierr))
      PetscCall(PetscBagSetName(bag,trim(name),"MEF90 Global options object",ierr))
      PetscCall(PetscBagSetOptionsPrefix(bag,trim(prefix),ierr))
      PetscCall(PetscBagRegisterInt (bag,MEF90CtxGlobalOptions%verbose,default%verbose,'verbose','Verbosity: level',ierr))
      PetscCall(PetscBagRegisterBool(bag,MEF90CtxGlobalOptions%dryrun,default%dryrun,'dryrun','Dry run in order to validate the options file. Use in combination with -h to print help or -verbose 1 to check input deck',ierr))
      PetscCall(PetscBagRegisterEnum(bag,MEF90CtxGlobalOptions%timeInterpolation,MEF90TimeInterpolationList,default%timeInterpolation,'time_interpolation','Time: interpolation type',ierr))
      PetscCall(PetscBagRegisterReal(bag,MEF90CtxGlobalOptions%timeMin,default%timeMin,'time_min','Time: min',ierr))
      PetscCall(PetscBagRegisterReal(bag,MEF90CtxGlobalOptions%timeMax,default%timeMax,'time_max','Time: max',ierr))
      PetscCall(PetscBagRegisterInt (bag,MEF90CtxGlobalOptions%timeNumStep,default%timeNumStep,'time_numstep','Time: number of time steps',ierr))
      PetscCall(PetscBagRegisterInt (bag,MEF90CtxGlobalOptions%timeSkip,   default%timeSkip,'time_skip','Time: number of time steps',ierr))
      PetscCall(PetscBagRegisterInt(bag,MEF90CtxGlobalOptions%timenumCycle,default%timenumCycle,'time_numCycle','Time: number of cycles',ierr))
      PetscCall(PetscBagRegisterEnum(bag,MEF90CtxGlobalOptions%ElementFamily,MEF90ElementFamilyList,default%ElementFamily,'element_family','Element family (possibly overridden in application contexts)',ierr))
      PetscCall(PetscBagRegisterInt(bag,MEF90CtxGlobalOptions%ElementOrder,default%ElementOrder,'element_order','Element order (possibly overridden in application contexts)',ierr))
   End Subroutine PetscBagRegisterMEF90CtxGlobalOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90CtxCreate"
!!!
!!!  
!!!  MEF90CtxCreate:
!!!  
!!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90CtxCreate(comm,MEF90Ctx,default,ierr)
      MPI_Comm,Intent(IN)                          :: comm
      Type(MEF90Ctx_type),Intent(OUT)              :: MEF90Ctx
      Type(MEF90CtxGlobalOptions_Type),Intent(IN)  :: default
      PetscErrorCode,Intent(INOUT)                 :: ierr
      
      Type(MEF90CtxGlobalOptions_Type),pointer     :: GlobalOptions
      Character(len=MEF90MXSTRLEN)                 :: IOBuffer,tmpPrefix
      PetscBool                                    :: hasPrefix,hasGeometry,hasResult

#ifdef PETSC_USE_DEBUG
      Character(len=MPI_MAX_PROCESSOR_NAME)        :: procName
      Integer                                      :: procNameLength
#endif      
   
   MEF90Ctx%comm = comm
   PetscCallMPI(MPI_COMM_RANK(MEF90Ctx%comm,MEF90Ctx%rank,ierr))
   PetscCallMPI(MPI_COMM_SIZE(MEF90Ctx%comm,MEF90Ctx%numProcs,ierr))
   PetscCall(PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-prefix',tmpPrefix,hasPrefix,ierr))
   PetscCall(PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-geometry',MEF90Ctx%geometryFile,hasGeometry,ierr))
   PetscCall(PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-result',MEF90Ctx%resultFile,hasResult,ierr))
   If (.NOT. (hasPrefix .NEQV. hasGeometry)) Then
      PetscCall(PetscPrintf(comm,"prefix or geometry must be given (-prefix or -geometry) \n",ierr))
      PetscCall(PetscFinalize(ierr))
      STOP
   End If
   If (hasPrefix .AND. hasResult) Then
      PetscCall(PetscPrintf(comm,"-prefix and -result options incompatible.\n",ierr))
      PetscCall(PetscFinalize(ierr))
      STOP
      !SETERRQ(comm,PETSC_ERR_FILE_OPEN,"no file prefix given\n")
   End If
   If (hasPrefix) Then
      !!! Old style calling sequence: geometryFile is <prefix>.gen, resultFile is <prefix>_out.gen
      MEF90Ctx%geometryFile = trim(tmpPrefix)//'.gen'
      MEF90Ctx%resultFile = trim(MEF90FilePrefix(MEF90Ctx%geometryFile))//'_out.gen'
   Else
      If (.NOT. hasResult) Then
         MEF90Ctx%resultFile = trim(MEF90FilePrefix(MEF90Ctx%geometryFile))//'_out.'//trim(MEF90FileExtension(MEF90Ctx%geometryFile))
      End If
   End If
   PetscCall(PetscBagCreate(comm,sizeofMEF90CtxGlobalOptions,MEF90Ctx%GlobalOptionsBag,ierr))
   PetscCall(PetscBagRegisterMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,'MEF90Ctx',PETSC_NULL_CHARACTER,default,ierr))

   PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,GlobalOptions,ierr))

#ifdef PETSC_USE_DEBUG
   call MPI_Get_processor_name(procName,procNameLength,ierr)
   write(*,200) MEF90Ctx%rank,MEF90Ctx%numProcs,trim(procName)
   ! Write(IOBuffer,200) MEF90Ctx%rank,MEF90Ctx%numProcs,trim(procName)
   ! Call PetscSynchronizedPrintf(MEF90Ctx%comm,IOBuffer,ierr);CHKERRQ(ierr)
   ! Call PetscSynchronizedFlush(MEF90Ctx%comm,ierr);CHKERRQ(ierr)
#endif

   If (GlobalOptions%verbose > 0) Then
      Write(IOBuffer,*) 'MEF90 Global Context: \n'
      PetscCall(PetscPrintf(comm,IOBuffer,ierr))
      Write(IOBuffer,100) trim(MEF90Ctx%geometryFile)
      PetscCall(PetscPrintf(comm,IOBuffer,ierr))
      Write(IOBuffer,101) trim(MEF90Ctx%resultFile)
      PetscCall(PetscPrintf(comm,IOBuffer,ierr))
      Write(IOBuffer,102) trim(MEF90FilePrefix(MEF90Ctx%resultFile))//'.log'
      PetscCall(PetscPrintf(comm,IOBuffer,ierr))
      PetscCall(PetscBagView(MEF90Ctx%GlobalOptionsBag,PETSC_VIEWER_STDOUT_WORLD,ierr))
      PetscCall(PetscPrintf(comm,"\n",ierr))
   End If
100 Format('  geometry file:       ',(A),'\n')
101 Format('  result file:         ',(A),'\n')
102 Format('  log file:            ',(A),'\n')
200 format(' # Task ',I6,'/',I6,' running on processor ',A,'\n')
   End Subroutine MEF90CtxCreate

#undef __FUNCT__
#define __FUNCT__ "MEF90CtxDestroy"
!!!
!!!  
!!!  MEF90CtxDestroy:
!!!  
!!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90CtxDestroy(MEF90Ctx,ierr)
      Type(MEF90Ctx_Type),Intent(INOUT)               :: MEF90Ctx
      PetscErrorCode,Intent(OUT)                      :: ierr

      Type(MEF90CtxGlobalOptions_Type),pointer        :: GlobalOptions

      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,GlobalOptions,ierr))
      If (GlobalOptions%dryrun) Then
         PetscCall(PetscOptionsLeft(PETSC_NULL_OPTIONS,ierr))
      End If
      PetscCall(PetscBagDestroy(MEF90Ctx%GlobalOptionsBag,ierr))
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
      !Real                                            :: dummyR
      !Character(len=1)                                :: dummyS
      !Integer                                         :: exoerr
      Type(MEF90CtxGlobalOptions_Type),pointer        :: GlobalOptions
      Integer                                         :: j,CycleLength 

      i = 0 ! silence gfortran silly warning
      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,GlobalOptions,ierr))
      Select Case(GlobalOptions%timeInterpolation)
      Case (MEF90TimeInterpolation_linear)
         Allocate(t(GlobalOptions%timeNumStep))
         dt = 0.0_Kr
         If (GlobalOptions%timeNumStep > 1) Then
            dt = (GlobalOptions%timeMax - GlobalOptions%timeMin) / Real(GlobalOptions%timeNumStep-1.0_Kr)
         End If
         t = [ ( GlobalOptions%timeMin + Real(i) * dt , i = 0,GlobalOptions%timeNumStep-1 ) ]
         t(GlobalOptions%timeNumStep) = GlobalOptions%timeMax

      Case (MEF90TimeInterpolation_Vcycle)
         cycleLength = GlobalOptions%timeNumStep / GlobalOptions%timeNumCycle
         Allocate(t(GlobalOptions%timeNumStep))
         dt = 0.0_Kr
         If (GlobalOptions%timeNumStep > 1) Then
            dt = (GlobalOptions%timeMax - GlobalOptions%timeMin) * GlobalOptions%timeNumCycle / Real((GlobalOptions%timeNumStep-1))
         End If
         Do i = 1, GlobalOptions%timeNumCycle
            Do j = 1, cycleLength 
               t((i-1)*cycleLength+j) = min( GlobalOptions%timeMin + 2.0_Kr * Real(j-1) * dt, 2.0_Kr * GlobalOptions%timeMax - GlobalOptions%timeMin - 2.0_Kr * Real(j-1) * dt )
            End Do
         End Do
         t(cycleLength*GlobalOptions%timeNumCycle+1:GlobalOptions%timeNumStep) = t(cycleLength*GlobalOptions%timeNumCycle)

      Case (MEF90TimeInterpolation_quadratic)
         !!! Natural time scale for the heat equation
         Allocate(t(GlobalOptions%timeNumStep))
         dt = 0.0_Kr
         If (GlobalOptions%timeNumStep > 1) Then
            dt = (sqrt(GlobalOptions%timeMax) - sqrt(GlobalOptions%timeMin)) / Real(GlobalOptions%timeNumStep-1.0_Kr)
         End If
         t = [ ((sqrt(GlobalOptions%timeMin) + Real(i) * dt)**2, i = 0,GlobalOptions%timeNumStep-1) ]
         t(GlobalOptions%timeNumStep) = GlobalOptions%timeMax

      ! Case (MEF90TimeInterpolation_exo)
      !    Select case(GlobalOptions%FileFormat)
      !    Case (MEF90FileFormat_EXOSingle)
      !       If (MEF90Ctx%rank == 0) Then
      !          If (MEF90Ctx%fileExoUnit /= 0) Then
      !             Call EXINQ(MEF90Ctx%fileExoUnit,EXTIMS,GlobalOptions%timeNumStep,dummyR,dummyS,exoerr)
      !             Allocate(t(GlobalOptions%timeNumStep))
      !             Call EXGATM(MEF90Ctx%fileExoUnit,t,exoerr)
      !             PetscCallMPI(MPI_Bcast(GlobalOptions%timeNumStep,1,MPIU_INTEGER,0,MEF90Ctx%comm,ierr))
      !             PetscCallMPI(MPI_Bcast(t,GlobalOptions%timeNumStep,MPIU_SCALAR,0,MEF90Ctx%comm,ierr))
      !          Else
      !             Call PetscPrintf(PETSC_COMM_SELF,"EXO input file must be open prior to calling MEF90Ctx_GetTime\n",ierr);
      !             SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"EXO input file must be open prior to calling MEF90Ctx_GetTime\n")
      !          End If
      !       Else
      !          PetscCallMPI(MPI_Bcast(GlobalOptions%timeNumStep,1,MPIU_INTEGER,0,MEF90Ctx%comm,ierr))
      !          Allocate(t(GlobalOptions%timeNumStep))
      !          PetscCallMPI(MPI_Bcast(t,GlobalOptions%timeNumStep,MPIU_SCALAR,0,MEF90Ctx%comm,ierr))
      !       End If            
      !       GlobalOptions%timeMin = t(1)
      !       GlobalOptions%timeMax = t(GlobalOptions%timeNumStep)
      !    Case (MEF90FileFormat_EXOSplit)
      !       If (MEF90Ctx%fileExoUnit /= 0) Then
      !          Call EXINQ(MEF90Ctx%fileExoUnit,EXTIMS,GlobalOptions%timeNumStep,dummyR,dummyS,exoerr)
      !          Allocate(t(GlobalOptions%timeNumStep))
      !          Call EXGATM(MEF90Ctx%fileExoUnit,t,exoerr)
      !       Else
      !          Call PetscPrintf(PETSC_COMM_SELF,"EXO input file must be open prior to calling MEF90Ctx_GetTime\n",ierr);
      !          SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"EXO input file must be open prior to calling MEF90Ctx_GetTime\n")
      !       End If
      !    End Select
      End Select
      If ((GlobalOptions%verbose > 0) .AND. (MEF90Ctx%rank == 0)) Then
         PetscCall(PetscPrintf(PETSC_COMM_SELF,"Time values array:\n",ierr))
         PetscCall(PetscRealView(GlobalOptions%timeNumStep,t,PETSC_VIEWER_STDOUT_SELF,ierr))
         PetscCall(PetscPrintf(PETSC_COMM_SELF,"===\n",ierr))
      End If
   End Subroutine MEF90CtxGetTime
End Module m_MEF90_Ctx
