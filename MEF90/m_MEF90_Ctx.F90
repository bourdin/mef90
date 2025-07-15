module m_MEF90_Ctx_Type
#include "petsc/finclude/petsc.h"
   use petscbag
   use m_MEF90_Parameters
   use m_MEF90_Utils
   use m_MEF90_LinAlg
   use m_MEF90_Elements
   use, intrinsic :: iso_c_binding
   implicit none(type, external)
   private
   public :: MEF90Ctx_Type
   public :: MEF90CtxGlobalOptions_Type

   type MEF90Ctx_Type
      MPI_Comm                                        :: comm
      PetscMPIInt                                     :: rank, numProcs
      character(len=MEF90MXSTRLEN, kind=c_char)        :: geometryfile, resultfile
      type(tPetscBag)                                 :: GlobalOptionsBag
      type(tPetscViewer)                              :: resultViewer
   end type MEF90Ctx_Type

   type MEF90CtxGlobalOptions_Type
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
   end type MEF90CtxGlobalOptions_Type
end module m_MEF90_Ctx_Type

module m_MEF90_Ctx
#include "petsc/finclude/petsc.h"
   use, intrinsic :: iso_c_binding
   use petscbag
   use m_MEF90_Parameters
   use m_MEF90_Utils
   use m_MEF90_LinAlg
   use m_MEF90_Elements
   use m_MEF90_Ctx_Type

   implicit none(type)

   public :: MEF90Ctx_Type
   public :: MEF90CtxGlobalOptions_Type
   public :: MEF90CtxInitialize_Private
   public :: PetscBagGetDataMEF90CtxGlobalOptions
   public :: MEF90CtxGetTime
   public :: sizeofMEF90CtxGlobalOptions

   private :: PetscBagGetData
   !!! Very important. PetscGetData must remain private to this module or others will not be able to declare their own interface
   !!! for other derived type

   PetscSizeT, protected    :: sizeofMEF90CtxGlobalOptions

   enum, bind(c)
      enumerator ::  MEF90Scaling_CST = 0, &
         MEF90Scaling_Linear, &
         MEF90Scaling_File, &
         MEF90Scaling_Expr, &
         MEF90Scaling_Null
   end enum
   character(len=MEF90MXSTRLEN), dimension(8), protected  :: MEF90ScalingList

   enum, bind(c)
      enumerator  :: MEF90TimeInterpolation_linear = 0, &
         MEF90TimeInterpolation_Vcycle, &
         MEF90TimeInterpolation_quadratic, &
         MEF90TimeInterpolation_exo
   end enum
   character(len=MEF90MXSTRLEN), dimension(7), protected  :: MEF90TimeInterpolationList

   interface PetscBagGetData
      subroutine PetscBagGetData(bag, data, ierr)
         use m_MEF90_Ctx_Type
         use petscbag
         implicit none(type)
         type(tPetscBag), intent(IN)                           :: bag
         type(MEF90CtxGlobalOptions_Type), pointer, intent(OUT) :: data
         PetscErrorCode, intent(INOUT)                         :: ierr
      end subroutine PetscBagGetData
   end interface PetscBagGetData

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90Ctx_InitializePrivate"
!!!
!!!
!!!  MEF90Ctx_InitializePrivate:
!!!
!!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90CtxInitialize_Private(ierr)
      PetscErrorCode, intent(OUT)                   :: ierr

      type(MEF90CtxGlobalOptions_Type)             :: MEF90CtxGlobalOptions
      character(len=1), pointer                     :: dummychar(:)
      PetscSizeT                                   :: sizeofchar

      PetscCall(PetscDataTypeGetSize(PETSC_CHAR, sizeofchar, ierr))
      sizeofMEF90CtxGlobalOptions = size(transfer(MEF90CtxGlobalOptions, dummychar)) * sizeofchar

      MEF90ScalingList(1) = 'constant'
      MEF90ScalingList(2) = 'linear'
      MEF90ScalingList(3) = 'file'
      MEF90ScalingList(4) = 'expression'
      MEF90ScalingList(5) = 'null'
      MEF90ScalingList(6) = 'MEF90scaling'
      MEF90ScalingList(7) = '_MEF90Scaling'
      MEF90ScalingList(8) = ''

      MEF90TimeInterpolationList(1) = 'linear'
      MEF90TimeInterpolationList(2) = 'Vcycle'
      MEF90TimeInterpolationList(3) = 'quadratic'
      MEF90TimeInterpolationList(4) = 'exo'
      MEF90TimeInterpolationList(5) = 'MEF90TimeInterpolation'
      MEF90TimeInterpolationList(6) = '_MEF90TimeInterpolation'
      MEF90TimeInterpolationList(7) = ''
   end subroutine MEF90CtxInitialize_Private

#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90CtxGlobalOptions"
!!!
!!!
!!!  PetscBagGetDataMEF90CtxGlobalOptions:
!!!
!!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine PetscBagGetDataMEF90CtxGlobalOptions(bag, data, ierr)
      type(tPetscBag), intent(IN)                           :: bag
      type(MEF90CtxGlobalOptions_Type), pointer, intent(OUT) :: data
      PetscErrorCode, intent(INOUT)                         :: ierr

      PetscCall(PetscBagGetData(bag, data, ierr))
   end subroutine PetscBagGetDataMEF90CtxGlobalOptions

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90CtxGlobalOptions"
!!!
!!!
!!!  PetscBagRegisterMEF90CtxGlobalOptions:
!!!
!!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine PetscBagRegisterMEF90CtxGlobalOptions(bag, name, prefix, default, ierr)
      type(tPetscBag), intent(IN)                      :: bag
      character(len=*), intent(IN)                     :: prefix, name
      type(MEF90CtxGlobalOptions_Type), intent(IN)     :: default
      PetscErrorCode, intent(INOUT)                    :: ierr

      type(MEF90CtxGlobalOptions_Type), pointer        :: MEF90CtxGlobalOptions

      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(bag, MEF90CtxGlobalOptions, ierr))
      PetscCall(PetscBagSetName(bag, trim(name), "MEF90 Global options object", ierr))
      PetscCall(PetscBagSetOptionsPrefix(bag, trim(prefix), ierr))
      PetscCall(PetscBagRegisterInt(bag, MEF90CtxGlobalOptions % verbose, default % verbose, 'verbose', 'Verbosity: level', ierr))
      PetscCall(PetscBagRegisterBool(bag, MEF90CtxGlobalOptions % dryrun, default % dryrun, 'dryrun', 'Dry run in order to validate the options file. Use in combination with -h to print help or -verbose 1 to check input deck', ierr))
      PetscCall(PetscBagRegisterEnum(bag, MEF90CtxGlobalOptions % timeInterpolation, MEF90TimeInterpolationList, default % timeInterpolation, 'time_interpolation', 'Time: interpolation type', ierr))
      PetscCall(PetscBagRegisterReal(bag, MEF90CtxGlobalOptions % timeMin, default % timeMin, 'time_min', 'Time: min', ierr))
      PetscCall(PetscBagRegisterReal(bag, MEF90CtxGlobalOptions % timeMax, default % timeMax, 'time_max', 'Time: max', ierr))
      PetscCall(PetscBagRegisterInt(bag, MEF90CtxGlobalOptions % timeNumStep, default % timeNumStep, 'time_numstep', 'Time: number of time steps', ierr))
      PetscCall(PetscBagRegisterInt(bag, MEF90CtxGlobalOptions % timeSkip, default % timeSkip, 'time_skip', 'Time: number of time steps', ierr))
      PetscCall(PetscBagRegisterInt(bag, MEF90CtxGlobalOptions % timenumCycle, default % timenumCycle, 'time_numCycle', 'Time: number of cycles', ierr))
      PetscCall(PetscBagRegisterEnum(bag, MEF90CtxGlobalOptions % ElementFamily, MEF90ElementFamilyList, default % ElementFamily, 'element_family', 'Element family (possibly overridden in application contexts)', ierr))
      PetscCall(PetscBagRegisterInt(bag, MEF90CtxGlobalOptions % ElementOrder, default % ElementOrder, 'element_order', 'Element order (possibly overridden in application contexts)', ierr))
   end subroutine PetscBagRegisterMEF90CtxGlobalOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90CtxCreate"
!!!
!!!
!!!  MEF90CtxCreate:
!!!
!!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90CtxCreate(comm, MEF90Ctx, default, ierr)
      MPI_Comm, intent(IN)                          :: comm
      type(MEF90Ctx_type), intent(OUT)              :: MEF90Ctx
      type(MEF90CtxGlobalOptions_Type), intent(IN)  :: default
      PetscErrorCode, intent(INOUT)                 :: ierr

      type(MEF90CtxGlobalOptions_Type), pointer     :: GlobalOptions
      character(len=MEF90MXSTRLEN)                 :: IOBuffer, tmpPrefix
      PetscBool                                    :: hasPrefix, hasGeometry, hasResult

#ifdef PETSC_USE_DEBUG
      character(len=MPI_MAX_PROCESSOR_NAME)        :: procName
      integer                                      :: procNameLength
#endif

      MEF90Ctx % comm = comm
      PetscCallMPI(MPI_COMM_RANK(MEF90Ctx % comm, MEF90Ctx % rank, ierr))
      PetscCallMPI(MPI_COMM_SIZE(MEF90Ctx % comm, MEF90Ctx % numProcs, ierr))
      PetscCall(PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-prefix', tmpPrefix, hasPrefix, ierr))
      PetscCall(PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-geometry', MEF90Ctx % geometryFile, hasGeometry, ierr))
      PetscCall(PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-result', MEF90Ctx % resultFile, hasResult, ierr))
      if (.not. (hasPrefix .neqv. hasGeometry)) then
         PetscCall(PetscPrintf(comm, "prefix or geometry must be given (-prefix or -geometry) \n", ierr))
         PetscCall(PetscFinalize(ierr))
         stop
      end if
      if (hasPrefix .and. hasResult) then
         PetscCall(PetscPrintf(comm, "-prefix and -result options incompatible.\n", ierr))
         PetscCall(PetscFinalize(ierr))
         stop
         !SETERRQ(comm,PETSC_ERR_FILE_OPEN,"no file prefix given\n")
      end if
      if (hasPrefix) then
         !!! Old style calling sequence: geometryFile is <prefix>.gen, resultFile is <prefix>_out.gen
         MEF90Ctx % geometryFile = trim(tmpPrefix)//'.gen'
         MEF90Ctx % resultFile = trim(MEF90FilePrefix(MEF90Ctx % geometryFile))//'_out.gen'
      else
         if (.not. hasResult) then
            MEF90Ctx % resultFile = trim(MEF90FilePrefix(MEF90Ctx % geometryFile))//'_out.gen'
         end if
      end if
      PetscCall(PetscBagCreate(comm, sizeofMEF90CtxGlobalOptions, MEF90Ctx % GlobalOptionsBag, ierr))
      PetscCall(PetscBagRegisterMEF90CtxGlobalOptions(MEF90Ctx % GlobalOptionsBag, 'MEF90Ctx', PETSC_NULL_CHARACTER, default, ierr))

      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx % GlobalOptionsBag, GlobalOptions, ierr))

#ifdef PETSC_USE_DEBUG
      PetscCallMPI(MPI_Get_processor_name(procName, procNameLength, ierr))
      write (IOBuffer, "(' # Task ',I6,'/',I6,' running on processor ',A,'\n')") MEF90Ctx % rank, MEF90Ctx % numProcs, trim(procName)
      PetscCall(PetscSynchronizedPrintf(MEF90Ctx % comm, IOBuffer, ierr))
      PetscCall(PetscSynchronizedFlush(MEF90Ctx % comm, PETSC_STDOUT, ierr))
#endif

      if (GlobalOptions % verbose > 0) then
         write (IOBuffer, *) 'MEF90 Global Context: \n'
         PetscCall(PetscPrintf(comm, IOBuffer, ierr))
         write (IOBuffer, "('  geometry file:       ',(A),'\n')") trim(MEF90Ctx % geometryFile)
         PetscCall(PetscPrintf(comm, IOBuffer, ierr))
         write (IOBuffer, "('  result file:         ',(A),'\n')") trim(MEF90Ctx % resultFile)
         PetscCall(PetscPrintf(comm, IOBuffer, ierr))
         write (IOBuffer, "('  log file:            ',(A),'\n')") trim(MEF90FilePrefix(MEF90Ctx % resultFile))//'.log'
         PetscCall(PetscPrintf(comm, IOBuffer, ierr))
         PetscCall(PetscBagView(MEF90Ctx % GlobalOptionsBag, PETSC_VIEWER_STDOUT_WORLD, ierr))
         PetscCall(PetscPrintf(comm, "\n", ierr))
      end if

      !!! Not sure if this should be there, but PETSc's gmsh reader defaults to ignoring vertex sets, which we defintely don't want...
      if (MEF90FileExtension(MEF90Ctx % geometryfile) == 'msh') then
         PetscCallA(PetscOptionsInsertString(PETSC_NULL_OPTIONS, "-dm_plex_gmsh_mark_vertices", ierr))
      end if
   end subroutine MEF90CtxCreate

#undef __FUNCT__
#define __FUNCT__ "MEF90CtxDestroy"
!!!
!!!
!!!  MEF90CtxDestroy:
!!!
!!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90CtxDestroy(MEF90Ctx, ierr)
      type(MEF90Ctx_Type), intent(INOUT)               :: MEF90Ctx
      PetscErrorCode, intent(OUT)                      :: ierr

      type(MEF90CtxGlobalOptions_Type), pointer        :: GlobalOptions

      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx % GlobalOptionsBag, GlobalOptions, ierr))
      if (GlobalOptions % dryrun) then
         PetscCall(PetscOptionsLeft(PETSC_NULL_OPTIONS, ierr))
      end if
      PetscCall(PetscBagDestroy(MEF90Ctx % GlobalOptionsBag, ierr))
      if (.not. PetscObjectIsNull(MEF90Ctx % resultViewer)) then
         PetscCall(PetscViewerDestroy(MEF90Ctx % resultViewer, ierr))
      end if
   end subroutine MEF90CtxDestroy

#undef __FUNCT__
#define __FUNCT__ "MEF90CtxGetTime"
!!!
!!!
!!!  MEF90CtxGetTime:
!!!
!!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!

   subroutine MEF90CtxGetTime(MEF90Ctx, t, ierr)
      type(MEF90Ctx_Type), intent(INOUT)               :: MEF90Ctx
      PetscReal, dimension(:), pointer                  :: t
      PetscErrorCode, intent(OUT)                      :: ierr

      PetscReal                                       :: dt
      integer                                         :: i
      real                                            :: dummyR
      character(len=1)                                :: dummyS
      PetscErrorCode                                  :: exoErr
      integer                                         :: exoUnit
      type(MEF90CtxGlobalOptions_Type), pointer        :: GlobalOptions
      integer                                         :: j, CycleLength
      character(len=MEF90MXSTRLEN)                    :: IOBuffer

      i = 0 ! silence gfortran silly warning
      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx % GlobalOptionsBag, GlobalOptions, ierr))
      select case (GlobalOptions % timeInterpolation)
      case (MEF90TimeInterpolation_linear)
         allocate (t(GlobalOptions % timeNumStep))
         dt = 0.0_kr
         if (GlobalOptions % timeNumStep > 1) then
            dt = (GlobalOptions % timeMax - GlobalOptions % timeMin) / real(GlobalOptions % timeNumStep - 1.0_kr)
         end if
         t = [(GlobalOptions % timeMin + real(i) * dt, i=0, GlobalOptions % timeNumStep - 1)]
         t(GlobalOptions % timeNumStep) = GlobalOptions % timeMax

      case (MEF90TimeInterpolation_Vcycle)
         cycleLength = GlobalOptions % timeNumStep / GlobalOptions % timeNumCycle
         allocate (t(GlobalOptions % timeNumStep))
         dt = 0.0_kr
         if (GlobalOptions % timeNumStep > 1) then
            dt = (GlobalOptions % timeMax - GlobalOptions % timeMin) * GlobalOptions % timeNumCycle / real((GlobalOptions % timeNumStep - 1))
         end if
         do i = 1, GlobalOptions % timeNumCycle
            do j = 1, cycleLength
               t((i - 1) * cycleLength + j) = min(GlobalOptions % timeMin + 2.0_kr * real(j - 1) * dt, 2.0_kr * GlobalOptions % timeMax - GlobalOptions % timeMin - 2.0_kr * real(j - 1) * dt)
            end do
         end do
         t(cycleLength * GlobalOptions % timeNumCycle + 1:GlobalOptions % timeNumStep) = t(cycleLength * GlobalOptions % timeNumCycle)

      case (MEF90TimeInterpolation_quadratic)
         !!! Natural time scale for the heat equation
         allocate (t(GlobalOptions % timeNumStep))
         dt = 0.0_kr
         if (GlobalOptions % timeNumStep > 1) then
            dt = (sqrt(GlobalOptions % timeMax) - sqrt(GlobalOptions % timeMin)) / real(GlobalOptions % timeNumStep - 1.0_kr)
         end if
         t = [((sqrt(GlobalOptions % timeMin) + real(i) * dt)**2, i=0, GlobalOptions % timeNumStep - 1)]
         t(GlobalOptions % timeNumStep) = GlobalOptions % timeMax

      case (MEF90TimeInterpolation_exo)
         PetscCall(PetscViewerExodusIIGetId(MEF90Ctx % resultViewer, exoUnit, ierr))
         if (.not. PetscObjectIsNull(MEF90Ctx % resultViewer)) then
            call EXINQ(exoUnit, EXTIMS, GlobalOptions % timeNumStep, dummyR, dummyS, exoErr)
            allocate (t(GlobalOptions % timeNumStep))
            call EXGATM(exoUnit, t, exoErr)
         else
            write (IOBuffer, "(A,'EXO input file must be open prior to calling MEF90Ctx_GetTime\n')") __FUNCT__
            SETERRQ(MEF90Ctx % Comm, PETSC_ERR_FILE_OPEN, IOBuffer)
         end if
      case Default
         write (IOBuffer, "(A,'Unimplemented time interpolation: ',I0,'\n')") __FUNCT__, GlobalOptions % timeInterpolation
         SETERRQ(MEF90Ctx % Comm, PETSC_ERR_ARG_OUTOFRANGE, IOBuffer)
      end select
      if ((GlobalOptions % verbose > 0) .and. (MEF90Ctx % rank == 0)) then
         PetscCall(PetscPrintf(MEF90Ctx % Comm, "Time values array:\n", ierr))
         PetscCall(PetscRealView(GlobalOptions % timeNumStep, t, PETSC_VIEWER_STDOUT_SELF, ierr))
         PetscCall(PetscPrintf(MEF90Ctx % Comm, "===\n", ierr))
      end if
   end subroutine MEF90CtxGetTime
end module m_MEF90_Ctx
