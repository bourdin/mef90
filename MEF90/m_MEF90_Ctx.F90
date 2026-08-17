#include "../MEF90/mef90.inc"
module m_MEF90_Ctx
#include "petsc/finclude/petsc.h"
   use, intrinsic :: iso_c_binding
   use petscsys
   use m_MEF90_Parameters
   use m_MEF90_Utils
   use m_MEF90_LinAlg
   use m_MEF90_Elements
   use m_MEF90_BaseClass

   implicit none(type)

   public :: MEF90Ctx_Type
   public :: MEF90CtxGlobalOptions_Type
   public :: MEF90CtxGetTime
   public :: MEF90CtxGlobalOptionsSetFromOptions

!!! author: Blaise Bourdin (2012-2014, bourdin@lsu.edu)
!!! author: Blaise Bourdin (2026, bourdin@mcmaster.ca)
!!!
!!!  MEF90Ctx_Type: The class holding the MEF90-wide state and options
!!!

   enum, bind(c)
      enumerator ::  MEF90Scaling_CST = 0, &
         MEF90Scaling_Linear, &
         MEF90Scaling_File, &
         MEF90Scaling_Expr, &
         MEF90Scaling_Null
   end enum
   character(len=MEF90MXSTRLEN), dimension(8), protected  :: MEF90ScalingList = [character(len=MEF90MXSTRLEN) :: &
      'constant', &
      'linear', &
      'file', &
      'expression', &
      'null', &
      'MEF90scaling', &
      '_MEF90Scaling', &
      '']

   enum, bind(c)
      enumerator  :: MEF90TimeInterpolation_linear = 0, &
         MEF90TimeInterpolation_Vcycle, &
         MEF90TimeInterpolation_quadratic, &
         MEF90TimeInterpolation_exo
   end enum
   character(len=MEF90MXSTRLEN), dimension(7), protected  :: MEF90TimeInterpolationList = [character(len=MEF90MXSTRLEN) :: &
      'linear', &
      'Vcycle', &
      'quadratic', &
      'exo', &
      'MEF90TimeInterpolation', &
      '_MEF90TimeInterpolation', &
      '']

!!!
!!!  MEF90CtxGlobalOptions_Type: the MEF90-wide options.
!!!                              The values given here are the defaults used by setFromOptions
!!!
   type MEF90CtxGlobalOptions_Type
      PetscInt                                        :: verbose = 0_ki
      PetscBool                                       :: dryrun = PETSC_FALSE
      PetscEnum                                       :: timeInterpolation = MEF90TimeInterpolation_linear
      PetscReal                                       :: timeMin = 0.0_kr
      PetscReal                                       :: timeMax = 1.0_kr
      PetscInt                                        :: timeNumStep = 11_ki
      PetscInt                                        :: timeSkip = 0_ki
      PetscInt                                        :: timeNumCycle = 1_ki
      PetscEnum                                       :: elementFamily = MEF90ElementFamilyLagrange
      PetscInt                                        :: elementOrder = 1_ki
   end type MEF90CtxGlobalOptions_Type

   type, extends(MEF90Object) :: MEF90Ctx_Type
      PetscMPIInt                                     :: rank, numProcs
      character(len=MEF90MXSTRLEN, kind=c_char)       :: geometryfile, resultfile
      !!! The options are not stored: they are read from the options database where they are
      !!! needed, with MEF90CtxGlobalOptionsSetFromOptions. Its options argument is intent(inout),
      !!! so a caller that wants a default of its own sets it on the record before the call.
      type(tPetscViewer)                              :: resultViewer
   contains
      procedure, pass(self) :: setFromOptions => MEF90CtxSetFromOptions
      procedure, pass(self) :: view_internal => MEF90CtxView
   end type MEF90Ctx_Type

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90CtxCreate"
!!! author: Blaise Bourdin (2012-2014, bourdin@lsu.edu)
!!! author: Blaise Bourdin (2022, bourdin@mcmaster.ca)
!!!
!!!  MEF90CtxCreate:
!!!

   subroutine MEF90CtxCreate(comm, MEF90Ctx, prefix, ierr)
      MPIU_Comm, intent(IN)                         :: comm
      type(MEF90Ctx_type), intent(OUT)              :: MEF90Ctx
      character(len=*), intent(IN)                  :: prefix
      PetscErrorCode, intent(INOUT)                 :: ierr

      character(len=MEF90MXSTRLEN)                  :: IOBuffer, tmpPrefix
      PetscBool                                     :: hasPrefix, hasGeometry, hasResult

#ifdef PETSC_USE_DEBUG
      character(len=MPI_MAX_PROCESSOR_NAME)         :: procName
      integer                                       :: procNameLength
#endif

      MEF90Ctx%comm = comm
      MEF90Ctx%prefix = prefix
      MEF90Ctx%name = trim(prefix)//"MEF90Ctx"
      PetscCallMPI(MPI_COMM_RANK(MEF90Ctx%comm, MEF90Ctx%rank, ierr))
      PetscCallMPI(MPI_COMM_SIZE(MEF90Ctx%comm, MEF90Ctx%numProcs, ierr))
      PetscCall(PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-prefix', tmpPrefix, hasPrefix, ierr))
      PetscCall(PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-geometry', MEF90Ctx%geometryFile, hasGeometry, ierr))
      PetscCall(PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-result', MEF90Ctx%resultFile, hasResult, ierr))
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
         !! Old style calling sequence: geometryFile is <prefix>.gen, resultFile is <prefix>_out.gen
         MEF90Ctx%geometryFile = trim(tmpPrefix)//'.gen'
         MEF90Ctx%resultFile = trim(MEF90FilePrefix(MEF90Ctx%geometryFile))//'_out.gen'
      else
         if (.not. hasResult) then
            MEF90Ctx%resultFile = trim(MEF90FilePrefix(MEF90Ctx%geometryFile))//'_out.gen'
         end if
      end if
#ifdef PETSC_USE_DEBUG
      PetscCallMPI(MPI_Get_processor_name(procName, procNameLength, ierr))
      write (IOBuffer, "(' # Task ',I6,'/',I6,' running on processor ',A,'\n')") MEF90Ctx%rank, MEF90Ctx%numProcs, trim(procName)
      PetscCall(PetscSynchronizedPrintf(MEF90Ctx%comm, IOBuffer, ierr))
      PetscCall(PetscSynchronizedFlush(MEF90Ctx%comm, PETSC_STDOUT, ierr))
#endif

      !! Not sure if this should be there, but PETSc's gmsh reader defaults to ignoring vertex sets, which we defintely don't want...
      if (MEF90FileExtension(MEF90Ctx%geometryfile) == 'msh') then
         PetscCallA(PetscOptionsInsertString(PETSC_NULL_OPTIONS, "-dm_plex_gmsh_mark_vertices", ierr))
      end if
   end subroutine MEF90CtxCreate

#undef __FUNCT__
#define __FUNCT__ "MEF90CtxDestroy"
!!! author: Blaise Bourdin (2012-2014, bourdin@lsu.edu)
!!! author: Blaise Bourdin (2022, bourdin@mcmaster.ca)
!!!
!!!  MEF90CtxDestroy:
!!!

   subroutine MEF90CtxDestroy(MEF90Ctx, ierr)
      type(MEF90Ctx_Type), intent(INOUT)               :: MEF90Ctx
      PetscErrorCode, intent(OUT)                      :: ierr

      type(MEF90CtxGlobalOptions_Type)                 :: options

      PetscCall(MEF90CtxGlobalOptionsSetFromOptions(MEF90Ctx%comm, trim(MEF90Ctx%prefix), options, ierr))
      if (options%dryrun) then
         PetscCall(PetscOptionsLeft(PETSC_NULL_OPTIONS, ierr))
      end if
      if (.not. PetscObjectIsNull(MEF90Ctx%resultViewer)) then
         PetscCall(PetscViewerDestroy(MEF90Ctx%resultViewer, ierr))
      end if
   end subroutine MEF90CtxDestroy

#undef __FUNCT__
#define __FUNCT__ "MEF90CtxGlobalOptionsSetFromOptions"
!!! author: Blaise Bourdin (2012-2014, bourdin@lsu.edu)
!!! author: Blaise Bourdin (2026, bourdin@mcmaster.ca)
!!!
!!!  MEF90CtxGlobalOptionsSetFromOptions: reads the MEF90-wide options.
!!!      options is intent(inout): whatever it holds on entry is used as the default, so that a
!!!      caller can override the defaults carried by MEF90CtxGlobalOptions_Type before calling.
!!!

   subroutine MEF90CtxGlobalOptionsSetFromOptions(comm, prefix, options, ierr)
      MPIU_Comm, intent(IN)                           :: comm
      character(len=*), intent(IN)                    :: prefix
      type(MEF90CtxGlobalOptions_Type), intent(INOUT) :: options
      PetscErrorCode, intent(INOUT)                   :: ierr

      PetscCall(PetscOptionsBegin(comm, prefix, "Options for MEF90Ctx_Type", "mef90", ierr))
         PetscCall(PetscOptionsInt('-verbose', 'Verbosity: level', 'mef90', options%verbose, options%verbose, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsBool('-dryrun', 'Dry run in order to validate the options file. Use in combination with -h to print help or -verbose 1 to check input deck', 'mef90', options%dryrun, options%dryrun, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsEnum('-time_interpolation', 'Time: interpolation type', 'mef90', MEF90TimeInterpolationList, options%timeInterpolation, options%timeInterpolation, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-time_min', 'Time: min', 'mef90', options%timeMin, options%timeMin, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-time_max', 'Time: max', 'mef90', options%timeMax, options%timeMax, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsInt('-time_numstep', 'Time: number of time steps', 'mef90', options%timeNumStep, options%timeNumStep, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsInt('-time_skip', 'Time: number of time steps', 'mef90', options%timeSkip, options%timeSkip, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsInt('-time_numCycle', 'Time: number of cycles', 'mef90', options%timeNumCycle, options%timeNumCycle, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsEnum('-element_family', 'Element family (possibly overridden in application contexts)', 'mef90', MEF90ElementFamilyList, options%elementFamily, options%elementFamily, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsInt('-element_order', 'Element order (possibly overridden in application contexts)', 'mef90', options%elementOrder, options%elementOrder, PETSC_NULL_BOOL, ierr))
      PetscCall(PetscOptionsEnd(ierr))
   end subroutine MEF90CtxGlobalOptionsSetFromOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90CtxSetFromOptions"
!!! author: Blaise Bourdin (2012-2014, bourdin@lsu.edu)
!!! author: Blaise Bourdin (2026, bourdin@mcmaster.ca)
!!!
!!!  MEF90CtxSetFromOptions: initializes a MEF90Ctx_Type from options
!!!

   subroutine MEF90CtxSetFromOptions(self, ierr)
      class(MEF90Ctx_Type), intent(INOUT)             :: self
      PetscErrorCode, intent(INOUT)                   :: ierr

      type(MEF90CtxGlobalOptions_Type)                :: options

      PetscCall(MEF90CtxGlobalOptionsSetFromOptions(self%comm, trim(self%prefix), options, ierr))

      if (options%verbose > 0) then
         call self%view(PETSC_VIEWER_STDOUT_WORLD, ierr)
      end if
   end subroutine MEF90CtxSetFromOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90CtxView"
!!! author: Blaise Bourdin (2026, bourdin@mcmaster.ca)
!!!
!!!  MEF90CtxView: the default viewer for a MEF90Ctx_Type
!!!

   subroutine MEF90CtxView(self, viewer, ierr)
      class(MEF90Ctx_Type), intent(IN)                :: self
      type(tPetscViewer), intent(IN)                  :: viewer
      PetscErrorCode, intent(INOUT)                   :: ierr

      character(len=MEF90MXSTRLEN, kind=c_char)       :: IOBuffer
      character(len=MEF90MXSTRLEN, kind=c_char)       :: viewerType
      type(MEF90CtxGlobalOptions_Type)                :: options

      PetscCall(PetscViewerGetType(viewer, viewerType, ierr))
      if (viewerType /= 'ascii') return

      PetscCall(MEF90CtxGlobalOptionsSetFromOptions(self%comm, trim(self%prefix), options, ierr))
      PetscCall(PetscViewerASCIIPrintf(viewer, "MEF90 Global Context: \n", ierr))
      write (IOBuffer, "('  geometry file:       ',(A),'\n')") trim(self%geometryFile)
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('  result file:         ',(A),'\n')") trim(self%resultFile)
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('  log file:            ',(A),'\n')") trim(MEF90FilePrefix(self%resultFile))//'.log'
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('  verbose / dryrun:    ',I4,' / ',L1,'\n')") options%verbose, options%dryrun
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('  time interpolation:  ',(A),'\n')") trim(MEF90TimeInterpolationList(options%timeInterpolation + 1))
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('  time min / max:      ',ES12.5,' / ',ES12.5,'\n')") options%timeMin, options%timeMax
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('  time numstep / skip / numCycle: ',3(I6,' '),'\n')") options%timeNumStep, options%timeSkip, options%timeNumCycle
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('  element family / order: ',(A),' / ',I4,'\n')") trim(MEF90ElementFamilyList(options%elementFamily + 1)), options%elementOrder
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      PetscCall(PetscViewerASCIIPrintf(viewer, "\n", ierr))
   end subroutine MEF90CtxView

#undef __FUNCT__
#define __FUNCT__ "MEF90CtxGetTime"
!!! author: Blaise Bourdin (2012-2014, bourdin@lsu.edu)
!!!
!!!  MEF90CtxGetTime:
!!!

   subroutine MEF90CtxGetTime(MEF90Ctx, t, ierr)
      !!! The number of time steps this returns is size(t). For the exo interpolation it comes from
      !!! the file rather than from -time_numstep, so callers must use size(t) and not re-read
      !!! -time_numstep to size their per-step arrays or to detect the last step.
      type(MEF90Ctx_Type), intent(INOUT)              :: MEF90Ctx
      PetscReal, dimension(:), pointer                :: t
      PetscErrorCode, intent(OUT)                     :: ierr

      PetscReal                                       :: dt
      real                                            :: dummyR
      character(len=1)                                :: dummyS
      PetscErrorCode                                  :: exoErr
      integer                                         :: exoUnit
      type(MEF90CtxGlobalOptions_Type)                :: GlobalOptions
      PetscInt                                        :: i, j, CycleLength
      character(len=MEF90MXSTRLEN)                    :: IOBuffer

      i = 0 ! silence gfortran silly warning
      PetscCall(MEF90CtxGlobalOptionsSetFromOptions(MEF90Ctx%comm, trim(MEF90Ctx%prefix), GlobalOptions, ierr))
      select case (GlobalOptions%timeInterpolation)
      case (MEF90TimeInterpolation_linear)
         allocate (t(GlobalOptions%timeNumStep))
         dt = 0.0_kr
         if (GlobalOptions%timeNumStep > 1) then
            dt = (GlobalOptions%timeMax - GlobalOptions%timeMin) / real(GlobalOptions%timeNumStep - 1)
         end if
         t = [(GlobalOptions%timeMin + i * dt, i = 0, GlobalOptions%timeNumStep - 1)]
         t(GlobalOptions%timeNumStep) = GlobalOptions%timeMax

      case (MEF90TimeInterpolation_Vcycle)
         cycleLength = GlobalOptions%timeNumStep / GlobalOptions%timeNumCycle
         allocate (t(GlobalOptions%timeNumStep))
         dt = 0.0_kr
         if (GlobalOptions%timeNumStep > 1) then
            dt = (GlobalOptions%timeMax - GlobalOptions%timeMin) * GlobalOptions%timeNumCycle / real((GlobalOptions%timeNumStep - 1))
         end if
         do i = 1, GlobalOptions%timeNumCycle
            do j = 1, cycleLength
               t((i - 1) * cycleLength + j) = min(GlobalOptions%timeMin + 2.0_kr * real(j - 1) * dt, 2.0_kr * GlobalOptions%timeMax - GlobalOptions%timeMin - 2.0_kr * real(j - 1) * dt)
            end do
         end do
         t(cycleLength * GlobalOptions%timeNumCycle + 1:GlobalOptions%timeNumStep) = t(cycleLength * GlobalOptions%timeNumCycle)

      case (MEF90TimeInterpolation_quadratic)
         !! Natural time scale for the heat equation
         allocate (t(GlobalOptions%timeNumStep))
         dt = 0.0_kr
         if (GlobalOptions%timeNumStep > 1) then
            dt = (sqrt(GlobalOptions%timeMax) - sqrt(GlobalOptions%timeMin)) / real(GlobalOptions%timeNumStep - 1.0_kr)
         end if
         t = [((sqrt(GlobalOptions%timeMin) + real(i) * dt)**2, i=0, GlobalOptions%timeNumStep - 1)]
         t(GlobalOptions%timeNumStep) = GlobalOptions%timeMax

      case (MEF90TimeInterpolation_exo)
         PetscCall(PetscViewerExodusIIGetId(MEF90Ctx%resultViewer, exoUnit, ierr))
         if (.not. PetscObjectIsNull(MEF90Ctx%resultViewer)) then
            call EXINQ(exoUnit, EXTIMS, GlobalOptions%timeNumStep, dummyR, dummyS, exoErr)
            allocate (t(GlobalOptions%timeNumStep))
            call EXGATM(exoUnit, t, exoErr)
         else
            write (IOBuffer, "(A,'EXO input file must be open prior to calling MEF90Ctx_GetTime\n')") __FUNCT__
            SETERRQ(MEF90Ctx%Comm, PETSC_ERR_FILE_OPEN, IOBuffer)
         end if
      case Default
         write (IOBuffer, "(A,'Unimplemented time interpolation: ',I0,'\n')") __FUNCT__, GlobalOptions%timeInterpolation
         SETERRQ(MEF90Ctx%Comm, PETSC_ERR_ARG_OUTOFRANGE, IOBuffer)
      end select
      if ((GlobalOptions%verbose > 0) .and. (MEF90Ctx%rank == 0)) then
         PetscCall(PetscPrintf(MEF90Ctx%Comm, "Time values array:\n", ierr))
         PetscCall(PetscRealView(GlobalOptions%timeNumStep, t, PETSC_VIEWER_STDOUT_SELF, ierr))
         PetscCall(PetscPrintf(MEF90Ctx%Comm, "===\n", ierr))
      end if
   end subroutine MEF90CtxGetTime
end module m_MEF90_Ctx
