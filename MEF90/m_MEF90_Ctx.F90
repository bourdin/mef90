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

!!!
!!!
!!!  MEF90Ctx_Type: The class holding the MEF90-wide state and options
!!!
!!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!      2026      Blaise Bourdin bourdin@mcmaster.ca
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
      type(MEF90CtxGlobalOptions_Type)                :: globalOptions
      type(tPetscViewer)                              :: resultViewer
   contains
      procedure, pass(self) :: setFromOptions => MEF90CtxSetFromOptions
      procedure, pass(self) :: view => MEF90CtxView
   end type MEF90Ctx_Type

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90CtxCreate"
!!!
!!!
!!!  MEF90CtxCreate:
!!!
!!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022      Blaise Bourdin bourdin@mcmaster.ca
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
         !!! Old style calling sequence: geometryFile is <prefix>.gen, resultFile is <prefix>_out.gen
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

      !!! Not sure if this should be there, but PETSc's gmsh reader defaults to ignoring vertex sets, which we defintely don't want...
      if (MEF90FileExtension(MEF90Ctx%geometryfile) == 'msh') then
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

      if (MEF90Ctx%globalOptions%dryrun) then
         PetscCall(PetscOptionsLeft(PETSC_NULL_OPTIONS, ierr))
      end if
      if (.not. PetscObjectIsNull(MEF90Ctx%resultViewer)) then
         PetscCall(PetscViewerDestroy(MEF90Ctx%resultViewer, ierr))
      end if
   end subroutine MEF90CtxDestroy

#undef __FUNCT__
#define __FUNCT__ "MEF90CtxSetFromOptions"
!!!
!!!
!!!  MEF90CtxSetFromOptions: initializes a MEF90Ctx_Type from options
!!!
!!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!      2026      Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90CtxSetFromOptions(self, ierr)
      class(MEF90Ctx_Type), intent(INOUT)             :: self
      PetscErrorCode, intent(INOUT)                   :: ierr

      PetscCall(PetscOptionsBegin(self%comm, trim(self%prefix), "Options for MEF90Ctx_Type", "mef90", ierr))
         PetscCall(PetscOptionsInt('-verbose', 'Verbosity: level', 'mef90', self%globalOptions%verbose, self%globalOptions%verbose, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsBool('-dryrun', 'Dry run in order to validate the options file. Use in combination with -h to print help or -verbose 1 to check input deck', 'mef90', self%globalOptions%dryrun, self%globalOptions%dryrun, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsEnum('-time_interpolation', 'Time: interpolation type', 'mef90', MEF90TimeInterpolationList, self%globalOptions%timeInterpolation, self%globalOptions%timeInterpolation, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-time_min', 'Time: min', 'mef90', self%globalOptions%timeMin, self%globalOptions%timeMin, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-time_max', 'Time: max', 'mef90', self%globalOptions%timeMax, self%globalOptions%timeMax, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsInt('-time_numstep', 'Time: number of time steps', 'mef90', self%globalOptions%timeNumStep, self%globalOptions%timeNumStep, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsInt('-time_skip', 'Time: number of time steps', 'mef90', self%globalOptions%timeSkip, self%globalOptions%timeSkip, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsInt('-time_numCycle', 'Time: number of cycles', 'mef90', self%globalOptions%timeNumCycle, self%globalOptions%timeNumCycle, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsEnum('-element_family', 'Element family (possibly overridden in application contexts)', 'mef90', MEF90ElementFamilyList, self%globalOptions%elementFamily, self%globalOptions%elementFamily, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsInt('-element_order', 'Element order (possibly overridden in application contexts)', 'mef90', self%globalOptions%elementOrder, self%globalOptions%elementOrder, PETSC_NULL_BOOL, ierr))
      PetscCall(PetscOptionsEnd(ierr))

      if (self%globalOptions%verbose > 0) then
         call self%view(PETSC_VIEWER_STDOUT_WORLD, ierr)
      end if
   end subroutine MEF90CtxSetFromOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90CtxView"
!!!
!!!
!!!  MEF90CtxView: the default viewer for a MEF90Ctx_Type
!!!
!!!  (c) 2026      Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90CtxView(self, viewer, ierr)
      class(MEF90Ctx_Type), intent(IN)                :: self
      type(tPetscViewer), intent(IN)                  :: viewer
      PetscErrorCode, intent(INOUT)                   :: ierr

      character(len=MEF90MXSTRLEN, kind=c_char)       :: IOBuffer
      character(len=MEF90MXSTRLEN, kind=c_char)       :: viewerType

      PetscCall(PetscViewerGetType(viewer, viewerType, ierr))
      if (viewerType /= 'ascii') return

      PetscCall(PetscViewerASCIIPrintf(viewer, "MEF90 Global Context: \n", ierr))
      write (IOBuffer, "('  geometry file:       ',(A),'\n')") trim(self%geometryFile)
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('  result file:         ',(A),'\n')") trim(self%resultFile)
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('  log file:            ',(A),'\n')") trim(MEF90FilePrefix(self%resultFile))//'.log'
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('  verbose / dryrun:    ',I4,' / ',L1,'\n')") self%globalOptions%verbose, self%globalOptions%dryrun
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('  time interpolation:  ',(A),'\n')") trim(MEF90TimeInterpolationList(self%globalOptions%timeInterpolation + 1))
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('  time min / max:      ',ES12.5,' / ',ES12.5,'\n')") self%globalOptions%timeMin, self%globalOptions%timeMax
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('  time numstep / skip / numCycle: ',3(I6,' '),'\n')") self%globalOptions%timeNumStep, self%globalOptions%timeSkip, self%globalOptions%timeNumCycle
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('  element family / order: ',(A),' / ',I4,'\n')") trim(MEF90ElementFamilyList(self%globalOptions%elementFamily + 1)), self%globalOptions%elementOrder
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      PetscCall(PetscViewerASCIIPrintf(viewer, "\n", ierr))
   end subroutine MEF90CtxView

#undef __FUNCT__
#define __FUNCT__ "MEF90CtxGetTime"
!!!
!!!
!!!  MEF90CtxGetTime:
!!!
!!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!

   subroutine MEF90CtxGetTime(MEF90Ctx, t, ierr)
      !!! MEF90Ctx has the target attribute since GlobalOptions below aliases its options, which this
      !!! routine updates when reading the time values from an exodus file
      type(MEF90Ctx_Type), target, intent(INOUT)      :: MEF90Ctx
      PetscReal, dimension(:), pointer                :: t
      PetscErrorCode, intent(OUT)                     :: ierr

      PetscReal                                       :: dt
      real                                            :: dummyR
      character(len=1)                                :: dummyS
      PetscErrorCode                                  :: exoErr
      integer                                         :: exoUnit
      type(MEF90CtxGlobalOptions_Type), pointer       :: GlobalOptions
      PetscInt                                        :: i, j, CycleLength
      character(len=MEF90MXSTRLEN)                    :: IOBuffer

      i = 0 ! silence gfortran silly warning
      GlobalOptions => MEF90Ctx%globalOptions
      select case (GlobalOptions%timeInterpolation)
      case (MEF90TimeInterpolation_linear)
         allocate (t(GlobalOptions%timeNumStep))
         dt = 0.0_kr
         if (GlobalOptions%timeNumStep > 1) then
            dt = (GlobalOptions%timeMax - GlobalOptions%timeMin) / real(GlobalOptions%timeNumStep - 1.0_kr)
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
         !!! Natural time scale for the heat equation
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
