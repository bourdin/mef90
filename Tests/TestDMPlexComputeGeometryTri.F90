program TestDMPlexComputeGeometry
#include <petsc/finclude/petsc.h>
   use m_MEF90
   use petsc
   implicit none(type, external)

   PetscErrorCode                      :: ierr
   type(MEF90Ctx_Type), target          :: MEF90Ctx
   type(MEF90CtxGlobalOptions_Type)    :: MEF90GlobalOptions_default
   type(MEF90CtxGlobalOptions_Type)    :: MEF90GlobalOptions
   type(tDM), target                    :: dm, dmDist
   PetscBool                           :: interpolate = PETSC_TRUE
   character(len=MEF90MXSTRLEN)       :: IOBuffer
   PetscInt                            :: set, cell
   type(tIS)                           :: CSIS, FSIS, CellIS
   PetscInt, dimension(:), pointer       :: setID, cellID
   PetscReal, dimension(:), pointer      :: v0
   PetscReal, dimension(:), pointer      :: B, Binv
   PetscReal                           :: detB
   PetscInt                            :: dim
   PetscMPIInt                         :: rank, numProc

   PetscReal                           :: vol
   PetscReal, dimension(:), pointer      :: centroid, normal

   MEF90GlobalOptions_default%verbose = 1
   MEF90GlobalOptions_default%dryrun = PETSC_FALSE
   MEF90GlobalOptions_default%timeInterpolation = MEF90TimeInterpolation_linear
   MEF90GlobalOptions_default%timeMin = 0.0_kr
   MEF90GlobalOptions_default%timeMax = 1.0_kr
   MEF90GlobalOptions_default%timeNumStep = 11
   MEF90GlobalOptions_default%timeSkip = 0
   MEF90GlobalOptions_default%timeNumCycle = 1
   MEF90GlobalOptions_default%fileFormat = MEF90FileFormat_EXOSingle

   PetscCall(PetscInitialize(PETSC_NULL_CHARACTER, ierr))
   PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr))
   PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, numProc, ierr))

   call MEF90Initialize(PETSC_COMM_WORLD, ierr)
   call MEF90CtxCreate(PETSC_COMM_WORLD, MEF90Ctx, "", ierr)
   MEF90GlobalOptions = MEF90GlobalOptions_default
   PetscCallA(MEF90CtxGlobalOptionsSetFromOptions(MEF90Ctx%comm, trim(MEF90Ctx%prefix), MEF90GlobalOptions, ierr))
   call MEF90Ctx%setFromOptions(ierr); CHKERRQ(ierr)

   PetscCall(PetscPrintf(PETSC_COMM_WORLD, MEF90Ctx%geometryfile, ierr))
   PetscCall(DMPlexCreateFromFile(PETSC_COMM_WORLD, MEF90Ctx%geometryfile, PETSC_NULL_CHARACTER, interpolate, dm, ierr))
   PetscCall(DMPlexDistributeSetDefault(dm, PETSC_FALSE, ierr))
   PetscCall(DMSetFromOptions(dm, ierr))
   PetscCall(DMGetDimension(dm, dim, ierr))
   PetscCall(DMViewFromOptions(dm, PETSC_NULL_OBJECT, "-dm_view", ierr))

   PetscCall(DMPlexDistribute(dm, 0, PETSC_NULL_SF, dmDist, ierr))
   if (numProc == 1) then
      dmDist = dm
   end if
   PetscCall(DMViewFromOptions(dmDist, PETSC_NULL_OBJECT, "-dm_view", ierr))
   PetscCall(DMGetDimension(dm, dim, ierr))

   allocate (v0(dim))
   allocate (B(dim**2))
   allocate (Binv(dim**2))
   PetscCall(DMGetLabelIdIS(dm, "Cell Sets", CSIS, ierr))
   PetscCall(ISGetIndices(CSIS, setID, ierr))
   do set = 1, size(setID)
      write (IOBuffer, *) 'Cell Set', setID(set), '\n'
      PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
      PetscCall(DMGetStratumIS(dm, "Cell Sets", setID(set), cellIS, ierr))
      PetscCall(ISGetIndices(cellIS, cellID, ierr))
      do cell = 1, size(cellID)
         PetscCall(DMPlexComputeCellGeometryAffineFEM(dm, cellID(cell), v0, B, Binv, detB, ierr))
         write (*, *) 'cell ', cellID(cell)
         write (*, *) '   v0      ', v0
         write (*, *) '   B       ', B
         write (*, *) '   Binv    ', Binv
         write (*, *) '   detB    ', detB
         !PetscCall(DMPlexComputeCellGeometryFVM(dm,cellID(cell),vol,centroid,normal,ierr))
         !write(*,*) '   vol      ', vol
         !write(*,*) '   centroid ', centroid
         !write(*,*) '   normal   ', normal
      end do
      PetscCall(PetscPrintf(PETSC_COMM_SELF, "\n", ierr))
      PetscCall(ISRestoreIndices(cellIS, cellID, ierr))
      PetscCall(ISDestroy(cellIS, ierr))
   end do
   PetscCall(ISRestoreIndices(CSIS, setID, ierr))
   PetscCall(ISDestroy(CSIS, ierr))
   deallocate (Binv)
   deallocate (B)
   deallocate (v0)

   allocate (centroid(dim))
   allocate (normal(dim))
   PetscCall(DMGetLabelIdIS(dm, "Face Sets", fsIS, ierr))

   PetscCall(ISGetIndices(FSIS, setID, ierr))
   do set = 1, size(setID)
      write (IOBuffer, *) 'Face Set', setID(set), '\n'
      PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
      PetscCall(DMGetStratumIS(dm, "Face Sets", setID(set), cellIS, ierr))
      PetscCall(ISGetIndices(cellIS, cellID, ierr))
      do cell = 1, size(cellID)
         !PetscCall(DMPlexComputeCellGeometryAffineFEM(dm,cellID(cell),v0,B,Binv,detB,ierr))
         !write(*,*) 'cell ',cellID(cell)
         !write(*,*) '   v0   ',v0
         !write(*,*) '   B    ',B
         !write(*,*) '   Binv ',Binv
         !write(*,*) '   detB ',detB
         PetscCall(DMPlexComputeCellGeometryFVM(dm, cellID(cell), vol, centroid, normal, ierr))
         write (*, *) '   vol         ', vol
         write (*, *) '   centroid    ', centroid
         write (*, *) '   normal      ', normal
      end do
      PetscCall(ISRestoreIndices(cellIS, cellID, ierr))
      PetscCall(ISDestroy(cellIS, ierr))
   end do
   PetscCall(ISRestoreIndices(FSIS, setID, ierr))
   PetscCall(ISDestroy(FSIS, ierr))

   call MEF90CtxDestroy(MEF90Ctx, ierr)
   call MEF90Finalize(ierr)
   call PetscFinalize(ierr)
end program TestDMPlexComputeGeometry
