module m_MEF90_EXO
#include "petsc/finclude/petsc.h"
   use m_MEF90_Parameters
   use m_MEF90_Utils
   use m_MEF90_Elements
   use m_MEF90_Ctx
   use m_MEF90_DMPlex

   implicit none(type)
#include "../mef90version.h"

   private
   PetscInt, public                                 :: exo_ver

   public :: MEF90CtxOpenEXO
   public :: MEF90CtxCloseEXO
   public :: MEF90EXOFormat
   public :: MEF90EXODMView
   public :: MEF90EXOVecView
   public :: MEF90EXOVecLoad

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90CtxOpenEXO"
!!!
!!!
!!!  MEF90CtxOpenEXO:
!!!
!!!  (c) 2022 Alexis Marboeuf marboeua@mcmaster.ca
!!!      2022 Blaise Bourdin  bourdin@mcmaster.ca
!!!
   subroutine MEF90CtxOpenEXO(MEF90Ctx, Viewer, mode, ierr)
      type(MEF90Ctx_Type), intent(IN)                  :: MEF90Ctx
      type(tPetscViewer), intent(INOUT)               :: Viewer
      type(ePetscFileMode), intent(IN)                 :: mode
      PetscErrorCode, intent(INOUT)                    :: ierr

      integer                                         :: opts

#ifdef PETSC_USE_DEBUG
      opts = EXVRBS + EXDEBG
#else
      opts = 0
#endif
      call exopts(opts, ierr)
      PetscCall(PetscViewerExodusIIOpen(MEF90Ctx % Comm, MEF90Ctx % resultFile, mode, Viewer, ierr))
   end subroutine MEF90CtxOpenEXO

#undef __FUNCT__
#define __FUNCT__ "MEF90CtxCloseEXO"
!!!
!!!
!!!  MEF90CtxCloseEXO:
!!!
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!      2022      Blaise Bourdin  bourdin@mcmaster.ca
!!!
   subroutine MEF90CtxCloseEXO(Viewer, ierr)
      type(tPetscViewer), intent(INOUT)                :: Viewer
      PetscErrorCode, intent(INOUT)                    :: ierr

      PetscCall(PetscViewerDestroy(Viewer, ierr))
   end subroutine MEF90CtxCloseEXO

#undef __FUNCT__
#define __FUNCT__ "MEF90EXOFormat"
!!!
!!!
!!!  MEF90EXOFormat:
!!!
!!!  (c) 2012-2022 Blaise Bourdin bourdin@lsu.edu
!!!           2022 Alexis Marboeuf marboeua@mcmaster.ca
!!!           2024 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine MEF90EXOFormat(Viewer, nameG, nameC, nameV, time, ierr)
      type(tPetscViewer), intent(IN)                         :: Viewer
      character(len=*), dimension(:), intent(IN)              :: nameG, nameC, nameV
      PetscReal, dimension(:), pointer                        :: time
      PetscErrorCode, intent(INOUT)                          :: ierr

      PetscInt                                              :: numCS, numC
      integer                                               :: i, exoid
      PetscInt                                              :: step
      character(len=MXSTLN)                                 :: sJunk
      PetscReal                                             :: rJunk
      logical, dimension(:, :), pointer                        :: truthtable

      if (size(nameV) > 0) then
         PetscCall(PetscViewerExodusIISetNodalVariable(Viewer, size(nameV), ierr))
         do i = 1, size(nameV)
            PetscCall(PetscViewerExodusIISetNodalVariableName(Viewer, i - 1, nameV(i), ierr))
         end do
      end if
      if (size(nameC) > 0) then
         PetscCall(PetscViewerExodusIISetZonalVariable(Viewer, size(nameC), ierr))
         do i = 1, size(nameC)
            PetscCall(PetscViewerExodusIISetZonalVariableName(Viewer, i - 1, nameC(i), ierr))
         end do
      end if

      if (.not. associated(time)) then
         SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_POINTER, "Time value must be allocated prior to calling MEF90EXOFormat")
         stop
      end if

      PetscCall(PetscViewerExodusIIGetId(Viewer, exoid, ierr))
      call exinq(exoid, EX_INQ_ELEM_BLK, numCS, rJunk, sjunk, ierr)

   !!! Write truth tables
      numC = size(nameC)
      if (numC > 0) then
         allocate (truthtable(numCS, numC))
         truthtable = .true.
         call expvtt(exoid, numCS, numC, truthtable, ierr)
         deallocate (truthtable)
      end if

      do step = 1, size(time)
         call exptim(exoid, step, time(step), ierr)
      end do
      PetscCall(PetscViewerFlush(Viewer, ierr))
   end subroutine MEF90EXOFormat

#undef __FUNCT__
#define __FUNCT__ "MEF90EXODMView"
!!!
!!!
!!!  MEF90EXODMView:
!!!
!!!
!!!  (c) 2022 Alexis Marboeuf marboeua@mcmaster.ca
!!!
   subroutine MEF90EXODMView(dm, Viewer, order, ierr)
      type(tPetscViewer), intent(IN)                      :: Viewer
      type(tDM), intent(IN)                               :: dm
      PetscInt, intent(IN)                                :: order
      PetscErrorCode, intent(INOUT)                       :: ierr

      character(len=PETSC_MAX_PATH_LEN)                  :: IOBuffer

      if ((order > 2) .or. (order < 1)) then
         write (IOBuffer, '("Unsupported polynomial order ", I2, " not in [1,2]")') order
         SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, IOBuffer)
      end if
      PetscCall(PetscViewerExodusIISetOrder(Viewer, order, ierr))
      PetscCall(DMView(dm, Viewer, ierr))
   end subroutine MEF90EXODMView

#undef __FUNCT__
#define __FUNCT__ "MEF90EXOVecView"
!!!
!!!
!!!  MEF90EXOVecView:
!!!
!!!
!!!  (c) 2022 Alexis Marboeuf marboeua@mcmaster.ca
!!!      2023 Blaise Bourdin  bourdin@mcmaster.ca
!!!
   subroutine MEF90EXOVecView(v, sf, invSF, Viewer, step, bs, ierr)
      type(tVec), intent(IN)                              :: v
      type(tPetscSF), intent(IN)                          :: sf, invSF
      type(tPetscViewer), intent(IN)                      :: Viewer
      PetscInt, intent(IN)                                :: step
      PetscInt, intent(IN)                                :: bs
      PetscErrorCode, intent(INOUT)                       :: ierr

      integer                                            :: exoid
      integer                                            :: offsetN, offsetZ
      type(tVec)                                         :: iov
      character(len=PETSC_MAX_PATH_LEN)                  :: vecname, IOBuffer

      offsetN = -1
      offsetZ = -1
      PetscCall(PetscObjectGetName(v, vecname, ierr))

      PetscCall(MEF90VecCreateIO(iov, bs, sf, ierr))
      PetscCall(PetscObjectSetName(iov, vecname, ierr))
      PetscCall(MEF90VecCopySF(v, iov, sf, ierr))

      PetscCall(PetscViewerExodusIIGetNodalVariableIndex(Viewer, vecname, offsetN, ierr))
      PetscCall(PetscViewerExodusIIGetZonalVariableIndex(Viewer, vecname, offsetZ, ierr))
      PetscCall(PetscViewerExodusIIGetId(Viewer, exoid, ierr))
      if (offsetN >= 0) then
         PetscCall(MEF90EXOVecViewNodal_Private(iov, exoid, step, offsetN + 1, ierr))
      else if (offsetZ >= 0) then
         PetscCall(MEF90EXOVecViewZonal_Private(iov, exoid, step, offsetZ + 1, ierr))
      else
         write (IOBuffer, '("Could not find nodal or zonal variable ", A, " in exodus file. ")') trim(vecname)
         SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED, IOBuffer)
      end if
      PetscCall(VecDestroy(iov, ierr))
   end subroutine MEF90EXOVecView

#undef __FUNCT__
#define __FUNCT__ "MEF90EXOVecLoad"
!!!
!!!
!!!  MEF90EXOVecLoad:
!!!
!!!
!!!  (c) 2022 Alexis Marboeuf marboeua@mcmaster.ca
!!!
   subroutine MEF90EXOVecLoad(v, sf, invSF, Viewer, step, bs, ierr)
      type(tVec), intent(INOUT)                           :: v
      type(tPetscSF), intent(IN)                          :: sf, invSF
      type(tPetscViewer), intent(IN)                      :: Viewer
      PetscInt, intent(IN)                                :: step
      PetscInt, intent(IN)                                :: bs
      PetscErrorCode, intent(INOUT)                       :: ierr

      integer                                            :: exoid
      integer                                            :: offsetN, offsetZ
      type(tVec)                                         :: iov
      type(tDM)                                          :: locDM, oDM
      type(tVec)                                         :: vGlob
      character(len=PETSC_MAX_PATH_LEN)                  :: vecname, IOBuffer

      offsetN = -1
      offsetZ = -1
      PetscCall(PetscObjectGetName(v, vecname, ierr))

      PetscCall(MEF90VecCreateIO(iov, bs, sf, ierr))
      PetscCall(PetscObjectSetName(iov, vecname, ierr))

      PetscCall(PetscViewerExodusIIGetNodalVariableIndex(Viewer, vecname, offsetN, ierr))
      PetscCall(PetscViewerExodusIIGetZonalVariableIndex(Viewer, vecname, offsetZ, ierr))
      PetscCall(PetscViewerExodusIIGetId(Viewer, exoid, ierr))
      if (offsetN >= 0) then
         PetscCall(MEF90EXOVecLoadNodal_Private(iov, exoid, step, offsetN + 1, ierr))
      else if (offsetZ >= 0) then
         PetscCall(MEF90EXOVecLoadZonal_Private(iov, exoid, step, offsetZ + 1, ierr))
      else
         write (IOBuffer, '("Could not find nodal or zonal variable ", A, " in exodus file. ")') trim(vecname)
         SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED, IOBuffer)
      end if
      PetscCall(MEF90VecCopySF(iov, v, invSF, ierr))
      PetscCall(VecDestroy(iov, ierr))

      !!! Make sure that the halo values in the local vector are updated by doing a L2G followed by a G2L
      !!! This is probably only needed for nodal vectors, unless we have distributed the mesh with an overlap.
      PetscCall(VecGetDM(V, locDM, ierr))
      PetscCall(DMGetOutputDM(locDM, oDM, ierr))
      PetscCall(DMGetGlobalVector(oDM, vGlob, ierr))
      PetscCall(DMLocalToGlobal(oDM, v, INSERT_ALL_VALUES, vGlob, ierr))
      PetscCall(DMGlobalToLocal(oDM, vGlob, INSERT_ALL_VALUES, v, ierr))
      PetscCall(DMRestoreGlobalVector(oDM, vGlob, ierr))

   end subroutine MEF90EXOVecLoad

#undef __FUNCT__
#define __FUNCT__ "MEF90EXOVecViewNodal_Private"
   subroutine MEF90EXOVecViewNodal_Private(v, exoid, step, offset, ierr)
      integer, intent(IN)               :: exoid
      PetscInt, intent(IN)              :: step
      integer, intent(IN)               :: offset
      type(tVec), intent(IN)            :: v
      PetscErrorCode, intent(INOUT)     :: ierr

      PetscInt                         :: xs, xe, bs
      integer                          :: c
      PetscScalar, dimension(:), pointer :: varray
      type(tVec)                       :: vComp
      type(tIS)                        :: compIS

      PetscCall(VecGetOwnershipRange(v, xs, xe, ierr))
      PetscCall(VecGetBlockSize(v, bs, ierr))

      if (bs == 1) then
         PetscCall(VecGetArrayRead(v, varray, ierr))
         call expnvs(exoid, step, offset, xs + 1, xe - xs, varray, ierr)
         PetscCall(VecRestoreArrayRead(v, varray, ierr))
      else
         PetscCall(ISCreateStride(PETSC_COMM_WORLD, (xe - xs) / bs, xs, bs, compIS, ierr))
         do c = 0, bs - 1
            PetscCall(ISStrideSetStride(compIS, (xe - xs) / bs, xs + c, bs, ierr))
            PetscCall(VecGetSubVector(v, compIS, vComp, ierr))
            PetscCall(VecGetArrayRead(vComp, varray, ierr))
            call expnvs(exoid, step, offset + c, xs / bs + 1, (xe - xs) / bs, varray, ierr)
            PetscCall(VecRestoreArrayRead(vComp, varray, ierr))
            PetscCall(VecRestoreSubVector(v, compIS, vComp, ierr))
         end do ! c
         PetscCall(ISDestroy(compIS, ierr))
      end if ! bs
   end subroutine MEF90EXOVecViewNodal_Private

#undef __FUNCT__
#define __FUNCT__ "MEF90EXOVecLoadNodal_Private"
   subroutine MEF90EXOVecLoadNodal_Private(v, exoid, step, offset, ierr)
      integer, intent(IN)               :: exoid
      PetscInt, intent(IN)              :: step
      integer, intent(IN)               :: offset
      type(tVec), intent(INOUT)         :: v
      PetscErrorCode, intent(INOUT)     :: ierr

      PetscInt                         :: xs, xe, bs
      integer                          :: c
      PetscScalar, dimension(:), pointer :: varray
      type(tVec)                       :: vComp
      type(tIS)                        :: compIS

      PetscCall(VecGetOwnershipRange(v, xs, xe, ierr))
      PetscCall(VecGetBlockSize(v, bs, ierr))
      if (bs == 1) then
         PetscCall(VecGetArrayRead(v, varray, ierr))
         call exgnnv(exoid, step, offset, xs + 1, xe - xs, varray, ierr)
         PetscCall(VecRestoreArrayRead(v, varray, ierr))
      else
         PetscCall(ISCreateStride(PETSC_COMM_WORLD, (xe - xs) / bs, xs, bs, compIS, ierr))
         do c = 0, bs - 1
            PetscCall(ISStrideSetStride(compIS, (xe - xs) / bs, xs + c, bs, ierr))
            PetscCall(VecGetSubVector(v, compIS, vComp, ierr))
            PetscCall(VecGetArrayRead(vComp, varray, ierr))
            call exgnnv(exoid, step, offset + c, xs / bs + 1, (xe - xs) / bs, varray, ierr)
            PetscCall(VecRestoreArrayRead(vComp, varray, ierr))
            PetscCall(VecISCopy(v, compIS, SCATTER_FORWARD, vComp, ierr))
            PetscCall(VecRestoreSubVector(v, compIS, vComp, ierr))
         end do
         PetscCall(ISDestroy(compIS, ierr))
      end if
   end subroutine MEF90EXOVecLoadNodal_Private

#undef __FUNCT__
#define __FUNCT__ "MEF90EXOVecViewZonal_Private"
   subroutine MEF90EXOVecViewZonal_Private(v, exoid, step, offset, ierr)
      integer, intent(IN)               :: exoid
      PetscInt, intent(IN)              :: step
      integer, intent(IN)               :: offset
      type(tVec), intent(IN)            :: v
      PetscErrorCode, intent(INOUT)     :: ierr

      PetscInt                         :: xs, xe, bs, numCS, set, csLocalSize, csxs
      integer                          :: c
      PetscScalar, dimension(:), pointer :: varray
      PetscInt, dimension(:), pointer    :: csID, csSize
      type(tVec)                       :: vComp
      type(tIS)                        :: compIS
      character(len=MXSTLN)            :: elemType
      PetscMPIInt                      :: rank

      csxs = 0_ki
      PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr))
      numCS = exinqi(exoid, EX_INQ_ELEM_BLK)
      allocate (csID(numCS))
      allocate (csSize(numCS))
      call exgebi(exoid, csID, ierr)
      do set = 1, numCS
         call exgelb(exoid, csID(set), elemType, csSize(set), PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)
      end do
      PetscCall(VecGetOwnershipRange(v, xs, xe, ierr))
      PetscCall(VecGetBlockSize(v, bs, ierr))
      if (bs > 1) then
         PetscCall(ISCreateStride(PETSC_COMM_WORLD, (xe - xs) / bs, xs, bs, compIS, ierr))
      end if
      do set = 1, numCS
         !  range of indices for set setID[set]: csxs:csxs + csSize[set]-1
         !  local slice of zonal values:         xs/bs,xm/bs-1
         !  intersection:                        max(xs/bs,csxs),min(xm/bs-1,csxs + csSize[set]-1)
         csLocalSize = max(0, min(xe / bs, csxs + csSize(set)) - max(xs / bs, csxs))
         if (bs == 1) then
            PetscCall(VecGetArrayRead(v, varray, ierr))
            call expevs(exoid, step, offset, csID(set), max(xs - csxs, 0) + 1, csLocalSize, varray, ierr)
            PetscCall(VecRestoreArrayRead(v, varray, ierr))
         else
            do c = 0, bs - 1
               PetscCall(ISStrideSetStride(compIS, (xe - xs) / bs, xs + c, bs, ierr))
               PetscCall(VecGetSubVector(v, compIS, vComp, ierr))
               PetscCall(VecGetArrayRead(vComp, varray, ierr))
               call expevs(exoid, step, offset + c, csID(set), max(xs / bs - csxs, 0) + 1, csLocalSize, varray(max(0, csxs - xs / bs) + 1:max(0, csxs - xs / bs) + csLocalSize), ierr)
               PetscCall(VecRestoreArrayRead(vComp, varray, ierr))
               PetscCall(VecRestoreSubVector(v, compIS, vComp, ierr))
            end do
         end if
         csxs = csxs + csSize(set)
      end do
      if (bs > 1) then
         PetscCall(ISDestroy(compIS, ierr))
      end if
      deallocate (csID)
      deallocate (csSize)
   end subroutine MEF90EXOVecViewZonal_Private

#undef __FUNCT__
#define __FUNCT__ "MEF90EXOVecLoadZonal_Private"
   subroutine MEF90EXOVecLoadZonal_Private(v, exoid, step, offset, ierr)
      integer, intent(IN)               :: exoid
      PetscInt, intent(IN)              :: step
      integer, intent(IN)               :: offset
      type(tVec), intent(INOUT)         :: v
      PetscErrorCode, intent(INOUT)     :: ierr

      PetscInt                         :: xs, xe, bs, numCS, set, csLocalSize, csxs
      integer                          :: c
      PetscScalar, dimension(:), pointer :: varray
      PetscInt, dimension(:), pointer    :: csID, csSize
      type(tVec)                       :: vComp
      type(tIS)                        :: compIS
      character(len=MXSTLN)            :: elemType
      PetscMPIInt                      :: rank

      csxs = 0_ki
      PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr))
      numCS = exinqi(exoid, EX_INQ_ELEM_BLK)
      allocate (csID(numCS))
      allocate (csSize(numCS))
      call exgebi(exoid, csID, ierr)
      do set = 1, numCS
         call exgelb(exoid, csID(set), elemType, csSize(set), PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)
      end do
      PetscCall(VecGetOwnershipRange(v, xs, xe, ierr))
      PetscCall(VecGetBlockSize(v, bs, ierr))
      if (bs > 1) then
         PetscCall(ISCreateStride(PETSC_COMM_WORLD, (xe - xs) / bs, xs, bs, compIS, ierr))
      end if
      do set = 1, numCS
         !  range of indices for set setID[set]: csxs:csxs + csSize[set]-1
         !  local slice of zonal values:         xs/bs,xm/bs-1
         !  intersection:                        max(xs/bs,csxs),min(xm/bs-1,csxs + csSize[set]-1)
         csLocalSize = max(0, min(xe / bs, csxs + csSize(set)) - max(xs / bs, csxs))
         if (bs == 1) then
            PetscCall(VecGetArray(v, varray, ierr))
            call exgnev(exoid, step, offset, csID(set), csSize(set), max(xs - csxs, 0) + 1, csLocalSize, varray, ierr)
            PetscCall(VecRestoreArray(v, varray, ierr))
         else
            do c = 0, bs - 1
               PetscCall(ISStrideSetStride(compIS, (xe - xs) / bs, xs + c, bs, ierr))
               PetscCall(VecGetSubVector(v, compIS, vComp, ierr))
               PetscCall(VecGetArray(vComp, varray, ierr))
               call exgnev(exoid, step, offset + c, csID(set), 0_ki, max(xs / bs - csxs, 0) + 1, csLocalSize, varray(max(0, csxs - xs / bs) + 1:max(0, csxs - xs / bs) + csLocalSize), ierr)
               ! the 5th argument of exgnev is unused
               PetscCall(VecRestoreArray(vComp, varray, ierr))
               PetscCall(VecISCopy(v, compIS, SCATTER_FORWARD, vComp, ierr))
               PetscCall(VecRestoreSubVector(v, compIS, vComp, ierr))
            end do
         end if
         csxs = csxs + csSize(set)
      end do
      if (bs > 1) then
         PetscCall(ISDestroy(compIS, ierr))
      end if
      deallocate (csID)
      deallocate (csSize)
   end subroutine MEF90EXOVecLoadZonal_Private

#undef __FUNCT__
#define __FUNCT__ "MEF90EXOVecViewSide_Private"
   subroutine MEF90EXOVecViewSide_Private(v, exoid, step, offset, ierr)
      integer, intent(IN)               :: exoid
      PetscInt, intent(IN)              :: step, offset
      type(tVec), intent(IN)            :: v
      PetscErrorCode, intent(INOUT)     :: ierr

      PetscInt                         :: xs, xe, bs, c, numSS, set, ssLocalSize, ssxs, sscs
      PetscScalar, dimension(:), pointer :: varray
      PetscInt, dimension(:), pointer    :: ssID, ssSize
      type(tVec)                       :: vComp
      type(tIS)                        :: compIS
      PetscMPIInt                      :: rank

      ssxs = 0_ki
      sscs = 0_ki
      PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr))
      numSS = exinqi(exoid, EX_INQ_SIDE_SETS)
      allocate (ssID(numSS))
      allocate (ssSize(numSS))
      call exgssi(exoid, ssID, ierr)
      do set = 1, numSS
         call exgsp(exoid, ssID(set), ssSize(set), PETSC_NULL_INTEGER, ierr)
      end do
      PetscCall(VecGetOwnershipRange(v, xs, xe, ierr))
      PetscCall(VecGetBlockSize(v, bs, ierr))
      do set = 1, numSS
         !  range of indices for set setID[set]: csxs:csxs + csSize[set]-1
         !  local slice of zonal values:         xs/bs,xm/bs-1
         !  intersection:                        max(xs/bs,csxs),min(xm/bs-1,csxs + csSize[set]-1)
         ssLocalSize = max(0, min(xe / bs, ssxs + ssSize(set)) - max(xs / bs, ssxs))
         if (bs == 1) then
            PetscCall(VecGetArrayRead(v, varray, ierr))
            call expssv(exoid, step, offset, ssID(set), ssLocalSize, varray, ierr)
            PetscCall(VecRestoreArrayRead(v, varray, ierr))
         else
            PetscCall(ISCreateStride(PETSC_COMM_WORLD, ssLocalSize, xs + sscs, bs, compIS, ierr))
            do c = 0, bs - 1
               PetscCall(ISStrideSetStride(compIS, ssLocalSize, xs + sscs + c, bs, ierr))
               PetscCall(VecGetSubVector(v, compIS, vComp, ierr))
               PetscCall(VecGetArrayRead(vComp, varray, ierr))
               call exppv(exoid, step, EX_SIDE_SET, offset + c, ssID(set), max(xs / bs - ssxs, 0) + 1, ssLocalSize, varray, ierr)
               PetscCall(VecRestoreArrayRead(vComp, varray, ierr))
               PetscCall(VecRestoreSubVector(v, compIS, vComp, ierr))
            end do
            PetscCall(ISDestroy(compIS, ierr))
         end if
         ssxs = ssxs + ssSize(set)
         sscs = sscs + bs * ssLocalSize
      end do
      deallocate (ssID)
      deallocate (ssSize)
   end subroutine MEF90EXOVecViewSide_Private

#undef __FUNCT__
#define __FUNCT__ "MEF90EXOVecLoadSide_Private"
   subroutine MEF90EXOVecLoadSide_Private(v, exoid, step, offset, ierr)
      integer, intent(IN)               :: exoid
      PetscInt, intent(IN)              :: step, offset
      type(tVec), intent(IN)            :: v
      PetscErrorCode, intent(INOUT)     :: ierr

      PetscInt                         :: xs, xe, bs, c, numSS, set, ssLocalSize, ssxs, sscs
      PetscScalar, dimension(:), pointer :: varray
      PetscInt, dimension(:), pointer    :: ssID, ssSize
      type(tVec)                       :: vComp
      type(tIS)                        :: compIS
      PetscMPIInt                      :: rank

      ssxs = 0_ki
      sscs = 0_ki
      PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr))
      numSS = exinqi(exoid, EX_INQ_SIDE_SETS)
      allocate (ssID(numSS))
      allocate (ssSize(numSS))
      call exgssi(exoid, ssID, ierr)
      do set = 1, numSS
         call exgsp(exoid, ssID(set), ssSize(set), PETSC_NULL_INTEGER, ierr)
      end do
      PetscCall(VecGetOwnershipRange(v, xs, xe, ierr))
      PetscCall(VecGetBlockSize(v, bs, ierr))
      do set = 1, numSS
         !  range of indices for set setID[set]: csxs:csxs + csSize[set]-1
         !  local slice of zonal values:         xs/bs,xm/bs-1
         !  intersection:                        max(xs/bs,csxs),min(xm/bs-1,csxs + csSize[set]-1)
         ssLocalSize = max(0, min(xe / bs, ssxs + ssSize(set)) - max(xs / bs, ssxs))
         if (bs == 1) then
            PetscCall(VecGetArray(v, varray, ierr))
            PetscCall(exgssv(exoid, step, offset, ssID(set), ssLocalSize, varray, ierr))
            PetscCall(VecRestoreArray(v, varray, ierr))
         else
            PetscCall(ISCreateStride(PETSC_COMM_WORLD, ssLocalSize, xs + sscs, bs, compIS, ierr))
            do c = 0, bs - 1
               PetscCall(ISStrideSetStride(compIS, ssLocalSize, xs + sscs + c, bs, ierr))
               PetscCall(VecGetSubVector(v, compIS, vComp, ierr))
               PetscCall(VecGetArray(vComp, varray, ierr))
               call exgpv(exoid, step, EX_SIDE_SET, offset + c, ssID(set), max(xs / bs - ssxs, 0) + 1, ssLocalSize, varray, ierr)
               PetscCall(VecRestoreArray(vComp, varray, ierr))
               PetscCall(VecISCopy(v, compIS, SCATTER_FORWARD, vComp, ierr))
               PetscCall(VecRestoreSubVector(v, compIS, vComp, ierr))
            end do
            PetscCall(ISDestroy(compIS, ierr))
         end if
         ssxs = ssxs + ssSize(set)
         sscs = sscs + bs * ssLocalSize
      end do
      deallocate (ssID)
      deallocate (ssSize)
   end subroutine MEF90EXOVecLoadSide_Private
end module m_MEF90_EXO
