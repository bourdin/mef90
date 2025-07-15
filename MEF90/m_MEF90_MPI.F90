module m_MEF90_MPI
#include "petsc/finclude/petsc.h"
   use m_MEF90_Parameters
   implicit none(type, external)

   private
   PetscMPIInt, public, protected   :: Vect2D_MPIType
   PetscMPIInt, public, protected   :: Vect3D_MPIType

   PetscMPIInt, public, protected   :: Mat2D_MPIType
   PetscMPIInt, public, protected   :: Mat3D_MPIType
   PetscMPIInt, public, protected   :: MatS2D_MPIType
   PetscMPIInt, public, protected   :: MatS3D_MPIType

   PetscMPIInt, public, protected   :: Tens4OS2D_MPIType
   PetscMPIInt, public, protected   :: Tens4OS3D_MPIType

   !PetscInt,Public,protected  :: MEF90_MyRank
   !PetscInt,Public,protected  :: MEF90_NumProcs

   public   :: MEF90MPIInitialize_Private
   public   :: MEF90MPIFinalize_Private

contains

!!!
!!!
!!!  MEF90MPIInitialize_Private:
!!!
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90MPIInitialize_Private"
   subroutine MEF90MPIInitialize_Private(ierr)
      PetscErrorCode, intent(INOUT)          :: ierr
      PetscMPIInt, dimension(:), pointer      :: BlkCounts, Offsets, DataTypes
      PetscMPIInt                           :: NumBlk
      !!! Vect2D,Vect3D,Mat2D,MatS2D,Mat3D,MatS3D,Tens4OS2D,Tens4OSD3D
      NumBlk = 1
      allocate (BlkCounts(0:NumBlk - 1))
      allocate (Offsets(0:NumBlk - 1))
      allocate (DataTypes(0:NumBlk - 1))

      Offsets(0) = 0
      DataTypes(0) = MPIU_SCALAR

      BlkCounts(0) = 2
      PetscCallMPI(MPI_TYPE_STRUCT(1, BlkCounts, Offsets, DataTypes, Vect2D_MPIType, ierr))
      PetscCallMPI(MPI_TYPE_COMMIT(Vect2D_MPIType, ierr))

      BlkCounts(0) = 3
      PetscCallMPI(MPI_TYPE_STRUCT(1, BlkCounts, Offsets, DataTypes, Vect3D_MPIType, ierr))
      PetscCallMPI(MPI_TYPE_COMMIT(Vect3D_MPIType, ierr))

      BlkCounts(0) = 4
      PetscCallMPI(MPI_TYPE_STRUCT(1, BlkCounts, Offsets, DataTypes, Mat2D_MPIType, ierr))
      PetscCallMPI(MPI_TYPE_COMMIT(Mat2D_MPIType, ierr))

      BlkCounts(0) = 3
      PetscCallMPI(MPI_TYPE_STRUCT(1, BlkCounts, Offsets, DataTypes, MatS2D_MPIType, ierr))
      PetscCallMPI(MPI_TYPE_COMMIT(MatS2D_MPIType, ierr))

      BlkCounts(0) = 9
      PetscCallMPI(MPI_TYPE_STRUCT(1, BlkCounts, Offsets, DataTypes, Mat3D_MPIType, ierr))
      PetscCallMPI(MPI_TYPE_COMMIT(Mat3D_MPIType, ierr))

      BlkCounts(0) = 6
      PetscCallMPI(MPI_TYPE_STRUCT(1, BlkCounts, Offsets, DataTypes, MatS3D_MPIType, ierr))
      PetscCallMPI(MPI_TYPE_COMMIT(MatS3D_MPIType, ierr))

      BlkCounts(0) = 6
      PetscCallMPI(MPI_TYPE_STRUCT(1, BlkCounts, Offsets, DataTypes, Tens4OS2D_MPIType, ierr))
      PetscCallMPI(MPI_TYPE_COMMIT(Tens4OS2D_MPIType, ierr))

      BlkCounts(0) = 21
      PetscCallMPI(MPI_TYPE_STRUCT(1, BlkCounts, Offsets, DataTypes, Tens4OS3D_MPIType, ierr))
      PetscCallMPI(MPI_TYPE_COMMIT(Tens4OS3D_MPIType, ierr))
      deallocate (BlkCounts, Offsets, DataTypes)
      ierr = 0
   end subroutine MEF90MPIInitialize_Private

!!!
!!!
!!!  MEF90MPIFinalize_Private:
!!!
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90MPIFinalize_Private"
   subroutine MEF90MPIFinalize_Private(ierr)
      PetscMPIInt, intent(INOUT)             :: ierr

      PetscCallMPI(MPI_TYPE_FREE(Vect2D_MPIType, ierr))
      PetscCallMPI(MPI_TYPE_FREE(Vect3D_MPIType, ierr))
      PetscCallMPI(MPI_TYPE_FREE(Mat2D_MPIType, ierr))
      PetscCallMPI(MPI_TYPE_FREE(MatS2D_MPIType, ierr))
      PetscCallMPI(MPI_TYPE_FREE(Mat3D_MPIType, ierr))
      PetscCallMPI(MPI_TYPE_FREE(MatS3D_MPIType, ierr))
      PetscCallMPI(MPI_TYPE_FREE(Tens4OS2D_MPIType, ierr))
      PetscCallMPI(MPI_TYPE_FREE(Tens4OS3D_MPIType, ierr))
      ierr = 0
   end subroutine MEF90MPIFinalize_Private
end module m_MEF90_MPI
