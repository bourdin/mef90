Module m_MEF90_MPI
#include "petsc/finclude/petsc.h"
   Use m_MEF90_Parameters
   Use petsc
   Implicit None

   Private
   PetscMPIInt,Public,protected   :: Vect2D_MPIType
   PetscMPIInt,Public,protected   :: Vect3D_MPIType
   
   PetscMPIInt,Public,protected   :: Mat2D_MPIType
   PetscMPIInt,Public,protected   :: Mat3D_MPIType
   PetscMPIInt,Public,protected   :: MatS2D_MPIType
   PetscMPIInt,Public,protected   :: MatS3D_MPIType
   
   PetscMPIInt,Public,protected   :: Tens4OS2D_MPIType
   PetscMPIInt,Public,protected   :: Tens4OS3D_MPIType
   
   !PetscInt,Public,protected  :: MEF90_MyRank
   !PetscInt,Public,protected  :: MEF90_NumProcs

   
   Public   :: MEF90MPIInitialize_Private
   Public   :: MEF90MPIFinalize_Private
   
Contains

!!!
!!!  
!!!  MEF90MPIInitialize_Private:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90MPIInitialize_Private"
   Subroutine MEF90MPIInitialize_Private(ierr)
      PetscErrorCode,Intent(INOUT)          :: ierr
      PetscMPIInt,Dimension(:),Pointer      :: BlkCounts,Offsets,DataTypes
      PetscMPIInt                           :: NumBlk
      !!! Vect2D,Vect3D,Mat2D,MatS2D,Mat3D,MatS3D,Tens4OS2D,Tens4OSD3D
      NumBlk=1
      Allocate(BlkCounts(0:NumBlk-1))
      Allocate(Offsets(0:NumBlk-1))
      Allocate(DataTypes(0:NumBlk-1))
      
      
      Offsets(0)   = 0
      DataTypes(0) = MPIU_SCALAR
      
      BlkCounts(0) = 2
      PetscCallMPI(MPI_TYPE_STRUCT(1,BlkCounts,Offsets,DataTypes,Vect2D_MPIType,ierr))
      PetscCallMPI(MPI_TYPE_COMMIT(Vect2D_MPIType,ierr))
      
      BlkCounts(0) = 3
      PetscCallMPI(MPI_TYPE_STRUCT(1,BlkCounts,Offsets,DataTypes,Vect3D_MPIType,ierr))
      PetscCallMPI(MPI_TYPE_COMMIT(Vect3D_MPIType,ierr))
      
      
      BlkCounts(0) = 4
      PetscCallMPI(MPI_TYPE_STRUCT(1,BlkCounts,Offsets,DataTypes,Mat2D_MPIType,ierr))
      PetscCallMPI(MPI_TYPE_COMMIT(Mat2D_MPIType,ierr))
      
      BlkCounts(0) = 3
      PetscCallMPI(MPI_TYPE_STRUCT(1,BlkCounts,Offsets,DataTypes,MatS2D_MPIType,ierr))
      PetscCallMPI(MPI_TYPE_COMMIT(MatS2D_MPIType,ierr))
      
      BlkCounts(0) = 9
      PetscCallMPI(MPI_TYPE_STRUCT(1,BlkCounts,Offsets,DataTypes,Mat3D_MPIType,ierr))
      PetscCallMPI(MPI_TYPE_COMMIT(Mat3D_MPIType,ierr))
      
      BlkCounts(0) = 6
      PetscCallMPI(MPI_TYPE_STRUCT(1,BlkCounts,Offsets,DataTypes,MatS3D_MPIType,ierr))
      PetscCallMPI(MPI_TYPE_COMMIT(MatS3D_MPIType,ierr))
      
      BlkCounts(0) = 6
      PetscCallMPI(MPI_TYPE_STRUCT(1,BlkCounts,Offsets,DataTypes,Tens4OS2D_MPIType,ierr))
      PetscCallMPI(MPI_TYPE_COMMIT(Tens4OS2D_MPIType,ierr))
      
      BlkCounts(0) = 21
      PetscCallMPI(MPI_TYPE_STRUCT(1,BlkCounts,Offsets,DataTypes,Tens4OS3D_MPIType,ierr))
      PetscCallMPI(MPI_TYPE_COMMIT(Tens4OS3D_MPIType,ierr))
      DeAllocate(BlkCounts,Offsets,DataTypes)
      ierr = 0
   End Subroutine MEF90MPIInitialize_Private
   
!!!
!!!  
!!!  MEF90MPIFinalize_Private:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90MPIFinalize_Private"
   Subroutine MEF90MPIFinalize_Private(ierr)
      PetscMPIInt,Intent(INOUT)             :: ierr
      
      PetscCallMPI(MPI_TYPE_FREE(Vect2D_MPIType,ierr))
      PetscCallMPI(MPI_TYPE_FREE(Vect3D_MPIType,ierr))
      PetscCallMPI(MPI_TYPE_FREE(Mat2D_MPIType,ierr))
      PetscCallMPI(MPI_TYPE_FREE(MatS2D_MPIType,ierr))
      PetscCallMPI(MPI_TYPE_FREE(Mat3D_MPIType,ierr))
      PetscCallMPI(MPI_TYPE_FREE(MatS3D_MPIType,ierr))
      PetscCallMPI(MPI_TYPE_FREE(Tens4OS2D_MPIType,ierr))
      PetscCallMPI(MPI_TYPE_FREE(Tens4OS3D_MPIType,ierr))
      ierr = 0
   End Subroutine MEF90MPIFinalize_Private
End Module m_MEF90_MPI
