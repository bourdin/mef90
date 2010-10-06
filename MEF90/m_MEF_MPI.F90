Module m_MEF_MPI
#include "finclude/petscdef.h"

   Use m_MEF_Parameters
   Use m_MEF_Types
   Use petsc
   Implicit None
   Private

   Public :: MEF90_Initialize
   Public :: MEF90_Finalize

   Integer, Public                      :: Vect2D_MPIType
   Integer, Public                      :: Vect3D_MPIType
   
   Integer, Public                      :: Mat2D_MPIType
   Integer, Public                      :: Mat3D_MPIType
   Integer, Public                      :: MatS2D_MPIType
   Integer, Public                      :: MatS3D_MPIType
   
   Integer, Public                      :: Tens4OS2D_MPIType
   Integer, Public                      :: Tens4OS3D_MPIType
   
 Contains
   Subroutine MEF90_Initialize()
       PetscInt                          :: iErr
       
      Call PetscInitialize(PETSC_NULL_CHARACTER, iErr); CHKERRQ(iErr)
      Call MPIType_Initialize()
      
      Call MPI_COMM_RANK(MPI_COMM_WORLD, MEF90_MyRank, iErr)
      Call MPI_COMM_SIZE(MPI_COMM_WORLD, MEF90_NumProcs, iErr)
   End Subroutine MEF90_Initialize
   
   Subroutine MEF90_Finalize()
      PetscInt                          :: iErr
      
      Call MPIType_Finalize()
      Call PetscFinalize(iErr)
   End Subroutine MEF90_Finalize

   Subroutine MPIType_Initialize()
      PetscInt, Dimension(:), Pointer   :: BlkCounts, Offsets, DataTypes
      PetscInt                          :: NumBlk, iErr
      
      !!! Vect2D, Vect3D, Mat2D, MatS2D, Mat3D, MatS3D, Tens4OS2D, Tens4OSD3D
      NumBlk=1
      Allocate(BlkCounts(0:NumBlk-1))
      Allocate(Offsets(0:NumBlk-1))
      Allocate(DataTypes(0:NumBlk-1))
      
      
      Offsets(0)   = 0
      DataTypes(0) = MPIU_SCALAR
      
      BlkCounts(0) = 2
      Call MPI_TYPE_STRUCT(1, BlkCounts, Offsets, DataTypes, Vect2D_MPIType, iErr)
      Call MPI_TYPE_COMMIT(Vect2D_MPIType, iErr)
      
      BlkCounts(0) = 3
      Call MPI_TYPE_STRUCT(1, BlkCounts, Offsets, DataTypes, Vect3D_MPIType, iErr)
      Call MPI_TYPE_COMMIT(Vect3D_MPIType, iErr)
      
      
      BlkCounts(0) = 4
      Call MPI_TYPE_STRUCT(1, BlkCounts, Offsets, DataTypes, Mat2D_MPIType, iErr)
      Call MPI_TYPE_COMMIT(Mat2D_MPIType, iErr)
      
      BlkCounts(0) = 3
      Call MPI_TYPE_STRUCT(1, BlkCounts, Offsets, DataTypes, MatS2D_MPIType, iErr)
      Call MPI_TYPE_COMMIT(MatS2D_MPIType, iErr)
      
      BlkCounts(0) = 9
      Call MPI_TYPE_STRUCT(1, BlkCounts, Offsets, DataTypes, Mat3D_MPIType, iErr)
      Call MPI_TYPE_COMMIT(Mat3D_MPIType, iErr)
      
      BlkCounts(0) = 6
      Call MPI_TYPE_STRUCT(1, BlkCounts, Offsets, DataTypes, MatS3D_MPIType, iErr)
      Call MPI_TYPE_COMMIT(MatS3D_MPIType, iErr)
      
      BlkCounts(0) = 6
      Call MPI_TYPE_STRUCT(1, BlkCounts, Offsets, DataTypes, Tens4OS2D_MPIType, iErr)
      Call MPI_TYPE_COMMIT(Tens4OS2D_MPIType, iErr)
      
      BlkCounts(0) = 21
      Call MPI_TYPE_STRUCT(1, BlkCounts, Offsets, DataTypes, Tens4OS3D_MPIType, iErr)
      Call MPI_TYPE_COMMIT(Tens4OS3D_MPIType, iErr)
      DeAllocate(BlkCounts, Offsets, DataTypes)
   End Subroutine MPIType_Initialize
   
   Subroutine MPIType_Finalize()
      PetscInt                          :: iErr
      
      Call MPI_TYPE_FREE(Vect2D_MPIType, iErr)
      Call MPI_TYPE_FREE(Vect3D_MPIType, iErr)
      Call MPI_TYPE_FREE(Mat2D_MPIType, iErr)
      Call MPI_TYPE_FREE(MatS2D_MPIType, iErr)
      Call MPI_TYPE_FREE(Mat3D_MPIType, iErr)
      Call MPI_TYPE_FREE(MatS3D_MPIType, iErr)
      Call MPI_TYPE_FREE(Tens4OS2D_MPIType, iErr)
      Call MPI_TYPE_FREE(Tens4OS3D_MPIType, iErr)
   End Subroutine MPIType_Finalize
End Module m_MEF_MPI
