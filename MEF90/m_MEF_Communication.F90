Module m_MEF_Communication
  Use m_MEF_Types
  Use m_MEF_EXO
!  Use MPI

  IMPLICIT NONE
  Private

#include "include/finclude/petsc.h"

  Public :: Broadcast_Geom
  Public :: Show_One_EXO_Geom_Info
  Public :: Show_All_EXO_Geom_Info
  Public :: MEF90_Initialize
  Public :: MEF90_Finalize
  Public :: Broadcast_Node
  Public :: Broadcast_Element

  Interface Broadcast_Node
     Module Procedure Broadcast_Node1D, Broadcast_Node2D, Broadcast_Node3D
  End Interface
  
  Interface Broadcast_Element
     Module Procedure Broadcast_Element2D_Scal, Broadcast_Element2D,          &
          & Broadcast_Element2D_Elast, Broadcast_Element3D_Scal,              &
          & Broadcast_Element3D, Broadcast_Element3D_Elast
  End Interface
  

  Integer, Public                                   :: Vect2D_MPIType
  Integer, Public                                   :: Vect3D_MPIType

  Integer, Public                                   :: Mat2D_MPIType
  Integer, Public                                   :: Mat3D_MPIType
  Integer, Public                                   :: MatS2D_MPIType
  Integer, Public                                   :: MatS3D_MPIType

  Integer, Public                                   :: Tens4OS2D_MPIType

  Integer, Public                                   :: Node1D_MPIType
  Integer, Public                                   :: Node2D_MPIType
  Integer, Public                                   :: Node3D_MPIType

  !!! Elements contain pointers, therefore it is not possible to create 
  !!! MPITypes for them...

  Integer, Parameter, Public                        :: MPI_MaxBuff = 2**19
  ! Slice of the largest vector we try to send in one call to MPI_BCAST

  Integer, Public                                   :: MEF90_MyRank
  Integer, Public                                   :: MEF90_NumProcs
  
Contains
  Subroutine MEF90_Initialize()
    Integer                           :: iErr
    
    Call PetscInitialize(PETSC_NULL_CHARACTER, iErr)
    Call MPIType_Initialize()
    Call MPI_BARRIER(MPI_COMM_WORLD, iErr)

    Call MPI_COMM_RANK(MPI_COMM_WORLD, MEF90_MyRank, iErr)
    Call MPI_COMM_SIZE(MPI_COMM_WORLD, MEF90_NumProcs, iErr)
  End Subroutine MEF90_Initialize
  
  Subroutine MEF90_Finalize()
    Integer                           :: iErr
    
    Call PetscFinalize(PETSC_NULL_CHARACTER, iErr)

  End Subroutine MEF90_Finalize

  Subroutine MPIType_Initialize()
    Integer, Dimension(:), Pointer    :: BlkCounts, Offsets, DataTypes

    Integer                           :: NumBlk, extent, iErr

    !!! Vect2D, Vect3D, Mat2D, MatS2D, Mat3D, MatS3D, Tens4OS2D
    NumBlk=1
    Allocate(BlkCounts(0:NumBlk-1))
    Allocate(Offsets(0:NumBlk-1))
    Allocate(DataTypes(0:NumBlk-1))


    Offsets(0)   = 0
    DataTypes(0) = MPI_DOUBLE_PRECISION

    BlkCounts(0) = 2
    Call MPI_TYPE_STRUCT(1, BlkCounts, Offsets, DataTypes, Vect2D_MPIType,    &
         & iErr)
    Call MPI_TYPE_COMMIT(Vect2D_MPIType, iErr)

    BlkCounts(0) = 3
    Call MPI_TYPE_STRUCT(1, BlkCounts, Offsets, DataTypes, Vect3D_MPIType,    &
         & iErr)
    Call MPI_TYPE_COMMIT(Vect3D_MPIType, iErr)


    BlkCounts(0) = 4
    Call MPI_TYPE_STRUCT(1, BlkCounts, Offsets, DataTypes, Mat2D_MPIType,     &
         & iErr)
    Call MPI_TYPE_COMMIT(Mat2D_MPIType, iErr)

    BlkCounts(0) = 3
    Call MPI_TYPE_STRUCT(1, BlkCounts, Offsets, DataTypes, MatS2D_MPIType,    &
         & iErr)
    Call MPI_TYPE_COMMIT(MatS2D_MPIType, iErr)

    BlkCounts(0) = 9
    Call MPI_TYPE_STRUCT(1, BlkCounts, Offsets, DataTypes, Mat3D_MPIType,     &
         & iErr)
    Call MPI_TYPE_COMMIT(Mat3D_MPIType, iErr)

    BlkCounts(0) = 6
    Call MPI_TYPE_STRUCT(1, BlkCounts, Offsets, DataTypes, MatS3D_MPIType,    &
         & iErr)
    Call MPI_TYPE_COMMIT(MatS2D_MPIType, iErr)

    BlkCounts(0) = 6
    Call MPI_TYPE_STRUCT(1, BlkCounts, Offsets, DataTypes, Tens4OS2D_MPIType, &
         & iErr)
    Call MPI_TYPE_COMMIT(Tens4OS2D_MPIType, iErr)
    
    DeAllocate(BlkCounts, Offsets, DataTypes)

    ! Node1D, Node2D, Node3D
    NumBlk=2
    Allocate(BlkCounts(0:NumBlk-1))
    Allocate(Offsets(0:NumBlk-1))
    Allocate(DataTypes(0:NumBlk-1))



    Call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, extent, iErr)
    BlkCounts(0) = 1
    Offsets(0)   = 0
    DataTypes(0) = MPI_DOUBLE_PRECISION
    BlkCounts(1) = 2
    Offsets(1)   = extent
    DataTypes(1) = MPI_INTEGER

    Call MPI_TYPE_STRUCT(2, BlkCounts, Offsets, DataTypes, Node1D_MPIType,    &
         & iErr)
    Call MPI_TYPE_COMMIT(Node1D_MPIType, iErr)
    

    Call MPI_TYPE_EXTENT(Vect2D_MPIType, extent, iErr)
    BlkCounts(0) = 1
    Offsets(0)   = 0
    DataTypes(0) = Vect2D_MPIType
    BlkCounts(1) = 2
    Offsets(1)   = extent
    DataTypes(1) = MPI_INTEGER

    Call MPI_TYPE_STRUCT(2, BlkCounts, Offsets, DataTypes, Node2D_MPIType,    &
         & iErr)
    Call MPI_TYPE_COMMIT(Node2D_MPIType, iErr)
    
    Call MPI_TYPE_EXTENT(Vect3D_MPIType, extent, iErr)
    BlkCounts(0) = 1
    Offsets(0)   = 0
    DataTypes(0) = Vect3D_MPIType
    BlkCounts(1) = 2
    Offsets(1)   = extent
    DataTypes(1) = MPI_INTEGER

    Call MPI_TYPE_STRUCT(2, BlkCounts, Offsets, DataTypes, Node3D_MPIType,    &
         & iErr)
    Call MPI_TYPE_COMMIT(Node3D_MPIType, iErr)
    
    DeAllocate(BlkCounts, Offsets, DataTypes)    
  End Subroutine MPIType_Initialize
  
  Subroutine MPIType_Finalize()
    Integer                           :: iErr

    Call MPI_TYPE_FREE(Vect2D_MPIType, iErr)
    Call MPI_TYPE_FREE(Vect3D_MPIType, iErr)
    Call MPI_TYPE_FREE(Mat2D_MPIType, iErr)
    Call MPI_TYPE_FREE(MatS2D_MPIType, iErr)
    Call MPI_TYPE_FREE(Mat3D_MPIType, iErr)
    Call MPI_TYPE_FREE(MatS3D_MPIType, iErr)
    Call MPI_TYPE_FREE(Tens4OS2D_MPIType, iErr)

    Call MPI_TYPE_FREE(Node1D_MPIType, iErr)
    Call MPI_TYPE_FREE(Node2D_MPIType, iErr)
    Call MPI_TYPE_FREE(Node3D_MPIType, iErr)
  End Subroutine MPIType_Finalize


  Subroutine Broadcast_Geom(Geom, Source, Communicator)
    Integer, Intent(IN)                              :: Source
    Type (EXO_Geom_Info), Intent(INOUT)              :: Geom
    Integer, Intent(IN), Optional                    :: Communicator

    Integer                                          :: COMM
    Integer                                          :: MyRank, iErr
    Integer                                          :: iBlk, iSet, iQA, i

!    Integer tag, i, j, k, iSet, RequestS, RequestR, iBlk, MPI_Size  

    If ( Present(Communicator) ) Then
       Comm = Communicator
    Else
       Comm = MPI_COMM_WORLD
    End If

    Call MPI_COMM_RANK(Comm, MyRank, iErr)
    If (MyRank /= Source) Then
       Call Destroy_EXO_Geom_Info(Geom)
    End If

    Call MPI_BCAST(Geom%filename, MXLNLN, MPI_CHARACTER, Source, Comm, ierr)
   
    Call MPI_BCAST(Geom%title, MXLNLN, MPI_CHARACTER, Source, Comm, ierr)
    Call MPI_BCAST(Geom%num_dim, 1, MPI_INTEGER, Source, Comm, iErr)
    Call MPI_BCAST(Geom%num_nodes, 1, MPI_INTEGER, Source, Comm, iErr)
    Call MPI_BCAST(Geom%num_elems, 1, MPI_INTEGER, Source, Comm, iErr)
    Call MPI_BCAST(Geom%Num_Elem_Blks, 1, MPI_INTEGER, Source, Comm, iErr)
    Call MPI_BCAST(Geom%Num_Node_Sets, 1, MPI_INTEGER, Source, Comm, iErr)
    Call MPI_BCAST(Geom%Num_Side_Sets, 1, MPI_INTEGER, Source, Comm, iErr)

    ! Broadcasting Elem Block
    If (MyRank /= Source) Then
       Allocate(Geom%Elem_Blk(Geom%Num_Elem_blks))
    End If
    Do iBlk = 1, Geom%Num_Elem_Blks
       Call MPI_BCAST(Geom%Elem_Blk(iBlk)%ID, 1, MPI_INTEGER, Source, Comm,   &
            & iErr)
       Call MPI_BCAST(Geom%Elem_Blk(iBlk)%Type, MXLNLN, MPI_CHARACTER,        &
            & Source, Comm, iErr)
       Call MPI_BCAST(Geom%Elem_Blk(iBlk)%Num_Elems, 1, MPI_INTEGER, Source,  &
            & Comm, iErr)
       Call MPI_BCAST(Geom%Elem_Blk(iBlk)%Num_Nodes_Per_Elem, 1, MPI_INTEGER, &
            & Source, Comm, iErr)
       Call MPI_BCAST(Geom%Elem_Blk(iBlk)%Num_Attr, 1, MPI_INTEGER, Source,   &
            & Comm, iErr)

       If (MyRank /= Source) Then
          Allocate( Geom%Elem_Blk(iBlk)%Elem_ID(                              &
               & Geom%Elem_Blk(iBlk)%Num_Elems) )
       End If
       Call MPI_BCAST(Geom%Elem_Blk(iBlk)%Elem_ID,                            &
            & Geom%elem_blk(iBlk)%Num_Elems, MPI_INTEGER, Source, Comm, iErr)
    End Do
    
    ! Broadcasting Node Sets
    If (MyRank /= Source) Then
       Allocate(Geom%Node_Set(Geom%Num_Node_Sets))
    End If
    Do iSet = 1, Geom%Num_Node_Sets
       Call MPI_BCAST(Geom%Node_Set(iSet)%ID, 1, MPI_INTEGER, Source, Comm,   &
            & iErr)
       Call MPI_BCAST(Geom%Node_Set(iSet)%Num_Nodes, 1, MPI_INTEGER, Source,  &
            & Comm, iErr)
       Call MPI_BCAST(Geom%Node_Set(iSet)%Num_Dist_Factors, 1, MPI_INTEGER,   &
            & Source, Comm, iErr)

       If (MyRank /= Source) Then
          Allocate( Geom%Node_Set(iSet)%Node_ID(                              &
               & Geom%Node_Set(iSet)%Num_Nodes) )
          Allocate( Geom%Node_Set(iSet)%Dist_Factor(                          &
               & Geom%Node_Set(iSet)%Num_Nodes) )
       End If
       Call MPI_BCAST(Geom%Node_Set(iSet)%Node_ID,                            &
            & Geom%Node_Set(iSet)%Num_Nodes, MPI_INTEGER, Source, Comm, iErr)
       Call MPI_BCAST(Geom%Node_Set(iSet)%Dist_Factor,                        &
            & Geom%Node_Set(iSet)%Num_Dist_Factors, MPI_DOUBLE_PRECISION,     &
            & Source, Comm, iErr)
    End Do
    
    ! Broadcasting QA
    Call MPI_BCAST(Geom%Num_QA, 1, MPI_INTEGER, Source, Comm, iErr)

    If (MyRank /= Source) Then
       Allocate (Geom%QA_Rec(4, Geom%Num_QA))
    End If
    Do iQA = 1, Geom%Num_QA
       Do i = 1, 4
          Call MPI_BCAST(Geom%QA_rec(i, iQA), MXSTLN, MPI_CHARACTER, Source,  &
               & Comm, iErr)
       End Do
    End Do
    
    Call MPI_BCAST(Geom%Numbering, 1, MPI_INTEGER, Source, Comm, iErr)
  End Subroutine Broadcast_Geom

  Subroutine Broadcast_Node1D(Node_db, Source, Communicator)
    Type(Node1D), Dimension(:), Pointer   :: Node_db
    Integer, Intent(IN)                   :: Source
    Integer, Intent(IN), Optional         :: Communicator
    
    Integer                               :: MyRank, iErr, NN, Comm

    If ( Present(Communicator) ) Then
       Comm = Communicator
    Else
       Comm = MPI_COMM_WORLD
    End If

    Call MPI_COMM_RANK(Comm, MyRank, iErr)
    If ( ( MyRank /= Source ) .AND. ( Associated(Node_db) ) ) Then
       DeAllocate (Node_db)
    End If

    If (MyRank == Source) Then
       NN = Size(Node_db)
    End If
    Call MPI_BCAST(NN, 1, MPI_INTEGER, Source, Comm, iErr)
    
    If (MyRank /= Source) Then
       Allocate (Node_db(NN))
    End If
    Call MPI_BCAST(Node_db, NN, Node1D_MPIType, Source, Comm, iErr)

  End Subroutine Broadcast_Node1D

  Subroutine Broadcast_Node2D(Node_db, Source, Communicator)
    Type(Node2D), Dimension(:), Pointer   :: Node_db
    Integer, Intent(IN)                   :: Source
    Integer, Intent(IN), Optional         :: Communicator
    
    Integer                               :: MyRank, iErr, NN, Comm

    If ( Present(Communicator) ) Then
       Comm = Communicator
    Else
       Comm = MPI_COMM_WORLD
    End If

    Call MPI_COMM_RANK(Comm, MyRank, iErr)
    If ( ( MyRank /= Source ) .AND. ( Associated(Node_db) ) ) Then
       DeAllocate (Node_db)
    End If

    If (MyRank == Source) Then
       NN = Size(Node_db)
    End If
    Call MPI_BCAST(NN, 1, MPI_INTEGER, Source, Comm, iErr)
    
    If (MyRank /= Source) Then
       Allocate (Node_db(NN))
    End If
    Call MPI_BCAST(Node_db, NN, Node2D_MPIType, Source, Comm, iErr)

  End Subroutine Broadcast_Node2D

  Subroutine Broadcast_Node3D(Node_db, Source, Communicator)
    Type(Node3D), Dimension(:), Pointer   :: Node_db
    Integer, Intent(IN)                   :: Source
    Integer, Intent(IN), Optional         :: Communicator
    
    Integer                               :: MyRank, iErr, NN, Comm

    If ( Present(Communicator) ) Then
       Comm = Communicator
    Else
       Comm = MPI_COMM_WORLD
    End If

    Call MPI_COMM_RANK(Comm, MyRank, iErr)
    If ( ( MyRank /= Source ) .AND. ( Associated(Node_db) ) ) Then
       DeAllocate (Node_db)
    End If

    If (MyRank == Source) Then
       NN = Size(Node_db)
    End If
    Call MPI_BCAST(NN, 1, MPI_INTEGER, Source, Comm, iErr)
    If (MyRank /= Source) Then
       Allocate (Node_db(NN))
    End If
    Call MPI_BCAST(Node_db, NN, Node3D_MPIType, Source, Comm, iErr)
  End Subroutine Broadcast_Node3D

!!!
!!! Broadcast element stuff
!!!


  Subroutine Broadcast_Element2D_Scal(Elem_db, Geom, Source, Communicator)
    Type(Element2D_Scal), Dimension(:), Pointer     :: Elem_db
    Type (EXO_Geom_Info), Intent(IN)                :: Geom
    Integer, Intent(IN)                             :: Source
    Integer, Intent(IN), Optional                   :: Communicator

    Integer                                         :: Comm, MyRank, iErr
    Integer                                         :: nE, eNum, iE, iBlk

    Integer                                         :: iSlice, SliceSize
    Integer                                         :: SliceIdxMin, SliceIdxMax
    Integer, Dimension(:), Pointer                  :: Buffer1
    Integer, Dimension(:,:), Pointer                :: Buffer2


    ! WARNING: We don't check if the elem_db are allocated or not
    ! We assume that Elem_db is allocated on source, but not on the 
    ! other cpu in te communicator.
    ! If this is not true, this subroutine is going to crash...

    ! We do NOT Broadcast the Gauss points informations. 
    ! It just does not make any sense to do so.

    ! Use the geom data and send element block after element block,
    ! since all element in a block must be of the same type

    ! We Assume that the Geom database has been broadcast, already

    If ( Present(Communicator) ) Then
       Comm = Communicator
    Else
       Comm = MPI_COMM_WORLD
    End If

    Call MPI_COMM_RANK(Comm, MyRank, iErr)
    If (MyRank /= Source) Then
       Allocate(Elem_db(Geom%Num_Elems))
    EndIf

    Do iSlice = 1, Geom%Num_Elems / MPI_MaxBuff + 1
       SliceIdxMin = (iSlice-1) * MPI_MaxBuff + 1
       SliceIdxMax = Min(Geom%Num_Elems, iSlice * MPI_MaxBuff)
       SliceSize   = SliceIdxMax - SliceIdxMin + 1
       Allocate (Buffer1(SliceSize))
       
       If (MyRank == Source) Then
          Buffer1 = Elem_db(SliceIdxMin:SliceIdxMax)%NB_DoF
       End If
       Call MPI_BCAST(Buffer1, SliceSize, MPI_INTEGER, Source, Comm, iErr)
       If (MyRank /= Source) Then
          Elem_db(SliceIdxMin:SliceIdxMax)%NB_DoF = Buffer1
       End If
       
       If (MyRank == Source) Then
          Buffer1 = Elem_db(SliceIdxMin:SliceIdxMax)%NB_Gauss
       End If
       Call MPI_BCAST(Buffer1, SliceSize, MPI_INTEGER, Source, Comm, iErr)
       If (MyRank /= Source) Then
          Elem_db(SliceIdxMin:SliceIdxMax)%NB_Gauss = Buffer1
       End If
       
       If (MyRank == Source) Then
          Buffer1 = Elem_db(SliceIdxMin:SliceIdxMax)%ID_EL
       End If
       Call MPI_BCAST(Buffer1, SliceSize, MPI_INTEGER, Source, Comm, iErr)
       If (MyRank /= Source) Then
          Elem_db(SliceIdxMin:SliceIdxMax)%ID_EL = Buffer1
       End If
       
       DeAllocate(Buffer1)
    End Do

    Do iBlk = 1, Geom%Num_Elem_Blks
       Do iSlice = 1, Geom%Elem_Blk(iBlk)%Num_Elems  *                        &
            & Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem / MPI_MaxBuff + 1
          SliceIdxMin = (iSlice-1) * MPI_MaxBuff /                            &
               & Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem + 1
          SliceIdxMax = Min(Geom%Elem_Blk(iBlk)%Num_Elems,                    &
               & iSlice * MPI_MaxBuff / Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem)
          SliceSize   = SliceIdxMax - SliceIdxMin + 1
          
          Allocate(Buffer2(SliceIdxMin:SliceIdxMax,                           &
               & Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem))

          If (MyRank == Source) Then
             Do iE = SliceIdxMin, SliceIdxMax
                Buffer2(iE,:) =                                              &
                     & Elem_db(Geom%Elem_Blk(iBlk)%Elem_ID(iE))%ID_DoF
             End Do
          End If

          Call MPI_BCast(Buffer2, SliceSize *                                 &
               & Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem, MPI_INTEGER, Source, &
               & Comm, iErr)

          If (MyRank /= Source) Then
             Do iE = SliceIdxMin, SliceIdxMax
                Allocate(Elem_db(Geom%Elem_Blk(iBlk)%Elem_ID(iE))%ID_DoF      &
                     & (Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem))
                Elem_db(Geom%Elem_Blk(iBlk)%Elem_ID(iE))%ID_DoF = Buffer2(iE,:)
             End Do
          End If
          DeAllocate(Buffer2)
       End Do
    End Do
!!$
!!$       If (MyRank == Source) Then
!!$          Do iE = 1, Geom%Elem_Blk(iBlk)%Num_Elems
!!$             Buffer2(iE, :) = Elem_db(Geom%Elem_Blk(iBlk)%Elem_ID(iE))%ID_DoF
!!$          End Do
!!$       End If
!!$
!!$       Call MPI_BCast(Buffer2, Geom%Elem_Blk(iBlk)%Num_Elems *             &
!!$            & Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem, MPI_INTEGER, Source,    &
!!$            & Comm, iErr)
!!$       If (MyRank /= Source) Then
!!$          Do iE = 1, Geom%Elem_Blk(iBlk)%Num_Elems
!!$             Allocate(Elem_db(Geom%Elem_Blk(iBlk)%Elem_ID(iE))%ID_DoF         &
!!$                  & (Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem))
!!$             Elem_db(Geom%Elem_Blk(iBlk)%Elem_ID(iE))%ID_DoF = Buffer2(iE, :)
!!$          End Do
!!$       End If
!!$       DeAllocate(Buffer2)
!!$    End Do
!!$
  End Subroutine Broadcast_Element2D_Scal


  Subroutine Broadcast_Element2D(Elem_db, Geom, Source, Communicator)
    Type(Element2D), Dimension(:), Pointer          :: Elem_db
    Type (EXO_Geom_Info), Intent(IN)                :: Geom
    Integer, Intent(IN)                             :: Source
    Integer, Intent(IN), Optional                   :: Communicator

    Integer                                         :: Comm, MyRank, iErr
    Integer                                         :: nE, eNum, iE, iBlk

    Integer, Dimension(:), Pointer                  :: Buffer1
    Integer, Dimension(:,:), Pointer                :: Buffer2


    ! WARNING: We don't check if the elem_db are allocated or not
    ! We assume that Elem_db is allocated on source, but not on the 
    ! other cpu in te communicator.
    ! If this is not true, this subroutine is going to crash...

    ! We do NOT Broadcast the Gauss points informations. 
    ! It just does not make any sense to do so.

    ! Use the geom data and send element block after element block,
    ! since all element in a block must be of the same type

    ! We Assume that the Geom database has been broadcast, already

    If ( Present(Communicator) ) Then
       Comm = Communicator
    Else
       Comm = MPI_COMM_WORLD
    End If

    Call MPI_COMM_RANK(Comm, MyRank, iErr)
    If (MyRank /= Source) Then
       Allocate(Elem_db(Geom%Num_Elems))
    EndIf

    Allocate (Buffer1(Geom%Num_Elems))

    If (MyRank == Source) Then
       Buffer1 = Elem_db(:)%NB_DoF
    End If
    Call MPI_BCAST(Buffer1, Geom%Num_Elems, MPI_INTEGER, Source, Comm, iErr)
    If (MyRank /= Source) Then
       Elem_db(:)%NB_DoF = Buffer1
    End If

    If (MyRank == Source) Then
       Buffer1 = Elem_db(:)%NB_Gauss
    End If
    Call MPI_BCAST(Buffer1, Geom%Num_Elems, MPI_INTEGER, Source, Comm, iErr)
    If (MyRank /= Source) Then
       Elem_db(:)%NB_Gauss = Buffer1
    End If

    If (MyRank == Source) Then
       Buffer1 = Elem_db(:)%ID_EL
    End If
    Call MPI_BCAST(Buffer1, Geom%Num_Elems, MPI_INTEGER, Source, Comm, iErr)
    If (MyRank /= Source) Then
       Elem_db(:)%ID_EL = Buffer1
    End If

    DeAllocate(Buffer1)

    Do iBlk = 1, Geom%Num_Elem_Blks
       Allocate(Buffer2(Geom%Elem_Blk(iBlk)%Num_Elems,                        &
            & 2 * Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem))

       If (MyRank == Source) Then
          Do iE = 1, Geom%Elem_Blk(iBlk)%Num_Elems
             Buffer2(iE, :) = Elem_db(Geom%Elem_Blk(iBlk)%Elem_ID(iE))%ID_DoF
          End Do
       End If

       Call MPI_BCast(Buffer2, Geom%Elem_Blk(iBlk)%Num_Elems * 2 *            &
            & Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem, MPI_INTEGER, Source,    &
            & Comm, iErr)
       If (MyRank /= Source) Then
          Do iE = 1, Geom%Elem_Blk(iBlk)%Num_Elems
             Allocate(Elem_db(Geom%Elem_Blk(iBlk)%Elem_ID(iE))%ID_DoF         &
                  & (2 * Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem))
             Elem_db(Geom%Elem_Blk(iBlk)%Elem_ID(iE))%ID_DoF = Buffer2(iE, :)
          End Do
       End If
       DeAllocate(Buffer2)
    End Do
  End Subroutine Broadcast_Element2D

  Subroutine Broadcast_Element2D_Elast(Elem_db, Geom, Source, Communicator)
    Type(Element2D_Elast), Dimension(:), Pointer    :: Elem_db
    Type (EXO_Geom_Info), Intent(IN)                :: Geom
    Integer, Intent(IN)                             :: Source
    Integer, Intent(IN), Optional                   :: Communicator

    Integer                                         :: Comm, MyRank, iErr
    Integer                                         :: nE, eNum, iE, iBlk

    Integer, Dimension(:), Pointer                  :: Buffer1
    Integer, Dimension(:,:), Pointer                :: Buffer2


    ! WARNING: We don't check if the elem_db are allocated or not
    ! We assume that Elem_db is allocated on source, but not on the 
    ! other cpu in te communicator.
    ! If this is not true, this subroutine is going to crash...

    ! We do NOT Broadcast the Gauss points informations. 
    ! It just does not make any sense to do so.

    ! Use the geom data and send element block after element block,
    ! since all element in a block must be of the same type

    ! We Assume that the Geom database has been broadcast, already

    If ( Present(Communicator) ) Then
       Comm = Communicator
    Else
       Comm = MPI_COMM_WORLD
    End If

    Call MPI_COMM_RANK(Comm, MyRank, iErr)
    If (MyRank /= Source) Then
       Allocate(Elem_db(Geom%Num_Elems))
    EndIf

    Allocate (Buffer1(Geom%Num_Elems))

    If (MyRank == Source) Then
       Buffer1 = Elem_db(:)%NB_DoF
    End If
    Call MPI_BCAST(Buffer1, Geom%Num_Elems, MPI_INTEGER, Source, Comm, iErr)
    If (MyRank /= Source) Then
       Elem_db(:)%NB_DoF = Buffer1
    End If

    If (MyRank == Source) Then
       Buffer1 = Elem_db(:)%NB_Gauss
    End If
    Call MPI_BCAST(Buffer1, Geom%Num_Elems, MPI_INTEGER, Source, Comm, iErr)
    If (MyRank /= Source) Then
       Elem_db(:)%NB_Gauss = Buffer1
    End If

    If (MyRank == Source) Then
       Buffer1 = Elem_db(:)%ID_EL
    End If
    Call MPI_BCAST(Buffer1, Geom%Num_Elems, MPI_INTEGER, Source, Comm, iErr)
    If (MyRank /= Source) Then
       Elem_db(:)%ID_EL = Buffer1
    End If

    DeAllocate(Buffer1)

    Do iBlk = 1, Geom%Num_Elem_Blks
       Allocate(Buffer2(Geom%Elem_Blk(iBlk)%Num_Elems,                        &
            & 2 * Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem))

       If (MyRank == Source) Then
          Do iE = 1, Geom%Elem_Blk(iBlk)%Num_Elems
             Buffer2(iE, :) = Elem_db(Geom%Elem_Blk(iBlk)%Elem_ID(iE))%ID_DoF
          End Do
       End If

       Call MPI_BCast(Buffer2, Geom%Elem_Blk(iBlk)%Num_Elems * 2 *            &
            & Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem, MPI_INTEGER, Source,    &
            & Comm, iErr)
       If (MyRank /= Source) Then
          Do iE = 1, Geom%Elem_Blk(iBlk)%Num_Elems
             Allocate(Elem_db(Geom%Elem_Blk(iBlk)%Elem_ID(iE))%ID_DoF         &
                  & (2 * Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem))
             Elem_db(Geom%Elem_Blk(iBlk)%Elem_ID(iE))%ID_DoF = Buffer2(iE, :)
          End Do
       End If
       DeAllocate(Buffer2)
    End Do

  End Subroutine Broadcast_Element2D_Elast

  Subroutine Broadcast_Element3D_Scal(Elem_db, Geom, Source, Communicator)
    Type(Element3D_Scal), Dimension(:), Pointer     :: Elem_db
    Type (EXO_Geom_Info), Intent(IN)                :: Geom
    Integer, Intent(IN)                             :: Source
    Integer, Intent(IN), Optional                   :: Communicator

    Integer                                         :: Comm, MyRank, iErr
    Integer                                         :: nE, eNum, iE, iBlk

    Integer                                         :: iSlice, SliceSize
    Integer                                         :: SliceIdxMin, SliceIdxMax
    Integer, Dimension(:), Pointer                  :: Buffer1
    Integer, Dimension(:,:), Pointer                :: Buffer2


    ! WARNING: We don't check if the elem_db are allocated or not
    ! We assume that Elem_db is allocated on source, but not on the 
    ! other cpu in te communicator.
    ! If this is not true, this subroutine is going to crash...

    ! We do NOT Broadcast the Gauss points informations. 
    ! It just does not make any sense to do so.

    ! Use the geom data and send element block after element block,
    ! since all element in a block must be of the same type

    ! We Assume that the Geom database has been broadcast, already

    If ( Present(Communicator) ) Then
       Comm = Communicator
    Else
       Comm = MPI_COMM_WORLD
    End If

    Call MPI_COMM_RANK(Comm, MyRank, iErr)
    If (MyRank /= Source) Then
       Allocate(Elem_db(Geom%Num_Elems))
    EndIf

    Do iSlice = 1, Geom%Num_Elems / MPI_MaxBuff + 1
       SliceIdxMin = (iSlice-1) * MPI_MaxBuff + 1
       SliceIdxMax = Min(Geom%Num_Elems, iSlice * MPI_MaxBuff)
       SliceSize   = SliceIdxMax - SliceIdxMin + 1
!!$       If (MyRank == Source) Then
!!$          Print*, 'Buffer1, Slice# ', iSlice, SliceIdxMin, SliceIdxMax,   &
!!$               & SliceSize
!!$       End If
       
       Allocate (Buffer1(SliceSize))

       If (MyRank == Source) Then
          Buffer1 = Elem_db(SliceIdxMin:SliceIdxMax)%NB_DoF
       End If
       Call MPI_BCAST(Buffer1, SliceSize, MPI_INTEGER, Source, Comm, iErr)
       If (MyRank /= Source) Then
          Elem_db(SliceIdxMin:SliceIdxMax)%NB_DoF = Buffer1
       End If
       
       If (MyRank == Source) Then
          Buffer1 = Elem_db(SliceIdxMin:SliceIdxMax)%NB_Gauss
       End If
       Call MPI_BCAST(Buffer1, SliceSize, MPI_INTEGER, Source, Comm, iErr)
       If (MyRank /= Source) Then
          Elem_db(SliceIdxMin:SliceIdxMax)%NB_Gauss = Buffer1
       End If
       
       If (MyRank == Source) Then
          Buffer1 = Elem_db(SliceIdxMin:SliceIdxMax)%ID_EL
       End If
       Call MPI_BCAST(Buffer1, SliceSize, MPI_INTEGER, Source, Comm, iErr)
       If (MyRank /= Source) Then
          Elem_db(SliceIdxMin:SliceIdxMax)%ID_EL = Buffer1
       End If
       
       DeAllocate(Buffer1)
    End Do

!!$    If (MyRank == Source) Then
!!$       Print*, 'OK Buffer 1'
!!$    End If

    Do iBlk = 1, Geom%Num_Elem_Blks
       Do iSlice = 1, Geom%Elem_Blk(iBlk)%Num_Elems  *                        &
            & Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem / MPI_MaxBuff + 1
          SliceIdxMin = (iSlice-1) * MPI_MaxBuff /                            &
               & Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem + 1
          SliceIdxMax = Min(Geom%Elem_Blk(iBlk)%Num_Elems,                    &
               & iSlice * MPI_MaxBuff / Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem)
          SliceSize   = SliceIdxMax - SliceIdxMin + 1
          
!!$          If (MyRank == Source) Then
!!$             Print*, 'Buffer2, Slice# ', iSlice, SliceIdxMin, SliceIdxMax,   &
!!$                  & SliceSize
!!$          End If
          Allocate(Buffer2(SliceIdxMin:SliceIdxMax,                           &
               & Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem))

          If (MyRank == Source) Then
             Do iE = SliceIdxMin, SliceIdxMax
                Buffer2(iE,:) =                                              &
                     & Elem_db(Geom%Elem_Blk(iBlk)%Elem_ID(iE))%ID_DoF
             End Do
          End If

          Call MPI_BCast(Buffer2, SliceSize *                                 &
               & Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem, MPI_INTEGER, Source, &
               & Comm, iErr)

          If (MyRank /= Source) Then
             Do iE = SliceIdxMin, SliceIdxMax
                Allocate(Elem_db(Geom%Elem_Blk(iBlk)%Elem_ID(iE))%ID_DoF      &
                     & (Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem))
                Elem_db(Geom%Elem_Blk(iBlk)%Elem_ID(iE))%ID_DoF = Buffer2(iE,:)
             End Do
          End If
          DeAllocate(Buffer2)
       End Do
    End Do
  End Subroutine Broadcast_Element3D_Scal


  Subroutine Broadcast_Element3D(Elem_db, Geom, Source, Communicator)
    Type(Element3D), Dimension(:), Pointer          :: Elem_db
    Type (EXO_Geom_Info), Intent(IN)                :: Geom
    Integer, Intent(IN)                             :: Source
    Integer, Intent(IN), Optional                   :: Communicator

    Integer                                         :: Comm, MyRank, iErr
    Integer                                         :: nE, eNum, iE, iBlk

    Integer, Dimension(:), Pointer                  :: Buffer1
    Integer, Dimension(:,:), Pointer                :: Buffer2


    ! WARNING: We don't check if the elem_db are allocated or not
    ! We assume that Elem_db is allocated on source, but not on the 
    ! other cpu in te communicator.
    ! If this is not true, this subroutine is going to crash...

    ! We do NOT Broadcast the Gauss points informations. 
    ! It just does not make any sense to do so.

    ! Use the geom data and send element block after element block,
    ! since all element in a block must be of the same type

    ! We Assume that the Geom database has been broadcast, already

    If ( Present(Communicator) ) Then
       Comm = Communicator
    Else
       Comm = MPI_COMM_WORLD
    End If

    Call MPI_COMM_RANK(Comm, MyRank, iErr)
    If (MyRank /= Source) Then
       Allocate(Elem_db(Geom%Num_Elems))
    EndIf

    Allocate (Buffer1(Geom%Num_Elems))

    If (MyRank == Source) Then
       Buffer1 = Elem_db(:)%NB_DoF
    End If
    Call MPI_BCAST(Buffer1, Geom%Num_Elems, MPI_INTEGER, Source, Comm, iErr)
    If (MyRank /= Source) Then
       Elem_db(:)%NB_DoF = Buffer1
    End If

    If (MyRank == Source) Then
       Buffer1 = Elem_db(:)%NB_Gauss
    End If
    Call MPI_BCAST(Buffer1, Geom%Num_Elems, MPI_INTEGER, Source, Comm, iErr)
    If (MyRank /= Source) Then
       Elem_db(:)%NB_Gauss = Buffer1
    End If

    If (MyRank == Source) Then
       Buffer1 = Elem_db(:)%ID_EL
    End If
    Call MPI_BCAST(Buffer1, Geom%Num_Elems, MPI_INTEGER, Source, Comm, iErr)
    If (MyRank /= Source) Then
       Elem_db(:)%ID_EL = Buffer1
    End If

    DeAllocate(Buffer1)

    Do iBlk = 1, Geom%Num_Elem_Blks
       Allocate(Buffer2(Geom%Elem_Blk(iBlk)%Num_Elems,                        &
            & 3 * Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem))

       If (MyRank == Source) Then
          Do iE = 1, Geom%Elem_Blk(iBlk)%Num_Elems
             Buffer2(iE, :) = Elem_db(Geom%Elem_Blk(iBlk)%Elem_ID(iE))%ID_DoF
          End Do
       End If

       Call MPI_BCast(Buffer2, Geom%Elem_Blk(iBlk)%Num_Elems * 3 *            &
            & Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem, MPI_INTEGER, Source,    &
            & Comm, iErr)
       If (MyRank /= Source) Then
          Do iE = 1, Geom%Elem_Blk(iBlk)%Num_Elems
             Allocate(Elem_db(Geom%Elem_Blk(iBlk)%Elem_ID(iE))%ID_DoF         &
                  & (3 * Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem))
             Elem_db(Geom%Elem_Blk(iBlk)%Elem_ID(iE))%ID_DoF = Buffer2(iE, :)
          End Do
       End If
       DeAllocate(Buffer2)
    End Do
  End Subroutine Broadcast_Element3D

  Subroutine Broadcast_Element3D_Elast(Elem_db, Geom, Source, Communicator)
    Type(Element3D_Elast), Dimension(:), Pointer    :: Elem_db
    Type (EXO_Geom_Info), Intent(IN)                :: Geom
    Integer, Intent(IN)                             :: Source
    Integer, Intent(IN), Optional                   :: Communicator

    Integer                                         :: Comm, MyRank, iErr
    Integer                                         :: nE, eNum, iE, iBlk

    Integer, Dimension(:), Pointer                  :: Buffer1
    Integer, Dimension(:,:), Pointer                :: Buffer2


    ! WARNING: We don't check if the elem_db are allocated or not
    ! We assume that Elem_db is allocated on source, but not on the 
    ! other cpu in te communicator.
    ! If this is not true, this subroutine is going to crash...

    ! We do NOT Broadcast the Gauss points informations. 
    ! It just does not make any sense to do so.

    ! Use the geom data and send element block after element block,
    ! since all element in a block must be of the same type

    ! We Assume that the Geom database has been broadcast, already

    If ( Present(Communicator) ) Then
       Comm = Communicator
    Else
       Comm = MPI_COMM_WORLD
    End If

    Call MPI_COMM_RANK(Comm, MyRank, iErr)
    If (MyRank /= Source) Then
       Allocate(Elem_db(Geom%Num_Elems))
    EndIf

    Allocate (Buffer1(Geom%Num_Elems))

    If (MyRank == Source) Then
       Buffer1 = Elem_db(:)%NB_DoF
    End If
    Call MPI_BCAST(Buffer1, Geom%Num_Elems, MPI_INTEGER, Source, Comm, iErr)
    If (MyRank /= Source) Then
       Elem_db(:)%NB_DoF = Buffer1
    End If

    If (MyRank == Source) Then
       Buffer1 = Elem_db(:)%NB_Gauss
    End If
    Call MPI_BCAST(Buffer1, Geom%Num_Elems, MPI_INTEGER, Source, Comm, iErr)
    If (MyRank /= Source) Then
       Elem_db(:)%NB_Gauss = Buffer1
    End If

    If (MyRank == Source) Then
       Buffer1 = Elem_db(:)%ID_EL
    End If
    Call MPI_BCAST(Buffer1, Geom%Num_Elems, MPI_INTEGER, Source, Comm, iErr)
    If (MyRank /= Source) Then
       Elem_db(:)%ID_EL = Buffer1
    End If

    DeAllocate(Buffer1)

    Do iBlk = 1, Geom%Num_Elem_Blks
       Allocate(Buffer2(Geom%Elem_Blk(iBlk)%Num_Elems,                        &
            & 3 * Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem))

       If (MyRank == Source) Then
          Do iE = 1, Geom%Elem_Blk(iBlk)%Num_Elems
             Buffer2(iE, :) = Elem_db(Geom%Elem_Blk(iBlk)%Elem_ID(iE))%ID_DoF
          End Do
       End If

       Call MPI_BCast(Buffer2, Geom%Elem_Blk(iBlk)%Num_Elems * 3 *            &
            & Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem, MPI_INTEGER, Source,    &
            & Comm, iErr)
       If (MyRank /= Source) Then
          Do iE = 1, Geom%Elem_Blk(iBlk)%Num_Elems
             Allocate(Elem_db(Geom%Elem_Blk(iBlk)%Elem_ID(iE))%ID_DoF         &
                  & (3 * Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem))
             Elem_db(Geom%Elem_Blk(iBlk)%Elem_ID(iE))%ID_DoF = Buffer2(iE, :)
          End Do
       End If
       DeAllocate(Buffer2)
    End Do

  End Subroutine Broadcast_Element3D_Elast


  Subroutine Show_One_EXO_Geom_Info(Geom, PrinterRank)
    Type (EXO_Geom_Info), Intent(IN)                :: Geom
    Integer, Intent(IN)                             :: PrinterRank
    Integer myrank, iErr
    Integer status(MPI_STATUS_SIZE), i 
    
    Call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, iErr)

    If (myrank == PrinterRank) then
       Write(*,100)
       Write(*,101) myrank, Trim(Geom%filename)
       Write(*,102) myrank, Trim(Geom%title)
       Write(*,103) myrank, Geom%num_dim
       Write(*,104) myrank, Geom%num_nodes
       Write(*,105) myrank, Geom%num_elems
       Write(*,106) myrank, Geom%Num_Elem_blks
       Write(*,107) myrank, Geom%Num_Node_Sets
       Write(*,108) myrank, Geom%Num_Side_Sets
       Write(*,*)
       Write(*,200) 
       Write(*,201) myrank, Geom%num_elem_blks
       Do i = 1, Geom%Num_Elem_blks
          Write(*,202) myrank, Geom%Elem_Blk(i)%ID, Geom%Elem_blk(i)%Type
          Write(*,203) myrank, Geom%Elem_Blk(i)%ID, Geom%Elem_blk(i)%Num_Elems
          Write(*,204) myrank, Geom%Elem_Blk(i)%ID,                           &
               & Geom%Elem_blk(i)%Num_Nodes_Per_Elem
          Write(*,205) myrank, Geom%Elem_Blk(i)%ID, Geom%Elem_blk(i)%Num_Attr
          Write(*,206, advance = 'no') myrank, Geom%Elem_Blk(i)%ID
          Write(*,*) Geom%Elem_blk(i)%Elem_ID
          
       End Do
       Write(*,300)
       Write(*,301) myrank, Geom%num_node_sets
       Do i = 1, Geom%num_node_sets
          Write(*,302) myrank, Geom%Node_Set(i)%ID, Geom%Node_Set(i)%Num_Nodes
          Write(*,303, advance = 'no') myrank, Geom%Node_Set(i)%ID
          Write(*,*) Geom%Node_Set(i)%Node_ID
          Write(*,304) myrank, Geom%Node_Set(i)%ID,                           &
               & Geom%Node_Set(i)%Num_Dist_Factors
       End Do
       Write(*,*)
       Write(*,400)
       Write(*,401) myrank, Geom%num_side_sets
       
       Write(*,*)
       Write(*,500)
       Write(*,501) myrank, Geom%num_QA
       Do i = 1, Geom%num_QA
          Write(*,502) myrank, i, Trim(Geom%QA_rec(1,i))
          Write(*,503) myrank, i, Trim(Geom%QA_rec(2,i))
          Write(*,504) myrank, i, Trim(Geom%QA_rec(3,i))
          Write(*,505) myrank, i, Trim(Geom%QA_rec(4,i))
       End Do
    End if

100 Format('rank: ', i2, ' IDs:      ', 3(i3,'  '), ' Block ', i3, ' Value',  &
         & 2(ES12.5, '  '))
101 Format('rank: ', I2,'    Filename ======================== ', A)
102 Format('rank: ', I2,'    Title =========================== ', A)
103 Format('rank: ', I2,'    Number of dimensions ============ ', I1)
104 Format('rank: ', I2,'    Number of nodes ================= ', I6)
105 Format('rank: ', I2,'    Number of elements ============== ', I6)
106 Format('rank: ', I2,'    Number of elements blocks ======= ', I6)
107 Format('rank: ', I2,'    Number of node sets ============= ', I6)
108 Format('rank: ', I2,'    Number of side sets ============= ', I6)
200 Format('rank: ', '*** ELEMENT BLOCKS ***')
201 Format('rank: ', I2,'    Number of blocks ================ ', I4)
202 Format('rank: ', I2,'    Block ', I3, ' Elements type ========= ', A)
203 Format('rank: ', I2,'    Block ', I3, ' Number of elements ==== ', I4)
204 Format('rank: ', I2,'    Block ', I3, ' Number of nodes per elt ', I4)
205 Format('rank: ', I2,'    Block ', I3, ' Number of attributes == ', I4)
206 Format('rank: ', I2,'    Block ', I3, ' IDs: ')
300 Format('*** NODE SETS ***')
301 Format('rank: ', I2,'    Number of sets ================== ', I4)
302 Format('rank: ', I2,'    Set ', I3, ' Number of nodes ========= ', I4)
303 Format('rank: ', I2,'    Set ', I3, ' IDs: ')
304 Format('rank: ', I2,'    Set ', I3, ' Number of dist. factors = ', I4)
    
400 Format('*** SIDE SETS ***')
401 Format('rank: ', I2,'    Number of side sets ============= ', I4)
    
500 Format('*** QA ***')
501 Format('rank: ', I2,'    Number of QA Records =========== ', I2)
502 Format('rank: ', I2,'    Rec ', I2, ' analysis code ============ ', A)
503 Format('rank: ', I2,'    Rec ', I2, ' analysis QA desc. ======== ', A)
504 Format('rank: ', I2,'    Rec ', I2, ' analysis date ============ ', A)
505 Format('rank: ', I2,'    Rec ', I2, ' analysis time ============ ', A)
  End Subroutine Show_One_EXO_Geom_Info
  
  Subroutine Show_All_EXO_Geom_Info(Out_Geom)
    Type (EXO_Geom_Info), Intent(IN)                 :: Out_Geom
    Integer myrank, iErr
    Integer status(MPI_STATUS_SIZE), nprocs, i, k
    
    Call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, iErr)
    Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, iErr)

    Do k = 0, nprocs-1
       If (k == myrank) then
          Write(*,101) myrank, Trim(Out_Geom%filename)
          Write(*,102) myrank, Trim(Out_Geom%title)	
          Write(*,103) myrank, Out_Geom%num_dim
          Write(*,104) myrank, Out_Geom%num_nodes
          Write(*,105) myrank, Out_Geom%num_elems
          Write(*,106) myrank, Out_Geom%Num_Elem_blks
          Write(*,107) myrank, Out_Geom%Num_Node_Sets
          Write(*,108) myrank, Out_Geom%Num_Side_Sets
          Do i = 1, Out_Geom%Num_Elem_blks
             Write(*,202) myrank, Out_Geom%Elem_Blk(i)%ID,                    &
                  & Out_Geom%Elem_blk(i)%Type
             Write(*,203) myrank, Out_Geom%Elem_Blk(i)%ID,                    &
                  & Out_Geom%Elem_blk(i)%Num_Elems
             Write(*,204) myrank, Out_Geom%Elem_Blk(i)%ID,                    &
                  & Out_Geom%Elem_blk(i)%Num_Nodes_Per_Elem
             Write(*,205) myrank, Out_Geom%Elem_Blk(i)%ID,                    &
                  & Out_Geom%Elem_blk(i)%Num_Attr
             Write(*,206, advance = 'no') myrank, Out_Geom%Elem_Blk(i)%ID
             Write(*,*) Out_Geom%Elem_blk(i)%Elem_ID
          End Do
          Write(*,301) myrank, Out_Geom%num_node_sets
          Do i = 1, Out_Geom%num_node_sets
             Write(*,302) myrank, Out_Geom%Node_Set(i)%ID,                    &
                  & Out_Geom%Node_Set(i)%Num_Nodes
             Write(*,303, advance = 'no') myrank, Out_Geom%Node_Set(i)%ID
             Write(*,*) Out_Geom%Node_Set(i)%Node_ID
             Write(*,304) myrank, Out_Geom%Node_Set(i)%ID,                   &
                  & Out_Geom%Node_Set(i)%Num_Dist_Factors
          End Do
          Write(*,401) myrank, Out_Geom%num_side_sets
          Write(*,501) myrank, Out_Geom%num_QA
          Do i = 1, Out_Geom%num_QA
             Write(*,502) myrank, i, Trim(Out_Geom%QA_rec(1,i))
             Write(*,503) myrank, i, Trim(Out_Geom%QA_rec(2,i))
             Write(*,504) myrank, i, Trim(Out_Geom%QA_rec(3,i))
             Write(*,505) myrank, i, Trim(Out_Geom%QA_rec(4,i))
          End Do
       End If
    End Do
    Call MPI_BARRIER(MPI_COMM_WORLD, iErr)
    
100 Format('rank: ', I2, ' IDs:      ', 3(i3,'  '), ' Block ', i3, ' Value',  &
         & 2(ES12.5, '  '))
101 Format('rank: ', I2,'    Filename ======================== ', A)
102 Format('rank: ', I2,'    Title =========================== ', A)
103 Format('rank: ', I2,'    Number of dimensions ============ ', I1)
104 Format('rank: ', I2,'    Number of nodes ================= ', I6)
105 Format('rank: ', I2,'    Number of elements ============== ', I6)
106 Format('rank: ', I2,'    Number of elements blocks ======= ', I6)
107 Format('rank: ', I2,'    Number of node sets ============= ', I6)
108 Format('rank: ', I2,'    Number of side sets ============= ', I6)
200 Format('rank: ', I2,'*** ELEMENT BLOCKS ***')
201 Format('rank: ', I2,'    Number of blocks ================ ', I4)
202 Format('rank: ', I2,'    Block ', I3, ' Elements type ========= ', A)
203 Format('rank: ', I2,'    Block ', I3, ' Number of elements ==== ', I4)
204 Format('rank: ', I2,'    Block ', I3, ' Number of nodes per elt ', I4)
205 Format('rank: ', I2,'    Block ', I3, ' Number of attributes == ', I4)
206 Format('rank: ', I2,'    Block ', I3, ' IDs: ')
300 Format('rank: ', I2,'*** NODE SETS ***')
301 Format('rank: ', I2,'    Number of node sets ================== ', I4)
302 Format('rank: ', I2,'    Set ', I3, ' Number of nodes ========= ', I4)
303 Format('rank: ', I2,'    Set ', I3, ' IDs: ')
304 Format('rank: ', I2,'    Set ', I3, ' Number of dist. factors = ', I4)

400 Format('rank: ', I2,'*** SIDE SETS ***')
401 Format('rank: ', I2,'    Number of side sets ================== ', I4)
    
500 Format('rank: ', I2,'*** QA ***')
501 Format('rank: ', I2,'    Number of QA Records ================= ', I4)
502 Format('rank: ', I2,'    Rec ', I2, ' analysis code ============ ', A)
503 Format('rank: ', I2,'    Rec ', I2, ' analysis QA desc. ======== ', A)
504 Format('rank: ', I2,'    Rec ', I2, ' analysis date ============ ', A)
505 Format('rank: ', I2,'    Rec ', I2, ' analysis time ============ ', A)

  End Subroutine Show_All_EXO_Geom_Info

  
End Module m_MEF_Communication
