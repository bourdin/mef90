Module m_MEF_EXO
#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmeshdef.h"

   Use m_MEF_LinAlg
   Use m_MEF_Parameters
   Use m_MEF_Types
   Use m_MEF_MPI
   Use m_MEF_Elements
   Use petsc
   Use petscvec
   Use petscmesh
   
   
   IMPLICIT NONE
   Private

   Integer, Parameter, Public                        :: exo_cpu_ws = 8
   Integer, Parameter, Public                        :: exo_io_ws = 8
   
   Public :: Write_MeshTopology
   Public :: Write_MeshTopologyGlobal

   Public :: Read_EXO_Result_Global
   Public :: Write_EXO_Result_Global   
   Public :: Read_EXO_Result_Vertex
   Public :: Write_EXO_Result_Vertex
   
   Interface Read_EXO_Result_Vertex
      Module Procedure Read_EXO_Result_VertexPtrInterlaced, Read_EXO_Result_VertexSection, Read_EXO_Result_VertexVec, Read_EXO_Result_VertexVect2D, Read_EXO_Result_VertexVect3D
   End Interface Read_EXO_Result_Vertex

   Interface Write_EXO_Result_Vertex
      Module Procedure Write_EXO_Result_VertexPtrInterlaced, Write_EXO_Result_VertexSection, Write_EXO_Result_VertexVec, Write_EXO_Result_VertexVect2D, Write_EXO_Result_VertexVect3D
   End Interface Write_EXO_Result_Vertex

!   Public :: Read_EXO_Result_Ptr_Nodes   
!   Public :: Write_EXO_Result_Ptr_Nodes   
!   Public :: Read_EXO_Result_Elems
!   Public :: Write_EXO_Result_Elems   
!   Public :: Read_EXO_Result_Ptr_Elems
!   Public :: Write_EXO_Result_Ptr_Elems   

!   Interface Read_EXO_Result_Elems
!      Module Procedure Read_EXO_Result_Scal_Elems, Read_EXO_Result_Vect2D_Elems, Read_EXO_Result_Mat2D_Elems, Read_EXO_Result_MatS2D_Elems, Read_EXO_Result_Vect3D_Elems, Read_EXO_Result_Mat3D_Elems, Read_EXO_Result_MatS3D_Elems
!   End Interface
!   
!   Interface Write_EXO_Result_Elems
!      Module Procedure Write_EXO_Result_Scal_Elems, Write_EXO_Result_Vect2D_Elems, Write_EXO_Result_Mat2D_Elems, Write_EXO_Result_MatS2D_Elems, Write_EXO_Result_Vect3D_Elems, Write_EXO_Result_Mat3D_Elems, Write_EXO_Result_MatS3D_Elems
!   End Interface

 Contains
 
   Subroutine Uniq(dComm, dMyVals, dVals)
      MPI_Comm                         :: dComm
      PetscInt, Dimension(:), Pointer  :: dMyVals, dVals
      
      Logical, Dimension(:), Pointer   :: ValCount
      PetscInt                         :: GlobMinVal, MyMinVal
      PetscInt                         :: GlobMaxVal, MyMaxVal
      PetscInt                         :: UniqCount
      PetscMPIInt                      :: rank
      PetscInt                         :: i, j, iErr

      Call MPI_Comm_Rank(PETSC_COMM_WORLD, rank, iErr)

      MyMinVal = MinVal(dMyVals)
      MyMaxVal = MaxVal(dMyVals)
      Call MPI_AllReduce(MyMinVal, GlobMinVal, 1, MPI_INTEGER, MPI_MIN, dComm, iErr)
      Call MPI_AllReduce(MyMaxVal, GlobMaxVal, 1, MPI_INTEGER, MPI_MAX, dComm, iErr)


      Allocate(ValCount(GlobMinVal:GlobMaxVal))
      ValCount = .FALSE.
      Do i = 1, Size(dMyVals)
         ValCount(dMyVals(i)) = .TRUE.
      End Do

      Call MPI_AllReduce(MPI_IN_PLACE, ValCount, GlobMaxVal-GlobMinVal+1, MPI_INTEGER, MPI_LOR, dComm, iErr)
      !!! This is suboptimal. I could gather only to CPU 0 and do everything else on CPU 0 before broadcasting
      
      UniqCount = Count(ValCount)

      Allocate(dVals(UniqCount))
      j = 1
      Do i = GlobMinVal, GlobMaxVal
         If (ValCount(i)) Then
            dVals(j) = i
            j = j+1
         End If
      End Do
      DeAllocate(ValCount)
   End Subroutine Uniq

   Subroutine Write_MeshTopology(dMeshTopology, dEXO)
      Type(MeshTopology_Info)                        :: dMeshTopology
      Type(EXO_Info)                                 :: dEXO
      PetscInt                                       :: vers
      PetscInt                                       :: iErr
      PetscInt                                       :: iDummy
      PetscReal                                      :: rDummy
      Character                                      :: cDummy
      
      PetscInt                                       :: iBlk, iSet, i
      PetscInt                                       :: iE, iELoc
      PetscInt                                       :: Num_Attr = 0         
      PetscInt                                       :: Num_Dist_Factor = 0
      
      Character(len=MXSTLN), Dimension(3)            :: Coord_Names, Elem_Type
      PetscReal, Dimension(:,:), Pointer             :: Coordinates
      PetscInt, Dimension(:,:), Pointer              :: ConnectMesh
      PetscInt, Dimension(:), Pointer                :: ConnectBlk
      PetscInt                                       :: offset
      
      Coord_Names(1) = 'X'
      Coord_Names(2) = 'Y'
      Coord_Names(3) = 'Z'
      
      Is_IO: If (dEXO%comm == PETSC_COMM_SELF) Then
         ! Open File
         dEXO%exoid = EXCRE (dEXO%filename, EXCLOB, exo_cpu_ws, exo_io_ws, iErr)
         dEXO%exoid = EXOPEN(dEXO%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers, ierr)
         
         ! Write Mesh Topology Info
         Call EXPINI(dEXO%exoid, dEXO%Title, dMeshTopology%Num_Dim, dMeshTopology%Num_Verts, dMeshTopology%Num_Elems, dMeshTopology%Num_Elem_Blks, dMeshTopology%Num_Node_Sets, dMeshTopology%Num_Side_Sets, iErr)
            
         ! Write Coordinate Names
         Call EXPCON(dEXO%exoid, Coord_Names, iErr)
         
         ! Write Elem blocks informations
         Do iBlk = 1, dMeshTopology%Num_Elem_Blks
            Select Case (dMeshTopology%elem_blk(iBlk)%Elem_Type)
            Case (MEF90_P1_Lagrange)
               Select Case (dMeshTopology%Num_Dim)
               Case (2)
                  Elem_Type = 'TRI3'
               Case (3)
                  Elem_Type = 'TETRA4'
               Case Default
                  SETERRQ(PETSC_ERR_SUP, 'Only 2 and 3 dimensional elements are supported', iErr)
               End Select
            Case (MEF90_P2_Lagrange)
               Select Case (dMeshTopology%Num_Dim)
               Case (2)
                  Elem_Type = 'TRI6'
               Case (3)
                  Elem_Type = 'TETRA10'
               Case Default
                  SETERRQ(PETSC_ERR_SUP, 'Only 2 and 3 dimensional elements are supported', iErr)
               End Select
            Case Default
                  SETERRQ(PETSC_ERR_SUP, 'Only MEF90_P1_Lagrange and MEF90_P2_Lagrange elements are supported', iErr)
            End Select
            Call EXPELB(dEXO%exoid, dMeshTopology%elem_blk(iBlk)%ID, Elem_Type, dMeshTopology%elem_blk(iBlk)%Num_Elems, dMeshTopology%elem_blk(iBlk)%Num_DoF, Num_Attr, iErr)
         End Do
   
         ! Write Side sets informations
         If (dMeshTopology%Num_Side_Sets > 0) Then
            SETERRQ(PETSC_ERR_SUP, 'Side sets not supported yet', iErr)
         End If
         
         ! Write Node sets informations
         Do iSet = 1, dMeshTopology%Num_Node_Sets
            Call EXPNP(dEXO%exoid, dMeshTopology%Node_Set(iSet)%ID, dMeshTopology%Node_Set(iSet)%Num_Nodes, Num_Dist_Factor, iErr)
            If (dMeshTopology%Node_Set(iSet)%Num_Nodes>0) Then
               Call EXPNS(dEXO%exoid, dMeshTopology%Node_Set(iSet)%ID, dMeshTopology%Node_Set(iSet)%Node_ID(:), iErr)
            End If
         End Do
   
         ! Write vertex coordinates
         Call MeshGetCoordinatesF90(dMeshTopology%mesh, Coordinates, iErr)
         Call EXPCOR(dEXO%exoid, Coordinates(:,1), Coordinates(:,2), Coordinates(:,3), iErr)
         Call MeshRestoreCoordinatesF90(dMeshTopology%mesh, Coordinates, iErr)
         
          ! Write Connectivity tables
         Call MeshGetElementsF90(dMeshTopology%mesh, ConnectMesh, iErr)
         Do iBlk = 1, dMeshTopology%Num_Elem_Blks
            If (dMeshTopology%Elem_Blk(iBlk)%Num_Elems > 0) Then
               Allocate (ConnectBlk(dMeshTopology%Elem_Blk(iBlk)%Num_Elems * dMeshTopology%Elem_Blk(iBlk)%Num_DoF))
               
               Do iELoc = 1, dMeshTopology%Elem_Blk(iBlk)%Num_Elems
                  iE = dMeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
                   offset = (iEloc-1) * dMeshTopology%Elem_Blk(iBlk)%Num_DoF
                   ConnectBlk(offset+1: offset+dMeshTopology%Elem_Blk(iBlk)%Num_DoF) = ConnectMesh(iE, :)
                End Do
                Call EXPELC (dEXO%exoid, iBlk, ConnectBlk, iErr)
                DeAllocate(ConnectBlk)
             End If
          End Do
          Call MeshRestoreElementsF90(dMeshTopology%mesh, ConnectMesh, iErr)

         Call EXCLOS(dEXO%exoid, iErr)
         dEXO%exoid = 0
      Else
         SETERRQ(PETSC_ERR_SUP, 'Synchronized I/O on PETSC_COMM_WORLD not supported yet', ierr)
      End If Is_IO
   End Subroutine Write_MeshTopology
   
   
   Subroutine Write_MeshTopologyGlobal(dMeshTopology, dEXO, dGlobalComm)
      !!! Reconstruct the topology informatins for of all meshes in dGlobalComm then write it
      Type(MeshTopology_Info)                        :: dMeshTopology
      Type(EXO_Info)                                 :: dEXO
      MPI_Comm                                       :: dGlobalComm
      
      PetscInt                                       :: iBlk, iSet, i, iErr
      PetscInt                                       :: iE, iELoc
      PetscInt                                       :: Num_Attr = 0         
      PetscInt                                       :: Num_Dist_Factor = 0
      
      Type(MeshTopology_Info)                        :: GlobalMeshTopology
      PetscInt, Dimension(:), Pointer                :: Tmp_GlobalID, Tmp_ID


      ! Gather Global Sizes
      GlobalMeshTopology%num_dim  = dMeshTopology%num_dim
      GlobalMeshTopology%num_verts = dMeshTopology%num_verts
      GlobalMeshTopology%num_elems = dMeshTopology%num_elems
      
      Allocate(Tmp_ID(dMeshTopology%num_elem_blks))
      Tmp_ID = dMeshTopology%elem_blk(:)%ID
      Call Uniq(dGlobalComm, Tmp_ID, Tmp_GlobalID)            
      GlobalMeshTopology%Num_elem_blks = Size(Tmp_GlobalID)
      Allocate(GlobalMeshTopology%elem_blk(GlobalMeshTopology%num_elem_blks))
      GlobalMeshTopology%elem_blk(:)%ID = Tmp_GlobalID
      DeAllocate(Tmp_ID)

      GlobalMeshTopology%num_side_sets = 0

      Allocate(Tmp_ID(dMeshTopology%num_node_sets))
      Tmp_ID = dMeshTopology%node_set(:)%ID
      Call Uniq(dGlobalComm, Tmp_ID, Tmp_GlobalID)            
      GlobalMeshTopology%Num_node_sets = Size(Tmp_GlobalID)
      Allocate(GlobalMeshTopology%node_set(GlobalMeshTopology%num_node_sets))
      GlobalMeshTopology%node_set(:)%ID = Tmp_GlobalID
      DeAllocate(Tmp_ID)
      
      ! Element Blocks
      Allocate(Tmp_ID(GlobalMeshTopology%num_elem_blks))
      Tmp_ID = 0
      Do i = 1, GlobalMeshTopology%num_elem_blks
         GlobalMeshTopology%elem_blk(i)%num_Elems = 0
         GlobalMeshTopology%elem_blk(i)%DoF_Location(:) = 0
         GlobalMeshTopology%elem_blk(i)%num_DoF = 0
         Do iBlk = 1, dMeshTopology%num_elem_blks
            If (GlobalMeshTopology%elem_blk(i)%ID == dMeshTopology%elem_blk(iBlk)%ID) Then
               Tmp_ID(i) = dMeshTopology%elem_blk(iBlk)%Elem_Type
               GlobalMeshTopology%elem_blk(i)%num_elems = dMeshTopology%elem_blk(iBlk)%num_elems
               Allocate(GlobalMeshTopology%elem_blk(i)%elem_id(GlobalMeshTopology%elem_blk(i)%num_elems))
               GlobalMeshTopology%elem_blk(i)%elem_id = dMeshTopology%elem_blk(iBlk)%elem_id
               GlobalMeshTopology%elem_blk(i)%DoF_Location = dMeshTopology%elem_blk(iBlk)%DoF_Location
               GlobalMeshTopology%elem_blk(i)%num_DoF = dMeshTopology%elem_blk(iBlk)%num_DoF
            End If
         End Do
      End Do
      Call MPI_AllReduce(MPI_IN_PLACE, Tmp_ID, GlobalMeshTopology%num_elem_blks, MPI_INTEGER, MPI_MAX, dGlobalComm, iErr)
      GlobalMeshTopology%elem_blk(:)%Elem_Type = Tmp_ID
      DeAllocate(Tmp_ID)
      Do iBlk = 1, GlobalMeshTopology%num_elem_blks
         Call MPI_AllReduce(MPI_IN_PLACE, GlobalMeshTopology%elem_blk(iBlk)%DoF_Location, 4, MPI_INTEGER, MPI_MAX, dGlobalComm, iErr)
         Call MPI_AllReduce(MPI_IN_PLACE, GlobalMeshTopology%elem_blk(iBlk)%Num_DoF, 1, MPI_INTEGER, MPI_MAX, dGlobalComm, iErr)
      End Do
      
      ! Node Sets
      GlobalMeshTopology%node_set(:)%num_nodes = 0
      Do i = 1, GlobalMeshTopology%num_node_sets
         Do iSet = 1, dMeshTopology%num_node_sets
            If (GlobalMeshTopology%node_set(i)%ID == dMeshTopology%node_set(iSet)%ID) Then
               GlobalMeshTopology%node_set(i)%num_nodes = dMeshTopology%node_set(iSet)%num_nodes
               Allocate(GlobalMeshTopology%node_set(i)%node_ID(GlobalMeshTopology%node_set(i)%num_nodes))
               GlobalMeshTopology%node_set(i)%node_ID = dMeshTopology%node_set(iSet)%node_ID
            End If
         End Do
      End Do
      Do iSet = 1, GlobalMeshTopology%num_node_sets
      End Do

      GlobalMeshTopology%mesh = dMeshTopology%mesh
      Call Write_MeshTopology(GlobalMeshTopology, dEXO)
      Call MeshTopologyDestroy(GlobalMeshTopology)
   End Subroutine Write_MeshTopologyGlobal


!!!
!!! RESULT FILES
!!!

!!! READ
!!! In the sequel, we always assume that we are doing distributed I/O when EXO%comm == PETSC_COMM_SELF (i.e. 1 file per CPU) and sequential (i.e. IO operation on a single file on CPU 0) when EXO%comm == PETSC_COM_WORLD
   
   Subroutine Read_EXO_Result_Global(dEXO, dIdx, dTS, dRes, dNumRec)
      Type (EXO_Info), Intent(INOUT)                 :: dEXO
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      PetscReal                                      :: dRes
      PetscInt                                       :: dNumRec
      
      PetscReal                                      :: MyRes
      PetscInt                                       :: iErr
      PetscReal                                      :: Vers
      PetscReal, Dimension(:), Pointer               :: Tmp_Res
      PetscInt                                       :: Num_Vars, Num_TS
      PetscReal                                      :: fDum
      Character                                      :: cDum
      
      dRes = 0.0_Kr
      If ( ((dEXO%comm == PETSC_COMM_WORLD) .AND. (MEF90_MyRank == 0)) .OR. (dEXO%comm == PETSC_COMM_SELF) ) Then
         dEXO%exoid = EXOPEN(dEXO%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers, ierr)
         ! Get the number of global variables stored in the database    
         Call EXGVP(dEXO%exoid, 'G', Num_Vars, iErr)
         Allocate(Tmp_Res(Num_Vars))
         
         ! Read All the global variables at time step TS and copy the right one
         ! into Res
         Call EXGGV(dEXO%exoid, dTS, Num_Vars, Tmp_Res, iErr)
         MyRes = Tmp_Res(dIdx)
         DeAllocate (Tmp_Res)
         Call EXCLOS(dEXO%exoid, iErr)
         dEXO%exoid = 0
      End If
            
      !!! Broacast if dEXO%comm == PETSC_COMM_WORLD
      If (dEXO%comm == PETSC_COMM_WORLD) Then
         Call MPI_BCast(MyRes, 1, MPIU_SCALAR, 0, dEXO%comm, iErr)
         dRes = MyRes
      End If
   End Subroutine Read_EXO_Result_Global

   Subroutine Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Info), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Info)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      PetscReal, Dimension(:), Pointer               :: dRes
      
      PetscInt                                       :: Num_Rec, iRec
      PetscInt                                       :: iErr
      PetscReal                                      :: Vers
      
      dEXO%exoid = EXOPEN(dEXO%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers, ierr)
      
      If (Mod(Size(dRes), dMeshTopology%Num_Verts) /= 0) Then
         SETERRQ(PETSC_ERR_ARG_SIZ, 'Read_EXO_Result_VertexPtrInterlaced: The argument does not match the number of vertices in the mesh', iErr)
      End If
      Num_Rec = Size(dRes) / dMeshTopology%Num_Verts
      
      Do iRec = 1, Num_Rec
         Call EXGNV(dEXO%exoid, dTS, dIdx + iRec-1, dMeshTopology%Num_Verts, dRes(iRec:dMeshTopology%Num_Verts*Num_Rec:Num_Rec), iErr); CHKERRQ(iErr)
      End Do
      Call EXCLOS(dEXO%exoid, iErr)
      dEXO%exoid = 0
   End Subroutine Read_EXO_Result_VertexPtrInterlaced

   Subroutine Read_EXO_Result_VertexSection(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Info), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Info)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (SectionReal)                             :: dRes
      
      Type (Vec)                                     :: Res_Vec
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr
      
      !!! We Assume that the section is initialized and has the proper size
      Call SectionRealCreateLocalVector(dRes, Res_Vec, iErr); CHKERRQ(iErr)
      Call VecGetArrayF90(Res_Vec, Res_Ptr, iErr); CHKERRQ(iErr)
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx, dTS, Res_Ptr)
      Call VecRestoreArrayF90(Res_Vec, Res_Ptr, iErr); CHKERRQ(iErr)
      Call VecDestroy(Res_Vec, iErr); CHKERRQ(iErr)
   End Subroutine Read_EXO_Result_VertexSection

   Subroutine Read_EXO_Result_VertexVec(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Info), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Info)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Vec)                                     :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr
      
      !!! We Assume that the section is initialized and has the proper size
      Call VecGetArrayF90(dRes, Res_Ptr, iErr); CHKERRQ(iErr)
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx, dTS, Res_Ptr)
      Call VecRestoreArrayF90(dRes, Res_Ptr, iErr); CHKERRQ(iErr)
   End Subroutine Read_EXO_Result_VertexVec

   Subroutine Read_EXO_Result_VertexVect2D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Info), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Info)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Vect2D), Dimension(:), Pointer           :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr
      
      If ( Size(dRes) /= dMeshTopology%Num_Verts) Then
         SETERRQ(PETSC_ERR_ARG_SIZ, 'Read_EXO_Result_VertexVect2D: The argument does not match the number of vertices in the mesh', iErr)
      End If

      Allocate(Res_Ptr(dMeshTopology%Num_Verts))
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      dRes(:)%X = Res_Ptr
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+1, dTS, Res_Ptr)
      dRes(:)%Y = Res_Ptr
      DeAllocate(Res_Ptr)
   End Subroutine Read_EXO_Result_VertexVect2D

   Subroutine Read_EXO_Result_VertexVect3D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Info), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Info)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Vect3D), Dimension(:), Pointer           :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr
      
      If ( Size(dRes) /= dMeshTopology%Num_Verts) Then
         SETERRQ(PETSC_ERR_ARG_SIZ, 'Read_EXO_Result_VertexVect3D: The argument does not match the number of vertices in the mesh', iErr)
      End If

      Allocate(Res_Ptr(dMeshTopology%Num_Verts))
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      dRes(:)%X = Res_Ptr
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+1, dTS, Res_Ptr)
      dRes(:)%Y = Res_Ptr
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+2, dTS, Res_Ptr)
      dRes(:)%Z = Res_Ptr
      DeAllocate(Res_Ptr)
   End Subroutine Read_EXO_Result_VertexVect3D

!!!##  Subroutine Read_EXO_Result_Mat2D_Nodes(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    PetscInt                                       :: Idx
!!!##    PetscInt                                       :: TS
!!!##    Type(Mat2D), Dimension(:), Pointer             :: Res
!!!##
!!!##    PetscInt                                       :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##
!!!##!!$    If (Geom%Num_Dim /= 2) Then 
!!!##!!$       Write(*,*) Trim(Geom%Filename), 'is not a 2D File.'
!!!##!!$       Write(*,*) 'Aborting'
!!!##!!$       STOP
!!!##!!$    End If
!!!##    If (Associated(Res)) Then
!!!##       DeAllocate(Res)
!!!##    End If
!!!##    Allocate(Res(Geom%Num_Nodes))
!!!##    Call EXGNV(Geom%exoid, TS, Idx,   Geom%Num_Nodes, Res%XX, iErr)
!!!##    Call EXGNV(Geom%exoid, TS, Idx+1, Geom%Num_Nodes, Res%XY, iErr)
!!!##    Call EXGNV(Geom%exoid, TS, Idx+2, Geom%Num_Nodes, Res%YX, iErr)
!!!##    Call EXGNV(Geom%exoid, TS, Idx+3, Geom%Num_Nodes, Res%YY, iErr)
!!!##    
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Read_EXO_Result_Mat2D_Nodes
!!!##
!!!##  Subroutine Read_EXO_Result_MatS2D_Nodes(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    PetscInt                                       :: Idx
!!!##    PetscInt                                       :: TS
!!!##    Type(MatS2D), Dimension(:), Pointer            :: Res
!!!##
!!!##    PetscInt                                       :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##
!!!##!!$    If (Geom%Num_Dim /= 2) Then 
!!!##!!$       Write(*,*) Trim(Geom%Filename), 'is not a 2D File.'
!!!##!!$       Write(*,*) 'Aborting'
!!!##!!$       STOP
!!!##!!$    End If
!!!##    If (Associated(Res)) Then
!!!##       DeAllocate(Res)
!!!##    End If
!!!##    Allocate(Res(Geom%Num_Nodes))
!!!##    Call EXGNV(Geom%exoid, TS, Idx,   Geom%Num_Nodes, Res%XX, iErr)
!!!##    Call EXGNV(Geom%exoid, TS, Idx+1, Geom%Num_Nodes, Res%YY, iErr)
!!!##    Call EXGNV(Geom%exoid, TS, Idx+2, Geom%Num_Nodes, Res%XY, iErr)
!!!##    
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Read_EXO_Result_MatS2D_Nodes
!!!##
!!!##  Subroutine Read_EXO_Result_Vect3D_Nodes(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    PetscInt                                       :: Idx
!!!##    PetscInt                                       :: TS
!!!##    Type(Vect3D), Dimension(:), Pointer            :: Res
!!!##
!!!##    PetscInt                                       :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##
!!!##!!$    If (Geom%Num_Dim /= 3) Then 
!!!##!!$       Write(*,*) Trim(Geom%Filename), 'is not a 3D File.'
!!!##!!$       Write(*,*) 'Aborting'
!!!##!!$       STOP
!!!##!!$    End If
!!!##    If (Associated(Res)) Then
!!!##       DeAllocate(Res)
!!!##    End If
!!!##    Allocate(Res(Geom%Num_Nodes))
!!!##    Call EXGNV(Geom%exoid, TS, Idx,   Geom%Num_Nodes, Res%X, iErr)
!!!##    Call EXGNV(Geom%exoid, TS, Idx+1, Geom%Num_Nodes, Res%Y, iErr)
!!!##    Call EXGNV(Geom%exoid, TS, Idx+2, Geom%Num_Nodes, Res%Z, iErr)
!!!##    
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Read_EXO_Result_Vect3D_Nodes
!!!##
!!!##  Subroutine Read_EXO_Result_Mat3D_Nodes(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    PetscInt                                       :: Idx
!!!##    PetscInt                                       :: TS
!!!##    Type(Mat3D), Dimension(:), Pointer             :: Res
!!!##
!!!##    PetscInt                                       :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##
!!!##!!$    If (Geom%Num_Dim /= 3) Then 
!!!##!!$       Write(*,*) Trim(Geom%Filename), 'is not a 3D File.'
!!!##!!$       Write(*,*) 'Aborting'
!!!##!!$       STOP
!!!##!!$    End If
!!!##    If (Associated(Res)) Then
!!!##       DeAllocate(Res)
!!!##    End If
!!!##    Allocate(Res(Geom%Num_Nodes))
!!!##    Call EXGNV(Geom%exoid, TS, Idx,   Geom%Num_Nodes, Res%XX, iErr)
!!!##    Call EXGNV(Geom%exoid, TS, Idx+1, Geom%Num_Nodes, Res%XY, iErr)
!!!##    Call EXGNV(Geom%exoid, TS, Idx+2, Geom%Num_Nodes, Res%XZ, iErr)
!!!##    Call EXGNV(Geom%exoid, TS, Idx+3, Geom%Num_Nodes, Res%YX, iErr)
!!!##    Call EXGNV(Geom%exoid, TS, Idx+4, Geom%Num_Nodes, Res%YY, iErr)
!!!##    Call EXGNV(Geom%exoid, TS, Idx+5, Geom%Num_Nodes, Res%YZ, iErr)
!!!##    Call EXGNV(Geom%exoid, TS, Idx+6, Geom%Num_Nodes, Res%ZX, iErr)
!!!##    Call EXGNV(Geom%exoid, TS, Idx+7, Geom%Num_Nodes, Res%ZY, iErr)
!!!##    Call EXGNV(Geom%exoid, TS, Idx+8, Geom%Num_Nodes, Res%ZZ, iErr)
!!!##    
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Read_EXO_Result_Mat3D_Nodes
!!!##
!!!##  Subroutine Read_EXO_Result_MatS3D_Nodes(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    PetscInt                                       :: Idx
!!!##    PetscInt                                       :: TS
!!!##    Type(MatS3D), Dimension(:), Pointer            :: Res
!!!##
!!!##    PetscInt                                       :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##
!!!##!!$    If (Geom%Num_Dim /= 3) Then 
!!!##!!$       Write(*,*) Trim(Geom%Filename), 'is not a 3D File.'
!!!##!!$       Write(*,*) 'Aborting'
!!!##!!$       STOP
!!!##!!$    End If
!!!##    If (Associated(Res)) Then
!!!##       DeAllocate(Res)
!!!##    End If
!!!##    Allocate(Res(Geom%Num_Nodes))
!!!##    Call EXGNV(Geom%exoid, TS, Idx,   Geom%Num_Nodes, Res%XX, iErr)
!!!##    Call EXGNV(Geom%exoid, TS, Idx+1, Geom%Num_Nodes, Res%YY, iErr)
!!!##    Call EXGNV(Geom%exoid, TS, Idx+2, Geom%Num_Nodes, Res%ZZ, iErr)
!!!##    Call EXGNV(Geom%exoid, TS, Idx+3, Geom%Num_Nodes, Res%XY, iErr)
!!!##    Call EXGNV(Geom%exoid, TS, Idx+4, Geom%Num_Nodes, Res%YZ, iErr)
!!!##    Call EXGNV(Geom%exoid, TS, Idx+5, Geom%Num_Nodes, Res%XZ, iErr)
!!!##    
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Read_EXO_Result_MatS3D_Nodes
!!!##
!!!##  Subroutine Read_EXO_Result_Scal_Elems(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    PetscInt                                       :: Idx
!!!##    PetscInt                                       :: TS
!!!##    Real(Kind = Kr), Dimension(:), Pointer         :: Res
!!!##
!!!##    PetscInt                                       :: iBlk
!!!##    PetscInt                                       :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##
!!!##    Do iBlk = 1, Geom%Num_Elem_Blks
!!!##       Call EXGEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID), iErr)
!!!##    End Do
!!!##
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Read_EXO_Result_Scal_Elems
!!!##
!!!##  Subroutine Read_EXO_Result_Vect2D_Elems(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    PetscInt                                       :: Idx
!!!##    PetscInt                                       :: TS
!!!##    Type(Vect2D), Dimension(:), Pointer            :: Res
!!!##
!!!##    PetscInt                                       :: iBlk
!!!##    PetscInt                                       :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##
!!!##    Do iBlk = 1, Geom%Num_Elem_Blks
!!!##       Call EXGEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%X, iErr)
!!!##       Call EXGEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%Y, iErr)
!!!##    End Do
!!!##    
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Read_EXO_Result_Vect2D_Elems
!!!##
!!!##  Subroutine Read_EXO_Result_Mat2D_Elems(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    PetscInt                                       :: Idx
!!!##    PetscInt                                       :: TS
!!!##    Type(Mat2D), Dimension(:), Pointer             :: Res
!!!##
!!!##    PetscInt                                       :: iBlk
!!!##    PetscInt                                       :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##
!!!##    Do iBlk = 1, Geom%Num_Elem_Blks
!!!##       Call EXGEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XX, iErr)
!!!##       Call EXGEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XY, iErr)
!!!##       Call EXGEV(Geom%exoid, TS, Idx+2, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YX, iErr)
!!!##       Call EXGEV(Geom%exoid, TS, Idx+3, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YY, iErr)
!!!##    End Do
!!!##    
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Read_EXO_Result_Mat2D_Elems
!!!##
!!!##  Subroutine Read_EXO_Result_MatS2D_Elems(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    PetscInt                                       :: Idx
!!!##    PetscInt                                       :: TS
!!!##    Type(MatS2D), Dimension(:), Pointer            :: Res
!!!##
!!!##    PetscInt                                       :: iBlk
!!!##    PetscInt                                       :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##
!!!##    Do iBlk = 1, Geom%Num_Elem_Blks
!!!##       Call EXGEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XX, iErr)
!!!##       Call EXGEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YY, iErr)
!!!##       Call EXGEV(Geom%exoid, TS, Idx+2, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XY, iErr)
!!!##    End Do
!!!##    
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Read_EXO_Result_MatS2D_Elems
!!!##
!!!##  Subroutine Read_EXO_Result_Vect3D_Elems(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    PetscInt                                       :: Idx
!!!##    PetscInt                                       :: TS
!!!##    Type(Vect3D), Dimension(:), Pointer            :: Res
!!!##
!!!##    PetscInt                                       :: iBlk
!!!##    PetscInt                                       :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##
!!!##    Do iBlk = 1, Geom%Num_Elem_Blks
!!!##       Call EXGEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%X, iErr)
!!!##       Call EXGEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%Y, iErr)
!!!##    End Do
!!!##    
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Read_EXO_Result_Vect3D_Elems
!!!##
!!!##  Subroutine Read_EXO_Result_Mat3D_Elems(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    PetscInt                                       :: Idx
!!!##    PetscInt                                       :: TS
!!!##    Type(Mat3D), Dimension(:), Pointer             :: Res
!!!##
!!!##    PetscInt                                       :: iBlk
!!!##    PetscInt                                       :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##
!!!##    Do iBlk = 1, Geom%Num_Elem_Blks
!!!##       Call EXGEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XX, iErr)
!!!##       Call EXGEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XY, iErr)
!!!##       Call EXGEV(Geom%exoid, TS, Idx+2, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XZ, iErr)
!!!##       Call EXGEV(Geom%exoid, TS, Idx+3, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YX, iErr)
!!!##       Call EXGEV(Geom%exoid, TS, Idx+4, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YY, iErr)
!!!##       Call EXGEV(Geom%exoid, TS, Idx+5, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YZ, iErr)
!!!##       Call EXGEV(Geom%exoid, TS, Idx+6, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%ZX, iErr)
!!!##       Call EXGEV(Geom%exoid, TS, Idx+7, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%ZY, iErr)
!!!##       Call EXGEV(Geom%exoid, TS, Idx+8, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%ZZ, iErr)
!!!##    End Do
!!!##    
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Read_EXO_Result_Mat3D_Elems
!!!##
!!!##  Subroutine Read_EXO_Result_MatS3D_Elems(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    PetscInt                                       :: Idx
!!!##    PetscInt                                       :: TS
!!!##    Type(MatS3D), Dimension(:), Pointer            :: Res
!!!##
!!!##    PetscInt                                       :: iBlk
!!!##    PetscInt                                       :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##
!!!##    Do iBlk = 1, Geom%Num_Elem_Blks
!!!##       Call EXGEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XX, iErr)
!!!##       Call EXGEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YY, iErr)
!!!##       Call EXGEV(Geom%exoid, TS, Idx+2, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%ZZ, iErr)
!!!##       Call EXGEV(Geom%exoid, TS, Idx+3, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XY, iErr)
!!!##       Call EXGEV(Geom%exoid, TS, Idx+4, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YZ, iErr)
!!!##       Call EXGEV(Geom%exoid, TS, Idx+5, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XZ, iErr)
!!!##    End Do
!!!##    
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Read_EXO_Result_MatS3D_Elems

!!! WRITE
   Subroutine Write_EXO_Result_Global(dEXO, dIdx, dTS, dRes)
      Type (EXO_Info), Intent(INOUT)                 :: dEXO
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      PetscReal                                      :: dRes
      
      PetscInt                                       :: iErr
      PetscReal                                      :: Vers
      PetscReal, Dimension(:), Pointer               :: Tmp_Res
      PetscInt                                       :: Num_Vars, Num_TS
      PetscReal                                      :: fDum
      Character                                      :: cDum
      
      If ( ((dEXO%comm == PETSC_COMM_WORLD) .AND. (MEF90_MyRank == 0)) .OR. (dEXO%comm == PETSC_COMM_SELF) ) Then
         dEXO%exoid = EXOPEN(dEXO%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers, ierr)
         ! Get the number of global variables stored in the database    
         Call EXGVP(dEXO%exoid, 'G', Num_Vars, iErr)
         Allocate(Tmp_Res(Num_Vars))
         Tmp_Res = 0.0_Kr
         
         ! Read All the global variables at time step TS into Tmp_Res
         ! Modify Tmp_Res(Idx) and write everything back...
         Call EXGGV(dEXO%exoid, dTS, Num_Vars, Tmp_Res, iErr)
         Tmp_Res(dIdx) = dRes
         
         Call EXPGV(dEXO%exoid, dTS, Num_Vars, Tmp_Res, iErr)
         DeAllocate (Tmp_Res)
         Call EXCLOS(dEXO%exoid, iErr)
         dEXO%exoid = 0
      End If
   End Subroutine Write_EXO_Result_Global

   Subroutine Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Info), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Info)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      PetscReal, Dimension(:), Pointer               :: dRes
      
      PetscInt                                       :: Num_Rec, iRec
      PetscInt                                       :: iErr
      PetscReal                                      :: Vers
      
      dEXO%exoid = EXOPEN(dEXO%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers, ierr)
      
      If (Mod(Size(dRes), dMeshTopology%Num_Verts) /= 0) Then
         SETERRQ(PETSC_ERR_ARG_SIZ, 'Write_EXO_Result_VertexPtrInterlaced: The argument does not match the number of vertices in the mesh', iErr)
      End If
      Num_Rec = Size(dRes) / dMeshTopology%Num_Verts
      Do iRec = 1, Num_Rec
         Call EXPNV(dEXO%exoid, dTS, dIdx + iRec-1, dMeshTopology%Num_Verts, dRes(iRec:dMeshTopology%Num_Verts*Num_Rec:Num_Rec), iErr); CHKERRQ(iErr)
      End Do
      Call EXCLOS(dEXO%exoid, iErr)
      dEXO%exoid = 0
   End Subroutine Write_EXO_Result_VertexPtrInterlaced

   Subroutine Write_EXO_Result_VertexSection(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Info), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Info)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (SectionReal)                             :: dRes
      
      Type (Vec)                                     :: Res_Vec
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr

      Call SectionRealCreateLocalVector(dRes, Res_Vec, iErr); CHKERRQ(iErr)
      Call VecGetArrayF90(Res_Vec, Res_Ptr, iErr); CHKERRQ(iErr)
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx, dTS, Res_Ptr)
      Call VecRestoreArrayF90(Res_Vec, Res_Ptr, iErr); CHKERRQ(iErr)
      Call VecDestroy(Res_Vec, iErr); CHKERRQ(iErr)
   End Subroutine Write_EXO_Result_VertexSection

   Subroutine Write_EXO_Result_VertexVec(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Info), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Info)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Vec)                                     :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr

      Call VecGetArrayF90(dRes, Res_Ptr, iErr); CHKERRQ(iErr)
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx, dTS, Res_Ptr)
      Call VecRestoreArrayF90(dRes, Res_Ptr, iErr); CHKERRQ(iErr)
   End Subroutine Write_EXO_Result_VertexVec

   Subroutine Write_EXO_Result_VertexVect2D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Info), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Info)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Vect2D), Dimension(:), Pointer           :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr

      If ( Size(dRes) /= dMeshTopology%Num_Verts ) Then
         SETERRQ(PETSC_ERR_ARG_SIZ, 'Write_EXO_Result_VertexVect2D: The argument size does not match the number of vertices in the mesh', iErr)
      End If
      Allocate(Res_Ptr(dMeshTopology%Num_Verts))
      Res_Ptr = dRes(:)%X
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%Y
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+1, dTS, Res_Ptr)
      DeAllocate(Res_Ptr)
   End Subroutine Write_EXO_Result_VertexVect2D

   Subroutine Write_EXO_Result_VertexVect3D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Info), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Info)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Vect3D), Dimension(:), Pointer           :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr

      If ( Size(dRes) /= dMeshTopology%Num_Verts ) Then
         SETERRQ(PETSC_ERR_ARG_SIZ, 'Write_EXO_Result_VertexVect2D: The argument size does not match the number of vertices in the mesh', iErr)
      End If
      Allocate(Res_Ptr(dMeshTopology%Num_Verts))
      Res_Ptr = dRes(:)%X
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%Y
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+1, dTS, Res_Ptr)
      Res_Ptr = dRes(:)%Z
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+2, dTS, Res_Ptr)
      DeAllocate(Res_Ptr)
   End Subroutine Write_EXO_Result_VertexVect3D


!!!##  Subroutine Write_EXO_Result_Mat2D_Nodes(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    PetscInt                                       :: Idx
!!!##    PetscInt                                       :: TS
!!!##    Type(Mat2D), Dimension(:), Pointer            :: Res
!!!##
!!!##    PetscInt                                       :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##
!!!##!!$    If (Geom%Num_Dim /= 2) Then 
!!!##!!$       Write(*,*) Trim(Geom%Filename), 'is not a 2D File.'
!!!##!!$       Write(*,*) 'Aborting'
!!!##!!$       STOP
!!!##!!$    End If
!!!##    Call EXPNV(Geom%exoid, TS, Idx,   Geom%Num_Nodes, Res%XX, iErr)
!!!##    Call EXPNV(Geom%exoid, TS, Idx+1, Geom%Num_Nodes, Res%XY, iErr)
!!!##    Call EXPNV(Geom%exoid, TS, Idx+2, Geom%Num_Nodes, Res%YX, iErr)
!!!##    Call EXPNV(Geom%exoid, TS, Idx+3, Geom%Num_Nodes, Res%YY, iErr)
!!!##    
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Write_EXO_Result_Mat2D_Nodes
!!!##
!!!##  Subroutine Write_EXO_Result_MatS2D_Nodes(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    PetscInt                                       :: Idx
!!!##    PetscInt                                       :: TS
!!!##    Type(MatS2D), Dimension(:), Pointer            :: Res
!!!##
!!!##    PetscInt                                       :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##
!!!##!!$    If (Geom%Num_Dim /= 2) Then 
!!!##!!$       Write(*,*) Trim(Geom%Filename), 'is not a 2D File.'
!!!##!!$       Write(*,*) 'Aborting'
!!!##!!$       STOP
!!!##!!$    End If
!!!##    Call EXPNV(Geom%exoid, TS, Idx,   Geom%Num_Nodes, Res%XX, iErr)
!!!##    Call EXPNV(Geom%exoid, TS, Idx+1, Geom%Num_Nodes, Res%YY, iErr)
!!!##    Call EXPNV(Geom%exoid, TS, Idx+2, Geom%Num_Nodes, Res%XY, iErr)
!!!##    
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Write_EXO_Result_MatS2D_Nodes
!!!##
!!!##  Subroutine Write_EXO_Result_Vect3D_Nodes(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    PetscInt                                       :: Idx
!!!##    PetscInt                                       :: TS
!!!##    Type(Vect3D), Dimension(:), Pointer            :: Res
!!!##
!!!##    PetscInt                                       :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##
!!!##!!$    If (Geom%Num_Dim /= 3) Then 
!!!##!!$       Write(*,*) Trim(Geom%Filename), 'is not a 3D File.'
!!!##!!$       Write(*,*) 'Aborting'
!!!##!!$       STOP
!!!##!!$    End If
!!!##    Call EXPNV(Geom%exoid, TS, Idx,   Geom%Num_Nodes, Res%X, iErr)
!!!##    Call EXPNV(Geom%exoid, TS, Idx+1, Geom%Num_Nodes, Res%Y, iErr)
!!!##    Call EXPNV(Geom%exoid, TS, Idx+2, Geom%Num_Nodes, Res%Z, iErr)
!!!##    
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Write_EXO_Result_Vect3D_Nodes
!!!##
!!!##  Subroutine Write_EXO_Result_Mat3D_Nodes(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    PetscInt                                       :: Idx
!!!##    PetscInt                                       :: TS
!!!##    Type(Mat3D), Dimension(:), Pointer             :: Res
!!!##
!!!##    PetscInt                                       :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##
!!!##!!$    If (Geom%Num_Dim /= 3) Then 
!!!##!!$       Write(*,*) Trim(Geom%Filename), 'is not a 3D File.'
!!!##!!$       Write(*,*) 'Aborting'
!!!##!!$       STOP
!!!##!!$    End If
!!!##    Call EXPNV(Geom%exoid, TS, Idx,   Geom%Num_Nodes, Res%XX, iErr)
!!!##    Call EXPNV(Geom%exoid, TS, Idx+1, Geom%Num_Nodes, Res%XY, iErr)
!!!##    Call EXPNV(Geom%exoid, TS, Idx+2, Geom%Num_Nodes, Res%XZ, iErr)
!!!##    Call EXPNV(Geom%exoid, TS, Idx+3, Geom%Num_Nodes, Res%YX, iErr)
!!!##    Call EXPNV(Geom%exoid, TS, Idx+4, Geom%Num_Nodes, Res%YY, iErr)
!!!##    Call EXPNV(Geom%exoid, TS, Idx+5, Geom%Num_Nodes, Res%YZ, iErr)
!!!##    Call EXPNV(Geom%exoid, TS, Idx+6, Geom%Num_Nodes, Res%ZX, iErr)
!!!##    Call EXPNV(Geom%exoid, TS, Idx+7, Geom%Num_Nodes, Res%ZY, iErr)
!!!##    Call EXPNV(Geom%exoid, TS, Idx+8, Geom%Num_Nodes, Res%ZZ, iErr)
!!!##
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Write_EXO_Result_Mat3D_Nodes
!!!##
!!!##  Subroutine Write_EXO_Result_MatS3D_Nodes(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    PetscInt                                       :: Idx
!!!##    PetscInt                                       :: TS
!!!##    Type(MatS3D), Dimension(:), Pointer            :: Res
!!!##
!!!##    PetscInt                                       :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##
!!!##!!$    If (Geom%Num_Dim /= 3) Then 
!!!##!!$       Write(*,*) Trim(Geom%Filename), 'is not a 3D File.'
!!!##!!$       Write(*,*) 'Aborting'
!!!##!!$       STOP
!!!##!!$    End If
!!!##    Call EXPNV(Geom%exoid, TS, Idx,   Geom%Num_Nodes, Res%XX, iErr)
!!!##    Call EXPNV(Geom%exoid, TS, Idx+1, Geom%Num_Nodes, Res%YY, iErr)
!!!##    Call EXPNV(Geom%exoid, TS, Idx+2, Geom%Num_Nodes, Res%ZZ, iErr)
!!!##    Call EXPNV(Geom%exoid, TS, Idx+3, Geom%Num_Nodes, Res%XY, iErr)
!!!##    Call EXPNV(Geom%exoid, TS, Idx+4, Geom%Num_Nodes, Res%YZ, iErr)
!!!##    Call EXPNV(Geom%exoid, TS, Idx+5, Geom%Num_Nodes, Res%XZ, iErr)
!!!##
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Write_EXO_Result_MatS3D_Nodes
!!!##
!!!##  Subroutine Write_EXO_Result_Scal_Elems(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    PetscInt                                       :: Idx
!!!##    PetscInt                                       :: TS
!!!##    Real(Kind = Kr), Dimension(:), Pointer         :: Res
!!!##
!!!##    Real(Kind = Kr), Dimension(:), Pointer         :: Tmp_Res
!!!##    PetscInt                                       :: iBlk, iE
!!!##    PetscInt                                       :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##
!!!##    Do iBlk = 1, Geom%Num_Elem_Blks
!!!##       Allocate (Tmp_Res(Geom%Elem_Blk(iBlk)%Num_Elems))
!!!##       Do iE = 1, Geom%Elem_Blk(iBlk)%Num_Elems
!!!##          Tmp_Res(iE) = Res(Geom%Elem_Blk(iBlk)%Elem_ID(iE))
!!!##       End Do
!!!##
!!!##       Call EXPEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
!!!##            & Tmp_Res, iErr)
!!!##       DeAllocate(Tmp_Res)
!!!##    End Do
!!!##
!!!##    
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Write_EXO_Result_Scal_Elems
!!!##
!!!##  Subroutine Write_EXO_Result_Vect2D_Elems(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    PetscInt                                       :: Idx
!!!##    PetscInt                                       :: TS
!!!##    Type(Vect2D), Dimension(:), Pointer            :: Res
!!!##
!!!##    PetscInt                                       :: iBlk
!!!##    PetscInt                                       :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##
!!!##    Do iBlk = 1, Geom%Num_Elem_Blks
!!!##       Call EXPEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%X, iErr)
!!!##       Call EXPEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%Y, iErr)
!!!##    End Do
!!!##    
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Write_EXO_Result_Vect2D_Elems
!!!##
!!!##  Subroutine Write_EXO_Result_Mat2D_Elems(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    PetscInt                                       :: Idx
!!!##    PetscInt                                       :: TS
!!!##    Type(Mat2D), Dimension(:), Pointer             :: Res
!!!##
!!!##    PetscInt                                       :: iBlk
!!!##    PetscInt                                       :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##
!!!##    Do iBlk = 1, Geom%Num_Elem_Blks
!!!##       Call EXPEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XX, iErr)
!!!##       Call EXPEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XY, iErr)
!!!##       Call EXPEV(Geom%exoid, TS, Idx+2, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YX, iErr)
!!!##       Call EXPEV(Geom%exoid, TS, Idx+3, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YY, iErr)
!!!##    End Do
!!!##    
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Write_EXO_Result_Mat2D_Elems
!!!##
!!!##  Subroutine Write_EXO_Result_MatS2D_Elems(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    PetscInt                                       :: Idx
!!!##    PetscInt                                       :: TS
!!!##    Type(MatS2D), Dimension(:), Pointer            :: Res
!!!##
!!!##    PetscInt                                       :: iBlk
!!!##    PetscInt                                       :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##
!!!##    Do iBlk = 1, Geom%Num_Elem_Blks
!!!##       Call EXPEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XX, iErr)
!!!##       Call EXPEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YY, iErr)
!!!##       Call EXPEV(Geom%exoid, TS, Idx+2, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XY, iErr)
!!!##    End Do
!!!##    
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Write_EXO_Result_MatS2D_Elems
!!!##
!!!##  Subroutine Write_EXO_Result_Vect3D_Elems(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    PetscInt                                       :: Idx
!!!##    PetscInt                                       :: TS
!!!##    Type(Vect3D), Dimension(:), Pointer            :: Res
!!!##
!!!##    Type(Vect3D), DImension(:), Pointer            :: Tmp_Res
!!!##    PetscInt                                       :: iBlk, iE
!!!##    PetscInt                                       :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##
!!!##    Do iBlk = 1, Geom%Num_Elem_Blks
!!!##       Allocate (Tmp_Res(Geom%Elem_Blk(iBlk)%Num_Elems))
!!!##       Do iE = 1, Geom%Elem_Blk(iBlk)%Num_Elems
!!!##          Tmp_Res(iE) = Res(Geom%Elem_Blk(iBlk)%Elem_ID(iE))
!!!##       End Do
!!!##       Call EXPEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
!!!##            & Tmp_Res%X, iErr)
!!!##       Call EXPEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Tmp_Res%Y, iErr)
!!!##       Call EXPEV(Geom%exoid, TS, Idx+2, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Tmp_Res%Z, iErr)
!!!##       DeAllocate(Tmp_Res)
!!!##    End Do
!!!##    
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Write_EXO_Result_Vect3D_Elems
!!!##
!!!##  Subroutine Write_EXO_Result_Mat3D_Elems(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    PetscInt                                       :: Idx
!!!##    PetscInt                                       :: TS
!!!##    Type(Mat3D), Dimension(:), Pointer             :: Res
!!!##
!!!##    PetscInt                                       :: iBlk
!!!##    PetscInt                                       :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##
!!!##    Do iBlk = 1, Geom%Num_Elem_Blks
!!!##       Call EXPEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XX, iErr)
!!!##       Call EXPEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XY, iErr)
!!!##       Call EXPEV(Geom%exoid, TS, Idx+2, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XZ, iErr)
!!!##       Call EXPEV(Geom%exoid, TS, Idx+3, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YX, iErr)
!!!##       Call EXPEV(Geom%exoid, TS, Idx+4, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YY, iErr)
!!!##       Call EXPEV(Geom%exoid, TS, Idx+5, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YZ, iErr)
!!!##       Call EXPEV(Geom%exoid, TS, Idx+6, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%ZX, iErr)
!!!##       Call EXPEV(Geom%exoid, TS, Idx+7, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%ZY, iErr)
!!!##       Call EXPEV(Geom%exoid, TS, Idx+8, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%ZZ, iErr)
!!!##    End Do
!!!##    
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Write_EXO_Result_Mat3D_Elems
!!!##
!!!##  Subroutine Write_EXO_Result_MatS3D_Elems(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    PetscInt                                       :: Idx
!!!##    PetscInt                                       :: TS
!!!##    Type(MatS3D), Dimension(:), Pointer            :: Res
!!!##
!!!##    Type(MatS3D), Dimension(:), Pointer            :: Tmp_Res
!!!##    PetscInt                                       :: iBlk, iE
!!!##    PetscInt                                       :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##
!!!##    Do iBlk = 1, Geom%Num_Elem_Blks
!!!##       Allocate (Tmp_Res(Geom%Elem_Blk(iBlk)%Num_Elems))
!!!##       Do iE = 1, Geom%Elem_Blk(iBlk)%Num_Elems
!!!##          Tmp_Res(iE) = Res(Geom%Elem_Blk(iBlk)%Elem_ID(iE))
!!!##       End Do
!!!##       Call EXPEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
!!!##            & Tmp_Res%XX, iErr)
!!!##       Call EXPEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Tmp_Res%YY, iErr)
!!!##       Call EXPEV(Geom%exoid, TS, Idx+2, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Tmp_Res%ZZ, iErr)
!!!##       Call EXPEV(Geom%exoid, TS, Idx+3, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Tmp_Res%XY, iErr)
!!!##       Call EXPEV(Geom%exoid, TS, Idx+4, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Tmp_Res%YZ, iErr)
!!!##       Call EXPEV(Geom%exoid, TS, Idx+5, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##            & Tmp_Res%XZ, iErr)
!!!##!!$
!!!##!!$       Call EXPEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
!!!##!!$            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XX, iErr)
!!!##!!$       Call EXPEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##!!$            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YY, iErr)
!!!##!!$       Call EXPEV(Geom%exoid, TS, Idx+2, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##!!$            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%ZZ, iErr)
!!!##!!$       Call EXPEV(Geom%exoid, TS, Idx+3, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##!!$            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XY, iErr)
!!!##!!$       Call EXPEV(Geom%exoid, TS, Idx+4, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##!!$            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YZ, iErr)
!!!##!!$       Call EXPEV(Geom%exoid, TS, Idx+5, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!!##!!$            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XZ, iErr)
!!!##       DeAllocate(Tmp_Res)
!!!##    End Do
!!!##    
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Write_EXO_Result_MatS3D_Elems
  
End Module m_MEF_EXO

