Module m_MEF_EXO
!!! ATTENTION:
!!! DoF ordering has changed. the order is now 1%X, 2%X, ..., 1%Y, 2%Y...
!!! Instead of 1%X, 1%Y, 2%X, 2%Y...
!!! I/O routines have to be changed accordingly...
!!! This is consistent with ensight ordering scheme, and easier to implement

!!! NOT ANYMORE...
!!! Uses Geom%Numbering to know wich scheme is used.
!!! Numbering_PerNodes -> 1%X, 1%Y, 1%Z, 2%X, 2%Y, 2%Z, ...
!!! Numbering_PerCoord -> 1%X, 2%X, ..., 1%Y, 2%Y, ..., 1%Z, 2%Z

   Use m_AlgebLin
   Use m_Constantes
   Use m_MEF_Types
   Use m_MEF_MPI
   Use m_MEF_Elements
   
   IMPLICIT NONE
   Private
   include "exodusII.inc"
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscsys.h"
#include "finclude/petscviewer.h"
#include "finclude/petscviewer.h90"
#include "finclude/petscmesh.h"
#include "finclude/petscmesh.h90"

   Integer, Parameter, Public                        :: exo_cpu_ws = 8
   Integer, Parameter, Public                        :: exo_io_ws = 8
   
   Public :: Write_MeshTopology
   Public :: Write_MeshTopologyGlobal

!   Public :: Show_MeshTopology_Info
!   Public :: Destroy_MeshTopology_Info
!   Public :: Read_EXO_Result_Nodes
!   Public :: Write_EXO_Result_Nodes
!   Public :: Read_EXO_Result_Ptr_Nodes   
!   Public :: Write_EXO_Result_Ptr_Nodes   
!   Public :: Read_EXO_Result_Elems
!   Public :: Write_EXO_Result_Elems   
!   Public :: Read_EXO_Result_Ptr_Elems
!   Public :: Write_EXO_Result_Ptr_Elems   

!   Public :: Read_EXO_Result_Global
!   Public :: Write_EXO_Result_Global   

   
!   Interface Read_MeshTopology_Info
!      Module Procedure Read_MeshTopology_Info_Seq, Read_MeshTopology_Info_Dist
!   End Interface
   
!   Interface Read_EXO_Node_Coord
!      Module Procedure Read_EXO_Node_Coord_1D, Read_EXO_Node_Coord_2D, Read_EXO_Node_Coord_3D
!   End Interface
   
!   Interface Write_EXO_Node_Coord
!      Module Procedure Write_EXO_Node_Coord_1D, Write_EXO_Node_Coord_2D, Write_EXO_Node_Coord_3D
!   End Interface
   
!   Interface Read_EXO_Connect
!      Module Procedure  Read_EXO_Connect_2D, Read_EXO_Connect_2D_Scal, Read_EXO_Connect_3D, Read_EXO_Connect_3D_Scal, Read_EXO_Connect_3D_Elast
!   End Interface
   
!   Interface Write_EXO_Connect
!      Module Procedure  Write_EXO_Connect_2D, Write_EXO_Connect_2D_Scal, Write_EXO_Connect_3D, Write_EXO_Connect_3D_Scal, Write_EXO_Connect_3D_Elast
!   End Interface
   
!   Interface Read_EXO_Result_Nodes
!      Module Procedure Read_EXO_Result_Scal_Nodes, Read_EXO_Result_Ptr_Nodes, Read_EXO_Result_Vect2D_Nodes, Read_EXO_Result_Mat2D_Nodes, Read_EXO_Result_MatS2D_Nodes, Read_EXO_Result_Vect3D_Nodes, Read_EXO_Result_Mat3D_Nodes, Read_EXO_Result_MatS3D_Nodes, Read_EXO_Result_Vec_Nodes
!   End Interface
!   
!   Interface Write_EXO_Result_Nodes
!      Module Procedure Write_EXO_Result_Ptr_Nodes, Write_EXO_Result_Vect2D_Nodes, Write_EXO_Result_Mat2D_Nodes, Write_EXO_Result_MatS2D_Nodes, Write_EXO_Result_Vect3D_Nodes, Write_EXO_Result_Mat3D_Nodes, Write_EXO_Result_MatS3D_Nodes, Write_EXO_Result_Vec_Nodes
!   End Interface
!   
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
      Integer                          :: UniqCount
      PetscMPIInt                      :: rank
      Integer                          :: i, j, iErr

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
      Integer                                        :: vers
      Integer                                        :: iErr
      Integer                                        :: iDummy
      PetscReal                                      :: rDummy
      Character                                      :: cDummy
      
      Integer                                        :: iBlk, iSet, i
      Integer                                        :: iE, iELoc
      Integer                                        :: Num_Attr = 0         
      Integer                                        :: Num_Dist_Factor = 0
      
      Character(len=MXSTLN), Dimension(3)            :: Coord_Names, Elem_Type
      PetscReal, Dimension(:,:), Pointer             :: Coordinates
      Integer, Dimension(:,:), Pointer               :: ConnectMesh
      Integer, Dimension(:), Pointer                 :: ConnectBlk
      Integer                                        :: offset
      
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
      
      Integer                                        :: iBlk, iSet, i, iErr
      Integer                                        :: iE, iELoc
      Integer                                        :: Num_Attr = 0         
      Integer                                        :: Num_Dist_Factor = 0
      
      Type(MeshTopology_Info)                        :: GlobalMeshTopology
      Integer, Dimension(:), Pointer                 :: Tmp_GlobalID, Tmp_ID


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
      Integer                                        :: dIdx
      Integer                                        :: dTS
      PetscReal                                      :: dRes
      Integer                                        :: dNumRec
      
      PetscReal                                      :: MyRes
      Integer                                        :: iErr
      PetscReal                                      :: Vers
      PetscReal, Dimension(:), Pointer               :: Tmp_Res
      Integer                                        :: Num_Vars, Num_TS
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

   Subroutine Read_EXO_Result_Section_Nodes(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Info), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Info)                       :: dMeshTopology
      Integer                                        :: dIdx
      Integer                                        :: dTS
      SectionReal                                    :: dRes
      
      Vec                                            :: Res_Vec
      PetscReal, Dimension(:), Pointer               :: Res_Ptr, Res_Comp_Ptr
      Integer                                        :: Num_Rec, iRec
      Integer                                        :: iErr
      PetscReal                                      :: Vers
      
      dEXO%exoid = EXOPEN(dEXO%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers, ierr)
      
      !!! We Assume that the section is initialized and has the proper size
      Call SectionRealCreateLocalVector(dRes, Res_Vec, iErr); CHKERRQ(iErr)
      Call VecGetArrayF90(Res_Vec, Res_Ptr, iErr); CHKERRQ(iErr)
      If (Mod(Size(Res_Ptr), dMeshTopology%Num_Verts) /= 0) Then
         SETERRQ(PETSC_ERR_ARG_SIZ, 'The Section does not match the number of dof in the mesh', iErr)
      End If
      Num_Rec = Size(Res_Ptr) / dMeshTopology%Num_Verts
      
      Allocate(Res_Comp_Ptr(dMeshTopology%Num_Verts))
      Do iRec = dIdx, dIdx+Num_Rec
         Call EXGNV(dEXO%exoid, dTS, iRec, dMeshTopology%Num_Verts, Res_Comp_Ptr, iErr)
         Res_Ptr((iRec-1)*dMeshTopology%Num_Verts: iRec*dMeshTopology%Num_Verts) = Res_Comp_Ptr
         !!! WRONG! this is not how the data is stored in Res_Ptr
         !!!        do as below
      End Do
      DeAllocate(Res_Comp_Ptr)      
      Call VecRestoreArrayF90(Res_Vec, Res_Ptr, iErr); CHKERRQ(iErr)
      !!! How do I release Res_Vec?
      
      Call EXCLOS(dEXO%exoid, iErr)
      dEXO%exoid = 0
   End Subroutine Read_EXO_Result_Section_Nodes


!!!##
!!!##   Subroutine Read_EXO_Result_Vec_Nodes(Geom, Layout, Idx, TS, Res)
!!!##      Type (MeshTopology_Info), Intent(INOUT)              :: Geom
!!!##      Type (Layout_Info), Intent(IN)                   :: Layout 
!!!##      Integer, Intent(IN)                              :: Idx
!!!##      Integer, Intent(IN)                              :: TS
!!!##      Vec                                              :: Res
!!!##                                                     
!!!##      Vec                                              :: Res_IO
!!!##      PetscReal, Dimension(:), Pointer                 :: Res_Array
!!!##      Integer                                          :: iErr
!!!##      Real(Kind = Kr)                                  :: Vers
!!!##      Integer                                          :: Num_Rec, Res_Size, iRec
!!!##      Integer                                          :: MyRank, NumProcs, IO_Size
!!!##
!!!##      Call VecGetSize(Res, Res_Size, iErr)
!!!##      If ( Mod(Res_Size, Geom%Num_Nodes) /= 0) Then
!!!##         Write(*,*)  '[ERROR] Cannot find the number of records to save'
!!!##         Write(*,*)  '        Mod(Res_Size, Geom%Num_Nodes) = ', Mod(Res_Size, Num_Rec) 
!!!##         Write(*,*)  '        Remember that Distributed IO requires the local form of the vector?'
!!!##         RETURN
!!!##      Else
!!!##         Num_Rec = Res_Size / Geom%Num_Nodes
!!!##      End If
!!!##      If (Num_Rec >1) Then
!!!##         Write(*,*) '[WARNING] not tested yet'
!!!##      End If      
!!!##         
!!!##      Call MPI_Comm_size(Geom%Comm, NumProcs, iErr)
!!!##      Select Case (NumProcs)
!!!##      Case(1)
!!!##         !!! We are doing Distributed IO: each CPU reads from its own EXO file
!!!##         Write(*,*) '[WARNING] not tested yet'
!!!##         Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers, ierr)
!!!##         Call VecGetArrayF90(Res, Res_Array, iErr)
!!!##         Do iRec = 0, Num_Rec-1
!!!##            Call EXGNV(Geom%exoid, TS, Idx + iRec, Res_Size, Res_Array(iRec:Geom%Num_Nodes*Num_Rec:Num_Rec), iErr)
!!!##         End Do
!!!##         Call VecRestoreArrayF90(Res, Res_Array, iErr)
!!!##         Call EXCLOS(Geom%exoid, iErr)
!!!##         Geom%exoid = 0
!!!##
!!!##      Case Default
!!!##         !!! We are doing Sequential IO for a MPI or Seq Vec
!!!##         Call MPI_Comm_rank(Geom%Comm, MyRank, iErr)
!!!##         If (MyRank == 0) Then
!!!##            IO_Size = Res_Size
!!!##         Else
!!!##            IO_Size = 0
!!!##         End If    
!!!##         Call VecCreateMPI(Geom%Comm, IO_Size, Res_Size, Res_IO, iErr)
!!!##
!!!##         !!! Read the fields, component by component
!!!##         If (MyRank == 0) Then
!!!##            Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers, ierr)
!!!##            Call VecGetArrayF90(Res_IO, Res_Array, iErr)
!!!##            Do iRec = 0, Num_Rec-1
!!!##               Call EXGNV(Geom%exoid, TS, Idx + iRec, Res_Size, Res_Array(1+iRec * Geom%Num_Nodes:(iRec+1) * Geom%Num_Nodes), iErr)
!!!##            End Do
!!!##            Call VecRestoreArrayF90(Res_IO, Res_Array, iErr)
!!!##            Call EXCLOS(Geom%exoid, iErr)
!!!##            Geom%exoid = 0
!!!##         End If
!!!##         
!!!##         !!! Scatter the Vec from CPU 0
!!!##         Call VecScatterBegin(Layout%ToIOSeq_N, Res_IO, Res, INSERT_VALUES, SCATTER_REVERSE, iErr) 
!!!##         Call VecScatterEnd  (Layout%ToIOSeq_N, Res_IO, Res, INSERT_VALUES, SCATTER_REVERSE, iErr) 
!!!##
!!!##
!!!##      End Select
!!!##   End Subroutine Read_EXO_Result_Vec_Nodes
!!!##
!!!##  Subroutine Read_EXO_Result_Ptr_Nodes(Geom, Idx, TS, Res, Num_Rec)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Real(Kind = Kr), Dimension(:), Pointer         :: Res
!!!##    Integer                                        :: Num_Rec
!!!##
!!!##    Integer                                        :: iErr, i
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Integer                                        :: iRec
!!!##    Real(Kind = Kr), Dimension(:), Pointer         :: Tmp_Res
!!!##
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##    
!!!##    If (Associated(Res)) Then
!!!##       DeAllocate(Res)
!!!##    End If
!!!##    Allocate(Res(Geom%Num_Nodes*Num_Rec))
!!!##
!!!##    Allocate(Tmp_Res(Geom%Num_Nodes))
!!!##    Select Case (Geom%Numbering)
!!!##    Case(Numbering_PerCoord)
!!!##       Do iRec = 0, Num_Rec-1
!!!##          Call EXGNV(Geom%exoid, TS, Idx + iRec, Geom%Num_Nodes,              &
!!!##               & Res(1+iRec * Geom%Num_Nodes:(iRec+1) * Geom%Num_Nodes), iErr)
!!!##       End Do
!!!##    Case(Numbering_PerNodes)
!!!##       Do iRec = 1, Num_Rec
!!!##          Call EXGNV(Geom%exoid, TS, Idx + iRec-1 , Geom%Num_Nodes,           &
!!!##               & Tmp_Res, iErr)
!!!##!               & Res(iRec:Num_Rec*Geom%Num_Nodes:Num_Rec), iErr)
!!!##          Res(iRec:Num_Rec*Geom%Num_Nodes:Num_Rec) = Tmp_Res
!!!##       End Do
!!!##
!!!##    Case Default
!!!##       Write(*,*) 'ERROR: Read_EXO_Result_Ptr_Nodes Numbering scheme unknown',&
!!!##            & Geom%Numbering
!!!##       STOP
!!!##    End Select
!!!##    DeAllocate(Tmp_Res)
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Read_EXO_Result_Ptr_Nodes
!!!##
!!!##  Subroutine Read_EXO_Result_Vect2D_Nodes(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Type(Vect2D), Dimension(:), Pointer            :: Res
!!!##
!!!##    Integer                                        :: iErr
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
!!!##    Call EXGNV(Geom%exoid, TS, Idx, Geom%Num_Nodes, Res%X, iErr)
!!!##    Call EXGNV(Geom%exoid, TS, Idx+1, Geom%Num_Nodes, Res%Y, iErr)
!!!##    
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Read_EXO_Result_Vect2D_Nodes
!!!##
!!!##  Subroutine Read_EXO_Result_Mat2D_Nodes(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Type(Mat2D), Dimension(:), Pointer             :: Res
!!!##
!!!##    Integer                                        :: iErr
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
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Type(MatS2D), Dimension(:), Pointer            :: Res
!!!##
!!!##    Integer                                        :: iErr
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
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Type(Vect3D), Dimension(:), Pointer            :: Res
!!!##
!!!##    Integer                                        :: iErr
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
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Type(Mat3D), Dimension(:), Pointer             :: Res
!!!##
!!!##    Integer                                        :: iErr
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
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Type(MatS3D), Dimension(:), Pointer            :: Res
!!!##
!!!##    Integer                                        :: iErr
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
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Real(Kind = Kr), Dimension(:), Pointer         :: Res
!!!##
!!!##    Integer                                        :: iBlk
!!!##    Integer                                        :: iErr
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
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Type(Vect2D), Dimension(:), Pointer            :: Res
!!!##
!!!##    Integer                                        :: iBlk
!!!##    Integer                                        :: iErr
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
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Type(Mat2D), Dimension(:), Pointer             :: Res
!!!##
!!!##    Integer                                        :: iBlk
!!!##    Integer                                        :: iErr
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
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Type(MatS2D), Dimension(:), Pointer            :: Res
!!!##
!!!##    Integer                                        :: iBlk
!!!##    Integer                                        :: iErr
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
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Type(Vect3D), Dimension(:), Pointer            :: Res
!!!##
!!!##    Integer                                        :: iBlk
!!!##    Integer                                        :: iErr
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
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Type(Mat3D), Dimension(:), Pointer             :: Res
!!!##
!!!##    Integer                                        :: iBlk
!!!##    Integer                                        :: iErr
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
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Type(MatS3D), Dimension(:), Pointer            :: Res
!!!##
!!!##    Integer                                        :: iBlk
!!!##    Integer                                        :: iErr
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
      Integer                                        :: dIdx
      Integer                                        :: dTS
      PetscReal                                      :: dRes
      
      Integer                                        :: iErr
      PetscReal                                      :: Vers
      PetscReal, Dimension(:), Pointer               :: Tmp_Res
      Integer                                        :: Num_Vars, Num_TS
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


!!!##!!$  Subroutine Write_EXO_Result_Scal_Nodes(Geom, Idx, TS, Res)
!!!##!!$    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##!!$    Integer                                        :: Idx
!!!##!!$    Integer                                        :: TS
!!!##!!$    Real(Kind = Kr), Dimension(:), Pointer         :: Res
!!!##!!$
!!!##!!$    Integer                                        :: iErr
!!!##!!$    Real(Kind = Kr)                                :: Vers
!!!##!!$    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
!!!##!!$         & ierr)
!!!##!!$
!!!##!!$    Call EXPNV(Geom%exoid, TS, Idx, Geom%Num_Nodes, Res, iErr)
!!!##!!$    
!!!##!!$    Call EXCLOS(Geom%exoid, iErr)
!!!##!!$    Geom%exoid = 0
!!!##!!$  End Subroutine Write_EXO_Result_Scal_Nodes
!!!##!!$
!!!##
!!!##   Subroutine Write_EXO_Result_Vec_Nodes(Geom, Layout, Idx, TS, Res)
!!!##      Type (MeshTopology_Info), Intent(INOUT)              :: Geom
!!!##      Type (Layout_Info), Intent(IN)                   :: Layout 
!!!##      Integer, Intent(IN)                              :: Idx
!!!##      Integer, Intent(IN)                              :: TS
!!!##      Vec                                              :: Res
!!!##                                                     
!!!##      Vec                                              :: Res_IO
!!!##      PetscReal, Dimension(:), Pointer                 :: Res_Array
!!!##      Integer                                          :: iErr
!!!##      Real(Kind = Kr)                                  :: Vers
!!!##      Integer                                          :: Num_Rec, Res_Size, iRec
!!!##      Integer                                          :: MyRank, NumProcs, IO_Size
!!!##
!!!##      Call VecGetSize(Res, Res_Size, iErr)
!!!##      If ( Mod(Res_Size, Geom%Num_Nodes) /= 0) Then
!!!##         Write(*,*)  '[ERROR] Cannot find the number of records to save'
!!!##         Write(*,*)  '        Mod(Res_Size, Geom%Num_Nodes) = ', Mod(Res_Size, Num_Rec) 
!!!##         Write(*,*)  '        Remember that Distributed IO requires the local form of the vector?'
!!!##         RETURN
!!!##      Else
!!!##         Num_Rec = Res_Size / Geom%Num_Nodes
!!!##      End If
!!!##      If (Num_Rec >1) Then
!!!##         Write(*,*) '[WARNING] not tested yet'
!!!##      End If      
!!!##         
!!!##      Call MPI_Comm_size(Geom%Comm, NumProcs, iErr)
!!!##      Select Case (NumProcs)
!!!##      Case(1)
!!!##         !!! We are doing Distributed IO: each CPU writes in its own EXO file
!!!##         Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers, ierr)
!!!##         Call VecGetArrayF90(Res, Res_Array, iErr)
!!!##         Do iRec = 0, Num_Rec-1
!!!##            Call EXPNV(Geom%exoid, TS, Idx + iRec, Res_Size, Res_Array(Geom%Num_Nodes*iRec+1:Geom%Num_Nodes*(iRec+1)), iErr)
!!!##         End Do
!!!##         Call VecRestoreArrayF90(Res, Res_Array, iErr)
!!!##         Call EXCLOS(Geom%exoid, iErr)
!!!##         Geom%exoid = 0
!!!##
!!!##      Case Default
!!!##         !!! We are doing Sequential IO for a MPI or Seq Vec
!!!##         Call MPI_Comm_rank(Geom%Comm, MyRank, iErr)
!!!##         If (MyRank == 0) Then
!!!##            IO_Size = Res_Size
!!!##         Else
!!!##            IO_Size = 0
!!!##         End If    
!!!##         Call VecCreateMPI(Geom%Comm, IO_Size, Res_Size, Res_IO, iErr)
!!!##
!!!##         !!! Scatter the Vec onto CPU 0
!!!##         Call VecScatterBegin(Layout%ToIOSeq_N, Res, Res_IO, INSERT_VALUES, SCATTER_FORWARD, iErr) 
!!!##         Call VecScatterEnd  (Layout%ToIOSeq_N, Res, Res_IO, INSERT_VALUES, SCATTER_FORWARD, iErr) 
!!!##
!!!##         !!! Save the fields, component by component
!!!##         If (MyRank == 0) Then
!!!##            Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers, ierr)
!!!##            Call VecGetArrayF90(Res_IO, Res_Array, iErr)
!!!##            Do iRec = 0, Num_Rec-1
!!!##               Call EXPNV(Geom%exoid, TS, Idx + iRec, Res_Size, Res_Array(1+iRec * Geom%Num_Nodes:(iRec+1) * Geom%Num_Nodes), iErr)
!!!##            End Do
!!!##            Call VecRestoreArrayF90(Res_IO, Res_Array, iErr)
!!!##            Call EXCLOS(Geom%exoid, iErr)
!!!##            Geom%exoid = 0
!!!##         End If
!!!##      End Select
!!!##   End Subroutine Write_EXO_Result_Vec_Nodes
!!!##    
!!!##  Subroutine Write_EXO_Result_Ptr_Nodes(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Real(Kind = Kr), Dimension(:), Pointer         :: Res
!!!##
!!!##    Real(Kind = Kr), Dimension(:), Pointer         :: Tmp_Res
!!!##    Integer                                        :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Integer                                        :: Num_Rec, iRec
!!!##
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##    Select Case (Geom%Numbering)
!!!##    Case(Numbering_PerCoord)
!!!##       If ( Mod(Size(Res), Geom%Num_Nodes) == 0) Then
!!!##          Num_Rec = Size(Res) / Geom%Num_Nodes
!!!##          Do iRec = 0, Num_Rec-1
!!!##             Call EXPNV(Geom%exoid, TS, Idx + iRec, Geom%Num_Nodes,           &
!!!##                  & Res(1+iRec * Geom%Num_Nodes:(iRec+1) * Geom%Num_Nodes),   &
!!!##                  &iErr)
!!!##          End Do
!!!##       Else
!!!##          Write(*,*) 'Error in Write_EXO_Result_Scal_Ptr_Nodes'
!!!##          Write(*,*) 'Can''t make any sense of the dimension of the arguments:'
!!!##          Write(*,*) 'Size(Res), Geom%Num_Dim ', Size(Res), Geom%Num_Dim
!!!##          Write(*,*) 'Not saving anything in', Trim(Geom%FileName)
!!!##       End If
!!!##    Case(Numbering_PerNodes)
!!!##       If ( Mod(Size(Res), Geom%Num_Nodes) == 0) Then
!!!##          Allocate(Tmp_Res(Geom%Num_Nodes))
!!!##          Num_Rec = Size(Res) / Geom%Num_Nodes
!!!##          Do iRec = 1, Num_Rec
!!!##             Tmp_Res = Res(iRec:Num_Rec * Geom%Num_Nodes:Num_Rec)
!!!##             Call EXPNV(Geom%exoid, TS, Idx + iRec - 1, Geom%Num_Nodes,       &
!!!##                  & Tmp_Res, iErr)
!!!##          End Do
!!!##          DeAllocate(Tmp_Res)
!!!##       Else
!!!##          Write(*,*) 'Error in Write_EXO_Result_Scal_Ptr_Nodes'
!!!##          Write(*,*) 'Can''t make any sense of the dimension of the arguments:'
!!!##          Write(*,*) 'Size(Res), Geom%Num_Dim ', Size(Res), Geom%Num_Dim
!!!##          Write(*,*) 'Not saving anything in', Trim(Geom%FileName)
!!!##       End If
!!!##    Case Default
!!!##       Write(*,*) 'ERROR:Write_EXO_Result_Ptr_Nodes Numbering scheme unknown',&
!!!##            & Geom%Numbering
!!!##       STOP
!!!##    End Select
!!!##
!!!##!!$    If ( Mod(Size(Res), Geom%Num_Nodes) == 0) Then
!!!##!!$       Num_Rec = Size(Res) / Geom%Num_Nodes
!!!##!!$       Do iRec = 0, Num_Rec-1
!!!##!!$          Call EXPNV(Geom%exoid, TS, Idx + iRec, Geom%Num_Nodes,              &
!!!##!!$               & Res(1+iRec * Geom%Num_Nodes:(iRec+1) * Geom%Num_Nodes), iErr)
!!!##!!$       End Do
!!!##!!$    Else
!!!##!!$       Write(*,*) 'Error in Write_EXO_Result_Scal_Ptr_Nodes'
!!!##!!$       Write(*,*) 'Can''t make any sense of the dimension of the arguments:'
!!!##!!$       Write(*,*) 'Size(Res), Geom%Num_Dim ', Size(Res), Geom%Num_Dim
!!!##!!$       Write(*,*) 'Not saving anything in', Trim(Geom%FileName)
!!!##!!$    End If
!!!##!!$       
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Write_EXO_Result_Ptr_Nodes
!!!##
!!!##  Subroutine Write_EXO_Result_Vect2D_Nodes(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Type(Vect2D), Dimension(:), Pointer            :: Res
!!!##
!!!##    Integer                                        :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##
!!!##!!$    If (Geom%Num_Dim /= 2) Then 
!!!##!!$       Write(*,*) Trim(Geom%Filename), 'is not a 2D File.'
!!!##!!$       Write(*,*) 'Aborting'
!!!##!!$       STOP
!!!##!!$    End If
!!!##    Call EXPNV(Geom%exoid, TS, Idx, Geom%Num_Nodes, Res%X, iErr)
!!!##    Call EXPNV(Geom%exoid, TS, Idx+1, Geom%Num_Nodes, Res%Y, iErr)
!!!##    
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Write_EXO_Result_Vect2D_Nodes
!!!##
!!!##  Subroutine Write_EXO_Result_Mat2D_Nodes(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Type(Mat2D), Dimension(:), Pointer            :: Res
!!!##
!!!##    Integer                                        :: iErr
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
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Type(MatS2D), Dimension(:), Pointer            :: Res
!!!##
!!!##    Integer                                        :: iErr
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
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Type(Vect3D), Dimension(:), Pointer            :: Res
!!!##
!!!##    Integer                                        :: iErr
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
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Type(Mat3D), Dimension(:), Pointer             :: Res
!!!##
!!!##    Integer                                        :: iErr
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
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Type(MatS3D), Dimension(:), Pointer            :: Res
!!!##
!!!##    Integer                                        :: iErr
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
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Real(Kind = Kr), Dimension(:), Pointer         :: Res
!!!##
!!!##    Real(Kind = Kr), Dimension(:), Pointer         :: Tmp_Res
!!!##    Integer                                        :: iBlk, iE
!!!##    Integer                                        :: iErr
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
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Type(Vect2D), Dimension(:), Pointer            :: Res
!!!##
!!!##    Integer                                        :: iBlk
!!!##    Integer                                        :: iErr
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
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Type(Mat2D), Dimension(:), Pointer             :: Res
!!!##
!!!##    Integer                                        :: iBlk
!!!##    Integer                                        :: iErr
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
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Type(MatS2D), Dimension(:), Pointer            :: Res
!!!##
!!!##    Integer                                        :: iBlk
!!!##    Integer                                        :: iErr
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
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Type(Vect3D), Dimension(:), Pointer            :: Res
!!!##
!!!##    Type(Vect3D), DImension(:), Pointer            :: Tmp_Res
!!!##    Integer                                        :: iBlk, iE
!!!##    Integer                                        :: iErr
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
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Type(Mat3D), Dimension(:), Pointer             :: Res
!!!##
!!!##    Integer                                        :: iBlk
!!!##    Integer                                        :: iErr
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
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Type(MatS3D), Dimension(:), Pointer            :: Res
!!!##
!!!##    Type(MatS3D), Dimension(:), Pointer            :: Tmp_Res
!!!##    Integer                                        :: iBlk, iE
!!!##    Integer                                        :: iErr
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

