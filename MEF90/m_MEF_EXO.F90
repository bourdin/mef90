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
      
      Character(len=MXSTLN), Dimension(3)            :: Coord_Names
      PetscReal, Dimension(:,:), Pointer             :: Coordinates
      Integer, Dimension(:,:), Pointer               :: ConnectMesh
      Integer, Dimension(:,:), Pointer               :: ConnectBlk
      
      Coord_Names(1) = 'X'
      Coord_Names(2) = 'Y'
      Coord_Names(3) = 'Z'
   
      Is_IO: If (dEXO%comm == PETSC_COMM_SELF) Then
         ! Open File
         dEXO%exoid = EXCRE (dEXO%filename, EXCLOB, exo_cpu_ws, exo_io_ws, iErr)
         dEXO%exoid = EXOPEN(dEXO%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers, ierr)
         
         ! Writes Global Geometric Parameters
         Call EXPINI(dEXO%exoid, dEXO%Title, dMeshTopology%Num_Dim, dMeshTopology%Num_Vert, dMeshTopology%Num_Elems, dMeshTopology%Num_Elem_Blks, dMeshTopology%Num_Node_Sets, dMeshTopology%Num_Side_Sets, iErr)
         Call EXPCON(dEXO%exoid, Coord_Names, iErr)
         
         ! Write Elem blocks informations
         Do iBlk = 1, dMeshTopology%Num_Elem_Blks
            Call EXPELB(dEXO%exoid, dMeshTopology%elem_blk(iBlk)%ID, dMeshTopology%elem_blk(iBlk)%Elem_Type, dMeshTopology%elem_blk(iBlk)%Num_Elems, dMeshTopology%elem_blk(iBlk)%Num_DoF, Num_Attr, iErr)
         End Do
   
         ! Write Side sets informations
         If (dMeshTopology%Num_Side_Sets > 0) Then
            SETERRQ(PETSC_ERR_SUP, 'Side sets not supported yet', iErr)
         End If
         
         ! Write Node sets informations
         !!! Check that the number of node sets has been compressed like the number of element blocks
         !!! It looks like it doesn;t
         Do iSet = 1, dMeshTopology%Num_Node_Sets
            Call EXPNP(dEXO%exoid, dMeshTopology%Node_Set(iSet)%ID, dMeshTopology%Node_Set(iSet)%Num_Nodes, Num_Dist_Factor, iErr)
            Call EXPNS(dEXO%exoid, dMeshTopology%Node_Set(iSet)%ID, dMeshTopology%Node_Set(iSet)%Node_ID(:), iErr)
         End Do
   
         ! Write vertex coordinates
         Call MeshGetCoordinatesF90(dMeshTopology%mesh, Coordinates, iErr)
         Call EXPCOR(dEXO%exoid, Coordinates(:,1), Coordinates(:,2), Coordinates(:,3), iErr)
         Call MeshRestoreCoordinatesF90(dMeshTopology%mesh, Coordinates, iErr)
         
          ! Write Connectivity tables
          Call MeshGetElementsF90(dMeshTopology%mesh, ConnectMesh, iErr)
          Do iBlk = 1, dMeshTopology%Num_Elem_Blks
             Allocate (ConnectBlk(dMeshTopology%Elem_Blk(iBlk)%Num_Elems, dMeshTopology%Elem_Blk(iBlk)%Num_DoF))
             Do iELoc = 1, dMeshTopology%Elem_Blk(iBlk)%Num_Elems
                iE = dMeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
               ConnectBlk(iELoc, :) = ConnectMesh(iE, :)
 !                There seems to be a problem with the EXPELC call. the block connectivity table looks fine
             End Do
             Call EXPELC (dEXO%exoid, iBlk, ConnectBlk, iErr)
             DeAllocate(ConnectBlk)
          End Do
          Call MeshRestoreElementsF90(dMeshTopology%mesh, ConnectMesh, iErr)

         Call EXCLOS(dEXO%exoid, iErr)
         dEXO%exoid = 0
      Else
         SETERRQ(PETSC_ERR_SUP, 'Synchronized I/O on PETSC_COMM_WORLD not supported yet', ierr)
      End If Is_IO
   End Subroutine Write_MeshTopology
   
   
!!!   Subroutine Write_MeshTopology_Info(Geom)
!!!      Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!      
!!!      Integer                                        :: vers
!!!      Integer                                        :: iErr
!!!      Integer                                        :: iDummy
!!!      Real(Kind = Kr)                                :: rDummy
!!!      Character                                      :: cDummy
!!!      
!!!      Integer                                        :: iBlk
!!!      Integer                                        :: iSet
!!!      Character(len=MXSTLN), Dimension(:,:), Pointer :: Tmp_QA
!!!      Character(len=MXSTLN), Dimension(3)            :: Coord_Names
!!!      Integer                                        :: MyRank
!!!      
!!!      Coord_Names(1) = 'X'
!!!      Coord_Names(2) = 'Y'
!!!      Coord_Names(3) = 'Z'
!!!      
!!!      Call MPI_Comm_rank(Geom%Comm, MyRank, iErr)
!!!      If (MyRank /= 0 ) Then
!!!         RETURN
!!!      End If   
!!!      
!!!      ! Open File
!!!      Geom%exoid = EXCRE (Geom%filename, EXCLOB, exo_cpu_ws, exo_io_ws, iErr)
!!!      Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers, ierr)
!!!    
!!!      ! Writes Global Geometric Parameters
!!!      Call EXPINI(Geom%exoid, Geom%Title, Geom%Num_Dim, Geom%Num_Nodes, Geom%Num_Elems, Geom%Num_Elem_Blks, Geom%Num_Node_Sets, Geom%Num_Side_Sets, iErr)
!!!      Call EXPCON(Geom%exoid, Coord_Names, iErr)
!!!      
!!!      If (Geom%Num_Side_Sets > 0) Then
!!!         Write(*,*) 'WARNING: Side_sets non implemented, continuing anyway'
!!!         Write(*,*) '         Got Geom%Num_Side_Sets = ', Geom%Num_Side_Sets
!!!      End If
!!!    
!!!      ! Write Elem blocks informations
!!!      If (Geom%Num_Elem_blks > 0) Then
!!!         Do iBlk = 1, Geom%Num_Elem_Blks
!!!            Call EXPELB(Geom%exoid, Geom%elem_blk(iBlk)%ID, Geom%elem_blk(iBlk)%Type, Geom%elem_blk(iBlk)%Num_Elems, Geom%elem_blk(iBlk)%Num_Nodes_Per_Elem, Geom%elem_blk(iBlk)%Num_Attr, iErr)
!!!         End Do
!!!      End If
!!!    
!!!      ! Write Node sets informations
!!!      If (Geom%Num_Node_Sets > 0) Then
!!!         Do iSet = 1, Geom%Num_node_sets
!!!            Call EXPNP(Geom%exoid, Geom%Node_Set(iSet)%ID, Geom%Node_Set(iSet)%Num_Nodes, Geom%Node_Set(iSet)%Num_Dist_Factors, iErr)
!!!            Call EXPNS(Geom%exoid, Geom%Node_Set(iSet)%ID, Geom%Node_Set(iSet)%Node_ID(:), iErr)
!!!         End Do
!!!      End If
!!!
!!!    ! Writes QA Records
!!!!    Geom%Num_QA = Geom%Num_QA+1
!!!!    Allocate (Tmp_QA(4, Geom%Num_QA))
!!!!    Tmp_QA(:,1:Geom%Num_QA-1) = Geom%QA_rec
!!!!    DeAllocate (Geom%QA_Rec)
!!!!    Allocate (Geom%QA_Rec(4, Geom%Num_QA))
!!!!    Geom%QA_rec = Tmp_QA
!!!!    DeAllocate (Tmp_QA)
!!!!    
!!!!    Geom%QA_Rec(1,Geom%Num_QA) = 'm_MEF_EXO'
!!!!    Geom%QA_Rec(2,Geom%Num_QA) = 'MEF90'
!!!!    Call Date_And_Time(date = Geom%QA_Rec(3,Geom%Num_QA))
!!!!    Call Date_And_Time(time = Geom%QA_Rec(4,Geom%Num_QA))
!!!!    Call EXPQA(Geom%exoid, Geom%num_QA, Geom%QA_Rec, iErr)
!!!    
!!!      Call EXCLOS(Geom%exoid, iErr)
!!!      Geom%exoid = 0
!!!   End Subroutine Write_MeshTopology_Info
!!!  
!!!   Subroutine Show_MeshTopology_Info(Geom, IO_Unit)
!!!      Type (MeshTopology_Info), Intent(IN)               :: Geom
!!!      Integer, Optional                              :: IO_Unit
!!!      
!!!      Integer                                        :: i
!!!    
!!!      If ( .NOT. Present(IO_UNit) ) Then
!!!         IO_Unit = 6 !!! STDOUT
!!!      End If  
!!!
!!!
!!!      Write(IO_Unit, 100) 
!!!      Write(IO_Unit, 101) Trim(Geom%filename)
!!!      Write(IO_Unit, 102) Trim(Geom%title)
!!!      Write(IO_Unit, 103) Geom%num_dim
!!!      Write(IO_Unit, 104) Geom%num_nodes
!!!      Write(IO_Unit, 105) Geom%num_elems
!!!      Write(IO_Unit, 106) Geom%Num_Elem_blks
!!!      Write(IO_Unit, 107) Geom%Num_Node_Sets
!!!      Write(IO_Unit, 108) Geom%Num_Side_Sets
!!!      
!!!      Write(IO_Unit, *)
!!!      Write(IO_Unit, 200)
!!!      Write(IO_Unit, 201) Geom%num_elem_blks
!!!      Do i = 1, Geom%Num_Elem_blks
!!!         Write(IO_Unit, 202) Geom%Elem_Blk(i)%ID, Geom%Elem_blk(i)%Type
!!!         Write(IO_Unit, 203) Geom%Elem_Blk(i)%ID, Geom%Elem_blk(i)%Num_Elems
!!!         Write(IO_Unit, 204) Geom%Elem_Blk(i)%ID, Geom%Elem_blk(i)%Num_Nodes_Per_Elem
!!!         Write(IO_Unit, 205) Geom%Elem_Blk(i)%ID, Geom%Elem_blk(i)%Num_Attr
!!!         Write(IO_Unit, 206, advance = 'no') Geom%Elem_Blk(i)%ID
!!!         Write(IO_Unit, *) Geom%Elem_blk(i)%Elem_ID
!!!      End Do
!!!      
!!!      
!!!      Write(IO_Unit, *)
!!!      Write(IO_Unit, 300)
!!!      Write(IO_Unit, 301) Geom%num_node_sets
!!!      Do i = 1, Geom%num_node_sets
!!!         Write(IO_Unit, 302) Geom%Node_Set(i)%ID, Geom%Node_Set(i)%Num_Nodes
!!!         Write(IO_Unit, 303, advance = 'no') Geom%Node_Set(i)%ID
!!!         Write(IO_Unit, *) Geom%Node_Set(i)%Node_ID
!!!         Write(IO_Unit, 304) Geom%Node_Set(i)%ID, Geom%Node_Set(i)%Num_Dist_Factors
!!!      End Do
!!!      
!!!      
!!!      Write(IO_Unit, *)
!!!      Write(IO_Unit, 400)
!!!      Write(IO_Unit, 401) Geom%num_side_sets
!!!      
!!!      Write(IO_Unit, *)
!!!      Write(IO_Unit, 500)
!!!      Write(IO_Unit, 501) Geom%num_QA
!!!      Do i = 1, Geom%num_QA
!!!         Write(IO_Unit, 502) i, Trim(Geom%QA_rec(1,i))
!!!         Write(IO_Unit, 503) i, Trim(Geom%QA_rec(2,i))
!!!         Write(IO_Unit, 504) i, Trim(Geom%QA_rec(3,i))
!!!         Write(IO_Unit, 505) i, Trim(Geom%QA_rec(4,i))
!!!      End Do
!!!
!!!100 Format('*** GLOBAL INFORMATIONS ***')
!!!101 Format('    Filename ======================== ', A)
!!!102 Format('    Title =========================== ', A)
!!!103 Format('    Number of dimensions ============ ', I1)
!!!104 Format('    Number of nodes ================= ', I6)
!!!105 Format('    Number of elements ============== ', I6)
!!!106 Format('    Number of elements blocks ======= ', I6)
!!!107 Format('    Number of node sets ============= ', I6)
!!!108 Format('    Number of side sets ============= ', I6)
!!!
!!!200 Format('*** ELEMENT BLOCKS ***')
!!!201 Format('    Number of blocks ================ ', I4)
!!!202 Format('    Block ', I3, ' Elements type ========= ', A)
!!!203 Format('    Block ', I3, ' Number of elements ==== ', I4)
!!!204 Format('    Block ', I3, ' Number of nodes per elt ', I4)
!!!205 Format('    Block ', I3, ' Number of attributes == ', I4)
!!!206 Format('    Block ', I3, ' IDs: ')
!!!
!!!300 Format('*** NODE SETS ***')
!!!301 Format('    Number of sets ================== ', I4)
!!!302 Format('    Set ', I3, ' Number of nodes ========= ', I4)
!!!303 Format('    Set ', I3, ' IDs: ')
!!!304 Format('    Set ', I3, ' Number of dist. factors = ', I4)
!!!
!!!400 Format('*** SIDE SETS ***')
!!!401 Format('    Number of side sets ============= ', I4)
!!!    
!!!500 Format('*** QA ***')
!!!501 Format('    Number of QA Records =========== ', I2)
!!!502 Format('    Rec ', I2, ' analysis code ============ ', A)
!!!503 Format('    Rec ', I2, ' analysis QA desc. ======== ', A)
!!!504 Format('    Rec ', I2, ' analysis date ============ ', A)
!!!505 Format('    Rec ', I2, ' analysis time ============ ', A)
!!!   End Subroutine Show_MeshTopology_Info




!!!
!!! RESULT FILES
!!!

!!! READ
!!! EXO stores values per coordinates
!!! V_PerNode = (/ (V_PerCoord(i:Size:Num_Nodes), i=1:NumDim /)
!!! V_PerCoord = (/ (V_PerNodes(i:Size:NumDim), i=1, Num_Nodes/)
!!! 

!!! In the sequel, we always assume that we are doing distributed I/O when EXO%comm == PETSC_COMM_SELF (i.e. 1 file per CPU) and sequential (i.e. IO operation on a single file on CPU 0) when EXO%comm == PETSC_COM_WORLD
   
!!!##   Subroutine Read_EXO_Result_Global(dEXO, dIdx, dTS, dRes, dNumRec)
!!!##      Type (EXO_Info), Intent(INOUT)                 :: dEXO
!!!##      Integer                                        :: dIdx
!!!##      Integer                                        :: dTS
!!!##      Real(Kind = Kr)                                :: dRes
!!!##      Integer                                        :: dNumRec
!!!##      
!!!##      Real(Kind = Kr)                                :: MyRes
!!!##      Integer                                        :: iErr
!!!##      Real(Kind = Kr)                                :: Vers
!!!##      Real(Kind = Kr), Dimension(:), Pointer         :: Tmp_Res
!!!##      Integer                                        :: Num_Vars, Num_TS
!!!##      Real(Kind = Kr)                                :: fDum
!!!##      Character                                      :: cDum
!!!##      
!!!##      If ( ((dEXO%comm == PETSC_COMM_WORLD) .AND. (MEF90_MyRank == 0)) .OR. (dEXO%comm == PETSC_COMM_SELF) ) Then
!!!##         dEXO%exoid = EXOPEN(dEXO%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers, ierr)
!!!##         ! Get the number of global variables stored in the database    
!!!##         Call EXGVP(dEXO%exoid, 'G', Num_Vars, iErr)
!!!##         Allocate(Tmp_Res(Num_Vars))
!!!##         
!!!##         ! Read All the global variables at time step TS and copy the right one
!!!##         ! into Res
!!!##         Call EXGGV(dEXO%exoid, dTS, Num_Vars, Tmp_Res, iErr)
!!!##         MyRes = Tmp_Res(dIdx)
!!!##         DeAllocate (Tmp_Res)
!!!##         Call EXCLOS(dEXO%exoid, iErr)
!!!##         dEXO%exoid = 0
!!!##         End If
!!!##      End If
!!!##            
!!!##      !!! Broacast if dEXO%comm == PETSC_COMM_WORLD, AllReduce if dEXO%comm == PETSC_COMM_SELF
!!!##      If (dEXO%comm == PETSC_COMM_WORLD) Then
!!!##         Call MPI_BCast(MyRes, 1, MPIU_SCALAR, 0, dEXO%comm, iErr)
!!!##         dRes = MyRes
!!!##      Else
!!!##         Call PetscGlobalSum(MyRes, dRes, dEXO%comm, iErr)
!!!##      End If
!!!##   End Subroutine Read_EXO_Result_Global
!!!##
!!!##   Subroutine Read_EXO_Result_Section_Nodes(dExo, dMeshTopology, dIdx, dTS, dRes)
!!!##      Type (EXO_Info), Intent(INOUT)                 :: dEXO
!!!##      Type (MeshTopology_Info)                       :: dMeshTopology
!!!##      Integer                                        :: dIdx
!!!##      Integer                                        :: dTS
!!!##      SectionReal                                    :: dRes
!!!##      
!!!##      Vec                                            :: Res_Vec
!!!##      Real(Kind = Kr), Dimension(:), Pointer         :: Res_Ptr
!!!##      Integer                                        :: iErr
!!!##      Real(Kind = Kr)                                :: Vers
!!!##      
!!!##      Geom%exoid = EXOPEN(dEXO%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers, ierr)
!!!##      
!!!##      !!! We Assume that the section is initialized and has the proper size
!!!##      
!!!##      Call EXGNV(Geom%exoid, TS, Idx, Geom%Num_Nodes, Res, iErr)
!!!##      
!!!##      Call EXCLOS(Geom%exoid, iErr)
!!!##      Geom%exoid = 0
!!!##   End Subroutine Read_EXO_Result_Section_Nodes
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
!!!##
!!!##!!! WRITE
!!!##  Subroutine Write_EXO_Result_Global(Geom, Idx, TS, Res)
!!!##    Type (MeshTopology_Info), Intent(INOUT)            :: Geom
!!!##    Integer                                        :: Idx
!!!##    Integer                                        :: TS
!!!##    Real(Kind = Kr)                                :: Res
!!!##
!!!##    Integer                                        :: iErr
!!!##    Real(Kind = Kr)                                :: Vers
!!!##    Real(Kind = Kr), Dimension(:), Pointer         :: Tmp_Res
!!!##    Integer                                        :: Num_Vars, Num_TS
!!!##    Real(Kind = Kr)                                :: fDum
!!!##    Character                                      :: cDum
!!!##
!!!##    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
!!!##         & ierr)
!!!##    ! Get the number of global variables stored in the database    
!!!##    Call EXGVP(Geom%exoid, 'G', Num_Vars, iErr)
!!!##    Allocate(Tmp_Res(Num_Vars))
!!!##    Tmp_Res = 0.0_Kr
!!!##
!!!##    ! Read All the global variables at time step TS into Tmp_Res
!!!##    ! Modify Tmp_Res(Idx) and write everything back...
!!!##    Call EXGGV(Geom%exoid, TS, Num_Vars, Tmp_Res, iErr)
!!!##    Tmp_Res(Idx) = Res
!!!##
!!!##    Call EXPGV(Geom%exoid, TS, Num_Vars, Tmp_Res, iErr)
!!!##    DeAllocate (Tmp_Res)
!!!##    Call EXCLOS(Geom%exoid, iErr)
!!!##    Geom%exoid = 0
!!!##  End Subroutine Write_EXO_Result_Global
!!!##
!!!##
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

