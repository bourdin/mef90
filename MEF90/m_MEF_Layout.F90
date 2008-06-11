Module m_MEF_Layout

   Use m_AlgebLin
   Use m_Constantes
   Use m_MEF_Types
   Use m_MEF_EXO
   IMPLICIT NONE
   Private

#include "include/finclude/petsc.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petscao.h"
  include "exodusII.inc"


   Public :: Layout_Info
   Public :: Init_Layout_TRI3
   
   Interface Init_Layout_TRI3
      Module Procedure Init_Layout_TRI3_Scal
   End Interface
   
 Contains
 
   Subroutine Init_Layout_TRI3_Scal(Geom, MyGeom, Layout, MyNode_db, MyElem_db)
      Type (EXO_Geom_Info), Intent(INOUT)           :: Geom
      Type (EXO_Geom_Info), Intent(OUT)             :: MyGeom
      Type (Layout_Info), Intent(OUT)               :: Layout
      Type (Node2D), Dimension(:), Pointer          :: MyNode_db
      TYpe (Element2D_Scal), Dimension(:), Pointer  :: MyElem_db
                                                    
      Integer, Dimension(:,:), Pointer              :: MyConnect
                                                    
      Integer, Dimension(:), Pointer                :: METIS_connect
      Integer, Dimension(:), Pointer                :: METIS_epart
      Integer, Dimension(:), Pointer                :: METIS_npart
      Integer                                       :: METIS_edgecut
      Integer                                       :: METIS_etype = 1
      Integer                                       :: METIS_numflag = 1
                                                    
      Integer                                       :: NumProcs, MyRank
                                                    
                                                    
      Integer                                       :: Offset = 0
      Integer                                       :: iBlk, iSet, iErr
      Integer                                       :: i, j, iE, iS
      Integer, Dimension(:), Pointer                :: MyNodes, MyElems
      Logical, Dimension(:), Pointer                :: IsGhost
      Real(Kind = Kr)                               :: Vers
      Vec                                           :: V_Dist, V_IO
      Integer                                       :: MyIdxMinN, MyIdxMaxN
      Integer                                       :: MyIdxMinE, MyIdxMaxE
      Integer, Dimension(:), Pointer                :: Tmp_Mapping_Array, Tmp_Connect, Tmp_Idx
      Real(Kind = Kr), Dimension(:), Pointer        :: Tmp_Coord_X, Tmp_Coord_Y, Tmp_Coord_Z
      AO                                            :: AO_N, AO_E

      !!! 1. Some initialization
      !!!    We Assume that Read_EXO_Geom_Info (or equivalent) has been called on each cpu 
      Call MPI_Comm_Size(Geom%comm, NumProcs, iErr)
      Call MPI_Comm_Rank(Geom%comm, MyRank, iErr)

      !!! IO node reads the mesh and gets the connectivity table, then BCasts it
      Allocate( METIS_connect(Geom%num_elems * 3) )
      
      If (MyRank == 0) Then
         Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers, ierr)

         Offset = 0
         Do iBlk = 1, Geom%num_elem_blks
            Call EXGELC(Geom%exoid, Geom%Elem_Blk(iBlk)%ID, METIS_connect(Offset+1:Offset+3*Geom%Elem_Blk(iBlk)%Num_Elems), iErr)
            Offset = Offset + 3*Geom%Elem_Blk(iBlk)%Num_Elems
         End Do
         Call EXCLOS(Geom%exoid, iErr)
         Geom%exoid = 0
      End If
      Call MPI_BCast(METIS_connect, Geom%num_elems * 3, MPI_INTEGER, 0, Geom%Comm, iErr)
      
      !!! Same for  coordinates
      Allocate (Tmp_Coord_X(Geom%Num_Nodes))
      Allocate (Tmp_Coord_Y(Geom%Num_Nodes))
      Allocate (Tmp_Coord_Z(Geom%Num_Nodes))
      If (MyRank == 0) Then
         Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers, ierr)

         Call EXGCOR(Geom%exoid, Tmp_Coord_X, Tmp_Coord_Y, Tmp_Coord_Z, iErr)
         Call EXCLOS(Geom%exoid, iErr)
         Geom%exoid = 0
      End If
      Call MPI_BCast(Tmp_Coord_X, Geom%num_nodes, MPIU_SCALAR, 0, Geom%Comm, iErr)
      Call MPI_BCast(Tmp_Coord_Y, Geom%num_nodes, MPIU_SCALAR, 0, Geom%Comm, iErr)

      Allocate(MyGeom%elem_blk(MyGeom%num_elem_blks))
      Allocate(MyGeom%Node_Set(MyGeom%Num_Node_Sets))

      ! Broadcasting Elem Block
      Do iBlk = 1, MyGeom%num_elem_blks
         If (MyRank == 0) Then
            MyGeom%Elem_blk(iBlk)%ID                 = Geom%Elem_blk(iBlk)%ID
            MyGeom%Elem_blk(iBlk)%Type               = Geom%Elem_blk(iBlk)%Type
            MyGeom%Elem_blk(iBlk)%Num_Nodes_Per_Elem = Geom%Elem_blk(iBlk)%Num_Nodes_Per_Elem
            MyGeom%Elem_blk(iBlk)%Num_Attr           = Geom%Elem_blk(iBlk)%Num_Attr
         End If
         Call MPI_BCAST(MyGeom%Elem_Blk(iBlk)%ID, 1, MPI_INTEGER, 0, Geom%Comm, iErr)
         Call MPI_BCAST(MyGeom%Elem_Blk(iBlk)%Type, MXLNLN, MPI_CHARACTER, 0, Geom%Comm, iErr)
         Call MPI_BCAST(MyGeom%Elem_Blk(iBlk)%Num_Elems, 1, MPI_INTEGER, 0, Geom%Comm, iErr)
         Call MPI_BCAST(MyGeom%Elem_Blk(iBlk)%Num_Nodes_Per_Elem, 1, MPI_INTEGER, 0, Geom%Comm, iErr)
         Call MPI_BCAST(MyGeom%Elem_Blk(iBlk)%Num_Attr, 1, MPI_INTEGER, 0, Geom%Comm, iErr)         
      End Do
         
      ! Broadcasting Node Sets
      Do iSet = 1, MyGeom%num_node_sets
         If (MyRank == 0) Then      
            MyGeom%Node_Set(iSet)%ID               = Geom%Node_Set(iSet)%ID
         End If
         Call MPI_BCAST(MyGeom%Node_Set(iSet)%ID, 1, MPI_INTEGER, 0, Geom%Comm, iErr)
      End Do
      
      Call MPI_BCAST(MyGeom%Numbering, 1, MPI_INTEGER, 0, Geom%Comm, iErr)
      MyGeom%Comm = PETSC_COMM_SELF
      !!! Geom%Numbering is obsolete

      !!! Call METIS and start filling in MyGeom and Layout
      Allocate(METIS_epart(Geom%num_elems))
      Allocate(METIS_npart(Geom%num_nodes))
      If (NumProcs > 1) Then
         Call METIS_PartMeshDual(Geom%num_elems, Geom%num_nodes, METIS_connect, METIS_etype, METIS_numflag, NumProcs, METIS_edgecut, METIS_epart, METIS_npart)
      Else
         !!! Metis crashes when called with npart=1
         METIS_npart = 1
         METIS_epart = 1
      End If

!!!   REMOVE
      
      Layout%num_local_dof = count(METIS_npart == MEF90_MyRank+1)
      Layout%num_local_elems = count(METIS_epart == MEF90_MyRank+1)
      
      !!! Find local and ghost nodes and elements
      Layout%num_ghost_elems = 0
      Allocate(MyNodes(Layout%num_local_dof))
      Offset = 1
      Do i = 1, Geom%num_nodes
         If(METIS_npart(i) == MyRank+1) Then
            MyNodes(Offset) = i
            Offset = Offset + 1
         End If
      End Do

      Allocate(MyElems(Layout%num_local_elems))
      Offset = 1
      Do i = 1, Geom%num_elems
         If(METIS_epart(i) == MyRank+1) Then
            MyElems(Offset) = i
            Offset = Offset + 1
         End If
      End Do
      
      Allocate(IsGhost(Geom%num_nodes))
      IsGhost = .FALSE.
      Do i = 1, Geom%num_elems
         If (METIS_epart(i) == MyRank+1) Then
            Do j = 1, 3
               If (METIS_npart( METIS_connect( 3*(i-1)+j) ) /= MyRank+1) Then
                  IsGhost( METIS_connect( 3*(i-1)+j) ) = .TRUE.
               End If
            End Do
         End If
      End Do
      Layout%num_ghost_dof = count(IsGhost)
      Allocate(Layout%ghost_dof(Layout%num_ghost_dof)) 
      Offset = 1
      Do i = 1, Geom%num_nodes
         If (IsGhost(i)) Then
            Layout%ghost_dof(Offset) = i     
            Offset = Offset + 1
         End If
      End Do 
      !!! At this point, the ghost number still refer to the mesh numbering. They will have to be renumbered as soon as AO_N is initialized
      
      MyGeom%num_nodes = Layout%num_local_dof + Layout%num_ghost_dof
      MyGeom%num_elems = Layout%num_local_elems + Layout%num_ghost_elems
      !!!! Write(100+MyRank,*) "MyNodes", Layout%num_local_dof
      !!!! Write(100+MyRank,*) MyNodes
      !!!! Write(100+MyRank,*) "MyElems", Layout%num_local_elems
      !!!! Write(100+MyRank,*) MyElems
      !!!! Write(100+MyRank,*) "Layout%ghost_dofs", Layout%num_ghost_dof
      !!!! Write(100+MyRank,*) Layout%ghost_dof
      

!!!   KEEP THIS PART WHICH REQUIRES ONLY MyNodes, MyElems, MyGhosts AND  NumDoFPerNode
      !!! IS and Scatters for Sequential IO and IdxMin/Max
      Call ISCreateGeneral(Geom%Comm, Layout%num_local_dof, MyNodes-1, Layout%IS_N, iErr)
      Call AOCreateBasicIS(Layout%IS_N, PETSC_NULL, AO_N, iErr)

      Call VecCreateMPI(Geom%Comm, Layout%num_local_dof, Geom%num_nodes, V_Dist, iErr)
      Call VecGetOwnershipRange(V_Dist, MyIdxMinN, MyIdxMaxN, iErr)
      If (MyRank == 0) Then
         Call VecCreateMPI(Geom%Comm, Geom%num_nodes, Geom%num_nodes, V_IO, iErr)
      Else
         Call VecCreateMPI(Geom%Comm, 0, Geom%num_nodes, V_IO, iErr)
      End If   
      Call VecScatterCreate(V_Dist, PETSC_NULL, V_IO, Layout%IS_N, Layout%ToIOSeq_N, iErr)
      Call VecDestroy(V_Dist, iErr)
      Call VecDestroy(V_IO, iErr)

      Call ISCreateGeneral(Geom%Comm, Layout%num_local_elems, MyElems-1, Layout%IS_E, iErr)
      Call AOCreateBasicIS(Layout%IS_E, PETSC_NULL, AO_E, iErr)

      Call VecCreateMPI(Geom%Comm, Layout%num_local_elems, Geom%num_elems, V_Dist, iErr)
      Call VecGetOwnershipRange(V_Dist, MyIdxMinE, MyIdxMaxE, iErr)
      If (MyRank == 0) Then
         Call VecCreateMPI(Geom%Comm, Geom%num_elems, Geom%num_elems, V_IO, iErr)
      Else
         Call VecCreateMPI(Geom%Comm, 0, Geom%num_elems, V_IO, iErr)
      End If   
      Call VecScatterCreate(V_Dist, PETSC_NULL, V_IO, Layout%IS_E, Layout%ToIOSeq_E, iErr)
      Call VecDestroy(V_Dist, iErr)
      Call VecDestroy(V_IO, iErr)
      
      !!! IS and Scatter for Dist IO
      Call VecCreateSeq(PETSC_COMM_SELF, Layout%num_local_dof + Layout%num_ghost_dof, V_IO, iErr)
      Allocate(Tmp_Idx((Layout%num_local_dof + Layout%num_ghost_dof)))
      DeAllocate(Tmp_Idx)
      Call VecDestroy(V_IO, iErr)
      
      
!!!   I NEED TO DECIDE IF  Layout%num_local_dof AND Layout%num_ghost_dof SHOULD REFER TO THE NUMBER OF VERTICES OR DOF
!!!   IN ORDER TO BE CONSISTENT IT SHOULD BE DOF

      !!!! Call PetscPrintf(PETSC_COMM_WORLD, "IS_N\n"c, iErr)
      !!!! Call ISView(Layout%IS_N, PETSC_VIEWER_STDOUT_WORLD, iErr)
      !!!! Call PetscPrintf(PETSC_COMM_WORLD, "IS_E\n"c, iErr)
      !!!! Call ISView(Layout%IS_E, PETSC_VIEWER_STDOUT_WORLD, iErr)
      !!!! 
      !!!! Call PetscPrintf(PETSC_COMM_WORLD, "AO_N\n"c, iErr)
      !!!! Call AOView(AO_N, PETSC_VIEWER_STDOUT_WORLD, iErr)
      !!!! Call PetscPrintf(PETSC_COMM_WORLD, "AO_E\n"c, iErr)
      !!!! Call AOView(AO_E, PETSC_VIEWER_STDOUT_WORLD, iErr)
      
      !!! LocalToGlobalMapping
      Allocate(Tmp_Mapping_Array(Layout%num_local_elems))
      Tmp_Mapping_Array = (/ (i+MyIdxMinE, i=0, Layout%num_local_elems - 1) /)
      Call ISLocalToGlobalMappingCreate(PETSC_COMM_SELF, Layout%num_local_elems, Tmp_Mapping_Array, Layout%Mapping_E, iErr)
      DeAllocate(Tmp_Mapping_Array)

      Allocate(Tmp_Mapping_Array(MyGeom%num_nodes))
      Tmp_Mapping_Array = 0
      Tmp_Mapping_Array(Layout%num_local_dof+1:MyGeom%num_nodes) = Layout%ghost_dof-1
      Call AOApplicationToPetsc(AO_N, MyGeom%num_nodes, Tmp_Mapping_Array, iErr)
      Tmp_Mapping_Array(1:Layout%num_local_dof) = (/ (i+MyIdxMinN, i=0, Layout%num_local_dof - 1) /)
      Call ISLocalToGlobalMappingCreate(PETSC_COMM_SELF, MyGeom%num_nodes, Tmp_Mapping_Array, Layout%Mapping_N, iErr)
      DeAllocate(Tmp_Mapping_Array)      
      
      !!!! Call PetscPrintf(PETSC_COMM_WORLD, "Layout_E\n"c, iErr)
      !!!! Call ISLocalToGlobalMappingView(Layout%Mapping_E, PETSC_VIEWER_STDOUT_WORLD, iErr)
      !!!! Call PetscPrintf(PETSC_COMM_WORLD, "Layout_N\n"c, iErr)
      !!!! Call ISLocalToGlobalMappingView(Layout%Mapping_N, PETSC_VIEWER_STDOUT_WORLD, iErr)
      
      !!! Renumbering the Element Blocks
      Do iBlk = 1, MyGeom%num_elem_blks
         Allocate(Tmp_Mapping_Array(Geom%elem_blk(iBlk)%num_elems))

         Tmp_Mapping_Array = Geom%elem_blk(iBlk)%Elem_ID-1
         Call AOApplicationToPetsc(AO_E, Geom%elem_blk(iBlk)%num_elems, Tmp_Mapping_Array, iErr)

         Allocate(MyGeom%elem_blk(iBlk)%Elem_ID(Geom%elem_blk(iBlk)%Num_Elems))

         Call ISGlobalToLocalMappingApply(Layout%Mapping_E, IS_GTOLM_DROP, Geom%elem_blk(iBlk)%Num_Elems, Tmp_Mapping_Array, MyGeom%elem_blk(iBlk)%Num_Elems, MyGeom%elem_blk(iBlk)%Elem_ID, iErr)
         DeAllocate(MyGeom%elem_blk(iBlk)%Elem_ID)

         Allocate(MyGeom%elem_blk(iBlk)%Elem_ID(MyGeom%elem_blk(iBlk)%Num_Elems))
         Call ISGlobalToLocalMappingApply(Layout%Mapping_E, IS_GTOLM_DROP, Geom%elem_blk(iBlk)%Num_Elems, Tmp_Mapping_Array, MyGeom%elem_blk(iBlk)%Num_Elems, MyGeom%elem_blk(iBlk)%Elem_ID, iErr)
         MyGeom%elem_blk(iBlk)%Elem_ID = MyGeom%elem_blk(iBlk)%Elem_ID+1

         !!!! Write(100+MyRank, *) "Geom%Blk", iBlk, Geom%elem_blk(iBlk)%Num_Elems
         !!!! Write(100+MyRank, *) Geom%elem_blk(iBlk)%ELem_ID
         !!!! Write(100+MyRank, *) "MyGeom%Blk", iBlk, MyGeom%elem_blk(iBlk)%Num_Elems
         !!!! Write(100+MyRank, *) MyGeom%elem_blk(iBlk)%ELem_ID
         DeAllocate(Tmp_Mapping_Array)
      End Do

      !!! Renumbering the node sets
      Do iSet = 1, MyGeom%num_node_sets
         MyGeom%node_set(iSet)%ID               = Geom%Node_Set(iSet)%ID
         
         Allocate(Tmp_Mapping_Array(Geom%node_set(iSet)%Num_Nodes))
         Tmp_Mapping_Array = Geom%node_set(iSet)%Node_ID-1
         Call AOApplicationToPetsc(AO_N, Geom%node_set(iSet)%num_nodes, Tmp_Mapping_Array, iErr)

         Allocate(MyGeom%node_set(iSet)%Node_ID(Geom%node_set(iSet)%Num_Nodes))
         Call ISGlobalToLocalMappingApply(Layout%Mapping_N, IS_GTOLM_DROP, Geom%node_set(iSet)%Num_Nodes, Tmp_Mapping_Array, MyGeom%node_set(iSet)%Num_Nodes, MyGeom%node_set(iSet)%Node_ID, iErr)
         DeAllocate(MyGeom%node_set(iSet)%Node_ID)

         Allocate(MyGeom%node_set(iSet)%Node_ID(MyGeom%node_set(iSet)%Num_Nodes))
         Call ISGlobalToLocalMappingApply(Layout%Mapping_N, IS_GTOLM_DROP, Geom%node_set(iSet)%Num_Nodes, Tmp_Mapping_Array, MyGeom%node_set(iSet)%Num_Nodes, MyGeom%node_set(iSet)%Node_ID, iErr)
         MyGeom%node_set(iSet)%Node_ID = MyGeom%node_set(iSet)%Node_ID+1
         MyGeom%Node_Set(iSet)%Num_Dist_Factors = MyGeom%node_set(iSet)%num_nodes
         
         !!!! Write(100+MyRank, *) "Geom%NodeSet", iSet, Geom%node_set(iSet)%Num_nodes
         !!!! Write(100+MyRank, *) Geom%node_set(iSet)%Node_ID
         !!!! Write(100+MyRank, *) "MyGeom%NodeSet", iSet, MyGeom%node_set(iSet)%Num_nodes
         !!!! Write(100+MyRank, *) MyGeom%node_set(iSet)%Node_ID
         DeAllocate(Tmp_Mapping_Array)
      End Do
      
!!! MOVE THIS INTO A SUBROUTINE NAME INIT_CONNECT{2D,2D_Scal,2D_Elast, 3D, 3D_Scal, 3D_Elast)
!!! RETURN THE PREALLOCATED Node and Elem OBJECT      
!!! DONE!

      !!! Construction of the local Nodes and Elem data structures
      Allocate (MyElem_db(MyGeom%Num_Elems))
      MyElem_db(:)%NB_DoF = 3
      
      Do iBlk = 1, MyGeom%num_elem_blks
         Do i = 1, MyGeom%elem_blk(iBlk)%num_Elems
            iE = MyGeom%elem_blk(iBlk)%Elem_ID(i)
            Allocate(MyElem_db(iE)%ID_DoF(MyElem_db(iE)%Nb_DoF))
            MyElem_db(iE)%ID_DoF = METIS_connect((MyElems(iE)-1) * 3 +1: MyElems(iE) * 3) -1
            Call AOApplicationToPetsc(AO_N, 3, MyElem_db(iE)%ID_DoF, iErr)
            Call ISGlobalToLocalMappingApply(Layout%Mapping_N, IS_GTOLM_MASK, 3, MyElem_db(iE)%ID_DoF, PETSC_NULL_INTEGER, MyElem_db(iE)%ID_DoF, iErr)
            MyElem_db(iE)%ID_DoF = MyElem_db(iE)%ID_DoF + 1
            MyElem_db(iE)%ID_EL  = iBlk
         End Do
      End Do

      !!! Construction of the local coordinates
      !!!??? What do I do with Node_db%ID ?
      Allocate(MyNode_db(MyGeom%num_nodes))
      Do i = 1, Layout%num_local_dof
         MyNode_db(i)%Coord%X = Tmp_Coord_X(MyNodes(i))
         MyNode_db(i)%Coord%Y = Tmp_Coord_Y(MyNodes(i))
      EndDo
      Do i = 1, Layout%num_ghost_dof
         MyNode_db(Layout%num_local_dof + i)%Coord%X = Tmp_Coord_X(Layout%ghost_dof(i))
         MyNode_db(Layout%num_local_dof + i)%Coord%Y = Tmp_Coord_Y(Layout%ghost_dof(i))
      EndDo

      !!! Renumbering the ghost dof
      !!! DANGER: THIS CAN NOT BE DONE BEFORE THE RECONTRUCTION OF THE CONNECTIVITY TABLE
      !!! I KNOW THIS IS UGLY
      Allocate(Tmp_Mapping_Array(Layout%num_ghost_dof))
      Tmp_Mapping_Array = Layout%ghost_dof - 1
      Call AOApplicationToPetsc(AO_N, Layout%num_ghost_dof, Tmp_Mapping_Array, iErr)
      Layout%ghost_dof = Tmp_Mapping_array+1
      DeAllocate(Tmp_Mapping_Array)



      !!!! Write(100+MyRank, *) "Coords"
      !!!! Do i = 1, Geom%Num_Nodes
      !!!!    Write(100+MyRank, *) i, Tmp_Coord_X(i), Tmp_Coord_Y(i)
      !!!! End Do
      !!!! Write(100+MyRank, *) "MyCoords"
      !!!! Do i = 1, MyGeom%Num_Nodes
      !!!!    Write(100+MyRank, *) i, MyCoords(i)%X, MyCoords(i)%Y
      !!!! End Do

      Call AODestroy(AO_N, iErr)
      Call AODestroy(AO_E, iErr)
      DeAllocate(METIS_npart)
      DeAllocate(METIS_epart)
      DeAllocate(METIS_connect)
      DeAllocate(MyNodes)
      DeAllocate(MyElems)
      DeAllocate(IsGhost)
      DeAllocate(Tmp_Coord_X)
      DeAllocate(Tmp_Coord_Y)
      DeAllocate(Tmp_Coord_Z)
   End Subroutine Init_Layout_TRI3_Scal
   
   
End Module m_MEF_Layout