Module m_MEF_EXO
!!! ATTENTION:
!!! DoF ordering has changed. the order is now 1%X, 2%X, ..., 1%Y, 2%Y...
!!! Instead of 1%X, 1%Y, 2%X, 2%Y...
!!! I/O routines have to be changed accordingly...
!!! This is consistent with ensight ordering scheme, and easier to implement

!!! NOT ANYMORE...
!!! Uses Geom%Numbering to know whci scheme is used.
!!! Numbering_PerNodes -> 1%X, 1%Y, 1%Z, 2%X, 2%Y, 2%Z, ...
!!! Numbering_PerCoord -> 1%X, 2%X, ..., 1%Y, 2%Y, ..., 1%Z, 2%Z
  Use m_MEF_Types
  IMPLICIT NONE
  

  Integer, Parameter                                :: exo_cpu_ws = 8
  Integer, Parameter                                :: exo_io_ws = 8

  Interface Read_EXO_Node_Coord
     Module Procedure Read_EXO_Node_Coord_1D, Read_EXO_Node_Coord_2D,         &
          & Read_EXO_Node_Coord_3D
  End Interface
  
  Interface Write_EXO_Node_Coord
     Module Procedure Write_EXO_Node_Coord_1D, Write_EXO_Node_Coord_2D,       &
          & Write_EXO_Node_Coord_3D
  End Interface

  Interface Read_EXO_Connect
     Module Procedure  Read_EXO_Connect_2D,  Read_EXO_Connect_2D_Scal,        &
          & Read_EXO_Connect_2D_Elast, Read_EXO_Connect_3D,                   &
          & Read_EXO_Connect_3D_Scal, Read_EXO_Connect_3D_Elast
  End Interface

  Interface Write_EXO_Connect
     Module Procedure  Write_EXO_Connect_2D,  Write_EXO_Connect_2D_Scal,      &
          & Write_EXO_Connect_2D_Elast, Write_EXO_Connect_3D,                 &
          & Write_EXO_Connect_3D_Scal, Write_EXO_Connect_3D_Elast
  End Interface

  Interface Read_EXO_Result_Nodes
     Module Procedure Read_EXO_Result_Scal_Nodes,                             &
          & Read_EXO_Result_Ptr_Nodes, Read_EXO_Result_Vect2D_Nodes,          &
          & Read_EXO_Result_Mat2D_Nodes, Read_EXO_Result_MatS2D_Nodes,        &
          & Read_EXO_Result_Vect3D_Nodes, Read_EXO_Result_Mat3D_Nodes,        &
          & Read_EXO_Result_MatS3D_Nodes
  End Interface

  Interface Write_EXO_Result_Nodes
     Module Procedure                                                         &
          & Write_EXO_Result_Ptr_Nodes, Write_EXO_Result_Vect2D_Nodes,        &
          & Write_EXO_Result_Mat2D_Nodes, Write_EXO_Result_MatS2D_Nodes,      &
          & Write_EXO_Result_Vect3D_Nodes, Write_EXO_Result_Mat3D_Nodes,      &
          & Write_EXO_Result_MatS3D_Nodes
  End Interface

  Interface Read_EXO_Result_Elems
     Module Procedure Read_EXO_Result_Scal_Elems,                             &
          & Read_EXO_Result_Vect2D_Elems,          &
          & Read_EXO_Result_Mat2D_Elems, Read_EXO_Result_MatS2D_Elems,        &
          & Read_EXO_Result_Vect3D_Elems, Read_EXO_Result_Mat3D_Elems,        &
          & Read_EXO_Result_MatS3D_Elems
  End Interface

  Interface Write_EXO_Result_Elems
     Module Procedure Write_EXO_Result_Scal_Elems,                            &
          & Write_EXO_Result_Vect2D_Elems,        &
          & Write_EXO_Result_Mat2D_Elems, Write_EXO_Result_MatS2D_Elems,      &
          & Write_EXO_Result_Vect3D_Elems, Write_EXO_Result_Mat3D_Elems,      &
          & Write_EXO_Result_MatS3D_Elems
  End Interface

Contains
  Subroutine Read_EXO_Geom_Info(Geom)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom

    Integer                                        :: vers
    Integer                                        :: iErr
    Integer                                        :: iDummy
    Real(Kind = Kr)                                :: rDummy
    Character                                      :: cDummy

    Integer                                        :: iBlk
    Integer                                        :: iSet
    Integer(Kind = Ki)                             :: Offset,i

    ! Open File
    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    ! Read Global Geometric Parameters
    Call EXGINI(Geom%exoid, Geom%Title, Geom%Num_Dim, Geom%Num_Nodes,         &
         & Geom%Num_Elems, Geom%Num_Elem_Blks, Geom%Num_Node_Sets,            &
         & Geom%Num_Side_Sets, iErr)
    If (Geom%Num_Side_Sets > 0) Then
       Write(*,*) 'WARNING: Side_sets non implemented, continuing anyway'
       Write(*,*) '         Got Geom%Num_Side_Sets = ', Geom%Num_Side_Sets
    End If
    
    ! Read Elem blocks informations
    Allocate(Geom%Elem_blk(Geom%Num_Elem_blks))
    Offset = 0
    If (Geom%Num_Elem_blks > 0) Then
       Call EXGEBI(Geom%exoid, Geom%Elem_Blk(:)%ID, iErr)
       Do iBlk = 1, Geom%Num_Elem_Blks
          Call EXGELB(Geom%exoid, Geom%elem_blk(iBlk)%ID,                     &
               & Geom%elem_blk(iBlk)%Type, Geom%elem_blk(iBlk)%Num_Elems,     &
               & Geom%elem_blk(iBlk)%Num_Nodes_Per_Elem,                      &
               & Geom%elem_blk(iBlk)%Num_Attr, iErr)
          Allocate(Geom%Elem_blk(iBlk)%Elem_ID(Geom%elem_blk(iBlk)%Num_Elems))
          Geom%Elem_Blk(iBlk)%Elem_ID = (/ (Offset+i,                         &
               &                           i=1,Geom%elem_blk(iBlk)%Num_Elems)/)
!          Do i = 1+Offset, Geom%Elem_Blk(iBlk)%Num_Elems+Offset
!             Geom%Elem_Blk(iBlk)%Elem_ID(i) = i
!          End Do
          Offset = Offset + Geom%elem_blk(iBlk)%Num_Elems
       End Do
    End If

    ! Read Node sets informations
    Allocate (Geom%Node_Set(Geom%Num_Node_Sets))
    If (Geom%Num_Node_Sets > 0) Then
       Call EXGNSI(Geom%exoid, Geom%Node_Set(:)%ID, iErr)
       Do iSet = 1, Geom%Num_node_sets
          Call EXGNP(Geom%exoid, Geom%Node_Set(iSet)%ID,                      &
               & Geom%Node_Set(iSet)%Num_Nodes,                               &
               & Geom%Node_Set(iSet)%Num_Dist_Factors, iErr)
          Allocate(Geom%Node_Set(iSet)%Node_ID(Geom%Node_Set(iSet)%Num_Nodes))
          Allocate(Geom%Node_Set(iSet)%Dist_Factor(                           &
               & Geom%Node_Set(iSet)%Num_Dist_Factors))
          Call EXGNS(Geom%exoid, Geom%Node_Set(iSet)%ID,                      &
               & Geom%Node_Set(iSet)%Node_ID(:), iErr)
       End Do
    End If

    ! Read QA Records
    Call EXINQ(Geom%exoid, EXQA, Geom%num_QA, rDummy, cDummy, iErr)
    Allocate (Geom%QA_Rec(4, Geom%Num_QA))
    Call EXGQA(Geom%exoid, Geom%QA_Rec, iErr)
    

    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Read_EXO_Geom_Info
  
  Subroutine Write_EXO_Geom_Info(Geom)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    
    Integer                                        :: vers
    Integer                                        :: iErr
    Integer                                        :: iDummy
    Real(Kind = Kr)                                :: rDummy
    Character                                      :: cDummy
    
    Integer                                        :: iBlk
    Integer                                        :: iSet
    Character(len=MXSTLN), Dimension(:,:), Pointer :: Tmp_QA
    Character(len=MXSTLN), Dimension(3)            :: Coord_Names

    Coord_Names(1) = 'X'
    Coord_Names(2) = 'Y'
    Coord_Names(3) = 'Z'
    
    
    ! Open File
    Geom%exoid = EXCRE (Geom%filename, EXCLOB, exo_cpu_ws, exo_io_ws, iErr)
    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)
    
    ! Writes Global Geometric Parameters
    Call EXPINI(Geom%exoid, Geom%Title, Geom%Num_Dim, Geom%Num_Nodes,         &
         & Geom%Num_Elems, Geom%Num_Elem_Blks, Geom%Num_Node_Sets,            &
         & Geom%Num_Side_Sets, iErr)
    Call EXPCON(Geom%exoid, Coord_Names, iErr)
    
    If (Geom%Num_Side_Sets > 0) Then
       Write(*,*) 'WARNING: Side_sets non implemented, continuing anyway'
       Write(*,*) '         Got Geom%Num_Side_Sets = ', Geom%Num_Side_Sets
    End If
    
    ! Write Elem blocks informations
    If (Geom%Num_Elem_blks > 0) Then
       Do iBlk = 1, Geom%Num_Elem_Blks
          Call EXPELB(Geom%exoid, Geom%elem_blk(iBlk)%ID,                     &
               & Geom%elem_blk(iBlk)%Type, Geom%elem_blk(iBlk)%Num_Elems,     &
               & Geom%elem_blk(iBlk)%Num_Nodes_Per_Elem,                      &
               & Geom%elem_blk(iBlk)%Num_Attr, iErr)
       End Do
    End If
    
    ! Write Node sets informations
    If (Geom%Num_Node_Sets > 0) Then
       Do iSet = 1, Geom%Num_node_sets
          Call EXPNP(Geom%exoid, Geom%Node_Set(iSet)%ID,                      &
               & Geom%Node_Set(iSet)%Num_Nodes,                               &
               & Geom%Node_Set(iSet)%Num_Dist_Factors, iErr)
          Call EXPNS(Geom%exoid, Geom%Node_Set(iSet)%ID,                      &
               & Geom%Node_Set(iSet)%Node_ID(:), iErr)
       End Do
    End If
    
    ! Writes QA Records
    Geom%Num_QA = Geom%Num_QA+1
    Allocate (Tmp_QA(4, Geom%Num_QA))
    Tmp_QA(:,1:Geom%Num_QA-1) = Geom%QA_rec
    DeAllocate (Geom%QA_Rec)
    Allocate (Geom%QA_Rec(4, Geom%Num_QA))
    Geom%QA_rec = Tmp_QA
    DeAllocate (Tmp_QA)
    
    Geom%QA_Rec(1,Geom%Num_QA) = 'm_MEF_EXO'
    Geom%QA_Rec(2,Geom%Num_QA) = 'MEF90-0.9.1'
    Call Date_And_Time(date = Geom%QA_Rec(3,Geom%Num_QA))
    Call Date_And_Time(time = Geom%QA_Rec(4,Geom%Num_QA))
    Call EXPQA(Geom%exoid, Geom%num_QA, Geom%QA_Rec, iErr)
    
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_EXO_Geom_Info
  
  Subroutine Show_EXO_Geom_Info(Geom)
    Type (EXO_Geom_Info), Intent(IN)               :: Geom

    Integer                                        :: i

    Write(*,100) 
    Write(*,101) Trim(Geom%filename)
    Write(*,102) Trim(Geom%title)
    Write(*,103) Geom%num_dim
    Write(*,104) Geom%num_nodes
    Write(*,105) Geom%num_elems
    Write(*,106) Geom%Num_Elem_blks
    Write(*,107) Geom%Num_Node_Sets
    Write(*,108) Geom%Num_Side_Sets

    Write(*,*)
    Write(*,200)
    Write(*,201) Geom%num_elem_blks
    Do i = 1, Geom%Num_Elem_blks
       Write(*,202) Geom%Elem_Blk(i)%ID, Geom%Elem_blk(i)%Type
       Write(*,203) Geom%Elem_Blk(i)%ID, Geom%Elem_blk(i)%Num_Elems
       Write(*,204) Geom%Elem_Blk(i)%ID, Geom%Elem_blk(i)%Num_Nodes_Per_Elem
       Write(*,205) Geom%Elem_Blk(i)%ID, Geom%Elem_blk(i)%Num_Attr
       Write(*,206, advance = 'no') Geom%Elem_Blk(i)%ID
       Write(*,*) Geom%Elem_blk(i)%Elem_ID
    End Do


    Write(*,*)
    Write(*,300)
    Write(*,301) Geom%num_node_sets
    Do i = 1, Geom%num_node_sets
       Write(*,302) Geom%Node_Set(i)%ID, Geom%Node_Set(i)%Num_Nodes
       Write(*,303, advance = 'no') Geom%Node_Set(i)%ID
       Write(*,*) Geom%Node_Set(i)%Node_ID
       Write(*,304) Geom%Node_Set(i)%ID, Geom%Node_Set(i)%Num_Dist_Factors
    End Do


    Write(*,*)
    Write(*,400)
    Write(*,401) Geom%num_side_sets

    Write(*,*)
    Write(*,500)
    Write(*,501) Geom%num_QA
    Do i = 1, Geom%num_QA
       Write(*,502) i, Trim(Geom%QA_rec(1,i))
       Write(*,503) i, Trim(Geom%QA_rec(2,i))
       Write(*,504) i, Trim(Geom%QA_rec(3,i))
       Write(*,505) i, Trim(Geom%QA_rec(4,i))
    End Do

100 Format('*** GLOBAL INFORMATIONS ***')
101 Format('    Filename ======================== ', A)
102 Format('    Title =========================== ', A)
103 Format('    Number of dimensions ============ ', I1)
104 Format('    Number of nodes ================= ', I6)
105 Format('    Number of elements ============== ', I6)
106 Format('    Number of elements blocks ======= ', I6)
107 Format('    Number of node sets ============= ', I6)
108 Format('    Number of side sets ============= ', I6)

200 Format('*** ELEMENT BLOCKS ***')
201 Format('    Number of blocks ================ ', I4)
202 Format('    Block ', I3, ' Elements type ========= ', A)
203 Format('    Block ', I3, ' Number of elements ==== ', I4)
204 Format('    Block ', I3, ' Number of nodes per elt ', I4)
205 Format('    Block ', I3, ' Number of attributes == ', I4)
206 Format('    Block ', I3, ' IDs: ')

300 Format('*** NODE SETS ***')
301 Format('    Number of sets ================== ', I4)
302 Format('    Set ', I3, ' Number of nodes ========= ', I4)
303 Format('    Set ', I3, ' IDs: ')
304 Format('    Set ', I3, ' Number of dist. factors = ', I4)

400 Format('*** SIDE SETS ***')
401 Format('    Number of side sets ============= ', I4)
    
500 Format('*** QA ***')
501 Format('    Number of QA Records =========== ', I2)
502 Format('    Rec ', I2, ' analysis code ============ ', A)
503 Format('    Rec ', I2, ' analysis QA desc. ======== ', A)
504 Format('    Rec ', I2, ' analysis date ============ ', A)
505 Format('    Rec ', I2, ' analysis time ============ ', A)
  End Subroutine Show_EXO_Geom_Info

  Subroutine Destroy_EXO_Geom_Info(Geom)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom

    Integer                                        :: iSet
    Integer                                        :: iBlk
    
    If (Associated(Geom%Elem_Blk)) Then
       Do iBlk = 1, Geom%Num_Elem_Blks
          If (Associated(Geom%Elem_Blk(iBlk)%Elem_ID)) Then
             DeAllocate(Geom%Elem_Blk(iBlk)%Elem_ID)
          End If
       End Do
       DeAllocate(Geom%Elem_Blk)
    End If
    
    If (Associated(Geom%Node_Set)) Then
       Do iSet = 1, Geom%Num_Node_Sets
          If (Associated(Geom%Node_Set(iSet)%Node_ID)) Then
             DeAllocate(Geom%Node_Set(iSet)%Node_ID)
          End If
       End Do
       DeAllocate(Geom%Node_Set)
    End If

    If (Associated(Geom%QA_Rec)) Then
       DeAllocate(Geom%QA_Rec)
    End If
  End Subroutine Destroy_EXO_Geom_Info

  Subroutine Read_EXO_Node_Coord_1D(Geom, Node_db, Num_DoF_per_node)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Type(Node1D), Dimension(:), Pointer            :: Node_db
    Integer, Intent(IN)                            :: Num_DoF_per_Node

    Integer                                        :: vers
    Real(Kind = Kr), Dimension(:), Pointer         :: X
    Real(Kind = Kr)                                :: Y, Z
    Integer                                        :: iDim, iN
    Integer                                        :: iSet
    Integer                                        :: iErr

    If (Geom%Num_Dim /= 1) Then
       Write(*,*) 'ERROR: not a 1D EXO file'
       Stop 
    End If

    Allocate (Node_db(Geom%Num_Nodes * Num_DoF_per_Node))
    
    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    Allocate (X(Geom%Num_Nodes))
    Call EXGCOR(Geom%exoid, X, Y, Z, iErr)
    Do iN = 1, Geom%Num_Nodes
       Node_db(Num_DoF_per_node*(iN-1)+1 : Num_DoF_per_node * iN)%Coord      &
            & = X(iN)
    End Do
    DeAllocate (X)

    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
    Geom%Numbering = Numbering_PerNodes
  End Subroutine Read_EXO_Node_Coord_1D

  Subroutine Read_EXO_Node_Coord_2D(Geom, Node_db, Num_DoF_per_node)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Type(Node2D), Dimension(:), Pointer            :: Node_db
    Integer, Intent(IN)                            :: Num_DoF_per_Node

    Integer                                        :: vers
    Real(Kind = Kr), Dimension(:), Pointer         :: X, Y
    Real                                           :: Z
    Integer                                        :: iDim, iN
    Integer                                        :: iSet
    Integer                                        :: iErr

    If (Geom%Num_Dim /= 2) Then
       Write(*,*) 'ERROR: not a 2D EXO file'
       Stop 
    End If

    Allocate (Node_db(Geom%Num_Nodes * Num_DoF_per_Node))
    
    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)
    Allocate (X(Geom%Num_Nodes))
    Allocate (Y(Geom%Num_Nodes))
    Call EXGCOR(Geom%exoid, X, Y, Z, iErr)
    Do iN = 1, Geom%Num_Nodes
       Node_db(Num_DoF_per_node*(iN-1) + 1 : Num_DoF_per_node * iN)%Coord%X = &
            & X(iN)
       Node_db(Num_DoF_per_node*(iN-1) + 1 : Num_DoF_per_node * iN)%Coord%Y = &
            & Y(iN)
    End Do
    DeAllocate (X,Y)

    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
    Geom%Numbering = Numbering_PerNodes
  End Subroutine Read_EXO_Node_Coord_2D

  Subroutine Read_EXO_Node_Coord_3D(Geom, Node_db, Num_DoF_per_node)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Type(Node3D), Dimension(:), Pointer            :: Node_db
    Integer, Intent(IN)                            :: Num_DoF_per_Node

    Integer                                        :: vers
    Real(Kind = Kr), Dimension(:), Pointer         :: X, Y, Z
    Integer                                        :: iDim, iN
    Integer                                        :: iSet
    Integer                                        :: iErr

    If (Geom%Num_Dim /= 3) Then
       Write(*,*) 'ERROR: not a 3D EXO file'
       Stop 
    End If

    Allocate (Node_db(Geom%Num_Nodes * Num_DoF_per_Node))
    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    Allocate (X(Geom%Num_Nodes))
    Allocate (Y(Geom%Num_Nodes))
    Allocate (Z(Geom%Num_Nodes))
    Call EXGCOR(Geom%exoid, X, Y, Z, iErr)
    Do iN = 1, Geom%Num_Nodes
       Node_db(Num_DoF_per_node*(iN-1) + 1 : Num_DoF_per_node * iN)%Coord%X = &
            & X(iN)
       Node_db(Num_DoF_per_node*(iN-1) + 1 : Num_DoF_per_node * iN)%Coord%Y = &
            & Y(iN)
       Node_db(Num_DoF_per_node*(iN-1) + 1 : Num_DoF_per_node * iN)%Coord%Z = &
            & Z(iN)
    End Do
    DeAllocate (X,Y,Z)

    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
    Geom%Numbering = Numbering_PerNodes
  End Subroutine Read_EXO_Node_Coord_3D

  Subroutine Write_EXO_Node_Coord_1D(Geom, Node_db, Num_DoF_per_node)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Type(Node1D), Dimension(:), Pointer            :: Node_db
    Integer, Intent(IN)                            :: Num_DoF_per_Node

    Real(Kind = Kr)                                :: Y, Z
    Real(Kind = Kr)                                :: Vers
    Integer                                        :: iErr

    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    Y = 0.0_Kr
    Z = 0.0_Kr

    Call EXPCOR(Geom%exoid,                                                   &
         & Node_db(1:Geom%Num_Nodes*Num_Dof_per_Node:Num_Dof_per_Node)%Coord, &
         & Y, Z, iErr)

    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_EXO_Node_Coord_1D

  Subroutine Write_EXO_Node_Coord_2D(Geom, Node_db, Num_DoF_per_node)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Type(Node2D), Dimension(:), Pointer            :: Node_db
    Integer, Intent(IN)                            :: Num_DoF_per_Node

    Real(Kind = Kr)                                :: Z 
    Real(Kind = Kr)                                :: Vers
    Integer                                        :: iErr

    Z = 0.0_Kr
    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    Call EXPCOR(Geom%exoid,                                                   &
         &Node_db(1:Geom%Num_Nodes*Num_Dof_per_Node:Num_DoF_per_Node)%Coord%X,&
         &Node_db(1:Geom%Num_Nodes*Num_Dof_per_Node:Num_DoF_per_Node)%Coord%Y,&
         & Z, iErr)

    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_EXO_Node_Coord_2D

  Subroutine Write_EXO_Node_Coord_3D(Geom, Node_db, Num_DoF_per_node)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Type(Node3D), Dimension(:), Pointer            :: Node_db
    Integer,Intent(IN)                             :: Num_DoF_per_Node

    Real(Kind = Kr)                                :: Vers
    Integer                                        :: iErr

    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    Call EXPCOR(Geom%exoid,                                                   &
         &Node_db(1:Geom%Num_Nodes*Num_Dof_per_Node:Num_DoF_per_Node)%Coord%X,&
         &Node_db(1:Geom%Num_Nodes*Num_Dof_per_Node:Num_DoF_per_Node)%Coord%Y,&
         &Node_db(1:Geom%Num_Nodes*Num_Dof_per_Node:Num_DoF_per_Node)%Coord%Z,&
         &iErr)

    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_EXO_Node_Coord_3D

  Subroutine Read_EXO_Connect_2D(Geom, Elem_db)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Type (Element2D), Dimension(:), Pointer        :: Elem_db

    Integer, Dimension(:,:), Pointer               :: Tmp_Connect
    Integer                                        :: iBlock, Offset, iE
    Real(Kind = Kr)                                :: Vers
    Integer                                        :: iErr
    Integer                                        :: Nb_DoF

    If (Geom%Num_Dim /= 2) Then
       Write(*,*) 'ERROR: Dimension mismatch in Read_EXO_Connect_2D'
       Stop
    End If

    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    Allocate (Elem_db(Geom%Num_Elems))
    Select Case (Geom%Numbering)
    Case(Numbering_PerCoord)       
       Offset = 0
       Do iBlock = 1, Geom%Num_Elem_Blks
          Allocate(Tmp_Connect(Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem,      &
               & Geom%Elem_Blk(iBlock)%Num_Elems))
          Call EXGELC(Geom%exoid, Geom%Elem_Blk(iBlock)%ID, Tmp_Connect, iErr)
          Do iE = 1, Geom%Elem_Blk(iBlock)%Num_Elems
             Elem_db(iE+Offset)%NB_DoF =                                      &
                  & 2 * Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem
             Elem_db(iE+Offset)%ID_EL = iBlock
             Allocate (Elem_db(iE+Offset)%ID_DoF(Elem_db(iE+Offset)%NB_DoF))
             Elem_db(iE+Offset)%ID_DoF = (/Tmp_Connect(:, iE),                &
                  & Tmp_Connect(:, iE) + Geom%Num_nodes /)
          End Do
          Offset = Offset +Geom%Elem_Blk(iBlock)%Num_Elems
          DeAllocate (Tmp_Connect)
       End Do
    Case(Numbering_PerNodes)
       Offset = 0
       Do iBlock = 1, Geom%Num_Elem_Blks
          Allocate(Tmp_Connect(Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem,      &
               & Geom%Elem_Blk(iBlock)%Num_Elems))
          Call EXGELC(Geom%exoid, Geom%Elem_Blk(iBlock)%ID, Tmp_Connect, iErr)
          Do iE = 1, Geom%Elem_Blk(iBlock)%Num_Elems
             Elem_db(iE+Offset)%NB_DoF =                                      &
                  & 2 * Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem
             Elem_db(iE+Offset)%ID_EL = iBlock
             Allocate (Elem_db(iE+Offset)%ID_DoF(Elem_db(iE+Offset)%NB_DoF))
             Elem_db(iE+Offset)%ID_DoF = (/ 2 * Tmp_Connect(:, iE) - 1,       &
                  & 2 * Tmp_Connect(:, iE) /)
          End Do
          Offset = Offset +Geom%Elem_Blk(iBlock)%Num_Elems
          DeAllocate (Tmp_Connect)
       End Do
    Case Default
       Write(*,*) 'ERROR: Read_EXO_Connect_2D Numbering scheme unknow: ',     &
            & Geom%Numbering
       STOP
    End Select
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Read_EXO_Connect_2D

  Subroutine Read_EXO_Connect_2D_Scal(Geom, Elem_db)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Type (Element2D_Scal), Dimension(:), Pointer   :: Elem_db

    Integer, Dimension(:,:), Pointer               :: Tmp_Connect
    Integer                                        :: iBlock, Offset, iE
    Real(Kind = Kr)                                :: Vers
    Integer                                        :: iErr
    Integer                                        :: Nb_DoF

    If (Geom%Num_Dim /= 2) Then
       Write(*,*) 'ERROR: Dimension mismatch in Read_EXO_Connect_2D_Scal'
    End If
    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    Allocate (Elem_db(Geom%Num_Elems))

    Select Case (Geom%Numbering)
    Case(Numbering_PerCoord, Numbering_PerNodes)       
!!! In the scalar case, Numbering_PerNodes and Numbering_PerCoord 
!!! are equivalent
       Offset = 0
       Do iBlock = 1, Geom%Num_Elem_Blks
          Allocate(Tmp_Connect(Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem,      &
               & Geom%Elem_Blk(iBlock)%Num_Elems))
          Call EXGELC (Geom%exoid, Geom%Elem_Blk(iBlock)%ID, Tmp_Connect, iErr)
          Do iE = 1, Geom%Elem_Blk(iBlock)%Num_Elems
             Elem_db(iE+Offset)%NB_DoF =                                      &
                  & Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem
             Elem_db(iE+Offset)%ID_EL = iBlock
             Allocate (Elem_db(iE+Offset)%ID_DoF(Elem_db(iE+Offset)%NB_DoF))
             Elem_db(iE+Offset)%ID_DoF = Tmp_Connect(:, iE)
          End Do
          Offset = Offset +Geom%Elem_Blk(iBlock)%Num_Elems
          DeAllocate (Tmp_Connect)
       End Do
    Case Default
       Write(*,*) 'ERROR: Read_EXO_Connect_2D_Scal Numbering scheme unknow: ',&
            & Geom%Numbering
       STOP
    End Select
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Read_EXO_Connect_2D_Scal

  Subroutine Read_EXO_Connect_2D_Elast(Geom, Elem_db)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Type (Element2D_Elast), Dimension(:), Pointer  :: Elem_db

    Integer, Dimension(:,:), Pointer               :: Tmp_Connect
    Integer                                        :: iBlock, Offset, iE
    Real(Kind = Kr)                                :: Vers
    Integer                                        :: iErr
    Integer                                        :: Nb_DoF

    If (Geom%Num_Dim /= 2) Then
       Write(*,*) 'ERROR: Dimension mismatch in Read_EXO_Connect_2D_Elast'
       Stop
    End If

    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    Allocate (Elem_db(Geom%Num_Elems))
    Offset = 0
    Select Case (Geom%Numbering)
    Case(Numbering_PerCoord)       
       Do iBlock = 1, Geom%Num_Elem_Blks
          Allocate(Tmp_Connect(Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem,      &
               & Geom%Elem_Blk(iBlock)%Num_Elems))
          Call EXGELC (Geom%exoid, Geom%Elem_Blk(iBlock)%ID, Tmp_Connect, iErr)
          Do iE = 1, Geom%Elem_Blk(iBlock)%Num_Elems
             Elem_db(iE+Offset)%NB_DoF =                                      &
                  & 2 * Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem
             Elem_db(iE+Offset)%ID_EL = iBlock
             Allocate (Elem_db(iE+Offset)%ID_DoF(Elem_db(iE+Offset)%NB_DoF))
             Elem_db(iE+Offset)%ID_DoF = (/Tmp_Connect(:, iE),                &
                  & Tmp_Connect(:, iE) + Geom%Num_nodes /)
          End Do
          Offset = Offset +Geom%Elem_Blk(iBlock)%Num_Elems
          DeAllocate (Tmp_Connect)
       End Do
    Case(Numbering_PerNodes)       
       Do iBlock = 1, Geom%Num_Elem_Blks
          Allocate(Tmp_Connect(Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem,      &
               & Geom%Elem_Blk(iBlock)%Num_Elems))
          Call EXGELC (Geom%exoid, Geom%Elem_Blk(iBlock)%ID, Tmp_Connect, iErr)
          Do iE = 1, Geom%Elem_Blk(iBlock)%Num_Elems
             Elem_db(iE+Offset)%NB_DoF =                                      &
                  & 2 * Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem
             Elem_db(iE+Offset)%ID_EL = iBlock
             Allocate (Elem_db(iE+Offset)%ID_DoF(Elem_db(iE+Offset)%NB_DoF))
             Elem_db(iE+Offset)%ID_DoF = (/ 2 * Tmp_Connect(:, iE) - 1,       &
                  & 2 * Tmp_Connect(:, iE) /)
          End Do
          Offset = Offset +Geom%Elem_Blk(iBlock)%Num_Elems
          DeAllocate (Tmp_Connect)
       End Do
    Case Default
       Write(*,*) 'ERROR: Read_EXO_Connect_2D_Elast Numbering scheme unknown',&
            & Geom%Numbering
       STOP
    End Select
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Read_EXO_Connect_2D_Elast

  Subroutine Read_EXO_Connect_3D(Geom, Elem_db)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Type (Element3D), Dimension(:), Pointer        :: Elem_db

    Integer, Dimension(:,:), Pointer               :: Tmp_Connect
    Integer                                        :: iBlock, Offset, iE
    Real(Kind = Kr)                                :: Vers
    Integer                                        :: iErr
    Integer                                        :: Nb_DoF

    If (Geom%Num_Dim /= 3) Then
       Write(*,*) 'ERROR: Dimension mismatch in Read_EXO_Connect_3D'
       Stop
    End If

    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    Allocate (Elem_db(Geom%Num_Elems))
    Offset = 0
    Select Case (Geom%Numbering)
    Case(Numbering_PerCoord)       
       Do iBlock = 1, Geom%Num_Elem_Blks
          Allocate(Tmp_Connect(Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem,      &
               & Geom%Elem_Blk(iBlock)%Num_Elems))
          Call EXGELC (Geom%exoid, Geom%Elem_Blk(iBlock)%ID, Tmp_Connect, iErr)
          Do iE = 1, Geom%Elem_Blk(iBlock)%Num_Elems
             Elem_db(iE+Offset)%NB_DoF = Geom%Num_Dim *                       &
                  & Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem
             Elem_db(iE+Offset)%ID_EL = iBlock
             Allocate (Elem_db(iE+Offset)%ID_DoF(Elem_db(iE+Offset)%NB_DoF))
             Elem_db(iE+Offset)%ID_DoF = (/Tmp_Connect(:, iE),                &
                  & Tmp_Connect(:, iE) + Geom%Num_nodes,                      &
                  & Tmp_Connect(:, iE) + 2 * Geom%Num_nodes /)
          End Do
          Offset = Offset +Geom%Elem_Blk(iBlock)%Num_Elems
          DeAllocate (Tmp_Connect)
       End Do
    Case(Numbering_PerNodes)       
       Do iBlock = 1, Geom%Num_Elem_Blks
          Allocate(Tmp_Connect(Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem,      &
               & Geom%Elem_Blk(iBlock)%Num_Elems))
          Call EXGELC (Geom%exoid, Geom%Elem_Blk(iBlock)%ID, Tmp_Connect, iErr)
          Do iE = 1, Geom%Elem_Blk(iBlock)%Num_Elems
             Elem_db(iE+Offset)%NB_DoF = Geom%Num_Dim *                       &
                  & Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem
             Elem_db(iE+Offset)%ID_EL = iBlock
             Allocate (Elem_db(iE+Offset)%ID_DoF(Elem_db(iE+Offset)%NB_DoF))
             Elem_db(iE+Offset)%ID_DoF = (/ 3 * Tmp_Connect(:, iE) - 2,       &
                  & 3 * Tmp_Connect(:, iE) - 1, 3 * Tmp_Connect(:, iE) /)
          End Do
          Offset = Offset +Geom%Elem_Blk(iBlock)%Num_Elems
          DeAllocate (Tmp_Connect)
       End Do
    Case Default
       Write(*,*) 'ERROR: Read_EXO_Connect_3D Numbering scheme unknown',&
            & Geom%Numbering
       STOP
    End Select

    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Read_EXO_Connect_3D

  Subroutine Read_EXO_Connect_3D_Scal(Geom, Elem_db)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Type (Element3D_Scal), Dimension(:), Pointer   :: Elem_db

    Integer, Dimension(:,:), Pointer               :: Tmp_Connect
    Integer                                        :: iBlock, Offset, iE
    Real(Kind = Kr)                                :: Vers
    Integer                                        :: iErr
    Integer                                        :: Nb_DoF

    If (Geom%Num_Dim /= 3) Then
       Write(*,*) 'ERROR: Dimension mismatch in Read_EXO_Connect_3D_Scal'
       Stop
    End If

    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    Allocate (Elem_db(Geom%Num_Elems))
    Offset = 0
    Select Case (Geom%Numbering)
    Case(Numbering_PerCoord, Numbering_PerNodes)       
!!! Scalar case Numbering_PerNode <=> Numbering_PerCoord
       Do iBlock = 1, Geom%Num_Elem_Blks
          Allocate(Tmp_Connect(Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem,      &
               & Geom%Elem_Blk(iBlock)%Num_Elems))
          Call EXGELC (Geom%exoid, Geom%Elem_Blk(iBlock)%ID, Tmp_Connect, iErr)
          Do iE = 1, Geom%Elem_Blk(iBlock)%Num_Elems
             Elem_db(iE+Offset)%NB_DoF =                                      &
                  & Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem
             Elem_db(iE+Offset)%ID_EL = iBlock
             Allocate (Elem_db(iE+Offset)%ID_DoF(Elem_db(iE+Offset)%NB_DoF))
             Elem_db(iE+Offset)%ID_DoF = Tmp_Connect(:, iE)
          End Do
          Offset = Offset +Geom%Elem_Blk(iBlock)%Num_Elems
          DeAllocate (Tmp_Connect)
       End Do
    Case Default
       Write(*,*) 'ERROR: Read_EXO_Connect_3D_Scal Numbering scheme unknown',&
            & Geom%Numbering
       STOP
    End Select
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Read_EXO_Connect_3D_Scal

  Subroutine Read_EXO_Connect_3D_Elast(Geom, Elem_db)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Type (Element3D_Elast), Dimension(:), Pointer  :: Elem_db

    Integer, Dimension(:,:), Pointer               :: Tmp_Connect
    Integer                                        :: iBlock, Offset, iE
    Real(Kind = Kr)                                :: Vers
    Integer                                        :: iErr
    Integer                                        :: Nb_DoF

    If (Geom%Num_Dim /= 3) Then
       Write(*,*) 'ERROR: Dimension mismatch in Read_EXO_Connect_3D_Elast'
       Stop
    End If

    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

!    Allocate (Elem_db(Geom%Num_Elems))
    Allocate (Elem_db(Geom%Num_Elems))
    Offset = 0
    Select Case (Geom%Numbering)
    Case(Numbering_PerCoord)       
       Do iBlock = 1, Geom%Num_Elem_Blks
          Allocate(Tmp_Connect(Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem,      &
               & Geom%Elem_Blk(iBlock)%Num_Elems))
          Call EXGELC (Geom%exoid, Geom%Elem_Blk(iBlock)%ID, Tmp_Connect, iErr)
          Do iE = 1, Geom%Elem_Blk(iBlock)%Num_Elems
             Elem_db(iE+Offset)%NB_DoF = Geom%Num_Dim *                       &
                  & Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem
             Elem_db(iE+Offset)%ID_EL = iBlock
             Allocate (Elem_db(iE+Offset)%ID_DoF(Elem_db(iE+Offset)%NB_DoF))
             Elem_db(iE+Offset)%ID_DoF = (/Tmp_Connect(:, iE),                &
                  & Tmp_Connect(:, iE) + Geom%Num_nodes,                      &
                  & Tmp_Connect(:, iE) + 2 * Geom%Num_nodes /)
          End Do
          Offset = Offset +Geom%Elem_Blk(iBlock)%Num_Elems
          DeAllocate (Tmp_Connect)
       End Do
    Case(Numbering_PerNodes)       
!       Print*, 'Read_Connect_3D_Elast, Numbering_PerNodes'
!       Print*, 'Geom%Num_Dim: ', Geom%Num_Dim
       Do iBlock = 1, Geom%Num_Elem_Blks
          Allocate(Tmp_Connect(Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem,      &
               & Geom%Elem_Blk(iBlock)%Num_Elems))
          Call EXGELC (Geom%exoid, Geom%Elem_Blk(iBlock)%ID, Tmp_Connect, iErr)
          Do iE = 1, Geom%Elem_Blk(iBlock)%Num_Elems
             Elem_db(iE+Offset)%NB_DoF = Geom%Num_Dim *                       &
                  & Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem
             Elem_db(iE+Offset)%ID_EL = iBlock
             Allocate (Elem_db(iE+Offset)%ID_DoF(Elem_db(iE+Offset)%NB_DoF))
             Elem_db(iE+Offset)%ID_DoF = (/ 3 * Tmp_Connect(:, iE) - 2,       &
                  & 3 * Tmp_Connect(:, iE) - 1, 3 * Tmp_Connect(:, iE) /)
!             Print*, 'Tmp_Connect: ', Tmp_Connect(:,iE)
!             Print*, 'ID_DoF: ', Elem_db(iE+Offset)%ID_DoF
          End Do
          Offset = Offset + Geom%Elem_Blk(iBlock)%Num_Elems
          DeAllocate (Tmp_Connect)
       End Do
    Case Default
       Write(*,*) 'ERROR: Read_EXO_Connect_3D Numbering scheme unknown',&
            & Geom%Numbering
       STOP
    End Select
  End Subroutine Read_EXO_Connect_3D_Elast

  Subroutine Write_EXO_Connect_2D(Geom, Elem_db)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Type (Element2D), Dimension(:), Pointer        :: Elem_db
    
    Integer, Dimension(:,:), Pointer               :: Tmp_Connect
    Integer                                        :: iBlock, iE, Offset
    Real(Kind = Kr)                                :: Vers
    Integer                                        :: iErr
    Integer                                        :: Nb_DoF, i
    
    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)
    
    Offset = 0
    Select Case (Geom%Numbering)
    Case(Numbering_PerCoord)       
       Do iBlock = 1, Geom%Num_Elem_Blks
          Nb_DoF = 2 * Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem
          Allocate(Tmp_Connect(Nb_DoF/2, Geom%Elem_Blk(iBlock)%Num_Elems))
          Do iE = 1, Geom%Elem_Blk(iBlock)%Num_Elems
             Tmp_Connect(:, iE) = Elem_db(iE+Offset)%ID_DoF(1:Nb_DoF/2)
          End Do
          
          Call EXPELC (Geom%exoid, Geom%Elem_Blk(iBlock)%ID, Tmp_Connect, iErr)
          DeAllocate (Tmp_Connect)
          Offset = Offset + Geom%Elem_Blk(iBlock)%Num_Elems
       End Do
    Case(Numbering_PerNodes)
       Do iBlock = 1, Geom%Num_Elem_Blks
          Nb_DoF = 2 * Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem 
          Allocate(Tmp_Connect(Nb_DoF/2, Geom%Elem_Blk(iBlock)%Num_Elems))
          Do iE = 1, Geom%Elem_Blk(iBlock)%Num_Elems
             Tmp_Connect(:, iE) =                                             &
                  & Elem_db(iE+Offset)%ID_DoF(1:Nb_DoF/2 )/2 +1
          End Do
          
          Call EXPELC (Geom%exoid, Geom%Elem_Blk(iBlock)%ID, Tmp_Connect, iErr)
          DeAllocate (Tmp_Connect)
          Offset = Offset + Geom%Elem_Blk(iBlock)%Num_Elems
       End Do
    Case Default
       Write(*,*) 'ERROR: Write_EXO_Connect_2D Numbering scheme unknown',&
            & Geom%Numbering
       STOP
    End Select
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_EXO_Connect_2D

  Subroutine Write_EXO_Connect_2D_Scal(Geom, Elem_db)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Type (Element2D_Scal), Dimension(:), Pointer   :: Elem_db

    Integer, Dimension(:,:), Pointer               :: Tmp_Connect
    Integer                                        :: iBlock, iE, Offset
    Real(Kind = Kr)                                :: Vers
    Integer                                        :: iErr
    Integer                                        :: Nb_DoF
    
    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)
    
    Offset = 0
    Select Case (Geom%Numbering)
    Case(Numbering_PerCoord, Numbering_PerNodes)       
!!! Scalar case Numbering_PerNode <=> Numbering_PerCoord
       Do iBlock = 1, Geom%Num_Elem_Blks
          Nb_DoF = Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem
          Allocate(Tmp_Connect(Nb_DoF, Geom%Elem_Blk(iBlock)%Num_Elems))
          Do iE = 1, Geom%Elem_Blk(iBlock)%Num_Elems
             Tmp_Connect(:, iE) = Elem_db(iE+Offset)%ID_DoF(1:Nb_DoF)
          End Do
          
          Call EXPELC (Geom%exoid, Geom%Elem_Blk(iBlock)%ID, Tmp_Connect, iErr)
          DeAllocate (Tmp_Connect)
          Offset = Offset + Geom%Elem_Blk(iBlock)%Num_Elems
       End Do
    Case Default
       Write(*,*) 'ERROR: Write_EXO_Connect_2D_Scal Numbering scheme unknown',&
            & Geom%Numbering
       STOP
    End Select

    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_EXO_Connect_2D_Scal

  Subroutine Write_EXO_Connect_2D_Elast(Geom, Elem_db)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Type (Element2D_Elast), Dimension(:), Pointer        :: Elem_db

    Integer, Dimension(:,:), Pointer               :: Tmp_Connect
    Integer                                        :: iBlock, iE, Offset
    Real(Kind = Kr)                                :: Vers
    Integer                                        :: iErr
    Integer                                        :: Nb_DoF
    
    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)
    
    Offset = 0
    Select Case (Geom%Numbering)
    Case(Numbering_PerCoord)       
       Do iBlock = 1, Geom%Num_Elem_Blks
          Nb_DoF = 2 * Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem
          Allocate(Tmp_Connect(Nb_DoF/2, Geom%Elem_Blk(iBlock)%Num_Elems))
          Do iE = 1, Geom%Elem_Blk(iBlock)%Num_Elems
             Tmp_Connect(:, iE) = Elem_db(iE+Offset)%ID_DoF(1:Nb_DoF/2)
          End Do
          
          Call EXPELC (Geom%exoid, Geom%Elem_Blk(iBlock)%ID, Tmp_Connect, iErr)
          DeAllocate (Tmp_Connect)
          Offset = Offset + Geom%Elem_Blk(iBlock)%Num_Elems
       End Do
    Case(Numbering_PerNodes)       
       Do iBlock = 1, Geom%Num_Elem_Blks
          Nb_DoF = 2 * Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem
          Allocate(Tmp_Connect(Nb_DoF/2, Geom%Elem_Blk(iBlock)%Num_Elems))
          Do iE = 1, Geom%Elem_Blk(iBlock)%Num_Elems
             Tmp_Connect(:, iE) =                                             &
                  & Elem_db(iE+Offset)%ID_DoF(1:Nb_DoF/2 )/2+1
          End Do
          
          Call EXPELC (Geom%exoid, Geom%Elem_Blk(iBlock)%ID, Tmp_Connect, iErr)
          DeAllocate (Tmp_Connect)
          Offset = Offset + Geom%Elem_Blk(iBlock)%Num_Elems
       End Do
    Case Default
       Write(*,*) 'ERROR:Write_EXO_Connect_2D_Elast Numbering scheme unknown',&
            & Geom%Numbering
       STOP
    End Select

    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_EXO_Connect_2D_Elast

  Subroutine Write_EXO_Connect_3D(Geom, Elem_db)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Type (Element3D), Dimension(:), Pointer        :: Elem_db
    
    Integer, Dimension(:,:), Pointer               :: Tmp_Connect
    Integer                                        :: iBlock, iE, Offset
    Real(Kind = Kr)                                :: Vers
    Integer                                        :: iErr
    Integer                                        :: Nb_DoF
    
    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)
    
    Offset = 0
    Select Case (Geom%Numbering)
    Case(Numbering_PerCoord)       
       Do iBlock = 1, Geom%Num_Elem_Blks
          Nb_DoF = 3 * Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem
          Allocate(Tmp_Connect(Nb_DoF/3, Geom%Elem_Blk(iBlock)%Num_Elems))
          Do iE = 1, Geom%Elem_Blk(iBlock)%Num_Elems
             Tmp_Connect(:, iE) = Elem_db(iE+Offset)%ID_DoF(1:Nb_DoF/3)
          End Do
          
          Call EXPELC (Geom%exoid, Geom%Elem_Blk(iBlock)%ID, Tmp_Connect, iErr)
          DeAllocate (Tmp_Connect)
          Offset = Offset + Geom%Elem_Blk(iBlock)%Num_Elems
       End Do
    Case(Numbering_PerNodes)       
       Do iBlock = 1, Geom%Num_Elem_Blks
          Nb_DoF = 3 * Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem
          Allocate(Tmp_Connect(Nb_DoF/3, Geom%Elem_Blk(iBlock)%Num_Elems))
          Do iE = 1, Geom%Elem_Blk(iBlock)%Num_Elems
             Tmp_Connect(:, iE) =                                             &
                  & Elem_db(iE+Offset)%ID_DoF(1:Nb_DoF/3 )/3+1
          End Do
          
          Call EXPELC (Geom%exoid, Geom%Elem_Blk(iBlock)%ID, Tmp_Connect, iErr)
          DeAllocate (Tmp_Connect)
          Offset = Offset + Geom%Elem_Blk(iBlock)%Num_Elems
       End Do
        Case Default
       Write(*,*) 'ERROR: Write_EXO_Connect_3D Numbering scheme unknown',&
            & Geom%Numbering
       STOP
    End Select

    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_EXO_Connect_3D

  Subroutine Write_EXO_Connect_3D_Scal(Geom, Elem_db)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Type (Element3D_Scal), Dimension(:), Pointer   :: Elem_db

    Integer, Dimension(:,:), Pointer               :: Tmp_Connect
    Integer                                        :: iBlock, iE, Offset
    Real(Kind = Kr)                                :: Vers
    Integer                                        :: iErr
    Integer                                        :: Nb_DoF
    
    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)
    
    Offset = 0
    Select Case (Geom%Numbering)
    Case(Numbering_PerCoord, Numbering_PerNodes)       
!!! Scalar case Numbering_PerNode <=> Numbering_PerCoord
       Do iBlock = 1, Geom%Num_Elem_Blks
          Nb_DoF = Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem
          Allocate(Tmp_Connect(Nb_DoF, Geom%Elem_Blk(iBlock)%Num_Elems))
          Do iE = 1, Geom%Elem_Blk(iBlock)%Num_Elems
             Tmp_Connect(:, iE) = Elem_db(iE+Offset)%ID_DoF(1:Nb_DoF)
          End Do
          
          Call EXPELC (Geom%exoid, Geom%Elem_Blk(iBlock)%ID, Tmp_Connect, iErr)
          DeAllocate (Tmp_Connect)
          Offset = Offset + Geom%Elem_Blk(iBlock)%Num_Elems
       End Do
    Case Default
       Write(*,*) 'ERROR: Write_EXO_Connect_3D_Scal Numbering scheme unknown',&
            & Geom%Numbering
       STOP
    End Select

    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_EXO_Connect_3D_Scal

  Subroutine Write_EXO_Connect_3D_Elast(Geom, Elem_db)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Type (Element3D_Elast), Dimension(:), Pointer  :: Elem_db

    Integer, Dimension(:,:), Pointer               :: Tmp_Connect
    Integer                                        :: iBlock, iE, Offset
    Real(Kind = Kr)                                :: Vers
    Integer                                        :: iErr
    Integer                                        :: Nb_DoF
    
    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)
    
    Offset = 0
    Select Case (Geom%Numbering)
    Case(Numbering_PerCoord)       
       Do iBlock = 1, Geom%Num_Elem_Blks
          Nb_DoF = 3 * Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem
          Allocate(Tmp_Connect(Nb_DoF/3, Geom%Elem_Blk(iBlock)%Num_Elems))
          Do iE = 1, Geom%Elem_Blk(iBlock)%Num_Elems
             Tmp_Connect(:, iE) = Elem_db(iE+Offset)%ID_DoF(1:Nb_DoF/3)
          End Do
          
          Call EXPELC (Geom%exoid, Geom%Elem_Blk(iBlock)%ID, Tmp_Connect, iErr)
          DeAllocate (Tmp_Connect)
          Offset = Offset + Geom%Elem_Blk(iBlock)%Num_Elems
       End Do
    Case(Numbering_PerNodes)       
       Do iBlock = 1, Geom%Num_Elem_Blks
          Nb_DoF = 3 * Geom%Elem_Blk(iBlock)%Num_Nodes_per_Elem
          Allocate(Tmp_Connect(Nb_DoF/3, Geom%Elem_Blk(iBlock)%Num_Elems))
          Do iE = 1, Geom%Elem_Blk(iBlock)%Num_Elems
             Tmp_Connect(:, iE) = Elem_db(iE+Offset)%ID_DoF(1:Nb_DoF/3)/3+1
          End Do
          
          Call EXPELC (Geom%exoid, Geom%Elem_Blk(iBlock)%ID, Tmp_Connect, iErr)
          DeAllocate (Tmp_Connect)
          Offset = Offset + Geom%Elem_Blk(iBlock)%Num_Elems
       End Do
    Case Default
       Write(*,*) 'ERROR: Read_EXO_Connect_3D Numbering scheme unknown',&
            & Geom%Numbering
       STOP
    End Select
       
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_EXO_Connect_3D_Elast

!!!
!!! RESULT FILES
!!!

!!! READ
!!! EXO stores values per coordinates
!!! V_PerNode = (/ (V_PerCoord(i:Size:Num_Nodes), i=1:NumDim /)
!!! V_PerCoord = (/ (V_PerNodes(i:Size:NumDim), i=1, Num_Nodes/)
!!! 
  Subroutine Read_EXO_Result_Global(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Real(Kind = Kr)                                :: Res

    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Real(Kind = Kr), Dimension(:), Pointer         :: Tmp_Res
    Integer                                        :: Num_Vars, Num_TS
    Real(Kind = Kr)                               :: fDum
    Character                                     :: cDum

    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)
    ! Get the number of global variables stored in the database    
    Call EXGVP(Geom%exoid, 'G', Num_Vars, iErr)
    Allocate(Tmp_Res(Num_Vars))

    ! Read All the global variables at time step TS and copy the right one
    ! into Res
    Call EXGGV(Geom%exoid, TS, Num_Vars, Tmp_Res, iErr)
    Res = Tmp_Res(Idx)
    DeAllocate (Tmp_Res)
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Read_EXO_Result_Global

  Subroutine Read_EXO_Result_Scal_Nodes(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Real(Kind = Kr), Dimension(:), Pointer         :: Res

    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    If (Associated(Res)) Then
       DeAllocate(Res)
    End If
    Allocate(Res(Geom%Num_Nodes))
    Call EXGNV(Geom%exoid, TS, Idx, Geom%Num_Nodes, Res, iErr)
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Read_EXO_Result_Scal_Nodes

  Subroutine Read_EXO_Result_Ptr_Nodes(Geom, Idx, TS, Res, Num_Rec)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Real(Kind = Kr), Dimension(:), Pointer         :: Res
    Integer                                        :: Num_Rec

    Integer                                        :: iErr, i
    Real(Kind = Kr)                                :: Vers
    Integer                                        :: iRec
    Real(Kind = Kr), Dimension(:), Pointer         :: Tmp_Res

    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)
    
    If (Associated(Res)) Then
       DeAllocate(Res)
    End If
    Allocate(Res(Geom%Num_Nodes*Num_Rec))

    Allocate(Tmp_Res(Geom%Num_Nodes))
    Select Case (Geom%Numbering)
    Case(Numbering_PerCoord)
       Do iRec = 0, Num_Rec-1
          Call EXGNV(Geom%exoid, TS, Idx + iRec, Geom%Num_Nodes,              &
               & Res(1+iRec * Geom%Num_Nodes:(iRec+1) * Geom%Num_Nodes), iErr)
       End Do
    Case(Numbering_PerNodes)
       Do iRec = 1, Num_Rec
          Call EXGNV(Geom%exoid, TS, Idx + iRec-1 , Geom%Num_Nodes,           &
               & Tmp_Res, iErr)
!               & Res(iRec:Num_Rec*Geom%Num_Nodes:Num_Rec), iErr)
          Res(iRec:Num_Rec*Geom%Num_Nodes:Num_Rec) = Tmp_Res

       End Do

    Case Default
       Write(*,*) 'ERROR: Read_EXO_Result_Ptr_Nodes Numbering scheme unknown',&
            & Geom%Numbering
       STOP
    End Select
    DeAllocate(Tmp_Res)
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Read_EXO_Result_Ptr_Nodes

  Subroutine Read_EXO_Result_Vect2D_Nodes(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Type(Vect2D), Dimension(:), Pointer            :: Res

    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

!!$    If (Geom%Num_Dim /= 2) Then 
!!$       Write(*,*) Trim(Geom%Filename), 'is not a 2D File.'
!!$       Write(*,*) 'Aborting'
!!$       STOP
!!$    End If
    If (Associated(Res)) Then
       DeAllocate(Res)
    End If
    Allocate(Res(Geom%Num_Nodes))
    Call EXGNV(Geom%exoid, TS, Idx, Geom%Num_Nodes, Res%X, iErr)
    Call EXGNV(Geom%exoid, TS, Idx+1, Geom%Num_Nodes, Res%Y, iErr)
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Read_EXO_Result_Vect2D_Nodes

  Subroutine Read_EXO_Result_Mat2D_Nodes(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Type(Mat2D), Dimension(:), Pointer             :: Res

    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

!!$    If (Geom%Num_Dim /= 2) Then 
!!$       Write(*,*) Trim(Geom%Filename), 'is not a 2D File.'
!!$       Write(*,*) 'Aborting'
!!$       STOP
!!$    End If
    If (Associated(Res)) Then
       DeAllocate(Res)
    End If
    Allocate(Res(Geom%Num_Nodes))
    Call EXGNV(Geom%exoid, TS, Idx,   Geom%Num_Nodes, Res%XX, iErr)
    Call EXGNV(Geom%exoid, TS, Idx+1, Geom%Num_Nodes, Res%XY, iErr)
    Call EXGNV(Geom%exoid, TS, Idx+2, Geom%Num_Nodes, Res%YX, iErr)
    Call EXGNV(Geom%exoid, TS, Idx+3, Geom%Num_Nodes, Res%YY, iErr)
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Read_EXO_Result_Mat2D_Nodes

  Subroutine Read_EXO_Result_MatS2D_Nodes(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Type(MatS2D), Dimension(:), Pointer            :: Res

    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

!!$    If (Geom%Num_Dim /= 2) Then 
!!$       Write(*,*) Trim(Geom%Filename), 'is not a 2D File.'
!!$       Write(*,*) 'Aborting'
!!$       STOP
!!$    End If
    If (Associated(Res)) Then
       DeAllocate(Res)
    End If
    Allocate(Res(Geom%Num_Nodes))
    Call EXGNV(Geom%exoid, TS, Idx,   Geom%Num_Nodes, Res%XX, iErr)
    Call EXGNV(Geom%exoid, TS, Idx+1, Geom%Num_Nodes, Res%YY, iErr)
    Call EXGNV(Geom%exoid, TS, Idx+2, Geom%Num_Nodes, Res%XY, iErr)
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Read_EXO_Result_MatS2D_Nodes

  Subroutine Read_EXO_Result_Vect3D_Nodes(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Type(Vect3D), Dimension(:), Pointer            :: Res

    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

!!$    If (Geom%Num_Dim /= 3) Then 
!!$       Write(*,*) Trim(Geom%Filename), 'is not a 3D File.'
!!$       Write(*,*) 'Aborting'
!!$       STOP
!!$    End If
    If (Associated(Res)) Then
       DeAllocate(Res)
    End If
    Allocate(Res(Geom%Num_Nodes))
    Call EXGNV(Geom%exoid, TS, Idx,   Geom%Num_Nodes, Res%X, iErr)
    Call EXGNV(Geom%exoid, TS, Idx+1, Geom%Num_Nodes, Res%Y, iErr)
    Call EXGNV(Geom%exoid, TS, Idx+2, Geom%Num_Nodes, Res%Z, iErr)
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Read_EXO_Result_Vect3D_Nodes

  Subroutine Read_EXO_Result_Mat3D_Nodes(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Type(Mat3D), Dimension(:), Pointer             :: Res

    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

!!$    If (Geom%Num_Dim /= 3) Then 
!!$       Write(*,*) Trim(Geom%Filename), 'is not a 3D File.'
!!$       Write(*,*) 'Aborting'
!!$       STOP
!!$    End If
    If (Associated(Res)) Then
       DeAllocate(Res)
    End If
    Allocate(Res(Geom%Num_Nodes))
    Call EXGNV(Geom%exoid, TS, Idx,   Geom%Num_Nodes, Res%XX, iErr)
    Call EXGNV(Geom%exoid, TS, Idx+1, Geom%Num_Nodes, Res%XY, iErr)
    Call EXGNV(Geom%exoid, TS, Idx+2, Geom%Num_Nodes, Res%XZ, iErr)
    Call EXGNV(Geom%exoid, TS, Idx+3, Geom%Num_Nodes, Res%YX, iErr)
    Call EXGNV(Geom%exoid, TS, Idx+4, Geom%Num_Nodes, Res%YY, iErr)
    Call EXGNV(Geom%exoid, TS, Idx+5, Geom%Num_Nodes, Res%YZ, iErr)
    Call EXGNV(Geom%exoid, TS, Idx+6, Geom%Num_Nodes, Res%ZX, iErr)
    Call EXGNV(Geom%exoid, TS, Idx+7, Geom%Num_Nodes, Res%ZY, iErr)
    Call EXGNV(Geom%exoid, TS, Idx+8, Geom%Num_Nodes, Res%ZZ, iErr)
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Read_EXO_Result_Mat3D_Nodes

  Subroutine Read_EXO_Result_MatS3D_Nodes(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Type(MatS3D), Dimension(:), Pointer            :: Res

    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

!!$    If (Geom%Num_Dim /= 3) Then 
!!$       Write(*,*) Trim(Geom%Filename), 'is not a 3D File.'
!!$       Write(*,*) 'Aborting'
!!$       STOP
!!$    End If
    If (Associated(Res)) Then
       DeAllocate(Res)
    End If
    Allocate(Res(Geom%Num_Nodes))
    Call EXGNV(Geom%exoid, TS, Idx,   Geom%Num_Nodes, Res%XX, iErr)
    Call EXGNV(Geom%exoid, TS, Idx+1, Geom%Num_Nodes, Res%YY, iErr)
    Call EXGNV(Geom%exoid, TS, Idx+2, Geom%Num_Nodes, Res%ZZ, iErr)
    Call EXGNV(Geom%exoid, TS, Idx+3, Geom%Num_Nodes, Res%XY, iErr)
    Call EXGNV(Geom%exoid, TS, Idx+4, Geom%Num_Nodes, Res%YZ, iErr)
    Call EXGNV(Geom%exoid, TS, Idx+5, Geom%Num_Nodes, Res%XZ, iErr)
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Read_EXO_Result_MatS3D_Nodes

  Subroutine Read_EXO_Result_Scal_Elems(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Real(Kind = Kr), Dimension(:), Pointer         :: Res

    Integer                                        :: iBlk
    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    Do iBlk = 1, Geom%Num_Elem_Blks
       Call EXGEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID), iErr)
    End Do

    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Read_EXO_Result_Scal_Elems

  Subroutine Read_EXO_Result_Vect2D_Elems(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Type(Vect2D), Dimension(:), Pointer            :: Res

    Integer                                        :: iBlk
    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    Do iBlk = 1, Geom%Num_Elem_Blks
       Call EXGEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%X, iErr)
       Call EXGEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%Y, iErr)
    End Do
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Read_EXO_Result_Vect2D_Elems

  Subroutine Read_EXO_Result_Mat2D_Elems(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Type(Mat2D), Dimension(:), Pointer             :: Res

    Integer                                        :: iBlk
    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    Do iBlk = 1, Geom%Num_Elem_Blks
       Call EXGEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XX, iErr)
       Call EXGEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XY, iErr)
       Call EXGEV(Geom%exoid, TS, Idx+2, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YX, iErr)
       Call EXGEV(Geom%exoid, TS, Idx+3, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YY, iErr)
    End Do
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Read_EXO_Result_Mat2D_Elems

  Subroutine Read_EXO_Result_MatS2D_Elems(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Type(MatS2D), Dimension(:), Pointer            :: Res

    Integer                                        :: iBlk
    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    Do iBlk = 1, Geom%Num_Elem_Blks
       Call EXGEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XX, iErr)
       Call EXGEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YY, iErr)
       Call EXGEV(Geom%exoid, TS, Idx+2, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XY, iErr)
    End Do
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Read_EXO_Result_MatS2D_Elems

  Subroutine Read_EXO_Result_Vect3D_Elems(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Type(Vect3D), Dimension(:), Pointer            :: Res

    Integer                                        :: iBlk
    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    Do iBlk = 1, Geom%Num_Elem_Blks
       Call EXGEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%X, iErr)
       Call EXGEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%Y, iErr)
    End Do
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Read_EXO_Result_Vect3D_Elems

  Subroutine Read_EXO_Result_Mat3D_Elems(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Type(Mat3D), Dimension(:), Pointer             :: Res

    Integer                                        :: iBlk
    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    Do iBlk = 1, Geom%Num_Elem_Blks
       Call EXGEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XX, iErr)
       Call EXGEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XY, iErr)
       Call EXGEV(Geom%exoid, TS, Idx+2, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XZ, iErr)
       Call EXGEV(Geom%exoid, TS, Idx+3, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YX, iErr)
       Call EXGEV(Geom%exoid, TS, Idx+4, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YY, iErr)
       Call EXGEV(Geom%exoid, TS, Idx+5, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YZ, iErr)
       Call EXGEV(Geom%exoid, TS, Idx+6, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%ZX, iErr)
       Call EXGEV(Geom%exoid, TS, Idx+7, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%ZY, iErr)
       Call EXGEV(Geom%exoid, TS, Idx+8, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%ZZ, iErr)
    End Do
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Read_EXO_Result_Mat3D_Elems

  Subroutine Read_EXO_Result_MatS3D_Elems(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Type(MatS3D), Dimension(:), Pointer            :: Res

    Integer                                        :: iBlk
    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    Do iBlk = 1, Geom%Num_Elem_Blks
       Call EXGEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XX, iErr)
       Call EXGEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YY, iErr)
       Call EXGEV(Geom%exoid, TS, Idx+2, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%ZZ, iErr)
       Call EXGEV(Geom%exoid, TS, Idx+3, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XY, iErr)
       Call EXGEV(Geom%exoid, TS, Idx+4, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YZ, iErr)
       Call EXGEV(Geom%exoid, TS, Idx+5, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XZ, iErr)
    End Do
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Read_EXO_Result_MatS3D_Elems

!!! WRITE
  Subroutine Write_EXO_Result_Global(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Real(Kind = Kr)                                :: Res

    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Real(Kind = Kr), Dimension(:), Pointer         :: Tmp_Res
    Integer                                        :: Num_Vars, Num_TS
    Real(Kind = Kr)                                :: fDum
    Character                                      :: cDum

    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)
    ! Get the number of global variables stored in the database    
    Call EXGVP(Geom%exoid, 'G', Num_Vars, iErr)
    Allocate(Tmp_Res(Num_Vars))
    Tmp_Res = 0.0_Kr

    ! Read All the global variables at time step TS into Tmp_Res
    ! Modify Tmp_Res(Idx) and write everything back...
    Call EXGGV(Geom%exoid, TS, Num_Vars, Tmp_Res, iErr)
    Tmp_Res(Idx) = Res

    Call EXPGV(Geom%exoid, TS, Num_Vars, Tmp_Res, iErr)
    DeAllocate (Tmp_Res)
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_EXO_Result_Global


!!$  Subroutine Write_EXO_Result_Scal_Nodes(Geom, Idx, TS, Res)
!!$    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
!!$    Integer                                        :: Idx
!!$    Integer                                        :: TS
!!$    Real(Kind = Kr), Dimension(:), Pointer         :: Res
!!$
!!$    Integer                                        :: iErr
!!$    Real(Kind = Kr)                                :: Vers
!!$    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
!!$         & ierr)
!!$
!!$    Call EXPNV(Geom%exoid, TS, Idx, Geom%Num_Nodes, Res, iErr)
!!$    
!!$    Call EXCLOS(Geom%exoid, iErr)
!!$    Geom%exoid = 0
!!$  End Subroutine Write_EXO_Result_Scal_Nodes
!!$
  Subroutine Write_EXO_Result_Ptr_Nodes(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Real(Kind = Kr), Dimension(:), Pointer         :: Res

    Real(Kind = Kr), Dimension(:), Pointer         :: Tmp_Res
    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Integer                                        :: Num_Rec, iRec

    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)
    Select Case (Geom%Numbering)
    Case(Numbering_PerCoord)
       If ( Mod(Size(Res), Geom%Num_Nodes) == 0) Then
          Num_Rec = Size(Res) / Geom%Num_Nodes
          Do iRec = 0, Num_Rec-1
             Call EXPNV(Geom%exoid, TS, Idx + iRec, Geom%Num_Nodes,           &
                  & Res(1+iRec * Geom%Num_Nodes:(iRec+1) * Geom%Num_Nodes),   &
                  &iErr)
          End Do
       Else
          Write(*,*) 'Error in Write_EXO_Result_Scal_Ptr_Nodes'
          Write(*,*) 'Can''t make any sense of the dimension of the arguments:'
          Write(*,*) 'Size(Res), Geom%Num_Dim ', Size(Res), Geom%Num_Dim
          Write(*,*) 'Not saving anything in', Trim(Geom%FileName)
       End If
    Case(Numbering_PerNodes)
       If ( Mod(Size(Res), Geom%Num_Nodes) == 0) Then
          Allocate(Tmp_Res(Geom%Num_Nodes))
          Num_Rec = Size(Res) / Geom%Num_Nodes
          Do iRec = 1, Num_Rec
             Tmp_Res = Res(iRec:Num_Rec * Geom%Num_Nodes:Num_Rec)
             Call EXPNV(Geom%exoid, TS, Idx + iRec - 1, Geom%Num_Nodes,       &
                  & Tmp_Res, iErr)
          End Do
          DeAllocate(Tmp_Res)
       Else
          Write(*,*) 'Error in Write_EXO_Result_Scal_Ptr_Nodes'
          Write(*,*) 'Can''t make any sense of the dimension of the arguments:'
          Write(*,*) 'Size(Res), Geom%Num_Dim ', Size(Res), Geom%Num_Dim
          Write(*,*) 'Not saving anything in', Trim(Geom%FileName)
       End If
    Case Default
       Write(*,*) 'ERROR:Write_EXO_Result_Ptr_Nodes Numbering scheme unknown',&
            & Geom%Numbering
       STOP
    End Select

!!$    If ( Mod(Size(Res), Geom%Num_Nodes) == 0) Then
!!$       Num_Rec = Size(Res) / Geom%Num_Nodes
!!$       Do iRec = 0, Num_Rec-1
!!$          Call EXPNV(Geom%exoid, TS, Idx + iRec, Geom%Num_Nodes,              &
!!$               & Res(1+iRec * Geom%Num_Nodes:(iRec+1) * Geom%Num_Nodes), iErr)
!!$       End Do
!!$    Else
!!$       Write(*,*) 'Error in Write_EXO_Result_Scal_Ptr_Nodes'
!!$       Write(*,*) 'Can''t make any sense of the dimension of the arguments:'
!!$       Write(*,*) 'Size(Res), Geom%Num_Dim ', Size(Res), Geom%Num_Dim
!!$       Write(*,*) 'Not saving anything in', Trim(Geom%FileName)
!!$    End If
!!$       
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_EXO_Result_Ptr_Nodes

  Subroutine Write_EXO_Result_Vect2D_Nodes(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Type(Vect2D), Dimension(:), Pointer            :: Res

    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

!!$    If (Geom%Num_Dim /= 2) Then 
!!$       Write(*,*) Trim(Geom%Filename), 'is not a 2D File.'
!!$       Write(*,*) 'Aborting'
!!$       STOP
!!$    End If
    Call EXPNV(Geom%exoid, TS, Idx, Geom%Num_Nodes, Res%X, iErr)
    Call EXPNV(Geom%exoid, TS, Idx+1, Geom%Num_Nodes, Res%Y, iErr)
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_EXO_Result_Vect2D_Nodes

  Subroutine Write_EXO_Result_Mat2D_Nodes(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Type(Mat2D), Dimension(:), Pointer            :: Res

    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

!!$    If (Geom%Num_Dim /= 2) Then 
!!$       Write(*,*) Trim(Geom%Filename), 'is not a 2D File.'
!!$       Write(*,*) 'Aborting'
!!$       STOP
!!$    End If
    Call EXPNV(Geom%exoid, TS, Idx,   Geom%Num_Nodes, Res%XX, iErr)
    Call EXPNV(Geom%exoid, TS, Idx+1, Geom%Num_Nodes, Res%XY, iErr)
    Call EXPNV(Geom%exoid, TS, Idx+2, Geom%Num_Nodes, Res%YX, iErr)
    Call EXPNV(Geom%exoid, TS, Idx+3, Geom%Num_Nodes, Res%YY, iErr)
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_EXO_Result_Mat2D_Nodes

  Subroutine Write_EXO_Result_MatS2D_Nodes(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Type(MatS2D), Dimension(:), Pointer            :: Res

    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

!!$    If (Geom%Num_Dim /= 2) Then 
!!$       Write(*,*) Trim(Geom%Filename), 'is not a 2D File.'
!!$       Write(*,*) 'Aborting'
!!$       STOP
!!$    End If
    Call EXPNV(Geom%exoid, TS, Idx,   Geom%Num_Nodes, Res%XX, iErr)
    Call EXPNV(Geom%exoid, TS, Idx+1, Geom%Num_Nodes, Res%YY, iErr)
    Call EXPNV(Geom%exoid, TS, Idx+2, Geom%Num_Nodes, Res%XY, iErr)
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_EXO_Result_MatS2D_Nodes

  Subroutine Write_EXO_Result_Vect3D_Nodes(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Type(Vect3D), Dimension(:), Pointer            :: Res

    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

!!$    If (Geom%Num_Dim /= 3) Then 
!!$       Write(*,*) Trim(Geom%Filename), 'is not a 3D File.'
!!$       Write(*,*) 'Aborting'
!!$       STOP
!!$    End If
    Call EXPNV(Geom%exoid, TS, Idx,   Geom%Num_Nodes, Res%X, iErr)
    Call EXPNV(Geom%exoid, TS, Idx+1, Geom%Num_Nodes, Res%Y, iErr)
    Call EXPNV(Geom%exoid, TS, Idx+2, Geom%Num_Nodes, Res%Z, iErr)
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_EXO_Result_Vect3D_Nodes

  Subroutine Write_EXO_Result_Mat3D_Nodes(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Type(Mat3D), Dimension(:), Pointer             :: Res

    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

!!$    If (Geom%Num_Dim /= 3) Then 
!!$       Write(*,*) Trim(Geom%Filename), 'is not a 3D File.'
!!$       Write(*,*) 'Aborting'
!!$       STOP
!!$    End If
    Call EXPNV(Geom%exoid, TS, Idx,   Geom%Num_Nodes, Res%XX, iErr)
    Call EXPNV(Geom%exoid, TS, Idx+1, Geom%Num_Nodes, Res%XY, iErr)
    Call EXPNV(Geom%exoid, TS, Idx+2, Geom%Num_Nodes, Res%XZ, iErr)
    Call EXPNV(Geom%exoid, TS, Idx+3, Geom%Num_Nodes, Res%YX, iErr)
    Call EXPNV(Geom%exoid, TS, Idx+4, Geom%Num_Nodes, Res%YY, iErr)
    Call EXPNV(Geom%exoid, TS, Idx+5, Geom%Num_Nodes, Res%YZ, iErr)
    Call EXPNV(Geom%exoid, TS, Idx+6, Geom%Num_Nodes, Res%ZX, iErr)
    Call EXPNV(Geom%exoid, TS, Idx+7, Geom%Num_Nodes, Res%ZY, iErr)
    Call EXPNV(Geom%exoid, TS, Idx+8, Geom%Num_Nodes, Res%ZZ, iErr)

    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_EXO_Result_Mat3D_Nodes

  Subroutine Write_EXO_Result_MatS3D_Nodes(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Type(MatS3D), Dimension(:), Pointer            :: Res

    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

!!$    If (Geom%Num_Dim /= 3) Then 
!!$       Write(*,*) Trim(Geom%Filename), 'is not a 3D File.'
!!$       Write(*,*) 'Aborting'
!!$       STOP
!!$    End If
    Call EXPNV(Geom%exoid, TS, Idx,   Geom%Num_Nodes, Res%XX, iErr)
    Call EXPNV(Geom%exoid, TS, Idx+1, Geom%Num_Nodes, Res%YY, iErr)
    Call EXPNV(Geom%exoid, TS, Idx+2, Geom%Num_Nodes, Res%ZZ, iErr)
    Call EXPNV(Geom%exoid, TS, Idx+3, Geom%Num_Nodes, Res%XY, iErr)
    Call EXPNV(Geom%exoid, TS, Idx+4, Geom%Num_Nodes, Res%YZ, iErr)
    Call EXPNV(Geom%exoid, TS, Idx+5, Geom%Num_Nodes, Res%XZ, iErr)

    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_EXO_Result_MatS3D_Nodes

  Subroutine Write_EXO_Result_Scal_Elems(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Real(Kind = Kr), Dimension(:), Pointer         :: Res

    Real(Kind = Kr), Dimension(:), Pointer         :: Tmp_Res
    Integer                                        :: iBlk, iE
    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    Do iBlk = 1, Geom%Num_Elem_Blks
       Allocate (Tmp_Res(Geom%Elem_Blk(iBlk)%Num_Elems))
       Do iE = 1, Geom%Elem_Blk(iBlk)%Num_Elems
          Tmp_Res(iE) = Res(Geom%Elem_Blk(iBlk)%Elem_ID(iE))
       End Do

       Call EXPEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
            & Tmp_Res, iErr)
       DeAllocate(Tmp_Res)
    End Do

    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_EXO_Result_Scal_Elems

  Subroutine Write_EXO_Result_Vect2D_Elems(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Type(Vect2D), Dimension(:), Pointer            :: Res

    Integer                                        :: iBlk
    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    Do iBlk = 1, Geom%Num_Elem_Blks
       Call EXPEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%X, iErr)
       Call EXPEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%Y, iErr)
    End Do
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_EXO_Result_Vect2D_Elems

  Subroutine Write_EXO_Result_Mat2D_Elems(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Type(Mat2D), Dimension(:), Pointer             :: Res

    Integer                                        :: iBlk
    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    Do iBlk = 1, Geom%Num_Elem_Blks
       Call EXPEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XX, iErr)
       Call EXPEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XY, iErr)
       Call EXPEV(Geom%exoid, TS, Idx+2, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YX, iErr)
       Call EXPEV(Geom%exoid, TS, Idx+3, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YY, iErr)
    End Do
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_EXO_Result_Mat2D_Elems

  Subroutine Write_EXO_Result_MatS2D_Elems(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Type(MatS2D), Dimension(:), Pointer            :: Res

    Integer                                        :: iBlk
    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    Do iBlk = 1, Geom%Num_Elem_Blks
       Call EXPEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XX, iErr)
       Call EXPEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YY, iErr)
       Call EXPEV(Geom%exoid, TS, Idx+2, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XY, iErr)
    End Do
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_EXO_Result_MatS2D_Elems

  Subroutine Write_EXO_Result_Vect3D_Elems(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Type(Vect3D), Dimension(:), Pointer            :: Res

    Type(Vect3D), DImension(:), Pointer            :: Tmp_Res
    Integer                                        :: iBlk, iE
    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    Do iBlk = 1, Geom%Num_Elem_Blks
       Allocate (Tmp_Res(Geom%Elem_Blk(iBlk)%Num_Elems))
       Do iE = 1, Geom%Elem_Blk(iBlk)%Num_Elems
          Tmp_Res(iE) = Res(Geom%Elem_Blk(iBlk)%Elem_ID(iE))
       End Do
       Call EXPEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
            & Tmp_Res%X, iErr)
       Call EXPEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Tmp_Res%Y, iErr)
       Call EXPEV(Geom%exoid, TS, Idx+2, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Tmp_Res%Z, iErr)
       DeAllocate(Tmp_Res)
    End Do
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_EXO_Result_Vect3D_Elems

  Subroutine Write_EXO_Result_Mat3D_Elems(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Type(Mat3D), Dimension(:), Pointer             :: Res

    Integer                                        :: iBlk
    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    Do iBlk = 1, Geom%Num_Elem_Blks
       Call EXPEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XX, iErr)
       Call EXPEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XY, iErr)
       Call EXPEV(Geom%exoid, TS, Idx+2, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XZ, iErr)
       Call EXPEV(Geom%exoid, TS, Idx+3, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YX, iErr)
       Call EXPEV(Geom%exoid, TS, Idx+4, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YY, iErr)
       Call EXPEV(Geom%exoid, TS, Idx+5, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YZ, iErr)
       Call EXPEV(Geom%exoid, TS, Idx+6, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%ZX, iErr)
       Call EXPEV(Geom%exoid, TS, Idx+7, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%ZY, iErr)
       Call EXPEV(Geom%exoid, TS, Idx+8, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%ZZ, iErr)
    End Do
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_EXO_Result_Mat3D_Elems

  Subroutine Write_EXO_Result_MatS3D_Elems(Geom, Idx, TS, Res)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Integer                                        :: Idx
    Integer                                        :: TS
    Type(MatS3D), Dimension(:), Pointer            :: Res

    Type(MatS3D), Dimension(:), Pointer            :: Tmp_Res
    Integer                                        :: iBlk, iE
    Integer                                        :: iErr
    Real(Kind = Kr)                                :: Vers
    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers,   &
         & ierr)

    Do iBlk = 1, Geom%Num_Elem_Blks
       Allocate (Tmp_Res(Geom%Elem_Blk(iBlk)%Num_Elems))
       Do iE = 1, Geom%Elem_Blk(iBlk)%Num_Elems
          Tmp_Res(iE) = Res(Geom%Elem_Blk(iBlk)%Elem_ID(iE))
       End Do
       Call EXPEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
            & Tmp_Res%XX, iErr)
       Call EXPEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Tmp_Res%YY, iErr)
       Call EXPEV(Geom%exoid, TS, Idx+2, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Tmp_Res%ZZ, iErr)
       Call EXPEV(Geom%exoid, TS, Idx+3, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Tmp_Res%XY, iErr)
       Call EXPEV(Geom%exoid, TS, Idx+4, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Tmp_Res%YZ, iErr)
       Call EXPEV(Geom%exoid, TS, Idx+5, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
            & Tmp_Res%XZ, iErr)
!!$
!!$       Call EXPEV(Geom%exoid, TS, Idx, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems,   &
!!$            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XX, iErr)
!!$       Call EXPEV(Geom%exoid, TS, Idx+1, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!$            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YY, iErr)
!!$       Call EXPEV(Geom%exoid, TS, Idx+2, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!$            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%ZZ, iErr)
!!$       Call EXPEV(Geom%exoid, TS, Idx+3, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!$            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XY, iErr)
!!$       Call EXPEV(Geom%exoid, TS, Idx+4, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!$            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%YZ, iErr)
!!$       Call EXPEV(Geom%exoid, TS, Idx+5, iBlk, Geom%Elem_Blk(iBlk)%Num_Elems, &
!!$            & Res(Geom%Elem_Blk(iBlk)%Elem_ID)%XZ, iErr)
       DeAllocate(Tmp_Res)
    End Do
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_EXO_Result_MatS3D_Elems


End Module m_MEF_EXO

