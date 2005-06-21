Module m_MEF_PU
!! Generating PU functions on each subdomains with the uniform size of overlap
  Use m_MEF_Types
  Use m_MEF_EXO
  IMPLICIT NONE

!  Integer, Parameter                                :: ovlp = 1



Contains

  Subroutine Gen_PU_2D_Scal(Geom, Elem_db, Node_db, ovlp, PU_function, Ovlp_Elem_Blks)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Type (Element2D_Scal), Dimension(:), Pointer        :: Elem_db
    Type (Node2D), Dimension(:), Pointer           :: Node_db
    Integer                                        :: ovlp
    Integer, Dimension(:,:), Pointer               :: Count_function
    Real, Dimension(:,:), Pointer               :: PU_function
    Integer, Dimension(:,:), Pointer               :: Count_ID
    Integer, Dimension(:), Pointer                 :: Total
    Integer, Dimension(:), Pointer                 :: Num_Elem_Ovlp_Elem_Blks
    Type(Elem_Blk_Info), Dimension(:), Pointer     :: Ovlp_Elem_Blks
    Integer                                        :: iBlock, Offset, iE, iO, iC
    Integer                                        :: IEx, iB
    Real(Kind = Kr)                                :: Vers
    Integer                                        :: iErr
    Integer                                        :: Nb_DoF

    Allocate (Count_function(Geom%Num_Elem_Blks,Geom%Num_Nodes))
    Allocate (PU_function(Geom%Num_Elem_Blks,Geom%Num_Nodes))
    Allocate (Count_ID(Geom%Num_Elem_Blks,Geom%Num_Elems))
    Allocate (Total(Geom%Num_Nodes))
    Allocate (Ovlp_Elem_Blks(Geom%Num_Elem_Blks))
    Allocate (Num_Elem_Ovlp_Elem_Blks(Geom%Num_Elem_Blks))
!! Initialize PU for nonoverlapping subdomains
    Offset=0
    Count_function = 0
    Count_ID = 0
    Total = 0
  Do iBlock = 1, Geom%Num_Elem_Blks
     Do iE = 1, Geom%Num_Elems  
        If (Elem_db(iE)%ID_EL == iBlock) Then
          Count_function(iBlock,Elem_db(iE+Offset)%ID_DoF) =1
        End If
      End Do
      Count_ID(iBlock,Geom%Elem_Blk(iBlock)%Elem_ID) = 1
    End Do

!! Extended Count_function with overlap 
   Do iB = 1, Geom%Num_Elem_Blks
   PU_function(iB,:)=Count_function(iB,:)
   End Do 
   If (ovlp >= 1) then
   Do iO = 1, ovlp
   Do iB = 1, Geom%Num_Elem_Blks

!! First check the one layer for extension; We need to restrict one layer
    Do iE = 1, Geom%Num_Elems
    If (any(Count_function(iB,Elem_db(iE)%ID_DoF)> 0) .and. Count_ID(iB,iE) == 0) Then
       Count_ID(iB,iE) = -1
    End If
    End Do
!! Change the count function for the extension
    Do iE = 1, Geom%Num_Elems 
    If (Count_ID(iB,iE) == -1) Then
       Count_function(iB,Elem_db(iE)%ID_DoF) = 1
       Count_ID(iB,iE) = 1
    End If
    End Do
      PU_function(iB,:) = Count_function(iB,:) + PU_function(iB,:)
    End Do
  End Do
  End If

  Do iB = 1, Geom%Num_Elem_Blks 
    Total(:) = Total(:) + PU_function(iB,:)
   Num_Elem_Ovlp_Elem_Blks(iB) = count(Count_ID(iB,:) ==1)
  End Do

!! Generate PU function Overlapping Element Blks
   Do iB = 1, Geom%Num_Elem_Blks
    Do iO = 1, Geom%Num_Nodes
    If (Total(iO) /= 0) Then
      PU_function(iB,iO) = PU_function(iB,iO)/Total(iO)
    Else
      PU_function(iB,iO) = 0
    End If
    End Do
   Ovlp_Elem_Blks(iB)%ID=iB
   Ovlp_Elem_Blks(iB)%Num_Elems =Num_Elem_Ovlp_Elem_Blks(iB)
   Allocate (Ovlp_Elem_Blks(iB)%Elem_ID(Ovlp_Elem_Blks(iB)%Num_Elems))
  End Do
   iEx =0
   Do iB = 1, Geom%Num_Elem_Blks
     iEx =0
    Do iE = 1, Geom%Num_Elems
     If (Count_ID(iB,iE) == 1) Then
       iEx = iEx+1
       Ovlp_Elem_Blks(iB)%Elem_ID(iEx)=iE
     End If 
    End Do
     iEx =0
   End Do

  End Subroutine Gen_PU_2D_Scal

  Subroutine Show_Ovlp_Info(Geom, PU_function, Ovlp_Elem_Blks)
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
    Real, Dimension(:,:), Pointer               :: PU_function
    Type(Elem_Blk_Info), Dimension(:), Pointer     :: Ovlp_Elem_Blks

    Integer                                        :: i
    Write(*,*)
    Write(*,100)
    Do i = 1, Geom%Num_Nodes
      Print *, PU_function(1,i), PU_function(2,i), PU_function(3,i)
    End Do


    Write(*,*)
    Write(*,200)
    Write(*,201) Geom%num_elem_blks
    Do i = 1, Geom%Num_Elem_blks
       Write(*,203) Ovlp_Elem_Blks(i)%ID, Ovlp_Elem_Blks(i)%Num_Elems
       Write(*,206, advance = 'no') Ovlp_Elem_Blks(i)%ID
       Write(*,*) Ovlp_Elem_Blks(i)%Elem_ID
    End Do

100 Format('*** PU Functions ***')
200 Format('*** OVERLAPPING ELEMENT BLOCKS ***')
201 Format('    Number of blocks ================ ', I4)
203 Format('    Block ', I3, ' Number of elements ==== ', I4)
206 Format('    Block ', I3, ' IDs: ')


  End Subroutine Show_Ovlp_Info



End Module m_MEF_PU

