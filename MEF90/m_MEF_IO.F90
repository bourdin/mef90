module m_MEF_IO
! Blaise Bourdin, 1996-1998
! Merci de faire parvenir toutes remarques et bugs 
! eventuels aux adresses suivantes :
! bourdin@lpmtm.univ-paris13.fr
! bourdin@mat.dtu.dk
!
! En cas de modifications, les noms de  modules et de 
! fichiers _DOIVENT_ etre renommes
!          ^^^^^^^^^
!
  Use m_MEF_Types
  IMPLICIT NONE

  Interface Read_Mesh_Tri_P1
     Module Procedure Read_Mesh_2D_Scal_Tri_P1, &
          & Read_Mesh_2D_Elast_Tri_P1
  End Interface

  Interface Write_Mesh_Tri_P1
     Module Procedure Write_Mesh_2D_Scal_Tri_P1, &
          & Write_Mesh_2D_Elast_Tri_P1
  End Interface

  Interface Read_bb
     Module Procedure Read_bb_1, Read_bb_2, Read_bb_3, Read_bb_Ptr
  End Interface

  Interface Write_bb
     Module Procedure Write_bb_Vect2D, Write_bb_Vect3D, &
          & Write_bb_Mat2D, Write_bb_Mat3D,             &
          & Write_bb_MatS2D, Write_bb_MatS3D,           &
          & Write_bb_Ptr
  End Interface


Contains
!$$ ============ MESH READING ============

  Subroutine Read_Mesh_2DA_Tri_P1(Mesh_fmt, Mesh_Str,Elem_bd, Nodes_bd)
    Character(len = *)                          :: Mesh_Fmt
    Character(len = *)                          :: Mesh_Str
    Type (Element2DA), Dimension(:), Pointer    :: Elem_bd
    Type (Node2D), Dimension(:), Pointer        :: Nodes_bd

    Integer                                     :: NE, NS
    Integer                                     :: iStat, iS, jS, iE, jE

    Open (Unit=F_In, file=Mesh_Str, status='old', err=100, iostat=iStat)
    Rewind(F_In)
    Msh_fmt : select case (Mesh_fmt)
    Case ('msh', 'MSH')
       Read(f_In,*) NS, NE
       Allocate (Nodes_bd(NS))
       Allocate (Elem_bd(NE))
       ! On dimensionne le Nb de DoF des elements de elem_bd : 
       !     Tri P1 => 3 DoF       
       Do iE=1, NE
          Elem_bd(iE)%NB_DoF = 3
          Allocate(Elem_bd(iE)%ID_DoF(3))
       EndDo

       Do iS=1, NS
          Read(f_In,*) Nodes_bd(iS)%Coord%X, Nodes_bd(iS)%Coord%Y, &
               &            Nodes_bd(iS)%ID
       EndDo
       Do iE=1, NE
          Read(f_IN,*) Elem_bd(iE)%ID_DoF(:), ELem_bd(iE)%ID_EL
       EndDo

    Case ('amdba', 'AMDBA')
       Read(f_In,*) NS, NE
       Allocate (Nodes_bd(NS))
       Allocate (Elem_bd(NE))
       ! On dimensionne le Nb de DoF des elements de elem_bd : 
       !     Tri P1 => 3 DoF       
       Do iE=1, NE
          Elem_bd(iE)%NB_DoF = 3
          Allocate(Elem_bd(iE)%ID_DoF(3))
       EndDo

       Do iS=1, NS
          Read(f_In,*) jS, Nodes_bd(jS)%Coord%X, Nodes_bd(jS)%Coord%Y, &
               &            Nodes_bd(jS)%ID
       EndDo
       Do iE=1, NE
          Read(f_IN,*) jE, Elem_bd(jE)%ID_DoF(:), ELem_bd(jE)%ID_EL
       EndDo

    End Select Msh_Fmt
    Close (f_In)

    Return

100 Print*, 'Erreur a l''ouverture du fichier '//Mesh_Str, iStat
    Stop
  End Subroutine Read_Mesh_2DA_Tri_P1

  Subroutine Read_Mesh_2D_Scal_Tri_P1(Mesh_fmt, Mesh_Str,Elem_bd, Nodes_bd)
    Character(*)                                        :: Mesh_Fmt
    Character(*)                                        :: Mesh_Str
    Type (Element2D_Scal), Dimension(:), Pointer       :: Elem_bd
    Type (Node2D), Dimension(:), Pointer                :: Nodes_bd

    Integer                                             :: NE, NS, iStat
    Integer                                             :: iS, jS, iE, jE

    Open (Unit=F_In, file=Mesh_Str, status='old', err=100, iostat=iStat)
    Rewind(F_In)
    Msh_fmt : select case (Mesh_fmt)
    Case ('msh', 'MSH')
       Read(f_In,*) NS, NE
       Allocate (Nodes_bd(NS))
       Allocate (Elem_bd(NE))
       ! On dimensionne le Nb de DoF des elements de elem_bd : 
       !     Tri P1 => 3 DoF       
       Do iE=1, NE
          Elem_bd(iE)%NB_DoF = 3
          Allocate(Elem_bd(iE)%ID_DoF(3))
       EndDo

       Do iS=1, NS
          Read(f_In,*) Nodes_bd(iS)%Coord%X, Nodes_bd(iS)%Coord%Y, &
               &            Nodes_bd(iS)%ID
       EndDo
       Do iE=1, NE
          Read(f_IN,*) Elem_bd(iE)%ID_DoF(:) 
       EndDo
       Elem_bd(:)%ID_EL = 0

    Case ('amdba', 'AMDBA')
       Read(f_In,*) NS, NE
       Allocate (Nodes_bd(NS))
       Allocate (Elem_bd(NE))
       ! On dimensionne le Nb de DoF des elements de elem_bd : 
       !     Tri P1 => 3 DoF       
       Do iE=1, NE
          Elem_bd(iE)%NB_DoF = 3
          Allocate(Elem_bd(iE)%ID_DoF(3))
       EndDo

       Do iS=1, NS
          Read(f_In,*) jS, Nodes_bd(jS)%Coord%X, Nodes_bd(jS)%Coord%Y, &
               &            Nodes_bd(jS)%ID
       EndDo
       Do iE=1, NE
          Read(f_IN,*) jE, Elem_bd(jE)%ID_DoF(:), Elem_bd(jE)%ID_EL
       EndDo

       
    End Select Msh_Fmt
    Close (f_In)

    Return

100 Print*, 'Erreur a l''ouverture du fichier ',Mesh_Str, iStat
    Stop
  End Subroutine Read_Mesh_2D_Scal_Tri_P1

  Subroutine Read_Mesh_2D_Elast_Tri_P1(Mesh_fmt, Mesh_Str,Elem_bd, Nodes_bd)
    Character(len=*)                                        :: Mesh_Fmt
    Character(len=*)                                        :: Mesh_Str
    Type (Element2D_Elast), Dimension(:), Pointer       :: Elem_bd
    Type (Node2D), Dimension(:), Pointer                :: Nodes_bd
    Integer, Dimension(3)                               :: Tmp_ID

    Integer                                             :: NE, NS, iDoF, iStat
    Integer                                             :: iS, jS, iE, jE


    Open (Unit=F_In, file=Mesh_Str, status='old', err=100, iostat=iStat)
    Rewind(F_In)
    Msh_fmt : select case (Mesh_fmt)
    Case ('msh', 'MSH')
       Read(f_In,*) NS, NE
       NS = 2*NS
       Allocate (Nodes_bd(NS))
       Allocate (Elem_bd(NE))
       ! On dimensionne le Nb de DoF des elements de elem_bd : 
       !     Tri P1 => 6 DoF       
       Do iE=1, NE
          Elem_bd(iE)%NB_DoF = 6
          Allocate(Elem_bd(iE)%ID_DoF(6))
       EndDo

       Do iS=1, NS,2
          Read(f_In,*) Nodes_bd(iS)%Coord%X, Nodes_bd(iS)%Coord%Y, &
               &            Nodes_bd(iS)%ID
          Nodes_bd(iS+1) = Nodes_bd(iS)
!!$          Nodes_bd(iS+1)%ID =  mod(Nodes_bd(iS)%ID, 10_Ki)
!!$          Nodes_bd(is)%ID = Nodes_bd(iS)%ID / 10_Ki
       EndDo
       Do iE=1, NE
          Read(f_IN,*) Tmp_ID, Elem_bd(iE)%ID_EL
          Do iDoF = 1,3
             Elem_bd(iE)%ID_DoF(2*iDoF-1) = 2 * Tmp_ID(iDoF)-1
             Elem_bd(iE)%ID_DoF(2*iDoF) = 2 * Tmp_ID(iDoF)
          EndDo
       EndDo

    Case ('amdba', 'AMDBA')
       Read(f_In,*) NS, NE
       NS = 2*NS
       Allocate (Nodes_bd(NS))
       Allocate (Elem_bd(NE))
       ! On dimensionne le Nb de DoF des elements de elem_bd : 
       !     Tri P1 =>  6 DoF      
       Do iE=1, NE
          Elem_bd(iE)%NB_DoF = 6
          Allocate(Elem_bd(iE)%ID_DoF(6))
       EndDo

       Do iS=1, NS, 2
          Read(f_In,*) jS, Nodes_bd(iS)%Coord%X, Nodes_bd(iS)%Coord%Y, &
               &            Nodes_bd(iS)%ID
          Nodes_bd(iS+1) = Nodes_bd(iS)
!!$          Nodes_bd(iS+1)%ID =  mod(Nodes_bd(iS)%ID, 10_Ki)
!!$          Nodes_bd(is)%ID = Nodes_bd(iS)%ID / 10_Ki
       EndDo
       Do iE=1, NE
          Read(f_IN,*) jE, Tmp_ID, Elem_bd(iE)%ID_EL
          Do iDoF = 1,3
             Elem_bd(iE)%ID_DoF(2*iDoF-1) = 2 * Tmp_ID(iDoF) -1
             Elem_bd(iE)%ID_DoF(2*iDoF) = 2 * Tmp_ID(iDoF)
          EndDo
       EndDo

    End Select Msh_Fmt
    Close (f_In)

    Return

100 Print*, 'Erreur a l''ouverture du fichier '//Mesh_Str, iStat
    Stop
  End Subroutine Read_Mesh_2D_Elast_Tri_P1

  Subroutine Read_Mesh_1D_P1(Mesh_fmt, Mesh_Str,Elem_bd, Nodes_bd)
    Character(*)                                :: Mesh_Fmt
    Character(*)                                :: Mesh_Str
    Type (Element1D), Dimension(:), Pointer     :: Elem_bd
    Type (Node1D), Dimension(:), Pointer        :: Nodes_bd

    Integer                                     :: NE, NS
    Integer                                     :: iStat, iS, jS, iE, jE

    Open (Unit=F_In, file=Mesh_Str, status='old', err=100, iostat=iStat)
    Rewind(F_In)
    Msh_fmt : select case (Mesh_fmt)
    Case ('msh', 'MSH')
       Read(f_In,*) NS, NE
       Allocate (Nodes_bd(NS))
       Allocate (Elem_bd(NE))
       ! On dimensionne le Nb de DoF des elements de elem_bd : 
       !     Tri P1 => 2 DoF       
       Do iE=1, NE
          Elem_bd(iE)%NB_DoF = 2
          Allocate(Elem_bd(iE)%ID_DoF(2))
       EndDo

       Do iS=1, NS
          Read(f_In,*) Nodes_bd(iS)%Coord, Nodes_bd(iS)%ID
       EndDo
       Do iE=1, NE
          Read(f_IN,*) Elem_bd(iE)%ID_DoF(:), Elem_bd(iE)%ID_EL
       EndDo

    Case ('amdba', 'AMDBA')
       Read(f_In,*) NS, NE
       Allocate (Nodes_bd(NS))
       Allocate (Elem_bd(NE))
       ! On dimensionne le Nb de DoF des elements de elem_bd : 
       !     Tri P1 => 2 DoF       
       Do iE=1, NE
          Elem_bd(iE)%NB_DoF = 2
          Allocate(Elem_bd(iE)%ID_DoF(2))
       EndDo

       Do iS=1, NS
          Read(f_In,*) jS, Nodes_bd(jS)%Coord, Nodes_bd(jS)%ID
       EndDo
       Do iE=1, NE
          Read(f_IN,*) jE, Elem_bd(jE)%ID_DoF(:), Elem_bd(jE)%ID_EL
       EndDo

    End Select Msh_Fmt
    Close (f_In)

    Return

100 Print*, 'Erreur a l''ouverture du fichier '//Mesh_Str, iStat
    Stop
  End Subroutine Read_Mesh_1D_P1

!$$ ============ MESH WRITE ============
  Subroutine Write_Mesh_2D_Scal_Tri_P1(Mesh_fmt, Mesh_Str,Elem_bd, Nodes_bd)
    Character(*)                                        :: Mesh_Fmt
    Character(*)                                        :: Mesh_Str
    Type (Element2D_Scal), Dimension(:), Pointer       :: Elem_bd
    Type (Node2D), Dimension(:), Pointer                :: Nodes_bd

    Integer                                             :: NE, NS, iDoF, iStat
    Integer                                             :: iS, iE

    Open (Unit=F_Out, file=Mesh_Str, status='unknown', err=100, iostat=iStat)
    Rewind(F_Out)
    Msh_fmt : select case (Mesh_fmt)
    Case ('msh', 'MSH')
       NS = size(Nodes_bd)
       NE = size(Elem_bd)
       Write(F_OUT,*) NS, NE
       Do iS=1, NS
          Write(f_OUT,*) Nodes_bd(iS)%Coord%X, Nodes_bd(iS)%Coord%Y, &
               &            Nodes_bd(iS)%ID
       EndDo
       Do iE=1, NE
          Read(f_IN,*) Elem_bd(iE)%ID_DoF, Elem_bd(iE)%ID_EL
       EndDo

    Case ('amdba', 'AMDBA')
       NS = size(Nodes_bd)
       NE = size(Elem_bd)
       Write(F_OUT, *) NS, NE
       Do iS=1, NS
          Write(F_OUT,*) iS, Nodes_bd(iS)%Coord%X, Nodes_bd(iS)%Coord%Y, &
               &            Nodes_bd(iS)%ID
       EndDo
       Do iE=1, NE
          Write(F_OUT,*) iE, Elem_bd(iE)%ID_DoF, Elem_bd(iE)%ID_EL
       End Do
    End Select Msh_Fmt
    Close (f_In)

    Return

100 Print*, 'Erreur a l''ouverture du fichier '//Mesh_Str, iStat
    Stop
  End Subroutine Write_Mesh_2D_Scal_Tri_P1

  Subroutine Write_Mesh_2D_Elast_Tri_P1(Mesh_fmt, Mesh_Str, Elem_bd, &
       & Nodes_bd, Disp)
    Character(*)                                        :: Mesh_Fmt
    Character(*)                                        :: Mesh_Str
    Type (Element2D_Elast), Dimension(:), Pointer       :: Elem_bd
    Type (Node2D), Dimension(:), Pointer                :: Nodes_bd
    Real(Kind = Kr), Dimension(:), Pointer, Optional    :: Disp

    Integer                                             :: NE, NS, iDoF, iStat
    Integer                                             :: iS, iE

    Open (Unit=F_Out, file=Mesh_Str, status='unknown', err=100, iostat=iStat)
    Rewind(F_Out)
    Msh_fmt : select case (Mesh_fmt)
    Case ('msh', 'MSH')
       NS = size(Nodes_bd)
       NE = size(Elem_bd)
       Write(F_OUT, *) NS, NE
       Do iS=1, NS,2
          Write(f_OUT,*) Nodes_bd(iS)%Coord%X, Nodes_bd(iS)%Coord%Y, &
               &            Nodes_bd(iS)%ID
       EndDo
       Do iE=1, NE
          Write(f_OUT,*) Elem_bd(iE)%ID_Dof(2)/2, &
               &           Elem_bd(iE)%ID_Dof(4)/2, &
               &           Elem_bd(iE)%ID_Dof(6)/2, &
               &           Elem_bd(iE)%ID_EL
       EndDo

    Case ('amdba', 'AMDBA')
       NS = size(Nodes_bd)
       NE = size(Elem_bd)
       Write(F_OUT, *) NS, NE
       Do iS=1, NS, 2
          Write(F_OUT,*) iS/2+1, Nodes_bd(iS)%Coord%X, Nodes_bd(iS)%Coord%Y, &
               &            Nodes_bd(iS)%ID
       EndDo
       Do iE=1, NE
          Write(F_OUT,*) iE, Elem_bd(iE)%ID_Dof(2)/2, &
               &           Elem_bd(iE)%ID_Dof(4)/2, &
               &           Elem_bd(iE)%ID_Dof(6)/2, &
               &           Elem_bd(iE)%ID_EL
       End Do
    Case ('mesh', 'MESH')
       Write(*,*) 'Warning: mesh format not complete'

       NS = size(Nodes_bd)
       NE = size(Elem_bd)
       
       Write(F_OUT, 200) 'MeshVersionFormatted 1'
       Write(F_OUT, 200) 'Dimension 2'
       Write(F_OUT, 200) '# Incomplete .mesh file created by m_MEF_IO.f'
       Write(F_OUT, 200) 'Vertices'
       Write(F_OUT, *)    NS/2
       If ( Present(Disp) ) Then
          Do iS = 1, NS, 2
             Write(F_OUT, *) Nodes_bd(iS)%Coord%X + Disp(iS),   &
                  &          Nodes_bd(iS)%Coord%Y + Disp(iS+1),         &
                  &          Nodes_bd(iS)%ID
          End Do
       Else
          Do iS = 2, NS, 2
             Write(F_OUT, 201) Nodes_bd(iS)%Coord%X,  Nodes_bd(iS)%Coord%Y , &
                  &          Nodes_bd(iS)%ID
          End Do
       End If
       Write(F_OUT, 200) 'Triangles'
       Write(F_OUT, *)    NE
       Do iE = 1, NE
          Write(F_OUT,202) Elem_bd(iE)%ID_Dof(2)/2, &
               &           Elem_bd(iE)%ID_Dof(4)/2, &
               &           Elem_bd(iE)%ID_Dof(6)/2, &
               &           Elem_bd(iE)%ID_EL
       End Do
    End Select Msh_Fmt

    Close (f_OUT) 
    Return 

100 Print*, 'Erreur a l''ouverture du fichier '//Mesh_Str, iStat
200 Format(A)
201 Format (2(ES13.5),I6)
202 Format(4(I6))
    Stop
  End Subroutine Write_Mesh_2D_Elast_Tri_P1

!$$ ============ Write_bb ============
  Subroutine Write_bb_Vect2D (bb_str,Vect, bb_Type)
    Character(*)                                        :: bb_str
    Type(Vect2D),  Dimension(:), Pointer                :: Vect
    Integer, intent(IN)                                 :: bb_Type

    Integer                                             :: NbVal, NbTerme
    Integer                                             :: iDim = 2
    Integer                                             :: iDum = 0
    Real(Kind = Kr)                                     :: rDum = 0._Kr
    Integer                                             :: iStat, i


    Open (Unit=F_Out, file=bb_Str, status='Unknown', err=100, iostat=iStat)
    Rewind(F_Out)

    NbVal = Size(Vect)
    NbTerme = 2

    Write(F_Out, *) iDim, NbTerme, NbVal, bb_Type

    Do i=1, NbVal
       Write(F_Out,101) Vect(i)%X, Vect(i)%Y
    EndDo
    Write(F_Out,*) iDum, rDum 
    Close(F_Out)
    Return

100 Print*, 'Erreur a l''ouverture du fichier '//bb_Str, iStat
    Stop
101 Format(ES13.5,'  ',ES13.5)
  End Subroutine Write_bb_Vect2D

  Subroutine Write_bb_Vect3D (bb_str,Vect, bb_Type)
    Character(*)                                        :: bb_str
    Type(Vect3D),  Dimension(:), Pointer                :: Vect
    Integer, intent(IN)                                 :: bb_Type

    Integer                                             :: NbVal, NbTerme
    Integer                                             :: iStat,i
    Integer                                             :: iDim = 3
    Integer                                             :: iDum = 0
    Real(Kind = Kr)                                     :: rDum = 0._Kr

    Open (Unit=F_Out, file=bb_Str, status='Unknown', err=100, iostat=iStat)
    Rewind(F_Out)

    NbVal = Size(Vect)
    NbTerme = 3

    Write(F_Out, *) iDim, NbTerme, NbVal, bb_Type

    Do i=1, NbVal
       Write(F_Out,101) Vect(i)%X, Vect(i)%Y, Vect(i)%Z
    EndDo
    Write(F_Out,*) iDum, rDum 
    Close(F_Out)
    Return

100 Print*, 'Erreur a l''ouverture du fichier '//bb_Str, iStat
    Stop
101 Format(ES13.5,'  ',ES13.5,'  ',ES13.5)
  End Subroutine Write_bb_Vect3D

  Subroutine Write_bb_Mat2D (bb_str,Mat, bb_Type)
    Character(*)                                        :: bb_str
    Type(Mat2D),  Dimension(:), Pointer                 :: Mat
    Integer, intent(IN)                                 :: bb_Type

    Integer                                             :: NbVal, NbTerme
    Integer                                             :: iStat,i
    Integer                                             :: iDim = 2
    Integer                                             :: iDum = 0
    Real(Kind = Kr)                                     :: rDum = 0._Kr

    Open (Unit=F_Out, file=bb_Str, status='Unknown', err=100, iostat=iStat)
    Rewind(F_Out)

    NbVal = Size(Mat)
    NbTerme = 4

    Write(F_Out, *) iDim, NbTerme, NbVal, bb_Type

    Do i=1, NbVal
       Write(F_Out,101) Mat(i)%XX, Mat(i)%XY, Mat(I)%YX, Mat(i)%YY
    EndDo
    Write(F_Out,*) iDum, rDum 
    Close(F_Out)
    Return

100 Print*, 'Erreur a l''ouverture du fichier '//bb_Str, iStat
    Stop
101 Format(ES13.5,'  ',ES13.5,'  ',ES13.5,'  ',ES13.5)
  End Subroutine Write_bb_Mat2D

  Subroutine Write_bb_Mat3D (bb_str,Mat, bb_Type)
    Character(*)                                        :: bb_str
    Type(Mat3D),  Dimension(:), Pointer                 :: Mat
    Integer, intent(IN)                                 :: bb_Type

    Integer                                             :: NbVal, NbTerme
    Integer                                             :: iStat,i
    Integer                                             :: iDim = 3
    Integer                                             :: iDum = 0
    Real(Kind = Kr)                                     :: rDum = 0._Kr

    Open (Unit=F_Out, file=bb_Str, status='Unknown', err=100, iostat=iStat)
    Rewind(F_Out)

    NbVal = Size(Mat)
    NbTerme = 9

    Write(F_Out, *) iDim, NbTerme, NbVal, bb_Type

    Do i=1, NbVal
       Write(F_Out,101) Mat(i)%XX, Mat(i)%XY, Mat(i)%XZ,        &
            &           Mat(i)%YX, Mat(i)%YY, Mat(i)%YZ,        &
            &           Mat(i)%ZX, Mat(i)%ZY, Mat(i)%ZZ
    EndDo
    Write(F_Out,*) iDum, rDum 
    Close(F_Out)
    Return

100 Print*, 'Erreur a l''ouverture du fichier '//bb_Str, iStat
    Stop
101 Format(8(ES13.5,'  '),ES13.5)
  End Subroutine Write_bb_Mat3D

  Subroutine Write_bb_MatS2D (bb_str, Mat, bb_Type)
    Character(*)                                        :: bb_str
    Type(MatS2D),  Dimension(:), Pointer                :: Mat
    Integer, intent(IN)                                 :: bb_Type

    Integer                                             :: NbVal, NbTerme
    Integer                                             :: iStat,i
    Integer                                             :: iDim = 2
    Integer                                             :: iDum = 0
    Real(Kind = Kr)                                     :: rDum = 0._Kr

    Open (Unit=F_Out, file=bb_Str, status='Unknown', err=100, iostat=iStat)
    Rewind(F_Out)

    NbVal = Size(Mat)
    NbTerme = 3

    Write(F_Out, *) iDim, NbTerme, NbVal, bb_Type

    Do i=1, NbVal
       Write(F_Out,101) Mat(i)%XX, Mat(i)%YY, Mat(I)%XY
    EndDo
    Write(F_Out,*) iDum, rDum 
    Close(F_Out)
    Return

100 Print*, 'Erreur a l''ouverture du fichier '//bb_Str, iStat
    Stop
101 Format(ES13.5,'  ',ES13.5,'  ',ES13.5)
  End Subroutine Write_bb_MatS2D
  
  Subroutine Write_bb_MatS3D (bb_str,Mat, bb_Type)
    Character(*)                                        :: bb_str
    Type(MatS3D),  Dimension(:), Pointer                :: Mat
    Integer, intent(IN)                                 :: bb_Type

    Integer                                             :: NbVal, NbTerme
    Integer                                             :: iStat,i
    Integer                                             :: iDim = 3
    Integer                                             :: iDum = 0
    Real(Kind = Kr)                                     :: rDum = 0._Kr

    Open (Unit=F_Out, file=bb_Str, status='Unknown', err=100, iostat=iStat)
    Rewind(F_Out)

    NbVal = Size(Mat)
    NbTerme = 6

    Write(F_Out, *) iDim, NbTerme, NbVal, bb_Type

    Do i=1, NbVal
       Write(F_Out,101) Mat(i)%XX, Mat(i)%YY, Mat(i)%ZZ,        &
            &           Mat(i)%YZ, Mat(i)%XZ, Mat(i)%XY
    EndDo
    Write(F_Out,*) iDum, rDum 
    Close(F_Out)
    Return

100 Print*, 'Erreur a l''ouverture du fichier '//bb_Str, iStat
    Stop
101 Format(5(ES13.5,'  '),ES13.5)
  End Subroutine Write_bb_MatS3D

  Subroutine Write_bb_Ptr(bb_str, Ptr, iDim, NbVal, bb_Type)
    Character(*)                                :: bb_str
    Real(Kind = Kr),  Dimension(:), Pointer     :: Ptr
    Integer(Kind = Ki)                          :: iDim, NbVal
    Integer, intent(IN)                         :: bb_Type

    Integer                                     :: NbRec, iStat,i
    Integer                                     :: iDum = 0
    Real(Kind = Kr)                             :: rDum = 0._Kr

!!$ iDim is the dimension of the space 
!!$ NbVal is the number of value for each record
!!$ bb_Type for Piecewise constant or piecewise linear

    Open (Unit=F_Out, file=bb_Str, status='Unknown', err=100, iostat=iStat)
    Rewind(F_Out)

    NbRec = Size(Ptr) / NbVal

    Write(F_Out, *) iDim, NbVal, NbRec, bb_Type
    Do i=1, NbRec
          Write(F_Out,*) Ptr((i-1)*NbVal + 1 : i*NbVal  )
    EndDo
    Write(F_Out,*) iDum, rDum 
    Close(F_Out)
    Return

100 Print*, 'Erreur a l''ouverture du fichier '//bb_Str, iStat
    Stop
  End Subroutine Write_bb_Ptr

!$$ ============ Read_bb ============
  Subroutine Read_bb_1 (bb_str,Vect)
    Character(*)                                :: bb_str
    Real(Kind = Kr), Dimension(:), Pointer      :: Vect

    Integer                                     :: NbVal, NbTerme, iStat, i
    Integer                                     :: iDum 

    Open (Unit=F_In, file=bb_Str, status='Old', err=100, iostat=iStat)
    Rewind(F_In)
    If (Associated(Vect)) Then
       DeAllocate(Vect)
    EndIf

    Read(F_In, *) iDum, NbTerme, NbVal, iDum
    If (NbTerme /= 1) Then
       Print*, 'Warning in Read_bb_1. File ', Trim(bb_Str), ' is not 1D'
    End If
    iDum = iDum
    ! This is only design to prevent from a compiler warning
    
    Allocate(Vect(NbVal))

    Do i=1, NbVal
       Read(F_In,*) Vect(i)
    EndDo

    Close(F_In)
    Return

100 Print*, 'Erreur a l''ouverture du fichier '//bb_Str, iStat
    Stop
  End Subroutine Read_bb_1

  Subroutine Read_bb_2 (bb_str,Vect)
    Character(*)                                :: bb_str
    Type(Vect2D),  Dimension(:), Pointer        :: Vect

    Integer                                     :: NbVal, NbTerme, iStat, i
    Integer                                     :: iDum 

    Open (Unit=F_In, file=bb_Str, status='Old', err=100, iostat=iStat)
    Rewind(F_In)
    If (Associated(Vect)) Then
       DeAllocate(Vect)
    EndIf

    Read(F_In, *) iDum, NbTerme, NbVal, iDum
    iDum = iDum
    If (NbTerme /= 2) Then
       Print*, 'Warning in Read_bb_2. File ', Trim(bb_Str), ' is not 2D'
    End If
    ! This is only designed to prevent from a compiler warning
           
    Allocate(Vect(NbVal))

    Do i=1, NbVal
       Read(F_In,*) Vect(i)%X, Vect(i)%Y
    EndDo

    Close(F_In)
    Return

100 Print*, 'Erreur a l''ouverture du fichier '//bb_Str, iStat
    Stop
  End Subroutine Read_bb_2

  Subroutine Read_bb_3 (bb_str,Vect)
    Character(*)                                :: bb_str
    Type(Vect3D),  Dimension(:), Pointer        :: Vect

    Integer                                     :: NbVal, NbTerme, iStat, i
    Integer                                     :: iDum 

    Open (Unit=F_In, file=bb_Str, status='Old', err=100, iostat=iStat)
    Rewind(F_In)
    If (Associated(Vect)) Then
       DeAllocate(Vect)
    EndIf

    Read(F_In, *) iDum, NbTerme, NbVal, iDum
    iDum = iDum
    ! This is only designed to prevent from a compiler warning
    If (NbTerme /= 3) Then
       Print*, 'Warning in Read_bb_3. File ', Trim(bb_Str), ' is not 3D'
    End If

    Allocate(Vect(NbVal))

    Do i=1, NbVal
       Read(F_In,*) Vect(i)%X, Vect(i)%Y, Vect(i)%Z
    EndDo

    Close(F_In)
    Return

100 Print*, 'Erreur a l''ouverture du fichier '//bb_Str, iStat
    Stop
  End Subroutine Read_bb_3


  Subroutine Read_bb_Ptr (bb_str,Ptr, iDim)
    Character(*)                                :: bb_str
    Real(Kind = Kr), Dimension(:), Pointer      :: Ptr
    Integer(Kind = Ki)                          :: iDim

    Integer                                     :: NbVal, NbTerme,  iStat, i, j
    Integer                                     :: iDum 

    Open (Unit=F_In, file=bb_Str, status='Old', err=100, iostat=iStat)
    Rewind(F_In)
    If (Associated(Ptr)) Then
       DeAllocate(Ptr)
    EndIf

    Read(F_In, *) iDum, NbTerme, NbVal, iDum
    If (NbTerme /= iDim) Then
       Print*, 'Warning in Read_bb_Ptr.'
       Print*, 'File ', Trim(bb_Str), ' is not' ,iDim,'D'
    End If
    iDum = iDum
    ! This is only design to prevent from a compiler warning

    Allocate(Ptr(iDIM * NbVal))

    Do i=1, NbVal
       Do j = 1, iDim
          Read(F_In,*) Ptr((i-1) * iDim + j)
       End Do
    EndDo

    Close(F_In)
    Return

100 Print*, 'Erreur a l''ouverture du fichier '//bb_Str, iStat
    Stop
  End Subroutine Read_bb_Ptr

End Module m_MEF_IO

