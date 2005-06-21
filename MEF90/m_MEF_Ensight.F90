Module m_MEF_Ensight
  Use m_MEF_Types
  Implicit None
  
  ! TODO: Be able to select data given per element or per nodes...

  Interface Write_Geo_ASCII
     Module Procedure Write_Geo_2D_Scal_ASCII, Write_Geo_2D_Elast_ASCII,     &
          &           Write_Geo_3D_Scal_ASCII
  End Interface


  Interface Write_Geo_Bin
     Module Procedure Write_Geo_2D_Scal_BIN, Write_Geo_2D_Elast_BIN,         &
          &           Write_Geo_3D_Scal_BIN
  End Interface

  Interface Write_Result_MatS_ASCII
     Module Procedure Write_Result_MatS2D_ASCII, Write_Result_MatS3D_ASCII
  End Interface

  Interface Write_Result_MatS_BIN
     Module Procedure Write_Result_MatS2D_BIN, Write_Result_MatS3D_BIN
  End Interface
  Contains

    Subroutine Write_Geo_2D_Scal_ASCII(Elem_db, Node_db, Geo_Str, Comment_Str)
      Type(Element2D_Scal), Dimension(:), Pointer      :: Elem_db
      Type(Node2D), Dimension(:), Pointer              :: Node_db
      Character(len = *)                               :: Geo_Str
      Character(len = *), Optional                     :: Comment_Str

      Integer(Kind = Ki)                               :: iS, iE
      Open(File = Geo_Str, Unit = F_OUT, Status = 'Unknown')
      Rewind(F_OUT)
      
      Write(F_OUT,200) Geo_Str
      Write(F_OUT,200) 'Created by Write_Geo_2S_Scal_ASCII'
      
      Write(F_OUT,100) 'node id off'
      Write(F_OUT,100) 'element id off'
      Write(F_OUT,100) 'part'
      Write(F_OUT,300) 1
      Write(F_OUT,200) Geo_Str
      Write(F_OUT,100) 'coordinates'
      Write(F_OUT,300) Size(Node_db)
      Write(F_OUT,500) Node_db%Coord%X
      Write(F_OUT,500) Node_db%Coord%Y
      Write(F_OUT,500) 0.0_Kr * Node_db%Coord%X
      Write(F_OUT,100) 'tria3'
      Write(F_OUT,300) Size(Elem_db)
      Do iE = 1, Size(Elem_db)
         Write(F_OUT,400) Elem_db(iE)%ID_DoF
      End Do

      Close(F_OUT)

100   Format(A)
200   Format(A80)
300   Format(I10)
400   Format(3I10)
500   Format(E12.5)
    End Subroutine Write_Geo_2D_Scal_ASCII

    Subroutine Write_Geo_2D_Scal_BIN(Elem_db, Node_db, Geo_Str, Comment_Str)
      Type(Element2D_Scal), Dimension(:), Pointer     :: Elem_db
      Type(Node2D), Dimension(:), Pointer              :: Node_db
      Character(len = *)                               :: Geo_Str
      Character(len = *), Optional                     :: Comment_Str

      Character(len = 80)                              :: Date_Str
      Integer(Kind = Ki)                               :: iS, iE
      Integer                                          :: NumRec
      Integer                                          :: iKind


      Open(File = Geo_Str, Unit = F_OUT, Status = 'Replace',             &
           & access = 'Direct', form = 'Unformatted', recl = 1)
      
      Write(F_OUT, rec = 1) 'C Binary'
      Write(F_OUT, rec = 81) Geo_Str
!      Call Date_And_Time(date = Date_Str)
!      Write(F_OUT, rec = 161) Date_Str
      Write(F_OUT, rec = 161) 'Created by Write_Geo_2S_Scal_ASCII'
      
      Write(F_OUT, rec = 241) 'node id off'
      Write(F_OUT, rec = 321) 'element id off'
      Write(F_OUT, rec = 401) 'part'
      Write(F_OUT, rec = 481) 1
      Write(F_OUT, rec = 485) Geo_Str
      Write(F_OUT, rec = 565) 'coordinates'
      Write(F_OUT, rec = 645) Size(Node_db)
      Write(F_OUT, rec = 649) real(Node_db%Coord%X, 4)
      NumRec = 649 + 4 * Size(Node_db)
      Write(F_OUT, rec = NumRec) real(Node_db%Coord%Y, 4)
      NumRec = NumRec + 4 * Size(Node_db)
      Write(F_OUT, rec = NumRec) 0.0_4 * real(Node_db%Coord%X,4)
      NumRec = NumRec + 4 * Size(Node_db)
      Write(F_OUT, rec = NumRec) 'tria3'
      NumRec = NumRec + 80
      Write(F_OUT, rec = NumRec) Size(Elem_db)
      NumRec = NumRec + 4
      Do iE = 1, Size(Elem_db)
         Write(F_OUT,rec = NumRec) Elem_db(iE)%ID_DoF
         NumRec = NumRec + 12
      End Do

      Close(F_OUT)

    End Subroutine Write_Geo_2D_Scal_BIN

    Subroutine Write_MatID_2D_Scal_BIN(Elem_db, Node_db, MatID_Str,          &
         & Comment_Str)
      Type(Element2D_Scal), Dimension(:), Pointer     :: Elem_db
      Type(Node2D), Dimension(:), Pointer              :: Node_db
      Character(len = *)                               :: MatID_Str
      Character(len = *), Optional                     :: Comment_Str

      Character(len = 80)                              :: Date_Str
      Integer(Kind = Ki)                               :: iS, iE
      Integer                                          :: NumRec
      Integer                                          :: iKind


      Open(File = MatID_Str, Unit = F_OUT, Status = 'Replace',                &
           & access = 'Direct', form = 'Unformatted', recl = 1)
      
      Write(F_OUT, rec = 1) 'C Binary'
      Write(F_OUT, rec = 81) MatID_Str
      Call Date_And_Time(date = Date_Str)
      Write(F_OUT, rec = 161) Date_Str
      
      Write(F_OUT, rec = 241) 'part'
      Write(F_OUT, rec = 321) 1
      Write(F_OUT, rec = 325) MatID_Str
      NumRec = 405
      Write(F_OUT, rec = NumRec) 'tria3'
      NumRec = NumRec + 80
      Write(F_OUT, rec = NumRec) Size(Elem_db)
      NumRec = NumRec + 4
      Do iE = 1, Size(Elem_db)
         Write(F_OUT,rec = NumRec) Elem_db(iE)%ID_EL
         NumRec = NumRec + 4
      End Do

      Close(F_OUT)

    End Subroutine Write_MatID_2D_Scal_BIN

    Subroutine Write_MatID_2D_Elast_BIN(Elem_db, Node_db, MatID_Str,          &
         & Comment_Str)
      Type(Element2D_Elast), Dimension(:), Pointer     :: Elem_db
      Type(Node2D), Dimension(:), Pointer              :: Node_db
      Character(len = *)                               :: MatID_Str
      Character(len = *), Optional                     :: Comment_Str

      Character(len = 80)                              :: Date_Str
      Integer(Kind = Ki)                               :: iS, iE
      Integer                                          :: NumRec
      Integer                                          :: iKind


      Open(File = MatID_Str, Unit = F_OUT, Status = 'Replace',                &
           & access = 'Direct', form = 'Unformatted', recl = 1)
      
!      Write(F_OUT, rec = 1) 'C Binary'
      Write(F_OUT, rec = 1) MatID_Str
!      Call Date_And_Time(date = Date_Str)
!      Write(F_OUT, rec = 161) Date_Str
      
      Write(F_OUT, rec = 81) 'part'
      Write(F_OUT, rec = 161) 1
      Write(F_OUT, rec = 165) MatID_Str
      NumRec = 245
      Write(F_OUT, rec = NumRec) 'tria3'
      NumRec = NumRec + 80
      Write(F_OUT, rec = NumRec) Size(Elem_db)
      NumRec = NumRec + 4
      Do iE = 1, Size(Elem_db)
         Write(F_OUT,rec = NumRec) Elem_db(iE)%ID_EL
         NumRec = NumRec + 4
      End Do

      Close(F_OUT)

    End Subroutine Write_MatID_2D_Elast_BIN

    Subroutine Write_MatID_2D_Elast_ASCII(Elem_db, Node_db, MatID_Str,        &
         & Comment_Str)
      Type(Element2D_Elast), Dimension(:), Pointer     :: Elem_db
      Type(Node2D), Dimension(:), Pointer              :: Node_db
      Character(len = *)                               :: MatID_Str
      Character(len = *), Optional                     :: Comment_Str

      Integer(Kind = Ki)                               :: iS, iE
      Integer                                          :: NumRec
      Integer                                          :: iKind


      Open(File = MatID_Str, Unit = F_OUT, Status = 'Unknown')
      
      Write(F_OUT, 200) MatID_Str
      Write(F_OUT, 100) 'part'
      Write(F_OUT, 300) 1
      Write(F_OUT, 100) 'tria3'
!      Write(F_OUT, rec = NumRec) Size(Elem_db)
 !     NumRec = NumRec + 4
      Do iE = 1, Size(Elem_db)
         Write(F_OUT, 300) Elem_db(iE)%ID_EL
!         NumRec = NumRec + 4
      End Do

      Close(F_OUT)

100   Format(A)
200   Format(A80)
300   Format(I10)
400   Format(3I10)
500   Format(E12.5)
    End Subroutine Write_MatID_2D_Elast_ASCII

    Subroutine Write_Geo_2D_Elast_ASCII(Elem_db, Node_db, Geo_Str, Comment_Str)
      Type(Element2D_Elast), Dimension(:), Pointer     :: Elem_db
      Type(Node2D), Dimension(:), Pointer              :: Node_db
      Character(len = *)                               :: Geo_Str
      Character(len = *), Optional                     :: Comment_Str

      Character(len = 80)                              :: Date_Str
      Integer(Kind = Ki)                               :: iS, iE
      Integer(Kind = Ki)                               :: Nb_DoF
      Integer(Kind = Ki)                               :: iStat
      Open(File = Geo_Str, Unit = F_OUT, Status = 'Unknown', err = 999,       &
           & iostat = iStat)
      Rewind(F_OUT)

      Write(F_OUT,200) Geo_Str
      Write(F_OUT,200) 'Created by Write_Geo_2D_Elast_ASCII'
      Nb_DoF = Size(Node_db)
      
      Write(F_OUT,100) 'node id off'
      Write(F_OUT,100) 'element id off'
      Write(F_OUT,100) 'part'
      Write(F_OUT,300) 1
      Write(F_OUT,200) Geo_Str
      Write(F_OUT,100) 'coordinates'
      Write(F_OUT,300) Nb_DoF/2
      Write(F_OUT,500) Node_db(1:Nb_DoF:2)%Coord%X
      Write(F_OUT,500) Node_db(1:Nb_DoF:2)%Coord%Y
      Write(F_OUT,500) 0.0_Kr * Node_db(1:Nb_DoF:2)%Coord%X
      Write(F_OUT,100) 'tria3'
      Write(F_OUT,300) Size(Elem_db)
      Do iE = 1, Size(Elem_db)
         Write(F_OUT,400) Elem_db(iE)%ID_DoF(2:6:2)/2
      End Do

      Close(F_OUT)

      Return

100   Format(A)
200   Format(A80)
300   Format(I10)
400   Format(3I10)
500   Format(E12.5)
999   Write(*,*) 'Can''t open file', Geo_Str, iStat
    End Subroutine Write_Geo_2D_Elast_ASCII

    Subroutine Write_Geo_2D_Elast_BIN(Elem_db, Node_db, Geo_Str, Comment_Str)
      Type(Element2D_Elast), Dimension(:), Pointer     :: Elem_db
      Type(Node2D), Dimension(:), Pointer              :: Node_db
      Character(len = *)                               :: Geo_Str
      Character(len = *), Optional                     :: Comment_Str

      Character(len = 80)                              :: Date_Str
      Integer(Kind = Ki)                               :: iS, iE
      Integer                                          :: NumRec
      Integer                                          :: iKind
      Integer(Kind = Ki)                               :: Nb_DoF
      Integer(Kind = Ki)                               :: iStat

      Nb_DoF = Size(Node_db)
      Open(File = Geo_Str, Unit = F_OUT, Status = 'Replace',                  &
           & access = 'Direct', form = 'Unformatted', recl = 1, err = 999,    &
           & iostat = iStat)
      
      Write(F_OUT, rec = 1) 'C Binary'
      Write(F_OUT, rec = 81) Geo_Str
      Call Date_And_Time(date = Date_Str)
      Write(F_OUT, rec = 161) Date_Str
      
      Write(F_OUT, rec = 241) 'node id off'
      Write(F_OUT, rec = 321) 'element id off'
      Write(F_OUT, rec = 401) 'part'
      Write(F_OUT, rec = 481) 1
      Write(F_OUT, rec = 485) Geo_Str
      Write(F_OUT, rec = 565) 'coordinates'
      Write(F_OUT, rec = 645) Nb_DoF/2
      Write(F_OUT, rec = 649) real(Node_db(1:Nb_DoF:2)%Coord%X, 4)
      NumRec = 649 + 4 * Nb_DoF / 2
      Write(F_OUT, rec = NumRec) real(Node_db(1:Nb_DoF:2)%Coord%Y, 4)
      NumRec = NumRec + 4 * Nb_DoF / 2
      Write(F_OUT, rec = NumRec) 0.0_4 * real(Node_db(1:Nb_DoF:2)%Coord%X,4)
      NumRec = NumRec + 4 * Nb_DoF / 2
      Write(F_OUT, rec = NumRec) 'tria3'
      NumRec = NumRec + 80
      Write(F_OUT, rec = NumRec) Size(Elem_db)
      NumRec = NumRec + 4
      Do iE = 1, Size(Elem_db)
         Write(F_OUT,rec = NumRec) Elem_db(iE)%ID_DoF(2:6:2)/2
         NumRec = NumRec + 12
      End Do

      Close(F_OUT)
      Return
999   Write(*,*) 'Can''t open file', Geo_Str, iStat

    End Subroutine Write_Geo_2D_Elast_BIN

    Subroutine Write_Geo_3D_Scal_ASCII(Elem_db, Node_db, Geo_Str, Comment_Str)
      Type(Element3D_Scal), Dimension(:), Pointer     :: Elem_db
      Type(Node3D), Dimension(:), Pointer              :: Node_db
      Character(len = *)                               :: Geo_Str
      Character(len = *), Optional                     :: Comment_Str

      Integer(Kind = Ki)                               :: iS, iE
      Open(File = Geo_Str, Unit = F_OUT, Status = 'Unknown')
      Rewind(F_OUT)

      Write(F_OUT,200) Geo_Str
      Write(F_OUT,200) 'Created by Write_Geo_3D_Elast_ASCII'
      
      Write(F_OUT,100) 'node id off'
      Write(F_OUT,100) 'element id off'
      Write(F_OUT,100) 'part'
      Write(F_OUT,300) 1
      Write(F_OUT,200) Geo_Str
      Write(F_OUT,100) 'coordinates'
      Write(F_OUT,300) Size(Node_db)
      Write(F_OUT,500) Node_db%Coord%X
      Write(F_OUT,500) Node_db%Coord%Y
      Write(F_OUT,500) Node_db%Coord%Z
      Write(F_OUT,100) 'tetra4'
      Write(F_OUT,300) Size(Elem_db)
      Do iE = 1, Size(Elem_db)
         Write(F_OUT,400) Elem_db(iE)%ID_DoF
      End Do

      Close(F_OUT)

100   Format(A)
200   Format(A80)
300   Format(I10)
400   Format(4I10)
500   Format(E12.5)
    End Subroutine Write_Geo_3D_Scal_ASCII

    Subroutine Write_Geo_3D_Scal_BIN(Elem_db, Node_db, Geo_Str, Comment_Str)
      Type(Element3D_Scal), Dimension(:), Pointer     :: Elem_db
      Type(Node3D), Dimension(:), Pointer              :: Node_db
      Character(len = *)                               :: Geo_Str
      Character(len = *), Optional                     :: Comment_Str

      Character(len = 80)                              :: Date_Str
      Integer(Kind = Ki)                               :: iS, iE
      Integer                                          :: NumRec
      Integer                                          :: iKind


      Open(File = Geo_Str, Unit = F_OUT, Status = 'Replace',             &
           & access = 'Direct', form = 'Unformatted', recl = 1)
      
      Write(F_OUT, rec = 1) 'C Binary'
      Write(F_OUT, rec = 81) Geo_Str
      Call Date_And_Time(date = Date_Str)
      Write(F_OUT, rec = 161) Date_Str
      
      Write(F_OUT, rec = 241) 'node id off'
      Write(F_OUT, rec = 321) 'element id off'
      Write(F_OUT, rec = 401) 'part'
      Write(F_OUT, rec = 481) 1
      Write(F_OUT, rec = 485) Geo_Str
      Write(F_OUT, rec = 565) 'coordinates'
      Write(F_OUT, rec = 645) Size(Node_db)
      Write(F_OUT, rec = 649) real(Node_db%Coord%X, 4)
      NumRec = 649 + 4 * Size(Node_db)
      Write(F_OUT, rec = NumRec) real(Node_db%Coord%Y, 4)
      NumRec = NumRec + 4 * Size(Node_db)
      Write(F_OUT, rec = NumRec) real(Node_db%Coord%z,4)
      NumRec = NumRec + 4 * Size(Node_db)
      Write(F_OUT, rec = NumRec) 'tetra4'
      NumRec = NumRec + 80
      Write(F_OUT, rec = NumRec) Size(Elem_db)
      NumRec = NumRec + 4
      Do iE = 1, Size(Elem_db)
         Write(F_OUT,rec = NumRec) Elem_db(iE)%ID_DoF
         NumRec = NumRec + 16
      End Do

      Close(F_OUT)

    End Subroutine Write_Geo_3D_Scal_BIN


    Subroutine Write_Result_Scal_ASCII(DV, Res_Str, Num_Dim, Comment_Str)
      Real(Kind = Kr), Dimension(:), Pointer           :: DV
      Character(len = *)                               :: Res_Str
      Integer, Optional                                :: Num_Dim
      Character(len = *), Optional                     :: Comment_Str

      Integer                                          :: iDim
      Open(File = Res_Str, Unit = F_OUT, Status = 'Unknown')
      Rewind(F_OUT)
      
      Write(F_OUT,200) Trim(Res_Str)
      Write(F_OUT,100) 'part'
      Write(F_OUT,300) 1
      Write(F_OUT,100) 'coordinates'
      If (Present(Num_Dim)) Then
         Do iDim = 1, Num_Dim
            Write(F_OUT,500) DV(iDim:size(DV):Num_Dim)
         End Do
      Else
         !!! Write the values in the order they were received
         !!! Will not make any sense to ensight in dim>1
         Write(F_OUT,500) DV
      End If

      Close(F_OUT)

100   Format(A)
200   Format(A80)
300   Format(I10)
500   Format(E12.5)
    End Subroutine Write_Result_Scal_ASCII

    Subroutine Write_Result_Scal_BIN(DV, Res_Str, Num_Dim, Comment_Str)
      Real(Kind = Kr), Dimension(:), Pointer           :: DV
      Character(len = *)                               :: Res_Str
      Integer, Optional                                :: Num_Dim
      Character(len = *), Optional                     :: Comment_Str

      Character(len = 80)                              :: Date_Str
      Integer                                          :: iDim, nVals

      Open(File = Res_Str, Unit = F_OUT, Status = 'Replace',             &
           & access = 'Direct', form = 'Unformatted', recl = 1)
      
      nVals = Size(DV)
!      Write(F_OUT, rec = 1) 'C Binary'
      Call Date_And_Time(date = Date_Str)
      Write(F_OUT, rec = 1) Trim(Res_Str) // ' ' // Trim(Date_Str)
      Write(F_OUT, rec = 81) 'part'
      Write(F_OUT, rec = 161) 1
      Write(F_OUT, rec = 165) 'coordinates'
      If (Present(Num_Dim)) Then
         Do iDim = 1, Num_Dim
            Write(F_OUT, rec=245+(iDim-1)*nVals/Num_Dim * 4)                  &
                 & real( DV(iDim:nVals:Num_Dim), 4) 
         End Do
      Else
         !!! Write the values in the order they were received
         !!! Will not make any sense to ensight in dim>1
         Write(F_OUT, rec = 245) real(DV, 4)
      End If

      Close(F_OUT)
    End Subroutine Write_Result_Scal_BIN
      
    Subroutine Write_Result_2D_ASCII(DV, Res_Str, Comment_Str)
      Real(Kind = Kr), Dimension(:), Pointer           :: DV
      Character(len = *)                               :: Res_Str
      Character(len = *), Optional                     :: Comment_Str

      Integer(Kind = Ki)                               :: iSize
      Open(File = Res_Str, Unit = F_OUT, Status = 'Unknown')
      Rewind(F_OUT)
      iSize = Size(DV)
      
      Write(F_OUT,200) Trim(Res_Str)
      Write(F_OUT,100) 'part'
      Write(F_OUT,300) 1
      Write(F_OUT,100) 'coordinates'
      Write(F_OUT,500) DV(1:iSize:2)
      Write(F_OUT,500) DV(2:iSize:2)
      Write(F_OUT,500) 0.0_Kr * DV(2:iSize:2)


      Close(F_OUT)

100   Format(A)
200   Format(A80)
300   Format(I10)
500   Format(E12.5)
    End Subroutine Write_Result_2D_ASCII

    Subroutine Write_Result_2D_BIN(DV, Res_Str, Comment_Str)
      Real(Kind = Kr), Dimension(:), Pointer           :: DV
      Character(len = *)                               :: Res_Str
      Character(len = *), Optional                     :: Comment_Str

      Integer(Kind = Ki)                               :: iSize
      Character(len = 80)                              :: Date_Str
      
      Open(File = Res_Str, Unit = F_OUT, Status = 'Replace',             &
           & access = 'Direct', form = 'Unformatted', recl = 1)
      
      iSize = Size(DV)
      Call Date_And_Time(date = Date_Str)
      Write(F_OUT, rec = 1) Trim(Res_Str) // ' ' // Trim(Date_Str)
      Write(F_OUT, rec = 81) 'part'
      Write(F_OUT, rec = 161) 1
      Write(F_OUT, rec = 165) 'coordinates'
      Write(F_OUT, rec = 245) real(DV(1:iSize:2), 4)
      Write(F_OUT, rec= 245+2*iSize) real(DV(2:iSize:2),4)
      Write(F_OUT, rec= 245+4*iSize) 0.0_4 * real(DV(2:iSize:2),4)

      Close(F_OUT)
    End Subroutine Write_Result_2D_BIN
      
    Subroutine Write_Result_MatS2D_ASCII(DV, Res_Str, Comment_Str)
      Type(MatS2D), Dimension(:), Pointer              :: DV
      Character(len = *)                               :: Res_Str
      Character(len = *), Optional                     :: Comment_Str

      Integer(Kind = Ki)                               :: iSize
      Real(Kind = Kr), Dimension(:), Pointer           :: Zeros

      Open(File = Res_Str, Unit = F_OUT, Status = 'Unknown')
      Rewind(F_OUT)
      iSize = Size(DV)
      Allocate(Zeros(iSize))
      Zeros = 0.0_Kr
      
      Write(F_OUT,200) Trim(Res_Str)
      Write(F_OUT,100) 'part'
      Write(F_OUT,300) 1
      Write(F_OUT,100) 'tria3'
      Write(F_OUT,500) DV%XX
      Write(F_OUT,500) DV%YY
      Write(F_OUT,500) Zeros
      Write(F_OUT,500) DV%XY
      Write(F_OUT,500) Zeros
      Write(F_OUT,500) Zeros

      Close(F_OUT)
      DeAllocate(Zeros)

100   Format(A)
200   Format(A80)
300   Format(I10)
500   Format(E12.5)
    End Subroutine Write_Result_MatS2D_ASCII

    Subroutine Write_Result_MatS2D_BIN(DV, Res_Str, Comment_Str)
      Type(MatS2D), Dimension(:), Pointer              :: DV
      Character(len = *)                               :: Res_Str
      Character(len = *), Optional                     :: Comment_Str

      Integer(Kind = Ki)                               :: iSize, iRec
      Character(len = 80)                              :: Date_Str
      Real(Kind = 4), Dimension(:), Pointer            :: Zeros
      
      Open(File = Res_Str, Unit = F_OUT, Status = 'Replace',             &
           & access = 'Direct', form = 'Unformatted', recl = 1)
      
      iSize = Size(DV)
      Allocate(Zeros(iSize))
      Zeros = 0.0_4

      Call Date_And_Time(date = Date_Str)
      Write(F_OUT, rec = 1) Trim(Res_Str) // ' ' // Trim(Date_Str)
      Write(F_OUT, rec = 81) 'part'
      Write(F_OUT, rec = 161) 1
      Write(F_OUT, rec = 165) 'tria3'
      Write(F_OUT, rec = 245) real(DV%XX, 4)
      iRec = 245 + 4 * iSize 
      Write(F_OUT, rec= iRec) real(DV%YY, 4)
      iRec = iRec + 4 * iSize
      Write(F_OUT, rec= iRec) Zeros
      iRec = iRec + 4 * iSize
      Write(F_OUT, rec= iRec) real(DV%XY, 4)
      iRec = iRec + 4 * iSize
      Write(F_OUT, rec= iRec) Zeros
      iRec = iRec + 4 * iSize
      Write(F_OUT, rec= iRec) Zeros
      Close(F_OUT)

      DeAllocate(Zeros)
    End Subroutine Write_Result_MatS2D_BIN


    Subroutine Write_Result_3D_ASCII(DV, Res_Str, Comment_Str)
      Real(Kind = Kr), Dimension(:), Pointer           :: DV
      Character(len = *)                               :: Res_Str
      Character(len = *), Optional                     :: Comment_Str

      Integer(Kind = Ki)                               :: iSize
      Open(File = Res_Str, Unit = F_OUT, Status = 'Unknown')
      Rewind(F_OUT)
      iSize = Size(DV)
      
      Write(F_OUT,200) Trim(Res_Str)
      Write(F_OUT,100) 'part'
      Write(F_OUT,300) 1
      Write(F_OUT,100) 'coordinates'
      Write(F_OUT,500) DV(1:iSize:3)
      Write(F_OUT,500) DV(2:iSize:3)
      Write(F_OUT,500) DV(3:iSize:3)


      Close(F_OUT)

100   Format(A)
200   Format(A80)
300   Format(I10)
500   Format(E12.5)
    End Subroutine Write_Result_3D_ASCII

    Subroutine Write_Result_3D_BIN(DV, Res_Str, Comment_Str)
      Real(Kind = Kr), Dimension(:), Pointer           :: DV
      Character(len = *)                               :: Res_Str
      Character(len = *), Optional                     :: Comment_Str

      Integer(Kind = Ki)                               :: iSize
      Character(len = 80)                              :: Date_Str
      
      Open(File = Res_Str, Unit = F_OUT, Status = 'Replace',             &
           & access = 'Direct', form = 'Unformatted', recl = 1)
      
      iSize = Size(DV)
      Call Date_And_Time(date = Date_Str)
      Write(F_OUT, rec = 1) Trim(Res_Str) // ' ' // Trim(Date_Str)
      Write(F_OUT, rec = 81) 'part'
      Write(F_OUT, rec = 161) 1
      Write(F_OUT, rec = 165) 'coordinates'
      Write(F_OUT, rec = 245) real(DV(1:iSize:3), 4)
      Write(F_OUT, rec= 245+2*iSize) real(DV(2:iSize:3),4)
      Write(F_OUT, rec= 245+4*iSize) real(DV(3:iSize:3),4)

      Close(F_OUT)
    End Subroutine Write_Result_3D_BIN
      
    Subroutine Write_Result_MatS3D_ASCII(DV, Res_Str, Comment_Str)
      Type(MatS3D), Dimension(:), Pointer              :: DV
      Character(len = *)                               :: Res_Str
      Character(len = *), Optional                     :: Comment_Str


      Open(File = Res_Str, Unit = F_OUT, Status = 'Unknown')
      Rewind(F_OUT)
      
      Write(F_OUT,200) Trim(Res_Str) 
      Write(F_OUT,100) 'part'
      Write(F_OUT,300) 1
      Write(F_OUT,100) 'tetra4'
      Write(F_OUT,500) DV%XX
      Write(F_OUT,500) DV%YY
      Write(F_OUT,500) DV%ZZ
      Write(F_OUT,500) DV%XY
      Write(F_OUT,500) DV%YZ
      Write(F_OUT,500) DV%XZ

      Close(F_OUT)

100   Format(A)
200   Format(A80)
300   Format(I10)
500   Format(E12.5)
    End Subroutine Write_Result_MatS3D_ASCII

    Subroutine Write_Result_MatS3D_BIN(DV, Res_Str, Comment_Str)
      Type(MatS3D), Dimension(:), Pointer              :: DV
      Character(len = *)                               :: Res_Str
      Character(len = *), Optional                     :: Comment_Str

      Integer(Kind = Ki)                               :: iSize, iRec
      Character(len = 80)                              :: Date_Str
      
      Open(File = Res_Str, Unit = F_OUT, Status = 'Replace',             &
           & access = 'Direct', form = 'Unformatted', recl = 1)
      
      iSize = Size(DV)

      Call Date_And_Time(date = Date_Str)
      Write(F_OUT, rec = 1) Trim(Res_Str) // ' ' // Trim(Date_Str)
      Write(F_OUT, rec = 81) 'part'
      Write(F_OUT, rec = 161) 1
      Write(F_OUT, rec = 165) 'tetra4'
      Write(F_OUT, rec = 245) real(DV%XX, 4)
      iRec = 245 + 4 * iSize 
      Write(F_OUT, rec= iRec) real(DV%YY, 4)
      iRec = iRec + 4 * iSize
      Write(F_OUT, rec= iRec) real(DV%ZZ, 4)
      iRec = iRec + 4 * iSize
      Write(F_OUT, rec= iRec) real(DV%XY, 4)
      iRec = iRec + 4 * iSize
      Write(F_OUT, rec= iRec) real(DV%YZ, 4)
      iRec = iRec + 4 * iSize
      Write(F_OUT, rec= iRec) real(DV%XZ, 4)
      Close(F_OUT)

    End Subroutine Write_Result_MatS3D_BIN


  End Module m_MEF_Ensight
