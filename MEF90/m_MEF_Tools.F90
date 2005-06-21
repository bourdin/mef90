module m_MEF_Tools
  ! Blaise Bourdin, 1996-1998
  ! Merci de faire parvenir toutes remarques et bugs 
  ! eventuels aux adresses suivantes :
  ! bourdin@lpmtm.univ-paris13.fr
  ! bourdin@mat.dtu.dk
  !
  ! En cas de modifications, les noms de  modules et de 
  ! fichiers _DOIVENT_ etre renommes
  !          ^^^^^^^^^
  Use m_MEF_Types

  IMPLICIT NONE

  Interface GC
     Module Procedure GC_Morse2, GC_Prof
  End Interface

  Interface PMV
     Module Procedure PMV_Morse, PMV_Prof
  End Interface

  Interface Elim_Dir
     Module Procedure Elim_Dir_Morse, Elim_Dir_Prof
  End Interface

  Interface Init_Matrix
     Module Procedure Init_Morse_2D_Scal, Init_Morse_2D_Elast,                &
          & Init_Morse_1D, Init_Profil_2D_Scal, Init_Profil_2D_Elast,         &
          & Init_Profil_1D
  End Interface

  Interface Get_Matrix_Size_Symm
     Module Procedure Get_Matrix_Size_2D_Scal_Symm,                           &
          & Get_Matrix_Size_2D_Elast_Symm,                                    &
          & Get_Matrix_Size_3D_Scal_Symm,                                     &
          & Get_Matrix_Size_3D_Elast_Symm
  End Interface

  Interface Get_Matrix_Size
     Module Procedure Get_Matrix_Size_2D_Scal, Get_Matrix_Size_2D_Elast,      &
          & Get_Matrix_Size_3D_Scal, Get_Matrix_Size_3D_Elast
  End Interface

  Interface LLT_Factor
     Module Procedure LLT_Factor_Prof
  End Interface

  Interface LLT_Solve
     Module Procedure LLT_Solve_Prof
  End Interface

Contains

  Subroutine GC_Morse (A, B, X, iPrec, eps, iStat)
    Type (Mat_Morse)                            :: A
    Real(Kind = Kr), Dimension(:), Pointer      :: B
    Real(Kind = Kr), Dimension(:), Pointer      :: X
    Integer                                     :: iPrec, iStat
    Real(Kind = Kr)                             :: eps

    Integer(Kind = Ki)                          :: iDim, iCpt, iL, iC
    Real(Kind = Kr)                             :: GCErr, Awn_Wn
    Real(Kind = Kr)                             :: AWn_Wn_I
    Real(Kind = Kr)                             :: Rhon, Gamman
    Real(Kind = Kr), Dimension(:), Pointer      :: Gn, Wn, AWn, GnP

    iDim = Ubound(B, dim=1)
    Allocate(Gn(iDim))
    Allocate(Wn(iDim))
    Allocate(AWn(iDim))
    Allocate(GnP(iDim))

    ! INITIALISATIONS :
    iCpt = 1

    Call PMV_Morse(A, X, Gn)
    Gn = Gn - B
    Wn = Gn
    GCErr = SQRT(Dot_Product(Gn, Gn))
    ! ITERATION
    Iteration : Do While ((iCPt < iDim+1) .AND. (GCErr > eps))
       Call PMV_Morse(A, Wn, AWn)
       AWn_Wn = Dot_Product(AWn,Wn)
!       TestDiv : If ( AWn_Wn < Epsilon(1.0_Kr) ) Then
       TestDiv : If ( AWn_Wn == 0.0_Kr ) Then
          Print*, 'Divergence du Gradient Conjugue', iCpt, AWn_Wn
          Stop
       EndIf TestDiv
       AWn_Wn_I = 1.0_Kr / AWn_Wn
       Rhon = Dot_Product(Wn,Gn) * AWn_Wn_I

       X = X - RhoN * Wn
       Gn = Gn - Rhon * AWn
       Precond : Select Case (iPrec)
       Case (1) ! Preconditionnement diagonal
          ! Find the diagonal term (This won't be really needed as soon 
          ! as the terms in the column will be sorted decreasingly
          Do iL = 1, iDim
             Do iC = 1, A%Row(iL)%NB_Col
                If (A%Row(iL)%Col(iC)%Pos == iL) Then
                   Exit
                End If
             End Do

             GnP(iL) = Gn(iL) / A%Row(iL)%Col(iC)%Val
          End Do
          Gamman = DbleM1 * Dot_Product(GnP,Awn) * AWn_Wn_I
          Wn = GnP + Gamman * Wn
       Case Default
          Gamman = DbleM1 * Dot_Product(Gn, AWn) * AWn_Wn_I
          Wn = Gn + Gamman * Wn

       End Select Precond
       GCErr = SQRT(Dot_Product(Gn, Gn))
       Print*, GCErr
       Print*, iCpt
       iCpt = iCpt + 1
    End Do Iteration
    If (GCErr > eps) Then
       iStat = -1
    Else
       iStat = iCpt
    EndIf

    DeAllocate(Gn)
    DeAllocate(Wn)
    DeAllocate(GnP)
    DeAllocate(AWn)
  End Subroutine GC_Morse

  Subroutine GC_Morse2 (A, B, X, iPrec, eps, iStat, Residuals)
    Type (Mat_Morse)                            :: A
    Real(Kind = Kr), Dimension(:), Pointer      :: B
    Real(Kind = Kr), Dimension(:), Pointer      :: X
    Integer                                     :: iPrec, iStat
    Real(Kind = Kr)                             :: eps
    Real(Kind = Kr), Dimension(:), Pointer, optional :: Residuals

    Integer(Kind = Ki)                          :: iDim, iCpt, iL, iC
    Real(Kind = Kr)                             :: GCErr, Awn_Wn
    Real(Kind = Kr)                             :: AWn_Wn_I
    Real(Kind = Kr)                             :: Rhon, Gamman
    Real(Kind = Kr), Dimension(:), Pointer      :: Gn, Wn, AWn, GnP
    Real(Kind = Kr), Dimension(:), Pointer      :: Tmp_Residuals

    Integer(Kind = Ki)                          :: iL_PMV, iC_PMV, Curr_pos
    Real(Kind = Kr)                             :: Curr_Val
    Real(Kind = Kr)                             :: NormB
    Real(Kind = Kr)                             :: NormRes

    iDim = Ubound(B, dim=1)
    Allocate(Gn(iDim))
    Allocate(Wn(iDim))
    Allocate(AWn(iDim))
    Allocate(GnP(iDim))
    Allocate(Tmp_Residuals(iDim))

    ! INITIALISATIONS :
    iCpt = 1
    NormB = sqrt(Dot_Product(B,B))

    AWn = 0.0_Kr
    Call PMV_Morse(A, X, Gn)
    Gn = Gn - B
    Wn = Gn
    GCErr = SQRT(Dot_Product(Gn, Gn))
    ! ITERATION
    Iteration : Do While ((iCPt < iDim+1) .AND. (GCErr > eps))
       Tmp_Residuals(iCpt) = GCErr
!       Print*, iCpt
!       Print*, GCErr
       Call PMV_Morse(A, Wn, AWn)
!!$!       The Matrix Vector is hard coded here:
!!$       AWn = 0.0_Kr
!!$       Do iL_PMV = 1, iDim
!!$          ! Scaning the Rows of the matrix
!!$          Do iC_PMV = 1, A%Row(iL_PMV)%NB_Col
!!$             ! Scaning the non zero terms of the current Row
!!$             Curr_pos = A%Row(iL_PMV)%Col(iC_PMV)%pos
!!$             !           Print*, Curr_Pos
!!$             ! Curr_pos is the column number of the current term of the matrix
!!$             Curr_val = A%Row(iL_PMV)%Col(iC_PMV)%Val
!!$             ! Curr_val is its value
!!$             
!!$             AWn(iL_PMV) = AWn(iL_PMV) + Curr_Val * Wn(Curr_Pos)
!!$             ! Contribution of the lower part
!!$             
!!$             If (iL_PMV /= Curr_Pos) Then
!!$                AWn(Curr_Pos) = AWn(Curr_Pos) + Curr_Val * Wn(iL_PMV)
!!$                ! Contribution of the upper part (So we don't have to
!!$                ! repeat the scan twice)
!!$             End If
!!$          End Do
!!$       End Do
!!$!      End of the Matrix-vector product

       AWn_Wn = Dot_Product(AWn,Wn)
       TestDiv : If ( AWn_Wn == 0.0_Kr ) Then
          Print*, 'Divergence du Gradient Conjugue', iCpt, AWn_Wn
          Stop
       EndIf TestDiv
       AWn_Wn_I = 1.0_Kr / AWn_Wn
       Rhon = Dot_Product(Wn,Gn) * AWn_Wn_I

       X = X - RhoN * Wn
!       X = X + RhoN * Wn
       Gn = Gn - Rhon * AWn
       Precond : Select Case (iPrec)
       Case (1) ! Preconditionnement diagonal
          ! Find the diagonal term (This won't be really needed as soon 
          ! as the terms in the column will be sorted decreasingly
          Do iL = 1, iDim
             Do iC = 1, A%Row(iL)%NB_Col
                If (A%Row(iL)%Col(iC)%Pos == iL) Then
                   Exit
                End If
             End Do

             GnP(iL) = Gn(iL) / A%Row(iL)%Col(iC)%Val
          End Do
          Gamman = DbleM1 * Dot_Product(GnP,Awn) * AWn_Wn_I
          Wn = GnP + Gamman * Wn
       Case Default
          Gamman = DbleM1 * Dot_Product(Gn, AWn) * AWn_Wn_I
          Wn = Gn + Gamman * Wn

       End Select Precond
       GCErr = SQRT(Dot_Product(Gn, Gn)) / NormB
       iCpt = iCpt + 1
    End Do Iteration
    If (GCErr > eps) Then
       iStat = -1
    Else
       iStat = iCpt
    EndIf
    If (Present(Residuals)) Then
       If (Associated(Residuals)) Then
          DeAllocate(Residuals)
       End If
       Allocate(Residuals(iCpt))
       Residuals = Tmp_Residuals(1:iCpt)
    End If

    DeAllocate(Tmp_Residuals)
    DeAllocate(Gn)
    DeAllocate(Wn)
    DeAllocate(GnP)
    DeAllocate(AWn)
  End Subroutine GC_Morse2


  Subroutine PMV_Morse(MM1, V1, VRes)
    Type (Mat_Morse)                            :: MM1
    Real(Kind=Kr), Dimension(:), Pointer        :: V1, VRes

    Integer(Kind = Ki)                          :: iL, iC, Curr_pos
    Integer(Kind = Ki)                          :: NL
    Real(Kind = Kr)                             :: Curr_Val

    NL = Size(MM1%Row)
    If (Associated (Vres)) Then
       Continue
    Else
       Allocate (VRes(NL))
    EndIf
    VRes = 0.0_Kr

    Do iL = 1, NL
       ! Scaning the Rows of the matrix
       Do iC = 1, MM1%Row(iL)%NB_Col
          ! Scaning the non zero terms of the current Row
          Curr_pos = MM1%Row(iL)%Col(iC)%pos
!           Print*, Curr_Pos
          ! Curr_pos is the column number of the current term of the matrix
          Curr_val = MM1%Row(iL)%Col(iC)%Val
          ! Curr_val is its value

!          Print*, 'iL, iC', iL, iC
!          Print*, ' Pos, Val, V1', Curr_pos, Curr_val, V1(Curr_Pos)
!          Print*, 'VRes, iL, iC', VRes(iL), iL, iC
!          Print*, 'V1    ', Curr_Pos, iC
!!$          Print*, '============================================='
          VRes(iL) = VRes(iL) + Curr_Val * V1(Curr_Pos)
          ! Contribution of the lower part

          If (iL /= Curr_Pos) Then
             VRes(Curr_Pos) = VRes(Curr_Pos) + Curr_Val * V1(iL)
             ! Contribution of the upper part (So we don't have to
             ! repeat the scan twice)
          End If
       End Do
    End Do
  End Subroutine PMV_Morse

  Subroutine Init_Morse_2D_Scal(MR, Elem_bd, Nod_bd)
    Type (Mat_Morse)                                    :: MR
    Type (Element2D_Scal), Dimension(:), Pointer       :: Elem_bd
    Type (Node2D), DImension(:), Pointer                :: Nod_bd

    Integer(Kind = Ki), Dimension(:), Pointer   :: TMP_Row
    Integer(Kind = Ki)                          :: NN, iL
    Integer(Kind = Ki)                          :: NE, iE
    Integer(Kind = Ki)                          :: iSL_L, iSL_C
    Integer(Kind = Ki)                          :: iSG_L, iSG_C
    Integer(Kind = Ki)                          :: NVals
    Integer(Kind = Ki)                          :: iCol
    Integer(Kind = Ki)                          :: iSwap
    Integer(Kind = Ki)                          :: Is_present
    Integer(Kind = Ki)                          :: Size_Inc = 10
    If (Associated (MR%Row)) Then
       Print*, 'MR should not be already allocated...'
!       stop
    Else
       NN = Size(Nod_bd)
       NE = Size(Elem_bd)
       Allocate (MR%Row(NN))
    End If
    !    Allocate (Row_Size(NN))
    !    Row_Size = 0_Ki
    Do iL = 1, NN
       MR%Row(iL)%NB_Col = 0_Ki
       ! reserve some space for the column to be filled
       Allocate (MR%Row(iL)%Col(Size_Inc))
    End Do

    Do iE = 1, NE
       Do iSL_L = 1, Elem_bd(iE)%NB_DoF
          Do iSL_C = 1, Elem_bd(iE)%NB_DoF
             iSG_L = Elem_bd(iE)%ID_DoF(iSL_L)
             iSG_C = Elem_bd(iE)%ID_DoF(iSL_C)
             ! Check if we are in the lower part 
             If (iSG_L >= iSG_C) Then
                ! Is iSG_C already present in the current Row
                Is_Present = 0

                Do iCol = 1, MR%Row(iSG_L)%NB_Col
                   If (MR%Row(iSG_L)%Col(iCol)%Pos == iSG_C) Then
                      Is_Present = 1
                      ! iSG_C has been found, exit the loop
                      Exit
                   End If
                End Do

                If (Is_Present == 0) Then
                   NVals = MR%Row(iSG_L)%NB_Col
                   ! Check if we need to reallocate the column
                   If (MR%Row(iSG_L)%NB_Col == Size(MR%Row(iSG_L)%Col)) Then
                      Allocate(TMP_Row(NVals + Size_Inc))
                      Do iSwap = 1, NVals
                         TMP_Row(iSwap) = MR%Row(iSG_L)%Col(iSwap)%Pos
                      End Do
                      DeAllocate(MR%Row(iSG_L)%Col)
                      Allocate(MR%Row(iSG_L)%Col(NVals + Size_Inc))
                      Do iSwap = 1, NVals
                         MR%Row(iSG_L)%Col(iSwap)%Pos = TMP_Row(iSwap)
                      End Do
                      DeAllocate(TMP_Row)
                   EndIf
                   MR%Row(iSG_L)%NB_Col = NVals + 1
                   MR%Row(iSG_L)%Col(NVals + 1)%Pos = iSG_C
                End If
             End If
          End Do
       End Do
    End Do

    ! resize the columns :
    Do iSG_L = 1, NN
       ! TODO :
       ! This is the perfect place for sorting the terms ...
       NVals = MR%Row(iSG_L)%NB_Col
       Allocate(TMP_Row(NVals))
       Do iSwap = 1, NVals
          TMP_Row(iSwap) = MR%Row(iSG_L)%Col(iSwap)%Pos
       End Do
       DeAllocate(MR%Row(iSG_L)%Col)
       Allocate(MR%Row(iSG_L)%Col(NVals))
       Do iSwap = 1, NVals
          MR%Row(iSG_L)%Col(iSwap)%Pos = TMP_Row(iSwap)
          MR%Row(iSG_L)%Col(iSwap)%Val = 0.0_Kr
       End Do
       DeAllocate(TMP_Row)
    End Do

  End Subroutine Init_Morse_2D_Scal


  Subroutine Init_Morse_2D_Elast(MR, Elem_bd, Nod_bd)
    Type (Mat_Morse)                                    :: MR
    Type (Element2D_Elast), Dimension(:), Pointer       :: Elem_bd
    Type (Node2D), Dimension(:), Pointer                :: Nod_bd

    Integer(Kind = Ki), Dimension(:), Pointer   :: TMP_Row
    Integer(Kind = Ki)                          :: NN, iL
    Integer(Kind = Ki)                          :: NE, iE
    Integer(Kind = Ki)                          :: iSL_L, iSL_C
    Integer(Kind = Ki)                          :: iSG_L, iSG_C
    Integer(Kind = Ki)                          :: NVals
    Integer(Kind = Ki)                          :: iCol
    Integer(Kind = Ki)                          :: iSwap
    Integer(Kind = Ki)                          :: Is_present
    Integer(Kind = Ki)                          :: Size_Inc = 10
!    Integer(Kind = Ki)                         :: iSL
!    If (Associated (MR%Row)) Then
!       Print*, 'MR should not be already allocated...'
!       stop
!    Else
       NN = Size(Nod_bd)
       NE = Size(Elem_bd)
       Allocate (MR%Row(NN))
!    End If
    !    Allocate (Row_Size(NN))
    !    Row_Size = 0_Ki
    Do iL = 1, NN
       MR%Row(iL)%NB_Col = 0_Ki
       ! reserve some space for the column to be filled
       Allocate (MR%Row(iL)%Col(Size_Inc))
    End Do

    Do_iE: Do iE = 1, NE
       Do_iSL_1: Do iSL_L = 1, Elem_bd(iE)%NB_DoF
          Do_iSLC: Do iSL_C = 1, Elem_bd(iE)%NB_DoF
             iSG_L = Elem_bd(iE)%ID_DoF(iSL_L)
             iSG_C = Elem_bd(iE)%ID_DoF(iSL_C)
             ! Check if we are in the lower part 
             Is_Lower: If (iSG_L >= iSG_C) Then
                ! Is iSG_C already present in the current Row
                Is_Present = 0

                Do_iClo: Do iCol = 1, MR%Row(iSG_L)%NB_Col
                   If (MR%Row(iSG_L)%Col(iCol)%Pos == iSG_C) Then
                      Is_Present = 1
                      ! iSG_C has been found, exit the loop
                      Exit
                   End If
                End Do Do_iClo

                If_Is_Present: If (Is_Present == 0) Then
                   NVals = MR%Row(iSG_L)%NB_Col
                   ! Check if we need to reallocate the column
                   If (MR%Row(iSG_L)%NB_Col == Size(MR%Row(iSG_L)%Col)) Then
                      Allocate(TMP_Row(NVals + Size_Inc))
                      Do iSwap = 1, NVals
                         TMP_Row(iSwap) = MR%Row(iSG_L)%Col(iSwap)%Pos
                      End Do
                      DeAllocate(MR%Row(iSG_L)%Col)
                      Allocate(MR%Row(iSG_L)%Col(NVals + Size_Inc))
                      Do iSwap = 1, NVals
                         MR%Row(iSG_L)%Col(iSwap)%Pos = TMP_Row(iSwap)
                      End Do
                      DeAllocate(TMP_Row)
                   EndIf
                   MR%Row(iSG_L)%NB_Col = NVals + 1
                   MR%Row(iSG_L)%Col(NVals + 1)%Pos = iSG_C
                End If If_Is_Present
             End If Is_Lower
          End Do Do_iSLC
       End Do Do_iSL_1
    End Do Do_iE
    ! resize the columns :
    Do_iL_2: Do iL = 1, NN
       ! TODO :
       ! This is the perfect place for sorting the terms ...
       NVals = MR%Row(iL)%NB_Col
       Allocate(TMP_Row(NVals))
       Do iSwap = 1, NVals
          TMP_Row(iSwap) = MR%Row(iL)%Col(iSwap)%Pos
       End Do
       DeAllocate(MR%Row(iL)%Col)
       Allocate(MR%Row(iL)%Col(NVals))
       Do iSwap = 1, NVals
          MR%Row(iL)%Col(iSwap)%Pos = TMP_Row(iSwap)
          MR%Row(iL)%Col(iSwap)%Val = 0.0_Kr
       End Do
       DeAllocate(TMP_Row)
    End Do Do_iL_2
  End Subroutine Init_Morse_2D_Elast

  Subroutine Init_Morse_1D(MR, Elem_bd, Nod_bd)
    Type (Mat_Morse)                            :: MR
    Type (Element1D), Dimension(:), Pointer     :: Elem_bd
    Type (Node2D), Dimension(:), Pointer        :: Nod_bd

    Integer(Kind = Ki), Dimension(:), Pointer   :: TMP_Row
    Integer(Kind = Ki)                          :: NN, iL
    Integer(Kind = Ki)                          :: NE, iE
    Integer(Kind = Ki)                          :: iSL_L, iSL_C
    Integer(Kind = Ki)                          :: iSG_L, iSG_C
    Integer(Kind = Ki)                          :: NVals
    Integer(Kind = Ki)                          :: iCol
    Integer(Kind = Ki)                          :: iSwap    
    Integer(Kind = Ki)                          :: Is_present
    Integer(Kind = Ki)                          :: Size_Inc = 10
    If (Associated (MR%Row)) Then
       Print*, 'MR should not be already allocated...'
!       stop
    Else
       NN = Size(Nod_bd)
       NE = Size(Elem_bd)
       Allocate (MR%Row(NN))
    End If
    !    Allocate (Row_Size(NN))
    !    Row_Size = 0_Ki
    Do iL = 1, NN
       MR%Row(iL)%NB_Col = 0_Ki
       ! reserve some space for the column to be filled
       Allocate (MR%Row(iL)%Col(Size_Inc))
    End Do

    Do iE = 1, NE
       Do iSL_L = 1, Elem_bd(iE)%NB_DoF
          Do iSL_C = 1, Elem_bd(iE)%NB_DoF
             iSG_L = Elem_bd(iE)%ID_DoF(iSL_L)
             iSG_C = Elem_bd(iE)%ID_DoF(iSL_C)
             ! Check if we are in the lower part 
             If (iSG_L >= iSG_C) Then
                ! Is iSG_C already present in the current Row
                Is_Present = 0

                Do iCol = 1, MR%Row(iSG_L)%NB_Col
                   If (MR%Row(iSG_L)%Col(iCol)%Pos == iSG_C) Then
                      Is_Present = 1
                      ! iSG_C has been found, exit the loop
                      Exit
                   End If
                End Do

                If (Is_Present == 0) Then
                   NVals = MR%Row(iSG_L)%NB_Col
                   ! Check if we need to reallocate the column
                   If (MR%Row(iSG_L)%NB_Col == Size(MR%Row(iSG_L)%Col)) Then
                      Allocate(TMP_Row(NVals + Size_Inc))
                      Do iSwap = 1, NVals
                         TMP_Row(iSwap) = MR%Row(iSG_L)%Col(iSwap)%Pos
                      End Do
                      DeAllocate(MR%Row(iSG_L)%Col)
                      Allocate(MR%Row(iSG_L)%Col(NVals + Size_Inc))
                      Do iSwap = 1, NVals
                         MR%Row(iSG_L)%Col(iSwap)%Pos = TMP_Row(iSwap)
                      End Do
                      DeAllocate(TMP_Row)
                   EndIf
                   MR%Row(iSG_L)%NB_Col = NVals + 1
                   MR%Row(iSG_L)%Col(NVals + 1)%Pos = iSG_C
                End If
             End If
          End Do
       End Do
    End Do
    ! resize the columns :
    Do iSL_L = 1, NN
       ! TODO :
       ! This is the perfect place for sorting the terms ...
       NVals = MR%Row(iSL_L)%NB_Col
       Allocate(TMP_Row(NVals))
       Do iSwap = 1, NVals
          TMP_Row(iSwap) = MR%Row(iSL_L)%Col(iSwap)%Pos
       End Do
       Deallocate(Mr%Row(iSL_L)%Col)
       Allocate(MR%Row(iSL_L)%Col(NVals))
       Do iSwap = 1, NVals
          MR%Row(iSL_L)%Col(iSwap)%Pos = TMP_Row(iSwap)
          MR%Row(iSL_L)%Col(iSwap)%Val = 0.0_Kr
       End Do
       DeAllocate(TMP_Row)
    End Do
  End Subroutine Init_Morse_1D

  Subroutine LLT_Factor_Prof(A)
    Type (Mat_Prof)             :: A

    Integer(Kind = Ki)          :: i, j, k
    Integer(Kind = Ki)          :: NE
    Real(Kind = Kr)             :: InvOfDiag

    NE = UBound(A%Profil,1)
    Do i = 1, NE-1
       If ( A%Val(A%Profil(i)) <= 0) Then
          Print*, 'Something went wrong in LLT_Factor:'
          Print*, 'Diagonal term ', i, '< 0'
          Stop
       End If
       
       !1. A(i,i) <-- sqrt(A(i,i))
       A%Val(A%Profil(i)) = sqrt(A%Val(A%Profil(i)))
       InvOfDiag = 1.0_Kr / A%Val(A%Profil(i))

       !2. A(i+1:n,i) <- A(i+1:n,i) / A(i,i)
       Do j = i+1, NE
          ! Is (j,i) in the profile ?
          If (i >= j+1 - A%Profil(j) + A%Profil(j-1)) Then
             ! A(j,i) = A%Val(A%Profil(j)+i-j)
             A%Val(A%Profil(j)+i-j) = A%Val(A%Profil(j)+i-j) * &
                  & InvOfDiag
          End If
       End Do

       !3. A(j,k) <- A(j,k) - A(j,i) * A(k,i) , j = i+1, NE, k = i+1,j
       Do j = i+1, NE
          Do k = j,NE
             If (k >= j+1 - A%Profil(j) + A%Profil(j-1)) Then
                If (i >= j+1 - A%Profil(j) + A%Profil(j-1)) Then
                   If (i >= k+1 - A%Profil(k) + A%Profil(k-1)) Then
                      A%Val(A%Profil(k)+j-k)  = A%Val(A%Profil(k)+j-k) - &
                           & A%Val(A%Profil(j)+i-j) * &
                           & A%Val(A%Profil(k)+i-k)
                   End If
                End If
             End If
          End Do
       End Do
    End Do
    A%Val(A%Profil(NE)) = sqrt(A%Val(A%Profil(NE)))
  End Subroutine LLT_Factor_Prof

  Subroutine LLT_Solve_Prof(L,X,b)
    Type (Mat_Prof)                             :: L
    Real(Kind = Kr), Dimension(:), Pointer      :: X,b
    
    Real(Kind = Kr), Dimension(:), Pointer      :: VTmp

    Allocate(VTmp(Size(b)))
    Call LLT_SolveInf_Prof(L,VTmp,b)
    Print*, 'LLT_Inf', MinVal(VTmp), MaxVal(VTmp)
    Call LLT_SolveSup_Prof(L,X,VTmp)
    Print*, 'LLT_Inf', MinVal(X), MaxVal(X)

    DeAllocate(VTmp)
  End Subroutine LLT_Solve_Prof

  Subroutine LLT_SolveInf_Prof(L,X,b)
    Type (Mat_Prof)                             :: L
    Real(Kind = Kr), Dimension(:), Pointer      :: X,b

    Integer(Kind = Ki)                          :: i, NE

    NE = Size(X)
    X = b
    X(1) = b(1) / L%Val(1)
    Do i = 2, NE
       X(i) = X(i) - Dot_Product(X(i-L%Profil(i)+L%Profil(i-1)+1:i-1), &
            & L%Val(L%Profil(i-1)+1:L%Profil(i)-1)) 
       X(i) = X(i) / L%Val(L%Profil(i))
    End Do
  End Subroutine LLT_SolveInf_Prof
  
  Subroutine LLT_SolveSup_Prof(L,X,b)
    Type (Mat_Prof)                             :: L
    Real(Kind = Kr), Dimension(:), Pointer      :: X,b

    Integer(Kind = Ki)                          :: i,j, NE

    NE = Size(X)
    X = b
    X(NE) = b(NE) / L%Val(L%Profil(NE))

    Do i = NE-1, 1, -1
       Do j = i+1, NE
          If (i >= j - L%Profil(j) + L%Profil(j-1) + 1) Then
             X(i) = X(i) - X(j) * L%Val(L%Profil(j)+i-j)
          End If
       End Do
       X(i) = X(i) / L%Val(L%Profil(i))
    End Do
  End Subroutine LLT_SolveSup_Prof

  Subroutine PMV_Prof(MP1, V1, VRes)
    Type (Mat_Prof)                             :: MP1
    Real(Kind = Kr), Pointer, Dimension(:)      :: V1
    Real(Kind = Kr), Pointer, Dimension(:)      :: Vres
    Integer                                     :: iDim, iC, iL
    Integer                                     :: Prof_iL
    iDim = Size(V1)
    If (Associated (Vres)) Then
       Continue
    Else
       Allocate (VRes(iDim))
    EndIf
    !Allocate (Vres(iDIM))

    VRes = 0.0_Kr
    Do iL = 1, iDim
       !       Vres(iL) = 0.0_Kr
       Prof_iL = MP1%Profil(iL)
       Do iC = iL-Prof_Il+MP1%Profil(iL-1)+1, iL-1
          ! On parcourt la partie tri. inf de la ligne
          Vres(iL) = Vres(iL) + &
               & MP1%Val(Prof_iL+iC-iL) * V1(iC)
          ! On rajoute les Vals correspondant a la partie tri sup
          Vres(iC) = Vres(iC) + &
               & MP1%Val(Prof_Il+iC-iL) * V1(iL)
       EndDo
       Vres(iL) = Vres(iL) + &
            & MP1%Val(Prof_Il) * V1(iL)
    EndDo
  End Subroutine PMV_Prof

  Subroutine GC_Prof (A, B, X, iPrec, eps, iStat, Residuals)
    Type (Mat_Prof)                             :: A
    Real(Kind = Kr), Dimension(:), Pointer      :: B
    Real(Kind = Kr), Dimension(:), Pointer      :: X
    Integer                                     :: iPrec, iStat
    Real(Kind = Kr)                             :: eps
    Real(Kind = Kr), Dimension(:), Pointer, Optional :: Residuals

    Integer                                     :: iDim, iCpt
    Real(Kind = Kr)                             :: GCErr, Awn_Wn
    Real(Kind = Kr)                             :: AWn_Wn_I
    Real(Kind = Kr)                             :: Rhon, Gamman
    Real(Kind = Kr), Dimension(:), Pointer      :: Gn, Wn, AWn, GnP
    Real(Kind = Kr), Dimension(:), Pointer      :: Tmp_Residuals

    iDim = Ubound(B, dim=1)
    Allocate(Gn(iDim))
    Allocate(Wn(iDim))
    Allocate(AWn(iDim))
    Allocate(GnP(iDim))
    Allocate(Tmp_Residuals(iDim))

    ! INITIALISATIONS :
    iCpt = 1

    Call PMV_Prof(A, X, Gn)
    Gn = Gn - B
    Wn = Gn
    GCErr = SQRT(Dot_Product(Gn, Gn))
    ! ITERATION
    Iteration : Do While ((iCPt < iDim+1) .AND. (GCErr > eps))
       Tmp_Residuals(iCpt) = GCErr
       Call PMV_Prof(A, Wn, AWn)
       AWn_Wn = Dot_Product(AWn,Wn)
       TestDiv : If ( AWn_Wn == 0.0_Kr ) Then
          Print*, 'Divergence du Gradient Conjugue'
          Stop
       EndIf TestDiv
       AWn_Wn_I = 1.0_Kr / AWn_Wn

       Rhon = Dot_Product(Wn,Gn) * AWn_Wn_I

       X = X - RhoN * Wn
       Gn = Gn - Rhon * AWn
       Precond : Select Case (iPrec)
       Case (1) ! Preconditionnement diagonal
          GnP(1:iDIM) = Gn(1:iDIM) / A%Val(A%Profil(1:iDIM))
          Gamman = DbleM1 * Dot_Product(GnP,Awn) * AWn_Wn_I
          Wn = GnP + Gamman * Wn
       Case Default
          Gamman = DbleM1*Dot_Product(Gn, AWn) * AWn_Wn_I
          Wn = Gn + Gamman * Wn

       End Select Precond
       GCErr = SQRT(Dot_Product(Gn, Gn))
       iCpt = iCpt + 1
    End Do Iteration
    If (GCErr > eps) Then
       iStat = -1
    Else
       iStat = iCpt
    EndIf
    If (Present(Residuals)) Then
       If(Associated(Residuals)) Then
          DeAllocate(Residuals)
       End If
       Allocate(Residuals(iCpt))
       Residuals = Tmp_Residuals(1:iCpt)
       DeAllocate(Tmp_Residuals)
    End If


    DeAllocate(Gn)
    DeAllocate(Wn)
    DeAllocate(GnP)
    DeAllocate(AWn)
  End Subroutine GC_Prof

  Subroutine Elim_Dir_Morse(MM1, Liste_Som)
    Type(Mat_Morse)                             :: MM1
    Integer(Kind = Ki), Dimension(:), Pointer   :: Liste_Som
!    Type(Node2D)                               :: Nod_bd

    Integer(Kind = Ki), Dimension(:), Pointer   :: Nod_IDs
    Integer(Kind = Ki)                          :: iL, iC
    Integer(Kind = Ki)                          :: Curr_Pos
!    Integer(Kind = Ki)                         :: iCod_CL
    Integer(Kind = Ki)                          :: NNod

    NNod = Size(MM1%Row)
    Allocate (Nod_IDs(NNod)) 
    Nod_Ids = 0_Ki
    Do iL = 1, Size(Liste_Som)
       Nod_Ids(Liste_Som(iL)) = 1_Ki
    End Do
       
    Do_iL: Do iL = 1, Size(MM1%Row)
       Do_iC: Do iC = 1, MM1%Row(iL)%NB_Col
          Curr_Pos = MM1%Row(iL)%Col(iC)%Pos
          Is_Cl: If ((Nod_Ids(Curr_Pos) == 1_Ki) .OR. &
               & (Nod_Ids(iL) == 1_Ki)) Then
             Is_Diag: If (Curr_Pos == iL) Then
! This is a diagonal term => set to 1
                MM1%Row(iL)%Col(iC)%Val = 1.0_Kr
             Else
                MM1%Row(iL)%Col(iC)%Val = 0.0_Kr
             End If Is_Diag
          End If Is_Cl
       End Do Do_iC
    End Do Do_iL
    DeAllocate (Nod_Ids)
    100 Format(I3, I3, I3)
  End Subroutine Elim_Dir_Morse
  


  Subroutine Elim_Dir_Prof(MP1, ListeSom)
    Type (Mat_Prof)                             :: MP1
    Integer (Kind = Ki), dimension(:), Pointer  :: ListeSom
    Integer                                     :: iDim, NbL, i, iL, iC

    iDim = ubound(MP1%Profil,1)
    NbL = Size(ListeSom)
    Do i = 1, NbL
       iL = ListeSom(i)
       ! Partie Tri_Inf :
       Do iC = MP1%Profil(iL-1)+1, MP1%Profil(iL)-1
          MP1%Val(iC) = Dble0
       EndDo
       ! Val Diagonal :
       MP1%Val(MP1%Profil(iL)) = Dble1
       ! Partie Tri_Sup :
       Do iC = iL + 1, iDim
          ! Le Val de coordonnees (iL, iC) est il dans le profil ?
          If (iL >= iC+1-MP1%Profil(iC)+MP1%Profil(iC-1)) Then
             ! La position du Val (iL,iC) est MP1%Profil(iC)+iL-iC
             MP1%Val(MP1%Profil(iC)+iL-iC) = Dble0
          EndIf
       EndDo
    EndDo
  End Subroutine Elim_Dir_Prof

  Subroutine Init_profil_2DA(MR, Elem_bd, Nodes_bd)
    Type (Mat_Prof)                             :: MR
    Type (Element2DA), Dimension(:), Pointer    :: Elem_bd
    Type (Node2D), Dimension(:), Pointer        :: Nodes_bd

    Integer                                     :: NS, NE, Nb_DoF
    Integer                                     :: iE, iSL1, ISG1
    Integer                                     :: iSL2, iSG2, iS

    NS = size(Nodes_bd)
    NE = Size(Elem_bd)
    Allocate(MR%Profil(0:NS))

    Do iS=0, NS
       MR%Profil(iS) = iS
    EndDo

    ! On determine la longueur de chaque demi ligne  
    Do iE = 1, NE
       Nb_DoF = Elem_bd(iE)%NB_DoF
       Do iSl1 = 1, Nb_DoF
          iSG1 = Elem_bd(iE)%ID_DoF(iSL1)
          Do iSL2 = 1, Nb_DoF
             iSG2 = Elem_bd(iE)%ID_DoF(iSL2)
             MR%Profil(iSG1) = min(MR%Profil(iSG1), iSG2)
          EndDo
       EndDo
    EndDo

    ! On calcule la position du Val diagonal dans MR%Val :
    MR%Profil(0) = 0
    MR%Profil(1) = 1
    Do iSG1 = 2, NS
       MR%Profil(iSG1) = MR%Profil(iSG1-1)-MR%Profil(iSG1)+iSG1+1
    EndDo

    ! On peut Maintenant Dimensionner MR%Val :
    Allocate (MR%Val(MR%Profil(NS)))  

  End Subroutine Init_Profil_2DA

  Subroutine Init_profil_2D_Scal(MR, Elem_bd, Nodes_bd)
    Type (Mat_Prof)                                     :: MR
    Type (Element2D_Scal), Dimension(:), Pointer       :: Elem_bd
    Type (Node2D), Dimension(:), Pointer                :: Nodes_bd

    Integer                                             :: NS, NE, Nb_DoF
    Integer                                             :: iE, iSL1, ISG1
    Integer                                             ::  iSL2, iSG2, iS

    NS = size(Nodes_bd)
    NE = Size(Elem_bd)
    Allocate(MR%Profil(0:NS))

    Do iS=0, NS
       MR%Profil(iS) = iS
    EndDo

    ! On determine la longueur de chaque demi ligne  
    Do iE = 1, NE
       Nb_DoF = Elem_bd(iE)%NB_DoF
       Do iSl1 = 1, Nb_DoF
          iSG1 = Elem_bd(iE)%ID_DoF(iSL1)
          Do iSL2 = 1, Nb_DoF
             iSG2 = Elem_bd(iE)%ID_DoF(iSL2)
             MR%Profil(iSG1) = min(MR%Profil(iSG1), iSG2)
          EndDo
       EndDo
    EndDo

    ! On calcule la position du Val diagonal dans MR%Val :
    MR%Profil(0) = 0
    MR%Profil(1) = 1
    Do iSG1 = 2, NS
       MR%Profil(iSG1) = MR%Profil(iSG1-1)-MR%Profil(iSG1)+iSG1+1
    EndDo

    ! On peut Maintenant Dimensionner MR%Val :
    Allocate (MR%Val(MR%Profil(NS)))  

  End Subroutine Init_Profil_2D_Scal

  Subroutine Init_profil_2D_Elast(MR, Elem_bd, Nodes_bd)
    Type (Mat_Prof)                                     :: MR
    Type (Element2D_Elast), Dimension(:), Pointer       :: Elem_bd
    Type (Node2D), Dimension(:), Pointer                :: Nodes_bd

    Integer                                             :: NS, NE, Nb_DoF
    Integer                                             :: iE, iSL1, ISG1
    Integer                                             :: iSL2, iSG2, iS
    !    Integer                                             :: NbDim

    !    NbDim = 2  
    NS = size(Nodes_bd)
    NE = Size(Elem_bd)
    Allocate(MR%Profil(0:NS))

    Do iS=0, NS
       MR%Profil(iS) = iS
    EndDo

    ! On determine la longueur de chaque demi ligne  
    Do iE = 1, NE
       Nb_DoF = Elem_bd(iE)%NB_DoF
       Do iSl1 = 1, Nb_DoF
          iSG1 = Elem_bd(iE)%ID_DoF(iSL1)
          Do iSL2 = 1, Nb_DoF
             iSG2 = Elem_bd(iE)%ID_DoF(iSL2)
             MR%Profil(iSG1) = min(MR%Profil(iSG1), iSG2)
          EndDo
       EndDo
    EndDo

    ! On calcule la position du Val diagonal dans MR%Val :
    MR%Profil(0) = 0
    MR%Profil(1) = 1
    Do iSG1 = 2, NS
       MR%Profil(iSG1) = MR%Profil(iSG1-1)-MR%Profil(iSG1)+iSG1+1
    EndDo

    ! On peut Maintenant Dimensionner MR%Val :
    Allocate (MR%Val(MR%Profil(NS)))  

  End Subroutine Init_Profil_2D_Elast

  Subroutine Init_profil_1D(MR, Elem_bd, Nodes_bd)
    Type (Mat_Prof)                             :: MR
    Type (Element1D), Dimension(:), Pointer     :: Elem_bd
    Type (Node1D), Dimension(:), Pointer        :: Nodes_bd

    Integer                                     :: NS, NE, Nb_DoF
    Integer                                     :: iE, iSL1, ISG1
    Integer                                     :: iSL2, iSG2, iS

    NS = size(Nodes_bd)
    NE = Size(Elem_bd)
    Allocate(MR%Profil(0:NS))

    Do iS=0, NS
       MR%Profil(iS) = iS
    EndDo

    ! On determine la longueur de chaque demi ligne  
    Do iE = 1, NE
       Nb_DoF = Elem_bd(iE)%NB_DoF
       Do iSl1 = 1, Nb_DoF
          iSG1 = Elem_bd(iE)%ID_DoF(iSL1)
          Do iSL2 = 1, Nb_DoF
             iSG2 = Elem_bd(iE)%ID_DoF(iSL2)
             MR%Profil(iSG1) = min(MR%Profil(iSG1), iSG2)
          EndDo
       EndDo
    EndDo

    ! On calcule la position du Val diagonal dans MR%Val :
    MR%Profil(0) = 0
    MR%Profil(1) = 1
    Do iSG1 = 2, NS
       MR%Profil(iSG1) = MR%Profil(iSG1-1)-MR%Profil(iSG1)+iSG1+1
    EndDo

    ! On peut Maintenant Dimensionner MR%Val :
    Allocate (MR%Val(MR%Profil(NS)))  

  End Subroutine Init_Profil_1D

  Subroutine Get_Matrix_Size_2D_Scal_Symm(NNZ, Elem_db, Node_db)
    
    Type Mat_AA_Entry
       Real(Kind = Kr)                            :: Val
       Integer(Kind = Ki)                         :: Pos
    End Type Mat_AA_Entry
    
    Type Mat_AA_Row
       Integer(Kind=Ki)                            :: NB_Col
       Type(Mat_AA_Entry), dimension(:), pointer   :: Col
    End Type Mat_AA_Row
    
    
    Integer, Dimension(:), Pointer                 :: NNZ
    Type(Element2D_Scal), Dimension(:), Pointer    :: Elem_db
    Type(Node2D), Dimension(:), Pointer            :: Node_db
    
    Integer(Kind = Ki)                             :: NS, NE, Nb_DoF
    Integer(Kind = Ki)                             :: iE, iSL1, ISG1
    Integer(Kind = Ki)                             :: iSL2, iSG2, iS, iT
    Type(Mat_AA_Row), Dimension(:), pointer        :: AA
    Integer(Kind = Ki), Dimension(:), Pointer      :: TMP_Row
    Integer(Kind = Ki)                             :: iL
    Integer(Kind = Ki)                             :: NVals
    Integer(Kind = Ki)                             :: iCol
    Integer(Kind = Ki)                             :: iSwap
    Integer(Kind = Ki)                             :: Is_present
    Integer(Kind = Ki)                             :: Size_Inc = 10



    NS = Size(Node_db)
    NE = Size(Elem_db)
    Allocate (NNZ(NS))
    NNZ = 0_Ki
    Allocate (AA(NS))
    
    Do iL = 1, NS
       AA(iL)%NB_Col = 0_Ki
       Allocate (AA(iL)%Col(Size_Inc))
    End Do
    
    Do_iE: Do iE = 1, NE
       Do_iSL_1: Do iSL1 = 1, Elem_db(iE)%NB_DoF
          iSG1 = Elem_db(iE)%ID_DoF(iSL1)
          Do_iSLC: Do iSL2 = 1, Elem_db(iE)%NB_DoF
             iSG2 = Elem_db(iE)%ID_DoF(iSL2)
             ! Check if we are in the lower part 
             Is_Lower: If (iSG1 <= iSG2) Then
                ! Is iSG2 already present in the current Row
                Is_Present = 0
                
                Do_iClo: Do iCol = 1, AA(iSG1)%NB_Col
                   If (AA(iSG1)%Col(iCol)%Pos == iSG2) Then
                      Is_Present = 1
                      ! iSG2 has been found, exit the loop
                      Exit
                   End If
                End Do Do_iClo
                
                If_Is_Present: If (Is_Present == 0) Then
                   NVals = AA(iSG1)%NB_Col
                   If (AA(iSG1)%NB_Col == Size(AA(iSG1)%Col)) Then
                      Allocate(TMP_Row(NVals + Size_Inc))
                      Do iSwap = 1, NVals
                         TMP_Row(iSwap) = AA(iSG1)%Col(iSwap)%Pos
                         If (AA(iSG1)%Col(iSwap)%Val /= 0.0_Kr) Then
                         End If
                      End Do
                      DeAllocate(AA(iSG1)%Col)
                      Allocate(AA(iSG1)%Col(NVals + Size_Inc))
                      Do iSwap = 1, NVals
                         AA(iSG1)%Col(iSwap)%Pos = TMP_Row(iSwap)
                      End Do
                      DeAllocate(TMP_Row)
                   EndIf
                   AA(iSG1)%NB_Col = NVals + 1
                   AA(iSG1)%Col(NVals + 1)%Pos = iSG2
                End If If_Is_Present
             End If Is_Lower
          End Do Do_iSLC
       End Do Do_iSL_1
    End Do Do_iE
    Do iS = 1, NS
       NNZ(iS)=AA(iS)%NB_Col
    End Do
  End Subroutine Get_Matrix_Size_2D_Scal_Symm

  Subroutine Get_Matrix_Size_2D_Scal(NNZ, Elem_db, Node_db)
    
    Type Mat_AA_Entry
       Real(Kind = Kr)                            :: Val
       Integer(Kind = Ki)                         :: Pos
    End Type Mat_AA_Entry
    
    Type Mat_AA_Row
       Integer(Kind=Ki)                            :: NB_Col
       Type(Mat_AA_Entry), dimension(:), pointer   :: Col
    End Type Mat_AA_Row
    
    
    Integer, Dimension(:), Pointer                 :: NNZ
    Type(Element2D_Scal), Dimension(:), Pointer    :: Elem_db
    Type(Node2D), Dimension(:), Pointer            :: Node_db
    
    Integer(Kind = Ki)                             :: NS, NE, Nb_DoF
    Integer(Kind = Ki)                             :: iE, iSL1, ISG1
    Integer(Kind = Ki)                             :: iSL2, iSG2, iS, iT
    Type(Mat_AA_Row), Dimension(:), pointer        :: AA
    Integer(Kind = Ki), Dimension(:), Pointer      :: TMP_Row
    Integer(Kind = Ki)                             :: iL
    Integer(Kind = Ki)                             :: NVals
    Integer(Kind = Ki)                             :: iCol
    Integer(Kind = Ki)                             :: iSwap
    Integer(Kind = Ki)                             :: Is_present
    Integer(Kind = Ki)                             :: Size_Inc = 10



    NS = Size(Node_db)
    NE = Size(Elem_db)
    Allocate (NNZ(NS))
    NNZ = 0_Ki
    Allocate (AA(NS))
    
    Do iL = 1, NS
       AA(iL)%NB_Col = 0_Ki
       Allocate (AA(iL)%Col(Size_Inc))
    End Do

    Do_iE: Do iE = 1, NE
       Do_iSL_1: Do iSL1 = 1, Elem_db(iE)%NB_DoF
          iSG1 = Elem_db(iE)%ID_DoF(iSL1)
          Do_iSLC: Do iSL2 = 1, Elem_db(iE)%NB_DoF
             iSG2 = Elem_db(iE)%ID_DoF(iSL2)
             ! Check if we are in the lower part 
!             Is_Lower: If (iSG1 <= iSG2) Then
                ! Is iSG2 already present in the current Row
                Is_Present = 0
                
                Do_iClo: Do iCol = 1, AA(iSG1)%NB_Col
                   If (AA(iSG1)%Col(iCol)%Pos == iSG2) Then
                      Is_Present = 1
                      ! iSG2 has been found, exit the loop
                      Exit
                   End If
                End Do Do_iClo
                
                If_Is_Present: If (Is_Present == 0) Then
                   NVals = AA(iSG1)%NB_Col
                   If (AA(iSG1)%NB_Col == Size(AA(iSG1)%Col)) Then
                      Allocate(TMP_Row(NVals + Size_Inc))
                      Do iSwap = 1, NVals
                         TMP_Row(iSwap) = AA(iSG1)%Col(iSwap)%Pos
                         If (AA(iSG1)%Col(iSwap)%Val /= 0.0_Kr) Then
                         End If
                      End Do
                      DeAllocate(AA(iSG1)%Col)
                      Allocate(AA(iSG1)%Col(NVals + Size_Inc))
                      Do iSwap = 1, NVals
                         AA(iSG1)%Col(iSwap)%Pos = TMP_Row(iSwap)
                      End Do
                      DeAllocate(TMP_Row)
                   EndIf
                   AA(iSG1)%NB_Col = NVals + 1
                   AA(iSG1)%Col(NVals + 1)%Pos = iSG2
                End If If_Is_Present
 !            End If Is_Lower
          End Do Do_iSLC
       End Do Do_iSL_1
    End Do Do_iE
    Do iS = 1, NS
       NNZ(iS)=AA(iS)%NB_Col
    End Do
  End Subroutine Get_Matrix_Size_2D_Scal

  Subroutine Get_Matrix_Size_2D_Elast_Symm(NNZ, Elem_db, Node_db)
    
    Type Mat_AA_Entry
       Real(Kind = Kr)                            :: Val
       Integer(Kind = Ki)                         :: Pos
    End Type Mat_AA_Entry
    
    Type Mat_AA_Row
       Integer(Kind=Ki)                            :: NB_Col
       Type(Mat_AA_Entry), dimension(:), pointer   :: Col
    End Type Mat_AA_Row
    
    
    Integer, Dimension(:), Pointer                 :: NNZ
    Type(Element2D_Elast), Dimension(:), Pointer   :: Elem_db
    Type(Node2D), Dimension(:), Pointer            :: Node_db
    
    Integer(Kind = Ki)                             :: NS, NE, Nb_DoF
    Integer(Kind = Ki)                             :: iE, iSL1, ISG1
    Integer(Kind = Ki)                             :: iSL2, iSG2, iS, iT
    Type(Mat_AA_Row), Dimension(:), pointer        :: AA
    Integer(Kind = Ki), Dimension(:), Pointer      :: TMP_Row
    Integer(Kind = Ki)                             :: iL
    Integer(Kind = Ki)                             :: NVals
    Integer(Kind = Ki)                             :: iCol
    Integer(Kind = Ki)                             :: iSwap
    Integer(Kind = Ki)                             :: Is_present
    Integer(Kind = Ki)                             :: Size_Inc = 10



    NS = Size(Node_db)
    NE = Size(Elem_db)
    Allocate (NNZ(NS))
    NNZ = 0_Ki
    Allocate (AA(NS))
    
    Do iL = 1, NS
       AA(iL)%NB_Col = 0_Ki
       Allocate (AA(iL)%Col(Size_Inc))
    End Do
    
    Do_iE: Do iE = 1, NE
       Do_iSL_1: Do iSL1 = 1, Elem_db(iE)%NB_DoF
          iSG1 = Elem_db(iE)%ID_DoF(iSL1)
          Do_iSLC: Do iSL2 = 1, Elem_db(iE)%NB_DoF
             iSG2 = Elem_db(iE)%ID_DoF(iSL2)
             ! Check if we are in the lower part 
             Is_Lower: If (iSG1 <= iSG2) Then
                ! Is iSG2 already present in the current Row
                Is_Present = 0
                
                Do_iClo: Do iCol = 1, AA(iSG1)%NB_Col
                   If (AA(iSG1)%Col(iCol)%Pos == iSG2) Then
                      Is_Present = 1
                      ! iSG2 has been found, exit the loop
                      Exit
                   End If
                End Do Do_iClo
                
                If_Is_Present: If (Is_Present == 0) Then
                   NVals = AA(iSG1)%NB_Col
                   If (AA(iSG1)%NB_Col == Size(AA(iSG1)%Col)) Then
                      Allocate(TMP_Row(NVals + Size_Inc))
                      Do iSwap = 1, NVals
                         TMP_Row(iSwap) = AA(iSG1)%Col(iSwap)%Pos
                         If (AA(iSG1)%Col(iSwap)%Val /= 0.0_Kr) Then
                         End If
                      End Do
                      DeAllocate(AA(iSG1)%Col)
                      Allocate(AA(iSG1)%Col(NVals + Size_Inc))
                      Do iSwap = 1, NVals
                         AA(iSG1)%Col(iSwap)%Pos = TMP_Row(iSwap)
                      End Do
                      DeAllocate(TMP_Row)
                   EndIf
                   AA(iSG1)%NB_Col = NVals + 1
                   AA(iSG1)%Col(NVals + 1)%Pos = iSG2
                End If If_Is_Present
             End If Is_Lower
          End Do Do_iSLC
       End Do Do_iSL_1
    End Do Do_iE
    Do iS = 1, NS
       NNZ(iS)=AA(iS)%NB_Col
    End Do
  End Subroutine Get_Matrix_Size_2D_Elast_Symm

  Subroutine Get_Matrix_Size_2D_Elast(NNZ, Elem_db, Node_db)
    
    Type Mat_AA_Entry
       Real(Kind = Kr)                            :: Val
       Integer(Kind = Ki)                         :: Pos
    End Type Mat_AA_Entry
    
    Type Mat_AA_Row
       Integer(Kind=Ki)                            :: NB_Col
       Type(Mat_AA_Entry), dimension(:), pointer   :: Col
    End Type Mat_AA_Row
    
    
    Integer, Dimension(:), Pointer                 :: NNZ
    Type(Element2D_Elast), Dimension(:), Pointer   :: Elem_db
    Type(Node2D), Dimension(:), Pointer            :: Node_db
    
    Integer(Kind = Ki)                             :: NS, NE, Nb_DoF
    Integer(Kind = Ki)                             :: iE, iSL1, ISG1
    Integer(Kind = Ki)                             :: iSL2, iSG2, iS, iT
    Type(Mat_AA_Row), Dimension(:), pointer        :: AA
    Integer(Kind = Ki), Dimension(:), Pointer      :: TMP_Row
    Integer(Kind = Ki)                             :: iL
    Integer(Kind = Ki)                             :: NVals
    Integer(Kind = Ki)                             :: iCol
    Integer(Kind = Ki)                             :: iSwap
    Integer(Kind = Ki)                             :: Is_present
    Integer(Kind = Ki)                             :: Size_Inc = 10



    NS = Size(Node_db)
    NE = Size(Elem_db)
    Allocate (NNZ(NS))
    NNZ = 0_Ki
    Allocate (AA(NS))
    
    Do iL = 1, NS
       AA(iL)%NB_Col = 0_Ki
       Allocate (AA(iL)%Col(Size_Inc))
    End Do

    Do_iE: Do iE = 1, NE
       Do_iSL_1: Do iSL1 = 1, Elem_db(iE)%NB_DoF
          iSG1 = Elem_db(iE)%ID_DoF(iSL1)
          Do_iSLC: Do iSL2 = 1, Elem_db(iE)%NB_DoF
             iSG2 = Elem_db(iE)%ID_DoF(iSL2)
             ! Check if we are in the lower part 
!             Is_Lower: If (iSG1 <= iSG2) Then
                ! Is iSG2 already present in the current Row
                Is_Present = 0
                
                Do_iClo: Do iCol = 1, AA(iSG1)%NB_Col
                   If (AA(iSG1)%Col(iCol)%Pos == iSG2) Then
                      Is_Present = 1
                      ! iSG2 has been found, exit the loop
                      Exit
                   End If
                End Do Do_iClo
                
                If_Is_Present: If (Is_Present == 0) Then
                   NVals = AA(iSG1)%NB_Col
                   If (AA(iSG1)%NB_Col == Size(AA(iSG1)%Col)) Then
                      Allocate(TMP_Row(NVals + Size_Inc))
                      Do iSwap = 1, NVals
                         TMP_Row(iSwap) = AA(iSG1)%Col(iSwap)%Pos
                         If (AA(iSG1)%Col(iSwap)%Val /= 0.0_Kr) Then
                         End If
                      End Do
                      DeAllocate(AA(iSG1)%Col)
                      Allocate(AA(iSG1)%Col(NVals + Size_Inc))
                      Do iSwap = 1, NVals
                         AA(iSG1)%Col(iSwap)%Pos = TMP_Row(iSwap)
                      End Do
                      DeAllocate(TMP_Row)
                   EndIf
                   AA(iSG1)%NB_Col = NVals + 1
                   AA(iSG1)%Col(NVals + 1)%Pos = iSG2
                End If If_Is_Present
 !            End If Is_Lower
          End Do Do_iSLC
       End Do Do_iSL_1
    End Do Do_iE
    Do iS = 1, NS
       NNZ(iS)=AA(iS)%NB_Col
    End Do
  End Subroutine Get_Matrix_Size_2D_Elast

  Subroutine Get_Matrix_Size_3D_Elast_Symm(NNZ, Elem_db, Node_db)
    
    Type Mat_AA_Entry
       Real(Kind = Kr)                            :: Val
       Integer(Kind = Ki)                         :: Pos
    End Type Mat_AA_Entry
    
    Type Mat_AA_Row
       Integer(Kind=Ki)                            :: NB_Col
       Type(Mat_AA_Entry), dimension(:), pointer   :: Col
    End Type Mat_AA_Row
    
    
    Integer, Dimension(:), Pointer                 :: NNZ
    Type(Element3D_Elast), Dimension(:), Pointer   :: Elem_db
    Type(Node3D), Dimension(:), Pointer            :: Node_db
    
    Integer(Kind = Ki)                             :: NS, NE, Nb_DoF
    Integer(Kind = Ki)                             :: iE, iSL1, ISG1
    Integer(Kind = Ki)                             :: iSL2, iSG2, iS, iT
    Type(Mat_AA_Row), Dimension(:), pointer        :: AA
    Integer(Kind = Ki), Dimension(:), Pointer      :: TMP_Row
    Integer(Kind = Ki)                             :: iL
    Integer(Kind = Ki)                             :: NVals
    Integer(Kind = Ki)                             :: iCol
    Integer(Kind = Ki)                             :: iSwap
    Integer(Kind = Ki)                             :: Is_present
    Integer(Kind = Ki)                             :: Size_Inc = 10



    NS = Size(Node_db)
    NE = Size(Elem_db)
    Allocate (NNZ(NS))
    NNZ = 0_Ki
    Allocate (AA(NS))
    
    Do iL = 1, NS
       AA(iL)%NB_Col = 0_Ki
       Allocate (AA(iL)%Col(Size_Inc))
    End Do
    
    Do_iE: Do iE = 1, NE
       Do_iSL_1: Do iSL1 = 1, Elem_db(iE)%NB_DoF
          iSG1 = Elem_db(iE)%ID_DoF(iSL1)
          Do_iSLC: Do iSL2 = 1, Elem_db(iE)%NB_DoF
             iSG2 = Elem_db(iE)%ID_DoF(iSL2)
             ! Check if we are in the lower part 
             Is_Lower: If (iSG1 <= iSG2) Then
                ! Is iSG2 already present in the current Row
                Is_Present = 0
                
                Do_iClo: Do iCol = 1, AA(iSG1)%NB_Col
                   If (AA(iSG1)%Col(iCol)%Pos == iSG2) Then
                      Is_Present = 1
                      ! iSG2 has been found, exit the loop
                      Exit
                   End If
                End Do Do_iClo
                
                If_Is_Present: If (Is_Present == 0) Then
                   NVals = AA(iSG1)%NB_Col
                   If (AA(iSG1)%NB_Col == Size(AA(iSG1)%Col)) Then
                      Allocate(TMP_Row(NVals + Size_Inc))
                      Do iSwap = 1, NVals
                         TMP_Row(iSwap) = AA(iSG1)%Col(iSwap)%Pos
                         If (AA(iSG1)%Col(iSwap)%Val /= 0.0_Kr) Then
                         End If
                      End Do
                      DeAllocate(AA(iSG1)%Col)
                      Allocate(AA(iSG1)%Col(NVals + Size_Inc))
                      Do iSwap = 1, NVals
                         AA(iSG1)%Col(iSwap)%Pos = TMP_Row(iSwap)
                      End Do
                      DeAllocate(TMP_Row)
                   EndIf
                   AA(iSG1)%NB_Col = NVals + 1
                   AA(iSG1)%Col(NVals + 1)%Pos = iSG2
                End If If_Is_Present
             End If Is_Lower
          End Do Do_iSLC
       End Do Do_iSL_1
    End Do Do_iE
    Do iS = 1, NS
       NNZ(iS)=AA(iS)%NB_Col
    End Do
  End Subroutine Get_Matrix_Size_3D_Elast_Symm

  Subroutine Get_Matrix_Size_3D_Elast(NNZ, Elem_db, Node_db)
    
    Type Mat_AA_Entry
       Real(Kind = Kr)                            :: Val
       Integer(Kind = Ki)                         :: Pos
    End Type Mat_AA_Entry
    
    Type Mat_AA_Row
       Integer(Kind=Ki)                            :: NB_Col
       Type(Mat_AA_Entry), dimension(:), pointer   :: Col
    End Type Mat_AA_Row
    
    
    Integer, Dimension(:), Pointer                 :: NNZ
    Type(Element3D_Elast), Dimension(:), Pointer   :: Elem_db
    Type(Node3D), Dimension(:), Pointer            :: Node_db
    
    Integer(Kind = Ki)                             :: NS, NE, Nb_DoF
    Integer(Kind = Ki)                             :: iE, iSL1, ISG1
    Integer(Kind = Ki)                             :: iSL2, iSG2, iS, iT
    Type(Mat_AA_Row), Dimension(:), pointer        :: AA
    Integer(Kind = Ki), Dimension(:), Pointer      :: TMP_Row
    Integer(Kind = Ki)                             :: iL
    Integer(Kind = Ki)                             :: NVals
    Integer(Kind = Ki)                             :: iCol
    Integer(Kind = Ki)                             :: iSwap
    Integer(Kind = Ki)                             :: Is_present
    Integer(Kind = Ki)                             :: Size_Inc = 10



    NS = Size(Node_db)
    NE = Size(Elem_db)
    Allocate (NNZ(NS))
    NNZ = 0_Ki
    Allocate (AA(NS))
    
    Do iL = 1, NS
       AA(iL)%NB_Col = 0_Ki
       Allocate (AA(iL)%Col(Size_Inc))
    End Do

    Do_iE: Do iE = 1, NE
       Do_iSL_1: Do iSL1 = 1, Elem_db(iE)%NB_DoF
          iSG1 = Elem_db(iE)%ID_DoF(iSL1)
          Do_iSLC: Do iSL2 = 1, Elem_db(iE)%NB_DoF
             iSG2 = Elem_db(iE)%ID_DoF(iSL2)
             ! Check if we are in the lower part 
!             Is_Lower: If (iSG1 <= iSG2) Then
                ! Is iSG2 already present in the current Row
                Is_Present = 0
                
                Do_iClo: Do iCol = 1, AA(iSG1)%NB_Col
                   If (AA(iSG1)%Col(iCol)%Pos == iSG2) Then
                      Is_Present = 1
                      ! iSG2 has been found, exit the loop
                      Exit
                   End If
                End Do Do_iClo
                
                If_Is_Present: If (Is_Present == 0) Then
                   NVals = AA(iSG1)%NB_Col
                   If (AA(iSG1)%NB_Col == Size(AA(iSG1)%Col)) Then
                      Allocate(TMP_Row(NVals + Size_Inc))
                      Do iSwap = 1, NVals
                         TMP_Row(iSwap) = AA(iSG1)%Col(iSwap)%Pos
                         If (AA(iSG1)%Col(iSwap)%Val /= 0.0_Kr) Then
                         End If
                      End Do
                      DeAllocate(AA(iSG1)%Col)
                      Allocate(AA(iSG1)%Col(NVals + Size_Inc))
                      Do iSwap = 1, NVals
                         AA(iSG1)%Col(iSwap)%Pos = TMP_Row(iSwap)
                      End Do
                      DeAllocate(TMP_Row)
                   EndIf
                   AA(iSG1)%NB_Col = NVals + 1
                   AA(iSG1)%Col(NVals + 1)%Pos = iSG2
                End If If_Is_Present
 !            End If Is_Lower
          End Do Do_iSLC
       End Do Do_iSL_1
    End Do Do_iE
    Do iS = 1, NS
       NNZ(iS)=AA(iS)%NB_Col
    End Do
  End Subroutine Get_Matrix_Size_3D_Elast

  Subroutine Get_Matrix_Size_3D_Scal_Symm(NNZ, Elem_db, Node_db)
    
    Type Mat_AA_Entry
       Real(Kind = Kr)                            :: Val
       Integer(Kind = Ki)                         :: Pos
    End Type Mat_AA_Entry
    
    Type Mat_AA_Row
       Integer(Kind=Ki)                            :: NB_Col
       Type(Mat_AA_Entry), dimension(:), pointer   :: Col
    End Type Mat_AA_Row
    
    
    Integer, Dimension(:), Pointer                 :: NNZ
    Type(Element3D_Scal), Dimension(:), Pointer   :: Elem_db
    Type(Node3D), Dimension(:), Pointer            :: Node_db
    
    Integer(Kind = Ki)                             :: NS, NE, Nb_DoF
    Integer(Kind = Ki)                             :: iE, iSL1, ISG1
    Integer(Kind = Ki)                             :: iSL2, iSG2, iS, iT
    Type(Mat_AA_Row), Dimension(:), pointer        :: AA
    Integer(Kind = Ki), Dimension(:), Pointer      :: TMP_Row
    Integer(Kind = Ki)                             :: iL
    Integer(Kind = Ki)                             :: NVals
    Integer(Kind = Ki)                             :: iCol
    Integer(Kind = Ki)                             :: iSwap
    Integer(Kind = Ki)                             :: Is_present
    Integer(Kind = Ki)                             :: Size_Inc = 10



    NS = Size(Node_db)
    NE = Size(Elem_db)
    Allocate (NNZ(NS))
    NNZ = 0_Ki
    Allocate (AA(NS))
    
    Do iL = 1, NS
       AA(iL)%NB_Col = 0_Ki
       Allocate (AA(iL)%Col(Size_Inc))
    End Do
    
    Do_iE: Do iE = 1, NE
       Do_iSL_1: Do iSL1 = 1, Elem_db(iE)%NB_DoF
          iSG1 = Elem_db(iE)%ID_DoF(iSL1)
          Do_iSLC: Do iSL2 = 1, Elem_db(iE)%NB_DoF
             iSG2 = Elem_db(iE)%ID_DoF(iSL2)
             ! Check if we are in the lower part 
             Is_Lower: If (iSG1 <= iSG2) Then
                ! Is iSG2 already present in the current Row
                Is_Present = 0
                
                Do_iClo: Do iCol = 1, AA(iSG1)%NB_Col
                   If (AA(iSG1)%Col(iCol)%Pos == iSG2) Then
                      Is_Present = 1
                      ! iSG2 has been found, exit the loop
                      Exit
                   End If
                End Do Do_iClo
                
                If_Is_Present: If (Is_Present == 0) Then
                   NVals = AA(iSG1)%NB_Col
                   If (AA(iSG1)%NB_Col == Size(AA(iSG1)%Col)) Then
                      Allocate(TMP_Row(NVals + Size_Inc))
                      Do iSwap = 1, NVals
                         TMP_Row(iSwap) = AA(iSG1)%Col(iSwap)%Pos
                         If (AA(iSG1)%Col(iSwap)%Val /= 0.0_Kr) Then
                         End If
                      End Do
                      DeAllocate(AA(iSG1)%Col)
                      Allocate(AA(iSG1)%Col(NVals + Size_Inc))
                      Do iSwap = 1, NVals
                         AA(iSG1)%Col(iSwap)%Pos = TMP_Row(iSwap)
                      End Do
                      DeAllocate(TMP_Row)
                   EndIf
                   AA(iSG1)%NB_Col = NVals + 1
                   AA(iSG1)%Col(NVals + 1)%Pos = iSG2
                End If If_Is_Present
             End If Is_Lower
          End Do Do_iSLC
       End Do Do_iSL_1
    End Do Do_iE
    Do iS = 1, NS
       NNZ(iS)=AA(iS)%NB_Col
    End Do
  End Subroutine Get_Matrix_Size_3D_Scal_Symm

  Subroutine Get_Matrix_Size_3D_Scal(NNZ, Elem_db, Node_db)
    
    Type Mat_AA_Entry
       Real(Kind = Kr)                            :: Val
       Integer(Kind = Ki)                         :: Pos
    End Type Mat_AA_Entry
    
    Type Mat_AA_Row
       Integer(Kind=Ki)                            :: NB_Col
       Type(Mat_AA_Entry), dimension(:), pointer   :: Col
    End Type Mat_AA_Row
    
    
    Integer, Dimension(:), Pointer                 :: NNZ
    Type(Element3D_Scal), Dimension(:), Pointer   :: Elem_db
    Type(Node3D), Dimension(:), Pointer            :: Node_db
    
    Integer(Kind = Ki)                             :: NS, NE, Nb_DoF
    Integer(Kind = Ki)                             :: iE, iSL1, ISG1
    Integer(Kind = Ki)                             :: iSL2, iSG2, iS, iT
    Type(Mat_AA_Row), Dimension(:), pointer        :: AA
    Integer(Kind = Ki), Dimension(:), Pointer      :: TMP_Row
    Integer(Kind = Ki)                             :: iL
    Integer(Kind = Ki)                             :: NVals
    Integer(Kind = Ki)                             :: iCol
    Integer(Kind = Ki)                             :: iSwap
    Integer(Kind = Ki)                             :: Is_present
    Integer(Kind = Ki)                             :: Size_Inc = 10



    NS = Size(Node_db)
    NE = Size(Elem_db)
    Allocate (NNZ(NS))
    NNZ = 0_Ki
    Allocate (AA(NS))
    
    Do iL = 1, NS
       AA(iL)%NB_Col = 0_Ki
       Allocate (AA(iL)%Col(Size_Inc))
    End Do

    Do_iE: Do iE = 1, NE
       Do_iSL_1: Do iSL1 = 1, Elem_db(iE)%NB_DoF
          iSG1 = Elem_db(iE)%ID_DoF(iSL1)
          Do_iSLC: Do iSL2 = 1, Elem_db(iE)%NB_DoF
             iSG2 = Elem_db(iE)%ID_DoF(iSL2)
             ! Check if we are in the lower part 
!             Is_Lower: If (iSG1 <= iSG2) Then
                ! Is iSG2 already present in the current Row
                Is_Present = 0
                
                Do_iClo: Do iCol = 1, AA(iSG1)%NB_Col
                   If (AA(iSG1)%Col(iCol)%Pos == iSG2) Then
                      Is_Present = 1
                      ! iSG2 has been found, exit the loop
                      Exit
                   End If
                End Do Do_iClo
                
                If_Is_Present: If (Is_Present == 0) Then
                   NVals = AA(iSG1)%NB_Col
                   If (AA(iSG1)%NB_Col == Size(AA(iSG1)%Col)) Then
                      Allocate(TMP_Row(NVals + Size_Inc))
                      Do iSwap = 1, NVals
                         TMP_Row(iSwap) = AA(iSG1)%Col(iSwap)%Pos
                         If (AA(iSG1)%Col(iSwap)%Val /= 0.0_Kr) Then
                         End If
                      End Do
                      DeAllocate(AA(iSG1)%Col)
                      Allocate(AA(iSG1)%Col(NVals + Size_Inc))
                      Do iSwap = 1, NVals
                         AA(iSG1)%Col(iSwap)%Pos = TMP_Row(iSwap)
                      End Do
                      DeAllocate(TMP_Row)
                   EndIf
                   AA(iSG1)%NB_Col = NVals + 1
                   AA(iSG1)%Col(NVals + 1)%Pos = iSG2
                End If If_Is_Present
 !            End If Is_Lower
          End Do Do_iSLC
       End Do Do_iSL_1
    End Do Do_iE
    Do iS = 1, NS
       NNZ(iS)=AA(iS)%NB_Col
    End Do
  End Subroutine Get_Matrix_Size_3D_Scal


  Subroutine Ext_Val_I(V_In, V_Out, Val)
    Integer(Kind = Ki), Dimension(:), Pointer   :: V_In
    Integer(Kind = Ki), Dimension(:), Pointer   :: V_Out
    Integer(Kind = Ki), Intent(IN)              :: Val

    Integer(Kind = Ki)                          :: Nb_Ext, iExt, i

!     Nullify(V_out)
    If (Associated(V_OUT)) Then
       DeAllocate(V_Out)
    EndIf

    Nb_Ext = Count(V_In == Val)
    Allocate(V_Out(Nb_Ext))
    iExt = 0
    Do i=1, Size(V_In)
       If ( V_In(i) == Val) Then
          iExt = iExt+1
          V_Out(iExt) = i
       EndIf
    EndDo
  End Subroutine Ext_Val_I



End Module m_MEF_Tools
