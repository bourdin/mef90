Program Prep_Poisson
  Use m_MEF90
  Use m_Poisson_Struct
  Implicit NONE

  Type(EXO_Geom_Info)                             :: Geom
  Integer                                         :: iS, iErr
  Type(Node2D), Dimension(:), Pointer             :: Node2D_db
  Type (Element2D_Scal), Dimension(:), Pointer    :: Elem2D_db
  Type(Node3D), Dimension(:), Pointer             :: Node3D_db
  Type (Element3D_Scal), Dimension(:), Pointer    :: Elem3D_db
  Type (Poisson_Params)                           :: Params
  Character(len = MXSTLN)                         :: EXO_Str
  Character(len=MXSTLN), Dimension(3)             :: Coord_Names

  Write(*,100, advance = 'no') 'ExodusII file name: '
  Read(*,100) EXO_Str

  Geom%filename = Trim(EXO_Str)//'.gen'
  Call Read_EXO_Geom_Info(Geom)

  Call Ask_Poisson_Params(Geom, Params)


  
  Select Case (Params%PB_Dim)
  Case (PB_2D)
     Call Read_EXO_Node_Coord(Geom, Node2D_db, 1)
     Call Read_EXO_Connect(Geom, Elem2D_db)
  Case (PB_3D)
     Call Read_EXO_Node_Coord(Geom, Node3D_db, 1)
     Call Read_EXO_Connect(Geom, Elem3D_db)
  Case Default
     Print*, 'Unknow Params%PB_Dim ', Params%PB_Dim
     STOP
  End Select

  Write(*,*) 'Updating geometrical information...'
!  Geom%exoid = EXCRE (Geom%filename, EXCLOB, exo_cpu_ws, exo_io_ws,iErr)
  Params%Sim_Str   = Trim(EXO_Str) // '_cst'
  Geom%filename    = Trim(Params%Sim_Str) // '.gen'
  Params%PARAM_Str = Trim(Params%Sim_Str) // '.PARAM'
  Call Write_EXO_Geom_Info(Geom)

  Select Case (Params%PB_Dim)
  Case (PB_2D)
     Write(*,*) 'Updating coordinates...'
     Call Write_EXO_Node_Coord(Geom, Node2D_db, 1)
     Write(*,*) 'Updating connectivity table...'
     Call Write_EXO_Connect(Geom, Elem2D_db)
  Case (PB_3D)
     Write(*,*) 'Updating coordinates...'
     Call Write_EXO_Node_Coord(Geom, Node3D_db, 1)
     Write(*,*) 'Updating connectivity table...'
     Call Write_EXO_Connect(Geom, Elem3D_db)
  End Select

  Write(*,*) 'Adding problem parameters...'
  Call Write_Poisson_EXO_Params(Geom, Params)

  Select Case (Params%PB_Dim)
  Case (PB_2D)
     Write(*,*) 'Adding Rhs...'
     Call Write_F2D(Geom, Params, Node2D_db)
  Case (PB_3D)
     Write(*,*) 'Adding Rhs...'
     Call Write_F3D(Geom, Params, Node3D_db)
  End Select

  Write(*,*) 'Writing problem parameters'
  Call Write_Poisson_DATA(Geom, Params)

  Write(*,*) 'DONE!'


100 Format(A)

Contains
  Subroutine Write_F2D(Geom, Params, Node_db)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Poisson_Params)                         :: Params
    Type (Node2D), Dimension(:), Pointer          :: Node_db
    
    Integer                                       :: iTS, iNode, iSet
    Real(Kind = Kr), Dimension(:), Pointer        :: Rhs_Result 

    Allocate (Rhs_Result(Geom%Num_Nodes))
    Rhs_Result=1.0_Kr

    Call Write_EXO_Result_Nodes(Geom, 1, 1, Rhs_Result)
    DeAllocate (Rhs_Result)
  End Subroutine Write_F2D

  Subroutine Write_F3D(Geom, Params, Node_db)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Poisson_Params)                         :: Params
    Type (Node3D), Dimension(:), Pointer          :: Node_db
    
    Real(Kind = Kr), Dimension(:), Pointer        :: Rhs_Result 
    
    Allocate (Rhs_Result(Geom%Num_Nodes))
    Rhs_Result=1.0_Kr 

    Call Write_EXO_Result_Nodes(Geom, 1, 1, Rhs_Result)
    DeAllocate (Rhs_Result)
  End Subroutine Write_F3D
  
  Subroutine Ask_Poisson_Params(Geom, Params)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Poisson_Params)                         :: Params
    Integer                                       :: i, iBlk, iSet, Set_ID
    Integer                                       :: Num_TS
    Real(Kind = Kr)                               :: TS_init, Tmp_BC
    
    Write(*,100) '===== Global Properties ======'
    Num_TS =1
    Params%Load = TS_INIT
    
    Write(*,100) '===== Node sets Properties ========='
    Allocate (Params%BC_Type(Geom%Num_Node_Sets))
    Allocate (Params%BC_DB(Geom%Num_Node_Sets)) 
    Do iSet = 1, Geom%Num_Node_Sets
       Params%BC_DB(iSet)%ID = Geom%Node_Set(iSet)%ID
       Allocate(Params%BC_DB(iSet)%Dis_F(Geom%Node_Set(iSet)%Num_Nodes))
       Write(*,400) Params%BC_DB(iSet)%ID
       Write(*,100, advance = 'no') 'BC type (None=0 DIRI = 1) '
       
       Params%BC_DB(iSet)%Dis_F = 0.0_Kr
       Read(*,*) Params%BC_Type(iSet) 
       If (Params%BC_Type(iSet) == 1) Then
          Write(*,100, advance = 'no') 'BC                    '
          Read(*,*) Tmp_BC
          Params%BC_DB(iSet)%Dis_F = Tmp_BC
       End If
       
    End Do

    Write(*,100) '===== Input Parameters ======'
    
    Write(*,100, advance = 'no') 'PB_Dim   [0=1D 1=2DA 2=2D 3=3D]   :'
    Read(*, 101) Params%PB_Dim
    Write(*,100, advance = 'no') 'Init_U   [0=PREV 1=0]   :'
    Read(*, 101) Params%Init_U
    Write(*,100, advance = 'no') 'MaxIterRelax    :'
    Read(*, 101) Params%MaxIterRelax
    Write(*,100, advance = 'no') 'TolRelax     :'
    Read(*, 102) Params%TolRelax
    Write(*,100, advance = 'no') 'TolKSP      :'
    Read(*, 102) Params%TolKSP
    

101 Format(I4)
102 Format(F12.5)

100 Format(A)
400 Format('=== Node set                  ', I3)
  End Subroutine Ask_Poisson_Params

End Program Prep_Poisson
