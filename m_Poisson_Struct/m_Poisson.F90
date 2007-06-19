Module m_Poisson_Struct
  Use m_MEF90
  Implicit NONE
  Private

#include "include/finclude/petsc.h"

  Public :: Write_Poisson_EXO_Params
  Public :: Read_Poisson_EXO_Params
  Public :: Read_Poisson_DATA
  Public :: Write_Poisson_DATA
  Public :: Broadcast_Poisson_Params

  Public :: Poisson_Params


  
  Integer, Parameter, Public                      :: BC_Type_NONE = 0
  Integer, Parameter, Public                      :: BC_Type_DIRI = 1

  Integer, Parameter, Public                      :: PB_2D  = 2
  Integer, Parameter, Public                      :: PB_3D  = 3

  Integer, Parameter, Public                      :: Num_Prop_NS  = 1
  Character(len=MXSTLN), Parameter  :: Prop_Name_NS = 'BC_Type  '

  Integer, Parameter                              :: Num_Res_G = 1
  Character(len=MXSTLN), Parameter  :: Res_Name_G = 'Bulk Energy    '

  Integer, Parameter                              :: Num_Res_N  = 3
  Character(len=MXSTLN), Dimension(3), Parameter  :: Res_Name_N =             &
      &      (/ 'Load      ',                                                 &
      &         'Solution  ',                                                 &
      &         'Node Owner' /)

  Integer, Parameter                              :: Num_Res_E  = 1
  Character(len=MXSTLN), Dimension(1), Parameter  :: Res_Name_E =             &
      &      (/ 'Element Owner'  /)

  Type Nd_DB
     Integer                                        :: ID
     Real(Kind = Kr), Dimension(:), Pointer         :: Dis_F
  End Type Nd_DB
  Type Poisson_Params
     !!! STORED IN THE EXODUS DB
     ! GLOBAL PROPERTIES
          
     Character(len = MXLNLN)                      :: Sim_Str
     Character(len = MXLNLN)                      :: PARAM_Str
     Character(len = MXLNLN)                      :: CST_Str

     
     ! NODE SETS PROPERTIES
     Integer, Dimension(:), Pointer               :: BC_Type
     Type(Nd_DB), Dimension(:), Pointer         :: BC_DB

     ! GLOBAL PARAMETERS (STORED AS RESULTS)
     Real(Kind = Kr), Dimension(:), Pointer       :: Load


     !!! STORED IN THE PARAM FILE
     ! GLOBAL DATAS
     ! The following parameter are NOT STORED in the Exodus db
     Integer                                      :: PB_Dim
     Integer                                      :: Init_U
     Integer                                      :: MaxIterRelax
     Real(Kind = Kr)                              :: TolRelax
     Real(Kind = Kr)                              :: TolKSP

  End Type Poisson_Params

Contains

  Subroutine Write_Poisson_EXO_Params(Geom, Params)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Poisson_Params)                         :: Params
    
    Character(len=MXSTLN), Dimension(:,:), Pointer:: Tmp_QA
    Integer                                       :: iErr
    Integer                                       :: iTs, iS, iSet
    Integer                                       :: exo_ver

    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, exo_ver,&
         & ierr)

    ! Object Properties
    Call EXPPN(Geom%exoid, EXNSET, Num_Prop_NS, Prop_Name_NS, iErr)
   Call EXPPA(Geom%exoid, EXNSET, Prop_Name_NS, Params%BC_Type, iErr)

    Do iSet = 1, Geom%Num_Node_Sets
    Call EXPNSD(Geom%exoid, Params%BC_DB(iSet)%ID, Params%BC_DB(iSet)%Dis_F, iErr)
    End Do   

    !Result and initial values names
    Call EXPVP (Geom%exoid, 'g', Num_Res_G, iErr)
    Call EXPVAN (Geom%exoid, 'g', Num_Res_G, Res_Name_G, iErr)
    Call EXPVP (Geom%exoid, 'n', Num_Res_N, iErr)
    Call EXPVAN (Geom%exoid, 'n', Num_Res_N, Res_Name_N, iErr)
    Call EXPVP (Geom%exoid, 'e', Num_Res_E, iErr)
    Call EXPVAN (Geom%exoid, 'e', Num_Res_E, Res_Name_E, iErr)
    
    Do iTS = 1, Size(Params%Load)
       Call EXPTIM(Geom%exoid, iTS, real(iTS, 8), iErr)
       Call EXPGV (Geom%exoid, iTS, 4, (/ 0.0_Kr, 0.0_Kr, 0.0_Kr,          &
            & Params%Load(iTS)/), iErr) ! This needs to be fixed
    End Do    

    Geom%Num_QA = Geom%Num_QA+1
    Allocate (Tmp_QA(4, Geom%Num_QA))
    Tmp_QA(:,1:Geom%Num_QA-1) = Geom%QA_rec
    DeAllocate (Geom%QA_Rec)
    Allocate (Geom%QA_Rec(4, Geom%Num_QA))
    Geom%QA_rec = Tmp_QA
    DeAllocate (Tmp_QA)
   
!    Geom%QA_Rec(1,Geom%Num_QA) = 'm_Poisson-Struct'
!    Geom%QA_Rec(2,Geom%Num_QA) = '0.0.1'
!    Call Date_And_Time(date = Geom%QA_Rec(3,Geom%Num_QA))
!    Call Date_And_Time(time = Geom%QA_Rec(4,Geom%Num_QA))
!    Call EXPQA(Geom%exoid, Geom%num_QA, Geom%QA_Rec, iErr)

    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_Poisson_EXO_Params

  Subroutine Read_Poisson_EXO_Params(Geom, Params)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Poisson_Params)                            :: Params
    
    Character(len=MXSTLN), Dimension(:,:), Pointer:: Tmp_QA
    Integer                                       :: iErr
    Integer                                       :: iTs, iS, iBlk, iSet
    Integer                                       :: Num_Vars, Num_TS
    Integer                                       :: exo_ver
    Real(Kind = Kr), Dimension(:), Pointer         :: Tmp_Res
    Real(Kind = Kr)                               :: fDum
    Character                                     :: cDum
  

    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, exo_ver,&
         & ierr)
    
    Allocate(Params%BC_Type(Geom%Num_Node_Sets))
    Allocate(Params%BC_DB(Geom%Num_Node_Sets))
    Do iSet = 1, Geom%Num_Node_Sets
       Call EXGP(Geom%exoid, EXNSET, iSet, Prop_Name_NS,                   &
            & Params%BC_Type(iSet), iErr)
       Allocate(Params%BC_DB(iSet)%Dis_F(Geom%Node_Set(iSet)%Num_Nodes))
       Call EXGNSD(Geom%exoid, iSet, Params%BC_DB(iSet)%Dis_F, iErr)
    End Do
    ! Time Step
    ! Time has to be monotonically increasing, so we store the number 
    ! of the time step instead, and store the displacement factor or the 
    ! temperature as a global result
    Call EXINQ (Geom%exoid, EXTIMS, iTS, fDum, cDum, iErr)    
    Allocate (Params%Load(iTS))
    !! We need to check this later
    Allocate(Tmp_Res(iTs))
    Call EXGVP(Geom%exoid, 'G', Num_Vars, iErr)
    Call EXGGV(Geom%exoid, iTs, Num_Vars, Tmp_Res, iErr)
    Params%Load(iTS) = Tmp_Res(4)
    DeAllocate (Tmp_Res)
    
    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0

!    Allocate (Params%Load(iTS))
!    Do iTS = 1, Size(Params%Load)
!       Call Read_EXO_Result_Global(Geom, iTS, 4, Params%Load(iTS))
!    End Do
!!   Call Read_EXO_Result_Global(Geom, iTS, 4, Params%Load)
!  Call EXCLOS(Geom%exoid, iErr)
!    Geom%exoid = 0
  End Subroutine Read_Poisson_EXO_Params


  Subroutine Write_Poisson_DATA(Geom, Params)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Poisson_Params)                            :: Params
    
    Integer                                       :: iBlk, Blk_ID

    Open(File = Params%PARAM_Str, Unit = F_OUT, Status = 'Unknown')
    Rewind(F_OUT)

    Write(F_OUT, 100) Params%PB_Dim,       'PB_Dim   [0=1D 1=2DA 2=2D 3=3D]'
    Write(F_OUT, 100) Params%Init_U,       'Init_U   [0=PREV 1=0]'
    Write(F_OUT, 100) Params%MaxIterRelax, 'MaxIterRelax'
    Write(F_OUT, 102) Params%TolRelax,     'TolRelax'
    Write(F_OUT, 102) Params%TolKSP,        'TolKSP'
    Close(F_OUT)


100 Format(I4,T16,'# ', A)
102 Format(ES12.5,T16,'# ', A)
110 Format(I4,T16,'# BLK_ID, Toughness, E, nu')
120 Format(I4, 3(ES12.5,' '))
  End Subroutine Write_Poisson_DATA

  Subroutine Read_Poisson_DATA(Geom, Params)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Poisson_Params)                            :: Params
    
    Integer                                       :: iBlk, Blk_ID, Num_Blk

    Open(File = Params%PARAM_Str, Unit = F_IN, Status = 'Old')
    Rewind(F_IN)

    Read(F_IN, 100) Params%PB_Dim
    Read(F_IN, 100) Params%Init_U
    Read(F_IN, 100) Params%MaxIterRelax
    Read(F_IN, 102) Params%TolRelax
    Read(F_IN, 102) Params%TolKSP
    Close(F_IN)
    

100 Format(I4)
102 Format(ES12.5)
  End Subroutine Read_Poisson_DATA

  Subroutine Broadcast_Poisson_Params(Params, Source, Communicator)
    Type (Poisson_Params), Intent(INOUT)                  :: Params
    Integer, Intent(IN)                                :: Source
    Integer, Intent(IN), Optional                      :: Communicator

    Integer                                            :: MyRank, iErr, Comm
    Integer                                            :: iBlk, iSet, iLoad
    Integer                                            :: nBlks, nSets, nLoads

    If ( Present(Communicator) ) Then
       Comm = Communicator
    Else
       Comm = MPI_COMM_WORLD
    End If

    Call MPI_COMM_RANK(Comm, MyRank, iErr)

    Call MPI_BCAST(Params%Sim_Str, MXSTLN, MPI_CHARACTER, Source, Comm, iErr)
    Call MPI_BCAST(Params%PARAM_Str, MXSTLN, MPI_CHARACTER, Source, Comm, iErr)

    If (MyRank == Source) Then
       nSets = Size(Params%BC_Type)
       nLoads = Size(Params%Load)
    End If
    Call MPI_BCAST(nSets, 1, MPI_INTEGER, Source, Comm, iErr)
    Call MPI_BCAST(nLoads, 1, MPI_INTEGER, Source, Comm, iErr)
    
    If (MyRank /= Source) Then
       If ( Associated (Params%BC_Type) ) Then
          DeAllocate(Params%BC_Type)
       End If
       If ( Associated (Params%Load) ) Then
          DeAllocate(Params%Load)
       End If
       Allocate(Params%BC_Type(nSets))      
       Allocate(Params%Load(nLoads))
    End If
    Call MPI_BCAST(Params%BC_Type, nSets, MPI_INTEGER, Source, Comm, iErr)
    Call MPI_BCAST(Params%Load, nLoads, MPI_DOUBLE_PRECISION, Source, Comm,   &
         & iErr)

    Call MPI_BCAST(Params%PB_Dim, 1, MPI_INTEGER, Source, Comm, iErr)
    Call MPI_BCAST(Params%Init_U, 1, MPI_INTEGER, Source, Comm, iErr)
    
    Call MPI_BCAST(Params%MaxIterRelax, 1, MPI_INTEGER, Source, Comm, iErr)
    Call MPI_BCAST(Params%TolRelax, 1, MPI_DOUBLE_PRECISION, Source, Comm,    &
         & iErr)
    Call MPI_BCAST(Params%TolKSP, 1, MPI_DOUBLE_PRECISION, Source, Comm,      &
         & iErr)
  End Subroutine Broadcast_Poisson_Params


End Module m_Poisson_Struct
