Module m_Rupt_Struct
  Use m_MEF90
  Implicit NONE
  Private

#include "include/finclude/petsc.h"

  Public :: Read_Rupt_EXO_Params
  Public :: Write_Rupt_EXO_Params
  Public :: Read_Rupt_DATA
  Public :: Write_Rupt_DATA

  Public :: Rupt_Params
  !!! OBSOLETE WILL BE REMOVED IN THE NEXT VERSION
  
  Public :: Rupt_Params2D
  Public :: Rupt_Params3D
  
  Public :: GenHL_Iso_LambdaMu
  Public :: GenHL_Iso3D_Enu
  Public :: GenHL_Iso2D_EnuPlaneStress
  Public :: GenHL_Iso2D_EnuPlaneStrain
  Public :: GenHL_Ortho2D_LambdaMu
  
  Interface Read_Rupt_EXO_Params
      Module Procedure Read_Rupt_EXO_Params_Iso, Read_Rupt_EXO_Params2D, Read_Rupt_EXO_Params3D
  End Interface

  Interface Write_Rupt_EXO_Params
      Module Procedure Write_Rupt_EXO_Params_Iso, Write_Rupt_EXO_Params2D, Write_Rupt_EXO_Params3D
  End Interface

  Interface Read_Rupt_DATA
      Module Procedure Read_Rupt_DATA_Iso, Read_Rupt_DATA2D, Read_Rupt_DATA3D
  End Interface
  
  Interface Write_Rupt_DATA
      Module Procedure Write_Rupt_DATA_Iso, Write_Rupt_DATA2D, Write_Rupt_DATA3D
  End Interface
  
  Interface GenHL_Iso_LambdaMu
      Module Procedure GenHL_Iso2D_LambdaMu, GenHL_Iso3D_LambdaMu
  End Interface

  Integer, Parameter, Public                      :: Init_V_PREV = 0
  Integer, Parameter, Public                      :: Init_V_ONE = 1
  Integer, Parameter, Public                      :: Init_V_RND = 2
  Integer, Parameter, Public                      :: Init_V_SPH = 3

  Integer, Parameter, Public                      :: Init_U_PREV = 0
  Integer, Parameter, Public                      :: Init_U_ZERO = 1
  
  Integer, Parameter, Public                      :: BC_Type_NONE = 0
  Integer, Parameter, Public                      :: BC_Type_DIRI = 1

  Integer, Parameter                              :: Num_Prop_EB  = 3
  Character(len=MXSTLN), Dimension(3), Parameter  :: Prop_Name_EB =           &
       &     (/ 'Is_Brittle',                                                 &
       &        'Is_Domain ',                                                 &
       &        'Has_Force ' /)

  Integer, Parameter                              :: Num_Prop_NS  = 3
  Character(len=MXSTLN), Dimension(3), Parameter  :: Prop_Name_NS =           &
       &     (/ 'BC_Type_X ',                                                 &
       &        'BC_Type_Y ',                                                 &
       &        'BC_Type_Z ' /)

  Integer, Parameter                              :: Num_Res_G = 4
  Character(len=MXSTLN), Dimension(4), Parameter  :: Res_Name_G =             &
       &     (/ 'Bulk energy   ',                                             &
       &        'Surface energy',                                             &
       &        'Total energy  ' ,                                            &
       &        'Load          '/)

  Integer, Parameter                              :: Num_Res_N  = 8
  Character(len=MXSTLN), Dimension(8), Parameter  :: Res_Name_N =             &
      &      (/ 'Fracture      ',                                             &
      &         'Displacement X',                                             &
      &         'Displacement Y',                                             &
      &         'Displacement Z',                                             &
      &         'Force X       ',                                             &
      &         'Force Y       ',                                             &
      &         'Force Z       ',                                             &
      &         'Temperature   ' /)

  Integer, Parameter                              :: Num_Res_E  = 12
  Character(len=MXSTLN), Dimension(12), Parameter :: Res_Name_E =             &
      &      (/ 'Strain XX',                                                  &
      &         'Strain YY',                                                  &
      &         'Strain ZZ',                                                  &
      &         'Strain XY',                                                  &
      &         'Strain YZ',                                                  &
      &         'Strain XZ',                                                  &
      &         'Stress XX',                                                  &
      &         'Stress YY',                                                  &
      &         'Stress ZZ',                                                  &
      &         'Stress XY',                                                  &
      &         'Stress YZ',                                                  &
      &         'Stress ZX' /)

  Type Rupt_Params
     !!! STORED IN THE EXODUS DB
     ! GLOBAL PROPERTIES
          
     Character(len = MXLNLN)                      :: Sim_Str
     Character(len = MXLNLN)                      :: PARAM_Str
     Character(len = MXLNLN)                      :: CST_Str

     ! ELEMENT BLOCK PROPERTIES
     Logical, Dimension(:), Pointer               :: Is_Brittle
     Logical, Dimension(:), Pointer               :: Is_Domain
     Logical, Dimension(:), Pointer               :: Has_Force
     
     ! NODE SETS PROPERTIES
     Integer, Dimension(:), Pointer               :: BC_Type_X
     Integer, Dimension(:), Pointer               :: BC_Type_Y
     Integer, Dimension(:), Pointer               :: BC_Type_Z

     ! GLOBAL PARAMETERS (STORED AS RESULTS)
     Real(Kind = Kr), Dimension(:), Pointer       :: Load

     ! ELEMENT PARAMETERS (STORED AS RESULTS)

     ! NODES PROPERTIES (STORES AS RESULTS)

     !!! STORED IN THE PARAM FILE
     ! GLOBAL DATAS
     ! The following parameter are NOT STORED in the Exodus db
!     Integer                                      :: PB_Type
!     Integer                                      :: PB_Dim
     Logical                                      :: Do_Irrev
     
     Logical                                      :: Do_BackTrack
     Real(Kind = Kr)                              :: Tol_Ener
               
     Integer                                      :: Init_U
     
     Integer                                      :: Init_V
     Integer                                      :: nbCracks
     Real(Kind = Kr)                              :: MaxCrackLength     

     Integer                                      :: MaxIterRelax
     Real(Kind = Kr)                              :: TolRelax
     Real(Kind = Kr)                              :: TolKSP

     Real(Kind = Kr)                              :: Epsilon
     Real(Kind = Kr)                              :: KEpsilon
     Real(Kind = Kr)                              :: TolIrrev

     !!! STORED IN THE CST FILE
     ! ELEMENT BLOCKS PARAMETERS
     Real(Kind = Kr), Dimension(:), Pointer       :: Toughness
     Real(Kind = Kr), Dimension(:), Pointer       :: Young_Mod
     Real(Kind = Kr), Dimension(:), Pointer       :: Poisson_Ratio
     Real(Kind = Kr), Dimension(:), Pointer       :: Therm_Exp
  End Type Rupt_Params

  Type Rupt_Params2D
     !!! STORED IN THE EXODUS DB
     ! GLOBAL PROPERTIES
          
     Character(len = MXLNLN)                      :: Sim_Str
     Character(len = MXLNLN)                      :: PARAM_Str
     Character(len = MXLNLN)                      :: CST_Str

     ! ELEMENT BLOCK PROPERTIES
     Logical, Dimension(:), Pointer               :: Is_Brittle
     Logical, Dimension(:), Pointer               :: Is_Domain
     Logical, Dimension(:), Pointer               :: Has_Force
     
     ! NODE SETS PROPERTIES
     Integer, Dimension(:), Pointer               :: BC_Type_X
     Integer, Dimension(:), Pointer               :: BC_Type_Y
     Integer, Dimension(:), Pointer               :: BC_Type_Z

     ! GLOBAL PARAMETERS (STORED AS RESULTS)
     Real(Kind = Kr), Dimension(:), Pointer       :: Load

     ! ELEMENT PARAMETERS (STORED AS RESULTS)

     ! NODES PROPERTIES (STORES AS RESULTS)

     !!! STORED IN THE PARAM FILE
     ! GLOBAL DATAS
     ! The following parameter are NOT STORED in the Exodus db
     Logical                                      :: Do_Irrev
     
     Logical                                      :: Do_BackTrack
     Real(Kind = Kr)                              :: Tol_Ener
               
     Integer                                      :: Init_U
     
     Integer                                      :: Init_V
     Integer                                      :: nbCracks
     Real(Kind = Kr)                              :: MaxCrackLength     

     Integer                                      :: MaxIterRelax
     Real(Kind = Kr)                              :: TolRelax
     Real(Kind = Kr)                              :: TolKSP

     Real(Kind = Kr)                              :: Epsilon
     Real(Kind = Kr)                              :: KEpsilon
     Real(Kind = Kr)                              :: TolIrrev

     !!! STORED IN THE CST FILE
     ! ELEMENT BLOCKS PARAMETERS
     Real(Kind = Kr), Dimension(:), Pointer       :: Toughness
     Type(Tens4OS2D), Dimension(:), Pointer       :: Hookes_Law
     Real(Kind = Kr), Dimension(:), Pointer       :: Therm_Exp
  End Type Rupt_Params2D

  Type Rupt_Params3D
     !!! STORED IN THE EXODUS DB
     ! GLOBAL PROPERTIES
          
     Character(len = MXLNLN)                      :: Sim_Str
     Character(len = MXLNLN)                      :: PARAM_Str
     Character(len = MXLNLN)                      :: CST_Str

     ! ELEMENT BLOCK PROPERTIES
     Logical, Dimension(:), Pointer               :: Is_Brittle
     Logical, Dimension(:), Pointer               :: Is_Domain
     Logical, Dimension(:), Pointer               :: Has_Force
     
     ! NODE SETS PROPERTIES
     Integer, Dimension(:), Pointer               :: BC_Type_X
     Integer, Dimension(:), Pointer               :: BC_Type_Y
     Integer, Dimension(:), Pointer               :: BC_Type_Z

     ! GLOBAL PARAMETERS (STORED AS RESULTS)
     Real(Kind = Kr), Dimension(:), Pointer       :: Load

     ! ELEMENT PARAMETERS (STORED AS RESULTS)

     ! NODES PROPERTIES (STORES AS RESULTS)

     !!! STORED IN THE PARAM FILE
     ! GLOBAL DATAS
     ! The following parameter are NOT STORED in the Exodus db
     Logical                                      :: Do_Irrev
     
     Logical                                      :: Do_BackTrack
     Real(Kind = Kr)                              :: Tol_Ener
               
     Integer                                      :: Init_U
     
     Integer                                      :: Init_V
     Integer                                      :: nbCracks
     Real(Kind = Kr)                              :: MaxCrackLength     

     Integer                                      :: MaxIterRelax
     Real(Kind = Kr)                              :: TolRelax
     Real(Kind = Kr)                              :: TolKSP

     Real(Kind = Kr)                              :: Epsilon
     Real(Kind = Kr)                              :: KEpsilon
     Real(Kind = Kr)                              :: TolIrrev

     !!! STORED IN THE CST FILE
     ! ELEMENT BLOCKS PARAMETERS
     Real(Kind = Kr), Dimension(:), Pointer       :: Toughness
     Type(Tens4OS3D), Dimension(:), Pointer       :: Hookes_Law
     Real(Kind = Kr), Dimension(:), Pointer       :: Therm_Exp
  End Type Rupt_Params3D

Contains

  Subroutine Read_Rupt_EXO_Params_Iso(Geom, Params)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params)                            :: Params
    
    Character(len=MXSTLN), Dimension(:,:), Pointer:: Tmp_QA
    Integer                                       :: iErr
    Integer                                       :: iTs, iS, iBlk, iSet
    Integer                                       :: Num_TS
    Integer                                       :: exo_ver

    Real(Kind = Kr)                               :: fDum
    Character                                     :: cDum
    Integer                                       :: Tmp_EB_Prop
  

    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, exo_ver,&
         & ierr)

    Allocate(Params%Is_Brittle(Geom%Num_Elem_Blks))
    Allocate(Params%Is_Domain(Geom%Num_Elem_Blks))
    Allocate(Params%Has_Force(Geom%Num_Elem_Blks))

    Do iBlk = 1, Geom%Num_Elem_Blks
       ! Object Properties
       Call EXGP(Geom%exoid, EXEBLK, iBlk, Prop_Name_EB(1),                   &
            & Tmp_EB_Prop, iErr)
       Params%Is_Brittle(iBlk) = Tmp_EB_Prop
       Call EXGP(Geom%exoid, EXEBLK, iBlk, Prop_Name_EB(2),                   &
            & Tmp_EB_Prop, iErr)
       Params%Is_Domain(iBlk) = Tmp_EB_Prop
       Call EXGP(Geom%exoid, EXEBLK, iBlk, Prop_Name_EB(3),                   &
            & Tmp_EB_Prop, iErr)
       Params%Has_Force(iBlk) = Tmp_EB_Prop
    End Do    
    
    Allocate(Params%BC_Type_X(Geom%Num_Node_Sets))
    Allocate(Params%BC_Type_Y(Geom%Num_Node_Sets))
    Allocate(Params%BC_Type_Z(Geom%Num_Node_Sets))
    Do iSet = 1, Geom%Num_Node_Sets
       Call EXGP(Geom%exoid, EXNSET, iSet, Prop_Name_NS(1),                   &
            & Params%BC_Type_X(iSet), iErr)
       Call EXGP(Geom%exoid, EXNSET, iSet, Prop_Name_NS(2),                   &
            & Params%BC_Type_Y(iSet), iErr)
       Call EXGP(Geom%exoid, EXNSET, iSet, Prop_Name_NS(3),                   &
            & Params%BC_Type_Z(iSet), iErr)   
    End Do

    ! Time Step
    ! Time has to be monotonically increasing, so we store the number 
    ! of the time step instead, and store the displacement factor or the 
    ! temperature as a global result
    Call EXINQ (Geom%exoid, EXTIMS, iTS, fDum, cDum, iErr)    

    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0

    Allocate (Params%Load(iTS))
    Do iTS = 1, Size(Params%Load)
       Call Read_EXO_Result_Global(Geom, 4, iTS, Params%Load(iTS))
    End Do
  End Subroutine Read_Rupt_EXO_Params_Iso

  Subroutine Read_Rupt_EXO_Params2D(Geom, Params)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params2D)                          :: Params
    
    Character(len=MXSTLN), Dimension(:,:), Pointer:: Tmp_QA
    Integer                                       :: iErr
    Integer                                       :: iTs, iS, iBlk, iSet
    Integer                                       :: Num_TS
    Integer                                       :: exo_ver

    Real(Kind = Kr)                               :: fDum
    Character                                     :: cDum
    Integer                                       :: Tmp_EB_Prop
  

    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, exo_ver, ierr)

    Allocate(Params%Is_Brittle(Geom%Num_Elem_Blks))
    Allocate(Params%Is_Domain(Geom%Num_Elem_Blks))
    Allocate(Params%Has_Force(Geom%Num_Elem_Blks))

    Do iBlk = 1, Geom%Num_Elem_Blks
       ! Object Properties
       Call EXGP(Geom%exoid, EXEBLK, iBlk, Prop_Name_EB(1),                   &
            & Tmp_EB_Prop, iErr)
       Params%Is_Brittle(iBlk) = Tmp_EB_Prop
       Call EXGP(Geom%exoid, EXEBLK, iBlk, Prop_Name_EB(2),                   &
            & Tmp_EB_Prop, iErr)
       Params%Is_Domain(iBlk) = Tmp_EB_Prop
       Call EXGP(Geom%exoid, EXEBLK, iBlk, Prop_Name_EB(3),                   &
            & Tmp_EB_Prop, iErr)
       Params%Has_Force(iBlk) = Tmp_EB_Prop
    End Do    
    
    Allocate(Params%BC_Type_X(Geom%Num_Node_Sets))
    Allocate(Params%BC_Type_Y(Geom%Num_Node_Sets))
    Allocate(Params%BC_Type_Z(Geom%Num_Node_Sets))
    Do iSet = 1, Geom%Num_Node_Sets
       Call EXGP(Geom%exoid, EXNSET, iSet, Prop_Name_NS(1),                   &
            & Params%BC_Type_X(iSet), iErr)
       Call EXGP(Geom%exoid, EXNSET, iSet, Prop_Name_NS(2),                   &
            & Params%BC_Type_Y(iSet), iErr)
       Call EXGP(Geom%exoid, EXNSET, iSet, Prop_Name_NS(3),                   &
            & Params%BC_Type_Z(iSet), iErr)   
    End Do

    ! Time Step
    ! Time has to be monotonically increasing, so we store the number 
    ! of the time step instead, and store the displacement factor or the 
    ! temperature as a global result
    Call EXINQ (Geom%exoid, EXTIMS, iTS, fDum, cDum, iErr)    

    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0

    Allocate (Params%Load(iTS))
    Do iTS = 1, Size(Params%Load)
       Call Read_EXO_Result_Global(Geom, 4, iTS, Params%Load(iTS))
    End Do
  End Subroutine Read_Rupt_EXO_Params2D

  Subroutine Read_Rupt_EXO_Params3D(Geom, Params)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params3D)                          :: Params
    
    Character(len=MXSTLN), Dimension(:,:), Pointer:: Tmp_QA
    Integer                                       :: iErr
    Integer                                       :: iTs, iS, iBlk, iSet
    Integer                                       :: Num_TS
    Integer                                       :: exo_ver

    Real(Kind = Kr)                               :: fDum
    Character                                     :: cDum
    Integer                                       :: Tmp_EB_Prop
  

    Geom%exoid = EXOPEN(Geom%filename, EXREAD, exo_cpu_ws, exo_io_ws, exo_ver,&
         & ierr)

    Allocate(Params%Is_Brittle(Geom%Num_Elem_Blks))
    Allocate(Params%Is_Domain(Geom%Num_Elem_Blks))
    Allocate(Params%Has_Force(Geom%Num_Elem_Blks))

    Do iBlk = 1, Geom%Num_Elem_Blks
       ! Object Properties
       Call EXGP(Geom%exoid, EXEBLK, iBlk, Prop_Name_EB(1),                   &
            & Tmp_EB_Prop, iErr)
       Params%Is_Brittle(iBlk) = Tmp_EB_Prop
       Call EXGP(Geom%exoid, EXEBLK, iBlk, Prop_Name_EB(2),                   &
            & Tmp_EB_Prop, iErr)
       Params%Is_Domain(iBlk) = Tmp_EB_Prop
       Call EXGP(Geom%exoid, EXEBLK, iBlk, Prop_Name_EB(3),                   &
            & Tmp_EB_Prop, iErr)
       Params%Has_Force(iBlk) = Tmp_EB_Prop
    End Do    
    
    Allocate(Params%BC_Type_X(Geom%Num_Node_Sets))
    Allocate(Params%BC_Type_Y(Geom%Num_Node_Sets))
    Allocate(Params%BC_Type_Z(Geom%Num_Node_Sets))
    Do iSet = 1, Geom%Num_Node_Sets
       Call EXGP(Geom%exoid, EXNSET, iSet, Prop_Name_NS(1),                   &
            & Params%BC_Type_X(iSet), iErr)
       Call EXGP(Geom%exoid, EXNSET, iSet, Prop_Name_NS(2),                   &
            & Params%BC_Type_Y(iSet), iErr)
       Call EXGP(Geom%exoid, EXNSET, iSet, Prop_Name_NS(3),                   &
            & Params%BC_Type_Z(iSet), iErr)   
    End Do

    ! Time Step
    ! Time has to be monotonically increasing, so we store the number 
    ! of the time step instead, and store the displacement factor or the 
    ! temperature as a global result
    Call EXINQ (Geom%exoid, EXTIMS, iTS, fDum, cDum, iErr)    

    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0

    Allocate (Params%Load(iTS))
    Do iTS = 1, Size(Params%Load)
       Call Read_EXO_Result_Global(Geom, 4, iTS, Params%Load(iTS))
    End Do
  End Subroutine Read_Rupt_EXO_Params3D

  Subroutine Write_Rupt_EXO_Params_Iso(Geom, Params)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params)                            :: Params
    
    Character(len=MXSTLN), Dimension(:,:), Pointer:: Tmp_QA
    Integer                                       :: iErr
    Integer                                       :: iTs, iS
    Integer                                       :: Num_TS
    Integer                                       :: exo_ver

    Integer, Dimension(:), Pointer                :: Tmp_EB_Prop
    
    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, exo_ver,&
         & ierr)
    ! Object Properties
    Call EXPPN(Geom%exoid, EXEBLK, Num_Prop_EB, Prop_Name_EB, iErr)
    Call EXPPN(Geom%exoid, EXNSET, Num_Prop_NS, Prop_Name_NS, iErr)


    Allocate (Tmp_EB_Prop(Geom%Num_Elem_Blks))
    Tmp_EB_Prop = Params%Is_Brittle
    Call EXPPA(Geom%exoid, EXEBLK, Prop_Name_EB(1), Tmp_EB_Prop, iErr)

    Tmp_EB_Prop = Params%Is_Domain
    Call EXPPA(Geom%exoid, EXEBLK, Prop_Name_EB(2), Tmp_EB_Prop, iErr)

    Tmp_EB_Prop = Params%Has_Force
    Call EXPPA(Geom%exoid, EXEBLK, Prop_Name_EB(3), Tmp_EB_Prop, iErr)
    DeAllocate (Tmp_EB_Prop)

    Call EXPPA(Geom%exoid, EXNSET, Prop_Name_NS(1), Params%BC_Type_X, iErr)
    Call EXPPA(Geom%exoid, EXNSET, Prop_Name_NS(2), Params%BC_Type_Y, iErr)
    Call EXPPA(Geom%exoid, EXNSET, Prop_Name_NS(3), Params%BC_Type_Z, iErr)

    !Result and initial values names
    Call EXPVP (Geom%exoid, 'g', Num_Res_G, iErr)
    Call EXPVAN (Geom%exoid, 'g', Num_Res_G, Res_Name_G, iErr)
    
    Call EXPVP (Geom%exoid, 'e', Num_Res_E, iErr)
    Call EXPVAN (Geom%exoid, 'e', Num_Res_E, Res_Name_E, iErr)
    
    Call EXPVP (Geom%exoid, 'n', Num_Res_N, iErr)
    Call EXPVAN (Geom%exoid, 'n', Num_Res_N, Res_Name_N, iErr)
    
    ! Time Step
    ! Time has to be monotonically increasing, so we store the number 
    ! of the time step instead, and store the displacement factor or the 
    ! temperature as a global result
    Do iTS = 1, Size(Params%Load)
       Call EXPTIM(Geom%exoid, iTS, real(iTS, 8), iErr)
       Call EXPGV (Geom%exoid, iTS, 4, (/ 0.0_Kr, 0.0_Kr, 0.0_Kr,          &
            & Params%Load(iTS)/), iErr)
    End Do

    Geom%Num_QA = Geom%Num_QA+1
    Allocate (Tmp_QA(4, Geom%Num_QA))
    Tmp_QA(:,1:Geom%Num_QA-1) = Geom%QA_rec
    DeAllocate (Geom%QA_Rec)
    Allocate (Geom%QA_Rec(4, Geom%Num_QA))
    Geom%QA_rec = Tmp_QA
    DeAllocate (Tmp_QA)
   
    Geom%QA_Rec(1,Geom%Num_QA) = 'm_Rupt-Struct'
    Geom%QA_Rec(2,Geom%Num_QA) = ''
    Call Date_And_Time(date = Geom%QA_Rec(3,Geom%Num_QA))
    Call Date_And_Time(time = Geom%QA_Rec(4,Geom%Num_QA))
    Call EXPQA(Geom%exoid, Geom%num_QA, Geom%QA_Rec, iErr)

    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_Rupt_EXO_Params_Iso
  
  Subroutine Write_Rupt_EXO_Params2D(Geom, Params)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params2D)                          :: Params
    
    Character(len=MXSTLN), Dimension(:,:), Pointer:: Tmp_QA
    Integer                                       :: iErr
    Integer                                       :: iTs, iS
    Integer                                       :: Num_TS
    Integer                                       :: exo_ver

    Integer, Dimension(:), Pointer                :: Tmp_EB_Prop
    
    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, exo_ver,&
         & ierr)
    ! Object Properties
    Call EXPPN(Geom%exoid, EXEBLK, Num_Prop_EB, Prop_Name_EB, iErr)
    Call EXPPN(Geom%exoid, EXNSET, Num_Prop_NS, Prop_Name_NS, iErr)


    Allocate (Tmp_EB_Prop(Geom%Num_Elem_Blks))
    Tmp_EB_Prop = Params%Is_Brittle
    Call EXPPA(Geom%exoid, EXEBLK, Prop_Name_EB(1), Tmp_EB_Prop, iErr)

    Tmp_EB_Prop = Params%Is_Domain
    Call EXPPA(Geom%exoid, EXEBLK, Prop_Name_EB(2), Tmp_EB_Prop, iErr)

    Tmp_EB_Prop = Params%Has_Force
    Call EXPPA(Geom%exoid, EXEBLK, Prop_Name_EB(3), Tmp_EB_Prop, iErr)
    DeAllocate (Tmp_EB_Prop)

    Call EXPPA(Geom%exoid, EXNSET, Prop_Name_NS(1), Params%BC_Type_X, iErr)
    Call EXPPA(Geom%exoid, EXNSET, Prop_Name_NS(2), Params%BC_Type_Y, iErr)
    Call EXPPA(Geom%exoid, EXNSET, Prop_Name_NS(3), Params%BC_Type_Z, iErr)

    !Result and initial values names
    Call EXPVP (Geom%exoid, 'g', Num_Res_G, iErr)
    Call EXPVAN (Geom%exoid, 'g', Num_Res_G, Res_Name_G, iErr)
    
    Call EXPVP (Geom%exoid, 'e', Num_Res_E, iErr)
    Call EXPVAN (Geom%exoid, 'e', Num_Res_E, Res_Name_E, iErr)
    
    Call EXPVP (Geom%exoid, 'n', Num_Res_N, iErr)
    Call EXPVAN (Geom%exoid, 'n', Num_Res_N, Res_Name_N, iErr)
    
    ! Time Step
    ! Time has to be monotonically increasing, so we store the number 
    ! of the time step instead, and store the displacement factor or the 
    ! temperature as a global result
    Do iTS = 1, Size(Params%Load)
       Call EXPTIM(Geom%exoid, iTS, real(iTS, 8), iErr)
       Call EXPGV (Geom%exoid, iTS, 4, (/ 0.0_Kr, 0.0_Kr, 0.0_Kr,          &
            & Params%Load(iTS)/), iErr)
    End Do

    Geom%Num_QA = Geom%Num_QA+1
    Allocate (Tmp_QA(4, Geom%Num_QA))
    Tmp_QA(:,1:Geom%Num_QA-1) = Geom%QA_rec
    DeAllocate (Geom%QA_Rec)
    Allocate (Geom%QA_Rec(4, Geom%Num_QA))
    Geom%QA_rec = Tmp_QA
    DeAllocate (Tmp_QA)
   
    Geom%QA_Rec(1,Geom%Num_QA) = 'm_Rupt-Struct'
    Geom%QA_Rec(2,Geom%Num_QA) = ''
    Call Date_And_Time(date = Geom%QA_Rec(3,Geom%Num_QA))
    Call Date_And_Time(time = Geom%QA_Rec(4,Geom%Num_QA))
    Call EXPQA(Geom%exoid, Geom%num_QA, Geom%QA_Rec, iErr)

    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_Rupt_EXO_Params2D
  
  Subroutine Write_Rupt_EXO_Params3D(Geom, Params)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params3D)                          :: Params
    
    Character(len=MXSTLN), Dimension(:,:), Pointer:: Tmp_QA
    Integer                                       :: iErr
    Integer                                       :: iTs, iS
    Integer                                       :: Num_TS
    Integer                                       :: exo_ver

    Integer, Dimension(:), Pointer                :: Tmp_EB_Prop

    Geom%exoid = EXOPEN(Geom%filename, EXWRIT, exo_cpu_ws, exo_io_ws, exo_ver,&
         & ierr)
    ! Object Properties
    Call EXPPN(Geom%exoid, EXEBLK, Num_Prop_EB, Prop_Name_EB, iErr)
    Call EXPPN(Geom%exoid, EXNSET, Num_Prop_NS, Prop_Name_NS, iErr)

    Allocate (Tmp_EB_Prop(Geom%Num_Elem_Blks))
    Tmp_EB_Prop = Params%Is_Brittle
    Call EXPPA(Geom%exoid, EXEBLK, Prop_Name_EB(1), Tmp_EB_Prop, iErr)
    Tmp_EB_Prop = Params%Is_Domain
    Call EXPPA(Geom%exoid, EXEBLK, Prop_Name_EB(2), Tmp_EB_Prop, iErr)
    Tmp_EB_Prop = Params%Has_Force
    Call EXPPA(Geom%exoid, EXEBLK, Prop_Name_EB(3), Tmp_EB_Prop, iErr)
    DeAllocate (Tmp_EB_Prop)
    Call EXPPA(Geom%exoid, EXNSET, Prop_Name_NS(1), Params%BC_Type_X, iErr)
    Call EXPPA(Geom%exoid, EXNSET, Prop_Name_NS(2), Params%BC_Type_Y, iErr)
    Call EXPPA(Geom%exoid, EXNSET, Prop_Name_NS(3), Params%BC_Type_Z, iErr)

    !Result and initial values names
    Call EXPVP (Geom%exoid, 'g', Num_Res_G, iErr)
    Call EXPVAN (Geom%exoid, 'g', Num_Res_G, Res_Name_G, iErr)
    
    Call EXPVP (Geom%exoid, 'e', Num_Res_E, iErr)
    Call EXPVAN (Geom%exoid, 'e', Num_Res_E, Res_Name_E, iErr)
    
    Call EXPVP (Geom%exoid, 'n', Num_Res_N, iErr)
    Call EXPVAN (Geom%exoid, 'n', Num_Res_N, Res_Name_N, iErr)
    
    ! Time Step
    ! Time has to be monotonically increasing, so we store the number 
    ! of the time step instead, and store the displacement factor or the 
    ! temperature as a global result
    Do iTS = 1, Size(Params%Load)
       Call EXPTIM(Geom%exoid, iTS, real(iTS, 8), iErr)
       Call EXPGV (Geom%exoid, iTS, 4, (/ 0.0_Kr, 0.0_Kr, 0.0_Kr,          &
            & Params%Load(iTS)/), iErr)
    End Do

!    Geom%Num_QA = Geom%Num_QA+1
!    Allocate (Tmp_QA(4, Geom%Num_QA))
!    Tmp_QA(:,1:Geom%Num_QA-1) = Geom%QA_rec
!    DeAllocate (Geom%QA_Rec)
!    Allocate (Geom%QA_Rec(4, Geom%Num_QA))
!    Geom%QA_rec = Tmp_QA
!    DeAllocate (Tmp_QA)
   
!    Geom%QA_Rec(1,Geom%Num_QA) = 'm_Rupt-Struct'
!    Geom%QA_Rec(2,Geom%Num_QA) = ''
!    Call Date_And_Time(date = Geom%QA_Rec(3,Geom%Num_QA))
!    Call Date_And_Time(time = Geom%QA_Rec(4,Geom%Num_QA))
!    Call EXPQA(Geom%exoid, Geom%num_QA, Geom%QA_Rec, iErr)

    Call EXCLOS(Geom%exoid, iErr)
    Geom%exoid = 0
  End Subroutine Write_Rupt_EXO_Params3D
  
  Subroutine Write_Rupt_DATA_Iso(Geom, Params)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params)                            :: Params
    
    Integer                                       :: iBlk, Blk_ID

    Open(File = Params%PARAM_Str, Unit = F_OUT, Status = 'Unknown')
    Rewind(F_OUT)

    Write(F_OUT, 101) Params%Do_Irrev,     'Do_Irrev [T/F]'
    If (Params%Do_BackTrack) Then
       Write(F_OUT, 103) Params%Do_BackTrack, Params%Tol_Ener ,               &
            & 'Do_BackTrack [T/F] &&  Tol_Ener '
    Else
       Write(F_OUT, 101) Params%Do_BackTrack,      'Do_BackTrack[T/F]'
    End If    
    Write(F_OUT, 100) Params%Init_U, 'Init_U   [0=PREV 1=0]'
    If (Params%Init_V==2) Then
       Write(F_OUT, 105) Params%Init_V, Params%nbCracks,                      &
            & Params%MaxCrackLength, 'Init_V   [0=PREV 1=1, 2=RND]  '//       &
            & 'nbCracks MaxCrackLength'    
    Else    
       Write(F_OUT, 100) Params%Init_V,       'Init_V   [0=PREV 1=1, 2=RND]'
    End If  
    Write(F_OUT, 100) Params%MaxIterRelax, 'MaxIterRelax'
    Write(F_OUT, 102) Params%TolRelax,     'TolRelax'
    Write(F_OUT, 102) Params%TolKSP,       'TolKSP'
    Write(F_OUT, 102) Params%Epsilon,      'Epsilon'
    Write(F_OUT, 102) Params%KEpsilon,     'KEpsilon'
    Write(F_OUT, 102) Params%TolIrrev,     'TolIrrev'
    Close(F_OUT)

    Open(File = Params%CST_Str, Unit = F_OUT, Status = 'Unknown')
    Rewind(F_OUT)
    Write(F_OUT, 110) Geom%Num_Elem_Blks
    Do iBlk = 1, Geom%Num_Elem_Blks
       Blk_ID = Geom%Elem_Blk(iBlk)%ID
       Write(F_OUT,120) Blk_ID, Params%Toughness(iBlk),                       &
            & Params%Young_Mod(iBlk), Params%Poisson_Ratio(iBlk),             &
            & Params%Therm_Exp(iBlk)
    End Do
    Close(F_OUT)

100 Format(I4,T30,'# ', A)
101 Format(L1,T30,'# ', A)
102 Format(ES12.5,T30,'# ', A)
103 Format(L1,ES12.5,T30,'# ', A)
105 Format(I4,I4,ES12.5,T30, '# ', A)
110 Format(I4,T30,'# BLK_ID, Toughness, E, nu, Alpha')
120 Format(I4, 4(ES12.5,' '))
  End Subroutine Write_Rupt_DATA_Iso
 
  Subroutine Write_Rupt_DATA2D(Geom, Params)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params2D)                          :: Params
    
    Integer                                       :: iBlk, Blk_ID

    Open(File = Params%PARAM_Str, Unit = F_OUT, Status = 'Unknown')
    Rewind(F_OUT)

    Write(F_OUT, 101) Params%Do_Irrev,     'Do_Irrev [T/F]'
    If (Params%Do_BackTrack) Then
       Write(F_OUT, 103) Params%Do_BackTrack, Params%Tol_Ener , 'Do_BackTrack [T/F] &&  Tol_Ener '
    Else
       Write(F_OUT, 101) Params%Do_BackTrack,      'Do_BackTrack[T/F]'
    End If    
    Write(F_OUT, 100) Params%Init_U, 'Init_U   [0=PREV 1=0]'
    If (Params%Init_V==2) Then
       Write(F_OUT, 105) Params%Init_V, Params%nbCracks, Params%MaxCrackLength, 'Init_V   [0=PREV 1=1, 2=RND]  '//                &
            & 'nbCracks MaxCrackLength'    
    Else    
       Write(F_OUT, 100) Params%Init_V,       'Init_V   [0=PREV 1=1, 2=RND]'
    End If  
    Write(F_OUT, 100) Params%MaxIterRelax, 'MaxIterRelax'
    Write(F_OUT, 102) Params%TolRelax,     'TolRelax'
    Write(F_OUT, 102) Params%TolKSP,       'TolKSP'
    Write(F_OUT, 102) Params%Epsilon,      'Epsilon'
    Write(F_OUT, 102) Params%KEpsilon,     'KEpsilon'
    Write(F_OUT, 102) Params%TolIrrev,     'TolIrrev'
    Close(F_OUT)

    Open(File = Params%CST_Str, Unit = F_OUT, Status = 'Unknown')
    Rewind(F_OUT)
    Write(F_OUT, 110) Geom%Num_Elem_Blks
    Do iBlk = 1, Geom%Num_Elem_Blks
       Blk_ID = Geom%Elem_Blk(iBlk)%ID
       Write(F_OUT,120) Blk_ID, Params%Toughness(iBlk), Params%Hookes_Law(iBlk), Params%Therm_Exp(iBlk)
    End Do
    Close(F_OUT)

100 Format(I4,T30,'# ', A)
101 Format(L1,T30,'# ', A)
102 Format(ES12.5,T30,'# ', A)
103 Format(L1,ES12.5,T30,'# ', A)
105 Format(I4,I4,ES12.5,T30, '# ', A)
110 Format(I6,'     Toughness    A1111        A1112        A1122        A1212        A1222        A2222        Alpha')
120 Format(I6, 8(ES12.5,' '))
  End Subroutine Write_Rupt_DATA2D

  Subroutine Write_Rupt_DATA3D(Geom, Params)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params3D)                          :: Params
    
    Integer                                       :: iBlk, Blk_ID

    Open(File = Params%PARAM_Str, Unit = F_OUT, Status = 'Unknown')
    Rewind(F_OUT)

    Write(F_OUT, 101) Params%Do_Irrev,     'Do_Irrev [T/F]'
    If (Params%Do_BackTrack) Then
       Write(F_OUT, 103) Params%Do_BackTrack, Params%Tol_Ener , 'Do_BackTrack [T/F] &&  Tol_Ener '
    Else
       Write(F_OUT, 101) Params%Do_BackTrack,      'Do_BackTrack[T/F]'
    End If    
    Write(F_OUT, 100) Params%Init_U, 'Init_U   [0=PREV 1=0]'
    If (Params%Init_V==2) Then
       Write(F_OUT, 105) Params%Init_V, Params%nbCracks, Params%MaxCrackLength, 'Init_V   [0=PREV 1=1, 2=RND]  '//                &
            & 'nbCracks MaxCrackLength'    
    Else    
       Write(F_OUT, 100) Params%Init_V,       'Init_V   [0=PREV 1=1, 2=RND]'
    End If  
    Write(F_OUT, 100) Params%MaxIterRelax, 'MaxIterRelax'
    Write(F_OUT, 102) Params%TolRelax,     'TolRelax'
    Write(F_OUT, 102) Params%TolKSP,       'TolKSP'
    Write(F_OUT, 102) Params%Epsilon,      'Epsilon'
    Write(F_OUT, 102) Params%KEpsilon,     'KEpsilon'
    Write(F_OUT, 102) Params%TolIrrev,     'TolIrrev'
    Close(F_OUT)

    Open(File = Params%CST_Str, Unit = F_OUT, Status = 'Unknown')
    Rewind(F_OUT)
    Write(F_OUT, 110) Geom%Num_Elem_Blks
    Do iBlk = 1, Geom%Num_Elem_Blks
       Blk_ID = Geom%Elem_Blk(iBlk)%ID
       Write(F_OUT,120) Blk_ID, Params%Toughness(iBlk), Params%Hookes_Law(iBlk), Params%Therm_Exp(iBlk)
    End Do
    Close(F_OUT)

100 Format(I4,T30,'# ', A)
101 Format(L1,T30,'# ', A)
102 Format(ES12.5,T30,'# ', A)
103 Format(L1,ES12.5,T30,'# ', A)
105 Format(I4,I4,ES12.5,T30, '# ', A)
110 Format(I6,' Toughness    A_1111       A_1112       A_1113       A_1122       A_1123       A_1133       A_1212       A_1213       A_1222       A_1223       A_12133      A_1313       A_1322       A_1323       A_1333       A_2222       A_2223       A_2233       A_2323       A_2333       A_3333       Alpha')
120 Format(I6, 23(ES12.5,' '))
  End Subroutine Write_Rupt_DATA3D

  Subroutine Read_Rupt_DATA_Iso(Geom, Params)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params)                            :: Params
    
    Integer                                       :: iBlk, Blk_ID, Num_Blk
    Real(Kind = Kr)                               :: Tmp_Toughness
    Real(Kind = Kr)                               :: Tmp_Young_Mod
    Real(Kind = Kr)                               :: Tmp_Poisson_Ratio
    Real(Kind = Kr)                               :: Tmp_Therm_Exp 

    Open(File = Params%PARAM_Str, Unit = F_IN, Status = 'Old')
    Rewind(F_IN)

    Read(F_IN, 101) Params%Do_Irrev
    Read(F_IN, 101, advance = 'no' ) Params%Do_BackTrack
    If (Params%Do_BackTrack) Then
       Read(F_IN, 102) Params%Tol_Ener
    Else
       Read(F_In,*)
    End If   
    Read(F_IN, 100) Params%Init_U
    Read(F_IN, 100, advance = 'no' ) Params%Init_V
    
    If ((Params%Init_V == Init_V_SPH) .OR. (Params%Init_V == Init_V_RND)) Then
      Read(F_IN, *) Params%nbCracks,Params%MaxCrackLength
    Else
       Read(F_In,*)
    End If
    Read(F_IN, 100) Params%MaxIterRelax     
    Read(F_IN, 102) Params%TolRelax
    Read(F_IN, 102) Params%TolKSP
    Read(F_IN, 102) Params%Epsilon
    Read(F_IN, 102) Params%KEpsilon
    Read(F_IN, 102) Params%TolIrrev
    Close(F_IN)
    
    Open(File = Params%CST_Str, Unit = F_IN, Status = 'Old')
    Rewind(F_IN)
    Read(F_IN, 110) Num_Blk
    If (Num_Blk /= Geom%Num_Elem_Blks) Then
       Write(*,*) 'ERROR number of elem blocks in CST file doesn''t',         &
            & ' match that of the exodus DB', Num_Blk, Geom%Num_Elem_Blks
       STOP
    End If
    
    Allocate (Params%Toughness(Num_Blk))
    Allocate (Params%Young_Mod(Num_Blk))
    Allocate (Params%Poisson_Ratio(Num_Blk))
    Allocate (Params%Therm_Exp(Num_Blk))

    Do iBlk = 1, Num_Blk
       Read(F_IN,*) Blk_ID, Tmp_Toughness, Tmp_Young_Mod, Tmp_Poisson_Ratio,  &
       	&           Tmp_Therm_Exp
       If (Blk_ID /= iBlk) Then
          Write(*,*) 'ERROR: Element Block must be numbered sequentially in', &
               & ' the CST file'
          STOP
       End If
       Params%Toughness(iBlk)     = Tmp_Toughness
       Params%Young_Mod(iBlk)     = Tmp_Young_Mod
       Params%Poisson_Ratio(iBlk) = Tmp_Poisson_ratio
       Params%Therm_Exp(iBlk)     = Tmp_Therm_Exp
    End Do
    Close(F_IN)

100 Format(I4)
101 Format(L1)
102 Format(ES12.5)
103 Format(I4,' ',ES12.5)
110 Format(I4)
120 Format(I4, 6(ES12.5,' '))
  End Subroutine Read_Rupt_DATA_Iso

  Subroutine Read_Rupt_DATA2D(Geom, Params)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params2D)                          :: Params
    
    Integer                                       :: iBlk, Blk_ID, Num_Blk
    Real(Kind = Kr)                               :: Tmp_Toughness
    Type(Tens4OS2D)                               :: Tmp_Hookes_Law
    Real(Kind = Kr)                               :: Tmp_Therm_Exp 

    Open(File = Params%PARAM_Str, Unit = F_IN, Status = 'Old')
    Rewind(F_IN)

    Read(F_IN, 101) Params%Do_Irrev
    Read(F_IN, 101, advance = 'no' ) Params%Do_BackTrack
    If (Params%Do_BackTrack) Then
       Read(F_IN, 102) Params%Tol_Ener
    Else
       Read(F_In,*)
    End If   
    Read(F_IN, 100) Params%Init_U
    Read(F_IN, 100, advance = 'no' ) Params%Init_V
    
    If ((Params%Init_V == Init_V_SPH) .OR. (Params%Init_V == Init_V_RND)) Then
      Read(F_IN, *) Params%nbCracks,Params%MaxCrackLength
    Else
       Read(F_In,*)
    End If
    Read(F_IN, 100) Params%MaxIterRelax     
    Read(F_IN, 102) Params%TolRelax
    Read(F_IN, 102) Params%TolKSP
    Read(F_IN, 102) Params%Epsilon
    Read(F_IN, 102) Params%KEpsilon
    Read(F_IN, 102) Params%TolIrrev
    Close(F_IN)
    
    Open(File = Params%CST_Str, Unit = F_IN, Status = 'Old')
    Rewind(F_IN)
    Read(F_IN, 110) Num_Blk
    If (Num_Blk /= Geom%Num_Elem_Blks) Then
       Write(*,*) 'ERROR number of elem blocks in CST file doesn''t',         &
            & ' match that of the exodus DB', Num_Blk, Geom%Num_Elem_Blks
       STOP
    End If
    
    Allocate (Params%Toughness(Num_Blk))
    Allocate (Params%Hookes_Law(Num_Blk))
    Allocate (Params%Therm_Exp(Num_Blk))

    Do iBlk = 1, Num_Blk
       Read(F_IN,*) Blk_ID, Tmp_Toughness, Tmp_Hookes_Law, Tmp_Therm_Exp
       If (Blk_ID /= iBlk) Then
          Write(*,*) 'ERROR: Element Block must be numbered sequentially in', &
               & ' the CST file'
          STOP
       End If
       Params%Toughness(iBlk)     = Tmp_Toughness
       Params%Hookes_Law(iBlk)    = Tmp_Hookes_Law
       Params%Therm_Exp(iBlk)     = Tmp_Therm_Exp
    End Do
    Close(F_IN)

100 Format(I4)
101 Format(L1)
102 Format(ES12.5)
103 Format(I6,' ',ES12.5)
!110 Format(I6)
  End Subroutine Read_Rupt_DATA2D
  
  Subroutine Read_Rupt_DATA3D(Geom, Params)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params3D)                          :: Params
    
    Integer                                       :: iBlk, Blk_ID, Num_Blk
    Real(Kind = Kr)                               :: Tmp_Toughness
    Type(Tens4OS3D)                               :: Tmp_Hookes_Law
    Real(Kind = Kr)                               :: Tmp_Therm_Exp 

    Open(File = Params%PARAM_Str, Unit = F_IN, Status = 'Old')
    Rewind(F_IN)

    Read(F_IN, 101) Params%Do_Irrev
    Read(F_IN, 101, advance = 'no' ) Params%Do_BackTrack
    If (Params%Do_BackTrack) Then
       Read(F_IN, 102) Params%Tol_Ener
    Else
       Read(F_In,*)
    End If   
    Read(F_IN, 100) Params%Init_U
    Read(F_IN, 100, advance = 'no' ) Params%Init_V
    
    If ((Params%Init_V == Init_V_SPH) .OR. (Params%Init_V == Init_V_RND)) Then
      Read(F_IN, *) Params%nbCracks,Params%MaxCrackLength
    Else
       Read(F_In,*)
    End If
    Read(F_IN, 100) Params%MaxIterRelax     
    Read(F_IN, 102) Params%TolRelax
    Read(F_IN, 102) Params%TolKSP
    Read(F_IN, 102) Params%Epsilon
    Read(F_IN, 102) Params%KEpsilon
    Read(F_IN, 102) Params%TolIrrev
    Close(F_IN)
    
    Open(File = Params%CST_Str, Unit = F_IN, Status = 'Old')
    Rewind(F_IN)
    Read(F_IN, 110) Num_Blk
    If (Num_Blk /= Geom%Num_Elem_Blks) Then
       Write(*,*) 'ERROR number of elem blocks in CST file doesn''t',         &
            & ' match that of the exodus DB', Num_Blk, Geom%Num_Elem_Blks
       STOP
    End If
    
    Allocate (Params%Toughness(Num_Blk))
    Allocate (Params%Hookes_Law(Num_Blk))
    Allocate (Params%Therm_Exp(Num_Blk))

    Do iBlk = 1, Num_Blk
       Read(F_IN,*) Blk_ID, Tmp_Toughness, Tmp_Hookes_Law, Tmp_Therm_Exp
       If (Blk_ID /= iBlk) Then
          Write(*,*) 'ERROR: Element Block must be numbered sequentially in', &
               & ' the CST file'
          STOP
       End If
       Params%Toughness(iBlk)     = Tmp_Toughness
       Params%Hookes_Law(iBlk)    = Tmp_Hookes_Law
       Params%Therm_Exp(iBlk)     = Tmp_Therm_Exp
    End Do
    Close(F_IN)

100 Format(I4)
101 Format(L1)
102 Format(ES12.5)
103 Format(I6,' ',ES12.5)
!110 Format(I6)
  End Subroutine Read_Rupt_DATA3D
  
   Subroutine GenHL_Iso2D_LambdaMu(lambda, mu, A) 
      Real(Kind = Kr), Intent(IN)         :: lambda, mu
      Type(Tens4OS2D), Intent(OUT)        :: A
      A = 0.0_Kr
      
      A%XXXX = lambda + 2.0_Kr * mu
      A%XXYY = lambda
      A%XYXY = mu
      A%YYYY = lambda + 2.0_Kr * mu
   End Subroutine GenHL_Iso2D_LambdaMu         

   Subroutine GenHL_Iso2D_EnuPlaneStress(E, nu, A) 
      Real(Kind = Kr), Intent(IN)         :: E, nu
      Type(Tens4OS2D), Intent(OUT)        :: A
      
      Real(Kind = Kr)                     :: Lambda, mu
      
      lambda = E * nu / (1.0_Kr - nu**2) 
      mu     = E / (1.0_Kr + nu) * .5_Kr
      A = 0.0_Kr
      A%XXXX = lambda + 2.0_Kr * mu
      A%XXYY = lambda
      A%XYXY = mu
      A%YYYY = lambda + 2.0_Kr * mu
   End Subroutine GenHL_Iso2D_EnuPlaneStress         
   
   Subroutine GenHL_Iso2D_EnuPlaneStrain(E, nu, A) 
      Real(Kind = Kr), Intent(IN)         :: E, nu
      Type(Tens4OS2D), Intent(OUT)        :: A
      
      Real(Kind = Kr)                     :: Lambda, mu
      
      lambda = E * nu / (1.0_Kr + nu) / (1.0_Kr - 2.0_Kr * nu)
      mu     = E / (1.0_Kr + nu) * .5_Kr      
      A = 0.0_Kr
      A%XXXX = lambda + 2.0_Kr * mu
      A%XXYY = lambda
      A%XYXY = mu
      A%YYYY = lambda + 2.0_Kr * mu
   End Subroutine GenHL_Iso2D_EnuPlaneStrain         

   Subroutine GenHL_Iso3D_LambdaMu(lambda, mu, A)
      Real(Kind = Kr), Intent(IN)         :: Lambda, Mu
      Type(Tens4OS3D), Intent(OUT)        :: A
   
      A = 0.0_Kr
      A%XXXX = lambda + mu * 2.0_Kr
      A%XXYY = lambda
      A%XXZZ = lambda
      
      A%XYXY = mu
      
      A%XZXZ = mu
      
      A%YYYY = lambda + mu * 2.0_Kr
      A%YYZZ = lambda
      
      A%YZYZ = mu
      
      A%ZZZZ = lambda + mu * 2.0_Kr
   End Subroutine GenHL_Iso3D_LambdaMu

   Subroutine GenHL_Iso3D_Enu(E, nu, A)
      Real(Kind = Kr), Intent(IN)         :: E, nu
      Type(Tens4OS3D), Intent(OUT)        :: A
      
      Real(Kind = Kr)                     :: Lambda, mu
   
      lambda = E * nu / (1.0_Kr + nu) / (1 - 2.0_Kr * nu)
      mu     = E / (1.0_Kr + nu) * .5_Kr      
   
      A = 0.0_Kr
      A%XXXX = lambda + mu * 2.0_Kr
      A%XXYY = lambda
      A%XXZZ = lambda
      
      A%XYXY = mu
      
      A%XZXZ = mu
      
      A%YYYY = lambda + mu * 2.0_Kr
      A%YYZZ = lambda
      
      A%YZYZ = mu
      
      A%ZZZZ = lambda + mu * 2.0_Kr
   End Subroutine GenHL_Iso3D_Enu

   Subroutine GenHL_Ortho2D_LambdaMu(lambda, mu1, mu2, theta, A)
      Real(Kind = Kr), Intent(IN)         :: Lambda, mu1, mu2, theta
      Type(Tens4OS2D), Intent(OUT)        :: A
      
      A = 0.0_Kr
      A%XXXX = lambda + mu1 * (1.0_Kr + (cos(theta))**2) +  mu2 * (sin(theta))**2
      A%XXXY = (mu1-mu2) * cos(theta) * sin(theta)
      A%XXYY = lambda + (mu1-mu2) * (sin(theta))**2

      A%XYXY = mu1 * (sin(theta))**2 + mu2 * (cos(theta))**2
      A%XYYY = -A%XXXY

      A%YYYY =  A%XXXX
   End Subroutine GenHL_Ortho2D_LambdaMu
   

End Module m_Rupt_Struct
