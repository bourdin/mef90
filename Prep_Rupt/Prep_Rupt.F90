Program Prep_Rupt
  Use m_MEF90
  Use m_Rupt_Struct
  Implicit NONE

  Type(EXO_Geom_Info)                             :: Geom

  Integer                                         :: iS
  Type(Node1D), Dimension(:), Pointer             :: Node1D_db
  Type(Node2D), Dimension(:), Pointer             :: Node2D_db
  Type(Node3D), Dimension(:), Pointer             :: Node3D_db

  Type (Element1D), Dimension(:), Pointer         :: Elem1D_db
  Type (Element2D_Scal), Dimension(:), Pointer    :: Elem2D_db
  Type (Element3D), Dimension(:), Pointer         :: Elem3D_db

  Character                                       :: JobType
  
  Type(Vect3D), Dimension(:), Pointer             :: BC_U_NS


  Type (Rupt_Params)                              :: Params
  Character(len = MXLNLN)                         :: EXO_Str
!  Type(Vect3D), Dimension(:), Pointer             :: BC_U
     

  Integer                                         :: Num_dof_per_nodes 
  Integer                                         :: iErr
!  Character(len=MXLNLN), Dimension(3)             :: Coord_Names


  Integer                                         :: iB, iE, iV, iG
  Integer                                         :: iSL1, iSL2, iSG1, iSG2
  

  Real(Kind = Kr), Dimension(:), Pointer          :: Res
  Real(Kind = Kr)                                 :: Vol
  
  Integer, Parameter                              :: PB_Gen         = 0
  Integer, Parameter                              :: PB_MixedMode   = 1
  Integer, Parameter                              :: PB_Antiplane   = 2
  Integer, Parameter                              :: PB_CylTwist_3D = 3
  Integer, Parameter                              :: PB_Dipping     = 4
  Integer, Parameter                              :: PB_Needleman   = 5
  Integer                                         :: PB_Type



  Write(*,100, advance = 'no') 'ExodusII file name: '
  Read(*,100) EXO_Str

  Geom%filename = Trim(EXO_Str)//'.gen'
  Call Read_EXO_Geom_Info(Geom) 
  Print*, 'OK Read_EXO_Geom_Info'
  Print*, 'Dimension: ', Geom%Num_Dim

  Print*, 'Node sets ', Geom%Num_Node_Sets
  Do iB = 1, Geom%Num_Node_Sets
     Print*, 'Node Set #', Geom%Node_Set(iB)%ID
     Print*, '   Size ', Geom%Node_Set(iB)%Num_Nodes
     Print*, '   Nodes', Geom%Node_Set(iB)%Node_ID
  End Do
  
  Do
     Write(*,100, advance = 'no') '[P]repare [C]heck, or [E]xit? '
     Read(*,100) JobType

     Select case (JobType)
     Case ('p', 'P')

  	 	Write(*,110) PB_Gen, 'Generic problem (CST fields, MIL)'
  	 	Write(*,110) PB_MixedMode, 'Mixed mode'
  	 	Write(*,110) PB_Antiplane, 'Antiplane problem (CST fields, MIL) currently broken see l 350'  	 	
  	 	Write(*,110) PB_CylTwist_3D, 'Torsion of a 3D cylinder along the Z-axis'
  	 	Write(*,110) PB_Dipping, 'Unstable crack propagation, thermal loads'
  	 	Write(*,110) PB_Needleman, 'Unstable crack propagation, elastodynamic'
  	 	Write(*,100, advance = 'no') 'Problem type: ' 
  	 	Read(*,*) PB_Type

		Write(*,100, advance = 'no') 'Reading Coordinates / Connectivity'
        Select Case (Geom%Num_Dim)
        Case(2)
           Call Read_EXO_Node_Coord(Geom, Node2D_db, 1)
           Call Read_EXO_Connect(Geom, Elem2D_db)
        Case(3)
           Call Read_EXO_Node_Coord(Geom, Node3D_db, 1)
           Call Read_EXO_Connect(Geom, Elem3D_db)
        End Select

        Params%Sim_Str   = Trim(EXO_Str) // '_out'
        Geom%filename    = Trim(Params%Sim_Str) // '.gen'
        Params%PARAM_Str = Trim(Params%Sim_Str) // '.PARAM'
        Params%CST_Str   = Trim(Params%Sim_Str) // '.CST'

        Geom%exoid = EXCRE (Geom%filename, EXCLOB, exo_cpu_ws, exo_io_ws,  &
                & iErr)

        Call Ask_Rupt_Params(Geom, Params, BC_U_NS)

        Write(*,*) 'Updating geometrical information...'
        Call Write_EXO_Geom_Info(Geom)
 
        Select Case (Geom%Num_Dim)
        Case(2)
	   Write(*,100, advance = 'no') 'Writing Coordinates / Connectivity'
           Call Write_EXO_Node_Coord(Geom, Node2D_db, 1)
           Call Write_EXO_Connect(Geom, Elem2D_db)
           Write(*,*) 'Writing problem parameters'
           Call Write_Rupt_EXO_Params(Geom, Params)
           Call Write_Rupt_DATA(Geom, Params)    
               
           Write(*,*) 'Writing Boundary conditions'
           Call Write_BC2D(Geom, Params, Node2D_db, PB_Type, BC_U_NS)
           Write(*,*) 'Writing Forces'
           Call Write_Forces2D(Geom, Params, Node2D_db, PB_Type)
           Write(*,*) 'Writing Temperature field'
           Call Write_Temp2D(Geom, Params, Node2D_db, PB_Type)
        Case(3)
           Write(*,100, advance = 'no') 'Writing Coordinates / Connectivity'
           Call Write_EXO_Node_Coord(Geom, Node3D_db, 1)
           Call Write_EXO_Connect(Geom, Elem3D_db)
           Write(*,*) 'Writing problem parameters...'
           Call Write_Rupt_EXO_Params(Geom, Params)
           Call Write_Rupt_DATA(Geom, Params)        

           Write(*,*) 'Writing Boundary conditions'
           Call Write_BC3D(Geom, Params, Node3D_db, PB_Type, BC_U_NS)
           Write(*,*) 'Writing Forces'
           Call Write_Forces3D(Geom, Params, Node3D_db, PB_Type)
           Write(*,*) 'Writing Temperature field'
           Call Write_Temp3D(Geom, Params, Node3D_db, PB_Type)
        End Select

  
     Case ('c', 'C')
        Call Show_EXO_Geom_Info(Geom)
!        Select Case (Geom%Num_Dim)
!        Case(2)
!           Call Read_EXO_Node_Coord(Geom, Node2D_db, 1)
!           Call Read_EXO_Connect(Geom, Elem2D_db)
!        Case(3)
!           Call Read_EXO_Node_Coord(Geom, Node3D_db, 1)
!           Call Read_EXO_Connect(Geom, Elem3D_db)
!        End Select
     Case ('e', 'E')
        EXIT
     End Select
  End Do


100 Format(A)
110 Format('   [',I1,'] ',A)
!101 Format('Node ', I3, ' ID ', I3,' X: ', ES10.3)
!102 Format('Node ', I3, ' ID ', I3, ' X: ', ES10.3, ' Y: ', ES10.3)
!103 Format('Node ', I3, ' ID ', I3, ' X: ', ES10.3, ' Y: ', ES10.3,           &
!         & ' Z: ', ES10.3)
!200 Format('Elem ', I3, ' ID ', I3, '  DoF ', 16(I4))

Contains
  Subroutine Write_BC3D(Geom, Params, Node_db, PB_Type, BC_NS)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params)                            :: Params
    Type (Node3D), Dimension(:), Pointer          :: Node_db
    Integer, Intent(IN)                           :: PB_Type
    Type (Vect3D), Dimension(:), Pointer          :: BC_NS
    
    Integer                                       :: iTS, iNode, iSet
    Type(Vect3D), Dimension(:), Pointer           :: BC_Result, BC_Unit
    Real(Kind = Kr), Dimension(:), Pointer        :: V_Init
    Real(Kind = Kr)                               :: Theta

    !!! Generates BC vector from the TS informations and the
    !!! Load factor. 
    !!! Saves them as the result for the displacement
    Allocate (BC_Result(Geom%Num_Nodes))
    Allocate (BC_Unit(Geom%Num_Nodes))
    Allocate (V_Init(Geom%Num_Nodes))
    V_Init = 1.0_Kr
    
    BC_Unit%X = 0.0_Kr
    BC_Unit%Y = 0.0_Kr
    BC_Unit%Z = 0.0_Kr

    Do iSet = 1, Geom%Num_Node_Sets
       Print*, Geom%Node_Set(iSet)%Num_Nodes
       If (Params%BC_Type_X(iSet) == 1) Then
          Do iNode = 1, Geom%Node_Set(iSet)%Num_Nodes
             BC_Unit(Geom%Node_Set(iSet)%Node_ID(iNode))%X = BC_NS(iSet)%X
          End Do
       End If
       If (Params%BC_Type_Y(iSet) == 1) Then
          Do iNode = 1, Geom%Node_Set(iSet)%Num_Nodes
             BC_Unit(Geom%Node_Set(iSet)%Node_ID(iNode))%Y = BC_NS(iSet)%Y
          End Do
       End If
       If (Params%BC_Type_Z(iSet) == 1) Then
          Do iNode = 1, Geom%Node_Set(iSet)%Num_Nodes
             BC_Unit(Geom%Node_Set(iSet)%Node_ID(iNode))%Z = BC_NS(iSet)%Z
          End Do
       End If
       Print*, 'Params%BC_Type_X(iSet)', Params%BC_Type_X(iSet)
       Print*, 'Params%BC_Type_Y(iSet)', Params%BC_Type_Y(iSet)
       Print*, 'Params%BC_Type_Z(iSet)', Params%BC_Type_Z(iSet)
      Print*, 'Done generating BC_Unit for nodeset', iSet
    End Do
    
    Select Case (PB_Type)
    Case (PB_Gen)
       Do iTS = 1, Size(Params%Load)
          BC_Result%X = Params%Load(iTS) * BC_Unit%X
    	  BC_Result%Y = Params%Load(iTS) * BC_Unit%Y
          BC_Result%Z = Params%Load(iTS) * BC_Unit%Z
          
          Call Write_EXO_Result_Nodes(Geom, 2, iTS, BC_Result)
          Call Write_EXO_Result_Nodes(Geom, 1, iTS, V_Init)
       End Do
    Case (PB_CylTwist_3D)    
       Do iTS = 1, Size(Params%Load)
          Theta = Params%Load(iTS) * Pi / 180.0_Kr 
          Do iS = 1, Geom%Num_Nodes          
             BC_Result(iS)%X =                                                &
                  &    ( Node_db(iS)%Coord%X * Cos(Theta * BC_Unit(iS)%X)     &
                  &    - Node_db(iS)%Coord%Y * Sin(Theta * BC_Unit(iS)%X)     &
                  &  - Node_db(iS)%Coord%X )
             BC_Result(iS)%Y =                                                &
                  &    ( Node_db(iS)%Coord%X * Sin(Theta * BC_Unit(iS)%Y)     &
                  &    + Node_db(iS)%Coord%Y * Cos(Theta * BC_Unit(iS)%Y)     &
                  &  - Node_db(iS)%Coord%Y )
             BC_Result(iS)%Z = Params%Load(iTS) * BC_Unit(iS)%Z
          End Do
          
          Call Write_EXO_Result_Nodes(Geom, 2, iTS, BC_Result)
          Call Write_EXO_Result_Nodes(Geom, 1, iTS, V_Init)
       End Do
    Case (PB_MixedMode)
       Write(*,*) 'Theta? '
       Read*, Theta
       Theta = Theta * Pi / 180.0_Kr 
       Do iTS = 1, Size(Params%Load)
          Do iS = 1, Geom%Num_Nodes          
             BC_Result(iS)%X = Cos(Theta) * BC_Unit(iS)%X * Params%Load(iTS)
             BC_Result(iS)%Y = Sin(Theta) * BC_Unit(iS)%Y * Params%Load(iTS)
             BC_Result(iS)%Z = 0.0_Kr
          End Do
          
          Call Write_EXO_Result_Nodes(Geom, 2, iTS, BC_Result)
          Call Write_EXO_Result_Nodes(Geom, 1, iTS, V_Init)
       End Do
    Case Default
       Write(*,*) 'Unknown dimension / PB_Type combination in Write_BC3D'
    End Select
    
    DeAllocate (BC_Result, BC_Unit, V_Init)
  End Subroutine Write_BC3D

  Subroutine Write_BC2D(Geom, Params, Node_db, PB_Type, BC_NS)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params)                            :: Params
    Type (Node2D), Dimension(:), Pointer          :: Node_db
    Integer, Intent(IN)                           :: PB_Type
    Type(Vect3D), Dimension(:), Pointer           :: BC_NS
    
    Integer                                       :: iTS, iNode, iSet
    Type(Vect3D), Dimension(:), Pointer           :: BC_Result, BC_Unit
    Real(Kind = Kr), Dimension(:), Pointer        :: BC_2DA
    Real(Kind = Kr), Dimension(:), Pointer        :: V_Init
    Real(Kind = Kr)                               :: Theta
    
!!! Generates BC vector from the TS informations and the
!!! Load factor. 
!!! Saves them as the result for the displacement
    Allocate (BC_Result(Geom%Num_Nodes))
    Allocate (BC_Unit(Geom%Num_Nodes))
    Allocate (BC_2DA(Geom%Num_Nodes))
    Allocate (V_Init(Geom%Num_Nodes))

    V_Init = 1.0_Kr
    
    BC_Unit%X = 0.0_Kr
    BC_Unit%Y = 0.0_Kr
    BC_Unit%Z = 0.0_Kr
    BC_Result%X = 0.0_Kr
    BC_Result%Y = 0.0_Kr
    BC_Result%Z = 0.0_Kr
       
    Do iSet = 1, Geom%Num_Node_Sets
       If (Params%BC_Type_X(iSet) == 1) Then
          Do iNode = 1, Geom%Node_Set(iSet)%Num_Nodes
             BC_Unit(Geom%Node_Set(iSet)%Node_ID(iNode))%X = BC_NS(iSet)%X
          End Do
       End If
       If (Params%BC_Type_Y(iSet) == 1) Then
          Do iNode = 1, Geom%Node_Set(iSet)%Num_Nodes
             BC_Unit(Geom%Node_Set(iSet)%Node_ID(iNode))%Y = BC_NS(iSet)%Y
          End Do
       End If
       If (Params%BC_Type_Z(iSet) == 1) Then
          Do iNode = 1, Geom%Node_Set(iSet)%Num_Nodes
             BC_Unit(Geom%Node_Set(iSet)%Node_ID(iNode))%Z = BC_NS(iSet)%Z
          End Do
       End If
    End Do
    
    Select Case (PB_Type)
    Case (PB_Gen, PB_Dipping)
       Do iTS = 1, Size(Params%Load)
          BC_Result%X = Params%Load(iTS) * BC_Unit%X
    	  BC_Result%Y = Params%Load(iTS) * BC_Unit%Y
          
          Call Write_EXO_Result_Nodes(Geom, 2, iTS, BC_Result)
          Call Write_EXO_Result_Nodes(Geom, 1, iTS, V_Init)
       End Do
    Case (PB_MixedMode)
       Write(*,*) 'Theta? '
       Read*, Theta
       Theta = Theta * Pi / 180.0_Kr 
       Do iTS = 1, Size(Params%Load)
          Do iS = 1, Geom%Num_Nodes          
             BC_Result(iS)%X = Cos(Theta) * BC_Unit(iS)%X * Params%Load(iTS)
             BC_Result(iS)%Y = Sin(Theta) * BC_Unit(iS)%Y * Params%Load(iTS)
          End Do
          
          Call Write_EXO_Result_Nodes(Geom, 2, iTS, BC_Result)
          Call Write_EXO_Result_Nodes(Geom, 1, iTS, V_Init)
       End Do
    Case (PB_Antiplane)
       Do iTS = 1, Size(Params%Load)
          BC_Result%Z = Params%Load(iTS) * BC_Unit%Z
!          BC_2DA = Params%Load(iTS) * BC_Unit%Z
          
          Call Write_EXO_Result_Nodes(Geom, 2, iTS, BC_Result)
          Call Write_EXO_Result_Nodes(Geom, 1, iTS, V_Init)
       End Do
    Case Default
       Write(*,*) 'Unknown dimension / PB_Type combination in Write_BC2D'
    End Select
    
    DeAllocate (BC_Result, BC_Unit, V_Init)
  End Subroutine Write_BC2D
  
!  Subroutine Write_BC2DA(Geom, Params, Node_db, PB_Type, BC_NS)
!    Type (EXO_Geom_Info)                          :: Geom
!    Type (Rupt_Params)                            :: Params
!    Type (Node2D), Dimension(:), Pointer          :: Node_db
!    Integer, Intent(IN)                           :: PB_Type
!    Type(Vect3D), Dimension(:), Pointer           :: BC_NS
!    
!    Integer                                       :: iTS, iNode, iSet
!    Real(Kind = Kr), Dimension(:), Pointer        :: BC_Result, BC_Unit
!    Real(Kind = Kr), Dimension(:), Pointer        :: V_Init
!    Real(Kind = Kr)                               :: Theta
!
!    !!! Generates BC vector from the TS informations and the
!    !!! Load factor. 
!    !!! Saves them as the result for the displacement
!    Allocate (BC_Result(Geom%Num_Nodes))
!    Allocate (BC_Unit(Geom%Num_Nodes))
!    Allocate (V_Init(Geom%Num_Nodes))
!    V_Init = 1.0_Kr
!    
!    BC_Unit = 0.0_Kr
!
!    Do iSet = 1, Geom%Num_Node_Sets
!       If (Params%BC_Type_Z(iSet) == 1) Then
!          Do iNode = 1, Geom%Node_Set(iSet)%Num_Nodes
!             BC_Unit(Geom%Node_Set(iSet)%Node_ID(iNode)) = BC_NS(iSet)%Z
!          End Do
!       End If
!    End Do
!    
!	Select Case (PB_Type)
!	Case (PB_Antiplane)
!	   Do iTS = 1, Size(Params%Load)
!	      BC_Result = Params%Load(iTS) * BC_Unit
!           
!          Call Write_EXO_Result_Nodes(Geom, 2, iTS, BC_Result)
!          Call Write_EXO_Result_Nodes(Geom, 1, iTS, V_Init)
!       End Do
!	Case Default
!	   Write(*,*) 'Unkown dimension / PB_Type combination in Write_BC2DA'
!	End Select
!
!    DeAllocate (BC_Result, BC_Unit, V_Init)
!  End Subroutine Write_BC2DA

  Subroutine Write_Forces3D(Geom, Params, Node_db, PB_Type)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params)                            :: Params
    Type (Node3D), Dimension(:), Pointer          :: Node_db
    Integer, Intent(IN)                           :: PB_Type
    
    Integer                                       :: iTS, iNode, iSet
    Real(Kind = Kr), Dimension(:), Pointer        :: F_Result, F_Unit
    Type(Vect3D)                                  :: Force
    Real(Kind = Kr), Dimension(:), Pointer        :: V_Init
    Real(Kind = Kr)                               :: Theta

    !!! Generates BC vector from the TS informations and the
    !!! Load factor. 
    !!! Saves them as the result for the displacement
    Allocate (F_Result(Geom%Num_Nodes * Geom%Num_Dim))
    Allocate (F_Unit(Geom%Num_Nodes * Geom%Num_Dim))
    F_Unit = 0.0_Kr

	Select Case (PB_Type)
	Case (PB_Gen)
	   Write (*,*) 'Force: '
	   Read (*,*) Force
	   F_Unit(1:Geom%Num_Nodes * Geom%Num_Dim:Geom%Num_Dim) = Force%X
	   F_Unit(2:Geom%Num_Nodes * Geom%Num_Dim:Geom%Num_Dim) = Force%Y
	   F_Unit(3:Geom%Num_Nodes * Geom%Num_Dim:Geom%Num_Dim) = Force%Z
    
       Do iTS = 1, Size(Params%Load)
		  F_Result = Params%Load(iTS) * F_Unit
          Call Write_EXO_Result_Nodes(Geom, 5, iTS, F_Result)
       End Do
    Case Default
	   Write(*,*) 'Unkown dimension / PB_Type combination in Write_Forces3D'
	End Select
    DeAllocate (F_Result, F_Unit)
  End Subroutine Write_Forces3D

  Subroutine Write_Forces2D(Geom, Params, Node_db, PB_Type)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params)                            :: Params
    Type (Node2D), Dimension(:), Pointer          :: Node_db
    Integer, Intent(IN)                           :: PB_Type
    
    Integer                                       :: iTS, iNode, iSet
    Real(Kind = Kr), Dimension(:), Pointer        :: F_Result, F_Unit
    Type(Vect2D)                                  :: Force
    Real(Kind = Kr), Dimension(:), Pointer        :: V_Init
    Real(Kind = Kr)                               :: Theta

    !!! Generates BC vector from the TS informations and the
    !!! Load factor. 
    !!! Saves them as the result for the displacement
    Allocate (F_Result(Geom%Num_Nodes * Geom%Num_Dim))
    Allocate (F_Unit(Geom%Num_Nodes * Geom%Num_Dim))
    F_Unit = 0.0_Kr

	Select Case (PB_Type)
	Case (PB_Gen)
	   Write (*,*) 'Force: '
	   Read (*,*) Force
	   F_Unit(1:Geom%Num_Nodes * Geom%Num_Dim:Geom%Num_Dim) = Force%X
	   F_Unit(2:Geom%Num_Nodes * Geom%Num_Dim:Geom%Num_Dim) = Force%Y
    
       Do iTS = 1, Size(Params%Load)
		  F_Result = Params%Load(iTS) * F_Unit
          Call Write_EXO_Result_Nodes(Geom, 5, iTS, F_Result)
       End Do
    Case Default
	   Write(*,*) 'Unkown dimension / PB_Type combination in Write_Force2D'
	End Select
    DeAllocate (F_Result, F_Unit)
  End Subroutine Write_Forces2D
  
  Subroutine Write_Temp3D(Geom, Params, Node_db, PB_Type)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params)                            :: Params
    Type (Node3D), Dimension(:), Pointer          :: Node_db
    Integer, Intent(IN)                           :: PB_Type
    
    Integer                                       :: iTS, iNode, iSet
    Real(Kind = Kr), Dimension(:), Pointer        :: Temp_Result, Temp_Unit
    Real(Kind = Kr)                               :: Temp
    Real(Kind = Kr), Dimension(:), Pointer        :: V_Init
    Real(Kind = Kr)                               :: Theta

    !!! Generates BC vector from the TS informations and the
    !!! Load factor. 
    !!! Saves them as the result for the displacement
    Allocate (Temp_Result(Geom%Num_Nodes))
    Allocate (Temp_Unit(Geom%Num_Nodes))
    Temp_Unit = 0.0_Kr

	Select Case (PB_Type)
	Case (PB_Gen)
	   Write (*,*) 'Temperature: '
	   Read (*,*) Temp
	   Temp_Unit = Temp
    
       Do iTS = 1, Size(Params%Load)
		  Temp_Result = Params%Load(iTS) * Temp_Unit
          Call Write_EXO_Result_Nodes(Geom, 8, iTS, Temp_Result)
       End Do
    Case Default
	   Write(*,*) 'Unkown dimension / PB_Type combination in Write_Temp3D'
	End Select
    DeAllocate (Temp_Result, Temp_Unit)
  End Subroutine Write_Temp3D

  Subroutine Write_Temp2D(Geom, Params, Node_db, PB_Type)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params)                            :: Params
    Type (Node2D), Dimension(:), Pointer          :: Node_db
    Integer, Intent(IN)                           :: PB_Type
    
    Integer                                       :: iTS, iNode, iSet
    Real(Kind = Kr), Dimension(:), Pointer        :: Temp_Result, Temp_Unit
    Real(Kind = Kr)                               :: Temp
    Real(Kind = Kr), Dimension(:), Pointer        :: V_Init
    Real(Kind = Kr)                               :: Theta

	Real(Kind = Kr)                               :: TCool, THot
	Real(Kind = Kr)                               :: b, V, D, P
	Real(Kind = Kr)                               :: Y
	

    !!! Generates BC vector from the TS informations and the
    !!! Load factor. 
    !!! Saves them as the result for the displacement
    Allocate (Temp_Result(Geom%Num_Nodes))
    Allocate (Temp_Unit(Geom%Num_Nodes))
    Temp_Unit = 0.0_Kr

	Select Case (PB_Type)
	Case (PB_Gen)
	   Write (*,*) 'Temperature: '
	   Read (*,*) Temp
	   Temp_Unit = Temp
    
       Do iTS = 1, Size(Params%Load)
		  Temp_Result = Params%Load(iTS) * Temp_Unit
          Call Write_EXO_Result_Nodes(Geom, 8, iTS, Temp_Result)
       End Do
    Case (PB_Dipping)
       Write(*, 100, advance = 'no') 'TCool, THot:                '
       Read(*,*) TCool, THot
       Write(*, 100, advance = 'no') 'Half plane width (b):       '
       Read(*,*) b
       Write(*, 100, advance = 'no') 'Heat diffusion coefficient: '
       Read(*,*) D
       Write(*, 100, advance = 'no') 'Dipping speed:              '
       Read(*,*) V
       P = V/D*b
       Do iTS = 1, Size(Params%Load)
          Temp_Result = TCool
          Do iNode = 1, Geom%Num_Nodes
             Y = Node_db(iNode)%Coord%Y
             If (Y >= Params%Load(iTS)) Then
!!$             	Temp_Result(iNode) = TCool + (THot - TCool) *                 &
!!$             		& (1.0_Kr - exp(-P*(Y-Params%Load(iTS) )) )
             	Temp_Result(iNode) = THot + (TCool - THot) *                  &
             		& exp( -P*(Y-Params%Load(iTS)) ) 
             End If
          End Do
          Call Write_EXO_Result_Nodes(Geom, 8, iTS, Temp_Result)
       End Do          
    Case Default
	   Write(*,*) 'Unkown dimension / PB_Type combination in Write_Temp2D'
	End Select
    DeAllocate (Temp_Result, Temp_Unit)
  100 Format(A)
  End Subroutine Write_Temp2D

  Subroutine Ask_Rupt_Params(Geom, Params, BC_NS)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params)                            :: Params
    Type(Vect3D), Dimension(:), Pointer           :: BC_NS
    
    Character                                     :: Char_Input
    Integer                                       :: Num_TS
    Real(Kind = Kr)                               :: TS_init, TS_Final
    Real(Kind = Kr)                               :: Tmp_DistFactor
    Integer                                       :: i, iBlk, iSet
    Integer                                       :: Blk_ID, Set_ID


    !!! The following are decent values.
    Params%nbCracks           = 0
    Params%MaxCrackLength     = 0.0             
    Params%MaxIterRelax       = 5000
    Params%TolKSP             = 1.0D-6
    Params%TolRelax           = 1.0D-3

    !!! This is totally meaningful in general...
    Params%Epsilon      = 1.0D-1
    Params%KEpsilon     = 1.0D-4
    Params%TolIrrev     = 1.0D-2
    Write(*,100) '===== Global Properties ======'
!    Do
!       Write(*,100, advance = 'no')                                           &
!            & 'Problem Type [D]isplacement, [T]ermal loads: '
!       Read(*,100) Char_Input
!       Select case (Char_Input)
!       Case('d', 'D')
!          Params%PB_Type = PB_DISPL
!          EXIT
!       Case('t', 'T')
!          Params%PB_Type = PB_THERM
!          EXIT
!       End Select
!    End Do

    Write(*, 100, advance = 'no') 'Number of time steps: '
    Read(*,*) Num_TS
    Allocate(Params%Load(Num_TS))
    Write(*, 100, advance = 'no') 'Initial load:         '
    Read(*,*) TS_Init
    If (Num_TS ==1) Then
       Params%Load = TS_INIT
    Else
       Write(*, 100, advance = 'no') 'Final load:           '
       Read(*,*) TS_Final
       Params%Load = (/ (TS_Init + (i-1.0)*(TS_Final - TS_Init)            &
            & / (Num_TS - 1.0_Kr), i = 1,Num_TS) /)
    End If!    Print*, 'Load vector: ', Params%Load
    Write(*, 100, advance = 'no')'Do_BackTrack [T/F]          '    
    Read(*,200) Params%Do_BackTrack
    If (Params%Do_BackTrack) Then
     	Write(*, 100, advance = 'no') 'BackTrack Tol : '
	Read(*,*) Params%Tol_Ener
    End If
    Write(*, 100, advance = 'no') 'Type of Crack initialization [0=PREV 1=1, 2=RND]: '
    Read(*,*) Params%Init_V   
    If (Params%Init_V==2) Then         
        Write(*, 100, advance = 'no') 'Number of Cracks : '
	Read(*,*) Params%nbCracks
	Write(*, 100, advance = 'no') 'MaxLength of Cracks : '
	Read(*,*) Params%MaxCrackLength
    End If
    Write(*,100) '===== Element block Properties ======'
    Allocate (Params%Is_Brittle(Geom%Num_Elem_Blks))
    Allocate (Params%Is_Domain(Geom%Num_Elem_Blks))
    Allocate (Params%Has_Force(Geom%Num_Elem_Blks))
    Allocate (Params%Toughness(Geom%Num_Elem_Blks))
    Allocate (Params%Young_Mod(Geom%Num_Elem_Blks))
    Allocate (Params%Poisson_Ratio(Geom%Num_Elem_Blks))
    Allocate (Params%Therm_Exp(Geom%Num_Elem_Blks))

    Do iBlk = 1, Geom%Num_Elem_Blks
       Blk_ID = Geom%Elem_Blk(iBlk)%ID
       Write(*,300) Blk_ID
       Write(*,100, advance = 'no') 'Is_Brittle [T/F]          '
       Read(*,200) Params%Is_Brittle(iBlk)
       Write(*,100, advance = 'no') 'Is_Domain [T/F]           '
       Read(*,200) Params%Is_Domain(iBlk)
       Write(*,100, advance = 'no') 'Has_Force [T/F]           '
       Read(*,200) Params%Has_Force(iBlk)
       If ( Params%Is_Brittle(iBLK) .AND. Params%Has_Force(iBlk)) Then
          Print*, 'WARNING: Element block has force and is Brittle'
       End If
!       If ( Params%Has_Force(iBlk) ) Then
!          Write(*,100, advance = 'no') 'Force X:                  '
!          Read(*,*) Params%Force(iBlk)%X
!          Write(*,100, advance = 'no') 'Force Y:                  '
!          Read(*,*) Params%Force(iBlk)%Y
!          Write(*,100, advance = 'no') 'Force Z:                  '
!          Read(*,*) Params%Force(iBlk)%Z
!       Else
!          Params%Force(iBlk)%X = 0.0_Kr
!          Params%Force(iBlk)%Y = 0.0_Kr
!          Params%Force(iBlk)%Z = 0.0_Kr
!       End If
       Write(*,100, advance = 'no') 'Toughness:                '
       Read(*,*) Params%Toughness(iBlk)
       Write(*,100, advance = 'no') 'Young''s modulus           '
       Read(*,*) Params%Young_Mod(iBlk)
       Write(*,100, advance = 'no') 'Poisson ratio             ' 
       Read(*,*) Params%Poisson_Ratio(iBlk)
       Write(*,100, advance = 'no') 'Thermal Expansion coef.   ' 
       Read(*,*) Params%Therm_Exp(iBlk)
    End Do

    Write(*,100) '===== Node sets Properties ========='
    Allocate (Params%BC_Type_X(Geom%Num_Node_Sets))
    Allocate (Params%BC_Type_Y(Geom%Num_Node_Sets))
    Allocate (Params%BC_Type_Z(Geom%Num_Node_Sets))
    Allocate (BC_NS(Geom%Num_Node_Sets))

    BC_NS%X = 0.0_Kr
    BC_NS%Y = 0.0_Kr
    BC_NS%Z = 0.0_Kr
    Do iSet = 1, Geom%Num_Node_Sets
       Set_ID = Geom%Node_Set(iSet)%ID
       Write(*,400) Set_ID
       Write(*,100, advance = 'no') 'BC U type, X direction (None=0 DIRI = 1) '

       Read(*,*) Params%BC_Type_X(iSet) 
       If (Params%BC_Type_X(iSet) == 1) Then
          Write(*,100, advance = 'no') 'BC U, X direction                     '
          Read(*,*) BC_NS(iSet)%X
       End If
       
       Write(*,100,advance = 'no') 'BC U type, Y direction (None=0 DIRI = 1) '
       Read(*,*) Params%BC_Type_Y(iSet) 
       If (Params%BC_Type_Y(iSet) == 1) Then
          Write(*,100, advance = 'no') 'BC U, Y direction                    '
          Read(*,*) BC_NS(iSet)%Y
       End If
       
!       If (Geom%Num_Dim == 3) Then
          Write(*,100,advance = 'no') 'BC U type, Z direction (None=0 DIRI = 1) '
          Read(*,*) Params%BC_Type_Z(iSet) 
          If (Params%BC_Type_Z(iSet) == 1) Then
             Write(*,100, advance = 'no') 'BC U, Z direction                    '
             Read(*,*) BC_NS(iSet)%Z 
          End If
!       End If
    End Do
    
100 Format(A)
200 Format(L1)
300 Format('=== Element Block             ', I3)
400 Format('=== Node set                  ', I3)


  End Subroutine Ask_Rupt_Params
    
!!$
!!$
!!$!!$  Call Show_EXO_Geom_Info(Geom)
!!$
!!$  Num_dof_per_nodes = 1
!!$  Num_dof_per_nodes = Geom%Num_Dim
!!$
!!$  Select Case (Geom%Num_Dim)
!!$  Case(1)
!!$     Print*, 'Calling Read_EXO_Node with Node1D_db'
!!$     Call Read_EXO_Node_Coord(Geom, Node1D_db, Num_dof_per_nodes)
!!$     Do iS = 1, Size(Node1D_db)
!!$        Write(*,101) iS, Node1D_db(iS)%ID, Node1D_db(iS)%Coord
!!$     End Do
!!$  Case(2)
!!$     Print*, 'Calling Read_EXO_Node with Node2D_db'
!!$     Call Read_EXO_Node_Coord(Geom, Node2D_db, Num_dof_per_nodes)
!!$     Call Read_Exo_Connect(Geom, Elem2D)
!!$     Write(*,*) 'Nodes'
!!$     Do iS = 1, Size(Node2D_db)
!!$        Write(*,102) iS, Node2D_db(iS)%ID, Node2D_db(iS)%Coord
!!$     End Do
!!$     Write(*,*) 'Elements'
!!$     Do iE = 1, Geom%Num_Elems
!!$        Write(*,200) iE, Elem2D(iE)%ID_EL, Elem2D(iE)%ID_DoF
!!$     End Do
!!$  Case(3)
!!$     Print*, 'Calling Read_EXO_Node with Node3D_db'
!!$     Call Read_EXO_Node_Coord(Geom, Node3D_db, Num_dof_per_nodes)
!!$     Print*, 'Calling Read_EXO_Connect with Elem3D_Elast'
!!$     Call Read_Exo_Connect(Geom, Elem3D_Elast)
!!$     Do iS = 1, Size(Node3D_db)
!!$        Write(*,103) iS, Node3D_db(iS)%ID, Node3D_db(iS)%Coord
!!$     End Do
!!$     Write(*,*) 'Elements'
!!$     Do iE = 1, Geom%Num_Elems
!!$        Write(*,200) iE, Elem3D_Elast(iE)%ID_EL, Elem3D_Elast(iE)%ID_DoF
!!$     End Do
!!$  End Select
!!$  
!!$
!!$  Geom%filename = Trim(Geom%Filename)//'.new'
!!$  Call Write_EXO_Geom_Info(Geom)
!!$  Call Destroy_EXO_Geom_Info(Geom)

End Program Prep_Rupt
