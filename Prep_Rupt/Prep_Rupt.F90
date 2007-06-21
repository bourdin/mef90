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

  Type (Rupt_Params2D)                            :: Params2D
  Type (Rupt_Params3D)                            :: Params3D
  
  Character(len = MXLNLN)                         :: EXO_Str

  Integer                                         :: Num_dof_per_nodes 
  Integer                                         :: iErr

  Integer                                         :: iB, iE, iV, iG
  Integer                                         :: iSL1, iSL2, iSG1, iSG2
  

  Real(Kind = Kr), Dimension(:), Pointer          :: Res
  Real(Kind = Kr)                                 :: Vol
  
  Integer, Parameter                              :: PB_Gen         = 0
  Integer, Parameter                              :: PB_MixedMode   = 1
  Integer, Parameter                              :: PB_Antiplane   = 2
  Integer, Parameter                              :: PB_CylTwist_3D = 3
  Integer, Parameter                              :: PB_Dipping     = 4
  Integer, Parameter                              :: PB_Dipping2    = 5
  Integer, Parameter                              :: PB_Dipping3    = 6
  Integer, Parameter                              :: PB_Needleman   = 7
  Integer, Parameter                              :: PB_KBMultiscale= 8
  Integer                                         :: PB_Type

  Write(*,100, advance = 'no') 'ExodusII file prefix: '
  Read(*,100) EXO_Str

  Geom%filename = Trim(EXO_Str)//'.gen'
  Call Read_EXO_Geom_Info(Geom) 
  
  Do
     Write(*,100, advance = 'no') '[P]repare [C]heck, or [E]xit?'
     Read(*,100) JobType

     Select case (JobType)
     Case ('p', 'P')

  	 	Write(*,110) PB_Gen,         'Generic problem (CST fields, MIL)'
  	 	Write(*,110) PB_MixedMode,   'Mixed mode'
  	 	Write(*,110) PB_Antiplane,   'Antiplane problem (CST fields, MIL)'  	 	
  	 	Write(*,110) PB_CylTwist_3D, 'Torsion of a 3D cylinder along the Z-axis'
  	 	Write(*,110) PB_Dipping,     'Unstable crack propagation, thermal loads'
  	 	Write(*,110) PB_Dipping2,    'Unstable crack propagation, thermal loads: Rescaled'
  	 	Write(*,110) PB_Dipping3,    'Unstable crack propagation, thermal loads: First time step'
  	 	Write(*,110) PB_Needleman,   'Unstable crack propagation, elastodynamic'
  	 	Write(*,110) PB_KBMultiScale,'Multi scale experiment on orthothropic crystals'
  	 	Write(*,100, advance = 'no') 'Problem type' 
  	 	Read(*,*) PB_Type

		Write(*,*) 'Reading Coordinates / Connectivity'
        Select Case (Geom%Num_Dim)
        Case(2)
           Call Read_EXO_Node_Coord(Geom, Node2D_db, 1)
           Call Read_EXO_Connect(Geom, Elem2D_db)
           Params2D%Sim_Str   = Trim(EXO_Str) // '_out'
           Geom%filename      = Trim(Params2D%Sim_Str) // '.gen'
           Params2D%PARAM_Str = Trim(Params2D%Sim_Str) // '.PARAM'
           Params2D%CST_Str   = Trim(Params2D%Sim_Str) // '.CST'
           Geom%exoid = EXCRE (Geom%filename, EXCLOB, exo_cpu_ws, exo_io_ws, iErr)
           Call Ask_Rupt_Params2D(Geom, Params2D, BC_U_NS)
           Write(*,*) 'Updating geometrical information...'
           Call Write_EXO_Geom_Info(Geom)

  	        Write(*,100, advance = 'no') 'Writing Coordinates / Connectivity'
           Call Write_EXO_Node_Coord(Geom, Node2D_db, 1)
           Call Write_EXO_Connect(Geom, Elem2D_db)
           Write(*,*) 'Writing problem parameters'
           Call Write_Rupt_EXO_Params(Geom, Params2D)
           Print*, 'OK Params'
           Call Write_Rupt_DATA(Geom, Params2D)    
            
           Write(*,*) 'Writing Boundary conditions'
           Call Write_BC2D(Geom, Params2D, Node2D_db, PB_Type, BC_U_NS)
           Write(*,*) 'Writing Forces'
           Call Write_Forces2D(Geom, Params2D, Node2D_db, PB_Type)
           Write(*,*) 'Writing Temperature field'
           Call Write_Temp2D(Geom, Params2D, Node2D_db, PB_Type)
        Case(3)
           Call Read_EXO_Node_Coord(Geom, Node3D_db, 1)
           Call Read_EXO_Connect(Geom, Elem3D_db)
           Params3D%Sim_Str   = Trim(EXO_Str) // '_out'
           Geom%filename      = Trim(Params3D%Sim_Str) // '.gen'
           Params3D%PARAM_Str = Trim(Params3D%Sim_Str) // '.PARAM'
           Params3D%CST_Str   = Trim(Params3D%Sim_Str) // '.CST'
           Geom%exoid = EXCRE (Geom%filename, EXCLOB, exo_cpu_ws, exo_io_ws, iErr)
           Call Ask_Rupt_Params3D(Geom, Params3D, BC_U_NS)
           Write(*,*) 'Updating geometrical information...'
           Call Write_EXO_Geom_Info(Geom)

           Write(*,100, advance = 'no') 'Writing Coordinates / Connectivity'
           Call Write_EXO_Node_Coord(Geom, Node3D_db, 1)
           Call Write_EXO_Connect(Geom, Elem3D_db)
           Write(*,*) 'Writing problem parameters...'
           Call Write_Rupt_EXO_Params(Geom, Params3D)
           Call Write_Rupt_DATA(Geom, Params3D)        

           Write(*,*) 'Writing Boundary conditions'
           Call Write_BC3D(Geom, Params3D, Node3D_db, PB_Type, BC_U_NS)
           Write(*,*) 'Writing Forces'
           Call Write_Forces3D(Geom, Params3D, Node3D_db, PB_Type)
           Write(*,*) 'Writing Temperature field'
           Call Write_Temp3D(Geom, Params3D, Node3D_db, PB_Type)
        End Select

     Case ('c', 'C')
        Call Show_EXO_Geom_Info(Geom)
     Case ('e', 'E')
        Write(*,*)
        EXIT
     End Select
  End Do


100 Format(A,T78,': ')
110 Format('   [',I1,'] ',A)

Contains
  Subroutine Write_BC3D(Geom, Params, Node_db, PB_Type, BC_NS)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params3D)                          :: Params
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
    Type (Rupt_Params2D)                          :: Params
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
    Case (PB_Gen, PB_Dipping, PB_Dipping2, PB_Dipping3, PB_KBMultiscale)
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
          
          Call Write_EXO_Result_Nodes(Geom, 2, iTS, BC_Result)
          Call Write_EXO_Result_Nodes(Geom, 1, iTS, V_Init)
       End Do
    Case Default
       Write(*,*) 'Unknown dimension / PB_Type combination in Write_BC2D'
    End Select
    
    DeAllocate (BC_Result, BC_Unit, V_Init)
  End Subroutine Write_BC2D
  
  Subroutine Write_Forces3D(Geom, Params, Node_db, PB_Type)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params3D)                          :: Params
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
	   Write (*,100, advance = 'no') 'Force%X: '
	   Read (*,*) Force%X
	   Write (*,100, advance = 'no') 'Force%Y: '
	   Read (*,*) Force%Y
	   Write (*,100, advance = 'no') 'Force%Z: '
	   Read (*,*) Force%Z
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
100 Format (A,T78,': ')
  End Subroutine Write_Forces3D

  Subroutine Write_Forces2D(Geom, Params, Node_db, PB_Type)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params2D)                          :: Params
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
	   Write (*,100, advance = 'no') 'Force%X: '
	   Read (*,*) Force%X
	   Write (*,100, advance = 'no') 'Force%Y: '
	   Read (*,*) Force%Y
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
100 Format (A,T78,': ')
  End Subroutine Write_Forces2D
  
  Subroutine Write_Temp3D(Geom, Params, Node_db, PB_Type)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params3D)                          :: Params
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
    Type (Rupt_Params2D)                          :: Params
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
	Real(Kind = Kr)                               :: l, DTheta, kappa
	

    !!! Generates BC vector from the TS informations and the
    !!! Load factor. 
    !!! Saves them as the result for the displacement
    Allocate (Temp_Result(Geom%Num_Nodes))
    Allocate (Temp_Unit(Geom%Num_Nodes))
    Temp_Unit = 0.0_Kr

	Select Case (PB_Type)
	Case (PB_Gen)
	   Write (*,100, advance = 'no') 'Temperature: '
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
             	Temp_Result(iNode) = THot + (TCool - THot) * exp( -P*(Y-Params%Load(iTS)) ) 
             End If
          End Do
          Call Write_EXO_Result_Nodes(Geom, 8, iTS, Temp_Result)
       End Do          
    Case (PB_Dipping2)
       Write(*, 100, advance = 'no') 'Temperature difference (Delta Theta) '
       Read(*,*) DTheta
       Write(*, 100, advance = 'no') 'Characteristic length (l)            '
       Read(*,*) l
       Write(*, 100, advance = 'no') 'Thermal conductivity (kappa)         '
       Read(*,*) kappa
       Write(*, 100, advance = 'no') 'Quenching speed (V)                  '
       Read(*,*) V
       Do iTS = 1, Size(Params%Load)
       Temp_Result = -DTheta
          Do iNode = 1, Geom%Num_Nodes
             Y = Node_db(iNode)%Coord%Y
             If (Y >= Params%Load(iTS)) Then
             	Temp_Result(iNode) = -DTheta * exp( -V/kappa*l*(Y-Params%Load(iTS)))
             End If
          End Do
          Call Write_EXO_Result_Nodes(Geom, 8, iTS, Temp_Result)
       End Do          
    Case (PB_Dipping3)
       Write(*, 100, advance = 'no') 'Characteristic length (l)            '
       Read(*,*) l
       Write(*, 100, advance = 'no') 'Thermal conductivity (kappa)         '
       Read(*,*) kappa
       Do iTS = 1, Size(Params%Load)
       Temp_Result = 0.0_Kr
          Do iNode = 1, Geom%Num_Nodes
             Y = Node_db(iNode)%Coord%Y
             If (Y == 0.0_Kr) Then
             	Temp_Result(iNode) = -Params%Load(iTS)
             End If
          End Do
          Call Write_EXO_Result_Nodes(Geom, 8, iTS, Temp_Result)
       End Do          
    Case Default
	   Write(*,*) 'Unkown dimension / PB_Type combination in Write_Temp2D'
	End Select
    DeAllocate (Temp_Result, Temp_Unit)
  100 Format(A,T78,': ')
  End Subroutine Write_Temp2D

  Subroutine Ask_Rupt_Params2D(Geom, Params, BC_NS)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params2D)                          :: Params
    Type (Vect3D), Dimension(:), Pointer          :: BC_NS
    
    Character                                     :: Char_Input
    Integer                                       :: Num_TS
    Real(Kind = Kr)                               :: TS_init, TS_Final
    Real(Kind = Kr)                               :: Tmp_DistFactor
    Integer                                       :: i, iBlk, iSet
    Integer                                       :: Blk_ID, Set_ID

    Real(Kind = Kr)                               :: Tmp_E, Tmp_nu
    Real(Kind = Kr)                               :: Tmp_Lambda, Tmp_mu1, Tmp_mu2
    Real(Kind = Kr)                               :: Tmp_Toughness, Tmp_Therm_Exp, Tmp_theta
    Integer                                       :: Num_Grains

    !!! The following are decent values.
    Params%nbCracks           = 0
    Params%MaxCrackLength     = 0.0             
    Params%MaxIterRelax       = 5000
    Params%TolKSP             = 1.0D-6
    Params%TolRelax           = 1.0D-3
    Params%Do_Irrev           = .TRUE.

    !!! This is totally meaningful in general...
    Params%Epsilon      = 1.0D-1
    Params%KEpsilon     = 1.0D-6
    Params%TolIrrev     = 1.0D-2
    Write(*,101) '===== Global Properties ======'

    Write(*, 100, advance = 'no') 'Number of time steps'
    Read(*,*) Num_TS
    Allocate(Params%Load(Num_TS))
    Write(*, 100, advance = 'no') 'Initial load'
    Read(*,*) TS_Init
    If (Num_TS ==1) Then
       Params%Load = TS_INIT
    Else
       Write(*, 100, advance = 'no') 'Final load'
       Read(*,*) TS_Final
       Params%Load = (/ (TS_Init + (i-1.0)*(TS_Final - TS_Init)  / (Num_TS - 1.0_Kr), i = 1,Num_TS) /)
    End If
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
    Write(*,101) '===== Element block Properties ======'
    Allocate (Params%Is_Brittle(Geom%Num_Elem_Blks))
    Allocate (Params%Is_Domain(Geom%Num_Elem_Blks))
    Allocate (Params%Has_Force(Geom%Num_Elem_Blks))
    Allocate (Params%Toughness(Geom%Num_Elem_Blks))
    Allocate (Params%Hookes_Law(Geom%Num_Elem_Blks))
    Allocate (Params%Therm_Exp(Geom%Num_Elem_Blks))
 
    Num_Grains = 0
    If (PB_Type == PB_KBMultiscale) Then
       Write(*,500)
       Write(*,100, advance = 'no') 'Number of grains:'
       Read(*,*) Num_Grains
       Write(*,100, advance = 'no') 'Toughness:                '
       Read(*,*) Tmp_Toughness
       Write(*,100, advance = 'no') 'Lambda'
       Read(*,*) tmp_Lambda      
       Write(*,100, advance = 'no') 'mu1'
       Read(*,*) tmp_mu1      
       Write(*,100, advance = 'no') 'mu2'
       Read(*,*) tmp_mu2      
       Write(*,100, advance = 'no') 'Thermal Expansion coef.   ' 
       Read(*,*) Tmp_Therm_Exp
       Call Random_Seed
       Do iBlk = 1, Num_Grains
          Params%Is_Brittle(iBlk) = .TRUE.
          Params%Is_Domain(iBlk)  = .TRUE.
          Params%Has_Force(iBlk)  = .FALSE.
          Call Random_Number(tmp_theta)
          tmp_theta = tmp_theta * pi
          Call GenHL_Ortho2D_LambdaMu(Tmp_Lambda, Tmp_mu1, Tmp_mu2, Tmp_theta, Params%Hookes_Law(iBlk))
          Params%Toughness(iBlk) = Tmp_Toughness
          Params%Therm_Exp(iBlk) = Tmp_Therm_Exp
          Write(*,600) iBlk, Tmp_Theta * 180.0_Kr / pi
       End Do
    End If
    Do iBlk = Num_Grains + 1, Geom%Num_Elem_Blks
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
       Write(*,100, advance = 'no') 'Toughness:                '
       Read(*,*) Params%Toughness(iBlk)
       Write(*,100, advance = 'no') 'Young''s modulus           '
       Read(*,*) Tmp_E
       Write(*,100, advance = 'no') 'Poisson ratio             ' 
       Read(*,*) Tmp_nu
       Write(*,100, advance = 'no') 'Thermal Expansion coef.   ' 
       Read(*,*) Params%Therm_Exp(iBlk)
       Call GenHL_Iso2D_EnuPlaneStress(Tmp_E, Tmp_nu, Params%Hookes_Law(iBlk))
    End Do

    Write(*,101) '===== Node sets Properties ========='
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
       
       Write(*,100,advance = 'no') 'BC U type, Z direction (None=0 DIRI = 1) '
       Read(*,*) Params%BC_Type_Z(iSet) 
       If (Params%BC_Type_Z(iSet) == 1) Then
          Write(*,100, advance = 'no') 'BC U, Z direction                    '
          Read(*,*) BC_NS(iSet)%Z 
       End If
    End Do
    
100 Format(A,T78,': ')
101 Format(A)
200 Format(L1)
300 Format('=== Element Block             ', I3)
400 Format('=== Node set                  ', I3)
500 Format('=== Grains')
600 Format('Grain ', I3,' : theta = ', F6.2) 
  End Subroutine Ask_Rupt_Params2D

  Subroutine Ask_Rupt_Params3D(Geom, Params, BC_NS)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params3D)                          :: Params
    Type (Vect3D), Dimension(:), Pointer          :: BC_NS
    
    Character                                     :: Char_Input
    Integer                                       :: Num_TS
    Real(Kind = Kr)                               :: TS_init, TS_Final
    Real(Kind = Kr)                               :: Tmp_DistFactor
    Integer                                       :: i, iBlk, iSet
    Integer                                       :: Blk_ID, Set_ID

    Real(Kind = Kr)                               :: Tmp_E, Tmp_nu

    !!! The following are decent values.
    Params%nbCracks           = 0
    Params%MaxCrackLength     = 0.0             
    Params%MaxIterRelax       = 5000
    Params%TolKSP             = 1.0D-6
    Params%TolRelax           = 1.0D-3
    Params%Do_Irrev           = .TRUE.

    !!! This is totally meaningful in general...
    Params%Epsilon      = 1.0D-1
    Params%KEpsilon     = 1.0D-6
    Params%TolIrrev     = 1.0D-2
    Write(*,100) '===== Global Properties ======'

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
       Params%Load = (/ (TS_Init + (i-1.0)*(TS_Final - TS_Init)  / (Num_TS - 1.0_Kr), i = 1,Num_TS) /)
    End If
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
    Allocate (Params%Hookes_Law(Geom%Num_Elem_Blks))
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
       Write(*,100, advance = 'no') 'Toughness:                '
       Read(*,*) Params%Toughness(iBlk)
       Write(*,100, advance = 'no') 'Young''s modulus           '
       Read(*,*) Tmp_E
       Write(*,100, advance = 'no') 'Poisson ratio             ' 
       Read(*,*) Tmp_nu
       Write(*,100, advance = 'no') 'Thermal Expansion coef.   ' 
       Read(*,*) Params%Therm_Exp(iBlk)
       Call GenHL_Iso3D_Enu(Tmp_E, Tmp_nu, Params%Hookes_Law(iBlk))
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
       
       Write(*,100,advance = 'no') 'BC U type, Z direction (None=0 DIRI = 1) '
       Read(*,*) Params%BC_Type_Z(iSet) 
       If (Params%BC_Type_Z(iSet) == 1) Then
          Write(*,100, advance = 'no') 'BC U, Z direction                    '
          Read(*,*) BC_NS(iSet)%Z 
       End If
    End Do
    
100 Format(A,T78,': ')
200 Format(L1)
300 Format('=== Element Block             ', I3)
400 Format('=== Node set                  ', I3)
  End Subroutine Ask_Rupt_Params3D
End Program Prep_Rupt
