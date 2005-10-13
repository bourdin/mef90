#if defined PB_2D
Module m_Rupt2D_Proc
#elif defined PB_3D
Module m_Rupt3D_Proc
#else 
Module m_Rupt2DA_Proc
#endif

  Use m_MEF90
  Use m_Rupt_Struct
#if defined PB_2D
  Use m_Rupt2D_Vars	
  Use m_Rupt2D_U	
  Use m_Rupt2D_V
#elif defined PB_3D
  Use m_Rupt3D_Vars	
  Use m_Rupt3D_U	
  Use m_Rupt3D_V
#else 
  Use m_Rupt2DA_Vars	
  Use m_Rupt2DA_U	
  Use m_Rupt2DA_V
#endif


  Implicit NONE
  PRIVATE


#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscsys.h"

  Public :: Init
  Public :: Export
  Public :: Export_V
  Public :: Finalize
  Public :: Solve_V
  Public :: Solve_U
  Public :: Init_BC_U
  Public :: Update_BC_V
  Public :: Update_BC_U
  Public :: Update_F
  Public :: Update_Temp

!  integer :: i

Contains
  Subroutine Init()
    PetscTruth                              :: Has_Sim_Str
    Real(Kind = Kr), Dimension(:), Pointer  :: Tmp_Ptr, U_Ptr, V_Ptr
    
    Call MEF90_Initialize()
    MEF90_GaussOrder = 2 
    
    
    Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-f', Params%Sim_Str,    &
         & Has_Sim_Str, iErr)

!!! Check for the -restart flag
!!! If found, we are restarting a computation at step Timestep
!!! -restart 1 is equivalent to making a new computation
    TimeStep = 1
    Call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-restart', TimeStep,       &
         & Is_Restarting, PETSC_NULL_INTEGER, iErr)

    If (.NOT. Has_Sim_Str) Then
       Write(CharBuffer, 100) 'Simulation name: \n'c
       Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
       If (MEF90_MyRank ==0) Then
          Read(*,100) Params%Sim_Str
       End If
       Call MPI_BCAST(Params%Sim_Str, MXSTLN, MPI_CHARACTER, 0, MPI_COMM_WORLD,  &
            & iErr)
       Is_Interactive = .TRUE.
       Log_Unit = 6 !!! messages sent to stdout
    Else
       Is_Interactive = .FALSE.
       Log_Str          = Trim(Params%Sim_Str) // '.log'
       If (Is_Restarting == PETSC_TRUE) Then
          Open (File = Log_Str, Unit = Log_Unit, Status = 'Old',              &
               & Position = 'Append')
       Else
          Open (File = Log_Str, Unit = Log_Unit, Status = 'Replace')
       End If
    End If
    
    Geom%Filename    = Trim(Params%Sim_Str) // '.gen'
    Params%PARAM_Str = Trim(Params%Sim_Str) // '.PARAM'
    Params%CST_Str   = Trim(Params%Sim_Str) // '.CST'
    Ener_Str         = Trim(Params%Sim_Str) // '.ener'
    
    If ((MEF90_MyRank == 0) .AND. (Is_Restarting == PETSC_FALSE) )Then
       Open (File = Ener_Str, Unit = Ener_Unit, Status = 'Replace')
       Rewind(Ener_Unit)
       Close(Ener_Unit)
       Else
       Open (File = Ener_Str, Unit = Ener_Unit, Position = 'Append')
       Write(Ener_Unit, *)
       Write(Ener_Unit, *)
       Close(Ener_Unit)
    End If
    
    Call Read_EXO_Geom_Info(Geom)
    
#ifdef PB_2DA
    Call Read_EXO_Node_Coord(Geom, Node_db_U, 1)
#else
    Call Read_EXO_Node_Coord(Geom, Node_db_U, Geom%Num_Dim)
#endif
    Call Read_EXO_Node_Coord(Geom, Node_db_V, 1)
    Call Read_EXO_Connect(Geom, Elem_db_U) 

    Call Read_EXO_Connect(Geom, Elem_db_V)
 
    Call Read_Rupt_EXO_Params(Geom, Params)
    Call Read_Rupt_DATA(Geom, Params)

    Call Init_BC_U(Geom, Params, Node_db_U)
    Node_db_V%BC = BC_Type_NONE

    If (MEF90_MyRank == 0) Then
       Write(Log_Unit,*) 'Number of nodes:      ', Geom%Num_Nodes
       Write(Log_Unit,*) 'Number of elements:  ', Geom%Num_Elems
    End	If

    Call Init_SD_NoOvlp(Geom, Elem_db_U, MySD_U, U_Dist, U_Loc, U_Master)
    Call VecDuplicate(U_Master, BCU_Master, iErr)
    Call VecDuplicate(U_Master, F_Master, iErr)
    Call VecDuplicate(U_Dist, BCU_Dist, iErr)
    Call VecDuplicate(U_Dist, RHS_U, iErr)
    Call VecDuplicate(U_Dist, F_Dist, iErr)
    Call VecGhostGetLocalForm(BCU_Dist, BCU_Loc, iErr)
    Call VecGhostGetLocalForm(F_Dist, F_Loc, iErr)

    Call Init_SD_NoOvlp(Geom, Elem_db_V, MySD_V, V_Dist, V_Loc, V_Master)
    Call VecDuplicate(V_Dist, V_Old, iErr)
    Call VecDuplicate(V_Dist, V_Change, iErr)
    Call VecDuplicate(V_Dist, RHS_V, iErr)
    
    Call VecDuplicate (V_Dist, Temp_Dist, iErr)
    Call VecDuplicate (V_Master, Temp_Master, iErr)

    Call VecGhostGetLocalForm(V_Dist, V_Loc, iErr)
    Call VecGhostGetLocalForm(Temp_Dist, Temp_Loc, iErr)

#ifdef PB_2DA
    Call MatCreateMPIAIJ(PETSC_COMM_WORLD, MySD_U%Num_Nodes, MySD_U%Num_Nodes,&
         & Geom%Num_Nodes, Geom%Num_Nodes, 24,                                &
         & PETSC_NULL_INTEGER, 24, PETSC_NULL_INTEGER, MR_U, iErr)
#else
    Call MatCreateMPIAIJ(PETSC_COMM_WORLD, MySD_U%Num_Nodes, MySD_U%Num_Nodes,&
         & Geom%Num_Nodes * Geom%Num_Dim, Geom%Num_Nodes * Geom%Num_Dim, 80,  &
         & PETSC_NULL_INTEGER, 80, PETSC_NULL_INTEGER, MR_U, iErr)
#endif
    Call MatSetOption(MR_U, MAT_SYMMETRIC, iErr)
    Call MatSetFromOptions(MR_U, iErr)

    Call MatCreateMPIAIJ(PETSC_COMM_WORLD, MySD_V%Num_Nodes, MySD_V%Num_Nodes,&
         & Geom%Num_Nodes, Geom%Num_Nodes, 24,                                &
         & PETSC_NULL_INTEGER, 24, PETSC_NULL_INTEGER, MR_V, iErr)
    Call MatSetOption(MR_V, MAT_SYMMETRIC, iErr)
    Call MatSetFromOptions(MR_V, iErr)       
     
     
    Allocate(Surf_Ener(0:Size(Params%Load)))
    Allocate(Tot_Ener(0:Size(Params%Load)))
    Allocate(Bulk_Ener(0:Size(Params%Load)))
    
    Surf_Ener(0) = 0.0_Kr
    Tot_Ener(0)  = 0.0_Kr
    Bulk_Ener(0) = 0.0_Kr      


    If (TimeStep == 1) Then
       
       Call VecSet(0.0_Kr, U_Dist, iErr)
       Call VecSet(0.0_Kr, U_Loc, iErr)
       
       Call VecSet(1.0_Kr, V_Dist, iErr)
       Call VecSet(1.0_Kr, V_Loc, iErr)
       
       Call VecSet(0.0_Kr, V_Old, iErr)
    Else
!!! Read and update V from TimeStep - 1
       If (MEF90_MyRank == 0) Then
          Call VecGetArrayF90(V_Master, V_Ptr, iErr)
          Call Read_EXO_Result_Nodes(Geom, 1, TimeStep-1, Tmp_Ptr, 1)
          V_Ptr = Tmp_Ptr
          Call VecRestoreArrayF90(V_Master, V_Ptr, iErr)
          DeAllocate(Tmp_Ptr)
       EndIf
       
!!! V_Master -> V_Dist
       Call VecScatterBegin(V_Master, V_Dist, INSERT_VALUES,                &
            & SCATTER_REVERSE, MySD_V%ToMaster, iErr)
       Call VecScatterEnd(V_Master, V_Dist, INSERT_VALUES, SCATTER_REVERSE, &
            & MySD_V%ToMaster, iErr)
       
!!! V_Dist -> V_Loc
       Call VecGhostUpdateBegin(V_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
       Call VecGhostUpdateEnd(V_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)

       Call VecCopy(V_Dist, V_Old, iErr)
       
!!! Read and update U from TimeStep - 1
       If (MEF90_MyRank == 0) Then
          Call VecGetArrayF90(U_Master, U_Ptr, iErr)
#ifdef PB_2DA
          Call Read_EXO_Result_Nodes(Geom, 4, TimeStep-1, Tmp_Ptr, 1)
#else
          Call Read_EXO_Result_Nodes(Geom, 2, TimeStep-1, Tmp_Ptr, Geom%Num_Dim)
#endif
          U_Ptr = Tmp_Ptr
          Call VecRestoreArrayF90(U_Master, U_Ptr, iErr)
          DeAllocate(Tmp_Ptr)
       EndIf
       
!!! U_Master -> U_Dist
       Call VecScatterBegin(U_Master, U_Dist, INSERT_VALUES,                &
            & SCATTER_REVERSE, MySD_U%ToMaster, iErr)
       Call VecScatterEnd(U_Master, U_Dist, INSERT_VALUES, SCATTER_REVERSE, &
            & MySD_U%ToMaster, iErr)
       
!!! U_Dist -> U_Loc
       Call VecGhostUpdateBegin(U_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
       Call VecGhostUpdateEnd(U_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)

!!! Reads the energy from the .gen file
!!! At this point, I let each CPU read the .gen file
!!! This is unefficient but is done only once 
       Do iTS = 1, TimeStep 
          Call Read_EXO_Result_Global(Geom, 1, iTS, Bulk_Ener(iTS))
          Call Read_EXO_Result_Global(Geom, 2, iTS, Surf_Ener(iTS))
          Call Read_EXO_Result_Global(Geom, 3, iTS, Tot_Ener(iTS))
       End Do
    End If
       
    Call Init_Ksps()
    
    Call Random_Seed
    
100 Format(A)
  End Subroutine Init
  
  Subroutine Init_BC_U(Geom, Params, Nodes_U)
    Type(EXO_Geom_Info), Intent(IN)                  :: Geom
    Type(Rupt_Params), Intent(IN)                    :: Params
#if defined PB_2D
    Type(Node2D), Dimension(:), Pointer              :: Nodes_U
#elif defined PB_3D
    Type(Node3D), Dimension(:), Pointer              :: Nodes_U
#else 
    Type(Node2D), Dimension(:), Pointer              :: Nodes_U
#endif

    Integer                                          :: iSet, iN
    
    If (Geom%Numbering /= Numbering_PerNodes) Then
       Print*, '[ERROR] Init_BC not implemented for this numbering scheme'
       STOP
    End If

    Nodes_U(:)%BC = BC_Type_NONE
    Do iSet = 1, Geom%Num_Node_Sets
       Do iN = 1, Geom%Node_Set(iSet)%Num_Nodes
#if defined PB_2DA
          Nodes_U(Geom%Node_Set(iSet)%Node_ID(iN))%BC =                       &
               & Params%BC_Type_Z(iSet)
#else          
          Nodes_U(Geom%Num_Dim * (Geom%Node_Set(iSet)%Node_ID(iN)-1)+1)%BC    &
               & = Params%BC_Type_X(iSet)
          Nodes_U(Geom%Num_Dim * (Geom%Node_Set(iSet)%Node_ID(iN)-1)+2)%BC    &
               & = Params%BC_Type_Y(iSet)
#ifdef PB_3D
          Nodes_U(Geom%Num_Dim * Geom%Node_Set(iSet)%Node_ID(iN) )%BC         &
               & = Params%BC_Type_Z(iSet)
#endif
#endif
       End Do
    End Do
  End Subroutine Init_BC_U

  Subroutine Update_BC_V(Geoms, Params, SD, Nodes, V, VOld, TS)
!!! CHECK THAT OUT. I am tired...
    Type(EXO_Geom_Info), Intent(IN)                  :: Geoms
    Type(Rupt_Params), Intent(IN)                    :: Params
    Type (SD_Info), Intent(IN)                       :: SD
#ifdef PB_3D
    Type(Node3D), Dimension(:), Pointer              :: Nodes
#else
    Type(Node2D), Dimension(:), Pointer              :: Nodes
#endif
    Vec                                              :: V
    Vec                                              :: VOld    
    Integer                                          :: TS

    PetscReal, Dimension(:), Pointer                 :: VPtr, VPtr_Old
    PetscReal, Dimension(:), Pointer                 :: BC_V_Ptr, Tmp_Ptr    
    Integer                                          :: i, iSloc, iS
    Integer, Dimension(:), Pointer                   :: Loc_Indices
      
    If(TS ==0 ) Then
       Do iSLoc = 1, SD%Num_Nodes
          iS = SD%Node(iSLoc)
          Nodes(iS)%BC = BC_Type_NONE
       End Do
       Do iSLoc = 1, SD%Num_GhostNodes
          iS = SD%GhostNode(iSLoc)
          Nodes(iS)%BC = BC_Type_NONE
       End Do             
!       Do iS = 1, Geom%Num_Nodes
!          If ( (ABS(Nodes(iS)%Coord%X) <= 1.0e-5) .AND.                        & 
!          & (Nodes(iS)%Coord%Y <= 0.0_Kr) ) Then
!             Nodes(iS)%BC = BC_TYPE_Diri
!          End If
!       End Do
    Else
       ! Read the V 
       If (MEF90_MyRank == 0) Then       
          Call VecGetArrayF90(V_Master, BC_V_Ptr, iErr)
          Call Read_EXO_Result_Nodes(Geom, 1, TS, Tmp_Ptr, 1)
          BC_V_Ptr = Tmp_Ptr
          Call VecRestoreArrayF90(V_Master, BC_V_Ptr, iErr)
          DeAllocate(Tmp_Ptr)
       EndIf   
       Call VecScatterBegin(V_Master, VOld, INSERT_VALUES,SCATTER_REVERSE, MySD_V%ToMaster, iErr)
       Call VecScatterEnd(V_Master, VOld, INSERT_VALUES, SCATTER_REVERSE,  MySD_V%ToMaster, iErr)
       Call VecGhostUpdateBegin(VOld, INSERT_VALUES, SCATTER_FORWARD, iErr)
       Call VecGhostUpdateEnd(VOld, INSERT_VALUES, SCATTER_FORWARD, iErr)
       
       ! V -> VOld
       Call VecGetArrayF90(VOld, VPtr_Old, iErr)
       
       Allocate(Loc_Indices(Geoms%Num_Nodes))
       Loc_Indices = (/ (i ,i = 0, Geoms%Num_Nodes - 1) /)
       Call AOApplicationToPETSc(SD%Loc_AO, Geoms%Num_Nodes, Loc_Indices, iErr)   
       Call VecGetArrayF90(V, VPtr, iErr)        
       Do iSLoc = 1, SD%Num_Nodes
          iS = SD%Node(iSLoc)
          If ( VPtr_Old(Loc_Indices(iS)+1) <= Params%TolIrrev ) Then 
             Nodes(iS)%BC = BC_Type_DIRI
             VPtr(Loc_Indices(iS)+1) = 0.0_Kr
          End If 
       End Do
       Do iSLoc = 1, SD%Num_GhostNodes
          iS = SD%GhostNode(iSLoc)
         If (. VPtr_Old(Loc_Indices(iS)+1) <= Params%TolIrrev) Then 
             Nodes(iS)%BC = BC_Type_DIRI
             VPtr(Loc_Indices(iS)+1) = 0.0_Kr
          End If 
       End Do
       Call VecRestoreArrayF90(V, VPtr, iErr)
       Call VecRestoreArrayF90(VOld, VPtr_Old, iErr)
       Call VecScatterBegin(V, V_Master, INSERT_VALUES, SCATTER_FORWARD, MySD_V%ToMaster, iErr)
       Call VecScatterEnd(V, V_Master, INSERT_VALUES, SCATTER_FORWARD, MySD_V%ToMaster, iErr)    
       DeAllocate(Loc_Indices)       
    End If      
    
  End Subroutine Update_BC_V

  Subroutine Init_KSPs()
    Call KSPCreate(PETSC_COMM_WORLD, KSP_U, iErr)
    Call KSPSetOperators(KSP_U, MR_U, MR_U, SAME_NONZERO_PATTERN, iErr)
    Call KSPGetPC(KSP_U, PC_U, iErr)
    Call PCSetType(PC_U, PCBJACOBI, iErr)
    Call KSPSetType(KSP_U, KSPCG, iErr)

    Call KSPSetInitialGuessNonzero(KSP_U, PETSC_TRUE, iErr)

!!$    Call KSPSetTolerances(KSP_U, Params%TolKSP,                              &
!!$           & PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION,  &
!!$           & 50000, iErr)
    Call KSPSetTolerances(KSP_U, Params%TolKSP,                              &
           & PETSC_DEFAULT_DOUBLE_PRECISION, 1.0e+9, 50000, iErr)
    Call KSPSetFromOptions(KSP_U, iErr)

    Call KSPCreate(PETSC_COMM_WORLD, KSP_V, iErr)
    Call KSPSetOperators(KSP_V, MR_V, MR_V, SAME_NONZERO_PATTERN, iErr)
    Call KSPGetPC(KSP_V, PC_V, iErr)
    Call PCSetType(PC_V, PCBJACOBI, iErr)
    Call KSPSetType(KSP_V, KSPCG, iErr)

    Call KSPSetInitialGuessNonzero(KSP_V, PETSC_TRUE, iErr)

!!$    Call KSPSetTolerances(KSP_V, Params%TolKSP,                              &
!!$           & PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION,  &
!!$           & PETSC_DEFAULT_INTEGER, iErr)
    Call KSPSetTolerances(KSP_V, Params%TolKSP,                              &
           & PETSC_DEFAULT_DOUBLE_PRECISION, 1.0e+9, PETSC_DEFAULT_INTEGER,  &
           & iErr)
    Call KSPSetFromOptions(KSP_V, iErr)
  End Subroutine Init_KSPs

  Subroutine Update_BC_U(TimeStep)
    Integer                                       :: TimeStep

    PetscReal, Dimension(:), Pointer              :: BC_Ptr, Tmp_Ptr

    If (MEF90_MyRank == 0) Then
       Call VecGetArrayF90(BCU_Master, BC_Ptr, iErr)
#ifdef PB_2DA
       Call Read_EXO_Result_Nodes(Geom, 4, TimeStep, Tmp_Ptr, 1)
#else
       Call Read_EXO_Result_Nodes(Geom, 2, TimeStep, Tmp_Ptr, Geom%Num_Dim)
#endif
       BC_Ptr = Tmp_Ptr
       Call VecRestoreArrayF90(BCU_Master, BC_Ptr, iErr)
       DeAllocate(Tmp_Ptr)
    EndIf

     !!! BC_Master -> BC_Dist
     Call VecScatterBegin(BCU_Master, BCU_Dist, INSERT_VALUES,                &
          & SCATTER_REVERSE, MySD_U%ToMaster, iErr)
     Call VecScatterEnd(BCU_Master, BCU_Dist, INSERT_VALUES, SCATTER_REVERSE, &
          & MySD_U%ToMaster, iErr)

     !!! BC_Dist -> BC_Loc
     Call VecGhostUpdateBegin(BCU_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
     Call VecGhostUpdateEnd(BCU_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
   End Subroutine Update_BC_U

  Subroutine Update_F(TimeStep)
    Integer                                       :: TimeStep

    PetscReal, Dimension(:), Pointer              :: F_Ptr, Tmp_Ptr

    If (MEF90_MyRank == 0) Then
       Call VecGetArrayF90(F_Master, F_Ptr, iErr)
#ifdef PB_2DA
       Call Read_EXO_Result_Nodes(Geom, 7, TimeStep, Tmp_Ptr, 1)
#else
       Call Read_EXO_Result_Nodes(Geom, 5, TimeStep, Tmp_Ptr, Geom%Num_Dim)
#endif
       F_Ptr = Tmp_Ptr
       Call VecRestoreArrayF90(F_Master, F_Ptr, iErr)
       DeAllocate(Tmp_Ptr)
    EndIf

     !!! F_Master -> F_Dist
     Call VecScatterBegin(F_Master, F_Dist, INSERT_VALUES,                &
          & SCATTER_REVERSE, MySD_U%ToMaster, iErr)
     Call VecScatterEnd(F_Master, F_Dist, INSERT_VALUES, SCATTER_REVERSE, &
          & MySD_U%ToMaster, iErr)

     !!! F_Dist -> F_Loc
     Call VecGhostUpdateBegin(F_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
     Call VecGhostUpdateEnd(F_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
   End Subroutine Update_F

  Subroutine Update_Temp(TimeStep)
  !!! Careful with variable names... Temp == temperature, tmp == temporary
    Integer                                       :: TimeStep

    PetscReal, Dimension(:), Pointer              :: Temp_Ptr, Tmp_Ptr

    If (MEF90_MyRank == 0) Then
       Call VecGetArrayF90(Temp_Master, Temp_Ptr, iErr)
       Call Read_EXO_Result_Nodes(Geom, 8, TimeStep, Tmp_Ptr, 1)
       Temp_Ptr = Tmp_Ptr
       Call VecRestoreArrayF90(Temp_Master, Temp_Ptr, iErr)
       DeAllocate(Tmp_Ptr)
    EndIf

     !!! Temp_Master -> Temp_Dist
     Call VecScatterBegin(Temp_Master, Temp_Dist, INSERT_VALUES,              &
          & SCATTER_REVERSE, MySD_V%ToMaster, iErr)
     Call VecScatterEnd(Temp_Master, Temp_Dist, INSERT_VALUES,                & 
          & SCATTER_REVERSE, MySD_V%ToMaster, iErr)

     !!! Temp_Dist -> Temp_Loc
     Call VecGhostUpdateBegin(Temp_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
     Call VecGhostUpdateEnd(Temp_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
   End Subroutine Update_Temp


  Subroutine Solve_U()
!    Call Update_BC_U(TimeStep)

    Call Assemb_MR_U(MR_U, V_Loc, Geom, Params, MySD_U, MySD_V,               &
         & Elem_db_U, Elem_db_V, Node_db_U, Node_db_V )
    Call Assemb_RHS_U(RHS_U, BCU_loc, Geom, Params, MySD_U, Elem_db_U,        &
         & Node_db_U, MySD_V, Elem_db_V, Node_db_V, V_Loc, F_Loc, Temp_Loc)
    Call PetscGetTime(SolveTS, iErr)
    Call KSPSolve(KSP_U, RHS_U, U_Dist, iErr)

#ifdef MEF90_TIMING
    Call PetscGetTime(SolveTF, iErr)
    If (MEF90_MyRank == 0) Then
       Write(CharBuffer,*) 'Total time in KSP_Solve:                 ',      &
            & SolveTF- SolveTS, '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
    End If
#endif

    Call VecGhostUpdateBegin(U_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
    Call VecGhostUpdateEnd(U_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
  End Subroutine Solve_U

  Subroutine Solve_V()
    
!!$    If (Params%Do_Irrev) Then
!!$       Call Apply_BC_V(Geom, Params, MySD_V, Node_db_V, V_Dist)
!!$       Call VecGhostUpdateBegin(V_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
!!$       Call VecGhostUpdateEnd(V_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
!!$    End If

    Call Assemb_MR_V(MR_V, U_Loc, Temp_Loc, Geom, Params, MySD_U, MySD_V,     &
         & Elem_db_U, Elem_db_V, Node_db_U, Node_db_V )
!    Call Assemb_RHS_V(RHS_V, Geom, Params, MySD_V, Elem_db_V, Node_db_V)

    Call PetscGetTime(SolveTS, iErr)
    Call KSPSolve(KSP_V, RHS_V, V_Dist, iErr)

#ifdef MEF90_TIMING
    Call PetscGetTime(SolveTF, iErr)
    If (MEF90_MyRank == 0) Then
       Write(CharBuffer,*) 'Total time in KSP_Solve:                 ',      &
            & SolveTF- SolveTS, '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
    End If
#endif

    Call VecGhostUpdateBegin(V_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
    Call VecGhostUpdateEnd(V_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
  End Subroutine Solve_V
  


  Subroutine Export(TS)
    Integer, Intent(IN)                           :: TS


    Real(Kind = Kr), Dimension(:), Pointer        :: SOL_Ptr

    Call VecScatterBegin(U_Dist, U_Master, INSERT_VALUES, SCATTER_FORWARD,    &
         & MySD_U%ToMaster, iErr)
    Call VecScatterEnd(U_Dist, U_Master, INSERT_VALUES, SCATTER_FORWARD,      &
         & MySD_U%ToMaster, iErr)
    
    If (MEF90_MyRank == 0) Then
       Call VecGetArrayF90(U_Master, Sol_Ptr, iErr)
#ifdef PB_2DA
       Call Write_EXO_Result_Ptr_Nodes(Geom, 4, TS, SOL_Ptr)
#else
       Call Write_EXO_Result_Ptr_Nodes(Geom, 2, TS, SOL_Ptr)
#endif
       Call VecRestoreArrayF90(U_Master, SOL_Ptr, iErr)
    End If
    Call VecScatterBegin(V_Dist, V_Master, INSERT_VALUES, SCATTER_FORWARD,    &
         & MySD_V%ToMaster, iErr)
    Call VecScatterEnd(V_Dist, V_Master, INSERT_VALUES, SCATTER_FORWARD,      &
         & MySD_V%ToMaster, iErr)
    
    If (MEF90_MyRank == 0) Then
       Call VecGetArrayF90(V_Master, Sol_Ptr, iErr)
       Call Write_EXO_Result_Ptr_Nodes(Geom, 1, TS, SOL_Ptr)
       Call VecRestoreArrayF90(V_Master, SOL_Ptr, iErr)
    
       Call Write_EXO_Result_Global(Geom, 1, TS, Bulk_Ener(TS))
       Call Write_EXO_Result_Global(Geom, 2, TS, Surf_Ener(TS))
       Call Write_EXO_Result_Global(Geom, 3, TS, Tot_Ener(TS))
    End If
  End Subroutine Export

  Subroutine Export_V(TS)
    Integer, Intent(IN)                           :: TS


    Real(Kind = Kr), Dimension(:), Pointer        :: SOL_Ptr

    Call VecScatterBegin(V_Dist, V_Master, INSERT_VALUES, SCATTER_FORWARD,    &
         & MySD_V%ToMaster, iErr)
    Call VecScatterEnd(V_Dist, V_Master, INSERT_VALUES, SCATTER_FORWARD,      &
         & MySD_V%ToMaster, iErr)
    
    If (MEF90_MyRank == 0) Then
       Call VecGetArrayF90(V_Master, Sol_Ptr, iErr)
       Call Write_EXO_Result_Ptr_Nodes(Geom, 1, TS, SOL_Ptr)
       Call VecRestoreArrayF90(V_Master, SOL_Ptr, iErr)
    End If
  End Subroutine Export_V


  Subroutine Finalize()
    DeAllocate(Node_db_U, Elem_db_U)
    DeAllocate(Node_db_V, Elem_db_V)
    
    Call MatDestroy(MR_U, iErr)
    Call MatDestroy(MR_V, iErr)

    Call VecDestroy(RHS_U, iErr)
    Call VecDestroy(U_Dist, iErr)
    Call VecDestroy(U_Loc, iErr)
    Call VecDestroy(U_Master, iErr)

    Call VecDestroy(RHS_V, iErr)
    Call VecDestroy(V_Dist, iErr)
    Call VecDestroy(V_Loc, iErr)
    Call VecDestroy(V_Master, iErr)


    Call VecDestroy(BCU_Dist, iErr)
    Call VecDestroy(BCU_Loc, iErr)
    Call VecDestroy(BCU_Master, iErr)

    Call VecDestroy(F_Dist, iErr)
    Call VecDestroy(F_Loc, iErr)
    Call VecDestroy(F_Master, iErr)

    Call VecDestroy(Temp_Dist, iErr)
    Call VecDestroy(Temp_Loc, iErr)
    Call VecDestroy(Temp_Master, iErr)
! Call PetscFinalize(PETSC_NULL_CHARACTER, iErr)
!Call MPI_FINALIZE(iErr)    
!    Call MEF90_Finalize()
  End Subroutine Finalize

#if defined PB_2D
End Module m_Rupt2D_Proc
#elif defined PB_3D
End Module m_Rupt3D_Proc
#else 
End Module m_Rupt2DA_Proc
#endif
