#ifdef PB_2D
Module m_OvSch2D_Procs
!  Use mpi
  Use m_MEF90
  Use m_Poisson_Struct
  Use m_OvSch2D_Vars
#else
Module m_OvSch3D_Procs
!  Use mpi
  Use m_MEF90
  Use m_Poisson_Struct
  Use m_OvSch3D_Vars
#endif


  Implicit NONE
  PRIVATE
#include "include/finclude/petsc.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petscao.h"

  Public :: Init_Broadcast
  Public :: Init_Node_Ownership
  Public :: Init_Size
  Public :: Init_Size_Local
  Public :: Assemb_RHS_Poisson_Local
  Public :: Init_PU
  Public :: G_PU_2D_Scal_Dist
  Public :: Communicate_PU
  Public :: Assemb_Coarse_Seq
  Public :: Coarse_Solve_Prep
  Public :: Coarse_Solve_Seq
  Public :: SET_OBDD
  Public :: SET_ASM
  Public :: Init_KSP_OLD
  Public :: Init_KSP
  Public :: Finalize
  Public :: CompL2
  Public :: CompSEMIH2
  Public :: Export
  Public :: SaveLayout
  Public :: Assemb_Mat_Poisson
  Public :: Assemb_Mat_Local_Poisson
  Public :: Assemb_Mat_Local_Poisson_Neu
  Public :: Assemb_RHS_Poisson

Contains


!!! Initialize Broadcast

  Subroutine Init_Broadcast()
    Integer                            :: iE, iN, Nb_DoF
    
    Integer, Dimension(:), Pointer                    :: METIS_elmnts
    Integer                                           :: METIS_etype
    Integer                                           :: METIS_numflag=1
    Integer                                           :: METIS_edgecut
    Integer, Dimension(:), Pointer                    :: METIS_epart
    Integer, Dimension(:), Pointer                    :: METIS_npart
 
    PetscLogDouble                                    :: metisTS, metisTF
    PetscTruth  :: flg

    Call MEF90_Initialize()
    Call MPI_COMM_RANK(MPI_COMM_WORLD, MyRank, iErr)
    Call MPI_COMM_SIZE(MPI_COMM_WORLD, NumProcs, iErr)
!!Read    If (MyRank ==0) Then
    Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-filename',             &
         & Params%Sim_Str,flg,ierr)
    
    If (MyRank ==0) Then
       Write(CharBuffer, *) 'Simulation name:   ', Params%Sim_Str, '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
    End If
    Geom%Filename    = Trim(Params%Sim_Str) // '.gen'
    Params%PARAM_Str = Trim(Params%Sim_Str) // '.PARAM'
    
    Call Read_EXO_Geom_Info(Geom)
    Call Read_EXO_Node_Coord(Geom, Node_db, 1)
    Call Read_EXO_Connect(Geom, Elem_db) 
    Call Read_Poisson_EXO_Params(Geom, Params)
    Call Read_Poisson_DATA(Geom, Params)
    Call Read_EXO_Result_Nodes(Geom, 1, 1, Load)
    
    If (MyRank == 0) Then
       Write(CharBuffer,*) 'Number of nodes:         ', Geom%Num_Nodes,       &
            & '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
       Write(CharBuffer,*) 'Number of elements:      ', Geom%Num_Elems, '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
       !!Change       Call Read_EXO_Result_Nodes(Geom, 1, 1, Load)
       Write(CharBuffer,*) 'Done reading... \n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
       !!Read    End If
       !!Change    If (MyRank == 0) Then
       Write(CharBuffer, *) 'OK Read ALL\n'c
    End If
    
    Call Init_BC(Geom, Params, Node_db)
    If (MyRank == 0) Then
       Write(CharBuffer, *) 'OK Init_BC\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr) 
    End If
    
    Allocate(Elem_Owner(Geom%Num_Elems))
    Allocate(Node_Owner(Geom%Num_Nodes))
    If (MyRank == 0) Then
       Write(CharBuffer, *) 'Calling METIS\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
       Call PetscGetTime(metisTS, iErr)
    End If
#ifdef PB_2D
    METIS_etype = 1
#else
    METIS_etype = 2
#endif
    
    Nb_DoF = Elem_db(1)%Nb_DoF
    Allocate (METIS_elmnts(Geom%Num_Elems*Nb_DoF))
    Do iE = 1, Geom%Num_Elems
       METIS_elmnts((iE-1)*Nb_DoF +1:iE*Nb_DoF) = Elem_db(iE)%ID_Dof        
    End Do
    
    Allocate(METIS_epart(Geom%Num_Elems))
    Allocate(METIS_npart(Geom%Num_Nodes))
    
    Call METIS_PartMeshNodal(Geom%Num_Elems, Geom%Num_Nodes,                  &
         & METIS_elmnts, METIS_etype, METIS_numflag, NumProcs,                &
         & METIS_edgecut, METIS_epart, METIS_npart)
    Elem_Owner = METIS_epart-1
    Node_Owner = METIS_npart-1
    DeAllocate(METIS_epart)
    DeAllocate(METIS_npart)
    If (MyRank == 0) Then
       Call PetscGetTime(metisTF, iErr)
       Write(CharBuffer,*) 'Total time in METIS                      ',       &
            & metisTF - metisTS, '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
       Write(CharBuffer,*) 'Total edge cuts                    ',             &
            & METIS_edgecut, '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
    End If
    
    Call MPI_BCAST(Elem_Owner, Geom%Num_Elems, MPI_INTEGER, 0,                &
         & PETSC_COMM_WORLD, iErr)
    Call MPI_BCAST(Node_Owner, Geom%Num_Nodes, MPI_INTEGER, 0,                &
         & PETSC_COMM_WORLD, iErr)
100 Format(A)
  End Subroutine Init_Broadcast
  
!!! Node Ownership
  Subroutine Init_Node_Ownership()
    Integer                            :: i, j, iE, ct
    
    MySize_Elem = count(Elem_Owner==MyRank)
    MySize_Node = count(Node_Owner==MyRank)
    Allocate(PU_node_size(NumProcs))
    Allocate(PU_node_init(NumProcs))
    ct=0
    Do i = 1, NumProcs
       PU_node_init(i)=ct
       PU_node_size(i)=count(Node_Owner==i-1)
       ct=ct+PU_node_size(i)
    End Do
    Call VecCreateMPI(PETSC_COMM_WORLD, MySize_Node, Geom%Num_Nodes, RHS,    &
         & ierr)
    Call VecGetOwnershipRange(RHS, MyIdxMin_Node, MyIdxMax_Node,ierr)
    Allocate(Renum(MySize_Node))
    Renum   = (/ (i, i = MyIdxMin_Node, MyIdxMax_Node-1 )/)
    If(MyRank ==0) Then       
       Write(CharBuffer,*) 'Data layout:\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
    End If
    Write(CharBuffer,300) MyRank, MySize_Elem, MySize_Node
    Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
    Call PetscSynchronizedFlush(PETSC_COMM_WORLD, CharBuffer, iErr)
    Allocate(MyElem(MySize_Elem))
    Allocate(MyNode(MySize_Node))    
    iE =0
    Do i=1, Geom%Num_Elems
       If (Elem_Owner(i)==MyRank) Then
          iE = iE+1
          MyElem(iE) = i-1
       End If
    End Do
    L_c =0
    Do i=1, Geom%Num_Nodes
       If ( Node_Owner(i)== MyRank) Then
          L_c = L_c +1
          MyNode(L_c) = i-1  
       End If
    End Do
    Call ISCreateGeneral(PETSC_COMM_WORLD, MySize_Node, MyNode, Mangled_IS, iErr)
    Call ISCreateGeneral(PETSC_COMM_WORLD, MySize_Node, Renum, Layout_IS, iErr)
    Call AOCreateBasicIS(Mangled_IS, Layout_IS, DeMangle_AO, iErr)
    Allocate(All_list(0:Geom%Num_Nodes-1))
    All_list = (/ (i, i=0, Geom%Num_Nodes-1) /)
    Allocate(PU_node_list(0:Geom%Num_Nodes-1))
    PU_node_list = (/ (i, i=0, Geom%Num_Nodes-1) /)
    Call AOApplicationToPetsc(DeMangle_AO, Geom%Num_Nodes, All_list, iErr)
    Call AOPetscToApplication(DeMangle_AO, Geom%Num_Nodes, PU_node_list, iErr)
    !    Call AOView(DeMangle_AO, PETSC_VIEWER_STDOUT_WORLD, iErr)
300 Format('Rank: ', I4,' # of Elements  : ', I7,' # of Nodes : ', I7,    &
         & ' \n'c)
  End Subroutine Init_Node_Ownership
  
  Subroutine Init_Size()
    
    Call VecDuplicate(RHS, SOL, iErr)
    
    !! For Export
    If (MyRank == 0) Then
       Call VecCreateMPI(PETSC_COMM_WORLD, Geom%Num_Nodes,     &
            & Geom%Num_Nodes, SOL_Master, iErr)
    Else
       Call VecCreateMPI(PETSC_COMM_WORLD, 0, Geom%Num_Nodes,  &
            & SOL_Master, iErr)
    End If
    Call VecScatterCreate(SOL, Layout_IS, SOL_Master, Mangled_IS,       &
         & SOL_DistToMaster, iErr)
    
    Call MatCreateMPIAIJ(PETSC_COMM_WORLD, MySize_Node, MySize_Node,    &
         & Geom%Num_Nodes, Geom%Num_Nodes, 25, PETSC_NUll_INTEGER, 25,    &
         & PETSC_NULL_INTEGER, MR, iErr)
    !    Call MatSetOption(MR, MAT_SYMMETRIC, iErr)
    Call MatSetFromOptions(MR, iErr)
    Allocate (Load(Geom%Num_Nodes))
    If (MyRank == 0) Then
       Call Read_EXO_Result_Nodes(Geom, 1, 1, Load)
    End If
    Call MPI_BCAST(Load, Geom%Num_Nodes, MPI_DOUBLE_PRECISION, 0,    &
         & PETSC_COMM_WORLD, iErr)
  End Subroutine Init_Size
  
!!! Initialize the data for local structures
  Subroutine Init_Size_Local()
    
    Call MatCreateSeqAIJ(PETSC_COMM_SELF, MySize_Node_Ovlp, &
         & MySize_Node_Ovlp, 25, PETSC_NUll_INTEGER, MyMR, iErr)
    Call MatSetFromOptions(MyMR, iErr)
    Call MatCreateSeqAIJ(PETSC_COMM_SELF, MySize_Node_Ovlp, &
         & MySize_Node_Ovlp, 25, PETSC_NUll_INTEGER, MyMR_Neu, iErr)
    Call MatSetFromOptions(MyMR_Neu, iErr)
    Call MatCreateSeqAIJ(PETSC_COMM_SELF, MySize_Node_Ovlp, &
         & MySize_Node_Ovlp, 25, PETSC_NUll_INTEGER, MyMR_Re, iErr)
    Call MatSetFromOptions(MyMR_Re, iErr)
    Call VecCreateGhost(PETSC_COMM_WORLD,MySize_Node,Geom%Num_Nodes,          &
         &    MySize_Node_Ghost_Metis,All_list(MyNode_Ghost),RHS_Ov,ierr)
    Call VecCreateGhost(PETSC_COMM_WORLD,MySize_Node,Geom%Num_Nodes,          &
         &    MySize_Node_Ghost_Metis,All_list(MyNode_Ghost),SOL_Ov,ierr)
    Call VecCreateGhost(PETSC_COMM_WORLD,MySize_Node,Geom%Num_Nodes,          &
         &    MySize_Node_Ghost_Metis,All_list(MyNode_Ghost),RES_Ov,ierr)
    Call VecCreateGhost(PETSC_COMM_WORLD,MySize_Node,Geom%Num_Nodes,          &
         &    MySize_Node_Ghost_Metis,All_list(MyNode_Ghost),IMP_Ov,ierr)
    Call VecCreateGhost(PETSC_COMM_WORLD,MySize_Node,Geom%Num_Nodes,          &
         &    MySize_Node_Ghost_Metis,All_list(MyNode_Ghost),CoPU_Ov,ierr)
    Call VecCreateGhost(PETSC_COMM_WORLD,MySize_Node,Geom%Num_Nodes,          &
         &    MySize_Node_Ghost_Metis,All_list(MyNode_Ghost),Update_Ov,ierr)

   End Subroutine Init_Size_local


!!! Assemble Local RHS

  Subroutine  Assemb_RHS_Poisson_Local(RHS_Ptr)
    PetscScalar, Dimension(:), Pointer                  :: RHS_Ptr
    Integer                                             :: i 
    Real(Kind = Kr), Dimension(:), Pointer                        :: value

    Call VecGhostGetLocalForm(RHS_Ov,MyRHS_Ov,ierr) 
    Call VecGetOwnershipRange(RHS_Ov, MyIdxMin_Node_Ov, MyIdxMax_Node_Ov,ierr)

    Allocate(Renum2(MySize_Node))
    Renum2   = (/ (i, i = MyIdxMin_Node_Ov, MyIdxMax_Node_Ov-1)/)
    Allocate(value(MySize_Node))
    value=RHS_Ptr(Renum2-MyIdxMin_Node_Ov+1)
    Call VecSetValues(RHS_Ov,MySize_Node,Renum2,value,INSERT_VALUES,ierr)
    Call VecAssemblyBegin(RHS_Ov,ierr) 
    Call VecAssemblyEnd(RHS_Ov,ierr)
!    Call VecView(RHS_Ov,PETSC_VIEWER_STDOUT_WORLD,ierr)
    Call VecGhostUpdateBegin(RHS_Ov,INSERT_VALUES,SCATTER_FORWARD,ierr) 
    Call VecGhostUpdateEnd(RHS_Ov,INSERT_VALUES,SCATTER_FORWARD,ierr)
    Call VecGhostRestoreLocalForm(RHS_Ov,MyRHS_Ov,ierr)

  End Subroutine  Assemb_RHS_Poisson_Local


!!! Initialize PU stuff
  Subroutine Init_PU()
  
    PetscTruth  :: flg

   If (MyRank ==0) Then
     Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-ovlp',ovlp,flg,ierr)
     Write(CharBuffer, *) 'Overlapping Size:   ', ovlp, '\n'c
     Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
   End If
     Call MPI_BCAST(ovlp, 1, MPI_INTEGER, 0,                &
         & PETSC_COMM_WORLD, iErr)
    100 Format(A)
  End Subroutine Init_PU


!!! Generating PU

  Subroutine G_PU_2D_Scal_Dist(Geom, Elem_db, Node_db, ovlp, MyPU, &
       & MyElem_Ovlp,  MyNode_Ovlp, MyElem_Ghost,  MyNode_Ghost, MyNode_Bd)
    
    Type (EXO_Geom_Info), Intent(INOUT)            :: Geom
#ifdef PB_2D
    Type (Element2D_Scal), Dimension(:), Pointer   :: Elem_db 
    Type (Node2D), Dimension(:), Pointer           :: Node_db 
#else
    Type (Element3D_Scal), Dimension(:), Pointer   :: Elem_db 
    Type (Node3D), Dimension(:), Pointer           :: Node_db 
#endif
    Integer                                        :: ovlp
    Real(Kind = Kr), Dimension(:), Pointer         :: MyPU
    Integer, Dimension(:), Pointer                 :: MyElem_Ovlp
    Integer, Dimension(:), Pointer                 :: MyNode_Ovlp
    Integer, Dimension(:), Pointer                 :: MyElem_Ghost
    Integer, Dimension(:), Pointer                 :: MyNode_Ghost
    Integer, Dimension(:), Pointer                 :: MyNode_Bd
    
!!! Local variables. Pointers MUST be Deallocated before exiting
    Integer, Dimension(:), Pointer                 :: Total
    Integer                                        :: iE, iO
    Integer                                        :: iFx, iEy, iEz
    Integer                                        :: i, iNe, iGr
    Integer                                        :: IEx, iB, iDB,  iDBc
    Integer                                        :: iErr
    Integer                                        :: Nb_DoF
    
    Integer(Kind = Ki)                             :: iBlk
    Integer                                        :: j
    
    Integer, Dimension(:), Pointer                 :: IS_Nbr
    Integer, Dimension(:,:), Pointer               :: Mr_Nbr
    
    
    !!! Allocation of global variables
    Allocate (MyPU(Geom%Num_Nodes))
    Allocate (MyPU2(Geom%Num_Nodes))
    Allocate (MyCount_function(Geom%Num_Nodes))
    Allocate (MyCount_ID_Elem(Geom%Num_Elems))
    Allocate (MyCount_ID_Node(Geom%Num_Nodes))
    Allocate (MyCount_ID_Elem_Ghost(Geom%Num_Elems))
    Allocate (MyCount_ID_Node_Ghost(Geom%Num_Nodes))
    Allocate (MyCount_ID_Node_METIS(Geom%Num_Nodes))
    Allocate (MyCount_ID_Node_Old(Geom%Num_Nodes))
    Allocate (MyCount_ID_Node_Bd(Geom%Num_Nodes))
    Allocate (MyCount_ID_Node_Ghost_Metis(Geom%Num_Nodes))

    !!! Allocation of local variables
    Allocate (Total(Geom%Num_Nodes))

!! Initialize PU for nonoverlapping subdomains
    MyCount_function = 0
    MyCount_ID_Elem = 0
    MyCount_ID_Node = 0
    MyCount_ID_Node_METIS = 0
    MyCount_ID_Elem_Ghost = 0
    MyCount_ID_Node_Ghost = 0
    MyCount_ID_Node_Ghost_METIS =0
    Total = 0

    Do iBlk = 1, Geom%Num_Elem_Blks
       Nb_DoF = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem
       Do iE = 1, Geom%Num_Elems  
          If (Elem_Owner(iE) /= MyRank) Then
             CYCLE
          End If
          MyCount_ID_Node(Elem_db(iE)%ID_DoF) =1     
          MyCount_ID_Elem(iE) = 1
          Do iFx =  1, Nb_DoF
             If (Node_Owner(Elem_db(iE)%ID_DoF(iFx))== MyRank) Then
                MyCount_ID_Node_METIS(Elem_db(iE)%ID_DoF(iFx)) =1
             End If
          End Do
       End Do
    End Do
    MySize_Node_METIS=count(MyCount_ID_Node_METIS == 1)
    MyPU=MyCount_ID_Node
    If (ovlp >= 1) then
       Do_iBlk: Do iBlk = 1, Geom%Num_Elem_Blks
          Nb_DoF = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem 
          Do iO = 1, ovlp
             MyCount_ID_Node_Bd=0
             MyCount_ID_Node_Old=0
             MyCount_ID_Node_Old = MyCount_ID_Node
             Do iE = 1, Geom%Num_Elems
                If (any(MyCount_ID_Node(Elem_db(iE)%ID_DoF)> 0) .and. MyCount_ID_Elem(iE) == 0) Then
                   MyCount_ID_Elem(iE) = -1
                   MyCount_ID_Elem_Ghost(iE) = 1
                End If
             End Do
             !! Change the count function for the extension
             Do iE = 1, Geom%Num_Elems 
                If (MyCount_ID_Elem(iE) == -1) Then
                   Do iFx =  1, Nb_DoF
                      If (MyCount_ID_Node(Elem_db(iE)%ID_DoF(iFx)) == 0) Then
                         MyCount_ID_Node_Ghost(Elem_db(iE)%ID_DoF(iFx)) = 1
                         MyCount_ID_Node(Elem_db(iE)%ID_DoF(iFx)) = 1
                      End If
                      If (MyCount_ID_Node_METIS(Elem_db(iE)%ID_DoF(iFx)) == 0) Then
                         MyCount_ID_Node_Ghost_METIS(Elem_db(iE)%ID_DoF(iFx)) = 1
                         MyCount_ID_Node_METIS(Elem_db(iE)%ID_DoF(iFx)) = 1 
                      End If
                   End Do
                   MyCount_ID_Elem(iE) = 1
                End If
             End Do
             MyCount_ID_Node_Bd(:)=MyCount_ID_Node_METIS(:)-MyCount_ID_Node_Old(:)
          End Do
          MyPU(:) = MyCount_ID_Node(:) + MyPU(:)
       End Do Do_iBlk
    End If
    Call MPI_REDUCE(Int(MyPU),Total,Geom%Num_Nodes,MPI_Integer,MPI_SUM,0,     &
         & PETSC_COMM_WORLD, iErr)
    Call MPI_BCAST(Total, Geom%Num_Nodes,MPI_Integer, 0, PETSC_COMM_WORLD, iErr)

    MySize_Elem_Ovlp = count(MyCount_ID_Elem==1)
    MySize_Node_Ovlp = count(MyCount_ID_Node==1)
    MySize_Elem_Ghost = count(MyCount_ID_Elem_Ghost==1)
    MySize_Node_Ghost = count(MyCount_ID_Node_Ghost==1)
    MySize_Node_Ovlp_Metis = count(MyCount_ID_Node_Metis==1)
    MySize_Node_Ghost_Metis = count(MyCount_ID_Node_Ghost_Metis==1)
    MySize_Node_Bd = count(MyCount_ID_Node_Bd==1)
    Size_Node_Diri_Bd = count(Node_db%BC /= BC_Type_NONE)
    Allocate(Node_Diri_Bd(Size_Node_Diri_Bd))
    iDBc=0  
    Do iDB= 1, Geom%Num_Nodes
       If (Node_db(iDB)%BC /= BC_Type_NONE) Then
          iDBc=iDBc+1
          Node_Diri_Bd(iDBc) = iDB
       End If
    End Do
    !! WE will use  MyPU for coarse space and MyPU2 for weight 
    MyPU =MyPU/Total
    MyPU2=MyPU
    MyPU(Node_Diri_Bd)=0.0_Kr
    Allocate (MyElem_Ovlp(MySize_Elem_Ovlp))
    Allocate (MyElem_Ghost(MySize_Elem_Ghost))
    Allocate (MyNode_Ovlp(MySize_Node_Ovlp_Metis))
    Allocate (MyNode_Ghost(MySize_Node_Ghost_Metis))
    Allocate (MyNode_Bd(MySize_Node_Bd))
    Allocate (MyNode_Bd_Local_list(MySize_Node_Bd))
    Allocate (MyNode_Ovlp_list(0:Geom%Num_Nodes-1))
    MyNode_Ovlp_list = -1
    iEx =0
    iEy =0
    Do iE = 1, Geom%Num_Elems
       If (MyCount_ID_Elem(iE) == 1) Then
          iEx = iEx+1
          MyElem_Ovlp(iEx)=iE
       End If
       If (MyCount_ID_Elem_Ghost(iE) == 1) Then
          iEy = iEy+1
          MyElem_Ghost(iEy)=iE
       End If
    End Do
    iEx =0
    iEy =0
    iEz=0
    Do iE = 1, Geom%Num_Nodes
       If (Node_Owner(iE)== MyRank) Then
          iEx = iEx+1
          MyNode_Ovlp_list(iE-1)=iEx-1
       End If
       If (MyCount_ID_Node_Bd(iE) == 1) Then
          iEz = iEz+1
          MyNode_Bd(iEz)=iE-1
       End If
       If (MyCount_ID_Node_Ghost_Metis(iE) == 1) Then
          iEy = iEy+1
          MyNode_Ghost(iEy)=iE-1
          MyNode_Ovlp_list(iE-1)=iEy-1+MySize_Node
       End If
    End Do
    MyNode_Ovlp(1:MySize_Node) = MyNode
    MyNode_Ovlp(MySize_Node+1:MySize_Node_Ovlp_Metis) =MyNode_Ghost
    
    
    Allocate(IS_Nbr(0:NumProcs-1))
    Allocate(Mr_Nbr(0:NumProcs-1, 0:NumProcs-1))
    
    IS_Nbr = -1
    IS_Nbr(MyRank) = 1
    
    Do i=1, MySize_Node_Ghost_Metis
       IS_Nbr(Node_Owner(MyNode_Ghost(i)+1))=1
    End Do
    
    Call MPI_AllGather(Is_Nbr, NumProcs, MPI_INTEGER, Mr_Nbr, NumProcs,       &
         & MPI_INTEGER,PETSC_COMM_WORLD, iErr)
    
    
!!! Symmetrizing the Mr_Nbr matrix
    Do i = 0, NumProcs-1
       Do j = 0, i
          If (Mr_Nbr(i,j) /= Mr_Nbr(j,i)) Then
             Mr_Nbr(i,j) = 1
             Mr_Nbr(j,i) = 1
          End If
       End Do
    End Do
    
    
    
    MySize_Neighbor = count(Mr_Nbr(MyRank,:)==1)-1
    Allocate(MyNeighbor(MySize_Neighbor))
    Allocate(MyNeighbor_list(NumProcs))
    MySize_Group=MySize_Neighbor+1
    
    Allocate(MyGroup(MySize_Group))
    Allocate(MyGroup_list(NumProcs))
    MyNeighbor_list=-2
    MyGroup_list =-2
    iNe =0
    iGr =0
    Do i=0, NumProcs-1
       If (MyRank == i) Then
          iGr =iGr+1
          MyGroup(iGr)=i
          MyGroup_list(i+1)=iGr
       End If
       If (MyRank /= i) Then
          If (Mr_Nbr(Myrank,i) == 1) Then
             iNe =iNe+1
             iGr =iGr+1
             Myneighbor(iNe)=i
             MyGroup(iGr)=i
             MyNeighbor_list(i+1)=iNe
             MyGroup_list(i+1)=iGr
          End If
       End If
    End Do
    
    DeAllocate(IS_Nbr)
    DeAllocate(Mr_Nbr)
    DeAllocate(Total)

  End Subroutine G_PU_2D_Scal_Dist
  
!!! PU function communication for neighbor
  Subroutine Communicate_PU()
    Integer                                       :: i
    
    Integer, Dimension(:), Pointer                :: Request
    Integer, Dimension(:,:), Pointer              :: Status
    Real(Kind = Kr), Dimension(:), Pointer        :: TmpPU

    Allocate(TmpPU(Geom%Num_Nodes*MySize_Neighbor))
    Allocate(Request(MySize_Neighbor))
    Allocate(Status(MPI_STATUS_SIZE,MySize_Neighbor))
    
    Do i = 1, MySize_Neighbor
       Call MPI_IRECV(TmpPU((i-1)*Geom%Num_Nodes+1),Geom%Num_Nodes, &
            & MPI_DOUBLE_PRECISION, Myneighbor(i), &
            & Myneighbor(i)*NumProcs+MyRank+1,PETSC_COMM_WORLD, &
            & Request(i),iErr)  
    EndDo
    
    Do i = 1, MySize_Neighbor
       Call MPI_SEND(MyPU,Geom%Num_Nodes, MPI_DOUBLE_PRECISION, &
            & Myneighbor(i), MyRank*NumProcs+Myneighbor(i)+1, &
            & PETSC_COMM_WORLD,iErr)
    End Do
    Call MPI_WAITALL(MySize_Neighbor,Request,Status,iErr)
    
    Allocate(GroupPU(Geom%Num_Nodes,MySize_Neighbor))
    Do i= 1, MySize_Neighbor
       GroupPU(:,MyNeighbor_list(MyNeighbor(i)+1))= &
            & TmpPU((i-1)*Geom%Num_Nodes+1:i*Geom%Num_Nodes)
    End Do
    DeAllocate (TmpPU)
    DeAllocate (Request)
    DeAllocate (Status)
  End Subroutine Communicate_PU
  
!! Assemble Coarse Matrix

  Subroutine Assemb_Coarse_Seq()
    Vec, Dimension(:), Pointer                       :: MyMRdotPU_Local
    Real(Kind = Kr), Dimension(:,:), Pointer         :: MyMRdotMyPU

    Vec                                              :: MyPUVec_Local
    Vec                                              :: MyMRdotMyPU_Local
    Integer                                          :: i, j

    PetscScalar, Dimension(:), Pointer               :: MyMRdotMyPU_Local_Ptr
    Integer, Dimension(:), Pointer                   :: Conum
    IS                                               :: MRCO_IS

    Allocate(MyMRdotMyPU(Geom%Num_Nodes,MySize_Group))
    Allocate(MyMRdotPU_Local(MySize_Group))
    MyMRdotMyPU=0.0_Kr
    Allocate(Conum(MySize_Node))
    Conum   = (/ (i, i = 0, MySize_Node-1)/)

    Call VecCreateSeq(PETSC_COMM_SELF, MySize_Node_Ovlp, MyPUVec_Local, iErr)
    Call VecCreateSeq(PETSC_COMM_SELF, MySize_Node_Ovlp, MyMRdotMyPU_Local,   &
         & iErr)
    Call VecSetValues(MyPUVec_Local, MySize_Node, Conum,MyPU(MyNode+1),       &
         & ADD_VALUES, iErr)
    Call VecAssemblyBegin(MyPUVec_Local, iErr)
    Call VecAssemblyEnd(MyPUVec_Local, iErr)
    Call MatMult(MyMR, MyPUVec_Local, MyMRdotMyPU_Local, iErr) 
    Call VecGetArrayF90(MyMRdotMyPU_Local, MyMRdotMyPU_Local_Ptr, iErr)
    MyMRdotMyPU(MyNode_Ovlp+1,MyGroup_list(MyRank+1))=MyMRdotMyPU_Local_Ptr
    Call VecRestoreArrayF90(MyMRdotMyPU_Local, MyMRdotMyPU_Local_Ptr, iErr)
    Call VecDestroy(MyPUVec_Local,iErr)
    Call VecDestroy(MyMRdotMyPU_Local,iErr)
    
    Do i= 1, MySize_Neighbor
       Call VecCreateSeq(PETSC_COMM_SELF,MySize_Node_Ovlp,MyPUVec_Local, iErr)
       Call VecCreateSeq(PETSC_COMM_SELF,MySize_Node_Ovlp,MyMRdotMyPU_Local, iErr)
       Call VecSetValues(MyPUVec_Local,MySize_Node,Conum,GroupPU(MyNode+1, &
            & MyNeighbor_list(MyNeighbor(i)+1)),ADD_VALUES, iErr)
       Call VecAssemblyBegin(MyPUVec_Local,iErr)
       Call VecAssemblyEnd(MyPUVec_Local,iErr)
       Call MatMult(MyMR,MyPUVec_Local,MyMRdotMyPU_Local, iErr) 
       Call VecGetArrayF90(MyMRdotMyPU_Local, MyMRdotMyPU_Local_Ptr, iErr)
       MyMRdotMyPU(MyNode_Ovlp+1,MyGroup_list(MyNeighbor(i)+1))=MyMRdotMyPU_Local_Ptr
       Call VecRestoreArrayF90(MyMRdotMyPU_Local, MyMRdotMyPU_Local_Ptr, iErr)
       Call VecDestroy(MyPUVec_Local,iErr)
       Call VecDestroy(MyMRdotMyPU_Local,iErr)
    End Do
!!$    If (MyRank == 1) Then
!!$       Print*, 'O.K 2'
!!$    End If
    If (MyRank == 0) Then
!!$       Print*, 'O.K 2.1'
       Call MatCreateMPIAIJ(PETSC_COMM_WORLD, NumProcs, NumProcs,  NumProcs,    &
            &  NumProcs, 30, PETSC_NULL_INTEGER, 15, PETSC_NULL_INTEGER,    &
            & MRCO, iErr)
!!$       Print*, 'O.K 2.2'
    Else
!!$       Print*, 'Rank :',MyRank ,'O.K 2.3'
       Call MatCreateMPIAIJ(PETSC_COMM_WORLD,  0,  0,  NumProcs,    &
            &  NumProcs, 0, PETSC_NULL_INTEGER, 0, PETSC_NULL_INTEGER,    &
            & MRCO, iErr)
!!$       Print*, 'Rank :',MyRank ,'O.K 2.4'
    End If
!!$    If (MyRank == 0) Then
!!$       Print*, 'O.K 3'
!!$    End If
    Call MatSetFromOptions(MRCO, iErr)
    Do j = 1, MySize_Group
       Call MatSetValue(MRCO,MyRank, MyGroup(j), Dot_Product(MyPU(MyNode_Ovlp+1),&
            & MyMRdotMyPU(MyNode_Ovlp+1,MyGroup_list(MyGroup(j)+1))), ADD_VALUES,&
            & iErr)
       Do i= 1, MySize_Neighbor
          Call MatSetValue(MRCO,MyNeighbor(i), MyGroup(j), &
          & Dot_Product(GroupPU(MyNode_Ovlp+1,MyNeighbor_list(MyNeighbor(i)+1)), &
          & MyMRdotMyPU(MyNode_Ovlp+1,MyGroup_list(MyGroup(j)+1))), ADD_VALUES,  &
          & iErr)
       End Do
    End Do
!!$    If (MyRank == 0) Then
!!$       Print*, 'O.K 4'
!!$    End If
    Call MatAssemblyBegin(MRCO, MAT_FINAL_ASSEMBLY, iErr)
    Call MatAssemblyEnd(MRCO, MAT_FINAL_ASSEMBLY,iErr)
    Allocate(MRCOnum(NumProcs))
    MRCOnum = (/ (i, i = 0, NumProcs-1)/)
    Call ISCreateGeneral(PETSC_COMM_WORLD, NumProcs, MRCOnum, MRCO_IS, iErr)

    If (MyRank==0) Then
       Call MatGetSubMatrices(MRCO, 1, MRCO_IS, MRCO_IS, MAT_INITIAL_MATRIX,  &
            & MRCO_S, iErr)
    Else
       Call MatGetSubMatrices(MRCO, 0, MRCO_IS, MRCO_IS, MAT_INITIAL_MATRIX,  &
            & MRCO_S, iErr)
    EndIf

    Call ISDestroy(MRCO_IS, iErr)
    DeAllocate(MyMRdotMyPU)
    DeAllocate(MyMRdotPU_Local)
    DeAllocate(Conum)
  End Subroutine Assemb_Coarse_Seq
  
  !! Initialize Data for coarse problem
  
  Subroutine Coarse_Solve_Prep()
    Integer                                  :: i, iErr
    Integer, Dimension(:), Pointer           :: MyPUSecIs 
    Integer, Dimension(:), Pointer           :: CoSOL_Master_IDx 

    Call VecCreateSeq(PETSC_COMM_SELF, MySize_Node_Ovlp_Metis, MyPUSec, iErr)

    Allocate(MyPUSecIs(MySize_Node_Ovlp_Metis))
    MyPUSecIs   = (/ (i, i = 0, MySize_Node_Ovlp_Metis-1)/)
    Call VecSetValues(MyPUSec, MySize_Node_Ovlp_Metis, MyPUSecIs,             &
         & MyPU(MyNode_Ovlp+1),INSERT_VALUES,iErr)
    Call VecAssemblyBegin(MyPUSec,ierr) 
    Call VecAssemblyEnd(MyPUSec,ierr)

    Allocate(Renum3(MySize_Node))
    Renum3   = (/ (i, i = MyIdxMin_Node_Ov, MyIdxMax_Node_Ov-1)/)

    Call VecGhostGetLocalForm(RHS_Ov,MyRHS_Ov,ierr) 
    Call VecGetOwnershipRange(RHS_Ov, MyIdxMin_Node_Ov, MyIdxMax_Node_Ov,ierr)
    Call VecGhostRestoreLocalForm(RHS_Ov,MyRHS_Ov,ierr)
    If (MyRank == 0) Then
       Allocate(CoSOL_Master_IDx(Geom%Num_Nodes))
       CoSOL_Master_IDx=(/ (i, i = 0, Geom%Num_Nodes-1)/)
    End If
    Call VecDuplicate(SOL_Master, CoSOL_Master,iErr)
    Call VecDuplicate(SOL, CoSOL,iErr)
    Call VecCreateMPI(PETSC_COMM_WORLD, 1,NumProcs, CoOut, iErr)

    Call VecDuplicate(CoOut, SOL_CoOut,iErr)

    If (MyRank == 0) Then
       Call VecCreateSeq(PETSC_COMM_SELF, NumProcs, CoIn_Seq, iErr)
    Else
       Call VecCreateSeq(PETSC_COMM_SELF, 0, CoIn_Seq, iErr)
    End If
    Call VecScatterCreateToZero(CoOut, CoVecScatter, CoIn_Seq, iErr)
    Call VecDuplicate(CoIn_Seq, SOL_CoIn_Seq,iErr)
    Call VecGhostGetLocalForm(CoPU_Ov,MyCoPU_Ov,ierr)

    DeAllocate(MyPUSecIs)
    If (MyRank == 0) Then
       DeAllocate(CoSOL_Master_IDx)
    End If
  End Subroutine Coarse_Solve_Prep
  
!! Coarse Solve Routine

  Subroutine Coarse_Solve_Seq(RHS, SOL)
    Vec                                        :: RHS, SOL 
    Integer                                    :: iErr
    PetscScalar                                :: vv, zz
    
    Call VecGetArrayF90(RHS, RHS_Ptr, iErr)
    !!    Call VecGhostGetLocalForm(RHS_Ov,MyRHS_Ov,ierr) 
    Call VecSetValues(RHS_Ov, MySize_Node, Renum3,                            &
         & RHS_Ptr(Renum3-MyIdxMin_Node_Ov+1),INSERT_VALUES,ierr)
    Call VecAssemblyBegin(RHS_Ov,ierr) 
    Call VecAssemblyEnd(RHS_Ov,ierr)
    Call VecRestoreArrayF90(RHS, RHS_Ptr, iErr)

!!$    Call VecGhostUpdateBegin(RHS_Ov,INSERT_VALUES,SCATTER_FORWARD,ierr) 
!!$    Call VecGhostUpdateEnd(RHS_Ov,INSERT_VALUES,SCATTER_FORWARD,ierr)
    Call VecDot(MyRHS_Ov,MyPUSec,vv,ierr)
    !!    Call VecGhostRestoreLocalForm(RHS_Ov,MyRHS_Ov,ierr)
    Call VecSetValue(CoOut, MyRank, vv, INSERT_VALUES,ierr)
    Call VecAssemblyBegin(CoOut,ierr) 
    Call VecAssemblyEnd(CoOut,ierr)

    Call VecScatterBegin(CoOut, CoIn_Seq, INSERT_VALUES, SCATTER_FORWARD,     &
         & CoVecScatter, iErr )
    Call VecScatterEnd(CoOut, CoIn_Seq, INSERT_VALUES, SCATTER_FORWARD,       &
         & CoVecScatter, iErr )
    If (MyRank == 0) Then
       Call KSPSetup(KSP_CO_Seq, iErr)
       Call KSPSolve(KSP_CO_Seq, CoIn_Seq, SOL_CoIn_Seq, iErr)
    End If
    Call VecScatterBegin(SOL_CoIn_Seq, SOL_CoOut, INSERT_VALUES,              &
         & SCATTER_REVERSE, CoVecScatter, iErr)
    Call VecScatterEnd(SOL_CoIn_Seq, SOL_CoOut, INSERT_VALUES,                &
         &SCATTER_REVERSE, CoVecScatter, iErr )
    !!    Call VecGhostGetLocalForm(CoPU_Ov,MyCoPU_Ov,ierr) 

    !!! I am sure that this block can also be much simplified
    Call VecGetArrayF90(MyCoPU_Ov, MyCoPU_Ov_Ptr, iErr)
    Call VecGetArrayF90(SOL_CoOut, SOL_CoOut_Ptr, iErr)
    MyCoPU_Ov_Ptr(MyNode_Ovlp_list(MyNode_Ovlp)+1) =                          &
         & MyPU(MyNode_Ovlp+1)*SOL_CoOut_Ptr(1)
    Call VecRestoreArrayF90(SOL_CoOut, SOL_CoOut_Ptr, iErr)
    Call VecRestoreArrayF90(MyCoPU_Ov, MyCoPU_Ov_Ptr, iErr)

!!$    Call VecGhostUpdateBegin(CoPU_Ov, ADD_VALUES,SCATTER_Reverse, ierr)
!!$    Call VecGhostUpdateEnd(CoPU_Ov, ADD_VALUES,SCATTER_Reverse, ierr)
    !!    Call VecGhostRestoreLocalForm(CoPU_Ov,MyCoPU_Ov,ierr)
    Call VecCopy(CoPU_Ov, SOL, iErr)
  End Subroutine Coarse_Solve_Seq
  
!! OBDD SETUP Routine
  Subroutine SET_OBDD(Dummy, InVec, OutVec, iErr)
    Integer                                  :: Dummy
    Vec                                      :: InVec, OutVec
    Integer                                  :: iErr


    Vec                                      :: MyInSOL_Ov
    PetscScalar, Dimension(:), Pointer       :: InVec_Ptr
    PetscLogDouble                           :: ASMTotTS, ASMTotTF, ASMTotT
    PetscLogDouble                           :: ASMSetTS, ASMSetTF, ASMSetT
    PetscLogDouble                           :: ASMTS, ASMTF
    Real(Kind = Kr), Dimension(:), Pointer   :: Coarse_SOL
    PetscScalar                              :: v1, v2, zz

    If (ASMiter == 0 ) Then
       Call PetscLogStagePush(stages(10), iErr)
    Else
       Call PetscLogStagePush(stages(11), iErr)
    End If
    !! STEP 1 : Coarse Step
    Call PetscGetTime(TTTS, iErr)
    Call PetscGetTime(TP1S, iErr)
    
    Call PetscGetTime(T1S, iErr)                   !! T1 --- Get Input
    Call PetscGetTime(ASMTS, iErr)
    Call VecGetArrayF90(InVec, InVec_Ptr, iErr) !! Get Input Data
    Call VecSetValues(RES_Ov,MySize_Node,Renum3,&
         & InVec_Ptr(Renum3-MyIdxMin_Node_Ov+1), &
         & INSERT_VALUES,ierr)                   !! Get Residual with ghost
    Call VecRestoreArrayF90(InVec, InVec_Ptr, iErr)
    Call VecAssemblyBegin(RES_Ov,ierr) 
    Call VecAssemblyEnd(RES_Ov,ierr)            !! Residual
    Call PetscGetTime(T1F, iErr)                  !! T1 End 
    
    Call PetscGetTime(T2S, iErr)                  !! T2 --- Get Res
    Call VecGhostUpdateBegin(RES_Ov,INSERT_VALUES,SCATTER_FORWARD,ierr) 
    Call VecGhostUpdateEnd(RES_Ov,INSERT_VALUES,SCATTER_FORWARD,ierr)
    Call PetscGetTime(T2F, iErr)                  !! T2 End
    
    !!*P1       Call PetscGetTime(T3S, iErr)                  !! T3 --- Coarse Rhs
    !!*P1        Call VecDot(MyRES_Ov,MyPUSec,v1,ierr)
    !!*P1            Call VecSetValue(CoOut,MyRank,v1,INSERT_VALUES,ierr)
    !!*P1    Call VecAssemblyBegin(CoOut,ierr) 
    !!*P1    Call VecAssemblyEnd(CoOut,ierr)
    !!*P1   Call PetscGetTime(T3F, iErr)                  !! T3 End
    
    !!*P1   Call PetscGetTime(T4S, iErr)                  !! T4 --- Collect in O
    !!*P1    Call VecScatterBegin(CoOut, CoIn,  INSERT_VALUES,SCATTER_FORWARD, CoVecScatter, iErr )
    !!*P1    Call VecScatterEnd(CoOut, CoIn,  INSERT_VALUES,SCATTER_FORWARD, CoVecScatter, iErr )
    !!*P1   Call PetscGetTime(T4F, iErr)                  !! T4 End
    
    !!*P1   Call PetscGetTime(T5S, iErr)                  !! T5 --- Coarse Sol
    !!*P1    If (MyRank == 0) Then
    !!*P1       Call VecGetArrayF90(CoIn, CoIn_Ptr, iErr)
    !!*P1       Call VecSetValues(CoIn_Seq,NumProcs,MRCOnum,CoIn_Ptr,INSERT_VALUES, iErr)
    !!*P1       Call VecRestoreArrayF90(CoIn, CoIn_Ptr, iErr)
    !!*P1       Call VecAssemblyBegin(CoIn_Seq,ierr) 
    !!*P1       Call VecAssemblyEnd(CoIn_Seq,ierr)
    !!*P1       Call KSPSolve(KSP_CO_Seq,CoIn_Seq, SOL_CoIn_Seq,iErr)
    !!*P1       Call VecGetArrayF90(SOL_CoIn_Seq, SOL_CoIn_Seq_Ptr, iErr)
    !!*P1       Call VecSetValues(SOL_CoIn,NumProcs,MRCOnum,SOL_CoIn_Seq_Ptr, INSERT_VALUES, iErr)
    !!*P1       Call VecRestoreArrayF90(SOL_CoIn_Seq, SOL_CoIn_Seq_Ptr, iErr)
    !!*P1     End If
    !!*P1     Call VecAssemblyBegin(SOL_CoIn,ierr) 
    !!*P1     Call VecAssemblyEnd(SOL_CoIn,ierr)           
    !!*P1   Call PetscGetTime(T5F, iErr)                   !! T5 End
    
    !!*P1   Call PetscGetTime(T6S, iErr)                   !! T6 --- To All
    !!*P1     Call VecScatterBegin(SOL_CoIn, SOL_CoOut,  INSERT_VALUES,SCATTER_REVERSE, CoVecScatter, iErr )
    !!*P1     Call VecScatterEnd(SOL_CoIn, SOL_CoOut,INSERT_VALUES,SCATTER_REVERSE, CoVecScatter, iErr )
    !!*P1   Call PetscGetTime(T6F, iErr)                   !! T6 End
    
    !!*P1   Call PetscGetTime(T7S, iErr)                   !! T7 --- Fine Space
    !!*P1     Call VecGetArrayF90(MyCoPU_Ov, MyCoPU_Ov_Ptr, iErr)
    !!*P1      Call VecGetArrayF90(SOL_CoOut, SOL_CoOut_Ptr, iErr)
    !!*P1     MyCoPU_Ov_Ptr(MyNode_Ovlp_list(MyNode_Ovlp)+1) = MyPU(MyNode_Ovlp+1)*SOL_CoOut_Ptr(1)
    !!*P1     Call VecRestoreArrayF90(SOL_CoOut, SOL_CoOut_Ptr, iErr)
    !!*P1     Call VecRestoreArrayF90(MyCoPU_Ov, MyCoPU_Ov_Ptr, iErr)
    !!*P1   Call PetscGetTime(T7F, iErr)                   !! T7 End
    
    !!*P1   Call PetscGetTime(T8S, iErr)                   !! T8 --- Fine Space Update
    !!*P1     Call VecGhostUpdateBegin(CoPU_Ov,ADD_VALUES,SCATTER_Reverse,ierr)
    !!*P1     Call VecGhostUpdateEnd(CoPU_Ov,ADD_VALUES,SCATTER_Reverse,ierr)
    !!*P1     Call VecGhostUpdateBegin(CoPU_Ov,INSERT_VALUES,SCATTER_FORWARD,ierr) !!New
    !!*P1     Call VecGhostUpdateEnd(CoPU_Ov,INSERT_VALUES,SCATTER_FORWARD,ierr) !!New
    !!*P1   Call PetscGetTime(T8F, iErr)                   !! T8 End
    
    !!*P1   Call PetscGetTime(T9S, iErr)                   !! T9 --- Update Res
    !!*P1     Call MatMult(MyMR_Neu, MyCoPU_Ov, MyUpdate_Ov, iErr)
    !!*P1     Call VecAXPY(none,MyUpdate_Ov, MyRES_Ov,iErr)
    !!    Call MatMult(MR, CoPU_Ov, MRCoSOL,iErr)
    !!    Call VecAXPY(none,MRCoSOL,RES_Ov,iErr)
    !!*P1   Call PetscGetTime(T9F, iErr)                   !! T9 End
    
    Call PetscGetTime(TP1F, iErr)
    
    !! Part 2 : Additive with PUs
    Call PetscGetTime(TP2S, iErr)
    
    Call PetscGetTime(T10S, iErr)                  !! T10 --- Local Res
    !!    Call VecGhostUpdateBegin(RES_Ov,INSERT_VALUES,SCATTER_FORWARD,ierr) 
    !!    Call VecGhostUpdateEnd(RES_Ov,INSERT_VALUES,SCATTER_FORWARD,ierr)
    Call VecGetArrayF90(MyRES_Ov, MyRES_Ov_Ptr, iErr)
    MyRES_Ov_Ptr(MyNode_Ovlp_list(MyNode_Ovlp)+1) = &
         & MyPU2(MyNode_Ovlp+1)*MyRES_Ov_Ptr(MyNode_Ovlp_list(MyNode_Ovlp)+1)
    Call VecRestoreArrayF90(MyRES_Ov, MyRES_Ov_Ptr, iErr)
    Call PetscGetTime(T10F, iErr)                  !! T10 End
    
    Call PetscGetTime(T11S, iErr)                  !! T11 --- Local Sol  
    If (Num_Method == 1 ) Then
       Call KSPSolve(KSP_MyMR, MyRES_Ov, MyIMP_Ov,iErr)
    End If
    If (Num_Method == 2 ) Then
       Call KSPSolve(KSP_MyMR_Neu, MyRES_Ov, MyIMP_Ov,iErr)
    End If
    Call PetscGetTime(T11F, iErr)                  !! T11 End
    
    Call PetscGetTime(T12S, iErr)                  !! T12 --- Balance
    Call VecGetArrayF90(MyIMP_Ov, MyIMP_Ov_Ptr, iErr)
    MyIMP_Ov_Ptr(MyNode_Ovlp_list(MyNode_Ovlp)+1) = &
         & MyPU2(MyNode_Ovlp+1)*MyIMP_Ov_Ptr(MyNode_Ovlp_list(MyNode_Ovlp)+1)
    Call VecRestoreArrayF90(MyIMP_Ov, MyIMP_Ov_Ptr, iErr)
    Call PetscGetTime(T12F, iErr)                  !! T12 End
    
    Call PetscGetTime(T13S, iErr)                  !! T13 --- Get Imp
    Call VecGhostUpdateBegin(IMP_Ov,ADD_VALUES,SCATTER_Reverse,ierr)
    Call VecGhostUpdateEnd(IMP_Ov,ADD_VALUES,SCATTER_Reverse,ierr)
    Call VecGhostUpdateBegin(IMP_Ov,INSERT_VALUES,SCATTER_FORWARD,ierr) !!New
    Call VecGhostUpdateEnd(IMP_Ov,INSERT_VALUES,SCATTER_FORWARD,ierr) !!New
    Call PetscGetTime(T13F, iErr)                  !! T13 End
    
    Call PetscGetTime(T14S, iErr)                  !! T14 --- Get Res 
    Call MatMult(MyMR_Neu, MyIMP_Ov, MyUpdate_Ov, iErr)
    Call VecAXPBY(MyRES_Ov,none, zz,MyUpdate_Ov,iErr)
    !!   Call MatMult(MR, IMP_Ov, MRIMP_Ov,iErr)
    !!   Call VecAXPBY(none, zz,MRIMP_Ov,MyRes_Ov,iErr) 
    !!   Call VecAXPY(pone, CoPU_Ov, IMP_Ov, iErr)
    Call PetscGetTime(T14F, iErr)                  !! T14 End
    Call PetscGetTime(TP2F, iErr) 
    
    
    !! Part 3 : Second Coarse Step
    Call PetscGetTime(TP3S, iErr)
    !!  Call PetscGetTime(T15S, iErr)                  !! T15 --- Update Ghost 
    !!    Call VecGhostUpdateBegin(RES_Ov,INSERT_VALUES,SCATTER_FORWARD,ierr) 
    !!    Call VecGhostUpdateEnd(RES_Ov,INSERT_VALUES,SCATTER_FORWARD,ierr)
    !!  Call PetscGetTime(T15F, iErr)                  !! T15 End
    
    Call PetscGetTime(T16S, iErr)                  !! T16 --- Get CoRhs
    Call VecDot(MyRES_Ov,MyPUSec,v2,ierr)
    Call VecSetValue(CoOut,MyRank,v2,INSERT_VALUES,ierr)
    Call VecAssemblyBegin(CoOut,ierr) 
    Call VecAssemblyEnd(CoOut,ierr)
    Call PetscGetTime(T16F, iErr)                  !! T16 End
    
    Call PetscGetTime(T17S, iErr)                  !! T17 --- Set to 0
    Call VecScatterBegin(CoOut, CoIn_Seq, INSERT_VALUES,SCATTER_FORWARD,        &
         & CoVecScatter, iErr )
    Call VecScatterEnd(CoOut, CoIn_Seq,  INSERT_VALUES,SCATTER_FORWARD,         &
         & CoVecScatter, iErr )
    Call PetscGetTime(T17F, iErr)                  !! T17 End
    
    Call PetscGetTime(T18S, iErr)                  !! T18 --- Coarse Sol
    If (MyRank == 0) Then
       Call KSPSolve(KSP_CO_Seq,CoIn_Seq, SOL_CoIn_Seq,iErr)
    End If
    Call PetscGetTime(T18F, iErr)                  !! T18 End
    
    Call PetscGetTime(T19S, iErr)                  !! T19 --- To All
    Call VecScatterBegin(SOL_CoIn_Seq, SOL_CoOut, INSERT_VALUES,                &
         & SCATTER_REVERSE, CoVecScatter, iErr )
    Call VecScatterEnd(SOL_CoIn_Seq, SOL_CoOut, INSERT_VALUES,                  &
         & SCATTER_REVERSE, CoVecScatter, iErr )
    Call PetscGetTime(T19F, iErr)                  !! T19 End
    
    Call PetscGetTime(T20S, iErr)                  !! T20 --- Fine Space
    Call VecGetArrayF90(MyCoPU_Ov, MyCoPU_Ov_Ptr, iErr)
    Call VecGetArrayF90(SOL_CoOut, SOL_CoOut_Ptr, iErr)
    MyCoPU_Ov_Ptr(MyNode_Ovlp_list(MyNode_Ovlp)+1)                            &
         & = MyPU(MyNode_Ovlp+1)*SOL_CoOut_Ptr(1)
    Call VecRestoreArrayF90(SOL_CoOut, SOL_CoOut_Ptr, iErr)
    Call VecRestoreArrayF90(MyCoPU_Ov, MyCoPU_Ov_Ptr, iErr)
    Call VecGhostUpdateBegin(CoPU_Ov,ADD_VALUES,SCATTER_Reverse,ierr)
    Call VecGhostUpdateEnd(CoPU_Ov,ADD_VALUES,SCATTER_Reverse,ierr)
    Call PetscGetTime(T20F, iErr)                  !! T20 End
    
    Call PetscGetTime(T21S, iErr)                  !! T21 --- Final
    Call VecAXPY(CoPU_Ov, pone, IMP_Ov, iErr)
    Call VecCopy(IMP_Ov, OutVec, iErr)
    Call PetscGetTime(T21F, iErr)                  !! T21 End
    
    Call PetscGetTime(TP3F, iErr)
    Call PetscGetTime(ASMTF, iErr)
    Call PetscGetTime(TTTF, iErr)
    If (ASMiter == 0 ) Then
       ASMSETUPT = ASMTF - ASMTS
       !!      Call PetscLogStagePop(iErr)
    ElseIf (ASMiter == 1) Then
       ASMT = ASMTF - ASMTS
       TT1= T1F-T1S
       TT2= T2F-T2S
       TT3= T3F-T3S
       TT4= T4F-T4S
       TT5= T5F-T5S
       TT6= T6F-T6S
       TT7= T7F-T7S
       TT8= T8F-T8S
       TT9= T9F-T9S
       TT10= T10F-T10S
       TT11= T11F-T11S
       TT12= T12F-T12S
       TT13= T13F-T13S
       TT14= T14F-T14S
       TT15= T15F-T15S
       TT16= T16F-T16S
       TT17= T17F-T17S
       TT18= T18F-T18S
       TT19= T19F-T19S
       TT20= T20F-T20S
       TT21= T21F-T21S
       TP1 = TP1F-TP1S
       TP2 = TP2F-TP2S
       TP3 = TP3F-TP3S
       TTTT= TTTF-TTTS
    Else
       ASMT = ASMTF - ASMTS+ASMT
    End If
    !%%      Print*, 'ASMiter:', ASMiter
    ASMiter=ASMiter+1
    Call PetscLogStagePop(iErr)
  End Subroutine SET_OBDD
  
  !! Setup ASM for pure CG case
  Subroutine SET_ASM(Dummy, InVec, OutVec, iErr)
    Vec                                      :: InVec, MyInSOL_Ov
    Vec                                      :: OutVec, MyOutSOL_Ov
    PetscScalar                              :: none
    PetscScalar                              :: pone
    Integer                                  :: i, Numiter, iErr
    Integer                                  :: Dummy
    PetscLogDouble              :: ASMTotTS, ASMTotTF, ASMTotT
    PetscLogDouble              :: ASMSetTS, ASMSetTF, ASMSetT
    PetscLogDouble              :: ASMTS, ASMTF

    none =-1.0
    pone = 1.0
    Call PetscGetTime(ASMTS, iErr)
    Call VecGhostGetLocalForm(RES_Ov,MyRES_Ov,ierr) 
    Call VecGhostGetLocalForm(IMP_Ov,MyIMP_Ov,ierr)
    Call VecCopy(InVec, RES_Ov, iErr)
    Call VecGhostUpdateBegin(RES_Ov,INSERT_VALUES,SCATTER_FORWARD,ierr) 
    Call VecGhostUpdateEnd(RES_Ov,INSERT_VALUES,SCATTER_FORWARD,ierr)
    Call VecGetArrayF90(MyRES_Ov, MyRES_Ov_Ptr, iErr)
    MyRES_Ov_Ptr(MyNode_Ovlp_list(MyNode_Bd))=0.0_Kr
    Call VecRestoreArrayF90(MyRES_Ov, MyRES_Ov_Ptr, iErr)
!!    Call KSPSetRHS(KSP_MyMR, MyRES_Ov, iErr)
!!    Call KSPSetSolution(KSP_MyMR, MyIMP_Ov, iErr)
    Call KSPSetup(KSP_MyMR, iErr)
    Call KSPSolve(KSP_MyMR, MyRES_Ov, MyIMP_Ov,iErr)
    Call VecGhostUpdateBegin(IMP_Ov,ADD_VALUES,SCATTER_Reverse,ierr)
    Call VecGhostUpdateEnd(IMP_Ov,ADD_VALUES,SCATTER_Reverse,ierr) 
    Call VecCopy(IMP_Ov, OutVec, iErr)
    Call VecGhostRestoreLocalForm(RES_Ov,MyRES_Ov,ierr) 
    Call VecGhostRestoreLocalForm(IMP_Ov,MyIMP_Ov,ierr)
    Call PetscGetTime(ASMTF, iErr)
    If (ASMiter == 0 ) Then
      ASMSETUPT = ASMTF - ASMTS
    ElseIF (ASMiter == 1) Then
      ASMT = ASMTF - ASMTS
    Else
      ASMT = ASMTF - ASMTS+ASMT
!     ASMT = ASMTF - ASMTS
    End If
      ASMiter=ASMiter+1
  End Subroutine SET_ASM

!! KSP Setup for CG
  Subroutine Init_KSP_OLD()
    integer                     :: user_defined_pc
    Real(Kind = Kr)             :: Tol 
    PetscTruth                  :: flg
                    

    Call KSPCreate(PETSC_COMM_WORLD, KSP_MR, iErr)
    Call KSPSetOperators(KSP_MR, MR, MR, SAME_NONZERO_PATTERN, iErr)
    Call KSPSetType(KSP_MR, KSPCG, iErr)
!    Call KSPSetInitialGuessNonzero(KSP_MR, PETSC_TRUE, iErr)
    Call KSPGetPC(KSP_MR, PC_MR, iErr)
!    Call PCSetType(PC_MR, PCLU, iErr)
    Call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-user_defined_pc',      &
      &                    user_defined_pc,ierr)
    Call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-pc_name',PC_Name,flg,ierr)
!    Call PCSetType(PC_MR, PCBJACOBI, iErr)
!    Call PCSetType(PC_MR, PCJACOBI, iErr)
    Tol = 0.000001
    Call KSPSetTolerances(KSP_MR, Tol,                               &
           & PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION,  &
           & PETSC_DEFAULT_INTEGER, iErr)
!    Call KSPSetTolerances(KSP_MR, Tol,                               &
!           & PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION,  &
!           & 1, iErr)
!    Tol = 0.00000000001
    Call KSPSetFromOptions(KSP_MR, iErr)
    Call KSPSetup(KSP_MR, iErr)
  End Subroutine Init_KSP_OLD

!! KSP Setup for DIRI and OBDD ---- I left Alternatives for referece
  
  Subroutine Init_KSP()
    integer                     :: user_defined_pc
    PetscInt                    :: levels
    Real(Kind = Kr)             :: Tol
    PetscScalar                 :: one, alpha
    PetscTruth                  :: flg
     one = 1.0

    Call KSPCreate(PETSC_COMM_WORLD, KSP_MR, iErr)
    Call KSPSetOperators(KSP_MR, MR, MR, SAME_NONZERO_PATTERN, iErr)
    Call KSPSetType(KSP_MR, KSPCG, iErr)
    Call KSPSetInitialGuessNonzero(KSP_MR, PETSC_TRUE, iErr)
    Call KSPGetPC(KSP_MR, PC_MR, iErr)
    Call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-user_defined_pc',      &
      &                    user_defined_pc,ierr)
 !!  Call PCSetType(PC_MR, PCJACOBI, iErr)
    Call PCSetType(PC_MR, PCSHELL, iErr)
 !!  Call PCShellSetApply(PC_MR,SET_ASM,PETSC_NULL_OBJECT,iErr)
    Call PCShellSetApply(PC_MR,SET_OBDD,PETSC_NULL_OBJECT,iErr)
    Tol = 0.000001
    Call KSPSetTolerances(KSP_MR, Tol,                               &
         & PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION,  &
         & PETSC_DEFAULT_INTEGER, iErr)
!    Call KSPSetTolerances(KSP_MR, Tol,                               &
!           & PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION,  &
!           & 4 , iErr)
    If (MyRank == 0 ) Then
       Call KSPCreate(PETSC_COMM_SELF, KSP_CO_Seq, iErr)
       Call KSPSetOperators(KSP_CO_Seq, MRCO_S, MRCO_S, SAME_NONZERO_PATTERN, &
            & iErr)
       Call KSPSetType(KSP_CO_Seq,  KSPPREONLY, iErr)
       !!     Call KSPSetType(KSP_CO_Seq,  KSPCG, iErr)
       Call KSPGetPC(KSP_CO_Seq, PC_CO_Seq, iErr)
       Call PCSetType(PC_CO_Seq, PCLU, iErr)
       Call PCFactorSetZeroPivot(PC_CO_Seq, 1.0D-30)

       !!     Tol = 0.1
!!!     Call KSPSetTolerances(KSP_CO_Seq, Tol,                               &
!!!           & PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION,  &
!!!           & PETSC_DEFAULT_INTEGER, iErr)
      Call KSPSetFromOptions(KSP_CO_Seq, iErr)
    End If


!! Dirichlet Interface Case
    If (Num_Method == 1 ) Then
      Call KSPCreate(PETSC_COMM_SELF, KSP_MyMR, iErr)
      Call KSPSetOperators(KSP_MyMR, MyMR, MyMR, SAME_NONZERO_PATTERN, iErr)

!! Direct Local Solver
      Call KSPSetType(KSP_MyMR, KSPPREONLY, iErr)
      Call KSPGetPC(KSP_MyMR, PC_MyMR, iErr)
      Call PCSetType(PC_MyMR, PCLU, iErr)
!!     Call PCSetType(PC_MyMR, PCICC, iErr)
!!     Call PCSetType(PC_MyMR, PCILU, iErr)
!!     levels =4
!!     Call PCILUSetLevels(PC_MyMR, levels, iErr)
!! Iterative Local Solver
!!*    Call KSPSetType(KSP_MyMR, KSPCG, iErr)
!!     Call KSPSetInitialGuessNonzero(KSP_MyMR, PETSC_TRUE, iErr)
!!*    Call KSPGetPC(KSP_MyMR, PC_MyMR, iErr)
!!*    Call PCSetType(PC_MyMR, PCJACOBI, iErr)
!!     Call KSPSetFromOptions(KSP_MyMR, iErr)
!!     Tol = 0.1
!!     Call KSPSetTolerances(KSP_MyMR, Tol,                               &
!!           & PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION,  &
!!          & PETSC_DEFAULT_INTEGER, iErr)
    End If

    If (Num_Method == 2 ) Then
!! Neumann Interface Case
      Call KSPCreate(PETSC_COMM_SELF, KSP_MyMR_Neu, iErr)
      Call KSPSetOperators(KSP_MyMR_Neu, MyMR_Neu, MyMR_Neu, SAME_NONZERO_PATTERN, iErr)
      Call KSPSetType(KSP_MyMR_Neu, KSPPREONLY, iErr)
!!     Call KSPSetType(KSP_MyMR_Neu, KSPCG, iErr)
      Call KSPGetPC(KSP_MyMR_Neu, PC_MyMR_Neu, iErr)
!!     Call PCSetType(PC_MyMR_Neu, PCJACOBI, iErr)
       Call PCSetType(PC_MyMR_Neu, PCLU, iErr)
!!*      Call PCSetType(PC_MyMR_Neu, PCILU, iErr)
      If (count(Node_db(MyNode_Ovlp+1)%BC /= BC_Type_NONE) == 0) Then
! Null Space Treatment 
! Construct orthonomal basis for Nullspace (1,1,1,1,...,1)^t
       Call VecCreateSeq(PETSC_COMM_SELF, MySize_Node_Ovlp, Orth, iErr) 
        Call VecSet(Orth,one,iErr)
        Call VecDot(Orth,Orth, alpha, iErr)
        alpha = 1.0/Sqrt(alpha)
        Call VecScale(Orth, alpha, iErr) 
        Call MatNullSpaceCreate(PETSC_COMM_SELF, PETSC_False, 1, Orth, NeuNullSpace, iErr)
        Call KSPSetNullSpace(KSP_MyMR_Neu, NeuNullSpace, iErr)
!!!        Call PCLUSetDamping(PC_MyMR_Neu,0.0000000001, iErr)
        Call PCFactorSetShiftNonzero(PC_MyMR_Neu, 1.0e-10, iErr)
!!       Call MatNullSpaceTest(NeuNullSpace, MyMR_Neu, iErr)
!       Call SLESGetKSP(SLES_M, KSP_M, iErr)
!       Call SLESGetPC(SLES_M, PC_M, iErr)
!       Call PCSetType(PC_M, PCJACOBI, iErr)
!       Call KSPSetType(KSP_M, KSPGMRES, iErr)
!       Call PCNullSpaceAttach(PC_M, nullSpace, iErr)
      End If
      Call KSPSetFromOptions(KSP_MyMR_Neu, iErr)
     Tol = 0.01
     Call KSPSetTolerances(KSP_MyMR_Neu, Tol,                               &
          & PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION,  &
         & PETSC_DEFAULT_INTEGER, iErr)
!          & 2, iErr)
    End If
    Call KSPSetFromOptions(KSP_MR, iErr)
    Call KSPSetup(KSP_MR, iErr)
  End Subroutine Init_KSP


!! Compute L2 norm

  Subroutine CompL2(D_Vec, Geom, Elem_db, Node_db, Value)
    Real(Kind = Kr), Dimension(:), Pointer              :: D_Vec
    Type (EXO_Geom_Info)                                :: Geom
#ifdef PB_2D
    Type (Node2D), Dimension(:), Pointer                :: Node_db 
    Type (Element2D_Scal), Dimension(:), Pointer        :: Elem_db 
#else
    Type (Node3D), Dimension(:), Pointer                :: Node_db 
    Type (Element3D_Scal), Dimension(:), Pointer        :: Elem_db 
#endif
    Integer(Kind = Ki), Dimension(:), Pointer           :: iSGSigIn
    Integer(Kind = Ki)                                  :: Nb_DoF
    Integer(Kind = Ki)                                  :: iSLSig, iSGSig
    Integer(Kind = Ki)                                  :: iSLEps, iSGEps
    Integer(Kind = Ki)                                  :: iE, iELoc
    Integer(Kind = Ki)                                  :: iBlk
    Integer(Kind = Ki)                                  :: iG, Nb_Gauss
    Real(Kind = Kr)                                     :: Value

       Value = 0.0_Kr
       Do_iBlk: Do iBlk = 1, Geom%Num_Elem_Blks
       Nb_DoF = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem 
       Do_iE: Do iELoc = 1, Geom%Num_Elems
          If (Elem_Owner(iELoc)/=MyRank) Then
             CYCLE
          End If
          iE = Geom%Elem_blk(iBlk)%Elem_ID(iELoc)
          Call Init_Gauss_EXO(Elem_db, Node_db, Geom, GaussOrder,    &
               & Elem=iE)
          Nb_Gauss = Elem_db(iE)%Nb_Gauss
          Do_iSLSig: Do iSLSig = 1, Nb_DoF
            ISGSig = Elem_db(iE)%ID_DoF(iSLSig)
            Do iSLEps = 1, Nb_DoF
              iSGEps = Elem_db(iE)%ID_DoF(iSLEps)
              Do iG = 1, Nb_Gauss
                  Value = Value+Elem_db(iE)%Gauss_C(iG) *  &
                     & Elem_db(iE)%BF(iSLSig,iG) * Elem_db(iE)%BF(iSLEps,iG)*D_Vec(iSGSig)*D_Vec(iSGEps)
               End Do
             End Do
           End Do Do_iSLSig
          Call Destroy_Gauss_EXO(Elem_db, Elem=iE)
       EndDo Do_iE
    End Do Do_iBlk
  End Subroutine CompL2

!! Compute SemiH1
  Subroutine CompSEMIH2(D_Vec, Geom, Elem_db, Node_db, Value)
    Real(Kind = Kr), Dimension(:), Pointer              :: D_Vec
    Type (EXO_Geom_Info)                                :: Geom
#ifdef PB_2D
    Type (Node2D), Dimension(:), Pointer                :: Node_db 
    Type (Element2D_Scal), Dimension(:), Pointer        :: Elem_db 
#else
    Type (Node3D), Dimension(:), Pointer                :: Node_db 
    Type (Element3D_Scal), Dimension(:), Pointer        :: Elem_db 
#endif
    Integer(Kind = Ki), Dimension(:), Pointer           :: iSGSigIn
    Integer(Kind = Ki)          :: Nb_DoF
    Integer(Kind = Ki)          :: iSLSig, iSGSig
    Integer(Kind = Ki)          :: iSLEps, iSGEps
    Integer(Kind = Ki)          :: iE, iELoc
    Integer(Kind = Ki)          :: iBlk
    Integer(Kind = Ki)          :: iG, Nb_Gauss
    Real(Kind = Kr)             :: Value

       Value = 0.0_Kr
       Do_iBlk: Do iBlk = 1, Geom%Num_Elem_Blks
       Nb_DoF = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem 
       Do_iE: Do iELoc = 1, Geom%Num_Elems
          If (Elem_Owner(iELoc)/=MyRank) Then
             CYCLE
          End If
          iE = Geom%Elem_blk(iBlk)%Elem_ID(iELoc)
          Call Init_Gauss_EXO(Elem_db, Node_db, Geom, GaussOrder,    &
               & Elem=iE)
          Nb_Gauss = Elem_db(iE)%Nb_Gauss
          Do_iSLSig: Do iSLSig = 1, Nb_DoF
            ISGSig = Elem_db(iE)%ID_DoF(iSLSig)
            Do iSLEps = 1, Nb_DoF
              iSGEps = Elem_db(iE)%ID_DoF(iSLEps)
              Do iG = 1, Nb_Gauss
                  Value = Value+ &
!                    &  Elem_db(iE)%Gauss_C(iG) *  &
!                     & Elem_db(iE)%BF(iSLSig,iG) * Elem_db(iE)%BF(iSLEps,iG)*D_Vec(iSGSig)*D_Vec(iSGEps) &
                     & + Elem_db(iE)%Gauss_C(iG) *  ( (D_Vec(iSGSig)*Elem_db(iE)%Grad_BF(iSLEps,iG)) .DotP.      &
                     &  (D_Vec(iSGEps)*Elem_db(iE)%Grad_BF(iSLSig,iG)) )
               End Do
             End Do
           End Do Do_iSLSig
          Call Destroy_Gauss_EXO(Elem_db, Elem=iE)
       EndDo Do_iE
    End Do Do_iBlk
  End Subroutine CompSEMIH2

!! Export

  Subroutine Export()
    Call VecScatterBegin(SOL, SOL_Master, INSERT_VALUES, SCATTER_FORWARD,     &
         & SOL_DistToMaster, iErr)
    Call VecScatterEnd(SOL, SOL_Master, INSERT_VALUES, SCATTER_FORWARD,       &
         & SOL_DistToMaster, iErr)
    Call VecGetArrayF90(SOL_Master, Disp_Sol, iErr)
    Allocate (Disp_L2(Geom%Num_Nodes))
    If (MyRank == 0) Then
      Write(*,90) MinVal(Disp_Sol), MaxVal(Disp_Sol)
      Disp_L2 = Disp_Sol
    End If
    Call VecRestoreArrayF90(SOL_Master, Disp_Sol, iErr)
    Call MPI_BCAST(Disp_L2, Geom%Num_Nodes, MPI_DOUBLE_PRECISION, 0,PETSC_COMM_WORLD, iErr)
90  Format('Solution Min / max: ', 2(ES10.3,' '))
  End Subroutine Export
  
  Subroutine SaveLayout()
    Real(Kind = Kr), Dimension(:), Pointer           :: Tmp_Ptr

    If (MyRank == 0) Then
       Allocate(Tmp_Ptr(Geom%Num_Nodes))
       Tmp_Ptr = real(Node_Owner)
       Call Write_EXO_Result_Ptr_Nodes(Geom, 3, 1, Tmp_Ptr)
       DeAllocate(Tmp_Ptr)
       
       Allocate(Tmp_Ptr(Geom%Num_Elems))
       Tmp_Ptr = real(Elem_Owner)
       Call Write_EXO_Result_Scal_Elems(Geom, 1, 1, Tmp_Ptr)
       DeAllocate(Tmp_Ptr)
    End If
  End Subroutine SaveLayout

  Subroutine Finalize()
    DeAllocate (Load)
    DeAllocate(Node_db, Elem_db)
    Call MatDestroy(MR, iErr)
    Call VecDestroy(RHS, iErr)
    Call PETScFinalize(iErr)
  End Subroutine Finalize

  Subroutine Init_BC(Geom, Params, Node_db)
    Type(EXO_Geom_Info), Intent(IN)                  :: Geom
    Type(Poisson_Params), Intent(IN)                 :: Params

#ifdef PB_2D
    Type(Node2D), Dimension(:), Pointer              :: Node_db
#else
    Type(Node3D), Dimension(:), Pointer              :: Node_db
#endif

    Integer                                          :: iN, iSet
    Node_db(:)%BC = BC_Type_NONE
    Do iSet = 1, Geom%Num_Node_Sets
     Do iN = 1, Geom%Node_Set(iSet)%Num_Nodes
      Node_db(Geom%Node_Set(iSet)%Node_ID(iN))%BC = Params%BC_Type(iSet)
     End Do
    End Do
  End Subroutine Init_BC

  Subroutine Assemb_Mat_Poisson(MR, Geom, Params, Elem_db, Node_db)
    Mat                                                 :: MR
    Type (EXO_Geom_Info)                                :: Geom
    Type (Poisson_Params)                               :: Params
#ifdef PB_2D
    Type (Node2D), Dimension(:), Pointer                :: Node_db 
    Type (Element2D_Scal), Dimension(:), Pointer        :: Elem_db 
#else
    Type (Node3D), Dimension(:), Pointer                :: Node_db 
    Type (Element3D_Scal), Dimension(:), Pointer        :: Elem_db 
#endif
    Integer(Kind = Ki)          :: Nb_Gauss, Nb_DoF
    Integer          :: iSLEps, iSLSig
    Integer          :: iSGEps, iSGSig
    Integer          :: iSGSigBIn 
    Integer, Dimension(:), Pointer :: iSGSigIn
    Integer(Kind = Ki)          :: iE, iG, iELoc
    Integer(Kind = Ki)          :: iBlk
    PetscScalar                 :: one
    PetscScalar, Dimension(:,:), Pointer              :: MR_Elem
    PetscLogDouble              :: TotTS, TotTF, TotT
    PetscLogDouble              :: SetTS, SetTF, SetT
    PetscLogDouble              :: AssembTS, AssembTF, AssembT
    PetscLogDouble              :: GaussTS, GaussTF, GaussT
    Integer                     :: SetN, SetNTot



    SetN = 0
    SetT = 0.0
    AssembT = 0.0
    GaussT = 0.0

    one = 1.0_Kr
    Call PetscGetTime(TotTS, iErr)
!!    Call MatZeroEntries(MR, iErr)

    Do_iBlk: Do iBlk = 1, Geom%Num_Elem_Blks
       Nb_DoF = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem 
       Allocate (MR_Elem(Nb_DoF, Nb_Dof))


       Do_iE: Do iELoc = 1, Geom%Num_Elems
          iE = Geom%Elem_blk(iBlk)%Elem_ID(iELoc)
          If (Elem_Owner(iE)/=MyRank) Then
            CYCLE
          End If

          Call PetscGetTime(GaussTS, iErr)
          Call Init_Gauss_EXO(Elem_db, Node_db, Geom, GaussOrder,    &
               & Elem=iE)
          Call PetscGetTime(GaussTF, iErr)
          GaussT = GaussT + GaussTF - GaussTS

          MR_Elem = 0.0_Kr
          Nb_Gauss = Elem_db(iE)%Nb_Gauss
          Allocate(ISGSigIn(Nb_DoF))
!          Allocate(ISGEpsIn(Nb_DoF))
          iSGSigIn = Elem_db(iE)%ID_DoF-1

          Do_iSLSig: Do iSLSig = 1, Nb_DoF

             iSGSig = Elem_db(iE)%ID_DoF(iSLSig)
             

!                   iSGSigIn(iSLSig) = iSGSig-1
!             Write(200+MyRank, *) 'Sig', iE, iSLSig, iSGSig, iSGSigIn(iSLSig)
!!nn               iSGSigIn(iSLSig) = All_list(iSGSig-1)
            Is_BCSig: If (Node_db(iSGSig)%BC /= BC_Type_NONE) Then
               CYCLE
             End If Is_BCSig
 !!x               If  (Any(iSGSigIn > Size(All_list))) Then
 !!x                  Print*, 'Proc:', MyRank, 'iSGSigIn', iSGSigIn
 !!x               End If
             Do_iSLEps: Do iSLEps = 1, Nb_DoF
                iSGEps = Elem_db(iE)%ID_DoF(iSLEps)
                
!                  iSGEpsIn(iSLEps) = iSGEps-1
                      
!!nn                iSGEpsIn(iSLEps) = All_list(iSGEps-1)
!             Write(200+MyRank, *) 'Eps', iE, iSLEps, iSGEps, iSGEpsIn(iSLEps)
                Is_BCEps: If (Node_db(iSGEps)%BC /= BC_Type_NONE) Then
                   CYCLE
                End If Is_BCEps
!!x                If (Any(iSGEpsIn > Size(All_list))) Then
!!x                   Print*, 'Proc:', MyRank, 'iSGEpsIn', iSGEpsIn
!!x                End If
                Do iG = 1, Nb_Gauss
                   MR_Elem(iSLSig, iSLEps) = MR_Elem(iSLSig, iSLEps)    &
                              & + Elem_db(iE)%Gauss_C(iG) *                   &
                              & ( Elem_db(iE)%Grad_BF(iSLEps,iG) .DotP.      &
                              &  Elem_db(iE)%Grad_BF(iSLSig,iG))
                End Do

             End Do Do_iSLEps
          End Do Do_iSLSig

          Call PetscGetTime(SetTS, iErr)
!               If  (Any(iSGSigIn > Size(All_list))) Then
!                   Print*, 'BF Proc:', MyRank, 'iSGSigIn', iSGSigIn
!                End If
!               If (Any(iSGEpsIn > Size(All_list))) Then
!                   Print*, 'BF Proc:', MyRank, 'iSGEpsIn', iSGEpsIn
!                   
!                End If
!          Call MatSetValues(MR, Nb_DoF, All_list(ISGSigIn), Nb_DoF,       &
!               & All_list(ISGEpsIn), MR_Elem, ADD_VALUES, iErr)

          Call MatSetValues(MR, Nb_DoF, All_list(ISGSigIn), Nb_DoF,       &
               & All_list(ISGSigIn), MR_Elem, ADD_VALUES, iErr)

!!nn          Call MatSetValues(MR, Nb_DoF, ISGSigIn, Nb_DoF,       &
!!nn               & ISGEpsIn, MR_Elem, ADD_VALUES, iErr)

          Call PetscGetTime(SetTF, iErr)
          SetN = SetN + 1
          SetT = SetT + SetTF - SetTS
          Call PetscGetTime(GaussTS, iErr)
          Call Destroy_Gauss_EXO(Elem_db, Elem=iE)
          Call PetscGetTime(GaussTF, iErr)
          GaussT = Gausst + GaussTF - GaussTS
          DeAllocate(ISGSigIn)
 !         DeAllocate(ISGEpsIn)
       EndDo Do_iE
       DeAllocate (MR_Elem)
    End Do Do_iBlk
    Call PetscGetTime(AssembTS, iErr)
    Call MatAssemblyBegin(MR, MAT_FINAL_ASSEMBLY, iErr)
    Call MatAssemblyEnd(MR, MAT_FINAL_ASSEMBLY, iErr)  
    Call PetscGetTime(AssembTF, iErr)
    AssembT = AssembT + AssembTF - AssembTS

    ! Assembly of the BC terms 
    Do_iBlk_BC: Do iBlk = 1, Geom%Num_Elem_Blks
       Do_iE_BC: Do iELoc = 1, Geom%Num_Elems
          iE = Geom%Elem_blk(iBlk)%Elem_ID(iELoc)
          If (Elem_Owner(iE)/=MyRank) Then
            CYCLE
          End If
          Do_iSLSig_BC: Do iSLSig = 1, Nb_DoF
             iSGSig = Elem_db(iE)%ID_DoF(iSLSig)
             Is_BC_BC: If (Node_db(iSGSig)%BC /= BC_Type_NONE) Then
                    iSGSigBIn= iSGSig-1
!!nn                 iSGSigBIn= All_list(iSGSig-1)
                Call PetscGetTime(SetTS, iErr)
              If  (iSGSigBIn > Size(All_list)) Then
                   Print*, 'BFB Proc:', MyRank, 'iSGSigBIn', iSGSigBIn
                End If
                Call MatSetValues(MR, 1, All_list(ISGSigBIn), 1,              &
                     & All_list(ISGSigBIn), one, INSERT_VALUES, iErr)
!!nn                Call MatSetValues(MR, 1, ISGSigBIn, 1,              &
!!nn                     & ISGSigBIn, one, INSERT_VALUES, iErr)
                Call PetscGetTime(SetTF, iErr)
                SetN = SetN+1
                SetT = SetT + SetTF - SetTS
             End If Is_BC_BC
          End Do Do_iSLSig_BC
       End Do Do_iE_BC
    End Do Do_iBlk_BC 

    Call PetscGetTime(AssembTS, iErr)
    Call MatAssemblyBegin(MR, MAT_FINAL_ASSEMBLY, iErr)
    Call MatAssemblyEnd(MR, MAT_FINAL_ASSEMBLY, iErr)
    Call PetscGetTime(AssembTF, iErr)
    AssembT = AssembT + AssembTF - AssembTS
    Call PetscGetTime(TotTF, iErr)
    TotT = TotTF - TotTS

    Call MPI_REDUCE(SetN, SetNTot, 1, MPI_INTEGER, MPI_SUM, 0,                &
         & PETSC_COMM_WORLD, iErr)

    If (MyRank ==0) Then
       Write(CharBuffer,*) 'Total time in Init/Destroy_Gauss:        ',       &
            & GaussT, '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
       Write(CharBuffer,*) 'Total time in MatSetValue:               ', SetT, &
            & '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
       Write(CharBuffer,*) 'Number of calls to MatSetValue:   ', SetNTot, '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
       Write(CharBuffer,*) 'Total time in MatAssembly:               ',       &
            & AssembT, '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
       Write(CharBuffer,*) 'Total time in Assemb_Mat_Poisson:        ',       &
            & TotT, '\n\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
    End If
  End Subroutine Assemb_Mat_Poisson





  Subroutine Assemb_Mat_Local_Poisson(MyMR, Geom, Params, Elem_db, Node_db)
    Mat                                                 :: MyMR
    Type (EXO_Geom_Info)                                :: Geom
    Type (Poisson_Params)                               :: Params
#ifdef PB_2D
    Type (Node2D), Dimension(:), Pointer                :: Node_db 
    Type (Element2D_Scal), Dimension(:), Pointer        :: Elem_db 
#else
    Type (Node3D), Dimension(:), Pointer                :: Node_db 
    Type (Element3D_Scal), Dimension(:), Pointer        :: Elem_db 
#endif
    Integer(Kind = Ki)          :: Nb_Gauss, Nb_DoF
    Integer(Kind = Ki)          :: iSLEps, iSLSig
    Integer(Kind = Ki)          :: iSGEps, iSGSig
    Integer(Kind = Ki)          :: iSGSigBIn 
    Integer, Dimension(:), Pointer :: iSGEpsIn, iSGSigIn
    Integer(Kind = Ki)          :: iE, iG, iELoc
    Integer(Kind = Ki)          :: iBlk
    PetscScalar                 :: one
    PetscScalar, Dimension(:,:), Pointer              :: MR_Elem
    PetscLogDouble              :: TotTS, TotTF, TotT
    PetscLogDouble              :: SetTS, SetTF, SetT
    PetscLogDouble              :: AssembTS, AssembTF, AssembT
    PetscLogDouble              :: GaussTS, GaussTF, GaussT
    Integer                     :: SetN, SetNTot



    SetN = 0
    SetT = 0.0
    AssembT = 0.0
    GaussT = 0.0

    one = 1.0_Kr
    Call PetscGetTime(TotTS, iErr)
!!    Call MatZeroEntries(MyMR, iErr)

    Do_iBlk: Do iBlk = 1, Geom%Num_Elem_Blks
       Nb_DoF = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem 
       Allocate (MR_Elem(Nb_DoF, Nb_Dof))
       Do_iE: Do iELoc = 1, MySize_Elem_Ovlp
          iE = MyElem_Ovlp(iELoc)
          Call PetscGetTime(GaussTS, iErr)
          Call Init_Gauss_EXO(Elem_db, Node_db, Geom, GaussOrder,    &
               & Elem=iE)
          Call PetscGetTime(GaussTF, iErr)
          GaussT = GaussT + GaussTF - GaussTS

          MR_Elem = 0.0_Kr
          Nb_Gauss = Elem_db(iE)%Nb_Gauss
          Allocate(ISGSigIn(Nb_DoF))
          Allocate(ISGEpsIn(Nb_DoF))
          Do_iSLSig: Do iSLSig = 1, Nb_DoF
             ISGSig = Elem_db(iE)%ID_DoF(iSLSig)
               ISGSigIn(iSLSig) = ISGSig-1
!             Is_BCSig: If ((Node_db(iSGSig)%BC /= BC_Type_NONE) .or. (MyCount_ID_Node_Bd(iSGSig)==1)) Then
!               CYCLE
!             End If Is_BCSig
             Do_iSLEps: Do iSLEps = 1, Nb_DoF
                iSGEps = Elem_db(iE)%ID_DoF(iSLEps)
                ISGEpsIn(iSLEps) = ISGEps-1
                Is_BCEps: If ((Node_db(iSGEps)%BC /= BC_Type_NONE) .or. (MyCount_ID_Node_Bd(iSGEps)==1)) Then
                   CYCLE
                End If Is_BCEps
                Is_BCSig: If ((Node_db(iSGSig)%BC /= BC_Type_NONE) .or. (MyCount_ID_Node_Bd(iSGSig)==1)) Then
                   CYCLE
                End If Is_BCSig
                Do iG = 1, Nb_Gauss
                   MR_Elem(iSLSig, iSLEps) = MR_Elem(iSLSig, iSLEps)    &
                              & + Elem_db(iE)%Gauss_C(iG) *                   &
                              & ( Elem_db(iE)%Grad_BF(iSLEps,iG) .DotP.      &
                              &  Elem_db(iE)%Grad_BF(iSLSig,iG))
                End Do
             End Do Do_iSLEps
          End Do Do_iSLSig
          Call PetscGetTime(SetTS, iErr)
          Call MatSetValues(MyMR, Nb_DoF, MyNode_Ovlp_list(ISGSigIn), Nb_DoF,       &
               &  MyNode_Ovlp_list(ISGEpsIn), MR_Elem, ADD_VALUES, iErr)
          Call PetscGetTime(SetTF, iErr)
          SetN = SetN + 1
          SetT = SetT + SetTF - SetTS
          Call PetscGetTime(GaussTS, iErr)
          Call Destroy_Gauss_EXO(Elem_db, Elem=iE)
          Call PetscGetTime(GaussTF, iErr)
          GaussT = Gausst + GaussTF - GaussTS
          DeAllocate(ISGSigIn)
          DeAllocate(ISGEpsIn)
       EndDo Do_iE
       DeAllocate (MR_Elem)
    End Do Do_iBlk

    Call PetscGetTime(AssembTS, iErr)
    Call MatAssemblyBegin(MyMR, MAT_FINAL_ASSEMBLY, iErr)
    Call MatAssemblyEnd(MyMR, MAT_FINAL_ASSEMBLY, iErr)  
    Call PetscGetTime(AssembTF, iErr)
    AssembT = AssembT + AssembTF - AssembTS

    ! Assembly of the BC terms 
    Do_iBlk_BC: Do iBlk = 1, Geom%Num_Elem_Blks
       Do_iE_BC: Do iELoc = 1, MySize_Elem_Ovlp
          iE = MyElem_Ovlp(iELoc)
!          If (Elem_Owner(iE)/=MyRank) Then
!            CYCLE
!          End If
          Do_iSLSig_BC: Do iSLSig = 1, Nb_DoF
             ISGSig = Elem_db(iE)%ID_DoF(iSLSig)
             Is_BC_BC: If ((Node_db(iSGSig)%BC /= BC_Type_NONE) .or. (MyCount_ID_Node_Bd(iSGSig)==1)) Then
                 ISGSigBIn= ISGSig-1
                Call PetscGetTime(SetTS, iErr)
                Call MatSetValues(MyMR, 1, MyNode_Ovlp_list(ISGSigBIn), 1,              &
                     &  MyNode_Ovlp_list(ISGSigBIn), one, INSERT_VALUES, iErr)
                Call PetscGetTime(SetTF, iErr)
                SetN = SetN+1
                SetT = SetT + SetTF - SetTS
             End If Is_BC_BC
          End Do Do_iSLSig_BC
       End Do Do_iE_BC
    End Do Do_iBlk_BC 

    Call PetscGetTime(AssembTS, iErr)
    Call MatAssemblyBegin(MyMR, MAT_FINAL_ASSEMBLY, iErr)
    Call MatAssemblyEnd(MyMR, MAT_FINAL_ASSEMBLY, iErr)
    Call PetscGetTime(AssembTF, iErr)
    AssembT = AssembT + AssembTF - AssembTS
    Call PetscGetTime(TotTF, iErr)
    TotT = TotTF - TotTS

    Call MPI_REDUCE(SetN, SetNTot, 1, MPI_INTEGER, MPI_SUM, 0,                &
         & PETSC_COMM_WORLD, iErr)

    If (MyRank ==0) Then
       Write(CharBuffer,*) 'Total time in Init/Destroy_Gauss:        ',       &
            & GaussT, '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
       Write(CharBuffer,*) 'Total time in MatSetValue:               ', SetT, &
            & '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
       Write(CharBuffer,*) 'Number of calls to MatSetValue:   ', SetNTot, '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
       Write(CharBuffer,*) 'Total time in MatAssembly:               ',       &
            & AssembT, '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
       Write(CharBuffer,*) 'Total time in Assemb_Mat_Poisson:        ',       &
            & TotT, '\n\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
    End If
  End Subroutine Assemb_Mat_Local_Poisson


  Subroutine Assemb_Mat_Local_Poisson_Neu(MyMR_Neu, Geom, Params, Elem_db, Node_db)
    Mat                                                 :: MyMR_Neu
    Type (EXO_Geom_Info)                                :: Geom
    Type (Poisson_Params)                               :: Params
#ifdef PB_2D
    Type (Node2D), Dimension(:), Pointer                :: Node_db 
    Type (Element2D_Scal), Dimension(:), Pointer        :: Elem_db 
#else
    Type (Node3D), Dimension(:), Pointer                :: Node_db 
    Type (Element3D_Scal), Dimension(:), Pointer        :: Elem_db 
#endif
    Integer(Kind = Ki)          :: Nb_Gauss, Nb_DoF
    Integer(Kind = Ki)          :: iSLEps, iSLSig
    Integer(Kind = Ki)          :: iSGEps, iSGSig
    Integer(Kind = Ki)          :: iSGSigBIn 
    Integer, Dimension(:), Pointer :: iSGEpsIn, iSGSigIn
    Integer(Kind = Ki)          :: iE, iG, iELoc
    Integer(Kind = Ki)          :: iBlk
    PetscScalar                 :: one
    PetscScalar, Dimension(:,:), Pointer              :: MR_Elem
    PetscLogDouble              :: TotTS, TotTF, TotT
    PetscLogDouble              :: SetTS, SetTF, SetT
    PetscLogDouble              :: AssembTS, AssembTF, AssembT
    PetscLogDouble              :: GaussTS, GaussTF, GaussT
    Integer                     :: SetN, SetNTot



    SetN = 0
    SetT = 0.0
    AssembT = 0.0
    GaussT = 0.0

    one = 1.0_Kr
    Call PetscGetTime(TotTS, iErr)
!!    Call MatZeroEntries(MyMR_Neu, iErr)

    Do_iBlk: Do iBlk = 1, Geom%Num_Elem_Blks
       Nb_DoF = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem 
       Allocate (MR_Elem(Nb_DoF, Nb_Dof))
       Do_iE: Do iELoc = 1, MySize_Elem_Ovlp
          iE = MyElem_Ovlp(iELoc)
          Call PetscGetTime(GaussTS, iErr)
          Call Init_Gauss_EXO(Elem_db, Node_db, Geom, GaussOrder,    &
               & Elem=iE)
          Call PetscGetTime(GaussTF, iErr)
          GaussT = GaussT + GaussTF - GaussTS

          MR_Elem = 0.0_Kr
          Nb_Gauss = Elem_db(iE)%Nb_Gauss
          Allocate(ISGSigIn(Nb_DoF))
          Allocate(ISGEpsIn(Nb_DoF))
          Do_iSLSig: Do iSLSig = 1, Nb_DoF
             ISGSig = Elem_db(iE)%ID_DoF(iSLSig)
               ISGSigIn(iSLSig) = ISGSig-1
!             Is_BCSig: If ((Node_db(iSGSig)%BC /= BC_Type_NONE) .or. (MyCount_ID_Node_Bd(iSGSig)==1)) Then
!               CYCLE
!             End If Is_BCSig
             Do_iSLEps: Do iSLEps = 1, Nb_DoF
                iSGEps = Elem_db(iE)%ID_DoF(iSLEps)
                ISGEpsIn(iSLEps) = ISGEps-1
                 Is_BCEps: If (Node_db(iSGEps)%BC /= BC_Type_NONE) Then          
                   CYCLE
                End If Is_BCEps
                Is_BCSig: If (Node_db(iSGSig)%BC /= BC_Type_NONE)  Then
                   CYCLE
                End If Is_BCSig
                Do iG = 1, Nb_Gauss
                   MR_Elem(iSLSig, iSLEps) = MR_Elem(iSLSig, iSLEps)    &
                              & + Elem_db(iE)%Gauss_C(iG) *                   &
                              & ( Elem_db(iE)%Grad_BF(iSLEps,iG) .DotP.      &
                              &  Elem_db(iE)%Grad_BF(iSLSig,iG))
                End Do
             End Do Do_iSLEps
          End Do Do_iSLSig
          Call PetscGetTime(SetTS, iErr)
          Call MatSetValues(MyMR_Neu, Nb_DoF, MyNode_Ovlp_list(ISGSigIn), Nb_DoF,       &
               &  MyNode_Ovlp_list(ISGEpsIn), MR_Elem, ADD_VALUES, iErr)
          Call PetscGetTime(SetTF, iErr)
          SetN = SetN + 1
          SetT = SetT + SetTF - SetTS
          Call PetscGetTime(GaussTS, iErr)
          Call Destroy_Gauss_EXO(Elem_db, Elem=iE)
          Call PetscGetTime(GaussTF, iErr)
          GaussT = Gausst + GaussTF - GaussTS
          DeAllocate(ISGSigIn)
          DeAllocate(ISGEpsIn)
       EndDo Do_iE
       DeAllocate (MR_Elem)
    End Do Do_iBlk

    Call PetscGetTime(AssembTS, iErr)
    Call MatAssemblyBegin(MyMR_Neu, MAT_FINAL_ASSEMBLY, iErr)
    Call MatAssemblyEnd(MyMR_Neu, MAT_FINAL_ASSEMBLY, iErr)  
    Call PetscGetTime(AssembTF, iErr)
    AssembT = AssembT + AssembTF - AssembTS

    ! Assembly of the BC terms 
    Do_iBlk_BC: Do iBlk = 1, Geom%Num_Elem_Blks
       Do_iE_BC: Do iELoc = 1, MySize_Elem_Ovlp
          iE = MyElem_Ovlp(iELoc)
!          If (Elem_Owner(iE)/=MyRank) Then
!            CYCLE
!          End If
          Do_iSLSig_BC: Do iSLSig = 1, Nb_DoF
             ISGSig = Elem_db(iE)%ID_DoF(iSLSig)
             Is_BC_BC: If (Node_db(iSGSig)%BC /= BC_Type_NONE) Then
                 ISGSigBIn= ISGSig-1
                Call PetscGetTime(SetTS, iErr)
                Call MatSetValues(MyMR_Neu, 1, MyNode_Ovlp_list(ISGSigBIn), 1,              &
                     &  MyNode_Ovlp_list(ISGSigBIn), one, INSERT_VALUES, iErr)
                Call PetscGetTime(SetTF, iErr)
                SetN = SetN+1
                SetT = SetT + SetTF - SetTS
             End If Is_BC_BC
          End Do Do_iSLSig_BC
       End Do Do_iE_BC
    End Do Do_iBlk_BC 

    Call PetscGetTime(AssembTS, iErr)
    Call MatAssemblyBegin(MyMR_Neu, MAT_FINAL_ASSEMBLY, iErr)
    Call MatAssemblyEnd(MyMR_Neu, MAT_FINAL_ASSEMBLY, iErr)
    Call PetscGetTime(AssembTF, iErr)
    AssembT = AssembT + AssembTF - AssembTS
    Call PetscGetTime(TotTF, iErr)
    TotT = TotTF - TotTS

    Call MPI_REDUCE(SetN, SetNTot, 1, MPI_INTEGER, MPI_SUM, 0,                &
         & PETSC_COMM_WORLD, iErr)

    If (MyRank ==0) Then
       Write(CharBuffer,*) 'Total time in Init/Destroy_Gauss:        ',       &
            & GaussT, '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
       Write(CharBuffer,*) 'Total time in MatSetValue:               ', SetT, &
            & '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
       Write(CharBuffer,*) 'Number of calls to MatSetValue:   ', SetNTot, '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
       Write(CharBuffer,*) 'Total time in MatAssembly:               ',       &
            & AssembT, '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
       Write(CharBuffer,*) 'Total time in Assemb_Mat_Poisson_Neu:        ',       &
            & TotT, '\n\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
    End If
  End Subroutine Assemb_Mat_Local_Poisson_Neu


  Subroutine Assemb_RHS_Poisson(RHS, Geom, Params, Elem_db, Node_db, Load)

    Vec                                                 :: RHS
    Type (EXO_Geom_Info)                                :: Geom
    Type (Poisson_Params)                               :: Params
#ifdef PB_2D
    Type(Element2D_Scal), Dimension(:), Pointer         :: Elem_db
    Type(Node2D), Dimension(:), Pointer                 :: Node_db
#else
    Type(Element3D_Scal), Dimension(:), Pointer         :: Elem_db
    Type(Node3D), Dimension(:), Pointer                 :: Node_db
#endif
    Real(Kind = Kr), Dimension(:), Pointer              :: Load
    PetscScalar, Dimension(:), Pointer                  :: RHS_Ptr
    Integer(Kind = Ki), Dimension(:), Pointer           :: iSGSigIn
    Integer(Kind = Ki)          :: Nb_DoF
    Integer(Kind = Ki)          :: iSLSig, iSGSig
    Integer(Kind = Ki)          :: iE, iELoc
    Integer(Kind = Ki)          :: iBlk
    Integer(Kind = Ki)          :: iG, Nb_Gauss

    PetscLogDouble              :: TotTS, TotTF, TotT
    PetscLogDouble              :: SetTS, SetTF, SetT
    PetscLogDouble              :: AssembTS, AssembTF, AssembT
    PetscLogDouble              :: GaussTS, GaussTF, GaussT
    Integer                     :: SetN, SetNTot


    SetN = 0
    SetT = 0.0
    AssembT = 0.0
    GaussT = 0.0

    Call PetscGetTime(TotTS, iErr)

    Do_iBlk: Do iBlk = 1, Geom%Num_Elem_Blks
       Nb_DoF = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem 
       Allocate(RHS_Ptr(Nb_DoF))
       Do_iE: Do iELoc = 1, Geom%Num_Elems
          iE = Geom%Elem_blk(iBlk)%Elem_ID(iELoc)
          If (Elem_Owner(iE)/=MyRank) Then
             CYCLE
          End If
          Call PetscGetTime(GaussTS, iErr)
          Call Init_Gauss_EXO(Elem_db, Node_db, Geom, GaussOrder, Elem=iE)
          Call PetscGetTime(GaussTF, iErr)
          GaussT = GaussT + GaussTF - GaussTS

          RHS_Ptr = 0.0_Kr
          Nb_Gauss = Elem_db(iE)%Nb_Gauss
          Allocate(ISGSigIn(Nb_DoF))
          
          Do_iSLSig: Do iSLSig = 1, Nb_DoF
             ISGSig = Elem_db(iE)%ID_DoF(iSLSig)
             ISGSigIn(iSLSig)=ISGSig-1
             Is_BCSig: If (Node_db(iSGSig)%BC /= BC_Type_NONE) Then
                Cycle
             End If Is_BCSig
             Do iG = 1, Nb_Gauss
#ifdef PB_2D
                RHS_Ptr(iSLSig) = RHS_Ptr(iSLSig)+Elem_db(iE)%Gauss_C(iG) *  &
                     & Elem_db(iE)%BF(iSLSig,iG) * Load(iSGSig) 

!!* Pi **2 /   & 
!!                     & 50.0_Kr
#else
                RHS_Ptr(iSLSig) = RHS_Ptr(iSLSig)+Elem_db(iE)%Gauss_C(iG) *  &
                     & Elem_db(iE)%BF(iSLSig,iG) * Load(iSGSig) 
!!* 3.0_Kr *   &
!!                     & Pi**2  / 100.0_Kr
#endif
             End Do
          End Do Do_iSLSig
          Call PetscGetTime(SetTS, iErr)
          Call VecSetValues(RHS, Nb_DoF, All_list(ISGSigIn), RHS_Ptr,         &
               & ADD_VALUES,iErr) 
          Call PetscGetTime(SetTF, iErr)
          SetN = SetN + 1
          SetT = SetT + SetTF - SetTS
          Call PetscGetTime(GaussTS, iErr)
          Call Destroy_Gauss_EXO(Elem_db, Elem=iE)
          Call PetscGetTime(GaussTF, iErr)
          GaussT = Gausst + GaussTF - GaussTS
          DeAllocate(ISGSigIn)
       EndDo Do_iE
       DeAllocate(RHS_Ptr)
    End Do Do_iBlk

    Call PetscGetTime(AssembTS, iErr)
    Call VecAssemblyBegin(RHS, iErr)
    Call VecAssemblyEnd(RHS, iErr)
    Call PetscGetTime(AssembTF, iErr)
    AssembT = AssembT + AssembTF - AssembTS

    Call PetscGetTime(TotTF, iErr)
    TotT = TotTF - TotTS
    Call MPI_REDUCE(SetN, SetNTot, 1, MPI_INTEGER, MPI_SUM, 0,                &
         & PETSC_COMM_WORLD, iErr)

    If (MyRank ==0) Then
       Write(CharBuffer,*) 'Total time in Init/Destroy_Gauss:        ',       &
            & GaussT, '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
       Write(CharBuffer,*) 'Total time in VecSetValue:               ', SetT, &
            & '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
       Write(CharBuffer,*) 'Number of calls to VecSetValue:   ', SetNTot, '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
       Write(CharBuffer,*) 'Total time in VecAssembly:               ',       &
            & AssembT, '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
       Write(CharBuffer,*) 'Total time in Assemb_RHS_Poisson:        ', TotT, &
            & '\n\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
    End If

  End Subroutine Assemb_RHS_Poisson


#ifdef PB_2D
End Module m_OVSch2D_Procs
#else
End Module m_OVSch3D_Procs
#endif
