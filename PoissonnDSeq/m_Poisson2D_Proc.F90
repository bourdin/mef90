
#if defined PB_2D
Module m_Poisson2D_Proc
#elif defined PB_3D
Module m_Poisson3D_Proc
#endif

  Use m_MEF90
  Use m_Poisson_Struct
#if defined PB_2D
  Use m_Poisson2D_Vars
#elif defined PB_3D
  Use m_Poisson3D_Vars
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

  Public :: Init
  Public :: SaveToEXO
  Public :: SaveToEnsight
  Public :: Finalize
  Public :: AssembMR_VLV
  Public :: AssembRHS_VLV

Contains
  Subroutine Init()
    Real(Kind = Kr), Dimension(:), Pointer       :: Tmp_Ptr
  
    !!! MEF90_Initialize initializes MPI, PETSc, and sets some constants.
    Call MEF90_Initialize()

    Write(CharBuffer, 100) 'Simulation name: \n'c
    Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
    If (MEF90_MyRank ==0) Then
       Read(*,100) Params%Sim_Str
    End If

!!! The names of the mesh (.gen) and parameter file (.PARAM)
    Geom%Filename    = Trim(Params%Sim_Str) // '.gen'
    Params%PARAM_Str = Trim(Params%Sim_Str) // '.PARAM'

!!! Read the mesh and the PARAM file
    Call Read_EXO_Geom_Info(Geom)
    Call Read_EXO_Node_Coord(Geom, Node_db, 1)
    Call Read_EXO_Connect(Geom, Elem_db) 

    Call Read_Poisson_EXO_Params(Geom, Params)
    Call Read_Poisson_DATA(Geom, Params)

!!! Print useless data to make the user happy
    Write(CharBuffer,*) 'Number of nodes:                          ',      &
         & Geom%Num_Nodes, '\n'c
    Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
    Write(CharBuffer,*) 'Number of elements:                       ',      &
         & Geom%Num_Elems, '\n'c
    Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)

!!! Initialize the data structure related to the Boundary Conditions
    Call Init_BC()

    Call VecCreateSeq(PETSC_COMM_WORLD, Geom%Num_Nodes, U, iErr) 
    Call VecDuplicate(U, Load, iErr)
    Call VecDuplicate(U, RHS_U, iErr)

!!! Initializes the stiffness matrix. See PETSc' doc for the meaning of the 24.
!!! It is a bad (as in wrong, stupid and lazy) choice
#ifdef PB_2D
    Call MatCreateSeqAIJ(PETSC_COMM_WORLD, Geom%Num_Nodes, Geom%Num_Nodes,    &
         & 24, PETSC_NULL_INTEGER, MR, iErr)
#else
    Call MatCreateSeqAIJ(PETSC_COMM_WORLD, Geom%Num_Nodes, Geom%Num_Nodes,    &
         & 60, PETSC_NULL_INTEGER, MR, iErr)
#endif
    Call MatSetOption(MR, MAT_SYMMETRIC, iErr)

!!! Takes command-line arguments related to the matrix into account
    Call MatSetFromOptions(MR, iErr)
    
!!! Read the load from the .gen file
!!! The load is the first Nodal Variable, and we make only one time step.
    If (MEF90_MyRank == 0) Then
       Call VecGetArrayF90(Load, Load_Ptr, iErr)
       !!! There is something subtle here: (actually, something stupid in 
       !!! Read_EXO_Result_Nodes) Tmp_Ptr is DeAllocated then Allocated 
       !!! in Read_EXO_Result_Nodes. PETSc forbids this.
       Call Read_EXO_Result_Nodes(Geom, 1, 1,  Tmp_Ptr)
       Load_Ptr(:) = Tmp_Ptr(:)
       DeAllocate(Tmp_Ptr)
       Call VecRestoreArrayF90(Load, Load_Ptr, iErr)
       !!! Another way to do so would be to call VecSetValues with argument
       !!! Load_Ptr or Tmp_Ptr.
       !!! That would be much smarter
    End If

    !!! Use the load as the initial value for U
    Call VecCopy(Load, U, iErr)

    !!! Create the KSP
    Call KSPCreate(PETSC_COMM_WORLD, KSP_U, iErr)

    !!! Both the matrix and the preconditionning matrix are MR
    Call KSPSetOperators(KSP_U, MR, MR, SAME_NONZERO_PATTERN, iErr)

    !!! Extract the preconditionner from KSP_U
    Call KSPGetPC(KSP_U, PC_U, iErr)

    !!! Use a conjugated gradient with ILU preconditioner
    !!! Pretty much the fastest solver for this type of problem.
    !!! For a parallel problem, ILU would not work. Replace with PCBJACOBI
    Call PCSetType(PC_U, PCILU, iErr)
    Call KSPSetType(KSP_U, KSPCG, iErr)

    !!! Tell PETSc to use the current value of U_DIst as the startting point of the
    !!! solver
    Call KSPSetInitialGuessNonzero(KSP_U, PETSC_TRUE, iErr)

    !!! Set the relative tolerance from that set in the PARAM file
    Call KSPSetTolerances(KSP_U, Params%TolKSP,                              &
         & PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION,  &
         & PETSC_DEFAULT_INTEGER, iErr)
    
    !!! Takes into account other options potentially given on the command line
    Call KSPSetFromOptions(KSP_U, iErr)
 
100 Format(A)
  End Subroutine Init

  Subroutine Init_BC()
    !!! Originally, I the type of boundary conditions at each node set 
    !!! However, it is more convenient to know for each degree of freedom
    !!! the type of BC (currenly BC_Type_NONE or BC_Type_DIRI)
    !!! this translation is done here

    Integer                        :: iBlk, iN

    Node_db(:)%BC = BC_Type_NONE
    Do iBlk = 1, Geom%Num_Node_Sets
       Do iN = 1, Geom%Node_Set(iBlk)%Num_Nodes
          Node_db(Geom%Node_Set(iBlk)%Node_ID(iN))%BC = Params%BC_Type(iBlk)
       End Do
    End Do
  End Subroutine Init_BC
        
  Subroutine Finalize()
    !!! Exits PETSc cleanly. (Mostly cosmetic)
    DeAllocate(Node_db, Elem_db)
    
    Call MatDestroy(MR, iErr)
    Call VecDestroy(U, iErr)
    Call VecDestroy(Load, iErr)

    Call VecDestroy(RHS_U, iErr)

    Call KSPDestroy(KSP_U, iErr)

    Call MEF90_Finalize()
  End Subroutine Finalize


  Subroutine SaveToEXO()
    !!! Saves and the energy in the .gen file 
    !!! There is space for more variables in this file.


    Call VecGetArrayF90(U, U_Ptr, iErr)
       !!! Saves U in the second nodal slot of the 1st time step of the file
       !!! See m_Poisson_Struct-0.0.1/m_Poisson_Struct.F90 for the structure
       !!! of a Poisson .gen file
    Call Write_EXO_Result_Ptr_Nodes(Geom, 2, 1, U_Ptr)
    Call VecRestoreArrayF90(U, U_Ptr, iErr)
       !!! Writes the energy in the first (and only) global solt of the 1st time 
       !!! step of the .gen file
    Call Write_EXO_Result_Global(Geom, 1, 1, Ener)
  End Subroutine SaveToEXO

  Subroutine SaveToEnsight()
    !!! Saves and the energy in the .gen file 
    !!! There is space for more variables in this file.

#ifdef PB_2D
    Call Write_Geo_2D_Scal_ASCII(Elem_db, Node_db, 'Ensight.geo')
#else
    Call Write_Geo_3D_Scal_ASCII(Elem_db, Node_db, 'Ensight.geo')
#endif

    !!! Usual stuff: U_Dist -> U_Master -> file
    Call VecGetArrayF90(U, U_Ptr, iErr)
       !!! Saves U in the second nodal slot of the 1st time step of the file
       !!! See m_Poisson_Struct-0.0.1/m_Poisson_Struct.F90 for the structure
       !!! of a Poisson .gen file
    Call Write_Result_Scal_ASCII(U_Ptr, 'U.res', 1)
    Call VecRestoreArrayF90(U, U_Ptr, iErr)
       !!! Writes the energy in the first (and only) global solt of the 1st time 
       !!! step of the .gen file

  End Subroutine SaveToEnsight

  Subroutine AssembMR_VLV(MR, Geom, Params, Elems, Nodes)
    Mat                                          :: MR
    Type (EXO_Geom_Info)                         :: Geom
    Type (Poisson_Params)                        :: Params
    
#ifdef PB_2D
    Type (Node2D), Dimension(:), Pointer         :: Nodes
    Type (Element2D_Scal), Dimension(:), Pointer :: Elems
#else
    Type (Node3D), Dimension(:), Pointer         :: Nodes
    Type (Element3D_Scal), Dimension(:), Pointer :: Elems
#endif

    Integer                                      :: Nb_Gauss, Nb_DoF
    Integer                                      :: iSLU1
    Integer                                      :: iSLU2
    Integer                                      :: iE, iG, i, iBlk

    PetscScalar, Dimension(:,:), Pointer         :: MR_Elem

    PetscTruth                                   :: ISAssembled

    Nb_DoF = Geom%Elem_blk(1)%Num_Nodes_per_elem
    
    Call MatAssembled(MR, IsAssembled, iErr)
    If (IsAssembled) Then
       Call MatZeroEntries(MR, iErr)
    End If
    !!! Non-BC Part of the matrix
    Allocate (MR_Elem(Nb_DoF, Nb_DoF))

    Do_iE: Do iE = 1, Geom%Num_Elems
       MR_Elem = 0.0_Kr

       !!! On large problems, it is not possible to store the base function 
       !!! informations for each element. Instead, I recompute them at assembly
       !!! time
       Call Init_Gauss_EXO(Elems, Nodes, Geom, GaussOrder, Elem=iE)

       Do_iGU: Do iG = 1, Elems(iE)%Nb_Gauss
          Do_iSLU1: Do iSLU1 = 1, Elems(iE)%Nb_DoF
             DoiSLU2: Do iSLU2 = 1, Elems(iE)%Nb_DoF
                   MR_Elem(iSLU1, iSLU2) = MR_Elem(iSLU1, iSLU2) + Elems(iE)%Gauss_C(iG) *                                         &
                                         ( Elems(iE)%Grad_BF(iSLU1, iG) .DotP. Elems(iE)%Grad_BF(iSLU2, iG) ) 
             End Do DoiSLU2
          End Do Do_iSLU1
       End Do Do_iGU

       Call MatSetValues(MR, Nb_DoF, Elems(iE)%Id_DoF-1, Nb_DoF, Elems(iE)%Id_DoF-1, MR_Elem, ADD_VALUES, iErr)
       Call Destroy_Gauss_EXO(Elems, Elem=iE)
    EndDo Do_iE


    Call MatAssemblyBegin(MR, MAT_FLUSH_ASSEMBLY, iErr)
    DeAllocate(MR_Elem)
    Call MatAssemblyEnd(MR, MAT_FLUSH_ASSEMBLY, iErr)
  
    !!! Dirichlet part
    Do i = 1, Geom%Num_Node_Sets
       iBlk = Geom%Node_Set(i)%ID
       If (Params%BC_Type(iBlk) == BC_Type_DIRI) Then
          Do iSLU1 = 1, Geom%Node_Set(i)%Num_Nodes
             iSLU2 = Geom%Node_Set(i)%Node_ID(iSLU1)
             Call MatSetValue(MR, iSLU2, iSLU2, VLV, INSERT_VALUES, iErr)
          End Do
       End If
    End Do

    Call MatAssemblyBegin(MR, MAT_FINAL_ASSEMBLY, iErr)
    Call MatAssemblyEnd(MR, MAT_FINAL_ASSEMBLY, iErr)
  End Subroutine AssembMR_VLV


  Subroutine AssembRHS_VLV(RHS, Geom, Params, Elems, Nodes)
    Vec                                          :: RHS
    Type (EXO_Geom_Info)                         :: Geom
    Type (Poisson_Params)                        :: Params

#ifdef PB_2D
    Type (Node2D), Dimension(:), Pointer         :: Nodes
    Type (Element2D_Scal), Dimension(:), Pointer :: Elems
#else
    Type (Node3D), Dimension(:), Pointer         :: Nodes
    Type (Element3D_Scal), Dimension(:), Pointer :: Elems
#endif

    Integer                                      :: Nb_Gauss, Nb_DoF
    Integer                                      :: iSL, iSG
    Integer                                      :: iE, iG, i
    Integer                                      :: iSet, iS
    Real(Kind = Kr), Dimension(:), Pointer       :: RHS_Elem
    Real(Kind = Kr)                              :: F

    PetscScalar                                  :: Tmp_Val

    Call VecSet(RHS, 0.0_Kr, iErr)

    Call VecGetArrayF90(Load, Load_Ptr, iErr)

    Nb_DoF = Geom%Elem_Blk(1)%Num_Nodes_Per_Elem
    Allocate (RHS_Elem(Nb_DoF))
    Do iE = 1, Geom%Num_Elems
       RHS_Elem = 0.0_Kr
       Call Init_Gauss_EXO(Elems, Nodes, Geom, GaussOrder, Elem=iE)
       Nb_Gauss = Elems(iE)%Nb_Gauss
       Do_iGU: Do iG = 1, Nb_Gauss
          F = 0.0_Kr
          Do_iSLU1: Do iSL = 1, Elems(iE)%Nb_DoF
             iSG = Elems(iE)%ID_DoF(iSL)
             F = F + Elems(iE)%BF(iSL, iG) * Load_Ptr(iSG)
          End Do Do_iSLU1
          Do_iSLU2: Do iSL = 1, Elems(iE)%Nb_DoF
             iSG = Elems(iE)%ID_DoF(iSL)
             RHS_Elem(iSL) = RHS_Elem(iSL) + Elems(iE)%Gauss_C(iG) * Elems(iE)%BF(iSL, iG) * F
          End Do Do_iSLU2
       End Do Do_iGU
       Call VecSetValues(RHS, Nb_DoF, Elems(iE)%ID_DoF-1, RHS_Elem, ADD_VALUES, iErr)
       Call Destroy_Gauss_EXO(Elems, Elem=iE)
    End Do
    DeAllocate (RHS_Elem)
                
    Call VecAssemblyBegin(RHS, iErr)
    Call VecRestoreArrayF90(Load, Load_Ptr, iErr)
    Call VecAssemblyEnd(RHS, iErr)
    
    Do iSet = 1, Geom%Num_Node_Sets
       If (Params%BC_Type(iSet) == BC_Type_NONE) Then
          Cycle
       End If
       Do iS = 1, Geom%Node_Set(iSet)%Num_Nodes
          iSG = Geom%Node_Set(iSet)%Node_ID(iS)
          Call VecSetValue(RHS, iSG, VLV * Params%BC_DB(iSet)%Dis_F(iS), INSERT_VALUES, iErr)
       End Do
    End Do
    Call VecAssemblyBegin(RHS, iErr)
    Call VecAssemblyEnd(RHS, iErr)

  End Subroutine AssembRHS_VLV

#if defined PB_2D
End Module m_Poisson2D_Proc
#elif defined PB_3D
End Module m_Poisson3D_Proc
#endif
