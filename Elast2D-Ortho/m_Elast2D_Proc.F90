Module m_Elast2D_Proc
  Use m_MEF90
  Use m_Rupt_Struct

  Use m_Elast2D_Vars

  Implicit NONE
  PRIVATE

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscao.h"

  Public :: Init
  Public :: Init_KSP
  Public :: Finalize
  Public :: Export

  Public :: Assemb_Mat_Elast
  Public :: Assemb_RHS_Elast
  Public :: Calc_Ener

  Public :: GenVect
Contains


  Subroutine Init()
    Integer, Dimension(:), Pointer     :: MyIndices
    Integer                            :: i, iS, iE

    PetscLogDouble                     :: BCastTS, BCastTF
    PetscLogDouble                     :: InitVectTS, InitVectTF
    PetscTruth                         :: Has_Sim_Str

    Call MEF90_Initialize()
    MEF90_GaussOrder = 2 


    Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-f', Params%Sim_Str, Has_Sim_Str, iErr)
    
    If (.NOT. Has_Sim_Str) Then
       Write(CharBuffer, 100) 'Simulation name: \n'c
       Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
       If (MEF90_MyRank ==0) Then
          Read(*,100) Params%Sim_Str
       End If
       Call MPI_BCAST(Params%Sim_Str, MXSTLN, MPI_CHARACTER, 0, MPI_COMM_WORLD, iErr)
    End If


    Geom%Filename    = Trim(Params%Sim_Str) // '.gen'
    Params%PARAM_Str = Trim(Params%Sim_Str) // '.PARAM'
    Params%CST_Str   = Trim(Params%Sim_Str) // '.CST'
    
    
    Call Read_EXO_Geom_Info(Geom)
    Call Read_EXO_Node_Coord(Geom, Node_db, Geom%Num_Dim)
    Call Read_EXO_Node_Coord(Geom, Node_Scal, 1)    
    Call Read_EXO_Connect(Geom, Elem_db) 
    Call Read_EXO_Connect(Geom, Elem_Scal)
    Call Read_Rupt_EXO_Params(Geom, Params)
    Call Read_Rupt_DATA(Geom, Params)

    Write(CharBuffer,*) 'Number of nodes:                          ', Geom%Num_Nodes, '\n'c
    Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
    Write(CharBuffer,*) 'Number of elements:                       ', Geom%Num_Elems, '\n'c
    Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
    Call Init_SD_NoOvlp(Geom, Elem_db, MySD_Vect, SOL_Dist, SOL_Loc, SOL_Master)
    Call Init_SD_NoOvlp(Geom, Elem_Scal, MySD_Scal, Temp_Dist, Temp_Loc, Temp_Master)

    Call Init_BC(Geom, Params, Node_db)

    Call VecDuplicate(Sol_Dist, RHS, iErr)
    Call VecGhostGetLocalForm(RHS, RHS_Loc, iErr)

    Call VecDuplicate(Sol_Master, BC_Master, iErr)
    Call VecDuplicate(Sol_Dist, BC_Dist, iErr)
    Call VecGhostGetLocalForm(BC_Dist, BC_Loc, iErr)

    Call VecDuplicate(Sol_Master, F_Master, iErr)
    Call VecDuplicate(Sol_Dist, F_Dist, iErr)
    Call VecGhostGetLocalForm(F_Dist, F_Loc, iErr)

    Call MatCreateMPIAIJ(PETSC_COMM_WORLD, MySD_Vect%Num_Nodes, MySD_Vect%Num_Nodes, Geom%Num_Nodes * Geom%Num_Dim,                &
    	 & Geom%Num_Nodes * Geom%Num_Dim, 60, PETSC_NULL_INTEGER, 60, PETSC_NULL_INTEGER, MR, iErr)

    Call MatSetOption(MR, MAT_SYMMETRIC, iErr)
    Call MatSetFromOptions(MR, iErr)

    Allocate (Stress_Sol(Geom%Num_Elems))
    Allocate (Strain_Sol(Geom%Num_Elems))


100 Format(A)
200 Format('Rank: ', I4,' Index Range: ', I7,' ',I7, '\n'c) 
  End Subroutine Init

  Subroutine Init_KSP()
    KSP                :: Sub_KSP_MR
    PC                 :: Sub_PC_MR
    PCType             :: PC_Type
    Call KSPCreate(PETSC_COMM_WORLD, KSP_MR, iErr)
    Call KSPSetOperators(KSP_MR, MR, MR, SAME_NONZERO_PATTERN, iErr)
    Call KSPGetPC(KSP_MR, PC_MR, iErr)
    Call PCSetType(PC_MR, PCBJACOBI, iErr)
    Call KSPSetType(KSP_MR, KSPCG, iErr)

    Call KSPSetInitialGuessNonzero(KSP_MR, PETSC_TRUE, iErr)

    Call KSPSetTolerances(KSP_MR, Params%TolKSP, PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION,                   &
                          PETSC_DEFAULT_INTEGER, iErr)
    Call KSPSetFromOptions(KSP_MR, iErr)

    Call PCGetType(PC_MR, PC_Type, iErr)
    If (PC_Type == PCBJACOBI) Then
       Call KSPSetUp(KSP_MR, iErr)
       Call PCBJacobiGetSubKSP(PC_MR, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, Sub_KSP_MR, iErr)
       Call KSPGetPC(Sub_KSP_MR, Sub_PC_MR, iErr)
       Call PCFactorSetZeroPivot(Sub_PC_MR, 1.0D-20, iErr)
       Call PCSetFromOptions(Sub_PC_MR, iErr)
    End If
  End Subroutine Init_KSP

  Subroutine Export()
    Real(Kind = Kr), Dimension(:), Pointer        :: SOL_Ptr
    
    Call VecScatterBegin(SOL_Dist, SOL_Master, INSERT_VALUES, SCATTER_FORWARD, MySD_Vect%ToMaster, iErr)
    Call VecScatterEnd(SOL_Dist, SOL_Master, INSERT_VALUES, SCATTER_FORWARD, MySD_Vect%ToMaster, iErr)
    
    If (MEF90_MyRank == 0) Then
       Call VecGetArrayF90(SOL_Master, Sol_Ptr, iErr)
       Call Write_EXO_Result_Ptr_Nodes(Geom, 2, TimeStep, SOL_Ptr)
       Write(*,90) MinVal(Sol_Ptr), MaxVal(SOL_Ptr)
       Call VecRestoreArrayF90(SOL_Master, SOL_Ptr, iErr)
    End If


90  Format('Displacement Min / max: ', 2(ES10.3,' '))
  End Subroutine Export

  Subroutine Finalize()
    DeAllocate(Node_db, Elem_db)
    DeAllocate(Stress_Sol, Strain_Sol)

    Call MatDestroy(MR, iErr)
    Call VecDestroy(RHS, iErr)
    Call VecDestroy(RHS_Loc, iErr)

    Call VecDestroy(SOL_Dist, iErr)
    Call VecDestroy(SOL_Loc, iErr)
    Call VecDestroy(SOL_Master, iErr)

    Call VecDestroy(BC_Dist, iErr)
    Call VecDestroy(BC_Loc, iErr)
    Call VecDestroy(BC_Master, iErr)
    
    Call VecDestroy(F_Dist, iErr)
    Call VecDestroy(F_Loc, iErr)
    Call VecDestroy(F_Master, iErr)

    Call PETScFinalize(iErr)
  End Subroutine Finalize

  Subroutine Init_BC(Geom, Params, Node_db)
    Type(EXO_Geom_Info), Intent(IN)                  :: Geom
    Type(Rupt_Params_Hooke2D), Intent(IN)            :: Params

    Type(Node2D), Dimension(:), Pointer              :: Node_db
    Integer                                          :: iN, iSet

    Select Case(Geom%Numbering)
    Case(Numbering_PerCoord)
       Node_db(:)%BC = BC_Type_NONE
       Do iSet = 1, Geom%Num_Node_Sets
          Do iN = 1, Geom%Node_Set(iSet)%Num_Nodes
             Node_db(Geom%Node_Set(iSet)%Node_ID(iN))%BC =                    &
                  & Params%BC_Type_X(iSet)
             Node_db(Geom%Node_Set(iSet)%Node_ID(iN)+Geom%Num_Nodes)%BC =     &
                  & Params%BC_Type_Y(iSet)
          End Do
       End Do
    Case(Numbering_PerNodes)
       Node_db(:)%BC = BC_Type_NONE
       Do iSet = 1, Geom%Num_Node_Sets
          Do iN = 1, Geom%Node_Set(iSet)%Num_Nodes
             Node_db(Geom%Num_Dim * (Geom%Node_Set(iSet)%Node_ID(iN)-1)+1)%BC &
                  & = Params%BC_Type_X(iSet)
             Node_db(Geom%Num_Dim * (Geom%Node_Set(iSet)%Node_ID(iN)-1)+2)%BC &
                  & = Params%BC_Type_Y(iSet)
          End Do
       End Do
    Case Default
       Write(*,*) 'Init_BC: Unknown numbering scheme: ', Geom%Numbering
       STOP
    End Select
  End Subroutine Init_BC


  Subroutine Assemb_Mat_Elast(MR, Geom, Params, Elem_db, Node_db)
    Mat                                                 :: MR
    Type (EXO_Geom_Info)                                :: Geom
    Type (Rupt_Params_Hooke2D)                          :: Params
    Type (Node2D), Dimension(:), Pointer                :: Node_db 
    Type (Element2D_Elast), Dimension(:), Pointer       :: Elem_db 
    Type (MatS2D)                                       :: Sigma, Epsilon

    Integer                                             :: Nb_Gauss, Nb_DoF
    Integer                                             :: iSL1, iSL2
    Integer                                             :: iSG1, iSG2
    Integer                                             :: iE, iG, iELoc
    Integer                                             :: iBlk
    Integer                                             :: i

    PetscScalar, Dimension(:,:), Pointer                :: MR_Elem
    Integer, Dimension(:), Pointer                      :: EXO_Indices
    PetscTruth                                          :: IsAssembled

    Call MatAssembled(MR, IsAssembled, iErr)
    If (IsAssembled) Then
       Call MatZeroEntries(MR, iErr)
    End If
    
    ! Assembly of the Non BC terms 
    Do_iBlk: Do iBlk = 1, Geom%Num_Elem_Blks

       Nb_DoF = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem * Geom%Num_Dim
       Allocate (MR_Elem(Nb_DoF, Nb_DoF))
       Allocate (EXO_Indices(Nb_DoF))

       Do_iE: Do iELoc = 1, Geom%Elem_Blk(iBlk)%Num_Elems
          iE = Geom%Elem_Blk(iBlk)%ELem_ID(iELoc)
          If (.NOT. MySD_Vect%IsLocal_Elem(iE)) Then
             CYCLE
          End If
          EXO_Indices = Elem_db(iE)%Id_DoF-1
          Call AOApplicationToPETSc(MySD_Vect%EXO_AO, Nb_DoF, EXO_Indices, iErr)

          Call Init_Gauss_EXO(Elem_db, Node_db, Geom, GaussOrder, Elem=iE)

          MR_Elem = 0.0_Kr
          Nb_Gauss = Elem_db(iE)%Nb_Gauss
          Do_iG: Do iG = 1, Nb_Gauss
             Do_iSL1: Do iSL1 = 1, Nb_DoF
                ISG1 = Elem_db(iE)%ID_DoF(iSL1)
                Epsilon = Elem_db(iE)%GradS_BF(iSL1,iG)
                Sigma = Params%Hookes_Law(iBlk) * Epsilon

                Do_iSL2: Do iSL2 = 1, Nb_DoF
                   iSG2 = Elem_db(iE)%ID_DoF(iSL2)
                   Epsilon = Elem_db(iE)%GradS_BF(iSL2,iG)
                   MR_Elem(iSL2, iSL1) = MR_Elem(iSL2, iSL1) + Elem_db(iE)%Gauss_C(iG) * (Sigma .DotP. Epsilon)
                End Do Do_iSL2
             End Do Do_iSL1
          End Do Do_iG

          Call MatSetValues(MR, Nb_DoF, EXO_Indices, Nb_DoF, EXO_Indices, MR_Elem, ADD_VALUES, iErr)
          Call Destroy_Gauss_EXO(Elem_db, Elem=iE)
       EndDo Do_iE

       DeAllocate (MR_Elem)
       DeAllocate (EXO_Indices)
    End Do Do_iBlk

    Call MatAssemblyBegin(MR, MAT_FLUSH_ASSEMBLY, iErr)  
    Call MatAssemblyEnd(MR, MAT_FLUSH_ASSEMBLY, iErr)  
    
    ! Assembly of the BC terms 
    Allocate(EXO_Indices(Geom%Num_Nodes * Geom%Num_Dim))
    EXO_Indices = (/ (i ,i = 0, Geom%Num_Nodes * Geom%Num_Dim - 1) /)
    Call AOApplicationToPETSc(MySD_Vect%EXO_AO, Geom%Num_Nodes * Geom%Num_Dim, EXO_Indices, iErr)
    Do_iBlk_BC: Do iBlk = 1, Geom%Num_Elem_Blks
       Nb_DoF = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem * Geom%Num_Dim

       Do_iE_BC: Do iELoc = 1, Geom%Elem_blk(iBlk)%Num_Elems
          iE = Geom%Elem_blk(iBlk)%Elem_ID(iELoc)
          If (.NOT. MySD_Vect%IsLocal_Elem(iE)) Then
             CYCLE
          End If
          
          Do_iSL1_BC: Do iSL1 = 1, Nb_DoF
             ISG1 = Elem_db(iE)%ID_DoF(iSL1)
             Is_BC_BC: If ( (Node_db(iSG1)%BC /= BC_Type_NONE) .AND. (MySD_Vect%IsLocal_Node(iSG1)) )Then
                Call MatSetValue(MR, EXO_Indices(iSG1), EXO_Indices(iSG1), VLV, INSERT_VALUES, iErr)
             End If Is_BC_BC
          End Do Do_iSL1_BC
       End Do Do_iE_BC
    End Do Do_iBlk_BC
    DeAllocate (EXO_Indices)

    Call MatAssemblyBegin(MR, MAT_FINAL_ASSEMBLY, iErr)
    Call MatAssemblyEnd(MR, MAT_FINAL_ASSEMBLY, iErr)  
  End Subroutine Assemb_Mat_Elast

  Subroutine Assemb_RHS_Elast(RHS, Geom, Params, Elem_Vect, Node_Vect, Elem_Scal, Node_Scal, BC_Loc, F_Loc, Temp_Loc)
    Vec                                                 :: RHS
    Type (EXO_Geom_Info)                                :: Geom
    Type (Rupt_Params_Hooke2D)                          :: Params
    Type(Element2D_Elast), Dimension(:), Pointer        :: Elem_Vect
    Type(Node2D), Dimension(:), Pointer                 :: Node_Vect
    Type(Element2D_Scal), Dimension(:), Pointer         :: Elem_Scal
    Type(Node2D), Dimension(:), Pointer                 :: Node_Scal
    Type(MatS2D)                                        :: ThetaId, AThetaId
    Type(Vect2D)                                        :: F
    Vec                                                 :: BC_Loc
    Vec                                                 :: F_Loc
    Vec                                                 :: Temp_Loc

    Real(Kind = Kr), Dimension(:), Pointer              :: BC_Ptr
    Real(Kind = Kr), Dimension(:), Pointer              :: RHS_Ptr
    Real(Kind = Kr), Dimension(:), Pointer              :: F_Ptr
    Real(Kind = Kr), Dimension(:), Pointer              :: Temp_Ptr
    

    Integer, Dimension(:), Pointer                      :: Loc_Indices_Vect
    Integer, Dimension(:), Pointer                      :: EXO_Indices_Vect
    Integer, Dimension(:), Pointer                      :: Loc_Indices_Scal

    Integer                                             :: Nb_Gauss, Nb_DoF_Vect
    Integer                                             :: iSL, iSG
    Integer                                             :: iE, iELoc, iG
    Integer                                             :: iBlk
    
    Integer                                             :: i, iS

    Real(Kind = Kr), Dimension(:), Pointer              :: RHS_Elem

    Allocate(Loc_Indices_Vect(Geom%Num_Nodes * Geom%Num_Dim))
    Loc_Indices_Vect = (/ (i ,i = 0, Geom%Num_Nodes * Geom%Num_Dim - 1) /)
    Call AOApplicationToPETSc(MySD_Vect%Loc_AO, Geom%Num_Nodes * Geom%Num_Dim, Loc_Indices_Vect, iErr)

    Allocate(Loc_Indices_Scal(Geom%Num_Nodes))
    Loc_Indices_Scal = (/ (i ,i = 0, Geom%Num_Nodes - 1) /)
    Call AOApplicationToPETSc(MySD_Scal%Loc_AO, Geom%Num_Nodes, Loc_Indices_Scal, iErr)

    Call VecGetArrayF90(F_Loc, F_Ptr, iErr)
    Call VecGetArrayF90(Temp_Loc, Temp_Ptr, iErr)


    Call VecSet(RHS, 0.0_Kr, iErr)
    Do_iBlk: Do iBlk = 1, Geom%Num_elem_blks
       Nb_DoF_Vect = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem * Geom%Num_Dim
       Allocate (RHS_Elem(Nb_DoF_Vect))
       Allocate (EXO_Indices_Vect(Nb_DoF_Vect))

       Do_iE: Do iELoc = 1, Geom%Elem_blk(iBlk)%Num_Elems
          iE = Geom%Elem_blk(iBlk)%Elem_ID(iELoc)
          If (.NOT. MySD_Vect%IsLocal_Elem(iE)) Then
             CYCLE
          End If
          EXO_Indices_Vect = Elem_Vect(iE)%Id_DoF-1
          Call AOApplicationToPETSc(MySD_Vect%EXO_AO, Nb_DoF_Vect, EXO_Indices_Vect, iErr)

          Call Init_Gauss_EXO(Elem_Vect, Node_Vect, Geom, GaussOrder, Elem=iE)
          Call Init_Gauss_EXO(Elem_Scal, Node_Scal, Geom, GaussOrder, Elem=iE)

          RHS_Elem = 0.0_Kr
          Do_iG: Do iG = 1, Elem_Vect(iE)%Nb_Gauss

!!! Need to think about thermal expansion with generalized Hookes laws          
!             ThetaId = 0.0_Kr
!             Do_iSL1: Do iSL = 1, Elem_Scal(iE)%Nb_DoF
!                iSG = Elem_Scal(iE)%ID_DoF(iSL)
!                ThetaId%XX = ThetaId%XX + Params%Therm_Exp(iBlk) * Temp_Ptr(Loc_Indices_Scal(iSG)+1) * Elem_Scal(iE)%BF(iSL, iG)
!                ThetaId%YY = ThetaId%YY + Params%Therm_Exp(iBlk) * Temp_Ptr(Loc_Indices_Scal(iSG)+1) * Elem_Scal(iE)%BF(iSL, iG)
!             End Do Do_iSL1
!             AThetaId = 2.0_Kr * K2 * ThetaId
!             AThetaId%XX = AThetaId%XX + K1 * Trace(ThetaId)
!             AThetaId%YY = AThetaId%YY + K1 * Trace(ThetaId)
             F = 0.0_Kr
             Do_iSL2: Do iSL = 1, Elem_Vect(iE)%Nb_DoF
                iSG = Elem_Vect(iE)%ID_DoF(iSL)
                F = F + F_Ptr(Loc_Indices_Vect(iSG)+1) * Elem_Vect(iE)%BF(iSL, iG)
             End Do Do_iSL2
             
             Do_iSL3: Do iSL = 1, Elem_Vect(iE)%Nb_DoF
                RHS_Elem(iSL) = RHS_Elem(iSL) + Elem_Vect(iE)%Gauss_C(iG) *                                                        &
                                ( (F .DotP. Elem_Vect(iE)%BF(iSL, iG)) + (AThetaId .DotP. Elem_Vect(iE)%GradS_BF(iSL, iG)) )
             End Do Do_iSL3
          End Do Do_iG
          Call VecSetValues(RHS, Elem_Vect(iE)%Nb_DoF, EXO_Indices_Vect, RHS_Elem, ADD_VALUES, iErr)
         
          Call Destroy_Gauss_EXO(Elem_Vect, Elem=iE)
          Call Destroy_Gauss_EXO(Elem_Scal, Elem=iE)
       End Do Do_iE
       DeAllocate(RHS_Elem)
       DeAllocate(EXO_Indices_Vect)
    End Do Do_iBlk

    Call VecAssemblyBegin(RHS, iErr)
    DeAllocate(Loc_Indices_Scal)
    Call VecAssemblyEnd(RHS, iErr)
    Call VecRestoreArrayF90(F_Loc, F_Ptr, iErr)
    Call VecRestoreArrayF90(Temp_Loc, Temp_Ptr, iErr)
    
    !!! BC Part, using VLV
    !!! There are two ways to do that:
    !!! Use VecSetValue and EXO_Indices or UseVecGetArray and Loc_Indices
    !!! I chose the second way for no special reason.

    Call VecGetArrayF90(BC_Loc, BC_Ptr, iErr)
    Call VecGetArrayF90(RHS, RHS_Ptr, iErr)

    Do_iS: Do iS = 1, Geom%Num_Nodes * Geom%Num_Dim
       Is_BC: If ((Node_Vect(iS)%BC /= BC_Type_NONE) .AND. (MySD_Vect%IsLocal_Node(iS)) ) Then
          RHS_Ptr(Loc_Indices_Vect(iS)+1) = BC_Ptr(Loc_Indices_Vect(iS)+1)*VLV
       End If Is_BC
    End Do Do_iS

    DeAllocate(Loc_Indices_Vect)
    Call VecRestoreArrayF90(BC_Loc, BC_Ptr, iErr)
    Call VecRestoreArrayF90(RHS, RHS_Ptr, iErr)
  End Subroutine Assemb_RHS_Elast

  Subroutine GenVect(Geom, SD, Nodes, Elems, Vect_Dist)
    Type (EXO_Geom_Info)                                :: Geom
    Type (SD_Info)                                      :: SD
    Type(Element2D_Elast), Dimension(:), Pointer        :: Elems
    Type(Node2D), Dimension(:), Pointer                 :: Nodes
    Vec                                                 :: Vect_Dist

    Integer                                             :: iS, i
    Integer, Dimension(:), Pointer                      :: EXO_Indices
    Real(Kind = Kr)                                     :: X, Y, Z
    

    EXO_Indices = (/ (i ,i = 0, Geom%Num_Nodes *Geom%Num_Dim - 1) /)
    Call AOApplicationToPETSc(MySD_Vect%EXO_AO, Geom%Num_Nodes* Geom%Num_Dim, EXO_Indices, iErr)

    Do iS = 0, Geom%Num_Nodes-1
       If ( .NOT. SD%IsLocal_Node(iS * Geom%Num_Dim +1) ) Then
          CYCLE
       EndIf
       X = Nodes(Geom%Num_Dim * iS +1)%Coord%X
       Y = Nodes(Geom%Num_Dim * iS +1)%Coord%Y
       Call VecSetValue(Vect_Dist, EXO_Indices(Geom%Num_Dim * iS +1), X, INSERT_VALUES, iErr)
    End Do

    Call VecAssemblyBegin(Vect_Dist, iErr)
    Call VecAssemblyEnd(Vect_Dist, iErr)

    DeAllocate(EXO_Indices)
  End Subroutine GenVect


  Subroutine Calc_Ener(DISP_Loc, Geom, Params, Elem_db_Vec, Node_db_Vec, SD_Vec, Elem_db_Scal, Node_db_Scal, SD_Scal,              & 
                       F_Loc, Temp_Loc, Ener) 

    Type(MatS2D)                                        :: Sigma
    Type(MatS2D)                                        :: Epsilon
    Type(Element2D_Elast), Dimension(:), Pointer        :: Elem_db_Vec 
    Type(Node2D), Dimension(:), Pointer                 :: Node_db_Vec
    Type(Vect2D)                                        :: F, U
    Type(Element2D_Scal), Dimension(:), Pointer         :: Elem_db_Scal 
    Type(Node2D), Dimension(:), Pointer                 :: Node_db_Scal
    Type (EXO_Geom_Info)                                :: Geom
    Type (Rupt_Params_Hooke2D)                          :: Params
    Type (SD_Info)                                      :: SD_Vec, SD_Scal
    
    Vec                                                 :: DISP_Loc
    Vec                                                 :: F_Loc
    Vec                                                 :: Temp_Loc
	    
    Real(Kind = Kr)                                     :: Ener
    Real(Kind = Kr)                                     :: MyEner

    Integer                                             :: Nb_Gauss, Nb_DoF_Vec, Nb_DoF_Scal
    Integer                                             :: iSL, iSG
    Integer                                             :: iE, iG, iELoc
    Integer(Kind = Ki)                                  :: iBlk
    Integer                                             :: i
    Integer, Dimension(:), Pointer                      :: Loc_Indices_Vec
    Integer, Dimension(:), Pointer                      :: Loc_Indices_Scal
    Real(Kind = Kr), Dimension(:), Pointer              :: DISP_Ptr
    Real(Kind = Kr), Dimension(:), Pointer              :: F_Ptr
    Real(Kind = Kr), Dimension(:), Pointer              :: Temp_Ptr
    
    Ener = 0.0_Kr
    MyEner = 0.0_Kr
    Allocate(Loc_Indices_Vec(Geom%Num_Nodes * Geom%Num_Dim))
    Loc_Indices_Vec = (/ (i ,i = 0, Geom%Num_Nodes * Geom%Num_Dim - 1) /)
    Call AOApplicationToPETSc(SD_Vec%Loc_AO, Geom%Num_Nodes * Geom%Num_Dim, Loc_Indices_Vec, iErr)

    Allocate(Loc_Indices_Scal(Geom%Num_Nodes))
    Loc_Indices_Scal = (/ (i ,i = 0, Geom%Num_Nodes - 1) /)
    Call AOApplicationToPETSc(SD_Scal%Loc_AO, Geom%Num_Nodes, Loc_Indices_Scal, iErr)

    Call VecGetArrayF90(DISP_Loc, DISP_Ptr, iErr)
    Call VecGetArrayF90(F_Loc, F_Ptr, iErr)
    Call VecGetArrayF90(Temp_Loc, Temp_Ptr, iErr)

    Do_iBlk: Do iBlk = 1, Geom%Num_Elem_Blks

       Nb_DoF_Vec  = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem * Geom%Num_Dim
       Nb_DoF_Scal = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem
       Do_iE: Do iELoc = 1, Geom%Elem_Blk(iBlk)%Num_Elems
          iE = Geom%Elem_Blk(iBlk)%Elem_ID(iELoc)
          If (.NOT. SD_Vec%IsLocal_Elem(iE)) Then
             CYCLE
          End If

          Call Init_Gauss_EXO(Elem_db_Vec, Node_db_Vec, Geom, GaussOrder, Elem=iE)
          Call Init_Gauss_EXO(Elem_db_Scal, Node_db_Scal, Geom, GaussOrder, Elem=iE)

          Nb_Gauss = Elem_db_Vec(iE)%Nb_Gauss
          Do iG = 1, Nb_Gauss
             Sigma   = 0.0_Kr
             Epsilon = 0.0_Kr
             F       = 0.0_Kr
             U       = 0.0_Kr
             Do_iSL1: Do iSL = 1, Nb_DoF_Vec
                iSG = Elem_db(iE)%ID_DoF(iSL)
                Epsilon = Epsilon + Elem_db_Vec(iE)%GradS_BF(iSL,iG) * Disp_Ptr(Loc_Indices_Vec(iSG)+1)
             End Do Do_iSL1

!!! Need to think about Thermal stuff
!             Do_iSL2: Do iSL = 1, Nb_DoF_Scal
!                iSG = Elem_db_Scal(iE)%ID_DoF(iSL)
!                Epsilon%XX = Epsilon%XX - Params%Therm_Exp(iBlk) * Elem_db_Scal(iE)%BF(iSL,iG) * Temp_Ptr(Loc_Indices_Scal(iSG)+1)
!                Epsilon%YY = Epsilon%YY - Params%Therm_Exp(iBlk) * Elem_db_Scal(iE)%BF(iSL,iG) * Temp_Ptr(Loc_Indices_Scal(iSG)+1)
!             End Do Do_iSL2

             Sigma = Params%Hookes_Law(iBlk) * Epsilon
             
             MyEner = MyEner + Elem_db(iE)%Gauss_C(iG) * ( Sigma .DotP. Epsilon ) * .5_Kr
            
!!! Forces stuff
             If (Params%Has_Force(iBlk)) Then
                Do_iSL3: Do iSL = 1, Nb_DoF_Vec
                   iSG = Elem_db(iE)%ID_DoF(iSL)
                   U = U + DISP_Ptr(Loc_Indices_Vec(iSG)+1) * Elem_db_Vec(iE)%BF(iSL,iG)
                   F = F + F_Ptr(Loc_Indices_Vec(iSG)+1) * Elem_db_Vec(iE)%BF(iSL,iG)
                End Do Do_iSL3
                MyEner = MyEner - Elem_db(iE)%Gauss_C(iG) * (F .DotP. U) 
              End If
          End Do

          Call Destroy_Gauss_EXO(Elem_db_Vec, Elem=iE)
          Call Destroy_Gauss_EXO(Elem_db_Scal, Elem=iE)
       EndDo Do_iE

    End Do Do_iBlk
    DeAllocate (Loc_Indices_Vec)
    DeAllocate (Loc_Indices_Scal)
    Call VecRestoreArrayF90(DISP_Loc, DISP_Ptr, iErr)
    Call VecRestoreArrayF90(F_Loc, F_Ptr, iErr)
    Call VecRestoreArrayF90(Temp_Loc, Temp_Ptr, iErr)

    Call MPI_AllReduce(MyEner, Ener, 1, MPI_DOUBLE_PRECISION, MPI_SUM, PETSC_COMM_WORLD, iErr)
  End Subroutine Calc_Ener


End Module m_Elast2D_Proc
