Module m_Elast2DA_Proc
  Use m_MEF90
  Use m_Rupt_Struct

  Use m_Elast2DA_Vars

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

  Public :: Assemb_Mat_Elast2DA
  Public :: Assemb_RHS_Elast2DA
!  Public :: Calc_Ener

Contains


  Subroutine Init()
    Integer, Dimension(:), Pointer     :: MyIndices
    Integer                            :: i, iS, iE

    PetscLogDouble                     :: BCastTS, BCastTF
    PetscLogDouble                     :: InitVectTS, InitVectTF
    PetscTruth                :: Has_Sim_Str

    Call MEF90_Initialize()
    MEF90_GaussOrder = 2 


    Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-f', Params%Sim_Str,     &
         & Has_Sim_Str, iErr)
    
    If (.NOT. Has_Sim_Str) Then
       Write(CharBuffer, 100) 'Simulation name: \n'c
       Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
       If (MEF90_MyRank ==0) Then
          Read(*,100) Params%Sim_Str
       End If
       Call MPI_BCAST(Params%Sim_Str, MXSTLN, MPI_CHARACTER, 0, MPI_COMM_WORLD,  &
            & iErr)
    End If


    Geom%Filename    = Trim(Params%Sim_Str) // '.gen'
    Params%PARAM_Str = Trim(Params%Sim_Str) // '.PARAM'
    Params%CST_Str   = Trim(Params%Sim_Str) // '.CST'
    
    
    Call Read_EXO_Geom_Info(Geom)
    Call Read_EXO_Node_Coord(Geom, Node_db, 1)
    Call Read_EXO_Connect(Geom, Elem_db) 
    Call Read_Rupt_EXO_Params(Geom, Params)
    Call Read_Rupt_DATA(Geom, Params)

    Write(CharBuffer,*) 'Number of nodes:                          ',      &
         & Geom%Num_Nodes, '\n'c
    Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
    Write(CharBuffer,*) 'Number of elements:                       ',      &
         & Geom%Num_Elems, '\n'c
    Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
    Call Init_SD_NoOvlp(Geom, Elem_db, MySD, SOL_Dist, SOL_Loc, SOL_Master)

    Call Init_BC(Geom, Params, Node_db)

    Call VecDuplicate(Sol_Dist, RHS, iErr)
    Call VecGhostGetLocalForm(RHS, RHS_Loc, iErr)

    Call VecDuplicate(Sol_Master, BC_Master, iErr)
    Call VecDuplicate(Sol_Dist, BC_Dist, iErr)
    Call VecGhostGetLocalForm(BC_Dist, BC_Loc, iErr)

    Call VecDuplicate(Sol_Master, F_Master, iErr)
    Call VecDuplicate(Sol_Dist, F_Dist, iErr)
    Call VecGhostGetLocalForm(F_Dist, F_Loc, iErr)

    Call MatCreateMPIAIJ(PETSC_COMM_WORLD, MySD%Num_Nodes,               &
    	 & MySD%Num_Nodes, Geom%Num_Nodes,                &
    	 & Geom%Num_Nodes, 24, PETSC_NULL_INTEGER, 24,         &
    	 & PETSC_NULL_INTEGER, MR, iErr)

    Call MatSetOption(MR, MAT_SYMMETRIC, iErr)
    Call MatSetFromOptions(MR, iErr)


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

    Call KSPSetTolerances(KSP_MR, Params%TolKSP,                              &
           & PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION,  &
           & PETSC_DEFAULT_INTEGER, iErr)
    Call KSPSetFromOptions(KSP_MR, iErr)

    Call PCGetType(PC_MR, PC_Type, iErr)
    If (PC_Type == PCBJACOBI) Then
       Call KSPSetUp(KSP_MR, iErr)
       Call PCBJacobiGetSubKSP(PC_MR, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
            & Sub_KSP_MR, iErr)
       Call KSPGetPC(Sub_KSP_MR, Sub_PC_MR, iErr)
       Call PCFactorSetZeroPivot(Sub_PC_MR, 1.0D-30, iErr)
       Call PCSetFromOptions(Sub_PC_MR, iErr)
    End If
  End Subroutine Init_KSP

  Subroutine Export()
    Real(Kind = Kr), Dimension(:), Pointer        :: SOL_Ptr
    
    Call VecScatterBegin(SOL_Dist, SOL_Master, INSERT_VALUES, SCATTER_FORWARD,&
         & MySD%ToMaster, iErr)
    Call VecScatterEnd(SOL_Dist, SOL_Master, INSERT_VALUES, SCATTER_FORWARD,  &
         & MySD%ToMaster, iErr)
    
    If (MEF90_MyRank == 0) Then
       Call VecGetArrayF90(SOL_Master, Sol_Ptr, iErr)
       Call Write_EXO_Result_Ptr_Nodes(Geom, 4, TimeStep, SOL_Ptr)
       Write(*,90) MinVal(Sol_Ptr), MaxVal(SOL_Ptr)
       Call VecRestoreArrayF90(SOL_Master, SOL_Ptr, iErr)
    End If

90  Format('Displacement Min / max: ', 2(ES10.3,' '))
  End Subroutine Export

  Subroutine Finalize()
    DeAllocate(Node_db, Elem_db)

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
    Type(Rupt_Params), Intent(IN)                    :: Params

    Type(Node2D), Dimension(:), Pointer              :: Node_db
    Integer                                          :: iN, iSet

    If (Geom%Numbering /= Numbering_PerNodes) Then
       Print*, '[ERROR] Init_BC not implemented for this numbering scheme'
       STOP
    End If

    Node_db(:)%BC = BC_Type_NONE
    Do iSet = 1, Geom%Num_Node_Sets
       Do iN = 1, Geom%Node_Set(iSet)%Num_Nodes
          Node_db(Geom%Node_Set(iSet)%Node_ID(iN))%BC =                       &
               & Params%BC_Type_Z(iSet)
        End Do
    End Do

  End Subroutine Init_BC


  Subroutine Assemb_Mat_Elast2DA(MR, Geom, Params, Elem_db, Node_db)
    Mat                                                 :: MR
    Type (EXO_Geom_Info)                                :: Geom
    Type (Rupt_Params)                                  :: Params
    Type (Node2D), Dimension(:), Pointer                :: Node_db 
    Type (Element2D_Scal), Dimension(:), Pointer        :: Elem_db 

    Integer(Kind = Ki)          :: Nb_Gauss, Nb_DoF
    Integer(Kind = Ki)          :: iSLEps, iSLSig
    Integer(Kind = Ki)          :: iSGEps, iSGSig
    Integer(Kind = Ki)          :: iE, iG, iELoc
    Real(Kind = Kr)             :: E, nu
    Real(Kind = Kr)             :: K2
    Integer(Kind = Ki)          :: iBlk
    Integer(Kind = Ki)          :: i

    PetscScalar                 :: value, one
    PetscScalar, Dimension(:,:), Pointer              :: MR_Elem
    Integer, Dimension(:), Pointer                    :: EXO_Indices
    PetscTruth                                        :: IsAssembled
    PetscLogDouble              :: SetTS, SetTF, SetT
    PetscLogDouble              :: GaussTS, GaussTF, GaussT

    Integer                     :: SetN 

    SetN = 0
    GaussT = 0.0

    one = 1.0_Kr

    Call MatAssembled(MR, IsAssembled, iErr)
    If (IsAssembled) Then
       Call MatZeroEntries(MR, iErr)
    End If
    
    ! Assembly of the Non BC terms 
    Do_iBlk: Do iBlk = 1, Geom%Num_Elem_Blks
       E  = Params%Young_Mod(iBlk)
       nu = Params%Poisson_Ratio(iBlk) 

!!! The isotropic Hooke's law is expressed as
!!! \sigma = K1 * trace(Epsilon) Id + 2*K2 * Epsilon - K3 Temp * Id
!!! K1, K2, K3 are computed in terms of E and nu
!!! in 3D, K1 = lambda, K2 = mu, K3 = E*alpha/(1-2nu) (= 3kappa alpha)
!!!       (alpha = therm exp coef).
!!! in 2D / plane stresses, the expressions are more complicated
!!!
!!! in 2DA
       K2 = E / (1.0_Kr + nu) * InvOf2
       
       Nb_DoF = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem
       Allocate (MR_Elem(Nb_DoF, Nb_DoF))
       Allocate (EXO_Indices(Nb_DoF))

       Do_iE: Do iELoc = 1, Geom%Elem_Blk(iBlk)%Num_Elems
          iE = Geom%Elem_Blk(iBlk)%ELem_ID(iELoc)
          If (.NOT. MySD%IsLocal_Elem(iE)) Then
             CYCLE
          End If
          EXO_Indices = Elem_db(iE)%Id_DoF-1
          Call AOApplicationToPETSc(MySD%EXO_AO, Nb_DoF, EXO_Indices, iErr)

          Call PetscGetTime(GaussTS, iErr)
          Call Init_Gauss_EXO(Elem_db, Node_db, Geom, GaussOrder, Elem=iE)
          Call PetscGetTime(GaussTF, iErr)
          GaussT = Gausst + GaussTF - GaussTS

          MR_Elem = 0.0_Kr
          Nb_Gauss = Elem_db(iE)%Nb_Gauss
          Do_iSLSig: Do iSLSig = 1, Nb_DoF
             ISGSig = Elem_db(iE)%ID_DoF(iSLSig)
                
             Do iG = 1, Nb_Gauss
             End Do
             Do_iSLEps: Do iSLEps = 1, Nb_DoF
                iSGEps = Elem_db(iE)%ID_DoF(iSLEps)
                Do iG = 1, Nb_Gauss
                   MR_Elem(iSLEps, iSLSig) = MR_Elem(iSLEps, iSLSig) +        &
                        & Elem_db(iE)%Gauss_C(iG) * (                         &
                        & (Elem_db(iE)%Grad_BF(iSLEps,iG) .DotP.              &
                        &  Elem_db(iE)%Grad_BF(iSLSig,iG)) ) * K2
                End Do
             End Do Do_iSLEps
          End Do Do_iSLSig

          Call MatSetValues(MR, Nb_DoF, EXO_Indices, Nb_DoF, EXO_Indices,     &
               & MR_Elem, ADD_VALUES, iErr)
          SetN = SetN + 1

          Call PetscGetTime(GaussTS, iErr)
          Call Destroy_Gauss_EXO(Elem_db, Elem=iE)
          Call PetscGetTime(GaussTF, iErr)
          GaussT = Gausst + GaussTF - GaussTS
       EndDo Do_iE

       DeAllocate (MR_Elem)
       DeAllocate (EXO_Indices)
    End Do Do_iBlk

    Call MatAssemblyBegin(MR, MAT_FLUSH_ASSEMBLY, iErr)  
    Call MatAssemblyEnd(MR, MAT_FLUSH_ASSEMBLY, iErr)  
    
    ! Assembly of the BC terms 
    Allocate(EXO_Indices(Geom%Num_Nodes))
    EXO_Indices = (/ (i ,i = 0, Geom%Num_Nodes - 1) /)
    Call AOApplicationToPETSc(MySD%EXO_AO, Geom%Num_Nodes,    &
         & EXO_Indices, iErr)
    Do_iBlk_BC: Do iBlk = 1, Geom%Num_Elem_Blks
       Nb_DoF = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem

       Do_iE_BC: Do iELoc = 1, Geom%Elem_blk(iBlk)%Num_Elems
          iE = Geom%Elem_blk(iBlk)%Elem_ID(iELoc)
          If (.NOT. MySD%IsLocal_Elem(iE)) Then
             CYCLE
          End If
          
          Do_iSLSig_BC: Do iSLSig = 1, Nb_DoF
             ISGSig = Elem_db(iE)%ID_DoF(iSLSig)
             Is_BC_BC: If ( (Node_db(iSGSig)%BC /= BC_Type_NONE) .AND.        &
                  & (MySD%IsLocal_Node(iSGSig)) )Then
                Call PetscGetTime(SetTS, iErr)
                Call MatSetValue(MR, EXO_Indices(iSGSig), EXO_Indices(iSGSig),&
                     &  VLV, INSERT_VALUES, iErr)
                Call PetscGetTime(SetTF, iErr)
                SetN = SetN+1
                SetT = SetT + SetTF - SetTS
             End If Is_BC_BC
          End Do Do_iSLSig_BC
       End Do Do_iE_BC
    End Do Do_iBlk_BC
    DeAllocate (EXO_Indices)

    Call MatAssemblyBegin(MR, MAT_FINAL_ASSEMBLY, iErr)
    Call MatAssemblyEnd(MR, MAT_FINAL_ASSEMBLY, iErr)  

#ifdef MEF90_TIMING
    If (MEF90_MyRank ==0) Then
       Write(CharBuffer,*) 'Total time in Init/Destroy_Gauss:         ',      &
            & GaussT, '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
       Write(CharBuffer,*) 'Total time in MatSetValue:                ',      &
            & SetT, '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
       Write(CharBuffer,*) 'Number of calls to MatSetValue:           ',      &
            & SetN, '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
    End If
#endif

  End Subroutine Assemb_Mat_Elast2DA

  Subroutine Assemb_RHS_Elast2DA(RHS, Geom, Params, Elem_db, Node_db, BC_Loc, F_Loc)
    Vec                                                 :: RHS
    Type (EXO_Geom_Info)                                :: Geom
    Type (Rupt_Params)                                  :: Params
    Type(Element2D_Scal), Dimension(:), Pointer         :: Elem_db
    Type(Node2D), Dimension(:), Pointer                 :: Node_db
    Vec                                                 :: BC_Loc
    Vec                                                 :: F_Loc

    Real(Kind = Kr), Dimension(:), Pointer              :: BC_Ptr
    Real(Kind = Kr), Dimension(:), Pointer              :: RHS_Ptr
    Real(Kind = Kr), Dimension(:), Pointer              :: F_Ptr
    

    Integer, Dimension(:), Pointer                      :: Loc_Indices
    Integer, Dimension(:), Pointer                      :: EXO_Indices

    Integer                                             :: Nb_Gauss, Nb_DoF
    Integer                                             :: iSL1, iSG1
    Integer                                             :: iSL2, iSG2
    Integer                                             :: iE, iELoc, iG
    Integer                                             :: iBlk
    
    Real(Kind = Kr)                                     :: K3
    Integer(Kind = Ki)                                  :: i, iS

    PetscLogDouble                                      :: TotTS, TotTF, TotT

    Real(Kind = Kr)                                     :: Tmp_Val


    Call PetscGetTime(TotTS, iErr)
    Call VecSet(RHS, 0.0_Kr, iErr)

    Allocate(EXO_Indices(Geom%Num_Nodes))
    EXO_Indices = (/ (i ,i = 0, Geom%Num_Nodes - 1) /)
    Call AOApplicationToPETSc(MySD%EXO_AO, Geom%Num_Nodes,      &
         & EXO_Indices, iErr)
    Allocate(Loc_Indices(Geom%Num_Nodes))
    Loc_Indices = (/ (i ,i = 0, Geom%Num_Nodes - 1) /)
    Call AOApplicationToPETSc(MySD%Loc_AO, Geom%Num_Nodes,      &
         & Loc_Indices, iErr)

    Call VecGetArrayF90(F_Loc, F_Ptr, iErr)

    Do_iBlk: Do iBlk = 1, Geom%Num_elem_blks

!!! The isotropic Hooke's law is expressed as
!!! \sigma = K1 * trace(Epsilon) Id + 2*K2 * Epsilon - K3 Temp * Id
!!! K1, K2, K3 are computed in terms of E and nu
!!! in 3D, K1 = lambda, K2 = mu, K3 = E*alpha/(1-2nu) (= 3kappa alpha)
!!!       (alpha = therm exp coef).
!!! in 2D / plane stresses, the expressions are more complicated
!!!

       Do_iE: Do iELoc = 1, Geom%Elem_blk(iBlk)%Num_Elems
          iE = Geom%Elem_blk(iBlk)%Elem_ID(iELoc)
          If (.NOT. MySD%IsLocal_Elem(iE)) Then
             CYCLE
          End If
          Call Init_Gauss_EXO(Elem_db, Node_db, Geom, GaussOrder, Elem=iE)

          Nb_Gauss = Elem_db(iE)%Nb_Gauss
          Do_iSL1: Do iSL1 = 1, Elem_db(iE)%Nb_DoF
             iSG1 = Elem_db(iE)%ID_DoF(iSL1)
             Tmp_Val = 0.0_Kr
             Do_iSL2: Do iSL2 = 1, Elem_db(iE)%Nb_DoF
             	iSG2 = Elem_db(iE)%ID_DoF(iSL2)
                Do_iG: Do iG = 1, Nb_Gauss
                   If (Params%Has_Force(iBlk)) Then
                      Tmp_Val = Tmp_Val + Elem_db(iE)%Gauss_C(iG) *           &
                   		   & F_Ptr(Loc_Indices(iSG2)+1) *                     &
                           & Elem_db(iE)%BF(iSL1, iG) *                       &
                           & Elem_db(iE)%BF(iSL2,iG) 
                   End If
                End Do Do_iG
             End Do Do_ISL2

             Call VecSetValue(RHS, EXO_Indices(iSG1), Tmp_Val, ADD_VALUES,   &
                  & iErr)
          End Do Do_iSL1
          Call Destroy_Gauss_EXO(Elem_db, Elem=iE)
       End Do Do_iE
    End Do Do_iBlk
    Call VecAssemblyBegin(RHS, iErr)
    DeAllocate(EXO_Indices)
    DeAllocate(Loc_Indices)
    Call VecAssemblyEnd(RHS, iErr)
    Call VecRestoreArrayF90(F_Loc, F_Ptr, iErr)
    
    !!! BC Part, using VLV
    !!! There are two ways to do that:
    !!! Use VecSetValue and EXO_Indices or UseVecGetArray and Loc_Indices
    !!! I chose the second way for no special reason.
    Allocate(Loc_Indices(Geom%Num_Nodes))
    Loc_Indices = (/ (i ,i = 0, Geom%Num_Nodes - 1) /)
    Call AOApplicationToPETSc(MySD%Loc_AO, Geom%Num_Nodes,    &
         & Loc_Indices, iErr)

    Call VecGetArrayF90(BC_Loc, BC_Ptr, iErr)
    Call VecGetArrayF90(RHS, RHS_Ptr, iErr)

    Do_iS: Do iS = 1, Geom%Num_Nodes
       Is_BC: If ((Node_db(iS)%BC /= BC_Type_NONE) .AND.            &
            & (MySD%IsLocal_Node(iS)) ) Then
          RHS_Ptr(Loc_Indices(iS)+1) = BC_Ptr(Loc_Indices(iS)+1)*VLV
       End If Is_BC
    End Do Do_iS


    Call PetscGetTime(TotTF, iErr)
    TotT = TotTF - TotTS
    DeAllocate(Loc_Indices)
    Call VecRestoreArrayF90(BC_Loc, BC_Ptr, iErr)
    Call VecRestoreArrayF90(RHS, RHS_Ptr, iErr)
  End Subroutine Assemb_RHS_Elast2DA

!!!  Subroutine Calc_Ener(DISP_Loc, Geom, Params, Elem_db, Node_db, SD, F_Loc,   &
!!!       & Ener) 
!!!
!!!#ifdef PB_2D
!!!    Type(MatS2D), Dimension(:), Pointer                 :: Sigma
!!!    Type(MatS2D), Dimension(:), Pointer                 :: Epsilon
!!!    Type(Element2D_Elast), Dimension(:), Pointer        :: Elem_db 
!!!    Type(Node2D), Dimension(:), Pointer                 :: Node_db 
!!!#else
!!!    Type(MatS3D), Dimension(:), Pointer                 :: Sigma
!!!    Type(MatS3D), Dimension(:), Pointer                 :: Epsilon
!!!    Type(Element3D_Elast), Dimension(:), Pointer        :: Elem_db 
!!!    Type(Node3D), Dimension(:), Pointer                 :: Node_db 
!!!#endif
!!!    Type (EXO_Geom_Info)                                :: Geom
!!!    Type (Rupt_Params)                                  :: Params
!!!    Type (SD_Info)                                      :: SD
!!!    Vec                                                 :: DISP_Loc
!!!	Vec                                                 :: F_Loc
!!!	    
!!!    Real(Kind = Kr)                                     :: Ener
!!!    Real(Kind = Kr)                                     :: MyEner
!!!
!!!    Integer(Kind = Ki)          :: Nb_Gauss, Nb_DoF
!!!    Integer(Kind = Ki)          :: iSLEps, iSLSig
!!!    Integer(Kind = Ki)          :: iSGEps, iSGSig
!!!    Integer(Kind = Ki)          :: iE, iG, iELoc
!!!    Real(Kind = Kr)             :: K1, K2, K3
!!!    Real(Kind = Kr)             :: E, nu
!!!    Integer(Kind = Ki)          :: iBlk
!!!    Integer(Kind = Ki)          :: i
!!!    Integer, Dimension(:), Pointer                    :: Loc_Indices
!!!    Real(Kind = Kr), Dimension(:), Pointer            :: DISP_Ptr
!!!    Real(Kind = Kr), Dimension(:), Pointer            :: F_Ptr
!!!    
!!!    Ener = 0.0_Kr
!!!    MyEner = 0.0_Kr
!!!    Allocate(Loc_Indices(Geom%Num_Nodes * Geom%Num_Dim))
!!!    Loc_Indices = (/ (i ,i = 0, Geom%Num_Nodes * Geom%Num_Dim -1) /)
!!!    Call AOApplicationToPETSc(MySD_Vect%Loc_AO, Geom%Num_Nodes * Geom%Num_Dim,    &
!!!         & Loc_Indices, iErr)
!!!
!!!    Call VecGetArrayF90(DISP_Loc, DISP_Ptr, iErr)
!!!    Call VecGetArrayF90(F_Loc, F_Ptr, iErr)
!!!
!!!    Do_iBlk: Do iBlk = 1, Geom%Num_Elem_Blks
!!!       E  = Params%Young_Mod(iBlk)
!!!       nu = Params%Poisson_Ratio(iBlk) 
!!!
!!!!!! The isotropic Hooke's law is expressed as
!!!!!! \sigma = K1 * trace(Epsilon) Id + 2*K2 * Epsilon - K3 Temp * Id
!!!!!! K1, K2, K3 are computed in terms of E and nu
!!!!!! in 3D, K1 = lambda, K2 = mu, K3 = E*alpha/(1-2nu) (= 3kappa alpha)
!!!!!!       (alpha = therm exp coef).
!!!!!! in 2D / plane stresses, the expressions are more complicated
!!!!!!
!!!#ifdef PB_2D
!!!       K1 = E * nu / (1.0_Kr - nu**2)
!!!       K2 = E / (1.0_Kr + nu) * InvOf2
!!!       K3 = Params%Therm_Exp(iBlk) * Params%Young_Mod(iBlk) /                 &
!!!            (1.0_Kr - Params%Poisson_Ratio(iBlk) )
!!!#else
!!!       K1 = E * nu / (1.0_Kr - 2.0_Kr * nu) / ( 1.0_Kr + nu)
!!!       K2 = E / (1.0_Kr + nu) * InvOf2
!!!       K3 = Params%Therm_Exp(iBlk) * Params%Young_Mod(iBlk) /                 &
!!!            (1.0_Kr - 2.0_Kr * Params%Poisson_Ratio(iBlk) )
!!!#endif
!!!
!!!       Nb_DoF = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem * Geom%Num_Dim
!!!       Do_iE: Do iELoc = 1, Geom%Elem_Blk(iBlk)%Num_Elems
!!!          iE = Geom%Elem_Blk(iBlk)%ELem_ID(iELoc)
!!!          If (.NOT. MySD_Vect%IsLocal_Elem(iE)) Then
!!!             CYCLE
!!!          End If
!!!
!!!          Call Init_Gauss_EXO(Elem_db, Node_db, Geom, GaussOrder, Elem=iE)
!!!
!!!          Nb_Gauss = Elem_db(iE)%Nb_Gauss
!!!          Allocate (Sigma(Nb_Gauss))
!!!          Do_iSLSig: Do iSLSig = 1, Nb_DoF
!!!             ISGSig = Elem_db(iE)%ID_DoF(iSLSig)
!!!
!!!             Do iG = 1, Nb_Gauss
!!!                Sigma(iG) = 2.0_Kr * K2 * Elem_db(iE)%GradS_BF(iSLSig,iG)
!!!                Sigma(iG)%XX = Sigma(iG)%XX +&
!!!                     & K1 * Trace(Elem_db(iE)%GradS_BF(iSLSig,iG))
!!!                Sigma(iG)%YY = Sigma(iG)%YY +&
!!!                     & K1 * Trace(Elem_db(iE)%GradS_BF(iSLSig,iG))
!!!#ifdef PB_3D
!!!                Sigma(iG)%ZZ = Sigma(iG)%ZZ +&
!!!                     & K1 * Trace(Elem_db(iE)%GradS_BF(iSLSig,iG))
!!!#endif
!!!             End Do
!!!             Do_iSLEps: Do iSLEps = 1, Nb_DoF
!!!                iSGEps = Elem_db(iE)%ID_DoF(iSLEps)
!!!                Do iG = 1, Nb_Gauss
!!!                   MyEner = MyEner + Elem_db(iE)%Gauss_C(iG) *                &
!!!                        & DISP_Ptr(Loc_Indices(iSGEps) + 1) *                 &
!!!                        & DISP_Ptr(Loc_Indices(iSGSig) + 1) *                 &
!!!                        & (Elem_db(iE)%GradS_BF(iSLEps,iG) .DotP. Sigma(iG))* &
!!!                        & InvOf2
!!!                End Do
!!!                If (Params%Has_Force(iBlk)) Then
!!!                   Do iG = 1, Nb_Gauss
!!!                      MyEner = MyEner - Elem_db(iE)%Gauss_C(iG) *             &
!!!                           & DISP_Ptr(Loc_Indices(iSGSig)+1) *                &
!!!                           & DISP_Ptr(Loc_Indices(iSGEps)+1) *                &
!!!                           & ( Elem_db(iE)%BF(iSLEps, iG) .DotP.              &
!!!                           &   Elem_db(iE)%BF(iSLSig, iG) )
!!!                   End Do
!!!                End If
!!!             End Do Do_iSLEps
!!!          End Do Do_iSLSig
!!!          DeAllocate(Sigma)
!!!
!!!
!!!          Call Destroy_Gauss_EXO(Elem_db, Elem=iE)
!!!       EndDo Do_iE
!!!
!!!    End Do Do_iBlk
!!!    DeAllocate (Loc_Indices)
!!!    Call VecRestoreArrayF90(DISP_Loc, DISP_Ptr, iErr)
!!!    Call VecRestoreArrayF90(F_Loc, F_Ptr, iErr)
!!!
!!!    Call MPI_AllReduce(MyEner, Ener, 1, MPI_DOUBLE_PRECISION, MPI_SUM,        &
!!!         & PETSC_COMM_WORLD, iErr)
!!!  End Subroutine Calc_Ener

End Module m_Elast2DA_Proc
