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
    End If
    Read (*,100) Params%Sim_Str

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

    Call Init_BC(Geom, Params, Node_db)

    Call VecCreateSeq(PETSC_COMM_WORLD, Geom%Num_Nodes, Sol, iErr)

    Call VecDuplicate(Sol, BC, iErr)
    Call VecDuplicate(Sol, F, iErr)
    Call VecDuplicate(Sol, RHS, iErr)

    Call MatCreateSeqAIJ(PETSC_COMM_WORLD, Geom%Num_Nodes,                &
    	 & Geom%Num_Nodes, 24, PETSC_NULL_INTEGER, MR, iErr)

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
    Call PCSetType(PC_MR, PCILU, iErr)
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
       Call PCFactorSetZeroPivot(Sub_PC_MR, 1.0D-20, iErr)
       Call PCSetFromOptions(Sub_PC_MR, iErr)
    End If
  End Subroutine Init_KSP

  Subroutine Export()
    Real(Kind = Kr), Dimension(:), Pointer        :: SOL_Ptr
    
    Call VecGetArrayF90(SOL, Sol_Ptr, iErr)
    Call Write_EXO_Result_Ptr_Nodes(Geom, 4, TimeStep, SOL_Ptr)
    Write(*,90) MinVal(Sol_Ptr), MaxVal(SOL_Ptr)
    Call VecRestoreArrayF90(SOL, SOL_Ptr, iErr)

90  Format('Displacement Min / max: ', 2(ES10.3,' '))
  End Subroutine Export

  Subroutine Finalize()
    DeAllocate(Node_db, Elem_db)

    Call MatDestroy(MR, iErr)
    Call VecDestroy(RHS, iErr)
    Call VecDestroy(SOL, iErr)
    Call VecDestroy(BC, iErr)
    Call VecDestroy(F, iErr)

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

       Do_iE: Do iELoc = 1, Geom%Elem_Blk(iBlk)%Num_Elems
          iE = Geom%Elem_Blk(iBlk)%ELem_ID(iELoc)

          Call Init_Gauss_EXO(Elem_db, Node_db, Geom, GaussOrder, Elem=iE)

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
                        &  Elem_db(iE)%Grad_BF(iSLSig,iG)) )
                End Do
             End Do Do_iSLEps
          End Do Do_iSLSig

          Call MatSetValues(MR, Nb_DoF, Elem_db(iE)%ID_DoF - 1, Nb_DoF, Elem_db(iE)%ID_DoF - 1,     &
               & MR_Elem, ADD_VALUES, iErr)

          Call Destroy_Gauss_EXO(Elem_db, Elem=iE)
       EndDo Do_iE

       DeAllocate (MR_Elem)
    End Do Do_iBlk

    Call MatAssemblyBegin(MR, MAT_FLUSH_ASSEMBLY, iErr)  
    Call MatAssemblyEnd(MR, MAT_FLUSH_ASSEMBLY, iErr)  
    
    ! Assembly of the BC terms 
    Do_iBlk_BC: Do iBlk = 1, Geom%Num_Elem_Blks
       Nb_DoF = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem

       Do_iE_BC: Do iELoc = 1, Geom%Elem_blk(iBlk)%Num_Elems
          iE = Geom%Elem_blk(iBlk)%Elem_ID(iELoc)
          
          Do_iSLSig_BC: Do iSLSig = 1, Nb_DoF
             ISGSig = Elem_db(iE)%ID_DoF(iSLSig)
             Is_BC_BC: If ( Node_db(iSGSig)%BC /= BC_Type_NONE ) Then
                Call MatSetValue(MR, iSGSig - 1, iSGSig - 1,&
                     &  VLV, INSERT_VALUES, iErr)
             End If Is_BC_BC
          End Do Do_iSLSig_BC
       End Do Do_iE_BC
    End Do Do_iBlk_BC

    Call MatAssemblyBegin(MR, MAT_FINAL_ASSEMBLY, iErr)
    Call MatAssemblyEnd(MR, MAT_FINAL_ASSEMBLY, iErr)  

  End Subroutine Assemb_Mat_Elast2DA

  Subroutine Assemb_RHS_Elast2DA(RHS, Geom, Params, Elem_db, Node_db, BC, F)
    Vec                                                 :: RHS
    Type (EXO_Geom_Info)                                :: Geom
    Type (Rupt_Params)                                  :: Params
    Type(Element2D_Scal), Dimension(:), Pointer         :: Elem_db
    Type(Node2D), Dimension(:), Pointer                 :: Node_db
    Vec                                                 :: BC, F

    Real(Kind = Kr), Dimension(:), Pointer              :: BC_Ptr
    Real(Kind = Kr), Dimension(:), Pointer              :: RHS_Ptr
    Real(Kind = Kr), Dimension(:), Pointer              :: F_Ptr
    

    Integer                                             :: Nb_Gauss, Nb_DoF
    Integer                                             :: iSL1, iSG1
    Integer                                             :: iSL2, iSG2
    Integer                                             :: iE, iELoc, iG
    Integer                                             :: iBlk
    
    Integer(Kind = Ki)                                  :: i, iS

    PetscLogDouble                                      :: TotTS, TotTF, TotT

    Real(Kind = Kr)                                     :: Tmp_Val


    Call VecSet(RHS, 0.0_Kr, iErr)
    Call VecGetArrayF90(F, F_Ptr, iErr)

    Do_iBlk: Do iBlk = 1, Geom%Num_elem_blks
       Do_iE: Do iELoc = 1, Geom%Elem_blk(iBlk)%Num_Elems
          iE = Geom%Elem_blk(iBlk)%Elem_ID(iELoc)
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
                   		   & F_Ptr(iSG2) *                     &
                           & Elem_db(iE)%BF(iSL1, iG) *                       &
                           & Elem_db(iE)%BF(iSL2,iG) 
                   End If
                End Do Do_iG
             End Do Do_ISL2

             Call VecSetValue(RHS, iSG1 - 1, Tmp_Val, ADD_VALUES, iErr)
          End Do Do_iSL1
          Call Destroy_Gauss_EXO(Elem_db, Elem=iE)
       End Do Do_iE
    End Do Do_iBlk
    Call VecAssemblyBegin(RHS, iErr)
    Call VecAssemblyEnd(RHS, iErr)
    Call VecRestoreArrayF90(F, F_Ptr, iErr)
    
    !!! BC Part, using VLV
    !!! There are two ways to do that:
    !!! Use VecSetValue and EXO_Indices or UseVecGetArray and Loc_Indices
    !!! I chose the second way for no special reason.

    Call VecGetArrayF90(BC, BC_Ptr, iErr)
    Call VecGetArrayF90(RHS, RHS_Ptr, iErr)

    Do_iS: Do iS = 1, Geom%Num_Nodes
       Is_BC: If ( Node_db(iS)%BC /= BC_Type_NONE ) Then
          RHS_Ptr(iS) = BC_Ptr(iS) * VLV
       End If Is_BC
    End Do Do_iS

    Call VecRestoreArrayF90(BC, BC_Ptr, iErr)
    Call VecRestoreArrayF90(RHS, RHS_Ptr, iErr)
  End Subroutine Assemb_RHS_Elast2DA

End Module m_Elast2DA_Proc