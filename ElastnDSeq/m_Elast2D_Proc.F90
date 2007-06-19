#ifdef PB_2D
Module m_Elast2D_Proc
#else
Module m_Elast3D_Proc
#endif
  Use m_MEF90
  Use m_Rupt_Struct

#ifdef PB_2D
  Use m_Elast2D_Vars
#else
  Use m_Elast3D_Vars
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
#include "include/finclude/petscao.h"

  Public :: Init
  Public :: Init_KSP
  Public :: Finalize
  Public :: Export

  Public :: Assemb_Mat_Elast
  Public :: Assemb_RHS_Elast
  Public :: Calc_Ener
Contains


  Subroutine Init()
    Integer                            :: i, iS, iE

    PetscTruth                :: Has_Sim_Str

    Call MEF90_Initialize()

    Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-f', Params%Sim_Str, Has_Sim_Str, iErr)
    
    If (.NOT. Has_Sim_Str) Then
       Write(CharBuffer, 100) 'Simulation name: \n'c
       Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
       Read(*,100) Params%Sim_Str
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

    Write(CharBuffer,*) 'Number of nodes:                          ',  Geom%Num_Nodes, '\n'c
    Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
    Write(CharBuffer,*) 'Number of elements:                       ', Geom%Num_Elems, '\n'c
    Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)

    Call Init_BC(Geom, Params, Node_db)

    Call VecCreateSeq(PETSC_COMM_WORLD, Geom%Num_Nodes * Geom%Num_Dim, Sol, iErr)
    Call VecDuplicate(Sol, RHS, iErr)
    Call VecDuplicate(Sol, F, iErr)
    Call VecDuplicate(Sol, BC, iErr)
    Call VecCreateSeq(PETSC_COMM_WORLD, Geom%Num_Nodes, Temp, iErr)

    Call MatCreateSeqAIJ(PETSC_COMM_WORLD, Geom%Num_Nodes * Geom%Num_Dim, Geom%Num_Nodes * Geom%Num_Dim, 60, PETSC_NULL_INTEGER,   &
                         MR, iErr)

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

    Call VecGetArrayF90(SOL, Sol_Ptr, iErr)
    Call Write_EXO_Result_Ptr_Nodes(Geom, 2, TimeStep, SOL_Ptr)
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

#ifdef PB_2D
    Type(Node2D), Dimension(:), Pointer              :: Node_db
#else
    Type(Node3D), Dimension(:), Pointer              :: Node_db
#endif
    Integer                                          :: iN, iSet

    Node_db(:)%BC = BC_Type_NONE
    Do iSet = 1, Geom%Num_Node_Sets
       Do iN = 1, Geom%Node_Set(iSet)%Num_Nodes
          Node_db(Geom%Num_Dim * (Geom%Node_Set(iSet)%Node_ID(iN)-1)+1)%BC  = Params%BC_Type_X(iSet)
          Node_db(Geom%Num_Dim * (Geom%Node_Set(iSet)%Node_ID(iN)-1)+2)%BC  = Params%BC_Type_Y(iSet)
#ifdef PB_3D
          Node_db(Geom%Num_Dim * Geom%Node_Set(iSet)%Node_ID(iN) )%BC       = Params%BC_Type_Z(iSet)
#endif
       End Do
    End Do
  End Subroutine Init_BC

  Subroutine Assemb_Mat_Elast(MR, Geom, Params, Elem_db, Node_db)
    Mat                                                 :: MR
    Type (EXO_Geom_Info)                                :: Geom
    Type (Rupt_Params)                                  :: Params
#ifdef PB_2D
    Type (Node2D), Dimension(:), Pointer                :: Node_db 
    Type (Element2D_Elast), Dimension(:), Pointer       :: Elem_db 
    Type (MatS2D)                                       :: Sigma, Epsilon
    Type (Tens4OS2D)                                    :: HookeLaw
#else
    Type (Node3D), Dimension(:), Pointer                :: Node_db 
    Type (Element3D_Elast), Dimension(:), Pointer       :: Elem_db 
    Type (MatS3D)                                       :: Sigma, Epsilon
    Type (Tens4OS3D)                                    :: HookeLaw
#endif

    Integer                                             :: Nb_Gauss, Nb_DoF
    Integer                                             :: iSL1, iSL2
    Integer                                             :: iSG1, iSG2
    Integer                                             :: iE, iG, iELoc
    Real(Kind = Kr)                                     :: E, nu
    Integer                                             :: iBlk
    Integer                                             :: i

    PetscScalar, Dimension(:,:), Pointer                :: MR_Elem
    PetscTruth                                          :: IsAssembled

    Call MatAssembled(MR, IsAssembled, iErr)
    If (IsAssembled) Then
       Call MatZeroEntries(MR, iErr)
    End If
    
    ! Assembly of the Non BC terms 
    Do_iBlk: Do iBlk = 1, Geom%Num_Elem_Blks
       E  = Params%Young_Mod(iBlk)
       nu = Params%Poisson_Ratio(iBlk) 

#ifdef PB_2D
       Call GenHL_Iso2D_EnuPlaneStress(E, nu, HookeLaw)
#else
       Call GenHL_Iso3D_Enu(E, nu, HookeLaw)
#endif

       Nb_DoF = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem * Geom%Num_Dim
       Allocate (MR_Elem(Nb_DoF, Nb_DoF))

       Do_iE: Do iELoc = 1, Geom%Elem_Blk(iBlk)%Num_Elems
          iE = Geom%Elem_Blk(iBlk)%ELem_ID(iELoc)
          Call Init_Gauss_EXO(Elem_db, Node_db, Geom, GaussOrder, Elem=iE)

          MR_Elem = 0.0_Kr
          Nb_Gauss = Elem_db(iE)%Nb_Gauss
          Do_iG: Do iG = 1, Nb_Gauss
             Do_iSL1: Do iSL1 = 1, Nb_DoF
                ISG1 = Elem_db(iE)%ID_DoF(iSL1)
                Epsilon = Elem_db(iE)%GradS_BF(iSL1,iG)
                Sigma = HookeLaw * Epsilon

                Do_iSL2: Do iSL2 = 1, Nb_DoF
                   iSG2 = Elem_db(iE)%ID_DoF(iSL2)
                   Epsilon = Elem_db(iE)%GradS_BF(iSL2,iG)
                   MR_Elem(iSL2, iSL1) = MR_Elem(iSL2, iSL1) + Elem_db(iE)%Gauss_C(iG) * (Sigma .DotP. Epsilon)
                End Do Do_iSL2
             End Do Do_iSL1
          End Do Do_iG

          Call MatSetValues(MR, Nb_DoF, Elem_db(iE)%ID_DoF-1, Nb_DoF, Elem_db(iE)%ID_DoF-1, MR_Elem, ADD_VALUES, iErr)
          Call Destroy_Gauss_EXO(Elem_db, Elem=iE)
       EndDo Do_iE

       DeAllocate (MR_Elem)
    End Do Do_iBlk

    Call MatAssemblyBegin(MR, MAT_FLUSH_ASSEMBLY, iErr)  
    Call MatAssemblyEnd(MR, MAT_FLUSH_ASSEMBLY, iErr)  
    
    ! Assembly of the BC terms 
    Do_iBlk_BC: Do iBlk = 1, Geom%Num_Elem_Blks
       Nb_DoF = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem * Geom%Num_Dim

       Do_iE_BC: Do iELoc = 1, Geom%Elem_blk(iBlk)%Num_Elems
          iE = Geom%Elem_blk(iBlk)%Elem_ID(iELoc)
          
          Do_iSL1_BC: Do iSL1 = 1, Nb_DoF
             ISG1 = Elem_db(iE)%ID_DoF(iSL1)
             Is_BC_BC: If (Node_db(iSG1)%BC /= BC_Type_NONE) Then
                Call MatSetValue(MR, iSG1-1, iSG1-1, VLV, INSERT_VALUES, iErr)
             End If Is_BC_BC
          End Do Do_iSL1_BC
       End Do Do_iE_BC
    End Do Do_iBlk_BC

    Call MatAssemblyBegin(MR, MAT_FINAL_ASSEMBLY, iErr)
    Call MatAssemblyEnd(MR, MAT_FINAL_ASSEMBLY, iErr)  
  End Subroutine Assemb_Mat_Elast

  Subroutine Assemb_RHS_Elast(RHS, Geom, Params, Elem_Vect, Node_Vect, Elem_Scal, Node_Scal, BC_Loc, F_Loc, Temp_Loc)
    Vec                                                 :: RHS
    Type (EXO_Geom_Info)                                :: Geom
    Type (Rupt_Params)                                  :: Params
#ifdef PB_2D
    Type(Element2D_Elast), Dimension(:), Pointer        :: Elem_Vect
    Type(Node2D), Dimension(:), Pointer                 :: Node_Vect
    Type(Element2D_Scal), Dimension(:), Pointer         :: Elem_Scal
    Type(Node2D), Dimension(:), Pointer                 :: Node_Scal
    Type(MatS2D)                                        :: ThetaId, AThetaId
    Type(Vect2D)                                        :: F
    Type(Tens4OS2D)                                     :: HookeLaw
#else
    Type(Element3D_Elast), Dimension(:), Pointer        :: Elem_Vect
    Type(Node3D), Dimension(:), Pointer                 :: Node_Vect
    Type(Element3D_Scal), Dimension(:), Pointer         :: Elem_Scal
    Type(Node3D), Dimension(:), Pointer                 :: Node_Scal
    Type(MatS3D)                                        :: ThetaId
    Type(Vect3D)                                        :: F
    Type(Tens4OS3D)                                     :: HookeLaw
#endif
    Vec                                                 :: BC_Loc
    Vec                                                 :: F_Loc
    Vec                                                 :: Temp_Loc

    Real(Kind = Kr), Dimension(:), Pointer              :: BC_Ptr
    Real(Kind = Kr), Dimension(:), Pointer              :: RHS_Ptr
    Real(Kind = Kr), Dimension(:), Pointer              :: F_Ptr
    Real(Kind = Kr), Dimension(:), Pointer              :: Temp_Ptr

    Integer                                             :: Nb_Gauss, Nb_DoF_Vect
    Integer                                             :: iSL, iSG
    Integer                                             :: iE, iELoc, iG
    Integer                                             :: iBlk
    
    Real(Kind = Kr)                                     :: E, nu
    Integer                                             :: i, iS

    Real(Kind = Kr), Dimension(:), Pointer              :: RHS_Elem


    Call VecGetArrayF90(F_Loc, F_Ptr, iErr)
    Call VecGetArrayF90(Temp_Loc, Temp_Ptr, iErr)


    Call VecSet(RHS, 0.0_Kr, iErr)
    Do_iBlk: Do iBlk = 1, Geom%Num_elem_blks
       Nb_DoF_Vect = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem * Geom%Num_Dim
       Allocate (RHS_Elem(Nb_DoF_Vect))

       E  = Params%Young_Mod(iBlk)
       nu = Params%Poisson_Ratio(iBlk) 

#ifdef PB_2D
       Call GenHL_ISO2D_EnuPlaneStress(E, nu, HookeLaw)
#else
       Call GenHL_ISO3D_Enu(E, nu, HookeLaw)
#endif


       Do_iE: Do iELoc = 1, Geom%Elem_blk(iBlk)%Num_Elems
          iE = Geom%Elem_blk(iBlk)%Elem_ID(iELoc)

          Call Init_Gauss_EXO(Elem_Vect, Node_Vect, Geom, GaussOrder, Elem=iE)
          Call Init_Gauss_EXO(Elem_Scal, Node_Scal, Geom, GaussOrder, Elem=iE)

          RHS_Elem = 0.0_Kr
          Do_iG: Do iG = 1, Elem_Vect(iE)%Nb_Gauss
          
             ThetaId = 0.0_Kr
             Do_iSL1: Do iSL = 1, Elem_Scal(iE)%Nb_DoF
                iSG = Elem_Scal(iE)%ID_DoF(iSL)
                ThetaId%XX = ThetaId%XX + Params%Therm_Exp(iBlk) * Temp_Ptr(iSG) * Elem_Scal(iE)%BF(iSL, iG)
                ThetaId%YY = ThetaId%YY + Params%Therm_Exp(iBlk) * Temp_Ptr(iSG) * Elem_Scal(iE)%BF(iSL, iG)
#ifdef PB_3D
                ThetaId%ZZ = ThetaId%ZZ + Params%Therm_Exp(iBlk) * Temp_Ptr(iSG) * Elem_Scal(iE)%BF(iSL, iG)
#endif
             End Do Do_iSL1

             F = 0.0_Kr
             Do_iSL2: Do iSL = 1, Elem_Vect(iE)%Nb_DoF
                iSG = Elem_Vect(iE)%ID_DoF(iSL)
                F = F + F_Ptr(iSG) * Elem_Vect(iE)%BF(iSL, iG)
             End Do Do_iSL2
             
             Do_iSL3: Do iSL = 1, Elem_Vect(iE)%Nb_DoF
                RHS_Elem(iSL) = RHS_Elem(iSL) + Elem_Vect(iE)%Gauss_C(iG) *                                                        &
                                ( (F .DotP. Elem_Vect(iE)%BF(iSL, iG)) +                                                           &
                                  (ThetaId .DotP.  (HookeLaw * Elem_Vect(iE)%GradS_BF(iSL, iG))) ) 
             End Do Do_iSL3
          End Do Do_iG
          Call VecSetValues(RHS, Elem_Vect(iE)%Nb_DoF, Elem_Vect(iE)%ID_DoF-1, RHS_Elem, ADD_VALUES, iErr)
         
          Call Destroy_Gauss_EXO(Elem_Vect, Elem=iE)
          Call Destroy_Gauss_EXO(Elem_Scal, Elem=iE)
       End Do Do_iE
       DeAllocate(RHS_Elem)
    End Do Do_iBlk

    Call VecAssemblyBegin(RHS, iErr)
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
       Is_BC: If (Node_Vect(iS)%BC /= BC_Type_NONE) Then
          RHS_Ptr(iS) = BC_Ptr(iS)*VLV
       End If Is_BC
    End Do Do_iS

    Call VecRestoreArrayF90(BC_Loc, BC_Ptr, iErr)
    Call VecRestoreArrayF90(RHS, RHS_Ptr, iErr)
  End Subroutine Assemb_RHS_Elast
  
  Subroutine Calc_Ener(DISP_Loc, Geom, Params, Elem_db_Vec, Node_db_Vec, Elem_db_Scal, Node_db_Scal, F_Loc, Temp_Loc, Ener) 

#ifdef PB_2D
    Type(MatS2D)                                        :: Sigma
    Type(MatS2D)                                        :: Epsilon
    Type(Element2D_Elast), Dimension(:), Pointer        :: Elem_db_Vec 
    Type(Node2D), Dimension(:), Pointer                 :: Node_db_Vec
    Type(Vect2D)                                        :: F, U
    Type(Element2D_Scal), Dimension(:), Pointer         :: Elem_db_Scal 
    Type(Node2D), Dimension(:), Pointer                 :: Node_db_Scal
    Type(Tens4OS2D)                                     :: HookeLaw
#else
    Type(MatS3D)                                        :: Sigma
    Type(MatS3D)                                        :: Epsilon
    Type(Element3D_Elast), Dimension(:), Pointer        :: Elem_db_Vec 
    Type(Node3D), Dimension(:), Pointer                 :: Node_db_Vec 
    Type(Vect3D)                                        :: F, U
    Type(Element3D_Scal), Dimension(:), Pointer         :: Elem_db_Scal 
    Type(Node3D), Dimension(:), Pointer                 :: Node_db_Scal
    Type(Tens4OS3D)                                     :: HookeLaw
#endif
    Type (EXO_Geom_Info)                                :: Geom
    Type (Rupt_Params)                                  :: Params
    
    Vec                                                 :: DISP_Loc
    Vec                                                 :: F_Loc
    Vec                                                 :: Temp_Loc
	    
    Real(Kind = Kr)                                     :: Ener

    Integer                                             :: Nb_Gauss, Nb_DoF_Vec, Nb_DoF_Scal
    Integer                                             :: iSL, iSG
    Integer                                             :: iE, iG, iELoc
    Real(Kind = Kr)                                     :: E, nu
    Integer(Kind = Ki)                                  :: iBlk
    Integer                                             :: i
    Real(Kind = Kr), Dimension(:), Pointer              :: DISP_Ptr
    Real(Kind = Kr), Dimension(:), Pointer              :: F_Ptr
    Real(Kind = Kr), Dimension(:), Pointer              :: Temp_Ptr
    
    Ener = 0.0_Kr

    Call VecGetArrayF90(DISP_Loc, DISP_Ptr, iErr)
    Call VecGetArrayF90(F_Loc, F_Ptr, iErr)
    Call VecGetArrayF90(Temp_Loc, Temp_Ptr, iErr)

    Do_iBlk: Do iBlk = 1, Geom%Num_Elem_Blks
       E  = Params%Young_Mod(iBlk)
       nu = Params%Poisson_Ratio(iBlk) 

#ifdef PB_2D
       Call GenHL_ISO2D_EnuPlaneStress(E, nu, HookeLaw)
#else
       Call GenHL_ISO3D_Enu(E, nu, HookeLaw)
#endif

       Nb_DoF_Vec  = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem * Geom%Num_Dim
       Nb_DoF_Scal = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem
       Do_iE: Do iELoc = 1, Geom%Elem_Blk(iBlk)%Num_Elems
          iE = Geom%Elem_Blk(iBlk)%Elem_ID(iELoc)

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
                Epsilon = Epsilon + Elem_db_Vec(iE)%GradS_BF(iSL,iG) * Disp_Ptr(iSG)
             End Do Do_iSL1

!!! Thermal stuff
             Do_iSL2: Do iSL = 1, Nb_DoF_Scal
                iSG = Elem_db_Scal(iE)%ID_DoF(iSL)
                Epsilon%XX = Epsilon%XX - Params%Therm_Exp(iBlk) * Elem_db_Scal(iE)%BF(iSL,iG) * Temp_Ptr(iSG)
                Epsilon%YY = Epsilon%YY - Params%Therm_Exp(iBlk) * Elem_db_Scal(iE)%BF(iSL,iG) * Temp_Ptr(iSG)
#ifdef PB_3D
                Epsilon%ZZ = Epsilon%ZZ - Params%Therm_Exp(iBlk) * Elem_db_Scal(iE)%BF(iSL,iG) * Temp_Ptr(iSG)
#endif
             End Do Do_iSL2
             Sigma = HookeLaw * Epsilon

             Ener = Ener + Elem_db(iE)%Gauss_C(iG) * ( Sigma .DotP. Epsilon ) * .5_Kr
            
!!! Forces stuff
             If (Params%Has_Force(iBlk)) Then
                Do_iSL3: Do iSL = 1, Nb_DoF_Vec
                   iSG = Elem_db(iE)%ID_DoF(iSL)
                   U = U + DISP_Ptr(iSG) * Elem_db_Vec(iE)%BF(iSL,iG)
                   F = F + F_Ptr(iSG) * Elem_db_Vec(iE)%BF(iSL,iG)
                End Do Do_iSL3
                Ener = Ener - Elem_db(iE)%Gauss_C(iG) * (F .DotP. U) 
              End If
          End Do

          Call Destroy_Gauss_EXO(Elem_db_Vec, Elem=iE)
          Call Destroy_Gauss_EXO(Elem_db_Scal, Elem=iE)
       EndDo Do_iE

    End Do Do_iBlk
    Call VecRestoreArrayF90(DISP_Loc, DISP_Ptr, iErr)
    Call VecRestoreArrayF90(F_Loc, F_Ptr, iErr)
    Call VecRestoreArrayF90(Temp_Loc, Temp_Ptr, iErr)

  End Subroutine Calc_Ener


#ifdef PB_2D
End Module m_Elast2D_Proc
#else
End Module m_Elast3D_Proc
#endif
