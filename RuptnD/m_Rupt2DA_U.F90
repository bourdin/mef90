#if defined PB_2D
Module m_Rupt2D_U
#elif defined PB_3D
Module m_Rupt3D_U
#else
Module m_Rupt2DA_U
#endif
  Use m_MEF90
  Use m_Rupt_Struct

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

  Public :: Assemb_MR_U
  Public :: Assemb_RHS_U

  Integer :: iErr
Contains

  Subroutine Assemb_MR_U(MR, V_Loc, Geom, Params, SD_U, SD_V,                 &
       & Elems_U, Elems_V, Nodes_U, Nodes_V)
    Mat                                           :: MR
    Vec                                           :: V_Loc
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params)                            :: Params
    Type (SD_Info)                                :: SD_U
    Type (SD_Info)                                :: SD_V


#if defined PB_2D
    Type (Node2D), Dimension(:), Pointer          :: Nodes_U
    Type (Node2D), Dimension(:), Pointer          :: Nodes_V
    Type (Element2D_Elast), Dimension(:), Pointer :: Elems_U
    Type (Element2D_Scal), Dimension(:), Pointer  :: Elems_V
    Type(MatS2D), Dimension(:), Pointer           :: Sigma
#elif defined PB_3D
    Type (Node3D), Dimension(:), Pointer          :: Nodes_U
    Type (Node3D), Dimension(:), Pointer          :: Nodes_V
    Type (Element3D_Elast), Dimension(:), Pointer :: Elems_U
    Type (Element3D_Scal), Dimension(:), Pointer  :: Elems_V
    Type(MatS3D), Dimension(:), Pointer           :: Sigma
#else
    Type (Node2D), Dimension(:), Pointer          :: Nodes_U
    Type (Node2D), Dimension(:), Pointer          :: Nodes_V
    Type (Element2D_Scal), Dimension(:), Pointer  :: Elems_U
    Type (Element2D_Scal), Dimension(:), Pointer  :: Elems_V
#endif

    Integer                                       :: Nb_Gauss, Nb_DoF_U
    Integer                                       :: iSLSig, iSGSig
    Integer                                       :: iSLEps
    Integer                                       :: iSLV1, iSGV1
    Integer                                       :: iSLV2, iSGV2
    Integer                                       :: iBlk, iEloc, iE, iG


    Real(Kind = Kr), Dimension(:), Pointer        :: V_Ptr
    Real(Kind = Kr), Dimension(:), Pointer        :: ContrV

    PetscScalar, Dimension(:,:), Pointer          :: MR_Elem
    Integer, Dimension(:), Pointer                :: EXO_Indices_U
    Integer, Dimension(:), Pointer                :: Loc_Indices_V

    PetscTruth                                    :: ISAssembled

    PetscLogDouble                                :: GaussTS, GaussTF, GaussT
    PetscLogDouble                                :: SetTS, SetTF, SetT

    Integer                                       :: SetN, i

    Real(Kind = Kr)                               :: E, Nu, Lambda, Mu

    SetN = 0
    GaussT = 0.0

    Call MatAssembled(MR, IsAssembled, iErr)
    If (IsAssembled) Then
       Call MatZeroEntries(MR, iErr)
    End If

    
    Allocate(Loc_Indices_V(Geom%Num_Nodes))
    Loc_Indices_V = (/ (i ,i = 0, Geom%Num_Nodes-1) /)
    Call AOApplicationToPETSc(SD_V%Loc_AO, Geom%Num_Nodes, Loc_Indices_V, iErr)

    Call VecGetArrayF90(V_Loc, V_Ptr, iErr) 

    Do_iBlk: Do iBlk = 1, Geom%Num_elem_blks
       E  = Params%Young_Mod(iBlk)
       nu = Params%Poisson_Ratio(iBlk) 

#if defined PB_2D
       Lambda = E * nu / (1.0_Kr - nu**2)
       Mu     = E / (1.0_Kr + nu) * InvOf2
#elif defined PB_3D
       Lambda = E * nu / (1.0_Kr - 2.0_Kr * nu) / ( 1.0_Kr + nu)
       Mu     = E / (1.0_Kr + nu) * InvOf2
#else
       Mu     = E / (1.0_Kr + nu) * InvOf2
#endif

#ifdef PB_2DA
       Nb_DoF_U = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem
#else
       Nb_DoF_U = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem * Geom%Num_Dim
#endif

       Allocate (MR_Elem(Nb_DoF_U, Nb_DoF_U))
       Allocate (EXO_Indices_U(Nb_DoF_U))

       Do_iE: Do iELoc = 1, Geom%Elem_Blk(iBlk)%Num_Elems
          iE = Geom%Elem_Blk(iBlk)%ELem_ID(iELoc)
          If ((.NOT. SD_U%IsLocal_Elem(iE)) .OR.                              &
               & (.NOT. Params%Is_Domain(iBlk))) Then
             CYCLE
          End If
          EXO_Indices_U = Elems_U(iE)%Id_DoF-1
          Call AOApplicationToPETSc(SD_U%EXO_AO, Nb_DoF_U, EXO_Indices_U, iErr)
       
          Call PetscGetTime(GaussTS, iErr)
          Call Init_Gauss_EXO(Elems_U, Nodes_U, Geom, MEF90_GaussOrder, Elem=iE)
          Call Init_Gauss_EXO(Elems_V, Nodes_V, Geom, MEF90_GaussOrder, Elem=iE)
          Call PetscGetTime(GaussTF, iErr)
          GaussT = Gausst + GaussTF - GaussTS
          MR_Elem = 0.0_Kr
          Nb_Gauss = Elems_U(iE)%Nb_Gauss
          Allocate(ContrV(Nb_Gauss))
#ifndef PB_2DA
          Allocate(Sigma(Nb_Gauss))
#endif
          Is_Brittle: If ( Params%Is_Brittle(iBlk)) Then
!             ContrV = Params%Kepsilon
             ContrV = 0.0_Kr
             Do_iSLV1: Do iSLV1 = 1, Elems_V(iE)%Nb_DoF
                iSGV1 = Elems_V(iE)%ID_DoF(iSLV1)
                DoiSLV2: Do iSLV2 = 1, Elems_V(iE)%Nb_DoF
                   iSGV2 = Elems_V(iE)%ID_DoF(iSLV2)
                   Do_iGV: Do iG = 1, Nb_Gauss
                      ContrV(iG) = ContrV(iG) +                               &
                           & V_Ptr(Loc_Indices_V(iSGV1)+1)                      &
                           &  * Elems_V(iE)%BF(iSLV1, iG)                     &
                           &  * V_Ptr(Loc_Indices_V(iSGV2)+1)                   &
                           &  * Elems_V(iE)%BF(iSLV2, iG)
                   End Do Do_iGV
                End Do DoiSLV2
              End Do Do_iSLV1
          Else
             ContrV = 1.0_Kr
          End If Is_Brittle

          Do_iSLSig: Do iSLSig = 1, Elems_U(iE)%Nb_DoF
#ifndef PB_2DA
             Do_iGSig: Do iG = 1, Nb_Gauss
                Sigma(iG) = 2.0_Kr * Mu * Elems_U(iE)%GradS_BF(iSLSig,iG)
                Sigma(iG)%XX = Sigma(iG)%XX +&
                     & Lambda * Trace(Elems_U(iE)%GradS_BF(iSLSig,iG))
                Sigma(iG)%YY = Sigma(iG)%YY +&
                     & Lambda * Trace(Elems_U(iE)%GradS_BF(iSLSig,iG))
#ifdef PB_3D
                Sigma(iG)%ZZ = Sigma(iG)%ZZ +&
                     & Lambda * Trace(Elems_U(iE)%GradS_BF(iSLSig,iG))
#endif
             End Do Do_iGSig
#endif
             Do_iSLEps: Do iSLEps = 1, Elems_U(iE)%Nb_DoF
                Do_iGUEps: Do iG = 1, Nb_Gauss
#if defined PB_2DA
                   MR_Elem(iSLEps, iSLSig) = MR_Elem(iSLEps, iSLSig) +        &
                        & Elems_U(iE)%Gauss_C(iG) *                           &
                        & ( ContrV(iG) + Params%Kepsilon ) *                  &
                        & ( Elems_U(iE)%Grad_BF(iSLEps,iG) .DotP.             &
                        &   Elems_U(iE)%Grad_BF(iSLSig,iG) ) * Mu
#else                   
                   MR_Elem(iSLEps, iSLSig) = MR_Elem(iSLEps, iSLSig) +        &
                        & Elems_U(iE)%Gauss_C(iG) * ( ContrV(iG) *            &
                        & (Elems_U(iE)%GradS_BF(iSLEps,iG) .DotP. Sigma(iG))  &
                        &  + Params%Kepsilon *                                &
                        &  (Elems_U(iE)%GradS_BF(iSLEps,iG) .DotP.            &
                        &   Elems_U(iE)%GradS_BF(iSLSig,iG)))
#endif
                End Do Do_iGUEps
             End Do Do_iSLEps
          End Do Do_iSLSig

          Call MatSetValues(MR, Nb_DoF_U, EXO_Indices_U, Nb_DoF_U,            &
               & EXO_Indices_U, MR_Elem, ADD_VALUES, iErr)
          SetN = SetN + 1

          Call PetscGetTime(GaussTS, iErr)
          Call Destroy_Gauss_EXO(Elems_U, Elem=iE)
          Call Destroy_Gauss_EXO(Elems_V, Elem=iE)
          Call PetscGetTime(GaussTF, iErr)
          GaussT = Gausst + GaussTF - GaussTS
          DeAllocate(ContrV)
#ifndef PB_2DA
          DeAllocate(Sigma)
#endif
       EndDo Do_iE
       DeAllocate(MR_Elem)
       DeAllocate(EXO_Indices_U)
    End Do Do_iBlk

    Call MatAssemblyBegin(MR, MAT_FLUSH_ASSEMBLY, iErr)
    DeAllocate(Loc_Indices_V)

    Call VecRestoreArrayF90(V_Loc, V_Ptr, iErr) 
    Call MatAssemblyEnd(MR, MAT_FLUSH_ASSEMBLY, iErr)  


    ! Assembly of the BC terms 
#ifdef PB_2DA
    Allocate(EXO_Indices_U(Geom%Num_Nodes))
    EXO_Indices_U = (/ (i ,i = 0, Geom%Num_Nodes - 1) /)
    Call AOApplicationToPETSc(SD_U%EXO_AO, Geom%Num_Nodes, EXO_Indices_U, iErr)

#else
    Allocate(EXO_Indices_U(Geom%Num_Nodes * Geom%Num_Dim))
    EXO_Indices_U = (/ (i ,i = 0, Geom%Num_Nodes * Geom%Num_Dim - 1) /)
    Call AOApplicationToPETSc(SD_U%EXO_AO, Geom%Num_Nodes * Geom%Num_Dim,    &
         & EXO_Indices_U, iErr)
#endif 

    Do_iBlk_BC: Do iBlk = 1, Geom%Num_Elem_Blks
#ifdef PB_2DA
       Nb_DoF_U = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem
#else
       Nb_DoF_U = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem * Geom%Num_Dim
#endif

       Do_iE_BC: Do iELoc = 1, Geom%Elem_blk(iBlk)%Num_Elems
          iE = Geom%Elem_blk(iBlk)%Elem_ID(iELoc)
          If (.NOT. SD_U%IsLocal_Elem(iE)) Then
             CYCLE
          End If
          
          Do_iSLSig_BC: Do iSLSig = 1, Nb_DoF_U
             ISGSig = Elems_U(iE)%ID_DoF(iSLSig)
             Is_BC_BC: If ( (Nodes_U(iSGSig)%BC /= BC_Type_NONE) .AND.        &
                  & (SD_U%IsLocal_Node(iSGSig)) )Then
                Call PetscGetTime(SetTS, iErr)
                Call MatSetValue(MR, EXO_Indices_U(iSGSig), EXO_Indices_U(iSGSig),&
                     &  MEF90_VLV, INSERT_VALUES, iErr)
                Call PetscGetTime(SetTF, iErr)
                SetN = SetN+1
                SetT = SetT + SetTF - SetTS
             End If Is_BC_BC
          End Do Do_iSLSig_BC
       End Do Do_iE_BC
    End Do Do_iBlk_BC
    DeAllocate (EXO_Indices_U)

    Call MatAssemblyBegin(MR, MAT_FINAL_ASSEMBLY, iErr)
    Call MatAssemblyEnd(MR, MAT_FINAL_ASSEMBLY, iErr)  

#ifdef MEF90_TIMING
    If (MEF90_MyRank ==0) Then
       Write(CharBuffer,*) 'Total time in Init/Destroy_Gauss:         ',      &
            & GaussT, '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
       Write(CharBuffer,*) 'Number of calls to MatSetValue:           ',      &
            & SetN, '\n'c
       Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
    End If
#endif
  End Subroutine Assemb_MR_U


  Subroutine Assemb_RHS_U(RHS, BC_U_Loc, Geom, Params, SD_U, Elems_U, Nodes_U,&
       & SD_V, Elems_V, Nodes_V, Vloc, FLoc, TempLoc)
    Vec                                           :: RHS
    Vec                                           :: BC_U_Loc
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params)                            :: Params
    Type (SD_Info)                                :: SD_U
    Type (SD_Info)                                :: SD_V

#if defined PB_2D
    Type (Node2D), Dimension(:), Pointer          :: Nodes_U
    Type (Element2D_Elast), Dimension(:), Pointer :: Elems_U
    Type (Node2D), Dimension(:), Pointer          :: Nodes_V
    Type (Element2D_Scal), Dimension(:), Pointer  :: Elems_V
#elif defined PB_3D
    Type (Node3D), Dimension(:), Pointer          :: Nodes_U
    Type (Element3D_Elast), Dimension(:), Pointer :: Elems_U
    Type (Node3D), Dimension(:), Pointer          :: Nodes_V
    Type (Element3D_Scal), Dimension(:), Pointer  :: Elems_V
#else
    Type (Node2D), Dimension(:), Pointer          :: Nodes_U
    Type (Element2D_Scal), Dimension(:), Pointer  :: Elems_U
    Type (Node2D), Dimension(:), Pointer          :: Nodes_V
    Type (Element2D_Scal), Dimension(:), Pointer  :: Elems_V
#endif

    Vec                                           :: VLoc
    Vec                                           :: FLoc
    Vec                                           :: TempLoc

    Real(Kind = Kr), Dimension(:), Pointer        :: F_Ptr
    Real(Kind = Kr), Dimension(:), Pointer        :: Temp_Ptr
    Real(Kind = Kr), Dimension(:), Pointer        :: V_Ptr

    Integer, Dimension(:), Pointer                :: Loc_Indices_Vect
    Integer, Dimension(:), Pointer                :: EXO_Indices_Vect
    Integer, Dimension(:), Pointer                :: Loc_Indices_Scal
    PetscReal, Dimension(:), Pointer              :: RHS_Ptr
    PetscReal, Dimension(:), Pointer              :: BC_U_Ptr
    Integer                                       :: i, iBlk
    Integer                                       :: iE, iELoc
    Integer                                       :: iSL1, iSG1
    Integer                                       :: iSL2, iSG2
    Integer                                       :: iSLV1, iSGV1
    Integer                                       :: iSLV2, iSGV2
    Integer                                       :: iS, iSLoc
    Integer                                       :: iG
    Real(Kind = Kr)                               :: Tmp_Val
    PetscLogDouble                                :: TotTS, TotTF, TotT
    Real(Kind = Kr), Dimension(:), Pointer        :: ContrV

    Call PetscGetTime(TotTS, iErr)
    
    Allocate(Loc_Indices_Scal(Geom%Num_Nodes))
    Loc_Indices_Scal = (/ (i ,i = 0, Geom%Num_Nodes - 1) /)
    Call AOApplicationToPETSc(SD_V%Loc_AO, Geom%Num_Nodes, Loc_Indices_Scal,  &
         & iErr)
    
#ifdef PB_2DA
    Allocate(Loc_Indices_Vect(Geom%Num_Nodes))
    Loc_Indices_Vect = (/ (i ,i = 0, Geom%Num_Nodes - 1) /)
    Call AOApplicationToPETSc(SD_U%Loc_AO, Geom%Num_Nodes, Loc_Indices_Vect,  &
         & iErr)
    
    Allocate(EXO_Indices_Vect(Geom%Num_Nodes))
    EXO_Indices_Vect = (/ (i ,i = 0, Geom%Num_Nodes - 1) /)
    Call AOApplicationToPETSc(SD_U%EXO_AO, Geom%Num_Nodes, EXO_Indices_Vect,  &
         & iErr)
#else
    Allocate(Loc_Indices_Vect(Geom%Num_Nodes * Geom%Num_Dim))
    Loc_Indices_Vect = (/ (i ,i = 0, Geom%Num_Nodes * Geom%Num_Dim -1) /)
    Call AOApplicationToPETSc(SD_U%Loc_AO, Geom%Num_Nodes * Geom%Num_Dim,    &
         & Loc_Indices_Vect, iErr)
    
    Allocate(EXO_Indices_Vect(Geom%Num_Nodes * Geom%Num_Dim))
    EXO_Indices_Vect = (/ (i ,i = 0, Geom%Num_Nodes * Geom%Num_Dim -1) /)
    Call AOApplicationToPETSc(SD_U%EXO_AO, Geom%Num_Nodes * Geom%Num_Dim,    &
         & EXO_Indices_Vect, iErr)
#endif
    
    Call VecGetArrayF90(FLoc, F_Ptr, iErr)
    Call VecGetArrayF90(TempLoc, Temp_Ptr, iErr)
    Call VecGetArrayF90(VLoc, V_Ptr, iErr)
    
    Call VecSet(0.0_Kr, RHS, iErr)
    Do_iBlk: Do iBlk = 1, Geom%Num_elem_blks
       !        If (.NOT. Params%Has_Force(iBlk)) Then
       !           CYCLE
       !        End If
!!! Let's see if this has a performance impact.
       
       Do_iE: Do iELoc = 1, Geom%Elem_blk(iBlk)%Num_Elems
          iE = Geom%Elem_blk(iBlk)%Elem_ID(iELoc)
          If (.NOT. SD_U%IsLocal_Elem(iE)) Then
             CYCLE
          End If
          Call Init_Gauss_EXO(Elems_U, Nodes_U, Geom, MEF90_GaussOrder,       &
               & Elem=iE)
          Call Init_Gauss_EXO(Elems_V, Nodes_V, Geom, MEF90_GaussOrder,       &
               & Elem=iE)
          
          Allocate(ContrV(Elems_V(iE)%Nb_Gauss))
          
#ifndef PB_2DA
          Is_Brittle: If ( Params%Is_Brittle(iBlk)) Then
!!! the antiplane hypothesis doesn;t make sense for thermoelasticity
             ContrV = Params%KEpsilon
             Do_iSLV1: Do iSLV1 = 1, Elems_V(iE)%Nb_DoF
                iSGV1 = Elems_V(iE)%ID_DoF(iSLV1)
                DoiSLV2: Do iSLV2 = 1, Elems_V(iE)%Nb_DoF
                   iSGV2 = Elems_V(iE)%ID_DoF(iSLV2)
                   Do_iGV: Do iG = 1, Elems_V(iE)%Nb_Gauss
                      ContrV(iG) = ContrV(iG) +                               &
                           & V_Ptr(Loc_Indices_Scal(iSGV1)+1)                 &
                           &  * Elems_V(iE)%BF(iSLV1, iG)                     &
                           &  * V_Ptr(Loc_Indices_Scal(iSGV2)+1)              &
                           &  * Elems_V(iE)%BF(iSLV2, iG)
                   End Do Do_iGV
                End Do DoiSLV2
             End Do Do_iSLV1
          Else
             ContrV = 1.0_Kr
          End If Is_Brittle
#endif
          
          Do_iSL1: Do iSL1 = 1, Elems_U(iE)%Nb_DoF
             iSG1 = Elems_U(iE)%ID_DoF(iSL1)
             Tmp_Val = 0.0_Kr
             
!!! Force part
             Has_Force: If (Params%Has_Force(iBlk)) Then
                Do_iSL2: Do iSL2 = 1, Elems_U(iE)%Nb_Gauss
                   iSG2 = Elems_U(iE)%ID_DoF(iSL2)
                   Do_iG: Do iG = 1, Elems_U(iE)%Nb_Gauss
#ifdef PB_2DA
                      Tmp_Val = Tmp_Val + Elems_U(iE)%Gauss_C(iG) * (         &
                           & F_Ptr(Loc_Indices_Vect(iSG2)+1) *                &
                           & ( Elems_U(iE)%BF(iSL1, iG) *                     &
                           &   Elems_U(iE)%BF(iSL2,iG) ))
#else
                      Tmp_Val = Tmp_Val + Elems_U(iE)%Gauss_C(iG) * (         &
                           & F_Ptr(Loc_Indices_Vect(iSG2)+1) *                &
                           & ( Elems_U(iE)%BF(iSL1, iG) .DotP.                &
                           &   Elems_U(iE)%BF(iSL2,iG) ))
#endif
                   End Do Do_iG
                End Do Do_iSL2
             EndIf Has_Force
             
!!! Thermal expansion part
#ifndef PB_2DA
             Do_iSL2_Temp: Do iSL2 = 1, Elems_V(iE)%Nb_Gauss
                iSG2 = Elems_V(iE)%ID_DoF(iSL2)
                Do_iG_Temp: Do iG = 1, Elems_U(iE)%Nb_Gauss
                   Tmp_Val = Tmp_Val + Elems_V(iE)%Gauss_C(iG) * (            &
                        & Temp_Ptr(Loc_Indices_Scal(iSG2)+1) *                &
                        & Params%Therm_Exp(iBlk) *                            &
                        & trace(Elems_U(iE)%GradS_BF(iSL1, iG)) *             &
                        & Elems_V(iE)% BF(iSL2, iG)) * ContrV(iG)
                End Do Do_iG_Temp
             End Do Do_iSL2_Temp
#endif
             
             Call VecSetValue(RHS, EXO_Indices_Vect(iSG1), Tmp_Val,           &
                  & ADD_VALUES, iErr)
          End Do Do_iSL1
          
          Call Destroy_Gauss_EXO(Elems_U, Elem=iE)
          Call Destroy_Gauss_EXO(Elems_V, Elem=iE)
          DeAllocate(ContrV)
       End Do Do_iE
    End Do Do_iBlk
    Call VecAssemblyBegin(RHS, iErr)
    DeAllocate(EXO_Indices_Vect)
    Call VecAssemblyEnd(RHS, iErr)
    
!!!
!!! Boundary Conditions, using VLV
!!!
    
    Call VecGetArrayF90(BC_U_Loc, BC_U_Ptr, iErr)
    Call VecGetArrayF90(RHS, RHS_Ptr, iErr)
    
    Do_iS: Do iSLoc = 1, SD_U%Num_Nodes
       iS = SD_U%Node(iSLoc)
       Is_BC: If (Nodes_U(iS)%BC /= BC_Type_NONE) Then
          RHS_Ptr(Loc_Indices_Vect(iS)+1) = BC_U_Ptr(Loc_Indices_Vect(iS)+1)  &
               & * MEF90_VLV
       End If Is_BC
    End Do Do_iS
    
    
    Call PetscGetTime(TotTF, iErr)
    TotT = TotTF - TotTS

    DeAllocate(Loc_Indices_Vect)
    DeAllocate(Loc_Indices_Scal)
    
    Call VecRestoreArrayF90(VLoc, V_Ptr, iErr)
    Call VecRestoreArrayF90(TempLoc, Temp_Ptr, iErr)
    Call VecRestoreArrayF90(FLoc, F_Ptr, iErr)
    Call VecRestoreArrayF90(BC_U_Loc, BC_U_Ptr, iErr)
    Call VecRestoreArrayF90(RHS, RHS_Ptr, iErr)
  End Subroutine Assemb_RHS_U
  


#if defined PB_2D
End Module m_Rupt2D_U
#elif defined PB_3D
End Module m_Rupt3D_U
#else
End Module m_Rupt2DA_U
#endif
