Module m_Rupt2DA_U

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

  Subroutine Assemb_MR_U(MR, V_Loc, Geom, Params, SD_U, SD_V, Elems_U, Elems_V, Nodes_U, Nodes_V)
    Mat                                           :: MR
    Vec                                           :: V_Loc
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params)                            :: Params
    Type (SD_Info)                                :: SD_U
    Type (SD_Info)                                :: SD_V



    Type (Node2D), Dimension(:), Pointer          :: Nodes_U
    Type (Node2D), Dimension(:), Pointer          :: Nodes_V
    Type (Element2D_Scal), Dimension(:), Pointer  :: Elems_U
    Type (Element2D_Scal), Dimension(:), Pointer  :: Elems_V


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

    Real(Kind = Kr)                               :: E, Nu
    Real(Kind = Kr)                               :: K2
    
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
       K2 = E / (1.0_Kr + nu) * InvOf2

!!! The isotropic Hooke's law is expressed as
!!! \sigma = K1 * trace(Epsilon) Id + 2*K2 * Epsilon - K3 Temp * Id
!!! K1, K2, K3 are computed in terms of E and nu
!!! in 3D, K1 = lambda, K2 = mu, K3 = E*alpha/(1-2nu) (= 3kappa alpha)
!!!       (alpha = therm exp coef).
!!! in 2D / plane stresses, the expressions are more complicated
!!!

       Nb_DoF_U = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem
       Allocate (MR_Elem(Nb_DoF_U, Nb_DoF_U))
       Allocate (EXO_Indices_U(Nb_DoF_U))

       Do_iE: Do iELoc = 1, Geom%Elem_Blk(iBlk)%Num_Elems
          iE = Geom%Elem_Blk(iBlk)%ELem_ID(iELoc)
          If ((.NOT. SD_U%IsLocal_Elem(iE)) .OR. (.NOT. Params%Is_Domain(iBlk))) Then
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

          Is_Brittle: If ( Params%Is_Brittle(iBlk)) Then
             ContrV = Params%Kepsilon
             Do_iSLV1: Do iSLV1 = 1, Elems_V(iE)%Nb_DoF
                iSGV1 = Elems_V(iE)%ID_DoF(iSLV1)
                DoiSLV2: Do iSLV2 = 1, Elems_V(iE)%Nb_DoF
                   iSGV2 = Elems_V(iE)%ID_DoF(iSLV2)
                   Do_iGV: Do iG = 1, Nb_Gauss
                      ContrV(iG) = ContrV(iG) +  V_Ptr(Loc_Indices_V(iSGV1)+1) * Elems_V(iE)%BF(iSLV1, iG) *                       &
                                                 V_Ptr(Loc_Indices_V(iSGV2)+1) * Elems_V(iE)%BF(iSLV2, iG)
                   End Do Do_iGV
                End Do DoiSLV2
              End Do Do_iSLV1
          Else
             ContrV = 1.0_Kr
          End If Is_Brittle

          Do_iSLSig: Do iSLSig = 1, Elems_U(iE)%Nb_DoF
             Do_iSLEps: Do iSLEps = 1, Elems_U(iE)%Nb_DoF
                Do_iGUEps: Do iG = 1, Nb_Gauss
                   MR_Elem(iSLEps, iSLSig) = MR_Elem(iSLEps, iSLSig) +  Elems_U(iE)%Gauss_C(iG) *                                  &
                                           ( ContrV(iG) + Params%Kepsilon ) *                                                      &
                                           ( Elems_U(iE)%Grad_BF(iSLEps,iG) .DotP. Elems_U(iE)%Grad_BF(iSLSig,iG) ) * K2

                End Do Do_iGUEps
             End Do Do_iSLEps
          End Do Do_iSLSig

          Call MatSetValues(MR, Nb_DoF_U, EXO_Indices_U, Nb_DoF_U, EXO_Indices_U, MR_Elem, ADD_VALUES, iErr)
          SetN = SetN + 1

          Call PetscGetTime(GaussTS, iErr)
          Call Destroy_Gauss_EXO(Elems_U, Elem=iE)
          Call Destroy_Gauss_EXO(Elems_V, Elem=iE)
          Call PetscGetTime(GaussTF, iErr)
          GaussT = Gausst + GaussTF - GaussTS
          DeAllocate(ContrV)
       EndDo Do_iE
       DeAllocate(MR_Elem)
       DeAllocate(EXO_Indices_U)
    End Do Do_iBlk

    Call MatAssemblyBegin(MR, MAT_FLUSH_ASSEMBLY, iErr)
    DeAllocate(Loc_Indices_V)

    Call VecRestoreArrayF90(V_Loc, V_Ptr, iErr) 
    Call MatAssemblyEnd(MR, MAT_FLUSH_ASSEMBLY, iErr)  


! Assembly of the BC terms

    Allocate(EXO_Indices_U(Geom%Num_Nodes))
    EXO_Indices_U = (/ (i ,i = 0, Geom%Num_Nodes - 1) /)
    Call AOApplicationToPETSc(SD_U%EXO_AO, Geom%Num_Nodes, EXO_Indices_U, iErr)

    Do_iBlk_BC: Do iBlk = 1, Geom%Num_Elem_Blks

       Nb_DoF_U = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem


       Do_iE_BC: Do iELoc = 1, Geom%Elem_blk(iBlk)%Num_Elems
          iE = Geom%Elem_blk(iBlk)%Elem_ID(iELoc)
          If (.NOT. SD_U%IsLocal_Elem(iE)) Then
             CYCLE
          End If
          
          Do_iSLSig_BC: Do iSLSig = 1, Nb_DoF_U
             ISGSig = Elems_U(iE)%ID_DoF(iSLSig)
             Is_BC_BC: If ( (Nodes_U(iSGSig)%BC /= BC_Type_NONE) .AND. (SD_U%IsLocal_Node(iSGSig)) )Then
                Call PetscGetTime(SetTS, iErr)
                Call MatSetValue(MR, EXO_Indices_U(iSGSig), EXO_Indices_U(iSGSig), MEF90_VLV, INSERT_VALUES, iErr)
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
  End Subroutine Assemb_MR_U


  Subroutine Assemb_RHS_U(RHS, BC_U_Loc, Geom, Params, SD_U, Elems_U, Nodes_U, SD_V, Elems_V, Nodes_V, Vloc, t)
    Vec                                           :: RHS
    Vec                                           :: BC_U_Loc
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params)                            :: Params
    Type (SD_Info)                                :: SD_U
    Type (SD_Info)                                :: SD_V
    Type (Node2D), Dimension(:), Pointer          :: Nodes_U
    Type (Element2D_Scal), Dimension(:), Pointer  :: Elems_U
    Type (Node2D), Dimension(:), Pointer          :: Nodes_V
    Type (Element2D_Scal), Dimension(:), Pointer  :: Elems_V
    Vec                                           :: VLoc
    Real(Kind = Kr), Intent(IN)                   :: t

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
    Integer                                       :: Nb_DoF_U
    Real(Kind = Kr)                               :: Tmp_Val, K2, E, nu
    PetscLogDouble                                :: TotTS, TotTF, TotT
    Real(Kind = Kr), Dimension(:), Pointer        :: ContrV
    Type(Vect2D)                                  :: Distorsion

    Call PetscGetTime(TotTS, iErr)
    
    Allocate(Loc_Indices_Scal(Geom%Num_Nodes))
    Loc_Indices_Scal = (/ (i ,i = 0, Geom%Num_Nodes - 1) /)
    Call AOApplicationToPETSc(SD_V%Loc_AO, Geom%Num_Nodes, Loc_Indices_Scal, iErr)

    Allocate(Loc_Indices_Vect(Geom%Num_Nodes))
    Loc_Indices_Vect = (/ (i ,i = 0, Geom%Num_Nodes - 1) /)
    Call AOApplicationToPETSc(SD_U%Loc_AO, Geom%Num_Nodes, Loc_Indices_Vect, iErr)
    
    Allocate(EXO_Indices_Vect(Geom%Num_Nodes))
    EXO_Indices_Vect = (/ (i ,i = 0, Geom%Num_Nodes - 1) /)
    Call AOApplicationToPETSc(SD_U%EXO_AO, Geom%Num_Nodes, EXO_Indices_Vect, iErr)
    
    Call VecGetArrayF90(VLoc, V_Ptr, iErr)
    
    Call VecSet(RHS, 0.0_Kr, iErr)
    Do_iBlk: Do iBlk = 1, Geom%Num_elem_blks
       E  = Params%Young_Mod(iBlk)
       nu = Params%Poisson_Ratio(iBlk) 
       K2 = E / (1.0_Kr + nu) * InvOf2

       Nb_DoF_U = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem

       Do_iE: Do iELoc = 1, Geom%Elem_blk(iBlk)%Num_Elems
          iE = Geom%Elem_blk(iBlk)%Elem_ID(iELoc)
          If (.NOT. SD_U%IsLocal_Elem(iE)) Then
             CYCLE
          End If
          Call Init_Gauss_EXO(Elems_U, Nodes_U, Geom, MEF90_GaussOrder, Elem=iE)
          Call Init_Gauss_EXO(Elems_V, Nodes_V, Geom, MEF90_GaussOrder, Elem=iE)
          
          Allocate(ContrV(Elems_V(iE)%Nb_Gauss))
          
          Is_Brittle: If ( Params%Is_Brittle(iBlk)) Then
             ContrV = Params%KEpsilon
             Do_iSLV1: Do iSLV1 = 1, Elems_V(iE)%Nb_DoF
                iSGV1 = Elems_V(iE)%ID_DoF(iSLV1)
                DoiSLV2: Do iSLV2 = 1, Elems_V(iE)%Nb_DoF
                   iSGV2 = Elems_V(iE)%ID_DoF(iSLV2)
                   Do_iGV: Do iG = 1, Elems_V(iE)%Nb_Gauss
                      ContrV(iG) = ContrV(iG) + V_Ptr(Loc_Indices_Scal(iSGV1)+1) * Elems_V(iE)%BF(iSLV1, iG) *                     &
                                                V_Ptr(Loc_Indices_Scal(iSGV2)+1) * Elems_V(iE)%BF(iSLV2, iG)
                   End Do Do_iGV
                End Do DoiSLV2
             End Do Do_iSLV1
          Else
             ContrV = 1.0_Kr
          End If Is_Brittle

             Do_iSL1: Do iSL1 = 1, Elems_U(iE)%Nb_DoF
                iSG1 = Elems_U(iE)%ID_DoF(iSL1)
                Tmp_Val = 0.0_Kr
                Do_iSL2_Temp: Do iSL2 = 1, Elems_V(iE)%Nb_DoF
                   iSG2 = Elems_V(iE)%ID_DoF(iSL2)
          Do_iG_Temp: Do iG = 1, Elems_U(iE)%Nb_Gauss
                   Distorsion%X =   Nodes_U(iSG2)%Coord%Y * Elems_U(iE)%BF(iSL2,iG)
	   			    Distorsion%Y = - Nodes_U(iSG2)%Coord%X * Elems_U(iE)%BF(iSL2,iG)
                   Tmp_Val = Tmp_Val + Elems_U(iE)%Gauss_C(iG) * ContrV(iG) * K2 * t *                                             &
                            ( Distorsion .DotP. Elems_U(iE)%Grad_BF(iSL1, iG) )
          End Do Do_iG_Temp
                End Do Do_iSL2_Temp             
                Call VecSetValue(RHS, EXO_Indices_Vect(iSG1), Tmp_Val, ADD_VALUES, iErr)
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
          RHS_Ptr(Loc_Indices_Vect(iS)+1) = BC_U_Ptr(Loc_Indices_Vect(iS)+1) * MEF90_VLV
       End If Is_BC
    End Do Do_iS
    
    
    Call PetscGetTime(TotTF, iErr)
    TotT = TotTF - TotTS

    DeAllocate(Loc_Indices_Vect)
    DeAllocate(Loc_Indices_Scal)
    
    Call VecRestoreArrayF90(VLoc, V_Ptr, iErr)
    Call VecRestoreArrayF90(BC_U_Loc, BC_U_Ptr, iErr)
    Call VecRestoreArrayF90(RHS, RHS_Ptr, iErr)
  End Subroutine Assemb_RHS_U
End Module m_Rupt2DA_U