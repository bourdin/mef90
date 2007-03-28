Module m_Rupt2DA_Ener

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

  Public :: Comp_Bulk_Ener
  Public :: Comp_Surf_Ener

  Integer  :: iErr

Contains
  Subroutine Comp_Bulk_Ener(Ener, ULoc, VLoc, Geom, Params, SD_U, SD_V, Elems_U, Elems_V, Nodes_U, Nodes_V, t)
    Type (EXO_Geom_Info)                         :: Geom
    Type (Rupt_Params)                           :: Params
    Type (SD_Info)                               :: SD_U
    Type (SD_Info)                               :: SD_V
    Real(Kind = Kr),Intent(IN)                   :: t
    

    Type (Node2D), Dimension(:), Pointer          :: Nodes_U
    Type (Node2D), Dimension(:), Pointer          :: Nodes_V
    Type (Element2D_Scal), Dimension(:), Pointer  :: Elems_U
    Type (Element2D_Scal), Dimension(:), Pointer  :: Elems_V


    Vec                                          :: ULoc
    Vec                                          :: VLoc
    Real(Kind = Kr), Intent(OUT)                 :: Ener

    PetscReal                                    :: E, Nu
    PetscReal                                    :: K2
    PetscReal                                    :: Toughness    
    Type(Vect2D)                                 :: DistorsionEps, DistorsionSig

    Real(Kind = Kr)                              :: MyEner
    Integer                                      :: Nb_Gauss, Nb_DoF
    Integer                                      :: iSLV1, iSGV1
    Integer                                      :: iSLV2, iSGV2
    Integer                                      :: iSLSig, iSGSig
    Integer                                      :: iSLEps, iSGEps
    Integer                                      :: iBlk, iELoc, iE, iG
    

    Real(Kind = Kr), Dimension(:), Pointer       :: UPtr, VPtr
    Real(Kind = Kr), Dimension(:), Pointer       :: ContrV
    Real(Kind = Kr)                              :: ContrU

    Integer, Dimension(:), Pointer               :: Loc_Indices_U
    Integer, Dimension(:), Pointer               :: Loc_Indices_V


    Integer                                      :: i

    Allocate(Loc_Indices_V(Geom%Num_Nodes))
    Loc_Indices_V = (/ (i ,i = 0, Geom%Num_Nodes-1) /)
    Call AOApplicationToPETSc(SD_V%Loc_AO, Geom%Num_Nodes, Loc_Indices_V, iErr)


    Allocate(Loc_Indices_U(Geom%Num_Nodes))
    Loc_Indices_U = (/ (i ,i = 0, Geom%Num_Nodes-1) /)
    Call AOApplicationToPETSc(SD_U%Loc_AO, Geom%Num_Nodes, Loc_Indices_U, iErr)


    Call VecGetArrayF90(VLoc, VPtr, iErr)
    Call VecGetArrayF90(ULoc, UPtr, iErr)

    MyEner = 0.0_Kr
    Ener = 0.0_Kr
    Do_iBlk: Do iBlk = 1, Geom%Num_elem_blks
       E         = Params%Young_Mod(iBlk)
       nu        = Params%Poisson_Ratio(iBlk) 
       Toughness = Params%Toughness(iBlk)

!!! The isotropic Hooke's law is expressed as
!!! \sigma = K1 * trace(Epsilon) Id + 2*K2 * Epsilon - K3 Temp * Id
!!! K1, K2, K3 are computed in terms of E and nu
!!! in 3D, K1 = lambda, K2 = mu, K3 = E*alpha/(1-2nu) (= 3kappa alpha)
!!!       (alpha = therm exp coef).
!!! in 2D / plane stresses, the expressions are more complicated
!!!

       K2 = E / (1.0_Kr + nu) * InvOf2

       Do_iE: Do iELoc = 1, Geom%Elem_Blk(iBlk)%Num_Elems
          iE = Geom%Elem_Blk(iBlk)%Elem_ID(iELoc)
          If ((.NOT. SD_U%IsLocal_Elem(iE)) .OR. (.NOT. Params%Is_Domain(iBlk))) Then
             CYCLE
          End If

          Call Init_Gauss_EXO(Elems_U, Nodes_U, Geom, MEF90_GaussOrder, Elem=iE)
          Call Init_Gauss_EXO(Elems_V, Nodes_V, Geom, MEF90_GaussOrder, Elem=iE)

          Nb_Gauss = Elems_U(iE)%Nb_Gauss

          Allocate(ContrV(Nb_Gauss))
          
! V^2 + K_\epsilon term
          Is_Brittle: If (Params%Is_Brittle(iBlk)) Then
             ContrV = Params%KEpsilon
             Do_iSLV1: Do iSLV1 = 1, Elems_V(iE)%Nb_DoF
                iSGV1 = Elems_V(iE)%ID_DoF(iSLV1)
                DoiSLV2: Do iSLV2 = 1, Elems_V(iE)%Nb_DoF
                   iSGV2 = Elems_V(iE)%ID_DoF(iSLV2)
                   Do_iGV: Do iG = 1, Nb_Gauss
                      ContrV(iG) = ContrV(iG) + Elems_V(iE)%BF(iSLV1, iG) * Elems_V(iE)%BF(iSLV2, iG) *                            &
                           &            ( VPtr(Loc_Indices_V(iSGV1)+1) * VPtr(Loc_Indices_V(iSGV2)+1) ) 
                   End Do Do_iGV
                End Do DoiSLV2
             End Do Do_iSLV1
          Else
             ContrV = 1.0_Kr
          End If Is_Brittle

! W(e(u)) Term
!!! part related to v^2+k_\e W(e(u))
          ContrU = 0.0_Kr
          Do_iGUEps: Do iG = 1, Nb_Gauss
             Do_iSLSig: Do iSLSig = 1, Elems_U(iE)%Nb_DoF
                iSGSig = Elems_U(iE)%ID_DoF(iSLSig)
                DistorsionSig%X =   Nodes_U(iSGSig)%Coord%Y * Elems_U(iE)%BF(iSLSig,iG)
                DistorsionSig%Y = - Nodes_U(iSGSig)%Coord%X * Elems_U(iE)%BF(iSLSig,iG)
                
                Do_iSLEps: Do iSLEps = 1, Elems_U(iE)%Nb_DoF
                   iSGEps = Elems_U(iE)%ID_DoF(iSLEps)
                   DistorsionEps%X =   Nodes_U(iSGEps)%Coord%Y * Elems_U(iE)%BF(iSLEps,iG)
                   DistorsionEps%Y = - Nodes_U(iSGEps)%Coord%X * Elems_U(iE)%BF(iSLEps,iG)
                   
!                   ContrU = ContrU + ( Elems_U(iE)%Grad_BF(iSLEps,iG) .DotP. Elems_U(iE)%Grad_BF(iSLSig,iG) ) *            &
!                        &                          UPtr(Loc_Indices_U(iSGSig)+1) * UPtr(Loc_Indices_U(iSGEps)+1) * K2 * ContrV(iG)

                   MyEner = MyEner + ( (UPtr(Loc_Indices_U(iSGEps)+1) * Elems_U(iE)%Grad_BF(iSLEps,iG) - t*DistorsionEps) .DotP.      &
                                    (UPtr(Loc_Indices_U(iSGSig)+1) * Elems_U(iE)%Grad_BF(iSLSig,iG) - t*DistorsionSig) ) *            &
                                    Elems_U(iE)%Gauss_C(iG) * ContrV(iG) * .5_Kr
                End Do Do_iSLEps
             End Do Do_iSLSig
          End Do Do_iGUEps
          Call Destroy_Gauss_EXO(Elems_U, Elem=iE)
          Call Destroy_Gauss_EXO(Elems_V, Elem=iE)
          DeAllocate(ContrV)

       End Do Do_iE
    End Do Do_iBlk

    Call MPI_Reduce(MyEner, Ener, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, iErr)
    Call MPI_BCAST(Ener, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, iErr)
    Call VecRestoreArrayF90(ULoc, UPtr, iErr)
    Call VecRestoreArrayF90(VLoc, VPtr, iErr)
    DeAllocate(Loc_Indices_U)
    DeAllocate(Loc_Indices_V)
  End Subroutine Comp_Bulk_Ener

  Subroutine Comp_Surf_Ener(Ener, VLoc, Geom, Params, SD, Elems, Nodes)
    Type (EXO_Geom_Info)                         :: Geom
    Type (Rupt_Params)                           :: Params
    Type (SD_Info)                               :: SD


    Type (Node2D), Dimension(:), Pointer         :: Nodes
    Type (Element2D_Scal), Dimension(:), Pointer :: Elems

    Vec                                          :: VLoc
    Real(Kind = Kr), Intent(OUT)                 :: Ener


    Real(Kind = Kr)                              :: MyEner1, MyEner2
    Real(Kind = Kr)                              :: Ener1, Ener2

    Integer                                      :: Nb_Gauss, Nb_DoF
    Integer                                      :: iSL1, iSG1
    Integer                                      :: iSL2, iSG2
    Integer                                      :: iBlk, iELoc, iE, iG

    Real(Kind = Kr), Dimension(:), Pointer       :: VPtr
    PetscReal                                    :: Toughness

    Integer, Dimension(:), Pointer               :: Loc_Indices

    Integer                                      :: SetN, i

    Allocate(Loc_Indices(Geom%Num_Nodes))
    Loc_Indices = (/ (i ,i = 0, Geom%Num_Nodes - 1) /)
    Call AOApplicationToPETSc(SD%Loc_AO, Geom%Num_Nodes, Loc_Indices, iErr)

    Call VecGetArrayF90(VLoc, VPtr, iErr)
    
    MyEner1 = 0.0_Kr
    MyEner2 = 0.0_Kr
    Ener = 0.0_Kr

    Do_iBlk: Do iBlk = 1, Geom%Num_elem_blks
       Toughness = Params%Toughness(iBlk)

       Do_iE: Do iELoc = 1, Geom%Elem_Blk(iBlk)%Num_Elems
          iE = Geom%Elem_Blk(iBlk)%Elem_ID(iELoc)
          If (.NOT. SD%IsLocal_Elem(iE)) Then
             CYCLE
          End If

          Call Init_Gauss_EXO(Elems, Nodes, Geom, MEF90_GaussOrder, Elem=iE)
          Nb_Gauss = Elems(iE)%Nb_Gauss
          Do_iSL1: Do iSL1 = 1, Elems(iE)%Nb_DoF
             iSG1 = Elems(iE)%ID_DoF(iSL1)
             DoiSL2: Do iSL2 = 1, Elems(iE)%Nb_DoF
                iSG2 = Elems(iE)%ID_DoF(iSL2)
                Do_iG: Do iG = 1, Nb_Gauss
                   MyEner1 = MyEner1 + Toughness * Elems(iE)%Gauss_C(iG) * ( Elems(iE)%BF(iSL1, iG) * Elems(iE)%BF(iSL2, iG) *      &
                         (1.0_Kr - VPtr(Loc_Indices(iSG1)+1)) * (1.0_Kr - VPtr(Loc_Indices(iSG2)+1)) / Params%Epsilon * .25_Kr )

                   MyEner2 = MyEner2 + Toughness * Elems(iE)%Gauss_C(iG) * ( (Elems(iE)%Grad_BF(iSL1, iG) .DotP.                    &
                             Elems(iE)%Grad_BF(iSL2, iG)) * VPtr(Loc_Indices(iSG1)+1) * VPtr(Loc_Indices(iSG2)+1) * Params%Epsilon )
                End Do Do_iG
             End Do DoiSL2
          End Do Do_iSL1
          Call Destroy_Gauss_EXO(Elems, Elem=iE)
       End Do Do_iE
    End Do Do_iBlk
    
    Call MPI_Reduce(MyEner1, Ener1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, iErr)
    Call MPI_BCAST(Ener1, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, iErr)
    Call MPI_Reduce(MyEner2, Ener2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, iErr)
    Call MPI_BCAST(Ener2, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, iErr)
    Ener = Ener1 + Ener2
    Call VecRestoreArrayF90(VLoc, VPtr, iErr)
    DeAllocate(Loc_Indices)
  End Subroutine Comp_Surf_Ener


End Module m_Rupt2DA_Ener

