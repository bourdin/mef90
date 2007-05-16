#if defined PB_2D
Module m_Rupt2D_Ener
#elif defined PB_3D
Module m_Rupt3D_Ener
#else
Module m_Rupt2DA_Ener
#endif
  Use m_MEF90
  Use m_Rupt_Struct

!!$#ifdef PB_2D
!!$  Use m_Rupt2D_Vars
!!$  Use m_Rupt2D_U
!!$  Use m_Rupt2D_V
!!$#else
!!$  Use m_Rupt3D_Vars
!!$  Use m_Rupt3D_U
!!$  Use m_Rupt3D_V
!!$#endif

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
  Subroutine Comp_Bulk_Ener(Ener, ULoc, VLoc, Geom, Params, SD_U, SD_V, Elems_U, Elems_V, Nodes_U, Nodes_V, FLoc, TempLoc)
    Type (EXO_Geom_Info)                         :: Geom
    Type (Rupt_Params)                           :: Params
    Type (SD_Info)                               :: SD_U
    Type (SD_Info)                               :: SD_V
    
#if defined PB_2D
    Type (Node2D), Dimension(:), Pointer          :: Nodes_U
    Type (Node2D), Dimension(:), Pointer          :: Nodes_V
    Type (Element2D_Elast), Dimension(:), Pointer :: Elems_U
    Type (Element2D_Scal), Dimension(:), Pointer  :: Elems_V
    Type(MatS2D)                                  :: Sigma, Epsilon
    Type(Vect2D)                                  :: F, U
#elif defined PB_3D
    Type (Node3D), Dimension(:), Pointer          :: Nodes_U
    Type (Node3D), Dimension(:), Pointer          :: Nodes_V
    Type (Element3D_Elast), Dimension(:), Pointer :: Elems_U
    Type (Element3D_Scal), Dimension(:), Pointer  :: Elems_V
    Type (MatS3D)                                 :: Sigma, Epsilon
    Type (Vect3D)                                 :: F, U
#else
    Type (Node2D), Dimension(:), Pointer          :: Nodes_U
    Type (Node2D), Dimension(:), Pointer          :: Nodes_V
    Type (Element2D_Scal), Dimension(:), Pointer  :: Elems_U
    Type (Element2D_Scal), Dimension(:), Pointer  :: Elems_V
    Type (Vect2D)                                 :: Sigma, Epsilon
#endif

    Vec                                          :: ULoc
    Vec                                          :: VLoc
    Vec                                          :: FLoc
    Vec                                          :: TempLoc
    Real(Kind = Kr), Intent(OUT)                 :: Ener

    PetscReal                                    :: E, Nu
    PetscReal                                    :: K1, K2, K3 

    Real(Kind = Kr)                              :: MyEner
    Integer                                      :: Nb_Gauss, Nb_DoF
    Integer                                      :: iSL, iSG
    Integer                                      :: iBlk, iELoc, iE, iG
    

    Real(Kind = Kr), Dimension(:), Pointer       :: UPtr, VPtr, FPtr, TempPtr
    Real(Kind = Kr)                              :: ContrV

    Integer, Dimension(:), Pointer               :: Loc_Indices_U
    Integer, Dimension(:), Pointer               :: Loc_Indices_V

    Integer                                      :: i

    Allocate(Loc_Indices_V(Geom%Num_Nodes))
    Loc_Indices_V = (/ (i ,i = 0, Geom%Num_Nodes-1) /)
    Call AOApplicationToPETSc(SD_V%Loc_AO, Geom%Num_Nodes, Loc_Indices_V, iErr)

#ifdef PB_2DA
    Allocate(Loc_Indices_U(Geom%Num_Nodes))
    Loc_Indices_U = (/ (i ,i = 0, Geom%Num_Nodes-1) /)
    Call AOApplicationToPETSc(SD_U%Loc_AO, Geom%Num_Nodes, Loc_Indices_U, iErr)
#else
    Allocate(Loc_Indices_U(Geom%Num_Nodes * Geom%Num_Dim))
    Loc_Indices_U = (/ (i ,i = 0, Geom%Num_Nodes * Geom%Num_Dim - 1) /)
    Call AOApplicationToPETSc(SD_U%Loc_AO, Geom%Num_Nodes * Geom%Num_Dim, Loc_Indices_U, iErr)
#endif

    Call VecGetArrayF90(VLoc, VPtr, iErr)
    Call VecGetArrayF90(ULoc, UPtr, iErr)
    Call VecGetArrayF90(FLoc, FPtr, iErr)
    Call VecGetArrayF90(TempLoc, TempPtr, iErr)

    MyEner = 0.0_Kr
    Ener = 0.0_Kr
    Do_iBlk: Do iBlk = 1, Geom%Num_elem_blks
       E         = Params%Young_Mod(iBlk)
       nu        = Params%Poisson_Ratio(iBlk) 

!!! The isotropic Hooke's law is expressed as
!!! \sigma = K1 * trace(Epsilon) Id + 2*K2 * Epsilon - K3 Temp * Id
!!! K1, K2, K3 are computed in terms of E and nu
!!! in 3D, K1 = lambda, K2 = mu, K3 = E*alpha/(1-2nu) (= 3kappa alpha)
!!!       (alpha = therm exp coef).
!!! in 2D / plane stresses, the expressions are more complicated
!!!
#ifdef PB_2D
       K1 = E * nu / (1.0_Kr - nu**2)
       K2 = E / (1.0_Kr + nu) * InvOf2
       K3 = Params%Therm_Exp(iBlk) * E / (1.0_Kr - nu )
#else
       K1 = E * nu / (1.0_Kr - 2.0_Kr * nu) / ( 1.0_Kr + nu)
       K2 = E / (1.0_Kr + nu) * InvOf2
       K3 = Params%Therm_Exp(iBlk) * E / (1.0_Kr - 2.0_Kr * nu)
#endif

#ifdef PB_2DA
       K2 = E / (1.0_Kr + nu) * InvOf2
#endif

       Do_iE: Do iELoc = 1, Geom%Elem_Blk(iBlk)%Num_Elems
          iE = Geom%Elem_Blk(iBlk)%Elem_ID(iELoc)
          If ((.NOT. SD_U%IsLocal_Elem(iE)) .OR. (.NOT. Params%Is_Domain(iBlk))) Then
             CYCLE
          End If

          Call Init_Gauss_EXO(Elems_U, Nodes_U, Geom, MEF90_GaussOrder, Elem=iE)
          Call Init_Gauss_EXO(Elems_V, Nodes_V, Geom, MEF90_GaussOrder, Elem=iE)

          Nb_Gauss = Elems_U(iE)%Nb_Gauss
          Do_iG: Do iG = 1, Nb_Gauss
 
          ! V^2 + K_\epsilon term
             Is_Brittle: If (Params%Is_Brittle(iBlk)) Then
                ContrV = 0.0
                Do_iSLV: Do iSL = 1, Elems_V(iE)%Nb_DoF
                   iSG = Elems_V(iE)%ID_DoF(iSL)
                   
                   ContrV = ContrV + VPtr(Loc_Indices_V(iSG)+1) * Elems_V(iE)%BF(iSL, iG)
                End Do Do_iSLV
                ContrV = ContrV**2 + Params%Kepsilon
             Else
                ContrV = 1.0_Kr
             End If Is_Brittle
             
             
             Sigma = 0.0_Kr
             Epsilon = 0.0_Kr
             Do_iSL1: Do iSL = 1, Elems_U(iE)%Nb_DoF
                iSG = Elems_U(iE)%ID_DoF(iSL)

#ifdef PB_2DA
                Sigma   = Sigma + K2 * Elems_U(iE)%Grad_BF(iSL,iG) * UPtr(Loc_Indices_U(iSG)+1)                
                Epsilon =    Epsilon + Elems_U(iE)%Grad_BF(iSL,iG) * UPtr(Loc_Indices_U(iSG)+1)                
#else
                Sigma    = Sigma + 2.0_Kr * K2 * Elems_U(iE)%GradS_BF(iSL,iG)  * UPtr(Loc_Indices_U(iSG)+1)
                Sigma%XX = Sigma%XX + K1 * Trace(Elems_U(iE)%GradS_BF(iSL,iG)) * UPtr(Loc_Indices_U(iSG)+1)
                Sigma%YY = Sigma%YY + K1 * Trace(Elems_U(iE)%GradS_BF(iSL,iG)) * UPtr(Loc_Indices_U(iSG)+1)
#ifdef PB_3D
                Sigma%ZZ = Sigma%ZZ + K1 * Trace(Elems_U(iE)%GradS_BF(iSL,iG)) * UPtr(Loc_Indices_U(iSG)+1)
#endif
                Epsilon  = Epsilon + Elems_U(iE)%GradS_BF(iSL,iG) * UPtr(Loc_Indices_U(iSG)+1)
#endif
             End Do Do_iSL1
             
#ifndef PB_2DA
!!! Thermal stuff
             Do_iSL2: Do iSL = 1, Elems_V(iE)%Nb_DoF
                iSG = Elems_V(iE)%ID_DoF(iSL)

                Sigma%XX   = Sigma%XX -   K3 * Elems_V(iE)%BF(iSL,iG) * TempPtr(Loc_Indices_V(iSG)+1)
                Sigma%YY   = Sigma%YY -   K3 * Elems_V(iE)%BF(iSL,iG) * TempPtr(Loc_Indices_V(iSG)+1)
                Epsilon%XX = Epsilon%XX - K3 * Elems_V(iE)%BF(iSL,iG) * TempPtr(Loc_Indices_V(iSG)+1)
                Epsilon%YY = Epsilon%YY - K3 * Elems_V(iE)%BF(iSL,iG) * TempPtr(Loc_Indices_V(iSG)+1)
#ifdef PB_3D
                Sigma%ZZ   = Sigma%ZZ -   K3 * Elems_V(iE)%BF(iSL,iG) * TempPtr(Loc_Indices_V(iSG)+1)
                Epsilon%ZZ = Epsilon%ZZ - K3 * Elems_V(iE)%BF(iSL,iG) * TempPtr(Loc_Indices_V(iSG)+1)
#endif
             End Do Do_iSL2

#endif
             MyEner = MyEner + Elems_U(iE)%Gauss_C(iG) * ContrV * (Sigma .DotP. Epsilon) * .5_Kr
             
             
!!! Forces stuff
#ifndef PB_2DA
             U = 0.0_Kr
             F = 0.0_Kr
             If (Params%Has_Force(iBlk)) Then
                Do_iSL3: Do iSL = 1, Elems_U(iE)%Nb_DoF
                   iSG = Elems_U(iE)%ID_DoF(iSL)
                   U = U + UPtr(Loc_Indices_U(iSG)+1) * Elems_U(iE)%BF(iSL,iG)
                   F = F + FPtr(Loc_Indices_U(iSG)+1) * Elems_U(iE)%BF(iSL,iG)
                End Do Do_iSL3
                MyEner = MyEner - Elems_U(iE)%Gauss_C(iG) * (F .DotP. U) 
              End If
#endif
          End Do Do_iG
          Call Destroy_Gauss_EXO(Elems_U, Elem=iE)
          Call Destroy_Gauss_EXO(Elems_V, Elem=iE)
       End Do Do_iE
    End Do Do_iBlk

    Call MPI_Reduce(MyEner, Ener, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, iErr)
    Call MPI_BCAST(Ener, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, iErr)
    Call VecRestoreArrayF90(ULoc, UPtr, iErr)
    Call VecRestoreArrayF90(VLoc, VPtr, iErr)
    Call VecRestoreArrayF90(FLoc, FPtr, iErr)
    Call VecRestoreArrayF90(TempLoc, TempPtr, iErr)
    DeAllocate(Loc_Indices_U)
    DeAllocate(Loc_Indices_V)
  End Subroutine Comp_Bulk_Ener

  Subroutine Comp_Surf_Ener(Ener, VLoc, Geom, Params, SD, Elems_V, Nodes_V)
    Type (EXO_Geom_Info)                         :: Geom
    Type (Rupt_Params)                           :: Params
    Type (SD_Info)                               :: SD

#ifdef PB_3D
    Type (Node3D), Dimension(:), Pointer         :: Nodes_V
    Type (Element3D_Scal), Dimension(:), Pointer :: Elems_V
    Type (Vect3D)                                :: GradV
#else
    Type (Node2D), Dimension(:), Pointer         :: Nodes_V
    Type (Element2D_Scal), Dimension(:), Pointer :: Elems_V
    Type (Vect2D)                                :: GradV
#endif
    Vec                                          :: VLoc
    Real(Kind = Kr), Intent(OUT)                 :: Ener


    Real(Kind = Kr)                              :: MyEner
    Real(Kind = Kr)                              :: V2

    Integer                                      :: Nb_Gauss, Nb_DoF
    Integer                                      :: iSL, iSG
    Integer                                      :: iBlk, iELoc, iE, iG

    Real(Kind = Kr), Dimension(:), Pointer       :: VPtr
    PetscReal                                    :: Toughness

    Integer, Dimension(:), Pointer               :: Loc_Indices_V

    Integer                                      :: SetN, i

    Allocate(Loc_Indices_V(Geom%Num_Nodes))
    Loc_Indices_V = (/ (i ,i = 0, Geom%Num_Nodes - 1) /)
    Call AOApplicationToPETSc(SD%Loc_AO, Geom%Num_Nodes, Loc_Indices_V, iErr)

    Call VecGetArrayF90(VLoc, VPtr, iErr)
    
    MyEner = 0.0_Kr
    Ener   = 0.0_Kr

    Do_iBlk: Do iBlk = 1, Geom%Num_elem_blks
       Toughness = Params%Toughness(iBlk)

       Do_iE: Do iELoc = 1, Geom%Elem_Blk(iBlk)%Num_Elems
          iE = Geom%Elem_Blk(iBlk)%Elem_ID(iELoc)
          If (.NOT. SD%IsLocal_Elem(iE)) Then
             CYCLE
          End If

          Call Init_Gauss_EXO(Elems_V, Nodes_V, Geom, MEF90_GaussOrder, Elem=iE)
          Nb_Gauss = Elems_V(iE)%Nb_Gauss
          
          Do_iG: Do iG = 1, Nb_Gauss
             
             V2 = 0.0_Kr
             GradV = 0.0_Kr
             Do_iSLV: Do iSL = 1, Elems_V(iE)%Nb_DoF
                iSG = Elems_V(iE)%ID_DoF(iSL)
                   
                V2 = V2 + (1.0_Kr - VPtr(Loc_Indices_V(iSG)+1)) * Elems_V(iE)%BF(iSL, iG)
                GradV = GradV + VPtr(Loc_Indices_V(iSG)+1) * Elems_V(iE)%Grad_BF(iSL, iG)
             End Do Do_iSLV
             MyEner = MyEner + Toughness * Elems_V(iE)%Gauss_C(iG) * ( V2**2 * .25_Kr / Params%Epsilon + Params%Epsilon * (GradV .DotP. GradV) )
          End Do Do_iG
          Call Destroy_Gauss_EXO(Elems_V, Elem=iE)
       End Do Do_iE
    End Do Do_iBlk
    
    Call MPI_Reduce(MyEner, Ener, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, iErr)
    Call MPI_BCAST(Ener, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, iErr)
    Call VecRestoreArrayF90(VLoc, VPtr, iErr)
    DeAllocate(Loc_Indices_V)
  End Subroutine Comp_Surf_Ener

#if defined PB_2D
End Module m_Rupt2D_Ener
#elif defined PB_3D
End Module m_Rupt3D_Ener
#else
End Module m_Rupt2DA_Ener
#endif
