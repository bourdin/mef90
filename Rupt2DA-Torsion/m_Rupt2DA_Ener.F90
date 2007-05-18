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
    Type (EXO_Geom_Info), Intent(IN)             :: Geom
    Type (Rupt_Params), Intent(IN)               :: Params
    Type (SD_Info), Intent(IN)                   :: SD_U
    Type (SD_Info), Intent(IN)                   :: SD_V
    Real(Kind = Kr), Intent(IN)                  :: t
    
    Type (Node2D), Dimension(:), Pointer          :: Nodes_U
    Type (Node2D), Dimension(:), Pointer          :: Nodes_V
    Type (Element2D_Scal), Dimension(:), Pointer  :: Elems_U
    Type (Element2D_Scal), Dimension(:), Pointer  :: Elems_V
    Type (Vect2D)                                 :: Sigma, Distortion

    Vec                                          :: ULoc
    Vec                                          :: VLoc
    Real(Kind = Kr), Intent(OUT)                 :: Ener

    Real(Kind = Kr)                              :: MyEner
    Integer                                      :: Nb_Gauss, Nb_DoF
    Integer                                      :: iSL, iSG
    Integer                                      :: iBlk, iELoc, iE, iG
    

    Real(Kind = Kr), Dimension(:), Pointer       :: UPtr, VPtr
    Real(Kind = Kr)                              :: ContrV

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
                ContrV = 0.0_Kr
                Do_iSLV: Do iSL = 1, Elems_V(iE)%Nb_DoF
                   iSG = Elems_V(iE)%ID_DoF(iSL)
                   
                   ContrV = ContrV + VPtr(Loc_Indices_V(iSG)+1) * Elems_V(iE)%BF(iSL, iG)
                End Do Do_iSLV
                ContrV = ContrV**2 + Params%Kepsilon
             Else
                ContrV = 1.0_Kr
             End If Is_Brittle
             
             
             Sigma = 0.0_Kr
             Distortion = 0.0_Kr
             Do_iSL1: Do iSL = 1, Elems_U(iE)%Nb_DoF
                iSG = Elems_U(iE)%ID_DoF(iSL)

                Sigma        =  Sigma + Elems_U(iE)%Grad_BF(iSL,iG) * UPtr(Loc_Indices_U(iSG)+1)                
                Distortion%X = Distortion%X + Elems_U(iE)%BF(iSL,iG) * Nodes_U(iSG)%Coord%Y 
                Distortion%Y = Distortion%Y - Elems_U(iE)%BF(iSL,iG) * Nodes_U(iSG)%Coord%X
                
             End Do Do_iSL1
             MyEner = MyEner + Elems_U(iE)%Gauss_C(iG) * ContrV * ( (Sigma - t * Distortion) .DotP. (sigma - t * Distortion) )     &
                      * .5_Kr
              
          End Do Do_iG
          Call Destroy_Gauss_EXO(Elems_U, Elem=iE)
          Call Destroy_Gauss_EXO(Elems_V, Elem=iE)
       End Do Do_iE
    End Do Do_iBlk

    Call MPI_Reduce(MyEner, Ener, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, iErr)
    Call MPI_BCAST(Ener, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, iErr)
    Call VecRestoreArrayF90(ULoc, UPtr, iErr)
    Call VecRestoreArrayF90(VLoc, VPtr, iErr)
    DeAllocate(Loc_Indices_U)
    DeAllocate(Loc_Indices_V)
  End Subroutine Comp_Bulk_Ener

  Subroutine Comp_Surf_Ener(Ener, VLoc, Geom, Params, SD, Elems_V, Nodes_V)
    Type (EXO_Geom_Info)                         :: Geom
    Type (Rupt_Params)                           :: Params
    Type (SD_Info)                               :: SD

    Type (Node2D), Dimension(:), Pointer         :: Nodes_V
    Type (Element2D_Scal), Dimension(:), Pointer :: Elems_V
    Type (Vect2D)                                :: GradV

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
             MyEner = MyEner + Toughness * Elems_V(iE)%Gauss_C(iG) *                                                               &
                      ( V2**2 * .25_Kr / Params%Epsilon + Params%Epsilon * (GradV .DotP. GradV) )
          End Do Do_iG
          Call Destroy_Gauss_EXO(Elems_V, Elem=iE)
       End Do Do_iE
    End Do Do_iBlk
    
    Call MPI_Reduce(MyEner, Ener, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, iErr)
    Call MPI_BCAST(Ener, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, iErr)
    Call VecRestoreArrayF90(VLoc, VPtr, iErr)
    DeAllocate(Loc_Indices_V)
  End Subroutine Comp_Surf_Ener

End Module m_Rupt2DA_Ener

