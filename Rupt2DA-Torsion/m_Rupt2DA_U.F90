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
    Integer                                       :: iSL, iSG
    Integer                                       :: iBlk, iEloc, iE, iG


    Real(Kind = Kr), Dimension(:), Pointer        :: VPtr
    Real(Kind = Kr)                               :: ContrV

    PetscScalar, Dimension(:,:), Pointer          :: MR_Elem
    Integer, Dimension(:), Pointer                :: EXO_Indices_Vect
    Integer, Dimension(:), Pointer                :: Loc_Indices_Scal

    PetscTruth                                    :: ISAssembled

    Integer                                       :: i

    Call MatAssembled(MR, IsAssembled, iErr)
    If (IsAssembled) Then
       Call MatZeroEntries(MR, iErr)
    End If

    
    Allocate(Loc_Indices_Scal(Geom%Num_Nodes))
    Loc_Indices_Scal = (/ (i ,i = 0, Geom%Num_Nodes-1) /)
    Call AOApplicationToPETSc(SD_V%Loc_AO, Geom%Num_Nodes, Loc_Indices_Scal, iErr)

    Call VecGetArrayF90(V_Loc, VPtr, iErr) 

    Do_iBlk: Do iBlk = 1, Geom%Num_elem_blks
       Nb_DoF_U = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem
       Allocate (MR_Elem(Nb_DoF_U, Nb_DoF_U))
       Allocate (EXO_Indices_Vect(Nb_DoF_U))

       Do_iE: Do iELoc = 1, Geom%Elem_Blk(iBlk)%Num_Elems
          iE = Geom%Elem_Blk(iBlk)%ELem_ID(iELoc)
          If ((.NOT. SD_U%IsLocal_Elem(iE)) .OR. (.NOT. Params%Is_Domain(iBlk))) Then
             CYCLE
          End If
          EXO_Indices_Vect = Elems_U(iE)%Id_DoF-1
          Call AOApplicationToPETSc(SD_U%EXO_AO, Nb_DoF_U, EXO_Indices_Vect, iErr)

          Call Init_Gauss_EXO(Elems_U, Nodes_U, Geom, MEF90_GaussOrder, Elem=iE)
          Call Init_Gauss_EXO(Elems_V, Nodes_V, Geom, MEF90_GaussOrder, Elem=iE)

          MR_Elem = 0.0_Kr
          Nb_Gauss = Elems_U(iE)%Nb_Gauss
          
          Do_iG: Do iG = 1, Nb_Gauss

          ! V^2 + K_\epsilon term
             Is_Brittle: If (Params%Is_Brittle(iBlk)) Then
                ContrV = 0.0_Kr
                Do_iSLV: Do iSL = 1, Elems_V(iE)%Nb_DoF
                   iSG = Elems_V(iE)%ID_DoF(iSL)
                   
                   ContrV = ContrV + VPtr(Loc_Indices_Scal(iSG)+1) * Elems_V(iE)%BF(iSL, iG)
                End Do Do_iSLV
                ContrV = ContrV**2 + Params%Kepsilon
             Else
                ContrV = 1.0_Kr
             End If Is_Brittle

             Do_iSLSig: Do iSLSig = 1, Elems_U(iE)%Nb_DoF
                Do_iSLEps: Do iSLEps = 1, Elems_U(iE)%Nb_DoF
                   MR_Elem(iSLEps, iSLSig) = MR_Elem(iSLEps, iSLSig) +  Elems_U(iE)%Gauss_C(iG) * ContrV     *                     &
                                           ( Elems_U(iE)%Grad_BF(iSLEps,iG) .DotP. Elems_U(iE)%Grad_BF(iSLSig,iG) )
                End Do Do_iSLEps
             End Do Do_iSLSig
   
          End Do Do_iG
          Call MatSetValues(MR, Nb_DoF_U, EXO_Indices_Vect, Nb_DoF_U, EXO_Indices_Vect, MR_Elem, ADD_VALUES, iErr)
          
          Call Destroy_Gauss_EXO(Elems_U, Elem=iE)
          Call Destroy_Gauss_EXO(Elems_V, Elem=iE)
       EndDo Do_iE
       DeAllocate(MR_Elem)
       DeAllocate(EXO_Indices_Vect)
    End Do Do_iBlk

    Call MatAssemblyBegin(MR, MAT_FLUSH_ASSEMBLY, iErr)
    DeAllocate(Loc_Indices_Scal)

    Call VecRestoreArrayF90(V_Loc, VPtr, iErr) 
    Call MatAssemblyEnd(MR, MAT_FLUSH_ASSEMBLY, iErr)  


! Assembly of the BC terms

    Allocate(EXO_Indices_Vect(Geom%Num_Nodes))
    EXO_Indices_Vect = (/ (i ,i = 0, Geom%Num_Nodes - 1) /)
    Call AOApplicationToPETSc(SD_U%EXO_AO, Geom%Num_Nodes, EXO_Indices_Vect, iErr)

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
                Call MatSetValue(MR, EXO_Indices_Vect(iSGSig), EXO_Indices_Vect(iSGSig), MEF90_VLV, INSERT_VALUES, iErr)
             End If Is_BC_BC
          End Do Do_iSLSig_BC
       End Do Do_iE_BC
    End Do Do_iBlk_BC
    DeAllocate (EXO_Indices_Vect)

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

    Real(Kind = Kr), Dimension(:), Pointer        :: VPtr

    Integer, Dimension(:), Pointer                :: Loc_Indices_Vect
    Integer, Dimension(:), Pointer                :: EXO_Indices_Vect
    Integer, Dimension(:), Pointer                :: Loc_Indices_Scal
    PetscReal, Dimension(:), Pointer              :: RHS_Ptr
    PetscReal, Dimension(:), Pointer              :: BC_U_Ptr
    Integer                                       :: i, iBlk
    Integer                                       :: iE, iELoc
    Integer                                       :: iSL, iSG
    Integer                                       :: iS, iSLoc
    Integer                                       :: iG
    Integer                                       :: Nb_DoF_U
    Real(Kind = Kr)                               :: ContrV
    Real(Kind = Kr), Dimension(:), Pointer        :: RHS_Elem
    Type(Vect2D)                                  :: Distortion

    Allocate(Loc_Indices_Scal(Geom%Num_Nodes))
    Loc_Indices_Scal = (/ (i ,i = 0, Geom%Num_Nodes - 1) /)
    Call AOApplicationToPETSc(SD_V%Loc_AO, Geom%Num_Nodes, Loc_Indices_Scal, iErr)

    Allocate(Loc_Indices_Vect(Geom%Num_Nodes))
    Loc_Indices_Vect = (/ (i ,i = 0, Geom%Num_Nodes - 1) /)
    Call AOApplicationToPETSc(SD_U%Loc_AO, Geom%Num_Nodes, Loc_Indices_Vect, iErr)
    
    Allocate(EXO_Indices_Vect(Geom%Num_Nodes))
    EXO_Indices_Vect = (/ (i ,i = 0, Geom%Num_Nodes - 1) /)
    Call AOApplicationToPETSc(SD_U%EXO_AO, Geom%Num_Nodes, EXO_Indices_Vect, iErr)
    
    Call VecGetArrayF90(VLoc, VPtr, iErr)
    
    Call VecSet(RHS, 0.0_Kr, iErr)
    Do_iBlk: Do iBlk = 1, Geom%Num_elem_blks
       Nb_DoF_U = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem
       Allocate (RHS_Elem(Nb_DoF_U))
       Allocate (EXO_Indices_Vect(Nb_DoF_U))


       Do_iE: Do iELoc = 1, Geom%Elem_blk(iBlk)%Num_Elems
          iE = Geom%Elem_blk(iBlk)%Elem_ID(iELoc)
          If (.NOT. SD_U%IsLocal_Elem(iE)) Then
             CYCLE
          End If
          EXO_Indices_Vect = Elems_U(iE)%Id_DoF-1
          Call AOApplicationToPETSc(SD_U%EXO_AO, Nb_DoF_U, EXO_Indices_Vect, iErr)
          
          Call Init_Gauss_EXO(Elems_U, Nodes_U, Geom, MEF90_GaussOrder, Elem=iE)
          Call Init_Gauss_EXO(Elems_V, Nodes_V, Geom, MEF90_GaussOrder, Elem=iE)
          
          RHS_Elem = 0.0_Kr
          Do_iG: Do iG = 1, Elems_U(iE)%Nb_Gauss
 
          ! V^2 + K_\epsilon term
             Is_Brittle: If (Params%Is_Brittle(iBlk)) Then
                ContrV = 0.0_Kr
                Do_iSLV: Do iSL = 1, Elems_V(iE)%Nb_DoF
                   iSG = Elems_V(iE)%ID_DoF(iSL)
                   
                   ContrV = ContrV + VPtr(Loc_Indices_Scal(iSG)+1) * Elems_V(iE)%BF(iSL, iG)
                End Do Do_iSLV
                ContrV = ContrV**2 + Params%Kepsilon
             Else
                ContrV = 1.0_Kr
             End If Is_Brittle
             
             Distortion = 0.0_Kr
             Do_iSL1: Do iSL = 1, Elems_U(iE)%Nb_DoF
                iSG = Elems_U(iE)%ID_DoF(iSL)

                Distortion%X = Distortion%X + Elems_U(iE)%BF(iSL,iG) * Nodes_U(iSG)%Coord%Y 
                Distortion%Y = Distortion%Y - Elems_U(iE)%BF(iSL,iG) * Nodes_U(iSG)%Coord%X
             End Do Do_iSL1

             Do_iSL2: Do iSL = 1, Elems_U(iE)%Nb_DoF
                RHS_Elem(iSL) = RHS_Elem(iSL) + Elems_U(iE)%Gauss_C(iG) * ContrV * t *                                       &
                            ( Distortion .DotP. Elems_U(iE)%Grad_BF(iSL, iG) )
             End Do Do_iSL2
          End Do Do_iG
          Call VecSetValues(RHS, Elems_U(iE)%Nb_DoF, EXO_Indices_Vect, RHS_Elem, ADD_VALUES, iErr)
          
          Call Destroy_Gauss_EXO(Elems_U, Elem=iE)
          Call Destroy_Gauss_EXO(Elems_V, Elem=iE)
       End Do Do_iE
       DeAllocate (RHS_Elem)
       DeAllocate(EXO_Indices_Vect)

    End Do Do_iBlk
    Call VecAssemblyBegin(RHS, iErr)
    Call VecAssemblyEnd(RHS, iErr)
    DeAllocate(Loc_Indices_Scal)
    
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
    
    DeAllocate(Loc_Indices_Vect)
    
    Call VecRestoreArrayF90(VLoc, VPtr, iErr)
    Call VecRestoreArrayF90(BC_U_Loc, BC_U_Ptr, iErr)
    Call VecRestoreArrayF90(RHS, RHS_Ptr, iErr)
  End Subroutine Assemb_RHS_U
  
End Module m_Rupt2DA_U
