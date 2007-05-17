Module m_Rupt2DA_V

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

  Public :: Assemb_MR_V
  Public :: Assemb_RHS_V
  Public :: Apply_BC_V
  Public :: Init_V_Cracks
  Public :: Init_V_Spheres
  
  Integer    :: iErr


Contains
  Function Distance(A, B, M)
! computes the distance from M to the segment [A,B]

    Type (Vect2D)                                :: A, B, M  

    Real(Kind = Kr)                              :: test, normAB

    Real(Kind = Kr)                              :: Distance
    
    normAB = sqrt((B-A).DotP.(B-A))
    test = ( (B-A).DotP.(M-A) ) / normAB
    If( test < 0.0_Kr) Then
       Distance = sqrt( (M-A) .DotP. (M-A) )
    Else If (test > normAB) Then
       Distance = sqrt( (M-B) .DotP. (M-B) )
    Else
       Distance = sqrt( ((M-A) .DotP. (M-A)) - ( ((M-A) .DotP. (B-A)) / normAB ) **2 )
    End If
    
  End Function Distance
  
  Function CrackProfile(t,Params)
    Real(Kind = Kr)                               :: t    
    Type (Rupt_Params)                            :: Params
    Real(Kind = Kr)                               :: CrackProfile
    
    CrackProfile = 1.0_Kr - EXP((Params%Epsilon**2 - t)/(2*Params%Epsilon) )
    CrackProfile= max(CrackProfile , 0.0_Kr)
  End Function  CrackProfile
  
  Subroutine  Init_V_Cracks(Geom, Params, SD, Elems_V,  Nodes_V, V)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params)                            :: Params
    Type (SD_Info)                                :: SD    

    Type (Node2D), Dimension(:), Pointer          :: Nodes_V
    Type (Element2D_Scal), Dimension(:), Pointer  :: Elems_V

    Vec                                           :: V
    
    Integer                                      :: Rand_Node
    Real(Kind = Kr)                              :: Rand_Theta ! colatitude
    Real(Kind = Kr)                              :: Rand_Phi   ! polar angle
    Real(Kind = Kr)                              :: Rand_Length, Length
    Integer                                      :: iTerCracks,Nb_DoF_V
    Integer                                      :: iBlk, iE, iEloc, iS, iSG ,i
    
    Real(Kind = Kr), Dimension(:), Pointer       :: V_Ptr
    Integer, Dimension(:), Pointer               :: Loc_Indices

    Type (Vect2D), Dimension(:),Pointer          :: Crack_Loc

    Real(Kind = Kr)                              :: Tmp_Node
    
    Allocate(Loc_Indices(Geom%Num_Nodes))
    Loc_Indices = (/ (i ,i = 0, Geom%Num_Nodes-1) /)
    Call AOApplicationToPETSc(SD%Loc_AO, Geom%Num_Nodes, Loc_Indices, iErr) 
    
    Allocate(Crack_Loc(2*Params%nbCracks))          
    If (MEF90_MyRank ==0) then  
       Do_Cracks : Do iTerCracks=1, Params%nbCracks
          Call Random_Number(Tmp_Node)
          Rand_Node = int( Tmp_Node * real(Geom%Num_Nodes) + 0.5_Kr)

          Call Random_Number(Rand_Theta)
          Rand_Theta = Rand_Theta * 2.0_Kr * Pi
          Call Random_Number(Rand_Length)
!!! @altech: fix the size of the initial cracks
!!!          Rand_Length = Rand_Length * Params%MaxCrackLength
!!!
          Rand_length = Params%MaxCrackLength

          Crack_Loc(iTerCracks*2 - 1) = Nodes_V(Rand_Node)%Coord
          Crack_Loc(iTerCracks*2 )    = Crack_Loc(iTerCracks*2 - 1) 

          Crack_Loc(iTerCracks*2 )%X  = Crack_Loc(iTerCracks*2)%X + Rand_Length * COS(Rand_Theta)  
          Crack_Loc(iTerCracks*2 )%Y  = Crack_Loc(iTerCracks*2)%Y + Rand_Length * SIN(Rand_Theta)

       End Do Do_Cracks
    End If
    

    Call MPI_Bcast(Crack_Loc, 2*Params%nbCracks, Vect2D_MPIType, 0, PETSC_COMM_WORLD,iErr)

    
    Do_iBlk: Do iBlk = 1, Geom%Num_elem_blks
!       Is_Brittle: If ( Params%Is_Brittle(iBlk)) Then
       Nb_Dof_V = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem
       Do_iE: Do iELoc = 1, Geom%Elem_Blk(iBlk)%Num_Elems
          iE = Geom%Elem_Blk(iBlk)%ELem_ID(iELoc)
          If (.NOT. SD%IsLocal_Elem(iE)) Then
             CYCLE
          End If
          Call VecGetArrayF90(V, V_Ptr, iErr)        
          Do_iS: Do iS = 1, Elems_V(iE)%Nb_DoF
             iSG = Elems_V(iE)%ID_DoF(iS)
             If (SD%IsLocal_Node(iSG)) Then
                Do_Crack : Do iTerCracks=1, Params%nbCracks
                   Length = Distance(Crack_Loc(iTerCracks*2 -1), Crack_Loc(iTerCracks*2 ), Nodes_V(iSG)%Coord)
                   V_Ptr(Loc_Indices(iSG)+1) = min(V_Ptr(Loc_Indices(iSG)+1), CrackProfile(Length,Params))
                End Do Do_Crack
             End IF
          End Do Do_iS
          Call VecRestoreArrayF90(V, V_Ptr, iErr)
       End Do Do_iE
!       End If  Is_Brittle
    End Do Do_iBlk
    
    Call VecGhostUpdateBegin(V, INSERT_VALUES, SCATTER_FORWARD, iErr)
    Call VecGhostUpdateEnd(V, INSERT_VALUES, SCATTER_FORWARD, iErr)
    
    DeAllocate(Loc_Indices)
    DeAllocate(Crack_Loc)
  End Subroutine Init_V_Cracks

  Subroutine  Init_V_Spheres(Geom, Params, SD, Elems_V,  Nodes_V, V)
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params)                            :: Params
    Type (SD_Info)                                :: SD    

    Type (Node2D), Dimension(:), Pointer          :: Nodes_V
    Type (Element2D_Scal), Dimension(:), Pointer  :: Elems_V

    Vec                                           :: V
    
    Integer                                      :: Rand_Node
    Real(Kind = Kr)                              :: Length
    Integer                                      :: iTerCracks,Nb_DoF_V
    Integer                                      :: iBlk, iE, iEloc, iS, iSG ,i
    
    Real(Kind = Kr), Dimension(:), Pointer       :: V_Ptr
    Integer, Dimension(:), Pointer               :: Loc_Indices

    Type (Vect2D), Dimension(:),Pointer          :: Crack_Loc

    Real(Kind = Kr)                              :: Tmp_Node
    
    Allocate(Loc_Indices(Geom%Num_Nodes))
    Loc_Indices = (/ (i ,i = 0, Geom%Num_Nodes-1) /)
    Call AOApplicationToPETSc(SD%Loc_AO, Geom%Num_Nodes, Loc_Indices, iErr) 
    
    Allocate(Crack_Loc(Params%nbCracks))          
    If (MEF90_MyRank ==0) then  
       Do_Cracks : Do iTerCracks=1, Params%nbCracks
          Call Random_Number(Tmp_Node)
          Rand_Node = int( Tmp_Node * real(Geom%Num_Nodes) + 0.5_Kr)
          Crack_Loc(iTerCracks) = Nodes_V(Rand_Node)%Coord
       End Do Do_Cracks
    End If
    

    Call MPI_Bcast(Crack_Loc, Params%nbCracks, Vect2D_MPIType, 0, PETSC_COMM_WORLD,iErr)

    
    Do_iBlk: Do iBlk = 1, Geom%Num_elem_blks
       Nb_Dof_V = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem
       Do_iE: Do iELoc = 1, Geom%Elem_Blk(iBlk)%Num_Elems
          iE = Geom%Elem_Blk(iBlk)%ELem_ID(iELoc)
          If (.NOT. SD%IsLocal_Elem(iE)) Then
             CYCLE
          End If
          Call VecGetArrayF90(V, V_Ptr, iErr)        
          Do_iS: Do iS = 1, Elems_V(iE)%Nb_DoF
             iSG = Elems_V(iE)%ID_DoF(iS)
             If (SD%IsLocal_Node(iSG)) Then
                Do_Crack : Do iTerCracks=1, Params%nbCracks
                   Length = sqrt((Crack_Loc(iTerCracks) - Nodes_V(iSG)%Coord) .DotP. (Crack_Loc(iTerCracks) - Nodes_V(iSG)%Coord))
                   Length = Max(0.0_Kr, Length - Params%MaxCrackLength)
                   V_Ptr(Loc_Indices(iSG)+1) = min(V_Ptr(Loc_Indices(iSG)+1), CrackProfile(Length,Params)+.1_Kr)
                End Do Do_Crack
             End IF
          End Do Do_iS
          Call VecRestoreArrayF90(V, V_Ptr, iErr)
       End Do Do_iE
    End Do Do_iBlk
    
    Call VecGhostUpdateBegin(V, INSERT_VALUES, SCATTER_FORWARD, iErr)
    Call VecGhostUpdateEnd(V, INSERT_VALUES, SCATTER_FORWARD, iErr)
    
    DeAllocate(Loc_Indices)
    DeAllocate(Crack_Loc)
  End Subroutine Init_V_Spheres



!!! USING VLV SEEMS TO MAKE THE KSP UNSTABLE, UNLESS V IS INITIALIZED WITH
!!! THE PROPER BC
  Subroutine Assemb_MR_V(MR, U_Loc, Geom, Params, SD_U, SD_V, Elems_U, Elems_V, Nodes_U, Nodes_V, t)
    Mat                                           :: MR
    Vec                                           :: U_Loc
    Type (EXO_Geom_Info)                          :: Geom
    Type (Rupt_Params)                            :: Params
    Type (SD_Info)                                :: SD_U
    Type (SD_Info)                                :: SD_V
	 Real(Kind = Kr), Intent(IN)                   :: t



    Type (Node2D), Dimension(:), Pointer          :: Nodes_U
    Type (Node2D), Dimension(:), Pointer          :: Nodes_V
    Type (Element2D_Scal), Dimension(:), Pointer  :: Elems_U
    Type (Element2D_Scal), Dimension(:), Pointer  :: Elems_V



    Integer                                       :: Nb_Gauss
    Integer                                       :: Nb_DoF_U, Nb_DoF_V
    Integer                                       :: iSLV1
    Integer                                       :: iSLV2
    Integer                                       :: iSL, iSG
    Integer                                       :: iSLSig, iSGSig
    Integer                                       :: iSLEps, iSGEps
    Integer                                       :: iBlk, iE, iEloc, iG
                                                  
    Real(Kind = Kr), Dimension(:), Pointer        :: UPtr
    Real(Kind = Kr)                               :: ContrU
                                                  
    PetscScalar, Dimension(:,:), Pointer          :: MR_Elem
    Integer, Dimension(:), Pointer                :: EXO_Indices_V
    Integer, Dimension(:), Pointer                :: Loc_Indices_U
    Integer, Dimension(:), Pointer                :: Loc_Indices_V
                                                  
    PetscTruth                                    :: ISAssembled
                                                  
    Real(Kind = Kr)                               :: Gc
    Type(Vect2D)                                  :: Sigma, Distortion
    Integer                                       :: i

    Call MatAssembled(MR, IsAssembled, iErr)
    If (IsAssembled) Then
       Call MatZeroEntries(MR, iErr)
    End If

    Allocate(Loc_Indices_U(Geom%Num_Nodes))
    Loc_Indices_U = (/ (i ,i = 0, Geom%Num_Nodes-1) /)
    Call AOApplicationToPETSc(SD_U%Loc_AO, Geom%Num_Nodes, Loc_Indices_U, iErr)

    Allocate(Loc_Indices_V(Geom%Num_Nodes))
    Loc_Indices_V = (/ (i ,i = 0, Geom%Num_Nodes-1) /)
    Call AOApplicationToPETSc(SD_V%Loc_AO, Geom%Num_Nodes, Loc_Indices_V, iErr)

    Call VecGetArrayF90(U_Loc, UPtr, iErr) 

    Do_iBlk: Do iBlk = 1, Geom%Num_elem_blks
       Gc  = Params%Toughness(iBlk)

       Nb_DoF_U = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem
	    Nb_Dof_V = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem

       Allocate (MR_Elem(Nb_DoF_V, Nb_DoF_V))
       Allocate (EXO_Indices_V(Nb_DoF_V))

       Do_iE: Do iELoc = 1, Geom%Elem_Blk(iBlk)%Num_Elems
          iE = Geom%Elem_Blk(iBlk)%ELem_ID(iELoc)
          If (.NOT. SD_V%IsLocal_Elem(iE)) Then
             CYCLE
          End If
          EXO_Indices_V = Elems_V(iE)%Id_DoF-1
          Call AOApplicationToPETSc(SD_V%EXO_AO, Nb_DoF_V, EXO_Indices_V, iErr)
       
          Call Init_Gauss_EXO(Elems_U, Nodes_U, Geom, MEF90_GaussOrder, Elem=iE)
          Call Init_Gauss_EXO(Elems_V, Nodes_V, Geom, MEF90_GaussOrder, Elem=iE)
          MR_Elem = 0.0_Kr
          
          Do_iG: Do iG = 1, Elems_U(iE)%Nb_Gauss

!!! part related to v^2+k_\e W(e(u))
             ContrU = 0.0_Kr
             Is_Brittle: If ( Params%Is_Brittle(iBlk)) Then
                Sigma = 0.0_Kr
                Distortion = 0.0_Kr
                Do_iSL1: Do iSL = 1, Elems_U(iE)%Nb_DoF
                   iSG = Elems_U(iE)%ID_DoF(iSL)
                   Sigma        =  Sigma + Elems_U(iE)%Grad_BF(iSL,iG) * UPtr(Loc_Indices_U(iSG)+1)                
                   Distortion%X = Distortion%X + Elems_U(iE)%BF(iSL,iG) * Nodes_U(iSG)%Coord%Y 
                   Distortion%Y = Distortion%Y - Elems_U(iE)%BF(iSL,iG) * Nodes_U(iSG)%Coord%X
                End Do Do_iSL1
                ContrU = (Sigma - t * Distortion) .DotP. (sigma - t * Distortion)
             
                DoiSLV1: Do iSLV1 = 1, Elems_V(iE)%Nb_DoF
                   DoiSLV2: Do iSLV2 = 1, Elems_V(iE)%Nb_DoF
                      MR_Elem(iSLV1, iSLV2) = MR_Elem(iSLV1, iSLV2) + Elems_V(iE)%Gauss_C(iG) * ContrU *                           &
   		                                     Elems_V(iE)%BF(iSLV1, iG) * Elems_V(iE)%BF(iSLV2, iG)
                   End Do DoiSLV2
                End Do DoiSLV1
             End If Is_Brittle

!!! Surface energy part
             DoiSLV1surf: Do iSLV1 = 1, Elems_V(iE)%Nb_DoF
                DoiSLV2surf: Do iSLV2 = 1, Elems_V(iE)%Nb_DoF
                   MR_Elem(iSLV1, iSLV2) = MR_Elem(iSLV1, iSLV2) + Elems_V(iE)%Gauss_C(iG) * Gc *                                  &
	                                        ( (Elems_V(iE)%BF(iSLV1, iG) * Elems_V(iE)%BF(iSLV2, iG)) / Params%Epsilon * .5_Kr +    &
				   			                     (Elems_V(iE)%Grad_BF(iSLV1, iG) .DotP. Elems_V(iE)%Grad_BF(iSLV2, iG))                &
				   			                      * Params%Epsilon * 2.0_Kr )
                End Do DoiSLV2surf
             End Do DoiSLV1surf
          End Do Do_iG
          
          If (Params%Do_Irrev) Then
             Do_iSLBC: Do iSL = 1, Elems_V(iE)%Nb_DoF
                iSG = Elems_V(iE)%ID_DoF(iSL)
                If (Nodes_V(iSG)%BC /= BC_TYPE_NONE) Then
                   MR_Elem(iSL,:)    = 0.0_Kr
                   MR_Elem(:,iSL)    = 0.0_Kr
                   MR_Elem(iSL, iSL) = 1.0_Kr
                End If
             End Do Do_iSLBC
          End If
          
          Call MatSetValues(MR, Nb_DoF_V, EXO_Indices_V, Nb_DoF_V, EXO_Indices_V, MR_Elem, ADD_VALUES, iErr)
          
          Call Destroy_Gauss_EXO(Elems_U, Elem=iE)
          Call Destroy_Gauss_EXO(Elems_V, Elem=iE)
       EndDo Do_iE
       DeAllocate(MR_Elem)
       DeAllocate(EXO_Indices_V)
    End Do Do_iBlk
    
    DeAllocate(Loc_Indices_U)
    DeAllocate(Loc_Indices_V)
    
    Call VecRestoreArrayF90(U_Loc, UPtr, iErr) 

    Call MatAssemblyBegin(MR, MAT_FINAL_ASSEMBLY, iErr)
    Call MatAssemblyEnd(MR, MAT_FINAL_ASSEMBLY, iErr)  
  End Subroutine Assemb_MR_V


  Subroutine Assemb_RHS_V(RHS, Geom, Params, SD, Elems_Scal, Nodes_Scal)
    Vec                                          :: RHS
    Type (EXO_Geom_Info)                         :: Geom
    Type (Rupt_Params)                           :: Params
    Type (SD_Info)                               :: SD


    Type (Node2D), Dimension(:), Pointer         :: Nodes_Scal
    Type (Element2D_Scal), Dimension(:), Pointer :: Elems_Scal


    PetscReal, Dimension(:), Pointer             :: RHS_Ptr


    Integer                                      :: Nb_DoF
    Integer                                      :: iSL
    Integer                                      :: iE, iELoc, iG
    Integer                                      :: iSLoc, iS, iBlk

    Integer, Dimension(:), Pointer               :: EXO_Indices
    Integer, Dimension(:), Pointer               :: Loc_Indices

    Integer                                      :: i

    Real(Kind = Kr)                              :: Gc
    Real(Kind = Kr), Dimension(:), Pointer       :: RHS_Elem
    
    
    Call VecSet(RHS, 0.0_Kr, iErr)

    Do_iBlk: Do iBlk = 1, Geom%Num_elem_blks
       Gc = Params%Toughness(iBlk)
       Nb_DoF = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem

       Allocate(RHS_Elem(Nb_DoF))
       Allocate(EXO_Indices(Nb_DoF))

       Do_iE: Do iELoc = 1, Geom%Elem_blk(iBlk)%Num_Elems
          iE = Geom%Elem_blk(iBlk)%Elem_ID(iELoc)
          If (.NOT. SD%IsLocal_Elem(iE)) Then
             CYCLE
          End If
          EXO_Indices = Elems_Scal(iE)%Id_DoF-1
          Call AOApplicationToPETSc(SD%EXO_AO, Nb_DoF, EXO_Indices, iErr)

          Call Init_Gauss_EXO(Elems_Scal, Nodes_Scal, Geom, MEF90_GaussOrder, Elem=iE)
          RHS_Elem = 0.0_Kr

          Do_iG: Do iG = 1, Elems_Scal(iE)%Nb_Gauss
             Do_iSL: Do iSL = 1, Nb_DoF
                RHS_Elem(iSL) = RHS_Elem(iSL) + Gc * Elems_Scal(iE)%Gauss_C(iG) * Elems_Scal(iE)%BF(iSL, iG)                       &
                                / Params%Epsilon * .5_Kr
             End Do Do_iSL
          End Do Do_iG
          Call VecSetValues(RHS, Nb_DoF, EXO_Indices, RHS_Elem, ADD_VALUES, iErr)
          Call Destroy_Gauss_EXO(Elems_Scal, Elem=iE)
       End Do Do_iE

       DeAllocate (RHS_Elem)
       DeAllocate (EXO_Indices)
    End Do Do_iBlk
    
    Call VecAssemblyBegin(RHS, iErr)
    Call VecAssemblyEnd(RHS, iErr)

! BC Due to irreversibility:
    Is_Irrev: If (Params%Do_Irrev) Then
       Allocate(Loc_Indices(Geom%Num_Nodes))
       Loc_Indices = (/ (i ,i = 0, Geom%Num_Nodes - 1) /)
       Call AOApplicationToPETSc(SD%Loc_AO, Geom%Num_Nodes, Loc_Indices, iErr)

       Call VecGetArrayF90(RHS, RHS_Ptr, iErr)
       Do_iS: Do iSLoc = 1, SD%Num_Nodes
          iS = SD%Node(iSLoc)
          Is_BC: If (Nodes_Scal(iS)%BC /= BC_Type_NONE) Then
             RHS_Ptr(Loc_Indices(iS)+1) = 0.0_Kr
          End If Is_BC
       End Do Do_iS
       Call VecRestoreArrayF90(RHS, RHS_Ptr, iErr)
       DeAllocate(Loc_Indices)
    End If Is_Irrev
  End Subroutine Assemb_RHS_V



  Subroutine Apply_BC_V(Geom, Params, SD, Nodes, V)
!!! CHECK THAT OUT. I am tired...
    Type(EXO_Geom_Info), Intent(IN)                  :: Geom
    Type(Rupt_Params), Intent(IN)                    :: Params
    Type (SD_Info), Intent(IN)                       :: SD

    Type(Node2D), Dimension(:), Pointer              :: Nodes

    Vec                                              :: V

    PetscReal, Dimension(:), Pointer                 :: VPtr
    Integer                                          :: i, iSloc, iS
    Integer, Dimension(:), Pointer                   :: Loc_Indices

    Allocate(Loc_Indices(Geom%Num_Nodes))
    Loc_Indices = (/ (i ,i = 0, Geom%Num_Nodes - 1) /)
    Call AOApplicationToPETSc(SD%Loc_AO, Geom%Num_Nodes, Loc_Indices, iErr)
    Call VecGetArrayF90(V, VPtr, iErr)

    Do iSLoc = 1, SD%Num_Nodes
       iS = SD%Node(iSLoc)
       If (Nodes(iS)%BC == BC_Type_DIRI) Then
          VPtr(Loc_Indices(iS)+1) = 0.0_Kr
       End If 
    End Do

    Call VecRestoreArrayF90(V, VPtr, iErr)
    DeAllocate(Loc_Indices)
  End Subroutine Apply_BC_V
End Module m_Rupt2DA_V

