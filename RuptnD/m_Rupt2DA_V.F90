#if defined PB_2D
Module m_Rupt2D_V
#elif defined PB_3D
Module m_Rupt3D_V
#else
Module m_Rupt2DA_V
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
#include "include/finclude/petscsys.h"

  Public :: Assemb_MR_V
  Public :: Assemb_RHS_V
  Public :: Apply_BC_V
  Public :: Init_V_Cracks
  
  Integer    :: iErr


Contains
!!$  Function ran1()  
!!$    !returns random number between 0 - 1
!!$    implicit none
!!$    Real(Kind = Kr)              ::  ran1, x
!!$    
!!$    call random_number(x) 
!!$    ran1=x
!!$    
!!$  End Function ran1
!!$  
!!$  Function Spread(min,max)  
!!$    !returns random number between min - max
!!$    implicit none
!!$    Real(Kind = Kr) spread
!!$    Real(Kind = Kr) min,max
!!$    
!!$    Spread=(max - min) * ran1() + min
!!$    
!!$  End Function Spread
!!$  
!!$  Function Spread_Int(min,max)  
!!$    !returns random Int_number between min - max
!!$    implicit none
!!$    Integer Spread_Int
!!$    Integer min,max
!!$    
!!$    Spread_Int=(max - min) * ran1() + min
!!$    
!!$  End Function Spread_Int
  
  Function Distance(A, B, M)
    ! computes the distance from M to the segment [A,B]
#if defined PB_3D
    Type (Vect3D)                                :: A, B, M
#else
    Type (Vect2D)                                :: A, B, M  
#endif    
    Real(Kind = Kr)                              :: test, normAB

    Real(Kind = Kr)                              :: Distance
    
    normAB = sqrt((B-A).DotP.(B-A))
    test = ( (B-A).DotP.(M-A) ) / normAB
    If( test < 0.0_Kr) Then
       Distance = sqrt( (M-A) .DotP. (M-A) )
    Else If (test > normAB) Then
       Distance = sqrt( (M-B) .DotP. (M-B) )
    Else
       Distance = sqrt( ((M-A) .DotP. (M-A)) -                             &
            & ( ((M-A) .DotP. (B-A)) / normAB ) **2 )
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
#if defined PB_2D
    Type (Node2D), Dimension(:), Pointer          :: Nodes_V
    Type (Element2D_Scal), Dimension(:), Pointer  :: Elems_V
#elif defined PB_3D
    Type (Node3D), Dimension(:), Pointer          :: Nodes_V
    Type (Element3D_Scal), Dimension(:), Pointer  :: Elems_V
#else
    Type (Node2D), Dimension(:), Pointer          :: Nodes_V
    Type (Element2D_Scal), Dimension(:), Pointer  :: Elems_V
#endif
    Vec                                           :: V
    
    Integer                                      :: Rand_Node
    Real(Kind = Kr)                              :: Rand_Theta ! colatitude
    Real(Kind = Kr)                              :: Rand_Phi   ! polar angle 
    Real(Kind = Kr)                              :: Rand_Length, Length
    Integer                                      :: iTerCracks,Nb_DoF_V
    Integer                                      :: iBlk, iE, iEloc, iS, iSG ,i
    
    Real(Kind = Kr), Dimension(:), Pointer       :: V_Ptr
    Integer, Dimension(:), Pointer               :: Loc_Indices
#if defined PB_3D
    Type (Vect3D), Dimension(:),Pointer          :: Crack_Loc
#else
    Type (Vect2D), Dimension(:),Pointer          :: Crack_Loc
#endif    
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
          Rand_Length = Rand_Length * Params%MaxCrackLength

          Crack_Loc(iTerCracks*2 - 1) = Nodes_V(Rand_Node)%Coord
          Crack_Loc(iTerCracks*2 )    = Crack_Loc(iTerCracks*2 - 1) 
#if defined PB_3D
          Crack_Loc(iTerCracks*2 )%X  = Crack_Loc(iTerCracks*2)%X             &
               & + Rand_Length * COS(Rand_Theta) * SIN(Rand_Phi)
          Crack_Loc(iTerCracks*2 )%Y  = Crack_Loc(iTerCracks*2)%Y             &
               & + Rand_Length * SIN(Rand_Theta) * SIN(Rand_Phi)
          Crack_Loc(iTerCracks*2 )%Z  = Crack_Loc(iTerCracks*2)%Z             &
               & + Rand_Length * COS(Rand_PHI)
#else
          Crack_Loc(iTerCracks*2 )%X  = Crack_Loc(iTerCracks*2)%X             &
               & + Rand_Length * COS(Rand_Theta)  
          Crack_Loc(iTerCracks*2 )%Y  = Crack_Loc(iTerCracks*2)%Y             &
               & + Rand_Length * SIN(Rand_Theta)
#endif
       End Do Do_Cracks
    End If
    
#if defined PB_3D
    Call MPI_Bcast(Crack_Loc, 2*Params%nbCracks, Vect3D_MPIType, 0,           &
         & PETSC_COMM_WORLD, iErr)
#else
    Call MPI_Bcast(Crack_Loc, 2*Params%nbCracks, Vect2D_MPIType, 0,           &
         & PETSC_COMM_WORLD,iErr)
#endif      
    
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
                   Length = Distance(Crack_Loc(iTerCracks*2 -1),              &
                        & Crack_Loc(iTerCracks*2 ), Nodes_V(iSG)%Coord)
                   V_Ptr(Loc_Indices(iSG)+1) = min(V_Ptr(Loc_Indices(iSG)+1), &
                        & CrackProfile(Length,Params))
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



!!! USING VLV SEEMS TO MAKE THE KSP UNSTABLE, UNLESS V IS INITIALIZED WITH 
!!! THE PROPER BC
  Subroutine Assemb_MR_V(MR, U_Loc, Temp_Loc, Geom, Params, SD_U, SD_V,       &
       & Elems_U, Elems_V, Nodes_U, Nodes_V)
    Mat                                           :: MR
    Vec                                           :: U_Loc
    Vec                                           :: Temp_Loc
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


    Integer                                      :: Nb_Gauss
    Integer                                      :: Nb_DoF_U, Nb_DoF_V
    Integer                                      :: iSLV1
    Integer                                      :: iSLV2
    Integer                                      :: iSL, iSG
    Integer                                      :: iSLSig, iSGSig
    Integer                                      :: iSLEps, iSGEps
    Integer                                      :: iBlk, iE, iEloc, iG

    Real(Kind = Kr), Dimension(:), Pointer       :: U_Ptr
    Real(Kind = Kr), Dimension(:), Pointer       :: Temp_Ptr
    Real(Kind = Kr), Dimension(:), Pointer       :: ContrU

    PetscScalar, Dimension(:,:), Pointer         :: MR_Elem
    Integer, Dimension(:), Pointer               :: EXO_Indices_V
    Integer, Dimension(:), Pointer               :: Loc_Indices_U
    Integer, Dimension(:), Pointer               :: Loc_Indices_V

    PetscTruth                                   :: ISAssembled

    PetscLogDouble                               :: GaussTS, GaussTF, GaussT
    PetscLogDouble                               :: SetTS, SetTF, SetT

    Real(Kind = Kr)                              :: E, Nu, Lambda, Mu, k

    Integer                                      :: SetN, i

    SetN = 0
    GaussT = 0.0

    Call MatAssembled(MR, IsAssembled, iErr)
    If (IsAssembled) Then
       Call MatZeroEntries(MR, iErr)
    End If

    
#ifdef PB_2DA
    Allocate(Loc_Indices_U(Geom%Num_Nodes))
    Loc_Indices_U = (/ (i ,i = 0, Geom%Num_Nodes-1) /)
    Call AOApplicationToPETSc(SD_U%Loc_AO, Geom%Num_Nodes, Loc_Indices_U, iErr)
#else
    Allocate(Loc_Indices_U(Geom%Num_Nodes * Geom%Num_Dim))
    Loc_Indices_U = (/ (i ,i = 0, Geom%Num_Nodes * Geom%Num_Dim - 1) /)
    Call AOApplicationToPETSc(SD_U%Loc_AO, Geom%Num_Nodes * Geom%Num_Dim,     &
         & Loc_Indices_U, iErr)
#endif
    Allocate(Loc_Indices_V(Geom%Num_Nodes))
    Loc_Indices_V = (/ (i ,i = 0, Geom%Num_Nodes-1) /)
    Call AOApplicationToPETSc(SD_V%Loc_AO, Geom%Num_Nodes, Loc_Indices_V, iErr)

    Call VecGetArrayF90(U_Loc, U_Ptr, iErr) 
    Call VecGetArrayF90(Temp_Loc, Temp_Ptr, iErr) 

    Do_iBlk: Do iBlk = 1, Geom%Num_elem_blks
       E  = Params%Young_Mod(iBlk)
       nu = Params%Poisson_Ratio(iBlk) 
       k  = Params%Toughness(iBlk)

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

       Nb_Dof_V = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem

       Allocate (MR_Elem(Nb_DoF_V, Nb_DoF_V))
       Allocate (EXO_Indices_V(Nb_DoF_V))

       Do_iE: Do iELoc = 1, Geom%Elem_Blk(iBlk)%Num_Elems
          iE = Geom%Elem_Blk(iBlk)%ELem_ID(iELoc)
!!$          If ((.NOT. SD_V%IsLocal_Elem(iE)) .OR.                              &
!!$               & (.NOT. Params%Is_Domain(iBlk))) Then
!!$             CYCLE
!!$          End If
          If (.NOT. SD_V%IsLocal_Elem(iE)) Then
             CYCLE
          End If
          EXO_Indices_V = Elems_V(iE)%Id_DoF-1
          Call AOApplicationToPETSc(SD_V%EXO_AO, Nb_DoF_V, EXO_Indices_V, iErr)
       
          Call PetscGetTime(GaussTS, iErr)
          Call Init_Gauss_EXO(Elems_U, Nodes_U, Geom, MEF90_GaussOrder,       &
               & Elem=iE)
          Call Init_Gauss_EXO(Elems_V, Nodes_V, Geom, MEF90_GaussOrder,       &
               & Elem=iE)
          Call PetscGetTime(GaussTF, iErr)
          GaussT = Gausst + GaussTF - GaussTS
          MR_Elem = 0.0_Kr
          Nb_Gauss = Elems_U(iE)%Nb_Gauss
          Allocate(ContrU(Nb_Gauss))
#ifndef PB_2DA
          Allocate(Sigma(Nb_Gauss))
#endif
          !!! part related to v^2+k_\e W(e(u))
          Is_Brittle: If ( Params%Is_Brittle(iBlk)) Then
             ContrU = 0.0_Kr
             Do_iSLSig: Do iSLSig = 1, Elems_U(iE)%Nb_DoF
                iSGSig = Elems_U(iE)%ID_DoF(iSLSig)
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
                   iSGEps = Elems_U(iE)%ID_DoF(iSLEps)
                   Do_iGUEps: Do iG = 1, Nb_Gauss
#if defined PB_2DA
                      ContrU(iG) = ContrU(iG) +                               &
                           & ( Elems_U(iE)%Grad_BF(iSLEps,iG) .DotP.          &
                           &   Elems_U(iE)%Grad_BF(iSLSig,iG) ) *             &
                           &  U_Ptr(Loc_Indices_U(iSGSig)+1) *                &
                           &  U_Ptr(Loc_Indices_U(iSGEps)+1) * Mu
#else                   
                      ContrU(iG) = ContrU(iG) +                               &
                           & ( Elems_U(iE)%GradS_BF(iSLEps,iG) .DotP.         &
                           &   Sigma(iG) ) *                                  &
                           &  U_Ptr(Loc_Indices_U(iSGSig)+1) *                &
                           &  U_Ptr(Loc_Indices_U(iSGEps)+1) 
#endif
                   End Do Do_iGUEps
                End Do Do_iSLEps

#ifndef PB_2DA
                Do_iSLEps_Temp: Do iSLEps = 1, Elems_V(iE)%Nb_DoF
                   iSGEps = Elems_V(iE)%ID_DoF(iSLEps)
                   Do_iGUEps_Temp: Do iG = 1, Nb_Gauss
                      ContrU(iG) = ContrU(iG) -                               &
                		& Params%Therm_Exp(iBlk) *                            &
                		& Temp_Ptr(Loc_Indices_V(iSGEps)+1) *                 &
                		& Elems_V(iE)% BF(iSLEps, iG) *                       &
                        & U_Ptr(Loc_Indices_U(iSGSig)+1) *                    &
                        & trace(Elems_U(iE)%GradS_BF(iSLSig, iG))
                   End Do Do_iGUEps_Temp
                End Do Do_iSLEps_Temp
#endif


             End Do Do_iSLSig
             
             DoiSLV1: Do iSLV1 = 1, Elems_V(iE)%Nb_DoF
                DoiSLV2: Do iSLV2 = 1, Elems_V(iE)%Nb_DoF
                   DoiGV: Do iG = 1, Nb_Gauss
                      MR_Elem(iSLV1, iSLV2) = MR_Elem(iSLV1, iSLV2) +         &
                           & Elems_V(iE)%Gauss_C(iG) * ContrU(iG) *           &
                           & Elems_V(iE)%BF(iSLV1, iG) *                      &
                           & Elems_V(iE)%BF(iSLV2, iG) * .5_Kr
                   End Do DoiGV
                End Do DoiSLV2
             End Do DoiSLV1
          End If Is_Brittle

       !!! Surface energy part
          DoiSLV1surf: Do iSLV1 = 1, Elems_V(iE)%Nb_DoF
             DoiSLV2surf: Do iSLV2 = 1, Elems_V(iE)%Nb_DoF
                DoiGVsurf: Do iG = 1, Nb_Gauss
                   MR_Elem(iSLV1, iSLV2) = MR_Elem(iSLV1, iSLV2) +            &
                        & Elems_V(iE)%Gauss_C(iG) * k *                       &
                        & (( Elems_V(iE)%BF(iSLV1, iG) *                      &
                        &    Elems_V(iE)%BF(iSLV2, iG) ) / Params%Epsilon     &
                        &    * InvOf4                                         &
                        &  +(Elems_V(iE)%Grad_BF(iSLV1, iG) .DotP.            &
                        &    Elems_V(iE)%Grad_BF(iSLV2, iG) )                 &
                        &    * Params%Epsilon)
                End Do DoiGVsurf
             End Do DoiSLV2surf
          End Do DoiSLV1surf
          
          If (Params%Do_Irrev) Then
             Do_iSLBC: Do iSL = 1, Elems_V(iE)%Nb_DoF
                iSG = Elems_V(iE)%ID_DoF(iSL)
                If (Nodes_V(iSG)%BC /= BC_TYPE_NONE) Then
                   MR_Elem(iSL,:)    = 0.0_Kr
                   MR_Elem(:,iSL)    = 0.0_Kr
                End If
             End Do Do_iSLBC
          End If
          
          Call MatSetValues(MR, Nb_DoF_V, EXO_Indices_V, Nb_DoF_V,            &
               & EXO_Indices_V, MR_Elem, ADD_VALUES, iErr)
          SetN = SetN + 1
          
          Call PetscGetTime(GaussTS, iErr)
          Call Destroy_Gauss_EXO(Elems_U, Elem=iE)
          Call Destroy_Gauss_EXO(Elems_V, Elem=iE)
          Call PetscGetTime(GaussTF, iErr)
          GaussT = Gausst + GaussTF - GaussTS
          DeAllocate(ContrU)
#ifndef PB_2DA
          DeAllocate(Sigma)
#endif
       EndDo Do_iE
       DeAllocate(MR_Elem)
       DeAllocate(EXO_Indices_V)
    End Do Do_iBlk
    
    Call MatAssemblyBegin(MR, MAT_FLUSH_ASSEMBLY, iErr)
    DeAllocate(Loc_Indices_U)
    
    Call VecRestoreArrayF90(U_Loc, U_Ptr, iErr) 
    Call VecRestoreArrayF90(Temp_Loc, Temp_Ptr, iErr) 
    Call MatAssemblyEnd(MR, MAT_FLUSH_ASSEMBLY, iErr)  
    

    Do_Irrev: If (Params%Do_Irrev) Then
       ! Assembly of the BC terms 
       Allocate(EXO_Indices_V(Geom%Num_Nodes))
       EXO_Indices_V = (/ (i, i=0, Geom%Num_Nodes - 1) /)
       Call AOApplicationToPETSc(SD_V%EXO_AO, Geom%Num_Nodes, EXO_Indices_V,  &
            & iErr)

       Do_iBlk_BC: Do iBlk = 1, Geom%Num_Elem_Blks
          Nb_DoF_V = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem

          Do_iE_BC: Do iELoc = 1, Geom%Elem_blk(iBlk)%Num_Elems
             iE = Geom%Elem_blk(iBlk)%Elem_ID(iELoc)
             If (.NOT. SD_V%IsLocal_Elem(iE)) Then
                CYCLE
             End If
             
             Do_iSLSig_BC: Do iSLSig = 1, Nb_DoF_V
                ISGSig = Elems_V(iE)%ID_DoF(iSLSig)
                Is_BC_BC: If ( (Nodes_V(iSGSig)%BC /= BC_Type_NONE) )Then
                   Call PetscGetTime(SetTS, iErr)
                   Call MatSetValue(MR, EXO_Indices_V(iSGSig),                &
                        & EXO_Indices_V(iSGSig), 1.0_Kr, INSERT_VALUES,    &
                        & iErr)
                   Call PetscGetTime(SetTF, iErr)
                   SetN = SetN+1
                   SetT = SetT + SetTF - SetTS
                End If Is_BC_BC
             End Do Do_iSLSig_BC
          End Do Do_iE_BC
       End Do Do_iBlk_BC
       DeAllocate (EXO_Indices_V)
       
    End If Do_Irrev
    Call MatAssemblyBegin(MR, MAT_FINAL_ASSEMBLY, iErr)
    Call MatAssemblyEnd(MR, MAT_FINAL_ASSEMBLY, iErr)  
  End Subroutine Assemb_MR_V

  Subroutine Assemb_RHS_V(RHS, Geom, Params, SD, Elems, Nodes)
    Vec                                          :: RHS
    Type (EXO_Geom_Info)                         :: Geom
    Type (Rupt_Params)                           :: Params
    Type (SD_Info)                               :: SD

#ifdef PB_3D
    Type (Node3D), Dimension(:), Pointer         :: Nodes
    Type (Element3D_Scal), Dimension(:), Pointer :: Elems
#else
    Type (Node2D), Dimension(:), Pointer         :: Nodes
    Type (Element2D_Scal), Dimension(:), Pointer :: Elems
#endif

    PetscReal, Dimension(:), Pointer             :: RHS_Ptr


    Integer                                      :: Nb_Gauss, Nb_DoF
    Integer                                      :: iSLV1, iSGV1
    Integer                                      :: iSLV2
    Integer                                      :: iE, iELoc, iG
    Integer                                      :: iSLoc, iS, iBlk

    Integer, Dimension(:), Pointer               :: EXO_Indices
    Integer, Dimension(:), Pointer               :: Loc_Indices

    Integer                                      :: SetN, i

    PetscReal                                    :: Tmp_Val
    PetscReal                                    :: Toughness

    Call VecSet(0.0_Kr, RHS, iErr)

    Allocate(EXO_Indices(Geom%Num_Nodes))
    EXO_Indices = (/ (i ,i = 0, Geom%Num_Nodes - 1) /)
    Call AOApplicationToPETSc(SD%EXO_AO, Geom%Num_Nodes, EXO_Indices, iErr)

    Do_iBlk: Do iBlk = 1, Geom%Num_elem_blks
       Toughness = Params%Toughness(iBlk)

       Do_iE: Do iELoc = 1, Geom%Elem_blk(iBlk)%Num_Elems
          iE = Geom%Elem_blk(iBlk)%Elem_ID(iELoc)
          If (.NOT. SD%IsLocal_Elem(iE)) Then
             CYCLE
          End If
          Call Init_Gauss_EXO(Elems, Nodes, Geom, MEF90_GaussOrder, Elem=iE)

          Nb_Gauss = Elems(iE)%Nb_Gauss
          Do_iSLV1: Do iSLV1 = 1, Elems(iE)%Nb_DoF
             iSGV1 = Elems(iE)%ID_DoF(iSLV1)
             Tmp_Val = 0.0_Kr
             DoiSLV2: Do iSLV2 = 1, Elems(iE)%Nb_DoF
                Do_iGV: Do iG = 1, Nb_Gauss
                   Tmp_Val = Tmp_Val + Toughness * Elems(iE)%Gauss_C(iG) *    &
                        & Elems(iE)%BF(iSLV1, iG) * Elems(iE)%BF(iSLV2, iG) / &
                        & Params%Epsilon / 4.0_Kr
                End Do Do_iGV
             End Do DoiSLV2
             Call VecSetValue(RHS, EXO_Indices(iSGV1), Tmp_Val, ADD_VALUES,   &
                  & iErr)
          End Do Do_iSLV1
          Call Destroy_Gauss_EXO(Elems, Elem=iE)
       End Do Do_iE
    End Do Do_iBlk
    Call VecAssemblyBegin(RHS, iErr)
    DeAllocate(EXO_Indices)
    Call VecAssemblyEnd(RHS, iErr)

    ! BC Due to irreversibility:
    Is_Irrev: If (Params%Do_Irrev) Then
       Allocate(Loc_Indices(Geom%Num_Nodes))
       Loc_Indices = (/ (i ,i = 0, Geom%Num_Nodes - 1) /)
       Call AOApplicationToPETSc(SD%Loc_AO, Geom%Num_Nodes, Loc_Indices, iErr)

       Call VecGetArrayF90(RHS, RHS_Ptr, iErr)
       Do_iS: Do iSLoc = 1, SD%Num_Nodes
          iS = SD%Node(iSLoc)
          Is_BC: If (Nodes(iS)%BC /= BC_Type_NONE) Then
             RHS_Ptr(Loc_Indices(iS)+1) = 0.0_Kr
!             RHS_Ptr(Loc_Indices(iS)+1) = -1.0_Kr * sqrt(MEF90_VLV)
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
#ifdef PB_3D
    Type(Node3D), Dimension(:), Pointer              :: Nodes
#else
    Type(Node2D), Dimension(:), Pointer              :: Nodes
#endif
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


#if defined PB_2D
End Module m_Rupt2D_V
#elif defined PB_3D
End Module m_Rupt3D_V
#else
End Module m_Rupt2DA_V
#endif
