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
      Integer, Dimension(:), Pointer                :: MyIndices
      Integer                                       :: i, iS, iE, iBlk
                                                    
      PetscLogDouble                                :: BCastTS, BCastTF
      PetscLogDouble                                :: InitVectTS, InitVectTF
      PetscTruth                                    :: Has_Sim_Str
      Integer, Dimension(:,:), Pointer              :: MyConnect
      Type(Vect2D), Dimension(:), Pointer           :: MyCoord
      
      
      Call MEF90_Initialize()
      MEF90_GaussOrder = 2 
      
      
      Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-f', Params%Sim_Str, Has_Sim_Str, iErr)
      
      If (.NOT. Has_Sim_Str) Then
         Write(CharBuffer, 100) 'Simulation name: \n'c
         Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
         If (MEF90_MyRank == 0) Then
            Read (*,100) Params%Sim_Str  
         End If
         Call MPI_BCAST(Params%Sim_Str, MXLNLN, MPI_CHARACTER, 0, MPI_COMM_WORLD, iErr)
      End If
      
      Geom%Filename    = Trim(Params%Sim_Str) // '.gen'
      Params%PARAM_Str = Trim(Params%Sim_Str) // '.PARAM'
      Params%CST_Str   = Trim(Params%Sim_Str) // '.CST'
      Geom%Comm        = PETSC_COMM_WORLD
      Geom%Numbering   = Numbering_PerNodes
      
      Write(CharBuffer, 99) MEF90_MyRank
 99   Format('-', I4.4, '.gen')
      MyGeom%FileName = Trim(Params%Sim_Str) // Trim(CharBuffer)
!      MyGeom%FileName = Trim(Params%Sim_Str) // Trim(adjustL(CharBuffer)) // '.gen'
      
      Call Read_EXO_Geom_Info(Geom, MyGeom)
      MyGeom%Numbering = Numbering_PerNodes
      Call Init_Layout_TRI3(Geom, MyGeom, Layout, Node_db, Elem_db)
      
!      Call MEF90_Finalize()
!      STOP
      Call Read_Rupt_EXO_Params(Geom, Params)
      Call Read_Rupt_DATA(MyGeom, Params)


      Call Write_EXO_Geom_Info(MyGeom)
      Call Write_EXO_Node_Coord(MyGeom, Node_db, 1)
      Call Write_EXO_Connect(MyGeom, Elem_db)
      Call Write_Rupt_EXO_Params(MyGeom, Params)
      
!      Call MEF90_Finalize()
!      STOP

      Write(CharBuffer,*) 'Total number of dof:        ', Geom%Num_Nodes, '\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
      Write(CharBuffer,*) 'Local number of dof/ghosts: ', Layout%Num_Local_dof, Layout%Num_Ghost_dof, '\n'c
      Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
      Call PetscSynchronizedFlush(PETSC_COMM_WORLD, iErr)
      
      Write(CharBuffer,*) 'Total number of elements:   ', Geom%Num_Elems, '\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
      Write(CharBuffer,*) 'Local number of elemens:    ', Layout%Num_Local_elems, '\n'c
      Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
      Call PetscSynchronizedFlush(PETSC_COMM_WORLD, iErr)

      Call Init_BC(MyGeom, Params, Node_db)
      
      Call VecCreateMPI(Geom%Comm, Layout%num_local_dof, Geom%num_nodes, Sol, iErr)
      Call VecDuplicate(Sol, RHS, iErr)
      Call VecSetLocalToGlobalMapping(RHS, Layout%Mapping_N, iErr)
      
      Call VecCreateGhost(Geom%Comm, Layout%num_local_dof, Geom%num_nodes, Layout%num_ghost_dof, Layout%ghost_dof-1, BC, iErr)
      Call VecCreateGhost(Geom%Comm, Layout%num_local_dof, Geom%num_nodes, Layout%num_ghost_dof, Layout%ghost_dof-1, F, iErr)
      
      Call VecGhostGetLocalForm(F, F_Local, iErr)
      Call VecGhostGetLocalForm(BC, BC_Local, iErr)
      
      If (MEF90_MyRank == 0) Then
         Call VecCreateMPI(Geom%Comm, Geom%num_nodes, Geom%num_nodes, IO_N, iErr)
      Else
         Call VecCreateMPI(Geom%Comm, 0, Geom%num_nodes, IO_N, iErr)
      End If                  

      Call MatCreateMPIAIJ(Geom%Comm, Layout%num_local_dof, Layout%num_local_dof, Geom%Num_Nodes, Geom%Num_Nodes, 24, PETSC_NULL_INTEGER, 24, PETSC_NULL_INTEGER, MR, iErr)
      Call MatSetLocalToGlobalMapping(MR, Layout%Mapping_N, iErr)

      Call MatSetOption(MR, MAT_SYMMETRIC, iErr)
      Call MatSetFromOptions(MR, iErr)

 100  Format(A)
 101  Format(I3.3)
 200  Format('Rank: ', I4,' Index Range: ', I7,' ',I7, '\n'c) 
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
      
      Call KSPSetTolerances(KSP_MR, Params%TolKSP, PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_INTEGER, iErr)
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
      Call Write_EXO_Result_Nodes(Geom, Layout, 4, TimeStep, Sol)
      !!! TEST
      Call VecCopy(Sol, F, iErr)
      Call VecGhostUpdateBegin(F, INSERT_VALUES, SCATTER_FORWARD, iErr)
      Call VecGhostUpdateEnd  (F, INSERT_VALUES, SCATTER_FORWARD, iErr)

      Call Write_EXO_Result_Nodes(MyGeom, Layout, 4, TimeStep, F_Local)
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
      Type(Rupt_Params2D), Intent(IN)                  :: Params
      
      Type(Node2D), Dimension(:), Pointer              :: Node_db
      Integer                                          :: iN, iSet
      
      If (Geom%Numbering /= Numbering_PerNodes) Then
         Print*, '[ERROR] Init_BC not implemented for this numbering scheme'
         STOP
      End If
      
      Node_db(:)%BC = BC_Type_NONE
      Do iSet = 1, Geom%Num_Node_Sets
         Do iN = 1, Geom%Node_Set(iSet)%Num_Nodes
            Node_db(Geom%Node_Set(iSet)%Node_ID(iN))%BC = Params%BC_Type_Z(iSet)
         End Do
      End Do
   End Subroutine Init_BC


   Subroutine Assemb_Mat_Elast2DA(MR, Geom, Params, Elem_db, Node_db)
      Mat                                                 :: MR
      Type (EXO_Geom_Info)                                :: Geom
      Type (Rupt_Params2D)                                :: Params
      Type (Node2D), Dimension(:), Pointer                :: Node_db 
      Type (Element2D_Scal), Dimension(:), Pointer        :: Elem_db 
      
      Integer                                             :: Nb_DoF
      Integer                                             :: iSL1, iSL2
      Integer                                             :: iSG1, iSG2
      Type (Vect2D)                                       :: Sigma, Epsilon
      Integer                                             :: iE, iG, iELoc
      Integer                                             :: iBlk
!      Integer                                             :: i
      
      PetscScalar, Dimension(:,:), Pointer                :: MR_Elem
      PetscTruth                                          :: IsAssembled
      Integer, Dimension(:), Pointer                      :: DoF
      
      Call MatAssembled(MR, IsAssembled, iErr)
      If (IsAssembled) Then
         Call MatZeroEntries(MR, iErr)
      End If
    
    ! Assembly of the Non BC terms 
      Do_iBlk: Do iBlk = 1, Geom%Num_Elem_Blks
         Nb_DoF = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem
         Allocate (MR_Elem(Nb_DoF, Nb_DoF))
         Allocate (DoF(Nb_DoF))
         
         Do_iE: Do iELoc = 1, Geom%Elem_Blk(iBlk)%Num_Elems
            iE = Geom%Elem_Blk(iBlk)%Elem_ID(iELoc)
            DoF = Elem_db(iE)%ID_DoF-1
            Call Init_Gauss_EXO(Elem_db, Node_db, Geom, GaussOrder, Elem=iE)
            Sigma = 0.0_Kr

            MR_Elem = 0.0_Kr
            Do_iG: Do iG = 1, Elem_db(iE)%Nb_Gauss
               Sigma = 0.0_Kr
               Do_iSL1: Do iSL1 = 1, Elem_db(iE)%Nb_DoF
                  iSG1 = Elem_db(iE)%ID_DoF(iSL1)
                  Sigma = 2.0_Kr * Params%Hookes_Law(iBlk)%XYXY * Elem_db(iE)%Grad_BF(iSL1,iG)
                  Do_iSL2: Do iSL2 = 1, Elem_db(iE)%Nb_DoF
                     iSG2 = Elem_db(iE)%ID_DoF(iSL2)
                     Epsilon = Elem_db(iE)%Grad_BF(iSL2,iG)
                     MR_Elem(iSL2, iSL1) = MR_Elem(iSL2, iSL1) + Elem_db(iE)%Gauss_C(iG)  * (Sigma .DotP. Epsilon)
                  End Do Do_iSL2
               End Do Do_iSL1
            End Do Do_iG
            Call MatSetValuesLocal(MR, Nb_DoF, DoF, Nb_DoF, DoF, MR_Elem, ADD_VALUES, iErr)
            Call Destroy_Gauss_EXO(Elem_db, Elem=iE)
         EndDo Do_iE
         DeAllocate (MR_Elem)
         DeAllocate(DoF)
      End Do Do_iBlk

      Call MatAssemblyBegin(MR, MAT_FLUSH_ASSEMBLY, iErr)  
      Call MatAssemblyEnd(MR, MAT_FLUSH_ASSEMBLY, iErr)  
    
      ! Assembly of the BC terms 
      Do_iBlk_BC: Do iBlk = 1, Geom%Num_Elem_Blks
         Nb_DoF = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem

         Do_iE_BC: Do iELoc = 1, Geom%Elem_blk(iBlk)%Num_Elems
         iE = Geom%Elem_blk(iBlk)%Elem_ID(iELoc)
         Do_iSLSig_BC: Do iSL1 = 1, Nb_DoF
            ISG1 = Elem_db(iE)%ID_DoF(iSL1)
               Is_BC_BC: If ( Node_db(iSG1)%BC /= BC_Type_NONE ) Then
                  Call MatSetValueLocal(MR, iSG1 - 1, iSG1 - 1, VLV, INSERT_VALUES, iErr)
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
      Type (Rupt_Params2D)                                :: Params
      Type(Element2D_Scal), Dimension(:), Pointer         :: Elem_db
      Type(Node2D), Dimension(:), Pointer                 :: Node_db
      Vec                                                 :: BC, F
      
      Real(Kind = Kr), Dimension(:), Pointer              :: BC_Ptr
      Real(Kind = Kr), Dimension(:), Pointer              :: RHS_Ptr
      Real(Kind = Kr), Dimension(:), Pointer              :: F_Ptr
      Real(Kind = Kr), Dimension(:), Pointer              :: RHS_Elem
      Real(Kind = Kr)                                     :: F_Elem
      
      
      Integer                                             :: Nb_DoF
      Integer                                             :: iSL1, iSG1
      Integer                                             :: iSL2, iSG2
      Integer                                             :: iE, iELoc, iG
      Integer                                             :: iBlk
      
      Integer                                             :: i, iS
      
      PetscLogDouble                                      :: TotTS, TotTF, TotT
      
      Real(Kind = Kr)                                     :: Tmp_Val


      Call VecSet(RHS, 0.0_Kr, iErr)
      Call VecGetArrayF90(F, F_Ptr, iErr)
      
      Do_iBlk: Do iBlk = 1, Geom%Num_elem_blks
         Nb_DoF = Geom%Elem_blk(iBlk)%Num_Nodes_per_elem
         Allocate(RHS_Elem(Nb_DoF))
         Do_iE: Do iELoc = 1, Geom%Elem_blk(iBlk)%Num_Elems
            iE = Geom%Elem_blk(iBlk)%Elem_ID(iELoc)
            Call Init_Gauss_EXO(Elem_db, Node_db, Geom, GaussOrder, Elem=iE)
            RHS_Elem = 0.0_Kr
            F_Elem   = 0.0_Kr
            
            Do_iG: Do iG = 1, Elem_db(iE)%Nb_Gauss
               Do_iSL1: Do iSL1 = 1, Nb_DoF
                  iSG1 = Elem_db(iE)%ID_DoF(iSL1)
                  F_Elem = F_Elem + F_Ptr(iSG1) * Elem_db(iE)%BF(iSL1, iG)
               End Do Do_iSL1
               Do_iSL2: Do iSL2 = 1, Nb_DoF
                  iSG2 = Elem_db(iE)%ID_DoF(iSL2)
                  RHS_Elem(iSL2) = RHS_Elem(iSL2) + F_Elem * Elem_db(iE)%Gauss_C(iG) * Elem_db(iE)%BF(iSL2, iG)
               End Do Do_iSL2
            End Do Do_iG
            Call VecSetValuesLocal(RHS, Elem_db(iE)%Nb_DoF, Elem_db(iE)%ID_DoF-1, RHS_Elem, ADD_VALUES, iErr)   
            Call Destroy_Gauss_EXO(Elem_db, Elem=iE)
         End Do Do_iE
         DeAllocate(RHS_Elem)
      End Do Do_iBlk
      Call VecAssemblyBegin(RHS, iErr)
      Call VecAssemblyEnd(RHS, iErr)
      Call VecRestoreArrayF90(F, F_Ptr, iErr)
      
      Call VecGetArrayF90(BC, BC_Ptr, iErr)
      Call VecGetArrayF90(RHS, RHS_Ptr, iErr)

      Do_iS: Do iS = 1, Geom%Num_Nodes !Layout%num_local_dof
         Is_BC: If ( Node_db(iS)%BC /= BC_Type_NONE ) Then
            RHS_Ptr(iS) = BC_Ptr(iS) * VLV
         End If Is_BC
      End Do Do_iS

      Call VecRestoreArrayF90(BC, BC_Ptr, iErr)
      Call VecRestoreArrayF90(RHS, RHS_Ptr, iErr)
  End Subroutine Assemb_RHS_Elast2DA

End Module m_Elast2DA_Proc
