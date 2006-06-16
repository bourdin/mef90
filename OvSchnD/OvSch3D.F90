Program Poisson

  Use m_MEF90
  Use m_Poisson_Struct

#ifdef PB_2D
  Use m_OvSch2D_Vars
  Use m_OvSch2D_Procs
#else
  Use m_OvSch3D_Vars
  Use m_OvSch3D_Procs
#endif

  Implicit NONE
#include "include/finclude/petsc.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petscao.h"

  Integer                                       :: NbIter
  PetscTruth                                    :: flg

Integer iBlk, iSL, iE, Nb_DoF
Integer, Dimension(:), Pointer     :: iSG
!! Initalization 
  Call Init_Broadcast()


  Allocate(stages(11))
  If (MyRank ==0) Then
    Write(CharBuffer, *) '====== Method Information  ========\n'c
    Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
    Call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-method',Method,flg,ierr)
    Write(CharBuffer, *) 'Method:   ', Method, '\n'c
    Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
    If ((Method == 'OBDD') .or. (Method =='DIRI') .or. (Method =='CG')) Then
      If (Method == 'OBDD') Then
        Num_Method = 2
      ElseIf (Method == 'DIRI') Then
        Num_Method = 1
      Else 
        Num_Method = 0
      End If
    Else
      Num_Method = -1
    End If
    Write(CharBuffer, *) '  \n'c
    Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
  End If
  Call MPI_BCAST(Num_Method, 1, MPI_INTEGER, 0,                &
         & PETSC_COMM_WORLD, iErr)
  Call Init_Node_Ownership()
  Call Init_PU()

!!! JH: Can you be more explicit about what the stages mean?
!!!     Do we need all 11 stages?

  Call PetscLogStageRegister(stages(1),"stage 1", iErr);
  Call PetscLogStageRegister(stages(2),"stage 2", iErr);
  Call PetscLogStageRegister(stages(3),"stage 3", iErr);
  Call PetscLogStageRegister(stages(4),"stage 4", iErr);
  Call PetscLogStageRegister(stages(5),"stage 5", iErr);
  Call PetscLogStageRegister(stages(6),"stage 6", iErr);
  Call PetscLogStageRegister(stages(7),"stage 7", iErr);
  Call PetscLogStageRegister(stages(8),"stage 8", iErr);
  Call PetscLogStageRegister(stages(9),"stage 9", iErr);  
  Call PetscLogStageRegister(stages(10),"stage 10", iErr);
  Call PetscLogStageRegister(stages(11),"stage 11", iErr);


  If (MyRank == 0) Then
     Write(CharBuffer, *) '====== PU Generation ========\n'c
     Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
  End If

!! PU Generation
  Call PetscGetTime(CoTS, iErr)

!****************** ST 1 ********
  Call PetscLogStagePush(stages(1), iErr)
  Call G_PU_2D_Scal_Dist(Geom, Elem_db, Node_db, ovlp, MyPU, MyElem_Ovlp, &
     &  MyNode_Ovlp, MyElem_Ghost, MyNode_Ghost, MyNode_Bd)
  Call PetscLogStagePop(iErr)
!*******************************

  Call PetscGetTime(CoTF, iErr)

  If (MyRank ==0) Then
    Write(CharBuffer,*) 'Total time in PUGenerating:                 ',      &
          & CoTF- CoTS, '\n'c
    Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
  End If
  If (MyRank == 0) Then
    Write(CharBuffer, *) '====== Initial Size and Save  ========\n'c
    Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
  End If
  Call PetscGetTime(CoTS, iErr)

!! Other Initializations

!****************** ST 2 ********  
  Call PetscLogStagePush(stages(2), iErr)
  Call Init_Size()
  Call PetscLogStagePop(iErr)
!*******************************

  Call SaveLayout()
  Call PetscGetTime(CoTF, iErr)
  If (MyRank ==0) Then
    Write(CharBuffer,*) 'Total time in Initial Size and Save:               ',&
          & CoTF- CoTS, '\n'c
    Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
  End If

  If (MyRank == 0) Then
    Write(CharBuffer, *) '========== MATRIX ==========\n'c
    Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
  End If

!! Accemble Mass and Local Matrices
 Call Assemb_Mat_Poisson(MR, Geom, Params, Elem_db, Node_db)



  If (MyRank == 0) Then
    Write(CharBuffer, *) '====== Local MATRIX ========\n'c
    Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
  End If

  !****************** ST 3 ********
  Call PetscLogStagePush(stages(3), iErr)
  Call Init_Size_Local()
  Call PetscLogStagePop(iErr)
  !*******************************
 
  !****************** ST 4 ********
   Call PetscLogStagePush(stages(4), iErr)
   Call Assemb_Mat_Local_Poisson(MyMR, Geom, Params, Elem_db, Node_db)
   Call PetscLogStagePop(iErr)
  !*******************************

  If (MyRank == 0) Then
    Write(CharBuffer, *) '====== Local Neumann MATRIX ========\n'c
    Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
  End If
  
  !****************** ST 5 ********
  Call PetscLogStagePush(stages(5), iErr)
  Call Assemb_Mat_Local_Poisson_Neu(MyMR_Neu, Geom, Params, Elem_db, Node_db)
  Call PetscLogStagePop(iErr)
  !*******************************

  If (MyRank == 0) Then
    Write(CharBuffer, *) '====== Coarse Matrix ========\n'c
    Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
  End If
  Call PetscGetTime(CoTS, iErr)

  !****************** ST 6 ********
  Call PetscLogStagePush(stages(6), iErr)
  Call Communicate_PU()
   Call PetscLogStagePop(iErr)
  !*******************************
    
   Call PetscGetTime(CoTF, iErr)
   If (MyRank ==0) Then
      Write(CharBuffer,*) 'Total time in Communicate_PU:                 ',    &
           & CoTF- CoTS, '\n'c
      Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
  End If

  Call PetscGetTime(CoTS, iErr)

  !****************** ST 7 ********
  Call PetscLogStagePush(stages(7), iErr)
  Call Assemb_Coarse_Seq()
  Call PetscLogStagePop(iErr)
  !*******************************

  Call PetscGetTime(CoTF, iErr)
  If (MyRank ==0) Then
    Write(CharBuffer,*) 'Total time in Coarse Seq Matrix:                 ', &
          & CoTF- CoTS, '\n'c
    Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
  End If

!! RHS 
  If (MyRank == 0) Then
    Write(CharBuffer, *) '========== VECTOR ==========\n'c
    Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
  End If
  Call Assemb_RHS_Poisson(RHS, Geom, Params, Elem_db, Node_db, Load)
  Call VecGetArrayF90(RHS, RHS_Ptr, iErr)
  Call Assemb_RHS_Poisson_Local(RHS_Ptr)
  Call VecRestoreArrayF90(RHS, RHS_Ptr, iErr)
  initial = 0.0
  Call VecSet(SOL_Ov, initial, iErr)
  Call VecGhostGetLocalForm(SOL_Ov,MySOL_Ov,ierr)

!! KSP Initialization
  If (Num_Method /= 0 ) Then
    Call Init_KSP()
  End If
  If (Num_Method == 0 ) Then
    Call Init_KSP_OLD()
  End If
  Call VecGhostRestoreLocalForm(SOL_Ov,MySOL_Ov,ierr) 
  If (MyRank == 0) Then
    Write(CharBuffer, *) '========== Coarse SOLVER ==========\n'c
    Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
  End If

!! Prepare Coarse Problem and Solve Routine Before CG
  Call PetscGetTime(CoSoTS, iErr)

  
  !****************** ST 8 ********
  Call PetscLogStagePush(stages(8), iErr)
  Call Coarse_Solve_Prep()
  Call PetscLogStagePop(iErr)
  !*******************************

  Call PetscGetTime(CoSoTF, iErr)
  If (MyRank ==0) Then
    Write(CharBuffer,*) 'Total time in Coarse Prep:                 ',  &
          & CoSoTF - CoSoTS, '\n'c
    Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
  End If
!!  Call PetscGetTime(CoSoTS, iErr)
   Call VecGhostGetLocalForm(RHS_Ov,MyRHS_Ov,ierr) 
  If (Num_Method /= 0 ) Then
!!  Call VecGhostGetLocalForm(RHS_Ov,MyRHS_Ov,ierr) 
   Call PetscGetTime(CoSoTS, iErr)
  Call Coarse_Solve_Seq(RHS,SOL_Ov)
   Call PetscGetTime(CoSoTF, iErr)
!!  Call VecGhostRestoreLocalForm(CoPU_Ov,MyCoPU_Ov,ierr)
  End If
    Call VecGhostRestoreLocalForm(RHS_Ov,MyRHS_Ov,ierr)
!!  Call PetscGetTime(CoSoTF, iErr)
  If (MyRank ==0) Then
    Write(CharBuffer,*) 'Total time in Coarse Solve:                 ',      &
          & CoSoTF - CoSoTS, '\n'c
    Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
  End If

!! Imformation for KSP and Coarse KSP
  If (MyRank == 0) Then
     Write(CharBuffer, *) '========== KSP SOLVER ==========\n'c
     Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
  End If
  Call KSPView(KSP_MR,PETSC_VIEWER_STDOUT_WORLD, iErr)
  If (MyRank == 0) Then
  Write(CharBuffer, *) '========== Local SOLVER ==========\n'c
  Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
  If (Num_Method == 1 ) Then
  Call KSPView(KSP_MyMR,PETSC_VIEWER_STDOUT_SELF, iErr)
  End If
  If (Num_Method == 2 ) Then
   Call KSPView(KSP_MyMR_Neu,PETSC_VIEWER_STDOUT_SELF, iErr)
  End If
  End If
  
!! Solver Using CG and PC
  If (MyRank == 0) Then
     Write(CharBuffer, *) '========== SOLVER  ==========\n'c
     Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
  End If
!  zz=0.0_Kr
!  vv=0.0_Kr
  none =-1.0
  pone = 1.0
  Call VecGhostGetLocalForm(RES_Ov,MyRES_Ov,ierr)
  Call VecGhostGetLocalForm(IMP_Ov,MyIMP_Ov,ierr)
  Call VecGhostGetLocalForm(Update_Ov,MyUpdate_Ov,ierr)
  Call PetscGetTime(SolveTS, iErr)

  !****************** ST 9 ********
  Call PetscLogStagePush(stages(9), iErr)
  Call KSPSolve(KSP_MR, RHS, SOL_Ov,iErr)
  Call PetscLogStagePop(iErr)
  !*******************************

  Call PetscGetTime(SolveTF, iErr)
  If (MyRank == 0) Then
     Call VecCreateMPI(PETSC_COMM_WORLD, Geom%Num_Nodes,     &
          & Geom%Num_Nodes, SOL_Ov_Master, iErr)
  Else
     Call VecCreateMPI(PETSC_COMM_WORLD, 0, Geom%Num_Nodes,  &
          & SOL_Ov_Master, iErr)
  End If
  Call VecScatterBegin(SOL_Ov, SOL_Ov_Master, INSERT_VALUES, SCATTER_FORWARD,     &
       & SOL_DistToMaster, iErr)
  Call VecScatterEnd(SOL_Ov, SOL_Ov_Master, INSERT_VALUES, SCATTER_FORWARD,       &
       & SOL_DistToMaster, iErr)
  Call KSPGetIterationNumber(KSP_MR, NbIter, iErr)
  Call KSPGetResidualNorm(KSP_MR, FinalRe, iErr)
  Call KSPGetConvergedReason(KSP_MR, Reason, iErr)
  If (MyRank ==0) Then
    Write(CharBuffer,*) 'Number of iterations:            ',      &
          & NbIter, '\n'c
    Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
    Write(CharBuffer,*) 'Norm of Residual:            ',      &
          & FinalRe, '\n'c
    Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
    Write(CharBuffer,*) 'Total time in KSP_Solve:                 ',      &
          & SolveTF- SolveTS, '\n'c
    Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
  End If
  If (MyRank== 0) Then
    Print*, 'Rank  : ' ,MyRank, 'ASM SetUp Time     :', ASMSETUPT 
    Print*, 'Rank  : ' ,MyRank, 'ToTal ASM Time     :', ASMT
    Print*, 'Rank  : ' ,MyRank, 'ASM Per iter       :', ASMT/NbIter
  End If

  Allocate (Disp_Ov_L2(Geom%Num_Nodes))
  If (MyRank == 0) Then
    Call VecGetArrayF90(SOL_Ov_Master, Disp_Sol_Ov, iErr)
    Call Write_EXO_Result_Ptr_Nodes(Geom, 2, 1, Disp_Sol_Ov)
    Write(*,90) MinVal(Disp_Sol_Ov), MaxVal(Disp_Sol_Ov)
    Disp_Ov_L2 = Disp_Sol_Ov
    Call VecRestoreArrayF90(SOL_Ov_Master, Disp_Sol_Ov, iErr)
  End If
  Call MPI_BCAST(Disp_Ov_L2, Geom%Num_Nodes, MPI_DOUBLE_PRECISION, 0, &
   & PETSC_COMM_WORLD, iErr)
  Allocate (D_Vec(Geom%Num_Nodes))
  D_vec = Load - Disp_Ov_L2
  Call CompL2(D_vec, Geom, Elem_db, Node_db, MyL2)
  Call MPI_REDUCE(MyL2,FinalL2,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,PETSC_COMM_WORLD, iErr)
  Call CompSEMIH2(D_vec, Geom, Elem_db, Node_db, MySEMIH2)
  Call MPI_REDUCE(MySEMIH2,FinalSEMIH2,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
   &PETSC_COMM_WORLD, iErr)

!! No Meaning but useful for Comparison
  If (MyRank ==0) Then
    Print*, 'L2 :   ', sqrt(FinalL2)
    Print*, 'SemiH1  : ' ,sqrt(FinalSEMIH2)
    Print*, 'H1 :   ', sqrt(FinalL2)+ sqrt(FinalSEMIH2)
    MyLIf = maxval(ABS(D_vec))
    Print*, 'Infty norm  : ', MyLif
  End If
  Call PetscLogPrintSummary(PETSC_COMM_WORLD, "petsc_log_summary.log", iErr)
  Call Finalize()
90 Format('Solution Min / max: ', 2(ES10.3,' '))
100 Format(A)

End Program Poisson
