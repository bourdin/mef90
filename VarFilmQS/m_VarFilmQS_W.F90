Module m_VarFilmQS_W
#include "finclude/petscdef.h"

Use m_VarFilmQS_Types
Use m_VarFilmQS_Post
Use m_MEF90
Use m_Film_Struct

Implicit NONE
Private 
   


Public :: Init_TS_W
Public :: Update_IrrevW
Public :: FW_Assembly
Public :: Step_W
Public :: W_Solve

Contains

Subroutine Init_TS_W(AppCtx)
   !!! Set the initial value of W at the beginning of each alternate minimizatins iterations
   Type(AppCtx_Type)                            :: AppCtx
   PetscInt                                     :: iErr
   Character(len=MEF90_MXSTRLEN)                :: IOBuffer   
   PetscInt                                     :: i
   PetscReal, Dimension(:,:), Pointer           :: CoordArray
   PetscReal, Dimension(:), Pointer             :: Coordlocal
   
     
   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) "Initializing W with ", AppCtx%VarFracSchemeParam%InitW, "\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If    
   Select Case(AppCtx%VarFracSchemeParam%InitW)
      Case(VarFrac_INIT_W_PREV) ! -initw 0 (default)
         Call FieldInsertVertexBoundaryValues(AppCtx%W, AppCtx%WBC, AppCtx%BCWFlag, AppCtx%MeshTopology)
         Call FieldInsertVertexBoundaryValues(AppCtx%W, AppCtx%WIrrev, AppCtx%WIrrevFlag, AppCtx%MeshTopology)
      Case(VarFrac_INIT_W_FILE) ! -initw 6
   If ( AppCtx%TimeStep==1 ) Then
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Delamination)%Offset, AppCtx%TimeStep, AppCtx%W)
   Else
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Delamination)%Offset, AppCtx%TimeStep-1, AppCtx%W)
   end if
   Call SectionRealToVec(AppCtx%W%Sec, AppCtx%W%Scatter, SCATTER_REVERSE, AppCtx%W%Vec, ierr); CHKERRQ(ierr)
      Case(VarFrac_INIT_W_ONE) ! -initw 4
         Call SectionRealSet(AppCtx%W%Sec, 1.0_Kr, iErr); CHKERRQ(iErr)      
         Call SectionRealToVec(AppCtx%W%Sec, AppCtx%W%Scatter, SCATTER_FORWARD, AppCtx%W%Vec, ierr); CHKERRQ(ierr)
         Call VecSet(AppCtx%W%Vec, 1.0_Kr, iErr); CHKERRQ(iErr)      
         
      Case Default   
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP, 'Not Implemented yet\n', iErr)
   End Select
     
End Subroutine Init_TS_W
   
Subroutine Update_IrrevW(AppCtx)
!!! Updates the WBCFlag field to account for irreversibilty 
Type(AppCtx_Type)                            :: AppCtx

PetscInt                                     :: i, iErr
Character(len=MEF90_MXSTRLEN)                :: IOBuffer      
PetscInt                                     :: MyIrrevEQ_Counter 
PetscInt                                     :: IrrevEQ_Counter
PetscReal, Dimension(:), Pointer             :: WIrrev_Ptr
PetscInt, Dimension(:), Pointer              :: WIrrevFlag_Ptr

MyIrrevEq_Counter = 0
IrrevEq_Counter   = 0

Select Case(AppCtx%VarFracSchemeParam%IrrevType)
Case Default
   If (AppCtx%AppParam%verbose > 1) Then
      Write(IOBuffer, *) "Irreversibility for W\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If

   Allocate(WIrrevFlag_Ptr(1))
   WIrrevFlag_Ptr = VarFrac_BC_Type_DIRI

   Do i = 1, AppCtx%MeshTopology%Num_Verts
      Call SectionRealRestrict(AppCtx%W%Sec, AppCtx%MeshTopology%Num_Elems + i-1, WIrrev_Ptr, iErr); CHKERRQ(ierr)      
      If (WIrrev_Ptr(1) < AppCtx%VarFracSchemeParam%IrrevTol) Then
!        If (AppCtx%AppParam%verbose > 1) Then
!           Write(IOBuffer, *) "WIrrev_Ptr(1)", WIrrev_Ptr(1), "\n"
!           Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
!        End If
         WIrrev_Ptr(1) = 0.0_Kr
         Call SectionIntUpdate(AppCtx%WIrrevFlag%Sec, AppCtx%MeshTopology%Num_Elems + i-1, WIrrevFlag_Ptr, INSERT_VALUES, iErr); CHKERRQ(iErr)
         Call SectionRealUpdate(AppCtx%WIrrev%Sec, AppCtx%MeshTopology%Num_Elems + i-1, WIrrev_Ptr, INSERT_VALUES, iErr); CHKERRQ(iErr)
         MyIrrevEQ_Counter = MyIrrevEQ_Counter + 1
      End If
      Call SectionRealRestore(AppCtx%W%Sec, AppCtx%MeshTopology%Num_Elems + i-1, WIrrev_Ptr, iErr); 
   End Do
   Call MPI_AllReduce(MyIrrevEQ_Counter, IrrevEQ_Counter, 1, MPIU_INTEGER, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) "Dirichlet nodes Wirrev", MyIrrevEQ_Counter, "\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If

   DeAllocate(WIrrevFlag_Ptr)

End Select
End Subroutine Update_IrrevW

!!!
!!! Global assembly functions
!!!


Subroutine FW_Assembly(AppCtx)
   !!! Global dispatch routine for F of the W-problem
      Type(AppCtx_Type)                            :: AppCtx

      PetscInt                                     :: iErr
      PetscInt                                     :: iBlk
      Type(Vec)                                    :: FW_Vec
   
       Call PetscLogStagePush(AppCtx%LogInfo%FWAssemblyW_Stage, iErr); CHKERRQ(iErr)

       Call SectionRealZero(AppCtx%FW%Sec, iErr); CHKERRQ(iErr)

      Do_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Call FW_AssemblyBlk(iBlk, AppCtx)
      End Do Do_iBlk

   Call SectionRealComplete(AppCtx%FW%Sec, iErr); CHKERRQ(iErr)

   If (AppCtx%AppParam%verbose > 2) Then
      Call VecView(AppCtx%FW%Vec, AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
   End If

   Call PetscLogStagePop(iErr); CHKERRQ(iErr)
End Subroutine FW_Assembly
   
Subroutine FW_AssemblyBlk(iBlk, AppCtx)
   PetscInt             :: iBlk
   Type(AppCtx_Type)          :: AppCtx

      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_ELem = \sum_i U_Loc(i) BF(i))
   PetscInt             :: iE, iEloc, iErr, iBlk_glob
   Type(Vect2D)               :: U_eff_elem
   PetscReal, Dimension(:), Pointer    :: F_loc
   PetscReal, Dimension(:), Pointer    :: U_loc, U0_loc


   PetscInt                                     :: NumDoFScal, NumDoFVect
   PetscInt                                     :: iDoF, iGauss
   PetscLogDouble, Parameter                    :: oneflop = 1.0
   
   NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
   NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim

   iBlk_glob=AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
   If (AppCtx%MyEXO%EBProperty(VarFrac_EBProp_IsDebondable)%Value(iBlk_glob) == 0) Then
      Call SectionRealZero(AppCtx%WBC%Sec, iErr); CHKERRQ(iErr)         ! When the block is non debondable, (min wrt W is trivial)
   End If
   Allocate(F_loc(NumDoFScal))
   Allocate(U_loc(NumDoFVect))
   Allocate(U0_loc(NumDoFVect))
   
   Do_iEloc: Do iEloc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
      iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
      
      Call SectionRealRestrictClosure(AppCtx%U%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U_loc, iErr); CHKERRQ(ierr)
      Call SectionRealRestrictClosure(AppCtx%U0%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U0_loc, iErr); CHKERRQ(ierr)
      Call SectionRealRestrictClosure(AppCtx%FW%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, F_loc, iErr); CHKERRQ(ierr)
    
      Do_iGauss: Do iGauss = 1, Size(AppCtx%ElemVect(iE)%Gauss_C)
   
		U_eff_elem = 0.0_Kr  
		
		! Compute U_elem vector
		Do iDof=1, NumDoFVect
			U_eff_elem=U_eff_elem+AppCtx%ElemVect(iE)%BF(iDoF, iGauss) * (U_loc(iDoF)-U0_loc(iDoF))
		End Do
		Do iDoF = 1, NumDoFScal
			F_loc(iDoF) = F_loc(iDoF) +  (2_Kr * AppCtx%MatProp(iBlk_glob)%DelamToughness / AppCtx%MatProp(iBlk_glob)%Ksubst - (U_eff_elem .DotP. U_eff_elem) )* AppCtx%ElemScal(iE)%Gauss_C(iGauss) *  AppCtx%ElemScal(iE)%BF(iDoF, iGauss) 
		End Do
	End Do Do_iGauss
      
	Call SectionRealUpdateClosure(AppCtx%FW%Sec, AppCtx%MeshTopology%Mesh, iE-1, F_loc, ADD_VALUES, iErr); CHKERRQ(iErr)
   End Do Do_iEloc


   DeAllocate(U_loc)
   DeAllocate(U0_loc)
   DeAllocate(F_loc)      
!       Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyLocalV_Event, iErr); CHKERRQ(iErr)

End Subroutine FW_AssemblyBlk

Subroutine W_Solve(AppCtx)
   Type(AppCtx_Type)          :: AppCtx
   
   PetscInt             :: i, delamnodes
   PetscReal, Dimension(:), Pointer    :: Fi_ptr
   PetscInt             :: iErr
   PetscReal, Dimension(:), Pointer    :: zero
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer      
   
   Allocate(zero(1))
   
   zero=0.0_Kr
   delamnodes=0

   Do i=1, AppCtx%MeshTopology%Num_Verts
      Call SectionRealRestrict(AppCtx%FW%Sec, AppCtx%MeshTopology%Num_Elems + i-1, Fi_ptr, iErr); CHKERRQ(iErr);
      If (Fi_ptr(1) .LE. 0.0_Kr) Then
         Call SectionRealUpdate(AppCtx%W%Sec, AppCtx%MeshTopology%Num_Elems + i-1, zero, INSERT_VALUES, iErr); CHKERRQ(iErr)
   delamnodes=delamnodes+1
      End If
      Call SectionRealRestore(AppCtx%W%Sec, AppCtx%MeshTopology%Num_Elems+i-1, Fi_ptr, iErr); CHKERRQ(iErr)
   End Do
   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) "Number of debonded nodes (W=0)", delamnodes, "\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If   
   Deallocate(zero)
End Subroutine W_Solve

#undef __FUNCT__
#define __FUNCT__ "Step_W"
Subroutine Step_W(AppCtx)
   !!! Do one W-step in the alternate minimization algorithm
   Type(AppCtx_Type)          :: AppCtx

   PetscInt             :: iBlk
   PetscInt             :: NumDoFScal, NumDoFVect
   PetscInt             :: iErr
   Type(Vec)               :: FW_Vec
   Character(len=MEF90_MXSTRLEN)       :: IOBuffer
   Call PetscLogStagePush(AppCtx%LogInfo%WStep_Stage, iErr); CHKERRQ(iErr)

! Minimization wrt to W: W(x)=1 iff FW=>0
!
!  for each block
!     assembleblk:
!        for each elem
!        SecRestClo USec -> array of values
!        for each dof
!           locally compute u_eff_elem
!           locally compute F_loc=G-k*u_eff_elem . u_eff_elem 
!           SecRealUpdateClo F_loc -> FW%Sec
!     Complete
!     Test
!  Impose BC's
!     !!! Set Dirichlet Boundary Values
!     Call FieldInsertVertexBoundaryValues(AppCtx%RHSV, AppCtx%VBC, AppCtx%BCVFlag, AppCtx%MeshTopology)
!  
   
   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Assembling the F for the W-subproblem \n' 
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
   End If
   
   Call VecSet(AppCtx%FW%Vec, 0.0_Kr, iErr); CHKERRQ(iErr)
   
   Call FW_Assembly(AppCtx)
   
   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Solving for the W-subproblem\n' 
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
   End If
   
   Call W_Solve(AppCtx)
!  Impose BC
!  Call FieldInsertVertexBoundaryValues(AppCtx%W, AppCtx%WBC, AppCtx%BCWFlag, AppCtx%MeshTopology)
   
!  Log and print

   Call PetscLogStagePop(iErr); CHKERRQ(iErr)

End Subroutine Step_W


End Module m_VarFilmQS_W
