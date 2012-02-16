Module m_VarFilmQS_NL_U

#include "finclude/petscdef.h"

   Use m_VarFilmQS_Types
   Use m_VarFilmQS_Post
   Use m_MEF90
   Use m_Film_Struct
   Use m_VarFilmQS_W

   Implicit NONE
   Private   

   Public :: Init_TS_U
   Public :: HessianU_Assembly
   Public :: GradientU_Assembly
   Public :: Step_NL_U
   
Contains
Subroutine Init_TS_U(AppCtx)
   Type(AppCtx_Type)                            :: AppCtx
   PetscInt                                     :: iErr
   Character(len=MEF90_MXSTRLEN)                :: IOBuffer   
   
   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) "Initializing U with ", AppCtx%VarFracSchemeParam%InitV, "\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If      
   
   Select Case(AppCtx%VarFracSchemeParam%InitV)
   Case(VarFrac_INIT_V_ONE, VarFrac_INIT_V_OSC)
      Call SectionRealSet(AppCtx%U%Sec, 0.0_Kr, iErr); CHKERRQ(iErr)      
      Call VecSet(AppCtx%U%Vec, 0.0_Kr, iErr); CHKERRQ(iErr)      
   Case(VarFrac_INIT_V_FILE)
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_DisplacementX)%Offset, AppCtx%TimeStep, AppCtx%U)
      Call SectionRealToVec(AppCtx%U%Sec, AppCtx%U%Scatter, SCATTER_REVERSE, AppCtx%U%Vec, ierr); CHKERRQ(ierr)
   End Select
   !!! Update boundary values
   Call FieldInsertVertexBoundaryValues(AppCtx%U, AppCtx%UBC, AppCtx%BCUFlag, AppCtx%MeshTopology)
End Subroutine Init_TS_U
!!!
!!! Global Assembly Functions
!!! 

Subroutine HessianU_Assembly(H, AppCtx)
	Type(Mat)                                    :: H
	PetscInt                                     :: iBlk, iBlkID, iErr
	
	! log
	Call PetscLogStagePush(AppCtx%LogInfo%HessianUAssembly_Stage, iErr); CHKERRQ(iErr)
	! initialize 
	Call MatZeroEntries(H, iErr); CHKERRQ(iErr)
	
	! bc's
	!!! MatInsertVertexBoundaryValues overwrites the entire block corresponding to all
	!!! dof of a point where a boundary condition is applied
	!!! it is to be called BEFORE assembling the matrix
	Call MatInsertVertexBoundaryValues(H, AppCtx%U, AppCtx%BCUFlag, AppCtx%MeshTopology)
	Call MatAssemblyBegin(H, MAT_FLUSH_ASSEMBLY, iErr); CHKERRQ(iErr)
	Call MatAssemblyEnd  (H, MAT_FLUSH_ASSEMBLY, iErr); CHKERRQ(iErr)
	Do_Elem_iBlk: Do iBlk=1, AppCtx%MeshTopology%Num_Elem_Blks
		iBlkID=AppCtx%MeshTopology%Elem_Blk(iBlkID)%ID
		
		If ( AppCtx%MyEXO%EBProperty(VarFrac_EBProp_IsBrittle)%Value(iBlkID) /=0 ) Then
			Call HessianU_AssemblyBlk_Brittle(K, iBlk, AppCtx, .TRUE.)
		Else
			Call HessianU_AssemblyBlk_NonBrittle(K, iBlk, AppCtx, .TRUE.)
		End If
	End Do Do_Elem_iBlk
	
	Call MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
	Call MatAssemblyEnd  (H, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
	Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   
End Subroutine Hessian


Subroutine HessianU_AssemblyBlk_Brittle(H, iBlkID, AppCtx, DoBC)
	Type(Mat)                                    :: H
	PetscInt                                     :: iBlkID
	Type(AppCtx_Type)                            :: AppCtx
	PetscBool                                    :: DoBC
	
	PetscInt                                     :: iE, iEloc
	PetscReal, Dimension(:), Pointer             :: U_loc
	PetscReal, Dimension(:), Pointer             :: U0_loc
	PetscReal, Dimension(:), Pointer             :: V_loc
	PetscReal, Dimension(:, :), Pointer          :: H_loc
	PetscReal                                    :: U_elem
	PetscReal                                    :: U0_elem
	PetscReal                                    :: V_elem
	
	PetscInt                                     :: NumDoFScal, NumDoFVect
	PetscInt, Dimension(:), Pointer              :: BCFlag_Loc
	PetscInt                                     :: iDoF1, iDoF2, iGauss

	!init
	NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
	NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
	iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
	
	Allocate(U_loc(NumDoFVect))
	Allocate(U0_loc(NumDoFVect))
	Allocate(V_loc(NumDoFScal))
	Allocate(BCFlag_Loc(NumDoFVect))
	Allocate(H_loc(NumDoFVect, NumDoFVect))
	
	! loop over el in blk
	Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
		iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
		H_Loc  = 0.0_Kr
		! get u field, compute u_loc
		
		SectionRealRestrictClosure(AppCtx%U%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U_Loc, iErr)
		SectionRealRestrictClosure(AppCtx%U0%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U0_Loc, iErr)
		SectionRealRestrictClosure(AppCtx%V%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, V_Loc, iErr)
		! get BC field for U if needed
		If (DoBC) Then
			Call SectionIntRestrictClosure(AppCtx%BCUFlag%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, BCUFlag_Loc, iErr); CHKERRQ(ierr)
		End If
		Do_iGauss: Do iGauss = 1, Size(AppCtx%ElemVect(iE)%Gauss_C)
			U_Elem = 0.0_Kr        
			U0_Elem = 0.0_Kr        
			V_Elem = 0.0_Kr        
		
			Do iDof1=1, NumDoFScal
				V_Elem = V_Elem + AppCtx%ElemScal(iE)%BF(iDoF1, iGauss) * V_Loc(iDoF1)
			End Do
			
			CoefV = V_Elem**2 + AppCtx%VarFracSchemeParam%KEpsilon
			
			Do iDof1=1, NumDoFVect
				U_Elem = U_Elem + AppCtx%ElemVect(iE)%BF(iDoF1, iGauss) * U_Loc(iDoF1)
				U0_Elem = U0_Elem + AppCtx%ElemVect(iE)%BF(iDoF1, iGauss) * U0_Loc(iDoF1)
				
				! assemble local H matrix
				
				Do iDof2=1, NumDoFVect
					H_loc(iDof1, iDof2) = H_loc(iDof1, iDof2) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * CoefV * ((AppCtx%MatProp(iBlkID)%Hookes_Law * AppCtx%ElemVect(iE)%GradS_BF(iDoF1, iGauss)) .DotP. AppCtx%ElemVect(iE)%GradS_BF(iDoF2, iGauss))
					If ((U_Elem-U0_Elem) *  (U_Elem-U0_Elem) .LE. 2_Kr * AppCtx%MatProp(iBlk_glob)%DelamToughness / AppCtx%MatProp(iBlk_glob)%Ksubst) Then
						H_loc(iDof1, iDof2) = H_loc(iDof1, iDof2) + 0.5_Kr * AppCtx%MatProp(iBlk_glob)%Ksubst
					End If
				End Do
			End Do 
	         
		End Do Do_iGauss
		! assemble global
		
		Call DMMeshAssembleMatrix(AppCtx%H, AppCtx%MeshTopology%mesh, AppCtx%U%Sec, iE-1, H_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)

	End Do Do_Elem_iE
	
	! clean
	DeAllocate(H_Loc)
	DeAllocate(BCFlag_Loc)
	DeAllocate(V_loc)
	DeAllocate(U0_loc)
	DeAllocate(U_loc)
   
	
End Subroutine HessianU_AssemblyBlk_Brittle

Subroutine HessianU_AssemblyBlk_NonBrittle(H, iBlkID, AppCtx, DoBC)
	Type(Mat)                                    :: H
	PetscInt                                     :: iBlkID
	Type(AppCtx_Type)                            :: AppCtx
	PetscBool                                    :: DoBC
	
	PetscInt                                     :: iE, iEloc
	PetscReal, Dimension(:), Pointer             :: U_loc
	PetscReal, Dimension(:), Pointer             :: U0_loc
	PetscReal, Dimension(:, :), Pointer          :: H_loc
	PetscReal                                    :: U_elem
	PetscReal                                    :: U0_elem
	
	PetscInt                                     :: NumDoFVect
	PetscInt, Dimension(:), Pointer              :: BCFlag_Loc
	PetscInt                                     :: iDoF1, iDoF2, iGauss

	!init
	NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
	iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
	
	Allocate(U_loc(NumDoFVect))
	Allocate(U0_loc(NumDoFVect))
	Allocate(BCFlag_Loc(NumDoFVect))
	Allocate(H_loc(NumDoFVect, NumDoFVect))
	
	! loop over el in blk
	Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
		iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
		H_Loc  = 0.0_Kr
		! get u field, compute u_loc
		
		SectionRealRestrictClosure(AppCtx%U%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U_Loc, iErr)
		SectionRealRestrictClosure(AppCtx%U0%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U0_Loc, iErr)
		! get BC field for U if needed
		If (DoBC) Then
			Call SectionIntRestrictClosure(AppCtx%BCUFlag%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, BCUFlag_Loc, iErr); CHKERRQ(ierr)
		End If
		Do_iGauss: Do iGauss = 1, Size(AppCtx%ElemVect(iE)%Gauss_C)
			U_Elem = 0.0_Kr        
			U0_Elem = 0.0_Kr        
			
			Do iDof1=1, NumDoFVect
				U_Elem = U_Elem + AppCtx%ElemVect(iE)%BF(iDoF1, iGauss) * U_Loc(iDoF1)
				U0_Elem = U0_Elem + AppCtx%ElemVect(iE)%BF(iDoF1, iGauss) * U0_Loc(iDoF1)
				
				! assemble local H matrix
				
				Do iDof2=1, NumDoFVect
					H_loc(iDof1, iDof2) = H_loc(iDof1, iDof2) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * ((AppCtx%MatProp(iBlkID)%Hookes_Law * AppCtx%ElemVect(iE)%GradS_BF(iDoF1, iGauss)) .DotP. AppCtx%ElemVect(iE)%GradS_BF(iDoF2, iGauss))
					If ((U_Elem-U0_Elem) *  (U_Elem-U0_Elem) .LE. 2_Kr * AppCtx%MatProp(iBlk_glob)%DelamToughness / AppCtx%MatProp(iBlk_glob)%Ksubst) Then
						H_loc(iDof1, iDof2) = H_loc(iDof1, iDof2) + 0.5_Kr * AppCtx%MatProp(iBlk_glob)%Ksubst
					End If
				End Do
			End Do 
	         
		End Do Do_iGauss
		! assemble global
		
		Call DMMeshAssembleMatrix(AppCtx%H, AppCtx%MeshTopology%mesh, AppCtx%U%Sec, iE-1, H_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)

	End Do Do_Elem_iE
	
	! clean
	DeAllocate(H_Loc)
	DeAllocate(BCFlag_Loc)
	DeAllocate(U0_loc)
	DeAllocate(U_loc)

End Subroutine HessianU_AssemblyBlk_NonBrittle

Subroutine MatU_Assembly(K, AppCtx)
   Type(Mat)                                    :: K
   Type(AppCtx_Type)                            :: AppCtx
   
   PetscInt                                     :: iBlk, iBlkID, iErr
   
   Call PetscLogStagePush(AppCtx%LogInfo%MatAssemblyU_Stage, iErr); CHKERRQ(iErr)
   
   Call MatZeroEntries(K, iErr); CHKERRQ(iErr)
   !!! MatInsertVertexBoundaryValues overwrites the entire block corresponding to all
   !!! dof of a point where a boundary condition is applied
   !!! it is to be called BEFORE assembling the matrix
   Call MatInsertVertexBoundaryValues(K, AppCtx%U, AppCtx%BCUFlag, AppCtx%MeshTopology)
   Call MatAssemblyBegin(K, MAT_FLUSH_ASSEMBLY, iErr); CHKERRQ(iErr)
   Call MatAssemblyEnd  (K, MAT_FLUSH_ASSEMBLY, iErr); CHKERRQ(iErr)
   
   Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
      If (AppCtx%MyEXO%EBProperty(VarFrac_EBProp_IsBrittle)%Value(iBlkID) /= 0) Then
         Call MatU_AssemblyBlk_Brittle(K, iBlk, AppCtx%V%Sec, .TRUE., AppCtx)
      Else
         Call MatU_AssemblyBlk_NonBrittle(K, iBlk, .TRUE., AppCtx)
      End If
   End Do Do_Elem_iBlk
   Call MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
   Call MatAssemblyEnd  (K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
   If (AppCtx%AppParam%verbose > 2) Then
      Call MatView(K, PETSC_VIEWER_STDOUT_WORLD, iErr)
   End If
   
   Call PetscLogStagePop(iErr); CHKERRQ(iErr)
End Subroutine MatU_Assembly

   Subroutine RHSU_Assembly(RHS_Vec, AppCtx)
      Type(Vec)                                    :: RHS_Vec
      Type(AppCtx_Type)                            :: AppCtx
      
      Type(SectionReal)                            :: RHSU_Sec
      PetscInt                                     :: iErr
      PetscInt                                     :: iBlk, iBlkID

      Call PetscLogStagePush(AppCtx%LogInfo%RHSAssemblyU_Stage, iErr); CHKERRQ(iErr)

      Call SectionRealZero(AppCtx%RHSU%Sec, iErr); CHKERRQ(iErr)
      Do_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
         If (AppCtx%MyEXO%EBProperty(VarFrac_EBProp_IsBrittle)%Value(iBlkID) /= 0) Then
            Call RHS_AssemblyBlock_ElastBrittle(AppCtx%RHSU%Sec, iBlk, AppCtx)
         Else
            Call RHS_AssemblyBlock_ElastNonBrittle(AppCtx%RHSU%Sec, iBlk, AppCtx)
         End If         

      End Do Do_iBlk

      Call SectionRealComplete(AppCtx%RHSU%Sec, iErr); CHKERRQ(iErr)

      !!! Set Dirichlet Boundary Values
      Call FieldInsertVertexBoundaryValues(AppCtx%RHSU, AppCtx%UBC, AppCtx%BCUFlag, AppCtx%MeshTopology)

      Call SectionRealToVec(AppCtx%RHSU%Sec, AppCtx%RHSU%Scatter, SCATTER_FORWARD, RHS_Vec, iErr); CHKERRQ(ierr)

      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine RHSU_Assembly
  
  
   !!! 
   !!! Block Assembly Routines
   !!!
Subroutine MatU_AssemblyBlk_Brittle(K, iBlk, V_Sec, DoBC, AppCtx)
   Type(Mat)                                    :: K
   PetscInt                                     :: iBlk
   Type(SectionReal)                            :: V_Sec
   PetscBool                                    :: DoBC
   Type(AppCtx_Type)                            :: AppCtx
   
   PetscInt                                     :: iBlkID
   PetscReal, Dimension(:,:), Pointer           :: Mat_Loc      
   PetscInt                                     :: iELoc, iE
   PetscInt                                     :: iErr
   
   PetscInt                                     :: NumDoFScal, NumDoFVect, NumGauss
   PetscInt, Dimension(:), Pointer              :: BCFlag_Loc
   PetscInt                                     :: iDoF1, iDoF2, iGauss
   
   !!!   _Loc are restriction of fields to local patch (the element)
   !!!   _Elem are local contribution over the element (u_Elem = \sum_i U_Loc(i) BF(i))
   PetscReal, Dimension(:), Pointer             :: V_Loc, W_Loc
   PetscReal                                    :: V_Elem, CoefV, W_Elem
   PetscLogDouble, Parameter                    :: oneflop = 1.0
   
   Call PetscLogEventBegin(AppCtx%LogInfo%MatAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
   
   NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
   NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
   iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
   
   Allocate(V_Loc(NumDoFScal))
   Allocate(BCFlag_Loc(NumDoFVect))
   Allocate(W_Loc(NumDoFScal))
   BCFlag_Loc = VarFrac_BC_Type_NONE
   Allocate(Mat_Loc(NumDoFVect, NumDoFVect))

   Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
      iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
      
      Mat_Loc  = 0.0_Kr
      Call SectionRealRestrictClosure(V_Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, V_Loc, iErr); CHKERRQ(ierr)
      Call SectionRealRestrictClosure(AppCtx%W%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, W_Loc, iErr); CHKERRQ(ierr)
      If (DoBC) Then
         Call SectionIntRestrictClosure(AppCtx%BCUFlag%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, BCFlag_Loc, iErr); CHKERRQ(ierr)
      End If
         
      Do iGauss = 1, Size(AppCtx%ElemVect(iE)%Gauss_C)
         !!! Compute the contribution of V to the stiffness matrix
         !!! CoefV = (1+\eta_\varepsilon)v^2 if Is_Brittle, 1 otherwise
         !! Calculate V at the gauss point
         V_Elem = 0.0_Kr        
         W_Elem = 0.0_Kr        
         Do iDoF1 = 1, NumDoFScal
            V_Elem = V_Elem + AppCtx%ElemScal(iE)%BF(iDoF1, iGauss) * V_Loc(iDoF1)
            W_Elem = W_Elem + AppCtx%ElemScal(iE)%BF(iDoF1, iGauss) * W_Loc(iDoF1)
            Call PetscLogFlops(4*oneflop, iErr); CHKERRQ(iErr)
         End Do
         CoefV = V_Elem**2 + AppCtx%VarFracSchemeParam%KEpsilon
         Call PetscLogFlops(2*oneflop, iErr); CHKERRQ(iErr)
         Do iDoF1 = 1, NumDoFVect
            If (BCFlag_Loc(iDoF1) == VarFrac_BC_Type_NONE) Then
               Do iDoF2 = 1, NumDoFVect                  
                  Mat_Loc(iDoF2, iDoF1) =  Mat_Loc(iDoF2, iDoF1) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * CoefV * ((AppCtx%MatProp(iBlkID)%Hookes_Law * AppCtx%ElemVect(iE)%GradS_BF(iDoF1, iGauss)) .DotP. AppCtx%ElemVect(iE)%GradS_BF(iDoF2, iGauss))
                  Mat_Loc(iDoF2, iDoF1) =  Mat_Loc(iDoF2, iDoF1) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * AppCtx%MatProp(iBlkId)%Ksubst * W_Elem * (AppCtx%ElemVect(iE)%BF(iDoF1, iGauss) .DotP. AppCtx%ElemVect(iE)%BF(iDoF2, iGauss))
                  Call PetscLogFlops(3*oneflop, iErr); CHKERRQ(iErr)
               End Do
            End If
         End Do
      End Do
      Call DMMeshAssembleMatrix(AppCtx%KU, AppCtx%MeshTopology%mesh, AppCtx%U%Sec, iE-1, Mat_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)
   End Do Do_Elem_iE
      
   DeAllocate(Mat_Loc)
   DeAllocate(V_Loc)
   DeAllocate(W_Loc)
   DeAllocate(BCFlag_Loc)
   Call PetscLogEventEnd(AppCtx%LogInfo%MatAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
End Subroutine MatU_AssemblyBlk_Brittle

Subroutine MatU_AssemblyBlk_NonBrittle(K, iBlk, DoBC, AppCtx)
   Type(Mat)                                    :: K
   PetscInt                                     :: iBlk
   PetscBool                                    :: DoBC
   Type(AppCtx_Type)                            :: AppCtx
   
   PetscInt                                     :: iBlkID
   PetscReal, Dimension(:,:), Pointer           :: Mat_Loc      
   PetscInt                                     :: iELoc, iE
   PetscInt                                     :: iErr
   
   PetscInt                                     :: NumDoFScal, NumDoFVect, NumGauss
   PetscInt, Dimension(:), Pointer              :: BCFlag_Loc, W_Loc
   PetscReal               :: W_Elem
   PetscInt                                     :: iDoF1, iDoF2, iGauss
   
   !!!   _Loc are restriction of fields to local patch (the element)
   !!!   _Elem are local contribution over the element (u_Elem = \sum_i U_Loc(i) BF(i))
   PetscLogDouble, Parameter                    :: oneflop = 1.0
   
   Call PetscLogEventBegin(AppCtx%LogInfo%MatAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
   
   NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
   NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
   iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
   
   Allocate(BCFlag_Loc(NumDoFVect))
   Allocate(W_Loc(NumDoFScal))
   BCFlag_Loc = VarFrac_BC_Type_NONE
   Allocate(Mat_Loc(NumDoFVect, NumDoFVect))
   
   Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
      iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
      
      Mat_Loc  = 0.0_Kr
      If (DoBC) Then
         Call SectionIntRestrictClosure(AppCtx%BCUFlag%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, BCFlag_Loc, iErr); CHKERRQ(ierr)
      End If
      Do iGauss = 1, Size(AppCtx%ElemVect(iE)%Gauss_C)
         W_Elem=0.0_Kr
         Do iDoF1=1, NumDoFScal
            W_Elem = AppCtx%ElemScal(iE)%BF(iDoF1, iGauss) * W_Loc(iDoF1) 
         End Do
         Do iDoF1 = 1, NumDoFVect
            If (BCFlag_Loc(iDoF1) == VarFrac_BC_Type_NONE) Then
               Do iDoF2 = 1, NumDoFVect
                  
                  Mat_Loc(iDoF2, iDoF1) =  Mat_Loc(iDoF2, iDoF1) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * ((AppCtx%MatProp(iBlkID)%Hookes_Law * AppCtx%ElemVect(iE)%GradS_BF(iDoF1, iGauss)) .DotP. AppCtx%ElemVect(iE)%GradS_BF(iDoF2, iGauss))
                  Mat_Loc(iDoF2, iDoF1) =  Mat_Loc(iDoF2, iDoF1) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * AppCtx%MatProp(iBlkId)%Ksubst * W_Elem * (AppCtx%ElemVect(iE)%BF(iDoF1, iGauss) .DotP. AppCtx%ElemVect(iE)%BF(iDoF2, iGauss))
                  Call PetscLogFlops(3*oneflop, iErr); CHKERRQ(iErr)
                  
               End Do
            End If
         End Do
      End Do
      Call DMMeshAssembleMatrix(AppCtx%KU, AppCtx%MeshTopology%mesh, AppCtx%U%Sec, iE-1, Mat_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)
   End Do Do_Elem_iE
   
   Deallocate(Mat_Loc)
   Deallocate(BCFlag_Loc)
   Deallocate(W_Loc)
   Call PetscLogEventEnd(AppCtx%LogInfo%MatAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
End Subroutine MatU_AssemblyBlk_NonBrittle
   

Subroutine RHS_AssemblyBlock_ElastBrittle(RHS_Sec, iBlk, AppCtx)
   PetscInt                                     :: iBlk
   Type(SectionReal)                            :: RHS_Sec
   Type(AppCtx_Type)                            :: AppCtx

      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_ELem = \sum_i U_Loc(i) BF(i))
   PetscReal, Dimension(:), Pointer             :: Theta_Loc, RHS_Loc, W_Loc
   PetscInt, Dimension(:), Pointer              :: BCFlag_Loc
   PetscReal, Dimension(:), Pointer              :: U0_Loc
   Type(Vect2D)               :: U0_Elem
   PetscReal                                    :: Theta_Elem, W_Elem
   PetscInt                                     :: iE, iEloc, iBlkId, iErr
   PetscInt                                     :: NumDoFScal, NumDoFVect
   PetscInt                                     :: iDoF, iGauss
   PetscReal, Dimension(:), Pointer             :: V_Loc
   PetscReal                                    :: V_Elem, CoefV
   PetscLogDouble, Parameter                    :: oneflop  = 1.0
   PetscBool                                    :: Has_ThermExp
      
   Call PetscLogEventBegin(AppCtx%LogInfo%RHSAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
      
   NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
   NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
   iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
   If (Norm(AppCtx%MatProp(iBlkId)%Therm_Exp) > 0.0_Kr) Then
      Has_ThermExp = PETSC_TRUE
   Else
      Has_ThermExp = PETSC_FALSE
   End If

   If (Has_ThermExp) Then
      Allocate(V_Loc(NumDoFScal))
      Allocate(BCFlag_Loc(NumDoFVect))
      Allocate(RHS_Loc(NumDoFVect))
      Allocate(U0_Loc(NumDoFVect))
      Allocate(Theta_Loc(NumDoFScal))
      Allocate(W_Loc(NumDoFScal))
         
      Do_iEloc: Do iEloc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         RHS_Loc = 0.0_Kr
         Call SectionRealRestrictClosure(AppCtx%V%Sec,      AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, V_Loc,      iErr); CHKERRQ(ierr)
         Call SectionIntRestrictClosure(AppCtx%BCUFlag%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, BCFlag_Loc, iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(AppCtx%Theta%Sec,  AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, Theta_Loc,  iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(AppCtx%W%Sec,  AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, W_Loc,  iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(AppCtx%U0%Sec,  AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U0_Loc,  iErr); CHKERRQ(ierr)
         Do_iGauss: Do iGauss = 1, size(AppCtx%ElemVect(iE)%Gauss_C)
            V_Elem = 0.0_Kr        
            Theta_Elem = 0.0_Kr
            W_Elem = 0.0_Kr
            U0_Elem = 0.0_Kr
            Do iDoF = 1, NumDoFScal
               V_Elem = V_Elem + AppCtx%ElemScal(iE)%BF(iDoF, iGauss) * V_Loc(iDoF)
               Theta_Elem = Theta_Elem + AppCtx%ElemScal(iE)%BF(iDoF, iGauss) * Theta_Loc(iDoF)
               W_Elem = W_Elem + AppCtx%ElemScal(iE)%BF(iDoF, iGauss) * W_Loc(iDoF)
               Call PetscLogFlops(3*oneflop, iErr); CHKERRQ(iErr)
            End Do

            CoefV = V_Elem**2 + AppCtx%VarFracSchemeParam%KEpsilon
            !!! CoefV = (1+\eta_\varepsilon)v^2 if Is_Brittle, 1 otherwise
            Call PetscLogFlops(2*oneflop, iErr); CHKERRQ(iErr)
            Do iDoF = 1, NumDoFVect
               If (BCFlag_Loc(iDoF) == 0) Then
                     ! RHS terms due to inelastic strains and substrate
                  U0_Elem = U0_Elem + (U0_Loc(iDoF)*AppCtx%ElemVect(iE)%BF(iDoF, iGauss))
                  RHS_Loc(iDoF) = RHS_Loc(iDoF) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * Theta_Elem * CoefV * ((AppCtx%MatProp(iBlkId)%Hookes_Law * AppCtx%ElemVect(iE)%GradS_BF(iDoF, iGauss)) .DotP. AppCtx%MatProp(iBlkId)%Therm_Exp)
                  RHS_Loc(iDoF) = RHS_Loc(iDoF) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * AppCtx%MatProp(iBlkId)%Ksubst * (U0_Elem .DotP. AppCtx%ElemVect(iE)%BF(iDoF, iGauss)) * W_Elem
               Call PetscLogFlops(8*oneflop, iErr); CHKERRQ(iErr)
               End If
            End Do
         End Do Do_iGauss
         Call SectionRealUpdateClosure(RHS_Sec, AppCtx%MeshTopology%Mesh, iE-1, RHS_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iEloc
      DeAllocate(BCFlag_Loc)
      DeAllocate(RHS_Loc)
      DeAllocate(Theta_Loc)
      DeAllocate(V_Loc)
      DeAllocate(W_Loc)
      DeAllocate(U0_Loc)
   End If
      
   Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
End Subroutine RHS_AssemblyBlock_ElastBrittle

Subroutine RHS_AssemblyBlock_ElastNonBrittle(RHS_Sec, iBlk, AppCtx)
      PetscInt                                     :: iBlk
      Type(SectionReal)                            :: RHS_Sec
      Type(AppCtx_Type)                            :: AppCtx

      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_ELem = \sum_i U_Loc(i) BF(i))
      PetscReal, Dimension(:), Pointer             :: Theta_Loc, RHS_Loc, U0_Loc, W_Loc
      PetscInt, Dimension(:), Pointer              :: BCFlag_Loc
   Type(Vect2D)               :: U0_Elem
      PetscReal                                    :: Theta_Elem, W_Elem
      PetscInt                                     :: iE, iEloc, iBlkId, iErr
      PetscInt                                     :: NumDoFScal, NumDoFVect
      PetscInt                                     :: iDoF, iGauss
      PetscLogDouble, Parameter                    :: oneflop = 1.0
      PetscBool                                    :: Has_ThermExp
      
      Call PetscLogEventBegin(AppCtx%LogInfo%RHSAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
      
      NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
      NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
      If (Norm(AppCtx%MatProp(iBlkId)%Therm_Exp) > 0.0_Kr) Then
         Has_ThermExp = PETSC_TRUE
      Else
         Has_ThermExp = PETSC_FALSE
      End If
 
   If (Has_THermExp) Then
      Allocate(BCFlag_Loc(NumDoFVect))
      Allocate(RHS_Loc(NumDoFVect))
      Allocate(U0_Loc(NumDoFVect))
      Allocate(Theta_Loc(NumDoFScal))
      Allocate(W_Loc(NumDoFScal))
 
   Do_iEloc: Do iEloc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
      iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
      RHS_Loc = 0.0_Kr
      Call SectionIntRestrictClosure(AppCtx%BCUFlag%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, BCFlag_Loc, iErr); CHKERRQ(ierr)
      Call SectionRealRestrictClosure(AppCtx%Theta%Sec , AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, Theta_Loc,  iErr); CHKERRQ(ierr)
      Call SectionRealRestrictClosure(AppCtx%W%Sec , AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, W_Loc,  iErr); CHKERRQ(ierr)
      Call SectionRealRestrictClosure(AppCtx%U0%Sec , AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U0_Loc,  iErr); CHKERRQ(ierr)
      Do_iGauss: Do iGauss = 1, size(AppCtx%ElemVect(iE)%Gauss_C)
         Theta_Elem = 0.0_Kr
         W_Elem = 0.0_Kr
         U0_Elem = 0.0_Kr
         Do iDoF = 1, NumDoFScal
            Theta_Elem = Theta_Elem + AppCtx%ElemScal(iE)%BF(iDoF, iGauss) * Theta_Loc(iDoF)
            W_Elem = Theta_Elem + AppCtx%ElemScal(iE)%BF(iDoF, iGauss) * W_Loc(iDoF)
            Call PetscLogFlops(4*oneflop, iErr); CHKERRQ(iErr)
         End Do
         Do iDoF = 1, NumDoFVect
            If (BCFlag_Loc(iDoF) == 0) Then
                       ! RHS terms due to inelastic strains and substrate
               U0_Elem = U0_Elem + (U0_Loc(iDoF)*AppCtx%ElemVect(iE)%BF(iDoF, iGauss))
               RHS_Loc(iDoF) = RHS_Loc(iDoF) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * Theta_Elem * ((AppCtx%MatProp(iBlkId)%Hookes_Law * AppCtx%ElemVect(iE)%GradS_BF(iDoF, iGauss)) .DotP. AppCtx%MatProp(iBlkId)%Therm_Exp)
               RHS_Loc(iDoF) = RHS_Loc(iDoF) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * AppCtx%MatProp(iBlkId)%Ksubst * (U0_Elem .DotP. AppCtx%ElemVect(iE)%BF(iDoF, iGauss) ) * W_Elem
               Call PetscLogFlops(8*oneflop, iErr); CHKERRQ(iErr)
            End If
         End Do
      End Do Do_iGauss
      Call SectionRealUpdateClosure(RHS_Sec, AppCtx%MeshTopology%Mesh, iE-1, RHS_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iEloc
      Deallocate(BCFlag_Loc)
      Deallocate(RHS_Loc)
      Deallocate(Theta_Loc)
      Deallocate(W_Loc)
      Deallocate(U0_Loc)
   End If
      
   Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
End Subroutine RHS_AssemblyBlock_ElastNonBrittle
   

Subroutine Step_U(AppCtx)
   Type(AppCtx_Type)                            :: AppCtx
   
   PetscInt                                     :: iErr
   KSPConvergedReason                           :: reason
   PetscInt                                     :: KSPNumIter
   Character(len=MEF90_MXSTRLEN)                :: IOBuffer
   PetscReal                                    :: VMin, VMax
   PetscReal                                    :: rDum
   PetscInt                                     :: iDum

   Call PetscLogStagePush(AppCtx%LogInfo%UStep_Stage, iErr); CHKERRQ(iErr)
  
   Call SectionRealToVec(AppCtx%U%Sec, AppCtx%U%Scatter, SCATTER_FORWARD, AppCtx%U%Vec, ierr); CHKERRQ(ierr)
   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Assembling the Matrix and RHS for the U-subproblem \n' 
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
   End If
   Call RHSU_Assembly(AppCtx%RHSU%Vec, AppCtx)
   If (AppCtx%AppParam%verbose > 2) Then
      Call VecView(AppCtx%RHSU%Vec, AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
   End If
   
   Call MatU_Assembly(AppCtx%KU, AppCtx)
   If (AppCtx%AppParam%verbose > 2) Then
      Call MatView(AppCtx%KU, AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
   End If
   
   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Calling KSPSolve for the U-subproblem\n' 
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
   End If
   
   Call KSPSolve(AppCtx%KSPU, AppCtx%RHSU%Vec, AppCtx%U%Vec, iErr); CHKERRQ(iErr)
   Call SectionRealToVec(AppCtx%U%Sec, AppCtx%U%Scatter, SCATTER_REVERSE, AppCtx%U%Vec, ierr); CHKERRQ(ierr)
   
   Call KSPGetConvergedReason(AppCtx%KSPU, reason, iErr); CHKERRQ(iErr)
   If ( reason > 0) Then
      Call KSPGetIterationNumber(AppCtx%KSPU, KSPNumIter, iErr); CHKERRQ(iErr)
      Write(IOBuffer, 100) KSPNumIter, reason
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Else
      Write(IOBuffer, 101) reason
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      If (AppCtx%AppParam%StopOnError) Then
         SETERRQ(PETSC_COMM_SELF,PETSC_ERR_CONV_FAILED, 'KSP failed to converge, aborting...\n', iErr)
      EndIf
   End If
   Call PetscLogStagePop(iErr); CHKERRQ(iErr)
100 Format('     KSP for U converged in  ', I5, ' iterations. KSPConvergedReason is    ', I5, '\n')
101 Format('[ERROR] KSP for U diverged. KSPConvergedReason is ', I2, '\n')
   End Subroutine Step_U   



End Module m_VarFilmQS_NL_U

