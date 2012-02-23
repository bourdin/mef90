Module m_VarFilmQS_NL_U

#include "finclude/petscdef.h"
      use petscsnes
   Use m_VarFilmQS_Types
   Use m_VarFilmQS_Post
   Use m_MEF90
   Use m_Film_Struct
   Use m_VarFilmQS_W

   Implicit NONE
   Private   

   Public :: Init_TS_NL_U
   Public :: HessianU_Assembly
   Public :: GradientU_Assembly
   
Contains
#undef __FUNC__ 
#define __FUNC__ "Init_TS_U"
Subroutine Init_TS_NL_U(AppCtx)
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
End Subroutine Init_TS_NL_U
!!!
!!! Global Assembly Functions
!!! 
#undef __FUNC__ 
#define __FUNC__ "HessianU_Assembly"
Subroutine HessianU_Assembly(SNESappU, U_Vec, H, HPC, flag, AppCtx)
	Type(SNES)                                   :: SNESappU
	Type(Vec)                                    :: U_Vec
	Type(Mat)                                    :: H, HPC
	Type(AppCtx_Type)                            :: AppCtx
	MatStructure                                 :: flag
      
   
	PetscInt                                     :: iBlk, iBlkID, iErr
	
	! log
! 	Call PetscLogStagePush(AppCtx%LogInfo%HessianUAssembly_Stage, iErr); CHKERRQ(iErr)
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
		
		If ( AppCtx%MyEXO%EBProperty(VarFrac_EBProp_IsBrittle)%Value(iBlkID) /= 0 ) Then
			Call HessianU_AssemblyBlk_Brittle(H, iBlk, AppCtx, .TRUE.)
		Else
			Call HessianU_AssemblyBlk_NonBrittle(H, iBlk, AppCtx, .TRUE.)
		End If
	End Do Do_Elem_iBlk
	
	Call MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
	Call MatAssemblyEnd  (H, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
	Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   
End Subroutine HessianU_Assembly

#undef __FUNC__ 
#define __FUNC__ "GradientU_Assembly"
Subroutine GradientU_Assembly(SNESappU, U_Vec, GradientU, AppCtx)
	Type(SNES)                                   :: SNESappU
	Type(Vec)                                    :: U_Vec
	Type(Vec)                                    :: GradientU
	Type(AppCtx_Type)                            :: AppCtx
	
	PetscInt                                     :: iBlk, iBlkID, iErr
	
	! log
	! 
	Do_Elem_iBlk: Do iBlk=1, AppCtx%MeshTopology%Num_Elem_Blks
		iBlkID=AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
		
		If ( AppCtx%MyEXO%EBProperty(VarFrac_EBProp_IsBrittle)%Value(iBlkID) /= 0 ) Then
			Call GradientU_AssemblyBlk_Brittle(GradientU, iBlk, AppCtx)
		Else
			Call GradientU_AssemblyBlk_NonBrittle(GradientU, iBlk, AppCtx)
		End If
	End Do Do_Elem_iBlk
	
End Subroutine GradientU_Assembly
 
   !!! 
   !!! Block Assembly Routines
   !!!
#undef __FUNC__ 
#define __FUNC__ "HessianU_AssemblyBlk_Brittle"
Subroutine HessianU_AssemblyBlk_Brittle(H, iBlkID, AppCtx, DoBC)
	Type(Mat)                                    :: H
	PetscInt                                     :: iBlk, iBlkID
	Type(AppCtx_Type)                            :: AppCtx
	PetscBool                                    :: DoBC
	PetscInt                                     :: iErr
   
	PetscInt                                     :: iE, iEloc
	PetscReal, Dimension(:), Pointer             :: U_loc, U0_loc, V_loc
	PetscInt, Dimension(:), Pointer              :: BCUFlag_Loc
	PetscReal, Dimension(:, :), Pointer          :: H_loc
	Type(Vect2D)                                  :: U_elem, U0_elem, U_eff_elem
	PetscReal                                    :: V_elem
	PetscReal                                    :: CoefV
	
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
	Allocate(BCUFlag_Loc(NumDoFVect))
	Allocate(H_loc(NumDoFVect, NumDoFVect))
	
	! loop over el in blk
	Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
		iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
		H_Loc  = 0.0_Kr
		! get u field, compute u_loc		
		Call SectionRealRestrictClosure(AppCtx%U%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U_Loc, iErr)
		Call SectionRealRestrictClosure(AppCtx%U0%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U0_Loc, iErr)
		Call SectionRealRestrictClosure(AppCtx%V%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, V_Loc, iErr)
		! get BC field for U if needed
		If (DoBC) Then
			Call SectionIntRestrictClosure(AppCtx%BCUFlag%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, BCUFlag_Loc, iErr); CHKERRQ(ierr)
		End If
		Do_iGauss: Do iGauss = 1, Size(AppCtx%ElemVect(iE)%Gauss_C)
			U_Elem = 0.0_Kr        
			U0_Elem = 0.0_Kr       
			U_eff_elem = 0.0_Kr
			
			V_Elem = 0.0_Kr        
		
			Do iDof1=1, NumDoFScal
				V_Elem = V_Elem + AppCtx%ElemScal(iE)%BF(iDoF1, iGauss) * V_Loc(iDoF1)
			End Do
			
			CoefV = V_Elem**2 + AppCtx%VarFracSchemeParam%KEpsilon
			
			Do iDof1=1, NumDoFVect
! 				U_Elem = U_Elem + U_Loc(iDoF1) * AppCtx%ElemVect(iE)%BF(iDoF1, iGauss) 
! 				U0_Elem = U0_Elem + U0_Loc(iDoF1) * AppCtx%ElemVect(iE)%BF(iDoF1, iGauss) 
				U_eff_elem = U_eff_elem + AppCtx%ElemVect(iE)%BF(iDoF1, iGauss) * (U_Loc(iDoF1) - U0_Loc(iDoF1))
				! assemble local H matrix
				
				Do iDof2=1, NumDoFVect
					H_loc(iDof1, iDof2) = H_loc(iDof1, iDof2) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * CoefV * ((AppCtx%MatProp(iBlkID)%Hookes_Law * AppCtx%ElemVect(iE)%GradS_BF(iDoF1, iGauss)) .DotP. AppCtx%ElemVect(iE)%GradS_BF(iDoF2, iGauss))
					If ((U_eff_elem .DotP. U_eff_elem) .LE. 2_Kr * AppCtx%MatProp(iBlkID)%DelamToughness / AppCtx%MatProp(iBlkID)%Ksubst) Then
						H_loc(iDof1, iDof2) = H_loc(iDof1, iDof2) + 0.5_Kr * AppCtx%MatProp(iBlkID)%Ksubst
					End If
				End Do
			End Do 
	         
		End Do Do_iGauss
		! assemble global
		
		Call DMMeshAssembleMatrix(AppCtx%KU, AppCtx%MeshTopology%mesh, AppCtx%U%Sec, iE-1, H_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)

	End Do Do_Elem_iE
	
	! clean
	DeAllocate(H_Loc)
	DeAllocate(BCFlag_Loc)
	DeAllocate(V_loc)
	DeAllocate(U0_loc)
	DeAllocate(U_loc)
   
	
End Subroutine HessianU_AssemblyBlk_Brittle
#undef __FUNC__ 
#define __FUNC__ "HessianU_AssemblyBlk_NonBrittle"
Subroutine HessianU_AssemblyBlk_NonBrittle(H, iBlkID, AppCtx, DoBC)
	Type(Mat)                                    :: H
	PetscInt                                     :: iBlk, iBlkID
	Type(AppCtx_Type)                            :: AppCtx
	PetscBool                                    :: DoBC
	PetscInt                                     :: iErr
   
	PetscInt                                     :: iE, iEloc
	PetscReal, Dimension(:), Pointer             :: U_loc, U0_loc
	PetscReal, Dimension(:, :), Pointer          :: H_loc
	PetscInt, Dimension(:), Pointer              :: BCUFlag_Loc
	
	Type(Vect2D)                                 :: U_elem, U0_elem, U_eff_elem
	
	PetscInt                                     :: NumDoFVect
	PetscInt, Dimension(:), Pointer              :: BCFlag_Loc
	PetscInt                                     :: iDoF1, iDoF2, iGauss

	!init
	NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
	iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
	
	Allocate(U_loc(NumDoFVect))
	Allocate(U0_loc(NumDoFVect))
	Allocate(BCUFlag_Loc(NumDoFVect))
	Allocate(H_loc(NumDoFVect, NumDoFVect))
	
	! loop over el in blk
	Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
		iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
		H_Loc  = 0.0_Kr
		! get u field, compute u_loc
		
		Call SectionRealRestrictClosure(AppCtx%U%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U_Loc, iErr)
		Call SectionRealRestrictClosure(AppCtx%U0%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U0_Loc, iErr)
		! get BC field for U if needed
		If (DoBC) Then
			Call SectionIntRestrictClosure(AppCtx%BCUFlag%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, BCUFlag_Loc, iErr); CHKERRQ(ierr)
			! do bcs
		End If
		Do_iGauss: Do iGauss = 1, Size(AppCtx%ElemVect(iE)%Gauss_C)
			U_Elem = 0.0_Kr        
			U0_Elem = 0.0_Kr        
			U_eff_Elem = 0.0_Kr        
			
			Do iDof1=1, NumDoFVect
! 				U_Elem = U_Elem + U_Loc(iDoF1) * AppCtx%ElemVect(iE)%BF(iDoF1, iGauss) 
! 				U0_Elem = U0_Elem + U0_Loc(iDoF1) * AppCtx%ElemVect(iE)%BF(iDoF1, iGauss) 
				U_eff_elem = U_eff_elem + AppCtx%ElemVect(iE)%BF(iDoF1, iGauss) * (U_Loc(iDoF1) - U0_Loc(iDoF1))
				! assemble local H matrix
				
				Do iDof2=1, NumDoFVect
					H_loc(iDof1, iDof2) = H_loc(iDof1, iDof2) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * ((AppCtx%MatProp(iBlkID)%Hookes_Law * AppCtx%ElemVect(iE)%GradS_BF(iDoF1, iGauss)) .DotP. AppCtx%ElemVect(iE)%GradS_BF(iDoF2, iGauss))
					If ((U_eff_elem .DotP. U_eff_elem) .LE. 2_Kr * AppCtx%MatProp(iBlkID)%DelamToughness / AppCtx%MatProp(iBlkID)%Ksubst) Then
						H_loc(iDof1, iDof2) = H_loc(iDof1, iDof2) + 0.5_Kr * AppCtx%MatProp(iBlkID)%Ksubst
					End If
				End Do
			End Do 
	         
		End Do Do_iGauss
		! assemble global
		
		Call DMMeshAssembleMatrix(AppCtx%KU, AppCtx%MeshTopology%mesh, AppCtx%U%Sec, iE-1, H_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)

	End Do Do_Elem_iE
	
	! clean
	DeAllocate(H_Loc)
	DeAllocate(BCFlag_Loc)
	DeAllocate(U0_loc)
	DeAllocate(U_loc)

End Subroutine HessianU_AssemblyBlk_NonBrittle
#undef __FUNC__ 
#define __FUNC__ "GradientU_AssemblyBlk_Brittle"
Subroutine GradientU_AssemblyBlk_Brittle(GradientU, iBlk, AppCtx)
	Type(Vec)                                    :: GradientU
	PetscInt                                     :: iBlk, iBlkID
	Type(AppCtx_Type)                            :: AppCtx
	PetscInt                                     :: iErr
   
	PetscInt                                     :: iE, iEloc
	PetscReal, Dimension(:), Pointer             :: U_loc, U0_loc, V_loc, Theta_loc
	Type(MatS2D)                                 :: Strain_Elem, EffectiveStrain_Elem, Stress_elem
	PetscReal, Dimension(:), Pointer             :: GradientU_loc
	Type(Vect2D)                                 :: U_elem, U_eff_elem, U0_elem
	PetscReal                                    :: V_elem, Theta_elem, CoefV
	PetscInt                                     :: NumDoFScal, NumDoFVect, iDoF, iGauss
	PetscBool                                    :: Has_ThermExp

	! log
	! init
	NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
	NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems * AppCtx%MeshTopology%Num_Dim

	If (Norm(AppCtx%MatProp(iBlkId)%Therm_Exp) > 0.0_Kr) Then
		Has_ThermExp = PETSC_TRUE
	Else                 
		Has_ThermExp = PETSC_FALSE
	End If	
	Allocate(U_loc(NumDoFVect))
	Allocate(U0_loc(NumDoFVect))
	Allocate(V_loc(NumDoFScal))
	Allocate(Theta_loc(NumDoFScal))
	Allocate(GradientU_loc(NumDoFVect))
	
	! loop over elems in block
	Do_Elem_iE: Do iEloc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
		iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
		GradientU_loc = 0.0_Kr
		! get values of the fields over the current patch (element)
		Call SectionRealRestrictClosure(AppCtx%U%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U_Loc, iErr)
		Call SectionRealRestrictClosure(AppCtx%U0%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U0_Loc, iErr)
		Call SectionRealRestrictClosure(AppCtx%V%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, V_Loc, iErr) 
		Call SectionRealRestrictClosure(AppCtx%Theta%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, Theta_Loc, iErr) 
		Do_iGauss: Do iGauss = 1, Size(AppCtx%ElemVect(iE)%Gauss_C)
			U_elem = 0.0_Kr
			U0_elem = 0.0_Kr
			V_elem = 0.0_Kr
			Strain_elem = 0.0_Kr
			Stress_elem = 0.0_Kr
			EffectiveStrain_elem = 0.0_Kr
	
			! compute V at current Gauss pt 
			Do iDof = 1, NumDoFScal
				V_elem = V_elem + V_Loc(iDoF) * AppCtx%ElemScal(iE)%BF(iDoF, iGauss)
				Theta_elem = Theta_elem + Theta_Loc(iDoF) * AppCtx%ElemScal(iE)%BF(iDoF, iGauss)
			End Do
			CoefV = V_elem**2 + AppCtx%VarFracSchemeParam%KEpsilon
			
			! compute U-U0 and Strain_elem and EffectiveStrain_elem at current Gauss pt 
			Do iDof = 1, NumDoFVect 
				U_eff_elem = U_elem + (U_Loc(iDoF)-U0_Loc(iDoF)) * AppCtx%ElemVect(iE)%BF(iDoF, iGauss)
				Strain_elem = Strain_elem + U_Loc(iDoF) * AppCtx%ElemVect(iE)%GradS_BF(iDoF, iGauss) 
				If(Has_ThermExp) Then
					EffectiveStrain_elem = Strain_elem - (Theta_Elem * AppCtx%MatProp(iBlkID)%Therm_Exp)
				Else
					EffectiveStrain_elem = Strain_elem
					Stress_elem = AppCtx%MatProp(iBlkID)%Hookes_Law * EffectiveStrain_elem 
				End If
			End Do
			
			! compute GradE_U vector
			Do iDof = 1, NumDoFVect
				GradientU_loc(iDof) = GradientU_loc(iDof) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * CoefV * (Stress_elem .DotP. AppCtx%ElemVect(iE)%GradS_BF(iDof, iGauss)) 
				If ((U_eff_elem .DotP. U_eff_elem ) .LT. 2.0_Kr * AppCtx%MatProp(iBlkID)%DelamToughness / AppCtx%MatProp(iBlkID)%Ksubst) Then
				GradientU_loc(iDof) = GradientU_loc(iDof) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * 0.5_Kr * AppCtx%MatProp(iBlkID)%Ksubst * (U_eff_elem .DotP. AppCtx%ElemVect(iE)%BF(iDof, iGauss) )
				End If
			End Do
		End Do Do_iGauss
	End Do Do_Elem_iE
	
End Subroutine GradientU_AssemblyBlk_Brittle
#undef __FUNC__ 
#define __FUNC__ "GradientU_AssemblyBlk_NonBrittle"
Subroutine GradientU_AssemblyBlk_NonBrittle(GradientU, iBlk, AppCtx)
	Type(Vec)                                    :: GradientU
	PetscInt                                     :: iBlk, iBlkID
	Type(AppCtx_Type)                            :: AppCtx
	PetscInt                                     :: iErr
   
	PetscInt                                     :: iE, iEloc
	PetscReal, Dimension(:), Pointer             :: U_loc, U0_loc, Theta_loc
	Type(MatS2D)                                 :: Strain_Elem, EffectiveStrain_Elem, Stress_elem
	PetscReal, Dimension(:), Pointer             :: GradientU_loc
	Type(Vect2D)                                  :: U_elem, U0_elem, U_eff_elem
	PetscReal                                    :: V_elem, Theta_elem
	PetscBool                                    :: Has_ThermExp
      
	PetscInt                                     :: NumDoFScal, NumDoFVect, iDoF, iGauss

	! log
	! init
	NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
	NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems * AppCtx%MeshTopology%Num_Dim

	If (Norm(AppCtx%MatProp(iBlkId)%Therm_Exp) > 0.0_Kr) Then
		Has_ThermExp = PETSC_TRUE
	Else
		Has_ThermExp = PETSC_FALSE
	End If	
	Allocate(U_loc(NumDoFVect))
	Allocate(U0_loc(NumDoFVect))
	Allocate(Theta_loc(NumDoFScal))
	Allocate(GradientU_loc(NumDoFVect))
	
	! loop over elems in block
	Do_Elem_iE: Do iEloc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
		iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
		GradientU_loc = 0.0_Kr
		! get values of the fields over the current patch (element)
		Call SectionRealRestrictClosure(AppCtx%U%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U_Loc, iErr)
		Call SectionRealRestrictClosure(AppCtx%U0%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U0_Loc, iErr)
		Call SectionRealRestrictClosure(AppCtx%Theta%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, Theta_Loc, iErr) 
		Do_iGauss: Do iGauss = 1, Size(AppCtx%ElemVect(iE)%Gauss_C)
			U_elem = 0.0_Kr
			U0_elem = 0.0_Kr
			Strain_elem = 0.0_Kr
			Stress_elem = 0.0_Kr
			EffectiveStrain_elem = 0.0_Kr
	
			! compute V at current Gauss pt 
			Do iDof = 1, NumDoFScal
				Theta_elem = Theta_elem + AppCtx%ElemScal(iE)%BF(iDoF, iGauss) * Theta_Loc(iDoF)
			End Do
			
			! compute U-U0 and Strain_elem and EffectiveStrain_elem at current Gauss pt 
			Do iDof = 1, NumDoFVect 
				U_eff_elem = U_elem + AppCtx%ElemVect(iE)%BF(iDoF, iGauss) * (U_Loc(iDoF)-U0_Loc(iDoF))
				Strain_elem = Strain_elem + AppCtx%ElemVect(iE)%GradS_BF(iDoF, iGauss) * U_Loc(iDoF)
				If(Has_ThermExp) Then
					EffectiveStrain_elem = Strain_elem - (Theta_Elem * AppCtx%MatProp(iBlkID)%Therm_Exp)
				Else
					EffectiveStrain_elem = Strain_elem
					Stress_elem = AppCtx%MatProp(iBlkID)%Hookes_Law * EffectiveStrain_elem
				End If
			End Do
			
			! compute GradE_U vector
			Do iDof = 1, NumDoFVect
				GradientU_loc = GradientU_loc + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * (Stress_elem .DotP. AppCtx%ElemVect(iE)%GradS_BF(iDof, iGauss))
				If ((U_eff_elem .DotP. U_eff_elem) .LT. 2.0_Kr * AppCtx%MatProp(iBlkID)%DelamToughness / AppCtx%MatProp(iBlkID)%Ksubst) Then
					GradientU_loc = GradientU_loc + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * 0.5_Kr * AppCtx%MatProp(iBlkID)%Ksubst *  (U_eff_elem .DotP. AppCtx%ElemVect(iE)%BF(iDof, iGauss)) 
				End If
			End Do
		End Do Do_iGauss
	End Do Do_Elem_iE
	
End Subroutine GradientU_AssemblyBlk_NonBrittle

End Module m_VarFilmQS_NL_U

