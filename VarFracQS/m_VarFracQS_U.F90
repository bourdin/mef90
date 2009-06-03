#if defined PB_2D
Module m_VarFracQS_U2D
#elif defined PB_3D
Module m_VarFracQS_U3D
#endif

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscmeshdef.h"
#include "finclude/petscviewerdef.h"

#if defined PB_2D
   Use m_VarFracQS_Types2D
#elif defined PB_3D
   Use m_VarFracQS_Types3D
#endif   
   Use m_MEF90
   Use m_VarFrac_Struct
   Use petsc
   Use petscvec
   Use petscmat
   Use petscksp
   Use petscmesh

   Implicit NONE   
   
Contains
   Subroutine Init_TS_U(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      Type(SectionReal)                            :: UBC
      PetscReal, Dimension(:), Pointer             :: U_Ptr, UBC_Ptr
      PetscInt, Dimension(:), Pointer              :: BCFlag
      PetscInt                                     :: i, j, iErr
      
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'UBC', AppCtx%MeshTopology%Num_Dim, UBC, iErr); CHKERRQ(iErr)
      !!! Using SectionRealDuplicate would make more sense
      
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_DisplacementX)%Offset, AppCtx%TimeStep, UBC)
      Allocate(BCFlag(AppCtx%MeshTopology%Num_Dim))
      Allocate(U_Ptr(AppCtx%MeshTopology%Num_Dim))
      Allocate(UBC_Ptr(AppCtx%MeshTopology%Num_Dim))

      Do i = 1, AppCtx%MeshTopology%Num_Verts
         Call MeshRestrictClosureInt(AppCtx%MeshTopology%mesh, AppCtx%BCUFlag, AppCtx%MeshTopology%Num_Elems + i-1, AppCtx%MeshTopology%Num_Dim, BCFlag, iErr); CHKERRQ(ierr)
         If (Sum(BCFlag) /= 0) Then
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%U, AppCtx%MeshTopology%Num_Elems + i-1, AppCtx%MeshTopology%Num_Dim, U_Ptr, iErr); CHKERRQ(ierr)      
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, UBC, AppCtx%MeshTopology%Num_Elems + i-1, AppCtx%MeshTopology%Num_Dim, UBC_Ptr, iErr); CHKERRQ(ierr)
            Do j = 1, AppCtx%MeshTopology%Num_Dim      
               If (BCFlag(j) /= 0) Then
                  U_Ptr(j) = UBC_Ptr(j)
               End If
            End Do
            Call MeshUpdateClosure(AppCtx%MeshTopology%mesh, AppCtx%U, AppCtx%MeshTopology%Num_Elems + i-1, UBC_Ptr, iErr); CHKERRQ(ierr)
         End If
      End Do
      DeAllocate(BCFlag)
      DeAllocate(U_Ptr)
      DeAllocate(UBC_Ptr)
      
      Call SectionRealDestroy(UBC, iErr); CHKERRQ(iErr)
   End Subroutine Init_TS_U

   !----------------------------------------------------------------------------------------!      
   ! MatAssembly (CM)  
   !----------------------------------------------------------------------------------------!      
   Subroutine MatU_Assembly(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iBlk, iErr
      
      Call PetscLogStagePush(AppCtx%LogInfo%MatAssemblyU_Stage, iErr); CHKERRQ(iErr)
      
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Call MatU_AssemblyBlk(iBlk, AppCtx)
      End Do Do_Elem_iBlk
      Call MatAssemblyBegin(AppCtx%KU, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (AppCtx%KU, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine MatU_Assembly

   Subroutine MatU_AssemblyBlk(iBlk, AppCtx)
      PetscInt                                     :: iBlk
      Type(AppCtx_Type)                            :: AppCtx

      PetscInt                                     :: iBlkID
      PetscReal, Dimension(:,:), Pointer           :: MatElem      
      PetscInt                                     :: iELoc, iE
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoFScal, NumDoFVect, NumGauss
      PetscInt, Dimension(:), Pointer              :: BCFlag
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      
      PetscReal, Dimension(:), Pointer             :: V
      PetscReal                                    :: V_Elem, CoefV
      PetscLogDouble                               :: flops = 0
      
      Call PetscLogEventBegin(AppCtx%LogInfo%MatAssemblyLocalU_Event, iErr); CHKERRQ(iErr)

      NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
      NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
      
      Allocate(V(NumDoFScal))
      Allocate(BCFlag(NumDoFVect))
      Allocate(MatElem(NumDoFVect, NumDoFVect))
      
      Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         MatElem  = 0.0_Kr
         Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%V, iE-1, NumDoFScal, V, iErr); CHKERRQ(ierr)
         Call MeshRestrictClosureInt(AppCtx%MeshTopology%mesh, AppCtx%BCUFlag, iE-1, NumDoFVect, BCFlag, iErr); CHKERRQ(ierr)
      
         Do iGauss = 1, Size(AppCtx%ElemVect(iE)%Gauss_C)
            !!! Compute the contribution of V to the stiffness matrix
            !!! CoefV = (1+\eta_\varepsilon)v^2 if Is_Brittle, 1 otherwise
            If (AppCtx%MyEXO%EBProperty(VarFrac_EBProp_IsBrittle)%Value(iBlkID) /= 0) Then
            !! Calculate V at the gauss point
               V_Elem = 0.0_Kr        
               Do iDoF1 = 1, NumDoFScal
                  V_Elem = V_Elem + AppCtx%ElemScal(iE)%BF(iDoF1, iGauss) * V(iDoF1)
                  flops = flops + 2
               End Do
               CoefV = V_Elem**2 + AppCtx%VarFracSchemeParam%KEpsilon
               flops = flops + 2
            Else
               CoefV = 1.0_Kr
            End If
            !! Assemble the element stiffness
            Do iDoF1 = 1, NumDoFVect
               If (BCFlag(iDoF1) == 0) Then
                  Do iDoF2 = 1, NumDoFVect
                     MatElem(iDoF2, iDoF1) =  MatElem(iDoF2, iDoF1) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * (CoefV * (AppCtx%MatProp(iBlkID)%Hookes_Law * AppCtx%ElemVect(iE)%GradS_BF(iDoF1, iGauss)) .DotP. AppCtx%ElemVect(iE)%GradS_BF(iDoF2, iGauss))
                     flops = flops + 3
                  End Do
               End If
            End Do
         End Do
         Call assembleMatrix(AppCtx%KU, AppCtx%MeshTopology%mesh, AppCtx%U, iE-1, MatElem, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_Elem_iE
      
      Call PetscLogFlops(flops, iErr); CHKERRQ(iErr)
      DeAllocate(MatElem)
      DeAllocate(V)
      DeAllocate(BCFlag)
      Call PetscLogEventEnd(AppCtx%LogInfo%MatAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
   End Subroutine MatU_AssemblyBlk

   !----------------------------------------------------------------------------------------!      
   ! RHSAssembly (CM)  
   !----------------------------------------------------------------------------------------!      
   Subroutine RHSU_Assembly(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      Type(SectionReal)                            :: RHSSec
      PetscInt                                     :: iBlk, iE, iELoc
      PetscReal, Dimension(:), Pointer             :: RHSElem
      PetscInt                                     :: NumDoFPerVertex 
      
      NumDoFPerVertex = AppCtx%MeshTopology%Num_Dim
      
      !!! Hopefully one day we will use assemble Vector instead of going through a section
      Call PetscLogStagePush(AppCtx%LogInfo%RHSAssemblyU_Stage, iErr); CHKERRQ(iErr)
      Call VecSet(AppCtx%RHSU, 0.0_Kr, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'RHSSec', NumDoFPerVertex, RHSSec, iErr); CHKERRQ(iErr)
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Allocate(RHSElem(AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim))
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            !--------------------------------------------------------
            Call RHSU_AssemblyLocal(iE, AppCtx%MatProp( AppCtx%MeshTopology%Elem_Blk(iBlk)%ID ),AppCtx, RHSElem)
            !--------------------------------------------------------
            Call MeshUpdateAddClosure(AppCtx%MeshTopology%Mesh, RHSSec, iE-1, RHSElem, iErr); CHKERRQ(iErr)
         End Do Do_Elem_iE
         DeAllocate(RHSElem)
      End Do Do_Elem_iBlk
      Call SectionRealComplete(RHSSec, iErr); CHKERRQ(iErr)
      !!! VERY important! This is the equivalent of a ghost update
      Call SectionRealToVec(RHSSec, AppCtx%ScatterVect, SCATTER_FORWARD, AppCtx%RHSU, ierr); CHKERRQ(ierr)
      Call SectionRealDestroy(RHSSec, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
  End Subroutine RHSU_Assembly
   
   !----------------------------------------------------------------------------------------!      
   ! RHSAssemblyLocal (CM)  
   !----------------------------------------------------------------------------------------!      
   Subroutine RHSU_AssemblyLocal(iE, MatProp, AppCtx, RHSElem)

#if defined PB_2D
      Type(MatProp2D_Type)                         :: MatProp
#elif defined PB_3D
      Type(MatProp3D_Type)                         :: MatProp
#endif
   
      PetscInt                                     :: iE
      Type(AppCtx_Type)                            :: AppCtx
      PetscReal, Dimension(:), Pointer             :: RHSElem
      PetscReal, Dimension(:), Pointer             :: F, V
      PetscReal, Dimension(:), Pointer             :: Theta
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoFScal, NumDoFVect, NumGauss
      PetscInt, Dimension(:), Pointer              :: BCFlag
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      PetscReal                                    :: Theta_Elem, V_Elem
#if defined PB_2D
      Type (Vect2D)             				         :: TmpRHS
#elif defined PB_3D  
      Type (Vect3D)             				         :: TmpRHS    
#endif
      PetscLogDouble                               :: flops = 0
      
      
      Call PetscLogEventBegin(AppCtx%LogInfo%RHSAssemblyLocalU_Event, iErr); CHKERRQ(iErr)

      RHSElem    = 0.0_Kr
      NumDoFVect = Size(AppCtx%ElemVect(iE)%BF,1)
      NumDoFScal = Size(AppCtx%ElemScal(iE)%BF,1)
      NumGauss   = Size(AppCtx%ElemVect(iE)%BF,2)
      Allocate(BCFlag(NumDoFVect))
      Call MeshRestrictClosureInt(AppCtx%MeshTopology%mesh, AppCtx%BCUFlag, iE-1, NumDoFVect, BCFlag, iErr); CHKERRQ(ierr)
      Allocate(F(NumDoFVect))
      Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%F, iE-1, NumDoFVect, F, iErr); CHKERRQ(ierr)
      Allocate(Theta(NumDoFScal))
      Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%Theta, iE-1, NumDoFScal, Theta, iErr); CHKERRQ(ierr)
      Allocate(V(NumDoFScal))
      Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%V, iE-1, NumDoFScal, V, iErr); CHKERRQ(ierr)
      Do_iGauss: Do iGauss = 1, NumGauss
         TmpRHS     = 0.0_Kr
         Theta_Elem = 0.0_Kr
         V_Elem     = 0.0_Kr
         Do iDoF2 = 1, NumDoFVect
            TmpRHS = TmpRHS + AppCtx%ElemVect(iE)%BF(iDoF2, iGauss) * F(iDoF2)
         End Do
         Do iDoF2 = 1, NumDoFScal
            Theta_Elem = Theta_Elem + AppCtx%ElemScal(iE)%BF(iDoF2, iGauss) * Theta(iDoF2)
            V_Elem = V_Elem + AppCtx%ElemScal(iE)%BF(iDoF2, iGauss) * V(iDoF2)
            flops = flops + 4
         End Do
         Do iDoF1 = 1, NumDoFVect
            If (BCFlag(iDoF1) == 0) Then
               ! RHS terms due to forces
               RHSElem(iDoF1) = RHSElem(iDoF1) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * ( AppCtx%ElemVect(iE)%BF(iDoF1, iGauss) .DotP. TmpRHS ) 
               ! RHS terms due to inelastic strains
   			   RHSElem(iDoF1) = RHSElem(iDoF1) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * Theta_Elem * ((V_Elem**2) * (MatProp%Hookes_Law * AppCtx%ElemVect(iE)%GradS_BF(iDoF1, iGauss)) .DotP. MatProp%Therm_Exp)
   			   flops = flops + 5
            End If
         End Do
      End Do Do_iGauss
      Call PetscLogFlops(flops, iErr); CHKERRQ(iErr)
      DeAllocate(BCFlag)
      DeAllocate(F)
      DeAllocate(Theta)
      DeAllocate(V)
      Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
   End Subroutine RHSU_AssemblyLocal
   
   !----------------------------------------------------------------------------------------!      
   ! Solve U (CM)   
   !----------------------------------------------------------------------------------------!      
   
   Subroutine Solve_U(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      KSPConvergedReason                           :: reason
      PetscInt                                     :: KSPNumIter
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
      Type(Vec)                                    :: U_Vec
      
      Call PetscLogStagePush(AppCtx%LogInfo%KSPSolveU_Stage, iErr); CHKERRQ(iErr)
      Call MeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%U, U_Vec, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(AppCtx%U, AppCtx%ScatterVect, SCATTER_FORWARD, U_Vec, ierr); CHKERRQ(ierr)
      Call KSPSolve(AppCtx%KSPU, AppCtx%RHSU, U_Vec, iErr); CHKERRQ(iErr)
      !!! Solve and store the solution in AppCtx%RHS
      
      Call SectionRealToVec(AppCtx%U, AppCtx%ScatterVect, SCATTER_REVERSE, U_Vec, ierr); CHKERRQ(ierr)
      !!! Scatter the solution from (Vec) AppCtx%RHS to (SectionReal) AppCtx%U
      
      Call KSPGetConvergedReason(AppCtx%KSPU, reason, iErr); CHKERRQ(iErr)
      If ( reason > 0) Then
         Call KSPGetIterationNumber(AppCtx%KSPU, KSPNumIter, iErr); CHKERRQ(iErr)
         Write(IOBuffer, 100) KSPNumIter
      Else
         Write(IOBuffer, 101) reason
      End If
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      
      Call VecDestroy(U_Vec, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
100 Format('     KSP for U converged in ', I5, ' iterations \n'c)
101 Format('[ERROR] KSP for U diverged. KSPConvergedReason is ', I2, '\n'c)
   End Subroutine Solve_U
      
#if defined PB_2D
End Module m_VarFracQS_U2D
#elif defined PB_3D
End Module m_VarFracQS_U3D
#endif
