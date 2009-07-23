Module m_VarFracFilm_U

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscmeshdef.h"
#include "finclude/petscviewerdef.h"

   Use m_VarFracFilm_Types
   Use m_MEF90
   Use m_VarFracFilm_Struct
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
      
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFracFilm_VertVar_DisplacementX)%Offset, AppCtx%TimeStep, UBC)
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
      
      PetscInt                                     :: iBlk, iE, iELoc, iErr
      PetscReal, Dimension(:,:), Pointer           :: MatElem
      
!      Call PetscLogStagePush(AppCtx%LogInfo%MatAssembly_Stage, iErr); CHKERRQ(iErr)
      
      Call MatZeroEntries(AppCtx%KU, iErr); CHKERRQ(iErr)
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Allocate(MatElem(AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim))
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)

            Call MatU_AssemblyLocal(iE, AppCtx%MatProp( AppCtx%MeshTopology%Elem_Blk(iBlk)%ID ), AppCtx, MatElem)
            Call assembleMatrix(AppCtx%KU, AppCtx%MeshTopology%mesh, AppCtx%U, iE-1, MatElem, ADD_VALUES, iErr); CHKERRQ(iErr)

         End Do Do_Elem_iE
         DeAllocate(MatElem)
      End Do Do_Elem_iBlk
      Call MatAssemblyBegin(AppCtx%KU, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (AppCtx%KU, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
!      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine MatU_Assembly

!----------------------------------------------------------------------------------------!      
! MatAssembly_U_Local (CM)
!----------------------------------------------------------------------------------------!      
   Subroutine MatU_AssemblyLocal(iE, MatProp, AppCtx, MatElem)
      Type(MatProp2D_Type)                         :: MatProp
      Type(AppCtx_Type)                            :: AppCtx
      PetscReal, Dimension(:,:), Pointer           :: MatElem      
      PetscInt                                     :: iE
   
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoFScal, NumDoFVect, NumGauss
      PetscInt, Dimension(:), Pointer              :: BCFlag
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      
      PetscReal, Dimension(:), Pointer             :: V
      PetscReal, Dimension(:), Pointer             :: Phi
      PetscReal                                    :: V_Elem
      
!      Call PetscLogEventBegin(AppCtx%LogInfo%MatAssemblyLocal_Event, iErr); CHKERRQ(iErr)

      MatElem  = 0.0_Kr
      NumDoFVect = Size(AppCtx%ElemVect(iE)%BF,1)
      NumDoFScal = Size(AppCtx%ElemScal(iE)%BF,1)
      NumGauss   = Size(AppCtx%ElemVect(iE)%BF,2)
      
      Allocate(V(NumDoFScal))
      Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%V, iE-1, NumDoFScal, V, iErr); CHKERRQ(ierr)
      Allocate(Phi(1))
      Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%Phi, iE-1, 1, Phi, iErr); CHKERRQ(ierr)
      Allocate(BCFlag(NumDoFVect))
      Call MeshRestrictClosureInt(AppCtx%MeshTopology%mesh, AppCtx%BCUFlag, iE-1, NumDoFVect, BCFlag, iErr); CHKERRQ(ierr)
      
      Do iGauss = 1, NumGauss
      !! Calculate V at the gauss points
         V_Elem = 0.0_Kr        
         Do iDoF1 = 1, NumDoFScal
            V_Elem = V_Elem + AppCtx%ElemScal(iE)%BF(iDoF1, iGauss)* V(iDoF1)
         End Do
    
     !! Assemble the element stiffness
         Do iDoF1 = 1, NumDoFVect
            If (BCFlag(iDoF1) == 0) Then
               Do iDoF2 = 1, NumDoFVect
               !Bulk contribution v^2(1/2*A(eps-eps0)^2)
                 MatElem(iDoF2, iDoF1) =  MatElem(iDoF2, iDoF1) +  AppCtx%ElemVect(iE)%Gauss_C(iGauss) * (V_Elem**2+AppCtx%VarFracFilmSchemeParam%KEpsilonV)*((MatProp%Hookes_Law * AppCtx%ElemVect(iE)%GradS_BF(iDoF1, iGauss)) .DotP. AppCtx%ElemVect(iE)%GradS_BF(iDoF2, iGauss))
               !Interface contribution phi*(1/2*K_interface*(u-u0)^2)
                   MatElem(iDoF2, iDoF1) =  MatElem(iDoF2, iDoF1) +  AppCtx%ElemVect(iE)%Gauss_C(iGauss) * ((Phi(1)+AppCtx%VarFracFilmSchemeParam%KEpsilonPhi)*(MatProp%K_interface * AppCtx%ElemVect(iE)%BF(iDoF1, iGauss)) .DotP. AppCtx%ElemVect(iE)%BF(iDoF2, iGauss))  
!			      Call PetscLogFlops(AppCtx%MeshTopology%num_dim * (AppCtx%MeshTopology%num_dim-1) +1 , iErr);CHKERRQ(iErr)
                  !!! Is that right?
               End Do
            End If
         End Do
      End Do
      DeAllocate(V)
      DeAllocate(Phi)
      DeAllocate(BCFlag)
!      Call PetscLogEventEnd(AppCtx%LogInfo%MatAssemblyLocal_Event, iErr); CHKERRQ(iErr)
   End Subroutine MatU_AssemblyLocal

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
!      Call PetscLogStagePush(AppCtx%LogInfo%RHSAssembly_Stage, iErr); CHKERRQ(iErr)
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
!      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
  End Subroutine RHSU_Assembly
   
   !----------------------------------------------------------------------------------------!      
   ! RHSAssemblyLocal (CM)  
   !----------------------------------------------------------------------------------------!      
   Subroutine RHSU_AssemblyLocal(iE, MatProp, AppCtx, RHSElem)

      Type(MatProp2D_Type)                         :: MatProp
   
      PetscInt                                     :: iE
      Type(AppCtx_Type)                            :: AppCtx
      PetscReal, Dimension(:), Pointer             :: RHSElem
      PetscReal, Dimension(:), Pointer             :: Phi
      PetscReal, Dimension(:), Pointer             :: U0
      PetscReal, Dimension(:), Pointer             :: V
      PetscReal, Dimension(:), Pointer             :: Theta
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoFScal, NumDoFVect, NumGauss
      PetscInt, Dimension(:), Pointer              :: BCFlag
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      PetscReal                                    :: Theta_Elem, V_Elem
      Type (Vect2D)                                :: U0_Elem

      RHSElem    = 0.0_Kr
      NumDoFVect = Size(AppCtx%ElemVect(iE)%BF,1)
      NumDoFScal = Size(AppCtx%ElemScal(iE)%BF,1)
      NumGauss   = Size(AppCtx%ElemVect(iE)%BF,2)
      Allocate(BCFlag(NumDoFVect))
      Call MeshRestrictClosureInt(AppCtx%MeshTopology%mesh, AppCtx%BCUFlag, iE-1, NumDoFVect, BCFlag, iErr); CHKERRQ(ierr)
      
      Allocate(U0(NumDoFVect))
      Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%U0, iE-1, NumDoFVect, U0, iErr); CHKERRQ(ierr)
      
      Allocate(Theta(NumDoFScal))
      Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%Theta, iE-1, NumDoFScal, Theta, iErr); CHKERRQ(ierr)
      
      Allocate(V(NumDoFScal))
      Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%V, iE-1, NumDoFScal, V, iErr); CHKERRQ(ierr)
      
      Allocate(Phi(1))
      Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%Phi, iE-1, 1, Phi, iErr); CHKERRQ(ierr)
      
      Do_iGauss: Do iGauss = 1, NumGauss
         U0_Elem     = 0.0_Kr
         Theta_Elem  = 0.0_Kr
         V_Elem      = 0.0_Kr
         Do iDoF2 = 1, NumDoFVect
            U0_Elem = U0_Elem + AppCtx%ElemVect(iE)%BF(iDoF2, iGauss) * U0(iDoF2)
         End Do
         Do iDoF2 = 1, NumDoFScal
            Theta_Elem = Theta_Elem + AppCtx%ElemScal(iE)%BF(iDoF2, iGauss) * Theta(iDoF2)
            V_Elem     = V_Elem     + AppCtx%ElemScal(iE)%BF(iDoF2, iGauss) * V(iDoF2)
         End Do
         Do iDoF1 = 1, NumDoFVect
            If (BCFlag(iDoF1) == 0) Then
               ! RHS terms due to U0
             !  RHSElem(iDoF1) = RHSElem(iDoF1) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) *  ( MatProp%K_interface * Phi(1) * (AppCtx%ElemVect(iE)%BF(iDoF1, iGauss)) .DotP. U0_Elem ) 
               ! RHS terms due to eps0
   			   RHSElem(iDoF1) = RHSElem(iDoF1) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * Theta_Elem * ((V_Elem**2) * (MatProp%Hookes_Law * AppCtx%ElemVect(iE)%GradS_BF(iDoF1, iGauss)) .DotP. MatProp%Therm_Exp)
!               Call PetscLogFlops(3 , iErr);CHKERRQ(iErr)
            End If
         End Do
      End Do Do_iGauss
      DeAllocate(BCFlag)
      DeAllocate(U0)
      DeAllocate(Theta)
      DeAllocate(V)
!      Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyLocal_Event, iErr); CHKERRQ(iErr)
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
      
!      Call PetscLogStagePush(AppCtx%LogInfo%KSPSolve_Stage, iErr); CHKERRQ(iErr)
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
!      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
100 Format('     KSP for U converged in ', I5, ' iterations \n')
101 Format('[ERROR] KSP for U diverged. KSPConvergedReason is ', I2, '\n')
   End Subroutine Solve_U
      
End Module m_VarFracFilm_U