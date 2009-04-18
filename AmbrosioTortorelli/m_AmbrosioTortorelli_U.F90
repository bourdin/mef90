#if defined PB_2D
Module m_AmbrosioTortorelli_U2D
#elif defined PB_3D
Module m_AmbrosioTortorelli_U3D
#endif

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscmeshdef.h"
#include "finclude/petscviewerdef.h"

#if defined PB_2D
   Use m_AmbrosioTortorelli_Types2D
#elif defined PB_3D
   Use m_AmbrosioTortorelli_Types3D
#endif   Use m_MEF90
   Use m_RuptStruct
   Use petsc
   Use petscvec
   Use petscmat
   Use petscksp
   Use petscmesh

   Implicit NONE   
   


   !----------------------------------------------------------------------------------------!      
   ! MatAssembly (CM)  
   !----------------------------------------------------------------------------------------!      
   Subroutine MatAssembly_U(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iBlk, iE, iELoc, iErr
      PetscReal, Dimension(:,:), Pointer           :: MatElem
      
!      Call PetscLogStagePush(AppCtx%LogInfo%MatAssembly_Stage, iErr); CHKERRQ(iErr)
      
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Allocate(MatElem(AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim))
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            !--------------------------------------------------------
            Call MatAssembly_U_Local(iE, AppCtx%MatProp( AppCtx%MeshTopology%Elem_Blk(iBlk)%ID ), AppCtx, MatElem)
            !--------------------------------------------------------
            Call assembleMatrix(AppCtx%KU, AppCtx%MeshTopology%mesh, AppCtx%U, iE-1, MatElem, ADD_VALUES, iErr); CHKERRQ(iErr)
         End Do Do_Elem_iE
         DeAllocate(MatElem)
      End Do Do_Elem_iBlk
      Call MatAssemblyBegin(AppCtx%KU, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (AppCtx%KU, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
!      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine MatAssembly

Subroutine MatAssemblyU(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iBlk, iE, iELoc, iErr
      PetscReal, Dimension(:,:), Pointer           :: MatElem
End Subroutine MatAssembly

!----------------------------------------------------------------------------------------!      
! MatAssembly_U_Local (CM)
!----------------------------------------------------------------------------------------!      
   Subroutine MatAssembly_U_Local(iE, MatProp, AppCtx, MatElem)
#if defined PB_2D
      Type(MatProp2D_Type)                         :: MatProp
#elif defined PB_3D
      Type(MatProp3D_Type)                         :: MatProp
#endif
      Type(AppCtx_Type)                            :: AppCtx
      PetscReal, Dimension(:,:), Pointer           :: MatElem      
      PetscInt                                     :: iE
   
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoFScal, NumDoFVect, NumGauss
      PetscInt, Dimension(:), Pointer              :: BCFlag
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      
      PetscReal, Dimension(:), Pointer             :: V
      PetscReal                                    :: V_Elem
      
!      Call PetscLogEventBegin(AppCtx%LogInfo%MatAssemblyLocal_Event, iErr); CHKERRQ(iErr)

      MatElem  = 0.0_Kr
      NumDoFVec  = Size(AppCtx%ElemVect(iE)%BF,1)
      NumDoFScal = Size(AppCtx%ElemScal(iE)%BF,1)
      NumGauss   = Size(AppCtx%ElemVect(iE)%BF,2)
      
      Allocate(V(NumDoFScal))
      Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%V, iE-1, NumDoFScal, V, iErr); CHKERRQ(ierr)

      Allocate(BCFlag(NumDoFVect))
      Call MeshRestrictClosureInt(AppCtx%MeshTopology%mesh, AppCtx%BCFlagU, iE-1, NumDoF, BCFlag, iErr); CHKERRQ(ierr)
      
      Do iGauss = 1, NumGauss
      !! Calculate V at the gauss point
         V_Elem = 0.0_Kr        
         Do iDoF1 = 1, NumDoFScal
            V_Elem = V_Elem + AppCtx%ElemScal(iE)%BF(iDoF1, iGauss)* V(iDoF1)
         End Do
      !! Assemble the element stiffness
         Do iDoF1 = 1, NumDoFVect
            If (BCFlag(iDoF1) == 0) Then
               Do iDoF2 = 1, NumDoFVect
                  MatElem(iDoF2, iDoF1) =  MatElem(iDoF2, iDoF1) +  AppCtx%ElemVect(iE)%Gauss_C(iGauss) * ((V_Elem**2)*(MatProp%Hookes_Law * AppCtx%ElemVect(iE)%GradS_BF(iDoF1, iGauss)) .DotP. AppCtx%ElemVect(iE)%GradS_BF(iDoF2, iGauss))
!			      Call PetscLogFlops(AppCtx%MeshTopology%num_dim * (AppCtx%MeshTopology%num_dim-1) +1 , iErr);CHKERRQ(iErr)
                  !!! Is that right?
               End Do
            End If
         End Do
      End Do
      DeAllocate(V)
      DeAllocate(BCFlag)
!      Call PetscLogEventEnd(AppCtx%LogInfo%MatAssemblyLocal_Event, iErr); CHKERRQ(iErr)
   End Subroutine MatAssembly_U_Local

   !----------------------------------------------------------------------------------------!      
   ! RHSAssembly (CM)  
   !----------------------------------------------------------------------------------------!      
   Subroutine RHS_U_Assembly(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      Type(SectionReal)                            :: RHSSec
      PetscInt                                     :: iBlk, iE, iELoc
      PetscReal, Dimension(:), Pointer             :: RHSElem
      PetscInt                                     :: NumDoFPerVertex 
      
      NumDoFPerVertex = AppCtx%MeshTopology%Num_Dim
      
      !!! Hopefully one day we will use assemble Vector instead of going through a section
!      Call PetscLogStagePush(AppCtx%LogInfo%RHSAssembly_Stage, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, NumDoFPerVertex, RHSSec, iErr); CHKERRQ(iErr)
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Allocate(RHSElem(AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim))
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            !--------------------------------------------------------
            Call RHS_U_AssemblyLocal(iE, AppCtx%MatProp( AppCtx%MeshTopology%Elem_Blk(iBlk)%ID ),AppCtx, RHSElem)
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
  End Subroutine RHS_U_Assembly
   
   !----------------------------------------------------------------------------------------!      
   ! RHSAssemblyLocal (CM)  
   !----------------------------------------------------------------------------------------!      
   Subroutine RHS_U_AssemblyLocal(iE, MatProp, AppCtx, RHSElem)

#if defined PB_2D
      Type(MatProp2D_Type)                         :: MatProp
#elif defined PB_3D
      Type(MatProp3D_Type)                         :: MatProp
#endif
   
      PetscInt                                     :: iE
      Type(AppCtx_Type)                            :: AppCtx
      PetscReal, Dimension(:), Pointer             :: RHSElem
      PetscReal, Dimension(:), Pointer             :: F
      PetscReal, Dimension(:), Pointer             :: Theta
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoFScal, NumDoFVect, NumGauss
      PetscInt, Dimension(:), Pointer              :: BCFlag
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      PetscReal                                    :: Theta_Elem
#if defined PB_2D
      Type (Vect2D)             				         :: TmpRHS
#elif defined PB_3D  
      Type (Vect3D)             				         :: TmpRHS    
#endif

      RHSElem    = 0.0_Kr
      NumDoFVect = Size(AppCtx%ElemVect(iE)%BF,1)
      NumDoFScal = Size(AppCtx%ElemScal(iE)%BF,1)
      NumGauss   = Size(AppCtx%ElemVect(iE)%BF,2)
      Allocate(BCFlag(NumDoFVect))
      Call MeshRestrictClosureInt(AppCtx%MeshTopology%mesh, AppCtx%BCFlagU, iE-1, NumDoFVect, BCFlag, iErr); CHKERRQ(ierr)
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
            Theta_Elem = Theta_Elem + AppCtx%ElemScal(iE)%BF(iDoF2, iGauss)* Theta(iDoF2)
            V_Elem = V_Elem + AppCtx%ElemVect(iE)%BF(iDoF2, iGauss) * F(iDoF2)
         End Do
         Do iDoF1 = 1, NumDoFVect
            If (BCFlag(iDoF1) == 0) Then
               ! RHS terms due to forces
               RHSElem(iDoF1) = RHSElem(iDoF1) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * ( AppCtx%ElemVect(iE)%BF(iDoF1, iGauss) .DotP. TmpRHS ) 
               ! RHS terms due to inelastic strains
   			   RHSElem(iDoF1) = RHSElem(iDoF1) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * Theta_Elem * ((V_Elem**2)*(MatProp%Hookes_Law * AppCtx%ElemVect(iE)%GradS_BF(iDoF1, iGauss)) .DotP. MatProp%Therm_Exp)
!               Call PetscLogFlops(3 , iErr);CHKERRQ(iErr)
            End If
         End Do
      End Do Do_iGauss
      DeAllocate(BCFlag)
      DeAllocate(F)
      DeAllocate(Theta)
      DeAllocate(V)
!      Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyLocal_Event, iErr); CHKERRQ(iErr)
   End Subroutine RHS_U_AssemblyLocal
   
   
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
      
      Call KSPGetIterationNumber(AppCtx%KSPU, KSPNumIter, iErr); CHKERRQ(iErr)
      Call KSPGetConvergedReason(AppCtx%KSPU, reason, iErr); CHKERRQ(iErr)
      Write(IOBuffer, 100) KSPNumIter, reason
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      
      Call VecDestroy(U_Vec, iErr); CHKERRQ(iErr)
!      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
100 Format('KSP converged in ', I5, ' iterations. KSPConvergedReason is ', I2, '\n'c)
   End Subroutine Solve_U
      
#if defined PB_2D
End Module m_ m_AmbrosioTortorelli_U2D
#elif defined PB_3D
End Module m_ m_AmbrosioTortorelli_U3D
#endif
