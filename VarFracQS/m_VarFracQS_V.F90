#if defined PB_2D
Module m_VarFracQS_V2D
#elif defined PB_3D
Module m_VarFracQS_V3D
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

!----------------------------------------------------------------------------------------!      
! MatAssembly V (CM)  
!----------------------------------------------------------------------------------------!      
   Subroutine MatV_Assembly(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iBlk, iE, iELoc, iErr
      PetscReal, Dimension(:,:), Pointer           :: MatElem
      
!      Call PetscLogStagePush(AppCtx%LogInfo%MatAssembly_Stage, iErr); CHKERRQ(iErr)

      Call MatZeroEntries(AppCtx%KV, iErr); CHKERRQ(iErr)
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Allocate(MatElem(AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF))
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            !--------------------------------------------------------
            ! The only line to change as a function of the problem.
            Call MatV_AssemblyLocal(iE, AppCtx%MatProp( AppCtx%MeshTopology%Elem_Blk(iBlk)%ID ), AppCtx, MatElem)
            !--------------------------------------------------------
            Call assembleMatrix(AppCtx%KV, AppCtx%MeshTopology%mesh, AppCtx%V, iE-1, MatElem, ADD_VALUES, iErr); CHKERRQ(iErr)
         End Do Do_Elem_iE
         DeAllocate(MatElem)
      End Do Do_Elem_iBlk
      Call MatAssemblyBegin(AppCtx%KV, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (AppCtx%KV, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)

!      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
End Subroutine MatV_Assembly

!----------------------------------------------------------------------------------------!      
! MatAssembly_V_Local (CM)
!----------------------------------------------------------------------------------------!      
   Subroutine MatV_AssemblyLocal(iE, MatProp, AppCtx, MatElem)
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
      
      PetscReal, Dimension(:), Pointer             :: U, Theta
      PetscReal                                    :: Theta_Elem, ElasEnDens_Elem
      PetscReal                                    :: C2_V, C2_GradV
#if defined PB_2D
      Type(MatS2D)                                 :: Strain_Elem, Effective_Strain_Elem
#elif defined PB_3D
      Type(MatS3D)                                 :: Strain_Elem, Effective_Strain_Elem
#endif      
!      Call PetscLogEventBegin(AppCtx%LogInfo%MatAssemblyLocal_Event, iErr); CHKERRQ(iErr),


      MatElem  = 0.0_Kr
      NumDoFVect = Size(AppCtx%ElemVect(iE)%BF,1)
      NumDoFScal = Size(AppCtx%ElemScal(iE)%BF,1)
      NumGauss   = Size(AppCtx%ElemVect(iE)%BF,2)
      
      !! Get the local nodal values of U, Theta, and BCs
      Allocate(U(NumDoFVect))
      Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%U, iE-1, NumDoFVect, U, iErr); CHKERRQ(ierr)
      Allocate(Theta(NumDoFScal))
      Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%Theta, iE-1, NumDoFScal, Theta, iErr); CHKERRQ(ierr)
      Allocate(BCFlag(NumDoFScal))
      Call MeshRestrictClosureInt(AppCtx%MeshTopology%mesh, AppCtx%BCFlagV, iE-1, NumDoFScal, BCFlag, iErr); CHKERRQ(ierr)
      
      Do iGauss = 1, NumGauss
      !! Calculate the strain at the gauss point
         Strain_Elem = 0.0_Kr
         Theta_Elem  = 0.0_Kr
         Do iDoF1 = 1, NumDoFVect
            Strain_Elem = Strain_Elem + (AppCtx%ElemVect(iE)%GradS_BF(iDoF1, iGauss) * U(iDoF1))
!             Call PetscLogFlops(3*AppCtx%MeshTopology%Num_Dim+2, iErr)
         End Do
      !! Calculate the temperature at the gauss point
         Do iDoF1 = 1, NumDoFScal
            Theta_Elem = Theta_Elem + AppCtx%ElemScal(iE)%BF(iDoF1, iGauss) * Theta(iDoF1)
         End Do
      !! Calculate the Effective Strain at the gauss point
         Effective_Strain_Elem  =  Strain_Elem - (Theta_Elem * MatProp%Therm_Exp)   
      !! Calculate the coefficients of the terms v^2 (C2_V) et GradV*GradV (C2_GradV) of the energy functional
         C2_V     = 0.5_Kr / AppCtx%VarFracSchemeParam%Epsilon * MatProp%Toughness + ((MatProp%Hookes_Law * Effective_Strain_Elem) .DotP. Effective_Strain_Elem)
         C2_GradV = 2.0_Kr * AppCtx%VarFracSchemeParam%Epsilon * MatProp%Toughness
      !! Assemble the element stiffness
         Do iDoF1 = 1, NumDoFScal
            If (BCFlag(iDoF1) == 0) Then
               Do iDoF2 = 1, NumDoFScal
               !! Terms in V^2
                  MatElem(iDoF2, iDoF1) =  MatElem(iDoF2, iDoF1) + AppCtx%ElemScal(iE)%Gauss_C(iGauss) * C2_V * AppCtx%ElemScal(iE)%BF(iDoF1, iGauss) * AppCtx%ElemScal(iE)%BF(iDoF2, iGauss)
               !! Terms in GradV^2
                  MatElem(iDoF2, iDoF1) =  MatElem(iDoF2, iDoF1) + AppCtx%ElemScal(iE)%Gauss_C(iGauss) * C2_GradV * (AppCtx%ElemScal(iE)%Grad_BF(iDoF1, iGauss) .DotP. AppCtx%ElemScal(iE)%Grad_BF(iDoF2, iGauss) )               
!              Call PetscLogFlops(AppCtx%MeshTopology%num_dim * (AppCtx%MeshTopology%num_dim-1) +1 , iErr);CHKERRQ(iErr)
                  !!! Is that right?
               End Do
            End If
         End Do
      End Do
      DeAllocate(U)      
      DeAllocate(Theta)
      DeAllocate(BCFlag)
!      Call PetscLogEventEnd(AppCtx%LogInfo%MatAssemblyLocal_Event, iErr); CHKERRQ(iErr)
   End Subroutine MatV_AssemblyLocal

!----------------------------------------------------------------------------------------!      
! RHSAssembly (CM)  
!----------------------------------------------------------------------------------------!      
   Subroutine RHSV_Assembly(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      Type(SectionReal)                            :: RHSSec
      PetscInt                                     :: iBlk, iE, iELoc
      PetscReal, Dimension(:), Pointer             :: RHSElem
      PetscInt                                     :: NumDoFPerVertex 
      
      NumDoFPerVertex = 1
      
      !!! Hopefully one day we will use assemble Vector instead of going through a section
!      Call PetscLogStagePush(AppCtx%LogInfo%RHSAssembly_Stage, iErr); CHKERRQ(iErr)
      Call VecSet(AppCtx%RHSV, 0.0_Kr, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'RHSSec', NumDoFPerVertex, RHSSec, iErr); CHKERRQ(iErr)
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Allocate(RHSElem(AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF))
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            Call RHSV_AssemblyLocal(iE, AppCtx%MatProp( AppCtx%MeshTopology%Elem_Blk(iBlk)%ID ), AppCtx, RHSElem)
            Write(MEF90_MyRank+100, *) iE, RHSElem
            Call MeshUpdateAddClosure(AppCtx%MeshTopology%Mesh, RHSSec, iE-1, RHSElem, iErr); CHKERRQ(iErr)
         End Do Do_Elem_iE
         DeAllocate(RHSElem)
      End Do Do_Elem_iBlk
      Call SectionRealComplete(RHSSec, iErr); CHKERRQ(iErr)
      !!! VERY important! This is the equivalent of a ghost update
      Call SectionRealToVec(RHSSec, AppCtx%ScatterScal, SCATTER_FORWARD, AppCtx%RHSV, ierr); CHKERRQ(ierr)
      Call SectionRealDestroy(RHSSec, iErr); CHKERRQ(iErr)
!      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
  End Subroutine RHSV_Assembly
   
!----------------------------------------------------------------------------------------!      
! RHSAssemblyLocal (CM)  
!----------------------------------------------------------------------------------------!      
   Subroutine RHSV_AssemblyLocal(iE, MatProp, AppCtx, RHSElem)

#if defined PB_2D
      Type(MatProp2D_Type)                         :: MatProp
#elif defined PB_3D
      Type(MatProp3D_Type)                         :: MatProp
#endif
   
      PetscInt                                     :: iE
      Type(AppCtx_Type)                            :: AppCtx
      PetscReal, Dimension(:), Pointer             :: RHSElem
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoFScal, NumGauss
      PetscInt, Dimension(:), Pointer              :: BCFlag
      PetscInt                                     :: iDoF1, iGauss
      PetscReal                                    :: C1_V

      RHSElem    = 0.0_Kr
      NumDoFScal = Size(AppCtx%ElemScal(iE)%BF,1)
      NumGauss   = Size(AppCtx%ElemVect(iE)%BF,2)
      Allocate(BCFlag(NumDoFScal))
      Call MeshRestrictClosureInt(AppCtx%MeshTopology%mesh, AppCtx%BCFlagV, iE-1, NumDoFScal, BCFlag, iErr); CHKERRQ(ierr)
      ! Calculate the coefficient of the term in V (C1_V) of the energy functional
      C1_V =  0.5_Kr / AppCtx%VarFracSchemeParam%Epsilon * MatProp%Toughness 
      Do_iGauss: Do iGauss = 1, NumGauss
          Do iDoF1 = 1, NumDoFScal
            If (BCFlag(iDoF1) == 0) Then
               ! RHS terms due to forces
               RHSElem(iDoF1) = RHSElem(iDoF1) + C1_V * AppCtx%ElemScal(iE)%Gauss_C(iGauss) *  AppCtx%ElemScal(iE)%BF(iDoF1, iGauss) 
!               Call PetscLogFlops(3 , iErr);CHKERRQ(iErr)
            End If
         End Do
      End Do Do_iGauss
      DeAllocate(BCFlag)
      Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyLocal_Event, iErr); CHKERRQ(iErr)
   End Subroutine RHSV_AssemblyLocal
      

   !----------------------------------------------------------------------------------------!      
   ! Solve V (CM)   
   !----------------------------------------------------------------------------------------!      
   
   Subroutine Solve_V(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      KSPConvergedReason                           :: reason
      PetscInt                                     :: KSPNumIter
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
      Type(Vec)                                    :: V_Vec, V_Old
      PetscReal                                    :: VMin, VMax
      
!      Call PetscLogStagePush(AppCtx%LogInfo%KSPSolve_Stage, iErr); CHKERRQ(iErr)
  
      Call MeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%V, V_Vec, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(AppCtx%V, AppCtx%ScatterScal, SCATTER_FORWARD, V_Vec, ierr); CHKERRQ(ierr)

      Call VecDuplicate(V_Vec, V_Old, iErr); CHKERRQ(iErr)
      Call VecCopy(V_Vec, V_Old, iErr); CHKERRQ(iErr)

      Call KSPSolve(AppCtx%KSPV, AppCtx%RHSV, V_Vec, iErr); CHKERRQ(iErr)
      !!! Solve and store the solution in AppCtx%RHS
      
      Call SectionRealToVec(AppCtx%V, AppCtx%ScatterScal, SCATTER_REVERSE, V_Vec, ierr); CHKERRQ(ierr)
       
      !!! Scatter the solution from (Vec) AppCtx%RHS to (SectionReal) AppCtx%V
      
      Call KSPGetConvergedReason(AppCtx%KSPV, reason, iErr); CHKERRQ(iErr)
      If ( reason > 0) Then
         Call KSPGetIterationNumber(AppCtx%KSPV, KSPNumIter, iErr); CHKERRQ(iErr)
         Write(IOBuffer, 100) KSPNumIter
      Else
         Write(IOBuffer, 101) reason
      End If
      
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

      Call VecMin(V_Vec, PETSC_NULL_INTEGER, VMin, iErr); CHKERRQ(iErr)
      Call VecMax(V_Vec, PETSC_NULL_INTEGER, VMax, iErr); CHKERRQ(iErr)
      Write(IOBuffer, 700) VMin, VMax
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

       Call VecAxPy(V_Old, -1.0_Kr, V_Vec, iErr)
       Call VecNorm(V_Old, NORM_INFINITY, AppCtx%ErrV, iErr)

      Write(IOBuffer, 800) AppCtx%ErrV
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

      
      Call VecDestroy(V_Old, iErr); CHKERRQ(iErr)
      Call VecDestroy(V_Vec, iErr); CHKERRQ(iErr)
!      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
100 Format('     KSP for V converged in ', I5, ' iterations \n'c)
101 Format('[ERROR] KSP for V diverged. KSPConvergedReason is ', I2, '\n'c)
700 Format('     VMin / Max:   ', T24, 2(ES12.5, '  '), '\n'c)
800 Format('     Max change V: ', T24, ES12.5, '\n'c)
   End Subroutine Solve_V

   
   
#if defined PB_2D
End Module m_VarFracQS_V2D
#elif defined PB_3D
End Module m_VarFracQS_V3D
#endif
