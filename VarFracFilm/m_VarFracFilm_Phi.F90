Module m_VarFracFilm_Phi

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
      

   !----------------------------------------------------------------------------------------!      
   ! Solve Phi (CM)   
   !----------------------------------------------------------------------------------------!      
   
   Subroutine Solve_Phi(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      KSPConvergedReason                           :: reason
      PetscInt                                     :: KSPNumIter
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
      Type(Vec)                                    :: Phi_Vec, Phi_Old
      PetscReal                                    :: PhiMin, PhiMax
      PetscReal, Dimension(:), Pointer             :: U, U0, DU
      PetscReal                                    :: NormKU
      PetscReal, Dimension(:), Pointer             :: Phi
      PetscReal                                    :: D_cond, ErrPhi
      PetscInt                                     :: iDoF, iGauss, iBlk, iE
      PetscInt                                     :: NumDoFVect, NumGauss
      PetscReal                                    :: NumberOfUpdates

!      Call PetscLogStagePush(AppCtx%LogInfo%KSPSolve_Stage, iErr); CHKERRQ(iErr)
!      !! Create the vector Phi_Vec and duplicate in Phi_Old
!      Call MeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%Phi, Phi_Vec, iErr); CHKERRQ(iErr)
!      Call SectionRealToVec(AppCtx%Phi, AppCtx%ScatterScal, SCATTER_FORWARD, Phi_Vec, ierr); CHKERRQ(ierr)
!      Call VecDuplicate(Phi_Vec, Phi_Old, iErr); CHKERRQ(iErr)
!      Call VecCopy(Phi_Vec, Phi_Old, iErr); CHKERRQ(iErr)     
      NumberOfUpdates = 0.0_Kr

      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Do_Elem_iE: Do iE = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
           NumDoFVect = Size(AppCtx%ElemVect(iE)%BF,1)
           NumGauss   = Size(AppCtx%ElemVect(iE)%BF,2)
           Allocate(U0(NumDoFVect))
           NormKU = 0.0_Kr
           Allocate(Phi(1))
           Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%Phi, iE-1, 1, Phi, iErr); CHKERRQ(ierr)
           Allocate(U(NumDoFVect))
           Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%U, iE-1, NumDoFVect, U, iErr); CHKERRQ(ierr)
           Allocate(U0(NumDoFVect))
           Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%U0, iE-1, NumDoFVect, U0, iErr); CHKERRQ(ierr)
           Do_iGauss: Do iGauss = 1, NumGauss
              DU = 0.0_Kr
              Do iDoF = 1, NumDoFVect
                   DU = DU + U(iDof)-U0(iDof)
                 ! Call PetscLogFlops(3 , iErr);CHKERRQ(iErr)
              End Do   
             ! NormKU = NormKU + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * ((AppCtx%MatProp( AppCtx%MeshTopology%Elem_Blk(iBlk)%ID ) * DU) .DotP. DU)
           End Do Do_iGauss
           D_cond =  NormKU - 2 * AppCtx%MatProp( AppCtx%MeshTopology%Elem_Blk(iBlk)%ID )%ToughnessD 
           If (D_cond .GT. 0) Then
           Call MeshUpdateClosure(AppCtx%MeshTopology%Mesh, AppCtx%Phi, iE-1, 0)
           ! SectionRealUpdate should work as well
           NumberOfUpdates = NumberOfUpdates + 1  !If there is one (or more) update(s) the error is greater than 1!!!!!!!!!
           End If
         End Do Do_Elem_iE
       End Do Do_Elem_iBlk
!      Call SectionRealToVec(AppCtx%U, AppCtx%ScatterPhi, SCATTER_FORWARD, Phi, ierr); CHKERRQ(ierr)       
       DeAllocate(U)
       DeAllocate(U0)
       DeAllocate(Phi) 
       AppCtx%ErrPhi=NumberOfUpdates   
      
!      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
    ! I do not how to write with Petsc, sigh !!!!
! 100 Format('There were ', I5, ' new delaminated elements \n'c)
   End Subroutine Solve_Phi

   
End Module m_VarFracFilm_Phi
