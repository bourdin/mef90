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
      PetscReal, Dimension(:), Pointer             :: U, U0
      Type(Vect2D)                                 :: DU
      PetscReal                                    :: NormKU
      PetscReal, Dimension(:), Pointer             :: Zero, One
      PetscReal                                    :: D_cond, ErrPhi
      PetscInt                                     :: iDoF, iGauss, iBlk, iE
      PetscInt                                     :: NumDoFVect, NumGauss
      PetscInt                                     :: numdec
      
      numdec = 0
      AppCtx%ErrPhi = 0
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Do_Elem_iE: Do iE = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
           NumDoFVect = Size(AppCtx%ElemVect(iE)%BF,1)
           NumGauss   = Size(AppCtx%ElemVect(iE)%BF,2)
           NormKU = 0.0_Kr
           Allocate(Zero(1))
           Zero = 0.0_Kr
           Allocate(One(1))
           One = 1.0_Kr
!           AppCtx%ErrPhi = 0
           Allocate(U(NumDoFVect))
           Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%U, iE-1, NumDoFVect, U, iErr); CHKERRQ(ierr)
           Allocate(U0(NumDoFVect))
           Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%U0, iE-1, NumDoFVect, U0, iErr); CHKERRQ(ierr)
           Do_iGauss: Do iGauss = 1, NumGauss
              DU = 0.0_Kr
              Do iDoF = 1, NumDoFVect
                   DU = DU + (U(iDof)-U0(iDof)) * AppCtx%ElemVect(iE)%BF(iE, iGauss)
                 ! Call PetscLogFlops(3 , iErr);CHKERRQ(iErr)
              End Do   
              NormKU = NormKU + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * ((AppCtx%MatProp( AppCtx%MeshTopology%Elem_Blk(iBlk)%ID )%K_Interface * DU ) .DotP. DU)
           End Do Do_iGauss
           D_cond =  NormKU - 2 * AppCtx%MatProp( AppCtx%MeshTopology%Elem_Blk(iBlk)%ID )%ToughnessD 
           If (D_cond >= 0.0_Kr) Then
              Call MeshUpdateClosure(AppCtx%MeshTopology%Mesh, AppCtx%Phi, iE-1, Zero, iErr); CHKERRQ(iErr)
              numdec = numdec+1
            Else
              Call MeshUpdateClosure(AppCtx%MeshTopology%Mesh, AppCtx%Phi, iE-1, One, iErr); CHKERRQ(iErr)
              ! SectionRealUpdate should work as well
!              AppCtx%ErrPhi = AppCtx%ErrPhi + 1  !If there is one (or more) update(s) the error is greater than 1!!!!!!!!!
           End If
         End Do Do_Elem_iE
       End Do Do_Elem_iBlk
       DeAllocate(U)
       DeAllocate(U0)
       DeAllocate(Zero) 
      
      Write(*,*) 'Number of debonded elements ', numdec 
!      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine Solve_Phi
End Module m_VarFracFilm_Phi
