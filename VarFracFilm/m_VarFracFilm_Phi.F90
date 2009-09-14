Module m_VarFracFilm_Phi

#include "finclude/petscdef.h"
!#include "finclude/petscvec.h"
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

  Subroutine Init_TS_PHI(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscReal, Dimension(:), Pointer             :: PHI_Ptr, Zero
      PetscInt                                     :: i, iErr
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer      
      !!! WTF? 
      !!! Why not using restrict / update instead of the mess?
            
          If (AppCtx%AppParam%verbose) Then
            Write(IOBuffer, *) "Initializing Phi with Init_Phi_One:", Init_PHI_Prev, "\n"c
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         End If
         Allocate(PHI_Ptr(1))
         Allocate(Zero(1))
         PHI_Ptr(1) = 1.0_Kr
         Zero(1) = 0.0_Kr               
         Do i = 1, AppCtx%MeshTopology%Num_Verts       
           Call MeshUpdateClosure(AppCtx%MeshTopology%Mesh, AppCtx%PHI, AppCtx%MeshTopology%Num_Elems + i-1, PHI_Ptr, iErr); CHKERRQ(iErr)
           Call MeshUpdateClosure(AppCtx%MeshTopology%Mesh,AppCtx%RHSPHI, AppCtx%MeshTopology%Num_Elems + i-1, Zero, iErr); CHKERRQ(iErr)
         End Do
         DeAllocate(PHI_Ptr)
         DeAllocate(Zero)
        End Subroutine Init_TS_PHI
        


!----------------------------------------------------------------------------------------!      
! RHSAssembly (CM)  
!----------------------------------------------------------------------------------------!      
   Subroutine RHSPHI_Assembly(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
!      Type(SectionReal)                            :: RHSPhiSec
      PetscInt                                     :: iBlk, iE, iELoc
      PetscReal, Dimension(:), Pointer             :: RHSElem
      PetscInt                                     :: NumDoFPerVertex 
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer

      NumDoFPerVertex = 1
      
      !!! Hopefully one day we will use assemble Vector instead of going through a section
!      Call PetscLogStagePush(AppCtx%LogInfo%RHSAssembly_Stage, iErr); CHKERRQ(iErr)
      !Call VecSet(AppCtx%RHSPHI, 0.0_Kr, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'RHSPhiSec', NumDoFPerVertex, AppCtx%RHSPHI, iErr); CHKERRQ(iErr)
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Allocate(RHSElem(AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF))
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            Call RHSPHI_AssemblyLocal(iE, AppCtx%MatProp( AppCtx%MeshTopology%Elem_Blk(iBlk)%ID ), AppCtx, RHSElem)
            Call MeshUpdateAddClosure(AppCtx%MeshTopology%Mesh, AppCtx%RHSPHI, iE-1, RHSElem, iErr); CHKERRQ(iErr)
         End Do Do_Elem_iE
         DeAllocate(RHSElem)
      End Do Do_Elem_iBlk
      Call SectionRealComplete(AppCtx%RHSPHI, iErr); CHKERRQ(iErr)
      !!! VERY important! This is the equivalent of a ghost update
      !Call SectionRealToVec(RHSPhiSec, AppCtx%ScatterScal, SCATTER_FORWARD, AppCtx%RHSPHI, ierr); CHKERRQ(ierr)
      !!! This update the vector RHSPHI scattering the values from RHSPhiSec
      !Call SectionRealDestroy(RHSPhiSec, iErr); CHKERRQ(iErr)
!      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
  
    
  End Subroutine RHSPHI_Assembly
   
!----------------------------------------------------------------------------------------!      
! RHSPHIAssemblyLocal (CM)  
!----------------------------------------------------------------------------------------!      
   Subroutine RHSPHI_AssemblyLocal(iE, MatProp, AppCtx, RHSPHIElem)

      Type(MatProp2D_Type)                         :: MatProp
   
      PetscInt                                     :: iE
      Type(AppCtx_Type)                            :: AppCtx
      PetscReal, Dimension(:), Pointer             :: RHSPHIElem
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoFScal, NumGauss,NumDoFVect
      PetscReal, Dimension(:), Pointer             :: U0 
      PetscReal, Dimension(:), Pointer             :: U
!      PetscInt, Dimension(:), Pointer              :: BCFlag
      PetscInt                                     :: iDoF1, iGauss, iDoF2
!      PetscReal                                    :: Vol
      Type (Vect2D)                                :: DU_Elem
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer



      NumDoFScal = Size(AppCtx%ElemScal(iE)%BF,1)
      NumDoFVect = Size(AppCtx%ElemVect(iE)%BF,1)
      NumGauss   = Size(AppCtx%ElemVect(iE)%BF,2)
      !! Get the local nodal values of U, U0

      RHSPHIElem    = 0.0_Kr

      Allocate(U0(NumDoFVect))
      Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%U0, iE-1, NumDoFVect, U0, iErr); CHKERRQ(ierr)
      Allocate(U(NumDoFVect))
      Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%U, iE-1, NumDoFVect, U, iErr); CHKERRQ(ierr)
!      Allocate(BCFlag(NumDoFScal))
!      Call MeshRestrictClosureInt(AppCtx%MeshTopology%mesh, AppCtx%BCVFlag, iE-1, NumDoFScal, BCFlag, iErr); CHKERRQ(ierr)
! 
      Do_iGauss: Do iGauss = 1, NumGauss
          DU_Elem     = 0.0_Kr
!        Vol = 0.0_Kr
          Do iDoF2 = 1, NumDoFVect
            DU_Elem = DU_Elem + AppCtx%ElemVect(iE)%BF(iDoF2, iGauss) * (U(iDoF2)-U0(iDoF2))
         End Do
          Do iDoF1 = 1, NumDoFScal
!           If (BCFlag(iDoF1) == 0) Then
               ! RHS terms due to forces 
               ! Gauss_C(iGauss) for Scalar or vectorial elements ????
               RHSPHIElem(iDoF1) = RHSPHIElem(iDoF1) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * AppCtx%ElemScal(iE)%BF(iDoF1, iGauss) * 0.5 * ((MatProp%K_Interface  * DU_Elem ) .DotP. DU_Elem) 
               RHSPHIElem(iDoF1) = RHSPHIElem(iDoF1) - AppCtx%ElemVect(iE)%Gauss_C(iGauss) * AppCtx%ElemScal(iE)%BF(iDoF1, iGauss) * MatProp%ToughnessD
!               Call PetscLogFlops(3 , iErr);CHKERRQ(iErr)
!            End If
         End Do
      End Do Do_iGauss
!      DeAllocate(BCFlag)
       DeAllocate(U0)
       DeAllocate(U)
  
     
       Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyLocal_Event, iErr); CHKERRQ(iErr)

   End Subroutine RHSPHI_AssemblyLocal


   !----------------------------------------------------------------------------------------!      
   ! Solve Phi (CM)   
   !----------------------------------------------------------------------------------------!      
   
   Subroutine Solve_Phi(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr
      PetscInt                                     :: i
      PetscReal, Dimension(:), Pointer             :: RHSPHI_Ptr, PHI0, PHI1
 !     Type(Vec)                                    :: Phi_Vec, Phi_Old
      PetscInt                                     :: MyNumDec, NumDec
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer

          

         MyNumDec = 0  
         
         Allocate(PHI0(1))
         PHI0(1) = 0.0_Kr         
         Allocate(PHI1(1))
         PHI1(1) = 1.0_Kr
         Allocate( RHSPhi_Ptr(1))
         Do i = 1, AppCtx%MeshTopology%Num_Verts
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%RHSPHI, AppCtx%MeshTopology%Num_Elems + i-1, AppCtx%MeshTopology%Num_Dim, RHSPHI_Ptr, iErr); CHKERRQ(ierr)      
            If (RHSPHI_Ptr(1) > 0) Then
               Call MeshUpdateClosure(AppCtx%MeshTopology%mesh, AppCtx%PHI, AppCtx%MeshTopology%Num_Elems + i-1, PHI0, iErr); CHKERRQ(iErr)
               MyNumDec = MyNumDec + 1
            !Else
            !   Call MeshUpdateClosureInt(AppCtx%MeshTopology%Mesh, AppCtx%PHI, AppCtx%MeshTopology%Num_Elems + i-1, PHI1, iErr); CHKERRQ(iErr)
            End If
         End Do
         DeAllocate(RHSPhi_Ptr)
         DeAllocate(PHI0)         
         DeAllocate(PHI1)

      Call PetscGlobalSum(MyNumDec, NumDec, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      !!! Count the total number of debonded elements

      Write(IOBuffer,*) 'Number of debonded elements ', numdec, '\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      !Call PetscLogStagePop(iErr); CHKERRQ(iErr)
      

   End Subroutine Solve_Phi

End Module m_VarFracFilm_Phi
