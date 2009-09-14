Module m_TestSieve
#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscmeshdef.h"

   Use m_MEF90
   Use petsc
   Use petscmesh
   
   Implicit NONE

   Public :: FormObjectiveFunction
!   Public :: FormGradient
!   Public :: FormHessian

 Contains
   Subroutine FormObjectiveFunction(dMyObjFunc, dMyMeshTopology, dMyElem, dU, dF, integrationEvent)
   !!! Compute the contribution of each CPU to the objective function
   !!! \int 1/2 |\nabla u|^2 - f.u \, dx
   !!! dU_Loc must be a LOCAL (ghosted) vector 
   
      PetscReal                                     :: dMyObjFunc
      Type(MeshTopology_Type)                       :: dMyMeshTopology
      Type(Element2D_Scal), Dimension(:), Pointer   :: dMyElem
      Type(SectionReal)                             :: dU, dF
      PetscLogEvent                                 :: integrationEvent
      
      PetscInt                                      :: iBlk, blkId, numElems, iELoc, iE, iG, iSL
      Character(len=256)                            :: CharBuffer
      PetscInt, Dimension(:), Pointer               :: blkIds, elemIds
      PetscScalar, Dimension(:), Pointer            :: U_Ptr, F_Ptr
      PetscReal                                     :: U_Elem, F_Elem
      Type(Vect2D)                                  :: Strain_Elem
      PetscErrorCode                                :: iErr

      call PetscLogEventBegin(integrationEvent, ierr); CHKERRQ(ierr)
      Allocate(U_Ptr(dMyMeshTopology%Elem_Blk(1)%Num_DoF))
      Allocate(F_Ptr(dMyMeshTopology%Elem_Blk(1)%Num_DoF))

      dMyObjFunc = 0.0_Kr
      CharBuffer = 'CellBlocks'
      Allocate(blkIds(dMyMeshTopology%Num_Elem_blks))
      call MeshGetLabelIds(dMyMeshTopology%mesh, CharBuffer, blkIds, ierr); CHKERRQ(ierr)
      Do_iBlk: Do iBlk = 1, dMyMeshTopology%Num_Elem_Blks
         blkId = blkIds(iBlk)
         call MeshGetStratumSize(dMyMeshTopology%mesh, CharBuffer, blkId, numElems, ierr); CHKERRQ(ierr)
         Allocate(elemIds(numElems))
         call MeshGetStratum(dMyMeshTopology%mesh, CharBuffer, blkId, elemIds, ierr); CHKERRQ(ierr)
         Do_iE: Do iELoc = 1, numElems
            iE = elemIds(iELoc)+1
      !Do_iBlk: Do iBlk = 1, dMyMeshTopology%Num_Elem_Blks
      !   Do_iE: Do iELoc = 1, dMyMeshTopology%Elem_Blk(iBlk)%Num_Elems
      !      iE = dMyMeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            call MeshRestrictClosure(dMyMeshTopology%mesh, dU, iE-1, Size(U_Ptr), U_Ptr, ierr); CHKERRQ(ierr)
            call MeshRestrictClosure(dMyMeshTopology%mesh, dF, iE-1, Size(F_Ptr), F_Ptr, ierr); CHKERRQ(ierr)
            Do_iG: Do iG = 1, size(dMyElem(iE)%BF,2)
               Strain_Elem = 0.0_Kr
               U_Elem      = 0.0_Kr
               F_Elem      = 0.0_Kr
              
               Do_iSL: Do iSL = 1, dMyMeshTopology%Elem_Blk(iBlk)%Num_DoF
                  Strain_Elem = Strain_Elem + dMyElem(iE)%Grad_BF(iSL, iG) * U_Ptr(iSL)
                  F_Elem      = F_Elem      + dMyElem(iE)%BF(iSL, iG)      * F_Ptr(iSL)
                  U_Elem      = U_Elem      + dMyElem(iE)%BF(iSL, iG)      * U_Ptr(iSL)
                  call PetscLogFlops(6._Kr, ierr)
               End Do Do_iSL
               dMyObjFunc = dMyObjFunc + dMyELem(iE)%Gauss_C(iG) * ( (Strain_Elem .DotP. Strain_Elem) * 0.5_Kr - F_Elem * U_Elem) 
               call PetscLogFlops(5 + dMyMeshTopology%num_dim*2 - 1, ierr)
            End Do Do_iG
         End Do Do_iE
         Deallocate(elemIds)
      End Do Do_iBlk
      Deallocate(blkIds)

      Deallocate(U_Ptr)
      Deallocate(F_Ptr)
      call PetscLogEventEnd(integrationEvent, ierr); CHKERRQ(ierr)
   End Subroutine FormObjectiveFunction
   
End Module m_TestSieve
