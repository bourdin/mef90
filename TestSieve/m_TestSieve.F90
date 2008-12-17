Module m_TestSieve
  Use m_MEF90

  Implicit NONE
  PRIVATE

#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscao.h"

   Public :: FormObjectiveFunction
!   Public :: FormGradient
!   Public :: FormHessian

 Contains
   Subroutine FormObjectiveFunction(dMyObjFunc, dMyMeshTopology, dMyElem, dU_Loc, dF_Loc)
   !!! Compute the contribution of each CPU to the objective function
   !!! \int 1/2 |\nabla u|^2 - f.u \, dx
   !!! dU_Loc must be a LOCAL (ghosted) vector 
   
      PetscReal                                     :: dMyObjFunc
      Type (MeshTopology_Info)                      :: dMyMeshTopology
      Type (Element2D_Scal), Dimension(:), Pointer  :: dMyElem
      Vec                                           :: dU_Loc, dF_Loc
      
      Integer                                       :: iBlk, iELoc, iE, iG, iSL, iSG
      PetscReal, Dimension(:), Pointer              :: U_Ptr, F_Ptr
      PetscReal                                     :: U_Elem, F_Elem
      Type (Vect2D)                                 :: Strain_Elem
      Integer                                       :: iErr
      
      dMyObjFunc = 0.0_Kr
      Call VecGetArrayF90(dU_Loc, U_Ptr, iErr)
      Call VecGetArrayF90(dF_Loc, F_Ptr, iErr)
      
      Do_iBlk: Do iBlk = 1, dMyMeshTopology%Num_Elem_Blks
         Do_iE: Do iELoc = 1, dMyMeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = dMyMeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            Do_iG: Do iG = 1, size(dMyElem(iE)%BF,2)
               Strain_Elem = 0.0_Kr
               U_Elem      = 0.0_Kr
               F_Elem      = 0.0_Kr
              
               Do_iSL: Do iSL = 1, dMyMeshTopology%Elem_Blk(iBlk)%Nb_DoF
                  iSG = dMyElem(iE)%ID_DoF(iSL)
                  Strain_Elem = Strain_Elem + dMyElem(iE)%Grad_BF(iSL, iG) * U_Ptr(iSG)
                  F_Elem      = F_Elem      + dMyElem(iE)%BF(iSL, iG)      * F_Ptr(iSG)
                  U_Elem      = U_Elem      + dMyElem(iE)%BF(iSL, iG)      * U_Ptr(iSG)
               End Do Do_iSL
               dMyObjFunc = dMyObjFunc + dMyELem(iE)%Gauss_C(iG) * ( (Strain_Elem .DotP. Strain_Elem) * 0.5_Kr - F_Elem * U_Elem) 
            End Do Do_iG
         End Do Do_iE
      End Do Do_iBlk
      Call VecRestoreArrayF90(dU_Loc, U_Ptr, iErr)
      Call VecRestoreArrayF90(dF_Loc, F_Ptr, iErr)
   End Subroutine FormObjectiveFunction
   
End Module m_TestSieve
