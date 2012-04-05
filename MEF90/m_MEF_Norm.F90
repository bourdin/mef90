Module m_MEF_Norm
#include "finclude/petscdef.h"
   Use m_MEF_Parameters
   Use m_MEF_Elements
   Use m_MEF_Utils
   Use petsc
   IMPLICIT NONE

   Private   
   Public :: SectionRealL2DotProduct
   
   Interface SectionRealL2DotProduct
      Module Procedure SectionRealL2DotProductSet_2DScal,SectionRealL2DotProduct_2DScal, &
                       SectionRealL2DotProductSet_3DScal,SectionRealL2DotProduct_3DScal
   End Interface SectionRealL2DotProduct

Contains
#undef __FUNCT__
#define __FUNCT__ "SectionRealL2DotProductSet_2DScal"
   Subroutine SectionRealL2DotProductSet_2DScal(Mesh,setID,U_sec,V_Sec,Elem,L2DotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec,V_Sec
      Type(Element2D_Scal), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: L2DotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: elemIS
      PetscInt,Dimension(:),Pointer                :: elemID
      PetscInt                                     :: iE
      PetscInt                                     :: numDof,iDof,numGauss,iGauss
      PetscReal, Dimension(:), Pointer             :: U_loc,V_Loc
      PetscReal                                    :: U_elem,V_elem
      PetscReal                                    :: myL2DotProduct,L2DotProductloc
      MPI_Comm                                     :: comm

      L2DotProduct = 0.0_Kr
      myL2DotProduct = 0.0_Kr
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID,elemIS,ierr); CHKERRQ(ierr)
      Call ISGetIndicesF90(elemIS,elemID,ierr);CHKERRQ(ierr)
      If (size(elemID) > 0) Then 
         !!! Since all elements in  a set have the same type, we get numDof and numGauss from the first element
         numDof   = size(elem(elemID(1)+1)%BF,1)
         numGauss = size(elem(elemID(1)+1)%BF,2)
         Allocate(U_loc(numDof))
         Allocate(V_loc(numDof))
         Do iE = 1,size(elemID)
            Call SectionRealRestrictClosure(U_Sec,mesh,elemID(iE),numDof,U_Loc,ierr); CHKERRQ(ierr)
            Call SectionRealRestrictClosure(V_Sec,mesh,elemID(iE),numDof,V_Loc,ierr); CHKERRQ(ierr)
            Do iGauss = 1, numGauss
               U_elem = 0.0_Kr
               V_elem = 0.0_Kr
               Do iDof = 1, numDof
                  U_elem = U_elem + U_loc(iDof) * elem(elemID(iE)+1)%BF(iDof,iGauss)
                  V_elem = V_elem + V_loc(iDof) * elem(elemID(iE)+1)%BF(iDof,iGauss)
               End Do !iDof
               myL2DotProduct = myL2DotProduct + U_elem * V_elem * elem(elemID(iE)+1)%Gauss_C(iGauss)
            End Do !iGauss
         End Do !elem
         DeAllocate(U_loc)
         DeAllocate(V_loc)
      End If !(size(elemID) > 0)
      Call ISRestoreIndicesF90(elemIS,elemID,ierr);CHKERRQ(ierr)
      !!!
      !!! Accumulate myL2DotProduct across all CPUs
      !!!
      Call PetscObjectGetComm(mesh,comm,ierr);CHKERRQ(ierr)
      Call MPI_AllReduce(myL2DotProduct,L2DotProduct,1,MPIU_SCALAR,MPI_SUM,comm,ierr);CHKERRQ(ierr)
   End Subroutine SectionRealL2DotProductSet_2DScal
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealL2DotProductSet_3DScal"
   Subroutine SectionRealL2DotProductSet_3DScal(Mesh,setID,U_sec,V_Sec,Elem,L2DotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec,V_Sec
      Type(Element3D_Scal), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: L2DotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: elemIS
      PetscInt,Dimension(:),Pointer                :: elemID
      PetscInt                                     :: iE
      PetscInt                                     :: numDof,iDof,numGauss,iGauss
      PetscReal, Dimension(:), Pointer             :: U_loc,V_Loc
      PetscReal                                    :: U_elem,V_elem
      PetscReal                                    :: myL2DotProduct,L2DotProductloc
      MPI_Comm                                     :: comm

      L2DotProduct = 0.0_Kr
      myL2DotProduct = 0.0_Kr
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID,elemIS,ierr); CHKERRQ(ierr)
      Call ISGetIndicesF90(elemIS,elemID,ierr);CHKERRQ(ierr)
      If (size(elemID) > 0) Then 
         !!! Since all elements in  a set have the same type, we get numDof and numGauss from the first element
         numDof   = size(elem(elemID(1)+1)%BF,1)
         numGauss = size(elem(elemID(1)+1)%BF,2)
         Allocate(U_loc(numDof))
         Allocate(V_loc(numDof))
         Do iE = 1,size(elemID)
            Call SectionRealRestrictClosure(U_Sec,mesh,elemID(iE),numDof,U_Loc,ierr); CHKERRQ(ierr)
            Call SectionRealRestrictClosure(V_Sec,mesh,elemID(iE),numDof,V_Loc,ierr); CHKERRQ(ierr)
            Do iGauss = 1, numGauss
               U_elem = 0.0_Kr
               V_elem = 0.0_Kr
               Do iDof = 1, numDof
                  U_elem = U_elem + U_loc(iDof) * elem(elemID(iE)+1)%BF(iDof,iGauss)
                  V_elem = V_elem + V_loc(iDof) * elem(elemID(iE)+1)%BF(iDof,iGauss)
               End Do !iDof
               myL2DotProduct = myL2DotProduct + U_elem * V_elem * elem(elemID(iE)+1)%Gauss_C(iGauss)
            End Do !iGauss
         End Do !elem
         DeAllocate(U_loc)
         DeAllocate(V_loc)
      End If !(size(elemID) > 0)
      Call ISRestoreIndicesF90(elemIS,elemID,ierr);CHKERRQ(ierr)
      !!!
      !!! Accumulate myL2DotProduct across all CPUs
      !!!
      Call PetscObjectGetComm(mesh,comm,ierr);CHKERRQ(ierr)
      Call MPI_AllReduce(myL2DotProduct,L2DotProduct,1,MPIU_SCALAR,MPI_SUM,comm,ierr);CHKERRQ(ierr)
   End Subroutine SectionRealL2DotProductSet_3DScal
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealL2DotProduct_2DScal"
   Subroutine SectionRealL2DotProduct_2DScal(Mesh,U_sec,V_Sec,Elem,L2DotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec,V_Sec
      Type(Element2D_Scal), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: L2DotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: L2DotProductSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      L2DotProductSet = 0.0_Kr
      L2DotProduct = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealL2DotProductSet_2DScal(Mesh,setID(set),U_sec,V_Sec,Elem,L2DotProductSet,ierr);CHKERRQ(ierr)
         L2DotProduct = L2DotProduct + L2DotProductSet
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   End Subroutine SectionRealL2DotProduct_2DScal

#undef __FUNCT__
#define __FUNCT__ "SectionRealL2DotProduct_3DScal"
   Subroutine SectionRealL2DotProduct_3DScal(Mesh,U_sec,V_Sec,Elem,L2DotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec,V_Sec
      Type(Element3D_Scal), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: L2DotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: L2DotProductSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      L2DotProductSet = 0.0_Kr
      L2DotProduct = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealL2DotProductSet_3DScal(Mesh,setID(set),U_sec,V_Sec,Elem,L2DotProductSet,ierr);CHKERRQ(ierr)
         L2DotProduct = L2DotProduct + L2DotProductSet
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   End Subroutine SectionRealL2DotProduct_3DScal
End Module m_MEF_Norm