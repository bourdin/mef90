Module m_MEF_Norm
#include "finclude/petscdef.h"
   Use m_MEF_LinAlg
   Use m_MEF_Parameters
   Use m_MEF_Elements
   Use m_MEF_Utils
   Use petsc
   IMPLICIT NONE

   Private   
   Public :: SectionRealL2DotProduct,SectionRealH1semiDotProduct
   
   Interface SectionRealL2DotProduct
      Module Procedure SectionRealL2DotProductSet_2DScal,SectionRealL2DotProduct_2DScal, &
                       SectionRealL2DotProductSet_2DVect,SectionRealL2DotProduct_2DVect, &
                       SectionRealL2DotProductSet_2DElast,SectionRealL2DotProduct_2DElast, &
                       SectionRealL2DotProductSet_3DScal,SectionRealL2DotProduct_3DScal, &
                       SectionRealL2DotProductSet_3DVect,SectionRealL2DotProduct_3DVect, &
                       SectionRealL2DotProductSet_3DElast,SectionRealL2DotProduct_3DElast
   End Interface SectionRealL2DotProduct

   Interface SectionRealL2Norm
      Module Procedure SectionRealL2NormSet_2DScal,SectionRealL2Norm_2DScal, &
                       SectionRealL2NormSet_2DVect,SectionRealL2Norm_2DVect, &
                       SectionRealL2NormSet_2DElast,SectionRealL2Norm_2DElast, &
                       SectionRealL2NormSet_3DScal,SectionRealL2Norm_3DScal, &
                       SectionRealL2NormSet_3DVect,SectionRealL2Norm_3DVect, &
                       SectionRealL2NormSet_3DElast,SectionRealL2Norm_3DElast
   End Interface SectionRealL2Norm

   Interface SectionRealH1semiDotProduct
      Module Procedure SectionRealH1semiDotProductSet_2DScal,SectionRealH1semiDotProduct_2DScal, &
                       SectionRealH1semiDotProductSet_2DVect,SectionRealH1semiDotProduct_2DVect, &
                       SectionRealH1semiDotProductSet_2DElast,SectionRealH1semiDotProduct_2DElast, &
                       SectionRealH1semiDotProductSet_3DScal,SectionRealH1semiDotProduct_3DScal, &
                       SectionRealH1semiDotProductSet_3DVect,SectionRealH1semiDotProduct_3DVect, &
                       SectionRealH1semiDotProductSet_3DElast,SectionRealH1semiDotProduct_3DElast
   End Interface SectionRealH1semiDotProduct

   Interface SectionRealH1semiNorm
      Module Procedure SectionRealH1semiNormSet_2DScal,SectionRealH1semiNorm_2DScal, &
                       SectionRealH1semiNormSet_2DVect,SectionRealH1semiNorm_2DVect, &
                       SectionRealH1semiNormSet_2DElast,SectionRealH1semiNorm_2DElast, &
                       SectionRealH1semiNormSet_3DScal,SectionRealH1semiNorm_3DScal, &
                       SectionRealH1semiNormSet_3DVect,SectionRealH1semiNorm_3DVect, &
                       SectionRealH1semiNormSet_3DElast,SectionRealH1semiNorm_3DElast
   End Interface SectionRealH1semiNorm

   Interface SectionRealH1DotProduct
      Module Procedure SectionRealH1DotProductSet_2DScal,SectionRealH1DotProduct_2DScal, &
                       SectionRealH1DotProductSet_2DVect,SectionRealH1DotProduct_2DVect, &
                       SectionRealH1DotProductSet_2DElast,SectionRealH1DotProduct_2DElast, &
                       SectionRealH1DotProductSet_3DScal,SectionRealH1DotProduct_3DScal, &
                       SectionRealH1DotProductSet_3DVect,SectionRealH1DotProduct_3DVect, &
                       SectionRealH1DotProductSet_3DElast,SectionRealH1DotProduct_3DElast
   End Interface SectionRealH1DotProduct

   Interface SectionRealH1Norm
      Module Procedure SectionRealH1NormSet_2DScal,SectionRealH1Norm_2DScal, &
                       SectionRealH1NormSet_2DVect,SectionRealH1Norm_2DVect, &
                       SectionRealH1NormSet_2DElast,SectionRealH1Norm_2DElast, &
                       SectionRealH1NormSet_3DScal,SectionRealH1Norm_3DScal, &
                       SectionRealH1NormSet_3DVect,SectionRealH1Norm_3DVect, &
                       SectionRealH1NormSet_3DElast,SectionRealH1Norm_3DElast
   End Interface SectionRealH1Norm

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
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID,elemIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(elemIS,elemID,ierr);CHKERRQ(ierr)
      If (size(elemID) > 0) Then 
         !!! Since all elements in  a set have the same type, we get numDof and numGauss from the first element
         numDof   = size(elem(elemID(1)+1)%BF,1)
         numGauss = size(elem(elemID(1)+1)%BF,2)
         Allocate(U_loc(numDof))
         Allocate(V_loc(numDof))
         Do iE = 1,size(elemID)
            Call SectionRealRestrictClosure(U_Sec,mesh,elemID(iE),numDof,U_Loc,ierr);CHKERRQ(ierr)
            Call SectionRealRestrictClosure(V_Sec,mesh,elemID(iE),numDof,V_Loc,ierr);CHKERRQ(ierr)
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
#define __FUNCT__ "SectionRealL2DotProductSet_2DVect"
   Subroutine SectionRealL2DotProductSet_2DVect(Mesh,setID,U_sec,V_Sec,Elem,L2DotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec,V_Sec
      Type(Element2D_Vect), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: L2DotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: elemIS
      PetscInt,Dimension(:),Pointer                :: elemID
      PetscInt                                     :: iE
      PetscInt                                     :: numDof,iDof,numGauss,iGauss
      PetscReal, Dimension(:), Pointer             :: U_loc,V_Loc
      Type(Vect2D)                                 :: U_elem,V_elem
      PetscReal                                    :: myL2DotProduct,L2DotProductloc
      MPI_Comm                                     :: comm

      L2DotProduct = 0.0_Kr
      myL2DotProduct = 0.0_Kr
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID,elemIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(elemIS,elemID,ierr);CHKERRQ(ierr)
      If (size(elemID) > 0) Then 
         !!! Since all elements in  a set have the same type, we get numDof and numGauss from the first element
         numDof   = size(elem(elemID(1)+1)%BF,1)
         numGauss = size(elem(elemID(1)+1)%BF,2)
         Allocate(U_loc(numDof))
         Allocate(V_loc(numDof))
         Do iE = 1,size(elemID)
            Call SectionRealRestrictClosure(U_Sec,mesh,elemID(iE),numDof,U_Loc,ierr);CHKERRQ(ierr)
            Call SectionRealRestrictClosure(V_Sec,mesh,elemID(iE),numDof,V_Loc,ierr);CHKERRQ(ierr)
            Do iGauss = 1, numGauss
               U_elem = 0.0_Kr
               V_elem = 0.0_Kr
               Do iDof = 1, numDof
                  U_elem = U_elem + elem(elemID(iE)+1)%BF(iDof,iGauss) * U_loc(iDof) 
                  V_elem = V_elem + elem(elemID(iE)+1)%BF(iDof,iGauss) * V_loc(iDof) 
               End Do !iDof
               myL2DotProduct = myL2DotProduct + (U_elem .DotP. V_elem) * elem(elemID(iE)+1)%Gauss_C(iGauss)
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
   End Subroutine SectionRealL2DotProductSet_2DVect
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealL2DotProductSet_2DElast"
   Subroutine SectionRealL2DotProductSet_2DElast(Mesh,setID,U_sec,V_Sec,Elem,L2DotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec,V_Sec
      Type(Element2D_Elast), Dimension(:), Pointer :: Elem
      PetscReal,Intent(OUT)                        :: L2DotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: elemIS
      PetscInt,Dimension(:),Pointer                :: elemID
      PetscInt                                     :: iE
      PetscInt                                     :: numDof,iDof,numGauss,iGauss
      PetscReal, Dimension(:), Pointer             :: U_loc,V_Loc
      Type(Vect2D)                                 :: U_elem,V_elem
      PetscReal                                    :: myL2DotProduct,L2DotProductloc
      MPI_Comm                                     :: comm

      L2DotProduct = 0.0_Kr
      myL2DotProduct = 0.0_Kr
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID,elemIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(elemIS,elemID,ierr);CHKERRQ(ierr)
      If (size(elemID) > 0) Then 
         !!! Since all elements in  a set have the same type, we get numDof and numGauss from the first element
         numDof   = size(elem(elemID(1)+1)%BF,1)
         numGauss = size(elem(elemID(1)+1)%BF,2)
         Allocate(U_loc(numDof))
         Allocate(V_loc(numDof))
         Do iE = 1,size(elemID)
            Call SectionRealRestrictClosure(U_Sec,mesh,elemID(iE),numDof,U_Loc,ierr);CHKERRQ(ierr)
            Call SectionRealRestrictClosure(V_Sec,mesh,elemID(iE),numDof,V_Loc,ierr);CHKERRQ(ierr)
            Do iGauss = 1, numGauss
               U_elem = 0.0_Kr
               V_elem = 0.0_Kr
               Do iDof = 1, numDof
                  U_elem = U_elem + elem(elemID(iE)+1)%BF(iDof,iGauss) * U_loc(iDof) 
                  V_elem = V_elem + elem(elemID(iE)+1)%BF(iDof,iGauss) * V_loc(iDof) 
               End Do !iDof
               myL2DotProduct = myL2DotProduct + (U_elem .DotP. V_elem) * elem(elemID(iE)+1)%Gauss_C(iGauss)
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
   End Subroutine SectionRealL2DotProductSet_2DElast
   
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
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID,elemIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(elemIS,elemID,ierr);CHKERRQ(ierr)
      If (size(elemID) > 0) Then 
         !!! Since all elements in  a set have the same type, we get numDof and numGauss from the first element
         numDof   = size(elem(elemID(1)+1)%BF,1)
         numGauss = size(elem(elemID(1)+1)%BF,2)
         Allocate(U_loc(numDof))
         Allocate(V_loc(numDof))
         Do iE = 1,size(elemID)
            Call SectionRealRestrictClosure(U_Sec,mesh,elemID(iE),numDof,U_Loc,ierr);CHKERRQ(ierr)
            Call SectionRealRestrictClosure(V_Sec,mesh,elemID(iE),numDof,V_Loc,ierr);CHKERRQ(ierr)
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
#define __FUNCT__ "SectionRealL2DotProductSet_3DVect"
   Subroutine SectionRealL2DotProductSet_3DVect(Mesh,setID,U_sec,V_Sec,Elem,L2DotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec,V_Sec
      Type(Element3D_Vect), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: L2DotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: elemIS
      PetscInt,Dimension(:),Pointer                :: elemID
      PetscInt                                     :: iE
      PetscInt                                     :: numDof,iDof,numGauss,iGauss
      PetscReal, Dimension(:), Pointer             :: U_loc,V_Loc
      Type(Vect3D)                                 :: U_elem,V_elem
      PetscReal                                    :: myL2DotProduct,L2DotProductloc
      MPI_Comm                                     :: comm

      L2DotProduct = 0.0_Kr
      myL2DotProduct = 0.0_Kr
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID,elemIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(elemIS,elemID,ierr);CHKERRQ(ierr)
      If (size(elemID) > 0) Then 
         !!! Since all elements in  a set have the same type, we get numDof and numGauss from the first element
         numDof   = size(elem(elemID(1)+1)%BF,1)
         numGauss = size(elem(elemID(1)+1)%BF,2)
         Allocate(U_loc(numDof))
         Allocate(V_loc(numDof))
         Do iE = 1,size(elemID)
            Call SectionRealRestrictClosure(U_Sec,mesh,elemID(iE),numDof,U_Loc,ierr);CHKERRQ(ierr)
            Call SectionRealRestrictClosure(V_Sec,mesh,elemID(iE),numDof,V_Loc,ierr);CHKERRQ(ierr)
            Do iGauss = 1, numGauss
               U_elem = 0.0_Kr
               V_elem = 0.0_Kr
               Do iDof = 1, numDof
                  U_elem = U_elem + elem(elemID(iE)+1)%BF(iDof,iGauss) * U_loc(iDof) 
                  V_elem = V_elem + elem(elemID(iE)+1)%BF(iDof,iGauss) * V_loc(iDof) 
               End Do !iDof
               myL2DotProduct = myL2DotProduct + (U_elem .DotP. V_elem) * elem(elemID(iE)+1)%Gauss_C(iGauss)
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
   End Subroutine SectionRealL2DotProductSet_3DVect
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealL2DotProductSet_3DElast"
   Subroutine SectionRealL2DotProductSet_3DElast(Mesh,setID,U_sec,V_Sec,Elem,L2DotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec,V_Sec
      Type(Element3D_Elast), Dimension(:), Pointer :: Elem
      PetscReal,Intent(OUT)                        :: L2DotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: elemIS
      PetscInt,Dimension(:),Pointer                :: elemID
      PetscInt                                     :: iE
      PetscInt                                     :: numDof,iDof,numGauss,iGauss
      PetscReal, Dimension(:), Pointer             :: U_loc,V_Loc
      Type(Vect3D)                                 :: U_elem,V_elem
      PetscReal                                    :: myL2DotProduct,L2DotProductloc
      MPI_Comm                                     :: comm

      L2DotProduct = 0.0_Kr
      myL2DotProduct = 0.0_Kr
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID,elemIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(elemIS,elemID,ierr);CHKERRQ(ierr)
      If (size(elemID) > 0) Then 
         !!! Since all elements in  a set have the same type, we get numDof and numGauss from the first element
         numDof   = size(elem(elemID(1)+1)%BF,1)
         numGauss = size(elem(elemID(1)+1)%BF,2)
         Allocate(U_loc(numDof))
         Allocate(V_loc(numDof))
         Do iE = 1,size(elemID)
            Call SectionRealRestrictClosure(U_Sec,mesh,elemID(iE),numDof,U_Loc,ierr);CHKERRQ(ierr)
            Call SectionRealRestrictClosure(V_Sec,mesh,elemID(iE),numDof,V_Loc,ierr);CHKERRQ(ierr)
            Do iGauss = 1, numGauss
               U_elem = 0.0_Kr
               V_elem = 0.0_Kr
               Do iDof = 1, numDof
                  U_elem = U_elem + elem(elemID(iE)+1)%BF(iDof,iGauss) * U_loc(iDof) 
                  V_elem = V_elem + elem(elemID(iE)+1)%BF(iDof,iGauss) * V_loc(iDof) 
               End Do !iDof
               myL2DotProduct = myL2DotProduct + (U_elem .DotP. V_elem) * elem(elemID(iE)+1)%Gauss_C(iGauss)
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
   End Subroutine SectionRealL2DotProductSet_3DElast

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
#define __FUNCT__ "SectionRealL2DotProduct_2DVect"
   Subroutine SectionRealL2DotProduct_2DVect(Mesh,U_sec,V_Sec,Elem,L2DotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec,V_Sec
      Type(Element2D_Vect), Dimension(:), Pointer  :: Elem
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
         Call SectionRealL2DotProductSet_2DVect(Mesh,setID(set),U_sec,V_Sec,Elem,L2DotProductSet,ierr);CHKERRQ(ierr)
         L2DotProduct = L2DotProduct + L2DotProductSet
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   End Subroutine SectionRealL2DotProduct_2DVect

#undef __FUNCT__
#define __FUNCT__ "SectionRealL2DotProduct_2DElast"
   Subroutine SectionRealL2DotProduct_2DElast(Mesh,U_sec,V_Sec,Elem,L2DotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec,V_Sec
      Type(Element2D_Elast), Dimension(:), Pointer :: Elem
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
         Call SectionRealL2DotProductSet_2DElast(Mesh,setID(set),U_sec,V_Sec,Elem,L2DotProductSet,ierr);CHKERRQ(ierr)
         L2DotProduct = L2DotProduct + L2DotProductSet
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   End Subroutine SectionRealL2DotProduct_2DElast

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

#undef __FUNCT__
#define __FUNCT__ "SectionRealL2DotProduct_3DVect"
   Subroutine SectionRealL2DotProduct_3DVect(Mesh,U_sec,V_Sec,Elem,L2DotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec,V_Sec
      Type(Element3D_Vect), Dimension(:), Pointer  :: Elem
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
         Call SectionRealL2DotProductSet_3DVect(Mesh,setID(set),U_sec,V_Sec,Elem,L2DotProductSet,ierr);CHKERRQ(ierr)
         L2DotProduct = L2DotProduct + L2DotProductSet
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   End Subroutine SectionRealL2DotProduct_3DVect

#undef __FUNCT__
#define __FUNCT__ "SectionRealL2DotProduct_3DElast"
   Subroutine SectionRealL2DotProduct_3DElast(Mesh,U_sec,V_Sec,Elem,L2DotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec,V_Sec
      Type(Element3D_Elast), Dimension(:), Pointer :: Elem
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
         Call SectionRealL2DotProductSet_3DElast(Mesh,setID(set),U_sec,V_Sec,Elem,L2DotProductSet,ierr);CHKERRQ(ierr)
         L2DotProduct = L2DotProduct + L2DotProductSet
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   End Subroutine SectionRealL2DotProduct_3DElast

#undef __FUNCT__
#define __FUNCT__ "SectionRealL2NormSet_2DScal"
   Subroutine SectionRealL2NormSet_2DScal(Mesh,setID,U_sec,Elem,L2NormSet,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec
      Type(Element2D_Scal), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: L2NormSet
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Call SectionRealL2DotProductSet_2DScal(Mesh,setID,U_sec,U_Sec,Elem,L2NormSet,ierr);CHKERRQ(ierr)
      L2NormSet = sqrt(L2NormSet)
   End Subroutine SectionRealL2NormSet_2DScal
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealL2NormSet_2DVect"
   Subroutine SectionRealL2NormSet_2DVect(Mesh,setID,U_sec,Elem,L2NormSet,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec
      Type(Element2D_Vect), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: L2NormSet
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Call SectionRealL2DotProductSet_2DVect(Mesh,setID,U_sec,U_Sec,Elem,L2NormSet,ierr);CHKERRQ(ierr)
      L2NormSet = sqrt(L2NormSet)
   End Subroutine SectionRealL2NormSet_2DVect
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealL2NormSet_2DElast"
   Subroutine SectionRealL2NormSet_2DElast(Mesh,setID,U_sec,Elem,L2NormSet,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec
      Type(Element2D_Elast), Dimension(:), Pointer :: Elem
      PetscReal,Intent(OUT)                        :: L2NormSet
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Call SectionRealL2DotProductSet_2DElast(Mesh,setID,U_sec,U_Sec,Elem,L2NormSet,ierr);CHKERRQ(ierr)
      L2NormSet = sqrt(L2NormSet)
   End Subroutine SectionRealL2NormSet_2DElast
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealL2NormSet_3DScal"
   Subroutine SectionRealL2NormSet_3DScal(Mesh,setID,U_sec,Elem,L2NormSet,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec
      Type(Element3D_Scal), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: L2NormSet
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Call SectionRealL3DotProductSet_3DScal(Mesh,setID,U_sec,U_Sec,Elem,L2NormSet,ierr);CHKERRQ(ierr)
      L2NormSet = sqrt(L2NormSet)
   End Subroutine SectionRealL2NormSet_3DScal
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealL2NormSet_3DVect"
   Subroutine SectionRealL2NormSet_3DVect(Mesh,setID,U_sec,Elem,L2NormSet,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec
      Type(Element3D_Vect), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: L2NormSet
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Call SectionRealL3DotProductSet_3DVect(Mesh,setID,U_sec,U_Sec,Elem,L2NormSet,ierr);CHKERRQ(ierr)
      L2NormSet = sqrt(L2NormSet)
   End Subroutine SectionRealL2NormSet_3DVect
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealL2NormSet_3DElast"
   Subroutine SectionRealL2NormSet_3DElast(Mesh,setID,U_sec,Elem,L2NormSet,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec
      Type(Element3D_Elast), Dimension(:), Pointer :: Elem
      PetscReal,Intent(OUT)                        :: L2NormSet
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Call SectionRealL3DotProductSet_3DElast(Mesh,setID,U_sec,U_Sec,Elem,L2NormSet,ierr);CHKERRQ(ierr)
      L2NormSet = sqrt(L2NormSet)
   End Subroutine SectionRealL2NormSet_3DElast
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealL2Norm_2DScal"
   Subroutine SectionRealL2Norm_2DScal(Mesh,U_sec,Elem,L2Norm,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec
      Type(Element2D_Scal), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: L2Norm
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: L2NormSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      L2NormSet = 0.0_Kr
      L2Norm = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealL2NormSet_2DScal(Mesh,setID(set),U_sec,Elem,L2NormSet,ierr);CHKERRQ(ierr)
         L2Norm = L2Norm + L2NormSet**2
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      L2Norm = sqrt(L2Norm)
   End Subroutine SectionRealL2Norm_2DScal
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealL2Norm_2DVect"
   Subroutine SectionRealL2Norm_2DVect(Mesh,U_sec,Elem,L2Norm,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec
      Type(Element2D_Vect), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: L2Norm
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: L2NormSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      L2NormSet = 0.0_Kr
      L2Norm = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealL2NormSet_2DVect(Mesh,setID(set),U_sec,Elem,L2NormSet,ierr);CHKERRQ(ierr)
         L2Norm = L2Norm + L2NormSet**2
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      L2Norm = sqrt(L2Norm)
   End Subroutine SectionRealL2Norm_2DVect
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealL2Norm_2DElast"
   Subroutine SectionRealL2Norm_2DElast(Mesh,U_sec,Elem,L2Norm,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec
      Type(Element2D_Elast), Dimension(:), Pointer :: Elem
      PetscReal,Intent(OUT)                        :: L2Norm
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: L2NormSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      L2NormSet = 0.0_Kr
      L2Norm = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealL2NormSet_2DElast(Mesh,setID(set),U_sec,Elem,L2NormSet,ierr);CHKERRQ(ierr)
         L2Norm = L2Norm + L2NormSet**2
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      L2Norm = sqrt(L2Norm)
   End Subroutine SectionRealL2Norm_2DElast
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealL2Norm_3DScal"
   Subroutine SectionRealL2Norm_3DScal(Mesh,U_sec,Elem,L2Norm,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec
      Type(Element3D_Scal), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: L2Norm
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: L2NormSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      L2NormSet = 0.0_Kr
      L2Norm = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealL2NormSet_3DScal(Mesh,setID(set),U_sec,Elem,L2NormSet,ierr);CHKERRQ(ierr)
         L2Norm = L2Norm + L2NormSet**2
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      L2Norm = sqrt(L2Norm)
   End Subroutine SectionRealL2Norm_3DScal
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealL2Norm_3DVect"
   Subroutine SectionRealL2Norm_3DVect(Mesh,U_sec,Elem,L2Norm,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec
      Type(Element3D_Vect), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: L2Norm
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: L2NormSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      L2NormSet = 0.0_Kr
      L2Norm = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealL2NormSet_3DVect(Mesh,setID(set),U_sec,Elem,L2NormSet,ierr);CHKERRQ(ierr)
         L2Norm = L2Norm + L2NormSet**2
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      L2Norm = sqrt(L2Norm)
   End Subroutine SectionRealL2Norm_3DVect
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealL2Norm_3DElast"
   Subroutine SectionRealL2Norm_3DElast(Mesh,U_sec,Elem,L2Norm,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec
      Type(Element3D_Elast), Dimension(:), Pointer :: Elem
      PetscReal,Intent(OUT)                        :: L2Norm
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: L2NormSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      L2NormSet = 0.0_Kr
      L2Norm = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealL2NormSet_3DElast(Mesh,setID(set),U_sec,Elem,L2NormSet,ierr);CHKERRQ(ierr)
         L2Norm = L2Norm + L2NormSet**2
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      L2Norm = sqrt(L2Norm)
   End Subroutine SectionRealL2Norm_3DElast

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1semiDotProductSet_2DScal"
   Subroutine SectionRealH1semiDotProductSet_2DScal(Mesh,setID,U_sec,V_Sec,Elem,H1semiDotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec,V_Sec
      Type(Element2D_Scal), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1semiDotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: elemIS
      PetscInt,Dimension(:),Pointer                :: elemID
      PetscInt                                     :: iE
      PetscInt                                     :: numDof,iDof,numGauss,iGauss
      PetscReal, Dimension(:), Pointer             :: U_loc,V_Loc
      Type(Vect2D)                                 :: GradU_elem,GradV_elem
      PetscReal                                    :: myH1semiDotProduct
      MPI_Comm                                     :: comm

      H1semiDotProduct = 0.0_Kr
      myH1semiDotProduct = 0.0_Kr
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID,elemIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(elemIS,elemID,ierr);CHKERRQ(ierr)
      If (size(elemID) > 0) Then 
         !!! Since all elements in  a set have the same type, we get numDof and numGauss from the first element
         numDof   = size(elem(elemID(1)+1)%BF,1)
         numGauss = size(elem(elemID(1)+1)%BF,2)
         Allocate(U_loc(numDof))
         Allocate(V_loc(numDof))
         Do iE = 1,size(elemID)
            Call SectionRealRestrictClosure(U_Sec,mesh,elemID(iE),numDof,U_Loc,ierr);CHKERRQ(ierr)
            Call SectionRealRestrictClosure(V_Sec,mesh,elemID(iE),numDof,V_Loc,ierr);CHKERRQ(ierr)
            Do iGauss = 1, numGauss
               GradU_elem = 0.0_Kr
               GradV_elem = 0.0_Kr
               Do iDof = 1, numDof
                  GradU_elem = GradU_elem + U_loc(iDof) * elem(elemID(iE)+1)%Grad_BF(iDof,iGauss)
                  GradV_elem = GradV_elem + V_loc(iDof) * elem(elemID(iE)+1)%Grad_BF(iDof,iGauss)
               End Do !iDof
               myH1semiDotProduct = myH1semiDotProduct + (GradU_elem .DotP. GradV_elem) * elem(elemID(iE)+1)%Gauss_C(iGauss)
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
      Call MPI_AllReduce(myH1semiDotProduct,H1semiDotProduct,1,MPIU_SCALAR,MPI_SUM,comm,ierr);CHKERRQ(ierr)
   End Subroutine SectionRealH1semiDotProductSet_2DScal

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1semiDotProductSet_2DVect"
   Subroutine SectionRealH1semiDotProductSet_2DVect(Mesh,setID,U_sec,V_Sec,Elem,H1semiDotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec,V_Sec
      Type(Element2D_Vect), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1semiDotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: elemIS
      PetscInt,Dimension(:),Pointer                :: elemID
      PetscInt                                     :: iE
      PetscInt                                     :: numDof,iDof,numGauss,iGauss
      PetscReal, Dimension(:), Pointer             :: U_loc,V_Loc
      Type(Mat2D)                                  :: GradU_elem,GradV_elem
      PetscReal                                    :: myH1semiDotProduct
      MPI_Comm                                     :: comm

      H1semiDotProduct = 0.0_Kr
      myH1semiDotProduct = 0.0_Kr
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID,elemIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(elemIS,elemID,ierr);CHKERRQ(ierr)
      If (size(elemID) > 0) Then 
         !!! Since all elements in  a set have the same type, we get numDof and numGauss from the first element
         numDof   = size(elem(elemID(1)+1)%BF,1)
         numGauss = size(elem(elemID(1)+1)%BF,2)
         Allocate(U_loc(numDof))
         Allocate(V_loc(numDof))
         Do iE = 1,size(elemID)
            Call SectionRealRestrictClosure(U_Sec,mesh,elemID(iE),numDof,U_Loc,ierr);CHKERRQ(ierr)
            Call SectionRealRestrictClosure(V_Sec,mesh,elemID(iE),numDof,V_Loc,ierr);CHKERRQ(ierr)
            Do iGauss = 1, numGauss
               GradU_elem = 0.0_Kr
               GradV_elem = 0.0_Kr
               Do iDof = 1, numDof
                  GradU_elem = GradU_elem + U_loc(iDof) * elem(elemID(iE)+1)%Der_BF(iDof,iGauss)
                  GradV_elem = GradV_elem + V_loc(iDof) * elem(elemID(iE)+1)%Der_BF(iDof,iGauss)
               End Do !iDof
               myH1semiDotProduct = myH1semiDotProduct + (GradU_elem .DotP. GradV_elem) * elem(elemID(iE)+1)%Gauss_C(iGauss)
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
      Call MPI_AllReduce(myH1semiDotProduct,H1semiDotProduct,1,MPIU_SCALAR,MPI_SUM,comm,ierr);CHKERRQ(ierr)
   End Subroutine SectionRealH1semiDotProductSet_2DVect

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1semiDotProductSet_2DElast"
   Subroutine SectionRealH1semiDotProductSet_2DElast(Mesh,setID,U_sec,V_Sec,Elem,H1semiDotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec,V_Sec
      Type(Element2D_Elast), Dimension(:), Pointer :: Elem
      PetscReal,Intent(OUT)                        :: H1semiDotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: elemIS
      PetscInt,Dimension(:),Pointer                :: elemID
      PetscInt                                     :: iE
      PetscInt                                     :: numDof,iDof,numGauss,iGauss
      PetscReal, Dimension(:), Pointer             :: U_loc,V_Loc
      Type(MatS2D)                                 :: GradSU_elem,GradSV_elem
      PetscReal                                    :: myH1semiDotProduct
      MPI_Comm                                     :: comm

      H1semiDotProduct = 0.0_Kr
      myH1semiDotProduct = 0.0_Kr
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID,elemIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(elemIS,elemID,ierr);CHKERRQ(ierr)
      If (size(elemID) > 0) Then 
         !!! Since all elements in  a set have the same type, we get numDof and numGauss from the first element
         numDof   = size(elem(elemID(1)+1)%BF,1)
         numGauss = size(elem(elemID(1)+1)%BF,2)
         Allocate(U_loc(numDof))
         Allocate(V_loc(numDof))
         Do iE = 1,size(elemID)
            Call SectionRealRestrictClosure(U_Sec,mesh,elemID(iE),numDof,U_Loc,ierr);CHKERRQ(ierr)
            Call SectionRealRestrictClosure(V_Sec,mesh,elemID(iE),numDof,V_Loc,ierr);CHKERRQ(ierr)
            Do iGauss = 1, numGauss
               GradSU_elem = 0.0_Kr
               GradSV_elem = 0.0_Kr
               Do iDof = 1, numDof
                  GradSU_elem = GradSU_elem + U_loc(iDof) * elem(elemID(iE)+1)%GradS_BF(iDof,iGauss)
                  GradSV_elem = GradSV_elem + V_loc(iDof) * elem(elemID(iE)+1)%GradS_BF(iDof,iGauss)
               End Do !iDof
               myH1semiDotProduct = myH1semiDotProduct + (GradSU_elem .DotP. GradSV_elem) * elem(elemID(iE)+1)%Gauss_C(iGauss)
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
      Call MPI_AllReduce(myH1semiDotProduct,H1semiDotProduct,1,MPIU_SCALAR,MPI_SUM,comm,ierr);CHKERRQ(ierr)
   End Subroutine SectionRealH1semiDotProductSet_2DElast

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1semiDotProductSet_3DScal"
   Subroutine SectionRealH1semiDotProductSet_3DScal(Mesh,setID,U_sec,V_Sec,Elem,H1semiDotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec,V_Sec
      Type(Element3D_Scal), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1semiDotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: elemIS
      PetscInt,Dimension(:),Pointer                :: elemID
      PetscInt                                     :: iE
      PetscInt                                     :: numDof,iDof,numGauss,iGauss
      PetscReal, Dimension(:), Pointer             :: U_loc,V_Loc
      Type(Vect3D)                                 :: GradU_elem,GradV_elem
      PetscReal                                    :: myH1semiDotProduct
      MPI_Comm                                     :: comm

      H1semiDotProduct = 0.0_Kr
      myH1semiDotProduct = 0.0_Kr
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID,elemIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(elemIS,elemID,ierr);CHKERRQ(ierr)
      If (size(elemID) > 0) Then 
         !!! Since all elements in  a set have the same type, we get numDof and numGauss from the first element
         numDof   = size(elem(elemID(1)+1)%BF,1)
         numGauss = size(elem(elemID(1)+1)%BF,2)
         Allocate(U_loc(numDof))
         Allocate(V_loc(numDof))
         Do iE = 1,size(elemID)
            Call SectionRealRestrictClosure(U_Sec,mesh,elemID(iE),numDof,U_Loc,ierr);CHKERRQ(ierr)
            Call SectionRealRestrictClosure(V_Sec,mesh,elemID(iE),numDof,V_Loc,ierr);CHKERRQ(ierr)
            Do iGauss = 1, numGauss
               GradU_elem = 0.0_Kr
               GradV_elem = 0.0_Kr
               Do iDof = 1, numDof
                  GradU_elem = GradU_elem + U_loc(iDof) * elem(elemID(iE)+1)%Grad_BF(iDof,iGauss)
                  GradV_elem = GradV_elem + V_loc(iDof) * elem(elemID(iE)+1)%Grad_BF(iDof,iGauss)
               End Do !iDof
               myH1semiDotProduct = myH1semiDotProduct + (GradU_elem .DotP. GradV_elem) * elem(elemID(iE)+1)%Gauss_C(iGauss)
            End Do !iGauss
         End Do !elem
         DeAllocate(U_loc)
         DeAllocate(V_loc)
      End If !(size(elemID) > 0)
      Call ISRestoreIndicesF90(elemIS,elemID,ierr);CHKERRQ(ierr)
      !!!
      !!! Accumulate myL3DotProduct across all CPUs
      !!!
      Call PetscObjectGetComm(mesh,comm,ierr);CHKERRQ(ierr)
      Call MPI_AllReduce(myH1semiDotProduct,H1semiDotProduct,1,MPIU_SCALAR,MPI_SUM,comm,ierr);CHKERRQ(ierr)
   End Subroutine SectionRealH1semiDotProductSet_3DScal

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1semiDotProductSet_3DVect"
   Subroutine SectionRealH1semiDotProductSet_3DVect(Mesh,setID,U_sec,V_Sec,Elem,H1semiDotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec,V_Sec
      Type(Element3D_Vect), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1semiDotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: elemIS
      PetscInt,Dimension(:),Pointer                :: elemID
      PetscInt                                     :: iE
      PetscInt                                     :: numDof,iDof,numGauss,iGauss
      PetscReal, Dimension(:), Pointer             :: U_loc,V_Loc
      Type(Mat3D)                                  :: GradU_elem,GradV_elem
      PetscReal                                    :: myH1semiDotProduct
      MPI_Comm                                     :: comm

      H1semiDotProduct = 0.0_Kr
      myH1semiDotProduct = 0.0_Kr
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID,elemIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(elemIS,elemID,ierr);CHKERRQ(ierr)
      If (size(elemID) > 0) Then 
         !!! Since all elements in  a set have the same type, we get numDof and numGauss from the first element
         numDof   = size(elem(elemID(1)+1)%BF,1)
         numGauss = size(elem(elemID(1)+1)%BF,2)
         Allocate(U_loc(numDof))
         Allocate(V_loc(numDof))
         Do iE = 1,size(elemID)
            Call SectionRealRestrictClosure(U_Sec,mesh,elemID(iE),numDof,U_Loc,ierr);CHKERRQ(ierr)
            Call SectionRealRestrictClosure(V_Sec,mesh,elemID(iE),numDof,V_Loc,ierr);CHKERRQ(ierr)
            Do iGauss = 1, numGauss
               GradU_elem = 0.0_Kr
               GradV_elem = 0.0_Kr
               Do iDof = 1, numDof
                  GradU_elem = GradU_elem + U_loc(iDof) * elem(elemID(iE)+1)%Der_BF(iDof,iGauss)
                  GradV_elem = GradV_elem + V_loc(iDof) * elem(elemID(iE)+1)%Der_BF(iDof,iGauss)
               End Do !iDof
               myH1semiDotProduct = myH1semiDotProduct + (GradU_elem .DotP. GradV_elem) * elem(elemID(iE)+1)%Gauss_C(iGauss)
            End Do !iGauss
         End Do !elem
         DeAllocate(U_loc)
         DeAllocate(V_loc)
      End If !(size(elemID) > 0)
      Call ISRestoreIndicesF90(elemIS,elemID,ierr);CHKERRQ(ierr)
      !!!
      !!! Accumulate myL3DotProduct across all CPUs
      !!!
      Call PetscObjectGetComm(mesh,comm,ierr);CHKERRQ(ierr)
      Call MPI_AllReduce(myH1semiDotProduct,H1semiDotProduct,1,MPIU_SCALAR,MPI_SUM,comm,ierr);CHKERRQ(ierr)
   End Subroutine SectionRealH1semiDotProductSet_3DVect

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1semiDotProductSet_3DElast"
   Subroutine SectionRealH1semiDotProductSet_3DElast(Mesh,setID,U_sec,V_Sec,Elem,H1semiDotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec,V_Sec
      Type(Element3D_Elast), Dimension(:), Pointer :: Elem
      PetscReal,Intent(OUT)                        :: H1semiDotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: elemIS
      PetscInt,Dimension(:),Pointer                :: elemID
      PetscInt                                     :: iE
      PetscInt                                     :: numDof,iDof,numGauss,iGauss
      PetscReal, Dimension(:), Pointer             :: U_loc,V_Loc
      Type(MatS3D)                                 :: GradSU_elem,GradSV_elem
      PetscReal                                    :: myH1semiDotProduct
      MPI_Comm                                     :: comm

      H1semiDotProduct = 0.0_Kr
      myH1semiDotProduct = 0.0_Kr
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID,elemIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(elemIS,elemID,ierr);CHKERRQ(ierr)
      If (size(elemID) > 0) Then 
         !!! Since all elements in  a set have the same type, we get numDof and numGauss from the first element
         numDof   = size(elem(elemID(1)+1)%BF,1)
         numGauss = size(elem(elemID(1)+1)%BF,2)
         Allocate(U_loc(numDof))
         Allocate(V_loc(numDof))
         Do iE = 1,size(elemID)
            Call SectionRealRestrictClosure(U_Sec,mesh,elemID(iE),numDof,U_Loc,ierr);CHKERRQ(ierr)
            Call SectionRealRestrictClosure(V_Sec,mesh,elemID(iE),numDof,V_Loc,ierr);CHKERRQ(ierr)
            Do iGauss = 1, numGauss
               GradSU_elem = 0.0_Kr
               GradSV_elem = 0.0_Kr
               Do iDof = 1, numDof
                  GradSU_elem = GradSU_elem + U_loc(iDof) * elem(elemID(iE)+1)%GradS_BF(iDof,iGauss)
                  GradSV_elem = GradSV_elem + V_loc(iDof) * elem(elemID(iE)+1)%GradS_BF(iDof,iGauss)
               End Do !iDof
               myH1semiDotProduct = myH1semiDotProduct + (GradSU_elem .DotP. GradSV_elem) * elem(elemID(iE)+1)%Gauss_C(iGauss)
            End Do !iGauss
         End Do !elem
         DeAllocate(U_loc)
         DeAllocate(V_loc)
      End If !(size(elemID) > 0)
      Call ISRestoreIndicesF90(elemIS,elemID,ierr);CHKERRQ(ierr)
      !!!
      !!! Accumulate myL3DotProduct across all CPUs
      !!!
      Call PetscObjectGetComm(mesh,comm,ierr);CHKERRQ(ierr)
      Call MPI_AllReduce(myH1semiDotProduct,H1semiDotProduct,1,MPIU_SCALAR,MPI_SUM,comm,ierr);CHKERRQ(ierr)
   End Subroutine SectionRealH1semiDotProductSet_3DElast

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1semiDotProduct_2DScal"
   Subroutine SectionRealH1semiDotProduct_2DScal(Mesh,U_sec,V_Sec,Elem,H1semiDotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec,V_Sec
      Type(Element2D_Scal), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1semiDotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: H1semiDotProductSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      H1semiDotProductSet = 0.0_Kr
      H1semiDotProduct = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealH1semiDotProductSet_2DScal(Mesh,setID(set),U_sec,V_Sec,Elem,H1semiDotProductSet,ierr);CHKERRQ(ierr)
         H1semiDotProduct = H1semiDotProduct + H1semiDotProductSet
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   End Subroutine SectionRealH1semiDotProduct_2DScal

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1semiDotProduct_2DVect"
   Subroutine SectionRealH1semiDotProduct_2DVect(Mesh,U_sec,V_Sec,Elem,H1semiDotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec,V_Sec
      Type(Element2D_Vect), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1semiDotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: H1semiDotProductSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      H1semiDotProductSet = 0.0_Kr
      H1semiDotProduct = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealH1semiDotProductSet_2DVect(Mesh,setID(set),U_sec,V_Sec,Elem,H1semiDotProductSet,ierr);CHKERRQ(ierr)
         H1semiDotProduct = H1semiDotProduct + H1semiDotProductSet
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   End Subroutine SectionRealH1semiDotProduct_2DVect

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1semiDotProduct_2DElast"
   Subroutine SectionRealH1semiDotProduct_2DElast(Mesh,U_sec,V_Sec,Elem,H1semiDotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec,V_Sec
      Type(Element2D_Elast), Dimension(:), Pointer :: Elem
      PetscReal,Intent(OUT)                        :: H1semiDotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: H1semiDotProductSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      H1semiDotProductSet = 0.0_Kr
      H1semiDotProduct = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealH1semiDotProductSet_2DElast(Mesh,setID(set),U_sec,V_Sec,Elem,H1semiDotProductSet,ierr);CHKERRQ(ierr)
         H1semiDotProduct = H1semiDotProduct + H1semiDotProductSet
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   End Subroutine SectionRealH1semiDotProduct_2DElast

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1semiDotProduct_3DScal"
   Subroutine SectionRealH1semiDotProduct_3DScal(Mesh,U_sec,V_Sec,Elem,H1semiDotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec,V_Sec
      Type(Element3D_Scal), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1semiDotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: H1semiDotProductSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      H1semiDotProductSet = 0.0_Kr
      H1semiDotProduct = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealH1semiDotProductSet_3DScal(Mesh,setID(set),U_sec,V_Sec,Elem,H1semiDotProductSet,ierr);CHKERRQ(ierr)
         H1semiDotProduct = H1semiDotProduct + H1semiDotProductSet
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   End Subroutine SectionRealH1semiDotProduct_3DScal

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1semiDotProduct_3DVect"
   Subroutine SectionRealH1semiDotProduct_3DVect(Mesh,U_sec,V_Sec,Elem,H1semiDotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec,V_Sec
      Type(Element3D_Vect), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1semiDotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: H1semiDotProductSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      H1semiDotProductSet = 0.0_Kr
      H1semiDotProduct = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealH1semiDotProductSet_3DVect(Mesh,setID(set),U_sec,V_Sec,Elem,H1semiDotProductSet,ierr);CHKERRQ(ierr)
         H1semiDotProduct = H1semiDotProduct + H1semiDotProductSet
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   End Subroutine SectionRealH1semiDotProduct_3DVect

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1semiDotProduct_3DElast"
   Subroutine SectionRealH1semiDotProduct_3DElast(Mesh,U_sec,V_Sec,Elem,H1semiDotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec,V_Sec
      Type(Element3D_Elast), Dimension(:), Pointer :: Elem
      PetscReal,Intent(OUT)                        :: H1semiDotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: H1semiDotProductSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      H1semiDotProductSet = 0.0_Kr
      H1semiDotProduct = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealH1semiDotProductSet_3DElast(Mesh,setID(set),U_sec,V_Sec,Elem,H1semiDotProductSet,ierr);CHKERRQ(ierr)
         H1semiDotProduct = H1semiDotProduct + H1semiDotProductSet
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   End Subroutine SectionRealH1semiDotProduct_3DElast

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1semiNormSet_2DScal"
   Subroutine SectionRealH1semiNormSet_2DScal(Mesh,setID,U_sec,Elem,H1semiNormSet,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec
      Type(Element2D_Scal), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1semiNormSet
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Call SectionRealH1semiDotProductSet_2DScal(Mesh,setID,U_sec,U_Sec,Elem,H1semiNormSet,ierr);CHKERRQ(ierr)
      H1semiNormSet = sqrt(H1semiNormSet)
   End Subroutine SectionRealH1semiNormSet_2DScal
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealH1semiNormSet_2DVect"
   Subroutine SectionRealH1semiNormSet_2DVect(Mesh,setID,U_sec,Elem,H1semiNormSet,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec
      Type(Element2D_Vect), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1semiNormSet
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Call SectionRealH1semiDotProductSet_2DVect(Mesh,setID,U_sec,U_Sec,Elem,H1semiNormSet,ierr);CHKERRQ(ierr)
      H1semiNormSet = sqrt(H1semiNormSet)
   End Subroutine SectionRealH1semiNormSet_2DVect
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealH1semiNormSet_2DElast"
   Subroutine SectionRealH1semiNormSet_2DElast(Mesh,setID,U_sec,Elem,H1semiNormSet,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec
      Type(Element2D_Elast), Dimension(:), Pointer :: Elem
      PetscReal,Intent(OUT)                        :: H1semiNormSet
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Call SectionRealH1semiDotProductSet_2DElast(Mesh,setID,U_sec,U_Sec,Elem,H1semiNormSet,ierr);CHKERRQ(ierr)
      H1semiNormSet = sqrt(H1semiNormSet)
   End Subroutine SectionRealH1semiNormSet_2DElast
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealH1semiNormSet_3DScal"
   Subroutine SectionRealH1semiNormSet_3DScal(Mesh,setID,U_sec,Elem,H1semiNormSet,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec
      Type(Element3D_Scal), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1semiNormSet
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Call SectionRealL3DotProductSet_3DScal(Mesh,setID,U_sec,U_Sec,Elem,H1semiNormSet,ierr);CHKERRQ(ierr)
      H1semiNormSet = sqrt(H1semiNormSet)
   End Subroutine SectionRealH1semiNormSet_3DScal
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealH1semiNormSet_3DVect"
   Subroutine SectionRealH1semiNormSet_3DVect(Mesh,setID,U_sec,Elem,H1semiNormSet,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec
      Type(Element3D_Vect), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1semiNormSet
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Call SectionRealL3DotProductSet_3DVect(Mesh,setID,U_sec,U_Sec,Elem,H1semiNormSet,ierr);CHKERRQ(ierr)
      H1semiNormSet = sqrt(H1semiNormSet)
   End Subroutine SectionRealH1semiNormSet_3DVect
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealH1semiNormSet_3DElast"
   Subroutine SectionRealH1semiNormSet_3DElast(Mesh,setID,U_sec,Elem,H1semiNormSet,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec
      Type(Element3D_Elast), Dimension(:), Pointer :: Elem
      PetscReal,Intent(OUT)                        :: H1semiNormSet
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Call SectionRealL3DotProductSet_3DElast(Mesh,setID,U_sec,U_Sec,Elem,H1semiNormSet,ierr);CHKERRQ(ierr)
      H1semiNormSet = sqrt(H1semiNormSet)
   End Subroutine SectionRealH1semiNormSet_3DElast
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealH1semiNorm_2DScal"
   Subroutine SectionRealH1semiNorm_2DScal(Mesh,U_sec,Elem,H1semiNorm,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec
      Type(Element2D_Scal), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1semiNorm
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: H1semiNormSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      H1semiNormSet = 0.0_Kr
      H1semiNorm = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealH1semiNormSet_2DScal(Mesh,setID(set),U_sec,Elem,H1semiNormSet,ierr);CHKERRQ(ierr)
         H1semiNorm = H1semiNorm + H1semiNormSet**2
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      H1semiNorm = sqrt(H1semiNorm)
   End Subroutine SectionRealH1semiNorm_2DScal
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealH1semiNorm_2DVect"
   Subroutine SectionRealH1semiNorm_2DVect(Mesh,U_sec,Elem,H1semiNorm,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec
      Type(Element2D_Vect), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1semiNorm
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: H1semiNormSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      H1semiNormSet = 0.0_Kr
      H1semiNorm = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealH1semiNormSet_2DVect(Mesh,setID(set),U_sec,Elem,H1semiNormSet,ierr);CHKERRQ(ierr)
         H1semiNorm = H1semiNorm + H1semiNormSet**2
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      H1semiNorm = sqrt(H1semiNorm)
   End Subroutine SectionRealH1semiNorm_2DVect
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealH1semiNorm_2DElast"
   Subroutine SectionRealH1semiNorm_2DElast(Mesh,U_sec,Elem,H1semiNorm,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec
      Type(Element2D_Elast), Dimension(:), Pointer :: Elem
      PetscReal,Intent(OUT)                        :: H1semiNorm
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: H1semiNormSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      H1semiNormSet = 0.0_Kr
      H1semiNorm = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealH1semiNormSet_2DElast(Mesh,setID(set),U_sec,Elem,H1semiNormSet,ierr);CHKERRQ(ierr)
         H1semiNorm = H1semiNorm + H1semiNormSet**2
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      H1semiNorm = sqrt(H1semiNorm)
   End Subroutine SectionRealH1semiNorm_2DElast
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealH1semiNorm_3DScal"
   Subroutine SectionRealH1semiNorm_3DScal(Mesh,U_sec,Elem,H1semiNorm,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec
      Type(Element3D_Scal), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1semiNorm
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: H1semiNormSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      H1semiNormSet = 0.0_Kr
      H1semiNorm = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealH1semiNormSet_3DScal(Mesh,setID(set),U_sec,Elem,H1semiNormSet,ierr);CHKERRQ(ierr)
         H1semiNorm = H1semiNorm + H1semiNormSet**2
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      H1semiNorm = sqrt(H1semiNorm)
   End Subroutine SectionRealH1semiNorm_3DScal
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealH1semiNorm_3DVect"
   Subroutine SectionRealH1semiNorm_3DVect(Mesh,U_sec,Elem,H1semiNorm,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec
      Type(Element3D_Vect), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1semiNorm
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: H1semiNormSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      H1semiNormSet = 0.0_Kr
      H1semiNorm = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealH1semiNormSet_3DVect(Mesh,setID(set),U_sec,Elem,H1semiNormSet,ierr);CHKERRQ(ierr)
         H1semiNorm = H1semiNorm + H1semiNormSet**2
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      H1semiNorm = sqrt(H1semiNorm)
   End Subroutine SectionRealH1semiNorm_3DVect
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealH1semiNorm_3DElast"
   Subroutine SectionRealH1semiNorm_3DElast(Mesh,U_sec,Elem,H1semiNorm,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec
      Type(Element3D_Elast), Dimension(:), Pointer :: Elem
      PetscReal,Intent(OUT)                        :: H1semiNorm
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: H1semiNormSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      H1semiNormSet = 0.0_Kr
      H1semiNorm = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealH1semiNormSet_3DElast(Mesh,setID(set),U_sec,Elem,H1semiNormSet,ierr);CHKERRQ(ierr)
         H1semiNorm = H1semiNorm + H1semiNormSet**2
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      H1semiNorm = sqrt(H1semiNorm)
   End Subroutine SectionRealH1semiNorm_3DElast

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1DotProductSet_2DScal"
   Subroutine SectionRealH1DotProductSet_2DScal(Mesh,setID,U_sec,V_Sec,Elem,H1DotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec,V_Sec
      Type(Element2D_Scal), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1DotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      PetscReal                                    :: L2DotProduct,H1semiDotProduct
      
      Call SectionRealL2DotProductSet_2DScal(Mesh,setID,U_sec,V_Sec,Elem,L2DotProduct,ierr)
      Call SectionRealH1semiDotProductSet_2DScal(Mesh,setID,U_sec,V_Sec,Elem,H1semiDotProduct,ierr)
      H1DotProduct = L2DotProduct + H1semiDotProduct
   End Subroutine SectionRealH1DotProductSet_2DScal

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1DotProductSet_2DVect"
   Subroutine SectionRealH1DotProductSet_2DVect(Mesh,setID,U_sec,V_Sec,Elem,H1DotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec,V_Sec
      Type(Element2D_Vect), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1DotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      PetscReal                                    :: L2DotProduct,H1semiDotProduct
      
      Call SectionRealL2DotProductSet_2DVect(Mesh,setID,U_sec,V_Sec,Elem,L2DotProduct,ierr)
      Call SectionRealH1semiDotProductSet_2DVect(Mesh,setID,U_sec,V_Sec,Elem,H1semiDotProduct,ierr)
      H1DotProduct = L2DotProduct + H1semiDotProduct
   End Subroutine SectionRealH1DotProductSet_2DVect

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1DotProductSet_2DElast"
   Subroutine SectionRealH1DotProductSet_2DElast(Mesh,setID,U_sec,V_Sec,Elem,H1DotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec,V_Sec
      Type(Element2D_Elast), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1DotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      PetscReal                                    :: L2DotProduct,H1semiDotProduct
      
      Call SectionRealL2DotProductSet_2DElast(Mesh,setID,U_sec,V_Sec,Elem,L2DotProduct,ierr)
      Call SectionRealH1semiDotProductSet_2DElast(Mesh,setID,U_sec,V_Sec,Elem,H1semiDotProduct,ierr)
      H1DotProduct = L2DotProduct + H1semiDotProduct
   End Subroutine SectionRealH1DotProductSet_2DElast

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1DotProductSet_3DScal"
   Subroutine SectionRealH1DotProductSet_3DScal(Mesh,setID,U_sec,V_Sec,Elem,H1DotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec,V_Sec
      Type(Element3D_Scal), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1DotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      PetscReal                                    :: L2DotProduct,H1semiDotProduct
      
      Call SectionRealL2DotProductSet_3DScal(Mesh,setID,U_sec,V_Sec,Elem,L2DotProduct,ierr)
      Call SectionRealH1semiDotProductSet_3DScal(Mesh,setID,U_sec,V_Sec,Elem,H1semiDotProduct,ierr)
      H1DotProduct = L2DotProduct + H1semiDotProduct
   End Subroutine SectionRealH1DotProductSet_3DScal

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1DotProductSet_3DVect"
   Subroutine SectionRealH1DotProductSet_3DVect(Mesh,setID,U_sec,V_Sec,Elem,H1DotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec,V_Sec
      Type(Element3D_Vect), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1DotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      PetscReal                                    :: L2DotProduct,H1semiDotProduct
      
      Call SectionRealL2DotProductSet_3DVect(Mesh,setID,U_sec,V_Sec,Elem,L2DotProduct,ierr)
      Call SectionRealH1semiDotProductSet_3DVect(Mesh,setID,U_sec,V_Sec,Elem,H1semiDotProduct,ierr)
      H1DotProduct = L2DotProduct + H1semiDotProduct
   End Subroutine SectionRealH1DotProductSet_3DVect

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1DotProductSet_3DElast"
   Subroutine SectionRealH1DotProductSet_3DElast(Mesh,setID,U_sec,V_Sec,Elem,H1DotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec,V_Sec
      Type(Element3D_Elast), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1DotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      PetscReal                                    :: L2DotProduct,H1semiDotProduct
      
      Call SectionRealL2DotProductSet_3DElast(Mesh,setID,U_sec,V_Sec,Elem,L2DotProduct,ierr)
      Call SectionRealH1semiDotProductSet_3DElast(Mesh,setID,U_sec,V_Sec,Elem,H1semiDotProduct,ierr)
      H1DotProduct = L2DotProduct + H1semiDotProduct
   End Subroutine SectionRealH1DotProductSet_3DElast

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1DotProductSet_2DScal"
   Subroutine SectionRealH1DotProduct_2DScal(Mesh,U_sec,V_Sec,Elem,H1DotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_Sec,V_Sec
      Type(Element2D_Scal), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1DotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      PetscReal                                    :: L2DotProduct,H1semiDotProduct
      
      Call SectionRealL2DotProduct_2DScal(Mesh,U_sec,V_Sec,Elem,L2DotProduct,ierr)
      Call SectionRealH1semiDotProduct_2DScal(Mesh,U_sec,V_Sec,Elem,H1semiDotProduct,ierr)
      H1DotProduct = L2DotProduct + H1semiDotProduct
   End Subroutine SectionRealH1DotProduct_2DScal

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1DotProduct_2DVect"
   Subroutine SectionRealH1DotProduct_2DVect(Mesh,U_sec,V_Sec,Elem,H1DotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_Sec,V_Sec
      Type(Element2D_Vect), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1DotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      PetscReal                                    :: L2DotProduct,H1semiDotProduct
      
      Call SectionRealL2DotProduct_2DVect(Mesh,U_sec,V_Sec,Elem,L2DotProduct,ierr)
      Call SectionRealH1semiDotProduct_2DVect(Mesh,U_sec,V_Sec,Elem,H1semiDotProduct,ierr)
      H1DotProduct = L2DotProduct + H1semiDotProduct
   End Subroutine SectionRealH1DotProduct_2DVect

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1DotProduct_2DElast"
   Subroutine SectionRealH1DotProduct_2DElast(Mesh,U_sec,V_Sec,Elem,H1DotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_Sec,V_Sec
      Type(Element2D_Elast), Dimension(:), Pointer :: Elem
      PetscReal,Intent(OUT)                        :: H1DotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      PetscReal                                    :: L2DotProduct,H1semiDotProduct
      
      Call SectionRealL2DotProduct_2DElast(Mesh,U_sec,V_Sec,Elem,L2DotProduct,ierr)
      Call SectionRealH1semiDotProduct_2DElast(Mesh,U_sec,V_Sec,Elem,H1semiDotProduct,ierr)
      H1DotProduct = L2DotProduct + H1semiDotProduct
   End Subroutine SectionRealH1DotProduct_2DElast

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1DotProduct_3DScal"
   Subroutine SectionRealH1DotProduct_3DScal(Mesh,U_sec,V_Sec,Elem,H1DotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_Sec,V_Sec
      Type(Element3D_Scal), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1DotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      PetscReal                                    :: L2DotProduct,H1semiDotProduct
      
      Call SectionRealL2DotProduct_3DScal(Mesh,U_sec,V_Sec,Elem,L2DotProduct,ierr)
      Call SectionRealH1semiDotProduct_3DScal(Mesh,U_sec,V_Sec,Elem,H1semiDotProduct,ierr)
      H1DotProduct = L2DotProduct + H1semiDotProduct
   End Subroutine SectionRealH1DotProduct_3DScal

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1DotProduct_3DVect"
   Subroutine SectionRealH1DotProduct_3DVect(Mesh,U_sec,V_Sec,Elem,H1DotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_Sec,V_Sec
      Type(Element3D_Vect), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1DotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      PetscReal                                    :: L2DotProduct,H1semiDotProduct
      
      Call SectionRealL2DotProduct_3DVect(Mesh,U_sec,V_Sec,Elem,L2DotProduct,ierr)
      Call SectionRealH1semiDotProduct_3DVect(Mesh,U_sec,V_Sec,Elem,H1semiDotProduct,ierr)
      H1DotProduct = L2DotProduct + H1semiDotProduct
   End Subroutine SectionRealH1DotProduct_3DVect

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1DotProduct_3DElast"
   Subroutine SectionRealH1DotProduct_3DElast(Mesh,U_sec,V_Sec,Elem,H1DotProduct,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_Sec,V_Sec
      Type(Element3D_Elast), Dimension(:), Pointer :: Elem
      PetscReal,Intent(OUT)                        :: H1DotProduct
      PetscErrorCode,intent(OUT)                   :: ierr
      
      PetscReal                                    :: L2DotProduct,H1semiDotProduct
      
      Call SectionRealL2DotProduct_3DElast(Mesh,U_sec,V_Sec,Elem,L2DotProduct,ierr)
      Call SectionRealH1semiDotProduct_3DElast(Mesh,U_sec,V_Sec,Elem,H1semiDotProduct,ierr)
      H1DotProduct = L2DotProduct + H1semiDotProduct
   End Subroutine SectionRealH1DotProduct_3DElast

#undef __FUNCT__
#define __FUNCT__ "SectionRealH1NormSet_2DScal"
   Subroutine SectionRealH1NormSet_2DScal(Mesh,setID,U_sec,Elem,H1NormSet,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec
      Type(Element2D_Scal), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1NormSet
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Call SectionRealH1DotProductSet_2DScal(Mesh,setID,U_sec,U_Sec,Elem,H1NormSet,ierr);CHKERRQ(ierr)
      H1NormSet = sqrt(H1NormSet)
   End Subroutine SectionRealH1NormSet_2DScal
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealH1NormSet_2DVect"
   Subroutine SectionRealH1NormSet_2DVect(Mesh,setID,U_sec,Elem,H1NormSet,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec
      Type(Element2D_Vect), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1NormSet
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Call SectionRealH1DotProductSet_2DVect(Mesh,setID,U_sec,U_Sec,Elem,H1NormSet,ierr);CHKERRQ(ierr)
      H1NormSet = sqrt(H1NormSet)
   End Subroutine SectionRealH1NormSet_2DVect
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealH1NormSet_2DElast"
   Subroutine SectionRealH1NormSet_2DElast(Mesh,setID,U_sec,Elem,H1NormSet,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec
      Type(Element2D_Elast), Dimension(:), Pointer :: Elem
      PetscReal,Intent(OUT)                        :: H1NormSet
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Call SectionRealH1DotProductSet_2DElast(Mesh,setID,U_sec,U_Sec,Elem,H1NormSet,ierr);CHKERRQ(ierr)
      H1NormSet = sqrt(H1NormSet)
   End Subroutine SectionRealH1NormSet_2DElast
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealH1NormSet_3DScal"
   Subroutine SectionRealH1NormSet_3DScal(Mesh,setID,U_sec,Elem,H1NormSet,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec
      Type(Element3D_Scal), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1NormSet
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Call SectionRealL3DotProductSet_3DScal(Mesh,setID,U_sec,U_Sec,Elem,H1NormSet,ierr);CHKERRQ(ierr)
      H1NormSet = sqrt(H1NormSet)
   End Subroutine SectionRealH1NormSet_3DScal
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealH1NormSet_3DVect"
   Subroutine SectionRealH1NormSet_3DVect(Mesh,setID,U_sec,Elem,H1NormSet,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec
      Type(Element3D_Vect), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1NormSet
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Call SectionRealL3DotProductSet_3DVect(Mesh,setID,U_sec,U_Sec,Elem,H1NormSet,ierr);CHKERRQ(ierr)
      H1NormSet = sqrt(H1NormSet)
   End Subroutine SectionRealH1NormSet_3DVect
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealH1NormSet_3DElast"
   Subroutine SectionRealH1NormSet_3DElast(Mesh,setID,U_sec,Elem,H1NormSet,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      PetscInt,Intent(IN)                          :: setID
      Type(SectionReal),Intent(IN)                 :: U_Sec
      Type(Element3D_Elast), Dimension(:), Pointer :: Elem
      PetscReal,Intent(OUT)                        :: H1NormSet
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Call SectionRealL3DotProductSet_3DElast(Mesh,setID,U_sec,U_Sec,Elem,H1NormSet,ierr);CHKERRQ(ierr)
      H1NormSet = sqrt(H1NormSet)
   End Subroutine SectionRealH1NormSet_3DElast
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealH1Norm_2DScal"
   Subroutine SectionRealH1Norm_2DScal(Mesh,U_sec,Elem,H1Norm,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec
      Type(Element2D_Scal), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1Norm
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: H1NormSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      H1NormSet = 0.0_Kr
      H1Norm = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealH1NormSet_2DScal(Mesh,setID(set),U_sec,Elem,H1NormSet,ierr);CHKERRQ(ierr)
         H1Norm = H1Norm + H1NormSet**2
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      H1Norm = sqrt(H1Norm)
   End Subroutine SectionRealH1Norm_2DScal
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealH1Norm_2DVect"
   Subroutine SectionRealH1Norm_2DVect(Mesh,U_sec,Elem,H1Norm,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec
      Type(Element2D_Vect), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1Norm
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: H1NormSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      H1NormSet = 0.0_Kr
      H1Norm = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealH1NormSet_2DVect(Mesh,setID(set),U_sec,Elem,H1NormSet,ierr);CHKERRQ(ierr)
         H1Norm = H1Norm + H1NormSet**2
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      H1Norm = sqrt(H1Norm)
   End Subroutine SectionRealH1Norm_2DVect
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealH1Norm_2DElast"
   Subroutine SectionRealH1Norm_2DElast(Mesh,U_sec,Elem,H1Norm,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec
      Type(Element2D_Elast), Dimension(:), Pointer :: Elem
      PetscReal,Intent(OUT)                        :: H1Norm
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: H1NormSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      H1NormSet = 0.0_Kr
      H1Norm = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealH1NormSet_2DElast(Mesh,setID(set),U_sec,Elem,H1NormSet,ierr);CHKERRQ(ierr)
         H1Norm = H1Norm + H1NormSet**2
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      H1Norm = sqrt(H1Norm)
   End Subroutine SectionRealH1Norm_2DElast
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealH1Norm_3DScal"
   Subroutine SectionRealH1Norm_3DScal(Mesh,U_sec,Elem,H1Norm,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec
      Type(Element3D_Scal), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1Norm
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: H1NormSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      H1NormSet = 0.0_Kr
      H1Norm = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealH1NormSet_3DScal(Mesh,setID(set),U_sec,Elem,H1NormSet,ierr);CHKERRQ(ierr)
         H1Norm = H1Norm + H1NormSet**2
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      H1Norm = sqrt(H1Norm)
   End Subroutine SectionRealH1Norm_3DScal
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealH1Norm_3DVect"
   Subroutine SectionRealH1Norm_3DVect(Mesh,U_sec,Elem,H1Norm,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec
      Type(Element3D_Vect), Dimension(:), Pointer  :: Elem
      PetscReal,Intent(OUT)                        :: H1Norm
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: H1NormSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      H1NormSet = 0.0_Kr
      H1Norm = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealH1NormSet_3DVect(Mesh,setID(set),U_sec,Elem,H1NormSet,ierr);CHKERRQ(ierr)
         H1Norm = H1Norm + H1NormSet**2
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      H1Norm = sqrt(H1Norm)
   End Subroutine SectionRealH1Norm_3DVect
   
#undef __FUNCT__
#define __FUNCT__ "SectionRealH1Norm_3DElast"
   Subroutine SectionRealH1Norm_3DElast(Mesh,U_sec,Elem,H1Norm,ierr)
      Type(DM),Intent(IN)                          :: Mesh
      Type(SectionReal),Intent(IN)                 :: U_sec
      Type(Element3D_Elast), Dimension(:), Pointer :: Elem
      PetscReal,Intent(OUT)                        :: H1Norm
      PetscErrorCode,intent(OUT)                   :: ierr
      
      Type(IS)                                     :: setIS,elemIS
      PetscInt,Dimension(:),Pointer                :: setID,elemID
      PetscInt                                     :: set
      PetscReal                                    :: H1NormSet

      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      H1NormSet = 0.0_Kr
      H1Norm = 0.0_Kr
      Do set = 1,size(setID)
         Call SectionRealH1NormSet_3DElast(Mesh,setID(set),U_sec,Elem,H1NormSet,ierr);CHKERRQ(ierr)
         H1Norm = H1Norm + H1NormSet**2
      End Do !set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      H1Norm = sqrt(H1Norm)
   End Subroutine SectionRealH1Norm_3DElast

End Module m_MEF_Norm