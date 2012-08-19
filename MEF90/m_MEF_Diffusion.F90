Module m_MEF_Diffusion
#include "finclude/petscdef.h"
   Use m_MEF_LinAlg
   Use m_MEF_Parameters
   Use m_MEF_Elements
   Use m_MEF_Sieve
   Use m_MEF_Utils
   Use petsc
   IMPLICIT NONE

   Private   
   Public :: MEF90_DiffusionBilinearFormSet
   Public :: MEF90_DiffusionOperatorSet
   Public :: MEF90_DiffusionOperatorAddTransientTermSet
   Public :: MEF90_DiffusionRHSSet
   Public :: MEF90_DiffusionEnergySet
   Public :: MEF90_DiffusionWorkSet
   
   Interface MEF90_DiffusionBilinearFormSet
      Module Procedure DiffusionBilinearFormSet_2D, DiffusionBilinearFormSet_3D
   End Interface MEF90_DiffusionBilinearFormSet
   
   Interface MEF90_DiffusionRHSSet
      Module procedure DiffusionRHSSet_2D,DiffusionRHSCellCstSet_2D,DiffusionRHSSet_3D,DiffusionRHSCellCstSet_3D
   End Interface MEF90_DiffusionRHSSet

   Interface MEF90_DiffusionOperatorSet
      Module procedure DiffusionOperatorSet_2D,DiffusionOperatorSet_3D
   End Interface MEF90_DiffusionOperatorSet

   Interface MEF90_DiffusionOperatorAddTransientTermSet
      Module procedure DiffusionOperatorAddTransientTermSet_2D,DiffusionOperatorAddTransientTermSet_3D
   End Interface MEF90_DiffusionOperatorAddTransientTermSet
   
   Interface MEF90_DiffusionEnergySet
      Module procedure DiffusionEnergySet_2D,DiffusionEnergySet_3D
   End Interface MEF90_DiffusionEnergySet

   Interface MEF90_DiffusionWorkSet
      Module procedure DiffusionWorkSet_2D,DiffusionWorkCellCstSet_2D,DiffusionWorkSet_3D,DiffusionWorkCellCstSet_3D
   End Interface MEF90_DiffusionWorkSet
   

!  Assembles all components required to solve a diffusion equation in the form
!
!       { -div[A\nabla v] + \lambda v = f in \Omega
!  (1)  {                           v = 0 on \partial \Omega_d
!       {                        Av.n = g on \partial \Omega_n
!
!  or equivalently to minimize
!
!  (2) E(v) := 1/2 \int A\nabla v \cdot v + \lambda v^2\, dx - \int_\Omega fv\, dx - \int_{\partial \Omega_n} gv\, dS

!  Nomenclature:
!     the "Bilinear Form" is K such that 
!        <Kv,w> := \int A\nabla v \cdot w + \lambda v w\, dx
!     the "Operator" is v -> G(v)  where  
!        <G(v),w> := \int A\nabla v \cdot w + \lambda v w\, dx
!     the "RHS" is F, where
!        <F,w> := \int_\Omega fw\, dx + \int_{\partial \Omega_n} gw\, dS
!     the "Energy" is E
!
!     when using a SNES to solve (1),
!        the "Function" is G (operator)
!        the "Jacobian" is K (bilinear form)
!        the "RHS" is F (RHS)
!        OR
!        the "Function" is G-F (operator - RHS)
!        the "Jacobian" is K (bilinear form)
!        the "RHS" is 0
!     when using a KSP to solve (1)
!        the "Matrix" is K (bilinear form)
!        the "RHS" is F (RHS)
!     when using TAO to solve (2):
!        the "Objective Function" is E (energy)
!        the "Gradient" is G-F (operator - RHS)
!        the "Hessian" is K (bilinear form)
!

Contains
#undef __FUNCT__
#define __FUNCT__ "DiffusionBilinearFormSet_2D"
   Subroutine DiffusionBilinearFormSet_2D(K,mesh,U,cellIS,A,lambda,elem,elemType,ierr)
      Type(Mat),Intent(IN)                         :: K 
      Type(DM),Intent(IN)                          :: mesh
      Type(SectionReal),Intent(IN)                 :: U        
      !!! only the layout part if U is used (like a PetscSection)
      Type(IS),Intent(IN)                          :: cellIS
      Type(MatS2D),Intent(IN)                      :: A
      PetscReal,Intent(IN)                         :: lambda
      Type(Element2D_Scal), Dimension(:), Pointer  :: elem
      Type(Element_Type),Intent(IN)                :: elemType
      PetscErrorCode,Intent(OUT)                   :: ierr
      
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscReal,Dimension(:,:),Pointer             :: MatElem
      PetscInt                                     :: iDoF1,iDoF2,iGauss
      PetscLogDouble                               :: flops
     
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(MatElem(elemType%numDof,elemType%numDof))
         Do cell = 1,size(cellID)   
            MatElem = 0.0_Kr
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               Do iDoF1 = 1,elemType%numDof
                  Do iDoF2 = 1,elemType%numDof
                     MatElem(iDoF2,iDoF1) = MatElem(iDoF2,iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                    (lambda * elem(cell)%BF(iDoF1,iGauss) * elem(cell)%BF(iDoF2,iGauss) + &
                                    ( (A * elem(cell)%Grad_BF(iDoF1,iGauss)) .DotP. elem(cell)%Grad_BF(iDoF2,iGauss)))
                  End Do
               End Do
            End Do
            Call DMmeshAssembleMatrix(K,mesh,U,cellID(cell),MatElem,ADD_VALUES,ierr);CHKERRQ(ierr)
         End Do
         flops = 5 * elemType%numDof**2 * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(MatElem)
      End If
   End Subroutine DiffusionBilinearFormSet_2D

#undef __FUNCT__
#define __FUNCT__ "DiffusionOperatorSet_2D"   
   Subroutine DiffusionOperatorSet_2D(G,mesh,V,cellIS,A,lambda,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                 :: G
      Type(DM),Intent(IN)                          :: mesh
      Type(SectionReal),Intent(IN)                 :: V
      Type(IS),Intent(IN)                          :: cellIS
      Type(MatS2D),Intent(IN)                      :: A
      PetscReal,Intent(IN)                         :: lambda
      Type(Element2D_Scal), Dimension(:), Pointer  :: elem
      Type(Element_Type),Intent(IN)                :: elemType
      PetscErrorCode,Intent(OUT)                   :: ierr
      
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscReal,Dimension(:),Pointer               :: Gloc
      PetscReal,Dimension(:),Pointer               :: Vloc
      PetscReal                                    :: Velem
      Type(Vect2D)                                 :: GradVelem
      PetscInt                                     :: iDoF1,iGauss
      PetscLogDouble                               :: flops
           
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(Vloc(elemType%numDof))
         Allocate(Gloc(elemType%numDof))
         Do cell = 1,size(cellID)      
            Gloc = 0.0_Kr
            Call SectionRealRestrictClosure(V,mesh,cellID(cell),elemType%numDof,Vloc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               Velem     = 0.0_Kr
               GradVElem = 0.0_Kr
               Do iDoF1 = 1,elemType%numDof
                  Velem     = Velem     + Vloc(iDof1) * elem(cell)%BF(iDoF1,iGauss)
                  GradVelem = GradVelem + Vloc(iDof1) * elem(cell)%Grad_BF(iDoF1,iGauss)
               End Do
               Do iDoF1 = 1,elemType%numDof
                  Gloc(iDoF1) = Gloc(iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                 (lambda * elem(cell)%BF(iDoF1,iGauss) * Velem + &
                                  ( (A * elem(cell)%Grad_BF(iDoF1,iGauss)) .DotP. GradVelem))
               End Do
            End Do
            Call SectionRealUpdateClosure(G,mesh,cellID(cell),Gloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do
      
         flops = 7 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(Gloc)
         DeAllocate(Vloc)
      End If
   End Subroutine DiffusionOperatorSet_2D

#undef __FUNCT__
#define __FUNCT__ "DiffusionOperatorAddTransientTermSet_2D"
!!!
!!!  
!!!  DiffusionOperatorAddTransientTermSet_2D:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine DiffusionOperatorAddTransientTermSet_2D(G,mesh,x,cellIS,alpha,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                 :: G
      Type(DM),Intent(IN)                          :: mesh
      Type(SectionReal),Intent(IN)                 :: x
      Type(IS),Intent(IN)                          :: cellIS
      PetscReal,Intent(IN)                         :: alpha
      Type(Element2D_Scal), Dimension(:), Pointer  :: elem
      Type(Element_Type),Intent(IN)                :: elemType
      PetscErrorCode,Intent(OUT)                   :: ierr

      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscReal,Dimension(:),Pointer               :: xloc,Gloc
      PetscReal                                    :: xelem
      PetscInt                                     :: iDoF1,iGauss
      PetscLogDouble                               :: flops
           
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(xloc(elemType%numDof))
         Allocate(Gloc(elemType%numDof))
         Do cell = 1,size(cellID)      
            Gloc = 0.0_Kr
            Call SectionRealRestrictClosure(x,mesh,cellID(cell),elemType%numDof,xloc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               xelem     = 0.0_Kr
               Do iDoF1 = 1,elemType%numDof
                  xelem = xelem + xloc(iDof1) * elem(cell)%BF(iDoF1,iGauss)
               End Do
               Do iDoF1 = 1,elemType%numDof
                  Gloc(iDoF1) = Gloc(iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                elem(cell)%BF(iDoF1,iGauss) * xelem
               End Do
            End Do
            Gloc = Gloc * alpha
            Call SectionRealUpdateClosure(G,mesh,cellID(cell),Gloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do
      
         flops = (5 * elemType%numDof * size(elem(1)%Gauss_C) + 1) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(xloc)
         DeAllocate(Gloc)
      End If
   End Subroutine DiffusionOperatorAddTransientTermSet_2D

#undef __FUNCT__
#define __FUNCT__ "DiffusionRHSSet_2D"
   Subroutine DiffusionRHSSet_2D(RHS,mesh,F,cellIS,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                 :: RHS
      Type(DM),Intent(IN)                          :: mesh
      Type(SectionReal),Intent(IN)                 :: F
      Type(IS),Intent(IN)                          :: cellIS
      Type(Element2D_Scal), Dimension(:), Pointer  :: elem
      Type(Element_Type),Intent(IN)                :: elemType
      PetscErrorCode,Intent(OUT)                   :: ierr
      
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscReal,Dimension(:),Pointer               :: Floc,RHSloc
      PetscReal                                    :: Felem
      PetscInt                                     :: iDoF1,iGauss
      PetscLogDouble                               :: flops
           
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(Floc(elemType%numDof))
         Allocate(RHSloc(elemType%numDof))
         Do cell = 1,size(cellID)      
            RHSloc = 0.0_Kr
            Call SectionRealRestrictClosure(F,mesh,cellID(cell),elemType%numDof,Floc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               Felem     = 0.0_Kr
               Do iDoF1 = 1,elemType%numDof
                  Felem = Felem + Floc(iDof1) * elem(cell)%BF(iDoF1,iGauss)
               End Do
               Do iDoF1 = 1,elemType%numDof
                  RHSloc(iDoF1) = RHSloc(iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                 elem(cell)%BF(iDoF1,iGauss) * Felem
               End Do
            End Do
            Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do
      
         flops = 5 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(RHSloc)
         DeAllocate(Floc)
      End If
   End Subroutine DiffusionRHSSet_2D

#undef __FUNCT__
#define __FUNCT__ "DiffusionRHSCellCstSet_2D"
   Subroutine DiffusionRHSCellCstSet_2D(RHS,mesh,F,cellIS,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                 :: RHS
      Type(DM),Intent(IN)                          :: mesh
      PetscReal,Intent(IN)                         :: F
      Type(IS),Intent(IN)                          :: cellIS
      Type(Element2D_Scal), Dimension(:), Pointer  :: elem
      Type(Element_Type),Intent(IN)                :: elemType
      PetscErrorCode,Intent(OUT)                   :: ierr
      
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscReal,Dimension(:),Pointer               :: RHSloc
      PetscInt                                     :: iDoF1,iGauss
      PetscLogDouble                               :: flops
           
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(RHSloc(elemType%numDof))
         Do cell = 1,size(cellID)      
            RHSloc = 0.0_Kr
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               Do iDoF1 = 1,elemType%numDof
                  RHSloc(iDoF1) = RHSloc(iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                 elem(cell)%BF(iDoF1,iGauss) * F
               End Do
            End Do
            Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do
      
         flops = 3 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(RHSloc)
      End If
   End Subroutine DiffusionRHSCellCstSet_2D

#undef __FUNCT__
#define __FUNCT__ "DiffusionEnergySet_2D"
   Subroutine DiffusionEnergySet_2D(energy,x,mesh,A,lambda,cellIS,elem,elemType,ierr)
      PetscReal,Intent(OUT)                        :: energy
      Type(SectionReal),Intent(IN)                 :: x
      Type(DM),Intent(IN)                          :: mesh
      Type(MatS2D),Intent(IN)                      :: A
      PetscReal,Intent(IN)                         :: lambda
      Type(IS),Intent(IN)                          :: cellIS
      Type(Element2D_Scal), Dimension(:), Pointer  :: elem
      Type(Element_Type),Intent(IN)                :: elemType
      PetscErrorCode,Intent(OUT)                   :: ierr

      Type(Vect2D)                                 :: stress,strain      
      PetscReal                                    :: xelem
      PetscReal,Dimension(:),Pointer               :: xloc
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscInt                                     :: iDoF1,iGauss
      PetscLogDouble                               :: flops
     
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(xloc(elemType%numDof))
         Do cell = 1,size(cellID)   
            Call SectionRealRestrictClosure(x,mesh,cellID(cell),elemType%numDof,xloc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               strain = 0.0_Kr   
               xelem  = 0.0_Kr
               Do iDoF1 = 1,elemType%numDof
                  strain = strain + elem(cell)%Grad_BF(iDoF1,iGauss) * xloc(iDoF1)
                  xelem = xelem + elem(cell)%BF(iDoF1,iGauss) * xloc(iDoF1)
               End Do
               stress = A * strain
               energy = energy + elem(cell)%Gauss_C(iGauss) * ( (stress .dotP. strain) + lambda * xelem **2) *.5_Kr
            End Do
         End Do
      
         flops = (2 * elemType%numDof + 6) * size(elem(1)%Gauss_C) * size(cellID) 
         
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(xloc)
      End If
   End Subroutine DiffusionEnergySet_2D

#undef __FUNCT__
#define __FUNCT__ "DiffusionWorkSet_2D"
   Subroutine DiffusionWorkSet_2D(work,x,mesh,F,cellIS,elem,elemType,ierr)
      PetscReal,Intent(OUT)                        :: work
      Type(SectionReal),Intent(IN)                 :: x
      Type(DM),Intent(IN)                          :: mesh
      Type(SectionReal),Intent(IN)                 :: F
      Type(IS),Intent(IN)                          :: cellIS
      Type(Element2D_Scal), Dimension(:), Pointer  :: elem
      Type(Element_Type),Intent(IN)                :: elemType
      PetscErrorCode,Intent(OUT)                   :: ierr

      PetscReal,Dimension(:),Pointer               :: xloc,floc
      PetscReal                                    :: xelem,felem
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscInt                                     :: iDoF1,iGauss
      PetscLogDouble                               :: flops
     
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(xloc(elemType%numDof))
         Allocate(floc(elemType%numDof))
         Do cell = 1,size(cellID)   
            Call SectionRealRestrictClosure(x,mesh,cellID(cell),elemType%numDof,xloc,ierr);CHKERRQ(ierr)
            Call SectionRealRestrictClosure(f,mesh,cellID(cell),elemType%numDof,floc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               xelem = 0.0_Kr
               felem = 0.0_Kr
               Do iDoF1 = 1,elemType%numDof
                  xelem = xelem + xloc(iDof1) * elem(cell)%BF(iDof1,iGauss)
                  felem = felem + floc(iDof1) * elem(cell)%BF(iDof1,iGauss)
               End Do
               work = work + elem(cell)%Gauss_C(iGauss) * felem * xelem
            End Do
         End Do
      
         flops = (4 * elemType%numDof + 3 )* size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(xloc)
         DeAllocate(floc)
      End If
   End Subroutine DiffusionWorkSet_2D

#undef __FUNCT__
#define __FUNCT__ "DiffusionWorkCellCstSet_2D"
   Subroutine DiffusionWorkCellCstSet_2D(work,x,mesh,F,cellIS,elem,elemType,ierr)
      PetscReal,Intent(OUT)                        :: work
      Type(SectionReal),Intent(IN)                 :: x
      Type(DM),Intent(IN)                          :: mesh
      PetscReal,Intent(IN)                         :: F
      Type(IS),Intent(IN)                          :: cellIS
      Type(Element2D_Scal), Dimension(:), Pointer  :: elem
      Type(Element_Type),Intent(IN)                :: elemType
      PetscErrorCode,Intent(OUT)                   :: ierr

      PetscReal,Dimension(:),Pointer               :: xloc
      PetscReal                                    :: xelem
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscInt                                     :: iDoF1,iGauss
      PetscLogDouble                               :: flops
     
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(xloc(elemType%numDof))
         Do cell = 1,size(cellID)   
            Call SectionRealRestrictClosure(x,mesh,cellID(cell),elemType%numDof,xloc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               xelem = 0.0_Kr
               Do iDoF1 = 1,elemType%numDof
                  xelem = xelem + xloc(iDof1) * elem(cell)%BF(iDof1,iGauss)
               End Do
               work = work + elem(cell)%Gauss_C(iGauss) * F * xelem
            End Do
         End Do
      
         flops = (2 * elemType%numDof + 3 )* size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(xloc)
      End If
   End Subroutine DiffusionWorkCellCstSet_2D

#undef __FUNCT__
#define __FUNCT__ "DiffusionBilinearFormSet_3D"
   Subroutine DiffusionBilinearFormSet_3D(K,mesh,U,cellIS,A,lambda,elem,elemType,ierr)
      Type(Mat),Intent(IN)                         :: K
      Type(DM),Intent(IN)                          :: mesh
      Type(SectionReal),Intent(IN)                 :: U
      Type(IS),Intent(IN)                          :: cellIS
      Type(MatS3D),Intent(IN)                      :: A
      PetscReal,Intent(IN)                         :: lambda
      Type(Element3D_Scal), Dimension(:), Pointer  :: elem
      Type(Element_Type),Intent(IN)                :: elemType
      PetscErrorCode,Intent(OUT)                   :: ierr
      
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscReal,Dimension(:,:),Pointer             :: MatElem
      PetscInt                                     :: iDoF1,iDoF2,iGauss
      PetscLogDouble                               :: flops
     
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(MatElem(elemType%numDof,elemType%numDof))
         Do cell = 1,size(cellID)      
            MatElem = 0.0_Kr
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               Do iDoF1 = 1,elemType%numDof
                  Do iDoF2 = 1,elemType%numDof
                     MatElem(iDoF2,iDoF1) = MatElem(iDoF2,iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                    (lambda * elem(cell)%BF(iDoF1,iGauss) * elem(cell)%BF(iDoF2,iGauss) + &
                                    ( (A * elem(cell)%Grad_BF(iDoF1,iGauss)) .DotP. elem(cell)%Grad_BF(iDoF2,iGauss)))
                  End Do
               End Do
            End Do
            Call DMmeshAssembleMatrix(K,mesh,U,cellID(cell),MatElem,ADD_VALUES,ierr);CHKERRQ(ierr)
         End Do
         flops =  5 * elemType%numDof**2 * size(elem(1)%Gauss_C) * size(cellID)   
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(MatElem)
      End If
   End Subroutine DiffusionBilinearFormSet_3D

#undef __FUNCT__
#define __FUNCT__ "DiffusionOperatorSet_3D"   
   Subroutine DiffusionOperatorSet_3D(G,mesh,V,cellIS,A,lambda,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                 :: G
      Type(DM),Intent(IN)                          :: mesh
      Type(SectionReal),Intent(IN)                 :: V
      Type(IS),Intent(IN)                          :: cellIS
      Type(MatS3D),Intent(IN)                      :: A
      PetscReal,Intent(IN)                         :: lambda
      Type(Element3D_Scal), Dimension(:), Pointer  :: elem
      Type(Element_Type),Intent(IN)                :: elemType
      PetscErrorCode,Intent(OUT)                   :: ierr
      
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscReal,Dimension(:),Pointer               :: Gloc
      PetscReal,Dimension(:),Pointer               :: Vloc
      PetscReal                                    :: Velem
      Type(Vect3D)                                 :: GradVelem
      PetscInt                                     :: iDoF1,iGauss
      PetscLogDouble                               :: flops
     
      
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(Vloc(elemType%numDof))
         Allocate(Gloc(elemType%numDof))
         Do cell = 1,size(cellID)      
            Gloc = 0.0_Kr
            Call SectionRealRestrictClosure(V,mesh,cellID(cell),elemType%numDof,Vloc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               Velem     = 0.0_Kr
               GradVElem = 0.0_Kr
               Do iDoF1 = 1,elemType%numDof
                  Velem     = Velem     + Vloc(iDof1) * elem(cell)%BF(iDoF1,iGauss)
                  GradVelem = GradVelem + Vloc(iDof1) * elem(cell)%Grad_BF(iDoF1,iGauss)
               End Do
               Do iDoF1 = 1,elemType%numDof
                  Gloc(iDoF1) = Gloc(iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                 (lambda * elem(cell)%BF(iDoF1,iGauss) * Velem + &
                                  ( (A * elem(cell)%Grad_BF(iDoF1,iGauss)) .DotP. GradVelem))
               End Do
            End Do
            Call SectionRealUpdateClosure(G,mesh,cellID(cell),Gloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do
      
         flops = 7 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(Gloc)
         DeAllocate(Vloc)
      End If
   End Subroutine DiffusionOperatorSet_3D

#undef __FUNCT__
#define __FUNCT__ "DiffusionOperatorAddTransientTermSet_3D"
!!!
!!!  
!!!  DiffusionOperatorAddTransientTermSet_3D:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine DiffusionOperatorAddTransientTermSet_3D(G,mesh,x,cellIS,alpha,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                 :: G
      Type(DM),Intent(IN)                          :: mesh
      Type(SectionReal),Intent(IN)                 :: x
      Type(IS),Intent(IN)                          :: cellIS
      PetscReal,Intent(IN)                         :: alpha
      Type(Element3D_Scal), Dimension(:), Pointer  :: elem
      Type(Element_Type),Intent(IN)                :: elemType
      PetscErrorCode,Intent(OUT)                   :: ierr

      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscReal,Dimension(:),Pointer               :: xloc,Gloc
      PetscReal                                    :: xelem
      PetscInt                                     :: iDoF1,iGauss
      PetscLogDouble                               :: flops
           
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(xloc(elemType%numDof))
         Allocate(Gloc(elemType%numDof))
         Do cell = 1,size(cellID)      
            Gloc = 0.0_Kr
            Call SectionRealRestrictClosure(x,mesh,cellID(cell),elemType%numDof,xloc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               xelem     = 0.0_Kr
               Do iDoF1 = 1,elemType%numDof
                  xelem = xelem + xloc(iDof1) * elem(cell)%BF(iDoF1,iGauss)
               End Do
               Do iDoF1 = 1,elemType%numDof
                  Gloc(iDoF1) = Gloc(iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                elem(cell)%BF(iDoF1,iGauss) * xelem
               End Do
            End Do
            Gloc = Gloc * alpha
            Call SectionRealUpdateClosure(G,mesh,cellID(cell),Gloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do
      
         flops = (5 * elemType%numDof * size(elem(1)%Gauss_C) + 1) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(xloc)
         DeAllocate(Gloc)
      End If
   End Subroutine DiffusionOperatorAddTransientTermSet_3D


#undef __FUNCT__
#define __FUNCT__ "DiffusionRHSSet_3D"
   Subroutine DiffusionRHSSet_3D(RHS,mesh,F,cellIS,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                 :: RHS
      Type(DM),Intent(IN)                          :: mesh
      Type(SectionReal),Intent(IN)                 :: F
      Type(IS),Intent(IN)                          :: cellIS
      Type(Element3D_Scal), Dimension(:), Pointer  :: elem
      Type(Element_Type),Intent(IN)                :: elemType
      PetscErrorCode,Intent(OUT)                   :: ierr
      
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscReal,Dimension(:),Pointer               :: Floc,RHSloc
      PetscReal                                    :: Felem
      PetscInt                                     :: iDoF1,iGauss
      PetscLogDouble                               :: flops
     
      
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(Floc(elemType%numDof))
         Allocate(RHSloc(elemType%numDof))
         Do cell = 1,size(cellID)      
            RHSloc = 0.0_Kr
            Call SectionRealRestrictClosure(F,mesh,cellID(cell),elemType%numDof,Floc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               Felem     = 0.0_Kr
               Do iDoF1 = 1,elemType%numDof
                  Felem = Felem + Floc(iDof1) * elem(cell)%BF(iDoF1,iGauss)
               End Do
               Do iDoF1 = 1,elemType%numDof
                  RHSloc(iDoF1) = RHSloc(iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                 elem(cell)%BF(iDoF1,iGauss) * Felem
               End Do
            End Do
            Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do
      
         flops = 5 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(RHSloc)
         DeAllocate(Floc)
      End If
   End Subroutine DiffusionRHSSet_3D

#undef __FUNCT__
#define __FUNCT__ "DiffusionRHSCellCstSet_3D"
   Subroutine DiffusionRHSCellCstSet_3D(RHS,mesh,F,cellIS,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                 :: RHS
      Type(DM),Intent(IN)                          :: mesh
      PetscReal,Intent(IN)                         :: F
      Type(IS),Intent(IN)                          :: cellIS
      Type(Element3D_Scal), Dimension(:), Pointer  :: elem
      Type(Element_Type),Intent(IN)                :: elemType
      PetscErrorCode,Intent(OUT)                   :: ierr
      
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscReal,Dimension(:),Pointer               :: RHSloc
      PetscInt                                     :: iDoF1,iGauss
      PetscLogDouble                               :: flops
     
      
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(RHSloc(elemType%numDof))
         Do cell = 1,size(cellID)      
            RHSloc = 0.0_Kr
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               Do iDoF1 = 1,elemType%numDof
                  RHSloc(iDoF1) = RHSloc(iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                 elem(cell)%BF(iDoF1,iGauss) * F
               End Do
            End Do
            Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do
      
         flops = 3 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(RHSloc)
      EndIf
   End Subroutine DiffusionRHSCellCstSet_3D

#undef __FUNCT__
#define __FUNCT__ "DiffusionEnergySet_3D"
   Subroutine DiffusionEnergySet_3D(energy,x,mesh,A,lambda,cellIS,elem,elemType,ierr)
      PetscReal,Intent(OUT)                        :: energy
      Type(SectionReal),Intent(IN)                 :: x
      Type(DM),Intent(IN)                          :: mesh
      Type(MatS3D),Intent(IN)                      :: A
      PetscReal,Intent(IN)                         :: lambda
      Type(IS),Intent(IN)                          :: cellIS
      Type(Element3D_Scal), Dimension(:), Pointer  :: elem
      Type(Element_Type),Intent(IN)                :: elemType
      PetscErrorCode,Intent(OUT)                   :: ierr

      Type(Vect3D)                                 :: stress,strain      
      PetscReal                                    :: xelem
      PetscReal,Dimension(:),Pointer               :: xloc
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscInt                                     :: iDoF1,iGauss
      PetscLogDouble                               :: flops
     
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(xloc(elemType%numDof))
         Do cell = 1,size(cellID)   
            Call SectionRealRestrictClosure(x,mesh,cellID(cell),elemType%numDof,xloc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               strain = 0.0_Kr   
               xelem  = 0.0_Kr
               Do iDoF1 = 1,elemType%numDof
                  strain = strain + elem(cell)%Grad_BF(iDoF1,iGauss) * xloc(iDoF1)
                  xelem = xelem + elem(cell)%BF(iDoF1,iGauss) * xloc(iDoF1)
               End Do
               stress = A * strain
               energy = energy + elem(cell)%Gauss_C(iGauss) * ( (stress .dotP. strain) + lambda * xelem **2) *.5_Kr
            End Do
         End Do
      
         flops = (2 * elemType%numDof + 6) * size(elem(1)%Gauss_C) * size(cellID) 
         
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(xloc)
      End If
   End Subroutine DiffusionEnergySet_3D

#undef __FUNCT__
#define __FUNCT__ "DiffusionWorkSet_3D"
   Subroutine DiffusionWorkSet_3D(work,x,mesh,F,cellIS,elem,elemType,ierr)
      PetscReal,Intent(OUT)                        :: work
      Type(SectionReal),Intent(IN)                 :: x
      Type(DM),Intent(IN)                          :: mesh
      Type(SectionReal),Intent(IN)                 :: F
      Type(IS),Intent(IN)                          :: cellIS
      Type(Element3D_Scal), Dimension(:), Pointer  :: elem
      Type(Element_Type),Intent(IN)                :: elemType
      PetscErrorCode,Intent(OUT)                   :: ierr

      PetscReal,Dimension(:),Pointer               :: xloc,floc
      PetscReal                                    :: xelem,felem
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscInt                                     :: iDoF1,iGauss
      PetscLogDouble                               :: flops
     
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(xloc(elemType%numDof))
         Allocate(floc(elemType%numDof))
         Do cell = 1,size(cellID)   
            Call SectionRealRestrictClosure(x,mesh,cellID(cell),elemType%numDof,xloc,ierr);CHKERRQ(ierr)
            Call SectionRealRestrictClosure(f,mesh,cellID(cell),elemType%numDof,floc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               xelem = 0.0_Kr
               felem = 0.0_Kr
               Do iDoF1 = 1,elemType%numDof
                  xelem = xelem + xloc(iDof1) * elem(cell)%BF(iDof1,iGauss)
                  felem = felem + floc(iDof1) * elem(cell)%BF(iDof1,iGauss)
               End Do
               work = work + elem(cell)%Gauss_C(iGauss) * felem * xelem
            End Do
         End Do
      
         flops = (4 * elemType%numDof + 3 )* size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(xloc)
         DeAllocate(floc)
      End If
   End Subroutine DiffusionWorkSet_3D


#undef __FUNCT__
#define __FUNCT__ "DiffusionWorkCellCstSet_3D"
   Subroutine DiffusionWorkCellCstSet_3D(work,x,mesh,F,cellIS,elem,elemType,ierr)
      PetscReal,Intent(OUT)                        :: work
      Type(SectionReal),Intent(IN)                 :: x
      Type(DM),Intent(IN)                          :: mesh
      PetscReal,Intent(IN)                         :: F
      Type(IS),Intent(IN)                          :: cellIS
      Type(Element3D_Scal), Dimension(:), Pointer  :: elem
      Type(Element_Type),Intent(IN)                :: elemType
      PetscErrorCode,Intent(OUT)                   :: ierr

      PetscReal,Dimension(:),Pointer               :: xloc
      PetscReal                                    :: xelem
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscInt                                     :: iDoF1,iGauss
      PetscLogDouble                               :: flops
     
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(xloc(elemType%numDof))
         Do cell = 1,size(cellID)   
            Call SectionRealRestrictClosure(x,mesh,cellID(cell),elemType%numDof,xloc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               xelem = 0.0_Kr
               Do iDoF1 = 1,elemType%numDof
                  xelem = xelem + xloc(iDof1) * elem(cell)%BF(iDof1,iGauss)
               End Do
               work = work + elem(cell)%Gauss_C(iGauss) * F * xelem
            End Do
         End Do
      
         flops = (2 * elemType%numDof + 3 )* size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(xloc)
      End If
   End Subroutine DiffusionWorkCellCstSet_3D
End Module m_MEF_Diffusion
