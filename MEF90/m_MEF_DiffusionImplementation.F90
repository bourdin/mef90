#include "mef90.inc"
Module MEF90_APPEND(m_MEF_DiffusionImplementation_,MEF90_DIM)D
#include "finclude/petscdef.h"
   Use m_MEF_LinAlg
   Use m_MEF_Parameters
   Use m_MEF_Elements
   Use m_MEF_Sieve
   Use m_MEF_Utils
   Use petsc
   IMPLICIT NONE

   Private   
   Public :: DiffusionBilinearFormSet
   Public :: DiffusionEnergySet
   Public :: DiffusionOperatorAddTransientTermSet
   Public :: DiffusionOperatorSet
   Public :: DiffusionRHSSetVertex
   Public :: DiffusionRHSSetCell
   Public :: DiffusionRHSSetCst
   Public :: DiffusionWorkSetVertex
   Public :: DiffusionWorkSetCell
   Public :: DiffusionWorkSetCst

   

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
#define __FUNCT__ "DiffusionBilinearFormSet"
   Subroutine DiffusionBilinearFormSet(K,mesh,cellIS,A,lambda,elem,elemType,ierr)
      Type(Mat),Intent(IN)                               :: K 
      Type(DM),Intent(IN)                                :: mesh
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_MATS),Intent(IN)                        :: A
      PetscReal,Intent(IN)                               :: lambda
      Type(MEF90_ELEMENTTYPE), Dimension(:), Pointer     :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscReal,Dimension(:,:),Pointer                   :: MatElem
      PetscInt                                           :: iDoF1,iDoF2,iGauss
      Type(SectionReal)                                  :: defaultSection
      PetscLogDouble                                     :: flops
     
      Call DMMeshGetSectionReal(mesh,'default',defaultSection,ierr);CHKERRQ(ierr)
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
            Call DMmeshAssembleMatrix(K,mesh,defaultSection,cellID(cell),MatElem,ADD_VALUES,ierr);CHKERRQ(ierr)
         End Do
         flops = 5 * elemType%numDof**2 * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(MatElem)
      End If
   End Subroutine DiffusionBilinearFormSet

#undef __FUNCT__
#define __FUNCT__ "DiffusionOperatorSet"   
   Subroutine DiffusionOperatorSet(G,mesh,V,cellIS,A,lambda,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                       :: G
      Type(DM),Intent(IN)                                :: mesh
      Type(SectionReal),Intent(IN)                       :: V
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_MATS),Intent(IN)                        :: A
      PetscReal,Intent(IN)                               :: lambda
      Type(MEF90_ELEMENTTYPE), Dimension(:), Pointer     :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscReal,Dimension(:),Pointer                     :: Gloc
      PetscReal,Dimension(:),Pointer                     :: Vloc
      PetscReal                                          :: Velem
      Type(MEF90_VECT)                                   :: GradVelem
      PetscInt                                           :: iDoF1,iGauss
      PetscLogDouble                                     :: flops
           
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
   End Subroutine DiffusionOperatorSet

#undef __FUNCT__
#define __FUNCT__ "DiffusionOperatorAddTransientTermSet"
!!!
!!!  
!!!  DiffusionOperatorAddTransientTermSet:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine DiffusionOperatorAddTransientTermSet(G,mesh,x,cellIS,alpha,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                       :: G
      Type(DM),Intent(IN)                                :: mesh
      Type(SectionReal),Intent(IN)                       :: x
      Type(IS),Intent(IN)                                :: cellIS
      PetscReal,Intent(IN)                               :: alpha
      Type(MEF90_ELEMENTTYPE), Dimension(:), Pointer     :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr

      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscReal,Dimension(:),Pointer                     :: xloc,Gloc
      PetscReal                                          :: xelem
      PetscInt                                           :: iDoF1,iGauss
      PetscLogDouble                                     :: flops
           
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
   End Subroutine DiffusionOperatorAddTransientTermSet

#undef __FUNCT__
#define __FUNCT__ "DiffusionRHSSetVertex"
   Subroutine DiffusionRHSSetVertex(RHS,mesh,F,cellIS,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                       :: RHS
      Type(DM),Intent(IN)                                :: mesh
      Type(SectionReal),Intent(IN)                       :: F
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENTTYPE), Dimension(:), Pointer     :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscReal,Dimension(:),Pointer                     :: Floc,RHSloc
      PetscReal                                          :: Felem
      PetscInt                                           :: iDoF1,iGauss
      PetscLogDouble                                     :: flops
           
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
   End Subroutine DiffusionRHSSetVertex

#undef __FUNCT__
#define __FUNCT__ "DiffusionRHSSetCell"
   Subroutine DiffusionRHSSetCell(RHS,mesh,F,cellIS,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                       :: RHS
      Type(DM),Intent(IN)                                :: mesh
      Type(SectionReal),Intent(IN)                       :: F
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENTTYPE), Dimension(:), Pointer     :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscReal,Dimension(:),Pointer                     :: Floc,RHSloc
      PetscReal                                          :: Felem
      PetscInt                                           :: iDoF1,iGauss
      PetscLogDouble                                     :: flops
           
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         !Allocate(Floc(1))
         Allocate(RHSloc(elemType%numDof))
         Do cell = 1,size(cellID)      
            RHSloc = 0.0_Kr
            Call SectionRealRestrict(F,cellID(cell),Floc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               Do iDoF1 = 1,elemType%numDof
                  RHSloc(iDoF1) = RHSloc(iDoF1) + (elem(cell)%Gauss_C(iGauss) * &
                                 elem(cell)%BF(iDoF1,iGauss) * Floc(1))
               End Do
            End Do
            Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
            Call SectionRealRestore(F,cellID(cell),Floc,ierr);CHKERRQ(ierr)
         End Do
      
         flops = 3 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(RHSloc)
      End If
   End Subroutine DiffusionRHSSetCell

#undef __FUNCT__
#define __FUNCT__ "DiffusionRHSSetCst"
   Subroutine DiffusionRHSSetCst(RHS,mesh,F,cellIS,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                       :: RHS
      Type(DM),Intent(IN)                                :: mesh
      PetscReal,Intent(IN)                               :: F
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENTTYPE), Dimension(:), Pointer     :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscReal,Dimension(:),Pointer                     :: RHSloc
      PetscInt                                           :: iDoF1,iGauss
      PetscLogDouble                                     :: flops
           
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
   End Subroutine DiffusionRHSSetCst

#undef __FUNCT__
#define __FUNCT__ "DiffusionEnergySet"
   Subroutine DiffusionEnergySet(energy,x,mesh,A,lambda,cellIS,elem,elemType,ierr)
      PetscReal,Intent(OUT)                              :: energy
      Type(SectionReal),Intent(IN)                       :: x
      Type(DM),Intent(IN)                                :: mesh
      Type(MEF90_MATS),Intent(IN)                        :: A
      PetscReal,Intent(IN)                               :: lambda
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENTTYPE), Dimension(:), Pointer     :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(MEF90_VECT)                                   :: stress,strain      
      PetscReal                                          :: xelem
      PetscReal,Dimension(:),Pointer                     :: xloc
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscInt                                           :: iDoF1,iGauss
      PetscLogDouble                                     :: flops
     
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
   End Subroutine DiffusionEnergySet

#undef __FUNCT__
#define __FUNCT__ "DiffusionWorkSetVertex"
   Subroutine DiffusionWorkSetVertex(work,x,mesh,F,cellIS,elem,elemType,ierr)
      PetscReal,Intent(OUT)                              :: work
      Type(SectionReal),Intent(IN)                       :: x
      Type(DM),Intent(IN)                                :: mesh
      Type(SectionReal),Intent(IN)                       :: F
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENTTYPE), Dimension(:), Pointer     :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr

      PetscReal,Dimension(:),Pointer                     :: xloc,floc
      PetscReal                                          :: xelem,felem
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscInt                                           :: iDoF1,iGauss
      PetscLogDouble                                     :: flops
     
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
   End Subroutine DiffusionWorkSetVertex

#undef __FUNCT__
#define __FUNCT__ "DiffusionWorkSetCell"
   Subroutine DiffusionWorkSetCell(work,x,mesh,F,cellIS,elem,elemType,ierr)
      PetscReal,Intent(OUT)                              :: work
      Type(SectionReal),Intent(IN)                       :: x
      Type(DM),Intent(IN)                                :: mesh
      Type(SectionReal),Intent(IN)                       :: F
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENTTYPE), Dimension(:), Pointer     :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr

      PetscReal,Dimension(:),Pointer                     :: xloc,floc
      PetscReal                                          :: xelem
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscInt                                           :: iDoF1,iGauss
      PetscLogDouble                                     :: flops
     
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(xloc(elemType%numDof))
         Allocate(floc(1))
         Do cell = 1,size(cellID)   
            Call SectionRealRestrictClosure(x,mesh,cellID(cell),elemType%numDof,xloc,ierr);CHKERRQ(ierr)
            Call SectionRealRestrict(f,cellID(cell),floc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               xelem = 0.0_Kr
               Do iDoF1 = 1,elemType%numDof
                  xelem = xelem + xloc(iDof1) * elem(cell)%BF(iDof1,iGauss)
               End Do
               work = work + elem(cell)%Gauss_C(iGauss) * xelem * Floc(1)
            End Do
         End Do
      
         flops = (4 * elemType%numDof + 3 )* size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(xloc)
         DeAllocate(floc)
      End If
   End Subroutine DiffusionWorkSetCell

#undef __FUNCT__
#define __FUNCT__ "DiffusionWorkSetCst"
   Subroutine DiffusionWorkSetCst(work,x,mesh,F,cellIS,elem,elemType,ierr)
      PetscReal,Intent(OUT)                              :: work
      Type(SectionReal),Intent(IN)                       :: x
      Type(DM),Intent(IN)                                :: mesh
      PetscReal,Intent(IN)                               :: F
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENTTYPE), Dimension(:), Pointer     :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr

      PetscReal,Dimension(:),Pointer                     :: xloc
      PetscReal                                          :: xelem
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscInt                                           :: iDoF1,iGauss
      PetscLogDouble                                     :: flops
     
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
   End Subroutine DiffusionWorkSetCst

End Module MEF90_APPEND(m_MEF_DiffusionImplementation_,MEF90_DIM)D