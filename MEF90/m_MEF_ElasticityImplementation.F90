#include "mef90.inc"
Module MEF90_APPEND(m_MEF_ElasticityImplementation_,MEF90_DIM)D
#include "finclude/petscdef.h"
   Use m_MEF_LinAlg
   Use m_MEF_Parameters
   Use m_MEF_Elements
   Use m_MEF_Sieve
   Use m_MEF_Utils
   Use petsc
   IMPLICIT NONE

   Private   
   Public :: ElasticityBilinearFormSet
   Public :: ElasticityEnergySet
   Public :: ElasticityOperatorAddTransientTermSet
   Public :: ElasticityOperatorSet
   Public :: ElasticityRHSSetVertex
   Public :: ElasticityRHSSetCell
   Public :: ElasticityRHSSetCst
   Public :: ElasticityWorkSetVertex
   Public :: ElasticityWorkSetCell
   Public :: ElasticityWorkSetCst

   

!  Assembles all components required to solve a Elasticity equation in the form
!
!       { -div[A\e(v)] = f   in \Omega
!  (1)  {            v = v_0 on \partial \Omega_d
!       {         Av.n = g   on \partial \Omega_n
!
!  or equivalently to minimize
!
!  (2) E(v) := 1/2 \int A\e(v):\e(v) \, dx - \int_\Omega fv\, dx - \int_{\partial \Omega_n} gv\, dS

!  Nomenclature:
!     the "Bilinear Form" is K such that 
!        <Kv,w> := \int A\e(v) : e(w) \, dx
!     the "Operator" is v -> G(v)  where  
!        <G(v),w> := \int A\e(v):\e(w)\, dx
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
#define __FUNCT__ "ElasticityBilinearFormSet"
   Subroutine ElasticityBilinearFormSet(K,mesh,cellIS,A,elem,elemType,ierr)
      Type(Mat),Intent(IN)                               :: K 
      Type(DM),Intent(IN)                                :: mesh
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_TENS4OS),Intent(IN)                     :: A
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
                                    ((A * elem(cell)%GradS_BF(iDoF1,iGauss)) .DotP. elem(cell)%GradS_BF(iDoF2,iGauss))
                  End Do
               End Do
            End Do
            Call DMmeshAssembleMatrix(K,mesh,defaultSection,cellID(cell),MatElem,ADD_VALUES,ierr);CHKERRQ(ierr)
         End Do
         !flops = 5 * elemType%numDof**2 * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(MatElem)
      End If
   End Subroutine ElasticityBilinearFormSet

#undef __FUNCT__
#define __FUNCT__ "ElasticityOperatorSet"   
   Subroutine ElasticityOperatorSet(G,mesh,V,cellIS,A,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                       :: G
      Type(DM),Intent(IN)                                :: mesh
      Type(SectionReal),Intent(IN)                       :: V
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_TENS4OS),Intent(IN)                     :: A
      Type(MEF90_ELEMENTTYPE), Dimension(:), Pointer     :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscReal,Dimension(:),Pointer                     :: Gloc
      PetscReal,Dimension(:),Pointer                     :: Vloc
      Type(MEF90_VECT)                                   :: Velem
      Type(MEF90_MATS)                                   :: GradSVelem
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
               Velem      = 0.0_Kr
               GradSVElem = 0.0_Kr
               Do iDoF1 = 1,elemType%numDof
                  Velem      = Velem      + Vloc(iDof1) * elem(cell)%BF(iDoF1,iGauss)
                  GradSVelem = GradSVelem + Vloc(iDof1) * elem(cell)%GradS_BF(iDoF1,iGauss)
               End Do
               Do iDoF1 = 1,elemType%numDof
                  Gloc(iDoF1) = Gloc(iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                 (( (A * elem(cell)%GradS_BF(iDoF1,iGauss)) .DotP. GradSVelem))
               End Do
            End Do
            Call SectionRealUpdateClosure(G,mesh,cellID(cell),Gloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do
      
         !flops = 7 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(Gloc)
         DeAllocate(Vloc)
      End If
   End Subroutine ElasticityOperatorSet

#undef __FUNCT__
#define __FUNCT__ "ElasticityOperatorAddTransientTermSet"
!!!
!!!  
!!!  ElasticityOperatorAddTransientTermSet:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine ElasticityOperatorAddTransientTermSet(G,mesh,x,cellIS,alpha,elem,elemType,ierr)
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
      Type(MEF90_VECT)                                   :: xelem
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
                                (elem(cell)%BF(iDoF1,iGauss) .DotP. xelem)
               End Do
            End Do
            Gloc = Gloc * alpha
            Call SectionRealUpdateClosure(G,mesh,cellID(cell),Gloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do
      
         !flops = (5 * elemType%numDof * size(elem(1)%Gauss_C) + 1) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(xloc)
         DeAllocate(Gloc)
      End If
   End Subroutine ElasticityOperatorAddTransientTermSet

#undef __FUNCT__
#define __FUNCT__ "ElasticityRHSSetVertex"
   Subroutine ElasticityRHSSetVertex(RHS,mesh,F,cellIS,elem,elemType,ierr)
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
      Type(MEF90_VECT)                                   :: Felem
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
                  Felem = Felem + (Floc(iDof1) * elem(cell)%BF(iDoF1,iGauss))
               End Do
               Do iDoF1 = 1,elemType%numDof
                  RHSloc(iDoF1) = RHSloc(iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                 (elem(cell)%BF(iDoF1,iGauss) .DotP. Felem)
               End Do
            End Do
            Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do
      
         !flops = 5 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(RHSloc)
         DeAllocate(Floc)
      End If
   End Subroutine ElasticityRHSSetVertex

#undef __FUNCT__
#define __FUNCT__ "ElasticityRHSSetCell"
   Subroutine ElasticityRHSSetCell(RHS,mesh,F,cellIS,elem,elemType,ierr)
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
      Type(MEF90_VECT)                                   :: fVect
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
            fVect = Floc(1)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               Do iDoF1 = 1,elemType%numDof
                  RHSloc(iDoF1) = RHSloc(iDoF1) + (elem(cell)%Gauss_C(iGauss) * &
                                 (elem(cell)%BF(iDoF1,iGauss) .DotP. fVect))
               End Do
            End Do
            Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
            Call SectionRealRestore(F,cellID(cell),Floc,ierr);CHKERRQ(ierr)
         End Do
      
         !flops = 3 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(RHSloc)
      End If
   End Subroutine ElasticityRHSSetCell

#undef __FUNCT__
#define __FUNCT__ "ElasticityRHSSetCst"
   Subroutine ElasticityRHSSetCst(RHS,mesh,F,cellIS,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                       :: RHS
      Type(DM),Intent(IN)                                :: mesh
      Type(MEF90_VECT),Intent(IN)                        :: F
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
                                 (elem(cell)%BF(iDoF1,iGauss) .DotP. F)
               End Do
            End Do
            Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do
      
         !flops = 3 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(RHSloc)
      End If
   End Subroutine ElasticityRHSSetCst

#undef __FUNCT__
#define __FUNCT__ "ElasticityEnergySet"
   Subroutine ElasticityEnergySet(energy,x,mesh,A,cellIS,elem,elemType,ierr)
      PetscReal,Intent(OUT)                              :: energy
      Type(SectionReal),Intent(IN)                       :: x
      Type(DM),Intent(IN)                                :: mesh
      Type(MEF90_TENS4OS),Intent(IN)                     :: A
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENTTYPE), Dimension(:), Pointer     :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(MEF90_MATS)                                   :: stress,strain      
      Type(MEF90_VECT)                                   :: xelem
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
                  strain = strain + elem(cell)%GradS_BF(iDoF1,iGauss) * xloc(iDoF1)
                  xelem = xelem + elem(cell)%BF(iDoF1,iGauss) * xloc(iDoF1)
               End Do
               stress = A * strain
               energy = energy + elem(cell)%Gauss_C(iGauss) * (stress .dotP. strain) *.5_Kr
            End Do
         End Do
      
         !flops = (2 * elemType%numDof + 6) * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(xloc)
      End If
   End Subroutine ElasticityEnergySet

#undef __FUNCT__
#define __FUNCT__ "ElasticityWorkSetVertex"
   Subroutine ElasticityWorkSetVertex(work,x,mesh,F,cellIS,elem,elemType,ierr)
      PetscReal,Intent(OUT)                              :: work
      Type(SectionReal),Intent(IN)                       :: x
      Type(DM),Intent(IN)                                :: mesh
      Type(SectionReal),Intent(IN)                       :: F
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENTTYPE), Dimension(:), Pointer     :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr

      PetscReal,Dimension(:),Pointer                     :: xloc,floc
      Type(MEF90_VECT)                                   :: xelem,felem
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
      
         !flops = (4 * elemType%numDof + 3 )* size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(xloc)
         DeAllocate(floc)
      End If
   End Subroutine ElasticityWorkSetVertex

#undef __FUNCT__
#define __FUNCT__ "ElasticityWorkSetCell"
   Subroutine ElasticityWorkSetCell(work,x,mesh,F,cellIS,elem,elemType,ierr)
      PetscReal,Intent(OUT)                              :: work
      Type(SectionReal),Intent(IN)                       :: x
      Type(DM),Intent(IN)                                :: mesh
      Type(SectionReal),Intent(IN)                       :: F
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENTTYPE), Dimension(:), Pointer     :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr

      PetscReal,Dimension(:),Pointer                     :: xloc,floc
      Type(MEF90_VECT)                                   :: fVect
      Type(MEF90_VECT)                                   :: xelem
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
            fVect = floc(1)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               xelem = 0.0_Kr
               Do iDoF1 = 1,elemType%numDof
                  xelem = xelem + xloc(iDof1) * elem(cell)%BF(iDof1,iGauss)
               End Do
               work = work + elem(cell)%Gauss_C(iGauss) * (xelem .DotP. fVect)
            End Do
         End Do
      
         !flops = (4 * elemType%numDof + 3 )* size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(xloc)
         DeAllocate(floc)
      End If
   End Subroutine ElasticityWorkSetCell

#undef __FUNCT__
#define __FUNCT__ "ElasticityWorkSetCst"
   Subroutine ElasticityWorkSetCst(work,x,mesh,F,cellIS,elem,elemType,ierr)
      PetscReal,Intent(OUT)                              :: work
      Type(SectionReal),Intent(IN)                       :: x
      Type(DM),Intent(IN)                                :: mesh
      Type(MEF90_VECT),Intent(IN)                        :: F
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENTTYPE), Dimension(:), Pointer     :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr

      PetscReal,Dimension(:),Pointer                     :: xloc
      Type(MEF90_VECT)                                   :: xelem
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
               work = work + elem(cell)%Gauss_C(iGauss) * (F .DotP. xelem)
            End Do
         End Do
      
         !flops = (2 * elemType%numDof + 3 )* size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(xloc)
      End If
   End Subroutine ElasticityWorkSetCst

End Module MEF90_APPEND(m_MEF_ElasticityImplementation_,MEF90_DIM)D