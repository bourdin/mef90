#include "mef90.inc"
Module MEF90_APPEND(m_MEF_ElasticityImplementation_,MEF90_DIM)D
#include "finclude/petscdef.h"
   Use m_MEF_LinAlg
   Use m_MEF_Parameters
   Use m_MEF_Elements
   Use m_MEF_Utils
   Use petsc
   IMPLICIT NONE

   Private   
   Public :: ElasticityBilinearFormSet
   Public :: ElasticityEnergySet
   Public :: ElasticityOperatorAddTransientTermSet
   Public :: ElasticityOperatorSet
   Public :: ElasticityInelasticStrainRHSSetCst
   Public :: ElasticityInelasticStrainRHSSetCell
   Public :: ElasticityInelasticStrainRHSSetCell2
   Public :: ElasticityInelasticStrainRHSSetVertex2
   Public :: ElasticityForceRHSSetCst
   Public :: ElasticityForceRHSSetCell
   Public :: ElasticityForceRHSSetVertex
   Public :: ElasticityWorkSetCst
   Public :: ElasticityWorkSetCell
   Public :: ElasticityWorkSetVertex

   

!  Assembles all components required to solve a Elasticity equation in the form
!
!       { -div[A\e(v)-\e_0] = f   in \Omega
!  (1)  {                 v = v_0 on \partial \Omega_d
!       {       A(v-\e_0).n = g   on \partial \Omega_n
!
!  or equivalently to minimize
!
!  (2) E(v) := 1/2 \int A(\e(v)-\e_0):(\e(v)-\e_0) \, dx - \int_\Omega fv\, dx - \int_{\partial \Omega_n} gv\, dS

!  Nomenclature:
!     the "Bilinear Form" is K such that 
!        <Kv,w> := \int A\e(v) : e(w) \, dx
!     the "Operator" is v -> G(v)  where  
!        <G(v),w> := \int A\e(v):\e(w)\, dx
!     the "RHS" is F, where
!        <F,w> := \int_\Omega A\e_0:\e(w)\, dx + \int_\Omega fw\, dx + \int_{\partial \Omega_n} gw\, dS
!                  ^^^^^^^^^^^^^^^^^^^^^^^^^^^    ^^^^^^^^^^^^^^^^^^     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!                      "InelasticStrainRHS"            "ForceRHS"                    "ForceRHS"
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
! an extra layer of complexity comes from the fact that \e_0 can be either
!  - a symmetric matrix defined in each block (ElasticityInelasticStrainRHSSetCst)
!  - a symmetric matrix defined in each cell of a block (ElasticityInelasticStrainRHSSetCell)
!  - f K_0, where f is a scalar defined on each cell and K_0 is a symmetric matrix defined given in each block (ElasticityInelasticStrainRHSSetCell2)
!  - f K_0, where f is a scalar defined on the finite element space and K_0 is a symmetric matrix defined given in each block (ElasticityInelasticStrainRHSSetVertex2)


Contains
#undef __FUNCT__
#define __FUNCT__ "ElasticityBilinearFormSet"
   Subroutine ElasticityBilinearFormSet(K,mesh,cellIS,A,elem,elemType,ierr)
      Type(Mat),Intent(IN)                               :: K 
      Type(DM),Intent(IN)                                :: mesh
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_TENS4OS),Intent(IN)                     :: A
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elem
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
   Subroutine ElasticityOperatorSet(G,mesh,U,cellIS,A,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                       :: G
      Type(DM),Intent(IN)                                :: mesh
      Type(SectionReal),Intent(IN)                       :: U
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_TENS4OS),Intent(IN)                     :: A
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscReal,Dimension(:),Pointer                     :: Gloc
      PetscReal,Dimension(:),Pointer                     :: Uloc
      Type(MEF90_MATS)                                   :: GradSUelem
      PetscInt                                           :: iDoF1,iGauss
      PetscLogDouble                                     :: flops
           
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(Uloc(elemType%numDof))
         Allocate(Gloc(elemType%numDof))
         Do cell = 1,size(cellID)      
            Gloc = 0.0_Kr
            Call SectionRealRestrictClosure(U,mesh,cellID(cell),elemType%numDof,Uloc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               GradSUelem = 0.0_Kr
               Do iDoF1 = 1,elemType%numDof
                  GradSUelem = GradSUelem + Uloc(iDof1) * elem(cell)%GradS_BF(iDoF1,iGauss)
               End Do
               Do iDoF1 = 1,elemType%numDof
                  Gloc(iDoF1) = Gloc(iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                 (( (A * elem(cell)%GradS_BF(iDoF1,iGauss)) .DotP. GradSUelem))
               End Do
            End Do
            Call SectionRealUpdateClosure(G,mesh,cellID(cell),Gloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do      
         !flops = 7 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(Gloc)
         DeAllocate(Uloc)
      End If
   End Subroutine ElasticityOperatorSet

#undef __FUNCT__
#define __FUNCT__ "ElasticityOperatorAddTransientTermSet"
!!!
!!!  
!!!  ElasticityOperatorAddTransientTermSet:
!!!  
!!!  (c) 2012-13 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine ElasticityOperatorAddTransientTermSet(G,mesh,x,cellIS,alpha,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                       :: G
      Type(DM),Intent(IN)                                :: mesh
      Type(SectionReal),Intent(IN)                       :: x
      Type(IS),Intent(IN)                                :: cellIS
      PetscReal,Intent(IN)                               :: alpha
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elem
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
#define __FUNCT__ "ElasticityInelasticStrainRHSSetCst"
   Subroutine ElasticityInelasticStrainRHSSetCst(RHS,mesh,e0,cellIS,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                       :: RHS
      Type(DM),Intent(IN)                                :: mesh
      Type(MEF90_MATS),Intent(IN)                        :: e0
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elem
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
                                 (elem(cell)%GradS_BF(iDoF1,iGauss) .DotP. e0)
               End Do
            End Do
            Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do
      
         !flops = 3 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(RHSloc)
      End If
   End Subroutine ElasticityInelasticStrainRHSSetCst

#undef __FUNCT__
#define __FUNCT__ "ElasticityInelasticStrainRHSSetCell"
   Subroutine ElasticityInelasticStrainRHSSetCell(RHS,mesh,meshMatS,e0,cellIS,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                       :: RHS
      Type(DM),Intent(IN)                                :: mesh,meshMatS
      Type(SectionReal),Intent(IN)                       :: e0
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      Type(MEF90_MATS)                                   :: e0_elem
      PetscReal,Dimension(:),Pointer                     :: e0_loc
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscReal,Dimension(:),Pointer                     :: RHSloc
      PetscInt                                           :: iDoF1,iGauss,i
      PetscLogDouble                                     :: flops
           
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(RHSloc(elemType%numDof))
         Do cell = 1,size(cellID)      
            RHSloc = 0.0_Kr
            Call SectionRealRestrict(e0,cellID(cell),e0_loc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               e0_elem = e0_loc
               Do iDoF1 = 1,elemType%numDof
                  RHSloc(iDoF1) = RHSloc(iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                 (elem(cell)%GradS_BF(iDoF1,iGauss) .DotP. e0_elem)
               End Do
            End Do
            Call SectionRealRestore(e0,cellID(cell),e0_loc,ierr);CHKERRQ(ierr)
            Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do
      
         !flops = 3 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(RHSloc)
      End If
   End Subroutine ElasticityInelasticStrainRHSSetCell

#undef __FUNCT__
#define __FUNCT__ "ElasticityInelasticStrainRHSSetCell2"
   Subroutine ElasticityInelasticStrainRHSSetCell2(RHS,mesh,meshMatS,e0,K0,cellIS,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                       :: RHS
      Type(DM),Intent(IN)                                :: mesh,meshMatS
      Type(SectionReal),Intent(IN)                       :: e0
      Type(MEF90_MATS),Intent(IN)                        :: K0 ! The inelastic strain is e0.K0
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscReal,Dimension(:),Pointer                     :: RHSloc,e0loc
      PetscInt                                           :: iDoF1,iGauss,i
      PetscLogDouble                                     :: flops
           
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(RHSloc(elemType%numDof))
         Do cell = 1,size(cellID)      
            RHSloc = 0.0_Kr
            Call SectionRealRestrict(e0,cellID(cell),e0loc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               Do iDoF1 = 1,elemType%numDof
                  RHSloc(iDoF1) = RHSloc(iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                 e0loc(1) * (elem(cell)%GradS_BF(iDoF1,iGauss) .DotP. K0)
               End Do
            End Do
            Call SectionRealRestore(e0,cellID(cell),e0loc,ierr);CHKERRQ(ierr)
            Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do
      
         !flops = 3 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(RHSloc)
      End If
   End Subroutine ElasticityInelasticStrainRHSSetCell2

#undef __FUNCT__
#define __FUNCT__ "ElasticityInelasticStrainRHSSetVertex2"
   Subroutine ElasticityInelasticStrainRHSSetVertex2(RHS,mesh,meshScal,e0,K0,cellIS,elem,elemType,elemScal,elemScalType,ierr)
      Type(SectionReal),Intent(IN)                       :: RHS
      Type(DM),Intent(IN)                                :: mesh,meshScal
      Type(SectionReal),Intent(IN)                       :: e0
      Type(MEF90_MATS),Intent(IN)                        :: K0 ! The inelastic strain is e0.K0
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elem
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elemScal
      Type(MEF90Element_Type),Intent(IN)                 :: elemType,elemScalType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscReal,Dimension(:),Pointer                     :: RHSloc,e0loc
      PetscReal                                          :: e0elem
      PetscInt                                           :: iDoF1,iGauss
      PetscLogDouble                                     :: flops
           
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(e0loc(elemScalType%numDof))
         Allocate(RHSloc(elemType%numDof))
         Do cell = 1,size(cellID)      
            RHSloc = 0.0_Kr
            Call SectionRealRestrictClosure(e0,mesh,cellID(cell),elemScalType%numDof,e0loc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               e0elem = 0.0_Kr
               Do iDoF1 = 1,elemScalType%numDof
                  e0elem = e0elem + (e0loc(iDof1) * elemScal(cell)%BF(iDoF1,iGauss))
               End Do
               Do iDoF1 = 1,elemType%numDof
                  RHSloc(iDoF1) = RHSloc(iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                 e0elem * (elem(cell)%GradS_BF(iDoF1,iGauss) .DotP. K0)
               End Do
            End Do
            Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do
      
         !flops = 5 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         DeAllocate(RHSloc)
         DeAllocate(e0loc)
      End If
   End Subroutine ElasticityInelasticStrainRHSSetVertex2

#undef __FUNCT__
#define __FUNCT__ "ElasticityForceRHSSetCst"
   Subroutine ElasticityForceRHSSetCst(RHS,mesh,F,cellIS,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                       :: RHS
      Type(DM),Intent(IN)                                :: mesh
      Type(MEF90_VECT),Intent(IN)                        :: F
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elem
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
   End Subroutine ElasticityFOrceRHSSetCst

#undef __FUNCT__
#define __FUNCT__ "ElasticityForceRHSSetCell"
   Subroutine ElasticityForceRHSSetCell(RHS,mesh,F,cellIS,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                       :: RHS
      Type(DM),Intent(IN)                                :: mesh
      Type(SectionReal),Intent(IN)                       :: F
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elem
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
   End Subroutine ElasticityForceRHSSetCell

#undef __FUNCT__
#define __FUNCT__ "ElasticityForceRHSSetVertex"
   Subroutine ElasticityForceRHSSetVertex(RHS,mesh,F,cellIS,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                       :: RHS
      Type(DM),Intent(IN)                                :: mesh
      Type(SectionReal),Intent(IN)                       :: F
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elem
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
   End Subroutine ElasticityForceRHSSetVertex

#undef __FUNCT__
#define __FUNCT__ "ElasticityEnergySet"
   Subroutine ElasticityEnergySet(energy,x,mesh,A,cellIS,elem,elemType,ierr)
      PetscReal,Intent(OUT)                              :: energy
      Type(SectionReal),Intent(IN)                       :: x
      Type(DM),Intent(IN)                                :: mesh
      Type(MEF90_TENS4OS),Intent(IN)                     :: A
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elem
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

#undef __FUNCT__
#define __FUNCT__ "ElasticityWorkSetCell"
   Subroutine ElasticityWorkSetCell(work,x,mesh,F,cellIS,elem,elemType,ierr)
      PetscReal,Intent(OUT)                              :: work
      Type(SectionReal),Intent(IN)                       :: x
      Type(DM),Intent(IN)                                :: mesh
      Type(SectionReal),Intent(IN)                       :: F
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elem
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
#define __FUNCT__ "ElasticityWorkSetVertex"
   Subroutine ElasticityWorkSetVertex(work,x,mesh,F,cellIS,elem,elemType,ierr)
      PetscReal,Intent(OUT)                              :: work
      Type(SectionReal),Intent(IN)                       :: x
      Type(DM),Intent(IN)                                :: mesh
      Type(SectionReal),Intent(IN)                       :: F
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elem
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

End Module MEF90_APPEND(m_MEF_ElasticityImplementation_,MEF90_DIM)D