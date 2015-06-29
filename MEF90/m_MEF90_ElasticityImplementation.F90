#include "mef90.inc"
Module MEF90_APPEND(m_MEF90_ElasticityImplementation_,MEF90_DIM)D
#include "finclude/petscdef.h"
   Use m_MEF90_LinAlg
   Use m_MEF90_Parameters
   Use m_MEF90_Elements
   Use m_MEF90_Utils
   Use m_MEF90_Materials
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
   Public :: ElasticityPressureForceRHSSetCst
   Public :: ElasticityPressureForceRHSSetCell
   Public :: ElasticityPressureForceRHSSetVertex
   Public :: ElasticityWorkSetCst
   Public :: ElasticityWorkSetCell
   Public :: ElasticityWorkSetVertex
   !Public :: ElasticityPressureWorkSetCst
   Public :: ElasticityPressureWorkSetCell
   !Public :: ElasticityPressureWorkSetVertex
   Public :: ElasticityStressSet

   

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
!!!
!!!  
!!!  ElasticityBilinearFormSet:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine ElasticityBilinearFormSet(K,mesh,cellIS,A,elem,elemType,ierr)
      Type(Mat),Intent(IN)                               :: K 
      Type(DM),Intent(IN)                                :: mesh
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_HOOKESLAW),Intent(IN)                   :: A
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
         DeAllocate(MatElem)
         flops = 2 * elemType%numDof**2 * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      End If
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(defaultSection,ierr);CHKERRQ(ierr)
   End Subroutine ElasticityBilinearFormSet

#undef __FUNCT__
#define __FUNCT__ "ElasticityOperatorSet"   
!!!
!!!  
!!!  ElasticityOperatorSet:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine ElasticityOperatorSet(G,mesh,U,cellIS,A,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                       :: G
      Type(DM),Intent(IN)                                :: mesh
      Type(SectionReal),Intent(IN)                       :: U
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_HOOKESLAW),Intent(IN)                   :: A
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
         DeAllocate(Gloc)
         DeAllocate(Uloc)
         flops = 2 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      End If
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine ElasticityOperatorSet

#undef __FUNCT__
#define __FUNCT__ "ElasticityOperatorAddTransientTermSet"
!!!
!!!  
!!!  ElasticityOperatorAddTransientTermSet: add the term corresponding to a \dot{x} to the operator
!!!  
!!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
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
      
         DeAllocate(xloc)
         DeAllocate(Gloc)
         flops = (2 * elemType%numDof * size(elem(1)%Gauss_C) + 1) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      End If
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine ElasticityOperatorAddTransientTermSet

#undef __FUNCT__
#define __FUNCT__ "ElasticityInelasticStrainRHSSetCst"
!!!
!!!  
!!!  ElasticityInelasticStrainRHSSetCst: Contribution of a constant inelastic strain e0
!!!                                      to a RHS in a cell set
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
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
           
      Call PetscPrintf(PETSC_COMM_WORLD,"Function not tested yet: "//__FUNCT__//"\n",ierr)
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Function not tested yet: "//__FUNCT__,ierr)
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(RHSloc(elemType%numDof))
         Do cell = 1,size(cellID)      
            RHSloc = 0.0_Kr
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               Do iDoF1 = 1,elemType%numDof
                  RHSloc(iDoF1) = RHSloc(iDoF1) - elem(cell)%Gauss_C(iGauss) * &
                                 (elem(cell)%GradS_BF(iDoF1,iGauss) .DotP. e0)
               End Do
            End Do
            Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do
      
         DeAllocate(RHSloc)
         flops = 2 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      End If
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine ElasticityInelasticStrainRHSSetCst

#undef __FUNCT__
#define __FUNCT__ "ElasticityInelasticStrainRHSSetCell"
!!!
!!!  
!!!  ElasticityInelasticStrainRHSSetCell: Contribution of a cell-based inelastic strain e0
!!!                                       to a RHS in a cell set
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine ElasticityInelasticStrainRHSSetCell(RHS,mesh,meshMatS,e0,HookesLaw,cellIS,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                       :: RHS
      Type(DM),Intent(IN)                                :: mesh,meshMatS
      Type(SectionReal),Intent(IN)                       :: e0
      Type(MEF90_HOOKESLAW),Intent(IN)                   :: HookesLaw
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
                  RHSloc(iDoF1) = RHSloc(iDoF1) - elem(cell)%Gauss_C(iGauss) * &
                                 (elem(cell)%GradS_BF(iDoF1,iGauss) .DotP. (HookesLaw * e0_elem))
               End Do
            End Do
            Call SectionRealRestore(e0,cellID(cell),e0_loc,ierr);CHKERRQ(ierr)
            Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do
      
         DeAllocate(RHSloc)
         flops = 2 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      End If
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine ElasticityInelasticStrainRHSSetCell

#undef __FUNCT__
#define __FUNCT__ "ElasticityInelasticStrainRHSSetCell2"
!!!
!!!  
!!!  ElasticityInelasticStrainRHSSetCell2: Contribution of a cell-based inelastic strain e0 K0
!!!                                        to a RHS in a cell set
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
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
           
      Call PetscPrintf(PETSC_COMM_WORLD,"Function not tested yet: "//__FUNCT__//"\n",ierr)
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Function not tested yet: "//__FUNCT__,ierr)
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(RHSloc(elemType%numDof))
         Do cell = 1,size(cellID)      
            RHSloc = 0.0_Kr
            Call SectionRealRestrict(e0,cellID(cell),e0loc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               Do iDoF1 = 1,elemType%numDof
                  RHSloc(iDoF1) = RHSloc(iDoF1) - elem(cell)%Gauss_C(iGauss) * &
                                 e0loc(1) * (elem(cell)%GradS_BF(iDoF1,iGauss) .DotP. K0)
               End Do
            End Do
            Call SectionRealRestore(e0,cellID(cell),e0loc,ierr);CHKERRQ(ierr)
            Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do
      
         DeAllocate(RHSloc)
         flops = 3 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      End If
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine ElasticityInelasticStrainRHSSetCell2

#undef __FUNCT__
#define __FUNCT__ "ElasticityInelasticStrainRHSSetVertex2"
!!!
!!!  
!!!  ElasticityInelasticStrainRHSSetVertex2: Contribution of a vertex-based inelastic strain e0 K0
!!!                                          to a RHS in a cell set
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
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
                  e0elem = e0elem + e0loc(iDof1) * elemScal(cell)%BF(iDoF1,iGauss)
               End Do
               Do iDoF1 = 1,elemType%numDof
                  RHSloc(iDoF1) = RHSloc(iDoF1) - elem(cell)%Gauss_C(iGauss) * &
                                 e0elem * (elem(cell)%GradS_BF(iDoF1,iGauss) .DotP. K0)
               End Do
            End Do
            Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do
      
         DeAllocate(RHSloc)
         DeAllocate(e0loc)
         flops = 5 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      End If
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine ElasticityInelasticStrainRHSSetVertex2

#undef __FUNCT__
#define __FUNCT__ "ElasticityForceRHSSetCst"
!!!
!!!  
!!!  ElasticityForceRHSSetCst:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
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
           
      Call PetscPrintf(PETSC_COMM_WORLD,"Function not tested yet: "//__FUNCT__//"\n",ierr)
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Function not tested yet: "//__FUNCT__,ierr)
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(RHSloc(elemType%numDof))
         Do cell = 1,size(cellID)      
            RHSloc = 0.0_Kr
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               Do iDoF1 = 1,elemType%numDof
                  RHSloc(iDoF1) = RHSloc(iDoF1) - elem(cell)%Gauss_C(iGauss) * &
                                 (elem(cell)%BF(iDoF1,iGauss) .DotP. F)
               End Do
            End Do
            Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do
      
         DeAllocate(RHSloc)
         flops = 2 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      End If
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine ElasticityFOrceRHSSetCst

#undef __FUNCT__
#define __FUNCT__ "ElasticityForceRHSSetCell"
!!!
!!!  
!!!  ElasticityForceRHSSetCell:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
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
      PetscInt                                           :: iDoF1,iGauss
      PetscLogDouble                                     :: flops
           
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(RHSloc(elemType%numDof))
         Do cell = 1,size(cellID)      
            RHSloc = 0.0_Kr
            Call SectionRealRestrict(F,cellID(cell),Floc,ierr);CHKERRQ(ierr)
            fVect = Floc
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               Do iDoF1 = 1,elemType%numDof
                  RHSloc(iDoF1) = RHSloc(iDoF1) - elem(cell)%Gauss_C(iGauss) * &
                                 (elem(cell)%BF(iDoF1,iGauss) .DotP. fVect)
               End Do
            End Do
            Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
            Call SectionRealRestore(F,cellID(cell),Floc,ierr);CHKERRQ(ierr)
         End Do
      
         flops = 2 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(RHSloc)
      End If
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine ElasticityForceRHSSetCell

#undef __FUNCT__
#define __FUNCT__ "ElasticityForceRHSSetVertex"
!!!
!!!  
!!!  ElasticityForceRHSSetVertex:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
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
           
      Call PetscPrintf(PETSC_COMM_WORLD,"Function not tested yet: "//__FUNCT__//"\n",ierr)
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Function not tested yet: "//__FUNCT__,ierr)
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
                  RHSloc(iDoF1) = RHSloc(iDoF1) - elem(cell)%Gauss_C(iGauss) * &
                                 (elem(cell)%BF(iDoF1,iGauss) .DotP. Felem)
               End Do
            End Do
            Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do
      
         flops = 2 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(RHSloc)
         DeAllocate(Floc)
      End If
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine ElasticityForceRHSSetVertex

#undef __FUNCT__
#define __FUNCT__ "ElasticityPressureForceRHSSetCst"
!!!
!!!   ElasticityPressureForceRHSSetCst contribution of a constant pressure force to a residual
!!!   (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine ElasticityPressureForceRHSSetCst(RHS,mesh,pressure,cellIS,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                       :: RHS
      Type(DM),Intent(IN)                                :: mesh
      PetscReal,Intent(IN)                               :: pressure
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscReal,Dimension(:),Pointer                     :: RHSloc
      PetscInt                                           :: iDoF1,iGauss,c,numCell
      PetscLogDouble                                     :: flops
           
      Call PetscPrintf(PETSC_COMM_WORLD,"Function not tested yet: "//__FUNCT__//"\n",ierr)
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Function not tested yet: "//__FUNCT__,ierr)
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Call DMMeshGetStratumSize(mesh,"height",0,numCell,ierr);CHKERRQ(ierr)
         Allocate(RHSloc(elemType%numDof))
         Do cell = 1,size(cellID)      
            RHSloc = 0.0_Kr
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               Do iDoF1 = 1,elemType%numDof
                  RHSloc(iDoF1) = RHSloc(iDoF1) - elem(cell)%Gauss_C(iGauss) * &
                                  pressure * (elem(cell)%BF(iDoF1,iGauss) .DotP. elem(cell)%outerNormal)
               End Do
            End Do
            Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do ! cell
      
         flops = 3 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(RHSloc)
      End If
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine ElasticityPressureForceRHSSetCst

#undef __FUNCT__
#define __FUNCT__ "ElasticityPressureForceRHSSetCell"
!!!
!!!   ElasticityPressureForceRHSSetCell contribution of a cell-based pressure force to a residual
!!!   (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine ElasticityPressureForceRHSSetCell(RHS,mesh,pressureSec,cellIS,elem,elemType,ierr)
      Type(SectionReal),Intent(IN)                       :: RHS
      Type(DM),Intent(IN)                                :: mesh
      Type(SectionReal),Intent(IN)                       :: pressureSec
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscReal,Dimension(:),Pointer                     :: pressureloc,RHSloc
      PetscInt                                           :: iDoF1,iGauss,c,numCell
      PetscLogDouble                                     :: flops
           
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(RHSloc(elemType%numDof))
         Do cell = 1,size(cellID)      
            RHSloc = 0.0_Kr
            Call SectionRealRestrict(pressureSec,cellID(cell),pressureloc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               Do iDoF1 = 1,elemType%numDof
                  RHSloc(iDoF1) = RHSloc(iDoF1) - elem(cell)%Gauss_C(iGauss) * &
                                  pressureloc(1) * (elem(cell)%BF(iDoF1,iGauss) .DotP. elem(cell)%outerNormal)
               End Do
            End Do
            Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
            Call SectionRealRestore(pressureSec,cellID(cell),pressureloc,ierr);CHKERRQ(ierr)
         End Do ! cell
      
         flops = 3 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(RHSloc)
      End If
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine ElasticityPressureForceRHSSetCell


#undef __FUNCT__
#define __FUNCT__ "ElasticityPressureForceRHSSetVertex"
!!!
!!!   ElasticityPressureForceRHSSetVertex contribution of a FE-based pressure force to a residual
!!!   (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine ElasticityPressureForceRHSSetVertex(RHS,meshDisplacement,meshPressure,pressureSec,cellIS,elemDisplacement,elemDisplacementType,elemPressure,elemPressureType,ierr)
      Type(SectionReal),Intent(IN)                       :: RHS
      Type(DM),Intent(IN)                                :: meshDisplacement,meshPressure
      Type(SectionReal),Intent(IN)                       :: pressureSec
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elemDisplacement
      Type(MEF90Element_Type),Intent(IN)                 :: elemDisplacementType
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elemPressure
      Type(MEF90Element_Type),Intent(IN)                 :: elemPressureType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscReal,Dimension(:),Pointer                     :: pressureloc,RHSloc
      PetscReal                                          :: pressureelem
      PetscInt                                           :: iDoF1,iGauss,c,numCell
      PetscLogDouble                                     :: flops
           
      Call PetscPrintf(PETSC_COMM_WORLD,"Function not tested yet: "//__FUNCT__//"\n",ierr)
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Function not tested yet: "//__FUNCT__,ierr)
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(pressureloc(elemPressureType%numDof))
         Allocate(RHSloc(elemDisplacementType%numDof))
         Do cell = 1,size(cellID)      
            RHSloc = 0.0_Kr
            Call SectionRealRestrictClosure(pressureSec,meshPressure,cellID(cell),elemPressureType%numDof,pressureloc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elemDisplacement(cell)%Gauss_C)
               pressureelem     = 0.0_Kr
               Do iDoF1 = 1,elemPressureType%numDof
                  pressureElem = pressureElem + pressureloc(iDof1) * elemPressure(cell)%BF(iDoF1,iGauss)
               End Do
               Do iDoF1 = 1,elemDisplacementType%numDof
                  RHSloc(iDoF1) = RHSloc(iDoF1) - elemDisplacement(cell)%Gauss_C(iGauss) * &
                                  pressureElem * (elemDisplacement(cell)%BF(iDoF1,iGauss) .DotP. elemDisplacement(cell)%outerNormal)
               End Do
            End Do
            Call SectionRealUpdateClosure(RHS,meshDisplacement,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do ! cell
      
         flops = (2 * elemPressureType%numDof + 3 * elemDisplacementType%numDof) * size(elemDisplacement(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(RHSloc)
         DeAllocate(pressureloc)
      End If
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine ElasticityPressureForceRHSSetVertex

#undef __FUNCT__
#define __FUNCT__ "ElasticityWorkSetCst"
!!!
!!!  
!!!  ElasticityWorkSetCst:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
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
     
      Call PetscPrintf(PETSC_COMM_WORLD,"Function not tested yet: "//__FUNCT__//"\n",ierr)
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Function not tested yet: "//__FUNCT__,ierr)
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
      
         flops = 2 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(xloc)
      End If
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine ElasticityWorkSetCst

#undef __FUNCT__
#define __FUNCT__ "ElasticityWorkSetCell"
!!!
!!!  
!!!  ElasticityWorkSetCell:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
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
            Call SectionRealRestore(f,cellID(cell),floc,ierr);CHKERRQ(ierr)
         End Do
      
         flops = 4 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(xloc)
      End If
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine ElasticityWorkSetCell

#undef __FUNCT__
#define __FUNCT__ "ElasticityWorkSetVertex"
!!!
!!!  
!!!  ElasticityWorkSetVertex:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
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
               work = work + elem(cell)%Gauss_C(iGauss) * (felem .dotP. xelem)
            End Do
         End Do
      
         flops = 2 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(xloc)
         DeAllocate(floc)
      End If
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine ElasticityWorkSetVertex

#undef __FUNCT__
#define __FUNCT__ "ElasticitypressureWorkSetCell"
!!!
!!!  
!!!  ElasticitypressureWorkSetCell:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine ElasticitypressureWorkSetCell(work,x,mesh,P,cellIS,elem,elemType,ierr)
      PetscReal,Intent(OUT)                              :: work
      Type(SectionReal),Intent(IN)                       :: x
      Type(DM),Intent(IN)                                :: mesh
      Type(SectionReal),Intent(IN)                       :: P
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr

      PetscReal,Dimension(:),Pointer                     :: xloc,ploc
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
            Call SectionRealRestrict(P,cellID(cell),ploc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               xelem = 0.0_Kr
               Do iDoF1 = 1,elemType%numDof
                  xelem = xelem + xloc(iDof1) * elem(cell)%BF(iDof1,iGauss)
               End Do
               work = work + ploc(1) * elem(cell)%Gauss_C(iGauss) * (xelem .DotP. elem(cell)%OuterNormal)
            End Do
            Call SectionRealRestore(P,cellID(cell),ploc,ierr);CHKERRQ(ierr)
         End Do
      
         flops = 3 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(xloc)
      End If
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine ElasticitypressureWorkSetCell

#undef __FUNCT__
#define __FUNCT__ "ElasticityEnergySet"
   !!!
   !!!  
   !!!  ElasticityEnergySet:  Contribution of a cell set to elastic energy. 
   !!!                        It is assumed that the temperature is interpolated on the FE space while the plastic strain 
   !!!                        is cell-based
   !!!  
   !!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
   !!!
   
   Subroutine ElasticityEnergySet(energy,x,plasticStrain,temperature,mesh,meshScal,cellIS,HookesLaw,ThermalExpansion,elemDisplacement,elemDisplacementType,elemTemperature,elemTemperatureType,ierr)
      PetscReal,Intent(OUT)                              :: energy
      Type(SectionReal),Intent(IN)                       :: x,plasticStrain,temperature
      Type(DM),Intent(IN)                                :: mesh,meshScal
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_HOOKESLAW),Intent(IN)                   :: HookesLaw
      Type(MEF90_MATS),Intent(IN)                        :: ThermalExpansion
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elemTemperature
      Type(MEF90Element_Type),Intent(IN)                 :: elemDisplacementType,elemTemperatureType
      PetscErrorCode,Intent(OUT)                         :: ierr

      PetscReal,Dimension(:),Pointer                     :: xloc,temperatureLoc,plasticStrainLoc
      PetscReal                                          :: temperatureElem
      Type(MEF90_MATS)                                   :: StrainElem,StressElem,plasticStrainElem
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscInt                                           :: iDoF1,iGauss
      PetscLogDouble                                     :: flops
     
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(xloc(elemDisplacementType%numDof))
         Allocate(temperatureloc(elemTemperatureType%numDof))
         Do cell = 1,size(cellID)   
            Call SectionRealRestrictClosure(x,mesh,cellID(cell),elemDisplacementType%numDof,xloc,ierr);CHKERRQ(ierr)
            If (temperature%v /= 0) Then
               Call SectionRealRestrictClosure(temperature,meshScal,cellID(cell),elemTemperatureType%numDof,temperatureLoc,ierr);CHKERRQ(ierr)
            End If
            If (plasticStrain%v /= 0) Then
               Call SectionRealRestrict(plasticStrain,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
            End If
            Do iGauss = 1,size(elemDisplacement(cell)%Gauss_C)
               temperatureElem   = 0.0_Kr
               plasticStrainElem = 0.0_Kr
               strainElem        = 0.0_Kr
               Do iDoF1 = 1,elemDisplacementType%numDof
                  strainElem = strainElem + xloc(iDof1) * elemDisplacement(cell)%GradS_BF(iDof1,iGauss)
               End Do
               If (temperature%v /= 0) Then
                  Do iDoF1 = 1,elemtemperatureType%numDof
                     temperatureElem = temperatureElem + temperatureLoc(iDof1) * elemTemperature(cell)%BF(iDof1,iGauss)
                  End Do
                  strainElem = strainElem - (temperatureElem * thermalExpansion)
               End If
               If (plasticStrain%v /= 0) Then
                  plasticStrainElem = plasticStrainLoc
                  strainElem = strainElem - plasticStrainElem
               End If
               stressElem = HookesLaw * strainElem
               energy = energy + (strainElem .dotP. stressElem) * elemDisplacement(cell)%Gauss_C(iGauss) * 0.5_Kr
            End Do ! Gauss
            If (plasticStrain%v /= 0) Then
               Call SectionRealRestore(plasticStrain,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
            End If
         End Do ! cell
         flops = 3 * size(elemDisplacement(1)%Gauss_C) * size(cellID) 
         If (temperature%v /= 0) Then
            flops = flops + 2 * elemtemperatureType%numDof * size(elemtemperature(1)%Gauss_C) * size(cellID)
         End If
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(xloc)
         DeAllocate(temperatureloc)
      End If   
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine ElasticityEnergySet
   
#undef __FUNCT__
#define __FUNCT__ "ElasticityStressSet"
!!!
!!!  
!!!  ElasticityStressSet: Compute the cell-averaged stress field associated with given displacement, temperature and plastic strain fields
!!!             
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine ElasticityStressSet(stress,x,plasticStrain,temperature,mesh,meshScal,cellIS,HookesLaw,ThermalExpansion,elemDisplacement,elemDisplacementType,elemTemperature,elemTemperatureType,ierr)
      Type(SectionReal),Intent(IN)                       :: stress,x,plasticStrain,temperature
      Type(DM),Intent(IN)                                :: mesh,meshScal
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_HOOKESLAW),Intent(IN)                   :: HookesLaw
      Type(MEF90_MATS),Intent(IN)                        :: ThermalExpansion
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elemTemperature
      Type(MEF90Element_Type),Intent(IN)                 :: elemDisplacementType,elemTemperatureType
      PetscErrorCode,Intent(OUT)                         :: ierr

      PetscReal,Dimension(:),Pointer                     :: xloc,temperatureLoc,plasticStrainLoc,stressPtr
      PetscReal                                          :: temperatureElem
      Type(MEF90_MATS)                                   :: StrainElem,plasticStrainElem,stressLoc
      PetscReal                                          :: cellSize
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscInt                                           :: iDoF1,iGauss
      PetscLogDouble                                     :: flops
      
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(xloc(elemDisplacementType%numDof))
         Allocate(temperatureloc(elemTemperatureType%numDof))
         Allocate(stressPtr(SIZEOFMEF90_MATS))
         Do cell = 1,size(cellID)  
            cellSize = 0.0_Kr 
            stressLoc = 0.0_Kr
            Call SectionRealRestrictClosure(x,mesh,cellID(cell),elemDisplacementType%numDof,xloc,ierr);CHKERRQ(ierr)
            If (temperature%v /= 0) Then
               Call SectionRealRestrictClosure(temperature,meshScal,cellID(cell),elemTemperatureType%numDof,temperatureLoc,ierr);CHKERRQ(ierr)
            End If
            If (plasticStrain%v /= 0) Then
               Call SectionRealRestrict(plasticStrain,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
            End If
            Do iGauss = 1,size(elemDisplacement(cell)%Gauss_C)
               temperatureElem   = 0.0_Kr
               plasticStrainElem = 0.0_Kr
               strainElem        = 0.0_Kr
               Do iDoF1 = 1,elemDisplacementType%numDof
                  strainElem = strainElem + xloc(iDof1) * elemDisplacement(cell)%GradS_BF(iDof1,iGauss)
               End Do
               If (temperature%v /= 0) Then
                  Do iDoF1 = 1,elemtemperatureType%numDof
                     temperatureElem = temperatureElem + temperatureLoc(iDof1) * elemTemperature(cell)%BF(iDof1,iGauss)
                  End Do
                  strainElem = strainElem - (temperatureElem * thermalExpansion)
               End If
               If (plasticStrain%v /= 0) Then
                  plasticStrainElem = plasticStrainLoc
                  strainElem = strainElem - plasticStrainElem
               End If
               stressLoc = stressLoc + elemDisplacement(cell)%Gauss_C(iGauss) * (HookesLaw * strainElem)
               cellSize = cellSize + elemDisplacement(cell)%Gauss_C(iGauss)
            End Do ! Gauss
            If (plasticStrain%v /= 0) Then
               Call SectionRealRestore(plasticStrain,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
            End If
            stressPtr = stressLoc / cellSize
            Call SectionRealUpdate(stress,cellID(cell),stressPtr,INSERT_VALUES,ierr);CHKERRQ(ierr)
         End Do ! cell
         flops = size(elemDisplacement(1)%Gauss_C) * size(cellID) 
         If (temperature%v /= 0) Then
            flops = flops + 2 * elemtemperatureType%numDof * size(elemtemperature(1)%Gauss_C) * size(cellID)
         End If
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(xloc)
         DeAllocate(temperatureloc)
         DeAllocate(stressPtr)
      End If   
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine ElasticityStressSet   




!!!      !!! InelasticStrainCell
!!!      
!!!      #undef __FUNCT__
!!!      #define __FUNCT__ "InelasticStrainCell"
!!!      !!!
!!!      !!!  
!!!      !!!  InelasticStrainCell: Compute the cell-averaged stress field associated with given displacement, temperature 
!!!      !!!             
!!!      !!!  
!!!      !!!  (c) 2014 Erwan TANNE erwan.tanne@gmail.com
!!!      !!!
!!!      
!!!         Subroutine InelasticStrainCell(stress,x,plasticStrain,temperature,mesh,meshScal,cell,HookesLaw,ThermalExpansion,elemDisplacement,elemDisplacementType,elemTemperature,elemTemperatureType,ierr)
!!!            Type(SectionReal),Intent(IN)                       :: stress,x,plasticStrain,temperature
!!!            Type(DM),Intent(IN)                                :: mesh,meshScal
!!!            !Type(IS),Intent(IN)                                :: cellIS
!!!            Type(MEF90_HOOKESLAW),Intent(IN)                   :: HookesLaw
!!!            Type(MEF90_MATS),Intent(IN)                        :: ThermalExpansion
!!!            Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elemDisplacement
!!!            Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elemTemperature
!!!            Type(MEF90Element_Type),Intent(IN)                 :: elemDisplacementType,elemTemperatureType
!!!            PetscErrorCode,Intent(OUT)                         :: ierr
!!!      
!!!            PetscReal,Dimension(:),Pointer                     :: xloc,temperatureLoc,plasticStrainLoc,stressPtr
!!!            PetscReal                                          :: temperatureElem
!!!            Type(MEF90_MATS)                                   :: StrainElem,plasticStrainElem,stressLoc
!!!            PetscReal                                          :: cellSize
!!!            PetscInt,Dimension(:),Pointer                      :: cellID
!!!            !!PetscInt                                           :: cell
!!!            PetscInt                                           :: iDoF1,iGauss
!!!            PetscLogDouble                                     :: flops
!!!            
!!!            !! Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
!!!            !!! If (Size(cellID) > 0) Then
!!!               Allocate(xloc(elemDisplacementType%numDof))
!!!               Allocate(temperatureloc(elemTemperatureType%numDof))
!!!               Allocate(stressPtr(SIZEOFMEF90_MATS))
!!!               !!!   Do cell = 1,size(cellID)  
!!!                  cellSize = 0.0_Kr 
!!!                  strainLoc = 0.0_Kr
!!!                  Call SectionRealRestrictClosure(x,mesh,cellID(cell),elemDisplacementType%numDof,xloc,ierr);CHKERRQ(ierr)
!!!                  If (temperature%v /= 0) Then
!!!                     Call SectionRealRestrictClosure(temperature,meshScal,cellID(cell),elemTemperatureType%numDof,temperatureLoc,ierr);CHKERRQ(ierr)
!!!                  End If
!!!                  !!! If (plasticStrain%v /= 0) Then
!!!                  !!!    Call SectionRealRestrict(plasticStrain,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
!!!                  !!! End If
!!!      
!!!                  Do iGauss = 1,size(elemDisplacement(cell)%Gauss_C)
!!!                     temperatureElem   = 0.0_Kr
!!!                     !!!   plasticStrainElem = 0.0_Kr
!!!                     strainElem        = 0.0_Kr
!!!                     Do iDoF1 = 1,elemDisplacementType%numDof
!!!                        strainElem = strainElem + xloc(iDof1) * elemDisplacement(cell)%GradS_BF(iDof1,iGauss)
!!!                     End Do
!!!                     If (temperature%v /= 0) Then
!!!                        Do iDoF1 = 1,elemtemperatureType%numDof
!!!                           temperatureElem = temperatureElem + temperatureLoc(iDof1) * elemTemperature(cell)%BF(iDof1,iGauss)
!!!                        End Do
!!!                        strainElem = strainElem - (temperatureElem * thermalExpansion)
!!!                     End If
!!!                     !!!   If (plasticStrain%v /= 0) Then
!!!                     !!!      plasticStrainElem = plasticStrainLoc
!!!                     !!!      strainElem = strainElem - plasticStrainElem
!!!                     End If
!!!                     strainLoc = strainLoc + elemDisplacement(cell)%Gauss_C(iGauss) * strainElem
!!!                     cellSize = cellSize + elemDisplacement(cell)%Gauss_C(iGauss)
!!!                  End Do ! Gauss
!!!      
!!!      
!!!      
!!!      
!!!                  !!!  If (plasticStrain%v /= 0) Then
!!!                  !!!     Call SectionRealRestore(plasticStrain,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
!!!                  !!!  End If
!!!      
!!!                  stressPtr = stressLoc / cellSize
!!!                  Call SectionRealUpdate(stress,cellID(cell),stressPtr,INSERT_VALUES,ierr);CHKERRQ(ierr)
!!!      
!!!      
!!!               !!!   End Do ! cell
!!!               flops = size(elemDisplacement(1)%Gauss_C) * size(cellID) 
!!!               If (temperature%v /= 0) Then
!!!                  flops = flops + 2 * elemtemperatureType%numDof * size(elemtemperature(1)%Gauss_C) * size(cellID)
!!!               !!!  End If
!!!               Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
!!!               DeAllocate(xloc)
!!!               DeAllocate(temperatureloc)
!!!               DeAllocate(stressPtr)
!!!            End If   
!!!            Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
!!!         End Subroutine InelasticStrainCell




End Module MEF90_APPEND(m_MEF90_ElasticityImplementation_,MEF90_DIM)D