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
   Public :: MEF90_DiffusionEnergySet
   
   Interface MEF90_DiffusionBilinearFormSet
      Module Procedure DiffusionBilinearFormSet_2D, DiffusionBilinearFormSet_3D
   End Interface MEF90_DiffusionBilinearFormSet
   
   Interface MEF90_DiffusionRHSSet
      Module procedure DiffusionRHSSet_2D
   End Interface MEF90_DiffusionRHSSet

   Interface MEF90_DiffusionOperatorSet
      Module procedure DiffusionOperatorSet_2D
   End Interface MEF90_DiffusionOperatorSet

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
   Subroutine DiffusionBilinearFormSet_2D(K,mesh,U,cellIS,A,lambda,elem,elemType,BC,ierr)
      Type(Mat),Intent(IN)                         :: K
      Type(DM),Intent(IN)                          :: mesh
      Type(SectionReal),Intent(IN)                 :: U
      Type(IS),Intent(IN)                          :: cellIS
      Type(MatS2D),Intent(IN)                      :: A
      PetscReal,Intent(IN)                         :: lambda
      Type(Element2D_Scal), Dimension(:), Pointer  :: elem
      Type(Element_Type),Intent(IN)                :: elemType
      Type(Flag),Intent(IN),Optional               :: BC
      PetscErrorCode,Intent(OUT)                   :: ierr
      
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscReal,Dimension(:,:),Pointer             :: MatElem
      PetscInt,Dimension(:),Pointer                :: BCFlag
      PetscInt                                     :: iDoF1,iDoF2,iGauss
      PetscLogDouble                               :: flops
     
      flops = 0
      Allocate(MatElem(elemType%numDof,elemType%numDof))
      Allocate(BCFlag(elemType%numDof))

      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Do cell = 1,size(cellID)      
         MatElem = 0.0_Kr
         BCFlag = 0
         If (present(BC)) Then
            Call SectionIntRestrictClosure(BC%Sec,mesh,cellID(cell),elemType%numDof,BCFlag,ierr);CHKERRQ(ierr)
         End If
         Do iGauss = 1,size(elem(cellID(cell)+1)%Gauss_C)
            Do iDoF1 = 1,elemType%numDof
               If (BCFlag(iDoF1) == 0) Then
                  Do iDoF2 = 1,elemType%numDof
                     If (BCFlag(iDoF2) == 0) Then
                        MatElem(iDoF2,iDoF1) = MatElem(iDoF2,iDoF1) + elem(cellID(cell)+1)%Gauss_C(iGauss) * &
                                       (lambda * elem(cellID(cell)+1)%BF(iDoF1,iGauss) * elem(cellID(cell)+1)%BF(iDoF2,iGauss) + &
                                       ( (A * elem(cellID(cell)+1)%Grad_BF(iDoF1,iGauss)) .DotP. elem(cellID(cell)+1)%Grad_BF(iDoF2,iGauss)))
                        flops = flops + 5
                     End If
                  End Do
               End If
            End Do
         End Do
         Call DMmeshAssembleMatrix(K,mesh,U,cellID(cell),MatElem,ADD_VALUES,ierr);CHKERRQ(ierr)
         !!
         !! Another way would be to get the point number through the cone of the element
         !! then calling MatSetValuesTopology
         !! The caveat is that this would only work if dof are only vertices
         !! but this we already make this asumption in many other locations
      End Do
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      DeAllocate(BCFlag)
      DeAllocate(MatElem)
   End Subroutine DiffusionBilinearFormSet_2D

   Subroutine DiffusionOperatorSet_2D(G,mesh,V,cellIS,A,lambda,elem,elemType,BC,ierr)
      Type(SectionReal),Intent(IN)                 :: G
      Type(DM),Intent(IN)                          :: mesh
      Type(SectionReal),Intent(IN)                 :: V
      Type(IS),Intent(IN)                          :: cellIS
      Type(MatS2D),Intent(IN)                      :: A
      PetscReal,Intent(IN)                         :: lambda
      Type(Element2D_Scal), Dimension(:), Pointer  :: elem
      Type(Element_Type),Intent(IN)                :: elemType
      Type(Flag),Intent(IN),Optional               :: BC
      PetscErrorCode,Intent(OUT)                   :: ierr
      
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscReal,Dimension(:),Pointer               :: Gloc
      PetscReal,Dimension(:),Pointer               :: Vloc
      PetscReal                                    :: Velem
      Type(Vect2D)                                 :: GradVelem
      PetscInt,Dimension(:),Pointer                :: BCFlag
      PetscInt                                     :: iDoF1,iDoF2,iGauss
      PetscLogDouble                               :: flops
     
      flops = 0
      Allocate(Vloc(elemType%numDof))
      Allocate(Gloc(elemType%numDof))
      Allocate(BCFlag(elemType%numDof))
      
      !sum_g cg  [ A (sum_i v_i  \nabla phi_i(xi_g))\cdot \nabla phi_j(xi_g) + \lambda (sum_i v_i phi_i(xi_g)) phi_j(xi_g)]
      ! GradVElem = (sum_i v_i  \nabla phi_i(xi_g))
      ! VElem     = (sum_i v_i  phi_i(xi_g))
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Do cell = 1,size(cellID)      
         Gloc = 0.0_Kr
         BCFlag = 0
         If (present(BC)) Then
            Call SectionIntRestrictClosure(BC%Sec,mesh,cellID(cell),elemType%numDof,BCFlag,ierr);CHKERRQ(ierr)
         End If
         Call SectionRealRestrictClosure(V,mesh,cellID(cell),elemType%numDof,Vloc,ierr);CHKERRQ(ierr)
         Do iGauss = 1,size(elem(cellID(cell)+1)%Gauss_C)
            Velem     = 0.0_Kr
            GradVElem = 0.0_Kr
            Do iDoF1 = 1,elemType%numDof
               Velem     = Velem     + Vloc(iDof1) * elem(cellID(cell)+1)%BF(iDoF1,iGauss)
               GradVelem = GradVelem + Vloc(iDof1) * elem(cellID(cell)+1)%Grad_BF(iDoF1,iGauss)
            End Do
            flops = flops + 2 * elemType%numDof
            Do iDoF1 = 1,elemType%numDof
               If (BCFlag(iDoF1) == 0) Then
                  Gloc(iDoF1) = Gloc(iDoF1) + elem(cellID(cell)+1)%Gauss_C(iGauss) * &
                                 (lambda * elem(cellID(cell)+1)%BF(iDoF1,iGauss) * Velem + &
                                  ( (A * elem(cellID(cell)+1)%Grad_BF(iDoF1,iGauss)) .DotP. GradVelem))
                  flops = flops + 5
               End If
            End Do
         End Do
         Call SectionRealUpdateClosure(G,mesh,cellID(cell),Gloc,ADD_VALUES,ierr);CHKERRQ(iErr)
      End Do
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      DeAllocate(BCFlag)
      DeAllocate(Gloc)
      DeAllocate(Vloc)
   End Subroutine DiffusionOperatorSet_2D

   Subroutine DiffusionRHSSet_2D(RHS,mesh,F,cellIS,elem,elemType,BC,ierr)
      Type(SectionReal),Intent(IN)                 :: RHS
      Type(DM),Intent(IN)                          :: mesh
      Type(SectionReal),Intent(IN)                 :: F
      Type(IS),Intent(IN)                          :: cellIS
      Type(Element2D_Scal), Dimension(:), Pointer  :: elem
      Type(Element_Type),Intent(IN)                :: elemType
      Type(Flag),Intent(IN),Optional               :: BC
      PetscErrorCode,Intent(OUT)                   :: ierr
      
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscReal,Dimension(:),Pointer               :: RHSloc
      PetscReal,Dimension(:),Pointer               :: Vloc
      PetscReal                                    :: Felem
      PetscInt,Dimension(:),Pointer                :: BCFlag
      PetscInt                                     :: iDoF1,iDoF2,iGauss
      PetscLogDouble                               :: flops
     
      flops = 0
      Allocate(Vloc(elemType%numDof))
      Allocate(Gloc(elemType%numDof))
      Allocate(BCFlag(elemType%numDof))
      
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Do cell = 1,size(cellID)      
         Gloc = 0.0_Kr
         BCFlag = 0
         If (present(BC)) Then
            Call SectionIntRestrictClosure(BC%Sec,mesh,cellID(cell),elemType%numDof,BCFlag,ierr);CHKERRQ(ierr)
         End If
         Call SectionRealRestrictClosure(F,mesh,cellID(cell),elemType%numDof,Vloc,ierr);CHKERRQ(ierr)
         Do iGauss = 1,size(elem(cellID(cell)+1)%Gauss_C)
            Felem     = 0.0_Kr
            Do iDoF1 = 1,elemType%numDof
               Felem = Felem + Floc(iDof1) * elem(cellID(cell)+1)%BF(iDoF1,iGauss)
            End Do
            flops = flops + 2 * elemType%numDof
            Do iDoF1 = 1,elemType%numDof
               If (BCFlag(iDoF1) == 0) Then
                  RHSloc(iDoF1) = RHSloc(iDoF1) + elem(cellID(cell)+1)%Gauss_C(iGauss) * &
                                 elem(cellID(cell)+1)%BF(iDoF1,iGauss) * Felem
                  flops = flops + 3
               End If
            End Do
         End Do
         Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
      End Do
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      DeAllocate(BCFlag)
      DeAllocate(RHSloc)
      DeAllocate(Floc)
   End Subroutine DiffusionRHSSet_2D

   Subroutine DiffusionBilinearFormSet_3D(K,mesh,U,cellIS,A,lambda,elem,elemType,BC,ierr)
      Type(Mat),Intent(IN)                         :: K
      Type(DM),Intent(IN)                          :: mesh
      Type(SectionReal),Intent(IN)                 :: U
      Type(IS),Intent(IN)                          :: cellIS
      Type(MatS3D),Intent(IN)                      :: A
      PetscReal,Intent(IN)                         :: lambda
      Type(Element3D_Scal), Dimension(:), Pointer  :: elem
      Type(Element_Type),Intent(IN)                :: elemType
      Type(Flag),Intent(IN),Optional               :: BC
      PetscErrorCode,Intent(OUT)                   :: ierr
      
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscReal,Dimension(:,:),Pointer             :: MatElem
      PetscInt,Dimension(:),Pointer                :: BCFlag
      PetscInt                                     :: iDoF1,iDoF2,iGauss
      PetscLogDouble                               :: flops
     
      flops = 0
      Allocate(MatElem(elemType%numDof,elemType%numDof))
      Allocate(BCFlag(elemType%numDof))

      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Do cell = 1,size(cellID)      
         MatElem = 0.0_Kr
         BCFlag = 0
         If (present(BC)) Then
            Call SectionIntRestrictClosure(BC%Sec,mesh,cellID(cell),elemType%numDof,BCFlag,ierr);CHKERRQ(ierr)
         End If
         Do iGauss = 1,size(elem(cellID(cell)+1)%Gauss_C)
            Do iDoF1 = 1,elemType%numDof
               If (BCFlag(iDoF1) == 0) Then
                  Do iDoF2 = 1,elemType%numDof
                     If (BCFlag(iDoF2) == 0) Then
                        MatElem(iDoF2,iDoF1) = MatElem(iDoF2,iDoF1) + elem(cellID(cell)+1)%Gauss_C(iGauss) * &
                                       (lambda * elem(cellID(cell)+1)%BF(iDoF1,iGauss) * elem(cellID(cell)+1)%BF(iDoF2,iGauss) + &
                                       ( (A * elem(cellID(cell)+1)%Grad_BF(iDoF1,iGauss)) .DotP. elem(cellID(cell)+1)%Grad_BF(iDoF2,iGauss)))
                        flops = flops + 5
                     End If
                  End Do
               End If
            End Do
         End Do
         Call DMmeshAssembleMatrix(K,mesh,U,cellID(cell),MatElem,ADD_VALUES,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      DeAllocate(BCFlag)
      DeAllocate(MatElem)
   End Subroutine DiffusionBilinearFormSet_3D
End Module m_MEF_Diffusion
