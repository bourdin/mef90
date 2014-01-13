#include "../MEF90/mef90.inc"
Module MEF90_APPEND(m_MEF90_GradDamageImplementation_,MEF90_DIM)D
#include "finclude/petscdef.h"
   Use m_MEF90_LinAlg
   Use m_MEF90_Parameters
   Use m_MEF90_Elements
   Use m_MEF90_Utils
   IMPLICIT NONE
   Private

   Public :: MEF90GradDamageDispBilinearFormSet
   Public :: MEF90GradDamageDispOperatorSet
   Public :: MEF90GradDamageDispInelasticStrainRHSSetVertex2
   Public :: MEF90GradDamageDispInelasticStrainRHSSetCell
Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageDispBilinearFormSet"
!!!
!!!  
!!!  MEF90GradDamageDispBilinearFormSet:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90GradDamageDispBilinearFormSet(K,mesh,meshScal,cellIS,A,alpha,eta,elem,elemType,elemScal,elemScalType,ierr)
      Type(Mat),Intent(IN)                               :: K 
      Type(DM),Intent(IN)                                :: mesh,meshScal
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_TENS4OS),Intent(IN)                     :: A
      Type(SectionReal),Intent(IN)                       :: alpha
      PetscReal,Intent(IN)                               :: eta
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elem
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elemScal
      Type(MEF90Element_Type),Intent(IN)                 :: elemType,elemScalType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscReal,Dimension(:,:),Pointer                   :: MatElem
      PetscReal,Dimension(:),Pointer                     :: alphaLoc
      PetscReal                                          :: alphaElem,SalphaElem
      PetscInt                                           :: iDoF1,iDoF2,iGauss
      Type(SectionReal)                                  :: defaultSection
      PetscLogDouble                                     :: flops
     
      Call DMMeshGetSectionReal(mesh,'default',defaultSection,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(alphaLoc(elemScalType%numDof))
         Allocate(MatElem(elemType%numDof,elemType%numDof))
         Do cell = 1,size(cellID)   
            MatElem = 0.0_Kr
            Call SectionRealRestrictClosure(alpha,meshScal,cellID(cell),elemScalType%numDof,alphaLoc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               alphaElem = 0.0_Kr
               Do iDoF1 = 1, elemScalType%numDof
                  alphaElem = alphaElem + alphaLoc(iDoF1) * elemScal(cell)%BF(iDoF1,iGauss)
               End Do ! iDoF1
               SalphaElem = eta + (1.0_Kr - alphaElem)**2
               Do iDoF1 = 1,elemType%numDof
                  Do iDoF2 = 1,elemType%numDof
                     MatElem(iDoF2,iDoF1) = MatElem(iDoF2,iDoF1) + elem(cell)%Gauss_C(iGauss) * SalphaElem * &
                                    ((A * elem(cell)%GradS_BF(iDoF1,iGauss)) .DotP. elem(cell)%GradS_BF(iDoF2,iGauss))
                  End Do
               End Do
            End Do
            Call DMmeshAssembleMatrix(K,mesh,defaultSection,cellID(cell),MatElem,ADD_VALUES,ierr);CHKERRQ(ierr)
         End Do
         !flops = 5 * elemType%numDof**2 * size(elem(1)%Gauss_C) * size(cellID) 
         !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(alphaLoc)
         DeAllocate(MatElem)
      End If   
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageDispBilinearFormSet
   
#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageDispOperatorSet"   
!!!
!!!  
!!!  MEF90GradDamageDispOperatorSet:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90GradDamageDispOperatorSet(G,mesh,meshScal,U,alpha,eta,cellIS,A,elem,elemType,elemScal,elemScalType,ierr)
      Type(SectionReal),Intent(IN)                       :: G
      Type(DM),Intent(IN)                                :: mesh,meshScal
      Type(SectionReal),Intent(IN)                       :: U,alpha
      PetscReal,Intent(IN)                               :: eta
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_TENS4OS),Intent(IN)                     :: A
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elem
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elemScal
      Type(MEF90Element_Type),Intent(IN)                 :: elemType,elemScalType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscReal,Dimension(:),Pointer                     :: Gloc
      PetscReal,Dimension(:),Pointer                     :: Uloc
      PetscReal,Dimension(:),Pointer                     :: alphaloc
      Type(MEF90_MATS)                                   :: GradSUelem
      PetscReal                                          :: alphaElem,SalphaElem
      PetscInt                                           :: iDoF1,iGauss
      PetscLogDouble                                     :: flops
           
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(Uloc(elemType%numDof))
         Allocate(Gloc(elemType%numDof))
         Allocate(alphaLoc(elemScalType%numDof))
         Do cell = 1,size(cellID)      
            Gloc = 0.0_Kr
            Call SectionRealRestrictClosure(U,mesh,cellID(cell),elemType%numDof,Uloc,ierr);CHKERRQ(ierr)
            Call SectionRealRestrictClosure(alpha,meshScal,cellID(cell),elemScalType%numDof,alphaLoc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               alphaElem = 0.0_Kr
               Do iDoF1 = 1, elemScalType%numDof
                  alphaElem = alphaElem + alphaLoc(iDoF1) * elemScal(cell)%BF(iDoF1,iGauss)
               End Do ! iDoF1
               SalphaElem = eta + (1.0_Kr - alphaElem)**2
               GradSUelem = 0.0_Kr
               Do iDoF1 = 1,elemType%numDof
                  GradSUelem = GradSUelem + Uloc(iDof1) * elem(cell)%GradS_BF(iDoF1,iGauss)
               End Do
               Do iDoF1 = 1,elemType%numDof
                  Gloc(iDoF1) = Gloc(iDoF1) + elem(cell)%Gauss_C(iGauss) * SalphaElem * &
                                 (( (A * elem(cell)%GradS_BF(iDoF1,iGauss)) .DotP. GradSUelem)) 
               End Do
            End Do
            Call SectionRealUpdateClosure(G,mesh,cellID(cell),Gloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do      
         DeAllocate(Gloc)
         DeAllocate(Uloc)
         DeAllocate(alphaLoc)
         !flops = 7 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      End If
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageDispOperatorSet

#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageDispInelasticStrainRHSSetVertex2"
!!!
!!!  
!!!  MEF90GradDamageDispInelasticStrainRHSSetVertex2: Contribution of a vertex-based inelastic strain e0 K0
!!!                                          to a RHS in a cell set
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90GradDamageDispInelasticStrainRHSSetVertex2(RHS,mesh,meshScal,e0,K0,alpha,eta,cellIS,elem,elemType,elemScal,elemScalType,ierr)
      Type(SectionReal),Intent(IN)                       :: RHS
      Type(DM),Intent(IN)                                :: mesh,meshScal
      Type(SectionReal),Intent(IN)                       :: e0
      Type(MEF90_MATS),Intent(IN)                        :: K0 ! The inelastic strain is e0.K0
      Type(SectionReal),Intent(IN)                       :: alpha
      PetscReal,Intent(IN)                               :: eta
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elem
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elemScal
      Type(MEF90Element_Type),Intent(IN)                 :: elemType,elemScalType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscReal,Dimension(:),Pointer                     :: RHSloc,e0loc,alphaloc
      PetscReal                                          :: e0elem
      PetscReal                                          :: alphaElem,SalphaElem
      PetscInt                                           :: iDoF1,iGauss
      PetscLogDouble                                     :: flops
           
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(e0loc(elemScalType%numDof))
         Allocate(RHSloc(elemType%numDof))
         Allocate(alphaloc(elemScalType%numDof))
         Do cell = 1,size(cellID)      
            RHSloc = 0.0_Kr
            Call SectionRealRestrictClosure(e0,mesh,cellID(cell),elemScalType%numDof,e0loc,ierr);CHKERRQ(ierr)
            Call SectionRealRestrictClosure(alpha,meshScal,cellID(cell),elemScalType%numDof,alphaLoc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               alphaElem = 0.0_Kr
               Do iDoF1 = 1, elemScalType%numDof
                  alphaElem = alphaElem + alphaLoc(iDoF1) * elemScal(cell)%BF(iDoF1,iGauss)
               End Do ! iDoF1
               SalphaElem = eta + (1.0_Kr - alphaElem)**2
               e0elem = 0.0_Kr
               Do iDoF1 = 1,elemScalType%numDof
                  e0elem = e0elem + (e0loc(iDof1) * elemScal(cell)%BF(iDoF1,iGauss))
               End Do
               Do iDoF1 = 1,elemType%numDof
                  RHSloc(iDoF1) = RHSloc(iDoF1) - elem(cell)%Gauss_C(iGauss) * SalphaElem * &
                                 e0elem * (elem(cell)%GradS_BF(iDoF1,iGauss) .DotP. K0)
               End Do
            End Do
            Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do
      
         DeAllocate(RHSloc)
         DeAllocate(e0loc)
         DeAllocate(alphaloc)
         !flops = 5 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      End If
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageDispInelasticStrainRHSSetVertex2

#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageDispInelasticStrainRHSSetCell"
!!!
!!!  
!!!  MEF90GradDamageDispInelasticStrainRHSSetCell: Contribution of a cell-based inelastic strain e0
!!!                                       to a RHS in a cell set
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90GradDamageDispInelasticStrainRHSSetCell(RHS,mesh,meshMatS,meshScal,e0,alpha,eta,cellIS,elem,elemType,elemScal,elemScalType,ierr)
      Type(SectionReal),Intent(IN)                       :: RHS
      Type(DM),Intent(IN)                                :: mesh,meshMatS,meshScal
      Type(SectionReal),Intent(IN)                       :: e0
      Type(SectionReal),Intent(IN)                       :: alpha
      PetscReal,Intent(IN)                               :: eta
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elem
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elemScal
      Type(MEF90Element_Type),Intent(IN)                 :: elemType,elemScalType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      Type(MEF90_MATS)                                   :: e0_elem
      PetscReal,Dimension(:),Pointer                     :: e0_loc
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscReal,Dimension(:),Pointer                     :: RHSloc,alphaloc
      PetscReal                                          :: alphaElem,SalphaElem
      PetscInt                                           :: iDoF1,iGauss,i
      PetscLogDouble                                     :: flops
           
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(RHSloc(elemType%numDof))
         Allocate(alphaloc(elemScalType%numDof))
         Do cell = 1,size(cellID)      
            RHSloc = 0.0_Kr
            Call SectionRealRestrict(e0,cellID(cell),e0_loc,ierr);CHKERRQ(ierr)
            Call SectionRealRestrictClosure(alpha,meshScal,cellID(cell),elemScalType%numDof,alphaLoc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               alphaElem = 0.0_Kr
               Do iDoF1 = 1, elemScalType%numDof
                  alphaElem = alphaElem + alphaLoc(iDoF1) * elemScal(cell)%BF(iDoF1,iGauss)
               End Do ! iDoF1
               SalphaElem = eta + (1.0_Kr - alphaElem)**2
               e0_elem = e0_loc
               Do iDoF1 = 1,elemType%numDof
                  RHSloc(iDoF1) = RHSloc(iDoF1) - elem(cell)%Gauss_C(iGauss) * SalphaElem * &
                                 (elem(cell)%GradS_BF(iDoF1,iGauss) .DotP. e0_elem)
               End Do
            End Do
            Call SectionRealRestore(e0,cellID(cell),e0_loc,ierr);CHKERRQ(ierr)
            Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do
         DeAllocate(alphaloc)
         DeAllocate(RHSloc)
         !flops = 3 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      End If
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageDispInelasticStrainRHSSetCell
End Module MEF90_APPEND(m_MEF90_GradDamageImplementation_,MEF90_DIM)D
