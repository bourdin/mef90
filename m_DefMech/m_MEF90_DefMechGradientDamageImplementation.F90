#include "../MEF90/mef90.inc"
Module MEF90_APPEND(m_MEF90_DefMechGradientDamageImplementation_,MEF90_DIM)D
#include "finclude/petscdef.h"
   Use m_MEF90
   Use m_MEF90_DefMechCtx
   IMPLICIT NONE
   Private

   Public :: MEF90DefMechGradientDamageBilinearFormSet

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechGradientDamageBilinearFormSet"
!!!
!!!  
!!!  MEF90DefMechGradientDamageBilinearFormSet:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechGradientDamageBilinearFormSet(K,mesh,meshDamage,cellIS,A,eta,alphaSec,elem,elemType,elemDamage,elemTypeDamage,ierr)
      Type(Mat),Intent(IN)                               :: K 
      Type(DM),Intent(IN)                                :: mesh,meshDamage
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_TENS4OS),Intent(IN)                     :: A
      Type(SectionReal),Intent(IN)                       :: alphaSec
      PetscReal,Intent(IN)                               :: eta
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elem
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elemDamage
      Type(MEF90Element_Type),Intent(IN)                 :: elemType,elemTypeDamage
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscReal,Dimension(:,:),Pointer                   :: MatElem
      PetscReal,Dimension(:),Pointer                     :: alphaLoc
      PetscReal                                          :: alphaElem,SElem
      PetscInt                                           :: iDoF1,iDoF2,iGauss
      Type(SectionReal)                                  :: defaultSection
      PetscLogDouble                                     :: flops
     
      Call DMMeshGetSectionReal(mesh,'default',defaultSection,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(alphaLoc(elemTypeDamage%numDof))
         Allocate(MatElem(elemType%numDof,elemType%numDof))
         Do cell = 1,size(cellID)   
            MatElem = 0.0_Kr
            Call SectionRealRestrictClosure(alphaSec,meshDamage,cellID(cell),elemTypeDamage%numDof,alphaLoc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               alphaElem = 0.0_Kr
               Do iDoF1 = 1, elemTypeDamage%numDof
                  alphaElem = alphaElem + alphaLoc(iDoF1) * elemDamage(cell)%BF(iDoF1,iGauss)
               End Do ! iDoF1
               SElem = eta + (1.0_Kr - alphaElem)**2
               Do iDoF1 = 1,elemType%numDof
                  Do iDoF2 = 1,elemType%numDof
                     MatElem(iDoF2,iDoF1) = MatElem(iDoF2,iDoF1) + elem(cell)%Gauss_C(iGauss) * SElem * &
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
   End Subroutine MEF90DefMechGradientDamageBilinearFormSet


End Module MEF90_APPEND(m_MEF90_DefMechGradientDamageImplementation_,MEF90_DIM)D
