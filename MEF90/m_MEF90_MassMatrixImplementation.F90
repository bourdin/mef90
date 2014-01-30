#include "mef90.inc"
Module MEF90_APPEND(m_MEF90_MassMatrixImplementation_,MEF90_ELEMENTTYPE)
#include "finclude/petscdef.h"
   Use m_MEF90_LinAlg
   Use m_MEF90_Parameters
   Use m_MEF90_Elements
   Use m_MEF90_Utils
   Use petsc
   IMPLICIT NONE

   Private   
   Public :: MEF90_MassMatrixAssembleSet   

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_MassMatrixAssembleSet"
   Subroutine MEF90_MassMatrixAssembleSet(M,mesh,cellIS,scaling,elem,elemType,ierr)
      Type(Mat),Intent(IN)                            :: M
      Type(DM),Intent(IN)                             :: mesh
      Type(IS),Intent(IN)                             :: cellIS
      PetscReal,Intent(IN)                            :: scaling
      Type(MEF90_ELEMENTTYPE), Dimension(:), Pointer  :: elem
      Type(MEF90Element_Type),Intent(IN)              :: elemType
      PetscErrorCode,Intent(OUT)                      :: ierr
      
      Type(SectionReal)                               :: defaultSection
      PetscInt,Dimension(:),Pointer                   :: cellID
      PetscInt                                        :: cell
      PetscReal,Dimension(:,:),Pointer                :: MatElem
      PetscInt                                        :: iDoF1,iDoF2,iGauss
      PetscLogDouble                                  :: flops
     
      flops = 0
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
                                           (elem(cell)%BF(iDoF1,iGauss) * elem(cell)%BF(iDoF2,iGauss) )
                  End Do ! iDoF2
               End Do ! iDoF1
            End Do ! iGauss
            MatElem = MatElem * scaling
            Call DMmeshAssembleMatrix(M,mesh,defaultSection,cellID(cell),MatElem,ADD_VALUES,ierr);CHKERRQ(ierr)
         End Do ! cell
         flops = 3 * elemType%numDof**2 * size(elem(1)%Gauss_C) * size(cellID)
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(MatElem,stat=ierr)
      End If 
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(defaultSection,ierr);CHKERRQ(ierr)
      
   End Subroutine MEF90_MassMatrixAssembleSet   
End Module MEF90_APPEND(m_MEF90_MassMatrixImplementation_,MEF90_ELEMENTTYPE)