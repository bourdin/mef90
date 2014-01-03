#include "mef90.inc"
Module MEF90_APPEND(m_MEF_MassMatrixImplementation_,MEF90_ELEMENTTYPE)
#include "finclude/petscdef.h"
   Use m_MEF_LinAlg
   Use m_MEF_Parameters
   Use m_MEF_Elements
   Use m_MEF_Utils
   Use petsc
   IMPLICIT NONE

   Private   
   Public :: MEF90_MassMatrixAssembleSet   

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_MassMatrixAssembleSet"
   Subroutine MEF90_MassMatrixAssembleSet(M,mesh,U,cellIS,scaling,elem,elemType,ierr)
      Type(Mat),Intent(IN)                            :: M
      Type(DM),Intent(IN)                             :: mesh
      Type(SectionReal),Intent(IN)                    :: U
      Type(IS),Intent(IN)                             :: cellIS
      PetscReal,Intent(IN)                            :: scaling
      Type(MEF90_ELEMENTTYPE), Dimension(:), Pointer  :: elem
      Type(MEF90Element_Type),Intent(IN)              :: elemType
      PetscErrorCode,Intent(OUT)                      :: ierr
      
      PetscInt,Dimension(:),Pointer                   :: cellID
      PetscInt                                        :: cell
      PetscReal,Dimension(:,:),Pointer                :: MatElem
      PetscInt                                        :: iDoF1,iDoF2,iGauss
      PetscLogDouble                                  :: flops
     
      flops = 0
      Allocate(MatElem(elemType%numDof,elemType%numDof),stat=ierr)

      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Do cell = 1,size(elem)      
         MatElem = 0.0_Kr
         Do iGauss = 1,size(elem(cell)%Gauss_C)
            Do iDoF1 = 1,elemType%numDof
               Do iDoF2 = 1,elemType%numDof
                  MatElem(iDoF2,iDoF1) = MatElem(iDoF2,iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                        (elem(cell)%BF(iDoF1,iGauss) * elem(cell)%BF(iDoF2,iGauss) )
               End Do
            End Do
         End Do
         MatElem = MatElem * scaling
         Call DMmeshAssembleMatrix(M,mesh,U,cellID(cell),MatElem,ADD_VALUES,ierr);CHKERRQ(ierr)
      End Do
      flops = elemType%numDof**2 * (3 * size(elem(1)%Gauss_C) + 1) * size(cellID)
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      DeAllocate(MatElem,stat=ierr)
   End Subroutine MEF90_MassMatrixAssembleSet   
End Module MEF90_APPEND(m_MEF_MassMatrixImplementation_,MEF90_ELEMENTTYPE)