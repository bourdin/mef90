#include "mef90.inc"
module MEF90_APPEND(m_MEF90_MassMatrixImplementation_,MEF90_ELEMENTTYPE)
#include "petsc/finclude/petsc.h"
use m_MEF90_Parameters
use m_MEF90_Utils
use m_MEF90_LinAlg
use m_MEF90_Elements
use m_MEF90_Ctx
use m_MEF90_DMPlex
implicit none(type, external)

private
public :: MEF90_MassMatrixAssembleSet

contains
!!!
!!!
!!!  MEF90_MassMatrixAssembleSet: Assemble the contribution of a cell / face / edge set to
!!!                               the mass matrix associated with the element type elemType
!!!
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022 Blaise Bourdin bourdin@mcmaster.ca
!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90_MassMatrixAssembleSet"

subroutine MEF90_MassMatrixAssembleSet(M, dm, setType, setID, elem, elemType, ierr)
   type(tMat), intent(IN)                           :: M
   type(tDM), intent(IN)                            :: dm
   PetscEnum, intent(IN)                            :: setType
   PetscInt                                         :: setID
   type(MEF90_ELEMENTTYPE), dimension(:), pointer   :: elem
   type(MEF90ElementType), intent(IN)               :: elemType
   PetscErrorCode, intent(OUT)                      :: ierr

   type(tIS)                                        :: setPointIS
   PetscInt, dimension(:), pointer                  :: setPointID
   PetscInt                                         :: point
   PetscReal, dimension(:), pointer                 :: MatElem
   PetscInt                                         :: iDoF1, iDoF2, numDof
   PetscInt                                         :: iGauss, numGauss

   PetscCall(DMGetStratumIS(dm, MEF90SetLabelName(setType), setID, setPointIS, ierr))
   PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
   if (size(setPointID) > 0) then
         !! This is really misleading: elemType doesn't know the number of component since we now use the
         !! same elemType for scalar and Vect elements, elem%numDof is NOT the number of dof...
         !! Maybe I need to change this to the old behaviour.
      numDof = size(elem(1)%BF(:, 1))
      numGauss = size(elem(1)%Gauss_C)
      allocate (MatElem(numDof**2))
      do point = 1, size(setPointID)
         MatElem = 0.0_kr
         do iGauss = 1, numGauss
            do iDoF1 = 1, numDof
               do iDoF2 = 1, numDof
                  MatElem((iDoF1 - 1) * numDof + iDof2) = MatElem((iDoF1 - 1) * numDof + iDof2) &
                                                          + elem(point)%Gauss_C(iGauss) * (elem(point)%BF(iDoF1, iGauss) * elem(point)%BF(iDoF2, iGauss))
               end do ! iDoF2
            end do ! iDoF1
         end do ! iGauss
         PetscCall(DMPlexMatSetClosure(dm, PETSC_NULL_SECTION, PETSC_NULL_SECTION, M, setPointID(point), MatElem, ADD_VALUES, ierr))
      end do ! point
      deallocate (MatElem)
   end if
   PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
   PetscCall(ISDestroy(setPointIS, ierr))
end subroutine MEF90_MassMatrixAssembleSet
end module MEF90_APPEND(m_MEF90_MassMatrixImplementation_,MEF90_ELEMENTTYPE)
