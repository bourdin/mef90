#include "mef90.inc"
module MEF90_APPEND(m_MEF90_NormsImplementation_,MEF90_ELEMENTTYPE)

#include "petsc/finclude/petsc.h"
use m_MEF90_Parameters
use m_MEF90_Utils
use m_MEF90_LinAlg
use m_MEF90_Elements
use m_MEF90_Ctx
use m_MEF90_DMPlex
implicit none(type, external)

private
#if MEF90_ELEMENTTYPE_SCALAR
public :: MEF90L2DotProductSet, MEF90H1DotProductSet, MEF90L2NormSet
#else
public :: MEF90L2DotProductSet, MEF90H1DotProductSet, MEF90H1symDotProductSet, MEF90L2NormSet
#endif

contains
!!! author: Blaise Bourdin (2022, bourdin@mcmaster.ca)
!!!
!!!  MEF90L2DotProductSet: Assemble and add the contribution one processor to the L2 dot products of Vec U and V
!!!                   on a cell / face / edge set, interpolated with element type elemType
!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90L2DotProductSet"

subroutine MEF90L2DotProductSet(myDotProductSet, U, V, setType, setID, elem, elemType, ierr)
   PetscReal, intent(OUT)                           :: myDotProductSet
   type(tVec), intent(IN)                           :: U, V
   PetscEnum, intent(IN)                            :: setType
   PetscInt                                        :: setID
   type(MEF90_ELEMENTTYPE), dimension(:), pointer  :: elem
   type(MEF90ElementType), intent(IN)               :: elemType
   PetscErrorCode, intent(INOUT)                    :: ierr

   PetscReal, dimension(:), pointer                  :: Uloc, Vloc
#if MEF90_ELEMENTTYPE_SCALAR
   PetscReal                                       :: UGauss, Vgauss
#else
   type(MEF90_VECT)                                :: UGauss, VGauss
#endif
   type(tDM)                                       :: dm
   type(tIS)                                       :: setPointIS
   PetscInt, dimension(:), pointer                   :: setPointID
   PetscInt                                        :: point
   PetscInt                                        :: iDoF1, numDof
   PetscInt                                        :: iGauss, numGauss

   PetscCall(VecGetDM(U, dm, ierr))
   PetscCall(DMGetStratumIS(dm, MEF90SetLabelName(setType), setID, setPointIS, ierr))
   PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
   if (size(setPointID) > 0) then
         !! This is really misleading: elemType doesn't know the number of component since we now use the
         !! same elemType for scalar and Vect elements, elem%numDof is NOT the number of dof...
         !! Maybe I need to change this to the old behaviour.
      numDof = size(elem(1)%BF(:, 1))
      numGauss = size(elem(1)%Gauss_C)
      allocate (Uloc(numDof), source=0.0_kr)
      allocate (Vloc(numDof), source=0.0_kr)
      do point = 1, size(setPointID)
         PetscCall(DMPlexVecGetClosure(dm, PETSC_NULL_SECTION, U, setPointID(point), PETSC_NULL_INTEGER, Uloc, ierr))
         PetscCall(DMPlexVecGetClosure(dm, PETSC_NULL_SECTION, V, setPointID(point), PETSC_NULL_INTEGER, Vloc, ierr))

         do iGauss = 1, size(elem(point)%Gauss_C)
            UGauss = 0.0_kr
            VGauss = 0.0_kr
            do iDoF1 = 1, numDof
               UGauss = UGauss + elem(point)%BF(iDoF1, iGauss) * Uloc(iDof1)
               VGauss = VGauss + elem(point)%BF(iDoF1, iGauss) * Vloc(iDof1)
            end do ! iDoF1
            myDotProductSet = myDotProductSet + (UGauss * VGauss) * elem(point)%Gauss_C(iGauss)
         end do ! iGauss
         PetscCall(DMPlexVecRestoreClosure(dm, PETSC_NULL_SECTION, V, setPointID(point), PETSC_NULL_INTEGER, Vloc, ierr))
         PetscCall(DMPlexVecRestoreClosure(dm, PETSC_NULL_SECTION, U, setPointID(point), PETSC_NULL_INTEGER, Uloc, ierr))
      end do ! point
      ! Flop computation is different for scalar and Vec
      deallocate (Vloc, stat=ierr)
      deallocate (Uloc, stat=ierr)
   end if
   PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
   PetscCall(ISDestroy(setPointIS, ierr))
end subroutine MEF90L2DotProductSet

!!! author: Blaise Bourdin (2022, bourdin@mcmaster.ca)
!!!
!!!  MEF90H1DotProductSet: Assemble and add the contribution one processor to the H1 dot products of Vec U and V
!!!                   on a cell / face / edge set, interpolated with element type elemType
!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90H1DotProductSet"

subroutine MEF90H1DotProductSet(myDotProductSet, U, V, setType, setID, elem, elemType, ierr)
   PetscReal, intent(OUT)                           :: myDotProductSet
   type(tVec), intent(IN)                           :: U, V
   PetscEnum, intent(IN)                            :: setType
   PetscInt                                        :: setID
   type(MEF90_ELEMENTTYPE), dimension(:), pointer  :: elem
   type(MEF90ElementType), intent(IN)               :: elemType
   PetscErrorCode, intent(INOUT)                    :: ierr

   PetscReal, dimension(:), pointer                  :: Uloc, Vloc
#if MEF90_ELEMENTTYPE_SCALAR
   type(MEF90_VECT)                                :: GradUGauss, GradVgauss
#else
   type(MEF90_MAT)                                 :: GradUGauss, GradVGauss
#endif
   type(tDM)                                       :: dm
   type(tIS)                                       :: setPointIS
   PetscInt, dimension(:), pointer                   :: setPointID
   PetscInt                                        :: point
   PetscInt                                        :: iDoF1, numDof
   PetscInt                                        :: iGauss, numGauss

   PetscCall(VecGetDM(U, dm, ierr))
   PetscCall(DMGetStratumIS(dm, MEF90SetLabelName(setType), setID, setPointIS, ierr))
   PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
   if (size(setPointID) > 0) then
         !! This is really misleading: elemType doesn't know the number of component since we now use the
         !! same elemType for scalar and Vect elements, elem%numDof is NOT the number of dof...
         !! Maybe I need to change this to the old behaviour.
      numDof = size(elem(1)%BF(:, 1))
      numGauss = size(elem(1)%Gauss_C)
      allocate (Uloc(numDof), source=0.0_kr)
      allocate (Vloc(numDof), source=0.0_kr)
      do point = 1, size(setPointID)
         PetscCall(DMPlexVecGetClosure(dm, PETSC_NULL_SECTION, U, setPointID(point), PETSC_NULL_INTEGER, Uloc, ierr))
         PetscCall(DMPlexVecGetClosure(dm, PETSC_NULL_SECTION, V, setPointID(point), PETSC_NULL_INTEGER, Vloc, ierr))
         do iGauss = 1, size(elem(point)%Gauss_C)
            GradUGauss = 0.0_kr
            GradVGauss = 0.0_kr
            do iDoF1 = 1, numDof
               GradUGauss = GradUGauss + elem(point)%Grad_BF(iDoF1, iGauss) * Uloc(iDof1)
               GradVGauss = GradVGauss + elem(point)%Grad_BF(iDoF1, iGauss) * Vloc(iDof1)
            end do ! iDoF1
            myDotProductSet = myDotProductSet + (GradUGauss.DotP.GradVGauss) * elem(point)%Gauss_C(iGauss)
         end do ! iGauss
         PetscCall(DMPlexVecRestoreClosure(dm, PETSC_NULL_SECTION, V, setPointID(point), PETSC_NULL_INTEGER, Vloc, ierr))
         PetscCall(DMPlexVecRestoreClosure(dm, PETSC_NULL_SECTION, U, setPointID(point), PETSC_NULL_INTEGER, Uloc, ierr))
      end do ! point
      ! Flop computation is different for scalar and Vec
      deallocate (Uloc, stat=ierr)
      deallocate (Vloc, stat=ierr)
   end if
   PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
   PetscCall(ISDestroy(setPointIS, ierr))
end subroutine MEF90H1DotProductSet

!!! author: Blaise Bourdin (2022, bourdin@mcmaster.ca)
!!!
!!!  MEF90H1SymDotProductSet: Assemble and add the contribution one processor to the H1-sym dot products of Vec U and V
!!!                   on a cell / face / edge set, interpolated with element type elemType
!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90H1SymDotProductSet"
#if MEF90_ELEMENTTYPE_VECT

subroutine MEF90H1SymDotProductSet(myDotProductSet, U, V, setType, setID, elem, elemType, ierr)
   PetscReal, intent(OUT)                           :: myDotProductSet
   type(tVec), intent(IN)                           :: U, V
   PetscEnum, intent(IN)                            :: setType
   PetscInt                                        :: setID
   type(MEF90_ELEMENTTYPE), dimension(:), pointer  :: elem
   type(MEF90ElementType), intent(IN)               :: elemType
   PetscErrorCode, intent(INOUT)                    :: ierr

   PetscReal, dimension(:), pointer                  :: Uloc, Vloc
   type(MEF90_MATS)                                :: GradSUGauss, GradSVGauss
   type(tDM)                                       :: dm
   type(tIS)                                       :: setPointIS
   PetscInt, dimension(:), pointer                   :: setPointID
   PetscInt                                        :: point
   PetscInt                                        :: iDoF1, numDof
   PetscInt                                        :: iGauss, numGauss

   PetscCall(VecGetDM(U, dm, ierr))
   PetscCall(DMGetStratumIS(dm, MEF90SetLabelName(setType), setID, setPointIS, ierr))
   PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
   if (size(setPointID) > 0) then
         !! This is really misleading: elemType doesn't know the number of component since we now use the
         !! same elemType for scalar and Vect elements, elem%numDof is NOT the number of dof...
         !! Maybe I need to change this to the old behaviour.
      numDof = size(elem(1)%BF(:, 1))
      numGauss = size(elem(1)%Gauss_C)
      allocate (Uloc(numDof), source=0.0_kr)
      do point = 1, size(setPointID)
         PetscCall(DMPlexVecGetClosure(dm, PETSC_NULL_SECTION, U, setPointID(point), PETSC_NULL_INTEGER, Uloc, ierr))
         PetscCall(DMPlexVecGetClosure(dm, PETSC_NULL_SECTION, V, setPointID(point), PETSC_NULL_INTEGER, Vloc, ierr))
         do iGauss = 1, size(elem(point)%Gauss_C)
            GradSUGauss = 0.0_kr
            GradSVGauss = 0.0_kr
            do iDoF1 = 1, numDof
               GradSUGauss = GradSUGauss + elem(point)%GradS_BF(iDoF1, iGauss) * Uloc(iDof1)
               GradSVGauss = GradSVGauss + elem(point)%GradS_BF(iDoF1, iGauss) * Vloc(iDof1)
            end do ! iDoF1
            myDotProductSet = myDotProductSet + (GradSUGauss.DotP.GradSVGauss) * elem(point)%Gauss_C(iGauss)
         end do ! iGauss
         PetscCall(DMPlexVecRestoreClosure(dm, PETSC_NULL_SECTION, V, setPointID(point), PETSC_NULL_INTEGER, Vloc, ierr))
         PetscCall(DMPlexVecRestoreClosure(dm, PETSC_NULL_SECTION, U, setPointID(point), PETSC_NULL_INTEGER, Uloc, ierr))
      end do ! point
      ! Flop computation is different for scalar and Vec
      deallocate (Uloc, stat=ierr)
   end if
   PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
   PetscCall(ISDestroy(setPointIS, ierr))
end subroutine MEF90H1SymDotProductSet
#endif
!!! author: Blaise Bourdin (2022, bourdin@mcmaster.ca)
!!!
!!!  MEF90L2NormSet: Assemble and add the contribution one processor to the L2 norm of a Vect U
!!!                   on a cell / face / edge set, interpolated with element type elemType
!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90L2NormSet"

subroutine MEF90L2NormSet(myNormSet, U, setType, setID, elem, elemType, ierr)
   PetscReal, intent(OUT)                           :: myNormSet
   type(tVec), intent(IN)                           :: U
   PetscEnum, intent(IN)                            :: setType
   PetscInt                                        :: setID
   type(MEF90_ELEMENTTYPE), dimension(:), pointer  :: elem
   type(MEF90ElementType), intent(IN)               :: elemType
   PetscErrorCode, intent(INOUT)                    :: ierr

   PetscReal, dimension(:), pointer                  :: Uloc
#if MEF90_ELEMENTTYPE_SCALAR
   PetscReal                                       :: UGauss
#else
   type(MEF90_VECT)                                :: UGauss
#endif
   type(tDM)                                       :: dm
   type(tIS)                                       :: setPointIS
   PetscInt, dimension(:), pointer                   :: setPointID
   PetscInt                                        :: point
   PetscInt                                        :: iDoF1, numDof
   PetscInt                                        :: iGauss, numGauss

   PetscCall(VecGetDM(U, dm, ierr))
   PetscCall(DMGetStratumIS(dm, MEF90SetLabelName(setType), setID, setPointIS, ierr))
   PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
   if (size(setPointID) > 0) then
            !! This is really misleading: elemType doesn't know the number of component since we now use the
            !! same elemType for scalar and Vect elements, elem%numDof is NOT the number of dof...
            !! Maybe I need to change this to the old behaviour.
      numDof = size(elem(1)%BF(:, 1))
      numGauss = size(elem(1)%Gauss_C)
      allocate (Uloc(numDof), source=0.0_kr)
      do point = 1, size(setPointID)
         PetscCall(DMPlexVecGetClosure(dm, PETSC_NULL_SECTION, U, setPointID(point), PETSC_NULL_INTEGER, Uloc, ierr))
         do iGauss = 1, size(elem(point)%Gauss_C)
            UGauss = 0.0_kr
            do iDoF1 = 1, numDof
               UGauss = UGauss + elem(point)%BF(iDoF1, iGauss) * Uloc(iDof1)
            end do ! iDoF1
            myNormSet = myNormSet + (UGauss * UGauss) * elem(point)%Gauss_C(iGauss)
         end do ! iGauss
         PetscCall(DMPlexVecRestoreClosure(dm, PETSC_NULL_SECTION, U, setPointID(point), PETSC_NULL_INTEGER, Uloc, ierr))
      end do ! point
      deallocate (Uloc, stat=ierr)
   end if
   PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
   PetscCall(ISDestroy(setPointIS, ierr))
end subroutine MEF90L2NormSet
end module MEF90_APPEND(m_MEF90_NormsImplementation_,MEF90_ELEMENTTYPE)
