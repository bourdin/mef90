Module m_MEF90_Elements
#include "petsc/finclude/petsc.h"
   Use m_MEF90_Elements_Scal
   Use m_MEF90_Elements_Vect
   Use m_MEF90_Elements_Type
   Use m_MEF90_LinAlg
   Use m_MEF90_Utils
   Use m_MEF90_Parameters
   Use petsc
   IMPLICIT NONE

   !Private   
   Public :: MEF90ElementCreate
   Public :: MEF90ElementGetType,MEF90ElementGetTypeBoundary
   Public :: MEF90ElementDestroy

   Interface MEF90ElementCreate
      Module Procedure Element2DScalInitSet,Element2DVectInitSet, &
                       Element3DScalInitSet,Element3DVectInitSet
   End Interface MEF90ElementCreate
   
   Interface MEF90ElementDestroy
      Module Procedure Element2DScalDestroy,Element2DVectDestroy,&
                       Element3DScalDestroy,Element3DVectDestroy,&
                       Element2DScalDestroySet,Element2DVectDestroySet,& 
                       Element3DScalDestroySet,Element3DVectDestroySet
   End Interface MEF90ElementDestroy
      
Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90ElementGetType"
!!!
!!!
!!!  MEF90ElementGetType: Return an element given a family and order
!!!
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90ElementGetType(elemFamily,order,cellType,elemType,ierr)
      PetscEnum,Intent(IN)                             :: elemFamily
      PetscInt,Intent(IN)                              :: order
      DMPolytopeType,Intent(IN)                        :: cellType
      Type(MEF90ElementType),Intent(OUT)               :: elemType
      PetscErrorCode,Intent(INOUT)                     :: ierr

      elemType = MEF90_NULL_ELEMENT
      Select Case(elemFamily)
      Case(MEF90ElementFamilyLagrange)
         Select Case(cellType)
         Case(DM_POLYTOPE_TRIANGLE)
            Select Case(order)
            Case(0)
               elemType = MEF90P0Lagrange2D
            Case(1)
               elemType = MEF90P1Lagrange2D
            Case(2)
               elemType = MEF90P2Lagrange2D
            ! Case Default
            !    Write(*,*) __FUNCT__,': Unimplemented order',order
            !    ierr = PETSC_ERR_SUP
            End Select ! order
         Case(DM_POLYTOPE_TETRAHEDRON)
            Select Case(order)
            Case(0)
               elemType = MEF90P0Lagrange3D
            Case(1)
               elemType = MEF90P1Lagrange3D
            Case(2)
               elemType = MEF90P2Lagrange3D
            ! Case Default
            !    Write(*,*) __FUNCT__,': Unimplemented order',order
            !    ierr = PETSC_ERR_SUP
            End Select ! order
         Case(DM_POLYTOPE_QUADRILATERAL)
            Select Case(order)
            Case(0)
               elemType = MEF90P0Lagrange2D
            Case(1)
               elemType = MEF90Q1Lagrange2D
            Case(2)
               elemType = MEF90Q2Lagrange2D
            ! Case Default
            !    Write(*,*) __FUNCT__,': Unimplemented order',order
            !    ierr = PETSC_ERR_SUP
            End Select ! order
         Case(DM_POLYTOPE_HEXAHEDRON)
            Select Case(order)
            Case(0)
               elemType = MEF90P0Lagrange3D
            Case(1)
               elemType = MEF90Q1Lagrange3D
            Case(2)
               elemType = MEF90Q2Lagrange3D
            ! Case Default
            !    Write(*,*) __FUNCT__,': Unimplemented order',order
            !    ierr = PETSC_ERR_SUP
            End Select ! order
         ! Case Default
         !    Write(*,*) __FUNCT__,': Unknown cell type (see $PETSC_DIR/srcdm/f90-mod/petscdm.h)',cellType
         !    ierr = PETSC_ERR_SUP
         End Select ! cellType
      ! Case Default
      ! Write(*,*) __FUNCT__,': Unknown element family',elemFamily
      ! ierr = PETSC_ERR_SUP
      End Select
   End Subroutine MEF90ElementGetType

#undef __FUNCT__
#define __FUNCT__ "MEF90ElementGetTypeBoundary"
!!!
!!!
!!!  MEF90ElementGetTypeBoundary: Return an element given a family and order
!!!
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90ElementGetTypeBoundary(elemFamily,order,cellType,elemType,ierr)
      PetscEnum,Intent(IN)                             :: elemFamily
      PetscInt,Intent(IN)                              :: order
      PetscEnum,Intent(IN)                             :: cellType
      Type(MEF90ElementType),Intent(OUT)               :: elemType
      PetscErrorCode,Intent(INOUT)                     :: ierr

      elemType = MEF90_NULL_ELEMENT
      Select Case(elemFamily)
      Case(MEF90ElementFamilyLagrange)
         Select Case(cellType)
         Case(DM_POLYTOPE_SEGMENT)
            Select Case(order)
            Case(0)
               elemType = MEF90P0Lagrange2DBoundary
            Case(1)
               elemType = MEF90P1Lagrange2DBoundary
            Case(2)
               elemType = MEF90P2Lagrange2DBoundary
            ! Case Default
            !    Write(*,*) __FUNCT__,': Unimplemented order',order
            !    ierr = PETSC_ERR_SUP
            End Select ! order
         Case(DM_POLYTOPE_TRIANGLE)
            Select Case(order)
            Case(0)
               elemType = MEF90P0Lagrange3DBoundary
            Case(1)
               elemType = MEF90P1Lagrange3DBoundary
            Case(2)
               elemType = MEF90P2Lagrange3DBoundary
            ! Case Default
            !    Write(*,*) __FUNCT__,': Unimplemented order',order
            !    ierr = PETSC_ERR_SUP
            End Select ! order
         Case(DM_POLYTOPE_QUADRILATERAL)
            Select Case(order)
            Case(0)
               elemType = MEF90P0Lagrange3DBoundary
            Case(1)
               elemType = MEF90Q1Lagrange3DBoundary
            Case(2)
               elemType = MEF90Q2Lagrange3DBoundary
            ! Case Default
            !    Write(*,*) __FUNCT__,': Unimplemented order',order
            !    ierr = PETSC_ERR_SUP
            End Select ! order
         ! Case Default
         !    Write(*,*) __FUNCT__,': Unknown cell type (see $PETSC_DIR/src/dm/f90-mod/petscdm.h)',cellType
         !    ierr = PETSC_ERR_SUP
         End Select ! cellType
      ! Case Default
      ! Write(*,*) __FUNCT__,': Unknown element family',elemFamily
      ! ierr = PETSC_ERR_SUP
      End Select
   End Subroutine MEF90ElementGetTypeBoundary
   
End Module m_MEF90_Elements
