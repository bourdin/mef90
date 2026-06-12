module m_MEF90_Elements
#include "petsc/finclude/petsc.h"
   use m_MEF90_Parameters
   use m_MEF90_LinAlg
   use petscdmplex
   use, intrinsic :: iso_c_binding
   implicit none(type, external)

   private
   public :: MEF90ElementCreate
   public :: MEF90ElementDestroy

   public :: MEF90ElementType
   public :: MEF90Element2DVect, MEF90Element2DScal
   public :: MEF90Element3DVect, MEF90Element3DScal
   public :: MEF90ElementGetType, MEF90ElementGetTypeBoundary
   public :: MEF90ElementFamilyLagrange
   public :: MEF90ElementFamilyList
   public :: MEF90ElementsInitialize_Private

   type MEF90ElementType
      ! name is the element name in english language
      character(len=MEF90MXSTRLEN)                 :: Name

      !!! Number of dof at each
      !!!   vertex edge face cell in 3D,
      !!!   vertex edge cell in 2D (last value is unused, i.e. always use depth to locate DoF)
      !!! i.e. point depth -1
      PetscInt, dimension(4)                        :: numDofs
      PetscInt                                     :: numDoF
      !!! numDof = numVertexDof * numVertex + numEdgeDof * numEdge
      !!!          + numFaceDof * numFace + numCellDof
      PetscInt                                     :: dim
      !!! The dimension of the cell associated with the element
      PetscInt                                     :: coDim
      !!! The co-dimension of the cell associated with the element
      PetscInt                                     :: order
      !!! The polynomial order
   end type MEF90ElementType

   enum, bind(c)
      enumerator :: &
         MEF90ElementFamilyLagrange = 0
   end enum

   character(kind=c_char, len=MEF90MXSTRLEN), dimension(4), parameter, public   :: MEF90ElementFamily = [ &
                                                                                   "Lagrange           ", &  ! 0
                                                                                   "MEF90ElementFamily ", &
                                                                                   "prefix_            ", &
                                                                                   c_null_char//"                  "]
   character(len=MEF90MXSTRLEN), dimension(4), protected  :: MEF90ElementFamilyList

   type(MEF90ElementType), parameter, public :: MEF90_NULL_ELEMENT = MEF90ElementType( &
                                                "NULL", &  ! name
                                                [0, 0, 0, 0], 0, &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
                                                0, 0, 0 &  ! dim,codim,order
                                                )

   !!!
   !!! P0 pseudo-elements. These are used for creating sections only
   !!!
   type(MEF90ElementType), parameter, public :: MEF90P0Lagrange2D = MEF90ElementType( &
                                                "MEF90P0Lagrange2D", &  ! name
                                                [0, 0, 1, 0], 1, &  ! numVertexDof,numEdgeDof,numCellDof,unused,numDof
                                                2, 0, 0 &  ! dim,codim,order
                                                )
   type(MEF90ElementType), parameter, public :: MEF90P0Lagrange3D = MEF90ElementType( &
                                                "MEF90P0Lagrange3D", &  ! name
                                                [0, 0, 0, 1], 1, &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
                                                3, 0, 0 &  ! dim,codim,order
                                                )
   type(MEF90ElementType), parameter, public :: MEF90P0Lagrange2DBoundary = MEF90ElementType( &
                                                "MEF90P0Lagrange2DBoundary", &  ! name
                                                [0, 1, 0, 0], 1, &  ! numVertexDof,numEdgeDof,numCellDof,unused,numDof
                                                2, 0, 0 &  ! dim,codim,order
                                                )
   type(MEF90ElementType), parameter, public :: MEF90P0Lagrange3DBoundary = MEF90ElementType( &
                                                "MEF90P0Lagrange3DBoundary", &  ! name
                                                [0, 0, 1, 0], 1, &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
                                                3, 0, 0 &  ! dim,codim,order
                                                )
   !!!
   !!! P1 Simplicial Lagrange elements
   !!!
   type(MEF90ElementType), parameter, public :: MEF90P1Lagrange2D = MEF90ElementType( &
                                                "MEF90P1Lagrange2D", &  ! name
                                                [1, 0, 0, 0], 3, &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
                                                2, 0, 1 &  ! dim,codim,order
                                                )
   type(MEF90ElementType), parameter, public :: MEF90P1Lagrange3D = MEF90ElementType( &
                                                "MEF90P1Lagrange3D", &  ! name
                                                [1, 0, 0, 0], 4, &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
                                                3, 0, 1 &  ! dim,codim,order
                                                )
   type(MEF90ElementType), parameter, public :: MEF90P1Lagrange2DBoundary = MEF90ElementType( &
                                                "MEF90P1Lagrange2DBoundary", &  ! name
                                                [1, 0, 0, 0], 2, &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
                                                2, 1, 1 &  ! dim,codim,order
                                                )
   type(MEF90ElementType), parameter, public :: MEF90P1Lagrange3DBoundary = MEF90ElementType( &
                                                "MEF90P1Lagrange3DBoundary", &  ! name
                                                [1, 0, 0, 0], 3, &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
                                                3, 1, 1 &  ! dim,codim,order
                                                )
   !!!
   !!! P2 Simplicial Lagrange elements
   !!!
   type(MEF90ElementType), parameter, public :: MEF90P2Lagrange2D = MEF90ElementType( &
                                                "MEF90P2Lagrange2D", &  ! name
                                                [1, 1, 0, 0], 6, &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
                                                2, 0, 2 &  ! dim,codim,order
                                                )
   type(MEF90ElementType), parameter, public :: MEF90P2Lagrange3D = MEF90ElementType( &
                                                "MEF90P2Lagrange3D", &  ! name
                                                [1, 1, 0, 0], 10, &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
                                                3, 0, 2 &  ! dim,codim,order
                                                )
   type(MEF90ElementType), parameter, public :: MEF90P2Lagrange2DBoundary = MEF90ElementType( &
                                                "MEF90P2Lagrange2DBoundary", &  ! name
                                                [1, 1, 0, 0], 3, &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
                                                2, 1, 2 &  ! dim,codim,order
                                                )
   type(MEF90ElementType), parameter, public :: MEF90P2Lagrange3DBoundary = MEF90ElementType( &
                                                "MEF90P2Lagrange3DBoundary", &  ! name
                                                [1, 1, 0, 0], 6, &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
                                                3, 1, 2 &  ! dim,codim,order
                                                )
   !!!
   !!! Q1 tensor product Lagrange elements
   !!!
   type(MEF90ElementType), parameter, public :: MEF90Q1Lagrange2D = MEF90ElementType( &
                                                "MEF90Q1Lagrange2D", &  ! name
                                                [1, 0, 0, 0], 4, &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
                                                2, 0, 1 &  ! dim,codim,order
                                                )
   type(MEF90ElementType), parameter, public :: MEF90Q1Lagrange3D = MEF90ElementType( &
                                                "MEF90Q1Lagrange3D", &  ! name
                                                [1, 0, 0, 0], 8, &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
                                                3, 0, 1 &  ! dim,codim,order
                                                )
   type(MEF90ElementType), parameter, public :: MEF90Q1Lagrange2DBoundary = MEF90ElementType( &
                                                "MEF90Q1Lagrange2DBoundary", &  ! name
                                                [1, 0, 0, 0], 2, &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
                                                2, 1, 1 &  ! dim,codim,order
                                                )
   type(MEF90ElementType), parameter, public :: MEF90Q1Lagrange3DBoundary = MEF90ElementType( &
                                                "MEF90Q1Lagrange3DBoundary", &  ! name
                                                [1, 0, 0, 0], 4, &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
                                                3, 1, 1 &  ! dim,codim,order
                                                )
   !!!
   !!! Q2 tensor product Lagrange elements
   !!!
   type(MEF90ElementType), parameter, public :: MEF90Q2Lagrange2D = MEF90ElementType( &
                                                "MEF90Q2Lagrange2D", &  ! name
                                                [1, 1, 0, 1], 9, &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
                                                2, 0, 1 &  ! dim,codim,order
                                                )
   type(MEF90ElementType), parameter, public :: MEF90Q2Lagrange3D = MEF90ElementType( &
                                                "MEF90Q2Lagrange3D", &  ! name
                                                [1, 1, 1, 1], 27, &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
                                                3, 0, 1 &  ! dim,codim,order
                                                )
   type(MEF90ElementType), parameter, public :: MEF90Q2Lagrange2DBoundary = MEF90ElementType( &
                                                "MEF90Q2Lagrange2DBoundary", &  ! name
                                                [1, 1, 0, 0], 3, &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
                                                2, 1, 1 &  ! dim,codim,order
                                                )
   type(MEF90ElementType), parameter, public :: MEF90Q2Lagrange3DBoundary = MEF90ElementType( &
                                                "MEF90Q2Lagrange3DBoundary", &  ! name
                                                [1, 1, 0, 1], 9, &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
                                                3, 1, 1 &  ! dim,codim,order
                                                )
   type MEF90Element2DScal
      PetscReal, dimension(:, :), pointer             :: BF => null()
      type(Vect2D), dimension(:, :), pointer          :: Grad_BF => null()
      PetscReal, dimension(:), pointer               :: Gauss_C => null()
      type(Vect2D)                                 :: outerNormal ! only makes sense for elements of codim 1
   end type MEF90Element2DScal

   type MEF90Element2DVect
      type(Vect2D), dimension(:, :), pointer         :: BF => null()
      type(Mat2D), dimension(:, :), pointer          :: Grad_BF => null()
      type(MatS2D), dimension(:, :), pointer         :: GradS_BF => null()
      PetscReal, dimension(:), pointer               :: Gauss_C => null()
      type(Vect2D)                                 :: outerNormal ! only makes sense for elements of codim 1
   end type MEF90Element2DVect

   type MEF90Element3DScal
      PetscReal, dimension(:, :), pointer             :: BF => null()
      type(Vect3D), dimension(:, :), pointer         :: Grad_BF => null()
      PetscReal, dimension(:), pointer               :: Gauss_C => null()
      type(Vect3D)                                 :: outerNormal ! only makes sense for elements of codim 1
   end type MEF90Element3DScal

   type MEF90Element3DVect
      type(Vect3D), dimension(:, :), pointer         :: BF => null()
      type(Mat3D), dimension(:, :), pointer          :: Grad_BF => null()
      type(MatS3D), dimension(:, :), pointer         :: GradS_BF => null()
      PetscReal, dimension(:), pointer               :: Gauss_C => null()
      type(Vect3D)                                 :: outerNormal ! only makes sense for elements of codim 1
   end type MEF90Element3DVect

   interface MEF90ElementCreate
      module procedure Element2DScalInitSet, Element2DVectInitSet, &
         Element3DScalInitSet, Element3DVectInitSet
   end interface MEF90ElementCreate

   interface MEF90ElementDestroy
      module procedure Element2DScalDestroy, Element2DVectDestroy, &
         Element3DScalDestroy, Element3DVectDestroy, &
         Element2DScalDestroySet, Element2DVectDestroySet, &
         Element3DScalDestroySet, Element3DVectDestroySet
   end interface MEF90ElementDestroy

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90ElementsInitialize_Private"
!!!
!!!
!!!  MEF90ElementsInitialize_Private:
!!!
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90ElementsInitialize_Private(ierr)
      PetscErrorCode, intent(OUT)                   :: ierr

      MEF90ElementFamilyList(1) = 'Lagrange'
      MEF90ElementFamilyList(2) = 'MEF90ElementFamily'
      MEF90ElementFamilyList(3) = '_MEF90ElementFamily'
      MEF90ElementFamilyList(4) = ''
      ierr = 0
   end subroutine MEF90ElementsInitialize_Private

#undef __FUNCT__
#define __FUNCT__ "MEF90ElementGetType"
!!!
!!!
!!!  MEF90ElementGetType: Return an element given a family and order
!!!
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90ElementGetType(elemFamily, order, cellType, elemType, ierr)
      PetscEnum, intent(IN)                             :: elemFamily
      PetscInt, intent(IN)                              :: order
      type(eDMPolytopeType), intent(IN)                 :: cellType
      type(MEF90ElementType), intent(OUT)               :: elemType
      PetscErrorCode, intent(INOUT)                     :: ierr

      elemType = MEF90_NULL_ELEMENT
      select case (elemFamily)
      case (MEF90ElementFamilyLagrange)
         select PetscEnumCase(cellType)
         PetscEnumCase(DM_POLYTOPE_TRIANGLE)
         select case (order)
         case (0)
            elemType = MEF90P0Lagrange2D
         case (1)
            elemType = MEF90P1Lagrange2D
         case (2)
            elemType = MEF90P2Lagrange2D
            ! Case Default
            !    Write(*,*) __FUNCT__,': Unimplemented order',order
            !    ierr = PETSC_ERR_SUP
         end select ! order
         PetscEnumCase(DM_POLYTOPE_TETRAHEDRON)
         select case (order)
         case (0)
            elemType = MEF90P0Lagrange3D
         case (1)
            elemType = MEF90P1Lagrange3D
         case (2)
            elemType = MEF90P2Lagrange3D
            ! Case Default
            !    Write(*,*) __FUNCT__,': Unimplemented order',order
            !    ierr = PETSC_ERR_SUP
         end select ! order
         PetscEnumCase(DM_POLYTOPE_QUADRILATERAL)
         select case (order)
         case (0)
            elemType = MEF90P0Lagrange2D
         case (1)
            elemType = MEF90Q1Lagrange2D
         case (2)
            elemType = MEF90Q2Lagrange2D
            ! Case Default
            !    Write(*,*) __FUNCT__,': Unimplemented order',order
            !    ierr = PETSC_ERR_SUP
         end select ! order
         PetscEnumCase(DM_POLYTOPE_HEXAHEDRON)
         select case (order)
         case (0)
            elemType = MEF90P0Lagrange3D
         case (1)
            elemType = MEF90Q1Lagrange3D
         case (2)
            elemType = MEF90Q2Lagrange3D
            ! Case Default
            !    Write(*,*) __FUNCT__,': Unimplemented order',order
            !    ierr = PETSC_ERR_SUP
         end select ! order
         ! Case Default
         !    Write(*,*) __FUNCT__,': Unknown cell type (see $PETSC_DIR/srcdm/f90-mod/petscdm.h)',cellType
         !    ierr = PETSC_ERR_SUP
      end select ! cellType
      ! Case Default
      ! Write(*,*) __FUNCT__,': Unknown element family',elemFamily
      ! ierr = PETSC_ERR_SUP
      end select
      end subroutine MEF90ElementGetType

#undef __FUNCT__
#define __FUNCT__ "MEF90ElementGetTypeBoundary"
!!!
!!!
!!!  MEF90ElementGetTypeBoundary: Return an element given a family and order
!!!
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

      subroutine MEF90ElementGetTypeBoundary(elemFamily, order, cellType, elemType, ierr)
         PetscEnum, intent(IN)                             :: elemFamily
         PetscInt, intent(IN)                              :: order
         type(eDMPolytopeType), intent(IN)                 :: cellType
         type(MEF90ElementType), intent(OUT)               :: elemType
         PetscErrorCode, intent(INOUT)                     :: ierr

         elemType = MEF90_NULL_ELEMENT
         select case (elemFamily)
         case (MEF90ElementFamilyLagrange)
            select PetscEnumCase(cellType)
            PetscEnumCase(DM_POLYTOPE_SEGMENT)
            select case (order)
            case (0)
               elemType = MEF90P0Lagrange2DBoundary
            case (1)
               elemType = MEF90P1Lagrange2DBoundary
            case (2)
               elemType = MEF90P2Lagrange2DBoundary
               ! Case Default
               !    Write(*,*) __FUNCT__,': Unimplemented order',order
               !    ierr = PETSC_ERR_SUP
            end select ! order
            PetscEnumCase(DM_POLYTOPE_TRIANGLE)
            select case (order)
            case (0)
               elemType = MEF90P0Lagrange3DBoundary
            case (1)
               elemType = MEF90P1Lagrange3DBoundary
            case (2)
               elemType = MEF90P2Lagrange3DBoundary
               ! Case Default
               !    Write(*,*) __FUNCT__,': Unimplemented order',order
               !    ierr = PETSC_ERR_SUP
            end select ! order
            PetscEnumCase(DM_POLYTOPE_QUADRILATERAL)
            select case (order)
            case (0)
               elemType = MEF90P0Lagrange3DBoundary
            case (1)
               elemType = MEF90Q1Lagrange3DBoundary
            case (2)
               elemType = MEF90Q2Lagrange3DBoundary
               ! Case Default
               !    Write(*,*) __FUNCT__,': Unimplemented order',order
               !    ierr = PETSC_ERR_SUP
            end select ! order
            ! Case Default
            !    Write(*,*) __FUNCT__,': Unknown cell type (see $PETSC_DIR/src/dm/f90-mod/petscdm.h)',cellType
            !    ierr = PETSC_ERR_SUP
         end select ! cellType
         ! Case Default
         ! Write(*,*) __FUNCT__,': Unknown element family',elemFamily
         ! ierr = PETSC_ERR_SUP
         end select
      end subroutine MEF90ElementGetTypeBoundary

#undef __FUNCT__
#define __FUNCT__ "Element2DScalInitSet"
!!!
!!!
!!!  Element2DScalInitSet:
!!!
!!!  (c) 2014-2022 Blaise Bourdin bourdin@lsu.edu
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

      subroutine Element2DScalInitSet(dm, cellIS, dElem, dQuadratureOrder, elemType, ierr)
         type(tDM), intent(IN)                             :: dm
         type(tIS), intent(IN)                             :: cellIS
         type(MEF90Element2DScal), dimension(:), pointer    :: dElem
         PetscInt, intent(IN)                              :: dQuadratureOrder
         type(MEF90ElementType), intent(IN)                :: elemType
         PetscErrorCode, intent(OUT)                       :: ierr

         PetscInt                                         :: iELoc
         PetscInt, dimension(:), pointer                    :: CellID
         PetscReal, dimension(:), pointer                 :: v0, BB, BBinv
         PetscReal                                        :: detBBinv
         type(Mat2D)                                      :: Bt
         PetscReal                                        :: length
         PetscReal, dimension(:), pointer                   :: centroid, innerNormal
         type(Vect2D)                                     :: outerNormal

         select case (elemType%name)
         case (MEF90P1Lagrange2D%name, MEF90P2Lagrange2D%name)
            PetscCall(ISGetIndices(CellIS, CellID, ierr))
            allocate (dElem(size(cellID)), stat=ierr)
            if (size(CellID) > 0) then
               allocate (v0(2))
               allocate (BB(4))
               allocate (BBinv(4))
               do iELoc = 1, size(CellID)
                  PetscCall(DMPlexComputeCellGeometryAffineFEM(dm, cellID(iELoc), v0, BB, BBinv, detBBinv, ierr))
                  !!! Petsc uses a reference simplex with vertices at (-1,-1), (1,-1) and (-1,1)
                  !!! Whereas MEF90 uses (0,0), (1,0), (0,1), so we need to rescale the affine transform
                  Bt%XX = BBinv(1) * 0.5_kr; Bt%XY = BBinv(3) * 0.5_kr
                  Bt%YX = BBinv(2) * 0.5_kr; Bt%YY = BBinv(4) * 0.5_kr
                  detBBinv = 4.0_kr * abs(detBBinv)
                  call ElementPLagrange2DScalInit(dElem(iELoc), Bt, detBBinv, elemType%order, dQuadratureOrder, ierr)
               end do
               deallocate (BBinv)
               deallocate (BB)
               deallocate (v0)
            end if
            PetscCall(ISRestoreIndices(CellIS, CellID, ierr))
         case (MEF90P1Lagrange2DBoundary%name, MEF90P2Lagrange2DBoundary%name)
            PetscCall(ISGetIndices(CellIS, CellID, ierr))
            allocate (dElem(size(cellID)), stat=ierr)
            if (size(CellID) > 0) then
               allocate (centroid(2))
               allocate (innerNormal(2))
               do iELoc = 1, size(CellID)
                  PetscCall(DMPlexComputeCellGeometryFVM(dm, cellID(iEloc), length, centroid, innerNormal, ierr))
                  outerNormal = -innerNormal
                  call ElementPLagrange2DBoundaryScalInit(dElem(iELoc), length, outerNormal, elemType%order, dQuadratureOrder, ierr)
               end do
               deallocate (innerNormal)
               deallocate (centroid)
            end if
            PetscCall(ISRestoreIndices(CellIS, CellID, ierr))
            !Case (MEF90Q1Lagrange2D%name,MEF90Q2Lagrange2D%name,MEF90Q1Lagrange2DBoundary%name,MEF90Q2Lagrange2DBoundary%name)
            !   !!! Get quadrature points for the current element using DMPlexComputeCellGeometryFEM
            !   !!! Initialize element
         case Default
            write (*, *) __FUNCT__, ': Element type not implemented yet', elemType%name
            ierr = PETSC_ERR_SUP
         end select
      end subroutine Element2DScalInitSet

#undef __FUNCT__
#define __FUNCT__ "Element2DVectInitSet"
      subroutine Element2DVectInitSet(dm, cellIS, dElem, dQuadratureOrder, elemType, ierr)
         type(tDM), intent(IN)                             :: dm
         type(tIS), intent(IN)                             :: cellIS
         type(MEF90Element2DVect), dimension(:), pointer    :: dElem
         PetscInt, intent(IN)                              :: dQuadratureOrder
         type(MEF90ElementType), intent(IN)                :: elemType
         PetscErrorCode, intent(OUT)                       :: ierr

         PetscInt                                         :: iELoc
         PetscInt, dimension(:), pointer                    :: CellID
         PetscReal, dimension(:), pointer                 :: v0, BB, BBinv
         PetscReal                                        :: detBBinv
         type(Mat2D)                                      :: Bt
         PetscReal                                        :: length
         PetscReal, dimension(:), pointer                   :: centroid, innerNormal
         type(Vect2D)                                     :: outerNormal

         select case (elemType%name)
         case (MEF90P1Lagrange2D%name, MEF90P2Lagrange2D%name)
            PetscCall(ISGetIndices(CellIS, CellID, ierr))
            allocate (dElem(size(cellID)), stat=ierr)
            if (size(CellID) > 0) then
               allocate (v0(2))
               allocate (BB(4))
               allocate (BBinv(4))
               Do_Elem_iE: do iELoc = 1, size(CellID)
                  PetscCall(DMPlexComputeCellGeometryAffineFEM(dm, cellID(iELoc), v0, BB, BBinv, detBBinv, ierr))
                  !!! Petsc uses a reference simplex with vertices at (-1,-1), (1,-1) and (-1,1)
                  !!! Whereas MEF90 uses (0,0), (1,0), (0,1), so we need to rescale the affine transform
                  Bt%XX = BBinv(1) * 0.5_kr; Bt%XY = BBinv(3) * 0.5_kr
                  Bt%YX = BBinv(2) * 0.5_kr; Bt%YY = BBinv(4) * 0.5_kr
                  detBBinv = 4.0_kr * abs(detBBinv)
                  call ElementPLagrange2DVectInit(dElem(iELoc), Bt, detBBinv, elemType%order, dQuadratureOrder, ierr)
               end do Do_Elem_iE
               deallocate (BBinv)
               deallocate (BB)
               deallocate (v0)
            end if
            PetscCall(ISRestoreIndices(CellIS, CellID, ierr))
         case (MEF90P1Lagrange2DBoundary%name, MEF90P2Lagrange2DBoundary%name)
            PetscCall(ISGetIndices(CellIS, CellID, ierr))
            allocate (dElem(size(cellID)), stat=ierr)
            if (size(CellID) > 0) then
               allocate (centroid(2))
               allocate (innerNormal(2))
               do iELoc = 1, size(CellID)
                  PetscCall(DMPlexComputeCellGeometryFVM(dm, cellID(iEloc), length, centroid, innerNormal, ierr))
                  outerNormal = -innerNormal
                  call ElementPLagrange2DBoundaryVectInit(dElem(iELoc), length, outerNormal, elemType%order, dQuadratureOrder, ierr)
               end do
               deallocate (innerNormal)
               deallocate (centroid)
            end if
            PetscCall(ISRestoreIndices(CellIS, CellID, ierr))
            !Case (MEF90Q1Lagrange2D%name,MEF90Q2Lagrange2D%name,MEF90Q1Lagrange2DBoundary%name,MEF90Q2Lagrange2DBoundary%name)
            !   !!! Get quadrature points for the current element using DMPlexComputeCellGeometryFEM
            !   !!! Initialize element
         case Default
            write (*, *) __FUNCT__, ': Element type not implemented yet', elemType%name
            ierr = PETSC_ERR_SUP
         end select
      end subroutine Element2DVectInitSet

#undef __FUNCT__
#define __FUNCT__ "Element3DScalInitSet"
      subroutine Element3DScalInitSet(dm, cellIS, dElem, dQuadratureOrder, elemType, ierr)
         type(tDM), intent(IN)                             :: dm
         type(tIS), intent(IN)                             :: cellIS
         type(MEF90Element3DScal), dimension(:), pointer    :: dElem
         PetscInt, intent(IN)                              :: dQuadratureOrder
         type(MEF90ElementType), intent(IN)                :: elemType
         PetscErrorCode, intent(OUT)                       :: ierr

         PetscInt                                         :: iELoc
         PetscInt, dimension(:), pointer                    :: CellID
         PetscReal, dimension(:), pointer                 :: v0, BB, BBinv
         PetscReal                                        :: detBBinv
         type(Mat3D)                                      :: Bt
         PetscReal                                        :: area
         PetscReal, dimension(:), pointer                   :: centroid, innerNormal
         type(Vect3D)                                     :: outerNormal

         select case (elemType%name)
         case (MEF90P1Lagrange3D%name, MEF90P2Lagrange3D%name)
            PetscCall(ISGetIndices(CellIS, CellID, ierr))
            allocate (dElem(size(cellID)), stat=ierr)
            if (size(CellID) > 0) then
               allocate (v0(3))
               allocate (BB(9))
               allocate (BBinv(9))
               Do_Elem_iE: do iELoc = 1, size(CellID)
                  PetscCall(DMPlexComputeCellGeometryAffineFEM(dm, cellID(iELoc), v0, BB, BBinv, detBBinv, ierr))
                  !!! Petsc uses a reference simplex with vertices at (-1,-1), (1,-1) and (-1,1)
                  !!! Whereas MEF90 uses (0,0), (1,0), (0,1), so we need to rescale the affine transform
                  Bt%XX = BBinv(1) * 0.5_kr; Bt%XY = BBinv(4) * 0.5_kr; Bt%XZ = BBinv(7) * 0.5_kr
                  Bt%YX = BBinv(2) * 0.5_kr; Bt%YY = BBinv(5) * 0.5_kr; Bt%YZ = BBinv(8) * 0.5_kr
                  Bt%ZX = BBinv(3) * 0.5_kr; Bt%ZY = BBinv(6) * 0.5_kr; Bt%ZZ = BBinv(9) * 0.5_kr
                  detBBinv = 8.0_kr * abs(detBBinv)
                  call ElementPLagrange3DScalInit(dElem(iELoc), Bt, detBBinv, elemType%order, dQuadratureOrder, ierr)
               end do Do_Elem_iE
               deallocate (BBinv)
               deallocate (BB)
               deallocate (v0)
            end if
            PetscCall(ISRestoreIndices(CellIS, CellID, ierr))
         case (MEF90P1Lagrange3DBoundary%name, MEF90P2Lagrange3DBoundary%name)
            PetscCall(ISGetIndices(CellIS, CellID, ierr))
            allocate (dElem(size(cellID)), stat=ierr)
            if (size(CellID) > 0) then
               allocate (centroid(3))
               allocate (innerNormal(3))
               do iELoc = 1, size(CellID)
                  PetscCall(DMPlexComputeCellGeometryFVM(dm, cellID(iEloc), area, centroid, innerNormal, ierr))
                  outerNormal = -innerNormal
                  call ElementPLagrange3DBoundaryScalInit(dElem(iELoc), area, outerNormal, elemType%order, dQuadratureOrder, ierr)
               end do
               deallocate (innerNormal)
               deallocate (centroid)
            end if
            PetscCall(ISRestoreIndices(CellIS, CellID, ierr))
            !Case (MEF90Q1Lagrange3D%name,MEF90Q2Lagrange3D%name,MEF90Q1Lagrange3DBoundary%name,MEF90Q2Lagrange3DBoundary%name)
            !   !!! Get quadrature points for the current element using DMPlexComputeCellGeometryFEM
            !   !!! Initialize element
         case Default
            write (*, *) __FUNCT__, ': Element type not implemented yet', elemType%name
            ierr = PETSC_ERR_SUP
         end select
      end subroutine Element3DScalInitSet

#undef __FUNCT__
#define __FUNCT__ "Element3DVectInitSet"
      subroutine Element3DVectInitSet(dm, cellIS, dElem, dQuadratureOrder, elemType, ierr)
         type(tDM), intent(IN)                             :: dm
         type(tIS), intent(IN)                             :: cellIS
         type(MEF90Element3DVect), dimension(:), pointer    :: dElem
         PetscInt, intent(IN)                              :: dQuadratureOrder
         type(MEF90ElementType), intent(IN)                :: elemType
         PetscErrorCode, intent(OUT)                       :: ierr

         PetscInt                                         :: iELoc
         PetscInt, dimension(:), pointer                    :: CellID
         PetscReal, dimension(:), pointer                 :: v0, BB, BBinv
         PetscReal                                        :: detBBinv
         type(Mat3D)                                      :: Bt
         PetscReal                                        :: area
         PetscReal, dimension(:), pointer                   :: centroid, innerNormal
         type(Vect3D)                                     :: outerNormal

         select case (elemType%name)
         case (MEF90P1Lagrange3D%name, MEF90P2Lagrange3D%name)
            PetscCall(ISGetIndices(CellIS, CellID, ierr))
            allocate (dElem(size(cellID)), stat=ierr)
            if (size(CellID) > 0) then
               allocate (v0(3))
               allocate (BB(9))
               allocate (BBinv(9))
               Do_Elem_iE: do iELoc = 1, size(CellID)
                  PetscCall(DMPlexComputeCellGeometryAffineFEM(dm, cellID(iELoc), v0, BB, BBinv, detBBinv, ierr))
                  !!! Petsc uses a reference simplex with vertices at (-1,-1), (1,-1) and (-1,1)
                  !!! Whereas MEF90 uses (0,0), (1,0), (0,1), so we need to rescale the affine transform
                  Bt%XX = BBinv(1) * 0.5_kr; Bt%XY = BBinv(4) * 0.5_kr; Bt%XZ = BBinv(7) * 0.5_kr
                  Bt%YX = BBinv(2) * 0.5_kr; Bt%YY = BBinv(5) * 0.5_kr; Bt%YZ = BBinv(8) * 0.5_kr
                  Bt%ZX = BBinv(3) * 0.5_kr; Bt%ZY = BBinv(6) * 0.5_kr; Bt%ZZ = BBinv(9) * 0.5_kr
                  detBBinv = 8.0_kr * abs(detBBinv)
                  call ElementPLagrange3DVectInit(dElem(iELoc), Bt, detBBinv, elemType%order, dQuadratureOrder, ierr)
               end do Do_Elem_iE
               deallocate (BBinv)
               deallocate (BB)
               deallocate (v0)
            end if
            PetscCall(ISRestoreIndices(CellIS, CellID, ierr))
         case (MEF90P1Lagrange3DBoundary%name, MEF90P2Lagrange3DBoundary%name)
            PetscCall(ISGetIndices(CellIS, CellID, ierr))
            allocate (dElem(size(cellID)), stat=ierr)
            if (size(CellID) > 0) then
               allocate (centroid(3))
               allocate (innerNormal(3))
               do iELoc = 1, size(CellID)
                  PetscCall(DMPlexComputeCellGeometryFVM(dm, cellID(iEloc), area, centroid, innerNormal, ierr))
                  outerNormal = -innerNormal
                  call ElementPLagrange3DBoundaryVectInit(dElem(iELoc), area, outerNormal, elemType%order, dQuadratureOrder, ierr)
               end do
               deallocate (innerNormal)
               deallocate (centroid)
            end if
            PetscCall(ISRestoreIndices(CellIS, CellID, ierr))
            !Case (MEF90Q1Lagrange3D%name,MEF90Q2Lagrange3D%name,MEF90Q1Lagrange3DBoundary%name,MEF90Q2Lagrange3DBoundary%name)
            !   !!! Get quadrature points for the current element using DMPlexComputeCellGeometryFEM
            !   !!! Initialize element
         case Default
            write (*, *) __FUNCT__, ': Element type not implemented yet', elemType%name
            ierr = PETSC_ERR_SUP
         end select
      end subroutine Element3DVectInitSet

#undef __FUNCT__
#define __FUNCT__ "ElementPLagrange2DScalInit"
      subroutine ElementPLagrange2DScalInit(dElem, Bt, DetBinv, dPolynomialOrder, dQuadratureOrder, ierr)
         ! Compute the quadrature weights and the value of the basis functions and their gradient
         ! at the quadrature points.
         ! Quadrature rules courtesy of Shawn Walker walker@math.lsu.edu
         type(MEF90Element2DScal), intent(INOUT) :: dElem
         type(Mat2D), intent(IN)                 :: Bt          ! The transposed of transformation matrix
         PetscReal, intent(IN)                   :: DetBinv     ! The determinant of B^{-1}
         PetscInt                               :: dPolynomialOrder, dQuadratureOrder
         PetscErrorCode, intent(OUT)             :: ierr

         PetscInt                               :: Nb_Gauss
         PetscInt                               :: Num_Dof
         PetscInt                               :: iDoF, iG

         PetscReal, dimension(:, :), pointer       :: PhiHat      ! PhiHat(i,k) The value of the ith basis function at the kth integration point
         type(Vect2D), dimension(:, :), pointer    :: GradPhiHat

         type(Vect2D), dimension(:), pointer     :: Xi ! The quadrature points coordinates in the reference element

         select case (dQuadratureOrder)
         case (0, 1)
            Nb_Gauss = 1
            allocate (Xi(Nb_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
            Xi(1)%X = 1.0_kr / 3.0_kr
            Xi(1)%Y = 1.0_kr / 3.0_kr
            dElem%Gauss_C = detBinv / 2.0_kr
         case (2)
            Nb_Gauss = 3
            allocate (Xi(Nb_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
            Xi(1) = [1.0_kr / 6.0_kr, 1.0_kr / 6.0_kr]
            Xi(2) = [2.0_kr / 3.0_kr, 1.0_kr / 6.0_kr]
            Xi(3) = [1.0_kr / 6.0_kr, 2.0_kr / 3.0_kr]
            dElem%Gauss_C(1:Nb_Gauss) = detBinv / 6.0_kr
            dElem%Gauss_C(Nb_Gauss) = detBinv / 2.0_kr - sum(dElem%Gauss_C(1:Nb_Gauss - 1))
         case (3)
            Nb_Gauss = 4
            allocate (Xi(Nb_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
            Xi(1) = [1.0_kr / 3.0_kr, 1.0_kr / 3.0_kr]
            Xi(2) = [3.0_kr / 5.0_kr, 1.0_kr / 5.0_kr]
            Xi(3) = [1.0_kr / 5.0_kr, 3.0_kr / 5.0_kr]
            Xi(4) = [1.0_kr / 5.0_kr, 1.0_kr / 5.0_kr]
            dElem%Gauss_C(1) = -detBinv * 9.0_kr / 32.0_kr
            dElem%Gauss_C(2:4) = detBinv * 25.0_kr / 96.0_kr
            dElem%Gauss_C(Nb_Gauss) = detBinv / 2.0_kr - sum(dElem%Gauss_C(1:Nb_Gauss - 1))
         case (4)
            Nb_Gauss = 6
            allocate (Xi(Nb_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
            Xi(1) = [0.816847572980459_kr, 0.091576213509771_kr]
            Xi(2) = [0.091576213509771_kr, 0.816847572980459_kr]
            Xi(3) = [0.091576213509771_kr, 0.091576213509771_kr]
            Xi(4) = [0.108103018168070_kr, 0.445948490915965_kr]
            Xi(5) = [0.445948490915965_kr, 0.108103018168070_kr]
            Xi(6) = [0.445948490915965_kr, 0.445948490915965_kr]
            dElem%Gauss_C(1:3) = 0.109951743655322 / 2.0_kr * detBinv
            dElem%Gauss_C(4:6) = 0.223381589678011 / 2.0_kr * detBinv
            dElem%Gauss_C(Nb_Gauss) = detBinv / 2.0_kr - sum(dElem%Gauss_C(1:Nb_Gauss - 1))
         case (5)
            Nb_Gauss = 7
            allocate (Xi(Nb_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
            Xi(1) = [0.3333333333333333_kr, 0.3333333333333333_kr]
            Xi(2) = [0.4701420641051150_kr, 0.4701420641051150_kr]
            Xi(3) = [0.0597158717897698_kr, 0.4701420641051150_kr]
            Xi(4) = [0.4701420641051150_kr, 0.0597158717897698_kr]
            Xi(5) = [0.1012865073234563_kr, 0.1012865073234563_kr]
            Xi(6) = [0.7974269853530872_kr, 0.1012865073234563_kr]
            Xi(7) = [0.1012865073234563_kr, 0.7974269853530872_kr]
            dElem%Gauss_C(1) = .225_kr / 2.0_kr * detBinv
            dElem%Gauss_C(2:4) = 0.13239415278850616 / 2.0_kr * detBinv
            dElem%Gauss_C(5:7) = 0.12593918054482717 / 2.0_kr * detBinv
            dElem%Gauss_C(Nb_Gauss) = detBinv / 2.0_kr - sum(dElem%Gauss_C(1:Nb_Gauss - 1))
         case (-6)
         !!! It seems to me that this is really a quadrature rule of order 5...
            Nb_Gauss = 9
            allocate (Xi(Nb_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
            Xi(1) = [0.124949503233232_kr, 0.437525248383384_kr]
            Xi(2) = [0.437525248383384_kr, 0.124949503233232_kr]
            Xi(3) = [0.437525248383384_kr, 0.437525248383384_kr]
            Xi(4) = [0.797112651860071_kr, 0.165409927389841_kr]
            Xi(5) = [0.797112651860071_kr, 0.037477420750088_kr]
            Xi(6) = [0.165409927389841_kr, 0.797112651860071_kr]
            Xi(7) = [0.165409927389841_kr, 0.037477420750088_kr]
            Xi(8) = [0.037477420750088_kr, 0.797112651860071_kr]
            Xi(9) = [0.037477420750088_kr, 0.165409927389841_kr]
            dElem%Gauss_C(1:3) = 0.205950504760887_kr / 2.0_kr * detBinv
            dElem%Gauss_C(4:9) = 0.063691414286223_kr / 2.0_kr * detBinv
            dElem%Gauss_C(Nb_Gauss) = detBinv / 2.0_kr - sum(dElem%Gauss_C(1:Nb_Gauss - 1))
         case (6)
            Nb_Gauss = 12
            allocate (Xi(Nb_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
            Xi(1) = [0.873821971016996_kr, 0.063089014491502_kr]
            Xi(2) = [0.063089014491502_kr, 0.873821971016996_kr]
            Xi(3) = [0.063089014491502_kr, 0.063089014491502_kr]
            Xi(4) = [0.501426509658179_kr, 0.249286745170910_kr]
            Xi(5) = [0.249286745170910_kr, 0.501426509658179_kr]
            Xi(6) = [0.249286745170910_kr, 0.249286745170910_kr]
            Xi(7) = [0.636502499121399_kr, 0.310352451033785_kr]
            Xi(8) = [0.636502499121399_kr, 0.053145049844816_kr]
            Xi(9) = [0.310352451033785_kr, 0.636502499121399_kr]
            Xi(10) = [0.310352451033785_kr, 0.053145049844816_kr]
            Xi(11) = [0.053145049844816_kr, 0.636502499121399_kr]
            Xi(12) = [0.053145049844816_kr, 0.310352451033785_kr]
            dElem%Gauss_C(1:3) = 0.050844906370207_kr / 2.0_kr * detBinv
            dElem%Gauss_C(4:6) = 0.116786275726379_kr / 2.0_kr * detBinv
            dElem%Gauss_C(7:12) = 0.082851075618374_kr / 2.0_kr * detBinv
            dElem%Gauss_C(Nb_Gauss) = detBinv / 2.0_kr - sum(dElem%Gauss_C(1:Nb_Gauss - 1))
         case (7)
            Nb_Gauss = 13
            allocate (Xi(Nb_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
            Xi(1) = [0.333333333333333_kr, 0.333333333333333_kr]
            Xi(2) = [0.479308067841923_kr, 0.260345966079038_kr]
            Xi(3) = [0.260345966079038_kr, 0.479308067841923_kr]
            Xi(4) = [0.260345966079038_kr, 0.260345966079038_kr]
            Xi(5) = [0.869739794195568_kr, 0.065130102902216_kr]
            Xi(6) = [0.065130102902216_kr, 0.869739794195568_kr]
            Xi(7) = [0.065130102902216_kr, 0.065130102902216_kr]
            Xi(8) = [0.638444188569809_kr, 0.312865496004875_kr]
            Xi(9) = [0.638444188569809_kr, 0.048690315425316_kr]
            Xi(10) = [0.312865496004875_kr, 0.638444188569809_kr]
            Xi(11) = [0.312865496004875_kr, 0.048690315425316_kr]
            Xi(12) = [0.048690315425316_kr, 0.638444188569809_kr]
            Xi(13) = [0.048690315425316_kr, 0.312865496004875_kr]
            dElem%Gauss_C(1) = -0.149570044467670_kr / 2.0_kr * detBinv
            dElem%Gauss_C(2:4) = 0.175615257433204_kr / 2.0_kr * detBinv
            dElem%Gauss_C(5:7) = 0.053347235608839_kr / 2.0_kr * detBinv
            dElem%Gauss_C(8:13) = 0.077113760890257_kr / 2.0_kr * detBinv
            dElem%Gauss_C(Nb_Gauss) = detBinv / 2.0_kr - sum(dElem%Gauss_C(1:Nb_Gauss - 1))
         case (8, 9)
            Nb_Gauss = 19
            allocate (Xi(Nb_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
            Xi(1) = [0.3333333333333333_kr, 0.3333333333333333_kr]
            Xi(2) = [0.0206349616025259_kr, 0.4896825191987370_kr]
            Xi(3) = [0.4896825191987370_kr, 0.0206349616025259_kr]
            Xi(4) = [0.4896825191987370_kr, 0.4896825191987370_kr]
            Xi(5) = [0.1258208170141290_kr, 0.4370895914929355_kr]
            Xi(6) = [0.4370895914929355_kr, 0.1258208170141290_kr]
            Xi(7) = [0.4370895914929355_kr, 0.4370895914929355_kr]
            Xi(8) = [0.6235929287619356_kr, 0.1882035356190322_kr]
            Xi(9) = [0.1882035356190322_kr, 0.6235929287619356_kr]
            Xi(10) = [0.1882035356190322_kr, 0.1882035356190322_kr]
            Xi(11) = [0.9105409732110941_kr, 0.0447295133944530_kr]
            Xi(12) = [0.0447295133944530_kr, 0.9105409732110941_kr]
            Xi(13) = [0.0447295133944530_kr, 0.0447295133944530_kr]
            Xi(14) = [0.7411985987844980_kr, 0.0368384120547363_kr]
            Xi(15) = [0.7411985987844980_kr, 0.2219628891607657_kr]
            Xi(16) = [0.0368384120547363_kr, 0.7411985987844980_kr]
            Xi(17) = [0.0368384120547363_kr, 0.2219628891607657_kr]
            Xi(18) = [0.2219628891607657_kr, 0.7411985987844980_kr]
            Xi(19) = [0.2219628891607657_kr, 0.0368384120547363_kr]
            dElem%Gauss_C(1) = 0.09713579628279610_kr / 2.0_kr * detBinv
            dElem%Gauss_C(2:4) = 0.03133470022713983_kr / 2.0_kr * detBinv
            dElem%Gauss_C(5:7) = 0.07782754100477543_kr / 2.0_kr * detBinv
            dElem%Gauss_C(8:10) = 0.07964773892720910_kr / 2.0_kr * detBinv
            dElem%Gauss_C(11:13) = 0.02557767565869810_kr / 2.0_kr * detBinv
            dElem%Gauss_C(14:19) = 0.04328353937728940_kr / 2.0_kr * detBinv
            dElem%Gauss_C(Nb_Gauss) = detBinv / 2.0_kr - sum(dElem%Gauss_C(1:Nb_Gauss - 1))
         case Default
            write (*, *) __FUNCT__, ': Unimplemented quadrature order', dQuadratureOrder
            stop
         end select

         select case (dPolynomialOrder)
         case (1)
            Num_DoF = 3
            allocate (PhiHat(Num_DoF, Nb_Gauss), stat=ierr)
            allocate (GradPhiHat(Num_DoF, Nb_Gauss), stat=ierr)
            PhiHat(1, :) = 1.0_kr - Xi%X - Xi%Y
            PhiHat(2, :) = Xi%X
            PhiHat(3, :) = Xi%Y

            GradPhiHat(1, :)%X = -1.0_kr; GradPhiHat(1, :)%Y = -1.0_kr
            GradPhiHat(2, :)%X = 1.0_kr; GradPhiHat(2, :)%Y = 0.0_kr
            GradPhiHat(3, :)%X = 0.0_kr; GradPhiHat(3, :)%Y = 1.0_kr

         case (2)
            Num_DoF = 6
            allocate (PhiHat(Num_DoF, Nb_Gauss), stat=ierr)
            allocate (GradPhiHat(Num_DoF, Nb_Gauss), stat=ierr)
         !!! See Dof local Ordering.md for the unusual dof ordering
            PhiHat(4, :) = (1.0_kr - Xi%X - Xi%Y) * (1.0_kr - 2.0_kr * Xi%X - 2.0_kr * Xi%Y)
            PhiHat(5, :) = Xi%X * (2.0_kr * Xi%X - 1.0_kr)
            PhiHat(6, :) = Xi%Y * (2.0_kr * Xi%Y - 1.0_kr)
            PhiHat(1, :) = 4.0_kr * Xi%X * (1.0_kr - Xi%X - Xi%Y)
            PhiHat(2, :) = 4.0_kr * Xi%X * Xi%Y
            PhiHat(3, :) = 4.0_kr * Xi%Y * (1.0_kr - Xi%X - Xi%Y)

            GradPhiHat(4, :)%X = 4.0_kr * Xi%X + 4.0_kr * Xi%y - 3.0_kr; GradPhiHat(4, :)%Y = 4.0_kr * Xi%X + 4.0_kr * Xi%y - 3.0_kr
            GradPhiHat(5, :)%X = 4.0_kr * Xi%X - 1.0_kr; GradPhiHat(5, :)%Y = 0.0_kr
            GradPhiHat(6, :)%X = 0.0_kr; GradPhiHat(6, :)%Y = 4.0_kr * Xi%Y - 1.0_kr
            GradPhiHat(1, :)%X = 4.0_kr * (1.0_kr - 2.0_kr * Xi%X - Xi%Y); GradPhiHat(1, :)%Y = -4.0_kr * Xi%X
            GradPhiHat(2, :)%X = 4.0_kr * Xi%Y; GradPhiHat(2, :)%Y = 4.0_kr * Xi%X
            GradPhiHat(3, :)%X = -4.0_kr * Xi%Y; GradPhiHat(3, :)%Y = 4.0_kr * (1.0_kr - Xi%X - 2.0_kr * Xi%Y)
         case Default
            Num_DoF = 0
            write (*, *) __FUNCT__, ': Unimplemented PolynomialOrder', dPolynomialOrder
            stop
         end select

         allocate (dElem%BF(Num_DoF, Nb_Gauss), stat=ierr)
         allocate (dElem%Grad_BF(Num_DoF, Nb_Gauss), stat=ierr)
         dElem%BF = PhiHat
         do iDoF = 1, Num_DoF
            do iG = 1, Nb_Gauss
               dElem%Grad_BF(iDoF, iG) = Bt * GradPhiHat(iDoF, iG)
            end do
         end do

         dElem%outerNormal = 0.0_kr

         deallocate (Xi, stat=ierr)
         deallocate (PhiHat, stat=ierr)
         deallocate (GradPhiHat, stat=ierr)
      end subroutine ElementPLagrange2DScalInit

#undef __FUNCT__
#define __FUNCT__ "ElementPLagrange2DBoundaryScalInit"
      subroutine ElementPLagrange2DBoundaryScalInit(dElem, l, outerNormal, dPolynomialOrder, dQuadratureOrder, ierr)
         type(MEF90Element2DScal), intent(INOUT) :: dElem
         PetscReal, intent(IN)                   :: l
         type(Vect2D), intent(IN)                :: outerNormal
         PetscInt, intent(IN)                    :: dPolynomialOrder, dQuadratureOrder
         PetscErrorCode, intent(OUT)             :: ierr

         PetscReal, dimension(:), pointer       :: Xi
         PetscInt                               :: iDoF, iG, Num_Gauss, Num_DoF

         num_Dof = 0
         num_Gauss = 0
         dElem%outerNormal = outerNormal
         select case (dQuadratureOrder)
         case (0, 1)
            Num_Gauss = 1
            allocate (Xi(Num_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Num_Gauss), stat=ierr)
            Xi(1) = 0.0_kr
            dElem%Gauss_C(1) = l
         case (2, 3)
            Num_Gauss = 2
            allocate (Xi(Num_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Num_Gauss), stat=ierr)
            Xi(1) = -1.0_kr / sqrt(3.0_kr)
            Xi(2) = 1.0_kr / sqrt(3.0_kr)
            dElem%Gauss_C = l*.5_kr
         case (4, 5)
            Num_Gauss = 3
            allocate (Xi(Num_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Num_Gauss), stat=ierr)
            Xi(1) = -sqrt(15.0_kr) / 5.0_kr
            Xi(2) = 0.0_kr
            Xi(3) = sqrt(15.0_kr) / 5.0_kr
            dElem%Gauss_C(1) = l * 5.0_kr / 18.0_kr
            dElem%Gauss_C(2) = l * 4.0_kr / 9.0_kr
            dElem%Gauss_C(3) = l * 5.0_kr / 18.0_kr
         case (6, 7)
            Num_Gauss = 4
            allocate (Xi(Num_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Num_Gauss), stat=ierr)
            Xi(1) = -sqrt((15.0_kr + 2.0_kr * sqrt(30.0_kr)) / 35.0_kr)
            Xi(2) = -sqrt((15.0_kr - 2.0_kr * sqrt(30.0_kr)) / 35.0_kr)
            Xi(3) = sqrt((15.0_kr - 2.0_kr * sqrt(30.0_kr)) / 35.0_kr)
            Xi(4) = sqrt((15.0_kr + 2.0_kr * sqrt(30.0_kr)) / 35.0_kr)
            dElem%Gauss_C(1) = l * (.25_kr - sqrt(5.0_kr / 864.0_kr))
            dElem%Gauss_C(2) = l * (.25_kr + sqrt(5.0_kr / 864.0_kr))
            dElem%Gauss_C(3) = l * (.25_kr + sqrt(5.0_kr / 864.0_kr))
            dElem%Gauss_C(4) = l * (.25_kr - sqrt(5.0_kr / 864.0_kr))
         case (8, 9, 10, 11, 12)
            Num_Gauss = 5
            allocate (Xi(Num_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Num_Gauss), stat=ierr)
            Xi(1) = -sqrt(5.0_kr + sqrt(40.0_kr / 7.0_kr)) / 3.0_kr
            Xi(2) = -sqrt(5.0_kr - sqrt(40.0_kr / 7.0_kr)) / 3.0_kr
            Xi(3) = 0.0_kr
            Xi(4) = sqrt(5.0_kr - sqrt(40.0_kr / 7.0_kr)) / 3.0_kr
            Xi(5) = sqrt(5.0_kr + sqrt(40.0_kr / 7.0_kr)) / 3.0_kr
            dElem%Gauss_C(1) = l * (322.0_kr - sqrt(11830.0_kr)) / 1800.0_kr
            dElem%Gauss_C(2) = l * (322.0_kr + sqrt(11830.0_kr)) / 1800.0_kr
            dElem%Gauss_C(3) = l * 128.0_kr / 450.0_kr
            dElem%Gauss_C(4) = l * (322.0_kr + sqrt(11830.0_kr)) / 1800.0_kr
            dElem%Gauss_C(5) = l * (322.0_kr - sqrt(11830.0_kr)) / 1800.0_kr
         case Default
            write (*, *) __FUNCT__, ': Unimplemented quadrature order', dQuadratureOrder
            stop
         end select
         select case (dPolynomialOrder)
         case (1)
            Num_DoF = 2
            allocate (dElem%BF(Num_DoF, Num_Gauss), stat=ierr)
            dElem%BF(1, :) = (1.0_kr - Xi)*.5_kr
            dElem%BF(2, :) = (1.0_kr + Xi)*.5_kr
         case (2)
            Num_DoF = 3
            allocate (dElem%BF(Num_DoF, Num_Gauss), stat=ierr)
         !!! See Dof local Ordering.md for the unusual dof ordering
            dElem%BF(2, :) = Xi * (Xi - 1.0_kr)*.5_kr
            dElem%BF(3, :) = Xi * (Xi + 1.0_kr)*.5_kr
            dElem%BF(1, :) = (1.0_kr - Xi) * (1.0_kr + Xi)
         case Default
            write (*, *) '[ERROR]: Polynomial order ', dPolynomialOrder, ' not implemented in ', __FUNCT__
         end select
         allocate (delem%Grad_BF(Num_DoF, Num_Gauss), stat=ierr)
         do iDof = 1, num_dof
            do iG = 1, num_Gauss
               delem%Grad_BF(iDof, iG) = 0.0_kr
            end do
         end do
         deallocate (Xi, stat=ierr)
      end subroutine ElementPLagrange2DBoundaryScalInit

#undef __FUNCT__
#define __FUNCT__ "ElementPLagrange2DVectInit"
      subroutine ElementPLagrange2DVectInit(dElem, Bt, DetBinv, dPolynomialOrder, dQuadratureOrder, ierr)
         type(MEF90Element2DVect), intent(INOUT) :: dElem
         type(Mat2D), intent(IN)                 :: Bt          ! The transposed of transformation matrix
         PetscReal, intent(IN)                   :: DetBinv     ! The determinant of B^{-1}
         PetscInt, intent(IN)                    :: dPolynomialOrder, dQuadratureOrder
         PetscErrorCode, intent(OUT)             :: ierr

         type(MEF90Element2DScal)              :: Elem_Scal
         PetscInt                               :: dim = 2
         PetscInt                               :: Num_DoF, Nb_Gauss, i, iDof, iG

         call ElementPLagrange2DScalInit(Elem_Scal, Bt, DetBinv, dPolynomialOrder, dQuadratureOrder, ierr)
         Num_DoF = size(Elem_Scal%BF, 1)
         Nb_Gauss = size(Elem_Scal%BF, 2)

         allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
         allocate (dElem%BF(Num_DoF * dim, Nb_Gauss), stat=ierr)
         allocate (dElem%Grad_BF(Num_DoF * dim, Nb_Gauss), stat=ierr)
         allocate (dElem%GradS_BF(Num_DoF * dim, Nb_Gauss), stat=ierr)

         dElem%Gauss_C = Elem_Scal%Gauss_C
         do iDof = 1, num_dof * dim
            do iG = 1, Nb_Gauss
               dElem%BF(iDof, iG) = 0.0_kr
               delem%Grad_BF(iDof, iG) = 0.0_kr
               delem%GradS_BF(iDof, iG) = 0.0_kr
            end do
         end do

         do i = 0, Num_DoF - 1
            dElem%BF(i * dim + 1, :)%X = Elem_Scal%BF(i + 1, :)
            dElem%BF(i * dim + 2, :)%Y = Elem_Scal%BF(i + 1, :)

            dElem%Grad_BF(i * dim + 1, :)%XX = Elem_Scal%Grad_BF(i + 1, :)%X
            dElem%Grad_BF(i * dim + 1, :)%XY = Elem_Scal%Grad_BF(i + 1, :)%Y
            dElem%Grad_BF(i * dim + 2, :)%YX = Elem_Scal%Grad_BF(i + 1, :)%X
            dElem%Grad_BF(i * dim + 2, :)%YY = Elem_Scal%Grad_BF(i + 1, :)%Y

            dElem%GradS_BF(i * dim + 1, :)%XX = Elem_Scal%Grad_BF(i + 1, :)%X
            dElem%GradS_BF(i * dim + 1, :)%XY = Elem_Scal%Grad_BF(i + 1, :)%Y / 2.0_kr
            dElem%GradS_BF(i * dim + 2, :)%XY = Elem_Scal%Grad_BF(i + 1, :)%X / 2.0_kr
            dElem%GradS_BF(i * dim + 2, :)%YY = Elem_Scal%Grad_BF(i + 1, :)%Y
         end do

         call MEF90ElementDestroy(Elem_Scal, ierr)
      end subroutine ElementPLagrange2DVectInit

#undef __FUNCT__
#define __FUNCT__ "ElementPLagrange2DBoundaryVectInit"
      subroutine ElementPLagrange2DBoundaryVectInit(dElem, length, outerNormal, dPolynomialOrder, dQuadratureOrder, ierr)
         type(MEF90Element2DVect), intent(INOUT) :: dElem
         PetscReal, intent(IN)                   :: length
         type(Vect2D), intent(IN)                :: outerNormal
         PetscInt, intent(IN)                    :: dPolynomialOrder, dQuadratureOrder
         PetscErrorCode, intent(OUT)             :: ierr

         type(MEF90Element2DScal)              :: Elem_Scal
         PetscInt                               :: dim = 2
         PetscInt                               :: Num_DoF, Nb_Gauss, iDof, iG

         call ElementPLagrange2DBoundaryScalInit(Elem_Scal, length, outerNormal, dPolynomialOrder, dQuadratureOrder, ierr)
         Num_DoF = size(Elem_Scal%BF, 1)
         Nb_Gauss = size(Elem_Scal%BF, 2)
         allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
         allocate (dElem%BF(Num_DoF * dim, Nb_Gauss), stat=ierr)
         allocate (dElem%Grad_BF(Num_DoF * dim, Nb_Gauss), stat=ierr)
         allocate (dElem%GradS_BF(Num_DoF * dim, Nb_Gauss), stat=ierr)

         dElem%Gauss_C = Elem_Scal%Gauss_C
         do iDof = 1, num_dof * dim
            do iG = 1, Nb_Gauss
               dElem%BF(iDof, iG) = 0.0_kr
               delem%Grad_BF(iDof, iG) = 0.0_kr
               delem%GradS_BF(iDof, iG) = 0.0_kr
            end do
         end do

         do idof = 0, Num_DoF - 1
            dElem%BF(iDof * dim + 1, :)%X = Elem_Scal%BF(iDof + 1, :)
            dElem%BF(iDof * dim + 2, :)%Y = Elem_Scal%BF(iDof + 1, :)
         end do

         dElem%outerNormal = Elem_Scal%outerNormal
         call MEF90ElementDestroy(Elem_Scal, ierr)
      end subroutine ElementPLagrange2DBoundaryVectInit

#undef __FUNCT__
#define __FUNCT__ "ElementPLagrange3DScalInit"
      subroutine ElementPLagrange3DScalInit(dElem, Bt, DetBinv, dPolynomialOrder, dQuadratureOrder, ierr)
         ! Compute the quadrature weights and the value of the basis functions and their gradient
         ! at the quadrature points.
         type(MEF90Element3DScal), intent(INOUT) :: dElem
         type(Mat3D), intent(IN)                 :: Bt          ! The transposed of transformation matrix
         PetscReal, intent(IN)                   :: DetBinv     ! The determinant of B^{-1}
         PetscInt, intent(IN)                    :: dPolynomialOrder, dQuadratureOrder
         PetscErrorCode, intent(OUT)             :: ierr

         PetscInt                               :: Nb_Gauss
         PetscInt                               :: Num_Dof
         PetscInt                               :: iDoF, iG

         PetscReal, dimension(:, :), pointer       :: PhiHat      ! PhiHat(i,k) The value of the ith basis function at the kth integration point
         type(Vect3D), dimension(:, :), pointer    :: GradPhiHat

         type(Vect3D), dimension(:), pointer      :: Xi          ! The quadrature points coordinates in the reference element
         PetscReal                              :: a, b         ! Location of integration points in Aiken p. 272 table 10.4

         select case (dQuadratureOrder)
         case (0, 1)
            Nb_Gauss = 1
            allocate (Xi(Nb_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
            Xi(1) = [.25_kr, .25_kr, .25_kr]
            dElem%Gauss_C(1) = 1.0_kr / 6.0_kr * detBinv
         case (2)
            a = (5.0_kr + 3.0_kr * sqrt(5.0_kr)) / 20.0_kr
            b = (5.0_kr - sqrt(5.0_kr)) / 20.0_kr
            Nb_Gauss = 4
            allocate (Xi(Nb_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
            Xi(1) = [a, b, b]
            Xi(2) = [b, a, b]
            Xi(3) = [b, b, a]
            Xi(4) = [b, b, b]
            dElem%Gauss_C(1:4) = 1.0_kr / 24.0_kr * detBinv
         case (3)
            Nb_Gauss = 5
            allocate (Xi(Nb_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
            Xi(1) = [.25_kr, .25_kr, .25_kr]
            Xi(2) = [.5_kr, 1._kr / 6._kr, 1._kr / 6._kr]
            Xi(3) = [1._kr / 6._kr, .5_kr, 1._kr / 6._kr]
            Xi(4) = [1._kr / 6._kr, 1._kr / 6._kr, .5_kr]
            Xi(5) = [1._kr / 6._kr, 1._kr / 6._kr, 1._kr / 6._kr]
            dElem%Gauss_C(1) = -2.0_kr / 15.0_kr * detBinv
            dElem%Gauss_C(2:5) = 3.0_kr / 40.0_kr * detBinv
         case (4)
            Nb_Gauss = 11
            allocate (Xi(Nb_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
            Xi(1) = [0.2500000000000000_kr, 0.2500000000000000_kr, 0.2500000000000000_kr]
            Xi(2) = [0.7857142857142857_kr, 0.0714285714285714_kr, 0.0714285714285714_kr]
            Xi(3) = [0.0714285714285714_kr, 0.0714285714285714_kr, 0.0714285714285714_kr]
            Xi(4) = [0.0714285714285714_kr, 0.0714285714285714_kr, 0.7857142857142857_kr]
            Xi(5) = [0.0714285714285714_kr, 0.7857142857142857_kr, 0.0714285714285714_kr]
            Xi(6) = [0.1005964238332008_kr, 0.3994035761667992_kr, 0.3994035761667992_kr]
            Xi(7) = [0.3994035761667992_kr, 0.1005964238332008_kr, 0.3994035761667992_kr]
            Xi(8) = [0.3994035761667992_kr, 0.3994035761667992_kr, 0.1005964238332008_kr]
            Xi(9) = [0.3994035761667992_kr, 0.1005964238332008_kr, 0.1005964238332008_kr]
            Xi(10) = [0.1005964238332008_kr, 0.3994035761667992_kr, 0.1005964238332008_kr]
            Xi(11) = [0.1005964238332008_kr, 0.1005964238332008_kr, 0.3994035761667992_kr]
            dElem%Gauss_C(1) = -0.0789333333333333_kr / 6.0_kr * detBinv
            dElem%Gauss_C(2:5) = 0.0457333333333333_kr / 6.0_kr * detBinv
            dElem%Gauss_C(6:11) = 0.1493333333333333_kr / 6.0_kr * detBinv
         case (5)
            Nb_Gauss = 15
            allocate (Xi(Nb_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
            Xi(1) = [0.2500000000000000_kr, 0.2500000000000000_kr, 0.2500000000000000_kr]
            Xi(2) = [0.0000000000000000_kr, 0.3333333333333333_kr, 0.3333333333333333_kr]
            Xi(3) = [0.3333333333333333_kr, 0.3333333333333333_kr, 0.3333333333333333_kr]
            Xi(4) = [0.3333333333333333_kr, 0.3333333333333333_kr, 0.0000000000000000_kr]
            Xi(5) = [0.3333333333333333_kr, 0.0000000000000000_kr, 0.3333333333333333_kr]
            Xi(6) = [0.7272727272727273_kr, 0.0909090909090909_kr, 0.0909090909090909_kr]
            Xi(7) = [0.0909090909090909_kr, 0.0909090909090909_kr, 0.0909090909090909_kr]
            Xi(8) = [0.0909090909090909_kr, 0.0909090909090909_kr, 0.7272727272727273_kr]
            Xi(9) = [0.0909090909090909_kr, 0.7272727272727273_kr, 0.0909090909090909_kr]
            Xi(10) = [0.4334498464263357_kr, 0.0665501535736643_kr, 0.0665501535736643_kr]
            Xi(11) = [0.0665501535736643_kr, 0.4334498464263357_kr, 0.0665501535736643_kr]
            Xi(12) = [0.0665501535736643_kr, 0.0665501535736643_kr, 0.4334498464263357_kr]
            Xi(13) = [0.0665501535736643_kr, 0.4334498464263357_kr, 0.4334498464263357_kr]
            Xi(14) = [0.4334498464263357_kr, 0.0665501535736643_kr, 0.4334498464263357_kr]
            Xi(15) = [0.4334498464263357_kr, 0.4334498464263357_kr, 0.0665501535736643_kr]
            dElem%Gauss_C(1) = 0.1817020685825351_kr / 6.0_kr * detBinv
            dElem%Gauss_C(2:5) = 0.0361607142857143_kr / 6.0_kr * detBinv
            dElem%Gauss_C(6:9) = 0.0698714945161738_kr / 6.0_kr * detBinv
            dElem%Gauss_C(10:15) = 0.0656948493683187_kr / 6.0_kr * detBinv
         case (6)
            Nb_Gauss = 24
            allocate (Xi(Nb_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
            Xi(1) = [0.3561913862225449_kr, 0.2146028712591517_kr, 0.2146028712591517_kr]
            Xi(2) = [0.2146028712591517_kr, 0.2146028712591517_kr, 0.2146028712591517_kr]
            Xi(3) = [0.2146028712591517_kr, 0.2146028712591517_kr, 0.3561913862225449_kr]
            Xi(4) = [0.2146028712591517_kr, 0.3561913862225449_kr, 0.2146028712591517_kr]
            Xi(5) = [0.8779781243961660_kr, 0.0406739585346113_kr, 0.0406739585346113_kr]
            Xi(6) = [0.0406739585346113_kr, 0.0406739585346113_kr, 0.0406739585346113_kr]
            Xi(7) = [0.0406739585346113_kr, 0.0406739585346113_kr, 0.8779781243961660_kr]
            Xi(8) = [0.0406739585346113_kr, 0.8779781243961660_kr, 0.0406739585346113_kr]
            Xi(9) = [0.0329863295731731_kr, 0.3223378901422757_kr, 0.3223378901422757_kr]
            Xi(10) = [0.3223378901422757_kr, 0.3223378901422757_kr, 0.3223378901422757_kr]
            Xi(11) = [0.3223378901422757_kr, 0.3223378901422757_kr, 0.0329863295731731_kr]
            Xi(12) = [0.3223378901422757_kr, 0.0329863295731731_kr, 0.3223378901422757_kr]
            Xi(13) = [0.2696723314583159_kr, 0.0636610018750175_kr, 0.0636610018750175_kr]
            Xi(14) = [0.0636610018750175_kr, 0.2696723314583159_kr, 0.0636610018750175_kr]
            Xi(15) = [0.0636610018750175_kr, 0.0636610018750175_kr, 0.2696723314583159_kr]
            Xi(16) = [0.6030056647916491_kr, 0.0636610018750175_kr, 0.0636610018750175_kr]
            Xi(17) = [0.0636610018750175_kr, 0.6030056647916491_kr, 0.0636610018750175_kr]
            Xi(18) = [0.0636610018750175_kr, 0.0636610018750175_kr, 0.6030056647916491_kr]
            Xi(19) = [0.0636610018750175_kr, 0.2696723314583159_kr, 0.6030056647916491_kr]
            Xi(20) = [0.2696723314583159_kr, 0.6030056647916491_kr, 0.0636610018750175_kr]
            Xi(21) = [0.6030056647916491_kr, 0.0636610018750175_kr, 0.2696723314583159_kr]
            Xi(22) = [0.0636610018750175_kr, 0.6030056647916491_kr, 0.2696723314583159_kr]
            Xi(23) = [0.2696723314583159_kr, 0.0636610018750175_kr, 0.6030056647916491_kr]
            Xi(24) = [0.6030056647916491_kr, 0.2696723314583159_kr, 0.0636610018750175_kr]
            dElem%Gauss_C(1:4) = 0.0399227502581679_kr / 6.0_kr * detBinv
            dElem%Gauss_C(5:8) = 0.0100772110553207_kr / 6.0_kr * detBinv
            dElem%Gauss_C(9:12) = 0.0553571815436544_kr / 6.0_kr * detBinv
            dElem%Gauss_C(13:24) = 0.0482142857142857_kr / 6.0_kr * detBinv
         case (7)
            Nb_Gauss = 31
            allocate (Xi(Nb_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
            Xi(1) = [0.2500000000000000_kr, 0.2500000000000000_kr, 0.2500000000000000_kr]
            Xi(2) = [0.7653604230090441_kr, 0.0782131923303186_kr, 0.0782131923303186_kr]
            Xi(3) = [0.0782131923303186_kr, 0.0782131923303186_kr, 0.0782131923303186_kr]
            Xi(4) = [0.0782131923303186_kr, 0.0782131923303186_kr, 0.7653604230090441_kr]
            Xi(5) = [0.0782131923303186_kr, 0.7653604230090441_kr, 0.0782131923303186_kr]
            Xi(6) = [0.6344703500082868_kr, 0.1218432166639044_kr, 0.1218432166639044_kr]
            Xi(7) = [0.1218432166639044_kr, 0.1218432166639044_kr, 0.1218432166639044_kr]
            Xi(8) = [0.1218432166639044_kr, 0.1218432166639044_kr, 0.6344703500082868_kr]
            Xi(9) = [0.1218432166639044_kr, 0.6344703500082868_kr, 0.1218432166639044_kr]
            Xi(10) = [0.0023825066607383_kr, 0.3325391644464206_kr, 0.3325391644464206_kr]
            Xi(11) = [0.3325391644464206_kr, 0.3325391644464206_kr, 0.3325391644464206_kr]
            Xi(12) = [0.3325391644464206_kr, 0.3325391644464206_kr, 0.0023825066607383_kr]
            Xi(13) = [0.3325391644464206_kr, 0.0023825066607383_kr, 0.3325391644464206_kr]
            Xi(14) = [0.0000000000000000_kr, 0.5000000000000000_kr, 0.5000000000000000_kr]
            Xi(15) = [0.5000000000000000_kr, 0.0000000000000000_kr, 0.5000000000000000_kr]
            Xi(16) = [0.5000000000000000_kr, 0.5000000000000000_kr, 0.0000000000000000_kr]
            Xi(17) = [0.5000000000000000_kr, 0.0000000000000000_kr, 0.0000000000000000_kr]
            Xi(18) = [0.0000000000000000_kr, 0.5000000000000000_kr, 0.0000000000000000_kr]
            Xi(19) = [0.0000000000000000_kr, 0.0000000000000000_kr, 0.5000000000000000_kr]
            Xi(20) = [0.2000000000000000_kr, 0.1000000000000000_kr, 0.1000000000000000_kr]
            Xi(21) = [0.1000000000000000_kr, 0.2000000000000000_kr, 0.1000000000000000_kr]
            Xi(22) = [0.1000000000000000_kr, 0.1000000000000000_kr, 0.2000000000000000_kr]
            Xi(23) = [0.6000000000000000_kr, 0.1000000000000000_kr, 0.1000000000000000_kr]
            Xi(24) = [0.1000000000000000_kr, 0.6000000000000000_kr, 0.1000000000000000_kr]
            Xi(25) = [0.1000000000000000_kr, 0.1000000000000000_kr, 0.6000000000000000_kr]
            Xi(26) = [0.1000000000000000_kr, 0.2000000000000000_kr, 0.6000000000000000_kr]
            Xi(27) = [0.2000000000000000_kr, 0.6000000000000000_kr, 0.1000000000000000_kr]
            Xi(28) = [0.6000000000000000_kr, 0.1000000000000000_kr, 0.2000000000000000_kr]
            Xi(29) = [0.1000000000000000_kr, 0.6000000000000000_kr, 0.2000000000000000_kr]
            Xi(30) = [0.2000000000000000_kr, 0.1000000000000000_kr, 0.6000000000000000_kr]
            Xi(31) = [0.6000000000000000_kr, 0.2000000000000000_kr, 0.1000000000000000_kr]
            dElem%Gauss_C(1) = 0.1095853407966528_kr / 6.0_kr * detBinv
            dElem%Gauss_C(2:5) = 0.0635996491464850_kr / 6.0_kr * detBinv
            dElem%Gauss_C(6:9) = -0.3751064406859797_kr / 6.0_kr * detBinv
            dElem%Gauss_C(10:13) = 0.0293485515784412_kr / 6.0_kr * detBinv
            dElem%Gauss_C(14:19) = 0.0058201058201058_kr / 6.0_kr * detBinv
            dElem%Gauss_C(20:31) = 0.1653439153439105_kr / 6.0_kr * detBinv
         case (8)
            Nb_Gauss = 45
            allocate (Xi(Nb_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
            Xi(1) = [0.2500000000000000_kr, 0.2500000000000000_kr, 0.2500000000000000_kr]
            Xi(2) = [0.6175871903000830_kr, 0.1274709365666390_kr, 0.1274709365666390_kr]
            Xi(3) = [0.1274709365666390_kr, 0.1274709365666390_kr, 0.1274709365666390_kr]
            Xi(4) = [0.1274709365666390_kr, 0.1274709365666390_kr, 0.6175871903000830_kr]
            Xi(5) = [0.1274709365666390_kr, 0.6175871903000830_kr, 0.1274709365666390_kr]
            Xi(6) = [0.9037635088221031_kr, 0.0320788303926323_kr, 0.0320788303926323_kr]
            Xi(7) = [0.0320788303926323_kr, 0.0320788303926323_kr, 0.0320788303926323_kr]
            Xi(8) = [0.0320788303926323_kr, 0.0320788303926323_kr, 0.9037635088221031_kr]
            Xi(9) = [0.0320788303926323_kr, 0.9037635088221031_kr, 0.0320788303926323_kr]
            Xi(10) = [0.4502229043567190_kr, 0.0497770956432810_kr, 0.0497770956432810_kr]
            Xi(11) = [0.0497770956432810_kr, 0.4502229043567190_kr, 0.0497770956432810_kr]
            Xi(12) = [0.0497770956432810_kr, 0.0497770956432810_kr, 0.4502229043567190_kr]
            Xi(13) = [0.0497770956432810_kr, 0.4502229043567190_kr, 0.4502229043567190_kr]
            Xi(14) = [0.4502229043567190_kr, 0.0497770956432810_kr, 0.4502229043567190_kr]
            Xi(15) = [0.4502229043567190_kr, 0.4502229043567190_kr, 0.0497770956432810_kr]
            Xi(16) = [0.3162695526014501_kr, 0.1837304473985499_kr, 0.1837304473985499_kr]
            Xi(17) = [0.1837304473985499_kr, 0.3162695526014501_kr, 0.1837304473985499_kr]
            Xi(18) = [0.1837304473985499_kr, 0.1837304473985499_kr, 0.3162695526014501_kr]
            Xi(19) = [0.1837304473985499_kr, 0.3162695526014501_kr, 0.3162695526014501_kr]
            Xi(20) = [0.3162695526014501_kr, 0.1837304473985499_kr, 0.3162695526014501_kr]
            Xi(21) = [0.3162695526014501_kr, 0.3162695526014501_kr, 0.1837304473985499_kr]
            Xi(22) = [0.0229177878448171_kr, 0.2319010893971509_kr, 0.2319010893971509_kr]
            Xi(23) = [0.2319010893971509_kr, 0.0229177878448171_kr, 0.2319010893971509_kr]
            Xi(24) = [0.2319010893971509_kr, 0.2319010893971509_kr, 0.0229177878448171_kr]
            Xi(25) = [0.5132800333608811_kr, 0.2319010893971509_kr, 0.2319010893971509_kr]
            Xi(26) = [0.2319010893971509_kr, 0.5132800333608811_kr, 0.2319010893971509_kr]
            Xi(27) = [0.2319010893971509_kr, 0.2319010893971509_kr, 0.5132800333608811_kr]
            Xi(28) = [0.2319010893971509_kr, 0.0229177878448171_kr, 0.5132800333608811_kr]
            Xi(29) = [0.0229177878448171_kr, 0.5132800333608811_kr, 0.2319010893971509_kr]
            Xi(30) = [0.5132800333608811_kr, 0.2319010893971509_kr, 0.0229177878448171_kr]
            Xi(31) = [0.2319010893971509_kr, 0.5132800333608811_kr, 0.0229177878448171_kr]
            Xi(32) = [0.0229177878448171_kr, 0.2319010893971509_kr, 0.5132800333608811_kr]
            Xi(33) = [0.5132800333608811_kr, 0.0229177878448171_kr, 0.2319010893971509_kr]
            Xi(34) = [0.7303134278075384_kr, 0.0379700484718286_kr, 0.0379700484718286_kr]
            Xi(35) = [0.0379700484718286_kr, 0.7303134278075384_kr, 0.0379700484718286_kr]
            Xi(36) = [0.0379700484718286_kr, 0.0379700484718286_kr, 0.7303134278075384_kr]
            Xi(37) = [0.1937464752488044_kr, 0.0379700484718286_kr, 0.0379700484718286_kr]
            Xi(38) = [0.0379700484718286_kr, 0.1937464752488044_kr, 0.0379700484718286_kr]
            Xi(39) = [0.0379700484718286_kr, 0.0379700484718286_kr, 0.1937464752488044_kr]
            Xi(40) = [0.0379700484718286_kr, 0.7303134278075384_kr, 0.1937464752488044_kr]
            Xi(41) = [0.7303134278075384_kr, 0.1937464752488044_kr, 0.0379700484718286_kr]
            Xi(42) = [0.1937464752488044_kr, 0.0379700484718286_kr, 0.7303134278075384_kr]
            Xi(43) = [0.0379700484718286_kr, 0.1937464752488044_kr, 0.7303134278075384_kr]
            Xi(44) = [0.7303134278075384_kr, 0.0379700484718286_kr, 0.1937464752488044_kr]
            Xi(45) = [0.1937464752488044_kr, 0.7303134278075384_kr, 0.0379700484718286_kr]
            dElem%Gauss_C(1) = -0.2359620398477557_kr / 6.0_kr * detBinv
            dElem%Gauss_C(2:5) = 0.0244878963560562_kr / 6.0_kr * detBinv
            dElem%Gauss_C(6:9) = 0.0039485206398261_kr / 6.0_kr * detBinv
            dElem%Gauss_C(10:15) = 0.0263055529507371_kr / 6.0_kr * detBinv
            dElem%Gauss_C(16:21) = 0.0829803830550589_kr / 6.0_kr * detBinv
            dElem%Gauss_C(22:33) = 0.0254426245481023_kr / 6.0_kr * detBinv
            dElem%Gauss_C(34:45) = 0.0134324384376852_kr / 6.0_kr * detBinv
         case Default
            write (*, *) __FUNCT__, ': Unimplemented quadrature order', dQuadratureOrder
            stop
         end select

         select case (dPolynomialOrder)
         case (1)
            Num_DoF = 4
            allocate (PhiHat(Num_DoF, Nb_Gauss), stat=ierr)
            allocate (GradPhiHat(Num_DoF, Nb_Gauss), stat=ierr)
            PhiHat(1, :) = 1.0_kr - Xi%X - Xi%Y - Xi%Z
            PhiHat(2, :) = Xi(:)%Y
            PhiHat(3, :) = Xi(:)%X
            PhiHat(4, :) = Xi(:)%Z

            GradPhiHat(1, :)%X = -1.0_kr; GradPhiHat(1, :)%Y = -1.0_kr; GradPhiHat(1, :)%Z = -1.0_kr
            GradPhiHat(2, :)%X = 0.0_kr; GradPhiHat(2, :)%Y = 1.0_kr; GradPhiHat(2, :)%Z = 0.0_kr
            GradPhiHat(3, :)%X = 1.0_kr; GradPhiHat(3, :)%Y = 0.0_kr; GradPhiHat(3, :)%Z = 0.0_kr
            GradPhiHat(4, :)%X = 0.0_kr; GradPhiHat(4, :)%Y = 0.0_kr; GradPhiHat(4, :)%Z = 1.0_kr
         case (2)
            Num_Dof = 10
            allocate (PhiHat(Num_DoF, Nb_Gauss), stat=ierr)
            allocate (GradPhiHat(Num_DoF, Nb_Gauss), stat=ierr)
            PhiHat(1, :) = 4.0_kr * (1.0_kr - Xi%X - Xi%Y - Xi%Z) * Xi%Y
            PhiHat(2, :) = 4.0_kr * Xi%X * Xi%Y
            PhiHat(3, :) = 4.0_kr * (1.0_kr - Xi%X - Xi%Y - Xi%Z) * Xi%X
            PhiHat(4, :) = 4.0_kr * (1.0_kr - Xi%X - Xi%Y - Xi%Z) * Xi%Z
            PhiHat(5, :) = 4.0_kr * Xi%Y * Xi%Z
            PhiHat(6, :) = 4.0_kr * Xi%X * Xi%Z
            PhiHat(7, :) = (1.0_kr - Xi%X - Xi%Y - Xi%Z) * (1.0_kr - 2.0_kr * Xi%X - 2.0_kr * Xi%Y - 2.0_kr * Xi%Z)
            PhiHat(8, :) = Xi%Y * (2.0_kr * Xi%Y - 1.0_kr)
            PhiHat(9, :) = Xi%X * (2.0_kr * Xi%X - 1.0_kr)
            PhiHat(10, :) = Xi%Z * (2.0_kr * Xi%Z - 1.0_kr)

            GradPhiHat(1, :)%X = -4.0_kr * Xi%Y
            GradPhiHat(1, :)%Y = 4.0_kr * (1.0_kr - Xi%X - 2.0_kr * Xi%Y - Xi%Z)
            GradPhiHat(1, :)%Z = -4.0_kr * Xi%Y

            GradPhiHat(2, :)%X = 4.0_kr * Xi%Y
            GradPhiHat(2, :)%Y = 4.0_kr * Xi%X
            GradPhiHat(2, :)%Z = 0.0_kr

            GradPhiHat(3, :)%X = 4.0_kr * (1.0_kr - 2.0_kr * Xi%X - Xi%Y - Xi%Z)
            GradPhiHat(3, :)%Y = -4.0_kr * Xi%X
            GradPhiHat(3, :)%Z = -4.0_kr * Xi%X

            GradPhiHat(4, :)%X = -4.0_kr * Xi%Z
            GradPhiHat(4, :)%Y = -4.0_kr * Xi%Z
            GradPhiHat(4, :)%Z = 4.0_kr * (1.0_kr - Xi%X - Xi%Y - 2.0_kr * Xi%Z)

            GradPhiHat(5, :)%X = 0.0_kr
            GradPhiHat(5, :)%Y = 4.0_kr * Xi%Z
            GradPhiHat(5, :)%Z = 4.0_kr * Xi%Y

            GradPhiHat(6, :)%X = 4.0_kr * Xi%Z
            GradPhiHat(6, :)%Y = 0.0_kr
            GradPhiHat(6, :)%Z = 4.0_kr * Xi%X

            GradPhiHat(7, :)%X = 4.0_kr * Xi%X + 4.0_kr * Xi%Y + 4.0_kr * Xi%Z - 3.0_kr
            GradPhiHat(7, :)%Y = 4.0_kr * Xi%X + 4.0_kr * Xi%Y + 4.0_kr * Xi%Z - 3.0_kr
            GradPhiHat(7, :)%Z = 4.0_kr * Xi%X + 4.0_kr * Xi%Y + 4.0_kr * Xi%Z - 3.0_kr

            GradPhiHat(8, :)%X = 0.0_kr
            GradPhiHat(8, :)%Y = 4.0_kr * Xi%Y - 1.0_kr
            GradPhiHat(8, :)%Z = 0.0_kr

            GradPhiHat(9, :)%X = 4.0_kr * Xi%X - 1.0_kr
            GradPhiHat(9, :)%Y = 0.0_kr
            GradPhiHat(9, :)%Z = 0.0_kr

            GradPhiHat(10, :)%X = 0.0_kr
            GradPhiHat(10, :)%Y = 0.0_kr
            GradPhiHat(10, :)%Z = 4.0_kr * Xi%Z - 1.0_kr

         case Default
            Num_DoF = 0
            write (*, *) __FUNCT__, ': Unimplemented PolynomialOrder', dPolynomialOrder
            stop
         end select

         allocate (dElem%BF(Num_DoF, Nb_Gauss), stat=ierr)
         allocate (dElem%Grad_BF(Num_DoF, Nb_Gauss), stat=ierr)
         dElem%BF = PhiHat
         do iDoF = 1, Num_DoF
            do iG = 1, Nb_Gauss
               dElem%Grad_BF(iDoF, iG) = Bt * GradPhiHat(iDoF, iG)
            end do
         end do

         deallocate (Xi, stat=ierr)
         deallocate (PhiHat, stat=ierr)
         deallocate (GradPhiHat, stat=ierr)
      end subroutine ElementPLagrange3DScalInit

#undef __FUNCT__
#define __FUNCT__ "ElementPLagrange3DBoundaryScalInit"
      subroutine ElementPLagrange3DBoundaryScalInit(dElem, area, outerNormal, dPolynomialOrder, dQuadratureOrder, ierr)
         ! Compute the quadrature weights and the value of the basis functions and their gradient
         ! at the quadrature points.
         type(MEF90Element3DScal), intent(INOUT) :: dElem
         PetscReal, intent(IN)                   :: area
         type(Vect3D), intent(IN)                :: outerNormal
         PetscInt, intent(IN)                    :: dPolynomialOrder, dQuadratureOrder
         PetscErrorCode, intent(OUT)             :: ierr

         PetscInt                               :: Nb_Gauss
         PetscInt                               :: Num_Dof
         PetscInt                               :: iDoF, iG

         type(Vect2D), dimension(:), pointer      :: Xi ! The quadrature points coordinates in the reference element

         num_Dof = 0
         nb_Gauss = 0

         dElem%outerNormal = outerNormal

         select case (dQuadratureOrder)
         case (0, 1)
            Nb_Gauss = 1
            allocate (Xi(Nb_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
            Xi(1)%X = 1.0_kr / 3.0_kr
            Xi(1)%Y = 1.0_kr / 3.0_kr
            dElem%Gauss_C = area

         case (2)
            Nb_Gauss = 3
            allocate (Xi(Nb_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
            Xi(1) = [1.0_kr / 6.0_kr, 1.0_kr / 6.0_kr]
            Xi(2) = [2.0_kr / 3.0_kr, 1.0_kr / 6.0_kr]
            Xi(3) = [1.0_kr / 6.0_kr, 2.0_kr / 3.0_kr]
            dElem%Gauss_C = area / 3.0_kr

         case (3)
            Nb_Gauss = 4
            allocate (Xi(Nb_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
            dElem%Gauss_C = area * 25.0_kr / 48.0_kr
            dElem%Gauss_C(1) = -area * 9.0_kr / 16.0_kr
            Xi(1) = [1.0_kr / 3.0_kr, 1.0_kr / 3.0_kr]
            Xi(2) = [3.0_kr / 5.0_kr, 1.0_kr / 5.0_kr]
            Xi(3) = [1.0_kr / 5.0_kr, 3.0_kr / 5.0_kr]
            Xi(4) = [1.0_kr / 5.0_kr, 1.0_kr / 5.0_kr]
         case (4)
            Nb_Gauss = 7
            allocate (Xi(Nb_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
            dElem%Gauss_C(1) = area / 20.0_kr
            dElem%Gauss_C(2) = area * 2.0_kr / 15.0_kr
            dElem%Gauss_C(3) = area / 20.0_kr
            dElem%Gauss_C(4) = area * 2.0_kr / 15.0_kr
            dElem%Gauss_C(5) = area / 20.0_kr
            dElem%Gauss_C(6) = area * 2.0_kr / 15.0_kr
            dElem%Gauss_C(7) = area * 9.0_kr / 20.0_kr
            Xi(1) = [0.0_kr, 0.0_kr]
            Xi(2) = [0.5_kr, 0.0_kr]
            Xi(3) = [1.0_kr, 0.0_kr]
            Xi(4) = [0.5_kr, 0.5_kr]
            Xi(5) = [0.0_kr, 1.0_kr]
            Xi(6) = [0.0_kr, 0.5_kr]
            Xi(7) = [1.0_kr / 3.0_kr, 1.0_kr / 3.0_kr]
         case (5, 6)
            Nb_Gauss = 9
            allocate (Xi(Nb_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
            Xi(1) = [0.124949503233232_kr, 0.437525248383384_kr]
            Xi(2) = [0.437525248383384_kr, 0.124949503233232_kr]
            Xi(3) = [0.437525248383384_kr, 0.437525248383384_kr]
            Xi(4) = [0.797112651860071_kr, 0.165409927389841_kr]
            Xi(5) = [0.797112651860071_kr, 0.037477420750088_kr]
            Xi(6) = [0.165409927389841_kr, 0.797112651860071_kr]
            Xi(7) = [0.165409927389841_kr, 0.037477420750088_kr]
            Xi(8) = [0.037477420750088_kr, 0.797112651860071_kr]
            Xi(9) = [0.037477420750088_kr, 0.165409927389841_kr]
            dElem%Gauss_C(1:3) = 0.205950504760887_kr * area
            dElem%Gauss_C(4:9) = 0.063691414286223_kr * area
         case (7)
            Nb_Gauss = 13
            allocate (Xi(Nb_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
            Xi(1) = [0.333333333333333_kr, 0.333333333333333_kr]
            Xi(2) = [0.479308067841923_kr, 0.260345966079038_kr]
            Xi(3) = [0.260345966079038_kr, 0.479308067841923_kr]
            Xi(4) = [0.260345966079038_kr, 0.260345966079038_kr]
            Xi(5) = [0.869739794195568_kr, 0.065130102902216_kr]
            Xi(6) = [0.065130102902216_kr, 0.869739794195568_kr]
            Xi(7) = [0.065130102902216_kr, 0.065130102902216_kr]
            Xi(8) = [0.638444188569809_kr, 0.312865496004875_kr]
            Xi(9) = [0.638444188569809_kr, 0.048690315425316_kr]
            Xi(10) = [0.312865496004875_kr, 0.638444188569809_kr]
            Xi(11) = [0.312865496004875_kr, 0.048690315425316_kr]
            Xi(12) = [0.048690315425316_kr, 0.638444188569809_kr]
            Xi(13) = [0.048690315425316_kr, 0.312865496004875_kr]
            dElem%Gauss_C(1) = -0.149570044467670_kr * area
            dElem%Gauss_C(2:4) = 0.175615257433204_kr * area
            dElem%Gauss_C(5:7) = 0.053347235608839_kr * area
            dElem%Gauss_C(8:13) = 0.077113760890257_kr * area
         case (8, 9)
            Nb_Gauss = 19
            allocate (Xi(Nb_Gauss), stat=ierr)
            allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
            Xi(1) = [0.3333333333333333_kr, 0.3333333333333333_kr]
            Xi(2) = [0.0206349616025259_kr, 0.4896825191987370_kr]
            Xi(3) = [0.4896825191987370_kr, 0.0206349616025259_kr]
            Xi(4) = [0.4896825191987370_kr, 0.4896825191987370_kr]
            Xi(5) = [0.1258208170141290_kr, 0.4370895914929355_kr]
            Xi(6) = [0.4370895914929355_kr, 0.1258208170141290_kr]
            Xi(7) = [0.4370895914929355_kr, 0.4370895914929355_kr]
            Xi(8) = [0.6235929287619356_kr, 0.1882035356190322_kr]
            Xi(9) = [0.1882035356190322_kr, 0.6235929287619356_kr]
            Xi(10) = [0.1882035356190322_kr, 0.1882035356190322_kr]
            Xi(11) = [0.9105409732110941_kr, 0.0447295133944530_kr]
            Xi(12) = [0.0447295133944530_kr, 0.9105409732110941_kr]
            Xi(13) = [0.0447295133944530_kr, 0.0447295133944530_kr]
            Xi(14) = [0.7411985987844980_kr, 0.0368384120547363_kr]
            Xi(15) = [0.7411985987844980_kr, 0.2219628891607657_kr]
            Xi(16) = [0.0368384120547363_kr, 0.7411985987844980_kr]
            Xi(17) = [0.0368384120547363_kr, 0.2219628891607657_kr]
            Xi(18) = [0.2219628891607657_kr, 0.7411985987844980_kr]
            Xi(19) = [0.2219628891607657_kr, 0.0368384120547363_kr]
            dElem%Gauss_C(1) = 0.09713579628279610_kr * area
            dElem%Gauss_C(2:4) = 0.03133470022713983_kr * area
            dElem%Gauss_C(5:7) = 0.07782754100477543_kr * area
            dElem%Gauss_C(8:10) = 0.07964773892720910_kr * area
            dElem%Gauss_C(11:13) = 0.02557767565869810_kr * area
            dElem%Gauss_C(14:19) = 0.04328353937728940_kr * area
         case Default
            write (*, *) __FUNCT__, ': Unimplemented quadrature order', dQuadratureOrder
            ierr = PETSC_ERR_SUP
            stop
         end select

         select case (dPolynomialOrder)
         case (1)
            Num_DoF = 3
            allocate (dElem%BF(Num_DoF, Nb_Gauss), stat=ierr)
            dElem%BF(1, :) = 1.0_kr - Xi%X - Xi%Y
            dElem%BF(2, :) = Xi(:)%X
            dElem%BF(3, :) = Xi(:)%Y
         case (2)
            Num_DoF = 6
            allocate (dElem%BF(Num_DoF, Nb_Gauss), stat=ierr)
            dElem%BF(5, :) = (1.0_kr - Xi%X - Xi%Y) * (1.0_kr - 2.0_kr * Xi%X - 2.0_kr * Xi%Y)
            dElem%BF(6, :) = Xi%X * (2.0_kr * Xi%X - 1.0_kr)
            dElem%BF(4, :) = Xi%Y * (2.0_kr * Xi%Y - 1.0_kr)
            dElem%BF(2, :) = 4.0_kr * Xi%X * (1.0_kr - Xi%X - Xi%Y)
            dElem%BF(3, :) = 4.0_kr * Xi%X * Xi%Y
            dElem%BF(1, :) = 4.0_kr * Xi%Y * (1.0_kr - Xi%X - Xi%Y)
         case Default
            write (*, *) __FUNCT__, ': Unimplemented PolynomialOrder', dPolynomialOrder
            ierr = PETSC_ERR_SUP
         end select

         allocate (delem%Grad_BF(Num_DoF, Nb_Gauss), stat=ierr)
         do iDof = 1, num_Dof
            do iG = 1, Nb_Gauss
               delem%Grad_BF(iDof, iG) = 0.0_kr
            end do
         end do

         deallocate (Xi, stat=ierr)
      end subroutine ElementPLagrange3DBoundaryScalInit

#undef __FUNCT__
#define __FUNCT__ "ElementPLagrange3DVectInit"
      subroutine ElementPLagrange3DVectInit(dElem, Bt, detBinv, dPolynomialOrder, dQuadratureOrder, ierr)
         type(MEF90Element3DVect), intent(INOUT) :: dElem
         type(Mat3D), intent(IN)                 :: Bt
         PetscReal, intent(IN)                   :: detBinv
         PetscInt, intent(IN)                    :: dPolynomialOrder, dQuadratureOrder
         PetscErrorCode, intent(OUT)             :: ierr

         type(MEF90Element3DScal)              :: Elem_Scal
         PetscInt                               :: dim = 3
         PetscInt                               :: Num_DoF, Nb_Gauss, i, iDof, iG

         call ElementPLagrange3DScalInit(Elem_Scal, Bt, detBinv, dPolynomialOrder, dQuadratureOrder, ierr)
         Num_DoF = size(Elem_Scal%BF, 1)
         Nb_Gauss = size(Elem_Scal%BF, 2)
         allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
         allocate (dElem%BF(Num_DoF * dim, Nb_Gauss), stat=ierr)
         allocate (dElem%Grad_BF(Num_DoF * dim, Nb_Gauss), stat=ierr)
         allocate (dElem%GradS_BF(Num_DoF * dim, Nb_Gauss), stat=ierr)

         dElem%Gauss_C = Elem_Scal%Gauss_C
         do iDof = 1, num_dof * dim
            do iG = 1, Nb_Gauss
               dElem%BF(iDof, iG) = 0.0_kr
               delem%Grad_BF(iDof, iG) = 0.0_kr
               delem%GradS_BF(iDof, iG) = 0.0_kr
            end do
         end do

         do i = 0, Num_DoF - 1
            dElem%BF(i * dim + 1, :)%X = Elem_Scal%BF(i + 1, :)
            dElem%BF(i * dim + 2, :)%Y = Elem_Scal%BF(i + 1, :)
            dElem%BF(i * dim + 3, :)%Z = Elem_Scal%BF(i + 1, :)
            dElem%Grad_BF(i * dim + 1, :)%XX = Elem_Scal%Grad_BF(i + 1, :)%X
            dElem%Grad_BF(i * dim + 1, :)%XY = Elem_Scal%Grad_BF(i + 1, :)%Y
            dElem%Grad_BF(i * dim + 1, :)%XZ = Elem_Scal%Grad_BF(i + 1, :)%Z
            dElem%Grad_BF(i * dim + 2, :)%YX = Elem_Scal%Grad_BF(i + 1, :)%X
            dElem%Grad_BF(i * dim + 2, :)%YY = Elem_Scal%Grad_BF(i + 1, :)%Y
            dElem%Grad_BF(i * dim + 2, :)%YZ = Elem_Scal%Grad_BF(i + 1, :)%Z
            dElem%Grad_BF(i * dim + 3, :)%ZX = Elem_Scal%Grad_BF(i + 1, :)%X
            dElem%Grad_BF(i * dim + 3, :)%ZY = Elem_Scal%Grad_BF(i + 1, :)%Y
            dElem%Grad_BF(i * dim + 3, :)%ZZ = Elem_Scal%Grad_BF(i + 1, :)%Z

            dElem%GradS_BF(i * dim + 1, :)%XX = Elem_Scal%Grad_BF(i + 1, :)%X
            dElem%GradS_BF(i * dim + 1, :)%XY = Elem_Scal%Grad_BF(i + 1, :)%Y * 0.5_kr
            dElem%GradS_BF(i * dim + 1, :)%XZ = Elem_Scal%Grad_BF(i + 1, :)%Z * 0.5_kr

            dElem%GradS_BF(i * dim + 2, :)%XY = Elem_Scal%Grad_BF(i + 1, :)%X * 0.5_kr
            dElem%GradS_BF(i * dim + 2, :)%YY = Elem_Scal%Grad_BF(i + 1, :)%Y
            dElem%GradS_BF(i * dim + 2, :)%YZ = Elem_Scal%Grad_BF(i + 1, :)%Z * 0.5_kr

            dElem%GradS_BF(i * dim + 3, :)%XZ = Elem_Scal%Grad_BF(i + 1, :)%X * 0.5_kr
            dElem%GradS_BF(i * dim + 3, :)%YZ = Elem_Scal%Grad_BF(i + 1, :)%Y * 0.5_kr
            dElem%GradS_BF(i * dim + 3, :)%ZZ = Elem_Scal%Grad_BF(i + 1, :)%Z
         end do
         call MEF90ElementDestroy(Elem_Scal, ierr)
      end subroutine ElementPLagrange3DVectInit

#undef __FUNCT__
#define __FUNCT__ "ElementPLagrange3DBoundaryVectInit"
      subroutine ElementPLagrange3DBoundaryVectInit(dElem, area, outerNormal, dPolynomialOrder, dQuadratureOrder, ierr)
         type(MEF90Element3DVect), intent(INOUT) :: dElem
         PetscReal, intent(IN)                   :: area
         type(Vect3D), intent(IN)                :: outerNormal
         PetscInt                               :: dPolynomialOrder, dQuadratureOrder
         PetscErrorCode, intent(OUT)             :: ierr

         type(MEF90Element3DScal)              :: Elem_Scal
         PetscInt                               :: dim = 3
         PetscInt                               :: Num_DoF, Nb_Gauss, iDof, iG

         call ElementPLagrange3DBoundaryScalInit(Elem_Scal, area, outerNormal, dPolynomialOrder, dQuadratureOrder, ierr)
         Num_DoF = size(Elem_Scal%BF, 1)
         Nb_Gauss = size(Elem_Scal%BF, 2)
         allocate (dElem%Gauss_C(Nb_Gauss), stat=ierr)
         allocate (dElem%BF(Num_DoF * dim, Nb_Gauss), stat=ierr)
         allocate (dElem%Grad_BF(Num_DoF * dim, Nb_Gauss), stat=ierr)
         allocate (dElem%GradS_BF(Num_DoF * dim, Nb_Gauss), stat=ierr)

         dElem%Gauss_C = Elem_Scal%Gauss_C
         do iDof = 1, num_dof * dim
            do iG = 1, Nb_Gauss
               dElem%BF(iDof, iG) = 0.0_kr
               delem%Grad_BF(iDof, iG) = 0.0_kr
               delem%GradS_BF(iDof, iG) = 0.0_kr
            end do
         end do
         do iDof = 0, Num_DoF - 1
            dElem%BF(iDof * dim + 1, :)%X = Elem_Scal%BF(iDof + 1, :)
            dElem%BF(iDof * dim + 2, :)%Y = Elem_Scal%BF(iDof + 1, :)
            dElem%BF(iDof * dim + 3, :)%Z = Elem_Scal%BF(iDof + 1, :)
         end do

         dElem%outerNormal = Elem_Scal%outerNormal
         call MEF90ElementDestroy(Elem_Scal, ierr)
      end subroutine ElementPLagrange3DBoundaryVectInit

#undef __FUNCT__
#define __FUNCT__ "Element2DScalDestroy"
      subroutine Element2DScalDestroy(dElem, ierr)
         type(MEF90Element2DScal), intent(INOUT) :: dElem
         PetscErrorCode, intent(OUT)             :: ierr

         if (associated(dElem%BF)) then
            deallocate (dElem%BF, stat=ierr)
         end if
         if (associated(dElem%Grad_BF)) then
            deallocate (dElem%Grad_BF, stat=ierr)
         end if
         if (associated(dElem%Gauss_C)) then
            deallocate (dElem%Gauss_C, stat=ierr)
         end if
      end subroutine Element2DScalDestroy

#undef __FUNCT__
#define __FUNCT__ "Element2DVectDestroy"
      subroutine Element2DVectDestroy(dElem, ierr)
         type(MEF90Element2DVect), intent(INOUT) :: dElem
         PetscErrorCode, intent(OUT)             :: ierr

         if (associated(dElem%BF)) then
            deallocate (dElem%BF, stat=ierr)
         end if
         if (associated(dElem%Grad_BF)) then
            deallocate (dElem%Grad_BF, stat=ierr)
         end if
         if (associated(dElem%GradS_BF)) then
            deallocate (dElem%GradS_BF, stat=ierr)
         end if
         if (associated(dElem%Gauss_C)) then
            deallocate (dElem%Gauss_C, stat=ierr)
         end if
      end subroutine Element2DVectDestroy

#undef __FUNCT__
#define __FUNCT__ "Element3DScalDestroy"
      subroutine Element3DScalDestroy(dElem, ierr)
         type(MEF90Element3DScal), intent(INOUT) :: dElem
         PetscErrorCode, intent(OUT)             :: ierr

         if (associated(dElem%BF)) then
            deallocate (dElem%BF, stat=ierr)
         end if
         if (associated(dElem%Grad_BF)) then
            deallocate (dElem%Grad_BF, stat=ierr)
         end if
         if (associated(dElem%Gauss_C)) then
            deallocate (dElem%Gauss_C, stat=ierr)
         end if
      end subroutine Element3DScalDestroy

#undef __FUNCT__
#define __FUNCT__ "Element3DVectDestroy"
      subroutine Element3DVectDestroy(dElem, ierr)
         type(MEF90Element3DVect), intent(INOUT) :: dElem
         PetscErrorCode, intent(OUT)             :: ierr
         if (associated(dElem%BF)) then
            deallocate (dElem%BF, stat=ierr)
         end if
         if (associated(dElem%Grad_BF)) then
            deallocate (dElem%Grad_BF, stat=ierr)
         end if
         if (associated(dElem%GradS_BF)) then
            deallocate (dElem%GradS_BF, stat=ierr)
         end if
         if (associated(dElem%Gauss_C)) then
            deallocate (dElem%Gauss_C, stat=ierr)
         end if
      end subroutine Element3DVectDestroy

#undef __FUNCT__
#define __FUNCT__ "Element2DScalDestroySet"
      subroutine Element2DScalDestroySet(dElem, ierr)
         type(MEF90Element2DScal), dimension(:), pointer    :: dElem
         PetscErrorCode, intent(OUT)                        :: ierr

         PetscInt                                          :: cell

         do cell = 1, size(dElem)
            call MEF90ElementDestroy(dElem(cell), ierr)
         end do
         deallocate (dElem, stat=ierr)
      end subroutine Element2DScalDestroySet

#undef __FUNCT__
#define __FUNCT__ "Element2DVectDestroySet"
      subroutine Element2DVectDestroySet(dElem, ierr)
         type(MEF90Element2DVect), dimension(:), pointer    :: dElem
         PetscErrorCode, intent(OUT)                        :: ierr

         PetscInt                                          :: cell

         do cell = 1, size(dElem)
            call MEF90ElementDestroy(dElem(cell), ierr)
         end do
         deallocate (dElem, stat=ierr)
      end subroutine Element2DVectDestroySet

#undef __FUNCT__
#define __FUNCT__ "Element3DScalDestroySet"
      subroutine Element3DScalDestroySet(dElem, ierr)
         type(MEF90Element3DScal), dimension(:), pointer    :: dElem
         PetscErrorCode, intent(OUT)                        :: ierr

         PetscInt                                          :: cell

         do cell = 1, size(dElem)
            call MEF90ElementDestroy(dElem(cell), ierr)
         end do
         deallocate (dElem, stat=ierr)
      end subroutine Element3DScalDestroySet

#undef __FUNCT__
#define __FUNCT__ "Element3DVectDestroySet"
      subroutine Element3DVectDestroySet(dElem, ierr)
         type(MEF90Element3DVect), dimension(:), pointer    :: dElem
         PetscErrorCode, intent(OUT)                        :: ierr

         PetscInt                                          :: cell

         do cell = 1, size(dElem)
            call MEF90ElementDestroy(dElem(cell), ierr)
         end do
         deallocate (dElem, stat=ierr)
      end subroutine Element3DVectDestroySet

   end module m_MEF90_Elements
