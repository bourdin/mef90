#include "../MEF90/mef90.inc"
program TestQuadrature
#include "petsc/finclude/petsc.h"
   use m_MEF90
   use petsc
   implicit none(type, external)

   PetscErrorCode                      :: ierr
   type(tDM), target                    :: dm, dmU
   type(tPetscSection), target          :: sectionU
   character(len=MEF90MXSTRLEN)        :: IOBuffer
   PetscInt                            :: dim
   PetscBool                           :: flg
   PetscInt                            :: i, j, k

   PetscReal                           :: scal, sol, xr = 1.0_kr, xl = 0.0_kr, yr = 1.0_kr, yl = 0.0_kr, zr = 1.0_kr, zl = 0.0_kr
   character(len=MEF90MXSTRLEN)        :: name
   type(tVec)                          :: locVecU
   type(Vect2D)                        :: v2D
   type(Vect3D)                        :: v3D
   PetscInt                            :: QuadOrderMax, QuadOrder

   type(MEF90Ctx_Type), target                         :: MEF90Ctx
   type(MEF90CtxGlobalOptions_Type), pointer           :: MEF90GlobalOptions
   type(MEF90CtxGlobalOptions_Type)                   :: MEF90GlobalOptions_default

   !!! Initialize MEF90
   PetscCallA(PetscInitialize(ierr))
   PetscCallA(MEF90Initialize(PETSC_COMM_WORLD, ierr))

   MEF90GlobalOptions_default%verbose = 1
   MEF90GlobalOptions_default%dryrun = PETSC_FALSE
   MEF90GlobalOptions_default%timeMin = 0.0_kr
   MEF90GlobalOptions_default%timeMax = 1.0_kr
   MEF90GlobalOptions_default%timeNumStep = 1
   MEF90GlobalOptions_default%timeSkip = 0
   MEF90GlobalOptions_default%timeNumCycle = 1
   MEF90GlobalOptions_default%timeInterpolation = MEF90TimeInterpolation_linear
   MEF90GlobalOptions_default%elementFamily = MEF90ElementFamilyLagrange
   MEF90GlobalOptions_default%elementOrder = 1

   !!! Get all MEF90-wide options
   PetscCallA(MEF90CtxCreate(PETSC_COMM_WORLD, MEF90Ctx, MEF90GlobalOptions_default, ierr))
   PetscCallA(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag, MEF90GlobalOptions, ierr))

   PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm, MEF90Ctx%geometryFile, PETSC_NULL_CHARACTER, PETSC_TRUE, dm, ierr))
   PetscCallA(DMSetFromOptions(dm, ierr))
   PetscCallA(DMViewFromOptions(dm, PETSC_NULL_OPTIONS, "-dm_view", ierr))
   PetscCall(DMGetDimension(dm, dim, ierr))

   QuadOrderMax = 4
   PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-order", QuadOrderMax, flg, ierr))
   PetscCallA(PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-xl", xl, flg, ierr))
   PetscCallA(PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-xr", xr, flg, ierr))
   PetscCallA(PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-yl", yl, flg, ierr))
   PetscCallA(PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-yr", yr, flg, ierr))
   PetscCallA(PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-zl", zl, flg, ierr))
   PetscCallA(PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-zr", zr, flg, ierr))

   ! Create nodal local Vec holding coordinates
   name = "U"
   PetscCallA(MEF90CreateLocalVector(dm, MEF90GlobalOptions%elementFamily, MEF90GlobalOptions%elementOrder, dim, name, locVecU, ierr))
   PetscCallA(VecGetDM(locVecU, dmU, ierr))
   PetscCallA(DMGetLocalSection(dmU, sectionU, ierr))
   PetscCall(project(locVecU, sectionU, ierr))

   do QuadOrder = QuadOrderMax, QuadOrderMax
      do k = 0, (dim - 2) * QuadOrderMax
         do j = 0, QuadOrderMax
            do i = 0, QuadOrderMax
               !!! Initialize a field
               if (i + j + k <= QuadOrder) then
                  !!! Integrate
                  if (dim == 2) then
                     sol = ((xr**(1.0_kr + i) - xl**(1.0_kr + i)) / (1.0_kr + i)) * ((yr**(1.0_kr + j) - yl**(1.0_kr + j)) / (1.0_kr + j))
                     PetscCallA(Integrate2D_Scal(MEF90Ctx, locVecU, i, j, QuadOrder, Scal, v2d, ierr))
                  else
                     sol = ((xr**(1.0_kr + i) - xl**(1.0_kr + i)) / (1.0_kr + i)) * ((yr**(1.0_kr + j) - yl**(1.0_kr + j)) / (1.0_kr + j)) * ((zr**(1.0_kr + k) - zl**(1.0_kr + k)) / (1.0_kr + k))
                     PetscCallA(Integrate3D_Scal(MEF90Ctx, locVecU, i, j, k, QuadOrder, Scal, v3d, ierr))
                  end if
                  write (IOBuffer, 100) i, j, k, QuadOrder, Scal, sol, abs(Scal - sol), abs(Scal - sol) / sol
                  PetscCallA(PetscPrintf(PETSC_COMM_WORLD, IOBuffer, ierr))
               end if
            end do
         end do
      end do
   end do

100 format('Integrating X^', I1, ' * Y^', I1, ' * Z^', I1, ' at order ', I4, ' : ', 2es12.5, ': absolute error ', ES12.5, ': relative error ', ES12.5, "\n")

   PetscCallA(VecDestroy(locVecU, ierr))
   PetscCall(DMDestroy(dm, ierr))
   PetscCall(MEF90CtxDestroy(MEF90Ctx, ierr))
   PetscCall(MEF90Finalize(ierr))
   PetscCall(PetscFinalize(ierr))

contains
   subroutine project(v, s, ierr)
      type(tVec), intent(IN)              :: v
      type(tPetscSection), intent(IN)     :: s
      PetscErrorCode, intent(INOUT)       :: ierr

      PetscInt                           :: pStart, pEnd, p, numDof, i
      type(tDM)                          :: dm
      type(tPetscSection)                :: coordSection
      type(tVec)                         :: coordVec
      PetscScalar, dimension(:), pointer   :: coordArray, vArray
      PetscScalar, dimension(3)           :: xyz
      PetscInt                           :: dim, pOffset

      PetscCallA(PetscSectionGetChart(s, pStart, pEnd, ierr))
      PetscCallA(VecGetDM(v, dm, ierr))
      PetscCallA(DMGetCoordinateSection(dm, coordSection, ierr))
      PetscCallA(DMGetCoordinatesLocal(dm, coordVec, ierr))
      PetscCallA(DMGetDimension(dm, dim, ierr))
      PetscCallA(VecGetArray(v, vArray, ierr))

      do p = pStart, pEnd - 1
         PetscCallA(PetscSectionGetDof(s, p, numDof, ierr))
         if (numDof > 0) then
            !!! trick: the coordinate of a point is the average of the coordinates of the points in its closure
            PetscCallA(DMPlexVecGetClosure(dm, coordSection, coordVec, p, coordArray, ierr))
            do i = 1, dim
               xyz(i) = sum(coordArray(i:size(coordArray):dim)) * dim / size(coordArray)
            end do
            PetscCallA(DMPlexVecRestoreClosure(dm, coordSection, coordVec, p, coordArray, ierr))

            PetscCallA(PetscSectionGetOffset(s, p, pOffset, ierr))
            do i = 1, numDof
               vArray(pOffset + i) = xyz(i)
            end do
         end if
      end do
      PetscCallA(VecRestoreArray(v, vArray, ierr))
      !!! Of course, this does not use informations from the section, so it does over-write constrained values
   end subroutine project

   subroutine Integrate3D_Scal(MEF90Ctx, v, i, j, k, QuadratureOrder, i1, i2, ierr)
      type(MEF90Ctx_Type), intent(IN)                     :: MEF90Ctx
      type(tVec), intent(IN)                              :: v
      PetscInt, intent(IN)                                :: QuadratureOrder, i, j, k
      PetscReal, intent(OUT)                              :: i1
      type(Vect3D), intent(OUT)                           :: i2
      PetscErrorCode, intent(OUT)                         :: ierr

      type(tDM)                                          :: dm
      type(MEF90CtxGlobalOptions_Type), pointer           :: MEF90CtxGlobalOptions
      type(tIS)                                          :: setIS, setPointIS
      type(MEF90_ELEMENT_SCAL), dimension(:), pointer      :: elem
      type(MEF90ElementType)                             :: elementType
      DMPolytopeType                                     :: cellType
      PetscInt                                           :: iDof, iGauss, cell, set, dim
      PetscInt, dimension(:), pointer                      :: setID, setPointID
      PetscReal                                          :: X, Y, Z
      PetscReal, dimension(:), pointer                     :: coordDof

      i1 = 0.0_kr
      i2 = 0.0_kr

      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag, MEF90CtxGlobalOptions, ierr))
      PetscCall(VecGetDM(v, dm, ierr))
      PetscCall(DMGetDimension(dm, dim, ierr))
      PetscCall(DMGetLabelIdIS(dm, MEF90CellSetLabelName, setIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90Ctx%comm, setIS, ierr))
      if (setIS /= PETSC_NULL_IS) then
         PetscCall(ISGetIndices(setIS, setID, ierr))
         do set = 1, size(setID)
            PetscCall(DMGetStratumIS(dm, MEF90CellSetLabelName, setID(set), setPointIS, ierr))
            if (setPointIS /= PETSC_NULL_IS) then

               PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
               PetscCall(DMPlexGetCellType(dm, setPointID(1), cellType, ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellType, elementType, ierr))
               PetscCall(MEF90ElementCreate(dm, setPointIS, elem, QuadratureOrder, elementType, ierr))

               do cell = 1, size(setPointID)
                  PetscCall(DMPlexVecGetClosure(dm, PETSC_NULL_SECTION, v, setPointID(cell), coordDof, ierr))
                  do iGauss = 1, size(elem(cell)%Gauss_C)
                     X = 0.0_kr
                     Y = 0.0_kr
                     Z = 0.0_kr
                     do iDof = 1, size(elem(cell)%BF(:, 1))
                        X = X + Elem(cell)%BF(iDoF, iGauss) * coordDof(3 * (iDof - 1) + 1)
                        Y = Y + Elem(cell)%BF(iDoF, iGauss) * coordDof(3 * (iDof - 1) + 2)
                        Z = Z + Elem(cell)%BF(iDoF, iGauss) * coordDof(3 * (iDof - 1) + 3)
                     end do
                     i1 = i1 + X**i * Y**j * Z**k * elem(cell)%Gauss_C(iGauss)
                  end do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dm, PETSC_NULL_SECTION, v, setPointID(cell), coordDof, ierr))
               end do ! cell

               PetscCall(MEF90ElementDestroy(elem, ierr))
               PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
            end if ! pointIS
            PetscCall(ISDestroy(setPointIS, ierr))
         end do ! set
         PetscCall(ISRestoreIndices(setIS, setID, ierr))
      end if ! setIS
      PetscCall(ISDestroy(setIS, ierr))
   end subroutine Integrate3D_Scal

   subroutine Integrate2D_Scal(MEF90Ctx, v, i, j, QuadratureOrder, i1, i2, ierr)
      type(MEF90Ctx_Type), intent(IN)                     :: MEF90Ctx
      type(tVec), intent(IN)                              :: v
      PetscInt, intent(IN)                                :: QuadratureOrder, i, j
      PetscReal, intent(OUT)                              :: i1
      type(Vect2D), intent(OUT)                           :: i2
      PetscErrorCode, intent(OUT)                         :: ierr

      type(tDM)                                          :: dm
      type(MEF90CtxGlobalOptions_Type), pointer           :: MEF90CtxGlobalOptions
      type(tIS)                                          :: setIS, setPointIS
      type(MEF90_ELEMENT_SCAL), dimension(:), pointer      :: elem
      type(MEF90ElementType)                             :: elementType
      DMPolytopeType                                     :: cellType
      PetscInt                                           :: iDof, iGauss, cell, set, dim
      PetscInt, dimension(:), pointer                      :: setID, setPointID
      PetscReal                                          :: X, Y
      PetscReal, dimension(:), pointer                     :: coordDof

      i1 = 0.0_kr
      i2 = 0.0_kr

      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag, MEF90CtxGlobalOptions, ierr))
      PetscCall(VecGetDM(v, dm, ierr))
      PetscCall(DMGetDimension(dm, dim, ierr))
      PetscCall(DMGetLabelIdIS(dm, MEF90CellSetLabelName, setIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90Ctx%comm, setIS, ierr))
      if (setIS /= PETSC_NULL_IS) then
         PetscCall(ISGetIndices(setIS, setID, ierr))
         do set = 1, size(setID)
            PetscCall(DMGetStratumIS(dm, MEF90CellSetLabelName, setID(set), setPointIS, ierr))
            if (setPointIS /= PETSC_NULL_IS) then

               PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
               PetscCall(DMPlexGetCellType(dm, setPointID(1), cellType, ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellType, elementType, ierr))
               PetscCall(MEF90ElementCreate(dm, setPointIS, elem, QuadratureOrder, elementType, ierr))

               do cell = 1, size(setPointID)
                  PetscCall(DMPlexVecGetClosure(dm, PETSC_NULL_SECTION, v, setPointID(cell), coordDof, ierr))
                  do iGauss = 1, size(elem(cell)%Gauss_C)
                     X = 0.0_kr
                     Y = 0.0_kr
                     do iDof = 1, size(elem(cell)%BF(:, 1))
                        X = X + Elem(cell)%BF(iDoF, iGauss) * coordDof(2 * (iDof - 1) + 1)
                        Y = Y + Elem(cell)%BF(iDoF, iGauss) * coordDof(2 * (iDof - 1) + 2)
                     end do
                     i1 = i1 + X**i * Y**j * elem(cell)%Gauss_C(iGauss)
                  end do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dm, PETSC_NULL_SECTION, v, setPointID(cell), coordDof, ierr))
               end do ! cell

               PetscCall(MEF90ElementDestroy(elem, ierr))
               PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
            end if ! pointIS
            PetscCall(ISDestroy(setPointIS, ierr))
         end do ! set
         PetscCall(ISRestoreIndices(setIS, setID, ierr))
      end if ! setIS
      PetscCall(ISDestroy(setIS, ierr))
   end subroutine Integrate2D_Scal
end program TestQuadrature
