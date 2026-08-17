#include "../MEF90/mef90.inc"
module MEF90_APPEND(m_MEF90_HeatXferAssembly,MEF90_DIM)D
#include "petsc/finclude/petsc.h"
use petscsnes
use petsctao
use m_MEF90_HeatXfer_class
implicit none(type, external)

private
public MEF90HeatXferOperator, &
   MEF90HeatXferBilinearForm, &
   MEF90HeatXFerEnergy, &
   MEF90HeatXFerIFunction, &
   MEF90HeatXferIJacobian

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferOperator"
!!! author: Blaise Bourdin (2012-14, bourdin@lsu.edu)
!!! author: Blaise Bourdin (2022, bourdin@mcmaster.ca)
!!! author: Alexis Marboeuf (2022, marboeua@mcmaster.ca)
!!!
!!!  MEF90HeatXferOperator: Build the operator. When called in SNES, the solution time should always match the target time,
!!!                         so there is no need for interpolation of the fluxes, external, and boundary values
!!!

subroutine MEF90HeatXferOperator(snesTemp, x, residual, MEF90HeatXferCtx, ierr)
   type(tSNES), intent(IN)                         :: snesTemp
   type(tVec), intent(IN)                          :: x
   type(tVec), intent(INOUT)                       :: residual
   type(MEF90HeatXfer_Type), intent(IN)            :: MEF90HeatXferCtx
   PetscErrorCode, intent(INOUT)                   :: ierr

   type(tDM)                                       :: dmTemperature, dmFlux, dmBoundaryFlux, dmExternalTemperature
   type(tPetscSection)                             :: sectionFlux, sectionBoundaryFlux, sectionExternalTemperature
   type(tVec)                                      :: locResidual
   type(tIS)                                       :: setIS, setPointIS
   PetscInt, dimension(:), pointer                 :: setID, setPointID
   PetscInt                                        :: set, QuadratureOrder, vecOffset
   PetscInt                                        :: cell, iDof, jDof, iGauss
   type(MEF90HeatXferCellSetOptions_Type)          :: cellSetOptions
   type(MEF90_MATS)          :: thermalConductivity
   type(MEF90HeatXferFaceSetOptions_Type)          :: faceSetOptions
   character(len=MEF90MXSTRLEN)                    :: setPrefix
   type(MEF90_ELEMENT_SCAL), dimension(:), pointer :: elem
   type(MEF90ElementType)                          :: elementType
   DMPolytopeType                                  :: cellType
   type(MEF90CtxGlobalOptions_Type)                :: MEF90CtxGlobalOptions
   type(MEF90HeatXferGlobalOptions_Type)           :: MEF90HeatXferGlobalOptions
   PetscReal, dimension(:), pointer                :: temperatureDof, fluxArray, boundaryFluxArray, externalTemperatureArray, residualDof
   type(MEF90_VECT)                                :: advectionVec

   PetscCall(MEF90CtxGlobalOptionsSetFromOptions(MEF90HeatXferCtx%MEF90Ctx%comm, trim(MEF90HeatXferCtx%MEF90Ctx%prefix), MEF90CtxGlobalOptions, ierr))
   MEF90HeatXferGlobalOptions = MEF90HeatXferCtx%globalOptions
   PetscCall(VecGetDM(MEF90HeatXferCtx%temperatureLocal, dmTemperature, ierr))
   PetscCall(VecGetDM(MEF90HeatXferCtx%fluxLocal, dmFlux, ierr))
   PetscCall(DMGetLocalSection(dmFlux, sectionFlux, ierr))
   PetscCall(VecGetDM(MEF90HeatXferCtx%boundaryFluxLocal, dmBoundaryFlux, ierr))
   PetscCall(DMGetLocalSection(dmBoundaryFlux, sectionBoundaryFlux, ierr))
   PetscCall(VecGetDM(MEF90HeatXferCtx%externalTemperatureLocal, dmExternalTemperature, ierr))
   PetscCall(DMGetLocalSection(dmExternalTemperature, sectionExternalTemperature, ierr))

   PetscCall(DMGetLocalVector(dmTemperature, locResidual, ierr))
   PetscCall(DMGlobalToLocal(dmTemperature, x, INSERT_VALUES, MEF90HeatXferCtx%TemperatureLocal, ierr))

   PetscCall(VecSet(residual, 0.0_kr, ierr))
   PetscCall(VecSet(locResidual, 0.0_kr, ierr))

      !! cell-based contributions
   PetscCall(DMGetLabelIdIS(dmTemperature, MEF90CellSetLabelName, setIS, ierr))
   PetscCall(MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm, setIS, ierr))
   if (.not. PetscObjectIsNull(setIS)) then
      PetscCall(ISGetIndices(setIS, setID, ierr))
      PetscCall(VecGetArray(MEF90HeatXferCtx%fluxLocal, fluxArray, ierr))
      do set = 1, size(setID)
         PetscCall(DMGetStratumIS(dmTemperature, MEF90CellSetLabelName, setID(set), setPointIS, ierr))
         if (.not. PetscObjectIsNull(setPointIS)) then
            write (setPrefix, '(A,"cs",I4.4,"_")') trim(MEF90HeatXferCtx%prefix), setID(set)
            PetscCall(MEF90HeatXferCellSetOptionsSetFromOptions(MEF90HeatXferCtx%comm, setPrefix, MEF90_DIM, cellSetOptions, ierr))
            select type (thermalConductivityMat => cellSetOptions%thermalConductivity)
            type is (MEF90_MATS)
               thermalConductivity = thermalConductivityMat
            end select

            PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
            PetscCall(DMPlexGetCellType(dmTemperature, setPointID(1), cellType, ierr))
            PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellType, elementType, ierr))
            QuadratureOrder = elementType%order * 2
            PetscCall(MEF90ElementCreate(dmTemperature, setPointIS, elem, QuadratureOrder, elementType, ierr))

            allocate (residualDof(size(elem(1)%BF(:, 1))))

            do cell = 1, size(setPointID)
               residualDof = 0.0_kr
                  !! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                  !! If this happens, we will need to protect this loop
               PetscCall(DMPlexVecGetClosure(dmTemperature, PETSC_NULL_SECTION, MEF90HeatXferCtx%TemperatureLocal, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
               do iGauss = 1, size(elem(cell)%Gauss_C)
                  do jDof = 1, size(elem(cell)%BF(:, 1))
                     do iDof = 1, size(elem(cell)%BF(:, 1))
                        residualDof(jDof) = residualDof(jDof) + (thermalConductivity * temperatureDof(iDof) * elem(cell)%Grad_BF(iDof, iGauss) .dotP.elem(cell)%Grad_BF(jDof, iGauss)) * elem(cell)%Gauss_C(iGauss)
                     end do ! iDof
                  end do ! jDof
               end do ! iGauss
               PetscCall(DMPlexVecRestoreClosure(dmTemperature, PETSC_NULL_SECTION, MEF90HeatXferCtx%TemperatureLocal, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
               PetscCall(DMPlexVecSetClosure(dmTemperature, PETSC_NULL_SECTION, locResidual, setPointID(cell), residualDof, ADD_VALUES, ierr))
            end do ! cell

            if (norm2(cellSetOptions%advectionVector) /= 0.0_kr) then
               advectionVec = cellSetOptions%advectionVector
               do cell = 1, size(setPointID)
                  residualDof = 0.0_kr
                     !! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !! If this happens, we will need to protect this loop
                  PetscCall(DMPlexVecGetClosure(dmTemperature, PETSC_NULL_SECTION, MEF90HeatXferCtx%TemperatureLocal, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
                  do iGauss = 1, size(elem(cell)%Gauss_C)
                     do jDof = 1, size(elem(cell)%BF(:, 1))
                        do iDof = 1, size(elem(cell)%BF(:, 1))
                           residualDof(jDof) = residualDof(jDof) - (cellSetOptions%density * cellSetOptions%specificHeat * advectionVec.dotP.temperatureDof(iDof) * elem(cell)%Grad_BF(iDof, iGauss)) * elem(cell)%BF(jDof, iGauss) * elem(cell)%Gauss_C(iGauss)
                        end do ! iDof
                     end do ! jDof
                  end do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmTemperature, PETSC_NULL_SECTION, MEF90HeatXferCtx%TemperatureLocal, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
                  PetscCall(DMPlexVecSetClosure(dmTemperature, PETSC_NULL_SECTION, locResidual, setPointID(cell), residualDof, ADD_VALUES, ierr))
               end do ! cell
            end if ! cellSetOptions%advectionVector

            if (cellSetOptions%flux /= 0.0_kr) then
               do cell = 1, size(setPointID)
                  residualDof = 0.0_kr
                     !! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !! If this happens, we will need to protect this loop
                  PetscCall(PetscSectionGetOffset(sectionFlux, setPointID(cell), vecOffset, ierr))
                  do iGauss = 1, size(elem(cell)%Gauss_C)
                     do iDof = 1, size(elem(cell)%BF(:, 1))
                        residualDof(iDof) = residualDof(iDof) - fluxArray(vecOffset + 1) * elem(cell)%BF(iDof, iGauss) * elem(cell)%Gauss_C(iGauss)
                     end do ! iDof
                  end do ! iGauss
                  PetscCall(DMPlexVecSetClosure(dmTemperature, PETSC_NULL_SECTION, locResidual, setPointID(cell), residualDof, ADD_VALUES, ierr))
               end do ! cell
            end if ! cellSetOptions%flux

            deallocate (residualDof)
            PetscCall(MEF90ElementDestroy(elem, ierr))
            PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
            PetscCall(ISDestroy(setPointIS, ierr))
         end if ! pointIS
      end do ! set
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(VecRestoreArray(MEF90HeatXferCtx%fluxLocal, fluxArray, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end if ! setIS

      !! face-based contributions
   PetscCall(DMGetLabelIdIS(dmTemperature, MEF90FaceSetLabelName, setIS, ierr))
   PetscCall(MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm, setIS, ierr))
   if (.not. PetscObjectIsNull(setIS)) then
      PetscCall(ISGetIndices(setIS, setID, ierr))
      PetscCall(VecGetArray(MEF90HeatXferCtx%boundaryFluxLocal, boundaryFluxArray, ierr))
      PetscCall(VecGetArray(MEF90HeatXferCtx%externalTemperatureLocal, externalTemperatureArray, ierr))
      do set = 1, size(setID)
         PetscCall(DMGetStratumIS(dmTemperature, MEF90FaceSetLabelName, setID(set), setPointIS, ierr))
         if (.not. PetscObjectIsNull(setPointIS)) then
            write (setPrefix, '(A,"fs",I4.4,"_")') trim(MEF90HeatXferCtx%prefix), setID(set)
            PetscCall(MEF90HeatXferFaceSetOptionsSetFromOptions(MEF90HeatXferCtx%comm, setPrefix, faceSetOptions, ierr))

            if (faceSetOptions%boundaryFlux /= 0.0_kr) then
               PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
               PetscCall(DMPlexGetCellType(dmTemperature, setPointID(1), cellType, ierr))
               PetscCall(MEF90ElementGetTypeBoundary(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellType, elementType, ierr))
               QuadratureOrder = elementType%order * 2
               PetscCall(MEF90ElementCreate(dmTemperature, setPointIS, elem, QuadratureOrder, elementType, ierr))
               allocate (residualDof(size(elem(1)%BF(:, 1))))
               do cell = 1, size(setPointID)
                  residualDof = 0.0_kr
                     !! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !! If this happens, we will need to protect this loop
                  PetscCall(PetscSectionGetOffset(sectionBoundaryFlux, setPointID(cell), vecOffset, ierr))
                  do iGauss = 1, size(elem(cell)%Gauss_C)
                     do iDof = 1, size(elem(cell)%BF(:, 1))
                        residualDof(iDof) = residualDof(iDof) - boundaryFluxArray(vecOffset + 1) * elem(cell)%BF(iDof, iGauss) * elem(cell)%Gauss_C(iGauss)
                     end do ! iDof
                  end do ! iGauss
                  PetscCall(DMPlexVecSetClosure(dmTemperature, PETSC_NULL_SECTION, locResidual, setPointID(cell), residualDof, ADD_VALUES, ierr))
               end do ! cell
               deallocate (residualDof)
               PetscCall(MEF90ElementDestroy(elem, ierr))
               PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
            end if ! faceSetOptions%boundaryFlux

            if (faceSetOptions%surfaceThermalConductivity /= 0.0_kr) then
               PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
               PetscCall(DMPlexGetCellType(dmTemperature, setPointID(1), cellType, ierr))
               PetscCall(MEF90ElementGetTypeBoundary(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellType, elementType, ierr))
               QuadratureOrder = elementType%order * 2
               PetscCall(MEF90ElementCreate(dmTemperature, setPointIS, elem, QuadratureOrder, elementType, ierr))
               allocate (residualDof(size(elem(1)%BF(:, 1))))
               do cell = 1, size(setPointID)
                  residualDof = 0.0_kr
                     !! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !! If this happens, we will need to protect this loop
                  PetscCall(DMPlexVecGetClosure(dmTemperature, PETSC_NULL_SECTION, MEF90HeatXferCtx%TemperatureLocal, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
                  PetscCall(PetscSectionGetOffset(sectionExternalTemperature, setPointID(cell), vecOffset, ierr))
                  do iGauss = 1, size(elem(cell)%Gauss_C)
                     do jDof = 1, size(elem(cell)%BF(:, 1))
                        do iDof = 1, size(elem(cell)%BF(:, 1))
                           residualDof(jDof) = residualDof(jDof) + faceSetOptions%surfaceThermalConductivity * temperatureDof(iDof) * elem(cell)%BF(iDof, iGauss) * elem(cell)%BF(jDof, iGauss) * elem(cell)%Gauss_C(iGauss)
                        end do ! iDof
                        residualDof(jDof) = residualDof(jDof) - faceSetOptions%surfaceThermalConductivity * externalTemperatureArray(vecOffset + 1) * elem(cell)%BF(jDof, iGauss) * elem(cell)%Gauss_C(iGauss)
                     end do ! jDof
                  end do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmTemperature, PETSC_NULL_SECTION, MEF90HeatXferCtx%TemperatureLocal, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
                  PetscCall(DMPlexVecSetClosure(dmTemperature, PETSC_NULL_SECTION, locResidual, setPointID(cell), residualDof, ADD_VALUES, ierr))
               end do ! cell
               deallocate (residualDof)
               PetscCall(MEF90ElementDestroy(elem, ierr))
               PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
            end if ! faceSetOptions%surfaceThermalConductivity
            PetscCall(ISDestroy(setPointIS, ierr))
         end if ! pointIS
      end do ! set
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(VecRestoreArray(MEF90HeatXferCtx%boundaryFluxLocal, boundaryFluxArray, ierr))
      PetscCall(VecRestoreArray(MEF90HeatXferCtx%externalTemperatureLocal, externalTemperatureArray, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end if ! setIS
   PetscCall(DMLocalToGlobalBegin(dmTemperature, locResidual, ADD_VALUES, residual, ierr))
   PetscCall(DMLocalToGlobalEnd(dmTemperature, locResidual, ADD_VALUES, residual, ierr))
   !PetscCall(DMRestoreLocalVector(dmTemperature,locTemperature,ierr))
   PetscCall(DMRestoreLocalVector(dmTemperature, locResidual, ierr))
end subroutine MEF90HeatXferOperator

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferBilinearForm"
!!! author: Blaise Bourdin (2012-14, bourdin@lsu.edu)
!!! author: Blaise Bourdin (2022, bourdin@mcmaster.ca)
!!! author: Alexis Marboeuf (2022, marboeua@mcmaster.ca)
!!!
!!!  MEF90HeatXferBilinearForm:
!!!

subroutine MEF90HeatXferBilinearForm(snesTemp, x, A, M, MEF90HeatXferCtx, ierr)
   type(tSNES), intent(IN)                          :: snesTemp
   type(tVec), intent(IN)                           :: x
   type(tMat), intent(INOUT)                        :: A, M
   type(MEF90HeatXfer_Type), intent(IN)             :: MEF90HeatXferCtx
   PetscErrorCode, intent(INOUT)                    :: ierr

   type(tDM)                                       :: dmTemperature
   type(tIS)                                       :: setIS, setPointIS
   PetscInt, dimension(:), pointer                   :: setID, setPointID
   PetscInt                                        :: set, QuadratureOrder
   PetscInt                                        :: cell, iDof, jDof, iGauss, nbDof
   type(MEF90HeatXferCellSetOptions_Type)          :: cellSetOptions
   type(MEF90_MATS)          :: thermalConductivity
   type(MEF90HeatXferFaceSetOptions_Type)           :: faceSetOptions
   character(len=MEF90MXSTRLEN)                     :: setPrefix
   type(MEF90_ELEMENT_SCAL), dimension(:), pointer   :: elem
   type(MEF90ElementType)                          :: elementType
   DMPolytopeType                                  :: cellType
   type(MEF90CtxGlobalOptions_Type)                 :: MEF90CtxGlobalOptions
   type(MEF90HeatXferGlobalOptions_Type)            :: MEF90HeatXferGlobalOptions
   PetscReal, dimension(:), pointer                  :: matDof
   type(MEF90_VECT)                                :: advectionVec

   PetscCall(MEF90CtxGlobalOptionsSetFromOptions(MEF90HeatXferCtx%MEF90Ctx%comm, trim(MEF90HeatXferCtx%MEF90Ctx%prefix), MEF90CtxGlobalOptions, ierr))
   MEF90HeatXferGlobalOptions = MEF90HeatXferCtx%globalOptions
   PetscCall(VecGetDM(MEF90HeatXferCtx%temperatureLocal, dmTemperature, ierr))

   PetscCall(MatZeroEntries(A, ierr))

      !! cell-based gradient contributions
   PetscCall(DMGetLabelIdIS(dmTemperature, MEF90CellSetLabelName, setIS, ierr))
   PetscCall(MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm, setIS, ierr))
   if (.not. PetscObjectIsNull(setIS)) then
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         PetscCall(DMGetStratumIS(dmTemperature, MEF90CellSetLabelName, setID(set), setPointIS, ierr))
         if (.not. PetscObjectIsNull(setPointIS)) then
            write (setPrefix, '(A,"cs",I4.4,"_")') trim(MEF90HeatXferCtx%prefix), setID(set)
            PetscCall(MEF90HeatXferCellSetOptionsSetFromOptions(MEF90HeatXferCtx%comm, setPrefix, MEF90_DIM, cellSetOptions, ierr))
            select type (thermalConductivityMat => cellSetOptions%thermalConductivity)
            type is (MEF90_MATS)
               thermalConductivity = thermalConductivityMat
            end select

            PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
            PetscCall(DMPlexGetCellType(dmTemperature, setPointID(1), cellType, ierr))
            PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellType, elementType, ierr))
            QuadratureOrder = elementType%order * 2
            PetscCall(MEF90ElementCreate(dmTemperature, setPointIS, elem, QuadratureOrder, elementType, ierr))

            nbDof = size(elem(1)%BF(:, 1))
            allocate (matDof(nbDof * nbDof))

            do cell = 1, size(setPointID)
               matDof = 0.0_kr
                  !! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                  !! If this happens, we will need to protect this loop
               do iGauss = 1, size(elem(cell)%Gauss_C)
                  do jDof = 0, nbDof - 1
                     do iDof = 1, nbDof
                        matDof(jDof * nbDof + iDof) = matDof(jDof * nbDof + iDof) + (thermalConductivity * elem(cell)%Grad_BF(iDof, iGauss) .dotP.elem(cell)%Grad_BF(jDof + 1, iGauss)) * elem(cell)%Gauss_C(iGauss)
                     end do ! iDof
                  end do ! jDof
               end do ! iGauss
               PetscCall(DMPlexMatSetClosure(dmTemperature, PETSC_NULL_SECTION, PETSC_NULL_SECTION, A, setPointID(cell), matDof, ADD_VALUES, ierr))
            end do ! cell

            if (norm2(cellSetOptions%advectionVector) /= 0.0_kr) then
               advectionVec = cellSetOptions%advectionVector
               do cell = 1, size(setPointID)
                  matDof = 0.0_kr
                     !! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !! If this happens, we will need to protect this loop
                  do iGauss = 1, size(elem(cell)%Gauss_C)
                     do jDof = 0, nbDof - 1
                        do iDof = 1, nbDof
                           matDof(jDof * nbDof + iDof) = matDof(jDof * nbDof + iDof) - (cellSetOptions%density * cellSetOptions%specificHeat * advectionVec.dotP.elem(cell)%Grad_BF(iDof, iGauss)) * elem(cell)%BF(jDof + 1, iGauss) * elem(cell)%Gauss_C(iGauss)
                        end do ! iDof
                     end do ! jDof
                  end do ! iGauss
                  PetscCall(DMPlexMatSetClosure(dmTemperature, PETSC_NULL_SECTION, PETSC_NULL_SECTION, A, setPointID(cell), matDof, ADD_VALUES, ierr))
               end do ! cell
            end if ! cellSetOptions%advectionVector

            deallocate (matDof)
            PetscCall(MEF90ElementDestroy(elem, ierr))
            PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
            PetscCall(ISDestroy(setPointIS, ierr))
         end if ! pointIS
      end do ! set
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end if ! setIS

      !! face-based energies
   PetscCall(DMGetLabelIdIS(dmTemperature, MEF90FaceSetLabelName, setIS, ierr))
   PetscCall(MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm, setIS, ierr))
   if (.not. PetscObjectIsNull(setIS)) then
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         PetscCall(DMGetStratumIS(dmTemperature, MEF90FaceSetLabelName, setID(set), setPointIS, ierr))
         if (.not. PetscObjectIsNull(setPointIS)) then
            write (setPrefix, '(A,"fs",I4.4,"_")') trim(MEF90HeatXferCtx%prefix), setID(set)
            PetscCall(MEF90HeatXferFaceSetOptionsSetFromOptions(MEF90HeatXferCtx%comm, setPrefix, faceSetOptions, ierr))

            if (faceSetOptions%surfaceThermalConductivity /= 0.0_kr) then
               PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
               PetscCall(DMPlexGetCellType(dmTemperature, setPointID(1), cellType, ierr))
               PetscCall(MEF90ElementGetTypeBoundary(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellType, elementType, ierr))
               QuadratureOrder = elementType%order * 2
               PetscCall(MEF90ElementCreate(dmTemperature, setPointIS, elem, QuadratureOrder, elementType, ierr))

               nbDof = size(elem(1)%BF(:, 1))
               allocate (matDof(nbDof * nbDof))

               do cell = 1, size(setPointID)
                  matDof = 0.0_kr
                     !! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !! If this happens, we will need to protect this loop
                  do iGauss = 1, size(elem(cell)%Gauss_C)
                     do jDof = 0, nbDof - 1
                        do iDof = 1, nbDof
                           matDof(jDof * nbDof + iDof) = matDof(jDof * nbDof + iDof) + faceSetOptions%surfaceThermalConductivity * elem(cell)%BF(iDof, iGauss) * elem(cell)%BF(jDof + 1, iGauss) * elem(cell)%Gauss_C(iGauss)
                        end do ! iDof
                     end do ! jDof
                  end do ! iGauss
                  PetscCall(DMPlexMatSetClosure(dmTemperature, PETSC_NULL_SECTION, PETSC_NULL_SECTION, A, setPointID(cell), matDof, ADD_VALUES, ierr))
               end do ! cell

               deallocate (matDof)
               PetscCall(MEF90ElementDestroy(elem, ierr))
               PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
            end if ! faceSetOptions%surfaceThermalConductivity
            PetscCall(ISDestroy(setPointIS, ierr))
         end if ! pointIS
      end do ! set
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end if ! setIS
   PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr))
   PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr))
   PetscCall(MatCopy(A, M, SAME_NONZERO_PATTERN, ierr))
end subroutine MEF90HeatXferBilinearForm

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXFerEnergy"
!!! author: Blaise Bourdin (2012-14, bourdin@lsu.edu)
!!! author: Blaise Bourdin (2022, bourdin@mcmaster.ca)
!!!
!!!  MEF90HeatXFerEnergy:
!!!

subroutine MEF90HeatXFerEnergy(MEF90HeatXferCtx, energy, bodyWork, surfaceWork, ierr)
   type(MEF90HeatXfer_Type), intent(IN)              :: MEF90HeatXferCtx
   PetscReal, dimension(:), pointer                  :: energy, bodyWork, surfaceWork
   PetscErrorCode, intent(INOUT)                    :: ierr

   type(tDM)                                       :: dmTemperature, dmFlux, dmBoundaryFlux
   type(tPetscSection)                             :: sectionFlux, sectionBoundaryFlux
   type(tIS)                                       :: setIS, setPointIS
   PetscInt, dimension(:), pointer                   :: setID, setPointID
   PetscInt                                        :: set, QuadratureOrder, vecOffset
   PetscInt                                        :: cell, iDof, iGauss
   type(MEF90HeatXferCellSetOptions_Type)          :: cellSetOptions
   type(MEF90_MATS)          :: thermalConductivity
   type(MEF90HeatXferFaceSetOptions_Type)           :: faceSetOptions
   character(len=MEF90MXSTRLEN)                     :: setPrefix
   type(MEF90_ELEMENT_SCAL), dimension(:), pointer   :: elem
   type(MEF90ElementType)                          :: elementType
   DMPolytopeType                                  :: cellType
   type(MEF90CtxGlobalOptions_Type)                 :: MEF90CtxGlobalOptions
   type(MEF90HeatXferGlobalOptions_Type)            :: MEF90HeatXferGlobalOptions
   PetscReal, dimension(:), pointer                  :: TemperatureDof, fluxArray, boundaryFluxArray
   type(MEF90_VECT)                                :: gradTemperatureCell
   PetscReal                                       :: bodyWorkCell, surfaceWorkCell
   PetscReal                                       :: myEnergy, myBodyWork, mySurfaceWork

   energy = 0.0_kr
   bodyWork = 0.0_kr
   surfaceWork = 0.0_kr
   PetscCall(MEF90CtxGlobalOptionsSetFromOptions(MEF90HeatXferCtx%MEF90Ctx%comm, trim(MEF90HeatXferCtx%MEF90Ctx%prefix), MEF90CtxGlobalOptions, ierr))
   MEF90HeatXferGlobalOptions = MEF90HeatXferCtx%globalOptions
   PetscCall(VecGetDM(MEF90HeatXferCtx%temperatureLocal, dmTemperature, ierr))

   PetscCall(VecGetDM(MEF90HeatXferCtx%fluxLocal, dmFlux, ierr))
   PetscCall(DMGetLocalSection(dmFlux, sectionFlux, ierr))
   PetscCall(VecGetDM(MEF90HeatXferCtx%boundaryFluxLocal, dmBoundaryFlux, ierr))
   PetscCall(DMGetLocalSection(dmBoundaryFlux, sectionBoundaryFlux, ierr))

      !! cell-based energies
   PetscCall(DMGetLabelIdIS(dmTemperature, MEF90CellSetLabelName, setIS, ierr))
   PetscCall(MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm, setIS, ierr))
   if (.not. PetscObjectIsNull(setIS)) then
      PetscCall(ISGetIndices(setIS, setID, ierr))
      PetscCall(VecGetArray(MEF90HeatXferCtx%fluxLocal, fluxArray, ierr))
      do set = 1, size(setID)
         !! myEnergy and myBodyWork are reduced over the whole communicator below, including over the ranks
         !! owning no point of this set, so they have to be initialized outside of the setPointIS test
         myEnergy = 0.0_kr
         myBodyWork = 0.0_kr
         PetscCall(DMGetStratumIS(dmTemperature, MEF90CellSetLabelName, setID(set), setPointIS, ierr))
         myenergy = 0.0_kr
         myBodyWork = 0.0_kr
         if (.not. PetscObjectIsNull(setPointIS)) then
            write (setPrefix, '(A,"cs",I4.4,"_")') trim(MEF90HeatXferCtx%prefix), setID(set)
            PetscCall(MEF90HeatXferCellSetOptionsSetFromOptions(MEF90HeatXferCtx%comm, setPrefix, MEF90_DIM, cellSetOptions, ierr))
            select type (thermalConductivityMat => cellSetOptions%thermalConductivity)
            type is (MEF90_MATS)
               thermalConductivity = thermalConductivityMat
            end select

            PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
            PetscCall(DMPlexGetCellType(dmTemperature, setPointID(1), cellType, ierr))
            PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellType, elementType, ierr))
            QuadratureOrder = elementType%order * 2
            PetscCall(MEF90ElementCreate(dmTemperature, setPointIS, elem, QuadratureOrder, elementType, ierr))

            do cell = 1, size(setPointID)
                  !! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                  !! If this happens, we will need to protect this loop
               PetscCall(DMPlexVecGetClosure(dmTemperature, PETSC_NULL_SECTION, MEF90HeatXferCtx%temperatureLocal, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
               do iGauss = 1, size(elem(cell)%Gauss_C)
                  gradTemperatureCell = 0.0_kr
                  do iDof = 1, size(elem(cell)%BF(:, 1))
                     gradTemperatureCell = gradTemperatureCell + temperatureDof(iDof) * elem(cell)%Grad_BF(iDof, iGauss)
                  end do ! iDof
                  myEnergy = myEnergy + ((thermalConductivity * gradTemperatureCell) .dotP.gradTemperatureCell) * elem(cell)%Gauss_C(iGauss)
               end do ! iGauss
               PetscCall(DMPlexVecRestoreClosure(dmTemperature, PETSC_NULL_SECTION, MEF90HeatXferCtx%temperatureLocal, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
            end do ! cell
            myEnergy = myEnergy * 0.5_kr

            if (cellSetOptions%flux /= 0.0_kr) then
               do cell = 1, size(setPointID)
                     !! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !! If this happens, we will need to protect this loop
                  PetscCall(DMPlexVecGetClosure(dmTemperature, PETSC_NULL_SECTION, MEF90HeatXferCtx%temperatureLocal, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
                  PetscCall(PetscSectionGetOffset(sectionFlux, setPointID(cell), vecOffset, ierr))
                  do iGauss = 1, size(elem(cell)%Gauss_C)
                     bodyWorkCell = 0.0_kr
                     do iDof = 1, size(elem(cell)%BF(:, 1))
                        bodyWorkCell = bodyWorkCell + fluxArray(vecOffset + 1) * temperatureDof(iDof) * elem(cell)%BF(iDof, iGauss)
                     end do ! iDof
                     myBodyWork = myBodyWork + bodyWorkCell * elem(cell)%Gauss_C(iGauss)
                  end do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmTemperature, PETSC_NULL_SECTION, MEF90HeatXferCtx%temperatureLocal, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
               end do ! cell
            end if ! cellSetOptions%flux
            PetscCall(MEF90ElementDestroy(elem, ierr))
            PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
            PetscCall(ISDestroy(setPointIS, ierr))
         end if ! pointIS
         PetscCallMPI(MPI_AllReduce(myEnergy, energy(set), 1, MPIU_SCALAR, MPI_SUM, MEF90HeatXferCtx%MEF90Ctx%comm, ierr))
         PetscCallMPI(MPI_AllReduce(myBodyWork, bodyWork(set), 1, MPIU_SCALAR, MPI_SUM, MEF90HeatXferCtx%MEF90Ctx%comm, ierr))
      end do ! set
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(VecRestoreArray(MEF90HeatXferCtx%fluxLocal, fluxArray, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end if ! setIS

      !! face-based energies
   PetscCall(DMGetLabelIdIS(dmTemperature, MEF90FaceSetLabelName, setIS, ierr))
   PetscCall(MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm, setIS, ierr))
   if (.not. PetscObjectIsNull(setIS)) then
      PetscCall(ISGetIndices(setIS, setID, ierr))
      PetscCall(VecGetArray(MEF90HeatXferCtx%boundaryFluxLocal, boundaryFluxArray, ierr))
      do set = 1, size(setID)
         !! mySurfaceWork is reduced over the whole communicator below, including over the ranks
         !! owning no point of this set, so it has to be initialized outside of the setPointIS test
         mySurfaceWork = 0.0_kr
         PetscCall(DMGetStratumIS(dmTemperature, MEF90FaceSetLabelName, setID(set), setPointIS, ierr))
         if (.not. PetscObjectIsNull(setPointIS)) then
            write (setPrefix, '(A,"fs",I4.4,"_")') trim(MEF90HeatXferCtx%prefix), setID(set)
            PetscCall(MEF90HeatXferFaceSetOptionsSetFromOptions(MEF90HeatXferCtx%comm, setPrefix, faceSetOptions, ierr))
            if (faceSetOptions%boundaryFlux /= 0.0_kr) then
               PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
               PetscCall(DMPlexGetCellType(dmTemperature, setPointID(1), cellType, ierr))
               PetscCall(MEF90ElementGetTypeBoundary(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellType, elementType, ierr))
               QuadratureOrder = elementType%order * 2
               PetscCall(MEF90ElementCreate(dmTemperature, setPointIS, elem, QuadratureOrder, elementType, ierr))

               do cell = 1, size(setPointID)
                     !! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !! If this happens, we will need to protect this loop
                  PetscCall(DMPlexVecGetClosure(dmTemperature, PETSC_NULL_SECTION, MEF90HeatXferCtx%temperatureLocal, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
                  PetscCall(PetscSectionGetOffset(sectionBoundaryFlux, setPointID(cell), vecOffset, ierr))
                  do iGauss = 1, size(elem(cell)%Gauss_C)
                     surfaceWorkCell = 0.0_kr
                     do iDof = 1, size(elem(cell)%BF(:, 1))
                        surfaceWorkCell = surfaceWorkCell + boundaryFluxArray(vecOffset + 1) * temperatureDof(iDof) * elem(cell)%BF(iDof, iGauss)
                     end do ! iDof
                     mySurfaceWork = mySurfaceWork + surfaceWorkCell * elem(cell)%Gauss_C(iGauss)
                  end do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmTemperature, PETSC_NULL_SECTION, MEF90HeatXferCtx%temperatureLocal, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
               end do ! cell

               PetscCall(MEF90ElementDestroy(elem, ierr))
               PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
            end if ! faceSetOptions%boundaryFlux
            PetscCall(ISDestroy(setPointIS, ierr))
         end if ! pointIS
         PetscCallMPI(MPI_AllReduce(mySurfaceWork, surfaceWork(set), 1, MPIU_SCALAR, MPI_SUM, MEF90HeatXferCtx%MEF90Ctx%comm, ierr))
      end do ! set
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(VecRestoreArray(MEF90HeatXferCtx%boundaryFluxLocal, boundaryFluxArray, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end if ! setIS
end subroutine MEF90HeatXFerEnergy

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXFerIFunction"
!!! author: Blaise Bourdin (2014, bourdin@lsu.edu)
!!! author: Blaise Bourdin (2022, bourdin@mcmaster.ca)
!!! author: Alexis Marboeuf (2022, marboeua@mcmaster.ca)
!!!
!!!  MEF90HeatXFerIFunction:
!!!

subroutine MEF90HeatXFerIFunction(tempTS, time, x, xdot, F, MEF90HeatXferCtx, ierr)
   type(tTS), intent(IN)                            :: tempTS
   PetscReal, intent(IN)                            :: time
   type(tVec), intent(IN)                           :: x, xdot
   type(tVec), intent(INOUT)                        :: F
   type(MEF90HeatXfer_Type), intent(IN)             :: MEF90HeatXferCtx
   PetscErrorCode, intent(INOUT)                    :: ierr

   type(tDM)                                       :: dmTemperature, dmFlux, dmBoundaryFlux, dmExternalTemperature
   type(tPetscSection)                             :: sectionFlux, sectionBoundaryFlux, sectionExternalTemperature
   type(tVec)                                      :: locTemperature, locTemperatureDot, locF
   type(tIS)                                       :: setIS, setPointIS
   PetscInt, dimension(:), pointer                   :: setID, setPointID
   PetscInt                                        :: set, QuadratureOrder, vecOffset
   PetscInt                                        :: cell, iDof, jDof, iGauss
   type(MEF90HeatXferCellSetOptions_Type)          :: cellSetOptions
   type(MEF90_MATS)          :: thermalConductivity
   type(MEF90HeatXferFaceSetOptions_Type)           :: faceSetOptions
   character(len=MEF90MXSTRLEN)                     :: setPrefix
   type(MEF90_ELEMENT_SCAL), dimension(:), pointer   :: elem
   type(MEF90ElementType)                          :: elementType
   DMPolytopeType                                  :: cellType
   type(MEF90CtxGlobalOptions_Type)                 :: MEF90CtxGlobalOptions
   type(MEF90HeatXferGlobalOptions_Type)            :: MEF90HeatXferGlobalOptions
   PetscReal, dimension(:), pointer                  :: temperatureDof, temperatureDotDof, fluxArray, boundaryFluxArray, externalTemperatureArray, residualDof
   type(MEF90_VECT)                                :: advectionVec

   PetscCall(MEF90CtxGlobalOptionsSetFromOptions(MEF90HeatXferCtx%MEF90Ctx%comm, trim(MEF90HeatXferCtx%MEF90Ctx%prefix), MEF90CtxGlobalOptions, ierr))
   MEF90HeatXferGlobalOptions = MEF90HeatXferCtx%globalOptions
   PetscCall(VecGetDM(MEF90HeatXferCtx%temperatureLocal, dmTemperature, ierr))
   PetscCall(VecGetDM(MEF90HeatXferCtx%fluxLocal, dmFlux, ierr))
   PetscCall(DMGetLocalSection(dmFlux, sectionFlux, ierr))
   PetscCall(VecGetDM(MEF90HeatXferCtx%boundaryFluxLocal, dmBoundaryFlux, ierr))
   PetscCall(DMGetLocalSection(dmBoundaryFlux, sectionBoundaryFlux, ierr))
   PetscCall(VecGetDM(MEF90HeatXferCtx%externalTemperatureLocal, dmExternalTemperature, ierr))
   PetscCall(DMGetLocalSection(dmExternalTemperature, sectionExternalTemperature, ierr))

   PetscCall(DMGetLocalVector(dmTemperature, locTemperature, ierr))
   PetscCall(DMGetLocalVector(dmTemperature, locTemperatureDot, ierr))
   PetscCall(DMGetLocalVector(dmTemperature, locF, ierr))
   PetscCall(MEF90VecGlobalToLocalConstraint(x, MEF90HeatXferCtx%temperatureLocal, locTemperature, ierr))
   PetscCall(VecSet(locTemperatureDot, 0.0_kr, ierr))
   PetscCall(DMGlobalToLocalBegin(dmTemperature, xdot, INSERT_VALUES, locTemperatureDot, ierr))
   PetscCall(DMGlobalToLocalEnd(dmTemperature, xdot, INSERT_VALUES, locTemperatureDot, ierr))

   PetscCall(VecSet(F, 0.0_kr, ierr))
   PetscCall(VecSet(locF, 0.0_kr, ierr))

      !! cell-based contributions
   PetscCall(DMGetLabelIdIS(dmTemperature, MEF90CellSetLabelName, setIS, ierr))
   PetscCall(MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm, setIS, ierr))
   if (.not. PetscObjectIsNull(setIS)) then
      PetscCall(ISGetIndices(setIS, setID, ierr))
      PetscCall(VecGetArray(MEF90HeatXferCtx%fluxLocal, fluxArray, ierr))
      do set = 1, size(setID)
         PetscCall(DMGetStratumIS(dmTemperature, MEF90CellSetLabelName, setID(set), setPointIS, ierr))
         if (.not. PetscObjectIsNull(setPointIS)) then
            write (setPrefix, '(A,"cs",I4.4,"_")') trim(MEF90HeatXferCtx%prefix), setID(set)
            PetscCall(MEF90HeatXferCellSetOptionsSetFromOptions(MEF90HeatXferCtx%comm, setPrefix, MEF90_DIM, cellSetOptions, ierr))
            select type (thermalConductivityMat => cellSetOptions%thermalConductivity)
            type is (MEF90_MATS)
               thermalConductivity = thermalConductivityMat
            end select

            PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
            PetscCall(DMPlexGetCellType(dmTemperature, setPointID(1), cellType, ierr))
            PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellType, elementType, ierr))
            QuadratureOrder = elementType%order * 2
            PetscCall(MEF90ElementCreate(dmTemperature, setPointIS, elem, QuadratureOrder, elementType, ierr))

            allocate (residualDof(size(elem(1)%BF(:, 1))))

            do cell = 1, size(setPointID)
               residualDof = 0.0_kr
                  !! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                  !! If this happens, we will need to protect this loop
               PetscCall(DMPlexVecGetClosure(dmTemperature, PETSC_NULL_SECTION, locTemperatureDot, setPointID(cell), PETSC_NULL_INTEGER, temperatureDotDof, ierr))
               do iGauss = 1, size(elem(cell)%Gauss_C)
                  do jDof = 1, size(elem(cell)%BF(:, 1))
                     do iDof = 1, size(elem(cell)%BF(:, 1))
                        residualDof(jDof) = residualDof(jDof) + cellSetOptions%density * cellSetOptions%specificHeat * temperatureDotDof(iDof) * elem(cell)%BF(iDof, iGauss) * elem(cell)%BF(jDof, iGauss) * elem(cell)%Gauss_C(iGauss)
                     end do ! iDof
                  end do ! jDof
               end do ! iGauss
               PetscCall(DMPlexVecRestoreClosure(dmTemperature, PETSC_NULL_SECTION, locTemperatureDot, setPointID(cell), PETSC_NULL_INTEGER, temperatureDotDof, ierr))
               PetscCall(DMPlexVecSetClosure(dmTemperature, PETSC_NULL_SECTION, locF, setPointID(cell), residualDof, ADD_VALUES, ierr))
            end do ! cell

            do cell = 1, size(setPointID)
               residualDof = 0.0_kr
                  !! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                  !! If this happens, we will need to protect this loop
               PetscCall(DMPlexVecGetClosure(dmTemperature, PETSC_NULL_SECTION, locTemperature, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
               do iGauss = 1, size(elem(cell)%Gauss_C)
                  do jDof = 1, size(elem(cell)%BF(:, 1))
                     do iDof = 1, size(elem(cell)%BF(:, 1))
                        residualDof(jDof) = residualDof(jDof) + (thermalConductivity * temperatureDof(iDof) * elem(cell)%Grad_BF(iDof, iGauss) .dotP.elem(cell)%Grad_BF(jDof, iGauss)) * elem(cell)%Gauss_C(iGauss)
                     end do ! iDof
                  end do ! jDof
               end do ! iGauss
               PetscCall(DMPlexVecRestoreClosure(dmTemperature, PETSC_NULL_SECTION, locTemperature, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
               PetscCall(DMPlexVecSetClosure(dmTemperature, PETSC_NULL_SECTION, locF, setPointID(cell), residualDof, ADD_VALUES, ierr))
            end do ! cell

            if (norm2(cellSetOptions%advectionVector) /= 0.0_kr) then
               advectionVec = cellSetOptions%advectionVector
               do cell = 1, size(setPointID)
                  residualDof = 0.0_kr
                     !! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !! If this happens, we will need to protect this loop
                  PetscCall(DMPlexVecGetClosure(dmTemperature, PETSC_NULL_SECTION, locTemperature, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
                  do iGauss = 1, size(elem(cell)%Gauss_C)
                     do jDof = 1, size(elem(cell)%BF(:, 1))
                        do iDof = 1, size(elem(cell)%BF(:, 1))
                           residualDof(jDof) = residualDof(jDof) - (cellSetOptions%density * cellSetOptions%specificHeat * advectionVec.dotP.temperatureDof(iDof) * elem(cell)%Grad_BF(iDof, iGauss)) * elem(cell)%BF(jDof, iGauss) * elem(cell)%Gauss_C(iGauss)
                        end do ! iDof
                     end do ! jDof
                  end do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmTemperature, PETSC_NULL_SECTION, locTemperature, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
                  PetscCall(DMPlexVecSetClosure(dmTemperature, PETSC_NULL_SECTION, locF, setPointID(cell), residualDof, ADD_VALUES, ierr))
               end do ! cell
            end if ! cellSetOptions%advectionVector

            if (cellSetOptions%flux /= 0.0_kr) then
               do cell = 1, size(setPointID)
                  residualDof = 0.0_kr
                     !! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !! If this happens, we will need to protect this loop
                  PetscCall(PetscSectionGetOffset(sectionFlux, setPointID(cell), vecOffset, ierr))
                  do iGauss = 1, size(elem(cell)%Gauss_C)
                     do iDof = 1, size(elem(cell)%BF(:, 1))
                        residualDof(iDof) = residualDof(iDof) - fluxArray(vecOffset + 1) * elem(cell)%BF(iDof, iGauss) * elem(cell)%Gauss_C(iGauss)
                     end do ! iDof
                  end do ! iGauss
                  PetscCall(DMPlexVecSetClosure(dmTemperature, PETSC_NULL_SECTION, locF, setPointID(cell), residualDof, ADD_VALUES, ierr))
               end do ! cell
            end if ! cellSetOptions%flux

            deallocate (residualDof)
            PetscCall(MEF90ElementDestroy(elem, ierr))
            PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
            PetscCall(ISDestroy(setPointIS, ierr))
         end if ! pointIS
      end do ! set
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(VecRestoreArray(MEF90HeatXferCtx%fluxLocal, fluxArray, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end if ! setIS

      !! face-based contributions
   PetscCall(DMGetLabelIdIS(dmTemperature, MEF90FaceSetLabelName, setIS, ierr))
   PetscCall(MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm, setIS, ierr))
   if (.not. PetscObjectIsNull(setIS)) then
      PetscCall(ISGetIndices(setIS, setID, ierr))
      PetscCall(VecGetArray(MEF90HeatXferCtx%boundaryFluxLocal, boundaryFluxArray, ierr))
      PetscCall(VecGetArray(MEF90HeatXferCtx%externalTemperatureLocal, externalTemperatureArray, ierr))
      do set = 1, size(setID)
         PetscCall(DMGetStratumIS(dmTemperature, MEF90FaceSetLabelName, setID(set), setPointIS, ierr))
         if (.not. PetscObjectIsNull(setPointIS)) then
            write (setPrefix, '(A,"fs",I4.4,"_")') trim(MEF90HeatXferCtx%prefix), setID(set)
            PetscCall(MEF90HeatXferFaceSetOptionsSetFromOptions(MEF90HeatXferCtx%comm, setPrefix, faceSetOptions, ierr))

            if (faceSetOptions%boundaryFlux /= 0.0_kr) then
               PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
               PetscCall(DMPlexGetCellType(dmTemperature, setPointID(1), cellType, ierr))
               PetscCall(MEF90ElementGetTypeBoundary(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellType, elementType, ierr))
               QuadratureOrder = elementType%order * 2
               PetscCall(MEF90ElementCreate(dmTemperature, setPointIS, elem, QuadratureOrder, elementType, ierr))
               allocate (residualDof(size(elem(1)%BF(:, 1))))
               do cell = 1, size(setPointID)
                  residualDof = 0.0_kr
                     !! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !! If this happens, we will need to protect this loop
                  PetscCall(PetscSectionGetOffset(sectionBoundaryFlux, setPointID(cell), vecOffset, ierr))
                  do iGauss = 1, size(elem(cell)%Gauss_C)
                     do iDof = 1, size(elem(cell)%BF(:, 1))
                        residualDof(iDof) = residualDof(iDof) - boundaryFluxArray(vecOffset + 1) * elem(cell)%BF(iDof, iGauss) * elem(cell)%Gauss_C(iGauss)
                     end do ! iDof
                  end do ! iGauss
                  PetscCall(DMPlexVecSetClosure(dmTemperature, PETSC_NULL_SECTION, locF, setPointID(cell), residualDof, ADD_VALUES, ierr))
               end do ! cell
               deallocate (residualDof)
               PetscCall(MEF90ElementDestroy(elem, ierr))
               PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
            end if ! faceSetOptions%boundaryFlux

            if (faceSetOptions%surfaceThermalConductivity /= 0.0_kr) then
               PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
               PetscCall(DMPlexGetCellType(dmTemperature, setPointID(1), cellType, ierr))
               PetscCall(MEF90ElementGetTypeBoundary(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellType, elementType, ierr))
               QuadratureOrder = elementType%order * 2
               PetscCall(MEF90ElementCreate(dmTemperature, setPointIS, elem, QuadratureOrder, elementType, ierr))
               allocate (residualDof(size(elem(1)%BF(:, 1))))
               do cell = 1, size(setPointID)
                  residualDof = 0.0_kr
                     !! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !! If this happens, we will need to protect this loop
                  PetscCall(DMPlexVecGetClosure(dmTemperature, PETSC_NULL_SECTION, locTemperature, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
                  PetscCall(PetscSectionGetOffset(sectionExternalTemperature, setPointID(cell), vecOffset, ierr))
                  do iGauss = 1, size(elem(cell)%Gauss_C)
                     do jDof = 1, size(elem(cell)%BF(:, 1))
                        do iDof = 1, size(elem(cell)%BF(:, 1))
                           residualDof(jDof) = residualDof(jDof) + faceSetOptions%surfaceThermalConductivity * temperatureDof(iDof) * elem(cell)%BF(iDof, iGauss) * elem(cell)%BF(jDof, iGauss) * elem(cell)%Gauss_C(iGauss)
                        end do ! iDof
                        residualDof(jDof) = residualDof(jDof) - faceSetOptions%surfaceThermalConductivity * externalTemperatureArray(vecOffset + 1) * elem(cell)%BF(jDof, iGauss) * elem(cell)%Gauss_C(iGauss)
                     end do ! jDof
                  end do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dmTemperature, PETSC_NULL_SECTION, locTemperature, setPointID(cell), PETSC_NULL_INTEGER, temperatureDof, ierr))
                  PetscCall(DMPlexVecSetClosure(dmTemperature, PETSC_NULL_SECTION, locF, setPointID(cell), residualDof, ADD_VALUES, ierr))
               end do ! cell
               deallocate (residualDof)
               PetscCall(MEF90ElementDestroy(elem, ierr))
               PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
            end if ! faceSetOptions%surfaceThermalConductivity
            PetscCall(ISDestroy(setPointIS, ierr))
         end if ! pointIS
      end do ! set
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(VecRestoreArray(MEF90HeatXferCtx%boundaryFluxLocal, boundaryFluxArray, ierr))
      PetscCall(VecRestoreArray(MEF90HeatXferCtx%externalTemperatureLocal, externalTemperatureArray, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end if ! setIS
   PetscCall(DMLocalToGlobalBegin(dmTemperature, locF, ADD_VALUES, F, ierr))
   PetscCall(DMLocalToGlobalEnd(dmTemperature, locF, ADD_VALUES, F, ierr))
   PetscCall(DMRestoreLocalVector(dmTemperature, locTemperature, ierr))
   PetscCall(DMRestoreLocalVector(dmTemperature, locTemperatureDot, ierr))
   PetscCall(DMRestoreLocalVector(dmTemperature, locF, ierr))
end subroutine MEF90HeatXFerIFunction

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferIJacobian"
!!! author: Blaise Bourdin (2014, bourdin@lsu.edu)
!!! author: Blaise Bourdin (2022, bourdin@mcmaster.ca)
!!! author: Alexis Marboeuf (2022, marboeua@mcmaster.ca)
!!!
!!!  MEF90HeatXferIJacobian:
!!!

subroutine MEF90HeatXferIJacobian(tempTS, t, x, xdot, shift, A, M, MEF90HeatXferCtx, ierr)
   type(tTS), intent(IN)                            :: tempTS
   PetscReal, intent(IN)                            :: t
   type(tVec), intent(IN)                           :: x, xdot
   PetscReal, intent(IN)                            :: shift
   type(tMat), intent(INOUT)                        :: A, M
   type(MEF90HeatXfer_Type), intent(IN)             :: MEF90HeatXferCtx
   PetscErrorCode, intent(INOUT)                    :: ierr

   type(tDM)                                       :: dmTemperature
   type(tIS)                                       :: setIS, setPointIS
   PetscInt, dimension(:), pointer                   :: setID, setPointID
   PetscInt                                        :: set, QuadratureOrder
   PetscInt                                        :: cell, iDof, jDof, iGauss, nbDof
   type(MEF90HeatXferCellSetOptions_Type)          :: cellSetOptions
   type(MEF90_MATS)          :: thermalConductivity
   type(MEF90HeatXferFaceSetOptions_Type)           :: faceSetOptions
   character(len=MEF90MXSTRLEN)                     :: setPrefix
   type(MEF90_ELEMENT_SCAL), dimension(:), pointer   :: elem
   type(MEF90ElementType)                          :: elementType
   DMPolytopeType                                  :: cellType
   type(MEF90CtxGlobalOptions_Type)                 :: MEF90CtxGlobalOptions
   type(MEF90HeatXferGlobalOptions_Type)            :: MEF90HeatXferGlobalOptions
   PetscReal, dimension(:), pointer                  :: matDof
   type(MEF90_VECT)                                :: advectionVec

   PetscCall(MEF90CtxGlobalOptionsSetFromOptions(MEF90HeatXferCtx%MEF90Ctx%comm, trim(MEF90HeatXferCtx%MEF90Ctx%prefix), MEF90CtxGlobalOptions, ierr))
   MEF90HeatXferGlobalOptions = MEF90HeatXferCtx%globalOptions
   PetscCall(VecGetDM(MEF90HeatXferCtx%temperatureLocal, dmTemperature, ierr))

   PetscCall(MatZeroEntries(A, ierr))

      !! cell-based gradient contributions
   PetscCall(DMGetLabelIdIS(dmTemperature, MEF90CellSetLabelName, setIS, ierr))
   PetscCall(MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm, setIS, ierr))
   if (.not. PetscObjectIsNull(setIS)) then
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         PetscCall(DMGetStratumIS(dmTemperature, MEF90CellSetLabelName, setID(set), setPointIS, ierr))
         if (.not. PetscObjectIsNull(setPointIS)) then
            write (setPrefix, '(A,"cs",I4.4,"_")') trim(MEF90HeatXferCtx%prefix), setID(set)
            PetscCall(MEF90HeatXferCellSetOptionsSetFromOptions(MEF90HeatXferCtx%comm, setPrefix, MEF90_DIM, cellSetOptions, ierr))
            select type (thermalConductivityMat => cellSetOptions%thermalConductivity)
            type is (MEF90_MATS)
               thermalConductivity = thermalConductivityMat
            end select

            PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
            PetscCall(DMPlexGetCellType(dmTemperature, setPointID(1), cellType, ierr))
            PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellType, elementType, ierr))
            QuadratureOrder = elementType%order * 2
            PetscCall(MEF90ElementCreate(dmTemperature, setPointIS, elem, QuadratureOrder, elementType, ierr))

            nbDof = size(elem(1)%BF(:, 1))
            allocate (matDof(nbDof * nbDof))

            do cell = 1, size(setPointID)
               matDof = 0.0_kr
                  !! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                  !! If this happens, we will need to protect this loop
               do iGauss = 1, size(elem(cell)%Gauss_C)
                  do jDof = 0, nbDof - 1
                     do iDof = 1, nbDof
                        matDof(jDof * nbDof + iDof) = matDof(jDof * nbDof + iDof) + shift * cellSetOptions%density * cellSetOptions%specificHeat * elem(cell)%BF(iDof, iGauss) * elem(cell)%BF(jDof + 1, iGauss) * elem(cell)%Gauss_C(iGauss)
                     end do ! iDof
                  end do ! jDof
               end do ! iGauss
               PetscCall(DMPlexMatSetClosure(dmTemperature, PETSC_NULL_SECTION, PETSC_NULL_SECTION, A, setPointID(cell), matDof, ADD_VALUES, ierr))
            end do ! cell

            do cell = 1, size(setPointID)
               matDof = 0.0_kr
                  !! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                  !! If this happens, we will need to protect this loop
               do iGauss = 1, size(elem(cell)%Gauss_C)
                  do jDof = 0, nbDof - 1
                     do iDof = 1, nbDof
                        matDof(jDof * nbDof + iDof) = matDof(jDof * nbDof + iDof) + (thermalConductivity * elem(cell)%Grad_BF(iDof, iGauss) .dotP.elem(cell)%Grad_BF(jDof + 1, iGauss)) * elem(cell)%Gauss_C(iGauss)
                     end do ! iDof
                  end do ! jDof
               end do ! iGauss
               PetscCall(DMPlexMatSetClosure(dmTemperature, PETSC_NULL_SECTION, PETSC_NULL_SECTION, A, setPointID(cell), matDof, ADD_VALUES, ierr))
            end do ! cell

            if (norm2(cellSetOptions%advectionVector) /= 0.0_kr) then
               advectionVec = cellSetOptions%advectionVector
               do cell = 1, size(setPointID)
                  matDof = 0.0_kr
                     !! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !! If this happens, we will need to protect this loop
                  do iGauss = 1, size(elem(cell)%Gauss_C)
                     do jDof = 0, nbDof - 1
                        do iDof = 1, nbDof
                           matDof(jDof * nbDof + iDof) = matDof(jDof * nbDof + iDof) - (cellSetOptions%density * cellSetOptions%specificHeat * advectionVec.dotP.elem(cell)%Grad_BF(iDof, iGauss)) * elem(cell)%BF(jDof + 1, iGauss) * elem(cell)%Gauss_C(iGauss)
                        end do ! iDof
                     end do ! jDof
                  end do ! iGauss
                  PetscCall(DMPlexMatSetClosure(dmTemperature, PETSC_NULL_SECTION, PETSC_NULL_SECTION, A, setPointID(cell), matDof, ADD_VALUES, ierr))
               end do ! cell
            end if ! cellSetOptions%advectionVector

            deallocate (matDof)
            PetscCall(MEF90ElementDestroy(elem, ierr))
            PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
            PetscCall(ISDestroy(setPointIS, ierr))
         end if ! pointIS
      end do ! set
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end if ! setIS

      !! face-based energies
   PetscCall(DMGetLabelIdIS(dmTemperature, MEF90FaceSetLabelName, setIS, ierr))
   PetscCall(MEF90ISAllGatherMerge(MEF90HeatXferCtx%MEF90Ctx%comm, setIS, ierr))
   if (.not. PetscObjectIsNull(setIS)) then
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         PetscCall(DMGetStratumIS(dmTemperature, MEF90FaceSetLabelName, setID(set), setPointIS, ierr))
         if (.not. PetscObjectIsNull(setPointIS)) then
            write (setPrefix, '(A,"fs",I4.4,"_")') trim(MEF90HeatXferCtx%prefix), setID(set)
            PetscCall(MEF90HeatXferFaceSetOptionsSetFromOptions(MEF90HeatXferCtx%comm, setPrefix, faceSetOptions, ierr))

            if (faceSetOptions%surfaceThermalConductivity /= 0.0_kr) then
               PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
               PetscCall(DMPlexGetCellType(dmTemperature, setPointID(1), cellType, ierr))
               PetscCall(MEF90ElementGetTypeBoundary(MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, cellType, elementType, ierr))
               QuadratureOrder = elementType%order * 2
               PetscCall(MEF90ElementCreate(dmTemperature, setPointIS, elem, QuadratureOrder, elementType, ierr))

               nbDof = size(elem(1)%BF(:, 1))
               allocate (matDof(nbDof * nbDof))

               do cell = 1, size(setPointID)
                  matDof = 0.0_kr
                     !! This could break if TemperatureLocal had no dof in any point in the closure of setPointID(set)
                     !! If this happens, we will need to protect this loop
                  do iGauss = 1, size(elem(cell)%Gauss_C)
                     do jDof = 0, nbDof - 1
                        do iDof = 1, nbDof
                           matDof(jDof * nbDof + iDof) = matDof(jDof * nbDof + iDof) + faceSetOptions%surfaceThermalConductivity * elem(cell)%BF(iDof, iGauss) * elem(cell)%BF(jDof + 1, iGauss) * elem(cell)%Gauss_C(iGauss)
                        end do ! iDof
                     end do ! jDof
                  end do ! iGauss
                  PetscCall(DMPlexMatSetClosure(dmTemperature, PETSC_NULL_SECTION, PETSC_NULL_SECTION, A, setPointID(cell), matDof, ADD_VALUES, ierr))
               end do ! cell

               deallocate (matDof)
               PetscCall(MEF90ElementDestroy(elem, ierr))
               PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
            end if ! faceSetOptions%surfaceThermalConductivity
            PetscCall(ISDestroy(setPointIS, ierr))
         end if ! pointIS
      end do ! set
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end if ! setIS
   PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr))
   PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr))
   PetscCall(MatCopy(A, M, SAME_NONZERO_PATTERN, ierr))
end subroutine MEF90HeatXferIJacobian
end module MEF90_APPEND(m_MEF90_HeatXferAssembly,MEF90_DIM)D
