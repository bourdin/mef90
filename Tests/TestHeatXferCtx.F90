program TestHeatXferCtx
#include "petsc/finclude/petsc.h"
   ! Use m_MEF90
   use m_MEF90_HeatXfer
   ! Use petsc
   implicit none(type, external)

   PetscErrorCode                                     :: ierr
   type(MEF90Ctx_Type), target                         :: MEF90Ctx
   type(MEF90CtxGlobalOptions_Type)                    :: MEF90GlobalOptions
   type(MEF90HeatXfer_Type), target                    :: MEF90HeatXferCtx
   type(MEF90HeatXferGlobalOptions_Type)              :: MEF90HeatXferGlobalOptions
   type(tDM)                                          :: dm

   type(tIS)                                          :: cellSetIS, faceSetIS
   PetscInt, dimension(:), pointer                      :: cellSetID, faceSetID
   PetscInt                                           :: step, dim, set, numCellSet, numFaceSet
   PetscReal, dimension(:), pointer                     :: time
   character(len=MEF90MXSTRLEN)                       :: IOBuffer
   PetscBool                                          :: flg
   PetscReal, dimension(:), pointer                     :: energy, bodyWork, surfaceWork

   !!! Initialize MEF90
   PetscCallA(PetscInitialize(ierr))
   PetscCallA(MEF90Initialize(PETSC_COMM_WORLD, ierr))

   !!! Get all MEF90-wide options
   PetscCallA(MEF90CtxCreate(PETSC_COMM_WORLD, MEF90Ctx, "", ierr))
   PetscCallA(MEF90Ctx%setFromOptions(ierr))
   PetscCallA(MEF90CtxGlobalOptionsSetFromOptions(MEF90Ctx%comm, trim(MEF90Ctx%prefix), MEF90GlobalOptions, ierr))
   PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm, MEF90Ctx%geometryFile, PETSC_NULL_CHARACTER, PETSC_TRUE, dm, ierr))
   PetscCallA(DMPlexDistributeSetDefault(dm, PETSC_FALSE, ierr))
   PetscCallA(DMSetUseNatural(dm, PETSC_TRUE, ierr))
   PetscCallA(DMSetFromOptions(dm, ierr))
   PetscCallA(DMViewFromOptions(dm, PETSC_NULL_OBJECT, "-heatXfer_dm_view", ierr))

   PetscCallA(MEF90CtxGetTime(MEF90Ctx, time, ierr))

   inquire (file=MEF90Ctx%resultFile, exist=flg)
   if (flg) then
      ! we assume that the output file exists and is formatted
      PetscCallA(MEF90CtxOpenEXO(MEF90Ctx, MEF90Ctx%resultViewer, FILE_MODE_APPEND, ierr))
   else
      ! we need to create the output file
      EXOFormat: block
         PetscInt                                            :: numNodalVar = 1, numCellVar = 1, numGVar = 0
         character(len=MEF90MXSTRLEN), dimension(:), pointer   :: nodalVarName, cellVarName, gVarName

         allocate (nodalVarName(numNodalVar))
         allocate (cellVarName(numCellVar))
         allocate (gVarName(numGVar))
         nodalVarName = ["Temperature        "]
         cellVarName = ["Flux               "]
         PetscCallA(MEF90CtxOpenEXO(MEF90Ctx, MEF90Ctx%resultViewer, FILE_MODE_WRITE, ierr))
         PetscCallA(MEF90EXODMView(dm, MEF90Ctx%resultViewer, MEF90GlobalOptions%elementOrder, ierr))
         PetscCallA(MEF90EXOFormat(MEF90Ctx%resultViewer, gVarName, cellVarName, nodalVarName, time, ierr))
         deallocate (nodalVarName)
         deallocate (cellVarName)
         deallocate (gVarName)
      end block EXOFormat
   end if
   distribute: block
      type(tDM), target                    :: dmDist
      PetscInt                            :: ovlp = 0_ki
      type(tPetscSF)                      :: naturalPointSF

      if (MEF90Ctx%NumProcs > 1) then
         PetscCallA(DMPlexDistribute(dm, ovlp, naturalPointSF, dmDist, ierr))
         PetscCallA(DMPlexSetMigrationSF(dmDist, naturalPointSF, ierr))
         PetscCallA(PetscSFDestroy(naturalPointSF, ierr))
         PetscCallA(DMDestroy(dm, ierr))
         dm = dmDist
      end if
   end block distribute
   PetscCallA(DMViewFromOptions(dm, PETSC_NULL_OBJECT, "-heatXfer_dm_view", ierr))

   !!! Create HeatXfer context, get all HeatXfer options
   PetscCallA(MEF90HeatXferCreate(MEF90HeatXferCtx, dm, MEF90Ctx, "", ierr))
   PetscCallA(MEF90HeatXferCtx%setFromOptions(ierr))
   !!! We no longer need the DM. We have the megaDM in MEF90HeatXferCtx
   PetscCallA(DMDestroy(dm, ierr))
   MEF90HeatXferGlobalOptions = MEF90HeatXferCtx%globalOptions

   PetscCallA(DMGetDimension(MEF90HeatXferCtx%megaDM, dim, ierr))

   PetscCallA(DMGetLabelIdIS(MEF90HeatXferCtx%megaDM, MEF90CellSetLabelName, cellSetIS, ierr))
   PetscCallA(MEF90ISAllGatherMerge(PETSC_COMM_WORLD, cellSetIS, ierr))
   PetscCallA(IsGetSize(cellSetIS, numCellSet, ierr))
   PetscCallA(ISGetIndices(cellSetIS, cellSetID, ierr))

   PetscCallA(DMGetLabelIdIS(MEF90HeatXferCtx%megaDM, MEF90FaceSetLabelName, faceSetIS, ierr))
   PetscCallA(MEF90ISAllGatherMerge(PETSC_COMM_WORLD, faceSetIS, ierr))
   PetscCallA(IsGetSize(faceSetIS, numFaceSet, ierr))
   PetscCallA(ISGetIndices(faceSetIS, faceSetID, ierr))

   allocate (energy(numCellSet))
   allocate (bodyWork(numCellSet))
   allocate (surfaceWork(numFaceSet))

   !!! Analysis loop:
   do step = 1, size(time)
      write (IOBuffer, '("Step: ",I4," Analysis time: ",ES12.5,"\n")') step, time(step)
      PetscCallA(PetscPrintf(PETSC_COMM_WORLD, IOBuffer, ierr))
      PetscCallA(VecSet(MEF90HeatXferCtx%temperatureLocal, time(step), ierr))
      PetscCallA(MEF90HeatXferSetTransients(MEF90HeatXferCtx, step, time(step), ierr))
      PetscCallA(MEF90HeatXFerEnergy(MEF90HeatXferCtx, energy, bodyWork, surfaceWork, ierr))
      do set = 1, numCellSet
         write (*, '("   Cell set ",I4," energy:       ",ES12.5," body work: ",ES12.5)') cellSetID(set), energy(set), bodyWork(set)
      end do
      do set = 1, numFaceSet
         write (*, '("   Face set ",I4," surface work: ",ES12.5)') faceSetID(set), surfaceWork(set)
      end do
      PetscCallA(MEF90HeatXferViewEXO(MEF90HeatXferCtx, step, ierr))
   end do ! step
   deallocate (energy)
   deallocate (bodyWork)
   deallocate (surfaceWork)
   PetscCallA(ISGetIndices(faceSetIS, faceSetID, ierr))
   PetscCallA(ISDestroy(faceSetIS, ierr))
   PetscCallA(ISGetIndices(cellSetIS, cellSetID, ierr))
   PetscCallA(ISDestroy(cellSetIS, ierr))

   PetscCallA(MEF90HeatXferDestroy(MEF90HeatXferCtx, ierr))
   PetscCallA(MEF90CtxDestroy(MEF90Ctx, ierr))
   PetscCallA(MEF90Finalize(ierr))
   PetscCallA(PetscFinalize(ierr))
end program TestHeatXferCtx
