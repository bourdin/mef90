program TestOrientation
#include <petsc/finclude/petsc.h>
   use m_MEF90
   use petsc
   implicit none(type, external)

   PetscErrorCode                            :: ierr
   type(MEF90Ctx_Type), target                :: MEF90Ctx
   type(MEF90CtxGlobalOptions_Type)          :: MEF90GlobalOptions_default
   type(tDM)                                 :: dm
   PetscBool                                 :: interpolate = PETSC_TRUE
   character(len=MEF90MXSTRLEN)              :: IOBuffer

   PetscInt                                  :: dim, pStart, pEnd
   PetscInt                                  :: cStart, cEnd, c, nFlip = 0
   PetscReal, dimension(:), pointer            :: v0, BB, BBinv
   PetscReal                                 :: detBBinv

   PetscBool                                 :: reversedCells = PETSC_FALSE
   type(MEF90CtxGlobalOptions_Type), pointer  :: MEF90GlobalOptions

   PetscCallA(PetscInitialize(ierr))
   PetscCallA(MEF90Initialize(PETSC_COMM_WORLD, ierr))

   MEF90GlobalOptions_default%verbose = 0
   MEF90GlobalOptions_default%dryrun = PETSC_FALSE
   MEF90GlobalOptions_default%timeMin = 0.0_kr
   MEF90GlobalOptions_default%timeMax = 1.0_kr
   MEF90GlobalOptions_default%timeNumStep = 11
   MEF90GlobalOptions_default%timeInterpolation = MEF90TimeInterpolation_linear
   MEF90GlobalOptions_default%timeSkip = 0
   MEF90GlobalOptions_default%timeNumCycle = 1
   MEF90GlobalOptions_default%elementFamily = MEF90ElementFamilyLagrange
   MEF90GlobalOptions_default%elementOrder = 1

   PetscCallA(MEF90CtxCreate(PETSC_COMM_WORLD, MEF90Ctx, MEF90GlobalOptions_default, ierr))
   PetscCallA(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag, MEF90GlobalOptions, ierr))

   PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm, MEF90Ctx%geometryfile, PETSC_NULL_CHARACTER, interpolate, dm, ierr))
   PetscCallA(DMPlexDistributeSetDefault(dm, PETSC_FALSE, ierr))
   PetscCallA(DMSetFromOptions(dm, ierr))
   PetscCallA(DMViewFromOptions(dm, PETSC_NULL_OBJECT, "-dm_view", ierr))

   distribute: block
      type(tDM), target                    :: dmDist
      PetscInt                            :: ovlp = 0_ki
      if (MEF90Ctx%NumProcs > 1) then
         PetscCallA(DMPlexDistribute(dm, ovlp, PETSC_NULL_SF, dmDist, ierr))
         PetscCallA(DMDestroy(dm, ierr))
         dm = dmDist
      end if
   end block distribute
   PetscCallA(DMViewFromOptions(dm, PETSC_NULL_OBJECT, "-mef90dm_view", ierr))

   PetscCallA(DMGetDimension(dm, dim, ierr))
   PetscCallA(DMPlexGetChart(dm, pStart, pEnd, ierr))
   PetscCallA(DMPlexGetHeightStratum(dm, 0_ki, cStart, cEnd, ierr))
   allocate (v0(2))
   allocate (BB(4))
   allocate (BBinv(4))
   do c = cStart, cEnd - 1
      PetscCallA(DMPlexComputeCellGeometryAffineFEM(dm, c, v0, BB, BBinv, detBBinv, ierr))
      if (detBBinv <= 0.0_kr) then
         reversedCells = PETSC_TRUE
         nFlip = nFlip + 1
         if (MEF90GlobalOptions%verbose > 0) then
            write (IOBuffer, '("Reversed cell: ",I0," Jacobian : ",ES12.5,"\n")') c, detBBinv
            PetscCallA(PetscPrintf(PETSC_COMM_WORLD, IOBuffer, ierr))
         end if
      end if
   end do
   if (reversedCells) then
      write (IOBuffer, '("WARNING: this mesh has ",I0," cells with negative Jacobian out of a total ",I0," cells. \nRun with -verbose 1 to see a list of flipped cells\n")') nFlip, pEnd - pStart
      PetscCallA(PetscPrintf(PETSC_COMM_WORLD, IOBuffer, ierr))
   else
      write (IOBuffer, *) "No reversed cells\n"
      PetscCallA(PetscPrintf(PETSC_COMM_WORLD, IOBuffer, ierr))
   end if
   PetscCallA(DMDestroy(dm, ierr))
   deallocate (v0)
   deallocate (BB)
   deallocate (BBinv)

   call MEF90CtxDestroy(MEF90Ctx, ierr)
   call MEF90Finalize(ierr)
   call PetscFinalize(ierr)
end program TestOrientation

