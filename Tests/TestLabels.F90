program TestLabels
#include <petsc/finclude/petsc.h>
   use m_MEF90
   use petsc
   implicit none(type, external)

   PetscErrorCode                      :: ierr
   type(MEF90Ctx_Type), target          :: MEF90Ctx
   type(MEF90CtxGlobalOptions_Type)    :: MEF90GlobalOptions_default
   type(MEF90CtxGlobalOptions_Type)    :: MEF90GlobalOptions
   type(tDM), target                    :: dm
   PetscBool                           :: interpolate = PETSC_TRUE

   PetscInt                            :: set
   type(tIS)                           :: setIS, setPointIS
   PetscInt, dimension(:), pointer       :: setID, pointID
   PetscInt                            :: labelSize, i
   character(len=MEF90MXSTRLEN)        :: IOBuffer

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

   PetscCallA(PetscInitialize(ierr))

   call MEF90Initialize(PETSC_COMM_WORLD, ierr)
   call MEF90CtxCreate(PETSC_COMM_WORLD, MEF90Ctx, "", ierr)
   MEF90GlobalOptions = MEF90GlobalOptions_default
   PetscCallA(MEF90CtxGlobalOptionsSetFromOptions(MEF90Ctx%comm, trim(MEF90Ctx%prefix), MEF90GlobalOptions, ierr))
   call MEF90Ctx%setFromOptions(ierr); CHKERRQ(ierr)

   PetscCallA(PetscPrintf(MEF90Ctx%Comm, MEF90Ctx%geometryfile//'\n', ierr))
   PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm, MEF90Ctx%geometryfile, PETSC_NULL_CHARACTER, interpolate, dm, ierr))
   PetscCallA(DMPlexDistributeSetDefault(dm, PETSC_FALSE, ierr))
   PetscCallA(DMSetFromOptions(dm, ierr))
   PetscCallA(DMViewFromOptions(dm, PETSC_NULL_OBJECT, "-mef90dm_view", ierr))

   distribute: block
      type(tDM), target                    :: dmDist
      PetscInt                            :: ovlp = 0_ki
      if (MEF90Ctx%NumProcs > 1) then
         PetscCallA(DMPlexDistribute(dm, ovlp, PETSC_NULL_SF, dmDist, ierr))
         PetscCallA(DMDestroy(dm, ierr))
         dm = dmDist
      end if
   end block distribute
   PetscCallA(DMViewFromOptions(dm, PETSC_NULL_OBJECT, "-dm_view", ierr))

   do i = 1, size(MEF90SetLabelName)
      PetscCallA(PetscPrintf(PETSC_COMM_WORLD, '=== '//trim(MEF90SetLabelName(i))//' ===\n', ierr))
      PetscCallA(DMGetLabelSize(dm, MEF90SetLabelName(i), labelSize, ierr))
      write (IOBuffer, '("label size",I5,"\n")') labelSize
      ! Note that label size refers only to the local dm, not the whole mesh
      PetscCallA(PetscSynchronizedPrintf(PETSC_COMM_WORLD, IOBuffer, ierr))
      PetscCallA(PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT, ierr))
      PetscCallA(DMGetLabelIdIS(dm, MEF90SetLabelName(i), setIS, ierr))
      if (setIS /= PETSC_NULL_IS) then
         PetscCallA(ISGetIndices(setIS, setID, ierr))
         write (IOBuffer, *) trim(MEF90SetLabelName(i))//' ID: ', setID, '\n'
         PetscCallA(PetscPrintf(PETSC_COMM_WORLD, IOBuffer, ierr))
         do set = 1, size(setID)
            PetscCallA(DMGetStratumIS(dm, MEF90SetLabelName(i), setID(set), setPointIS, ierr))
            if (setPointIS /= PETSC_NULL_IS) then
               PetscCallA(ISGetIndices(setPointIS, pointID, ierr))
               write (*, *) '   points'
               write (*, *) pointID
               PetscCallA(ISRestoreIndices(setPointIS, pointID, ierr))
            end if ! setPointIS
            PetscCallA(ISDestroy(setPointIS, ierr))
         end do
         PetscCallA(ISRestoreIndices(setIS, setID, ierr))
      end if ! setIS
      PetscCallA(ISDestroy(setIS, ierr))
   end do
   PetscCallA(DMDestroy(dm, ierr))

   call MEF90CtxDestroy(MEF90Ctx, ierr)
   call MEF90Finalize(ierr)
   call PetscFinalize(ierr)
end program TestLabels
