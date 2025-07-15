program TestMEF90Ctx
#include <petsc/finclude/petsc.h>
   use m_MEF90
   use petsc
   implicit none(type, external)

   PetscErrorCode                      :: ierr
   type(MEF90Ctx_Type), target          :: MEF90Ctx
   type(MEF90CtxGlobalOptions_Type)    :: MEF90GlobalOptions_default
   type(tDM), target                    :: dm
   PetscBool                           :: flg
   PetscBool                           :: interpolate = PETSC_FALSE
   character(len=MEF90MXSTRLEN)        :: IOBuffer

   MEF90GlobalOptions_default % verbose = 1
   MEF90GlobalOptions_default % dryrun = PETSC_FALSE
   MEF90GlobalOptions_default % timeInterpolation = MEF90TimeInterpolation_linear
   MEF90GlobalOptions_default % timeMin = 0.0_kr
   MEF90GlobalOptions_default % timeMax = 1.0_kr
   MEF90GlobalOptions_default % timeNumStep = 11
   MEF90GlobalOptions_default % timeSkip = 0
   MEF90GlobalOptions_default % timeNumCycle = 1
   MEF90GlobalOptions_default % elementFamily = MEF90ElementFamilyLagrange
   MEF90GlobalOptions_default % elementOrder = 1

   PetscCallA(PetscInitialize(ierr))
   PetscCallA(MEF90Initialize(PETSC_COMM_WORLD, ierr))
   PetscCallA(MEF90CtxCreate(PETSC_COMM_WORLD, MEF90Ctx, MEF90GlobalOptions_default, ierr))

   PetscCallA(DMPlexCreateFromFile(MEF90Ctx % Comm, MEF90Ctx % geometryfile, PETSC_NULL_CHARACTER, interpolate, dm, ierr))
   PetscCallA(DMPlexDistributeSetDefault(dm, PETSC_FALSE, ierr))
   PetscCallA(DMSetFromOptions(dm, ierr))
   PetscCallA(DMViewFromOptions(dm, PETSC_NULL_OBJECT, "-dm_view", ierr))

   distribute: block
      type(tDM), target                    :: dmDist
      PetscInt                            :: ovlp = 0_ki
      if (MEF90Ctx % NumProcs > 1) then
         PetscCallA(DMPlexDistribute(dm, ovlp, PETSC_NULL_SF, dmDist, ierr))
         PetscCallA(DMDestroy(dm, ierr))
         dm = dmDist
      end if
   end block distribute
   PetscCallA(DMViewFromOptions(dm, PETSC_NULL_OBJECT, "-mef90dm_view", ierr))

   ! PetscCallA(MEF90CtxOpenEXO(MEF90Ctx,dm,ierr))
   ! PetscCallA(MEF90CtxGetTime(MEF90Ctx,time,ierr))

   PetscCallA(MEF90CtxDestroy(MEF90Ctx, ierr))

   PetscCallA(PetscOptionsHasName(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-verbose", flg, ierr))
   if (flg) then
      write (IOBuffer, *) "verbose is set\n"
   else
      write (IOBuffer, *) "verbose is NOT set\n"
   end if
   PetscCallA(PetscPrintf(PETSC_COMM_WORLD, IOBuffer, ierr))
   PetscCallA(PetscOptionsHasName(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-time_min", flg, ierr))
   if (flg) then
      write (IOBuffer, *) "time_min is set\n"
   else
      write (IOBuffer, *) "time_min is NOT set\n"
   end if
   PetscCallA(PetscPrintf(PETSC_COMM_WORLD, IOBuffer, ierr))

   PetscCallA(MEF90Finalize(ierr))
   PetscCallA(PetscFinalize(ierr))
end program TestMEF90Ctx
