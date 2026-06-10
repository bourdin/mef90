program TestMEF90Ctx
#include <petsc/finclude/petsc.h>
   use m_MEF90
   use petsc
   implicit none(type, external)

   PetscErrorCode                      :: ierr
   type(MEF90Ctx_Type), target         :: MEF90Ctx
   type(MEF90CtxGlobalOptions_Type)    :: MEF90GlobalOptions_default
   type(tDM), target                   :: dm
   PetscBool                           :: flg
   PetscBool                           :: interpolate = PETSC_TRUE
   character(len=MEF90MXSTRLEN)        :: IOBuffer
   character(len=MEF90MXSTRLEN),dimension(:), pointer :: vnameG, vnameZ, vnameN
   integer                             :: i
   PetscReal, dimension(:), pointer    :: time
   PetscViewer :: viewer


   MEF90GlobalOptions_default%verbose = 1
   MEF90GlobalOptions_default%dryrun = PETSC_FALSE
   MEF90GlobalOptions_default%timeInterpolation = MEF90TimeInterpolation_linear
   MEF90GlobalOptions_default%timeMin = 0.0_kr
   MEF90GlobalOptions_default%timeMax = 1.0_kr
   MEF90GlobalOptions_default%timeNumStep = 11
   MEF90GlobalOptions_default%timeSkip = 0
   MEF90GlobalOptions_default%timeNumCycle = 1
   MEF90GlobalOptions_default%elementFamily = MEF90ElementFamilyLagrange
   MEF90GlobalOptions_default%elementOrder = 1

   PetscCallA(PetscInitialize(ierr))
   PetscCallA(MEF90Initialize(PETSC_COMM_WORLD, ierr))
   PetscCallA(MEF90CtxCreate(PETSC_COMM_WORLD, MEF90Ctx, MEF90GlobalOptions_default, ierr))

   PetscCallA(PetscPrintf(PETSC_COMM_WORLD, "Reading geometry file ...", ierr))
   PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm, MEF90Ctx%geometryfile, PETSC_NULL_CHARACTER, interpolate, dm, ierr))
   PetscCallA(PetscPrintf(PETSC_COMM_WORLD, "Done \n", ierr))
   PetscCallA(DMPlexDistributeSetDefault(dm, PETSC_FALSE, ierr))
   PetscCallA(DMSetUseNatural(dm, PETSC_TRUE, ierr))
   PetscCallA(DMSetFromOptions(dm, ierr))
   PetscCallA(DMViewFromOptions(dm, PETSC_NULL_OBJECT, "-dm_view", ierr))

   PetscCallA(PetscPrintf(PETSC_COMM_WORLD, "Saving result file ...", ierr))
   PetscCallA(MEF90CtxOpenEXO(MEF90Ctx, MEF90Ctx%resultViewer, FILE_MODE_WRITE, ierr))
   PetscCallA(MEF90EXODMView(dm, MEF90Ctx%resultViewer, MEF90GlobalOptions_default%elementOrder, ierr))
   PetscCallA(PetscPrintf(PETSC_COMM_WORLD, "Done \n", ierr))


   distribute: block
      type(tDM), target                    :: dmDist
      PetscInt                             :: ovlp = 0_ki

      PetscCallA(PetscPrintf(PETSC_COMM_WORLD, "Distributing mesh ...", ierr))
      if (MEF90Ctx%NumProcs > 1) then
         PetscCallA(DMPlexDistribute(dm, ovlp, PETSC_NULL_SF, dmDist, ierr))
         PetscCallA(DMDestroy(dm, ierr))
         dm = dmDist
      end if
      PetscCallA(PetscPrintf(PETSC_COMM_WORLD, "Done \n", ierr))
   end block distribute
   PetscCallA(DMViewFromOptions(dm, PETSC_NULL_OBJECT, "-mef90dm_view", ierr))

   allocate(vnameG(0))
   allocate(vnameN(3))
   do i = 1, 3
      write (vnameN(i), '(A,I0)') "Nodal_", i
   end do
   allocate(vnameZ(6))
   do i = 1, 6
      write (vnameZ(i), '(A,I0)') "Zonal_", i
   end do

   allocate(time(101))
   time = (/(i, i = 1, 101)/)


   ! ! PetscCallA(PetscViewerExodusIIOpen(PETSC_COMM_WORLD, MEF90Ctx%resultfile, FILE_MODE_WRITE, viewer, ierr))
   PetscCallA(PetscPrintf(PETSC_COMM_WORLD, "Formatting result file ...", ierr))
   Call MEF90EXOFormat(MEF90Ctx%resultViewer, vnameG, vnameZ, vnameN, time, ierr)
   PetscCallA(PetscPrintf(PETSC_COMM_WORLD, "Done \n", ierr))

   PetscCallA(MEF90CtxDestroy(MEF90Ctx, ierr))

   PetscCallA(MEF90Finalize(ierr))
   PetscCallA(PetscFinalize(ierr))
end program TestMEF90Ctx
