program TestSection
#include <petsc/finclude/petsc.h>
   use m_MEF90
   use petsc
   implicit none(type, external)

   PetscErrorCode                      :: ierr
   type(MEF90Ctx_Type), target          :: MEF90Ctx
   type(MEF90CtxGlobalOptions_Type)    :: MEF90GlobalOptions_default
   type(tDM)                           :: dm, dmU
   PetscBool                           :: interpolate = PETSC_TRUE
   character(len=MEF90MXSTRLEN)        :: name

   PetscInt                            :: dim
   type(tPetscSection)                 :: section
   type(tVec)                          :: U, U0, V, V0

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

   call MEF90Initialize(PETSC_COMM_WORLD, ierr)
   call MEF90CtxCreate(PETSC_COMM_WORLD, MEF90Ctx, "", ierr)
   MEF90Ctx%globalOptions = MEF90GlobalOptions_default
   call MEF90Ctx%setFromOptions(ierr); CHKERRQ(ierr)

   PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm, MEF90Ctx%geometryfile, PETSC_NULL_CHARACTER, interpolate, dm, ierr))
   PetscCallA(DMPlexDistributeSetDefault(dm, PETSC_FALSE, ierr))
   PetscCallA(DMSetFromOptions(dm, ierr))
   PetscCallA(DMViewFromOptions(dm, PETSC_NULL_OBJECT, "-mef90dm_view", ierr))
   distribute: block
      type(tDM), target                    :: dmDist
      PetscInt                            :: ovlp = 2
      if (MEF90Ctx%NumProcs > 1) then
         PetscCallA(DMPlexDistribute(dm, ovlp, PETSC_NULL_SF, dmDist, ierr))
         PetscCallA(DMDestroy(dm, ierr))
         dm = dmDist
      end if
   end block distribute
   PetscCallA(DMViewFromOptions(dm, PETSC_NULL_OBJECT, "-mef90dm_view", ierr))

   PetscCallA(DMGetDimension(dm, dim, ierr))

   name = "U"
   PetscCallA(MEF90CreateLocalVector(dm, MEF90GlobalOptions_default%elementFamily, MEF90GlobalOptions_default%elementOrder, dim, name, U, ierr))
   name = "U0"
   PetscCallA(MEF90CreateBoundaryLocalVector(dm, MEF90GlobalOptions_default%elementFamily, MEF90GlobalOptions_default%elementOrder, dim, name, U0, ierr))

   name = "V"
   PetscCallA(MEF90CreateCellVector(dm, dim, name, V, ierr))
   name = "V0"
   PetscCallA(MEF90CreateBoundaryCellVector(dm, dim, name, V0, ierr))

   PetscCallA(VecGetDM(U, dmU, ierr))
   PetscCallA(DMGetLocalSection(dmU, section, ierr))
   PetscCallA(PetscSectionViewFromOptions(section, PETSC_NULL_OBJECT, "-sectionU_view", ierr))
   PetscCallA(VecViewFromOptions(U, PETSC_NULL_OBJECT, "-U_view", ierr))

   PetscCallA(VecGetDM(U0, dmU, ierr))
   PetscCallA(DMGetLocalSection(dmU, section, ierr))
   PetscCallA(PetscSectionViewFromOptions(section, PETSC_NULL_OBJECT, "-sectionU0_view", ierr))
   PetscCallA(VecViewFromOptions(U0, PETSC_NULL_OBJECT, "-U0_view", ierr))
   ! dmU and section were obtained by a "Get" without a matching "Restore" so they do not need to be destroyed

   PetscCallA(VecGetDM(V, dmU, ierr))
   PetscCallA(DMGetLocalSection(dmU, section, ierr))
   PetscCallA(PetscSectionViewFromOptions(section, PETSC_NULL_OBJECT, "-sectionV_view", ierr))
   PetscCallA(VecViewFromOptions(V, PETSC_NULL_OBJECT, "-V_view", ierr))

   PetscCallA(VecGetDM(V0, dmU, ierr))
   PetscCallA(DMGetLocalSection(dmU, section, ierr))
   PetscCallA(PetscSectionViewFromOptions(section, PETSC_NULL_OBJECT, "-sectionV0_view", ierr))
   PetscCallA(VecViewFromOptions(V0, PETSC_NULL_OBJECT, "-V0_view", ierr))

   PetscCallA(VecDestroy(U, ierr))
   PetscCallA(VecDestroy(U0, ierr))
   PetscCallA(VecDestroy(V, ierr))
   PetscCallA(VecDestroy(V0, ierr))
   PetscCallA(DMDestroy(dm, ierr))

   call MEF90CtxDestroy(MEF90Ctx, ierr)
   call MEF90Finalize(ierr)
   call PetscFinalize(ierr)
end program TestSection

! mpirun -np 3 ./TestSection -geometry ../TestMeshes/SquareFaceSetCubit2CS.gen -result test.exo -mef90dm_view -log_view