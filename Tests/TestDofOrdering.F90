program TestDofOrdering
#include <petsc/finclude/petsc.h>
   use m_MEF90
   use petsc
   implicit none(type, external)

   PetscErrorCode                                 :: ierr
   type(MEF90Ctx_Type), target                     :: MEF90Ctx
   type(MEF90CtxGlobalOptions_Type)               :: MEF90GlobalOptions_default
   type(tDM)                                      :: dm, dmU
   type(tPetscSection)                            :: sectionU
   PetscBool                                      :: interpolate = PETSC_TRUE
   character(len=MEF90MXSTRLEN)                   :: name
   type(MEF90CtxGlobalOptions_Type)                :: MEF90GlobalOptions

   PetscInt                                       :: nVal, i, dim
   type(tVec)                                     :: U
   PetscReal, dimension(:), pointer                 :: UArray

   PetscCallA(PetscInitialize(ierr))
   PetscCallA(MEF90Initialize(PETSC_COMM_WORLD, ierr))

   MEF90GlobalOptions_default%verbose = 1
   MEF90GlobalOptions_default%dryrun = PETSC_FALSE
   MEF90GlobalOptions_default%timeMin = 0.0_kr
   MEF90GlobalOptions_default%timeMax = 1.0_kr
   MEF90GlobalOptions_default%timeNumStep = 11
   MEF90GlobalOptions_default%timeInterpolation = MEF90TimeInterpolation_linear
   MEF90GlobalOptions_default%timeSkip = 0
   MEF90GlobalOptions_default%timeNumCycle = 1
   MEF90GlobalOptions_default%elementFamily = MEF90ElementFamilyLagrange
   MEF90GlobalOptions_default%elementOrder = 1

   PetscCallA(MEF90CtxCreate(PETSC_COMM_WORLD, MEF90Ctx, "", ierr))
   MEF90GlobalOptions = MEF90GlobalOptions_default
   PetscCallA(MEF90Ctx%setFromOptions(ierr))
   PetscCallA(MEF90CtxGlobalOptionsSetFromOptions(MEF90Ctx%comm, trim(MEF90Ctx%prefix), MEF90GlobalOptions, ierr))
   PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm, MEF90Ctx%geometryfile, PETSC_NULL_CHARACTER, interpolate, dm, ierr))
   PetscCallA(DMPlexDistributeSetDefault(dm, PETSC_FALSE, ierr))
   PetscCallA(DMSetFromOptions(dm, ierr))
   PetscCallA(DMViewFromOptions(dm, PETSC_NULL_OBJECT, "-dm_view", ierr))
   PetscCallA(DMGetDimension(dm, dim, ierr))

   distribute: block
      type(tDM), target                    :: dmDist
      PetscInt                            :: ovlp = 0
      if (MEF90Ctx%NumProcs > 1) then
         PetscCallA(DMPlexDistribute(dm, ovlp, PETSC_NULL_SF, dmDist, ierr))
         PetscCallA(DMDestroy(dm, ierr))
         dm = dmDist
         PetscCallA(DMSetUseNatural(dm, PETSC_TRUE, ierr))
      end if
   end block distribute
   PetscCallA(DMViewFromOptions(dm, PETSC_NULL_OBJECT, "-mef90dm_view", ierr))

   ! Create nodal local Vec holding coordinates
   name = "U"
   PetscCallA(MEF90CreateLocalVector(dm, MEF90GlobalOptions%elementFamily, MEF90GlobalOptions%elementOrder, 1_ki, name, U, ierr))
   PetscCallA(VecGetDM(U, dmU, ierr))
   PetscCallA(DMGetLocalSection(dmU, sectionU, ierr))
   PetscCallA(PetscSectionViewFromOptions(sectionU, PETSC_NULL_OBJECT, "-mef90sectionU_view", ierr))
   PetscCallA(VecSet(U, 0.0_kr, ierr))

   PetscCallA(VecGetLocalSize(U, nVal, ierr))
   do i = 1, nVal
      PetscCallA(VecGetArray(U, UArray, ierr))
      UArray = 0.0_kr
      UArray(i) = 1.0_kr
      PetscCallA(VecRestoreArray(U, UArray, ierr))
      PetscCallA(DMPlexVecGetClosure(dmU, sectionU, U, 0_ki, UArray, ierr))
      write (*, *) i, UArray
      PetscCallA(DMPlexVecRestoreClosure(dmU, sectionU, U, 0_ki, UArray, ierr))
   end do

   PetscCallA(VecDestroy(U, ierr))
   PetscCallA(DMDestroy(dm, ierr))

   call MEF90CtxDestroy(MEF90Ctx, ierr)
   call MEF90Finalize(ierr)
   call PetscFinalize(ierr)
end program TestDofOrdering

