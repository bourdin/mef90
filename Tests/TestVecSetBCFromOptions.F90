program TestVecSetBCFromOptions
#include <petsc/finclude/petsc.h>
   use petsc
   use m_MEF90
   implicit none(type, external)

   PetscErrorCode                                     :: ierr
   type(MEF90Ctx_Type), target                         :: MEF90Ctx
   type(MEF90CtxGlobalOptions_Type)                   :: MEF90GlobalOptions_default
   type(MEF90CtxGlobalOptions_Type)                    :: MEF90GlobalOptions
   type(tDM)                                          :: dm
   PetscBool                                          :: flg, interpolate = PETSC_TRUE
   character(len=MEF90MXSTRLEN)                       :: name
   PetscInt                                           :: sDim = 1

   type(tVec)                                         :: v
   PetscReal                                          :: scalingFactor

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
   MEF90GlobalOptions = MEF90GlobalOptions_default
   call MEF90Ctx%setFromOptions(ierr); CHKERRQ(ierr)
   PetscCallA(MEF90CtxGlobalOptionsSetFromOptions(MEF90Ctx%comm, trim(MEF90Ctx%prefix), MEF90GlobalOptions, ierr))
   PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm, MEF90Ctx%geometryfile, PETSC_NULL_CHARACTER, interpolate, dm, ierr))
   PetscCallA(DMPlexDistributeSetDefault(dm, PETSC_FALSE, ierr))
   PetscCallA(DMSetFromOptions(dm, ierr))
   PetscCallA(DMViewFromOptions(dm, PETSC_NULL_OBJECT, "-mef90dm_view", ierr))

   distribute: block
      type(tDM), target                    :: dmDist
      PetscInt                            :: ovlp = 0
      if (MEF90Ctx%NumProcs > 1) then
         PetscCallA(DMPlexDistribute(dm, ovlp, PETSC_NULL_SF, dmDist, ierr))
         PetscCallA(DMDestroy(dm, ierr))
         dm = dmDist
      end if
   end block distribute

   name = "Temperature"
   PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-sdim', sdim, flg, ierr))
   PetscCallA(MEF90CreateLocalVector(dm, MEF90GlobalOptions%elementFamily, MEF90GlobalOptions%elementOrder, sdim, name, V, ierr))
   PetscCallA(VecSet(v, -1.234_kr, ierr))
   scalingFactor = 1.0_kr
   PetscCallA(MEF90VecSetBCValuesFromOptions(V, scalingFactor, ierr))
   PetscCallA(VecViewFromOptions(V, PETSC_NULL_OBJECT, "-temperature_vec_view", ierr))

   PetscCallA(VecSet(v, 5.678_kr, ierr))
   scalingFactor = 2.0_kr
   PetscCallA(MEF90VecSetValuesFromOptions(V, scalingFactor, ierr))
   PetscCallA(VecViewFromOptions(V, PETSC_NULL_OBJECT, "-temperature_vec_view", ierr))

   ViewSec: block
      type(tPetscSection)     :: sectionV
      type(tDM)               :: dmV

      PetscCallA(VecGetDM(V, dmV, ierr))
      PetscCallA(DMGetLocalSection(dmV, sectionV, ierr))
      PetscCallA(PetscSectionViewFromOptions(sectionV, PETSC_NULL_OBJECT, "-temperature_section_view", ierr))
   end block ViewSec

   PetscCallA(VecDestroy(v, ierr))
   PetscCallA(DMDestroy(dm, ierr))
   call MEF90CtxDestroy(MEF90Ctx, ierr)
   call MEF90Finalize(ierr)
   call PetscFinalize(ierr)
end program TestVecSetBCFromOptions

!!! TEST
!!! ./TestVecSetBCFromOptions -geometry ../TestMeshes/SquareFaceSet.gen -temperature_section_view    \
!!!     -temperature_vec_view -sdim 2 -fs0021_temperatureBC yes,no -fs0021_boundaryTemperature 10,20 \
!!!     -vs0010_temperatureBC yes,yes -vs0010_boundaryTemperature 100,101 -cs0001_temperature 1000   \
!!!     -vs0010_temperature -2 -fs0020_temperature 20 -fs0021_temperature 21 -vs0010_temperature 100
