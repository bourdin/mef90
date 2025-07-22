program TestOffsets
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
   type(tVec)                          :: U

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
   call MEF90CtxCreate(PETSC_COMM_WORLD, MEF90Ctx, MEF90GlobalOptions_default, ierr)

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

   PetscCallA(VecGetDM(U, dmU, ierr))
   PetscCallA(DMGetLocalSection(dmU, section, ierr))
   PetscCallA(PetscSectionViewFromOptions(section, PETSC_NULL_OBJECT, "-sectionU_view", ierr))
   PetscCallA(VecViewFromOptions(U, PETSC_NULL_OBJECT, "-U_view", ierr))

   PetscCallA(VecDestroy(U, ierr))

   testRemote: block
      type(tPetscSection)                     :: locsection, gsection
      type(tPetscSF)                          :: overlapSF, idSF, sf
      PetscInt                                :: pStart, pEnd, n, p
      type(sPetscSFNode), dimension(:), pointer :: remote
      PetscInt, dimension(:), pointer           :: remoteOffsets

      PetscCallA(DMGetLocalSection(dm, locSection, ierr))
      PetscCall(DMGetPointSF(dm, overlapSF, ierr))
      PetscCallA(PetscSectionCreateGlobalSection(locSection, overlapSF, PETSC_TRUE, PETSC_TRUE, PETSC_TRUE, gSection, ierr))

      PetscCall(PetscSectionGetChart(locSection, pStart, pEnd, ierr))
      n = pEnd - pStart
      allocate (remote(n))
      do p = 1, n
         remote(p)%rank = MEF90Ctx%rank
         remote(p)%index = p - 1
      end do

      PetscCall(PetscSFCreate(MEF90Ctx%Comm, idSF, ierr))
      PetscCall(PetscSFSetFromOptions(idSF, ierr))
      PetscCall(PetscSFSetGraph(idSF, n, n, PETSC_NULL_INTEGER_ARRAY, PETSC_COPY_VALUES, remote, PETSC_COPY_VALUES, ierr))
      PetscCall(PetscSFSetUp(idSF, ierr))
      PetscCall(PetscSFCreateRemoteOffsets(idSF, locSection, gSection, remoteOffsets, ierr))
      PetscCallA(PetscSectionViewFromOptions(locsection, PETSC_NULL_OBJECT, "-locsection_view", ierr))
      PetscCallA(PetscSectionViewFromOptions(gsection, PETSC_NULL_OBJECT, "-gsection_view", ierr))
      if (associated(remoteOffsets)) then
         PetscCall(PetscSFCreateSectionSF(idSF, locSection, remoteOffsets, gSection, sf, ierr))
         PetscCall(PetscSFDestroyRemoteOffsets(remoteOffsets, ierr))
      end if
   end block testRemote

   PetscCallA(DMDestroy(dm, ierr))

   call MEF90CtxDestroy(MEF90Ctx, ierr)
   call MEF90Finalize(ierr)
   call PetscFinalize(ierr)
end program TestOffsets

! mpirun -np 3 ./TestOffsets -geometry ../TestMeshes/SquareFaceSetCubit2CS.gen -result test.exo -mef90dm_view -log_view