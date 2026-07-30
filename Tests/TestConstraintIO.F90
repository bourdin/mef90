module localFunctions
#include <petsc/finclude/petsc.h>
   use petsc
   use m_MEF90
   implicit none(type, external)

contains

#undef __FUNCT__
#define __FUNCT__ "project"

   subroutine project(v, s, ierr)
      type(tVec), intent(IN)              :: v
      type(tPetscSection), intent(IN)     :: s
      PetscErrorCode, intent(INOUT)       :: ierr

      PetscInt                           :: pStart, pEnd, p, numDof, i
      type(tDM)                          :: dm
      type(tPetscSection)                :: coordSection
      type(tVec)                         :: coordVec
      PetscScalar, dimension(:), pointer   :: coordArray, vArray
      PetscScalar, dimension(3)           :: xyz
      PetscInt                           :: dim, pOffset

      PetscCallA(PetscSectionGetChart(s, pStart, pEnd, ierr))
      PetscCallA(VecGetDM(v, dm, ierr))
      PetscCallA(DMGetCoordinateSection(dm, coordSection, ierr))
      PetscCallA(DMGetCoordinatesLocal(dm, coordVec, ierr))
      PetscCallA(DMGetDimension(dm, dim, ierr))
      PetscCallA(VecGetArray(v, vArray, ierr))

      do p = pStart, pEnd - 1
         PetscCallA(PetscSectionGetDof(s, p, numDof, ierr))
         if (numDof > 0) then
                !!! trick: the coordinate of a point is the average of the coordinates of the points in its closure
            PetscCallA(DMPlexVecGetClosure(dm, coordSection, coordVec, p, coordArray, ierr))
            do i = 1, dim
               xyz(i) = sum(coordArray(i:size(coordArray):dim)) * dim / size(coordArray)
            end do
            PetscCallA(DMPlexVecRestoreClosure(dm, coordSection, coordVec, p, coordArray, ierr))

            PetscCallA(PetscSectionGetOffset(s, p, pOffset, ierr))
            do i = 1, numDof
               vArray(pOffset + i) = xyz(i)
            end do
         end if
      end do
      PetscCallA(VecRestoreArray(v, vArray, ierr))
        !!! Of course, this does not use informations from the section, so it does over-write constrained values
   end subroutine project

end module localFunctions

program TestConstraintIO
#include <petsc/finclude/petsc.h>
   use m_MEF90
   use petsc
   use localFunctions
   implicit none(type, external)

   PetscErrorCode                                      :: ierr
   type(MEF90Ctx_Type), target                          :: MEF90Ctx
   type(MEF90CtxGlobalOptions_Type)                    :: MEF90GlobalOptions_default
   type(MEF90CtxGlobalOptions_Type)                     :: MEF90GlobalOptions
   type(tDM), target                                    :: dm, dmU, dmU0, dmSigma, dmSigma0
   PetscBool                                           :: interpolate = PETSC_TRUE

   PetscInt                                            :: numNodalVar = 3, numCellVar = 3, numGVar = 0
   PetscInt                                            :: numFS
   character(len=MEF90MXSTRLEN), dimension(:), pointer   :: nodalVarName, cellVarName, gVarName
   character(len=MEF90MXSTRLEN)                        :: name
   type(tIS)                                           :: cellIS, csIS, faceIS, ssIS
   PetscInt, dimension(:), pointer                       :: cellID, csID, faceID, ssID
   PetscInt                                            :: dim, pStart, pEnd, size, cell, c, set, face
   type(tPetscSection)                                 :: sectionU, sectionU0, sectionSigma, sectionSigma0, coordSection, coordSection0
   type(tPetscSF)                                      :: naturalPointSF, lcSF, clSF, lcSSF, clSSF, lioSF, iolSF, lioBSF, iolBSF, lioSSF, iolSSF, lioBSSF, iolBSSF
   type(tVec)                                          :: locCoord, locVecU, locVecU0, locVecSigma, locVecSigma0
   PetscReal, dimension(:), pointer                      :: cval, xyz
   PetscReal, dimension(:), pointer                      :: time
   PetscInt                                            :: step = 1_ki

   PetscCallA(PetscInitialize(ierr))
   PetscCallA(MEF90Initialize(PETSC_COMM_WORLD, ierr))

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

   PetscCallA(MEF90CtxCreate(PETSC_COMM_WORLD, MEF90Ctx, "", ierr))
   MEF90GlobalOptions = MEF90GlobalOptions_default
   PetscCallA(MEF90Ctx%setFromOptions(ierr))
   PetscCallA(MEF90CtxGlobalOptionsSetFromOptions(MEF90Ctx%comm, trim(MEF90Ctx%prefix), MEF90GlobalOptions, ierr))
   PetscCallA(MEF90CtxGetTime(MEF90Ctx, time, ierr))

   allocate (nodalVarName(numNodalVar))
   allocate (cellVarName(numCellVar))
   allocate (gVarName(numGVar))
   nodalVarName = ["U_X", "U_Y", "U_Z"]
   cellVarName = ["Sigma_11", "Sigma_22", "Sigma_33"]

   ! Create DM from file
   PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm, MEF90Ctx%geometryfile, PETSC_NULL_CHARACTER, interpolate, dm, ierr))
   PetscCallA(DMPlexDistributeSetDefault(dm, PETSC_FALSE, ierr))
   PetscCallA(DMSetUseNatural(dm, PETSC_TRUE, ierr))
   PetscCallA(DMSetFromOptions(dm, ierr))
   PetscCallA(DMViewFromOptions(dm, PETSC_NULL_OBJECT, "-mef90dm_view", ierr))

   ! Open exodus file + write geometry + format the file
   PetscCallA(MEF90CtxOpenEXO(MEF90Ctx, MEF90Ctx%resultViewer, FILE_MODE_WRITE, ierr))
   PetscCallA(MEF90EXODMView(dm, MEF90Ctx%resultViewer, MEF90GlobalOptions%elementOrder, ierr))
   PetscCallA(MEF90EXOFormat(MEF90Ctx%resultViewer, gVarName, cellVarName, nodalVarName, time, ierr))

   deallocate (nodalVarName)
   deallocate (cellVarName)
   deallocate (gVarName)

   ! Distribute DM
   distribute: block
      type(tDM), target                    :: dmDist
      PetscInt                            :: ovlp = 0
      if (MEF90Ctx%NumProcs > 1) then
         PetscCallA(DMPlexDistribute(dm, ovlp, naturalPointSF, dmDist, ierr))
         PetscCallA(DMPlexSetMigrationSF(dmDist, naturalPointSF, ierr))
         PetscCallA(PetscSFDestroy(naturalPointSF, ierr))
         PetscCallA(DMDestroy(dm, ierr))
         dm = dmDist
         PetscCallA(DMSetUseNatural(dm, PETSC_TRUE, ierr))
      end if
   end block distribute
   PetscCallA(DMViewFromOptions(dm, PETSC_NULL_OBJECT, "-mef90dm_view", ierr))

   PetscCallA(DMGetDimension(dm, dim, ierr))
   PetscCallA(DMPlexGetChart(dm, pStart, pEnd, ierr))

   ! Create nodal local Vec holding coordinates
   name = "U"
   PetscCallA(MEF90CreateLocalVector(dm, MEF90GlobalOptions%elementFamily, MEF90GlobalOptions%elementOrder, dim, name, locVecU, ierr))
   PetscCallA(VecGetDM(locVecU, dmU, ierr))
   PetscCallA(DMGetLocalSection(dmU, sectionU, ierr))

   ! Create nodal local Vec holding constraints
   name = "U0"
   PetscCallA(MEF90CreateBoundaryLocalVector(dm, MEF90GlobalOptions%elementFamily, MEF90GlobalOptions%elementOrder, dim, name, locVecU0, ierr))
   PetscCallA(VecGetDM(locVecU0, dmU0, ierr))
   PetscCallA(DMGetLocalSection(dmU0, sectionU0, ierr))

   ! create cell Vec holding sigma
   name = "Sigma"
   PetscCallA(MEF90CreateCellVector(dm, 3_ki, name, locVecSigma, ierr))
   PetscCallA(VecGetDM(locVecSigma, dmSigma, ierr))
   PetscCallA(DMGetLocalSection(dmSigma, sectionSigma, ierr))

   ! create cell Vec holding sigma
   name = "Sigma0"
   PetscCallA(MEF90CreateBoundaryCellVector(dm, 3_ki, name, locVecSigma0, ierr))
   PetscCallA(VecGetDM(locVecSigma0, dmSigma0, ierr))
   PetscCallA(DMGetLocalSection(dmSigma0, sectionSigma0, ierr))

   ! View all sections
   PetscCallA(PetscSectionViewFromOptions(SectionU, PETSC_NULL_OBJECT, "-sectionU_view", ierr))
   PetscCallA(PetscSectionViewFromOptions(SectionU0, PETSC_NULL_OBJECT, "-sectionU0_view", ierr))
   PetscCallA(PetscSectionViewFromOptions(SectionSigma, PETSC_NULL_OBJECT, "-sectionSigma_view", ierr))
   PetscCallA(PetscSectionViewFromOptions(SectionSigma0, PETSC_NULL_OBJECT, "-sectionSigma0_view", ierr))

   ! Create SFs for copying from/into IO coordinates Vec
   PetscCallA(MEF90IOSFCreate(MEF90Ctx, locVecU, lioSF, iolSF, ierr))
   ! Create SFs for copying from/into IO Constraint Vec
   PetscCallA(MEF90IOSFCreate(MEF90Ctx, locVecU0, lioBSF, iolBSF, ierr))
   ! Create SFs for copying constrained dofs from/into Constraint Vec
   PetscCallA(MEF90ConstraintSFCreate(MEF90Ctx, locVecU, locVecU0, lcSF, clSF, ierr))
   ! Create SFs for copying from/into IO cell Vec
   PetscCallA(MEF90IOSFCreate(MEF90Ctx, locVecSigma, lioSSF, iolSSF, ierr))
   ! Create SFs for copying from/into IO cell Vec
   PetscCallA(MEF90FaceSetIOSFCreate(MEF90Ctx, locVecSigma0, lioBSSF, iolBSSF, ierr))
   ! Create SFs for copying constrained dofs from/into Constraint Vec
   PetscCallA(MEF90ConstraintSFCreate(MEF90Ctx, locVecSigma, locVecSigma0, lcSSF, clSSF, ierr))

   ! Fill locVecU holding coordinates and copy constraint values = 10 from locVecU0
   PetscCallA(VecSet(locVecU0, 10.0_kr, ierr))
   ! locCoord is obtained from dmSigma but all DMs have the same coordinates
   PetscCallA(DMGetCoordinatesLocal(dmSigma, locCoord, ierr))
   ! Fill locVecU with coordinates
   PetscCallA(project(locVecU, sectionU, ierr))
   PetscCallA(MEF90VecCopySF(locVecU0, locVecU, clSF, ierr))

   ! Reorder locVecU into ioVec and write ioVec
   PetscCallA(VecViewFromOptions(locVecU, PETSC_NULL_OBJECT, "-iovec_view", ierr))
   PetscCallA(MEF90EXOVecView(locVecU, lioSF, iolSF, MEF90Ctx%resultViewer, step, dim, ierr))

   ! Create Vec with 3 components on each cell: (i) rank, (ii) x geometric center,
   ! (iii) y geometric center
   PetscCallA(DMGetCoordinateSection(dmSigma0, coordSection0, ierr))
   PetscCallA(DMGetLabelIdIS(dmSigma0, "Face Sets", ssIS, ierr))
   PetscCallA(ISGetIndices(ssIS, ssID, ierr))
   numFS = size(ssID)
   do set = 1, numFS
      PetscCallA(DMGetStratumIS(dmSigma0, "Face Sets", ssID(set), faceIS, ierr))
      if (faceIS /= PETSC_NULL_IS) then
         PetscCallA(ISGetIndices(faceIS, faceID, ierr))
         do face = 1, size(faceID)
            PetscCallA(DMPlexVecGetClosure(dmSigma0, PETSC_NULL_SECTION, locVecSigma0, faceID(face), cval, ierr))
            PetscCallA(DMPlexVecGetClosure(dmSigma0, coordSection0, locCoord, faceID(face), xyz, ierr))
            cval(1) = MEF90Ctx%rank
            cval(2) = 0.0_kr
            do c = 1, size(xyz), dim
               cval(2) = cval(2) + xyz(c)
            end do
            cval(2) = cval(2) * dim / size(xyz)
            cval(3) = 0.0_kr
            do c = 0, size(xyz), dim
               cval(3) = cval(3) + xyz(c)
            end do
            cval(3) = cval(3) * dim / size(xyz)
            !write(*,*) 'RANK = ',MEF90Ctx%rank,' set = ',set,' set ID = ',ssID(set),' face = ',face,' faceID = ',faceID(face),' x = ',cval(2),' y = ',cval(3)
            PetscCallA(DMPlexVecRestoreClosure(dmSigma0, coordSection0, locCoord, faceID(face), xyz, ierr))
            PetscCallA(DMPlexVecSetClosure(dmSigma0, PETSC_NULL_SECTION, locVecSigma0, faceID(face), cval, INSERT_VALUES, ierr))
            PetscCallA(DMPlexVecRestoreClosure(dmSigma0, PETSC_NULL_SECTION, locVecSigma0, faceID(face), cval, ierr))
         end do
         PetscCallA(ISRestoreIndices(faceIS, faceID, ierr))
      end if ! cellIS
      PetscCallA(ISDestroy(faceIS, ierr))
   end do
   PetscCallA(ISRestoreIndices(ssIS, ssID, ierr))
   PetscCallA(ISDestroy(ssIS, ierr))
   PetscCallA(DMGetCoordinateSection(dmSigma, coordSection, ierr))
   PetscCallA(DMGetLabelIdIS(dmSigma, "Cell Sets", csIS, ierr))
   PetscCallA(ISGetIndices(csIS, csID, ierr))
   do set = 1, size(csID)
      PetscCallA(DMGetStratumIS(dmSigma, "Cell Sets", csID(set), cellIS, ierr))
      if (cellIS /= PETSC_NULL_IS) then
         PetscCallA(ISGetIndices(cellIS, cellID, ierr))
         do cell = 1, size(cellID)
            PetscCallA(DMPlexVecGetClosure(dmSigma, PETSC_NULL_SECTION, locVecSigma, cellID(cell), cval, ierr))
            PetscCallA(DMPlexVecGetClosure(dmSigma, coordSection, locCoord, cellID(cell), xyz, ierr))
            cval(1) = MEF90Ctx%rank
            cval(2) = 0.0_kr
            do c = 1, size(xyz), dim - 1
               cval(2) = cval(2) + xyz(c)
            end do
            cval(2) = cval(2) * dim / size(xyz)
            cval(3) = 0.0_kr
            do c = 1, size(xyz), dim
               cval(3) = cval(3) + xyz(c)
            end do
            cval(3) = cval(3) * dim / size(xyz)
            PetscCallA(DMPlexVecRestoreClosure(dmSigma, coordSection, locCoord, cellID(cell), xyz, ierr))
            PetscCallA(DMPlexVecSetClosure(dmSigma, PETSC_NULL_SECTION, locVecSigma, cellID(cell), cval, INSERT_ALL_VALUES, ierr))
            PetscCallA(DMPlexVecRestoreClosure(dmSigma, PETSC_NULL_SECTION, locVecSigma, cellID(cell), cval, ierr))
         end do
         PetscCallA(ISRestoreIndices(cellIS, cellID, ierr))
      end if ! cellIS
      PetscCallA(ISDestroy(cellIS, ierr))
   end do
   PetscCallA(ISRestoreIndices(csIS, csID, ierr))
   PetscCallA(ISDestroy(csIS, ierr))

   ! Reorder and write ioS
   PetscCallA(VecViewFromOptions(locVecSigma, PETSC_NULL_OBJECT, "-ios_view", ierr))
   PetscCallA(MEF90EXOVecView(locVecSigma, lioSSF, iolSSF, MEF90Ctx%resultViewer, step, dim * (dim + 1) / 2, ierr))

   ! Reorder and write ioS0
   if (numFS > 0) then
      PetscCallA(VecViewFromOptions(locVecSigma0, PETSC_NULL_OBJECT, "-ios0_view", ierr))
      PetscCallA(MEF90EXOVecView(locVecSigma0, lioBSSF, iolBSSF, MEF90Ctx%resultViewer, step, dim * (dim + 1) / 2, ierr))
   end if

   ! Test read ioVecRead and ioSRead
   PetscCallA(VecSet(locVecU, 1000.0_kr, ierr))
   PetscCallA(MEF90EXOVecLoad(locVecU, lioSF, iolSF, MEF90Ctx%resultViewer, step, dim, ierr))
   PetscCallA(VecViewFromOptions(locVecU, PETSC_NULL_OBJECT, "-iovec_view", ierr))
   PetscCallA(VecSet(locVecSigma, 1000.0_kr, ierr))
   PetscCallA(MEF90EXOVecLoad(locVecSigma, lioSSF, iolSSF, MEF90Ctx%resultViewer, step, dim * (dim + 1) / 2_ki, ierr))
   PetscCallA(VecViewFromOptions(locVecSigma, PETSC_NULL_OBJECT, "-ios_view", ierr))
   PetscCallA(VecSet(locVecSigma0, 1000.0_kr, ierr))
   if (numFS > 0) then
      PetscCallA(MEF90EXOVecLoad(locVecSigma0, lioBSSF, iolBSSF, MEF90Ctx%resultViewer, step, dim * (dim + 1) / 2_ki, ierr))
      PetscCallA(VecViewFromOptions(locVecSigma0, PETSC_NULL_OBJECT, "-ios0_view", ierr))
   end if

   ! Cleanup Vec
   PetscCallA(VecDestroy(locVecU, ierr))
   PetscCallA(VecDestroy(locVecU0, ierr))
   PetscCallA(VecDestroy(locVecSigma, ierr))
   PetscCallA(VecDestroy(locVecSigma0, ierr))

   ! Cleanup SF
   PetscCallA(PetscSFDestroy(lcSF, ierr))
   PetscCallA(PetscSFDestroy(clSF, ierr))
   PetscCallA(PetscSFDestroy(lcSSF, ierr))
   PetscCallA(PetscSFDestroy(clSSF, ierr))
   PetscCallA(PetscSFDestroy(lioSF, ierr))
   PetscCallA(PetscSFDestroy(iolSF, ierr))
   PetscCallA(PetscSFDestroy(lioBSF, ierr))
   PetscCallA(PetscSFDestroy(iolBSF, ierr))
   PetscCallA(PetscSFDestroy(lioSSF, ierr))
   PetscCallA(PetscSFDestroy(iolSSF, ierr))
   PetscCallA(PetscSFDestroy(lioBSSF, ierr))
   PetscCallA(PetscSFDestroy(iolBSSF, ierr))

   ! Cleanup DMs
   deallocate (time)
   ! Note that I would need to manually destroy these DM no matter what
   PetscCallA(DMDestroy(dm, ierr))

   ! Exit nicely
   PetscCallA(MEF90CtxDestroy(MEF90Ctx, ierr))
   PetscCallA(MEF90Finalize(ierr))
   PetscCallA(PetscFinalize(ierr))
end program TestConstraintIO
