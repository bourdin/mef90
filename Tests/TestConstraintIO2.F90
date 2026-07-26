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

program TestConstraintIO2
#include <petsc/finclude/petsc.h>
   use m_MEF90
   use petsc
   use localFunctions
   implicit none(type, external)

   PetscErrorCode                                      :: ierr
   type(MEF90Ctx_Type), target                          :: MEF90Ctx
   type(MEF90CtxGlobalOptions_Type)                    :: MEF90GlobalOptions_default
   type(MEF90CtxGlobalOptions_Type), pointer            :: MEF90GlobalOptions
   type(tDM), target                                    :: dm, dmU
   PetscBool                                           :: interpolate = PETSC_TRUE

   PetscInt                                            :: numNodalVar = 3, numCellVar = 0, numGVar = 0
   character(len=MEF90MXSTRLEN), dimension(:), pointer   :: nodalVarName, cellVarName, gVarName
   character(len=MEF90MXSTRLEN)                        :: name, IOBuffer
   PetscInt                                            :: dim, pStart, pEnd
   type(tPetscSection)                                 :: sectionU
   type(tPetscSF)                                      :: naturalPointSF, lcSF, clSF, lcSSF, clSSF, lioSF, iolSF, lioBSF, iolBSF, lioSSF, iolSSF, lioBSSF, iolBSSF
   type(tVec)                                          :: locVecU, U, U0, locVecV, V
   PetscReal, dimension(:), pointer                      :: time
   PetscInt                                            :: step = 1_ki
   PetscReal                                           :: myerr, err

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
   MEF90Ctx%globalOptions = MEF90GlobalOptions_default
   PetscCallA(MEF90Ctx%setFromOptions(ierr))
   MEF90GlobalOptions => MEF90Ctx%globalOptions
   PetscCallA(MEF90CtxGetTime(MEF90Ctx, time, ierr))

   allocate (nodalVarName(numNodalVar))
   allocate (cellVarName(numCellVar))
   allocate (gVarName(numGVar))
   nodalVarName = ["U_X", "U_Y", "U_Z"]

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
   PetscCallA(DMCreateGlobalVector(dmU, U, ierr))

   ! Fill locVecU with coordinates
   PetscCallA(project(locVecU, sectionU, ierr))
   PetscCallA(DMLocalToGlobal(dmU, locVecU, INSERT_VALUES, U, ierr))
   PetscCallA(VecViewFromOptions(locVecU, PETSC_NULL_OBJECT, "-Uloc_view", ierr))
   PetscCallA(VecViewFromOptions(U, PETSC_NULL_OBJECT, "-U_view", ierr))

   ! Save locVecU in exo file
   PetscCallA(MEF90EXOVecView(locVecU, lioSF, iolSF, MEF90Ctx%resultViewer, step, dim, ierr))

   ! Create new vectors and read them from the file
   PetscCallA(VecDuplicate(locVecU, locVecV, ierr))
   PetscCallA(VecDuplicate(U, V, ierr))
   PetscCallA(VecSet(locVecV, -1000.0_kr, ierr))
   PetscCallA(PetscObjectSetName(locVecV, "U", ierr))
   PetscCallA(VecSet(locVecV, 0.0_kr, ierr))
   PetscCallA(MEF90EXOVecLoad(locVecV, lioSF, iolSF, MEF90Ctx%resultViewer, step, dim, ierr))

   PetscCallA(VecViewFromOptions(locVecV, PETSC_NULL_OBJECT, "-Vloc_view", ierr))
   PetscCallA(VecViewFromOptions(V, PETSC_NULL_OBJECT, "-V_view", ierr))

   ! Save it again
   PetscCallA(MEF90EXOVecView(locVecU, lioSF, iolSF, MEF90Ctx%resultViewer, step + 1, dim, ierr))

   ! Compute the difference between the LOCAL vector we wrote and the one we read
   PetscCallA(VecAXPY(locVecV, -1.0_kr, locVecU, ierr))
   PetscCallA(VecNorm(locVecV, NORM_INFINITY, err, ierr))
   write (IOBuffer, '("Local vector L^infty error:  ",ES12.5,"\n")') err
   PetscCallA(PetscPrintf(PETSC_COMM_WORLD, IOBuffer, ierr))
   PetscCallA(VecViewFromOptions(locVecV, PETSC_NULL_OBJECT, "-diffloc_view", ierr))

   ! Compute the difference between the LOCAL vector we wrote and the one we read
   PetscCallA(VecAXPY(V, -1.0_kr, U, ierr))
   PetscCallA(VecNorm(V, NORM_INFINITY, myerr, ierr))
   PetscCallMPIA(MPI_AllReduce(myerr, err, 1, MPIU_SCALAR, MPI_MAX, PETSC_COMM_WORLD, ierr))
   write (IOBuffer, '("Global vector L^infty error: ",ES12.5,"\n")') err
   PetscCallA(PetscPrintf(PETSC_COMM_WORLD, IOBuffer, ierr))
   PetscCallA(VecViewFromOptions(V, PETSC_NULL_OBJECT, "-diff_view", ierr))

   ! Cleanup Vec
   PetscCallA(VecDestroy(locVecU, ierr))
   PetscCallA(VecDestroy(locVecV, ierr))

   PetscCallA(VecDestroy(U, ierr))
   PetscCallA(VecDestroy(U0, ierr))
   PetscCallA(VecDestroy(V, ierr))

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
end program TestConstraintIO2
