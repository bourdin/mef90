! mpirun -np 6 ./TestConstraintIO3 -geometry ../TestMeshes/TallRectangle-0.1.gen -result out.exo -options_file ~/Desktop/Debug/TensionInternalBC_FS.yaml

program TestConstraintIO3
#include <petsc/finclude/petsc.h>
   use m_MEF90
   use m_MEF90_DefMech
   use petsc
   implicit none(type, external)

   PetscErrorCode                                      :: ierr
   type(MEF90Ctx_Type), target                          :: MEF90Ctx
   type(MEF90CtxGlobalOptions_Type)                    :: MEF90GlobalOptions_default
   type(MEF90CtxGlobalOptions_Type), pointer            :: MEF90GlobalOptions
   type(tDM), target                                    :: dm, dmU
   PetscBool                                           :: interpolate = PETSC_TRUE

   !!! Defect mechanics contexts
   type(MEF90DefMech_Type), target                     :: MEF90DefMechCtx
   type(MEF90DefMechGlobalOptions_Type)                :: MEF90DefMechGlobalOptions

   PetscInt                                            :: numNodalVar = 3, numCellVar = 0, numGVar = 0
   character(len=MEF90MXSTRLEN), dimension(:), pointer   :: nodalVarName, cellVarName, gVarName
   character(len=MEF90MXSTRLEN)                        :: filename
   PetscInt                                            :: dim
   type(tPetscSection)                                 :: sectionU
   type(tVec)                                          :: U, V, locVecV, tmpVec
   PetscReal, dimension(:), pointer                      :: time
   PetscInt                                            :: i, step = 1_ki
   PetscReal, dimension(:), pointer                      :: DisplacementArray, locVecVArray
   ! PetscReal                                           :: myerr,err

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

   ! Create DM from file
   PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm, MEF90Ctx%geometryfile, PETSC_NULL_CHARACTER, interpolate, dm, ierr))
   PetscCallA(DMPlexDistributeSetDefault(dm, PETSC_FALSE, ierr))
   PetscCallA(DMSetUseNatural(dm, PETSC_TRUE, ierr))
   PetscCallA(DMSetFromOptions(dm, ierr))
   PetscCallA(DMViewFromOptions(dm, PETSC_NULL_OBJECT, "-mef90dm_view", ierr))

   PetscCallA(DMGetDimension(dm, dim, ierr))

   ! Open exodus file + write geometry + format the file
   numNodalVar = dim
   numCellVar = 0
   numGVar = 0

   allocate (nodalVarName(numNodalVar))
   allocate (cellVarName(numCellVar))
   allocate (gVarName(numGVar))
   if (dim == 2) then
      nodalVarName = ["Displacement_X", "Displacement_Y"]
   else
      nodalVarName = ["Displacement_X", "Displacement_Y", "Displacement_Z"]
   end if

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
      type(tPetscSF)                      :: naturalPointSF
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

   PetscCallA(MEF90DefMechCreate(MEF90DefMechCtx, dm, MEF90Ctx, "", ierr))
   PetscCallA(MEF90DefMechCtx%setFromOptions(ierr))
   PetscCallA(MEF90DefMechGlobalOptionsSetFromOptions(MEF90DefMechCtx%comm, trim(MEF90DefMechCtx%prefix), MEF90DefMechGlobalOptions, ierr))
   PetscCallA(DMDestroy(dm, ierr))

   PetscCallA(VecGetDM(MEF90DefMechCtx%displacementLocal, dmU, ierr))
   PetscCallA(DMGetLocalSection(dmU, sectionU, ierr))
   PetscCallA(DMGetGlobalVector(dmU, U, ierr))
   PetscCallA(DMGetGlobalVector(dmU, V, ierr))

   ! Initialize boundary values of MEF90DefMechCtx%displacementLocal with values from the command line
   PetscCallA(VecSet(MEF90DefMechCtx%displacementLocal, 0.0_kr, ierr))
   PetscCallA(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%displacementLocal, 1.0_kr, ierr))

   ! Save MEF90DefMechCtx%displacementLocal in exo file
   PetscCallA(MEF90EXOVecView(MEF90DefMechCtx%displacementLocal, MEF90DefMechCtx%displacementToIOSF, MEF90DefMechCtx%IOTodisplacementSF, MEF90Ctx%resultViewer, step, dim, ierr))

   ! ! Create new vectors and read them from the file
   PetscCallA(VecDuplicate(MEF90DefMechCtx%displacementLocal, locVecV, ierr))
   PetscCallA(VecSet(locVecV, 0.0_kr, ierr))
   PetscCallA(PetscObjectSetName(locVecV, "Displacement", ierr))

   ! PetscCallA(DMGetLocalVector(dmDisplacement,tmpVec,ierr))
   PetscCallA(VecDuplicate(MEF90DefMechCtx%displacementLocal, tmpVec, ierr))
   PetscCallA(PetscObjectSetName(tmpVec, "Displacement", ierr))
   PetscCallA(MEF90EXOVecLoad(tmpVec, MEF90DefMechCtx%displacementToIOSF, MEF90DefMechCtx%IOToDisplacementSF, MEF90DefMechCtx%MEF90Ctx%resultViewer, step, MEF90DefMechCtx%dim, ierr))

   ! PetscCallA(MEF90VecCopySF(tmpVec,locVecV,MEF90DefMechCtx%displacementConstraintsSF,ierr))
   PetscCallA(VecCopy(tmpVec, locVecV, ierr))
        !!! VecCopy and MEF90VecCopySF give the same error so the problem is in the reading
        !!! I assume that the issue is in copying from IO to local?

   PetscCallA(VecDestroy(tmpVec, ierr))

   ! PetscCallA(MEF90EXOVecLoad(locVecV,MEF90DefMechCtx%displacementToIOSF,MEF90DefMechCtx%IOTodisplacementSF,MEF90Ctx%resultViewer,step,dim,ierr))

   write (filename, '("out-",I4.4,".txt")') MEF90Ctx%rank
   open (file=filename, unit=99)
   write (*, *) 'Opening ', filename

   PetscCallA(VecGetArrayRead(MEF90DefMechCtx%displacementLocal, DisplacementArray, ierr))
   PetscCallA(VecGetArrayRead(locVecV, locVecVArray, ierr))

   write (99, *) 'Max diff', maxval(DisplacementArray - locVecVArray), maxloc(DisplacementArray - locVecVArray)
   do i = 1, size(locVecVArray)
      if (abs(DisplacementArray(i) - locVecVArray(i)) > 1.e-5) then
         write (99, *) i, DisplacementArray(i), locVecVArray(i)
      end if
   end do
   PetscCallA(VecRestoreArrayRead(MEF90DefMechCtx%displacementLocal, DisplacementArray, ierr))
   PetscCallA(VecRestoreArrayRead(locVecV, locVecVArray, ierr))

   PetscCallA(VecSet(U, -1.0_kr, ierr))
   PetscCallA(VecSet(V, 1.0_kr, ierr))
   PetscCallA(DMLocalToGlobal(dmU, MEF90DefMechCtx%displacementLocal, INSERT_VALUES, U, ierr))
   PetscCallA(DMLocalToGlobal(dmU, locVecV, INSERT_VALUES, V, ierr))
   PetscCallA(VecGetArrayRead(U, DisplacementArray, ierr))
   PetscCallA(VecGetArrayRead(V, locVecVArray, ierr))

   write (99, *) 'Max diff', maxval(DisplacementArray - locVecVArray), maxloc(DisplacementArray - locVecVArray)
   do i = 1, size(locVecVArray)
      if (abs(DisplacementArray(i) - locVecVArray(i)) > 1.e-5) then
         write (99, *) i, DisplacementArray(i), locVecVArray(i)
      end if
   end do
   PetscCallA(VecRestoreArrayRead(U, DisplacementArray, ierr))
   PetscCallA(VecRestoreArrayRead(V, locVecVArray, ierr))
   close (99)

   ! Cleanup Vec
   PetscCallA(DMRestoreGlobalVector(dmU, U, ierr))
   PetscCallA(DMRestoreGlobalVector(dmU, V, ierr))

   ! Cleanup DMs
   deallocate (time)

   ! Exit nicely
   PetscCallA(MEF90CtxDestroy(MEF90Ctx, ierr))
   PetscCallA(MEF90Finalize(ierr))
   PetscCallA(PetscFinalize(ierr))
end program TestConstraintIO3
