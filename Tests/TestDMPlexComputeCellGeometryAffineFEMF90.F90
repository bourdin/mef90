program TestDMPlexComputeCellGeometryAffineFEMF90
#include <petsc/finclude/petsc.h>
   use petsc
   implicit none(type, external)

   type(tDM)                      :: dm
   character(len=2048)            :: filename, IOBuffer
   type(tVec)                     :: coordLoc
   PetscReal, dimension(:), pointer :: v0, J, invJ; 
   PetscReal                      :: detJ
   PetscInt                       :: dim, i
   PetscErrorCode                 :: ierr
   PetscBool                      :: flg
   PetscInt                       :: c = 0

   PetscCallA(PetscInitialize(ierr))
   PetscCallA(PetscOptionsGetString(PETSC_NULL_OPTIONS, '', '-i', filename, flg, ierr))
   write (IOBuffer, '("Filename :",A,"\n")') trim(filename)
   PetscCallA(PetscPrintf(PETSC_COMM_WORLD, IOBuffer, ierr))
   PetscCallA(DMPlexCreateFromFile(PETSC_COMM_WORLD, filename, PETSC_NULL_CHARACTER, PETSC_FALSE, dm, ierr))

   PetscCallA(DMGetCoordinatesLocal(dm, coordLoc, ierr))
   PetscCallA(VecView(coordLoc, PETSC_VIEWER_STDOUT_WORLD, ierr))

   PetscCallA(DMGetDimension(dm, dim, ierr))
   allocate (v0(dim))
   allocate (J(dim * dim))
   allocate (invJ(dim * dim))
   PetscCallA(DMPlexComputeCellGeometryAffineFEM(dm, c, v0, J, invJ, detJ, ierr))
   write (IOBuffer, *) "V0: ", v0, "\n"
   PetscCallA(PetscPrintf(PETSC_COMM_WORLD, IOBuffer, ierr))
   PetscCallA(PetscPrintf(PETSC_COMM_WORLD, "J\n", ierr))
   do i = 0, dim - 1
      write (IOBuffer, *) J(i * dim + 1:(i + 1) * dim)
      PetscCallA(PetscPrintf(PETSC_COMM_WORLD, trim(IOBuffer)//"\n", ierr))
   end do
   PetscCallA(PetscPrintf(PETSC_COMM_WORLD, "invJ\n", ierr))
   do i = 0, dim - 1
      write (IOBuffer, *) invJ(i * dim + 1:(i + 1) * dim)
      PetscCallA(PetscPrintf(PETSC_COMM_WORLD, trim(IOBuffer)//"\n", ierr))
   end do
   write (IOBuffer, '("detJ: ",ES12.5,"\n")') detJ
   PetscCallA(PetscPrintf(PETSC_COMM_WORLD, IOBuffer, ierr))

   deallocate (v0)
   deallocate (J)
   deallocate (invJ)

   PetscCallA(DMDestroy(dm, ierr))
   PetscCallA(PetscFinalize(ierr))
end program TestDMPlexComputeCellGeometryAffineFEMF90
