Program YAMLValidator
#include "finclude/petscdef.h"
   use petsc
   IMPLICIT NONE

   PetscInt                         :: ierr
   Integer                          :: rank

   Call MPI_Init(ierr)
   Call MPI_Comm_Rank(MPI_COMM_WORLD,rank,ierr)
   if (rank == 0) then
      write(*,*) "Parsing options. If this takes more than a few seconds, "
      write(*,*) "there is probably a problem with the options file"
      end if
   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   Call PetscOptionsView(PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr);

   Call PetscFinalize()
End Program YAMLValidator