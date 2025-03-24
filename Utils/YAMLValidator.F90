Program YAMLValidator
#include "petsc/finclude/petsc.h"
   use petsc
   IMPLICIT NONE

   PetscErrorCode                   :: ierr
   Integer                          :: rank

   PetscCallA(MPI_Init(ierr))
   PetscCallA(MPI_Comm_Rank(MPI_COMM_WORLD,rank,ierr))
   if (rank == 0) then
      write(*,*) "Parsing options. If this takes more than a few seconds, "
      write(*,*) "there is probably a problem with the options file"
      end if
   PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
   PetscCallA(PetscOptionsView(PETSC_NULL_OPTIONS,PETSC_VIEWER_STDOUT_WORLD,ierr))

   PetscCallA(PetscFinalize(ierr))
End Program YAMLValidator