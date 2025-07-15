program TestLog
#include "petsc/finclude/petsc.h"
   use m_MEF90
   use petsc

   implicit none(type, external)

   type(tPetscViewer)                   :: logViewer
   PetscErrorCode                       :: ierr

   PetscCallA(PetscInitialize(ierr))
   PetscCallA(MEF90Initialize(PETSC_COMM_WORLD, ierr))

   PetscCallA(PetscLogDefaultBegin(ierr))
   PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD, 'test.log', logViewer, ierr))
   PetscCallA(PetscLogView(logViewer, ierr))
   PetscCallA(PetscViewerDestroy(logViewer, ierr))
   PetscCallA(MEF90Finalize(ierr))
   PetscCallA(PetscFinalize(ierr))
end program TestLog
