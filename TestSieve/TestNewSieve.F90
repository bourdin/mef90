Program TestNewSieve
#include "finclude/petscdef.h"

   Use petsc
   Implicit NONE
   
   Type(DM)                                     :: mesh
   Character(len=256)                           :: IOBuffer, filename
   PetscBool                                    :: flg
	PetscErrorCode                               :: ierr
   
   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr);CHKERRQ(ierr)
	Call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-f',filename,flg,ierr);CHKERRQ(ierr)
   Call DMMeshCreateExodus(PETSC_COMM_WORLD,filename,mesh,ierr);CHKERRQ(ierr)
   Call DMView(mesh,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
   Call DMDestroy(mesh,ierr);CHKERRQ(ierr)
   Call PetscFinalize(ierr)
End Program TestNewSieve
