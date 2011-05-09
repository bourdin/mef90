Program TestNewSieve
#include "finclude/petscdef.h"

   Use petsc
   Implicit NONE
   
   Type(DM)                                     :: mesh
   Character(len=256)                           :: IOBuffer, filename
   PetscBool                                    :: flg
   PetscErrorCode                               :: ierr
   Type(SectionReal)                            :: S
   PetscReal, Dimension(:,:), Pointer           :: Coord
   PetscInt, Dimension(:,:), Pointer            :: connect
   
   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr);CHKERRQ(ierr)
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-f',filename,flg,ierr);CHKERRQ(ierr)
   Call DMMeshCreateExodus(PETSC_COMM_WORLD,filename,mesh,ierr);CHKERRQ(ierr)
   Call DMView(mesh,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)

   Call DMMeshGetVertexSectionReal(mesh,'U',1,S,iErr);CHKERRQ(iErr)
   Call SectionRealDestroy(S,ierr);CHKERRQ(ierr)
   
   Call DMMeshGetCoordinatesF90(mesh,Coord,ierr);CHKERRQ(ierr)
   Write(*,*) 'Coordinates: '
   Write(*,*) Coord
   Call DMMeshRestoreCoordinatesF90(mesh,Coord,ierr);CHKERRQ(ierr)
   
   Call DMMeshGetElementsF90(mesh,connect,ierr);CHKERRQ(ierr)
   write(*,*) 'Connect: '
   write(*,*) connect
   Call DMMeshRestoreElementsF90(mesh,connect,ierr);CHKERRQ(ierr)
   Call DMDestroy(mesh,ierr);CHKERRQ(ierr)
   Call PetscFinalize(ierr)
End Program TestNewSieve
