Program TestNewSieve
#include "finclude/petscdef.h"

	Use petsc
	Implicit NONE
	
	Type(DM)                                     :: mesh
	Character(len=256)                           :: IOBuffer,infile,outfile
	Type(PetscViewer)                            :: meshViewer
	PetscBool                                    :: flg
	PetscErrorCode                               :: ierr
	PetscReal, Dimension(:,:), Pointer           :: Coord
	PetscInt, Dimension(:,:), Pointer            :: connect
	
	Call PetscInitialize(PETSC_NULL_CHARACTER,ierr);CHKERRQ(ierr)

	Call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-i',infile,flg,ierr);CHKERRQ(ierr)
	 	
#ifdef PETSC_HAVE_EXODUSII
	Call DMMeshCreateExodus(PETSC_COMM_WORLD,infile,mesh,ierr);CHKERRQ(ierr)
	
        Call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-o',outfile,flg,ierr);CHKERRQ(ierr)
	If (flg) Then
		Call PetscViewerBinaryOpen(PETSC_COMM_WORLD,outfile,FILE_MODE_WRITE,meshViewer,ierr);CHKERRQ(ierr)
		Call DMView(mesh,meshViewer,ierr); CHKERRQ(ierr)
		Call PetscViewerDestroy(meshViewer,ierr);CHKERRQ(ierr)
	End If
#else
	Call PetscViewerBinaryOpen(PETSC_COMM_WORLD,infile,FILE_MODE_READ,meshViewer,ierr);CHKERRQ(ierr)
	Call DMMeshLoad(mesh,meshViewer,ierr); CHKERRQ(ierr)
	Call PetscViewerDestroy(meshViewer,ierr);CHKERRQ(ierr)
#endif
	
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
