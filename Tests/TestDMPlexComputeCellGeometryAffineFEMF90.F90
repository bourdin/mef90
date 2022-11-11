Program  TestDMPlexComputeCellGeometryAffineFEMF90
#include <petsc/finclude/petsc.h>
  Use petsc
  Implicit NONE   
  
  Type(tDM)                      :: dm
  Character(len=2048)            :: filename,IOBuffer
  Type(tVec)                     :: coordLoc
  PetscReal,Dimension(:),Pointer :: v0,J,invJ;
  PetscReal                      :: detJ
  PetscInt                       :: dim,i
  PetscErrorCode                 :: ierr
  PetscBool                      :: flg
  PetscInt                       :: c = 0


  PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
  PetscCallA(PetscOptionsGetString(PETSC_NULL_OPTIONS,'','-i',filename,flg,ierr))
  Write(IOBuffer,'("Filename :",A,"\n")') trim(filename)
  PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
  PetscCallA(DMPlexCreateFromFile(PETSC_COMM_WORLD,filename,PETSC_NULL_CHARACTER,PETSC_FALSE,dm,ierr))

  PetscCallA(DMGetCoordinatesLocal(dm,coordLoc,ierr))
  PetscCallA(VecView(coordLoc,PETSC_VIEWER_STDOUT_WORLD,ierr))

  PetscCallA(DMGetDimension(dm,dim,ierr))
  Allocate(v0(dim))
  Allocate(J(dim*dim))
  Allocate(invJ(dim*dim))
  PetscCallA(DMPlexComputeCellGeometryAffineFEM(dm,c,v0,J,invJ,detJ,ierr))
  Write(IOBuffer,*) "V0: ", v0,"\n"
  PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
  PetscCallA(PetscPrintf(PETSC_COMM_WORLD,"J\n",ierr))
  Do i = 0, dim-1
    Write(IOBuffer,*) J(i*dim+1:(i+1)*dim)
    PetscCallA(PetscPrintf(PETSC_COMM_WORLD,trim(IOBuffer)//"\n",ierr))
  End Do
  PetscCallA(PetscPrintf(PETSC_COMM_WORLD,"invJ\n",ierr))
  Do i = 0, dim-1
    Write(IOBuffer,*) invJ(i*dim+1:(i+1)*dim)
    PetscCallA(PetscPrintf(PETSC_COMM_WORLD,trim(IOBuffer)//"\n",ierr))
  End Do
  Write(IOBuffer,'("detJ: ",ES12.5,"\n")') detJ
  PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))

  DeAllocate(v0)
  DeAllocate(J)
  DeAllocate(invJ)

  PetscCallA(DMDestroy(dm,ierr))
  PetscCallA(PetscFinalize(ierr))
End program TestDMPlexComputeCellGeometryAffineFEMF90