Program TestNewSieve
#include "finclude/petscdef.h"
   Use petsc
   Implicit NONE   
#include "exodusII.inc"
  
   Type(DM)                                     :: SeqMesh,DistMesh
   Character(len=256)                           :: IOBuffer,infile,outfile
   Type(PetscViewer)                            :: meshViewer
   PetscBool                                    :: flg
   PetscErrorCode                               :: ierr
   PetscReal,Dimension(:,:),Pointer             :: Coord
   PetscInt,Dimension(:,:),Pointer              :: connect
   Integer                                      :: rank,numproc
   Integer                                      :: exoid
   Integer                                      :: cpu_ws = 0
   Integer                                      :: io_ws = 0
   PetscReal                                    :: vers
  
   Type(SectionReal)                            :: SR
   Type(SectionInt)                             :: SI
   Type(Vec)                                    :: V
   Type(Mat)                                    :: M
   Type(VecScatter)                             :: SR2V
   PetscReal                                    :: val = -12.34
   
   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr);CHKERRQ(ierr)
   Call MPI_Comm_size(PETSC_COMM_WORLD,numproc,iErr);CHKERRQ(iErr)
   Call MPI_Comm_rank(PETSC_COMM_WORLD,rank,iErr);CHKERRQ(iErr)
   
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-i',infile,flg,ierr);CHKERRQ(ierr)
   
   Write(IOBuffer,*) "Reading the mesh and distributing if necessary\n"
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr);CHKERRQ(iErr)
   
   If (rank == 0) Then
      exoid = exopen(infile,EXREAD,cpu_ws,io_ws,vers,ierr)
   End If
   If (numproc == 1) Then
      Call DMMeshCreateExodusNG(PETSC_COMM_WORLD,exoid,DistMesh,ierr);CHKERRQ(ierr)
   Else
      Call DMMeshCreateExodusNG(PETSC_COMM_WORLD,exoid,SeqMesh,ierr);CHKERRQ(ierr)
      Call DMMeshDistribute(SeqMesh,PETSC_NULL_CHARACTER,DistMesh,ierr);CHKERRQ(iErr)
      Call DMDestroy(SeqMesh,ierr);CHKERRQ(iErr)
   End If
   If (rank == 0) Then
      Call EXCLOS(exoid,ierr)
   End If
   
   
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-o',outfile,flg,ierr);CHKERRQ(ierr)
   If (flg) Then
     Call PetscViewerHDF5Open(PETSC_COMM_WORLD,outfile,FILE_MODE_WRITE,meshViewer,ierr);CHKERRQ(ierr)
     !Call DMView(DistMesh,meshViewer,ierr);CHKERRQ(ierr)
     !Call PetscViewerDestroy(meshViewer,ierr);CHKERRQ(ierr)
   End If
   
   
   Write(IOBuffer,*) 'Testing DMMeshGetCoordinatesF90...'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr);CHKERRQ(iErr)
   Call DMMeshGetCoordinatesF90(DistMesh,Coord,ierr);CHKERRQ(ierr)
   Write(100+rank,*) 'Coordinates: '
   Write(100+rank,*) Coord
   Call DMMeshRestoreCoordinatesF90(DistMesh,Coord,ierr);CHKERRQ(ierr)
   Write(IOBuffer,*) 'done.\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr);CHKERRQ(iErr)
   
   
   Write(IOBuffer,*) 'Testing DMMeshGetElementsF90...'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr);CHKERRQ(iErr)
   Call DMMeshGetElementsF90(DistMesh,connect,ierr);CHKERRQ(ierr)
   write(100+rank,*) 'Connect: '
   write(100+rank,*) connect
   Call DMMeshRestoreElementsF90(DistMesh,connect,ierr);CHKERRQ(ierr)
   Write(IOBuffer,*) 'done.\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr);CHKERRQ(iErr)
   
   Write(IOBuffer,*) 'Testing SectionReal stuff.\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr);CHKERRQ(iErr)
   
   Write(IOBuffer,*) '   SectionReal Creation.\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr);CHKERRQ(iErr)
   Call DMMeshGetVertexSectionReal(DistMesh,'Test SectionReal',1,SR,iErr);CHKERRQ(iErr)
   Call SectionRealSet(SR,val,iErr);CHKERRQ(iErr)
   Call SectionRealView(SR,PETSC_VIEWER_STDOUT_WORLD,iErr);CHKERRQ(iErr)   
   
   
   Write(IOBuffer,*) '   Vec Creation from SectionReal.\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr);CHKERRQ(iErr)
   Call DMMeshCreateVector(DistMesh,SR,V,iErr);CHKERRQ(iErr)
   
   Write(IOBuffer,*) '   Scatter from SectionReal to Vec.\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr);CHKERRQ(iErr)
   Call DMMeshCreateGlobalScatter(DistMesh,SR,SR2V,iErr);CHKERRQ(iErr)
   !Call VecScatterView(SR2V,PETSC_VIEWER_STDOUT_WORLD,iErr);CHKERRQ(iErr)
   Call SectionRealToVec(SR,SR2V,SCATTER_FORWARD,V,iErr);CHKERRQ(iErr)
   Call VecView(V,PETSC_VIEWER_STDOUT_WORLD,iErr);CHKERRQ(iErr)   
   If (flg) Then
      Call VecView(V,meshViewer,iErr);CHKERRQ(iErr)   
   End If
   
   Write(IOBuffer,*) 'Matrix creation.\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr);CHKERRQ(iErr)
   Call DMMeshCreateMatrix(DistMesh,SR,MATMPIAIJ,M,iErr);CHKERRQ(iErr)
   Call MatZeroEntries(M,iErr);CHKERRQ(iErr)
   Call MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
   Call MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
   Call MatView(M,PETSC_VIEWER_STDOUT_WORLD,iErr);CHKERRQ(iErr)
   
   
   Write(IOBuffer,*) 'Destroying objects.\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr);CHKERRQ(iErr)
   
   Call SectionRealDestroy(SR,iErr);CHKERRQ(iErr)
   Call VecDestroy(V,iErr);CHKERRQ(iErr)
   Call VecScatterDestroy(SR2V,iErr);CHKERRQ(iErr)
   Call MatDestroy(M,iErr);CHKERRQ(iErr)
   Call DMDestroy(DistMesh,iErr);CHKERRQ(iErr)
   
   Write(IOBuffer,*) 'Done.\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr);CHKERRQ(iErr)
   
   If (flg) Then
      Call PetscViewerDestroy(meshViewer,ierr);CHKERRQ(ierr)
   End If
   
   Call PetscFinalize(iErr)
End Program TestNewSieve
