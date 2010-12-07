Program TestLabel

#include "finclude/petscdef.h"

   Use petsc

   Implicit NONE   

   Type(Mesh)                                   :: Seq_Mesh,Dist_Mesh
   
   PetscBool                                    :: HasPrefix,flg
   PetscErrorCode                               :: iErr
   Character(len=256)                           :: prefix,filename, buffer, IOBuffer
   Type(PetscViewer)                            :: MeshViewer
   Integer                                      :: numproc, rank
   PetscInt                                     :: numids
   PetscInt, Dimension(:), Pointer              :: ids
   
   
   Call PetscInitialize(PETSC_NULL_CHARACTER,iErr); CHKERRQ(iErr)
   Call MPI_Comm_size(PETSC_COMM_WORLD,numproc,iErr); CHKERRQ(iErr)
   Call MPI_Comm_rank(PETSC_COMM_WORLD,rank,iErr); CHKERRQ(iErr)
   
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-p',prefix,HasPrefix,ierr);CHKERRQ(ierr);
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD,"No input file prefix given\n",ierr);CHKERRQ(ierr);
      Call PetscFinalize()
      STOP
   End If

   !!! Read Mesh
   filename = Trim(prefix)//'.gen'
   Call MeshCreateExodus(PETSC_COMM_WORLD,filename,Seq_Mesh,ierr);CHKERRQ(ierr)

   buffer = 'CellBlocks'
   Call MeshGetLabelSize(Seq_Mesh, buffer, numIds, ierr); CHKERRQ(ierr)  
   Allocate(ids(numids))
   Call MeshGetLabelIds(Seq_Mesh, buffer, ids, ierr); CHKERRQ(ierr)

   Write(IOBuffer,*) 'Label CellBlocks size:', numids, '\n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(ierr)
   Write(IOBuffer,*) 'ids: ',ids, '\n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(ierr)
   DeAllocate(ids)



   !!! Distribute mesh if numproc > 1
   If (numproc > 1) Then
      Call MeshDistribute(Seq_Mesh,PETSC_NULL_CHARACTER,Dist_Mesh,ierr);CHKERRQ(ierr)
      Call MeshDestroy(Seq_Mesh,ierr); CHKERRQ(ierr)

      Call MeshGetLabelSize(Dist_Mesh, buffer, numIds, ierr); CHKERRQ(ierr)  
      Allocate(ids(numids))
      Call MeshGetLabelIds(Dist_Mesh, buffer, ids, ierr); CHKERRQ(ierr)

      Write(IOBuffer,*) rank, ' Label CellBlocks size:', numids, '\n'
      Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(ierr)
      Call PetscSynchronizedFlush(PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      Write(IOBuffer,*) rank, ' ids: ',ids, '\n'
      Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(ierr)
      Call PetscSynchronizedFlush(PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      DeAllocate(ids)


   End If

   If (numproc > 1) Then
      Call MeshDestroy(Dist_Mesh,ierr);CHKERRQ(ierr);
   Else
      Call MeshDestroy(Seq_Mesh,ierr);CHKERRQ(ierr);
   End If

   Call PetscFinalize()
End Program TestLabel