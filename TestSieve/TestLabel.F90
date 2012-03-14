Program TestLabel

#include "finclude/petscdef.h"

   Use m_mef90
   Use petsc

   Implicit NONE   

   Type(DM)                                     :: dmBody,tmpDM
   
   PetscBool                                    :: HasPrefix,flg
   PetscErrorCode                               :: ierr
   Character(len=256)                           :: prefix,filename,buffer,IOBuffer
   Type(PetscViewer)                            :: MeshViewer
   PetscInt                                     :: num_set
   PetscInt,Dimension(:),Pointer                :: set_ids
   Type(IS)                                     :: set_IS
   Integer                                      :: cpu_ws = 0
   Integer                                      :: io_ws = 0
   Integer                                      :: numCellSet,numDim
   Integer                                      :: exoid
   PetscReal                                    :: vers
   
   Call MEF90_Initialize()
   
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-p',prefix,HasPrefix,ierr);CHKERRQ(ierr);
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD,"No input file prefix given\n",ierr);CHKERRQ(ierr);
      Call PetscFinalize()
      STOP
   End If

   !!! Read Mesh
   filename = Trim(prefix)//'.gen'
   If (MEF90_Myrank == 0) Then
      exoid = EXOPEN(filename,EXREAD,cpu_ws,io_ws,vers,ierr)
   End If
   If (MEF90_numprocs == 1) Then
      Call DMMeshCreateExodusNG(PETSC_COMM_WORLD,exoid,dmBody,ierr);CHKERRQ(ierr)
   Else
      Call DMMeshCreateExodusNG(PETSC_COMM_WORLD,exoid,tmpDM,ierr);CHKERRQ(ierr)
      Call DMMeshDistribute(tmpDM,PETSC_NULL_CHARACTER,dmBody,ierr);CHKERRQ(iErr)
      Call DMDestroy(tmpDM,ierr);CHKERRQ(iErr)
   End If
   If (MEF90_Myrank == 0) Then
      Call EXCLOS(exoid,ierr)
   End If

   Call DMMeshGetLabelIdIS(dmBody,'Cell Sets',set_IS,ierr);CHKERRQ(ierr)
   Call DMMeshGetLabelSize(dmBody,'Cell Sets',num_set,ierr);CHKERRQ(ierr)
   Write(IOBuffer,*) MEF90_MyRank,': Label Cell Sets size:',num_set,'\n'
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)


   Call PetscPrintf(PETSC_COMM_WORLD,'Raw from the DM\n',ierr);CHKERRQ(ierr)

   Call ISGetLocalSize(set_IS,num_set,ierr);CHKERRQ(ierr)
   Write(IOBuffer,*) MEF90_MyRank,': label_IS local size:',num_set,'\n'
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)

   Call ISGetSize(set_IS,num_set,ierr);CHKERRQ(ierr)
   Write(IOBuffer,*) MEF90_MyRank,': label_IS size:      ',num_set,'\n'
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)

   Call ISGetIndicesF90(set_IS,set_ids,ierr);CHKERRQ(ierr)
   Write(IOBuffer,*) MEF90_MyRank,': labels              ',set_ids,'\n'
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
   Call ISRestoreIndicesF90(set_IS,set_ids,ierr);CHKERRQ(ierr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)
   
   Call PetscPrintf(PETSC_COMM_WORLD,'Synchronized version\n',ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,set_IS)
   Call ISGetLocalSize(set_IS,num_set,ierr);CHKERRQ(ierr)
   Write(IOBuffer,*) MEF90_MyRank,': label_IS local size:',num_set,'\n'
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)

   Call ISGetSize(set_IS,num_set,ierr);CHKERRQ(ierr)
   Write(IOBuffer,*) MEF90_MyRank,': label_IS size:      ',num_set,'\n'
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)

   Call ISGetIndicesF90(set_IS,set_ids,ierr);CHKERRQ(ierr)
   Write(IOBuffer,*) MEF90_MyRank,': labels              ',set_ids,'\n'
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
   Call ISRestoreIndicesF90(set_IS,set_ids,ierr);CHKERRQ(ierr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)


   Call DMDestroy(dmBody,ierr);CHKERRQ(ierr);
   Call MEF90_Finalize()
End Program TestLabel