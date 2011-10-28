Program Partitioner

#include "finclude/petscdef.h"

   Use m_MEF90
   Use petsc

   Implicit NONE   

   Type(MeshTopology_Type)                      :: MeshTopology
   Type(DM)                                     :: Tmp_Mesh
   Type(EXO_Type)                               :: EXO, MyEXO
   
   PetscBool                                    :: HasPrefix, flg
   PetscInt                                     :: verbose
   PetscErrorCode                               :: iErr, iBlk
   Character(len=256)                           :: CharBuffer, IOBuffer, filename
   Character(len=256)                           :: prefix
   Type(PetscViewer)                            :: LogViewer, MyLogViewer, MeshViewer
   Character(len=MEF90_MXSTRLEN), Dimension(4)  :: stagename
   PetscInt                                     :: iDebug
   
   
   
   Call MEF90_Initialize()


   Call PetscMemorySetGetMaximumUsage(iErr); CHKERRQ(iErr)
   iDebug = 0
   Write(stagename(1), "(A)") "Everything"
   Write(stagename(2), "(A)") "MeshCreate"
   Write(stagename(3), "(A)") "MeshDistribute"
   Write(stagename(4), "(A)") "MEF90 stuff"
   Call ALEStagePush(stagename(1), iDebug, iErr); CHKERRQ(iErr)

   verbose = 0



   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-verbose', verbose, flg, iErr)    
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p', prefix, HasPrefix, iErr)    
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No input file prefix given\n", iErr)
      Call MEF90_Finalize()
      STOP
   End If

   If (verbose > 1) Then
      Write(filename, 101) Trim(prefix), MEF90_MyRank
      Call PetscViewerASCIIOpen(PETSC_COMM_SELF, filename, MyLogViewer, iErr); CHKERRQ(iErr);   
      Write(IOBuffer, 102) MEF90_MyRank, Trim(filename)
      Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call PetscSynchronizedFlush (PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)

      Write(filename, 103) Trim(prefix)
      Call PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, LogViewer, iErr); CHKERRQ(iErr);   
      Write(IOBuffer, 104) Trim(filename)
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
101 Format(A, '-', I4.4, '.log')
102 Format('Output from processor ', I4.4, ' redirected to file ', A, '\n')
103 Format(A,'.log')
104 Format('Collective output redirected to file ', A, '\n')

   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'


   If (MEF90_NumProcs == 1) Then
      If (verbose > 0) Then
         Write(IOBuffer, *) "Reading the mesh\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Call ALEStagePush(stagename(2), iDebug, iErr); CHKERRQ(iErr)
      End If
      Call DMMeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, MeshTopology%mesh, ierr); CHKERRQ(iErr)
      If (verbose > 0) Then
         Call ALEStagePrintMemory(stagename(2), iErr); CHKERRQ(iErr)
         Call ALEStagePop(iDebug, iErr); CHKERRQ(iErr)
         Call PetscMemoryShowUsage(PETSC_VIEWER_STDOUT_WORLD, "PetscMemoryShowUsage AfterMeshCreate: \n", iErr); CHKERRQ(iErr)
         Write(IOBuffer, *) "\n\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
   Else
      If (verbose > 0) Then
         Write(IOBuffer, *) "Calling DMMeshCreateExodus\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Call ALEStagePush(stagename(2), iDebug, iErr); CHKERRQ(iErr)
      End If
      Call DMMeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
      If (verbose > 0) Then
         Call ALEStagePrintMemory(stagename(2), iErr); CHKERRQ(iErr)
         Call ALEStagePop(iDebug, iErr); CHKERRQ(iErr)
         Call PetscMemoryShowUsage(PETSC_VIEWER_STDOUT_WORLD, "PetscMemoryShowUsage AfterMeshCreate: ", iErr); CHKERRQ(iErr)
         Write(IOBuffer, *) "\n\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If

      If (verbose > 0) Then
         Write(IOBuffer, *) "Calling MeshDistribute\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Call ALEStagePush(stagename(3), iDebug, iErr); CHKERRQ(iErr)
      End If
      Call DMMeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, MeshTopology%mesh, ierr); CHKERRQ(iErr)
      If (verbose > 0) Then
         Call ALEStagePrintMemory(stagename(3), iErr); CHKERRQ(iErr)
         Call ALEStagePop(iDebug, iErr); CHKERRQ(iErr)
         Call PetscMemoryShowUsage(PETSC_VIEWER_STDOUT_WORLD, "PetscMemoryShowUsage After MeshCreate: ", iErr); CHKERRQ(iErr)
         Write(IOBuffer, *) "\n\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If

      If (verbose > 0) Then
         Write(IOBuffer, *) "Calling MeshDestroy\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      Call DMDestroy(Tmp_mesh, ierr); CHKERRQ(iErr)
      If (verbose > 0) Then
         Call ALEStagePrintMemory(stagename(1), iErr); CHKERRQ(iErr)
         Call PetscMemoryShowUsage(PETSC_VIEWER_STDOUT_WORLD, "PetscMemoryShowUsage After MeshDestroy: ", iErr); CHKERRQ(iErr)
         Write(IOBuffer, *) "\n\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
   End If
   
   If (verbose > 0) Then
      Write(IOBuffer, *) "Initializing MeshTopology object\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call ALEStagePush(stagename(4), iDebug, iErr); CHKERRQ(iErr)
   End If
   Call MeshTopologyGetInfo(MeshTopology, PETSC_COMM_WORLD)
   
   If (verbose > 0) Then
      Write(IOBuffer, *) "Initializing Element types\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   MeshTopology%Elem_Blk%Elem_Type    = MEF90_P1_Lagrange
   Do iBlk = 1, MeshTopology%Num_Elem_Blks
      Call Init_Elem_Blk_Type(MeshTopology%Elem_Blk(iBlk), MeshTopology%num_dim)
   End Do
   If (verbose > 0) Then
      Call ALEStagePrintMemory(stagename(4), iErr); CHKERRQ(iErr)
      Call PetscMemoryShowUsage(PETSC_VIEWER_STDOUT_WORLD, "PetscMemoryShowUsage After MEF90 Initializations: ", iErr); CHKERRQ(iErr)
      Write(IOBuffer, *) "\n\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   MyEXO%comm = PETSC_COMM_SELF
   MyEXO%exoid = EXO%exoid
   Write(MyEXO%filename, 200) trim(prefix), MEF90_MyRank
200 Format(A, '-', I4.4, '.gen')
 
   MyEXO%title = trim(EXO%title)
   MyEXO%Num_QA = EXO%Num_QA
   
   If (verbose > 0) Then
      Write(IOBuffer, *) "Writing EXO files\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   Call MeshTopologyWriteGlobal(MeshTopology, MyEXO, PETSC_COMM_WORLD)
   Call Write_EXO_Case(prefix, '%0.4d', MEF90_NumProcs)
   
   !!! Print mesh
   !Call MeshView(MeshTopology%Mesh, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   !!! Flushing the parallel mesh onto the disk
   !filename = trim(prefix) // '.mesh'
   !Call PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, MeshViewer, iErr); CHKERRQ(iErr)
   !Call MeshView(MeshTopology%Mesh, MeshViewer, iErr); CHKERRQ(iErr)
   !Call PetscViewerDestroy(MeshViewer, iErr); CHKERRQ(iErr)
   
   Call MeshTopologyDestroy(MeshTopology)
   If (verbose > 0) Then
      Write(IOBuffer, *) "Cleaning up\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call ALEStagePop(iDebug, iErr); CHKERRQ(iErr)
      Call ALEStagePrintMemory(stagename(1), iErr); CHKERRQ(iErr)
      Call PetscMemoryShowUsage(PETSC_VIEWER_STDOUT_WORLD, "PetscMemoryShowUsage After MeshCreate: ", iErr); CHKERRQ(iErr)
      Write(IOBuffer, *) "\n\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   If (verbose>1) Then
      Call PetscViewerFlush(MyLogViewer, iErr); CHKERRQ(iErr)
      Call PetscViewerDestroy(MyLogViewer, iErr); CHKERRQ(iErr)
      Call PetscViewerFlush(LogViewer, iErr); CHKERRQ(iErr)
      Call PetscViewerDestroy(LogViewer, iErr); CHKERRQ(iErr)
   End If
   
300 Format(A)   
   Call MEF90_Finalize()
End Program Partitioner
