Program TestMeshTopologyGetInfo

#include "finclude/petscdef.h"

   Use m_MEF90
   Use petsc

   Implicit NONE   
   

!!!startregion VARIABLES

   Character(len=MEF90_MXSTRLEN)                :: prefix, IOBuffer, myfilename
   Type(MeshTopology_Type)                      :: MeshTopology
   Type(EXO_Type)                               :: EXO
   Type(PetscViewer)                            :: myviewer
   PetscInt                                     :: verbose = 0
   PetscInt                                     :: iErr
   PetscBool                                    :: HasPrefix
   PetscInt                                     :: iBlock




   Call MEF90_Initialize()
   
   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-verbose', verbose, HasPrefix, iErr)    
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p', prefix, HasPrefix, iErr); CHKERRQ(iErr)
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No mesh prefix given\n", iErr)
      Call MEF90_Finalize()
      STOP
   End If
   
   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'

   Call DMMeshCreateExodusNG(PETSC_COMM_WORLD, EXO%filename, MeshTopology%mesh, MeshTopology%meshFS,ierr); CHKERRQ(iErr)
   
   If (verbose > 0) Then
      Write(IOBuffer, *) "Done reading and partitioning the mesh\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   write(myfilename,100) MEF90_MyRank
   Write(IOBuffer,101) trim(myfilename)
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, IOBuffer,ierr);CHKERRQ(ierr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)
   Call PetscViewerASCIIOpen(PETSC_COMM_SELF, myfilename, myviewer, iErr); CHKERRQ(iErr); 
   Call MeshTopologyGetInfo(MeshTopology, PETSC_COMM_WORLD)

   Call MeshTopologyView(MeshTopology, myviewer)
   Write(IOBuffer, *) "MeshTopology./.num_elem_blk,s ", MeshTopology%num_elem_blks, "\n"
   Call PetscViewerASCIIPrintf(myviewer, IOBuffer, iErr); CHKERRQ(iErr)
   Write(IOBuffer, *) "MeshTopology./.num_side_sets, ", MeshTopology%num_side_sets, "\n"
   Call PetscViewerASCIIPrintf(myviewer, IOBuffer, iErr); CHKERRQ(iErr)
   Write(IOBuffer, *) "MeshTopology./.num_node_sets, ", MeshTopology%num_node_sets, "\n"
   Call PetscViewerASCIIPrintf(myviewer, IOBuffer, iErr); CHKERRQ(iErr)

   Call MeshTopologyDestroy(MeshTopology)
   
   Call PetscViewerDestroy(myviewer,ierr);CHKERRQ(ierr)

   Call MEF90_Finalize()
100 Format("TestMeshTopologyGetInfo-",I0.4,".txt")
101 Format("output in file ",A,"\n")
End Program TestMeshTopologyGetInfo
