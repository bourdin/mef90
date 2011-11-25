Program TestMeshTopologyGetInfo

#include "finclude/petscdef.h"

   Use m_MEF90
   Use petsc

   Implicit NONE   
   

!!!startregion VARIABLES

   Character(len=MEF90_MXSTRLEN)                :: prefix, IOBuffer, filename
   Type(MeshTopology_Type)                      :: MeshTopology
   Type(EXO_Type)                               :: EXO
   PetscInt                                     :: verbose = 0
   PetscInt                                     :: iErr
   PetscBool                                    :: HasPrefix
   Type(DM)                                     :: Tmp_Mesh
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
	EXO%exoid = EXOPEN(EXO%filename, EXREAD, exo_cpu_ws, exo_io_ws, PETSC_NULL_INTEGER, ierr)
	
! 	Call EXO_Check_Numbering(EXO, iErr)
! 	If (iErr /= 0) Then
! 		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, 'Unsupported numbering of the element blocks, side sets or node sets\n', iErr)
! 	End If

	If (MEF90_NumProcs == 1) Then
	   Call DMMeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, MeshTopology%mesh, ierr); CHKERRQ(iErr)
	Else    
	   Call DMMeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
	   Call DMMeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, MeshTopology%mesh, ierr); CHKERRQ(iErr)
	   Call DMDestroy(Tmp_mesh, ierr); CHKERRQ(iErr)
	End If
	
	Call MeshTopologyGetInfo(MeshTopology, PETSC_COMM_WORLD)
	If (verbose > 0) Then
	   Write(IOBuffer, *) "Done reading and partitioning the mesh\n"
	   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
	End If
	
   	Call MeshTopologyGetInfo(MeshTopology, PETSC_COMM_WORLD)
	Write(IOBuffer, *) "MeshTopology%num_elem_blk, ", MeshTopology%num_elem_blks, "\n"
	Call PetscPrintf(PETSC_COMM_SELF, IOBuffer, iErr); CHKERRQ(iErr)
	Write(IOBuffer, *) "MeshTopology%num_elem_blk_global, ", MeshTopology%num_elem_blks_global, "\n"
	Call PetscPrintf(PETSC_COMM_SELF, IOBuffer, iErr); CHKERRQ(iErr)
	Write(IOBuffer, *) "MeshTopology%num_node_sets, ", MeshTopology%num_node_sets, "\n"
	Call PetscPrintf(PETSC_COMM_SELF, IOBuffer, iErr); CHKERRQ(iErr)
	Write(IOBuffer, *) "MeshTopology%num_node_sets_global, ", MeshTopology%num_node_sets_global, "\n"
	Call PetscPrintf(PETSC_COMM_SELF, IOBuffer, iErr); CHKERRQ(iErr)

	Call MeshTopologyDestroy(MeshTopology)
	Call EXCLOS(EXO%exoid, iErr)
	EXO%exoid = 0

	Call MEF90_Finalize()

End Program TestMeshTopologyGetInfo
