Program Partitioner

#include "finclude/petscdef.h"

   Use m_MEF90
   Use petsc
   Use m_TestCorruption

   Implicit NONE   

   Type(MeshTopology_Type)                      :: MeshTopology
   Type(Mesh)                                   :: Tmp_Mesh
   Type(EXO_Type)                               :: EXO
   
   PetscBool                                    :: HasPrefix, flg
   PetscErrorCode                               :: ierr
   Character(len=256)                           :: CharBuffer,IOBuffer, filename
   Character(len=256)                           :: prefix
   PetscInt                                     :: numIds
   PetscInt, Dimension(:), Pointer              :: blkIds
   
   
   
   Call MEF90_Initialize()
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p', prefix, HasPrefix, ierr)    
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No input file prefix given\n", ierr)
      Call MEF90_Finalize()
      STOP
   End If

   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'


   Write(IOBuffer, *) "Reading the mesh\n"
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, ierr); CHKERRQ(ierr)
   If (MEF90_NumProcs == 1) Then
      Call MeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, MeshTopology%mesh, ierr); CHKERRQ(ierr)
   Else
      Call MeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, Tmp_mesh, ierr); CHKERRQ(ierr)
      Call MeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, MeshTopology%mesh, ierr); CHKERRQ(ierr)
      Call MeshDestroy(Tmp_mesh, ierr); CHKERRQ(ierr)
   End If
   
   CHKMEMQ
   !Call MeshTopologyReadEXO(MeshTopology, EXO)
   !!! trying step by step:
   Call MeshExodusGetInfo(MeshTopology%mesh, MeshTopology%Num_Dim, MeshTopology%Num_Verts, MeshTopology%Num_Elems, MeshTopology%Num_Elem_Blks, MeshTopology%Num_Node_Sets, ierr); CHKERRQ(ierr)
   Write(*,*) 'MeshTopology%Num_Dim            ',MeshTopology%Num_Dim
   Write(*,*) 'MeshTopology%Num_Verts          ',MeshTopology%Num_Verts
   Write(*,*) 'MeshTopology%Num_Elems          ',MeshTopology%Num_Elems
   Write(*,*) 'MeshTopology%Num_Elem_Blk       ',MeshTopology%Num_Elem_Blks
   Write(*,*) 'MeshTopology%Num_Node_Sets      ',MeshTopology%Num_Node_Sets
   !!! Extracts sizes from the Mesh oject

   ! Read Elem block information
   CharBuffer = 'CellBlocks'
   Call MeshGetLabelSize(MeshTopology%mesh, CharBuffer, numIds, ierr); CHKERRQ(ierr)  
   Write(*,*) 'numIds =                        ', numIds
   Allocate(blkIds(numIds))
   Call MeshGetLabelIds(MeshTopology%mesh, CharBuffer, blkIds, ierr); CHKERRQ(ierr)
   Write(*,*) 'blkIds                          ', blkIds
   
   Call MyMeshTopologyReadEXO(MeshTopology, EXO)
   Call ModMeshTopologyReadEXO(MeshTopology, EXO)

   Write(*,*) '============ MEF90 MeshTopologyReadEXO ============'
   Call MeshTopologyReadEXO(MeshTopology, EXO)
   Write(*,*) '============ MEF90 MeshTopologyReadEXO ============'
   
   Call MeshTopologyDestroy(MeshTopology)
   Call MEF90_Finalize()

Contains
   Subroutine MyMeshTopologyReadEXO(dMeshTopology, dEXO)
      !!! Remove the element and coordinate stuff and move in separate functions
      Type (MeshTopology_Type)                     :: dMeshTopology
      Type (EXO_Type)                              :: dEXO
      PetscInt                                     :: iErr, iBlk, iSet
      Character(len=256)                           :: CharBuffer
      
      PetscReal, Dimension(:,:), Pointer           :: array
      PetscInt, Dimension(:,:), Pointer            :: arrayCon
      PetscInt                                     :: embedDim
      PetscInt                                     :: iE, iElem, numIds, blkId, setId
      PetscInt, Dimension(:), Pointer              :: setIds
      PetscInt, Dimension(:), Pointer              :: blkIds
      PetscInt, Dimension(:), Pointer              :: Tmp_ID, Tmp_GlobalID


      
      ! Read Global Geometric Parameters
      Call MeshExodusGetInfo(dMeshTopology%mesh, dMeshTopology%Num_Dim, dMeshTopology%Num_Verts, dMeshTopology%Num_Elems, dMeshTopology%Num_Elem_Blks, dMeshTopology%Num_Node_Sets, iErr); CHKERRQ(iErr)
      Write(*,*) '============ MyMeshTopologyReadEXO ============'
      Write(*,*) 'dMeshTopology%Num_Dim            ',dMeshTopology%Num_Dim
      Write(*,*) 'dMeshTopology%Num_Verts          ',dMeshTopology%Num_Verts
      Write(*,*) 'dMeshTopology%Num_Elems          ',dMeshTopology%Num_Elems
      Write(*,*) 'dMeshTopology%Num_Elem_Blk       ',dMeshTopology%Num_Elem_Blks
      Write(*,*) 'dMeshTopology%Num_Node_Sets      ',dMeshTopology%Num_Node_Sets
      !!! Extracts sizes from the Mesh oject

      ! Read Elem block information
      CharBuffer = 'CellBlocks'
      Call MeshGetLabelSize(dMeshTopology%mesh, CharBuffer, numIds, ierr); CHKERRQ(ierr)  
      !!! Get the number of labels of type 'CellBlocks' in the mesh
      If (numIds .ne. dMeshTopology%Num_Elem_blks) Then
         SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, 'Invalid number of element ids', ierr)
      End If
      !!! Compare to the number initialized in MeshTopology
      
      Allocate(blkIds(numIds))
      Write(*,*) 'numIds =                        ', numIds
      
      Call MeshGetLabelIds(dMeshTopology%mesh, CharBuffer, blkIds, ierr); CHKERRQ(ierr)
      Write(*,*) 'blkIds                          ', blkIds
      
      Allocate(dMeshTopology%Elem_blk(dMeshTopology%Num_Elem_blks))
      If (dMeshTopology%Num_Elem_blks > 0) Then
         Do iBlk = 1, dMeshTopology%Num_Elem_Blks
            blkId = blkIds(iBlk)
            dMeshTopology%Elem_blk(iBlk)%ID = blkId
            Call MeshGetStratumSize(dMeshTopology%mesh, CharBuffer, blkId, dMeshTopology%elem_blk(iBlk)%Num_Elems, ierr); CHKERRQ(iErr)
            !!! Get the size of the layer (stratum) 'CellBlock' of Mesh
            Allocate(dMeshTopology%Elem_blk(iBlk)%Elem_ID(dMeshTopology%elem_blk(iBlk)%Num_Elems))
            Call MeshGetStratum(dMeshTopology%mesh, CharBuffer, blkId, dMeshTopology%Elem_blk(iBlk)%Elem_ID, ierr); CHKERRQ(iErr)
            !!! Get the layer (stratum) 'CellBlock' of Mesh in C numbering
            dMeshTopology%Elem_blk(iBlk)%Elem_ID = dMeshTopology%Elem_blk(iBlk)%Elem_ID + 1
            !!! Converts to Fortran style indexing
         End Do
      End If
      Deallocate(blkIds)
      Write(*,*) '============ MyMeshTopologyReadEXO ============'
   End Subroutine MyMeshTopologyReadEXO

End Program Partitioner
