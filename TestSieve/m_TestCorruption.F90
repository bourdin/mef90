Module m_TestCorruption

#include "finclude/petscdef.h"

   Use m_MEF90
   Use petsc

   Implicit NONE   
   Public :: ModMeshTopologyReadEXO 

Contains
   Subroutine ModMeshTopologyReadEXO(dMeshTopology, dEXO)
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
      Write(*,*) '============ ModMeshTopologyReadEXO ============'
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
      Write(*,*) '============ ModMeshTopologyReadEXO ============'
   End Subroutine ModMeshTopologyReadEXO


End Module m_TestCorruption