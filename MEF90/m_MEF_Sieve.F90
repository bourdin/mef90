Module m_MEF_Sieve
#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmeshdef.h"

   Use m_MEF_LinAlg
   Use m_MEF_Types
   Use petsc
   Use petscmesh
   
   IMPLICIT NONE
   Private
   
   Public :: MeshTopologyReadEXO
   Public :: MeshInitElementConnectivity
   Public :: MeshInitCoordinates
   
   Interface MeshInitCoordinates
      Module Procedure MeshInitCoordinatesVect2D, MeshInitCoordinatesVect3D
   End Interface MeshInitCoordinates
   
   Interface MeshInitElementConnectivity
      Module Procedure MeshInitElementConnectivity2D_Scal, MeshInitElementConnectivity2D, MeshInitElementConnectivity2D_Elast, MeshInitElementConnectivity3D_Scal, MeshInitElementConnectivity3D, MeshInitElementConnectivity3D_Elast
   End Interface MeshInitElementConnectivity
      
Contains
   Subroutine MeshTopologyReadEXO(dMeshTopology, dEXO)
      !!! Remove the element and coordinate stuff and move in separate functions
      Type (MeshTopology_Info)                     :: dMeshTopology
      Type (EXO_Info)                              :: dEXO
      PetscInt                                     :: iErr, iBlk, iSet
      Character(len=256)                           :: CharBuffer
      
!      Type(Mesh)                                   :: Tmp_mesh
      PetscReal, Dimension(:,:), Pointer           :: array
      PetscInt, Dimension(:,:), Pointer            :: arrayCon
      PetscInt                                     :: embedDim
      PetscInt                                     :: iE, iElem, numIds, blkId, setId
      PetscInt, Dimension(:), Pointer              :: blkIds
      PetscInt, Dimension(:), Pointer              :: setIds
      
      ! Open File
!      Call MeshCreateExodus(PETSC_COMM_WORLD, dEXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
      !!! reads exo file, stores all information in a Mesh

!      Call MeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, dMeshTopology%mesh, ierr); CHKERRQ(iErr)
      !!! Partitions using a partitioner (currently PETSC_NULL_CHARACTER) 
     
!      Call MeshDestroy(Tmp_mesh, ierr); CHKERRQ(iErr)

      ! Read Global Geometric Parameters
      Call MeshExodusGetInfo(dMeshTopology%mesh, dMeshTopology%Num_Dim, dMeshTopology%Num_Verts, dMeshTopology%Num_Elems, dMeshTopology%Num_Elem_Blks, dMeshTopology%Num_Node_Sets, iErr); CHKERRQ(iErr)
      !!! Extracts sizes from the Mesh oject

      ! Read Elem block information
      CharBuffer = 'CellBlocks'
      Call MeshGetLabelSize(dMeshTopology%mesh, CharBuffer, numIds, ierr); CHKERRQ(ierr)  
      !!! Get the number of labels of type 'CellBlocks' in the mesh
      If (numIds .ne. dMeshTopology%Num_Elem_blks) Then
         SETERRQ(PETSC_ERR_ARG_SIZ, 'Invalid number of element ids', ierr)
      End If
      !!! Compare to the number initialized in MeshTopology
      
      Allocate(dMeshTopology%Elem_blk(dMeshTopology%Num_Elem_blks))
      Allocate(blkIds(numIds))
      Call MeshGetLabelIds(dMeshTopology%mesh, CharBuffer, blkIds, ierr); CHKERRQ(ierr)
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
      
      ! Read Node set information
      CharBuffer = 'VertexSets'
      Call MeshGetLabelSize(dMeshTopology%mesh, CharBuffer, numIds, ierr); CHKERRQ(ierr)
      If (numIds .ne. dMeshTopology%Num_node_sets) Then
         SETERRQ(PETSC_ERR_ARG_SIZ, 'Invalid number of node ids', ierr)
      End If
      Allocate(dMeshTopology%Node_Set(dMeshTopology%Num_Node_Sets))
      Allocate(setIds(numIds))
      Call MeshGetLabelIds(dMeshTopology%mesh, CharBuffer, setIds, ierr); CHKERRQ(ierr)
      If (dMeshTopology%Num_Node_Sets > 0) Then
         Do iSet = 1, dMeshTopology%Num_node_sets
            setId = setIds(iSet)
            dMeshTopology%Node_Set(iSet)%ID = setId
            Call MeshGetStratumSize(dMeshTopology%mesh, CharBuffer, setId, dMeshTopology%Node_Set(iSet)%Num_Nodes, ierr); CHKERRQ(iErr)
            Allocate(dMeshTopology%Node_Set(iSet)%Node_ID(dMeshTopology%Node_Set(iSet)%Num_Nodes))
            Call MeshGetStratum(dMeshTopology%mesh, CharBuffer, setId, dMeshTopology%Node_Set(iSet)%Node_ID, ierr); CHKERRQ(iErr)
            dMeshTopology%Node_Set(iSet)%Node_ID = dMeshTopology%Node_Set(iSet)%Node_ID - dMeshTopology%Num_Elems + 1
         End Do
      End If
      Deallocate(setIds)
   End Subroutine MeshTopologyReadEXO
    
   Subroutine MeshInitCoordinatesVect2D(dMeshTopology, dCoords)
      Type(MeshTopology_Info)                      :: dMeshTopology
      Type(Vect2D), Dimension(:), Pointer          :: dCoords
      
      PetscReal, Dimension(:,:), Pointer           :: array
      PetscInt                                     :: iErr

      Call MeshGetCoordinatesF90(dMeshTopology%mesh, array, iErr); CHKERRQ(iErr)
      dCoords%X = array(:,1)
      dCoords%Y = array(:,2)
      Call MeshRestoreCoordinatesF90(dMeshTopology%mesh, array, iErr); CHKERRQ(iErr)
   End Subroutine MeshInitCoordinatesVect2D 

   Subroutine MeshInitCoordinatesVect3D(dMeshTopology, dCoords)
      Type(MeshTopology_Info)                      :: dMeshTopology
      Type(Vect3D), Dimension(:), Pointer          :: dCoords

      PetscReal, Dimension(:,:), Pointer           :: array
      PetscInt                                     :: iErr

      Call MeshGetCoordinatesF90(dMeshTopology%mesh, array, iErr); CHKERRQ(iErr)
      dCoords%X = array(:,1)
      dCoords%Y = array(:,2)
      dCoords%Z = array(:,3)
      Call MeshRestoreCoordinatesF90(dMeshTopology%mesh, array, iErr); CHKERRQ(iErr)
   End Subroutine MeshInitCoordinatesVect3D 

   Subroutine MeshInitElementConnectivity2D_Scal(dMeshTopology, dElem)
      Type(MeshTopology_Info)                      :: dMeshTopology
      Type(Element2D_Scal), Dimension(:), Pointer  :: dElem
      
      PetscInt, Dimension(:,:), Pointer            :: arrayCon
      PetscInt                                     :: iBlk, iE, iEloc
      PetscInt                                     :: iErr
      
      Call MeshGetElementsF90(dMeshTopology%mesh, arrayCon, iErr); CHKERRQ(iErr)
      !!! What would happen if the number of dof per  element wasn't constant?
      Do iBlk = 1, dMeshTopology%Num_Elem_Blks
         Do iEloc = 1, dMeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = dMeshTopology%Elem_Blk(iBlk)%Elem_ID(iEloc)
            Allocate(dElem(iE)%ID_DoF(dMeshTopology%Elem_Blk(iBlk)%Num_DoF))
            dElem(iE)%ID_DoF = arrayCon(iE,:)
         End Do
      End Do
      Call MeshRestoreElementsF90(dMeshTopology%mesh, arrayCon, iErr); CHKERRQ(iErr)
    End Subroutine MeshInitElementConnectivity2D_Scal

   Subroutine MeshInitElementConnectivity2D(dMeshTopology, dElem)
      Type(MeshTopology_Info)                      :: dMeshTopology
      Type(Element2D), Dimension(:), Pointer       :: dElem
      
      PetscInt, Dimension(:,:), Pointer            :: arrayCon
      PetscInt                                     :: iBlk, iE, iEloc
      PetscInt                                     :: iErr
      
      Call MeshGetElementsF90(dMeshTopology%mesh, arrayCon, iErr); CHKERRQ(iErr)
      !!! What would happen if the number of dof per  element wasn't constant?
      Do iBlk = 1, dMeshTopology%Num_Elem_Blks
         Do iEloc = 1, dMeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = dMeshTopology%Elem_Blk(iBlk)%Elem_ID(iEloc)
            Allocate(dElem(iE)%ID_DoF(dMeshTopology%Elem_Blk(iBlk)%Num_DoF))
            dElem(iE)%ID_DoF = arrayCon(iE,:)
         End Do
      End Do
      Call MeshRestoreElementsF90(dMeshTopology%mesh, arrayCon, iErr); CHKERRQ(iErr)
    End Subroutine MeshInitElementConnectivity2D

   Subroutine MeshInitElementConnectivity2D_Elast(dMeshTopology, dElem)
      Type(MeshTopology_Info)                      :: dMeshTopology
      Type(Element2D_Elast), Dimension(:), Pointer :: dElem
      
      PetscInt, Dimension(:,:), Pointer            :: arrayCon
      PetscInt                                     :: iBlk, iE, iEloc
      PetscInt                                     :: iErr
      
      Call MeshGetElementsF90(dMeshTopology%mesh, arrayCon, iErr); CHKERRQ(iErr)
      !!! What would happen if the number of dof per  element wasn't constant?
      Do iBlk = 1, dMeshTopology%Num_Elem_Blks
         Do iEloc = 1, dMeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = dMeshTopology%Elem_Blk(iBlk)%Elem_ID(iEloc)
            Allocate(dElem(iE)%ID_DoF(dMeshTopology%Elem_Blk(iBlk)%Num_DoF))
            dElem(iE)%ID_DoF = arrayCon(iE,:)
         End Do
      End Do
      Call MeshRestoreElementsF90(dMeshTopology%mesh, arrayCon, iErr); CHKERRQ(iErr)
    End Subroutine MeshInitElementConnectivity2D_Elast

   Subroutine MeshInitElementConnectivity3D_Scal(dMeshTopology, dElem)
      Type(MeshTopology_Info)                      :: dMeshTopology
      Type(Element3D_Scal), Dimension(:), Pointer  :: dElem
      
      PetscInt, Dimension(:,:), Pointer            :: arrayCon
      PetscInt                                     :: iBlk, iE, iEloc
      PetscInt                                     :: iErr
      
      Call MeshGetElementsF90(dMeshTopology%mesh, arrayCon, iErr); CHKERRQ(iErr)
      !!! What would happen if the number of dof per  element wasn't constant?
      Do iBlk = 1, dMeshTopology%Num_Elem_Blks
         Do iEloc = 1, dMeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = dMeshTopology%Elem_Blk(iBlk)%Elem_ID(iEloc)
            Allocate(dElem(iE)%ID_DoF(dMeshTopology%Elem_Blk(iBlk)%Num_DoF))
            dElem(iE)%ID_DoF = arrayCon(iE,:)
         End Do
      End Do
      Call MeshRestoreElementsF90(dMeshTopology%mesh, arrayCon, iErr); CHKERRQ(iErr)
    End Subroutine MeshInitElementConnectivity3D_Scal

   Subroutine MeshInitElementConnectivity3D(dMeshTopology, dElem)
      Type(MeshTopology_Info)                      :: dMeshTopology
      Type(Element3D), Dimension(:), Pointer       :: dElem
      
      PetscInt, Dimension(:,:), Pointer            :: arrayCon
      PetscInt                                     :: iBlk, iE, iEloc
      PetscInt                                     :: iErr
      
      Call MeshGetElementsF90(dMeshTopology%mesh, arrayCon, iErr); CHKERRQ(iErr)
      !!! What would happen if the number of dof per  element wasn't constant?
      Do iBlk = 1, dMeshTopology%Num_Elem_Blks
         Do iEloc = 1, dMeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = dMeshTopology%Elem_Blk(iBlk)%Elem_ID(iEloc)
            Allocate(dElem(iE)%ID_DoF(dMeshTopology%Elem_Blk(iBlk)%Num_DoF))
            dElem(iE)%ID_DoF = arrayCon(iE,:)
         End Do
      End Do
      Call MeshRestoreElementsF90(dMeshTopology%mesh, arrayCon, iErr); CHKERRQ(iErr)
    End Subroutine MeshInitElementConnectivity3D

   Subroutine MeshInitElementConnectivity3D_Elast(dMeshTopology, dElem)
      Type(MeshTopology_Info)                      :: dMeshTopology
      Type(Element3D_Elast), Dimension(:), Pointer :: dElem
      
      PetscInt, Dimension(:,:), Pointer            :: arrayCon
      PetscInt                                     :: iBlk, iE, iEloc
      PetscInt                                     :: iErr
      
      Call MeshGetElementsF90(dMeshTopology%mesh, arrayCon, iErr); CHKERRQ(iErr)
      !!! What would happen if the number of dof per  element wasn't constant?
      Do iBlk = 1, dMeshTopology%Num_Elem_Blks
         Do iEloc = 1, dMeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = dMeshTopology%Elem_Blk(iBlk)%Elem_ID(iEloc)
            Allocate(dElem(iE)%ID_DoF(dMeshTopology%Elem_Blk(iBlk)%Num_DoF))
            dElem(iE)%ID_DoF = arrayCon(iE,:)
         End Do
      End Do
      Call MeshRestoreElementsF90(dMeshTopology%mesh, arrayCon, iErr); CHKERRQ(iErr)
    End Subroutine MeshInitElementConnectivity3D_Elast

End Module m_MEF_Sieve