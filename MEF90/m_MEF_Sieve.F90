Module m_MEF_Sieve

   Use m_AlgebLin
   Use m_MEF_Types
   
   IMPLICIT NONE
   Private
   
#include "finclude/petsc.h"
#include "finclude/petscviewer.h"
#include "finclude/petscviewer.h90"
#include "finclude/petscmesh.h"
#include "finclude/petscmesh.h90"

   include "exodusII.inc"

   Public :: MeshTopologyReadEXO
   
   Interface MeshTopologyReadEXO
      Module Procedure MeshTopologyReadEXO_2DScal
   End Interface MeshTopologyReadEXO
   
!   Interface MeshCreateCoordinates
!      Module Procedure MeshCreateCoordinatesPtr MeshCreateCoordinatesVect2D MeshCreateCoordinatesVect3D
!   End Interface MeshCreateCoordinates
   
!   Interface Interface MeshInitElemConnectivity
!      Module Procedure Interface MeshInitElemConnectivity2D_Scal, MeshInitElemConnectivity2D, MeshInitElemConnectivity2D_Elast, Interface MeshInitElemConnectivity3D_Scal, MeshInitElemConnectivity3D, MeshInitElemConnectivity3D_Elast
!   End Interface MeshInitElemConnectivity
      
Contains

!   Subroutine MeshTopologyReadGeometryEXO(dMeshTopology, dEXO)
!   End Subroutine MeshTopologyReadGeometryEXO
   
!   Subroutine MeshGeometryGetConnectivity_2DScal(dMeshTopology, Elem2DA)
!   End Subroutine MeshGeometryGetConnectivity_2DScal
   
   
   Subroutine MeshTopologyReadEXO_2DScal(dMeshTopology, Coords, Elem2DA, dEXO)
      !!! Remove the element and coordinate stuff and move in separate functions
      Type (MeshTopology_Info)                     :: dMeshTopology
      Type (Vect3D), Dimension(:), Pointer         :: Coords
      Type (Element2D_Scal), Dimension(:), Pointer :: Elem2DA
      Type (EXO_Info)                              :: dEXO
      Integer                                      :: iErr, iBlk, iSet
      Character(len=256)                           :: CharBuffer
      
      Mesh                                         :: mesh
      PetscReal, Dimension(:,:), Pointer           :: array
      PetscInt, Dimension(:,:), Pointer            :: arrayCon
      Integer                                      :: embedDim
      Integer                                      :: iE, iElem, numIds, blkId, setId
      PetscInt, Dimension(:), Pointer              :: blkIds
      PetscInt, Dimension(:), Pointer              :: setIds
      
      ! Open File
      call MeshCreateExodus(PETSC_COMM_WORLD, dEXO%filename, mesh, ierr)
      !!! reads exo file, stores all information in a Mesh

      call MeshDistribute(mesh, PETSC_NULL_CHARACTER, dMeshTopology%mesh, ierr)
      !!! Partitions using a partitioner (currently PETSC_NULL_CHARACTER) 
      
      call MeshDestroy(mesh, ierr)

      ! Read Global Geometric Parameters
      call MeshExodusGetInfo(dMeshTopology%mesh, dMeshTopology%Num_Dim, dMeshTopology%Num_Verts, dMeshTopology%Num_Elems, dMeshTopology%Num_Elem_Blks, dMeshTopology%Num_Node_Sets, iErr)
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
      call MeshGetLabelIds(dMeshTopology%mesh, CharBuffer, blkIds, ierr); CHKERRQ(ierr)
      If (dMeshTopology%Num_Elem_blks > 0) Then
         Do iBlk = 1, dMeshTopology%Num_Elem_Blks
            blkId = blkIds(iBlk)
            dMeshTopology%Elem_blk(iBlk)%ID = blkId
            call MeshGetStratumSize(dMeshTopology%mesh, CharBuffer, blkId, dMeshTopology%elem_blk(iBlk)%Num_Elems, ierr)
            !!! Get the size of the layer (stratum) 'CellBlock' of Mesh
            Allocate(dMeshTopology%Elem_blk(iBlk)%Elem_ID(dMeshTopology%elem_blk(iBlk)%Num_Elems))
            call MeshGetStratum(dMeshTopology%mesh, CharBuffer, blkId, dMeshTopology%Elem_blk(iBlk)%Elem_ID, ierr)
            !!! Get the layer (stratum) 'CellBlock' of Mesh in C numbering
            dMeshTopology%Elem_blk(iBlk)%Elem_ID = dMeshTopology%Elem_blk(iBlk)%Elem_ID + 1
            !!! Converts to Fortran style indexing
         End Do
      End If
      Deallocate(blkIds)
      
      ! Read Node set information
      CharBuffer = 'VertexSets'
      call MeshGetLabelSize(dMeshTopology%mesh, CharBuffer, numIds, ierr); CHKERRQ(ierr)
      If (numIds .ne. dMeshTopology%Num_node_sets) Then
         SETERRQ(PETSC_ERR_ARG_SIZ, 'Invalid number of node ids', ierr)
      End If
      Allocate(dMeshTopology%Node_Set(dMeshTopology%Num_Node_Sets))
      Allocate(setIds(numIds))
      call MeshGetLabelIds(dMeshTopology%mesh, CharBuffer, setIds, ierr); CHKERRQ(ierr)
      If (dMeshTopology%Num_Node_Sets > 0) Then
         Do iSet = 1, dMeshTopology%Num_node_sets
            setId = setIds(iSet)
            dMeshTopology%Node_Set(iSet)%ID = setId
            call MeshGetStratumSize(dMeshTopology%mesh, CharBuffer, setId, dMeshTopology%Node_Set(iSet)%Num_Nodes, ierr)
            Allocate(dMeshTopology%Node_Set(iSet)%Node_ID(dMeshTopology%Node_Set(iSet)%Num_Nodes))
            call MeshGetStratum(dMeshTopology%mesh, CharBuffer, setId, dMeshTopology%Node_Set(iSet)%Node_ID, ierr)
            dMeshTopology%Node_Set(iSet)%Node_ID = dMeshTopology%Node_Set(iSet)%Node_ID - dMeshTopology%Num_Elems + 1
         End Do
      End If
      Deallocate(setIds)

      ! Read the vertices coordinates
      Allocate(Coords(dMeshTopology%Num_Verts))
      call MeshGetCoordinatesF90(dMeshTopology%mesh, array, iErr)
      embedDim = size(array,2)
      Coords%X = array(:,1)
      If (embedDim > 1) Then
         Coords%Y = array(:,2)
      Else
         Coords%Y = 0.0
      EndIf
      If (embedDim > 2) Then
         Coords%Z = array(:,3)
      Else
         Coords%z = 0.0
      EndIf
      call MeshRestoreCoordinatesF90(dMeshTopology%mesh, array, iErr)
 
      ! Read the connectivity table
      Allocate(Elem2DA(dMeshTopology%Num_Elems))
      call MeshGetElementsF90(dMeshTopology%mesh, arrayCon, iErr)
      Do iBlk = 1, dMeshTopology%Num_Elem_Blks
         Do iE = 1, dMeshTopology%Elem_Blk(iBlk)%Num_Elems
            iElem = dMeshTopology%Elem_Blk(iBlk)%Elem_ID(iE)
            Allocate(Elem2DA(iElem)%ID_DoF(3))
            Elem2DA(iElem)%ID_DoF = arrayCon(iElem,:)
         End Do
      End Do
      call MeshRestoreElementsF90(dMeshTopology%mesh, arrayCon, iErr)

      dEXO%exoid = 0
    End Subroutine MeshTopologyReadEXO_2DScal

End Module m_MEF_Sieve