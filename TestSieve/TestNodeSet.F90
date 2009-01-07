Program TestSieve

   Use m_MEF90
   Implicit NONE
    
#include "finclude/petsc.h"
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscviewer.h"
#include "finclude/petscviewer.h90"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscao.h"
#include "finclude/petscmesh.h"
#include "finclude/petscmesh.h90"

   Type (MeshTopology_Info)                     :: MeshTopology, GlobalMeshTopology
   Type (EXO_Info)                              :: EXO, MyEXO
   Type (Element2D_Scal), Dimension(:), Pointer :: Elem2DA
   Type (Vect3D), Dimension(:), Pointer         :: Coords
!   PetscReal, Dimension(:,:), Pointer           :: Vertices
   
   PetscTruth                                   :: HasPrefix
   PetscTruth                                   :: verbose
   PetscErrorCode                               :: iErr
   Integer                                      :: iBlk
   Character(len=256)                           :: CharBuffer, filename
   Character(len=256)                           :: prefix
   PetscViewer                                  :: viewer, myviewer
     
   Call MEF90_Initialize()
   Call PetscOptionsHasName(PETSC_NULL_CHARACTER, '-verbose', verbose, iErr)    
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p', prefix, HasPrefix, iErr)    
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No input file prefix given\n", iErr)
      Call MEF90_Finalize()
      STOP
   End If

   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'


   !!! Clean that later
!   If(MEF90_MyRank == 0) Then
      Call MeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, GlobalMeshTopology%mesh, ierr)
      Call MeshExodusGetInfo(GlobalMeshTopology%mesh, GlobalMeshTopology%Num_Dim, GlobalMeshTopology%Num_Vert, GlobalMeshTopology%Num_Elems, GlobalMeshTopology%Num_Elem_Blks, GlobalMeshTopology%Num_Node_Sets, iErr)
      Call MeshDestroy(GlobalMeshTopology%mesh, ierr)
!   End If
   Call MPI_BCast(GlobalMeshTopology%Num_Dim, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, iErr)
   Call MPI_BCast(GlobalMeshTopology%Num_Vert, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, iErr)
   Call MPI_BCast(GlobalMeshTopology%Num_Elems, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, iErr)
   Call MPI_BCast(GlobalMeshTopology%Num_Elem_Blks, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, iErr)
   Call MPI_BCast(GlobalMeshTopology%Num_Node_Sets, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, iErr)
   !!!

   Call MeshTopologyReadEXO(MeshTopology, Coords, Elem2DA, EXO)
   
   MeshTopology%Elem_Blk%Elem_Type    = MEF90_P1_Lagrange
   Do iBlk = 1, MeshTopology%Num_Elem_Blks
      Call Init_Elem_Blk_Info(MeshTopology%Elem_Blk(iBlk), MeshTopology%num_dim)
   End Do
   
   If (verbose) Then
      Write(filename, 102) trim(prefix), MEF90_MyRank
      Call PetscViewerASCIIOpen(PETSC_COMM_SELF, filename, myviewer, iErr); CHKERRQ(iErr);   
      Write(CharBuffer, 103) MEF90_MyRank, trim(filename)
      Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, CharBuffer, iErr); CHKERRQ(iErr)
      Call PetscSynchronizedFlush(PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
   End If

102 Format(A, '-', I4.4, '.log')
103 Format('Output from processor ', I4.4, ' redirected to ', A, '\n'c)

   MyEXO%comm = PETSC_COMM_SELF
   MyEXO%exoid = EXO%exoid
   Write(MyEXO%filename, 200) trim(prefix), MEF90_MyRank
 200 Format(A, '-', I4.4, '.gen')
 
   MyEXO%title = trim(EXO%title)
   MyEXO%Num_QA = EXO%Num_QA
   

   MeshTopology%elem_blk(:)%Elem_type = MEF90_P1_Lagrange
   Call Write_MeshTopologyGlobal(MeshTopology, MyEXO, PETSC_COMM_WORLD)
!   Call Write_MeshTopology(MeshTopology, MyEXO)

   If (verbose) Then
      Write(CharBuffer, 300) 'EXO_Info\n'c
      Call PetscViewerASCIIPrintf(myviewer, CharBuffer, iErr); CHKERRQ(iErr)
      Call EXOView(EXO, myviewer)

      Write(CharBuffer, 300) '\n\nMeshTopology\n'c
      Call PetscViewerASCIIPrintf(myviewer, CharBuffer, iErr); CHKERRQ(iErr)
      Call MeshTopologyView(MeshTopology, myviewer)
      
      Write(CharBuffer, 300) '\n\nNode Sets\n'c
      Call PetscViewerASCIIPrintf(myviewer, CharBuffer, iErr); CHKERRQ(iErr)
      Call ShowNodesets(MeshTopology%mesh, myviewer)
   End If
300 Format(A)

   Call MeshTopologyDestroy(MeshTopology)

   If (verbose) Then
      Call PetscViewerFlush(myviewer); CHKERRQ(iErr)
      Call PetscViewerDestroy(myviewer, iErr); CHKERRQ(iErr)
   End If
   
   
   Call MEF90_Finalize()

 Contains
   Subroutine ShowNodeSets(dmesh, viewer)
      Mesh                                         :: dmesh
      PetscViewer                                  :: viewer
      
      Type(MeshTopology_Info)                      :: lMeshTopology
      Character(len=512)                           :: CharBuffer, IOBuffer
      Integer                                      :: iElem, numIds, blkId, setId, iSet, i
      Integer                                      :: iErr
      PetscInt, Dimension(:), Pointer              :: blkIds
      PetscInt, Dimension(:), Pointer              :: setIds

      lMeshTopology%mesh = dmesh
      ! Read Global Geometric Parameters
      call MeshExodusGetInfo(lMeshTopology%mesh, lMeshTopology%Num_Dim, lMeshTopology%Num_Vert, lMeshTopology%Num_Elems, lMeshTopology%Num_Elem_Blks, lMeshTopology%Num_Node_Sets, iErr)
      
      ! Read Node set information
      CharBuffer = 'VertexSets'
      call MeshGetLabelSize(lMeshTopology%mesh, CharBuffer, numIds, ierr); CHKERRQ(ierr)
      If (numIds .ne. lMeshTopology%Num_node_sets) Then
         SETERRQ(PETSC_ERR_ARG_SIZ, 'Invalid number of node ids', ierr)
      End If
      
      Write(IOBuffer, 100) numIds
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      
      Allocate(lMeshTopology%Node_Set(lMeshTopology%Num_Node_Sets))
      Allocate(setIds(numIds))
      call MeshGetLabelIds(lMeshTopology%mesh, CharBuffer, setIds, ierr); CHKERRQ(ierr)
      Write(IOBuffer, 101) 
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Do iSet = 1, lMeshTopology%Num_node_sets
         Write(IOBuffer, 102) setIds(iSet)
         Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      End Do
      Write(IOBuffer, 200) '\n'c
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
         
      Do iSet = 1, lMeshTopology%Num_node_sets
         setId = setIds(iSet)
         Write(IOBuffer, 103) setID
         Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
         
         lMeshTopology%Node_Set(iSet)%ID = setId
         call MeshGetStratumSize(lMeshTopology%mesh, CharBuffer, setId, lMeshTopology%Node_Set(iSet)%Num_Nodes, ierr)
         Write(IOBuffer, 104) lMeshTopology%Node_Set(iSet)%Num_Nodes
         Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
         
         Allocate(lMeshTopology%Node_Set(iSet)%Node_ID(lMeshTopology%Node_Set(iSet)%Num_Nodes))
         Write(IOBuffer, 200) '   Vertices:'
         Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
         Do i = 1, lMeshTopology%Node_Set(iSet)%Num_Nodes
            call MeshGetStratum(lMeshTopology%mesh, CharBuffer, setId, lMeshTopology%Node_Set(iSet)%Node_ID, ierr)
            Write(IOBuffer, 105) lMeshTopology%Node_Set(iSet)%Node_ID(i) - lMeshTopology%Num_Elems + 1
            Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(IOBuffer, 200) '\n'c
         Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
         DeAllocate(lMeshTopology%Node_Set(iSet)%Node_ID)
      End Do
   Deallocate(setIds)
   DeAllocate(lMeshTopology%Node_Set)
   
100 Format('Number of Vertex Sets', I4, '\n'c)
101 Format('Ids of VertexSets:')
102 Format('   ', I3)
103 Format('VertexSet ', I3, '\n'c)
104 Format('   Number of Vertices in VertexSet', I5, '\n'c) 
105 Format('  ', I5)
200 Format(A)
   End Subroutine ShowNodeSets
End Program TestSieve
