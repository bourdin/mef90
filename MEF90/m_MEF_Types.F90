Module m_MEF_Types
   Use m_AlgebLin

   IMPLICIT NONE
   Private
   
   Public :: Element1D
   Public :: Element2D, Element2D_Scal, Element2D_Elast 
   Public :: Element3D, Element3D_Scal, Element3D_Elast 

   Public :: Elem_Blk_Info, Node_Set_Info, MeshTopology_Info, EXO_Info
   
   Public :: EXOView
   Public :: MeshTopologyDestroy, MeshTopologyView
      
#include "finclude/petsc.h"
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscao.h"
#include "finclude/petscmesh.h"
   include "exodusII.inc"

   ! Defines the basid bata structures for nodes, elements and geometry
   
   ! Basic element data structure:
   ! NB_DoF              Number of degree of freedoms in the element
   ! NB_GAUSS            Number of integration points
   ! ID_EL               Element Type
   ! ID_DoF              Index of teh DoF (Dim = NB_DoF)
   ! BF                  Value of the Basis Fuctions at the integration
   !                     Type depends on the element
   !                             BF(i,j) = i^th BF at the  j^th integration point
   ! DER_BF / rad_BF / GradS_BF Derivative / gradient / symmetrized gradient of the BF


!!! mef90-sieve:
!!! Do I still need to sequence my data structures?

!!! REMOVE
!!! Type Node{1,2,3}D 
!!! Layout_Info ???

!!! MODIFY      
!!! Type element:
!!!      Remove NB_DoF, NB_Gauss, ID_EL and add an element type to blocks
!!!      Add parent block information?
!!!      Can we possibly need to have elements owned by more than one block? I am assuming that we don't
!!! 
!!! Elem_Blk_Info
!!!      Remove Type, replace with Element_Info

!!! ADD
!!! Type Element_Info
!!!
!!! Various element names

!!! RENAME
!!! EXO_Geom_Info is now Geom_Info

   Type Element1D
!      Integer                                    :: NB_DoF
!      Integer                                    :: NB_Gauss     
!      Integer                                    :: ID_EL
      Integer, Dimension(:), pointer             :: ID_DoF
      Real(Kind = Kr), Dimension(:,:), pointer   :: BF
      Real(Kind = Kr), Dimension(:,:), pointer   :: Der_BF
      Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
!      Integer                                    :: Parent_Block
   End Type Element1D
 
   Type Element2D_Scal
!      Integer                                    :: NB_DoF       
!      Integer                                    :: NB_Gauss     
!      Integer                                    :: ID_EL
      Integer, Dimension(:), Pointer             :: ID_DoF
      Real(Kind = Kr), Dimension(:,:), pointer   :: BF
      Type(Vect2D), Dimension(:,:), pointer      :: Grad_BF
      Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
!      Integer                                    :: Parent_Block
   End Type Element2D_Scal
 
   Type Element2D
!      Integer                                    :: NB_DoF       
!      Integer                                    :: NB_Gauss     
!      Integer                                    :: ID_EL
      Integer, Dimension(:), pointer             :: ID_DoF
      Type (Vect2D), Dimension(:,:), pointer     :: BF
      Type (Mat2D), Dimension(:,:), pointer      :: Der_BF
      Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
!      Integer                                    :: Parent_Block
   End Type Element2D
 
   Type Element2D_Elast
!      Integer                                    :: NB_DoF       
!      Integer                                    :: NB_Gauss     
!      Integer                                    :: ID_EL
      Integer, Dimension(:), pointer             :: ID_DoF
      Type (Vect2D), Dimension(:,:), pointer     :: BF
      Type (MatS2D), Dimension(:,:), pointer     :: GradS_BF
      Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
!      Integer                                    :: Parent_Block
   End Type Element2D_Elast
 
   Type Element3D
!      Integer                                    :: NB_DoF       
!      Integer                                    :: NB_Gauss     
!      Integer                                    :: ID_EL
      Integer, Dimension(:), pointer             :: ID_DoF
      Type (Vect3D), Dimension(:,:), pointer     :: BF
      Type (Mat3D), Dimension(:,:), pointer      :: Der_BF
      Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
!      Integer                                    :: Parent_Block
   End Type Element3D
 
   Type Element3D_Scal
!      Integer                                    :: NB_DoF       
!      Integer                                    :: NB_Gauss     
!      Integer                                    :: ID_EL
      Integer, Dimension(:), pointer             :: ID_DoF
      Real(Kind = Kr), Dimension(:,:), pointer   :: BF
      Type (Vect3D), Dimension(:,:), pointer     :: Grad_BF
      Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
!      Integer                                    :: Parent_Block
   End Type Element3D_Scal
 
   Type Element3D_Elast
!      Integer                                    :: NB_DoF       
!      Integer                                    :: NB_Gauss     
!      Integer                                    :: ID_EL
      Integer, Dimension(:), pointer             :: ID_DoF
      Type (Vect3D), Dimension(:,:), pointer     :: BF
      Type (MatS3D), Dimension(:,:), pointer     :: GradS_BF
      Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
!      Integer                                    :: Parent_Block
   End Type Element3D_Elast
 
!   Type Node1D
!      Sequence
!      Real(Kind = Kr)                            :: Coord
!      Integer                                    :: ID
!      Integer                                    :: BC
!      Integer                                    :: Parent_Block
!   End Type Node1D
 
!   Type Node2D
!      Sequence
!      Type (Vect2D)                              :: Coord
!      Integer                                    :: ID
!      Integer                                    :: BC
!   End Type Node2D
 
!   Type Node3D
!      Sequence
!      Type (Vect3D)                              :: Coord
!      Integer                                    :: ID
!      Integer                                    :: BC
!   End Type Node3D
 
   Type Elem_Blk_Info
      Sequence
      Integer                                        :: ID
      Integer                                        :: Elem_Type
      Integer, Dimension(4)                          :: DoF_Location
      !!! Edge location is Cells, faces, edges, vertices in 3D and
      !!!                  Cell, unused, edges, vertices in2D
      Integer                                        :: Num_DoF !! = sum(DoF_Location)
!      Integer                                        :: NB_Gauss
      Integer                                        :: Num_Elems
      Integer, Dimension(:), Pointer                 :: Elem_ID
   End Type Elem_Blk_Info
 
   Type Node_Set_Info
      Sequence
      Integer                                        :: ID
      Integer                                        :: Num_Nodes
!      Integer                                        :: Num_Dist_Factors
      Integer, Dimension(:), Pointer                 :: Node_ID
!      Real(Kind = Kr), Dimension(:), Pointer         :: Dist_Factor
   End Type Node_Set_Info
 
   Type MeshTopology_Info
      Sequence
      ! Global datas
      Integer                                        :: num_dim
      Integer                                        :: num_verts
      Integer                                        :: num_elems
!      Integer                                        :: num_ghost_nodes
!      Integer                                        :: num_ghost_elems
      ! Element Blocks datas
      Integer                                        :: num_elem_blks
      Type(Elem_Blk_Info), Dimension(:), Pointer     :: elem_blk
      ! Node sets datas
      Integer                                        :: num_node_sets 
      Type(Node_Set_Info), Dimension(:), Pointer     :: node_set
      ! Side Sets DATAS
      Integer                                        :: num_side_sets
      Mesh                                           :: mesh
   End Type MeshTopology_Info
   
   Type EXO_Info
      MPI_Comm                                       :: comm      
      Integer                                        :: exoid
      Character(len=MXLNLN)                          :: filename  
      Character(len=MXLNLN)                          :: title     
      ! QA DATAS
      Integer                                        :: num_QA    
      Character(len=MXSTLN), Dimension(:,:), Pointer :: QA_rec    
   End Type EXO_Info
   
Contains
   Subroutine EXOView(dEXO, viewer)
      Type(EXO_Info)              :: dEXO
      PetscViewer                 :: viewer
      
      Integer                     :: iErr
      Character(len=512)          :: CharBuffer
   
      If (dEXO%comm == PETSC_COMM_WORLD) Then
         Write(CharBuffer, 100) 'PETSC_COMM_WORLD'
      ElseIf (dEXO%comm == PETSC_COMM_SELF) Then
         Write(CharBuffer, 100) 'PETSC_COMM_SELF'
      Else  
         Write(CharBuffer, 105) 'unknown', dEXO%comm
      End If
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      
      Write(CharBuffer, 101) dEXO%exoid
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 102) dEXO%filename
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 103) dEXO%title
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
!      Do i = 1, dEXO%Num_QA
!         Write(CharBuffer, 104) i, dEXO%QA_rec(i,:)
!         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
!      End Do
      
 100 Format('Communicator:       ', A, '\n'c)
 101 Format('exo ID:             ', I3, '\n'c)
 102 Format('filename:           ', A, '\n'c)
 103 Format('title:              ', A, '\n'c)
 104 Format('QR_rec ', I2.2, '         ', A, '\n'c)
 105 Format('Communicator:       ', A, I3, '\n'c)
   End Subroutine EXOView


   Subroutine MeshTopologyDestroy(dMeshTopology)
     Type (MeshTopology_Info)        :: dMeshTopology
     PetscInt                        :: iSet, iBlk

     If (dMeshTopology%Num_Node_Sets > 0) Then
        Do iSet = 1, dMeshTopology%Num_node_sets
           Deallocate(dMeshTopology%Node_Set(iSet)%Node_ID)
        End Do
     End If
     Deallocate (dMeshTopology%Node_Set)
     If (dMeshTopology%Num_Elem_blks > 0) Then
        Do iBlk = 1, dMeshTopology%Num_Elem_Blks
           Deallocate(dMeshTopology%Elem_blk(iBlk)%Elem_ID)
        End Do
     End If
     Deallocate(dMeshTopology%Elem_blk)
   End Subroutine MeshTopologyDestroy

   Subroutine MeshTopologyView(dMeshTopology, viewer)
      Type (MeshTopology_Info), Intent(IN)           :: dMeshTopology
      PetscViewer                                    :: viewer
      
      Integer                                        :: iErr
      Character(len=512)                             :: CharBuffer
      Integer                                        :: i, j

      Write(CharBuffer, 103) dMeshTopology%num_dim
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 104) dMeshTopology%num_Verts
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 105) dMeshTopology%num_elems
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 106) dMeshTopology%Num_Elem_blks
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 107) dMeshTopology%Num_Node_Sets
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 108) dMeshTopology%Num_Side_Sets
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
   
      Write(CharBuffer, 600) '\n'c
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 200)
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 201) dMeshTopology%num_elem_blks
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Do i = 1, dMeshTopology%Num_Elem_blks
         Write(CharBuffer, 203) dMeshTopology%Elem_Blk(i)%ID, dMeshTopology%Elem_blk(i)%Num_Elems
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Write(CharBuffer, 204) dMeshTopology%Elem_Blk(i)%ID, dMeshTopology%Elem_blk(i)%Elem_Type
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Write(CharBuffer, 205) dMeshTopology%Elem_Blk(i)%ID, dMeshTopology%Elem_blk(i)%DoF_Location
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Write(CharBuffer, 207) dMeshTopology%Elem_Blk(i)%ID
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do j = 1, dMeshTopology%Elem_blk(i)%Num_DoF
            Write(CharBuffer, 208) dMeshTopology%Elem_blk(i)%Elem_ID(j)
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, 600) '\n'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      End Do
      
      
      Write(CharBuffer, 600) '\n'c
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 300)
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 301) dMeshTopology%num_node_sets
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Do i = 1, dMeshTopology%num_node_sets
         Write(CharBuffer, 302) dMeshTopology%Node_Set(i)%ID, dMeshTopology%Node_Set(i)%Num_Nodes
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Write(CharBuffer, 303) dMeshTopology%Node_Set(i)%ID
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do j = 1, dMeshTopology%Node_Set(i)%Num_Nodes
            Write(CharBuffer, 500) dMeshTopology%Node_Set(i)%Node_ID(j)
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
      End Do
      
      
      Write(CharBuffer, 600) '\n'c
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 400)
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 401) dMeshTopology%num_side_sets
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      
103 Format('    Number of dimensions ============ ', I6, '\n'c)
104 Format('    Number of vertices ============== ', I6, '\n'c)
105 Format('    Number of elements ============== ', I6, '\n'c)
106 Format('    Number of elements blocks ======= ', I6, '\n'c)
107 Format('    Number of node sets ============= ', I6, '\n'c)
108 Format('    Number of side sets ============= ', I6, '\n'c)

200 Format('*** ELEMENT BLOCKS ***', '\n'c)
201 Format('    Number of blocks ================ ', I4, '\n'c)
203 Format('    Block ', I3, ' Number of elements ==== ', I4, '\n'c)
204 Format('    Block ', I3, ' Element type ========== ', I4, '\n'c)
205 Format('    Block ', I3, ' DoF location ========== ', 4(I4, ' '), '\n'c)
207 Format('    Block ', I3, ' IDs: ')
208 Format('   ', I4)

300 Format('*** NODE SETS ***', '\n'c)
301 Format('    Number of sets ================== ', I4, '\n'c)
302 Format('    Set ', I3, ' Number of nodes ========= ', I4, '\n'c)
303 Format('    Set ', I3, ' IDs: ')
!304 Format('    Set ', I3, ' Number of dist. factors = ', I4)

400 Format('*** SIDE SETS ***', '\n'c)
401 Format('    Number of side sets ============= ', I4, '\n'c)
    
500 Format(I4) 
600 Format(A)
   End Subroutine MeshTopologyView


End Module m_MEF_Types
