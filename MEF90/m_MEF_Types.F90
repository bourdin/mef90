Module m_MEF_Types
   Use m_AlgebLin

   IMPLICIT NONE
   Private
   
   Public :: Element1D
   Public :: Element2D, Element2D_Scal, Element2D_Elast 
   Public :: Element3D, Element3D_Scal, Element3D_Elast 

!   Public :: Node1D, Node2D, Node3D
   Public :: Elem_Blk_Info, Node_Set_Info, MeshTopology_Info, EXO_Info
      
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
      Integer                                    :: Parent_Block
   End Type Element1D
 
   Type Element2D_Scal
!      Integer                                    :: NB_DoF       
!      Integer                                    :: NB_Gauss     
!      Integer                                    :: ID_EL
      Integer, Dimension(:), Pointer             :: ID_DoF
      Real(Kind = Kr), Dimension(:,:), pointer   :: BF
      Type(Vect2D), Dimension(:,:), pointer      :: Grad_BF
      Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
      Integer                                    :: Parent_Block
   End Type Element2D_Scal
 
   Type Element2D
!      Integer                                    :: NB_DoF       
!      Integer                                    :: NB_Gauss     
!      Integer                                    :: ID_EL
      Integer, Dimension(:), pointer             :: ID_DoF
      Type (Vect2D), Dimension(:,:), pointer     :: BF
      Type (Mat2D), Dimension(:,:), pointer      :: Der_BF
      Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
      Integer                                    :: Parent_Block
   End Type Element2D
 
   Type Element2D_Elast
!      Integer                                    :: NB_DoF       
!      Integer                                    :: NB_Gauss     
!      Integer                                    :: ID_EL
      Integer, Dimension(:), pointer             :: ID_DoF
      Type (Vect2D), Dimension(:,:), pointer     :: BF
      Type (MatS2D), Dimension(:,:), pointer     :: GradS_BF
      Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
      Integer                                    :: Parent_Block
   End Type Element2D_Elast
 
   Type Element3D
!      Integer                                    :: NB_DoF       
!      Integer                                    :: NB_Gauss     
!      Integer                                    :: ID_EL
      Integer, Dimension(:), pointer             :: ID_DoF
      Type (Vect3D), Dimension(:,:), pointer     :: BF
      Type (Mat3D), Dimension(:,:), pointer      :: Der_BF
      Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
      Integer                                    :: Parent_Block
   End Type Element3D
 
   Type Element3D_Scal
!      Integer                                    :: NB_DoF       
!      Integer                                    :: NB_Gauss     
!      Integer                                    :: ID_EL
      Integer, Dimension(:), pointer             :: ID_DoF
      Real(Kind = Kr), Dimension(:,:), pointer   :: BF
      Type (Vect3D), Dimension(:,:), pointer     :: Grad_BF
      Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
      Integer                                    :: Parent_Block
   End Type Element3D_Scal
 
   Type Element3D_Elast
!      Integer                                    :: NB_DoF       
!      Integer                                    :: NB_Gauss     
!      Integer                                    :: ID_EL
      Integer, Dimension(:), pointer             :: ID_DoF
      Type (Vect3D), Dimension(:,:), pointer     :: BF
      Type (MatS3D), Dimension(:,:), pointer     :: GradS_BF
      Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
      Integer                                    :: Parent_Block
   End Type Element3D_Elast
 
   Type Node1D
      Sequence
      Real(Kind = Kr)                            :: Coord
      Integer                                    :: ID
      Integer                                    :: BC
      Integer                                    :: Parent_Block
   End Type Node1D
 
   Type Node2D
      Sequence
      Type (Vect2D)                              :: Coord
      Integer                                    :: ID
      Integer                                    :: BC
   End Type Node2D
 
   Type Node3D
      Sequence
      Type (Vect3D)                              :: Coord
      Integer                                    :: ID
      Integer                                    :: BC
   End Type Node3D
 
   Type Elem_Blk_Info
      Sequence
      Integer                                        :: ID
      Integer                                        :: Elem_Type
      Integer, Dimension(4)                          :: DoF_Location
      !!! Edge location is Cells, faces, edges, vertices in 3D and
      !!!                  Cell, unused, edges, vertices in2D
      Integer                                        :: Nb_DoF !! = sum(DoF_Location)
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
      Integer                                        :: num_vert 
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
      MPI_Comm                                       :: comm      ! Move somewhere else
      Integer                                        :: exoid
      Character(len=MXLNLN)                          :: filename  ! Move to EXO_DATA
      Character(len=MXLNLN)                          :: title     ! Move to EXO_Info
      ! QA DATAS
      Integer                                        :: num_QA    ! move to EXO_Info
      Character(len=MXSTLN), Dimension(:,:), Pointer :: QA_rec    ! move to EXO_Info   
   End Type EXO_Info
   
!!!   Type Layout_Info
!!!      IS                                         :: IS_N, IS_E
!!!      ISLocalToGlobalMapping                     :: Mapping_N, Mapping_E
!!!      VecScatter                                 :: ToIOSeq_N, ToIOSeq_E
!!!      !!! Scatters everything onto the IO node ordered component after component
!!!      !!! Collective on Geom%Comm
!!!      VecScatter                                 :: ToIODist_N
!!!      !!! Scatters locally ordered component by component
!!!      !!! Collective on PETSC_COMM_SELF (because of ghost points)
!!!      
!!!      Integer                                    :: num_local_dof, num_ghost_dof
!!!      Integer, Dimension(:), Pointer             :: ghost_dof
!!!      Integer                                    :: num_local_elems, num_ghost_elems
!!!      Integer, Dimension(:), Pointer             :: ghost_elem
!!!   End Type Layout_Info      
End Module m_MEF_Types
