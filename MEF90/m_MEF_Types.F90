Module m_MEF_Types
   Use m_AlgebLin

   IMPLICIT NONE
   Private
   
   Public :: Element1D
   Public :: Element2D, Element2D_Scal, Element2D_Elast 
   Public :: Element3D, Element3D_Scal, Element3D_Elast 

   Public :: Node1D, Node2D, Node3D
   Public :: Elem_Blk_Info, Node_Set_Info, EXO_Geom_Info
   Public :: Layout_Info
   
#include "include/finclude/petsc.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petscao.h"
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

   Type Element1D
      Integer                                    :: NB_DoF
      Integer                                    :: NB_Gauss     
      Integer                                    :: ID_EL
      Integer, Dimension(:), pointer             :: ID_DoF
      Real(Kind = Kr), Dimension(:,:), pointer   :: BF
      Real(Kind = Kr), Dimension(:,:), pointer   :: Der_BF
      Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
   End Type Element1D
 
   Type Element2D_Scal
      Integer                                    :: NB_DoF       
      Integer                                    :: NB_Gauss     
      Integer                                    :: ID_EL
      Integer, Dimension(:), Pointer             :: ID_DoF
      Real(Kind = Kr), Dimension(:,:), pointer   :: BF
      Type(Vect2D), Dimension(:,:), pointer      :: Grad_BF
      Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
   End Type Element2D_Scal
 
   Type Element2D
      Integer                                    :: NB_DoF       
      Integer                                    :: NB_Gauss     
      Integer                                    :: ID_EL
      Integer, Dimension(:), pointer             :: ID_DoF
      Type (Vect2D), Dimension(:,:), pointer     :: BF
      Type (Mat2D), Dimension(:,:), pointer      :: Der_BF
      Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
   End Type Element2D
 
   Type Element2D_Elast
      Integer                                    :: NB_DoF       
      Integer                                    :: NB_Gauss     
      Integer                                    :: ID_EL
      Integer, Dimension(:), pointer             :: ID_DoF
      Type (Vect2D), Dimension(:,:), pointer     :: BF
      Type (MatS2D), Dimension(:,:), pointer     :: GradS_BF
      Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
   End Type Element2D_Elast
 
   Type Element3D
      Integer                                    :: NB_DoF       
      Integer                                    :: NB_Gauss     
      Integer                                    :: ID_EL
      Integer, Dimension(:), pointer             :: ID_DoF
      Type (Vect3D), Dimension(:,:), pointer     :: BF
      Type (Mat3D), Dimension(:,:), pointer      :: Der_BF
      Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
   End Type Element3D
 
   Type Element3D_Scal
      Integer                                    :: NB_DoF       
      Integer                                    :: NB_Gauss     
      Integer                                    :: ID_EL
      Integer, Dimension(:), pointer             :: ID_DoF
      Real(Kind = Kr), Dimension(:,:), pointer   :: BF
      Type (Vect3D), Dimension(:,:), pointer     :: Grad_BF
      Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
   End Type Element3D_Scal
 
   Type Element3D_Elast
      Integer                                    :: NB_DoF       
      Integer                                    :: NB_Gauss     
      Integer                                    :: ID_EL
      Integer, Dimension(:), pointer             :: ID_DoF
      Type (Vect3D), Dimension(:,:), pointer     :: BF
      Type (MatS3D), Dimension(:,:), pointer     :: GradS_BF
      Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
   End Type Element3D_Elast
 
   Type Node1D
      Sequence
      Real(Kind = Kr)                            :: Coord
      Integer                                    :: ID
      Integer                                    :: BC
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
      Character(len=MXSTLN)                          :: Type
      Integer                                        :: Num_Elems
      Integer                                        :: Num_Nodes_Per_Elem
      Integer                                        :: Num_Attr
      Integer, Dimension(:), Pointer                 :: Elem_ID
   End Type Elem_Blk_Info
 
   Type Node_Set_Info
      Sequence
      Integer                                        :: ID
      Integer                                        :: Num_Nodes
      Integer                                        :: Num_Dist_Factors
      Integer, Dimension(:), Pointer                 :: Node_ID
      Real(Kind = Kr), Dimension(:), Pointer         :: Dist_Factor
   End Type Node_Set_Info
 
   Type EXO_Geom_Info
      Sequence
      ! Global datas
      MPI_Comm                                       :: comm
      Character(len=MXLNLN)                          :: filename
      Character(len=MXLNLN)                          :: title
      Integer                                        :: exoid
      Integer                                        :: num_dim
      Integer                                        :: num_nodes 
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
      ! QA DATAS
      Integer                                        :: num_QA
      Character(len=MXSTLN), Dimension(:,:), Pointer :: QA_rec
      Integer                                        :: Numbering
   End Type EXO_Geom_Info
   
   Type Layout_Info
      IS                                         :: IS_N, IS_E
      ISLocalToGlobalMapping                     :: Mapping_N, Mapping_E
      VecScatter                                 :: ToIOSeq_N, ToIOSeq_E
      !!! Scatters everything onto the IO node ordered component after component
      !!! Collective on Geom%Comm
      VecScatter                                 :: ToIODist_N
      !!! Scatters locally ordered component by component
      !!! Collective on PETSC_COMM_SELF (because of ghost points)
      
      Integer                                    :: num_local_dof, num_ghost_dof
      Integer, Dimension(:), Pointer             :: ghost_dof
      Integer                                    :: num_local_elems, num_ghost_elems
      Integer, Dimension(:), Pointer             :: ghost_elem
   End Type Layout_Info

   
End Module m_MEF_Types
