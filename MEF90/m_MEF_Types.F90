Module m_Mef_types
Use m_AlgebLin

  IMPLICIT NONE

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

  Type Element2DA
     Integer                                    :: NB_DoF       
     Integer                                    :: NB_Gauss     
     Integer                                    :: ID_EL
     Integer, Dimension(:), Pointer             :: ID_DoF
     Real(Kind = Kr), Dimension(:,:), pointer   :: BF
     Type(Vect2D), Dimension(:,:), pointer      :: Der_BF
     Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
  End Type Element2DA

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
     Character(len=MXLNLN)                          :: filename
     Character(len=MXLNLN)                          :: title
     Integer                                        :: exoid
     Integer                                        :: num_dim
     Integer                                        :: num_nodes 
     Integer                                        :: num_elems
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
End Module m_Mef_types
