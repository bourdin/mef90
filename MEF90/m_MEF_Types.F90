Module m_Mef_types
! Blaise Bourdin, 1996-1998
! Merci de faire parvenir toutes remarques et bugs 
! eventuels aux adresses suivantes :
! bourdin@lpmtm.univ-paris13.fr
! bourdin@mat.dtu.dk
!
! En cas de modifications, les noms de  modules et de 
! fichiers _DOIVENT_ etre renommes
!          ^^^^^^^^^

Use m_AlgebLin

  IMPLICIT NONE

  include "exodusII.inc"

  ! Definition de types elements :
  ! NB_DoF              Nombre de degres de libertes de l'element
  ! NB_GAUSS    Nombre de points de Gauss
  ! ID_EL               Type d'element
  ! ID_DoF              Numeros globaux des DL (Dim = NB_DoF)
  ! BF          Valeur des fonctions de base locales (BF)
  !                     aux points de Gauss :
  !                             BF(i,j) = i^eme BF au j^eme point de Gauss
  !                             (Type Vect2D, Vect3D ou Double)
  ! DER_BF      Valeur des derivees des BF aux point de Gauss :
  !                             DER_BF(i,j) = derivee de la i^eme BF au
  !                              j^eme point de Gauss
  !                             (Type Double, Vect2D, Mat2D ou Mat3D)

  !     DEFINITION DES TYPES VECT* ET MAT* dans m_algeblin


  ! Definition d'un type matrice profil (Mat_Prof):
  ! XXXX%profil Pointeur de compactage
  ! XXXX%Val    Termes stockes de la matrice
  ! Redefinition des operations arithmetiques sur Mat_Prof

  ! Attention : Risque de pbs avec les affectations de profils
  ! On risque de creer des infos redondantes

  ! Blaise Bourdin   04-97
  ! Version 3      09-97
  ! Version 0.7.0  09-99


  Type Element1D
     Integer(Kind = Ki)                         :: NB_DoF
     Integer(Kind = Ki)                         :: NB_Gauss     
     Integer(Kind = Ki)                         :: ID_EL
     Integer(Kind = Ki), Dimension(:), pointer  :: ID_DoF
     Real(Kind = Kr), Dimension(:,:), pointer   :: BF
     Real(Kind = Kr), Dimension(:,:), pointer   :: Der_BF
     Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
  End Type Element1D

  Type Element2DA
     Integer(Kind = Ki)                         :: NB_DoF       
     Integer(Kind = Ki)                         :: NB_Gauss     
     Integer(Kind = Ki)                         :: ID_EL
     Integer(Kind = Ki), Dimension(:), Pointer  :: ID_DoF
     Real(Kind = Kr), Dimension(:,:), pointer   :: BF
     Type(Vect2D), Dimension(:,:), pointer      :: Der_BF
     Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
  End Type Element2DA

  Type Element2D_Scal
     Integer(Kind = Ki)                         :: NB_DoF       
     Integer(Kind = Ki)                         :: NB_Gauss     
     Integer(Kind = Ki)                         :: ID_EL
     Integer(Kind = Ki), Dimension(:), Pointer  :: ID_DoF
     Real(Kind = Kr), Dimension(:,:), pointer   :: BF
     Type(Vect2D), Dimension(:,:), pointer      :: Grad_BF
     Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
  End Type Element2D_Scal

  Type Element2D
     Integer(Kind = Ki)                         :: NB_DoF       
     Integer(Kind = Ki)                         :: NB_Gauss     
     Integer(Kind = Ki)                         :: ID_EL
     Integer(Kind = Ki), Dimension(:), pointer  :: ID_DoF
     Type (Vect2D), Dimension(:,:), pointer     :: BF
     Type (Mat2D), Dimension(:,:), pointer      :: Der_BF
     Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
  End Type Element2D

  Type Element2D_Elast
     Integer(Kind = Ki)                         :: NB_DoF       
     Integer(Kind = Ki)                         :: NB_Gauss     
     Integer(Kind = Ki)                         :: ID_EL
     Integer(Kind = Ki), Dimension(:), pointer  :: ID_DoF
     Type (Vect2D), Dimension(:,:), pointer     :: BF
     Type (MatS2D), Dimension(:,:), pointer     :: GradS_BF
     Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
  End Type Element2D_Elast

  Type Element3D
     Integer(Kind = Ki)                         :: NB_DoF       
     Integer(Kind = Ki)                         :: NB_Gauss     
     Integer(Kind = Ki)                         :: ID_EL
     Integer(Kind = Ki), Dimension(:), pointer  :: ID_DoF
     Type (Vect3D), Dimension(:,:), pointer     :: BF
     Type (Mat3D), Dimension(:,:), pointer      :: Der_BF
     Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
  End Type Element3D

  Type Element3D_Scal
     Integer(Kind = Ki)                         :: NB_DoF       
     Integer(Kind = Ki)                         :: NB_Gauss     
     Integer(Kind = Ki)                         :: ID_EL
     Integer(Kind = Ki), Dimension(:), pointer  :: ID_DoF
     Real(Kind = Kr), Dimension(:,:), pointer   :: BF
     Type (Vect3D), Dimension(:,:), pointer     :: Grad_BF
     Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
  End Type Element3D_Scal

  Type Element3D_Elast
     Integer                                    :: NB_DoF       
     Integer                                    :: NB_Gauss     
     Integer                                    :: ID_EL
!!$     Integer(Kind = Ki)                         :: NB_DoF       
!!$     Integer(Kind = Ki)                         :: NB_Gauss     
!!$     Integer(Kind = Ki)                         :: ID_EL
     Integer(Kind = Ki), Dimension(:), pointer  :: ID_DoF
     Type (Vect3D), Dimension(:,:), pointer     :: BF
     Type (MatS3D), Dimension(:,:), pointer     :: GradS_BF
     Real(Kind = Kr), Dimension(:), Pointer     :: Gauss_C
  End Type Element3D_Elast

  Type Node1D
     Sequence
     Real(Kind = Kr)                            :: Coord
     Integer(Kind = Ki)                         :: ID
     Integer(Kind = Ki)                         :: BC
  End Type Node1D

  Type Node2D
     Sequence
     Type (Vect2D)                              :: Coord
     Integer(Kind = Ki)                         :: ID
     Integer(Kind = Ki)                         :: BC
  End Type Node2D

  Type Node3D
     Sequence
     Type (Vect3D)                              :: Coord
     Integer(Kind = Ki)                         :: ID
     Integer(Kind = Ki)                         :: BC
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

  Type Mat_Prof
     Integer(Kind = Ki), Dimension(:), Pointer      :: Profil
     Real(Kind = Kr),Dimension(:), Pointer          :: Val
  End Type Mat_Prof

! Data structures for a nice morse coding  (Compressed Row Coding)
  
  Type Mat_Morse_Entry
     Real(Kind = Kr)                                :: Val
     Integer(Kind = Ki)                             :: Pos
  End Type Mat_Morse_Entry

  Type Mat_Morse_Row
     Integer(Kind=Ki)                               :: Nb_Col
     Type(Mat_Morse_Entry), dimension(:), pointer   :: Col
  End Type Mat_Morse_Row

!!$  Type Mat_Morse_Col
!!$     Integer(Kind=Ks)                                        :: Nb_Row
!!$     Type(Mat_Morse_Row_ptr), Dimension(:), Pointer  :: Row_Ptr
!!$  End Type Mat_Morse_Col
!!$
!!$  Type Mat_Morse_Row_Ptr
!!$     ! Gives a fast way of scanning the matrix by column
!!$     Integer                                         :: Row_Number
!!$     Integer(Kind=Ks)                                        :: Row_Index
!!$  End Type Mat_Morse_Row_Ptr
  
  Type Mat_Morse
     Type(Mat_Morse_Row), dimension(:), pointer     :: Row
!!$     Type(Mat_Morse_Col_Ptr), Dimension(:), Pointer  :: Col
  End Type Mat_Morse

! Surcharge des operations arithmetiques
  Interface assignment (=)
     Module Procedure Mat_Prof_EQ
  End Interface

Contains
  Subroutine Mat_Prof_EQ(MP1,MP2)
    Type (Mat_Prof), Intent(OUT)                    :: MP1
    Type (Mat_Prof), Intent(IN)                     :: MP2

    MP1%Profil = MP2%Profil
    MP1%Val    = MP2%Val
  End Subroutine Mat_Prof_EQ
End Module m_Mef_types
