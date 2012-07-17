Module m_MEF_Elements
#include "finclude/petscdef.h"
   Use m_MEF_LinAlg
   Use m_MEF_Utils
   Use m_MEF_Parameters
   Use petsc
   use,intrinsic :: iso_c_binding
   IMPLICIT NONE

   !Private   
   ! not sure how to make the enumerator public short of listing them one after each other.
   Public :: ElementInit
   Public :: ElementDestroy
   Public :: ElementView
   Public :: Element_Type
   Public :: Element1D
   Public :: Element2D_Vect,Element2D_Scal,Element2D_Elast 
   Public :: Element3D_Vect,Element3D_Scal,Element3D_Elast 
   Public :: Element_TypeFindByID,Element_TypeFindByName
   Public :: EXO2Element_TypeScal,EXO2Element_TypeVect,EXO2Element_TypeElast
   
   Type Element_Type
      ! name is the element name in english language
      Character(len=MEF90_MXSTRLEN)                :: Name
      ! shortID is a numerical shortcut used in the prep, for instance
      PetscInt                                     :: ShortID
      
      !!! Geometry of the cell associated with the element
      PetscInt                                     :: numVertex
      PetscInt                                     :: numEdge
      PetscInt                                     :: numFaces
      
      !!! Number of dof at each vertex edge face cell
      PetscInt                                     :: numVertexDof
      PetscInt                                     :: numEdgeDof
      PetscInt                                     :: numFaceDof
      PetscInt                                     :: numCellDof
      PetscInt                                     :: numDoF 
      !!! numDof = numVertexDof * numVertex + numEdgeDof * numEdge 
      !!!          + numFaceDof * numFace + numCellDof
      PetscInt                                     :: dim
      !!! The dimension of the cell associated with the element
      PetscInt                                     :: coDim
      !!! The co-dimension of the cell associated with the element
      PetscInt                                     :: order
      !!! The polynomial order
   End Type Element_Type
 
   Enum,bind(c)
      enumerator  :: &
         MEF90_P1_Lagrange_2D_Scal_ShortID = 1,       &  ! 1
         MEF90_P1_Lagrange_3D_Scal_ShortID,           &  ! 2
         MEF90_P1_Lagrange_2D_Elast_ShortID,          &  ! 3
         MEF90_P1_Lagrange_3D_Elast_ShortID,          &  ! 4
         MEF90_P1_Lagrange_2D_Vect_ShortID,           &  ! 5
         MEF90_P1_Lagrange_3D_Vect_ShortID,           &  ! 6
         MEF90_P1_Lagrange_2DBoundary_Scal_ShortID,   &  ! 7
         MEF90_P1_Lagrange_3DBoundary_Scal_ShortID,   &  ! 8
         MEF90_P1_Lagrange_2DBoundary_Elast_ShortID,  &  ! 9
         MEF90_P1_Lagrange_3DBoundary_Elast_ShortID,  &  ! 10
         MEF90_P1_Lagrange_2DBoundary_Vect_ShortID,   &  ! 11 
         MEF90_P1_Lagrange_3DBoundary_Vect_ShortID,   &  ! 12 
         MEF90_P2_Lagrange_2D_Scal_ShortID,           &  ! 13
         MEF90_P2_Lagrange_3D_Scal_ShortID,           &  ! 14
         MEF90_P2_Lagrange_2D_Elast_ShortID,          &  ! 15
         MEF90_P2_Lagrange_3D_Elast_ShortID,          &  ! 16
         MEF90_P2_Lagrange_2D_Vect_ShortID,           &  ! 17
         MEF90_P2_Lagrange_3D_Vect_ShortID,           &  ! 18
         MEF90_P2_Lagrange_2DBoundary_Scal_ShortID,   &  ! 19
         MEF90_P2_Lagrange_3DBoundary_Scal_ShortID,   &  ! 20
         MEF90_P2_Lagrange_2DBoundary_Elast_ShortID,  &  ! 21
         MEF90_P2_Lagrange_3DBoundary_Elast_ShortID,  &  ! 22
         MEF90_P2_Lagrange_2DBoundary_Vect_ShortID,   &  ! 23 
         MEF90_P2_Lagrange_3DBoundary_Vect_ShortID       ! 24 
   End Enum      

   !!! 
   !!! Linear Lagrange elements
   !!!
   Type(Element_Type),Parameter,Public :: MEF90_P1_Lagrange_2D_Scal = Element_Type(   &
      "MEF90_P1_Lagrange_2D_Scal",        &  ! name
      MEF90_P1_Lagrange_2D_Scal_ShortID,  &  ! shortID
      3,3,0,                              &  ! numVertex,numEdge,numFace
      1,0,0,0,3,                          &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
      2,0,1                               &  ! dim,codim,order                             
   )
   Type(Element_Type),Parameter,Public :: MEF90_P1_Lagrange_3D_Scal = Element_Type(   &
      "MEF90_P1_Lagrange_3D_Scal",        &  ! name
      MEF90_P1_Lagrange_3D_Scal_ShortID,  &  ! shortID
      4,6,4,                              &  ! numVertex,numEdge,numFace
      1,0,0,0,4,                          &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
      3,0,1                               &  ! dim,codim,order                         
   )
   Type(Element_Type),Parameter,Public :: MEF90_P1_Lagrange_2D_Elast = Element_Type(   &
      "MEF90_P1_Lagrange_2D_Elast",       &  ! name
      MEF90_P1_Lagrange_2D_Elast_ShortID, &  ! shortID
      3,3,0,                              &  ! numVertex,numEdge,numFace
      2,0,0,0,6,                          &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
      2,0,1                               &  ! dim,codim,order                             
   )
   Type(Element_Type),Parameter,Public :: MEF90_P1_Lagrange_3D_Elast = Element_Type(   &
      "MEF90_P1_Lagrange_3D_Elast",       &  ! name
      MEF90_P1_Lagrange_3D_Elast_ShortID, &  ! shortID
      4,6,4,                              &  ! numVertex,numEdge,numFace
      3,0,0,0,12,                         &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
      3,0,1                               &  ! dim,codim,order                             
   )
   Type(Element_Type),Parameter,Public :: MEF90_P1_Lagrange_2D_Vect = Element_Type(   &
      "MEF90_P1_Lagrange_2D_Vect",        &  ! name
      MEF90_P1_Lagrange_2D_Vect_ShortID,  &  ! shortID
      3,3,0,                              &  ! numVertex,numEdge,numFace
      2,0,0,0,6,                          &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
      2,0,1                               &  ! dim,codim,order                             
   )
   Type(Element_Type),Parameter,Public :: MEF90_P1_Lagrange_3D_Vect = Element_Type(   &
      "MEF90_P1_Lagrange_3D_Vect",        &  ! name
      MEF90_P1_Lagrange_2D_Vect_ShortID,  &  ! shortID
      4,6,4,                              &  ! numVertex,numEdge,numFace
      3,0,0,0,12,                         &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
      3,0,1                               &  ! dim,codim,order                             
   )
   Type(Element_Type),Parameter,Public :: MEF90_P1_Lagrange_2DBoundary_Scal = Element_Type(   &
      "MEF90_P1_Lagrange_2DBoundary_Scal",         &  ! name
      MEF90_P1_Lagrange_2DBoundary_Scal_ShortID,   &  ! shortID
      2,1,0,                                       &  ! numVertex,numEdge,numFace
      1,0,0,0,2,                                   &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
      2,1,1                                        &  ! dim,codim,order                             
   )
   Type(Element_Type),Parameter,Public :: MEF90_P1_Lagrange_3DBoundary_Scal = Element_Type(   &
      "MEF90_P1_Lagrange_3DBoundary_Scal",         &  ! name
      MEF90_P1_Lagrange_3DBoundary_Scal_ShortID,   &  ! shortID
      3,3,0,                                       &  ! numVertex,numEdge,numFace
      1,0,0,0,3,                                   &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
      3,1,1                                        &  ! dim,codim,order                             
   )
   Type(Element_Type),Parameter,Public :: MEF90_P1_Lagrange_2DBoundary_Vect = Element_Type(   &
      "MEF90_P1_Lagrange_2DBoundary_Vect",         &  ! name
      MEF90_P1_Lagrange_2DBoundary_Vect_ShortID,   &  ! shortID
      2,1,0,                                       &  ! numVertex,numEdge,numFace
      2,0,0,0,4,                                   &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
      2,1,1                                        &  ! dim,codim,order                             
   )
   Type(Element_Type),Parameter,Public :: MEF90_P1_Lagrange_3DBoundary_Vect = Element_Type(   &
      "MEF90_P1_Lagrange_3DBoundary_Vect",         &  ! name
      MEF90_P1_Lagrange_3DBoundary_Vect_ShortID,   &  ! shortID
      3,3,0,                                       &  ! numVertex,numEdge,numFace
      3,0,0,0,9,                                   &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
      3,1,1                                        &  ! dim,codim,order                             
   )
   Type(Element_Type),Parameter,Public :: MEF90_P1_Lagrange_2DBoundary_Elast = Element_Type(   &
      "MEF90_P1_Lagrange_2DBoundary_Elast",        &  ! name
      MEF90_P1_Lagrange_2DBoundary_Elast_ShortID,  &  ! shortID
      2,1,0,                                       &  ! numVertex,numEdge,numFace
      2,0,0,0,4,                                   &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
      2,1,1                                        &  ! dim,codim,order                             
   )
   Type(Element_Type),Parameter,Public :: MEF90_P1_Lagrange_3DBoundary_Elast = Element_Type(   &
      "MEF90_P1_Lagrange_3DBoundary_Elast",        &  ! name
      MEF90_P1_Lagrange_3DBoundary_Elast_ShortID,  &  ! shortID
      3,3,0,                                       &  ! numVertex,numEdge,numFace
      3,0,0,0,9,                                   &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
      3,1,1                                        &  ! dim,codim,order                             
   )
   !!!
   !!! Quadratic Lagrange elements
   !!!
   Type(Element_Type),Parameter,Public :: MEF90_P2_Lagrange_2D_Scal = Element_Type(   &
      "MEF90_P2_Lagrange_2D_Scal",        &  ! name
      MEF90_P2_Lagrange_2D_Scal_ShortID,  &  ! shortID
      3,3,0,                              &  ! numVertex,numEdge,numFace
      3,3,0,0,6,                          &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
      2,0,2                               &  ! dim,codim,order                             
   )
   Type(Element_Type),Parameter,Public :: MEF90_P2_Lagrange_3D_Scal = Element_Type(   &
      "MEF90_P2_Lagrange_3D_Scal",        &  ! name
      MEF90_P2_Lagrange_3D_Scal_ShortID,  &  ! shortID
      4,6,4,                              &  ! numVertex,numEdge,numFace
      4,6,0,0,10,                         &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
      3,0,2                               &  ! dim,codim,order                             
   )
   Type(Element_Type),Parameter,Public :: MEF90_P2_Lagrange_2D_Elast = Element_Type(   &
      "MEF90_P2_Lagrange_2D_Elast",       &  ! name
      MEF90_P2_Lagrange_2D_Elast_ShortID, &  ! shortID
      3,3,0,                              &  ! numVertex,numEdge,numFace
      6,6,0,0,12,                         &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
      2,0,2                               &  ! dim,codim,order                             
   )
   Type(Element_Type),Parameter,Public :: MEF90_P2_Lagrange_3D_Elast = Element_Type(   &
      "MEF90_P2_Lagrange_3D_Elast",       &  ! name
      MEF90_P2_Lagrange_3D_Elast_ShortID, &  ! shortID
      4,6,4,                              &  ! numVertex,numEdge,numFace
      8,12,0,0,20,                        &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
      3,0,2                               &  ! dim,codim,order                             
   )
   Type(Element_Type),Parameter,Public :: MEF90_P2_Lagrange_2D_Vect = Element_Type(   &
      "MEF90_P2_Lagrange_2D_Vect",        &  ! name
      MEF90_P2_Lagrange_2D_Vect_ShortID,  &  ! shortID
      3,3,0,                              &  ! numVertex,numEdge,numFace
      6,6,0,0,12,                         &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
      2,0,2                               &  ! dim,codim,order                             
   )
   Type(Element_Type),Parameter,Public :: MEF90_P2_Lagrange_3D_Vect = Element_Type(   &
      "MEF90_P2_Lagrange_3D_Vect",        &  ! name
      MEF90_P2_Lagrange_3D_Vect_ShortID,  &  ! shortID
      4,6,4,                              &  ! numVertex,numEdge,numFace
      8,12,0,0,20,                        &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
      3,0,2                               &  ! dim,codim,order                             
   )
   Type(Element_Type),Parameter,Public :: MEF90_P2_Lagrange_2DBoundary_Scal = Element_Type(   &
      "MEF90_P2_Lagrange_2DBoundary_Scal",         &  ! name
      MEF90_P2_Lagrange_2DBoundary_Scal_ShortID,   &  ! shortID
      2,1,0,                                       &  ! numVertex,numEdge,numFace
      1,1,0,0,3,                                   &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
      2,1,2                                        &  ! dim,codim,order                             
   )
   Type(Element_Type),Parameter,Public :: MEF90_P2_Lagrange_3DBoundary_Scal = Element_Type(   &
      "MEF90_P2_Lagrange_3DBoundary_Scal",         &  ! name
      MEF90_P2_Lagrange_3DBoundary_Scal_ShortID,   &  ! shortID
      3,3,0,                                       &  ! numVertex,numEdge,numFace
      1,1,0,0,6,                                   &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
      3,1,2                                        &  ! dim,codim,order                             
   )
   Type(Element_Type),Parameter,Public :: MEF90_P2_Lagrange_2DBoundary_Vect = Element_Type(   &
      "MEF90_P2_Lagrange_2DBoundary_Vect",         &  ! name
      MEF90_P2_Lagrange_2DBoundary_Vect_ShortID,   &  ! shortID
      2,1,0,                                       &  ! numVertex,numEdge,numFace
      2,2,0,0,6,                                   &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
      2,1,2                                        &  ! dim,codim,order                             
   )
   Type(Element_Type),Parameter,Public :: MEF90_P2_Lagrange_3DBoundary_Vect = Element_Type(   &
      "MEF90_P2_Lagrange_3DBoundary_Vect",         &  ! name
      MEF90_P2_Lagrange_3DBoundary_Vect_ShortID,   &  ! shortID
      3,3,0,                                       &  ! numVertex,numEdge,numFace
      3,3,0,0,18,                                  &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
      3,1,2                                        &  ! dim,codim,order                             
   )
   Type(Element_Type),Parameter,Public :: MEF90_P2_Lagrange_2DBoundary_Elast = Element_Type(   &
      "MEF90_P2_Lagrange_2DBoundary_Elast",        &  ! name
      MEF90_P2_Lagrange_2DBoundary_Elast_ShortID,  &  ! shortID
      2,1,0,                                       &  ! numVertex,numEdge,numFace
      2,2,0,0,6,                                   &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
      2,1,2                                        &  ! dim,codim,order                             
   )
   Type(Element_Type),Parameter,Public :: MEF90_P2_Lagrange_3DBoundary_Elast = Element_Type(   &
      "MEF90_P2_Lagrange_3DBoundary_Elast",        &  ! name
      MEF90_P2_Lagrange_3DBoundary_Elast_ShortID,  &  ! shortID
      3,3,0,                                       &  ! numVertex,numEdge,numFace
      3,3,0,0,18,                                  &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
      3,1,2                                        &  ! dim,codim,order                             
   )

   Integer,Parameter,Public :: MEF90_numKnownElements = 24       
   Type(Element_Type),dimension(MEF90_numKnownElements),Parameter,Public   :: MEF90_knownElements = (/ &
      MEF90_P1_Lagrange_2D_Scal,          &
      MEF90_P1_Lagrange_3D_Scal,          &
      MEF90_P1_Lagrange_2D_Elast,         & 
      MEF90_P1_Lagrange_3D_Elast,         &
      MEF90_P1_Lagrange_2D_Vect,          &
      MEF90_P1_Lagrange_3D_Vect,          &
      MEF90_P1_Lagrange_2DBoundary_Scal,  &
      MEF90_P1_Lagrange_3DBoundary_Scal,  &
      MEF90_P1_Lagrange_2DBoundary_Elast, &
      MEF90_P1_Lagrange_3DBoundary_Elast, &
      MEF90_P1_Lagrange_2DBoundary_Vect,  &
      MEF90_P1_Lagrange_3DBoundary_Vect,  &
      MEF90_P2_Lagrange_2D_Scal,          &
      MEF90_P2_Lagrange_3D_Scal,          &
      MEF90_P2_Lagrange_2D_Vect,          &
      MEF90_P2_Lagrange_3D_Vect,          &
      MEF90_P2_Lagrange_2D_Elast,         &
      MEF90_P2_Lagrange_3D_Elast,         &
      MEF90_P2_Lagrange_2DBoundary_Scal,  &
      MEF90_P2_Lagrange_3DBoundary_Scal,  &
      MEF90_P2_Lagrange_2DBoundary_Elast, &
      MEF90_P2_Lagrange_3DBoundary_Elast, &
      MEF90_P2_Lagrange_2DBoundary_Vect,  &
      MEF90_P2_Lagrange_3DBoundary_Vect   &
   /)

   Character(kind=c_char,len=MEF90_MXSTRLEN),dimension(MEF90_numKnownElements+3),Parameter,Public   :: MEF90_knownElementNames = (/ &
      "P1_Lagrange_2D_Scal         ",     &
      "P1_Lagrange_3D_Scal         ",     &
      "P1_Lagrange_2D_Elast        ",     &
      "P1_Lagrange_3D_Elast        ",     &
      "P1_Lagrange_2D_Vect         ",     &
      "P1_Lagrange_3D_Vect         ",     &
      "P1_Lagrange_2DBoundary_Scal ",     &
      "P1_Lagrange_3DBoundary_Scal ",     &
      "P1_Lagrange_2DBoundary_Elast",     &
      "P1_Lagrange_3DBoundary_Elast",     &
      "P1_Lagrange_2DBoundary_Vect ",     &
      "P1_Lagrange_3DBoundary_Vect ",     &
      "P2_Lagrange_2D_Scal         ",     &
      "P2_Lagrange_3D_Scal         ",     &
      "P2_Lagrange_2D_Elast        ",     &
      "P2_Lagrange_3D_Elast        ",     &
      "P2_Lagrange_2D_Vect         ",     &
      "P2_Lagrange_3D_Vect         ",     &
      "P2_Lagrange_2DBoundary_Scal ",     &
      "P2_Lagrange_3DBoundary_Scal ",     &
      "P2_Lagrange_2DBoundary_Elast",     &
      "P2_Lagrange_3DBoundary_Elast",     &
      "P2_Lagrange_2DBoundary_Vect ",     &
      "P2_Lagrange_3DBoundary_Vect ",     &
      "MEF90_knownElementNames     ",     &
      "prefix_                     ",     &
      C_NULL_CHAR//"                           "/)

      
   Type Element1D
      PetscReal,Dimension(:,:),Pointer             :: BF
      PetscReal,Dimension(:,:),Pointer             :: Der_BF
      PetscReal,Dimension(:),Pointer               :: Gauss_C
   End Type Element1D
 
   Type Element2D_Scal
      PetscReal,Dimension(:,:),Pointer             :: BF
      Type(Vect2D),Dimension(:,:),Pointer          :: Grad_BF
      PetscReal,Dimension(:),Pointer               :: Gauss_C
   End Type Element2D_Scal
 
   Type Element2D_Vect
      Type (Vect2D),Dimension(:,:),Pointer         :: BF
      Type (Mat2D),Dimension(:,:),Pointer          :: Der_BF
      PetscReal,Dimension(:),Pointer               :: Gauss_C
   End Type Element2D_Vect
 
   Type Element2D_Elast
      Type (Vect2D),Dimension(:,:),Pointer         :: BF
      Type (MatS2D),Dimension(:,:),Pointer         :: GradS_BF
      PetscReal,Dimension(:),Pointer               :: Gauss_C
   End Type Element2D_Elast
 
   Type Element3D_Vect
      Type (Vect3D),Dimension(:,:),Pointer         :: BF
      Type (Mat3D),Dimension(:,:),Pointer          :: Der_BF
      PetscReal,Dimension(:),Pointer               :: Gauss_C
   End Type Element3D_Vect
 
   Type Element3D_Scal
      PetscReal,Dimension(:,:),Pointer             :: BF
      Type (Vect3D),Dimension(:,:),Pointer         :: Grad_BF
      PetscReal,Dimension(:),Pointer               :: Gauss_C
   End Type Element3D_Scal
 
   Type Element3D_Elast
      Type (Vect3D),Dimension(:,:),Pointer         :: BF
      Type (MatS3D),Dimension(:,:),Pointer         :: GradS_BF
      PetscReal,Dimension(:),Pointer               :: Gauss_C
   End Type Element3D_Elast

   Interface ElementInit
      Module Procedure Element2D_Scal_Init,Element2D_Vect_Init,Element2D_Elast_Init, &
                       Element3D_Scal_Init,Element3D_Vect_Init,Element3D_Elast_Init, &
                       Element2D_Scal_InitSet,Element2D_Vect_InitSet,Element2D_Elast_InitSet, &
                       Element3D_Scal_InitSet,Element3D_Vect_InitSet,Element3D_Elast_InitSet, &
                       Element2D_Scal_InitSet_ByShortID,Element2D_Vect_InitSet_ByShortID,Element2D_Elast_InitSet_ByShortID,&
                       Element3D_Scal_InitSet_ByShortID,Element3D_Vect_InitSet_ByShortID,Element3D_Elast_InitSet_ByShortID
   End Interface ElementInit
   
   Interface ElementDestroy
      Module Procedure Element2D_Scal_Destroy,Element2D_Vect_Destroy,Element2D_Elast_Destroy,&
                       Element3D_Scal_Destroy,Element3D_Vect_Destroy,Element3D_Elast_Destroy
   End Interface ElementDestroy
   
   Interface ElementView
      Module Procedure Element2D_Scal_View,Element2D_Vect_View,Element2D_Elast_View, &
                       Element3D_Scal_View,Element3D_Vect_View,Element3D_Elast_View
   End Interface ElementView
   
Contains
#undef __FUNCT__
#define __FUNCT__ "EXO2Element_TypeScal"
   Subroutine EXO2Element_TypeScal(exoName,dim,elemType)
      Character(len=*),Intent(IN)                 :: EXOName
      PetscInt,Intent(IN)                         :: dim
      Type(Element_Type),Intent(OUT)              :: elemType

      Select Case(trim(exoName))
         Case("TETRA","TETRA4")
            elemType = MEF90_P1_Lagrange_3D_Scal
         Case("TETRA10")
            elemType = MEF90_P2_Lagrange_3D_Scal
         Case("TRI","TRI3","TRISHELL","TRISHELL3")
            If (dim == 2) Then
               elemType = MEF90_P1_Lagrange_2D_Scal
            Else
               elemType = MEF90_P1_Lagrange_3DBoundary_Scal
            End If
         Case("TRI6","TRISHELL6")
            If (dim == 2) Then
               elemType = MEF90_P2_Lagrange_2D_Scal
            Else
               elemType = MEF90_P2_Lagrange_3DBoundary_Scal
            End If
         !Case("QUAD","QUAD4","SHELL","SHELL4")
         !   If (dim == 2) Then
         !      elemType = MEF90_Q1_Lagrange_2D_Scal
         !   Else
         !      elemType = MEF90_Q1_Lagrange_3DBoundary_Scal
         !   End If
         !Case("QUAD9","SHELL9")
         !   If (dim == 2) Then
         !      elemType = MEF90_Q2_Lagrange_2D_Scal
         !   Else
         !      elemType = MEF90_Q2_Lagrange_3DBoundary_Scal
         !   End If
         Case("BAR","BAR2")
            elemType = MEF90_P1_Lagrange_2DBoundary_Scal
         Case("BAR3")
            elemType = MEF90_P2_Lagrange_2DBoundary_Scal
         Case default
            Write(*,*),__FUNCT__,': Element ',trim(exoName),'not recognized. Set type manually.'
      End Select
   End Subroutine EXO2Element_TypeScal

#undef __FUNCT__
#define __FUNCT__ "EXO2Element_TypeVect"
   Subroutine EXO2Element_TypeVect(exoName,dim,elemType)
      Character(len=*),Intent(IN)                 :: EXOName
      PetscInt,Intent(IN)                         :: dim
      Type(Element_Type),Intent(OUT)              :: elemType

      Select Case(trim(exoName))
         Case("TETRA","TETRA4")
            elemType = MEF90_P1_Lagrange_3D_Vect
         Case("TETRA10")
            elemType = MEF90_P2_Lagrange_3D_Vect
         Case("TRI","TRI3","TRISHELL","TRISHELL3")
            If (dim == 2) Then
               elemType = MEF90_P1_Lagrange_2D_Vect
            Else
               elemType = MEF90_P1_Lagrange_3DBoundary_Vect
            End If
         Case("TRI6","TRISHELL6")
            If (dim == 2) Then
               elemType = MEF90_P2_Lagrange_2D_Vect
            Else
               elemType = MEF90_P2_Lagrange_3DBoundary_Vect
            End If
         !Case("QUAD","QUAD4","SHELL","SHELL4")
            !If (dim == 2) Then
            !   elemType = MEF90_Q1_Lagrange_2D_Vect
            !Else
            !   elemType = MEF90_Q1_Lagrange_3DBoundary_Vect
            !End If
         !Case("QUAD9","SHELL9")
            !If (dim == 2) Then
            !   elemType = MEF90_Q2_Lagrange_2D_Vect
            !Else
            !   elemType = MEF90_Q2_Lagrange_3DBoundary_Vect
            !End If
         Case("BAR","BAR2")
            elemType = MEF90_P1_Lagrange_2DBoundary_Vect
         Case("BAR3")
            elemType = MEF90_P2_Lagrange_2DBoundary_Vect
         Case default
            Write(*,*),__FUNCT__,': Element ',trim(exoName),'not recognized. Set type manually.'
      End Select
   End Subroutine EXO2Element_TypeVect

#undef __FUNCT__
#define __FUNCT__ "EXO2Element_TypeElast"
   Subroutine EXO2Element_TypeElast(exoName,dim,elemType)
      Character(len=*),Intent(IN)                 :: EXOName
      PetscInt,Intent(IN)                         :: dim
      Type(Element_Type),Intent(OUT)              :: elemType

      Select Case(trim(exoName))
         Case("TETRA","TETRA4")
            elemType = MEF90_P1_Lagrange_3D_Elast
         Case("TETRA10")
            elemType = MEF90_P2_Lagrange_3D_Elast
         Case("TRI","TRI3","TRISHELL","TRISHELL3")
            If (dim == 2) Then
               elemType = MEF90_P1_Lagrange_2D_Elast
            Else
               elemType = MEF90_P1_Lagrange_3DBoundary_Elast
            End If
         Case("TRI6","TRISHELL6")
            If (dim == 2) Then
               elemType = MEF90_P2_Lagrange_2D_Elast
            Else
               elemType = MEF90_P2_Lagrange_3DBoundary_Elast
            End If
         !Case("QUAD","QUAD4","SHELL","SHELL4")
            !If (dim == 2) Then
            !   elemType = MEF90_Q1_Lagrange_2D_Elast
            !Else
            !   elemType = MEF90_Q1_Lagrange_3DBoundary_Elast
            !End If
         !Case("QUAD9","SHELL9")
            !If (dim == 2) Then
            !   elemType = MEF90_Q2_Lagrange_2D_Elast
            !Else
            !   elemType = MEF90_Q2_Lagrange_3DBoundary_Elast
            !End If
         Case("BAR","BAR2")
            elemType = MEF90_P1_Lagrange_2DBoundary_Elast
         Case("BAR3")
            elemType = MEF90_P2_Lagrange_2DBoundary_Elast
         Case default
            Write(*,*),__FUNCT__,': Element ',trim(exoName),'not recognized. Set type manually.'
      End Select
   End Subroutine EXO2Element_TypeElast

#undef __FUNCT__
#define __FUNCT__ "Element_TypeFindByID"
   Subroutine Element_TypeFindByID(elemID,elemType)
      PetscInt, Intent(IN)                        :: elemID
      Type(Element_Type),Intent(OUT)              :: elemType
      
      Integer                                     :: i
      PetscBool                                   :: knownID = PETSC_FALSE      
      Do i = 1, size(MEF90_knownElements)
         If (MEF90_knownElements(i)%shortID == elemID) Then
            elemType = MEF90_knownElements(i)
            knownID  = PETSC_TRUE
            EXIT
         End If
      End Do
      If (.NOT. knownID) Then
         Write(*,*) "[ERROR]: Unknown element ID", elemID
      End If
   End Subroutine Element_TypeFindByID
   
#undef __FUNCT__
#define __FUNCT__ "Element_TypeFindByname"
   Subroutine Element_TypeFindByname(elemName,elemType)
      Character(len=*), Intent(IN)                :: elemName
      Type(Element_Type),Intent(OUT)              :: elemType

      Integer                                     :: i
      PetscBool                                   :: knownID = PETSC_FALSE      
      Do i = 1, size(MEF90_knownElements)
         If (trim(MEF90_knownElements(i)%name) == trim(elemName)) Then
            elemType = MEF90_knownElements(i)
            knownID  = PETSC_TRUE
            EXIT
         End If
      End Do
      If (.NOT. knownID) Then
         Write(*,*) "[ERROR]: Unknown element ", trim(elemName)
      End If
   End Subroutine Element_TypeFindByname

#undef __FUNCT__
#define __FUNCT__ "Element2D_Scal_InitSet_ByName"
   Subroutine Element2D_Scal_InitSet_ByName(mesh,cellIS,dElem,dQuadratureOrder,elemName)
      Type(DM),intent(IN)                         :: mesh
      Type(IS),intent(IN)                         :: cellIS
      Type(Element2D_Scal),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                         :: dQuadratureOrder
      Character(len=*), Intent(IN)                :: elemName
      
      Type(Element_Type)                          :: elemType

      Call Element_TypeFindByName(elemName,elemType)
      Call Element2D_Scal_InitSet(mesh,cellIS,dElem,dQuadratureOrder,elemType)
   End Subroutine Element2D_Scal_InitSet_ByName

#undef __FUNCT__
#define __FUNCT__ "Element2D_Vect_InitSet_ByName"
   Subroutine Element2D_Vect_InitSet_ByName(mesh,cellIS,dElem,dQuadratureOrder,elemName)
      Type(DM),intent(IN)                         :: mesh
      Type(IS),intent(IN)                         :: cellIS
      Type(Element2D_Vect),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                         :: dQuadratureOrder
      Character(len=*), Intent(IN)                :: elemName
      
      Type(Element_Type)                          :: elemType

      Call Element_TypeFindByName(elemName,elemType)
      Call Element2D_Vect_InitSet(mesh,cellIS,dElem,dQuadratureOrder,elemType)
   End Subroutine Element2D_Vect_InitSet_ByName

#undef __FUNCT__
#define __FUNCT__ "Element2D_Elast_InitSet_ByName"
   Subroutine Element2D_Elast_InitSet_ByName(mesh,cellIS,dElem,dQuadratureOrder,elemName)
      Type(DM),intent(IN)                         :: mesh
      Type(IS),intent(IN)                         :: cellIS
      Type(Element2D_Elast),Dimension(:),Pointer  :: dElem
      PetscInt,Intent(IN)                         :: dQuadratureOrder
      Character(len=*), Intent(IN)                :: elemName
      
      Type(Element_Type)                          :: elemType

      Call Element_TypeFindByName(elemName,elemType)
      Call Element2D_Elast_InitSet(mesh,cellIS,dElem,dQuadratureOrder,elemType)
   End Subroutine Element2D_Elast_InitSet_ByName

#undef __FUNCT__
#define __FUNCT__ "Element3D_Scal_InitSet_ByName"
   Subroutine Element3D_Scal_InitSet_ByName(mesh,cellIS,dElem,dQuadratureOrder,elemName)
      Type(DM),intent(IN)                         :: mesh
      Type(IS),intent(IN)                         :: cellIS
      Type(Element3D_Scal),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                         :: dQuadratureOrder
      Character(len=*), Intent(IN)                :: elemName
      
      Type(Element_Type)                          :: elemType

      Call Element_TypeFindByName(elemName,elemType)
      Call Element3D_Scal_InitSet(mesh,cellIS,dElem,dQuadratureOrder,elemType)
   End Subroutine Element3D_Scal_InitSet_ByName

#undef __FUNCT__
#define __FUNCT__ "Element3D_Vect_InitSet_ByName"
   Subroutine Element3D_Vect_InitSet_ByName(mesh,cellIS,dElem,dQuadratureOrder,elemName)
      Type(DM),intent(IN)                         :: mesh
      Type(IS),intent(IN)                         :: cellIS
      Type(Element3D_Vect),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                         :: dQuadratureOrder
      Character(len=*), Intent(IN)                :: elemName
      
      Type(Element_Type)                          :: elemType

      Call Element_TypeFindByName(elemName,elemType)
      Call Element3D_Vect_InitSet(mesh,cellIS,dElem,dQuadratureOrder,elemType)
   End Subroutine Element3D_Vect_InitSet_ByName

#undef __FUNCT__
#define __FUNCT__ "Element3D_Elast_InitSet_ByName"
   Subroutine Element3D_Elast_InitSet_ByName(mesh,cellIS,dElem,dQuadratureOrder,elemName)
      Type(DM),intent(IN)                         :: mesh
      Type(IS),intent(IN)                         :: cellIS
      Type(Element3D_Elast),Dimension(:),Pointer  :: dElem
      PetscInt,Intent(IN)                         :: dQuadratureOrder
      Character(len=*), Intent(IN)                :: elemName
      
      Type(Element_Type)                          :: elemType

      Call Element_TypeFindByName(elemName,elemType)
      Call Element3D_Elast_InitSet(mesh,cellIS,dElem,dQuadratureOrder,elemType)
   End Subroutine Element3D_Elast_InitSet_ByName

#undef __FUNCT__
#define __FUNCT__ "Element2D_Scal_InitSet_ByShortID"
   Subroutine Element2D_Scal_InitSet_ByShortID(mesh,cellIS,dElem,dQuadratureOrder,shortID)
      Type(DM),intent(IN)                         :: mesh
      Type(IS),intent(IN)                         :: cellIS
      Type(Element2D_Scal),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                         :: dQuadratureOrder,shortID
      
      Type(Element_Type)                          :: elemType

      Call Element_TypeFindByID(shortID,elemType)
      Call Element2D_Scal_InitSet(mesh,cellIS,dElem,dQuadratureOrder,elemType)
   End Subroutine Element2D_Scal_InitSet_ByShortID

#undef __FUNCT__
#define __FUNCT__ "Element2D_Vect_InitSet_ByShortID"
   Subroutine Element2D_Vect_InitSet_ByShortID(mesh,cellIS,dElem,dQuadratureOrder,shortID)
      Type(DM),intent(IN)                         :: mesh
      Type(IS),intent(IN)                         :: cellIS
      Type(Element2D_Vect),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                         :: dQuadratureOrder,shortID
      
      Type(Element_Type)                          :: elemType

      Call Element_TypeFindByID(shortID,elemType)
      Call Element2D_Vect_InitSet(mesh,cellIS,dElem,dQuadratureOrder,elemType)
   End Subroutine Element2D_Vect_InitSet_ByShortID

#undef __FUNCT__
#define __FUNCT__ "Element2D_Elast_InitSet_ByShortID"
   Subroutine Element2D_Elast_InitSet_ByShortID(mesh,cellIS,dElem,dQuadratureOrder,shortID)
      Type(DM),intent(IN)                         :: mesh
      Type(IS),intent(IN)                         :: cellIS
      Type(Element2D_Elast),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                         :: dQuadratureOrder,shortID
      
      Type(Element_Type)                          :: elemType

      Call Element_TypeFindByID(shortID,elemType)
      Call Element2D_Elast_InitSet(mesh,cellIS,dElem,dQuadratureOrder,elemType)
   End Subroutine Element2D_Elast_InitSet_ByShortID

#undef __FUNCT__
#define __FUNCT__ "Element3D_Scal_InitSet_ByShortID"
   Subroutine Element3D_Scal_InitSet_ByShortID(mesh,cellIS,dElem,dQuadratureOrder,shortID)
      Type(DM),intent(IN)                         :: mesh
      Type(IS),intent(IN)                         :: cellIS
      Type(Element3D_Scal),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                         :: dQuadratureOrder,shortID
      
      Type(Element_Type)                          :: elemType

      Call Element_TypeFindByID(shortID,elemType)
      Call Element3D_Scal_InitSet(mesh,cellIS,dElem,dQuadratureOrder,elemType)
   End Subroutine Element3D_Scal_InitSet_ByShortID

#undef __FUNCT__
#define __FUNCT__ "Element3D_Vect_InitSet_ByShortID"
   Subroutine Element3D_Vect_InitSet_ByShortID(mesh,cellIS,dElem,dQuadratureOrder,shortID)
      Type(DM),intent(IN)                         :: mesh
      Type(IS),intent(IN)                         :: cellIS
      Type(Element3D_Vect),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                         :: dQuadratureOrder,shortID
      
      Type(Element_Type)                          :: elemType

      Call Element_TypeFindByID(shortID,elemType)
      Call Element3D_Vect_InitSet(mesh,cellIS,dElem,dQuadratureOrder,elemType)
   End Subroutine Element3D_Vect_InitSet_ByShortID

#undef __FUNCT__
#define __FUNCT__ "Element3D_Elast_InitSet_ByShortID"
   Subroutine Element3D_Elast_InitSet_ByShortID(mesh,cellIS,dElem,dQuadratureOrder,shortID)
      Type(DM),intent(IN)                         :: mesh
      Type(IS),intent(IN)                         :: cellIS
      Type(Element3D_Elast),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                         :: dQuadratureOrder,shortID
      
      Type(Element_Type)                          :: elemType

      Call Element_TypeFindByID(shortID,elemType)
      Call Element3D_Elast_InitSet(mesh,cellIS,dElem,dQuadratureOrder,elemType)
   End Subroutine Element3D_Elast_InitSet_ByShortID

#undef __FUNCT__
#define __FUNCT__ "Element2D_Scal_InitSet"
   Subroutine Element2D_Scal_InitSet(mesh,cellIS,dElem,dQuadratureOrder,elemType)
      Type(DM),intent(IN)                         :: mesh
      Type(IS),intent(IN)                         :: cellIS
      Type(Element2D_Scal),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                         :: dQuadratureOrder
      Type(Element_Type),intent(IN)               :: elemType
      
      PetscInt                                    :: conesize,iELoc,ierr
      PetscInt                                    :: i,j,k
      PetscInt,Dimension(:),Pointer               :: CellID
      PetscReal,Dimension(:,:),Pointer            :: Coord
      PetscReal,Dimension(:),Pointer              :: TmpCoord
      Type(SectionReal)                           :: CoordSection
      
      Call DMMeshGetSectionReal(mesh,'coordinates',CoordSection,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellIS,CellID,ierr);CHKERRQ(ierr)
      
      If (size(CellID) > 0) Then
         iELoc = CellID(1)
         Call DMMeshGetConeSize(mesh,iELoc,coneSize,ierr);CHKERRQ(ierr)    
         !!! 
         !!! conesize*numdim is an upper bound on the size of the restriction
         !!! of the coordinate section to a cell (interpolated mesh) 
         Allocate(TmpCoord(conesize * elemType%dim))
         Allocate(Coord(elemType%dim,elemType%numVertex))
         Do_Elem_iE: Do iELoc = 1,size(CellID)
            Call SectionRealRestrictClosure(CoordSection,mesh,CellID(iELoc),Size(TmpCoord),TmpCoord,ierr);CHKERRQ(ierr)
            k = 1
            Do i = 1, elemType%numVertex
               Do j = 1, elemType%dim
                  Coord(j,i) = TmpCoord(k)
                  k = k+1
               End Do
            End Do
            Call Element2D_Scal_Init(dElem(CellID(iELoc)+1),Coord,dQuadratureOrder,elemType)
         End Do Do_Elem_iE
         DeAllocate(TmpCoord)
         DeAllocate(Coord)
      End If
      Call SectionRealDestroy(CoordSection,ierr);CHKERRQ(ierr)
      Call ISRestoreIndicesF90(CellIS,CellID,ierr);CHKERRQ(ierr)
   End Subroutine Element2D_Scal_InitSet

#undef __FUNCT__
#define __FUNCT__ "Element2D_Vect_InitSet"
   Subroutine Element2D_Vect_InitSet(mesh,cellIS,dElem,dQuadratureOrder,elemType)
      Type(DM),intent(IN)                         :: mesh
      Type(IS),intent(IN)                         :: cellIS
      Type(Element2D_Vect),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                         :: dQuadratureOrder
      Type(Element_Type),intent(IN)               :: elemType
      
      PetscInt                                    :: conesize,iELoc,ierr
      PetscInt                                    :: i,j,k
      PetscInt,Dimension(:),Pointer               :: CellID
      PetscReal,Dimension(:,:),Pointer            :: Coord
      PetscReal,Dimension(:),Pointer              :: TmpCoord
      Type(SectionReal)                           :: CoordSection
      
      Call DMMeshGetSectionReal(mesh,'coordinates',CoordSection,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellIS,CellID,ierr);CHKERRQ(ierr)

      Call DMMeshGetConeSize(mesh,CellID(1),coneSize,ierr);CHKERRQ(ierr)    
      !!! 
      !!! conesize*numdim is an upper bound on the size of the restriction
      !!! of the coordinate section to a cell (interpolated mesh) 
      Allocate(TmpCoord(conesize * elemType%dim))
      Allocate(Coord(elemType%dim,elemType%numVertex))
      Do_Elem_iE: Do iELoc = 1,size(CellID)
         Call SectionRealRestrictClosure(CoordSection,mesh,CellID(iELoc),Size(TmpCoord),TmpCoord,ierr);CHKERRQ(ierr)
         k = 1
         Do i = 1, elemType%numVertex
            Do j = 1, elemType%dim
               Coord(j,i) = TmpCoord(k)
               k = k+1
            End Do
         End Do
         Call Element2D_Vect_Init(dElem(cellID(iELoc)+1),Coord,dQuadratureOrder,elemType)
      End Do Do_Elem_iE
      DeAllocate(TmpCoord)
      DeAllocate(Coord)

      Call SectionRealDestroy(CoordSection,ierr);CHKERRQ(ierr)
      Call ISRestoreIndicesF90(CellIS,CellID,ierr);CHKERRQ(ierr)
   End Subroutine Element2D_Vect_InitSet

#undef __FUNCT__
#define __FUNCT__ "Element2D_Elast_InitSet"
   Subroutine Element2D_Elast_InitSet(mesh,cellIS,dElem,dQuadratureOrder,elemType)
      Type(DM),intent(IN)                         :: mesh
      Type(IS),intent(IN)                         :: cellIS
      Type(Element2D_Elast),Dimension(:),Pointer  :: dElem
      PetscInt,Intent(IN)                         :: dQuadratureOrder
      Type(Element_Type),intent(IN)               :: elemType
      
      PetscInt                                    :: conesize,iELoc,ierr
      PetscInt                                    :: i,j,k
      PetscInt,Dimension(:),Pointer               :: CellID
      PetscReal,Dimension(:,:),Pointer            :: Coord
      PetscReal,Dimension(:),Pointer              :: TmpCoord
      Type(SectionReal)                           :: CoordSection
      
      Call DMMeshGetSectionReal(mesh,'coordinates',CoordSection,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellIS,CellID,ierr);CHKERRQ(ierr)

      Call DMMeshGetConeSize(mesh,CellID(1),coneSize,ierr);CHKERRQ(ierr)    
      !!! 
      !!! conesize*numdim is an upper bound on the size of the restriction
      !!! of the coordinate section to a cell (interpolated mesh) 
      Allocate(TmpCoord(conesize * elemType%dim))
      Allocate(Coord(elemType%dim,elemType%numVertex))
      Do_Elem_iE: Do iELoc = 1,size(CellID)
         Call SectionRealRestrictClosure(CoordSection,mesh,CellID(iELoc),Size(TmpCoord),TmpCoord,ierr);CHKERRQ(ierr)
         k = 1
         Do i = 1, elemType%numVertex
            Do j = 1, elemType%dim
               Coord(j,i) = TmpCoord(k)
               k = k+1
            End Do
         End Do
         Call Element2D_Elast_Init(dElem(cellID(iELoc)+1),Coord,dQuadratureOrder,elemType)
      End Do Do_Elem_iE
      DeAllocate(TmpCoord)
      DeAllocate(Coord)

      Call SectionRealDestroy(CoordSection,ierr);CHKERRQ(ierr)
      Call ISRestoreIndicesF90(CellIS,CellID,ierr);CHKERRQ(ierr)
   End Subroutine Element2D_Elast_InitSet

#undef __FUNCT__
#define __FUNCT__ "Element3D_Scal_InitSet"
   Subroutine Element3D_Scal_InitSet(mesh,cellIS,dElem,dQuadratureOrder,elemType)
      Type(DM),intent(IN)                         :: mesh
      Type(IS),intent(IN)                         :: cellIS
      Type(Element3D_Scal),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                         :: dQuadratureOrder
      Type(Element_Type),intent(IN)               :: elemType
      
      PetscInt                                    :: conesize,iELoc,ierr
      PetscInt                                    :: i,j,k
      PetscInt,Dimension(:),Pointer               :: CellID
      PetscReal,Dimension(:,:),Pointer            :: Coord
      PetscReal,Dimension(:),Pointer              :: TmpCoord
      Type(SectionReal)                           :: CoordSection
      
      Call DMMeshGetSectionReal(mesh,'coordinates',CoordSection,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellIS,CellID,ierr);CHKERRQ(ierr)

      Call DMMeshGetConeSize(mesh,CellID(1),coneSize,ierr);CHKERRQ(ierr)    
      !!! 
      !!! conesize*numdim is an upper bound on the size of the restriction
      !!! of the coordinate section to a cell (interpolated mesh)
      Allocate(TmpCoord(conesize * elemType%dim))
      Allocate(Coord(elemType%dim,elemType%numVertex))
      Do_Elem_iE: Do iELoc = 1,size(CellID)
         Call SectionRealRestrictClosure(CoordSection,mesh,CellID(iELoc),Size(TmpCoord),TmpCoord,ierr);CHKERRQ(ierr)
         k = 1
         Do i = 1, elemType%numVertex
            Do j = 1, elemType%dim
               Coord(j,i) = TmpCoord(k)
               k = k+1
            End Do
         End Do
         Call Element3D_Scal_Init(dElem(CellID(iELoc)+1),Coord,dQuadratureOrder,elemType)
      End Do Do_Elem_iE
      DeAllocate(TmpCoord)
      DeAllocate(Coord)

      Call SectionRealDestroy(CoordSection,ierr);CHKERRQ(ierr)
      Call ISRestoreIndicesF90(CellIS,CellID,ierr);CHKERRQ(ierr)
   End Subroutine Element3D_Scal_InitSet

#undef __FUNCT__
#define __FUNCT__ "Element3D_Vect_InitSet"
   Subroutine Element3D_Vect_InitSet(mesh,cellIS,dElem,dQuadratureOrder,elemType)
      Type(DM),intent(IN)                         :: mesh
      Type(IS),intent(IN)                         :: cellIS
      Type(Element3D_Vect),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                         :: dQuadratureOrder
      Type(Element_Type),intent(IN)               :: elemType
      
      PetscInt                                    :: conesize,iELoc,ierr
      PetscInt                                    :: i,j,k
      PetscInt,Dimension(:),Pointer               :: CellID
      PetscReal,Dimension(:,:),Pointer            :: Coord
      PetscReal,Dimension(:),Pointer              :: TmpCoord
      Type(SectionReal)                           :: CoordSection
      
      Call DMMeshGetSectionReal(mesh,'coordinates',CoordSection,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellIS,CellID,ierr);CHKERRQ(ierr)

      Call DMMeshGetConeSize(mesh,CellID(1),coneSize,ierr);CHKERRQ(ierr)    
      !!! 
      !!! conesize*numdim is an upper bound on the size of the restriction
      !!! of the coordinate section to a cell (interpolated mesh) 
      Allocate(TmpCoord(conesize * elemType%dim))
      Allocate(Coord(elemType%dim,elemType%numVertex))
      Do_Elem_iE: Do iELoc = 1,size(CellID)
         Call SectionRealRestrictClosure(CoordSection,mesh,CellID(iELoc),Size(TmpCoord),TmpCoord,ierr);CHKERRQ(ierr)
         k = 1
         Do i = 1, elemType%numVertex
            Do j = 1, elemType%dim
               Coord(j,i) = TmpCoord(k)
               k = k+1
            End Do
         End Do
         Call Element3D_Vect_Init(dElem(cellID(iELoc)+1),Coord,dQuadratureOrder,elemType)
      End Do Do_Elem_iE
      DeAllocate(TmpCoord)
      DeAllocate(Coord)

      Call SectionRealDestroy(CoordSection,ierr);CHKERRQ(ierr)
      Call ISRestoreIndicesF90(CellIS,CellID,ierr);CHKERRQ(ierr)
   End Subroutine Element3D_Vect_InitSet

#undef __FUNCT__
#define __FUNCT__ "Element3D_Elast_InitSet"
   Subroutine Element3D_Elast_InitSet(mesh,cellIS,dElem,dQuadratureOrder,elemType)
      Type(DM),intent(IN)                         :: mesh
      Type(IS),intent(IN)                         :: cellIS
      Type(Element3D_Elast),Dimension(:),Pointer  :: dElem
      PetscInt,Intent(IN)                         :: dQuadratureOrder
      Type(Element_Type),intent(IN)               :: elemType
      
      PetscInt                                    :: conesize,iELoc,ierr
      PetscInt                                    :: i,j,k
      PetscInt,Dimension(:),Pointer               :: CellID
      PetscReal,Dimension(:,:),Pointer            :: Coord
      PetscReal,Dimension(:),Pointer              :: TmpCoord
      Type(SectionReal)                           :: CoordSection
      
      Call DMMeshGetSectionReal(mesh,'coordinates',CoordSection,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellIS,CellID,ierr);CHKERRQ(ierr)

      Call DMMeshGetConeSize(mesh,CellID(1),coneSize,ierr);CHKERRQ(ierr)    
      !!! 
      !!! conesize*numdim is an upper bound on the size of the restriction
      !!! of the coordinate section to a cell (interpolated mesh) 
      Allocate(TmpCoord(conesize * elemType%dim))
      Allocate(Coord(elemType%dim,elemType%numVertex))
      Do_Elem_iE: Do iELoc = 1,size(CellID)
         Call SectionRealRestrictClosure(CoordSection,mesh,CellID(iELoc),Size(TmpCoord),TmpCoord,ierr);CHKERRQ(ierr)
         k = 1
         Do i = 1, elemType%numVertex
            Do j = 1, elemType%dim
               Coord(j,i) = TmpCoord(k)
               k = k+1
            End Do
         End Do
         Call Element3D_Elast_Init(dElem(cellID(iELoc)+1),Coord,dQuadratureOrder,elemType)
      End Do Do_Elem_iE
      DeAllocate(TmpCoord)
      DeAllocate(Coord)

      Call SectionRealDestroy(CoordSection,ierr);CHKERRQ(ierr)
      Call ISRestoreIndicesF90(CellIS,CellID,ierr);CHKERRQ(ierr)
   End Subroutine Element3D_Elast_InitSet


#undef __FUNCT__
#define __FUNCT__ "Element2D_Scal_Init"
   Subroutine Element2D_Scal_Init(dElem,dCoord,QuadratureOrder,elemType)
      Type(Element2D_Scal)                   :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord
      PetscInt,Intent(IN)                    :: QuadratureOrder
      Type(Element_Type),intent(IN)          :: elemType
      
      Select Case (elemType%shortID)
         Case (MEF90_P1_Lagrange_2D_Scal%shortID)
            Call Element_P_Lagrange_2D_Scal_Init(dElem,dCoord,1,QuadratureOrder)
         Case (MEF90_P2_Lagrange_2D_Scal%shortID)
            Call Element_P_Lagrange_2D_Scal_Init(dElem,dCoord,2,QuadratureOrder)
         Case (MEF90_P1_Lagrange_2DBoundary_Scal%shortID)
            Call Element_P_Lagrange_2DBoundary_Scal_Init(dElem,dCoord,1,QuadratureOrder)
         Case (MEF90_P2_Lagrange_2DBoundary_Scal%shortID)
            Call Element_P_Lagrange_2DBoundary_Scal_Init(dElem,dCoord,2,QuadratureOrder)
!         Case (MEF90_Q1_Lagrange_2D_Scal%shortID)
!            Call Element_Q_Lagrange_2D_Scal_Init(dElem,dCoord,1,QuadratureOrder)
!         Case (MEF90_Q2_Lagrange_2D_Scal%shortID)
!            Call Element_Q_Lagrange_2D_Scal_Init(dElem,dCoord,2,QuadratureOrder)
         Case Default
            Print*,__FUNCT__,': Element type not implemented yet',elemType%name,elemType%shortID
      End Select
   End Subroutine Element2D_Scal_Init                                
   
#undef __FUNCT__
#define __FUNCT__ "Element2D_Vect_Init"
   Subroutine Element2D_Vect_Init(dElem,dCoord,QuadratureOrder,elemType)
      Type(Element2D_Vect)                   :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord
      PetscInt,Intent(IN)                    :: QuadratureOrder
      Type(Element_Type),intent(IN)          :: elemType
      
      Select Case (elemType%shortID)
         Case (MEF90_P1_Lagrange_2D_Vect%shortID)
            Call Element_P_Lagrange_2D_Vect_Init(dElem,dCoord,1,QuadratureOrder)
         Case (MEF90_P2_Lagrange_2D_Vect%shortID)
            Call Element_P_Lagrange_2D_Vect_Init(dElem,dCoord,2,QuadratureOrder)
         Case (MEF90_P1_Lagrange_2DBoundary_Vect%shortID)
            Call Element_P_Lagrange_2DBoundary_Vect_Init(dElem,dCoord,1,QuadratureOrder)
         Case (MEF90_P2_Lagrange_2DBoundary_Vect%shortID)
            Call Element_P_Lagrange_2DBoundary_Vect_Init(dElem,dCoord,2,QuadratureOrder)
!         Case (MEF90_Q1_Lagrange_2D_Vect%shortID)
!            Call Element_Q_Lagrange_2D_Vect_Init(dElem,dCoord,1,QuadratureOrder)
!         Case (MEF90_Q2_Lagrange_2D_Vect%shortID)
!            Call Element_Q_Lagrange_2D_Vect_Init(dElem,dCoord,2,QuadratureOrder)
         Case Default
            Print*,__FUNCT__,': Element type not implemented yet',elemType%name,elemType%shortID
      End Select
   End Subroutine Element2D_Vect_Init                                

#undef __FUNCT__
#define __FUNCT__ "Element2D_Elast_Init"
   Subroutine Element2D_Elast_Init(dElem,dCoord,QuadratureOrder,elemType)
      Type(Element2D_Elast)                   :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord
      PetscInt,Intent(IN)                    :: QuadratureOrder
      Type(Element_Type),intent(IN)          :: elemType
      
      Select Case (elemType%shortID)
         Case (MEF90_P1_Lagrange_2D_Elast%shortID)
            Call Element_P_Lagrange_2D_Elast_Init(dElem,dCoord,1,QuadratureOrder)
         Case (MEF90_P2_Lagrange_2D_Elast%shortID)
            Call Element_P_Lagrange_2D_Elast_Init(dElem,dCoord,2,QuadratureOrder)
         Case (MEF90_P1_Lagrange_2DBoundary_Elast%shortID)
            Call Element_P_Lagrange_2DBoundary_Elast_Init(dElem,dCoord,1,QuadratureOrder)
         Case (MEF90_P2_Lagrange_2DBoundary_Elast%shortID)
            Call Element_P_Lagrange_2DBoundary_Elast_Init(dElem,dCoord,2,QuadratureOrder)
!         Case (MEF90_Q1_Lagrange_2D_Elast%shortID)
!            Call Element_Q_Lagrange_2D_Elast_Init(dElem,dCoord,1,QuadratureOrder)
!         Case (MEF90_Q2_Lagrange_2D_Elast%shortID)
!            Call Element_Q_Lagrange_2D_Elast_Init(dElem,dCoord,2,QuadratureOrder)
         Case Default
            Print*,__FUNCT__,': Element type not implemented yet',elemType%name,elemType%shortID
      End Select
   End Subroutine Element2D_Elast_Init                                
   
#undef __FUNCT__
#define __FUNCT__ "Element3D_Scal_Init"
   Subroutine Element3D_Scal_Init(dElem,dCoord,QuadratureOrder,elemType)
      Type(Element3D_Scal)                   :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord
      PetscInt,Intent(IN)                    :: QuadratureOrder
      Type(Element_Type),intent(IN)          :: elemType
      
      Select Case (elemType%shortID)
         Case (MEF90_P1_Lagrange_3D_Scal%shortID)
            Call Element_P_Lagrange_3D_Scal_Init(dElem,dCoord,1,QuadratureOrder)
         Case (MEF90_P2_Lagrange_3D_Scal%shortID)
            Call Element_P_Lagrange_3D_Scal_Init(dElem,dCoord,2,QuadratureOrder)
         Case (MEF90_P1_Lagrange_3DBoundary_Scal%shortID)
            Call Element_P_Lagrange_3DBoundary_Scal_Init(dElem,dCoord,1,QuadratureOrder)
         Case (MEF90_P2_Lagrange_3DBoundary_Scal%shortID)
            Call Element_P_Lagrange_3DBoundary_Scal_Init(dElem,dCoord,2,QuadratureOrder)
!         Case (MEF90_Q1_Lagrange_3D_Scal%shortID)
!            Call Element_Q_Lagrange_3D_Scal_Init(dElem,dCoord,1,QuadratureOrder)
!         Case (MEF90_Q2_Lagrange_3D_Scal%shortID)
!            Call Element_Q_Lagrange_3D_Scal_Init(dElem,dCoord,2,QuadratureOrder)
         Case Default
            Print*,__FUNCT__,': Element type not implemented yet',elemType%name,elemType%shortID
      End Select
   End Subroutine Element3D_Scal_Init                                
   
#undef __FUNCT__
#define __FUNCT__ "Element3D_Vect_Init"
   Subroutine Element3D_Vect_Init(dElem,dCoord,QuadratureOrder,elemType)
      Type(Element3D_Vect)                   :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord
      PetscInt,Intent(IN)                    :: QuadratureOrder
      Type(Element_Type),intent(IN)          :: elemType
      
      Select Case (elemType%shortID)
         Case (MEF90_P1_Lagrange_3D_Vect%shortID)
            Call Element_P_Lagrange_3D_Vect_Init(dElem,dCoord,1,QuadratureOrder)
         Case (MEF90_P2_Lagrange_3D_Vect%shortID)
            Call Element_P_Lagrange_3D_Vect_Init(dElem,dCoord,2,QuadratureOrder)
         Case (MEF90_P1_Lagrange_3DBoundary_Vect%shortID)
            Call Element_P_Lagrange_3DBoundary_Vect_Init(dElem,dCoord,1,QuadratureOrder)
         Case (MEF90_P2_Lagrange_3DBoundary_Vect%shortID)
            Call Element_P_Lagrange_3DBoundary_Vect_Init(dElem,dCoord,2,QuadratureOrder)
!         Case (MEF90_Q1_Lagrange_3D_Vect)
!            Call Element_Q_Lagrange_3D_Vect_Init(dElem,dCoord,1,QuadratureOrder)
!         Case (MEF90_Q2_Lagrange_3D_Vect)
!            Call Element_Q_Lagrange_3D_Vect_Init(dElem,dCoord,2,QuadratureOrder)
         Case Default
            Print*,__FUNCT__,': Element type not implemented yet',elemType%name,elemType%shortID
      End Select
   End Subroutine Element3D_Vect_Init                                

#undef __FUNCT__
#define __FUNCT__ "Element3D_Elast_Init"
   Subroutine Element3D_Elast_Init(dElem,dCoord,QuadratureOrder,elemType)
      Type(Element3D_Elast)                   :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord
      PetscInt,Intent(IN)                    :: QuadratureOrder
      Type(Element_Type),intent(IN)          :: elemType
      
      Select Case (elemType%shortID)
         Case (MEF90_P1_Lagrange_3D_Elast%shortID)
            Call Element_P_Lagrange_3D_Elast_Init(dElem,dCoord,1,QuadratureOrder)
         Case (MEF90_P2_Lagrange_3D_Elast%shortID)
            Call Element_P_Lagrange_3D_Elast_Init(dElem,dCoord,2,QuadratureOrder)
         Case (MEF90_P1_Lagrange_3DBoundary_Elast%shortID)
            Call Element_P_Lagrange_3DBoundary_Elast_Init(dElem,dCoord,1,QuadratureOrder)
         Case (MEF90_P2_Lagrange_3DBoundary_Elast%shortID)
            Call Element_P_Lagrange_3DBoundary_Elast_Init(dElem,dCoord,2,QuadratureOrder)
!         Case (MEF90_Q1_Lagrange_3D_Elast%shortID)
!            Call Element_Q_Lagrange_3D_Elast_Init(dElem,dCoord,1,QuadratureOrder)
!         Case (MEF90_Q2_Lagrange_3D_Elast%shortID)
!            Call Element_Q_Lagrange_3D_Elast_Init(dElem,dCoord,2,QuadratureOrder)
         Case Default
            Print*,__FUNCT__,': Element type not implemented yet',elemType%name,elemType%shortID
      End Select
   End Subroutine Element3D_Elast_Init                                



#undef __FUNCT__
#define __FUNCT__ "Element_P_Lagrange_2D_Scal_Init"
   Subroutine Element_P_Lagrange_2D_Scal_Init(dElem,dCoord,dPolynomialOrder,dQuadratureOrder)
      ! Compute the quadrature weights and the value of the basis functions and their gradient 
      ! at the quadrature points.
      Type(Element2D_Scal)                   :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord      ! coord(i,j)=ith coord of jth vertice
      PetscInt                               :: dPolynomialOrder,dQuadratureOrder
      
      PetscInt                               :: Nb_Gauss
      PetscInt                               :: Num_Dof
      PetscInt                               :: iDoF,iG
      Type(Mat2D)                            :: Bt          ! The transposed of transformation matrix
      PetscReal                              :: DetBinv     ! The determinant of B^{-1}

      PetscReal,Dimension(:,:),Pointer       :: PhiHat      ! PhiHat(i,k) The value of the ith basis function at the kth integration point
      Type(Vect2D),Dimension(:,:),Pointer    :: GradPhiHat
      
      
      Type(Vect2D),Dimension(:),Pointer     :: Xi ! The quadrature points coordinates in the reference element
      
      !!! The transformation matrix and the determinant of its inverse
      Bt%XX = dCoord(1,2) - dCoord(1,1) 
      Bt%XY = dCoord(2,2) - dCoord(2,1)
      Bt%YX = dCoord(1,3) - dCoord(1,1)
      Bt%YY = dCoord(2,3) - dCoord(2,1)
      
      DetBinv = Det(Bt)
      Bt = Invert(Bt)

      Select Case (dQuadratureOrder)
      Case(1)
         Nb_Gauss = 1
         Allocate(Xi(Nb_Gauss))
         Allocate(dElem%Gauss_C(Nb_Gauss))
         Xi(1)%X = 1.0_Kr / 3.0_Kr
         Xi(1)%Y = 1.0_Kr / 3.0_Kr
         dElem%Gauss_C = detBinv / 2.0_Kr

      Case(2)
         Nb_Gauss = 3
         Allocate(Xi(Nb_Gauss))
         Allocate(dElem%Gauss_C(Nb_Gauss))
         Xi(1) = (/ 1.0_Kr / 6.0_Kr,1.0_Kr / 6.0_Kr /)
         Xi(2) = (/ 2.0_Kr / 3.0_Kr,1.0_Kr / 6.0_Kr /)
         Xi(3) = (/ 1.0_Kr / 6.0_Kr,2.0_Kr / 3.0_Kr /)
         dElem%Gauss_C = detBinv / 6.0_Kr

      Case(3)
         Nb_Gauss = 4
         Allocate(Xi(Nb_Gauss))
         Allocate(dElem%Gauss_C(Nb_Gauss))
         dElem%Gauss_C    =  detBinv * 25.0_Kr / 96.0_Kr
         dElem%Gauss_C(1) = -detBinv * 9.0_Kr / 32.0_Kr
         Xi(1) = (/ 1.0_Kr / 3.0_Kr,1.0_Kr / 3.0_Kr /)
         Xi(2) = (/ 3.0_Kr / 5.0_Kr,1.0_Kr / 5.0_Kr /)
         Xi(3) = (/ 1.0_Kr / 5.0_Kr,3.0_Kr / 5.0_Kr /)
         Xi(4) = (/ 1.0_Kr / 5.0_Kr,1.0_Kr / 5.0_Kr /)
      Case(4)
         Nb_Gauss = 7
         Allocate(Xi(Nb_Gauss))
         Allocate(dElem%Gauss_C(Nb_Gauss))
         dElem%Gauss_C(1) = detBinv / 40.0_Kr
         dElem%Gauss_C(2) = detBinv / 15.0_Kr
         dElem%Gauss_C(3) = detBinv / 40.0_Kr
         dElem%Gauss_C(4) = detBinv / 15.0_Kr
         dElem%Gauss_C(5) = detBinv / 40.0_Kr
         dElem%Gauss_C(6) = detBinv / 15.0_Kr
         dElem%Gauss_C(7) = detBinv * 9.0_Kr / 40.0_Kr
         Xi(1) = (/ 0.0_Kr         ,0.0_Kr          /)
         Xi(2) = (/ 1.0_Kr / 2.0_Kr,0.0_Kr          /)
         Xi(3) = (/ 1.0_Kr         ,0.0_Kr          /)
         Xi(4) = (/ 1.0_Kr / 2.0_Kr,1.0_Kr / 2.0_Kr /)
         Xi(5) = (/ 0.0_Kr         ,1.0_Kr          /)
         Xi(6) = (/ 0.0_Kr         ,1.0_Kr / 2.0_Kr /)
         Xi(7) = (/ 1.0_Kr / 3.0_Kr,1.0_Kr / 3.0_Kr /)
      Case(6)
         Nb_Gauss = 9
         Allocate(Xi(Nb_Gauss))
         Allocate(dElem%Gauss_C(Nb_Gauss))
         Xi(1) = (/ 0.124949503233232_Kr, 0.437525248383384_Kr /)
         Xi(2) = (/ 0.437525248383384_Kr, 0.124949503233232_Kr /)
         Xi(3) = (/ 0.437525248383384_Kr, 0.437525248383384_Kr /)
         Xi(4) = (/ 0.797112651860071_Kr, 0.165409927389841_Kr /)
         Xi(5) = (/ 0.797112651860071_Kr, 0.037477420750088_Kr /)
         Xi(6) = (/ 0.165409927389841_Kr, 0.797112651860071_Kr /)
         Xi(7) = (/ 0.165409927389841_Kr, 0.037477420750088_Kr /)
         Xi(8) = (/ 0.037477420750088_Kr, 0.797112651860071_Kr /)
         Xi(9) = (/ 0.037477420750088_Kr, 0.165409927389841_Kr /)
         dElem%Gauss_C(1:3) = 0.205950504760887_Kr / 2.0_Kr * detBinv
         dElem%Gauss_C(4:9) = 0.063691414286223_Kr / 2.0_Kr * detBinv
      Case Default
         Print*,__FUNCT__,': Unimplemented quadrature order',dQuadratureOrder
         STOP
      End Select
      
      Select Case (dPolynomialOrder)
      Case(1)
         Num_DoF = 3
         Allocate(PhiHat(Num_DoF,Nb_Gauss))
         Allocate(GradPhiHat(Num_DoF,Nb_Gauss))
         PhiHat(1,:) = 1.0_Kr - Xi%X - Xi%Y
         PhiHat(2,:) = Xi(:)%X
         PhiHat(3,:) = Xi(:)%Y
         
         GradPhiHat(1,:)%X = -1.0_Kr;GradPhiHat(1,:)%Y = -1.0_Kr 
         GradPhiHat(2,:)%X =  1.0_Kr;GradPhiHat(2,:)%Y =  0.0_Kr 
         GradPhiHat(3,:)%X =  0.0_Kr;GradPhiHat(3,:)%Y =  1.0_Kr
          
      Case(2)
         Num_DoF = 6
         Allocate(PhiHat(Num_DoF,Nb_Gauss))
         Allocate(GradPhiHat(Num_DoF,Nb_Gauss))
         PhiHat(1,:) = (1.0_Kr - Xi%X - Xi%Y) * (1.0_Kr - 2.0_Kr * Xi%X - 2.0_Kr * Xi%Y)      
         PhiHat(2,:) = Xi%X * (2.0_Kr * Xi%X - 1.0_Kr)
         PhiHat(3,:) = Xi%Y * (2.0_Kr * Xi%Y - 1.0_Kr)
         PhiHat(4,:) = 4.0_Kr * Xi%X * (1.0_Kr - Xi%X - Xi%Y)
         PhiHat(5,:) = 4.0_Kr * Xi%X * Xi%Y
         PhiHat(6,:) = 4.0_Kr * Xi%X * (1.0_Kr - Xi%X - Xi%Y)
         
         GradPhiHat(1,:)%X = 4.0_Kr * Xi%X + 4.0_Kr * Xi%y -3.0_Kr;     GradPhiHat(1,:)%Y = 4.0_Kr * Xi%X + 4.0_Kr * Xi%y -3.0_Kr
         GradPhiHat(2,:)%X = 4.0_Kr * Xi%X - 1.0_Kr;                    GradPhiHat(2,:)%Y = 0.0_Kr
         GradPhiHat(3,:)%X = 0.0_Kr;                                    GradPhiHat(3,:)%Y = 4.0_Kr * Xi%Y - 1.0_Kr
         GradPhiHat(4,:)%X = 4.0_Kr * (1.0_Kr - 2.0_Kr * Xi%X - Xi%Y);  GradPhiHat(4,:)%Y = -4.0_Kr * Xi%X
         GradPhiHat(5,:)%X = 4.0_Kr * Xi%Y;                             GradPhiHat(5,:)%Y = 4.0_Kr * Xi%X
         GradPhiHat(6,:)%X = -4.0_Kr * Xi%Y;                            GradPhiHat(6,:)%Y = 4.0_Kr * (1.0_Kr - Xi%X - 2.0_Kr * Xi%Y);
      Case Default
         Print*,__FUNCT__,': Unimplemented PolynomialOrder',dPolynomialOrder
      End Select
      
      Allocate (dElem%BF(Num_DoF,Nb_Gauss)) 
      Allocate (dElem%Grad_BF(Num_DoF,Nb_Gauss))
      dElem%BF = PhiHat
      Do iDoF = 1,Num_DoF
         Do iG = 1,Nb_Gauss
            dElem%Grad_BF(iDoF,iG) = Bt * GradPhiHat(iDoF,iG) 
         End Do
      End Do
     
      DeAllocate(Xi)
      DeAllocate(PhiHat)
      DeAllocate(GradPhiHat)
   End Subroutine Element_P_Lagrange_2D_Scal_Init
   
#undef __FUNCT__
#define __FUNCT__ "Element_P_Lagrange_2DBoundary_Scal_Init"
   Subroutine Element_P_Lagrange_2DBoundary_Scal_Init(dElem,dCoord,dPolynomialOrder,dQuadratureOrder)
      Type(Element2D_Scal),intent(INOUT)     :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord      ! coord(i,j)=ith coord of jth vertice
      PetscInt,Intent(IN)                    :: dPolynomialOrder,dQuadratureOrder
      
      Type(Element2D_Scal)                   :: tmpElem
      PetscReal,Dimension(:,:),Pointer       :: tmpCoord
      PetscInt                               :: i,j,iDoF,iG,Num_Gauss,Num_DoF
      Type(Vect2D)                           :: NormalVector
      PetscReal                              :: InnerBF
      
      !!! Create a bogus tri element with unit height by adding a 3rd vertex
      Allocate(tmpCoord(2,3))
      Do i = 1,2
         Do j = 1,2
            tmpCoord(i,j) = dCoord(i,j)
         End Do
      End Do
      NormalVector%X = tmpCoord(2,2) - tmpCoord(2,1)
      NormalVector%Y = tmpCoord(1,1) - tmpCoord(1,2)
      NormalVector = NormalVector / Norm(NormalVector)
      tmpCoord(1,3) = (dCoord(1,1) + dCoord(1,2))*.5_Kr - NormalVector%X 
      tmpCoord(2,3) = (dCoord(2,1) + dCoord(2,2))*.5_Kr - NormalVector%Y
      tmpCoord(1,3) = dCoord(1,1) - NormalVector%X 
      tmpCoord(2,3) = dCoord(2,1) - NormalVector%Y
      
      Call Element_P_Lagrange_2D_Scal_Init(tmpElem,tmpCoord,dPolynomialOrder,dQuadratureOrder)
      Select Case (dPolynomialOrder)
         Case (1)
            Num_DoF  = 2
            Num_Gauss = size(tmpElem%BF,2)
            Allocate(dElem%Gauss_C(Num_Gauss))
            Allocate(dElem%BF(Num_DoF,Num_Gauss))
            dElem%Gauss_C = tmpElem%Gauss_C * 2.0_Kr
            Do iDoF = 1,Num_doF
               Do iG = 1,Num_Gauss
                  dElem%BF(iDoF,iG) = tmpElem%BF(iDoF,iG) + tmpElem%BF(Num_DoF+1,iG) * .5_Kr
               End Do
            End Do
         Case (2)
            !!! dof corrrespondance between a BEAM3 and  a TRI6 element (exodus.pdf figure 4 p.  17)
            !!! BEAM3   TRI6
            !!! 1       1
            !!! 2       2
            !!! 3       6
            !!! Inner dof: 3,4,5
            Num_DoF  = 3
            Num_Gauss = size(tmpElem%BF,2)
            Allocate(dElem%Gauss_C(Num_Gauss))
            Allocate(dElem%BF(Num_DoF,Num_Gauss))
            dElem%Gauss_C = tmpElem%Gauss_C * 2.0_Kr
            Do iG = 1,Num_Gauss
               InnerBF = (tmpElem%BF(3,iG) + tmpElem%BF(4,iG) + tmpElem%BF(5,iG)) / 3.0_Kr
               dElem%BF(1,iG) = tmpElem%BF(1,iG) + InnerBF
               dElem%BF(2,iG) = tmpElem%BF(2,iG) + InnerBF
               dElem%BF(3,iG) = tmpElem%BF(6,iG) + InnerBF
            End Do
         Case Default
            Print*,'[ERROR]: Polynomial order ',dPolynomialOrder,' not implemented in ',__FUNCT__
      End Select
      Call Element2D_Scal_Destroy(tmpElem)
      deAllocate(tmpCoord)
      Allocate(delem%Grad_BF(0,0))
   End Subroutine Element_P_Lagrange_2DBoundary_Scal_Init                          
      
#undef __FUNCT__
#define __FUNCT__ "Element_P_Lagrange_2D_Vect_Init"
   Subroutine Element_P_Lagrange_2D_Vect_Init(dElem,dCoord,dPolynomialOrder,dQuadratureOrder)
      Type(Element2D_Vect)                   :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord      ! coord(i,j)=ith coord of jth vertice
      PetscInt                               :: dPolynomialOrder,dQuadratureOrder
   
      Type(Element2D_Scal)                   :: Elem_Scal
      PetscInt                               :: dim = 2 
      PetscInt                               :: Num_DoF,Nb_Gauss,i
      
      
      Call Element_P_Lagrange_2D_Scal_Init(Elem_Scal,dCoord,dPolynomialOrder,dQuadratureOrder)
      Num_DoF  = Size(Elem_Scal%BF,1) 
      Nb_Gauss = Size(Elem_Scal%BF,2)
      Allocate(dElem%Gauss_C(Nb_Gauss))
      Allocate(dElem%BF(Num_DoF * dim,Nb_Gauss))
      Allocate(dElem%Der_BF(Num_DoF * dim,Nb_Gauss))
         
      dElem%Gauss_C = Elem_Scal%Gauss_C
      dElem%BF(:,:)%X = 0.0_Kr
      dElem%BF(:,:)%Y = 0.0_Kr
      dElem%Der_BF(:,:)%XX = 0.0_Kr
      dElem%Der_BF(:,:)%XY = 0.0_Kr
      dElem%Der_BF(:,:)%YX = 0.0_Kr
      dElem%Der_BF(:,:)%YY = 0.0_Kr
      
      Do i = 0,Num_DoF-1
         dElem%BF(i*dim+1,:)%X = Elem_Scal%BF(i+1,:)
         dElem%BF(i*dim+2,:)%Y = Elem_Scal%BF(i+1,:)
         dElem%Der_BF(i*dim+1,:)%XX = Elem_Scal%Grad_BF(i+1,:)%X
         dElem%Der_BF(i*dim+1,:)%XY = Elem_Scal%Grad_BF(i+1,:)%Y
         dElem%Der_BF(i*dim+2,:)%YX = Elem_Scal%Grad_BF(i+1,:)%X
         dElem%Der_BF(i*dim+2,:)%YY = Elem_Scal%Grad_BF(i+1,:)%Y
      End Do
      Call ElementDestroy(Elem_Scal)
   End Subroutine Element_P_Lagrange_2D_Vect_Init

#undef __FUNCT__
#define __FUNCT__ "Element_P_Lagrange_2DBoundary_Vect_Init"
   Subroutine Element_P_Lagrange_2DBoundary_Vect_Init(dElem,dCoord,dPolynomialOrder,dQuadratureOrder)
      Type(Element2D_Vect)                   :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord      ! coord(i,j)=ith coord of jth vertice
      PetscInt                               :: dPolynomialOrder,dQuadratureOrder
   
      Type(Element2D_Scal)                   :: Elem_Scal
      PetscInt                               :: dim = 2 
      PetscInt                               :: Num_DoF,Nb_Gauss,i
      
      
      Call Element_P_Lagrange_2DBoundary_Scal_Init(Elem_Scal,dCoord,dPolynomialOrder,dQuadratureOrder)
      Num_DoF  = Size(Elem_Scal%BF,1) 
      Nb_Gauss = Size(Elem_Scal%BF,2)
      Allocate(dElem%Gauss_C(Nb_Gauss))
      Allocate(dElem%BF(Num_DoF * dim,Nb_Gauss))
      Allocate(dElem%Der_BF(0,Nb_Gauss))
         
      dElem%Gauss_C = Elem_Scal%Gauss_C
      dElem%BF(:,:)%X = 0.0_Kr
      dElem%BF(:,:)%Y = 0.0_Kr
      
      Do i = 0,Num_DoF-1
         dElem%BF(i*dim+1,:)%X = Elem_Scal%BF(i+1,:)
         dElem%BF(i*dim+2,:)%Y = Elem_Scal%BF(i+1,:)
      End Do
      Call ElementDestroy(Elem_Scal)
   End Subroutine Element_P_Lagrange_2DBoundary_Vect_Init

#undef __FUNCT__
#define __FUNCT__ "Element_P_Lagrange_2D_Elast_Init"
   Subroutine Element_P_Lagrange_2D_Elast_Init(dElem,dCoord,dPolynomialOrder,dQuadratureOrder)
      Type(Element2D_Elast)                  :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord      ! coord(i,j)=ith coord of jth vertice
      PetscInt                               :: dPolynomialOrder,dQuadratureOrder
   
      Type(Element2D_Scal)                   :: Elem_Scal
      PetscInt                               :: dim = 2 
      PetscInt                               :: Num_DoF,Nb_Gauss,i
      
      
      Call Element_P_Lagrange_2D_Scal_Init(Elem_Scal,dCoord,dPolynomialOrder,dQuadratureOrder)
      Num_DoF   = Size(Elem_Scal%BF,1) 
      Nb_Gauss = Size(Elem_Scal%BF,2)
      Allocate(dElem%Gauss_C(Nb_Gauss))
      Allocate(dElem%BF(Num_DoF * dim,Nb_Gauss))
      Allocate(dElem%GradS_BF(Num_DoF * dim,Nb_Gauss))
         
      dElem%Gauss_C = Elem_Scal%Gauss_C
      dElem%BF(:,:)%X = 0.0_Kr
      dElem%BF(:,:)%Y = 0.0_Kr
      dElem%GradS_BF(:,:)%XX = 0.0_Kr
      dElem%GradS_BF(:,:)%XY = 0.0_Kr
      dElem%GradS_BF(:,:)%YY = 0.0_Kr
      Do i = 0,Num_DoF-1
         dElem%BF(i*dim+1,:)%X = Elem_Scal%BF(i+1,:)
         dElem%BF(i*dim+2,:)%Y = Elem_Scal%BF(i+1,:)
         dElem%GradS_BF(i*dim+1,:)%XX = Elem_Scal%Grad_BF(i+1,:)%X
         dElem%GradS_BF(i*dim+1,:)%XY = Elem_Scal%Grad_BF(i+1,:)%Y / 2.0_Kr
         dElem%GradS_BF(i*dim+2,:)%XY = Elem_Scal%Grad_BF(i+1,:)%X / 2.0_Kr
         dElem%GradS_BF(i*dim+2,:)%YY = Elem_Scal%Grad_BF(i+1,:)%Y
      End Do
      Call ElementDestroy(Elem_Scal)
   End Subroutine Element_P_Lagrange_2D_Elast_Init

#undef __FUNCT__
#define __FUNCT__ "Element_P_Lagrange_2DBoundary_Elast_Init"
   Subroutine Element_P_Lagrange_2DBoundary_Elast_Init(dElem,dCoord,dPolynomialOrder,dQuadratureOrder)
      Type(Element2D_Elast)                  :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord      ! coord(i,j)=ith coord of jth vertice
      PetscInt                               :: dPolynomialOrder,dQuadratureOrder
   
      Type(Element2D_Scal)                   :: Elem_Scal
      PetscInt                               :: dim = 2 
      PetscInt                               :: Num_DoF,Nb_Gauss,i
      
      
      Call Element_P_Lagrange_2DBoundary_Scal_Init(Elem_Scal,dCoord,dPolynomialOrder,dQuadratureOrder)
      Num_DoF   = Size(Elem_Scal%BF,1) 
      Nb_Gauss = Size(Elem_Scal%BF,2)
      Allocate(dElem%Gauss_C(Nb_Gauss))
      Allocate(dElem%BF(Num_DoF * dim,Nb_Gauss))
      Allocate(dElem%GradS_BF(0,Nb_Gauss))
         
      dElem%Gauss_C = Elem_Scal%Gauss_C
      dElem%BF(:,:)%X = 0.0_Kr
      dElem%BF(:,:)%Y = 0.0_Kr
      Do i = 0,Num_DoF-1
         dElem%BF(i*dim+1,:)%X = Elem_Scal%BF(i+1,:)
         dElem%BF(i*dim+2,:)%Y = Elem_Scal%BF(i+1,:)
      End Do
      Call ElementDestroy(Elem_Scal)
   End Subroutine Element_P_Lagrange_2DBoundary_Elast_Init

#undef __FUNCT__
#define __FUNCT__ "Element_P_Lagrange_3D_Scal_Init"
   Subroutine Element_P_Lagrange_3D_Scal_Init(dElem,dCoord,dPolynomialOrder,dQuadratureOrder)
      ! Compute the quadrature weights and the value of the basis functions and their gradient 
      ! at the quadrature points.
      Type(Element3D_Scal)                   :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord      ! coord(i,j)=ith coord of jth vertice
      PetscInt                               :: dPolynomialOrder,dQuadratureOrder
      
      PetscInt                               :: Nb_Gauss
      PetscInt                               :: Num_Dof
      PetscInt                               :: iDoF,iG
      Type(Mat3D)                            :: Bt          ! The transposed of transformation matrix
      PetscReal                              :: DetBinv     ! The determinant of B^{-1}

      PetscReal,Dimension(:,:),Pointer       :: PhiHat      ! PhiHat(i,k) The value of the ith basis function at the kth integration point
      Type(Vect3D),Dimension(:,:),Pointer    :: GradPhiHat
      
      
      Type(Vect3D),Dimension(:),Pointer      :: Xi          ! The quadrature points coordinates in the reference element
      PetscReal                              :: a,b         ! Location of integration points in Aiken p. 272 table 10.4
      
      !!! The transformation matrix and the determinant of its inverse
      Bt%XX = dCoord(1,2) - dCoord(1,1) 
      Bt%XY = dCoord(2,2) - dCoord(2,1)
      Bt%XZ = dCoord(3,2) - dCoord(3,1)
      
      Bt%YX = dCoord(1,3) - dCoord(1,1)
      Bt%YY = dCoord(2,3) - dCoord(2,1)
      Bt%YZ = dCoord(3,3) - dCoord(3,1)
      
      Bt%ZX = dCoord(1,4) - dCoord(1,1)
      Bt%ZY = dCoord(2,4) - dCoord(2,1)
      Bt%ZZ = dCoord(3,4) - dCoord(3,1)
      
      DetBinv = Det(Bt)
      Bt = Invert(Bt)
      
      Select Case (dQuadratureOrder)
      Case (1)
         Nb_Gauss = 1
         Allocate(Xi(Nb_Gauss))
         Allocate(dElem%Gauss_C(Nb_Gauss))
         Xi(1) = (/ .25_Kr,.25_Kr,.25_Kr /)
         dElem%Gauss_C(1) = 1.0_Kr / 6.0_Kr * detBinv
         
      Case(2)
         a = (5.0_Kr + 3.0_Kr * sqrt(5.0_Kr)) / 20.0_Kr
         b = (5.0_Kr - sqrt(5.0_Kr)) / 20.0_Kr
         Nb_Gauss = 4
         Allocate(Xi(Nb_Gauss))
         Allocate(dElem%Gauss_C(Nb_Gauss))
         Xi(1) = (/ a,b,b /)
         Xi(2) = (/ b,a,b /)
         Xi(3) = (/ b,b,a /)
         Xi(4) = (/ b,b,b /)
         dElem%Gauss_C(1:4) = 1.0_Kr / 24.0_Kr * detBinv
            
      Case(3)
         Nb_Gauss = 5
         Allocate(Xi(Nb_Gauss))
         Allocate(dElem%Gauss_C(Nb_Gauss))
         Xi(1) = (/ .25_Kr,.25_Kr,.25_Kr /)
         Xi(2) = (/ .5_Kr,1._Kr / 6._Kr,1._Kr / 6._Kr /)
         Xi(3) = (/ 1._Kr / 6._Kr,.5_Kr,1._Kr / 6._Kr /)
         Xi(4) = (/ 1._Kr / 6._Kr,1._Kr / 6._Kr,.5_Kr /)
         Xi(5) = (/ 1._Kr / 6._Kr,1._Kr / 6._Kr,1._Kr / 6._Kr /)
         dElem%Gauss_C(1) = -2.0_Kr / 15.0_Kr * detBinv
         dElem%Gauss_C(2:5) = 3.0_Kr / 40.0_Kr * detBinv
      Case(4)
         Nb_Gauss = 11
         Allocate(Xi(Nb_Gauss))
         Allocate(dElem%Gauss_C(Nb_Gauss))
         Xi(1)  = (/ 0.2500000000000000_Kr, 0.2500000000000000_Kr, 0.2500000000000000_Kr /)
         Xi(2)  = (/ 0.7857142857142857_Kr, 0.0714285714285714_Kr, 0.0714285714285714_Kr /)
         Xi(3)  = (/ 0.0714285714285714_Kr, 0.0714285714285714_Kr, 0.0714285714285714_Kr /)
         Xi(4)  = (/ 0.0714285714285714_Kr, 0.0714285714285714_Kr, 0.7857142857142857_Kr /)
         Xi(5)  = (/ 0.0714285714285714_Kr, 0.7857142857142857_Kr, 0.0714285714285714_Kr /)
         Xi(6)  = (/ 0.1005964238332008_Kr, 0.3994035761667992_Kr, 0.3994035761667992_Kr /)
         Xi(7)  = (/ 0.3994035761667992_Kr, 0.1005964238332008_Kr, 0.3994035761667992_Kr /)
         Xi(8)  = (/ 0.3994035761667992_Kr, 0.3994035761667992_Kr, 0.1005964238332008_Kr /)
         Xi(9)  = (/ 0.3994035761667992_Kr, 0.1005964238332008_Kr, 0.1005964238332008_Kr /)
         Xi(10) = (/ 0.1005964238332008_Kr, 0.3994035761667992_Kr, 0.1005964238332008_Kr /)
         Xi(11) = (/ 0.1005964238332008_Kr, 0.1005964238332008_Kr, 0.3994035761667992_Kr /)
         dElem%Gauss_C(1)    = -0.0789333333333333_Kr / 6.0_Kr * detBinv
         dElem%Gauss_C(2:5)  =  0.0457333333333333_Kr / 6.0_Kr * detBinv
         dElem%Gauss_C(6:11) =  0.1493333333333333_Kr / 6.0_Kr * detBinv
      Case(6)
         Nb_Gauss = 24
         Allocate(Xi(Nb_Gauss))
         Allocate(dElem%Gauss_C(Nb_Gauss))
         Xi(1)  = (/ 0.3561913862225449_Kr, 0.2146028712591517_Kr, 0.2146028712591517_Kr /)
         Xi(2)  = (/ 0.2146028712591517_Kr, 0.2146028712591517_Kr, 0.2146028712591517_Kr /)
         Xi(3)  = (/ 0.2146028712591517_Kr, 0.2146028712591517_Kr, 0.3561913862225449_Kr /)
         Xi(4)  = (/ 0.2146028712591517_Kr, 0.3561913862225449_Kr, 0.2146028712591517_Kr /)
         Xi(5)  = (/ 0.8779781243961660_Kr, 0.0406739585346113_Kr, 0.0406739585346113_Kr /)
         Xi(6)  = (/ 0.0406739585346113_Kr, 0.0406739585346113_Kr, 0.0406739585346113_Kr /)
         Xi(7)  = (/ 0.0406739585346113_Kr, 0.0406739585346113_Kr, 0.8779781243961660_Kr /)
         Xi(8)  = (/ 0.0406739585346113_Kr, 0.8779781243961660_Kr, 0.0406739585346113_Kr /)
         Xi(9)  = (/ 0.0329863295731731_Kr, 0.3223378901422757_Kr, 0.3223378901422757_Kr /)
         Xi(10) = (/ 0.3223378901422757_Kr, 0.3223378901422757_Kr, 0.3223378901422757_Kr /)
         Xi(11) = (/ 0.3223378901422757_Kr, 0.3223378901422757_Kr, 0.0329863295731731_Kr /)
         Xi(12) = (/ 0.3223378901422757_Kr, 0.0329863295731731_Kr, 0.3223378901422757_Kr /)
         Xi(13) = (/ 0.2696723314583159_Kr, 0.0636610018750175_Kr, 0.0636610018750175_Kr /)
         Xi(14) = (/ 0.0636610018750175_Kr, 0.2696723314583159_Kr, 0.0636610018750175_Kr /)
         Xi(15) = (/ 0.0636610018750175_Kr, 0.0636610018750175_Kr, 0.2696723314583159_Kr /)
         Xi(16) = (/ 0.6030056647916491_Kr, 0.0636610018750175_Kr, 0.0636610018750175_Kr /)
         Xi(17) = (/ 0.0636610018750175_Kr, 0.6030056647916491_Kr, 0.0636610018750175_Kr /)
         Xi(18) = (/ 0.0636610018750175_Kr, 0.0636610018750175_Kr, 0.6030056647916491_Kr /)
         Xi(19) = (/ 0.0636610018750175_Kr, 0.2696723314583159_Kr, 0.6030056647916491_Kr /)
         Xi(20) = (/ 0.2696723314583159_Kr, 0.6030056647916491_Kr, 0.0636610018750175_Kr /)
         Xi(11) = (/ 0.6030056647916491_Kr, 0.0636610018750175_Kr, 0.2696723314583159_Kr /)
         Xi(22) = (/ 0.0636610018750175_Kr, 0.6030056647916491_Kr, 0.2696723314583159_Kr /)
         Xi(23) = (/ 0.2696723314583159_Kr, 0.0636610018750175_Kr, 0.6030056647916491_Kr /)
         Xi(24) = (/ 0.6030056647916491_Kr, 0.2696723314583159_Kr, 0.0636610018750175_Kr /)
         dElem%Gauss_C(1:4)   = 0.0399227502581679_Kr / 6.0_Kr * detBinv
         dElem%Gauss_C(5:8)   = 0.0100772110553207_Kr / 6.0_Kr * detBinv
         dElem%Gauss_C(9:12)  = 0.0553571815436544_Kr / 6.0_Kr * detBinv
         dElem%Gauss_C(13:24) = 0.0482142857142857_Kr / 6.0_Kr * detBinv
      Case Default
         Print*,__FUNCT__,': Unimplemented quadrature order',dQuadratureOrder
         STOP
      End Select
      
      Select Case (dPolynomialOrder)
      Case(1)
         Num_DoF = 4
         Allocate(PhiHat(Num_DoF,Nb_Gauss))
         Allocate(GradPhiHat(Num_DoF,Nb_Gauss))
         PhiHat(1,:) = 1.0_Kr - Xi%X - Xi%Y - Xi%Z
         PhiHat(2,:) = Xi(:)%X
         PhiHat(3,:) = Xi(:)%Y
         PhiHat(4,:) = Xi(:)%Z
         
         GradPhiHat(1,:)%X = -1.0_Kr;GradPhiHat(1,:)%Y = -1.0_Kr;GradPhiHat(1,:)%Z = -1.0_Kr;
         GradPhiHat(2,:)%X =  1.0_Kr;GradPhiHat(2,:)%Y =  0.0_Kr;GradPhiHat(2,:)%Z =  0.0_Kr;
         GradPhiHat(3,:)%X =  0.0_Kr;GradPhiHat(3,:)%Y =  1.0_Kr;GradPhiHat(3,:)%Z =  0.0_Kr;
         GradPhiHat(4,:)%X =  0.0_Kr;GradPhiHat(4,:)%Y =  0.0_Kr;GradPhiHat(4,:)%Z =  1.0_Kr;
      Case(2)
         Num_Dof = 10
         Allocate(PhiHat(Num_DoF,Nb_Gauss))
         Allocate(GradPhiHat(Num_DoF,Nb_Gauss))
         PhiHat(1,:)  = (1.0_Kr - Xi%X - Xi%Y - Xi%Z) * (1.0_Kr - 2.0_Kr*Xi%X - 2.0_Kr*Xi%Y - 2.0_Kr*Xi%Z)
         PhiHat(2,:)  = Xi%X * (2.0_Kr *  Xi%X - 1.0_Kr)
         PhiHat(3,:)  = Xi%Y * (2.0_Kr * Xi%Y - 1.0_Kr)
         PhiHat(4,:)  = Xi%Z * (2.0_Kr * Xi%Z - 1.0_Kr)
         PhiHat(5,:)  = 4.0_Kr * (1.0_Kr - Xi%X - Xi%Y - Xi%Z) * Xi%X
         PhiHat(6,:)  = 4.0_Kr *  Xi%X * Xi%Y
         PhiHat(7,:)  = 4.0_Kr * (1.0_Kr - Xi%X - Xi%Y - Xi%Z) * Xi%Y
!         PhiHat(8,:)  = 4.0_Kr *  Xi%X * Xi%Z
!         PhiHat(9,:)  = 4.0_Kr * Xi%Y * Xi%Z
!         PhiHat(10,:) = 4.0_Kr * (1.0_Kr - Xi%X - Xi%Y - Xi%Z) * Xi%Z
         PhiHat(8,:)  = 4.0_Kr * (1.0_Kr - Xi%X - Xi%Y - Xi%Z) * Xi%Z
         PhiHat(9,:)  = 4.0_Kr * Xi%X * Xi%Z
         PhiHat(10,:) = 4.0_Kr * Xi%Y * Xi%Z
         
         GradPhiHat(1,:)%X = 4.0_Kr * Xi%X + 4.0_Kr * Xi%Y + 4.0_Kr * Xi%Z - 3.0_Kr
         GradPhiHat(1,:)%Y = 4.0_Kr * Xi%X + 4.0_Kr * Xi%Y + 4.0_Kr * Xi%Z - 3.0_Kr
         GradPhiHat(1,:)%Z = 4.0_Kr * Xi%X + 4.0_Kr * Xi%Y + 4.0_Kr * Xi%Z - 3.0_Kr

         GradPhiHat(2,:)%X = 4.0_Kr * Xi%X - 1.0_Kr
         GradPhiHat(2,:)%Y = 0.0_Kr
         GradPhiHat(2,:)%Z = 0.0_Kr
         
         GradPhiHat(3,:)%X = 0.0_Kr
         GradPhiHat(3,:)%Y = 4.0_Kr * Xi%Y - 1.0_Kr
         GradPhiHat(3,:)%Z = 0.0_Kr
         
         GradPhiHat(4,:)%X = 0.0_Kr
         GradPhiHat(4,:)%Y = 0.0_Kr
         GradPhiHat(4,:)%Z = 4.0_Kr * Xi%Z - 1.0_Kr
         
         GradPhiHat(5,:)%X = 4.0_Kr * (1.0_Kr - 2.0_Kr * Xi%X - Xi%Y - Xi%Z)
         GradPhiHat(5,:)%Y = -4.0_Kr * Xi%X
         GradPhiHat(5,:)%Z = -4.0_Kr * Xi%X

         GradPhiHat(6,:)%X = 4.0_Kr * Xi%Y
         GradPhiHat(6,:)%Y = 4.0_Kr * Xi%X
         GradPhiHat(6,:)%Z = 0.0_Kr

         GradPhiHat(7,:)%X = -4.0_Kr * Xi%Y
         GradPhiHat(7,:)%Y = 4.0_Kr * (1.0_Kr - Xi%X - 2.0_Kr * Xi%Y - Xi%Z)
         GradPhiHat(7,:)%Z = -4.0_Kr * Xi%Y

         GradPhiHat(8,:)%X = -4.0_Kr * Xi%Z
         GradPhiHat(8,:)%Y = -4.0_Kr * Xi%Z
         GradPhiHat(8,:)%Z = 4.0_Kr * (1.0_Kr - Xi%X - Xi%Y -2.0_Kr * Xi%Z)

         GradPhiHat(9,:)%X = 4.0_Kr * Xi%Z
         GradPhiHat(9,:)%Y = 0.0_Kr
         GradPhiHat(9,:)%Z = 4.0_Kr * Xi%X
          
         GradPhiHat(10,:)%X = 0.0_Kr
         GradPhiHat(10,:)%Y = 4.0_Kr * Xi%Z
         GradPhiHat(10,:)%Z = 4.0_Kr * Xi%Y
      Case Default
         Print*,__FUNCT__,': Unimplemented PolynomialOrder',dPolynomialOrder
      End Select
      
      Allocate (dElem%BF(Num_DoF,Nb_Gauss)) 
      Allocate (dElem%Grad_BF(Num_DoF,Nb_Gauss))
      dElem%BF = PhiHat
      Do iDoF = 1,Num_DoF
         Do iG = 1,Nb_Gauss
            dElem%Grad_BF(iDoF,iG) = Bt * GradPhiHat(iDoF,iG) 
         End Do
      End Do
     
      DeAllocate(Xi)
      DeAllocate(PhiHat)
      DeAllocate(GradPhiHat)
   End Subroutine Element_P_Lagrange_3D_Scal_Init

#undef __FUNCT__
#define __FUNCT__ "Element_P_Lagrange_3DBoundary_Scal_Init"
   Subroutine Element_P_Lagrange_3DBoundary_Scal_Init(dElem,dCoord,dPolynomialOrder,dQuadratureOrder)
      Type(Element3D_Scal)                   :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord      ! coord(i,j)=ith coord of jth vertice
      PetscInt,Intent(IN)                    :: dPolynomialOrder,dQuadratureOrder
      
      Type(Element3D_Scal)                   :: tmpElem
      PetscReal,Dimension(:,:),Pointer       :: tmpCoord
      PetscInt                               :: i,j,iDoF,iG,Num_Gauss,Num_DoF
      Type(Vect3D)                           :: Edge1,Edge2,NormalVector
      PetscReal                              :: InnerBF
      
      !!!
      !!! Create a bogus tet element of unit volume by adding a 4th vertex along the normal of the
      !!! face, so that if a function is constant along the normal
      !!! direction of the face, one has \int_face fdx = \int_tet fdx 
      !!!
      Allocate(tmpCoord(3,4))
      Do i = 1,3
         Do j = 1,3
            tmpCoord(i,j) = dCoord(i,j)
         End Do
      End Do
      Edge1 = dCoord(:,2) - dCoord(:,1)
      Edge2 = dCoord(:,3) - dCoord(:,1)
      NormalVector = CrossP3D(Edge1,Edge2)
      NormalVector = NormalVector / Norm(NormalVector)
      tmpCoord(1,4) = dCoord(1,1) + NormalVector%X 
      tmpCoord(2,4) = dCoord(2,1) + NormalVector%Y 
      tmpCoord(3,4) = dCoord(3,1) + NormalVector%Z 
      

      Call Element_P_Lagrange_3D_Scal_Init(tmpElem,tmpCoord,dPolynomialOrder,dQuadratureOrder)
      Select Case (dPolynomialOrder)
         Case (1)
            Num_DoF  = 3
            Num_Gauss = size(tmpElem%BF,2)
            Allocate(dElem%Gauss_C(Num_Gauss))
            Allocate(dElem%BF(Num_DoF,Num_Gauss))
            dElem%Gauss_C = tmpElem%Gauss_C * 3.0_Kr
            Do iDoF = 1,Num_doF
               Do iG = 1,Num_Gauss
                  dElem%BF(iDoF,iG) = (tmpElem%BF(iDoF,iG) + tmpElem%BF(Num_DoF+1,iG) / 3.0_Kr)
               End Do
            End Do
         Case (2)
            !!! dof corrrespondance between a TRI6 and  a TETRA10 element (exodus.pdf figure 4 p.  17)
            !!! TRI6 TETRA10
            !!! 1    1
            !!! 2    2
            !!! 3    3
            !!! 4    6
            !!! 5    7
            !!! 6    5
            !!! Inner dof: 4,8,9,10
            Num_DoF  = 6
            Num_Gauss = size(tmpElem%BF,2)
            Allocate(dElem%Gauss_C(Num_Gauss))
            Allocate(dElem%BF(Num_DoF,Num_Gauss))
            dElem%Gauss_C = tmpElem%Gauss_C * 3.0_Kr
            Do iG = 1,Num_Gauss
               InnerBF = (tmpElem%BF(4,iG) + tmpElem%BF(8,iG) + tmpElem%BF(9,iG) + tmpElem%BF(10,iG)) *.25_Kr
               !!! This does not seem right,
               !!! I the inner BF must be distributed another way, but how...
               !!! I suspect 
               !!!      4,9 -> 2
               !!!      8 ->5
               !!!      10 -> 6
               dElem%BF(1,iG) = tmpElem%BF(1,iG) + InnerBF
               dElem%BF(2,iG) = tmpElem%BF(2,iG) + InnerBF
               dElem%BF(3,iG) = tmpElem%BF(3,iG) + InnerBF
               dElem%BF(4,iG) = tmpElem%BF(6,iG) + InnerBF
               dElem%BF(5,iG) = tmpElem%BF(7,iG) + InnerBF
               dElem%BF(6,iG) = tmpElem%BF(5,iG) + InnerBF
            End Do
         Case Default
            Print*,'[ERROR]: Polynomial order ',dPolynomialOrder,' not implemented in ',__FUNCT__
      End Select
      Call Element3D_Scal_Destroy(tmpElem)
      deAllocate(tmpCoord)
      Allocate(dElem%Grad_BF(0,0))
   End Subroutine Element_P_Lagrange_3DBoundary_Scal_Init
   
#undef __FUNCT__
#define __FUNCT__ "Element_P_Lagrange_3D_Vect_Init"
   Subroutine Element_P_Lagrange_3D_Vect_Init(dElem,dCoord,dPolynomialOrder,dQuadratureOrder)
      Type(Element3D_Vect)                   :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord      ! coord(i,j)=ith coord of jth vertice
      PetscInt                               :: dPolynomialOrder,dQuadratureOrder
   
      Type(Element3D_Scal)                   :: Elem_Scal
      PetscInt                               :: dim = 3 
      PetscInt                               :: Num_DoF,Nb_Gauss,i
      
      
      Call Element_P_Lagrange_3D_Scal_Init(Elem_Scal,dCoord,dPolynomialOrder,dQuadratureOrder)
      Num_DoF   = Size(Elem_Scal%BF,1) 
      Nb_Gauss = Size(Elem_Scal%BF,2)
      Allocate(dElem%Gauss_C(Nb_Gauss))
      Allocate(dElem%BF(Num_DoF * dim,Nb_Gauss))
      Allocate(dElem%Der_BF(Num_DoF * dim,Nb_Gauss))
         
      dElem%Gauss_C = Elem_Scal%Gauss_C
      dElem%BF(:,:)%X = 0.0_Kr
      dElem%BF(:,:)%Y = 0.0_Kr
      dElem%BF(:,:)%Z = 0.0_Kr
      dElem%Der_BF(:,:)%XX = 0.0_Kr
      dElem%Der_BF(:,:)%XY = 0.0_Kr
      dElem%Der_BF(:,:)%XZ = 0.0_Kr
      dElem%Der_BF(:,:)%YX = 0.0_Kr
      dElem%Der_BF(:,:)%YY = 0.0_Kr
      dElem%Der_BF(:,:)%YZ = 0.0_Kr
      dElem%Der_BF(:,:)%ZX = 0.0_Kr
      dElem%Der_BF(:,:)%ZY = 0.0_Kr
      dElem%Der_BF(:,:)%ZZ = 0.0_Kr
      
      Do i = 0,Num_DoF-1
         dElem%BF(i*dim+1,:)%X = Elem_Scal%BF(i+1,:)
         dElem%BF(i*dim+2,:)%Y = Elem_Scal%BF(i+1,:)
         dElem%BF(i*dim+3,:)%Z = Elem_Scal%BF(i+1,:)
         dElem%Der_BF(i*dim+1,:)%XX = Elem_Scal%Grad_BF(i+1,:)%X
         dElem%Der_BF(i*dim+1,:)%XY = Elem_Scal%Grad_BF(i+1,:)%Y
         dElem%Der_BF(i*dim+1,:)%XZ = Elem_Scal%Grad_BF(i+1,:)%Z
         dElem%Der_BF(i*dim+2,:)%YX = Elem_Scal%Grad_BF(i+1,:)%X
         dElem%Der_BF(i*dim+2,:)%YY = Elem_Scal%Grad_BF(i+1,:)%Y
         dElem%Der_BF(i*dim+2,:)%YZ = Elem_Scal%Grad_BF(i+1,:)%Z
         dElem%Der_BF(i*dim+3,:)%ZX = Elem_Scal%Grad_BF(i+1,:)%X
         dElem%Der_BF(i*dim+3,:)%ZY = Elem_Scal%Grad_BF(i+1,:)%Y
         dElem%Der_BF(i*dim+3,:)%ZZ = Elem_Scal%Grad_BF(i+1,:)%Z
      End Do
      Call ElementDestroy(Elem_Scal)
   End Subroutine Element_P_Lagrange_3D_Vect_Init

#undef __FUNCT__
#define __FUNCT__ "Element_P_Lagrange_3DBoundary_Vect_Init"
   Subroutine Element_P_Lagrange_3DBoundary_Vect_Init(dElem,dCoord,dPolynomialOrder,dQuadratureOrder)
      Type(Element3D_Vect)                   :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord      ! coord(i,j)=ith coord of jth vertice
      PetscInt                               :: dPolynomialOrder,dQuadratureOrder
   
      Type(Element3D_Scal)                   :: Elem_Scal
      PetscInt                               :: dim = 3 
      PetscInt                               :: Num_DoF,Nb_Gauss,i
      
      
      Call Element_P_Lagrange_3DBoundary_Scal_Init(Elem_Scal,dCoord,dPolynomialOrder,dQuadratureOrder)
      Num_DoF   = Size(Elem_Scal%BF,1) 
      Nb_Gauss = Size(Elem_Scal%BF,2)
      Allocate(dElem%Gauss_C(Nb_Gauss))
      Allocate(dElem%BF(Num_DoF * dim,Nb_Gauss))
      Allocate(dElem%Der_BF(0,Nb_Gauss))
         
      dElem%Gauss_C = Elem_Scal%Gauss_C
      dElem%BF(:,:)%X = 0.0_Kr
      dElem%BF(:,:)%Y = 0.0_Kr
      dElem%BF(:,:)%Z = 0.0_Kr
      
      Do i = 0,Num_DoF-1
         dElem%BF(i*dim+1,:)%X = Elem_Scal%BF(i+1,:)
         dElem%BF(i*dim+2,:)%Y = Elem_Scal%BF(i+1,:)
         dElem%BF(i*dim+3,:)%Z = Elem_Scal%BF(i+1,:)
      End Do
      Call ElementDestroy(Elem_Scal)
   End Subroutine Element_P_Lagrange_3DBoundary_Vect_Init

#undef __FUNCT__
#define __FUNCT__ "Element_P_Lagrange_3D_Elast_Init"
   Subroutine Element_P_Lagrange_3D_Elast_Init(dElem,dCoord,dPolynomialOrder,dQuadratureOrder)
      Type(Element3D_Elast)                  :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord      ! coord(i,j)=ith coord of jth vertice
      PetscInt                               :: dPolynomialOrder,dQuadratureOrder
   
      Type(Element3D_Scal)                   :: Elem_Scal
      PetscInt                               :: dim = 3 
      PetscInt                               :: Num_DoF,Nb_Gauss,i
      
      
      Call Element_P_Lagrange_3D_Scal_Init(Elem_Scal,dCoord,dPolynomialOrder,dQuadratureOrder)
      Num_DoF  = Size(Elem_Scal%BF,1) 
      Nb_Gauss = Size(Elem_Scal%BF,2)
      Allocate(dElem%Gauss_C(Nb_Gauss))
      Allocate(dElem%BF(Num_DoF * dim,Nb_Gauss))
      Allocate(dElem%GradS_BF(Num_DoF * dim,Nb_Gauss))
         
      dElem%Gauss_C = Elem_Scal%Gauss_C
      dElem%BF(:,:)%X = 0.0_Kr
      dElem%BF(:,:)%Y = 0.0_Kr
      dElem%BF(:,:)%Z = 0.0_Kr
      dElem%GradS_BF(:,:)%XX = 0.0_Kr
      dElem%GradS_BF(:,:)%YY = 0.0_Kr
      dElem%GradS_BF(:,:)%ZZ = 0.0_Kr
      dElem%GradS_BF(:,:)%YZ = 0.0_Kr
      dElem%GradS_BF(:,:)%XZ = 0.0_Kr
      dElem%GradS_BF(:,:)%XY = 0.0_Kr
      
      Do i = 0,Num_DoF-1
         dElem%BF(i*dim+1,:)%X = Elem_Scal%BF(i+1,:)
         dElem%BF(i*dim+2,:)%Y = Elem_Scal%BF(i+1,:)
         dElem%BF(i*dim+3,:)%Z = Elem_Scal%BF(i+1,:)
         dElem%GradS_BF(i*dim+1,:)%XX = Elem_Scal%Grad_BF(i+1,:)%X
         dElem%GradS_BF(i*dim+1,:)%XY = Elem_Scal%Grad_BF(i+1,:)%Y * 0.5_Kr
         dElem%GradS_BF(i*dim+1,:)%XZ = Elem_Scal%Grad_BF(i+1,:)%Z * 0.5_Kr

         dElem%GradS_BF(i*dim+2,:)%XY = Elem_Scal%Grad_BF(i+1,:)%X * 0.5_Kr
         dElem%GradS_BF(i*dim+2,:)%YY = Elem_Scal%Grad_BF(i+1,:)%Y 
         dElem%GradS_BF(i*dim+2,:)%YZ = Elem_Scal%Grad_BF(i+1,:)%Z * 0.5_Kr
         
         dElem%GradS_BF(i*dim+3,:)%XZ = Elem_Scal%Grad_BF(i+1,:)%X * 0.5_Kr
         dElem%GradS_BF(i*dim+3,:)%YZ = Elem_Scal%Grad_BF(i+1,:)%Y * 0.5_Kr
         dElem%GradS_BF(i*dim+3,:)%ZZ = Elem_Scal%Grad_BF(i+1,:)%Z
      End Do
      Call ElementDestroy(Elem_Scal)
   End Subroutine Element_P_Lagrange_3D_Elast_Init
   
#undef __FUNCT__
#define __FUNCT__ "Element_P_Lagrange_3DBoundary_Elast_Init"
   Subroutine Element_P_Lagrange_3DBoundary_Elast_Init(dElem,dCoord,dPolynomialOrder,dQuadratureOrder)
      Type(Element3D_Elast)                  :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord      ! coord(i,j)=ith coord of jth vertice
      PetscInt                               :: dPolynomialOrder,dQuadratureOrder
   
      Type(Element3D_Scal)                   :: Elem_Scal
      PetscInt                               :: dim = 3 
      PetscInt                               :: Num_DoF,Nb_Gauss,i
      
      
      Call Element_P_Lagrange_3DBoundary_Scal_Init(Elem_Scal,dCoord,dPolynomialOrder,dQuadratureOrder)
      Num_DoF   = Size(Elem_Scal%BF,1) 
      Nb_Gauss = Size(Elem_Scal%BF,2)
      Allocate(dElem%Gauss_C(Nb_Gauss))
      Allocate(dElem%BF(Num_DoF * dim,Nb_Gauss))
      Allocate(dElem%GradS_BF(0,Nb_Gauss))
         
      dElem%Gauss_C = Elem_Scal%Gauss_C
      dElem%BF(:,:)%X = 0.0_Kr
      dElem%BF(:,:)%Y = 0.0_Kr
      dElem%BF(:,:)%Z = 0.0_Kr
      
      Do i = 0,Num_DoF-1
         dElem%BF(i*dim+1,:)%X = Elem_Scal%BF(i+1,:)
         dElem%BF(i*dim+2,:)%Y = Elem_Scal%BF(i+1,:)
         dElem%BF(i*dim+3,:)%Z = Elem_Scal%BF(i+1,:)
      End Do
      Call ElementDestroy(Elem_Scal)
   End Subroutine Element_P_Lagrange_3DBoundary_Elast_Init

#undef __FUNCT__
#define __FUNCT__ "Element2D_Scal_Destroy"
   Subroutine Element2D_Scal_Destroy(dElem)
      Type(Element2D_Scal)                   :: dElem
      
      If (Associated(dElem%BF)) Then
         DeAllocate(dElem%BF)
      End If
      If (Associated(dElem%Grad_BF)) Then
         DeAllocate(dElem%Grad_BF)
      End If
      If (Associated(dElem%Gauss_C)) Then
         DeAllocate(dElem%Gauss_C)
      End If
   End Subroutine Element2D_Scal_Destroy

#undef __FUNCT__
#define __FUNCT__ "Element2D_Destroy"
   Subroutine Element2D_Vect_Destroy(dElem)
      Type(Element2D_Vect)                        :: dElem
      If (Associated(dElem%BF)) Then
         DeAllocate(dElem%BF)
      End If
      If (Associated(dElem%Der_BF)) Then
         DeAllocate(dElem%Der_BF)
      End If
      If (Associated(dElem%Gauss_C)) Then
         DeAllocate(dElem%Gauss_C)
      End If
   End Subroutine Element2D_Vect_Destroy
   
#undef __FUNCT__
#define __FUNCT__ "Element2D_Elast_Destroy"
   Subroutine Element2D_Elast_Destroy(dElem)
      Type(Element2D_Elast)                  :: dElem
      If (Associated(dElem%BF)) Then
         DeAllocate(dElem%BF)
      End If
      If (Associated(dElem%GradS_BF)) Then
         DeAllocate(dElem%GradS_BF)
      End If
      If (Associated(dElem%Gauss_C)) Then
         DeAllocate(dElem%Gauss_C)
      End If
   End Subroutine Element2D_Elast_Destroy
   
#undef __FUNCT__
#define __FUNCT__ "Element3D_Scal_Destroy"
   Subroutine Element3D_Scal_Destroy(dElem)
      Type(Element3D_Scal)                   :: dElem
      If (Associated(dElem%BF)) Then
         DeAllocate(dElem%BF)
      End If
      If (Associated(dElem%Grad_BF)) Then
         DeAllocate(dElem%Grad_BF)
      End If
      If (Associated(dElem%Gauss_C)) Then
         DeAllocate(dElem%Gauss_C)
      End If
   End Subroutine Element3D_Scal_Destroy

#undef __FUNCT__
#define __FUNCT__ "Element3D_Vect_Destroy"
   Subroutine Element3D_Vect_Destroy(dElem)
      Type(Element3D_Vect)                        :: dElem
      If (Associated(dElem%BF)) Then
         DeAllocate(dElem%BF)
      End If
      If (Associated(dElem%Der_BF)) Then
         DeAllocate(dElem%Der_BF)
      End If
      If (Associated(dElem%Gauss_C)) Then
         DeAllocate(dElem%Gauss_C)
      End If
   End Subroutine Element3D_Vect_Destroy
   
#undef __FUNCT__
#define __FUNCT__ "Element3D_Elast_Destroy"
   Subroutine Element3D_Elast_Destroy(dElem)
      Type(Element3D_Elast)                  :: dElem
      If (Associated(dElem%BF)) Then
         DeAllocate(dElem%BF)
      End If
      If (Associated(dElem%GradS_BF)) Then
         DeAllocate(dElem%GradS_BF)
      End If
      If (Associated(dElem%Gauss_C)) Then
         DeAllocate(dElem%Gauss_C)
      End If
   End Subroutine Element3D_Elast_Destroy

!!!
!!! ELEMENT VIEWERS
!!!

#undef __FUNCT__
#define __FUNCT__ "Element2D_Scal_View"
   Subroutine Element2D_Scal_View(dElem,viewer)
      Type(Element2D_Scal)                            :: dElem
      Type(PetscViewer)                               :: viewer
      
      PetscInt                                        :: Nb_Gauss,Nb_DoF,iDoF,iG,ierr
      Character(len=512)                              :: CharBuffer
                
      Nb_DoF   = size(dElem%BF,1)
      Nb_Gauss = size(dElem%BF,2)
      Write(CharBuffer,102) Nb_DoF
      Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
      Write(CharBuffer,103) Nb_Gauss
      Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)

      Nb_Gauss = Size(dElem%BF,2)
      Do iDoF = 1,Nb_DoF
         Write(CharBuffer,104) iDoF
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)

         Write(CharBuffer,200) '        BF \n    '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%BF(iDoF,iG)
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n'
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)

         Write(CharBuffer,200) '        Grad_BF (X,Y)\n    '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%Grad_BF(iDoF,iG)%X
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n    '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%Grad_BF(iDoF,iG)%Y
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n'
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
      End Do
102 Format('    Nb_DoF   ',I9,'\n')
103 Format('    Nb_Gauss ',I9,'\n')
104 Format('    *** DoF  ',I9,'\n')
200 Format(A)
201 Format('   ',F5.2)
   End Subroutine Element2D_Scal_View
   
#undef __FUNCT__
#define __FUNCT__ "Element2D_Vect_View"
   Subroutine Element2D_Vect_View(dElem,viewer)
      Type(Element2D_Vect)                            :: dElem
      Type(PetscViewer)                               :: viewer
      
      PetscInt                                        :: Nb_Gauss,Nb_DoF,iDoF,iG,ierr
      Character(len=512)                              :: CharBuffer

      Nb_DoF   = size(dElem%BF,1)
      Nb_Gauss = size(dElem%BF,2)
      Write(CharBuffer,102) Nb_DoF
      Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
      Write(CharBuffer,103) Nb_Gauss
      Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)

      Nb_Gauss = Size(dElem%BF,2)
      Do iDoF = 1,Nb_DoF
         Write(CharBuffer,104) iDoF
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)

         Write(CharBuffer,200) '        BF (X,Y)\n    '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%BF(iDoF,iG)%X
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n   '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%BF(iDoF,iG)%Y
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n'
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)

         Write(CharBuffer,200) '        Der_BF (XX,XY,YX,YY)\n    '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%Der_BF(iDoF,iG)%XX
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n    '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%Der_BF(iDoF,iG)%XY
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n'
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%Der_BF(iDoF,iG)%YX
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n    '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%Der_BF(iDoF,iG)%YY
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n    '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
      End Do
102 Format('    Nb_DoF   ',I9,'\n')
103 Format('    Nb_Gauss ',I9,'\n')
104 Format('    *** DoF  ',I9,'\n')
200 Format(A)
201 Format('   ',F5.2)
   End Subroutine Element2D_Vect_View
   
   Subroutine Element2D_Elast_View(dElem,viewer)
      Type(Element2D_Elast)                           :: dElem
      Type(PetscViewer)                               :: viewer
      
      PetscInt                                        :: Nb_Gauss,Nb_DoF,iDoF,iG,ierr
      Character(len=512)                              :: CharBuffer

      Nb_DoF   = size(dElem%BF,1)
      Nb_Gauss = size(dElem%BF,2)
      Write(CharBuffer,102) Nb_DoF
      Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
      Write(CharBuffer,103) Nb_Gauss
      Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)

      Nb_Gauss = Size(dElem%BF,2)
      Do iDoF = 1,Nb_DoF
         Write(CharBuffer,104) iDoF
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)

         Write(CharBuffer,200) '        BF (X,Y)\n    '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%BF(iDoF,iG)%X
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n   '
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%BF(iDoF,iG)%Y
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n   '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)

         Write(CharBuffer,200) '        GradS_BF (XX,YY,XY)\n    '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%GradS_BF(iDoF,iG)%XX
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n    '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%GradS_BF(iDoF,iG)%YY
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n'
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%GradS_BF(iDoF,iG)%XY
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n    '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
      End Do
102 Format('    Nb_DoF   ',I9,'\n')
103 Format('    Nb_Gauss ',I9,'\n')
104 Format('    *** DoF  ',I9,'\n')
200 Format(A)
201 Format('   ',F5.2)
   End Subroutine Element2D_Elast_View

#undef __FUNCT__
#define __FUNCT__ "Element3D_Scal_View"
   Subroutine Element3D_Scal_View(dElem,viewer)
      Type(Element3D_Scal)                            :: dElem
      Type(PetscViewer)                               :: viewer
      
      PetscInt                                        :: Nb_Gauss,Nb_DoF,iDoF,iG,ierr
      Character(len=512)                              :: CharBuffer
                
      Nb_DoF   = size(dElem%BF,1)
      Nb_Gauss = size(dElem%BF,2)
      Write(CharBuffer,102) Nb_DoF
      Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
      Write(CharBuffer,103) Nb_Gauss
      Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)

      Nb_Gauss = Size(dElem%BF,2)
      Do iDoF = 1,Nb_DoF
         Write(CharBuffer,104) iDoF
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)

         Write(CharBuffer,200) '        BF \n    '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%BF(iDoF,iG)
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n'
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)

         Write(CharBuffer,200) '        Grad_BF (X,Y,Z)\n    '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%Grad_BF(iDoF,iG)%X
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n    '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%Grad_BF(iDoF,iG)%Y
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n    '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%Grad_BF(iDoF,iG)%Z
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n    '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)         
      End Do
102 Format('    Nb_DoF   ',I9,'\n')
103 Format('    Nb_Gauss ',I9,'\n')
104 Format('    *** DoF  ',I9,'\n')
200 Format(A)
201 Format('   ',F5.2)
   End Subroutine Element3D_Scal_View
   
#undef __FUNCT__
#define __FUNCT__ "Element3D_View"
   Subroutine Element3D_Vect_View(dElem,viewer)
      Type(Element3D_Vect)                            :: dElem
      Type(PetscViewer)                               :: viewer
      
      PetscInt                                        :: Nb_Gauss,Nb_DoF,iDoF,iG,ierr
      Character(len=512)                              :: CharBuffer

      Nb_DoF   = size(dElem%BF,1)
      Nb_Gauss = size(dElem%BF,2)
      Write(CharBuffer,102) Nb_DoF
      Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
      Write(CharBuffer,103) Nb_Gauss
      Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)

      Nb_Gauss = Size(dElem%BF,2)
      Do iDoF = 1,Nb_DoF
         Write(CharBuffer,104) iDoF
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)

         Write(CharBuffer,200) '        BF (X,Y,Z)\n    '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%BF(iDoF,iG)%X
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n   '
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%BF(iDoF,iG)%Y
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n   '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%BF(iDoF,iG)%Z
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n'
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)

         Write(CharBuffer,200) '        Der_BF (XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ)\n    '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%Der_BF(iDoF,iG)%XX
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n    '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%Der_BF(iDoF,iG)%XY
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n'
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%Der_BF(iDoF,iG)%XZ
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n'
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%Der_BF(iDoF,iG)%YX
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n    '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%Der_BF(iDoF,iG)%YY
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n    '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%Der_BF(iDoF,iG)%YZ
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n'
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%Der_BF(iDoF,iG)%ZX
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n'
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%Der_BF(iDoF,iG)%ZY
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n'
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%Der_BF(iDoF,iG)%ZZ
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n'
      End Do
102 Format('    Nb_DoF   ',I9,'\n')
103 Format('    Nb_Gauss ',I9,'\n')
104 Format('    *** DoF  ',I9,'\n')
200 Format(A)
201 Format('   ',F5.2)
   End Subroutine Element3D_Vect_View
   
   Subroutine Element3D_Elast_View(dElem,viewer)
      Type(Element3D_Elast)                           :: dElem
      Type(PetscViewer)                               :: viewer
      
      PetscInt                                        :: Nb_Gauss,Nb_DoF,iDoF,iG,ierr
      Character(len=512)                              :: CharBuffer

      Nb_DoF   = size(dElem%BF,1)
      Nb_Gauss = size(dElem%BF,2)
      Write(CharBuffer,102) Nb_DoF
      Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
      Write(CharBuffer,103) Nb_Gauss
      Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)

      Nb_Gauss = Size(dElem%BF,2)
      Do iDoF = 1,Nb_DoF
         Write(CharBuffer,104) iDoF
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)

         Write(CharBuffer,200) '        BF (X,Y,Z)\n    '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%BF(iDoF,iG)%X
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n   '
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%BF(iDoF,iG)%Y
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n   '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%BF(iDoF,iG)%Z
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n'
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)

         Write(CharBuffer,200) '        GradS_BF (XX,YY,ZZ,YZ,XZ,XY)\n    '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%GradS_BF(iDoF,iG)%XX
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n    '
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%GradS_BF(iDoF,iG)%YY
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n'
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%GradS_BF(iDoF,iG)%ZZ
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n    '
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%GradS_BF(iDoF,iG)%YZ
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n    '
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%GradS_BF(iDoF,iG)%XZ
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n    '
         Do iG = 1,Nb_Gauss
            Write(CharBuffer,201) dElem%GradS_BF(iDoF,iG)%XY
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(CharBuffer,*) '\n    '
      End Do
102 Format('    Nb_DoF   ',I9,'\n')
103 Format('    Nb_Gauss ',I9,'\n')
104 Format('    *** DoF  ',I9,'\n')
200 Format(A)
201 Format('   ',F5.2)
   End Subroutine Element3D_Elast_View
End Module m_MEF_Elements
