!!! TODO
!!! - Add a validation of the shortID and names
Module m_MEF90_Elements
#include "petsc/finclude/petsc.h"
   Use m_MEF90_LinAlg
   Use m_MEF90_Utils
   Use m_MEF90_Parameters
   Use petsc
   Use,intrinsic :: iso_c_binding
   IMPLICIT NONE

   !Private   
   ! not sure how to make the enumerator public short of listing them one after each other.
   Public :: MEF90Element_Create
   Public :: MEF90Element_Destroy

   Public :: MEF90Element_Type
   Public :: MEF90Element2D_Vect,MEF90Element2D_Scal
   Public :: MEF90Element3D_Vect,MEF90Element3D_Scal
   Public :: MEF90Element_TypeFindByID,MEF90Element_TypeFindByName
   
   Type MEF90Element_Type
      ! name is the element name in english language
      Character(len=MEF90_MXSTRLEN)                :: Name
      ! shortID is a numerical shortcut used in the prep, for instance
      PetscInt                                     :: ShortID
      
      !!! Geometry of the cell associated with the element
      PetscInt                                     :: numVertex
      PetscInt                                     :: numEdge
      PetscInt                                     :: numFaces
      
      !!! Number of dof at each 
      !!!   vertex edge face cell in 3D,
      !!!   vertex edge cell in 2D (last value is unused, i.e. always use depth to locate DoF)
      !!! i.e. point depth -1
      PetscInt,Dimension(4)                        :: numDofs
      ! !!!
      ! PetscInt                                     :: numVertexDof
      ! PetscInt                                     :: numEdgeDof
      ! PetscInt                                     :: numFaceDof
      ! PetscInt                                     :: numCellDof
      PetscInt                                     :: numDoF 
      !!! numDof = numVertexDof * numVertex + numEdgeDof * numEdge 
      !!!          + numFaceDof * numFace + numCellDof
      PetscInt                                     :: dim
      !!! The dimension of the cell associated with the element
      PetscInt                                     :: coDim
      !!! The co-dimension of the cell associated with the element
      PetscInt                                     :: order
      !!! The polynomial order
   End Type MEF90Element_Type
 
   Enum,bind(c)
      enumerator  :: &
      MEF90_P1_Lagrange_2D_ShortID = 1,       &  ! 1
      MEF90_P1_Lagrange_3D_ShortID,           &  ! 2
      MEF90_P1_Lagrange_2DBoundary_ShortID,   &  ! 3
      MEF90_P1_Lagrange_3DBoundary_ShortID,   &  ! 4
      MEF90_P2_Lagrange_2D_ShortID,           &  ! 5
      MEF90_P2_Lagrange_3D_ShortID,           &  ! 6
      MEF90_P2_Lagrange_2DBoundary_ShortID,   &  ! 7
      MEF90_P2_Lagrange_3DBoundary_ShortID,   &  ! 8
      MEF90_Q1_Lagrange_2D_ShortID,           &  ! 9
      MEF90_Q1_Lagrange_3D_ShortID,           &  ! 10
      MEF90_Q1_Lagrange_2DBoundary_ShortID,   &  ! 11
      MEF90_Q1_Lagrange_3DBoundary_ShortID,   &  ! 12
      MEF90_Q2_Lagrange_2D_ShortID,           &  ! 13
      MEF90_Q2_Lagrange_3D_ShortID,           &  ! 14
      MEF90_Q2_Lagrange_2DBoundary_ShortID,   &  ! 15
      MEF90_Q2_Lagrange_3DBoundary_ShortID       ! 16
End Enum      
   !!! 
   !!! P1 Simplicial Lagrange elements
   !!!
Type(MEF90Element_Type),Parameter,Public :: MEF90_P1_Lagrange_2D = MEF90Element_Type(   &
"MEF90_P1_Lagrange_2D",                      &  ! name
MEF90_P1_Lagrange_2D_ShortID,                &  ! shortID
3,3,0,                                       &  ! numVertex,numEdge,numFace
[1,0,0,0],3,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
2,0,1                                        &  ! dim,codim,order
)
Type(MEF90Element_Type),Parameter,Public :: MEF90_P1_Lagrange_3D = MEF90Element_Type(   &
"MEF90_P1_Lagrange_3D",                      &  ! name
MEF90_P1_Lagrange_3D_ShortID,                &  ! shortID
4,6,4,                                       &  ! numVertex,numEdge,numFace
[1,0,0,0],4,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
3,0,1                                        &  ! dim,codim,order
)
Type(MEF90Element_Type),Parameter,Public :: MEF90_P1_Lagrange_2DBoundary = MEF90Element_Type(   &
"MEF90_P1_Lagrange_2DBoundary",              &  ! name
MEF90_P1_Lagrange_2DBoundary_ShortID,        &  ! shortID
2,1,0,                                       &  ! numVertex,numEdge,numFace
[1,0,0,0],2,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
2,1,1                                        &  ! dim,codim,order
)
Type(MEF90Element_Type),Parameter,Public :: MEF90_P1_Lagrange_3DBoundary = MEF90Element_Type(   &
"MEF90_P1_Lagrange_3DBoundary",              &  ! name
MEF90_P1_Lagrange_3DBoundary_ShortID,        &  ! shortID
3,3,0,                                       &  ! numVertex,numEdge,numFace
[1,0,0,0],3,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
3,1,1                                        &  ! dim,codim,order
)
!!!
!!! P2 Simplicial Lagrange elements
!!!
Type(MEF90Element_Type),Parameter,Public :: MEF90_P2_Lagrange_2D = MEF90Element_Type(   &
"MEF90_P2_Lagrange_2D",                      &  ! name
MEF90_P2_Lagrange_2D_ShortID,                &  ! shortID
3,3,0,                                       &  ! numVertex,numEdge,numFace
[1,1,0,0],6,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
2,0,2                                        &  ! dim,codim,order
)
Type(MEF90Element_Type),Parameter,Public :: MEF90_P2_Lagrange_3D = MEF90Element_Type(   &
"MEF90_P2_Lagrange_3D",                      &  ! name
MEF90_P2_Lagrange_3D_ShortID,                &  ! shortID
4,6,4,                                       &  ! numVertex,numEdge,numFace
[1,1,0,0],10,                                &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
3,0,2                                        &  ! dim,codim,order
)
Type(MEF90Element_Type),Parameter,Public :: MEF90_P2_Lagrange_2DBoundary = MEF90Element_Type(   &
"MEF90_P2_Lagrange_2DBoundary",              &  ! name
MEF90_P2_Lagrange_2DBoundary_ShortID,        &  ! shortID
2,1,0,                                       &  ! numVertex,numEdge,numFace
[1,1,0,0],3,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
2,1,2                                        &  ! dim,codim,order
)
Type(MEF90Element_Type),Parameter,Public :: MEF90_P2_Lagrange_3DBoundary = MEF90Element_Type(   &
"MEF90_P2_Lagrange_3DBoundary",              &  ! name
MEF90_P2_Lagrange_3DBoundary_ShortID,        &  ! shortID
3,3,0,                                       &  ! numVertex,numEdge,numFace
[1,1,0,0],6,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
3,1,2                                        &  ! dim,codim,order
)
!!! 
!!! Q1 tensor product Lagrange elements
!!!
Type(MEF90Element_Type),Parameter,Public :: MEF90_Q1_Lagrange_2D = MEF90Element_Type(   &
"MEF90_Q1_Lagrange_2D",                      &  ! name
MEF90_Q1_Lagrange_2D_ShortID,                &  ! shortID
4,4,0,                                       &  ! numVertex,numEdge,numFace
[1,0,0,0],4,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
2,0,1                                        &  ! dim,codim,order
)
Type(MEF90Element_Type),Parameter,Public :: MEF90_Q1_Lagrange_3D = MEF90Element_Type(   &
"MEF90_Q1_Lagrange_3D",                      &  ! name
MEF90_Q1_Lagrange_3D_ShortID,                &  ! shortID
8,12,6,                                      &  ! numVertex,numEdge,numFace
[1,0,0,0],8,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
3,0,1                                        &  ! dim,codim,order
)
Type(MEF90Element_Type),Parameter,Public :: MEF90_Q1_Lagrange_2DBoundary = MEF90Element_Type(   &
"MEF90_Q1_Lagrange_2DBoundary",              &  ! name
MEF90_Q1_Lagrange_2DBoundary_ShortID,        &  ! shortID
2,1,0,                                       &  ! numVertex,numEdge,numFace
[1,0,0,0],2,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
2,1,1                                        &  ! dim,codim,order
)
Type(MEF90Element_Type),Parameter,Public :: MEF90_Q1_Lagrange_3DBoundary = MEF90Element_Type(   &
"MEF90_Q1_Lagrange_3DBoundary",              &  ! name
MEF90_Q1_Lagrange_3DBoundary_ShortID,        &  ! shortID
4,4,0,                                       &  ! numVertex,numEdge,numFace
[1,0,0,0],4,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
3,1,1                                        &  ! dim,codim,order
)
!!! 
!!! Q2 tensor product Lagrange elements
!!!
Type(MEF90Element_Type),Parameter,Public :: MEF90_Q2_Lagrange_2D = MEF90Element_Type(   &
"MEF90_Q2_Lagrange_2D",                      &  ! name
MEF90_Q2_Lagrange_2D_ShortID,                &  ! shortID
4,4,0,                                       &  ! numVertex,numEdge,numFace
[1,1,0,1],9,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
2,0,1                                        &  ! dim,codim,order
)
Type(MEF90Element_Type),Parameter,Public :: MEF90_Q2_Lagrange_3D = MEF90Element_Type(   &
"MEF90_Q2_Lagrange_3D",                      &  ! name
MEF90_Q2_Lagrange_3D_ShortID,                &  ! shortID
8,12,6,                                      &  ! numVertex,numEdge,numFace
[1,1,1,1],27,                                &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
3,0,1                                        &  ! dim,codim,order
)
Type(MEF90Element_Type),Parameter,Public :: MEF90_Q2_Lagrange_2DBoundary = MEF90Element_Type(   &
"MEF90_Q2_Lagrange_2DBoundary",              &  ! name
MEF90_Q2_Lagrange_2DBoundary_ShortID,        &  ! shortID
2,1,0,                                       &  ! numVertex,numEdge,numFace
[1,1,0,0],3,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
2,1,1                                        &  ! dim,codim,order
)
Type(MEF90Element_Type),Parameter,Public :: MEF90_Q2_Lagrange_3DBoundary = MEF90Element_Type(   &
"MEF90_Q2_Lagrange_3DBoundary",              &  ! name
MEF90_Q2_Lagrange_3DBoundary_ShortID,        &  ! shortID
4,4,0,                                       &  ! numVertex,numEdge,numFace
[1,1,0,1],9,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
3,1,1                                        &  ! dim,codim,order
)
Type(MEF90Element_Type),Parameter,Public :: MEF90_NULL_ELEMENT = MEF90Element_Type(   &
"NULL",                                      &  ! name
0,                                           &  ! shortID
0,0,0,                                       &  ! numVertex,numEdge,numFace
[0,0,0,0],0,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
0,0,0                                        &  ! dim,codim,order
)

Integer,Parameter,Public :: MEF90_numKnownElements = 16       
Type(MEF90Element_Type),dimension(MEF90_numKnownElements),Parameter,Public   :: MEF90knownElements = [ &
   MEF90_P1_Lagrange_2D,          &  ! 1
   MEF90_P1_Lagrange_3D,          &  ! 2
   MEF90_P1_Lagrange_2DBoundary,  &  ! 3
   MEF90_P1_Lagrange_3DBoundary,  &  ! 4
   MEF90_P2_Lagrange_2D,          &  ! 5
   MEF90_P2_Lagrange_3D,          &  ! 6
   MEF90_P2_Lagrange_2DBoundary,  &  ! 7
   MEF90_P2_Lagrange_3DBoundary,  &  ! 8
   MEF90_Q1_Lagrange_2D,          &  ! 9
   MEF90_Q1_Lagrange_3D,          &  ! 10
   MEF90_Q1_Lagrange_2DBoundary,  &  ! 11
   MEF90_Q1_Lagrange_3DBoundary,  &  ! 12
   MEF90_Q2_Lagrange_2D,          &  ! 13
   MEF90_Q2_Lagrange_3D,          &  ! 14
   MEF90_Q2_Lagrange_2DBoundary,  &  ! 15
   MEF90_Q2_Lagrange_3DBoundary   &  ! 16
]

Character(kind=c_char,len=MEF90_MXSTRLEN),dimension(MEF90_numKnownElements+3),Parameter,Public   :: MEF90_knownElementNames = [ &
   "P1_Lagrange_2D         ",     &  ! 1
   "P1_Lagrange_3D         ",     &  ! 2
   "P1_Lagrange_2DBoundary ",     &  ! 3
   "P1_Lagrange_3DBoundary ",     &  ! 4
   "P2_Lagrange_2D         ",     &  ! 5
   "P2_Lagrange_3D         ",     &  ! 6
   "P2_Lagrange_2DBoundary ",     &  ! 7
   "P2_Lagrange_3DBoundary ",     &  ! 8
   "Q1_Lagrange_2D         ",     &  ! 9
   "Q1_Lagrange_3D         ",     &  ! 10
   "Q1_Lagrange_2DBoundary ",     &  ! 11
   "Q1_Lagrange_3DBoundary ",     &  ! 12
   "Q2_Lagrange_2D         ",     &  ! 13
   "Q2_Lagrange_3D         ",     &  ! 14
   "Q2_Lagrange_2DBoundary ",     &  ! 15
   "Q2_Lagrange_3DBoundary ",     &  ! 16
   "MEF90_knownElementNames",     &
   "prefix_                ",     &
   C_NULL_CHAR//"                      "]

      
   Type MEF90Element2D_Scal
      PetscReal,Dimension(:,:),Pointer             :: BF
      Type(Vect2D),Dimension(:,:),Pointer          :: Grad_BF
      PetscReal,Dimension(:),Pointer               :: Gauss_C
      Type(Vect2D)                                 :: outerNormal ! only makes sense for elements of codim 1
   End Type MEF90Element2D_Scal
 
   Type MEF90Element2D_Vect
      Type (Vect2D),Dimension(:,:),Pointer         :: BF
      Type (Mat2D),Dimension(:,:),Pointer          :: Grad_BF
      Type (MatS2D),Dimension(:,:),Pointer         :: GradS_BF
      PetscReal,Dimension(:),Pointer               :: Gauss_C
      Type(Vect2D)                                 :: outerNormal ! only makes sense for elements of codim 1
   End Type MEF90Element2D_Vect
 
   Type MEF90Element3D_Scal
      PetscReal,Dimension(:,:),Pointer             :: BF
      Type (Vect3D),Dimension(:,:),Pointer         :: Grad_BF
      PetscReal,Dimension(:),Pointer               :: Gauss_C
      Type(Vect3D)                                 :: outerNormal ! only makes sense for elements of codim 1
   End Type MEF90Element3D_Scal
 
   Type MEF90Element3D_Vect
      Type (Vect3D),Dimension(:,:),Pointer         :: BF
      Type (Mat3D),Dimension(:,:),Pointer          :: Grad_BF
      Type (MatS3D),Dimension(:,:),Pointer         :: GradS_BF
      PetscReal,Dimension(:),Pointer               :: Gauss_C
      Type(Vect3D)                                 :: outerNormal ! only makes sense for elements of codim 1
   End Type MEF90Element3D_Vect
 
   Interface MEF90Element_Create
      Module Procedure Element2D_Scal_InitSet,Element2D_Vect_InitSet, &
                       Element3D_Scal_InitSet,Element3D_Vect_InitSet, &
                       Element2D_Scal_InitSet_ByShortID,Element2D_Vect_InitSet_ByShortID,&
                       Element3D_Scal_InitSet_ByShortID,Element3D_Vect_InitSet_ByShortID
   End Interface MEF90Element_Create
   
   Interface MEF90Element_Destroy
      Module Procedure Element2D_Scal_Destroy,Element2D_Vect_Destroy,&
                       Element3D_Scal_Destroy,Element3D_Vect_Destroy,&
                       Element2D_Scal_DestroySet,Element2D_Vect_DestroySet,& 
                       Element3D_Scal_DestroySet,Element3D_Vect_DestroySet
   End Interface MEF90Element_Destroy
      
Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90Element_TypeFindByID"
   Subroutine MEF90Element_TypeFindByID(elemID,elemType,ierr)
      PetscInt, Intent(IN)                        :: elemID
      Type(MEF90Element_Type),Intent(OUT)         :: elemType
      PetscErrorCode,Intent(OUT)                  :: ierr
      
      Integer                                     :: i
      PetscBool                                   :: knownID = PETSC_FALSE      
      Do i = 1, size(MEF90knownElements)
         If (MEF90knownElements(i)%shortID == elemID) Then
            elemType = MEF90knownElements(i)
            knownID  = PETSC_TRUE
            EXIT
         End If
      End Do
      If (.NOT. knownID) Then
         Write(*,*) "[ERROR]: Unknown element ID", elemID
         ierr = PETSC_ERR_ARG_WRONG
      End If
   End Subroutine MEF90Element_TypeFindByID
   
#undef __FUNCT__
#define __FUNCT__ "MEF90Element_TypeFindByname"
   Subroutine MEF90Element_TypeFindByname(elemName,elemType,ierr)
      Character(len=*), Intent(IN)                :: elemName
      Type(MEF90Element_Type),Intent(OUT)         :: elemType
      PetscErrorCode,Intent(OUT)                  :: ierr

      Integer                                     :: i
      PetscBool                                   :: knownID = PETSC_FALSE      
      Do i = 1, size(MEF90knownElements)
         If (trim(MEF90knownElements(i)%name) == trim(elemName)) Then
            elemType = MEF90knownElements(i)
            knownID  = PETSC_TRUE
            EXIT
         End If
      End Do
      If (.NOT. knownID) Then
         Write(*,*) "[ERROR]: Unknown element ", trim(elemName)
         ierr = PETSC_ERR_ARG_WRONG
      End If
   End Subroutine MEF90Element_TypeFindByname

#undef __FUNCT__
#define __FUNCT__ "Element2D_Scal_InitSet_ByName"
   Subroutine Element2D_Scal_InitSet_ByName(mesh,cellIS,dElem,dQuadratureOrder,elemName,ierr)
      Type(tDM),intent(IN)                             :: mesh
      Type(tIS),intent(IN)                             :: cellIS
      Type(MEF90Element2D_Scal),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                              :: dQuadratureOrder
      Character(len=*), Intent(IN)                     :: elemName
      PetscErrorCode,Intent(OUT)                       :: ierr
      
      Type(MEF90Element_Type)                          :: elemType

      Call MEF90Element_TypeFindByName(elemName,elemType,ierr)
      Call Element2D_Scal_InitSet(mesh,cellIS,dElem,dQuadratureOrder,elemType,ierr)
   End Subroutine Element2D_Scal_InitSet_ByName

#undef __FUNCT__
#define __FUNCT__ "Element2D_Vect_InitSet_ByName"
   Subroutine Element2D_Vect_InitSet_ByName(mesh,cellIS,dElem,dQuadratureOrder,elemName,ierr)
      Type(tDM),intent(IN)                             :: mesh
      Type(tIS),intent(IN)                             :: cellIS
      Type(MEF90Element2D_Vect),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                              :: dQuadratureOrder
      Character(len=*), Intent(IN)                     :: elemName
      PetscErrorCode,Intent(OUT)                       :: ierr
      
      Type(MEF90Element_Type)                          :: elemType

      Call MEF90Element_TypeFindByName(elemName,elemType,ierr)
      Call Element2D_Vect_InitSet(mesh,cellIS,dElem,dQuadratureOrder,elemType,ierr)
   End Subroutine Element2D_Vect_InitSet_ByName

#undef __FUNCT__
#define __FUNCT__ "Element3D_Scal_InitSet_ByName"
   Subroutine Element3D_Scal_InitSet_ByName(mesh,cellIS,dElem,dQuadratureOrder,elemName,ierr)
      Type(tDM),intent(IN)                             :: mesh
      Type(tIS),intent(IN)                             :: cellIS
      Type(MEF90Element3D_Scal),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                              :: dQuadratureOrder
      Character(len=*), Intent(IN)                     :: elemName
      PetscErrorCode,Intent(OUT)                       :: ierr
      
      Type(MEF90Element_Type)                          :: elemType

      Call MEF90Element_TypeFindByName(elemName,elemType,ierr)
      Call Element3D_Scal_InitSet(mesh,cellIS,dElem,dQuadratureOrder,elemType,ierr)
   End Subroutine Element3D_Scal_InitSet_ByName

#undef __FUNCT__
#define __FUNCT__ "Element3D_Vect_InitSet_ByName"
   Subroutine Element3D_Vect_InitSet_ByName(mesh,cellIS,dElem,dQuadratureOrder,elemName,ierr)
      Type(tDM),intent(IN)                             :: mesh
      Type(tIS),intent(IN)                             :: cellIS
      Type(MEF90Element3D_Vect),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                              :: dQuadratureOrder
      Character(len=*), Intent(IN)                     :: elemName
      PetscErrorCode,Intent(OUT)                       :: ierr
      
      Type(MEF90Element_Type)                          :: elemType

      Call MEF90Element_TypeFindByName(elemName,elemType,ierr)
      Call Element3D_Vect_InitSet(mesh,cellIS,dElem,dQuadratureOrder,elemType,ierr)
   End Subroutine Element3D_Vect_InitSet_ByName

#undef __FUNCT__
#define __FUNCT__ "Element2D_Scal_InitSet_ByShortID"
   Subroutine Element2D_Scal_InitSet_ByShortID(mesh,cellIS,dElem,dQuadratureOrder,shortID,ierr)
      Type(tDM),intent(IN)                             :: mesh
      Type(tIS),intent(IN)                             :: cellIS
      Type(MEF90Element2D_Scal),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                              :: dQuadratureOrder,shortID
      PetscErrorCode,Intent(OUT)                       :: ierr
      
      Call Element2D_Scal_InitSet(mesh,cellIS,dElem,dQuadratureOrder,MEF90knownElements(ShortID),ierr)
   End Subroutine Element2D_Scal_InitSet_ByShortID

#undef __FUNCT__
#define __FUNCT__ "Element2D_Vect_InitSet_ByShortID"
   Subroutine Element2D_Vect_InitSet_ByShortID(mesh,cellIS,dElem,dQuadratureOrder,shortID,ierr)
      Type(tDM),intent(IN)                             :: mesh
      Type(tIS),intent(IN)                             :: cellIS
      Type(MEF90Element2D_Vect),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                              :: dQuadratureOrder,shortID
      PetscErrorCode,Intent(OUT)                       :: ierr
      
      Call Element2D_Vect_InitSet(mesh,cellIS,dElem,dQuadratureOrder,MEF90knownElements(ShortID),ierr)
   End Subroutine Element2D_Vect_InitSet_ByShortID

#undef __FUNCT__
#define __FUNCT__ "Element3D_Scal_InitSet_ByShortID"
   Subroutine Element3D_Scal_InitSet_ByShortID(mesh,cellIS,dElem,dQuadratureOrder,shortID,ierr)
      Type(tDM),intent(IN)                             :: mesh
      Type(tIS),intent(IN)                             :: cellIS
      Type(MEF90Element3D_Scal),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                              :: dQuadratureOrder,shortID
      PetscErrorCode,Intent(OUT)                       :: ierr
      
      Call Element3D_Scal_InitSet(mesh,cellIS,dElem,dQuadratureOrder,MEF90knownElements(ShortID),ierr)
   End Subroutine Element3D_Scal_InitSet_ByShortID

#undef __FUNCT__
#define __FUNCT__ "Element3D_Vect_InitSet_ByShortID"
   Subroutine Element3D_Vect_InitSet_ByShortID(mesh,cellIS,dElem,dQuadratureOrder,shortID,ierr)
      Type(tDM),intent(IN)                             :: mesh
      Type(tIS),intent(IN)                             :: cellIS
      Type(MEF90Element3D_Vect),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                              :: dQuadratureOrder,shortID
      PetscErrorCode,Intent(OUT)                       :: ierr
      
      Call Element3D_Vect_InitSet(mesh,cellIS,dElem,dQuadratureOrder,MEF90knownElements(ShortID),ierr)
   End Subroutine Element3D_Vect_InitSet_ByShortID

#undef __FUNCT__
#define __FUNCT__ "Element2D_Scal_InitSet"
   Subroutine Element2D_Scal_InitSet(dm,cellIS,dElem,dQuadratureOrder,elemType,ierr)
      Type(tDM),intent(IN)                             :: dm
      Type(tIS),intent(IN)                             :: cellIS
      Type(MEF90Element2D_Scal),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                              :: dQuadratureOrder
      Type(MEF90Element_Type),intent(IN)               :: elemType
      PetscErrorCode,Intent(OUT)                       :: ierr
      
      PetscInt                                         :: iELoc
      PetscInt,Dimension(:),Pointer                    :: CellID
      PetscReal, Dimension(:), Pointer                 :: v0,BB,BBinv
      PetscReal                                        :: detBBinv
      Type(Mat2D)                                      :: Bt
      PetscReal                                        :: length
      PetscReal,dimension(:),pointer                   :: centroid,innerNormal
      Type(Vect2D)                                     :: outerNormal

      Select Case (elemType%shortID)
         Case (MEF90_P1_Lagrange_2D%shortID,MEF90_P2_Lagrange_2D%shortID)
            PetscCall(ISGetIndicesF90(CellIS,CellID,ierr))
            Allocate(dElem(size(cellID)),stat=ierr)
            If (size(CellID) > 0) Then
               Allocate(v0(2))
               Allocate(BB(4))
               Allocate(BBinv(4))
               Do iELoc = 1,size(CellID)
                  PetscCall(DMPlexComputeCellGeometryAffineFEM(dm,cellID(iELoc),v0,BB,BBinv,detBBinv,ierr))
                  !!! Petsc uses a reference simplex with vertices at (-1,-1), (1,-1) and (-1,1)
                  !!! Whereas MEF90 uses (0,0), (1,0), (0,1), so we need to rescale the affine transform
                  Bt%XX = BBinv(1)*0.5_Kr; Bt%XY = BBinv(3)*0.5_Kr
                  Bt%YX = BBinv(2)*0.5_Kr; Bt%YY = BBinv(4)*0.5_Kr
                  detBBinv = 4.0_Kr * detBBinv
                  Call Element_P_Lagrange_2D_Scal_Init(dElem(iELoc),Bt,detBBinv,elemType%order,dQuadratureOrder,ierr)
               End Do
               DeAllocate(BBinv)
               DeAllocate(BB)
               DeAllocate(v0)
            End If
            PetscCall(ISRestoreIndicesF90(CellIS,CellID,ierr))
         Case (MEF90_P1_Lagrange_2DBoundary%shortID,MEF90_P2_Lagrange_2DBoundary%shortID)
            PetscCall(ISGetIndicesF90(CellIS,CellID,ierr))
            Allocate(dElem(size(cellID)),stat=ierr)
            If (size(CellID) > 0) Then
               allocate(centroid(2))
               allocate(innerNormal(2))
               Do iELoc = 1,size(CellID)
                  PetscCall(DMPlexComputeCellGeometryFVM(dm,cellID(iEloc),length,centroid,innerNormal,ierr))
                  outerNormal = -innerNormal
                  Call Element_P_Lagrange_2DBoundary_Scal_Init(dElem(iELoc),length,outerNormal,elemType%order,dQuadratureOrder,ierr)
               End Do
               DeAllocate(innerNormal)
               DeAllocate(centroid)
            End If
            PetscCall(ISRestoreIndicesF90(CellIS,CellID,ierr))
         !Case (MEF90_Q1_Lagrange_2D%shortID,MEF90_Q2_Lagrange_2D%shortID,MEF90_Q1_Lagrange_2DBoundary%shortID,MEF90_Q2_Lagrange_2DBoundary%shortID)
         !   !!! Get quadrature points for the current element using DMPlexComputeCellGeometryFEM
         !   !!! Initialize element
         Case Default
            Write(*,*) __FUNCT__,': Element type not implemented yet',elemType%name,elemType%shortID
            ierr = PETSC_ERR_SUP
      End Select
   End Subroutine Element2D_Scal_InitSet

#undef __FUNCT__
#define __FUNCT__ "Element2D_Vect_InitSet"
   Subroutine Element2D_Vect_InitSet(dm,cellIS,dElem,dQuadratureOrder,elemType,ierr)
      Type(tDM),intent(IN)                             :: dm
      Type(tIS),intent(IN)                             :: cellIS
      Type(MEF90Element2D_Vect),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                              :: dQuadratureOrder
      Type(MEF90Element_Type),intent(IN)               :: elemType
      PetscErrorCode,Intent(OUT)                       :: ierr
      
      PetscInt                                         :: iELoc
      PetscInt,Dimension(:),Pointer                    :: CellID
      PetscReal, Dimension(:), Pointer                 :: v0,BB,BBinv
      PetscReal                                        :: detBBinv
      Type(Mat2D)                                      :: Bt
      PetscReal                                        :: length
      PetscReal,dimension(:),pointer                   :: centroid,innerNormal
      Type(Vect2D)                                     :: outerNormal

      Select Case (elemType%shortID)
         Case (MEF90_P1_Lagrange_2D%shortID,MEF90_P2_Lagrange_2D%shortID)
            PetscCall(ISGetIndicesF90(CellIS,CellID,ierr))
            Allocate(dElem(size(cellID)),stat=ierr)
            If (size(CellID) > 0) Then
               Allocate(v0(2))
               Allocate(BB(4))
               Allocate(BBinv(4))
               Do_Elem_iE: Do iELoc = 1,size(CellID)
                  PetscCall(DMPlexComputeCellGeometryAffineFEM(dm,cellID(iELoc),v0,BB,BBinv,detBBinv,ierr))
                  !!! Petsc uses a reference simplex with vertices at (-1,-1), (1,-1) and (-1,1)
                  !!! Whereas MEF90 uses (0,0), (1,0), (0,1), so we need to rescale the affine transform
                  Bt%XX = BBinv(1)*0.5_Kr; Bt%XY = BBinv(3)*0.5_Kr
                  Bt%YX = BBinv(2)*0.5_Kr; Bt%YY = BBinv(4)*0.5_Kr
                  detBBinv = 4.0_Kr * detBBinv
                  Call Element_P_Lagrange_2D_Vect_Init(dElem(iELoc),Bt,detBBinv,elemType%order,dQuadratureOrder,ierr)
               End Do Do_Elem_iE
               DeAllocate(BBinv)
               DeAllocate(BB)
               DeAllocate(v0)
            End If
            PetscCall(ISRestoreIndicesF90(CellIS,CellID,ierr))
         Case (MEF90_P1_Lagrange_2DBoundary%shortID,MEF90_P2_Lagrange_2DBoundary%shortID)
            PetscCall(ISGetIndicesF90(CellIS,CellID,ierr))
            Allocate(dElem(size(cellID)),stat=ierr)
            If (size(CellID) > 0) Then
               allocate(centroid(2))
               allocate(innerNormal(2))
               Do iELoc = 1,size(CellID)
                  PetscCall(DMPlexComputeCellGeometryFVM(dm,cellID(iEloc),length,centroid,innerNormal,ierr))
                  outerNormal = -innerNormal
                  Call Element_P_Lagrange_2DBoundary_Vect_Init(dElem(iELoc),length,outerNormal,elemType%order,dQuadratureOrder,ierr)
               End Do
               DeAllocate(innerNormal)
               DeAllocate(centroid)
            End If
            PetscCall(ISRestoreIndicesF90(CellIS,CellID,ierr))
         !Case (MEF90_Q1_Lagrange_2D%shortID,MEF90_Q2_Lagrange_2D%shortID,MEF90_Q1_Lagrange_2DBoundary%shortID,MEF90_Q2_Lagrange_2DBoundary%shortID)
         !   !!! Get quadrature points for the current element using DMPlexComputeCellGeometryFEM
         !   !!! Initialize element
         Case Default
            Write(*,*) __FUNCT__,': Element type not implemented yet',elemType%name,elemType%shortID
            ierr = PETSC_ERR_SUP
      End Select
   End Subroutine Element2D_Vect_InitSet

#undef __FUNCT__
#define __FUNCT__ "Element3D_Scal_InitSet"
   Subroutine Element3D_Scal_InitSet(dm,cellIS,dElem,dQuadratureOrder,elemType,ierr)
      Type(tDM),intent(IN)                             :: dm
      Type(tIS),intent(IN)                             :: cellIS
      Type(MEF90Element3D_Scal),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                              :: dQuadratureOrder
      Type(MEF90Element_Type),intent(IN)               :: elemType
      PetscErrorCode,Intent(OUT)                       :: ierr
      
      PetscInt                                         :: iELoc
      PetscInt,Dimension(:),Pointer                    :: CellID
      PetscReal, Dimension(:), Pointer                 :: v0,BB,BBinv
      PetscReal                                        :: detBBinv
      Type(Mat3D)                                      :: Bt
      PetscReal                                        :: area
      PetscReal,dimension(:),pointer                   :: centroid,innerNormal
      Type(Vect3D)                                     :: outerNormal

      Select Case (elemType%shortID)
         Case (MEF90_P1_Lagrange_3D%shortID,MEF90_P2_Lagrange_3D%shortID)
            PetscCall(ISGetIndicesF90(CellIS,CellID,ierr))
            Allocate(dElem(size(cellID)),stat=ierr)
            If (size(CellID) > 0) Then
               Allocate(v0(3))
               Allocate(BB(9))
               Allocate(BBinv(9))
               Do_Elem_iE: Do iELoc = 1,size(CellID)
                  PetscCall(DMPlexComputeCellGeometryAffineFEM(dm,cellID(iELoc),v0,BB,BBinv,detBBinv,ierr))
                  !!! Petsc uses a reference simplex with vertices at (-1,-1), (1,-1) and (-1,1)
                  !!! Whereas MEF90 uses (0,0), (1,0), (0,1), so we need to rescale the affine transform
                  Bt%XX = BBinv(1)*0.5_Kr; Bt%XY = BBinv(4)*0.5_Kr; Bt%XZ = BBinv(7)*0.5_Kr
                  Bt%YX = BBinv(2)*0.5_Kr; Bt%YY = BBinv(5)*0.5_Kr; Bt%YZ = BBinv(8)*0.5_Kr
                  Bt%ZX = BBinv(3)*0.5_Kr; Bt%ZY = BBinv(6)*0.5_Kr; Bt%ZZ = BBinv(9)*0.5_Kr
                  detBBinv = 8.0_Kr * detBBinv
                  Call Element_P_Lagrange_3D_Scal_Init(dElem(iELoc),Bt,detBBinv,elemType%order,dQuadratureOrder,ierr)
               End Do Do_Elem_iE
               DeAllocate(BBinv)
               DeAllocate(BB)
               DeAllocate(v0)
            End If
            PetscCall(ISRestoreIndicesF90(CellIS,CellID,ierr))
         Case (MEF90_P1_Lagrange_3DBoundary%shortID,MEF90_P2_Lagrange_3DBoundary%shortID)
            PetscCall(ISGetIndicesF90(CellIS,CellID,ierr))
            Allocate(dElem(size(cellID)),stat=ierr)
            If (size(CellID) > 0) Then
               allocate(centroid(3))
               allocate(innerNormal(3))
               Do iELoc = 1,size(CellID)
                  PetscCall(DMPlexComputeCellGeometryFVM(dm,cellID(iEloc),area,centroid,innerNormal,ierr))
                  outerNormal = -innerNormal
                  Call Element_P_Lagrange_3DBoundary_Scal_Init(dElem(iELoc),area,outerNormal,elemType%order,dQuadratureOrder,ierr)
               End Do
               DeAllocate(innerNormal)
               DeAllocate(centroid)
            End If
            PetscCall(ISRestoreIndicesF90(CellIS,CellID,ierr))
         !Case (MEF90_Q1_Lagrange_3D%shortID,MEF90_Q2_Lagrange_3D%shortID,MEF90_Q1_Lagrange_3DBoundary%shortID,MEF90_Q2_Lagrange_3DBoundary%shortID)
         !   !!! Get quadrature points for the current element using DMPlexComputeCellGeometryFEM
         !   !!! Initialize element
         Case Default
            Write(*,*) __FUNCT__,': Element type not implemented yet',elemType%name,elemType%shortID
            ierr = PETSC_ERR_SUP
      End Select
   End Subroutine Element3D_Scal_InitSet

#undef __FUNCT__
#define __FUNCT__ "Element3D_Vect_InitSet"
   Subroutine Element3D_Vect_InitSet(dm,cellIS,dElem,dQuadratureOrder,elemType,ierr)
      Type(tDM),intent(IN)                             :: dm
      Type(tIS),intent(IN)                             :: cellIS
      Type(MEF90Element3D_Vect),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                              :: dQuadratureOrder
      Type(MEF90Element_Type),intent(IN)               :: elemType
      PetscErrorCode,Intent(OUT)                       :: ierr
      
      PetscInt                                         :: iELoc
      PetscInt,Dimension(:),Pointer                    :: CellID
      PetscReal, Dimension(:), Pointer                 :: v0,BB,BBinv
      PetscReal                                        :: detBBinv
      Type(Mat3D)                                      :: Bt
      PetscReal                                        :: area
      PetscReal,dimension(:),pointer                   :: centroid,innerNormal
      Type(Vect3D)                                     :: outerNormal

      Select Case (elemType%shortID)
         Case (MEF90_P1_Lagrange_3D%shortID,MEF90_P2_Lagrange_3D%shortID)
            PetscCall(ISGetIndicesF90(CellIS,CellID,ierr))
            Allocate(dElem(size(cellID)),stat=ierr)
            If (size(CellID) > 0) Then
               Allocate(v0(3))
               Allocate(BB(9))
               Allocate(BBinv(9))
               Do_Elem_iE: Do iELoc = 1,size(CellID)
                  PetscCall(DMPlexComputeCellGeometryAffineFEM(dm,cellID(iELoc),v0,BB,BBinv,detBBinv,ierr))
                  !!! Petsc uses a reference simplex with vertices at (-1,-1), (1,-1) and (-1,1)
                  !!! Whereas MEF90 uses (0,0), (1,0), (0,1), so we need to rescale the affine transform
                  Bt%XX = BBinv(1)*0.5_Kr; Bt%XY = BBinv(4)*0.5_Kr; Bt%XZ = BBinv(7)*0.5_Kr
                  Bt%YX = BBinv(2)*0.5_Kr; Bt%YY = BBinv(5)*0.5_Kr; Bt%YZ = BBinv(8)*0.5_Kr
                  Bt%ZX = BBinv(3)*0.5_Kr; Bt%ZY = BBinv(6)*0.5_Kr; Bt%ZZ = BBinv(9)*0.5_Kr
                  detBBinv = 8.0_Kr * detBBinv
                  Call Element_P_Lagrange_3D_Vect_Init(dElem(iELoc),Bt,detBBinv,elemType%order,dQuadratureOrder,ierr)
               End Do Do_Elem_iE
               DeAllocate(BBinv)
               DeAllocate(BB)
               DeAllocate(v0)
            End If
            PetscCall(ISRestoreIndicesF90(CellIS,CellID,ierr))
         Case (MEF90_P1_Lagrange_3DBoundary%shortID,MEF90_P2_Lagrange_3DBoundary%shortID)
            PetscCall(ISGetIndicesF90(CellIS,CellID,ierr))
            Allocate(dElem(size(cellID)),stat=ierr)
            If (size(CellID) > 0) Then
               allocate(centroid(3))
               allocate(innerNormal(3))
               Do iELoc = 1,size(CellID)
                  PetscCall(DMPlexComputeCellGeometryFVM(dm,cellID(iEloc),area,centroid,innerNormal,ierr))
                  outerNormal = -innerNormal
                  Call Element_P_Lagrange_3DBoundary_Vect_Init(dElem(iELoc),area,outerNormal,elemType%order,dQuadratureOrder,ierr)
               End Do
               DeAllocate(innerNormal)
               DeAllocate(centroid)
            End If
            PetscCall(ISRestoreIndicesF90(CellIS,CellID,ierr))
         !Case (MEF90_Q1_Lagrange_3D%shortID,MEF90_Q2_Lagrange_3D%shortID,MEF90_Q1_Lagrange_3DBoundary%shortID,MEF90_Q2_Lagrange_3DBoundary%shortID)
         !   !!! Get quadrature points for the current element using DMPlexComputeCellGeometryFEM
         !   !!! Initialize element
         Case Default
            Write(*,*) __FUNCT__,': Element type not implemented yet',elemType%name,elemType%shortID
            ierr = PETSC_ERR_SUP
      End Select
   End Subroutine Element3D_Vect_InitSet

#undef __FUNCT__
#define __FUNCT__ "Element_P_Lagrange_2D_Scal_Init"
   Subroutine Element_P_Lagrange_2D_Scal_Init(dElem,Bt,DetBinv,dPolynomialOrder,dQuadratureOrder,ierr)
      ! Compute the quadrature weights and the value of the basis functions and their gradient 
      ! at the quadrature points.
      ! Quadrature rules courtesy of Shawn Walker walker@math.lsu.edu
      Type(MEF90Element2D_Scal)              :: dElem
      Type(Mat2D),Intent(IN)                 :: Bt          ! The transposed of transformation matrix
      PetscReal,Intent(IN)                   :: DetBinv     ! The determinant of B^{-1}
      PetscInt                               :: dPolynomialOrder,dQuadratureOrder
      PetscErrorCode,Intent(OUT)             :: ierr
      
      PetscInt                               :: Nb_Gauss
      PetscInt                               :: Num_Dof
      PetscInt                               :: iDoF,iG

      PetscReal,Dimension(:,:),Pointer       :: PhiHat      ! PhiHat(i,k) The value of the ith basis function at the kth integration point
      Type(Vect2D),Dimension(:,:),Pointer    :: GradPhiHat
      
      Type(Vect2D),Dimension(:),Pointer     :: Xi ! The quadrature points coordinates in the reference element


      Select Case (dQuadratureOrder)
      Case(0,1)
         Nb_Gauss = 1
         Allocate(Xi(Nb_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
         Xi(1)%X = 1.0_Kr / 3.0_Kr
         Xi(1)%Y = 1.0_Kr / 3.0_Kr
         dElem%Gauss_C = detBinv / 2.0_Kr
      Case(2)
         Nb_Gauss = 3
         Allocate(Xi(Nb_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
         Xi(1) = [ 1.0_Kr / 6.0_Kr,1.0_Kr / 6.0_Kr ]
         Xi(2) = [ 2.0_Kr / 3.0_Kr,1.0_Kr / 6.0_Kr ]
         Xi(3) = [ 1.0_Kr / 6.0_Kr,2.0_Kr / 3.0_Kr ]
         dElem%Gauss_C = detBinv / 6.0_Kr
      Case(3)
         Nb_Gauss = 4
         Allocate(Xi(Nb_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
         Xi(1) = [ 1.0_Kr / 3.0_Kr,1.0_Kr / 3.0_Kr ]
         Xi(2) = [ 3.0_Kr / 5.0_Kr,1.0_Kr / 5.0_Kr ]
         Xi(3) = [ 1.0_Kr / 5.0_Kr,3.0_Kr / 5.0_Kr ]
         Xi(4) = [ 1.0_Kr / 5.0_Kr,1.0_Kr / 5.0_Kr ]
         dElem%Gauss_C(1)   = -detBinv * 9.0_Kr / 32.0_Kr
         dElem%Gauss_C(2:4) =  detBinv * 25.0_Kr / 96.0_Kr
      Case(4)
         Nb_Gauss = 6
         Allocate(Xi(Nb_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
         Xi(1) = [ 0.816847572980459_Kr,  0.091576213509771_Kr ]
         Xi(2) = [ 0.091576213509771_Kr,  0.816847572980459_Kr ]
         Xi(3) = [ 0.091576213509771_Kr,  0.091576213509771_Kr ]
         Xi(4) = [ 0.108103018168070_Kr,  0.445948490915965_Kr ]
         Xi(5) = [ 0.445948490915965_Kr,  0.108103018168070_Kr ]
         Xi(6) = [ 0.445948490915965_Kr,  0.445948490915965_Kr ]
         dElem%Gauss_C(1:3) = 0.109951743655322 / 2.0_Kr * detBinv
         dElem%Gauss_C(4:6) = 0.223381589678011 / 2.0_Kr * detBinv
      Case(5)
         Nb_Gauss = 7
         Allocate(Xi(Nb_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
         Xi(1) = [ 0.3333333333333333_Kr, 0.3333333333333333_Kr] 
         Xi(2) = [ 0.4701420641051150_Kr, 0.4701420641051150_Kr] 
         Xi(3) = [ 0.0597158717897698_Kr, 0.4701420641051150_Kr] 
         Xi(4) = [ 0.4701420641051150_Kr, 0.0597158717897698_Kr] 
         Xi(5) = [ 0.1012865073234563_Kr, 0.1012865073234563_Kr] 
         Xi(6) = [ 0.7974269853530872_Kr, 0.1012865073234563_Kr] 
         Xi(7) = [ 0.1012865073234563_Kr, 0.7974269853530872_Kr] 
         dElem%Gauss_C(1) = .225_Kr / 2.0_Kr * detBinv
         dElem%Gauss_C(2:4) = 0.13239415278850616 / 2.0_Kr * detBinv
         dElem%Gauss_C(5:7) = 0.12593918054482717 / 2.0_Kr * detBinv
      Case(-6)
         !!! It seems to me that this is really a quadrature rule of order 5...
         Nb_Gauss = 9
         Allocate(Xi(Nb_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
         Xi(1) = [ 0.124949503233232_Kr,  0.437525248383384_Kr ]
         Xi(2) = [ 0.437525248383384_Kr,  0.124949503233232_Kr ]
         Xi(3) = [ 0.437525248383384_Kr,  0.437525248383384_Kr ]
         Xi(4) = [ 0.797112651860071_Kr,  0.165409927389841_Kr ]
         Xi(5) = [ 0.797112651860071_Kr,  0.037477420750088_Kr ]
         Xi(6) = [ 0.165409927389841_Kr,  0.797112651860071_Kr ]
         Xi(7) = [ 0.165409927389841_Kr,  0.037477420750088_Kr ]
         Xi(8) = [ 0.037477420750088_Kr,  0.797112651860071_Kr ]
         Xi(9) = [ 0.037477420750088_Kr,  0.165409927389841_Kr ]
         dElem%Gauss_C(1:3) = 0.205950504760887_Kr / 2.0_Kr * detBinv
         dElem%Gauss_C(4:9) = 0.063691414286223_Kr / 2.0_Kr * detBinv
      Case(6)
         Nb_Gauss = 12
         Allocate(Xi(Nb_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
         Xi(1)  = [ 0.873821971016996_Kr,  0.063089014491502_Kr ]
         Xi(2)  = [ 0.063089014491502_Kr,  0.873821971016996_Kr ]
         Xi(3)  = [ 0.063089014491502_Kr,  0.063089014491502_Kr ]
         Xi(4)  = [ 0.501426509658179_Kr,  0.249286745170910_Kr ]
         Xi(5)  = [ 0.249286745170910_Kr,  0.501426509658179_Kr ]
         Xi(6)  = [ 0.249286745170910_Kr,  0.249286745170910_Kr ]
         Xi(7)  = [ 0.636502499121399_Kr,  0.310352451033785_Kr ]
         Xi(8)  = [ 0.636502499121399_Kr,  0.053145049844816_Kr ]
         Xi(9)  = [ 0.310352451033785_Kr,  0.636502499121399_Kr ]
         Xi(10) = [ 0.310352451033785_Kr,  0.053145049844816_Kr ]
         Xi(11) = [ 0.053145049844816_Kr,  0.636502499121399_Kr ]
         Xi(12) = [ 0.053145049844816_Kr,  0.310352451033785_Kr ]
         dElem%Gauss_C(1:3)  = 0.050844906370207_Kr / 2.0_Kr * detBinv
         dElem%Gauss_C(4:6)  = 0.116786275726379_Kr / 2.0_Kr * detBinv
         dElem%Gauss_C(7:12) = 0.082851075618374_Kr / 2.0_Kr * detBinv
      Case(7)
         Nb_Gauss = 13
         Allocate(Xi(Nb_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
         Xi(1)  = [ 0.333333333333333_Kr, 0.333333333333333_Kr ]
         Xi(2)  = [ 0.479308067841923_Kr, 0.260345966079038_Kr ]
         Xi(3)  = [ 0.260345966079038_Kr, 0.479308067841923_Kr ]
         Xi(4)  = [ 0.260345966079038_Kr, 0.260345966079038_Kr ]
         Xi(5)  = [ 0.869739794195568_Kr, 0.065130102902216_Kr ]
         Xi(6)  = [ 0.065130102902216_Kr, 0.869739794195568_Kr ]
         Xi(7)  = [ 0.065130102902216_Kr, 0.065130102902216_Kr ]
         Xi(8)  = [ 0.638444188569809_Kr, 0.312865496004875_Kr ]
         Xi(9)  = [ 0.638444188569809_Kr, 0.048690315425316_Kr ]
         Xi(10) = [ 0.312865496004875_Kr, 0.638444188569809_Kr ]
         Xi(11) = [ 0.312865496004875_Kr, 0.048690315425316_Kr ]
         Xi(12) = [ 0.048690315425316_Kr, 0.638444188569809_Kr ]
         Xi(13) = [ 0.048690315425316_Kr, 0.312865496004875_Kr ]
         dElem%Gauss_C(1)    = -0.149570044467670_Kr / 2.0_Kr * detBinv
         dElem%Gauss_C(2:4)  =  0.175615257433204_Kr / 2.0_Kr * detBinv
         dElem%Gauss_C(5:7)  =  0.053347235608839_Kr / 2.0_Kr * detBinv
         dElem%Gauss_C(8:13) =  0.077113760890257_Kr / 2.0_Kr * detBinv
      Case(8,9)
         Nb_Gauss = 19
         Allocate(Xi(Nb_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
         Xi(1)  = [ 0.3333333333333333_Kr,  0.3333333333333333_Kr ]
         Xi(2)  = [ 0.0206349616025259_Kr,  0.4896825191987370_Kr ]
         Xi(3)  = [ 0.4896825191987370_Kr,  0.0206349616025259_Kr ]
         Xi(4)  = [ 0.4896825191987370_Kr,  0.4896825191987370_Kr ]
         Xi(5)  = [ 0.1258208170141290_Kr,  0.4370895914929355_Kr ]
         Xi(6)  = [ 0.4370895914929355_Kr,  0.1258208170141290_Kr ]
         Xi(7)  = [ 0.4370895914929355_Kr,  0.4370895914929355_Kr ]
         Xi(8)  = [ 0.6235929287619356_Kr,  0.1882035356190322_Kr ]
         Xi(9)  = [ 0.1882035356190322_Kr,  0.6235929287619356_Kr ]
         Xi(10) = [ 0.1882035356190322_Kr,  0.1882035356190322_Kr ]
         Xi(11) = [ 0.9105409732110941_Kr,  0.0447295133944530_Kr ]
         Xi(12) = [ 0.0447295133944530_Kr,  0.9105409732110941_Kr ]
         Xi(13) = [ 0.0447295133944530_Kr,  0.0447295133944530_Kr ]
         Xi(14) = [ 0.7411985987844980_Kr,  0.0368384120547363_Kr ]
         Xi(15) = [ 0.7411985987844980_Kr,  0.2219628891607657_Kr ]
         Xi(16) = [ 0.0368384120547363_Kr,  0.7411985987844980_Kr ]
         Xi(17) = [ 0.0368384120547363_Kr,  0.2219628891607657_Kr ]
         Xi(18) = [ 0.2219628891607657_Kr,  0.7411985987844980_Kr ]
         Xi(19) = [ 0.2219628891607657_Kr,  0.0368384120547363_Kr ]
         dElem%Gauss_C(1)     = 0.09713579628279610_Kr / 2.0_Kr * detBinv
         dElem%Gauss_C(2:4)   = 0.03133470022713983_Kr / 2.0_Kr * detBinv
         dElem%Gauss_C(5:7)   = 0.07782754100477543_Kr / 2.0_Kr * detBinv
         dElem%Gauss_C(8:10)  = 0.07964773892720910_Kr / 2.0_Kr * detBinv
         dElem%Gauss_C(11:13) = 0.02557767565869810_Kr / 2.0_Kr * detBinv
         dElem%Gauss_C(14:19) = 0.04328353937728940_Kr / 2.0_Kr * detBinv
      Case Default
         Write(*,*) __FUNCT__,': Unimplemented quadrature order',dQuadratureOrder
         STOP
      End Select
      
      Select Case (dPolynomialOrder)
      Case(1)
         Num_DoF = 3
         Allocate(PhiHat(Num_DoF,Nb_Gauss),stat=ierr)
         Allocate(GradPhiHat(Num_DoF,Nb_Gauss),stat=ierr)
         PhiHat(1,:) = 1.0_Kr - Xi%X - Xi%Y
         PhiHat(2,:) = Xi(:)%X
         PhiHat(3,:) = Xi(:)%Y
         
         GradPhiHat(1,:)%X = -1.0_Kr;GradPhiHat(1,:)%Y = -1.0_Kr 
         GradPhiHat(2,:)%X =  1.0_Kr;GradPhiHat(2,:)%Y =  0.0_Kr 
         GradPhiHat(3,:)%X =  0.0_Kr;GradPhiHat(3,:)%Y =  1.0_Kr
          
      Case(2)
         Num_DoF = 6
         Allocate(PhiHat(Num_DoF,Nb_Gauss),stat=ierr)
         Allocate(GradPhiHat(Num_DoF,Nb_Gauss),stat=ierr)
         PhiHat(1,:) = (1.0_Kr - Xi%X - Xi%Y) * (1.0_Kr - 2.0_Kr * Xi%X - 2.0_Kr * Xi%Y)      
         PhiHat(2,:) = Xi%X * (2.0_Kr * Xi%X - 1.0_Kr)
         PhiHat(3,:) = Xi%Y * (2.0_Kr * Xi%Y - 1.0_Kr)
         PhiHat(4,:) = 4.0_Kr * Xi%X * (1.0_Kr - Xi%X - Xi%Y)
         PhiHat(5,:) = 4.0_Kr * Xi%X * Xi%Y
         PhiHat(6,:) = 4.0_Kr * Xi%Y * (1.0_Kr - Xi%X - Xi%Y)
         
         GradPhiHat(1,:)%X = 4.0_Kr * Xi%X + 4.0_Kr * Xi%y -3.0_Kr;     GradPhiHat(1,:)%Y = 4.0_Kr * Xi%X + 4.0_Kr * Xi%y -3.0_Kr
         GradPhiHat(2,:)%X = 4.0_Kr * Xi%X - 1.0_Kr;                    GradPhiHat(2,:)%Y = 0.0_Kr
         GradPhiHat(3,:)%X = 0.0_Kr;                                    GradPhiHat(3,:)%Y = 4.0_Kr * Xi%Y - 1.0_Kr
         GradPhiHat(4,:)%X = 4.0_Kr * (1.0_Kr - 2.0_Kr * Xi%X - Xi%Y);  GradPhiHat(4,:)%Y = -4.0_Kr * Xi%X
         GradPhiHat(5,:)%X = 4.0_Kr * Xi%Y;                             GradPhiHat(5,:)%Y = 4.0_Kr * Xi%X
         GradPhiHat(6,:)%X = -4.0_Kr * Xi%Y;                            GradPhiHat(6,:)%Y = 4.0_Kr * (1.0_Kr - Xi%X - 2.0_Kr * Xi%Y);
      Case Default
         Num_DoF = 0
         Write(*,*) __FUNCT__,': Unimplemented PolynomialOrder',dPolynomialOrder
         STOP
      End Select
      
      Allocate (dElem%BF(Num_DoF,Nb_Gauss),stat=ierr) 
      Allocate (dElem%Grad_BF(Num_DoF,Nb_Gauss),stat=ierr)
      dElem%BF = PhiHat
      Do iDoF = 1,Num_DoF
         Do iG = 1,Nb_Gauss
            dElem%Grad_BF(iDoF,iG) = Bt * GradPhiHat(iDoF,iG) 
         End Do
      End Do
      
      dElem%outerNormal = 0.0_Kr
     
      DeAllocate(Xi,stat=ierr)
      DeAllocate(PhiHat,stat=ierr)
      DeAllocate(GradPhiHat,stat=ierr)
   End Subroutine Element_P_Lagrange_2D_Scal_Init
   
#undef __FUNCT__
#define __FUNCT__ "Element_P_Lagrange_2DBoundary_Scal_Init"
   Subroutine Element_P_Lagrange_2DBoundary_Scal_Init(dElem,l,outerNormal,dPolynomialOrder,dQuadratureOrder,ierr)
      Type(MEF90Element2D_Scal),intent(INOUT):: dElem
      PetscReal,Intent(IN)                   :: l
      Type(Vect2D),Intent(IN)                :: outerNormal
      PetscInt,Intent(IN)                    :: dPolynomialOrder,dQuadratureOrder
      PetscErrorCode,Intent(OUT)             :: ierr
      
      PetscReal, Dimension(:), Pointer       :: Xi
      PetscInt                               :: iDoF,iG,Num_Gauss,Num_DoF

      num_Dof = 0
      num_Gauss = 0
      dElem%outerNormal = outerNormal
      Select Case(dQuadratureOrder)
      Case(0,1)
         Num_Gauss = 1
         Allocate(Xi(Num_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Num_Gauss),stat=ierr)
         Xi(1) = 0.0_Kr
         dElem%Gauss_C(1) = l
      Case(2,3)
         Num_Gauss = 2
         Allocate(Xi(Num_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Num_Gauss),stat=ierr)
         Xi(1) = -1.0_Kr / sqrt(3.0_Kr)
         Xi(2) =  1.0_Kr / sqrt(3.0_Kr)
         dElem%Gauss_C = l * .5_Kr
      Case(4,5)
         Num_Gauss = 3
         Allocate(Xi(Num_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Num_Gauss),stat=ierr)
         Xi(1) = -sqrt(15.0_Kr) / 5.0_Kr
         Xi(2) = 0.0_Kr
         Xi(3) = sqrt(15.0_Kr) / 5.0_Kr
         dElem%Gauss_C(1) = l * 5.0_Kr / 18.0_Kr
         dElem%Gauss_C(2) = l * 4.0_Kr / 9.0_Kr
         dElem%Gauss_C(3) = l * 5.0_Kr / 18.0_Kr
       Case(6,7)
          Num_Gauss = 4
          Allocate(Xi(Num_Gauss),stat=ierr)
          Allocate(dElem%Gauss_C(Num_Gauss),stat=ierr)
          Xi(1) = -sqrt((15.0_Kr + 2.0_Kr * sqrt(30.0_Kr))/35.0_Kr)
          Xi(2) = -sqrt((15.0_Kr - 2.0_Kr * sqrt(30.0_Kr))/35.0_Kr)
          Xi(3) =  sqrt((15.0_Kr - 2.0_Kr * sqrt(30.0_Kr))/35.0_Kr)
          Xi(4) =  sqrt((15.0_Kr + 2.0_Kr * sqrt(30.0_Kr))/35.0_Kr)
          dElem%Gauss_C(1) = l * (.25_kr - sqrt(5.0_Kr / 864.0_Kr)) 
          dElem%Gauss_C(2) = l * (.25_kr + sqrt(5.0_Kr / 864.0_Kr))
          dElem%Gauss_C(3) = l * (.25_kr + sqrt(5.0_Kr / 864.0_Kr))
          dElem%Gauss_C(4) = l * (.25_kr - sqrt(5.0_Kr / 864.0_Kr))
      Case(8,9,10,11,12)
         Num_Gauss = 5
         Allocate(Xi(Num_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Num_Gauss),stat=ierr)
         Xi(1) = -sqrt(5.0_Kr + sqrt(40.0_Kr / 7.0_Kr)) / 3.0_Kr
         Xi(2) = -sqrt(5.0_Kr - sqrt(40.0_Kr / 7.0_Kr)) / 3.0_Kr
         Xi(3) =  0.0_Kr
         Xi(4) =  sqrt(5.0_Kr - sqrt(40.0_Kr / 7.0_Kr)) / 3.0_Kr
         Xi(5) =  sqrt(5.0_Kr + sqrt(40.0_Kr / 7.0_Kr)) / 3.0_Kr
         dElem%Gauss_C(1) = l * (322.0_Kr - sqrt(11830.0_Kr)) / 1800.0_Kr
         dElem%Gauss_C(2) = l * (322.0_Kr + sqrt(11830.0_Kr)) / 1800.0_Kr
         dElem%Gauss_C(3) = l * 128.0_Kr / 450.0_Kr
         dElem%Gauss_C(4) = l * (322.0_Kr + sqrt(11830.0_Kr)) / 1800.0_Kr
         dElem%Gauss_C(5) = l * (322.0_Kr - sqrt(11830.0_Kr)) / 1800.0_Kr
      Case Default
         Write(*,*) __FUNCT__,': Unimplemented quadrature order',dQuadratureOrder
         STOP
      End Select         
      Select Case (dPolynomialOrder)
         Case (1)
            Num_DoF  = 2
            Allocate(dElem%BF(Num_DoF,Num_Gauss),stat=ierr)
            dElem%BF(1,:) = (1.0_Kr - Xi) * .5_Kr
            dElem%BF(2,:) = (1.0_Kr + Xi) * .5_Kr
         Case (2)
            Num_DoF  = 3
            Allocate(dElem%BF(Num_DoF,Num_Gauss),stat=ierr)
            dElem%BF(1,:) = Xi * (Xi - 1.0_Kr) * .5_Kr
            dElem%BF(2,:) = Xi * (Xi + 1.0_Kr) * .5_Kr
            dElem%BF(3,:) = 1.0_Kr - dElem%BF(1,:)  - dElem%BF(2,:) !(1.0_Kr - Xi) * (1.0_Kr + Xi)
         Case Default
            Write(*,*) '[ERROR]: Polynomial order ',dPolynomialOrder,' not implemented in ',__FUNCT__
      End Select
      Allocate(delem%Grad_BF(Num_DoF,Num_Gauss),stat=ierr)
      Do iDof = 1, num_dof
         Do iG = 1, num_Gauss
            delem%Grad_BF(iDof,iG) = 0.0_Kr
         End Do
      End Do
      DeAllocate(Xi,stat=ierr)
   End Subroutine Element_P_Lagrange_2DBoundary_Scal_Init                          
      
#undef __FUNCT__
#define __FUNCT__ "Element_P_Lagrange_2D_Vect_Init"
   Subroutine Element_P_Lagrange_2D_Vect_Init(dElem,Bt,DetBinv,dPolynomialOrder,dQuadratureOrder,ierr)
      Type(MEF90Element2D_Vect)              :: dElem
      Type(Mat2D),Intent(IN)                 :: Bt          ! The transposed of transformation matrix
      PetscReal,Intent(IN)                   :: DetBinv     ! The determinant of B^{-1}
      PetscInt,Intent(IN)                    :: dPolynomialOrder,dQuadratureOrder
      PetscErrorCode,Intent(OUT)             :: ierr
   
      Type(MEF90Element2D_Scal)              :: Elem_Scal
      PetscInt                               :: dim = 2 
      PetscInt                               :: Num_DoF,Nb_Gauss,i,iDof,iG
      
      
      Call Element_P_Lagrange_2D_Scal_Init(Elem_Scal,Bt,DetBinv,dPolynomialOrder,dQuadratureOrder,ierr)
      Num_DoF  = Size(Elem_Scal%BF,1) 
      Nb_Gauss = Size(Elem_Scal%BF,2)

      Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
      Allocate(dElem%BF(Num_DoF * dim,Nb_Gauss),stat=ierr)
      Allocate(dElem%Grad_BF(Num_DoF * dim,Nb_Gauss),stat=ierr)
      Allocate(dElem%GradS_BF(Num_DoF * dim,Nb_Gauss),stat=ierr)
         
      dElem%Gauss_C = Elem_Scal%Gauss_C
      Do iDof = 1, num_dof * dim
         Do iG = 1, Nb_Gauss
            dElem%BF(iDof,iG)       = 0.0_Kr
            delem%Grad_BF(iDof,iG)  = 0.0_Kr
            delem%GradS_BF(iDof,iG) = 0.0_Kr
         End Do
      End Do
      
      Do i = 0,Num_DoF-1
         dElem%BF(i*dim+1,:)%X = Elem_Scal%BF(i+1,:)
         dElem%BF(i*dim+2,:)%Y = Elem_Scal%BF(i+1,:)

         dElem%Grad_BF(i*dim+1,:)%XX = Elem_Scal%Grad_BF(i+1,:)%X
         dElem%Grad_BF(i*dim+1,:)%XY = Elem_Scal%Grad_BF(i+1,:)%Y
         dElem%Grad_BF(i*dim+2,:)%YX = Elem_Scal%Grad_BF(i+1,:)%X
         dElem%Grad_BF(i*dim+2,:)%YY = Elem_Scal%Grad_BF(i+1,:)%Y

         dElem%GradS_BF(i*dim+1,:)%XX = Elem_Scal%Grad_BF(i+1,:)%X
         dElem%GradS_BF(i*dim+1,:)%XY = Elem_Scal%Grad_BF(i+1,:)%Y / 2.0_Kr
         dElem%GradS_BF(i*dim+2,:)%XY = Elem_Scal%Grad_BF(i+1,:)%X / 2.0_Kr
         dElem%GradS_BF(i*dim+2,:)%YY = Elem_Scal%Grad_BF(i+1,:)%Y
      End Do
      
      Call MEF90Element_Destroy(Elem_Scal,ierr)
   End Subroutine Element_P_Lagrange_2D_Vect_Init

#undef __FUNCT__
#define __FUNCT__ "Element_P_Lagrange_2DBoundary_Vect_Init"
   Subroutine Element_P_Lagrange_2DBoundary_Vect_Init(dElem,length,outerNormal,dPolynomialOrder,dQuadratureOrder,ierr)
      Type(MEF90Element2D_Vect)              :: dElem
      PetscReal,Intent(IN)                   :: length
      Type(Vect2D),Intent(IN)                :: outerNormal 
      PetscInt,Intent(IN)                    :: dPolynomialOrder,dQuadratureOrder
      PetscErrorCode,Intent(OUT)             :: ierr
   
      Type(MEF90Element2D_Scal)              :: Elem_Scal
      PetscInt                               :: dim = 2 
      PetscInt                               :: Num_DoF,Nb_Gauss,iDof,iG
      
      
      Call Element_P_Lagrange_2DBoundary_Scal_Init(Elem_Scal,length,outerNormal,dPolynomialOrder,dQuadratureOrder,ierr)
      Num_DoF  = Size(Elem_Scal%BF,1) 
      Nb_Gauss = Size(Elem_Scal%BF,2)
      Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
      Allocate(dElem%BF(Num_DoF * dim,Nb_Gauss),stat=ierr)
      Allocate(dElem%Grad_BF(Num_DoF * dim,Nb_Gauss),stat=ierr)
      Allocate(dElem%GradS_BF(Num_DoF * dim,Nb_Gauss),stat=ierr)
         
      dElem%Gauss_C = Elem_Scal%Gauss_C
      Do iDof = 1, num_dof * dim
         Do iG = 1, Nb_Gauss
            dElem%BF(iDof,iG)       = 0.0_Kr
            delem%Grad_BF(iDof,iG)  = 0.0_Kr
            delem%GradS_BF(iDof,iG) = 0.0_Kr
         End Do
      End Do

      Do idof = 0,Num_DoF-1
         dElem%BF(iDof*dim+1,:)%X = Elem_Scal%BF(iDof+1,:)
         dElem%BF(iDof*dim+2,:)%Y = Elem_Scal%BF(iDof+1,:)
      End Do

      dElem%outerNormal = Elem_Scal%outerNormal
      Call MEF90Element_Destroy(Elem_Scal,ierr)
   End Subroutine Element_P_Lagrange_2DBoundary_Vect_Init


#undef __FUNCT__
#define __FUNCT__ "Element_P_Lagrange_3D_Scal_Init"
   Subroutine Element_P_Lagrange_3D_Scal_Init(dElem,Bt,DetBinv,dPolynomialOrder,dQuadratureOrder,ierr)
      ! Compute the quadrature weights and the value of the basis functions and their gradient 
      ! at the quadrature points.
      Type(MEF90Element3D_Scal)              :: dElem
      Type(Mat3D),Intent(IN)                 :: Bt          ! The transposed of transformation matrix
      PetscReal,Intent(IN)                   :: DetBinv     ! The determinant of B^{-1}
      PetscInt,Intent(IN)                    :: dPolynomialOrder,dQuadratureOrder
      PetscErrorCode,Intent(OUT)             :: ierr
      
      PetscInt                               :: Nb_Gauss
      PetscInt                               :: Num_Dof
      PetscInt                               :: iDoF,iG

      PetscReal,Dimension(:,:),Pointer       :: PhiHat      ! PhiHat(i,k) The value of the ith basis function at the kth integration point
      Type(Vect3D),Dimension(:,:),Pointer    :: GradPhiHat
      
      
      Type(Vect3D),Dimension(:),Pointer      :: Xi          ! The quadrature points coordinates in the reference element
      PetscReal                              :: a,b         ! Location of integration points in Aiken p. 272 table 10.4
      
      Select Case (dQuadratureOrder)
      Case (0,1)
         Nb_Gauss = 1
         Allocate(Xi(Nb_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
         Xi(1) = [ .25_Kr,.25_Kr,.25_Kr ]
         dElem%Gauss_C(1) = 1.0_Kr / 6.0_Kr * detBinv
      Case(2)
         a = (5.0_Kr + 3.0_Kr * sqrt(5.0_Kr)) / 20.0_Kr
         b = (5.0_Kr - sqrt(5.0_Kr)) / 20.0_Kr
         Nb_Gauss = 4
         Allocate(Xi(Nb_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
         Xi(1) = [ a,b,b ]
         Xi(2) = [ b,a,b ]
         Xi(3) = [ b,b,a ]
         Xi(4) = [ b,b,b ]
         dElem%Gauss_C(1:4) = 1.0_Kr / 24.0_Kr * detBinv
      Case(3)
         Nb_Gauss = 5
         Allocate(Xi(Nb_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
         Xi(1) = [ .25_Kr,.25_Kr,.25_Kr ]
         Xi(2) = [ .5_Kr,1._Kr / 6._Kr,1._Kr / 6._Kr ]
         Xi(3) = [ 1._Kr / 6._Kr,.5_Kr,1._Kr / 6._Kr ]
         Xi(4) = [ 1._Kr / 6._Kr,1._Kr / 6._Kr,.5_Kr ]
         Xi(5) = [ 1._Kr / 6._Kr,1._Kr / 6._Kr,1._Kr / 6._Kr ]
         dElem%Gauss_C(1) = -2.0_Kr / 15.0_Kr * detBinv
         dElem%Gauss_C(2:5) = 3.0_Kr / 40.0_Kr * detBinv
      Case(4)
         Nb_Gauss = 11
         Allocate(Xi(Nb_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
         Xi(1)  = [ 0.2500000000000000_Kr, 0.2500000000000000_Kr, 0.2500000000000000_Kr ]
         Xi(2)  = [ 0.7857142857142857_Kr, 0.0714285714285714_Kr, 0.0714285714285714_Kr ]
         Xi(3)  = [ 0.0714285714285714_Kr, 0.0714285714285714_Kr, 0.0714285714285714_Kr ]
         Xi(4)  = [ 0.0714285714285714_Kr, 0.0714285714285714_Kr, 0.7857142857142857_Kr ]
         Xi(5)  = [ 0.0714285714285714_Kr, 0.7857142857142857_Kr, 0.0714285714285714_Kr ]
         Xi(6)  = [ 0.1005964238332008_Kr, 0.3994035761667992_Kr, 0.3994035761667992_Kr ]
         Xi(7)  = [ 0.3994035761667992_Kr, 0.1005964238332008_Kr, 0.3994035761667992_Kr ]
         Xi(8)  = [ 0.3994035761667992_Kr, 0.3994035761667992_Kr, 0.1005964238332008_Kr ]
         Xi(9)  = [ 0.3994035761667992_Kr, 0.1005964238332008_Kr, 0.1005964238332008_Kr ]
         Xi(10) = [ 0.1005964238332008_Kr, 0.3994035761667992_Kr, 0.1005964238332008_Kr ]
         Xi(11) = [ 0.1005964238332008_Kr, 0.1005964238332008_Kr, 0.3994035761667992_Kr ]
         dElem%Gauss_C(1)    = -0.0789333333333333_Kr / 6.0_Kr * detBinv
         dElem%Gauss_C(2:5)  =  0.0457333333333333_Kr / 6.0_Kr * detBinv
         dElem%Gauss_C(6:11) =  0.1493333333333333_Kr / 6.0_Kr * detBinv
      Case(5)
         Nb_Gauss = 15
         Allocate(Xi(Nb_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
         Xi(1)  = [0.2500000000000000_Kr,  0.2500000000000000_Kr,  0.2500000000000000_Kr]
         Xi(2)  = [0.0000000000000000_Kr,  0.3333333333333333_Kr,  0.3333333333333333_Kr]
         Xi(3)  = [0.3333333333333333_Kr,  0.3333333333333333_Kr,  0.3333333333333333_Kr]
         Xi(4)  = [0.3333333333333333_Kr,  0.3333333333333333_Kr,  0.0000000000000000_Kr]
         Xi(5)  = [0.3333333333333333_Kr,  0.0000000000000000_Kr,  0.3333333333333333_Kr]
         Xi(6)  = [0.7272727272727273_Kr,  0.0909090909090909_Kr,  0.0909090909090909_Kr]
         Xi(7)  = [0.0909090909090909_Kr,  0.0909090909090909_Kr,  0.0909090909090909_Kr]
         Xi(8)  = [0.0909090909090909_Kr,  0.0909090909090909_Kr,  0.7272727272727273_Kr]
         Xi(9)  = [0.0909090909090909_Kr,  0.7272727272727273_Kr,  0.0909090909090909_Kr]
         Xi(10) = [0.4334498464263357_Kr,  0.0665501535736643_Kr,  0.0665501535736643_Kr]
         Xi(11) = [0.0665501535736643_Kr,  0.4334498464263357_Kr,  0.0665501535736643_Kr]
         Xi(12) = [0.0665501535736643_Kr,  0.0665501535736643_Kr,  0.4334498464263357_Kr]
         Xi(13) = [0.0665501535736643_Kr,  0.4334498464263357_Kr,  0.4334498464263357_Kr]
         Xi(14) = [0.4334498464263357_Kr,  0.0665501535736643_Kr,  0.4334498464263357_Kr]
         Xi(15) = [0.4334498464263357_Kr,  0.4334498464263357_Kr,  0.0665501535736643_Kr]
         dElem%Gauss_C(1)  = 0.1817020685825351_Kr / 6.0_Kr * detBinv
         dElem%Gauss_C(2:5)  = 0.0361607142857143_Kr / 6.0_Kr * detBinv
         dElem%Gauss_C(6:9)  = 0.0698714945161738_Kr / 6.0_Kr * detBinv
         dElem%Gauss_C(10:15) = 0.0656948493683187_Kr / 6.0_Kr * detBinv
      Case(6)
         Nb_Gauss = 24
         Allocate(Xi(Nb_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
         Xi(1)  = [0.3561913862225449_Kr,  0.2146028712591517_Kr,  0.2146028712591517_Kr]
         Xi(2)  = [0.2146028712591517_Kr,  0.2146028712591517_Kr,  0.2146028712591517_Kr]
         Xi(3)  = [0.2146028712591517_Kr,  0.2146028712591517_Kr,  0.3561913862225449_Kr]
         Xi(4)  = [0.2146028712591517_Kr,  0.3561913862225449_Kr,  0.2146028712591517_Kr]
         Xi(5)  = [0.8779781243961660_Kr,  0.0406739585346113_Kr,  0.0406739585346113_Kr]
         Xi(6)  = [0.0406739585346113_Kr,  0.0406739585346113_Kr,  0.0406739585346113_Kr]
         Xi(7)  = [0.0406739585346113_Kr,  0.0406739585346113_Kr,  0.8779781243961660_Kr]
         Xi(8)  = [0.0406739585346113_Kr,  0.8779781243961660_Kr,  0.0406739585346113_Kr]
         Xi(9)  = [0.0329863295731731_Kr,  0.3223378901422757_Kr,  0.3223378901422757_Kr]
         Xi(10) = [0.3223378901422757_Kr,  0.3223378901422757_Kr,  0.3223378901422757_Kr]
         Xi(11) = [0.3223378901422757_Kr,  0.3223378901422757_Kr,  0.0329863295731731_Kr]
         Xi(12) = [0.3223378901422757_Kr,  0.0329863295731731_Kr,  0.3223378901422757_Kr]
         Xi(13) = [0.2696723314583159_Kr,  0.0636610018750175_Kr,  0.0636610018750175_Kr]
         Xi(14) = [0.0636610018750175_Kr,  0.2696723314583159_Kr,  0.0636610018750175_Kr]
         Xi(15) = [0.0636610018750175_Kr,  0.0636610018750175_Kr,  0.2696723314583159_Kr]
         Xi(16) = [0.6030056647916491_Kr,  0.0636610018750175_Kr,  0.0636610018750175_Kr]
         Xi(17) = [0.0636610018750175_Kr,  0.6030056647916491_Kr,  0.0636610018750175_Kr]
         Xi(18) = [0.0636610018750175_Kr,  0.0636610018750175_Kr,  0.6030056647916491_Kr]
         Xi(19) = [0.0636610018750175_Kr,  0.2696723314583159_Kr,  0.6030056647916491_Kr]
         Xi(20) = [0.2696723314583159_Kr,  0.6030056647916491_Kr,  0.0636610018750175_Kr]
         Xi(21) = [0.6030056647916491_Kr,  0.0636610018750175_Kr,  0.2696723314583159_Kr]
         Xi(22) = [0.0636610018750175_Kr,  0.6030056647916491_Kr,  0.2696723314583159_Kr]
         Xi(23) = [0.2696723314583159_Kr,  0.0636610018750175_Kr,  0.6030056647916491_Kr]
         Xi(24) = [0.6030056647916491_Kr,  0.2696723314583159_Kr,  0.0636610018750175_Kr]
         dElem%Gauss_C(1:4)   = 0.0399227502581679_Kr / 6.0_Kr * detBinv
         dElem%Gauss_C(5:8)   = 0.0100772110553207_Kr / 6.0_Kr * detBinv
         dElem%Gauss_C(9:12)  = 0.0553571815436544_Kr / 6.0_Kr * detBinv
         dElem%Gauss_C(13:24) = 0.0482142857142857_Kr / 6.0_Kr * detBinv
      Case(7)
         Nb_Gauss = 31
         Allocate(Xi(Nb_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
         Xi(1)  = [ 0.2500000000000000_Kr, 0.2500000000000000_Kr, 0.2500000000000000_Kr ]
         Xi(2)  = [ 0.7653604230090441_Kr, 0.0782131923303186_Kr, 0.0782131923303186_Kr ]
         Xi(3)  = [ 0.0782131923303186_Kr, 0.0782131923303186_Kr, 0.0782131923303186_Kr ]
         Xi(4)  = [ 0.0782131923303186_Kr, 0.0782131923303186_Kr, 0.7653604230090441_Kr ]
         Xi(5)  = [ 0.0782131923303186_Kr, 0.7653604230090441_Kr, 0.0782131923303186_Kr ]
         Xi(6)  = [ 0.6344703500082868_Kr, 0.1218432166639044_Kr, 0.1218432166639044_Kr ]
         Xi(7)  = [ 0.1218432166639044_Kr, 0.1218432166639044_Kr, 0.1218432166639044_Kr ]
         Xi(8)  = [ 0.1218432166639044_Kr, 0.1218432166639044_Kr, 0.6344703500082868_Kr ]
         Xi(9)  = [ 0.1218432166639044_Kr, 0.6344703500082868_Kr, 0.1218432166639044_Kr ]
         Xi(10) = [ 0.0023825066607383_Kr, 0.3325391644464206_Kr, 0.3325391644464206_Kr ]
         Xi(11) = [ 0.3325391644464206_Kr, 0.3325391644464206_Kr, 0.3325391644464206_Kr ]
         Xi(12) = [ 0.3325391644464206_Kr, 0.3325391644464206_Kr, 0.0023825066607383_Kr ]
         Xi(13) = [ 0.3325391644464206_Kr, 0.0023825066607383_Kr, 0.3325391644464206_Kr ]
         Xi(14) = [ 0.0000000000000000_Kr, 0.5000000000000000_Kr, 0.5000000000000000_Kr ]
         Xi(15) = [ 0.5000000000000000_Kr, 0.0000000000000000_Kr, 0.5000000000000000_Kr ]
         Xi(16) = [ 0.5000000000000000_Kr, 0.5000000000000000_Kr, 0.0000000000000000_Kr ]
         Xi(17) = [ 0.5000000000000000_Kr, 0.0000000000000000_Kr, 0.0000000000000000_Kr ]
         Xi(18) = [ 0.0000000000000000_Kr, 0.5000000000000000_Kr, 0.0000000000000000_Kr ]
         Xi(19) = [ 0.0000000000000000_Kr, 0.0000000000000000_Kr, 0.5000000000000000_Kr ]
         Xi(20) = [ 0.2000000000000000_Kr, 0.1000000000000000_Kr, 0.1000000000000000_Kr ]
         Xi(21) = [ 0.1000000000000000_Kr, 0.2000000000000000_Kr, 0.1000000000000000_Kr ]
         Xi(22) = [ 0.1000000000000000_Kr, 0.1000000000000000_Kr, 0.2000000000000000_Kr ]
         Xi(23) = [ 0.6000000000000000_Kr, 0.1000000000000000_Kr, 0.1000000000000000_Kr ]
         Xi(24) = [ 0.1000000000000000_Kr, 0.6000000000000000_Kr, 0.1000000000000000_Kr ]
         Xi(25) = [ 0.1000000000000000_Kr, 0.1000000000000000_Kr, 0.6000000000000000_Kr ]
         Xi(26) = [ 0.1000000000000000_Kr, 0.2000000000000000_Kr, 0.6000000000000000_Kr ]
         Xi(27) = [ 0.2000000000000000_Kr, 0.6000000000000000_Kr, 0.1000000000000000_Kr ]
         Xi(28) = [ 0.6000000000000000_Kr, 0.1000000000000000_Kr, 0.2000000000000000_Kr ]
         Xi(29) = [ 0.1000000000000000_Kr, 0.6000000000000000_Kr, 0.2000000000000000_Kr ]
         Xi(30) = [ 0.2000000000000000_Kr, 0.1000000000000000_Kr, 0.6000000000000000_Kr ]
         Xi(31) = [ 0.6000000000000000_Kr, 0.2000000000000000_Kr, 0.1000000000000000_Kr ]
         dElem%Gauss_C(1)     =  0.1095853407966528_Kr / 6.0_Kr * detBinv
         dElem%Gauss_C(2:5)   =  0.0635996491464850_Kr / 6.0_Kr * detBinv
         dElem%Gauss_C(6:9)   = -0.3751064406859797_Kr / 6.0_Kr * detBinv
         dElem%Gauss_C(10:13) =  0.0293485515784412_Kr / 6.0_Kr * detBinv
         dElem%Gauss_C(14:19) =  0.0058201058201058_Kr / 6.0_Kr * detBinv
         dElem%Gauss_C(20:31) =  0.1653439153439105_Kr / 6.0_Kr * detBinv
      Case(8)
         Nb_Gauss = 45
         Allocate(Xi(Nb_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
         Xi(1)  = [ 0.2500000000000000_Kr, 0.2500000000000000_Kr, 0.2500000000000000_Kr ]
         Xi(2)  = [ 0.6175871903000830_Kr, 0.1274709365666390_Kr, 0.1274709365666390_Kr ]
         Xi(3)  = [ 0.1274709365666390_Kr, 0.1274709365666390_Kr, 0.1274709365666390_Kr ]
         Xi(4)  = [ 0.1274709365666390_Kr, 0.1274709365666390_Kr, 0.6175871903000830_Kr ]
         Xi(5)  = [ 0.1274709365666390_Kr, 0.6175871903000830_Kr, 0.1274709365666390_Kr ]
         Xi(6)  = [ 0.9037635088221031_Kr, 0.0320788303926323_Kr, 0.0320788303926323_Kr ]
         Xi(7)  = [ 0.0320788303926323_Kr, 0.0320788303926323_Kr, 0.0320788303926323_Kr ]
         Xi(8)  = [ 0.0320788303926323_Kr, 0.0320788303926323_Kr, 0.9037635088221031_Kr ]
         Xi(9)  = [ 0.0320788303926323_Kr, 0.9037635088221031_Kr, 0.0320788303926323_Kr ]
         Xi(10) = [ 0.4502229043567190_Kr, 0.0497770956432810_Kr, 0.0497770956432810_Kr ]
         Xi(11) = [ 0.0497770956432810_Kr, 0.4502229043567190_Kr, 0.0497770956432810_Kr ]
         Xi(12) = [ 0.0497770956432810_Kr, 0.0497770956432810_Kr, 0.4502229043567190_Kr ]
         Xi(13) = [ 0.0497770956432810_Kr, 0.4502229043567190_Kr, 0.4502229043567190_Kr ]
         Xi(14) = [ 0.4502229043567190_Kr, 0.0497770956432810_Kr, 0.4502229043567190_Kr ]
         Xi(15) = [ 0.4502229043567190_Kr, 0.4502229043567190_Kr, 0.0497770956432810_Kr ]
         Xi(16) = [ 0.3162695526014501_Kr, 0.1837304473985499_Kr, 0.1837304473985499_Kr ]
         Xi(17) = [ 0.1837304473985499_Kr, 0.3162695526014501_Kr, 0.1837304473985499_Kr ]
         Xi(18) = [ 0.1837304473985499_Kr, 0.1837304473985499_Kr, 0.3162695526014501_Kr ]
         Xi(19) = [ 0.1837304473985499_Kr, 0.3162695526014501_Kr, 0.3162695526014501_Kr ]
         Xi(20) = [ 0.3162695526014501_Kr, 0.1837304473985499_Kr, 0.3162695526014501_Kr ]
         Xi(21) = [ 0.3162695526014501_Kr, 0.3162695526014501_Kr, 0.1837304473985499_Kr ]
         Xi(22) = [ 0.0229177878448171_Kr, 0.2319010893971509_Kr, 0.2319010893971509_Kr ]
         Xi(23) = [ 0.2319010893971509_Kr, 0.0229177878448171_Kr, 0.2319010893971509_Kr ]
         Xi(24) = [ 0.2319010893971509_Kr, 0.2319010893971509_Kr, 0.0229177878448171_Kr ]
         Xi(25) = [ 0.5132800333608811_Kr, 0.2319010893971509_Kr, 0.2319010893971509_Kr ]
         Xi(26) = [ 0.2319010893971509_Kr, 0.5132800333608811_Kr, 0.2319010893971509_Kr ]
         Xi(27) = [ 0.2319010893971509_Kr, 0.2319010893971509_Kr, 0.5132800333608811_Kr ]
         Xi(28) = [ 0.2319010893971509_Kr, 0.0229177878448171_Kr, 0.5132800333608811_Kr ]
         Xi(29) = [ 0.0229177878448171_Kr, 0.5132800333608811_Kr, 0.2319010893971509_Kr ]
         Xi(30) = [ 0.5132800333608811_Kr, 0.2319010893971509_Kr, 0.0229177878448171_Kr ]
         Xi(31) = [ 0.2319010893971509_Kr, 0.5132800333608811_Kr, 0.0229177878448171_Kr ]
         Xi(32) = [ 0.0229177878448171_Kr, 0.2319010893971509_Kr, 0.5132800333608811_Kr ]
         Xi(33) = [ 0.5132800333608811_Kr, 0.0229177878448171_Kr, 0.2319010893971509_Kr ]
         Xi(34) = [ 0.7303134278075384_Kr, 0.0379700484718286_Kr, 0.0379700484718286_Kr ]
         Xi(35) = [ 0.0379700484718286_Kr, 0.7303134278075384_Kr, 0.0379700484718286_Kr ]
         Xi(36) = [ 0.0379700484718286_Kr, 0.0379700484718286_Kr, 0.7303134278075384_Kr ]
         Xi(37) = [ 0.1937464752488044_Kr, 0.0379700484718286_Kr, 0.0379700484718286_Kr ]
         Xi(38) = [ 0.0379700484718286_Kr, 0.1937464752488044_Kr, 0.0379700484718286_Kr ]
         Xi(39) = [ 0.0379700484718286_Kr, 0.0379700484718286_Kr, 0.1937464752488044_Kr ]
         Xi(40) = [ 0.0379700484718286_Kr, 0.7303134278075384_Kr, 0.1937464752488044_Kr ]
         Xi(41) = [ 0.7303134278075384_Kr, 0.1937464752488044_Kr, 0.0379700484718286_Kr ]
         Xi(42) = [ 0.1937464752488044_Kr, 0.0379700484718286_Kr, 0.7303134278075384_Kr ]
         Xi(43) = [ 0.0379700484718286_Kr, 0.1937464752488044_Kr, 0.7303134278075384_Kr ]
         Xi(44) = [ 0.7303134278075384_Kr, 0.0379700484718286_Kr, 0.1937464752488044_Kr ]
         Xi(45) = [ 0.1937464752488044_Kr, 0.7303134278075384_Kr, 0.0379700484718286_Kr ]
         dElem%Gauss_C(1)     = -0.2359620398477557_Kr / 6.0_Kr * detBinv
         dElem%Gauss_C(2:5)   =  0.0244878963560562_Kr / 6.0_Kr * detBinv
         dElem%Gauss_C(6:9)   =  0.0039485206398261_Kr / 6.0_Kr * detBinv
         dElem%Gauss_C(10:15) =  0.0263055529507371_Kr / 6.0_Kr * detBinv
         dElem%Gauss_C(16:21) =  0.0829803830550589_Kr / 6.0_Kr * detBinv
         dElem%Gauss_C(22:33) =  0.0254426245481023_Kr / 6.0_Kr * detBinv
         dElem%Gauss_C(34:45) =  0.0134324384376852_Kr / 6.0_Kr * detBinv
      Case Default
         Write(*,*) __FUNCT__,': Unimplemented quadrature order',dQuadratureOrder
         STOP
      End Select
      
      Select Case (dPolynomialOrder)
      Case(1)
         Num_DoF = 4
         Allocate(PhiHat(Num_DoF,Nb_Gauss),stat=ierr)
         Allocate(GradPhiHat(Num_DoF,Nb_Gauss),stat=ierr)
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
         Allocate(PhiHat(Num_DoF,Nb_Gauss),stat=ierr)
         Allocate(GradPhiHat(Num_DoF,Nb_Gauss),stat=ierr)
         PhiHat(1,:)  = (1.0_Kr - Xi%X - Xi%Y - Xi%Z) * (1.0_Kr - 2.0_Kr*Xi%X - 2.0_Kr*Xi%Y - 2.0_Kr*Xi%Z)
         PhiHat(2,:)  = Xi%X * (2.0_Kr * Xi%X - 1.0_Kr)
         PhiHat(3,:)  = Xi%Y * (2.0_Kr * Xi%Y - 1.0_Kr)
         PhiHat(4,:)  = Xi%Z * (2.0_Kr * Xi%Z - 1.0_Kr)
         PhiHat(5,:)  = 4.0_Kr * (1.0_Kr - Xi%X - Xi%Y - Xi%Z) * Xi%X
         PhiHat(6,:)  = 4.0_Kr *  Xi%X * Xi%Y
         PhiHat(7,:)  = 4.0_Kr * (1.0_Kr - Xi%X - Xi%Y - Xi%Z) * Xi%Y
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
         Num_DoF = 0
         Write(*,*) __FUNCT__,': Unimplemented PolynomialOrder',dPolynomialOrder
         STOP
      End Select
      
      Allocate (dElem%BF(Num_DoF,Nb_Gauss),stat=ierr) 
      Allocate (dElem%Grad_BF(Num_DoF,Nb_Gauss),stat=ierr)
      dElem%BF = PhiHat
      Do iDoF = 1,Num_DoF
         Do iG = 1,Nb_Gauss
            dElem%Grad_BF(iDoF,iG) = Bt * GradPhiHat(iDoF,iG) 
         End Do
      End Do
     
      DeAllocate(Xi,stat=ierr)
      DeAllocate(PhiHat,stat=ierr)
      DeAllocate(GradPhiHat,stat=ierr)
   End Subroutine Element_P_Lagrange_3D_Scal_Init

#undef __FUNCT__
#define __FUNCT__ "Element_P_Lagrange_3DBoundary_Scal_Init"
   Subroutine Element_P_Lagrange_3DBoundary_Scal_Init(dElem,area,outerNormal,dPolynomialOrder,dQuadratureOrder,ierr)
      ! Compute the quadrature weights and the value of the basis functions and their gradient 
      ! at the quadrature points.
      Type(MEF90Element3D_Scal)              :: dElem
      PetscReal,Intent(IN)                   :: area
      Type(Vect3D),Intent(IN)                :: outerNormal
      PetscInt,Intent(IN)                    :: dPolynomialOrder,dQuadratureOrder
      PetscErrorCode,Intent(OUT)             :: ierr
      
      PetscInt                               :: Nb_Gauss
      PetscInt                               :: Num_Dof
      PetscInt                               :: iDoF,iG

      Type(Vect2D),Dimension(:),Pointer      :: Xi ! The quadrature points coordinates in the reference element
      
      num_Dof = 0
      nb_Gauss = 0
      
      dElem%outerNormal = outerNormal

      Select Case (dQuadratureOrder)
      Case(0,1)
         Nb_Gauss = 1
         Allocate(Xi(Nb_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
         Xi(1)%X = 1.0_Kr / 3.0_Kr
         Xi(1)%Y = 1.0_Kr / 3.0_Kr
         dElem%Gauss_C = area

      Case(2)
         Nb_Gauss = 3
         Allocate(Xi(Nb_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
         Xi(1) = [ 1.0_Kr / 6.0_Kr,1.0_Kr / 6.0_Kr ]
         Xi(2) = [ 2.0_Kr / 3.0_Kr,1.0_Kr / 6.0_Kr ]
         Xi(3) = [ 1.0_Kr / 6.0_Kr,2.0_Kr / 3.0_Kr ]
         dElem%Gauss_C = area / 3.0_Kr

      Case(3)
         Nb_Gauss = 4
         Allocate(Xi(Nb_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
         dElem%Gauss_C    =  area * 25.0_Kr / 48.0_Kr
         dElem%Gauss_C(1) = -area * 9.0_Kr / 16.0_Kr
         Xi(1) = [ 1.0_Kr / 3.0_Kr,1.0_Kr / 3.0_Kr ]
         Xi(2) = [ 3.0_Kr / 5.0_Kr,1.0_Kr / 5.0_Kr ]
         Xi(3) = [ 1.0_Kr / 5.0_Kr,3.0_Kr / 5.0_Kr ]
         Xi(4) = [ 1.0_Kr / 5.0_Kr,1.0_Kr / 5.0_Kr ]
      Case(4)
         Nb_Gauss = 7
         Allocate(Xi(Nb_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
         dElem%Gauss_C(1) = area / 20.0_Kr
         dElem%Gauss_C(2) = area * 2.0_Kr / 15.0_Kr
         dElem%Gauss_C(3) = area / 20.0_Kr
         dElem%Gauss_C(4) = area * 2.0_Kr / 15.0_Kr
         dElem%Gauss_C(5) = area / 20.0_Kr
         dElem%Gauss_C(6) = area * 2.0_Kr / 15.0_Kr
         dElem%Gauss_C(7) = area * 9.0_Kr / 20.0_Kr
         Xi(1) = [ 0.0_Kr,0.0_Kr ]
         Xi(2) = [ 0.5_Kr,0.0_Kr ]
         Xi(3) = [ 1.0_Kr,0.0_Kr ]
         Xi(4) = [ 0.5_Kr,0.5_Kr ]
         Xi(5) = [ 0.0_Kr,1.0_Kr ]
         Xi(6) = [ 0.0_Kr,0.5_Kr ]
         Xi(7) = [ 1.0_Kr / 3.0_Kr,1.0_Kr / 3.0_Kr ]
      Case(5,6)
         Nb_Gauss = 9
         Allocate(Xi(Nb_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
         Xi(1) = [ 0.124949503233232_Kr, 0.437525248383384_Kr ]
         Xi(2) = [ 0.437525248383384_Kr, 0.124949503233232_Kr ]
         Xi(3) = [ 0.437525248383384_Kr, 0.437525248383384_Kr ]
         Xi(4) = [ 0.797112651860071_Kr, 0.165409927389841_Kr ]
         Xi(5) = [ 0.797112651860071_Kr, 0.037477420750088_Kr ]
         Xi(6) = [ 0.165409927389841_Kr, 0.797112651860071_Kr ]
         Xi(7) = [ 0.165409927389841_Kr, 0.037477420750088_Kr ]
         Xi(8) = [ 0.037477420750088_Kr, 0.797112651860071_Kr ]
         Xi(9) = [ 0.037477420750088_Kr, 0.165409927389841_Kr ]
         dElem%Gauss_C(1:3) = 0.205950504760887_Kr * area
         dElem%Gauss_C(4:9) = 0.063691414286223_Kr * area
      Case(7)
         Nb_Gauss = 13
         Allocate(Xi(Nb_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
         Xi(1)  = [ 0.333333333333333_Kr, 0.333333333333333_Kr ]
         Xi(2)  = [ 0.479308067841923_Kr, 0.260345966079038_Kr ]
         Xi(3)  = [ 0.260345966079038_Kr, 0.479308067841923_Kr ]
         Xi(4)  = [ 0.260345966079038_Kr, 0.260345966079038_Kr ]
         Xi(5)  = [ 0.869739794195568_Kr, 0.065130102902216_Kr ]
         Xi(6)  = [ 0.065130102902216_Kr, 0.869739794195568_Kr ]
         Xi(7)  = [ 0.065130102902216_Kr, 0.065130102902216_Kr ]
         Xi(8)  = [ 0.638444188569809_Kr, 0.312865496004875_Kr ]
         Xi(9)  = [ 0.638444188569809_Kr, 0.048690315425316_Kr ]
         Xi(10) = [ 0.312865496004875_Kr, 0.638444188569809_Kr ]
         Xi(11) = [ 0.312865496004875_Kr, 0.048690315425316_Kr ]
         Xi(12) = [ 0.048690315425316_Kr, 0.638444188569809_Kr ]
         Xi(13) = [ 0.048690315425316_Kr, 0.312865496004875_Kr ]
         dElem%Gauss_C(1)    = -0.149570044467670_Kr * area
         dElem%Gauss_C(2:4)  =  0.175615257433204_Kr * area
         dElem%Gauss_C(5:7)  =  0.053347235608839_Kr * area
         dElem%Gauss_C(8:13) =  0.077113760890257_Kr * area
      Case(8,9)
         Nb_Gauss = 19
         Allocate(Xi(Nb_Gauss),stat=ierr)
         Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
         Xi(1)  = [ 0.3333333333333333_Kr, 0.3333333333333333_Kr ]
         Xi(2)  = [ 0.0206349616025259_Kr, 0.4896825191987370_Kr ]
         Xi(3)  = [ 0.4896825191987370_Kr, 0.0206349616025259_Kr ]
         Xi(4)  = [ 0.4896825191987370_Kr, 0.4896825191987370_Kr ]
         Xi(5)  = [ 0.1258208170141290_Kr, 0.4370895914929355_Kr ]
         Xi(6)  = [ 0.4370895914929355_Kr, 0.1258208170141290_Kr ]
         Xi(7)  = [ 0.4370895914929355_Kr, 0.4370895914929355_Kr ]
         Xi(8)  = [ 0.6235929287619356_Kr, 0.1882035356190322_Kr ]
         Xi(9)  = [ 0.1882035356190322_Kr, 0.6235929287619356_Kr ]
         Xi(10) = [ 0.1882035356190322_Kr, 0.1882035356190322_Kr ]
         Xi(11) = [ 0.9105409732110941_Kr, 0.0447295133944530_Kr ]
         Xi(12) = [ 0.0447295133944530_Kr, 0.9105409732110941_Kr ]
         Xi(13) = [ 0.0447295133944530_Kr, 0.0447295133944530_Kr ]
         Xi(14) = [ 0.7411985987844980_Kr, 0.0368384120547363_Kr ]
         Xi(15) = [ 0.7411985987844980_Kr, 0.2219628891607657_Kr ]
         Xi(16) = [ 0.0368384120547363_Kr, 0.7411985987844980_Kr ]
         Xi(17) = [ 0.0368384120547363_Kr, 0.2219628891607657_Kr ]
         Xi(18) = [ 0.2219628891607657_Kr, 0.7411985987844980_Kr ]
         Xi(19) = [ 0.2219628891607657_Kr, 0.0368384120547363_Kr ]
         dElem%Gauss_C(1)     = 0.09713579628279610_Kr * area
         dElem%Gauss_C(2:4)   = 0.03133470022713983_Kr * area
         dElem%Gauss_C(5:7)   = 0.07782754100477543_Kr * area
         dElem%Gauss_C(8:10)  = 0.07964773892720910_Kr * area
         dElem%Gauss_C(11:13) = 0.02557767565869810_Kr * area
         dElem%Gauss_C(14:19) = 0.04328353937728940_Kr * area
      Case Default
         Write(*,*) __FUNCT__,': Unimplemented quadrature order',dQuadratureOrder
         ierr = PETSC_ERR_SUP
         STOP
      End Select
      
      Select Case (dPolynomialOrder)
      Case(1)
         Num_DoF = 3
         Allocate(dElem%BF(Num_DoF,Nb_Gauss),stat=ierr) 
         dElem%BF(1,:) = 1.0_Kr - Xi%X - Xi%Y
         dElem%BF(2,:) = Xi(:)%X
         dElem%BF(3,:) = Xi(:)%Y
      Case(2)
         Num_DoF = 6
         Allocate(dElem%BF(Num_DoF,Nb_Gauss),stat=ierr) 
         dElem%BF(1,:) = (1.0_Kr - Xi%X - Xi%Y) * (1.0_Kr - 2.0_Kr * Xi%X - 2.0_Kr * Xi%Y)      
         dElem%BF(2,:) = Xi%X * (2.0_Kr * Xi%X - 1.0_Kr)
         dElem%BF(3,:) = Xi%Y * (2.0_Kr * Xi%Y - 1.0_Kr)
         dElem%BF(4,:) = 4.0_Kr * Xi%X * (1.0_Kr - Xi%X - Xi%Y)
         dElem%BF(5,:) = 4.0_Kr * Xi%X * Xi%Y
         dElem%BF(6,:) = 4.0_Kr * Xi%Y * (1.0_Kr - Xi%X - Xi%Y)
      Case Default
         Write(*,*) __FUNCT__,': Unimplemented PolynomialOrder',dPolynomialOrder
         ierr = PETSC_ERR_SUP
      End Select
      
      Allocate(delem%Grad_BF(Num_DoF,Nb_Gauss),stat=ierr)
      Do iDof = 1, num_Dof
         Do iG = 1, Nb_Gauss
            delem%Grad_BF(iDof,iG) = 0.0_Kr
         End Do
      End Do
     
      DeAllocate(Xi,stat=ierr)
   End Subroutine Element_P_Lagrange_3DBoundary_Scal_Init
   
#undef __FUNCT__
#define __FUNCT__ "Element_P_Lagrange_3D_Vect_Init"
   Subroutine Element_P_Lagrange_3D_Vect_Init(dElem,Bt,detBinv,dPolynomialOrder,dQuadratureOrder,ierr)
      Type(MEF90Element3D_Vect)              :: dElem
      Type(Mat3D),Intent(IN)                 :: Bt
      PetscReal,Intent(IN)                   :: detBinv
      PetscInt,Intent(IN)                    :: dPolynomialOrder,dQuadratureOrder
      PetscErrorCode,Intent(OUT)             :: ierr
   
      Type(MEF90Element3D_Scal)              :: Elem_Scal
      PetscInt                               :: dim = 3 
      PetscInt                               :: Num_DoF,Nb_Gauss,i,iDof,iG
      
      
      Call Element_P_Lagrange_3D_Scal_Init(Elem_Scal,Bt,detBinv,dPolynomialOrder,dQuadratureOrder,ierr)
      Num_DoF   = Size(Elem_Scal%BF,1) 
      Nb_Gauss = Size(Elem_Scal%BF,2)
      Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
      Allocate(dElem%BF(Num_DoF * dim,Nb_Gauss),stat=ierr)
      Allocate(dElem%Grad_BF(Num_DoF * dim,Nb_Gauss),stat=ierr)
      Allocate(dElem%GradS_BF(Num_DoF * dim,Nb_Gauss),stat=ierr)
         
      dElem%Gauss_C = Elem_Scal%Gauss_C
      Do iDof = 1, num_dof * dim
         Do iG = 1, Nb_Gauss
            dElem%BF(iDof,iG)       = 0.0_Kr
            delem%Grad_BF(iDof,iG)  = 0.0_Kr
            delem%GradS_BF(iDof,iG) = 0.0_Kr
         End Do
      End Do
      
      Do i = 0,Num_DoF-1
         dElem%BF(i*dim+1,:)%X = Elem_Scal%BF(i+1,:)
         dElem%BF(i*dim+2,:)%Y = Elem_Scal%BF(i+1,:)
         dElem%BF(i*dim+3,:)%Z = Elem_Scal%BF(i+1,:)
         dElem%Grad_BF(i*dim+1,:)%XX = Elem_Scal%Grad_BF(i+1,:)%X
         dElem%Grad_BF(i*dim+1,:)%XY = Elem_Scal%Grad_BF(i+1,:)%Y
         dElem%Grad_BF(i*dim+1,:)%XZ = Elem_Scal%Grad_BF(i+1,:)%Z
         dElem%Grad_BF(i*dim+2,:)%YX = Elem_Scal%Grad_BF(i+1,:)%X
         dElem%Grad_BF(i*dim+2,:)%YY = Elem_Scal%Grad_BF(i+1,:)%Y
         dElem%Grad_BF(i*dim+2,:)%YZ = Elem_Scal%Grad_BF(i+1,:)%Z
         dElem%Grad_BF(i*dim+3,:)%ZX = Elem_Scal%Grad_BF(i+1,:)%X
         dElem%Grad_BF(i*dim+3,:)%ZY = Elem_Scal%Grad_BF(i+1,:)%Y
         dElem%Grad_BF(i*dim+3,:)%ZZ = Elem_Scal%Grad_BF(i+1,:)%Z

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
      Call MEF90Element_Destroy(Elem_Scal,ierr)
   End Subroutine Element_P_Lagrange_3D_Vect_Init

#undef __FUNCT__
#define __FUNCT__ "Element_P_Lagrange_3DBoundary_Vect_Init"
   Subroutine Element_P_Lagrange_3DBoundary_Vect_Init(dElem,area,outerNormal,dPolynomialOrder,dQuadratureOrder,ierr)
      Type(MEF90Element3D_Vect)              :: dElem
      PetscReal,Intent(IN)                   :: area
      Type(Vect3D),Intent(IN)                :: outerNormal
      PetscInt                               :: dPolynomialOrder,dQuadratureOrder
      PetscErrorCode,Intent(OUT)             :: ierr
   
      Type(MEF90Element3D_Scal)              :: Elem_Scal
      PetscInt                               :: dim = 3 
      PetscInt                               :: Num_DoF,Nb_Gauss,iDof,iG
      
      
      Call Element_P_Lagrange_3DBoundary_Scal_Init(Elem_Scal,area,outerNormal,dPolynomialOrder,dQuadratureOrder,ierr)
      Num_DoF   = Size(Elem_Scal%BF,1) 
      Nb_Gauss = Size(Elem_Scal%BF,2)
      Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
      Allocate(dElem%BF(Num_DoF * dim,Nb_Gauss),stat=ierr)
      Allocate(dElem%Grad_BF(Num_DoF * dim,Nb_Gauss),stat=ierr)
      Allocate(dElem%GradS_BF(Num_DoF * dim,Nb_Gauss),stat=ierr)
         
      dElem%Gauss_C = Elem_Scal%Gauss_C
      Do iDof = 1, num_dof * dim
         Do iG = 1, Nb_Gauss
            dElem%BF(iDof,iG)       = 0.0_Kr
            delem%Grad_BF(iDof,iG)  = 0.0_Kr
            delem%GradS_BF(iDof,iG) = 0.0_Kr
         End Do
      End Do
      Do iDof = 0,Num_DoF-1
         dElem%BF(iDof*dim+1,:)%X = Elem_Scal%BF(iDof+1,:)
         dElem%BF(iDof*dim+2,:)%Y = Elem_Scal%BF(iDof+1,:)
         dElem%BF(iDof*dim+3,:)%Z = Elem_Scal%BF(iDof+1,:)
      End Do

      dElem%outerNormal = Elem_Scal%outerNormal
      Call MEF90Element_Destroy(Elem_Scal,ierr)
   End Subroutine Element_P_Lagrange_3DBoundary_Vect_Init

#undef __FUNCT__
#define __FUNCT__ "Element2D_Scal_Destroy"
   Subroutine Element2D_Scal_Destroy(dElem,ierr)
      Type(MEF90Element2D_Scal)              :: dElem
      PetscErrorCode,Intent(OUT)             :: ierr
      
      If (Associated(dElem%BF)) Then
         DeAllocate(dElem%BF,stat=ierr)
      End If
      If (Associated(dElem%Grad_BF)) Then
         DeAllocate(dElem%Grad_BF,stat=ierr)
      End If
      If (Associated(dElem%Gauss_C)) Then
         DeAllocate(dElem%Gauss_C,stat=ierr)
      End If
   End Subroutine Element2D_Scal_Destroy

#undef __FUNCT__
#define __FUNCT__ "Element2D_Vect_Destroy"
   Subroutine Element2D_Vect_Destroy(dElem,ierr)
      Type(MEF90Element2D_Vect)              :: dElem
      PetscErrorCode,Intent(OUT)             :: ierr
      
      If (Associated(dElem%BF)) Then
         DeAllocate(dElem%BF,stat=ierr)
      End If
      If (Associated(dElem%Grad_BF)) Then
         DeAllocate(dElem%Grad_BF,stat=ierr)
      End If
      If (Associated(dElem%GradS_BF)) Then
         DeAllocate(dElem%GradS_BF,stat=ierr)
      End If
      If (Associated(dElem%Gauss_C)) Then
         DeAllocate(dElem%Gauss_C,stat=ierr)
      End If
   End Subroutine Element2D_Vect_Destroy
      
#undef __FUNCT__
#define __FUNCT__ "Element3D_Scal_Destroy"
   Subroutine Element3D_Scal_Destroy(dElem,ierr)
      Type(MEF90Element3D_Scal)              :: dElem
      PetscErrorCode,Intent(OUT)             :: ierr
      
      If (Associated(dElem%BF)) Then
         DeAllocate(dElem%BF,stat=ierr)
      End If
      If (Associated(dElem%Grad_BF)) Then
         DeAllocate(dElem%Grad_BF,stat=ierr)
      End If
      If (Associated(dElem%Gauss_C)) Then
         DeAllocate(dElem%Gauss_C,stat=ierr)
      End If
   End Subroutine Element3D_Scal_Destroy

#undef __FUNCT__
#define __FUNCT__ "Element3D_Vect_Destroy"
   Subroutine Element3D_Vect_Destroy(dElem,ierr)
      Type(MEF90Element3D_Vect)              :: dElem
      PetscErrorCode,Intent(OUT)             :: ierr
      If (Associated(dElem%BF)) Then
         DeAllocate(dElem%BF,stat=ierr)
      End If
      If (Associated(dElem%Grad_BF)) Then
         DeAllocate(dElem%Grad_BF,stat=ierr)
      End If
      If (Associated(dElem%GradS_BF)) Then
         DeAllocate(dElem%GradS_BF,stat=ierr)
      End If
      If (Associated(dElem%Gauss_C)) Then
         DeAllocate(dElem%Gauss_C,stat=ierr)
      End If
   End Subroutine Element3D_Vect_Destroy
   
#undef __FUNCT__
#define __FUNCT__ "Element2D_Scal_DestroySet"
   Subroutine Element2D_Scal_DestroySet(dElem,ierr)
      Type(MEF90Element2D_Scal),dimension(:),Pointer    :: dElem
      PetscErrorCode,Intent(OUT)                        :: ierr
      
      PetscInt                                          :: cell
      
      Do cell = 1, size(dElem)
         Call MEF90Element_Destroy(dElem(cell),ierr)
      End Do
      DeAllocate(dElem,stat=ierr)
   End Subroutine Element2D_Scal_DestroySet

#undef __FUNCT__
#define __FUNCT__ "Element2D_Vect_DestroySet"
   Subroutine Element2D_Vect_DestroySet(dElem,ierr)
      Type(MEF90Element2D_Vect),dimension(:),Pointer    :: dElem
      PetscErrorCode,Intent(OUT)                        :: ierr
      
      PetscInt                                          :: cell
      
      Do cell = 1, size(dElem)
         Call MEF90Element_Destroy(dElem(cell),ierr)
      End Do
      DeAllocate(dElem,stat=ierr)
   End Subroutine Element2D_Vect_DestroySet

#undef __FUNCT__
#define __FUNCT__ "Element3D_Scal_DestroySet"
   Subroutine Element3D_Scal_DestroySet(dElem,ierr)
      Type(MEF90Element3D_Scal),dimension(:),Pointer    :: dElem
      PetscErrorCode,Intent(OUT)                        :: ierr
      
      PetscInt                                          :: cell
      
      Do cell = 1, size(dElem)
         Call MEF90Element_Destroy(dElem(cell),ierr)
      End Do
      DeAllocate(dElem,stat=ierr)
   End Subroutine Element3D_Scal_DestroySet

#undef __FUNCT__
#define __FUNCT__ "Element3D_Vect_DestroySet"
   Subroutine Element3D_Vect_DestroySet(dElem,ierr)
      Type(MEF90Element3D_Vect),dimension(:),Pointer    :: dElem
      PetscErrorCode,Intent(OUT)                        :: ierr
      
      PetscInt                                          :: cell
        
      Do cell = 1, size(dElem)
         Call MEF90Element_Destroy(dElem(cell),ierr)
      End Do
      DeAllocate(dElem,stat=ierr)
   End Subroutine Element3D_Vect_DestroySet
   
End Module m_MEF90_Elements
