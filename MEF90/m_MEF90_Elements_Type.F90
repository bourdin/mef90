Module m_MEF90_Elements_Type
#include "petsc/finclude/petsc.h"
    Use m_MEF90_LinAlg
    Use m_MEF90_Utils
    Use m_MEF90_Parameters
    Use petsc
    IMPLICIT NONE
   
    Public :: MEF90ElementType
    
    Type MEF90ElementType
        ! name is the element name in english language
        Character(len=MEF90MXSTRLEN)                 :: Name
        
        !!! Number of dof at each 
        !!!   vertex edge face cell in 3D,
        !!!   vertex edge cell in 2D (last value is unused, i.e. always use depth to locate DoF)
        !!! i.e. point depth -1
        PetscInt,Dimension(4)                        :: numDofs
        PetscInt                                     :: numDoF 
        !!! numDof = numVertexDof * numVertex + numEdgeDof * numEdge 
        !!!          + numFaceDof * numFace + numCellDof
        PetscInt                                     :: dim
        !!! The dimension of the cell associated with the element
        PetscInt                                     :: coDim
        !!! The co-dimension of the cell associated with the element
        PetscInt                                     :: order
        !!! The polynomial order
    End Type MEF90ElementType
    
    Enum,bind(c)
        enumerator :: &
        MEF90ElementFamilyLagrange = 0
    End Enum

    Character(kind=c_char,len=MEF90MXSTRLEN),dimension(4),Parameter,Public   :: MEF90ElementFamily = [ &
        "Lagrange           ",     &  ! 0
        "MEF90ElementFamily ",     &
        "prefix_            ",     &
        C_NULL_CHAR//"                  "]
        Character(len=MEF90MXSTRLEN),dimension(4),protected  :: MEF90ElementFamilyList


    Type(MEF90ElementType),Parameter,Public :: MEF90_NULL_ELEMENT = MEF90ElementType(   &
    "NULL",                                      &  ! name
    [0,0,0,0],0,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
    0,0,0                                        &  ! dim,codim,order
    )

    !!! 
    !!! P0 pseudo-elements. These are used for creating sections only
    !!!
    Type(MEF90ElementType),Parameter,Public :: MEF90P0Lagrange2D = MEF90ElementType(   &
    "MEF90P0Lagrange2D",                         &  ! name
    [0,0,1,0],1,                                 &  ! numVertexDof,numEdgeDof,numCellDof,unused,numDof
    2,0,0                                        &  ! dim,codim,order
    )
    Type(MEF90ElementType),Parameter,Public :: MEF90P0Lagrange3D = MEF90ElementType(   &
    "MEF90P0Lagrange3D",                         &  ! name
    [0,0,0,1],1,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
    3,0,0                                        &  ! dim,codim,order
    )
    Type(MEF90ElementType),Parameter,Public :: MEF90P0Lagrange2DBoundary = MEF90ElementType(   &
    "MEF90P0Lagrange2DBoundary",                 &  ! name
    [0,1,0,0],1,                                 &  ! numVertexDof,numEdgeDof,numCellDof,unused,numDof
    2,0,0                                        &  ! dim,codim,order
    )
    Type(MEF90ElementType),Parameter,Public :: MEF90P0Lagrange3DBoundary = MEF90ElementType(   &
    "MEF90P0Lagrange3DBoundary",                 &  ! name
    [0,0,1,0],1,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
    3,0,0                                        &  ! dim,codim,order
    )
    !!! 
    !!! P1 Simplicial Lagrange elements
    !!!
    Type(MEF90ElementType),Parameter,Public :: MEF90P1Lagrange2D = MEF90ElementType(   &
    "MEF90P1Lagrange2D",                         &  ! name
    [1,0,0,0],3,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
    2,0,1                                        &  ! dim,codim,order
    )
    Type(MEF90ElementType),Parameter,Public :: MEF90P1Lagrange3D = MEF90ElementType(   &
    "MEF90P1Lagrange3D",                         &  ! name
    [1,0,0,0],4,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
    3,0,1                                        &  ! dim,codim,order
    )
    Type(MEF90ElementType),Parameter,Public :: MEF90P1Lagrange2DBoundary = MEF90ElementType(   &
    "MEF90P1Lagrange2DBoundary",                 &  ! name
    [1,0,0,0],2,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
    2,1,1                                        &  ! dim,codim,order
    )
    Type(MEF90ElementType),Parameter,Public :: MEF90P1Lagrange3DBoundary = MEF90ElementType(   &
    "MEF90P1Lagrange3DBoundary",                 &  ! name
    [1,0,0,0],3,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
    3,1,1                                        &  ! dim,codim,order
    )
    !!!
    !!! P2 Simplicial Lagrange elements
    !!!
    Type(MEF90ElementType),Parameter,Public :: MEF90P2Lagrange2D = MEF90ElementType(   &
    "MEF90P2Lagrange2D",                         &  ! name
    [1,1,0,0],6,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
    2,0,2                                        &  ! dim,codim,order
    )
    Type(MEF90ElementType),Parameter,Public :: MEF90P2Lagrange3D = MEF90ElementType(   &
    "MEF90P2Lagrange3D",                         &  ! name
    [1,1,0,0],10,                                &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
    3,0,2                                        &  ! dim,codim,order
    )
    Type(MEF90ElementType),Parameter,Public :: MEF90P2Lagrange2DBoundary = MEF90ElementType(   &
    "MEF90P2Lagrange2DBoundary",                 &  ! name
    [1,1,0,0],3,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
    2,1,2                                        &  ! dim,codim,order
    )
    Type(MEF90ElementType),Parameter,Public :: MEF90P2Lagrange3DBoundary = MEF90ElementType(   &
    "MEF90P2Lagrange3DBoundary",                 &  ! name
    [1,1,0,0],6,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
    3,1,2                                        &  ! dim,codim,order
    )
    !!! 
    !!! Q1 tensor product Lagrange elements
    !!!
    Type(MEF90ElementType),Parameter,Public :: MEF90Q1Lagrange2D = MEF90ElementType(   &
    "MEF90Q1Lagrange2D",                         &  ! name
    [1,0,0,0],4,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
    2,0,1                                        &  ! dim,codim,order
    )
    Type(MEF90ElementType),Parameter,Public :: MEF90Q1Lagrange3D = MEF90ElementType(   &
    "MEF90Q1Lagrange3D",                         &  ! name
    [1,0,0,0],8,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
    3,0,1                                        &  ! dim,codim,order
    )
    Type(MEF90ElementType),Parameter,Public :: MEF90Q1Lagrange2DBoundary = MEF90ElementType(   &
    "MEF90Q1Lagrange2DBoundary",                 &  ! name
    [1,0,0,0],2,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
    2,1,1                                        &  ! dim,codim,order
    )
    Type(MEF90ElementType),Parameter,Public :: MEF90Q1Lagrange3DBoundary = MEF90ElementType(   &
    "MEF90Q1Lagrange3DBoundary",                 &  ! name
    [1,0,0,0],4,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
    3,1,1                                        &  ! dim,codim,order
    )
    !!! 
    !!! Q2 tensor product Lagrange elements
    !!!
    Type(MEF90ElementType),Parameter,Public :: MEF90Q2Lagrange2D = MEF90ElementType(   &
    "MEF90Q2Lagrange2D",                         &  ! name
    [1,1,0,1],9,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
    2,0,1                                        &  ! dim,codim,order
    )
    Type(MEF90ElementType),Parameter,Public :: MEF90Q2Lagrange3D = MEF90ElementType(   &
    "MEF90Q2Lagrange3D",                         &  ! name
    [1,1,1,1],27,                                &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
    3,0,1                                        &  ! dim,codim,order
    )
    Type(MEF90ElementType),Parameter,Public :: MEF90Q2Lagrange2DBoundary = MEF90ElementType(   &
    "MEF90Q2Lagrange2DBoundary",                 &  ! name
    [1,1,0,0],3,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
    2,1,1                                        &  ! dim,codim,order
    )
    Type(MEF90ElementType),Parameter,Public :: MEF90Q2Lagrange3DBoundary = MEF90ElementType(   &
    "MEF90Q2Lagrange3DBoundary",                 &  ! name
    [1,1,0,1],9,                                 &  ! numVertexDof,numEdgeDof,numFaceDof,numCellDof,numDof
    3,1,1                                        &  ! dim,codim,order
    )       
Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90ElementsInitialize_Private"
    !!!
    !!!  
    !!!  MEF90ElementsInitialize_Private:
    !!!  
    !!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
    !!!

    Subroutine MEF90ElementsInitialize_Private(ierr)
        PetscErrorCode,Intent(OUT)                   :: ierr

        MEF90ElementFamilyList(1) = 'Lagrange'
        MEF90ElementFamilyList(2) = 'MEF90ElementFamily'
        MEF90ElementFamilyList(3) = '_MEF90ElementFamily'
        MEF90ElementFamilyList(4) = ''
        ierr = 0
    End Subroutine MEF90ElementsInitialize_Private
End Module m_MEF90_Elements_Type
