Module m_MEF_Elements
#include "finclude/petscdef.h"
   Use m_MEF_LinAlg
   Use m_MEF_Types
   Use m_MEF_Utils
   Use m_MEF_Parameters
   Use petsc
   IMPLICIT NONE

   Private   
   Public :: ElementInit
   Public :: ElementDestroy
   Public :: ElementView
   Public :: cellSetElemTypeInit


   Interface ElementInit
      Module Procedure Element2D_Scal_Init,Element2D_Init,Element2D_Elast_Init,Element3D_Scal_Init,Element3D_Init, &
      Element3D_Elast_Init,Element2D_Scal_InitAll,Element2D_InitAll,Element2D_Elast_InitAll,Element3D_Scal_InitAll, &
      Element3D_InitAll,Element3D_Elast_InitAll
   End Interface ElementInit
   
   Interface ElementDestroy
      Module Procedure Element2D_Scal_Destroy,Element2D_Destroy,Element2D_Elast_Destroy,Element3D_Scal_Destroy, &
      Element3D_Destroy,Element3D_Elast_Destroy
   End Interface ElementDestroy
   
   Interface ElementView
      Module Procedure Element2D_Scal_View,Element2D_View,Element2D_Elast_View, &
      Element3D_Scal_View,Element3D_View,Element3D_Elast_View
   End Interface ElementView
   
 Contains
#undef __FUNCT__
#define __FUNCT__ "cellSetElemTypeInit"
   Subroutine cellSetElemTypeInit(set,dDim)
      Type(CellSet_Type)                     :: set
      PetscInt                               :: dDim
      
      Select Case (dDim)
      Case (2)
         Select Case (set%ElemType)
         Case (MEF90_P1_Lagrange)         
            set%dofLocation = (/ 0,0,0,3 /)
            set%Codimension = 0
         Case (MEF90_P1_Lagrange_Boundary)         
            set%dofLocation = (/ 0,0,0,2 /)
            set%Codimension = 1
         Case (MEF90_P2_Lagrange)
            set%dofLocation = (/ 0,0,3,3 /)
            set%Codimension = 0
         Case (MEF90_P2_Lagrange_Boundary)         
            set%dofLocation = (/ 0,0,1,2 /)
            set%Codimension = 1
         Case Default
            Print*,__FUNCT__,': Unknown element type',set%ElemType
            STOP
         End Select
      Case (3)
         Select Case (set%ElemType)
         Case (MEF90_P1_Lagrange)         
            set%dofLocation = (/ 0,0,0,4 /)
            set%Codimension = 0
         Case (MEF90_P1_Lagrange_Boundary)         
            set%dofLocation = (/ 0,0,0,3 /)
            set%Codimension = 1
         Case (MEF90_P2_Lagrange)
            set%dofLocation = (/ 0,0,6,4 /)
            set%Codimension = 0
         Case (MEF90_P2_Lagrange_Boundary)
            set%dofLocation = (/ 0,0,3,3 /)
            set%Codimension = 1
         Case Default
            Print*,__FUNCT__,': Unknown element type',set%ElemType
            STOP
         End Select
      Case Default
         Print*,__FUNCT__,': Unknown dimension',dDim
         STOP
      End Select
      set%numDoF = sum(set%dofLocation)
   End Subroutine cellSetElemTypeInit     
   
#undef __FUNCT__
#define __FUNCT__ "Element2D_Scal_InitAll"
!!!
!!! This routine does not use the values cached in dMeshtopology that are also directly accessible from the mesh
!!!
   Subroutine Element2D_Scal_InitAll(dMeshTopology,dElem,dQuadratureOrder)
      Type(MeshTopology_Type)                     :: dMeshTopology
      Type(Element2D_Scal),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                         :: dQuadratureOrder
      
      PetscInt                                    :: set,setID,iELoc,iE,ierr
      PetscInt                                    :: numDim,numCells,numVertices
      PetscInt                                    :: numVertexinCell,point
      PetscInt,Dimension(:),Pointer               :: cone
      Type(IS)                                    :: CellIS,CellSetIS
      PetscInt,Dimension(:),Pointer               :: CellID,CellSetID
      PetscReal,Dimension(:,:),Pointer            :: Coords
      PetscReal,Dimension(:),Pointer              :: TmpCoords
      Type(SectionReal)                           :: CoordSection
      
      Call DMMeshGetDimension(dMeshTopology%mesh,numDim,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumSize(dMeshTopology%mesh,"height",0,numCells,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumSize(dMeshTopology%mesh,"depth",0,numVertices,ierr);CHKERRQ(ierr)

      Allocate(dElem(numCells))
      Call DMMeshGetLabelIdIS(dMeshTopology%mesh,'Cell Sets',CellSetIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetIS,CellSetID,ierr);CHKERRQ(ierr)

      !!! Initialize the Basis Functions in each element
      Call DMMeshGetSectionReal(dMeshTopology%mesh,'coordinates',CoordSection,ierr);CHKERRQ(ierr)
      Do_CellSet: Do set = 1,size(CellSetID)
         setID = CellSetID(set)
         Call DMMeshGetStratumIS(dMeshTopology%mesh,'Cell Sets',setId,CellIS,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(CellIS,CellID,ierr);CHKERRQ(ierr)

         !!! Get the number of vertices per cell for the cells in the set
         !!! We always assume that all cells in a set are of the same type so
         !!! we only need to do this once.
         !!! This is a bit crazy but will get much better when we switch to DMComplex
         Call DMMeshGetConeSize(dMeshTopology%mesh,CellID(1),numVertexinCell,ierr);CHKERRQ(ierr)    
         !!! The cone may contain faces and edges in the case of an interpolated mesh
         !!! so we need to filter out based on the vertex point number:
         Allocate(cone(numVertexinCell))
         numVertexinCell = 0
         Do point = 1, size(cone)
            If ( (cone(point) >= numCells) .AND. (cone(point) < numCells + numVertices) ) Then
               numVertexinCell = numVertexinCell + 1
            End If
         End Do
         DeAllocate(cone)
         
         Allocate(TmpCoords(numDim * numVertexinCell))
         Allocate(Coords   (numDim,  numVertexinCell))

         Do_Elem_iE: Do iELoc = 1,size(CellID)
            iE = CellID(iELoc)
            Call SectionRealRestrictClosure(CoordSection,dMeshTopology%mesh,iE-1,Size(TmpCoords),TmpCoords,ierr);CHKERRQ(ierr)
            Coords = Reshape(TmpCoords,(/numDim,numVertexinCell /) )
            !!! WTF? why not reshaping the arguments in Init_Element? 
            Call ElementInit(dElem(iE),Coords,dQuadratureOrder,dMeshTopology%cellSet(set)%ElemType)
         End Do Do_Elem_iE
         Call ISRestoreIndicesF90(CellIS,CellID,ierr);CHKERRQ(ierr)
         DeAllocate(TmpCoords)
         DeAllocate(Coords)
      End Do Do_CellSet
      Call SectionRealDestroy(CoordSection,ierr);CHKERRQ(ierr)
      Call ISRestoreIndicesF90(CellSetIS,CellSetID,ierr);CHKERRQ(ierr)
   End Subroutine Element2D_Scal_InitAll

#undef __FUNCT__
#define __FUNCT__ "Element2D_Scal_Init"
   Subroutine Element2D_Scal_Init(dElem,dCoord,QuadratureOrder,Element_Type)
      Type(Element2D_Scal)                   :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord
      PetscInt,Intent(IN)                    :: QuadratureOrder
      PetscInt,Intent(IN)                    :: Element_Type
      
      Select Case (Element_Type)
         Case (MEF90_P1_Lagrange)
            Call Element_P_Lagrange_2D_Scal_Init(dElem,dCoord,1,QuadratureOrder)

         Case (MEF90_P2_Lagrange)
            Call Element_P_Lagrange_2D_Scal_Init(dElem,dCoord,2,QuadratureOrder)

         Case (MEF90_P1_Lagrange_Boundary)
            !Call Element_P_Lagrange_2D_Scal_Init(dElem,dCoord,1,QuadratureOrder)

         Case (MEF90_P2_Lagrange_Boundary)
            !Call Element_P_Lagrange_2D_Scal_Init(dElem,dCoord,2,QuadratureOrder)

!         Case (MEF90_Q1_Lagrange)
!            Call Element_Q_Lagrange_2D_Scal_Init(dElem,dCoord,1,QuadratureOrder)
!         Case (MEF90_Q2_Lagrange)
!            Call Element_Q_Lagrange_2D_Scal_Init(dElem,dCoord,2,QuadratureOrder)
         Case Default
            Print*,__FUNCT__,': Element type not implemented yet',Element_Type
      End Select
   End Subroutine Element2D_Scal_Init                                
   
#undef __FUNCT__
#define __FUNCT__ "Element2D_InitAll"
!!!
!!! This routine does not use the values cached in dMeshtopology that are also directly accessible from the mesh
!!!
   Subroutine Element2D_InitAll(dMeshTopology,dElem,dQuadratureOrder)
      Type(MeshTopology_Type)                     :: dMeshTopology
      Type(Element2D),Dimension(:),Pointer        :: dElem
      PetscInt,Intent(IN)                         :: dQuadratureOrder
      
      PetscInt                                    :: set,setID,iELoc,iE,ierr
      PetscInt                                    :: numDim,numCells,numVertices
      PetscInt                                    :: numVertexinCell,point
      PetscInt,Dimension(:),Pointer               :: cone
      Type(IS)                                    :: CellIS,CellSetIS
      PetscInt,Dimension(:),Pointer               :: CellID,CellSetID
      PetscReal,Dimension(:,:),Pointer            :: Coords
      PetscReal,Dimension(:),Pointer              :: TmpCoords
      Type(SectionReal)                           :: CoordSection
      
      Call DMMeshGetDimension(dMeshTopology%mesh,numDim,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumSize(dMeshTopology%mesh,"height",0,numCells,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumSize(dMeshTopology%mesh,"depth",0,numVertices,ierr);CHKERRQ(ierr)

      Allocate(dElem(numCells))
      Call DMMeshGetLabelIdIS(dMeshTopology%mesh,'Cell Sets',CellSetIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetIS,CellSetID,ierr);CHKERRQ(ierr)

      !!! Initialize the Basis Functions in each element
      Call DMMeshGetSectionReal(dMeshTopology%mesh,'coordinates',CoordSection,ierr);CHKERRQ(ierr)
      Do_CellSet: Do set = 1,size(CellSetID)
         setID = CellSetID(set)
         Call DMMeshGetStratumIS(dMeshTopology%mesh,'Cell Sets',setId,CellIS,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(CellIS,CellID,ierr);CHKERRQ(ierr)

         !!! Get the number of vertices per cell for the cells in the set
         !!! We always assume that all cells in a set are of the same type so
         !!! we only need to do this once.
         !!! This is a bit crazy but will get much better when we switch to DMComplex
         Call DMMeshGetConeSize(dMeshTopology%mesh,CellID(1),numVertexinCell,ierr);CHKERRQ(ierr)    
         !!! The cone may contain faces and edges in the case of an interpolated mesh
         !!! so we need to filter out based on the vertex point number:
         Allocate(cone(numVertexinCell))
         numVertexinCell = 0
         Do point = 1, size(cone)
            If ( (cone(point) >= numCells) .AND. (cone(point) < numCells + numVertices) ) Then
               numVertexinCell = numVertexinCell + 1
            End If
         End Do
         DeAllocate(cone)
         
         Allocate(TmpCoords(numDim * numVertexinCell))
         Allocate(Coords   (numDim,  numVertexinCell))

         Do_Elem_iE: Do iELoc = 1,size(CellID)
            iE = CellID(iELoc)
            Call SectionRealRestrictClosure(CoordSection,dMeshTopology%mesh,iE-1,Size(TmpCoords),TmpCoords,ierr);CHKERRQ(ierr)
            Coords = Reshape(TmpCoords,(/numDim,numVertexinCell /) )
            !!! WTF? why not reshaping the arguments in Init_Element? 
            Call ElementInit(dElem(iE),Coords,dQuadratureOrder,dMeshTopology%cellSet(set)%ElemType)
         End Do Do_Elem_iE
         Call ISRestoreIndicesF90(CellIS,CellID,ierr);CHKERRQ(ierr)
         DeAllocate(TmpCoords)
         DeAllocate(Coords)
      End Do Do_CellSet
      Call SectionRealDestroy(CoordSection,ierr);CHKERRQ(ierr)
      Call ISRestoreIndicesF90(CellSetIS,CellSetID,ierr);CHKERRQ(ierr)
   End Subroutine Element2D_InitAll

#undef __FUNCT__
#define __FUNCT__ "Element2D_Init"
   Subroutine Element2D_Init(dElem,dCoord,QuadratureOrder,Element_Type)
      Type(Element2D)                        :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord
      PetscInt,Intent(IN)                    :: QuadratureOrder
      PetscInt,Intent(IN)                    :: Element_Type
      
      Select Case (Element_Type)
         Case (MEF90_P1_Lagrange)
            Call Element_P_Lagrange_2D_Init(dElem,dCoord,1,QuadratureOrder)

         Case (MEF90_P2_Lagrange)
            Call Element_P_Lagrange_2D_Init(dElem,dCoord,2,QuadratureOrder)

!         Case (MEF90_Q1_Lagrange)
!            Call Element_Q_Lagrange_2D_Scal_Init(dElem,dCoord,1,QuadratureOrder)
!         Case (MEF90_Q2_Lagrange)
!            Call Element_Q_Lagrange_2D_Scal_Init(dElem,dCoord,2,QuadratureOrder)
         Case Default
            Print*,__FUNCT__,': Element type not implemented yet',Element_Type
      End Select
   End Subroutine Element2D_Init                                

#undef __FUNCT__
#define __FUNCT__ "Element2D_Elast_InitAll"
!!!
!!! This routine does not use the values cached in dMeshtopology that are also directly accessible from the mesh
!!!
   Subroutine Element2D_Elast_InitAll(dMeshTopology,dElem,dQuadratureOrder)
      Type(MeshTopology_Type)                     :: dMeshTopology
      Type(Element2D_Elast),Dimension(:),Pointer  :: dElem
      PetscInt,Intent(IN)                         :: dQuadratureOrder
      
      PetscInt                                    :: set,setID,iELoc,iE,ierr
      PetscInt                                    :: numDim,numCells,numVertices
      PetscInt                                    :: numVertexinCell,point
      PetscInt,Dimension(:),Pointer               :: cone
      Type(IS)                                    :: CellIS,CellSetIS
      PetscInt,Dimension(:),Pointer               :: CellID,CellSetID
      PetscReal,Dimension(:,:),Pointer            :: Coords
      PetscReal,Dimension(:),Pointer              :: TmpCoords
      Type(SectionReal)                           :: CoordSection
      
      Call DMMeshGetDimension(dMeshTopology%mesh,numDim,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumSize(dMeshTopology%mesh,"height",0,numCells,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumSize(dMeshTopology%mesh,"depth",0,numVertices,ierr);CHKERRQ(ierr)

      Allocate(dElem(numCells))
      Call DMMeshGetLabelIdIS(dMeshTopology%mesh,'Cell Sets',CellSetIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetIS,CellSetID,ierr);CHKERRQ(ierr)

      !!! Initialize the Basis Functions in each element
      Call DMMeshGetSectionReal(dMeshTopology%mesh,'coordinates',CoordSection,ierr);CHKERRQ(ierr)
      Do_CellSet: Do set = 1,size(CellSetID)
         setID = CellSetID(set)
         Call DMMeshGetStratumIS(dMeshTopology%mesh,'Cell Sets',setId,CellIS,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(CellIS,CellID,ierr);CHKERRQ(ierr)

         !!! Get the number of vertices per cell for the cells in the set
         !!! We always assume that all cells in a set are of the same type so
         !!! we only need to do this once.
         !!! This is a bit crazy but will get much better when we switch to DMComplex
         Call DMMeshGetConeSize(dMeshTopology%mesh,CellID(1),numVertexinCell,ierr);CHKERRQ(ierr)    
         !!! The cone may contain faces and edges in the case of an interpolated mesh
         !!! so we need to filter out based on the vertex point number:
         Allocate(cone(numVertexinCell))
         numVertexinCell = 0
         Do point = 1, size(cone)
            If ( (cone(point) >= numCells) .AND. (cone(point) < numCells + numVertices) ) Then
               numVertexinCell = numVertexinCell + 1
            End If
         End Do
         DeAllocate(cone)
         
         Allocate(TmpCoords(numDim * numVertexinCell))
         Allocate(Coords   (numDim,  numVertexinCell))

         Do_Elem_iE: Do iELoc = 1,size(CellID)
            iE = CellID(iELoc)
            Call SectionRealRestrictClosure(CoordSection,dMeshTopology%mesh,iE-1,Size(TmpCoords),TmpCoords,ierr);CHKERRQ(ierr)
            Coords = Reshape(TmpCoords,(/numDim,numVertexinCell /) )
            !!! WTF? why not reshaping the arguments in Init_Element? 
            Call ElementInit(dElem(iE),Coords,dQuadratureOrder,dMeshTopology%cellSet(set)%ElemType)
         End Do Do_Elem_iE
         Call ISRestoreIndicesF90(CellIS,CellID,ierr);CHKERRQ(ierr)
         DeAllocate(TmpCoords)
         DeAllocate(Coords)
      End Do Do_CellSet
      Call SectionRealDestroy(CoordSection,ierr);CHKERRQ(ierr)
      Call ISRestoreIndicesF90(CellSetIS,CellSetID,ierr);CHKERRQ(ierr)
   End Subroutine Element2D_Elast_InitAll

#undef __FUNCT__
#define __FUNCT__ "Element2D_Elast_Init"
   Subroutine Element2D_Elast_Init(dElem,dCoord,QuadratureOrder,Element_Type)
      Type(Element2D_Elast)                  :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord
      PetscInt,Intent(IN)                    :: QuadratureOrder
      PetscInt,Intent(IN)                    :: Element_Type
      
      Select Case (Element_Type)
         Case (MEF90_P1_Lagrange)
            Call Element_P_Lagrange_2D_Elast_Init(dElem,dCoord,1,QuadratureOrder)

         Case (MEF90_P2_Lagrange)
            Call Element_P_Lagrange_2D_Elast_Init(dElem,dCoord,2,QuadratureOrder)

!         Case (MEF90_Q1_Lagrange)
!            Call Element_Q_Lagrange_2D_Scal_Init(dElem,dCoord,1,QuadratureOrder)
!         Case (MEF90_Q2_Lagrange)
!            Call Element_Q_Lagrange_2D_Scal_Init(dElem,dCoord,2,QuadratureOrder)
         Case Default
            Print*,__FUNCT__,': Element type not implemented yet',Element_Type
      End Select
   End Subroutine Element2D_Elast_Init


#undef __FUNCT__
#define __FUNCT__ "Element3D_Scal_InitAll"
!!!
!!! This routine does not use the values cached in dMeshtopology that are also directly accessible from the mesh
!!!
   Subroutine Element3D_Scal_InitAll(dMeshTopology,dElem,dQuadratureOrder)
      Type(MeshTopology_Type)                     :: dMeshTopology
      Type(Element3D_Scal),Dimension(:),Pointer   :: dElem
      PetscInt,Intent(IN)                         :: dQuadratureOrder
      
      PetscInt                                    :: set,setID,iELoc,iE,ierr
      PetscInt                                    :: numDim,numCells,numVertices
      PetscInt                                    :: numVertexinCell,point
      PetscInt,Dimension(:),Pointer               :: cone
      Type(IS)                                    :: CellIS,CellSetIS
      PetscInt,Dimension(:),Pointer               :: CellID,CellSetID
      PetscReal,Dimension(:,:),Pointer            :: Coords
      PetscReal,Dimension(:),Pointer              :: TmpCoords
      Type(SectionReal)                           :: CoordSection
      
      Call DMMeshGetDimension(dMeshTopology%mesh,numDim,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumSize(dMeshTopology%mesh,"height",0,numCells,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumSize(dMeshTopology%mesh,"depth",0,numVertices,ierr);CHKERRQ(ierr)

      Allocate(dElem(numCells))
      Call DMMeshGetLabelIdIS(dMeshTopology%mesh,'Cell Sets',CellSetIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetIS,CellSetID,ierr);CHKERRQ(ierr)

      !!! Initialize the Basis Functions in each element
      Call DMMeshGetSectionReal(dMeshTopology%mesh,'coordinates',CoordSection,ierr);CHKERRQ(ierr)
      Do_CellSet: Do set = 1,size(CellSetID)
         setID = CellSetID(set)
         Call DMMeshGetStratumIS(dMeshTopology%mesh,'Cell Sets',setId,CellIS,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(CellIS,CellID,ierr);CHKERRQ(ierr)

         !!! Get the number of vertices per cell for the cells in the set
         !!! We always assume that all cells in a set are of the same type so
         !!! we only need to do this once.
         !!! This is a bit crazy but will get much better when we switch to DMComplex
         Call DMMeshGetConeSize(dMeshTopology%mesh,CellID(1),numVertexinCell,ierr);CHKERRQ(ierr)    
         !!! The cone may contain faces and edges in the case of an interpolated mesh
         !!! so we need to filter out based on the vertex point number:
         Allocate(cone(numVertexinCell))
         numVertexinCell = 0
         Do point = 1, size(cone)
            If ( (cone(point) >= numCells) .AND. (cone(point) < numCells + numVertices) ) Then
               numVertexinCell = numVertexinCell + 1
            End If
         End Do
         DeAllocate(cone)
         
         Allocate(TmpCoords(numDim * numVertexinCell))
         Allocate(Coords   (numDim,  numVertexinCell))

         Do_Elem_iE: Do iELoc = 1,size(CellID)
            iE = CellID(iELoc)
            Call SectionRealRestrictClosure(CoordSection,dMeshTopology%mesh,iE-1,Size(TmpCoords),TmpCoords,ierr);CHKERRQ(ierr)
            Coords = Reshape(TmpCoords,(/numDim,numVertexinCell /) )
            !!! WTF? why not reshaping the arguments in Init_Element? 
            Call ElementInit(dElem(iE),Coords,dQuadratureOrder,dMeshTopology%cellSet(set)%ElemType)
         End Do Do_Elem_iE
         Call ISRestoreIndicesF90(CellIS,CellID,ierr);CHKERRQ(ierr)
         DeAllocate(TmpCoords)
         DeAllocate(Coords)
      End Do Do_CellSet
      Call SectionRealDestroy(CoordSection,ierr);CHKERRQ(ierr)
      Call ISRestoreIndicesF90(CellSetIS,CellSetID,ierr);CHKERRQ(ierr)
   End Subroutine Element3D_Scal_InitAll

#undef __FUNCT__
#define __FUNCT__ "Element3D_Scal_Init"
   Subroutine Element3D_Scal_Init(dElem,dCoord,QuadratureOrder,Element_Type)
      Type(Element3D_Scal)                   :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord
      PetscInt,Intent(IN)                    :: QuadratureOrder
      PetscInt,Intent(IN)                    :: Element_Type
      
      Select Case (Element_Type)
         Case (MEF90_P1_Lagrange)
            Call Element_P_Lagrange_3D_Scal_Init(dElem,dCoord,1,QuadratureOrder)

         Case (MEF90_P2_Lagrange)
            Call Element_P_Lagrange_3D_Scal_Init(dElem,dCoord,2,QuadratureOrder)

!         Case (MEF90_Q1_Lagrange)
!            Call Element_Q_Lagrange_3D_Scal_Init(dElem,dCoord,1,QuadratureOrder)
!         Case (MEF90_Q2_Lagrange)
!            Call Element_Q_Lagrange_3D_Scal_Init(dElem,dCoord,2,QuadratureOrder)
         Case Default
            Print*,__FUNCT__,': Element type not implemented yet',Element_Type
      End Select
   End Subroutine Element3D_Scal_Init                                

#undef __FUNCT__
#define __FUNCT__ "Element3D_InitAll"
!!!
!!! This routine does not use the values cached in dMeshtopology that are also directly accessible from the mesh
!!!
   Subroutine Element3D_InitAll(dMeshTopology,dElem,dQuadratureOrder)
      Type(MeshTopology_Type)                     :: dMeshTopology
      Type(Element3D),Dimension(:),Pointer        :: dElem
      PetscInt,Intent(IN)                         :: dQuadratureOrder
      
      PetscInt                                    :: set,setID,iELoc,iE,ierr
      PetscInt                                    :: numDim,numCells,numVertices
      PetscInt                                    :: numVertexinCell,point
      PetscInt,Dimension(:),Pointer               :: cone
      Type(IS)                                    :: CellIS,CellSetIS
      PetscInt,Dimension(:),Pointer               :: CellID,CellSetID
      PetscReal,Dimension(:,:),Pointer            :: Coords
      PetscReal,Dimension(:),Pointer              :: TmpCoords
      Type(SectionReal)                           :: CoordSection
      
      Call DMMeshGetDimension(dMeshTopology%mesh,numDim,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumSize(dMeshTopology%mesh,"height",0,numCells,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumSize(dMeshTopology%mesh,"depth",0,numVertices,ierr);CHKERRQ(ierr)

      Allocate(dElem(numCells))
      Call DMMeshGetLabelIdIS(dMeshTopology%mesh,'Cell Sets',CellSetIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetIS,CellSetID,ierr);CHKERRQ(ierr)

      !!! Initialize the Basis Functions in each element
      Call DMMeshGetSectionReal(dMeshTopology%mesh,'coordinates',CoordSection,ierr);CHKERRQ(ierr)
      Do_CellSet: Do set = 1,size(CellSetID)
         setID = CellSetID(set)
         Call DMMeshGetStratumIS(dMeshTopology%mesh,'Cell Sets',setId,CellIS,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(CellIS,CellID,ierr);CHKERRQ(ierr)

         !!! Get the number of vertices per cell for the cells in the set
         !!! We always assume that all cells in a set are of the same type so
         !!! we only need to do this once.
         !!! This is a bit crazy but will get much better when we switch to DMComplex
         Call DMMeshGetConeSize(dMeshTopology%mesh,CellID(1),numVertexinCell,ierr);CHKERRQ(ierr)    
         !!! The cone may contain faces and edges in the case of an interpolated mesh
         !!! so we need to filter out based on the vertex point number:
         Allocate(cone(numVertexinCell))
         numVertexinCell = 0
         Do point = 1, size(cone)
            If ( (cone(point) >= numCells) .AND. (cone(point) < numCells + numVertices) ) Then
               numVertexinCell = numVertexinCell + 1
            End If
         End Do
         DeAllocate(cone)
         
         Allocate(TmpCoords(numDim * numVertexinCell))
         Allocate(Coords   (numDim,  numVertexinCell))

         Do_Elem_iE: Do iELoc = 1,size(CellID)
            iE = CellID(iELoc)
            Call SectionRealRestrictClosure(CoordSection,dMeshTopology%mesh,iE-1,Size(TmpCoords),TmpCoords,ierr);CHKERRQ(ierr)
            Coords = Reshape(TmpCoords,(/numDim,numVertexinCell /) )
            !!! WTF? why not reshaping the arguments in Init_Element? 
            Call ElementInit(dElem(iE),Coords,dQuadratureOrder,dMeshTopology%cellSet(set)%ElemType)
         End Do Do_Elem_iE
         Call ISRestoreIndicesF90(CellIS,CellID,ierr);CHKERRQ(ierr)
         DeAllocate(TmpCoords)
         DeAllocate(Coords)
      End Do Do_CellSet
      Call SectionRealDestroy(CoordSection,ierr);CHKERRQ(ierr)
      Call ISRestoreIndicesF90(CellSetIS,CellSetID,ierr);CHKERRQ(ierr)
   End Subroutine Element3D_InitAll

#undef __FUNCT__
#define __FUNCT__ "Element3D_Init"
   Subroutine Element3D_Init(dElem,dCoord,QuadratureOrder,Element_Type)
      Type(Element3D)                        :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord
      PetscInt,Intent(IN)                    :: QuadratureOrder
      PetscInt,Intent(IN)                    :: Element_Type
      
      Select Case (Element_Type)
         Case (MEF90_P1_Lagrange)
            Call Element_P_Lagrange_3D_Init(dElem,dCoord,1,QuadratureOrder)

         Case (MEF90_P2_Lagrange)
            Call Element_P_Lagrange_3D_Init(dElem,dCoord,2,QuadratureOrder)

!         Case (MEF90_Q1_Lagrange)
!            Call Element_Q_Lagrange_3D_Init(dElem,dCoord,1,QuadratureOrder)
!         Case (MEF90_Q2_Lagrange)
!            Call Element_Q_Lagrange_3D_Init(dElem,dCoord,2,QuadratureOrder)
         Case Default
            Print*,__FUNCT__,': Element type not implemented yet',Element_Type
      End Select
   End Subroutine Element3D_Init                                

#undef __FUNCT__
#define __FUNCT__ "Element3D_Elast_InitAll"
!!!
!!! This routine does not use the values cached in dMeshtopology that are also directly accessible from the mesh
!!!
   Subroutine Element3D_Elast_InitAll(dMeshTopology,dElem,dQuadratureOrder)
      Type(MeshTopology_Type)                     :: dMeshTopology
      Type(Element3D_Elast),Dimension(:),Pointer  :: dElem
      PetscInt,Intent(IN)                         :: dQuadratureOrder
      
      PetscInt                                    :: set,setID,iELoc,iE,ierr
      PetscInt                                    :: numDim,numCells,numVertices
      PetscInt                                    :: numVertexinCell,point
      PetscInt,Dimension(:),Pointer               :: cone
      Type(IS)                                    :: CellIS,CellSetIS
      PetscInt,Dimension(:),Pointer               :: CellID,CellSetID
      PetscReal,Dimension(:,:),Pointer            :: Coords
      PetscReal,Dimension(:),Pointer              :: TmpCoords
      Type(SectionReal)                           :: CoordSection
      
      Call DMMeshGetDimension(dMeshTopology%mesh,numDim,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumSize(dMeshTopology%mesh,"height",0,numCells,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumSize(dMeshTopology%mesh,"depth",0,numVertices,ierr);CHKERRQ(ierr)

      Allocate(dElem(numCells))
      Call DMMeshGetLabelIdIS(dMeshTopology%mesh,'Cell Sets',CellSetIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetIS,CellSetID,ierr);CHKERRQ(ierr)

      !!! Initialize the Basis Functions in each element
      Call DMMeshGetSectionReal(dMeshTopology%mesh,'coordinates',CoordSection,ierr);CHKERRQ(ierr)
      Do_CellSet: Do set = 1,size(CellSetID)
         setID = CellSetID(set)
         Call DMMeshGetStratumIS(dMeshTopology%mesh,'Cell Sets',setId,CellIS,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(CellIS,CellID,ierr);CHKERRQ(ierr)

         !!! Get the number of vertices per cell for the cells in the set
         !!! We always assume that all cells in a set are of the same type so
         !!! we only need to do this once.
         !!! This is a bit crazy but will get much better when we switch to DMComplex
         Call DMMeshGetConeSize(dMeshTopology%mesh,CellID(1),numVertexinCell,ierr);CHKERRQ(ierr)    
         !!! The cone may contain faces and edges in the case of an interpolated mesh
         !!! so we need to filter out based on the vertex point number:
         Allocate(cone(numVertexinCell))
         numVertexinCell = 0
         Do point = 1, size(cone)
            If ( (cone(point) >= numCells) .AND. (cone(point) < numCells + numVertices) ) Then
               numVertexinCell = numVertexinCell + 1
            End If
         End Do
         DeAllocate(cone)
         
         Allocate(TmpCoords(numDim * numVertexinCell))
         Allocate(Coords   (numDim,  numVertexinCell))

         Do_Elem_iE: Do iELoc = 1,size(CellID)
            iE = CellID(iELoc)
            Call SectionRealRestrictClosure(CoordSection,dMeshTopology%mesh,iE-1,Size(TmpCoords),TmpCoords,ierr);CHKERRQ(ierr)
            Coords = Reshape(TmpCoords,(/numDim,numVertexinCell /) )
            !!! WTF? why not reshaping the arguments in Init_Element? 
            Call ElementInit(dElem(iE),Coords,dQuadratureOrder,dMeshTopology%cellSet(set)%ElemType)
         End Do Do_Elem_iE
         Call ISRestoreIndicesF90(CellIS,CellID,ierr);CHKERRQ(ierr)
         DeAllocate(TmpCoords)
         DeAllocate(Coords)
      End Do Do_CellSet
      Call SectionRealDestroy(CoordSection,ierr);CHKERRQ(ierr)
      Call ISRestoreIndicesF90(CellSetIS,CellSetID,ierr);CHKERRQ(ierr)
   End Subroutine Element3D_Elast_InitAll

#undef __FUNCT__
#define __FUNCT__ "Element3D_Elast_Init"
   Subroutine Element3D_Elast_Init(dElem,dCoord,QuadratureOrder,Element_Type)
      Type(Element3D_Elast)                  :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord
      PetscInt,Intent(IN)                    :: QuadratureOrder
      PetscInt,Intent(IN)                    :: Element_Type
      
      Select Case (Element_Type)
         Case (MEF90_P1_Lagrange)
            Call Element_P_Lagrange_3D_Elast_Init(dElem,dCoord,1,QuadratureOrder)

         Case (MEF90_P2_Lagrange)
            Call Element_P_Lagrange_3D_Elast_Init(dElem,dCoord,2,QuadratureOrder)

!         Case (MEF90_Q1_Lagrange)
!            Call Element_Q_Lagrange_3D_Elast_Init(dElem,dCoord,1,QuadratureOrder)
!         Case (MEF90_Q2_Lagrange)
!            Call Element_Q_Lagrange_3D_Elast_Init(dElem,dCoord,2,QuadratureOrder)
         Case Default
            Print*,__FUNCT__,': Element type not implemented yet',Element_Type
      End Select
   End Subroutine Element3D_Elast_Init                                

#undef __FUNCT__
#define __FUNCT__ "Element_P_Lagrange_2D_Scal_Init"
   Subroutine Element_P_Lagrange_2D_Scal_Init(dElem,dCoord,dPolynomialOrder,dQuadratureOrder)
      ! Compute the quadrature weights and the value of the basis functions and their gradient 
      ! at the quadrature points.
      ! One day when I am smart I will use FIAT for that...
      !
      ! Assumes that the elements connectivity is known
      !
      
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
         PhiHat(1,:) = 4.0_Kr * (1.0_Kr - Xi%X - Xi%Y) * Xi%X      
         PhiHat(2,:) = 4.0_Kr * Xi%X * Xi%Y     
         PhiHat(3,:) = 4.0_Kr * (1.0_Kr - Xi%X - Xi%Y) * Xi%Y     
         PhiHat(4,:) = 2.0_Kr * (1.0_Kr - Xi%X - Xi%Y)**2  - (1.0_Kr - Xi%X - Xi%Y)
         PhiHat(5,:) = 2.0_Kr * Xi%X**2 - Xi%X
         PhiHat(6,:) = 2.0_Kr * Xi%Y**2 - Xi%Y
         
         GradPhiHat(1,:)%X = 4.0_Kr * (1.0_Kr - 2.0_Kr * Xi%X - Xi%Y);GradPhiHat(1,:)%Y =-4.0_Kr * Xi%X
         GradPhiHat(2,:)%X = 4.0_Kr * Xi%Y;                           GradPhiHat(2,:)%Y = 4.0_Kr * Xi%X
         GradPhiHat(3,:)%X =-4.0_Kr * Xi%Y;                           GradPhiHat(3,:)%Y = 4.0_Kr * (1.0_Kr - Xi%X - 2.0_Kr * Xi%Y)
         GradPhiHat(4,:)%X = 1.0_Kr - 4.0_Kr * (1.0_Kr - Xi%X - Xi%Y);GradPhiHat(4,:)%Y = 1.0_Kr - 4.0_Kr * (1.0_Kr - Xi%X - Xi%Y)
         GradPhiHat(5,:)%X = 4.0_Kr * Xi%X - 1.0_Kr;                  GradPhiHat(5,:)%Y = 0.0_Kr
         GradPhiHat(6,:)%X = 0.0_Kr;                                  GradPhiHat(6,:)%Y = 4.0_Kr * Xi%Y - 1.0_Kr
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
#define __FUNCT__ "Element_P_Lagrange_2D_Init"
   Subroutine Element_P_Lagrange_2D_Init(dElem,dCoord,dPolynomialOrder,dQuadratureOrder)
      Type(Element2D)                        :: dElem
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
   End Subroutine Element_P_Lagrange_2D_Init

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
#define __FUNCT__ "Element_P_Lagrange_3D_Scal_Init"
   Subroutine Element_P_Lagrange_3D_Scal_Init(dElem,dCoord,dPolynomialOrder,dQuadratureOrder)
      ! Compute the quadrature weights and the value of the basis functions and their gradient 
      ! at the quadrature points.
      ! One day when I am smart I will use FIAT for that...
      !
      ! Assumes that the elements connectivity is known
      !
      
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
#define __FUNCT__ "Element_P_Lagrange_3D_Init"
   Subroutine Element_P_Lagrange_3D_Init(dElem,dCoord,dPolynomialOrder,dQuadratureOrder)
      Type(Element3D)                        :: dElem
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
   End Subroutine Element_P_Lagrange_3D_Init

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
!      If (Associated(dElem%ID_DoF)) Then
!         DeAllocate(dElem%ID_DoF)
!      End If
   End Subroutine Element2D_Scal_Destroy

#undef __FUNCT__
#define __FUNCT__ "Element2D_Destroy"
   Subroutine Element2D_Destroy(dElem)
      Type(Element2D)                        :: dElem
      If (Associated(dElem%BF)) Then
         DeAllocate(dElem%BF)
      End If
      If (Associated(dElem%Der_BF)) Then
         DeAllocate(dElem%Der_BF)
      End If
      If (Associated(dElem%Gauss_C)) Then
         DeAllocate(dElem%Gauss_C)
      End If
!      If (Associated(dElem%ID_DoF)) Then
!         DeAllocate(dElem%ID_DoF)
!      End If
   End Subroutine Element2D_Destroy
   
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
!      If (Associated(dElem%ID_DoF)) Then
!         DeAllocate(dElem%ID_DoF)
!      End If
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
!      If (Associated(dElem%ID_DoF)) Then
!         DeAllocate(dElem%ID_DoF)
!      End If
   End Subroutine Element3D_Scal_Destroy

#undef __FUNCT__
#define __FUNCT__ "Element3D_Destroy"
   Subroutine Element3D_Destroy(dElem)
      Type(Element3D)                        :: dElem
      If (Associated(dElem%BF)) Then
         DeAllocate(dElem%BF)
      End If
      If (Associated(dElem%Der_BF)) Then
         DeAllocate(dElem%Der_BF)
      End If
      If (Associated(dElem%Gauss_C)) Then
         DeAllocate(dElem%Gauss_C)
      End If
!      If (Associated(dElem%ID_DoF)) Then
!         DeAllocate(dElem%ID_DoF)
!      End If
   End Subroutine Element3D_Destroy
   
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
!      If (Associated(dElem%ID_DoF)) Then
!         DeAllocate(dElem%ID_DoF)
!      End If
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
#define __FUNCT__ "Element2D_View"
   Subroutine Element2D_View(dElem,viewer)
      Type(Element2D)                                 :: dElem
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
   End Subroutine Element2D_View
   
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
   Subroutine Element3D_View(dElem,viewer)
      Type(Element3D)                                 :: dElem
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
   End Subroutine Element3D_View
   
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
