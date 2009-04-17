Module m_MEF_Elements

#include "finclude/petscdef.h"
#include "finclude/petscmeshdef.h"
#include "finclude/petscviewerdef.h"

   Use m_MEF_LinAlg
   Use m_MEF_Types
   Use m_MEF_Utils
   Use petsc
   Use petscmesh
      
   IMPLICIT NONE
   Private
   
   Public :: ElementInit
   Public :: ElementDestroy
   Public :: ElementView
   Public :: Init_Elem_Blk_Type


   Interface ElementInit
      Module Procedure Init_Element2D_Scal, Init_Element2D, Init_Element2D_Elast, Init_Element3D_Scal, Init_Element3D, Init_Element3D_Elast, Init_AllElement2D_Scal, Init_AllElement2D, Init_AllElement2D_Elast, Init_AllElement3D_Scal, Init_AllElement3D, Init_AllElement3D_Elast
   End Interface ElementInit
   
   Interface ElementDestroy
      Module Procedure Destroy_Element2D_Scal, Destroy_Element2D, Destroy_Element2D_Elast, Destroy_Element3D_Scal, Destroy_Element3D, Destroy_Element3D_Elast
   End Interface ElementDestroy
   
   Interface ElementView
      Module Procedure Element2D_ScalView, Element2D_ScalPtrView, Element2DView, Element2DPtrView, Element2D_ElastView, Element2D_ElastPtrView, Element3D_ScalView, Element3D_ScalPtrView, Element3DView, Element3DPtrView, Element3D_ElastView, Element3D_ElastPtrView
   End Interface ElementView
   
   PetscInt, Parameter, Public                   :: MEF90_P1_Lagrange = 1
   PetscInt, Parameter, Public                   :: MEF90_P2_Lagrange = 2
   
!   PetscInt, Parameter, Public                   :: MEF90_Q1_Lagrange = 3
!   PetscInt, Parameter, Public                   :: MEF90_Q2_Lagrange = 4

 Contains
   Subroutine Init_Elem_Blk_Type(dBlk, dDim)
      Type (Elem_Blk_Type)                   :: dBlk
      PetscInt                               :: dDim
      
      Select Case (dDim)
      Case (2)
         Select Case (dBlk%Elem_Type)
         Case (MEF90_P1_Lagrange)         
            dBlk%DoF_Location = (/ 0, 0, 0, 3/)
            dBlk%Num_Face = 0
            dBlk%Num_Edge = 3
            dBlk%Num_Vert = 3
         Case (MEF90_P2_Lagrange)
            dBlk%DoF_Location = (/ 0, 0, 3, 3/)
            dBlk%Num_Face = 0
            dBlk%Num_Edge = 3
            dBlk%Num_Vert = 3
         Case Default
            Print*, 'Unknown element type', dBlk%Elem_Type
            STOP
         End Select
      Case (3)
         Select Case (dBlk%Elem_Type)
         Case (MEF90_P1_Lagrange)         
            dBlk%DoF_Location = (/ 0, 0, 0, 4/)
            dBlk%Num_Face = 4
            dBlk%Num_Edge = 6
            dBlk%Num_Vert = 4
         Case (MEF90_P2_Lagrange)
            dBlk%DoF_Location = (/ 0, 0, 6, 4/)
            dBlk%Num_Face = 4
            dBlk%Num_Edge = 6
            dBlk%Num_Vert = 4
         Case Default
            Print*, 'Unknown element type', dBlk%Elem_Type
            STOP
         End Select
      Case Default
         Print*, 'Unknown dimension', dDim
         STOP
      End Select
      dBlk%Num_DoF = sum(dBlk%DoF_Location)

   End Subroutine Init_Elem_Blk_Type     
   
   Subroutine Init_AllElement2D_Scal(dMeshTopology, dElem, dQuadratureOrder)
      Type(MeshTopology_Type)                     :: dMeshTopology
      Type(Element2D_Scal), Dimension(:), Pointer :: dElem
      PetscInt, Intent(IN)                        :: dQuadratureOrder
      
      PetscInt                                    :: iBlk, iELoc, iE, iErr
      PetscReal, Dimension(:,:), Pointer          :: Coords
      PetscReal, Dimension(:), Pointer            :: TmpCoords
      Type(SectionReal)                           :: CoordSection
      
      Allocate(dElem(dMeshTopology%Num_Elems))
      !!! Initialize the Basis Functions in each element
      Call MeshGetSectionReal(dMeshTopology%mesh, 'coordinates', CoordSection, iErr); CHKERRQ(iErr)
      Do_Elem_iBlk: Do iBlk = 1, dMeshTopology%Num_Elem_Blks
         Allocate(TmpCoords(dMeshTopology%Num_Dim * dMeshTopology%Elem_Blk(iBlk)%Num_Vert))
         Allocate(Coords   (dMeshTopology%Num_Dim,  dMeshTopology%Elem_Blk(iBlk)%Num_Vert))

         Do_Elem_iE: Do iELoc = 1, dMeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = dMeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            Call MeshRestrictClosure(dMeshTopology%mesh, CoordSection, iE-1, Size(TmpCoords), TmpCoords, iErr); CHKERRQ(iErr)
             Coords = Reshape(TmpCoords, (/dMeshTopology%Num_Dim, dMeshTopology%Elem_Blk(iBlk)%Num_Vert /) )
             !!! WTF? why not reshaping the arguments in Init_Element? 
            Call ElementInit(dElem(iE), Coords, dQuadratureOrder, dMeshTopology%Elem_Blk(iBlk)%Elem_Type)
         End Do Do_Elem_iE
         DeAllocate(TmpCoords)
         DeAllocate(Coords)
      End Do Do_Elem_iBlk
      Call SectionRealDestroy(CoordSection, iErr); CHKERRQ(iErr)
   End Subroutine Init_AllElement2D_Scal
      
   
   Subroutine Init_Element2D_Scal(dElem, dCoord, QuadratureOrder, Element_Type)
      Type (Element2D_Scal)                  :: dElem
      PetscReal, Dimension(:,:), Pointer     :: dCoord
      PetscInt, Intent(IN)                   :: QuadratureOrder
      PetscInt, Intent(IN)                   :: Element_Type
      
      Select Case (Element_Type)
         Case (MEF90_P1_Lagrange)
            Call Init_Element_P_Lagrange_2D_Scal(dElem, dCoord, 1, QuadratureOrder)

         Case (MEF90_P2_Lagrange)
            Call Init_Element_P_Lagrange_2D_Scal(dElem, dCoord, 2, QuadratureOrder)

!         Case (MEF90_Q1_Lagrange)
!            Call Init_Element_Q_Lagrange_2D_Scal(dElem, dCoord, 1, QuadratureOrder)
!         Case (MEF90_Q2_Lagrange)
!            Call Init_Element_Q_Lagrange_2D_Scal(dElem, dCoord, 2, QuadratureOrder)
         Case Default
            Print*, 'Element type not implemented yet', Element_Type
      End Select
   End Subroutine Init_Element2D_Scal                                
   
   Subroutine Init_AllElement2D(dMeshTopology, dElem, dQuadratureOrder)
      Type(MeshTopology_Type)                     :: dMeshTopology
      Type(Element2D), Dimension(:), Pointer      :: dElem
      PetscInt, Intent(IN)                        :: dQuadratureOrder
      
      PetscInt                                    :: iBlk, iELoc, iE, iErr
      PetscReal, Dimension(:,:), Pointer          :: Coords
      PetscReal, Dimension(:), Pointer            :: TmpCoords
      Type(SectionReal)                           :: CoordSection
      
      Allocate(dElem(dMeshTopology%Num_Elems))
      !!! Initialize the Basis Functions in each element
      Call MeshGetSectionReal(dMeshTopology%mesh, 'coordinates', CoordSection, iErr); CHKERRQ(iErr)
      Do_Elem_iBlk: Do iBlk = 1, dMeshTopology%Num_Elem_Blks
         Allocate(TmpCoords(dMeshTopology%Num_Dim * dMeshTopology%Elem_Blk(iBlk)%Num_Vert))
         Allocate(Coords   (dMeshTopology%Num_Dim,  dMeshTopology%Elem_Blk(iBlk)%Num_Vert))

         Do_Elem_iE: Do iELoc = 1, dMeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = dMeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            Call MeshRestrictClosure(dMeshTopology%mesh, CoordSection, iE-1, Size(TmpCoords), TmpCoords, iErr); CHKERRQ(iErr)
             Coords = Reshape(TmpCoords, (/dMeshTopology%Num_Dim, dMeshTopology%Elem_Blk(iBlk)%Num_Vert /) )
             !!! WTF? why not reshaping the arguments in Init_Element? 
            Call ElementInit(dElem(iE), Coords, dQuadratureOrder, dMeshTopology%Elem_Blk(iBlk)%Elem_Type)
         End Do Do_Elem_iE
         DeAllocate(TmpCoords)
         DeAllocate(Coords)
      End Do Do_Elem_iBlk
      Call SectionRealDestroy(CoordSection, iErr); CHKERRQ(iErr)
   End Subroutine Init_AllElement2D

   Subroutine Init_Element2D(dElem, dCoord, QuadratureOrder, Element_Type)
      Type (Element2D)                       :: dElem
      PetscReal, Dimension(:,:), Pointer     :: dCoord
      PetscInt, Intent(IN)                   :: QuadratureOrder
      PetscInt, Intent(IN)                   :: Element_Type
      
      Select Case (Element_Type)
         Case (MEF90_P1_Lagrange)
            Call Init_Element_P_Lagrange_2D(dElem, dCoord, 1, QuadratureOrder)

         Case (MEF90_P2_Lagrange)
            Call Init_Element_P_Lagrange_2D(dElem, dCoord, 2, QuadratureOrder)

!         Case (MEF90_Q1_Lagrange)
!            Call Init_Element_Q_Lagrange_2D_Scal(dElem, dCoord, 1, QuadratureOrder)
!         Case (MEF90_Q2_Lagrange)
!            Call Init_Element_Q_Lagrange_2D_Scal(dElem, dCoord, 2, QuadratureOrder)
         Case Default
            Print*, 'Element type not implemented yet', Element_Type
      End Select
   End Subroutine Init_Element2D                                

   Subroutine Init_AllElement2D_Elast(dMeshTopology, dElem, dQuadratureOrder)
      Type(MeshTopology_Type)                      :: dMeshTopology
      Type(Element2D_Elast), Dimension(:), Pointer :: dElem
      PetscInt, Intent(IN)                         :: dQuadratureOrder
      
      PetscInt                                     :: iBlk, iELoc, iE, iErr
      PetscReal, Dimension(:,:), Pointer           :: Coords
      PetscReal, Dimension(:), Pointer             :: TmpCoords
      Type(SectionReal)                            :: CoordSection
      
      Allocate(dElem(dMeshTopology%Num_Elems))
      !!! Initialize the Basis Functions in each element
      Call MeshGetSectionReal(dMeshTopology%mesh, 'coordinates', CoordSection, iErr); CHKERRQ(iErr)
      Do_Elem_iBlk: Do iBlk = 1, dMeshTopology%Num_Elem_Blks
         Allocate(TmpCoords(dMeshTopology%Num_Dim * dMeshTopology%Elem_Blk(iBlk)%Num_Vert))
         Allocate(Coords   (dMeshTopology%Num_Dim,  dMeshTopology%Elem_Blk(iBlk)%Num_Vert))

         Do_Elem_iE: Do iELoc = 1, dMeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = dMeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            Call MeshRestrictClosure(dMeshTopology%mesh, CoordSection, iE-1, Size(TmpCoords), TmpCoords, iErr); CHKERRQ(iErr)
             Coords = Reshape(TmpCoords, (/dMeshTopology%Num_Dim, dMeshTopology%Elem_Blk(iBlk)%Num_Vert /) )
             !!! WTF? why not reshaping the arguments in Init_Element? 
            Call ElementInit(dElem(iE), Coords, dQuadratureOrder, dMeshTopology%Elem_Blk(iBlk)%Elem_Type)
         End Do Do_Elem_iE
         DeAllocate(TmpCoords)
         DeAllocate(Coords)
      End Do Do_Elem_iBlk
      Call SectionRealDestroy(CoordSection, iErr); CHKERRQ(iErr)
   End Subroutine Init_AllElement2D_Elast

   Subroutine Init_Element2D_Elast(dElem, dCoord, QuadratureOrder, Element_Type)
      Type (Element2D_Elast)                 :: dElem
      PetscReal, Dimension(:,:), Pointer     :: dCoord
      PetscInt, Intent(IN)                   :: QuadratureOrder
      PetscInt, Intent(IN)                   :: Element_Type
      
      Select Case (Element_Type)
         Case (MEF90_P1_Lagrange)
            Call Init_Element_P_Lagrange_2D_Elast(dElem, dCoord, 1, QuadratureOrder)

         Case (MEF90_P2_Lagrange)
            Call Init_Element_P_Lagrange_2D_Elast(dElem, dCoord, 2, QuadratureOrder)

!         Case (MEF90_Q1_Lagrange)
!            Call Init_Element_Q_Lagrange_2D_Scal(dElem, dCoord, 1, QuadratureOrder)
!         Case (MEF90_Q2_Lagrange)
!            Call Init_Element_Q_Lagrange_2D_Scal(dElem, dCoord, 2, QuadratureOrder)
         Case Default
            Print*, 'Element type not implemented yet', Element_Type
      End Select
   End Subroutine Init_Element2D_Elast


   Subroutine Init_AllElement3D_Scal(dMeshTopology, dElem, dQuadratureOrder)
      Type(MeshTopology_Type)                      :: dMeshTopology
      Type(Element3D_Scal), Dimension(:), Pointer  :: dElem
      PetscInt, Intent(IN)                         :: dQuadratureOrder
      
      PetscInt                                     :: iBlk, iELoc, iE, iErr
      PetscReal, Dimension(:,:), Pointer           :: Coords
      PetscReal, Dimension(:), Pointer             :: TmpCoords
      Type(SectionReal)                            :: CoordSection
      
      Allocate(dElem(dMeshTopology%Num_Elems))
      !!! Initialize the Basis Functions in each element
      Call MeshGetSectionReal(dMeshTopology%mesh, 'coordinates', CoordSection, iErr); CHKERRQ(iErr)
      Do_Elem_iBlk: Do iBlk = 1, dMeshTopology%Num_Elem_Blks
         Allocate(TmpCoords(dMeshTopology%Num_Dim * dMeshTopology%Elem_Blk(iBlk)%Num_Vert))
         Allocate(Coords   (dMeshTopology%Num_Dim,  dMeshTopology%Elem_Blk(iBlk)%Num_Vert))

         Do_Elem_iE: Do iELoc = 1, dMeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = dMeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            Call MeshRestrictClosure(dMeshTopology%mesh, CoordSection, iE-1, Size(TmpCoords), TmpCoords, iErr); CHKERRQ(iErr)
             Coords = Reshape(TmpCoords, (/dMeshTopology%Num_Dim, dMeshTopology%Elem_Blk(iBlk)%Num_Vert /) )
             !!! WTF? why not reshaping the arguments in Init_Element? 
            Call ElementInit(dElem(iE), Coords, dQuadratureOrder, dMeshTopology%Elem_Blk(iBlk)%Elem_Type)
         End Do Do_Elem_iE
         DeAllocate(TmpCoords)
         DeAllocate(Coords)
      End Do Do_Elem_iBlk
      Call SectionRealDestroy(CoordSection, iErr); CHKERRQ(iErr)
   End Subroutine Init_AllElement3D_Scal

   Subroutine Init_Element3D_Scal(dElem, dCoord, QuadratureOrder, Element_Type)
      Type (Element3D_Scal)                  :: dElem
      PetscReal, Dimension(:,:), Pointer     :: dCoord
      PetscInt, Intent(IN)                   :: QuadratureOrder
      PetscInt, Intent(IN)                   :: Element_Type
      
      Select Case (Element_Type)
         Case (MEF90_P1_Lagrange)
            Call Init_Element_P_Lagrange_3D_Scal(dElem, dCoord, 1, QuadratureOrder)

         Case (MEF90_P2_Lagrange)
            Call Init_Element_P_Lagrange_3D_Scal(dElem, dCoord, 2, QuadratureOrder)

!         Case (MEF90_Q1_Lagrange)
!            Call Init_Element_Q_Lagrange_3D_Scal(dElem, dCoord, 1, QuadratureOrder)
!         Case (MEF90_Q2_Lagrange)
!            Call Init_Element_Q_Lagrange_3D_Scal(dElem, dCoord, 2, QuadratureOrder)
         Case Default
            Print*, 'Element type not implemented yet', Element_Type
      End Select
   End Subroutine Init_Element3D_Scal                                

   Subroutine Init_AllElement3D(dMeshTopology, dElem, dQuadratureOrder)
      Type(MeshTopology_Type)                      :: dMeshTopology
      Type(Element3D), Dimension(:), Pointer       :: dElem
      PetscInt, Intent(IN)                         :: dQuadratureOrder
      
      PetscInt                                     :: iBlk, iELoc, iE, iErr
      PetscReal, Dimension(:,:), Pointer           :: Coords
      PetscReal, Dimension(:), Pointer             :: TmpCoords
      Type(SectionReal)                            :: CoordSection
      
      Allocate(dElem(dMeshTopology%Num_Elems))
      !!! Initialize the Basis Functions in each element
      Call MeshGetSectionReal(dMeshTopology%mesh, 'coordinates', CoordSection, iErr); CHKERRQ(iErr)
      Do_Elem_iBlk: Do iBlk = 1, dMeshTopology%Num_Elem_Blks
         Allocate(TmpCoords(dMeshTopology%Num_Dim * dMeshTopology%Elem_Blk(iBlk)%Num_Vert))
         Allocate(Coords   (dMeshTopology%Num_Dim,  dMeshTopology%Elem_Blk(iBlk)%Num_Vert))

         Do_Elem_iE: Do iELoc = 1, dMeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = dMeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            Call MeshRestrictClosure(dMeshTopology%mesh, CoordSection, iE-1, Size(TmpCoords), TmpCoords, iErr); CHKERRQ(iErr)
             Coords = Reshape(TmpCoords, (/dMeshTopology%Num_Dim, dMeshTopology%Elem_Blk(iBlk)%Num_Vert /) )
             !!! WTF? why not reshaping the arguments in Init_Element? 
            Call ElementInit(dElem(iE), Coords, dQuadratureOrder, dMeshTopology%Elem_Blk(iBlk)%Elem_Type)
         End Do Do_Elem_iE
         DeAllocate(TmpCoords)
         DeAllocate(Coords)
      End Do Do_Elem_iBlk
      Call SectionRealDestroy(CoordSection, iErr); CHKERRQ(iErr)
   End Subroutine Init_AllElement3D

   Subroutine Init_Element3D(dElem, dCoord, QuadratureOrder, Element_Type)
      Type (Element3D)                       :: dElem
      PetscReal, Dimension(:,:), Pointer     :: dCoord
      PetscInt, Intent(IN)                   :: QuadratureOrder
      PetscInt, Intent(IN)                   :: Element_Type
      
      Select Case (Element_Type)
         Case (MEF90_P1_Lagrange)
            Call Init_Element_P_Lagrange_3D(dElem, dCoord, 1, QuadratureOrder)

         Case (MEF90_P2_Lagrange)
            Call Init_Element_P_Lagrange_3D(dElem, dCoord, 2, QuadratureOrder)

!         Case (MEF90_Q1_Lagrange)
!            Call Init_Element_Q_Lagrange_3D(dElem, dCoord, 1, QuadratureOrder)
!         Case (MEF90_Q2_Lagrange)
!            Call Init_Element_Q_Lagrange_3D(dElem, dCoord, 2, QuadratureOrder)
         Case Default
            Print*, 'Element type not implemented yet', Element_Type
      End Select
   End Subroutine Init_Element3D                                

   Subroutine Init_AllElement3D_Elast(dMeshTopology, dElem, dQuadratureOrder)
      Type(MeshTopology_Type)                      :: dMeshTopology
      Type(Element3D_Elast), Dimension(:), Pointer :: dElem
      PetscInt, Intent(IN)                         :: dQuadratureOrder
      
      PetscInt                                     :: iBlk, iELoc, iE, iErr
      PetscReal, Dimension(:,:), Pointer           :: Coords
      PetscReal, Dimension(:), Pointer             :: TmpCoords
      Type(SectionReal)                            :: CoordSection
      
      Allocate(dElem(dMeshTopology%Num_Elems))
      !!! Initialize the Basis Functions in each element
      Call MeshGetSectionReal(dMeshTopology%mesh, 'coordinates', CoordSection, iErr); CHKERRQ(iErr)
      Do_Elem_iBlk: Do iBlk = 1, dMeshTopology%Num_Elem_Blks
         Allocate(TmpCoords(dMeshTopology%Num_Dim * dMeshTopology%Elem_Blk(iBlk)%Num_Vert))
         Allocate(Coords   (dMeshTopology%Num_Dim,  dMeshTopology%Elem_Blk(iBlk)%Num_Vert))

         Do_Elem_iE: Do iELoc = 1, dMeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = dMeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            Call MeshRestrictClosure(dMeshTopology%mesh, CoordSection, iE-1, Size(TmpCoords), TmpCoords, iErr); CHKERRQ(iErr)
             Coords = Reshape(TmpCoords, (/dMeshTopology%Num_Dim, dMeshTopology%Elem_Blk(iBlk)%Num_Vert /) )
             !!! WTF? why not reshaping the arguments in Init_Element? 
            Call ElementInit(dElem(iE), Coords, dQuadratureOrder, dMeshTopology%Elem_Blk(iBlk)%Elem_Type)
         End Do Do_Elem_iE
         DeAllocate(TmpCoords)
         DeAllocate(Coords)
      End Do Do_Elem_iBlk
      Call SectionRealDestroy(CoordSection, iErr); CHKERRQ(iErr)
   End Subroutine Init_AllElement3D_Elast

   Subroutine Init_Element3D_Elast(dElem, dCoord, QuadratureOrder, Element_Type)
      Type (Element3D_Elast)                 :: dElem
      PetscReal, Dimension(:,:), Pointer     :: dCoord
      PetscInt, Intent(IN)                   :: QuadratureOrder
      PetscInt, Intent(IN)                   :: Element_Type
      
      Select Case (Element_Type)
         Case (MEF90_P1_Lagrange)
            Call Init_Element_P_Lagrange_3D_Elast(dElem, dCoord, 1, QuadratureOrder)

         Case (MEF90_P2_Lagrange)
            Call Init_Element_P_Lagrange_3D_Elast(dElem, dCoord, 2, QuadratureOrder)

!         Case (MEF90_Q1_Lagrange)
!            Call Init_Element_Q_Lagrange_3D_Elast(dElem, dCoord, 1, QuadratureOrder)
!         Case (MEF90_Q2_Lagrange)
!            Call Init_Element_Q_Lagrange_3D_Elast(dElem, dCoord, 2, QuadratureOrder)
         Case Default
            Print*, 'Element type not implemented yet', Element_Type
      End Select
   End Subroutine Init_Element3D_Elast                                

   Subroutine Init_Element_P_Lagrange_2D_Scal(dElem, dCoord, dPolynomialOrder, dQuadratureOrder)
      ! Compute the quadrature weights and the value of the basis functions and their gradient 
      ! at the quadrature points.
      ! One day when I am smart I will use FIAT for that...
      !
      ! Assumes that the elements connectivity is known
      !
      
      Type (Element2D_Scal)                  :: dElem
      PetscReal, Dimension(:,:), Pointer     :: dCoord      ! coord(i,j)=ith coord of jth vertice
      PetscInt                               :: dPolynomialOrder, dQuadratureOrder
      
      PetscInt                               :: Nb_Gauss
      PetscInt                               :: Num_Dof
      PetscInt                               :: iDoF, iG
      Type (Mat2D)                           :: Bt          ! The transposed of transformation matrix
      PetscReal                              :: DetBinv     ! The determinant of B^{-1}

      PetscReal, Dimension(:,:), Pointer     :: PhiHat      ! PhiHat(i,k) The value of the ith basis function at the kth integration point
      Type (Vect2D), Dimension(:,:), Pointer :: GradPhiHat
      
      
      Type (Vect2D), Dimension(:), Pointer   :: Xi ! The quadrature points coordinates in the reference element
      
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
         Xi(1) = (/ 1.0_Kr / 6.0_Kr, 1.0_Kr / 6.0_Kr /)
         Xi(2) = (/ 2.0_Kr / 3.0_Kr, 1.0_Kr / 6.0_Kr /)
         Xi(3) = (/ 1.0_Kr / 6.0_Kr, 2.0_Kr / 3.0_Kr /)
         dElem%Gauss_C = detBinv / 6.0_Kr

      Case(3)
         Nb_Gauss = 4
         Allocate(Xi(Nb_Gauss))
         Allocate(dElem%Gauss_C(Nb_Gauss))
         dElem%Gauss_C    =  detBinv * 25.0_Kr / 96.0_Kr
         dElem%Gauss_C(1) = -detBinv * 9.0_Kr / 32.0_Kr
         Xi(1) = (/ 1.0_Kr / 3.0_Kr, 1.0_Kr / 3.0_Kr /)
         Xi(2) = (/ 3.0_Kr / 5.0_Kr, 1.0_Kr / 5.0_Kr /)
         Xi(3) = (/ 1.0_Kr / 5.0_Kr, 3.0_Kr / 5.0_Kr /)
         Xi(4) = (/ 1.0_Kr / 5.0_Kr, 1.0_Kr / 5.0_Kr /)
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
         Xi(1) = (/ 0.0_Kr         , 0.0_Kr          /)
         Xi(2) = (/ 1.0_Kr / 2.0_Kr, 0.0_Kr          /)
         Xi(3) = (/ 1.0_Kr         , 0.0_Kr          /)
         Xi(4) = (/ 1.0_Kr / 2.0_Kr, 1.0_Kr / 2.0_Kr /)
         Xi(5) = (/ 0.0_Kr         , 1.0_Kr          /)
         Xi(6) = (/ 0.0_Kr         , 1.0_Kr / 2.0_Kr /)
         Xi(7) = (/ 1.0_Kr / 3.0_Kr, 1.0_Kr / 3.0_Kr /)
         
      Case Default
         Print*, 'Unimplemented quadrature order', dQuadratureOrder
         STOP
      End Select
      
      Select Case (dPolynomialOrder)
      Case(1)
         Num_DoF = 3
         Allocate(PhiHat(Num_DoF, Nb_Gauss))
         Allocate(GradPhiHat(Num_DoF, Nb_Gauss))
         PhiHat(1,:) = 1.0_Kr - Xi%X - Xi%Y
         PhiHat(2,:) = Xi(:)%X
         PhiHat(3,:) = Xi(:)%Y
         
         GradPhiHat(1,:)%X = -1.0_Kr; GradPhiHat(1,:)%Y = -1.0_Kr 
         GradPhiHat(2,:)%X =  1.0_Kr; GradPhiHat(2,:)%Y =  0.0_Kr 
         GradPhiHat(3,:)%X =  0.0_Kr; GradPhiHat(3,:)%Y =  1.0_Kr
          
      Case(2)
         Num_DoF = 6
         Allocate(PhiHat(Num_DoF, Nb_Gauss))
         Allocate(GradPhiHat(Num_DoF, Nb_Gauss))
         PhiHat(1,:) = 4.0_Kr * (1.0_Kr - Xi%X - Xi%Y) * Xi%X      
         PhiHat(2,:) = 4.0_Kr * Xi%X * Xi%Y     
         PhiHat(3,:) = 4.0_Kr * (1.0_Kr - Xi%X - Xi%Y) * Xi%Y     
         PhiHat(4,:) = 2.0_Kr * (1.0_Kr - Xi%X - Xi%Y)**2  - (1.0_Kr - Xi%X - Xi%Y)
         PhiHat(5,:) = 2.0_Kr * Xi%X**2 - Xi%X
         PhiHat(6,:) = 2.0_Kr * Xi%Y**2 - Xi%Y
         
         GradPhiHat(1,:)%X = 4.0_Kr * (1.0_Kr - 2.0_Kr * Xi%X - Xi%Y); GradPhiHat(1,:)%Y =-4.0_Kr * Xi%X
         GradPhiHat(2,:)%X = 4.0_Kr * Xi%Y;                            GradPhiHat(2,:)%Y = 4.0_Kr * Xi%X
         GradPhiHat(3,:)%X =-4.0_Kr * Xi%Y;                            GradPhiHat(3,:)%Y = 4.0_Kr * (1.0_Kr - Xi%X - 2.0_Kr * Xi%Y)
         GradPhiHat(4,:)%X = 1.0_Kr - 4.0_Kr * (1.0_Kr - Xi%X - Xi%Y); GradPhiHat(4,:)%Y = 1.0_Kr - 4.0_Kr * (1.0_Kr - Xi%X - Xi%Y)
         GradPhiHat(5,:)%X = 4.0_Kr * Xi%X - 1.0_Kr;                   GradPhiHat(5,:)%Y = 0.0_Kr
         GradPhiHat(6,:)%X = 0.0_Kr;                                   GradPhiHat(6,:)%Y = 4.0_Kr * Xi%Y - 1.0_Kr
      Case Default
         Print*, 'Unimplemented PolynomialOrder', dPolynomialOrder
      End Select
      
      Allocate (dElem%BF(Num_DoF, Nb_Gauss)) 
      Allocate (dElem%Grad_BF(Num_DoF, Nb_Gauss))
      dElem%BF = PhiHat
      Do iDoF = 1, Num_DoF
         Do iG = 1, Nb_Gauss
            dElem%Grad_BF(iDoF, iG) = Bt * GradPhiHat(iDoF, iG) 
         End Do
      End Do
     
      DeAllocate(Xi)
      DeAllocate(PhiHat)
      DeAllocate(GradPhiHat)
   End Subroutine Init_Element_P_Lagrange_2D_Scal
   
   Subroutine Init_Element_P_Lagrange_2D(dElem, dCoord, dPolynomialOrder, dQuadratureOrder)
      Type (Element2D)                       :: dElem
      PetscReal, Dimension(:,:), Pointer     :: dCoord      ! coord(i,j)=ith coord of jth vertice
      PetscInt                               :: dPolynomialOrder, dQuadratureOrder
   
      Type (Element2D_Scal)                  :: Elem_Scal
      PetscInt                               :: dim = 2 
      PetscInt                               :: Num_DoF, Nb_Gauss, i
      
      
      Call Init_Element_P_Lagrange_2D_Scal(Elem_Scal, dCoord, dPolynomialOrder, dQuadratureOrder)
      Num_DoF  = Size(Elem_Scal%BF, 1) 
      Nb_Gauss = Size(Elem_Scal%BF, 2)
      Allocate(dElem%Gauss_C(Nb_Gauss))
      Allocate(dElem%BF(Num_DoF * dim, Nb_Gauss))
      Allocate(dElem%Der_BF(Num_DoF * dim, Nb_Gauss))
         
      dElem%Gauss_C = Elem_Scal%Gauss_C
      dElem%BF(:,:)%X = 0.0_Kr
      dElem%BF(:,:)%Y = 0.0_Kr
      dElem%Der_BF(:,:)%XX = 0.0_Kr
      dElem%Der_BF(:,:)%XY = 0.0_Kr
      dElem%Der_BF(:,:)%YX = 0.0_Kr
      dElem%Der_BF(:,:)%YY = 0.0_Kr
      
      Do i = 0, Num_DoF-1
         dElem%BF(i*dim+1,:)%X = Elem_Scal%BF(i+1,:)
         dElem%BF(i*dim+2,:)%Y = Elem_Scal%BF(i+1,:)
         dElem%Der_BF(i*dim+1,:)%XX = Elem_Scal%Grad_BF(i+1,:)%X
         dElem%Der_BF(i*dim+1,:)%XY = Elem_Scal%Grad_BF(i+1,:)%Y
         dElem%Der_BF(i*dim+2,:)%YX = Elem_Scal%Grad_BF(i+1,:)%X
         dElem%Der_BF(i*dim+2,:)%YY = Elem_Scal%Grad_BF(i+1,:)%Y
      End Do
      Call ElementDestroy(Elem_Scal)
   End Subroutine Init_Element_P_Lagrange_2D

   Subroutine Init_Element_P_Lagrange_2D_Elast(dElem, dCoord, dPolynomialOrder, dQuadratureOrder)
      Type (Element2D_Elast)                 :: dElem
      PetscReal, Dimension(:,:), Pointer     :: dCoord      ! coord(i,j)=ith coord of jth vertice
      PetscInt                               :: dPolynomialOrder, dQuadratureOrder
   
      Type (Element2D_Scal)                  :: Elem_Scal
      PetscInt                               :: dim = 2 
      PetscInt                               :: Num_DoF, Nb_Gauss, i
      
      
      Call Init_Element_P_Lagrange_2D_Scal(Elem_Scal, dCoord, dPolynomialOrder, dQuadratureOrder)
      Num_DoF   = Size(Elem_Scal%BF, 1) 
      Nb_Gauss = Size(Elem_Scal%BF, 2)
      Allocate(dElem%Gauss_C(Nb_Gauss))
      Allocate(dElem%BF(Num_DoF * dim, Nb_Gauss))
      Allocate(dElem%GradS_BF(Num_DoF * dim, Nb_Gauss))
         
      dElem%Gauss_C = Elem_Scal%Gauss_C
      dElem%BF(:,:)%X = 0.0_Kr
      dElem%BF(:,:)%Y = 0.0_Kr
      dElem%GradS_BF(:,:)%XX = 0.0_Kr
      dElem%GradS_BF(:,:)%XY = 0.0_Kr
      dElem%GradS_BF(:,:)%YY = 0.0_Kr
      Do i = 0, Num_DoF-1
         dElem%BF(i*dim+1,:)%X = Elem_Scal%BF(i+1,:)
         dElem%BF(i*dim+2,:)%Y = Elem_Scal%BF(i+1,:)
         dElem%GradS_BF(i*dim+1,:)%XX = Elem_Scal%Grad_BF(i+1,:)%X
         dElem%GradS_BF(i*dim+1,:)%XY = Elem_Scal%Grad_BF(i+1,:)%Y / 2.0_Kr
         dElem%GradS_BF(i*dim+2,:)%XY = Elem_Scal%Grad_BF(i+1,:)%X / 2.0_Kr
         dElem%GradS_BF(i*dim+2,:)%YY = Elem_Scal%Grad_BF(i+1,:)%Y
      End Do
      Call ElementDestroy(Elem_Scal)
   End Subroutine Init_Element_P_Lagrange_2D_Elast

   Subroutine Init_Element_P_Lagrange_3D_Scal(dElem, dCoord, dPolynomialOrder, dQuadratureOrder)
      ! Compute the quadrature weights and the value of the basis functions and their gradient 
      ! at the quadrature points.
      ! One day when I am smart I will use FIAT for that...
      !
      ! Assumes that the elements connectivity is known
      !
      
      Type (Element3D_Scal)                  :: dElem
      PetscReal, Dimension(:,:), Pointer     :: dCoord      ! coord(i,j)=ith coord of jth vertice
      PetscInt                               :: dPolynomialOrder, dQuadratureOrder
      
      PetscInt                               :: Nb_Gauss
      PetscInt                               :: Num_Dof
      PetscInt                               :: iDoF, iG
      Type (Mat3D)                           :: Bt          ! The transposed of transformation matrix
      PetscReal                              :: DetBinv     ! The determinant of B^{-1}

      PetscReal, Dimension(:,:), Pointer     :: PhiHat      ! PhiHat(i,k) The value of the ith basis function at the kth integration point
      Type (Vect3D), Dimension(:,:), Pointer :: GradPhiHat
      
      
      Type (Vect3D), Dimension(:), Pointer   :: Xi          ! The quadrature points coordinates in the reference element
      
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
      Case(1,2,3)
         !!! 3rd order cubature on a tetrahedron forula from J.E. Akin' book
         Nb_Gauss = 5
         Allocate(Xi(Nb_Gauss))
         Allocate(dElem%Gauss_C(Nb_Gauss))
         Xi(1) = (/ .25_Kr, .25_Kr, .25_Kr /)
         Xi(2) = (/ .5_Kr, 1._Kr / 6._Kr, 1._Kr / 6._Kr /)
         Xi(3) = (/ 1._Kr / 6._Kr, .5_Kr, 1._Kr / 6._Kr /)
         Xi(4) = (/ 1._Kr / 6._Kr, 1._Kr / 6._Kr, .5_Kr /)
         Xi(5) = (/ 1._Kr / 6._Kr, 1._Kr / 6._Kr, 1._Kr / 6._Kr /)
         dElem%Gauss_C(1) = -2.0_Kr / 15.0_Kr * detBinv
         dElem%Gauss_C(2:5) = 3.0_Kr / 40.0_Kr * detBinv

      Case Default
         Print*, 'Unimplemented quadrature order', dQuadratureOrder
         STOP
      End Select
      
      Select Case (dPolynomialOrder)
      Case(1)
         Num_DoF = 4
         Allocate(PhiHat(Num_DoF, Nb_Gauss))
         Allocate(GradPhiHat(Num_DoF, Nb_Gauss))
         PhiHat(1,:) = 1.0_Kr - Xi%X - Xi%Y - Xi%Z
         PhiHat(2,:) = Xi(:)%X
         PhiHat(3,:) = Xi(:)%Y
         PhiHat(4,:) = Xi(:)%Z
         
         GradPhiHat(1,:)%X = -1.0_Kr; GradPhiHat(1,:)%Y = -1.0_Kr; GradPhiHat(1,:)%Z = -1.0_Kr; 
         GradPhiHat(2,:)%X =  1.0_Kr; GradPhiHat(2,:)%Y =  0.0_Kr; GradPhiHat(2,:)%Z =  0.0_Kr; 
         GradPhiHat(3,:)%X =  0.0_Kr; GradPhiHat(3,:)%Y =  1.0_Kr; GradPhiHat(3,:)%Z =  0.0_Kr;
         GradPhiHat(4,:)%X =  0.0_Kr; GradPhiHat(4,:)%Y =  0.0_Kr; GradPhiHat(4,:)%Z =  1.0_Kr;
          
      Case Default
         Print*, 'Unimplemented PolynomialOrder', dPolynomialOrder
      End Select
      
      Allocate (dElem%BF(Num_DoF, Nb_Gauss)) 
      Allocate (dElem%Grad_BF(Num_DoF, Nb_Gauss))
      dElem%BF = PhiHat
      Do iDoF = 1, Num_DoF
         Do iG = 1, Nb_Gauss
            dElem%Grad_BF(iDoF, iG) = Bt * GradPhiHat(iDoF, iG) 
         End Do
      End Do
     
      DeAllocate(Xi)
      DeAllocate(PhiHat)
      DeAllocate(GradPhiHat)
   End Subroutine Init_Element_P_Lagrange_3D_Scal

   Subroutine Init_Element_P_Lagrange_3D(dElem, dCoord, dPolynomialOrder, dQuadratureOrder)
      Type (Element3D)                       :: dElem
      PetscReal, Dimension(:,:), Pointer     :: dCoord      ! coord(i,j)=ith coord of jth vertice
      PetscInt                               :: dPolynomialOrder, dQuadratureOrder
   
      Type (Element3D_Scal)                  :: Elem_Scal
      PetscInt                               :: dim = 3 
      PetscInt                               :: Num_DoF, Nb_Gauss, i
      
      
      Call Init_Element_P_Lagrange_3D_Scal(Elem_Scal, dCoord, dPolynomialOrder, dQuadratureOrder)
      Num_DoF   = Size(Elem_Scal%BF, 1) 
      Nb_Gauss = Size(Elem_Scal%BF, 2)
      Allocate(dElem%Gauss_C(Nb_Gauss))
      Allocate(dElem%BF(Num_DoF * dim, Nb_Gauss))
      Allocate(dElem%Der_BF(Num_DoF * dim, Nb_Gauss))
         
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
      
      Do i = 0, Num_DoF-1
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
   End Subroutine Init_Element_P_Lagrange_3D

   Subroutine Init_Element_P_Lagrange_3D_Elast(dElem, dCoord, dPolynomialOrder, dQuadratureOrder)
      Type (Element3D_Elast)                 :: dElem
      PetscReal, Dimension(:,:), Pointer     :: dCoord      ! coord(i,j)=ith coord of jth vertice
      PetscInt                               :: dPolynomialOrder, dQuadratureOrder
   
      Type (Element3D_Scal)                  :: Elem_Scal
      PetscInt                               :: dim = 3 
      PetscInt                               :: Num_DoF, Nb_Gauss, i
      
      
      Call Init_Element_P_Lagrange_3D_Scal(Elem_Scal, dCoord, dPolynomialOrder, dQuadratureOrder)
      Num_DoF  = Size(Elem_Scal%BF, 1) 
      Nb_Gauss = Size(Elem_Scal%BF, 2)
      Allocate(dElem%Gauss_C(Nb_Gauss))
      Allocate(dElem%BF(Num_DoF * dim, Nb_Gauss))
      Allocate(dElem%GradS_BF(Num_DoF * dim, Nb_Gauss))
         
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
      
      Do i = 0, Num_DoF-1
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
   End Subroutine Init_Element_P_Lagrange_3D_Elast

   Subroutine Destroy_Element2D_Scal(dElem)
      Type (Element2D_Scal)                  :: dElem
      If (Associated(dElem%BF)) Then
         DeAllocate(dElem%BF)
      End If
      If (Associated(dElem%Grad_BF)) Then
         DeAllocate(dElem%Grad_BF)
      End If
      If (Associated(dElem%Gauss_C)) Then
         DeAllocate(dElem%Gauss_C)
      End If
      If (Associated(dElem%ID_DoF)) Then
         DeAllocate(dElem%ID_DoF)
      End If
   End Subroutine Destroy_Element2D_Scal

   Subroutine Destroy_Element2D(dElem)
      Type (Element2D)                       :: dElem
      If (Associated(dElem%BF)) Then
         DeAllocate(dElem%BF)
      End If
      If (Associated(dElem%Der_BF)) Then
         DeAllocate(dElem%Der_BF)
      End If
      If (Associated(dElem%Gauss_C)) Then
         DeAllocate(dElem%Gauss_C)
      End If
      If (Associated(dElem%ID_DoF)) Then
         DeAllocate(dElem%ID_DoF)
      End If
   End Subroutine Destroy_Element2D
   
   Subroutine Destroy_Element2D_Elast(dElem)
      Type (Element2D_Elast)                 :: dElem
      If (Associated(dElem%BF)) Then
         DeAllocate(dElem%BF)
      End If
      If (Associated(dElem%GradS_BF)) Then
         DeAllocate(dElem%GradS_BF)
      End If
      If (Associated(dElem%Gauss_C)) Then
         DeAllocate(dElem%Gauss_C)
      End If
      If (Associated(dElem%ID_DoF)) Then
         DeAllocate(dElem%ID_DoF)
      End If
   End Subroutine Destroy_Element2D_Elast
   
   Subroutine Destroy_Element3D_Scal(dElem)
      Type (Element3D_Scal)                  :: dElem
      If (Associated(dElem%BF)) Then
         DeAllocate(dElem%BF)
      End If
      If (Associated(dElem%Grad_BF)) Then
         DeAllocate(dElem%Grad_BF)
      End If
      If (Associated(dElem%Gauss_C)) Then
         DeAllocate(dElem%Gauss_C)
      End If
      If (Associated(dElem%ID_DoF)) Then
         DeAllocate(dElem%ID_DoF)
      End If
   End Subroutine Destroy_Element3D_Scal

   Subroutine Destroy_Element3D(dElem)
      Type (Element3D)                       :: dElem
      If (Associated(dElem%BF)) Then
         DeAllocate(dElem%BF)
      End If
      If (Associated(dElem%Der_BF)) Then
         DeAllocate(dElem%Der_BF)
      End If
      If (Associated(dElem%Gauss_C)) Then
         DeAllocate(dElem%Gauss_C)
      End If
      If (Associated(dElem%ID_DoF)) Then
         DeAllocate(dElem%ID_DoF)
      End If
   End Subroutine Destroy_Element3D
   
   Subroutine Destroy_Element3D_Elast(dElem)
      Type (Element3D_Elast)                 :: dElem
      If (Associated(dElem%BF)) Then
         DeAllocate(dElem%BF)
      End If
      If (Associated(dElem%GradS_BF)) Then
         DeAllocate(dElem%GradS_BF)
      End If
      If (Associated(dElem%Gauss_C)) Then
         DeAllocate(dElem%Gauss_C)
      End If
      If (Associated(dElem%ID_DoF)) Then
         DeAllocate(dElem%ID_DoF)
      End If
   End Subroutine Destroy_Element3D_Elast

!!!
!!! ELEMENT VIEWERS
!!!

   Subroutine Element2D_ScalView(dElem, viewer)
      Type(Element2D_Scal)                           :: dElem
      Type(PetscViewer)                              :: viewer
      
      PetscInt                                       :: Nb_Gauss, Nb_DoF, iDoF, iG, iErr
      Character(len=512)                             :: CharBuffer
                
      Nb_DoF   = size(dElem%BF,1)
      Nb_Gauss = size(dElem%BF,2)
      Write(CharBuffer, 102) Nb_DoF
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 103) Nb_Gauss
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)

      Nb_Gauss = Size(dElem%BF,2)
      Do iDoF = 1, Nb_DoF
         Write(CharBuffer, 104) iDoF
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)

         Write(CharBuffer, 200) '        BF \n    'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%BF(iDoF, iG)
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)

         Write(CharBuffer, 200) '        Grad_BF (X,Y)\n    'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%Grad_BF(iDoF, iG)%X
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n    'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%Grad_BF(iDoF, iG)%Y
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      End Do
101 Format('*** Element  ', I9, '\n'c)
102 Format('    Nb_DoF   ', I9, '\n'c)
103 Format('    Nb_Gauss ', I9, '\n'c)
104 Format('    *** DoF  ', I9, '\n'c)
200 Format(A)
201 Format('   ', F5.2)
   End Subroutine Element2D_ScalView
   
   Subroutine Element2D_ScalPtrView(dElems, viewer)
      Type(Element2D_Scal), Dimension(:), Pointer    :: dElems
      Type(PetscViewer)                              :: viewer

      PetscInt                                       :: iE, iErr
      Character(len=512)                             :: CharBuffer
      
      Write(CharBuffer, 100) Size(dElems)
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Do iE = 1, Size(dElems)
         Write(CharBuffer, 101) iE
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         
         Call ElementView(dElems(iE), viewer)
      End Do
100 Format('    Number of elements =================== ', I9, '\n'c)
101 Format('*** Element  ', I9, '\n'c)
   End Subroutine Element2D_ScalPtrView

   Subroutine Element2DView(dElem, viewer)
      Type(Element2D)                                :: dElem
      Type(PetscViewer)                              :: viewer
      
      PetscInt                                       :: Nb_Gauss, Nb_DoF, iDoF, iG, iErr
      Character(len=512)                             :: CharBuffer

      Nb_DoF   = size(dElem%BF,1)
      Nb_Gauss = size(dElem%BF,2)
      Write(CharBuffer, 102) Nb_DoF
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 103) Nb_Gauss
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)

      Nb_Gauss = Size(dElem%BF,2)
      Do iDoF = 1, Nb_DoF
         Write(CharBuffer, 104) iDoF
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)

         Write(CharBuffer, 200) '        BF (X,Y)\n    'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%BF(iDoF, iG)%X
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n   'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%BF(iDoF, iG)%Y
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)

         Write(CharBuffer, 200) '        Der_BF (XX, XY, YX, YY)\n    'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%Der_BF(iDoF, iG)%XX
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n    'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%Der_BF(iDoF, iG)%XY
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%Der_BF(iDoF, iG)%YX
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n    'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%Der_BF(iDoF, iG)%YY
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n    'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      End Do
101 Format('*** Element  ', I9, '\n'c)
102 Format('    Nb_DoF   ', I9, '\n'c)
103 Format('    Nb_Gauss ', I9, '\n'c)
104 Format('    *** DoF  ', I9, '\n'c)
200 Format(A)
201 Format('   ', F5.2)
   End Subroutine Element2DView
   
   Subroutine Element2DPtrView(dElems, viewer)
      Type(Element2D), Dimension(:), Pointer         :: dElems
      Type(PetscViewer)                              :: viewer

      PetscInt                                       :: iE, iErr
      Character(len=512)                             :: CharBuffer
      
      Write(CharBuffer, 100) Size(dElems)
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Do iE = 1, Size(dElems)
         Write(CharBuffer, 101) iE
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         
         Call ElementView(dElems(iE), viewer)
      End Do
100 Format('    Number of elements =================== ', I9, '\n'c)
101 Format('*** Element  ', I9, '\n'c)
   End Subroutine Element2DPtrView

   Subroutine Element2D_ElastView(dElem, viewer)
      Type(Element2D_Elast)                          :: dElem
      Type(PetscViewer)                              :: viewer
      
      PetscInt                                       :: Nb_Gauss, Nb_DoF, iDoF, iG, iErr
      Character(len=512)                             :: CharBuffer

      Nb_DoF   = size(dElem%BF,1)
      Nb_Gauss = size(dElem%BF,2)
      Write(CharBuffer, 102) Nb_DoF
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 103) Nb_Gauss
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)

      Nb_Gauss = Size(dElem%BF,2)
      Do iDoF = 1, Nb_DoF
         Write(CharBuffer, 104) iDoF
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)

         Write(CharBuffer, 200) '        BF (X,Y)\n    'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%BF(iDoF, iG)%X
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n   'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%BF(iDoF, iG)%Y
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)

         Write(CharBuffer, 200) '        GradS_BF (XX, YY, XY)\n    'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%GradS_BF(iDoF, iG)%XX
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n    'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%GradS_BF(iDoF, iG)%YY
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%GradS_BF(iDoF, iG)%XY
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n    'c
      End Do
101 Format('*** Element  ', I9, '\n'c)
102 Format('    Nb_DoF   ', I9, '\n'c)
103 Format('    Nb_Gauss ', I9, '\n'c)
104 Format('    *** DoF  ', I9, '\n'c)
200 Format(A)
201 Format('   ', F5.2)
   End Subroutine Element2D_ElastView
   
   Subroutine Element2D_ElastPtrView(dElems, viewer)
      Type(Element2D_Elast), Dimension(:), Pointer   :: dElems
      Type(PetscViewer)                              :: viewer

      PetscInt                                       :: iE, iErr
      Character(len=512)                             :: CharBuffer
      
      Write(CharBuffer, 100) Size(dElems)
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Do iE = 1, Size(dElems)
         Write(CharBuffer, 101) iE
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         
         Call ElementView(dElems(iE), viewer)
      End Do
100 Format('    Number of elements =================== ', I9, '\n'c)
101 Format('*** Element  ', I9, '\n'c)
   End Subroutine Element2D_ElastPtrView

   Subroutine Element3D_ScalView(dElem, viewer)
      Type(Element3D_Scal)                           :: dElem
      Type(PetscViewer)                              :: viewer
      
      PetscInt                                       :: Nb_Gauss, Nb_DoF, iDoF, iG, iErr
      Character(len=512)                             :: CharBuffer
                
      Nb_DoF   = size(dElem%BF,1)
      Nb_Gauss = size(dElem%BF,2)
      Write(CharBuffer, 102) Nb_DoF
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 103) Nb_Gauss
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)

      Nb_Gauss = Size(dElem%BF,2)
      Do iDoF = 1, Nb_DoF
         Write(CharBuffer, 104) iDoF
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)

         Write(CharBuffer, 200) '        BF \n    'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%BF(iDoF, iG)
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)

         Write(CharBuffer, 200) '        Grad_BF (X,Y,Z)\n    'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%Grad_BF(iDoF, iG)%X
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n    'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%Grad_BF(iDoF, iG)%Y
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n    'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%Grad_BF(iDoF, iG)%Z
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n    'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)         
      End Do
101 Format('*** Element  ', I9, '\n'c)
102 Format('    Nb_DoF   ', I9, '\n'c)
103 Format('    Nb_Gauss ', I9, '\n'c)
104 Format('    *** DoF  ', I9, '\n'c)
200 Format(A)
201 Format('   ', F5.2)
   End Subroutine Element3D_ScalView
   
   Subroutine Element3D_ScalPtrView(dElems, viewer)
      Type(Element3D_Scal), Dimension(:), Pointer    :: dElems
      Type(PetscViewer)                              :: viewer

      PetscInt                                       :: iE, iErr
      Character(len=512)                             :: CharBuffer
      
      Write(CharBuffer, 100) Size(dElems)
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Do iE = 1, Size(dElems)
         Write(CharBuffer, 101) iE
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         
!! MGK: Something is broken here
!         Call ElementView(dElems(iE), viewer)
      End Do
100 Format('    Number of elements =================== ', I9, '\n'c)
101 Format('*** Element  ', I9, '\n'c)
   End Subroutine Element3D_ScalPtrView

   Subroutine Element3DView(dElem, viewer)
      Type(Element3D)                                :: dElem
      Type(PetscViewer)                              :: viewer
      
      PetscInt                                       :: Nb_Gauss, Nb_DoF, iDoF, iG, iErr
      Character(len=512)                             :: CharBuffer

      Nb_DoF   = size(dElem%BF,1)
      Nb_Gauss = size(dElem%BF,2)
      Write(CharBuffer, 102) Nb_DoF
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 103) Nb_Gauss
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)

      Nb_Gauss = Size(dElem%BF,2)
      Do iDoF = 1, Nb_DoF
         Write(CharBuffer, 104) iDoF
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)

         Write(CharBuffer, 200) '        BF (X,Y,Z)\n    'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%BF(iDoF, iG)%X
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n   'c
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%BF(iDoF, iG)%Y
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n   'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%BF(iDoF, iG)%Z
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)

         Write(CharBuffer, 200) '        Der_BF (XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ)\n    'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%Der_BF(iDoF, iG)%XX
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n    'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%Der_BF(iDoF, iG)%XY
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%Der_BF(iDoF, iG)%XZ
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%Der_BF(iDoF, iG)%YX
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n    'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%Der_BF(iDoF, iG)%YY
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n    'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%Der_BF(iDoF, iG)%YZ
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%Der_BF(iDoF, iG)%ZX
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%Der_BF(iDoF, iG)%ZY
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%Der_BF(iDoF, iG)%ZZ
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n'c
      End Do
101 Format('*** Element  ', I9, '\n'c)
102 Format('    Nb_DoF   ', I9, '\n'c)
103 Format('    Nb_Gauss ', I9, '\n'c)
104 Format('    *** DoF  ', I9, '\n'c)
200 Format(A)
201 Format('   ', F5.2)
   End Subroutine Element3DView
   
   Subroutine Element3DPtrView(dElems, viewer)
      Type(Element3D), Dimension(:), Pointer         :: dElems
      Type(PetscViewer)                              :: viewer

      PetscInt                                       :: iE, iErr
      Character(len=512)                             :: CharBuffer
      
      Write(CharBuffer, 100) Size(dElems)
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Do iE = 1, Size(dElems)
         Write(CharBuffer, 101) iE
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         
!! MGK: Something is broken here
!         Call ElementView(dElems(iE), viewer)
      End Do
100 Format('    Number of elements =================== ', I9, '\n'c)
101 Format('*** Element  ', I9, '\n'c)
   End Subroutine Element3DPtrView

   Subroutine Element3D_ElastView(dElem, viewer)
      Type(Element3D_Elast)                          :: dElem
      Type(PetscViewer)                              :: viewer
      
      PetscInt                                       :: Nb_Gauss, Nb_DoF, iDoF, iG, iErr
      Character(len=512)                             :: CharBuffer

      Nb_DoF   = size(dElem%BF,1)
      Nb_Gauss = size(dElem%BF,2)
      Write(CharBuffer, 102) Nb_DoF
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 103) Nb_Gauss
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)

      Nb_Gauss = Size(dElem%BF,2)
      Do iDoF = 1, Nb_DoF
         Write(CharBuffer, 104) iDoF
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)

         Write(CharBuffer, 200) '        BF (X,Y,Z)\n    'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%BF(iDoF, iG)%X
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n   'c
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%BF(iDoF, iG)%Y
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n   'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%BF(iDoF, iG)%Z
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)

         Write(CharBuffer, 200) '        GradS_BF (XX, YY, ZZ, YZ, XZ, XY)\n    'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%GradS_BF(iDoF, iG)%XX
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n    'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%GradS_BF(iDoF, iG)%YY
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n'c
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%GradS_BF(iDoF, iG)%ZZ
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n    'c
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%GradS_BF(iDoF, iG)%YZ
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n    'c
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%GradS_BF(iDoF, iG)%XZ
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n    'c
         Do iG = 1, Nb_Gauss
            Write(CharBuffer, 201) dElem%GradS_BF(iDoF, iG)%XY
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, *) '\n    'c
      End Do
101 Format('*** Element  ', I9, '\n'c)
102 Format('    Nb_DoF   ', I9, '\n'c)
103 Format('    Nb_Gauss ', I9, '\n'c)
104 Format('    *** DoF  ', I9, '\n'c)
200 Format(A)
201 Format('   ', F5.2)
   End Subroutine Element3D_ElastView
   
   Subroutine Element3D_ElastPtrView(dElems, viewer)
      Type(Element3D_Elast), Dimension(:), Pointer   :: dElems
      Type(PetscViewer)                              :: viewer

      PetscInt                                       :: iE, iErr
      Character(len=512)                             :: CharBuffer
      
      Write(CharBuffer, 100) Size(dElems)
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Do iE = 1, Size(dElems)
         Write(CharBuffer, 101) iE
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         
!! MGK: Something is broken here
!         Call ElementView(dElems(iE), viewer)
      End Do
100 Format('    Number of elements =================== ', I9, '\n'c)
101 Format('*** Element  ', I9, '\n'c)
   End Subroutine Element3D_ElastPtrView
End Module m_MEF_Elements
