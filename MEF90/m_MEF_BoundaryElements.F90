Module m_MEF_BoundaryElements
#include "finclude/petscdef.h"

   Use m_MEF_LinAlg
   Use m_MEF_Types
   Use m_MEF_Utils
   Use m_MEF_Elements
   Use petsc
      
   IMPLICIT NONE
   Private
   
   Public :: BoundaryElement_Init
   Public :: BoundaryElement_Destroy
   Public :: BoundaryElement_View
   Public :: Init_Side_Set_Type


   Interface BoundaryElement_Init
      Module Procedure BoundaryElement2D_Init,BoundaryElement3D_Init
   End Interface BoundaryElement_Init
   
   Interface BoundaryElement_Destroy
      Module Procedure BoundaryElement2D_Destroy,BoundaryElement3D_Destroy
   End Interface BoundaryElement_Destroy
   
   Interface BoundaryElement_View
      Module Procedure BoundaryElement2D_View,BoundaryElement3D_View
   End Interface BoundaryElement_View

Contains
!!!
!!! 2D
!!!
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
      tmpCoord(1,3) = dCoord(1,1) - NormalVector%X 
      tmpCoord(2,3) = dCoord(2,1) - NormalVector%Y
      
      Call Element_P_Lagrange_2D_Scal_Init(tmpElem,tmpCoord,dQuadratureOrder)
      Select Case (dPolynomialOrder)
         Case (1)
            Num_DoF  = 2
            Num_Gauss = size(tmpElem%BF,2)
            Allocate(dElem%Gauss_C(Num_Gauss))
            Allocate(dElem%BF(Num_DoF,Num_Gauss))
            dElem%Gauss_C = tmpElem%Gauss_C * 2.0_Kr
            Do iDoF = 1,Num_doF
               Do iG = 1,Num_Gauss
                  dElem%BF(iDoF,iG) = tmpElem%BF(iDoF,iG) + tmpElem%BF(Num_DoF+1,iG) / 2.0_Kr
               End Do
            End Do
         Case Default
            Print*,'[ERROR]: Polynomial order ',dPolynomialOrder,' not implemented in ',__FUNCT__
      End Select
      Call Element2D_Scal_Destroy(tmpElem)
      deAllocate(tmpCoord)
   End Subroutine Element_P_Lagrange_2DBoundary_Scal_Init                          
   
#undef __FUNCT__
#define __FUNCT__ "Element_P_Lagrange_3DBoundary_Scal_Init"
   Subroutine Element_P_Lagrange_3DBoundary_Scal_Init(dElem,dCoord,dPolynomialOrder,dQuadratureOrder)
      Type(BoundaryElement3D)                :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord      ! coord(i,j)=ith coord of jth vertice
      PetscInt,Intent(IN)                    :: dPolynomialOrder,dQuadratureOrder
      
      Type(Element3D_Scal)                   :: tmpElem
      PetscReal,Dimension(:,:),Pointer       :: tmpCoord
      PetscInt                               :: i,j,iDoF,iG,Num_Gauss,Num_DoF
      Type(Vect3D)                           :: Edge1,Edge2,NormalVector
      
      !!!
      !!! Create a bogus tet element by adding a 4th vertex along the normal of the
      !!! face,at distance XXX so that if a function is constant along the normal
      !!! direction of the face,one has \int_face fdx = \int_tet fdx 
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
                  dElem%BF(iDoF,iG) = (tmpElem%BF(iDoF,iG) + tmpElem%BF(Num_DoF+1,iG) / 3.0_Kr) * NormalVector
               End Do
            End Do
            !dElem%BF(3,:) = dElem%BF(3,:) + tmpElem%BF(4,:)

         !Case (2)
            !!! I need to think about DoF ordering in this case... 
            !!! This is going to work out the same way. The mid-edge dof with 
            !!! xi3>0 are linear combinations of the other dof  
         Case Default
            Print*,'[ERROR]: Polynomial order ',dPolynomialOrder,' not implemented in ',__FUNCT__
      End Select
      Call Element3D_Scal_Destroy(tmpElem)
      deAllocate(tmpCoord)
   End Subroutine Element_P_Lagrange_3DBoundary_Scal_Init
   
End Module m_MEF_BoundaryElements

