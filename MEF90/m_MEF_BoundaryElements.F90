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
#define __FUNCT__ "BoundaryElement2D_Init"
   Subroutine BoundaryElement2D_Init(dElem,dCoord,QuadratureOrder,Element_Type)
      Type(BoundaryElement2D)                :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord
      PetscInt,Intent(IN)                    :: QuadratureOrder
      PetscInt,Intent(IN)                    :: Element_Type
      
      Select Case (Element_Type)
         Case (MEF90_P1_Lagrange)
            Call BoundaryElement_P_Lagrange_2D_Init(dElem,dCoord,1,QuadratureOrder,Element_Type)

         Case (MEF90_P2_Lagrange)
            Call BoundaryElement_P_Lagrange_2D_Init(dElem,dCoord,2,QuadratureOrder,Element_Type)

         !!!Case (MEF90_Q1_Lagrange)
         !!!   Call BoundaryElement_Q_Lagrange_2D_Init(dElem,dCoord,1,QuadratureOrder)
         !!!Case (MEF90_Q2_Lagrange)
         !!!   Call BoundaryElement_Q_Lagrange_2D_Init(dElem,dCoord,2,QuadratureOrder)
         Case Default
            Print*,__FUNCT__,'Element type not implemented yet',Element_Type
      End Select
   End Subroutine BoundaryElement2D_Init         
   
#undef __FUNCT__
#define __FUNCT__ "BoundaryElement_P_Lagrange_2D_Init"
   Subroutine BoundaryElement_P_Lagrange_2D_Init(dElem,dCoord,dPolynomialOrder,dQuadratureOrder,Element_Type)
      Type(BoundaryElement2D)                :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord      ! coord(i,j)=ith coord of jth vertice
      PetscInt                               :: dPolynomialOrder,dQuadratureOrder
      PetscInt,Intent(IN)                    :: Element_Type
      
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
      
      Call ElementInit(tmpElem,tmpCoord,dQuadratureOrder,Element_Type)
      Select Case (Element_Type)
         Case (MEF90_P1_Lagrange)
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

         !Case (MEF90_P2_Lagrange)
         Case Default
            Print*,__FUNCT__,'Element type not implemented yet',Element_Type
      End Select


      Call ElementDestroy(tmpElem)
      deAllocate(tmpCoord)
      
   End Subroutine BoundaryElement_P_Lagrange_2D_Init                          
   
!!! 
!!! 3D
!!!
#undef __FUNCT__
#define __FUNCT__ "BoundaryElement3D_Init"
   Subroutine BoundaryElement3D_Init(dElem,dCoord,QuadratureOrder,Element_Type)
      Type(BoundaryElement3D)                :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord
      PetscInt,Intent(IN)                    :: QuadratureOrder
      PetscInt,Intent(IN)                    :: Element_Type
      
      Select Case (Element_Type)
         Case (MEF90_P1_Lagrange,MEF90_P2_Lagrange)
            Call BoundaryElement_P_Lagrange_3D_Init(dElem,dCoord,QuadratureOrder,Element_Type)

         !!!Case (MEF90_Q1_Lagrange,MEF90_Q2_Lagrange)
         !!!   Call BoundaryElement_Q_Lagrange_3D_Init(dElem,dCoord,QuadratureOrder,Element_Type)
         Case Default
            Print*,__FUNCT__,'Element type not implemented yet',Element_Type
      End Select
   End Subroutine BoundaryElement3D_Init         
   

#undef __FUNCT__
#define __FUNCT__ "BoundaryElement_P_Lagrange_3D_Init"
   Subroutine BoundaryElement_P_Lagrange_3D_Init(dElem,dCoord,dQuadratureOrder,Element_Type)
      Type(BoundaryElement3D)                :: dElem
      PetscReal,Dimension(:,:),Pointer       :: dCoord      ! coord(i,j)=ith coord of jth vertice
      PetscInt,Intent(IN)                    :: dQuadratureOrder
      PetscInt,Intent(IN)                    :: Element_Type
      
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
      
            
      Call ElementInit(tmpElem,tmpCoord,dQuadratureOrder,Element_Type)
      Select Case (Element_Type)
         Case (MEF90_P1_Lagrange)
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

         !Case (MEF90_P2_Lagrange)
            !!! I need to think about DoF ordering in this case... 
            !!! This is going to work out the same way. The mid-edge dof with 
            !!! xi3>0 are linear combinations of the other dof  
         Case Default
            Print*,__FUNCT__,'Element type not implemented yet',Element_Type
      End Select


      Call ElementDestroy(tmpElem)
      deAllocate(tmpCoord)
   End Subroutine BoundaryElement_P_Lagrange_3D_Init
   
End Module m_MEF_BoundaryElements

