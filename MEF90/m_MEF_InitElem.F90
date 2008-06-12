Module m_MEF_InitElem

   Use m_AlgebLin
   Use m_MEF_Types
   Use m_Utils
   
   IMPLICIT NONE
   Private
   
#include "include/finclude/petsc.h"

   Public :: Init_Element
   Public :: Init_Element_2D_Scal
!   Public :: Destroy_Element


   Public :: Init_Element_P_Lagrange_2D_Scal

   Interface Init_Element
      Module Procedure Init_Element_2D_Scal
   End Interface Init_Element
      
 !  Interface Destroy_Element
 !     Module Procedure Destroy_Element2D_Scal
 !  End Interface Destroy_Element
 
   Integer, Parameter, Public                    :: MEF90_P1_Lagrange_2D_Scal = 1
!   Integer, Parameter, Public                    :: MEF90_P1_Lagrange_3D_Scal = 2
!   Integer, Parameter, Public                    :: MEF90_P1_Lagrange_2D_Vect = 3
!   Integer, Parameter, Public                    :: MEF90_P1_Lagrange_3D_Vect = 4
   
   Integer, Parameter, Public                    :: MEF90_P2_Lagrange_2D_Scal = 5
!   Integer, Parameter, Public                    :: MEF90_P2_Lagrange_3D_Scal = 6
!   Integer, Parameter, Public                    :: MEF90_P2_Lagrange_2D_Vect = 7
!   Integer, Parameter, Public                    :: MEF90_P2_Lagrange_3D_Vect = 8
   
!   Integer, Parameter, Public                    :: MEF90_Q1_Lagrange_2D_Scal = 9
!   Integer, Parameter, Public                    :: MEF90_Q1_Lagrange_3D_Scal = 10
!   Integer, Parameter, Public                    :: MEF90_Q1_Lagrange_2D_Vect = 11
!   Integer, Parameter, Public                    :: MEF90_Q1_Lagrange_3D_Vect = 12

!   Integer, Parameter, Public                    :: MEF90_Q2_Lagrange_2D_Scal = 13
!   Integer, Parameter, Public                    :: MEF90_Q2_Lagrange_3D_Scal = 14
!   Integer, Parameter, Public                    :: MEF90_Q2_Lagrange_2D_Vect = 15
!   Integer, Parameter, Public                    :: MEF90_Q2_Lagrange_3D_Vect = 16

 Contains
   Subroutine Init_Element_2D_Scal(dElem, dCoord, QuadratureOrder, Element_Type)
      Type (Element2D_Scal)                  :: dElem
      PetscReal, Dimension(:,:), Pointer     :: dCoord
      Integer, Intent(IN)                    :: QuadratureOrder
      Integer, Intent(IN)                    :: Element_Type
      
      Select Case (Element_Type)
         Case (MEF90_P1_Lagrange_2D_Scal)
            Call Init_Element_P_Lagrange_2D_Scal(dElem, dCoord, 1, QuadratureOrder)

         Case (MEF90_P2_Lagrange_2D_Scal)
            Call Init_Element_P_Lagrange_2D_Scal(dElem, dCoord, 2, QuadratureOrder)

!         Case (MEF90_Q1_Lagrange_2D_Scal)
!            Call Init_Element_Q_Lagrange_2D_Scal(dElem, dCoord, 1, QuadratureOrder)
!         Case (MEF90_Q2_Lagrange_2D_Scal)
!            Call Init_Element_Q_Lagrange_2D_Scal(dElem, dCoord, 2, QuadratureOrder)
         Case Default
            Print*, 'Element type not implemented yet', Element_Type
      End Select
   End Subroutine Init_Element_2D_Scal                                
   
   Subroutine Init_Element_P_Lagrange_2D_Scal(dElem, dCoord, PolynomialOrder, QuadratureOrder)
      ! Compute the quadrature weights and the value of the basis functions and their gradient 
      ! at the quadrature points.
      ! One day when I am smart I will use FIAT for that...
      !
      ! Assumes that the elements connectivity is known
      !
      
      Type (Element2D_Scal)                  :: dElem
      PetscReal, Dimension(:,:), Pointer     :: dCoord      ! coord(i,j)=ith coord of jth vertice
      Integer                                :: PolynomialOrder, QuadratureOrder
      
      Integer                                :: Nb_Gauss
      Integer                                :: Nb_Dof
      Integer                                :: iDoF, iG
      Type (Mat2D)                           :: Bt          ! The transposed of transformation matrix
      PetscReal                              :: DetBinv     ! The determinant of B^{-1}

      PetscReal, Dimension(:,:), Pointer     :: PhiHat      ! PhiHat(i,k) The value of the ith basis function at the kth integration point
      Type (Vect2D), Dimension(:,:), Pointer :: GradPhiHat
      
      
      Type (Vect2D), Dimension(:), Pointer   :: Xi ! The quadrature points coordinates in the reference element
      
      !!! The transformation matrix and the determinant of its inverse
      Bt%XX = dCoord(2,3) - dCoord(2,1)
      Bt%XY = dCoord(2,1) - dCoord(2,2)
      Bt%YX = dCoord(1,1) - dCoord(1,3)
      Bt%YY = dCoord(1,2) - dCoord(1,1)
      
      detBinv = Bt%XX * Bt%YY - Bt%XY * Bt%YX
      Bt = Bt / detBinv
      
      Select Case (QuadratureOrder)
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
         Xi(1)%X = 1.0_Kr / 6.0_Kr; Xi(2)%X = 2.0_Kr / 3.0_Kr; Xi(3)%X = 1.0_Kr / 6.0_Kr
         Xi(1)%Y = 1.0_Kr / 6.0_Kr; Xi(2)%Y = 1.0_Kr / 6.0_Kr; Xi(3)%Y = 2.0_Kr / 3.0_Kr
         dElem%Gauss_C = detBinv / 6.0_Kr

      Case(3)
         Nb_Gauss = 4
         Allocate(Xi(Nb_Gauss))
         Allocate(dElem%Gauss_C(Nb_Gauss))
      
      Case Default
         Print*, 'Unimplemented quadrature order', QuadratureOrder
         STOP
      End Select
      
      Select Case (PolynomialOrder)
      Case(1)
         Nb_DoF = 3
         Allocate(PhiHat(Nb_DoF, Nb_Gauss))
         Allocate(GradPhiHat(Nb_DoF, Nb_Gauss))
         PhiHat(1,:) = 1.0_Kr - Xi%X - Xi%Y
         PhiHat(2,:) = Xi(:)%X
         PhiHat(3,:) = Xi(:)%Y
         
         GradPhiHat(1,:)%X = -1.0_Kr; GradPhiHat(1,:)%Y = -1.0_Kr 
         GradPhiHat(2,:)%X =  1.0_Kr; GradPhiHat(2,:)%Y =  0.0_Kr 
         GradPhiHat(3,:)%X =  0.0_Kr; GradPhiHat(3,:)%Y =  1.0_Kr
          
      Case(2)
         Nb_DoF = 6
         Allocate(PhiHat(Nb_DoF, Nb_Gauss))
         Allocate(GradPhiHat(Nb_DoF, Nb_Gauss))
      !!! Needs to be checked
         PhiHat(1,:) = 4.0_Kr * (1.0_Kr - Xi%X - Xi%Y) * Xi%X      
         PhiHat(2,:) = 4.0_Kr * Xi%X * Xi%Y     
         PhiHat(3,:) = 4.0_Kr * (1.0_Kr - Xi%X - Xi%Y) * Xi%Y     
         PhiHat(4,:) = 2.0_Kr * Xi%Y**2 - Xi%Y
         PhiHat(5,:) = 2.0_Kr * (1.0_Kr - Xi%X - Xi%Y)**2  - (1.0_Kr - Xi%X - Xi%Y)
         PhiHat(6,:) = 2.0_Kr * Xi%X**2 - Xi%X
         
         GradPhiHat(1,:)%X = 4.0_Kr * (1.0_Kr - 2.0_Kr * Xi%X - Xi%Y); GradPhiHat(1,:)%Y =-4.0_Kr * Xi%X
         GradPhiHat(2,:)%X = 4.0_Kr * Xi%Y;                            GradPhiHat(2,:)%Y = 4.0_Kr * Xi%X
         GradPhiHat(3,:)%X =-4.0_Kr * Xi%Y;                            GradPhiHat(3,:)%Y = 4.0_Kr * (1.0_Kr - Xi%X - 2.0_Kr * Xi%Y)
         GradPhiHat(4,:)%X = 0.0_Kr;                                   GradPhiHat(4,:)%Y = 4.0_Kr * Xi%Y - 1.0_Kr
         GradPhiHat(5,:)%X = 1.0_Kr - 4.0_Kr * (1.0_Kr - Xi%X - Xi%Y); GradPhiHat(4,:)%Y = 1.0_Kr - 4.0_Kr * (1.0_Kr - Xi%X - Xi%Y)
         GradPhiHat(6,:)%X = 4.0_Kr * Xi%X - 1.0_Kr;                   GradPhiHat(6,:)%Y = 0.0_Kr
      Case Default
         Print*, 'Unimplemented PolynomialOrder', PolynomialOrder
      End Select
      
      Allocate (dElem%BF(Nb_DoF, Nb_Gauss)) 
      Allocate (dElem%Grad_BF(Nb_DoF, Nb_Gauss))
      dElem%BF = PhiHat
      Do iDoF = 1, Nb_DoF
         Do iG = 1, Nb_Gauss
            dElem%Grad_BF(iDoF, iG) = Bt * GradPhiHat(iDoF, iG) 
         End Do
      End Do
     
      DeAllocate(Xi)
      DeAllocate(PhiHat)
      DeAllocate(GradPhiHat)
   End Subroutine Init_Element_P_Lagrange_2D_Scal

End Module m_MEF_InitElem