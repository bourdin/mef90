Module m_MEF_InitElem

   Use m_AlgebLin
   Use m_MEF_Types
   Use m_Utils
   
   IMPLICIT NONE
   Private
   
#include "include/finclude/petsc.h"

   Public :: Init_Element
   Public :: Destroy_Element
   Public :: Init_Elem_Blk_Info


   Interface Init_Element
      Module Procedure Init_Element2D_Scal, Init_Element2D, Init_Element2D_Elast
   End Interface Init_Element
   
   Interface Destroy_Element
      Module Procedure Destroy_Element2D_Scal, Destroy_Element2D, Destroy_Element2D_Elast, Destroy_Element3D_Scal, Destroy_Element3D, Destroy_Element3D_Elast
   End Interface Destroy_Element
   
   Integer, Parameter, Public                    :: MEF90_P1_Lagrange = 1
   Integer, Parameter, Public                    :: MEF90_P2_Lagrange = 2
   
!   Integer, Parameter, Public                    :: MEF90_Q1_Lagrange = 3
!   Integer, Parameter, Public                    :: MEF90_Q2_Lagrange = 4

 Contains
   Subroutine Init_Elem_Blk_Info(dBlk, dDim)
      Type (Elem_Blk_Info)                   :: dBlk
      Integer                                :: dDim
      
      Select Case (dDim)
      Case (2)
         Select Case (dBlk%Elem_Type)
         Case (MEF90_P1_Lagrange)         
               dBlk%DoF_Location = (/ 0, 0, 0, 3/)
         Case (MEF90_P2_Lagrange)
            dBlk%DoF_Location = (/ 0, 0, 3, 3/)
         Case Default
            Print*, 'Unknown element type', dBlk%Elem_Type
            STOP
         End Select
      Case (3)
         Select Case (dBlk%Elem_Type)
         Case (MEF90_P1_Lagrange)         
               dBlk%DoF_Location = (/ 0, 0, 0, 4/)
         Case (MEF90_P2_Lagrange)
            dBlk%DoF_Location = (/ 0, 0, 6, 4/)
         Case Default
            Print*, 'Unknown element type', dBlk%Elem_Type
            STOP
         End Select
      Case Default
         Print*, 'Unknown dimension', dDim
         STOP
      End Select
      dBlk%Nb_DoF = sum(dBlk%DoF_Location)

   End Subroutine Init_Elem_Blk_Info      
      
   
   Subroutine Init_Element2D_Scal(dElem, dCoord, QuadratureOrder, Element_Type)
      Type (Element2D_Scal)                  :: dElem
      PetscReal, Dimension(:,:), Pointer     :: dCoord
      Integer, Intent(IN)                    :: QuadratureOrder
      Integer, Intent(IN)                    :: Element_Type
      
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
   
   Subroutine Init_Element2D(dElem, dCoord, QuadratureOrder, Element_Type)
      Type (Element2D)                       :: dElem
      PetscReal, Dimension(:,:), Pointer     :: dCoord
      Integer, Intent(IN)                    :: QuadratureOrder
      Integer, Intent(IN)                    :: Element_Type
      
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

   Subroutine Init_Element2D_Elast(dElem, dCoord, QuadratureOrder, Element_Type)
      Type (Element2D_Elast)                 :: dElem
      PetscReal, Dimension(:,:), Pointer     :: dCoord
      Integer, Intent(IN)                    :: QuadratureOrder
      Integer, Intent(IN)                    :: Element_Type
      
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


   Subroutine Init_Element_P_Lagrange_2D_Scal(dElem, dCoord, dPolynomialOrder, dQuadratureOrder)
      ! Compute the quadrature weights and the value of the basis functions and their gradient 
      ! at the quadrature points.
      ! One day when I am smart I will use FIAT for that...
      !
      ! Assumes that the elements connectivity is known
      !
      
      Type (Element2D_Scal)                  :: dElem
      PetscReal, Dimension(:,:), Pointer     :: dCoord      ! coord(i,j)=ith coord of jth vertice
      Integer                                :: dPolynomialOrder, dQuadratureOrder
      
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
         Print*, 'Unimplemented PolynomialOrder', dPolynomialOrder
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
   
   Subroutine Init_Element_P_Lagrange_2D(dElem, dCoord, dPolynomialOrder, dQuadratureOrder)
      Type (Element2D)                       :: dElem
      PetscReal, Dimension(:,:), Pointer     :: dCoord      ! coord(i,j)=ith coord of jth vertice
      Integer                                :: dPolynomialOrder, dQuadratureOrder
   
      Type (Element2D_Scal)                  :: Elem_Scal
      Integer                                :: dim = 2 
      Integer                                :: Nb_DoF, Nb_Gauss, i
      
      
      Call Init_Element_P_Lagrange_2D_Scal(Elem_Scal, dCoord, dPolynomialOrder, dQuadratureOrder)
      Nb_DoF   = Size(Elem_Scal%BF, 1) 
      Nb_Gauss = Size(Elem_Scal%BF, 2)
      Allocate(dElem%Gauss_C(Nb_Gauss))
      Allocate(dElem%BF(Nb_DoF * dim, Nb_Gauss))
      Allocate(dElem%Der_BF(Nb_DoF * dim, Nb_Gauss))
         
      dElem%Gauss_C = Elem_Scal%Gauss_C
      dElem%BF(:,:)%X = 0.0_Kr
      dElem%BF(:,:)%Y = 0.0_Kr
      dElem%Der_BF(:,:)%XX = 0.0_Kr
      dElem%Der_BF(:,:)%XY = 0.0_Kr
      dElem%Der_BF(:,:)%YX = 0.0_Kr
      dElem%Der_BF(:,:)%YY = 0.0_Kr
      
      Do i = 0, Nb_DoF-1
         dElem%BF(i*dim+1,:)%X = Elem_Scal%BF(i+1,:)
         dElem%BF(i*dim+2,:)%Y = Elem_Scal%BF(i+1,:)
         dElem%Der_BF(i*dim+1,:)%XX = Elem_Scal%Grad_BF(i+1,:)%X
         dElem%Der_BF(i*dim+1,:)%XY = Elem_Scal%Grad_BF(i+1,:)%Y
         dElem%Der_BF(i*dim+2,:)%YX = Elem_Scal%Grad_BF(i+1,:)%X
         dElem%Der_BF(i*dim+2,:)%YY = Elem_Scal%Grad_BF(i+1,:)%Y
      End Do
      Call Destroy_Element(Elem_Scal)
   End Subroutine Init_Element_P_Lagrange_2D

   Subroutine Init_Element_P_Lagrange_2D_Elast(dElem, dCoord, dPolynomialOrder, dQuadratureOrder)
      Type (Element2D_Elast)                 :: dElem
      PetscReal, Dimension(:,:), Pointer     :: dCoord      ! coord(i,j)=ith coord of jth vertice
      Integer                                :: dPolynomialOrder, dQuadratureOrder
   
      Type (Element2D_Scal)                  :: Elem_Scal
      Integer                                :: dim = 2 
      Integer                                :: Nb_DoF, Nb_Gauss, i
      
      
      Call Init_Element_P_Lagrange_2D_Scal(Elem_Scal, dCoord, dPolynomialOrder, dQuadratureOrder)
      Nb_DoF   = Size(Elem_Scal%BF, 1) 
      Nb_Gauss = Size(Elem_Scal%BF, 2)
      Allocate(dElem%Gauss_C(Nb_Gauss))
      Allocate(dElem%BF(Nb_DoF * dim, Nb_Gauss))
      Allocate(dElem%GradS_BF(Nb_DoF * dim, Nb_Gauss))
         
      dElem%Gauss_C = Elem_Scal%Gauss_C
      dElem%BF(:,:)%X = 0.0_Kr
      dElem%BF(:,:)%Y = 0.0_Kr
      dElem%GradS_BF(:,:)%XX = 0.0_Kr
      dElem%GradS_BF(:,:)%XY = 0.0_Kr
      dElem%GradS_BF(:,:)%YY = 0.0_Kr
      Do i = 0, Nb_DoF-1
         dElem%BF(i*dim+1,:)%X = Elem_Scal%BF(i+1,:)
         dElem%BF(i*dim+2,:)%Y = Elem_Scal%BF(i+1,:)
         dElem%GradS_BF(i*dim+1,:)%XX = Elem_Scal%Grad_BF(i+1,:)%X
         dElem%GradS_BF(i*dim+1,:)%XY = (Elem_Scal%Grad_BF(i+1,:)%Y + Elem_Scal%Grad_BF(i+1,:)%X) / 2.0_Kr
         dElem%GradS_BF(i*dim+2,:)%YY = Elem_Scal%Grad_BF(i+1,:)%Y
      End Do
      Call Destroy_Element(Elem_Scal)
   End Subroutine Init_Element_P_Lagrange_2D_Elast

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
End Module m_MEF_InitElem