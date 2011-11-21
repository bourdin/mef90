Program TestBoundaryElement3D

#include "finclude/petscdef.h"


   Use m_MEF90
   Use petsc

   Implicit NONE   

   Character(len=MEF90_MXSTRLEN)                :: exoname
   PetscInt, Dimension(:), Pointer              :: tmpconnect
   Integer, Dimension(:,:), Pointer             :: connect
   PetscInt                                     :: num_elem,num_elem_blk
   PetscInt                                     :: num_vert,num_vert_set
   PetscInt                                     :: num_side,num_side_set
   PetscInt                                     :: num_dim
   PetscReal, Dimension(:), Pointer             :: U1
   Type(Vect3D), Dimension(:), Pointer          :: U2
   PetscReal                                    :: average2,flux
   Real, Dimension(:), Pointer                  :: X,Y,Z
   Real, Dimension(3)                           :: V1,V2,N
   
   PetscReal, Dimension(:,:), Pointer           :: Coord3, Coord2
   PetscInt                                     :: i,j,iE,iG,iDoF
   PetscInt                                     :: num_Gauss,Num_DoF
   Type(BoundaryElement3d)                      :: bElem
   Type(Element2d_Scal)                         :: Elem
   Type(PetscViewer)                            :: viewer
   PetscInt                                     :: IntegOrder = 3
   PetscReal                                    :: area2,area3
   
   Integer                                      :: cpu_ws,io_ws,mod_sz,exoid
   Real                                         :: vers
   Character(len=MXLNLN)                        :: titl
   PetscBool                                    :: flg
   PetscErrorCode                               :: ierr
   
   
   Call MEF90_Initialize()
   viewer = PetscViewer(PETSC_VIEWER_STDOUT_WORLD)
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-i',exoname,flg,ierr);CHKERRQ(ierr);
   cpu_ws = 0
   io_ws = 0
   exoid = exopen (exoname, EXREAD, cpu_ws, io_ws, vers, ierr)
   Call exgini (exoid,titl,num_dim,num_vert,num_elem,num_elem_blk,num_vert_set,num_side_set,ierr)
   Write(*,*) 'num_dim: ',num_dim
   Write(*,*) 'num_vert: ',num_vert
   Write(*,*) 'num_elem: ', num_elem


   Allocate(X(num_vert))
   Allocate(Y(num_vert))
   Allocate(Z(num_vert))

   Allocate(Coord3(3,3)) ! 3 coordinates, 3 vertices
   Allocate(Coord2(2,3)) ! 3 coordinates, 3 vertices
   Allocate(tmpconnect(3*num_elem))
   Allocate(connect(3,num_elem))

   Call exgcor (exoid,X,Y,Z,ierr)
   !Write(*,*) 'X: ',X
   !Write(*,*) 'Y: ',Y
   !Write(*,*) 'Z: ',Z

   Call exgelc (exoid,1,tmpconnect,ierr)
   connect = reshape(tmpconnect,(/3,num_elem/))
   !Write(*,*) 'Connect: '
   !Do i = 1, num_Elem
   !   Write(*,*) i, connect(:,i)
   !End Do
   
   !!! Check that the elements orientation is consistent
   Do iE = 1, num_elem
      V1(1) = X(connect(2,iE)) - X(connect(1,iE))
      V2(1) = X(connect(3,iE)) - X(connect(1,iE))
      V1(2) = Y(connect(2,iE)) - Y(connect(1,iE))
      V2(2) = Y(connect(3,iE)) - Y(connect(1,iE))
      V1(3) = Z(connect(2,iE)) - Z(connect(1,iE))
      V2(3) = Z(connect(3,iE)) - Z(connect(1,iE))
      N(1) = V1(2) * V2(3) - V1(3) * V2(2)
      N(2) = V1(3) * V2(1) - V1(1) * V2(3)
      N(3) = V1(1) * V2(2) - V1(2) * V2(1)
      Write(*,*) 'Normal vector is: ', N
      If (N(3) < 0.) Then
         Write(*,*) 'Flipping orientation of element ',iE
         i = connect(1,iE)
         connect(1,iE) = connect(2,iE)
         connect(2,iE) = i
      End If
   End Do


   Allocate(U1(num_vert))
   Allocate(U2(num_vert))
   Do i = 1, num_vert
      U1(i) = X(i)**2 + Y(i)**2 + Z(i)**2
      U2(i)%X = 0.0_Kr!X(i)
      U2(i)%Y = 0.0_Kr!Y(i)
      U2(i)%Z = 1.0_Kr!Z(i)
      
   End Do
   
   average2 = 0.0_Kr
   flux    = 0.0_Kr
   
   iE = 1
   Coord3(1,:) = X(connect(:,iE))
   Coord3(2,:) = Y(connect(:,iE))
   Coord3(3,:) = Z(connect(:,iE))
   Call BoundaryElement_Init(bElem,Coord3,IntegOrder,MEF90_P1_Lagrange)
   Call BoundaryElement_View(bElem,PetscViewer(PETSC_VIEWER_STDOUT_WORLD))
   Call BoundaryElement_Destroy(bElem)
   
   Do iE = 1, num_Elem
      Coord3(1,:) = X(connect(:,iE))
      Coord3(2,:) = Y(connect(:,iE))
      Coord3(3,:) = Z(connect(:,iE))
      Call BoundaryElement_Init(bElem,Coord3,IntegOrder,MEF90_P1_Lagrange)
      Num_DoF   = size(bElem%BF,1)
      Num_Gauss = size(bElem%BF,2)
      Do iG = 1, Num_Gauss
         area3  = area3 + bElem%Gauss_C(iG)
         Do iDof = 1, Num_DoF
            flux     = flux + bElem%Gauss_C(iG) * (bElem%BF(iDoF,iG) .DotP.  U2(connect(iDoF,iE)))
         End Do
      End Do
      Call BoundaryElement_Destroy(bElem)
      
      Coord2(1,:) = X(connect(:,iE))
      Coord2(2,:) = Y(connect(:,iE))
      Call ElementInit(Elem,Coord2,IntegOrder,MEF90_P1_Lagrange)
      Num_DoF   = size(Elem%BF,1)
      Num_Gauss = size(Elem%BF,2)
      Do iG = 1, Num_Gauss
         area2 = area2 + Elem%Gauss_C(iG)
         Do iDof = 1, Num_DoF
            average2  = average2 + Elem%Gauss_C(iG) * (Elem%BF(iDoF,iG) * U1(connect(iDoF,iE)))
         End Do
      End Do
      Call ElementDestroy(Elem)
   End Do
   Write(*,*) 'Area2 (using 2D body elements) is     ', area2
   Write(*,*) 'Area3 (using 3D boundary elements) is ', area3
   Write(*,*) 'Average2 (using 2D body elements) is     ', average2
   Write(*,*) 'Flux is ', flux   
   DeAllocate(X)
   DeAllocate(Y)
   DeAllocate(Z)
   DeAllocate(Coord2)
   DeAllocate(Coord3)
   DeAllocate(tmpconnect)
   DeAllocate(connect)
   DeAllocate(U1)
   DeAllocate(U2)
   Call exclos(exoid,ierr)   
   Call MEF90_Finalize()
End Program TestBoundaryElement3D
