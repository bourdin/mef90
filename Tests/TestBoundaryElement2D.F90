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
   Type(Vect2D), Dimension(:), Pointer          :: U2
   PetscReal                                    :: area,flux
   Real, Dimension(:), Pointer                  :: X,Y,Z
   Type(Vect2D)                                 :: N
   
   PetscReal, Dimension(:,:), Pointer           :: Coord2
   PetscInt                                     :: i,j,iE,iG,iDoF
   PetscInt                                     :: num_Gauss,Num_DoF
   Type(BoundaryElement2d)                      :: bElem
   Type(Element2d_Scal)                         :: Elem
   Type(PetscViewer)                            :: viewer
   PetscInt                                     :: IntegOrder = 3

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

   Allocate(Coord2(2,2)) ! 3 coordinates, 3 vertices
   Allocate(tmpconnect(2*num_elem))
   Allocate(connect(2,num_elem))

   Call exgcor (exoid,X,Y,Z,ierr)
   !Write(*,*) 'X: ',X
   !Write(*,*) 'Y: ',Y
   !Write(*,*) 'Z: ',Z

   Call exgelc (exoid,1,tmpconnect,ierr)
   connect = reshape(tmpconnect,(/2,num_elem/))
   !Write(*,*) 'Connect: '
   !Do i = 1, num_Elem
   !   Write(*,*) i, connect(:,i)
   !End Do
   
   !!! Check that the elements orientation is consistent
   Do iE = 1, num_elem
      N%X = Y(connect(2,iE)) - Y(connect(1,iE))
      N%Y = X(connect(1,iE)) - X(connect(2,iE))
      N = N / Norm(N)
      If (N%X * X(connect(1,iE)) + N%Y * Y(connect(1,iE)) < 0.) Then
         Write(*,*) 'Flipping orientation of element ',iE
         i = connect(1,iE)
         connect(1,iE) = connect(2,iE)
         connect(2,iE) = i
      End If
   End Do


   Allocate(U2(num_vert))
   Do i = 1, num_vert
      U2(i)%X = X(i)
      U2(i)%Y = Y(i)
   End Do
   
   area = 0.0_Kr
   flux = 0.0_Kr
   
   iE = 1
   Coord2(1,:) = X(connect(:,iE))
   Coord2(2,:) = Y(connect(:,iE))
   Call BoundaryElement_Init(bElem,Coord2,IntegOrder,MEF90_P1_Lagrange)
   Call BoundaryElement_View(bElem,PetscViewer(PETSC_VIEWER_STDOUT_WORLD))
   Call BoundaryElement_Destroy(bElem)
   
   Do iE = 1, num_Elem
      Coord2(1,:) = X(connect(:,iE))
      Coord2(2,:) = Y(connect(:,iE))
      Call BoundaryElement_Init(bElem,Coord2,IntegOrder,MEF90_P1_Lagrange)
      Num_DoF   = size(bElem%BF,1)
      Num_Gauss = size(bElem%BF,2)
      Do iG = 1, Num_Gauss
         area  = area + bElem%Gauss_C(iG)
         Do iDof = 1, Num_DoF
            flux     = flux + bElem%Gauss_C(iG) * (bElem%BF(iDoF,iG) .DotP.  U2(connect(iDoF,iE)))
         End Do
      End Do
      Call BoundaryElement_Destroy(bElem)
      
   End Do
   Write(*,*) 'Area is ', area
   Write(*,*) 'Flux is ', flux   
   DeAllocate(X)
   DeAllocate(Y)
   DeAllocate(Z)
   DeAllocate(Coord2)
   DeAllocate(tmpconnect)
   DeAllocate(connect)
   DeAllocate(U2)
   Call exclos(exoid,ierr)   
   Call MEF90_Finalize()
End Program TestBoundaryElement3D
