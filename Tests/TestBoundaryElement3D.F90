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
   PetscReal, Dimension(:), Pointer             :: UScal
   Type(Vect3D), Dimension(:), Pointer          :: UVect
   PetscReal                                    :: average2,average3,flux
   Real, Dimension(:), Pointer                  :: X,Y,Z
   PetscReal, Dimension(:,:), Pointer           :: Coord3, Coord2
   PetscInt                                     :: i,iE,iG,iDoF
   PetscInt                                     :: num_Gauss,Num_DoF
   Type(BoundaryElement3d_Scal)                 :: bElem
   Type(Element2d_Scal)                         :: Elem
   Type(PetscViewer)                            :: viewer
   PetscInt                                     :: IntegOrder = 2
   PetscReal                                    :: area2,area3
   PetscReal                                    :: inc2,inc3
   
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
   Do i = 1, num_Elem
      Write(*,*) i, connect(:,i)
   End Do
   
   Allocate(UScal(num_vert))
   Allocate(UVect(num_vert))
   Do i = 1, num_vert
      UScal(i) = 1.0_Kr
      UVect(i)%X = 0.0_Kr
      UVect(i)%Y = 0.0_Kr
      UVect(i)%Z = 1.0_Kr
   End Do
   
   average2 = 0.0_Kr
   average3 = 0.0_Kr
   flux    = 0.0_Kr
   Do iE = 1, num_Elem
      area3 = 0.0_Kr
      inc3 = 0.0_Kr
      Coord3(1,:) = X(connect(:,iE))
      Coord3(2,:) = Y(connect(:,iE))
      Coord3(3,:) = Z(connect(:,iE))
      Call BoundaryElement_Init(bElem,Coord3,IntegOrder,MEF90_P1_Lagrange)
      Num_DoF   = size(bElem%BF,1)
      Num_Gauss = size(bElem%BF,2)
      Do iG = 1, Num_Gauss
         area3  = area3 + bElem%Gauss_C(iG)
         Do iDof = 1, Num_DoF
            inc3 = inc3 + bElem%Gauss_C(iG) * (bElem%BF(iDoF,iG) * UScal(connect(iE,iDoF)))
         End Do
      End Do
      average3 = average3 + inc3
      !Write(*,*) 'Area (boundary)', sum(bElem%Gauss_C)
      
      area2 = 0.0_Kr
      inc2 = 0.0_Kr
      Coord2(1,:) = X(connect(:,iE))
      Coord2(2,:) = Y(connect(:,iE))
      Call ElementInit(Elem,Coord2,IntegOrder,MEF90_P1_Lagrange)
      Num_DoF   = size(Elem%BF,1)
      Num_Gauss = size(Elem%BF,2)
      Do iG = 1, Num_Gauss
         area2 = area2 + Elem%Gauss_C(iG)
         Do iDof = 1, Num_DoF
            inc2  = inc2 + Elem%Gauss_C(iG) * (Elem%BF(iDoF,iG) * UScal(connect(iE,iDoF)))
         End Do
      End Do
      average2 = average2 + inc2
      !Write(*,*) 'Area (Body)', sum(Elem%Gauss_C)
      !Write(*,*) 'inc ratio ', inc2 / inc3
      !Write(*,*) 'area ratio ', area2 / area3
      Call BoundaryElement_Destroy(bElem)
      Call ElementDestroy(Elem)
   End Do
   Write(*,*) 'Average2 (using 2D body elements) is     ', average2
   Write(*,*) 'Average3 (using 3D boundary elements) is ', average3
   
   DeAllocate(X)
   DeAllocate(Y)
   DeAllocate(Z)
   DeAllocate(Coord2)
   DeAllocate(Coord3)
   DeAllocate(tmpconnect)
   DeAllocate(connect)
   DeAllocate(UScal)
   DeAllocate(UVect)
   Call exclos(exoid,ierr)   
   Call MEF90_Finalize()
End Program TestBoundaryElement3D
