Module m_MEF90_Elements_Vect
#include "petsc/finclude/petsc.h"
   Use m_MEF90_Elements_Type
   Use m_MEF90_Elements_Scal
   Use m_MEF90_LinAlg
   Use m_MEF90_Utils
   Use m_MEF90_Parameters
   Use petsc
   IMPLICIT NONE

   Public :: MEF90Element2DVect
   Public :: MEF90Element3DVect

   Type MEF90Element2DVect
      Type (Vect2D),Dimension(:,:),Pointer         :: BF
      Type (Mat2D),Dimension(:,:),Pointer          :: Grad_BF
      Type (MatS2D),Dimension(:,:),Pointer         :: GradS_BF
      PetscReal,Dimension(:),Pointer               :: Gauss_C
      Type(Vect2D)                                 :: outerNormal ! only makes sense for elements of codim 1
   End Type MEF90Element2DVect
   
   Type MEF90Element3DVect
      Type (Vect3D),Dimension(:,:),Pointer         :: BF
      Type (Mat3D),Dimension(:,:),Pointer          :: Grad_BF
      Type (MatS3D),Dimension(:,:),Pointer         :: GradS_BF
      PetscReal,Dimension(:),Pointer               :: Gauss_C
      Type(Vect3D)                                 :: outerNormal ! only makes sense for elements of codim 1
   End Type MEF90Element3DVect
 
Contains
#undef __FUNCT__
#define __FUNCT__ "Element2DVectInitSet"
   Subroutine Element2DVectInitSet(dm,cellIS,dElem,dQuadratureOrder,elemType,ierr)
      Type(tDM),intent(IN)                             :: dm
      Type(tIS),intent(IN)                             :: cellIS
      Type(MEF90Element2DVect),Dimension(:),Pointer    :: dElem
      PetscInt,Intent(IN)                              :: dQuadratureOrder
      Type(MEF90ElementType),intent(IN)                :: elemType
      PetscErrorCode,Intent(OUT)                       :: ierr
      
      PetscInt                                         :: iELoc
      PetscInt,Dimension(:),Pointer                    :: CellID
      PetscReal, Dimension(:), Pointer                 :: v0,BB,BBinv
      PetscReal                                        :: detBBinv
      Type(Mat2D)                                      :: Bt
      PetscReal                                        :: length
      PetscReal,dimension(:),pointer                   :: centroid,innerNormal
      Type(Vect2D)                                     :: outerNormal

      Select Case (elemType%name)
         Case (MEF90P1Lagrange2D%name,MEF90P2Lagrange2D%name)
            PetscCall(ISGetIndicesF90(CellIS,CellID,ierr))
            Allocate(dElem(size(cellID)),stat=ierr)
            If (size(CellID) > 0) Then
               Allocate(v0(2))
               Allocate(BB(4))
               Allocate(BBinv(4))
               Do_Elem_iE: Do iELoc = 1,size(CellID)
                  PetscCall(DMPlexComputeCellGeometryAffineFEM(dm,cellID(iELoc),v0,BB,BBinv,detBBinv,ierr))
                  !!! Petsc uses a reference simplex with vertices at (-1,-1), (1,-1) and (-1,1)
                  !!! Whereas MEF90 uses (0,0), (1,0), (0,1), so we need to rescale the affine transform
                  Bt%XX = BBinv(1)*0.5_Kr; Bt%XY = BBinv(3)*0.5_Kr
                  Bt%YX = BBinv(2)*0.5_Kr; Bt%YY = BBinv(4)*0.5_Kr
                  detBBinv = 4.0_Kr * detBBinv
                  Call ElementPLagrange2DVectInit(dElem(iELoc),Bt,detBBinv,elemType%order,dQuadratureOrder,ierr)
               End Do Do_Elem_iE
               DeAllocate(BBinv)
               DeAllocate(BB)
               DeAllocate(v0)
            End If
            PetscCall(ISRestoreIndicesF90(CellIS,CellID,ierr))
         Case (MEF90P1Lagrange2DBoundary%name,MEF90P2Lagrange2DBoundary%name)
            PetscCall(ISGetIndicesF90(CellIS,CellID,ierr))
            Allocate(dElem(size(cellID)),stat=ierr)
            If (size(CellID) > 0) Then
               allocate(centroid(2))
               allocate(innerNormal(2))
               Do iELoc = 1,size(CellID)
                  PetscCall(DMPlexComputeCellGeometryFVM(dm,cellID(iEloc),length,centroid,innerNormal,ierr))
                  outerNormal = -innerNormal
                  Call ElementPLagrange2DBoundaryVectInit(dElem(iELoc),length,outerNormal,elemType%order,dQuadratureOrder,ierr)
               End Do
               DeAllocate(innerNormal)
               DeAllocate(centroid)
            End If
            PetscCall(ISRestoreIndicesF90(CellIS,CellID,ierr))
         !Case (MEF90Q1Lagrange2D%name,MEF90Q2Lagrange2D%name,MEF90Q1Lagrange2DBoundary%name,MEF90Q2Lagrange2DBoundary%name)
         !   !!! Get quadrature points for the current element using DMPlexComputeCellGeometryFEM
         !   !!! Initialize element
         Case Default
            Write(*,*) __FUNCT__,': Element type not implemented yet',elemType%name
            ierr = PETSC_ERR_SUP
      End Select
   End Subroutine Element2DVectInitSet

#undef __FUNCT__
#define __FUNCT__ "Element3DVectInitSet"
   Subroutine Element3DVectInitSet(dm,cellIS,dElem,dQuadratureOrder,elemType,ierr)
      Type(tDM),intent(IN)                             :: dm
      Type(tIS),intent(IN)                             :: cellIS
      Type(MEF90Element3DVect),Dimension(:),Pointer    :: dElem
      PetscInt,Intent(IN)                              :: dQuadratureOrder
      Type(MEF90ElementType),intent(IN)                :: elemType
      PetscErrorCode,Intent(OUT)                       :: ierr
      
      PetscInt                                         :: iELoc
      PetscInt,Dimension(:),Pointer                    :: CellID
      PetscReal, Dimension(:), Pointer                 :: v0,BB,BBinv
      PetscReal                                        :: detBBinv
      Type(Mat3D)                                      :: Bt
      PetscReal                                        :: area
      PetscReal,dimension(:),pointer                   :: centroid,innerNormal
      Type(Vect3D)                                     :: outerNormal

      Select Case (elemType%name)
         Case (MEF90P1Lagrange3D%name,MEF90P2Lagrange3D%name)
            PetscCall(ISGetIndicesF90(CellIS,CellID,ierr))
            Allocate(dElem(size(cellID)),stat=ierr)
            If (size(CellID) > 0) Then
               Allocate(v0(3))
               Allocate(BB(9))
               Allocate(BBinv(9))
               Do_Elem_iE: Do iELoc = 1,size(CellID)
                  PetscCall(DMPlexComputeCellGeometryAffineFEM(dm,cellID(iELoc),v0,BB,BBinv,detBBinv,ierr))
                  !!! Petsc uses a reference simplex with vertices at (-1,-1), (1,-1) and (-1,1)
                  !!! Whereas MEF90 uses (0,0), (1,0), (0,1), so we need to rescale the affine transform
                  Bt%XX = BBinv(1)*0.5_Kr; Bt%XY = BBinv(4)*0.5_Kr; Bt%XZ = BBinv(7)*0.5_Kr
                  Bt%YX = BBinv(2)*0.5_Kr; Bt%YY = BBinv(5)*0.5_Kr; Bt%YZ = BBinv(8)*0.5_Kr
                  Bt%ZX = BBinv(3)*0.5_Kr; Bt%ZY = BBinv(6)*0.5_Kr; Bt%ZZ = BBinv(9)*0.5_Kr
                  detBBinv = 8.0_Kr * detBBinv
                  Call ElementPLagrange3DVectInit(dElem(iELoc),Bt,detBBinv,elemType%order,dQuadratureOrder,ierr)
               End Do Do_Elem_iE
               DeAllocate(BBinv)
               DeAllocate(BB)
               DeAllocate(v0)
            End If
            PetscCall(ISRestoreIndicesF90(CellIS,CellID,ierr))
         Case (MEF90P1Lagrange3DBoundary%name,MEF90P2Lagrange3DBoundary%name)
            PetscCall(ISGetIndicesF90(CellIS,CellID,ierr))
            Allocate(dElem(size(cellID)),stat=ierr)
            If (size(CellID) > 0) Then
               allocate(centroid(3))
               allocate(innerNormal(3))
               Do iELoc = 1,size(CellID)
                  PetscCall(DMPlexComputeCellGeometryFVM(dm,cellID(iEloc),area,centroid,innerNormal,ierr))
                  outerNormal = -innerNormal
                  Call ElementPLagrange3DBoundaryVectInit(dElem(iELoc),area,outerNormal,elemType%order,dQuadratureOrder,ierr)
               End Do
               DeAllocate(innerNormal)
               DeAllocate(centroid)
            End If
            PetscCall(ISRestoreIndicesF90(CellIS,CellID,ierr))
         !Case (MEF90Q1Lagrange3D%name,MEF90Q2Lagrange3D%name,MEF90Q1Lagrange3DBoundary%name,MEF90Q2Lagrange3DBoundary%name)
         !   !!! Get quadrature points for the current element using DMPlexComputeCellGeometryFEM
         !   !!! Initialize element
         Case Default
            Write(*,*) __FUNCT__,': Element type not implemented yet',elemType%name
            ierr = PETSC_ERR_SUP
      End Select
   End Subroutine Element3DVectInitSet

      
#undef __FUNCT__
#define __FUNCT__ "ElementPLagrange2DVectInit"
   Subroutine ElementPLagrange2DVectInit(dElem,Bt,DetBinv,dPolynomialOrder,dQuadratureOrder,ierr)
      Type(MEF90Element2DVect),Intent(INOUT) :: dElem
      Type(Mat2D),Intent(IN)                 :: Bt          ! The transposed of transformation matrix
      PetscReal,Intent(IN)                   :: DetBinv     ! The determinant of B^{-1}
      PetscInt,Intent(IN)                    :: dPolynomialOrder,dQuadratureOrder
      PetscErrorCode,Intent(OUT)             :: ierr
   
      Type(MEF90Element2DScal)              :: Elem_Scal
      PetscInt                               :: dim = 2 
      PetscInt                               :: Num_DoF,Nb_Gauss,i,iDof,iG
      
      
      Call ElementPLagrange2DScalInit(Elem_Scal,Bt,DetBinv,dPolynomialOrder,dQuadratureOrder,ierr)
      Num_DoF  = Size(Elem_Scal%BF,1) 
      Nb_Gauss = Size(Elem_Scal%BF,2)

      Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
      Allocate(dElem%BF(Num_DoF * dim,Nb_Gauss),stat=ierr)
      Allocate(dElem%Grad_BF(Num_DoF * dim,Nb_Gauss),stat=ierr)
      Allocate(dElem%GradS_BF(Num_DoF * dim,Nb_Gauss),stat=ierr)
         
      dElem%Gauss_C = Elem_Scal%Gauss_C
      Do iDof = 1, num_dof * dim
         Do iG = 1, Nb_Gauss
            dElem%BF(iDof,iG)       = 0.0_Kr
            delem%Grad_BF(iDof,iG)  = 0.0_Kr
            delem%GradS_BF(iDof,iG) = 0.0_Kr
         End Do
      End Do
      
      Do i = 0,Num_DoF-1
         dElem%BF(i*dim+1,:)%X = Elem_Scal%BF(i+1,:)
         dElem%BF(i*dim+2,:)%Y = Elem_Scal%BF(i+1,:)

         dElem%Grad_BF(i*dim+1,:)%XX = Elem_Scal%Grad_BF(i+1,:)%X
         dElem%Grad_BF(i*dim+1,:)%XY = Elem_Scal%Grad_BF(i+1,:)%Y
         dElem%Grad_BF(i*dim+2,:)%YX = Elem_Scal%Grad_BF(i+1,:)%X
         dElem%Grad_BF(i*dim+2,:)%YY = Elem_Scal%Grad_BF(i+1,:)%Y

         dElem%GradS_BF(i*dim+1,:)%XX = Elem_Scal%Grad_BF(i+1,:)%X
         dElem%GradS_BF(i*dim+1,:)%XY = Elem_Scal%Grad_BF(i+1,:)%Y / 2.0_Kr
         dElem%GradS_BF(i*dim+2,:)%XY = Elem_Scal%Grad_BF(i+1,:)%X / 2.0_Kr
         dElem%GradS_BF(i*dim+2,:)%YY = Elem_Scal%Grad_BF(i+1,:)%Y
      End Do
      
      Call MEF90ElementDestroy_2DScal(Elem_Scal,ierr)
   End Subroutine ElementPLagrange2DVectInit

#undef __FUNCT__
#define __FUNCT__ "ElementPLagrange2DBoundaryVectInit"
   Subroutine ElementPLagrange2DBoundaryVectInit(dElem,length,outerNormal,dPolynomialOrder,dQuadratureOrder,ierr)
      Type(MEF90Element2DVect),Intent(INOUT) :: dElem
      PetscReal,Intent(IN)                   :: length
      Type(Vect2D),Intent(IN)                :: outerNormal 
      PetscInt,Intent(IN)                    :: dPolynomialOrder,dQuadratureOrder
      PetscErrorCode,Intent(OUT)             :: ierr
   
      Type(MEF90Element2DScal)              :: Elem_Scal
      PetscInt                               :: dim = 2 
      PetscInt                               :: Num_DoF,Nb_Gauss,iDof,iG
      
      
      Call ElementPLagrange2DBoundaryScalInit(Elem_Scal,length,outerNormal,dPolynomialOrder,dQuadratureOrder,ierr)
      Num_DoF  = Size(Elem_Scal%BF,1) 
      Nb_Gauss = Size(Elem_Scal%BF,2)
      Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
      Allocate(dElem%BF(Num_DoF * dim,Nb_Gauss),stat=ierr)
      Allocate(dElem%Grad_BF(Num_DoF * dim,Nb_Gauss),stat=ierr)
      Allocate(dElem%GradS_BF(Num_DoF * dim,Nb_Gauss),stat=ierr)
         
      dElem%Gauss_C = Elem_Scal%Gauss_C
      Do iDof = 1, num_dof * dim
         Do iG = 1, Nb_Gauss
            dElem%BF(iDof,iG)       = 0.0_Kr
            delem%Grad_BF(iDof,iG)  = 0.0_Kr
            delem%GradS_BF(iDof,iG) = 0.0_Kr
         End Do
      End Do

      Do idof = 0,Num_DoF-1
         dElem%BF(iDof*dim+1,:)%X = Elem_Scal%BF(iDof+1,:)
         dElem%BF(iDof*dim+2,:)%Y = Elem_Scal%BF(iDof+1,:)
      End Do

      dElem%outerNormal = Elem_Scal%outerNormal
      Call MEF90ElementDestroy_2DScal(Elem_Scal,ierr)
   End Subroutine ElementPLagrange2DBoundaryVectInit
   
#undef __FUNCT__
#define __FUNCT__ "ElementPLagrange3DVectInit"
   Subroutine ElementPLagrange3DVectInit(dElem,Bt,detBinv,dPolynomialOrder,dQuadratureOrder,ierr)
      Type(MEF90Element3DVect),Intent(INOUT) :: dElem
      Type(Mat3D),Intent(IN)                 :: Bt
      PetscReal,Intent(IN)                   :: detBinv
      PetscInt,Intent(IN)                    :: dPolynomialOrder,dQuadratureOrder
      PetscErrorCode,Intent(OUT)             :: ierr
   
      Type(MEF90Element3DScal)              :: Elem_Scal
      PetscInt                               :: dim = 3 
      PetscInt                               :: Num_DoF,Nb_Gauss,i,iDof,iG
      
      
      Call ElementPLagrange3DScalInit(Elem_Scal,Bt,detBinv,dPolynomialOrder,dQuadratureOrder,ierr)
      Num_DoF   = Size(Elem_Scal%BF,1) 
      Nb_Gauss = Size(Elem_Scal%BF,2)
      Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
      Allocate(dElem%BF(Num_DoF * dim,Nb_Gauss),stat=ierr)
      Allocate(dElem%Grad_BF(Num_DoF * dim,Nb_Gauss),stat=ierr)
      Allocate(dElem%GradS_BF(Num_DoF * dim,Nb_Gauss),stat=ierr)
         
      dElem%Gauss_C = Elem_Scal%Gauss_C
      Do iDof = 1, num_dof * dim
         Do iG = 1, Nb_Gauss
            dElem%BF(iDof,iG)       = 0.0_Kr
            delem%Grad_BF(iDof,iG)  = 0.0_Kr
            delem%GradS_BF(iDof,iG) = 0.0_Kr
         End Do
      End Do
      
      Do i = 0,Num_DoF-1
         dElem%BF(i*dim+1,:)%X = Elem_Scal%BF(i+1,:)
         dElem%BF(i*dim+2,:)%Y = Elem_Scal%BF(i+1,:)
         dElem%BF(i*dim+3,:)%Z = Elem_Scal%BF(i+1,:)
         dElem%Grad_BF(i*dim+1,:)%XX = Elem_Scal%Grad_BF(i+1,:)%X
         dElem%Grad_BF(i*dim+1,:)%XY = Elem_Scal%Grad_BF(i+1,:)%Y
         dElem%Grad_BF(i*dim+1,:)%XZ = Elem_Scal%Grad_BF(i+1,:)%Z
         dElem%Grad_BF(i*dim+2,:)%YX = Elem_Scal%Grad_BF(i+1,:)%X
         dElem%Grad_BF(i*dim+2,:)%YY = Elem_Scal%Grad_BF(i+1,:)%Y
         dElem%Grad_BF(i*dim+2,:)%YZ = Elem_Scal%Grad_BF(i+1,:)%Z
         dElem%Grad_BF(i*dim+3,:)%ZX = Elem_Scal%Grad_BF(i+1,:)%X
         dElem%Grad_BF(i*dim+3,:)%ZY = Elem_Scal%Grad_BF(i+1,:)%Y
         dElem%Grad_BF(i*dim+3,:)%ZZ = Elem_Scal%Grad_BF(i+1,:)%Z

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
      Call MEF90ElementDestroy_3DScal(Elem_Scal,ierr)
   End Subroutine ElementPLagrange3DVectInit

#undef __FUNCT__
#define __FUNCT__ "ElementPLagrange3DBoundaryVectInit"
   Subroutine ElementPLagrange3DBoundaryVectInit(dElem,area,outerNormal,dPolynomialOrder,dQuadratureOrder,ierr)
      Type(MEF90Element3DVect),Intent(INOUT) :: dElem
      PetscReal,Intent(IN)                   :: area
      Type(Vect3D),Intent(IN)                :: outerNormal
      PetscInt                               :: dPolynomialOrder,dQuadratureOrder
      PetscErrorCode,Intent(OUT)             :: ierr
   
      Type(MEF90Element3DScal)              :: Elem_Scal
      PetscInt                               :: dim = 3 
      PetscInt                               :: Num_DoF,Nb_Gauss,iDof,iG
      
      
      Call ElementPLagrange3DBoundaryScalInit(Elem_Scal,area,outerNormal,dPolynomialOrder,dQuadratureOrder,ierr)
      Num_DoF   = Size(Elem_Scal%BF,1) 
      Nb_Gauss = Size(Elem_Scal%BF,2)
      Allocate(dElem%Gauss_C(Nb_Gauss),stat=ierr)
      Allocate(dElem%BF(Num_DoF * dim,Nb_Gauss),stat=ierr)
      Allocate(dElem%Grad_BF(Num_DoF * dim,Nb_Gauss),stat=ierr)
      Allocate(dElem%GradS_BF(Num_DoF * dim,Nb_Gauss),stat=ierr)
         
      dElem%Gauss_C = Elem_Scal%Gauss_C
      Do iDof = 1, num_dof * dim
         Do iG = 1, Nb_Gauss
            dElem%BF(iDof,iG)       = 0.0_Kr
            delem%Grad_BF(iDof,iG)  = 0.0_Kr
            delem%GradS_BF(iDof,iG) = 0.0_Kr
         End Do
      End Do
      Do iDof = 0,Num_DoF-1
         dElem%BF(iDof*dim+1,:)%X = Elem_Scal%BF(iDof+1,:)
         dElem%BF(iDof*dim+2,:)%Y = Elem_Scal%BF(iDof+1,:)
         dElem%BF(iDof*dim+3,:)%Z = Elem_Scal%BF(iDof+1,:)
      End Do

      dElem%outerNormal = Elem_Scal%outerNormal
      Call MEF90ElementDestroy_3DScal(Elem_Scal,ierr)
   End Subroutine ElementPLagrange3DBoundaryVectInit


#undef __FUNCT__
#define __FUNCT__ "Element2DVectDestroy"
   Subroutine Element2DVectDestroy(dElem,ierr)
      Type(MEF90Element2DVect),Intent(INOUT) :: dElem
      PetscErrorCode,Intent(OUT)             :: ierr
      
      If (Associated(dElem%BF)) Then
         DeAllocate(dElem%BF,stat=ierr)
      End If
      If (Associated(dElem%Grad_BF)) Then
         DeAllocate(dElem%Grad_BF,stat=ierr)
      End If
      If (Associated(dElem%GradS_BF)) Then
         DeAllocate(dElem%GradS_BF,stat=ierr)
      End If
      If (Associated(dElem%Gauss_C)) Then
         DeAllocate(dElem%Gauss_C,stat=ierr)
      End If
   End Subroutine Element2DVectDestroy
      
#undef __FUNCT__
#define __FUNCT__ "Element3DVectDestroy"
   Subroutine Element3DVectDestroy(dElem,ierr)
      Type(MEF90Element3DVect),Intent(INOUT) :: dElem
      PetscErrorCode,Intent(OUT)             :: ierr
      If (Associated(dElem%BF)) Then
         DeAllocate(dElem%BF,stat=ierr)
      End If
      If (Associated(dElem%Grad_BF)) Then
         DeAllocate(dElem%Grad_BF,stat=ierr)
      End If
      If (Associated(dElem%GradS_BF)) Then
         DeAllocate(dElem%GradS_BF,stat=ierr)
      End If
      If (Associated(dElem%Gauss_C)) Then
         DeAllocate(dElem%Gauss_C,stat=ierr)
      End If
   End Subroutine Element3DVectDestroy
   
#undef __FUNCT__
#define __FUNCT__ "Element2DVectDestroySet"
   Subroutine Element2DVectDestroySet(dElem,ierr)
      Type(MEF90Element2DVect),dimension(:),Pointer    :: dElem
      PetscErrorCode,Intent(OUT)                        :: ierr
      
      PetscInt                                          :: cell
      
      Do cell = 1, size(dElem)
         Call Element2DVectDestroy(dElem(cell),ierr)
      End Do
      DeAllocate(dElem,stat=ierr)
   End Subroutine Element2DVectDestroySet

#undef __FUNCT__
#define __FUNCT__ "Element3DVectDestroySet"
   Subroutine Element3DVectDestroySet(dElem,ierr)
      Type(MEF90Element3DVect),dimension(:),Pointer    :: dElem
      PetscErrorCode,Intent(OUT)                        :: ierr
      
      PetscInt                                          :: cell
        
      Do cell = 1, size(dElem)
         Call Element3DVectDestroy(dElem(cell),ierr)
      End Do
      DeAllocate(dElem,stat=ierr)
   End Subroutine Element3DVectDestroySet

End Module m_MEF90_Elements_Vect
