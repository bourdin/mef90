Module m_MEF_MassMatrix
#include "finclude/petscdef.h"
   Use m_MEF_LinAlg
   Use m_MEF_Parameters
   Use m_MEF_Elements
   Use m_MEF_Sieve
   Use m_MEF_Utils
   Use petsc
   IMPLICIT NONE

   Private   
   Public :: MEF90_MassMatrixAssembleSet
   
   Interface MEF90_MassMatrixAssembleSet
      Module Procedure MassMatrixAssembleSet_2DScal, MassMatrixAssembleSet_2DVect, MassMatrixAssembleSet_2DElast, &
                       MassMatrixAssembleSet_3DScal, MassMatrixAssembleSet_3DVect, MassMatrixAssembleSet_3DElast 
   End Interface MEF90_MassMatrixAssembleSet

Contains
#undef __FUNCT__
#define __FUNCT__ "MassMatrixAssembleSet_2DScal"
   Subroutine MassMatrixAssembleSet_2DScal(M,mesh,U,cellIS,scaling,elem,elemType,ierr)
      Type(Mat),Intent(IN)                              :: M
      Type(DM),Intent(IN)                               :: mesh
      Type(SectionReal),Intent(IN)                      :: U
      Type(IS),Intent(IN)                               :: cellIS
      PetscReal,Intent(IN),optional                     :: scaling
      Type(MEF90Element2D_Scal), Dimension(:), Pointer  :: elem
      Type(MEF90Element_Type),Intent(IN)                :: elemType
      PetscErrorCode,Intent(OUT)                        :: ierr
      
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscReal,Dimension(:,:),Pointer             :: MatElem
      PetscInt                                     :: iDoF1,iDoF2,iGauss
      PetscLogDouble                               :: flops
     
      flops = 0
      Allocate(MatElem(elemType%numDof,elemType%numDof),stat=ierr)

      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Do cell = 1,size(elem)      
         MatElem = 0.0_Kr
         Do iGauss = 1,size(elem(cell)%Gauss_C)
            Do iDoF1 = 1,elemType%numDof
               Do iDoF2 = 1,elemType%numDof
                  MatElem(iDoF2,iDoF1) = MatElem(iDoF2,iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                        (elem(cell)%BF(iDoF1,iGauss) * elem(cell)%BF(iDoF2,iGauss) )
               End Do
            End Do
         End Do
         If (present(scaling)) Then
            MatElem = MatElem * scaling
         End If
         Call DMmeshAssembleMatrix(M,mesh,U,cellID(cell),MatElem,ADD_VALUES,ierr);CHKERRQ(ierr)
      End Do
      flops = 3 * elemType%numDof**2 * size(elem(1)%Gauss_C) * size(cellID) 
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      If (present(scaling)) Then
         flops = size(cellID) * elemType%numDof**2
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      End If
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      DeAllocate(MatElem,stat=ierr)
   End Subroutine MassMatrixAssembleSet_2DScal
   
#undef __FUNCT__
#define __FUNCT__ "MassMatrixAssembleSet_2DVect"
   Subroutine MassMatrixAssembleSet_2DVect(M,mesh,U,iBlk,elem,elemType,BC,ierr)
      Type(Mat),Intent(IN)                         :: M
      Type(DM),Intent(IN)                          :: mesh
      Type(SectionReal),Intent(IN)                 :: U
      PetscInt,Intent(IN)                          :: iBlk
      Type(MEF90Element2D_Vect), Dimension(:), Pointer  :: elem
      Type(MEF90Element_Type),Intent(IN)                :: elemType
      Type(Flag),Intent(IN),Optional               :: BC
      
      Type(IS)                                     :: cellIS
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscInt                                     :: numDof
      PetscInt                                     :: ierr
      PetscReal,Dimension(:,:),Pointer             :: MatElem
      PetscInt,Dimension(:),Pointer                :: BCFlag
      PetscInt                                     :: iDoF1,iDoF2,iGauss
      PetscLogDouble                               :: flops = 0
     
      numDof = elemType%numDof
      Allocate(MatElem(numDof,numDof))
      Allocate(BCFlag(numDof))
      BCFlag = 0

      Call DMmeshGetStratumIS(mesh,'Cell Sets',iBlk,cellIS,ierr); CHKERRQ(ierr)
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Do cell = 1,size(cellID)      
         MatElem = 0.0_Kr
         BCFlag = 0
         If (present(BC)) Then
            Call SectionIntRestrictClosure(BC%Sec,mesh,cellID(cell),numDof,BCFlag,ierr);CHKERRQ(ierr)
         End If
         Do iGauss = 1,size(elem(cellID(cell)+1)%Gauss_C)
            Do iDoF1 = 1,numDof
               If (BCFlag(iDoF1) == 0) Then
                  Do iDoF2 = 1,numDof
                     MatElem(iDoF2,iDoF1) = MatElem(iDoF2,iDoF1) + elem(cellID(cell)+1)%Gauss_C(iGauss) * &
                                           (elem(cellID(cell)+1)%BF(iDoF1,iGauss) .DotP. elem(cellID(cell)+1)%BF(iDoF2,iGauss) )
                  End Do
                  flops = flops + 2 * numDof
               End If
            End Do
         End Do
         Call DMmeshAssembleMatrix(M,mesh,U,cellID(cell),MatElem,ADD_VALUES,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Call ISDestroy(cellIS,ierr);CHKERRQ(ierr)
   
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      DeAllocate(BCFlag)
      DeAllocate(MatElem)
   End Subroutine MassMatrixAssembleSet_2DVect
   
#undef __FUNCT__
#define __FUNCT__ "MassMatrixAssembleSet_2DElast"
   Subroutine MassMatrixAssembleSet_2DElast(M,mesh,U,iBlk,elem,elemType,BC)
      Type(Mat),Intent(IN)                         :: M
      Type(DM),Intent(IN)                          :: mesh
      Type(SectionReal),Intent(IN)                 :: U
      PetscInt,Intent(IN)                          :: iBlk
      Type(MEF90Element2D_Elast), Dimension(:), Pointer :: elem
      Type(MEF90Element_Type),Intent(IN)                :: elemType
      Type(Flag),Intent(IN),Optional               :: BC
      
      Type(IS)                                     :: cellIS
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscInt                                     :: numDof
      PetscInt                                     :: ierr
      PetscReal,Dimension(:,:),Pointer             :: MatElem
      PetscInt,Dimension(:),Pointer                :: BCFlag
      PetscInt                                     :: iDoF1,iDoF2,iGauss
      PetscLogDouble                               :: flops = 0
     
      numDof = elemType%numDof
      Allocate(MatElem(numDof,numDof))
      Allocate(BCFlag(numDof))
      BCFlag = 0

      Call DMmeshGetStratumIS(mesh,'Cell Sets',iBlk,cellIS,ierr); CHKERRQ(ierr)
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Do cell = 1,size(cellID)      
         MatElem = 0.0_Kr
         BCFlag = 0
         If (present(BC)) Then
            Call SectionIntRestrictClosure(BC%Sec,mesh,cellID(cell),numDof,BCFlag,ierr);CHKERRQ(ierr)
         End If
         Do iGauss = 1,size(elem(cellID(cell)+1)%Gauss_C)
            Do iDoF1 = 1,numDof
               If (BCFlag(iDoF1) == 0) Then
                  Do iDoF2 = 1,numDof
                     MatElem(iDoF2,iDoF1) = MatElem(iDoF2,iDoF1) + elem(cellID(cell)+1)%Gauss_C(iGauss) * &
                                           (elem(cellID(cell)+1)%BF(iDoF1,iGauss) .DotP. elem(cellID(cell)+1)%BF(iDoF2,iGauss) )
                  End Do
                  flops = flops + 2 * numDof
               End If
            End Do
         End Do
         Call DMmeshAssembleMatrix(M,mesh,U,cellID(cell),MatElem,ADD_VALUES,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Call ISDestroy(cellIS,ierr);CHKERRQ(ierr)
   
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      DeAllocate(BCFlag)
      DeAllocate(MatElem)
   End Subroutine MassMatrixAssembleSet_2DElast
   
#undef __FUNCT__
#define __FUNCT__ "MassMatrixAssembleSet_3DScal"
   Subroutine MassMatrixAssembleSet_3DScal(M,mesh,U,iBlk,elem,elemType,BC)
      Type(Mat),Intent(IN)                         :: M
      Type(DM),Intent(IN)                          :: mesh
      Type(SectionREal),Intent(IN)                 :: U
      PetscInt,Intent(IN)                          :: iBlk
      Type(MEF90Element3D_Scal), Dimension(:), Pointer  :: elem
      Type(MEF90Element_Type),Intent(IN)                :: elemType
      Type(Flag),Intent(IN),Optional               :: BC
      
      Type(IS)                                     :: cellIS
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscInt                                     :: numDof
      PetscInt                                     :: ierr
      PetscReal,Dimension(:,:),Pointer             :: MatElem
      PetscInt,Dimension(:),Pointer                :: BCFlag
      PetscInt                                     :: iDoF1,iDoF2,iGauss
      PetscLogDouble                               :: flops = 0
     
      numDof = elemType%numDof
      Allocate(MatElem(numDof,numDof))
      Allocate(BCFlag(numDof))
      BCFlag = 0

      Call DMmeshGetStratumIS(mesh,'Cell Sets',iBlk,cellIS,ierr); CHKERRQ(ierr)
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Do cell = 1,size(cellID)      
         MatElem = 0.0_Kr
         BCFlag = 0
         If (present(BC)) Then
            Call SectionIntRestrictClosure(BC%Sec,mesh,cellID(cell),numDof,BCFlag,ierr);CHKERRQ(ierr)
         End If
         Do iGauss = 1,size(elem(cellID(cell)+1)%Gauss_C)
            Do iDoF1 = 1,numDof
               If (BCFlag(iDoF1) == 0) Then
                  Do iDoF2 = 1,numDof
                     MatElem(iDoF2,iDoF1) = MatElem(iDoF2,iDoF1) + elem(cellID(cell)+1)%Gauss_C(iGauss) * &
                                           (elem(cellID(cell)+1)%BF(iDoF1,iGauss) * elem(cellID(cell)+1)%BF(iDoF2,iGauss) )
                  End Do
                  flops = flops + 3 * numDof
               End If
            End Do
         End Do
         Call DMmeshAssembleMatrix(M,mesh,U,cellID(cell),MatElem,ADD_VALUES,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Call ISDestroy(cellIS,ierr);CHKERRQ(ierr)
   
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      DeAllocate(BCFlag)
      DeAllocate(MatElem)
   End Subroutine MassMatrixAssembleSet_3DScal
   
#undef __FUNCT__
#define __FUNCT__ "MassMatrixAssembleSet_3DVect"
   Subroutine MassMatrixAssembleSet_3DVect(M,mesh,U,iBlk,elem,elemType,BC)
      Type(Mat),Intent(IN)                         :: M
      Type(DM),Intent(IN)                          :: mesh
      Type(SectionReal),Intent(IN)                 :: U
      PetscInt,Intent(IN)                          :: iBlk
      Type(MEF90Element3D_Vect), Dimension(:), Pointer  :: elem
      Type(MEF90Element_Type),Intent(IN)                :: elemType
      Type(Flag),Intent(IN),Optional               :: BC
      
      Type(IS)                                     :: cellIS
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscInt                                     :: numDof
      PetscInt                                     :: ierr
      PetscReal,Dimension(:,:),Pointer             :: MatElem
      PetscInt,Dimension(:),Pointer                :: BCFlag
      PetscInt                                     :: iDoF1,iDoF2,iGauss
      PetscLogDouble                               :: flops = 0
     
      numDof = elemType%numDof
      Allocate(MatElem(numDof,numDof))
      Allocate(BCFlag(numDof))
      BCFlag = 0

      Call DMmeshGetStratumIS(mesh,'Cell Sets',iBlk,cellIS,ierr); CHKERRQ(ierr)
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Do cell = 1,size(cellID)      
         MatElem = 0.0_Kr
         BCFlag = 0
         If (present(BC)) Then
            Call SectionIntRestrictClosure(BC%Sec,mesh,cellID(cell),numDof,BCFlag,ierr);CHKERRQ(ierr)
         End If
         Do iGauss = 1,size(elem(cellID(cell)+1)%Gauss_C)
            Do iDoF1 = 1,numDof
               If (BCFlag(iDoF1) == 0) Then
                  Do iDoF2 = 1,numDof
                     MatElem(iDoF2,iDoF1) = MatElem(iDoF2,iDoF1) + elem(cellID(cell)+1)%Gauss_C(iGauss) * &
                                           (elem(cellID(cell)+1)%BF(iDoF1,iGauss) .DotP. elem(cellID(cell)+1)%BF(iDoF2,iGauss) )
                  End Do
                  flops = flops + 2 * numDof
               End If
            End Do
         End Do
         Call DMmeshAssembleMatrix(M,mesh,U,cellID(cell),MatElem,ADD_VALUES,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Call ISDestroy(cellIS,ierr);CHKERRQ(ierr)
   
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      DeAllocate(BCFlag)
      DeAllocate(MatElem)
   End Subroutine MassMatrixAssembleSet_3DVect
   
#undef __FUNCT__
#define __FUNCT__ "MassMatrixAssembleSet_3DElast"
   Subroutine MassMatrixAssembleSet_3DElast(M,mesh,U,iBlk,elem,elemType,BC)
      Type(Mat),Intent(IN)                         :: M
      Type(DM),Intent(IN)                          :: mesh
      Type(SectionReal),Intent(IN)                 :: U
      PetscInt,Intent(IN)                          :: iBlk
      Type(MEF90Element3D_Elast), Dimension(:), Pointer :: elem
      Type(MEF90Element_Type),Intent(IN)                :: elemType
      Type(Flag),Intent(IN),Optional               :: BC
      
      Type(IS)                                     :: cellIS
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscInt                                     :: numDof
      PetscInt                                     :: ierr
      PetscReal,Dimension(:,:),Pointer             :: MatElem
      PetscInt,Dimension(:),Pointer                :: BCFlag
      PetscInt                                     :: iDoF1,iDoF2,iGauss
      PetscLogDouble                               :: flops = 0
     
      numDof = elemType%numDof
      Allocate(MatElem(numDof,numDof))
      Allocate(BCFlag(numDof))
      BCFlag = 0

      Call DMmeshGetStratumIS(mesh,'Cell Sets',iBlk,cellIS,ierr); CHKERRQ(ierr)
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Do cell = 1,size(cellID)      
         MatElem = 0.0_Kr
         BCFlag = 0
         If (present(BC)) Then
            Call SectionIntRestrictClosure(BC%Sec,mesh,cellID(cell),numDof,BCFlag,ierr);CHKERRQ(ierr)
         End If
         Do iGauss = 1,size(elem(cellID(cell)+1)%Gauss_C)
            Do iDoF1 = 1,numDof
               If (BCFlag(iDoF1) == 0) Then
                  Do iDoF2 = 1,numDof
                     MatElem(iDoF2,iDoF1) = MatElem(iDoF2,iDoF1) + elem(cellID(cell)+1)%Gauss_C(iGauss) * &
                                           (elem(cellID(cell)+1)%BF(iDoF1,iGauss) .DotP. elem(cellID(cell)+1)%BF(iDoF2,iGauss) )
                  End Do
                  flops = flops + 2 * numDof
               End If
            End Do
         End Do
         Call DMmeshAssembleMatrix(M,mesh,U,cellID(cell),MatElem,ADD_VALUES,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Call ISDestroy(cellIS,ierr);CHKERRQ(ierr)
   
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      DeAllocate(BCFlag)
      DeAllocate(MatElem)
   End Subroutine MassMatrixAssembleSet_3DElast
End Module m_MEF_MassMatrix
