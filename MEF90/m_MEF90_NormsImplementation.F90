#include "mef90.inc"
Module MEF90_APPEND(m_MEF90_NormsImplementation_,MEF90_ELEMENTTYPE)

#include "petsc/finclude/petsc.h"
   Use m_MEF90_LinAlg
   Use m_MEF90_Parameters
   Use m_MEF90_Elements
   Use m_MEF90_Utils
   Use m_MEF90_DMPlex
   Use petsc
   IMPLICIT NONE

   Private   
#if MEF90_ELEMENTTYPE_SCALAR
   Public :: MEF90_L2DotProductSet,MEF90_H1DotProductSet,MEF90_L2NormSet
#else
   Public :: MEF90_L2DotProductSet,MEF90_H1DotProductSet,MEF90_H1symDotProductSet,MEF90_L2NormSet
#endif

Contains
!!!
!!!  
!!!  MEF90_L2DotProductSet: Assemble and add the contribution one processor to the L2 dot products of Vec U and V 
!!!                   on a cell / face / edge set, interpolated with element type elemType
!!!  
!!!  (c) 2022 Blaise Bourdin bourdin@mcmaster.ca
!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90_L2DotProductSet"

   Subroutine MEF90_L2DotProductSet(myDotProductSet,U,V,setType,setID,elem,elemType,ierr)
      PetscReal,Intent(OUT)                           :: myDotProductSet
      Type(tVec),Intent(IN)                           :: U,V
      PetscEnum,Intent(IN)                            :: setType
      PetscInt                                        :: setID
      Type(MEF90_ELEMENTTYPE), Dimension(:), Pointer  :: elem
      Type(MEF90Element_Type),Intent(IN)              :: elemType
      PetscErrorCode,Intent(INOUT)                    :: ierr

      PetscReal,Dimension(:),Pointer                  :: Uloc,Vloc
#if MEF90_ELEMENTTYPE_SCALAR
      PetscReal                                       :: UGauss,Vgauss
#else
      Type(MEF90_VECT)                                :: UGauss,VGauss
#endif
      Type(tDM)                                       :: dm
      Type(tIS)                                       :: setPointIS
      PetscInt,Dimension(:),Pointer                   :: setPointID
      PetscInt                                        :: point
      PetscInt                                        :: iDoF1,numDof
      PetscInt                                        :: iGauss,numGauss
      !PetscLogDouble                                  :: flops

      !flops = 0.0_pflop
      PetscCall(VecGetDM(U,dm,ierr))
      PetscCall(DMGetStratumIS(dm,MEF90_DMPlexSetLabelName(setType),setID,setPointIS,ierr))
      PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
      If (size(setPointID) > 0) Then
         !!! This is really misleading: elemType doesn't know the number of component since we now use the 
         !!! same elemType for scalar and Vect elements, elem%numDof is NOT the number of dof...
         !!! Maybe I need to change this to the old behaviour.
         numDof = size(elem(1)%BF(:,1))
         numGauss = size(elem(1)%Gauss_C)
         Allocate(Uloc(numDof),source = 0.0_Kr)
         Do point = 1,size(setPointID)
            PetscCall(DMPlexVecGetClosure(dm,PETSC_NULL_SECTION,U,point,Uloc,ierr))
            PetscCall(DMPlexVecGetClosure(dm,PETSC_NULL_SECTION,V,point,Vloc,ierr))
            Do iGauss = 1,size(elem(point)%Gauss_C)
               UGauss = 0.0_Kr
               VGauss = 0.0_Kr
               Do iDoF1 = 1,numDof
                     UGauss = UGauss + elem(point)%BF(iDoF1,iGauss) * Uloc(iDof1)
                     VGauss = VGauss + elem(point)%BF(iDoF1,iGauss) * Vloc(iDof1)
               End Do ! iDoF1
               myDotProductSet = myDotProductSet + (UGauss * VGauss) * elem(point)%Gauss_C(iGauss)
            End Do ! iGauss
            PetscCall(DMPlexVecRestoreClosure(dm,PETSC_NULL_SECTION,V,point,Vloc,ierr))
            PetscCall(DMPlexVecRestoreClosure(dm,PETSC_NULL_SECTION,U,point,Uloc,ierr))
         End Do ! point
         !flops = 3 * numDof**2 * numGauss * size(setPointID)
         !PetscCall(PetscLogFlops(flops,ierr))
         ! Flop computation is different for scalar and Vec
         DeAllocate(Uloc,stat=ierr)
      End If 
      PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
      PetscCall(ISDestroy(setPointIS,ierr))
   End Subroutine MEF90_L2DotProductSet   

!!!
!!!  
!!!  MEF90_H1DotProductSet: Assemble and add the contribution one processor to the H1 dot products of Vec U and V 
!!!                   on a cell / face / edge set, interpolated with element type elemType
!!!  
!!!  (c) 2022 Blaise Bourdin bourdin@mcmaster.ca
!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90_H1DotProductSet"
   
   Subroutine MEF90_H1DotProductSet(myDotProductSet,U,V,setType,setID,elem,elemType,ierr)
      PetscReal,Intent(OUT)                           :: myDotProductSet
      Type(tVec),Intent(IN)                           :: U,V
      PetscEnum,Intent(IN)                            :: setType
      PetscInt                                        :: setID
      Type(MEF90_ELEMENTTYPE), Dimension(:), Pointer  :: elem
      Type(MEF90Element_Type),Intent(IN)              :: elemType
      PetscErrorCode,Intent(INOUT)                    :: ierr

      PetscReal,Dimension(:),Pointer                  :: Uloc,Vloc
#if MEF90_ELEMENTTYPE_SCALAR
      Type(MEF90_VECT)                                :: GradUGauss,GradVgauss
#else
      Type(MEF90_MAT)                                 :: GradUGauss,GradVGauss
#endif
      Type(tDM)                                       :: dm
      Type(tIS)                                       :: setPointIS
      PetscInt,Dimension(:),Pointer                   :: setPointID
      PetscInt                                        :: point
      PetscInt                                        :: iDoF1,numDof
      PetscInt                                        :: iGauss,numGauss
      !PetscLogDouble                                  :: flops

      !flops = 0.0_pflop
      PetscCall(VecGetDM(U,dm,ierr))
      PetscCall(DMGetStratumIS(dm,MEF90_DMPlexSetLabelName(setType),setID,setPointIS,ierr))
      PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
      If (size(setPointID) > 0) Then
         !!! This is really misleading: elemType doesn't know the number of component since we now use the 
         !!! same elemType for scalar and Vect elements, elem%numDof is NOT the number of dof...
         !!! Maybe I need to change this to the old behaviour.
         numDof = size(elem(1)%BF(:,1))
         numGauss = size(elem(1)%Gauss_C)
         Allocate(Uloc(numDof),source = 0.0_Kr)
         Do point = 1,size(setPointID)
            PetscCall(DMPlexVecGetClosure(dm,PETSC_NULL_SECTION,U,point,Uloc,ierr))
            PetscCall(DMPlexVecGetClosure(dm,PETSC_NULL_SECTION,V,point,Vloc,ierr))
            Do iGauss = 1,size(elem(point)%Gauss_C)
               GradUGauss = 0.0_Kr
               GradVGauss = 0.0_Kr
               Do iDoF1 = 1,numDof
                     GradUGauss = GradUGauss + elem(point)%Grad_BF(iDoF1,iGauss) * Uloc(iDof1)
                     GradVGauss = GradVGauss + elem(point)%Grad_BF(iDoF1,iGauss) * Vloc(iDof1)
               End Do ! iDoF1
               myDotProductSet = myDotProductSet + (GradUGauss .DotP. GradVGauss) * elem(point)%Gauss_C(iGauss)
            End Do ! iGauss
            PetscCall(DMPlexVecRestoreClosure(dm,PETSC_NULL_SECTION,V,point,Vloc,ierr))
            PetscCall(DMPlexVecRestoreClosure(dm,PETSC_NULL_SECTION,U,point,Uloc,ierr))
         End Do ! point
         !flops = 3 * numDof**2 * numGauss * size(setPointID)
         !PetscCall(PetscLogFlops(flops,ierr))
         ! Flop computation is different for scalar and Vec
         DeAllocate(Uloc,stat=ierr)
      End If 
      PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
      PetscCall(ISDestroy(setPointIS,ierr))
   End Subroutine MEF90_H1DotProductSet   
   
!!!
!!!  
!!!  MEF90_H1SymDotProductSet: Assemble and add the contribution one processor to the H1-sym dot products of Vec U and V 
!!!                   on a cell / face / edge set, interpolated with element type elemType
!!!  
!!!  (c) 2022 Blaise Bourdin bourdin@mcmaster.ca
!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90_H1SymDotProductSet"
#if MEF90_ELEMENTTYPE_VECT
   
   Subroutine MEF90_H1SymDotProductSet(myDotProductSet,U,V,setType,setID,elem,elemType,ierr)
      PetscReal,Intent(OUT)                           :: myDotProductSet
      Type(tVec),Intent(IN)                           :: U,V
      PetscEnum,Intent(IN)                            :: setType
      PetscInt                                        :: setID
      Type(MEF90_ELEMENTTYPE), Dimension(:), Pointer  :: elem
      Type(MEF90Element_Type),Intent(IN)              :: elemType
      PetscErrorCode,Intent(INOUT)                    :: ierr

      PetscReal,Dimension(:),Pointer                  :: Uloc,Vloc
      Type(MEF90_MATS)                                :: GradSUGauss,GradSVGauss
      Type(tDM)                                       :: dm
      Type(tIS)                                       :: setPointIS
      PetscInt,Dimension(:),Pointer                   :: setPointID
      PetscInt                                        :: point
      PetscInt                                        :: iDoF1,numDof
      PetscInt                                        :: iGauss,numGauss
      !PetscLogDouble                                  :: flops

      !flops = 0.0_pflop
      PetscCall(VecGetDM(U,dm,ierr))
      PetscCall(DMGetStratumIS(dm,MEF90_DMPlexSetLabelName(setType),setID,setPointIS,ierr))
      PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
      If (size(setPointID) > 0) Then
         !!! This is really misleading: elemType doesn't know the number of component since we now use the 
         !!! same elemType for scalar and Vect elements, elem%numDof is NOT the number of dof...
         !!! Maybe I need to change this to the old behaviour.
         numDof = size(elem(1)%BF(:,1))
         numGauss = size(elem(1)%Gauss_C)
         Allocate(Uloc(numDof),source = 0.0_Kr)
         Do point = 1,size(setPointID)
            PetscCall(DMPlexVecGetClosure(dm,PETSC_NULL_SECTION,U,point,Uloc,ierr))
            PetscCall(DMPlexVecGetClosure(dm,PETSC_NULL_SECTION,V,point,Vloc,ierr))
            Do iGauss = 1,size(elem(point)%Gauss_C)
               GradSUGauss = 0.0_Kr
               GradSVGauss = 0.0_Kr
               Do iDoF1 = 1,numDof
                     GradSUGauss = GradSUGauss + elem(point)%GradS_BF(iDoF1,iGauss) * Uloc(iDof1)
                     GradSVGauss = GradSVGauss + elem(point)%GradS_BF(iDoF1,iGauss) * Vloc(iDof1)
               End Do ! iDoF1
               myDotProductSet = myDotProductSet + (GradSUGauss .DotP. GradSVGauss) * elem(point)%Gauss_C(iGauss)
            End Do ! iGauss
            PetscCall(DMPlexVecRestoreClosure(dm,PETSC_NULL_SECTION,V,point,Vloc,ierr))
            PetscCall(DMPlexVecRestoreClosure(dm,PETSC_NULL_SECTION,U,point,Uloc,ierr))
         End Do ! point
         !flops = 3 * numDof**2 * numGauss * size(setPointID)
         !PetscCall(PetscLogFlops(flops,ierr))
         ! Flop computation is different for scalar and Vec
         DeAllocate(Uloc,stat=ierr)
      End If 
      PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
      PetscCall(ISDestroy(setPointIS,ierr))
   End Subroutine MEF90_H1SymDotProductSet   
#endif
!!!
!!!  
!!!  MEF90_L2NormSet: Assemble and add the contribution one processor to the L2 norm of a Vect U 
!!!                   on a cell / face / edge set, interpolated with element type elemType
!!!  
!!!  (c) 2022 Blaise Bourdin bourdin@mcmaster.ca
!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90_L2NormSet"
   
      Subroutine MEF90_L2NormSet(myNormSet,U,setType,setID,elem,elemType,ierr)
         PetscReal,Intent(OUT)                           :: myNormSet
         Type(tVec),Intent(IN)                           :: U
         PetscEnum,Intent(IN)                            :: setType
         PetscInt                                        :: setID
         Type(MEF90_ELEMENTTYPE), Dimension(:), Pointer  :: elem
         Type(MEF90Element_Type),Intent(IN)              :: elemType
         PetscErrorCode,Intent(INOUT)                    :: ierr
   
         PetscReal,Dimension(:),Pointer                  :: Uloc
#if MEF90_ELEMENTTYPE_SCALAR
         PetscReal                                       :: UGauss
#else
         Type(MEF90_VECT)                                :: UGauss
#endif
         Type(tDM)                                       :: dm
         Type(tIS)                                       :: setPointIS
         PetscInt,Dimension(:),Pointer                   :: setPointID
         PetscInt                                        :: point
         PetscInt                                        :: iDoF1,numDof
         PetscInt                                        :: iGauss,numGauss
         !PetscLogDouble                                  :: flops
   
         !flops = 0.0_pflop
         PetscCall(VecGetDM(U,dm,ierr))
         PetscCall(DMGetStratumIS(dm,MEF90_DMPlexSetLabelName(setType),setID,setPointIS,ierr))
         PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
         If (size(setPointID) > 0) Then
            !!! This is really misleading: elemType doesn't know the number of component since we now use the 
            !!! same elemType for scalar and Vect elements, elem%numDof is NOT the number of dof...
            !!! Maybe I need to change this to the old behaviour.
            numDof = size(elem(1)%BF(:,1))
            numGauss = size(elem(1)%Gauss_C)
            Allocate(Uloc(numDof),source = 0.0_Kr)
            Do point = 1,size(setPointID)
               PetscCall(DMPlexVecGetClosure(dm,PETSC_NULL_SECTION,U,point,Uloc,ierr))
               Do iGauss = 1,size(elem(point)%Gauss_C)
                  UGauss = 0.0_Kr
                  Do iDoF1 = 1,numDof
                        UGauss = UGauss + elem(point)%BF(iDoF1,iGauss) * Uloc(iDof1)
                  End Do ! iDoF1
                  myNormSet = myNormSet + (UGauss*UGauss) * elem(point)%Gauss_C(iGauss)
               End Do ! iGauss
               PetscCall(DMPlexVecRestoreClosure(dm,PETSC_NULL_SECTION,U,point,Uloc,ierr))
            End Do ! point
            !flops = 3 * numDof**2 * numGauss * size(setPointID)
            !PetscCall(PetscLogFlops(flops,ierr))
            ! Flop computation is different for scalar and Vec
            DeAllocate(Uloc,stat=ierr)
         End If 
         PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
         PetscCall(ISDestroy(setPointIS,ierr))
      End Subroutine MEF90_L2NormSet   
End Module MEF90_APPEND(m_MEF90_NormsImplementation_,MEF90_ELEMENTTYPE)