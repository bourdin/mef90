#include "mef90.inc"
Module MEF90_APPEND(m_MEF90_MassMatrixImplementation_,MEF90_ELEMENTTYPE)
#include "petsc/finclude/petsc.h"
   Use m_MEF90_LinAlg
   Use m_MEF90_Parameters
   Use m_MEF90_Elements
   Use m_MEF90_Utils
   Use petsc
   IMPLICIT NONE

   Private   
   Public :: MEF90_MassMatrixAssembleSet   

Contains
!!!
!!!  
!!!  MEF90_MassMatrixAssembleSet: Assemble the contribution of a cell / face / edge set to 
!!!                               the mass matrix associated with the element type elemType
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022 Blaise Bourdin bourdin@mcmaster.ca
!!!
#undef __FUNCT__
#define __FUNCT__ "MEF90_MassMatrixAssembleSet"

   Subroutine MEF90_MassMatrixAssembleSet(M,dm,setType,setID,elem,elemType,ierr)
      Type(tMat),Intent(IN)                           :: M
      Type(tDM),Intent(IN)                            :: dm
      Character(len=*),Intent(IN)                     :: setType
      PetscInt                                        :: setID
      Type(MEF90_ELEMENTTYPE), Dimension(:), Pointer  :: elem
      Type(MEF90Element_Type),Intent(IN)              :: elemType
      PetscErrorCode,Intent(OUT)                      :: ierr
      
      Type(tIS)                                       :: setPointIS
      PetscInt,Dimension(:),Pointer                   :: setPointID
      PetscInt                                        :: point
      PetscReal,Dimension(:),Pointer                  :: MatElem
      PetscInt                                        :: iDoF1,iDoF2,iGauss
      PetscLogDouble                                  :: flops
     
      flops = 0.0_pflop
      PetscCall(DMGetStratumIS(dm,setType,setID,setPointIS,ierr))
      PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
      If (size(setPointID) > 0) Then
         Allocate(MatElem(elemType%numDof**2))
         Do point = 1,size(setPointID)   
            MatElem = 0.0_Kr
            Do iGauss = 1,size(elem(point)%Gauss_C)
               Do iDoF1 = 1,elemType%numDof
                  Do iDoF2 = 1,elemType%numDof
                     MatElem((iDoF1-1)*elemType%numDof+iDof2) = MatElem((iDoF1-1)*elemType%numDof+iDof2) &
                        + elem(point)%Gauss_C(iGauss) * ( elem(point)%BF(iDoF1,iGauss) * elem(point)%BF(iDoF2,iGauss) )
                  End Do ! iDoF2
               End Do ! iDoF1
            End Do ! iGauss
            PetscCall(DMPlexMatSetClosure(dm,PETSC_NULL_SECTION,PETSC_NULL_SECTION,M,setPointID(point),MatElem,ADD_VALUES,ierr))
         End Do ! point
         flops = 3 * elemType%numDof**2 * size(elem(1)%Gauss_C) * size(setPointID)
         PetscCall(PetscLogFlops(flops,ierr))
         DeAllocate(MatElem,stat=ierr)
      End If 
      PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
      PetscCall(ISDestroy(setPointIS,ierr))
   End Subroutine MEF90_MassMatrixAssembleSet   
End Module MEF90_APPEND(m_MEF90_MassMatrixImplementation_,MEF90_ELEMENTTYPE)