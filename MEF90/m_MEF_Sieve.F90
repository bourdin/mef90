Module m_MEF_Sieve
#include "finclude/petscdef.h"

   Use m_MEF_LinAlg
   Use m_MEF_Elements
   Use petsc
         
   IMPLICIT NONE
   Private
   
   Public :: SectionRealSetVertexFiberDimensionSet,SectionRealSetCellFiberDimensionSet
   Public :: Field
   Public :: Flag
   Public :: FieldCreateVertex
   Public :: FlagCreateVertex
   Public :: FieldDestroy
   Public :: FlagDestroy
   Public :: MatInsertVertexBoundaryValues
   Public :: FieldInsertVertexBoundaryValues   
   
   Type Field
      !!! A Field is a general container containing a section,a vector of sections for its components
      !!! a (global) vec,a local Vec,and the associated scatter.
      !!! All the Vec and SectionReal share the same memory space for their data storage
      Type(SectionReal)                            :: Sec
      Type(SectionReal),Dimension(:),Pointer       :: Component_Sec
      Type(Vec)                                    :: Vec
      Type(Vec)                                    :: LocalVec
      Type(VecScatter)                             :: Scatter
      PetscInt                                     :: num_components
      PetscInt,Dimension(:),Pointer                :: component_size
   End Type Field
   
   Type Flag
      !!! a Flag is to field what a SectionInt is to a SectionReal...
      !!! except that there is not VecInt,so the whole Vec and scatter part is ignored.
      Type(SectionInt)                             :: Sec
      Type(SectionInt),Dimension(:),Pointer        :: Component_Sec
      PetscInt                                     :: num_components
      PetscInt,Dimension(:),Pointer                :: component_size
   End Type Flag

Contains
#undef __FUNCT__
#define __FUNCT__ "SectionRealSetVertexFiberDimensionSet"
!!!
!!!  
!!!  SectionRealSetVertexFiberDimensionSet:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine SectionRealSetVertexFiberDimensionSet(Mesh,section,setIS,numDoF,ierr)
   Type(DM),Intent(IN)                             :: Mesh
   Type(SectionReal),Intent(IN)                    :: section
   Type(IS),Intent(IN)                             :: setIS
   PetscInt,Intent(IN)                             :: numDoF
   PetscErrorCode,Intent(OUT)                      :: ierr

   PetscInt                                        :: cell,point
   PetscInt,Dimension(:),Pointer                   :: cellID,cone
   
   Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
   Do cell = 1, size(cellID)
      Call DMMeshGetConeF90(mesh,cellID(cell),Cone,ierr);CHKERRQ(ierr)
      Do point = 1, size(cone)
         Call SectionRealSetFiberDimension(section,Cone(point),numDof,ierr);CHKERRQ(ierr)
      End Do
      Call DMMeshRestoreConeF90(mesh,cellID(cell),Cone,ierr);CHKERRQ(ierr)
   End Do
   Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
End Subroutine SectionRealSetVertexFiberDimensionSet

#undef __FUNCT__
#define __FUNCT__ "SectionRealSetCellFiberDimensionSet"
!!!
!!!  
!!!  SectionRealSetCellFiberDimensionSet:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine SectionRealSetCellFiberDimensionSet(Mesh,section,setIS,numDoF,ierr)
   Type(DM),Intent(IN)                             :: Mesh
   Type(SectionReal),Intent(IN)                    :: section
   Type(IS),Intent(IN)                             :: setIS
   PetscInt,Intent(IN)                             :: numDoF
   PetscErrorCode,Intent(OUT)                      :: ierr

   PetscInt                                        :: cell
   PetscInt,Dimension(:),Pointer                   :: cellID
   
   Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
   Do cell = 1, size(cellID)
      Call SectionRealSetFiberDimension(section,cellID(cell),numDof,ierr);CHKERRQ(ierr)
   End Do
   Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
End Subroutine SectionRealSetCellFiberDimensionSet

#undef __FUNCT__
#define __FUNCT__ "FieldCreateVertex"
   Subroutine FieldCreateVertex(F,Fname,mesh,component_size)
      Type(Field)                                  :: F
      Character(len=*)                             :: Fname
      Type(DM)                                     :: mesh
      PetscInt,Dimension(:),Pointer                :: component_size
      Character(len=256)                           :: component_name
      PetscInt                                     :: numCells,numVertices
      PetscInt                                     :: i,j,ierr
      
      Call DMMeshGetStratumSize(mesh,"height",0,numCells,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumSize(mesh,"depth",0,numVertices,ierr);CHKERRQ(ierr)

      F%num_components = Size(component_size)

      !!! Create the Main Section
      Call DMMeshGetVertexSectionReal(mesh,Fname,sum(component_size),F%Sec,ierr);CHKERRQ(ierr)
      
      !!! Add space for each of the individual component section
      Allocate(F%Component_size(F%num_components))
      Do i = 1,F%num_components
         Call SectionRealAddSpace(F%Sec,ierr);CHKERRQ(ierr)
         F%component_size(i) = component_size(i)
      End Do 
      Call SectionRealAllocate(F%Sec,ierr);CHKERRQ(ierr)
      
      !!! Set the fibration size
      Do i = 1,numVertices
         Do j = 1,F%num_components
            Call SectionRealSetFiberDimensionField(F%Sec,i+numCells-1,component_size(j),j-1,ierr);CHKERRQ(ierr)
         End Do
      End Do 
      
      !!! Create the individual component sections
      Allocate(F%Component_Sec(F%num_Components))
      Do i = 1,F%num_Components      
         Call SectionRealGetFibration(F%Sec,i-1,F%Component_sec(i),ierr);CHKERRQ(ierr)
         Write(component_name,"(A,'.',I3.3)") Trim(Fname),i
         Call PetscObjectSetName(F%Component_Sec(i),component_name,ierr);CHKERRQ(ierr)
      End Do
      !!! Create the Scatter and global and local Vec
!      F%Has_Vec  = .TRUE.
      Call DMMeshCreateGlobalScatter(mesh,F%Sec,F%Scatter,ierr);CHKERRQ(ierr)
      Call DMMeshCreateVector(mesh,F%Sec,F%Vec,ierr);CHKERRQ(ierr)
      Call SectionRealCreateLocalVector(F%Sec,F%LocalVec,ierr);CHKERRQ(ierr)
      Call VecSetBlockSize(F%LocalVec,sum(component_size),ierr);CHKERRQ(ierr)
   End Subroutine FieldCreateVertex
   
#undef __FUNCT__
#define __FUNCT__ "FieldDestroy"
   Subroutine FieldDestroy(F)
      Type(Field)                                  :: F
      PetscInt                                     :: i,ierr
      
      DeAllocate(F%Component_size)
      Call SectionRealDestroy(F%Sec,ierr);CHKERRQ(ierr)
      Do i = 1,F%Num_Components   
         Call SectionRealDestroy(F%Component_Sec(i),ierr);CHKERRQ(ierr)   
      End Do
      DeAllocate(F%Component_Sec)
      
      Call VecDestroy(F%Vec,ierr);CHKERRQ(ierr)
      Call VecDestroy(F%LocalVec,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(F%Scatter,ierr);CHKERRQ(ierr)
   End Subroutine FieldDestroy

#undef __FUNCT__
#define __FUNCT__ "FlagCreateVertex"
   Subroutine FlagCreateVertex(F,Fname,mesh,component_size)
      !!! creates a Flag derived type with sum(component_size) values at each vertex
      !!! add a way to scatter it into size(component_size) component sections
      !!! I am not sure this actualy works and does anything meaningful if component_size(i) /= 1...
      Type(Flag)                                   :: F
      Character(len=*)                             :: Fname
      Type(DM)                                     :: mesh
      PetscInt,Dimension(:),Pointer                :: component_size
      Character(len=256)                           :: component_name
      PetscInt                                     :: numCells,numVertices
      PetscInt                                     :: i,j,ierr
      
      Call DMMeshGetStratumSize(mesh,"height",0,numCells,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumSize(mesh,"depth",0,numVertices,ierr);CHKERRQ(ierr)

      F%num_components = Size(component_size)

      !!! Create the Main Section
      Call DMMeshGetVertexSectionInt(Mesh,Fname,sum(component_size),F%Sec,ierr);CHKERRQ(ierr)
      
      !!! Add space for each of the individual component section
      Allocate(F%Component_size(F%num_components))
      Do i = 1,F%num_components
         Call SectionIntAddSpace(F%Sec,ierr);CHKERRQ(ierr)
         F%component_size(i) = component_size(i)
      End Do 
      Call SectionIntAllocate(F%Sec,ierr);CHKERRQ(ierr)
      
      !!! Set the fibration size
      Do i = 1,numVertices
         Do j = 1,F%num_components
            Call SectionIntSetFiberDimensionField(F%Sec,i+numCells-1,component_size(j),j-1,ierr);CHKERRQ(ierr)
         End Do
      End Do 
      
      !!! Create the individual component sections
      Allocate(F%Component_Sec(F%num_Components))
      Do i = 1,F%num_Components      
         Call SectionIntGetFibration(F%Sec,i-1,F%Component_sec(i),ierr);CHKERRQ(ierr)
         Write(component_name,"(A,'.',I3.3)") Trim(Fname),i
         Call PetscObjectSetName(F%Component_Sec(i),component_name,ierr);CHKERRQ(ierr)
      End Do
   End Subroutine FlagCreateVertex

#undef __FUNCT__
#define __FUNCT__ "FlagDestroy"
   Subroutine FlagDestroy(F)
      Type(Flag)                                   :: F
      PetscInt                                     :: i,ierr
      
      Call SectionIntDestroy(F%Sec,ierr);CHKERRQ(ierr)
      DeAllocate(F%Component_size)
      Do i = 1,F%Num_Components   
         Call SectionIntDestroy(F%Component_Sec(i),ierr);CHKERRQ(ierr)   
      End Do
      DeAllocate(F%Component_Sec)
   End Subroutine FlagDestroy
   
#undef __FUNCT__
#define __FUNCT__ "MatInsertVertexBoundaryValues"
   Subroutine MatInsertVertexBoundaryValues(M,U,BCFlag,mesh)
      Type(Mat)                                    :: M
      Type(Field)                                  :: U
      Type(Flag)                                   :: BCFlag
      Type(DM)                                     :: mesh
      
      PetscInt                                     :: ierr
      PetscInt,Dimension(:),Pointer                :: BCFlag_Ptr
      PetscReal,Dimension(:,:),Pointer             :: MatElem
      PetscInt                                     :: i,j,num_dof,zero,numCells,numVertices

      !!!
      !!! DMMeshAssembleMatrix does not work with fibrated sections at this point
      !!! in order to INSERT boundary values,I need to insert a block for ALL dof associated to a given point
      !!! therefre erasing exsting values....
      !!! MatInsertBoundaryValues needs to be called BEFORE building the hessian or stiffness matrix
      !!!
      Call DMMeshGetStratumSize(mesh,"height",0,numCells,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumSize(mesh,"depth",0,numVertices,ierr);CHKERRQ(ierr)
      zero = 0
      num_dof = sum(BCFlag%component_size)
      Allocate(MatElem(num_dof,num_dof))  
      Do i = 1,numVertices
         Call SectionIntRestrict(BCFlag%Sec,numCells+i-1,BCFlag_Ptr,ierr);CHKERRQ(ierr)
         If (Sum(BCFlag_Ptr) /= 0) Then
            MatElem = 0.0_Kr
            Do j = 1,num_dof
               If (BCFlag_Ptr(j) /= zero) Then
                  MatElem(j,j) = 1.0_Kr
               End If  
            End Do
            Call DMMeshAssembleMatrix(M,mesh,U%Sec,numCells+i-1,MatElem,INSERT_VALUES,ierr);CHKERRQ(ierr)
         End If
      End Do
      DeAllocate(MatElem)
   End Subroutine MatInsertVertexBoundaryValues

#undef __FUNCT__
#define __FUNCT__ "FieldInsertVertexBoundaryValues"
   Subroutine FieldInsertVertexBoundaryValues(F,FBC,BCFlag,mesh)
      !!! Insert in F the value of FBC at each point where BCFlag is non 0
      Type(Field)                                  :: F,FBC
      Type(Flag)                                   :: BCFlag
      Type(DM)                                     :: mesh
      
      PetscInt                                     :: ierr,zero
      PetscInt,Dimension(:),Pointer                :: BCFlag_Ptr
      PetscReal,Dimension(:),Pointer               :: FBC_Ptr
      PetscInt                                     :: i,j,numCells,numVertices

      Call DMMeshGetStratumSize(mesh,"height",0,numCells,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumSize(mesh,"depth",0,numVertices,ierr);CHKERRQ(ierr)
      zero = 0
      Do j = 1,BCFlag%num_components 
         If (BCFlag%Component_size(j) /= 1 ) Then
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,'FieldInsertVertexBoundaryValues requires scalar components',ierr)
         End If
         Do i = 1,numVertices
            Call SectionIntRestrict(BCFlag%Component_Sec(j),numCells+i-1,BCFlag_Ptr,ierr);CHKERRQ(ierr)
            If (BCFlag_Ptr(1) /= zero) Then
               Call SectionRealRestrict(FBC%Component_Sec(j),numCells+i-1,FBC_Ptr,ierr);CHKERRQ(ierr)
               Call SectionRealUpdate  (F%Component_Sec(j),numCells+i-1,FBC_Ptr,INSERT_VALUES,ierr);CHKERRQ(ierr)
               Call SectionRealRestore (FBC%Component_Sec(j),numCells+i-1,FBC_Ptr,ierr);CHKERRQ(ierr)
            End If
            Call SectionIntRestore(BCFlag%Component_Sec(j),numCells+i-1,BCFlag_Ptr,ierr);CHKERRQ(ierr)
         End Do
      End Do
   End Subroutine FieldInsertVertexBoundaryValues
End Module m_MEF_Sieve
