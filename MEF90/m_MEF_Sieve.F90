Module m_MEF_Sieve
#include "finclude/petscdef.h"

   Use m_MEF_LinAlg
   Use m_MEF_Types
   Use m_MEF_Utils
   Use petsc
         
   IMPLICIT NONE
   Private
   
   Public :: MeshTopologyGetInfo
   !Public :: MeshTopologyGetInfo2

   Public :: FieldCreateVertex
   Public :: FlagCreateVertex
   Public :: FieldDestroy
   Public :: FlagDestroy
   Public :: SectionIntAddNSProperty
   Public :: MatInsertVertexBoundaryValues
   Public :: FieldInsertVertexBoundaryValues   
   
Contains
#undef __FUNCT__
#define __FUNCT__ "MeshTopologyGetInfo"
!!! Same as before except that each cpu has the same number of blocks and sets
!!! with potentially 0 size. 
!!! 
   Subroutine MeshTopologyGetInfo(dMeshTopology)
      Type (MeshTopology_Type)                     :: dMeshTopology

      MPI_Comm                                     :: comm
      PetscErrorCode                               :: ierr
      Type(IS)                                     :: cellSetIS
      PetscInt,Dimension(:),Pointer                :: cellSetID
      PetscInt                                     :: set
      
      Call PetscObjectGetComm(dMeshTopology%mesh,comm,ierr);CHKERRQ(ierr)      
      !!!
      !!! Get global IS for vertex and cell sets.
      !!! This is required for synchronous loops over sets (like I/O)
      !!!
      Call DMMeshGetLabelIdIS(dMeshTopology%mesh,'Cell Sets',dMeshTopology%cellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(comm,dMeshTopology%cellSetGlobalIS)
      Call DMMeshGetLabelIdIS(dMeshTopology%mesh,'Vertex Sets',dMeshTopology%vertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(comm,dMeshTopology%vertexSetGlobalIS)
      
      Call DMMeshGetLabelIdIS(dMeshTopology%mesh,'Cell Sets',cellSetIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(cellSetIS,cellSetID,ierr);CHKERRQ(ierr)
      Allocate(dMeshTopology%cellSet(size(cellSetID)))
      Do set = 1, size(cellSetID)
         dMeshTopology%cellSet(set)%ElemType    = -1
         dMeshTopology%cellSet(set)%dofLocation = -1
         dMeshTopology%cellSet(set)%numDof      = -1
         dMeshTopology%cellSet(set)%coDimension = -1
      End Do
      Call ISRestoreIndicesF90(cellSetIS,cellSetID,ierr);CHKERRQ(ierr)
   End Subroutine MeshTopologyGetInfo

!!!
!!! Replace with an IS version
!!!
#undef __FUNCT__
#define __FUNCT__ "FieldCreateVertex"
   Subroutine FieldCreateVertex(F,Fname,MeshTopology,component_size)
      Type(Field)                                  :: F
      Character(len=*)                             :: Fname
      Type(MeshTopology_Type)                      :: MeshTopology
      PetscInt,Dimension(:),Pointer                :: component_size
      Character(len=256)                           :: component_name
      PetscInt                                     :: numCells,numVertices
      PetscInt                                     :: i,j,ierr
      
      Call DMMeshGetStratumSize(MeshTopology%mesh,"height",0,numCells,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumSize(MeshTopology%mesh,"depth",0,numVertices,ierr);CHKERRQ(ierr)

      F%num_components = Size(component_size)

      !!! Create the Main Section
      Call DMMeshGetVertexSectionReal(MeshTopology%Mesh,Fname,sum(component_size),F%Sec,ierr);CHKERRQ(ierr)
      
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
      Call DMMeshCreateGlobalScatter(MeshTopology%mesh,F%Sec,F%Scatter,ierr);CHKERRQ(ierr)
      Call DMMeshCreateVector(MeshTopology%mesh,F%Sec,F%Vec,ierr);CHKERRQ(ierr)
      Call SectionRealCreateLocalVector(F%Sec,F%LocalVec,ierr);CHKERRQ(ierr)
   End Subroutine FieldCreateVertex
   
#undef __FUNCT__
#define __FUNCT__ "FieldDestroy"
   Subroutine FieldDestroy(F)
      Type(Field)                                  :: F
      PetscInt                                     :: i,ierr
      
      Call SectionRealDestroy(F%Sec,ierr);CHKERRQ(ierr)
      DeAllocate(F%Component_size)
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
   Subroutine FlagCreateVertex(F,Fname,MeshTopology,component_size)
      !!! creates a Flag derived type with sum(component_size) values at each vertex
      !!! add a way to scatter it into size(component_size) component sections
      !!! I am not sure this actualy works and does anything meaningful if component_size(i) /= 1...
      Type(Flag)                                   :: F
      Character(len=*)                             :: Fname
      Type(MeshTopology_Type)                      :: MeshTopology
      PetscInt,Dimension(:),Pointer                :: component_size
      Character(len=256)                           :: component_name
      PetscInt                                     :: numCells,numVertices
      PetscInt                                     :: i,j,ierr
      
      Call DMMeshGetStratumSize(MeshTopology%mesh,"height",0,numCells,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumSize(MeshTopology%mesh,"depth",0,numVertices,ierr);CHKERRQ(ierr)

      F%num_components = Size(component_size)

      !!! Create the Main Section
      Call DMMeshGetVertexSectionInt(MeshTopology%Mesh,Fname,sum(component_size),F%Sec,ierr);CHKERRQ(ierr)
      
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
#define __FUNCT__ "SectionIntAddNSProperty"
   Subroutine SectionIntAddNSProperty(Sec,NSProperty,MeshTopology)
      !!!
      !!! initialize a SectionInt by adding at each point of a nodeset the value of the
      !!! NSProperty associated to the nodeset
      Type(SectionInt)                             :: Sec 
      Type(EXO_Property_Type)                      :: NSProperty
      Type(MeshTopology_Type)                      :: MeshTopology
      
      PetscErrorCode                               :: ierr
      PetscInt                                     :: set,vertex
      Type(IS)                                     :: setIS,vertexIS
      PetscInt,Dimension(:),Pointer                :: setID,vertexID
      PetscInt                                     :: numCells
      PetscInt,Dimension(:),Pointer                :: Sec_Ptr
      
      Call DMMeshGetStratumSize(MeshTopology%mesh,"height",0,numCells,ierr);CHKERRQ(ierr)
      Call DMMeshGetLabelIdIS(MeshTopology%mesh,'Vertex Sets',setIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Allocate(Sec_Ptr(1))
         Sec_Ptr = NSProperty%Value(setID(set))
         Call DMMeshGetStratumIS(MeshTopology%mesh,'Vertex Sets',setID(set),vertexIS,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(vertexIS,vertexID,ierr);CHKERRQ(ierr)
         Do vertex = 1,size(vertexID)
            Call SectionIntUpdate(Sec,vertexID(vertex) + numCells - 1,Sec_Ptr,ADD_VALUES,ierr);CHKERRQ(ierr)
         End Do
         DeAllocate(Sec_Ptr)
         Call ISRestoreIndicesF90(vertexIS,vertexID,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
   End Subroutine SectionIntAddNSProperty
   
#undef __FUNCT__
#define __FUNCT__ "MatInsertVertexBoundaryValues"
   Subroutine MatInsertVertexBoundaryValues(M,U,BCFlag,MeshTopology)
      Type(Mat)                                    :: M
      Type(Field)                                  :: U
      Type(Flag)                                   :: BCFlag
      Type(MeshTopology_Type)                      :: MeshTopology
      
      PetscInt                                     :: ierr
      PetscInt,Dimension(:),Pointer                :: BCFlag_Ptr
      PetscReal,Dimension(:,:),Pointer             :: MatElem
      PetscInt                                     :: i,j,num_dof,zero,numCells,numVertices
      !!! As soon as I can get access to the layout data of a SectionReal,I won't need the MeshTopology and the specific Vertex case

      !!!
      !!! DMMeshAssembleMatrix does not work with fibrated sections at this point
      !!! in order to INSERT boundary values,I need to insert a block for ALL dof associated to a given point
      !!! therefre erasing exsting values....
      !!! MatInsertBoundaryValues needs to be called BEFORE building the hessian or stiffness matrix
      !!!
      Call DMMeshGetStratumSize(MeshTopology%mesh,"height",0,numCells,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumSize(MeshTopology%mesh,"depth",0,numVertices,ierr);CHKERRQ(ierr)
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
            Call DMMeshAssembleMatrix(M,MeshTopology%mesh,U%Sec,numCells+i-1,MatElem,INSERT_VALUES,ierr);CHKERRQ(ierr)
         End If
      End Do
      DeAllocate(MatElem)
   End Subroutine MatInsertVertexBoundaryValues

#undef __FUNCT__
#define __FUNCT__ "FieldInsertVertexBoundaryValues"
   Subroutine FieldInsertVertexBoundaryValues(F,FBC,BCFlag,MeshTopology)
      !!! Insert in F the value of FBC at each point where BCFlag is non 0
      Type(Field)                                  :: F,FBC
      Type(Flag)                                   :: BCFlag
      Type(MeshTopology_Type)                      :: MeshTopology
      
      PetscInt                                     :: ierr,zero
      PetscInt,Dimension(:),Pointer                :: BCFlag_Ptr
      PetscReal,Dimension(:),Pointer               :: FBC_Ptr
      PetscInt                                     :: i,j,numCells,numVertices
      !!! As soon as I can get access to the layout data of a SectionReal,I won't need the MeshTopology and the specific Vertex case

      Call DMMeshGetStratumSize(MeshTopology%mesh,"height",0,numCells,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumSize(MeshTopology%mesh,"depth",0,numVertices,ierr);CHKERRQ(ierr)
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
