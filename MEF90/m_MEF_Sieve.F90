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
   Subroutine MeshTopologyGetInfo(dMeshTopology,Comm)
      Type (MeshTopology_Type)                     :: dMeshTopology
      MPI_Comm                                     :: Comm
      PetscInt                                     :: iErr,iBlk,iSet
      
      PetscReal,Dimension(:,:),Pointer             :: array
      PetscInt,Dimension(:,:),Pointer              :: arrayCon
      PetscInt                                     :: embedDim
      PetscInt                                     :: iE,iElem,blkId,setId,i
      Type(IS)                                     :: labels
      PetscInt,Dimension(:),Pointer                :: labels_ptr,set_ptr


      !!! Zero out all informations
      dMeshTopology%num_dim = 0
      dMeshTopology%num_verts = 0
      dMeshTopology%num_elems = 0
      dMeshTopology%num_elem_blks = 0
      dMeshTopology%num_side_sets = 0
      dMeshTopology%num_node_sets = 0
      
      ! Read Global Geometric Parameters
      !!! Extracts sizes from the Mesh oject
      Call DMMeshGetDimension(dMeshTopology%mesh,dMeshTopology%num_dim,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumSize(dMeshTopology%mesh,"height",0,dMeshTopology%Num_Elems,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumSize(dMeshTopology%meshFS,"height",0,dMeshTopology%Num_Faces,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumSize(dMeshTopology%mesh,"depth",0,dMeshTopology%Num_Verts,ierr);CHKERRQ(ierr)
      
      ! Read Elem block information
      !!! Get the number of labels of type 'CellBlocks' in the mesh
      !!! Compare to the number initialized in MeshTopology
      Call DMMeshGetLabelIdIS(dMeshTopology%mesh,'Cell Sets',labels,iErr);CHKERRQ(iErr)
      Call MEF90_ISAllGatherMerge(Comm,labels)
      Call ISGetLocalSize(labels,dMeshTopology%Num_Elem_blks,iErr);CHKERRQ(iErr)
      Call ISGetIndicesF90(labels,labels_ptr,ierr);CHKERRQ(ierr)

      Allocate(dMeshTopology%Elem_blk(dMeshTopology%Num_elem_blks))
      Do iBlk = 1,dMeshTopology%Num_elem_blks
         blkId = labels_ptr(iBlk)
         dMeshTopology%elem_blk(iBlk)%elem_type = -1
         dMeshTopology%elem_blk(iBlk)%dof_location = -1
         dMeshTopology%Elem_blk(iBlk)%ID = blkId
         Call DMMeshGetStratumSize(dMeshTopology%mesh,'Cell Sets',blkId,dMeshTopology%elem_blk(iBlk)%Num_Elems,ierr); CHKERRQ(iErr)
         Call DMMeshGetStratumIS(dMeshTopology%mesh,'Cell Sets',blkId,dMeshTopology%elem_blk(iBlk)%Cell_IS,ierr); CHKERRQ(iErr)
         !!!
         !!! Remove when  getting rid of Elem_ID in Elem_Blk_Type
         !!! vvvvvv
         If (dMeshTopology%elem_blk(iBlk)%Num_Elems > 0) Then
            Allocate(dMeshTopology%Elem_blk(iBlk)%Elem_ID(dMeshTopology%elem_blk(iBlk)%Num_Elems))
            Call ISGetIndicesF90(dMeshTopology%elem_blk(iBlk)%Cell_IS,set_ptr,iErr);CHKERRQ(iErr)
            !!! Get the layer (stratum) 'CellBlock' of Mesh in C numbering
            Do i = 1,dMeshTopology%elem_blk(iBlk)%Num_Elems
               dMeshTopology%Elem_blk(iBlk)%Elem_ID(i) = set_ptr(i) + 1 
            End Do   
            !!! Converts to Fortran style indexing
            Call ISRestoreIndicesF90(dMeshTopology%elem_blk(iBlk)%Cell_IS,set_ptr,iErr);CHKERRQ(iErr)
         End If
         !!! ^^^^^^
         !!! Remove when  getting rid of Elem_ID in Elem_Blk_Type
         !!! 
      End Do
      Call ISRestoreIndicesF90(labels,labels_ptr,iErr);CHKERRQ(iErr)
      
      ! Read Node set information
      Call DMMeshGetLabelIdIS(dMeshTopology%mesh,'Vertex Sets',labels,iErr);CHKERRQ(iErr)
      Call MEF90_ISAllGatherMerge(Comm,labels)
      Call ISGetLocalSize(labels,dMeshTopology%Num_node_sets,iErr);CHKERRQ(iErr)
      Call ISGetIndicesF90(labels,labels_ptr,ierr);CHKERRQ(ierr)

      Allocate(dMeshTopology%node_set(dMeshTopology%Num_node_sets))
      Do iBlk = 1,dMeshTopology%Num_node_sets
         blkId = labels_ptr(iBlk)
         dMeshTopology%node_set(iBlk)%ID = blkId
         Call DMMeshGetStratumSize(dMeshTopology%mesh,'Vertex Sets',blkId,dMeshTopology%node_set(iBlk)%Num_Nodes,ierr); CHKERRQ(iErr)
         Call DMMeshGetStratumIS(dMeshTopology%mesh,'Vertex Sets',blkId,dMeshTopology%node_set(iBlk)%Vertex_IS,ierr); CHKERRQ(iErr)
         !!!
         !!! Remove everything when geting rid of node_IS in Node_Set_Type
         !!! vvvvvv
         If (dMeshTopology%node_set(iBlk)%num_nodes > 0) Then
            Allocate(dMeshTopology%node_set(iBlk)%Node_ID(dMeshTopology%node_set(iBlk)%num_nodes))
         
            Call ISGetIndicesF90(dMeshTopology%node_set(iBlk)%Vertex_IS,set_ptr,iErr);CHKERRQ(iErr)
            Do i = 1,dMeshTopology%node_set(iBlk)%num_nodes
               dMeshTopology%node_set(iBlk)%node_ID(i) = set_ptr(i) + 1 - dMeshTopology%num_elems
            End Do   
            !!! Converts to Fortran style indexing
            Call ISRestoreIndicesF90(dMeshTopology%node_set(iBlk)%Vertex_IS,set_ptr,iErr);CHKERRQ(iErr)
         End If
         !!! ^^^^^^
         !!! Remove everything when geting rid of node_IS in Node_Set_Type
         !!! 
      End Do
      Call ISRestoreIndicesF90(labels,labels_ptr,iErr);CHKERRQ(iErr)
   End Subroutine MeshTopologyGetInfo


#undef __FUNCT__
#define __FUNCT__ "FieldCreateVertex"
   Subroutine FieldCreateVertex(F,Fname,MeshTopology,component_size)
      Type(Field)                                  :: F
      Character(len=*)                             :: Fname
      Type(MeshTopology_Type)                      :: MeshTopology
      PetscInt,Dimension(:),Pointer              :: component_size
      Character(len=256)                           :: component_name

      PetscInt                                     :: i,j,iErr
      
      F%num_components = Size(component_size)

      !!! Create the Main Section
      Call DMMeshGetVertexSectionReal(MeshTopology%Mesh,Fname,sum(component_size),F%Sec,iErr); CHKERRQ(iErr)
      
      !!! Add space for each of the individual component section
      Allocate(F%Component_size(F%num_components))
      Do i = 1,F%num_components
         Call SectionRealAddSpace(F%Sec,iErr); CHKERRQ(iErr)
         F%component_size(i) = component_size(i)
      End Do 
      Call SectionRealAllocate(F%Sec,iErr); CHKERRQ(iErr)
      
      !!! Set the fibration size
      Do i = 1,MeshTopology%num_verts
         Do j = 1,F%num_components
            Call SectionRealSetFiberDimensionField(F%Sec,i+MeshTopology%Num_Elems-1,component_size(j),j-1,iErr); CHKERRQ(iErr)
         End Do
      End Do 
      
      !!! Create the individual component sections
      Allocate(F%Component_Sec(F%num_Components))
      Do i = 1,F%num_Components      
         Call SectionRealGetFibration(F%Sec,i-1,F%Component_sec(i),iErr); CHKERRQ(iErr)
         Write(component_name,"(A,'.',I3.3)") Trim(Fname),i
         Call PetscObjectSetName(F%Component_Sec(i),component_name,iErr); CHKERRQ(iErr)
      End Do
      !!! Create the Scatter and global and local Vec
!      F%Has_Vec  = .TRUE.
      Call DMMeshCreateGlobalScatter(MeshTopology%mesh,F%Sec,F%Scatter,iErr); CHKERRQ(iErr)
      Call DMMeshCreateVector(MeshTopology%mesh,F%Sec,F%Vec,iErr); CHKERRQ(iErr)
      Call SectionRealCreateLocalVector(F%Sec,F%LocalVec,iErr); CHKERRQ(iErr)
   End Subroutine FieldCreateVertex
   
#undef __FUNCT__
#define __FUNCT__ "FieldDestroy"
   Subroutine FieldDestroy(F)
      Type(Field)                                  :: F
      PetscInt                                     :: i,iErr
      
      Call SectionRealDestroy(F%Sec,iErr); CHKERRQ(iErr)
      DeAllocate(F%Component_size)
      Do i = 1,F%Num_Components   
         Call SectionRealDestroy(F%Component_Sec(i),iErr); CHKERRQ(iErr)   
      End Do
      DeAllocate(F%Component_Sec)
      
      Call VecDestroy(F%Vec,iErr); CHKERRQ(iErr)
      Call VecDestroy(F%LocalVec,iErr); CHKERRQ(iErr)
      Call VecScatterDestroy(F%Scatter,iErr); CHKERRQ(iErr)
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
      PetscInt,Dimension(:),Pointer              :: component_size
      Character(len=256)                           :: component_name

      PetscInt                                     :: i,j,iErr
      
      F%num_components = Size(component_size)

      !!! Create the Main Section
      Call DMMeshGetVertexSectionInt(MeshTopology%Mesh,Fname,sum(component_size),F%Sec,iErr); CHKERRQ(iErr)
      
      !!! Add space for each of the individual component section
      Allocate(F%Component_size(F%num_components))
      Do i = 1,F%num_components
         Call SectionIntAddSpace(F%Sec,iErr); CHKERRQ(iErr)
         F%component_size(i) = component_size(i)
      End Do 
      Call SectionIntAllocate(F%Sec,iErr); CHKERRQ(iErr)
      
      !!! Set the fibration size
      Do i = 1,MeshTopology%num_verts
         Do j = 1,F%num_components
            Call SectionIntSetFiberDimensionField(F%Sec,i+MeshTopology%Num_Elems-1,component_size(j),j-1,iErr); CHKERRQ(iErr)
         End Do
      End Do 
      
      !!! Create the individual component sections
      Allocate(F%Component_Sec(F%num_Components))
      Do i = 1,F%num_Components      
         Call SectionIntGetFibration(F%Sec,i-1,F%Component_sec(i),iErr); CHKERRQ(iErr)
         Write(component_name,"(A,'.',I3.3)") Trim(Fname),i
         Call PetscObjectSetName(F%Component_Sec(i),component_name,iErr); CHKERRQ(iErr)
      End Do
   End Subroutine FlagCreateVertex

#undef __FUNCT__
#define __FUNCT__ "FlagDestroy"
   Subroutine FlagDestroy(F)
      Type(Flag)                                   :: F
      PetscInt                                     :: i,iErr
      
      Call SectionIntDestroy(F%Sec,iErr); CHKERRQ(iErr)
      DeAllocate(F%Component_size)
      Do i = 1,F%Num_Components   
         Call SectionIntDestroy(F%Component_Sec(i),iErr); CHKERRQ(iErr)   
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
      
      PetscInt                                     :: iErr,i,j
      PetscInt,Dimension(:),Pointer              :: Sec_Ptr
      
      Do i = 1,MeshTopology%Num_Node_Sets
         Allocate(Sec_Ptr(1))
         Sec_Ptr = NSProperty%Value( MeshTopology%Node_Set(i)%ID )
         Do j = 1,MeshTopology%Node_Set(i)%Num_Nodes
            Call SectionIntUpdate(Sec,MeshTopology%Node_Set(i)%Node_ID(j) + MeshTopology%Num_Elems-1,Sec_Ptr,ADD_VALUES,iErr); CHKERRQ(iErr)
         End Do
         DeAllocate(Sec_Ptr)
      End Do
   End Subroutine SectionIntAddNSProperty
   
#undef __FUNCT__
#define __FUNCT__ "MatInsertVertexBoundaryValues"
   Subroutine MatInsertVertexBoundaryValues(M,U,BCFlag,MeshTopology)
      Type(Mat)                                    :: M
      Type(Field)                                  :: U
      Type(Flag)                                   :: BCFlag
      Type(MeshTopology_Type)                      :: MeshTopology
      
      PetscInt                                     :: iErr
      PetscInt,Dimension(:),Pointer              :: BCFlag_Ptr
      PetscReal,Dimension(:,:),Pointer           :: MatElem
      PetscInt                                     :: i,j,num_dof,zero
      !!! As soon as I can get access to the layout data of a SectionReal,I won't need the MeshTopology and the specific Vertex case

      !!!
      !!! DMMeshAssembleMatrix does not work with fibrated sections at this point
      !!! in order to INSERT boundary values,I need to insert a block for ALL dof associated to a given point
      !!! therefre erasing exsting values....
      !!! MatInsertBoundaryValues needs to be called BEFORE building the hessian or stiffness matrix
      !!!
      zero = 0
      num_dof = sum(BCFlag%component_size)
      Allocate(MatElem(num_dof,num_dof))  
      Do i = 1,MeshTopology%num_verts
         Call SectionIntRestrict(BCFlag%Sec,MeshTopology%Num_Elems+i-1,BCFlag_Ptr,iErr); CHKERRQ(iErr)
         If (Sum(BCFlag_Ptr) /= 0) Then
            MatElem = 0.0_Kr
            Do j = 1,num_dof
               If (BCFlag_Ptr(j) /= zero) Then
                  MatElem(j,j) = 1.0_Kr
               End If  
            End Do
            Call DMMeshAssembleMatrix(M,MeshTopology%mesh,U%Sec,MeshTopology%Num_Elems+i-1,MatElem,INSERT_VALUES,iErr); CHKERRQ(iErr)
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
      
      PetscInt                                     :: iErr,zero
      PetscInt,Dimension(:),Pointer              :: BCFlag_Ptr
      PetscReal,Dimension(:),Pointer             :: FBC_Ptr
      PetscInt                                     :: i,j
      !!! As soon as I can get access to the layout data of a SectionReal,I won't need the MeshTopology and the specific Vertex case
      zero = 0
      Do j = 1,BCFlag%num_components 
         If (BCFlag%Component_size(j) /= 1 ) Then
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,'FieldInsertVertexBoundaryValues requires scalar components',ierr)
         End If
         Do i = 1,MeshTopology%Num_Verts
            Call SectionIntRestrict(BCFlag%Component_Sec(j),MeshTopology%Num_Elems+i-1,BCFlag_Ptr,iErr); CHKERRQ(iErr)
            If (BCFlag_Ptr(1) /= zero) Then
               Call SectionRealRestrict(FBC%Component_Sec(j),MeshTopology%Num_Elems+i-1,FBC_Ptr,iErr); CHKERRQ(iErr)
               Call SectionRealUpdate  (F%Component_Sec(j),  MeshTopology%Num_Elems+i-1,FBC_Ptr,INSERT_VALUES,iErr); CHKERRQ(iErr)
               Call SectionRealRestore (FBC%Component_Sec(j),MeshTopology%Num_Elems+i-1,FBC_Ptr,iErr); CHKERRQ(iErr)
            End If
            Call SectionIntRestore(BCFlag%Component_Sec(j),MeshTopology%Num_Elems+i-1,BCFlag_Ptr,iErr); CHKERRQ(iErr)
         End Do
      End Do
   End Subroutine FieldInsertVertexBoundaryValues
End Module m_MEF_Sieve
