Module m_MEF_Sieve
#include "finclude/petscdef.h"

   Use m_MEF_LinAlg
   Use m_MEF_Types
   Use m_MEF_Utils
   Use petsc
         
   IMPLICIT NONE
   Private
   
   Public :: MeshTopologyReadEXO
   Public :: MeshInitCoordinates
   Public :: FieldCreateVertex
   Public :: FlagCreateVertex
   Public :: FieldDestroy
   Public :: FlagDestroy
   Public :: SectionIntAddNSProperty
   Public :: MatInsertVertexBoundaryValues
   Public :: FieldInsertVertexBoundaryValues   
   
   Interface MeshInitCoordinates
      Module Procedure MeshInitCoordinatesVect2D, MeshInitCoordinatesVect3D
   End Interface MeshInitCoordinates
   
Contains
#undef __FUNCT__
#define __FUNCT__ "MeshTopologyReadEXO"
   Subroutine MeshTopologyReadEXO(dMeshTopology, dEXO)
      !!! Remove the element and coordinate stuff and move in separate functions
      Type (MeshTopology_Type)                     :: dMeshTopology
      Type (EXO_Type)                              :: dEXO
      PetscInt                                     :: iErr, iBlk, iSet
      Character(len=256)                           :: CharBuffer
      
      PetscReal, Dimension(:,:), Pointer           :: array
      PetscInt, Dimension(:,:), Pointer            :: arrayCon
      PetscInt                                     :: embedDim
      PetscInt                                     :: iE, iElem, numIds, blkId, setId
      PetscInt, Dimension(:), Pointer              :: setIds
      PetscInt, Dimension(:), Pointer              :: blkIds
      PetscInt, Dimension(:), Pointer              :: Tmp_ID, Tmp_GlobalID


      !!! Zero out all informations
      dMeshTopology%num_dim = 0
      dMeshTopology%num_verts = 0
      dMeshTopology%num_elems = 0
      dMeshTopology%num_elem_blks = 0
      dMeshTopology%num_elem_blks_global = 0
      dMeshTopology%num_side_sets = 0
      dMeshTopology%num_side_sets_global = 0
      dMeshTopology%num_node_sets = 0
      dMeshTopology%num_node_sets_global = 0
      
      ! Read Global Geometric Parameters
      !!! Extracts sizes from the Mesh oject
      Call DMMeshExodusGetInfo(dMeshTopology%mesh, dMeshTopology%Num_Dim, dMeshTopology%Num_Verts, dMeshTopology%Num_Elems, dMeshTopology%Num_Elem_Blks, dMeshTopology%Num_Node_Sets, iErr); CHKERRQ(iErr)

      ! Read Elem block information
      CharBuffer = 'CellBlocks'
      Call DMMeshGetLabelSize(dMeshTopology%mesh, CharBuffer, numIds, ierr); CHKERRQ(ierr)  
      !!! Get the number of labels of type 'CellBlocks' in the mesh
      If (numIds .ne. dMeshTopology%Num_Elem_blks) Then
         SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, 'Invalid number of element ids', ierr)
      End If
      !!! Compare to the number initialized in MeshTopology
      
      Allocate(blkIds(numIds))
      Call DMMeshGetLabelIds(dMeshTopology%mesh, CharBuffer, blkIds, ierr); CHKERRQ(ierr)
      Allocate(dMeshTopology%Elem_blk(dMeshTopology%Num_Elem_blks))

      If (dMeshTopology%Num_Elem_blks > 0) Then
         Do iBlk = 1, dMeshTopology%Num_Elem_Blks
            blkId = blkIds(iBlk)
            dMeshTopology%Elem_blk(iBlk)%ID = blkId
            Call DMMeshGetStratumSize(dMeshTopology%mesh, CharBuffer, blkId, dMeshTopology%elem_blk(iBlk)%Num_Elems, ierr); CHKERRQ(iErr)
            !!! Get the size of the layer (stratum) 'CellBlock' of Mesh
            Allocate(dMeshTopology%Elem_blk(iBlk)%Elem_ID(dMeshTopology%elem_blk(iBlk)%Num_Elems))
            Call DMMeshGetStratum(dMeshTopology%mesh, CharBuffer, blkId, dMeshTopology%Elem_blk(iBlk)%Elem_ID, ierr); CHKERRQ(iErr)
            !!! Get the layer (stratum) 'CellBlock' of Mesh in C numbering
            dMeshTopology%Elem_blk(iBlk)%Elem_ID = dMeshTopology%Elem_blk(iBlk)%Elem_ID + 1
            !!! Converts to Fortran style indexing
         End Do
      End If
      Deallocate(blkIds)
      
      ! Get the overall number of element blocks
      Allocate(Tmp_ID(dMeshTopology%num_elem_blks))
      Tmp_ID = dMeshTopology%elem_blk(:)%ID
      Call Uniq(dEXO%Comm, Tmp_ID, Tmp_GlobalID)            
      dMeshTopology%Num_elem_blks_global = Size(Tmp_GlobalID)
      DeAllocate(Tmp_ID)
      DeAllocate(Tmp_GlobalID)
      
      ! Read Node set information
      CharBuffer = 'VertexSets'
      Call DMMeshGetLabelSize(dMeshTopology%mesh, CharBuffer, numIds, ierr); CHKERRQ(ierr)
      If (numIds .ne. dMeshTopology%Num_node_sets) Then
         SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, 'Invalid number of node set ids', ierr)
      End If
      Allocate(dMeshTopology%Node_Set(dMeshTopology%Num_Node_Sets))
      Allocate(setIds(numIds))
      Call DMMeshGetLabelIds(dMeshTopology%mesh, CharBuffer, setIds, ierr); CHKERRQ(ierr)
      If (dMeshTopology%Num_Node_Sets > 0) Then
         Do iSet = 1, dMeshTopology%Num_node_sets
            setId = setIds(iSet)
            dMeshTopology%Node_Set(iSet)%ID = setId
            Call DMMeshGetStratumSize(dMeshTopology%mesh, CharBuffer, setId, dMeshTopology%Node_Set(iSet)%Num_Nodes, ierr); CHKERRQ(iErr)
            Allocate(dMeshTopology%Node_Set(iSet)%Node_ID(dMeshTopology%Node_Set(iSet)%Num_Nodes))
            Call DMMeshGetStratum(dMeshTopology%mesh, CharBuffer, setId, dMeshTopology%Node_Set(iSet)%Node_ID, ierr); CHKERRQ(iErr)
            dMeshTopology%Node_Set(iSet)%Node_ID = dMeshTopology%Node_Set(iSet)%Node_ID - dMeshTopology%Num_Elems + 1
         End Do
      End If
      Deallocate(setIds)

      ! Get the overall number of Node Sets
      Allocate(Tmp_ID(dMeshTopology%num_node_sets))
      Tmp_ID = dMeshTopology%node_set(:)%ID
      Call Uniq(dEXO%Comm, Tmp_ID, Tmp_GlobalID)            
      dMeshTopology%Num_node_sets_global = Size(Tmp_GlobalID)
      DeAllocate(Tmp_ID)
      DeAllocate(Tmp_GlobalID)
   End Subroutine MeshTopologyReadEXO
    
#undef __FUNCT__
#define __FUNCT__ "MeshInitCoordinatesVect2D"
   Subroutine MeshInitCoordinatesVect2D(dMeshTopology, dCoords)
      Type(MeshTopology_Type)                      :: dMeshTopology
      Type(Vect2D), Dimension(:), Pointer          :: dCoords
      
      PetscReal, Dimension(:,:), Pointer           :: array
      PetscInt                                     :: iErr

      Call DMMeshGetCoordinatesF90(dMeshTopology%mesh, array, iErr); CHKERRQ(iErr)
      dCoords%X = array(:,1)
      dCoords%Y = array(:,2)
      Call DMMeshRestoreCoordinatesF90(dMeshTopology%mesh, array, iErr); CHKERRQ(iErr)
   End Subroutine MeshInitCoordinatesVect2D 

#undef __FUNCT__
#define __FUNCT__ "MeshInitCoordinatesVect3D"
   Subroutine MeshInitCoordinatesVect3D(dMeshTopology, dCoords)
      Type(MeshTopology_Type)                      :: dMeshTopology
      Type(Vect3D), Dimension(:), Pointer          :: dCoords

      PetscReal, Dimension(:,:), Pointer           :: array
      PetscInt                                     :: iErr

      Call DMMeshGetCoordinatesF90(dMeshTopology%mesh, array, iErr); CHKERRQ(iErr)
      dCoords%X = array(:,1)
      dCoords%Y = array(:,2)
      dCoords%Z = array(:,3)
      Call DMMeshRestoreCoordinatesF90(dMeshTopology%mesh, array, iErr); CHKERRQ(iErr)
   End Subroutine MeshInitCoordinatesVect3D 


#undef __FUNCT__
#define __FUNCT__ "FieldCreateVertex"
   Subroutine FieldCreateVertex(F, Fname, MeshTopology, component_size)
      Type(Field)                                  :: F
      Character(len=*)                             :: Fname
      Type(MeshTopology_Type)                      :: MeshTopology
      PetscInt, Dimension(:), Pointer              :: component_size
      Character(len=256)                           :: component_name

      PetscInt                                     :: i, j, iErr
      
      F%num_components = Size(component_size)

      !!! Create the Main Section
      Call DMMeshGetVertexSectionReal(MeshTopology%Mesh, Fname, sum(component_size), F%Sec, iErr); CHKERRQ(iErr)
      
      !!! Add space for each of the individual component section
      Allocate(F%Component_size(F%num_components))
      Do i = 1, F%num_components
         Call SectionRealAddSpace(F%Sec, iErr); CHKERRQ(iErr)
         F%component_size(i) = component_size(i)
      End Do 
      Call SectionRealAllocate(F%Sec, iErr); CHKERRQ(iErr)
      
      !!! Set the fibration size
      Do i = 1, MeshTopology%num_verts
         Do j = 1, F%num_components
            Call SectionRealSetFiberDimensionField(F%Sec, i+MeshTopology%Num_Elems-1, component_size(j), j-1, iErr); CHKERRQ(iErr)
         End Do
      End Do 
      
      !!! Create the individual component sections
      Allocate(F%Component_Sec(F%num_Components))
      Do i = 1, F%num_Components      
         Call SectionRealGetFibration(F%Sec, i-1, F%Component_sec(i), iErr); CHKERRQ(iErr)
         Write(component_name, "(A, '.', I3.3)") Trim(Fname), i
         Call PetscObjectSetName(F%Component_Sec(i), component_name, iErr); CHKERRQ(iErr)
      End Do
      !!! Create the Scatter and global and local Vec
!      F%Has_Vec  = .TRUE.
      Call DMMeshCreateGlobalScatter(MeshTopology%mesh, F%Sec, F%Scatter, iErr); CHKERRQ(iErr)
      Call DMMeshCreateVector(MeshTopology%mesh, F%Sec, F%Vec, iErr); CHKERRQ(iErr)
      Call SectionRealCreateLocalVector(F%Sec, F%LocalVec, iErr); CHKERRQ(iErr)
   End Subroutine FieldCreateVertex
   
#undef __FUNCT__
#define __FUNCT__ "FieldDestroy"
   Subroutine FieldDestroy(F)
      Type(Field)                                  :: F
      PetscInt                                     :: i, iErr
      
      Call SectionRealDestroy(F%Sec, iErr); CHKERRQ(iErr)
      DeAllocate(F%Component_size)
      Do i = 1, F%Num_Components   
         Call SectionRealDestroy(F%Component_Sec(i), iErr); CHKERRQ(iErr)   
      End Do
      DeAllocate(F%Component_Sec)
      
      Call VecDestroy(F%Vec, iErr); CHKERRQ(iErr)
      Call VecDestroy(F%LocalVec, iErr); CHKERRQ(iErr)
      Call VecScatterDestroy(F%Scatter, iErr); CHKERRQ(iErr)
   End Subroutine FieldDestroy

#undef __FUNCT__
#define __FUNCT__ "FlagCreateVertex"
   Subroutine FlagCreateVertex(F, Fname, MeshTopology, component_size)
      !!! creates a Flag derived type with sum(component_size) values at each vertex
      !!! add a way to scatter it into size(component_size) component sections
      !!! I am not sure this actualy works and does anything meaningful if component_size(i) /= 1...
      Type(Flag)                                   :: F
      Character(len=*)                             :: Fname
      Type(MeshTopology_Type)                      :: MeshTopology
      PetscInt, Dimension(:), Pointer              :: component_size
      Character(len=256)                           :: component_name

      PetscInt                                     :: i, j, iErr
      
      F%num_components = Size(component_size)

      !!! Create the Main Section
      Call DMMeshGetVertexSectionInt(MeshTopology%Mesh, Fname, sum(component_size), F%Sec, iErr); CHKERRQ(iErr)
      
      !!! Add space for each of the individual component section
      Allocate(F%Component_size(F%num_components))
      Do i = 1, F%num_components
         Call SectionIntAddSpace(F%Sec, iErr); CHKERRQ(iErr)
         F%component_size(i) = component_size(i)
      End Do 
      Call SectionIntAllocate(F%Sec, iErr); CHKERRQ(iErr)
      
      !!! Set the fibration size
      Do i = 1, MeshTopology%num_verts
         Do j = 1, F%num_components
            Call SectionIntSetFiberDimensionField(F%Sec, i+MeshTopology%Num_Elems-1, component_size(j), j-1, iErr); CHKERRQ(iErr)
         End Do
      End Do 
      
      !!! Create the individual component sections
      Allocate(F%Component_Sec(F%num_Components))
      Do i = 1, F%num_Components      
         Call SectionIntGetFibration(F%Sec, i-1, F%Component_sec(i), iErr); CHKERRQ(iErr)
         Write(component_name, "(A, '.', I3.3)") Trim(Fname), i
         Call PetscObjectSetName(F%Component_Sec(i), component_name, iErr); CHKERRQ(iErr)
      End Do
   End Subroutine FlagCreateVertex

#undef __FUNCT__
#define __FUNCT__ "FlagDestroy"
   Subroutine FlagDestroy(F)
      Type(Flag)                                   :: F
      PetscInt                                     :: i, iErr
      
      Call SectionIntDestroy(F%Sec, iErr); CHKERRQ(iErr)
      DeAllocate(F%Component_size)
      Do i = 1, F%Num_Components   
         Call SectionIntDestroy(F%Component_Sec(i), iErr); CHKERRQ(iErr)   
      End Do
      DeAllocate(F%Component_Sec)
   End Subroutine FlagDestroy

#undef __FUNCT__
#define __FUNCT__ "SectionIntAddNSProperty"
   Subroutine SectionIntAddNSProperty(Sec, NSProperty, MeshTopology)
      !!!
      !!! initialize a SectionInt by adding at each point of a nodeset the value of the
      !!! NSProperty associated to the nodeset
      Type(SectionInt)                             :: Sec 
      Type(EXO_Property_Type)                      :: NSProperty
      Type(MeshTopology_Type)                      :: MeshTopology
      
      PetscInt                                     :: iErr, i, j
      PetscInt, Dimension(:), Pointer              :: Sec_Ptr
      
      Do i = 1, MeshTopology%Num_Node_Sets
         Allocate(Sec_Ptr(1))
         Sec_Ptr = NSProperty%Value( MeshTopology%Node_Set(i)%ID )
         Do j = 1, MeshTopology%Node_Set(i)%Num_Nodes
            Call SectionIntUpdate(Sec, MeshTopology%Node_Set(i)%Node_ID(j) + MeshTopology%Num_Elems-1, Sec_Ptr, ADD_VALUES, iErr); CHKERRQ(iErr)
         End Do
         DeAllocate(Sec_Ptr)
      End Do
   End Subroutine SectionIntAddNSProperty
   
#undef __FUNCT__
#define __FUNCT__ "MatInsertVertexBoundaryValues"
   Subroutine MatInsertVertexBoundaryValues(M, U, BCFlag, MeshTopology)
      Type(Mat)                                    :: M
      Type(Field)                                  :: U
      Type(Flag)                                   :: BCFlag
      Type(MeshTopology_Type)                      :: MeshTopology
      
      PetscInt                                     :: iErr
      PetscInt, Dimension(:), Pointer              :: BCFlag_Ptr
      PetscReal, Dimension(:,:), Pointer           :: MatElem
      PetscInt                                     :: i, j, num_dof, zero
      !!! As soon as I can get access to the layout data of a SectionReal, I won't need the MeshTopology and the specific Vertex case

      !!!
      !!! DMMeshAssembleMatrix does not work with fibrated sections at this point
      !!! in order to INSERT boundary values, I need to insert a block for ALL dof associated to a given point
      !!! therefre erasing exsting values....
      !!! MatInsertBoundaryValues needs to be called BEFORE building the hessian or stiffness matrix
      !!!
      zero = 0
      num_dof = sum(BCFlag%component_size)
      Allocate(MatElem(num_dof, num_dof))  
      Do i = 1, MeshTopology%num_verts
         Call SectionIntRestrict(BCFlag%Sec, MeshTopology%Num_Elems+i-1, BCFlag_Ptr, iErr); CHKERRQ(iErr)
         If (Sum(BCFlag_Ptr) /= 0) Then
            MatElem = 0.0_Kr
            Do j = 1, num_dof
               If (BCFlag_Ptr(j) /= zero) Then
                  MatElem(j,j) = 1.0_Kr
               End If  
            End Do
            Call DMMeshAssembleMatrix(M, MeshTopology%mesh, U%Sec, MeshTopology%Num_Elems+i-1, MatElem, INSERT_VALUES, iErr); CHKERRQ(iErr)
         End If
      End Do
      DeAllocate(MatElem)
   End Subroutine MatInsertVertexBoundaryValues

#undef __FUNCT__
#define __FUNCT__ "FieldInsertVertexBoundaryValues"
   Subroutine FieldInsertVertexBoundaryValues(F, FBC, BCFlag, MeshTopology)
      !!! Insert in F the value of FBC at each point where BCFlag is non 0
      Type(Field)                                  :: F, FBC
      Type(Flag)                                   :: BCFlag
      Type(MeshTopology_Type)                      :: MeshTopology
      
      PetscInt                                     :: iErr, zero
      PetscInt, Dimension(:), Pointer              :: BCFlag_Ptr
      PetscReal, Dimension(:), Pointer             :: FBC_Ptr
      PetscInt                                     :: i, j
      !!! As soon as I can get access to the layout data of a SectionReal, I won't need the MeshTopology and the specific Vertex case
      zero = 0
      Do j = 1, BCFlag%num_components 
         If (BCFlag%Component_size(j) /= 1 ) Then
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, 'FieldInsertVertexBoundaryValues requires scalar components', ierr)
         End If
         Do i = 1, MeshTopology%Num_Verts
            Call SectionIntRestrict(BCFlag%Component_Sec(j), MeshTopology%Num_Elems+i-1, BCFlag_Ptr, iErr); CHKERRQ(iErr)
            If (BCFlag_Ptr(1) /= zero) Then
               Call SectionRealRestrict(FBC%Component_Sec(j), MeshTopology%Num_Elems+i-1, FBC_Ptr, iErr); CHKERRQ(iErr)
               Call SectionRealUpdate  (F%Component_Sec(j),   MeshTopology%Num_Elems+i-1, FBC_Ptr, INSERT_VALUES, iErr); CHKERRQ(iErr)
               Call SectionRealRestore (FBC%Component_Sec(j), MeshTopology%Num_Elems+i-1, FBC_Ptr, iErr); CHKERRQ(iErr)
            End If
            Call SectionIntRestore(BCFlag%Component_Sec(j), MeshTopology%Num_Elems+i-1, BCFlag_Ptr, iErr); CHKERRQ(iErr)
         End Do
      End Do
   End Subroutine FieldInsertVertexBoundaryValues
End Module m_MEF_Sieve
