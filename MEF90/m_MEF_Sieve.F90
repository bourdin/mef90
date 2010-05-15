Module m_MEF_Sieve
#include "finclude/petscdef.h"
#include "finclude/petscmeshdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"

   Use m_MEF_LinAlg
   Use m_MEF_Types
   Use m_MEF_Utils

   Use petsc
   Use petscmesh
   Use petscvec
   Use petscmat
         
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
      PetscInt, Dimension(:), Pointer              :: blkIds
      PetscInt, Dimension(:), Pointer              :: setIds
      PetscInt, Dimension(:), Pointer              :: Tmp_ID, Tmp_GlobalID


      
      ! Read Global Geometric Parameters
      Call MeshExodusGetInfo(dMeshTopology%mesh, dMeshTopology%Num_Dim, dMeshTopology%Num_Verts, dMeshTopology%Num_Elems, dMeshTopology%Num_Elem_Blks, dMeshTopology%Num_Node_Sets, iErr); CHKERRQ(iErr)
      !!! Extracts sizes from the Mesh oject

      ! Read Elem block information
      CharBuffer = 'CellBlocks'
      Call MeshGetLabelSize(dMeshTopology%mesh, CharBuffer, numIds, ierr); CHKERRQ(ierr)  
      !!! Get the number of labels of type 'CellBlocks' in the mesh
      If (numIds .ne. dMeshTopology%Num_Elem_blks) Then
         SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, 'Invalid number of element ids', ierr)
      End If
      !!! Compare to the number initialized in MeshTopology
      
      Allocate(dMeshTopology%Elem_blk(dMeshTopology%Num_Elem_blks))
      Allocate(blkIds(numIds))
      Call MeshGetLabelIds(dMeshTopology%mesh, CharBuffer, blkIds, ierr); CHKERRQ(ierr)
      If (dMeshTopology%Num_Elem_blks > 0) Then
         Do iBlk = 1, dMeshTopology%Num_Elem_Blks
            blkId = blkIds(iBlk)
            dMeshTopology%Elem_blk(iBlk)%ID = blkId
            Call MeshGetStratumSize(dMeshTopology%mesh, CharBuffer, blkId, dMeshTopology%elem_blk(iBlk)%Num_Elems, ierr); CHKERRQ(iErr)
            !!! Get the size of the layer (stratum) 'CellBlock' of Mesh
            Allocate(dMeshTopology%Elem_blk(iBlk)%Elem_ID(dMeshTopology%elem_blk(iBlk)%Num_Elems))
            Call MeshGetStratum(dMeshTopology%mesh, CharBuffer, blkId, dMeshTopology%Elem_blk(iBlk)%Elem_ID, ierr); CHKERRQ(iErr)
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
      Call MeshGetLabelSize(dMeshTopology%mesh, CharBuffer, numIds, ierr); CHKERRQ(ierr)
      If (numIds .ne. dMeshTopology%Num_node_sets) Then
         SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, 'Invalid number of node set ids', ierr)
      End If
      Allocate(dMeshTopology%Node_Set(dMeshTopology%Num_Node_Sets))
      Allocate(setIds(numIds))
      Call MeshGetLabelIds(dMeshTopology%mesh, CharBuffer, setIds, ierr); CHKERRQ(ierr)
      If (dMeshTopology%Num_Node_Sets > 0) Then
         Do iSet = 1, dMeshTopology%Num_node_sets
            setId = setIds(iSet)
            dMeshTopology%Node_Set(iSet)%ID = setId
            Call MeshGetStratumSize(dMeshTopology%mesh, CharBuffer, setId, dMeshTopology%Node_Set(iSet)%Num_Nodes, ierr); CHKERRQ(iErr)
            Allocate(dMeshTopology%Node_Set(iSet)%Node_ID(dMeshTopology%Node_Set(iSet)%Num_Nodes))
            Call MeshGetStratum(dMeshTopology%mesh, CharBuffer, setId, dMeshTopology%Node_Set(iSet)%Node_ID, ierr); CHKERRQ(iErr)
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
    
   Subroutine MeshInitCoordinatesVect2D(dMeshTopology, dCoords)
      Type(MeshTopology_Type)                      :: dMeshTopology
      Type(Vect2D), Dimension(:), Pointer          :: dCoords
      
      PetscReal, Dimension(:,:), Pointer           :: array
      PetscInt                                     :: iErr

      Call MeshGetCoordinatesF90(dMeshTopology%mesh, array, iErr); CHKERRQ(iErr)
      dCoords%X = array(:,1)
      dCoords%Y = array(:,2)
      Call MeshRestoreCoordinatesF90(dMeshTopology%mesh, array, iErr); CHKERRQ(iErr)
   End Subroutine MeshInitCoordinatesVect2D 

   Subroutine MeshInitCoordinatesVect3D(dMeshTopology, dCoords)
      Type(MeshTopology_Type)                      :: dMeshTopology
      Type(Vect3D), Dimension(:), Pointer          :: dCoords

      PetscReal, Dimension(:,:), Pointer           :: array
      PetscInt                                     :: iErr

      Call MeshGetCoordinatesF90(dMeshTopology%mesh, array, iErr); CHKERRQ(iErr)
      dCoords%X = array(:,1)
      dCoords%Y = array(:,2)
      dCoords%Z = array(:,3)
      Call MeshRestoreCoordinatesF90(dMeshTopology%mesh, array, iErr); CHKERRQ(iErr)
   End Subroutine MeshInitCoordinatesVect3D 

!!!$   Subroutine MeshInitElementConnectivity2D_Scal(dMeshTopology, dElem)
!!!$      Type(MeshTopology_Type)                      :: dMeshTopology
!!!$      Type(Element2D_Scal), Dimension(:), Pointer  :: dElem
!!!$      
!!!$      PetscInt, Dimension(:,:), Pointer            :: arrayCon
!!!$      PetscInt                                     :: iBlk, iE, iEloc
!!!$      PetscInt                                     :: iErr
!!!$      
!!!$      Call MeshGetElementsF90(dMeshTopology%mesh, arrayCon, iErr); CHKERRQ(iErr)
!!!$      !!! What would happen if the number of dof per  element wasn't constant?
!!!$      Do iBlk = 1, dMeshTopology%Num_Elem_Blks
!!!$         Do iEloc = 1, dMeshTopology%Elem_Blk(iBlk)%Num_Elems
!!!$            iE = dMeshTopology%Elem_Blk(iBlk)%Elem_ID(iEloc)
!!!$            Allocate(dElem(iE)%ID_DoF(dMeshTopology%Elem_Blk(iBlk)%Num_DoF))
!!!$            dElem(iE)%ID_DoF = arrayCon(iE,:)
!!!$         End Do
!!!$      End Do
!!!$      Call MeshRestoreElementsF90(dMeshTopology%mesh, arrayCon, iErr); CHKERRQ(iErr)
!!!$    End Subroutine MeshInitElementConnectivity2D_Scal
!!!$
!!!$   Subroutine MeshInitElementConnectivity2D(dMeshTopology, dElem)
!!!$      Type(MeshTopology_Type)                      :: dMeshTopology
!!!$      Type(Element2D), Dimension(:), Pointer       :: dElem
!!!$      
!!!$      PetscInt, Dimension(:,:), Pointer            :: arrayCon
!!!$      PetscInt                                     :: iBlk, iE, iEloc
!!!$      PetscInt                                     :: iErr
!!!$      
!!!$      Call MeshGetElementsF90(dMeshTopology%mesh, arrayCon, iErr); CHKERRQ(iErr)
!!!$      !!! What would happen if the number of dof per  element wasn't constant?
!!!$      Do iBlk = 1, dMeshTopology%Num_Elem_Blks
!!!$         Do iEloc = 1, dMeshTopology%Elem_Blk(iBlk)%Num_Elems
!!!$            iE = dMeshTopology%Elem_Blk(iBlk)%Elem_ID(iEloc)
!!!$            Allocate(dElem(iE)%ID_DoF(dMeshTopology%Elem_Blk(iBlk)%Num_DoF))
!!!$            dElem(iE)%ID_DoF = arrayCon(iE,:)
!!!$         End Do
!!!$      End Do
!!!$      Call MeshRestoreElementsF90(dMeshTopology%mesh, arrayCon, iErr); CHKERRQ(iErr)
!!!$    End Subroutine MeshInitElementConnectivity2D
!!!$
!!!$   Subroutine MeshInitElementConnectivity2D_Elast(dMeshTopology, dElem)
!!!$      Type(MeshTopology_Type)                      :: dMeshTopology
!!!$      Type(Element2D_Elast), Dimension(:), Pointer :: dElem
!!!$      
!!!$      PetscInt, Dimension(:,:), Pointer            :: arrayCon
!!!$      PetscInt                                     :: iBlk, iE, iEloc
!!!$      PetscInt                                     :: iErr
!!!$      
!!!$      Call MeshGetElementsF90(dMeshTopology%mesh, arrayCon, iErr); CHKERRQ(iErr)
!!!$      !!! What would happen if the number of dof per  element wasn't constant?
!!!$      Do iBlk = 1, dMeshTopology%Num_Elem_Blks
!!!$         Do iEloc = 1, dMeshTopology%Elem_Blk(iBlk)%Num_Elems
!!!$            iE = dMeshTopology%Elem_Blk(iBlk)%Elem_ID(iEloc)
!!!$            Allocate(dElem(iE)%ID_DoF(dMeshTopology%Elem_Blk(iBlk)%Num_DoF))
!!!$            dElem(iE)%ID_DoF = arrayCon(iE,:)
!!!$         End Do
!!!$      End Do
!!!$      Call MeshRestoreElementsF90(dMeshTopology%mesh, arrayCon, iErr); CHKERRQ(iErr)
!!!$    End Subroutine MeshInitElementConnectivity2D_Elast
!!!$
!!!$   Subroutine MeshInitElementConnectivity3D_Scal(dMeshTopology, dElem)
!!!$      Type(MeshTopology_Type)                      :: dMeshTopology
!!!$      Type(Element3D_Scal), Dimension(:), Pointer  :: dElem
!!!$      
!!!$      PetscInt, Dimension(:,:), Pointer            :: arrayCon
!!!$      PetscInt                                     :: iBlk, iE, iEloc
!!!$      PetscInt                                     :: iErr
!!!$      
!!!$      Call MeshGetElementsF90(dMeshTopology%mesh, arrayCon, iErr); CHKERRQ(iErr)
!!!$      !!! What would happen if the number of dof per  element wasn't constant?
!!!$      Do iBlk = 1, dMeshTopology%Num_Elem_Blks
!!!$         Do iEloc = 1, dMeshTopology%Elem_Blk(iBlk)%Num_Elems
!!!$            iE = dMeshTopology%Elem_Blk(iBlk)%Elem_ID(iEloc)
!!!$            Allocate(dElem(iE)%ID_DoF(dMeshTopology%Elem_Blk(iBlk)%Num_DoF))
!!!$            dElem(iE)%ID_DoF = arrayCon(iE,:)
!!!$         End Do
!!!$      End Do
!!!$      Call MeshRestoreElementsF90(dMeshTopology%mesh, arrayCon, iErr); CHKERRQ(iErr)
!!!$    End Subroutine MeshInitElementConnectivity3D_Scal
!!!$
!!!$   Subroutine MeshInitElementConnectivity3D(dMeshTopology, dElem)
!!!$      Type(MeshTopology_Type)                      :: dMeshTopology
!!!$      Type(Element3D), Dimension(:), Pointer       :: dElem
!!!$      
!!!$      PetscInt, Dimension(:,:), Pointer            :: arrayCon
!!!$      PetscInt                                     :: iBlk, iE, iEloc
!!!$      PetscInt                                     :: iErr
!!!$      
!!!$      Call MeshGetElementsF90(dMeshTopology%mesh, arrayCon, iErr); CHKERRQ(iErr)
!!!$      !!! What would happen if the number of dof per  element wasn't constant?
!!!$      Do iBlk = 1, dMeshTopology%Num_Elem_Blks
!!!$         Do iEloc = 1, dMeshTopology%Elem_Blk(iBlk)%Num_Elems
!!!$            iE = dMeshTopology%Elem_Blk(iBlk)%Elem_ID(iEloc)
!!!$            Allocate(dElem(iE)%ID_DoF(dMeshTopology%Elem_Blk(iBlk)%Num_DoF))
!!!$            dElem(iE)%ID_DoF = arrayCon(iE,:)
!!!$         End Do
!!!$      End Do
!!!$      Call MeshRestoreElementsF90(dMeshTopology%mesh, arrayCon, iErr); CHKERRQ(iErr)
!!!$    End Subroutine MeshInitElementConnectivity3D
!!!$
!!!$   Subroutine MeshInitElementConnectivity3D_Elast(dMeshTopology, dElem)
!!!$      Type(MeshTopology_Type)                      :: dMeshTopology
!!!$      Type(Element3D_Elast), Dimension(:), Pointer :: dElem
!!!$      
!!!$      PetscInt, Dimension(:,:), Pointer            :: arrayCon
!!!$      PetscInt                                     :: iBlk, iE, iEloc
!!!$      PetscInt                                     :: iErr
!!!$      
!!!$      Call MeshGetElementsF90(dMeshTopology%mesh, arrayCon, iErr); CHKERRQ(iErr)
!!!$      !!! What would happen if the number of dof per  element wasn't constant?
!!!$      Do iBlk = 1, dMeshTopology%Num_Elem_Blks
!!!$         Do iEloc = 1, dMeshTopology%Elem_Blk(iBlk)%Num_Elems
!!!$            iE = dMeshTopology%Elem_Blk(iBlk)%Elem_ID(iEloc)
!!!$            Allocate(dElem(iE)%ID_DoF(dMeshTopology%Elem_Blk(iBlk)%Num_DoF))
!!!$            dElem(iE)%ID_DoF = arrayCon(iE,:)
!!!$         End Do
!!!$      End Do
!!!$      Call MeshRestoreElementsF90(dMeshTopology%mesh, arrayCon, iErr); CHKERRQ(iErr)
!!!$    End Subroutine MeshInitElementConnectivity3D_Elast

   Subroutine FieldCreateVertex(F, Fname, MeshTopology, component_size)
      Type(Field)                                  :: F
      Character(len=*)                             :: Fname
      Type(MeshTopology_Type)                      :: MeshTopology
      PetscInt, Dimension(:), Pointer              :: component_size
      Character(len=256)                           :: component_name

      PetscInt                                     :: i, j, iErr
      
      F%num_components = Size(component_size)

      !!! Create the Main Section
      Call MeshGetVertexSectionReal(MeshTopology%Mesh, Fname, sum(component_size), F%Sec, iErr); CHKERRQ(iErr)
      
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
      Call MeshCreateGlobalScatter(MeshTopology%mesh, F%Sec, F%Scatter, iErr); CHKERRQ(iErr)
      Call MeshCreateVector(MeshTopology%mesh, F%Sec, F%Vec, iErr); CHKERRQ(iErr)
      Call SectionRealCreateLocalVector(F%Sec, F%LocalVec, iErr); CHKERRQ(iErr)
   End Subroutine FieldCreateVertex
   
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

   Subroutine FlagCreateVertex(F, Fname, MeshTopology, component_size)
      Type(Flag)                                   :: F
      Character(len=*)                             :: Fname
      Type(MeshTopology_Type)                      :: MeshTopology
      PetscInt, Dimension(:), Pointer              :: component_size
      Character(len=256)                           :: component_name

      PetscInt                                     :: i, j, iErr
      
      F%num_components = Size(component_size)

      !!! Create the Main Section
      Call MeshGetVertexSectionInt(MeshTopology%Mesh, Fname, sum(component_size), F%Sec, iErr); CHKERRQ(iErr)
      
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

   Subroutine SectionIntAddNSProperty(Flag, NSProperty, MeshTopology)
      Type(SectionInt)                             :: Flag 
      Type(EXO_Property_Type)                      :: NSProperty
      Type(MeshTopology_Type)                      :: MeshTopology
      
      PetscInt                                     :: iErr, i, j
      PetscInt, Dimension(:), Pointer              :: Flag_Ptr
      
      Do i = 1, MeshTopology%Num_Node_Sets
         Allocate(Flag_Ptr(1))
         Flag_Ptr = NSProperty%Value( MeshTopology%Node_Set(i)%ID )
         Do j = 1, MeshTopology%Node_Set(i)%Num_Nodes
            Call SectionIntUpdate(Flag, MeshTopology%Node_Set(i)%Node_ID(j) + MeshTopology%Num_Elems-1, Flag_Ptr, ADD_VALUES, iErr); CHKERRQ(iErr)
         End Do
         DeAllocate(Flag_Ptr)
      End Do
   End Subroutine SectionIntAddNSProperty
   
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
      !!! assembleMatrix does not work with fibrated sections at this point
      !!! in order to INSERT boundary values, I need to insert a block for ALL dof associated to a given point
      !!! therefre erasing exsting values....
      !!! MatInsertBoundaryValues needs to be called BEFORE building the hessian of stifness matrix
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
            Call assembleMatrix(M, MeshTopology%mesh, U%Sec, MeshTopology%Num_Elems+i-1, MatElem, INSERT_VALUES, iErr); CHKERRQ(iErr)
         End If
      End Do
      DeAllocate(MatElem)
!!!$      Allocate(MatElem(1,1))
!!!$      MatElem = 1.0_Kr
!!!$      Do j = 1, BCFlag%num_components 
!!!$         Write(*,*) 'Doing component', j
!!!$         If (BCFlag%Component_size(j) /= 1 ) Then
!!!$            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, 'MatInsertVertexBoundaryValues requires scalar components', ierr)
!!!$         End If
!!!$         Do i = 1, MeshTopology%Num_Verts
!!!$            Call SectionIntRestrict(BCFlag%Component_Sec(j), MeshTopology%Num_Elems+i-1, BCFlag_Ptr, iErr); CHKERRQ(iErr)
!!!$            If (BCFlag_Ptr(1) /= 0) Then
!!!$               Call assembleMatrix(M, MeshTopology%mesh, U%Component_Sec(j), MeshTopology%Num_Elems+i-1, MatElem, ADD_VALUES, iErr); CHKERRQ(iErr)
!!!$            End If
!!!$            Call SectionIntRestore(BCFlag%Component_Sec(j), MeshTopology%Num_Elems+i-1, BCFlag_Ptr, iErr); CHKERRQ(iErr)
!!!$         End Do
!!!$      End Do
!!!$      DeAllocate(MatElem)
   End Subroutine MatInsertVertexBoundaryValues

   Subroutine FieldInsertVertexBoundaryValues(F, FBC, BCFlag, MeshTopology)
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
