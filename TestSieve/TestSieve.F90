Program TestSieve

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscmeshdef.h"

   Use m_MEF90
   Use petsc
   Use petscvec
   Use petscmat
   Use petscmesh
   Use m_TestSieve
   
   Implicit NONE
    
   Type(MeshTopology_Info)                      :: MeshTopology
   Type(EXO_Info)                               :: EXO
   Type(Element2D_Scal), Dimension(:), Pointer  :: Elem2DA
   Type(Vect3D), Dimension(:), Pointer          :: Coords
   PetscReal, Dimension(:,:), Pointer           :: Vertices
   
   PetscReal                                    :: MyObjectiveFunction, ObjectiveFunction
   Type(SectionReal)                            :: U, F, coordSection
   PetscReal, Dimension(:), Pointer             :: values
   Type(Mat)                                    :: K
   Type(Vec)                                    :: V, V_Loc
   PetscTruth                                   :: HasF
   PetscInt                                     :: dof
   PetscLogEvent                                :: integrationEvent
   PetscTruth                                   :: verbose
   PetscErrorCode                               :: iErr
   PetscInt                                     :: iBlk, iELoc, iE, iV, iX
   Character(len=256)                           :: CharBuffer
   Type(VecScatter)                             :: scatter
   
     
   Call MEF90_Initialize()
   dof = 1
   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-dof', dof, HasF, iErr)    
   Call PetscOptionsHasName(PETSC_NULL_CHARACTER, '-verbose', verbose, iErr)    
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-f', EXO%filename, HasF, iErr)    
   If (.NOT. HasF) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No input file given\n", iErr)
      Call MEF90_Finalize()
      STOP
   End If

   call PetscLogEventRegister('ElemInteg', 0, integrationEvent, ierr); CHKERRQ(ierr)

   EXO%Comm = PETSC_COMM_WORLD
   
   Call Read_MeshTopology_Info_EXO(MeshTopology, Coords, Elem2DA, EXO)
   
   MeshTopology%Elem_Blk%Elem_Type    = MEF90_P1_Lagrange
   Do iBlk = 1, MeshTopology%Num_Elem_Blks
      Call Init_Elem_Blk_Info(MeshTopology%Elem_Blk(iBlk), MeshTopology%num_dim)
   End Do
   If (verbose) Then
      Call MeshTopologyView(MeshTopology, PetscViewer(PETSC_VIEWER_STDOUT_SELF))
   End If

   !!! Initialize the element   
   Allocate(Vertices(2,3))
   Allocate(values(6))
   Call MeshGetSectionReal(MeshTopology%mesh, 'coordinates', coordSection, iErr); CHKERRQ(ierr)
   !!! Gets the Section 'coordinate'
   Do iBlk = 1, MeshTopology%Num_Elem_Blks
      Do iELoc = 1, MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         call MeshRestrictClosure(MeshTopology%mesh, coordSection, iE-1, Size(values), values, ierr)
         !!! iE-1 is a point name
         !!! puts all the "values" for element iE (iE-1 in C) of the section coordsection into the vector 
         !!! "values"
         !Vertices(1,:) = Coords(Elem2DA(iE)%ID_DoF(:))%X
         !Vertices(2,:) = Coords(Elem2DA(iE)%ID_DoF(:))%Y
         Do iX = 1, 3
            Vertices(1,iX) = values(2*iX-1)
            Vertices(2,iX) = values(2*iX)
         End Do
         Call Init_Element(Elem2DA(iE), Vertices, 4, MeshTopology%Elem_Blk(iBlk)%Elem_Type)
      End Do
   End Do
   Deallocate(Vertices)
   Deallocate(values)
   
!   Call Show_Elem2D_Scal(Elem2DA)

   Call MeshGetVertexSectionReal(MeshTopology%mesh, dof, U, ierr); CHKERRQ(iErr)
   Call MeshGetVertexSectionReal(MeshTopology%mesh, dof, F, ierr); CHKERRQ(iErr)
   !!! Unlike the coordinate section, it is not stored in the mesh hence doesn't have a name
   !!! It would not be redistributed if we'd repartition the mesh whereas the coordinate section does

   !!! In order to add the section to the mesh (and get it automatically repartitioned)
   
   !!! MeshSetSectionInt(MeshTopology%mesh, 'BCFlag', flag, ierr);
   !!! Or Use MeshGetSectionInt which would not create it
   !!! At this point the section would still be empty
   !!! The calls would be 
   !!!    SetFiberDimension (flag, iV, dof) (loop through all points iV) 
   !!!    (note that iV can be of different types, and dof can be 0


   Call MeshCreateMatrix(MeshTopology%mesh, U, MATMPIAIJ, K, iErr); CHKERRQ(iErr)
   Call MeshCreateVector(MeshTopology%mesh, U, V, iErr); CHKERRQ(iErr)
   Call MatZeroEntries(K, iErr); CHKERRQ(ierr)
   Call MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(ierr)
   Call MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(ierr)
   Call PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INFO_DETAIL, ierr); CHKERRQ(ierr)
   Call MatView(K, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(ierr)

   Call MeshCreateGlobalScatter(MeshTopology%mesh, U, scatter, iErr); CHKERRQ(iErr)
   !!! This Scatter maps from unassembled section storage to global section storage
   !!! To be used in SectionRealToVec

   !Call VecScatterView(scatter, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)

   Call SectionRealSet(U, 1.0_Kr, iErr); CHKERRQ(iErr)
   Call SectionRealSet(F, 1.0_Kr, iErr); CHKERRQ(iErr)

   
   Allocate(values(dof))
   Do iV = 1, MeshTopology%Num_Verts
      values = 1.0+Coords(iV)%Y
      call MeshUpdateClosure(MeshTopology%mesh, U, MeshTopology%Num_Elems+iV-1, values, ierr)
      !!! Internal storage of points in the mesh is cell then vertices then everything else
      !!! So vertex iV is at offset MeshTopology%Num_Elems+iV-1
      !!! This is purely local (does not update "ghost values")
      !!! This sets dof values

      !!! Same thing for a Section defined over elements would be
      !!! call MeshUpdateClosure(MeshTopology%mesh, U, iE, value, ierr)

      !!! The inverse function MeshRestrictClosure
      !!! MeshUpdateClosure(MeshTopology%mesh, 'coordinates', iE, value, ierr)
      !!! will give the coordinates of the closure of iE which can be a cell an edge a vertex etc
   End Do
   Deallocate(values)


   Call SectionRealView(U, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(ierr)
   Call SectionRealView(F, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(ierr)

   Call FormObjectiveFunction(MyObjectiveFunction, MeshTopology, Elem2DA, U, F, integrationEvent)

   Call PetscGlobalSum(MyObjectiveFunction, ObjectiveFunction, PETSC_COMM_WORLD, ierr); CHKERRQ(iErr)
   Write(CharBuffer,*) MEF90_MyRank, ' My Objective Function: ', MyObjectiveFunction, '\n'
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, CharBuffer, ierr); CHKERRQ(iErr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)

   Write(CharBuffer,*) '               Objective Function: ', ObjectiveFunction, '\n'
   Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, ierr); CHKERRQ(iErr)

   Call SectionRealToVec(U, scatter, SCATTER_FORWARD, V, ierr); CHKERRQ(ierr)
   !!! U is the equivalent of a ghosted vector
   !!! V is the equivalent of the global vector
   !!! They are both collective on PETSC_COMM_WORLD
   !!! The ghost update is INSERT by default
   !!! Doing a SCATTER_REVERSE after that would result in a ghost update

   !!! Ghost update can also be done using SectionComplete (equivalent of DALocalToLocal)

   Write(CharBuffer,*) 'V\n'c
   Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, ierr); CHKERRQ(iErr)
   Call VecView(V, PETSC_VIEWER_STDOUT_WORLD, ierr); CHKERRQ(ierr)

   !!! Destroy the element   
   Do iBlk = 1, MeshTopology%Num_Elem_Blks
      Do iELoc = 1, MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         Call Destroy_Element(Elem2DA(iE))
      End Do
   End Do
   Deallocate(Elem2DA)
   Deallocate(Coords)
   Call MeshTopologyDestroy(MeshTopology)
   Call MEF90_Finalize()

 Contains
   Subroutine Read_MeshTopology_Info_EXO(dMeshTopology, Coords, Elem2DA, dEXO)
      Type(MeshTopology_Info)                      :: dMeshTopology
      Type(Vect3D), Dimension(:), Pointer          :: Coords
      Type(Element2D_Scal), Dimension(:), Pointer  :: Elem2DA
      Type(EXO_Info)                               :: dEXO
      PetscInt                                     :: iErr, iBlk, iSet
      Character(len=256)                           :: CharBuffer
      
      PetscReal, Dimension(:,:), Pointer           :: array
      PetscInt, Dimension(:,:), Pointer            :: arrayCon
      PetscInt                                     :: embedDim
      PetscInt                                     :: iElem, numIds, blkId, setId
      PetscInt, Dimension(:), Pointer              :: blkIds
      PetscInt, Dimension(:), Pointer              :: setIds
      Type(Mesh)                                   :: mesh
      Type(PetscViewer)                            :: viewer

      ! Open File
      call MeshCreateExodus(PETSC_COMM_WORLD, dEXO%filename, mesh, ierr)
      !!! reads exo file, stores all information in a Mesh

      call MeshDistribute(mesh, PETSC_NULL_CHARACTER, dMeshTopology%mesh, ierr)
      !!! Partitions using a partitioner (currently PETSC_NULL_CHARACTER) 
      call MeshDestroy(mesh, ierr)

      If (verbose) Then
         call PetscViewerASCIIOpen(PETSC_COMM_WORLD, PETSC_NULL_CHARACTER, viewer, ierr)
         call PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INFO_DETAIL, ierr)
         call MeshView(dMeshTopology%mesh, viewer, ierr)
         call PetscViewerDestroy(viewer, ierr)
      End If

      ! Read Global Geometric Parameters
      call MeshExodusGetInfo(dMeshTopology%mesh, dMeshTopology%Num_Dim, dMeshTopology%Num_Verts, dMeshTopology%Num_Elems, dMeshTopology%Num_Elem_Blks, dMeshTopology%Num_Node_Sets, iErr)
      !!! Extracts sizes from the Mesh oject


      ! Read Elem block information
      CharBuffer = 'CellBlocks'
      Call MeshGetLabelSize(dMeshTopology%mesh, CharBuffer, numIds, ierr); CHKERRQ(ierr)  
      !!! Get the number of labels of type 'CellBlocks' in the mesh
      If (numIds .ne. dMeshTopology%Num_Elem_blks) Then
         SETERRQ(PETSC_ERR_ARG_SIZ, 'Invalid number of element ids', ierr)
      End If
      !!! Compare to the number initialized in MeshTopology
      
      Allocate(dMeshTopology%Elem_blk(dMeshTopology%Num_Elem_blks))
      Allocate(blkIds(numIds))
      call MeshGetLabelIds(dMeshTopology%mesh, CharBuffer, blkIds, ierr); CHKERRQ(ierr)
      If (dMeshTopology%Num_Elem_blks > 0) Then
         Do iBlk = 1, dMeshTopology%Num_Elem_Blks
            blkId = blkIds(iBlk)
            dMeshTopology%Elem_blk(iBlk)%ID = blkId
            call MeshGetStratumSize(dMeshTopology%mesh, CharBuffer, blkId, dMeshTopology%elem_blk(iBlk)%Num_Elems, ierr)
            !!! Get the size of the layer (stratum) 'CellBlock' of Mesh
            Allocate(dMeshTopology%Elem_blk(iBlk)%Elem_ID(dMeshTopology%elem_blk(iBlk)%Num_Elems))
            call MeshGetStratum(dMeshTopology%mesh, CharBuffer, blkId, dMeshTopology%Elem_blk(iBlk)%Elem_ID, ierr)
            !!! Get the layer (stratum) 'CellBlock' of Mesh in C numbering
            dMeshTopology%Elem_blk(iBlk)%Elem_ID = dMeshTopology%Elem_blk(iBlk)%Elem_ID + 1
            !!! Converts to Fortran style indexing
         End Do
      End If
      Deallocate(blkIds)
      
      ! Read Node set information
      CharBuffer = 'VertexSets'
      call MeshGetLabelSize(dMeshTopology%mesh, CharBuffer, numIds, ierr); CHKERRQ(ierr)
      If (numIds .ne. dMeshTopology%Num_node_sets) Then
         SETERRQ(PETSC_ERR_ARG_SIZ, 'Invalid number of node ids', ierr)
      End If
      Allocate(dMeshTopology%Node_Set(dMeshTopology%Num_Node_Sets))
      Allocate(setIds(numIds))
      call MeshGetLabelIds(dMeshTopology%mesh, CharBuffer, setIds, ierr); CHKERRQ(ierr)
      If (dMeshTopology%Num_Node_Sets > 0) Then
         Do iSet = 1, dMeshTopology%Num_node_sets
            setId = setIds(iSet)
            dMeshTopology%Node_Set(iSet)%ID = setId
            call MeshGetStratumSize(dMeshTopology%mesh, CharBuffer, setId, dMeshTopology%Node_Set(iSet)%Num_Nodes, ierr)
            Allocate(dMeshTopology%Node_Set(iSet)%Node_ID(dMeshTopology%Node_Set(iSet)%Num_Nodes))
            call MeshGetStratum(dMeshTopology%mesh, CharBuffer, setId, dMeshTopology%Node_Set(iSet)%Node_ID, ierr)
            dMeshTopology%Node_Set(iSet)%Node_ID = dMeshTopology%Node_Set(iSet)%Node_ID - dMeshTopology%Num_Elems + 1
         End Do
      End If
      Deallocate(setIds)

      ! Read the vertices coordinates
      Allocate(Coords(MeshTopology%Num_Verts))
      call MeshGetCoordinatesF90(dMeshTopology%mesh, array, iErr)
      embedDim = size(array,2)
      Coords%X = array(:,1)
      If (embedDim > 1) Then
         Coords%Y = array(:,2)
      Else
         Coords%Y = 0.0
      EndIf
      If (embedDim > 2) Then
         Coords%Z = array(:,3)
      Else
         Coords%z = 0.0
      EndIf
      call MeshRestoreCoordinatesF90(dMeshTopology%mesh, array, iErr)
 
      ! Read the connectivity table
      Allocate(Elem2DA(MeshTopology%Num_Elems))
      call MeshGetElementsF90(dMeshTopology%mesh, arrayCon, iErr)
      Do iBlk = 1, MeshTopology%Num_Elem_Blks
         Do iE = 1, MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iElem = MeshTopology%Elem_Blk(iBlk)%Elem_ID(iE)
            Allocate(Elem2DA(iElem)%ID_DoF(3))
            Elem2DA(iElem)%ID_DoF = arrayCon(iElem,:)
         End Do
      End Do
      call MeshRestoreElementsF90(dMeshTopology%mesh, arrayCon, iErr)

      dEXO%exoid = 0
    End Subroutine Read_MeshTopology_Info_EXO



End Program TestSieve
