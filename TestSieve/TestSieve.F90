Program TestSieve

   Use m_MEF90
   Use m_TestSieve
   Implicit NONE
    
#include "finclude/petsc.h"
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscviewer.h"
#include "finclude/petscviewer.h90"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscao.h"
#include "finclude/petscmesh.h"
#include "finclude/petscmesh.h90"

   Type (MeshTopology_Info)                     :: MeshTopology
   Type (EXO_Info)                              :: EXO
   Type (Element2D_Scal), Dimension(:), Pointer :: Elem2DA
   Type (Vect3D), Dimension(:), Pointer         :: Coords
   PetscReal, Dimension(:,:), Pointer           :: Vertices
   
   PetscReal                                    :: MyObjectiveFunction, ObjectiveFunction
   SectionReal                                  :: U, F
   PetscReal, Dimension(:), Pointer             :: values
   Mat                                          :: K
   Vec                                          :: V
   PetscTruth                                   :: HasF
   PetscInt                                     :: dof
   PetscLogEvent                                :: integrationEvent
   PetscTruth                                   :: verbose
   PetscErrorCode                               :: iErr
   Integer                                      :: iBlk, iELoc, iE, iV
   Character(len=256)                           :: CharBuffer
   VecScatter                                   :: scatter
   
   Integer, Parameter                           :: exo_cpu_ws = 8
   Integer, Parameter                           :: exo_io_ws = 8
     
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
      Call Show_MeshTopology_Info(MeshTopology)
   End If
   Call Show_MeshTopology_Info(MeshTopology, MEF90_MyRank+100)

   !!! Initialize the element   
   Allocate(Vertices(2,3))
   Do iBlk = 1, MeshTopology%Num_Elem_Blks
      Do iELoc = 1, MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         Vertices(1,:) = Coords(Elem2DA(iE)%ID_DoF(:))%X
         Vertices(2,:) = Coords(Elem2DA(iE)%ID_DoF(:))%Y
         Call Init_Element(Elem2DA(iE), Vertices, 4, MeshTopology%Elem_Blk(iBlk)%Elem_Type)
      End Do
   End Do
   Deallocate(Vertices)

!   Call Show_Elem2D_Scal(Elem2DA)
   
   Call MeshGetVertexSectionReal(MeshTopology%mesh, dof, U, ierr); CHKERRQ(iErr)
   Call MeshGetVertexSectionReal(MeshTopology%mesh, dof, F, ierr); CHKERRQ(iErr)
   Call MeshCreateMatrix(MeshTopology%mesh, U, MATMPIAIJ, K, iErr); CHKERRQ(iErr)
   Call MeshCreateVector(MeshTopology%mesh, U, V, iErr); CHKERRQ(iErr)
   Call MatZeroEntries(K, iErr); CHKERRQ(ierr)
   Call MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(ierr)
   Call MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(ierr)
   Call PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INFO_DETAIL, ierr); CHKERRQ(ierr)
   Call MatView(K, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(ierr)

   Call MeshCreateGlobalScatter(MeshTopology%mesh, U, scatter, iErr); CHKERRQ(iErr)
   !Call VecScatterView(scatter, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)

   Call SectionRealSet(U, 1.0_Kr, iErr); CHKERRQ(iErr)
   Call SectionRealSet(F, 1.0_Kr, iErr); CHKERRQ(iErr)
   Allocate(values(dof))
   Do iV = 1, MeshTopology%Num_Vert
      values(1) = 1.0+Coords(iV)%Y
      call MeshUpdateClosure(MeshTopology%mesh, U, MeshTopology%Num_Elems+iV-1, values, ierr)
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

   !!! Destroy the element   
   Do iBlk = 1, MeshTopology%Num_Elem_Blks
      Do iELoc = 1, MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         Call Destroy_Element(Elem2DA(iE))
      End Do
   End Do
   Deallocate(Elem2DA)
   Deallocate(Coords)
   Call Destroy_MeshTopology_Info(MeshTopology)
   Call MEF90_Finalize()

 Contains
   Subroutine Show_Elem2D_Scal(dElems, Unit)
      Type (Element2D_Scal), DImension(:), Pointer   :: dElems
      Integer, Optional                              :: Unit
      
      Integer                                        :: iE, Nb_Gauss, Nb_DoF, iDoF
      Integer                                        :: IO_Unit
    
      If ( Present(Unit) ) Then
         IO_Unit = Unit
      Else
         IO_Unit = 6 !!! STDOUT
      End If 
 
      
      Write(IO_Unit, 100) Size(dElems)
      
      Do iE = 1, Size(dElems)
         Write(IO_Unit, 101) iE
         Nb_DoF   = size(dElems(iE)%BF,1)
         Nb_Gauss = size(dElems(iE)%BF,2)
         Write(IO_Unit, 102) Nb_DoF
         Write(IO_Unit, 103) Nb_Gauss

         Do iDoF = 1, Nb_DoF
            Write(IO_Unit, 201) iDoF
            Write(IO_Unit, 200, advance = 'no') 'BF'
            Write(IO_Unit, *) dElems(iE)%BF(iDoF, :)
            Write(IO_Unit, 200, advance = 'no') 'Grad_BF%X'
            Write(IO_Unit, *) dElems(iE)%Grad_BF(iDoF, :)%X
            Write(IO_Unit, 200, advance = 'no') 'Grad_BF%Y'
            Write(IO_Unit, *) dElems(iE)%Grad_BF(iDoF, :)%Y
         End Do
      End Do
100 Format('    Number of elements =================== ', I9)
101 Format('*** Element  ', I9)
102 Format('    Nb_DoF   ', I9)
103 Format('    Nb_Gauss ', I9)
200 Format(A)
201 Format('    *** DoF  ', I9)
   End Subroutine Show_Elem2D_Scal
   
   
   Subroutine Read_MeshTopology_Info_EXO(dMeshTopology, Coords, Elem2DA, dEXO)
      Type (MeshTopology_Info)                     :: dMeshTopology
      Type (Vect3D), Dimension(:), Pointer         :: Coords
      Type (Element2D_Scal), Dimension(:), Pointer :: Elem2DA
      Type (EXO_Info)                              :: dEXO
      Integer                                      :: iErr, iBlk, iSet
      Character(len=256)                           :: CharBuffer
      
      PetscReal, Dimension(:,:), Pointer           :: array
      PetscInt, Dimension(:,:), Pointer            :: arrayCon
      Integer                                      :: embedDim
      Integer                                      :: iElem, numIds, blkId, setId
      PetscInt, Dimension(:), Pointer              :: blkIds
      PetscInt, Dimension(:), Pointer              :: setIds
      Mesh                                         :: mesh
      PetscViewer                                  :: viewer

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
      call MeshExodusGetInfo(dMeshTopology%mesh, dMeshTopology%Num_Dim, dMeshTopology%Num_Vert, dMeshTopology%Num_Elems, &
           & dMeshTopology%Num_Elem_Blks, dMeshTopology%Num_Node_Sets, iErr)
      !!! Extracts sizes from the Mesh oject


      ! Read Elem block information
      CharBuffer = 'CellBlocks'
      call MeshGetLabelSize(dMeshTopology%mesh, CharBuffer, numIds, ierr); CHKERRQ(ierr)
      If (numIds .ne. dMeshTopology%Num_Elem_blks) Then
         SETERRQ(PETSC_ERR_ARG_SIZ, 'Invalid number of element ids', ierr)
      End If
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
            !!! COnverts to Fortran style indexing
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
            dMeshTopology%Node_Set(iSet)%Node_ID = dMeshTopology%Node_Set(iSet)%Node_ID - MeshTopology%Num_Elems
         End Do
      End If
      Deallocate(setIds)

      ! Read the vertices coordinates
      Allocate(Coords(MeshTopology%Num_Vert))
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

   Subroutine Show_MeshTopology_Info(dMeshTopology, Unit)
      Type (MeshTopology_Info), Intent(IN)           :: dMeshTopology
      Integer, Optional                              :: Unit
      
      Integer                                        :: IO_Unit, i
    
      If ( Present(Unit) ) Then
         IO_Unit = Unit
      Else
         IO_Unit = 6 !!! STDOUT
      End If  


      Write(IO_Unit, 103) MeshTopology%num_dim
      Write(IO_Unit, 104) MeshTopology%num_Vert
      Write(IO_Unit, 105) MeshTopology%num_elems
      Write(IO_Unit, 106) MeshTopology%Num_Elem_blks
      Write(IO_Unit, 107) MeshTopology%Num_Node_Sets
      Write(IO_Unit, 108) MeshTopology%Num_Side_Sets
      
      Write(IO_Unit, *)
      Write(IO_Unit, 200)
      Write(IO_Unit, 201) MeshTopology%num_elem_blks
      Do i = 1, dMeshTopology%Num_Elem_blks
         Write(IO_Unit, 203) dMeshTopology%Elem_Blk(i)%ID, dMeshTopology%Elem_blk(i)%Num_Elems
         Write(IO_Unit, 204) dMeshTopology%Elem_Blk(i)%ID, dMeshTopology%Elem_blk(i)%Elem_Type
         Write(IO_Unit, 205) dMeshTopology%Elem_Blk(i)%ID, dMeshTopology%Elem_blk(i)%DoF_Location
!         Write(IO_Unit, 206) dMeshTopology%Elem_Blk(i)%ID, dMeshTopology%Elem_blk(i)%Nb_Gauss
         Write(IO_Unit, 207, advance = 'no') dMeshTopology%Elem_Blk(i)%ID
         Write(IO_Unit, *) dMeshTopology%Elem_blk(i)%Elem_ID
      End Do
      
      
      Write(IO_Unit, *)
      Write(IO_Unit, 300)
      Write(IO_Unit, 301) dMeshTopology%num_node_sets
      Do i = 1, dMeshTopology%num_node_sets
         Write(IO_Unit, 302) dMeshTopology%Node_Set(i)%ID, dMeshTopology%Node_Set(i)%Num_Nodes
         Write(IO_Unit, 303, advance = 'no') dMeshTopology%Node_Set(i)%ID
         Write(IO_Unit, *) dMeshTopology%Node_Set(i)%Node_ID
      End Do
      
      
      Write(IO_Unit, *)
      Write(IO_Unit, 400)
      Write(IO_Unit, 401) dMeshTopology%num_side_sets
      
103 Format('    Number of dimensions ============ ', I6)
104 Format('    Number of vertices ============== ', I6)
105 Format('    Number of elements ============== ', I6)
106 Format('    Number of elements blocks ======= ', I6)
107 Format('    Number of node sets ============= ', I6)
108 Format('    Number of side sets ============= ', I6)

200 Format('*** ELEMENT BLOCKS ***')
201 Format('    Number of blocks ================ ', I4)
!202 Format('    Block ', I3, ' Elements type ========= ', A)
203 Format('    Block ', I3, ' Number of elements ==== ', I4)
204 Format('    Block ', I3, ' Element type ========== ', I4)
205 Format('    Block ', I3, ' DoF location ========== ', 4(I4, ' '))
!206 Format('    Block ', I3, ' Nb_Gauss==== ========== ', I4)
207 Format('    Block ', I3, ' IDs: ')

300 Format('*** NODE SETS ***')
301 Format('    Number of sets ================== ', I4)
302 Format('    Set ', I3, ' Number of nodes ========= ', I4)
303 Format('    Set ', I3, ' IDs: ')
!304 Format('    Set ', I3, ' Number of dist. factors = ', I4)

400 Format('*** SIDE SETS ***')
401 Format('    Number of side sets ============= ', I4)
    
   End Subroutine Show_MeshTopology_Info

   Subroutine Destroy_MeshTopology_Info(dMeshTopology)
     Type (MeshTopology_Info) :: dMeshTopology
     PetscInt                 :: iSet, iBlk

     If (dMeshTopology%Num_Node_Sets > 0) Then
        Do iSet = 1, dMeshTopology%Num_node_sets
           Deallocate(dMeshTopology%Node_Set(iSet)%Node_ID)
        End Do
     End If
     Deallocate (dMeshTopology%Node_Set)
     If (dMeshTopology%Num_Elem_blks > 0) Then
        Do iBlk = 1, dMeshTopology%Num_Elem_Blks
           Deallocate(dMeshTopology%Elem_blk(iBlk)%Elem_ID)
        End Do
     End If
     Deallocate(dMeshTopology%Elem_blk)
   End Subroutine Destroy_MeshTopology_Info

End Program TestSieve
