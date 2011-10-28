Program TestSieve3D

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
   Use m_TestSieve3D
   
   Implicit NONE
    
   Type(MeshTopology_Type)                      :: MeshTopology
   Type(Mesh)                                   :: Tmp_Mesh
   Type(EXO_Type)                               :: EXO
   Type(Element3D_Scal), Dimension(:), Pointer  :: Elem
   Type(Vect3D), Dimension(:), Pointer          :: Coords
   PetscReal, Dimension(:,:), Pointer           :: Vertices
   
   PetscReal                                    :: MyObjectiveFunction, ObjectiveFunction
   Type(SectionReal)                            :: U, F, coordSection
   PetscReal, Dimension(:), Pointer             :: values
   Type(Mat)                                    :: K
   Type(Vec)                                    :: V, V_Loc
   PetscBool                                   :: HasF
   PetscInt                                     :: dof
   PetscLogEvent                                :: integrationEvent
   PetscBool                                   :: verbose
   PetscErrorCode                               :: iErr
   PetscInt                                     :: iBlk, iELoc, iE, iV, iX
   Character(len=256)                           :: CharBuffer
   Type(VecScatter)                             :: scatter
   Type(PetscViewer)                            :: MyViewer
   
     
   Call MEF90_Initialize()
   
   Write(CharBuffer, 100) MEF90_MyRank
100 Format('Output-',  I4.4, '.log')
   Call PetscViewerASCIIOpen(PETSC_COMM_SELF, CharBuffer, MyViewer, iErr); CHKERRQ(iErr);   

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
   
   Call MeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
   Call MeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, MeshTopology%mesh, ierr); CHKERRQ(iErr)
   Call MeshDestroy(Tmp_Mesh)
   
   Call MeshTopologyGetInfo(MeshTopology, PETSC_COMM_WORLD)
   MeshTopology%Elem_Blk%Elem_Type    = MEF90_P1_Lagrange
   Do iBlk = 1, MeshTopology%Num_Elem_Blks
      Call Init_Elem_Blk_Type(MeshTopology%Elem_Blk(iBlk), MeshTopology%num_dim)
   End Do
   If (verbose) Then
      Call MeshTopologyView(MeshTopology, PetscViewer(PETSC_VIEWER_STDOUT_SELF))
   End If

   !!! Initialize the element   
   Allocate(Elem(MeshTopology%Num_Elems))
!   Allocate(Vertices(2,3))
   Allocate(Vertices(3,4))
!   Allocate(values(6))
   Allocate(values(12))
   
   Call MeshGetSectionReal(MeshTopology%mesh, 'coordinates', coordSection, iErr); CHKERRQ(ierr)
   !!! Gets the Section 'coordinate'
   Do iBlk = 1, MeshTopology%Num_Elem_Blks
      Do iELoc = 1, MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         call MeshRestrictClosure(MeshTopology%mesh, coordSection, iE-1, Size(values), values, ierr)
         !!! iE-1 is a point name
         !!! puts all the "values" for element iE (iE-1 in C) of the section coordsection into the vector 
         !!! "values"
         !Vertices(1,:) = Coords(Elem(iE)%ID_DoF(:))%X
         !Vertices(2,:) = Coords(Elem(iE)%ID_DoF(:))%Y
         Vertices = Reshape(values, (/MeshTopology%Num_Dim, MeshTopology%Elem_Blk(iBlk)%Num_DoF /) )

         Call Init_Element(Elem(iE), Vertices, 2, MeshTopology%Elem_Blk(iBlk)%Elem_Type)
         Call ElementView(Elem(iE), MyViewer)
      End Do
   End Do
   Deallocate(Vertices)
   Deallocate(values)
   
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

   
   Allocate(values(dof))
   Allocate(Coords(MeshTopology%Num_Verts))
   Call MeshInitCoordinates(MeshTopology, Coords)
   Do iV = 1, MeshTopology%Num_Verts
      values = 1.0+Coords(iV)%Z
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
   Deallocate(Coords)


!   Call SectionRealView(U, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(ierr)
!   Call SectionRealView(F, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(ierr)

   Call FormObjectiveFunction(MyObjectiveFunction, MeshTopology, Elem, U, F, integrationEvent)

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

!   Write(CharBuffer,*) 'V\n'
!   Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, ierr); CHKERRQ(iErr)
!   Call VecView(V, PETSC_VIEWER_STDOUT_WORLD, ierr); CHKERRQ(ierr)

   !!! Destroy the element   
   Do iBlk = 1, MeshTopology%Num_Elem_Blks
      Do iELoc = 1, MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
!         Call Destroy_Element(Elem(iE))
         !!! DestroyELement crashes when the element has not been fully allocated
         !!! Remove the deallocation of the ID_DoF later if I indeed don;t use the 
         !!! connectivity table anymore
      End Do
   End Do
   Deallocate(Elem)
   Call MeshTopologyDestroy(MeshTopology)
   Call MEF90_Finalize()

End Program TestSieve3D
