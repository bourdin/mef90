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
    
   Type(MeshTopology_Type)                      :: MeshTopology
   Type(EXO_Type)                               :: EXO
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
   
   Call MeshTopologyReadEXO(MeshTopology, EXO)
   MeshTopology%Elem_Blk%Elem_Type    = MEF90_P1_Lagrange
   Do iBlk = 1, MeshTopology%Num_Elem_Blks
      Call Init_Elem_Blk_Type(MeshTopology%Elem_Blk(iBlk), MeshTopology%num_dim)
   End Do
   If (verbose) Then
      Call MeshTopologyView(MeshTopology, PetscViewer(PETSC_VIEWER_STDOUT_SELF))
   End If

   !!! Initialize the element   
   Allocate(Elem2DA(MeshTopology%Num_Elems))
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
   Allocate(Coords(MeshTopology%Num_Verts))
   Call MeshInitCoordinates(MeshTopology, Coords)
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
!         Call Destroy_Element(Elem2DA(iE))
         !!! DestroyELement crashes when the element has not been fully allocated
         !!! Remove the deallocation of the ID_DoF later if I indeed don;t use the 
         !!! connectivity table anymore
      End Do
   End Do
   Deallocate(Elem2DA)
   Deallocate(Coords)
   Call MeshTopologyDestroy(MeshTopology)
   Call MEF90_Finalize()

End Program TestSieve
