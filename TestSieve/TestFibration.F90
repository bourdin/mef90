Program TestFibration


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

   Implicit NONE   

   Type (MeshTopology_Type)                     :: MeshTopology
   Type (EXO_Type)                              :: EXO, MyEXO
   Type(Element2D_Elast), Dimension(:), Pointer :: ElemVect
   Type(Element2D_Scal), Dimension(:), Pointer  :: ElemScal
   Type(Field)                                  :: Field1, Field2
   Type(Flag)                                   :: Flag1, Flag2
   
   PetscTruth                                   :: HasPrefix, flg
   PetscInt                                     :: verbose
   PetscErrorCode                               :: iErr, i, j, iBlk
   Character(len=256)                           :: CharBuffer, IOBuffer, filename
   Character(len=256)                           :: prefix
   Type(PetscViewer)                            :: LogViewer, MyLogViewer
   Type(Mesh)                                   :: Tmp_Mesh
   PetscInt                                     :: num_components, num_dof
   PetscInt, Dimension(:), Pointer              :: component_length 

   Call MEF90_Initialize()
   verbose = 0
   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-verbose', verbose, flg, iErr)    
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p', prefix, HasPrefix, iErr)    
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No input file prefix given\n", iErr)
      Call MEF90_Finalize()
      STOP
   End If

   If (verbose > 1) Then
      Write(filename, 101) Trim(prefix), MEF90_MyRank
      Call PetscViewerASCIIOpen(PETSC_COMM_SELF, filename, MyLogViewer, iErr); CHKERRQ(iErr);   
      Write(IOBuffer, 102) MEF90_MyRank, Trim(filename)
      Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call PetscSynchronizedFlush (PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)

      Write(filename, 103) Trim(prefix)
      Call PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, LogViewer, iErr); CHKERRQ(iErr);   
      Write(IOBuffer, 104) Trim(filename)
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
101 Format(A, '-', I4.4, '.log')
102 Format('Output from processor ', I4.4, ' redirected to file ', A, '\n')
103 Format(A,'.log')
104 Format('Collective output redirected to file ', A, '\n')

   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'


   If (MEF90_NumProcs == 1) Then
      If (verbose > 0) Then
         Write(IOBuffer, *) "Reading the mesh\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      Call MeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, MeshTopology%mesh, ierr); CHKERRQ(iErr)
   Else
      If (verbose > 0) Then
         Write(IOBuffer, *) "Calling MeshCreateExodus\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      Call MeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
      If (verbose > 0) Then
         Write(IOBuffer, *) "Calling MeshDistribute\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      Call MeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, MeshTopology%mesh, ierr); CHKERRQ(iErr)
      Call MeshDestroy(Tmp_mesh, ierr); CHKERRQ(iErr)
   End If
   
   If (verbose > 0) Then
      Write(IOBuffer, *) "Initializing MeshTopology object\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   Call MeshTopologyReadEXO(MeshTopology, EXO)
   
   If (verbose > 0) Then
      Write(IOBuffer, *) "Initializing Element types\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   MeshTopology%Elem_Blk%Elem_Type    = MEF90_P1_Lagrange
   Do iBlk = 1, MeshTopology%Num_Elem_Blks
      Call Init_Elem_Blk_Type(MeshTopology%Elem_Blk(iBlk), MeshTopology%num_dim)
   End Do
   
   MyEXO%comm = PETSC_COMM_SELF
   MyEXO%exoid = EXO%exoid
   Write(MyEXO%filename, 200) trim(prefix), MEF90_MyRank
200 Format(A, '-', I4.4, '.gen')
 

!!! Creating the Sec component of the field and flags
   If (verbose > 0) Then
      Write(IOBuffer, *) "Creating Sections\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If

   num_components=2
   Allocate (component_length(2))
   component_length(1) = 1
   component_length(2) = 2
   num_dof = sum(component_length)

   If (verbose > 0) Then
      Write(IOBuffer, *) "Field1.Sec\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   Call MeshGetSectionReal(MeshTopology%Mesh, 'Field1.Sec', Field1%Sec, iErr); CHKERRQ(iErr)
   Do i = 1, MeshTopology%num_verts
      Call SectionRealSetFiberDimension(Field1%Sec, i+MeshTopology%Num_Elems-1, num_dof, iErr); CHKERRQ(iErr)
   End Do 
   Call SectionRealAllocate(Field1%Sec, iErr); CHKERRQ(iErr)
   Do j = 1, num_dof
   Call SectionRealAddSpace(Field1%Sec, iErr); CHKERRQ(iErr)
   End Do 
   Do i = 1, MeshTopology%num_verts
      Do j = 1, num_components
         Call SectionRealSetFiberDimensionField(Field1%Sec, i+MeshTopology%Num_Elems-1, component_length(j), j, iErr); CHKERRQ(iErr)
      End Do
   End Do 

   Allocate(Field1%Component_Sec(num_Components))
   Do i = 1, Num_Components
      Call SectionRealGetFibration(Field1%Sec, i, Field1%Component_sec(i), iErr); CHKERRQ(iErr)
   End Do

   Call SectionRealSet(Field1%Sec, 1.0_Kr, iErr); CHKERRQ(iErr)

   Write(IOBuffer, *) "Field1.Sec\n"
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call SectionRealView(Field1%Sec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   
   Write(IOBuffer, *) "Field1%Component_Sec(1)\n"
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call SectionRealView(Field1%Component_Sec(1), PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   
   Write(IOBuffer, *) "Field1%Component_Sec(2)\n"
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call SectionRealView(Field1%Component_Sec(2), PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   
!!$   Call FieldDestroy(Field1)
!!$   
!!$   Call MeshTopologyDestroy(MeshTopology)
!!$   If (verbose > 0) Then
!!$      Write(IOBuffer, *) "Cleaning up\n"
!!$      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
!!$   End If
!!$   If (verbose>1) Then
!!$      Call PetscViewerFlush(MyLogViewer, iErr); CHKERRQ(iErr)
!!$      Call PetscViewerDestroy(MyLogViewer, iErr); CHKERRQ(iErr)
!!$      Call PetscViewerFlush(LogViewer, iErr); CHKERRQ(iErr)
!!$      Call PetscViewerDestroy(LogViewer, iErr); CHKERRQ(iErr)
!!$   End If

300 Format(A)   
   Call MEF90_Finalize()

End Program TestFibration
