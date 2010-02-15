Program TestDuplicate

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

   Type(MeshTopology_Type)                      :: MeshTopology
   Type(Mesh)                                   :: Tmp_Mesh
   Type(EXO_Type)                               :: EXO, MyEXO
   Type(Element2D_Scal), Dimension(:), Pointer  :: Elem2DA
   
   PetscTruth                                   :: HasPrefix
   PetscTruth                                   :: verbose
   PetscErrorCode                               :: iErr
   Character(len=256)                           :: filename
   Character(len=256)                           :: prefix
   Type(SectionReal)                            :: U_Sec, V_Sec
   Type(VecScatter)                             :: scatter
   PetscInt                                     :: i
   PetscInt                                     :: dof
   PetscInt, Dimension(:), Pointer              :: iRows
   PetscReal                                    :: one = 1.0_Kr
   Character(len=256)                           :: IO_Buffer
   Type(PetscViewer)                            :: MeshViewer
     
   Call MEF90_Initialize()
   dof=1
   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-dof', dof, HasPrefix, iErr)    
   Call PetscOptionsHasName(PETSC_NULL_CHARACTER, '-verbose', verbose, iErr)    
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p', prefix, HasPrefix, iErr)    
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No input file prefix given\n", iErr)
      Call MEF90_Finalize()
      STOP
   End If

   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'


   If (MEF90_NumProcs == 1) Then
      Call MeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, MeshTopology%mesh, ierr); CHKERRQ(iErr)
   Else
      Call MeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
      Call MeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, MeshTopology%mesh, ierr); CHKERRQ(iErr)
      Call MeshDestroy(Tmp_mesh, ierr); CHKERRQ(iErr)
   End If
   
   Call MeshTopologyReadEXO(MeshTopology, EXO)
   
   Call MeshSetMaxDof(MeshTopology%Mesh, dof, iErr); CHKERRQ(iErr) 

   Call MeshGetVertexSectionReal(MeshTopology%mesh, 'U', dof, U_Sec, iErr); CHKERRQ(iErr)

   Call MeshCreateGlobalScatter(MeshTopology%mesh, U_Sec, scatter, iErr); CHKERRQ(iErr)



   Call PetscPrintf(PETSC_COMM_WORLD, "\n\nSection U_Sec: \n", iErr); CHKERRQ(iErr)
   Call SectionRealView(U_Sec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   Call SectionRealGetSize(U_Sec, i, iErr)
   Write(IO_Buffer, *) 'Size of U_Sec: ', i, '\n'
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, IO_Buffer, iErr); CHKERRQ(iErr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)   
!   Call SectionRealDuplicate(U_Sec, V_Sec, iErr); CHKERRQ(iErr)
!   Call PetscPrintf(PETSC_COMM_WORLD, "\n\nSection V_Sec: \n", iErr); CHKERRQ(iErr)
!   Call SectionRealView(V_Sec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)

   Call PetscViewerBinaryOpen(PETSC_COMM_WORLD, Trim(Prefix)//'.dat', MeshViewer, iErr); CHKERRQ(iErr)
!   Call PetscViewerASCIIPrintf(MeshViewer, "MeshTopology%Mesh: \n", iErr); CHKERRQ(iErr)
   Call MeshView(MeshTopology%Mesh, MeshViewer, iErr); CHKERRQ(iErr)
   

   Call MeshTopologyDestroy(MeshTopology)
   Call SectionRealDestroy(U_Sec, iErr); CHKERRQ(iErr)
   
   
   Call VecScatterDestroy(scatter, iErr); CHKERRQ(iErr)
   
300 Format(A)   
   Call MEF90_Finalize()
End Program TestDuplicate
