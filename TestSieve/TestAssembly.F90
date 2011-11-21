Program TestAssembly

#include "finclude/petscdef.h"

   Use m_MEF90
   Use petsc

   Implicit NONE   

   Type(MeshTopology_Type)                      :: MeshTopology
   Type(DM)                                     :: Tmp_Mesh
   Type(EXO_Type)                               :: EXO, MyEXO
   Type(Element2D_Scal), Dimension(:), Pointer  :: Elem2DA
   
   PetscBool                                   :: HasPrefix
   PetscBool                                   :: verbose
   PetscErrorCode                               :: iErr
   Character(len=256)                           :: filename
   Character(len=256)                           :: prefix
   Type(SectionReal)                            :: U_Sec
   PetscReal, Dimension(:), Pointer             :: val
   PetscReal, Dimension(:,:), Pointer           :: Kelem
   Type(Mat)                                    :: K
   Type(VecScatter)                             :: scatter
   PetscInt                                     :: i
   PetscInt                                     :: dof
   PetscInt, Dimension(:), Pointer              :: iRows
   PetscReal                                    :: one = 1.0_Kr
     
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
      Call DMMeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, MeshTopology%mesh, ierr); CHKERRQ(iErr)
   Else
      Call DMMeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
      Call DMMeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, MeshTopology%mesh, ierr); CHKERRQ(iErr)
      Call DMDestroy(Tmp_mesh, ierr); CHKERRQ(iErr)
   End If
   
   Call MeshTopologyGetInfo(MeshTopology, PETSC_COMM_WORLD)
   
   Call DMMeshSetMaxDof(MeshTopology%Mesh, dof, iErr); CHKERRQ(iErr) 

   Call DMMeshGetVertexSectionReal(MeshTopology%mesh, 'U', dof, U_Sec, iErr); CHKERRQ(iErr)

   Call DMMeshCreateGlobalScatter(MeshTopology%mesh, U_Sec, scatter, iErr); CHKERRQ(iErr)

   Call DMMeshCreateMatrix(MeshTopology%mesh, U_Sec, MATMPIAIJ, K, iErr); CHKERRQ(iErr)

   Allocate(val(3*dof))
   Allocate(Kelem(3*dof,3*dof))
   val = 1.0_Kr
   Kelem = 1.0_Kr
   Do i = 0, MeshTopology%Num_Elems-1
      Call DMMeshAssembleMatrix(K, MeshTopology%mesh, U_Sec, i, Kelem, ADD_VALUES, iErr); CHKERRQ(iErr)
      Call SectionRealUpdateClosure(U_Sec, MeshTopology%mesh, i, val, ADD_VALUES, iErr); CHKERRQ(iErr)  
   End Do
   DeAllocate(val)
   DeAllocate(Kelem)
   Call MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
   Call MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)

   Call PetscPrintf(PETSC_COMM_WORLD, "Matrix K: \n", iErr); CHKERRQ(iErr)
   Call MatView(K, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   
   !MatSetZeroRowsLocal does not work....
   !Allocate(iRows(1))
   !iRows = 2
   !Call MatZeroRowsLocal(K, 1, iRows, 1, one, iErr); CHKERRQ(iErr)
   !DeAllocate(iRows)
   !Call PetscPrintf(PETSC_COMM_WORLD, "Matrix K: \n", iErr); CHKERRQ(iErr)
   !Call MatView(K, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)

   Call PetscPrintf(PETSC_COMM_WORLD, "\n\nSection U_Sec: \n", iErr); CHKERRQ(iErr)
   Call SectionRealView(U_Sec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)

   Call MeshTopologyDestroy(MeshTopology)
   Call SectionRealDestroy(U_Sec, iErr); CHKERRQ(iErr)
   
   
   Call VecScatterDestroy(scatter, iErr); CHKERRQ(iErr)
   
300 Format(A)   
   Call MEF90_Finalize()
End Program TestAssembly
