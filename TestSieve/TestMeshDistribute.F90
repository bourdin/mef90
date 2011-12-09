Program TestMeshDistribute


#include "finclude/petscdef.h"

   Use m_MEF90
   Use petsc
   Use petscvec
   Use petscmat
   Use petscdm
   
   Implicit NONE   
   
   Type(MeshTopology_Type)                      :: MeshTopology
   Type(EXO_Type)                               :: EXO, MyEXO
   Type(Vec)                                    :: LocalVec
   PetscReal, Dimension(:), Pointer             :: LocalPtr
         
   PetscBool                                    :: HasPrefix
   PetscErrorCode                               :: iErr
   Character(len=256)                           :: CharBuffer, IOBuffer, filename
   Character(len=256)                           :: prefix
   Type(DM)                                     :: Tmp_Mesh
   PetscInt, Dimension(:), Pointer              :: component_length 
   PetscInt                                     :: MySize
   Character(len=MEF90_MXSTRLEN)                :: stagename
   PetscReal                                    :: val   
   
   Call MEF90_Initialize()
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p', prefix, HasPrefix, iErr)    
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No input file prefix given\n", iErr)
      Call MEF90_Finalize()
      STOP
   End If
   
   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'
   
   
   If (MEF90_NumProcs == 1) Then
      Write(IOBuffer, *) "Test me on multiple procs\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call MEF90_Finalize()
      STOP
   Else
      Write(IOBuffer, *) "Calling MeshCreateExodus\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call DMMeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
      Write(IOBuffer, *) "Calling MeshDistribute\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call DMMeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, MeshTopology%mesh, ierr); CHKERRQ(iErr)
      Call DMDestroy(Tmp_mesh, ierr); CHKERRQ(iErr)
   End If
   
   Call MEF90_Finalize()

End Program TestMeshDistribute