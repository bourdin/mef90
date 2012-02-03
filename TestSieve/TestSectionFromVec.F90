Program TestSectionFromVec

#include "finclude/petscdef.h"

   Use m_MEF90
   Use petsc

   Implicit NONE   

   Type (MeshTopology_Type)                     :: MeshTopology
   Type (EXO_Type)                              :: EXO, MyEXO
   Type (Element2D_Scal), Dimension(:), Pointer :: Elem2DA
   Type (Vect3D), Dimension(:), Pointer         :: Coords
   
   Type(DM)                                     :: Tmp_Mesh
   PetscBool                                    :: HasPrefix
   PetscBool                                    :: verbose
   PetscErrorCode                               :: iErr
   PetscInt                                     :: i
   Character(len=256)                           :: CharBuffer
   Character(len=256)                           :: prefix
   Type(PetscViewer)                            :: viewer, myviewer
   Type(SectionReal)                            :: Sscal
   Type(Vec)                                    :: Vscal, VVec
   Type(VecScatter)                             :: ScatterScal
   Type(Mat)                                    :: KScal, KVect
   PetscInt                                     :: dof = 2
   PetscReal, Dimension(:), Pointer             :: ValScal
   PetscInt                                     :: nloc, nglob
   PetscReal                                    :: one = 1.0
   Integer                                      :: n
     
   Call MEF90_Initialize()
   Call PetscOptionsHasName(PETSC_NULL_CHARACTER, '-verbose', verbose, iErr); CHKERRQ(iErr)
   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-dof', dof,  HasPrefix, iErr); CHKERRQ(iErr) 
   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-n', n,  HasPrefix, iErr); CHKERRQ(iErr) 
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p', prefix, HasPrefix, iErr); CHKERRQ(iErr)    
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No input file prefix given\n", iErr); CHKERRQ(iErr)
      Call MEF90_Finalize()
      STOP
   End If

   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'


   Call DMMeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
   Call DMMeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, MeshTopology%mesh, ierr); CHKERRQ(iErr)
   Call DMDestroy(Tmp_mesh, iErr); CHKERRQ(iErr)
   Call MeshTopologyGetInfo(MeshTopology, PETSC_COMM_WORLD)

   Call PetscMemorySetGetMaximumUsage(iErr); CHKERRQ(iErr)
   Call DMMeshGetVertexSectionReal(MeshTopology%mesh, 'scal', 1, SScal, ierr); CHKERRQ(iErr)
   Call DMMeshCreateGlobalScatter(MeshTopology%mesh, SScal, ScatterScal, iErr); CHKERRQ(iErr)
   Call DMMeshCreateGlobalVector(MeshTopology%mesh, SScal, VScal, iErr); CHKERRQ(iErr)

   Call DMMeshSetMaxDof(MeshTopology%Mesh, dof+1, iErr); CHKERRQ(iErr) 


   Do i = 1, N
      If (mod(i, 1000) == 0 ) Then 
         Write(CharBuffer, *) i, "\n"
         Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr); CHKERRQ(iErr)
         Call PetscMemoryShowUsage(PETSC_VIEWER_STDOUT_WORLD, " ", iErr); CHKERRQ(iErr)
      End If
      
      Call SectionRealDestroy(SScal, iErr); CHKERRQ(iErr)
      Call VecDestroy(VScal, iErr); CHKERRQ(iErr)
   
      Call DMMeshGetVertexSectionReal(MeshTopology%mesh, 'scal', 1, SScal, ierr); CHKERRQ(iErr)
      Call DMMeshCreateGlobalVector(MeshTopology%mesh, SScal, VScal, iErr); CHKERRQ(iErr)
      
      
      Call SectionRealSet(SScal, one, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(SScal, ScatterScal, SCATTER_FORWARD,  VScal, ierr); CHKERRQ(ierr)

      Call VecSet(VScal, one, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(SScal, ScatterScal, SCATTER_REVERSE,  VScal, ierr); CHKERRQ(ierr)
   End Do   
   
   Call SectionRealDestroy(SScal, iErr); CHKERRQ(iErr)
   Call VecDestroy(VScal, iErr); CHKERRQ(iErr)
   Call DMDestroy(MeshTopology%Mesh, iErr); CHKERRQ(iErr)
   Call MEF90_Finalize()
End Program TestSectionFromVec
