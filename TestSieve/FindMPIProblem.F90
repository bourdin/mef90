Program FindMPIProblem


#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscdmmeshdef.h"

   Use m_MEF90
   Use petsc
   Use petscvec
   Use petscmat
   Use petscdmmesh

   Implicit NONE   

   Type (MeshTopology_Type)                     :: MeshTopology
   Type (EXO_Type)                              :: EXO, MyEXO
   Type(Field)                                  :: Field1
   Type(Vec)                                    :: LocalVec
   PetscReal, Dimension(:), Pointer             :: LocalPtr
         
   PetscBool                                   :: HasPrefix, flg
   PetscErrorCode                               :: iErr, i, j, iBlk
   Character(len=256)                           :: CharBuffer, IOBuffer, filename
   Character(len=256)                           :: prefix
   Type(Mesh)                                   :: Tmp_Mesh
   PetscInt                                     :: n, iDebug, num_components, num_dof
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
   n=5000
   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-n', n, HasPrefix, iErr)    

   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'


   If (MEF90_NumProcs == 1) Then
      Write(IOBuffer, *) "Reading the mesh\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call MeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, MeshTopology%mesh, ierr); CHKERRQ(iErr)
   Else
      Write(IOBuffer, *) "Calling MeshCreateExodus\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call MeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
      Write(IOBuffer, *) "Calling MeshDistribute\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call MeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, MeshTopology%mesh, ierr); CHKERRQ(iErr)
      Call MeshDestroy(Tmp_mesh, ierr); CHKERRQ(iErr)
   End If
   
   Write(IOBuffer, *) "Initializing MeshTopology object\n"
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call MeshTopologyGetInfo(MeshTopology, PETSC_COMM_WORLD)
   
   Write(IOBuffer, *) "Initializing Element types\n"
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   MeshTopology%Elem_Blk%Elem_Type    = MEF90_P1_Lagrange
   Do iBlk = 1, MeshTopology%Num_Elem_Blks
      Call Init_Elem_Blk_Type(MeshTopology%Elem_Blk(iBlk), MeshTopology%num_dim)
   End Do
   
   MyEXO%comm = PETSC_COMM_SELF
   MyEXO%exoid = EXO%exoid
   Write(MyEXO%filename, 200) trim(prefix), MEF90_MyRank
200 Format(A, '-', I4.4, '.gen')
 

!!! Creating the Sec component of the field and flags
   Write(IOBuffer, *) "Creating Sections\n"
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

   num_components=2
   Allocate (component_length(2))
   component_length(1) = 1
   component_length(2) = 2
   num_dof = sum(component_length)

   Write(IOBuffer, *) "Field1.Sec\n"
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   
   Call FieldCreateVertex(Field1, 'Field1', MeshTopology, component_length)
   
   
   
   Call SectionRealSet(Field1%Sec, 1.0_Kr, iErr); CHKERRQ(iErr)
   Call SectionRealSet(Field1%Component_Sec(1), 2.23_Kr, iErr); CHKERRQ(iErr)
   Call SectionRealSet(Field1%Component_Sec(2), 4.43_Kr, iErr); CHKERRQ(iErr)

   Call PetscMemorySetGetMaximumUsage(iErr); CHKERRQ(iErr)
   iDebug = 0
   Write(stagename, "(A)") "Complete test"
   Call ALEStagePush(stagename, iDebug, iErr); CHKERRQ(iErr)
   Do i = 1, n

      val = i+1.23_Kr
      Call SectionRealSet(Field1%Sec, val, ierr); CHKERRQ(iErr)
      Call SectionRealComplete(Field1%Sec, iErr); CHKERRQ(iErr)

      If (mod(i, 1000)==0) Then
         Call ALEStagePrintMemory(stagename, iErr); CHKERRQ(iErr)

         Call PetscMemoryShowUsage(PETSC_VIEWER_STDOUT_WORLD, "PetscMemoryShowUsage output for PETSC_COMM_WORLD: ", iErr); CHKERRQ(iErr)
      End If

   End Do
   Call ALEStagePop(iDebug, iErr); CHKERRQ(iErr)


   Call FieldDestroy(Field1)
300 Format(A)   
   Call MEF90_Finalize()

End Program FindMPIProblem
