Program TestSectionToArray


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
   Type(Field)                                  :: Field1
   Type(Vec)                                    :: LocalVec
         
   PetscTruth                                   :: HasPrefix, flg
   PetscErrorCode                               :: iErr, i, j, iBlk
   Character(len=256)                           :: CharBuffer, IOBuffer, filename
   Character(len=256)                           :: prefix
   Type(Mesh)                                   :: Tmp_Mesh
   PetscInt                                     :: num_components, num_dof
   PetscInt, Dimension(:), Pointer              :: component_length 
   PetscInt                                     :: MySize

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
   Call MeshTopologyReadEXO(MeshTopology, EXO)
   
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

   Write(IOBuffer, *) "Field1.Sec\n"
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call SectionRealView(Field1%Sec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   
   
   Write(IOBuffer, *) "Field1.Component_Sec(1)\n"
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call SectionRealView(Field1%Component_Sec(1), PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   
   Write(IOBuffer, *) "Field1.Component_Sec(2)\n"
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call SectionRealView(Field1%Component_Sec(2), PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   
   Call SectionRealToVec(Field1%Sec, Field1%Scatter, SCATTER_FORWARD, Field1%Vec, iErr); CHKERRQ(iErr)
   Write(IOBuffer, *) "Field1.Vec\n"
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call VecView(Field1%Vec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   
   Call SectionRealGetSize(Field1%Sec, MySize, iErr); CHKERRQ(iErr)
   Write(IOBuffer, *) 'local size of Field1./.Sec: ', MySize, '\n'
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD, IErr); CHKERRQ(iErr)

   Call SectionRealGetSize(Field1%Component_Sec(1), MySize, iErr); CHKERRQ(iErr)
   Write(IOBuffer, *) 'local size of Field1./.Component_Sec(1): ', MySize, '\n'
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD, IErr); CHKERRQ(iErr)

   Call SectionRealGetSize(Field1%Component_Sec(2), MySize, iErr); CHKERRQ(iErr)
   Write(IOBuffer, *) 'local size of Field1./.Component_Sec(2): ', MySize, '\n'
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD, IErr); CHKERRQ(iErr)
   
   Call SectionRealCreateLocalVector(Field1%Sec, LocalVec, iErr); CHKERRQ(iErr)
   Call VecGetSize(LocalVec, MySize, iErr); CHKERRQ(iErr)
   Write(IOBuffer, *) 'local size of LocalVec obtained from Field1./.Sec: ', MySize, '\n'
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD, IErr); CHKERRQ(iErr)
!   Call VecView(LocalVec, PETSC_VIEWER_STDOUT_SELF, iErr); CHKERRQ(iErr)
   Call VecDestroy(LocalVec, iErr); CHKERRQ(iErr)

   Call SectionRealCreateLocalVector(Field1%Component_Sec(1), LocalVec, iErr); CHKERRQ(iErr)
   Call VecGetSize(LocalVec, MySize, iErr); CHKERRQ(iErr)
   Write(IOBuffer, *) 'local size of LocalVec obtained from Field1./.Component_Sec(1): ', MySize, '\n'
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD, IErr); CHKERRQ(iErr)
!   Call VecView(LocalVec, PETSC_VIEWER_STDOUT_SELF, iErr); CHKERRQ(iErr)
   Call VecDestroy(LocalVec, iErr); CHKERRQ(iErr)

   Call SectionRealCreateLocalVector(Field1%Component_Sec(2), LocalVec, iErr); CHKERRQ(iErr)
   Call VecGetSize(LocalVec, MySize, iErr); CHKERRQ(iErr)
   Write(IOBuffer, *) 'local size of LocalVec obtained from Field1./.Component_Sec(2): ', MySize, '\n'
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD, IErr); CHKERRQ(iErr)
!   Call VecView(LocalVec, PETSC_VIEWER_STDOUT_SELF, iErr); CHKERRQ(iErr)
   Call VecDestroy(LocalVec, iErr); CHKERRQ(iErr)
   
   !!! Now trying to gain access to the section directly through a fortran array
   Call VecGetSize(Field1%LocalVec, MySize, iErr); CHKERRQ(iErr)
   Write(IOBuffer, *) 'local size of Field1./.LocalVec: ', MySize, '\n'
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD, IErr); CHKERRQ(iErr)
!   Call VecView(Field1%LocalVec, PETSC_VIEWER_STDOUT_SELF, iErr); CHKERRQ(iErr)


   Call FieldDestroy(Field1)
300 Format(A)   
   Call MEF90_Finalize()

End Program TestSectionToArray