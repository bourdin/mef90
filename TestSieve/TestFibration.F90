Program TestFibration


#include "finclude/petscdef.h"

   Use m_MEF90
   Use petsc

   Implicit NONE   

   Type (MeshTopology_Type)                     :: MeshTopology
   Type (EXO_Type)                              :: EXO, MyEXO
   Type(Element2D_Elast), Dimension(:), Pointer :: ElemVect
   Type(Element2D_Scal), Dimension(:), Pointer  :: ElemScal
   Type(Field)                                  :: Field1
   Type(Flag)                                   :: Flag1
   Type(Mat)                                    :: M
   PetscReal, Dimension(:,:), Pointer           :: MatElem
      
   PetscBool                                    :: HasPrefix, flg
   PetscErrorCode                               :: iErr, i, j, iBlk
   Character(len=256)                           :: CharBuffer, IOBuffer, filename
   Character(len=256)                           :: prefix
   PetscInt                                     :: num_components, num_dof
   PetscInt, Dimension(:), Pointer              :: component_length 

   Call MEF90_Initialize()
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p', prefix, HasPrefix, iErr)    
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No input file prefix given\n", iErr)
      Call MEF90_Finalize()
      STOP
   End If

   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'


   Call DMMeshCreateExodusNG(PETSC_COMM_WORLD, EXO%filename, MeshTopology%mesh, MeshTopology%meshFS,ierr); CHKERRQ(iErr)
   
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
   

   Call FlagCreateVertex(Flag1, 'Flag1', MeshTopology, component_length)
   
   Call SectionIntSet(Flag1%Sec, 1, iErr); CHKERRQ(iErr)
   Call SectionIntSet(Flag1%Component_Sec(1), 2, iErr); CHKERRQ(iErr)
   Call SectionIntSet(Flag1%Component_Sec(2), 4, iErr); CHKERRQ(iErr)

   Write(IOBuffer, *) "Flag1.Sec\n"
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call SectionIntView(Flag1%Sec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   
   
   Write(IOBuffer, *) "Flag1.Component_Sec(1)\n"
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call SectionIntView(Flag1%Component_Sec(1), PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   
   Write(IOBuffer, *) "Flag1.Component_Sec(2)\n"
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call SectionIntView(Flag1%Component_Sec(2), PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)


!   Call MeshSetMaxDof(MeshTopology%Mesh, num_dof, iErr); CHKERRQ(iErr) 
!   Call MeshCreateMatrix(MeshTopology%mesh, Field1%Sec, MATMPIAIJ, M, iErr); CHKERRQ(iErr)
!
!   j = 5
!   Do i = 1, num_components
!      Allocate(MatElem(component_length(i),component_length(i)))
!      MatElem = 100.0_Kr * i 
!      Call DMMeshAssembleMatrix(M, MeshTopology%mesh, Field1%Component_Sec(i), MeshTopology%Num_Elems+j-1, MatElem, ADD_VALUES, iErr); CHKERRQ(iErr)
!      DeAllocate(MatElem)
!   End Do
!   
!   Call MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
!   Call MatAssemblyEnd  (M, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
!
!   Call MatView(M, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
!   
!   Call MatDestroy(M, iErr); CHKERRQ(iErr)
   Call FieldDestroy(Field1)
   Call FlagDestroy(Flag1)
300 Format(A)   
   Call MEF90_Finalize()

End Program TestFibration
