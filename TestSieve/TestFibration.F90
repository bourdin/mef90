Program TestFibration


#include "finclude/petscdef.h"

   Use m_MEF90
   Use petsc

   Implicit NONE   

   Type (MeshTopology_Type)                     :: MeshTopology
   Type (EXO_Type)                              :: EXO
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
   Integer                                      :: cpu_ws = 0
   Integer                                      :: io_ws = 0
   Integer                                      :: rank,numproc
   Type(DM)                                     :: tmpDM
   Integer                                      :: numCellSet,numDim
   Integer                                      :: exoid
   PetscReal                                    :: vers

   Call MEF90_Initialize()
   Call MPI_Comm_size(PETSC_COMM_WORLD,numproc,iErr);CHKERRQ(iErr)
   Call MPI_Comm_rank(PETSC_COMM_WORLD,rank,iErr);CHKERRQ(iErr)
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p', prefix, HasPrefix, iErr)    
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No input file prefix given\n", iErr)
      Call MEF90_Finalize()
      STOP
   End If

   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'


   If (rank == 0) Then
      exoid = EXOPEN(EXO%filename,EXREAD,cpu_ws,io_ws,vers,ierr)
   End If
   If (numproc == 1) Then
      Call DMMeshCreateExodusNG(PETSC_COMM_WORLD,exoid,MeshTopology%mesh,ierr);CHKERRQ(ierr)
   Else
      Call DMMeshCreateExodusNG(PETSC_COMM_WORLD,exoid,tmpDM,ierr);CHKERRQ(ierr)
      Call DMMeshDistribute(tmpDM,PETSC_NULL_CHARACTER,MeshTopology%mesh,ierr);CHKERRQ(iErr)
      Call DMDestroy(tmpDM,ierr);CHKERRQ(iErr)
   End If
   
   Write(IOBuffer, *) "Initializing MeshTopology object\n"
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call MeshTopologyGetInfo(MeshTopology)
   Call DMMeshGetLabelSize(MeshTopology%mesh,"Cell Sets",numCellSet,ierr);CHKERRQ(ierr)
   Call DMMeshGetDimension(MeshTopology%Mesh,numDim,ierr);CHKERRQ(ierr)
   
   Write(IOBuffer, *) "Initializing Element types\n"
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   MeshTopology%cellSet%ElemType    = MEF90_P1_Lagrange
   Do iBlk = 1, numCellSet
      Call cellSetElemTypeInit(MeshTopology%cellSet(iBlk), numDim)
   End Do
    

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

   If (rank == 0) Then
      Call EXCLOS(exoid,ierr)
   End If
   Call FieldDestroy(Field1)
   Call FlagDestroy(Flag1)
   Call MEF90_Finalize()

End Program TestFibration
