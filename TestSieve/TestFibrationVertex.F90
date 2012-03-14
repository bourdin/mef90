Program TestFibration


#include "finclude/petscdef.h"

   Use m_MEF90
   Use petsc

   Implicit NONE   

   Type(MeshTopology_Type)                      :: MeshTopology
   Type(EXO_Type)                               :: EXO, MyEXO
   Type(Element2D_Elast), Dimension(:), Pointer :: ElemVect
   Type(Element2D_Scal), Dimension(:), Pointer  :: ElemScal
   Type(Field)                                  :: Field1, Field2
   
   PetscBool                                    :: HasPrefix, flg
   PetscInt                                     :: verbose
   PetscErrorCode                               :: iErr, i, j, iBlk
   Character(len=256)                           :: CharBuffer, IOBuffer, filename
   Character(len=256)                           :: prefix
   Type(PetscViewer)                            :: LogViewer, MyLogViewer
   PetscInt                                     :: num_components, num_dof
   PetscInt, Dimension(:), Pointer              :: component_length 
   Type(SectionReal)                            :: Comp1, Comp2
   Integer                                      :: cpu_ws = 0
   Integer                                      :: io_ws = 0
   Integer                                      :: rank,numproc
   Type(DM)                                     :: tmpDM
   Integer                                      :: numCellSet,numDim,numVertices,numCells
   Integer                                      :: exoid
   PetscReal                                    :: vers

   Call MEF90_Initialize()
   Call MPI_Comm_size(PETSC_COMM_WORLD,numproc,iErr);CHKERRQ(iErr)
   Call MPI_Comm_rank(PETSC_COMM_WORLD,rank,iErr);CHKERRQ(iErr)
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
   
   If (verbose > 0) Then
      Write(IOBuffer, *) "Initializing MeshTopology object\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   Call MeshTopologyGetInfo(MeshTopology)
   Call DMMeshGetLabelSize(MeshTopology%mesh,"Cell Sets",numCellSet,ierr);CHKERRQ(ierr)
   Call DMMeshGetDimension(MeshTopology%Mesh,numDim,ierr);CHKERRQ(ierr)
   
   If (verbose > 0) Then
      Write(IOBuffer, *) "Initializing Element types\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   MeshTopology%cellSet%ElemType    = MEF90_P1_Lagrange
   Do iBlk = 1, numCellSet
      Call cellSetElemTypeInit(MeshTopology%cellSet(iBlk), numDim)
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
   
   Call DMMeshGetVertexSectionReal(MeshTopology%mesh, 'Field1.Sec', num_dof, Field1%Sec, iErr); CHKERRQ(iErr)
   Do i = 1, num_components
      Call SectionRealAddSpace(Field1%Sec, iErr); CHKERRQ(iErr)
   End Do 
!!! For some reason, it looks like one doesn't need to call SectionRealAllocate if the Section was created with MeshGetVertexSectionReal
   Call SectionRealAllocate(Field1%Sec, iErr); CHKERRQ(iErr)

   Call DMMeshGetStratumSize(MeshTopology%mesh,"depth",0,numVertices,ierr);CHKERRQ(ierr)
   Call DMMeshGetStratumSize(MeshTopology%mesh,"height",0,numCells,ierr);CHKERRQ(ierr)
   Do i = 1, numVertices
      Do j = 1, num_components
         Call SectionRealSetFiberDimensionField(Field1%Sec, i+numCells-1, component_length(j), j-1, iErr); CHKERRQ(iErr)
      End Do
   End Do 

   Allocate(Field1%Component_Sec(Num_Components))
   Do i = 1, Num_Components
      Call SectionRealGetFibration(Field1%Sec, i-1, Field1%Component_sec(i), iErr); CHKERRQ(iErr)
   End Do

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
   
   Call PetscPrintf(PETSC_COMM_WORLD, "Now trying to scatter into a Vec",iErr);CHKERRQ(iErr)
   Call DMMeshCreateGlobalScatter(MeshTopology%mesh,Field1%Sec,Field1%Scatter,ierr);CHKERRQ(ierr)
   Call DMMeshCreateVector(MeshTopology%mesh,Field1%Sec,Field1%Vec,ierr);CHKERRQ(ierr)
   Call SectionRealToVec(Field1%Sec,Field1%Scatter,SCATTER_FORWARD,Field1%Vec,iErr);CHKERRQ(iErr);
   Call VecView(Field1%Vec,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(iErr)

   If (rank == 0) Then
      Call EXCLOS(exoid,ierr)
   End If
   Call MEF90_Finalize()
End Program TestFibration
