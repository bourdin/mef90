Program TestSectionToArray
#include "finclude/petscdef.h"
   Use m_MEF90
   Use petsc
   Implicit NONE   

   Type(MeshTopology_Type)                      :: MeshTopology
   Type(EXO_Type)                               :: EXO
   Type(Field)                                  :: Field1
   Type(Vec)                                    :: LocalVec
         
   PetscBool                                    :: HasPrefix,flg
   PetscErrorCode                               :: iErr,i,j,iBlk
   Character(len=256)                           :: CharBuffer,IOBuffer,filename
   Character(len=256)                           :: prefix
   PetscInt                                     :: num_components,num_dof
   PetscInt,Dimension(:),Pointer                :: component_length 
   Type(DM)                                     :: tmpDM
   Integer                                      :: cpu_ws = 0
   Integer                                      :: io_ws = 0
   PetscReal                                    :: vers
   Integer                                      :: numCellSet,numVertexSet,numDim,numVertices,numCells    
   Integer                                      :: mySize

   Call MEF90_Initialize()
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-p',prefix,HasPrefix,iErr)    
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD,"No input file prefix given\n",iErr)
      Call MEF90_Finalize()
      STOP
   End If

   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'


   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'
   If (MEF90_Myrank == 0) Then
      EXO%exoid = EXOPEN(EXO%filename,EXREAD,cpu_ws,io_ws,vers,ierr)
   End If
   If (MEF90_numprocs == 1) Then
      Call DMMeshCreateExodusNG(PETSC_COMM_WORLD,EXO%exoid,MeshTopology%mesh,ierr);CHKERRQ(ierr)
   Else
      Call DMMeshCreateExodusNG(PETSC_COMM_WORLD,EXO%exoid,tmpDM,ierr);CHKERRQ(ierr)
      Call DMMeshDistribute(tmpDM,PETSC_NULL_CHARACTER,MeshTopology%mesh,ierr);CHKERRQ(iErr)
      Call DMDestroy(tmpDM,ierr);CHKERRQ(iErr)
   End If
   If (MEF90_Myrank == 0) Then
      Call EXCLOS(EXO%exoid,ierr)
   End If
   Call MeshTopologyGetInfo(MeshTopology)
   
   !!! Get various local sizes
   Call DMMeshGetStratumSize(MeshTopology%mesh,"height",0,numCells,ierr);CHKERRQ(ierr)
   Call DMMeshGetStratumSize(MeshTopology%mesh,"depth",0,numVertices,ierr);CHKERRQ(ierr)
   Call DMMeshGetLabelSize(MeshTopology%mesh,"Cell Sets",numCellSet,ierr);CHKERRQ(ierr)
   Call DMMeshGetLabelSize(MeshTopology%mesh,"Vertex Sets",numVertexSet,ierr);CHKERRQ(ierr)
   Call DMMeshGetDimension(MeshTopology%Mesh,numDim,ierr);CHKERRQ(ierr)

   Write(IOBuffer, *) "Initializing Element types\n"
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   MeshTopology%cellSet%ElemType    = MEF90_P1_Lagrange
   Do iBlk = 1, numCellSet
      Call cellSetElemTypeInit(MeshTopology%cellSet(iBlk), numDim)
   End Do     

!!! Creating the Sec component of the field and flags
   Write(IOBuffer,*) "Creating Sections\n"
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)

   num_components=2
   Allocate (component_length(2))
   component_length(1) = 1
   component_length(2) = 2
   num_dof = sum(component_length)

   Write(IOBuffer,*) "Field1.Sec\n"
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   
   Call FieldCreateVertex(Field1,'Field1',MeshTopology,component_length)
   
   Call SectionRealSet(Field1%Sec,1.0_Kr,iErr); CHKERRQ(iErr)
   Call SectionRealSet(Field1%Component_Sec(1),2.23_Kr,iErr); CHKERRQ(iErr)
   Call SectionRealSet(Field1%Component_Sec(2),4.43_Kr,iErr); CHKERRQ(iErr)

   Write(IOBuffer,*) "Field1.Sec\n"
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   Call SectionRealView(Field1%Sec,PETSC_VIEWER_STDOUT_WORLD,iErr); CHKERRQ(iErr)
   
   
   Write(IOBuffer,*) "Field1.Component_Sec(1)\n"
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   Call SectionRealView(Field1%Component_Sec(1),PETSC_VIEWER_STDOUT_WORLD,iErr); CHKERRQ(iErr)
   
   Write(IOBuffer,*) "Field1.Component_Sec(2)\n"
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   Call SectionRealView(Field1%Component_Sec(2),PETSC_VIEWER_STDOUT_WORLD,iErr); CHKERRQ(iErr)
   
   Call SectionRealToVec(Field1%Sec,Field1%Scatter,SCATTER_FORWARD,Field1%Vec,iErr); CHKERRQ(iErr)
   Write(IOBuffer,*) "Field1.Vec\n"
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   Call VecView(Field1%Vec,PETSC_VIEWER_STDOUT_WORLD,iErr); CHKERRQ(iErr)
   
   Call SectionRealGetSize(Field1%Sec,mySize,iErr); CHKERRQ(iErr)
   Write(IOBuffer,*) 'local size of Field1./.Sec: ',mySize,'\n'
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD,IErr); CHKERRQ(iErr)

   Call SectionRealGetSize(Field1%Component_Sec(1),mySize,iErr); CHKERRQ(iErr)
   Write(IOBuffer,*) 'local size of Field1./.Component_Sec(1): ',mySize,'\n'
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD,IErr); CHKERRQ(iErr)

   Call SectionRealGetSize(Field1%Component_Sec(2),mySize,iErr); CHKERRQ(iErr)
   Write(IOBuffer,*) 'local size of Field1./.Component_Sec(2): ',mySize,'\n'
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD,IErr); CHKERRQ(iErr)
   
   Call SectionRealCreateLocalVector(Field1%Sec,LocalVec,iErr); CHKERRQ(iErr)
   Call VecGetSize(LocalVec,mySize,iErr); CHKERRQ(iErr)
   Write(IOBuffer,*) 'local size of LocalVec obtained from Field1./.Sec: ',mySize,'\n'
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD,IErr); CHKERRQ(iErr)
   Call VecDestroy(LocalVec,iErr); CHKERRQ(iErr)

   Call SectionRealCreateLocalVector(Field1%Component_Sec(1),LocalVec,iErr); CHKERRQ(iErr)
   Call VecGetSize(LocalVec,mySize,iErr); CHKERRQ(iErr)
   Write(IOBuffer,*) 'local size of LocalVec obtained from Field1./.Component_Sec(1): ',mySize,'\n'
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD,IErr); CHKERRQ(iErr)
   Call VecDestroy(LocalVec,iErr); CHKERRQ(iErr)

   Call SectionRealCreateLocalVector(Field1%Component_Sec(2),LocalVec,iErr); CHKERRQ(iErr)
   Call VecGetSize(LocalVec,mySize,iErr); CHKERRQ(iErr)
   Write(IOBuffer,*) 'local size of LocalVec obtained from Field1./.Component_Sec(2): ',mySize,'\n'
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD,IErr); CHKERRQ(iErr)
   Call VecDestroy(LocalVec,iErr); CHKERRQ(iErr)

   Call FieldDestroy(Field1)
   Call MEF90_Finalize()
End Program TestSectionToArray