Program TestUpdateRestrict
#include "finclude/petscdef.h"
   Use m_MEF90
   Use petsc
   Implicit NONE   

   Type(MeshTopology_Type)                      :: MeshTopology
   Type(DM)                                     :: Tmp_Mesh
   Type(EXO_Type)                               :: EXO,MyEXO
   Type(Element2D_Scal),Dimension(:),Pointer    :: Elem2DA
   
   PetscBool                                    :: HasPrefix
   PetscBool                                    :: verbose
   PetscErrorCode                               :: iErr,iBlk,iE,i
   Character(len=256)                           :: CharBuffer,IOBuffer,filename
   Character(len=256)                           :: prefix
   Type(SectionReal)                            :: rSec
   Type(SectionInt)                             :: iSec

   PetscInt                                     :: Point
   PetscInt                                     :: SecSize

   PetscReal                                    :: rVal
   PetscReal,Dimension(:),Pointer               :: rVals
   PetscInt                                     :: iVal
   PetscInt,Dimension(:),Pointer                :: iVals
   
   PetscInt                                     :: dof = 1
   Integer                                      :: numCellSet,numVertexSet,numDim,numVertices,numCells    
   Type(DM)                                     :: tmpDM
   Integer                                      :: cpu_ws = 0
   Integer                                      :: io_ws = 0
   PetscReal                                    :: vers

     
   Call MEF90_Initialize()
   Call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-verbose',verbose,iErr)    
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-p',prefix,HasPrefix,iErr)    
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD,"No input file prefix given\n",iErr)
      Call MEF90_Finalize()
      STOP
   End If

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
   
   Write(IOBuffer, *) "Initializing MeshTopology object\n"
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call MeshTopologyGetInfo(MeshTopology)
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

   Call DMMeshGetVertexSectionReal(MeshTopology%Mesh,'rSec',dof,rSec,iErr); CHKERRQ(iErr)
   Call DMMeshGetVertexSectionInt (MeshTopology%Mesh,'iSec',dof,iSec,iErr); CHKERRQ(iErr)
   
   Write(IOBuffer,100) MEF90_MyRank,numCells,numVertices,numCellSet,numVertexSet
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)
100 Format("[",I2,"]: number of cells:        ",I5,"\n      number of vertices:     ",I5,&
    &      "\n      number of cell sets:   ",I5,"\n      number of vertex sets: ",I5,"\n")


!!! Initializing Sections
   Write(IOBuffer,*) '\n\n === Initializing Sections ===\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   rVal = 1.00_Kr
   Call SectionRealSet(rSec,rVal,iErr); CHKERRQ(iErr)
   Call SectionIntZero (iSec,iErr); CHKERRQ(iErr)

   Write(IOBuffer,*) 'rSec after SectionRealSet: \n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   Call SectionRealView(rSec,PETSC_VIEWER_STDOUT_WORLD,iErr); CHKERRQ(iErr)

   Write(IOBuffer,*) 'iSec after SectionIntSet:  \n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   Call SectionIntView(iSec,PETSC_VIEWER_STDOUT_WORLD,iErr); CHKERRQ(iErr)

!!! Testing SectionRealUpdate / SectionIntUpdate
   Write(IOBuffer,*) '\n\n === Testing SectionRealUpdate / SectionIntUpdate ===\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   
   Point = numCells
   Allocate(rVals(dof))
   rVals = -2.00_Kr
   Call SectionRealUpdate(rSec,Point,rVals,INSERT_VALUES,iErr); CHKERRQ(iErr)

   Point = numCells+1
   rVals = 2.00_Kr
   Call SectionRealUpdate(rSec,Point,rVals,ADD_VALUES,iErr); CHKERRQ(iErr)

   Write(IOBuffer,*) 'rSec after SectionRealUpdate:  \n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   Call SectionRealView(rSec,PETSC_VIEWER_STDOUT_WORLD,iErr); CHKERRQ(iErr)
   DeAllocate(rVals)

   Point = numCells
   Allocate(iVals(dof))
   iVals = -7
   Call SectionIntUpdate(iSec,Point,iVals,INSERT_VALUES,iErr); CHKERRQ(iErr)

   Point = numCells+1
   iVals = 4
   Call SectionIntUpdate(iSec,Point,iVals,ADD_VALUES,iErr); CHKERRQ(iErr)
   Write(IOBuffer,*) 'iSec after SectionIntUpdate:  \n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   Call SectionIntView(iSec,PETSC_VIEWER_STDOUT_WORLD,iErr); CHKERRQ(iErr)
   DeAllocate(iVals)
   
!!! Testing SectionRealRestrict / SectionIntRestrict
   Write(IOBuffer,*) '\n\n === SectionRealRestrict / SectionIntRestrict ===\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)

   Point = numCells
   Call SectionRealRestrict(rSec,Point,rVals,iErr); CHKERRQ(iErr)
   Write(IOBuffer,*) 'rSec(0)=',rVals,'\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   Call SectionRealRestore(rSec,Point,rVals,iErr); CHKERRQ(iErr)

   Point = numCells+1
   Call SectionRealRestrict(rSec,Point,rVals,iErr); CHKERRQ(iErr)
   Write(IOBuffer,*) ' rSec(1)=',rVals,'\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   Call SectionRealRestore(rSec,Point,rVals,iErr); CHKERRQ(iErr)
   

   Point = numCells
   Call SectionIntRestrict(iSec,Point,iVals,iErr); CHKERRQ(iErr)
   Write(IOBuffer,*) 'iSec(0)=',iVals,'\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   Call SectionIntRestore(iSec,Point,iVals,iErr); CHKERRQ(iErr)

   Point = numCells+1
   Call SectionIntRestrict(iSec,Point,iVals,iErr); CHKERRQ(iErr)
   Write(IOBuffer,*) ' iSec(1)=',iVals,'\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   Call SectionIntRestore(iSec,Point,iVals,iErr); CHKERRQ(iErr)


!!! Testing SectionRealUpdateClosure / SectionIntUpdateClosure
   Write(IOBuffer,*) '\n\n === SectionRealUpdateClosure / SectionIntUpdateClosure ===\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)

   Allocate(rVals(3*dof))
   Point = 2
   rVals = 6.5_Kr
   Call SectionRealUpdateClosure(rSec,MeshTopology%Mesh,Point,rVals,INSERT_VALUES,iErr); CHKERRQ(iErr)

   Point = 18
   rVals = 3.0_Kr
   Call SectionRealUpdateClosure(rSec,MeshTopology%Mesh,Point,rVals,ADD_VALUES,iErr); CHKERRQ(iErr)
   DeAllocate(rVals)

   Write(IOBuffer,*) 'rSec after SectionRealUpdateClosure:  \n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   Call SectionRealView(rSec,PETSC_VIEWER_STDOUT_WORLD,iErr); CHKERRQ(iErr)

   Allocate(iVals(3*dof))
   Point = 2
   iVals = -9
   Call SectionIntUpdateClosure(iSec,MeshTopology%Mesh,Point,iVals,INSERT_VALUES,iErr); CHKERRQ(iErr)

   Point = 18
   iVals = 4
   Call SectionIntUpdateClosure(iSec,MeshTopology%Mesh,Point,iVals,ADD_VALUES,iErr); CHKERRQ(iErr)
   Write(IOBuffer,*) 'iSec after SectionIntUpdate:  \n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   Call SectionIntView(iSec,PETSC_VIEWER_STDOUT_WORLD,iErr); CHKERRQ(iErr)
   DeAllocate(iVals)

!!! Testing SectionRealRestrictClosure / SectionIntRestrictClosure
   Write(IOBuffer,*) '\n\n === SectionRealRestrictClosure / SectionIntRestrictClosure ===\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)

   Allocate(rVals(3*dof))
   Point = 0
   Call SectionRealRestrictClosure(rSec,MeshTopology%Mesh,Point,3*dof,rVals,iErr); CHKERRQ(iErr)
   Write(IOBuffer,*) 'rSec(0)=',rVals,'\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)

   Point = 14
   Call SectionRealRestrictClosure(rSec,MeshTopology%Mesh,Point,3*dof,rVals,iErr); CHKERRQ(iErr)
   Write(IOBuffer,*) ' rSec(1)=',rVals,'\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   DeAllocate(rVals)
   

   Allocate(iVals(3*dof))
   Point = 0
   Call SectionIntRestrictClosure(iSec,MeshTopology%Mesh,Point,3*dof,iVals,iErr); CHKERRQ(iErr)
   Write(IOBuffer,*) 'iSec(0)=',iVals   ,'\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)

   Point = 14
   Call SectionIntRestrictClosure(iSec,MeshTopology%Mesh,Point,3*dof,iVals,iErr); CHKERRQ(iErr)
   Write(IOBuffer,*) ' iSec(1)=',iVals,'\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   DeAllocate(iVals)
   
!!! Testing SectionRealGetSize
   Call SectionRealGetSize(rSec,SecSize,iErr); CHKERRQ(iErr)
   Write(IOBuffer,*) 'size(rSec)=',SecSize,'\n'   
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,iErr); CHKERRQ(iErr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD,iErr); CHKERRQ(iErr)
   

   Call MeshTopologyDestroy(MeshTopology)
   Call SectionRealDestroy(rSec,iErr); CHKERRQ(iErr)
   Call SectionIntDestroy (iSec,iErr); CHKERRQ(iErr)
   Call MEF90_Finalize()
End Program TestUpdateRestrict
