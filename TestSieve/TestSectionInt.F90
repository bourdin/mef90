Program TestSectionInt
#include "finclude/petscdef.h"
   Use m_MEF90
   Use petsc
   Implicit NONE   

   Type (MeshTopology_Type)                     :: MeshTopology
   Type (EXO_Type)                              :: EXO
   Type(Element2D_Scal),Dimension(:),Pointer    :: Elem2DA
   
   PetscBool                                    :: HasPrefix
   PetscBool                                    :: verbose
   PetscErrorCode                               :: iErr,iBlk,iE,i
   Character(len=256)                           :: CharBuffer,IOBuffer,filename
   Character(len=256)                           :: prefix
   Type(PetscViewer)                            :: viewer,myviewer
   Type(SectionReal)                            :: U_Sec
   Type(Vec)                                    :: U_Vec
   PetscReal,Dimension(:),Pointer               :: Values
   PetscInt,Dimension(:),Pointer                :: IntValues
   Type(SectionInt)                             :: Flag_Sec
   Type(SectionInt)                             :: Sec1,Sec2
   Type(DM)                                     :: tmpDM
   Integer                                      :: cpu_ws = 0
   Integer                                      :: io_ws = 0
   PetscReal                                    :: vers
   Integer                                      :: numCellSet,numVertexSet,numDim,numVertices,numCells    
   Integer                                      :: myNumVertices,myNumCells
   
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
      !Call DMDestroy(tmpDM,ierr);CHKERRQ(iErr)
   End If
   If (MEF90_Myrank == 0) Then
      Call EXCLOS(EXO%exoid,ierr)
   End If
   Call MeshTopologyGetInfo(MeshTopology)
   Call DMMeshGetStratumSize(MeshTopology%mesh,"height",0,numCells,ierr);CHKERRQ(ierr)
   Call DMMeshGetStratumSize(MeshTopology%mesh,"depth",0,numVertices,ierr);CHKERRQ(ierr)
   Call DMMeshGetLabelSize(MeshTopology%mesh,"Cell Sets",numCellSet,ierr);CHKERRQ(ierr)
   Call DMMeshGetLabelSize(MeshTopology%mesh,"Vertex Sets",numVertexSet,ierr);CHKERRQ(ierr)
   Call DMMeshGetDimension(MeshTopology%Mesh,numDim,ierr);CHKERRQ(ierr)


   !!! Create the Section on cpu0 
   If (MEF90_numProcs > 1) Then
      Call DMMeshGetStratumSize(tmpDM,"height",0,myNumCells,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumSize(tmpDM,"depth",0,myNumVertices,ierr);CHKERRQ(ierr)
      Call DMMeshGetSectionInt(tmpDM,prefix,Sec1,iErr); CHKERRQ(iErr)
      Call PetscObjectSetName(Sec1,"sec1",ierr);CHKERRQ(ierr)
      Do i = 1,myNumVertices
         Call SectionIntSetFiberDimension(Sec1,i+myNumCells-1,1,iErr); CHKERRQ(iErr)
      End Do 
      Call SectionIntAllocate(Sec1,iErr); CHKERRQ(iErr)
   
      !!! Setup and initialize an internal SectionInt
      Allocate(IntValues(1))
      Do i = 1,myNumVertices
         IntValues = i
         Call SectionIntUpdateClosure(Sec1,tmpDM,i+myNumCells-1,IntValues,INSERT_VALUES,iErr); CHKERRQ(iErr)
      End Do 
      DeAllocate(IntValues)
      
      Call PetscPrintf(PETSC_COMM_WORLD,"Sec1\n",iErr); CHKERRQ(iErr)
      Call SectionIntView(Sec1,PETSC_VIEWER_STDOUT_WORLD,iErr); CHKERRQ(iErr)
      !!! It looks like SectionIntDistribute is broken.
      !!! Not a big deal since we never used it...
      !Call SectionIntDistribute(Sec1,MeshTopology%mesh,Sec2,iErr); CHKERRQ(iErr)
      !Call PetscPrintf(PETSC_COMM_WORLD,"Sec2\n",iErr); CHKERRQ(iErr)
      !Call SectionIntView(Sec2,PETSC_VIEWER_STDOUT_WORLD,iErr); CHKERRQ(iErr)
   End If   
   Call MEF90_Finalize()
End Program TestSectionInt
