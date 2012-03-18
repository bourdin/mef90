Program TestLocal

#include "finclude/petscdef.h"
   Use m_MEF90
   Use petsc
   Implicit NONE   

   Type (MeshTopology_Type)                     :: MeshTopology
   Type (EXO_Type)                              :: EXO
   Type (Element2D_Scal),Dimension(:),Pointer   :: Elem2DA
   Type (Vect3D),Dimension(:),Pointer           :: Coords
   
   PetscBool                                    :: HasPrefix
   PetscBool                                    :: verbose
   PetscErrorCode                               :: ierr
   PetscInt                                     :: iBlk
   Character(len=256)                           :: CharBuffer,IOBuffer,filename
   Character(len=256)                           :: prefix
   Type(PetscViewer)                            :: viewer,myviewer
   Type(SectionReal)                            :: Coords_Sec,U_Sec
   Type(Vec)                                    :: Coords_VecG,Coords_VecL
   PetscReal,Dimension(:),Pointer               :: Coords_Ptr
   PetscInt                                     :: VSize
   Type(VecScatter)                             :: scatter
   PetscInt                                     :: i
   PetscReal                                    :: T
   Type(Vect2D),Dimension(:),Pointer            :: V2D
   Type(Vect3D),Dimension(:),Pointer            :: V3D
   PetscReal,Dimension(:),Pointer               :: V_Ptr
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
   Call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-verbose',verbose,ierr)    
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-p',prefix,HasPrefix,ierr)    
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD,"No input file prefix given\n",ierr)
      Call MEF90_Finalize()
      STOP
   End If

   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'
   EXO%exoid = EXOPEN(EXO%filename,EXREAD,cpu_ws,io_ws,vers,ierr)

   If (numproc == 1) Then
      Call DMMeshCreateExodusNG(PETSC_COMM_WORLD,EXO%exoid,MeshTopology%mesh,ierr);CHKERRQ(ierr)
   Else
      Call DMMeshCreateExodusNG(PETSC_COMM_WORLD,EXO%exoid,tmpDM,ierr);CHKERRQ(ierr)
      Call DMMeshDistribute(tmpDM,PETSC_NULL_CHARACTER,MeshTopology%mesh,ierr);CHKERRQ(iErr)
      Call DMDestroy(tmpDM,ierr);CHKERRQ(iErr)
   End If
   Call MeshTopologyGetInfo(MeshTopology)
   Call DMMeshGetLabelSize(MeshTopology%mesh,"Cell Sets",numCellSet,ierr);CHKERRQ(ierr)
   Call DMMeshGetDimension(MeshTopology%Mesh,numDim,ierr);CHKERRQ(ierr)
   Call EXOProperty_Read(EXO)
   Call EXOVariable_Read(EXO)
   
   MeshTopology%cellSet%ElemType    = MEF90_P1_Lagrange
   Do iBlk = 1, numCellSet
      Call cellSetElemTypeInit(MeshTopology%cellSet(iBlk), numDim)
   End Do
   
   If (verbose) Then
      Write(filename,102) Trim(prefix),MEF90_MyRank
      Call PetscViewerASCIIOpen(PETSC_COMM_SELF,filename,myviewer,ierr);CHKERRQ(ierr);  
      Write(IOBuffer,103) MEF90_MyRank,Trim(filename)
      Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Call PetscSynchronizedFlush(PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)

      Write(filename,101) Trim(prefix)
      Call PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,viewer,ierr);CHKERRQ(ierr);  
      Write(IOBuffer,104) Trim(filename)
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
   End If
   
 
   If (verbose) Then
      Write(IOBuffer,200) 'EXO\n'
      Call PetscViewerASCIIPrintf(myviewer,IOBuffer,ierr);CHKERRQ(ierr)
      Call EXOView(EXO,myviewer)
   End If
   
   Call EXCLOS(EXO%exoid,ierr)
   Call MEF90_Finalize()
101 Format(A,'.log')
102 Format(A,'-',I4.4,'.log')
103 Format('Output from processor ',I4.4,' redirected to ',A,'\n')
104 Format('Collective output redirected to ',A,'\n')
200 Format(A)
End Program TestLocal
