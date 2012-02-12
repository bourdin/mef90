Program TestLocal
#include "finclude/petscdef.h"
   Use m_MEF90
   Use petsc
   Implicit NONE   

   Type (MeshTopology_Type)                     :: MeshTopology
   Type (EXO_Type)                              :: EXO,MyEXO
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
   PetscInt                                     :: vers

     
   Call MEF90_Initialize()
   Call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-verbose',verbose,ierr)    
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-p',prefix,HasPrefix,ierr)    
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD,"No input file prefix given\n",ierr)
      Call MEF90_Finalize()
      STOP
   End If

   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'
   EXO%exoid = EXOPEN(EXO%filename,EXREAD,exo_cpu_ws,exo_io_ws,vers,ierr)


   Call DMMeshCreateExodusNG(EXO%Comm,EXO%filename,MeshTopology%mesh,MeshTopology%meshFS,ierr);CHKERRQ(ierr)
   Call MeshTopologyGetInfo(MeshTopology,PETSC_COMM_WORLD)
   

   Call EXOProperty_Read(EXO)
   Call EXOVariable_Read(EXO)
   
   MeshTopology%Elem_Blk%Elem_Type = MEF90_P1_Lagrange
   Do iBlk = 1,MeshTopology%Num_Elem_Blks
      Call Init_Elem_Blk_Type(MeshTopology%Elem_Blk(iBlk),MeshTopology%num_dim)
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
   
101 Format(A,'.log')
102 Format(A,'-',I4.4,'.log')
103 Format('Output from processor ',I4.4,' redirected to ',A,'\n')
104 Format('Collective output redirected to ',A,'\n')
   MyEXO%comm = PETSC_COMM_SELF
   MyEXO%exoid = EXO%exoid
   Write(MyEXO%filename,200) trim(prefix),MEF90_MyRank
200 Format(A,'-',I4.4,'.gen')
 
   MyEXO%title  = trim(EXO%title)
   MyEXO%Num_QA = EXO%Num_QA
   MyEXO%exoid  = EXOPEN(MyEXO%filename,EXWRIT,exo_cpu_ws,exo_io_ws,vers,ierr)

   Call EXOProperty_Copy(EXO,MyEXO)
   Call EXOVariable_Copy(EXO,MyEXO)

   Call MeshTopologyWriteGlobal(MeshTopology,MyEXO,PETSC_COMM_WORLD)
   Call EXOProperty_Write(MyEXO)
   Call EXOVariable_Write(MyEXO)

   If (verbose) Then
      Write(IOBuffer,300) 'EXO\n'
      Call PetscViewerASCIIPrintf(myviewer,IOBuffer,ierr);CHKERRQ(ierr)
      Call EXOView(EXO,myviewer)
      Write(IOBuffer,300) 'MyEXO\n'
      Call PetscViewerASCIIPrintf(myviewer,IOBuffer,ierr);CHKERRQ(ierr)
      Call EXOView(MyEXO,myviewer)

      Write(IOBuffer,300) '\n\nMeshTopology\n'
      Call PetscViewerASCIIPrintf(myviewer,IOBuffer,ierr);CHKERRQ(ierr)
      Call MeshTopologyView(MeshTopology,myviewer);CHKERRQ(ierr)
   End If
   
      Call MEF90_Finalize()
      STOP
   
   !!! ****************** PLAYING WITH SECTIONS AND VECS ******************
   Call DMMeshGetSectionReal(MeshTopology%mesh,'coordinates',Coords_Sec,ierr);CHKERRQ(ierr) 
   If (verbose) Then
      Write(IOBuffer,300) '\n\nSECTION STUFF\n'
      Call PetscViewerASCIIPrintf(viewer,IOBuffer,ierr);CHKERRQ(ierr)
      Call PetscViewerASCIIPrintf(myviewer,IOBuffer,ierr);CHKERRQ(ierr)

      Call SectionRealView(Coords_Sec,viewer,ierr);CHKERRQ(ierr)
   End If


   Write(IOBuffer,400) MeshTopology%Num_Verts,MeshTopology%Num_Elems
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)
   
   If (verbose) Then
      Call PetscViewerASCIIPrintf(myviewer,IOBuffer,ierr);CHKERRQ(ierr)
   End If
400 Format('Number of Vertices: ',I3,' Number of elements: ',I3,'\n')


   Call DMMeshCreateVector(MeshTopology%mesh,Coords_Sec,Coords_VecG,ierr);CHKERRQ(ierr)

   Write(IOBuffer,*) '\n\nVec obtained with MeshCreateVector have no \"ghost values\" \n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)

   Call VecGetSize(Coords_VecG,VSize,ierr);CHKERRQ(ierr)
   Write(IOBuffer,401) VSize
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)

   Call VecGetLocalSize(Coords_VecG,VSize,ierr);CHKERRQ(ierr)
   Write(IOBuffer,402) VSize
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)
401 Format('   Global size of Coords_VecG: ',I3,'\n')
402 Format('   Local size of Coords_VecG: ',I3,'\n')

   If (verbose) Then
      !!! Scatter the values
      Write(IOBuffer,*) 'Coords_VecG before scattering\n'
      Call PetscViewerASCIIPrintf(viewer,IOBuffer,ierr);CHKERRQ(ierr)
      Call VecView(Coords_VecG,viewer,ierr);CHKERRQ(ierr)
      
      Call DMMeshCreateGlobalScatter(MeshTopology%mesh,Coords_Sec,scatter,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(Coords_Sec,scatter,SCATTER_FORWARD,Coords_VecG,ierr);CHKERRQ(ierr)
      Write(IOBuffer,*) 'Coords_VecG after scattering\n'
      Call PetscViewerASCIIPrintf(viewer,IOBuffer,ierr);CHKERRQ(ierr)
      Call VecView(Coords_VecG,viewer,ierr);CHKERRQ(ierr)
   End If

   Call SectionRealCreateLocalVector(Coords_Sec,Coords_VecL,ierr);CHKERRQ(ierr)
   If (verbose) Then
      !!! Scatter the values
      Write(IOBuffer,*) 'Coords_VecL\n'
      Call PetscViewerASCIIPrintf(myviewer,IOBuffer,ierr);CHKERRQ(ierr)
      Call VecView(Coords_VecL,myviewer,ierr);CHKERRQ(ierr)
   End If


   Write(IOBuffer,*) '\n\nVec obtained with SectionRealCreateLocalVector do have \"ghost values\" \n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)

   Call VecGetLocalSize(Coords_VecL,VSize,ierr);CHKERRQ(ierr)
   Write(IOBuffer,404) VSize
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)
403 Format('   Global size of Coords_VecL: ',I3,'\n')
404 Format('   Local size of Coords_VecL: ',I3,'\n')

   Call VecDestroy(Coords_VecG,ierr);CHKERRQ(ierr)
   Call VecDestroy(Coords_VecL,ierr);CHKERRQ(ierr)

   Call DMMeshGetVertexSectionReal(MeshTopology%mesh,2,U_Sec,ierr);CHKERRQ(ierr)

!!! Format the output mesh with variables
      MyEXO%exoid = EXOPEN(MyEXO%filename,EXWRIT,exo_cpu_ws,exo_io_ws,exo_ver,ierr)
      
      Call EXPVP (MyEXO%exoid,'n',4,ierr)
      Call EXPVAN (MyEXO%exoid,'n',4,(/'nRes1','nRes2','nRes3','nRes4'/),ierr)
      Call EXPVP (MyEXO%exoid,'e',4,ierr)
      Call EXPVAN (MyEXO%exoid,'e',4,(/'eRes1','eRes2','eRes3','eRes4'/),ierr)
      Call EXPVP (MyEXO%exoid,'g',2,ierr)
      Call EXPVAN (MyEXO%exoid,'g',2,(/'gRes1','gRes2'/),ierr)
      Do i = 1,10
         T = 1.0_Kr + i
         Call EXPTIM (MyEXO%exoid,i,T,ierr)
      End Do
      Call EXCLOS(MyEXO%exoid,ierr)
      
      !!! TEST GLOBAL IO
      T = 1.41_Kr + MEF90_MyRank
      Call Write_EXO_Result_Global(MyExo,2,1,T)
      T = 3.14_Kr + MEF90_MyRank
      Call Write_EXO_Result_Global(MyExo,1,1,T)

      Call Read_EXO_Result_Global(MyExo,2,1,T)

      Write(IOBuffer,*) 'Read ',T,'from EXO file ',Trim(MyEXO%filename),'\n'
      Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Call PetscSynchronizedFlush(PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)

      !!! TEST NODAL IO
      Call Write_EXO_Result_Vertex(MyExo,MeshTopology,2,1,Coords_Sec)
      Call SectionRealCreateLocalVector(Coords_Sec,Coords_VecL,ierr);CHKERRQ(ierr)
      Call Write_EXO_Result_Vertex(MyExo,MeshTopology,2,2,Coords_VecL)
      Call VecDestroy(Coords_VecL,ierr);CHKERRQ(ierr)
      Call Read_EXO_Result_Vertex(MyEXO,MeshTopology,2,1,U_Sec)
      Call SectionRealView(U_Sec,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
      
      Allocate(V2D(MeshTopology%Num_Verts))
      V2D(:)%X = -1.0_Kr
      V2D(:)%Y = -3.14_Kr
      Call Read_EXO_Result_Vertex(MyEXO,MeshTopology,2,1,V2D)
      If (verbose) Then
         Do i = 1,MeshTopology%Num_Verts
            Write(IOBuffer,*) i,V2D(i)%X,V2D(i)%Y,'\n'
            Call PetscViewerASCIIPrintf(myviewer,IOBuffer,ierr);CHKERRQ(ierr)
         End Do
      End If
      
      
      !!! TEST ELEMENTAL IO
      Allocate(V_Ptr(MeshTopology%Num_Elems))
      V_Ptr = MEF90_MyRank + 1.0_Kr
      Call Write_EXO_Result_Cell(MyExo,MeshTopology,2,2,V_Ptr)
      V_Ptr = 0.0_Kr
      Call Read_EXO_Result_Cell(MyExo,MeshTopology,2,2,V_Ptr)      
      Write(IOBuffer,*) MEF90_MyRank,'V_Ptr Min/Max:',MinVal(V_Ptr),MaxVal(V_Ptr),'\n'
      Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Call PetscSynchronizedFlush(PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)
      
      DeAllocate(V_Ptr)

      
   Call MeshTopologyDestroy(MeshTopology)
   Call SectionRealDestroy(Coords_Sec,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(U_Sec,ierr)
   DeAllocate(V2D)
   Call VecScatterDestroy(scatter,ierr);CHKERRQ(ierr)
   
   
   If (verbose) Then
      Call PetscViewerFlush(myviewer,ierr);CHKERRQ(ierr)
      Call PetscViewerDestroy(myviewer,ierr);CHKERRQ(ierr)
      Call PetscViewerFlush(viewer,ierr);CHKERRQ(ierr)
      Call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
   End If
   
300 Format(A)   
   Call MEF90_Finalize()
End Program TestLocal
