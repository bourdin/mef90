Program TestLocal

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscmeshdef.h"

   Use m_MEF90
   Use petsc
   Use petscvec
   Use petscmesh

   Implicit NONE   

   Type (MeshTopology_Type)                     :: MeshTopology
   Type (EXO_Type)                              :: EXO, MyEXO
   Type (Element2D_Scal), Dimension(:), Pointer :: Elem2DA
   Type (Vect3D), Dimension(:), Pointer         :: Coords
   
   Type(Mesh)                                   :: Tmp_Mesh
   PetscBool                                   :: HasPrefix
   PetscBool                                   :: verbose
   PetscErrorCode                               :: iErr
   PetscInt                                     :: iBlk
   Character(len=256)                           :: CharBuffer, IOBuffer, filename
   Character(len=256)                           :: prefix
   Type(PetscViewer)                            :: viewer, myviewer
   Type(SectionReal)                            :: Coords_Sec, U_Sec
   Type(Vec)                                    :: Coords_VecG, Coords_VecL
   PetscReal, Dimension(:), Pointer             :: Coords_Ptr
   PetscInt                                     :: VSize
   Type(VecScatter)                             :: scatter
   PetscInt                                     :: i
   PetscReal                                    :: T
   Type(Vect2D), Dimension(:), Pointer          :: V2D
   Type(Vect3D), Dimension(:), Pointer          :: V3D
   PetscReal, Dimension(:), Pointer             :: V_Ptr
     
   Call MEF90_Initialize()
   Call PetscOptionsHasName(PETSC_NULL_CHARACTER, '-verbose', verbose, iErr)    
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p', prefix, HasPrefix, iErr)    
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No input file prefix given\n", iErr)
      Call MEF90_Finalize()
      STOP
   End If

   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'


   Call MeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
   Call MeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, MeshTopology%mesh, ierr); CHKERRQ(iErr)
   Call MeshDestroy(Tmp_mesh, iErr); CHKERRQ(iErr)
   Call MeshTopologyGetInfo(MeshTopology, PETSC_COMM_WORLD)
   Call EXOProperty_Read(EXO)
   Call EXOVariable_Read(EXO)
   
   MeshTopology%Elem_Blk%Elem_Type    = MEF90_P1_Lagrange
   Do iBlk = 1, MeshTopology%Num_Elem_Blks
      Call Init_Elem_Blk_Type(MeshTopology%Elem_Blk(iBlk), MeshTopology%num_dim)
   End Do
   
   If (verbose) Then
      Write(filename, 102) Trim(prefix), MEF90_MyRank
      Call PetscViewerASCIIOpen(PETSC_COMM_SELF, filename, myviewer, iErr); CHKERRQ(iErr);   
      Write(IOBuffer, 103) MEF90_MyRank, Trim(filename)
      Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call PetscSynchronizedFlush(PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)

      Write(filename, 101) Trim(prefix)
      Call PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, viewer, iErr); CHKERRQ(iErr);   
      Write(IOBuffer, 104) Trim(filename)
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
101 Format(A,'.log')
102 Format(A, '-', I4.4, '.log')
103 Format('Output from processor ', I4.4, ' redirected to ', A, '\n')
104 Format('Collective output redirected to ', A, '\n')
   MyEXO%comm = PETSC_COMM_SELF
   MyEXO%exoid = EXO%exoid
   Write(MyEXO%filename, 200) trim(prefix), MEF90_MyRank
200 Format(A, '-', I4.4, '.gen')
 
   MyEXO%title = trim(EXO%title)
   MyEXO%Num_QA = EXO%Num_QA
   Call EXOProperty_Copy(EXO, MyEXO)
   Call EXOVariable_Copy(EXO, MyEXO)

   Call Write_MeshTopologyGlobal(MeshTopology, MyEXO, PETSC_COMM_WORLD)
   Call EXOProperty_Write(MyEXO)
   Call EXOVariable_Write(MyEXO)

   If (verbose) Then
      Write(IOBuffer, 300) 'EXO\n'
      Call PetscViewerASCIIPrintf(myviewer, IOBuffer, iErr); CHKERRQ(iErr)
      Call EXOView(EXO, myviewer)
      Write(IOBuffer, 300) 'MyEXO\n'
      Call PetscViewerASCIIPrintf(myviewer, IOBuffer, iErr); CHKERRQ(iErr)
      Call EXOView(MyEXO, myviewer)

      Write(IOBuffer, 300) '\n\nMeshTopology\n'
      Call PetscViewerASCIIPrintf(myviewer, IOBuffer, iErr); CHKERRQ(iErr)
      Call MeshTopologyView(MeshTopology, myviewer); CHKERRQ(iErr)
   End If
   
      Call MEF90_Finalize()
      STOP
   
   !!! ****************** PLAYING WITH SECTIONS AND VECS ******************
   Call MeshGetSectionReal(MeshTopology%mesh, 'coordinates', Coords_Sec, iErr); CHKERRQ(ierr) 
   If (verbose) Then
      Write(IOBuffer, 300) '\n\nSECTION STUFF\n'
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Call PetscViewerASCIIPrintf(myviewer, IOBuffer, iErr); CHKERRQ(iErr)

      Call SectionRealView(Coords_Sec, viewer, iErr); CHKERRQ(iErr)
   End If


   Write(IOBuffer, 400) MeshTopology%Num_Verts, MeshTopology%Num_Elems
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
   
   If (verbose) Then
      Call PetscViewerASCIIPrintf(myviewer, IOBuffer, iErr); CHKERRQ(iErr)
   End If
400 Format('Number of Vertices: ', I3, ' Number of elements: ', I3, '\n')


   Call MeshCreateVector(MeshTopology%mesh, Coords_Sec, Coords_VecG, iErr); CHKERRQ(iErr)

   Write(IOBuffer, *) '\n\nVec obtained with MeshCreateVector have no \"ghost values\" \n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

   Call VecGetSize(Coords_VecG, VSize, iErr); CHKERRQ(iErr)
   Write(IOBuffer, 401) VSize
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

   Call VecGetLocalSize(Coords_VecG, VSize, iErr); CHKERRQ(iErr)
   Write(IOBuffer, 402) VSize
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
401 Format('   Global size of Coords_VecG: ', I3, '\n')
402 Format('   Local size of Coords_VecG: ', I3, '\n')

   If (verbose) Then
      !!! Scatter the values
      Write(IOBuffer, *) 'Coords_VecG before scattering\n'
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Call VecView(Coords_VecG, viewer, iErr); CHKERRQ(iErr)
      
      Call MeshCreateGlobalScatter(MeshTopology%mesh, Coords_Sec, scatter, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(Coords_Sec, scatter, SCATTER_FORWARD, Coords_VecG, ierr); CHKERRQ(ierr)
      Write(IOBuffer, *) 'Coords_VecG after scattering\n'
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Call VecView(Coords_VecG, viewer, iErr); CHKERRQ(iErr)
   End If

   Call SectionRealCreateLocalVector(Coords_Sec, Coords_VecL, iErr); CHKERRQ(iErr)
   If (verbose) Then
      !!! Scatter the values
      Write(IOBuffer, *) 'Coords_VecL\n'
      Call PetscViewerASCIIPrintf(myviewer, IOBuffer, iErr); CHKERRQ(iErr)
      Call VecView(Coords_VecL, myviewer, iErr); CHKERRQ(iErr)
   End If


   Write(IOBuffer, *) '\n\nVec obtained with SectionRealCreateLocalVector do have \"ghost values\" \n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

   Call VecGetLocalSize(Coords_VecL, VSize, iErr); CHKERRQ(iErr)
   Write(IOBuffer, 404) VSize
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
403 Format('   Global size of Coords_VecL: ', I3, '\n')
404 Format('   Local size of Coords_VecL: ', I3, '\n')

   Call VecDestroy(Coords_VecG, iErr); CHKERRQ(iErr)
   Call VecDestroy(Coords_VecL, iErr); CHKERRQ(iErr)

   Call MeshGetVertexSectionReal(MeshTopology%mesh, 2, U_Sec, iErr); CHKERRQ(iErr)

!!! Format the output mesh with variables
      MyEXO%exoid = EXOPEN(MyEXO%filename, EXWRIT, exo_cpu_ws, exo_io_ws, exo_ver, ierr)
      
      Call EXPVP (MyEXO%exoid, 'n', 4, iErr)
      Call EXPVAN (MyEXO%exoid, 'n', 4, (/'nRes1', 'nRes2', 'nRes3', 'nRes4'/), iErr)
      Call EXPVP (MyEXO%exoid, 'e', 4, iErr)
      Call EXPVAN (MyEXO%exoid, 'e', 4, (/'eRes1', 'eRes2', 'eRes3', 'eRes4'/), iErr)
      Call EXPVP (MyEXO%exoid, 'g', 2, iErr)
      Call EXPVAN (MyEXO%exoid, 'g', 2, (/'gRes1', 'gRes2'/), iErr)
      Do i = 1, 10
         T = 1.0_Kr + i
         Call EXPTIM (MyEXO%exoid, i, T, iErr)
      End Do
      Call EXCLOS(MyEXO%exoid, iErr)
      
      !!! TEST GLOBAL IO
      T = 1.41_Kr + MEF90_MyRank
      Call Write_EXO_Result_Global(MyExo, 2, 1, T)
      T = 3.14_Kr + MEF90_MyRank
      Call Write_EXO_Result_Global(MyExo, 1, 1, T)

      Call Read_EXO_Result_Global(MyExo, 2, 1, T)

      Write(IOBuffer, *) 'Read ', T, 'from EXO file ', Trim(MyEXO%filename), '\n'
      Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call PetscSynchronizedFlush(PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)

      !!! TEST NODAL IO
      Call Write_EXO_Result_Vertex(MyExo, MeshTopology, 2, 1, Coords_Sec)
      Call SectionRealCreateLocalVector(Coords_Sec, Coords_VecL, iErr); CHKERRQ(iErr)
      Call Write_EXO_Result_Vertex(MyExo, MeshTopology, 2, 2, Coords_VecL)
      Call VecDestroy(Coords_VecL, iErr); CHKERRQ(iErr)
      Call Read_EXO_Result_Vertex(MyEXO, MeshTopology, 2, 1, U_Sec)
      Call SectionRealView(U_Sec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
      
      Allocate(V2D(MeshTopology%Num_Verts))
      V2D(:)%X = -1.0_Kr
      V2D(:)%Y = -3.14_Kr
      Call Read_EXO_Result_Vertex(MyEXO, MeshTopology, 2, 1, V2D)
      If (verbose) Then
         Do i = 1, MeshTopology%Num_Verts
            Write(IOBuffer, *) i, V2D(i)%X, V2D(i)%Y, '\n'
            Call PetscViewerASCIIPrintf(myviewer, IOBuffer, iErr); CHKERRQ(iErr)
         End Do
      End If
      
      
      !!! TEST ELEMENTAL IO
      Allocate(V_Ptr(MeshTopology%Num_Elems))
      V_Ptr = MEF90_MyRank + 1.0_Kr
      Call Write_EXO_Result_Cell(MyExo, MeshTopology, 2, 2, V_Ptr)
      V_Ptr = 0.0_Kr
      Call Read_EXO_Result_Cell(MyExo, MeshTopology, 2, 2, V_Ptr)      
      Write(IOBuffer, *) MEF90_MyRank, 'V_Ptr Min/Max:', MinVal(V_Ptr), MaxVal(V_Ptr), '\n'
      Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call PetscSynchronizedFlush(PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      
      DeAllocate(V_Ptr)

      
   Call MeshTopologyDestroy(MeshTopology)
   Call SectionRealDestroy(Coords_Sec, iErr); CHKERRQ(iErr)
   Call SectionRealDestroy(U_Sec, iErr)
   DeAllocate(V2D)
   Call VecScatterDestroy(scatter, iErr); CHKERRQ(iErr)
   
   
   If (verbose) Then
      Call PetscViewerFlush(myviewer, iErr); CHKERRQ(iErr)
      Call PetscViewerDestroy(myviewer, iErr); CHKERRQ(iErr)
      Call PetscViewerFlush(viewer, iErr); CHKERRQ(iErr)
      Call PetscViewerDestroy(viewer, iErr); CHKERRQ(iErr)
   End If
   
300 Format(A)   
   Call MEF90_Finalize()
End Program TestLocal
