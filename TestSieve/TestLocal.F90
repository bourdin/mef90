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

   Type (MeshTopology_Info)                     :: MeshTopology
   Type (EXO_Info)                              :: EXO, MyEXO
   Type (Element2D_Scal), Dimension(:), Pointer :: Elem2DA
   Type (Vect3D), Dimension(:), Pointer         :: Coords
   
   PetscTruth                                   :: HasPrefix
   PetscTruth                                   :: verbose
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
   PetscInt                                     :: exo_ver
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


!   Call MeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, GlobalMeshTopology%mesh, ierr)

   Call MeshTopologyReadEXO(MeshTopology, Coords, Elem2DA, EXO)
   
   MeshTopology%Elem_Blk%Elem_Type    = MEF90_P1_Lagrange
   Do iBlk = 1, MeshTopology%Num_Elem_Blks
      Call Init_Elem_Blk_Info(MeshTopology%Elem_Blk(iBlk), MeshTopology%num_dim)
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
103 Format('Output from processor ', I4.4, ' redirected to ', A, '\n'c)
104 Format('Collective output redirected to ', A, '\n'c)
   MyEXO%comm = PETSC_COMM_SELF
   MyEXO%exoid = EXO%exoid
   Write(MyEXO%filename, 200) trim(prefix), MEF90_MyRank
200 Format(A, '-', I4.4, '.gen')
 
   MyEXO%title = trim(EXO%title)
   MyEXO%Num_QA = EXO%Num_QA
   

   Call Write_MeshTopologyGlobal(MeshTopology, MyEXO, PETSC_COMM_WORLD)

   If (verbose) Then
      Write(IOBuffer, 300) 'EXO_Info\n'c
      Call PetscViewerASCIIPrintf(myviewer, IOBuffer, iErr); CHKERRQ(iErr)
      Call EXOView(EXO, myviewer)

      Write(IOBuffer, 300) '\n\nMeshTopology\n'c
      Call PetscViewerASCIIPrintf(myviewer, IOBuffer, iErr); CHKERRQ(iErr)
      Call MeshTopologyView(MeshTopology, myviewer); CHKERRQ(iErr)
   End If
   
   !!! ****************** PLAYING WITH SECTIONS AND VECS ******************
   Call MeshGetSectionReal(MeshTopology%mesh, 'coordinates', Coords_Sec, iErr); CHKERRQ(ierr) 
   If (verbose) Then
      Write(IOBuffer, 300) '\n\nSECTION STUFF\n'c
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
400 Format('Number of Vertices: ', I3, ' Number of elements: ', I3, '\n'c)


   Call MeshCreateVector(MeshTopology%mesh, Coords_Sec, Coords_VecG, iErr); CHKERRQ(iErr)

   Write(IOBuffer, *) '\n\nVec obtained with MeshCreateVector have no \"ghost values\" \n'c
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

   Call VecGetSize(Coords_VecG, VSize, iErr); CHKERRQ(iErr)
   Write(IOBuffer, 401) VSize
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

   Call VecGetLocalSize(Coords_VecG, VSize, iErr); CHKERRQ(iErr)
   Write(IOBuffer, 402) VSize
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
401 Format('   Global size of Coords_VecG: ', I3, '\n'c)
402 Format('   Local size of Coords_VecG: ', I3, '\n'c)

   If (verbose) Then
      !!! Scatter the values
      Write(IOBuffer, *) 'Coords_VecG before scattering\n'c
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Call VecView(Coords_VecG, viewer, iErr); CHKERRQ(iErr)
      
      Call MeshCreateGlobalScatter(MeshTopology%mesh, Coords_Sec, scatter, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(Coords_Sec, scatter, SCATTER_FORWARD, Coords_VecG, ierr); CHKERRQ(ierr)
      Write(IOBuffer, *) 'Coords_VecG after scattering\n'c
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Call VecView(Coords_VecG, viewer, iErr); CHKERRQ(iErr)
   End If

   Call SectionRealCreateLocalVector(Coords_Sec, Coords_VecL, iErr); CHKERRQ(iErr)
   If (verbose) Then
      !!! Scatter the values
      Write(IOBuffer, *) 'Coords_VecL\n'c
      Call PetscViewerASCIIPrintf(myviewer, IOBuffer, iErr); CHKERRQ(iErr)
      Call VecView(Coords_VecL, myviewer, iErr); CHKERRQ(iErr)
   End If


   Write(IOBuffer, *) '\n\nVec obtained with SectionRealCreateLocalVector do have \"ghost values\" \n'c
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

   Call VecGetLocalSize(Coords_VecL, VSize, iErr); CHKERRQ(iErr)
   Write(IOBuffer, 404) VSize
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
403 Format('   Global size of Coords_VecL: ', I3, '\n'c)
404 Format('   Local size of Coords_VecL: ', I3, '\n'c)

   Call VecDestroy(Coords_VecG, iErr); CHKERRQ(iErr)
   Call VecDestroy(Coords_VecL, iErr); CHKERRQ(iErr)

   Call MeshGetVertexSectionReal(MeshTopology%mesh, 2, U_Sec, iErr); CHKERRQ(iErr)

!      Write(*,*) 'Writing in ', trim(MyEXO%filename)
      MyEXO%exoid = EXOPEN(MyEXO%filename, EXWRIT, exo_cpu_ws, exo_io_ws, exo_ver, ierr)
      
      Call EXPVP (MyEXO%exoid, 'n', 4, iErr)
      Call EXPVAN (MyEXO%exoid, 'n', 4, (/'res1', 'Res2', 'Res3', 'Res4'/), iErr)
      Do i = 1, 10
         T = 1.0_Kr + i
         Call EXPTIM (MyEXO%exoid, i, T, iErr)
      End Do
      Call EXCLOS(MyEXO%exoid, iErr)

      Call Write_EXO_Result_Vertex(MyExo, MeshTopology, 2, 1, Coords_Sec)
      Call Read_EXO_Result_Vertex(MyEXO, MeshTopology, 2, 1, U_Sec)
      Call SectionRealView(U_Sec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
      
      Allocate(V2D(MeshTopology%Num_Verts))
      V2D(:)%X = -1.0_Kr
      V2D(:)%Y = -3.14_Kr
      Call Read_EXO_Result_Vertex(MyEXO, MeshTopology, 2, 1, V2D)
      If (verbose) Then
         Do i = 1, MeshTopology%Num_Verts
            Write(IOBuffer, *) i, V2D(i)%X, V2D(i)%Y, '\n'c
            Call PetscViewerASCIIPrintf(myviewer, IOBuffer, iErr); CHKERRQ(iErr)
         End Do
      End If
   Call MeshTopologyDestroy(MeshTopology)
   Call SectionRealDestroy(Coords_Sec, iErr); CHKERRQ(iErr)
   If (verbose) Then
      Call PetscViewerFlush(myviewer, iErr); CHKERRQ(iErr)
      Call PetscViewerDestroy(myviewer, iErr); CHKERRQ(iErr)
      Call PetscViewerFlush(viewer, iErr); CHKERRQ(iErr)
      Call PetscViewerDestroy(viewer, iErr); CHKERRQ(iErr)
      Call VecScatterDestroy(scatter, iErr); CHKERRQ(iErr)
   End If
   
300 Format(A)   
   Call MEF90_Finalize()
End Program TestLocal
