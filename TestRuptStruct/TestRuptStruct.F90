Program TestRuptStruct

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscmeshdef.h"

   Use m_MEF90
   Use m_RuptStruct

   Use petsc
   Use petscvec
   Use petscmat
   Use petscksp
   Use petscmesh

   Implicit NONE   
   
   Character(len=MEF90_MXSTRLEN)                :: prefix, IOBuffer, filename
   Type(MeshTopology_Type)                      :: MeshTopology
   Type(EXO_Type)                               :: EXO, MyEXO
   Type(Mesh)                                   :: Tmp_Mesh
   PetscTruth                                   :: HasPrefix, verbose
   PetscInt                                     :: iErr, i
   Type(PetscViewer)                            :: viewer, myviewer
   Type(SectionInt)                             :: BCFlag
   
   Call MEF90_Initialize()
   Call PetscOptionsHasName(PETSC_NULL_CHARACTER, '-verbose', verbose, iErr)    
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p',       prefix, HasPrefix, iErr); CHKERRQ(iErr)
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No input file prefix given\n", iErr)
      Call MEF90_Finalize()
      STOP
   End If

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


   Call Write_EXO_Case(prefix, '%0.4d', MEF90_NumProcs)

   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'
   
   Call EXO_Check_Numbering(EXO, iErr)
   If (iErr /= 0) Then
      SETERRQ(PETSC_ERR_SUP, 'Unsupported numbering of the element blocks, side sets or node sets\n'c, iErr)
   End If

   Call MeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
   If (verbose) Then
      Write(IOBuffer, '(A)') 'Done reading Sequential EXO file\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If

   Call MeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, MeshTopology%mesh, ierr); CHKERRQ(iErr)
   If (verbose) Then
      Write(IOBuffer, '(A)') 'Done distributing mesh\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   Call MeshDestroy(Tmp_mesh, iErr); CHKERRQ(iErr)

   Call MeshTopologyReadEXO(MeshTopology, EXO)
   !!! Sets the type of elements for each block
   Do i = 1, MeshTopology%Num_Elem_Blks
      MeshTopology%Elem_Blk(i)%Elem_Type = MEF90_P1_Lagrange
      Call Init_Elem_Blk_Type(MeshTopology%Elem_Blk(i), MeshTopology%num_dim)
   End Do
   If (verbose) Then
      Write(IOBuffer, '(A)') 'Done Initializing EXO data structure\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If


      
   MyEXO%comm = PETSC_COMM_SELF
   MyEXO%exoid = EXO%exoid
   Write(MyEXO%filename, 200) trim(prefix), MEF90_MyRank
   200 Format(A, '-', I4.4, '.gen')
   Call Write_MeshTopologyGlobal(MeshTopology, MyEXO, PETSC_COMM_WORLD)
   
   
   Call RuptEXOProperty_Init(MyEXO)   
   If (verbose) Then
      Write(IOBuffer, '(A)') 'Done with RuptEXOProperty_Init\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
!!!   Do i = 1, MeshTopology%Num_Elem_Blks
!!!      MyEXO%EBProperty(Rupt_EBProp_BCTypeX)%Value(i) = 100*i
!!!      MyEXO%EBProperty(Rupt_EBProp_BCTypeY)%Value(i) = 100*i
!!!      MyEXO%EBProperty(Rupt_EBProp_BCTypeZ)%Value(i) = 100*i
!!!   End Do
!!!   
!!!!   Do i = 1, size(
!!!      MyEXO%NSProperty(Rupt_NSProp_BCTypeX)%Value(:) = 1
!!!      MyEXO%NSProperty(Rupt_NSProp_BCTypeY)%Value(:) = 1
!!!      MyEXO%NSProperty(Rupt_NSProp_BCTypeZ)%Value(:) = 1
!!!!   End Do

   Call RuptEXOVariable_Init(MyEXO)   
   If (verbose) Then
      Write(IOBuffer, '(A)') 'Done with RuptEXOVariable_Init\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If

   If (verbose) Then
      Write(IOBuffer, '(A)') '\n\nMeshTopology\n'c
      Call PetscViewerASCIIPrintf(myviewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, '(A)') 'Wrote MeshTopology\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call MeshTopologyView(MeshTopology, myviewer); CHKERRQ(iErr)

      Write(IOBuffer, '(A)') '\n\nEXO\n'c
      Call PetscViewerASCIIPrintf(myviewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, '(A)') 'Wrote EXO\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call EXOView(EXO, myviewer)

      Write(IOBuffer, '(A)') '\n\nMyEXO\n'c
      Call PetscViewerASCIIPrintf(myviewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, '(A)') 'Wrote MyEXO\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call EXOView(MyEXO, myviewer)

   End If

   Call EXO_Variable_Write(MyEXO)
   If (verbose) Then
      Write(IOBuffer, '(A)') 'Done with EXO_Variable_Write\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   Call EXO_Property_Write(MyEXO)
   If (verbose) Then
      Write(IOBuffer, '(A)') 'Done with EXO_Property_Write\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   

   Call MeshGetVertexSectionInt(MeshTopology%mesh, 1, BCFlag, iErr); CHKERRQ(iErr)
   Call SectionIntZero(BCFlag, iErr); CHKERRQ(iErr)
   Call EXOProperty_InitBCFlag2DA(MyEXO, MeshTopology, BCFlag)
   Call SectionIntView(BCFlag, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)



   !!! Done playing, cleaning up
   Call MeshTopologyDestroy(MeshTopology)
   If (verbose) Then
      Call PetscViewerFlush(myviewer, iErr); CHKERRQ(iErr)
      Call PetscViewerDestroy(myviewer, iErr); CHKERRQ(iErr)
      Call PetscViewerFlush(viewer, iErr); CHKERRQ(iErr)
      Call PetscViewerDestroy(viewer, iErr); CHKERRQ(iErr)
   End If

   Call MEF90_Finalize()

End Program  TestRuptStruct
