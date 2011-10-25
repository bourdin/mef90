#if defined PB_2D
Module m_TransientHeat2D
#elif defined PB_3D
Module m_TransientHeat3D
#endif

#include "finclude/petscdef.h"

   Use m_MEF90

#if defined PB_2D
   Use m_Poisson2D
#elif defined PB_3D 
   Use m_Poisson3D
#endif   

   Implicit NONE   

   
Type TestCase_Type
   PetscInt                                  :: Index
   Character(len=MEF90_MXSTRLEN)             :: Description
End Type

Contains

    Subroutine AskInt(val, msg, ArgUnit, IsBatch)
      PetscInt                                  :: Val
      Character(len=*)                          :: msg 
      PetscInt                                  :: argunit
      PetscBool                                 :: IsBatch

      Character(len=MEF90_MXSTRLEN)             :: prefix, IOBuffer   
      PetscInt                                  :: iErr   
      If (IsBatch) Then
         If (MEF90_MyRank == 0) Then
            Read(ArgUnit,*) Val
         End If
         Call MPI_BCast(Val, 1, MPIU_INTEGER, 0, PETSC_COMM_WORLD, iErr)
      Else
         Write(IOBuffer, "(A, t60,':  ')") Trim(msg)
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         If (MEF90_MyRank == 0) Then
            Read(*,*) Val
            Write(ArgUnit, "(I4, t60, A)") val, Trim(msg)
         End If
         Call MPI_BCast(Val, 1, MPIU_INTEGER, 0, PETSC_COMM_WORLD, iErr)
      End If
   End Subroutine AskInt   

  Subroutine AskReal(val, msg, ArgUnit, IsBatch)
      PetscReal                                 :: Val
      Character(len=*)                          :: msg 
      PetscInt                                  :: argunit
      PetscBool                                 :: IsBatch

      Character(len=MEF90_MXSTRLEN)             :: prefix, IOBuffer      
      PetscInt                                  :: iErr
      If (IsBatch) Then
         If (MEF90_MyRank == 0) Then
            Read(ArgUnit,*) Val
         End If
         Call MPI_BCast(Val, 1, MPIU_SCALAR, 0, PETSC_COMM_WORLD, iErr)
      Else
         Write(IOBuffer, "(A, t60,':  ')") Trim(msg)
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         If (MEF90_MyRank == 0) Then
            Read(*,*) Val
            Write(ArgUnit, "(ES12.5, t60, A)") val, Trim(msg)
         End If
         Call MPI_BCast(Val, 1, MPIU_SCALAR, 0, PETSC_COMM_WORLD, iErr)
      End If
   End Subroutine AskReal

 Subroutine EXONSProperty_AskWithBatch(dEXO, dMeshTopology, BatchUnit, IsBatch)
      Type(EXO_Type)                                 :: dEXO
      Type(MeshTopology_Type)                        :: dMeshTopology
      PetscInt                                       :: BatchUnit
      PetscBool                                      :: IsBatch

      PetscInt                                       :: iErr
      PetscInt                                       :: i, j, IntBuffer

      PetscInt                                       :: NumNS
      PetscInt                                       :: EXO_MyRank
      Character(len=MEF90_MXSTRLEN)                  :: IOBuffer

      Do i = 1, dMeshTopology%Num_Node_Sets_Global
         Write(IOBuffer, 102) i
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Do j = 1, dEXO%Num_NSProperties
            Write(IOBuffer, 202) i, Trim(dEXO%NSProperty(j)%Name)
            Call AskInt(dEXO%NSProperty(j)%Value(i), IOBuffer, BatchUnit, IsBatch)
         End Do
         If ((.NOT. IsBatch) .AND. (MEF90_MyRank == 0)) Then
            Write(BatchUnit, *)
         End If
      End Do
102 Format('    Node Set      ', T24, I3, '\n')
202 Format('NS', I4.4, ': ', A)
   End Subroutine EXONSProperty_AskWithBatch


   Subroutine EXOEBProperty_AskWithBatch(dEXO, dMeshTopology, BatchUnit, IsBatch)
      Type(EXO_Type)                                 :: dEXO
      Type(MeshTopology_Type)                        :: dMeshTopology
      PetscInt                                       :: BatchUnit
      PetscBool                                      :: IsBatch

      PetscInt                                       :: iErr
      PetscInt                                       :: i, j, IntBuffer

      PetscInt                                       :: NumEB
      PetscInt                                       :: EXO_MyRank
      Character(len=MEF90_MXSTRLEN)                  :: IOBuffer
      PetscReal                                      :: TmpEBProperty

      Do i = 1, dMeshTopology%Num_Elem_Blks_Global
         Write(IOBuffer, 100) i
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Do j = 1, dEXO%Num_EBProperties
            Write(IOBuffer, 200) i, Trim(dEXO%EBProperty(j)%Name)
            Call AskInt(dEXO%EBProperty(j)%Value(i), IOBuffer, BatchUnit, IsBatch)
         End Do
         If ((.NOT. IsBatch) .AND. (MEF90_MyRank == 0)) Then
            Write(BatchUnit, *)
         End If
      End Do
100 Format('    Element Block ', T24, I3, '\n')
200 Format('EB', I4.4, ': ', A)
   End Subroutine EXOEBProperty_AskWithBatch
   
   Subroutine EXOSSProperty_AskWithBatch(dEXO, dMeshTopology, BatchUnit, IsBatch)
      Type(EXO_Type)                                 :: dEXO
      Type(MeshTopology_Type)                        :: dMeshTopology
      PetscInt                                       :: BatchUnit
      PetscBool                                      :: IsBatch

      PetscInt                                       :: iErr
      PetscInt                                       :: i, j, IntBuffer

      PetscInt                                       :: NumSS
      PetscInt                                       :: EXO_MyRank
      Character(len=MEF90_MXSTRLEN)                  :: IOBuffer

      Do i = 1, dMeshTopology%Num_Side_Sets_Global
         Write(IOBuffer, 101) i
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Do j = 1, dEXO%Num_SSProperties
            Write(IOBuffer, 201) i, Trim(dEXO%SSProperty(j)%Name)
            Call AskInt(dEXO%SSProperty(j)%Value(i), IOBuffer, BatchUnit, IsBatch)
         End Do
         If ((.NOT. IsBatch) .AND. (MEF90_MyRank == 0)) Then
            Write(BatchUnit, *)
         End If
      End Do
101 Format('    Side Set      ', T24, I3, '\n')
201 Format('SS', I4.4, ': ', A)
   End Subroutine EXOSSProperty_AskWithBatch


#undef __FUNCT__
#define __FUNCT__ "PoissonInit"
   Subroutine PoissonInit(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr
      PetscBool                                    :: HasPrefix, Flag
      PetscInt                                     :: iBlk, iDoF, i 
      Character(len=MEF90_MXSTRLEN)                :: BatchFileName
      PetscInt                                     :: BatchUnit=99
      PetscInt                                     :: NumTestCase
      Type(TestCase_Type), Dimension(:) , Pointer  :: TestCase
      PetscBool                                    :: EraseBatch
      PetscBool                                    :: IsBatch, HasBatchFile
      Character(len=MEF90_MXSTRLEN)                :: prefix, IOBuffer, filename   
      Type(DM)                                     :: Tmp_Mesh
      PetscReal, Dimension(:), Pointer             :: TmpCoords, ValPtr
      PetscReal, Dimension(:,:), Pointer           :: Coords
      PetscReal                                    :: ValU, ValF
      PetscInt, Dimension(:), Pointer              :: SizeVect, SizeScal
      PetscInt, Dimension(:), Pointer              :: TmpFlag
      PetscInt                                     :: TmpPoint
      PetscInt, Parameter                          :: VarFrac_NSProp_BCT =1
      PetscReal, Dimension(:), Pointer             :: U, Uelem


      Call MEF90_Initialize()
      
      NumTestCase = 2
      Allocate(TestCase(NumTestCase))
      Do i = 1, NumTestCase
         TestCase(i)%Index = i
      End Do
      TestCase(1)%Description = "Simple Poisson \Delta u = f"
      TestCase(2)%Description = "Heat equation u,t - \Delta u = f"
      
      AppCtx%AppParam%verbose = 0
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-verbose', AppCtx%AppParam%verbose, Flag, iErr); CHKERRQ(iErr)
      
      !AppCtx%AppParam%Restart = PETSC_FALSE
      !Call PetscOptionsGetBool(PETSC_NULL_CHARACTER, '-restart', AppCtx%AppParam%restart, Flag, iErr); CHKERRQ(iErr)
      
      Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p',       AppCtx%AppParam%prefix, HasPrefix, iErr); CHKERRQ(iErr)
      If (.NOT. HasPrefix) Then
         Call PetscPrintf(PETSC_COMM_WORLD, "No mesh prefix given\n", iErr)
         Call MEF90_Finalize()
         STOP
      End If
      Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-i', BatchFileName, IsBatch, iErr);CHKERRQ(iErr)

      EraseBatch=.False.
      Call PetscOptionsGetBool(PETSC_NULL_CHARACTER, '-force', EraseBatch, HasPrefix, iErr)    
      If (MEF90_MyRank==0) Then
         If (IsBatch) Then
            Write(IOBuffer, *) "\nProcessing batch file ", Trim(BatchFileName), "\n"
            Call PetscPrintf(PETSC_COMM_SELF, IOBuffer, iErr); CHKERRQ(iErr)
            Open(Unit=BatchUnit, File=BatchFileName, Status='Old', Action='Read')
            Rewind(BatchUnit)
         Else
            BatchFileName = Trim(prefix)//'.args'
            Inquire(File=BatchFileName, EXIST=HasBatchFile)
            If (HasBatchFile .AND. (.NOT. EraseBatch)) Then
               Write(IOBuffer, *) "Batch file ", trim(BatchFileName), " already exists. Erase it or use -force flag\n"
               Call PetscPrintf(PETSC_COMM_SELF, IOBuffer, iErr); CHKERRQ(iErr)
               SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_UNEXPECTED, IOBuffer, iErr)
            Else
               Write(IOBuffer, *) "Running interactively and generating batch file ", trim(BatchFileName), "\n"
               Call PetscPrintf(PETSC_COMM_SELF, IOBuffer, iErr); CHKERRQ(iErr)
               Open(Unit=BatchUnit, File=BatchFileName, Status='Unknown')
               Rewind(BatchUnit)
            End If
         End If
      End If

      AppCtx%AppParam%TestCase = 1
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-test',       AppCtx%AppParam%TestCase, Flag, iErr); CHKERRQ(iErr)
      Call InitLog(AppCtx)
      Call PetscLogStagePush(AppCtx%LogInfo%Setup_Stage, iErr); CHKERRQ(iErr)
      If (AppCtx%AppParam%verbose > 1) Then
         Write(filename, 101) Trim(AppCtx%AppParam%prefix), MEF90_MyRank
         Call PetscViewerASCIIOpen(PETSC_COMM_SELF, filename, AppCtx%AppParam%MyLogViewer, iErr); CHKERRQ(iErr);   
         Write(IOBuffer, 102) MEF90_MyRank, Trim(filename)
         Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Call PetscSynchronizedFlush (PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
   
         Write(filename, 103) Trim(AppCtx%AppParam%prefix)
         Call PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr);   
         Write(IOBuffer, 104) Trim(filename)
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
   
101 Format(A, '-', I4.4, '.log')
102 Format('Output from processor ', I4.4, ' redirected to file ', A, '\n')
103 Format(A,'.log')
104 Format('Collective output redirected to file ', A, '\n')

      AppCtx%EXO%Comm = PETSC_COMM_WORLD
      AppCtx%EXO%num_nsproperties = 0
      AppCtx%EXO%num_ssproperties = 0
      AppCtx%EXO%num_ebproperties = 0
      AppCtx%EXO%filename = Trim(AppCtx%AppParam%prefix)//'.gen'
      !!! Read and partition the mesh
      If (MEF90_NumProcs == 1) Then
         Call PetscLogStagePush(AppCtx%LogInfo%MeshCreateExodus_Stage, iErr); CHKERRQ(iErr)
         Call DMMeshCreateExodus(PETSC_COMM_WORLD, AppCtx%EXO%filename, AppCtx%MeshTopology%mesh, ierr); CHKERRQ(iErr)
         Call PetscLogStagePop(iErr); CHKERRQ(iErr)
      Else
         Call PetscLogStagePush(AppCtx%LogInfo%MeshCreateExodus_Stage, iErr); CHKERRQ(iErr)
         Call DMMeshCreateExodus(PETSC_COMM_WORLD, AppCtx%EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
         Call PetscLogStagePop(iErr); CHKERRQ(iErr)
      
         Call PetscLogStagePush(AppCtx%LogInfo%MeshDistribute_Stage, iErr); CHKERRQ(iErr)
         Call DMMeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, AppCtx%MeshTopology%mesh, ierr); CHKERRQ(iErr)
         Call DMDestroy(Tmp_mesh, ierr); CHKERRQ(iErr)
         Call PetscLogStagePop(iErr); CHKERRQ(iErr)
      End If

      Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
      Call MeshTopologyReadEXO(AppCtx%MeshTopology, AppCtx%EXO)

      !!! Sets the type of elements for each block
      Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_Type = MEF90_P1_Lagrange
         Call Init_Elem_Blk_Type(AppCtx%MeshTopology%Elem_Blk(iBlk), AppCtx%MeshTopology%num_dim)
      End Do
   
      Call ElementInit(AppCtx%MeshTopology, AppCtx%Elem, 2)
      !!! Allocate the Section for U and F
      Allocate(SizeScal(1))
      SizeScal=1

      Call FieldCreateVertex(AppCtx%U,     'U',         AppCtx%MeshTopology, SizeScal)
      Call FieldCreateVertex(AppCtx%F,     'F',         AppCtx%MeshTopology, SizeScal)
      Call FieldCreateVertex(AppCtx%RHS,     'RHS',       AppCtx%MeshTopology, SizeScal)

      
      !!! Allocate and initialize the Section for the flag
      Call FlagCreateVertex(AppCtx%BCFlag, 'BC', AppCtx%MeshTopology, SizeScal)
      DeAllocate(SizeScal)

      !Allocate(TmpFlag(1))
      !Do iBlk = 1, AppCtx%MeshTopology%num_node_sets  
      !   Do iDoF = 1, AppCtx%MeshTopology%node_set(iBlk)%Num_Nodes     
      !      TmpPoint = AppCtx%MeshTopology%Num_Elems + AppCtx%MeshTopology%Node_Set(iBlk)%Node_ID(iDoF) - 1
      !      Call SectionIntUpdate(AppCtx%BCFlag%Sec, TmpPoint, TmpFlag, INSERT_VALUES, iErr); CHKERRQ(iErr)
      !   End Do
      !End Do
      !DeAllocate(TmpFlag)
      
      AppCtx%MyEXO%comm = PETSC_COMM_SELF
      AppCtx%MyEXO%exoid = AppCtx%EXO%exoid
      Write(AppCtx%MyEXO%filename, 200) trim(AppCtx%AppParam%prefix), MEF90_MyRank
200 Format(A, '-', I4.4, '.gen')
      AppCtx%MyEXO%title = trim(AppCtx%EXO%title)
      AppCtx%MyEXO%Num_QA = AppCtx%EXO%Num_QA
      !If (AppCtx%AppParam%Restart) Then
      !   Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
      !   Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, 1, 1, AppCtx%U) 
      !   Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, 2, 1, AppCtx%F) 
      !   Call PetscLogStagePop(iErr); CHKERRQ(iErr)
      !Else
         !!! Prepare and format the output mesh   
      Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
      Call MeshTopologyWriteGlobal(AppCtx%MeshTopology, AppCtx%MyEXO, PETSC_COMM_WORLD)
      Call EXOFormat_SimplePoisson(AppCtx)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)

      Write(IOBuffer, *) '\nTest Case:\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)      
      Do i = 1, NumTestCase
         Write(IOBuffer, "('   [',I2.2,'] ',A)"), TestCase(i)%Index,   Trim(TestCase(i)%Description)//'\n'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr);   CHKERRQ(iErr)
      End Do




         
      Call AskInt(AppCtx%AppParam%TestCase, 'Test Case', BatchUnit, IsBatch)
!pk est ce que les SectionRealSet suivant fonctionnent ?????         
! TODO pb we are just reading the first and seconf lines of the files !!!! 
         !Setting initiale value
      Call AskReal(ValU, 'Initial value in U ', BatchUnit, IsBatch)
      Call SectionRealSet(AppCtx%U, ValU, iErr); CHKERRQ(iErr);
         !Setting force
      Call AskReal(ValF, 'RHS F', BatchUnit, IsBatch)
      Call SectionRealSet(AppCtx%F, ValF, iErr); CHKERRQ(iErr);
      print *, ValU
      print *, ValF

!         Select Case (AppCtx%AppParam%TestCase)
!         Case(1)
!            If (AppCtx%AppParam%verbose > 0) Then
!               Write(IOBuffer, *) 'Setting U to 0 and F to 1\n'
!               Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
!            End If
!            !!! U = 0, F=1
!            Val = 1.0_Kr
!            Call SectionRealSet(AppCtx%F, Val, iErr); CHKERRQ(iErr);
!            Val = 0.0_Kr
!            Call SectionRealSet(AppCtx%U, Val, iErr); CHKERRQ(iErr);
!         Case(2)
!            If (AppCtx%AppParam%verbose > 0) Then
!               Write(IOBuffer, *) 'Solving Test Case 2: pure dirichlet problem, no force\n'
!               Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
!            End If

!            !!! Test of non homogeneous Dirichlet BC
!            Allocate(ValPtr(1))
!            ValPtr = 0.0_Kr
!            Call SectionRealSet(AppCtx%F, 0.0_Kr, iErr); CHKERRQ(iErr);
!            Call DMMeshGetCoordinatesF90(AppCtx%MeshTopology%Mesh, Coords, iErr); CHKERRQ(iErr)
!               
!            Do iDoF = 1, Size(Coords,1)
!#if defined PB_2D
!               ValPtr = (Coords(iDoF,1)-0.5_Kr)**2 + (Coords(iDoF,2)+0.5_Kr)**2
!               ValPtr = Coords(iDoF,1)**2 + Coords(iDoF,2)**2
!               ValPtr = Coords(iDoF,1)**3
!#elif defined PB_3D
!               ValPtr = (Coords(iDoF,1)-0.5_Kr)**2 + (Coords(iDoF,2)+0.5_Kr)**2 +  (Coords(iDoF,3)+.5_Kr)**2
!#endif
!               Call SectionRealUpdate(AppCtx%U%Sec, AppCtx%MeshTopology%Num_Elems+iDoF-1, ValPtr, INSERT_VALUES, iErr); CHKERRQ(iErr)
!            End Do
!            DeAllocate(ValPtr)
!            Call DMMeshRestoreCoordinatesF90(AppCtx%MeshTopology%Mesh, Coords, iErr); CHKERRQ(iErr)
!         Case(3)
!            If (AppCtx%AppParam%verbose > 0) Then
!               Write(IOBuffer, *) 'Solving Test Case 3: F=sgn(x) . sgn(y) [. sgn(z)] \n'
!               Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
!            End If

            !!! Test of non homogeneous Dirichlet BC
!            Allocate(ValPtr(1))
!            ValPtr = 0.0_Kr
!            Call SectionRealSet(AppCtx%U, 1.0_Kr, iErr); CHKERRQ(iErr);
!            Call DMMeshGetCoordinatesF90(AppCtx%MeshTopology%Mesh, Coords, iErr); CHKERRQ(iErr)
!            
!            Do iDoF = 1, Size(Coords,1)
!#if defined PB_2D
!               If ( Coords(iDoF,1) * Coords(iDoF,2) < 0.0_Kr ) Then
!                  ValPtr = -1.0_Kr
!               Else
!                  ValPtr = 1.0_Kr
!               End If
!#elif defined PB_3D
!               If ( Coords(iDoF,1) * Coords(iDoF,2) * Coords(iDoF,3) < 0.0_Kr ) Then
!                  ValPtr = -1.0_Kr
!               Else
!                  ValPtr = 1.0_Kr
!               End If
!#endif
!               Call SectionRealUpdate(AppCtx%F%Sec, AppCtx%MeshTopology%Num_Elems+iDoF-1, ValPtr, INSERT_VALUES, iErr); CHKERRQ(iErr)
!            End Do
!            Call DMMeshRestoreCoordinatesF90(AppCtx%MeshTopology%Mesh, Coords, iErr); CHKERRQ(iErr)
!         End Select            
      Call SectionIntAddNSProperty(AppCtx%BCFlag%Sec,  AppCtx%MyEXO%NSProperty(VarFrac_NSProp_BCT),  AppCtx%MeshTopology)
!lines 535-550 PrepVarFracNG
      Allocate(U(AppCtx%MeshTopology%Num_Node_Sets)) 
      U = 0.0_Kr
      Do i = 1, AppCtx%MeshTopology%Num_Node_Sets_Global
         Write(IOBuffer, 202) i
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         If (AppCtx%MyEXO%NSProperty(VarFrac_NSProp_BCT)%Value(i) /= 0 ) Then
            Write(IOBuffer, 302) i, 'Ux'
            Call AskReal(U(i), IOBuffer, BatchUnit, IsBatch)
         End If
      End Do
      !End If

202 Format('    Node Set      ', T24, I3, '\n')
302 Format('NS', I4.4, ': ', A)


   End Subroutine PoissonInit

#undef __FUNCT__
#define __FUNCT__ "KSPSetUp"
   Subroutine KSPSetUp(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr
      PetscInt                                     :: iBlk, iDoF      
      PetscInt, Dimension(:), Pointer              :: TmpFlag
      PetscInt                                     :: iE, iELoc
      Character(len=MEF90_MXSTRLEN)                :: prefix, IOBuffer, filename   
      PetscLogDouble                               :: TS, TF
      Type(Vec)                                    :: F
      PetscReal                                    :: Val, tol

      

      !!! Initialize the matrix and vector for the linear system
      Call DMMeshSetMaxDof(AppCtx%MeshTopology%Mesh, 1, iErr); CHKERRQ(iErr) 
      Call DMMeshCreateMatrix(AppCtx%MeshTopology%mesh, AppCtx%U%Sec, MATMPIAIJ, AppCtx%K, iErr); CHKERRQ(iErr)
     

      Call KSPCreate(PETSC_COMM_WORLD, AppCtx%KSP, iErr); CHKERRQ(iErr)
      Call KSPSetOperators(AppCtx%KSP, AppCtx%K, AppCtx%K, SAME_NONZERO_PATTERN, iErr); CHKERRQ(iErr)
      Call KSPSetInitialGuessNonzero(AppCtx%KSP, PETSC_TRUE, iErr); CHKERRQ(iErr)
      Call KSPGetPC(AppCtx%KSP, AppCtx%PC, iErr); CHKERRQ(iErr)
      Call PCSetType(AppCtx%PC, PCBJACOBI, iErr); CHKERRQ(iErr)
      Call KSPSetType(AppCtx%KSP, KSPCG, iErr); CHKERRQ(iErr)
      Call KSPAppendOptionsPrefix(AppCtx%KSP, "U_", iErr); CHKERRQ(iErr)
      Call KSPSetFromOptions(AppCtx%KSP, iErr); CHKERRQ(iErr)
      
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
      
      
      !!! Read Force and BC from Data file or reformat it
      !call PoissonBC(AppCtx) 
      !!! Create the EXO case file
      Call Write_EXO_Case(AppCtx%AppParam%prefix, '%0.4d', MEF90_NumProcs)
      
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine KSPSetUP
   
#undef __FUNCT__
#define __FUNCT__ "TSSetUp"
   Subroutine TSSetUP(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      PetscInt                                     :: iBlk, iDoF      
      PetscBool                                    :: HasPrefix, Flag
      PetscInt                                     :: TmpPoint
      
      PetscInt                                     :: iE, iELoc
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer, filename   
      PetscLogDouble                               :: TS, TF
      Type(DM)                                     :: Tmp_Mesh
      Type(Vec)                                    :: F
      PetscReal                                    :: tol
      PetscInt, Dimension(:), Pointer              :: SizeVect, SizeScal

      

      !!! Allocate the Section for U and F
!      Call DMMeshGetCellSectionReal(AppCtx%MeshTopology%mesh,   'GradU',  AppCtx%MeshTopology%Num_Dim, AppCtx%GradU, iErr); CHKERRQ(iErr)
      

      !!! Initialize the matrix and vector for the linear system
      Call DMMeshSetMaxDof(AppCtx%MeshTopology%Mesh, 1, iErr); CHKERRQ(iErr) 
      Call DMMeshCreateMatrix(AppCtx%MeshTopology%mesh, AppCtx%U%sec, MATMPIAIJ, AppCtx%K, iErr); CHKERRQ(iErr)
      Call DMMeshCreateMatrix(AppCtx%MeshTopology%mesh, AppCtx%U%sec, MATMPIAIJ, AppCtx%M, iErr); CHKERRQ(iErr)
      Call DMMeshCreateMatrix(AppCtx%MeshTopology%mesh, AppCtx%U%sec, MATMPIAIJ, AppCtx%Jac, iErr); CHKERRQ(iErr)
      
      Call TSCreate(PETSC_COMM_WORLD, AppCtx%TS, iErr); CHKERRQ(iErr)
      
      !Call TSSetProblemType(AppCtx%TS, TS_LINEAR, ierr); CHKERRQ(iErr)
      !Call TSSetType(AppCtx%TS, TSBEULER, ierr);  CHKERRQ(iErr)
      call TSSetType(AppCtx%TS, TSARKIMEX, ierr); CHKERRQ(iErr)  


      Call DMMeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'U_0', 1,   AppCtx%U_0, iErr); CHKERRQ(iErr)
      call SectionRealSet(AppCtx%U_0, 0, iErr); CHKERRQ(iErr) !This should be replaced by evaluating the initial solution or most probabaly reading from .args file
      Call DMMeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%U_0, AppCtx%U_0_Vec, iErr); CHKERRQ(iErr)
!TODO check that U_0_vec est bien défini 
      Call TSSetSolution(AppCtx%TS, AppCtx%U_0_Vec,  ierr);   CHKERRQ(iErr) 

      Call TSSetInitialTimeStep(AppCtx%TS, 0.0, 1.,  ierr); CHKERRQ(iErr)
!TODO set these parameters in the PrepPoisson from the.args file 
      AppCtx%maxsteps = 1000
      AppCtx%maxtime =10.0
      Call TSSetDuration(AppCtx%TS, AppCtx%maxsteps, AppCtx%maxtime, iErr); CHKERRQ(iErr)
      Call TSSetFromOptions(AppCtx%TS,  iErr); CHKERRQ(iErr)

      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
       
      
      !!! Read Force and BC from Data file or reformat it
!      call PoissonBC(AppCtx) 
      
      !!! Create the EXO case file
      Call Write_EXO_Case(AppCtx%AppParam%prefix, '%0.4d', MEF90_NumProcs)

      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine TSSetUp


#undef __FUNCT__
#define __FUNCT__ "IFunctionPoisson"
   SubRoutine IFunctionPoisson(dummyTS, t, U, Udot, GlobalOut, AppCtx, iErr)
      PetscReal                                    :: t 
      Type(Vec)                                    :: U, Udot, GlobalOut, dummyvec
      Type(AppCtx_Type)                            :: AppCtx
      Type(TS)                                     :: dummyTS
      PetscInt                                     :: iErr
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer   
         
  
      Write(IOBuffer, *) 'Begin IFunctionPoisson \n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
          
       !Call VecView(global_in,  PETSC_VIEWER_STDOUT_SELF, iErr); CHKERRQ(iErr) 
       !F = AppCtx%M * Udot + AppCtx%K * U
      !Call VecView(U,  PETSC_VIEWER_STDOUT_SELF, iErr); CHKERRQ(iErr)
      Call MatMult(AppCtx%K, U, GlobalOut, iErr);CHKERRQ(iErr)
      !Call VecView(GlobalOut,  PETSC_VIEWER_STDOUT_SELF, iErr); CHKERRQ(iErr)
    !  Call MatMultAdd(AppCtx%M, Udot, GlobalOut, GlobalOut, iErr);CHKERRQ(iErr)
!TODO Can the constant and result factors be the same ? 
      !Call VecView(GlobalOut,  PETSC_VIEWER_STDOUT_SELF, iErr); CHKERRQ(iErr)
       
   End Subroutine IFunctionPoisson


#undef __FUNCT__
#define __FUNCT__ "IJacobianPoisson"
   SubRoutine IJacobianPoisson(dummyTS, t, U, Udot, a, Jac, PreJac, AppCtx, iErr)
      Type(TS)                                     :: dummyTS
      PetscReal                                    :: t, a 
      Type(Vec)                                    :: U, Udot
      Type(Mat)                                    :: Jac, PreJac  
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer   

! a*AppCtx%M + AppCtx%K
      Write(IOBuffer, *) 'Begin IJacobianPoisson \n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

      call MatAXPY(Jac, a, AppCtx%M, AppCtx%K, iErr); CHKERRQ(iErr)

   End Subroutine IJacobianPoisson
 
#undef __FUNCT__
#define __FUNCT__ "RHSPoisson"
   SubRoutine RHSPoisson(dummyTS, t, U, GlobalOut, AppCtx, iErr)
      Type(TS)                                     :: dummyTS
      PetscReal                                    :: t 
      Type(Vec)                                    :: U, GlobalOut
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer   
      Type(VecScatter)                             :: dummyScatter

      Write(IOBuffer, *) 'Begin RHSPoisson \n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
!AppCtx%RHS Is a section real we want to return a vector!
!Return AppCtx%RHS
      Call DMMeshCreateGlobalScatter(AppCtx%MeshTopology%mesh, AppCtx%RHS,  dummyScatter, iErr); CHKERRQ(iErr) 
      Call SectionRealToVec(AppCtx%RHS, dummyScatter, SCATTER_FORWARD, GlobalOut, iErr); CHKERRQ(iErr)
      Call VecScatterDestroy(dummyScatter, iErr); CHKERRQ(iErr)

   End Subroutine RHSPoisson
   
   
#undef __FUNCT__
#define __FUNCT__ "SolveTransient"
   Subroutine SolveTransient(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
      PetscInt                                     :: iDum !, steps, maxsteps
!      PetscReal                                     ::ftime 
      
      Call PetscLogStagePush(AppCtx%LogInfo%KSPSolve_Stage, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(AppCtx%U%Sec, AppCtx%U%Scatter, SCATTER_FORWARD, AppCtx%U%Vec, iErr); CHKERRQ(iErr)

      ! Using IMEX see section 6.1.2 PetSc documentation 
      Call TSSetIFunction(AppCtx%TS, PETSC_NULL_OBJECT, IFunctionPoisson, AppCtx, ierr); CHKERRQ(iErr)
      Call TSSetIJacobian(AppCtx%TS, AppCtx%Jac, AppCtx%Jac, IJacobianPoisson, AppCtx, ierr); CHKERRQ(iErr)
      Call TSSetRHSFunction(AppCtx%TS, PETSC_NULL_OBJECT, RHSPoisson, AppCtx, ierr); CHKERRQ(iErr)
!      call TSSetFromOptions(AppCtx%TS,ierr) 
      call TSSolve(AppCtx%TS, AppCtx%U%Vec, AppCtx%maxtime, iErr); CHKERRQ(iErr)
       

  !    Write(IOBuffer, 100) KSPNumIter, KSPreason
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
!      Call VecDestroy(RHS_Vec, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(AppCtx%U%Sec, AppCtx%U%Scatter, SCATTER_REVERSE, AppCtx%U%Vec, ierr); CHKERRQ(ierr)
 
!      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
!100 Format('KSP converged in ', I5, ' iterations. KSPConvergedReason is ', I2, '\n')
   End Subroutine SolveTransient
   
 
     
#undef __FUNCT__
#define __FUNCT__ "MatMassAssembly"
   Subroutine MatMassAssembly(AppCtx)   
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iBlk, iErr
      
      Call PetscLogStagePush(AppCtx%LogInfo%MatAssembly_Stage, iErr); CHKERRQ(iErr)
      Do_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Call MatMassAssemblyBlock(iBlk, AppCtx)
      End Do Do_iBlk
      Call MatAssemblyBegin(AppCtx%M, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (AppCtx%M, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)

      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine MatMassAssembly
      
#undef __FUNCT__
#define __FUNCT__ "MatMassAssemblyBlock"
   Subroutine MatMassAssemblyBlock(iBlk, AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iBlk
      
      PetscInt                                     :: iE, iELoc, iErr
      PetscReal, Dimension(:,:), Pointer           :: MatElem
      PetscInt, Dimension(:), Pointer              :: BCFlag
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      PetscLogDouble                               :: flops = 0
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
      
      Call PetscLogEventBegin(AppCtx%LogInfo%MatAssemblyBlock_Event, iErr); CHKERRQ(iErr)
      
      Allocate(MatElem(AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF))
      Allocate(BCFlag(AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF))
   
      Do_iELoc: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         MatElem = 0.0_Kr
         BCFlag = 0
!TODO ERROR here : Mass matrice !!!!!!!
         Call SectionIntRestrictClosure(AppCtx%BCFlag%Sec, AppCtx%MeshTopology%mesh, iE-1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF, BCFlag, iErr); CHKERRQ(ierr)
         Write(IOBuffer, *) 'valeur pour le flag BC', BcFlag, '\n'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

         Do iGauss = 1, size(AppCtx%Elem(iE)%Gauss_C)
            Do iDoF1 = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
               If (BCFlag(iDoF1) == 0) Then
                  Do iDoF2 = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
                     MatElem(iDoF2, iDoF1) = MatElem(iDoF2, iDoF1) + AppCtx%Elem(iE)%Gauss_C(iGauss) * ( AppCtx%Elem(iE)%BF(iDoF1, iGauss) *  AppCtx%Elem(iE)%BF(iDoF2, iGauss) )
                     flops = flops + 1
                  End Do
               End If
            End Do
         End Do
         Call DMMeshAssembleMatrix(AppCtx%M, AppCtx%MeshTopology%mesh, AppCtx%U%sec, iE-1, MatElem, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iELoc
   
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
      DeAllocate(BCFlag)
      DeAllocate(MatElem)
      Call PetscLogEventEnd(AppCtx%LogInfo%MatAssemblyBlock_Event, iErr); CHKERRQ(iErr)
   End Subroutine MatMassAssemblyBlock

 
#undef __FUNCT__
#define __FUNCT__ "TSPoissonFinalize"
   Subroutine TSPoissonFinalize(AppCtx)   
      Type(AppCtx_Type)                            :: AppCtx

      PetscInt                                     :: iErr
      Character(len=MEF90_MXSTRLEN)                :: filename

      Call FieldDestroy(AppCtx%U)
      Call FieldDestroy(AppCtx%F)
      Call FieldDestroy(AppCtx%RHS)
      Call SectionRealDestroy(AppCtx%U_0, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(AppCtx%GradU, iErr); CHKERRQ(iErr)
      Call SectionIntDestroy(AppCtx%BCFlag, iErr); CHKERRQ(iErr)
      Call MatDestroy(AppCtx%K, iErr); CHKERRQ(iErr)
      Call MatDestroy(AppCtx%M, iErr); CHKERRQ(iErr)
      Call MatDestroy(AppCtx%Jac, iErr); CHKERRQ(iErr)
      Call VecDestroy(AppCtx%U_0_Vec, iErr); CHKERRQ(iErr)

      Call TSDestroy(AppCtx%TS, iErr); CHKERRQ(iErr)
      Call DMDestroy(AppCtx%MeshTopology%Mesh, iErr); CHKERRQ(ierr)
      
      If (AppCtx%AppParam%verbose > 1) Then
         Call PetscViewerFlush(AppCtx%AppParam%MyLogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerFlush(AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerDestroy(AppCtx%AppParam%MyLogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerDestroy(AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
      End If
      Write(filename, 103) Trim(AppCtx%AppParam%prefix)
!      Call PetscLogPrintSummary(PETSC_COMM_WORLD, filename, iErr); CHKERRQ(iErr)
      Call MEF90_Finalize()
103 Format(A,'.log.txt')
   End Subroutine TSPoissonFinalize

#if defined PB_2D
End Module m_TransientHeat2D
#elif defined PB_3D
End Module m_TransientHeat3D
#endif
