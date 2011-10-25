Program  SimplePoisson

#include "finclude/petscdef.h"

   Use m_MEF90
#if defined PB_2D
   Use m_Poisson2D
#elif defined PB_3D 
   Use m_Poisson3D
#endif

   Implicit NONE   


   Type(AppCtx_Type)                            :: AppCtx
   PetscInt                                     :: iErr
   Character(len=MEF90_MXSTRLEN)                :: IOBuffer
   !PetscReal                                    :: one 
   Call SimplePoissonInit(AppCtx)
   
   If (AppCtx%AppParam%verbose > 4) Then
      Call EXOView(AppCtx%EXO, AppCtx%AppParam%LogViewer)
      Call EXOView(AppCtx%MyEXO, AppCtx%AppParam%MyLogViewer)
      Call MeshTopologyView(AppCtx%MeshTopology, AppCtx%AppParam%MyLogViewer)
   End If   

   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Assembling the matrix\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   Call MatAssembly(AppCtx)
   If (AppCtx%AppParam%verbose > 3) Then
      Write(IOBuffer, *) 'Matrix\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call MatView(AppCtx%K, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   End If

   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Assembling the RHS\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   Call RHSAssembly(AppCtx)
   If (AppCtx%AppParam%verbose > 2) Then
      Write(IOBuffer, *) 'RHS\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call SectionRealView(AppCtx%RHS, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   End If
   
   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Calling Solve\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   Call Solve(AppCtx)
   
   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Computing Energies\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   Call ComputeEnergy(AppCtx)

   Write(IOBuffer, 100) AppCtx%ElasticEnergy, AppCtx%ExtForcesWork, AppCtx%TotalEnergy
100 Format('Elastic energy: ', ES12.5, ' Forces Work: ', ES12.5, ' Total: ', ES12.5, '\n')    
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

   Call ComputeGradU(AppCtx)

   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Saving results\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If

   Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
   Call Write_EXO_Result_Global(AppCtx%MyExo, 1, 1, AppCtx%ElasticEnergy)
   Call Write_EXO_Result_Global(AppCtx%MyExo, 2, 1, AppCtx%ExtForcesWork)
   Call Write_EXO_Result_Global(AppCtx%MyExo, 3, 1, AppCtx%TotalEnergy)
   Call Write_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, 1, 1, AppCtx%U%Sec) 
   Call Write_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, 2, 1, AppCtx%F%Sec) 
   Call Write_EXO_Result_Cell(AppCtx%MyEXO, AppCtx%MeshTopology, 1, 1, AppCtx%GradU) 
   Call PetscLogStagePop (AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
   
   Call SimplePoissonFinalize(AppCtx)

Contains 
#undef __FUNCT__
#define __FUNCT__ "SimplePoissonInit"
   Subroutine SimplePoissonInit(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      PetscInt                                     :: iBlk, iDoF      
      PetscBool                                    :: HasPrefix, Flag
      PetscInt, Dimension(:), Pointer              :: TmpFlag
      PetscInt                                     :: TmpPoint
      
!      Type(SectionReal)                            :: CoordSection
      PetscReal, Dimension(:), Pointer             ::  ValPtr
      PetscReal, Dimension(:,:), Pointer           :: Coords
!      PetscInt                                     :: iE, iELoc
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer, filename   
!      PetscLogDouble                               :: TS, TF
      Type(DM)                                     :: Tmp_Mesh
!      Type(Vec)                                    :: F
      PetscReal                                    :: Val
      PetscInt, Dimension(:), Pointer              :: SizeScal

      Call MEF90_Initialize()
      AppCtx%AppParam%verbose = 0
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-verbose', AppCtx%AppParam%verbose, Flag, iErr); CHKERRQ(iErr)
      
      AppCtx%AppParam%Restart = PETSC_FALSE
      Call PetscOptionsGetBool(PETSC_NULL_CHARACTER, '-restart', AppCtx%AppParam%restart, Flag, iErr); CHKERRQ(iErr)
      
      Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p',       AppCtx%AppParam%prefix, HasPrefix, iErr); CHKERRQ(iErr)
      If (.NOT. HasPrefix) Then
         Call PetscPrintf(PETSC_COMM_WORLD, "No mesh prefix given\n", iErr)
         Call MEF90_Finalize()
         STOP
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
!could we use    Call SectionIntAddNsProperty()  yes 
      
      Allocate(TmpFlag(1))
      Do iBlk = 1, AppCtx%MeshTopology%num_node_sets  
         Do iDoF = 1, AppCtx%MeshTopology%node_set(iBlk)%Num_Nodes     
            TmpPoint = AppCtx%MeshTopology%Num_Elems + AppCtx%MeshTopology%Node_Set(iBlk)%Node_ID(iDoF) - 1
            Call SectionIntUpdate(AppCtx%BCFlag%Sec, TmpPoint, TmpFlag, INSERT_VALUES, iErr); CHKERRQ(iErr)
         End Do
      End Do
      DeAllocate(TmpFlag)

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
      AppCtx%MyEXO%comm = PETSC_COMM_SELF
      AppCtx%MyEXO%exoid = AppCtx%EXO%exoid
      Write(AppCtx%MyEXO%filename, 200) trim(AppCtx%AppParam%prefix), MEF90_MyRank
   200 Format(A, '-', I4.4, '.gen')
      AppCtx%MyEXO%title = trim(AppCtx%EXO%title)
      AppCtx%MyEXO%Num_QA = AppCtx%EXO%Num_QA
      If (AppCtx%AppParam%Restart) Then
         Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
         Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, 1, 1, AppCtx%U) 
         Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, 2, 1, AppCtx%F) 
         Call PetscLogStagePop(iErr); CHKERRQ(iErr)
      Else
         !!! Prepare and format the output mesh   
         Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
         Call MeshTopologyWriteGlobal(AppCtx%MeshTopology, AppCtx%MyEXO, PETSC_COMM_WORLD)
         Call EXOFormat_SimplePoisson(AppCtx)
         Call PetscLogStagePop(iErr); CHKERRQ(iErr)
         
         Select Case (AppCtx%AppParam%TestCase)
         Case(1)
            If (AppCtx%AppParam%verbose > 0) Then
               Write(IOBuffer, *) 'Setting U to 0 and F to 1\n'
               Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
            End If
            !!! U = 0, F=1
            Val = 1.0_Kr
            Call SectionRealSet(AppCtx%F, Val, iErr); CHKERRQ(iErr);
            Val = 0.0_Kr
            Call SectionRealSet(AppCtx%U, Val, iErr); CHKERRQ(iErr);
         Case(2)
            If (AppCtx%AppParam%verbose > 0) Then
               Write(IOBuffer, *) 'Solving Test Case 2: pure dirichlet problem, no force\n'
               Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
            End If

            !!! Test of non homogeneous Dirichlet BC
            Allocate(ValPtr(1))
            ValPtr = 0.0_Kr
            Call SectionRealSet(AppCtx%F, 0.0_Kr, iErr); CHKERRQ(iErr);
            Call DMMeshGetCoordinatesF90(AppCtx%MeshTopology%Mesh, Coords, iErr); CHKERRQ(iErr)
               
            Do iDoF = 1, Size(Coords,1)
#if defined PB_2D
               ValPtr = (Coords(iDoF,1)-0.5_Kr)**2 + (Coords(iDoF,2)+0.5_Kr)**2
               ValPtr = Coords(iDoF,1)**2 + Coords(iDoF,2)**2
               ValPtr = Coords(iDoF,1)**3
#elif defined PB_3D
               ValPtr = (Coords(iDoF,1)-0.5_Kr)**2 + (Coords(iDoF,2)+0.5_Kr)**2 +  (Coords(iDoF,3)+.5_Kr)**2
#endif
               Call SectionRealUpdate(AppCtx%U%Sec, AppCtx%MeshTopology%Num_Elems+iDoF-1, ValPtr, INSERT_VALUES, iErr); CHKERRQ(iErr)
            End Do
            DeAllocate(ValPtr)
            Call DMMeshRestoreCoordinatesF90(AppCtx%MeshTopology%Mesh, Coords, iErr); CHKERRQ(iErr)
         Case(3)
            If (AppCtx%AppParam%verbose > 0) Then
               Write(IOBuffer, *) 'Solving Test Case 3: F=sgn(x) . sgn(y) [. sgn(z)] \n'
               Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
            End If

            !!! Test of non homogeneous Dirichlet BC
            Allocate(ValPtr(1))
            ValPtr = 0.0_Kr
            Call SectionRealSet(AppCtx%U, 1.0_Kr, iErr); CHKERRQ(iErr);
            Call DMMeshGetCoordinatesF90(AppCtx%MeshTopology%Mesh, Coords, iErr); CHKERRQ(iErr)
            
            Do iDoF = 1, Size(Coords,1)
#if defined PB_2D
               If ( Coords(iDoF,1) * Coords(iDoF,2) < 0.0_Kr ) Then
                  ValPtr = -1.0_Kr
               Else
                  ValPtr = 1.0_Kr
               End If
#elif defined PB_3D
               If ( Coords(iDoF,1) * Coords(iDoF,2) * Coords(iDoF,3) < 0.0_Kr ) Then
                  ValPtr = -1.0_Kr
               Else
                  ValPtr = 1.0_Kr
               End If
#endif
               Call SectionRealUpdate(AppCtx%F%Sec, AppCtx%MeshTopology%Num_Elems+iDoF-1, ValPtr, INSERT_VALUES, iErr); CHKERRQ(iErr)
            End Do
            Call DMMeshRestoreCoordinatesF90(AppCtx%MeshTopology%Mesh, Coords, iErr); CHKERRQ(iErr)
         End Select            
      End If
      
      !!! Create the EXO case file
      Call Write_EXO_Case(AppCtx%AppParam%prefix, '%0.4d', MEF90_NumProcs)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine SimplePoissonInit


End Program  SimplePoisson
