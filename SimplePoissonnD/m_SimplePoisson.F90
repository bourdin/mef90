#if defined PB_2D
Module m_SimplePoisson2D
#elif defined PB_3D
Module m_SimplePoisson3D
#endif

#include "finclude/petscdef.h"

   Use m_MEF90
   
   Implicit NONE   

#if defined WITH_TAO
#include "include/finclude/tao_solver.h"
#endif
   
   Type LogInfo_Type
      PetscLogStage                                :: IO_Stage
      PetscLogStage                                :: Setup_Stage 
      PetscLogStage                                :: MeshCreateExodus_Stage
      PetscLogStage                                :: MeshDistribute_Stage
      PetscLogStage                                :: MatAssembly_Stage    
      PetscLogStage                                :: RHSAssembly_Stage
      PetscLogStage                                :: KSPSolve_Stage
      PetscLogStage                                :: PostProc_Stage
      
      PetscLogEvent                                :: MatAssemblyBlock_Event
      PetscLogEvent                                :: RHSAssemblyBlock_Event
      PetscLogEvent                                :: PostProc_Event
   End Type LogInfo_Type

   Type AppParam_Type
      PetscBool                                    :: Restart, Use_tao
      PetscInt                                     :: Verbose
      PetscInt                                     :: TestCase
      Character(len=MEF90_MXSTRLEN)                :: prefix
      Type(PetscViewer)                            :: LogViewer, MyLogViewer
   End Type AppParam_Type

   Type AppCtx_Type
      Type (MeshTopology_Type)                     :: MeshTopology
      Type (EXO_Type)                              :: EXO, MyEXO
#if defined PB_2D
      Type(Element2D_Scal), Dimension(:), Pointer  :: Elem
#elif defined PB_3D
      Type(Element3D_Scal), Dimension(:), Pointer  :: Elem
#endif
      Type(SectionReal)                            :: U
      Type(Vec)                                    :: U_Vec
      Type(SectionReal)                            :: GradU
      Type(SectionReal)                            :: F
      PetscReal                                    :: ElasticEnergy
      PetscReal                                    :: ExtForcesWork
      PetscReal                                    :: TotalEnergy
      Type(VecScatter)                             :: Scatter
      Type(SectionInt)                             :: BCFlag
      Type(Mat)                                    :: K
      Type(SectionReal)                            :: RHS
      Type(KSP)                                    :: KSP
      Type(PC)                                     :: PC
      Type(LogInfo_Type)                           :: LogInfo
      Type(AppParam_Type)                          :: AppParam
#if defined WITH_TAO
      TAO_SOLVER                                   :: taoU
      TAO_APPLICATION                              :: taoappU
      Type(Vec)                                    :: Uplus, Uminus
#endif
   End Type AppCtx_Type
   
   
Contains

   Subroutine SimplePoissonInit(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      PetscInt                                     :: iBlk, iDoF      
      PetscBool                                    :: HasPrefix, Flag
      PetscInt, Dimension(:), Pointer              :: TmpFlag
      PetscInt                                     :: TmpPoint
      
      Type(SectionReal)                            :: CoordSection
      PetscReal, Dimension(:), Pointer             :: TmpCoords, ValPtr
      PetscReal, Dimension(:,:), Pointer           :: Coords
      PetscInt                                     :: iE, iELoc
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer, filename   
      PetscLogDouble                               :: TS, TF
      Type(Mesh)                                   :: Tmp_Mesh
      Type(Vec)                                    :: F
      PetscReal                                    :: Val, tol

      Call MEF90_Initialize()
#if defined WITH_TAO
      Call TaoInitialize(PETSC_NULL_CHARACTER, iErr); CHKERRQ(iErr)
#endif
      AppCtx%AppParam%verbose = 0
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-verbose', AppCtx%AppParam%verbose, Flag, iErr); CHKERRQ(iErr)
      
      AppCtx%AppParam%Restart = PETSC_FALSE
      Call PetscOptionsGetBool(PETSC_NULL_CHARACTER, '-restart', AppCtx%AppParam%restart, Flag, iErr); CHKERRQ(iErr)
      
      AppCtx%AppParam%Use_Tao = PETSC_FALSE
      Call PetscOptionsGetBool(PETSC_NULL_CHARACTER, '-tao', AppCtx%AppParam%Use_Tao, Flag, iErr); CHKERRQ(iErr)


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
      AppCtx%EXO%filename = Trim(AppCtx%AppParam%prefix)//'.gen'
      !!! Read and partition the mesh
      If (MEF90_NumProcs == 1) Then
         Call PetscLogStagePush(AppCtx%LogInfo%MeshCreateExodus_Stage, iErr); CHKERRQ(iErr)
         Call MeshCreateExodus(PETSC_COMM_WORLD, AppCtx%EXO%filename, AppCtx%MeshTopology%mesh, ierr); CHKERRQ(iErr)
         Call PetscLogStagePop(iErr); CHKERRQ(iErr)
      Else
         Call PetscLogStagePush(AppCtx%LogInfo%MeshCreateExodus_Stage, iErr); CHKERRQ(iErr)
         Call MeshCreateExodus(PETSC_COMM_WORLD, AppCtx%EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
         Call PetscLogStagePop(iErr); CHKERRQ(iErr)
      
         Call PetscLogStagePush(AppCtx%LogInfo%MeshDistribute_Stage, iErr); CHKERRQ(iErr)
         Call MeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, AppCtx%MeshTopology%mesh, ierr); CHKERRQ(iErr)
         Call MeshDestroy(Tmp_mesh, ierr); CHKERRQ(iErr)
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
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'U', 1,   AppCtx%U, iErr); CHKERRQ(iErr)
      Call MeshGetCellSectionReal(AppCtx%MeshTopology%mesh,   'GradU',  AppCtx%MeshTopology%Num_Dim, AppCtx%GradU, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'F', 1,   AppCtx%F, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'RHS', 1, AppCtx%RHS, iErr); CHKERRQ(iErr)
      Call MeshCreateGlobalScatter(AppCtx%MeshTopology%mesh, AppCtx%U, AppCtx%Scatter, iErr); CHKERRQ(iErr)
      Call MeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%U, AppCtx%U_Vec, iErr); CHKERRQ(iErr)
#if defined WITH_TAO
      If (AppCtx%AppParam%Use_Tao) Then
         Call MeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%U, AppCtx%Uplus, iErr); CHKERRQ(iErr)
         Call MeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%U, AppCtx%Uminus, iErr); CHKERRQ(iErr)
      End If
#endif

      !!! Allocate and initialize the Section for the flag
      Call MeshGetVertexSectionInt(AppCtx%MeshTopology%mesh, 'BCFlag', 1, AppCtx%BCFlag, iErr); CHKERRQ(iErr)
      Allocate(TmpFlag(1))
      Do iBlk = 1, AppCtx%MeshTopology%num_node_sets  
         Do iDoF = 1, AppCtx%MeshTopology%node_set(iBlk)%Num_Nodes     
            TmpPoint = AppCtx%MeshTopology%Num_Elems + AppCtx%MeshTopology%Node_Set(iBlk)%Node_ID(iDoF) - 1
            Call SectionIntUpdate(AppCtx%BCFlag, TmpPoint, TmpFlag, INSERT_VALUES, iErr); CHKERRQ(iErr)
         End Do
      End Do
      DeAllocate(TmpFlag)

      !!! Initialize the matrix and vector for the linear system
      Call MeshSetMaxDof(AppCtx%MeshTopology%Mesh, 1, iErr); CHKERRQ(iErr) 
      Call MeshCreateMatrix(AppCtx%MeshTopology%mesh, AppCtx%U, MATMPIAIJ, AppCtx%K, iErr); CHKERRQ(iErr)
      
      If (AppCtx%AppParam%Use_Tao) Then
#if defined WITH_TAO
         Call TaoCreate(PETSC_COMM_WORLD, 'tao_tron', AppCTx%taoU, iErr); CHKERRQ(iErr)      
         Call TaoApplicationCreate(PETSC_COMM_WORLD, AppCTx%taoappU, iErr); CHKERRQ(iErr)      
         Call TaoAppGetKSP(AppCtx%taoappU, AppCtx%KSP, iErr); CHKERRQ(iErr)
         Call KSPAppendOptionsPrefix(AppCtx%KSP, "U_", iErr); CHKERRQ(iErr)
         Call KSPSetType(AppCtx%KSP, KSPCG, iErr); CHKERRQ(iErr)
         Call KSPSetFromOptions(AppCtx%KSP, iErr); CHKERRQ(iErr)

         Call KSPGetPC(AppCtx%KSP, AppCtx%PC, iErr); CHKERRQ(iErr)
         Call PCSetType(AppCtx%PC, PCBJACOBI, iErr); CHKERRQ(iErr)
         Call PCSetFromOptions(AppCtx%PC, iErr); CHKERRQ(iErr)

         Call TaoAppSetObjectiveAndGradientRoutine(AppCtx%taoappU, FormFunctionAndGradient, AppCtx, iErr); CHKERRQ(iErr)
         Call TaoAppSetHessianRoutine(AppCtx%taoappU, FormHessian, AppCtx, iErr); CHKERRQ(iErr)

         Call TaoAppSetHessianMat(AppCtx%taoappU, AppCtx%K, AppCtx%K, iErr); CHKERRQ(iErr)
         Call TaoAppSetInitialSolutionVec(AppCtx%taoappU, AppCtx%U_Vec, iErr); CHKERRQ(iErr)
         
         tol = 1.0D-8
         Call TaoSetTolerances(AppCtx%taoU, tol, tol, tol, tol, iErr); CHKERRQ(iErr)
         Call TaoSetGradientTolerances(AppCtx%taoU, tol, tol, tol, iErr); CHKERRQ(iErr)

         Call TaoSetOptions(AppCtx%taoappU, AppCtx%taoU, iErr); CHKERRQ(iErr)
#endif
      Else      
      !!! Create the KSP and PCs
         Call KSPCreate(PETSC_COMM_WORLD, AppCtx%KSP, iErr); CHKERRQ(iErr)
         Call KSPSetOperators(AppCtx%KSP, AppCtx%K, AppCtx%K, SAME_NONZERO_PATTERN, iErr); CHKERRQ(iErr)
         Call KSPSetInitialGuessNonzero(AppCtx%KSP, PETSC_TRUE, iErr); CHKERRQ(iErr)
         Call KSPGetPC(AppCtx%KSP, AppCtx%PC, iErr); CHKERRQ(iErr)
         Call PCSetType(AppCtx%PC, PCBJACOBI, iErr); CHKERRQ(iErr)
         Call KSPSetType(AppCtx%KSP, KSPCG, iErr); CHKERRQ(iErr)
         Call KSPAppendOptionsPrefix(AppCtx%KSP, "U_", iErr); CHKERRQ(iErr)
         
         
         Call KSPSetFromOptions(AppCtx%KSP, iErr); CHKERRQ(iErr)
      End If
      
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
         Call Write_MeshTopologyGlobal(AppCtx%MeshTopology, AppCtx%MyEXO, PETSC_COMM_WORLD)
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
            Call MeshGetCoordinatesF90(AppCtx%MeshTopology%Mesh, Coords, iErr); CHKERRQ(iErr)
               
            Do iDoF = 1, Size(Coords,1)
#if defined PB_2D
               ValPtr = (Coords(iDoF,1)-0.5_Kr)**2 + (Coords(iDoF,2)+0.5_Kr)**2
               ValPtr = Coords(iDoF,1)**2 + Coords(iDoF,2)**2
               ValPtr = Coords(iDoF,1)**3
#elif defined PB_3D
               ValPtr = (Coords(iDoF,1)-0.5_Kr)**2 + (Coords(iDoF,2)+0.5_Kr)**2 +  (Coords(iDoF,3)+.5_Kr)**2
#endif
               Call SectionRealUpdate(AppCtx%U, AppCtx%MeshTopology%Num_Elems+iDoF-1, ValPtr, INSERT_VALUES, iErr); CHKERRQ(iErr)
            End Do
            DeAllocate(ValPtr)
            Call MeshRestoreCoordinatesF90(AppCtx%MeshTopology%Mesh, Coords, iErr); CHKERRQ(iErr)
         Case(3)
            If (AppCtx%AppParam%verbose > 0) Then
               Write(IOBuffer, *) 'Solving Test Case 3: F=sgn(x) . sgn(y) [. sgn(z)] \n'
               Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
            End If

            !!! Test of non homogeneous Dirichlet BC
            Allocate(ValPtr(1))
            ValPtr = 0.0_Kr
            Call SectionRealSet(AppCtx%U, 1.0_Kr, iErr); CHKERRQ(iErr);
            Call MeshGetCoordinatesF90(AppCtx%MeshTopology%Mesh, Coords, iErr); CHKERRQ(iErr)
            
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
               Call SectionRealUpdate(AppCtx%F, AppCtx%MeshTopology%Num_Elems+iDoF-1, ValPtr, INSERT_VALUES, iErr); CHKERRQ(iErr)
            End Do
            Call MeshRestoreCoordinatesF90(AppCtx%MeshTopology%Mesh, Coords, iErr); CHKERRQ(iErr)
         End Select            
      End If
      
      !!! Create the EXO case file
      Call Write_EXO_Case(AppCtx%AppParam%prefix, '%0.4d', MEF90_NumProcs)

#if defined WITH_TAO
      If (AppCtx%AppParam%Use_Tao) Then   
         Call SetupBCBounds(AppCtx, iErr)
      End If
#endif      

      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine SimplePoissonInit
   
   Subroutine EXOFormat_SimplePoisson(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr
   
      AppCtx%MyEXO%exoid = EXOPEN(AppCtx%MyEXO%filename, EXWRIT, exo_cpu_ws, exo_io_ws, exo_ver, iErr)
      Call EXPVP (AppCtx%MyEXO%exoid, 'g', 3, iErr)
      Call EXPVAN(AppCtx%MyEXO%exoid, 'g', 3, (/'Elastic Energy ', 'Ext Forces work', 'Total Energy   '/), iErr)
      Call EXPVP (AppCtx%MyEXO%exoid, 'n', 2, iErr)
      Call EXPVAN(AppCtx%MyEXO%exoid, 'n', 2, (/'U', 'F'/), iErr)
#if defined PB_2D
      Call EXPVP (AppCtx%MyEXO%exoid, 'e', 2, iErr)
      Call EXPVAN(AppCtx%MyEXO%exoid, 'e', 2, (/'Grad U_X', 'Grad U_Y'/), iErr)
#elif defined PB_3D
      Call EXPVP (AppCtx%MyEXO%exoid, 'e', 3, iErr)
      Call EXPVAN(AppCtx%MyEXO%exoid, 'e', 3, (/'Grad U_X', 'Grad U_Y', 'Grad U_Z'/), iErr)
#endif
      Call EXPTIM(AppCtx%MyEXO%exoid, 1, 1.0_Kr, iErr)

      Call EXCLOS(AppCtx%MyEXO%exoid, iErr)
   End Subroutine EXOFormat_SimplePoisson
   
   Subroutine InitLog(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr
      
      Call PetscLogEventRegister('MatAssembly Block', 0, AppCtx%LogInfo%MatAssemblyBlock_Event, ierr); CHKERRQ(ierr)
      Call PetscLogEventRegister('RHSAssembly Block', 0, AppCtx%LogInfo%RHSAssemblyBlock_Event, ierr); CHKERRQ(ierr)
      Call PetscLogEventRegister('Post Processing',   0, AppCtx%LogInfo%PostProc_Event,         ierr); CHKERRQ(ierr)

      Call PetscLogStageRegister("MeshCreateExodus", AppCtx%LogInfo%MeshCreateExodus_Stage, iErr)
      Call PetscLogStageRegister("MeshDistribute",   AppCtx%LogInfo%MeshDistribute_Stage,   iErr)
      Call PetscLogStageRegister("IO Stage",         AppCtx%LogInfo%IO_Stage,               iErr)
      Call PetscLogStageRegister("Setup",            AppCtx%LogInfo%Setup_Stage,            iErr)
      Call PetscLogStageRegister("Mat Assembly",     AppCtx%LogInfo%MatAssembly_Stage,      iErr)
      Call PetscLogStageRegister("RHS Assembly",     AppCtx%LogInfo%RHSAssembly_Stage,      iErr)
      Call PetscLogStageRegister("KSP Solve",        AppCtx%LogInfo%KSPSolve_Stage,         iErr)
      Call PetscLogStageRegister("Post Proc",        AppCtx%LogInfo%PostProc_Stage,         iErr)
   End Subroutine InitLog
   
   Subroutine Solve(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      KSPConvergedReason                           :: KSPreason
#if defined WITH_TAO
      TaoTerminateReason                           :: TaoReason
#endif
      PetscReal                                    :: TaoResidual
      PetscInt                                     :: KSPNumIter
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
      PetscInt                                     :: iDum
      Type(Vec)                                    :: RHS_Vec
      
      Call PetscLogStagePush(AppCtx%LogInfo%KSPSolve_Stage, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(AppCtx%U, AppCtx%Scatter, SCATTER_FORWARD, AppCtx%U_Vec, iErr); CHKERRQ(iErr)

      If (AppCtx%AppParam%Use_Tao) Then
#if defined WITH_TAO
         Call TaoSolveApplication(AppCtx%taoappU, AppCtx%taoU, iErr); CHKERRQ(iErr)
         Call TaoGetTerminationReason(AppCtx%TaoU, TaoReason, iErr); CHKERRQ(iErr)
         If ( TaoReason > 0) Then
            Write(IOBuffer, *) 'Tao convergence status: ', TaoReason, '\n'
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Else
            Write(IOBuffer, 103) TaoReason
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         End If
103 Format('[ERROR] TaoSolveApplication did not converge. ', I2, '\n')      
#else
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "TAO not supported in this build.", iErr); CHKERRQ(iErr)
#endif
      Else
         Call MeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%RHS, RHS_Vec, iErr); CHKERRQ(iErr)
         Call SectionRealToVec(AppCtx%RHS, AppCtx%Scatter, SCATTER_FORWARD, RHS_Vec, iErr); CHKERRQ(iErr)

         Call KSPSolve(AppCtx%KSP, RHS_Vec, AppCtx%U_Vec, iErr); CHKERRQ(iErr)
         !!! Solve and store the solution in AppCtx%RHS
         
         Call KSPGetIterationNumber(AppCtx%KSP, KSPNumIter, iErr); CHKERRQ(iErr)
         Call KSPGetConvergedReason(AppCtx%KSP, KSPreason, iErr); CHKERRQ(iErr)
         Write(IOBuffer, 100) KSPNumIter, KSPreason
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
!         Call SectionRealToVec(AppCtx%RHS, AppCtx%Scatter, SCATTER_REVERSE, RHS_Vec, iErr); CHKERRQ(iErr)
         Call VecDestroy(RHS_Vec, iErr); CHKERRQ(iErr)
      End If
      Call SectionRealToVec(AppCtx%U, AppCtx%Scatter, SCATTER_REVERSE, AppCtx%U_Vec, ierr); CHKERRQ(ierr)
 
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
100 Format('KSP converged in ', I5, ' iterations. KSPConvergedReason is ', I2, '\n')
   End Subroutine Solve
   
      
   Subroutine MatAssembly(AppCtx)   
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iBlk, iErr
      
      Call PetscLogStagePush(AppCtx%LogInfo%MatAssembly_Stage, iErr); CHKERRQ(iErr)
      Do_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Call MatAssemblyBlock(iBlk, AppCtx)
      End Do Do_iBlk
      Call MatAssemblyBegin(AppCtx%K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (AppCtx%K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)

      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine MatAssembly
      
   Subroutine MatAssemblyBlock(iBlk, AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iBlk
      
      PetscInt                                     :: iE, iELoc, iErr
      PetscReal, Dimension(:,:), Pointer           :: MatElem
      PetscInt, Dimension(:), Pointer              :: BCFlag
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      PetscLogDouble                               :: flops = 0
      
      Call PetscLogEventBegin(AppCtx%LogInfo%MatAssemblyBlock_Event, iErr); CHKERRQ(iErr)
      
      Allocate(MatElem(AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF))
      Allocate(BCFlag(AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF))
   
      Do_iELoc: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         MatElem = 0.0_Kr
         BCFlag = 0
         Call SectionIntRestrictClosure(AppCtx%BCFlag, AppCtx%MeshTopology%mesh, iE-1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF, BCFlag, iErr); CHKERRQ(ierr)
         Write(MEF90_MyRank+300, *) iE, BCFlag
         Do iGauss = 1, size(AppCtx%Elem(iE)%Gauss_C)
            Do iDoF1 = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
               If (BCFlag(iDoF1) == 0) Then
                  Do iDoF2 = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
                     MatElem(iDoF2, iDoF1) = MatElem(iDoF2, iDoF1) + AppCtx%Elem(iE)%Gauss_C(iGauss) * ( AppCtx%Elem(iE)%Grad_BF(iDoF1, iGauss) .DotP. AppCtx%Elem(iE)%Grad_BF(iDoF2, iGauss) )
                     flops = flops + 1
                  End Do
               End If
            End Do
         End Do
         Call assembleMatrix(AppCtx%K, AppCtx%MeshTopology%mesh, AppCtx%U, iE-1, MatElem, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iELoc
   
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
      DeAllocate(BCFlag)
      DeAllocate(MatElem)
      Call PetscLogEventEnd(AppCtx%LogInfo%MatAssemblyBlock_Event, iErr); CHKERRQ(iErr)
   End Subroutine MatAssemblyBlock

   Subroutine RHSAssembly(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      PetscInt                                     :: iBlk

      Call PetscLogStagePush(AppCtx%LogInfo%RHSAssembly_Stage, iErr); CHKERRQ(iErr)

      Call SectionRealZero(AppCtx%RHS, iErr); CHKERRQ(iErr)
      
      Do_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Call RHSAssemblyBlock(iBlk, AppCtx)
      End Do Do_iBlk

      Call SectionRealComplete(AppCtx%RHS, iErr); CHKERRQ(iErr)
      !!! VERY important! This is the equivalent of a ghost update
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine RHSAssembly



   Subroutine RHSAssemblyBlock(iBlk, AppCtx)
      PetscInt                                     :: iBlk
      Type(AppCtx_Type)                            :: AppCtx

      PetscInt                                     :: iErr
      PetscInt                                     :: iBlkId
      PetscInt                                     :: Num_DoF, iDoF
      PetscInt                                     :: iEloc, iE, iGauss
      PetscInt, Dimension(:), Pointer              :: BCFlag_Loc
      PetscReal, Dimension(:), Pointer             :: RHS_Loc, F_Loc
      PetscReal                                    :: F_Elem
      PetscLogDouble                               :: flops 
      
      Call PetscLogEventBegin(AppCtx%LogInfo%RHSAssemblyBlock_Event, iErr); CHKERRQ(iErr)
      flops = 0.0

      Num_DoF = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
      Allocate(F_Loc(Num_DoF))
      Allocate(RHS_Loc(Num_DoF))
      Allocate(BCFlag_Loc(Num_DoF))

      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
      Do_iEloc: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         RHS_Loc = 0.0_Kr
         Call SectionRealRestrictClosure(AppCtx%F, AppCtx%MeshTopology%mesh, iE-1, Num_DoF, F_Loc, iErr); CHKERRQ(ierr)
         Call SectionIntRestrictClosure(AppCtx%BCFlag, AppCtx%MeshTopology%mesh, iE-1, Num_DoF, BCFlag_Loc, iErr); CHKERRQ(ierr)
         Do iGauss = 1, Size(AppCtx%Elem(iE)%Gauss_C)
            F_Elem = 0.0_Kr
            Do iDoF = 1, Num_DoF
               F_Elem = F_Elem + AppCtx%Elem(iE)%BF(iDoF, iGauss) * F_Loc(iDoF)
               flops = flops + 2.0
            End Do
            Do iDoF = 1, Num_DoF
               If (BCFlag_Loc(iDoF) == 0) Then
                  RHS_Loc(iDoF) = RHS_Loc(iDoF) + AppCtx%Elem(iE)%Gauss_C(iGauss) * ( F_Elem * AppCtx%Elem(iE)%BF(iDoF, iGauss) )
                  flops = flops + 3.0
               End If
            End Do
         End Do
         Call SectionRealUpdateClosure(AppCtx%RHS, AppCtx%MeshTopology%Mesh, iE-1, RHS_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iEloc

      DeAllocate(BCFlag_Loc)
      DeAllocate(RHS_Loc)
      DeAllocate(F_Loc)
      Call PetscLogFlops(flops , iErr);CHKERRQ(iErr)
      Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyBlock_Event, iErr); CHKERRQ(iErr)
   End Subroutine RHSAssemblyBlock

   Subroutine ComputeEnergy(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoF, NumGauss
      PetscReal, Dimension(:), Pointer             :: F, U
      PetscInt                                     :: iBlk, iELoc, iE
      PetscInt                                     :: iDoF, iGauss
#if defined PB_2D
      Type(Vect2D)                                 :: Strain_Elem, Stress_Elem      
#elif defined PB_3D
      Type(Vect3D)                                 :: Strain_Elem, Stress_Elem      
#endif
      PetscReal                                    :: F_Elem, U_Elem
      PetscReal                                    :: MyElasticEnergy, MyExtForcesWork
      PetscLogDouble                               :: flops

      Call PetscLogStagePush(AppCtx%LogInfo%PostProc_Stage, iErr); CHKERRQ(iErr)
      Call PetscLogEventBegin(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      
      MyElasticEnergy = 0.0_Kr
      MyExtForcesWork = 0.0_Kr
      flops = 0.0_Kr
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            NumDoF   = Size(AppCtx%Elem(iE)%BF,1)
            NumGauss = Size(AppCtx%Elem(iE)%BF,2)
            Allocate(F(NumDoF))
            Call SectionRealRestrictClosure(AppCtx%F, AppCtx%MeshTopology%mesh, iE-1, NumDoF, F, iErr); CHKERRQ(ierr)
            Allocate(U(NumDoF))
            Call SectionRealRestrictClosure(AppCtx%U, AppCtx%MeshTopology%mesh, iE-1, NumDoF, U, iErr); CHKERRQ(ierr)
            Do iGauss = 1, NumGauss
               Strain_Elem = 0.0_Kr
               Stress_Elem = 0.0_Kr
               F_Elem      = 0.0_Kr
               U_Elem      = 0.0_Kr
               Do iDoF = 1, NumDoF
                  Stress_Elem = Stress_Elem + AppCtx%Elem(iE)%Grad_BF(iDoF, iGauss) * U(iDoF)
                  Strain_Elem = Strain_Elem + AppCtx%Elem(iE)%Grad_BF(iDoF, iGauss) * U(iDoF)
                  F_Elem = F_Elem + AppCtx%Elem(iE)%BF(iDoF, iGauss) * F(iDoF)
                  U_Elem = U_Elem + AppCtx%Elem(iE)%BF(iDoF, iGauss) * U(iDoF)
                  flops = flops + 4
               End Do
               MyElasticEnergy = MyElasticEnergy + AppCtx%Elem(iE)%Gauss_C(iGauss) * ( (Strain_Elem .DotP. Strain_Elem) * 0.5_Kr)
               MyExtForcesWork = MyExtForcesWork - AppCtx%Elem(iE)%Gauss_C(iGauss) * F_Elem * U_Elem
               flops = flops + 5
            End Do
            DeAllocate(U)
            DeAllocate(F)
         End Do Do_Elem_iE
      End Do Do_Elem_iBlk
      Call MPI_AllReduce(MyElasticEnergy, AppCtx%ElasticEnergy, 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      Call MPI_AllReduce(MyExtForcesWork, AppCtx%ExtForcesWork, 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      AppCtx%TotalEnergy = AppCtx%ElasticEnergy + AppCtx%ExtForcesWork
      Call PetscLogFlops(flops, iErr)
      Call PetscLogEventEnd  (AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine ComputeEnergy
   
   Subroutine ComputeGradU(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
#if defined PB_2D
      Type(Vect2D)                                 :: Grad
#elif defined PB_3D
      Type(Vect3D)                                 :: Grad
#endif
      PetscReal                                    :: Vol
      PetscInt                                     :: NumDoF, NumGauss
      PetscReal, Dimension(:), Pointer             :: U
      PetscInt                                     :: iBlk, iELoc, iE
      PetscInt                                     :: iDoF, iGauss
      PetscReal, Dimension(:), Pointer             :: Grad_Ptr
      PetscLogDouble                               :: flops = 0

      Call PetscLogStagePush (AppCtx%LogInfo%PostProc_Stage, iErr); CHKERRQ(iErr)
      Call PetscLogEventBegin(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      Allocate(Grad_Ptr(AppCtx%MeshTopology%Num_Dim))
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            NumDoF   = Size(AppCtx%Elem(iE)%BF,1)
            NumGauss = Size(AppCtx%Elem(iE)%BF,2)
            Allocate(U(NumDoF))
            Call SectionRealRestrictClosure(AppCtx%U, AppCtx%MeshTopology%mesh, iE-1, NumDoF, U, iErr); CHKERRQ(ierr)
            Grad = 0.0_Kr
            Vol  = 0.0_Kr
            Do iGauss = 1, NumGauss
               Do iDoF = 1, NumDoF
                  Grad = Grad + (AppCtx%Elem(iE)%Gauss_C(iGauss) * U(iDoF) ) * AppCtx%Elem(iE)%Grad_BF(iDoF, iGauss) 
                  Vol = Vol + AppCtx%Elem(iE)%Gauss_C(iGauss) * AppCtx%Elem(iE)%BF(iDoF, iGauss)
                  flops = flops + 3
               End Do
            End Do
            Grad = Grad / Vol
#if defined PB_2D
            Grad_Ptr = (/ Grad%X, Grad%Y /)
#elif defined PB_3D
            Grad_Ptr = (/ Grad%X, Grad%Y, Grad%Z /)
#endif
            Call SectionRealUpdateClosure(AppCtx%GradU, AppCtx%MeshTopology%Mesh, iE-1, Grad_Ptr, INSERT_VALUES, iErr)
            DeAllocate(U)
         End Do Do_Elem_iE
      End Do Do_Elem_iBlk
      DeAllocate(Grad_Ptr)
      Call PetscLogFlops(flops, iErr)
      Call PetscLogEventEnd  (AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine ComputeGradU
  
  
#if defined WITH_TAO 
   Subroutine FormHessian(tao, X, H, Hpre, flg, AppCtx, iErr)
      TAO_SOLVER         :: tao
      Type(Vec)          :: X
      Type(Mat)          :: H, Hpre
      PetscInt           :: iErr
      MatStructure       :: flg
      Type(AppCtx_Type)  :: AppCtx
      
      PetscInt           :: iBlk
      
      Call PetscLogStagePush(AppCtx%LogInfo%MatAssembly_Stage, iErr); CHKERRQ(iErr)

      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Call HessianAssemblyBlock(iBlk, H, AppCtx)
         Call MatAssemblyBegin(H, MAT_FLUSH_ASSEMBLY, iErr); CHKERRQ(iErr)
         Call MatAssemblyEnd  (H, MAT_FLUSH_ASSEMBLY, iErr); CHKERRQ(iErr)
      End Do Do_Elem_iBlk
      Call MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (H, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)

      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine FormHessian
#endif

#if defined WITH_TAO
   Subroutine HessianAssemblyBlock(iBlk, H, AppCtx)
      PetscInt                                     :: iBlk
      Type(Mat)                                    :: H
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iE, iELoc, iErr
      PetscReal, Dimension(:,:), Pointer           :: MatElem
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      PetscLogDouble                               :: flops = 0
      
      Call PetscLogEventBegin(AppCtx%LogInfo%MatAssemblyBlock_Event, iErr); CHKERRQ(iErr)
      
      Allocate(MatElem(AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF))
   
      Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         MatElem = 0.0_Kr
         Do iGauss = 1, size(AppCtx%Elem(iE)%Gauss_C)
            Do iDoF1 = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
                  Do iDoF2 = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
                     MatElem(iDoF2, iDoF1) = MatElem(iDoF2, iDoF1) + AppCtx%Elem(iE)%Gauss_C(iGauss) * ( AppCtx%Elem(iE)%Grad_BF(iDoF1, iGauss) .DotP. AppCtx%Elem(iE)%Grad_BF(iDoF2, iGauss) )
                     flops = flops + 1
                  End Do
            End Do
         End Do
         Call assembleMatrix(H, AppCtx%MeshTopology%mesh, AppCtx%U, iE-1, MatElem, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_Elem_iE
   
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
      DeAllocate(MatElem)
      Call PetscLogEventEnd(AppCtx%LogInfo%MatAssemblyBlock_Event, iErr); CHKERRQ(iErr)
   End Subroutine HessianAssemblyBlock
#endif

#if defined WITH_TAO
   Subroutine FormFunctionandGradient(tao, X, func, Gradient, AppCtx, iErr)
      TAO_SOLVER                                   :: tao
      Type(Vec)                                    :: X, Gradient
      PetscInt                                     :: iErr
      PetscReal                                    :: func
      Type(AppCtx_Type)                            :: AppCtx
      
      Type(SectionReal)                            :: Gradient_Sec, X_Sec
      PetscInt                                     :: iBlk 
      PetscReal                                    :: myfunc

      Call PetscLogStagePush(AppCtx%LogInfo%RHSAssembly_Stage, iErr); CHKERRQ(iErr)

      Call VecZeroEntries(Gradient, iErr); CHKERRQ(iErr)

      !!! Create Section to temporarily hold the surrent point and gradient
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'Gradient', 1, Gradient_Sec, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'X', 1, X_Sec, iErr); CHKERRQ(iErr)
      
      !!! Scatter values from the Vec to the Sections
      Call SectionRealToVec(X_Sec, AppCtx%Scatter, SCATTER_REVERSE, X, iErr); CHKERRQ(iErr)
      
      myfunc = 0.0_Kr
      Do_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Call FormFunctionAndGradientBlock(iBlk, X_Sec, myfunc, Gradient_Sec, AppCtx)
      End Do Do_iBlk
      Call PetscGlobalSum(Myfunc, func, PETSC_COMM_WORLD,iErr); CHKERRQ(iErr)

      Call SectionRealComplete(Gradient_Sec, iErr); CHKERRQ(iErr)
      !!! Scatter values from the Sections back to the Vec
      Call SectionRealToVec(Gradient_Sec, AppCtx%Scatter, SCATTER_FORWARD, Gradient, iErr); CHKERRQ(iErr)

      Call SectionRealDestroy(Gradient_Sec, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(X_Sec, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine FormFunctionandGradient
#endif

#if defined WITH_TAO
   Subroutine FormFunctionAndGradientBlock(iBlk, X_Sec, func, Gradient_Sec, AppCtx)
      PetscInt                                     :: iBlk
      Type(SectionReal)                            :: X_Sec, Gradient_Sec
      PetscReal                                    :: func
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      PetscInt                                     :: iELoc, iE, iBlkID, iGauss
      PetscInt                                     :: Num_DoF, iDoF
      PetscReal, Dimension(:), Pointer             :: X_Loc, F_Loc, Gradient_Loc
      PetscReal                                    :: X_Elem, F_Elem
#if defined PB_2D
      Type (Vect2D)                                :: Strain_Elem
#elif defined PB_3D  
      Type (Vect3D)                                :: Strain_Elem    
#endif
      PetscLogDouble                               :: flops
      
      Call PetscLogEventBegin(AppCtx%LogInfo%RHSAssemblyBlock_Event, iErr); CHKERRQ(iErr)
      flops = 0.0
      
      Num_DoF = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
      Allocate(X_Loc(Num_DoF))
      Allocate(F_Loc(Num_DoF))
      Allocate(Gradient_Loc(Num_DoF))

      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID

      Do_iEloc: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         Gradient_Loc = 0.0_Kr
         Call SectionRealRestrictClosure(X_Sec, AppCtx%MeshTopology%mesh, iE-1, Num_DoF, X_Loc, iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(AppCtx%F, AppCtx%MeshTopology%mesh, iE-1, Num_DoF, F_Loc, iErr); CHKERRQ(ierr)
         Do iGauss = 1, Size(AppCtx%Elem(iE)%Gauss_C)
            Strain_Elem = 0.0_Kr
            F_Elem = 0.0_Kr
            X_Elem = 0.0_Kr
            Do iDoF = 1, Num_DoF
               F_Elem = F_Elem + AppCtx%Elem(iE)%BF(iDoF, iGauss) * F_Loc(iDoF)
               X_Elem = X_Elem + AppCtx%Elem(iE)%BF(iDoF, iGauss) * X_Loc(iDoF)
               Strain_Elem = Strain_Elem + AppCtx%Elem(iE)%Grad_BF(iDoF, iGauss) * X_Loc(iDoF)
               flops = flops + 2.0
            End Do
            Do iDoF = 1, Num_DoF
               Gradient_Loc(iDoF) = Gradient_Loc(iDoF) + AppCtx%Elem(iE)%Gauss_C(iGauss) * ( (Strain_Elem .DotP. AppCtx%Elem(iE)%Grad_BF(iDoF, iGauss)) - F_Elem * AppCtx%Elem(iE)%BF(iDoF, iGauss) )
               flops = flops + 4.0
            End Do
            func = func + AppCtx%Elem(iE)%Gauss_C(iGauss) * ( (Strain_Elem .DotP. Strain_Elem) * 0.5_Kr - F_Elem * X_Elem )
            flops = flops + 5.0
         End Do
         Call SectionRealUpdateClosure(Gradient_Sec, AppCtx%MeshTopology%Mesh, iE-1, Gradient_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iEloc

      DeAllocate(Gradient_Loc)
      DeAllocate(F_Loc)
      DeAllocate(X_Loc)
      Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyBlock_Event, iErr); CHKERRQ(iErr)
   End Subroutine FormFunctionAndGradientBlock
#endif

#if defined WITH_TAO
   Subroutine SetupBCBounds(AppCtx, iErr)
      Type(AppCtx_Type)                            :: AppCtx

      PetscInt                                     :: i, iErr
      PetscInt, Dimension(:), Pointer              :: BCFlagPtr
      PetscReal                                    :: Val
      PetscReal, Dimension(:), Pointer             :: UPtr, UplusPtr, UminusPtr
      Type(SectionReal)                            :: UplusSec, UminusSec

      Val = 1.0D20
      Call VecSet(AppCtx%Uplus, Val, iErr); CHKERRQ(iErr)
      Val = -1.0D20
      Call VecSet(AppCtx%Uminus, Val, iErr); CHKERRQ(iErr)
      
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'Uplus',  1, UplusSec , iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'Uminus', 1, UminusSec, iErr); CHKERRQ(iErr)

      Call SectionRealToVec(UplusSec,  AppCtx%Scatter, SCATTER_REVERSE, AppCtx%Uplus, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(UminusSec, AppCtx%Scatter, SCATTER_REVERSE, AppCtx%Uminus, iErr); CHKERRQ(iErr)
      
      Allocate(UplusPtr(1))
      Allocate(UminusPtr(1))
      Do i = 1, AppCtx%MeshTopology%Num_Verts
         Call SectionIntRestrict(AppCtx%BCFlag, AppCtx%MeshTopology%Num_Elems+i-1, BCFlagPtr, iErr);
         If ( BCFlagPtr(1) /= 0 ) Then
            Call SectionRealRestrict(AppCtx%U,     AppCtx%MeshTopology%Num_Elems+i-1, UPtr,      iErr);
            UplusPtr  = UPtr
            UminusPtr = UPtr
            Call SectionRealRestore(AppCtx%U,     AppCtx%MeshTopology%Num_Elems+i-1, UPtr,      iErr);
            Call SectionRealUpdate(UplusSec,  AppCtx%MeshTopology%Num_Elems+i-1, UplusPtr,  INSERT_VALUES, iErr);
            Call SectionRealUpdate(UminusSec, AppCtx%MeshTopology%Num_Elems+i-1, UminusPtr, INSERT_VALUES, iErr);
         End If
         Call SectionIntRestore(AppCtx%BCFlag, AppCtx%MeshTopology%Num_Elems+i-1, BCFlagPtr, iErr);
      End Do
      DeAllocate(UplusPtr)
      DeAllocate(UminusPtr)

      Call SectionRealToVec(UplusSec,  AppCtx%Scatter, SCATTER_FORWARD, AppCtx%Uplus, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(UminusSec, AppCtx%Scatter, SCATTER_FORWARD, AppCtx%Uminus, iErr); CHKERRQ(iErr)
      Call TaoAppSetVariableBounds(AppCtx%TaoAppU, AppCtx%Uminus, AppCtx%Uplus, iErr); CHKERRQ(iErr)
      
      Call SectionRealDestroy(UplusSec,  iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(UminusSec, iErr); CHKERRQ(iErr)
   End Subroutine SetupBCBounds   
#endif

   Subroutine SimplePoissonFinalize(AppCtx)   
      Type(AppCtx_Type)                            :: AppCtx

      PetscInt                                     :: iErr
      Character(len=MEF90_MXSTRLEN)                :: filename

      Call SectionRealDestroy(AppCtx%U, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(AppCtx%GradU, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(AppCtx%F, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(AppCtx%RHS, iErr); CHKERRQ(iErr)
      Call VecScatterDestroy(AppCtx%Scatter, iErr); CHKERRQ(iErr)
      Call SectionIntDestroy(AppCtx%BCFlag, iErr); CHKERRQ(iErr)
      Call MatDestroy(AppCtx%K, iErr); CHKERRQ(iErr)
      Call VecDestroy(AppCtx%U_Vec, iErr); CHKERRQ(iErr)
      If (AppCtx%AppParam%Use_Tao) Then
#if defined WITH_TAO
         Call TaoApplicationDestroy(AppCtx%taoappU, iErr); CHKERRQ(iErr)
         Call TaoDestroy(AppCtx%taoU, iErr); CHKERRQ(iErr)
         Call VecDestroy(AppCtx%Uplus, iErr); CHKERRQ(iErr)
         Call VecDestroy(AppCtx%Uminus, iErr); CHKERRQ(iErr)
#endif
      Else
         Call KSPDestroy(AppCtx%KSP, iErr); CHKERRQ(iErr)
      End If
      Call MeshDestroy(AppCtx%MeshTopology%Mesh, iErr); CHKERRQ(ierr)
      
      If (AppCtx%AppParam%verbose > 1) Then
         Call PetscViewerFlush(AppCtx%AppParam%MyLogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerFlush(AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerDestroy(AppCtx%AppParam%MyLogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerDestroy(AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
      End If
      Write(filename, 103) Trim(AppCtx%AppParam%prefix)
      Call PetscLogPrintSummary(PETSC_COMM_WORLD, filename, iErr); CHKERRQ(iErr)
#if defined WITH_TAO
      Call TaoFinalize(iErr); CHKERRQ(iErr)
#endif
      Call MEF90_Finalize()
103 Format(A,'.log.txt')
   End Subroutine SimplePoissonFinalize

#if defined PB_2D
End Module m_SimplePoisson2D
#elif defined PB_3D
End Module m_SimplePoisson3D
#endif
