#if defined PB_2D
Module m_Poisson2D
#elif defined PB_3D
Module m_Poisson3D
#endif

#include "finclude/petscdef.h"

   Use m_MEF90
   
   Implicit NONE   

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
      PetscBool                                    :: Restart
      PetscInt                                     :: Verbose
      PetscInt                                     :: TestCase
      Character(len=MEF90_MXSTRLEN)                :: prefix
      Type(PetscViewer)                            :: LogViewer,MyLogViewer
   End Type AppParam_Type


   Type Heat_AppCtx_Type
      Type (dm)                                    :: mesh
      Type (EXO_Type)                              :: EXO
#if defined PB_2D

      Type(Element2D_Scal),Dimension(:),Pointer    :: Elem
#elif defined PB_3D
      Type(Element3D_Scal),Dimension(:),Pointer    :: Elem
#endif
      Type(IS)                                     :: CellSetGlobalIS
      Type(IS)                                     :: VertexSetGlobalIS
      Type(Elem_Type),Dimension(:),Pointer         :: cellSet
      !!! Since everything is done block after block, it is simpler to have one set of elements per block
      Type(SectionReal)                            :: GradU
      PetscReal                                    :: ElasticEnergy
      Type(Field)                                  :: U
      Type(Field)                                  :: UBC
      Type(Field)                                  :: F
      PetscReal                                    :: ExtForcesWork
      PetscReal                                    :: TotalEnergy
      Type(Flag)                                   :: BCFlag
      Type(Mat)                                    :: K
      Type(Mat)                                    :: M ! Mass matrix 
      Type(Mat)                                    :: Jac ! jacobian
      Type(Field)                                  :: RHS
      Type(KSP)                                    :: KSP
      Type(PC)                                     :: PC
      Type(LogInfo_Type)                           :: LogInfo
      Type(AppParam_Type)                          :: AppParam
   !For TS
      Type(TS)                                     :: TS
      PetscInt                                     :: maxsteps,NumSteps
      PetscReal                                    :: maxtime
      PetscInt                                     :: VertVar_Temperature 
      PetscReal,Dimension(:),Pointer               :: Diff,B_Mensi
   End Type Heat_AppCtx_Type
   
   
Contains
#undef __FUNCT__
#define __FUNCT__ "SimplePoissonInit"
   Subroutine SimplePoissonInit(AppCtx)
      Type(Heat_AppCtx_Type)                       :: AppCtx
      
      PetscInt                                     :: iErr
      PetscInt                                     :: iBlk,iDoF      
      PetscBool                                    :: HasPrefix,Flag
      PetscInt,Dimension(:),Pointer                :: TmpFlag
      
      PetscReal,Dimension(:),Pointer               :: ValPtr
      PetscReal,Dimension(:,:),Pointer             :: Coords
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer,filename   
      Type(DM)                                     :: Tmp_Mesh
      PetscReal                                    :: Val
      PetscInt,Dimension(:),Pointer                :: SizeScal
      PetscInt                                     :: numCellSet,numVertexSet,numDim
      Type(IS)                                     :: setIS,vertexIS
      PetscInt,Dimension(:),Pointer                :: setID,vertexID
      PetscInt                                     :: set,vertex
      Type(SectionReal)                            :: coordSection

      Call MEF90_Initialize()
      AppCtx%AppParam%verbose = 0
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-verbose',AppCtx%AppParam%verbose,Flag,iErr);CHKERRQ(iErr)
      
      AppCtx%AppParam%Restart = PETSC_FALSE
      Call PetscOptionsGetBool(PETSC_NULL_CHARACTER,'-restart',AppCtx%AppParam%restart,Flag,iErr);CHKERRQ(iErr)
      
      Call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-p',      AppCtx%AppParam%prefix,HasPrefix,iErr);CHKERRQ(iErr)
      If (.NOT. HasPrefix) Then
         Call PetscPrintf(PETSC_COMM_WORLD,"No mesh prefix given\n",iErr)
         Call MEF90_Finalize()
         STOP
      End If

      AppCtx%AppParam%TestCase = 1
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-test',      AppCtx%AppParam%TestCase,Flag,iErr);CHKERRQ(iErr)
      
      Call InitLog(AppCtx)
      Call PetscLogStagePush(AppCtx%LogInfo%Setup_Stage,iErr);CHKERRQ(iErr)
      If (AppCtx%AppParam%verbose > 1) Then
         Write(filename,101) Trim(AppCtx%AppParam%prefix),MEF90_MyRank
         Call PetscViewerASCIIOpen(PETSC_COMM_SELF,filename,AppCtx%AppParam%MyLogViewer,iErr);CHKERRQ(iErr);  
         Write(IOBuffer,102) MEF90_MyRank,Trim(filename)
         Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,iErr);CHKERRQ(iErr)
         Call PetscSynchronizedFlush (PETSC_COMM_WORLD,iErr);CHKERRQ(iErr)
   
         Write(filename,103) Trim(AppCtx%AppParam%prefix)
         Call PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,AppCtx%AppParam%LogViewer,iErr);CHKERRQ(iErr);  
         Write(IOBuffer,104) Trim(filename)
         Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr);CHKERRQ(iErr)
      End If
   
101 Format(A,'-',I4.4,'.log')
102 Format('Output from processor ',I4.4,' redirected to file ',A,'\n')
103 Format(A,'.log')
104 Format('Collective output redirected to file ',A,'\n')

      AppCtx%EXO%Comm = PETSC_COMM_WORLD
      AppCtx%EXO%num_nsproperties = 0
      AppCtx%EXO%num_ebproperties = 0
      AppCtx%EXO%filename = Trim(AppCtx%AppParam%prefix)//'.gen'
      If (MEF90_MyRank == 0) Then
         AppCtx%EXO%exoid = EXOPEN(AppCtx%EXO%filename,EXWRIT,exo_cpu_ws,exo_io_ws,PETSC_NULL_INTEGER,ierr)
      End If
      !!! Read and partition the mesh
      If (MEF90_NumProcs == 1) Then
         Call PetscLogStagePush(AppCtx%LogInfo%MeshCreateExodus_Stage,iErr);CHKERRQ(iErr)
         Call DMMeshCreateExodusNG(PETSC_COMM_WORLD,AppCtx%EXO%exoid,AppCtx%mesh,ierr);CHKERRQ(iErr)
         Call PetscLogStagePop(iErr);CHKERRQ(iErr)
      Else
         Call PetscLogStagePush(AppCtx%LogInfo%MeshCreateExodus_Stage,iErr);CHKERRQ(iErr)
         Call DMMeshCreateExodusNG(PETSC_COMM_WORLD,AppCtx%EXO%exoid,Tmp_mesh,ierr);CHKERRQ(iErr)
         Call PetscLogStagePop(iErr);CHKERRQ(iErr)
      
         Call PetscLogStagePush(AppCtx%LogInfo%MeshDistribute_Stage,iErr);CHKERRQ(iErr)
         Call DMMeshDistribute(Tmp_mesh,PETSC_NULL_CHARACTER,AppCtx%mesh,ierr);CHKERRQ(iErr)
         Call DMDestroy(Tmp_mesh,ierr);CHKERRQ(iErr)
         Call PetscLogStagePop(iErr);CHKERRQ(iErr)
      End If

      Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage,iErr);CHKERRQ(iErr)

      
      Call DMMeshGetLabelIdIS(AppCtx%mesh,'Cell Sets',AppCtx%CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,AppCtx%CellSetGlobalIS)
      Call DMMeshGetLabelIdIS(AppCtx%mesh,'Vertex Sets',AppCtx%VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,AppCtx%VertexSetGlobalIS)
      
      Call DMMeshGetStratumSize(AppCtx%mesh,"height",0,numCells,ierr);CHKERRQ(ierr)
      Allocate(AppCtx%Elem(numCells))

      Call DMMeshGetLabelIdIS(AppCtx%mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)

      Allocate(AppCtx%ElementType(size(setID)))
      !!! Sets the type of elements for each block
      !!! NO THis is a GLOBAL property and the size should be that of the GLOBAL # cell sets
      !!!
#if defined PB_2D
      AppCtx%ElementType(:) = MEF90_P1_Lagrange_2D_Scal
#elif defined PB_3D
      AppCtx%ElementType(:) = MEF90_P1_Lagrange_2D_Scal
#endif

      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call DMMeshGetStratumIS(AppCtx%mesh,'Cell Sets',setID(set),cellIS,ierr); CHKERRQ(iErr)
#if defined PB_2D
         Call ElementInit(AppCtx%mesh,cellIS,AppCtx%Elem,2,MEF90_P1_Lagrange_2D_Scal)
#elif defined PB_3D
         Call ElementInit(AppCtx%mesh,cellIS,AppCtx%Elem,2,MEF90_P1_Lagrange_3D_Scal)
#endif
         Call ISDestroy(cellIS)
      End Do    
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   
      !!! Allocate the Section for U and F
      Allocate(SizeScal(1))
      SizeScal=1

      Call FieldCreateVertex(AppCtx%U,    'U',        AppCtx%mesh,SizeScal)
      Call FieldCreateVertex(AppCtx%F,    'F',        AppCtx%mesh,SizeScal)
      Call FieldCreateVertex(AppCtx%RHS,  'RHS',      AppCtx%mesh,SizeScal)


      !!! Allocate and initialize the Section for the flag
      Call FlagCreateVertex(AppCtx%BCFlag,'BC',AppCtx%mesh,SizeScal)
!could we use    Call SectionIntAddNsProperty()  yes 
      
      Allocate(TmpFlag(1))
      TmpFlag = 1
      Call DMMeshGetLabelIdIS(AppCtx%mesh,'Vertex Sets',setIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call DMMeshGetStratumIS(AppCtx%mesh,'Vertex Sets',setID(set),vertexIS,ierr); CHKERRQ(iErr)
         Call ISGetIndicesF90(vertexIS,vertexID,iErr);CHKERRQ(iErr)
         Do vertex = 1, size(vertexID)
            Call SectionIntUpdate(AppCtx%BCFlag%Sec,vertexID(vertex),TmpFlag,INSERT_VALUES,iErr);CHKERRQ(iErr)
         End Do
         Call ISRestoreIndicesF90(vertexIS,vertexID,iErr);CHKERRQ(iErr)
         Call ISDestroy(vertexIS,ierr);CHKERRQ(ierr)
      End Do
      DeAllocate(TmpFlag)
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)

      !!! Initialize the matrix and vector for the linear system
      Call DMMeshSetMaxDof(AppCtx%mesh,1,iErr);CHKERRQ(iErr) 
      Call DMMeshCreateMatrix(AppCtx%mesh,AppCtx%U%Sec,MATMPIAIJ,AppCtx%K,iErr);CHKERRQ(iErr)
     

      Call KSPCreate(PETSC_COMM_WORLD,AppCtx%KSP,iErr);CHKERRQ(iErr)
      Call KSPSetOperators(AppCtx%KSP,AppCtx%K,AppCtx%K,SAME_NONZERO_PATTERN,iErr);CHKERRQ(iErr)
      Call KSPSetInitialGuessNonzero(AppCtx%KSP,PETSC_TRUE,iErr);CHKERRQ(iErr)
      Call KSPGetPC(AppCtx%KSP,AppCtx%PC,iErr);CHKERRQ(iErr)
      Call PCSetType(AppCtx%PC,PCBJACOBI,iErr);CHKERRQ(iErr)
      Call KSPSetType(AppCtx%KSP,KSPCG,iErr);CHKERRQ(iErr)
      Call KSPAppendOptionsPrefix(AppCtx%KSP,"U_",iErr);CHKERRQ(iErr)
      Call KSPSetFromOptions(AppCtx%KSP,iErr);CHKERRQ(iErr)
      
      Call PetscLogStagePop(iErr);CHKERRQ(iErr)
      
      
      !!! Read Force and BC from Data file or reformat it
      
      !!!
      !!! Create local vectors, then call 
      !!! VecLoadExodusVertex(dm,res,IOcomm,exoidout,step,offset)
      !!!
      If (AppCtx%AppParam%Restart) Then
         Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage,iErr);CHKERRQ(iErr)
         Call VecLoadExodusVertex(AppCtx%mesh,AppCtx%U%LocalVec,PETSC_COMM_WORLD,AppCtx%EXO%exoid,1,1,ierr);CHKERRQ(ierr)
         Call VecLoadExodusVertex(AppCtx%mesh,AppCtx%F%LocalVec,PETSC_COMM_WORLD,AppCtx%EXO%exoid,1,2,ierr);CHKERRQ(ierr)
         Call PetscLogStagePop(iErr);CHKERRQ(iErr)
      Else
         !!! Prepare and format the output mesh   
         Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage,iErr);CHKERRQ(iErr)
         !!! Check. We probably don't want to overwrite the input file
         Call EXOFormat_SimplePoisson(AppCtx)
         Call PetscLogStagePop(iErr);CHKERRQ(iErr)
         
         Select Case (AppCtx%AppParam%TestCase)
         Case(1)
            If (AppCtx%AppParam%verbose > 0) Then
               Write(IOBuffer,*) 'Setting U to 0 and F to 1\n'
               Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr);CHKERRQ(iErr)
            End If
            !!! U = 0,F=1
            Val = 1.0_Kr
            Call SectionRealSet(AppCtx%F,Val,iErr);CHKERRQ(iErr);
            Val = 0.0_Kr
            Call SectionRealSet(AppCtx%U,Val,iErr);CHKERRQ(iErr);
         Case(2)
            If (AppCtx%AppParam%verbose > 0) Then
               Write(IOBuffer,*) 'Solving Test Case 2: pure dirichlet problem,no force\n'
               Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr);CHKERRQ(iErr)
            End If


      !!! OLD WAY
! Test of non homogeneous Dirichlet BC
!            Allocate(ValPtr(1))
!            ValPtr = 0.0_Kr
!            Call SectionRealSet(AppCtx%F,0.0_Kr,iErr);CHKERRQ(iErr);
!            Call DMMeshGetCoordinatesF90(AppCtx%mesh,Coords,iErr);CHKERRQ(iErr)
!               
!            Do iDoF = 1,Size(Coords,1)
!#if defined PB_2D
!               ValPtr = (Coords(iDoF,1)-0.5_Kr)**2 + (Coords(iDoF,2)+0.5_Kr)**2
!               ValPtr = Coords(iDoF,1)**2 + Coords(iDoF,2)**2
!               ValPtr = Coords(iDoF,1)**3
!#elif defined PB_3D
!               ValPtr = (Coords(iDoF,1)-0.5_Kr)**2 + (Coords(iDoF,2)+0.5_Kr)**2 +  (Coords(iDoF,3)+.5_Kr)**2
!#endif
!               Call SectionRealUpdate(AppCtx%U%Sec,AppCtx%MeshTopology%Num_Elems+iDoF-1,ValPtr,INSERT_VALUES,iErr);CHKERRQ(iErr)
!            End Do
!            DeAllocate(ValPtr)
!            Call DMMeshRestoreCoordinatesF90(AppCtx%mesh,Coords,iErr);CHKERRQ(iErr)
      !!! NEW WAY
            Allocate(ValPtr(1))
            Call DMMeshGetSectionReal(AppCtx%mesh,'coordinates',coordSection,iErr); CHKERRQ(ierr)
            Call DMMeshGetStratumIS(AppCtx%mesh,'height',0,vertexIS,ierr);CHKERRQ(ierr)
            Call ISGetIndicesF90(vertexIS,vertexID,iErr);CHKERRQ(iErr)
            Do vertex = 1, size(vertexID)
               Call SectionRealRestrict(coordSection,vertexID(vertex),coord,ierr);CHKERRQ(ierr)
               ValPtr = 0
               Do c = 1,numDim
                  ValPtr = ValPtr + coord(c)**2
               End Do
               Call SectionRealUpdate(AppCtx%U%Sec,vertexID(vertex),ValPtr,INSERT_VALUES,iErr);CHKERRQ(iErr)
               Call SectionRealRestore(coordSection,vertexID(vertex),coord,ierr);CHKERRQ(ierr)
            End Do
            Call ISRestoreIndicesF90(vertexIS,vertexID,iErr);CHKERRQ(iErr)
            Call ISDestroy(vertexIS,ierr);CHKERRQ(ierr)
            Call SectionRealDestroy(coordSection,ierr);CHKERRQ(ierr)
            DeAllocate(ValPtr)
!         Case(3)
!            If (AppCtx%AppParam%verbose > 0) Then
!               Write(IOBuffer,*) 'Solving Test Case 3: F=sgn(x) . sgn(y) [. sgn(z)] \n'
!               Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr);CHKERRQ(iErr)
!            End If
!
!            !!! Test of non homogeneous Dirichlet BC
!            Allocate(ValPtr(1))
!            ValPtr = 0.0_Kr
!            Call SectionRealSet(AppCtx%U,1.0_Kr,iErr);CHKERRQ(iErr);
!            Call DMMeshGetCoordinatesF90(AppCtx%mesh,Coords,iErr);CHKERRQ(iErr)
!            
!            Do iDoF = 1,Size(Coords,1)
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
!               Call SectionRealUpdate(AppCtx%F%Sec,AppCtx%MeshTopology%Num_Elems+iDoF-1,ValPtr,INSERT_VALUES,iErr);CHKERRQ(iErr)
!            End Do
!            Call DMMeshRestoreCoordinatesF90(AppCtx%mesh,Coords,iErr);CHKERRQ(iErr)
         End Select            
      End If
      
      !!! Create the EXO case file
      Call Write_EXO_Case(AppCtx%AppParam%prefix,'%0.4d',MEF90_NumProcs)
      Call PetscLogStagePop(iErr);CHKERRQ(iErr)
   End Subroutine SimplePoissonInit


   
#undef __FUNCT__
#define __FUNCT__ "ExoFormat_SimplePoisson"
   Subroutine EXOFormat_SimplePoisson(EXO)
      Type(EXO_Type)                               :: EXO
      PetscInt                                     :: iErr
   
      Call EXPVP (EXO%exoid,'g',3,iErr)
      Call EXPVAN(EXO%exoid,'g',3,(/'Elastic Energy ','Ext Forces work','Total Energy   '/),iErr)
      Call EXPVP (EXO%exoid,'n',2,iErr)
      Call EXPVAN(EXO%exoid,'n',2,(/'U','F'/),iErr)
#if defined PB_2D
      Call EXPVP (EXO%exoid,'e',2,iErr)
      Call EXPVAN(EXO%exoid,'e',2,(/'Grad U_X','Grad U_Y'/),iErr)
#elif defined PB_3D
      Call EXPVP (EXO%exoid,'e',3,iErr)
      Call EXPVAN(EXO%exoid,'e',3,(/'Grad U_X','Grad U_Y','Grad U_Z'/),iErr)
#endif
      Call EXPTIM(EXO%exoid,1,1.0_Kr,iErr)
   End Subroutine EXOFormat_SimplePoisson
   
#undef __FUNCT__
#define __FUNCT__ "InitLog"
   Subroutine InitLog(AppCtx)
      Type(Heat_AppCtx_Type)                       :: AppCtx
      PetscInt                                     :: iErr
      
      Call PetscLogEventRegister('MatAssembly Block',0,AppCtx%LogInfo%MatAssemblyBlock_Event,ierr);CHKERRQ(ierr)
      Call PetscLogEventRegister('RHSAssembly Block',0,AppCtx%LogInfo%RHSAssemblyBlock_Event,ierr);CHKERRQ(ierr)
      Call PetscLogEventRegister('Post Processing',  0,AppCtx%LogInfo%PostProc_Event,        ierr);CHKERRQ(ierr)

      Call PetscLogStageRegister("MeshCreateExodus",AppCtx%LogInfo%MeshCreateExodus_Stage,iErr)
      Call PetscLogStageRegister("MeshDistribute",  AppCtx%LogInfo%MeshDistribute_Stage,  iErr)
      Call PetscLogStageRegister("IO Stage",        AppCtx%LogInfo%IO_Stage,              iErr)
      Call PetscLogStageRegister("Setup",           AppCtx%LogInfo%Setup_Stage,           iErr)
      Call PetscLogStageRegister("Mat Assembly",    AppCtx%LogInfo%MatAssembly_Stage,     iErr)
      Call PetscLogStageRegister("RHS Assembly",    AppCtx%LogInfo%RHSAssembly_Stage,     iErr)
      Call PetscLogStageRegister("KSP Solve",       AppCtx%LogInfo%KSPSolve_Stage,        iErr)
      Call PetscLogStageRegister("Post Proc",       AppCtx%LogInfo%PostProc_Stage,        iErr)
   End Subroutine InitLog
   
#undef __FUNCT__
#define __FUNCT__ "Solve"
   Subroutine Solve(AppCtx)
      Type(Heat_AppCtx_Type)                       :: AppCtx
      
      PetscInt                                     :: iErr
      KSPConvergedReason                           :: KSPreason
      PetscInt                                     :: KSPNumIter
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
      
      Call PetscLogStagePush(AppCtx%LogInfo%KSPSolve_Stage,iErr);CHKERRQ(iErr)
      Call SectionRealToVec(AppCtx%U%Sec,AppCtx%U%Scatter,SCATTER_FORWARD,AppCtx%U%Vec,iErr);CHKERRQ(iErr)

      Call DMMeshCreateVector(AppCtx%mesh,AppCtx%RHS,AppCtx%RHS%Vec,iErr);CHKERRQ(iErr)
      Call SectionRealToVec(AppCtx%RHS%Sec,AppCtx%RHS%Scatter,SCATTER_FORWARD,AppCtx%RHS%Vec,iErr);CHKERRQ(iErr)

      Call KSPSolve(AppCtx%KSP,AppCtx%RHS%Vec,AppCtx%U%Vec,iErr);CHKERRQ(iErr)
         
      Call KSPGetIterationNumber(AppCtx%KSP,KSPNumIter,iErr);CHKERRQ(iErr)
      Call KSPGetConvergedReason(AppCtx%KSP,KSPreason,iErr);CHKERRQ(iErr)
      Write(IOBuffer,100) KSPNumIter,KSPreason
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,iErr);CHKERRQ(iErr)
      Call SectionRealToVec(AppCtx%U%Sec,AppCtx%U%Scatter,SCATTER_REVERSE,AppCtx%U%Vec,ierr);CHKERRQ(ierr)
 
      Call PetscLogStagePop(iErr);CHKERRQ(iErr)
100 Format('KSP converged in ',I5,' iterations. KSPConvergedReason is ',I2,'\n')
   End Subroutine Solve
   
 
#undef __FUNCT__
#define __FUNCT__ "MatAssembly"
   Subroutine HeatMatAssembly(AppCtx)   
      Type(Heat_AppCtx_Type)                       :: AppCtx
      PetscInt                                     :: iBlk,iErr
      
      Type(IS)                                     :: setIS
      PetscInt,Dimension(:),Pointer                :: setID
      PetscInt                                     :: set
      
      Call MatInsertVertexBoundaryValues(AppCtx%K,AppCtx%U,AppCtx%BCFlag,AppCtx%mesh)
      Call MatAssemblyBegin(AppCtx%K,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
      Call MatAssemblyEnd  (AppCtx%K,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)

!      Call PetscLogStagePush(AppCtx%LogInfo%MatAssembly_Stage,iErr);CHKERRQ(iErr)
      Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call MatAssemblyBlock(set,AppCtx)
         !!! set or setID(set)?
      End Do
      Call ISRestoreIndices(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      Call MatAssemblyBegin(AppCtx%K,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
      Call MatAssemblyEnd  (AppCtx%K,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)

!      Call PetscLogStagePop(iErr);CHKERRQ(iErr)
   End Subroutine HeatMatAssembly
      
   Subroutine MatAssemblyBlock(iBlk,AppCtx)
      Type(Heat_AppCtx_Type)                       :: AppCtx
      PetscInt                                     :: iBlk
      
      PetscInt                                     :: iE,iELoc,iErr,i
      PetscReal,Dimension(:,:),Pointer             :: MatElem
      PetscInt,Dimension(:),Pointer                :: BCFlag
      PetscInt                                     :: iDoF1,iDoF2,iGauss
      PetscInt                                     :: NumDoFScal
      PetscLogDouble                               :: flops = 0
      PetscReal                                    :: lDiff
      PetscReal,Dimension(:),Pointer               :: T_Loc
      PetscReal                                    :: T_Elem
!      Call PetscLogEventBegin(AppCtx%LogInfo%MatAssemblyBlock_Event,iErr);CHKERRQ(iErr)
     
     !!! iBlk needs to be the GLOBAL set ID so that I can recover the proper material properties
     !!! get Diff as a PETScOption
     !!! 
      NumDoFScal = MeshTopology%Elem_Blk(iBlk)%Num_DoF
      Allocate(MatElem(MeshTopology%Elem_Blk(iBlk)%Num_DoF,MeshTopology%Elem_Blk(iBlk)%Num_DoF))
      Allocate(BCFlag(MeshTopology%Elem_Blk(iBlk)%Num_DoF))
      Allocate(T_Loc(NumDoFScal))

      Do_iELoc: Do iELoc = 1,MeshTopology%Elem_Blk(iBlk)%Num_Elems
         i = MeshTopology%Elem_Blk(iBlk)%ID
         iE = MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         MatElem = 0.0_Kr
         BCFlag = 0
         Call SectionIntRestrictClosure(AppCtx%BCFlag%Sec,mesh,iE-1,MeshTopology%Elem_Blk(iBlk)%Num_DoF,BCFlag,iErr);CHKERRQ(ierr)
         Do iGauss = 1,size(AppCtx%Elem(iE)%Gauss_C)
            Select Case(AppCtx%AppParam%TestCase)
            Case(2)
               lDiff = AppCtx%Diff(i) 
            Case(3)
               ! Mensi Law D(C) = A exp (B*C) 
               ! 1988 Mensi-Acker-Attolou Mater Struct
               ! C : water content, A and B material parameters
               Call SectionRealRestrictClosure(AppCtx%U%Sec,mesh, iE-1,NumDoFScal,T_Loc,iErr);CHKERRQ(ierr)
               T_Elem = 0.0_Kr
               Do iDoF1 = 1,NumDoFScal
                  T_Elem = T_Elem + AppCtx%Elem(iE)%BF(iDoF1,iGauss) * T_Loc(iDoF1)
               End DO
               lDiff = AppCtx%Diff(i)*exp(AppCtx%B_Mensi(i)*T_Elem) 
            Case Default
               ldiff =1 
            End Select
            Do iDoF1 = 1,MeshTopology%Elem_Blk(iBlk)%Num_DoF
               If (BCFlag(iDoF1) == 0) Then
                  Do iDoF2 = 1,MeshTopology%Elem_Blk(iBlk)%Num_DoF
                    ! MatElem(iDoF1,iDoF1) = 1./2.
                     MatElem(iDoF2,iDoF1) = MatElem(iDoF2,iDoF1) + lDiff * AppCtx%Elem(iE)%Gauss_C(iGauss) * ( AppCtx%Elem(iE)%Grad_BF(iDoF1,iGauss) .DotP. AppCtx%Elem(iE)%Grad_BF(iDoF2,iGauss) )
                     flops = flops + 1
                  End Do
               End If
            End Do
         End Do
         Call DMMeshAssembleMatrix(AppCtx%K,mesh,AppCtx%U%Sec,iE-1,MatElem,ADD_VALUES,iErr);CHKERRQ(iErr)
      End Do Do_iELoc
   
      Call PetscLogFlops(flops,iErr);CHKERRQ(iErr)
      DeAllocate(BCFlag)
      DeAllocate(MatElem)
!      Call PetscLogEventEnd(AppCtx%LogInfo%MatAssemblyBlock_Event,iErr);CHKERRQ(iErr)
   End Subroutine MatAssemblyBlock


#undef __FUNCT__
#define __FUNCT__ "RHSAssembly"
   Subroutine RHSAssembly(AppCtx,MeshTopology,MyExo)
      Type(Heat_AppCtx_Type)                       :: AppCtx
      Type (MeshTopology_Type)                     :: MeshTopology
      Type (EXO_Type)                              :: MyEXO

      PetscInt                                     :: iErr
      PetscInt                                     :: iBlk
!      Call PetscLogStagePush(AppCtx%LogInfo%RHSAssembly_Stage,iErr);CHKERRQ(iErr)

      Call SectionRealZero(AppCtx%RHS%Sec,iErr);CHKERRQ(iErr)
      
      Do_iBlk: Do iBlk = 1,MeshTopology%Num_Elem_Blks
         Call RHSAssemblyBlock(iBlk,AppCtx,MeshTopology)
      End Do Do_iBlk

      Call SectionRealComplete(AppCtx%RHS%Sec,iErr);CHKERRQ(iErr)

      !!! Set Dirichlet Boundary Values
!Suppose that loading is contant (replace second to last by AppCtx%Timestep otherwise)
      Call Read_EXO_Result_Vertex(MyEXO,MeshTopology,AppCtx%VertVar_Temperature,1,AppCtx%UBC)
      Call FieldInsertVertexBoundaryValues(AppCtx%RHS,AppCtx%UBC,AppCtx%BCFlag,MeshTopology)

      Call SectionRealToVec(AppCtx%RHS%Sec,AppCtx%RHS%Scatter,SCATTER_FORWARD,AppCtx%RHS%Vec,iErr);CHKERRQ(iErr)
      Call FieldInsertVertexBoundaryValues(AppCtx%RHS,AppCtx%UBC,AppCtx%BCFlag,MeshTopology)
      !!! VERY important! This is the equivalent of a ghost update
!      Call PetscLogStagePop(iErr);CHKERRQ(iErr)
   End Subroutine RHSAssembly

#undef __FUNCT__
#define __FUNCT__ "RHSAssemblyBlock"
   Subroutine RHSAssemblyBlock(iBlk,AppCtx,MeshTopology)
      PetscInt                                     :: iBlk
      Type(Heat_AppCtx_Type)                       :: AppCtx
      Type (MeshTopology_Type)                     :: MeshTopology

      PetscInt                                     :: iErr
      PetscInt                                     :: iBlkId
      PetscInt                                     :: Num_DoF,iDoF
      PetscInt                                     :: iEloc,iE,iGauss
      PetscInt,Dimension(:),Pointer                :: BCFlag_Loc
      PetscReal,Dimension(:),Pointer               :: RHS_Loc,F_Loc
      PetscReal                                    :: F_Elem
      PetscLogDouble                               :: flops 
      
      Call PetscLogEventBegin(AppCtx%LogInfo%RHSAssemblyBlock_Event,iErr);CHKERRQ(iErr)
      flops = 0.0

      Num_DoF = MeshTopology%Elem_Blk(iBlk)%Num_DoF
      Allocate(F_Loc(Num_DoF))
      Allocate(RHS_Loc(Num_DoF))
      Allocate(BCFlag_Loc(Num_DoF))

      iBlkID = MeshTopology%Elem_Blk(iBlk)%ID
      Do_iEloc: Do iELoc = 1,MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         RHS_Loc = 0.0_Kr
         Call SectionRealRestrictClosure(AppCtx%F%Sec,mesh,iE-1,Num_DoF,F_Loc,iErr);CHKERRQ(ierr)
         Call SectionIntRestrictClosure(AppCtx%BCFlag%Sec,mesh,iE-1,Num_DoF,BCFlag_Loc,iErr);CHKERRQ(ierr)
         Do iGauss = 1,Size(AppCtx%Elem(iE)%Gauss_C)
            F_Elem = 0.0_Kr
            Do iDoF = 1,Num_DoF
               F_Elem = F_Elem + AppCtx%Elem(iE)%BF(iDoF,iGauss) * F_Loc(iDoF)
               flops = flops + 2.0
            End Do
            Do iDoF = 1,Num_DoF
               If (BCFlag_Loc(iDoF) == 0) Then
                  RHS_Loc(iDoF) = RHS_Loc(iDoF) + AppCtx%Elem(iE)%Gauss_C(iGauss) * ( F_Elem * AppCtx%Elem(iE)%BF(iDoF,iGauss) )
                  flops = flops + 3.0
               End If
            End Do
         End Do
         Call SectionRealUpdateClosure(AppCtx%RHS%Sec,mesh,iE-1,RHS_Loc,ADD_VALUES,iErr);CHKERRQ(iErr)
      End Do Do_iEloc

      DeAllocate(BCFlag_Loc)
      DeAllocate(RHS_Loc)
      DeAllocate(F_Loc)
      Call PetscLogFlops(flops ,iErr);CHKERRQ(iErr)
      Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyBlock_Event,iErr);CHKERRQ(iErr)
   End Subroutine RHSAssemblyBlock

#undef __FUNCT__
#define __FUNCT__ "ComputeEnergy"
   Subroutine ComputeEnergy(AppCtx)
      Type(Heat_AppCtx_Type)                       :: AppCtx
      
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoF,NumGauss
      PetscReal,Dimension(:),Pointer               :: F,U
      PetscInt                                     :: iBlk,iELoc,iE
      PetscInt                                     :: iDoF,iGauss
#if defined PB_2D
      Type(Vect2D)                                 :: Strain_Elem,Stress_Elem      
#elif defined PB_3D
      Type(Vect3D)                                 :: Strain_Elem,Stress_Elem      
#endif
      PetscReal                                    :: F_Elem,U_Elem
      PetscReal                                    :: MyElasticEnergy,MyExtForcesWork
      PetscLogDouble                               :: flops

      Call PetscLogStagePush(AppCtx%LogInfo%PostProc_Stage,iErr);CHKERRQ(iErr)
      Call PetscLogEventBegin(AppCtx%LogInfo%PostProc_Event,iErr);CHKERRQ(iErr)
      
      MyElasticEnergy = 0.0_Kr
      MyExtForcesWork = 0.0_Kr
      flops = 0.0_Kr
      Do_Elem_iBlk: Do iBlk = 1,AppCtx%MeshTopology%Num_Elem_Blks
         Do_Elem_iE: Do iELoc = 1,AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            NumDoF   = Size(AppCtx%Elem(iE)%BF,1)
            NumGauss = Size(AppCtx%Elem(iE)%BF,2)
            Allocate(F(NumDoF))
            Call SectionRealRestrictClosure(AppCtx%F%Sec,AppCtx%mesh,iE-1,NumDoF,F,iErr);CHKERRQ(ierr)
            Allocate(U(NumDoF))
            Call SectionRealRestrictClosure(AppCtx%U%Sec,AppCtx%mesh,iE-1,NumDoF,U,iErr);CHKERRQ(ierr)
            Do iGauss = 1,NumGauss
               Strain_Elem = 0.0_Kr
               Stress_Elem = 0.0_Kr
               F_Elem      = 0.0_Kr
               U_Elem      = 0.0_Kr
               Do iDoF = 1,NumDoF
                  Stress_Elem = Stress_Elem + AppCtx%Elem(iE)%Grad_BF(iDoF,iGauss) * U(iDoF)
                  Strain_Elem = Strain_Elem + AppCtx%Elem(iE)%Grad_BF(iDoF,iGauss) * U(iDoF)
                  F_Elem = F_Elem + AppCtx%Elem(iE)%BF(iDoF,iGauss) * F(iDoF)
                  U_Elem = U_Elem + AppCtx%Elem(iE)%BF(iDoF,iGauss) * U(iDoF)
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
      Call MPI_AllReduce(MyElasticEnergy,AppCtx%ElasticEnergy,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD,iErr);CHKERRQ(iErr)
      Call MPI_AllReduce(MyExtForcesWork,AppCtx%ExtForcesWork,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD,iErr);CHKERRQ(iErr)
      AppCtx%TotalEnergy = AppCtx%ElasticEnergy + AppCtx%ExtForcesWork
      Call PetscLogFlops(flops,iErr)
      Call PetscLogEventEnd  (AppCtx%LogInfo%PostProc_Event,iErr);CHKERRQ(iErr)
      Call PetscLogStagePop(iErr);CHKERRQ(iErr)
   End Subroutine ComputeEnergy
 
#undef __FUNCT__
#define __FUNCT__ "ComputeGradU"
   Subroutine ComputeGradU(AppCtx)
      Type(Heat_AppCtx_Type)                       :: AppCtx
      
      PetscInt                                     :: iErr
#if defined PB_2D
      Type(Vect2D)                                 :: Grad
#elif defined PB_3D
      Type(Vect3D)                                 :: Grad
#endif
      PetscReal                                    :: Vol
      PetscInt                                     :: NumDoF,NumGauss
      PetscReal,Dimension(:),Pointer               :: U
      PetscInt                                     :: iBlk,iELoc,iE
      PetscInt                                     :: iDoF,iGauss
      PetscReal,Dimension(:),Pointer               :: Grad_Ptr
      PetscLogDouble                               :: flops = 0

      Call PetscLogStagePush (AppCtx%LogInfo%PostProc_Stage,iErr);CHKERRQ(iErr)
      Call PetscLogEventBegin(AppCtx%LogInfo%PostProc_Event,iErr);CHKERRQ(iErr)
      Call DMMeshGetCellSectionReal(AppCtx%mesh,  'GradU',AppCtx%MeshTopology%Num_Dim,AppCtx%GradU,iErr);CHKERRQ(iErr)

      Allocate(Grad_Ptr(AppCtx%MeshTopology%Num_Dim))
      Do_Elem_iBlk: Do iBlk = 1,AppCtx%MeshTopology%Num_Elem_Blks
         Do_Elem_iE: Do iELoc = 1,AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            NumDoF   = Size(AppCtx%Elem(iE)%BF,1)
            NumGauss = Size(AppCtx%Elem(iE)%BF,2)
            Allocate(U(NumDoF))
            Call SectionRealRestrictClosure(AppCtx%U%Sec,AppCtx%mesh,iE-1,NumDoF,U,iErr);CHKERRQ(ierr)
            Grad = 0.0_Kr
            Vol  = 0.0_Kr
            Do iGauss = 1,NumGauss
               Do iDoF = 1,NumDoF
                  Grad = Grad + (AppCtx%Elem(iE)%Gauss_C(iGauss) * U(iDoF) ) * AppCtx%Elem(iE)%Grad_BF(iDoF,iGauss) 
                  Vol = Vol + AppCtx%Elem(iE)%Gauss_C(iGauss) * AppCtx%Elem(iE)%BF(iDoF,iGauss)
                  flops = flops + 3
               End Do
            End Do
            Grad = Grad / Vol
#if defined PB_2D
            Grad_Ptr = (/ Grad%X,Grad%Y /)
#elif defined PB_3D
            Grad_Ptr = (/ Grad%X,Grad%Y,Grad%Z /)
#endif
            Call SectionRealUpdateClosure(AppCtx%GradU,AppCtx%mesh,iE-1,Grad_Ptr,INSERT_VALUES,iErr)
            DeAllocate(U)
         End Do Do_Elem_iE
      End Do Do_Elem_iBlk
      DeAllocate(Grad_Ptr)
      Call PetscLogFlops(flops,iErr)
      Call PetscLogEventEnd  (AppCtx%LogInfo%PostProc_Event,iErr);CHKERRQ(iErr)
      Call PetscLogStagePop(iErr);CHKERRQ(iErr)
   End Subroutine ComputeGradU
 
#undef __FUNCT__
#define __FUNCT__ "SimplePoissonFinalize"
   Subroutine SimplePoissonFinalize(AppCtx)   
      Type(Heat_AppCtx_Type)                       :: AppCtx

      PetscInt                                     :: iErr
      Character(len=MEF90_MXSTRLEN)                :: filename
      Type(PetscViewer)                            :: LogViewer

      Call SectionRealDestroy(AppCtx%GradU,iErr);CHKERRQ(iErr)
      Call FieldDestroy(AppCtx%U)
      !Call FieldDestroy(AppCtx%GradU)
      Call FieldDestroy(AppCtx%F)
      Call FieldDestroy(AppCtx%RHS)

      Call SectionIntDestroy(AppCtx%BCFlag,iErr);CHKERRQ(iErr)
      Call MatDestroy(AppCtx%K,iErr);CHKERRQ(iErr)
      Call KSPDestroy(AppCtx%KSP,iErr);CHKERRQ(iErr)
      Call DMDestroy(AppCtx%mesh,iErr);CHKERRQ(ierr)
      
      If (AppCtx%AppParam%verbose > 1) Then
         Call PetscViewerFlush(AppCtx%AppParam%MyLogViewer,iErr);CHKERRQ(iErr)
         Call PetscViewerFlush(AppCtx%AppParam%LogViewer,iErr);CHKERRQ(iErr)
         Call PetscViewerDestroy(AppCtx%AppParam%MyLogViewer,iErr);CHKERRQ(iErr)
         Call PetscViewerDestroy(AppCtx%AppParam%LogViewer,iErr);CHKERRQ(iErr)
      End If
      Write(filename,103) Trim(AppCtx%AppParam%prefix)
103 Format(A,'-logsummary.txt')
      Call PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,LogViewer,ierr);CHKERRQ(ierr)
      Call PetscViewerASCIIPrintf(LogViewer,filename,ierr);CHKERRQ(ierr)
      Call PetscLogView(LogViewer,ierr);CHKERRQ(ierr)
      Call PetscViewerDestroy(LogViewer,ierr);CHKERRQ(ierr)
!       Call PetscLogPrintSummary(PETSC_COMM_WORLD,filename,iErr);CHKERRQ(iErr)
      Call MEF90_Finalize()
   End Subroutine SimplePoissonFinalize

#if defined PB_2D
End Module m_Poisson2D
#elif defined PB_3D
End Module m_Poisson3D
#endif
