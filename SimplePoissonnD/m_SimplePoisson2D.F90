#if defined PB_2D
Module m_SimplePoisson2D
#elif defined PB_3D
Module m_SimplePoisson3D
#endif

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscmeshdef.h"
#include "finclude/petscviewerdef.h"

   Use m_MEF90
   Use petsc
   Use petscvec
   Use petscmat
   Use petscksp
   Use petscmesh

   Implicit NONE   
   
   Type LogInfo_Type
      PetscLogStage               :: IO_Stage
      PetscLogStage               :: Distribute_Stage
      PetscLogStage               :: DataSetup_Stage
      PetscLogStage               :: MatAssembly_Stage
      PetscLogStage               :: RHSAssembly_Stage
      PetscLogStage               :: KSPSolve_Stage
      PetscLogStage               :: EnergyEval_Stage
      
      PetscLogEvent               :: MatAssemblyLocal_Event
      PetscLogEvent               :: RHSAssemblyLocal_Event
      PetscLogEvent               :: EnergyEval_Event
   End Type LogInfo_Type

   Type AppParam_Type
      PetscTruth                                   :: Verbose
      Character(len=MEF90_MXSTRLEN)                :: prefix
      Type(PetscViewer)                            :: LogViewer, MyLogViewer
   End Type AppParam_Type

   Type AppCtx_Type
      Type (MeshTopology_Info)                     :: MeshTopology
      Type (EXO_Info)                              :: EXO, MyEXO
#if defined PB_2D
      Type(Element2D_Scal), Dimension(:), Pointer  :: Elem
#elif defined PB_3D
      Type(Element3D_Scal), Dimension(:), Pointer  :: Elem
#endif
      Type(SectionReal)                            :: U
      Type(SectionReal)                            :: F
      PetscReal                                    :: Energy
      Type(VecScatter)                             :: Scatter
      Type(SectionInt)                             :: BCFlag
      Type(Mat)                                    :: K
      Type(Vec)                                    :: RHS
      Type(KSP)                                    :: KSP
      Type(PC)                                     :: PC
      Type(LogInfo_Type)                           :: LogInfo
      Type(AppParam_Type)                          :: AppParam
   End Type AppCtx_Type
   
   
Contains

   Subroutine SimplePoissonInit(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      PetscInt                                     :: iBlk, iDoF      
      PetscTruth                                   :: HasPrefix
      PetscInt, Dimension(:), Pointer              :: TmpFlag
      PetscInt                                     :: TmpPoint
      
      Type(SectionReal)                            :: CoordSection
      PetscReal, Dimension(:), Pointer             :: TmpCoords
      PetscReal, Dimension(:,:), Pointer           :: Coords
      PetscInt                                     :: iE, iELoc
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer, filename   
      PetscLogDouble                               :: TS, TF
      Type(Mesh)                                   :: Tmp_Mesh


      Call MEF90_Initialize()
      Call PetscOptionsHasName(PETSC_NULL_CHARACTER, '-verbose', AppCtx%AppParam%verbose, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p', AppCtx%AppParam%prefix, HasPrefix, iErr); CHKERRQ(iErr)
      
      Call InitLog(AppCtx)
      If (AppCtx%AppParam%verbose) Then
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
102 Format('Output from processor ', I4.4, ' redirected to file ', A, '\n'c)
103 Format(A,'.log')
104 Format('Collective output redirected to file ', A, '\n'c)


      !!! Read and partition the mesh
      Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
      AppCtx%EXO%Comm = PETSC_COMM_WORLD
      AppCtx%EXO%filename = Trim(AppCtx%AppParam%prefix)//'.gen'

      
      Call MeshCreateExodus(PETSC_COMM_WORLD, AppCtx%EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
      !!! reads exo file, stores all information in a Mesh
      Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)

      Call PetscLogStagePush(AppCtx%LogInfo%Distribute_Stage, iErr); CHKERRQ(iErr)
      Call MeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, AppCtx%MeshTopology%mesh, ierr); CHKERRQ(iErr)
      Call MeshDestroy(Tmp_mesh, ierr); CHKERRQ(iErr)
      Call PetscLogStagePop(AppCtx%LogInfo%Distribute_Stage, iErr); CHKERRQ(iErr)

      Call MeshTopologyReadEXO(AppCtx%MeshTopology, AppCtx%EXO)
         
      Call PetscLogStagePush(AppCtx%LogInfo%DataSetup_Stage, iErr); CHKERRQ(iErr)
      !!! Sets the type of elements for each block
      Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_Type = MEF90_P1_Lagrange
         Call Init_Elem_Blk_Info(AppCtx%MeshTopology%Elem_Blk(iBlk), AppCtx%MeshTopology%num_dim)
      End Do
      Call PetscLogStagePop(AppCtx%LogInfo%Distribute_Stage, iErr); CHKERRQ(iErr)
   
      !!! Allocate the elements
      Allocate(AppCtx%Elem(AppCtx%MeshTopology%Num_Elems))

      !!! Initialize the Basis Functions in each element
      Call MeshGetSectionReal(AppCtx%MeshTopology%mesh, 'coordinates', CoordSection, iErr); CHKERRQ(iErr)
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Allocate(TmpCoords(AppCtx%MeshTopology%Num_Dim * AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF))
         Allocate(Coords(AppCtx%MeshTopology%Num_Dim, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF))
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, CoordSection, iE-1, Size(TmpCoords), TmpCoords, iErr); CHKERRQ(iErr)
             Coords = Reshape(TmpCoords, (/AppCtx%MeshTopology%Num_Dim, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF /) )
            Call Init_Element(AppCtx%Elem(iE), Coords, 2, AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_Type)
         End Do Do_Elem_iE
         DeAllocate(TmpCoords)
         DeAllocate(Coords)
      End Do Do_Elem_iBlk
      Call SectionRealDestroy(CoordSection, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(AppCtx%LogInfo%DataSetup_Stage, iErr); CHKERRQ(iErr)

      Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
      !!! Prepare and format the output mesh
      !!! 1. Geometry
      AppCtx%MyEXO%comm = PETSC_COMM_SELF
      AppCtx%MyEXO%exoid = AppCtx%EXO%exoid
      Write(AppCtx%MyEXO%filename, 200) trim(AppCtx%AppParam%prefix), MEF90_MyRank
   200 Format(A, '-', I4.4, '.gen')
      AppCtx%MyEXO%title = trim(AppCtx%EXO%title)
      AppCtx%MyEXO%Num_QA = AppCtx%EXO%Num_QA
      Call Write_MeshTopologyGlobal(AppCtx%MeshTopology, AppCtx%MyEXO, PETSC_COMM_WORLD)
      
      !!! 2. Variables
      AppCtx%MyEXO%exoid = EXOPEN(AppCtx%MyEXO%filename, EXWRIT, exo_cpu_ws, exo_io_ws, exo_ver, iErr)
      Call EXPVP (AppCtx%MyEXO%exoid, 'g', 1, iErr)
      Call EXPVAN(AppCtx%MyEXO%exoid, 'g', 1, (/'Energy'/), iErr)
      Call EXPVP (AppCtx%MyEXO%exoid, 'n', 1, iErr)
      Call EXPVAN(AppCtx%MyEXO%exoid, 'n', 1, (/'U'/), iErr)
#if defined PB_2D
      Call EXPVP (AppCtx%MyEXO%exoid, 'e', 2, iErr)
      Call EXPVAN(AppCtx%MyEXO%exoid, 'e', 2, (/'Grad U_X', 'Grad U_Y'/), iErr)
#elif defined PB_3D
      Call EXPVP (AppCtx%MyEXO%exoid, 'e', 3, iErr)
      Call EXPVAN(AppCtx%MyEXO%exoid, 'e', 3, (/'Grad U_X', 'Grad U_Y', 'Grad U_Z'/), iErr)
#endif
      Call EXPTIM(AppCtx%MyEXO%exoid, 1, 1.0_Kr, iErr)

      Call EXCLOS(AppCtx%MyEXO%exoid, iErr)

      Call PetscLogStagePop(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)

      Call PetscLogStagePush(AppCtx%LogInfo%DataSetup_Stage, iErr); CHKERRQ(iErr)
      !!! Allocate the Section for U and F
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 1, AppCtx%U, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 1, AppCtx%F, iErr); CHKERRQ(iErr)
      Call MeshCreateGlobalScatter(AppCtx%MeshTopology%mesh, AppCtx%U, AppCtx%Scatter, iErr); CHKERRQ(iErr)
      
      !!! Allocate and initialize the Section for the flag
      Call MeshGetVertexSectionInt(AppCtx%MeshTopology%mesh, 1, AppCtx%BCFlag, iErr); CHKERRQ(iErr)
      Allocate(TmpFlag(1))
      Do iBlk = 1, AppCtx%MeshTopology%num_node_sets  
         Do iDoF = 1, AppCtx%MeshTopology%node_set(iBlk)%Num_Nodes     
            TmpPoint = AppCtx%MeshTopology%Num_Elems + AppCtx%MeshTopology%Node_Set(iBlk)%Node_ID(iDoF)-1
            TmpFlag = 1
            Call MeshUpdateClosureInt(AppCtx%MeshTopology%Mesh, AppCtx%BCFlag, TmpPoint, TmpFlag, iErr); CHKERRQ(iErr)
         End Do
      End Do
      DeAllocate(TmpFlag)

      !!! Initialize the matrix and vector for the linear system
      Call MeshSetMaxDof(AppCtx%MeshTopology%Mesh, 1, iErr); CHKERRQ(iErr) 
      !Max DoF per point is 1 (Should it be 3?)
      Call MeshCreateMatrix(AppCtx%MeshTopology%mesh, AppCtx%U, MATMPIAIJ, AppCtx%K, iErr); CHKERRQ(iErr)
      Call MeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%U, AppCtx%RHS, iErr); CHKERRQ(iErr)
      
      !!! Create the KSP and PCs
      Call KSPCreate(PETSC_COMM_WORLD, AppCtx%KSP, iErr); CHKERRQ(iErr)
      Call KSPSetOperators(AppCtx%KSP, AppCtx%K, AppCtx%K, SAME_NONZERO_PATTERN, iErr); CHKERRQ(iErr)
      
      Call KSPGetPC(AppCtx%KSP, AppCtx%PC, iErr); CHKERRQ(iErr)
      Call PCSetType(AppCtx%PC, PCBJACOBI, iErr); CHKERRQ(iErr)
      Call KSPSetType(AppCtx%KSP, KSPCG, iErr); CHKERRQ(iErr)
      
      Call KSPSetInitialGuessNonzero(AppCtx%KSP, PETSC_TRUE, iErr); CHKERRQ(iErr)
      
      Call KSPSetFromOptions(AppCtx%KSP, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(AppCtx%LogInfo%DataSetup_Stage, iErr); CHKERRQ(iErr)
      
      !!! Create the EXO cae file
      Call Write_EXO_Case(AppCtx%AppParam%prefix, '%0.4d', MEF90_NumProcs)
   End Subroutine SimplePoissonInit
   
   Subroutine InitLog(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr
      
      Call PetscLogEventRegister('MatAssembly Local', 0, AppCtx%LogInfo%MatAssemblyLocal_Event, ierr); CHKERRQ(ierr)
      Call PetscLogEventRegister('RHSAssembly Local', 0, AppCtx%LogInfo%RHSAssemblyLocal_Event, ierr); CHKERRQ(ierr)
      Call PetscLogEventRegister('Energy Eval',       0, AppCtx%LogInfo%EnergyEval_Event,       ierr); CHKERRQ(ierr)

      Call PetscLogStageRegister("IO Stage",          AppCtx%LogInfo%IO_Stage,          iErr)
      Call PetscLogStageRegister("Distribution",      AppCtx%LogInfo%Distribute_Stage,  iErr)
      Call PetscLogStageRegister("Objects Setup",     AppCtx%LogInfo%DataSetup_Stage,   iErr)
      Call PetscLogStageRegister("Mat Assembly",      AppCtx%LogInfo%MatAssembly_Stage, iErr)
      Call PetscLogStageRegister("RHS Assembly",      AppCtx%LogInfo%RHSAssembly_Stage, iErr)
      Call PetscLogStageRegister("KSP Solve",         AppCtx%LogInfo%KSPSolve_Stage,    iErr)
      Call PetscLogStageRegister("Energy Eval",       AppCtx%LogInfo%EnergyEval_Stage,  iErr)
   End Subroutine InitLog
   
   Subroutine Solve(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      KSPConvergedReason                           :: reason
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
      
      Call PetscLogStagePush(AppCtx%LogInfo%KSPSolve_Stage, iErr); CHKERRQ(iErr)
      Call KSPSolve(AppCtx%KSP, AppCtx%RHS, AppCtx%RHS, iErr); CHKERRQ(iErr)
      !!! Solve and store the solution in AppCtx%RHS
      
      Call SectionRealToVec(AppCtx%U, AppCtx%Scatter, SCATTER_REVERSE, AppCtx%RHS, ierr); CHKERRQ(ierr)
      !!! Scatter the solution from (Vec) AppCtx%RHS to (SectionReal) AppCtx%U
      
      Call KSPGetConvergedReason(AppCtx%KSP, reason, iErr); CHKERRQ(iErr)
      Write(IOBuffer, *) 'KSPGetConvergedReason returned ', reason, '\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      
      Call PetscLogStagePop (AppCtx%LogInfo%KSPSolve_Stage, iErr); CHKERRQ(iErr)
   End Subroutine Solve
   
      
   Subroutine MatAssembly(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iBlk, iE, iELoc, iErr
      PetscReal, Dimension(:,:), Pointer           :: MatElem
      
      Call PetscLogStagePush(AppCtx%LogInfo%MatAssembly_Stage, iErr); CHKERRQ(iErr)
      
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Allocate(MatElem(AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF))
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            Call MatAssemblyLocal(iE, AppCtx, MatElem)
            Call assembleMatrix(AppCtx%K, AppCtx%MeshTopology%mesh, AppCtx%U, iE-1, MatElem, ADD_VALUES, iErr); CHKERRQ(iErr)
         End Do Do_Elem_iE
         DeAllocate(MatElem)
      End Do Do_Elem_iBlk
      Call MatAssemblyBegin(AppCtx%K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (AppCtx%K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)

      Call PetscLogStagePop(AppCtx%LogInfo%MatAssembly_Stage, iErr); CHKERRQ(iErr)
   End Subroutine MatAssembly
   
   
   Subroutine MatAssemblyLocal(iE, AppCtx, MatElem)
      Type(AppCtx_Type)                            :: AppCtx
      PetscReal, Dimension(:,:), Pointer           :: MatElem 
      PetscInt                                     :: iE
   
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoF, NumGauss
      PetscInt, Dimension(:), Pointer              :: BCFlag
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      
      Call PetscLogEventBegin(AppCtx%LogInfo%MatAssemblyLocal_Event, iErr); CHKERRQ(iErr)
      
      MatElem  = 0.0_Kr
      NumDoF   = Size(AppCtx%Elem(iE)%BF,1)
      NumGauss = Size(AppCtx%Elem(iE)%BF,2)
      Allocate(BCFlag(NumDoF))
      Call MeshRestrictClosureInt(AppCtx%MeshTopology%mesh, AppCtx%BCFlag, iE-1, NumDoF, BCFlag, iErr); CHKERRQ(ierr)
      Do iGauss = 1, NumGauss
         Do iDoF1 = 1, NumDoF
            If (BCFlag(iDoF1) == 0) Then
               Do iDoF2 = 1, NumDoF
                  MatElem(iDoF2, iDoF1) = AppCtx%Elem(iE)%Gauss_C(iGauss) * AppCtx%Elem(iE)%Grad_BF(iDoF1, iGauss) .DotP. AppCtx%Elem(iE)%Grad_BF(iDoF2, iGauss) 
                  Call PetscLogFlops(AppCtx%MeshTopology%num_dim * (AppCtx%MeshTopology%num_dim-1) +1 , iErr);CHKERRQ(iErr)
                  !!! Is that right?
               End Do
            End If
         End Do
      End Do
      
      DeAllocate(BCFlag)
      Call PetscLogEventEnd(AppCtx%LogInfo%MatAssemblyLocal_Event, iErr); CHKERRQ(iErr)
   End Subroutine MatAssemblyLocal
   
   Subroutine RHSAssembly(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      Type(SectionReal)                            :: RHSSec
      PetscInt                                     :: iBlk, iE, iELoc
      PetscReal, Dimension(:), Pointer             :: RHSElem
      PetscInt                                     :: NumDoFPerVertex = 1
      
      !!! Hopefully one day we will use assembleVector instead of going through a section
      Call PetscLogStagePush(AppCtx%LogInfo%RHSAssembly_Stage, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, NumDoFPerVertex, RHSSec, iErr); CHKERRQ(iErr)
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Allocate(RHSElem(AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF))
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            Call RHSAssemblyLocal(iE, AppCtx, RHSElem)
            Call MeshUpdateAddClosure(AppCtx%MeshTopology%Mesh, RHSSec, iE-1, RHSElem, iErr)
         End Do Do_Elem_iE
         DeAllocate(RHSElem)
      End Do Do_Elem_iBlk
      Call SectionRealComplete(RHSSec, iErr); CHKERRQ(iErr)
      !!! VERY important! This is the equivalent of a ghost update
      Call SectionRealToVec(RHSSec, AppCtx%Scatter, SCATTER_FORWARD, AppCtx%RHS, ierr); CHKERRQ(ierr)
      Call SectionRealDestroy(RHSSec, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(AppCtx%LogInfo%RHSAssembly_Stage, iErr); CHKERRQ(iErr)
   End Subroutine RHSAssembly



   Subroutine RHSAssemblyLocal(iE, AppCtx, RHSElem)
      PetscInt                                     :: iE
      Type(AppCtx_Type)                            :: AppCtx
      PetscReal, Dimension(:), Pointer             :: RHSElem
      
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoF, NumGauss
      PetscInt, Dimension(:), Pointer              :: BCFlag
      PetscReal, Dimension(:), Pointer             :: F
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      PetscReal                                    :: TmpRHS
      
      Call PetscLogEventBegin(AppCtx%LogInfo%RHSAssemblyLocal_Event, iErr); CHKERRQ(iErr)
      
      RHSElem  = 0.0_Kr
      NumDoF   = Size(AppCtx%Elem(iE)%BF,1)
      NumGauss = Size(AppCtx%Elem(iE)%BF,2)
      Allocate(BCFlag(NumDoF))
      Call MeshRestrictClosureInt(AppCtx%MeshTopology%mesh, AppCtx%BCFlag, iE-1, NumDoF, BCFlag, iErr); CHKERRQ(ierr)
      Allocate(F(NumDoF))
      Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%F, iE-1, NumDoF, F, iErr); CHKERRQ(ierr)
      Do iGauss = 1, NumGauss
         TmpRHS = 0.0_Kr
         Do iDoF2 = 1, NumDoF
            TmpRHS = TmpRHS + AppCtx%Elem(iE)%BF(iDoF2, iGauss)! * F(iDoF2)
            Call PetscLogFlops(2 , iErr);CHKERRQ(iErr)
         End Do
         Do iDoF1 = 1, NumDoF
            If (BCFlag(iDoF1) == 0) Then
               RHSElem(iDoF1) = RHSElem(iDoF1) + AppCtx%Elem(iE)%Gauss_C(iGauss) * AppCtx%Elem(iE)%BF(iDoF1, iGauss) * TmpRHS
               Call PetscLogFlops(2 , iErr);CHKERRQ(iErr)
            End If
         End Do
      End Do
      
      DeAllocate(BCFlag)
      DeAllocate(F)
      Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyLocal_Event, iErr); CHKERRQ(iErr)
   End Subroutine RHSAssemblyLocal

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
      PetscReal                                    :: MyEnergy

      Call PetscLogStagePush(AppCtx%LogInfo%EnergyEval_Stage, iErr); CHKERRQ(iErr)
      Call PetscLogEventBegin(AppCtx%LogInfo%EnergyEval_Event, iErr); CHKERRQ(iErr)
      
      MyEnergy = 0.0_Kr
      
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            NumDoF   = Size(AppCtx%Elem(iE)%BF,1)
            NumGauss = Size(AppCtx%Elem(iE)%BF,2)
            Allocate(F(NumDoF))
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%F, iE-1, NumDoF, F, iErr); CHKERRQ(ierr)
            Allocate(U(NumDoF))
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%U, iE-1, NumDoF, U, iErr); CHKERRQ(ierr)
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
                  Call PetscLogFlops(0, iErr) !!! FIX THAT!
               End Do
               MyEnergy = MyEnergy + AppCtx%Elem(iE)%Gauss_C(iGauss) * ( (Strain_Elem .DotP. Strain_Elem) * 0.5_Kr - F_Elem * U_Elem)
            End Do
         End Do Do_Elem_iE
      End Do Do_Elem_iBlk
      Call PetscGlobalSum(MyEnergy, AppCtx%Energy, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
               
      Call PetscLogEventEnd  (AppCtx%LogInfo%EnergyEval_Event, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(AppCtx%LogInfo%EnergyEval_Stage, iErr); CHKERRQ(iErr)
   End Subroutine ComputeEnergy
   
   Subroutine SimplePoissonFinalize(AppCtx)   
      Type(AppCtx_Type)                            :: AppCtx

      PetscInt                                     :: iErr

      Call SectionRealDestroy(AppCtx%U, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(AppCtx%F, iErr); CHKERRQ(iErr)
      Call VecScatterDestroy(AppCtx%Scatter, iErr); CHKERRQ(iErr)
      Call SectionIntDestroy(AppCtx%BCFlag, iErr); CHKERRQ(iErr)
      Call MatDestroy(AppCtx%K, iErr); CHKERRQ(iErr)
      Call VecDestroy(AppCtx%RHS, iErr); CHKERRQ(iErr)
      Call KSPDestroy(AppCtx%KSP, iErr); CHKERRQ(iErr)
      Call MeshDestroy(AppCtx%MeshTopology%Mesh, iErr); CHKERRQ(ierr)
      
      If (AppCtx%AppParam%verbose) Then
         Call PetscViewerFlush(AppCtx%AppParam%MyLogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerFlush(AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerDestroy(AppCtx%AppParam%MyLogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerDestroy(AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
      End If
      Call MEF90_Finalize()
   End Subroutine SimplePoissonFinalize

#if defined PB_2D
End Module m_SimplePoisson2D
#elif defined PB_3D
End Module m_SimplePoisson3D
#endif
