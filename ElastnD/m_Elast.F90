#if defined PB_2D
Module m_Elast2D
#elif defined PB_3D
Module m_Elast3D
#endif

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscmeshdef.h"
#include "finclude/petscviewerdef.h"

   Use m_MEF90
   Use m_VarFrac_Struct
   Use petsc
   Use petscvec
   Use petscmat
   Use petscksp
   Use petscmesh

   Implicit NONE   
#include "include/finclude/tao_solver.h"
   
   Type LogInfo_Type
      PetscLogStage               :: IO_Stage
      PetscLogStage               :: Setup_Stage      
      PetscLogStage               :: MeshDistribute_Stage
      PetscLogStage               :: MatAssemblyU_Stage
      PetscLogStage               :: RHSAssemblyU_Stage
      PetscLogStage               :: SolveU_Stage
      PetscLogStage               :: PostProc_Stage
      
      PetscLogEvent               :: MatAssemblyLocalU_Event
      PetscLogEvent               :: RHSAssemblyLocalU_Event
      PetscLogEvent               :: PostProc_Event
   End Type LogInfo_Type

   Type AppParam_Type
      PetscTruth                                   :: Restart
      PetscInt                                     :: Verbose
      Character(len=MEF90_MXSTRLEN)                :: prefix
      Type(PetscViewer)                            :: LogViewer, MyLogViewer
      Type(PetscViewer)                            :: EnergyViewer
   End Type AppParam_Type

   Type AppCtx_Type
      Type (MeshTopology_Type)                     :: MeshTopology
      Type (EXO_Type)                              :: EXO, MyEXO
#if defined PB_2D
      Type(Element2D_Elast), Dimension(:), Pointer :: ElemVect
      Type(Element2D_Scal), Dimension(:), Pointer  :: ElemScal
#elif defined PB_3D
      Type(Element3D_Elast), Dimension(:), Pointer :: ElemVect
      Type(Element3D_Scal), Dimension(:), Pointer  :: ElemScal
#endif
      Type(SectionReal)                            :: U
      Type(Vec)                                    :: U_Vec
      Type(SectionReal)                            :: StressU
      Type(SectionReal)                            :: StrainU
      Type(SectionReal)                            :: F
      Type(SectionReal)                            :: Theta
      PetscInt                                     :: NumTimeSteps
      PetscInt                                     :: TimeStep
      PetscReal                                    :: Load
      PetscReal                                    :: Time
      PetscReal                                    :: ElasticEnergy
      PetscReal                                    :: ExtForcesWork
      PetscReal                                    :: TotalEnergy
      Type(VecScatter)                             :: ScatterVect
      Type(VecScatter)                             :: ScatterScal
      Type(SectionInt)                             :: BCFlagU
      Type(Mat)                                    :: KU
      Type(SectionReal)                            :: RHSU
      Type(KSP)                                    :: KSPU
      Type(PC)                                     :: PCU
!!! TAO stuff      
      TAO_SOLVER                                   :: taoU
      TAO_APPLICATION                              :: taoappU

      
      Type(LogInfo_Type)                           :: LogInfo
#if defined PB_2D
      Type(MatProp2D_Type), Dimension(:), Pointer  :: MatProp      
#elif defined PB_3D
      Type(MatProp3D_Type), Dimension(:), Pointer  :: MatProp
#endif
      Type(AppParam_Type)                          :: AppParam
      Type(VarFracSchemeParam_Type)                :: VarFracSchemeParam
   End Type AppCtx_Type
   
   
Contains

   Subroutine ElastInit(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx

      PetscInt                                     :: iErr, i
      Type(Mesh)                                   :: Tmp_Mesh
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer, filename
      PetscInt                                     :: NumComponents
      PetscTruth                                   :: HasPrefix, flag
      
      PetscReal                                    :: rDummy
      Character                                    :: cDummy
      PetscInt                                     :: vers
      PetscReal                                    :: fatol, frtol, catol, crtol, gatol, grtol, gttol
      Type(Vec)                                    :: LowerBoundU, UpperBoundU
      
      Call MEF90_Initialize()
      Call TaoInitialize(PETSC_NULL_CHARACTER, iErr); CHKERRQ(iErr)
      Call InitLog(AppCtx)

      Call PetscLogStagePush(AppCtx%LogInfo%Setup_Stage, iErr); CHKERRQ(iErr)

      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-verbose', AppCtx%AppParam%Verbose, flag, iErr)    
      Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p',     AppCtx%AppParam%prefix, HasPrefix, iErr); CHKERRQ(iErr)
      If (.NOT. HasPrefix) Then
         Call PetscPrintf(PETSC_COMM_WORLD, "No input file prefix given\n", iErr)
         Call MEF90_Finalize()
         STOP
      End If
      Call VarFracSchemeParam_GetFromOptions(AppCtx%VarFracSchemeParam)
      
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
      Write(filename, "(A,'.ener')") Trim(AppCtx%AppParam%prefix)
      Call PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, AppCtx%AppParam%EnergyViewer, iErr); CHKERRQ(iErr);   
      Write(IOBuffer, 105) Trim(filename)
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
105 Format('Energies saved in file ', A, '\n')

      Call Write_EXO_Case(AppCtx%AppParam%prefix, '%0.4d', MEF90_NumProcs)
      AppCtx%EXO%Comm = PETSC_COMM_WORLD
      AppCtx%EXO%filename = Trim(AppCtx%AppParam%prefix)//'.gen'

      !!! Reading and distributing sequential mesh
      If (MEF90_NumProcs == 1) Then
         Call MeshCreateExodus(PETSC_COMM_WORLD, AppCtx%EXO%filename, AppCtx%MeshTopology%mesh, ierr); CHKERRQ(iErr)
      Else
         Call MeshCreateExodus(PETSC_COMM_WORLD, AppCtx%EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
         Call PetscLogStagePush(AppCtx%LogInfo%MeshDistribute_Stage, iErr); CHKERRQ(iErr)
         Call MeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, AppCtx%MeshTopology%mesh, ierr); CHKERRQ(iErr)
         Call PetscLogStagePop(iErr); CHKERRQ(iErr)
         Call MeshDestroy(Tmp_mesh, ierr); CHKERRQ(iErr)
      End If
   
      Call MeshTopologyReadEXO(AppCtx%MeshTopology, AppCtx%EXO)
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) "Done reading and partitioning the mesh\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
   
      AppCtx%MyEXO%comm = PETSC_COMM_SELF
      AppCtx%MyEXO%exoid = AppCtx%EXO%exoid
      Write(AppCtx%MyEXO%filename, 99) trim(AppCtx%AppParam%prefix), MEF90_MyRank
 99  Format(A, '-', I4.4, '.gen')
   
      !!! Initializes the values and names of the properties and variables
      Call VarFracEXOVariable_Init(AppCtx%MyEXO)
      Call EXOProperty_Read(AppCtx%MyEXO)   
      If (AppCtx%AppParam%verbose > 1) Then
         Write(IOBuffer, *) "Done with VarFracEXOVariable_Init and VarFracEXOProperty_Read\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      
      !!! Read Mat Properties from the CST file
      Call MatProp_Read(AppCtx%MeshTopology, AppCtx%MatProp, trim(AppCtx%AppParam%prefix)//'.CST')
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) "Done with MatProp_Read\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If

      !!! Set the element type for each block so that we can call ElementInit
      Do i = 1, AppCtx%MeshTopology%num_elem_blks
         AppCtx%MeshTopology%elem_blk(i)%Elem_Type = AppCtx%MyEXO%EBProperty( VarFrac_EBProp_Elem_Type )%Value( AppCtx%MeshTopology%elem_blk(i)%ID )
         Call Init_Elem_Blk_Type(AppCtx%MeshTopology%Elem_Blk(i), AppCtx%MeshTopology%num_dim)
      End Do
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) "Done with Init_Elem_Blk_Type\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      
      
      Call ElementInit(AppCtx%MeshTopology, AppCtx%ElemVect, AppCtx%VarFracSchemeParam%IntegOrder)
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) "Done with ElementInit Vect\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      Call ElementInit(AppCtx%MeshTopology, AppCtx%ElemScal, AppCtx%VarFracSchemeParam%IntegOrder)
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) "Done with ElementInit Scal\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If

      If (AppCtx%AppParam%verbose > 1) Then
         Call MeshTopologyView(AppCtx%MeshTopology, AppCtx%AppParam%MyLogViewer)
         Call EXOView(AppCtx%MyEXO, AppCtx%AppParam%MyLogViewer)
      End If

      !!! Create the Sections for the variables
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'U', AppCtx%MeshTopology%Num_Dim, AppCtx%U, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'RHSU', AppCtx%MeshTopology%Num_Dim, AppCtx%RHSU, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'F', AppCtx%MeshTopology%Num_Dim, AppCtx%F, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'Theta', 1, AppCtx%Theta, iErr); CHKERRQ(iErr)

      NumComponents = AppCtx%MeshTopology%Num_Dim * (AppCtx%MeshTopology%Num_Dim + 1) / 2
      If ( (AppCtx%VarFracSchemeParam%SaveStress) .OR. ( AppCtx%VarFracSchemeParam%SaveStrain) ) Then
         Call MeshGetCellSectionReal(AppCtx%MeshTopology%mesh, 'Strain', NumComponents, AppCtx%StrainU, iErr); CHKERRQ(iErr)
         Call MeshGetCellSectionReal(AppCtx%MeshTopology%mesh, 'Stress', NumComponents, AppCtx%StressU, iErr); CHKERRQ(iErr)
      End If


      !!! Create the Scatter, Vec, Mat, KSP, PC
      Call MeshCreateGlobalScatter(AppCtx%MeshTopology%mesh, AppCtx%U, AppCtx%ScatterVect, iErr); CHKERRQ(iErr)
      Call MeshCreateGlobalScatter(AppCtx%MeshTopology%mesh, AppCtx%Theta, AppCtx%ScatterScal, iErr); CHKERRQ(iErr)
      Call MeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%U, AppCtx%U_Vec, iErr); CHKERRQ(iErr)

      If (AppCtx%VarFracSchemeParam%U_UseTao) Then
         Call MeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%U, LowerBoundU, iErr); CHKERRQ(iErr)
         Call MeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%U, UpperBoundU, iErr); CHKERRQ(iErr)
      End If
   
      Call MeshSetMaxDof(AppCtx%MeshTopology%Mesh, AppCtx%MeshTopology%Num_Dim, iErr); CHKERRQ(iErr) 
      Call MeshCreateMatrix(AppCtx%MeshTopology%mesh, AppCtx%U, MATMPIAIJ, AppCtx%KU, iErr); CHKERRQ(iErr)
      
      If (AppCtx%VarFracSchemeParam%U_UseTao) Then
         Call TaoCreate(PETSC_COMM_WORLD, 'tao_tron', AppCtx%taoU, iErr); CHKERRQ(iErr)
         Call TaoApplicationCreate(PETSC_COMM_WORLD, AppCtx%taoappU, iErr); CHKERRQ(iErr)

         Call TaoAppSetObjectiveAndGradientRoutine(AppCtx%taoappU, FormFunctionAndGradient, AppCtx, iErr); CHKERRQ(iErr)
         Call TaoAppSetHessianRoutine(AppCtx%taoappU, FormHessian, AppCtx, iErr); CHKERRQ(iErr)
         Call TaoAppSetVariableBoundsRoutine(AppCtx%taoappU, FormBoundsTao, AppCtx, iErr); CHKERRQ(iErr)

         Call TaoAppSetHessianMat(AppCtx%taoappU, AppCtx%KU, AppCtx%KU, iErr); CHKERRQ(iErr)
         Call TaoAppSetInitialSolutionVec(AppCtx%taoappU, AppCtx%U_Vec, iErr); CHKERRQ(iErr)


         Call TaoSetOptions(AppCtx%taoappU, AppCtx%taoU, iErr); CHKERRQ(iErr)
         Call TaoAppSetFromOptions(AppCtx%taoappU, iErr); CHKERRQ(iErr)
         Call TaoAppGetKSP(AppCtx%taoappU, AppCtx%KSPU, iErr); CHKERRQ(iErr)
      Else
         Call KSPCreate(PETSC_COMM_WORLD, AppCtx%KSPU, iErr); CHKERRQ(iErr)
         Call KSPSetOperators(AppCtx%KSPU, AppCtx%KU, AppCtx%KU, SAME_NONZERO_PATTERN, iErr); CHKERRQ(iErr)
         Call KSPSetInitialGuessNonzero(AppCtx%KSPU, PETSC_TRUE, iErr); CHKERRQ(iErr)
      End If

      Call KSPAppendOptionsPrefix(AppCtx%KSPU, "U_", iErr); CHKERRQ(iErr)
      Call KSPSetType(AppCtx%KSPU, KSPCG, iErr); CHKERRQ(iErr)
      Call KSPSetFromOptions(AppCtx%KSPU, iErr); CHKERRQ(iErr)

      Call KSPGetPC(AppCtx%KSPU, AppCtx%PCU, iErr); CHKERRQ(iErr)
      Call PCSetType(AppCtx%PCU, PCBJACOBI, iErr); CHKERRQ(iErr)
      Call PCSetFromOptions(AppCtx%PCU, iErr); CHKERRQ(iErr)

      !!! Create the Section for the BC
      Call MeshGetVertexSectionInt(AppCtx%MeshTopology%mesh, 'BCU', AppCtx%MeshTopology%Num_Dim, AppCtx%BCFlagU, iErr); CHKERRQ(iErr)
#if defined PB_2D
      Call EXOProperty_InitBCUFlag2D(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%BCFlagU)
#elif defined PB_3D
      Call EXOProperty_InitBCUFlag3D(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%BCFlagU)
#endif

      !!! Get the number of time steps
      AppCtx%MyEXO%exoid = EXOPEN(AppCtx%MyEXO%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers, iErr)
      Call EXINQ(AppCtx%MyEXO%exoid, EXTIMS, AppCtx%NumTimeSteps, rDummy, cDummy, iErr)
      Call EXCLOS(AppCtx%MyEXO%exoid, iErr)
      AppCtx%MyEXO%exoid = 0

      AppCtx%TimeStep = 1
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine ElastInit


   Subroutine InitLog(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr
      
      Call PetscLogEventRegister('MatAssembly Local U', 0, AppCtx%LogInfo%MatAssemblyLocalU_Event, ierr); CHKERRQ(ierr)
      Call PetscLogEventRegister('RHSAssembly Local U', 0, AppCtx%LogInfo%RHSAssemblyLocalU_Event, ierr); CHKERRQ(ierr)
      Call PetscLogEventRegister('Post Processing',     0, AppCtx%LogInfo%PostProc_Event,          ierr); CHKERRQ(ierr)

      Call PetscLogStageRegister("Setup",          AppCtx%LogInfo%Setup_Stage,          iErr)
      Call PetscLogStageRegister("MeshDistribute", AppCtx%LogInfo%MeshDistribute_Stage, iErr)
      Call PetscLogStageRegister("Mat Assembly U", AppCtx%LogInfo%MatAssemblyU_Stage,   iErr)
      Call PetscLogStageRegister("RHS Assembly U", AppCtx%LogInfo%RHSAssemblyU_Stage,   iErr)
      Call PetscLogStageRegister("Solve U",        AppCtx%LogInfo%SolveU_Stage,         iErr)
      Call PetscLogStageRegister("IO Stage",       AppCtx%LogInfo%IO_Stage,             iErr)
      Call PetscLogStageRegister("Post Proc",      AppCtx%LogInfo%PostProc_Stage,       iErr)
   End Subroutine InitLog


!----------------------------------------------------------------------------------------!      
! Solve (CM)   
!----------------------------------------------------------------------------------------!      
   Subroutine FormBoundsTao(TaoApp, LowerBound, UpperBound, AppCtx, iErr)
      TAO_APPLICATION                              :: TaoApp
      Type(Vec)                                    :: LowerBound, UpperBound
      Type(AppCtx_Type)                            :: AppCtx
      PetscErrorCode                               :: iErr

      Type(SectionReal)                            :: UBC_Sec, LowerBoundU_Sec, UpperBoundU_Sec
      PetscInt, Dimension(:), Pointer              :: BCFlag_Ptr
      PetscReal, Dimension(:), Pointer             :: UBC_Ptr, LowerBoundU_Ptr, UpperBoundU_Ptr
      PetscInt                                     :: iDoF1, iDoF2
      Character(len=256)                           :: IOBuffer
      
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) "Updating TAO bounds\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If

      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'UBC',  AppCtx%MeshTopology%Num_Dim, UBC_Sec , iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'LowerBoundU_Sec',  AppCtx%MeshTopology%Num_Dim, LowerBoundU_Sec , iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'LowerBoundU_Sec',  AppCtx%MeshTopology%Num_Dim, UpperBoundU_Sec , iErr); CHKERRQ(iErr)
      Allocate(LowerBoundU_Ptr(AppCtx%MeshTopology%Num_Dim))
      Allocate(UpperBoundU_Ptr(AppCtx%MeshTopology%Num_Dim))

      Call SectionRealSet(LowerBoundU_Sec, -1.0D20, iErr); CHKERRQ(iErr)
      Call SectionRealSet(UpperBoundU_Sec,  1.0D20, iErr); CHKERRQ(iErr)
         
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_DisplacementX)%Offset, AppCtx%TimeStep, UBC_Sec) 

      Do iDoF1 = 1, AppCtx%MeshTopology%Num_Verts !!! WRONG!
         Call SectionIntRestrict(AppCtx%BCFlagU, AppCtx%MeshTopology%Num_Elems+iDoF1-1, BCFlag_Ptr, iErr); CHKERRQ(iErr)
         If (Sum(BCFlag_Ptr) /= 0) Then
            LowerBoundU_Ptr = -1.0D20
            UpperBoundU_Ptr =  1.0D20
            Call SectionRealRestrict(UBC_Sec, AppCtx%MeshTopology%Num_Elems+iDoF1-1, UBC_Ptr, iErr); CHKERRQ(iErr)
            Do iDof2 = 1, AppCtx%MeshTopology%Num_Dim
               If (BCFlag_Ptr(iDoF2) /= 0) Then
                  LowerBoundU_Ptr(iDoF2) = UBC_Ptr(iDoF2)
                  UpperBoundU_Ptr(iDoF2) = UBC_Ptr(iDoF2)
               End If
            End Do
            Call SectionRealRestore(UBC_Sec, AppCtx%MeshTopology%Num_Elems+iDoF1-1, UBC_Ptr, iErr); CHKERRQ(iErr)
            Call SectionRealUpdate(LowerBoundU_Sec, AppCtx%MeshTopology%Num_Elems+iDoF1-1, LowerBoundU_Ptr, INSERT_VALUES, iErr); CHKERRQ(iErr)
            Call SectionRealUpdate(UpperBoundU_Sec, AppCtx%MeshTopology%Num_Elems+iDoF1-1, UpperBoundU_Ptr, INSERT_VALUES, iErr); CHKERRQ(iErr)
         EndIf
         Call SectionIntRestore(AppCtx%BCFlagU, AppCtx%MeshTopology%Num_Elems+iDoF1-1, BCFlag_Ptr, iErr); CHKERRQ(iErr)
      End Do
      Call SectionRealDestroy(UBC_Sec, iErr); CHKERRQ(iErr)

      Call SectionRealToVec(LowerBoundU_Sec, AppCtx%ScatterVect, SCATTER_FORWARD, LowerBound, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(UpperBoundU_Sec, AppCtx%ScatterVect, SCATTER_FORWARD, UpperBound, iErr); CHKERRQ(iErr)

      Call SectionRealDestroy(LowerBoundU_Sec, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(UpperBoundU_Sec, iErr); CHKERRQ(iErr)
      DeAllocate(LowerBoundU_Ptr)
      DeAllocate(UpperBoundU_Ptr)
   End Subroutine FormBoundsTao   
   
   Subroutine InitU(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscErrorCode                               :: iErr

      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_DisplacementX)%Offset, AppCtx%TimeStep, AppCtx%U) 
   End Subroutine InitU
   
   Subroutine InitLoads(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscErrorCode                               :: iErr

      !!! Read F, and Temperature
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_ForceX)%Offset, AppCtx%TimeStep, AppCtx%F) 
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Temperature)%Offset, AppCtx%TimeStep, AppCtx%Theta) 
   End Subroutine InitLoads
   
   Subroutine Solve(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      KSPConvergedReason                           :: KSPreason
      TaoTerminateReason                           :: TaoReason
      PetscReal                                    :: TaoResidual
      PetscInt                                     :: KSPNumIter
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
      Integer                                      :: iDum
      PetscTruth                                   :: flg
      Type(Vec)                                    :: RHSU_Vec
      Type(SectionReal)                            :: UBC_Sec, LowerBoundU_Sec, UpperBoundU_Sec
      PetscInt, Dimension(:), Pointer              :: BCFlag_Ptr
      PetscReal, Dimension(:), Pointer             :: UBC_Ptr, U_Ptr, LowerBoundU_Ptr, UpperBoundU_Ptr
      PetscInt                                     :: iDoF1, iDoF2
      
      PetscReal      fatol, frtol, gatol, grtol, gttol, radius, catol, crtol
   
      Call InitLoads(AppCtx)
      Call InitU(AppCtx)      

      !!! KSPSolve and Tao work with vectors, not sections
      Call SectionRealToVec(AppCtx%U, AppCtx%ScatterVect, SCATTER_FORWARD, AppCtx%U_Vec, ierr); CHKERRQ(ierr)

      If (AppCtx%VarFracSchemeParam%U_UseTao) Then
         If (AppCtx%AppParam%verbose > 0) Then
            Write(IOBuffer, *) 'Calling TaoSolveApplication\n'
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         End If

         Call TaoSolveApplication(AppCtx%taoappU, AppCtx%taoU, iErr); CHKERRQ(iErr)
         Call TaoGetSolutionStatus(AppCtx%taoU, KSPNumIter, AppCtx%TotalEnergy, TaoResidual, iDum, iDum, TaoReason, iErr); CHKERRQ(iErr)
         If ( TaoReason > 0) Then
            Write(IOBuffer, 102) KSPNumiter, TAOReason
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Else
            Write(IOBuffer, 103) TaoReason
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         End If
      Else
         If (AppCtx%AppParam%verbose > 0) Then
            Write(IOBuffer, *) 'Assembling the RHS\n'
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         End If
         Call RHSAssembly(AppCtx)
         If (AppCtx%AppParam%verbose > 1) Then
            Call SectionRealView(AppCtx%RHSU, AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
         End If
         Call MeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%RHSU, RHSU_Vec, iErr); CHKERRQ(iErr)
         Call SectionRealToVec(AppCtx%RHSU, AppCtx%ScatterVect, SCATTER_FORWARD, RHSU_Vec, iErr); CHKERRQ(ierr)
      
         If (AppCtx%AppParam%verbose > 0) Then
            Write(IOBuffer, *) 'Calling KSPSolve\n'
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         End If
   
         Call PetscLogStagePush(AppCtx%LogInfo%SolveU_Stage, iErr); CHKERRQ(iErr)
         Call KSPSolve(AppCtx%KSPU, RHSU_Vec, AppCtx%U_Vec, iErr); CHKERRQ(iErr)
         Call PetscLogStagePop(iErr); CHKERRQ(iErr)
         
         Call KSPGetConvergedReason(AppCtx%KSPU, KSPreason, iErr); CHKERRQ(iErr)
         If ( KSPreason > 0) Then
            Call KSPGetIterationNumber(AppCtx%KSPU, KSPNumIter, iErr); CHKERRQ(iErr)
            Write(IOBuffer, 100) KSPNumIter
         Else
            Write(IOBuffer, 101) KSPreason
         End If
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Call VecDestroy(RHSU_Vec, iErr); CHKERRQ(iErr)
      End If
      Call SectionRealToVec(AppCtx%U, AppCtx%ScatterVect, SCATTER_REVERSE, AppCtx%U_Vec, iErr); CHKERRQ(ierr)
      
100 Format('     KSP for U converged in ', I5, ' iterations \n')
101 Format('[ERROR] KSP for U diverged. KSPConvergedReason is ', I2, '\n')      
102 Format('     TAO solver converged in ', I5, ' iterations. Tao termination reason is ', I2, '\n')
103 Format('[ERROR] TaoSolveApplication did not converge. ', I2, '\n')      
   End Subroutine Solve
   
!----------------------------------------------------------------------------------------!      
!             Global assembly routines
!----------------------------------------------------------------------------------------!      
   Subroutine MatAssembly(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iBlk, iErr
      PetscReal, Dimension(:,:), Pointer           :: MatElem
      
      Call PetscLogStagePush(AppCtx%LogInfo%MatAssemblyU_Stage, iErr); CHKERRQ(iErr)

      Call MatZeroEntries(AppCtx%KU, iErr); CHKERRQ(iErr)
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Call HessianAssemblyBlock(iBlk, AppCtx%KU, .TRUE.,  AppCtx)
      End Do Do_Elem_iBlk

      Call MatAssemblyBegin(AppCtx%KU, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (AppCtx%KU, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)

      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine MatAssembly
   
   Subroutine FormHessian(tao, X, H, Hpre, flg, AppCtx, iErr)
      TAO_SOLVER         :: tao
      Type(Vec)          :: X
      Type(Mat)          :: H, Hpre
      PetscInt           :: iErr
      MatStructure       :: flg
      Type(AppCtx_Type)  :: AppCtx
      
      PetscInt           :: iBlk
      
      Call PetscLogStagePush(AppCtx%LogInfo%MatAssemblyU_Stage, iErr); CHKERRQ(iErr)

      Call MatZeroEntries(H, iErr); CHKERRQ(iErr)

      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Call HessianAssemblyBlock(iBlk, H, .FALSE.,  AppCtx)
      End Do Do_Elem_iBlk
      Call MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (H, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      If (AppCtx%AppParam%verbose > 1) Then
         Call MatView(H, PETSC_VIEWER_STDOUT_WORLD, iErr)
      End If

      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine FormHessian
   
   
   Subroutine FormFunctionAndGradient(tao, X, func, Gradient, AppCtx, iErr)
      TAO_SOLVER                                   :: tao
      Type(Vec)                                    :: X, Gradient
      PetscInt                                     :: iErr
      PetscReal                                    :: func
      Type(AppCtx_Type)                            :: AppCtx
      
      Type(SectionReal)                            :: Gradient_Sec, X_Sec
      PetscInt                                     :: iBlk     
      PetscReal                                    :: myfunc, ElasticEnergyBlock, ExtForcesWorkBlock 
      
      
      Call PetscLogStagePush(AppCtx%LogInfo%RHSAssemblyU_Stage, iErr); CHKERRQ(iErr)

      func = 0.0_Kr
            
      Call VecZeroEntries(Gradient, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'Gradient', AppCtx%MeshTopology%Num_Dim, Gradient_Sec, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'X', AppCtx%MeshTopology%Num_Dim, X_Sec, iErr); CHKERRQ(iErr)

      Call SectionRealToVec(X_Sec, AppCtx%ScatterVect, SCATTER_REVERSE, X, iErr); CHKERRQ(ierr)

      myfunc = 0.0_Kr
      Do_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Call ComputeEnergiesBlock(iBlk, X_Sec, ElasticEnergyBlock, ExtForcesWorkBlock, AppCtx)
         myfunc = myfunc + ElasticEnergyBlock - ExtForcesWorkBlock
         Call FormGradientBlock(iBlk, X_Sec, Gradient_Sec, AppCtx)
      End Do Do_iBlk
      Call PetscGlobalSum(myfunc, func, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)

      Call SectionRealComplete(Gradient_Sec, iErr); CHKERRQ(iErr)
      !!! Scatter values from the Sections back to the Vec
      Call SectionRealToVec(Gradient_Sec, AppCtx%ScatterVect, SCATTER_FORWARD, Gradient, iErr); CHKERRQ(iErr)

      Call SectionRealDestroy(Gradient_Sec, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(X_Sec, iErr); CHKERRQ(iErr)
      
      If (AppCtx%AppParam%verbose > 1) Then
         Call VecView(Gradient, PETSC_VIEWER_STDOUT_WORLD, iErr)
      End If
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine FormFunctionAndGradient
   
   Subroutine ComputeEnergies(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iBlk, iErr     
      PetscReal                                    :: ElasticEnergyBlock, ExtForcesWorkBlock 
      PetscReal                                    :: MyElasticEnergy, MyExtForcesWork 
      
      
      Call PetscLogStagePush(AppCtx%LogInfo%RHSAssemblyU_Stage, iErr); CHKERRQ(iErr)

      MyElasticEnergy = 0.0_Kr
      MyExtForcesWork = 0.0_Kr
            
      Do_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Call ComputeEnergiesBlock(iBlk, AppCtx%U, ElasticEnergyBlock, ExtForcesWorkBlock, AppCtx)
         MyElasticEnergy = MyElasticEnergy + ElasticEnergyBlock
         MyExtForcesWork = MyExtForcesWork + ExtForcesWorkBlock
      End Do Do_iBlk
      Call PetscGlobalSum(MyElasticEnergy, AppCtx%ElasticEnergy, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      Call PetscGlobalSum(MyExtForcesWork, AppCtx%ExtForcesWork, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)

      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine ComputeEnergies

   Subroutine RHSAssembly(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      PetscInt                                     :: iBlk

      Call PetscLogStagePush(AppCtx%LogInfo%RHSAssemblyU_Stage, iErr); CHKERRQ(iErr)

      Call SectionRealZero(AppCtx%RHSU, iErr); CHKERRQ(iErr)
      
      Do_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Call RHSAssemblyBlock(iBlk, AppCtx)
      End Do Do_iBlk

      Call SectionRealComplete(AppCtx%RHSU, iErr); CHKERRQ(iErr)
      !!! VERY important! This is the equivalent of a ghost update
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine RHSAssembly

!----------------------------------------------------------------------------------------!      
!             Block assembly routines
!----------------------------------------------------------------------------------------!      
   Subroutine HessianAssemblyBlock(iBlk, H, DoBC, AppCtx)
      PetscInt                                     :: iBlk
      Type(Mat)                                    :: H
      Logical                                      :: DoBC
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iBlkID
      PetscReal, Dimension(:,:), Pointer           :: MatElem      
      PetscInt                                     :: iELoc, iE
   
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoF
      PetscInt, Dimension(:), Pointer              :: BCFlag
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      PetscLogDouble                               :: flops       
      
      Call PetscLogEventBegin(AppCtx%LogInfo%MatAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
      flops = 0.0
      NumDoF   =  AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
      Allocate(MatElem(NumDoF, NumDoF))
      Allocate(BCFlag(NumDoF))
      BCFlag = 0

      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID

      Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         If (DoBC) Then
            Call SectionIntRestrictClosure(AppCtx%BCFlagU, AppCtx%MeshTopology%mesh, iE-1, NumDoF, BCFlag, iErr); CHKERRQ(ierr)      
         End If
         MatElem  = 0.0_Kr
         Do iGauss = 1, Size(AppCtx%ElemVect(iE)%Gauss_C)
            Do iDoF1 = 1, NumDoF
               If (BCFlag(iDOF1) == 0) Then
                  Do iDoF2 = 1, NumDoF
                     MatElem(iDoF2, iDoF1) =  MatElem(iDoF2, iDoF1) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * ( (AppCtx%MatProp(iBlkID)%Hookes_Law * AppCtx%ElemVect(iE)%GradS_BF(iDoF1, iGauss) ) .DotP. AppCtx%ElemVect(iE)%GradS_BF(iDoF2, iGauss))        
                     flops = flops + 2.0     
                  End Do
               End If
            End Do
         End Do
         Call assembleMatrix(H, AppCtx%MeshTopology%mesh, AppCtx%U, iE-1, MatElem, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_Elem_iE
            
      DeAllocate(MatElem)
      DeAllocate(BCFlag)
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
      Call PetscLogEventEnd(AppCtx%LogInfo%MatAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
   End Subroutine HessianAssemblyBlock

   Subroutine FormGradientBlock(iBlk, X_Sec, Gradient_Sec, AppCtx)
      PetscInt                                     :: iBlk
      Type(SectionReal)                            :: X_Sec
      Type(SectionReal)                            :: Gradient_Sec
      Type(AppCtx_Type)                            :: AppCtx
      
      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_ELem = \sum_i U_Loc(i) BF(i))
      PetscReal, Dimension(:), Pointer             :: X_Loc, F_Loc, Theta_Loc, Gradient_Loc
#if defined PB_2D
      Type(Vect2D)                                 :: X_Elem, F_Elem
      Type(Mats2D)                                 :: Strain_Elem, EffectiveStrain_Elem
#elif defined PB_3D  
      Type(Vect3D)                                 :: X_ElemF_Elem    
      Type(Mats3D)                                 :: Strain_Elem, EffectiveStrain_Elem
#endif
      PetscReal                                    :: Theta_Elem
      PetscInt                                     :: iE, iEloc, iBlkId, iErr
      PetscInt                                     :: NumDoFScal, NumDoFVect
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      PetscLogDouble                               :: flops       

      Call PetscLogEventBegin(AppCtx%LogInfo%RHSAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
      flops        = 0.0

      NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
      NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF

      Allocate(X_Loc(NumDoFVect))
      Allocate(F_Loc(NumDoFVect))
      Allocate(Gradient_Loc(NumDoFVect))
      Allocate(Theta_Loc(NumDoFScal))

      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID

      Do_iEloc: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         Gradient_Loc = 0.0_Kr
         Call SectionRealRestrictClosure(X_Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, X_Loc, iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(AppCtx%F, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, F_Loc, iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(AppCtx%Theta, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, Theta_Loc, iErr); CHKERRQ(ierr)
         Do_iGauss: Do iGauss = 1, size(AppCtx%ElemVect(iE)%Gauss_C)
            F_Elem = 0.0_Kr
            X_Elem = 0.0_Kr
            Strain_Elem = 0.0_Kr
            Do iDoF2 = 1, NumDoFVect
               F_Elem      = F_Elem + AppCtx%ElemVect(iE)%BF(iDoF2, iGauss) * F_Loc(iDoF2)
               X_Elem      = X_Elem + X_Loc(iDoF2) * AppCtx%ElemVect(iE)%BF(iDoF2, iGauss)
               Strain_Elem = Strain_Elem + X_Loc(iDoF2) * AppCtx%ElemVect(iE)%GradS_BF(iDoF2, iGauss)
            End Do
            Theta_Elem = 0.0_Kr
            Do iDoF2 = 1, NumDoFScal
               Theta_Elem = Theta_Elem + AppCtx%ElemScal(iE)%BF(iDoF2, iGauss) * Theta_Loc(iDoF2)
               flops = flops + 2.0
            End Do
            EffectiveStrain_Elem = Strain_Elem - AppCtx%MatProp(iBlkId)%Therm_Exp * Theta_Elem
            
            Do iDoF1 = 1, NumDofVect
               Gradient_Loc(iDoF1) = Gradient_Loc(iDoF1) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * ( ((AppCtx%MatProp(iBlkId)%Hookes_Law * EffectiveStrain_Elem) .DotP. AppCtx%ElemVect(iE)%GradS_BF(iDoF1, iGauss)) - (F_Elem .DotP. AppCtx%ElemVect(iE)%BF(iDoF1, iGauss)) )
            End Do
         End Do Do_iGauss
         Call SectionRealUpdateClosure(Gradient_Sec, AppCtx%MeshTopology%Mesh, iE-1, Gradient_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iEloc

      DeAllocate(X_Loc)
      DeAllocate(F_Loc)
      DeAllocate(Gradient_Loc)
      DeAllocate(Theta_Loc)
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
      Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
   End Subroutine FormGradientBlock

   Subroutine ComputeEnergiesBlock(iBlk, X_Sec, ElasticEnergy, ExtForcesWork, AppCtx)
      PetscInt                                     :: iBlk
      Type(SectionReal)                            :: X_Sec
      PetscReal                                    :: ElasticEnergy, ExtForcesWork
      Type(AppCtx_Type)                            :: AppCtx
      
      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_ELem = \sum_i U_Loc(i) BF(i))
      PetscReal, Dimension(:), Pointer             :: X_Loc, F_Loc, Theta_Loc, Gradient_Loc
#if defined PB_2D
      Type(Vect2D)                                 :: X_Elem, F_Elem
      Type(Mats2D)                                 :: Strain_Elem, EffectiveStrain_Elem
#elif defined PB_3D  
      Type(Vect3D)                                 :: X_ElemF_Elem    
      Type(Mats3D)                                 :: Strain_Elem, EffectiveStrain_Elem
#endif
      PetscReal                                    :: Theta_Elem
      PetscInt                                     :: iE, iEloc, iBlkId, iErr
      PetscInt                                     :: NumDoFScal, NumDoFVect
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      PetscLogDouble                               :: flops       

      Call PetscLogEventBegin(AppCtx%LogInfo%RHSAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
      flops  = 0.0
      ElasticEnergy = 0.0_Kr
      ExtFOrcesWork = 0.0_Kr

      NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
      NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF

      Allocate(X_Loc(NumDoFVect))
      Allocate(F_Loc(NumDoFVect))
      Allocate(Theta_Loc(NumDoFScal))

      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID

      Do_iEloc: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         Call SectionRealRestrictClosure(X_Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, X_Loc, iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(AppCtx%F, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, F_Loc, iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(AppCtx%Theta, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, Theta_Loc, iErr); CHKERRQ(ierr)
         Do_iGauss: Do iGauss = 1, size(AppCtx%ElemVect(iE)%Gauss_C)
            F_Elem = 0.0_Kr
            X_ELem = 0.0_Kr
            Strain_Elem = 0.0_Kr
            Do iDoF2 = 1, NumDoFVect
               F_Elem      = F_Elem + AppCtx%ElemVect(iE)%BF(iDoF2, iGauss) * F_Loc(iDoF2)
               X_Elem      = X_Elem + X_Loc(iDoF2) * AppCtx%ElemVect(iE)%BF(iDoF2, iGauss)
               Strain_Elem = Strain_Elem + X_Loc(iDoF2) * AppCtx%ElemVect(iE)%GradS_BF(iDoF2, iGauss)
            End Do
            Theta_Elem = 0.0_Kr
            Do iDoF2 = 1, NumDoFScal
               Theta_Elem = Theta_Elem + AppCtx%ElemScal(iE)%BF(iDoF2, iGauss) * Theta_Loc(iDoF2)
               flops = flops + 2.0
            End Do
            EffectiveStrain_Elem = Strain_Elem - AppCtx%MatProp(iBlkId)%Therm_Exp * Theta_Elem
            
            ElasticEnergy = ElasticEnergy + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * ((AppCtx%MatProp(iBlkId)%Hookes_Law * EffectiveStrain_Elem) .DotP. EffectiveStrain_Elem ) * .5_Kr 
            ExtForcesWork = ExtForcesWork + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * (F_Elem .DotP. X_Elem)
         End Do Do_iGauss
      End Do Do_iEloc

      DeAllocate(X_Loc)
      DeAllocate(F_Loc)
      DeAllocate(Theta_Loc)
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
      Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
   End Subroutine ComputeEnergiesBlock

   Subroutine RHSAssemblyBlock(iBlk, AppCtx)
      PetscInt                                     :: iBlk
      Type(AppCtx_Type)                            :: AppCtx

      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_ELem = \sum_i U_Loc(i) BF(i))
      PetscReal, Dimension(:), Pointer             :: F_Loc, Theta_Loc, RHS_Loc
      PetscInt, Dimension(:), Pointer              :: BCFlag_Loc
      PetscReal                                    :: Theta_Elem
#if defined PB_2D
      Type (Vect2D)                                :: F_Elem
#elif defined PB_3D  
      Type (Vect3D)                                :: F_Elem    
#endif
      PetscInt                                     :: iE, iEloc, iBlkId, iErr
      PetscInt                                     :: NumDoFScal, NumDoFVect
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      PetscLogDouble                               :: flops       

      Call PetscLogEventBegin(AppCtx%LogInfo%RHSAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
      flops      = 0.0
      
      NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
      NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF

      Allocate(BCFlag_Loc(NumDoFVect))
      Allocate(F_Loc(NumDoFVect))
      Allocate(RHS_Loc(NumDoFVect))
      Allocate(Theta_Loc(NumDoFScal))

      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
      Do_iEloc: Do iEloc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         RHS_Loc = 0.0_Kr
         Call SectionIntRestrictClosure(AppCtx%BCFlagU, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, BCFlag_Loc, iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(AppCtx%F, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, F_Loc, iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(AppCtx%Theta, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, Theta_Loc, iErr); CHKERRQ(ierr)
         Do_iGauss: Do iGauss = 1, size(AppCtx%ElemVect(iE)%Gauss_C)
            F_Elem = 0.0_Kr
            Do iDoF2 = 1, NumDoFVect
               F_Elem = F_Elem + AppCtx%ElemVect(iE)%BF(iDoF2, iGauss) * F_Loc(iDoF2)
            End Do
            Theta_Elem = 0.0_Kr
            Do iDoF2 = 1, NumDoFScal
               Theta_Elem = Theta_Elem + AppCtx%ElemScal(iE)%BF(iDoF2, iGauss) * Theta_Loc(iDoF2)
               flops = flops + 2.0
            End Do
            Do iDoF1 = 1, NumDoFVect
               If (BCFlag_Loc(iDoF1) == 0) Then
                  ! RHS terms due to forces
                  RHS_Loc(iDoF1) = RHS_Loc(iDoF1) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * ( AppCtx%ElemVect(iE)%BF(iDoF1, iGauss) .DotP. F_Elem ) 
                  ! RHS terms due to inelastic strains
                  RHS_Loc(iDoF1) = RHS_Loc(iDoF1) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * Theta_Elem * ((AppCtx%MatProp(iBlkId)%Hookes_Law * AppCtx%ElemVect(iE)%GradS_BF(iDoF1, iGauss)) .DotP. AppCtx%MatProp(iBlkId)%Therm_Exp)
                  flops = flops + 5.0
               End If
            End Do
         End Do Do_iGauss
         Call SectionRealUpdateClosure(AppCtx%RHSU, AppCtx%MeshTopology%Mesh, iE-1, RHS_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iEloc
      
      DeAllocate(BCFlag_Loc)
      DeAllocate(F_Loc)
      DeAllocate(RHS_Loc)
      DeAllocate(Theta_Loc)
      
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
      Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
   End Subroutine RHSAssemblyBLock
   
!----------------------------------------------------------------------------------------!      
! ComputeStrainStress (CM)  
!----------------------------------------------------------------------------------------!      
   Subroutine ComputeStrainStress(AppCtx)
     Type(AppCtx_Type)                             :: AppCtx
      
      PetscInt                                     :: iErr
#if defined PB_2D
     Type(MatS2D)                                  :: Strain_Elem, Stress_Elem 
     Type(MatS2D)                                  :: Effective_Strain_Elem
     Type(Vect2D)                                  :: F_Elem, U_Elem  
#elif defined PB_3D 
     Type(MatS3D)                                  :: Strain_Elem, Stress_Elem 
     Type(MatS3D)                                  :: Effective_Strain_Elem
     Type(Vect3D)                                  :: F_Elem, U_Elem  
#endif
      PetscReal                                    :: Theta_Elem
      PetscReal                                    :: Vol
      PetscInt                                     :: NumDoFScal, NumDoFVect, NumGauss
      PetscReal, Dimension(:), Pointer             :: U, Theta
      PetscInt                                     :: iBlk, iELoc, iE
      PetscInt                                     :: iDoF, iGauss
      PetscReal, Dimension(:), Pointer             :: Stress_Ptr, Strain_Ptr
      PetscLogDouble                               :: flops       
       
        
      Call PetscLogStagePush (AppCtx%LogInfo%PostProc_Stage, iErr); CHKERRQ(iErr)
      Call PetscLogEventBegin(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      flops = 0.0
      Allocate(Stress_Ptr( AppCtx%MeshTopology%Num_Dim * ( AppCtx%MeshTopology%Num_Dim+1 ) / 2))
      Allocate(Strain_Ptr( AppCtx%MeshTopology%Num_Dim * ( AppCtx%MeshTopology%Num_Dim+1 ) / 2))
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            NumDoFVect = Size(AppCtx%ElemVect(iE)%BF,1)
            NumDoFScal = Size(AppCtx%ElemScal(iE)%BF,1)
            NumGauss   = Size(AppCtx%ElemVect(iE)%BF,2)
            Allocate(U(NumDoFVect))
            Call SectionRealRestrictClosure(AppCtx%U, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U, iErr); CHKERRQ(ierr)
            Allocate(Theta(NumDoFScal))
            Call SectionRealRestrictClosure(AppCtx%Theta, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, Theta, iErr); CHKERRQ(ierr)

            Strain_Elem           = 0.0_Kr
            Stress_Elem           = 0.0_Kr
            Theta_Elem            = 0.0_Kr
            Effective_Strain_Elem = 0.0_Kr
            Vol  = 0.0_Kr
            Do iGauss = 1, NumGauss
               Do iDoF = 1, NumDoFVect
                  Strain_Elem = Strain_Elem + AppCtx%ElemScal(iE)%Gauss_C(iGauss) * AppCtx%ElemVect(iE)%GradS_BF(iDoF, iGauss) * U(iDoF)
               End Do
               Do iDoF = 1, NumDoFScal
                  Theta_Elem  = Theta_Elem +  AppCtx%ElemScal(iE)%Gauss_C(iGauss) * AppCtx%ElemScal(iE)%BF(iDoF, iGauss) * Theta(iDoF)
                  Vol = Vol + AppCtx%ElemScal(iE)%Gauss_C(iGauss) * AppCtx%ElemScal(iE)%BF(iDoF, iGauss)
                  flops = flops + 5.0
               End Do
            End Do
            Strain_Elem = Strain_Elem / Vol
            Theta_Elem  = Theta_Elem / Vol
            Effective_Strain_Elem  = Strain_Elem - Theta_Elem * (AppCtx%MatProp(iBlk)%Therm_Exp) 
            Stress_Elem            = AppCtx%MatProp( AppCtx%MeshTopology%Elem_Blk(iBlk)%ID )%Hookes_Law * Effective_Strain_Elem
            flops = flops + 2.0
            
#if defined PB_2D
            Stress_Ptr = (/ Stress_Elem%XX, Stress_Elem%YY, Stress_Elem%XY /)
            Strain_Ptr = (/ Strain_Elem%XX, Strain_Elem%YY, Strain_Elem%XY /)
#elif defined PB_3D
            Stress_Ptr = (/ Stress_Elem%XX, Stress_Elem%YY, Stress_Elem%ZZ, Stress_Elem%YZ, Stress_Elem%XZ, Stress_Elem%XY  /)
            Strain_Ptr = (/ Strain_Elem%XX, Strain_Elem%YY, Strain_Elem%ZZ, Strain_Elem%YZ, Strain_Elem%XZ, Strain_Elem%XY  /)
#endif
            ! Update the Sections with the local values
            Call SectionRealUpdateClosure(AppCtx%StressU, AppCtx%MeshTopology%Mesh, iE-1, Stress_Ptr, INSERT_VALUES, iErr)
            Call SectionRealUpdateClosure(AppCtx%StrainU, AppCtx%MeshTopology%Mesh, iE-1, Strain_Ptr, INSERT_VALUES, iErr)
            DeAllocate(U)
            DeAllocate(Theta)
         End Do Do_Elem_iE
      End Do Do_Elem_iBlk
      DeAllocate(Stress_Ptr)
      DeAllocate(Strain_Ptr)
      Call SectionRealComplete(AppCtx%StressU, iErr); CHKERRQ(iErr)
      Call SectionRealComplete(AppCtx%StrainU, iErr); CHKERRQ(iErr)

      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
      Call PetscLogEventEnd(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine ComputeStrainStress

   
   Subroutine Save(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr

      Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
      Call Write_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_DisplacementX)%Offset, AppCtx%TimeStep, AppCtx%U) 
      If (AppCtx%VarFracSchemeParam%SaveStress) Then
         Call Write_EXO_Result_Cell(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%CellVariable(VarFrac_CellVar_StressXX)%Offset, AppCtx%TimeStep, AppCtx%StressU) 
      End If
      If (AppCtx%VarFracSchemeParam%SaveStrain) Then
         Call Write_EXO_Result_Cell(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%CellVariable(VarFrac_CellVar_StrainXX)%Offset, AppCtx%TimeStep, AppCtx%StrainU) 
      End If
      Call Write_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_ElasticEnergy)%Offset, AppCtx%TimeStep, AppCtx%ElasticEnergy)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine Save
   
   Subroutine ElastFinalize(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx

      PetscInt                                     :: iErr
      Character(len=MEF90_MXSTRLEN)                :: filename
   

      If (AppCtx%VarFracSchemeParam%U_UseTao) Then
         Call TaoDestroy(AppCtx%taoU, iErr); CHKERRQ(iErr)
         Call TaoApplicationDestroy(AppCtx%taoappU, iErr); CHKERRQ(iErr)
      Else 
         Call KSPDestroy(AppCtx%KSPU, iErr); CHKERRQ(iErr)
      End If
      
      Call VecDestroy(AppCtx%U_Vec, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(AppCtx%U, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(AppCtx%F, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(AppCtx%Theta, iErr); CHKERRQ(iErr)
      If ( (AppCtx%VarFracSchemeParam%SaveStrain) .OR. (AppCtx%VarFracSchemeParam%SaveStress) ) Then
         Call SectionRealDestroy(AppCtx%StrainU, iErr); CHKERRQ(iErr)
         Call SectionRealDestroy(AppCtx%StressU, iErr); CHKERRQ(iErr)
      End If
      
      Call VecScatterDestroy(AppCtx%ScatterVect, iErr); CHKERRQ(iErr)
      Call VecScatterDestroy(AppCtx%ScatterScal, iErr); CHKERRQ(iErr)
      Call SectionIntDestroy(AppCtx%BCFlagU, iErr); CHKERRQ(iErr)
      Call MatDestroy(AppCtx%KU, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(AppCtx%RHSU, iErr); CHKERRQ(iErr)
      Call MeshDestroy(AppCtx%MeshTopology%Mesh, iErr); CHKERRQ(ierr)

      If (AppCtx%AppParam%verbose > 1) Then
         Call PetscViewerFlush(AppCtx%AppParam%MyLogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerFlush(AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerDestroy(AppCtx%AppParam%MyLogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerDestroy(AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
      End If
      Call PetscViewerDestroy(AppCtx%AppParam%EnergyViewer, iErr); CHKERRQ(iErr)
      Write(filename, 103) Trim(AppCtx%AppParam%prefix)
      Call PetscLogPrintSummary(PETSC_COMM_WORLD, filename, iErr); CHKERRQ(iErr)
      Call TaoFinalize(iErr); CHKERRQ(iErr)
      Call MEF90_Finalize()
103 Format(A,'-logsummary.txt')
   End Subroutine ElastFinalize
   
   
#if defined PB_2D
End Module m_Elast2D
#elif defined PB_3D
End Module m_Elast3D
#endif
