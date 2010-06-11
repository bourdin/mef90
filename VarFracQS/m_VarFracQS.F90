#if defined PB_2D
Module m_VarFracQS2D
#elif defined PB_3D
Module m_VarFracQS3D
#endif
#include "finclude/petscdef.h"

#if defined PB_2D
   Use m_VarFracQS_Types2D
   Use m_VarFracQS_U2D
   Use m_VarFracQS_V2D
   Use m_VarFracQS_Post2D
#elif defined PB_3D
   Use m_VarFracQS_Types3D   
   Use m_VarFracQS_U3D
   Use m_VarFracQS_V3D
   Use m_VarFracQS_Post3D
#endif   
   Use m_MEF90
   Use m_VarFrac_Struct

   Implicit NONE   
   
Contains

   Subroutine VarFracQSInit(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx

      PetscInt                                     :: iErr, i, iBlk
      Type(Mesh)                                   :: Tmp_Mesh
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer, filename
      PetscInt                                     :: NumComponents
      PetscTruth                                   :: HasPrefix
      PetscReal                                    :: rDummy
      Character                                    :: cDummy
      PetscInt                                     :: vers
      PetscTruth                                   :: flag
      
      PetscReal                                    :: KSP_Default_rtol
      PetscReal                                    :: KSP_Default_atol
      PetscInt                                     :: KSP_Default_MaxIt
      PetscReal                                    :: TAO_Default_fatol
      PetscReal                                    :: TAO_Default_frtol
      PetscReal                                    :: TAO_Default_gatol
      PetscReal                                    :: TAO_Default_grtol
      PetscReal                                    :: TAO_Default_gttol
      PetscReal                                    :: TAO_Default_catol
      PetscReal                                    :: TAO_Default_crtol
      Type(PetscViewer)                            :: flgviewer
      PetscInt, Dimension(:), Pointer              :: SizeVect, SizeScal
      
      Call MEF90_Initialize()
#if defined WITH_TAO
      Call TaoInitialize(PETSC_NULL_CHARACTER, iErr); CHKERRQ(iErr)
#endif
      Call PetscMemorySetGetMaximumUsage(iErr); CHKERRQ(iErr)

      KSP_Default_rtol  = 1.0D-6
      KSP_Default_atol  = 1.0D-6
      KSP_Default_MaxIt = 50000
      TAO_Default_fatol = 1.0D-10
      TAO_Default_frtol = 1.0D-8
      TAO_Default_gatol = 0.
      TAO_Default_grtol = 0.
      TAO_Default_gttol = 0.
      TAO_Default_catol = 0.
      TAO_Default_crtol = 0.
      
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-verbose', AppCtx%AppParam%Verbose, flag, iErr); CHKERRQ(iErr)
      Call InitFileNames(AppCtx)    
      AppCtx%AppParam%StopOnError = PETSC_FALSE
      Call PetscOptionsGetTruth(PETSC_NULL_CHARACTER, '-stop_on_error', AppCtx%AppParam%StopOnError, flag, iErr) ; CHKERRQ(iErr)

      Call VarFracSchemeParam_GetFromOptions(AppCtx%VarFracSchemeParam)
      If (AppCtx%AppParam%verbose > 0) Then
         Call VarFracSchemeParam_View(AppCtx%VarFracSchemeParam, PetscViewer(PETSC_VIEWER_STDOUT_WORLD))
      End If
      Call InitLog(AppCtx)
      Call PetscLogStagePush(AppCtx%LogInfo%Setup_Stage, iErr); CHKERRQ(iErr)
      
      Write(filename, 100) Trim(AppCtx%AppParam%prefix)
      Call PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, flgviewer, iErr); CHKERRQ(iErr);   
      Call VarFracSchemeParam_View(AppCtx%VarFracSchemeParam, flgviewer)
      Call PetscViewerFlush(flgviewer, iErr); CHKERRQ(iErr)
      Call PetscViewerDestroy(flgviewer, iErr); CHKERRQ(iErr)
      
      If (AppCtx%AppParam%verbose > 0) Then
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
100 Format(A, '.flg')      
101 Format(A, '-', I4.4, '.log')
102 Format('Output from processor ', I4.4, ' redirected to file ', A, '\n')
103 Format(A,'.log')
104 Format('Collective output redirected to file ', A, '\n')

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
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) "Done with VarFracQSEXOVariable_Init and VarFracQSEXOProperty_Read\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      If (AppCtx%AppParam%verbose > 1) Then
         Call MeshTopologyView(AppCtx%MeshTopology, AppCtx%AppParam%MyLogViewer)
         Call EXOView(AppCtx%MyEXO, AppCtx%AppParam%MyLogViewer)
      End If
      
      !!! Read Mat Properties from the CST file
      Call MatProp_Read(AppCtx%MeshTopology, AppCtx%MatProp, AppCtx%AppParam%CST_FileName)
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) "Done with MatProp_Read\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      If (AppCtx%AppParam%verbose > 2) Then
         Call EXOView(AppCtx%MyEXO, AppCtx%AppParam%MyLogViewer)
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

      !!! Create the Fields for the variables
      Allocate(SizeVect(AppCtx%MeshTopology%Num_dim))
      SizeVect= 1
      Allocate(SizeScal(1))
      SizeScal=1
      Call FieldCreateVertex(AppCtx%U,      'U',      AppCtx%MeshTopology, SizeVect)
      Call FieldCreateVertex(AppCtx%UBC,    'UBC',    AppCtx%MeshTopology, SizeVect)
      Call FieldCreateVertex(AppCtx%F,      'F',      AppCtx%MeshTopology, SizeVect)
      Call FieldCreateVertex(AppCtx%V,      'V',      AppCtx%MeshTopology, SizeScal)
      Call FieldCreateVertex(AppCtx%VBC,    'VBC',    AppCtx%MeshTopology, SizeScal)
      Call FieldCreateVertex(AppCtx%VIrrev, 'VIrrev', AppCtx%MeshTopology, SizeScal)
      Call FieldCreateVertex(AppCtx%Theta,  'Theta',  AppCtx%MeshTopology, SizeScal)
      If (AppCtx%VarFracSchemeParam%U_UseTao) Then
         Call FieldCreateVertex(AppCtx%GradientU,   'GradientU',   AppCtx%MeshTopology, SizeVect)
         Call FieldCreateVertex(AppCtx%LowerBoundU, 'LowerBoundU', AppCtx%MeshTopology, SizeVect)
         Call FieldCreateVertex(AppCtx%UpperBoundU, 'UpperBoundU', AppCtx%MeshTopology, SizeVect)
      Else
         Call FieldCreateVertex(AppCtx%RHSU, 'RHSU', AppCtx%MeshTopology, SizeVect)
      End If
      If (AppCtx%VarFracSchemeParam%V_UseTao) Then
         Call FieldCreateVertex(AppCtx%GradientV,   'GradientV',   AppCtx%MeshTopology, SizeScal)
         Call FieldCreateVertex(AppCtx%LowerBoundV, 'LowerBoundV', AppCtx%MeshTopology, SizeScal)
         Call FieldCreateVertex(AppCtx%UpperBoundV, 'UpperBoundV', AppCtx%MeshTopology, SizeScal)
      Else
         Call FieldCreateVertex(AppCtx%RHSV, 'RHSV', AppCtx%MeshTopology, SizeScal)
      End If
      Call FlagCreateVertex(AppCtx%BCUFlag,   'BCUFlag',   AppCtx%MeshTopology, SizeVect)
      Call FlagCreateVertex(AppCtx%BCVFlag,   'BCVFlag',   AppCtx%MeshTopology, SizeScal)
      Call FlagCreateVertex(AppCtx%IrrevFlag, 'IrrevFlag', AppCtx%MeshTopology, SizeScal)
      DeAllocate(SizeVect)
      DeAllocate(SizeScal)
      Call VecDuplicate(AppCtx%V%Vec, AppCtx%V_Old, iErr); CHKERRQ(iErr)

      If ( (AppCtx%VarFracSchemeParam%SaveStress) .OR. (AppCtx%VarFracSchemeParam%SaveStrain) ) Then
         NumComponents = AppCtx%MeshTopology%Num_Dim * (AppCtx%MeshTopology%Num_Dim + 1) / 2
         Call MeshGetCellSectionReal(AppCtx%MeshTopology%mesh, 'Strain', NumComponents, AppCtx%StrainU, iErr); CHKERRQ(iErr)
         Call MeshGetCellSectionReal(AppCtx%MeshTopology%mesh, 'Stress', NumComponents, AppCtx%StressU, iErr); CHKERRQ(iErr)
      End If

      !!! Create the Mat, KSP, PC
      Call MeshSetMaxDof(AppCtx%MeshTopology%Mesh, AppCtx%MeshTopology%Num_Dim, iErr); CHKERRQ(iErr) 
      Call MeshCreateMatrix(AppCtx%MeshTopology%mesh, AppCtx%U%Sec, MATMPIAIJ, AppCtx%KU, iErr); CHKERRQ(iErr)
      Call MeshCreateMatrix(AppCtx%MeshTopology%mesh, AppCtx%V%Sec, MATMPIAIJ, AppCtx%KV, iErr); CHKERRQ(iErr)
      
      !! Solver context for U      
      If (AppCtx%VarFracSchemeParam%U_UseTao) Then
#if defined WITH_TAO
         Call TaoCreate(PETSC_COMM_WORLD, 'tao_tron', AppCtx%taoU, iErr); CHKERRQ(iErr)
         Call TaoApplicationCreate(PETSC_COMM_WORLD, AppCtx%taoappU, iErr); CHKERRQ(iErr)
         Call TaoAppendOptionsPrefix(AppCtx%taoU, "U_", iErr); CHKERRQ(iErr)

         Call TaoAppSetObjectiveAndGradientRoutine(AppCtx%taoappU, FormFunctionAndGradientU, AppCtx, iErr); CHKERRQ(iErr)
         Call TaoAppSetHessianRoutine(AppCtx%taoappU, HessianU_Assembly, AppCtx, iErr); CHKERRQ(iErr)
         Call TaoAppSetVariableBoundsRoutine(AppCtx%taoappU, InitTaoBoundsU, AppCtx, iErr); CHKERRQ(iErr)

         Call TaoAppSetHessianMat(AppCtx%taoappU, AppCtx%KU, AppCtx%KU, iErr); CHKERRQ(iErr)

         Call TaoAppSetDefaultSolutionVec(AppCtx%taoappU, AppCtx%U%Vec, iErr); CHKERRQ(iErr)
         
         Call TaoSetOptions(AppCtx%taoappU, AppCtx%taoU, iErr); CHKERRQ(iErr)
         Call TaoAppSetFromOptions(AppCtx%taoappU, iErr); CHKERRQ(iErr)
         Call TaoSetTolerances(AppCtx%taoU, TAO_Default_fatol, TAO_Default_frtol, TAO_Default_catol, TAO_Default_crtol, iErr); CHKERRQ(iErr)
         Call TaoSetGradientTolerances(AppCtx%taoU, TAO_Default_gatol, TAO_Default_grtol, TAO_Default_gttol, iErr); CHKERRQ(iErr)
         Call TaoAppGetKSP(AppCtx%taoappU, AppCtx%KSPU, iErr); CHKERRQ(iErr)
         Call KSPSetType(AppCtx%KSPU, KSPSTCG, iErr); CHKERRQ(iErr)
         Call PetscPrintf(PETSC_COMM_WORLD, "TAO U Solver:\n", iErr); CHKERRQ(iErr)
         Call TaoView(AppCtx%taoU, iErr); CHKERRQ(iErr)
#endif
      Else
         Call KSPCreate(PETSC_COMM_WORLD, AppCtx%KSPU, iErr); CHKERRQ(iErr)
         Call KSPSetOperators(AppCtx%KSPU, AppCtx%KU, AppCtx%KU, SAME_NONZERO_PATTERN, iErr); CHKERRQ(iErr)
         Call KSPSetType(AppCtx%KSPU, KSPCG, iErr); CHKERRQ(iErr)
      End If
      Call KSPSetInitialGuessNonzero(AppCtx%KSPU, PETSC_TRUE, iErr); CHKERRQ(iErr)
      Call KSPAppendOptionsPrefix(AppCtx%KSPU, "U_", iErr); CHKERRQ(iErr)
      Call KSPSetTolerances(AppCtx%KSPU, KSP_Default_rtol, KSP_Default_atol, PETSC_DEFAULT_DOUBLE_PRECISION, KSP_Default_MaxIt, iErr)
      Call KSPSetFromOptions(AppCtx%KSPU, iErr); CHKERRQ(iErr)

      Call KSPGetPC(AppCtx%KSPU, AppCtx%PCU, iErr); CHKERRQ(iErr)
      Call PCSetType(AppCtx%PCU, PCBJACOBI, iErr); CHKERRQ(iErr)
      Call PCSetFromOptions(AppCtx%PCU, iErr); CHKERRQ(iErr)
!      Call KSPView(AppCtx%KSPU, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
         
      !! Solver context for V      
      If (AppCtx%VarFracSchemeParam%V_UseTao) Then
#if defined WITH_TAO
         Call TaoCreate(PETSC_COMM_WORLD, 'tao_tron', AppCtx%taoV, iErr); CHKERRQ(iErr)
         Call TaoApplicationCreate(PETSC_COMM_WORLD, AppCtx%taoappV, iErr); CHKERRQ(iErr)
         Call TaoAppendOptionsPrefix(AppCtx%taoV, "V_", iErr); CHKERRQ(iErr)

         Call TaoAppSetObjectiveAndGradientRoutine(AppCtx%taoappV, FormFunctionAndGradientV, AppCtx, iErr); CHKERRQ(iErr)
         Call TaoAppSetHessianRoutine(AppCtx%taoappV, HessianV_Assembly, AppCtx, iErr); CHKERRQ(iErr)
         Call TaoAppSetVariableBoundsRoutine(AppCtx%taoappV, InitTaoBoundsV, AppCtx, iErr); CHKERRQ(iErr)

         Call TaoAppSetHessianMat(AppCtx%taoappV, AppCtx%KV, AppCtx%KV, iErr); CHKERRQ(iErr)

         Call TaoAppSetDefaultSolutionVec(AppCtx%taoappV, AppCtx%V%Vec, iErr); CHKERRQ(iErr)
         
         Call TaoSetOptions(AppCtx%taoappV, AppCtx%taoV, iErr); CHKERRQ(iErr)
         Call TaoAppSetFromOptions(AppCtx%taoappV, iErr); CHKERRQ(iErr)
         Call TaoSetTolerances(AppCtx%taoV, TAO_Default_fatol, TAO_Default_frtol, TAO_Default_catol, TAO_Default_crtol, iErr); CHKERRQ(iErr)
         Call TaoSetGradientTolerances(AppCtx%taoV, TAO_Default_gatol, TAO_Default_grtol, TAO_Default_gttol, iErr); CHKERRQ(iErr)
         Call TaoAppGetKSP(AppCtx%taoappV, AppCtx%KSPV, iErr); CHKERRQ(iErr)
         Call KSPSetType(AppCtx%KSPV, KSPSTCG, iErr); CHKERRQ(iErr)
         Call PetscPrintf(PETSC_COMM_WORLD, "TAO V Solver:\n", iErr); CHKERRQ(iErr)
         Call TaoView(AppCtx%taoV, iErr); CHKERRQ(iErr)
#endif
      Else
         Call KSPCreate(PETSC_COMM_WORLD, AppCtx%KSPV, iErr); CHKERRQ(iErr)
         Call KSPSetOperators(AppCtx%KSPV, AppCtx%KV, AppCtx%KV, SAME_NONZERO_PATTERN, iErr); CHKERRQ(iErr)
         Call KSPSetType(AppCtx%KSPV, KSPCG, iErr); CHKERRQ(iErr)
      End If
      Call KSPSetInitialGuessNonzero(AppCtx%KSPV, PETSC_TRUE, iErr); CHKERRQ(iErr)
      Call KSPAppendOptionsPrefix(AppCtx%KSPV, "V_", iErr); CHKERRQ(iErr)
      Call KSPSetTolerances(AppCtx%KSPV, KSP_Default_rtol, KSP_Default_atol, PETSC_DEFAULT_DOUBLE_PRECISION, KSP_Default_MaxIt, iErr)
      Call KSPSetFromOptions(AppCtx%KSPV, iErr); CHKERRQ(iErr)

      Call KSPGetPC(AppCtx%KSPV, AppCtx%PCV, iErr); CHKERRQ(iErr)
      Call PCSetType(AppCtx%PCV, PCBJACOBI, iErr); CHKERRQ(iErr)
      Call PCSetFromOptions(AppCtx%PCV, iErr); CHKERRQ(iErr)
!      Call KSPView(AppCtx%KSPV, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)

      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) "Done Creating fields Section, Vec, KSP and Mat\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If

      !!! Initialize flags
      Call SectionIntZero(AppCtx%BCUFlag%Sec, iErr); CHKERRQ(iErr)
      Call SectionIntAddNSProperty(AppCtx%BCUFlag%Component_Sec(1), AppCtx%MyEXO%NSProperty(VarFrac_NSProp_BCUTypeX), AppCtx%MeshTopology)
      Call SectionIntAddNSProperty(AppCtx%BCUFlag%Component_Sec(2), AppCtx%MyEXO%NSProperty(VarFrac_NSProp_BCUTypeY), AppCtx%MeshTopology)
      If (AppCtx%MeshTopology%Num_Dim == 3) Then
         Call SectionIntAddNSProperty(AppCtx%BCUFlag%Component_Sec(3), AppCtx%MyEXO%NSProperty(VarFrac_NSProp_BCUTypeZ), AppCtx%MeshTopology)
      End If

      Call SectionIntZero(AppCtx%BCVFlag%Sec, iErr); CHKERRQ(iErr)
      Call SectionIntAddNSProperty(AppCtx%BCVFlag%Sec, AppCtx%MyEXO%NSProperty(VarFrac_NSProp_BCVType), AppCtx%MeshTopology)

      Call SectionIntZero(AppCtx%IrrevFlag%Sec, iErr); CHKERRQ(iErr)

      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) "Done Initializing BC Sections\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If

      !!! Get the number of time steps
      AppCtx%MyEXO%exoid = EXOPEN(AppCtx%MyEXO%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers, iErr)
      Call EXINQ(AppCtx%MyEXO%exoid, EXTIMS, AppCtx%NumTimeSteps, rDummy, cDummy, iErr)
      Call EXCLOS(AppCtx%MyEXO%exoid, iErr)
      AppCtx%MyEXO%exoid = 0
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Total Number of Time Steps', AppCtx%NumTimeSteps, '\n'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      
      !!! Allocate energies, and open matching files
      Allocate(AppCtx%SurfaceEnergy(AppCtx%NumTimeSteps))
      Allocate(AppCtx%ElasticEnergy(AppCtx%NumTimeSteps))
      Allocate(AppCtx%ExtForcesWork(AppCtx%NumTimeSteps))
      Allocate(AppCtx%TotalEnergy(AppCtx%NumTimeSteps))
      Allocate(AppCtx%Load(AppCtx%NumTimeSteps))
      
      Allocate(AppCtx%SurfaceEnergyBlock(AppCtx%MeshTopology%Num_Elem_Blks_Global))
      Allocate(AppCtx%ElasticEnergyBlock(AppCtx%MeshTopology%Num_Elem_Blks_Global))
      Allocate(AppCtx%ExtForcesWorkBlock(AppCtx%MeshTopology%Num_Elem_Blks_Global))
      Allocate(AppCtx%TotalEnergyblock  (AppCtx%MeshTopology%Num_Elem_Blks_Global))
      
      If (MEF90_MyRank == 0) Then
         AppCtx%AppParam%Ener_Unit = 71
         Allocate(AppCtx%AppParam%EnerBlock_Unit(AppCtx%MeshTopology%Num_Elem_Blks_Global))
         Open(File = AppCtx%AppParam%Ener_FileName, Unit = AppCtx%AppParam%Ener_Unit, Status = 'Unknown')
         Rewind(AppCtx%AppParam%Ener_Unit)
         Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks_Global
            AppCtx%AppParam%EnerBlock_Unit(iBlk) = 170+iBlk
            Write(AppCtx%AppParam%EnerBlock_FileName, 110) Trim(AppCtx%AppParam%Prefix), iBlk, Trim(AppCtx%AppParam%EnerBlock_Suffix)
            Open(File = AppCtx%AppParam%EnerBlock_FileName, Unit = AppCtx%AppParam%EnerBlock_Unit(iBlk), Status = 'Unknown')
            Rewind(AppCtx%AppParam%EnerBlock_Unit(iBlk))
         End Do
 110 Format(A, '-', I4.4, '.', A)
      End If
      
      !!! Set V=1
      Call SectionRealSet(AppCtx%VIrrev%Sec, 0.0_Kr, iErr); CHKERRQ(iErr)
      Call SectionRealSet(AppCtx%V%Sec, 1.0_Kr, iErr); CHKERRQ(iErr)
      Call VecSet(AppCtx%V%Vec, 1.0_Kr, iErr); CHKERRQ(iErr)
      
      !!! We are not backTracking
      AppCtx%IsBT = PETSC_FALSE

      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine VarFracQSInit
   
   Subroutine InitFileNames(dAppCtx) 
      Type(AppCtx_Type)                            :: dAppCtx
      PetscInt                                     :: iErr
      PetscTruth                                   :: HasPrefix
      Character(len=MEF90_MXSTRLEN)                :: tmpStr
      Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p', dAppCtx%AppParam%prefix, HasPrefix, iErr); CHKERRQ(iErr)
      If (.NOT. HasPrefix) Then
         Call PetscPrintf(PETSC_COMM_WORLD, "No input file prefix given\n", iErr)
         Call MEF90_Finalize()
         STOP
      End If

      !!!
      !!! CST file
      !!!
      Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-cst', dAppCtx%AppParam%CST_FileName, HasPrefix, iErr); CHKERRQ(iErr)
      If (.NOT. HasPrefix) Then
         dAppCtx%AppParam%CST_FileName = trim(dAppCtx%AppParam%prefix)//'.CST'
      End If

      !!! 
      !!! global energy file
      !!!
      Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-ener', dAppCtx%AppParam%Ener_FileName, HasPrefix, iErr); CHKERRQ(iErr)
      If (.NOT. HasPrefix) Then
         dAppCtx%AppParam%Ener_FileName = trim(dAppCtx%AppParam%prefix)//'.ener'
      End If
      
      !!!
      !!! blockwise energy files
      !!!
      Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-enerblk', dAppCtx%AppParam%EnerBlock_Suffix, HasPrefix, iErr); CHKERRQ(iErr)
      If (.NOT. HasPrefix) Then
         dAppCtx%AppParam%EnerBlock_Suffix = 'enerblk'
      End If
      
   End Subroutine InitFileNames




   Subroutine InitLog(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr
      
      If (AppCtx%VarFracSchemeParam%U_UseTao) Then
         Call PetscLogEventRegister('Hessian Local U',  0, AppCtx%LogInfo%MatAssemblyLocalU_Event, ierr); CHKERRQ(ierr)
         Call PetscLogEventRegister('Obj+Grad Local U', 0, AppCtx%LogInfo%RHSAssemblyLocalU_Event, ierr); CHKERRQ(ierr)
      Else
         Call PetscLogEventRegister('MatAssembly Local U', 0, AppCtx%LogInfo%MatAssemblyLocalU_Event, ierr); CHKERRQ(ierr)
         Call PetscLogEventRegister('RHSAssembly Local U', 0, AppCtx%LogInfo%RHSAssemblyLocalU_Event, ierr); CHKERRQ(ierr)
      End If
      If (AppCtx%VarFracSchemeParam%V_UseTao) Then
         Call PetscLogEventRegister('Hessian Local V',  0, AppCtx%LogInfo%MatAssemblyLocalV_Event, ierr); CHKERRQ(ierr)
         Call PetscLogEventRegister('Obj+Grad Local V', 0, AppCtx%LogInfo%RHSAssemblyLocalV_Event, ierr); CHKERRQ(ierr)
      Else
         Call PetscLogEventRegister('MatAssembly Local V', 0, AppCtx%LogInfo%MatAssemblyLocalV_Event, ierr); CHKERRQ(ierr)
         Call PetscLogEventRegister('RHSAssembly Local V', 0, AppCtx%LogInfo%RHSAssemblyLocalV_Event, ierr); CHKERRQ(ierr)
      End If
      Call PetscLogEventRegister('Post Processing',     0, AppCtx%LogInfo%PostProc_Event,          ierr); CHKERRQ(ierr)

      Call PetscLogStageRegister("Setup",            AppCtx%LogInfo%Setup_Stage,            iErr)
      Call PetscLogStageRegister("MeshDistribute",   AppCtx%LogInfo%MeshDistribute_Stage,   iErr)
      If (AppCtx%VarFracSchemeParam%U_UseTao) Then
         Call PetscLogStageRegister("Hessian U",  AppCtx%LogInfo%MatAssemblyU_Stage,     iErr)
         Call PetscLogStageRegister("Obj+Grad U", AppCtx%LogInfo%RHSAssemblyU_Stage,     iErr)
      Else
         Call PetscLogStageRegister("Mat Assembly U",   AppCtx%LogInfo%MatAssemblyU_Stage,     iErr)
         Call PetscLogStageRegister("RHS Assembly U",   AppCtx%LogInfo%RHSAssemblyU_Stage,     iErr)
      End If      
      Call PetscLogStageRegister("U-Step",      AppCtx%LogInfo%UStep_Stage,        iErr)
      If (AppCtx%VarFracSchemeParam%V_UseTao) Then
         Call PetscLogStageRegister("Hessian V",  AppCtx%LogInfo%MatAssemblyV_Stage,     iErr)
         Call PetscLogStageRegister("Obj+Grad V", AppCtx%LogInfo%RHSAssemblyV_Stage,     iErr)
      Else
         Call PetscLogStageRegister("Mat Assembly V",   AppCtx%LogInfo%MatAssemblyV_Stage,     iErr)
         Call PetscLogStageRegister("RHS Assembly V",   AppCtx%LogInfo%RHSAssemblyV_Stage,     iErr)
      End If
      Call PetscLogStageRegister("V-Step",      AppCtx%LogInfo%VStep_Stage,        iErr)
      Call PetscLogStageRegister("IO Stage",         AppCtx%LogInfo%IO_Stage,               iErr)
      Call PetscLogStageRegister("Post Proc",        AppCtx%LogInfo%PostProc_Stage,         iErr)
   End Subroutine InitLog
   
   Subroutine Save_U(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr

      Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
      Call Write_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_DisplacementX)%Offset, AppCtx%TimeStep, AppCtx%U) 
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine Save_U
   
   Subroutine Save_V(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr

      Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
      Call Write_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Fracture)%Offset, AppCtx%TimeStep, AppCtx%V) 
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine Save_V

   Subroutine Save_StrainStress(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr
   
      Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
      If (AppCtx%VarFracSchemeParam%SaveStress) Then
         Call Write_EXO_Result_Cell(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%CellVariable(VarFrac_CellVar_StressXX)%Offset, AppCtx%TimeStep, AppCtx%StressU) 
      End If
      If (AppCtx%VarFracSchemeParam%SaveStrain) Then
         Call Write_EXO_Result_Cell(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%CellVariable(VarFrac_CellVar_StrainXX)%Offset, AppCtx%TimeStep, AppCtx%StrainU) 
      End If
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine Save_StrainStress

   Subroutine ComputeEnergies(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx

      Call ElasticEnergy_Assembly(AppCtx%ElasticEnergy(AppCtx%TimeStep), AppCtx%ElasticEnergyBlock, AppCtx)
      Call ExtForcesWork_Assembly(AppCtx%ExtForcesWork(AppCtx%TimeStep), AppCtx%ExtForcesWorkBlock, AppCtx)
      Call SurfaceEnergy_Assembly(AppCtx%SurfaceEnergy(AppCtx%TimeStep), AppCtx%SurfaceEnergyBlock, AppCtx)
      AppCtx%TotalEnergy(AppCtx%TimeStep) = AppCtx%ElasticEnergy(AppCtx%TimeStep) - AppCtx%ExtForcesWork(AppCtx%TimeStep) + AppCtx%SurfaceEnergy(AppCtx%TimeStep)
      AppCtx%TotalEnergyBlock = AppCtx%ElasticEnergyBlock - AppCtx%ExtForcesWorkBlock + AppCtx%SurfaceEnergyBlock
   End Subroutine ComputeEnergies

   Subroutine Save_Ener(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr, iBlk

      Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
      Call Write_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_SurfaceEnergy)%Offset, AppCtx%TimeStep, AppCtx%SurfaceEnergy(AppCtx%TimeStep))
      Call Write_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_ElasticEnergy)%Offset, AppCtx%TimeStep, AppCtx%ElasticEnergy(AppCtx%TimeStep))
      Call Write_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_ExtForcesWork)%Offset, AppCtx%TimeStep, AppCtx%ExtForcesWork(AppCtx%TimeStep))
      Call Write_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_TotalEnergy)%Offset, AppCtx%TimeStep, AppCtx%TotalEnergy(AppCtx%TimeStep))
      Call Write_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_Load)%Offset, AppCtx%TimeStep, AppCtx%Load(AppCtx%TimeStep))
      
      If (MEF90_MyRank == 0) Then
         Write(AppCtx%AppParam%Ener_Unit, 100) AppCtx%TimeStep, AppCtx%Load(AppCtx%TimeStep), AppCtx%ElasticEnergy(AppCtx%TimeStep), AppCtx%ExtForcesWork(AppCtx%TimeStep), AppCtx%SurfaceEnergy(AppCtx%TimeStep), AppCtx%TotalEnergy(AppCtx%TimeStep)
         Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks_Global
            Write(AppCtx%AppParam%EnerBlock_Unit(iBlk), 100) AppCtx%TimeStep, AppCtx%Load(AppCtx%TimeStep), AppCtx%ElasticEnergyBlock(iBlk), AppCtx%ExtForcesWorkBlock(iBlk), AppCtx%SurfaceEnergyBlock(iBlk), AppCtx%TotalEnergyBlock(iBlk)
         End Do    
      End If
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
100   Format(I6, 5(ES13.5,'  '))  
   End Subroutine Save_Ener
   
   Subroutine Init_TS_Loads(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr

      Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_ForceX)%Offset, AppCtx%TimeStep, AppCtx%F) 
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Temperature)%Offset, AppCtx%TimeStep, AppCtx%Theta) 
      
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Fracture)%Offset, AppCtx%TimeStep, AppCtx%VBC) 
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_DisplacementX)%Offset, AppCtx%TimeStep, AppCtx%UBC) 

      Call Read_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_Load)%Offset, AppCtx%TimeStep, AppCtx%Load(AppCtx%TimeStep))
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine Init_TS_Loads   
   
   Subroutine VarFracQSFinalize(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx

      PetscInt                                     :: iErr, iBlk
      Character(len=MEF90_MXSTRLEN)                :: filename
   
      Call PetscLogStagePush(AppCtx%LogInfo%Setup_Stage, iErr); CHKERRQ(iErr)
      Call FieldDestroy(AppCtx%U);      CHKERRQ(iErr)
      Call FieldDestroy(AppCtx%UBC);    CHKERRQ(iErr)
      Call FieldDestroy(AppCtx%V);      CHKERRQ(iErr)
      Call FieldDestroy(AppCtx%VBC);    CHKERRQ(iErr)
      Call FieldDestroy(AppCtx%VIrrev); CHKERRQ(iErr)
      Call FieldDestroy(AppCtx%F);      CHKERRQ(iErr)
      Call FieldDestroy(AppCtx%Theta);  CHKERRQ(iErr)
      If (AppCtx%VarFracSchemeParam%U_UseTao) Then
         Call FieldDestroy(AppCtx%GradientU);   CHKERRQ(iErr)
         Call FieldDestroy(AppCtx%LowerBoundU); CHKERRQ(iErr)
         Call FieldDestroy(AppCtx%UpperBoundU); CHKERRQ(iErr)
      Else
         Call FieldDestroy(AppCtx%RHSU); CHKERRQ(iErr)
      End If
      If (AppCtx%VarFracSchemeParam%V_UseTao) Then
         Call FieldDestroy(AppCtx%GradientV);   CHKERRQ(iErr)
         Call FieldDestroy(AppCtx%LowerBoundV); CHKERRQ(iErr)
         Call FieldDestroy(AppCtx%UpperBoundV); CHKERRQ(iErr)
      Else
         Call FieldDestroy(AppCtx%RHSV); CHKERRQ(iErr)
      End If
      Call VecDestroy(AppCtx%V_Old, iErr); CHKERRQ(iErr)

      If ( (AppCtx%VarFracSchemeParam%SaveStress) .OR. (AppCtx%VarFracSchemeParam%SaveStrain) ) Then
         Call SectionRealDestroy(AppCtx%StrainU, iErr); CHKERRQ(iErr)
         Call SectionRealDestroy(AppCtx%StressU, iErr); CHKERRQ(iErr)
      End If
      
      Call FlagDestroy(AppCtx%BCUFlag);   CHKERRQ(iErr)
      Call FlagDestroy(AppCtx%BCVFlag);   CHKERRQ(iErr)
      Call FlagDestroy(AppCtx%IrrevFlag); CHKERRQ(iErr)

      Call MatDestroy(AppCtx%KU, iErr); CHKERRQ(iErr)
      Call MatDestroy(AppCtx%KV, iErr); CHKERRQ(iErr)
      Call MeshDestroy(AppCtx%MeshTopology%Mesh, iErr); CHKERRQ(ierr)

      DeAllocate(AppCtx%SurfaceEnergy)
      DeAllocate(AppCtx%ElasticEnergy)
      DeAllocate(AppCtx%ExtForcesWork)
      DeAllocate(AppCtx%TotalEnergy)
      DeAllocate(AppCtx%Load)
      DeAllocate(AppCtx%SurfaceEnergyBlock)
      DeAllocate(AppCtx%ElasticEnergyBlock)
      DeAllocate(AppCtx%ExtForcesWorkBlock)
      DeAllocate(AppCtx%TotalEnergyBlock)
      If (AppCtx%AppParam%verbose > 1) Then
         Call PetscViewerFlush(AppCtx%AppParam%MyLogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerFlush(AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerDestroy(AppCtx%AppParam%MyLogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerDestroy(AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
      End If
      
      If (MEF90_MyRank == 0) Then
         Close(AppCtx%AppParam%Ener_Unit)
         Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks_Global
            Close(AppCtx%AppParam%EnerBlock_Unit(iBlk))
         End Do
         DeAllocate(AppCtx%AppParam%EnerBlock_Unit)
      End If
      
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
      Write(filename, 103) Trim(AppCtx%AppParam%prefix)
      Call PetscLogPrintSummary(PETSC_COMM_WORLD, filename, iErr); CHKERRQ(iErr)
      
      If (AppCtx%VarFracSchemeParam%U_UseTao) Then
#if defined WITH_TAO
         Call TaoDestroy(AppCtx%taoU, iErr); CHKERRQ(iErr)
         Call TaoApplicationDestroy(AppCtx%taoAppU, iErr); CHKERRQ(iErr)
#endif
      Else
         Call KSPDestroy(AppCtx%KSPU, iErr); CHKERRQ(iErr)
      End If
      If (AppCtx%VarFracSchemeParam%V_UseTao) Then
#if defined WITH_TAO
         Call TaoDestroy(AppCtx%taoV, iErr); CHKERRQ(iErr)
         Call TaoApplicationDestroy(AppCtx%taoAppV, iErr); CHKERRQ(iErr)
#endif
      Else
         Call KSPDestroy(AppCtx%KSPV, iErr); CHKERRQ(iErr)
      End If

#if defined WITH_TAO
      Call TaoFinalize(iErr)
#endif
      Call MEF90_Finalize()
103 Format(A,'-logsummary.txt')
   End Subroutine VarFracQSFinalize


   ! Bactracking subroutine
   !    - Assumes that the energies have been computed
   !    - Returns iBTStep
   Subroutine BackTracking(AppCtx, iBTStep)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt, Intent(OUT)                        :: iBTStep
      
      PetscInt                                     :: iErr
      PetscReal                                    :: EnerBT, EnerRef
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
   
      !!! Check the BT condition
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Doing BackTracking\n' 
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
      End If
      Do iBTStep = max(1, AppCtx%TimeStep-AppCtx%VarFracSchemeParam%BTScope), AppCtx%TimeStep-1
         EnerBT  = AppCtx%Load(iBTStep)**2 * (AppCtx%ElasticEnergy(AppCtx%TimeStep) - AppCtx%ExtForcesWork(AppCtx%TimeStep)) + AppCtx%Load(AppCtx%TimeStep)**2 * AppCtx%SurfaceEnergy(AppCtx%TimeStep)
         EnerRef = AppCtx%Load(AppCtx%TimeStep)**2 * (AppCtx%TotalEnergy(iBTStep) - AppCtx%ExtForcesWork(AppCtx%TimeStep))
         If (AppCtx%AppParam%verbose > 0) Then
            Write(IOBuffer, *) 'Checking against timestep', iBTStep, ':', EnerBT, EnerRef, (1.0_Kr - AppCtx%VarFracSchemeParam%BTTol) * EnerRef, '\n'
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
         End If
         
         If (EnerBT < (1.0_Kr - AppCtx%VarFracSchemeParam%BTTol) * EnerRef) Then
            If (AppCtx%AppParam%verbose > 0) Then
               Write(IOBuffer, *) 'BackTracking to step', iBTStep, '\n' 
               Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
            End If
            EXIT
         End If
      End Do
   End Subroutine BackTracking   
   
#if defined PB_2D
End Module m_VarFracQS2D
#elif defined PB_3D
End Module m_VarFracQS3D
#endif
