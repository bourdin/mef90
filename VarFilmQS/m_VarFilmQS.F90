Module m_VarFilmQS
#include "finclude/petscdef.h"

   Use m_VarFilmQS_Types
   Use m_VarFilmQS_U
   Use m_VarFilmQS_NL_U
   Use m_VarFilmQS_V
   Use m_VarFilmQS_W
   Use m_VarFilmQS_Post
   
   Use m_MEF90
   Use m_Film_Struct

   Implicit NONE   
   
Contains

#undef __FUNC__ 
#define __FUNC__ "VarFilmQSInit"
   Subroutine VarFilmQSInit(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx

      PetscInt                                     :: iErr, i, iBlk, iTS
      Type(DM)                                   :: Tmp_Mesh
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer, filename
      PetscInt                                     :: NumComponents
      PetscBool                                    :: HasPrefix
      PetscReal                                    :: rDummy
      Character                                    :: cDummy
      PetscInt                                     :: vers
      PetscBool                                    :: flag
      
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
      Call PetscOptionsGetBool(PETSC_NULL_CHARACTER, '-stop_on_error', AppCtx%AppParam%StopOnError, flag, iErr) ; CHKERRQ(iErr)
      AppCtx%TimeStep = 1
      AppCtx%AppParam%Restart = PETSC_FALSE
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER, "-restart", AppCtx%TimeStep, AppCtx%AppParam%Restart, iErr); CHKERRQ(iErr)

      Call VarFracSchemeParam_GetFromOptions(AppCtx%VarFracSchemeParam)
      If (AppCtx%AppParam%verbose > 0) Then
         Call VarFracSchemeParam_View(AppCtx%VarFracSchemeParam, PetscViewer(PETSC_VIEWER_STDOUT_WORLD))
      End If
      Call InitLog(AppCtx)
      Call PetscLogStagePush(AppCtx%LogInfo%Setup_Stage, iErr); CHKERRQ(iErr)
      
      Write(filename, 100) Trim(AppCtx%AppParam%prefix)
      Call PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, flgviewer, iErr); CHKERRQ(iErr);          ! ???
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
      AppCtx%EXO%exoid = EXOPEN(AppCtx%EXO%filename, EXREAD, exo_cpu_ws, exo_io_ws, PETSC_NULL_INTEGER, ierr) 

      !!! Reading and distributing sequential mesh
      If (MEF90_NumProcs == 1) Then
         Call DMMeshCreateExodus(PETSC_COMM_WORLD, AppCtx%EXO%filename, AppCtx%MeshTopology%mesh, ierr); CHKERRQ(iErr)
      Else
         Call DMMeshCreateExodus(PETSC_COMM_WORLD, AppCtx%EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
         Call PetscLogStagePush(AppCtx%LogInfo%MeshDistribute_Stage, iErr); CHKERRQ(iErr)
         Call DMMeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, AppCtx%MeshTopology%mesh, ierr); CHKERRQ(iErr)
         Call PetscLogStagePop(iErr); CHKERRQ(iErr)
         Call DMDestroy(Tmp_mesh, ierr); CHKERRQ(iErr)
      End If
   
   Call MeshTopologyGetInfo(AppCtx%MeshTopology, PETSC_COMM_WORLD)
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) "Done reading and partitioning the mesh\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
   
      AppCtx%MyEXO%comm = PETSC_COMM_SELF
      AppCtx%MyEXO%exoid = AppCtx%EXO%exoid
      Write(AppCtx%MyEXO%filename, 99) trim(AppCtx%AppParam%prefix), MEF90_MyRank
      AppCtx%MyEXO%exoid = EXOPEN(AppCtx%MyEXO%filename, EXWRIT, exo_cpu_ws, exo_io_ws, PETSC_NULL_INTEGER, ierr) 
 99  Format(A, '-', I4.4, '.gen')

!!! Initializes the values and names of the properties and variables
      If ( (AppCtx%VarFracSchemeParam%SaveStress) .OR. (AppCtx%VarFracSchemeParam%SaveStrain) ) Then
         Call VarFracEXOVariable_Init(AppCtx%MyEXO,.TRUE.)
      Else
         Call VarFracEXOVariable_Init(AppCtx%MyEXO,.FALSE.)
      End If
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
	If (AppCtx%VarFracSchemeParam%U_UseSNES) Then
		Call FieldCreateVertex(AppCtx%GradientU,   'GradientU',   AppCtx%MeshTopology, SizeVect)
	End If

   Call FieldCreateVertex(AppCtx%U0,      'U0',    AppCtx%MeshTopology, SizeVect)
      
   Call FieldCreateVertex(AppCtx%V,      'V',      AppCtx%MeshTopology, SizeScal)
   Call FieldCreateVertex(AppCtx%VBC,    'VBC',    AppCtx%MeshTopology, SizeScal)
   Call FieldCreateVertex(AppCtx%VIrrev, 'VIrrev', AppCtx%MeshTopology, SizeScal)
   
   Call FieldCreateVertex(AppCtx%W,      'W',      AppCtx%MeshTopology, SizeScal)
   Call FieldCreateVertex(AppCtx%WBC,    'WBC',    AppCtx%MeshTopology, SizeScal)
   Call FieldCreateVertex(AppCtx%FW,     'FW',     AppCtx%MeshTopology, SizeScal)
   Call FieldCreateVertex(AppCtx%WIrrev, 'WIrrev', AppCtx%MeshTopology, SizeScal)
   
   Call FieldCreateVertex(AppCtx%Theta,  'Theta',  AppCtx%MeshTopology, SizeScal)
   Call FieldCreateVertex(AppCtx%RHSU,   'RHSU',   AppCtx%MeshTopology, SizeVect)

   If (AppCtx%VarFracSchemeParam%V_UseTao) Then
      Call FieldCreateVertex(AppCtx%GradientV,   'GradientV',   AppCtx%MeshTopology, SizeScal)
      Call FieldCreateVertex(AppCtx%LowerBoundV, 'LowerBoundV', AppCtx%MeshTopology, SizeScal)
      Call FieldCreateVertex(AppCtx%UpperBoundV, 'UpperBoundV', AppCtx%MeshTopology, SizeScal)
   Else
      Call FieldCreateVertex(AppCtx%RHSV, 'RHSV', AppCtx%MeshTopology, SizeScal)
   End If
   Call FlagCreateVertex(AppCtx%BCUFlag,   'BCUFlag',   AppCtx%MeshTopology, SizeVect)
   Call FlagCreateVertex(AppCtx%BCVFlag,   'BCVFlag',   AppCtx%MeshTopology, SizeScal)
   Call FlagCreateVertex(AppCtx%BCWFlag,   'BCWFlag',   AppCtx%MeshTopology, SizeScal)

	Call FlagCreateVertex(AppCtx%IrrevFlag, 'IrrevFlag', AppCtx%MeshTopology, SizeScal)
	Call FlagCreateVertex(AppCtx%WIrrevFlag, 'WIrrevFlag', AppCtx%MeshTopology, SizeScal)

      DeAllocate(SizeVect)
      DeAllocate(SizeScal)
      Call VecDuplicate(AppCtx%V%Vec, AppCtx%V_Old, iErr); CHKERRQ(iErr)
      Call VecDuplicate(AppCtx%U%Vec, AppCtx%U_Old, iErr); CHKERRQ(iErr)
      Call VecDuplicate(AppCtx%W%Vec, AppCtx%W_Old, iErr); CHKERRQ(iErr)

   If ( (AppCtx%VarFracSchemeParam%SaveStress) .OR. (AppCtx%VarFracSchemeParam%SaveStrain) ) Then
      NumComponents = AppCtx%MeshTopology%Num_Dim * (AppCtx%MeshTopology%Num_Dim + 1) / 2
      Call DMMeshGetCellSectionReal(AppCtx%MeshTopology%mesh, 'Strain', NumComponents, AppCtx%StrainU, iErr); CHKERRQ(iErr)
      Call DMMeshGetCellSectionReal(AppCtx%MeshTopology%mesh, 'Stress', NumComponents, AppCtx%StressU, iErr); CHKERRQ(iErr)
   End If

      !!! Create the Mat, KSP, PC
   Call DMMeshSetMaxDof(AppCtx%MeshTopology%Mesh, AppCtx%MeshTopology%Num_Dim, iErr); CHKERRQ(iErr) 
   Call DMMeshCreateMatrix(AppCtx%MeshTopology%mesh, AppCtx%U%Sec, MATMPIAIJ, AppCtx%KU, iErr); CHKERRQ(iErr)
   Call DMMeshCreateMatrix(AppCtx%MeshTopology%mesh, AppCtx%V%Sec, MATMPIAIJ, AppCtx%KV, iErr); CHKERRQ(iErr)


!!!startregion Solver context for U
If (AppCtx%VarFracSchemeParam%U_UseSNES) Then
! SNES Solver Ctx for U
	Call SNESCreate(PETSC_COMM_WORLD, AppCtx%snesU, iErr); CHKERRQ(iErr)
	Call SNESSetFunction(AppCtx%snesU, AppCtx%GradientU%Vec, GradientU_Assembly, AppCtx, iErr); CHKERRQ(iErr)
	Call SNESSetJacobian(AppCtx%snesU, AppCtx%KU, AppCtx%KU, HessianU_Assembly, AppCtx, iErr); CHKERRQ(iErr)
	Call SNESSetFromOptions(AppCtx%snesU, iErr); CHKERRQ(iErr)
	Call SNESGetKSP(AppCtx%snesU, AppCtx%KSPU, iErr); CHKERRQ(iErr)
	Call KSPSetType(AppCtx%KSPU, KSPCG, iErr); CHKERRQ(iErr)
	Call KSPSetFromOptions(AppCtx%KSPU, iErr); CHKERRQ(iErr)
	Call KSPGetPC(AppCtx%KSPU, AppCtx%PCU, iErr); CHKERRQ(iErr) 
	Call PCSetType(AppCtx%PCU, PCBJACOBI, iErr); CHKERRQ(iErr) 
	Call KSPSetTolerances(AppCtx%KSPU,1.e-4, KSP_Default_rtol, KSP_Default_atol, PETSC_DEFAULT_DOUBLE_PRECISION, KSP_Default_MaxIt, iErr); CHKERRQ(ierr);
Else

   Call KSPCreate(PETSC_COMM_WORLD, AppCtx%KSPU, iErr); CHKERRQ(iErr)
   Call KSPSetOperators(AppCtx%KSPU, AppCtx%KU, AppCtx%KU, SAME_NONZERO_PATTERN, iErr); CHKERRQ(iErr)
   Call KSPSetType(AppCtx%KSPU, KSPCG, iErr); CHKERRQ(iErr)
      Call KSPSetInitialGuessNonzero(AppCtx%KSPU, PETSC_TRUE, iErr); CHKERRQ(iErr)
      Call KSPAppendOptionsPrefix(AppCtx%KSPU, "U_", iErr); CHKERRQ(iErr)
      Call KSPSetTolerances(AppCtx%KSPU, KSP_Default_rtol, KSP_Default_atol, PETSC_DEFAULT_DOUBLE_PRECISION, KSP_Default_MaxIt, iErr)
      Call KSPSetFromOptions(AppCtx%KSPU, iErr); CHKERRQ(iErr)

      Call KSPGetPC(AppCtx%KSPU, AppCtx%PCU, iErr); CHKERRQ(iErr)
      Call PCSetType(AppCtx%PCU, PCBJACOBI, iErr); CHKERRQ(iErr)
      Call PCSetFromOptions(AppCtx%PCU, iErr); CHKERRQ(iErr)
!      Call KSPView(AppCtx%KSPU, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
End If
!!!endregion Solver context for U      
         
!!!startregion Solver context for V      
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
!!!endregion Solver context for V      

      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) "Done Creating fields Section, Vec, KSP and Mat\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If

!!!startregion BC Flags
      Call SectionIntZero(AppCtx%BCUFlag%Sec, iErr); CHKERRQ(iErr)
      Call SectionIntAddNSProperty(AppCtx%BCUFlag%Component_Sec(1), AppCtx%MyEXO%NSProperty(VarFrac_NSProp_BCUTypeX), AppCtx%MeshTopology)
      Call SectionIntAddNSProperty(AppCtx%BCUFlag%Component_Sec(2), AppCtx%MyEXO%NSProperty(VarFrac_NSProp_BCUTypeY), AppCtx%MeshTopology)


   Call SectionIntZero(AppCtx%BCVFlag%Sec, iErr); CHKERRQ(iErr)
   Call SectionIntAddNSProperty(AppCtx%BCVFlag%Sec, AppCtx%MyEXO%NSProperty(VarFrac_NSProp_BCVType), AppCtx%MeshTopology)

   Call SectionIntZero(AppCtx%BCWFlag%Sec, iErr); CHKERRQ(iErr)
   Call SectionIntAddNSProperty(AppCtx%BCWFlag%Sec, AppCtx%MyEXO%NSProperty(VarFrac_NSProp_BCWType), AppCtx%MeshTopology)

   Call SectionIntZero(AppCtx%IrrevFlag%Sec, iErr); CHKERRQ(iErr)

      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) "Done Initializing BC Sections\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
!!!endregion BC Flags

      !!! Get the number of time steps
      Call EXINQ(AppCtx%MyEXO%exoid, EXTIMS, AppCtx%NumTimeSteps, rDummy, cDummy, iErr)
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Total Number of Time Steps', AppCtx%NumTimeSteps, '\n'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If

      !!! Allocate energies, and open matching files
   Allocate(AppCtx%FractureEnergy(AppCtx%NumTimeSteps))
   Allocate(AppCtx%DelaminationEnergy(AppCtx%NumTimeSteps))
   Allocate(AppCtx%BondingLayerEnergy(AppCtx%NumTimeSteps))
   Allocate(AppCtx%FilmEnergy(AppCtx%NumTimeSteps))
   Allocate(AppCtx%ElasticEnergy(AppCtx%NumTimeSteps))
   Allocate(AppCtx%ExtForcesWork(AppCtx%NumTimeSteps))
   Allocate(AppCtx%TotalEnergy(AppCtx%NumTimeSteps))
   Allocate(AppCtx%Load(AppCtx%NumTimeSteps))
      
   Allocate(AppCtx%FractureEnergyBlock(AppCtx%MeshTopology%Num_Elem_Blks_Global))
   Allocate(AppCtx%DelaminationEnergyBlock(AppCtx%MeshTopology%Num_Elem_Blks_Global))
   Allocate(AppCtx%BondingLayerEnergyBlock(AppCtx%MeshTopology%Num_Elem_Blks_Global))
   Allocate(AppCtx%FilmEnergyBlock(AppCtx%MeshTopology%Num_Elem_Blks_Global))
   Allocate(AppCtx%ElasticEnergyBlock(AppCtx%MeshTopology%Num_Elem_Blks_Global))
   Allocate(AppCtx%ExtForcesWorkBlock(AppCtx%MeshTopology%Num_Elem_Blks_Global))
   Allocate(AppCtx%TotalEnergyblock  (AppCtx%MeshTopology%Num_Elem_Blks_Global))
      
      if_rank0: If (MEF90_MyRank == 0) Then
         AppCtx%AppParam%Ener_Unit = 71
         if_restart: If (AppCtx%AppParam%Restart) Then
            Open(File = AppCtx%AppParam%Ener_FileName, Unit = AppCtx%AppParam%Ener_Unit, Status = 'old', Position='Append')
            Rewind(AppCtx%AppParam%Ener_Unit)
            Do 
               Read(AppCtx%AppParam%Ener_Unit, *, end=999, err=999) iTS, rdummy, AppCtx%ElasticEnergy(iTS), AppCtx%ExtForcesWork(iTS), AppCtx%FractureEnergy(iTS), AppCtx%TotalEnergy(iTS)
               CYCLE
 999            EXIT               
            End Do
            Write(AppCtx%AppParam%Ener_Unit, *)
            Write(AppCtx%AppParam%Ener_Unit, *)
            Write(AppCtx%AppParam%Ener_Unit, *)
            If (AppCtx%TimeStep == 0) Then
               AppCtx%TimeStep = iTS+1
            End If
            If (AppCtx%AppParam%verbose > 0) Then
               Write(IOBuffer, *) 'Restarting from step', AppCtx%TimeStep, '\n'
               Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
            End If
         Else
            Open(File = AppCtx%AppParam%Ener_FileName, Unit = AppCtx%AppParam%Ener_Unit, Status = 'Unknown')
            Rewind(AppCtx%AppParam%Ener_Unit)
         End If if_restart
         
         if_saveblk: If (AppCtx%VarFracSchemeParam%SaveBlk) Then
            Allocate(AppCtx%AppParam%EnerBlock_Unit(AppCtx%MeshTopology%Num_Elem_Blks_Global))
            Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks_Global
               AppCtx%AppParam%EnerBlock_Unit(iBlk) = 170+iBlk
               Write(AppCtx%AppParam%EnerBlock_FileName, 110) Trim(AppCtx%AppParam%Prefix), iBlk, Trim(AppCtx%AppParam%EnerBlock_Suffix)
               If (AppCtx%AppParam%Restart) Then
                  Open(File = AppCtx%AppParam%EnerBlock_FileName, Unit = AppCtx%AppParam%EnerBlock_Unit(iBlk), Status = 'Old', Position='Append')
                  Write(AppCtx%AppParam%EnerBlock_Unit(iBlk),*)
                  Write(AppCtx%AppParam%EnerBlock_Unit(iBlk),*)
                  Write(AppCtx%AppParam%EnerBlock_Unit(iBlk),*)
               Else
                  Open(File = AppCtx%AppParam%EnerBlock_FileName, Unit = AppCtx%AppParam%EnerBlock_Unit(iBlk), Status = 'Unknown')
                  Rewind(AppCtx%AppParam%EnerBlock_Unit(iBlk))
               End If
            End Do
         End If if_saveblk
 110 Format(A, '-', I4.4, '.', A)
      End If if_rank0
      
      !!! Broacasting Energies in case we are restarting and it is needed for backtracking
      If (AppCtx%AppParam%Restart) Then
         Call MPI_BCast(AppCtx%TimeStep, 1, MPIU_INTEGER, 0, PETSC_COMM_WORLD, iErr);CHKERRQ(iErr)
         Do iTS = 1, AppCtx%TimeStep
            Call MPI_BCast(AppCtx%ElasticEnergy(iTS), 1, MPIU_SCALAR, 0, PETSC_COMM_WORLD, iErr);CHKERRQ(iErr)
            Call MPI_BCast(AppCtx%ExtForcesWork(iTS), 1, MPIU_SCALAR, 0, PETSC_COMM_WORLD, iErr);CHKERRQ(iErr)
            Call MPI_BCast(AppCtx%ElasticEnergy(iTS), 1, MPIU_SCALAR, 0, PETSC_COMM_WORLD, iErr);CHKERRQ(iErr)
            Call MPI_BCast(AppCtx%TotalEnergy(iTS), 1, MPIU_SCALAR, 0, PETSC_COMM_WORLD, iErr);CHKERRQ(iErr)
         End Do
      End If
      
      !!! Set V=1 if we are not restarting, and rekoad last time step if not
      Call SectionRealSet(AppCtx%WIrrev%Sec, 0.0_Kr, iErr); CHKERRQ(iErr)
      Call SectionRealSet(AppCtx%VIrrev%Sec, 0.0_Kr, iErr); CHKERRQ(iErr)

	If (AppCtx%AppParam%Restart .AND. (AppCtx%TimeSTep > 1) ) Then
		Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Fracture)%Offset, AppCtx%TimeStep-1, AppCtx%V)
		Call SectionRealToVec(AppCtx%V%Sec, AppCtx%V%Scatter, SCATTER_REVERSE, AppCtx%V%Vec, ierr); CHKERRQ(ierr)
		Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_DisplacementX)%Offset, AppCtx%TimeStep-1, AppCtx%U)
		Call SectionRealToVec(AppCtx%U%Sec, AppCtx%U%Scatter, SCATTER_REVERSE, AppCtx%U%Vec, ierr); CHKERRQ(ierr)
	Else
		Call SectionRealSet(AppCtx%V%Sec, 1.0_Kr, iErr); CHKERRQ(iErr)
		Call SectionRealSet(AppCtx%W%Sec, 1.0_Kr, iErr); CHKERRQ(iErr)
		Call VecSet(AppCtx%V%Vec, 1.0_Kr, iErr); CHKERRQ(iErr)
		Call VecSet(AppCtx%W%Vec, 1.0_Kr, iErr); CHKERRQ(iErr)
	End If
      
      !!! We are not backTracking
      AppCtx%IsBT = PETSC_FALSE

      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine VarFilmQSInit
   
#undef __FUNC__ 
#define __FUNC__ "InitFileNames"
   Subroutine InitFileNames(dAppCtx) 
      Type(AppCtx_Type)                            :: dAppCtx
      PetscInt                                     :: iErr
      PetscBool                                    :: HasPrefix
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

#undef __FUNC__ 
#define __FUNC__ "InitLog"
   Subroutine InitLog(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr
      
   Call PetscLogEventRegister('MatAssembly Local U',    0, AppCtx%LogInfo%MatAssemblyLocalU_Event, ierr); CHKERRQ(ierr)
   Call PetscLogEventRegister('RHSAssembly Local U',    0, AppCtx%LogInfo%RHSAssemblyLocalU_Event, ierr); CHKERRQ(ierr)
   If (AppCtx%VarFracSchemeParam%V_UseTao) Then
         Call PetscLogEventRegister('Hessian Local V',  0, AppCtx%LogInfo%MatAssemblyLocalV_Event, ierr); CHKERRQ(ierr)
         Call PetscLogEventRegister('Obj+Grad Local V', 0, AppCtx%LogInfo%RHSAssemblyLocalV_Event, ierr); CHKERRQ(ierr)
   Else
      Call PetscLogEventRegister('MatAssembly Local V', 0, AppCtx%LogInfo%MatAssemblyLocalV_Event, ierr); CHKERRQ(ierr)
      Call PetscLogEventRegister('RHSAssembly Local V', 0, AppCtx%LogInfo%RHSAssemblyLocalV_Event, ierr); CHKERRQ(ierr)
   End If
   Call PetscLogEventRegister('FWAssembly Local W',     0, AppCtx%LogInfo%FWAssemblyLocalW_Event, ierr); CHKERRQ(ierr)
   Call PetscLogEventRegister('Post Processing',        0, AppCtx%LogInfo%PostProc_Event,         ierr); CHKERRQ(ierr)

   Call PetscLogStageRegister("Setup",            AppCtx%LogInfo%Setup_Stage,            iErr)
   Call PetscLogStageRegister("MeshDistribute",   AppCtx%LogInfo%MeshDistribute_Stage,   iErr)
   Call PetscLogStageRegister("Mat Assembly U",   AppCtx%LogInfo%MatAssemblyU_Stage,     iErr)
   Call PetscLogStageRegister("RHS Assembly U",   AppCtx%LogInfo%RHSAssemblyU_Stage,     iErr)
   Call PetscLogStageRegister("U-Step",           AppCtx%LogInfo%UStep_Stage,            iErr)
   If (AppCtx%VarFracSchemeParam%V_UseTao) Then
      Call PetscLogStageRegister("Hessian V",     AppCtx%LogInfo%MatAssemblyV_Stage,     iErr)
      Call PetscLogStageRegister("Obj+Grad V",    AppCtx%LogInfo%RHSAssemblyV_Stage,     iErr)
   Else
      Call PetscLogStageRegister("Mat Assembly V",AppCtx%LogInfo%MatAssemblyV_Stage,     iErr)
      Call PetscLogStageRegister("RHS Assembly V",AppCtx%LogInfo%RHSAssemblyV_Stage,     iErr)
   End If
   Call PetscLogStageRegister("FW Assembly W",    AppCtx%LogInfo%FWAssemblyW_Stage,      iErr)
   Call PetscLogStageRegister("V-Step",           AppCtx%LogInfo%VStep_Stage,            iErr)
   Call PetscLogStageRegister("W-Step",           AppCtx%LogInfo%WStep_Stage,            iErr)
   Call PetscLogStageRegister("IO Stage",         AppCtx%LogInfo%IO_Stage,               iErr)
   Call PetscLogStageRegister("Post Proc",        AppCtx%LogInfo%PostProc_Stage,         iErr)
   End Subroutine InitLog
   
#undef __FUNC__ 
#define __FUNC__ "Save_U"
   Subroutine Save_U(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr

      Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
      Call Write_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_DisplacementX)%Offset, AppCtx%TimeStep, AppCtx%U) 
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine Save_U
   
#undef __FUNC__ 
#define __FUNC__ "Save_V"
   Subroutine Save_V(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr

      Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
      Call Write_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Fracture)%Offset, AppCtx%TimeStep, AppCtx%V) 
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine Save_V

#undef __FUNC__ 
#define __FUNC__ "Save_W"
Subroutine Save_W(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr
      
      Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
      Call Write_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Delamination)%Offset, AppCtx%TimeStep, AppCtx%W) 
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
End Subroutine Save_W


#undef __FUNC__ 
#define __FUNC__ "Save_StrainStress"
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

#undef __FUNC__ 
#define __FUNC__ "ComputeEnergies"
   Subroutine ComputeEnergies(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx

      Call ExtForcesWork_Assembly(AppCtx%ExtForcesWork(AppCtx%TimeStep), AppCtx%ExtForcesWorkBlock, AppCtx)
      Call FractureEnergy_Assembly(AppCtx%FractureEnergy(AppCtx%TimeStep), AppCtx%FractureEnergyBlock, AppCtx)
      Call DelaminationEnergy_Assembly(AppCtx%DelaminationEnergy(AppCtx%TimeStep), AppCtx%DelaminationEnergyBlock, AppCtx)
      Call BondingLayerEnergy_Assembly(AppCtx%BondingLayerEnergy(AppCtx%TimeStep), AppCtx%BondingLayerEnergyBlock, AppCtx)
      Call FilmEnergy_Assembly(AppCtx%FilmEnergy(AppCtx%TimeStep), AppCtx%FilmEnergyBlock, AppCtx)

      AppCtx%ElasticEnergy(AppCtx%TimeStep) = AppCtx%FilmEnergy(AppCtx%TimeStep) + AppCtx%BondingLayerEnergy(AppCtx%TimeStep)
      AppCtx%ElasticEnergyBlock = AppCtx%FilmEnergyBlock + AppCtx%BondingLayerEnergyBlock 

      AppCtx%TotalEnergy(AppCtx%TimeStep) = AppCtx%ElasticEnergy(AppCtx%TimeStep) - AppCtx%ExtForcesWork(AppCtx%TimeStep) + AppCtx%FractureEnergy(AppCtx%TimeStep) + AppCtx%DelaminationEnergy(AppCtx%TimeStep)
      AppCtx%TotalEnergyBlock = AppCtx%ElasticEnergyBlock - AppCtx%ExtForcesWorkBlock + AppCtx%FractureEnergyBlock + AppCtx%DelaminationEnergyBlock

   End Subroutine ComputeEnergies

#undef __FUNC__ 
#define __FUNC__ "Save_Ener"
   Subroutine Save_Ener(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr, iBlk

      Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
      Call Write_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_FractureEnergy)%Offset, AppCtx%TimeStep, AppCtx%FractureEnergy(AppCtx%TimeStep))
      Call Write_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_ElasticEnergy)%Offset, AppCtx%TimeStep, AppCtx%ElasticEnergy(AppCtx%TimeStep))
      Call Write_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_DelaminationEnergy)%Offset, AppCtx%TimeStep, AppCtx%DelaminationEnergy(AppCtx%TimeStep))
      Call Write_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_BondingLayerEnergy)%Offset, AppCtx%TimeStep, AppCtx%BondingLayerEnergy(AppCtx%TimeStep))
      Call Write_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_FilmEnergy)%Offset, AppCtx%TimeStep, AppCtx%FilmEnergy(AppCtx%TimeStep))
      Call Write_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_TotalEnergy)%Offset, AppCtx%TimeStep, AppCtx%TotalEnergy(AppCtx%TimeStep))
      Call Write_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_Load)%Offset, AppCtx%TimeStep, AppCtx%Load(AppCtx%TimeStep))
      
      If (MEF90_MyRank == 0) Then
         Write(AppCtx%AppParam%Ener_Unit, 100) AppCtx%TimeStep, AppCtx%Load(AppCtx%TimeStep), AppCtx%ElasticEnergy(AppCtx%TimeStep), AppCtx%ExtForcesWork(AppCtx%TimeStep), AppCtx%FractureEnergy(AppCtx%TimeStep), AppCtx%DelaminationEnergy(AppCtx%TimeStep),  AppCtx%FilmEnergy(AppCtx%TimeStep), AppCtx%BondingLayerEnergy(AppCtx%TimeStep), AppCtx%TotalEnergy(AppCtx%TimeStep)
         If (AppCtx%VarFracSchemeParam%SaveBlk) Then
            Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks_Global
               Write(AppCtx%AppParam%EnerBlock_Unit(iBlk), 100) AppCtx%TimeStep, AppCtx%Load(AppCtx%TimeStep), AppCtx%ExtForcesWorkBlock(iBlk),  AppCtx%ElasticEnergyBlock(iBlk), AppCtx%FilmEnergyBlock(iBlk), AppCtx%BondingLayerEnergyBlock(iBlk), AppCtx%FractureEnergyBlock(iBlk),  AppCtx%DelaminationEnergyBlock(iBlk),   AppCtx%TotalEnergyBlock(iBlk)
            End Do    
      End If
      End If
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
100   Format(I6, 8(ES13.5,'  '))  
   End Subroutine Save_Ener
   
#undef __FUNC__ 
#define __FUNC__ "Load_Ener"
   Subroutine Load_Ener(AppCtx, TimeStep)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: TimeStep, iErr

      Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)

      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine Load_Ener
   
#undef __FUNC__ 
#define __FUNC__ "Init_TS_Loads"
   Subroutine Init_TS_Loads(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr

      Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Temperature)%Offset, AppCtx%TimeStep, AppCtx%Theta) 
      
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Fracture)%Offset, AppCtx%TimeStep, AppCtx%VBC) 
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Delamination)%Offset, AppCtx%TimeStep, AppCtx%WBC) 
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_DisplacementX)%Offset, AppCtx%TimeStep, AppCtx%UBC) 

      Call Read_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_Load)%Offset, AppCtx%TimeStep, AppCtx%Load(AppCtx%TimeStep))
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine Init_TS_Loads   
   
#undef __FUNC__ 
#define __FUNC__ "VarFracQSFinalize"
   Subroutine VarFracQSFinalize(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr, iBlk
      Character(len=MEF90_MXSTRLEN)                :: filename
      Type(PetscViewer)                            :: LogViewer
   
      Call PetscLogStagePush(AppCtx%LogInfo%Setup_Stage, iErr); CHKERRQ(iErr)
      Call FieldDestroy(AppCtx%U);      CHKERRQ(iErr)
      Call FieldDestroy(AppCtx%UBC);    CHKERRQ(iErr)
	If (AppCtx%VarFracSchemeParam%U_UseSNES) Then
		Call FieldDestroy(AppCtx%GradientU); CHKERRQ(iErr)
	End If
      Call FieldDestroy(AppCtx%V);      CHKERRQ(iErr)
      Call FieldDestroy(AppCtx%VBC);    CHKERRQ(iErr)
   Call FieldDestroy(AppCtx%W);      CHKERRQ(iErr)
   Call FieldDestroy(AppCtx%WBC);    CHKERRQ(iErr)
      Call FieldDestroy(AppCtx%VIrrev); CHKERRQ(iErr)

      Call FieldDestroy(AppCtx%Theta);  CHKERRQ(iErr)
   Call FieldDestroy(AppCtx%RHSU); CHKERRQ(iErr)
      If (AppCtx%VarFracSchemeParam%V_UseTao) Then
         Call FieldDestroy(AppCtx%GradientV);   CHKERRQ(iErr)
         Call FieldDestroy(AppCtx%LowerBoundV); CHKERRQ(iErr)
         Call FieldDestroy(AppCtx%UpperBoundV); CHKERRQ(iErr)
      Else
         Call FieldDestroy(AppCtx%RHSV); CHKERRQ(iErr)
      End If
      Call VecDestroy(AppCtx%V_Old, iErr); CHKERRQ(iErr)
	Call VecDestroy(AppCtx%U_Old, iErr); CHKERRQ(iErr)
	Call VecDestroy(AppCtx%W_Old, iErr); CHKERRQ(iErr)

      If ( (AppCtx%VarFracSchemeParam%SaveStress) .OR. (AppCtx%VarFracSchemeParam%SaveStrain) ) Then
         Call SectionRealDestroy(AppCtx%StrainU, iErr); CHKERRQ(iErr)
         Call SectionRealDestroy(AppCtx%StressU, iErr); CHKERRQ(iErr)
      End If
      
      Call FlagDestroy(AppCtx%BCUFlag);   CHKERRQ(iErr)
      Call FlagDestroy(AppCtx%BCVFlag);   CHKERRQ(iErr)
      Call FlagDestroy(AppCtx%IrrevFlag); CHKERRQ(iErr)

      Call MatDestroy(AppCtx%KU, iErr); CHKERRQ(iErr)
      Call MatDestroy(AppCtx%KV, iErr); CHKERRQ(iErr)
      Call DMDestroy(AppCtx%MeshTopology%Mesh, iErr); CHKERRQ(ierr)

      DeAllocate(AppCtx%FractureEnergy)
      DeAllocate(AppCtx%ElasticEnergy)
      DeAllocate(AppCtx%BondingLayerEnergy)
      DeAllocate(AppCtx%FilmEnergy)
      DeAllocate(AppCtx%ExtForcesWork)
      DeAllocate(AppCtx%TotalEnergy)
      DeAllocate(AppCtx%Load)
      DeAllocate(AppCtx%FractureEnergyBlock)
      DeAllocate(AppCtx%ElasticEnergyBlock)
      DeAllocate(AppCtx%BondingLayerEnergyBlock)
      DeAllocate(AppCtx%FilmEnergyBlock)
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
         If (AppCtx%VarFracSchemeParam%SaveBlk) Then
            Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks_Global
               Close(AppCtx%AppParam%EnerBlock_Unit(iBlk))
            End Do
         DeAllocate(AppCtx%AppParam%EnerBlock_Unit)
      End If
      End If
      
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
      Write(filename, 103) Trim(AppCtx%AppParam%prefix)
      Call PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, LogViewer, ierr);CHKERRQ(ierr)
      Call PetscViewerASCIIPrintf(LogViewer,filename,ierr);CHKERRQ(ierr)
      Call PetscLogView(LogViewer,ierr);CHKERRQ(ierr)
      Call PetscViewerDestroy(LogViewer,ierr);CHKERRQ(ierr)
      
      Call KSPDestroy(AppCtx%KSPU, iErr); CHKERRQ(iErr)

      If (AppCtx%VarFracSchemeParam%V_UseTao) Then
#if defined WITH_TAO
         Call TaoDestroy(AppCtx%taoV, iErr); CHKERRQ(iErr)
         Call TaoApplicationDestroy(AppCtx%taoAppV, iErr); CHKERRQ(iErr)
#endif
      Else
         Call KSPDestroy(AppCtx%KSPV, iErr); CHKERRQ(iErr)
      End If

      Call EXCLOS(AppCtx%EXO%exoid, iErr)
      AppCtx%EXO%exoid = 0
      Call EXCLOS(AppCtx%MyEXO%exoid, iErr)
      AppCtx%MyEXO%exoid = 0
#if defined WITH_TAO
      Call TaoFinalize(iErr)
#endif
      Call MEF90_Finalize()
103 Format(A,'-logsummary.txt')
   End Subroutine VarFracQSFinalize

#undef __FUNC__ 
#define __FUNC__ "Step_UW"
   Subroutine Step_UW(AppCtx)
	Type(AppCtx_Type)                            :: AppCtx

	PetscInt                                     :: iErr, UWiter
	PetscReal                                    :: ErrW
	Character(len=MEF90_MXSTRLEN)                :: IOBuffer   
	
	UWiter=1
	
	Call Step_U(AppCtx)
	Call Step_W(AppCtx)

	Call VecCopy(AppCtx%W%Vec, AppCtx%W_Old, iErr); CHKERRQ(iErr)
	
! 	Call VecSet(AppCtx%U_Old, 0.0_Kr, iErr); CHKERRQ(iErr)
	
	Do While (ErrW == 1  .OR. UWiter == 1)

		If (AppCtx%AppParam%verbose > 0) Then
			Write(IOBuffer, *) '      UWiter: ', UWiter, '\n' 
			Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
		End If

		Call Step_U(AppCtx)
		Call Step_W(AppCtx)

		Call VecAxPy(AppCtx%W_Old, -1.0_Kr, AppCtx%W%Vec, iErr) ! W_old <- -W + W_old
! 		If (AppCtx%AppParam%verbose > 0) Then
! 			Call PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INDEX , iErr)
! 			Call VecView(AppCtx%W_Old, PETSC_VIEWER_STDOUT_WORLD, iErr)
! 		End If
		Call VecNorm(AppCtx%W_Old, NORM_INFINITY, ErrW, iErr)

		Write(IOBuffer, 800) ErrW
		Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
		Call VecCopy(AppCtx%W%Vec, AppCtx%W_Old, iErr); CHKERRQ(iErr)

		Call Save_U(AppCtx)
		Call Save_W(AppCtx)

		UWiter = UWiter + 1
	End Do

	800 Format('     Max change W: ', T24, ES12.5, '\n')

End Subroutine Step_UW

End Module m_VarFilmQS

