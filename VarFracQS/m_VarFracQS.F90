#if defined PB_2D
Module m_VarFracQS2D
#elif defined PB_3D
Module m_VarFracQS3D
#endif

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscmeshdef.h"
#include "finclude/petscviewerdef.h"

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
   Use petsc
   Use petscvec
   Use petscmat
   Use petscksp
   Use petscmesh

   Implicit NONE   
   
Contains

   Subroutine VarFracQSInit(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx

      PetscInt                                     :: iErr, i
      Type(Mesh)                                   :: Tmp_Mesh
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer, filename
      PetscInt                                     :: NumComponents
      PetscTruth                                   :: HasPrefix
      PetscReal                                    :: rDummy
      Character                                    :: cDummy
      PetscInt                                     :: vers

      
      Call MEF90_Initialize()
      Call PetscOptionsHasName(PETSC_NULL_CHARACTER, '-verbose', AppCtx%AppParam%Verbose, iErr)    
      Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p',     AppCtx%AppParam%prefix, HasPrefix, iErr); CHKERRQ(iErr)
      If (.NOT. HasPrefix) Then
         Call PetscPrintf(PETSC_COMM_WORLD, "No input file prefix given\n", iErr)
         Call MEF90_Finalize()
         STOP
      End If
      Call VarFracSchemeParam_GetFromOptions(AppCtx%VarFracSchemeParam)
      
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

      Call Write_EXO_Case(AppCtx%AppParam%prefix, '%0.4d', MEF90_NumProcs)
      AppCtx%EXO%Comm = PETSC_COMM_WORLD
      AppCtx%EXO%filename = Trim(AppCtx%AppParam%prefix)//'.gen'

         !!! Reading and distributing sequential mesh
      If (MEF90_NumProcs == 1) Then
         Call MeshCreateExodus(PETSC_COMM_WORLD, AppCtx%EXO%filename, AppCtx%MeshTopology%mesh, ierr); CHKERRQ(iErr)
      Else
         Call MeshCreateExodus(PETSC_COMM_WORLD, AppCtx%EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
         Call MeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, AppCtx%MeshTopology%mesh, ierr); CHKERRQ(iErr)
         Call MeshDestroy(Tmp_mesh, ierr); CHKERRQ(iErr)
      End If
   
      Call MeshTopologyReadEXO(AppCtx%MeshTopology, AppCtx%EXO)
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) "Done reading and partitioning the mesh\n"c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
   
      AppCtx%MyEXO%comm = PETSC_COMM_SELF
      AppCtx%MyEXO%exoid = AppCtx%EXO%exoid
      Write(AppCtx%MyEXO%filename, 99) trim(AppCtx%AppParam%prefix), MEF90_MyRank
 99  Format(A, '-', I4.4, '.gen')
   
      !!! Initializes the values and names of the properties and variables
      Call VarFracEXOVariable_Init(AppCtx%MyEXO)
      Call EXOProperty_Read(AppCtx%MyEXO)   
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) "Done with VarFracQSEXOVariable_Init and VarFracQSEXOProperty_Read\n"c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Call MeshTopologyView(AppCtx%MeshTopology, AppCtx%AppParam%MyLogViewer)
         Call EXOView(AppCtx%MyEXO, AppCtx%AppParam%MyLogViewer)
      End If
      
      !!! Read Mat Properties from the CST file
      Call MatProp_Read(AppCtx%MeshTopology, AppCtx%MatProp, trim(AppCtx%AppParam%prefix)//'.CST')
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) "Done with MatProp_Read\n"c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
!         Call EXOView(AppCtx%MyEXO, AppCtx%AppParam%MyLogViewer)
      End If

      !!! Set the element type for each block so that we can call ElementInit
      Do i = 1, AppCtx%MeshTopology%num_elem_blks
         AppCtx%MeshTopology%elem_blk(i)%Elem_Type = AppCtx%MyEXO%EBProperty( VarFrac_EBProp_Elem_Type )%Value( AppCtx%MeshTopology%elem_blk(i)%ID )
         Call Init_Elem_Blk_Type(AppCtx%MeshTopology%Elem_Blk(i), AppCtx%MeshTopology%num_dim)
      End Do
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) "Done with Init_Elem_Blk_Type\n"c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      
      
      Call ElementInit(AppCtx%MeshTopology, AppCtx%ElemVect, AppCtx%VarFracSchemeParam%IntegOrder)
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) "Done with ElementInit Vect\n"c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      Call ElementInit(AppCtx%MeshTopology, AppCtx%ElemScal, AppCtx%VarFracSchemeParam%IntegOrder)
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) "Done with ElementInit Scal\n"c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If

      !!! Create the Sections for the variables
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'U', AppCtx%MeshTopology%Num_Dim, AppCtx%U, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'F', AppCtx%MeshTopology%Num_Dim, AppCtx%F, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'Theta', 1, AppCtx%Theta, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'V', 1, AppCtx%V, iErr); CHKERRQ(iErr)

      NumComponents = AppCtx%MeshTopology%Num_Dim * (AppCtx%MeshTopology%Num_Dim + 1) / 2
      Call MeshGetCellSectionReal(AppCtx%MeshTopology%mesh, 'Strain', NumComponents, AppCtx%StrainU, iErr); CHKERRQ(iErr)
      Call MeshGetCellSectionReal(AppCtx%MeshTopology%mesh, 'Stress', NumComponents, AppCtx%StressU, iErr); CHKERRQ(iErr)
      Call MeshGetCellSectionReal(AppCtx%MeshTopology%mesh, 'GradV', AppCtx%MeshTopology%Num_Dim, AppCtx%GradV, iErr); CHKERRQ(iErr)


      !!! Create the Scatter, Vec, Mat, KSP, PC
      Call MeshCreateGlobalScatter(AppCtx%MeshTopology%mesh, AppCtx%U, AppCtx%ScatterVect, iErr); CHKERRQ(iErr)
      Call MeshCreateGlobalScatter(AppCtx%MeshTopology%mesh, AppCtx%V, AppCtx%ScatterScal, iErr); CHKERRQ(iErr)

      Call MeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%U, AppCtx%RHSU, iErr); CHKERRQ(iErr)
      Call MeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%V, AppCtx%RHSV, iErr); CHKERRQ(iErr)

      Call MeshSetMaxDof(AppCtx%MeshTopology%Mesh, AppCtx%MeshTopology%Num_Dim, iErr); CHKERRQ(iErr) 
      Call MeshCreateMatrix(AppCtx%MeshTopology%mesh, AppCtx%U, MATMPIAIJ, AppCtx%KU, iErr); CHKERRQ(iErr)
      Call MeshCreateMatrix(AppCtx%MeshTopology%mesh, AppCtx%V, MATMPIAIJ, AppCtx%KV, iErr); CHKERRQ(iErr)
      
      !! Solver context for U      
      Call KSPCreate(PETSC_COMM_WORLD, AppCtx%KSPU, iErr); CHKERRQ(iErr)
      Call KSPSetOperators(AppCtx%KSPU, AppCtx%KU, AppCtx%KU, SAME_NONZERO_PATTERN, iErr); CHKERRQ(iErr)
      Call KSPSetType(AppCtx%KSPU, KSPCG, iErr); CHKERRQ(iErr)
      Call KSPSetInitialGuessNonzero(AppCtx%KSPU, PETSC_TRUE, iErr); CHKERRQ(iErr)
      Call KSPSetFromOptions(AppCtx%KSPU, iErr); CHKERRQ(iErr)

      Call KSPGetPC(AppCtx%KSPU, AppCtx%PCU, iErr); CHKERRQ(iErr)
      Call PCSetType(AppCtx%PCU, PCBJACOBI, iErr); CHKERRQ(iErr)
   
      !! Solver context for V      
      Call KSPCreate(PETSC_COMM_WORLD, AppCtx%KSPV, iErr); CHKERRQ(iErr)
      Call KSPSetOperators(AppCtx%KSPV, AppCtx%KV, AppCtx%KV, SAME_NONZERO_PATTERN, iErr); CHKERRQ(iErr)
      Call KSPSetType(AppCtx%KSPV, KSPCG, iErr); CHKERRQ(iErr)
      Call KSPSetInitialGuessNonzero(AppCtx%KSPV, PETSC_TRUE, iErr); CHKERRQ(iErr)
      Call KSPSetFromOptions(AppCtx%KSPV, iErr); CHKERRQ(iErr)

      Call KSPGetPC(AppCtx%KSPV, AppCtx%PCV, iErr); CHKERRQ(iErr)
      Call PCSetType(AppCtx%PCV, PCBJACOBI, iErr); CHKERRQ(iErr)
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) "Done Creating fields Section, Vec, KSP and Mat\n"c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If


      !!! Create the Section for the BC
      Call MeshGetVertexSectionInt(AppCtx%MeshTopology%mesh, 'BCUFlag', AppCtx%MeshTopology%Num_Dim, AppCtx%BCUFlag, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionInt(AppCtx%MeshTopology%mesh, 'BCVFlag', 1, AppCtx%BCVFlag, iErr); CHKERRQ(iErr)
#if defined PB_2D
      Call EXOProperty_InitBCUFlag2D(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%BCUFlag)
#elif defined PB_3D
      Call EXOProperty_InitBCUFlag3D(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%BCUFlag)
#endif
      Call EXOProperty_InitBCVFlag(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%BCVFlag)

      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) "Done Initializing BC Sections\n"c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If

      !!! Get the number of time steps
      AppCtx%MyEXO%exoid = EXOPEN(AppCtx%MyEXO%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers, iErr)
      Call EXINQ(AppCtx%MyEXO%exoid, EXTIMS, AppCtx%NumTimeSteps, rDummy, cDummy, iErr)
      Call EXCLOS(AppCtx%MyEXO%exoid, iErr)
      AppCtx%MyEXO%exoid = 0
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) 'Total Number of Time Steps', AppCtx%NumTimeSteps, '\n'c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      
      !!! Set V=1
      Call SectionRealSet(AppCtx%V, 1.0_Kr, iErr); CHKERRQ(iErr)
   End Subroutine VarFracQSInit
         
   
   Subroutine Save_U(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx

      Call Write_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_DisplacementX)%Offset, AppCtx%TimeStep, AppCtx%U) 
   End Subroutine Save_U

   
   Subroutine Save_V(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx

      Call Write_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Fracture)%Offset, AppCtx%TimeStep, AppCtx%V) 
   End Subroutine Save_V

   Subroutine Save_StrainStress(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
   
      Call Write_EXO_Result_Cell(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%CellVariable(VarFrac_CellVar_StrainXX)%Offset, AppCtx%TimeStep, AppCtx%StrainU) 
      Call Write_EXO_Result_Cell(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%CellVariable(VarFrac_CellVar_StressXX)%Offset, AppCtx%TimeStep, AppCtx%StressU) 
   End Subroutine Save_StrainStress


   Subroutine Save_Ener(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx

      Call Write_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_SurfaceEnergy)%Offset, AppCtx%TimeStep, AppCtx%SurfaceEnergy)
      Call Write_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_ElasticEnergy)%Offset, AppCtx%TimeStep, AppCtx%ElasticEnergy)
      Call Write_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_ExtForcesWork)%Offset, AppCtx%TimeStep, AppCtx%ExtForcesWork)
      Call Write_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_TotalEnergy)%Offset, AppCtx%TimeStep, AppCtx%TotalEnergy)
      Call Write_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_Load)%Offset, AppCtx%TimeStep, AppCtx%Load)
   End Subroutine Save_Ener
   
   Subroutine Init_TS_Loads(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr

      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_ForceX)%Offset, AppCtx%TimeStep, AppCtx%F) 
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Temperature)%Offset, AppCtx%TimeStep, AppCtx%Theta) 

      Call Read_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_Load)%Offset, AppCtx%TimeStep, AppCtx%Load)
   End Subroutine Init_TS_Loads   
   
   Subroutine VarFracQSFinalize(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx

      PetscInt                                     :: iErr
      Character(len=MEF90_MXSTRLEN)                :: filename
   
      Call SectionRealDestroy(AppCtx%U, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(AppCtx%V, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(AppCtx%F, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(AppCtx%Theta, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(AppCtx%StrainU, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(AppCtx%StressU, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(AppCtx%GradV, iErr); CHKERRQ(iErr)
      
      Call VecScatterDestroy(AppCtx%ScatterVect, iErr); CHKERRQ(iErr)
      Call VecScatterDestroy(AppCtx%ScatterScal, iErr); CHKERRQ(iErr)
      Call SectionIntDestroy(AppCtx%BCUFlag, iErr); CHKERRQ(iErr)
      Call SectionIntDestroy(AppCtx%BCVFlag, iErr); CHKERRQ(iErr)
      Call MatDestroy(AppCtx%KU, iErr); CHKERRQ(iErr)
      Call VecDestroy(AppCtx%RHSU, iErr); CHKERRQ(iErr)
      Call KSPDestroy(AppCtx%KSPU, iErr); CHKERRQ(iErr)
      Call MatDestroy(AppCtx%KV, iErr); CHKERRQ(iErr)
      Call VecDestroy(AppCtx%RHSV, iErr); CHKERRQ(iErr)
      Call KSPDestroy(AppCtx%KSPV, iErr); CHKERRQ(iErr)
      Call MeshDestroy(AppCtx%MeshTopology%Mesh, iErr); CHKERRQ(ierr)

      If (AppCtx%AppParam%verbose) Then
         Call PetscViewerFlush(AppCtx%AppParam%MyLogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerFlush(AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerDestroy(AppCtx%AppParam%MyLogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerDestroy(AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
      End If
      Write(filename, 103) Trim(AppCtx%AppParam%prefix)
      Call PetscLogPrintSummary(PETSC_COMM_WORLD, filename, iErr); CHKERRQ(iErr)
      Call MEF90_Finalize()
103 Format(A,'.log')
   End Subroutine VarFracQSFinalize
   
   
#if defined PB_2D
End Module m_VarFracQS2D
#elif defined PB_3D
End Module m_VarFracQS3D
#endif
