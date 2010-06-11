Program  Elast
#include "finclude/petscdef.h"

   Use m_MEF90
#if defined PB_2D
   Use m_Elast2D
#elif defined PB_3D 
   Use m_Elast3D
#endif

   Implicit NONE   


   Type(AppCtx_Type)                            :: AppCtx
   PetscInt                                     :: i, iErr
   Character(len=MEF90_MXSTRLEN)                :: IOBuffer
   PetscInt                                     :: vers
   
   PetscReal                                    :: VolChange

   Call ElastInit(AppCtx)
   
   If (AppCtx%AppParam%verbose > 1) Then
      Call EXOView(AppCtx%EXO, AppCtx%AppParam%LogViewer)
      Call EXOView(AppCtx%MyEXO, AppCtx%AppParam%MyLogViewer)
      Call MeshTopologyView(AppCtx%MeshTopology, AppCtx%AppParam%MyLogViewer)
   End If   

   If (.NOT. AppCtx%VarFracSchemeParam%U_UseTao) Then
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Assembling the matrix\n'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If   
      Call MatAssembly(AppCtx)
      If (AppCtx%AppParam%verbose > 1) Then
         Call MatView(AppCtx%KU, AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
      End If
   End If
   
   Do i = 1, AppCtx%NumTimeSteps
      AppCtx%TimeStep = i
      
      AppCtx%MyEXO%exoid = EXOPEN(AppCtx%MyEXO%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers, iErr)
      Call EXGTIM(AppCtx%MyEXO%exoid, AppCtx%TimeStep, AppCtx%Time, iErr)
      Call EXCLOS(AppCtx%MyEXO%exoid, iErr)
      AppCtx%MyEXO%exoid = 0
      Call Read_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_Load)%Offset, AppCtx%TimeStep, AppCtx%Load)

      Call Solve(AppCtx)
   
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Computing Elastic energy, strains and stresses\n'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      Call ComputeEnergies(AppCtx)
      Write(IOBuffer, 108) AppCtx%TimeStep, AppCtx%Time, AppCtx%Load, AppCtx%ElasticEnergy, AppCtx%ExtForcesWork, AppCtx%ElasticEnergy - AppCtx%ExtForcesWork
108 Format('TS ',I4, ' Time:', ES10.3, ' Load:', ES10.3, ' Elast:', ES10.3, ' Work:', ES10.3, ' Total:', ES10.3, '\n')
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, 110) AppCtx%TimeStep, AppCtx%Load, AppCtx%ElasticEnergy, AppCtx%ExtForcesWork, 0.0_Kr, AppCtx%ElasticEnergy - AppCtx%ExtForcesWork
110 Format(I4, 5(ES13.5, '   '), '\n')
      Call PetscViewerASCIIPrintf(AppCtx%AppParam%EnergyViewer, IOBuffer, iErr); CHKERRQ(iErr)

      Call ComputeVolumeChange(VolChange, AppCtx)
      Write(IOBuffer, 120) VolChange / AppCtx%Load
120 Format('Volume change: ', ES12.5, ' \n')
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

      If ( (AppCtx%VarFracSchemeParam%SaveStress) .OR. ( AppCtx%VarFracSchemeParam%SaveStrain) ) Then
         Call ComputeStrainStress(AppCtx)
      End If   
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Saving results\n'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If      
      Call Save(AppCtx)
   End Do

   Call ElastFinalize(AppCtx)
End Program  Elast
