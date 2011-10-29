Program  TransientHeat

#include "finclude/petscdef.h"

   Use m_MEF90
#if defined PB_2D
   Use m_Poisson2D
   Use m_TransientHeat2D
#elif defined PB_3D 
   Use m_Poisson3D
   Use m_TransientHeat3D
#endif

   Implicit NONE   


   Type(Heat_AppCtx_Type)                       :: AppCtx
   PetscInt                                     :: iErr
   Character(len=MEF90_MXSTRLEN)                :: IOBuffer
   Type(Vec)                                    :: W1, W2, W3 
   PetscScalar                                  :: prodMassUnit
   PetscReal                                    :: Mass 
   Call PoissonInit(AppCtx)

   Select Case (AppCtx%AppParam%TestCase)
   Case (1)
      Call KSPSetUp(AppCtx)
   Case(2)
      Call Poisson_TSSetUp(AppCtx, AppCtx%MeshTopology)
   End Select


   If (AppCtx%AppParam%verbose > 4) Then
      Call EXOView(AppCtx%EXO, AppCtx%AppParam%LogViewer)
      Call EXOView(AppCtx%MyEXO, AppCtx%AppParam%MyLogViewer)
      Call MeshTopologyView(AppCtx%MeshTopology, AppCtx%AppParam%MyLogViewer)
   End If   

   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Assembling the matrix\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   Call HeatMatAssembly(AppCtx, AppCtx%MeshTopology)
   If (AppCtx%AppParam%verbose > 3) Then
      Write(IOBuffer, *) 'Matrix\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call MatView(AppCtx%K, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   End If

   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Assembling the RHS\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   Call RHSAssembly(AppCtx, AppCtx%MeshTopology)
   If (AppCtx%AppParam%verbose > 2) Then
      Write(IOBuffer, *) 'RHS\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call SectionRealView(AppCtx%RHS, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   End If
   
   
   Select Case (AppCtx%AppParam%TestCase)
   Case (1) 
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Calling Solve\n'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      Call Solve(AppCtx)
   Case(2)
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Assembling the Mass - Variational  Identity   matrix\n'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      Call MatMassAssembly(AppCtx, AppCtx%MeshTopology)

      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Calling Solve\n'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      Call SolveTransient(AppCtx, AppCtx%MyEXO, AppCtx%MeshTopology)
   End Select 

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

!! Computing water loss
!   Call ComputeWaterLoss(AppCtx%MeshTopology, AppCtx%Exo, 'U')
!   Call IntegrateScal(AppCtx%MeshTopology,1,AppCtx%Elem,AppCtx%U%Sec,Mass)

   Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
   Call Write_EXO_Result_Global(AppCtx%MyExo, 1, 1, AppCtx%ElasticEnergy)
   Call Write_EXO_Result_Global(AppCtx%MyExo, 2, 1, AppCtx%ExtForcesWork)
   Call Write_EXO_Result_Global(AppCtx%MyExo, 3, 1, AppCtx%TotalEnergy)
!   Call Write_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, 1, 1, AppCtx%U%Sec) 
   Call Write_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, 2, 1, AppCtx%F%Sec) 
   Call Write_EXO_Result_Cell(AppCtx%MyEXO, AppCtx%MeshTopology, 1, 1, AppCtx%GradU) 
   Call PetscLogStagePop (AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)

   Call EXCLOS(AppCtx%EXO%exoid, iErr)
   AppCtx%EXO%exoid = 0
   Call EXCLOS(AppCtx%MyEXO%exoid, iErr)
   AppCtx%MyEXO%exoid = 0

   Select Case (AppCtx%AppParam%TestCase)
   Case (1)
      Call SimplePoissonFinalize(AppCtx)
   Case (2)
      Call TSPoissonFinalize(AppCtx)
   End Select 
End Program  TransientHeat
