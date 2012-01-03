Program  TransientHeat

#include "finclude/petscdef.h"

   Use m_MEF90
#if defined PB_2D
   Use m_TransientHeat2D
#elif defined PB_3D 
   Use m_TransientHeat3D
#endif

   Implicit NONE   


   Type(Heat_AppCtx_Type)                       :: AppCtx
   PetscInt                                     :: iErr, iStep
   Character(len=MEF90_MXSTRLEN)                :: IOBuffer
   Type(Vec)                                    :: W1, W2, W3 
   PetscScalar                                  :: prodMassUnit
   PetscReal                                    :: Mass
   PetscReal, Dimension(:), Pointer             :: CurTime
   Type(Field)                                  :: ExtraField


   Call PoissonInit(AppCtx)

   Call Poisson_TSSetUp(AppCtx, AppCtx%MeshTopology)


   If (AppCtx%AppParam%verbose > 4) Then
      Call EXOView(AppCtx%EXO, AppCtx%AppParam%LogViewer)
      Call EXOView(AppCtx%MyEXO, AppCtx%AppParam%MyLogViewer)
      Call MeshTopologyView(AppCtx%MeshTopology, AppCtx%AppParam%MyLogViewer)
   End If   

   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Assembling the matrix\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Assembling the RHS\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   Call RHSAssembly(AppCtx, AppCtx%MeshTopology, AppCtx%MyExo)
   If (AppCtx%AppParam%verbose > 2) Then
      Write(IOBuffer, *) 'RHS\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call SectionRealView(AppCtx%RHS, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   End If
  

   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Assembling the Mass - Variational  Identity   matrix\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   Call MatMassAssembly(AppCtx, AppCtx%MeshTopology)

   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Calling Solve\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   Allocate(CurTime(AppCtx%NumSteps))
   !CurTime(1) = 0
   TimeStep :DO iStep = 1, AppCtx%NumSteps-1
       !! Non uniform time stepping adapted to the time scale of the  thermal problem in tau=sqrt(t)
      Call MatZeroEntries(AppCtx%K, iErr); CHKERRQ(iErr)
      Call HeatMatAssembly(AppCtx, AppCtx%MeshTopology, ExtraField)
      If (AppCtx%AppParam%verbose > 3) Then
         Write(IOBuffer, *) 'Matrix\n'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Call MatView(AppCtx%K, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
      End If
!Next 2 lines how it should be : 
!      Call Read_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(Heat_GlobVar_Load)%Offset, iStep, CurTime(iStep))
!      Call Read_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(Heat_GlobVar_Load)%Offset, iStep+1, CurTime(iStep+1))
      Call Read_EXO_Result_Global(AppCtx%MyEXO, Heat_GlobVar_Load, iStep, CurTime(iStep))
      Call Read_EXO_Result_Global(AppCtx%MyEXO, Heat_GlobVar_Load, iStep+1, CurTime(iStep+1))
      Call SolveTransientStep(AppCtx, AppCtx%MyEXO, AppCtx%MeshTopology, CurTime(iStep), CurTime(iStep+1), iStep)
      Call Write_EXO_Result_Global(AppCtx%MyEXO, Heat_GlobVar_Load, iStep+1, CurTime(iStep+1))
      Call EXPTIM(AppCtx%MyEXO%exoid, iStep+1, CurTime(iStep+1), iErr)
!Save the end of the time step in case it has been changed 
   End Do TimeStep
   DeAllocate(CurTime) 

   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Computing Energies\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   Call ComputeEnergy(AppCtx)

   Write(IOBuffer, 100) AppCtx%ElasticEnergy, AppCtx%ExtForcesWork, AppCtx%TotalEnergy
100 Format('Elastic energy: ', ES12.5, ' Forces Work: ', ES12.5, ' Total: ', ES12.5, '\n')    
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Saving results\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If

!! Computing water loss
   Call ComputeWaterMass(AppCtx%MeshTopology, AppCtx%MyExo, 1, AppCtx%Elem, AppCtx%NumSteps)

   Call Write_EXO_Result_Global(AppCtx%MyExo, Heat_GlobVar_ElasticEnergy,  1, AppCtx%ElasticEnergy)
   Call Write_EXO_Result_Global(AppCtx%MyExo, Heat_GlobVar_ExtForcesWork,  1, AppCtx%ExtForcesWork)
   Call Write_EXO_Result_Global(AppCtx%MyExo, Heat_GlobVar_TotalEnergy,    1, AppCtx%TotalEnergy)
   Call Write_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, 2, 1, AppCtx%F%Sec) 

   Call EXCLOS(AppCtx%EXO%exoid, iErr)
   AppCtx%EXO%exoid = 0
   Call EXCLOS(AppCtx%MyEXO%exoid, iErr)
   AppCtx%MyEXO%exoid = 0

   Call TSPoissonFinalize(AppCtx)
End Program  TransientHeat
