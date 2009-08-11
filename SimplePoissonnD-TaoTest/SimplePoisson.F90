Program  SimplePoisson

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscmeshdef.h"

   Use m_MEF90
#if defined PB_2D
   Use m_SimplePoisson2D
!   Use m_SimplePoissonTao2D
#elif defined PB_3D 
   Use m_SimplePoisson3D
!   Use m_SimplePoissonTao3D
#endif

   Use petsc
   Use petscvec
   Use petscmat
   Use petscksp
   Use petscmesh

   Implicit NONE   


   Type(AppCtx_Type)                            :: AppCtx
   PetscInt                                     :: iErr
   Character(len=MEF90_MXSTRLEN)                :: IOBuffer
   PetscInt                                     :: Step
   PetscReal                                    :: t

   Call SimplePoissonInit(AppCtx)
   
   If (AppCtx%AppParam%verbose > 1) Then
      Call EXOView(AppCtx%EXO, AppCtx%AppParam%LogViewer)
      Call EXOView(AppCtx%MyEXO, AppCtx%AppParam%MyLogViewer)
      Call MeshTopologyView(AppCtx%MeshTopology, AppCtx%AppParam%MyLogViewer)
   End If   

   If (.NOT. AppCtx%AppParam%Use_Tao) Then
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Assembling the matrix\n'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      Call MatAssembly(AppCtx)
      If (AppCtx%AppParam%verbose > 1) Then
         Write(IOBuffer, *) 'Matrix\n'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Call MatView(AppCtx%K, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
      End If
   End If
   
   
   Do step = 0, AppCtx%NumSteps-1
      If (AppCtx%NumSteps == 1) Then
         t = AppCtx%tmin
      Else
         t = AppCtx%tmin + Real(step) * (AppCtx%Tmax - AppCtx%Tmin) / (AppCtx%NumSteps - 1.0)
      End If
      Write(*,*) 'Step ', step, ' load: ', t

      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Calling InitLoads\n'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      Call InitLoads(AppCtx, t, iErr)

      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Calling Solve\n'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      Call Solve(AppCtx)
   
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Computing Energies\n'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      Call ComputeEnergy(AppCtx)
   
      Write(IOBuffer, 100) AppCtx%ElasticEnergy, AppCtx%ExtForcesWork, AppCtx%TotalEnergy
100   Format('Elastic energy: ', ES12.5, ' Forces Work: ', ES12.5, ' Total: ', ES12.5, '\n')    
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   
      Call ComputeGradU(AppCtx)
   
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Saving results\n'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
   
      Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
!      Call Write_EXO_Result_Global(AppCtx%MyExo, 1, step+1, AppCtx%ElasticEnergy)
!      Call Write_EXO_Result_Global(AppCtx%MyExo, 2, step+1, AppCtx%ExtForcesWork)
!      Call Write_EXO_Result_Global(AppCtx%MyExo, 3, step+1, AppCtx%TotalEnergy)
      Call Write_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, 1, step+1, AppCtx%U) 
      Call Write_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, 2, step+1, AppCtx%F) 
      Call Write_EXO_Result_Cell(AppCtx%MyEXO, AppCtx%MeshTopology, 1, step+1, AppCtx%GradU) 
      Call PetscLogStagePop (AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
   End Do
      
   Call SimplePoissonFinalize(AppCtx)
End Program  SimplePoisson
