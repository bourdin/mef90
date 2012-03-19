Program  SimplePoisson

#include "finclude/petscdef.h"

   Use m_MEF90
#if defined PB_2D
   Use m_Poisson2D
#elif defined PB_3D 
   Use m_Poisson3D
#endif

   Implicit NONE   


   Type(Poisson_AppCtx_Type)                    :: AppCtx
   PetscInt                                     :: iErr
   Character(len=MEF90_MXSTRLEN)                :: IOBuffer
   Call SimplePoissonInit(AppCtx)
   
   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Assembling the matrix\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   Call HeatMatAssembly(AppCtx, AppCtx%mesh)
   If (AppCtx%AppParam%verbose > 3) Then
      Write(IOBuffer, *) 'Matrix\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call MatView(AppCtx%K, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   End If

   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Assembling the RHS\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   Call RHSAssembly(AppCtx%mesh,RHS,V,AppCtx)
   If (AppCtx%AppParam%verbose > 2) Then
      Write(IOBuffer, *) 'RHS\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call SectionRealView(AppCtx%RHS, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   End If
   
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
100 Format('Elastic energy: ', ES12.5, ' Forces Work: ', ES12.5, ' Total: ', ES12.5, '\n')    
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

   Call ComputeGradU(AppCtx)

   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Saving results\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If

   Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
   !!! 
   Call VecViewExodusVertex(AppCtx%mesh,AppCtx%U%LocalVec,IO_COMM,AppCtx%EXO%exoid,1,1)
   Call VecViewExodusVertex(AppCtx%mesh,AppCtx%F%LocalVec,IO_COMM,AppCtx%EXO%exoid,1,2)
   Call VecViewExodusCell(AppCtx%mesh,AppCtx%GradU%LocalVec,IO_COMM,AppCtx%EXO%exoid,1,1)
   
   Call Write_EXO_Result_Global(AppCtx%MyExo, 1, 1, AppCtx%ElasticEnergy)
   Call Write_EXO_Result_Global(AppCtx%MyExo, 2, 1, AppCtx%ExtForcesWork)
   Call Write_EXO_Result_Global(AppCtx%MyExo, 3, 1, AppCtx%TotalEnergy)
   Call PetscLogStagePop (AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
  
    Call EXCLOS(AppCtx%EXO%exoid, iErr)
    AppCtx%EXO%exoid = 0

   Call SimplePoissonFinalize(AppCtx)
End Program  SimplePoisson
