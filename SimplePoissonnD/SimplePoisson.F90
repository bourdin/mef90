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
#elif defined PB_3D 
   Use m_SimplePoisson3D
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
   Type(Vec)                                    :: V

   Call SimplePoissonInit(AppCtx)
   
   If (AppCtx%AppParam%verbose > 1) Then
      Call EXOView(AppCtx%EXO, AppCtx%AppParam%LogViewer)
      Call EXOView(AppCtx%MyEXO, AppCtx%AppParam%MyLogViewer)
      Call MeshTopologyView(AppCtx%MeshTopology, AppCtx%AppParam%MyLogViewer)
   End If   

   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Assembling the matrix\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   Call MatAssembly_Elementwise(AppCtx)
   If (AppCtx%AppParam%verbose > 1) Then
      Write(IOBuffer, *) 'Matrix\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call MatView(AppCtx%K, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   End If
   
   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Assembling the RHS\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   Call RHSAssembly(AppCtx)
   If (AppCtx%AppParam%verbose > 1) Then
      Write(IOBuffer, *) 'RHS\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call VecView(AppCtx%RHS, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   End If
   
   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Calling KSPSolve\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   Call Solve(AppCtx)
   
   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Computing energy and gradient\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   Call ComputeEnergy(AppCtx)

   Write(IOBuffer, 100) AppCtx%Energy
100 Format('Total energy: ', ES12.5, '\n'c)    
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

   Call ComputeGradient(AppCtx)

   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Saving results\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If

   Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
   Call Write_EXO_Result_Global(AppCtx%MyExo, 1, 1, AppCtx%Energy)
   Call Write_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, 1, 1, AppCtx%U) 
   Call Write_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, 2, 1, AppCtx%F) 
   Call Write_EXO_Result_Cell(AppCtx%MyEXO, AppCtx%MeshTopology, 1, 1, AppCtx%GradU) 
   Call PetscLogStagePop (AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
   
   Call SimplePoissonFinalize(AppCtx)
End Program  SimplePoisson
