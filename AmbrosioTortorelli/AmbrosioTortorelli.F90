Program  Elast

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscmeshdef.h"

#if defined PB_2D
   Use m_AmbrosioTortorelli2D
#elif defined PB_3D
   Use m_AmbrosioTortorelli3D
#endif   
   Use m_MEF90
   Use m_RuptStruct
   Use petsc
   Use petscvec
   Use petscmat
   Use petscksp
   Use petscmesh
   
   Implicit NONE   


   Type(AppCtx_Type)                            :: AppCtx
   PetscInt                                     :: iErr
   Character(len=MEF90_MXSTRLEN)                :: IOBuffer
   PetscInt                                     :: iter

   Call AmbrosioTortorelliInit(AppCtx)
   
   If (AppCtx%AppParam%verbose) Then
      Call EXOView(AppCtx%EXO, AppCtx%AppParam%LogViewer)
      Call EXOView(AppCtx%MyEXO, AppCtx%AppParam%MyLogViewer) 
      Call MeshTopologyView(AppCtx%MeshTopology, AppCtx%AppParam%MyLogViewer) 
   End If   
   
   AppCtx%TimeStep = 1
   Do iter=1, AppCtx%RuptSchemeParam%RelaxMaxIter
      Write(IOBuffer, "('Iteration ', I4,A)") iter, '\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

      !------------------------------------------------------------------- 
      ! Problem for U
      !-------------------------------------------------------------------
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) 'Assembling the Matrix and RHS for  the U-subproblem \n'c 
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      Call RHSU_Assembly(AppCtx)
   !   Call VecView(AppCtx%RHSU, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
      Call MatU_Assembly(AppCtx)
   !   Call MatView(AppCtx%KU, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
      
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) 'Calling KSPSolve for the U-subproblem\n'c 
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
      End If
      Call Solve_U(AppCtx)
      
      !------------------------------------------------------------------- 
      ! Problem for V
      !-------------------------------------------------------------------
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) 'Assembling the Matrix and RHS for the V-subproblem \n'c 
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
      End If
      Call RHSV_Assembly(AppCtx)
   !   Call VecView(AppCtx%RHSV, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
      Call MatV_Assembly(AppCtx)
   !   Call MatView(AppCtx%KV, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
      
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) 'Calling KSPSolve for the V-subproblem\n'c 
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
      End If
      Call Solve_V(AppCtx)

      !------------------------------------------------------------------- 
      ! Check the exit condition: tolerance on the error in V 
      !------------------------------------------------------------------- 
      If (AppCtx.ErrV.LT.AppCtx%RuptSchemeParam%RelaxTol) then 
         EXIT 
      End If
      If (Mod(iter, AppCtx%RuptSchemeParam%SaveInt) == 0) Then
         If (AppCtx%AppParam%verbose) Then
            Write(IOBuffer, *) 'Saving U and V\n'c
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
         End If   
         Call Save_U(AppCtx)
         Call Save_V(AppCtx)
         Call ComputeEnergy(AppCtx)
      
         Write(IOBuffer, 100) AppCtx%ElasticEnergy
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Write(IOBuffer, 101) AppCtx%ExtForcesWork
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Write(IOBuffer, 102) AppCtx%SurfaceEnergy
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Write(IOBuffer, 103) AppCtx%TotalEnergy
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

         Call ComputeStrainStress(AppCtx)
         Call Save_StrainStress(AppCtx)
      End If
   End Do
   
   !------------------------------------------------------------------- ! Save the results
   Call Save_U(AppCtx)
   Call Save_V(AppCtx)
   !-------------------------------------------------------------------
   
   If (AppCtx%AppParam%verbose) Then
      Write(IOBuffer, *) 'Computing bulk energy, strains and stresses\n'c 
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
   End If
   
   Call ComputeEnergy(AppCtx)
   Call Save_Ener(AppCtx)

   Write(IOBuffer, 100) AppCtx%ElasticEnergy
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Write(IOBuffer, 101) AppCtx%ExtForcesWork
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Write(IOBuffer, 102) AppCtx%SurfaceEnergy
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Write(IOBuffer, 103) AppCtx%TotalEnergy
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

100 Format('Elastic energy:       ', ES12.5, '\n'c)    
101 Format('External Forces Work: ', ES12.5, '\n'c)    
102 Format('Surface energy:       ', ES12.5, '\n'c)    
103 Format('Total energy:         ', ES12.5, '\n'c)    

   Call ComputeStrainStress(AppCtx)

   Call AmbrosioTortorelliFinalize(AppCtx)
End Program  Elast
