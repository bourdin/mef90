Program  VarFracQS

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscmeshdef.h"

#if defined PB_2D
   Use m_VarFracQS2D
#elif defined PB_3D
   Use m_VarFracQS3D
#endif   
   Use m_MEF90
   Use m_VarFrac_Struct
   Use petsc
   Use petscvec
   Use petscmat
   Use petscksp
   Use petscmesh
   
   Implicit NONE   


   Type(AppCtx_Type)                            :: AppCtx
   PetscInt                                     :: iErr
   Character(len=MEF90_MXSTRLEN)                :: IOBuffer
   PetscInt                                     :: iter, iTS

   Call VarFracQSInit(AppCtx)
   
   If (AppCtx%AppParam%verbose > 1) Then
      Call EXOView(AppCtx%EXO, AppCtx%AppParam%LogViewer)
      Call EXOView(AppCtx%MyEXO, AppCtx%AppParam%MyLogViewer) 
      Call MeshTopologyView(AppCtx%MeshTopology, AppCtx%AppParam%MyLogViewer) 
   End If   
   
   TimeStep: Do iTS = 1, AppCtx%NumTimeSteps
      AppCtx%TimeStep = iTS
      Write(IOBuffer, 99) AppCtx%TimeStep
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
99    Format('\n=== Solving time step ', I4, '\n\n'c)

      !!! Init the fields:
      Call Init_TS_Loads(AppCtx)      
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Done with Init_TS_Loads \n'c 
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If

      Call Init_TS_U(AppCtx)
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Done with Init_TS_U \n'c 
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If

      Call Init_TS_V(AppCtx)
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Done with Init_TS_V \n'c 
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      
      AltMin: Do !iter=1, AppCtx%VarFracSchemeParam%AltMinMaxIter
         Write(IOBuffer, "('Iteration ', I4, ' /', I4, A)") iTS, iter,'\n'c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   
         !------------------------------------------------------------------- 
         ! Problem for U
         !-------------------------------------------------------------------
         If (AppCtx%AppParam%verbose > 0) Then
            Write(IOBuffer, *) 'Assembling the Matrix and RHS for  the U-subproblem \n'c 
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         End If
         Call RHSU_Assembly(AppCtx)
         If (AppCtx%AppParam%verbose > 2) Then
            Call VecView(AppCtx%RHSU, AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
         End If
         
         
         Call MatU_Assembly(AppCtx)
         If (AppCtx%AppParam%verbose > 2) Then
            Call MatView(AppCtx%KU, AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
         End If
               
         If (AppCtx%AppParam%verbose > 0) Then
            Write(IOBuffer, *) 'Calling KSPSolve for the U-subproblem\n'c 
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
         End If
         Call Solve_U(AppCtx)
         
         !------------------------------------------------------------------- 
         ! Problem for V
         !-------------------------------------------------------------------
         If (AppCtx%AppParam%verbose > 0) Then
            Write(IOBuffer, *) 'Assembling the Matrix and RHS for the V-subproblem \n'c 
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
         End If
         Call RHSV_Assembly(AppCtx)
         If (AppCtx%AppParam%verbose > 2) Then
            Call VecView(AppCtx%RHSV, AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
         End If
         Call MatV_Assembly(AppCtx)
         If (AppCtx%AppParam%verbose > 2) Then
            Call MatView(AppCtx%KV, AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
         End If
         
         If (AppCtx%AppParam%verbose > 0) Then
            Write(IOBuffer, *) 'Calling KSPSolve for the V-subproblem\n'c 
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
         End If
         Call Solve_V(AppCtx)
   
         !------------------------------------------------------------------- 
         ! Check the exit condition: tolerance on the error in V 
         !------------------------------------------------------------------- 
         If ( (Mod(iter, AppCtx%VarFracSchemeParam%AltMinSaveInt) == 0) .OR. (AppCtx%ErrV < AppCtx%VarFracSchemeParam%AltMinTol)) Then
            If (AppCtx%AppParam%verbose > 0) Then
               Write(IOBuffer, *) 'Saving U and V\n'c
               Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
            End If   
            Call Save_U(AppCtx)
            Call Save_V(AppCtx)

            If (AppCtx%AppParam%verbose > 0) Then
               Write(IOBuffer, *) 'Computing bulk energy, strains and stresses and saving\n'c 
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
            If ( (AppCtx%VarFracSchemeParam%SaveStress) .OR. (AppCtx%VarFracSchemeParam%SaveStrain) ) Then
               Call ComputeStrainStress(AppCtx)
               Call Save_StrainStress(AppCtx)
            End If
         End If
         If (AppCtx%ErrV < AppCtx%VarFracSchemeParam%AltMinTol) then 
            EXIT 
         End If
         iter = iter + 1
      End Do AltMin
      
      !------------------------------------------------------------------- 
      ! Save the results
      !-------------------------------------------------------------------
      
      

100   Format('Elastic energy:       ', ES12.5, '\n'c)    
101   Format('External Forces Work: ', ES12.5, '\n'c)    
102   Format('Surface energy:       ', ES12.5, '\n'c)    
103   Format('Total energy:         ', ES12.5, '\n'c)    

   End Do TimeStep

   Call VarFracQSFinalize(AppCtx)
End Program  VarFracQS
