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
   PetscInt                                     :: iBTStep
   Character(len=MEF90_MXSTRLEN)                :: IOBuffer
   PetscInt                                     :: AltMinIter
   Character(len=MEF90_MXSTRLEN)                :: filename

   Call VarFracQSInit(AppCtx)
   
   If (AppCtx%AppParam%verbose > 1) Then
      Call EXOView(AppCtx%EXO, AppCtx%AppParam%LogViewer)
      Call EXOView(AppCtx%MyEXO, AppCtx%AppParam%MyLogViewer) 
      Call MeshTopologyView(AppCtx%MeshTopology, AppCtx%AppParam%MyLogViewer) 
   End If   
   
   AppCtx%TimeStep = 1
   TimeStep: Do 
      Write(IOBuffer, 99) AppCtx%TimeStep
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
99    Format('\n=== Solving time step ', I4, '\n\n')

      !!! Init the fields:
      Call Init_TS_Loads(AppCtx)      
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Done with Init_TS_Loads \n' 
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If

      Call Init_TS_U(AppCtx)
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Done with Init_TS_U \n' 
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If

      Call Init_TS_V(AppCtx)
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Done with Init_TS_V \n' 
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      
      AltMinIter = 1
      AltMin: Do 
         Write(IOBuffer, "('Iteration ', I4, ' /', I4, A)") AppCtx%TimeStep, AltMinIter,'\n'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   
         !------------------------------------------------------------------- 
         ! Problem for U
         !-------------------------------------------------------------------
         If (AppCtx%AppParam%verbose > 0) Then
            Write(IOBuffer, *) 'Assembling the Matrix and RHS for  the U-subproblem \n' 
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
            Write(IOBuffer, *) 'Calling KSPSolve for the U-subproblem\n' 
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
         End If
         Call Solve_U(AppCtx)
         
         !------------------------------------------------------------------- 
         ! Problem for V
         !-------------------------------------------------------------------
         If (AppCtx%AppParam%verbose > 0) Then
            Write(IOBuffer, *) 'Assembling the Matrix and RHS for the V-subproblem \n' 
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
            Write(IOBuffer, *) 'Calling KSPSolve for the V-subproblem\n' 
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
         End If
         Call Solve_V(AppCtx)
         
         !------------------------------------------------------------------- 
         ! Check For BackTracking 
         !------------------------------------------------------------------- 
         AppCtx%IsBT = PETSC_FALSE
         If ((AppCtx%VarFracSchemeParam%DoBT) .AND. (Mod(AltMinIter, AppCtx%VarFracSchemeParam%BTInt) == 0) ) Then
            Call ComputeEnergy(AppCtx)
            Call BackTracking(AppCtx, iBTStep)
            
            If (iBTStep < AppCtx%TimeStep) Then
               AppCtx%IsBT = PETSC_TRUE
               AppCtx%TimeStep = iBTStep - 1

               !!! Insert 2 blank lines in the energy file so that gnuplot breaks lines
               Write(AppCtx%AppParam%Ener_Unit, *)
               Write(AppCtx%AppParam%Ener_Unit, *)
               
               !!! Exit the AltMin loop
               EXIT 
            End If
         End If
   
         !------------------------------------------------------------------- 
         ! Check the exit condition: tolerance on the error in V 
         !------------------------------------------------------------------- 
         If ( (Mod(AltMinIter, AppCtx%VarFracSchemeParam%AltMinSaveInt) == 0) .OR. (AppCtx%ErrV < AppCtx%VarFracSchemeParam%AltMinTol) .OR. (AltMinIter == AppCtx%VarFracSchemeParam%AltMinMaxIter)) Then
            If (AppCtx%AppParam%verbose > 0) Then
               Write(IOBuffer, *) 'Saving U and V\n'
               Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
            End If   
            Call Save_U(AppCtx)
            Call Save_V(AppCtx)

            If (AppCtx%AppParam%verbose > 0) Then
               Write(IOBuffer, *) 'Computing bulk energy, strains and stresses and saving\n' 
               Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
            End If
            Call ComputeEnergy(AppCtx)
            Write(IOBuffer, 104) AppCtx%Load(AppCtx%TimeStep)
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
            Write(IOBuffer, 100) AppCtx%ElasticEnergy(AppCtx%TimeStep)
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
            Write(IOBuffer, 101) AppCtx%ExtForcesWork(AppCtx%TimeStep)
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
            Write(IOBuffer, 102) AppCtx%SurfaceEnergy(AppCtx%TimeStep)
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
            Write(IOBuffer, 103) AppCtx%TotalEnergy(AppCtx%TimeStep)
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
            If ( (AppCtx%VarFracSchemeParam%SaveStress) .OR. (AppCtx%VarFracSchemeParam%SaveStrain) ) Then
               Call ComputeStrainStress(AppCtx)
               Call Save_StrainStress(AppCtx)
            End If
         End If
         If ( (AppCtx%ErrV < AppCtx%VarFracSchemeParam%AltMinTol) .OR. (AltMinIter == AppCtx%VarFracSchemeParam%AltMinMaxIter) ) then 
            Call Save_Ener(AppCtx)
            EXIT 
         End If
         AltMinIter = AltMinIter + 1
      End Do AltMin
      !------------------------------------------------------------------- 
      ! Check For BackTracking again
      !------------------------------------------------------------------- 
      If ((AppCtx%VarFracSchemeParam%DoBT) .AND. (.NOT. AppCtx%IsBT) ) Then
         Call ComputeEnergy(AppCtx)
         Call BackTracking(AppCtx, iBTStep)
         
         If (iBTStep < AppCtx%TimeStep) Then
            AppCtx%IsBT = PETSC_TRUE
            AppCtx%TimeStep = iBTStep - 1

            !!! Insert 2 blank lines in the energy file so that gnuplot breaks lines
            Write(AppCtx%AppParam%Ener_Unit, *)
            Write(AppCtx%AppParam%Ener_Unit, *)
         Else
            AppCtx%IsBT = PETSC_FALSE
         End If
      End If
      
      !------------------------------------------------------------------- 
      ! Save the results
      !-------------------------------------------------------------------
      
      If ( AppCtx%TimeStep == AppCtx%NumTimeSteps ) Then
         EXIT
      End If
      AppCtx%TimeStep = AppCtx%TimeStep + 1
      Write(filename, 105) Trim(AppCtx%AppParam%prefix)
      Call PetscLogPrintSummary(PETSC_COMM_WORLD, filename, iErr); CHKERRQ(iErr)
   End Do TimeStep

100   Format('Elastic energy:       ', ES12.5, '\n')    
101   Format('External Forces Work: ', ES12.5, '\n')    
102   Format('Surface energy:       ', ES12.5, '\n')    
103   Format('Total energy:         ', ES12.5, '\n')    
104   Format('Load:                 ', ES12.5, '\n')    

105   Format(A,'-logsummary.txt')

   Call VarFracQSFinalize(AppCtx)
End Program  VarFracQS
