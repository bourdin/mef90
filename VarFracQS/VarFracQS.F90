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
   PetscInt                                     :: AltMinIter, iBTStep
   PetscReal                                    :: EnerBT

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
      
      AltMinIter = 1
      AltMin: Do 
         Write(IOBuffer, "('Iteration ', I4, ' /', I4, A)") AppCtx%TimeStep, AltMinIter,'\n'c
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
         ! Check For BackTracking 
         !------------------------------------------------------------------- 
         If ((AppCtx%VarFracSchemeParam%DoBT) .AND. (Mod(AltMinIter, AppCtx%VarFracSchemeParam%BTInt) == 0) ) Then
            AppCtx%IsBT = PETSC_FALSE
            !!! Check the BT condition
            If (AppCtx%AppParam%verbose > 0) Then
               Write(IOBuffer, *) 'Doing BackTracking\n'c 
               Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
            End If
            Do iBTStep = 1, AppCtx%TimeStep-1
               EnerBT = AppCtx%Load(iBTStep)**2 / AppCtx%Load(AppCtx%TimeStep)**2 * AppCtx%ElasticEnergy(AppCtx%TimeStep) + AppCtx%Load(iBTStep) / AppCtx%Load(AppCtx%TimeStep) * AppCtx%ExtForcesWork(AppCtx%TimeStep) + AppCtx%SurfaceEnergy(AppCtx%TimeStep)
               If ( (AppCtx%TotalEnergy(iBTStep) - EnerBT ) / EnerBT > AppCtx%VarFracSchemeParam%BTTol ) Then
                  Write(IOBuffer, *) '*********** BackTracking to step', iBTStep, '\n'c 
                  Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
   
                  !!! Insert 2 blank lines in the energy file so that gnuplot breaks lines
                  Write(AppCtx%AppParam%Ener_Unit, *)
                  Write(AppCtx%AppParam%Ener_Unit, *)

                  AppCtx%TimeStep = iBTStep - 1
                  ! Since AppCtx%TimeStep is going to be incremented
                  AppCtx%IsBT = PETSC_TRUE
                  !!! Exit the BT loop
                  EXIT
               End If
            End Do
            If (AppCtx%IsBT) Then
               !!! Exit the AltMin loop
               EXIT
            End If
         End If
   
         !------------------------------------------------------------------- 
         ! Check the exit condition: tolerance on the error in V 
         !------------------------------------------------------------------- 
         If ( (Mod(AltMinIter, AppCtx%VarFracSchemeParam%AltMinSaveInt) == 0) .OR. (AppCtx%ErrV < AppCtx%VarFracSchemeParam%AltMinTol) .OR. (AltMinIter == AppCtx%VarFracSchemeParam%AltMinMaxIter)) Then
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
            EXIT 
         End If
         AltMinIter = AltMinIter + 1
      End Do AltMin
      
      !------------------------------------------------------------------- 
      ! Save the results
      !-------------------------------------------------------------------
      
      If ( AppCtx%TimeStep == AppCtx%NumTimeSteps ) Then
         EXIT
      End If
      AppCtx%TimeStep = AppCtx%TimeStep + 1
   End Do TimeStep

100   Format('Elastic energy:       ', ES12.5, '\n'c)    
101   Format('External Forces Work: ', ES12.5, '\n'c)    
102   Format('Surface energy:       ', ES12.5, '\n'c)    
103   Format('Total energy:         ', ES12.5, '\n'c)    
104   Format('Load:                 ', ES12.5, '\n'c)    


   Call VarFracQSFinalize(AppCtx)
End Program  VarFracQS
