Program  VarFracQS

#include "finclude/petscdef.h"

#if defined PB_2D
   Use m_VarFracQS2D
#elif defined PB_3D
   Use m_VarFracQS3D
#endif   
   Use m_MEF90
   Use m_VarFrac_Struct
   
   Implicit NONE   

   Type(AppCtx_Type)                            :: AppCtx
   PetscInt                                     :: iErr, iBlk
   PetscInt                                     :: iDebug
   PetscInt                                     :: iBTStep
   Character(len=MEF90_MXSTRLEN)                :: IOBuffer
   PetscInt                                     :: AltMinIter
   Character(len=MEF90_MXSTRLEN)                :: filename
   Character(len=MEF90_MXSTRLEN), Dimension(4)  :: stagename
   PetscLogDouble                               :: CurrentMemoryUsage, MaximumMemoryUsage
   PetscTruth                                   :: restart

   Call VarFracQSInit(AppCtx)
   
   If (AppCtx%AppParam%verbose > 1) Then
      Call EXOView(AppCtx%EXO, AppCtx%AppParam%LogViewer)
      Call EXOView(AppCtx%MyEXO, AppCtx%AppParam%MyLogViewer) 
      Call MeshTopologyView(AppCtx%MeshTopology, AppCtx%AppParam%MyLogViewer) 
   End If   
   
   iDebug = 0
   Write(stagename(1), "(A)") "Outer loop"
   Write(stagename(2), "(A)") "AltMin loop"
   Write(stagename(3), "(A)") "U-step"
   Write(stagename(4), "(A)") "V-step"
   TimeStep: Do 
      Call ALEStagePush(stagename(1), iDebug, iErr); CHKERRQ(iErr)

      Write(IOBuffer, 99) AppCtx%TimeStep
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
99    Format('\n=== Solving time step ', I4, '\n\n')

      !!! Init the fields:
      Call Init_TS_Loads(AppCtx)      
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Done with Init_TS_Loads \n' 
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If


      !!! Update U at fixed nodes
      Call FieldInsertVertexBoundaryValues(AppCtx%U, AppCtx%UBC, AppCtx%BCUFlag, AppCtx%MeshTopology)
      
     !!! Rebuild AppCtx%VIrrev and AppCtx%IrrevFlag
      Call Update_Irrev(AppCtx)

      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Done with Update_Irrev \n' 
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If

      !!! Update V at fixed nodes
      Call Init_TS_V(AppCtx)

      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Done with Init_TS_V \n' 
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      
      
      AltMinIter = 1
      AltMin: Do 
         Call ALEStagePush(stagename(2), iDebug, iErr); CHKERRQ(iErr)
         If (AppCtx%AppParam%verbose > 0) Then
            Call PetscMemoryGetCurrentUsage(CurrentMemoryUsage,iErr); CHKERRQ(iErr)
            Call PetscMemoryGetMaximumUsage(MaximumMemoryUsage,iErr); CHKERRQ(iErr)
            Write(MEF90_MyRank+100, *) AppCtx%TimeStep, AltMinIter, CurrentMemoryUsage, MaximumMemoryUsage
            Call PetscMemoryShowUsage(PETSC_VIEWER_STDOUT_WORLD, "PetscMemoryShowUsage output for PETSC_COMM_WORLD: ", iErr); CHKERRQ(iErr)
         End If
         Write(IOBuffer, "('Iteration ', I4, ' /', I4, A)") AppCtx%TimeStep, AltMinIter,'\n'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   
         !------------------------------------------------------------------- 
         ! Problem for U
         !-------------------------------------------------------------------
         Call ALEStagePush(stagename(3), iDebug, iErr); CHKERRQ(iErr)
         Call Step_U(AppCtx)
         Call ALEStagePop(iDebug, iErr); CHKERRQ(iErr)
         If (AppCtx%AppParam%verbose > 0) Then
            Call ALEStagePrintMemory(stagename(3), iErr); CHKERRQ(iErr)
         End If
         !------------------------------------------------------------------- 
         ! Problem for V
         !-------------------------------------------------------------------
         Call ALEStagePush(stagename(4), iDebug, iErr); CHKERRQ(iErr)
         Call Step_V(AppCtx)
         Call ALEStagePop(iDebug, iErr); CHKERRQ(iErr)
         If (AppCtx%AppParam%verbose > 0) Then
            Call ALEStagePrintMemory(stagename(4), iErr); CHKERRQ(iErr)
         EndIf      
         !------------------------------------------------------------------- 
         ! Check For BackTracking 
         !------------------------------------------------------------------- 
         !
         If ((AppCtx%VarFracSchemeParam%DoBT) .AND. (Mod(AltMinIter, AppCtx%VarFracSchemeParam%BTInt) == 0) .AND. (.NOT. AppCtx%IsBT)) Then
            Call ComputeEnergies(AppCtx)
            Call BackTracking(AppCtx, iBTStep)
            AppCtx%IsBT = PETSC_FALSE   
            If (iBTStep < AppCtx%TimeStep) Then
               AppCtx%IsBT = PETSC_TRUE
               AppCtx%TimeStep = max(1, iBTStep-1)

               If (MEF90_MyRank ==0) Then
                  !!! Insert 2 blank lines in the energy files so that gnuplot breaks lines
                  Write(AppCtx%AppParam%Ener_Unit, *)
                  Write(AppCtx%AppParam%Ener_Unit, *)
                  
                  If (AppCtx%VarFracSchemeParam%SaveBlk) Then
                     Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks_Global
                        Write(AppCtx%AppParam%EnerBlock_Unit(iBlk), *)
                        Write(AppCtx%AppParam%EnerBlock_Unit(iBlk), *)
                     End Do
                  End If
               End If
               
               !!! Exit the AltMin loop
               Call ALEStagePop(iDebug, iErr); CHKERRQ(iErr)
               If (AppCtx%AppParam%verbose > 0) Then
                  Call ALEStagePrintMemory(stagename(2), iErr); CHKERRQ(iErr)
               End If
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
            Call ComputeEnergies(AppCtx)

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
            Call ALEStagePop(iDebug, iErr); CHKERRQ(iErr)
            If (AppCtx%AppParam%verbose > 0) Then
               Call ALEStagePrintMemory(stagename(2), iErr); CHKERRQ(iErr)
            End If
            EXIT 
         End If
         AltMinIter = AltMinIter + 1
         Call ALEStagePop(iDebug, iErr); CHKERRQ(iErr)
         If (AppCtx%AppParam%verbose > 0) Then
            Call ALEStagePrintMemory(stagename(2), iErr); CHKERRQ(iErr)
         End If
      End Do AltMin
      !------------------------------------------------------------------- 
      ! Check For BackTracking again
      !------------------------------------------------------------------- 
      If ((AppCtx%VarFracSchemeParam%DoBT) .AND. (.NOT. AppCtx%IsBT)) Then
         Call ComputeEnergies(AppCtx)
         Call BackTracking(AppCtx, iBTStep)
         
         If (iBTStep < AppCtx%TimeStep) Then
            AppCtx%IsBT = PETSC_TRUE
            AppCtx%TimeStep = max(1, iBTStep-1)
            !!! Insert 2 blank lines in the energy file so that gnuplot breaks lines
            If (MEF90_MyRank ==0) Then
               !!! Insert 2 blank lines in the energy files so that gnuplot breaks lines
               Write(AppCtx%AppParam%Ener_Unit, *)
               Write(AppCtx%AppParam%Ener_Unit, *)
            
               If (AppCtx%VarFracSchemeParam%SaveBlk) Then
                  Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks_Global
                     Write(AppCtx%AppParam%EnerBlock_Unit(iBlk), *)
                     Write(AppCtx%AppParam%EnerBlock_Unit(iBlk), *)
                  End Do
               End If
            End If
         Else
            AppCtx%IsBT = PETSC_FALSE
         End If
      Else 
         AppCtx%IsBT = PETSC_FALSE
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
      Call ALEStagePop(iDebug, iErr); CHKERRQ(iErr)
      If (AppCtx%AppParam%verbose > 0) Then
         Call ALEStagePrintMemory(stagename(1), iErr); CHKERRQ(iErr)
      End If
   End Do TimeStep

100   Format('Elastic energy:       ', ES12.5, '\n')    
101   Format('External Forces Work: ', ES12.5, '\n')    
102   Format('Surface energy:       ', ES12.5, '\n')    
103   Format('Total energy:         ', ES12.5, '\n')    
104   Format('Load:                 ', ES12.5, '\n')    

105   Format(A,'-logsummary.txt')
106   Format('Time Loop')

   Call VarFracQSFinalize(AppCtx)
End Program  VarFracQS
