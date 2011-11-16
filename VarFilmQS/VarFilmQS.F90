Program  VarFilmQS

#include "finclude/petscdef.h"

Use m_VarFilmQS
Use m_MEF90
Use m_Film_Struct
Use m_VarFilmQS_W

Implicit NONE   

Type(AppCtx_Type)                            :: AppCtx
PetscInt                                     :: iErr, iBlk
PetscInt                                     :: iDebug
PetscInt                                     :: iBTStep
Character(len=MEF90_MXSTRLEN)                :: IOBuffer
PetscInt                                     :: AltMinIter
Character(len=MEF90_MXSTRLEN)                :: filename
Character(len=MEF90_MXSTRLEN), Dimension(5)  :: stagename
PetscLogDouble                               :: CurrentMemoryUsage, MaximumMemoryUsage
PetscBool                                    :: restart

Call VarFilmQSInit(AppCtx)

If (AppCtx%AppParam%verbose > 1) Then
   Call EXOView(AppCtx%EXO, AppCtx%AppParam%LogViewer)
   Call EXOView(AppCtx%MyEXO, AppCtx%AppParam%MyLogViewer) 
   Call MeshTopologyView(AppCtx%MeshTopology, AppCtx%AppParam%MyLogViewer) 
End If   

iDebug = 0
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
   
   !!! Update U at fixed nodes
   Call Init_TS_U(AppCtx)
   
   !! Rebuild AppCtx%VIrrev and AppCtx%IrrevFlag
   Call Update_Irrev(AppCtx)
   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Done with Update_Irrev \n' 
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If

! 
! 	Write(IOBuffer, *) 'W Section prior to Update Irrev W \n' 
! 	Call SectionRealView(AppCtx%V%Sec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)

   !! Update WBCFlag accounting for irreversibility
   Call Update_IrrevW(AppCtx)
   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Done with Update_IrrevW \n' 
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If

   !!! Update V at fixed nodes
   Call Init_TS_V(AppCtx)
   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Done with Init_TS_V \n' 
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   !!! Update W at fixed nodes
   Call Init_TS_W(AppCtx)
   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Done with Init_TS_W \n' 
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If

   
   
   AltMinIter = 1
   AltMin: Do 
      If (AppCtx%AppParam%verbose > 0) Then
         Call PetscMemoryGetCurrentUsage(CurrentMemoryUsage,iErr); CHKERRQ(iErr)
         Call PetscMemoryGetMaximumUsage(MaximumMemoryUsage,iErr); CHKERRQ(iErr)
         Write(MEF90_MyRank+100, *) AppCtx%TimeStep, AltMinIter, CurrentMemoryUsage, MaximumMemoryUsage
         Call PetscMemoryShowUsage(PETSC_VIEWER_STDOUT_WORLD, "PetscMemoryShowUsage output for PETSC_COMM_WORLD: ", iErr); CHKERRQ(iErr)
      End If
      Write(IOBuffer, "('Iteration ', I4, ' /', I4, A)") AppCtx%TimeStep, AltMinIter,'\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

	Select Case(AppCtx%VarFracSchemeParam%CoupledUW)
		Case(VarFrac_UW_Coupled)
		If (AppCtx%AppParam%verbose > 0) Then
			Write(IOBuffer, *) 'Alternate UW \n' 
			Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
		End If
			Call Step_UW(AppCtx)
		Case(VarFrac_UW_Uncoupled)
			Call Step_U(AppCtx)
			Call Step_W(AppCtx)
	End Select
      !------------------------------------------------------------------- 
      ! Problem for U
      !-------------------------------------------------------------------
!       Call Step_U(AppCtx)
      !------------------------------------------------------------------- 
      ! Problem for V
      !-------------------------------------------------------------------
      Call Step_V(AppCtx)
      !------------------------------------------------------------------- 
      ! Problem for W
      !-------------------------------------------------------------------
!       Call Step_W(AppCtx)
      
      !------------------------------------------------------------------- 
      ! Check the exit condition: tolerance on the error in V 
      !------------------------------------------------------------------- 
      If ( (Mod(AltMinIter, AppCtx%VarFracSchemeParam%AltMinSaveInt) == 0) .OR. (AppCtx%ErrV < AppCtx%VarFracSchemeParam%AltMinTol) .OR. (AltMinIter == AppCtx%VarFracSchemeParam%AltMinMaxIter)) Then
         If (AppCtx%AppParam%verbose > 0) Then
            Write(IOBuffer, *) 'Saving U, V and W\n'
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
         End If   
         Call Save_U(AppCtx)
         Call Save_V(AppCtx)
         Call Save_W(AppCtx)
         
         If (AppCtx%AppParam%verbose > 0) Then
            Write(IOBuffer, *) 'Computing bulk energy, strains and stresses and saving\n' 
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
         End If
         Call ComputeEnergies(AppCtx)
         
	Write(IOBuffer, 104) AppCtx%Load(AppCtx%TimeStep)
	Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

	Write(IOBuffer, 100) AppCtx%ElasticEnergy(AppCtx%TimeStep)
	Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

	Write(IOBuffer, 109) AppCtx%FilmEnergy(AppCtx%TimeStep)
	Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
	
	Write(IOBuffer, 108) AppCtx%BondingLayerEnergy(AppCtx%TimeStep)
	Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

	Write(IOBuffer, 102) AppCtx%FractureEnergy(AppCtx%TimeStep)
	Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
	
	Write(IOBuffer, 107) AppCtx%DelaminationEnergy(AppCtx%TimeStep)
	Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

	Write(IOBuffer, 103) AppCtx%TotalEnergy(AppCtx%TimeStep)
	Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
	
	Write(IOBuffer, 101) AppCtx%ExtForcesWork(AppCtx%TimeStep)
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
   ! Save the results
   !-------------------------------------------------------------------
   
   If ( AppCtx%TimeStep == AppCtx%NumTimeSteps ) Then
      EXIT
   End If
   AppCtx%TimeStep = AppCtx%TimeStep + 1
   Write(filename, 105) Trim(AppCtx%AppParam%prefix)
  ! Call PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, LogViewer, ierr);CHKERRQ(ierr)
   ! Call PetscViewerASCIIPrintf(LogViewer,filename,ierr);CHKERRQ(ierr)
!   Call PetscLogView(LogViewer,ierr);CHKERRQ(ierr)
  !  Call PetscViewerDestroy(LogViewer,ierr);CHKERRQ(ierr)
End Do TimeStep

100   Format('Elastic energy:           ', ES12.5, '\n')    
101   Format('External Forces Work:     ', ES12.5, '\n')    
102   Format('Surface energy:           ', ES12.5, '\n')    
103   Format('Total energy:             ', ES12.5, '\n')    
104   Format('Load:                     ', ES12.5, '\n')    
107   Format('Delamination energy:      ', ES12.5, '\n')    
108   Format('Bonding Layer energy:     ', ES12.5, '\n')    
109   Format('Film energy:              ', ES12.5, '\n')    
110   Format('Elastic energy:           ', ES12.5, '\n')    

105   Format(A,'-logsummary.txt')


   Call VarFracQSFinalize(AppCtx)
End Program  VarFilmQS
