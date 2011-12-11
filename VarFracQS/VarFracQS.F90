Program  VarFracQS

#include "finclude/petscdef.h"

#if defined PB_2D
   Use m_VarFracQS2D
   Use m_Poisson2D
   Use m_TransientHeat2D
#elif defined PB_3D
   Use m_VarFracQS3D
   Use m_Poisson3D
   Use m_TransientHeat3D
#endif   
   Use m_MEF90
   Use m_VarFrac_Struct
   Use m_Heat_Struct
   
   Implicit NONE   

   Type(AppCtx_Type)                            :: AppCtx
   PetscInt                                     :: iErr, iBlk
   PetscInt                                     :: iDebug
   PetscInt                                     :: iBTStep
   Character(len=MEF90_MXSTRLEN)                :: IOBuffer
   PetscInt                                     :: AltMinIter
   Character(len=MEF90_MXSTRLEN)                :: filename
   Type(PetscViewer)                            :: LogViewer
   Character(len=MEF90_MXSTRLEN), Dimension(4)  :: stagename
   PetscLogDouble                               :: CurrentMemoryUsage, MaximumMemoryUsage
   PetscBool                                    :: restart
   
   PetscInt                                     :: StepOUT
   !PetscBool                                    :: AppCtx%IsBT

   Character(len=MEF90_MXSTRLEN)                :: TCST_FileName
   PetscInt, Dimension(:), Pointer              :: SizeScal
   Type(Heat_AppCtx_Type)                       :: HeatAppCtx

   Call VarFracQSInit(AppCtx)
   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Done with VarFracQSInit'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   If (AppCtx%AppParam%verbose > 1) Then
      Call EXOView(AppCtx%EXO, AppCtx%AppParam%LogViewer)
      Call EXOView(AppCtx%MyEXO, AppCtx%AppParam%MyLogViewer) 
      Call MeshTopologyView(AppCtx%MeshTopology, AppCtx%AppParam%MyLogViewer) 
   End If   

   !Read Heat material parameters from *.TSCT file 
   !!! Read Mat Properties from the CST file
   TCST_FileName = trim(AppCtx%AppParam%prefix)//'.TCST'
   Call MatHeat_Read(AppCtx%MeshTopology, HeatAppCtx%MatProp, TCST_FileName)
   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) "Done with MatHeat_Read\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If


   !Init For HEAT
   !HeatAppCtx%NumSteps = NumSteps
   HeatAppCtx%NumSteps = 100
   HeatAppCtx%maxsteps = 1000 
   HeatAppCtx%VertVar_Temperature = AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Temperature)%Offset
   Call ElementInit(AppCtx%MeshTopology, HeatAppCtx%Elem, 2)
  
   Allocate(SizeScal(1)) 
   SizeScal=1
   Call FieldCreateVertex(HeatAppCtx%U,     'T',   AppCtx%MeshTopology, SizeScal)
   Call FieldCreateVertex(HeatAppCtx%F,     'F',  AppCtx%MeshTopology,  SizeScal)
   Call FieldCreateVertex(HeatAppCtx%RHS,   'RHS', AppCtx%MeshTopology, SizeScal)
   Call FlagCreateVertex(HeatAppCtx%BCFlag, 'BC',   AppCtx%MeshTopology, SizeScal)
   Call FieldCreateVertex(HeatAppCtx%UBC,    'TBC',      AppCtx%MeshTopology, SizeScal)
   DeAllocate(SizeScal)

   !Read Initial Temerature Field 
   Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology,  AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Temperature)%Offset, 1, HeatAppCtx%U)
   Call SectionRealToVec(HeatAppCtx%U%Sec, HeatAppCtx%U%Scatter,  SCATTER_REVERSE, HeatAppCtx%U%Vec, ierr); CHKERRQ(ierr)
   
   !      Call HeatSetInitial(HeatAppCtx, MeshTopology, ValT_Init,ValT_F)
      !Set BC in Temperature
!      Call HeatSetBC(HeatAppCtx, T_BC, MyExo, MeshTopology,VarFrac_NSProp_HasPForce)
   Call SectionIntZero(HeatAppCtx%BCFlag%Sec, iErr); CHKERRQ(iErr)
   Call SectionIntAddNSProperty(HeatAppCtx%BCFlag%Sec,   AppCtx%MyEXO%NSProperty(VarFrac_NSProp_HasPForce), AppCtx%MeshTopology)

   Call Poisson_TSSetUp(HeatAppCtx, AppCtx%MeshTopology)
   Call HeatMatAssembly(HeatAppCtx, AppCtx%MeshTopology)
   Call RHSAssembly(HeatAppCtx, AppCtx%MeshTopology, AppCtx%MyExo)
   Call MatMassAssembly(HeatAppCtx, AppCtx%MeshTopology)


   iDebug = 0
   Write(stagename(1), "(A)") "Outer loop"
   Write(stagename(2), "(A)") "AltMin loop"
   Write(stagename(3), "(A)") "U-step"
   Write(stagename(4), "(A)") "V-step"
   TimeStep: Do 
      !Call ALEStagePush(stagename(1), iDebug, iErr); CHKERRQ(iErr)

!       Write(IOBuffer, 99) AppCtx%TimeStep
!       Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
! 99    Format('\n=== Solving time step ', I4, '\n\n')

      !!! Init the fields:
      Call Init_TS_Loads(AppCtx) 
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_ForceTemp)%Offset, AppCtx%TimeStep, HeatAppCtx%F)
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Temperature)%Offset,  AppCtx%TimeStep, HeatAppCtx%UBC) 

      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Done with Init_TS_Loads \n' 
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If

      !!Compute Heat Field 
      If (AppCtx%TimeStep < 2) Then
!         Call SolveTransientStep(HeatAppCtx, AppCtx%MyEXO, AppCtx%MeshTopology, 0.0_Kr , AppCtx%Load(AppCtx%TimeStep), AppCtx%TimeStep)
      else 
         Call SolveTransientStep(HeatAppCtx, AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%Load(AppCtx%TimeStep-1), AppCtx%Load(AppCtx%TimeStep), AppCtx%TimeStep)
      End If 
!AppCtx%Load is the list of time steps. Quasi-Static ..... 
!Trasnfert HeatAppCtx%U to AppCtx%Theta

      !!! Update U at fixed nodes
      Call Init_TS_U(AppCtx)
      
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
      
      AppCtx%IsBT = .False.
      AltMinIter = 1
      AltMin: Do 
         !Call ALEStagePush(stagename(2), iDebug, iErr); CHKERRQ(iErr)
         If (AppCtx%AppParam%verbose > 0) Then
            Call PetscMemoryGetCurrentUsage(CurrentMemoryUsage,iErr); CHKERRQ(iErr)
            Call PetscMemoryGetMaximumUsage(MaximumMemoryUsage,iErr); CHKERRQ(iErr)
            Call PetscMemoryShowUsage(PETSC_VIEWER_STDOUT_WORLD, "PetscMemoryShowUsage output for PETSC_COMM_WORLD: ", iErr); CHKERRQ(iErr)
         End If
         Write(IOBuffer, "('Iteration ', I4, ' /', I4, A)") AppCtx%TimeStep, AltMinIter,'\n'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   
         !------------------------------------------------------------------- 
         ! Problem for U
         !-------------------------------------------------------------------
         !Call ALEStagePush(stagename(3), iDebug, iErr); CHKERRQ(iErr)
         Call Step_U(AppCtx)
         !Call ALEStagePop(iDebug, iErr); CHKERRQ(iErr)
         !If (AppCtx%AppParam%verbose > 0) Then
         !   Call ALEStagePrintMemory(stagename(3), iErr); CHKERRQ(iErr)
         !End If
         !------------------------------------------------------------------- 
         ! Problem for V
         !-------------------------------------------------------------------
         !Call ALEStagePush(stagename(4), iDebug, iErr); CHKERRQ(iErr)
         Call Step_V(AppCtx)
         !Call ALEStagePop(iDebug, iErr); CHKERRQ(iErr)
         !If (AppCtx%AppParam%verbose > 0) Then
         !   Call ALEStagePrintMemory(stagename(4), iErr); CHKERRQ(iErr)
         !EndIf      
         !------------------------------------------------------------------- 
         ! Check For BackTracking 
         !------------------------------------------------------------------- 
         !
         TestBT1: If ((AppCtx%VarFracSchemeParam%DoBT) .AND.                                                &
             (Mod(AltMinIter, AppCtx%VarFracSchemeParam%BTInt) == 0) .AND.                         &
             (.NOT. AppCtx%IsBT)) Then
            Call ComputeEnergies(AppCtx)
            Call Backtracking(AppCtx,AppCtx%TimeStep,StepOUT,AppCtx%IsBT)
            BTFound1: If (AppCtx%IsBT) Then
               !AppCtx%IsBT = PETSC_TRUE
               AppCtx%TimeStep = StepOUT

               !!! Insert 2 blank lines in the energy files so that gnuplot breaks lines
               If (MEF90_MyRank ==0) Then
                  Write(AppCtx%AppParam%Ener_Unit, *)
                  Write(AppCtx%AppParam%Ener_Unit, *)
                  
                  If (AppCtx%VarFracSchemeParam%SaveBlk) Then
                     Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks_Global
                        Write(AppCtx%AppParam%EnerBlock_Unit(iBlk), *)
                        Write(AppCtx%AppParam%EnerBlock_Unit(iBlk), *)
                     End Do
                  End If
               End If
               !Call ALEStagePop(iDebug, iErr); CHKERRQ(iErr)
               !If (AppCtx%AppParam%verbose > 0) Then
               !   Call ALEStagePrintMemory(stagename(2), iErr); CHKERRQ(iErr)
               !End If
               !!! Exit AltMin loop
               EXIT
            End If BTFound1
         End If TestBT1

         !------------------------------------------------------------------- 
         ! Check the exit condition: tolerance on the error in V 
         !------------------------------------------------------------------- 
         If ((Mod(AltMinIter, AppCtx%VarFracSchemeParam%AltMinSaveInt) == 0) .OR.                  &
             (AppCtx%ErrV < AppCtx%VarFracSchemeParam%AltMinTol) .OR.                              &
             (AltMinIter == AppCtx%VarFracSchemeParam%AltMinMaxIter)) Then
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

            If ((AppCtx%VarFracSchemeParam%SaveStress) .OR.                                        &
                (AppCtx%VarFracSchemeParam%SaveStrain) ) Then
               Call ComputeStrainStress(AppCtx)
               Call Save_StrainStress(AppCtx)
            End If
         End If
         If ((AppCtx%ErrV < AppCtx%VarFracSchemeParam%AltMinTol) .OR.                              &
             (AltMinIter == AppCtx%VarFracSchemeParam%AltMinMaxIter) ) then 
            Call Save_Ener(AppCtx)
            !Call ALEStagePop(iDebug, iErr); CHKERRQ(iErr)
            !If (AppCtx%AppParam%verbose > 0) Then
            !   Call ALEStagePrintMemory(stagename(2), iErr); CHKERRQ(iErr)
            !End If
            EXIT 
         End If
         AltMinIter = AltMinIter + 1
         !Call ALEStagePop(iDebug, iErr); CHKERRQ(iErr)
         !If (AppCtx%AppParam%verbose > 0) Then
         !   Call ALEStagePrintMemory(stagename(2), iErr); CHKERRQ(iErr)
         !End If
      End Do AltMin
      !------------------------------------------------------------------- 
      ! Check For BackTracking again
      !------------------------------------------------------------------- 
      TestBTFinal: If ((AppCtx%VarFracSchemeParam%DoBT) .AND. (.NOT. AppCtx%IsBT)) Then
         Call Backtracking(AppCtx,AppCtx%TimeSTep,StepOUT,AppCtx%IsBT)
         BTFoundFinal: If (AppCtx%IsBT) Then
            !AppCtx%IsBT = PETSC_TRUE
            AppCtx%TimeStep = StepOUT

            !!! Insert 2 blank lines in the energy files so that gnuplot breaks lines
            If (MEF90_MyRank ==0) Then
               Write(AppCtx%AppParam%Ener_Unit, *)
               Write(AppCtx%AppParam%Ener_Unit, *)
               
               If (AppCtx%VarFracSchemeParam%SaveBlk) Then
                  Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks_Global
                     Write(AppCtx%AppParam%EnerBlock_Unit(iBlk), *)
                     Write(AppCtx%AppParam%EnerBlock_Unit(iBlk), *)
                  End Do
               End If
            End If
            !Call ALEStagePop(iDebug, iErr); CHKERRQ(iErr)
            !If (AppCtx%AppParam%verbose > 0) Then
            !   Call ALEStagePrintMemory(stagename(2), iErr); CHKERRQ(iErr)
            !End If
         Else 
            AppCtx%TimeStep = AppCtx%TimeStep +1
         End If BTFoundFinal
      Else
         If (.NOT. AppCtx%IsBT) Then
            AppCtx%TimeStep = AppCtx%TimeStep +1      
         End If
      End If TestBTFinal

      !------------------------------------------------------------------- 
      ! Save the results
      !-------------------------------------------------------------------
      Write(filename, 105) Trim(AppCtx%AppParam%prefix)
      Call PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, LogViewer, ierr);CHKERRQ(ierr)
      Call PetscViewerASCIIPrintf(LogViewer,filename,ierr);CHKERRQ(ierr)
      Call PetscLogView(LogViewer,ierr);CHKERRQ(ierr)
      Call PetscViewerDestroy(LogViewer,ierr);CHKERRQ(ierr)
      !Call ALEStagePop(iDebug, iErr); CHKERRQ(iErr)
      !If (AppCtx%AppParam%verbose > 0) Then
      !   Call ALEStagePrintMemory(stagename(1), iErr); CHKERRQ(iErr)
      !End If
      If ( AppCtx%TimeStep > AppCtx%NumTimeSteps ) Then
         EXIT
      End If
   End Do TimeStep
   

100   Format('Elastic energy:       ', ES12.5, '\n')    
101   Format('External Forces Work: ', ES12.5, '\n')    
102   Format('Surface energy:       ', ES12.5, '\n')    
103   Format('Total energy:         ', ES12.5, '\n')    
104   Format('Load:                 ', ES12.5, '\n')    

105   Format(A,'-logsummary.txt')
   Call VarFracQSFinalize(AppCtx)
End Program  VarFracQS
