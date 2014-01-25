#include "../MEF90/mef90.inc"
Program ThermoElasticity
#include <finclude/petscdef.h>
   Use m_MEF90
   Use m_MEF90_DefMechCtx
   Use m_MEF90_DefMech
   Use m_MEF90_HeatXfer
   Use m_MEF90_HeatXferCtx
   Use petsc
   Implicit NONE   

   PetscErrorCode                                     :: ierr
   Type(MEF90Ctx_Type),target                         :: MEF90Ctx
   Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions
   Type(MEF90CtxGlobalOptions_Type),Parameter         :: MEF90DefaultGlobalOptions = MEF90CtxGlobalOptions_Type( &
                                                         1,                             & ! verbose
                                                         PETSC_FALSE,                   & ! helponly
                                                         MEF90TimeInterpolation_linear, & ! timeInterpolation
                                                         0.0_Kr,                        & ! timeMin
                                                         1.0_Kr,                        & ! timeMax
                                                         11,                            & ! timeNumStep
                                                         MEF90FileFormat_EXOSingle)       ! fileFormat

   !!! Defect mechanics contexts
   Type(MEF90DefMechCtx_Type)                         :: MEF90DefMechCtx
   Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
   Type(MEF90DefMechGlobalOptions_Type),Parameter     :: MEF90DefMechDefaultGlobalOptions2D = MEF90DefMechGlobalOptions_Type( &
                                                         MEF90DefMech_ModeQuasiStatic, & ! mode
                                                         PETSC_TRUE,              & ! disp_addNullSpace
                                                         3,                       & ! DisplacementOffset
                                                         2,                       & ! DamageOffset
                                                         0,                       & ! boundaryDisplacementOffset
                                                         0,                       & ! boundaryDamageOffset
                                                         1,                       & ! temperatureOffset
                                                         4,                       & ! ForceOffset
                                                         3,                       & ! pressureForceOffset
                                                         0,                       & ! plasticStrainOffset
                                                         0,                       & ! StressOffset
                                                         MEF90Scaling_Linear,     & ! boundaryDisplacementScaling
                                                         MEF90Scaling_CST,        & ! boundaryDamageScaling
                                                         MEF90Scaling_Linear,     & ! ForceScaling
                                                         MEF90Scaling_Linear,     & ! pressureForceScaling
                                                         1e-4,                    & ! damage_atol
                                                         1000,                    & ! maxit
                                                         0.,                      & ! irrevThres 
                                                         MEF90DefMech_BTTypeNULL, & ! BTType
                                                         -1,                      & ! BTInt
                                                         -1,                      & ! BTScope
                                                         1.0e-2)                    ! BTTol

   Type(MEF90DefMechGlobalOptions_Type),Parameter     :: MEF90DefMechDefaultGlobalOptions3D = MEF90DefMechGlobalOptions_Type( &
                                                         MEF90DefMech_ModeQuasiStatic, & ! mode
                                                         PETSC_TRUE,              & ! disp_addNullSpace
                                                         3,                       & ! DisplacementOffset
                                                         2,                       & ! DamageOffset
                                                         0,                       & ! boundaryDisplacementOffset
                                                         0,                       & ! boundaryDamageOffset
                                                         1,                       & ! temperatureOffset
                                                         4,                       & ! ForceOffset
                                                         3,                       & ! pressureForceOffset
                                                         0,                       & ! plasticStrainOffset
                                                         0,                       & ! StressOffset
                                                         MEF90Scaling_Linear,     & ! boundaryDisplacementScaling
                                                         MEF90Scaling_CST,        & ! boundaryDamageScaling
                                                         MEF90Scaling_Linear,     & ! ForceScaling
                                                         MEF90Scaling_Linear,     & ! pressureForceScaling
                                                         1e-4,                    & ! damage_atol
                                                         1000,                    & ! maxit
                                                         0.,                      & ! irrevThres 
                                                         MEF90DefMech_BTTypeNULL, & ! BTType
                                                         -1,                      & ! BTInt
                                                         -1,                      & ! BTScope
                                                         1.0e-2)                    ! BTTol


   Type(MEF90DefMechCellSetOptions_Type),Parameter    :: MEF90DefMechDefaultCellSetOptions = MEF90DefMechCellSetOptions_Type( &
                                                         -1,                                      & ! elemTypeShortIDDispl will be overriden
                                                         -1,                                      & ! elemTypeShortIDDamage will be overriden
                                                         [0.0_Kr,0.0_Kr,0.0_Kr],                  & ! force
                                                         0.0_Kr,                                  & ! pressureForce
                                                         MEF90DefMech_defectLawElasticity,        & ! defect law
                                                         MEF90DefMech_defectLawGradientDamageAT1, & ! gradientDamageLaw
                                                         MEF90DefMech_defectLawPlasticityVonMises,& ! plasticityLaw
                                                         [PETSC_FALSE,PETSC_FALSE,PETSC_FALSE],   & ! Has Displacement BC
                                                         0.0_Kr,                                  & ! boundary Displacement
                                                         PETSC_FALSE,                             & ! Has Damage BC
                                                         0.0_Kr,                                  & ! Boundary Damage
                                                         1.0D-9)                                    ! residualStiffness
   Type(MEF90DefMechVertexSetOptions_Type),Parameter  :: MEF90DefMechDefaultVertexSetOptions = MEF90DefMechVertexSetOptions_Type( &
                                                         [PETSC_FALSE,PETSC_FALSE,PETSC_FALSE],   & ! Has Displacement BC
                                                         0.0_Kr,                                  & ! boundary Displacement
                                                         PETSC_FALSE,                             & ! Has Damage BC
                                                         0.0_Kr)                                    ! boundary Damage

   !!! HeatXfer contexts
   Type(MEF90HeatXferCtx_Type)                        :: MEF90HeatXferCtx
   Type(MEF90HeatXferGlobalOptions_Type),Pointer      :: MEF90HeatXferGlobalOptions
   Type(MEF90HeatXferGlobalOptions_Type),Parameter    :: MEF90HeatXferDefaultGlobalOptions = MEF90HeatXferGlobalOptions_Type( &
                                                         MEF90HeatXFer_ModeSteadyState, & ! mode
                                                         PETSC_FALSE,         & ! addNullSpace
                                                         1,                   & ! tempOffset
                                                         0.,                  & ! initialTemperature
                                                         MEF90Scaling_Linear, & ! boundaryTempScaling
                                                         0,                   & ! boundaryTempOffset
                                                         MEF90Scaling_Linear, & ! externalTempScaling
                                                         2,                   & ! externalTempOffset
                                                         MEF90Scaling_Linear, & ! fluxScaling
                                                         1)                     ! fluxOffset
   Type(MEF90HeatXferCellSetOptions_Type),Parameter   :: MEF90HeatXferDefaultCellSetOptions = MEF90HeatXferCellSetOptions_Type( &
                                                         -1,            & ! elemTypeShortID will be overriden
                                                         0.0_Kr,        & ! flux
                                                         0.0_Kr,        & ! surfaceThermalConductivity
                                                         0.0_Kr,        & ! externalTemp
                                                         PETSC_FALSE,   & ! Has BC
                                                         0.0_Kr)          ! boundaryTemp
                                                         
   Type(MEF90HeatXferVertexSetOptions_Type),Parameter :: MEF90HeatXferDefaultVertexSetOptions = MEF90HeatXferVertexSetOptions_Type( &
                                                         PETSC_FALSE,   & ! Has BC
                                                         0.0_Kr)          ! boundaryTemp

   
   Type(DM),target                                    :: Mesh
   Type(IS)                                           :: setIS,CellSetGlobalIS
   PetscInt,Dimension(:),Pointer                      :: setID
   PetscInt                                           :: numset,set
   PetscReal,Dimension(:),Pointer                     :: time,energy,work

   Type(SNES)                                         :: snesDisp
   SNESConvergedReason                                :: snesDispConvergedReason
   Type(Vec)                                          :: residualDisp

   Type(SNES)                                         :: snesTemp
   Type(TS)                                           :: tsTemp
   Type(TSAdapt)                                      :: tsAdaptTemp
   Type(Vec)                                          :: residualTemp

   PetscReal                                          :: tsTempInitialStep,tsTempInitialTime
   PetscInt                                           :: tsTempmaxIter
   PetscReal                                          :: t
          
   PetscBool                                          :: flg
   Character(len=MEF90_MXSTRLEN)                      :: IOBuffer
   Type(PetscViewer)                                  :: logViewer
   Integer                                            :: numfield
   
   Integer                                            :: step
   PetscInt                                           :: dim
      
   !!! Initialize MEF90
   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   Call MEF90Initialize(ierr)

   !!! Get all MEF90-wide options
   Call MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90DefaultGlobalOptions,ierr);CHKERRQ(ierr)
   Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)

   !!! Get DM from mesh
   Call MEF90CtxGetDMMeshEXO(MEF90Ctx,Mesh,ierr);CHKERRQ(ierr)
   Call DMMeshGetDimension(Mesh,dim,ierr);CHKERRQ(ierr)
   Call DMMeshSetMaxDof(Mesh,dim,ierr);CHKERRQ(ierr) 
   Call DMSetBlockSize(Mesh,dim,ierr);CHKERRQ(ierr)
   
   !!! Open output file
   Call MEF90CtxOpenEXO(MEF90Ctx,Mesh,ierr)
   
   !!! Create DefMech context, get all DefMech options
   Call MEF90DefMechCtxCreate(MEF90DefMechCtx,Mesh,MEF90Ctx,ierr);CHKERRQ(ierr)
   If (dim == 2) Then
      Call MEF90DefMechCtxSetFromOptions(MEF90DefMechCtx,PETSC_NULL_CHARACTER,MEF90DefMechDefaultGlobalOptions2D, &
                                         MEF90DefMechDefaultCellSetOptions,MEF90DefMechDefaultVertexSetOptions,ierr)
   Else
      Call MEF90DefMechCtxSetFromOptions(MEF90DefMechCtx,PETSC_NULL_CHARACTER,MEF90DefMechDefaultGlobalOptions3D, &
                                         MEF90DefMechDefaultCellSetOptions,MEF90DefMechDefaultVertexSetOptions,ierr)
   End If
   Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)
   
   !!! Create HeatXfer context, get all HeatXfer options
   Call MEF90HeatXferCtxCreate(MEF90HeatXferCtx,Mesh,MEF90Ctx,ierr);CHKERRQ(ierr)
   Call MEF90HeatXferCtxSetFromOptions(MEF90HeatXferCtx,PETSC_NULL_CHARACTER,MEF90HeatXferDefaultGlobalOptions, &
                                       MEF90HeatXferDefaultCellSetOptions,MEF90HeatXferDefaultVertexSetOptions,ierr)
   Call PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr);CHKERRQ(ierr)

   !!! Get material properties bags
   If (dim == 2) Then
      Call MEF90MatPropBagSetFromOptions(MEF90DefMechCtx%MaterialPropertiesBag,MEF90DefMechCtx%DMVect,MEF90Mathium2D,MEF90Ctx,ierr)
   Else
      Call MEF90MatPropBagSetFromOptions(MEF90DefMechCtx%MaterialPropertiesBag,MEF90DefMechCtx%DMVect,MEF90Mathium3D,MEF90Ctx,ierr)
   End If   
   MEF90HeatXferCtx%MaterialPropertiesBag => MEF90DefMechCtx%MaterialPropertiesBag

   !!! Create time array from global options
   Call MEF90CtxGetTime(MEF90Ctx,time,ierr)

   !!! Set the data layout
   Call MEF90DefMechCtxSetSections(MEF90DefMechCtx,ierr)
   Call MEF90HeatXferCtxSetSections(MEF90HeatXferCtx,ierr)

   !!! Create vectors
   Call MEF90DefMechCtxCreateVectors(MEF90DefMechCtx,ierr)
   Call MEF90HeatXferCtxCreateVectors(MEF90HeatXferCtx,ierr)
   
   !!! Link the temperature field from teh DefMech context to that of the HeatXfer context
   DeAllocate(MEF90DefMechCtx%temperature)
   MEF90DefMechCtx%temperature => MEF90HeatXferCtx%temperature

   !!! 
   !!! Create SNES or TS, Mat and set KSP default options
   !!!
   Call VecDuplicate(MEF90DefMechCtx%displacement,residualDisp,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(residualDisp,"residualDisp",ierr);CHKERRQ(ierr)
   Call MEF90DefMechCreateSolversDisp(MEF90DefMechCtx,snesDisp,residualDisp,ierr)

   Call VecDuplicate(MEF90HeatXferCtx%temperature,residualTemp,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(residualTemp,"residualTemp",ierr);CHKERRQ(ierr)
   If (MEF90HeatXferGlobalOptions%mode == MEF90HeatXFer_ModeSteadyState) Then
      Call MEF90HeatXferCreateSNES(MEF90HeatXferCtx,snesTemp,residualTemp,ierr)
   Else
      Call MEF90HeatXferCreateTS(MEF90HeatXferCtx,tsTemp,residualTemp,ierr)
      tsTempInitialStep = (time(size(time))-time(1)) / (size(time) + 0.0_Kr) / 10.0_Kr
      tsTempInitialTime = time(1)
      Call TSSetInitialTimeStep(tsTemp,tsTempInitialTime,tsTempInitialStep,ierr);CHKERRQ(ierr)
      Call TSGetAdapt(tsTemp,tsAdaptTemp,ierr);CHKERRQ(ierr)
      Call TSAdaptSetFromOptions(tsAdaptTemp,ierr);CHKERRQ(ierr)
   End If

   !!! 
   !!! Allocate array of works and energies
   !!!
   Allocate(energy(size(MEF90DefMechCtx%CellSetOptionsBag)))
   energy = 0.0_Kr
   Allocate(work(size(MEF90DefMechCtx%CellSetOptionsBag)))
   work = 0.0_Kr

   !!!
   !!! Try to figure out if the file was formatted
   !!!
   If (MEF90Ctx%rank == 0) Then
      Call EXGVP(MEF90Ctx%fileExoUnit,"N",numfield,ierr)
   End If
   Call MPI_Bcast(numfield,1,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)   
   If (numfield == 0) Then
      Call MEF90DefMechFormatEXO(MEF90DefMechCtx,ierr)
      !!! Will have to figure out this one
   End If
   
   !!!
   !!! Actual computations / time stepping
   !!!
   If (MEF90DefMechGlobalOptions%mode == MEF90DefMech_ModeQuasiStatic) Then
      Do step = 1,MEF90GlobalOptions%timeNumStep
         Write(IOBuffer,100) step,time(step)
         Call PetscPrintf(MEF90Ctx%comm,IOBuffer,ierr);CHKERRQ(ierr)

         !!! Solve for temperature
         Select Case (MEF90HeatXferGlobalOptions%mode)
         Case (MEF90HeatXFer_ModeSteadyState) 
            !!! Update fields
            Call MEF90HeatXferSetTransients(MEF90HeatXferCtx,step,time(step),ierr)
            !!! Solve SNES
            Call MEF90HeatXferUpdateboundaryTemperature(MEF90HeatXferCtx%temperature,MEF90HeatXferCtx,ierr);
            Call SNESSolve(snesTemp,PETSC_NULL_OBJECT,MEF90HeatXferCtx%temperature,ierr);CHKERRQ(ierr)

            !!! Compute thermal energy
            Call MEF90HeatXFerEnergy(MEF90HeatXferCtx%temperature,time(step),MEF90HeatXferCtx,energy,work,ierr);CHKERRQ(ierr)
            Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
            Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
            Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
            Do set = 1, size(setID)
               Write(IOBuffer,101) setID(set),energy(set),work(set),energy(set)-work(set)
               Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
            End Do
            Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
            Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
            Write(IOBuffer,102) sum(energy),sum(work),sum(energy)-sum(work)
            Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
            !!! Save results
            Call MEF90HeatXferViewEXO(MEF90HeatXferCtx,step,ierr)
         Case (MEF90HeatXFer_ModeTransient)
            If (step > 1) Then
               !!! Update fields
               Call MEF90HeatXferSetTransients(MEF90HeatXferCtx,step,time(step),ierr)
               Call MEF90HeatXferUpdateboundaryTemperature(MEF90HeatXferCtx%temperature,MEF90HeatXferCtx,ierr);
               !!! Make sure TS does not overstep
               Call TSGetTime(tsTemp,t,ierr);CHKERRQ(ierr)
               If (t < time(step)) Then
                  Call TSAdaptSetStepLimits(tsAdaptTemp,PETSC_DECIDE,(time(step)-time)/2.0_Kr,ierr);CHKERRQ(ierr)
                  !!! Something is up here. 
                  !!! replacing the constant 10000 with a variable leads to divergence of TSAdapt
                  !!! when using gcc
                  Call TSSetDuration(tsTemp,10000,time(step),ierr);CHKERRQ(ierr)
                  Call TSSolve(tsTemp,MEF90HeatXferCtx%temperature,time(step),ierr);CHKERRQ(ierr)
                  Call TSGetTime(tsTemp,t,ierr);CHKERRQ(ierr)
                  time(step) = t
               Else
                  Write(IOBuffer,*) 'TS exceeded analysis time. Skipping step\n'
                  Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
               End If
            End If

            !!! Compute thermal energy
            Call MEF90HeatXFerEnergy(MEF90HeatXferCtx%temperature,time(step),MEF90HeatXferCtx,energy,work,ierr);CHKERRQ(ierr)
            Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
            Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
            Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
            Do set = 1, size(setID)
               Write(IOBuffer,101) setID(set),energy(set),work(set),energy(set)-work(set)
               Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
            End Do
            Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
            Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
            Write(IOBuffer,102) sum(energy),sum(work),sum(energy)-sum(work)
            Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
            !!! Save results
            Call MEF90HeatXferViewEXO(MEF90HeatXferCtx,step,ierr)
         End Select

         !!! Solve for displacement
         Select case(MEF90DefMechGlobalOptions%mode)
         Case (MEF90DefMech_ModeQuasiStatic)
            !!! Update fields
            Call MEF90DefMechSetTransients(MEF90DefMechCtx,step,time(step),ierr)
            Call MEF90DefMechUpdateboundaryDisplacement(MEF90DefMechCtx%displacement,MEF90DefMechCtx,ierr)

            !!! Solve SNES
            Call SNESSolve(snesDisp,PETSC_NULL_OBJECT,MEF90DefMechCtx%displacement,ierr);CHKERRQ(ierr)
            Call SNESGetConvergedReason(snesDisp,snesDispConvergedReason,ierr);CHKERRQ(ierr)
            Write(IOBuffer,*) "SNESConvergedReason returned ",snesDispConvergedReason,"\n"
            Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
         
            !!! Compute energies
            energy = 0.0_Kr
            work = 0.0_Kr
            Call MEF90DefMechWork(MEF90DefMechCtx%displacement,MEF90DefMechCtx,work,ierr);CHKERRQ(ierr)
            Call MEF90DefMechElasticEnergy(MEF90DefMechCtx%displacement,MEF90DefMechCtx,energy,ierr);CHKERRQ(ierr)
            Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
            Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
            Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
            Do set = 1, size(setID)
               Write(IOBuffer,201) setID(set),energy(set),work(set),energy(set)-work(set)
               Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
            End Do
            Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
            Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
            Write(IOBuffer,202) sum(energy),sum(work),sum(energy)-sum(work)
            Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
     
            !!! Save results and boundary Values
            Call MEF90DefMechViewEXO(MEF90DefMechCtx,step,ierr)
         End Select
      End Do
   End If
100 Format("Solving steady state step ",I4,", t=",ES12.5,"\n")
101 Format("cell set ",I4," thermal energy: ",ES12.5," fluxes work: ",ES12.5," total: ",ES12.5,"\n")
102 Format("======= Total thermal energy: ",ES12.5," fluxes work: ",ES12.5," total: ",ES12.5,"\n")
201 Format("cell set ",I4," elastic energy: ",ES12.5," work: ",ES12.5," total: ",ES12.5,"\n")
202 Format("======= Total elastic energy: ",ES12.5," work: ",ES12.5," total: ",ES12.5,"\n")

   !!! Clean up and exit nicely
   Select case(MEF90DefMechGlobalOptions%mode)
   Case (MEF90DefMech_ModeQuasiStatic)
      Call SNESDestroy(snesDisp,ierr);CHKERRQ(ierr)
      Call VecDestroy(residualDisp,ierr);CHKERRQ(ierr)
   End Select
   
   Select Case (MEF90HeatXferGlobalOptions%mode)
   Case (MEF90HeatXFer_ModeSteadyState) 
      Call SNESDestroy(snesTemp,ierr);CHKERRQ(ierr)
   Case (MEF90HeatXFer_ModeTransient) 
      Call TSDestroy(tsTemp,ierr);CHKERRQ(ierr)
   End Select

   Call MEF90DefMechCtxDestroyVectors(MEF90DefMechCtx,ierr)
   Nullify(MEF90HeatXferCtx%temperature)
   Call MEF90HeatXferCtxDestroyVectors(MEF90HeatXferCtx,ierr)
   Call VecDestroy(residualDisp,ierr);CHKERRQ(ierr)
   Call VecDestroy(residualTemp,ierr);CHKERRQ(ierr)

   Call DMDestroy(Mesh,ierr);CHKERRQ(ierr)

   DeAllocate(time)
   DeAllocate(energy)
   DeAllocate(work)
   Call MEF90DefMechCtxDestroy(MEF90DefMechCtx,ierr);CHKERRQ(ierr)
   Call MEF90HeatXferCtxDestroy(MEF90HeatXferCtx,ierr);CHKERRQ(ierr)
   Call MEF90CtxCloseEXO(MEF90Ctx,ierr)

   Call PetscViewerASCIIOpen(MEF90Ctx%comm,trim(MEF90Ctx%prefix)//'.log',logViewer, ierr);CHKERRQ(ierr)
   Call PetscLogView(logViewer,ierr);CHKERRQ(ierr)
   Call PetscViewerDestroy(logViewer,ierr);CHKERRQ(ierr)
   Call MEF90CtxDestroy(MEF90Ctx,ierr)
   Call MEF90Finalize(ierr)
   Call PetscFinalize(ierr)
End Program ThermoElasticity
