#include "../MEF90/mef90.inc"
Program CoupledPlasticityDamage
#include <finclude/petscdef.h>

   Use petsc
   Use m_MEF90
   Use m_MEF90_DefMechCtx
   Use m_MEF90_DefMech
   Use m_MEF90_HeatXferCtx
   Use m_MEF90_HeatXfer
   Implicit NONE   

   PetscErrorCode                                     :: ierr
   Type(MEF90Ctx_Type),target                         :: MEF90Ctx
   Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions

   !!! Defect mechanics contexts
   Type(MEF90DefMechCtx_Type)                         :: MEF90DefMechCtx
   Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions

   !!! HeatXfer contexts
   Type(MEF90HeatXferCtx_Type)                        :: MEF90HeatXferCtx
   Type(MEF90HeatXferGlobalOptions_Type),Pointer      :: MEF90HeatXferGlobalOptions

   
   Type(DM),target                                    :: Mesh
   Type(IS)                                           :: setIS,CellSetGlobalIS
   PetscInt,Dimension(:),Pointer                      :: setID
   PetscInt                                           :: numset,set
   PetscReal,Dimension(:),Pointer                     :: time
   PetscReal,Dimension(:),Pointer                     :: thermalEnergySet,heatFluxWorkSet
   PetscReal,Dimension(:),Pointer                     :: elasticEnergySet,surfaceEnergySet,forceWorkSet,cohesiveEnergySet
   PetscReal,Dimension(:),Pointer                     :: elasticEnergy,surfaceEnergy,forceWork,cohesiveEnergy,totalMechanicalEnergy

   Type(SNES)                                         :: snesDisp
   SNESConvergedReason                                :: snesDispConvergedReason
   Type(Vec)                                          :: residualDisp
   Type(SNES)                                         :: snesDamage
   SNESConvergedReason                                :: snesDamageConvergedReason
   Type(Vec)                                          :: residualDamage,damageOld
   PetscInt                                           :: AltMinIter
   PetscReal                                          :: damageMaxChange

   Type(SNES)                                         :: snesTemp
   SNESConvergedReason                                :: snesTempConvergedReason
   Type(TS)                                           :: tsTemp
   TSConvergedReason                                  :: tsTempConvergedReason
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
   
   PetscReal                                          :: alphaMaxChange,alphaMin,alphaMax
   
   PetscBool                                          :: BTActive = Petsc_False
   !PetscInt                                           :: BTStep,BTminStep,BTMaxSTep,BTDirection

   !!! plastic variables
   PetscReal,Dimension(:),Pointer                     :: plasticDissipationSet,plasticDissipationvariationSet
   PetscReal,Dimension(:),Pointer                     :: plasticDissipation,plasticDissipationvariation
   Type(Vec)                                          :: plasticStrainOld
   PetscInt                                           :: AltProjIter
   PetscReal                                          :: PlasticStrainMaxChange
   Type(Vec)                                          :: plasticstrainerror
   Type(Vec)                                          :: plasticStrainPrevious

   !!! WorkControl variables
   Type(Vec)                                          :: pressureForce_1
   PetscReal                                          :: pressure
   PetscReal                                          :: Work



!!! Default values of the contexts
   Type(MEF90CtxGlobalOptions_Type),Parameter         :: vDefDefaultGlobalOptions = MEF90CtxGlobalOptions_Type( &
                                                         1,                             & ! verbose
                                                         PETSC_FALSE,                   & ! validate
                                                         MEF90TimeInterpolation_linear, & ! timeInterpolation
                                                         0.0_Kr,                        & ! timeMin
                                                         1.0_Kr,                        & ! timeMax
                                                         11,                            & ! timeNumStep
                                                         MEF90FileFormat_EXOSingle)       ! fileFormat

   Type(MEF90DefMechGlobalOptions_Type),Parameter     :: vDefDefMechDefaultGlobalOptions2D = MEF90DefMechGlobalOptions_Type( &
                                                         MEF90DefMech_ModeQuasiStatic, & ! mode
                                                         PETSC_TRUE,              & ! disp_addNullSpace
                                                         3,                       & ! DisplacementOffset
                                                         2,                       & ! DamageOffset
                                                         3,                       & ! boundaryDisplacementOffset
                                                         0,                       & ! boundaryDamageOffset
                                                         1,                       & ! temperatureOffset
                                                         4,                       & ! ForceOffset
                                                         3,                       & ! pressureForceOffset
                                                         6,                       & ! plasticStrainOffset
                                                         6,                       & ! StressOffset
                                                         MEF90Scaling_Linear,     & ! boundaryDisplacementScaling
                                                         MEF90Scaling_CST,        & ! boundaryDamageScaling
                                                         MEF90Scaling_Linear,     & ! ForceScaling
                                                         MEF90Scaling_Linear,     & ! pressureForceScaling
                                                         1e-4,                    & ! damage_atol
                                                         1000,                    & ! maxit
                                                         10,                      & ! PCLag
                                                         0.,                      & ! irrevThres 
                                                         MEF90DefMech_BTTypeNULL, & ! BTType
                                                         -1,                      & ! BTInt
                                                         -1,                      & ! BTScope
                                                         1.0e-2,                  & ! BTTol
                                                         1.0e-4,                  & ! plasticStrainAtol
                                                         0)                         ! bloacknumberworkcontrolled
   Type(MEF90DefMechGlobalOptions_Type),Parameter     :: vDefDefMechDefaultGlobalOptions3D = MEF90DefMechGlobalOptions_Type( &
                                                         MEF90DefMech_ModeQuasiStatic, & ! mode
                                                         PETSC_TRUE,              & ! disp_addNullSpace
                                                         3,                       & ! DisplacementOffset
                                                         2,                       & ! DamageOffset
                                                         3,                       & ! boundaryDisplacementOffset
                                                         0,                       & ! boundaryDamageOffset
                                                         1,                       & ! temperatureOffset
                                                         4,                       & ! ForceOffset
                                                         3,                       & ! pressureForceOffset
                                                         0,                       & ! plasticStrainOffset
                                                         7,                       & ! StressOffset
                                                         MEF90Scaling_Linear,     & ! boundaryDisplacementScaling
                                                         MEF90Scaling_CST,        & ! boundaryDamageScaling
                                                         MEF90Scaling_Linear,     & ! ForceScaling
                                                         MEF90Scaling_Linear,     & ! pressureForceScaling
                                                         1e-4,                    & ! damage_atol
                                                         1000,                    & ! maxit
                                                         10,                      & ! PCLag
                                                         0.,                      & ! irrevThres 
                                                         MEF90DefMech_BTTypeNULL, & ! BTType
                                                         -1,                      & ! BTInt
                                                         -1,                      & ! BTScope
                                                         1.0e-2,                  & ! BTTol
                                                         1.0e-4,                  & ! plasticStrainAtol
                                                         0)                         ! bloacknumberworkcontrolled

   Type(MEF90DefMechCellSetOptions_Type),Parameter    :: vDefDefMechDefaultCellSetOptions = MEF90DefMechCellSetOptions_Type( &
                                                         -1,                                      & ! elemTypeShortIDDispl will be overriden
                                                         -1,                                      & ! elemTypeShortIDDamage will be overriden
                                                         [0.0_Kr,0.0_Kr,0.0_Kr],                  & ! force
                                                         0.0_Kr,                                  & ! pressureForce
                                                         MEF90DefMech_damageTypeAT1,              & ! damageType
                                                         MEF90DefMech_plasticityTypeNone,         & ! plasticityType
                                                         MEF90DefMech_unilateralContactTypeNone,  & ! unilateralContactType
                                                         [PETSC_FALSE,PETSC_FALSE,PETSC_FALSE],   & ! Has Displacement BC
                                                         [0.0_Kr,0.0_Kr,0.0_Kr],                  & ! boundary Displacement
                                                         PETSC_FALSE,                             & ! Has Damage BC
                                                         0._Kr)                                     ! Boundary Damage
   Type(MEF90DefMechVertexSetOptions_Type),Parameter  :: vDefDefMechDefaultVertexSetOptions = MEF90DefMechVertexSetOptions_Type( &
                                                         [PETSC_FALSE,PETSC_FALSE,PETSC_FALSE],   & ! Has Displacement BC
                                                         [0.0_Kr,0.0_Kr,0.0_Kr],                  & ! boundary Displacement
                                                         PETSC_FALSE,                             & ! Has Damage BC
                                                         0.0_Kr)                                    ! boundary Damage

   Type(MEF90HeatXferGlobalOptions_Type),Parameter    :: vDefHeatXferDefaultGlobalOptions = MEF90HeatXferGlobalOptions_Type( &
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
   Type(MEF90HeatXferCellSetOptions_Type),Parameter   :: vDefHeatXferDefaultCellSetOptions = MEF90HeatXferCellSetOptions_Type( &
                                                         -1,            & ! elemTypeShortID will be overriden
                                                         0.0_Kr,        & ! flux
                                                         0.0_Kr,        & ! surfaceThermalConductivity
                                                         0.0_Kr,        & ! externalTemp
                                                         PETSC_FALSE,   & ! Has BC
                                                         0.0_Kr)          ! boundaryTemp
                                                         
   Type(MEF90HeatXferVertexSetOptions_Type),Parameter :: vDefHeatXferDefaultVertexSetOptions = MEF90HeatXferVertexSetOptions_Type( &
                                                         PETSC_FALSE,   & ! Has BC
                                                         0.0_Kr)          ! boundaryTemp


   !!! Initialize MEF90
   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   Call MEF90Initialize(ierr)
   Call PetscPrintf(PETSC_COMM_WORLD," # vDucDef: numerical implementation of variational models of Ductile Defect Mechanics\n",ierr);CHKERRQ(ierr)
   
   
   !!! Get all MEF90-wide options
   Call MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,vDefDefaultGlobalOptions,ierr);CHKERRQ(ierr)
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
      Call MEF90DefMechCtxSetFromOptions(MEF90DefMechCtx,PETSC_NULL_CHARACTER,vDefDefMechDefaultGlobalOptions2D, &
                                         vDefDefMechDefaultCellSetOptions,vDefDefMechDefaultVertexSetOptions,ierr)
   Else
      Call MEF90DefMechCtxSetFromOptions(MEF90DefMechCtx,PETSC_NULL_CHARACTER,vDefDefMechDefaultGlobalOptions3D, &
                                         vDefDefMechDefaultCellSetOptions,vDefDefMechDefaultVertexSetOptions,ierr)
   End If
   Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)
   
   !!! Create HeatXfer context, get all HeatXfer options
   Call MEF90HeatXferCtxCreate(MEF90HeatXferCtx,Mesh,MEF90Ctx,ierr);CHKERRQ(ierr)
   Call MEF90HeatXferCtxSetFromOptions(MEF90HeatXferCtx,PETSC_NULL_CHARACTER,vDefHeatXferDefaultGlobalOptions, &
                                       vDefHeatXferDefaultCellSetOptions,vDefHeatXferDefaultVertexSetOptions,ierr)
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


   !!! Create sections, vectors, and solvers for DefMech Context
   Call MEF90DefMechCtxSetSections(MEF90DefMechCtx,ierr)
   Call MEF90DefMechCtxCreateVectors(MEF90DefMechCtx,ierr)
   Call VecDuplicate(MEF90DefMechCtx%damage,damageOld,ierr);CHKERRQ(ierr)
   Call VecDuplicate(MEF90DefMechCtx%displacement,residualDisp,ierr);CHKERRQ(ierr)

   Call PetscObjectSetName(residualDisp,"residualDisp",ierr);CHKERRQ(ierr)
   Call MEF90DefMechCreateSNESDisplacement(MEF90DefMechCtx,snesDisp,residualDisp,ierr)
   Call VecDuplicate(MEF90DefMechCtx%damage,residualDamage,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(residualDamage,"residualDamage",ierr);CHKERRQ(ierr)
   Call MEF90DefMechCreateSNESDamage(MEF90DefMechCtx,snesDamage,residualDamage,ierr)
   DeAllocate(MEF90DefMechCtx%temperature)
   
   Call VecDuplicate(MEF90DefMechCtx%plasticStrain,plasticStrainOld,ierr);CHKERRQ(ierr)
   Call VecDuplicate(MEF90DefMechCtx%plasticStrain,plasticStrainPrevious,ierr);CHKERRQ(ierr)

   
   !!! As long as plasticity is not implemented, there is no point in keeping the pastic strain around
   !!!DeAllocate(MEF90DefMechCtx%plasticStrain)
   Call VecDuplicate(MEF90DefMechCtx%plasticStrain,plasticStrainOld,ierr);CHKERRQ(ierr)
   Call VecDuplicate(MEF90DefMechCtx%plasticStrain,plasticStrainPrevious,ierr);CHKERRQ(ierr)

   !!! Create sections, vectors, and solvers for HeatXfer Context
   If (MEF90HeatXferGlobalOptions%mode /= MEF90HeatXfer_ModeNULL) Then
      Call MEF90HeatXferCtxSetSections(MEF90HeatXferCtx,ierr)
      Call MEF90HeatXferCtxCreateVectors(MEF90HeatXferCtx,ierr)
      Call VecDuplicate(MEF90HeatXferCtx%temperature,residualTemp,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(residualTemp,"residualTemp",ierr);CHKERRQ(ierr)
      Select Case(MEF90HeatXferGlobalOptions%mode)
      Case (MEF90HeatXFer_ModeSteadyState)
         Call MEF90HeatXferCreateSNES(MEF90HeatXferCtx,snesTemp,residualTemp,ierr)
      Case (MEF90HeatXFer_ModeTransient)
         Call MEF90HeatXferCreateTS(MEF90HeatXferCtx,tsTemp,residualTemp,ierr)
         tsTempInitialStep = (time(size(time))-time(1)) / (size(time) + 0.0_Kr) / 10.0_Kr
         tsTempInitialTime = time(1)
         Call TSSetInitialTimeStep(tsTemp,tsTempInitialTime,tsTempInitialStep,ierr);CHKERRQ(ierr)
         Call TSGetAdapt(tsTemp,tsAdaptTemp,ierr);CHKERRQ(ierr)
         Call TSAdaptSetFromOptions(tsAdaptTemp,ierr);CHKERRQ(ierr)
      End Select

      !!! Link the temperature field in the DefMechContext with that of the HeatXfer
      MEF90DefMechCtx%temperature => MEF90HeatXferCtx%temperature
   End If   
   !!! 
   !!! Allocate array of works and energies
   !!!
   Allocate(elasticEnergySet(size(MEF90DefMechCtx%CellSetOptionsBag)))
   elasticEnergySet = 0.0_Kr
   Allocate(surfaceEnergySet(size(MEF90DefMechCtx%CellSetOptionsBag)))
   surfaceEnergySet = 0.0_Kr
   Allocate(forceWorkSet(size(MEF90DefMechCtx%CellSetOptionsBag)))
   forceWorkSet = 0.0_Kr
   Allocate(cohesiveEnergySet(size(MEF90DefMechCtx%CellSetOptionsBag)))
   cohesiveEnergySet = 0.0_Kr
   Allocate(thermalEnergySet(size(MEF90DefMechCtx%CellSetOptionsBag)))
   thermalEnergySet = 0.0_Kr
   Allocate(heatFluxWorkSet(size(MEF90DefMechCtx%CellSetOptionsBag)))
   heatFluxWorkSet = 0.0_Kr

   Allocate(plasticDissipationvariationSet(size(MEF90DefMechCtx%CellSetOptionsBag)))
   plasticDissipationvariationSet = 0.0_Kr
   Allocate(plasticDissipationSet(size(MEF90DefMechCtx%CellSetOptionsBag)))
   plasticDissipationSet = 0.0_Kr

   
   Allocate(elasticEnergy(MEF90GlobalOptions%timeNumStep))
   elasticEnergy = 0.0_Kr
   Allocate(surfaceEnergy(MEF90GlobalOptions%timeNumStep))
   surfaceEnergy = 0.0_Kr
   Allocate(forceWork(MEF90GlobalOptions%timeNumStep))
   forceWork = 0.0_Kr
   Allocate(cohesiveEnergy(MEF90GlobalOptions%timeNumStep))
   cohesiveEnergy = 0.0_Kr
   Allocate(totalMechanicalEnergy(MEF90GlobalOptions%timeNumStep))
   totalMechanicalEnergy = 0.0_Kr

   Allocate(plasticDissipation(MEF90GlobalOptions%timeNumStep))
   plasticDissipation = 0.0_Kr

   Allocate(plasticDissipationvariation(MEF90GlobalOptions%timeNumStep))
   plasticDissipationvariation = 0.0_Kr
   
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
   
   !
   !!! Actual computations / time stepping
   !!!
   If (.NOT. MEF90GlobalOptions%dryrun) Then
      step = 1
      mainloopQS: Do
         BTActive = PETSC_FALSE
         !!! Solve for temperature
         Select Case (MEF90HeatXferGlobalOptions%mode)
         Case (MEF90HeatXFer_ModeSteadyState) 
            Write(IOBuffer,100) step,time(step)
            Call PetscPrintf(MEF90Ctx%comm,IOBuffer,ierr);CHKERRQ(ierr)

            !!! Update fields
            Call MEF90HeatXferSetTransients(MEF90HeatXferCtx,step,time(step),ierr)
            !!! Solve SNES
            Call MEF90HeatXferUpdateboundaryTemperature(MEF90HeatXferCtx%temperature,MEF90HeatXferCtx,ierr);
            Call SNESSolve(snesTemp,PETSC_NULL_OBJECT,MEF90HeatXferCtx%temperature,ierr);CHKERRQ(ierr)
            Call SNESGetConvergedReason(snesTemp,snesTempConvergedReason,ierr);CHKERRQ(ierr)
            If (snesTempConvergedReason < 0) Then  
               Write(IOBuffer,400) "temperature",snesTempConvergedReason
               Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
            End If

            !!! Compute thermal energy
            Call MEF90HeatXFerEnergy(MEF90HeatXferCtx%temperature,time(step),MEF90HeatXferCtx,thermalEnergySet,heatFluxWorkSet,ierr);CHKERRQ(ierr)
            Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
            Call MEF90ISAllGatherMerge(MEF90Ctx%Comm,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
            Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
            Call PetscPrintf(MEF90Ctx%Comm,"\nThermal energies: \n",ierr);CHKERRQ(ierr)
            Do set = 1, size(setID)
               Write(IOBuffer,101) setID(set),thermalEnergySet(set),heatFluxWorkSet(set),thermalEnergySet(set)-heatFluxWorkSet(set)
               Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
            End Do
            Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
            Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
            Write(IOBuffer,102) sum(thermalEnergySet),sum(heatFluxWorkSet),sum(thermalEnergySet)-sum(heatFluxWorkSet)
            Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)

            !!! Save results
            Call MEF90HeatXferViewEXO(MEF90HeatXferCtx,step,ierr)
         Case (MEF90HeatXFer_ModeTransient)
            If (step > 1) Then
               Write(IOBuffer,110) step,time(step)
               Call PetscPrintf(MEF90Ctx%comm,IOBuffer,ierr);CHKERRQ(ierr)
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
                  Call TSGetConvergedReason(tsTemp,tsTempConvergedReason,ierr);CHKERRQ(ierr)
                  If (tsTempConvergedReason < 0) Then  
                     Write(IOBuffer,410) "temperature",tsTempConvergedReason
                     Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
                  End If
                  Call TSGetTime(tsTemp,t,ierr);CHKERRQ(ierr)
                  time(step) = t
               Else
                  Write(IOBuffer,*) 'TS exceeded analysis time. Skipping step\n'
                  Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
               End If
            End If

            !!! Compute thermal energy
            Call MEF90HeatXFerEnergy(MEF90HeatXferCtx%temperature,time(step),MEF90HeatXferCtx,thermalEnergySet,heatFluxWorkSet,ierr);CHKERRQ(ierr)
            Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
            Call MEF90ISAllGatherMerge(MEF90Ctx%Comm,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
            Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
            Call PetscPrintf(MEF90Ctx%Comm,"\nThermal energies: \n",ierr);CHKERRQ(ierr)
            Do set = 1, size(setID)
               Write(IOBuffer,101) setID(set),thermalEnergySet(set),heatFluxWorkSet(set),thermalEnergySet(set)-heatFluxWorkSet(set)
               Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
            End Do
            Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
            Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
            Write(IOBuffer,102) sum(thermalEnergySet),sum(heatFluxWorkSet),sum(thermalEnergySet)-sum(heatFluxWorkSet)
            Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
            !!! Save results
            Call MEF90HeatXferViewEXO(MEF90HeatXferCtx,step,ierr)
         Case (MEF90HeatXfer_ModeNULL)
            Continue
         Case default
            Write(IOBuffer,*) "Implemented HeatXfer mode: ", MEF90HeatXferGlobalOptions%mode, "\n"
            Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
            STOP
         End Select
         
         !!! Solve for displacement and damage
         Select case(MEF90DefMechGlobalOptions%mode)
         Case (MEF90DefMech_ModeQuasiStatic)
            Write(IOBuffer,200) step,time(step) 
            Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
            damageMaxChange = 1.0D+20


            Call MEF90DefMechUpdateDamageBounds(MEF90DefMechCtx,snesDamage,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)


            !!! Update fields
            Call MEF90DefMechSetTransients(MEF90DefMechCtx,step,time(step),ierr)
            Call MEF90DefMechUpdateboundaryDisplacement(MEF90DefMechCtx%displacement,MEF90DefMechCtx,ierr)
            Call MEF90DefMechUpdateboundaryDamage(MEF90DefMechCtx%damage,MEF90DefMechCtx,ierr)

            Call SNESSetLagPreconditioner(snesDamage,1,ierr);CHKERRQ(ierr)
            Call SNESSetLagPreconditioner(snesDisp,1,ierr);CHKERRQ(ierr)
            


            !!! Save the pressure equal to one
            If (MEF90DefMechGlobalOptions%BlockNumberWorkControlled /= 0) Then 
               Call VecDuplicate(MEF90DefMechCtx%pressureForce,pressureForce_1,ierr);CHKERRQ(ierr)
               Call VecCopy(MEF90DefMechCtx%pressureForce,pressureForce_1,ierr);CHKERRQ(ierr)
            End If

            AltMin: Do AltMinIter = 1, MEF90DefMechGlobalOptions%maxit 
               Write(IObuffer,208) AltMinIter
               !Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
               AltProj: Do AltProjIter = 1, MEF90DefMechGlobalOptions%maxit 
                  Write(IObuffer,308) AltProjIter
                  !Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)

                  !!! Solve SNES

                  Call SNESSolve(snesDisp,PETSC_NULL_OBJECT,MEF90DefMechCtx%displacement,ierr);CHKERRQ(ierr)
                  Call SNESGetConvergedReason(snesDisp,snesDispConvergedReason,ierr);CHKERRQ(ierr)
                  If (snesDispConvergedReason < 0) Then
                     Write(IOBuffer,400) "displacement",snesDispConvergedReason
                     Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
                  End If


                  !!! Work Controled part 
                  If (MEF90DefMechGlobalOptions%BlockNumberWorkControlled /= 0) Then 
                     forceWorkSet      = 0.0_Kr
                     Call MEF90DefMechWork(MEF90DefMechCtx%displacement,MEF90DefMechCtx,forceWorkSet,ierr);CHKERRQ(ierr)
                     work = forceWorkSet(MEF90DefMechGlobalOptions%BlockNumberWorkControlled)
                     pressure = sqrt(time(step)/work)
                     Call VecScale(MEF90DefMechCtx%pressureForce,pressure)
                     Call VecScale(MEF90DefMechCtx%displacement,pressure)
                  Else
                     pressure=0
                  End If


                  !!! Solve PlasticProjection
                  !!! Save the PlasticStrainPrevious, add damage in DefMechPlasticStrainUpdate, error on the norm PlasticStrain

                  Call VecCopy(MEF90DefMechCtx%plasticStrain,plasticStrainPrevious,ierr);CHKERRQ(ierr)
                  Call MEF90DefMechPlasticStrainUpdate(MEF90DefMechCtx,MEF90DefMechCtx%plasticStrain,MEF90DefMechCtx%displacement,plasticStrainOld,plasticStrainPrevious,ierr);CHKERRQ(ierr)
                  Call VecAxPy(plasticStrainPrevious,-1.0_Kr,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)
                  Call VecNorm(plasticStrainPrevious,NORM_INFINITY,PlasticStrainMaxChange,ierr);CHKERRQ(ierr)


                  Write(IOBuffer,300) PlasticStrainMaxChange
                  Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)

                  If (PlasticStrainMaxChange <= MEF90DefMechGlobalOptions%plasticStrainATol) Then
                     EXIT
                  End If
               End Do AltProj

               Call VecCopy(MEF90DefMechCtx%damage,damageOld,ierr);CHKERRQ(ierr)
               Call SNESSolve(snesDamage,PETSC_NULL_OBJECT,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)
               Call SNESGetConvergedReason(snesDamage,snesDamageConvergedReason,ierr);CHKERRQ(ierr)
               If (snesDamageConvergedReason < 0) Then
                  Write(IOBuffer,400) "damage field",snesDamageConvergedReason
                  Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
               End If

               Call VecMin(MEF90DefMechCtx%damage,PETSC_NULL_INTEGER,alphaMin,ierr);CHKERRQ(ierr)
               Call VecMax(MEF90DefMechCtx%damage,PETSC_NULL_INTEGER,alphaMax,ierr);CHKERRQ(ierr)
               Call VecAxPy(damageOld,-1.0_Kr,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)
               Call VecNorm(damageOld,NORM_INFINITY,damageMaxChange,ierr);CHKERRQ(ierr)
               Write(IOBuffer,209) alphamin,alphamax,damageMaxChange
               Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)

               Call VecNorm(damageOld,NORM_INFINITY,damageMaxChange,ierr);CHKERRQ(ierr)


               !!! Evaluation of the Work pressure and compare with the Work input (time)
               If (MEF90DefMechGlobalOptions%BlockNumberWorkControlled /= 0) Then 
                  forceWorkSet      = 0.0_Kr
                  Call MEF90DefMechWork(MEF90DefMechCtx%displacement,MEF90DefMechCtx,forceWorkSet,ierr);CHKERRQ(ierr)
                  work = forceWorkSet(MEF90DefMechGlobalOptions%BlockNumberWorkControlled)
                  If (damageMaxChange <= MEF90DefMechGlobalOptions%damageATol  .and. abs(time(step)-work)<= 1e-4) Then
                     EXIT
                  End If
               Else
                  If (damageMaxChange <= MEF90DefMechGlobalOptions%damageATol) Then
                     EXIT
                  End If
               End If

               If (mod(AltMinIter,25) == 0) Then
                  Call MEF90DefMechViewEXO(MEF90DefMechCtx,step,ierr)
               End If
            End Do AltMin

            !!! Work controlled destroy Vec
            If (MEF90DefMechGlobalOptions%BlockNumberWorkControlled /= 0) Then 
               Call VecDestroy(pressureForce_1,ierr);CHKERRQ(ierr)
            End If


            !!! Compute energies
            elasticEnergySet  = 0.0_Kr
            forceWorkSet      = 0.0_Kr
            surfaceEnergySet  = 0.0_Kr
            cohesiveEnergySet = 0.0_Kr
            plasticDissipationvariationSet = 0.0_Kr


            Call MEF90DefMechElasticEnergy(MEF90DefMechCtx%displacement,MEF90DefMechCtx,elasticEnergySet,ierr);CHKERRQ(ierr)
            Call MEF90DefMechWork(MEF90DefMechCtx%displacement,MEF90DefMechCtx,forceWorkSet,ierr);CHKERRQ(ierr)
            Call MEF90DefMechSurfaceEnergy(MEF90DefMechCtx%damage,MEF90DefMechCtx,surfaceEnergySet,ierr);CHKERRQ(ierr)
            Call MEF90DefMechCohesiveEnergy(MEF90DefMechCtx%displacement,MEF90DefMechCtx,cohesiveEnergySet,ierr);CHKERRQ(ierr)
            Call MEF90DefMechPlasticDissipation(MEF90DefMechCtx%displacement,MEF90DefMechCtx,plasticStrainOld,plasticDissipationvariationSet,ierr);CHKERRQ(ierr)
            
            plasticDissipationvariation(step) = sum(plasticDissipationvariationSet)
            plasticDissipation(step) = plasticDissipation(step-1) + plasticDissipationvariation(step)

            elasticEnergy(step)  = sum(elasticEnergySet)
            forceWork(step)      = sum(forceWorkSet)
            surfaceEnergy(step)  = sum(surfaceEnergySet)
            cohesiveEnergy(step) = sum(cohesiveEnergySet)
            totalMechanicalEnergy(step) = elasticEnergy(step) - forceWork(step) + cohesiveEnergy(step) + surfaceEnergy(step) + plasticDissipation(step)
            !!!
            !!! Print and save energies
            !!!
            Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
            Call MEF90ISAllGatherMerge(MEF90Ctx%Comm,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
            Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
            Call PetscPrintf(MEF90Ctx%Comm,"\nMechanical energies: \n",ierr);CHKERRQ(ierr)
            Do set = 1, size(setID)

               plasticDissipationSet(set)= plasticDissipationSet(set) + plasticDissipationvariationSet(set)
               Write(IOBuffer,201) setID(set),elasticEnergySet(set),forceWorkSet(set),cohesiveEnergySet(set),surfaceEnergySet(set),elasticEnergySet(set) - forceWorkSet(set) + cohesiveEnergySet(set) + surfaceEnergySet(set) + plasticDissipationSet(set), plasticDissipationSet(set)
               Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)

               Write(IOBuffer,500) step,time(step),elasticEnergySet(set),forceWorkSet(set),cohesiveEnergySet(set),surfaceEnergySet(set),elasticEnergySet(set) - forceWorkSet(set) + cohesiveEnergySet(set) + surfaceEnergySet(set) + plasticDissipationSet(set) , plasticDissipationSet(set) 
               Call PetscViewerASCIIPrintf(MEF90DefMechCtx%setEnergyViewer(set),IOBuffer,ierr);CHKERRQ(ierr)
               Call PetscViewerFlush(MEF90DefMechCtx%setEnergyViewer(set),ierr);CHKERRQ(ierr)
            End Do
            Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
            Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
            Write(IOBuffer,202) elasticEnergy(step),forceWork(step),cohesiveEnergy(step),surfaceEnergy(step),totalMechanicalEnergy(step),plasticDissipation(step)
            Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
            Write(IOBuffer,500) step,time(step),elasticEnergy(step),cohesiveEnergy(step),forceWork(step),surfaceEnergy(step),totalMechanicalEnergy(step),plasticDissipation(step)
            Call PetscViewerASCIIPrintf(MEF90DefMechCtx%globalEnergyViewer,IOBuffer,ierr);CHKERRQ(ierr)
            Call PetscViewerFlush(MEF90DefMechCtx%globalEnergyViewer,ierr);CHKERRQ(ierr)



         Case (MEF90DefMech_ModeNULL)
            Continue
         Case default
            Write(IOBuffer,*) "Implemented DefMech mode: ", MEF90DefMechGlobalOptions%mode, "\n"
            Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
            STOP
         End Select
         !!!
         !!! Save results and boundary Values
         !!!
         If (MEF90DefMechGlobalOptions%stressOffset > 0) Then
            Call MEF90DefMechStress(MEF90DefMechCtx%displacement,MEF90DefMechCtx,MEF90DefMechCtx%stress,ierr)
         End If
         Call MEF90DefMechViewEXO(MEF90DefMechCtx,step,ierr)

         !!! Update plasticstrainold
         Call VecCopy(MEF90DefMechCtx%plasticStrain,plasticStrainOld,ierr);CHKERRQ(ierr)

         !!!
         !!! Save performance log file
         !!!
         Call PetscViewerASCIIOpen(MEF90Ctx%comm,trim(MEF90Ctx%prefix)//'.log',logViewer, ierr);CHKERRQ(ierr)
         Call PetscLogView(logViewer,ierr);CHKERRQ(ierr)
         Call PetscViewerDestroy(logViewer,ierr);CHKERRQ(ierr)

         If (step == MEF90GlobalOptions%timeNumStep) Then
            EXIT
         Else
            step = step + 1
         End If
      End Do MainloopQS
   End If
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
   Call VecDestroy(damageOld,ierr);CHKERRQ(ierr)
   Call VecDestroy(plasticStrainPrevious,ierr);CHKERRQ(ierr)

   Call DMDestroy(Mesh,ierr);CHKERRQ(ierr)

   DeAllocate(elasticEnergySet)
   DeAllocate(surfaceEnergySet)
   DeAllocate(forceWorkSet)
   DeAllocate(cohesiveEnergySet)
   DeAllocate(thermalEnergySet)
   DeAllocate(heatFluxWorkSet)
   DeAllocate(plasticDissipationSet)
   DeAllocate(plasticDissipationvariationSet)
   
   DeAllocate(elasticEnergy)
   DeAllocate(forceWork)
   DeAllocate(cohesiveEnergy)
   DeAllocate(surfaceEnergy)
   DeAllocate(totalMechanicalEnergy)
   DeAllocate(plasticDissipation)
   DeAllocate(plasticDissipationvariation)

   Call MEF90DefMechCtxDestroy(MEF90DefMechCtx,ierr);CHKERRQ(ierr)
   Call MEF90HeatXferCtxDestroy(MEF90HeatXferCtx,ierr);CHKERRQ(ierr)
   Call MEF90CtxCloseEXO(MEF90Ctx,ierr)

   Call PetscViewerASCIIOpen(MEF90Ctx%comm,trim(MEF90Ctx%prefix)//'.log',logViewer, ierr);CHKERRQ(ierr)
   Call PetscLogView(logViewer,ierr);CHKERRQ(ierr)
   Call PetscViewerDestroy(logViewer,ierr);CHKERRQ(ierr)
   Call MEF90CtxDestroy(MEF90Ctx,ierr)
   Call MEF90Finalize(ierr)
   Call PetscFinalize(ierr)
100 Format("\nHeat transfer: solving steady state step ",I4,", t=",ES12.5,"\n")
101 Format("cell set ",I4," thermal energy: ",ES12.5," fluxes work: ",ES12.5," total: ",ES12.5,"\n")
102 Format("======= Total thermal energy: ",ES12.5," fluxes work: ",ES12.5," total: ",ES12.5,"\n")
110 Format("\nHeat transfer: step ",I4,", until t=",ES12.5,"\n")
200 Format("\nMechanics: step ",I4,", t=",ES12.5,"\n")
201 Format("cell set ",I4,"  elastic energy: ",ES12.5," work: ",ES12.5," cohesive: ",ES12.5," surface: ",ES12.5," total: ",ES12.5," plastic dissipation: ",ES12.5,"\n")
202 Format("======= Total: elastic energy: ",ES12.5," work: ",ES12.5," cohesive: ",ES12.5," surface: ",ES12.5," total: ",ES12.5, " plastic dissipation: ",ES12.5,"\n")
208 Format("   Alt. Min. step ",I5," ")
209 Format(" alpha min / max", ES12.5, " / ", ES12.5, ", max change ", ES12.5,"\n")
300 Format("                 plastic strain max change: ", ES12.5, "\n")
308 Format("   Alt. Proj. step ",I5," ")
400 Format(" [ERROR]: ",A," SNESSolve failed with SNESConvergedReason ",I2,". \n Check http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESConvergedReason.html for error code meaning.\n")
410 Format(" [ERROR]: ",A," TSSolve failed with TSConvergedReason ",I2,". \n Check http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESConvergedReason.html for error code meaning.\n")
500 Format(I6, 7(ES16.5), "\n")
End Program CoupledPlasticityDamage
