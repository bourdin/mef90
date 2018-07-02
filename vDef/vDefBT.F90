#include "../MEF90/mef90.inc"
Program vDef
#include <finclude/petscdef.h>
   Use petsc
   Use m_MEF90
   !Use m_vDef
   Use m_MEF90_DefMechCtx
   Use m_MEF90_DefMech
   Use m_MEF90_HeatXferCtx
   Use m_MEF90_HeatXfer
   use m_vDefDefault
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

   !Type(SNES)                                         :: snesDisp
   SNESConvergedReason                                :: snesDispConvergedReason
   Type(Vec)                                          :: residualDisp
   Type(Vec)                                          :: damageAltMinOld,displacementAltMinOld
   Type(Vec)                                          :: damageLB,damageUB   
   PetscReal,Dimension(:),Pointer                     :: damageArray,damageAltMinOldArray,damageLBArray,damageUBArray
   PetscInt                                           :: iDof
   PetscReal                                          :: SOROmega,mySOROmega
   Type(Vec)                                          :: plasticStrainOld,plasticStrainPrevious
   !Type(SNES)                                         :: snesDamage
   SNESConvergedReason                                :: snesDamageConvergedReason
   Type(Vec)                                          :: residualDamage,localVec
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
   
   Integer                                            :: step,stepold
   PetscInt                                           :: dim
   PetscInt                                           :: AltMinStep = 0
   
   PetscReal                                          :: alphaMaxChange,alphaMin,alphaMax
   
   PetscBool                                          :: BTActive = Petsc_False
   PetscInt                                           :: BTStep,BTminStep,BTMaxSTep,BTDirection

   Type(Vec)                                          :: cumulatedDissipatedPlasticEnergyOld
   Type(Vec)                                          :: cumulatedDissipatedPlasticEnergyVariation

      


   !!! Initialize MEF90
   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   Call MEF90Initialize(ierr)
   Call PetscPrintf(PETSC_COMM_WORLD," # vDefBT: numerical implementation of variational models of Defect Mechanics\n",ierr);CHKERRQ(ierr)
   
   
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
   Call VecDuplicate(MEF90DefMechCtx%displacement,displacementAltMinOld,ierr);CHKERRQ(ierr)
   Call VecDuplicate(MEF90DefMechCtx%displacement,residualDisp,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(residualDisp,"residualDisp",ierr);CHKERRQ(ierr)
   Call MEF90DefMechCreateSNESDisplacement(MEF90DefMechCtx,MEF90DefMechCtx%snesDisp,residualDisp,ierr)
   Call VecDuplicate(MEF90DefMechCtx%damage,damageAltMinOld,ierr);CHKERRQ(ierr)
   Call VecDuplicate(MEF90DefMechCtx%damage,residualDamage,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(residualDamage,"residualDamage",ierr);CHKERRQ(ierr)
   Call MEF90DefMechCreateSNESDamage(MEF90DefMechCtx,MEF90DefMechCtx%snesDamage,residualDamage,ierr)
   DeAllocate(MEF90DefMechCtx%temperature)
   
   !!!cumulatedDissipatedPlasticEnergy Vectors
   Call VecDuplicate(MEF90DefMechCtx%cumulatedDissipatedPlasticEnergy,cumulatedDissipatedPlasticEnergyOld,ierr);CHKERRQ(ierr)
   Call VecDuplicate(MEF90DefMechCtx%cumulatedDissipatedPlasticEnergy,cumulatedDissipatedPlasticEnergyVariation,ierr);CHKERRQ(ierr)
   !Call VecCopy(MEF90DefMechCtx%cumulatedDissipatedPlasticEnergy,cumulatedDissipatedPlasticEnergyOld,ierr);CHKERRQ(ierr)

   
   Call VecDuplicate(MEF90DefMechCtx%plasticStrain,plasticStrainOld,ierr);CHKERRQ(ierr)
   Call VecDuplicate(MEF90DefMechCtx%PlasticStrain,plasticStrainPrevious,ierr);CHKERRQ(ierr)
   
   !!! Create sections, vectors, and solvers for HeatXfer Context
   If (MEF90HeatXferGlobalOptions%timeSteppingType /= MEF90HeatXfer_timeSteppingTypeNULL) Then
      Call MEF90HeatXferCtxSetSections(MEF90HeatXferCtx,ierr)
      Call MEF90HeatXferCtxCreateVectors(MEF90HeatXferCtx,ierr)
      Call VecDuplicate(MEF90HeatXferCtx%temperature,residualTemp,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(residualTemp,"residualTemp",ierr);CHKERRQ(ierr)
      Select Case(MEF90HeatXferGlobalOptions%timeSteppingType)
      Case (MEF90HeatXfer_timeSteppingTypeSteadyState)
         Call MEF90HeatXferCreateSNES(MEF90HeatXferCtx,snesTemp,residualTemp,ierr)
      Case (MEF90HeatXfer_timeSteppingTypeTransient)
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
   
   !!!
   !!! Try to figure out if the file was formatted
   !!!
   If (MEF90Ctx%rank == 0) Then
      Call EXGVP(MEF90Ctx%fileExoUnit,"N",numfield,ierr)
   End If
   Call MPI_Bcast(numfield,1,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)   
   If (numfield == 0) Then
      Call MEF90DefMechFormatEXO(MEF90DefMechCtx,time,ierr)
      !!! Will have to figure out this one
   End If
   
   !
   !!! Actual computations / time stepping
   !!!
   If (.NOT. MEF90GlobalOptions%dryrun) Then
      step = MEF90GlobalOptions%timeSkip+1
      stepold = step+1
      mainloopQS: Do
         BTActive = PETSC_FALSE
         !!! Solve for temperature
         Select Case (MEF90HeatXferGlobalOptions%timeSteppingType)
         Case (MEF90HeatXfer_timeSteppingTypeSteadyState) 
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
         Case (MEF90HeatXfer_timeSteppingTypeTransient)
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
         Case (MEF90HeatXfer_timeSteppingTypeNULL)
            Continue
         Case default
            Write(IOBuffer,*) "Implemented HeatXfer mode: ", MEF90HeatXferGlobalOptions%timeSteppingType, "\n"
            Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
            STOP
         End Select
         
         !!! Solve for displacement and damage
         Select case(MEF90DefMechGlobalOptions%timeSteppingType)
         Case (MEF90DefMech_timeSTeppingTypeQuasiStatic)
            Write(IOBuffer,200) step,time(step) 
            Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
            damageMaxChange = 1.0D+20
            If (stepOld > step) Then
               If (step == 1) Then
                  call VecSet(damageAltMinOld,0.0_Kr,ierr);CHKERRQ(ierr)
               Else
                  !!! We are going back to this step from the future.
                  !!! We need to reload damageAltMinOld in order to recompute bounds properly
                  Call DMGetLocalVector(MEF90DefMechCtx%DMScal,localVec,ierr);CHKERRQ(ierr)
                  Call VecLoadExodusVertex(MEF90DefMechCtx%DMScal,localVec,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
                                           MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step-1,MEF90DefMechGlobalOptions%damageOffset,ierr);CHKERRQ(ierr)
                  Call DMLocalToGlobalBegin(MEF90DefMechCtx%DMScal,localVec,INSERT_VALUES,damageAltMinOld,ierr);CHKERRQ(ierr)
                  Call DMLocalToGlobalEnd(MEF90DefMechCtx%DMScal,localVec,INSERT_VALUES,damageAltMinOld,ierr);CHKERRQ(ierr)
                  Call DMRestoreLocalVector(MEF90DefMechCtx%DMScal,localVec,ierr);CHKERRQ(ierr)
                  Call MEF90DefMechUpdateDamageBounds(MEF90DefMechCtx,MEF90DefMechCtx%snesDamage,damageAltMinOld,ierr);CHKERRQ(ierr)
               End If
            Else
               Call MEF90DefMechUpdateDamageBounds(MEF90DefMechCtx,MEF90DefMechCtx%snesDamage,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)
            EndIf

            !!! Update fields
            Call MEF90DefMechSetTransients(MEF90DefMechCtx,step,time(step),ierr)
            Call MEF90DefMechUpdateboundaryDisplacement(MEF90DefMechCtx%displacement,MEF90DefMechCtx,ierr)
            Call MEF90DefMechUpdateboundaryDamage(MEF90DefMechCtx%damage,MEF90DefMechCtx,ierr)

            Call SNESSetLagPreconditioner(MEF90DefMechCtx%snesDamage,1,ierr);CHKERRQ(ierr)
            Call SNESSetLagPreconditioner(MEF90DefMechCtx%snesDisp,1,ierr);CHKERRQ(ierr)
            AltMin: Do AltMinIter = 1, MEF90DefMechGlobalOptions%maxit
               AltMinStep = altminstep + 1
               Write(IObuffer,208) AltMinIter
               Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)

               !!! Solve SNES
               If (mod(AltMinIter-1,MEF90DefMechGlobalOptions%PCLag) == 0) Then
                  Call SNESSetLagPreconditioner(MEF90DefMechCtx%snesDisp,-2,ierr);CHKERRQ(ierr)
                  Call SNESSetLagPreconditioner(MEF90DefMechCtx%snesDamage,-2,ierr);CHKERRQ(ierr)
               End If 

               Call VecCopy(MEF90DefMechCtx%displacement,displacementAltMinOld,ierr);CHKERRQ(ierr)
               Call SNESSolve(MEF90DefMechCtx%snesDisp,PETSC_NULL_OBJECT,MEF90DefMechCtx%displacement,ierr);CHKERRQ(ierr)
               !!! Over relaxation on u, skipping the first alternate minimization step
               If (AltMinIter > 1) Then
                  SOROmega = abs(MEF90DefMechGlobalOptions%SOROmega)
                  Call VecAXPBY(MEF90DefMechCtx%displacement,1.0_Kr - SOROmega,SOROmega,displacementAltMinOld,ierr);CHKERRQ(ierr)
               EndIf

               Call SNESGetConvergedReason(MEF90DefMechCtx%snesDisp,snesDispConvergedReason,ierr);CHKERRQ(ierr)
               If (snesDispConvergedReason < 0) Then  
                  Write(IOBuffer,400) "displacement",snesDispConvergedReason
                  Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
               End If
               Call VecCopy(MEF90DefMechCtx%damage,damageAltMinOld,ierr);CHKERRQ(ierr)
               Call SNESSolve(MEF90DefMechCtx%snesDamage,PETSC_NULL_OBJECT,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)
               Call SNESGetConvergedReason(MEF90DefMechCtx%snesDamage,snesDamageConvergedReason,ierr);CHKERRQ(ierr)
               If (snesDamageConvergedReason < 0) Then
                  Write(IOBuffer,400) "damage field",snesDamageConvergedReason
                  Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
               End If

               !!! Over relaxation of the damage variable
               If (AltMinIter > 1) Then
                  If (MEF90DefMechGlobalOptions%SOROmega > 1.0_Kr) Then
                     mySOROmega = MEF90DefMechGlobalOptions%SOROmega
                     !!! LIMITED SOR
                     Call SNESVIGetVariableBounds(MEF90DefMechCtx%snesDamage,damageLB,damageUB,ierr);CHKERRQ(ierr)
                     Call VecGetArrayF90(damageLB,damageLBArray,ierr);CHKERRQ(ierr)
                     Call VecGetArrayF90(damageUB,damageUBArray,ierr);CHKERRQ(ierr)
                     Call VecGetArrayF90(damageAltMinOld,damageAltMinOldArray,ierr);CHKERRQ(ierr)
                     Call VecGetArrayF90(MEF90DefMechCtx%damage,damageArray,ierr);CHKERRQ(ierr)
                     Do iDof = 1, size(damageArray)
                        If (damageArray(iDof) > damageAltMinOldArray(iDof)) Then
                           mySOROmega = min(mySOROmega,(damageUBArray(iDof)-damageAltMinOldArray(iDof)) / (damageArray(iDof) - damageAltMinOldArray(iDof)))
                        Else If (damageArray(iDof) < damageAltMinOldArray(iDof)) Then
                           mySOROmega = min(mySOROmega,(damageLBArray(iDof)-damageAltMinOldArray(iDof)) / (damageArray(iDof) - damageAltMinOldArray(iDof)))
                        End If
                     End Do
                     Call MPI_AllReduce(mySOROmega,SOROmega,1,MPIU_SCALAR,MPI_MIN,PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)
                     Call VecRestoreArrayF90(MEF90DefMechCtx%damage,damageArray,ierr);CHKERRQ(ierr)
                     Call VecRestoreArrayF90(damageAltMinOld,damageAltMinOldArray,ierr);CHKERRQ(ierr)
                     Call VecRestoreArrayF90(damageUB,damageUBArray,ierr);CHKERRQ(ierr)
                     Call VecRestoreArrayF90(damageLB,damageLBArray,ierr);CHKERRQ(ierr)
                     Call VecAXPBY(MEF90DefMechCtx%damage,1.0_Kr - SOROmega,SOROmega,damageAltMinOld,ierr);CHKERRQ(ierr)
                  Else If (MEF90DefMechGlobalOptions%SOROmega < -1.0_Kr) Then
                     !!! PROJECTED SOR
                     SOROmega = -MEF90DefMechGlobalOptions%SOROmega
                     Call VecAXPBY(MEF90DefMechCtx%damage,1.0_Kr - SOROmega,SOROmega,damageAltMinOld,ierr);CHKERRQ(ierr)
                     Call SNESVIGetVariableBounds(MEF90DefMechCtx%snesDamage,damageLB,damageUB,ierr);CHKERRQ(ierr)
                     Call VecPointwiseMax(MEF90DefMechCtx%damage,MEF90DefMechCtx%damage,damageLB,ierr);CHKERRQ(ierr)
                     Call VecPointwiseMin(MEF90DefMechCtx%damage,MEF90DefMechCtx%damage,damageUB,ierr);CHKERRQ(ierr)
                  EndIf
               End If

               Call VecMin(MEF90DefMechCtx%damage,PETSC_NULL_INTEGER,alphaMin,ierr);CHKERRQ(ierr)
               Call VecMax(MEF90DefMechCtx%damage,PETSC_NULL_INTEGER,alphaMax,ierr);CHKERRQ(ierr)
               Call VecAxPy(damageAltMinOld,-1.0_Kr,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)
               ! damageAltMinOld = damageAltMinOld - damage
               Call VecNorm(damageAltMinOld,NORM_INFINITY,damageMaxChange,ierr);CHKERRQ(ierr)
               Write(IOBuffer,209) alphamin,alphamax,damageMaxChange
               Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)

               Call VecCopy(MEF90DefMechCtx%PlasticStrain,plasticStrainPrevious,ierr);CHKERRQ(ierr)
               Call MEF90DefMechPlasticStrainUpdate(MEF90DefMechCtx,MEF90DefMechCtx%PlasticStrain,MEF90DefMechCtx%displacement,PlasticStrainOld,plasticStrainPrevious,cumulatedDissipatedPlasticEnergyVariation,cumulatedDissipatedPlasticEnergyOld,ierr);CHKERRQ(ierr)
               Call VecWAXPY(MEF90DefMechCtx%cumulatedDissipatedPlasticEnergy,1.0_Kr,cumulatedDissipatedPlasticEnergyOld,cumulatedDissipatedPlasticEnergyVariation,ierr);CHKERRQ(ierr)

               ! Check for BT if necessary
               BTCheck: If ((MEF90DefMechGlobalOptions%BTInterval > 0) .AND. &
                   (mod(AltMinIter,MEF90DefMechGlobalOptions%BTInterval) == 0) .AND. &
                   (MEF90DefMechGlobalOptions%BTType /= MEF90DefMech_BTTypeNULL)) Then
                   !!!
                   !!! Recompute all energies
                   !!!
                  elasticEnergySet  = 0.0_Kr
                  forceWorkSet      = 0.0_Kr
                  surfaceEnergySet  = 0.0_Kr
                  cohesiveEnergySet = 0.0_Kr
                  Call MEF90DefMechElasticEnergy(MEF90DefMechCtx%displacement,MEF90DefMechCtx,elasticEnergySet,ierr);CHKERRQ(ierr)
                  Call MEF90DefMechWork(MEF90DefMechCtx%displacement,MEF90DefMechCtx,forceWorkSet,ierr);CHKERRQ(ierr)
                  Call MEF90DefMechSurfaceEnergy(MEF90DefMechCtx%damage,MEF90DefMechCtx,surfaceEnergySet,ierr);CHKERRQ(ierr)
                  Call MEF90DefMechCohesiveEnergy(MEF90DefMechCtx%displacement,MEF90DefMechCtx,cohesiveEnergySet,ierr);CHKERRQ(ierr)
                  elasticEnergy(step)  = sum(elasticEnergySet)
                  forceWork(step)      = sum(forceWorkSet)
                  surfaceEnergy(step)  = sum(surfaceEnergySet)
                  cohesiveEnergy(step) = sum(cohesiveEnergySet)
                  totalMechanicalEnergy(step) = elasticEnergy(step) - forceWork(step) + cohesiveEnergy(step) + surfaceEnergy(step)
                  
                  !!!
                  !!! Check for a BT step
                  !!!
                  If (MEF90DefMechGlobalOptions% BTType == MEF90DefMech_BTTypeForward) Then
                     If (MEF90DefMechGlobalOptions%BTScope > 0) Then
                        BTMinStep = max(1,step - MEF90DefMechGlobalOptions%BTScope)
                     Else
                         BTMinStep = 1
                     End If 
                     BTMaxStep   = step - 1
                     BTDirection = 1
                  Else
                     BTMinStep   = step - 1
                     If (MEF90DefMechGlobalOptions%BTScope > 0) Then
                        BTMaxStep = max(1,step - MEF90DefMechGlobalOptions%BTScope)
                     Else
                         BTMaxStep = step
                     End If 
                     BTDirection = -1
                  End If
                  Do BTStep = BTminStep,BTMaxSTep,BTDirection
                     If (time(step)**2 * (totalMechanicalEnergy(BTStep) - surfaceEnergy(step)) - time(BTStep)**2 * (elasticEnergy(step) - forceWork(step) + cohesiveEnergy(step)) & 
                         > time(step)**2 * abs(totalMechanicalEnergy(step)) * MEF90DefMechGlobalOptions%BTtol ) Then
                         !!! current solution is a better test field for step BTstep, backtracking
                        BTActive = PETSC_TRUE
                        Call PetscViewerASCIIPrintf(MEF90DefMechCtx%globalEnergyViewer,"\n\n",ierr);CHKERRQ(ierr)
                        Write(IOBuffer,450) BTStep
                        Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
                        EXIT
                     End If
                  End Do
               End If BTCheck
               If (BTActive) Then
                  EXIT
               End If
               If (damageMaxChange <= MEF90DefMechGlobalOptions%damageATol) Then
                  EXIT
               End If
               If (mod(AltMinIter,25) == 0) Then
                  Call MEF90DefMechViewEXO(MEF90DefMechCtx,step,ierr)
               End If
            End Do AltMin

            Call VecCopy(MEF90DefMechCtx%cumulatedDissipatedPlasticEnergy,cumulatedDissipatedPlasticEnergyOld,ierr);CHKERRQ(ierr)

            EndStep: If (.NOT. BTActive) Then
               !!! Compute energies
               elasticEnergySet  = 0.0_Kr
               forceWorkSet      = 0.0_Kr
               surfaceEnergySet  = 0.0_Kr
               cohesiveEnergySet = 0.0_Kr
               Call MEF90DefMechElasticEnergy(MEF90DefMechCtx%displacement,MEF90DefMechCtx,elasticEnergySet,ierr);CHKERRQ(ierr)
               Call MEF90DefMechWork(MEF90DefMechCtx%displacement,MEF90DefMechCtx,forceWorkSet,ierr);CHKERRQ(ierr)
               Call MEF90DefMechSurfaceEnergy(MEF90DefMechCtx%damage,MEF90DefMechCtx,surfaceEnergySet,ierr);CHKERRQ(ierr)
               Call MEF90DefMechCohesiveEnergy(MEF90DefMechCtx%displacement,MEF90DefMechCtx,cohesiveEnergySet,ierr);CHKERRQ(ierr)
               elasticEnergy(step)  = sum(elasticEnergySet)
               forceWork(step)      = sum(forceWorkSet)
               surfaceEnergy(step)  = sum(surfaceEnergySet)
               cohesiveEnergy(step) = sum(cohesiveEnergySet)
               totalMechanicalEnergy(step) = elasticEnergy(step) - forceWork(step) + cohesiveEnergy(step) + surfaceEnergy(step)
               !!!
               !!! Print and save energies
               !!!
               Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
               Call MEF90ISAllGatherMerge(MEF90Ctx%Comm,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
               Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
               Call PetscPrintf(MEF90Ctx%Comm,"\nMechanical energies: \n",ierr);CHKERRQ(ierr)
               Do set = 1, size(setID)
                  Write(IOBuffer,201) setID(set),elasticEnergySet(set),forceWorkSet(set),cohesiveEnergySet(set),surfaceEnergySet(set),elasticEnergySet(set) - forceWorkSet(set) + cohesiveEnergySet(set) + surfaceEnergySet(set)
                  Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
                  Write(IOBuffer,500) step,time(step),elasticEnergySet(set),forceWorkSet(set),cohesiveEnergySet(set),surfaceEnergySet(set),elasticEnergySet(set) - forceWorkSet(set) + cohesiveEnergySet(set) + surfaceEnergySet(set)
                  Call PetscViewerASCIIPrintf(MEF90DefMechCtx%setEnergyViewer(set),IOBuffer,ierr);CHKERRQ(ierr)
                  Call PetscViewerFlush(MEF90DefMechCtx%setEnergyViewer(set),ierr);CHKERRQ(ierr)
               End Do
               Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
               Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
               Write(IOBuffer,202) elasticEnergy(step),forceWork(step),cohesiveEnergy(step),surfaceEnergy(step),totalMechanicalEnergy(step)
               Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)

               Write(IOBuffer,500) step,time(step),elasticEnergy(step),forceWork(step),cohesiveEnergy(step),surfaceEnergy(step),totalMechanicalEnergy(step)
               Call PetscViewerASCIIPrintf(MEF90DefMechCtx%globalEnergyViewer,IOBuffer,ierr);CHKERRQ(ierr)
               Call PetscViewerFlush(MEF90DefMechCtx%globalEnergyViewer,ierr);CHKERRQ(ierr)

               !!!
               !!! Check for a BT step
               !!!
               BTCheck2: If ((MEF90DefMechGlobalOptions%BTInterval == 0) .AND. &
                   (MEF90DefMechGlobalOptions%BTType /= MEF90DefMech_BTTypeNULL)) Then
                  !!!
                  !!! Compute a BT step
                  !!!
                  If (MEF90DefMechGlobalOptions% BTType == MEF90DefMech_BTTypeForward) Then
                     If (MEF90DefMechGlobalOptions%BTScope > 0) Then
                        BTMinStep = max(1,step - MEF90DefMechGlobalOptions%BTScope)
                     Else
                         BTMinStep = 1
                     End If 
                     BTMaxStep   = step - 1
                     BTDirection = 1
                  Else
                     BTMinStep   = step - 1
                     BTMaxStep   = max(1,step - MEF90DefMechGlobalOptions%BTScope)
                     If (MEF90DefMechGlobalOptions%BTScope > 0) Then
                        BTMaxStep = max(1,step - MEF90DefMechGlobalOptions%BTScope)
                     Else
                         BTMaxStep = step
                     End If 
                     BTDirection = -1
                  End If
                  Do BTStep = BTminStep,BTMaxSTep,BTDirection
                     If (time(step)**2 * (totalMechanicalEnergy(BTStep) - surfaceEnergy(step)) - time(BTStep)**2 * (elasticEnergy(step) - forceWork(step) + cohesiveEnergy(step)) & 
                         > time(step)**2 * abs(totalMechanicalEnergy(step)) * MEF90DefMechGlobalOptions%BTtol ) Then
                         !!! current solution is a better test field for step BTstep, backtracking
                        BTActive = PETSC_TRUE
                        Call PetscViewerASCIIPrintf(MEF90DefMechCtx%globalEnergyViewer,"\n\n",ierr);CHKERRQ(ierr)
                        Write(IOBuffer,450) BTStep
                        Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
                        EXIT
                     End If
                  End Do
               End If BTCheck2
            End If EndStep
         Case (MEF90DefMech_timeSteppingTypeNULL)
            Continue
         Case default
            Write(IOBuffer,*) "Unimplemented DefMech time stepping type: ", MEF90DefMechGlobalOptions%timeSteppingType, "\n"
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
         !!!
         !!! Save performance log file
         !!!
         Call PetscViewerASCIIOpen(MEF90Ctx%comm,trim(MEF90FilePrefix(MEF90Ctx%resultFile))//'.log',logViewer, ierr);CHKERRQ(ierr)
         Call PetscLogView(logViewer,ierr);CHKERRQ(ierr)
         Call PetscViewerFlush(logViewer,ierr);CHKERRQ(ierr)
         !Call PetscViewerDestroy(logViewer,ierr);CHKERRQ(ierr)

         If (step == MEF90GlobalOptions%timeNumStep) Then
            EXIT
         ElseIf (BTActive) Then
            stepOld = step
            step = BTStep
            BTActive = PETSC_FALSE
         Else
            stepOld = step
            step = step + 1
         End If
      End Do MainloopQS
   End If
   !!! Clean up and exit nicely
   Select case(MEF90DefMechGlobalOptions%timeSteppingType)
   Case (MEF90DefMech_timeSteppingTypeQuasiStatic)
      Call SNESDestroy(MEF90DefMechCtx%snesDisp,ierr);CHKERRQ(ierr)
      Call SNESDestroy(MEF90DefMechCtx%snesDamage,ierr);CHKERRQ(ierr)
      Call VecDestroy(residualDisp,ierr);CHKERRQ(ierr)
      Call VecDestroy(residualDamage,ierr);CHKERRQ(ierr)
   End Select
   
   Select Case (MEF90HeatXferGlobalOptions%timeSteppingType)
   Case (MEF90HeatXfer_timeSteppingTypeSteadyState) 
      Call SNESDestroy(snesTemp,ierr);CHKERRQ(ierr)
      call VecDestroy(residualTemp,ierr);CHKERRQ(ierr)
   Case (MEF90HeatXfer_timeSteppingTypeTransient) 
      Call TSDestroy(tsTemp,ierr);CHKERRQ(ierr)
   End Select

   Write(IOBuffer,*) 'Total number of alternate minimizations:',AltMinStep,'\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)

   Call MEF90DefMechCtxDestroyVectors(MEF90DefMechCtx,ierr)
   Nullify(MEF90HeatXferCtx%temperature)
   Call MEF90HeatXferCtxDestroyVectors(MEF90HeatXferCtx,ierr)
   Call VecDestroy(damageAltMinOld,ierr);CHKERRQ(ierr)
   Call VecDestroy(displacementAltMinOld,ierr);CHKERRQ(ierr)

   Call VecDestroy(plasticStrainOld,ierr);CHKERRQ(ierr)
   Call VecDestroy(plasticStrainPrevious,ierr);CHKERRQ(ierr)
   Call VecDestroy(cumulatedDissipatedPlasticEnergyOld,ierr);CHKERRQ(ierr)
   Call VecDestroy(cumulatedDissipatedPlasticEnergyOld,ierr);CHKERRQ(ierr)
   Call VecDestroy(cumulatedDissipatedPlasticEnergyVariation,ierr);CHKERRQ(ierr)


   DeAllocate(elasticEnergySet)
   DeAllocate(surfaceEnergySet)
   DeAllocate(forceWorkSet)
   DeAllocate(cohesiveEnergySet)
   DeAllocate(thermalEnergySet)
   DeAllocate(heatFluxWorkSet)
   
   DeAllocate(elasticEnergy)
   DeAllocate(forceWork)
   DeAllocate(cohesiveEnergy)
   DeAllocate(surfaceEnergy)
   DeAllocate(totalMechanicalEnergy)

   DeAllocate(time)
   Call MEF90DefMechCtxDestroy(MEF90DefMechCtx,ierr);CHKERRQ(ierr)
   Call MEF90HeatXferCtxDestroy(MEF90HeatXferCtx,ierr);CHKERRQ(ierr)
   Call MEF90CtxCloseEXO(MEF90Ctx,ierr)
   Call DMDestroy(Mesh,ierr);CHKERRQ(ierr)

   !!Call PetscViewerASCIIOpen(MEF90Ctx%comm,trim(MEF90Ctx%prefix)//'.log',logViewer, ierr);CHKERRQ(ierr)
   If (.NOT. MEF90GlobalOptions%dryrun) Then
      Call PetscLogView(logViewer,ierr);CHKERRQ(ierr)
   End If
   Call PetscViewerDestroy(logViewer,ierr);CHKERRQ(ierr)
   Call MEF90CtxDestroy(MEF90Ctx,ierr)
   Call MEF90Finalize(ierr)
   Call PetscFinalize(ierr)
100 Format("\nHeat transfer: solving steady state step ",I4,", t=",ES12.5,"\n")
101 Format("cell set ",I4," thermal energy: ",ES12.5," fluxes work: ",ES12.5," total: ",ES12.5,"\n")
102 Format("======= Total thermal energy: ",ES12.5," fluxes work: ",ES12.5," total: ",ES12.5,"\n")
110 Format("\nHeat transfer: step ",I4,", until t=",ES12.5,"\n")
200 Format("\nMechanics: step ",I4,", t=",ES12.5,"\n")
201 Format("cell set ",I4,"  elastic energy: ",ES12.5," work: ",ES12.5," cohesive: ",ES12.5," surface: ",ES12.5," total: ",ES12.5,"\n")
202 Format("======= Total: elastic energy: ",ES12.5," work: ",ES12.5," cohesive: ",ES12.5," surface: ",ES12.5," total: ",ES12.5,"\n")
208 Format("   Alt. Min. step ",I5," ")
209 Format(" alpha min / max", ES12.5, " / ", ES12.5, ", max change ", ES12.5,"\n")
400 Format(" [ERROR]: ",A," SNESSolve failed with SNESConvergedReason ",I2,". \n Check http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESConvergedReason.html for error code meaning.\n")
410 Format(" [ERROR]: ",A," TSSolve failed with TSConvergedReason ",I2,". \n Check http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESConvergedReason.html for error code meaning.\n")
450 Format("BT: going back to step ",I4,"\n")
500 Format(I6, 6(ES16.5),"\n")

End Program vDef
