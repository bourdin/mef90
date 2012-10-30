#include "../MEF90/mef90.inc"
Program TestHeatXfer
#include <finclude/petscdef.h>
   Use m_MEF90
   Use m_MEF90_HeatXfer
   Use m_MEF90_HeatXferAssembly2D
   Use petsc
   Implicit NONE   

   PetscErrorCode                                     :: ierr
   Type(MEF90HeatXferCtx_Type)                        :: MEF90HeatXferCtx
   Type(MEF90HeatXferGlobalOptions_Type),Parameter    :: MEF90HeatXferDefaultGlobalOptions = MEF90HeatXferGlobalOptions_Type( &
                                                         MEF90HeatXFer_ModeSteadyState, & ! mode
                                                         PETSC_FALSE,         & ! addNullSpace
                                                         1,                   & ! tempOffset
                                                         MEF90Scaling_Linear, & ! boundaryTempScaling
                                                         2,                   & ! boundaryTempOffset
                                                         MEF90Scaling_Linear, & ! externalTempScaling
                                                         1,                   & ! externalTempOffset
                                                         MEF90Scaling_CST, & ! fluxScaling
                                                         2)                     ! fluxOffset
   Type(MEF90HeatXferCellSetOptions_Type),Parameter   :: MEF90HeatXferDefaultCellSetOptions = MEF90HeatXferCellSetOptions_Type( &
                                                         -1,            & ! elemTypeShortID will be overriden
                                                         0.0_Kr,        & ! flux
                                                         0.0_Kr,        & ! surfaceThermalConductivity
                                                         0.0_Kr)          ! externalTemp
   Type(MEF90HeatXferVertexSetOptions_Type),Parameter :: MEF90HeatXferDefaultVertexSetOptions = MEF90HeatXferVertexSetOptions_Type( &
                                                         PETSC_TRUE,    & ! Has BC
                                                         0.0_Kr)          ! boundaryTemp
   Type(MEF90HeatXferGlobalOptions_Type),pointer      :: MEF90HeatXferGlobalOptions
                                                         
   Type(MEF90Ctx_Type),target                         :: MEF90Ctx
   Type(MEF90CtxGlobalOptions_Type),Parameter         :: MEF90DefaultGlobalOptions = MEF90CtxGlobalOptions_Type( &
                                                         1,             & ! verbose
                                                         MEF90TimeInterpolation_linear, & ! timeInterpolation
                                                         0.0_Kr,        & ! timeMin
                                                         1.0_Kr,        & ! timeMax
                                                         11,            & ! timeNumStep
                                                         MEF90FileFormat_EXOSingle) ! fileFormat
   Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions
   PetscBag,dimension(:),pointer                      :: MEF90MatPropBag

   Type(DM),target                                    :: Mesh
   Type(IS)                                           :: setIS,cellIS
   PetscInt,Dimension(:),Pointer                      :: setID
   PetscInt                                           :: numset,set
   Type(SectionReal)                                  :: defaultSection
   Type(Vec),target                                   :: temperature
   Type(Vec),target                                   :: fluxPrevious,fluxTarget,flux
   Type(Vec),target                                   :: boundaryTemperaturePrevious,boundaryTemperatureTarget,boundaryTemperature
   Type(Vec),target                                   :: externalTemperaturePrevious,externalTemperatureTarget,externalTemperature
   Type(Vec)                                          :: residual,RHS
   PetscReal,Dimension(:),Pointer                     :: time

   Type(SNES)                                         :: snesTemp
   Type(TS)                                           :: tsTemp
   Type(KSP)                                          :: kspTemp
   Type(PC)                                           :: pcTemp
   Type(Mat)                                          :: matTemp
   Type(MatNullSpace)                                 :: nspTemp
   PetscReal                                          :: rtol,atol,dtol
   PetscInt                                           :: maxits
          
   PetscBool                                          :: flg
   Character(len=MEF90_MXSTRLEN)                      :: IOBuffer
   Character(len=MEF90_MXSTRLEN)                      :: setName,setprefix
   Character(len=MXSTLN),Dimension(:),Pointer         :: nameG,nameC,nameV
   Integer                                            :: numfield
   
   Integer                                            :: step
   Type(Vec)                                          :: localVec
   
   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   Call MEF90_Initialize(ierr)

   !!! Initialize MEF90
   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   Call MEF90_Initialize(ierr)

   !!! Get all MEF90-wide options
   Call MEF90Ctx_Create(PETSC_COMM_WORLD,MEF90Ctx,MEF90DefaultGlobalOptions,ierr);CHKERRQ(ierr)
   Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)

   !!! Get DM from mesh
   Call MEF90Ctx_GetDMMeshEXO(MEF90Ctx,Mesh,ierr);CHKERRQ(ierr)
   Call DMMeshSetMaxDof(mesh,1,ierr);CHKERRQ(ierr) 
   
   !!! Open output file
   Call MEF90Ctx_OpenEXO(MEF90Ctx,Mesh,ierr)

   !!! Create HeatXfer context, get all HeatXfer options
   Call MEF90HeatXferCtx_Create(MEF90HeatXferCtx,Mesh,MEF90Ctx,ierr);CHKERRQ(ierr)
   Call MEF90HeatXferCtx_SetFromOptions(MEF90HeatXferCtx,PETSC_NULL_CHARACTER,MEF90HeatXferDefaultGlobalOptions, &
                                        MEF90HeatXferDefaultCellSetOptions,MEF90HeatXferDefaultVertexSetOptions,ierr)
   Call PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr);CHKERRQ(ierr)

   !!! Get material properties bags
   Call MEF90MatPropBag_SetFromOptions(MEF90MatPropBag,MEF90HeatXferCtx%DM,MEF90_Mathium2D,MEF90Ctx,ierr)
   MEF90HeatXferCtx%MaterialPropertiesBag => MEF90MatPropBag

   Call MEF90Ctx_GetTime(MEF90Ctx,time,ierr)

   !!! Create default section matching element type
   Call DMMeshGetVertexSectionReal(MEF90HeatXferCtx%DM,"default",1,defaultSection,ierr);CHKERRQ(ierr)
   Call DMMeshSetSectionReal(MEF90HeatXferCtx%DM,"default",defaultSection,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(defaultSection,ierr);CHKERRQ(ierr)

   Call DMMeshGetCellSectionReal(MEF90HeatXferCtx%cellDM,"default",1,defaultSection,ierr);CHKERRQ(ierr)
   Call DMMeshSetSectionReal(MEF90HeatXferCtx%cellDM,"default",defaultSection,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(defaultSection,ierr);CHKERRQ(ierr)
      
   !!! Create vectors
   Call DMCreateGlobalVector(MEF90HeatXferCtx%DM,temperature,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(temperature,"temperature",ierr);CHKERRQ(ierr)
   Call DMCreateGlobalVector(MEF90HeatXferCtx%DM,residual,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(residual,"residual",ierr);CHKERRQ(ierr)
   Call DMCreateGlobalVector(MEF90HeatXferCtx%DM,RHS,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(MEF90HeatXferCtx%DM,"RHS",ierr);CHKERRQ(ierr)

   Call DMCreateGlobalVector(MEF90HeatXferCtx%DM,boundaryTemperatureTarget,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(boundaryTemperatureTarget,"boundary Temperature",ierr);CHKERRQ(ierr)
   MEF90HeatXferCtx%boundaryTemperatureTarget => boundaryTemperatureTarget
   Call DMCreateGlobalVector(MEF90HeatXferCtx%DM,boundaryTemperaturePrevious,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(boundaryTemperaturePrevious,"boundary Temperature (previous time step)",ierr);CHKERRQ(ierr)
   MEF90HeatXferCtx%boundaryTemperaturePrevious => boundaryTemperaturePrevious
   
   Call DMCreateGlobalVector(MEF90HeatXferCtx%cellDM,externalTemperatureTarget,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(externalTemperatureTarget,"external Temperature",ierr);CHKERRQ(ierr)
   MEF90HeatXferCtx%externalTemperatureTarget => externalTemperatureTarget
   Call DMCreateGlobalVector(MEF90HeatXferCtx%cellDM,externalTemperaturePrevious,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(externalTemperaturePrevious,"external Temperature (previous time step)",ierr);CHKERRQ(ierr)
   MEF90HeatXferCtx%externalTemperaturePrevious => externalTemperaturePrevious

   Call DMCreateGlobalVector(MEF90HeatXferCtx%cellDM,fluxTarget,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(fluxTarget,"Flux",ierr);CHKERRQ(ierr)
   MEF90HeatXferCtx%fluxTarget => fluxTarget
   Call DMCreateGlobalVector(MEF90HeatXferCtx%cellDM,fluxPrevious,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(fluxPrevious,"Flux (previous time step)",ierr);CHKERRQ(ierr)
   MEF90HeatXferCtx%fluxPrevious => fluxPrevious

   !!! 
   !!! Create SNES or TS, Mat and set KSP default options
   !!!
   Call DMCreateMatrix(MEF90HeatXferCtx%DM,MATAIJ,matTemp,iErr);CHKERRQ(iErr)
   Call MatSetOption(matTemp,MAT_SPD,PETSC_TRUE,ierr);CHKERRQ(ierr)
   Call MatSetOption(matTemp,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr);CHKERRQ(ierr)
   Call MatSetOption(matTemp,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr);CHKERRQ(ierr)
   Call MatSetFromOptions(matTemp,ierr);CHKERRQ(ierr)

   If (MEF90HeatXferGlobalOptions%mode == MEF90HeatXFer_ModeSteadyState) Then
      Call SNESCreate(PETSC_COMM_WORLD,snesTemp,ierr);CHKERRQ(ierr)
      Call SNESSetDM(snesTemp,MEF90HeatXferCtx%DM,ierr);CHKERRQ(ierr)
      Call SNESSetOptionsPrefix(snesTemp,'temp_',ierr);CHKERRQ(ierr)

      Call SNESSetFunction(snesTemp,residual,MEF90HeatXferOperator,MEF90HeatXferCtx,ierr);CHKERRQ(ierr)
      Call SNESSetJacobian(snesTemp,matTemp,matTemp,MEF90HeatXferBilinearForm,MEF90HeatXferCtx,ierr);CHKERRQ(ierr)
      Call SNESSetFromOptions(snesTemp,ierr);CHKERRQ(ierr)
      If (MEF90GlobalOptions%verbose > 0) Then
         Call SNESView(snesTemp,PETSC_VIEWER_STDOUT_WORLD,ierr)
      End If
   Else
      Call TSCreate(PETSC_COMM_WORLD,tsTemp,ierr);CHKERRQ(ierr)
      Call TSSetDM(tsTemp,MEF90HeatXferCtx%DM,ierr);CHKERRQ(ierr)
      Call TSSetOptionsPrefix(tsTemp,'temp_',ierr);CHKERRQ(ierr)
      Call TSGetSNES(tsTemp,snesTemp,ierr);CHKERRQ(ierr)

      !Call TSSetIFunction(tsTemp,resTemp,SimplePoissonTSIFunction,MEF90Ctx,ierr);CHKERRQ(ierr)
      !Call TSSetIJacobian(tsTemp,matTemp,matTemp,SimplePoissonTSIJacobian,MEF90Ctx,ierr);CHKERRQ(ierr)
      !If ((GlobalProperties%LoadingType == Poisson_MIL) .OR. (GlobalProperties%LoadingType == Poisson_CST)) Then
      !   Call TSSetRHSFunction(tsTemp,PETSC_NULL_OBJECT,SimplePoissonTSRHS_Cst,MEF90Ctx,ierr);CHKERRQ(ierr)
      !Else
      !   !Call TSSetRHSFunction(tsTemp,PETSC_NULL_OBJECT,SimplePoissonTSRHS,MEF90Ctx,ierr);CHKERRQ(ierr)
      !End If
      Call TSSetType(tsTemp,'rosw',ierr);CHKERRQ(ierr)
      Call TSRosWSetType(tsTemp,'ra3pw',ierr);CHKERRQ(ierr)
      !TSinitialStep = (time(size(time))-time(0)) / (size(time) + 0.0_Kr) / 10.0_Kr
      !TSinitialTime = time(1)
      !Call TSSetInitialTimeStep(tsTemp,TSinitialTime,TSinitialStep,ierr);CHKERRQ(ierr)
      !Call TSSetProblemType(tsTemp,TS_LINEAR,ierr);CHKERRQ(ierr)
      !Call VecSet(solTemp,GlobalProperties%initialTemp,ierr);CHKERRQ(ierr)
      !Call TSSetSolution(tsTemp,solTemp,ierr);CHKERRQ(ierr)
      !Call TSSetFromOptions(tsTemp,ierr);CHKERRQ(ierr)
      If (MEF90GlobalOptions%verbose > 0) Then
         Call TSView(tsTemp,PETSC_VIEWER_STDOUT_WORLD,ierr)
      End If
   End If
   If (MEF90HeatXferGlobalOptions%addNullSpace) Then
      Call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,PETSC_NULL_OBJECT,nspTemp,ierr);CHKERRQ(ierr)
      Call MatSetNullSpace(matTemp,nspTemp,ierr);CHKERRQ(ierr)
   End If
   !!! 
   !!! Set some KSP options
   !!!
   Call SNESGetKSP(snesTemp,kspTemp,ierr);CHKERRQ(ierr)
   Call KSPSetType(kspTemp,KSPCG,ierr);CHKERRQ(ierr)
   Call KSPSetInitialGuessNonzero(kspTemp,PETSC_TRUE,ierr);CHKERRQ(ierr)
   If (MEF90HeatXferGlobalOptions%addNullSpace) Then
      Call KSPSetNullSpace(kspTemp,nspTemp,ierr);CHKERRQ(ierr)
   End If
   Call KSPGetTolerances(kspTemp,rtol,atol,dtol,maxits,ierr);CHKERRQ(ierr)
   rtol = 1.0D-8
   Call KSPSetTolerances(kspTemp,rtol,atol,dtol,maxits,ierr);CHKERRQ(ierr)
   Call KSPSetFromOptions(kspTemp,ierr);CHKERRQ(ierr)
   


   !!!
   !!! Try to figure out if the file was formatted
   !!!
   If (MEF90Ctx%rank == 0) Then
      Call EXGVP(MEF90Ctx%fileExoUnit,"N",numfield,ierr)
   End If
   Call MPI_Bcast(numfield,1,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)

   If (numfield == 0) Then
      Allocate(nameG(2))
      nameG(1) = "Energy"
      nameG(2) = "work"
   
      numfield = max(MEF90HeatXferGlobalOptions%tempOffset, &
                     MEF90HeatXferGlobalOptions%boundaryTempOffset)
      Allocate(nameV(numfield))
      nameV = "empty"
      nameV(MEF90HeatXferGlobalOptions%tempOffset) = "temperature"
      nameV(MEF90HeatXferGlobalOptions%boundaryTempOffset) = "boundary temperature"
                     
      numfield = max(MEF90HeatXferGlobalOptions%externalTempOffset, &
                     MEF90HeatXferGlobalOptions%fluxOffset)
      Allocate(nameC(numfield))
      nameC = "empty"
      nameC(MEF90HeatXferGlobalOptions%externalTempOffset) = "external temperature"
      nameC(MEF90HeatXferGlobalOptions%fluxOffset) = "heat flux"

      Call MEF90EXOFormat(MEF90Ctx%fileEXOUNIT,nameG,nameC,nameV,ierr)
   End If
   
   !!!
   !!! Actual computations / time stepping
   !!!
   If (MEF90HeatXferGlobalOptions%mode == MEF90HeatXFer_ModeSteadyState) Then
      Call MEF90HeatXferGetTransients(MEF90HeatXferCtx,1,time(1),ierr)
      Do step = 1,MEF90GlobalOptions%timeNumStep
         Write(IOBuffer,100) step,time(step)
         Call PetscPrintf(MEF90Ctx%comm,IOBuffer,ierr);CHKERRQ(ierr)
         !!! Update fields
         Call MEF90HeatXferGetTransients(MEF90HeatXferCtx,step,time(step),ierr)

         !!! Solve SNES
         Call SNESSolve(snesTemp,PETSC_NULL_OBJECT,temperature,ierr);CHKERRQ(ierr)
         
         !!! Compute energies
         
         
         !!! Save results
         Call DMGetLocalVector(MEF90HeatXferCtx%cellDM,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalBegin(MEF90HeatXferCtx%cellDM,MEF90HeatXferCtx%fluxTarget,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalEnd(MEF90HeatXferCtx%cellDM,MEF90HeatXferCtx%fluxTarget,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call VecViewExodusCell(MEF90HeatXferCtx%cellDM,localVec,MEF90HeatXferCtx%MEF90Ctx%IOcomm, &
                                MEF90HeatXferCtx%MEF90Ctx%fileExoUnit,step,MEF90HeatXferGlobalOptions%fluxOffset,ierr);CHKERRQ(ierr)

         Call DMGlobalToLocalBegin(MEF90HeatXferCtx%cellDM,MEF90HeatXferCtx%externalTemperatureTarget,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalEnd(MEF90HeatXferCtx%cellDM,MEF90HeatXferCtx%externalTemperatureTarget,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call VecViewExodusCell(MEF90HeatXferCtx%cellDM,localVec,MEF90HeatXferCtx%MEF90Ctx%IOcomm, &
                                MEF90HeatXferCtx%MEF90Ctx%fileExoUnit,step,MEF90HeatXferGlobalOptions%externalTempOffset,ierr);CHKERRQ(ierr)
         Call DMGetLocalVector(MEF90HeatXferCtx%DM,localVec,ierr);CHKERRQ(ierr)

         Call DMGlobalToLocalBegin(MEF90HeatXferCtx%DM,temperature,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalEnd(MEF90HeatXferCtx%DM,temperature,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call VecViewExodusVertex(MEF90HeatXferCtx%DM,localVec,MEF90HeatXferCtx%MEF90Ctx%IOcomm, &
                                  MEF90HeatXferCtx%MEF90Ctx%fileExoUnit,step,MEF90HeatXferGlobalOptions%tempOffset,ierr);CHKERRQ(ierr)

         Call DMGlobalToLocalBegin(MEF90HeatXferCtx%DM,MEF90HeatXferCtx%boundaryTemperatureTarget,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalEnd(MEF90HeatXferCtx%DM,MEF90HeatXferCtx%boundaryTemperatureTarget,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call VecViewExodusVertex(MEF90HeatXferCtx%DM,localVec,MEF90HeatXferCtx%MEF90Ctx%IOcomm, &
                                  MEF90HeatXferCtx%MEF90Ctx%fileExoUnit,step,MEF90HeatXferGlobalOptions%boundaryTempOffset,ierr);CHKERRQ(ierr)
         Call DMRestoreLocalVector(MEF90HeatXferCtx%DM,localVec,ierr);CHKERRQ(ierr)
      End Do
   !Else
   End If
100 Format("Solving steady state step ",I4,", t=",ES12.5,"\n")
   !!! Clean up and exit nicely
   If (MEF90HeatXferGlobalOptions%mode == MEF90HeatXFer_ModeSteadyState) Then
      Call SNESDestroy(snesTemp,ierr);CHKERRQ(ierr)
   Else
      Call TSDestroy(tsTemp,ierr);CHKERRQ(ierr)
   End If

   Call VecDestroy(Temperature,ierr);CHKERRQ(ierr)
   Call VecDestroy(residual,ierr);CHKERRQ(ierr)
   Call VecDestroy(RHS,ierr);CHKERRQ(ierr)
   
   If (Associated(MEF90HeatXferCtx%boundaryTemperaturePrevious)) Then 
      Call VecDestroy(MEF90HeatXferCtx%boundaryTemperaturePrevious,ierr);CHKERRQ(ierr)
      Nullify(MEF90HeatXferCtx%boundaryTemperaturePrevious)
   End If
   If (Associated(MEF90HeatXferCtx%boundaryTemperatureTarget)) Then 
      Call VecDestroy(MEF90HeatXferCtx%boundaryTemperatureTarget,ierr);CHKERRQ(ierr)
      Nullify(MEF90HeatXferCtx%boundaryTemperatureTarget)
   End If
   
   If (Associated(MEF90HeatXferCtx%externalTemperaturePrevious)) Then 
      Call VecDestroy(MEF90HeatXferCtx%externalTemperaturePrevious,ierr);CHKERRQ(ierr)
      Nullify(MEF90HeatXferCtx%externalTemperaturePrevious)
   End If
   If (Associated(MEF90HeatXferCtx%externalTemperatureTarget)) Then 
      Call VecDestroy(MEF90HeatXferCtx%externalTemperatureTarget,ierr);CHKERRQ(ierr)
      Nullify(MEF90HeatXferCtx%externalTemperatureTarget)
   End If

   If (Associated(MEF90HeatXferCtx%fluxPrevious)) Then 
      Call VecDestroy(MEF90HeatXferCtx%fluxPrevious,ierr);CHKERRQ(ierr)
      Nullify(MEF90HeatXferCtx%fluxPrevious)
   End If
   If (Associated(MEF90HeatXferCtx%fluxTarget)) Then 
      Call VecDestroy(MEF90HeatXferCtx%fluxTarget,ierr);CHKERRQ(ierr)
      Nullify(MEF90HeatXferCtx%fluxTarget)
   End If

   DeAllocate(time)
   Call MEF90HeatXferCtx_Destroy(MEF90HeatXferCtx,ierr);CHKERRQ(ierr)
   Call MEF90Ctx_CloseEXO(MEF90Ctx,ierr)
   Call MEF90_Finalize(ierr)
   Call PetscFinalize()
End Program TestHeatXfer