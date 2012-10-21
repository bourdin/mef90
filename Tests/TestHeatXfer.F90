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
                                                         PETSC_FALSE,   & ! addNullSpace
                                                         1,             & ! tempOffset
                                                         PETSC_FALSE,   & ! boundaryTempCst
                                                         2,             & ! boundaryTempOffset
                                                         PETSC_TRUE,    & ! externalTempCst
                                                         3,             & ! externalTempOffset
                                                         PETSC_TRUE,    & ! fluxCst
                                                         4)               ! fluxOffset
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
   Type(PetscBag),dimension(:),pointer                :: MEF90MatPropBag

   Type(DM),target                                    :: Mesh
   Type(IS)                                           :: setIS,cellIS
   PetscInt,Dimension(:),Pointer                      :: setID
   PetscInt                                           :: numset,set
   Type(SectionReal)                                  :: defaultSection
   Type(Vec),target                                   :: temperature,flux,boundaryTemperature,externalTemperature
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
   
   PetscInt :: point,numCell,numVertex

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
   Call DMView(Mesh,PETSC_VIEWER_STDOUT_WORLD,ierr)
   
   !!! Open output file
   Call MEF90Ctx_OpenEXO(MEF90Ctx,Mesh,ierr)

   !!! Create HeatXfer context, get all HeatXfer options
   Call MEF90HeatXferCtx_Create(MEF90HeatXferCtx,Mesh,MEF90Ctx,ierr);CHKERRQ(ierr)
   Call MEF90HeatXferCtx_SetFromOptions(MEF90HeatXferCtx,PETSC_NULL_CHARACTER,Mesh,MEF90Ctx,MEF90HeatXferDefaultGlobalOptions, &
                                        MEF90HeatXferDefaultCellSetOptions,MEF90HeatXferDefaultVertexSetOptions,ierr)
   Call PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr);CHKERRQ(ierr)

   !!! Get material properties bags
   Call MEF90MatPropBag_SetFromOptions(MEF90MatPropBag,Mesh,MEF90_Mathium2D,MEF90Ctx,ierr)
   MEF90HeatXferCtx%MaterialPropertiesBag => MEF90MatPropBag

   Call MEF90Ctx_GetTime(MEF90Ctx,time,ierr)

   !!! Create default section matching element type
   Call DMMeshGetSectionReal(mesh,'default',defaultSection,ierr);CHKERRQ(ierr)
   Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
   Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1, size(setID)
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),cellIS,ierr);CHKERRQ(iErr)
      Call SectionRealSetFiberDimensionSet(Mesh,defaultSection,cellIS,1,ierr)
      Call ISDestroy(cellIS,ierr);CHKERRQ(ierr)
   End Do
   Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
   Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   Call SectionRealAllocate(defaultSection,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(defaultSection,ierr);CHKERRQ(ierr)
      
   Call DMCreateGlobalVector(Mesh,temperature,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(temperature,"temperature",ierr);CHKERRQ(ierr)
   Call DMCreateGlobalVector(Mesh,residual,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(residual,"residual",ierr);CHKERRQ(ierr)
   Call DMCreateGlobalVector(Mesh,RHS,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(RHS,"RHS",ierr);CHKERRQ(ierr)

   If (.NOT. MEF90HeatXferGlobalOptions%boundaryTempCst) Then 
      Call DMCreateGlobalVector(Mesh,boundaryTemperature,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(boundaryTemperature,"boundary Temperature",ierr);CHKERRQ(ierr)
      MEF90HeatXferCtx%boundaryValue => boundaryTemperature
   End If
   If (.NOT. MEF90HeatXferGlobalOptions%fluxCst) Then 
      Call DMCreateGlobalVector(Mesh,flux,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(flux,"Heat flux",ierr);CHKERRQ(ierr)
      MEF90HeatXferCtx%flux => flux
   End If
   If (.NOT. MEF90HeatXferGlobalOptions%externalTempCst) Then 
      Call DMCreateGlobalVector(Mesh,externalTemperature,ierr);CHKERRQ(ierr)
      Call PetscObjectSetName(externalTemperature,"external Temperature",ierr);CHKERRQ(ierr)
      MEF90HeatXferCtx%externalValue => externalTemperature
   End If
   
   !!! 
   !!! Create SNES or TS, Mat and set KSP default options
   !!!
   Call DMCreateMatrix(mesh,MATAIJ,matTemp,iErr);CHKERRQ(iErr)
   Call MatSetOption(matTemp,MAT_SPD,PETSC_TRUE,ierr);CHKERRQ(ierr)
   Call MatSetOption(matTemp,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr);CHKERRQ(ierr)
   Call MatSetOption(matTemp,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr);CHKERRQ(ierr)
   Call MatSetFromOptions(matTemp,ierr);CHKERRQ(ierr)

   If (MEF90HeatXferGlobalOptions%mode == MEF90HeatXFer_ModeSteadyState) Then
      Call SNESCreate(PETSC_COMM_WORLD,snesTemp,ierr);CHKERRQ(ierr)
      Call SNESSetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)
      Call SNESSetOptionsPrefix(snesTemp,'temp_',ierr);CHKERRQ(ierr)

      !Call SNESSetFunction(snesTemp,resTemp,MEF90HeatXferOperator,MEF90HeatXferCtx,ierr);CHKERRQ(ierr)
      Call SNESSetJacobian(snesTemp,matTemp,matTemp,MEF90HeatXferBilinearForm,MEF90HeatXferCtx,ierr);CHKERRQ(ierr)
      Call SNESSetFromOptions(snesTemp,ierr);CHKERRQ(ierr)
      If (MEF90GlobalOptions%verbose > 0) Then
         Call SNESView(snesTemp,PETSC_VIEWER_STDOUT_WORLD,ierr)
      End If
   Else
      Call TSCreate(PETSC_COMM_WORLD,tsTemp,ierr);CHKERRQ(ierr)
      Call TSSetDM(tsTemp,mesh,ierr);CHKERRQ(ierr)
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

   
   

   !!! Clean up and exit nicely
   If (MEF90HeatXferGlobalOptions%mode == MEF90HeatXFer_ModeSteadyState) Then
      Call SNESDestroy(snesTemp,ierr);CHKERRQ(ierr)
   Else
      Call TSDestroy(tsTemp,ierr);CHKERRQ(ierr)
   End If

   Call VecDestroy(Temperature,ierr);CHKERRQ(ierr)
   Call VecDestroy(residual,ierr);CHKERRQ(ierr)
   Call VecDestroy(RHS,ierr);CHKERRQ(ierr)
   If (.NOT. MEF90HeatXferGlobalOptions%boundaryTempCst) Then 
      Call VecDestroy(boundaryTemperature,ierr);CHKERRQ(ierr)
   End If
   If (.NOT. MEF90HeatXferGlobalOptions%fluxCst) Then 
      Call VecDestroy(flux,ierr);CHKERRQ(ierr)
   End If
   If (.NOT. MEF90HeatXferGlobalOptions%externalTempCst) Then 
      Call VecDestroy(externalTemperature,ierr);CHKERRQ(ierr)
   End If

   DeAllocate(time)
   Call MEF90HeatXferCtx_Destroy(MEF90HeatXferCtx,ierr);CHKERRQ(ierr)
   Call MEF90Ctx_CloseEXO(MEF90Ctx,ierr)
   Call MEF90_Finalize(ierr)
   Call PetscFinalize()
End Program TestHeatXfer