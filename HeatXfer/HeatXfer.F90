#include "../MEF90/mef90.inc"
Program TestHeatXfer
#include <finclude/petscdef.h>
   Use m_MEF90
   Use m_MEF90_HeatXfer
   Use m_MEF90_HeatXferCtx
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
                                                         MEF90Scaling_Linear, & ! fluxScaling
                                                         2)                     ! fluxOffset
   Type(MEF90HeatXferCellSetOptions_Type),Parameter   :: MEF90HeatXferDefaultCellSetOptions = MEF90HeatXferCellSetOptions_Type( &
                                                         -1,            & ! elemTypeShortID will be overriden
                                                         0.0_Kr,        & ! flux
                                                         0.0_Kr,        & ! surfaceThermalConductivity
                                                         0.0_Kr,        & ! externalTemp
                                                         PETSC_FALSE,   & ! Has BC
                                                         0.0_Kr)          ! boundaryTemp
                                                         
   Type(MEF90HeatXferVertexSetOptions_Type),Parameter :: MEF90HeatXferDefaultVertexSetOptions = MEF90HeatXferVertexSetOptions_Type( &
                                                         PETSC_FALSE,    & ! Has BC
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
   Type(IS)                                           :: setIS,cellIS,CellSetGlobalIS
   PetscInt,Dimension(:),Pointer                      :: setID
   PetscInt                                           :: numset,set
   Type(SectionReal)                                  :: defaultSection
   Type(Vec),target                                   :: temperature
   Type(Vec),target                                   :: flux
   Type(Vec),target                                   :: boundaryTemperature
   Type(Vec),target                                   :: externalTemperature
   Type(Vec)                                          :: residual,RHS
   PetscReal,Dimension(:),Pointer                     :: time,energy,work

   Type(SNES)                                         :: snesTemp
   Type(TS)                                           :: tsTemp
   Type(KSP)                                          :: kspTemp
   Type(PC)                                           :: pcTemp
   Type(Mat)                                          :: matTemp
   Type(MatNullSpace)                                 :: nspTemp
   PetscReal                                          :: rtol,dtol
          
   PetscBool                                          :: flg
   Character(len=MEF90_MXSTRLEN)                      :: IOBuffer
   Character(len=MEF90_MXSTRLEN)                      :: setName,setprefix
   Character(len=MXSTLN),Dimension(:),Pointer         :: nameG,nameC,nameV
   Integer                                            :: numfield
   
   Integer                                            :: step
   Type(Vec)                                          :: localVec
   PetscInt                                           :: dim
      
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
   Call DMMeshGetDimension(Mesh,dim,ierr);CHKERRQ(ierr)
   If (dim == 2) Then
      Call MEF90MatPropBag_SetFromOptions(MEF90MatPropBag,MEF90HeatXferCtx%DM,MEF90_Mathium2D,MEF90Ctx,ierr)
   Else
      Call MEF90MatPropBag_SetFromOptions(MEF90MatPropBag,MEF90HeatXferCtx%DM,MEF90_Mathium3D,MEF90Ctx,ierr)
   End If   
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
   Call PetscObjectSetName(RHS,"RHS",ierr);CHKERRQ(ierr)

   Call DMCreateGlobalVector(MEF90HeatXferCtx%DM,boundaryTemperature,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(boundaryTemperature,"boundary Temperature",ierr);CHKERRQ(ierr)
   MEF90HeatXferCtx%boundaryTemperature => boundaryTemperature
   
   Call DMCreateGlobalVector(MEF90HeatXferCtx%cellDM,externalTemperature,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(externalTemperature,"external Temperature",ierr);CHKERRQ(ierr)
   MEF90HeatXferCtx%externalTemperature => externalTemperature

   Call DMCreateGlobalVector(MEF90HeatXferCtx%cellDM,flux,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(flux,"Flux",ierr);CHKERRQ(ierr)
   MEF90HeatXferCtx%flux => flux

   !!! 
   !!! Create SNES or TS, Mat and set KSP default options
   !!!
   Call DMCreateMatrix(MEF90HeatXferCtx%DM,MATAIJ,matTemp,iErr);CHKERRQ(iErr)
   Call MatSetOptionsPrefix(matTemp,"temp_",ierr);CHKERRQ(ierr)
   Call MatSetOption(matTemp,MAT_SPD,PETSC_TRUE,ierr);CHKERRQ(ierr)
   Call MatSetOption(matTemp,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr);CHKERRQ(ierr)
   Call MatSetOption(matTemp,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr);CHKERRQ(ierr)
   If (MEF90HeatXferGlobalOptions%addNullSpace) Then
      Call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,PETSC_NULL_OBJECT,nspTemp,ierr);CHKERRQ(ierr)
      Call MatSetNullSpace(matTemp,nspTemp,ierr);CHKERRQ(ierr)
   End If
   Call MatSetFromOptions(matTemp,ierr);CHKERRQ(ierr)

   If (MEF90HeatXferGlobalOptions%mode == MEF90HeatXFer_ModeSteadyState) Then
      Call SNESCreate(PETSC_COMM_WORLD,snesTemp,ierr);CHKERRQ(ierr)
      Call SNESSetApplicationContext(snesTemp,MEF90HeatXferCtx,ierr);CHKERRQ(ierr)
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
   !!! 
   !!! Set some KSP options
   !!!
   Call SNESGetKSP(snesTemp,kspTemp,ierr);CHKERRQ(ierr)
   Call KSPSetType(kspTemp,KSPCG,ierr);CHKERRQ(ierr)
   Call KSPSetInitialGuessNonzero(kspTemp,PETSC_TRUE,ierr);CHKERRQ(ierr)
   If (MEF90HeatXferGlobalOptions%addNullSpace) Then
      Call KSPSetNullSpace(kspTemp,nspTemp,ierr);CHKERRQ(ierr)
   End If
   rtol = 1.0D-8
   dtol = 1.0D+10
   Call KSPSetTolerances(kspTemp,rtol,PETSC_DEFAULT_DOUBLE_PRECISION,dtol,PETSC_DEFAULT_INTEGER,ierr);CHKERRQ(ierr)
   Call KSPSetFromOptions(kspTemp,ierr);CHKERRQ(ierr)
   
   !!! 
   !!! Allocate array of works and energies
   !!!
   Allocate(energy(size(MEF90HeatXferCtx%CellSetOptionsBag)))
   Allocate(work(size(MEF90HeatXferCtx%CellSetOptionsBag)))

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
      Do step = 1,MEF90GlobalOptions%timeNumStep
         Write(IOBuffer,100) step,time(step)
         Call PetscPrintf(MEF90Ctx%comm,IOBuffer,ierr);CHKERRQ(ierr)

         !!! Update fields
         Call MEF90HeatXferSetTransients(MEF90HeatXferCtx,step,time(step),ierr)

         !!! Solve SNES
         Call MEF90HeatXferUpdateboundaryTemperature(temperature,MEF90HeatXferCtx,ierr);
         Call SNESSolve(snesTemp,PETSC_NULL_OBJECT,temperature,ierr);CHKERRQ(ierr)
         
         !!! Compute energies
         Call MEF90HeatXFerEnergy(temperature,time(step),MEF90HeatXferCtx,energy,work,ierr);CHKERRQ(ierr)
         Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
         Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
         Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
         Do set = 1, size(setID)
            Write(IOBuffer,101) setID(set),energy(set),work(set)
            Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
         End Do
         Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
         Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
         Write(IOBuffer,102) sum(energy),sum(work)
         Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
     
         
         !!! Save results
         Call DMGetLocalVector(MEF90HeatXferCtx%cellDM,localVec,ierr);CHKERRQ(ierr)
         If (MEF90HeatXferGlobalOptions%fluxOffset > 0) Then
            Call DMGlobalToLocalBegin(MEF90HeatXferCtx%cellDM,MEF90HeatXferCtx%flux,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
            Call DMGlobalToLocalEnd(MEF90HeatXferCtx%cellDM,MEF90HeatXferCtx%flux,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
            Call VecViewExodusCell(MEF90HeatXferCtx%cellDM,localVec,MEF90HeatXferCtx%MEF90Ctx%IOcomm, &
                                   MEF90HeatXferCtx%MEF90Ctx%fileExoUnit,step,MEF90HeatXferGlobalOptions%fluxOffset,ierr);CHKERRQ(ierr)
         End If
         
         If (MEF90HeatXferGlobalOptions%externalTempOffset > 0) Then
            Call DMGlobalToLocalBegin(MEF90HeatXferCtx%cellDM,MEF90HeatXferCtx%externalTemperature,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
            Call DMGlobalToLocalEnd(MEF90HeatXferCtx%cellDM,MEF90HeatXferCtx%externalTemperature,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
            Call VecViewExodusCell(MEF90HeatXferCtx%cellDM,localVec,MEF90HeatXferCtx%MEF90Ctx%IOcomm, &
                                   MEF90HeatXferCtx%MEF90Ctx%fileExoUnit,step,MEF90HeatXferGlobalOptions%externalTempOffset,ierr);CHKERRQ(ierr)
         End If
         Call DMRestoreLocalVector(MEF90HeatXferCtx%cellDM,localVec,ierr);CHKERRQ(ierr)
         
         Call DMGetLocalVector(MEF90HeatXferCtx%DM,localVec,ierr);CHKERRQ(ierr)
         If (MEF90HeatXferGlobalOptions%tempOffset > 0) Then
            Call DMGlobalToLocalBegin(MEF90HeatXferCtx%DM,temperature,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
            Call DMGlobalToLocalEnd(MEF90HeatXferCtx%DM,temperature,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
            Call VecViewExodusVertex(MEF90HeatXferCtx%DM,localVec,MEF90HeatXferCtx%MEF90Ctx%IOcomm, &
                                     MEF90HeatXferCtx%MEF90Ctx%fileExoUnit,step,MEF90HeatXferGlobalOptions%tempOffset,ierr);CHKERRQ(ierr)
         End If

         If (MEF90HeatXferGlobalOptions%boundaryTempOffset > 0) Then
            Call DMGlobalToLocalBegin(MEF90HeatXferCtx%DM,MEF90HeatXferCtx%boundaryTemperature,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
            Call DMGlobalToLocalEnd(MEF90HeatXferCtx%DM,MEF90HeatXferCtx%boundaryTemperature,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
            Call VecViewExodusVertex(MEF90HeatXferCtx%DM,localVec,MEF90HeatXferCtx%MEF90Ctx%IOcomm, &
                                     MEF90HeatXferCtx%MEF90Ctx%fileExoUnit,step,MEF90HeatXferGlobalOptions%boundaryTempOffset,ierr);CHKERRQ(ierr)
            End If
         Call DMRestoreLocalVector(MEF90HeatXferCtx%DM,localVec,ierr);CHKERRQ(ierr)
      End Do
   !Else
   End If
100 Format("Solving steady state step ",I4,", t=",ES12.5,"\n")
101 Format("cell set ",I4," thermal energy: ",ES12.5," fluxes work: ",ES12.5,"\n")
102 Format("======= Total thermal energy: ",ES12.5," fluxes work: ",ES12.5,"\n")
   !!! Clean up and exit nicely
   If (MEF90HeatXferGlobalOptions%mode == MEF90HeatXFer_ModeSteadyState) Then
      Call SNESDestroy(snesTemp,ierr);CHKERRQ(ierr)
   Else
      Call TSDestroy(tsTemp,ierr);CHKERRQ(ierr)
   End If

   Call VecDestroy(Temperature,ierr);CHKERRQ(ierr)
   Call VecDestroy(residual,ierr);CHKERRQ(ierr)
   Call VecDestroy(RHS,ierr);CHKERRQ(ierr)
   
   If (Associated(MEF90HeatXferCtx%boundaryTemperature)) Then 
      Call VecDestroy(MEF90HeatXferCtx%boundaryTemperature,ierr);CHKERRQ(ierr)
      Nullify(MEF90HeatXferCtx%boundaryTemperature)
   End If
   
   If (Associated(MEF90HeatXferCtx%externalTemperature)) Then 
      Call VecDestroy(MEF90HeatXferCtx%externalTemperature,ierr);CHKERRQ(ierr)
      Nullify(MEF90HeatXferCtx%externalTemperature)
   End If

   If (Associated(MEF90HeatXferCtx%flux)) Then 
      Call VecDestroy(MEF90HeatXferCtx%flux,ierr);CHKERRQ(ierr)
      Nullify(MEF90HeatXferCtx%flux)
   End If

   DeAllocate(time)
   DeAllocate(energy)
   DeAllocate(work)
   Call PetscLogView(MEF90Ctx%logViewer,ierr);CHKERRQ(ierr)
   Call MEF90HeatXferCtx_Destroy(MEF90HeatXferCtx,ierr);CHKERRQ(ierr)
   Call MEF90Ctx_CloseEXO(MEF90Ctx,ierr)
   Call MEF90_Finalize(ierr)
   Call PetscFinalize()
End Program TestHeatXfer
