Program  TestHeatXferCtx
#include <petsc/finclude/petsc.h>
    Use m_MEF90
    Use m_MEF90_HeatXfer
    Use petsc
    Implicit NONE   
    
    PetscErrorCode                                     :: ierr
    Type(MEF90Ctx_Type),target                         :: MEF90Ctx
    Type(MEF90CtxGlobalOptions_Type)                   :: MEF90GlobalOptions_default
    Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions
    Type(MEF90HeatXferCtx_Type)                        :: MEF90HeatXferCtx
    Type(tDM)                                          :: dm
    Type(tPetscSF)                                     :: naturalPointSF

    PetscInt                                           :: step
    PetscReal,Dimension(:),Pointer                     :: time
    Character(len=MEF90MXSTRLEN)                       :: IOBuffer
    PetscInt                                           :: numNodalVar = 1, numCellVar = 3, numGVar = 0
    Character(len=MEF90MXSTRLEN),Dimension(:),Pointer  :: nodalVarName, cellVarName, gVarName

    Type(MEF90HeatXferGlobalOptions_Type),Parameter    :: MEF90HeatXferDefaultGlobalOptions = MEF90HeatXferGlobalOptions_Type( &
        MEF90HeatXFer_timeSteppingTypeSteadyState, & ! timeSteppingType
        PETSC_FALSE,         & ! addNullSpace
        0.,                  & ! initialTemperature
        MEF90Scaling_Linear, & ! boundaryTemperatureScaling
        MEF90Scaling_Linear, & ! externalTemperatureScaling
        MEF90Scaling_Linear, & ! fluxScaling
        MEF90Scaling_Linear)   ! boundaryFluxScaling

    Type(MEF90HeatXferCellSetOptions_Type),Parameter   :: MEF90HeatXferDefaultCellSetOptions = MEF90HeatXferCellSetOptions_Type( &
        0.0_Kr,        & ! flux
        0.0_Kr,        & ! surfaceThermalConductivity
        0.0_Kr,        & ! externalTemperature
        PETSC_FALSE,   & ! Has BC
        0.0_Kr,        & ! boundaryTemperature
        [0.0_Kr,0.0_Kr,0.0_Kr]) ! AdvectionVector
        
    Type(MEF90HeatXferVertexSetOptions_Type),Parameter :: MEF90HeatXferDefaultVertexSetOptions = MEF90HeatXferVertexSetOptions_Type( &
        PETSC_FALSE,   & ! Has BC
        0.0_Kr)          ! boundaryTemperature


    PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
    PetscCallA(MEF90Initialize(ierr))

    MEF90GlobalOptions_default%verbose           = 1
    MEF90GlobalOptions_default%dryrun            = PETSC_FALSE
    MEF90GlobalOptions_default%timeInterpolation = MEF90TimeInterpolation_linear
    MEF90GlobalOptions_default%timeMin           = 0.0_Kr
    MEF90GlobalOptions_default%timeMax           = 1.0_Kr
    MEF90GlobalOptions_default%timeNumStep       = 11
    MEF90GlobalOptions_default%timeSkip          = 0
    MEF90GlobalOptions_default%timeNumCycle      = 1
    MEF90GlobalOptions_default%elementFamily     = MEF90ElementFamilyLagrange
    MEF90GlobalOptions_default%elementOrder      = 1
    
    PetscCallA(MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90GlobalOptions_default,ierr))
    PetscCallA(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr))

    PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm,MEF90Ctx%geometryfile,PETSC_NULL_CHARACTER,PETSC_TRUE,dm,ierr))
    PetscCallA(DMPlexDistributeSetDefault(dm,PETSC_FALSE,ierr))
    PetscCallA(DMSetFromOptions(dm,ierr))
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-mef90dm_view",ierr))
    PetscCallA(MEF90CtxOpenEXO(MEF90Ctx,MEF90Ctx%resultViewer,ierr))
    PetscCallA(MEF90EXODMView(dm,MEF90Ctx%resultViewer,MEF90GlobalOptions%elementOrder,ierr))

    distribute: Block 
        Type(tDM),target                    :: dmDist
        PetscInt                            :: ovlp = 0
        If (MEF90Ctx%NumProcs > 1) Then
            PetscCallA(DMPlexDistribute(dm,ovlp,naturalPointSF,dmDist,ierr))
            PetscCallA(DMPlexSetMigrationSF(dmDist,naturalPointSF, ierr))
            PetscCallA(PetscSFDestroy(naturalPointSF,ierr))
            PetscCallA(DMDestroy(dm,ierr))
            dm = dmDist
        End If
    End Block distribute
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-mef90dm_view",ierr))
    PetscCallA(MEF90CtxGetTime(MEF90Ctx,time,ierr))

    PetscCallA(MEF90HeatXferCtxCreate(MEF90HeatXferCtx,dm,MEF90Ctx,ierr))
    PetscCallA(MEF90HeatXferCtxSetFromOptions(MEF90HeatXferCtx,PETSC_NULL_CHARACTER,MEF90HeatXferDefaultGlobalOptions,MEF90HeatXferDefaultCellSetOptions,MEF90HeatXferDefaultVertexSetOptions,ierr))

    Allocate(nodalVarName(numNodalVar))
    Allocate(cellVarName(numCellVar))
    Allocate(gVarName(numGVar))
    nodalVarName = ["Temperature        "]
    cellVarName  = ["ExternalTemperature","Flux               ","BoundaryFlux       "]
    PetscCallA(MEF90EXOFormat(MEF90Ctx%resultViewer,gVarName,cellVarName,nodalVarName,time,ierr))
    DeAllocate(nodalVarName)
    DeAllocate(cellVarName)
    DeAllocate(gVarName)
    PetscCallA(DMDestroy(dm,ierr))

    !!! Analysis loop:
    Do step = 1, size(time)
        Write(IOBuffer,'("Step: ",I4," Analysis time: ",ES12.5,"\n")') step,time(step)
        PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
        PetscCallA(VecSet(MEF90HeatXferCtx%temperatureLocal,time(step),ierr))
        PetscCallA(MEF90HeatXferUpdateTransients(MEF90HeatXferCtx,step,time(step),ierr))
        PetscCallA(MEF90HeatXferViewEXO(MEF90HeatXferCtx,step,ierr))
    End Do ! step

    PetscCallA(MEF90HeatXferCtxDestroy(MEF90HeatXferCtx,ierr))
    PetscCallA(MEF90CtxDestroy(MEF90Ctx,ierr))
    PetscCallA(MEF90Finalize(ierr))
    PetscCallA(PetscFinalize(ierr))
End Program  TestHeatXferCtx