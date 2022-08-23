Program  TestHeatXferCtx
#include <petsc/finclude/petsc.h>
    Use m_MEF90
    Use m_MEF90_HeatXfer
    Use petsc
    Implicit NONE   
    
    PetscErrorCode                        :: ierr
    Type(MEF90Ctx_Type),target            :: MEF90Ctx
    Type(MEF90CtxGlobalOptions_Type)      :: MEF90GlobalOptions_default
    Type(MEF90HeatXferCtx_Type)           :: MEF90HeatXferCtx
    Type(tDM)                             :: dm

    Type(tVec)                            :: temperatureGlobal


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
    
    PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
    
    PetscCallA(MEF90Initialize(ierr))
    PetscCallA(MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90GlobalOptions_default,ierr))

    PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm,MEF90Ctx%geometryfile,PETSC_NULL_CHARACTER,PETSC_TRUE,dm,ierr))
    PetscCallA(DMPlexDistributeSetDefault(dm,PETSC_FALSE,ierr))
    PetscCallA(DMSetFromOptions(dm,ierr))
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-mef90dm_view",ierr))
    distribute: Block 
        Type(tDM),target                    :: dmDist
        PetscInt                            :: ovlp = 0
        If (MEF90Ctx%NumProcs > 1) Then
            PetscCallA(DMPlexDistribute(dm,ovlp,PETSC_NULL_SF,dmDist,ierr))
            PetscCallA(DMDestroy(dm,ierr))
            dm = dmDist
        End If
    End Block distribute
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-mef90dm_view",ierr))

    PetscCallA(MEF90HeatXferCtxCreate(MEF90HeatXferCtx,dm,MEF90Ctx,ierr))
    PetscCallA(MEF90HeatXferCtxSetFromOptions(MEF90HeatXferCtx,PETSC_NULL_CHARACTER,MEF90HeatXferDefaultGlobalOptions,MEF90HeatXferDefaultCellSetOptions,MEF90HeatXferDefaultVertexSetOptions,ierr))

    
    PetscCallA(MEF90HeatXferCtxDestroy(MEF90HeatXferCtx,ierr))
    PetscCallA(MEF90CtxDestroy(MEF90Ctx,ierr))
    PetscCallA(MEF90Finalize(ierr))
    PetscCallA(PetscFinalize(ierr))
End Program  TestHeatXferCtx