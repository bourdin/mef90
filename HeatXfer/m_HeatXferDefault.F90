Module m_MEF90_HeatXferDefault
#include <petsc/finclude/petsc.h>
    Use m_MEF90_HeatXfer
    Implicit NONE   

    Type(MEF90CtxGlobalOptions_Type),Parameter         :: MEF90DefaultGlobalOptions = MEF90CtxGlobalOptions_Type( &
                                                            1,                               & ! verbose
                                                            PETSC_FALSE,                     & ! dryrun
                                                            MEF90TimeInterpolation_linear,   & ! timeInterpolation
                                                            0.0_Kr,                          & ! timeMin
                                                            1.0_Kr,                          & ! timeMax
                                                            11,                              & ! timeNumStep
                                                            0,                               & ! timeSkip
                                                            1.0_Kr,                          & ! timeNumCycle
                                                            MEF90ElementFamilyLagrange,      & ! elementFamily
                                                            1)                                 ! elementOrder
     
    Type(MEF90HeatXferGlobalOptions_Type),Parameter    :: MEF90HeatXferDefaultGlobalOptions = MEF90HeatXferGlobalOptions_Type( &
                                                        MEF90HeatXFer_timeSteppingTypeSteadyState, & ! timeSteppingType
                                                        PETSC_FALSE,         & ! addNullSpace
                                                        0.,                  & ! initialTemperature
                                                        MEF90Scaling_Linear, & ! boundaryTempScaling
                                                        MEF90Scaling_Linear, & ! externalTempScaling
                                                        MEF90Scaling_Linear, & ! fluxScaling
                                                        MEF90Scaling_Linear)   ! boundaryFluxScaling

    Type(MEF90HeatXferCellSetOptions_Type),Parameter   :: MEF90HeatXferDefaultCellSetOptions = MEF90HeatXferCellSetOptions_Type( &
                                                        0.0_Kr,        & ! flux
                                                        PETSC_FALSE,   & ! Has BC
                                                        0.0_Kr,        & ! boundaryTemperature
                                                        [0.0_Kr,0.0_Kr,0.0_Kr]) ! AdvectionVector
                                                        
    Type(MEF90HeatXferFaceSetOptions_Type),Parameter   :: MEF90HeatXferDefaultFaceSetOptions = MEF90HeatXferFaceSetOptions_Type( &
                                                        0.0_Kr,        & ! boundaryFlux
                                                        0.0_Kr,        & ! surfaceThermalConductivity
                                                        0.0_Kr,        & ! externalTemp
                                                        PETSC_FALSE,   & ! Has BC
                                                        0.0_Kr)          ! boundaryTemperature
                                                        
    Type(MEF90HeatXferVertexSetOptions_Type),Parameter :: MEF90HeatXferDefaultVertexSetOptions = MEF90HeatXferVertexSetOptions_Type( &
                                                        PETSC_FALSE,   & ! Has BC
                                                        0.0_Kr)          ! boundaryTemp
   
End Module m_MEF90_HeatXferDefault
