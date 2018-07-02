#include "../MEF90/mef90.inc"
module m_vDefDefault
#include <finclude/petscdef.h>
   Use petsc
   Use m_MEF90
   Use m_MEF90_DefMechCtx
   Use m_MEF90_HeatXferCtx
   Implicit NONE   


!!! Default values of the contexts
   Type(MEF90CtxGlobalOptions_Type),Parameter         :: vDefDefaultGlobalOptions = MEF90CtxGlobalOptions_Type( &
                                                         1,                             & ! verbose
                                                         PETSC_FALSE,                   & ! validate
                                                         MEF90TimeInterpolation_linear, & ! timeInterpolation
                                                         0.0_Kr,                        & ! timeMin
                                                         1.0_Kr,                        & ! timeMax
                                                         11,                            & ! timeNumStep
                                                         0,                             & ! timeSkip
                                                         1.0_Kr,                        & ! frequency
                                                         MEF90FileFormat_EXOSingle)       ! fileFormat

   Type(MEF90DefMechGlobalOptions_Type),Parameter     :: vDefDefMechDefaultGlobalOptions2D = MEF90DefMechGlobalOptions_Type( &
                                                         MEF90DefMech_TimeSteppingTypeQuasiStatic, & ! solverType
                                                         MEF90DefMech_SolverTypeAltMin,            & ! timeSteppingType
                                                         PETSC_TRUE,              & ! disp_addNullSpace
                                                         3,                       & ! DisplacementOffset
                                                         2,                       & ! DamageOffset
                                                         3,                       & ! boundaryDisplacementOffset
                                                         0,                       & ! boundaryDamageOffset
                                                         1,                       & ! temperatureOffset
                                                         4,                       & ! ForceOffset
                                                         3,                       & ! pressureForceOffset
                                                         0,                       & ! CrackPressureOffset
                                                         0,                       & ! plasticStrainOffset
                                                         6,                       & ! StressOffset
                                                         MEF90Scaling_Linear,     & ! boundaryDisplacementScaling
                                                         MEF90Scaling_CST,        & ! boundaryDamageScaling
                                                         MEF90Scaling_Linear,     & ! ForceScaling
                                                         MEF90Scaling_Linear,     & ! pressureForceScaling
                                                         MEF90Scaling_Linear,     & ! CrackPressureScaling
                                                         1e-3,                    & ! damage_atol
                                                         1000,                    & ! maxit
                                                         10,                      & ! PCLag
                                                         1.0_Kr,                  & ! SOROmega
                                                         0.,                      & ! irrevThres 
                                                         MEF90DefMech_BTTypeNULL, & ! BTType
                                                         -1,                      & ! BTInt
                                                         -1,                      & ! BTScope
                                                         1.0e-2,                  & ! BTTol
                                                         1.0e-4,                  & ! plasticStrainAtol
                                                         1,                       & ! cumulatedDissipatedPlasticEnergyOffset
                                                         1.0e-3,                  & ! InjectedVolumeAtol
                                                         0.0_Kr,                  & ! dampingCoefficientDisplacement
                                                         0.0_Kr)                    ! dampingCoefficientDamage

   Type(MEF90DefMechGlobalOptions_Type),Parameter     :: vDefDefMechDefaultGlobalOptions3D = MEF90DefMechGlobalOptions_Type( &
                                                         MEF90DefMech_SolverTypeAltMin,            & ! timeSteppingType
                                                         MEF90DefMech_TimeSteppingTypeQuasiStatic, & ! mode
                                                         PETSC_TRUE,              & ! disp_addNullSpace
                                                         3,                       & ! DisplacementOffset
                                                         2,                       & ! DamageOffset
                                                         3,                       & ! boundaryDisplacementOffset
                                                         0,                       & ! boundaryDamageOffset
                                                         1,                       & ! temperatureOffset
                                                         4,                       & ! ForceOffset
                                                         3,                       & ! pressureForceOffset
                                                         0,                       & ! CrackPressureOffset
                                                         0,                       & ! plasticStrainOffset
                                                         7,                       & ! StressOffset
                                                         MEF90Scaling_Linear,     & ! boundaryDisplacementScaling
                                                         MEF90Scaling_CST,        & ! boundaryDamageScaling
                                                         MEF90Scaling_Linear,     & ! ForceScaling
                                                         MEF90Scaling_Linear,     & ! pressureForceScaling
                                                         MEF90Scaling_Linear,     & ! CrackPressureScaling
                                                         1e-4,                    & ! damage_atol
                                                         1000,                    & ! maxit
                                                         10,                      & ! PCLag
                                                         1.0_Kr,                  & ! SOROmega
                                                         0.,                      & ! irrevThres 
                                                         MEF90DefMech_BTTypeNULL, & ! BTType
                                                         -1,                      & ! BTInt
                                                         -1,                      & ! BTScope
                                                         1.0e-2,                  & ! BTTol
                                                         1.0e-4,                  & ! plasticStrainAtol
                                                         1,                       & ! cumulatedDissipatedPlasticEnergyOffset
                                                         1.0e-3,                  & ! InjectedVolumeAtol
                                                         0.0_Kr,                  & ! dampingCoefficientDisplacement
                                                         0.0_Kr)                    ! dampingCoefficientDamage

   
   Type(MEF90DefMechCellSetOptions_Type),Parameter    :: vDefDefMechDefaultCellSetOptions = MEF90DefMechCellSetOptions_Type( &
                                                         -1,                                      & ! elemTypeShortIDDispl will be overriden
                                                         -1,                                      & ! elemTypeShortIDDamage will be overriden
                                                         [0.0_Kr,0.0_Kr,0.0_Kr],                  & ! force
                                                         0.0_Kr,                                  & ! pressureForce
                                                         0.0_Kr,                                  & ! CrackPressure
                                                         MEF90DefMech_damageTypeAT1,              & ! damageType
                                                         MEF90DefMech_plasticityTypeNone,         & ! plasticityType
                                                         MEF90DefMech_unilateralContactTypeNone,  & ! unilateralContactType
                                                         [PETSC_FALSE,PETSC_FALSE,PETSC_FALSE],   & ! Has Displacement BC
                                                         [0.0_Kr,0.0_Kr,0.0_Kr],                  & ! boundary Displacement
                                                         PETSC_FALSE,                             & ! Has Damage BC
                                                         PETSC_FALSE,                             & ! CrackVolumeControlled
                                                         PETSC_FALSE,                             & ! WorkControlled
                                                         0._Kr)                                     ! Boundary Damage
   Type(MEF90DefMechVertexSetOptions_Type),Parameter  :: vDefDefMechDefaultVertexSetOptions = MEF90DefMechVertexSetOptions_Type( &
                                                         [PETSC_FALSE,PETSC_FALSE,PETSC_FALSE],   & ! Has Displacement BC
                                                         [0.0_Kr,0.0_Kr,0.0_Kr],                  & ! boundary Displacement
                                                         PETSC_FALSE,                             & ! Has Damage BC
                                                         0.0_Kr)                                    ! boundary Damage

   Type(MEF90HeatXferGlobalOptions_Type),Parameter    :: vDefHeatXferDefaultGlobalOptions = MEF90HeatXferGlobalOptions_Type( &
                                                         MEF90HeatXFer_timeSteppingTypeSteadyState, & ! timeSteppingType
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


end module m_vDefDefault