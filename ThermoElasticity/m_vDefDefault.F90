#include "../MEF90/mef90.inc"
module m_vDefDefault
#include "petsc/finclude/petsc.h"
   use petsc
   use m_MEF90
   use m_MEF90_DefMechCtx
   use m_MEF90_DefMechAT, only: MEF90DefMech_damageTypeAT1, MEF90DefMech_damageTypeAT1exp, MEF90DefMech_damageTypeAT2, MEF90DefMech_damageTypeKKL
   use m_MEF90_HeatXferCtx
   implicit none(type, external)

!!! Default values of the contexts
   type(MEF90CtxGlobalOptions_Type), parameter         :: MEF90CtxDefaultGlobalOptions = MEF90CtxGlobalOptions_Type( &
                                                          0, & ! verbose
                                                          PETSC_FALSE, & ! dryrun
                                                          MEF90TimeInterpolation_linear, & ! timeInterpolation
                                                          0.0_kr, & ! timeMin
                                                          1.0_kr, & ! timeMax
                                                          11_ki, & ! timeNumStep
                                                          0_ki, & ! timeSkip
                                                          1_ki, & ! numCycle
                                                          MEF90ElementFamilyLagrange, & ! elementFamily
                                                          1_ki)                            ! elementOrder

   type(MEF90DefMechGlobalOptions_Type), parameter     :: DefMechDefaultGlobalOptions = MEF90DefMechGlobalOptions_Type( &
                                                          MEF90DefMech_TimeSteppingTypeQuasiStatic, & ! timeSteppingType
                                                          MEF90DefMech_SolverTypeAltMin, & ! solverType
                                                          MEF90DefMech_DamageSolverTypeSNES, & ! damageSolverType
                                                          MEF90Scaling_Linear, & ! boundaryDisplacementScaling
                                                          MEF90Scaling_CST, & ! displacementLowerBoundScaling
                                                          MEF90Scaling_CST, & ! displacementUpperBoundScaling
                                                          MEF90Scaling_CST, & ! cohesiveDisplacementScaling
                                                          MEF90Scaling_CST, & ! boundaryDamageScaling
                                                          MEF90Scaling_Linear, & ! bodyForceScaling
                                                          MEF90Scaling_Linear, & ! boundaryForceScaling
                                                          MEF90Scaling_Linear, & ! pressureForceScaling
                                                          MEF90Scaling_Linear, & ! crackPressureScaling
                                                          1e-3, & ! damage_atol
                                                          1000_ki, & ! maxit
                                                          10_ki, & ! PCLag
                                                          1.0_kr, & ! SOROmega
                                                          0., & ! irrevThres
                                                          MEF90DefMech_BTTypeNULL, & ! BTType
                                                          -1_ki, & ! BTInt
                                                          -1_ki, & ! BTScope
                                                          1.0e-2, & ! BTTol
                                                          1.0e-4, & ! plasticStrainAtol
                                                          1.0e-3, & ! InjectedVolumeAtol
                                                          0.0_kr, & ! dampingCoefficientDisplacement
                                                          0.0_kr, & ! dampingCoefficientDamage
                                                          PETSC_FALSE, & ! temperatureExport
                                                          PETSC_TRUE, & ! displacementExport
                                                          PETSC_TRUE, & ! damageExport
                                                          PETSC_TRUE, & ! stressExport
                                                          PETSC_FALSE, & ! plasticStrainExport
                                                          PETSC_FALSE)               ! cumulatedPlasticDissipationExport

   type(MEF90DefMechCellSetOptions_Type), parameter    :: DefMechDefaultCellSetOptions = MEF90DefMechCellSetOptions_Type( &
                                                          [0.0_kr, 0.0_kr, 0.0_kr], & ! bodyForce
                                                          0.0_kr, & ! crackPressure
                                                          MEF90DefMech_damageTypeAT1, & ! damageType
                                                          MEF90DefMech_plasticityTypeNone, & ! plasticityType
                                                          MEF90DefMech_unilateralContactTypeNone, & ! unilateralContactType
                                                          1.0e-5, & ! unilateralContactHydrostaticDeviatoricGamma
                                                          PETSC_FALSE, & ! unilateralContactHybrid
                                                          1.0_kr, & ! damageATLinSoftk
                                                          1.25_kr, & ! damageAT1expb
                                                          MEF90DefMech_drivingForceTypeNone, & ! drivingForceType
                                                          [0.0_kr, 0.0_kr, 0.0_kr], & ! cohesiveDisplacement
                                                          [PETSC_FALSE, PETSC_FALSE, PETSC_FALSE], & ! hasDisplacementBC
                                                          [0.0_kr, 0.0_kr, 0.0_kr], & ! boundaryDisplacement
                                                          [MEF90NINFINITY, MEF90NINFINITY, MEF90NINFINITY], & ! displacementLowerBound
                                                          [MEF90INFINITY, MEF90INFINITY, MEF90INFINITY], & ! displacementUpperBound
                                                          PETSC_FALSE, & ! hasDamageBC
                                                          0.0_kr, & ! boundaryDamage
                                                          PETSC_FALSE, & ! crackVolumeControlled
                                                          PETSC_FALSE)                                         ! workControlled

   type(MEF90DefMechFaceSetOptions_Type), parameter    :: DefMechDefaultFaceSetOptions = MEF90DefMechFaceSetOptions_Type( &
                                                          [0.0_kr, 0.0_kr, 0.0_kr], & ! boundaryForce
                                                          0.0_kr, & ! pressureForce
                                                          [PETSC_FALSE, PETSC_FALSE, PETSC_FALSE], & ! hasDisplacementBC
                                                          [0.0_kr, 0.0_kr, 0.0_kr], & ! boundaryDisplacement
                                                          [MEF90NINFINITY, MEF90NINFINITY, MEF90NINFINITY], & ! displacementLowerBound
                                                          [MEF90INFINITY, MEF90INFINITY, MEF90INFINITY], & ! displacementUpperBound
                                                          PETSC_FALSE, & ! hasDamageBC
                                                          0.0_kr)                                              ! boundaryDamage

   type(MEF90DefMechVertexSetOptions_Type), parameter  :: DefMechDefaultVertexSetOptions = MEF90DefMechVertexSetOptions_Type( &
                                                          [PETSC_FALSE, PETSC_FALSE, PETSC_FALSE], & ! hasDisplacementBC
                                                          [0.0_kr, 0.0_kr, 0.0_kr], & ! boundaryDisplacement
                                                          [MEF90NINFINITY, MEF90NINFINITY, MEF90NINFINITY], & ! displacementLowerBound
                                                          [MEF90INFINITY, MEF90INFINITY, MEF90INFINITY], & ! displacementUpperBound
                                                          PETSC_FALSE, & ! hasDamageBC
                                                          0.0_kr)                                              ! boundaryDamage

   type(MEF90HeatXferGlobalOptions_Type), parameter    :: HeatXferDefaultGlobalOptions = MEF90HeatXferGlobalOptions_Type( &
                                                          MEF90HeatXFer_timeSteppingTypeSteadyState, & ! timeSteppingType
                                                          PETSC_FALSE, & ! addNullSpace
                                                          0.0_kr, & ! initialTemperature
                                                          MEF90Scaling_Linear, & ! boundaryTempScaling
                                                          MEF90Scaling_Linear, & ! externalTempScaling
                                                          MEF90Scaling_Linear, & ! fluxScaling
                                                          MEF90Scaling_Linear, & ! boundaryFluxScaling
                                                          PETSC_TRUE)            ! temperatureExport

   type(MEF90HeatXferCellSetOptions_Type), parameter   :: HeatXferDefaultCellSetOptions = MEF90HeatXferCellSetOptions_Type( &
                                                          0.0_kr, & ! flux
                                                          PETSC_FALSE, & ! hasTemperatureBC
                                                          0.0_kr, & ! boundaryTemperature
                                                          [0.0_kr, 0.0_kr, 0.0_kr]) ! advectionVector

   type(MEF90HeatXferFaceSetOptions_Type), parameter   :: HeatXferDefaultFaceSetOptions = MEF90HeatXferFaceSetOptions_Type( &
                                                          0.0_kr, & ! boundaryFlux
                                                          0.0_kr, & ! surfaceThermalConductivity
                                                          0.0_kr, & ! externalTemp
                                                          PETSC_FALSE, & ! hasTemperatureBC
                                                          0.0_kr)          ! boundaryTemperature

   type(MEF90HeatXferVertexSetOptions_Type), parameter ::HeatXferDefaultVertexSetOptions = MEF90HeatXferVertexSetOptions_Type( &
                                                          PETSC_FALSE, & ! hasTemperatureBC
                                                          0.0_kr)          ! boundaryTemperature
end module m_vDefDefault
