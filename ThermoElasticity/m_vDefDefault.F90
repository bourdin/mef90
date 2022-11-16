#include "../MEF90/mef90.inc"
module m_vDefDefault
#include "petsc/finclude/petsc.h"
   Use petsc
   Use m_MEF90
   Use m_MEF90_DefMechCtx
   Use m_MEF90_HeatXferCtx
   Implicit NONE   


!!! Default values of the contexts
   Type(MEF90CtxGlobalOptions_Type),Parameter         :: MEF90CtxDefaultGlobalOptions = MEF90CtxGlobalOptions_Type( &
                                                         0,                             & ! verbose
                                                         PETSC_FALSE,                   & ! dryrun
                                                         MEF90TimeInterpolation_linear, & ! timeInterpolation
                                                         0.0_Kr,                        & ! timeMin
                                                         1.0_Kr,                        & ! timeMax
                                                         11_Ki,                         & ! timeNumStep
                                                         0_Ki,                          & ! timeSkip
                                                         1_Ki,                          & ! numCycle
                                                         MEF90ElementFamilyLagrange,    & ! elementFamily
                                                         1_Ki)                            ! elementOrder

   Type(MEF90DefMechGlobalOptions_Type),Parameter     :: DefMechDefaultGlobalOptions = MEF90DefMechGlobalOptions_Type( &
                                                         MEF90DefMech_TimeSteppingTypeQuasiStatic, & ! solverType
                                                         MEF90DefMech_SolverTypeAltMin,            & ! timeSteppingType
                                                         PETSC_TRUE,              & ! disp_addNullSpace
                                                         MEF90Scaling_Linear,     & ! boundaryDisplacementScaling
                                                         MEF90Scaling_CST,        & ! displacementLowerBoundScaling
                                                         MEF90Scaling_CST,        & ! displacementUpperBoundScaling
                                                         MEF90Scaling_CST,        & ! boundaryDamageScaling
                                                         MEF90Scaling_Linear,     & ! bodyForceScaling
                                                         MEF90Scaling_Linear,     & ! boundaryForceScaling
                                                         MEF90Scaling_Linear,     & ! pressureForceScaling
                                                         MEF90Scaling_Linear,     & ! crackPressureScaling
                                                         1e-3,                    & ! damage_atol
                                                         1000_Ki,                 & ! maxit
                                                         10_Ki,                   & ! PCLag
                                                         1.0_Kr,                  & ! SOROmega
                                                         0.,                      & ! irrevThres 
                                                         MEF90DefMech_BTTypeNULL, & ! BTType
                                                         -1_Ki,                   & ! BTInt
                                                         -1_Ki,                   & ! BTScope
                                                         1.0e-2,                  & ! BTTol
                                                         1.0e-4,                  & ! plasticStrainAtol
                                                         1.0e-3,                  & ! InjectedVolumeAtol
                                                         0.0_Kr,                  & ! dampingCoefficientDisplacement
                                                         0.0_Kr)                    ! dampingCoefficientDamage

   Type(MEF90DefMechCellSetOptions_Type),Parameter    :: DefMechDefaultCellSetOptions = MEF90DefMechCellSetOptions_Type( &
                                                         [0.0_Kr,0.0_Kr,0.0_Kr],                            & ! bodyForce
                                                         0.0_Kr,                                            & ! CrackPressure
                                                         MEF90DefMech_damageTypeAT1,                        & ! damageType
                                                         MEF90DefMech_plasticityTypeNone,                   & ! plasticityType
                                                         MEF90DefMech_unilateralContactTypeNone,            & ! unilateralContactType
                                                         1.0e-5,                                            & ! unilateralContactHydrostatocDeviatoricGamma
                                                         PETSC_FALSE,                                       & ! unilateralContactHybrid
                                                         1.0_Kr,                                            & ! DamageATLinSoftk
                                                         1.25_Kr,                                           & ! DamageAT1expb
                                                         MEF90DefMech_drivingForceTypeNone,                 & ! drivingForceType
                                                         [PETSC_FALSE,PETSC_FALSE,PETSC_FALSE],             & ! Has Displacement BC
                                                         [0.0_Kr,0.0_Kr,0.0_Kr],                            & ! boundary Displacement
                                                         [MEF90NINFINITY,MEF90NINFINITY,MEF90NINFINITY],    & ! displacement Lower Bound
                                                         [MEF90INFINITY,MEF90INFINITY,MEF90INFINITY],       & ! displacement Upper Bound
                                                         PETSC_FALSE,                                       & ! Has Damage BC
                                                         0.0_Kr,                                            & ! boundaryDamage
                                                         PETSC_FALSE,                                       & ! CrackVolumeControlled
                                                         PETSC_FALSE)                                         ! WorkControlled

   Type(MEF90DefMechFaceSetOptions_Type),Parameter    :: DefMechDefaultFaceSetOptions = MEF90DefMechFaceSetOptions_Type( &
                                                         [0.0_Kr,0.0_Kr,0.0_Kr],                            & ! boundaryForce
                                                         0.0_Kr,                                            & ! pressureForce
                                                         [PETSC_FALSE,PETSC_FALSE,PETSC_FALSE],             & ! Has Displacement BC
                                                         [0.0_Kr,0.0_Kr,0.0_Kr],                            & ! boundary Displacement
                                                         [MEF90NINFINITY,MEF90NINFINITY,MEF90NINFINITY],    & ! displacement Lower Bound
                                                         [MEF90INFINITY,MEF90INFINITY,MEF90INFINITY],       & ! displacement Upper Bound
                                                         PETSC_FALSE,                                       & ! Has Damage BC
                                                         0.0_Kr)                                              ! boundaryDamage

   Type(MEF90DefMechVertexSetOptions_Type),Parameter  :: DefMechDefaultVertexSetOptions = MEF90DefMechVertexSetOptions_Type( &
                                                         [PETSC_FALSE,PETSC_FALSE,PETSC_FALSE],             & ! Has Displacement BC
                                                         [0.0_Kr,0.0_Kr,0.0_Kr],                            & ! boundary Displacement
                                                         [MEF90NINFINITY,MEF90NINFINITY,MEF90NINFINITY],    & ! displacement Lower Bound
                                                         [MEF90INFINITY,MEF90INFINITY,MEF90INFINITY],       & ! displacement Upper Bound
                                                         PETSC_FALSE,                                       & ! Has Damage BC
                                                         0.0_Kr)                                              ! boundary Damage


   Type(MEF90HeatXferGlobalOptions_Type),Parameter    :: HeatXferDefaultGlobalOptions = MEF90HeatXferGlobalOptions_Type( &
                                                         MEF90HeatXFer_timeSteppingTypeSteadyState, & ! timeSteppingType
                                                         PETSC_FALSE,         & ! addNullSpace
                                                         0.0_Kr,              & ! initialTemperature
                                                         MEF90Scaling_Linear, & ! boundaryTempScaling
                                                         MEF90Scaling_Linear, & ! externalTempScaling
                                                         MEF90Scaling_Linear, & ! fluxScaling
                                                         MEF90Scaling_Linear)   ! boundaryFluxScaling
 
     Type(MEF90HeatXferCellSetOptions_Type),Parameter :: HeatXferDefaultCellSetOptions = MEF90HeatXferCellSetOptions_Type( &
                                                         0.0_Kr,        & ! flux
                                                         PETSC_FALSE,   & ! Has BC
                                                         0.0_Kr,        & ! boundaryTemperature
                                                         [0.0_Kr,0.0_Kr,0.0_Kr]) ! AdvectionVector
                                                         
     Type(MEF90HeatXferFaceSetOptions_Type),Parameter :: HeatXferDefaultFaceSetOptions = MEF90HeatXferFaceSetOptions_Type( &
                                                         0.0_Kr,        & ! boundaryFlux
                                                         0.0_Kr,        & ! surfaceThermalConductivity
                                                         0.0_Kr,        & ! externalTemp
                                                         PETSC_FALSE,   & ! Has BC
                                                         0.0_Kr)          ! boundaryTemperature
                                                         
     Type(MEF90HeatXferVertexSetOptions_Type),Parameter::HeatXferDefaultVertexSetOptions = MEF90HeatXferVertexSetOptions_Type( &
                                                         PETSC_FALSE,   & ! Has BC
                                                         0.0_Kr)          ! boundaryTemp
end module m_vDefDefault