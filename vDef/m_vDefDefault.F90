#include "../MEF90/mef90.inc"
module m_vDefDefault
#include "petsc/finclude/petsc.h"
   use petsc
   use m_MEF90
   use m_MEF90_DefMech_class
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

   type(MEF90HeatXferGlobalOptions_Type), parameter    :: HeatXferDefaultGlobalOptions = MEF90HeatXferGlobalOptions_Type( &
                                                          MEF90HeatXFer_timeSteppingTypeSteadyState, & ! timeSteppingType
                                                          PETSC_FALSE, & ! addNullSpace
                                                          0.0_kr, & ! initialTemperature
                                                          MEF90Scaling_Linear, & ! boundaryTempScaling
                                                          MEF90Scaling_Linear, & ! externalTempScaling
                                                          MEF90Scaling_Linear, & ! fluxScaling
                                                          MEF90Scaling_Linear, & ! boundaryFluxScaling
                                                          PETSC_FALSE)           ! temperatureExport

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
