#include "../MEF90/mef90.inc"
module m_vDefDefault
#include "petsc/finclude/petsc.h"
   use petsc
   use m_MEF90
   use m_MEF90_DefMech_class
   use m_MEF90_HeatXfer_class
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
end module m_vDefDefault
