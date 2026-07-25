Module m_MEF90_HeatXferDefault
#include <petsc/finclude/petsc.h>
    Use m_MEF90_HeatXfer
    implicit none (type, external)   

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
   
End Module m_MEF90_HeatXferDefault
