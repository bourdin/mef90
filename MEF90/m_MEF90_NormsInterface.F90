Module m_MEF90_NormsInterface
#include "petsc/finclude/petsc.h"
    Use m_MEF90_NormsImplementation_MEF90Element2DScal, MEF90L2DotProductSet2DScal    => MEF90L2DotProductSet, &
                                                        MEF90H1DotProductSet2DScal    => MEF90H1DotProductSet, &
                                                        MEF90L2NormSet2DScal          => MEF90L2NormSet 
    Use m_MEF90_NormsImplementation_MEF90Element2DVect, MEF90L2DotProductSet2DVect    => MEF90L2DotProductSet, &
                                                        MEF90H1DotProductSet2DVect    => MEF90H1DotProductSet, &
                                                        MEF90H1SymDotProductSet2DVect => MEF90H1SymDotProductSet, &
                                                        MEF90L2NormSet2DVect          => MEF90L2NormSet
    Use m_MEF90_NormsImplementation_MEF90Element3DScal, MEF90L2DotProductSet3DScal    => MEF90L2DotProductSet, &
                                                        MEF90H1DotProductSet3DScal    => MEF90H1DotProductSet, &
                                                        MEF90L2NormSet3DScal          => MEF90L2NormSet
    Use m_MEF90_NormsImplementation_MEF90Element3DVect, MEF90L2DotProductSet3DVect    => MEF90L2DotProductSet, &
                                                        MEF90H1DotProductSet3DVect    => MEF90H1DotProductSet, &
                                                        MEF90H1SymDotProductSet3DVect => MEF90H1SymDotProductSet, &
                                                        MEF90L2NormSet3DVect          => MEF90L2NormSet

    IMPLICIT NONE

    Private

    Public :: MEF90L2DotProductSet,MEF90H1DotProductSet,MEF90H1symDotProductSet,MEF90L2NormSet

    Interface MEF90L2DotProductSet
        Module Procedure MEF90L2DotProductSet2DScal, MEF90L2DotProductSet2DVect, &
                         MEF90L2DotProductSet3DScal, MEF90L2DotProductSet3DVect 
    End Interface MEF90L2DotProductSet

    Interface MEF90H1DotProductSet
        Module Procedure MEF90H1DotProductSet2DScal, MEF90H1DotProductSet2DVect, &
                         MEF90H1DotProductSet3DScal, MEF90H1DotProductSet3DVect 
    End Interface MEF90H1DotProductSet

    Interface MEF90H1symDotProductSet
        Module Procedure MEF90H1symDotProductSet2DVect, &
                         MEF90H1symDotProductSet3DVect 
    End Interface MEF90H1symDotProductSet

    Interface MEF90L2NormSet
        Module Procedure MEF90L2NormSet2DScal, MEF90L2NormSet2DVect, &
                        MEF90L2NormSet3DScal, MEF90L2NormSet3DVect 
    End Interface MEF90L2NormSet
End Module m_MEF90_NormsInterface
    
    