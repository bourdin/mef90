Module m_MEF90_NormsInterface
#include "petsc/finclude/petsc.h"
    Use m_MEF90_NormsImplementation_MEF90Element2D_Scal,  MEF90_L2DotProductSet_2DScal    => MEF90_L2DotProductSet, &
                                                          MEF90_H1DotProductSet_2DScal    => MEF90_H1DotProductSet, &
                                                          MEF90_L2NormSet_2DScal          => MEF90_L2NormSet 
    Use m_MEF90_NormsImplementation_MEF90Element2D_Vect,  MEF90_L2DotProductSet_2DVect    => MEF90_L2DotProductSet, &
                                                          MEF90_H1DotProductSet_2DVect    => MEF90_H1DotProductSet, &
                                                          MEF90_H1SymDotProductSet_2DVect => MEF90_H1SymDotProductSet, &
                                                          MEF90_L2NormSet_2DVect          => MEF90_L2NormSet
    Use m_MEF90_NormsImplementation_MEF90Element3D_Scal,  MEF90_L2DotProductSet_3DScal    => MEF90_L2DotProductSet, &
                                                          MEF90_H1DotProductSet_3DScal    => MEF90_H1DotProductSet, &
                                                          MEF90_L2NormSet_3DScal          => MEF90_L2NormSet
    Use m_MEF90_NormsImplementation_MEF90Element3D_Vect,  MEF90_L2DotProductSet_3DVect    => MEF90_L2DotProductSet, &
                                                          MEF90_H1DotProductSet_3DVect    => MEF90_H1DotProductSet, &
                                                          MEF90_H1SymDotProductSet_3DVect => MEF90_H1SymDotProductSet, &
                                                          MEF90_L2NormSet_3DVect          => MEF90_L2NormSet

    IMPLICIT NONE

    Private

    Public :: MEF90_L2DotProductSet,MEF90_H1DotProductSet,MEF90_H1symDotProductSet,MEF90_L2NormSet

    Interface MEF90_L2DotProductSet
        Module Procedure MEF90_L2DotProductSet_2DScal, MEF90_L2DotProductSet_2DVect, &
                         MEF90_L2DotProductSet_3DScal, MEF90_L2DotProductSet_3DVect 
    End Interface MEF90_L2DotProductSet

    Interface MEF90_H1DotProductSet
        Module Procedure MEF90_H1DotProductSet_2DScal, MEF90_H1DotProductSet_2DVect, &
                         MEF90_H1DotProductSet_3DScal, MEF90_H1DotProductSet_3DVect 
    End Interface MEF90_H1DotProductSet

    Interface MEF90_H1symDotProductSet
        Module Procedure MEF90_H1symDotProductSet_2DVect, &
                         MEF90_H1symDotProductSet_3DVect 
    End Interface MEF90_H1symDotProductSet

    Interface MEF90_L2NormSet
        Module Procedure MEF90_L2NormSet_2DScal, MEF90_L2NormSet_2DVect, &
                        MEF90_L2NormSet_3DScal, MEF90_L2NormSet_3DVect 
    End Interface MEF90_L2NormSet
End Module m_MEF90_NormsInterface
    
    