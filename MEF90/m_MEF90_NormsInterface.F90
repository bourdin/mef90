Module m_MEF90_NormsInterface
#include "petsc/finclude/petsc.h"
       Use m_MEF90_NormsImplementation_MEF90Element2D_Scal,  MEF90_L2DotProductSet_2DScal => MEF90_L2DotProductSet, &
                                                             MEF90_L2NormSet_2DScal  => MEF90_L2NormSet 
       Use m_MEF90_NormsImplementation_MEF90Element2D_Vect,  MEF90_L2DotProductSet_2DVect => MEF90_L2DotProductSet, &
                                                             MEF90_L2NormSet_2DVect  => MEF90_L2NormSet
       Use m_MEF90_NormsImplementation_MEF90Element3D_Scal,  MEF90_L2DotProductSet_3DScal => MEF90_L2DotProductSet, &
                                                             MEF90_L2NormSet_3DScal  => MEF90_L2NormSet
       Use m_MEF90_NormsImplementation_MEF90Element3D_Vect,  MEF90_L2DotProductSet_3DVect => MEF90_L2DotProductSet, &
                                                             MEF90_L2NormSet_3DVect  => MEF90_L2NormSet
    
       IMPLICIT NONE
    
       Private
    
       Public :: MEF90_L2DotProductSet,MEF90_L2NormSet

       Interface MEF90_L2DotProductSet
          Module Procedure MEF90_L2DotProductSet_2DScal, MEF90_L2DotProductSet_2DVect, &
                           MEF90_L2DotProductSet_3DScal, MEF90_L2DotProductSet_3DVect 
       End Interface MEF90_L2DotProductSet
       Interface MEF90_L2NormSet
          Module Procedure MEF90_L2NormSet_2DScal, MEF90_L2NormSet_2DVect, &
                           MEF90_L2NormSet_3DScal, MEF90_L2NormSet_3DVect 
       End Interface MEF90_L2NormSet
    End Module m_MEF90_NormsInterface
    
    