Module m_MEF90_DiffusionInterface
#include "finclude/petscdef.h"
   Use m_MEF90_LinAlg
   Use m_MEF90_Parameters
   Use m_MEF90_Elements
   Use m_MEF90_Utils
   Use m_MEF90_DiffusionImplementation_2D, &
         MEF90DiffusionBilinearFormSet_2D => MEF90DiffusionBilinearFormSet, &
         MEF90DiffusionEnergySet_2D => MEF90DiffusionEnergySet, &
         MEF90DiffusionOperatorAddTransientTermSet_2D => MEF90DiffusionOperatorAddTransientTermSet, &
         MEF90DiffusionOperatorSet_2D => MEF90DiffusionOperatorSet, &
         MEF90DiffusionRHSSetVertex_2D => MEF90DiffusionRHSSetVertex, &
         MEF90DiffusionRHSSetCell_2D => MEF90DiffusionRHSSetCell, &
         MEF90DiffusionRHSSetCst_2D => MEF90DiffusionRHSSetCst, &
         MEF90DiffusionWorkSetVertex_2D => MEF90DiffusionWorkSetVertex, &
         MEF90DiffusionWorkSetCell_2D => MEF90DiffusionWorkSetCell, &
         MEF90DiffusionWorkSetCst_2D => MEF90DiffusionWorkSetCst
   Use m_MEF90_DiffusionImplementation_3D, &
         MEF90DiffusionBilinearFormSet_3D => MEF90DiffusionBilinearFormSet, &
         MEF90DiffusionEnergySet_3D => MEF90DiffusionEnergySet, &
         MEF90DiffusionOperatorAddTransientTermSet_3D => MEF90DiffusionOperatorAddTransientTermSet, &
         MEF90DiffusionOperatorSet_3D => MEF90DiffusionOperatorSet, &
         MEF90DiffusionRHSSetVertex_3D => MEF90DiffusionRHSSetVertex, &
         MEF90DiffusionRHSSetCell_3D => MEF90DiffusionRHSSetCell, &
         MEF90DiffusionRHSSetCst_3D => MEF90DiffusionRHSSetCst, &
         MEF90DiffusionWorkSetVertex_3D => MEF90DiffusionWorkSetVertex, &
         MEF90DiffusionWorkSetCell_3D => MEF90DiffusionWorkSetCell, &
         MEF90DiffusionWorkSetCst_3D => MEF90DiffusionWorkSetCst
   Use petsc

   IMPLICIT NONE

   Private   
   Public :: MEF90DiffusionBilinearFormSet
   Public :: MEF90DiffusionEnergySet
   Public :: MEF90DiffusionOperatorAddTransientTermSet
   Public :: MEF90DiffusionOperatorSet
   Public :: MEF90DiffusionRHSSetVertex
   Public :: MEF90DiffusionRHSSetCell
   Public :: MEF90DiffusionRHSSetCst
   Public :: MEF90DiffusionWorkSetVertex
   Public :: MEF90DiffusionWorkSetCell
   Public :: MEF90DiffusionWorkSetCst
   
   Interface MEF90DiffusionBilinearFormSet
      Module Procedure MEF90DiffusionBilinearFormSet_2D, MEF90DiffusionBilinearFormSet_3D
   End Interface MEF90DiffusionBilinearFormSet
   
   Interface MEF90DiffusionEnergySet
      Module procedure MEF90DiffusionEnergySet_2D,MEF90DiffusionEnergySet_3D
   End Interface MEF90DiffusionEnergySet

   Interface MEF90DiffusionOperatorAddTransientTermSet
      Module procedure MEF90DiffusionOperatorAddTransientTermSet_2D,MEF90DiffusionOperatorAddTransientTermSet_3D
   End Interface MEF90DiffusionOperatorAddTransientTermSet
   
   Interface MEF90DiffusionOperatorSet
      Module procedure MEF90DiffusionOperatorSet_2D,MEF90DiffusionOperatorSet_3D
   End Interface MEF90DiffusionOperatorSet
   
   Interface MEF90DiffusionRHSSetVertex
      Module Procedure MEF90DiffusionRHSSetVertex_2D,MEF90DiffusionRHSSetVertex_3D
   End Interface MEF90DiffusionRHSSetVertex

   Interface MEF90DiffusionRHSSetCell
      Module Procedure MEF90DiffusionRHSSetCell_2D,MEF90DiffusionRHSSetCell_3D
   End Interface MEF90DiffusionRHSSetCell

   Interface MEF90DiffusionRHSSetCst
      Module Procedure MEF90DiffusionRHSSetCst_2D,MEF90DiffusionRHSSetCst_3D
   End Interface MEF90DiffusionRHSSetCst

   Interface MEF90DiffusionWorkSetVertex
      Module Procedure MEF90DiffusionWorkSetVertex_2D,MEF90DiffusionWorkSetVertex_3D
   End Interface MEF90DiffusionWorkSetVertex

   Interface MEF90DiffusionWorkSetCell
      Module Procedure MEF90DiffusionWorkSetCell_2D,MEF90DiffusionWorkSetCell_3D
   End Interface MEF90DiffusionWorkSetCell

   Interface MEF90DiffusionWorkSetCst
      Module Procedure MEF90DiffusionWorkSetCst_2D,MEF90DiffusionWorkSetCst_3D
   End Interface MEF90DiffusionWorkSetCst
End Module m_MEF90_DiffusionInterface
