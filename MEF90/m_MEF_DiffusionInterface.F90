Module m_MEF_DiffusionInterface
#include "finclude/petscdef.h"
   Use m_MEF_LinAlg
   Use m_MEF_Parameters
   Use m_MEF_Elements
   Use m_MEF_Sieve
   Use m_MEF_Utils
   Use m_MEF_DiffusionImplementation_2D, &
         DiffusionBilinearFormSet_2D => DiffusionBilinearFormSet, &
         DiffusionEnergySet_2D => DiffusionEnergySet, &
         DiffusionOperatorAddTransientTermSet_2D => DiffusionOperatorAddTransientTermSet, &
         DiffusionOperatorSet_2D => DiffusionOperatorSet, &
         DiffusionRHSSetVertex_2D => DiffusionRHSSetVertex, &
         DiffusionRHSSetCell_2D => DiffusionRHSSetCell, &
         DiffusionRHSSetCst_2D => DiffusionRHSSetCst, &
         DiffusionWorkSetVertex_2D => DiffusionWorkSetVertex, &
         DiffusionWorkSetCell_2D => DiffusionWorkSetCell, &
         DiffusionWorkSetCst_2D => DiffusionWorkSetCst
   Use m_MEF_DiffusionImplementation_3D, &
         DiffusionBilinearFormSet_3D => DiffusionBilinearFormSet, &
         DiffusionEnergySet_3D => DiffusionEnergySet, &
         DiffusionOperatorAddTransientTermSet_3D => DiffusionOperatorAddTransientTermSet, &
         DiffusionOperatorSet_3D => DiffusionOperatorSet, &
         DiffusionRHSSetVertex_3D => DiffusionRHSSetVertex, &
         DiffusionRHSSetCell_3D => DiffusionRHSSetCell, &
         DiffusionRHSSetCst_3D => DiffusionRHSSetCst, &
         DiffusionWorkSetVertex_3D => DiffusionWorkSetVertex, &
         DiffusionWorkSetCell_3D => DiffusionWorkSetCell, &
         DiffusionWorkSetCst_3D => DiffusionWorkSetCst
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
      Module Procedure DiffusionBilinearFormSet_2D, DiffusionBilinearFormSet_3D
   End Interface MEF90DiffusionBilinearFormSet
   
   Interface MEF90DiffusionEnergySet
      Module procedure DiffusionEnergySet_2D,DiffusionEnergySet_3D
   End Interface MEF90DiffusionEnergySet

   Interface MEF90DiffusionOperatorAddTransientTermSet
      Module procedure DiffusionOperatorAddTransientTermSet_2D,DiffusionOperatorAddTransientTermSet_3D
   End Interface MEF90DiffusionOperatorAddTransientTermSet
   
   Interface MEF90DiffusionOperatorSet
      Module procedure DiffusionOperatorSet_2D,DiffusionOperatorSet_3D
   End Interface MEF90DiffusionOperatorSet
   
   Interface MEF90DiffusionRHSSetVertex
      Module Procedure DiffusionRHSSetVertex_2D,DiffusionRHSSetVertex_3D
   End Interface MEF90DiffusionRHSSetVertex

   Interface MEF90DiffusionRHSSetCell
      Module Procedure DiffusionRHSSetCell_2D,DiffusionRHSSetCell_3D
   End Interface MEF90DiffusionRHSSetCell

   Interface MEF90DiffusionRHSSetCst
      Module Procedure DiffusionRHSSetCst_2D,DiffusionRHSSetCst_3D
   End Interface MEF90DiffusionRHSSetCst

   Interface MEF90DiffusionWorkSetVertex
      Module Procedure DiffusionWorkSetVertex_2D,DiffusionWorkSetVertex_3D
   End Interface MEF90DiffusionWorkSetVertex

   Interface MEF90DiffusionWorkSetCell
      Module Procedure DiffusionWorkSetCell_2D,DiffusionWorkSetCell_3D
   End Interface MEF90DiffusionWorkSetCell

   Interface MEF90DiffusionWorkSetCst
      Module Procedure DiffusionWorkSetCst_2D,DiffusionWorkSetCst_3D
   End Interface MEF90DiffusionWorkSetCst
End Module m_MEF_DiffusionInterface
