Module m_MEF_ElasticityInterface
#include "finclude/petscdef.h"
   Use m_MEF_LinAlg
   Use m_MEF_Parameters
   Use m_MEF_Elements
   Use m_MEF_Utils
   Use m_MEF_ElasticityImplementation_2D, &
         ElasticityOperatorSet_2D => ElasticityOperatorSet, &
         ElasticityBilinearFormSet_2D => ElasticityBilinearFormSet, &
         ElasticityEnergySet_2D => ElasticityEnergySet, &
         ElasticityOperatorAddTransientTermSet_2D => ElasticityOperatorAddTransientTermSet, &
         ElasticityInelasticStrainRHSSetCst_2D => ElasticityInelasticStrainRHSSetCst, &
         ElasticityInelasticStrainRHSSetCell_2D => ElasticityInelasticStrainRHSSetCell, &
         ElasticityInelasticStrainRHSSetCell2_2D => ElasticityInelasticStrainRHSSetCell2, &
         ElasticityInelasticStrainRHSSetVertex2_2D => ElasticityInelasticStrainRHSSetVertex2, &
         ElasticityForceRHSSetCst_2D => ElasticityForceRHSSetCst, &
         ElasticityForceRHSSetCell_2D => ElasticityForceRHSSetCell, &
         ElasticityForceRHSSetVertex_2D => ElasticityForceRHSSetVertex, &
         ElasticityPressureForceRHSSetCst_2D => ElasticityPressureForceRHSSetCst, &
         ElasticityPressureForceRHSSetCell_2D => ElasticityPressureForceRHSSetCell, &
         ElasticityPressureForceRHSSetVertex_2D => ElasticityPressureForceRHSSetVertex, &
         ElasticityWorkSetCst_2D => ElasticityWorkSetCst, &
         ElasticityWorkSetCell_2D => ElasticityWorkSetCell, &
         ElasticityWorkSetVertex_2D => ElasticityWorkSetVertex
   Use m_MEF_ElasticityImplementation_3D, &
         ElasticityOperatorSet_3D => ElasticityOperatorSet, &
         ElasticityBilinearFormSet_3D => ElasticityBilinearFormSet, &
         ElasticityEnergySet_3D => ElasticityEnergySet, &
         ElasticityOperatorAddTransientTermSet_3D => ElasticityOperatorAddTransientTermSet, &
         ElasticityInelasticStrainRHSSetCst_3D => ElasticityInelasticStrainRHSSetCst, &
         ElasticityInelasticStrainRHSSetCell_3D => ElasticityInelasticStrainRHSSetCell, &
         ElasticityInelasticStrainRHSSetCell2_3D => ElasticityInelasticStrainRHSSetCell2, &
         ElasticityInelasticStrainRHSSetVertex2_3D => ElasticityInelasticStrainRHSSetVertex2, &
         ElasticityForceRHSSetCst_3D => ElasticityForceRHSSetCst, &
         ElasticityForceRHSSetCell_3D => ElasticityForceRHSSetCell, &
         ElasticityForceRHSSetVertex_3D => ElasticityForceRHSSetVertex, &
         ElasticityPressureForceRHSSetCst_3D => ElasticityPressureForceRHSSetCst, &
         ElasticityPressureForceRHSSetCell_3D => ElasticityPressureForceRHSSetCell, &
         ElasticityPressureForceRHSSetVertex_3D => ElasticityPressureForceRHSSetVertex, &
         ElasticityWorkSetCst_3D => ElasticityWorkSetCst, &
         ElasticityWorkSetCell_3D => ElasticityWorkSetCell, &
         ElasticityWorkSetVertex_3D => ElasticityWorkSetVertex
   Use petsc

   IMPLICIT NONE

   Private   
   Public :: MEF90ElasticityOperatorSet
   Public :: MEF90ElasticityBilinearFormSet
   Public :: MEF90ElasticityEnergySet
   Public :: MEF90ElasticityOperatorAddTransientTermSet
   Public :: MEF90ElasticityInelasticStrainRHSSetCst
   Public :: MEF90ElasticityInelasticStrainRHSSetCell
   Public :: MEF90ElasticityInelasticStrainRHSSetVertex
   Public :: MEF90ElasticityForceRHSSetCst
   Public :: MEF90ElasticityForceRHSSetCell
   Public :: MEF90ElasticityForceRHSSetVertex
   Public :: MEF90ElasticityPressureForceRHSSetCst
   Public :: MEF90ElasticityPressureForceRHSSetCell
   Public :: MEF90ElasticityPressureForceRHSSetVertex
   Public :: MEF90ElasticityWorkSetCst
   Public :: MEF90ElasticityWorkSetCell
   Public :: MEF90ElasticityWorkSetVertex
   
   Interface MEF90ElasticityBilinearFormSet
      Module Procedure ElasticityBilinearFormSet_2D, ElasticityBilinearFormSet_3D
   End Interface MEF90ElasticityBilinearFormSet
   
   Interface MEF90ElasticityEnergySet
      Module procedure ElasticityEnergySet_2D,ElasticityEnergySet_3D
   End Interface MEF90ElasticityEnergySet

   Interface MEF90ElasticityOperatorAddTransientTermSet
      Module procedure ElasticityOperatorAddTransientTermSet_2D,ElasticityOperatorAddTransientTermSet_3D
   End Interface MEF90ElasticityOperatorAddTransientTermSet
   
   Interface MEF90ElasticityOperatorSet
      Module procedure ElasticityOperatorSet_2D,ElasticityOperatorSet_3D
   End Interface MEF90ElasticityOperatorSet
   
   Interface MEF90ElasticityInelasticStrainRHSSetCst
      Module Procedure ElasticityInelasticStrainRHSSetCst_2D, ElasticityInelasticStrainRHSSetCst_3D
   End Interface MEF90ElasticityInelasticStrainRHSSetCst

   Interface MEF90ElasticityInelasticStrainRHSSetCell
      Module Procedure ElasticityInelasticStrainRHSSetCell_2D, ElasticityInelasticStrainRHSSetCell_3D, &
                       ElasticityInelasticStrainRHSSetCell2_2D, ElasticityInelasticStrainRHSSetCell2_3D
   End Interface MEF90ElasticityInelasticStrainRHSSetCell

   Interface MEF90ElasticityInelasticStrainRHSSetVertex
      Module Procedure ElasticityInelasticStrainRHSSetVertex2_2D, ElasticityInelasticStrainRHSSetVertex2_3D
   End Interface MEF90ElasticityInelasticStrainRHSSetVertex

   Interface MEF90ElasticityForceRHSSetCst
      Module procedure ElasticityForceRHSSetCst_2D, ElasticityForceRHSSetCst_3D
   End Interface MEF90ElasticityForceRHSSetCst

   Interface MEF90ElasticityForceRHSSetCell
      Module procedure ElasticityForceRHSSetCell_2D, ElasticityForceRHSSetCell_3D
   End Interface MEF90ElasticityForceRHSSetCell

   Interface MEF90ElasticityForceRHSSetVertex
      Module procedure ElasticityForceRHSSetVertex_2D, ElasticityForceRHSSetVertex_3D
   End Interface MEF90ElasticityForceRHSSetVertex

   Interface MEF90ElasticityPressureForceRHSSetCst
      Module procedure ElasticityPressureForceRHSSetCst_2D, ElasticityPressureForceRHSSetCst_3D
   End Interface MEF90ElasticityPressureForceRHSSetCst

   Interface MEF90ElasticityPressureForceRHSSetCell
      Module procedure ElasticityPressureForceRHSSetCell_2D, ElasticityPressureForceRHSSetCell_3D
   End Interface MEF90ElasticityPressureForceRHSSetCell

   Interface MEF90ElasticityPressureForceRHSSetVertex
      Module procedure ElasticityPressureForceRHSSetVertex_2D, ElasticityPressureForceRHSSetVertex_3D
   End Interface MEF90ElasticityPressureForceRHSSetVertex

   Interface MEF90ElasticityWorkSetCst
      Module procedure ElasticityWorkSetCst_2D, ElasticityWorkSetCst_3D
   End Interface MEF90ElasticityWorkSetCst

   Interface MEF90ElasticityWorkSetCell
      Module procedure ElasticityWorkSetCell_2D, ElasticityWorkSetCell_3D
   End Interface MEF90ElasticityWorkSetCell

   Interface MEF90ElasticityWorkSetVertex
      Module procedure ElasticityWorkSetVertex_2D, ElasticityWorkSetVertex_3D
   End Interface MEF90ElasticityWorkSetVertex
End Module m_MEF_ElasticityInterface
