Module m_MEF90_ElasticityInterface
#include "finclude/petscdef.h"
   Use m_MEF90_LinAlg
   Use m_MEF90_Parameters
   Use m_MEF90_Elements
   Use m_MEF90_Utils
   Use m_MEF90_ElasticityImplementation_2D, &
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
         ElasticityWorkSetVertex_2D => ElasticityWorkSetVertex, &
         ElasticityCohesiveEnergySet_2D => ElasticityCohesiveEnergySet, &
         !ElasticityPressureWorkSetCst_2D => ElasticityPressureWorkSetCst, &
         ElasticityPressureWorkSetCell_2D => ElasticityPressureWorkSetCell, &
         !ElasticityWorkPressureSetVertex_2D => ElasticityPressureWorkSetVertex
         ElasticityStressSet_2D => ElasticityStressSet, &
         InelasticStrainSet_2D => InelasticStrainSet, &
         PlasticityEnergySet_2D => PlasticityEnergySet



   Use m_MEF90_ElasticityImplementation_3D, &
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
         ElasticityWorkSetVertex_3D => ElasticityWorkSetVertex, &
         ElasticityCohesiveEnergySet_3D => ElasticityCohesiveEnergySet, &
         !ElasticityPressureWorkSetCst_3D => ElasticityPressureWorkSetCst, &
         ElasticityPressureWorkSetCell_3D => ElasticityPressureWorkSetCell, &
         !ElasticityWorkPressureSetVertex_3D => ElasticityPressureWorkSetVertex
         ElasticityStressSet_3D => ElasticityStressSet, &
         InelasticStrainSet_3D => InelasticStrainSet, &
         PlasticityEnergySet_3D => PlasticityEnergySet
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
   Public :: MEF90ElasticityCohesiveEnergySet
   !Public :: MEF90ElasticityPressureWorkSetCst
   Public :: MEF90ElasticityPressureWorkSetCell
   !Public :: MEF90ElasticityPressureWorkSetVertex
   Public :: MEF90ElasticityStressSet
   Public :: MEF90InelasticStrainSet
   Public :: MEF90PlasticityEnergySet
   

   
   Interface MEF90ElasticityBilinearFormSet
      Module Procedure ElasticityBilinearFormSet_2D, ElasticityBilinearFormSet_3D
   End Interface 
   
   Interface MEF90ElasticityEnergySet
      Module procedure ElasticityEnergySet_2D,ElasticityEnergySet_3D
   End Interface 

   Interface MEF90ElasticityOperatorAddTransientTermSet
      Module procedure ElasticityOperatorAddTransientTermSet_2D,ElasticityOperatorAddTransientTermSet_3D
   End Interface 
   
   Interface MEF90ElasticityOperatorSet
      Module procedure ElasticityOperatorSet_2D,ElasticityOperatorSet_3D
   End Interface 
   
   Interface MEF90ElasticityInelasticStrainRHSSetCst
      Module Procedure ElasticityInelasticStrainRHSSetCst_2D, ElasticityInelasticStrainRHSSetCst_3D
   End Interface 

   Interface MEF90ElasticityInelasticStrainRHSSetCell
      Module Procedure ElasticityInelasticStrainRHSSetCell_2D, ElasticityInelasticStrainRHSSetCell_3D, &
                       ElasticityInelasticStrainRHSSetCell2_2D, ElasticityInelasticStrainRHSSetCell2_3D
   End Interface 

   Interface MEF90ElasticityInelasticStrainRHSSetVertex
      Module Procedure ElasticityInelasticStrainRHSSetVertex2_2D, ElasticityInelasticStrainRHSSetVertex2_3D
   End Interface 

   Interface MEF90ElasticityForceRHSSetCst
      Module procedure ElasticityForceRHSSetCst_2D, ElasticityForceRHSSetCst_3D
   End Interface 

   Interface MEF90ElasticityForceRHSSetCell
      Module procedure ElasticityForceRHSSetCell_2D, ElasticityForceRHSSetCell_3D
   End Interface 

   Interface MEF90ElasticityForceRHSSetVertex
      Module procedure ElasticityForceRHSSetVertex_2D, ElasticityForceRHSSetVertex_3D
   End Interface 

   Interface MEF90ElasticityPressureForceRHSSetCst
      Module procedure ElasticityPressureForceRHSSetCst_2D, ElasticityPressureForceRHSSetCst_3D
   End Interface 

   Interface MEF90ElasticityPressureForceRHSSetCell
      Module procedure ElasticityPressureForceRHSSetCell_2D, ElasticityPressureForceRHSSetCell_3D
   End Interface 

   Interface MEF90ElasticityPressureForceRHSSetVertex
      Module procedure ElasticityPressureForceRHSSetVertex_2D, ElasticityPressureForceRHSSetVertex_3D
   End Interface 

   Interface MEF90ElasticityWorkSetCst
      Module procedure ElasticityWorkSetCst_2D, ElasticityWorkSetCst_3D
   End Interface 

   Interface MEF90ElasticityWorkSetCell
      Module procedure ElasticityWorkSetCell_2D, ElasticityWorkSetCell_3D
   End Interface 

   Interface MEF90ElasticityWorkSetVertex
      Module procedure ElasticityWorkSetVertex_2D, ElasticityWorkSetVertex_3D
   End Interface 
   
   Interface MEF90ElasticityCohesiveEnergySet
      Module Procedure ElasticityCohesiveEnergySet_2D, ElasticityCohesiveEnergySet_3D
   End Interface

   !Interface MEF90ElasticityPressureWorkSetCst
   !   Module procedure ElasticityPressureWorkSetCst_2D, ElasticityPressureWorkSetCst_3D
   !End Interface MEF90ElasticityPressureWorkSetCst
   
   Interface MEF90ElasticityPressureWorkSetCell
      Module procedure ElasticityPressureWorkSetCell_2D, ElasticityPressureWorkSetCell_3D
   End Interface 

   !Interface MEF90ElasticityPressureWorkSetVertex
   !   Module procedure ElasticityPressureWorkSetVertex_2D, ElasticityPressureWorkSetVertex_3D
   !End Interface MEF90ElasticityPressureWorkSetVertex

   Interface MEF90ElasticityStressSet
      Module procedure ElasticityStressSet_2D, ElasticityStressSet_3D
   End Interface 

   Interface MEF90InelasticStrainSet
      Module procedure InelasticStrainSet_2D, InelasticStrainSet_3D
   End Interface 

   Interface MEF90PlasticityEnergySet
      Module procedure PlasticityEnergySet_2D,PlasticityEnergySet_3D
   End Interface

End Module m_MEF90_ElasticityInterface
