#include "../MEF90/mef90.inc"
Module m_MEF90_GradDamageInterface
#include "finclude/petscdef.h"
   Use m_MEF90_LinAlg
   Use m_MEF90_Parameters
   Use m_MEF90_Elements
   Use m_MEF90_Utils

   Use m_MEF90_GradDamageImplementation_2D , &
      MEF90GradDamageDispBilinearFormSet_2D                => MEF90GradDamageDispBilinearFormSet, &
      MEF90GradDamageDispOperatorSet_2D                    => MEF90GradDamageDispOperatorSet, &
      MEF90GradDamageDispInelasticStrainRHSSetVertex2_2D   => MEF90GradDamageDispInelasticStrainRHSSetVertex2, &
      MEF90GradDamageDispInelasticStrainRHSSetCell_2D      => MEF90GradDamageDispInelasticStrainRHSSetCell, &
      MEF90GradDamageElasticEnergySet_2D                   => MEF90GradDamageElasticEnergySet, &
      MEF90GradDamageDamageBilinearFormSetAT1Elastic_2D    => MEF90GradDamageDamageBilinearFormSetAT1Elastic, &
      MEF90GradDamageDamageOperatorSetAT1Elastic_2D        => MEF90GradDamageDamageOperatorSetAT1Elastic, &
      MEF90GradDamageDamageRHSSetAT1Elastic_2D             => MEF90GradDamageDamageRHSSetAT1Elastic, &
      MEF90GradDamageDamageBilinearFormSetAT1_2D           => MEF90GradDamageDamageBilinearFormSetAT1, &
      MEF90GradDamageDamageBilinearFormSetAT1Unilateral_2D => MEF90GradDamageDamageBilinearFormSetAT1Unilateral, &
      MEF90GradDamageDamageOperatorSetAT1_2D               => MEF90GradDamageDamageOperatorSetAT1, &
      MEF90GradDamageDamageOperatorSetAT1Unilateral_2D     => MEF90GradDamageDamageOperatorSetAT1Unilateral, &
      MEF90GradDamageDamageRHSSetAT1_2D                    => MEF90GradDamageDamageRHSSetAT1, &
      MEF90GradDamageDamageRHSSetAT1Unilateral_2D          => MEF90GradDamageDamageRHSSetAT1Unilateral, &
      MEF90GradDamageSurfaceEnergySetAT1_2D                => MEF90GradDamageSurfaceEnergySetAT1, &
      MEF90GradDamageDamageBilinearFormSetAT2Elastic_2D    => MEF90GradDamageDamageBilinearFormSetAT2Elastic, &
      MEF90GradDamageDamageOperatorSetAT2Elastic_2D        => MEF90GradDamageDamageOperatorSetAT2Elastic, &
      MEF90GradDamageDamageBilinearFormSetAT2_2D           => MEF90GradDamageDamageBilinearFormSetAT2, &
      MEF90GradDamageDamageBilinearFormSetAT2Unilateral_2D => MEF90GradDamageDamageBilinearFormSetAT2Unilateral, &
      MEF90GradDamageDamageOperatorSetAT2_2D               => MEF90GradDamageDamageOperatorSetAT2, &
      MEF90GradDamageDamageOperatorSetAT2Unilateral_2D     => MEF90GradDamageDamageOperatorSetAT2Unilateral, &
      MEF90GradDamageDamageRHSSetAT2_2D                    => MEF90GradDamageDamageRHSSetAT2, &
      MEF90GradDamageDamageRHSSetAT2Unilateral_2D          => MEF90GradDamageDamageRHSSetAT2Unilateral, &
      MEF90GradDamageSurfaceEnergySetAT2_2D                => MEF90GradDamageSurfaceEnergySetAT2
      
   Use m_MEF90_GradDamageImplementation_3D , &
      MEF90GradDamageDispBilinearFormSet_3D                => MEF90GradDamageDispBilinearFormSet, &
      MEF90GradDamageDispOperatorSet_3D                    => MEF90GradDamageDispOperatorSet, &
      MEF90GradDamageDispInelasticStrainRHSSetVertex2_3D   => MEF90GradDamageDispInelasticStrainRHSSetVertex2, &
      MEF90GradDamageDispInelasticStrainRHSSetCell_3D      => MEF90GradDamageDispInelasticStrainRHSSetCell, &
      MEF90GradDamageElasticEnergySet_3D                   => MEF90GradDamageElasticEnergySet, &
      MEF90GradDamageDamageBilinearFormSetAT1Elastic_3D    => MEF90GradDamageDamageBilinearFormSetAT1Elastic, &
      MEF90GradDamageDamageOperatorSetAT1Elastic_3D        => MEF90GradDamageDamageOperatorSetAT1Elastic, &
      MEF90GradDamageDamageRHSSetAT1Elastic_3D             => MEF90GradDamageDamageRHSSetAT1Elastic, &
      MEF90GradDamageDamageBilinearFormSetAT1_3D           => MEF90GradDamageDamageBilinearFormSetAT1, &
      MEF90GradDamageDamageBilinearFormSetAT1Unilateral_3D => MEF90GradDamageDamageBilinearFormSetAT1Unilateral, &
      MEF90GradDamageDamageOperatorSetAT1_3D               => MEF90GradDamageDamageOperatorSetAT1, &
      MEF90GradDamageDamageOperatorSetAT1Unilateral_3D     => MEF90GradDamageDamageOperatorSetAT1Unilateral, &
      MEF90GradDamageDamageRHSSetAT1_3D                    => MEF90GradDamageDamageRHSSetAT1, &
      MEF90GradDamageDamageRHSSetAT1Unilateral_3D          => MEF90GradDamageDamageRHSSetAT1Unilateral, &
      MEF90GradDamageSurfaceEnergySetAT1_3D                => MEF90GradDamageSurfaceEnergySetAT1, &
      MEF90GradDamageDamageBilinearFormSetAT2Elastic_3D    => MEF90GradDamageDamageBilinearFormSetAT2Elastic, &
      MEF90GradDamageDamageOperatorSetAT2Elastic_3D        => MEF90GradDamageDamageOperatorSetAT2Elastic, &
      MEF90GradDamageDamageBilinearFormSetAT2_3D           => MEF90GradDamageDamageBilinearFormSetAT2, &
      MEF90GradDamageDamageBilinearFormSetAT2Unilateral_3D => MEF90GradDamageDamageBilinearFormSetAT2Unilateral, &
      MEF90GradDamageDamageOperatorSetAT2_3D               => MEF90GradDamageDamageOperatorSetAT2, &
      MEF90GradDamageDamageOperatorSetAT2Unilateral_3D     => MEF90GradDamageDamageOperatorSetAT2Unilateral, &
      MEF90GradDamageDamageRHSSetAT2_3D                    => MEF90GradDamageDamageRHSSetAT2, &
      MEF90GradDamageDamageRHSSetAT2Unilateral_3D          => MEF90GradDamageDamageRHSSetAT2Unilateral, &
      MEF90GradDamageSurfaceEnergySetAT2_3D                => MEF90GradDamageSurfaceEnergySetAT2
      
   IMPLICIT NONE
   Private
   
   Public :: MEF90GradDamageDispBilinearFormSet
   Public :: MEF90GradDamageDispOperatorSet
   Public :: MEF90GradDamageDispInelasticStrainRHSSetVertex
   Public :: MEF90GradDamageDispInelasticStrainRHSSetCell
   Public :: MEF90GradDamageElasticEnergySet
   Public :: MEF90GradDamageDamageBilinearFormSetAT1Elastic
   Public :: MEF90GradDamageDamageOperatorSetAT1Elastic
   Public :: MEF90GradDamageDamageRHSSetAT1Elastic
   Public :: MEF90GradDamageDamageBilinearFormSetAT1
   Public :: MEF90GradDamageDamageBilinearFormSetAT1Unilateral
   Public :: MEF90GradDamageDamageOperatorSetAT1
   Public :: MEF90GradDamageDamageOperatorSetAT1Unilateral
   Public :: MEF90GradDamageDamageRHSSetAT1
   Public :: MEF90GradDamageDamageRHSSetAT1Unilateral
   Public :: MEF90GradDamageSurfaceEnergySetAT1
   Public :: MEF90GradDamageDamageBilinearFormSetAT2Elastic
   Public :: MEF90GradDamageDamageOperatorSetAT2Elastic
   Public :: MEF90GradDamageDamageBilinearFormSetAT2
   Public :: MEF90GradDamageDamageBilinearFormSetAT2Unilateral
   Public :: MEF90GradDamageDamageOperatorSetAT2
   Public :: MEF90GradDamageDamageOperatorSetAT2Unilateral
   Public :: MEF90GradDamageDamageRHSSetAT2
   Public :: MEF90GradDamageDamageRHSSetAT2Unilateral
   Public :: MEF90GradDamageSurfaceEnergySetAT2
   
   Interface MEF90GradDamageDispBilinearFormSet
      Module Procedure MEF90GradDamageDispBilinearFormSet_2D, MEF90GradDamageDispBilinearFormSet_3D
   End Interface
   
   Interface MEF90GradDamageDispOperatorSet
      Module Procedure MEF90GradDamageDispOperatorSet_2D, MEF90GradDamageDispOperatorSet_3D
   End Interface
   
   Interface MEF90GradDamageDispInelasticStrainRHSSetVertex
      Module Procedure MEF90GradDamageDispInelasticStrainRHSSetVertex2_2D, MEF90GradDamageDispInelasticStrainRHSSetVertex2_3D
   End Interface 

   Interface MEF90GradDamageDispInelasticStrainRHSSetCell
      Module Procedure MEF90GradDamageDispInelasticStrainRHSSetCell_2D, MEF90GradDamageDispInelasticStrainRHSSetCell_3D
   End Interface 

   Interface MEF90GradDamageElasticEnergySet
      Module Procedure MEF90GradDamageElasticEnergySet_2D, MEF90GradDamageElasticEnergySet_3D
   End Interface
   
   Interface MEF90GradDamageDamageBilinearFormSetAT1Elastic
      Module Procedure MEF90GradDamageDamageBilinearFormSetAT1Elastic_2D, MEF90GradDamageDamageBilinearFormSetAT1Elastic_3D
   End Interface 

   Interface MEF90GradDamageDamageOperatorSetAT1Elastic
      Module Procedure MEF90GradDamageDamageOperatorSetAT1Elastic_2D, MEF90GradDamageDamageOperatorSetAT1Elastic_3D
   End Interface 

   Interface MEF90GradDamageDamageRHSSetAT1Elastic
      Module Procedure MEF90GradDamageDamageRHSSetAT1Elastic_2D, MEF90GradDamageDamageRHSSetAT1Elastic_3D
   End Interface 

   Interface MEF90GradDamageDamageBilinearFormSetAT1
      Module Procedure MEF90GradDamageDamageBilinearFormSetAT1_2D, MEF90GradDamageDamageBilinearFormSetAT1_3D
   End Interface 

   Interface MEF90GradDamageDamageBilinearFormSetAT1Unilateral
      Module Procedure MEF90GradDamageDamageBilinearFormSetAT1Unilateral_2D, MEF90GradDamageDamageBilinearFormSetAT1Unilateral_3D
   End Interface 

   Interface MEF90GradDamageDamageOperatorSetAT1
      Module Procedure MEF90GradDamageDamageOperatorSetAT1_2D, MEF90GradDamageDamageOperatorSetAT1_3D
   End Interface 

   Interface MEF90GradDamageDamageOperatorSetAT1Unilateral
      Module Procedure MEF90GradDamageDamageOperatorSetAT1Unilateral_2D, MEF90GradDamageDamageOperatorSetAT1Unilateral_3D
   End Interface 

   Interface MEF90GradDamageDamageRHSSetAT1
      Module Procedure MEF90GradDamageDamageRHSSetAT1_2D, MEF90GradDamageDamageRHSSetAT1_3D
   End Interface 

   Interface MEF90GradDamageDamageRHSSetAT1Unilateral
      Module Procedure MEF90GradDamageDamageRHSSetAT1Unilateral_2D, MEF90GradDamageDamageRHSSetAT1Unilateral_3D
   End Interface 

   Interface MEF90GradDamageSurfaceEnergySetAT1
      Module Procedure MEF90GradDamageSurfaceEnergySetAT1_2D, MEF90GradDamageSurfaceEnergySetAT1_3D
   End Interface 
   
   Interface MEF90GradDamageDamageBilinearFormSetAT2Elastic
      Module Procedure MEF90GradDamageDamageBilinearFormSetAT2Elastic_2D, MEF90GradDamageDamageBilinearFormSetAT2Elastic_3D
   End Interface 

   Interface MEF90GradDamageDamageOperatorSetAT2Elastic
      Module Procedure MEF90GradDamageDamageOperatorSetAT2Elastic_2D, MEF90GradDamageDamageOperatorSetAT2Elastic_3D
   End Interface 

   Interface MEF90GradDamageDamageBilinearFormSetAT2
      Module Procedure MEF90GradDamageDamageBilinearFormSetAT2_2D, MEF90GradDamageDamageBilinearFormSetAT2_3D
   End Interface 

   Interface MEF90GradDamageDamageBilinearFormSetAT2Unilateral
      Module Procedure MEF90GradDamageDamageBilinearFormSetAT2Unilateral_2D, MEF90GradDamageDamageBilinearFormSetAT2Unilateral_3D
   End Interface 

   Interface MEF90GradDamageDamageOperatorSetAT2
      Module Procedure MEF90GradDamageDamageOperatorSetAT2_2D, MEF90GradDamageDamageOperatorSetAT2_3D
   End Interface 

   Interface MEF90GradDamageDamageOperatorSetAT2Unilateral
      Module Procedure MEF90GradDamageDamageOperatorSetAT2Unilateral_2D, MEF90GradDamageDamageOperatorSetAT2Unilateral_3D
   End Interface 

   Interface MEF90GradDamageDamageRHSSetAT2
      Module Procedure MEF90GradDamageDamageRHSSetAT2_2D, MEF90GradDamageDamageRHSSetAT2_3D
   End Interface 

   Interface MEF90GradDamageDamageRHSSetAT2Unilateral
      Module Procedure MEF90GradDamageDamageRHSSetAT2Unilateral_2D, MEF90GradDamageDamageRHSSetAT2Unilateral_3D
   End Interface 

   Interface MEF90GradDamageSurfaceEnergySetAT2
      Module Procedure MEF90GradDamageSurfaceEnergySetAT2_2D, MEF90GradDamageSurfaceEnergySetAT2_3D
   End Interface 
End Module m_MEF90_GradDamageInterface