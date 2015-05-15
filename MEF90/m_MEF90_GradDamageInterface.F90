#include "../MEF90/mef90.inc"
Module m_MEF90_GradDamageInterface
#include "finclude/petscdef.h"
   Use m_MEF90_LinAlg
   Use m_MEF90_Parameters
   Use m_MEF90_Elements
   Use m_MEF90_Utils

   Use m_MEF90_GradDamageImplementation_2D , &
      MEF90GradDamageDispBilinearFormSet_2D                  => MEF90GradDamageDispBilinearFormSet, &
      MEF90GradDamageDispOperatorSet_2D                      => MEF90GradDamageDispOperatorSet, &
      MEF90GradDamageDispInelasticStrainRHSSetVertex2_2D     => MEF90GradDamageDispInelasticStrainRHSSetVertex2, &
      MEF90GradDamageDispInelasticStrainRHSSetCell_2D        => MEF90GradDamageDispInelasticStrainRHSSetCell, &
      MEF90GradDamageElasticEnergySet_2D                     => MEF90GradDamageElasticEnergySet, &
      MEF90GradDamageDamageBilinearFormSetAT1Elastic_2D      => MEF90GradDamageDamageBilinearFormSetAT1Elastic, &
      MEF90GradDamageDamageOperatorSetAT1Elastic_2D          => MEF90GradDamageDamageOperatorSetAT1Elastic, &
      MEF90GradDamageDamageRHSSetAT1Elastic_2D               => MEF90GradDamageDamageRHSSetAT1Elastic, &
      MEF90GradDamageDamageBilinearFormSetAT1_2D             => MEF90GradDamageDamageBilinearFormSetAT1, &
      MEF90GradDamageDamageOperatorSetAT1_2D                 => MEF90GradDamageDamageOperatorSetAT1, &
      MEF90GradDamageDamageRHSSetAT1_2D                      => MEF90GradDamageDamageRHSSetAT1, &
      MEF90GradDamageSurfaceEnergySetAT1_2D                  => MEF90GradDamageSurfaceEnergySetAT1, &
      MEF90GradDamageDamageBilinearFormSetAT2Elastic_2D      => MEF90GradDamageDamageBilinearFormSetAT2Elastic, &
      MEF90GradDamageDamageOperatorSetAT2Elastic_2D          => MEF90GradDamageDamageOperatorSetAT2Elastic, &
      MEF90GradDamageDamageBilinearFormSetAT2_2D             => MEF90GradDamageDamageBilinearFormSetAT2, &
      MEF90GradDamageDamageOperatorSetAT2_2D                 => MEF90GradDamageDamageOperatorSetAT2, &
      MEF90GradDamageDamageRHSSetAT2_2D                      => MEF90GradDamageDamageRHSSetAT2, &
      MEF90GradDamageSurfaceEnergySetAT2_2D                  => MEF90GradDamageSurfaceEnergySetAT2, &

!!erwan-->!!
!      MEF90GradDamageDamageBilinearFormSetATkElastic_2D      => MEF90GradDamageDamageBilinearFormSetATkElastic, &
!      MEF90GradDamageDamageOperatorSetATkElastic_2D          => MEF90GradDamageDamageOperatorSetATkElastic, &
!      MEF90GradDamageDamageRHSSetATkElastic_2D               => MEF90GradDamageDamageRHSSetATkElastic, &
      MEF90GradDamageDamageBilinearFormSetATk_2D             => MEF90GradDamageDamageBilinearFormSetATk, &
      MEF90GradDamageDamageOperatorSetATk_2D                 => MEF90GradDamageDamageOperatorSetATk, &
      MEF90GradDamageDamageRHSSetATk_2D                      => MEF90GradDamageDamageRHSSetATk, &
      MEF90GradDamageSurfaceEnergySetATk_2D                  => MEF90GradDamageSurfaceEnergySetATk
!!<--erwan!!

        
   Use m_MEF90_GradDamageImplementation_3D , &  
      MEF90GradDamageDispBilinearFormSet_3D                  => MEF90GradDamageDispBilinearFormSet, &
      MEF90GradDamageDispOperatorSet_3D                      => MEF90GradDamageDispOperatorSet, &
      MEF90GradDamageDispInelasticStrainRHSSetVertex2_3D     => MEF90GradDamageDispInelasticStrainRHSSetVertex2, &
      MEF90GradDamageDispInelasticStrainRHSSetCell_3D        => MEF90GradDamageDispInelasticStrainRHSSetCell, &
      MEF90GradDamageElasticEnergySet_3D                     => MEF90GradDamageElasticEnergySet, &
      MEF90GradDamageDamageBilinearFormSetAT1Elastic_3D      => MEF90GradDamageDamageBilinearFormSetAT1Elastic, &
      MEF90GradDamageDamageOperatorSetAT1Elastic_3D          => MEF90GradDamageDamageOperatorSetAT1Elastic, &
      MEF90GradDamageDamageRHSSetAT1Elastic_3D               => MEF90GradDamageDamageRHSSetAT1Elastic, &
      MEF90GradDamageDamageBilinearFormSetAT1_3D             => MEF90GradDamageDamageBilinearFormSetAT1, &
      MEF90GradDamageDamageOperatorSetAT1_3D                 => MEF90GradDamageDamageOperatorSetAT1, &
      MEF90GradDamageDamageRHSSetAT1_3D                      => MEF90GradDamageDamageRHSSetAT1, &
      MEF90GradDamageSurfaceEnergySetAT1_3D                  => MEF90GradDamageSurfaceEnergySetAT1, &
      MEF90GradDamageDamageBilinearFormSetAT2Elastic_3D      => MEF90GradDamageDamageBilinearFormSetAT2Elastic, &
      MEF90GradDamageDamageOperatorSetAT2Elastic_3D          => MEF90GradDamageDamageOperatorSetAT2Elastic, &
      MEF90GradDamageDamageBilinearFormSetAT2_3D             => MEF90GradDamageDamageBilinearFormSetAT2, &
      MEF90GradDamageDamageOperatorSetAT2_3D                 => MEF90GradDamageDamageOperatorSetAT2, &
      MEF90GradDamageDamageRHSSetAT2_3D                      => MEF90GradDamageDamageRHSSetAT2, &
      MEF90GradDamageSurfaceEnergySetAT2_3D                  => MEF90GradDamageSurfaceEnergySetAT2, &
!!erwan-->!!
!      MEF90GradDamageDamageBilinearFormSetATkElastic_3D      => MEF90GradDamageDamageBilinearFormSetATkElastic, &
!      MEF90GradDamageDamageOperatorSetATkElastic_3D          => MEF90GradDamageDamageOperatorSetATkElastic, &
!      MEF90GradDamageDamageRHSSetATkElastic_3D               => MEF90GradDamageDamageRHSSetATkElastic, &
      MEF90GradDamageDamageBilinearFormSetATk_3D             => MEF90GradDamageDamageBilinearFormSetATk, &
      MEF90GradDamageDamageOperatorSetATk_3D                 => MEF90GradDamageDamageOperatorSetATk, &
      MEF90GradDamageDamageRHSSetATk_3D                      => MEF90GradDamageDamageRHSSetATk, &
      MEF90GradDamageSurfaceEnergySetATk_3D                  => MEF90GradDamageSurfaceEnergySetATk
!!<--erwan!!


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
   Public :: MEF90GradDamageDamageOperatorSetAT1
   Public :: MEF90GradDamageDamageRHSSetAT1
   Public :: MEF90GradDamageSurfaceEnergySetAT1
   Public :: MEF90GradDamageDamageBilinearFormSetAT2Elastic
   Public :: MEF90GradDamageDamageOperatorSetAT2Elastic
   Public :: MEF90GradDamageDamageBilinearFormSetAT2
   Public :: MEF90GradDamageDamageOperatorSetAT2
   Public :: MEF90GradDamageDamageRHSSetAT2
   Public :: MEF90GradDamageSurfaceEnergySetAT2

!!erwan-->!!
!   Public :: MEF90GradDamageDamageBilinearFormSetATkElastic
!   Public :: MEF90GradDamageDamageOperatorSetATkElastic
!   Public :: MEF90GradDamageDamageRHSSetATkElastic
   Public :: MEF90GradDamageDamageBilinearFormSetATk
   Public :: MEF90GradDamageDamageOperatorSetATk
   Public :: MEF90GradDamageDamageRHSSetATk
   Public :: MEF90GradDamageSurfaceEnergySetATk
!!<--erwan!!
   
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

   Interface MEF90GradDamageDamageOperatorSetAT1
      Module Procedure MEF90GradDamageDamageOperatorSetAT1_2D, MEF90GradDamageDamageOperatorSetAT1_3D
   End Interface 

   Interface MEF90GradDamageDamageRHSSetAT1
      Module Procedure MEF90GradDamageDamageRHSSetAT1_2D, MEF90GradDamageDamageRHSSetAT1_3D
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

   Interface MEF90GradDamageDamageOperatorSetAT2
      Module Procedure MEF90GradDamageDamageOperatorSetAT2_2D, MEF90GradDamageDamageOperatorSetAT2_3D
   End Interface 

   Interface MEF90GradDamageDamageRHSSetAT2
      Module Procedure MEF90GradDamageDamageRHSSetAT2_2D, MEF90GradDamageDamageRHSSetAT2_3D
   End Interface 

   Interface MEF90GradDamageSurfaceEnergySetAT2
      Module Procedure MEF90GradDamageSurfaceEnergySetAT2_2D, MEF90GradDamageSurfaceEnergySetAT2_3D
   End Interface 

!!erwan-->!!

   Interface MEF90GradDamageDamageBilinearFormSetATk
      Module Procedure MEF90GradDamageDamageBilinearFormSetATk_2D, MEF90GradDamageDamageBilinearFormSetATk_3D
   End Interface 

   Interface MEF90GradDamageDamageOperatorSetATk
      Module Procedure MEF90GradDamageDamageOperatorSetATk_2D, MEF90GradDamageDamageOperatorSetATk_3D
   End Interface 

   Interface MEF90GradDamageDamageRHSSetATk
      Module Procedure MEF90GradDamageDamageRHSSetATk_2D, MEF90GradDamageDamageRHSSetATk_3D
   End Interface 

   Interface MEF90GradDamageSurfaceEnergySetATk
      Module Procedure MEF90GradDamageSurfaceEnergySetATk_2D, MEF90GradDamageSurfaceEnergySetATk_3D
   End Interface 
!!<--erwan!!


End Module m_MEF90_GradDamageInterface