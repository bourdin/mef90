#include "../MEF90/mef90.inc"
Module m_MEF90_GradDamageInterface
#include "finclude/petscdef.h"
   Use m_MEF90_LinAlg
   Use m_MEF90_Parameters
   Use m_MEF90_Elements
   Use m_MEF90_Utils

   Use m_MEF90_GradDamageImplementation_2D , &
      MEF90GradDamageDispBilinearFormSet_2D              => MEF90GradDamageDispBilinearFormSet, &
      MEF90GradDamageDispOperatorSet_2D                  => MEF90GradDamageDispOperatorSet, &
      MEF90GradDamageDispInelasticStrainRHSSetVertex2_2D => MEF90GradDamageDispInelasticStrainRHSSetVertex2, &
      MEF90GradDamageDispInelasticStrainRHSSetCell_2D    => MEF90GradDamageDispInelasticStrainRHSSetCell, &
      MEF90GradDamageDamageBilinearFormSetAT1_2D         => MEF90GradDamageDamageBilinearFormSetAT1, &
      MEF90GradDamageDamageBilinearFormSetAT2_2D         => MEF90GradDamageDamageBilinearFormSetAT2
      
   Use m_MEF90_GradDamageImplementation_3D , &
      MEF90GradDamageDispBilinearFormSet_3D              => MEF90GradDamageDispBilinearFormSet, &
      MEF90GradDamageDispOperatorSet_3D                  => MEF90GradDamageDispOperatorSet, &
      MEF90GradDamageDispInelasticStrainRHSSetVertex2_3D => MEF90GradDamageDispInelasticStrainRHSSetVertex2, &
      MEF90GradDamageDispInelasticStrainRHSSetCell_3D    => MEF90GradDamageDispInelasticStrainRHSSetCell, &
      MEF90GradDamageDamageBilinearFormSetAT1_3D         => MEF90GradDamageDamageBilinearFormSetAT1, &
      MEF90GradDamageDamageBilinearFormSetAT2_3D         => MEF90GradDamageDamageBilinearFormSetAT2
      
   IMPLICIT NONE
   Private
   
   Public :: MEF90GradDamageDispBilinearFormSet
   Public :: MEF90GradDamageDispOperatorSet
   Public :: MEF90GradDamageDispInelasticStrainRHSSetVertex
   Public :: MEF90GradDamageDispInelasticStrainRHSSetCell
   Public :: MEF90GradDamageDamageBilinearFormSetAT1
   Public :: MEF90GradDamageDamageBilinearFormSetAT2
   
   Interface MEF90GradDamageDispBilinearFormSet
      Module Procedure MEF90GradDamageDispBilinearFormSet_2D, MEF90GradDamageDispBilinearFormSet_3D
   End Interface
   
   Interface MEF90GradDamageDispOperatorSet
      Module Procedure MEF90GradDamageDispOperatorSet_2D, MEF90GradDamageDispOperatorSet_3D
   End Interface
   
   Interface MEF90GradDamageDispInelasticStrainRHSSetVertex
      Module Procedure MEF90GradDamageDispInelasticStrainRHSSetVertex2_2D, MEF90GradDamageDispInelasticStrainRHSSetVertex2_3D
   End Interface MEF90GradDamageDispInelasticStrainRHSSetVertex

   Interface MEF90GradDamageDispInelasticStrainRHSSetCell
      Module Procedure MEF90GradDamageDispInelasticStrainRHSSetCell_2D, MEF90GradDamageDispInelasticStrainRHSSetCell_3D
   End Interface MEF90GradDamageDispInelasticStrainRHSSetCell

   Interface MEF90GradDamageDamageBilinearFormSetAT1
      Module Procedure MEF90GradDamageDamageBilinearFormSetAT1_2D, MEF90GradDamageDamageBilinearFormSetAT1_3D
   End Interface MEF90GradDamageDamageBilinearFormSetAT1

   Interface MEF90GradDamageDamageBilinearFormSetAT2
      Module Procedure MEF90GradDamageDamageBilinearFormSetAT2_2D, MEF90GradDamageDamageBilinearFormSetAT2_3D
   End Interface MEF90GradDamageDamageBilinearFormSetAT2
End Module m_MEF90_GradDamageInterface