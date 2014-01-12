#include "../MEF90/mef90.inc"
Module m_MEF90_DefMechGradientDamageInterface
#include "finclude/petscdef.h"
   Use m_MEF90
   Use m_MEF90_DefMechCtx

   Use m_MEF90_DefMechGradientDamageImplementation_2D , &
      MEF90DefMechGradientDamageBilinearFormSet_2D => MEF90DefMechGradientDamageBilinearFormSet
   Use m_MEF90_DefMechGradientDamageImplementation_3D , &
      MEF90DefMechGradientDamageBilinearFormSet_3D => MEF90DefMechGradientDamageBilinearFormSet
      
   IMPLICIT NONE
   Private
   
   Public :: MEF90DefMechGradientDamageBilinearFormSet
   
   Interface MEF90DefMechGradientDamageBilinearFormSet
      Module Procedure MEF90DefMechGradientDamageBilinearFormSet_2D, MEF90DefMechGradientDamageBilinearFormSet_3D
   End Interface
End Module m_MEF90_DefMechGradientDamageInterface