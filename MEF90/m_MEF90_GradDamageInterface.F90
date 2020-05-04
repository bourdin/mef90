#include "../MEF90/mef90.inc"
Module m_MEF90_GradDamageInterface
#include "finclude/petscdef.h"
   Use m_MEF90_LinAlg
   Use m_MEF90_Parameters
   Use m_MEF90_Elements
   Use m_MEF90_Utils

   Use m_MEF90_GradDamageImplementation_2D , &
      MEF90GradDamageCellAverage_2D                          => MEF90GradDamageCellAverage
        
   Use m_MEF90_GradDamageImplementation_3D , &  
      MEF90GradDamageCellAverage_3D                          => MEF90GradDamageCellAverage

   IMPLICIT NONE
   Private
   
   Public :: MEF90GradDamageCellAverage
   Interface MEF90GradDamageCellAverage
      Module Procedure MEF90GradDamageCellAverage_2D, MEF90GradDamageCellAverage_3D
   End Interface 
End Module m_MEF90_GradDamageInterface