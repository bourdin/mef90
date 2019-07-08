#include "../MEF90/mef90.inc"
Module m_MEF90_GradDamageInterface
#include "finclude/petscdef.h"
   Use m_MEF90_LinAlg
   Use m_MEF90_Parameters
   Use m_MEF90_Elements
   Use m_MEF90_Utils

   Use m_MEF90_GradDamageImplementation_2D , &
      MEF90GradDamageElasticEnergySet_2D                     => MEF90GradDamageElasticEnergySet, &
      MEF90GradDamageSurfaceEnergySetAT1_2D                  => MEF90GradDamageSurfaceEnergySetAT1, &
      MEF90GradDamageSurfaceEnergySetAT2_2D                  => MEF90GradDamageSurfaceEnergySetAT2, &
      MEF90GradDamageSurfaceEnergySetLinSoft_2D              => MEF90GradDamageSurfaceEnergySetLinSoft, &
      MEF90GradDamageCellAverage_2D                          => MEF90GradDamageCellAverage
        
   Use m_MEF90_GradDamageImplementation_3D , &  
      MEF90GradDamageElasticEnergySet_3D                     => MEF90GradDamageElasticEnergySet, &
      MEF90GradDamageSurfaceEnergySetAT1_3D                  => MEF90GradDamageSurfaceEnergySetAT1, &
      MEF90GradDamageSurfaceEnergySetAT2_3D                  => MEF90GradDamageSurfaceEnergySetAT2, &
      MEF90GradDamageSurfaceEnergySetLinSoft_3D              => MEF90GradDamageSurfaceEnergySetLinSoft, &
      MEF90GradDamageCellAverage_3D                          => MEF90GradDamageCellAverage

   IMPLICIT NONE
   Private
   
   Public :: MEF90GradDamageElasticEnergySet
   Public :: MEF90GradDamageSurfaceEnergySetAT1
   Public :: MEF90GradDamageSurfaceEnergySetAT2
   Public :: MEF90GradDamageSurfaceEnergySetLinSoft
   Public :: MEF90GradDamageCellAverage

   Interface MEF90GradDamageElasticEnergySet
      Module Procedure MEF90GradDamageElasticEnergySet_2D, MEF90GradDamageElasticEnergySet_3D
   End Interface
   
   Interface MEF90GradDamageSurfaceEnergySetAT1
      Module Procedure MEF90GradDamageSurfaceEnergySetAT1_2D, MEF90GradDamageSurfaceEnergySetAT1_3D
   End Interface 
   
   Interface MEF90GradDamageSurfaceEnergySetAT2
      Module Procedure MEF90GradDamageSurfaceEnergySetAT2_2D, MEF90GradDamageSurfaceEnergySetAT2_3D
   End Interface 

   Interface MEF90GradDamageSurfaceEnergySetLinSoft
      Module Procedure MEF90GradDamageSurfaceEnergySetLinSoft_2D, MEF90GradDamageSurfaceEnergySetLinSoft_3D
   End Interface 

   Interface MEF90GradDamageCellAverage
      Module Procedure MEF90GradDamageCellAverage_2D, MEF90GradDamageCellAverage_3D
   End Interface 
End Module m_MEF90_GradDamageInterface