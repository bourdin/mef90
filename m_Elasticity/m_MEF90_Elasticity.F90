#include "../MEF90/mef90.inc"
Module m_MEF90_Elasticity
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_MEF90_ElasticityCtx
   Use m_MEF90_ElasticityAssembly2D, &
      MEF90ElasticityOperator2D     => MEF90ElasticityOperator, &
      MEF90ElasticityBilinearForm2D => MEF90ElasticityBilinearForm      
   Use m_MEF90_ElasticityAssembly3D, &
      MEF90ElasticityOperator3D     => MEF90ElasticityOperator, &
      MEF90ElasticityBilinearForm3D => MEF90ElasticityBilinearForm

   Implicit none
   Private
   Public MEF90ElasticityOperator
   Public MEF90ElasticityBilinearForm
   
   Interface MEF90ElasticityOperator
      Module Procedure MEF90ElasticityOperator2D,MEF90ElasticityOperator3D
   End Interface

   Interface MEF90ElasticityBilinearForm
      Module Procedure MEF90ElasticityBilinearForm2D,MEF90ElasticityBilinearForm3D
   End Interface
   !Public MEF90ElasticitySetTransients

   
!Contains
End Module m_MEF90_Elasticity
