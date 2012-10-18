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
         DiffusionRHSSet_2D => DiffusionRHSSet, &
         DiffusionRHSCellCstSet_2D => DiffusionRHSCellCstSet, &
         DiffusionWorkSet_2D => DiffusionWorkSet, &
         DiffusionWorkCellCstSet_2D => DiffusionWorkCellCstSet
   Use m_MEF_DiffusionImplementation_3D, &
         DiffusionBilinearFormSet_3D => DiffusionBilinearFormSet, &
         DiffusionEnergySet_3D => DiffusionEnergySet, &
         DiffusionOperatorAddTransientTermSet_3D => DiffusionOperatorAddTransientTermSet, &
         DiffusionOperatorSet_3D => DiffusionOperatorSet, &
         DiffusionRHSSet_3D => DiffusionRHSSet, &
         DiffusionRHSCellCstSet_3D => DiffusionRHSCellCstSet, &
         DiffusionWorkSet_3D => DiffusionWorkSet, &
         DiffusionWorkCellCstSet_3D => DiffusionWorkCellCstSet
   Use petsc

   IMPLICIT NONE

   Private   
   Public :: MEF90Diffusion_BilinearFormSet
   Public :: MEF90Diffusion_EnergySet
   Public :: MEF90Diffusion_OperatorAddTransientTermSet
   Public :: MEF90Diffusion_OperatorSet
   Public :: MEF90Diffusion_RHSSet
   Public :: MEF90Diffusion_WorkSet
   
   Interface MEF90Diffusion_BilinearFormSet
      Module Procedure DiffusionBilinearFormSet_2D, DiffusionBilinearFormSet_3D
   End Interface MEF90Diffusion_BilinearFormSet
   
   Interface MEF90Diffusion_EnergySet
      Module procedure DiffusionEnergySet_2D,DiffusionEnergySet_3D
   End Interface MEF90Diffusion_EnergySet

   Interface MEF90Diffusion_OperatorAddTransientTermSet
      Module procedure DiffusionOperatorAddTransientTermSet_2D,DiffusionOperatorAddTransientTermSet_3D
   End Interface MEF90Diffusion_OperatorAddTransientTermSet
   
   Interface MEF90Diffusion_OperatorSet
      Module procedure DiffusionOperatorSet_2D,DiffusionOperatorSet_3D
   End Interface MEF90Diffusion_OperatorSet

   Interface MEF90Diffusion_RHSSet
      Module procedure DiffusionRHSSet_2D,DiffusionRHSCellCstSet_2D,DiffusionRHSSet_3D,DiffusionRHSCellCstSet_3D
   End Interface MEF90Diffusion_RHSSet

   Interface MEF90Diffusion_WorkSet
      Module procedure DiffusionWorkSet_2D,DiffusionWorkCellCstSet_2D,DiffusionWorkSet_3D,DiffusionWorkCellCstSet_3D
   End Interface MEF90Diffusion_WorkSet

End Module m_MEF_DiffusionInterface
