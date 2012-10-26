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
         DiffusionRHSSetVertex_2D => DiffusionRHSSetVertex, &
         DiffusionRHSSetCell_2D => DiffusionRHSSetCell, &
         DiffusionRHSSetCst_2D => DiffusionRHSSetCst, &
         DiffusionWorkSetVertex_2D => DiffusionWorkSetVertex, &
         DiffusionWorkSetCell_2D => DiffusionWorkSetCell, &
         DiffusionWorkSetCst_2D => DiffusionWorkSetCst
   Use m_MEF_DiffusionImplementation_3D, &
         DiffusionBilinearFormSet_3D => DiffusionBilinearFormSet, &
         DiffusionEnergySet_3D => DiffusionEnergySet, &
         DiffusionOperatorAddTransientTermSet_3D => DiffusionOperatorAddTransientTermSet, &
         DiffusionOperatorSet_3D => DiffusionOperatorSet, &
         DiffusionRHSSetVertex_3D => DiffusionRHSSetVertex, &
         DiffusionRHSSetCell_3D => DiffusionRHSSetCell, &
         DiffusionRHSSetCst_3D => DiffusionRHSSetCst, &
         DiffusionWorkSetVertex_3D => DiffusionWorkSetVertex, &
         DiffusionWorkSetCell_3D => DiffusionWorkSetCell, &
         DiffusionWorkSetCst_3D => DiffusionWorkSetCst
   Use petsc

   IMPLICIT NONE

   Private   
   Public :: MEF90Diffusion_BilinearFormSet
   Public :: MEF90Diffusion_EnergySet
   Public :: MEF90Diffusion_OperatorAddTransientTermSet
   Public :: MEF90Diffusion_OperatorSet
   Public :: MEF90Diffusion_RHSSetVertex
   Public :: MEF90Diffusion_RHSSetCell
   Public :: MEF90Diffusion_RHSSetCst
   Public :: MEF90Diffusion_WorkSetVertex
   Public :: MEF90Diffusion_WorkSetCell
   Public :: MEF90Diffusion_WorkSetCst
   
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
   
   Interface MEF90Diffusion_RHSSetVertex
      Module Procedure DiffusionRHSSetVertex_2D,DiffusionRHSSetVertex_3D
   End Interface MEF90Diffusion_RHSSetVertex

   Interface MEF90Diffusion_RHSSetCell
      Module Procedure DiffusionRHSSetCell_2D,DiffusionRHSSetCell_3D
   End Interface MEF90Diffusion_RHSSetCell

   Interface MEF90Diffusion_RHSSetCst
      Module Procedure DiffusionRHSSetCst_2D,DiffusionRHSSetCst_3D
   End Interface MEF90Diffusion_RHSSetCst

   Interface MEF90Diffusion_WorkSetVertex
      Module Procedure DiffusionWorkSetVertex_2D,DiffusionWorkSetVertex_3D
   End Interface MEF90Diffusion_WorkSetVertex

   Interface MEF90Diffusion_WorkSetCell
      Module Procedure DiffusionWorkSetCell_2D,DiffusionWorkSetCell_3D
   End Interface MEF90Diffusion_WorkSetCell

   Interface MEF90Diffusion_WorkSetCst
      Module Procedure DiffusionWorkSetCst_2D,DiffusionWorkSetCst_3D
   End Interface MEF90Diffusion_WorkSetCst
End Module m_MEF_DiffusionInterface
