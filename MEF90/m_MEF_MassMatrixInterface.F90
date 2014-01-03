Module m_MEF_MassMatrixInterface
#include "finclude/petscdef.h"
   Use m_MEF_LinAlg
   Use m_MEF_Parameters
   Use m_MEF_Elements
   Use m_MEF_Utils
   Use m_MEF_MassMatrixImplementation_MEF90Element2D_Scal,  MassMatrixAssembleSet_2DScal  => MEF90_MassMatrixAssembleSet
   Use m_MEF_MassMatrixImplementation_MEF90Element2D_Vect,  MassMatrixAssembleSet_2DVect  => MEF90_MassMatrixAssembleSet
   Use m_MEF_MassMatrixImplementation_MEF90Element2D_Elast, MassMatrixAssembleSet_2DElast => MEF90_MassMatrixAssembleSet
   Use m_MEF_MassMatrixImplementation_MEF90Element3D_Scal,  MassMatrixAssembleSet_3DScal  => MEF90_MassMatrixAssembleSet
   Use m_MEF_MassMatrixImplementation_MEF90Element3D_Vect,  MassMatrixAssembleSet_3DVect  => MEF90_MassMatrixAssembleSet
   Use m_MEF_MassMatrixImplementation_MEF90Element3D_Elast, MassMatrixAssembleSet_3DElast => MEF90_MassMatrixAssembleSet
   Use petsc

   IMPLICIT NONE

   Interface MEF90_MassMatrixAssembleSet
      Module Procedure MassMatrixAssembleSet_2DScal, MassMatrixAssembleSet_2DVect, MassMatrixAssembleSet_2DElast, &
                       MassMatrixAssembleSet_3DScal, MassMatrixAssembleSet_3DVect, MassMatrixAssembleSet_3DElast 
   End Interface MEF90_MassMatrixAssembleSet
End Module m_MEF_MassMatrixInterface

