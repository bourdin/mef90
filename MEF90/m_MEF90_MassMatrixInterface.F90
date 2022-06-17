Module m_MEF90_MassMatrixInterface
#include "petsc/finclude/petsc.h"
   Use m_MEF90_MassMatrixImplementation_MEF90Element2D_Scal,  MassMatrixAssembleSet_2DScal  => MEF90_MassMatrixAssembleSet
   Use m_MEF90_MassMatrixImplementation_MEF90Element2D_Vect,  MassMatrixAssembleSet_2DVect  => MEF90_MassMatrixAssembleSet
   Use m_MEF90_MassMatrixImplementation_MEF90Element3D_Scal,  MassMatrixAssembleSet_3DScal  => MEF90_MassMatrixAssembleSet
   Use m_MEF90_MassMatrixImplementation_MEF90Element3D_Vect,  MassMatrixAssembleSet_3DVect  => MEF90_MassMatrixAssembleSet

   IMPLICIT NONE

   Private

   Public :: MEF90_MassMatrixAssembleSet
   Interface MEF90_MassMatrixAssembleSet
      Module Procedure MassMatrixAssembleSet_2DScal, MassMatrixAssembleSet_2DVect, &
                       MassMatrixAssembleSet_3DScal, MassMatrixAssembleSet_3DVect 
   End Interface MEF90_MassMatrixAssembleSet

End Module m_MEF90_MassMatrixInterface

