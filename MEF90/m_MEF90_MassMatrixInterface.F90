Module m_MEF90_MassMatrixInterface
#include "petsc/finclude/petsc.h"
   Use m_MEF90_MassMatrixImplementation_MEF90Element2DScal,  MassMatrixAssembleSet2DScal  => MEF90_MassMatrixAssembleSet
   Use m_MEF90_MassMatrixImplementation_MEF90Element2DVect,  MassMatrixAssembleSet2DVect  => MEF90_MassMatrixAssembleSet
   Use m_MEF90_MassMatrixImplementation_MEF90Element3DScal,  MassMatrixAssembleSet3DScal  => MEF90_MassMatrixAssembleSet
   Use m_MEF90_MassMatrixImplementation_MEF90Element3DVect,  MassMatrixAssembleSet3DVect  => MEF90_MassMatrixAssembleSet

   IMPLICIT NONE

   Private

   Public :: MEF90_MassMatrixAssembleSet
   Interface MEF90_MassMatrixAssembleSet
      Module Procedure MassMatrixAssembleSet2DScal, MassMatrixAssembleSet2DVect, &
                       MassMatrixAssembleSet3DScal, MassMatrixAssembleSet3DVect 
   End Interface MEF90_MassMatrixAssembleSet

End Module m_MEF90_MassMatrixInterface

