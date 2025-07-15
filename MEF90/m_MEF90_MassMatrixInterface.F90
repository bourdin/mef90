module m_MEF90_MassMatrixInterface
#include "petsc/finclude/petsc.h"
   use m_MEF90_MassMatrixImplementation_MEF90Element2DScal, MassMatrixAssembleSet2DScal => MEF90_MassMatrixAssembleSet
   use m_MEF90_MassMatrixImplementation_MEF90Element2DVect, MassMatrixAssembleSet2DVect => MEF90_MassMatrixAssembleSet
   use m_MEF90_MassMatrixImplementation_MEF90Element3DScal, MassMatrixAssembleSet3DScal => MEF90_MassMatrixAssembleSet
   use m_MEF90_MassMatrixImplementation_MEF90Element3DVect, MassMatrixAssembleSet3DVect => MEF90_MassMatrixAssembleSet

   implicit none(type, external)

   private

   public :: MEF90_MassMatrixAssembleSet
   interface MEF90_MassMatrixAssembleSet
      module procedure MassMatrixAssembleSet2DScal, MassMatrixAssembleSet2DVect, &
         MassMatrixAssembleSet3DScal, MassMatrixAssembleSet3DVect
   end interface MEF90_MassMatrixAssembleSet

end module m_MEF90_MassMatrixInterface

