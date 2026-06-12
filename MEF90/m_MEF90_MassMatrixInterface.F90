module m_MEF90_MassMatrixInterface
#include "petsc/finclude/petsc.h"
   use m_MEF90_MassMatrixImplementation_MEF90Element2DScal, only: MassMatrixAssembleSet2DScal => MEF90_MassMatrixAssembleSet
   use m_MEF90_MassMatrixImplementation_MEF90Element2DVect, only: MassMatrixAssembleSet2DVect => MEF90_MassMatrixAssembleSet
   use m_MEF90_MassMatrixImplementation_MEF90Element3DScal, only: MassMatrixAssembleSet3DScal => MEF90_MassMatrixAssembleSet
   use m_MEF90_MassMatrixImplementation_MEF90Element3DVect, only: MassMatrixAssembleSet3DVect => MEF90_MassMatrixAssembleSet

   implicit none(type, external)

   private

   public :: MEF90_MassMatrixAssembleSet
   interface MEF90_MassMatrixAssembleSet
      module procedure MassMatrixAssembleSet2DScal, MassMatrixAssembleSet2DVect, &
         MassMatrixAssembleSet3DScal, MassMatrixAssembleSet3DVect
   end interface MEF90_MassMatrixAssembleSet

end module m_MEF90_MassMatrixInterface

