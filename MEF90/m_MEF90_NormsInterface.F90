module m_MEF90_NormsInterface
#include "petsc/finclude/petsc.h"
   use m_MEF90_NormsImplementation_MEF90Element2DScal, MEF90L2DotProductSet2DScal => MEF90L2DotProductSet, &
      MEF90H1DotProductSet2DScal => MEF90H1DotProductSet, &
      MEF90L2NormSet2DScal => MEF90L2NormSet
   use m_MEF90_NormsImplementation_MEF90Element2DVect, MEF90L2DotProductSet2DVect => MEF90L2DotProductSet, &
      MEF90H1DotProductSet2DVect => MEF90H1DotProductSet, &
      MEF90H1SymDotProductSet2DVect => MEF90H1SymDotProductSet, &
      MEF90L2NormSet2DVect => MEF90L2NormSet
   use m_MEF90_NormsImplementation_MEF90Element3DScal, MEF90L2DotProductSet3DScal => MEF90L2DotProductSet, &
      MEF90H1DotProductSet3DScal => MEF90H1DotProductSet, &
      MEF90L2NormSet3DScal => MEF90L2NormSet
   use m_MEF90_NormsImplementation_MEF90Element3DVect, MEF90L2DotProductSet3DVect => MEF90L2DotProductSet, &
      MEF90H1DotProductSet3DVect => MEF90H1DotProductSet, &
      MEF90H1SymDotProductSet3DVect => MEF90H1SymDotProductSet, &
      MEF90L2NormSet3DVect => MEF90L2NormSet

   implicit none(type, external)

   private

   public :: MEF90L2DotProductSet, MEF90H1DotProductSet, MEF90H1symDotProductSet, MEF90L2NormSet

   interface MEF90L2DotProductSet
      module procedure MEF90L2DotProductSet2DScal, MEF90L2DotProductSet2DVect, &
         MEF90L2DotProductSet3DScal, MEF90L2DotProductSet3DVect
   end interface MEF90L2DotProductSet

   interface MEF90H1DotProductSet
      module procedure MEF90H1DotProductSet2DScal, MEF90H1DotProductSet2DVect, &
         MEF90H1DotProductSet3DScal, MEF90H1DotProductSet3DVect
   end interface MEF90H1DotProductSet

   interface MEF90H1symDotProductSet
      module procedure MEF90H1symDotProductSet2DVect, &
         MEF90H1symDotProductSet3DVect
   end interface MEF90H1symDotProductSet

   interface MEF90L2NormSet
      module procedure MEF90L2NormSet2DScal, MEF90L2NormSet2DVect, &
         MEF90L2NormSet3DScal, MEF90L2NormSet3DVect
   end interface MEF90L2NormSet
end module m_MEF90_NormsInterface

