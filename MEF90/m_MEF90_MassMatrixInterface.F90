Module m_MEF90_MassMatrixInterface
#include "finclude/petscdef.h"
   Use m_MEF90_LinAlg
   Use m_MEF90_Parameters
   Use m_MEF90_Elements
   Use m_MEF90_Utils
   Use m_MEF90_MassMatrixImplementation_MEF90Element2D_Scal,  MassMatrixAssembleSet_2DScal  => MEF90_MassMatrixAssembleSet, &
                                                              MassMatrixAssembleLoc_2DScal  => MEF90_MassMatrixAssembleLoc
   Use m_MEF90_MassMatrixImplementation_MEF90Element2D_Vect,  MassMatrixAssembleSet_2DVect  => MEF90_MassMatrixAssembleSet, &
                                                              MassMatrixAssembleLoc_2DVect  => MEF90_MassMatrixAssembleLoc
   Use m_MEF90_MassMatrixImplementation_MEF90Element2D_Elast, MassMatrixAssembleSet_2DElast => MEF90_MassMatrixAssembleSet, &
                                                              MassMatrixAssembleLoc_2DElast  => MEF90_MassMatrixAssembleLoc
   Use m_MEF90_MassMatrixImplementation_MEF90Element3D_Scal,  MassMatrixAssembleSet_3DScal  => MEF90_MassMatrixAssembleSet, &
                                                              MassMatrixAssembleLoc_3DScal  => MEF90_MassMatrixAssembleLoc
   Use m_MEF90_MassMatrixImplementation_MEF90Element3D_Vect,  MassMatrixAssembleSet_3DVect  => MEF90_MassMatrixAssembleSet, &
                                                              MassMatrixAssembleLoc_3DVect  => MEF90_MassMatrixAssembleLoc
   Use m_MEF90_MassMatrixImplementation_MEF90Element3D_Elast, MassMatrixAssembleSet_3DElast => MEF90_MassMatrixAssembleSet, &
                                                              MassMatrixAssembleLoc_3DElast  => MEF90_MassMatrixAssembleLoc
   Use petsc

   IMPLICIT NONE

   Interface MEF90_MassMatrixAssembleSet
      Module Procedure MassMatrixAssembleSet_2DScal, MassMatrixAssembleSet_2DVect, MassMatrixAssembleSet_2DElast, &
                       MassMatrixAssembleSet_3DScal, MassMatrixAssembleSet_3DVect, MassMatrixAssembleSet_3DElast 
   End Interface MEF90_MassMatrixAssembleSet

   Interface MEF90_MassMatrixAssembleLoc
      Module Procedure MassMatrixAssembleLoc_2DScal, MassMatrixAssembleLoc_2DVect, MassMatrixAssembleLoc_2DElast, &
                       MassMatrixAssembleLoc_3DScal, MassMatrixAssembleLoc_3DVect, MassMatrixAssembleLoc_3DElast 
   End Interface MEF90_MassMatrixAssembleLoc
End Module m_MEF90_MassMatrixInterface

