#if defined PB_2D
Module m_AmbrosioTortorelli_Types2D
#elif defined PB_3D
Module m_AmbrosioTortorelli_Types3D
#endif

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscmeshdef.h"
#include "finclude/petscviewerdef.h"

   Use m_MEF90
   Use m_RuptStruct
   Use petsc
   Use petscvec
   Use petscmat
   Use petscksp
   Use petscmesh

   Implicit NONE   
   
   Type LogInfo_Type
      PetscLogStage               :: IO_Stage
      PetscLogStage               :: Distribute_Stage
      PetscLogStage               :: DataSetup_Stage
      PetscLogStage               :: MatAssembly_Stage
      PetscLogStage               :: RHSAssembly_Stage
      PetscLogStage               :: KSPSolve_Stage
      PetscLogStage               :: PostProc_Stage
      
      PetscLogEvent               :: MatAssemblyLocal_Event
      PetscLogEvent               :: RHSAssemblyLocal_Event
      PetscLogEvent               :: PostProc_Event
   End Type LogInfo_Type

   Type AppParam_Type
      PetscTruth                                   :: Restart
      PetscTruth                                   :: Verbose
      Character(len=MEF90_MXSTRLEN)                :: prefix
      Type(PetscViewer)                            :: LogViewer, MyLogViewer
   End Type AppParam_Type

   Type AppCtx_Type
      Type (MeshTopology_Type)                     :: MeshTopology
      Type (EXO_Type)                              :: EXO, MyEXO
#if defined PB_2D
      Type(Element2D_Elast), Dimension(:), Pointer :: ElemVect
      Type(Element2D_Scal), Dimension(:), Pointer  :: ElemScal
#elif defined PB_3D
      Type(Element3D_Elast), Dimension(:), Pointer :: ElemVect
      Type(Element3D_Scal), Dimension(:), Pointer  :: ElemScal
#endif
      Type(SectionReal)                            :: U
      Type(SectionReal)                            :: StressU
      Type(SectionReal)                            :: StrainU
      Type(SectionReal)                            :: F
      Type(SectionReal)                            :: V
      Type(SectionReal)                            :: GradV
      Type(SectionReal)                            :: Theta
      PetscReal                                    :: Load
      PetscInt                                     :: TimeStep
      PetscReal                                    :: SurfaceEnergy
      PetscReal                                    :: ElasticEnergy
      PetscReal                                    :: ExtForcesWork
      PetscReal                                    :: TotalEnergy
      PetscReal                                    :: ErrV
      Type(VecScatter)                             :: ScatterVect
      Type(VecScatter)                             :: ScatterScal
      Type(SectionInt)                             :: BCFlagU, BCFlagV
      Type(Mat)                                    :: KU, KV
      Type(Vec)                                    :: RHSU, RHSV
      Type(KSP)                                    :: KSPU, KSPV
      Type(PC)                                     :: PCU, PCV
      Type(LogInfo_Type)                           :: LogInfo
#if defined PB_2D
      Type(MatProp2D_Type), Dimension(:), Pointer  :: MatProp      
#elif defined PB_3D
      Type(MatProp3D_Type), Dimension(:), Pointer  :: MatProp
#endif
      Type(AppParam_Type)                          :: AppParam
      Type(RuptSchemeParam_Type)                   :: RuptSchemeParam
   End Type AppCtx_Type
      
   
#if defined PB_2D
End Module m_AmbrosioTortorelli_Types2D
#elif defined PB_3D
End Module m_AmbrosioTortorelli_Types3D
#endif
