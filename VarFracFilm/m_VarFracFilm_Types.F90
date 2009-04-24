Module m_VarFracFilm_Types

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscmeshdef.h"
#include "finclude/petscviewerdef.h"

   Use m_MEF90
   Use m_VarFracFilm_Struct
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
      Type(Element2D_Elast), Dimension(:), Pointer :: ElemVect
      Type(Element2D_Scal), Dimension(:), Pointer  :: ElemScal
      Type(SectionReal)                            :: U
      Type(SectionReal)                            :: StressU
      Type(SectionReal)                            :: StrainU
      Type(SectionReal)                            :: V
      Type(SectionReal)                            :: PHI
      Type(SectionReal)                            :: GradV
      Type(SectionReal)                            :: Theta
      Type(SectionReal)                            :: U0
      PetscReal                                    :: Load
      PetscInt                                     :: TimeStep
      PetscReal                                    :: SurfaceEnergyT
      PetscReal                                    :: SurfaceEnergyD
      PetscReal                                    :: ElasticBulkEnergy
      PetscReal                                    :: ElasticInterEnergy
      PetscReal                                    :: TotalEnergy
      PetscReal                                    :: ErrV
      PetscInt                                     :: ErrPHI
      Type(VecScatter)                             :: ScatterVect
      Type(VecScatter)                             :: ScatterScal
      Type(SectionInt)                             :: BCFlagU, BCFlagV
      Type(Mat)                                    :: KU, KV
      Type(Vec)                                    :: RHSU, RHSV
      Type(KSP)                                    :: KSPU, KSPV
      Type(PC)                                     :: PCU, PCV
      Type(LogInfo_Type)                           :: LogInfo
      Type(MatProp2D_Type), Dimension(:), Pointer  :: MatProp      
      Type(AppParam_Type)                          :: AppParam
      Type(VarFracFilmSchemeParam_Type)            :: VarFracFilmSchemeParam
   End Type AppCtx_Type
      
   
End Module m_VarFracFilm_Types