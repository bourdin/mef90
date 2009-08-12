#if defined PB_2D
Module m_VarFracQS_Types2D
#elif defined PB_3D
Module m_VarFracQS_Types3D
#endif

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscmeshdef.h"
#include "finclude/petscviewerdef.h"

   Use m_MEF90
   Use m_VarFrac_Struct
   Use petsc
   Use petscvec
   Use petscmat
   Use petscksp
   Use petscmesh

   Implicit NONE   
   Private
#if defined WITH_TAO
#include "include/finclude/tao_solver.h"
#endif

   Public :: LogInfo_Type
   Public :: AppParam_Type
   Public :: AppCtx_Type

   Type LogInfo_Type
      PetscLogStage               :: IO_Stage
      PetscLogStage               :: Setup_Stage      
      PetscLogStage               :: MeshDistribute_Stage
      PetscLogStage               :: MatAssemblyU_Stage
      PetscLogStage               :: RHSAssemblyU_Stage
      PetscLogStage               :: KSPSolveU_Stage
      PetscLogStage               :: MatAssemblyV_Stage
      PetscLogStage               :: RHSAssemblyV_Stage
      PetscLogStage               :: KSPSolveV_Stage
      PetscLogStage               :: PostProc_Stage
      
      PetscLogEvent               :: MatAssemblyLocalU_Event
      PetscLogEvent               :: RHSAssemblyLocalU_Event
      PetscLogEvent               :: MatAssemblyLocalV_Event
      PetscLogEvent               :: RHSAssemblyLocalV_Event
      PetscLogEvent               :: PostProc_Event
   End Type LogInfo_Type

   Type AppParam_Type
      PetscInt                                     :: Restart
      PetscInt                                     :: Verbose
      Character(len=MEF90_MXSTRLEN)                :: prefix
      Type(PetscViewer)                            :: LogViewer, MyLogViewer
      Integer                                      :: Ener_Unit
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
      Type(SectionReal)                            :: Theta
      PetscInt                                     :: NumTimeSteps
      PetscInt                                     :: TimeStep
      PetscReal, Dimension(:), Pointer             :: Load
      PetscReal, Dimension(:), Pointer             :: SurfaceEnergy
      PetscReal, Dimension(:), Pointer             :: ElasticEnergy
      PetscReal, Dimension(:), Pointer             :: ExtForcesWork
      PetscReal, Dimension(:), Pointer             :: TotalEnergy
      PetscReal                                    :: ErrV
      Type(VecScatter)                             :: ScatterVect
      Type(VecScatter)                             :: ScatterScal
      Type(SectionInt)                             :: BCUFlag, BCVFlag, IrrevFlag
      Type(Mat)                                    :: KU, KV
      Type(SectionReal)                            :: RHSV
      Type(Vec)                            :: RHSU
      Type(KSP)                                    :: KSPU, KSPV
      Type(PC)                                     :: PCU, PCV
      Type(LogInfo_Type)                           :: LogInfo
      PetscTruth                                   :: IsBT
#if defined WITH_TAO
      TAO_SOLVER                                   :: taoV
      TAO_APPLICATION                              :: taoAppV
#endif
      
#if defined PB_2D
      Type(MatProp2D_Type), Dimension(:), Pointer  :: MatProp      
#elif defined PB_3D
      Type(MatProp3D_Type), Dimension(:), Pointer  :: MatProp
#endif
      Type(AppParam_Type)                          :: AppParam
      Type(VarFracSchemeParam_Type)                :: VarFracSchemeParam
   End Type AppCtx_Type
      
   
#if defined PB_2D
End Module m_VarFracQS_Types2D
#elif defined PB_3D
End Module m_VarFracQS_Types3D
#endif
