#if defined PB_2D
Module m_VarFracQS_Types2D
#elif defined PB_3D
Module m_VarFracQS_Types3D
#endif
#include "finclude/petscdef.h"

   Use m_MEF90
   Use m_VarFrac_Struct

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
      PetscLogStage               :: UStep_Stage
      PetscLogStage               :: MatAssemblyV_Stage
      PetscLogStage               :: RHSAssemblyV_Stage
      PetscLogStage               :: VStep_Stage
      PetscLogStage               :: PostProc_Stage
      
      PetscLogEvent               :: MatAssemblyLocalU_Event
      PetscLogEvent               :: RHSAssemblyLocalU_Event
      PetscLogEvent               :: MatAssemblyLocalV_Event
      PetscLogEvent               :: RHSAssemblyLocalV_Event
      PetscLogEvent               :: PostProc_Event
   End Type LogInfo_Type

   Type AppParam_Type
      PetscTruth                                   :: Restart
      PetscInt                                     :: Verbose
      PetscBool                                    :: StopOnError
      Character(len=MEF90_MXSTRLEN)                :: prefix
      Character(len=MEF90_MXSTRLEN)                :: CST_FileName
      Character(len=MEF90_MXSTRLEN)                :: Ener_FileName
      Integer                                      :: Ener_Unit
      Type(PetscViewer)                            :: LogViewer, MyLogViewer
      Character(len=MEF90_MXSTRLEN)                :: Log_FileName, MyLog_FileName
      Integer, Dimension(:), Pointer               :: EnerBlock_Unit
      Character(len=MEF90_MXSTRLEN)                :: EnerBlock_Suffix
      Character(len=MEF90_MXSTRLEN)                :: EnerBlock_Filename
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
      Type(Field)                                  :: U
      Type(Field)                                  :: UBC
      Type(Field)                                  :: F
      Type(Field)                                  :: V
      Type(Field)                                  :: VBC, VIrrev
      Type(Field)                                  :: Theta
      Type(Field)                                  :: RHSU, GradientU, LowerBoundU, UpperBoundU
      Type(Field)                                  :: RHSV, GradientV, LowerBoundV, UpperBoundV
      Type(Vec)                                    :: V_Old
      Type(Flag)                                   :: BCUFlag, BCVFlag, IrrevFlag
      Type(SectionReal)                            :: StressU
      Type(SectionReal)                            :: StrainU
      PetscInt                                     :: NumTimeSteps
      PetscInt                                     :: TimeStep
      PetscReal, Dimension(:), Pointer             :: Load                 ! All Time Steps
      PetscReal, Dimension(:), Pointer             :: SurfaceEnergy        ! All Time Steps
      PetscReal, Dimension(:), Pointer             :: ElasticEnergy        ! All Time Steps
      PetscReal, Dimension(:), Pointer             :: ExtForcesWork        ! All Time Steps
      PetscReal, Dimension(:), Pointer             :: TotalEnergy
      PetscReal, Dimension(:), Pointer             :: SurfaceEnergyBlock   ! Current TS, All Blocks
      PetscReal, Dimension(:), Pointer             :: ElasticEnergyBlock   ! Current TS, All Blocks
      PetscReal, Dimension(:), Pointer             :: ExtForcesWorkBlock   ! Current TS, All Blocks
      PetscReal, Dimension(:), Pointer             :: TotalEnergyBlock     ! Current TS, All Blocks
      PetscReal                                    :: ErrV
      Type(Mat)                                    :: KU, KV
      Type(KSP)                                    :: KSPU, KSPV
      Type(PC)                                     :: PCU, PCV
      Type(LogInfo_Type)                           :: LogInfo
      PetscBool                                    :: IsBT
#if defined WITH_TAO
      TAO_SOLVER                                   :: taoV
      TAO_APPLICATION                              :: taoAppV
      TAO_SOLVER                                   :: taoU
      TAO_APPLICATION                              :: taoAppU
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
