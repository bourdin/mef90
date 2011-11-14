Module m_VarFilmQS_Types
#include "finclude/petscdef.h"

   Use m_MEF90
   Use m_Film_Struct

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
      PetscLogStage               :: FWAssemblyW_Stage
      PetscLogStage               :: WStep_Stage
      PetscLogStage               :: PostProc_Stage
   
      PetscLogEvent               :: MatAssemblyLocalU_Event
      PetscLogEvent               :: RHSAssemblyLocalU_Event
      PetscLogEvent               :: MatAssemblyLocalV_Event
      PetscLogEvent               :: RHSAssemblyLocalV_Event
      PetscLogEvent               :: FWAssemblyLocalW_Event
      PetscLogEvent               :: PostProc_Event
   End Type LogInfo_Type

   Type AppParam_Type
      PetscBool                                    :: Restart
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
      Type(Element2D_Elast), Dimension(:), Pointer :: ElemVect
      Type(Element2D_Scal), Dimension(:), Pointer  :: ElemScal
      Type(Field)                                  :: U
      Type(Field)                                  :: U0
      Type(Field)                                  :: UBC
      Type(Field)                                  :: V
      Type(Field)                                  :: W
      Type(Field)                                  :: WBC, WIrrev
      Type(Field)                                  :: VBC, VIrrev
      Type(Field)                                  :: Theta
      Type(Field)                                  :: RHSU, GradientU, LowerBoundU, UpperBoundU
      Type(Field)                                  :: RHSV, GradientV, LowerBoundV, UpperBoundV
   Type(Field)                                  :: FW
   Type(Vec)                                    :: V_Old
   Type(Vec)                                    :: W_Old
   Type(Flag)                                   :: BCUFlag, BCVFlag, BCWFlag, IrrevFlag, WIrrevFlag
      Type(SectionReal)                            :: StressU
      Type(SectionReal)                            :: StrainU
      PetscInt                                     :: NumTimeSteps
      PetscInt                                     :: TimeStep
      PetscReal, Dimension(:), Pointer             :: Load                 ! All Time Steps
   PetscReal, Dimension(:), Pointer             :: FractureEnergy        ! All Time Steps
   PetscReal, Dimension(:), Pointer             :: DelaminationEnergy        ! All Time Steps
   PetscReal, Dimension(:), Pointer             :: CohesiveEnergy        ! All Time Steps
      PetscReal, Dimension(:), Pointer             :: ElasticEnergy        ! All Time Steps
      PetscReal, Dimension(:), Pointer             :: ExtForcesWork        ! All Time Steps
      PetscReal, Dimension(:), Pointer             :: TotalEnergy
   PetscReal, Dimension(:), Pointer             :: FractureEnergyBlock   ! Current TS, All Blocks
   PetscReal, Dimension(:), Pointer             :: DelaminationEnergyBlock   ! Current TS, All Blocks
   PetscReal, Dimension(:), Pointer             :: CohesiveEnergyBlock   ! Current TS, All Blocks
      PetscReal, Dimension(:), Pointer             :: ElasticEnergyBlock   ! Current TS, All Blocks
      PetscReal, Dimension(:), Pointer             :: ExtForcesWorkBlock   ! Current TS, All Blocks
      PetscReal, Dimension(:), Pointer             :: TotalEnergyBlock     ! Current TS, All Blocks
      PetscReal                                    :: ErrV
      Type(Mat)                                    :: KU, KV, KW
      Type(KSP)                                    :: KSPU, KSPV, KSPW
      Type(PC)                                     :: PCU, PCV, PCW
      Type(LogInfo_Type)                           :: LogInfo
      PetscBool                                    :: IsBT
#if defined WITH_TAO
      TAO_SOLVER                                   :: taoV
      TAO_APPLICATION                              :: taoAppV
#endif
      
      Type(MatProp2D_Type), Dimension(:), Pointer  :: MatProp
      Type(AppParam_Type)                          :: AppParam
      Type(VarFracSchemeParam_Type)                :: VarFracSchemeParam
   End Type AppCtx_Type
      
   
End Module m_VarFilmQS_Types
