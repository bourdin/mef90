#if defined PB_2D
Module m_SimplePoisson2D
#elif defined PB_3D
Module m_SimplePoisson3D
#endif

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscmeshdef.h"
#include "finclude/petscviewerdef.h"

   Use m_MEF90
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
      PetscInt                                     :: TestCase
      Character(len=MEF90_MXSTRLEN)                :: prefix
      Type(PetscViewer)                            :: LogViewer, MyLogViewer
   End Type AppParam_Type

   Type AppCtx_Type
      Type (MeshTopology_Type)                     :: MeshTopology
      Type (EXO_Type)                              :: EXO, MyEXO
#if defined PB_2D
      Type(Element2D_Scal), Dimension(:), Pointer  :: Elem
#elif defined PB_3D
      Type(Element3D_Scal), Dimension(:), Pointer  :: Elem
#endif
      Type(SectionReal)                            :: U
      Type(SectionReal)                            :: GradU
      Type(SectionReal)                            :: F
      PetscReal                                    :: Energy
      Type(VecScatter)                             :: Scatter
      Type(SectionInt)                             :: BCFlag
      Type(Mat)                                    :: K
      Type(Vec)                                    :: RHS
      Type(KSP)                                    :: KSP
      Type(PC)                                     :: PC
      Type(LogInfo_Type)                           :: LogInfo
      Type(AppParam_Type)                          :: AppParam
   End Type AppCtx_Type
   
   
Contains

   Subroutine ElastInit(AppCtx)
      !!! BB
   End Subroutine ElastInit
   
   Subroutine Solve(AppCtx)
      !!! CM
   End Subroutine Solve
   
   Subroutine MatAssembly(AppCtx)
      !!! CM
   End Subroutine MatAssembly
   
   Subroutine MatAssemblyLocal(iE, Mat_Prop, AppCtx, MatElem)
      !!! CM
   End Subroutine MatAssemblyLocal

   Subroutine RHSAssembly(AppCtx)
      !!! CM
   End Subroutine RHSAssembly
   
   Subroutine RHSAssemblyLocal(iE, Map_Prop, AppCtx, RHSElem)
      !!! CM
   End Subroutine RHSAssemblyLocal
   
   Subroutine ComputeEnergy(AppCtx)
      !!! CM
   End Subroutine ComputeEnergy
   
   Subroutine Save(AppCtx)
      !!! BB
   End Subroutine Save
   
   Subroutine ElastFinalize(AppCtx)
      !!! BB
   End Subroutine ELastFinalize
   
   
#if defined PB_2D
End Module m_SimplePoisson2D
#elif defined PB_3D
End Module m_SimplePoisson3D
#endif
