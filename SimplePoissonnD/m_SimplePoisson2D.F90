#if defined PB_2D
Module m_SimplePoisson2D
#elif defined PB_3D
Module m_SimplePoisson3D
#endif

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscmeshdef.h"

   Use m_MEF90
   Use petsc
   Use petscvec
   Use petscmat
   Use petscmesh

   Implicit NONE   
   
   Type LogInfo_Type
      PetscLogStage               :: IO_Stage
      PetscLogStage               :: Distribute_Stage
      PetscLogStage               :: MatAssembly_Stage
      PetscLogStage               :: RHSAssembly_Stage
      PetscLogStage               :: KSPSolve_Stage
      PetscLogStage               :: EnergyEval_Stage
      
      PetscLogEvent               :: MatAssemblyLocal_Event
      PetscLogEvent               :: RHSAssemblyLocal_Event
      PetscLogEvent               :: EnergyEval_Event
   End Type LogInfo_Type

   Type AppCtx_Type
      Type (MeshTopology_Info)                     :: MeshTopology
      Type (EXO_Info)                              :: EXO, MyEXO
#if defined PB_2D
      Type(Element2D_Scal), Dimension(:), Pointer  :: Elem
#elif defined PB_3D
      Type(Element3D_Scal), Dimension(:), Pointer  :: Elem
#endif
      Type(SectionReal)                            :: U
      Type(SectionInt)                             :: BCFlag
      Type(LogInfo_Type)                           :: LogInfo
      Type(Mat)                                    :: K
      Type(Vec)                                    :: RHS
   End Type AppCtx_Type
   
Contains

   Subroutine SimplePoissonInit(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      PetscInt                                     :: iBlk, iDoF      
      PetscTruth                                   :: HasPrefix
      PetscTruth                                   :: verbose
      Character(len=MXSTLN)                        :: prefix
      PetscInt, Dimension(:), Pointer              :: TmpFlag
      PetscInt                                     :: TmpPoint
      
      Type(SectionReal)                            :: CoordSection
      PetscReal, Dimension(:), Pointer             :: TmpCoords
      PetscReal, Dimension(:,:), Pointer           :: Coords
      PetscInt                                     :: iE, iELoc


      Call MEF90_Initialize()
      Call PetscOptionsHasName(PETSC_NULL_CHARACTER, '-verbose', verbose, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p', prefix, HasPrefix, iErr); CHKERRQ(iErr)
      
      Call InitLog(AppCtx)

      !!! Read and partition the mesh
      Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
      AppCtx%EXO%Comm = PETSC_COMM_WORLD
      AppCtx%EXO%filename = Trim(prefix)//'.gen'

      Call MeshTopologyReadEXO(AppCtx%MeshTopology, AppCtx%EXO)
      !!! Split in order to be able to take the distribute out if necessary
   
      !!! Sets the type of elements for each block
      Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_Type = MEF90_P1_Lagrange
         Call Init_Elem_Blk_Info(AppCtx%MeshTopology%Elem_Blk(iBlk), AppCtx%MeshTopology%num_dim)
      End Do
   
      !!! Allocate the elements
      Allocate(AppCtx%Elem(AppCtx%MeshTopology%Num_Elems))

      !!! Initialize the Basis Functions in each element
      Call MeshGetSectionReal(AppCtx%MeshTopology%mesh, 'coordinates', CoordSection, iErr); CHKERRQ(iErr)
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Allocate(TmpCoords(AppCtx%MeshTopology%Num_Dim * AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF))
         Allocate(Coords(AppCtx%MeshTopology%Num_Dim, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF))
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, CoordSection, iE-1, Size(TmpCoords), TmpCoords, iErr); CHKERRQ(iErr)
            Coords = Reshape(TmpCoords, (/AppCtx%MeshTopology%Num_Dim, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF /) )
            Call Init_Element(AppCtx%Elem(iE), Coords, 4, AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_Type)
         End Do Do_Elem_iE
         DeAllocate(TmpCoords)
         DeAllocate(Coords)
      End Do Do_Elem_iBlk
      Call SectionRealDestroy(CoordSection, iErr); CHKERRQ(iErr)

      !!! Prepare and format the output mesh
      !!! 1. Geometry
      AppCtx%MyEXO%comm = PETSC_COMM_SELF
      AppCtx%MyEXO%exoid = AppCtx%EXO%exoid
      Write(AppCtx%MyEXO%filename, 200) trim(prefix), MEF90_MyRank
   200 Format(A, '-', I4.4, '.gen')
      AppCtx%MyEXO%title = trim(AppCtx%EXO%title)
      AppCtx%MyEXO%Num_QA = AppCtx%EXO%Num_QA
      Call Write_MeshTopologyGlobal(AppCtx%MeshTopology, AppCtx%MyEXO, PETSC_COMM_WORLD)
      
      !!! 2. Variables
      AppCtx%MyEXO%exoid = EXOPEN(AppCtx%MyEXO%filename, EXWRIT, exo_cpu_ws, exo_io_ws, exo_ver, iErr)
      Call EXPVP (AppCtx%MyEXO%exoid, 'g', 1, iErr)
      Call EXPVAN(AppCtx%MyEXO%exoid, 'g', 1, (/'Energy'/), iErr)
      Call EXPVP (AppCtx%MyEXO%exoid, 'n', 1, iErr)
      Call EXPVAN(AppCtx%MyEXO%exoid, 'n', 1, (/'U'/), iErr)
#if defined PB_2D
      Call EXPVP (AppCtx%MyEXO%exoid, 'e', 2, iErr)
      Call EXPVAN(AppCtx%MyEXO%exoid, 'e', 2, (/'Grad U_X', 'Grad U_Y'/), iErr)
#elif defined PB_3D
      Call EXPVP (AppCtx%MyEXO%exoid, 'e', 3, iErr)
      Call EXPVAN(AppCtx%MyEXO%exoid, 'e', 3, (/'Grad U_X', 'Grad U_Y', 'Grad U_Z'/), iErr)
#endif
      Call EXCLOS(AppCtx%MyEXO%exoid, iErr)

      Call PetscLogStagePop(iErr); CHKERRQ(iErr)

      !!! Allocate the Section for U
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 1, AppCtx%U, iErr); CHKERRQ(iErr)
      
      !!! Allocate and initialize the Section for the flag
      Call MeshGetVertexSectionInt(AppCtx%MeshTopology%mesh, 1, AppCtx%BCFlag, iErr); CHKERRQ(iErr)
      Allocate(TmpFlag(1))
      Do iBlk = 1, AppCtx%MeshTopology%num_node_sets  
         Do iDoF = 1, AppCtx%MeshTopology%node_set(iBlk)%Num_Nodes     
            TmpPoint = AppCtx%MeshTopology%Num_Elems + AppCtx%MeshTopology%Node_Set(iBlk)%Node_ID(iDoF)-1
            TmpFlag = 1
            Call MeshUpdateClosureInt(AppCtx%MeshTopology%Mesh, AppCtx%BCFlag, TmpPoint, TmpFlag, iErr); CHKERRQ(iErr)
         End Do
      End Do
      DeAllocate(TmpFlag)

      !!! Initialize the matrix and vector for the linear system
      Call MeshSetMaxDof(AppCtx%MeshTopology%Mesh, 1, iErr); CHKERRQ(iErr) 
      !Max DoF per point is 1 (Should it be 3?)
      Call MeshCreateMatrix(AppCtx%MeshTopology%mesh, AppCtx%U, MATMPIAIJ, AppCtx%K, iErr); CHKERRQ(iErr)
      Call MeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%U, AppCtx%RHS, iErr); CHKERRQ(iErr)
   End Subroutine SimplePoissonInit
   
   Subroutine InitLog(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr
      
      Call PetscLogEventRegister('MatAssembly Local', 0, AppCtx%LogInfo%MatAssemblyLocal_Event, ierr); CHKERRQ(ierr)
      Call PetscLogEventRegister('RHSAssembly Local', 0, AppCtx%LogInfo%RHSAssemblyLocal_Event, ierr); CHKERRQ(ierr)
      Call PetscLogEventRegister('Energy Eval',       0, AppCtx%LogInfo%EnergyEval_Event,       ierr); CHKERRQ(ierr)

      Call PetscLogStageRegister("IO Stage",          AppCtx%LogInfo%IO_Stage,     iErr)
      Call PetscLogStageRegister("Mesh Distribution", AppCtx%LogInfo%Distribute_Stage,  iErr)
      Call PetscLogStageRegister("Mat Assembly",      AppCtx%LogInfo%MatAssembly_Stage, iErr)
      Call PetscLogStageRegister("RHS Assembly",      AppCtx%LogInfo%RHSAssembly_Stage, iErr)
      Call PetscLogStageRegister("KSP Solve",         AppCtx%LogInfo%KSPSolve_Stage,    iErr)
      Call PetscLogStageRegister("Energy Evaluation", AppCtx%LogInfo%EnergyEval_Stage,  iErr)
   End Subroutine InitLog
   
   Subroutine MatAssembly(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iBlk, iE, iELoc, iErr
      PetscReal, Dimension(:,:), Pointer           :: MatElem
      
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Allocate(MatElem(AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF))
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            Call MatAssemblyLocal(iE, AppCtx, MatElem)
            Call assembleMatrix(AppCtx%K, AppCtx%MeshTopology%mesh, AppCtx%U, iE-1, MatElem, ADD_VALUES, iErr); CHKERRQ(iErr)
         End Do Do_Elem_iE
         DeAllocate(MatElem)
      End Do Do_Elem_iBlk
      Call MatAssemblyBegin(AppCtx%K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (AppCtx%K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)

   End Subroutine MatAssembly
   
   
   Subroutine MatAssemblyLocal(iE, AppCtx, MatElem)
      Type(AppCtx_Type)                            :: AppCtx
      PetscReal, Dimension(:,:), Pointer           :: MatElem 
      PetscInt                                     :: iE
   
      PetscInt                                     :: NumDoF, NumGauss
      PetscInt, Dimension(:), Pointer              :: BCFlag
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      PetscInt                                     :: iErr
      
      Call PetscLogEventBegin(AppCtx%LogInfo%MatAssemblyLocal_Event, iErr); CHKERRQ(iErr)
      
      MatElem  = 0.0_Kr
      NumDoF   = Size(AppCtx%Elem(iE)%BF,1)
      NumGauss = Size(AppCtx%Elem(iE)%BF,2)
      Allocate(BCFlag(NumDoF))
      Call MeshRestrictClosureInt(AppCtx%MeshTopology%mesh, AppCtx%BCFlag, iE-1, NumDoF, BCFlag, iErr); CHKERRQ(ierr)
      Do iGauss = 1, NumGauss
         Do iDoF1 = 1, NumDoF
            If (BCFlag(iDoF1) == 0) Then
               Do iDoF2 = 1, NumDoF
                  MatElem(iDoF2, iDoF1) = AppCtx%Elem(iE)%Gauss_C(iGauss) * AppCtx%Elem(iE)%Grad_BF(iDoF1, iGauss) .DotP. AppCtx%Elem(iE)%Grad_BF(iDoF2, iGauss) 
               End Do
            End If
         End Do
      End Do
      Call PetscLogFlops(1, iErr);CHKERRQ(iErr)
      
      DeAllocate(BCFlag)
      Call PetscLogEventEnd(AppCtx%LogInfo%MatAssemblyLocal_Event, iErr); CHKERRQ(iErr)
   End Subroutine MatAssemblyLocal
   
   Subroutine SimplePoissonFinalize(AppCtx)   
      Type(AppCtx_Type)                            :: AppCtx

      PetscInt                                     :: iErr

      Call SectionRealDestroy(AppCtx%U, iErr); CHKERRQ(iErr)
      Call SectionIntDestroy(AppCtx%BCFlag, iErr); CHKERRQ(iErr)
      Call MeshDestroy(AppCtx%MeshTopology%Mesh, iErr); CHKERRQ(ierr)
      Call MEF90_Finalize()
   End Subroutine SimplePoissonFinalize

#if defined PB_2D
End Module m_SimplePoisson2D
#elif defined PB_3D
End Module m_SimplePoisson3D
#endif
