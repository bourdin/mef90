#if defined PB_2D  
Module m_VarFracQS_V2D
#elif defined PB_3D
Module m_VarFracQS_V3D
#endif

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscmeshdef.h"
#include "finclude/petscviewerdef.h"

#if defined PB_2D
   Use m_VarFracQS_Types2D
   Use m_VarFracQS_Post2D
#elif defined PB_3D
   Use m_VarFracQS_Types3D
   Use m_VarFracQS_Post3D
#endif   
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
   Public :: HessianV_Assembly
   Public :: FormFunctionAndGradientV
   Public :: InitTaoBoundsV
#endif

   Public :: Init_TS_V
   Public :: Init_TS_Irrev
   Public :: MatV_Assembly
   Public :: RHSV_Assembly
   Public :: Step_V

Contains

#if defined WITH_TAO
   Subroutine InitTaoBoundsV(TaoApp, LowerBoundV_Vec, UpperBoundV_Vec, AppCtx, iErr)
      TAO_APPLICATION                              :: taoapp
      Type(Vec)                                    :: LowerBoundV_Vec, UpperBoundV_Vec
      Type(AppCtx_Type)                            :: AppCtx
      PetscErrorCode                               :: iErr
      
      Type(SectionReal)                            :: LowerBoundV_Sec, UpperBoundV_Sec
      PetscInt, Dimension(:), Pointer              :: BCVFlag_Ptr, IrrevFlag_Ptr
      PetscReal, Dimension(:), Pointer             :: Bound_Ptr, V_Ptr
      PetscInt                                     :: iDoF, NumDoF
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer      
      
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) "Updating bounds for V\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'LowerBoundV_Sec',  1, LowerBoundV_Sec , iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'UpperBoundV_Sec',  1, UpperBoundV_Sec , iErr); CHKERRQ(iErr) 
      Call SectionRealSet(LowerBoundV_Sec, 0.0_Kr, iErr); CHKERRQ(iErr)
      Call SectionRealSet(UpperBoundV_Sec, 1.0_Kr, iErr); CHKERRQ(iErr)
      
      Select Case(AppCtx%VarFracSchemeParam%IrrevType)
      Case (VarFrac_Irrev_Eq)
         !!! Set Lower bound and upper bounds to 0 at each point where BCFlag == VarFrac_BC_Type_DIRI and to V whereever IrrevFlag = VarFrac_BC_Type_DIRI
         NumDoF = AppCtx%MeshTopology%Num_Verts !!! Change that!
         Allocate(Bound_Ptr(1))
         Bound_Ptr = 0.0_Kr
         Do iDof = 1, NumDoF
            Call SectionIntRestrict(AppCtx%IrrevFlag, AppCtx%MeshTopology%Num_Elems+iDoF-1, IrrevFlag_Ptr, iErr); CHKERRQ(iErr)
            If (IrrevFlag_Ptr(1) /= VarFrac_BC_Type_NONE) Then
               Call SectionRealUpdate(LowerBoundV_Sec, AppCtx%MeshTopology%Num_Elems+iDoF-1, Bound_Ptr, INSERT_VALUES, iErr); CHKERRQ(iErr)
               Call SectionRealUpdate(UpperBoundV_Sec, AppCtx%MeshTopology%Num_Elems+iDoF-1, Bound_Ptr, INSERT_VALUES, iErr); CHKERRQ(iErr)
            End If
            Call SectionIntRestore(AppCtx%IrrevFlag, AppCtx%MeshTopology%Num_Elems+iDoF-1, IrrevFlag_Ptr, iErr); CHKERRQ(iErr)

            Call SectionIntRestrict(AppCtx%BCVFlag, AppCtx%MeshTopology%Num_Elems+iDoF-1, BCVFlag_Ptr, iErr); CHKERRQ(iErr)
            If (BCVFlag_Ptr(1) /= VarFrac_BC_Type_NONE) Then
               Call SectionRealRestrict(AppCtx%V,      AppCtx%MeshTopology%Num_Elems+iDoF-1, V_Ptr, iErr); CHKERRQ(iErr)
               Call SectionRealUpdate(LowerBoundV_Sec, AppCtx%MeshTopology%Num_Elems+iDoF-1, V_Ptr, INSERT_VALUES, iErr); CHKERRQ(iErr)
               Call SectionRealUpdate(UpperBoundV_Sec, AppCtx%MeshTopology%Num_Elems+iDoF-1, V_Ptr, INSERT_VALUES, iErr); CHKERRQ(iErr)
               Call SectionRealRestore(AppCtx%V,       AppCtx%MeshTopology%Num_Elems+iDoF-1, V_Ptr, iErr); CHKERRQ(iErr)
            End If
            Call SectionIntRestore(AppCtx%BCVFlag, AppCtx%MeshTopology%Num_Elems+iDoF-1, BCVFlag_Ptr, iErr); CHKERRQ(iErr)
         End Do  
         DeAllocate(Bound_Ptr)   
         
      Case (VarFrac_Irrev_InEq)      
         !!! Set Set lower bound to 0 and Upper bound to V, THEN to 0 at each point where at each point where BCFlag == VarFrac_BC_Type_DIRI
         SETERRQ(PETSC_ERR_SUP, 'Not Implemented yet\n', iErr)
      End Select
      
      Call SectionRealToVec(LowerBoundV_Sec, AppCtx%ScatterScal, SCATTER_FORWARD, LowerBoundV_Vec, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(UpperBoundV_Sec, AppCtx%ScatterScal, SCATTER_FORWARD, UpperBoundV_Vec, iErr); CHKERRQ(iErr)

      Call SectionRealDestroy(LowerBoundV_Sec, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(UpperBoundV_Sec, iErr); CHKERRQ(iErr)
      CHKMEMQ
   End Subroutine InitTaoBoundsV
#endif

   Subroutine Init_TS_V(AppCtx)
   !!! Split into Init_V and Init_IrrevBC_V
      Type(AppCtx_Type)                            :: AppCtx
      Type(SectionReal)                            :: VBC, VBT
      PetscReal, Dimension(:), Pointer             :: V_Ptr, VBC_Ptr
      PetscInt, Dimension(:), Pointer              :: BCVFlag, IrrevFlag
      PetscInt                                     :: i, iErr
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer      
      
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'VBC', 1, VBC, iErr); CHKERRQ(iErr)
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Fracture)%Offset, AppCtx%TimeStep, VBC)
            
      Select Case(AppCtx%VarFracSchemeParam%InitV)
      Case(VarFrac_INIT_V_PREV)
         If (AppCtx%AppParam%verbose > 0) Then
            Write(IOBuffer, *) "Initializing V with ", AppCtx%VarFracSchemeParam%InitV, "\n"
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         End If      
         Allocate(V_Ptr(1))
         Do i = 1, AppCtx%MeshTopology%Num_Verts  !!! NO!      
            !!! Take care of potential BC on V
            Call SectionIntRestrict(AppCtx%BCVFlag, AppCtx%MeshTopology%Num_Elems + i-1, BCVFlag, iErr); CHKERRQ(ierr)
            
            If (AppCtx%VarFracSchemeParam%IrrevType == VarFrac_Irrev_Eq ) Then
               !!! Take care of Irreversibility 
               Call SectionIntRestrict(AppCtx%IrrevFlag, AppCtx%MeshTopology%Num_Elems + i-1, IrrevFlag, iErr); CHKERRQ(ierr)
               If (IrrevFlag(1) /= VarFrac_BC_TYPE_NONE) Then
                  V_Ptr = 0.0_Kr
                  Call SectionRealUpdate(AppCtx%V, AppCtx%MeshTopology%Num_Elems + i-1, V_Ptr, INSERT_VALUES, iErr); CHKERRQ(iErr)
               End If
               Call SectionIntRestore(AppCtx%IrrevFlag, AppCtx%MeshTopology%Num_Elems + i-1, IrrevFlag, iErr); CHKERRQ(ierr)
            End If
            Call SectionIntRestore(AppCtx%BCVFlag, AppCtx%MeshTopology%Num_Elems + i-1, BCVFlag, iErr); CHKERRQ(ierr)
         End Do      
         DeAllocate(V_Ptr)
      Case default   
         SETERRQ(PETSC_ERR_SUP, 'Not Implemented yet\n', iErr)
      End Select
      
      Call SectionRealDestroy(VBC, iErr); CHKERRQ(iErr)
      CHKMEMQ
   End Subroutine Init_TS_V

   Subroutine Init_TS_Irrev(AppCtx)
   !!! Split into Init_V and Init_IrrevBC_V
      Type(AppCtx_Type)                            :: AppCtx
      Type(SectionReal)                            :: VBT
      PetscReal, Dimension(:), Pointer             :: V_Ptr
      PetscInt, Dimension(:), Pointer              :: BCVFlag, IrrevFlag
      PetscInt                                     :: i, iErr
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer      
      PetscReal                                    :: MyIrrevEQ_Counter 
      PetscReal                                    :: IrrevEQ_Counter
      
      MyIrrevEq_Counter = 0.0_Kr
      IrrevEq_Counter   = 0.0_Kr

      !!! Update the BC Flag from the current V
      Select Case(AppCtx%VarFracSchemeParam%IrrevType)
      Case(VarFrac_Irrev_Eq, VarFrac_Irrev_InEq)
         If (AppCtx%AppParam%verbose > 0) Then
            Write(IOBuffer, *) "Irreversibility with Irrev_EQ\n"
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         End If
         !!! Update the IrrevFlag Section
         Allocate(IrrevFlag(1))
         IrrevFlag = VarFrac_BC_Type_DIRI
         If (AppCtx%IsBT) Then
            If (AppCtx%AppParam%verbose > 0) Then
               Write(IOBuffer, *) "Backtracking, so reading timestep if TS>1\n"
               Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
            End If
            Call SectionIntZero(AppCtx%IrrevFlag, iErr); CHKERRQ(iErr)
            If (AppCtx%TimeStep > 1) Then
               Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'VBT', 1, VBT, iErr); CHKERRQ(iErr)      
               Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Fracture)%Offset, AppCtx%TimeStep-1, VBT)
               Do i = 1, AppCtx%MeshTopology%Num_Verts
                  Call SectionRealRestrict(VBT, AppCtx%MeshTopology%Num_Elems + i-1, V_Ptr, iErr); CHKERRQ(iErr)      
                  If (V_Ptr(1) < AppCtx%VarFracSchemeParam%IrrevTol) Then
                     Call SectionIntUpdate(AppCtx%IrrevFlag, AppCtx%MeshTopology%Num_Elems + i-1, IrrevFlag, INSERT_VALUES, iErr); CHKERRQ(iErr)
                     MyIrrevEQ_Counter = MyIrrevEQ_Counter + 1.0
                  End If
                  Call SectionRealRestore(VBT, AppCtx%MeshTopology%Num_Elems + i-1, V_Ptr, iErr); CHKERRQ(iErr)
               End Do
               Call SectionRealDestroy(VBT, iErr); CHKERRQ(iErr)
            End If
         Else
            If (AppCtx%AppParam%verbose > 0) Then
               Write(IOBuffer, *) "Not Backtracking, so getting blocked nodes from V\n"
               Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
            End If
            Do i = 1, AppCtx%MeshTopology%Num_Verts
               Call SectionRealRestrict(AppCtx%V, AppCtx%MeshTopology%Num_Elems + i-1, V_Ptr, iErr); CHKERRQ(ierr)      
               If (V_Ptr(1) < AppCtx%VarFracSchemeParam%IrrevTol) Then
                  Call SectionIntUpdate(AppCtx%IrrevFlag, AppCtx%MeshTopology%Num_Elems + i-1, IrrevFlag, INSERT_VALUES, iErr); CHKERRQ(iErr)
                  MyIrrevEQ_Counter = MyIrrevEQ_Counter + 1.0
               End If
               Call SectionRealRestore(AppCtx%V, AppCtx%MeshTopology%Num_Elems + i-1, V_Ptr, iErr);             End Do
         End If
         DeAllocate(IrrevFlag)
         If (AppCtx%AppParam%verbose > 0) Then
            Call PetscGlobalSum(MyIrrevEQ_Counter, IrrevEQ_Counter, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
            Write(IOBuffer, *) "Number of blocked nodes for V: ", IrrevEQ_Counter, "\n"
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         End If      
      Case default   
         SETERRQ(PETSC_ERR_SUP, 'Not Implemented yet\n', iErr)
      End Select
      CHKMEMQ
   End Subroutine Init_TS_Irrev

!!!
!!! Global assembly functions
!!!
   
   Subroutine MatV_Assembly(K, AppCtx)
      Type(Mat)                                    :: K
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iBlk, iBlkID, iErr
      
      Call PetscLogStagePush(AppCtx%LogInfo%MatAssemblyV_Stage, iErr); CHKERRQ(iErr)

      Call MatZeroEntries(K, iErr); CHKERRQ(iErr)
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
         Select Case (AppCtx%VarFracSchemeParam%AtNum)
         Case(1)
            SETERRQ(PETSC_ERR_SUP, 'AT1 requires V_tao', iErr)
         Case(2)
            If (AppCtx%MyEXO%EBProperty(VarFrac_EBProp_IsBrittle)%Value(iBlkID) /= 0) Then
               If (AppCtx%VarFracSchemeParam%Unilateral /= 0) Then
                  Call MatV_AssemblyBlk_ElastBrittleUnilateral(K, iBlk, .TRUE., AppCtx)
               Else
                  Call MatV_AssemblyBlk_ElastBrittle(K, iBlk, .TRUE., AppCtx)
               End If
            End If
            Call MatV_AssemblyBlk_SurfaceAT2(K, iBlk, .TRUE., AppCtx)
         Case Default
            SETERRQ(PETSC_ERR_SUP, 'Only AT1 and AT2 are implemented\n', iErr)
      End Select
      End Do Do_Elem_iBlk
      Call MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      If (AppCtx%AppParam%verbose > 1) Then
         Call MatView(K, PETSC_VIEWER_STDOUT_WORLD, iErr)
      End If

      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
      CHKMEMQ
   End Subroutine MatV_Assembly

#if defined WITH_TAO
   Subroutine HessianV_Assembly(tao, X, H, Hpre, flg, AppCtx, iErr)
      TAO_SOLVER         :: tao
      Type(Vec)          :: X
      Type(Mat)          :: H, Hpre
      PetscInt           :: iErr
      MatStructure       :: flg
      Type(AppCtx_Type)  :: AppCtx
      
      PetscInt           :: iBlk, iBlkID
      
      Call PetscLogStagePush(AppCtx%LogInfo%MatAssemblyV_Stage, iErr); CHKERRQ(iErr)

      Call MatZeroEntries(H, iErr); CHKERRQ(iErr)
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
         Select Case (AppCtx%VarFracSchemeParam%AtNum)
         Case(1)
            If (AppCtx%MyEXO%EBProperty(VarFrac_EBProp_IsBrittle)%Value(iBlkID) /= 0) Then
               If (AppCtx%VarFracSchemeParam%Unilateral /= 0) Then
                  Call MatV_AssemblyBlk_ElastBrittleUnilateral(H, iBlk, .FALSE., AppCtx)
               Else
                  Call MatV_AssemblyBlk_ElastBrittle(H, iBlk, .FALSE., AppCtx)
               End If
            End If
            Call MatV_AssemblyBlk_SurfaceAT1(H, iBlk, .FALSE., AppCtx)
         Case(2)
            If (AppCtx%MyEXO%EBProperty(VarFrac_EBProp_IsBrittle)%Value(iBlkID) /= 0) Then
               Call MatV_AssemblyBlk_ElastBrittle(H, iBlk, .FALSE., AppCtx)
            End If
            Call MatV_AssemblyBlk_SurfaceAT2(H, iBlk, .FALSE., AppCtx)
         Case Default
            SETERRQ(PETSC_ERR_SUP, 'Only AT1 and AT2 are implemented\n', iErr)
         End Select
      End Do Do_Elem_iBlk
      Call MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (H, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      If (AppCtx%AppParam%verbose > 1) Then
         Call MatView(H, PETSC_VIEWER_STDOUT_WORLD, iErr)
      End If

      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
      CHKMEMQ
   End Subroutine HessianV_Assembly
#endif

#if defined WITH_TAO
   Subroutine FormFunctionAndGradientV(tao, V_Vec, ObjFunc, GradientV_Vec, AppCtx, iErr)
      TAO_SOLVER                                   :: tao
      Type(Vec)                                    :: V_Vec, GradientV_Vec
      PetscInt                                     :: iErr
      PetscReal                                    :: ObjFunc
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iBlk, iBlkId
      Type(SectionReal)                            :: GradientV_Sec
      PetscReal                                    :: MyElasticEnergyBlock, MySurfaceEnergyBlock
      PetscReal                                    :: MyObjFunc
      
      !!! Objective function is ElasticEnergy + SurfaceEnergy
      Call SectionRealToVec(AppCtx%V, AppCtx%ScatterScal, SCATTER_REVERSE, V_Vec, iErr); CHKERRQ(ierr)

      MyObjFunc = 0.0_Kr
      Do_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
         If (AppCtx%MyEXO%EBProperty(VarFrac_EBProp_IsBrittle)%Value(iBlkID) /= 0) Then
            If (AppCtx%VarFracSchemeParam%Unilateral /= 0) Then
               Call ElasticEnergy_AssemblyBlk_BrittleUnilateral(MyElasticEnergyBlock, iBlk, AppCtx%U, AppCtx%Theta, AppCtx%V, AppCtx)
            Else
               Call ElasticEnergy_AssemblyBlk_Brittle(MyElasticEnergyBlock, iBlk, AppCtx%U, AppCtx%Theta, AppCtx%V, AppCtx)
            End If
            MyObjFunc = MyObjFunc + MyElasticEnergyBlock
         End If

         Select Case (AppCtx%VarFracSchemeParam%AtNum)
         Case(1)
            Call SurfaceEnergy_AssemblyBlk_AT1(MySurfaceEnergyBlock, iBlk, AppCtx%V, AppCtx)
         Case(2)
            Call SurfaceEnergy_AssemblyBlk_AT2(MySurfaceEnergyBlock, iBlk, AppCtx%V, AppCtx)
         Case Default
            SETERRQ(PETSC_ERR_SUP, 'Only AT1 and AT2 are implemented\n', iErr)
         End Select
         MyObjFunc = MyObjFunc + MySurfaceEnergyBlock
      End Do Do_iBlk
      Call PetscGlobalSum(MyObjFunc, ObjFunc, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      
      !!! Gradient
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'GradientV', 1, GradientV_Sec, iErr); CHKERRQ(iErr)
      Call SectionRealZero(GradientV_Sec, iErr); CHKERRQ(iErr)

      Do_iBlk2: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
         If (AppCtx%MyEXO%EBProperty(VarFrac_EBProp_IsBrittle)%Value(iBlkID) /= 0) Then
            !!! Contribution of the bulk term 1/2 \int v^2 Ae(u):e(u)
            If (AppCtx%VarFracSchemeParam%Unilateral /= 0) Then
               Call GradientV_AssemblyBlk_ElastBrittleUnilateral(GradientV_Sec, iBlk, AppCtx%U, AppCtx%Theta, AppCtx%V, AppCtx)
            Else
               Call GradientV_AssemblyBlk_ElastBrittle(GradientV_Sec, iBlk, AppCtx%U, AppCtx%Theta, AppCtx%V, AppCtx)
            End If
            
         End If
         !!! Contribution of the surface term
         Select Case (AppCtx%VarFracSchemeParam%AtNum)
         Case(1)
            Call GradientV_AssemblyBlk_SurfaceAT1(GradientV_Sec, iBlk, AppCtx%V, AppCtx)
         Case(2)
            Call GradientV_AssemblyBlk_SurfaceAT2(GradientV_Sec, iBlk, AppCtx%V, AppCtx)
         Case Default
            SETERRQ(PETSC_ERR_SUP, 'Only AT1 and AT2 are implemented\n', iErr)
         End Select
      End Do Do_iBlk2
      Call SectionRealComplete(GradientV_Sec, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(GradientV_Sec, AppCtx%ScatterScal, SCATTER_FORWARD, GradientV_Vec, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(GradientV_Sec, iErr); CHKERRQ(iErr)
      CHKMEMQ
   End Subroutine FormFunctionAndGradientV
#endif

   Subroutine RHSV_Assembly(RHSV_Vec, AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      Type(Vec)                                    :: RHSV_Vec

      Type(SectionReal)                            :: RHSV_Sec
      
      PetscInt                                     :: iErr
      PetscInt                                     :: iBlk
            
      Call PetscLogStagePush(AppCtx%LogInfo%RHSAssemblyV_Stage, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'RHSV_Sec', 1, RHSV_Sec, iErr); CHKERRQ(iErr)
      
      Call SectionRealZero(RHSV_Sec, iErr); CHKERRQ(iErr)
      Do_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Select Case (AppCtx%VarFracSchemeParam%AtNum)
         Case(2)
            Call RHSV_AssemblyBlk_AT2(RHSV_Sec, iBlk, AppCtx)
         Case Default
            SETERRQ(PETSC_ERR_SUP, 'Only AT2 is implemented with KSP\n', iErr)
      End Select
      End Do Do_iBlk

      Call SectionRealComplete(RHSV_Sec, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(RHSV_Sec, AppCtx%ScatterScal, SCATTER_FORWARD, RHSV_Vec, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(RHSV_Sec, iErr); CHKERRQ(iErr)
      
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
      CHKMEMQ
   End Subroutine RHSV_Assembly
   
   Subroutine MatV_AssemblyBlk_ElastBrittle(H, iBlk, DoBC, AppCtx)
      Type(Mat)                                    :: H
      PetscInt                                     :: iBlk
      PetscTruth                                   :: DoBC
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iBlkID
      PetscReal, Dimension(:,:), Pointer           :: MatElem      
      PetscInt                                     :: iELoc, iE
   
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoFScal, NumDoFVect
      PetscInt, Dimension(:), Pointer              :: BCFlag, IrrevFlag
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      
      PetscReal, Dimension(:), Pointer             :: U, Theta
      PetscReal                                    :: Theta_Elem, ElasEnDens_Elem
#if defined PB_2D
      Type(MatS2D)                                 :: Strain_Elem, EffectiveStrain_Elem
#elif defined PB_3D
      Type(MatS3D)                                 :: Strain_Elem, EffectiveStrain_Elem
#endif      
      PetscReal                                    :: ElasticEnergyDensity
      PetscLogDouble                               :: flops
      
      Call PetscLogEventBegin(AppCtx%LogInfo%MatAssemblyLocalV_Event, iErr); CHKERRQ(iErr)

      flops = 0.0
      NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
      NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID

      Allocate(U(NumDoFVect))
      Allocate(Theta(NumDoFScal))
      Allocate(BCFlag(NumDoFScal))
      BCFlag = VarFrac_BC_Type_NONE
      Allocate(IrrevFlag(NumDoFScal))
      IrrevFlag = VarFrac_BC_Type_NONE
      Allocate(MatElem(NumDoFScal, NumDoFScal))
      
      Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         MatElem  = 0.0_Kr
         !! Get the local nodal values of U, Theta, and BCs
         Call SectionRealRestrictClosure(AppCtx%U, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U, iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(AppCtx%Theta, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, Theta, iErr); CHKERRQ(ierr)
         If (DoBC) Then
            Call SectionIntRestrictClosure(AppCtx%BCVFlag, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, BCFlag, iErr); CHKERRQ(ierr)
            Call SectionIntRestrictClosure(AppCtx%IrrevFlag, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, IrrevFlag, iErr); CHKERRQ(ierr)
         End If
      
         Do iGauss = 1, Size(AppCtx%ElemScal(iE)%Gauss_C)
            !! Calculate the strain at the gauss point
            Strain_Elem = 0.0_Kr
            Theta_Elem  = 0.0_Kr
            Do iDoF1 = 1, NumDoFVect
               Strain_Elem = Strain_Elem + (AppCtx%ElemVect(iE)%GradS_BF(iDoF1, iGauss) * U(iDoF1))
            End Do
            !! Calculate the temperature at the gauss point
            Do iDoF1 = 1, NumDoFScal
               Theta_Elem = Theta_Elem + AppCtx%ElemScal(iE)%BF(iDoF1, iGauss) * Theta(iDoF1)
               flops = flops + 2.0
            End Do
            !! Calculate the Effective Strain at the gauss point
            EffectiveStrain_Elem  = Strain_Elem - (Theta_Elem * AppCtx%MatProp(iBlkID)%Therm_Exp)   
            ElasticEnergyDensity  = (AppCtx%MatProp(iBlkID)%Hookes_Law * EffectiveStrain_Elem) .DotP. EffectiveStrain_Elem
            !! Assemble the element stiffness
            Do iDoF1 = 1, NumDoFScal
               If ( (BCFlag(iDoF1) == VarFrac_BC_Type_NONE) .AND. (IrrevFlag(iDoF1) == VarFrac_BC_Type_NONE) ) Then
                  Do iDoF2 = 1, NumDoFScal
                  !! Terms in V^2
                     MatElem(iDoF2, iDoF1) =  MatElem(iDoF2, iDoF1) + AppCtx%ElemScal(iE)%Gauss_C(iGauss) * ElasticEnergyDensity * AppCtx%ElemScal(iE)%BF(iDoF1, iGauss) * AppCtx%ElemScal(iE)%BF(iDoF2, iGauss)
                     flops = flops + 4.0
                  End Do
               End If
            End Do
         End Do
         Call assembleMatrix(H, AppCtx%MeshTopology%mesh, AppCtx%V, iE-1, MatElem, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_Elem_iE
      Call PetscLogFlops(flops, iErr); CHKERRQ(iErr)
      DeAllocate(MatElem)
      DeAllocate(U)      
      DeAllocate(Theta)
      DeAllocate(BCFlag)
      DeAllocate(IrrevFlag)
      Call PetscLogEventEnd(AppCtx%LogInfo%MatAssemblyLocalV_Event, iErr); CHKERRQ(iErr)
      CHKMEMQ
   End Subroutine MatV_AssemblyBlk_ElastBrittle

   Subroutine MatV_AssemblyBlk_ElastBrittleUnilateral(H, iBlk, DoBC, AppCtx)
      Type(Mat)                                    :: H
      PetscInt                                     :: iBlk
      PetscTruth                                   :: DoBC
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iBlkID
      PetscReal, Dimension(:,:), Pointer           :: MatElem      
      PetscInt                                     :: iELoc, iE
   
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoFScal, NumDoFVect
      PetscInt, Dimension(:), Pointer              :: BCFlag, IrrevFlag
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      
      PetscReal, Dimension(:), Pointer             :: U, Theta
      PetscReal                                    :: Theta_Elem
#if defined PB_2D
      Type(MatS2D)                                 :: Strain_Elem, EffectiveStrain_Elem
      Type(MatS2D)                                 :: EffectiveStrain_Elem_D
#elif defined PB_3D
      Type(MatS3D)                                 :: Strain_Elem, EffectiveStrain_Elem
      Type(MatS3D)                                 :: EffectiveStrain_Elem_D
#endif      
      PetscReal                                    :: EffectiveStrain_Trace
      PetscReal                                    :: ElasticEnergyDensity
      PetscLogDouble                               :: flops
      
      Call PetscLogEventBegin(AppCtx%LogInfo%MatAssemblyLocalV_Event, iErr); CHKERRQ(iErr)

      flops = 0.0
      NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
      NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID

      Allocate(U(NumDoFVect))
      Allocate(Theta(NumDoFScal))
      Allocate(BCFlag(NumDoFScal))
      BCFlag = VarFrac_BC_Type_NONE
      Allocate(IrrevFlag(NumDoFScal))
      IrrevFlag = VarFrac_BC_Type_NONE
      Allocate(MatElem(NumDoFScal, NumDoFScal))
      
      Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         MatElem  = 0.0_Kr
         !! Get the local nodal values of U, Theta, and BCs
         Call SectionRealRestrictClosure(AppCtx%U, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U, iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(AppCtx%Theta, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, Theta, iErr); CHKERRQ(ierr)
         If (DoBC) Then
            Call SectionIntRestrictClosure(AppCtx%BCVFlag, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, BCFlag, iErr); CHKERRQ(ierr)
            Call SectionIntRestrictClosure(AppCtx%IrrevFlag, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, IrrevFlag, iErr); CHKERRQ(ierr)
         End If
      
         Do iGauss = 1, Size(AppCtx%ElemScal(iE)%Gauss_C)
            !! Calculate the strain at the gauss point
            Strain_Elem = 0.0_Kr
            Theta_Elem  = 0.0_Kr
            Do iDoF1 = 1, NumDoFVect
               Strain_Elem = Strain_Elem + (AppCtx%ElemVect(iE)%GradS_BF(iDoF1, iGauss) * U(iDoF1))
            End Do
            !! Calculate the temperature at the gauss point
            Do iDoF1 = 1, NumDoFScal
               Theta_Elem = Theta_Elem + AppCtx%ElemScal(iE)%BF(iDoF1, iGauss) * Theta(iDoF1)
               flops = flops + 2.0
            End Do
            !! Calculate the Effective Strain at the gauss point
            EffectiveStrain_Elem  = Strain_Elem - (Theta_Elem * AppCtx%MatProp(iBlkID)%Therm_Exp)   
            EffectiveStrain_Trace = Trace(EffectiveStrain_Elem)

            If (EffectiveStrain_Trace >= 0.0_Kr) Then
               ElasticEnergyDensity = (AppCtx%MatProp(iBlkID)%Hookes_Law * EffectiveStrain_Elem) .DotP. EffectiveStrain_Elem
            Else
               EffectiveStrain_Elem_D = DeviatoricPart(EffectiveStrain_Elem)
               ElasticEnergyDensity = (AppCtx%MatProp(iBlkID)%Hookes_Law * EffectiveStrain_Elem_D) .DotP. EffectiveStrain_Elem_D
            End If

            !! Assemble the element stiffness
            Do iDoF1 = 1, NumDoFScal
               If ( (BCFlag(iDoF1) == VarFrac_BC_Type_NONE) .AND. (IrrevFlag(iDoF1) == VarFrac_BC_Type_NONE) ) Then
                  Do iDoF2 = 1, NumDoFScal
                  !! Terms in V^2
                     MatElem(iDoF2, iDoF1) =  MatElem(iDoF2, iDoF1) + AppCtx%ElemScal(iE)%Gauss_C(iGauss) * ElasticEnergyDensity * AppCtx%ElemScal(iE)%BF(iDoF1, iGauss) * AppCtx%ElemScal(iE)%BF(iDoF2, iGauss)
                     flops = flops + 4.0
                  End Do
               End If
            End Do
         End Do
         Call assembleMatrix(H, AppCtx%MeshTopology%mesh, AppCtx%V, iE-1, MatElem, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_Elem_iE
      Call PetscLogFlops(flops, iErr); CHKERRQ(iErr)
      DeAllocate(MatElem)
      DeAllocate(U)      
      DeAllocate(Theta)
      DeAllocate(BCFlag)
      DeAllocate(IrrevFlag)
      Call PetscLogEventEnd(AppCtx%LogInfo%MatAssemblyLocalV_Event, iErr); CHKERRQ(iErr)
      CHKMEMQ
   End Subroutine MatV_AssemblyBlk_ElastBrittleUnilateral

   Subroutine MatV_AssemblyBlk_SurfaceAT2(H, iBlk, DoBC, AppCtx)
      Type(Mat)                                    :: H
      PetscInt                                     :: iBlk
      PetscTruth                                   :: DoBC
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iBlkID
      PetscReal, Dimension(:,:), Pointer           :: MatElem      
      PetscInt                                     :: iELoc, iE
   
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoFScal, NumDoFVect
      PetscInt, Dimension(:), Pointer              :: BCFlag, IrrevFlag
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      
      PetscReal                                    :: C2_V, C2_GradV
      PetscLogDouble                               :: flops
      
      Call PetscLogEventBegin(AppCtx%LogInfo%MatAssemblyLocalV_Event, iErr); CHKERRQ(iErr)

      flops = 0.0
      NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
      NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID

      Allocate(BCFlag(NumDoFScal))
      BCFlag = VarFrac_BC_Type_NONE
      Allocate(IrrevFlag(NumDoFScal))
      IrrevFlag = VarFrac_BC_Type_NONE
      Allocate(MatElem(NumDoFScal, NumDoFScal))
      
      Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         MatElem  = 0.0_Kr
         !! Get the local nodal values of BC sections if necessary
         If (DoBC) Then
            Call SectionIntRestrictClosure(AppCtx%BCVFlag, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, BCFlag, iErr); CHKERRQ(ierr)
            Call SectionIntRestrictClosure(AppCtx%IrrevFlag, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, IrrevFlag, iErr); CHKERRQ(ierr)
         End If
      
         Do iGauss = 1, Size(AppCtx%ElemScal(iE)%Gauss_C)
            C2_V = AppCtx%MatProp(iBlkID)%Toughness / AppCtx%VarFracSchemeParam%Epsilon / AppCtx%VarFracSchemeParam%ATCv * 0.5_Kr
            flops = flops + 2.0
            C2_GradV = AppCtx%MatProp(iBlkID)%Toughness * AppCtx%VarFracSchemeParam%Epsilon / AppCtx%VarFracSchemeParam%ATCv * 0.5_Kr
            flops = flops + 2.0
            !! Assemble the element stiffness
            Do iDoF1 = 1, NumDoFScal
               If ( (BCFlag(iDoF1) == VarFrac_BC_Type_NONE) .AND. (IrrevFlag(iDoF1) == VarFrac_BC_Type_NONE) ) Then
                  Do iDoF2 = 1, NumDoFScal
                  !! Terms in V^2
                     MatElem(iDoF2, iDoF1) =  MatElem(iDoF2, iDoF1) + AppCtx%ElemScal(iE)%Gauss_C(iGauss) * C2_V * AppCtx%ElemScal(iE)%BF(iDoF1, iGauss) * AppCtx%ElemScal(iE)%BF(iDoF2, iGauss)
                     flops = flops + 4.0
                  !! Terms in GradV^2
                     MatElem(iDoF2, iDoF1) =  MatElem(iDoF2, iDoF1) + AppCtx%ElemScal(iE)%Gauss_C(iGauss) * C2_GradV * (AppCtx%ElemScal(iE)%Grad_BF(iDoF1, iGauss) .DotP. AppCtx%ElemScal(iE)%Grad_BF(iDoF2, iGauss) )               
                     flops = flops + 3.0
                  End Do
               End If
            End Do
         End Do
         Call assembleMatrix(H, AppCtx%MeshTopology%mesh, AppCtx%V, iE-1, MatElem, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_Elem_iE
      Call PetscLogFlops(flops, iErr); CHKERRQ(iErr)
      DeAllocate(MatElem)
      DeAllocate(BCFlag)
      DeAllocate(IrrevFlag)
      Call PetscLogEventEnd(AppCtx%LogInfo%MatAssemblyLocalV_Event, iErr); CHKERRQ(iErr)
      CHKMEMQ
   End Subroutine MatV_AssemblyBlk_SurfaceAT2

   Subroutine MatV_AssemblyBlk_SurfaceAT1(H, iBlk, DoBC, AppCtx)
      Type(Mat)                                    :: H
      PetscInt                                     :: iBlk
      PetscTruth                                   :: DoBC
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iBlkID
      PetscReal, Dimension(:,:), Pointer           :: MatElem      
      PetscInt                                     :: iELoc, iE
   
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoFScal, NumDoFVect
      PetscInt, Dimension(:), Pointer              :: BCFlag, IrrevFlag
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      
      PetscReal                                    :: C2_GradV
      PetscLogDouble                               :: flops
      
      Call PetscLogEventBegin(AppCtx%LogInfo%MatAssemblyLocalV_Event, iErr); CHKERRQ(iErr)

      flops = 0.0
      NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
      NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID

      Allocate(BCFlag(NumDoFScal))
      BCFlag = VarFrac_BC_Type_NONE
      Allocate(IrrevFlag(NumDoFScal))
      IrrevFlag = VarFrac_BC_Type_NONE
      Allocate(MatElem(NumDoFScal, NumDoFScal))
      
      Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         MatElem  = 0.0_Kr
         !! Get the local nodal values of BC sections if necessary
         If (DoBC) Then
            Call SectionIntRestrictClosure(AppCtx%BCVFlag, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, BCFlag, iErr); CHKERRQ(ierr)
            Call SectionIntRestrictClosure(AppCtx%IrrevFlag, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, IrrevFlag, iErr); CHKERRQ(ierr)
         End If
      
         Do iGauss = 1, Size(AppCtx%ElemScal(iE)%Gauss_C)
            C2_GradV = AppCtx%MatProp(iBlkID)%Toughness * AppCtx%VarFracSchemeParam%Epsilon / AppCtx%VarFracSchemeParam%ATCv * 0.5_Kr
            flops = flops + 2.0
            !! Assemble the element stiffness
            Do iDoF1 = 1, NumDoFScal
               If ( (BCFlag(iDoF1) == VarFrac_BC_Type_NONE) .AND. (IrrevFlag(iDoF1) == VarFrac_BC_Type_NONE) ) Then
                  Do iDoF2 = 1, NumDoFScal
                  !! Terms in GradV^2
                     MatElem(iDoF2, iDoF1) =  MatElem(iDoF2, iDoF1) + AppCtx%ElemScal(iE)%Gauss_C(iGauss) * C2_GradV * (AppCtx%ElemScal(iE)%Grad_BF(iDoF1, iGauss) .DotP. AppCtx%ElemScal(iE)%Grad_BF(iDoF2, iGauss) )               
                     flops = flops + 3.0
                  End Do
               End If
            End Do
         End Do
         Call assembleMatrix(H, AppCtx%MeshTopology%mesh, AppCtx%V, iE-1, MatElem, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_Elem_iE
      Call PetscLogFlops(flops, iErr); CHKERRQ(iErr)
      DeAllocate(MatElem)
      DeAllocate(BCFlag)
      DeAllocate(IrrevFlag)
      Call PetscLogEventEnd(AppCtx%LogInfo%MatAssemblyLocalV_Event, iErr); CHKERRQ(iErr)
      CHKMEMQ
   End Subroutine MatV_AssemblyBlk_SurfaceAT1

   Subroutine RHSV_AssemblyBlk_AT2(RHSV_Sec, iBlk, AppCtx)
      Type(SectionReal)                            :: RHSV_Sec
      PetscInt                                     :: iBlk
      Type(AppCtx_Type)                            :: AppCtx

      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_ELem = \sum_i U_Loc(i) BF(i))
      PetscReal, Dimension(:), Pointer             :: RHS_Loc
      PetscInt, Dimension(:), Pointer              :: BCFlag_Loc, IrrevFlag_Loc
      PetscInt                                     :: iE, iEloc, iBlkId, iErr

      PetscInt                                     :: NumDoFScal
      PetscInt                                     :: iDoF1, iGauss
      PetscReal                                    :: C1_V
      PetscLogDouble                               :: flops

      Call PetscLogEventBegin(AppCtx%LogInfo%RHSAssemblyLocalV_Event, iErr); CHKERRQ(iErr)
      flops      = 0.0
      
      NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
      Allocate(BCFlag_Loc(NumDoFScal))
      Allocate(IrrevFlag_Loc(NumDoFScal))      
      
      Allocate(RHS_Loc(NumDoFScal))

      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
      C1_V =  AppCtx%MatProp(iBlkId)%Toughness / AppCtx%VarFracSchemeParam%Epsilon / AppCtx%VarFracSchemeParam%ATCv * 0.5_Kr
      flops = flops + 2.0

      Do_iEloc: Do iEloc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         RHS_Loc = 0.0_Kr

         Call SectionIntRestrictClosure(AppCtx%BCVFlag, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, BCFlag_Loc, iErr); CHKERRQ(ierr)
         Call SectionIntRestrictClosure(AppCtx%IrrevFlag, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, IrrevFlag_Loc, iErr); CHKERRQ(ierr)

         Do_iGauss: Do iGauss = 1, Size(AppCtx%ElemScal(iE)%Gauss_C)
            Do iDoF1 = 1, NumDoFScal
               If  ( (BCFlag_Loc(iDoF1) == VarFrac_BC_Type_NONE) .AND. (IrrevFlag_Loc(iDoF1) == VarFrac_BC_Type_NONE) ) Then
               ! RHS terms due to forces
                  RHS_Loc(iDoF1) = RHS_Loc(iDoF1) + C1_V * AppCtx%ElemScal(iE)%Gauss_C(iGauss) *  AppCtx%ElemScal(iE)%BF(iDoF1, iGauss) 
                  flops = flops + 3.0
               End If
            End Do
         End Do Do_iGauss
         Call SectionRealUpdateClosure(RHSV_Sec, AppCtx%MeshTopology%Mesh, iE-1, RHS_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iEloc

      DeAllocate(BCFlag_Loc)
      DeAllocate(IrrevFlag_Loc)      
      DeAllocate(RHS_Loc)
      Call PetscLogFlops(flops, iErr); CHKERRQ(iErr)
      Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyLocalV_Event, iErr); CHKERRQ(iErr)
      CHKMEMQ
   End Subroutine RHSV_AssemblyBlk_AT2
   
#if defined WITH_TAO   
   Subroutine GradientV_AssemblyBlk_ElastBrittle(GradientV_Sec, iBlk, U_Sec, Theta_Sec, V_Sec, AppCtx)
      Type(SectionReal)                            :: GradientV_Sec
      PetscInt                                     :: iBlk
      Type(SectionReal)                            :: U_Sec, Theta_SeC, V_Sec
      Type(AppCtx_Type)                            :: AppCtx

      PetscReal, Dimension(:), Pointer             :: U_Loc, Theta_Loc, V_Loc, GradientV_Loc
      PetscReal                                    :: Theta_Elem, V_Elem
#if defined PB_2D
      Type(MatS2D)                                 :: Strain_Elem, EffectiveStrain_Elem
#elif defined PB_3D
      Type(MatS3D)                                 :: Strain_Elem, EffectiveStrain_Elem
#endif      
      PetscReal                                    :: ElasticEnergyDensity
      PetscInt                                     :: iE, iEloc, iBlkId, iErr
      PetscInt                                     :: NumDoFScal, NumDoFVect
      PetscInt                                     :: iDoF, iGauss
      PetscLogDouble                               :: flops

      flops        = 0.0

      NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
      NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF

      Allocate(U_Loc(NumDoFVect))
      Allocate(Theta_Loc(NumDoFScal))
      Allocate(GradientV_Loc(NumDoFScal))
      Allocate(V_Loc(NumDoFScal))

      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
      Do_iEloc: Do iEloc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         GradientV_Loc = 0.0_Kr
         Call SectionRealRestrictClosure(U_Sec,     AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U_Loc,     iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(Theta_Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, Theta_Loc, iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(V_Sec,     AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, V_Loc,     iErr); CHKERRQ(ierr)
         Do_iGauss: Do iGauss = 1, size(AppCtx%ElemVect(iE)%Gauss_C)
            Strain_Elem = 0.0_Kr
            Do iDoF = 1, NumDoFVect
               Strain_Elem = Strain_Elem + U_Loc(iDoF) * AppCtx%ElemVect(iE)%GradS_BF(iDoF, iGauss)
            End Do
            Theta_Elem = 0.0_Kr
            V_Elem     = 0.0_Kr
            Do iDoF = 1, NumDoFScal
               Theta_Elem = Theta_Elem + Theta_Loc(iDoF) * AppCtx%ElemScal(iE)%BF(iDoF, iGauss)
               V_Elem     = V_Elem     + V_Loc(iDoF)     * AppCtx%ElemScal(iE)%BF(iDoF, iGauss)
               flops = flops + 4.0
            End Do
            EffectiveStrain_Elem  = Strain_Elem - AppCtx%MatProp(iBlkId)%Therm_Exp * Theta_Elem
            ElasticEnergyDensity  = (AppCtx%MatProp(iBlkId)%Hookes_Law * EffectiveStrain_Elem) .DotP. EffectiveStrain_Elem
            Do iDoF = 1, NumDofScal
               GradientV_Loc(iDoF) = GradientV_Loc(iDoF) + AppCtx%ElemScal(iE)%Gauss_C(iGauss) * ElasticEnergyDensity * V_Elem * AppCtx%ElemScal(iE)%BF(iDoF, iGauss)
            End Do
         End Do Do_iGauss
         Call SectionRealUpdateClosure(GradientV_Sec, AppCtx%MeshTopology%Mesh, iE-1, GradientV_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iEloc

      DeAllocate(U_Loc)
      DeAllocate(Theta_Loc)
      DeAllocate(GradientV_Loc)
      DeAllocate(V_Loc)
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
      CHKMEMQ
   End Subroutine GradientV_AssemblyBlk_ElastBrittle
   
#if defined WITH_TAO   
   Subroutine GradientV_AssemblyBlk_ElastBrittleUnilateral(GradientV_Sec, iBlk, U_Sec, Theta_Sec, V_Sec, AppCtx)
      Type(SectionReal)                            :: GradientV_Sec
      PetscInt                                     :: iBlk
      Type(SectionReal)                            :: U_Sec, Theta_Sec, V_Sec
      Type(AppCtx_Type)                            :: AppCtx

      PetscReal, Dimension(:), Pointer             :: U_Loc, Theta_Loc, V_Loc, GradientV_Loc
      PetscReal                                    :: Theta_Elem, V_Elem
#if defined PB_2D
      Type(MatS2D)                                 :: Strain_Elem, EffectiveStrain_Elem
      Type(MatS2D)                                 :: EffectiveStrain_Elem_D
#elif defined PB_3D
      Type(MatS3D)                                 :: Strain_Elem, EffectiveStrain_Elem
      Type(MatS3D)                                 :: EffectiveStrain_Elem_D
#endif      
      PetscReal                                    :: EffectiveStrain_Trace
      PetscReal                                    :: ElasticEnergyDensity
      PetscInt                                     :: iE, iEloc, iBlkId, iErr
      PetscInt                                     :: NumDoFScal, NumDoFVect
      PetscInt                                     :: iDoF, iGauss
      PetscLogDouble                               :: flops

      flops        = 0.0

      NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
      NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF

      Allocate(U_Loc(NumDoFVect))
      Allocate(Theta_Loc(NumDoFScal))
      Allocate(GradientV_Loc(NumDoFScal))
      Allocate(V_Loc(NumDoFScal))

      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
      Do_iEloc: Do iEloc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         GradientV_Loc = 0.0_Kr
         Call SectionRealRestrictClosure(U_Sec,     AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U_Loc,     iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(Theta_Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, Theta_Loc, iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(V_Sec,     AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, V_Loc,     iErr); CHKERRQ(ierr)
         Do_iGauss: Do iGauss = 1, size(AppCtx%ElemVect(iE)%Gauss_C)
            Strain_Elem = 0.0_Kr
            Do iDoF = 1, NumDoFVect
               Strain_Elem = Strain_Elem + U_Loc(iDoF) * AppCtx%ElemVect(iE)%GradS_BF(iDoF, iGauss)
            End Do
            Theta_Elem = 0.0_Kr
            V_Elem     = 0.0_Kr
            Do iDoF = 1, NumDoFScal
               Theta_Elem = Theta_Elem + Theta_Loc(iDoF) * AppCtx%ElemScal(iE)%BF(iDoF, iGauss)
               V_Elem     = V_Elem     + V_Loc(iDoF)     * AppCtx%ElemScal(iE)%BF(iDoF, iGauss)
               flops = flops + 4.0
            End Do
            EffectiveStrain_Elem  = Strain_Elem - AppCtx%MatProp(iBlkId)%Therm_Exp * Theta_Elem
            EffectiveStrain_Trace = Trace(EffectiveStrain_Elem)
            
            If (EffectiveStrain_Trace >= 0.0_Kr) Then
               ElasticEnergyDensity = (AppCtx%MatProp(iBlkID)%Hookes_Law * EffectiveStrain_Elem) .DotP. EffectiveStrain_Elem
            Else
               EffectiveStrain_Elem_D = DeviatoricPart(EffectiveStrain_Elem)
               ElasticEnergyDensity = (AppCtx%MatProp(iBlkID)%Hookes_Law * EffectiveStrain_Elem_D) .DotP. EffectiveStrain_Elem_D
            End If
            Do iDoF = 1, NumDofScal
               GradientV_Loc(iDoF) = GradientV_Loc(iDoF) + AppCtx%ElemScal(iE)%Gauss_C(iGauss) * ElasticEnergyDensity * V_Elem * AppCtx%ElemScal(iE)%BF(iDoF, iGauss)
            End Do
         End Do Do_iGauss
         Call SectionRealUpdateClosure(GradientV_Sec, AppCtx%MeshTopology%Mesh, iE-1, GradientV_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iEloc

      DeAllocate(U_Loc)
      DeAllocate(Theta_Loc)
      DeAllocate(GradientV_Loc)
      DeAllocate(V_Loc)
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
      CHKMEMQ
   End Subroutine GradientV_AssemblyBlk_ElastBrittleUnilateral
#endif
   
   Subroutine GradientV_AssemblyBlk_SurfaceAT1(GradientV_Sec, iBlk, V_Sec, AppCtx)
      Type(SectionReal)                            :: GradientV_Sec
      PetscInt                                     :: iBlk
      Type(SectionReal)                            :: V_Sec
      Type(AppCtx_Type)                            :: AppCtx

      PetscReal, Dimension(:), Pointer             :: V_Loc, GradientV_Loc
#if defined PB_2D
      Type(Vect2D)                                 :: GradV_Elem
#elif defined PB_3D
      Type(Vect3D)                                 :: GradV_Elem
#endif      
      PetscInt                                     :: iE, iEloc, iBlkId, iErr
      PetscInt                                     :: NumDoFScal
      PetscInt                                     :: iDoF, iGauss
      PetscLogDouble                               :: flops

      flops        = 0.0

      NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF

      Allocate(V_Loc(NumDoFScal))
      Allocate(GradientV_Loc(NumDoFScal))

      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
      Do_iEloc: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         GradientV_Loc = 0.0_Kr
         Call SectionRealRestrictClosure(V_Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, V_Loc, iErr); CHKERRQ(ierr)
         Do_iGauss: Do iGauss = 1, size(AppCtx%ElemVect(iE)%Gauss_C)
            GradV_Elem = 0.0_Kr
            Do iDoF = 1, NumDoFScal
               GradV_Elem = GradV_Elem + V_Loc(iDoF) * AppCtx%ElemScal(iE)%Grad_BF(iDoF, iGauss)
            End Do
            Do iDoF = 1, NumDofScal
               GradientV_Loc(iDoF) = GradientV_Loc(iDoF) + AppCtx%ElemScal(iE)%Gauss_C(iGauss) *  (-AppCtx%ElemScal(iE)%BF(iDoF, iGauss) / AppCtx%VarFracSchemeParam%Epsilon +  2.0_Kr * AppCtx%VarFracSchemeParam%Epsilon * (GradV_Elem .DotP. AppCtx%ElemScal(iE)%Grad_BF(iDoF, iGauss)) )
            End Do
         End Do Do_iGauss
         GradientV_Loc = GradientV_Loc * AppCtx%MatProp(iBlkID)%Toughness / AppCtx%VarFracSchemeParam%ATCv * 0.25_Kr
         Call SectionRealUpdateClosure(GradientV_Sec, AppCtx%MeshTopology%Mesh, iE-1, GradientV_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iEloc

      DeAllocate(V_Loc)
      DeAllocate(GradientV_Loc)
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
      CHKMEMQ
   End Subroutine GradientV_AssemblyBlk_SurfaceAT1

   Subroutine GradientV_AssemblyBlk_SurfaceAT2(GradientV_Sec, iBlk, V_Sec, AppCtx)
      Type(SectionReal)                            :: GradientV_Sec
      PetscInt                                     :: iBlk
      Type(SectionReal)                            :: V_Sec
      Type(AppCtx_Type)                            :: AppCtx

      PetscReal, Dimension(:), Pointer             :: V_Loc, GradientV_Loc
      PetscReal                                    :: V_Elem
#if defined PB_2D
      Type(Vect2D)                                 :: GradV_Elem
#elif defined PB_3D
      Type(Vect3D)                                 :: GradV_Elem
#endif      
      PetscInt                                     :: iE, iEloc, iBlkId, iErr
      PetscInt                                     :: NumDoFScal
      PetscInt                                     :: iDoF, iGauss
      PetscLogDouble                               :: flops

      flops        = 0.0

      NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF

      Allocate(V_Loc(NumDoFScal))
      Allocate(GradientV_Loc(NumDoFScal))

      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
      Do_iEloc: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         GradientV_Loc = 0.0_Kr
         Call SectionRealRestrictClosure(V_Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, V_Loc, iErr); CHKERRQ(ierr)
         Do_iGauss: Do iGauss = 1, size(AppCtx%ElemVect(iE)%Gauss_C)
            V_Elem = 0.0_Kr
            GradV_Elem = 0.0_Kr
            Do iDoF = 1, NumDoFScal
               V_Elem     = V_Elem     + V_Loc(iDoF) * AppCtx%ElemScal(iE)%BF(iDoF, iGauss)
               GradV_Elem = GradV_Elem + V_Loc(iDoF) * AppCtx%ElemScal(iE)%Grad_BF(iDoF, iGauss)
               flops = flops + 2.0
            End Do
            Do iDoF = 1, NumDofScal
               GradientV_Loc(iDoF) = GradientV_Loc(iDoF) + AppCtx%ElemScal(iE)%Gauss_C(iGauss) * ( (V_Elem-1.0_Kr) * AppCtx%ElemScal(iE)%BF(iDoF, iGauss) / AppCtx%VarFracSchemeParam%Epsilon + AppCtx%VarFracSchemeParam%Epsilon * (GradV_Elem .DotP. AppCtx%ElemScal(iE)%Grad_BF(iDoF, iGauss)) )
            End Do
         End Do Do_iGauss
         GradientV_Loc = GradientV_Loc * AppCtx%MatProp(iBlkID)%Toughness / AppCtx%VarFracSchemeParam%ATCv * 0.5_Kr
         Call SectionRealUpdateClosure(GradientV_Sec, AppCtx%MeshTopology%Mesh, iE-1, GradientV_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iEloc

      DeAllocate(V_Loc)
      DeAllocate(GradientV_Loc)
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
      CHKMEMQ
   End Subroutine GradientV_AssemblyBlk_SurfaceAT2
#endif

   
   Subroutine Step_V(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      KSPConvergedReason                           :: reason
      PetscInt                                     :: KSPNumIter
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
      Type(Vec)                                    :: RHSV_Vec, V_Vec, V_Old
      PetscReal                                    :: VMin, VMax
      PetscReal                                    :: rDum
      PetscInt                                     :: iDum
#if defined WITH_TAO
      TaoTerminateReason                           :: TaoReason
      PetscReal                                    :: TaoResidual
#endif
      
      Call PetscLogStagePush(AppCtx%LogInfo%KSPSolveV_Stage, iErr); CHKERRQ(iErr)
  
      !!! Create Vectors for V and RHSV

      If (AppCtx%VarFracSchemeParam%V_UseTao) Then
#if defined WITH_TAO
         Call TaoAppGetSolutionVec(AppCtx%taoappV, V_Vec, iErr); CHKERRQ(iErr)
         Call SectionRealToVec(AppCtx%V, AppCtx%ScatterScal, SCATTER_FORWARD, V_Vec, ierr); CHKERRQ(ierr)
         Call VecDuplicate(V_Vec, V_Old, iErr); CHKERRQ(iErr)
         Call VecCopy(V_Vec, V_Old, iErr); CHKERRQ(iErr)
         If (AppCtx%AppParam%verbose > 0) Then
            Write(IOBuffer, *) 'Calling TaoSolveApplication\n'
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         End If
         Call TaoSolveApplication(AppCtx%taoappV, AppCtx%taoV, iErr); CHKERRQ(iErr)
         Call TaoGetSolutionStatus(AppCtx%taoV, KSPNumIter, rDum, TaoResidual, rDum, rDum, TaoReason, iErr); CHKERRQ(iErr)
         If ( TaoReason > 0) Then
            Write(IOBuffer, 102) KSPNumiter, TAOReason
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Else
            Write(IOBuffer, 103) TaoReason
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         End If     
         Call KSPGetConvergedReason(AppCtx%KSPV, reason, iErr); CHKERRQ(iErr)
         If ( reason > 0) Then
            Call KSPGetIterationNumber(AppCtx%KSPV, KSPNumIter, iErr); CHKERRQ(iErr)
            Write(IOBuffer, 100) KSPNumIter, reason
         Else
            Write(IOBuffer, 101) reason
         End If
#endif      
      Else
         Call MeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%V, V_Vec, iErr); CHKERRQ(iErr)
         Call MeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%V, RHSV_Vec, iErr); CHKERRQ(iErr)
         Call SectionRealToVec(AppCtx%V, AppCtx%ScatterScal, SCATTER_FORWARD, V_Vec, ierr); CHKERRQ(ierr)
         Call VecDuplicate(V_Vec, V_Old, iErr); CHKERRQ(iErr)
         Call VecCopy(V_Vec, V_Old, iErr); CHKERRQ(iErr)

         If (AppCtx%AppParam%verbose > 0) Then
            Write(IOBuffer, *) 'Assembling the Matrix and RHS for the V-subproblem \n' 
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
         End If
         Call RHSV_Assembly(RHSV_Vec, AppCtx)
         If (AppCtx%AppParam%verbose > 2) Then
            Call VecView(RHSV_Vec, AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
         End If

         Call MatV_Assembly(AppCtx%KV, AppCtx)
         If (AppCtx%AppParam%verbose > 2) Then
            Call MatView(AppCtx%KV, AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
         End If
         
         If (AppCtx%AppParam%verbose > 0) Then
            Write(IOBuffer, *) 'Calling KSPSolve for the V-subproblem\n' 
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
         End If

         Call KSPSolve(AppCtx%KSPV, RHSV_Vec, V_Vec, iErr); CHKERRQ(iErr)
      
         Call KSPGetConvergedReason(AppCtx%KSPV, reason, iErr); CHKERRQ(iErr)
         If ( reason > 0) Then
            Call KSPGetIterationNumber(AppCtx%KSPV, KSPNumIter, iErr); CHKERRQ(iErr)
            Write(IOBuffer, 100) KSPNumIter, reason
         Else
            Write(IOBuffer, 101) reason
         End If
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Call VecDestroy(RHSV_Vec, iErr); CHKERRQ(iErr)
      End If      
   
      !!! Compute ||v-v_old||_\infty
      Call VecMin(V_Vec, PETSC_NULL_INTEGER, VMin, iErr); CHKERRQ(iErr)
      Call VecMax(V_Vec, PETSC_NULL_INTEGER, VMax, iErr); CHKERRQ(iErr)
      Write(IOBuffer, 700) VMin, VMax
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

       Call VecAxPy(V_Old, -1.0_Kr, V_Vec, iErr)
       Call VecNorm(V_Old, NORM_INFINITY, AppCtx%ErrV, iErr)

      Write(IOBuffer, 800) AppCtx%ErrV
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

      !!! Scatter V back into a SectionReal
      Call SectionRealToVec(AppCtx%V, AppCtx%ScatterScal, SCATTER_REVERSE, V_Vec, ierr); CHKERRQ(ierr)
      
      Call VecDestroy(V_Old, iErr); CHKERRQ(iErr)
      If (.NOT. AppCtx%VarFracSchemeParam%V_UseTao) Then
         Call VecDestroy(V_Vec, iErr); CHKERRQ(iErr)
      End If
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
100 Format('     KSP for V converged in  ', I5, ' iterations. KSPConvergedReason is     ', I5, '\n')
101 Format('[ERROR] KSP for V diverged. KSPConvergedReason is ', I2, '\n')
102 Format('     TAO solver converged in ', I5, ' iterations. Tao termination reason is ', I5, '\n')
103 Format('[ERROR] TaoSolveApplication did not converge. ', I2, '\n')      
700 Format('     VMin / Max:   ', T24, 2(ES12.5, '  '), '\n')
800 Format('     Max change V: ', T24, ES12.5, '\n')
      CHKMEMQ
   End Subroutine Step_V

   
   
#if defined PB_2D
End Module m_VarFracQS_V2D
#elif defined PB_3D
End Module m_VarFracQS_V3D
#endif
