#if defined PB_2D
Module m_VarFracQS_U2D
#elif defined PB_3D
Module m_VarFracQS_U3D
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
   Public :: HessianU_Assembly
   Public :: FormFunctionAndGradientU
   Public :: InitTaoBoundsU
#endif
   Public :: Init_TS_U
   Public :: MatU_Assembly
   Public :: RHSU_Assembly
   Public :: Step_U
   
Contains
#if defined WITH_TAO
   Subroutine InitTaoBoundsU(TaoApp, LowerBoundU_Vec, UpperBoundU_Vec, AppCtx, iErr)
      TAO_APPLICATION                              :: taoapp
      Type(Vec)                                    :: LowerBoundU_Vec, UpperBoundU_Vec
      Type(AppCtx_Type)                            :: AppCtx
      PetscErrorCode                               :: iErr
      
      Type(SectionReal)                            :: LowerBoundU_Sec, UpperBoundU_Sec
      PetscInt, Dimension(:), Pointer              :: BCUFlag_Ptr
      PetscReal, Dimension(:), Pointer             :: LowerBoundU_Ptr, UpperBoundU_Ptr, U_Ptr
      PetscInt                                     :: i, j
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
      PetscReal                                    :: DummyLowerBound = -1.0D+20      
      PetscReal                                    :: DummyUpperBound =  1.0D+20      
      
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) "Updating bounds for U\n"c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'LowerBoundU_Sec', AppCtx%MeshTopology%Num_Dim, LowerBoundU_Sec , iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'UpperBoundU_Sec', AppCtx%MeshTopology%Num_Dim,  UpperBoundU_Sec , iErr); CHKERRQ(iErr) 
      Call SectionRealSet(LowerBoundU_Sec, DummyLowerBound, iErr); CHKERRQ(iErr)
      Call SectionRealSet(UpperBoundU_Sec, DummyUpperBound, iErr); CHKERRQ(iErr)
      
      Allocate(LowerBoundU_Ptr(AppCtx%MeshTopology%Num_Dim))
      Allocate(UpperBoundU_Ptr(AppCtx%MeshTopology%Num_Dim))
      Do i = 1, AppCtx%MeshTopology%Num_Verts !!! Baaad
         Call SectionIntRestrict(AppCtx%BCUFlag, AppCtx%MeshTopology%Num_Elems + i-1, BCUFlag_Ptr, iErr); CHKERRQ(ierr)
         If (Sum(BCUFlag_Ptr) /= 0) Then
            LowerBoundU_Ptr = DummyLowerBound
            UpperBoundU_Ptr = DummyUpperBound
            Call SectionRealRestrict(AppCtx%U, AppCtx%MeshTopology%Num_Elems + i-1, U_Ptr, iErr); CHKERRQ(iErr)      
            Do j = 1, AppCtx%MeshTopology%Num_Dim      
               If (BCUFlag_Ptr(j) /= VarFrac_BC_Type_NONE) Then
                  LowerBoundU_Ptr(j) = U_Ptr(j)
                  UpperBoundU_Ptr(j) = U_Ptr(j)
               End If
            End Do
            Call SectionRealUpdate(LowerBoundU_Sec, AppCtx%MeshTopology%Num_Elems + i-1, LowerBoundU_Ptr, INSERT_VALUES, iErr); CHKERRQ(iErr)
            Call SectionRealUpdate(UpperBoundU_Sec, AppCtx%MeshTopology%Num_Elems + i-1, UpperBoundU_Ptr, INSERT_VALUES, iErr); CHKERRQ(iErr)
            Call SectionRealRestore(AppCtx%U, AppCtx%MeshTopology%Num_Elems + i-1, U_Ptr, iErr); CHKERRQ(iErr)      
         End If
         Call SectionIntRestore(AppCtx%BCUFlag, AppCtx%MeshTopology%Num_Elems + i-1, BCUFlag_Ptr, iErr); CHKERRQ(iErr)
      End Do
      DeAllocate(LowerBoundU_Ptr)   
      DeAllocate(UpperBoundU_Ptr)   
         
      Call SectionRealComplete(LowerBoundU_Sec, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(LowerBoundU_Sec, AppCtx%ScatterVect, SCATTER_FORWARD, LowerBoundU_Vec, iErr); CHKERRQ(iErr)
      Call SectionRealComplete(UpperBoundU_Sec, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(UpperBoundU_Sec, AppCtx%ScatterVect, SCATTER_FORWARD, UpperBoundU_Vec, iErr); CHKERRQ(iErr)

      Call SectionRealDestroy(LowerBoundU_Sec, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(UpperBoundU_Sec, iErr); CHKERRQ(iErr)
      CHKMEMQ
   End Subroutine InitTaoBoundsU
#endif

   Subroutine Init_TS_U(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      Type(SectionReal)                            :: UBC
      PetscReal, Dimension(:), Pointer             :: U_Ptr, UBC_Ptr
      PetscInt, Dimension(:), Pointer              :: BCFlag
      PetscInt                                     :: i, j, iErr
      
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'UBC', AppCtx%MeshTopology%Num_Dim, UBC, iErr); CHKERRQ(iErr)
      !!! Using SectionRealDuplicate would make more sense
      
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_DisplacementX)%Offset, AppCtx%TimeStep, UBC)

      Do i = 1, AppCtx%MeshTopology%Num_Verts
         Call SectionIntRestrict(AppCtx%BCUFlag, AppCtx%MeshTopology%Num_Elems + i-1, BCFlag, iErr); CHKERRQ(ierr)
         If (Sum(BCFlag) /= 0) Then
            Call SectionRealRestrict(AppCtx%U, AppCtx%MeshTopology%Num_Elems + i-1, U_Ptr, iErr); CHKERRQ(iErr)      
            Call SectionRealRestrict(UBC, AppCtx%MeshTopology%Num_Elems + i-1, UBC_Ptr, iErr); CHKERRQ(iErr)
            Do j = 1, AppCtx%MeshTopology%Num_Dim      
               If (BCFlag(j) /= 0) Then
                  U_Ptr(j) = UBC_Ptr(j)
               End If
            End Do
            Call SectionRealUpdate(AppCtx%U, AppCtx%MeshTopology%Num_Elems + i-1, UBC_Ptr, INSERT_VALUES, iErr); CHKERRQ(iErr)
            Call SectionRealRestore(AppCtx%U, AppCtx%MeshTopology%Num_Elems + i-1, U_Ptr, iErr); CHKERRQ(iErr)      
            Call SectionRealRestore(UBC, AppCtx%MeshTopology%Num_Elems + i-1, UBC_Ptr, iErr); CHKERRQ(iErr)
         End If
         Call SectionIntRestore(AppCtx%BCUFlag, AppCtx%MeshTopology%Num_Elems + i-1, BCFlag, iErr); CHKERRQ(iErr)
      End Do
      Call SectionRealDestroy(UBC, iErr); CHKERRQ(iErr)
   End Subroutine Init_TS_U
   
!!!
!!! Global Assembly Functions
!!! 

   Subroutine MatU_Assembly(K, AppCtx)
      Type(Mat)                                    :: K
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iBlk, iBlkID, iErr
      
      Call PetscLogStagePush(AppCtx%LogInfo%MatAssemblyV_Stage, iErr); CHKERRQ(iErr)

      Call MatZeroEntries(K, iErr); CHKERRQ(iErr)
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
         If (AppCtx%MyEXO%EBProperty(VarFrac_EBProp_IsBrittle)%Value(iBlkID) /= 0) Then
            Call MatU_AssemblyBlk_Brittle(K, iBlk, .TRUE., AppCtx)
         Else
            Call MatU_AssemblyBlk_NonBrittle(K, iBlk, .TRUE., AppCtx)
         EndIf
      End Do Do_Elem_iBlk
      Call MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      If (AppCtx%AppParam%verbose > 1) Then
         Call MatView(K, PETSC_VIEWER_STDOUT_WORLD, iErr)
      End If

      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
      CHKMEMQ
   End Subroutine MatU_Assembly

#if defined WITH_TAO
   Subroutine HessianU_Assembly(tao, X, H, Hpre, flg, AppCtx, iErr)
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
         If (AppCtx%MyEXO%EBProperty(VarFrac_EBProp_IsBrittle)%Value(iBlkID) /= 0) Then
            Call MatU_AssemblyBlk_Brittle(H, iBlk, .FALSE., AppCtx)
         Else
            Call MatU_AssemblyBlk_NonBrittle(H, iBlk, .FALSE., AppCtx)
         EndIf
      End Do Do_Elem_iBlk
      Call MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (H, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      If (AppCtx%AppParam%verbose > 1) Then
         Call MatView(H, PETSC_VIEWER_STDOUT_WORLD, iErr)
      End If

      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
      CHKMEMQ
   End Subroutine HessianU_Assembly
#endif
  
#if defined WITH_TAO
   Subroutine FormFunctionAndGradientU(tao, U_Vec, ObjFunc, GradientU_Vec, AppCtx, iErr)
      TAO_SOLVER                                   :: tao
      Type(Vec)                                    :: U_Vec, GradientU_Vec
      PetscInt                                     :: iErr
      PetscReal                                    :: ObjFunc
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iBlk, iBlkId
      Type(SectionReal)                            :: GradientU_Sec
      PetscReal                                    :: MyElasticEnergyBlock, MyExtForcesWorkBlock
      PetscReal                                    :: MyObjFunc
      
      !!! Objective function is ElasticEnergy + ExtForcesWork
      Call SectionRealToVec(AppCtx%U, AppCtx%ScatterVect, SCATTER_REVERSE, U_Vec, iErr); CHKERRQ(ierr)

      MyObjFunc = 0.0_Kr
      Do_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
         If (AppCtx%MyEXO%EBProperty(VarFrac_EBProp_IsBrittle)%Value(iBlkID) /= 0) Then
            Call ElasticEnergy_AssemblyBlk_Brittle(MyElasticEnergyBlock, iBlk, AppCtx%U, AppCtx%Theta, AppCtx%V, AppCtx)
         Else
            Call ElasticEnergy_AssemblyBlk_NonBrittle(MyElasticEnergyBlock, iBlk, AppCtx%U, AppCtx%Theta, AppCtx)
         End If
         MyObjFunc = MyObjFunc + MyElasticEnergyBlock

         If (AppCtx%MyEXO%EBProperty(VarFrac_EBProp_HasBForce)%Value(iBlkID) /= 0) Then
            Call ElasticEnergy_AssemblyBlk_Brittle(MyExtForcesWorkBlock, iBlk, AppCtx%U, AppCtx%Theta, AppCtx%F, AppCtx)
            MyObjFunc = MyObjFunc + MyExtForcesWorkBlock
         End If
      End Do Do_iBlk
      Call PetscGlobalSum(MyObjFunc, ObjFunc, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      
      !!! Gradient
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'GradientU', AppCtx%MeshTopology%Num_Dim, GradientU_Sec, iErr); CHKERRQ(iErr)
      Call SectionRealZero(GradientU_Sec, iErr); CHKERRQ(iErr)

      Do_iBlk2: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
         If (AppCtx%MyEXO%EBProperty(VarFrac_EBProp_IsBrittle)%Value(iBlkID) /= 0) Then
            Call GradientU_AssemblyBlk_ElastBrittle(GradientU_Sec, iBlk, AppCtx%U, AppCtx%Theta, AppCtx%V, AppCtx)
         Else
            Call GradientU_AssemblyBlk_ElastNonBrittle(GradientU_Sec, iBlk, AppCtx%U, AppCtx%Theta, AppCtx)
         End If
         If (AppCtx%MyEXO%EBProperty(VarFrac_EBProp_IsBrittle)%Value(iBlkID) /= 0) Then
            Call GradientU_AssemblyBlk_ExtForcesWork(GradientU_Sec, iBlk, AppCtx%U, AppCtx%F, AppCtx)
         End If
      End Do Do_iBlk2
      
      Call SectionRealComplete(GradientU_Sec, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(GradientU_Sec, AppCtx%ScatterVect, SCATTER_FORWARD, GradientU_Vec, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(GradientU_Sec, iErr); CHKERRQ(iErr)
      CHKMEMQ
   End Subroutine FormFunctionAndGradientU
#endif


   Subroutine RHSU_Assembly(RHS_Vec, AppCtx)
      Type(Vec)                                    :: RHS_Vec
      Type(AppCtx_Type)                            :: AppCtx
      
      Type(SectionReal)                            :: RHSU_Sec
      PetscInt                                     :: iErr
      PetscInt                                     :: iBlk, iBlkID

      Call PetscLogStagePush(AppCtx%LogInfo%RHSAssemblyU_Stage, iErr); CHKERRQ(iErr)

      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'RHSU_Sec',  AppCtx%MeshTopology%Num_Dim, RHSU_Sec , iErr); CHKERRQ(iErr)
      Call SectionRealZero(RHSU_Sec, iErr); CHKERRQ(iErr)
      
      Do_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
         If (AppCtx%MyEXO%EBProperty(VarFrac_EBProp_IsBrittle)%Value(iBlkID) /= 0) Then
            Call RHSAssemblyBlock_ElastBrittle(RHSU_Sec, iBlk, AppCtx)
         Else
            Call RHSAssemblyBlock_ElastNonBrittle(RHSU_Sec, iBlk, AppCtx)
         End If         
         If (AppCtx%MyEXO%EBProperty(VarFrac_EBProp_HasBForce)%Value(iBlkID) /= 0) Then
            Call RHSAssemblyBlock_Force(RHSU_Sec, iBlk, AppCtx)
         End If
      End Do Do_iBlk

      Call SectionRealComplete(RHSU_Sec, iErr); CHKERRQ(iErr)
      !!! VERY important! This is the equivalent of a ghost update
      Call SectionRealToVec(RHSU_Sec, AppCtx%ScatterVect, SCATTER_FORWARD, RHS_Vec, iErr); CHKERRQ(ierr)
      Call SectionRealDestroy(RHSU_Sec, iErr); CHKERRQ(iErr)

      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine RHSU_Assembly
  
  
   !!! 
   !!! Block Assembly Routines
   !!!
   Subroutine MatU_AssemblyBlk_Brittle(K, iBlk, DoBC, AppCtx)
      Type(Mat)                                    :: K
      PetscInt                                     :: iBlk
      PetscTruth                                   :: DoBC
      Type(AppCtx_Type)                            :: AppCtx

      PetscInt                                     :: iBlkID
      PetscReal, Dimension(:,:), Pointer           :: Mat_Loc      
      PetscInt                                     :: iELoc, iE
      PetscInt                                     :: iErr

      PetscInt                                     :: NumDoFScal, NumDoFVect, NumGauss
      PetscInt, Dimension(:), Pointer              :: BCFlag_Loc
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      
      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_Elem = \sum_i U_Loc(i) BF(i))
      PetscReal, Dimension(:), Pointer             :: V_Loc
      PetscReal                                    :: V_Elem, CoefV
      PetscLogDouble                               :: flops
      
      Call PetscLogEventBegin(AppCtx%LogInfo%MatAssemblyLocalU_Event, iErr); CHKERRQ(iErr)

      flops = 0.0
      NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
      NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
      
      Allocate(V_Loc(NumDoFScal))
      Allocate(BCFlag_Loc(NumDoFVect))
      BCFlag_Loc = VarFrac_BC_Type_NONE
      Allocate(Mat_Loc(NumDoFVect, NumDoFVect))
      
      Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         
         Mat_Loc  = 0.0_Kr
         Call SectionRealRestrictClosure(AppCtx%V, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, V_Loc, iErr); CHKERRQ(ierr)
         If (DoBC) Then
            Call SectionIntRestrictClosure(AppCtx%BCUFlag, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, BCFlag_Loc, iErr); CHKERRQ(ierr)
         End If
            
         Do iGauss = 1, Size(AppCtx%ElemVect(iE)%Gauss_C)
            !!! Compute the contribution of V to the stiffness matrix
            !!! CoefV = (1+\eta_\varepsilon)v^2 if Is_Brittle, 1 otherwise
            !! Calculate V at the gauss point
            V_Elem = 0.0_Kr        
            Do iDoF1 = 1, NumDoFScal
               V_Elem = V_Elem + AppCtx%ElemScal(iE)%BF(iDoF1, iGauss) * V_Loc(iDoF1)
               flops = flops + 2.0
            End Do
            CoefV = V_Elem**2 + AppCtx%VarFracSchemeParam%KEpsilon
            flops = flops + 2.0
            Do iDoF1 = 1, NumDoFVect
               If (BCFlag_Loc(iDoF1) == VarFrac_BC_Type_NONE) Then
                  Do iDoF2 = 1, NumDoFVect
                     Mat_Loc(iDoF2, iDoF1) =  Mat_Loc(iDoF2, iDoF1) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * CoefV * ((AppCtx%MatProp(iBlkID)%Hookes_Law * AppCtx%ElemVect(iE)%GradS_BF(iDoF1, iGauss)) .DotP. AppCtx%ElemVect(iE)%GradS_BF(iDoF2, iGauss))
                     flops = flops + 3.0
                  End Do
               End If
            End Do
         End Do
         Call assembleMatrix(AppCtx%KU, AppCtx%MeshTopology%mesh, AppCtx%U, iE-1, Mat_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_Elem_iE
      
      Call PetscLogFlops(flops, iErr); CHKERRQ(iErr)
      DeAllocate(Mat_Loc)
      DeAllocate(V_Loc)
      DeAllocate(BCFlag_Loc)
      Call PetscLogEventEnd(AppCtx%LogInfo%MatAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
   End Subroutine MatU_AssemblyBlk_Brittle

   Subroutine MatU_AssemblyBlk_NonBrittle(K, iBlk, DoBC, AppCtx)
      Type(Mat)                                    :: K
      PetscInt                                     :: iBlk
      PetscTruth                                   :: DoBC
      Type(AppCtx_Type)                            :: AppCtx

      PetscInt                                     :: iBlkID
      PetscReal, Dimension(:,:), Pointer           :: Mat_Loc      
      PetscInt                                     :: iELoc, iE
      PetscInt                                     :: iErr

      PetscInt                                     :: NumDoFScal, NumDoFVect, NumGauss
      PetscInt, Dimension(:), Pointer              :: BCFlag_Loc
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      
      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_Elem = \sum_i U_Loc(i) BF(i))
      PetscLogDouble                               :: flops
      
      Call PetscLogEventBegin(AppCtx%LogInfo%MatAssemblyLocalU_Event, iErr); CHKERRQ(iErr)

      flops = 0.0
      NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
      NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
      
      Allocate(BCFlag_Loc(NumDoFVect))
      BCFlag_Loc = VarFrac_BC_Type_NONE
      Allocate(Mat_Loc(NumDoFVect, NumDoFVect))
      
      Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         
         Mat_Loc  = 0.0_Kr
         If (DoBC) Then
            Call SectionIntRestrictClosure(AppCtx%BCUFlag, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, BCFlag_Loc, iErr); CHKERRQ(ierr)
         End If
            
         Do iGauss = 1, Size(AppCtx%ElemVect(iE)%Gauss_C)
            Do iDoF1 = 1, NumDoFVect
               If (BCFlag_Loc(iDoF1) == VarFrac_BC_Type_NONE) Then
                  Do iDoF2 = 1, NumDoFVect
                     Mat_Loc(iDoF2, iDoF1) =  Mat_Loc(iDoF2, iDoF1) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * ((AppCtx%MatProp(iBlkID)%Hookes_Law * AppCtx%ElemVect(iE)%GradS_BF(iDoF1, iGauss)) .DotP. AppCtx%ElemVect(iE)%GradS_BF(iDoF2, iGauss))
                     flops = flops + 3.0
                  End Do
               End If
            End Do
         End Do
         Call assembleMatrix(AppCtx%KU, AppCtx%MeshTopology%mesh, AppCtx%U, iE-1, Mat_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_Elem_iE
      
      Call PetscLogFlops(flops, iErr); CHKERRQ(iErr)
      DeAllocate(Mat_Loc)
      DeAllocate(BCFlag_Loc)
      Call PetscLogEventEnd(AppCtx%LogInfo%MatAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
   End Subroutine MatU_AssemblyBlk_NonBrittle
   
   Subroutine GradientU_AssemblyBlk_ElastBrittle(Gradient_Sec, iBlk, X_Sec, Theta_Sec, V_Sec, AppCtx)
      Type(SectionReal)                            :: Gradient_Sec
      PetscInt                                     :: iBlk
      Type(SectionReal)                            :: X_Sec, Theta_Sec, V_Sec
      Type(AppCtx_Type)                            :: AppCtx
      
      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_ELem = \sum_i U_Loc(i) BF(i))
      PetscReal, Dimension(:), Pointer             :: X_Loc, Theta_Loc, V_Loc, Gradient_Loc
#if defined PB_2D
      Type(Vect2D)                                 :: X_Elem
      Type(Mats2D)                                 :: Strain_Elem, EffectiveStrain_Elem
#elif defined PB_3D  
      Type(Vect3D)                                 :: X_Elem
      Type(Mats3D)                                 :: Strain_Elem, EffectiveStrain_Elem
#endif
      PetscReal                                    :: Theta_Elem, V_Elem, CoefV
      PetscInt                                     :: iE, iEloc, iBlkId, iErr
      PetscInt                                     :: NumDoFScal, NumDoFVect
      PetscInt                                     :: iDoF, iGauss
      PetscLogDouble                               :: flops       

      Call PetscLogEventBegin(AppCtx%LogInfo%RHSAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
      flops        = 0.0

      NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
      NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF

      Allocate(X_Loc(NumDoFVect))
      Allocate(Gradient_Loc(NumDoFVect))
      Allocate(Theta_Loc(NumDoFScal))
      Allocate(V_Loc(NumDoFScal))

      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID

      Do_iEloc: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         Gradient_Loc = 0.0_Kr
         Call SectionRealRestrictClosure(X_Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, X_Loc, iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(V_Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, V_Loc, iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(Theta_Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, Theta_Loc, iErr); CHKERRQ(ierr)
         Do_iGauss: Do iGauss = 1, size(AppCtx%ElemVect(iE)%Gauss_C)
            X_Elem = 0.0_Kr
            Strain_Elem = 0.0_Kr
            Do iDoF = 1, NumDoFVect
               X_Elem      = X_Elem + X_Loc(iDoF) * AppCtx%ElemVect(iE)%BF(iDoF, iGauss)
               Strain_Elem = Strain_Elem + X_Loc(iDoF) * AppCtx%ElemVect(iE)%GradS_BF(iDoF, iGauss)
            End Do
            Theta_Elem = 0.0_Kr
            V_Elem = 0.0_Kr        
            Do iDoF = 1, NumDoFScal
               V_Elem = V_Elem + AppCtx%ElemScal(iE)%BF(iDoF, iGauss) * V_Loc(iDoF)
               Theta_Elem = Theta_Elem + AppCtx%ElemScal(iE)%BF(iDoF, iGauss) * Theta_Loc(iDoF)
               flops = flops + 4.0
            End Do
            CoefV = V_Elem**2 + AppCtx%VarFracSchemeParam%KEpsilon
            !!! CoefV = (1+\eta_\varepsilon)v^2 if Is_Brittle, 1 otherwise
            flops = flops + 2.0

            EffectiveStrain_Elem = Strain_Elem - AppCtx%MatProp(iBlkId)%Therm_Exp * Theta_Elem
            
            Do iDoF = 1, NumDofVect
               Gradient_Loc(iDoF) = Gradient_Loc(iDoF) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * CoefV * ((AppCtx%MatProp(iBlkId)%Hookes_Law * EffectiveStrain_Elem) .DotP. AppCtx%ElemVect(iE)%GradS_BF(iDoF, iGauss))
            End Do
         End Do Do_iGauss
         Call SectionRealUpdateClosure(Gradient_Sec, AppCtx%MeshTopology%Mesh, iE-1, Gradient_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iEloc

      DeAllocate(X_Loc)
      DeAllocate(Gradient_Loc)
      DeAllocate(Theta_Loc)
      DeAllocate(V_Loc)
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
      Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
   End Subroutine GradientU_AssemblyBlk_ElastBrittle

   Subroutine GradientU_AssemblyBlk_ElastNonBrittle(Gradient_Sec, iBlk, X_Sec, Theta_Sec, AppCtx)
      Type(SectionReal)                            :: Gradient_Sec
      PetscInt                                     :: iBlk
      Type(SectionReal)                            :: X_Sec, Theta_Sec
      Type(AppCtx_Type)                            :: AppCtx
      
      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_ELem = \sum_i U_Loc(i) BF(i))
      PetscReal, Dimension(:), Pointer             :: X_Loc, Theta_Loc, Gradient_Loc
#if defined PB_2D
      Type(Vect2D)                                 :: X_Elem
      Type(Mats2D)                                 :: Strain_Elem, EffectiveStrain_Elem
#elif defined PB_3D  
      Type(Vect3D)                                 :: X_Elem
      Type(Mats3D)                                 :: Strain_Elem, EffectiveStrain_Elem
#endif
      PetscReal                                    :: Theta_Elem
      PetscInt                                     :: iE, iEloc, iBlkId, iErr
      PetscInt                                     :: NumDoFScal, NumDoFVect
      PetscInt                                     :: iDoF, iGauss
      PetscLogDouble                               :: flops       

      Call PetscLogEventBegin(AppCtx%LogInfo%RHSAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
      flops        = 0.0

      NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
      NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF

      Allocate(X_Loc(NumDoFVect))
      Allocate(Gradient_Loc(NumDoFVect))
      Allocate(Theta_Loc(NumDoFScal))

      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID

      Do_iEloc: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         Gradient_Loc = 0.0_Kr
         Call SectionRealRestrictClosure(X_Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, X_Loc, iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(Theta_Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, Theta_Loc, iErr); CHKERRQ(ierr)
         Do_iGauss: Do iGauss = 1, size(AppCtx%ElemVect(iE)%Gauss_C)
            X_Elem = 0.0_Kr
            Strain_Elem = 0.0_Kr
            Do iDoF = 1, NumDoFVect
               X_Elem      = X_Elem + X_Loc(iDoF) * AppCtx%ElemVect(iE)%BF(iDoF, iGauss)
               Strain_Elem = Strain_Elem + X_Loc(iDoF) * AppCtx%ElemVect(iE)%GradS_BF(iDoF, iGauss)
            End Do
            Theta_Elem = 0.0_Kr
            Do iDoF = 1, NumDoFScal
               Theta_Elem = Theta_Elem + AppCtx%ElemScal(iE)%BF(iDoF, iGauss) * Theta_Loc(iDoF)
               flops = flops + 2.0
            End Do
            EffectiveStrain_Elem = Strain_Elem - AppCtx%MatProp(iBlkId)%Therm_Exp * Theta_Elem
            
            Do iDoF = 1, NumDofVect
               Gradient_Loc(iDoF) = Gradient_Loc(iDoF) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * ((AppCtx%MatProp(iBlkId)%Hookes_Law * EffectiveStrain_Elem) .DotP. AppCtx%ElemVect(iE)%GradS_BF(iDoF, iGauss))
            End Do
         End Do Do_iGauss
         Call SectionRealUpdateClosure(Gradient_Sec, AppCtx%MeshTopology%Mesh, iE-1, Gradient_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iEloc

      DeAllocate(X_Loc)
      DeAllocate(Gradient_Loc)
      DeAllocate(Theta_Loc)
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
      Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
   End Subroutine GradientU_AssemblyBlk_ElastNonBrittle

   Subroutine GradientU_AssemblyBlk_ExtForcesWork(Gradient_Sec, iBlk, X_Sec, F_Sec, AppCtx)
      Type(SectionReal)                            :: Gradient_Sec
      PetscInt                                     :: iBlk
      Type(SectionReal)                            :: X_Sec, F_Sec
      Type(AppCtx_Type)                            :: AppCtx
      
      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_ELem = \sum_i U_Loc(i) BF(i))
      PetscReal, Dimension(:), Pointer             :: X_Loc, F_Loc, Gradient_Loc
#if defined PB_2D
      Type(Vect2D)                                 :: X_Elem, F_Elem
#elif defined PB_3D  
      Type(Vect3D)                                 :: X_Elem, F_Elem    
#endif
      PetscInt                                     :: iE, iEloc, iBlkId, iErr
      PetscInt                                     :: NumDoFVect
      PetscInt                                     :: iDoF, iGauss
      PetscLogDouble                               :: flops       

      Call PetscLogEventBegin(AppCtx%LogInfo%RHSAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
      flops        = 0.0

      NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim

      Allocate(X_Loc(NumDoFVect))
      Allocate(F_Loc(NumDoFVect))
      Allocate(Gradient_Loc(NumDoFVect))

      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID

      Do_iEloc: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         Gradient_Loc = 0.0_Kr
         Call SectionRealRestrictClosure(X_Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, X_Loc, iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(F_Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, F_Loc, iErr); CHKERRQ(ierr)
         Do_iGauss: Do iGauss = 1, size(AppCtx%ElemVect(iE)%Gauss_C)
            F_Elem = 0.0_Kr
            X_Elem = 0.0_Kr
            Do iDoF = 1, NumDoFVect
               F_Elem      = F_Elem + AppCtx%ElemVect(iE)%BF(iDoF, iGauss) * F_Loc(iDoF)
               X_Elem      = X_Elem + X_Loc(iDoF) * AppCtx%ElemVect(iE)%BF(iDoF, iGauss)
            End Do
            
            Do iDoF = 1, NumDofVect
               Gradient_Loc(iDoF) = Gradient_Loc(iDoF) - AppCtx%ElemVect(iE)%Gauss_C(iGauss) * (F_Elem .DotP. AppCtx%ElemVect(iE)%BF(iDoF, iGauss))
            End Do
         End Do Do_iGauss
         Call SectionRealUpdateClosure(Gradient_Sec, AppCtx%MeshTopology%Mesh, iE-1, Gradient_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iEloc

      DeAllocate(X_Loc)
      DeAllocate(F_Loc)
      DeAllocate(Gradient_Loc)
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
      Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
   End Subroutine GradientU_AssemblyBlk_ExtForcesWork

   Subroutine RHSAssemblyBlock_ElastBrittle(RHS_Sec, iBlk, AppCtx)
      PetscInt                                     :: iBlk
      Type(SectionReal)                            :: RHS_Sec
      Type(AppCtx_Type)                            :: AppCtx

      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_ELem = \sum_i U_Loc(i) BF(i))
      PetscReal, Dimension(:), Pointer             :: Theta_Loc, RHS_Loc
      PetscInt, Dimension(:), Pointer              :: BCFlag_Loc
      PetscReal                                    :: Theta_Elem
      PetscInt                                     :: iE, iEloc, iBlkId, iErr
      PetscInt                                     :: NumDoFScal, NumDoFVect
      PetscInt                                     :: iDoF, iGauss
      PetscReal, Dimension(:), Pointer             :: V_Loc
      PetscReal                                    :: V_Elem, CoefV
      PetscLogDouble                               :: flops       

      Call PetscLogEventBegin(AppCtx%LogInfo%RHSAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
      flops      = 0.0
      
      NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
      NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF

      Allocate(V_Loc(NumDoFScal))
      Allocate(BCFlag_Loc(NumDoFVect))
      Allocate(RHS_Loc(NumDoFVect))
      Allocate(Theta_Loc(NumDoFScal))

      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
      Do_iEloc: Do iEloc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         RHS_Loc = 0.0_Kr
         Call SectionRealRestrictClosure(AppCtx%V, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, V_Loc, iErr); CHKERRQ(ierr)
         Call SectionIntRestrictClosure(AppCtx%BCUFlag, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, BCFlag_Loc, iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(AppCtx%Theta, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, Theta_Loc, iErr); CHKERRQ(ierr)
         Do_iGauss: Do iGauss = 1, size(AppCtx%ElemVect(iE)%Gauss_C)
            V_Elem = 0.0_Kr        
            Theta_Elem = 0.0_Kr
            Do iDoF = 1, NumDoFScal
               V_Elem = V_Elem + AppCtx%ElemScal(iE)%BF(iDoF, iGauss) * V_Loc(iDoF)
               Theta_Elem = Theta_Elem + AppCtx%ElemScal(iE)%BF(iDoF, iGauss) * Theta_Loc(iDoF)
               flops = flops + 4.0
            End Do
            CoefV = V_Elem**2 + AppCtx%VarFracSchemeParam%KEpsilon
            !!! CoefV = (1+\eta_\varepsilon)v^2 if Is_Brittle, 1 otherwise
            flops = flops + 2.0
            Do iDoF = 1, NumDoFVect
               If (BCFlag_Loc(iDoF) == 0) Then
                  ! RHS terms due to inelastic strains
                  RHS_Loc(iDoF) = RHS_Loc(iDoF) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * Theta_Elem * CoefV * ((AppCtx%MatProp(iBlkId)%Hookes_Law * AppCtx%ElemVect(iE)%GradS_BF(iDoF, iGauss)) .DotP. AppCtx%MatProp(iBlkId)%Therm_Exp)
                  flops = flops + 5.0
               End If
            End Do
         End Do Do_iGauss
         Call SectionRealUpdateClosure(RHS_Sec, AppCtx%MeshTopology%Mesh, iE-1, RHS_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iEloc
      
      DeAllocate(BCFlag_Loc)
      DeAllocate(RHS_Loc)
      DeAllocate(Theta_Loc)
      
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
      Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
   End Subroutine RHSAssemblyBlock_ElastBrittle

   Subroutine RHSAssemblyBlock_ElastNonBrittle(RHS_Sec, iBlk, AppCtx)
      PetscInt                                     :: iBlk
      Type(SectionReal)                            :: RHS_Sec
      Type(AppCtx_Type)                            :: AppCtx

      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_ELem = \sum_i U_Loc(i) BF(i))
      PetscReal, Dimension(:), Pointer             :: Theta_Loc, RHS_Loc
      PetscInt, Dimension(:), Pointer              :: BCFlag_Loc
      PetscReal                                    :: Theta_Elem
      PetscInt                                     :: iE, iEloc, iBlkId, iErr
      PetscInt                                     :: NumDoFScal, NumDoFVect
      PetscInt                                     :: iDoF, iGauss
      PetscLogDouble                               :: flops       

      Call PetscLogEventBegin(AppCtx%LogInfo%RHSAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
      flops      = 0.0
      
      NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
      NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF

      Allocate(BCFlag_Loc(NumDoFVect))
      Allocate(RHS_Loc(NumDoFVect))
      Allocate(Theta_Loc(NumDoFScal))

      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
      Do_iEloc: Do iEloc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         RHS_Loc = 0.0_Kr
         Call SectionIntRestrictClosure(AppCtx%BCUFlag, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, BCFlag_Loc, iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(AppCtx%Theta, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, Theta_Loc, iErr); CHKERRQ(ierr)
         Do_iGauss: Do iGauss = 1, size(AppCtx%ElemVect(iE)%Gauss_C)
            Theta_Elem = 0.0_Kr
            Do iDoF = 1, NumDoFScal
               Theta_Elem = Theta_Elem + AppCtx%ElemScal(iE)%BF(iDoF, iGauss) * Theta_Loc(iDoF)
               flops = flops + 2.0
            End Do
            Do iDoF = 1, NumDoFVect
               If (BCFlag_Loc(iDoF) == 0) Then
                  ! RHS terms due to inelastic strains
                  RHS_Loc(iDoF) = RHS_Loc(iDoF) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * Theta_Elem * ((AppCtx%MatProp(iBlkId)%Hookes_Law * AppCtx%ElemVect(iE)%GradS_BF(iDoF, iGauss)) .DotP. AppCtx%MatProp(iBlkId)%Therm_Exp)
                  flops = flops + 5.0
               End If
            End Do
         End Do Do_iGauss
         Call SectionRealUpdateClosure(RHS_Sec, AppCtx%MeshTopology%Mesh, iE-1, RHS_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iEloc
      
      DeAllocate(BCFlag_Loc)
      DeAllocate(RHS_Loc)
      DeAllocate(Theta_Loc)
      
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
      Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
   End Subroutine RHSAssemblyBlock_ElastNonBrittle
   
      Subroutine RHSAssemblyBlock_Force(RHS_Sec, iBlk, AppCtx)
      PetscInt                                     :: iBlk
      Type(SectionReal)                            :: RHS_Sec
      Type(AppCtx_Type)                            :: AppCtx

      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_ELem = \sum_i U_Loc(i) BF(i))
      PetscReal, Dimension(:), Pointer             :: F_Loc, RHS_Loc
      PetscInt, Dimension(:), Pointer              :: BCFlag_Loc
#if defined PB_2D
      Type (Vect2D)                                :: F_Elem
#elif defined PB_3D  
      Type (Vect3D)                                :: F_Elem    
#endif
      PetscInt                                     :: iE, iEloc, iBlkId, iErr
      PetscInt                                     :: NumDoFVect
      PetscInt                                     :: iDoF, iGauss
      PetscLogDouble                               :: flops       

      Call PetscLogEventBegin(AppCtx%LogInfo%RHSAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
      flops      = 0.0
      
      NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim

      Allocate(BCFlag_Loc(NumDoFVect))
      Allocate(F_Loc(NumDoFVect))
      Allocate(RHS_Loc(NumDoFVect))

      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
      Do_iEloc: Do iEloc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         RHS_Loc = 0.0_Kr
         Call SectionIntRestrictClosure(AppCtx%BCUFlag, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, BCFlag_Loc, iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(AppCtx%F, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, F_Loc, iErr); CHKERRQ(ierr)
         Do_iGauss: Do iGauss = 1, size(AppCtx%ElemVect(iE)%Gauss_C)
            F_Elem = 0.0_Kr
            Do iDoF = 1, NumDoFVect
               F_Elem = F_Elem + AppCtx%ElemVect(iE)%BF(iDoF, iGauss) * F_Loc(iDoF)
            End Do
            Do iDoF = 1, NumDoFVect
               If (BCFlag_Loc(iDoF) == 0) Then
                  ! RHS terms due to forces
                  RHS_Loc(iDoF) = RHS_Loc(iDoF) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * ( AppCtx%ElemVect(iE)%BF(iDoF, iGauss) .DotP. F_Elem ) 
               End If
            End Do
         End Do Do_iGauss
         Call SectionRealUpdateClosure(RHS_Sec, AppCtx%MeshTopology%Mesh, iE-1, RHS_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iEloc
      
      DeAllocate(BCFlag_Loc)
      DeAllocate(F_Loc)
      DeAllocate(RHS_Loc)
      
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
      Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyLocalU_Event, iErr); CHKERRQ(iErr)
   End Subroutine RHSAssemblyBlock_Force
   


   Subroutine Step_U(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      KSPConvergedReason                           :: reason
      PetscInt                                     :: KSPNumIter
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
      Type(Vec)                                    :: U_Vec, RHSU_Vec
      PetscReal                                    :: VMin, VMax
      PetscReal                                    :: rDum
      PetscInt                                     :: iDum
#if defined WITH_TAO
      TaoTerminateReason                           :: TaoReason
      PetscReal                                    :: TaoResidual
#endif
      If (AppCtx%VarFracSchemeParam%U_UseTao) Then
#if defined WITH_TAO
         Call TaoAppGetSolutionVec(AppCtx%taoappU, U_Vec, iErr); CHKERRQ(iErr)
         Call SectionRealToVec(AppCtx%U, AppCtx%ScatterVect, SCATTER_FORWARD, U_Vec, ierr); CHKERRQ(ierr)
         If (AppCtx%AppParam%verbose > 0) Then
            Write(IOBuffer, *) 'Calling TaoSolveApplication\n'
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         End If
         Call TaoSolveApplication(AppCtx%taoappU, AppCtx%taoU, iErr); CHKERRQ(iErr)
         Call TaoGetSolutionStatus(AppCtx%taoU, KSPNumIter, rDum, TaoResidual, iDum, iDum, TaoReason, iErr); CHKERRQ(iErr)
         If ( TaoReason > 0) Then
            Write(IOBuffer, 102) KSPNumiter, TAOReason
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Else
            Write(IOBuffer, 103) TaoReason
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         End If     
#endif      
      Else
         Call MeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%U, RHSU_Vec, iErr); CHKERRQ(iErr)
         Call MeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%U, U_Vec, iErr); CHKERRQ(iErr)
         Call SectionRealToVec(AppCtx%U, AppCtx%ScatterVect, SCATTER_FORWARD, U_Vec, ierr); CHKERRQ(ierr)
         If (AppCtx%AppParam%verbose > 0) Then
            Write(IOBuffer, *) 'Assembling the Matrix and RHS for the U-subproblem \n' 
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
         End If
         Call RHSU_Assembly(RHSU_Vec, AppCtx)
         If (AppCtx%AppParam%verbose > 2) Then
            Call VecView(RHSU_Vec, AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
         End If

         Call MatU_Assembly(AppCtx%KU, AppCtx)
         If (AppCtx%AppParam%verbose > 2) Then
            Call MatView(AppCtx%KU, AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
         End If
         
         If (AppCtx%AppParam%verbose > 0) Then
            Write(IOBuffer, *) 'Calling KSPSolve for the U-subproblem\n' 
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
         End If

         Call KSPSolve(AppCtx%KSPU, RHSU_Vec, U_Vec, iErr); CHKERRQ(iErr)
         Call SectionRealToVec(AppCtx%U, AppCtx%ScatterVect, SCATTER_REVERSE, U_Vec, ierr); CHKERRQ(ierr)
      
         Call KSPGetConvergedReason(AppCtx%KSPU, reason, iErr); CHKERRQ(iErr)
         If ( reason > 0) Then
            Call KSPGetIterationNumber(AppCtx%KSPU, KSPNumIter, iErr); CHKERRQ(iErr)
            Write(IOBuffer, 100) KSPNumIter, reason
         Else
            Write(IOBuffer, 101) reason
         End If


         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Call VecDestroy(RHSU_Vec, iErr); CHKERRQ(iErr)
         Call VecDestroy(U_Vec, iErr); CHKERRQ(iErr)
      End If
      CHKMEMQ
100 Format('     KSP for U converged in  ', I5, ' iterations. KSPConvergedReason is     ', I3, '\n')
101 Format('[ERROR] KSP for U diverged. KSPConvergedReason is ', I2, '\n')
102 Format('     TAO solver converged in ', I5, ' iterations. Tao termination reason is ', I3, '\n')
103 Format('[ERROR] TaoSolveApplication did not converge. ', I2, '\n')      
   End Subroutine Step_U   


#if defined PB_2D
End Module m_VarFracQS_U2D
#elif defined PB_3D
End Module m_VarFracQS_U3D
#endif
