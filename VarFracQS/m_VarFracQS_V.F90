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
#elif defined PB_3D
   Use m_VarFracQS_Types3D
#endif   
   Use m_MEF90
   Use m_VarFrac_Struct
   Use petsc
   Use petscvec
   Use petscmat
   Use petscksp
   Use petscmesh

   Implicit NONE   

Contains
   Subroutine Init_TS_V(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      Type(SectionReal)                            :: VBC, VBT
      PetscReal, Dimension(:), Pointer             :: V_Ptr, VBC_Ptr
      PetscInt, Dimension(:), Pointer              :: BCVFlag, IrrevFlag
      PetscInt                                     :: i, iErr
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer      
      PetscReal                                    :: MyIrrevEQ_Counter  = 0.0
      PetscReal                                    :: IrrevEQ_Counter
      
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'VBC', 1, VBC, iErr); CHKERRQ(iErr)
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Fracture)%Offset, AppCtx%TimeStep, VBC)
            
      !!! Update the BC Flag from the current V
      Select Case(AppCtx%VarFracSchemeParam%IrrevType)
      Case(VarFrac_Irrev_Eq)
         If (AppCtx%AppParam%verbose > 0) Then
            Write(IOBuffer, *) "Irreversibility with Irrev_EQ\n"c
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         End If
         !!! Update the IrrevFlag Section
         !!! If I'd understand how to use SectionGetArrayF90 I'd use it
         Allocate(IrrevFlag(1))
         IrrevFlag = VarFrac_BC_Type_DIRI
         If (AppCtx%IsBT) Then
            If (AppCtx%AppParam%verbose > 0) Then
               Write(IOBuffer, *) "Backtracking, so reading timestep if TS>1\n"c
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
               Write(IOBuffer, *) "Not Backtracking, so getting blocked nodes from V\n"c
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
            Write(IOBuffer, *) "Number of blocked nodes for V: ", IrrevEQ_Counter, "\n"c
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         End If      
      Case(VarFrac_Irrev_Ineq)   
         SETERRQ(PETSC_ERR_SUP, 'NotImplemented yet\n', iErr)
      End Select
      
      Select Case(AppCtx%VarFracSchemeParam%InitV)
      Case(VarFrac_INIT_V_PREV)
         If (AppCtx%AppParam%verbose > 0) Then
            Write(IOBuffer, *) "Initializing V with ", AppCtx%VarFracSchemeParam%InitV, "\n"c
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         End If      
         Allocate(V_Ptr(1))
         Do i = 1, AppCtx%MeshTopology%Num_Verts       
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
   End Subroutine Init_TS_V

   
!----------------------------------------------------------------------------------------!      
! MatAssembly V (CM)  
!----------------------------------------------------------------------------------------!      
   Subroutine MatV_Assembly(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iBlk, iBlkID, iErr
      
      Call PetscLogStagePush(AppCtx%LogInfo%MatAssemblyV_Stage, iErr); CHKERRQ(iErr)

      Call MatZeroEntries(AppCtx%KV, iErr); CHKERRQ(iErr)
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
         Select Case (AppCtx%VarFracSchemeParam%AtNum)
         !!! Case(1)
         !!!    Call MatV_AssemblyBlk_AT1(iBlk, AppCtx%KV, .TRUE., AppCtx)
         Case(2)
            If (AppCtx%MyEXO%EBProperty(VarFrac_EBProp_IsBrittle)%Value(iBlkID) /= 0) Then
               Call MatV_AssemblyBlk_Brittle_AT2(iBlk, AppCtx%KV, .TRUE., AppCtx)
            Else
               Call MatV_AssemblyBlk_NonBrittle_AT2(iBlk, AppCtx%KV, .TRUE., AppCtx)
            End If
         Case Default
            SETERRQ(PETSC_ERR_SUP, 'Only AT1 and AT2 are implemented\n', iErr)
      End Select
      End Do Do_Elem_iBlk
      Call MatAssemblyBegin(AppCtx%KV, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (AppCtx%KV, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)

      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine MatV_Assembly

#if defined WITH_TAO
   Subroutine FormHessianV(tao, X, H, Hpre, flg, AppCtx, iErr)
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
               Call MatV_AssemblyBlk_Brittle_AT1(iBlk, AppCtx%KV, .FALSE., AppCtx)
            Else
               Call MatV_AssemblyBlk_NonBrittle_AT1(iBlk, AppCtx%KV, .FALSE., AppCtx)
            End If
         Case(2)
            If (AppCtx%MyEXO%EBProperty(VarFrac_EBProp_IsBrittle)%Value(iBlkID) /= 0) Then
               Call MatV_AssemblyBlk_Brittle_AT2(iBlk, AppCtx%KV, .FALSE., AppCtx)
            Else
               Call MatV_AssemblyBlk_NonBrittle_AT2(iBlk, AppCtx%KV, .FALSE., AppCtx)
            End If
         Case Default
            SETERRQ(PETSC_ERR_SUP, 'Only AT1 and AT2 are implemented\n', iErr)
      End Do Do_Elem_iBlk
      Call MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (H, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      If (AppCtx%AppParam%verbose > 1) Then
         Call MatView(H, PETSC_VIEWER_STDOUT_WORLD, iErr)
      End If

      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine FormHessianV
#endif

   Subroutine MatV_AssemblyBlk_Brittle_AT2(iBlk, H, DoBC, AppCtx)
      PetscInt                                     :: iBlk
      Type(Mat)                                    :: H
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
      PetscReal                                    :: C2_V, C2_GradV
#if defined PB_2D
      Type(MatS2D)                                 :: Strain_Elem, Effective_Strain_Elem
#elif defined PB_3D
      Type(MatS3D)                                 :: Strain_Elem, Effective_Strain_Elem
#endif      
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
            Effective_Strain_Elem  =  Strain_Elem - (Theta_Elem * AppCtx%MatProp(iBlkID)%Therm_Exp)   
            !! Calculate the coefficients of the terms v^2 (C2_V) et GradV*GradV (C2_GradV) of the energy functional
            C2_V = 2.0_Kr * AppCtx%VarFracSchemeParam%ATCv * AppCtx%MatProp(iBlkID)%Toughness / AppCtx%VarFracSchemeParam%Epsilon  + ((AppCtx%MatProp(iBlkID)%Hookes_Law * Effective_Strain_Elem) .DotP. Effective_Strain_Elem)
            flops = flops + 3.0
            C2_GradV = 2.0_Kr * AppCtx%VarFracSchemeParam%ATCv * AppCtx%MatProp(iBlkID)%Toughness * AppCtx%VarFracSchemeParam%Epsilon 
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
      DeAllocate(U)      
      DeAllocate(Theta)
      DeAllocate(BCFlag)
      DeAllocate(IrrevFlag)
      Call PetscLogEventEnd(AppCtx%LogInfo%MatAssemblyLocalV_Event, iErr); CHKERRQ(iErr)
   End Subroutine MatV_AssemblyBlk_Brittle_AT2

   Subroutine MatV_AssemblyBlk_NonBrittle_AT2(iBlk, H, DoBC, AppCtx)
      PetscInt                                     :: iBlk
      Type(Mat)                                    :: H
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
            C2_V = 2.0_Kr * AppCtx%VarFracSchemeParam%ATCv * AppCtx%MatProp(iBlkID)%Toughness / AppCtx%VarFracSchemeParam%Epsilon  
            flops = flops + 2.0
            C2_GradV = 2.0_Kr * AppCtx%VarFracSchemeParam%ATCv * AppCtx%MatProp(iBlkID)%Toughness * AppCtx%VarFracSchemeParam%Epsilon 
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
   End Subroutine MatV_AssemblyBlk_NonBrittle_AT2

   Subroutine MatV_AssemblyBlk_Brittle_AT1(iBlk, H, DoBC, AppCtx)
      PetscInt                                     :: iBlk
      Type(Mat)                                    :: H
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
      PetscReal                                    :: C2_V, C2_GradV
#if defined PB_2D
      Type(MatS2D)                                 :: Strain_Elem, Effective_Strain_Elem
#elif defined PB_3D
      Type(MatS3D)                                 :: Strain_Elem, Effective_Strain_Elem
#endif      
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
            Effective_Strain_Elem  =  Strain_Elem - (Theta_Elem * AppCtx%MatProp(iBlkID)%Therm_Exp)   
            !! Calculate the coefficients of the terms v^2 (C2_V) et GradV*GradV (C2_GradV) of the energy functional
            C2_V = (AppCtx%MatProp(iBlkID)%Hookes_Law * Effective_Strain_Elem) .DotP. Effective_Strain_Elem
            C2_GradV = 2.0_Kr * AppCtx%VarFracSchemeParam%ATCv * AppCtx%MatProp(iBlkID)%Toughness * AppCtx%VarFracSchemeParam%Epsilon 
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
      DeAllocate(U)      
      DeAllocate(Theta)
      DeAllocate(BCFlag)
      DeAllocate(IrrevFlag)
      Call PetscLogEventEnd(AppCtx%LogInfo%MatAssemblyLocalV_Event, iErr); CHKERRQ(iErr)
   End Subroutine MatV_AssemblyBlk_Brittle_AT1

   Subroutine MatV_AssemblyBlk_NonBrittle_AT1(iBlk, H, DoBC, AppCtx)
      PetscInt                                     :: iBlk
      Type(Mat)                                    :: H
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
            C2_GradV = 2.0_Kr * AppCtx%VarFracSchemeParam%ATCv * AppCtx%MatProp(iBlkID)%Toughness * AppCtx%VarFracSchemeParam%Epsilon 
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
   End Subroutine MatV_AssemblyBlk_NonBrittle_AT1

!----------------------------------------------------------------------------------------!      
! RHSAssembly (CM)  
!----------------------------------------------------------------------------------------!      
   Subroutine RHSV_Assembly(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      PetscInt                                     :: iBlk
            
      !!! Hopefully one day we will use assemble Vector instead of going through a section
      Call PetscLogStagePush(AppCtx%LogInfo%RHSAssemblyV_Stage, iErr); CHKERRQ(iErr)
      
      Call SectionRealZero(AppCtx%RHSV, iErr); CHKERRQ(iErr)
      Do_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Select Case (AppCtx%VarFracSchemeParam%AtNum)
         Case(2)
            Call RHSV_AssemblyBlock_AT2(iBlk, AppCtx)
         Case Default
            SETERRQ(PETSC_ERR_SUP, 'Only AT1 and AT2 are implemented\n', iErr)
      End Select
      End Do Do_iBlk

      Call SectionRealComplete(AppCtx%RHSV, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
  End Subroutine RHSV_Assembly
   
!----------------------------------------------------------------------------------------!      
! RHSAssemblyLocal (CM)  
!----------------------------------------------------------------------------------------!      
   Subroutine RHSV_AssemblyBlock_AT2(iBlk, AppCtx)
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
      C1_V =  2.0_Kr * AppCtx%VarFracSchemeParam%ATCv * AppCtx%MatProp(iBlkId)%Toughness / AppCtx%VarFracSchemeParam%Epsilon  
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
         Call SectionRealUpdateClosure(AppCtx%RHSV, AppCtx%MeshTopology%Mesh, iE-1, RHS_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iEloc

      DeAllocate(BCFlag_Loc)
      DeAllocate(IrrevFlag_Loc)      
      DeAllocate(RHS_Loc)
      Call PetscLogFlops(flops, iErr); CHKERRQ(iErr)
      Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyLocalV_Event, iErr); CHKERRQ(iErr)
   End Subroutine RHSV_AssemblyBlock_AT2
      

   !----------------------------------------------------------------------------------------!      
   ! Solve V (CM)   
   !----------------------------------------------------------------------------------------!      
   
   Subroutine Solve_V(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      KSPConvergedReason                           :: reason
      PetscInt                                     :: KSPNumIter
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
      Type(Vec)                                    :: V_Vec, V_Old, RHSV_Vec
      PetscReal                                    :: VMin, VMax
      
      Call PetscLogStagePush(AppCtx%LogInfo%KSPSolveV_Stage, iErr); CHKERRQ(iErr)
  
      Call MeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%V, V_Vec, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(AppCtx%V, AppCtx%ScatterScal, SCATTER_FORWARD, V_Vec, ierr); CHKERRQ(ierr)

      Call MeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%RHSV, RHSV_Vec, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(AppCtx%RHSV, AppCtx%ScatterScal, SCATTER_FORWARD, RHSV_Vec, ierr); CHKERRQ(ierr)

      Call VecDuplicate(V_Vec, V_Old, iErr); CHKERRQ(iErr)
      Call VecCopy(V_Vec, V_Old, iErr); CHKERRQ(iErr)

      Call KSPSolve(AppCtx%KSPV, RHSV_Vec, V_Vec, iErr); CHKERRQ(iErr)
      !!! Solve and store the solution in AppCtx%RHS
      
      Call SectionRealToVec(AppCtx%V, AppCtx%ScatterScal, SCATTER_REVERSE, V_Vec, ierr); CHKERRQ(ierr)
       
      !!! Scatter the solution from (Vec) AppCtx%RHS to (SectionReal) AppCtx%V
      
      Call KSPGetConvergedReason(AppCtx%KSPV, reason, iErr); CHKERRQ(iErr)
      If ( reason > 0) Then
         Call KSPGetIterationNumber(AppCtx%KSPV, KSPNumIter, iErr); CHKERRQ(iErr)
         Write(IOBuffer, 100) KSPNumIter, reason
      Else
         Write(IOBuffer, 101) reason
      End If
      
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

      Call VecMin(V_Vec, PETSC_NULL_INTEGER, VMin, iErr); CHKERRQ(iErr)
      Call VecMax(V_Vec, PETSC_NULL_INTEGER, VMax, iErr); CHKERRQ(iErr)
      Write(IOBuffer, 700) VMin, VMax
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

       Call VecAxPy(V_Old, -1.0_Kr, V_Vec, iErr)
       Call VecNorm(V_Old, NORM_INFINITY, AppCtx%ErrV, iErr)

      Write(IOBuffer, 800) AppCtx%ErrV
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

      
      Call VecDestroy(V_Old, iErr); CHKERRQ(iErr)
      Call VecDestroy(V_Vec, iErr); CHKERRQ(iErr)
      Call VecDestroy(RHSV_Vec, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
100 Format('     KSP for V converged in ', I5, ' iterations. KSPConvergedReason is ', I3, '\n')
101 Format('[ERROR] KSP for V diverged. KSPConvergedReason is ', I2, '\n')
700 Format('     VMin / Max:   ', T24, 2(ES12.5, '  '), '\n')
800 Format('     Max change V: ', T24, ES12.5, '\n')
   End Subroutine Solve_V

   
   
#if defined PB_2D
End Module m_VarFracQS_V2D
#elif defined PB_3D
End Module m_VarFracQS_V3D
#endif
