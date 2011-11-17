#if defined PB_2D
Module m_Poisson2D
#elif defined PB_3D
Module m_Poisson3D
#endif

#include "finclude/petscdef.h"

   Use m_MEF90
   
   Implicit NONE   

   Type LogInfo_Type
      PetscLogStage                                :: IO_Stage
      PetscLogStage                                :: Setup_Stage 
      PetscLogStage                                :: MeshCreateExodus_Stage
      PetscLogStage                                :: MeshDistribute_Stage
      PetscLogStage                                :: MatAssembly_Stage    
      PetscLogStage                                :: RHSAssembly_Stage
      PetscLogStage                                :: KSPSolve_Stage
      PetscLogStage                                :: PostProc_Stage
      
      PetscLogEvent                                :: MatAssemblyBlock_Event
      PetscLogEvent                                :: RHSAssemblyBlock_Event
      PetscLogEvent                                :: PostProc_Event
   End Type LogInfo_Type

   Type AppParam_Type
      PetscBool                                    :: Restart
      PetscInt                                     :: Verbose
      PetscInt                                     :: TestCase
      Character(len=MEF90_MXSTRLEN)                :: prefix
      Type(PetscViewer)                            :: LogViewer, MyLogViewer
   End Type AppParam_Type


   Type Heat_AppCtx_Type
      Type (MeshTopology_Type)                     :: MeshTopology
      Type (EXO_Type)                              :: EXO, MyEXO
#if defined PB_2D
      Type(Element2D_Scal), Dimension(:), Pointer  :: Elem
#elif defined PB_3D
      Type(Element3D_Scal), Dimension(:), Pointer  :: Elem
#endif
      Type(SectionReal)                            :: GradU
      PetscReal                                    :: ElasticEnergy
      Type(Field)                                  :: U
      Type(Field)                                  :: UBC
      Type(Field)                                  :: F
      PetscReal                                    :: ExtForcesWork
      PetscReal                                    :: TotalEnergy
      Type(Flag)                                   :: BCFlag
      Type(Mat)                                    :: K
      Type(Mat)                                    :: M ! Mass matrix 
      Type(Mat)                                    :: Jac ! jacobian
      Type(Field)                                  :: RHS
      Type(KSP)                                    :: KSP
      Type(PC)                                     :: PC
      Type(LogInfo_Type)                           :: LogInfo
      Type(AppParam_Type)                          :: AppParam
   !For TS
      Type(TS)                                     :: TS
      PetscInt                                     :: maxsteps, NumSteps
      PetscReal                                    :: maxtime
      PetscInt                                     :: VertVar_Temperature 
      PetscReal, Dimension(:), Pointer             :: Diff, B_Mensi
   End Type Heat_AppCtx_Type
   
   
Contains

   
#undef __FUNCT__
#define __FUNCT__ "ExoFormat_SimplePoisson"
   Subroutine EXOFormat_SimplePoisson(AppCtx)
      Type(Heat_AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr
   
      Call EXPVP (AppCtx%MyEXO%exoid, 'g', 3, iErr)
      Call EXPVAN(AppCtx%MyEXO%exoid, 'g', 3, (/'Elastic Energy ', 'Ext Forces work', 'Total Energy   '/), iErr)
      Call EXPVP (AppCtx%MyEXO%exoid, 'n', 2, iErr)
      Call EXPVAN(AppCtx%MyEXO%exoid, 'n', 2, (/'U', 'F'/), iErr)
#if defined PB_2D
      Call EXPVP (AppCtx%MyEXO%exoid, 'e', 2, iErr)
      Call EXPVAN(AppCtx%MyEXO%exoid, 'e', 2, (/'Grad U_X', 'Grad U_Y'/), iErr)
#elif defined PB_3D
      Call EXPVP (AppCtx%MyEXO%exoid, 'e', 3, iErr)
      Call EXPVAN(AppCtx%MyEXO%exoid, 'e', 3, (/'Grad U_X', 'Grad U_Y', 'Grad U_Z'/), iErr)
#endif
      Call EXPTIM(AppCtx%MyEXO%exoid, 1, 1.0_Kr, iErr)
   End Subroutine EXOFormat_SimplePoisson
   
#undef __FUNCT__
#define __FUNCT__ "InitLog"
   Subroutine InitLog(AppCtx)
      Type(Heat_AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr
      
      Call PetscLogEventRegister('MatAssembly Block', 0, AppCtx%LogInfo%MatAssemblyBlock_Event, ierr); CHKERRQ(ierr)
      Call PetscLogEventRegister('RHSAssembly Block', 0, AppCtx%LogInfo%RHSAssemblyBlock_Event, ierr); CHKERRQ(ierr)
      Call PetscLogEventRegister('Post Processing',   0, AppCtx%LogInfo%PostProc_Event,         ierr); CHKERRQ(ierr)

      Call PetscLogStageRegister("MeshCreateExodus", AppCtx%LogInfo%MeshCreateExodus_Stage, iErr)
      Call PetscLogStageRegister("MeshDistribute",   AppCtx%LogInfo%MeshDistribute_Stage,   iErr)
      Call PetscLogStageRegister("IO Stage",         AppCtx%LogInfo%IO_Stage,               iErr)
      Call PetscLogStageRegister("Setup",            AppCtx%LogInfo%Setup_Stage,            iErr)
      Call PetscLogStageRegister("Mat Assembly",     AppCtx%LogInfo%MatAssembly_Stage,      iErr)
      Call PetscLogStageRegister("RHS Assembly",     AppCtx%LogInfo%RHSAssembly_Stage,      iErr)
      Call PetscLogStageRegister("KSP Solve",        AppCtx%LogInfo%KSPSolve_Stage,         iErr)
      Call PetscLogStageRegister("Post Proc",        AppCtx%LogInfo%PostProc_Stage,         iErr)
   End Subroutine InitLog
   
#undef __FUNCT__
#define __FUNCT__ "Solve"
   Subroutine Solve(AppCtx)
      Type(Heat_AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      KSPConvergedReason                           :: KSPreason
      PetscInt                                     :: KSPNumIter
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
      
      Call PetscLogStagePush(AppCtx%LogInfo%KSPSolve_Stage, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(AppCtx%U%Sec, AppCtx%U%Scatter, SCATTER_FORWARD, AppCtx%U%Vec, iErr); CHKERRQ(iErr)

      Call DMMeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%RHS, AppCtx%RHS%Vec, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(AppCtx%RHS%Sec, AppCtx%RHS%Scatter, SCATTER_FORWARD, AppCtx%RHS%Vec, iErr); CHKERRQ(iErr)

      Call KSPSolve(AppCtx%KSP, AppCtx%RHS%Vec, AppCtx%U%Vec, iErr); CHKERRQ(iErr)
         
      Call KSPGetIterationNumber(AppCtx%KSP, KSPNumIter, iErr); CHKERRQ(iErr)
      Call KSPGetConvergedReason(AppCtx%KSP, KSPreason, iErr); CHKERRQ(iErr)
      Write(IOBuffer, 100) KSPNumIter, KSPreason
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(AppCtx%U%Sec, AppCtx%U%Scatter, SCATTER_REVERSE, AppCtx%U%Vec, ierr); CHKERRQ(ierr)
 
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
100 Format('KSP converged in ', I5, ' iterations. KSPConvergedReason is ', I2, '\n')
   End Subroutine Solve
   
 
#undef __FUNCT__
#define __FUNCT__ "MatAssembly"
   Subroutine HeatMatAssembly(AppCtx, MeshTopology)   
      Type(Heat_AppCtx_Type)                            :: AppCtx
      Type (MeshTopology_Type)                     :: MeshTopology
      PetscInt                                     :: iBlk, iErr
      
      Call MatInsertVertexBoundaryValues(AppCtx%K, AppCtx%U, AppCtx%BCFlag, MeshTopology)
      Call MatAssemblyBegin(AppCtx%K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (AppCtx%K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)

!      Call PetscLogStagePush(AppCtx%LogInfo%MatAssembly_Stage, iErr); CHKERRQ(iErr)
      Do_iBlk: Do iBlk = 1, MeshTopology%Num_Elem_Blks
         Call MatAssemblyBlock(iBlk, AppCtx, MeshTopology)
      End Do Do_iBlk
      Call MatAssemblyBegin(AppCtx%K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (AppCtx%K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)

!      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine HeatMatAssembly
      
   Subroutine MatAssemblyBlock(iBlk, AppCtx, MeshTopology)
      Type(Heat_AppCtx_Type)                       :: AppCtx
      Type (MeshTopology_Type)                     :: MeshTopology
      PetscInt                                     :: iBlk
      
      PetscInt                                     :: iE, iELoc, iErr, i
      PetscReal, Dimension(:,:), Pointer           :: MatElem
      PetscInt, Dimension(:), Pointer              :: BCFlag
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      PetscInt                                     :: NumDoFScal
      PetscLogDouble                               :: flops = 0
      PetscReal                                    :: lDiff
      PetscReal, Dimension(:), Pointer             :: T_Loc
      PetscReal                                    :: T_Elem
!      Call PetscLogEventBegin(AppCtx%LogInfo%MatAssemblyBlock_Event, iErr); CHKERRQ(iErr)
     
      NumDoFScal = MeshTopology%Elem_Blk(iBlk)%Num_DoF
      Allocate(MatElem(MeshTopology%Elem_Blk(iBlk)%Num_DoF, MeshTopology%Elem_Blk(iBlk)%Num_DoF))
      Allocate(BCFlag(MeshTopology%Elem_Blk(iBlk)%Num_DoF))
      Allocate(T_Loc(NumDoFScal))

      Do_iELoc: Do iELoc = 1, MeshTopology%Elem_Blk(iBlk)%Num_Elems
         i = MeshTopology%Elem_Blk(iBlk)%ID
         iE = MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         MatElem = 0.0_Kr
         BCFlag = 0
         Call SectionIntRestrictClosure(AppCtx%BCFlag%Sec, MeshTopology%mesh, iE-1, MeshTopology%Elem_Blk(iBlk)%Num_DoF, BCFlag, iErr); CHKERRQ(ierr)
         Do iGauss = 1, size(AppCtx%Elem(iE)%Gauss_C)
            Select Case(AppCtx%AppParam%TestCase)
            Case(2)
               lDiff = AppCtx%Diff(i) 
            Case(3)
               ! Mensi Law D(C) = A exp (B*C) 
               ! 1988 Mensi-Acker-Attolou Mater Struct
               ! C : water content,  A and B material parameters
               Call SectionRealRestrictClosure(AppCtx%U%Sec, MeshTopology%mesh,  iE-1, NumDoFScal, T_Loc, iErr); CHKERRQ(ierr)
               T_Elem = 0.0_Kr
               Do iDoF1 = 1, NumDoFScal
                  T_Elem = T_Elem + AppCtx%Elem(iE)%BF(iDoF1, iGauss) * T_Loc(iDoF1)
               End DO
               lDiff = AppCtx%Diff(i)*exp(AppCtx%B_Mensi(i)*T_Elem) 
            Case Default
               ldiff =1 
            End Select
            Do iDoF1 = 1, MeshTopology%Elem_Blk(iBlk)%Num_DoF
               If (BCFlag(iDoF1) == 0) Then
                  Do iDoF2 = 1, MeshTopology%Elem_Blk(iBlk)%Num_DoF
                    ! MatElem(iDoF1, iDoF1) = 1./2.
                     MatElem(iDoF2, iDoF1) = MatElem(iDoF2, iDoF1) + lDiff * AppCtx%Elem(iE)%Gauss_C(iGauss) * ( AppCtx%Elem(iE)%Grad_BF(iDoF1, iGauss) .DotP. AppCtx%Elem(iE)%Grad_BF(iDoF2, iGauss) )
                     flops = flops + 1
                  End Do
               End If
            End Do
         End Do
         Call DMMeshAssembleMatrix(AppCtx%K, MeshTopology%mesh, AppCtx%U%Sec, iE-1, MatElem, ADD_VALUES, iErr); CHKERRQ(iErr)
         !Call DMMeshAssembleMatrix(AppCtx%K, MeshTopology%mesh, AppCtx%U%Sec, iE-1, MatElem, INSERT_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iELoc
   
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
      DeAllocate(BCFlag)
      DeAllocate(MatElem)
!      Call PetscLogEventEnd(AppCtx%LogInfo%MatAssemblyBlock_Event, iErr); CHKERRQ(iErr)
   End Subroutine MatAssemblyBlock


#undef __FUNCT__
#define __FUNCT__ "RHSAssembly"
   Subroutine RHSAssembly(AppCtx, MeshTopology, MyExo)
      Type(Heat_AppCtx_Type)                       :: AppCtx
      Type (MeshTopology_Type)                     :: MeshTopology
      Type (EXO_Type)                              :: MyEXO

      PetscInt                                     :: iErr
      PetscInt                                     :: iBlk
!      Call PetscLogStagePush(AppCtx%LogInfo%RHSAssembly_Stage, iErr); CHKERRQ(iErr)

      Call SectionRealZero(AppCtx%RHS%Sec, iErr); CHKERRQ(iErr)
      
      Do_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Call RHSAssemblyBlock(iBlk, AppCtx, MeshTopology)
      End Do Do_iBlk

      Call SectionRealComplete(AppCtx%RHS%Sec, iErr); CHKERRQ(iErr)

      !!! Set Dirichlet Boundary Values
!Suppose that loading is contant (replace second to last by AppCtx%Timestep otherwise)
      Call Read_EXO_Result_Vertex(MyEXO, MeshTopology, 1, 1, AppCtx%UBC)
      Call FieldInsertVertexBoundaryValues(AppCtx%RHS, AppCtx%UBC, AppCtx%BCFlag, MeshTopology)

      Call SectionRealToVec(AppCtx%RHS%Sec, AppCtx%RHS%Scatter, SCATTER_FORWARD, AppCtx%RHS%Vec, iErr); CHKERRQ(iErr)
      Call FieldInsertVertexBoundaryValues(AppCtx%RHS, AppCtx%UBC, AppCtx%BCFlag, MeshTopology)
      !!! VERY important! This is the equivalent of a ghost update
!      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine RHSAssembly

#undef __FUNCT__
#define __FUNCT__ "RHSAssemblyBlock"
   Subroutine RHSAssemblyBlock(iBlk, AppCtx, MeshTopology)
      PetscInt                                     :: iBlk
      Type(Heat_AppCtx_Type)                       :: AppCtx
      Type (MeshTopology_Type)                     :: MeshTopology

      PetscInt                                     :: iErr
      PetscInt                                     :: iBlkId
      PetscInt                                     :: Num_DoF, iDoF
      PetscInt                                     :: iEloc, iE, iGauss
      PetscInt, Dimension(:), Pointer              :: BCFlag_Loc
      PetscReal, Dimension(:), Pointer             :: RHS_Loc, F_Loc
      PetscReal                                    :: F_Elem
      PetscLogDouble                               :: flops 
      
      Call PetscLogEventBegin(AppCtx%LogInfo%RHSAssemblyBlock_Event, iErr); CHKERRQ(iErr)
      flops = 0.0

      Num_DoF = MeshTopology%Elem_Blk(iBlk)%Num_DoF
      Allocate(F_Loc(Num_DoF))
      Allocate(RHS_Loc(Num_DoF))
      Allocate(BCFlag_Loc(Num_DoF))

      iBlkID = MeshTopology%Elem_Blk(iBlk)%ID
      Do_iEloc: Do iELoc = 1, MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         RHS_Loc = 0.0_Kr
         Call SectionRealRestrictClosure(AppCtx%F%Sec, MeshTopology%mesh, iE-1, Num_DoF, F_Loc, iErr); CHKERRQ(ierr)
         Call SectionIntRestrictClosure(AppCtx%BCFlag%Sec, MeshTopology%mesh, iE-1, Num_DoF, BCFlag_Loc, iErr); CHKERRQ(ierr)
         Do iGauss = 1, Size(AppCtx%Elem(iE)%Gauss_C)
            F_Elem = 0.0_Kr
            Do iDoF = 1, Num_DoF
               F_Elem = F_Elem + AppCtx%Elem(iE)%BF(iDoF, iGauss) * F_Loc(iDoF)
               flops = flops + 2.0
            End Do
            Do iDoF = 1, Num_DoF
               If (BCFlag_Loc(iDoF) == 0) Then
                  RHS_Loc(iDoF) = RHS_Loc(iDoF) + AppCtx%Elem(iE)%Gauss_C(iGauss) * ( F_Elem * AppCtx%Elem(iE)%BF(iDoF, iGauss) )
                  flops = flops + 3.0
               End If
            End Do
         End Do
         Call SectionRealUpdateClosure(AppCtx%RHS%Sec, MeshTopology%Mesh, iE-1, RHS_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iEloc

      DeAllocate(BCFlag_Loc)
      DeAllocate(RHS_Loc)
      DeAllocate(F_Loc)
      Call PetscLogFlops(flops , iErr);CHKERRQ(iErr)
      Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyBlock_Event, iErr); CHKERRQ(iErr)
   End Subroutine RHSAssemblyBlock

#undef __FUNCT__
#define __FUNCT__ "ComputeEnergy"
   Subroutine ComputeEnergy(AppCtx)
      Type(Heat_AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoF, NumGauss
      PetscReal, Dimension(:), Pointer             :: F, U
      PetscInt                                     :: iBlk, iELoc, iE
      PetscInt                                     :: iDoF, iGauss
#if defined PB_2D
      Type(Vect2D)                                 :: Strain_Elem, Stress_Elem      
#elif defined PB_3D
      Type(Vect3D)                                 :: Strain_Elem, Stress_Elem      
#endif
      PetscReal                                    :: F_Elem, U_Elem
      PetscReal                                    :: MyElasticEnergy, MyExtForcesWork
      PetscLogDouble                               :: flops

      Call PetscLogStagePush(AppCtx%LogInfo%PostProc_Stage, iErr); CHKERRQ(iErr)
      Call PetscLogEventBegin(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      
      MyElasticEnergy = 0.0_Kr
      MyExtForcesWork = 0.0_Kr
      flops = 0.0_Kr
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            NumDoF   = Size(AppCtx%Elem(iE)%BF,1)
            NumGauss = Size(AppCtx%Elem(iE)%BF,2)
            Allocate(F(NumDoF))
            Call SectionRealRestrictClosure(AppCtx%F%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoF, F, iErr); CHKERRQ(ierr)
            Allocate(U(NumDoF))
            Call SectionRealRestrictClosure(AppCtx%U%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoF, U, iErr); CHKERRQ(ierr)
            Do iGauss = 1, NumGauss
               Strain_Elem = 0.0_Kr
               Stress_Elem = 0.0_Kr
               F_Elem      = 0.0_Kr
               U_Elem      = 0.0_Kr
               Do iDoF = 1, NumDoF
                  Stress_Elem = Stress_Elem + AppCtx%Elem(iE)%Grad_BF(iDoF, iGauss) * U(iDoF)
                  Strain_Elem = Strain_Elem + AppCtx%Elem(iE)%Grad_BF(iDoF, iGauss) * U(iDoF)
                  F_Elem = F_Elem + AppCtx%Elem(iE)%BF(iDoF, iGauss) * F(iDoF)
                  U_Elem = U_Elem + AppCtx%Elem(iE)%BF(iDoF, iGauss) * U(iDoF)
                  flops = flops + 4
               End Do
               MyElasticEnergy = MyElasticEnergy + AppCtx%Elem(iE)%Gauss_C(iGauss) * ( (Strain_Elem .DotP. Strain_Elem) * 0.5_Kr)
               MyExtForcesWork = MyExtForcesWork - AppCtx%Elem(iE)%Gauss_C(iGauss) * F_Elem * U_Elem
               flops = flops + 5
            End Do
            DeAllocate(U)
            DeAllocate(F)
         End Do Do_Elem_iE
      End Do Do_Elem_iBlk
      Call MPI_AllReduce(MyElasticEnergy, AppCtx%ElasticEnergy, 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      Call MPI_AllReduce(MyExtForcesWork, AppCtx%ExtForcesWork, 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      AppCtx%TotalEnergy = AppCtx%ElasticEnergy + AppCtx%ExtForcesWork
      Call PetscLogFlops(flops, iErr)
      Call PetscLogEventEnd  (AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine ComputeEnergy
 
#undef __FUNCT__
#define __FUNCT__ "ComputeGradU"
   Subroutine ComputeGradU(AppCtx)
      Type(Heat_AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
#if defined PB_2D
      Type(Vect2D)                                 :: Grad
#elif defined PB_3D
      Type(Vect3D)                                 :: Grad
#endif
      PetscReal                                    :: Vol
      PetscInt                                     :: NumDoF, NumGauss
      PetscReal, Dimension(:), Pointer             :: U
      PetscInt                                     :: iBlk, iELoc, iE
      PetscInt                                     :: iDoF, iGauss
      PetscReal, Dimension(:), Pointer             :: Grad_Ptr
      PetscLogDouble                               :: flops = 0

      Call PetscLogStagePush (AppCtx%LogInfo%PostProc_Stage, iErr); CHKERRQ(iErr)
      Call PetscLogEventBegin(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      Call DMMeshGetCellSectionReal(AppCtx%MeshTopology%mesh,   'GradU', AppCtx%MeshTopology%Num_Dim, AppCtx%GradU, iErr); CHKERRQ(iErr)

      Allocate(Grad_Ptr(AppCtx%MeshTopology%Num_Dim))
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            NumDoF   = Size(AppCtx%Elem(iE)%BF,1)
            NumGauss = Size(AppCtx%Elem(iE)%BF,2)
            Allocate(U(NumDoF))
            Call SectionRealRestrictClosure(AppCtx%U%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoF, U, iErr); CHKERRQ(ierr)
            Grad = 0.0_Kr
            Vol  = 0.0_Kr
            Do iGauss = 1, NumGauss
               Do iDoF = 1, NumDoF
                  Grad = Grad + (AppCtx%Elem(iE)%Gauss_C(iGauss) * U(iDoF) ) * AppCtx%Elem(iE)%Grad_BF(iDoF, iGauss) 
                  Vol = Vol + AppCtx%Elem(iE)%Gauss_C(iGauss) * AppCtx%Elem(iE)%BF(iDoF, iGauss)
                  flops = flops + 3
               End Do
            End Do
            Grad = Grad / Vol
#if defined PB_2D
            Grad_Ptr = (/ Grad%X, Grad%Y /)
#elif defined PB_3D
            Grad_Ptr = (/ Grad%X, Grad%Y, Grad%Z /)
#endif
            Call SectionRealUpdateClosure(AppCtx%GradU, AppCtx%MeshTopology%Mesh, iE-1, Grad_Ptr, INSERT_VALUES, iErr)
            DeAllocate(U)
         End Do Do_Elem_iE
      End Do Do_Elem_iBlk
      DeAllocate(Grad_Ptr)
      Call PetscLogFlops(flops, iErr)
      Call PetscLogEventEnd  (AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine ComputeGradU
 
#undef __FUNCT__
#define __FUNCT__ "SimplePoissonFinalize"
   Subroutine SimplePoissonFinalize(AppCtx)   
      Type(Heat_AppCtx_Type)                            :: AppCtx

      PetscInt                                     :: iErr
      Character(len=MEF90_MXSTRLEN)                :: filename
      Type(PetscViewer)                            :: LogViewer

      Call SectionRealDestroy(AppCtx%GradU, iErr); CHKERRQ(iErr)
      Call FieldDestroy(AppCtx%U)
      !Call FieldDestroy(AppCtx%GradU)
      Call FieldDestroy(AppCtx%F)
      Call FieldDestroy(AppCtx%RHS)

      Call SectionIntDestroy(AppCtx%BCFlag, iErr); CHKERRQ(iErr)
      Call MatDestroy(AppCtx%K, iErr); CHKERRQ(iErr)
      Call KSPDestroy(AppCtx%KSP, iErr); CHKERRQ(iErr)
      Call DMDestroy(AppCtx%MeshTopology%Mesh, iErr); CHKERRQ(ierr)
      
      If (AppCtx%AppParam%verbose > 1) Then
         Call PetscViewerFlush(AppCtx%AppParam%MyLogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerFlush(AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerDestroy(AppCtx%AppParam%MyLogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerDestroy(AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
      End If
      Write(filename, 103) Trim(AppCtx%AppParam%prefix)
103 Format(A,'-logsummary.txt')
      Call PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, LogViewer, ierr);CHKERRQ(ierr)
      Call PetscViewerASCIIPrintf(LogViewer,filename,ierr);CHKERRQ(ierr)
      Call PetscLogView(LogViewer,ierr);CHKERRQ(ierr)
      Call PetscViewerDestroy(LogViewer,ierr);CHKERRQ(ierr)
!       Call PetscLogPrintSummary(PETSC_COMM_WORLD, filename, iErr); CHKERRQ(iErr)
      Call MEF90_Finalize()
   End Subroutine SimplePoissonFinalize

#if defined PB_2D
End Module m_Poisson2D
#elif defined PB_3D
End Module m_Poisson3D
#endif
