#if defined PB_2D
Module m_Elast2D
#elif defined PB_3D
Module m_Elast3D
#endif

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscmeshdef.h"
#include "finclude/petscviewerdef.h"

   Use m_MEF90
   Use m_RuptStruct
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
#if defined PB_2D
      Type(Element2D_Elast), Dimension(:), Pointer  :: ElemU
#elif defined PB_3D
      Type(Element3D_Elast), Dimension(:), Pointer  :: ElemU
#endif
      Type(SectionReal)                            :: U
      Type(SectionReal)                            :: StressU
      Type(SectionReal)                            :: StrainU
      Type(SectionReal)                            :: F
      Type(SectionReal)                            :: Theta
      PetscReal                                    :: Load
      PetscInt                                     :: TimeStep
      PetscReal                                    :: BulkEnergy
      Type(VecScatter)                             :: ScatterU
      Type(SectionInt)                             :: BCFlagU
      Type(Mat)                                    :: KU
      Type(Vec)                                    :: RHSU
      Type(KSP)                                    :: KSPU
      Type(PC)                                     :: PCU
      Type(LogInfo_Type)                           :: LogInfo
#if defined PB_2D
      Type(MatProp2D_Type)                         :: MatProp      
#elif defined PB_3D
      Type(MatProp3D_Type)                         :: MatProp
#endif
      Type(AppParam_Type)                          :: AppParam
      Type(RuptSchemeParam_Type)                   :: RuptSchemeParam
   End Type AppCtx_Type
   
   
Contains

   Subroutine ElastInit(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx

      PetscInt                                     :: iErr, i
      Type(Mesh)                                   :: Tmp_Mesh
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer, filename
      PetscInt                                     :: NumComponents
      PetscTruth                                   :: HasPrefix

      
      Call MEF90_Initialize()
      Call PetscOptionsHasName(PETSC_NULL_CHARACTER, '-verbose', AppCtx%AppParam%Verbose, iErr)    
      Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p',     AppCtx%AppParam%prefix, HasPrefix, iErr); CHKERRQ(iErr)
      If (.NOT. HasPrefix) Then
         Call PetscPrintf(PETSC_COMM_WORLD, "No input file prefix given\n", iErr)
         Call MEF90_Finalize()
         STOP
      End If
      
      If (AppCtx%AppParam%verbose) Then
         Write(filename, 101) Trim(AppCtx%AppParam%prefix), MEF90_MyRank
         Call PetscViewerASCIIOpen(PETSC_COMM_SELF, filename, AppCtx%AppParam%MyLogViewer, iErr); CHKERRQ(iErr);   
         Write(IOBuffer, 102) MEF90_MyRank, Trim(filename)
         Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Call PetscSynchronizedFlush (PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
   
         Write(filename, 103) Trim(AppCtx%AppParam%prefix)
         Call PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr);   
         Write(IOBuffer, 104) Trim(filename)
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
   
101 Format(A, '-', I4.4, '.log')
102 Format('Output from processor ', I4.4, ' redirected to file ', A, '\n'c)
103 Format(A,'.log')
104 Format('Collective output redirected to file ', A, '\n'c)

      Call Write_EXO_Case(AppCtx%AppParam%prefix, '%0.4d', MEF90_NumProcs)
      AppCtx%EXO%Comm = PETSC_COMM_WORLD
      AppCtx%EXO%filename = Trim(AppCtx%AppParam%prefix)//'.gen'

         !!! Reading and distributing sequential mesh
      If (MEF90_NumProcs == 1) Then
         Call MeshCreateExodus(PETSC_COMM_WORLD, AppCtx%EXO%filename, AppCtx%MeshTopology%mesh, ierr); CHKERRQ(iErr)
      Else
         Call MeshCreateExodus(PETSC_COMM_WORLD, AppCtx%EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
         Call MeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, AppCtx%MeshTopology%mesh, ierr); CHKERRQ(iErr)
         Call MeshDestroy(Tmp_mesh, ierr); CHKERRQ(iErr)
      End If
   
      Call MeshTopologyReadEXO(AppCtx%MeshTopology, AppCtx%EXO)
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) "Done reading and partitioning the mesh\n"c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
   
      AppCtx%MyEXO%comm = PETSC_COMM_SELF
      AppCtx%MyEXO%exoid = AppCtx%EXO%exoid
      Write(AppCtx%MyEXO%filename, 99) trim(AppCtx%AppParam%prefix), MEF90_MyRank
 99  Format(A, '-', I4.4, '.gen')
   
      !!! Initializes the values and names of the properties and variables
      Call RuptEXOVariable_Init(AppCtx%MyEXO)
      Call RuptEXOProperty_Init(AppCtx%MyEXO, AppCtx%MeshTopology)   
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) "Done with RuptEXOVariable_Init and RuptEXOProperty_Init\n"c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If

      !!! Set the element type for each block so that we can call ElementInit
      Do i = 1, AppCtx%MeshTopology%num_elem_blks
         AppCtx%MeshTopology%elem_blk(i)%Elem_Type = AppCtx%MyEXO%EBProperty( Rupt_EBProp_Elem_Type )%Value( AppCtx%MeshTopology%elem_blk(i)%ID )
         Call Init_Elem_Blk_Type(AppCtx%MeshTopology%Elem_Blk(i), AppCtx%MeshTopology%num_dim)
      End Do
      Call ElementInit(AppCtx%MeshTopology, AppCtx%ElemU, AppCtx%RuptSchemeParam%IntegOrder)
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) "Done with ElementInit\n"c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If

      !!! Create the Sections for the variables
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, AppCtx%MeshTopology%Num_Dim, AppCtx%U, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, AppCtx%MeshTopology%Num_Dim, AppCtx%F, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 1, AppCtx%Theta, iErr); CHKERRQ(iErr)

      NumComponents = AppCtx%MeshTopology%Num_Dim * (AppCtx%MeshTopology%Num_Dim + 1) / 2
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, NumComponents, AppCtx%StrainU, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, NumComponents, AppCtx%StressU, iErr); CHKERRQ(iErr)

      !!! Create the Scatter, Vec, Mat, KSP, PC
      Call MeshCreateGlobalScatter(AppCtx%MeshTopology%mesh, AppCtx%U, AppCtx%ScatterU, iErr); CHKERRQ(iErr)
      Call MeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%U, AppCtx%RHSU, iErr); CHKERRQ(iErr)

      Call MeshSetMaxDof(AppCtx%MeshTopology%Mesh, 1, iErr); CHKERRQ(iErr) 
      Call MeshCreateMatrix(AppCtx%MeshTopology%mesh, AppCtx%U, MATMPIAIJ, AppCtx%KU, iErr); CHKERRQ(iErr)
      
      Call KSPCreate(PETSC_COMM_WORLD, AppCtx%KSPU, iErr); CHKERRQ(iErr)
      Call KSPSetOperators(AppCtx%KSPU, AppCtx%KU, AppCtx%KU, SAME_NONZERO_PATTERN, iErr); CHKERRQ(iErr)
      Call KSPSetType(AppCtx%KSPU, KSPCG, iErr); CHKERRQ(iErr)
      Call KSPSetInitialGuessNonzero(AppCtx%KSPU, PETSC_TRUE, iErr); CHKERRQ(iErr)
      Call KSPSetFromOptions(AppCtx%KSPU, iErr); CHKERRQ(iErr)

      Call KSPGetPC(AppCtx%KSPU, AppCtx%PCU, iErr); CHKERRQ(iErr)
      Call PCSetType(AppCtx%PCU, PCBJACOBI, iErr); CHKERRQ(iErr)

      !!! Create the Section for the BC
      Call MeshGetVertexSectionInt(AppCtx%MeshTopology%mesh, AppCtx%MeshTopology%Num_Dim, AppCtx%BCFlagU, iErr); CHKERRQ(iErr)
#if defined PB_2D
      Call EXOProperty_InitBCFlag2D(AppCtx%EXO, AppCtx%MeshTopology, AppCtx%BCFlagU)
#elif defined PB_3D
      Call EXOProperty_InitBCFlag3D(AppCtx%EXO, AppCtx%MeshTopology, AppCtx%BCFlagU)
#endif

      !!! Read U, F, and Temperature
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(Rupt_VertVar_DisplacementX)%Offset, AppCtx%TimeStep, AppCtx%U) 
      Call Read_EXO_Result_Cell(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(Rupt_CellVar_StrainXX)%Offset, AppCtx%TimeStep, AppCtx%StrainU) 
      Call Read_EXO_Result_Cell(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(Rupt_CellVar_StressXX)%Offset, AppCtx%TimeStep, AppCtx%StressU) 


      AppCtx%TimeStep = 1
      !!! BB
   End Subroutine ElastInit
!----------------------------------------------------------------------------------------!      
! Solve (CM)   
! No changes wrt VectPoisson
!----------------------------------------------------------------------------------------!      
   
   Subroutine Solve(AppCtx)
      !!! CM
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      KSPConvergedReason                           :: reason
      PetscInt                                     :: KSPNumIter
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
      Type(Vec)                                    :: U_Vec
      
      Call PetscLogStagePush(AppCtx%LogInfo%KSPSolve_Stage, iErr); CHKERRQ(iErr)
      Call MeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%U, U_Vec, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(AppCtx%U, AppCtx%Scatter, SCATTER_FORWARD, U_Vec, ierr); CHKERRQ(ierr)
      Call KSPSolve(AppCtx%KSP, AppCtx%RHS, U_Vec, iErr); CHKERRQ(iErr)
      !!! Solve and store the solution in AppCtx%RHS
      
      Call SectionRealToVec(AppCtx%U, AppCtx%Scatter, SCATTER_REVERSE, U_Vec, ierr); CHKERRQ(ierr)
      !!! Scatter the solution from (Vec) AppCtx%RHS to (SectionReal) AppCtx%U
      
      Call KSPGetIterationNumber(AppCtx%KSP, KSPNumIter, iErr); CHKERRQ(iErr)
      Call KSPGetConvergedReason(AppCtx%KSP, reason, iErr); CHKERRQ(iErr)
      Write(IOBuffer, 100) KSPNumIter, reason
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      
      Call VecDestroy(U_Vec, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
100 Format('KSP converged in ', I5, ' iterations. KSPConvergedReason is ', I2, '\n'c)
   End Subroutine Solve
!----------------------------------------------------------------------------------------!      
! MatAssembly (CM)  
!----------------------------------------------------------------------------------------!      
   
   Subroutine MatAssembly(AppCtx)
      !!! CM
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iBlk, iE, iELoc, iErr
      PetscReal, Dimension(:,:), Pointer           :: MatElem
      
      Call PetscLogStagePush(AppCtx%LogInfo%MatAssembly_Stage, iErr); CHKERRQ(iErr)
      
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
		 Allocate(MatElem(AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim))
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            !---------------------------------------------------
			Call MatAssemblyLocal(iE, AppCtx%Mat_Prop(iBlk), AppCtx, MatElem)
			! Is it ok? AppCtx%Mat_Prop(iBlk)
            !---------------------------------------------------
			Write(MEF90_MyRank+100, *) 'ELement ', iE
            Write(MEF90_MyRank+100, *)  MatElem
            Call assembleMatrix(AppCtx%K, AppCtx%MeshTopology%mesh, AppCtx%U, iE-1, MatElem, ADD_VALUES, iErr); CHKERRQ(iErr)
         End Do Do_Elem_iE
         DeAllocate(MatElem)
      End Do Do_Elem_iBlk
      Call MatAssemblyBegin(AppCtx%K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (AppCtx%K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)

      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine MatAssembly
!----------------------------------------------------------------------------------------!      
! MatAsemblyLocal (CM)  
!----------------------------------------------------------------------------------------!      
   
   Subroutine MatAssemblyLocal(iE, Mat_Prop, AppCtx, MatElem)
      !!! CM
      PetscInt                                     :: iE
#if defined PB_2D
      Type(MatProp2D_Type)                         :: Mat_Prop
#elif
      Type(MatProp3D_Type)                         :: Mat_Prop
#endif
      Type(AppCtx_Type)                            :: AppCtx
      PetscReal, Dimension(:,:), Pointer           :: MatElem      
      PetscInt                                     :: iE
   
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoF, NumGauss
      PetscInt, Dimension(:), Pointer              :: BCFlag
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      
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
			   !------------------------- 
				  MatElem(iDoF2, iDoF1) =  MatElem(iDoF2, iDoF1)+ AppCtx%Elem(iE)%Gauss_C(iGauss) * ((AppCtx%MatProp%Hookes_Law * AppCtx%Elem(iE)%GradS_BF(iDoF1, iGauss)) .DotP. AppCtx%Elem(iE)%GradS_BF(iDoF2, iGauss))				  
               !-------------------------
			      Call PetscLogFlops(AppCtx%MeshTopology%num_dim * (AppCtx%MeshTopology%num_dim-1) +1 , iErr);CHKERRQ(iErr)
                  !!! Is that right?
               End Do
            End If
         End Do
      End Do
      
      DeAllocate(BCFlag)
      Call PetscLogEventEnd(AppCtx%LogInfo%MatAssemblyLocal_Event, iErr); CHKERRQ(iErr)
   End Subroutine MatAssemblyLocal

!----------------------------------------------------------------------------------------!      
! RHSAssembly (CM)  
!----------------------------------------------------------------------------------------!      
   Subroutine RHSAssembly(AppCtx)
      !!! CM
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      Type(SectionReal)                            :: RHSSec
      PetscInt                                     :: iBlk, iE, iELoc
      PetscReal, Dimension(:), Pointer             :: RHSElem
      PetscInt                                     :: NumDoFPerVertex 
      
      NumDoFPerVertex = AppCtx%MeshTopology%Num_Dim
      
      !!! Hopefully one day we will use assemble Vector instead of going through a section
      Call PetscLogStagePush(AppCtx%LogInfo%RHSAssembly_Stage, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, NumDoFPerVertex, RHSSec, iErr); CHKERRQ(iErr)
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Allocate(RHSElem(AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim))
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            Call RHSAssemblyLocal(iE, AppCtx, RHSElem)
            Call MeshUpdateAddClosure(AppCtx%MeshTopology%Mesh, RHSSec, iE-1, RHSElem, iErr); CHKERRQ(iErr)
         End Do Do_Elem_iE
         DeAllocate(RHSElem)
      End Do Do_Elem_iBlk
      Call SectionRealComplete(RHSSec, iErr); CHKERRQ(iErr)
      !!! VERY important! This is the equivalent of a ghost update
      Call SectionRealToVec(RHSSec, AppCtx%Scatter, SCATTER_FORWARD, AppCtx%RHS, ierr); CHKERRQ(ierr)
      Call SectionRealDestroy(RHSSec, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
  End Subroutine RHSAssembly
   
!----------------------------------------------------------------------------------------!      
! RHSAssemblyLocal (CM)  
!----------------------------------------------------------------------------------------!      
   Subroutine RHSAssemblyLocal(iE, Mat_Prop, AppCtx, RHSElem)
      !!! CM
      PetscInt                                     :: iE
      Type(AppCtx_Type)                            :: AppCtx
      PetscReal, Dimension(:), Pointer             :: RHSElem
      PetscReal, Dimension(:), Pointer             :: F
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoF, NumGauss
      PetscInt, Dimension(:), Pointer              :: BCFlag
      PetscInt                                     :: iDoF1, iDoF2, iGauss
#if defined PB_2D
	  Type (Vect2D)             				   :: TmpRHS
#elif defined PB_3D  
	  Type (Vect3D)             				   :: TmpRHS    
#endif
      RHSElem  = 0.0_Kr
      NumDoF   = Size(AppCtx%Elem(iE)%BF,1)
      NumGauss = Size(AppCtx%Elem(iE)%BF,2)
      Allocate(BCFlag(NumDoF))
      Call MeshRestrictClosureInt(AppCtx%MeshTopology%mesh, AppCtx%BCFlag, iE-1, NumDoF, BCFlag, iErr); CHKERRQ(ierr)
      Allocate(F(NumDoF))
      Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%F, iE-1, NumDoF, F, iErr); CHKERRQ(ierr)
      Do iGauss = 1, NumGauss
         TmpRHS = 0.0_Kr
         Do iDoF2 = 1, NumDoF
	         TmpRHS = TmpRHS + AppCtx%Elem(iE)%BF(iDoF2, iGauss) * F(iDoF2)
             Call PetscLogFlops(2 , iErr);CHKERRQ(iErr)
         End Do
         Do iDoF1 = 1, NumDoF
            If (BCFlag(iDoF1) == 0) Then
               RHSElem(iDoF1) = RHSElem(iDoF1) + AppCtx%Elem(iE)%Gauss_C(iGauss) * ( AppCtx%Elem(iE)%BF(iDoF1, iGauss) .DotP. TmpRHS )
               Call PetscLogFlops(3 , iErr);CHKERRQ(iErr)
            End If
         End Do
      End Do
      DeAllocate(BCFlag)
      DeAllocate(F)
      Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyLocal_Event, iErr); CHKERRQ(iErr)
   End Subroutine RHSAssemblyLocal
   
!----------------------------------------------------------------------------------------!      
! ComputeEnergy (CM)  
!----------------------------------------------------------------------------------------!      
   Subroutine ComputeEnergy(AppCtx)
      !!! CM
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoF, NumGauss
      PetscReal, Dimension(:), Pointer             :: F, U
      PetscInt                                     :: iBlk, iELoc, iE
      PetscInt                                     :: iDoF, iGauss
#if defined PB_2D
	  Type(Mat2D)                                  :: Strain_Elem, Stress_Elem    
	  Type(Vect2D)								   :: F_Elem, U_Elem  
#elif defined PB_3D
	  Type(Mat3D)                                 :: Strain_Elem, Stress_Elem
	  Type(Vect3D)								   :: F_Elem, U_Elem  
#endif
      PetscReal                                    :: MyEnergy

      Call PetscLogStagePush(AppCtx%LogInfo%PostProc_Stage, iErr); CHKERRQ(iErr)
      Call PetscLogEventBegin(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      
      MyEnergy = 0.0_Kr
      
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            NumDoF   = Size(AppCtx%Elem(iE)%BF,1)
            NumGauss = Size(AppCtx%Elem(iE)%BF,2)
            Allocate(F(NumDoF))
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%F, iE-1, NumDoF, F, iErr); CHKERRQ(ierr)
            Allocate(U(NumDoF))
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%U, iE-1, NumDoF, U, iErr); CHKERRQ(ierr)
            Do iGauss = 1, NumGauss
               Strain_Elem = 0.0_Kr
               Stress_Elem = 0.0_Kr
               F_Elem      = 0.0_Kr
               U_Elem      = 0.0_Kr
               Do iDoF = 1, NumDoF
                  Strain_Elem = Strain_Elem + AppCtx%Elem(iE)%GradS_BF(iDoF1, iGauss) * U(iDoF) 
                  F_Elem = F_Elem + F(iDoF)* AppCtx%Elem(iE)%BF(iDoF, iGauss) 
                  U_Elem = U_Elem + U(iDoF)*AppCtx%Elem(iE)%BF(iDoF, iGauss) 
                  Call PetscLogFlops(4*AppCtx%MeshTopology%Num_Dim+4, iErr)
               End Do
			   Stress_Elem = Stress_Elem + (AppCtx%MatProp%Hookes_Law * Strain_Elem) 
               MyEnergy = MyEnergy + AppCtx%Elem(iE)%Gauss_C(iGauss) * ( (Strain_Elem .DotP. Strain_Elem) * 0.5_Kr - (F_Elem .DotP. U_Elem))
               Call PetscLogFlops(AppCtx%MeshTopology%Num_Dim+4, iErr)
            End Do
         End Do Do_Elem_iE
      End Do Do_Elem_iBlk
      Call PetscGlobalSum(MyEnergy, AppCtx%Energy, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
               
      Call PetscLogEventEnd  (AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine ComputeEnergy
   
!----------------------------------------------------------------------------------------!      
! ComputeStrainStress (CM)  
!----------------------------------------------------------------------------------------!      
   Subroutine ComputeStrainStress(AppCtx)
     Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
#if defined PB_2D
      Type(Mat2D)                                 :: Grad
#elif defined PB_3D
      Type(Mat3D)                                 :: Grad
#endif
      PetscReal                                    :: Vol
      PetscInt                                     :: NumDoF, NumGauss
      PetscReal, Dimension(:), Pointer             :: U
      PetscInt                                     :: iBlk, iELoc, iE
      PetscInt                                     :: iDoF, iGauss
	  PetscReal, Dimension(:), Pointer             :: Grad_Ptr
	  	  
      Call PetscLogEventBegin(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      Call PetscLogStagePush (AppCtx%LogInfo%PostProc_Stage, iErr); CHKERRQ(iErr)
      Allocate(Grad_Ptr(AppCtx%MeshTopology%Num_Dim**2))
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            NumDoF   = Size(AppCtx%Elem(iE)%BF,1)
            NumGauss = Size(AppCtx%Elem(iE)%BF,2)
            Allocate(U(NumDoF))
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%U, iE-1, NumDoF, U, iErr); CHKERRQ(ierr)
            Grad = 0.0_Kr
            Vol  = 0.0_Kr
            Do iGauss = 1, NumGauss
               Do iDoF = 1, NumDoF
                  Grad = Grad + AppCtx%Elem(iE)%Gauss_C(iGauss) * AppCtx%Elem(iE)%Der_BF(iDoF, iGauss) * U(iDoF)
                  Vol = Vol + AppCtx%Elem(iE)%Gauss_C(iGauss) * (AppCtx%Elem(iE)%BF(iDoF, iGauss).DotP.AppCtx%Elem(iE)%BF(iDoF, iGauss))**0.5_Kr
                  Call PetscLogFlops(3*AppCtx%MeshTopology%Num_Dim+2, iErr)
               End Do
            End Do
            Grad = Grad / Vol
            Call PetscLogFlops(AppCtx%MeshTopology%Num_Dim, iErr)
#if defined PB_2D
            Grad_Ptr = (/ Grad%XX, Grad%XY, Grad%YX, Grad%YY /)
#elif defined PB_3D
            Grad_Ptr = (/ Grad%XX, Grad%XY, Grad%XZ, Grad%YX, Grad%YY, Grad%YZ,Grad%ZX, Grad%ZY, Grad%ZZ /)
#endif
            Call MeshUpdateClosure(AppCtx%MeshTopology%Mesh, AppCtx%GradU, iE-1, Grad_Ptr, iErr)
         End Do Do_Elem_iE
      End Do Do_Elem_iBlk
      DeAllocate(Grad_Ptr)
      Call PetscLogEventEnd  (AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End ComputeStrainStress
!----------------------------------------------------------------------------------------!      
   
   Subroutine Save(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx

      Call Write_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(Rupt_VertVar_DisplacementX)%Offset, AppCtx%TimeStep, AppCtx%U) 
      Call Write_EXO_Result_Cell(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(Rupt_CellVar_StrainXX)%Offset, AppCtx%TimeStep, AppCtx%StrainU) 
      Call Write_EXO_Result_Cell(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(Rupt_CellVar_StressXX)%Offset, AppCtx%TimeStep, AppCtx%StressU) 
      Call Write_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(Rupt_GlobVar_Load)%Offset, AppCtx%TimeStep, AppCtx%BulkEnergy)
   End Subroutine Save
   
   Subroutine ElastFinalize(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx

      PetscInt                                     :: iErr
      Character(len=MEF90_MXSTRLEN)                :: filename
   
      Call SectionRealDestroy(AppCtx%U, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(AppCtx%F, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(AppCtx%Theta, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(AppCtx%StrainU, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(AppCtx%STressU, iErr); CHKERRQ(iErr)
      
      Call VecScatterDestroy(AppCtx%ScatterU, iErr); CHKERRQ(iErr)
      Call SectionIntDestroy(AppCtx%BCFlagU, iErr); CHKERRQ(iErr)
      Call MatDestroy(AppCtx%KU, iErr); CHKERRQ(iErr)
      Call VecDestroy(AppCtx%RHSU, iErr); CHKERRQ(iErr)
      Call KSPDestroy(AppCtx%KSPU, iErr); CHKERRQ(iErr)
      Call MeshDestroy(AppCtx%MeshTopology%Mesh, iErr); CHKERRQ(ierr)

      If (AppCtx%AppParam%verbose) Then
         Call PetscViewerFlush(AppCtx%AppParam%MyLogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerFlush(AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerDestroy(AppCtx%AppParam%MyLogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerDestroy(AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
      End If
      Write(filename, 103) Trim(AppCtx%AppParam%prefix)
      Call PetscLogPrintSummary(PETSC_COMM_WORLD, filename, iErr); CHKERRQ(iErr)
      Call MEF90_Finalize()
103 Format(A,'.log')
   End Subroutine ElastFinalize
   
   
#if defined PB_2D
End Module m_Elast2D
#elif defined PB_3D
End Module m_Elast3D
#endif
