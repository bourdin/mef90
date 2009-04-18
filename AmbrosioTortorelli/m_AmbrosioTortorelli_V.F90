#if defined PB_2D
Module m_AmbrosioTortorelli_U2D
#elif defined PB_3D
Module m_AmbrosioTortorelli_U3D
#endif

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscmeshdef.h"
#include "finclude/petscviewerdef.h"

#if defined PB_2D
   Use m_AmbrosioTortorelli_Types2D
#elif defined PB_3D
   Use m_AmbrosioTortorelli_Types3D
#endif   Use m_AmbrosioTortorelli_Types
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
      Call RuptSchemeParam_GetFromOptions(AppCtx%RuptSchemeParam)
      
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
      Call EXOProperty_Read(AppCtx%MyEXO)   
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) "Done with RuptEXOVariable_Init and RuptEXOProperty_Read\n"c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Call MeshTopologyView(AppCtx%MeshTopology, AppCtx%AppParam%MyLogViewer)
         Call EXOView(AppCtx%MyEXO, AppCtx%AppParam%MyLogViewer)
      End If
      
      !!! Read Mat Properties from the CST file
      Call MatProp_Read(AppCtx%MeshTopology, AppCtx%MatProp, trim(AppCtx%AppParam%prefix)//'.CST')
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) "Done with MatProp_Read\n"c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
!         Call EXOView(AppCtx%MyEXO, AppCtx%AppParam%MyLogViewer)
      End If

      !!! Set the element type for each block so that we can call ElementInit
      Do i = 1, AppCtx%MeshTopology%num_elem_blks
         AppCtx%MeshTopology%elem_blk(i)%Elem_Type = AppCtx%MyEXO%EBProperty( Rupt_EBProp_Elem_Type )%Value( AppCtx%MeshTopology%elem_blk(i)%ID )
         Call Init_Elem_Blk_Type(AppCtx%MeshTopology%Elem_Blk(i), AppCtx%MeshTopology%num_dim)
      End Do
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) "Done with Init_Elem_Blk_Type\n"c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      
      
      Call ElementInit(AppCtx%MeshTopology, AppCtx%ElemVect, AppCtx%RuptSchemeParam%IntegOrder)
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) "Done with ElementInit Vect\n"c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      Call ElementInit(AppCtx%MeshTopology, AppCtx%ElemScal, AppCtx%RuptSchemeParam%IntegOrder)
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) "Done with ElementInit Scal\n"c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If

      !!! Create the Sections for the variables
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, AppCtx%MeshTopology%Num_Dim, AppCtx%U, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, AppCtx%MeshTopology%Num_Dim, AppCtx%F, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 1, AppCtx%Theta, iErr); CHKERRQ(iErr)

      NumComponents = AppCtx%MeshTopology%Num_Dim * (AppCtx%MeshTopology%Num_Dim + 1) / 2
      Call MeshGetCellSectionReal(AppCtx%MeshTopology%mesh, NumComponents, AppCtx%StrainU, iErr); CHKERRQ(iErr)
      Call MeshGetCellSectionReal(AppCtx%MeshTopology%mesh, NumComponents, AppCtx%StressU, iErr); CHKERRQ(iErr)

      !!! Create the Scatter, Vec, Mat, KSP, PC
      Call MeshCreateGlobalScatter(AppCtx%MeshTopology%mesh, AppCtx%U, AppCtx%ScatterVect, iErr); CHKERRQ(iErr)
      Call MeshCreateGlobalScatter(AppCtx%MeshTopology%mesh, AppCtx%Theta, AppCtx%ScatterScal, iErr); CHKERRQ(iErr)
      Call MeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%U, AppCtx%RHSU, iErr); CHKERRQ(iErr)

      Call MeshSetMaxDof(AppCtx%MeshTopology%Mesh, AppCtx%MeshTopology%Num_Dim, iErr); CHKERRQ(iErr) 
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
      Call EXOProperty_InitBCFlag2D(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%BCFlagU)
#elif defined PB_3D
      Call EXOProperty_InitBCFlag3D(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%BCFlagU)
#endif

      AppCtx%TimeStep = 1
      !!! Read U, F, and Temperature
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(Rupt_VertVar_DisplacementX)%Offset, AppCtx%TimeStep, AppCtx%U) 
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(Rupt_VertVar_ForceX)%Offset, AppCtx%TimeStep, AppCtx%F) 
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(Rupt_VertVar_Temperature)%Offset, AppCtx%TimeStep, AppCtx%Theta) 
   End Subroutine ElastInit
!----------------------------------------------------------------------------------------!      
! Solve U (CM)   
!----------------------------------------------------------------------------------------!      
   
   Subroutine Solve_U(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      KSPConvergedReason                           :: reason
      PetscInt                                     :: KSPNumIter
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
      Type(Vec)                                    :: U_Vec
      
!      Call PetscLogStagePush(AppCtx%LogInfo%KSPSolve_Stage, iErr); CHKERRQ(iErr)
      Call MeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%U, U_Vec, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(AppCtx%U, AppCtx%ScatterVect, SCATTER_FORWARD, U_Vec, ierr); CHKERRQ(ierr)
      Call KSPSolve(AppCtx%KSPU, AppCtx%RHSU, U_Vec, iErr); CHKERRQ(iErr)
      !!! Solve and store the solution in AppCtx%RHS
      
      Call SectionRealToVec(AppCtx%U, AppCtx%ScatterVect, SCATTER_REVERSE, U_Vec, ierr); CHKERRQ(ierr)
      !!! Scatter the solution from (Vec) AppCtx%RHS to (SectionReal) AppCtx%U
      
      Call KSPGetIterationNumber(AppCtx%KSPU, KSPNumIter, iErr); CHKERRQ(iErr)
      Call KSPGetConvergedReason(AppCtx%KSPU, reason, iErr); CHKERRQ(iErr)
      Write(IOBuffer, 100) KSPNumIter, reason
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      
      Call VecDestroy(U_Vec, iErr); CHKERRQ(iErr)
!      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
100 Format('KSP converged in ', I5, ' iterations. KSPConvergedReason is ', I2, '\n'c)
   End Subroutine Solve_U

!----------------------------------------------------------------------------------------!      
! MatAssembly (CM)  
!----------------------------------------------------------------------------------------!      
   Subroutine MatAssembly_U(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iBlk, iE, iELoc, iErr
      PetscReal, Dimension(:,:), Pointer           :: MatElem
      
!      Call PetscLogStagePush(AppCtx%LogInfo%MatAssembly_Stage, iErr); CHKERRQ(iErr)
      
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Allocate(MatElem(AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim))
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            !--------------------------------------------------------
            ! The only line to change as a function of the problem.
            Call MatAssembly_U_Local(iE, AppCtx%MatProp( AppCtx%MeshTopology%Elem_Blk(iBlk)%ID ), AppCtx, MatElem)
            !--------------------------------------------------------
            Call assembleMatrix(AppCtx%KU, AppCtx%MeshTopology%mesh, AppCtx%U, iE-1, MatElem, ADD_VALUES, iErr); CHKERRQ(iErr)
         End Do Do_Elem_iE
         DeAllocate(MatElem)
      End Do Do_Elem_iBlk
      Call MatAssemblyBegin(AppCtx%KU, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (AppCtx%KU, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)

!      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine MatAssembly

Subroutine MatAssemblyU(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iBlk, iE, iELoc, iErr
      PetscReal, Dimension(:,:), Pointer           :: MatElem
End Subroutine MatAssembly

!----------------------------------------------------------------------------------------!      
! MatAssembly_U_Local (CM)
!----------------------------------------------------------------------------------------!      
   Subroutine MatAssembly_U_Local(iE, MatProp, AppCtx, MatElem)
#if defined PB_2D
      Type(MatProp2D_Type)                         :: MatProp
#elif defined PB_3D
      Type(MatProp3D_Type)                         :: MatProp
#endif
      Type(AppCtx_Type)                            :: AppCtx
      PetscReal, Dimension(:,:), Pointer           :: MatElem      
      PetscInt                                     :: iE
   
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoFScal, NumDoFVect, NumGauss
      PetscInt, Dimension(:), Pointer              :: BCFlag
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      
      PetscReal, Dimension(:), Pointer             :: V
      PetscReal                                    :: V_Elem
      
!      Call PetscLogEventBegin(AppCtx%LogInfo%MatAssemblyLocal_Event, iErr); CHKERRQ(iErr)

      MatElem  = 0.0_Kr
      NumDoFVec  = Size(AppCtx%ElemVect(iE)%BF,1)
      NumDoFScal = Size(AppCtx%ElemScal(iE)%BF,1)
      NumGauss   = Size(AppCtx%ElemVect(iE)%BF,2)
      
      Allocate(V(NumDoFScal))
      Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%V, iE-1, NumDoFScal, V, iErr); CHKERRQ(ierr)

      Allocate(BCFlag(NumDoFVect))
      Call MeshRestrictClosureInt(AppCtx%MeshTopology%mesh, AppCtx%BCFlagU, iE-1, NumDoF, BCFlag, iErr); CHKERRQ(ierr)
      
      Do iGauss = 1, NumGauss
      !! Calculate V at the gauss point
         V_Elem = 0.0_Kr        
         Do iDoF1 = 1, NumDoFScal
            V_Elem = V_Elem + AppCtx%ElemScal(iE)%BF(iDoF1, iGauss)* V(iDoF1)
         End Do
      !! Assemble the element stiffness
         Do iDoF1 = 1, NumDoFVect
            If (BCFlag(iDoF1) == 0) Then
               Do iDoF2 = 1, NumDoFVect
                  MatElem(iDoF2, iDoF1) =  MatElem(iDoF2, iDoF1) +  AppCtx%ElemVect(iE)%Gauss_C(iGauss) * ((V_Elem**2)*(MatProp%Hookes_Law * AppCtx%ElemVect(iE)%GradS_BF(iDoF1, iGauss)) .DotP. AppCtx%ElemVect(iE)%GradS_BF(iDoF2, iGauss))
!			      Call PetscLogFlops(AppCtx%MeshTopology%num_dim * (AppCtx%MeshTopology%num_dim-1) +1 , iErr);CHKERRQ(iErr)
                  !!! Is that right?
               End Do
            End If
         End Do
      End Do
      DeAllocate(V)
      DeAllocate(BCFlag)
!      Call PetscLogEventEnd(AppCtx%LogInfo%MatAssemblyLocal_Event, iErr); CHKERRQ(iErr)
   End Subroutine MatAssembly_U_Local

!----------------------------------------------------------------------------------------!      
! RHSAssembly (CM)  
!----------------------------------------------------------------------------------------!      
   Subroutine RHS_U_Assembly(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      Type(SectionReal)                            :: RHSSec
      PetscInt                                     :: iBlk, iE, iELoc
      PetscReal, Dimension(:), Pointer             :: RHSElem
      PetscInt                                     :: NumDoFPerVertex 
      
      NumDoFPerVertex = AppCtx%MeshTopology%Num_Dim
      
      !!! Hopefully one day we will use assemble Vector instead of going through a section
!      Call PetscLogStagePush(AppCtx%LogInfo%RHSAssembly_Stage, iErr); CHKERRQ(iErr)
      Call MeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, NumDoFPerVertex, RHSSec, iErr); CHKERRQ(iErr)
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Allocate(RHSElem(AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim))
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            Call RHS_U_AssemblyLocal(iE, AppCtx%MatProp( AppCtx%MeshTopology%Elem_Blk(iBlk)%ID ),AppCtx, RHSElem)
            Call MeshUpdateAddClosure(AppCtx%MeshTopology%Mesh, RHSSec, iE-1, RHSElem, iErr); CHKERRQ(iErr)
         End Do Do_Elem_iE
         DeAllocate(RHSElem)
      End Do Do_Elem_iBlk
      Call SectionRealComplete(RHSSec, iErr); CHKERRQ(iErr)
      !!! VERY important! This is the equivalent of a ghost update
      Call SectionRealToVec(RHSSec, AppCtx%ScatterVect, SCATTER_FORWARD, AppCtx%RHSU, ierr); CHKERRQ(ierr)
      Call SectionRealDestroy(RHSSec, iErr); CHKERRQ(iErr)
!      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
  End Subroutine RHS_U_Assembly
   
!----------------------------------------------------------------------------------------!      
! RHSAssemblyLocal (CM)  
!----------------------------------------------------------------------------------------!      
   Subroutine RHS_U_AssemblyLocal(iE, MatProp, AppCtx, RHSElem)

#if defined PB_2D
      Type(MatProp2D_Type)                         :: MatProp
#elif defined PB_3D
      Type(MatProp3D_Type)                         :: MatProp
#endif
   
      PetscInt                                     :: iE
      Type(AppCtx_Type)                            :: AppCtx
      PetscReal, Dimension(:), Pointer             :: RHSElem
      PetscReal, Dimension(:), Pointer             :: F
      PetscReal, Dimension(:), Pointer             :: Theta
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoFScal, NumDoFVect, NumGauss
      PetscInt, Dimension(:), Pointer              :: BCFlag
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      PetscReal                                    :: Theta_Elem
#if defined PB_2D
      Type (Vect2D)             				         :: TmpRHS
#elif defined PB_3D  
      Type (Vect3D)             				         :: TmpRHS    
#endif

      RHSElem    = 0.0_Kr
      NumDoFVect = Size(AppCtx%ElemVect(iE)%BF,1)
      NumDoFScal = Size(AppCtx%ElemScal(iE)%BF,1)
      NumGauss   = Size(AppCtx%ElemVect(iE)%BF,2)
      Allocate(BCFlag(NumDoFVect))
      Call MeshRestrictClosureInt(AppCtx%MeshTopology%mesh, AppCtx%BCFlagU, iE-1, NumDoFVect, BCFlag, iErr); CHKERRQ(ierr)
      Allocate(F(NumDoFVect))
      Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%F, iE-1, NumDoFVect, F, iErr); CHKERRQ(ierr)
      Allocate(Theta(NumDoFScal))
      Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%Theta, iE-1, NumDoFScal, Theta, iErr); CHKERRQ(ierr)
      Allocate(V(NumDoFScal))
      Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%V, iE-1, NumDoFScal, V, iErr); CHKERRQ(ierr)
      Do_iGauss: Do iGauss = 1, NumGauss
         TmpRHS     = 0.0_Kr
         Theta_Elem = 0.0_Kr
         V_Elem     = 0.0_Kr
         Do iDoF2 = 1, NumDoFVect
            TmpRHS = TmpRHS + AppCtx%ElemVect(iE)%BF(iDoF2, iGauss) * F(iDoF2)
         End Do
         Do iDoF2 = 1, NumDoFScal
            Theta_Elem = Theta_Elem + AppCtx%ElemScal(iE)%BF(iDoF2, iGauss)* Theta(iDoF2)
            V_Elem = V_Elem + AppCtx%ElemVect(iE)%BF(iDoF2, iGauss) * F(iDoF2)
         End Do
         Do iDoF1 = 1, NumDoFVect
            If (BCFlag(iDoF1) == 0) Then
               ! RHS terms due to forces
               RHSElem(iDoF1) = RHSElem(iDoF1) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * ( AppCtx%ElemVect(iE)%BF(iDoF1, iGauss) .DotP. TmpRHS ) 
               ! RHS terms due to inelastic strains
   			   RHSElem(iDoF1) = RHSElem(iDoF1) + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * Theta_Elem * ((V_Elem**2)*(MatProp%Hookes_Law * AppCtx%ElemVect(iE)%GradS_BF(iDoF1, iGauss)) .DotP. MatProp%Therm_Exp)
!               Call PetscLogFlops(3 , iErr);CHKERRQ(iErr)
            End If
         End Do
      End Do Do_iGauss
      DeAllocate(BCFlag)
      DeAllocate(F)
      DeAllocate(Theta)
      DeAllocate(V)
!      Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyLocal_Event, iErr); CHKERRQ(iErr)
   End Subroutine RHS_U_AssemblyLocal
      
   
   Subroutine Save(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx

      Call Write_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(Rupt_VertVar_DisplacementX)%Offset, AppCtx%TimeStep, AppCtx%U) 
      Write(*,*) 'Wrote Result at offset ', AppCtx%MyEXO%VertVariable(Rupt_VertVar_DisplacementX)%Offset
!!!      Call Write_EXO_Result_Cell(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(Rupt_CellVar_StrainXX)%Offset, AppCtx%TimeStep, AppCtx%StrainU) 
!!!      Write(*,*) 'Wrote Strains at offset ', AppCtx%MyEXO%VertVariable(Rupt_CellVar_StrainXX)%Offset
!!!      Call Write_EXO_Result_Cell(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(Rupt_CellVar_StressXX)%Offset, AppCtx%TimeStep, AppCtx%StressU) 
!!!      Write(*,*) 'Wrote Stresses at offset ', AppCtx%MyEXO%VertVariable(Rupt_CellVar_StressXX)%Offset
!!!      Call Write_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(Rupt_GlobVar_Load)%Offset, AppCtx%TimeStep, AppCtx%BulkEnergy)
!!!      Write(*,*) 'Wrote Energy at offset ', AppCtx%MyEXO%GlobVariable(Rupt_GlobVar_Load)%Offset
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
      
      Call VecScatterDestroy(AppCtx%ScatterVect, iErr); CHKERRQ(iErr)
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
End Module m_ m_AmbrosioTortorelli_U2D
#elif defined PB_3D
End Module m_ m_AmbrosioTortorelli_U3D
#endif
