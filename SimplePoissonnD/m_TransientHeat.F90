#if defined PB_2D
Module m_TransientHeat2D
#elif defined PB_3D
Module m_TransientHeat3D
#endif

#include "finclude/petscdef.h"

   Use m_MEF90

#if defined PB_2D
   Use m_Poisson2D
#elif defined PB_3D 
   Use m_Poisson3D
#endif   

   Implicit NONE   

   PetscInt, Parameter, Public                     :: Poisson_Num_EBProperties  = 2
   PetscInt, Parameter, Public                     :: VarFrac_EBProp_HasBForce  = 1
   PetscInt, Parameter, Public                     :: VarFrac_EBProp_Elem_Type  = 2
   
   PetscInt, Parameter, Public                     :: Poisson_Num_SSProperties  = 2
   PetscInt, Parameter, Public                     :: VarFrac_SSProp_HasSForce  = 1
   PetscInt, Parameter, Public                     :: VarFrac_SSProp_Elem_Type  = 2

   PetscInt, Parameter, Public                     :: Poisson_Num_NSProperties  = 1
   PetscInt, Parameter, Public                     :: VarFrac_NSProp_HasPForce  = 1
   
Type TestCase_Type
   PetscInt                                  :: Index
   Character(len=MEF90_MXSTRLEN)             :: Description
End Type

Contains


#undef __FUNCT__
#define __FUNCT__ "PoissonEXOProperty_Init"
   Subroutine PoissonEXOProperty_Init(dEXO, dMeshTopology)
      Type(EXO_Type)                      :: dEXO
      Type(MeshTopology_Type)             :: dMeshTopology
      PetscInt                            :: i, iErr
      PetscInt                            :: NumEB, NumSS, NumNS
      Integer                             :: EXO_MyRank
          
      Call MPI_COMM_RANK(dEXO%Comm, EXO_MyRank, iErr)

      NumEB = dMeshTopology%Num_Elem_Blks_Global
      NumSS = dMeshTopology%Num_Side_Sets_Global
      NumNS = dMeshTopology%Num_Node_Sets_Global
      
      If ( (NumEB == 0) .AND. (NumSS == 0) .AND. (NumSS ==0) ) Then
         Call PetscPrintf(PETSC_COMM_WORLD, '[WARNING]: The EXO file contains no EB, SS or NS is this right?\n', iErr); CHKERRQ(iErr)
         Call PetscPrintf(PETSC_COMM_WORLD, '           Was Write_MeshTopologyGlobal called before VarFracEXOProperty_Init?\n', iErr); CHKERRQ(iErr)
      End If

      dEXO%Num_EBProperties = Poisson_Num_EBProperties
      Allocate(dEXO%EBProperty(dEXO%Num_EBProperties))
      dEXO%EBProperty(VarFrac_EBProp_HasBForce)%Name = 'Has_BForce'
      dEXO%EBProperty(VarFrac_EBProp_Elem_Type)%Name = 'Elem_Type'
      Do i = 1, dEXO%Num_EBProperties
         Allocate(dEXO%EBProperty(i)%Value(NumEB))
         dEXO%EBProperty(i)%Value = 0
      End Do
      
      dEXO%Num_SSProperties = Poisson_Num_SSProperties
      Allocate(dEXO%SSProperty(dEXO%Num_SSProperties))
      dEXO%SSProperty(VarFrac_SSProp_HasSForce)%Name = 'Has_SForce'
      dEXO%SSProperty(VarFrac_SSProp_Elem_Type)%Name = 'Elem_Type'
      Do i = 1, dEXO%Num_SSProperties
         Allocate(dEXO%SSProperty(i)%Value(NumSS))
         dEXO%SSProperty(i)%Value = 0
      End Do
      
      dEXO%Num_NSProperties = Poisson_Num_NSProperties
      Allocate(dEXO%NSProperty(dEXO%Num_NSProperties))
      dEXO%NSProperty(VarFrac_NSProp_HasPForce)%Name = 'Has_PForce'
      Do i = 1, dEXO%Num_NSProperties
         Allocate(dEXO%NSProperty(i)%Value(NumNS))
         dEXO%NSProperty(i)%Value = 0
      End Do
   End Subroutine PoissonEXOProperty_Init   

#undef __FUNCT__
#define __FUNCT__ "PoissonInit"
   Subroutine PoissonInit(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr
      PetscBool                                    :: HasPrefix, Flag
      PetscInt                                     :: iBlk, i, j, iloc 
      Character(len=MEF90_MXSTRLEN)                :: BatchFileName
      PetscInt                                     :: BatchUnit=99
      PetscInt                                     :: NumTestCase
      Type(TestCase_Type), Dimension(:) , Pointer  :: TestCase
      PetscBool                                    :: EraseBatch
      PetscBool                                    :: IsBatch, HasBatchFile
      Character(len=MEF90_MXSTRLEN)                :: prefix, IOBuffer, filename   
      Type(DM)                                     :: Tmp_Mesh
      PetscReal                                    :: ValU, ValF
      PetscInt, Dimension(:), Pointer              :: SizeScal
      PetscInt, Parameter                          :: VarFrac_NSProp_BCT =1
      PetscReal, Dimension(:), Pointer             :: U, Uelem


      Call MEF90_Initialize()
      
      NumTestCase = 2
      Allocate(TestCase(NumTestCase))
      Do i = 1, NumTestCase
         TestCase(i)%Index = i
      End Do
      TestCase(1)%Description = "Simple Poisson \Delta u = f"
      TestCase(2)%Description = "Heat equation u,t - \Delta u = f"
      
      AppCtx%AppParam%verbose = 0
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-verbose', AppCtx%AppParam%verbose, Flag, iErr); CHKERRQ(iErr)
      
      Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p',       AppCtx%AppParam%prefix, HasPrefix, iErr); CHKERRQ(iErr)
      If (.NOT. HasPrefix) Then
         Call PetscPrintf(PETSC_COMM_WORLD, "No mesh prefix given\n", iErr)
         Call MEF90_Finalize()
         STOP
      End If
      Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-i', BatchFileName, IsBatch, iErr);CHKERRQ(iErr)

      EraseBatch=.False.
      Call PetscOptionsGetBool(PETSC_NULL_CHARACTER, '-force', EraseBatch, HasPrefix, iErr)    
      If (MEF90_MyRank==0) Then
         If (IsBatch) Then
            Write(IOBuffer, *) "\nProcessing batch file ", Trim(BatchFileName), "\n"
            Call PetscPrintf(PETSC_COMM_SELF, IOBuffer, iErr); CHKERRQ(iErr)
            Open(Unit=BatchUnit, File=BatchFileName, Status='Old', Action='Read')
            Rewind(BatchUnit)
         Else
            BatchFileName = Trim(prefix)//'.args'
            Inquire(File=BatchFileName, EXIST=HasBatchFile)
            If (HasBatchFile .AND. (.NOT. EraseBatch)) Then
               Write(IOBuffer, *) "Batch file ", trim(BatchFileName), " already exists. Erase it or use -force flag\n"
               Call PetscPrintf(PETSC_COMM_SELF, IOBuffer, iErr); CHKERRQ(iErr)
               SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_UNEXPECTED, IOBuffer, iErr)
            Else
               Write(IOBuffer, *) "Running interactively and generating batch file ", trim(BatchFileName), "\n"
               Call PetscPrintf(PETSC_COMM_SELF, IOBuffer, iErr); CHKERRQ(iErr)
               Open(Unit=BatchUnit, File=BatchFileName, Status='Unknown')
               Rewind(BatchUnit)
            End If
         End If
      End If

      AppCtx%AppParam%TestCase = 1
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-test',       AppCtx%AppParam%TestCase, Flag, iErr); CHKERRQ(iErr)
      Call InitLog(AppCtx)
      Call PetscLogStagePush(AppCtx%LogInfo%Setup_Stage, iErr); CHKERRQ(iErr)
      If (AppCtx%AppParam%verbose > 1) Then
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
102 Format('Output from processor ', I4.4, ' redirected to file ', A, '\n')
103 Format(A,'.log')
104 Format('Collective output redirected to file ', A, '\n')

      AppCtx%EXO%Comm = PETSC_COMM_WORLD
      AppCtx%EXO%filename = Trim(AppCtx%AppParam%prefix)//'.gen'
      !!! Read and partition the mesh
      If (MEF90_NumProcs == 1) Then
         Call PetscLogStagePush(AppCtx%LogInfo%MeshCreateExodus_Stage, iErr); CHKERRQ(iErr)
         Call DMMeshCreateExodus(PETSC_COMM_WORLD, AppCtx%EXO%filename, AppCtx%MeshTopology%mesh, ierr); CHKERRQ(iErr)
         Call PetscLogStagePop(iErr); CHKERRQ(iErr)
      Else
         Call PetscLogStagePush(AppCtx%LogInfo%MeshCreateExodus_Stage, iErr); CHKERRQ(iErr)
         Call DMMeshCreateExodus(PETSC_COMM_WORLD, AppCtx%EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
         Call PetscLogStagePop(iErr); CHKERRQ(iErr)
      
         Call PetscLogStagePush(AppCtx%LogInfo%MeshDistribute_Stage, iErr); CHKERRQ(iErr)
         Call DMMeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, AppCtx%MeshTopology%mesh, ierr); CHKERRQ(iErr)
         Call DMDestroy(Tmp_mesh, ierr); CHKERRQ(iErr)
         Call PetscLogStagePop(iErr); CHKERRQ(iErr)
      End If

      Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
      Call MeshTopologyReadEXO(AppCtx%MeshTopology, AppCtx%EXO)

      !!! Sets the type of elements for each block
      Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_Type = MEF90_P1_Lagrange
         Call Init_Elem_Blk_Type(AppCtx%MeshTopology%Elem_Blk(iBlk), AppCtx%MeshTopology%num_dim)
      End Do
   
      Call ElementInit(AppCtx%MeshTopology, AppCtx%Elem, 2)
      
      !!! Allocate the Section for U and F
      Allocate(SizeScal(1))
      SizeScal=1
      Call FieldCreateVertex(AppCtx%U,     'U',         AppCtx%MeshTopology, SizeScal)
      Call FieldCreateVertex(AppCtx%F,     'F',         AppCtx%MeshTopology, SizeScal)
      Call FieldCreateVertex(AppCtx%RHS,   'RHS',       AppCtx%MeshTopology, SizeScal)
      Call FlagCreateVertex(AppCtx%BCFlag, 'BC',        AppCtx%MeshTopology, SizeScal)
      Call FieldCreateVertex(AppCtx%U_0,   'U_0',        AppCtx%MeshTopology, SizeScal)
      DeAllocate(SizeScal)

      
      AppCtx%MyEXO%comm = PETSC_COMM_SELF
      AppCtx%MyEXO%exoid = AppCtx%EXO%exoid
      Write(AppCtx%MyEXO%filename, 200) trim(AppCtx%AppParam%prefix), MEF90_MyRank
200 Format(A, '-', I4.4, '.gen')
      AppCtx%MyEXO%title = trim(AppCtx%EXO%title)
      AppCtx%MyEXO%Num_QA = AppCtx%EXO%Num_QA
!!! Prepare and format the output mesh   
      Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
      Call MeshTopologyWriteGlobal(AppCtx%MeshTopology, AppCtx%MyEXO, PETSC_COMM_WORLD)
      Call EXOFormat_SimplePoisson(AppCtx)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)

      Write(IOBuffer, *) '\nTest Case:\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)      
      Do i = 1, NumTestCase
         Write(IOBuffer, "('   [',I2.2,'] ',A)"), TestCase(i)%Index,   Trim(TestCase(i)%Description)//'\n'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr);   CHKERRQ(iErr)
      End Do

   !!!
   !!! EB, SS, NS Properties
   !!!
      Call PoissonEXOProperty_Init(AppCtx%MyEXO, AppCtx%MeshTopology) 
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) "Done with PoissonEXOProperty_Init\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
   
      Write(IOBuffer, *) '\nElement Block and Node Set Properties\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)      

      Call EXOEBProperty_AskWithBatch(AppCtx%MyEXO, AppCtx%MeshTopology, BatchUnit, IsBatch)
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) "Done with EXOEBProperty_AskWithBatch\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      Call EXOSSProperty_AskWithBatch(AppCtx%MyEXO, AppCtx%MeshTopology, BatchUnit, IsBatch)
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) "Done with EXOSSProperty_AskWithBatch\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      Call EXONSProperty_AskWithBatch(AppCtx%MyEXO, AppCtx%MeshTopology, BatchUnit, IsBatch)
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) "Done with EXONSProperty_AskWithBatch\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If


      Call MEF90_AskInt(AppCtx%AppParam%TestCase, 'Test Case', BatchUnit, IsBatch)
!TODO set these parameters in the PrepPoisson from the.args file 
      Select Case(AppCtx%AppParam%TestCase)
      Case(2)
      !AppCtx%maxsteps = 1000
      !AppCtx%maxtime =10.0
         Call MEF90_AskInt(AppCtx%maxsteps, 'Max number of steps for TS computation', BatchUnit, IsBatch)
         Call MEF90_AskReal(AppCtx%maxtime,  'Max time for TS computation', BatchUnit, IsBatch)
      End Select

!  Set EB Properties : U, F         

         ! This has no impact for TestCase 1
      Call MEF90_AskReal(ValU, 'Initial value in U ', BatchUnit, IsBatch)
      Call SectionRealSet(AppCtx%U%Sec, ValU, iErr); CHKERRQ(iErr);
      Call DMMeshGetVertexSectionReal(AppCtx%MeshTopology%mesh, 'U_0', 1,   AppCtx%U_0%Sec, iErr); CHKERRQ(iErr)
      call SectionRealSet(AppCtx%U_0%Sec, ValU, iErr); CHKERRQ(iErr) !This should be replaced by evaluating the initial solution or most probabaly reading from .args file
      Call SectionRealToVec(AppCtx%U_0, AppCtx%U_0%Scatter, SCATTER_FORWARD, AppCtx%U_0%Vec,ierr)
! TODO ligne précedente pour U et F ?? 
      !Setting force
      Call MEF90_AskReal(ValF, 'RHS F', BatchUnit, IsBatch)
      Call SectionRealSet(AppCtx%F%Sec, ValF, iErr); CHKERRQ(iErr);

!Get BC values 
!lines 535-550 PrepVarFracNG
      Allocate(U(AppCtx%MeshTopology%Num_Node_Sets)) 
      U = 0.0_Kr
      Do i = 1, AppCtx%MeshTopology%Num_Node_Sets_Global
         Write(IOBuffer, 202) i
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         If (AppCtx%MyEXO%NSProperty(VarFrac_NSProp_HasPForce)%Value(i) /= 0 ) Then
            Write(IOBuffer, 302) i, 'Ux'
            Call MEF90_AskReal(U(i), IOBuffer, BatchUnit, IsBatch)
         End If
      End Do

202 Format('    Node Set      ', T24, I3, '\n')
302 Format('NS', I4.4, ': ', A)

! Set BC on NS 
! lines 635-647 PrepVarfracNG
      Allocate(Uelem(1))
      Do iloc = 1, AppCtx%MeshTopology%Num_Node_Sets         
         i = AppCtx%MeshTopology%Node_Set(iloc)%ID
            Uelem(1) = U(i)
            Do j = 1, AppCtx%MeshTopology%Node_Set(iloc)%Num_Nodes
               Call SectionRealUpdate(AppCtx%U%Sec, AppCtx%MeshTopology%Num_Elems + AppCtx%MeshTopology%Node_Set(iloc)%Node_ID(j)-1, Uelem, INSERT_VALUES, iErr); CHKERRQ(iErr)
            End Do
      End Do
      DeAllocate(Uelem)
      Call SectionIntAddNSProperty(AppCtx%BCFlag%Sec,  AppCtx%MyEXO%NSProperty(VarFrac_NSProp_BCT),  AppCtx%MeshTopology)

   End Subroutine PoissonInit

#undef __FUNCT__
#define __FUNCT__ "KSPSetUp"
   Subroutine KSPSetUp(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr
      

      !!! Initialize the matrix and vector for the linear system
      Call DMMeshSetMaxDof(AppCtx%MeshTopology%Mesh, 1, iErr); CHKERRQ(iErr) 
      Call DMMeshCreateMatrix(AppCtx%MeshTopology%mesh, AppCtx%U%Sec, MATMPIAIJ, AppCtx%K, iErr); CHKERRQ(iErr)
     

      Call KSPCreate(PETSC_COMM_WORLD, AppCtx%KSP, iErr); CHKERRQ(iErr)
      Call KSPSetOperators(AppCtx%KSP, AppCtx%K, AppCtx%K, SAME_NONZERO_PATTERN, iErr); CHKERRQ(iErr)
      Call KSPSetInitialGuessNonzero(AppCtx%KSP, PETSC_TRUE, iErr); CHKERRQ(iErr)
      Call KSPGetPC(AppCtx%KSP, AppCtx%PC, iErr); CHKERRQ(iErr)
      Call PCSetType(AppCtx%PC, PCBJACOBI, iErr); CHKERRQ(iErr)
      Call KSPSetType(AppCtx%KSP, KSPCG, iErr); CHKERRQ(iErr)
      Call KSPAppendOptionsPrefix(AppCtx%KSP, "U_", iErr); CHKERRQ(iErr)
      Call KSPSetFromOptions(AppCtx%KSP, iErr); CHKERRQ(iErr)
      
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
      
      
      !!! Read Force and BC from Data file or reformat it
      !call PoissonBC(AppCtx) 
      !!! Create the EXO case file
      Call Write_EXO_Case(AppCtx%AppParam%prefix, '%0.4d', MEF90_NumProcs)
      
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine KSPSetUP
   
#undef __FUNCT__
#define __FUNCT__ "TSSetUp"
   Subroutine TSSetUP(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr

      !!! Initialize the matrix and vector for the linear system
      Call DMMeshSetMaxDof(AppCtx%MeshTopology%Mesh, 1, iErr); CHKERRQ(iErr) 
      Call DMMeshCreateMatrix(AppCtx%MeshTopology%mesh, AppCtx%U%sec, MATMPIAIJ, AppCtx%K, iErr); CHKERRQ(iErr)
      Call DMMeshCreateMatrix(AppCtx%MeshTopology%mesh, AppCtx%U%sec, MATMPIAIJ, AppCtx%M, iErr); CHKERRQ(iErr)
      Call DMMeshCreateMatrix(AppCtx%MeshTopology%mesh, AppCtx%U%sec, MATMPIAIJ, AppCtx%Jac, iErr); CHKERRQ(iErr)
      
      Call TSCreate(PETSC_COMM_WORLD, AppCtx%TS, iErr); CHKERRQ(iErr)
      
      
      Select Case(AppCtx%AppParam%TestCase)
      Case(2)
         Call TSSetType(AppCtx%TS, TSARKIMEX, ierr); CHKERRQ(iErr)  
      !Call TSSetProblemType(AppCtx%TS, TS_LINEAR, ierr); CHKERRQ(iErr)
      !Call TSSetType(AppCtx%TS, TSBEULER, ierr);  CHKERRQ(iErr)
      End Select

      !call DMMeshSetSectionReal(AppCtx%MeshTopology%mesh, "default", AppCtx%U%sec).
      call TSSetDM(AppCtx%TS, AppCtx%MeshTopology%mesh, ierr);   CHKERRQ(iErr) 


      !Call DMMeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%U_0%Sec, AppCtx%U_0%Vec, iErr); CHKERRQ(iErr)
      Call TSSetSolution(AppCtx%TS, AppCtx%U_0%Vec,  ierr);   CHKERRQ(iErr) 

      Call TSSetInitialTimeStep(AppCtx%TS, 0.0, .01,  ierr); CHKERRQ(iErr)
      Call TSSetDuration(AppCtx%TS, AppCtx%maxsteps, AppCtx%maxtime, iErr); CHKERRQ(iErr)
      
      Call TSSetIFunction(AppCtx%TS, PETSC_NULL_OBJECT, IFunctionPoisson, AppCtx, ierr); CHKERRQ(iErr)
      Call TSSetIJacobian(AppCtx%TS, AppCtx%Jac, AppCtx%Jac, IJacobianPoisson, AppCtx, ierr); CHKERRQ(iErr)
      Call TSSetRHSFunction(AppCtx%TS, PETSC_NULL_OBJECT, RHSPoisson, AppCtx, ierr); CHKERRQ(iErr)
      
      Call TSSetFromOptions(AppCtx%TS,  iErr); CHKERRQ(iErr)

      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
       
      !!! Create the EXO case file
      Call Write_EXO_Case(AppCtx%AppParam%prefix, '%0.4d', MEF90_NumProcs)

      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine TSSetUp


#undef __FUNCT__
#define __FUNCT__ "IFunctionPoisson"
   SubRoutine IFunctionPoisson(dummyTS, t, U, Udot, GlobalOut, AppCtx, iErr)
      PetscReal                                    :: t 
      Type(Vec)                                    :: U, Udot, GlobalOut, dummyVec
      Type(AppCtx_Type)                            :: AppCtx
      Type(TS)                                     :: dummyTS
      PetscInt                                     :: iErr
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer   
         
  
      Write(IOBuffer, *) 'Begin IFunctionPoisson \n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call DMMeshCreateVector(AppCtx%MeshTopology%mesh, AppCtx%U_0%Sec, dummyVec, iErr); CHKERRQ(iErr)

       !GlobalOut = AppCtx%M * Udot + AppCtx%K * U
      !Call VecView(U,  PETSC_VIEWER_STDOUT_SELF, iErr); CHKERRQ(iErr)
      Call MatMult(AppCtx%K, U, GlobalOut, iErr);CHKERRQ(iErr)
!Udot is non at the second iteration; why ?? 
      Call VecView(Udot,  PETSC_VIEWER_STDOUT_SELF, iErr); CHKERRQ(iErr)
      !Call VecView(dummyVec,  PETSC_VIEWER_STDOUT_SELF, iErr); CHKERRQ(iErr)
      Call MatMultAdd(AppCtx%M, Udot, GlobalOut, GlobalOut, iErr);CHKERRQ(iErr)
!TODO Can the constant and result factors be the same ? 
!      Call VecView(GlobalOut,  PETSC_VIEWER_STDOUT_SELF, iErr); CHKERRQ(iErr)
       
   End Subroutine IFunctionPoisson


#undef __FUNCT__
#define __FUNCT__ "IJacobianPoisson"
   SubRoutine IJacobianPoisson(dummyTS, t, U, Udot, a, Jac, PreJac, AppCtx, iErr)
      Type(TS)                                     :: dummyTS
      PetscReal                                    :: t, a 
      Type(Vec)                                    :: U, Udot
      Type(Mat)                                    :: Jac, PreJac  
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer   
!Pk est ce que Ijacobian n'est pas appelé ????????? 
! a*AppCtx%M + AppCtx%K
      Write(IOBuffer, *) 'Begin IJacobianPoisson \n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
!Initialisation de la matrice Jac ?
!!! Att manque un copy ! la on a Y = aX + Y 
!! Condition aux limites ?? 
      call MatAXPY(Jac, a, AppCtx%M, AppCtx%K, iErr); CHKERRQ(iErr)

   End Subroutine IJacobianPoisson
 
#undef __FUNCT__
#define __FUNCT__ "RHSPoisson"
   SubRoutine RHSPoisson(dummyTS, t, U, GlobalOut, AppCtx, iErr)
      Type(TS)                                     :: dummyTS
      PetscReal                                    :: t 
      Type(Vec)                                    :: U, GlobalOut
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer   
!      Type(VecScatter)                             :: dummyScatter
!TODO est ce que Globalout est correctement initialise ? 
      Write(IOBuffer, *) 'Begin RHSPoisson \n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
!AppCtx%RHS Is a section real we want to return a vector!
!Return AppCtx%RHS
!      Call DMMeshCreateGlobalScatter(AppCtx%MeshTopology%mesh, AppCtx%RHS%Sec,  AppCtx%RHS%Scatter, iErr); CHKERRQ(iErr) 
!      Call SectionRealToVec(AppCtx%RHS%Sec, AppCtx%RHS%Scatter, SCATTER_FORWARD, GlobalOut, iErr); CHKERRQ(iErr)
      Call VecCopy(AppCtx%RHS%Vec, GlobalOut, iErr); CHKERRQ(iErr)
     
   End Subroutine RHSPoisson
   
   
#undef __FUNCT__
#define __FUNCT__ "SolveTransient"
   Subroutine SolveTransient(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
      
      Call PetscLogStagePush(AppCtx%LogInfo%KSPSolve_Stage, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(AppCtx%U%Sec, AppCtx%U%Scatter, SCATTER_FORWARD, AppCtx%U%Vec, iErr); CHKERRQ(iErr)

      ! Using IMEX see section 6.1.2 PetSc documentation 
!      call TSSetFromOptions(AppCtx%TS,ierr)
      Call TSView(AppCtx%TS,  PETSC_VIEWER_STDOUT_WORLD)

      call TSSolve(AppCtx%TS, AppCtx%U%Vec, AppCtx%maxtime, iErr); CHKERRQ(iErr)
       

  !    Write(IOBuffer, 100) KSPNumIter, KSPreason
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      
      Call SectionRealToVec(AppCtx%U%Sec, AppCtx%U%Scatter, SCATTER_REVERSE, AppCtx%U%Vec, ierr); CHKERRQ(ierr)
 
!      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
!100 Format('KSP converged in ', I5, ' iterations. KSPConvergedReason is ', I2, '\n')
   End Subroutine SolveTransient
   
 
     
#undef __FUNCT__
#define __FUNCT__ "MatMassAssembly"
   Subroutine MatMassAssembly(AppCtx)   
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iBlk, iErr
!TOTO Set M to 0 before each iteration before calling assembly      
      Call PetscLogStagePush(AppCtx%LogInfo%MatAssembly_Stage, iErr); CHKERRQ(iErr)
      Do_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Call MatMassAssemblyBlock(iBlk, AppCtx)
      End Do Do_iBlk
      Call MatAssemblyBegin(AppCtx%M, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (AppCtx%M, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)

      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine MatMassAssembly
      
#undef __FUNCT__
#define __FUNCT__ "MatMassAssemblyBlock"
   Subroutine MatMassAssemblyBlock(iBlk, AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iBlk
      
      PetscInt                                     :: iE, iELoc, iErr
      PetscReal, Dimension(:,:), Pointer           :: MatElem
      PetscInt, Dimension(:), Pointer              :: BCFlag
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      PetscLogDouble                               :: flops = 0
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
      
      Call PetscLogEventBegin(AppCtx%LogInfo%MatAssemblyBlock_Event, iErr); CHKERRQ(iErr)
      
      Allocate(MatElem(AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF))
   
      Do_iELoc: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         MatElem = 0.0_Kr
         Do iGauss = 1, size(AppCtx%Elem(iE)%Gauss_C)
            Do iDoF1 = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
               Do iDoF2 = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
                  MatElem(iDoF2, iDoF1) = MatElem(iDoF2, iDoF1) + AppCtx%Elem(iE)%Gauss_C(iGauss) * ( AppCtx%Elem(iE)%BF(iDoF1, iGauss) *  AppCtx%Elem(iE)%BF(iDoF2, iGauss) )
                  flops = flops + 2 !Check 
               End Do
            End Do
         End Do
         Call DMMeshAssembleMatrix(AppCtx%M, AppCtx%MeshTopology%mesh, AppCtx%U%sec, iE-1, MatElem, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iELoc
   
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
      DeAllocate(MatElem)
      Call PetscLogEventEnd(AppCtx%LogInfo%MatAssemblyBlock_Event, iErr); CHKERRQ(iErr)
   End Subroutine MatMassAssemblyBlock

 
#undef __FUNCT__
#define __FUNCT__ "TSPoissonFinalize"
   Subroutine TSPoissonFinalize(AppCtx)   
      Type(AppCtx_Type)                            :: AppCtx

      PetscInt                                     :: iErr
      Character(len=MEF90_MXSTRLEN)                :: filename

      Call FieldDestroy(AppCtx%U)
      Call FieldDestroy(AppCtx%F)
      Call FieldDestroy(AppCtx%RHS)
      Call FieldDestroy(AppCtx%U_0)
!      Call SectionRealDestroy(AppCtx%U_0, iErr); CHKERRQ(iErr)
      Call SectionRealDestroy(AppCtx%GradU, iErr); CHKERRQ(iErr)
      Call SectionIntDestroy(AppCtx%BCFlag, iErr); CHKERRQ(iErr)
      Call MatDestroy(AppCtx%K, iErr); CHKERRQ(iErr)
      Call MatDestroy(AppCtx%M, iErr); CHKERRQ(iErr)
      Call MatDestroy(AppCtx%Jac, iErr); CHKERRQ(iErr)
!      Call VecDestroy(AppCtx%U_0_Vec, iErr); CHKERRQ(iErr)

      Call TSDestroy(AppCtx%TS, iErr); CHKERRQ(iErr)
      Call DMDestroy(AppCtx%MeshTopology%Mesh, iErr); CHKERRQ(ierr)
      
      If (AppCtx%AppParam%verbose > 1) Then
         Call PetscViewerFlush(AppCtx%AppParam%MyLogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerFlush(AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerDestroy(AppCtx%AppParam%MyLogViewer, iErr); CHKERRQ(iErr)
         Call PetscViewerDestroy(AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
      End If
      Write(filename, 103) Trim(AppCtx%AppParam%prefix)
!      Call PetscLogPrintSummary(PETSC_COMM_WORLD, filename, iErr); CHKERRQ(iErr)
      Call MEF90_Finalize()
103 Format(A,'.log.txt')
   End Subroutine TSPoissonFinalize

#if defined PB_2D
End Module m_TransientHeat2D
#elif defined PB_3D
End Module m_TransientHeat3D
#endif
