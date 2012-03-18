#if defined PB_2D
Module m_TransientHeat2D
#elif defined PB_3D
Module m_TransientHeat3D
#endif

#include "finclude/petscdef.h"

   Use m_MEF90
   Use m_MEF_Integrate
   use m_Heat_Struct

   Implicit NONE   

   PetscInt, Parameter, Public                     :: Poisson_Num_EBProperties   = 2
   PetscInt, Parameter, Public                     :: Heat_EBProp_HasBForce      = 1
   PetscInt, Parameter, Public                     :: Heat_EBProp_Elem_Type      = 2
   
   PetscInt, Parameter, Public                     :: Poisson_Num_SSProperties   = 2
   PetscInt, Parameter, Public                     :: Heat_SSProp_HasSForce      = 1
   PetscInt, Parameter, Public                     :: Heat_SSProp_Elem_Type      = 2

   PetscInt, Parameter, Public                     :: Poisson_Num_NSProperties   = 1
   PetscInt, Parameter, Public                     :: Heat_NSProp_HasPForce      = 1
   
   PetscInt, Parameter, Public                     :: Heat_Num_GlobVar           = 4
   PetscInt, Parameter, Public                     :: Heat_GlobVar_ElasticEnergy = 1
   PetscInt, Parameter, Public                     :: Heat_GlobVar_ExtForcesWork = 2
   PetscInt, Parameter, Public                     :: Heat_GlobVar_TotalEnergy   = 3
   PetscInt, Parameter, Public                     :: Heat_GlobVar_Load          = 4
   
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
      Type(Field)                                  :: T
      Type(Field)                                  :: TBC
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
      Type(AppParam_Type)                          :: AppParam
   !For TS
      Type(TS)                                     :: TS
      PetscInt                                     :: NumSteps
      PetscInt                                     :: VertVar_Temperature 
      Type(MatHeat_Type), Dimension(:), Pointer    :: MatProp
      Type(HeatSchemeParam_Type)                   :: HeatSchemeParam 
   End Type Heat_AppCtx_Type

Contains

#undef __FUNCT__
#define __FUNCT__ "ExoFormat_SimplePoisson"
   Subroutine EXOFormat_SimplePoisson(AppCtx)
      Type(Heat_AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr
   
      Call EXPVP (AppCtx%MyEXO%exoid, 'g', 4, iErr)
      Call EXPVAN(AppCtx%MyEXO%exoid, 'g', 4, (/'Elastic Energy ', 'Ext Forces work', 'Total Energy   ', 'Time           '/), iErr)
      Call EXPVP (AppCtx%MyEXO%exoid, 'n', 2, iErr)
      Call EXPVAN(AppCtx%MyEXO%exoid, 'n', 2, (/'T', 'F'/), iErr)
      Call EXPTIM(AppCtx%MyEXO%exoid, 1, 1.0_Kr, iErr)
   End Subroutine EXOFormat_SimplePoisson

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
      dEXO%EBProperty(Heat_EBProp_HasBForce)%Name = 'Has_BForce'
      dEXO%EBProperty(Heat_EBProp_Elem_Type)%Name = 'Elem_Type'
      Do i = 1, dEXO%Num_EBProperties
         Allocate(dEXO%EBProperty(i)%Value(NumEB))
         dEXO%EBProperty(i)%Value = 0
      End Do
      
      dEXO%Num_SSProperties = Poisson_Num_SSProperties
      Allocate(dEXO%SSProperty(dEXO%Num_SSProperties))
      dEXO%SSProperty(Heat_SSProp_HasSForce)%Name = 'Has_SForce'
      dEXO%SSProperty(Heat_SSProp_Elem_Type)%Name = 'Elem_Type'
      Do i = 1, dEXO%Num_SSProperties
         Allocate(dEXO%SSProperty(i)%Value(NumSS))
         dEXO%SSProperty(i)%Value = 0
      End Do
      
      dEXO%Num_NSProperties = Poisson_Num_NSProperties
      Allocate(dEXO%NSProperty(dEXO%Num_NSProperties))
      dEXO%NSProperty(Heat_NSProp_HasPForce)%Name = 'Has_PForce'
      Do i = 1, dEXO%Num_NSProperties
         Allocate(dEXO%NSProperty(i)%Value(NumNS))
         dEXO%NSProperty(i)%Value = 0
      End Do
   End Subroutine PoissonEXOProperty_Init   

#undef __FUNCT__
#define __FUNCT__ "PoissonInit"
   Subroutine PoissonInit(AppCtx)
      Type(Heat_AppCtx_Type)                       :: AppCtx
      PetscInt                                     :: iErr, iDiff
      PetscBool                                    :: HasPrefix, Flag
      PetscInt                                     :: iBlk, i, j, iloc, iEloc 
      Character(len=MEF90_MXSTRLEN)                :: BatchFileName
      PetscInt                                     :: BatchUnit=99
      PetscInt                                     :: NumTestCase
      PetscBool                                    :: EraseBatch
      PetscBool                                    :: IsBatch, HasBatchFile
      Character(len=MEF90_MXSTRLEN)                :: prefix, IOBuffer, filename   
      Type(DM)                                     :: Tmp_Mesh
      PetscInt, Parameter                          :: VarFrac_NSProp_BCT =1
      PetscReal, Dimension(:), Pointer             :: U
      PetscReal, Dimension(:), Pointer             :: ValU, ValF

      PetscReal                                    :: Tmin
      PetscReal                                    :: Tmax
      PetscReal, Dimension(:), Pointer             :: Time
      PetscReal, Dimension(:), Pointer             :: GlobVars
      
      Call MEF90_Initialize()
      
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
            BatchFileName = Trim(AppCtx%AppParam%prefix)//'.args'
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

      Call HeatSchemeParam_GetFromOptions(AppCtx%HeatSchemeParam)
      If (AppCtx%AppParam%verbose > 0) Then
         Call HeatSchemeParam_View(AppCtx%HeatSchemeParam, PetscViewer(PETSC_VIEWER_STDOUT_WORLD))
      End If

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
      AppCtx%EXO%exoid = EXOPEN(AppCtx%EXO%filename, EXREAD, exo_cpu_ws, exo_io_ws, PETSC_NULL_INTEGER, ierr)
      !!! Read and partition the mesh
      If (MEF90_NumProcs == 1) Then
         Call DMMeshCreateExodus(PETSC_COMM_WORLD, AppCtx%EXO%filename, AppCtx%MeshTopology%mesh, ierr); CHKERRQ(iErr)
      Else
         Call DMMeshCreateExodus(PETSC_COMM_WORLD, AppCtx%EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
         Call DMMeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, AppCtx%MeshTopology%mesh, ierr); CHKERRQ(iErr)
         Call DMDestroy(Tmp_mesh, ierr); CHKERRQ(iErr)
      End If

      Call MeshTopologyGetInfo(AppCtx%MeshTopology,PETSC_COMM_WORLD)

      !!! Sets the type of elements for each block
      Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_Type = MEF90_P1_Lagrange
         Call Init_Elem_Blk_Type(AppCtx%MeshTopology%Elem_Blk(iBlk), AppCtx%MeshTopology%num_dim)
      End Do
   
      Call ElementInit(AppCtx%MeshTopology, AppCtx%Elem, 2)
      Call HeatFieldCreate(AppCtx, AppCtx%MeshTopology)

      
      AppCtx%MyEXO%comm = PETSC_COMM_SELF
      AppCtx%MyEXO%exoid = AppCtx%EXO%exoid
      Write(AppCtx%MyEXO%filename, 200) trim(AppCtx%AppParam%prefix), MEF90_MyRank
200 Format(A, '-', I4.4, '.gen')
      AppCtx%MyEXO%exoid = EXCRE (AppCtx%MyEXO%filename, EXCLOB, exo_cpu_ws, exo_io_ws, iErr)
      AppCtx%MyEXO%title = trim(AppCtx%EXO%title)
      AppCtx%MyEXO%Num_QA = AppCtx%EXO%Num_QA
!!! Prepare and format the output mesh   
      Call MeshTopologyWriteGlobal(AppCtx%MeshTopology, AppCtx%MyEXO, PETSC_COMM_WORLD)
      Call EXOFormat_SimplePoisson(AppCtx)

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


      Call MEF90_AskReal(Tmin,  'Min time for TS computation',    BatchUnit, IsBatch)
      Call MEF90_AskReal(Tmax,  'Max time for TS computation',    BatchUnit, IsBatch)
      Call MEF90_AskInt(AppCtx%NumSteps,  'Number of time steps', BatchUnit, IsBatch)

      Allocate(GlobVars(Heat_Num_GlobVar))
      GlobVars = 0.0_Kr
      Allocate(Time(AppCtx%NumSteps))
      Do i = 1, AppCtx%NumSteps-1
         !! Non uniform time stepping adapted to the time scale of the thermal problem in tau=sqrt(t)
         !! Pay attention in visualization softwares
         Time(i) = Tmin + (Real(i-1) * (sqrt(Tmax-TMin))/Real(AppCtx%NumSteps-1))**2

         GlobVars(Heat_GlobVar_Load) = Time(i)
         Call Write_EXO_AllResult_Global(AppCtx%MyEXO, i, GlobVars)
         Call EXPTIM(AppCtx%MyEXO%exoid, i, Time(i), iErr)
      End Do
   
      Time(AppCtx%NumSteps) = Tmax
      GlobVars(Heat_GlobVar_Load) = Time(AppCtx%NumSteps)
      Call Write_EXO_AllResult_Global(AppCtx%MyEXO, AppCtx%NumSteps, GlobVars)
      Call EXPTIM(AppCtx%MyEXO%exoid, AppCtx%NumSteps, Time(AppCtx%NumSteps), iErr)
      DeAllocate(GlobVars)

      Call View_Available_Diffusion_Laws 
!  Set EB Properties : U, F     
      Allocate(ValU(AppCtx%MeshTopology%Num_Elem_Blks_Global)) 
      Allocate(ValF(AppCtx%MeshTopology%Num_Elem_Blks_Global)) 
      Allocate(AppCtx%MatProp(AppCtx%MeshTopology%Num_Elem_Blks_Global)) 
      Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks_Global
         Write(IOBuffer, 300) iBlk, 'Initial Value in U'
         Call MEF90_AskReal(ValU(iBlk), IOBuffer, BatchUnit, IsBatch)
         Write(IOBuffer, 300) iBlk, 'RHS F'
         Call MEF90_AskReal(ValF(iBlk),IOBuffer, BatchUnit, IsBatch)
         Write(IOBuffer, 300) iBlk, 'Diffusion Law'
         Call MEF90_AskInt(AppCtx%MatProp(iBlk)%Type_Law,IOBuffer, BatchUnit, IsBatch)
         Allocate(AppCtx%MatProp(iBlk)%Diffusivity(Heat_Num_Param(AppCtx%MatProp(iBlk)%Type_Law)))
         Do iDiff = 1, Heat_Num_Param(AppCtx%MatProp(iBlk)%Type_Law) 
            Write(IOBuffer, 301) i, 'Diffusivity', iDiff 
            Call MEF90_AskReal(AppCtx%MatProp(iBlk)%Diffusivity(iDiff), IOBuffer, BatchUnit, IsBatch)
         End Do 
      End Do 
      Call HeatSetInitial(AppCtx, AppCtx%MeshTopology, ValU, ValF)
      DeAllocate(ValU)
      DeAllocate(ValF)


!Get BC values 
      Allocate(U(AppCtx%MeshTopology%Num_Node_Sets_Global)) 
      U = 0.0_Kr 
      Do i = 1, AppCtx%MeshTopology%Num_Node_Sets_Global
         Write(IOBuffer, 202) i
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         If (AppCtx%MyEXO%NSProperty(Heat_NSProp_HasPForce)%Value(i) /= 0 ) Then
            Write(IOBuffer, 302) i, 'Temperature'
            Call MEF90_AskReal(U(i), IOBuffer, BatchUnit, IsBatch)
         End If
      End Do

202 Format('    Node Set      ', T24, I3, '\n')
300 Format('EB', I4.4, ': ', A) 
301 Format('EB', I4.4, ': ', A, I3)
302 Format('NS', I4.4, ': ', A)
      AppCtx%VertVar_Temperature = 1
! Set BC on NS
      Call HeatSetBC(AppCtx, U, AppCtx%MyExo, AppCtx%MeshTopology, VarFrac_NSProp_BCT)
      !!! Create the EXO case file
      Call Write_EXO_Case(AppCtx%AppParam%prefix, '%0.4d', MEF90_NumProcs)

   End Subroutine PoissonInit


#undef __FUNCT__
#define __FUNCT__ "HeatSetInitial"
   Subroutine HeatSetInitial(AppCtx, MeshTopology, ValU, ValF)
      Type(Heat_AppCtx_Type)                             :: AppCtx
      Type (MeshTopology_Type)                           :: MeshTopology 
      PetscReal, Dimension(:), Pointer                   :: ValU, ValF
      PetscReal, Dimension(:), Pointer                   :: Uelem, Felem
      PetscInt                                           :: iBlk, iEloc, Num_Dof, iE
      PetscInt                                           :: iErr, i

      Do_Elem_IBlk: Do iBlk = 1, MeshTopology%Num_Elem_Blks
         i = MeshTopology%Elem_Blk(iBlk)%ID
         Num_DoF = MeshTopology%Elem_Blk(iBlk)%Num_DoF
         Allocate(Uelem(Num_DoF))
         Allocate(Felem(Num_DoF))
         Uelem = ValU(i)
         Felem = ValF(i)
         Do_Elem_iE: Do iELoc = 1, MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            Call SectionRealUpdateClosure(AppCtx%T%Sec, MeshTopology%Mesh, iE-1, Uelem, INSERT_VALUES, ierr);CHKERRQ(iErr)
! attention les forces c est pour tous les pas de temps !!!! ici que premier  
            Call SectionRealUpdateClosure(AppCtx%F%Sec, MeshTopology%Mesh, iE-1, Felem, INSERT_VALUES, ierr);CHKERRQ(iErr) 
         End Do Do_Elem_iE 
         DeAllocate(Uelem)
         DeAllocate(Felem)
      End Do Do_Elem_Iblk
   End Subroutine HeatSetInitial

#undef __FUNCT__
#define __FUNCT__ "HeatSetBC"
   Subroutine HeatSetBC(AppCtx, U, MyExo, MeshTopology, NS_Offset)
      Type(Heat_AppCtx_Type)                             :: AppCtx
      Type (MeshTopology_Type)                           :: MeshTopology 
      Type (EXO_Type)                                    :: MyEXO
      PetscReal, Dimension(:), Pointer                   :: U, Uelem
      PetscInt                                           :: NS_Offset 
      PetscInt                                           :: iErr, iLoc, i, j  

      Allocate(Uelem(1))
      Do iloc = 1, MeshTopology%Num_Node_Sets         
         i = MeshTopology%Node_Set(iloc)%ID
         If (MyEXO%NSProperty(NS_Offset)%Value( MeshTopology%Node_Set(iloc)%ID)== 1) then 
            Uelem(1) =  U(i) !Dirichlet BC 
            Do j = 1, MeshTopology%Node_Set(iloc)%Num_Nodes
               Call SectionRealUpdate(AppCtx%T%Sec, MeshTopology%Num_Elems + MeshTopology%Node_Set(iloc)%Node_ID(j)-1, Uelem, INSERT_VALUES, iErr); CHKERRQ(iErr)
            End Do
         End If 
      End Do
      DeAllocate(Uelem)
!TODO write for all timesteps
      Call Write_EXO_Result_Vertex(MyEXO, MeshTopology,  AppCtx%VertVar_Temperature, 1, AppCtx%T%Sec)
      Call SectionIntAddNSProperty(AppCtx%BCFlag%Sec,  MyEXO%NSProperty(NS_Offset),  MeshTopology)
   End Subroutine HeatSetBC


#undef __FUNCT__
#define __FUNCT__ "HeatFieldCreate"
   Subroutine HeatFieldCreate(AppCtx, MeshTopology)
      Type(Heat_AppCtx_Type)                             :: AppCtx
      Type (MeshTopology_Type)                     :: MeshTopology
      
      PetscInt, Dimension(:), Pointer                    :: SizeScal
      
      !!! Allocate the Section for U and F
      Allocate(SizeScal(1))
      SizeScal=1
      Call FieldCreateVertex(AppCtx%T,     'T',         MeshTopology, SizeScal)
      Call FieldCreateVertex(AppCtx%F,     'F',         MeshTopology, SizeScal)
      Call FieldCreateVertex(AppCtx%RHS,   'RHS',       MeshTopology, SizeScal)
      Call FlagCreateVertex(AppCtx%BCFlag, 'BC',        MeshTopology, SizeScal)
      Call FieldCreateVertex(AppCtx%TBC,    'TBC',      MeshTopology, SizeScal)
      DeAllocate(SizeScal)
      
   End Subroutine HeatFieldCreate

#undef __FUNCT__
#define __FUNCT__ "Poisson_TSSetUp"
   Subroutine Poisson_TSSetUp(AppCtx, MeshTopology)
      Type(Heat_AppCtx_Type)                            :: AppCtx
      Type (MeshTopology_Type)                     :: MeshTopology

      Type(SectionReal)                            :: TmpSec
      Character(len=256)                           :: TmpSecName
      PetscInt                                     :: iErr

      !!! Initialize the matrix and vector for the linear system
      Call DMMeshSetMaxDof(MeshTopology%Mesh, MeshTopology%Num_Dim, iErr); CHKERRQ(iErr) 
      Call DMMeshCreateMatrix(MeshTopology%mesh, AppCtx%T%sec, MATMPIAIJ, AppCtx%K, iErr); CHKERRQ(iErr)
      Call PetscObjectSetName(AppCtx%K,"K matrix",iErr);CHKERRQ(iErr)
      Call DMMeshCreateMatrix(MeshTopology%mesh, AppCtx%T%sec, MATMPIAIJ, AppCtx%M, iErr); CHKERRQ(iErr)
      Call PetscObjectSetName(AppCtx%M,"M matrix",iErr);CHKERRQ(iErr)
      Call DMMeshCreateMatrix(MeshTopology%mesh, AppCtx%T%sec, MATMPIAIJ, AppCtx%Jac, iErr); CHKERRQ(iErr)
      Call PetscObjectSetName(AppCtx%Jac,"Jac matrix",iErr);CHKERRQ(iErr)
      
      Call TSCreate(PETSC_COMM_WORLD, AppCtx%TS, iErr); CHKERRQ(iErr)

!!! Could be simpler with SectionRealDuplicate (when Fortran Binding exists)
      Call SectionRealDuplicate(AppCtx%T%Sec,TmpSec,iErr);CHKERRQ(iErr)
      Call PetscObjectSetName(TmpSec,"default",iErr);CHKERRQ(iErr)
      call DMMeshSetSectionReal(MeshTopology%mesh,trim('default'),TmpSec,iErr);CHKERRQ(iErr)
      !Destroy section tmpsec!!?
      
      call TSSetDM(AppCtx%TS, MeshTopology%mesh, ierr);   CHKERRQ(iErr) 
      Call TSSetIFunction(AppCtx%TS, PETSC_NULL_OBJECT, IFunctionPoisson, AppCtx, ierr); CHKERRQ(iErr)
      Call TSSetIJacobian(AppCtx%TS, AppCtx%Jac, AppCtx%Jac, IJacobianPoisson, AppCtx, ierr); CHKERRQ(iErr)
      Call TSSetRHSFunction(AppCtx%TS, PETSC_NULL_OBJECT, RHSPoisson, AppCtx, ierr); CHKERRQ(iErr)

! Setting time discretization scheme for TS
      Call TSSetType(AppCtx%TS, AppCtx%HeatSchemeParam%HeatTsType, iErr); CHKERRQ(iErr)
      Call TSRosWSetType(AppCtx%TS, AppCtx%HeatSchemeParam%TsTypeOpt1, iErr); CHKERRQ(iErr)
      Call TSSetFromOptions(AppCtx%TS,  iErr); CHKERRQ(iErr) ! overwritting if cli arguments
       
   End Subroutine Poisson_TSSetUp


#undef __FUNCT__
#define __FUNCT__ "IFunctionPoisson"
   SubRoutine IFunctionPoisson(dummyTS, t, U, Udot, GlobalOut, AppCtx, iErr)
      PetscReal                                    :: t 
      Type(Vec)                                    :: U, Udot, GlobalOut
      Type(Heat_AppCtx_Type)                            :: AppCtx
      Type(TS)                                     :: dummyTS
      PetscInt                                     :: iErr
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer   
      If (AppCtx%AppParam%verbose > 0) Then    
         Write(IOBuffer, *) 'Begin IFunctionPoisson \n'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If 
      
!GlobalOut = AppCtx%M * Udot + AppCtx%K * U
      Call MatMult(AppCtx%K, U, GlobalOut, iErr);CHKERRQ(iErr)
      Call MatMultAdd(AppCtx%M, Udot, GlobalOut, GlobalOut, iErr);CHKERRQ(iErr)
       
   End Subroutine IFunctionPoisson


#undef __FUNCT__
#define __FUNCT__ "IJacobianPoisson"
   SubRoutine IJacobianPoisson(dummyTS, t, U, Udot, a, Jac, PreJac, mStruc, AppCtx, iErr)
      Type(TS)                                     :: dummyTS
      PetscReal                                    :: t, a 
      Type(Vec)                                    :: U, Udot
      Type(Mat)                                    :: Jac, PreJac  
      Type(Heat_AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer  
      MatStructure                           :: mStruc
! a*AppCtx%M + AppCtx%K
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Begin IJacobianPoisson \n'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End IF

      Call MatCopy(AppCtx%K, Jac, mStruc, iErr); CHKERRQ(iErr)
      Call MatAXPY(Jac, a, AppCtx%M, mStruc, iErr); CHKERRQ(iErr)

   End Subroutine IJacobianPoisson
 
#undef __FUNCT__
#define __FUNCT__ "RHSPoisson"
   SubRoutine RHSPoisson(dummyTS, t, U, GlobalOut, AppCtx, iErr)
      Type(TS)                                     :: dummyTS
      PetscReal                                    :: t 
      Type(Vec)                                    :: U, GlobalOut
      Type(Heat_AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer   
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Begin RHSPoisson \n'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If

      Call VecCopy(AppCtx%RHS%Vec, GlobalOut, iErr); CHKERRQ(iErr)
     
   End Subroutine RHSPoisson
   
   
#undef __FUNCT__
#define __FUNCT__ "HeatMatAssembly"
   Subroutine HeatMatAssembly(AppCtx, MeshTopology, ExtraField)   
      Type(Heat_AppCtx_Type)                             :: AppCtx
      Type (MeshTopology_Type)                           :: MeshTopology
      PetscInt                                           :: iBlk, iErr
      Type(Field)                                        :: ExtraField
      
      Call MatInsertVertexBoundaryValues(AppCtx%K, AppCtx%T, AppCtx%BCFlag, MeshTopology)
      Call MatAssemblyBegin(AppCtx%K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (AppCtx%K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)

      Do_iBlk: Do iBlk = 1, MeshTopology%Num_Elem_Blks
         Call HeatMatAssemblyBlock(iBlk, AppCtx, MeshTopology, ExtraField)
      End Do Do_iBlk
      Call MatAssemblyBegin(AppCtx%K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (AppCtx%K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)

   End Subroutine HeatMatAssembly
      
   Subroutine HeatMatAssemblyBlock(iBlk, AppCtx, MeshTopology, ExtraField)
      Type(Heat_AppCtx_Type)                       :: AppCtx
      Type (MeshTopology_Type)                     :: MeshTopology
      PetscInt                                     :: iBlk
      Type(Field)                                  :: ExtraField
      
      PetscInt                                     :: iE, iELoc, iErr, i
      PetscReal, Dimension(:,:), Pointer           :: MatElem
      PetscInt, Dimension(:), Pointer              :: BCFlag
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      PetscInt                                     :: NumDoFScal
      PetscLogDouble                               :: flops = 0
      PetscReal                                    :: lDiff
      PetscReal, Dimension(:), Pointer             :: T_Loc, Extra_Loc
      PetscReal                                    :: T_Elem, Extra_Elem
     
      NumDoFScal = MeshTopology%Elem_Blk(iBlk)%Num_DoF
      Allocate(MatElem(MeshTopology%Elem_Blk(iBlk)%Num_DoF, MeshTopology%Elem_Blk(iBlk)%Num_DoF))
      Allocate(BCFlag(MeshTopology%Elem_Blk(iBlk)%Num_DoF))
      Allocate(T_Loc(NumDoFScal))
      Allocate(Extra_Loc(NumDoFScal))

      Do_iELoc: Do iELoc = 1, MeshTopology%Elem_Blk(iBlk)%Num_Elems
         i = MeshTopology%Elem_Blk(iBlk)%ID
         iE = MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         MatElem = 0.0_Kr
         BCFlag = 0
         Call SectionIntRestrictClosure(AppCtx%BCFlag%Sec,  MeshTopology%mesh, iE-1, NumDoFScal, BCFlag,    iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(AppCtx%T%Sec,      MeshTopology%mesh, iE-1, NumDoFScal, T_Loc,     iErr); CHKERRQ(ierr)
         Do iGauss = 1, size(AppCtx%Elem(iE)%Gauss_C)
         Select Case(AppCtx%MatProp(i)%Type_Law)
            Case(Heat_Constant_ID)
            ! Constant Diffusion
               lDiff = AppCtx%MatProp(i)%Diffusivity(1) 
            Case(Heat_Increasing_ID )
            !Increasing diffusion 
               ! Mensi Law D(C) = A exp (B*C) 
               ! 1988 Mensi-Acker-Attolou Mater Struct
               ! C : water content,  A and B material parameters
               T_Elem = 0.0_Kr
               Do iDoF1 = 1, NumDoFScal
                  T_Elem = T_Elem + AppCtx%Elem(iE)%BF(iDoF1, iGauss) * T_Loc(iDoF1)
               End DO
               lDiff = AppCtx%MatProp(i)%Diffusivity(1)*exp(AppCtx%MatProp(i)%Diffusivity(2)*T_Elem) 
            Case(Heat_Decreasing_ID)
      !Diffusion is decreasing with the variable
               T_Elem = 0.0_Kr
               Do iDoF1 = 1, NumDoFScal
                  T_Elem = T_Elem + AppCtx%Elem(iE)%BF(iDoF1, iGauss) * T_Loc(iDoF1)
               End DO
               lDiff =  AppCtx%MatProp(i)%Diffusivity(1)*(AppCtx%MatProp(i)%Diffusivity(2)-T_Elem)**AppCtx%MatProp(i)%Diffusivity(3) 
            Case(Heat_Non_Monotonic_ID)
      !Diffusion is non monotonic with the variable
               T_Elem = 0.0_Kr
               Do iDoF1 = 1, NumDoFScal
                  T_Elem = T_Elem + AppCtx%Elem(iE)%BF(iDoF1, iGauss) * T_Loc(iDoF1)
               End DO
      !2007 Cement Concrete Research_baroghel-bouny_I and II         
               lDiff = AppCtx%MatProp(i)%Diffusivity(1)*(AppCtx%MatProp(i)%Diffusivity(2)-T_Elem)**AppCtx%MatProp(i)%Diffusivity(3)+AppCtx%MatProp(i)%Diffusivity(4)*exp(AppCtx%MatProp(i)%Diffusivity(5)*T_Elem) 
           Case(Heat_Damage_ID)
      ! Diffusion Depends on the damaging variable
               Call SectionRealRestrictClosure(ExtraField%Sec, MeshTopology%mesh,  iE-1, NumDoFScal, Extra_Loc, iErr); CHKERRQ(ierr)
               Extra_Elem = 0.001_Kr
               Do iDoF1 = 1, NumDoFScal
                  Extra_Elem = Extra_Elem + AppCtx%Elem(iE)%BF(iDoF1, iGauss) * Extra_Loc(iDoF1)
               End DO
               lDiff = min(AppCtx%MatProp(i)%Diffusivity(1)/(Extra_Elem+0.001), AppCtx%MatProp(i)%Diffusivity(2))
      ! TODO : Test that ExtraField is defined to evoid Segment Fault error
            Case(Heat_Crack_Open_ID)
      ! Diffusion depends on crack opening
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
         Call DMMeshAssembleMatrix(AppCtx%K, MeshTopology%mesh, AppCtx%T%Sec, iE-1, MatElem, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iELoc
   
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
      DeAllocate(BCFlag)
      DeAllocate(MatElem)
   End Subroutine HeatMatAssemblyBlock


#undef __FUNCT__
#define __FUNCT__ "RHSAssembly"
   Subroutine RHSAssembly(AppCtx, MeshTopology, MyExo)
      Type(Heat_AppCtx_Type)                       :: AppCtx
      Type (MeshTopology_Type)                     :: MeshTopology
      Type (EXO_Type)                              :: MyEXO

      PetscInt                                     :: iErr
      PetscInt                                     :: iBlk

      Call SectionRealZero(AppCtx%RHS%Sec, iErr); CHKERRQ(iErr)
      
      Do_iBlk: Do iBlk = 1, MeshTopology%Num_Elem_Blks
         Call RHSAssemblyBlock(iBlk, AppCtx, MeshTopology)
      End Do Do_iBlk

      Call SectionRealComplete(AppCtx%RHS%Sec, iErr); CHKERRQ(iErr)

      !!! Set Dirichlet Boundary Values
!Suppose that loading is contant (replace second to last by AppCtx%Timestep otherwise)
      Call Read_EXO_Result_Vertex(MyEXO, MeshTopology, AppCtx%VertVar_Temperature, 1, AppCtx%TBC)
      Call FieldInsertVertexBoundaryValues(AppCtx%RHS, AppCtx%TBC, AppCtx%BCFlag, MeshTopology)

      Call SectionRealToVec(AppCtx%RHS%Sec, AppCtx%RHS%Scatter, SCATTER_FORWARD, AppCtx%RHS%Vec, iErr); CHKERRQ(iErr)
      Call FieldInsertVertexBoundaryValues(AppCtx%RHS, AppCtx%TBC, AppCtx%BCFlag, MeshTopology)
      !!! VERY important! This is the equivalent of a ghost update
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
      
      flops = 0.0

      Num_DoF = MeshTopology%Elem_Blk(iBlk)%Num_DoF
      Allocate(F_Loc(Num_DoF))
      Allocate(RHS_Loc(Num_DoF))
      Allocate(BCFlag_Loc(Num_DoF))

      iBlkID = MeshTopology%Elem_Blk(iBlk)%ID
      Do_iEloc: Do iELoc = 1, MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
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
   End Subroutine RHSAssemblyBlock
#undef __FUNCT__
#define __FUNCT__ "SolveTransientStep"
   Subroutine SolveTransientStep(AppCtx, MyEXO, MeshTopology, TimeStepIni, TimeStepFinal, iStep)
      Type(Heat_AppCtx_Type)                            :: AppCtx
      Type (MeshTopology_Type)                     :: MeshTopology
      Type(EXO_Type)                               :: MyEXO
      PetscInt                                     :: TSTimeSteps, iErr, iStep
      TSConvergedReason                            :: TSreason 
      PetscReal                                    :: TimeStepIni, TimeStepFinal
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
     
      Write(IOBuffer, *) "Warning : TSSolve does not assure that computation stops as the exact asked time \n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr);   CHKERRQ(iErr)
!TODO see TSSetExactFinalTime, TSGetTimeStep
!TSSetExactFinalTime can not be used : 'TSRosW 2p does not have an interpolation formula'
!      Call TSSetExactFinalTime(AppCtx%TS, PETSC_TRUE, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(AppCtx%T%Sec, AppCtx%T%Scatter, SCATTER_FORWARD, AppCtx%T%Vec, iErr); CHKERRQ(iErr)

! Using IMEX see section 6.1.2 PetSc documentation 
      if (AppCtx%AppParam%verbose > 1) Then                                                                                                                                                
         Call TSView(AppCtx%TS,  PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
      End If

      Write(IOBuffer, 400) iStep, TimeStepIni, TimeStepFinal 
      Call FieldInsertVertexBoundaryValues(AppCtx%T, AppCtx%TBC, AppCtx%BCFlag, MeshTopology)
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr);   CHKERRQ(iErr)
      Call TSSetSolution(AppCtx%TS, AppCtx%T%Vec,  ierr);   CHKERRQ(iErr)

      If (AppCtx%AppParam%verbose > 3) Then
         Call VecView(AppCtx%T%Vec,  PETSC_VIEWER_STDOUT_WORLD, iErr)
         Call VecView(AppCtx%TBC%Vec,  PETSC_VIEWER_STDOUT_WORLD, iErr)
      End If

      Call TSSetInitialTimeStep(AppCtx%TS, TimeStepIni, (TimeStepFinal - TimeStepIni)/100.,  ierr); CHKERRQ(iErr)
      Call TSSetDuration(AppCtx%TS, AppCtx%HeatSchemeParam%HeatMaxiter, TimeStepFinal, iErr); CHKERRQ(iErr)
      Call TSSolve(AppCtx%TS, AppCtx%T%Vec, TimeStepFinal, iErr); CHKERRQ(iErr)

      Call SectionRealToVec(AppCtx%T%Sec, AppCtx%T%Scatter, SCATTER_REVERSE, AppCtx%T%Vec, ierr); CHKERRQ(ierr)
      Call Write_EXO_Result_Vertex(MyEXO, MeshTopology, AppCtx%VertVar_Temperature , iStep+1, AppCtx%T%Sec)
      Call TSGetTimeStepNumber(AppCtx%TS, TSTimeSteps, iErr); CHKERRQ(iErr) 
      Call TSGetConvergedReason(AppCtx%TS, TSreason, iErr); CHKERRQ(iErr)
      Write(IOBuffer, 400) iStep, TimeStepIni, TimeStepFinal 
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, 100) TSTimeSteps, TSreason
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      
100 Format('TS', I5, ' TimeSteps. TSConvergedReason is ', I2, '\n')
400 Format('Solving TS Time step : ', I5, '    Ini Time : ', ES12.5,  ",    Final Time  :", ES12.5, '\n')
   End Subroutine SolveTransientStep

     
#undef __FUNCT__
#define __FUNCT__ "MatMassAssembly"
   Subroutine MatMassAssembly(AppCtx, MeshTopology)   
      Type(Heat_AppCtx_Type)                            :: AppCtx
      Type (MeshTopology_Type)                     :: MeshTopology
      PetscInt                                     :: iBlk, iErr
      !Call PetscLogStagePush(AppCtx%LogInfo%MatAssembly_Stage, iErr); CHKERRQ(iErr)
      Do_iBlk: Do iBlk = 1, MeshTopology%Num_Elem_Blks
         Call MatMassAssemblyBlock(iBlk, AppCtx, MeshTopology)
      End Do Do_iBlk
      Call MatAssemblyBegin(AppCtx%M, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (AppCtx%M, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)

   End Subroutine MatMassAssembly
      
#undef __FUNCT__
#define __FUNCT__ "MatMassAssemblyBlock"
   Subroutine MatMassAssemblyBlock(iBlk, AppCtx, MeshTopology)
      Type(Heat_AppCtx_Type)                            :: AppCtx
      Type (MeshTopology_Type)                     :: MeshTopology
      PetscInt                                     :: iBlk
      
      PetscInt                                     :: iE, iELoc, iErr
      PetscReal, Dimension(:,:), Pointer           :: MatElem
      PetscInt, Dimension(:), Pointer              :: BCFlag
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      PetscLogDouble                               :: flops = 0
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
      
      Allocate(MatElem(MeshTopology%Elem_Blk(iBlk)%Num_DoF, MeshTopology%Elem_Blk(iBlk)%Num_DoF))
      Allocate(BCFlag(MeshTopology%Elem_Blk(iBlk)%Num_DoF))

      Do_iELoc: Do iELoc = 1, MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         MatElem = 0.0_Kr
         BCFlag = 0 
         Call SectionIntRestrictClosure(AppCtx%BCFlag%Sec, MeshTopology%mesh, iE-1, MeshTopology%Elem_Blk(iBlk)%Num_DoF, BCFlag, iErr); CHKERRQ(ierr)
         Do iGauss = 1, size(AppCtx%Elem(iE)%Gauss_C)
            Do iDoF1 = 1, MeshTopology%Elem_Blk(iBlk)%Num_DoF
             !     MatElem(iDoF1, iDoF1) = 1
          !  If (BCFlag(iDoF1) == 0) Then
               Do iDoF2 = 1, MeshTopology%Elem_Blk(iBlk)%Num_DoF
                  MatElem(iDoF2, iDoF1) = MatElem(iDoF2, iDoF1) + AppCtx%Elem(iE)%Gauss_C(iGauss) * ( AppCtx%Elem(iE)%BF(iDoF1, iGauss) *  AppCtx%Elem(iE)%BF(iDoF2, iGauss) )
                  flops = flops + 2 !Check 
               End Do
         !  End If
            End Do
         End Do
         Call DMMeshAssembleMatrix(AppCtx%M, MeshTopology%mesh, AppCtx%T%sec, iE-1, MatElem, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iELoc
   
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
      DeAllocate(MatElem)
   End Subroutine MatMassAssemblyBlock

 
#undef __FUNCT__
#define __FUNCT__ "Heat_Field_Mat_Finalize"
   Subroutine Heat_Field_Mat_Finalize(AppCtx)
      Type(Heat_AppCtx_Type)                       :: AppCtx
      PetscInt                                     :: iErr

      Call FieldDestroy(AppCtx%T)
      Call FieldDestroy(AppCtx%F)
      Call FieldDestroy(AppCtx%RHS)
      Call FieldDestroy(AppCtx%TBC)
      Call SectionIntDestroy(AppCtx%BCFlag, iErr); CHKERRQ(iErr)
      Call MatDestroy(AppCtx%K, iErr); CHKERRQ(iErr)
      Call MatDestroy(AppCtx%M, iErr); CHKERRQ(iErr)
      Call MatDestroy(AppCtx%Jac, iErr); CHKERRQ(iErr)

      Call TSDestroy(AppCtx%TS, iErr); CHKERRQ(iErr)

   End Subroutine Heat_Field_Mat_Finalize 

#undef __FUNCT__
#define __FUNCT__ "TSPoissonFinalize"
   Subroutine TSPoissonFinalize(AppCtx)   
      Type(Heat_AppCtx_Type)                            :: AppCtx

      PetscInt                                     :: iErr
      Character(len=MEF90_MXSTRLEN)                :: filename

      Call Heat_Field_Mat_Finalize(AppCtx)
      Call SectionRealDestroy(AppCtx%GradU, iErr); CHKERRQ(iErr)

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
            Call SectionRealRestrictClosure(AppCtx%T%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoF, U, iErr); CHKERRQ(ierr)
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
   End Subroutine ComputeEnergy

   
#undef __FUNCT__
#define __FUNCT__ "ComputeWaterMass"
   Subroutine ComputeWaterMass(MeshTopology, Exo, Id_Vertex, Elem, NumTimeSteps)
      Type(SectionReal)                            :: V_Sec_0, V_Sec_i
      Type(MeshTopology_Type)                      :: MeshTopology
      Type (EXO_Type)                              :: EXO
#if defined PB_2D
      Type(Element2D_Scal), Dimension(:), Pointer  :: Elem
#elif defined PB_3D 
      Type(Element3D_Scal), Dimension(:), Pointer  :: Elem
#endif   
      PetscInt                                     :: iErr, iStep, iBlk 
      PetscInt                                     :: NumTimeSteps
      PetscInt                                     :: MassFileUnit
      PetscScalar                                  :: neg
      PetscReal, Dimension(:), Pointer             :: WaterMass
      PetscInt                                     :: Id_Vertex
!Add a density material parameter      
      Allocate(WaterMass(MeshTopology%Num_Elem_Blks_Global))
      Call DMMeshGetVertexSectionReal(MeshTopology%mesh, 'V_Sec_0', 1, V_Sec_0, iErr); CHKERRQ(iErr)
      Call DMMeshGetVertexSectionReal(MeshTopology%mesh, 'V_Sec_i', 1, V_Sec_i, iErr); CHKERRQ(iErr)
      Call Read_EXO_Result_Vertex(EXO, MeshTopology, Id_Vertex, 1, V_Sec_0)
      neg = -1.
      MassFileUnit = 80
      WaterMass = 0.0_Kr
      Open(File = "water.mass", Unit = MassFileUnit, Status = 'Unknown')
      Rewind(MassFileUnit)
      Write(MassFileUnit, 400, advance='no') 0
      Call IntegrateScalLp(MeshTopology,Elem,V_Sec_0,WaterMass,1)
      Write(MassFileUnit, *)  WaterMass
      Do istep = 1, NumTimeSteps-1
         Write(MassFileUnit, 400, advance='no') istep
         Call Read_EXO_Result_Vertex(EXO, MeshTopology, Id_Vertex, 1, V_Sec_0)
         Call Read_EXO_Result_Vertex(EXO, MeshTopology, Id_Vertex, iStep+1, V_Sec_i)
         Call SectionRealAXPY(V_Sec_0, MeshTopology%mesh, neg, V_Sec_i, iErr); CHKERRQ(iErr)
         Call IntegrateScalLp(MeshTopology,Elem,V_Sec_i,WaterMass,1)
         Write(MassFileUnit, *)  WaterMass
      End DO 
400 Format(I6) 
!401 Format(ES13.5, '   ') 
      Close(MassFileUnit)
      DeAllocate(WaterMass)
   End Subroutine ComputeWaterMass

#if defined PB_2D
End Module m_TransientHeat2D
#elif defined PB_3D
End Module m_TransientHeat3D
#endif
