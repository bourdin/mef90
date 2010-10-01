Program PrepVarFrac

#include "finclude/petscdef.h"

   Use m_MEF90
   Use m_VarFrac_Struct
   Use m_PrepVarFrac
   Use petsc

   Implicit NONE   
   
   Type TestCase_Type
      PetscInt                                  :: Index
      Character(len=MEF90_MXSTRLEN)             :: Description
   End Type
   
   
   Character(len=MEF90_MXSTRLEN)                :: prefix, IOBuffer, filename
   Type(MeshTopology_Type)                      :: MeshTopology, GlobalMeshTopology
   Type(EXO_Type)                               :: EXO, MyEXO
   Type(Mesh)                                   :: Tmp_Mesh
   Type(Element2D_Scal), Dimension(:), Pointer  :: Elem2D
   Type(Element3D_Scal), Dimension(:), Pointer  :: Elem3D
   PetscBool                                    :: HasPrefix
   PetscInt                                     :: verbose = 0
   PetscInt                                     :: iErr, iloc, i, j, k, iCase, iStep
   PetscInt                                     :: iBlock

   PetscInt                                     :: NumTestCase
   Type(TestCase_Type), Dimension(:) , Pointer  :: TestCase
   PetscInt, Parameter                          :: QuadOrder=2
   Type(SectionReal)                            :: USec, FSec, VSec, ThetaSec, CoordSec
   Type(Vect3D), Dimension(:), Pointer          :: U, F
   Type(Vect3D)                                 :: FGrain
   PetscReal, Dimension(:), Pointer             :: V, Theta
   PetscReal                                    :: ThetaGrain
   PetscReal, Dimension(:), Pointer             :: Uelem, Felem, Velem, Thetaelem, Coordelem
   PetscReal                                    :: Tmin, Tmax
   PetscReal, Dimension(:), Pointer             :: T
   PetscInt                                     :: NumSteps
   PetscInt                                     :: Num_DoF
   PetscReal                                    :: RealBuffer
   Type(MatProp2D_Type), Dimension(:), Pointer  :: MatProp2D
   Type(MatProp3D_Type), Dimension(:), Pointer  :: MatProp3D
   Type(Tens4OS2D)                              :: TmpHL_2D
   Type(Tens4OS3D)                              :: TmpHL_3D   
   PetscReal                                    :: E, nu, Toughness, Therm_ExpScal
   PetscReal, Dimension(:), Pointer             :: GlobVars
   PetscInt                                     :: vers
   PetscBool                                    :: IsBatch, HasBatchFile
   PetscInt                                     :: BatchUnit=99
   Character(len=MEF90_MXSTRLEN)                :: BatchFileName
   PetscBool                                    :: EraseBatch
   
   PetscInt                                     :: NumGrains
   PetscReal                                    :: ToughnessGrain, ThermExpScalGrain
   PetscReal                                    :: BGrain, CGrain, CpGrain
   PetscReal                                    :: alphaGrain, alphamin, alphamax
   PetscReal                                    :: Eeff, nueff, Kappa, mu
   PetscRandom                                  :: RandomCtx
   PetscInt                                     :: Seed
   PetscLogDouble                               :: Time
   PetscBool                                    :: Has_Seed, Has_n
   
   PetscReal                                    :: R, Ctheta, CTheta2, Stheta, STheta2

   Call MEF90_Initialize()

   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-verbose', verbose, HasPrefix, iErr)    
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p', prefix, HasPrefix, iErr); CHKERRQ(iErr)
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No mesh prefix given\n", iErr)
      Call MEF90_Finalize()
      STOP
   End If
   EraseBatch=.False.
   Call PetscOptionsGetTruth(PETSC_NULL_CHARACTER, '-force', EraseBatch, HasPrefix, iErr)    
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-i', BatchFileName, IsBatch, iErr); CHKERRQ(iErr)
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
   
   !!! Initialize random number generator
   Call PetscRandomCreate(PETSC_COMM_WORLD, RandomCtx, iErr); CHKERRQ(iErr)
   Call PetscRandomSetFromOptions(RandomCtx, iErr); CHKERRQ(iErr)
   Call PetscOptionsGetReal(PETSC_NULL_CHARACTER, '-seed', Seed, Has_Seed, iErr); CHKERRQ(iErr)
   If (.NOT. Has_Seed) Then
      Call PetscGetTime(Time, iErr); CHKERRQ(iErr)
      Seed =  Time * (Time - Int(Time))
   End If
   !Write(IOBuffer, *) 'Seeding random number generator with seed ', Seed, '\n'
   !Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call PetscRandomSetSeed(RandomCtx, Seed, iErr); CHKERRQ(iErr)
   Call PetscRandomSeed(RandomCtx, iErr); CHKERRQ(iErr)
   Call PetscRandomGetSeed(RandomCtx, Seed, iErr); CHKERRQ(iErr)
   NumTestCase = 4
   Allocate(TestCase(NumTestCase))
   Do i = 1, NumTestCase
      TestCase(i)%Index = i
   End Do
   TestCase(1)%Description = "Polycrystal, MIL"
   TestCase(2)%Description = "Polycrystal, Mode-I scaling experiment"
   TestCase(3)%Description = "Polycrystal, Mode-I surfing experiment"
   TestCase(4)%Description = "Polycrystal, Mode-I isotropic scaling experiment"
   

   Call Write_EXO_Case(prefix, '%0.4d', MEF90_NumProcs)
   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'
   
   Call EXO_Check_Numbering(EXO, iErr)
   If (iErr /= 0) Then
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, 'Unsupported numbering of the element blocks, side sets or node sets\n', iErr)
   End If

   !!! Reading and distributing sequential mesh
   If (MEF90_NumProcs == 1) Then
      Call MeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, MeshTopology%mesh, ierr); CHKERRQ(iErr)
   Else
      Call MeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
      Call MeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, MeshTopology%mesh, ierr); CHKERRQ(iErr)
      Call MeshDestroy(Tmp_mesh, ierr); CHKERRQ(iErr)
   End If

   Call MeshTopologyReadEXO(MeshTopology, EXO)
   If (verbose > 0) Then
      Write(IOBuffer, *) "Done reading and partitioning the mesh\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If

   MyEXO%comm = PETSC_COMM_SELF
   MyEXO%exoid = EXO%exoid
   Write(MyEXO%filename, 99) trim(prefix), MEF90_MyRank
 99  Format(A, '-', I4.4, '.gen')
   Call MeshGetSectionReal(MeshTopology%mesh, 'coordinates', CoordSec, iErr); CHKERRQ(ierr)
   
   Call VarFracEXOProperty_Init(MyEXO, MeshTopology)   
   If (verbose > 0) Then
      Write(IOBuffer, *) "Done with VarFracEXOProperty_Init\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   Write(IOBuffer, *) '\nElement Block and Node Set Properties\n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)      

   !!! Get number of grains from command line
   numGrains = MeshTopology%Num_Elem_Blks_Global
   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-n', NumGrains, has_n, iErr) 
   If (NumGrains < 0) Then
      NumGrains = MeshTopology%Num_Elem_Blks_Global + NumGrains
   End If
   If (NumGrains > MeshTopology%Num_Elem_Blks_Global) Then
      Write(IOBuffer, *) "Number of grains cannot be greater than number of blocks! ", NumGrains, MeshTopology%Num_Elem_Blks_Global, "\n"
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ, IOBuffer, iErr)
   End If
   Call EXOEBProperty_AskWithBatchGrains(MyEXO, MeshTopology, BatchUnit, IsBatch, NumGrains)
   Call EXOSSProperty_AskWithBatchGrains(MyEXO, MeshTopology, BatchUnit, IsBatch, NumGrains)
   Call EXONSProperty_AskWithBatchGrains(MyEXO, MeshTopology, BatchUnit, IsBatch, NumGrains)
   
   Do i = 1, MeshTopology%num_elem_blks
      MeshTopology%elem_blk(i)%Elem_Type = MyEXO%EBProperty( VarFrac_EBProp_Elem_Type )%Value( MeshTopology%elem_blk(i)%ID )
      Call Init_Elem_Blk_Type(MeshTopology%Elem_Blk(i), MeshTopology%num_dim)
   End Do

   Call Write_MeshTopologyGlobal(MeshTopology, MyEXO, PETSC_COMM_WORLD)
   If (verbose > 0) Then
      Write(IOBuffer, *) "Done with Write_MeshTopologyGlobal\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If

   Call EXOProperty_Write(MyEXO)
   If (verbose > 0) Then
      Write(IOBuffer, '(A)') 'Done with EXOProperty_Write\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If

   Call VarFracEXOVariable_Init(MyEXO)
   Call EXOVariable_Write(MyEXO)
   If (verbose > 0) Then
      Write(IOBuffer, '(A)') 'Done with EXOVariable_Write\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If


!!! Now we are done with the geometry, the variable definition and properties, so we can pre-compute the loads, displacement, time steps etc
   Select Case (MeshTopology%Num_Dim)
   Case (2) 
      Call ElementInit(MeshTopology, Elem2D, QuadOrder)
   Case(3)
      Call ElementInit(MeshTopology, Elem3D, QuadOrder)   
   Case Default
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, 'Only 2 and 3 dimensional elements are supported', iErr)
   End Select
   If (verbose > 0) Then
      Write(IOBuffer, '(A)') 'Done with ElementInit\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
!!! List all test cases and wait for case number 
   Do i = 1, NumTestCase
      Write(IOBuffer, "('[',I2.2,'] ',A)"), TestCase(i)%Index, Trim(TestCase(i)%Description)//'\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End Do
   
   Call AskInt(iCase, 'Test Case', BatchUnit, IsBatch)

!!!
!!! Global Variables: Time Steps and Load
!!!
   Write(IOBuffer, *) '\nGlobal Variables\n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)      
   Call AskReal(TMin, 'TMin', BatchUnit, IsBatch)
   Call AskReal(TMax, 'TMax', BatchUnit, IsBatch)
   Call AskInt(NumSteps, 'NumSteps', BatchUnit, IsBatch)

   Allocate(GlobVars(VarFrac_Num_GlobVar))
   GlobVars = 0.0_Kr
   Allocate(T(NumSteps))
   Select Case(iCase)
      Case Default
      Do i = 1, NumSteps-1
         T(i) = Tmin + Real(i-1) * (Tmax-TMin)/Real(NumSteps-1)
         GlobVars(VarFrac_GlobVar_Load) = T(i)
         Call Write_EXO_AllResult_Global(MyEXO, i, GlobVars)
   
         MyEXO%exoid = EXOPEN(MyEXO%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers, iErr)
         Call EXPTIM(MyEXO%exoid, i, T(i), iErr)
         Call EXCLOS(MyEXO%exoid, iErr)
         MyEXO%exoid = 0
      End Do
   
      T(NumSteps) = Tmax
      GlobVars(VarFrac_GlobVar_Load) = T(NumSteps)
      Call Write_EXO_AllResult_Global(MyEXO, NumSteps, GlobVars)
      MyEXO%exoid = EXOPEN(MyEXO%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers, iErr)
      Call EXPTIM(MyEXO%exoid, NumSteps, T(NumSteps), iErr)
      Call EXCLOS(MyEXO%exoid, iErr)
      MyEXO%exoid = 0
   End Select
   DeAllocate(GlobVars)
   
!!!
!!! EB Properties
!!!
   Select Case(MeshTopology%Num_Dim)
   Case(2)
      Allocate(MatProp2D(MeshTopology%Num_Elem_Blks_Global))
   Case(3)
      Allocate(MatProp3D(MeshTopology%Num_Elem_Blks_Global))
   End Select
   
   Write(IOBuffer, *) '\nMaterial Properties\n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)      
   Select Case(iCase)
   !!! Write special cases here
   Case Default
      If (NumGrains > 0) Then 
         !!! Default behavior is that the grains are the first element block and share the same properties (but not orientation)
         Call AskReal(ToughnessGrain,    'Grains: toughness',   BatchUnit, IsBatch)
         Call AskReal(BGrain,            'Grains: Bulk modulus B',         BatchUnit, IsBatch)
         Call AskReal(CGrain,            'Grains: Shear modulus C',        BatchUnit, IsBatch)
         Call AskReal(CpGrain,           'Grains: Shear modulus Cp',       BatchUnit, IsBatch)
         Call AskReal(ThermExpScalGrain, 'Grains: therm. exp.', BatchUnit, IsBatch)
      End If
      alphamin = 0.0_Kr
      alphamax = PETSC_PI
      Call PetscRandomSetInterval(RandomCtx, alphamin, alphamax, iErr); CHKERRQ(iErr)
      Do iBlock = 1, NumGrains
         Select Case(MeshTopology%Num_Dim)
         Case(2)
            MatProp2D(iBlock)%Toughness = ToughnessGrain
            Call PetscRandomGetValue(RandomCtx, alphaGrain, iErr); CHKERRQ(iErr)
            Write(IOBuffer, 90) iBlock, alphaGrain * 180.0_kr / PETSC_PI
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
            Call GenHL_Cubic2DPlaneStress_Zener(BGrain, CGrain, CpGrain, alphaGrain, MatProp2D(iBlock)%Hookes_Law)
            MatProp2D(iBlock)%Therm_Exp    = 0.0_Kr
            MatProp2D(iBlock)%Therm_Exp%XX = ThermExpScalGrain
            MatProp2D(iBlock)%Therm_Exp%YY = ThermExpScalGrain
         Case(3)
            MatProp3D(iBlock)%Toughness = ToughnessGrain
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, 'Polycrystal only supported in 2D at this point\n', iErr)
            MatProp3D(iBlock)%Therm_Exp    = 0.0_Kr
            MatProp3D(iBlock)%Therm_Exp%XX = ThermExpScalGrain
            MatProp3D(iBlock)%Therm_Exp%YY = ThermExpScalGrain
            MatProp3D(iBlock)%Therm_Exp%ZZ = ThermExpScalGrain
         End Select 
      End Do
      !!! Non grain EB if necessary
      Do iBlock = NumGrains+1, MeshTopology%Num_Elem_Blks_Global
         Write(IOBuffer, 100) iBlock
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         
         Write(IOBuffer, 300) iBlock, 'Toughness'
         Call AskReal(Toughness, IOBuffer, BatchUnit, IsBatch)
         Write(IOBuffer, 300) iBlock, 'Young Modulus'
         Call AskReal(E, IOBuffer, BatchUnit, IsBatch)
         Write(IOBuffer, 300) iBlock, 'Poisson Ratio'
         Call AskReal(nu, IOBuffer, BatchUnit, IsBatch)
         Write(IOBuffer, 300) iBlock, 'Therm Exp'
         Call AskReal(Therm_ExpScal, IOBuffer, BatchUnit, IsBatch)

         Select Case(MeshTopology%Num_Dim)
         Case(2)
            MatProp2D(iBlock)%Toughness = Toughness
            Call GenHL_Iso2D_EnuPlaneStress(E, nu, MatProp2D(iBlock)%Hookes_Law)
            MatProp2D(iBlock)%Therm_Exp    = 0.0_Kr
            MatProp2D(iBlock)%Therm_Exp%XX = Therm_ExpScal
            MatProp2D(iBlock)%Therm_Exp%YY = Therm_ExpScal
         Case(3)
            MatProp3D(iBlock)%Toughness = Toughness
            Call GenHL_Iso3D_Enu(E, nu, MatProp3D(iBlock)%Hookes_Law)
            MatProp3D(iBlock)%Therm_Exp    = 0.0_Kr
            MatProp3D(iBlock)%Therm_Exp%XX = Therm_ExpScal
            MatProp3D(iBlock)%Therm_Exp%YY = Therm_ExpScal
            MatProp3D(iBlock)%Therm_Exp%ZZ = Therm_ExpScal
         End Select 
      End Do
   Case(4)
      If (NumGrains > 0) Then 
         !!! Default behavior is that the grains are the first element block and share the same properties (but not orientation)
         Call AskReal(ToughnessGrain,    'Grains: toughness', BatchUnit, IsBatch)
         Call AskReal(E,                 'Grains: Young''s modulus E', BatchUnit, IsBatch)
         Call AskReal(nu,                'Grains: Poisson Ration nu', BatchUnit, IsBatch)
         Call AskReal(ThermExpScalGrain, 'Grains: therm. exp.', BatchUnit, IsBatch)
      End If
      Do iBlock = 1, NumGrains
         Select Case(MeshTopology%Num_Dim)
         Case(2)
            MatProp2D(iBlock)%Toughness = ToughnessGrain
            Call GenHL_Iso2D_EnuPlaneStress(E, nu, MatProp2D(iBlock)%Hookes_Law)
            MatProp2D(iBlock)%Therm_Exp    = 0.0_Kr
            MatProp2D(iBlock)%Therm_Exp%XX = ThermExpScalGrain
            MatProp2D(iBlock)%Therm_Exp%YY = ThermExpScalGrain
         Case(3)
            MatProp3D(iBlock)%Toughness = ToughnessGrain
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, 'Polycrystal only supported in 2D at this point\n', iErr)
            MatProp3D(iBlock)%Therm_Exp    = 0.0_Kr
            MatProp3D(iBlock)%Therm_Exp%XX = ThermExpScalGrain
            MatProp3D(iBlock)%Therm_Exp%YY = ThermExpScalGrain
            MatProp3D(iBlock)%Therm_Exp%ZZ = ThermExpScalGrain
         End Select 
      End Do
      !!! Non grain EB if necessary
      Do iBlock = NumGrains+1, MeshTopology%Num_Elem_Blks_Global
         Write(IOBuffer, 100) iBlock
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      
         Write(IOBuffer, 300) iBlock, 'Toughness'
         Call AskReal(Toughness, IOBuffer, BatchUnit, IsBatch)
         Write(IOBuffer, 300) iBlock, 'Young Modulus'
         Call AskReal(E, IOBuffer, BatchUnit, IsBatch)
         Write(IOBuffer, 300) iBlock, 'Poisson Ratio'
         Call AskReal(nu, IOBuffer, BatchUnit, IsBatch)
         Write(IOBuffer, 300) iBlock, 'Therm Exp'
         Call AskReal(Therm_ExpScal, IOBuffer, BatchUnit, IsBatch)

         Select Case(MeshTopology%Num_Dim)
         Case(2)
            MatProp2D(iBlock)%Toughness = Toughness
            Call GenHL_Iso2D_EnuPlaneStress(E, nu, MatProp2D(iBlock)%Hookes_Law)
            MatProp2D(iBlock)%Therm_Exp    = 0.0_Kr
            MatProp2D(iBlock)%Therm_Exp%XX = Therm_ExpScal
            MatProp2D(iBlock)%Therm_Exp%YY = Therm_ExpScal
         Case(3)
            MatProp3D(iBlock)%Toughness = Toughness
            Call GenHL_Iso3D_Enu(E, nu, MatProp3D(iBlock)%Hookes_Law)
            MatProp3D(iBlock)%Therm_Exp    = 0.0_Kr
            MatProp3D(iBlock)%Therm_Exp%XX = Therm_ExpScal
            MatProp3D(iBlock)%Therm_Exp%YY = Therm_ExpScal
            MatProp3D(iBlock)%Therm_Exp%ZZ = Therm_ExpScal
         End Select 
      End Do
   End Select
   
   !!! Write EB Properties
   If (MEF90_MyRank == 0) Then
      Select Case(MeshTopology%Num_Dim)
      Case (2)
         Call MatProp_Write(MeshTopology, MatProp2D, Trim(prefix)//'.CST')
         DeAllocate(MatProp2D)
      Case(3)
         Call MatProp_Write(MeshTopology, MatProp3D, Trim(prefix)//'.CST')
         DeAllocate(MatProp3D)
      End Select
   End If
90 Format('Grain ', I4.4, ' angle ', F6.2, '\n')

!!!
!!! Temperature and Forces
!!!
   Write(IOBuffer, *) '\nTemperature and Forces\n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)      
   !!! Variable initialized on EB: F and Theta
   Allocate(F(MeshTopology%Num_Elem_Blks_Global))
   Allocate(Theta(MeshTopology%Num_Elem_Blks_Global))
   F%X = 0.0_Kr
   F%Y = 0.0_Kr
   F%Z = 0.0_Kr
   Theta = 0.0_Kr
   If (NumGrains>0) Then
      If (MyEXO%EBProperty(VarFrac_EBProp_HasBForce)%Value(1) /= 0) Then
         !!! Force
         Write(IOBuffer, 310) 'Fx'
         Call AskReal(FGrain%X, IOBuffer, BatchUnit, IsBatch)
         Write(IOBuffer, 310) 'Fy'
         Call AskReal(FGrain%Y, IOBuffer, BatchUnit, IsBatch)
         Write(IOBuffer, 310) 'Fz'
         Call AskReal(FGrain%Z, IOBuffer, BatchUnit, IsBatch)
      End If
      F(1:NumGrains)=F

      !!! Temperature
      Write(IOBuffer, 310) 'Theta'
      Call AskReal(ThetaGrain, IOBuffer, BatchUnit, IsBatch)
      Theta(1:NumGrains) = ThetaGrain
   End If
   
   Do i = Numgrains+1, MeshTopology%Num_Elem_Blks_Global
      Write(IOBuffer, 100) i
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      !!! Force
      If (MyEXO%EBProperty(VarFrac_EBProp_HasBForce)%Value(i) /= 0 ) Then
         Write(IOBuffer, 300) i, 'Fx'
         Call AskReal(F(i)%X, IOBuffer, BatchUnit, IsBatch)
         Write(IOBuffer, 300) i, 'Fy'
         Call AskReal(F(i)%Y, IOBuffer, BatchUnit, IsBatch)
         Write(IOBuffer, 300) i, 'Fz'
         Call AskReal(F(i)%Z, IOBuffer, BatchUnit, IsBatch)
      End If
      !!! Temperature
      Write(IOBuffer, 300) i, 'Theta'
      Call AskReal(Theta(i), IOBuffer, BatchUnit, IsBatch)
      !!!
      If (.NOT. IsBatch) Then
         Write(BatchUnit, *)
      End If
      
   End Do
   !!!
   !!! Compute the value of the force field at the vertices
   !!! Here is the place to request additional parameters if needed
   !!!
   Call MeshGetVertexSectionReal(MeshTopology%mesh, 'F', 3, FSec, iErr); CHKERRQ(iErr)
   Do_Step_F: Do iStep = 1, NumSteps
      Call SectionRealSet(FSec, 0.0_Kr, iErr); CHKERRQ(iErr)
      Do_Blk_F: Do iloc = 1, MeshTopology%Num_Elem_Blks
         i = MeshTopology%Elem_Blk(iLoc)%ID
         !!! Initialize the Section
         Num_DoF = MeshTopology%Elem_Blk(iloc)%Num_DoF
         Allocate(Felem(3*Num_DoF))
         Select Case(iCase)
         !!! Write special cases here
         Case default
            !!! Default is MIL
            If ( MyEXO%EBProperty(VarFrac_EBProp_HasBForce)%Value(i) /= 0 ) Then
               Do k = 0, Num_DoF-1
                  Felem(3*k+1) = T(iStep) * F(i)%X
                  Felem(3*k+2) = T(iStep) * F(i)%Y
                  Felem(3*k+3) = T(iStep) * F(i)%Z
               End Do
               Do j = 1, MeshTopology%Elem_Blk(iloc)%Num_Elems
                  Call SectionRealUpdateClosure(FSec, MeshTopology%Mesh, MeshTopology%Elem_Blk(iloc)%Elem_ID(j)-1, Felem, INSERT_VALUES, iErr); CHKERRQ(iErr)            
               End Do
            End If
         End Select      
         DeAllocate(Felem)
      End Do Do_Blk_F
      Call Write_EXO_Result_Vertex(MyEXO, MeshTopology, MyEXO%VertVariable(VarFrac_VertVar_ForceX)%Offset, iStep, FSec) 
   End Do Do_Step_F
   Call SectionRealDestroy(FSec, iErr); CHKERRQ(iErr)
   DeAllocate(F)
   !!!
   !!! Compute the value of the force field at the vertices
   !!! Here is the place to request additional parameters if needed
   !!!
   Call MeshGetVertexSectionReal(MeshTopology%mesh, 'Theta', 1, ThetaSec, iErr); CHKERRQ(iErr)
   Do_Step_Theta: Do iStep = 1, NumSteps
      Call SectionRealSet(ThetaSec, 0.0_Kr, iErr); CHKERRQ(iErr)
      Do_Blk_Theta: Do iloc = 1, MeshTopology%Num_Elem_Blks
         i = MeshTopology%Elem_Blk(iLoc)%ID
         !!! Initialize the Section
         Num_DoF = MeshTopology%Elem_Blk(iloc)%Num_DoF
         Allocate(Thetaelem(Num_DoF))
         Select Case(iCase)
         !!! Write special cases here
         Case(2,3,4)
            Call SectionRealSet(ThetaSec, 0.0_Kr, iErr); CHKERRQ(iErr)
         Case default
            !!! Default is MIL
            Do k = 1, Num_DoF
               Thetaelem(k) = T(iStep) * Theta(i)
            End Do
            Do j = 1, MeshTopology%Elem_Blk(iloc)%Num_Elems
               Call SectionRealUpdateClosure(ThetaSec, MeshTopology%Mesh, MeshTopology%Elem_Blk(iloc)%Elem_ID(j)-1, Thetaelem, INSERT_VALUES, iErr); CHKERRQ(iErr)            
            End Do
         End Select      
         DeAllocate(Thetaelem)
      End Do Do_Blk_Theta
      Call Write_EXO_Result_Vertex(MyEXO, MeshTopology, MyEXO%VertVariable(VarFrac_VertVar_Temperature)%Offset, iStep, ThetaSec) 
   End Do Do_Step_Theta
   Call SectionRealDestroy(ThetaSec, iErr); CHKERRQ(iErr)
   DeAllocate(Theta)


   !!!
   !!! U and V
   !!!
   Write(IOBuffer, *) '\nU and V\n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)      
   !!! Variable initialized on EB: F and Theta
   Allocate(U(MeshTopology%Num_Node_Sets_Global))
   Allocate(V(MeshTopology%Num_Node_Sets_Global))
   U%X = 0.0_Kr
   U%Y = 0.0_Kr
   U%Z = 0.0_Kr
   V   = 1.0_Kr
   Do i = 1, MeshTopology%Num_Node_Sets_Global
      Write(IOBuffer, 102) i
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      If (MyEXO%NSProperty(VarFrac_NSProp_BCUTypeX)%Value(i) /= 0 ) Then
         Write(IOBuffer, 302) i, 'Ux'
         Call AskReal(U(i)%X, IOBuffer, BatchUnit, IsBatch)
      End If
      If (MyEXO%NSProperty(VarFrac_NSProp_BCUTypeY)%Value(i) /= 0 ) Then
         Write(IOBuffer, 302) i, 'Uy'
         Call AskReal(U(i)%Y, IOBuffer, BatchUnit, IsBatch)
      End If
      If (MyEXO%NSProperty(VarFrac_NSProp_BCUTypeZ)%Value(i) /= 0 ) Then
         Write(IOBuffer, 302) i, 'Uz'
         Call AskReal(U(i)%Z, IOBuffer, BatchUnit, IsBatch)
      End If
      If (MyEXO%NSProperty(VarFrac_NSProp_BCVType)%Value(i) /= 0 ) Then
         Write(IOBuffer, 302) i, 'V'
         Call AskReal(V(i), IOBuffer, BatchUnit, IsBatch)
      End If
      If (.NOT. IsBatch) Then
         Write(BatchUnit, *)
      End If
   End Do
   Write(IOBuffer, *) '   \n\n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   !!!
   !!! Compute the value of the Displacement at the vertices
   !!! Here is the place to request additional parameters if needed
   !!!
   If ((iCase == 2) .OR. (iCase==3) .OR. (iCase==4)) Then
      Call AskReal(Eeff,  'E effective (for displacement field)',  BatchUnit, IsBatch)
      Call AskReal(nueff, 'nu effective (for displacement field)', BatchUnit, IsBatch)
      Kappa = (3.0-nu)/(1.0+nu)
      Mu = E / (1.0_Kr + nu) * .5_Kr
   End If
   Call MeshGetVertexSectionReal(MeshTopology%mesh, 'U', 3, USec, iErr); CHKERRQ(iErr)
   Allocate(Uelem(3))
   Do iStep = 1, NumSteps
      Call SectionRealSet(USec, 0.0_Kr, iErr); CHKERRQ(iErr)
      Do iloc = 1, MeshTopology%Num_Node_Sets         
         i = MeshTopology%Node_Set(iloc)%ID
         Select Case(iCase)
         !!! Write special cases here 
         Case(2,4)
            Do j = 1, MeshTopology%Node_Set(iloc)%Num_Nodes
               Call SectionRealRestrict(CoordSec, MeshTopology%Num_Elems + MeshTopology%Node_Set(iloc)%Node_ID(j)-1, Coordelem, iErr); CHKERRQ(iErr)
               R        = sqrt(CoordElem(1)**2 + CoordElem(2)**2)
               CTheta   = CoordElem(1) / R
               Ctheta2 = sqrt((1.0_Kr + CTheta) * .5_Kr)
               If (CoordElem(2) > 0.) Then
                  Stheta2 = sqrt((1.0_Kr - CTheta) * .5_Kr)
               Else
                  Stheta2 = -sqrt((1.0_Kr - CTheta) * .5_Kr)
               End If
               Uelem(1) = T(iStep) * sqrt(R / PETSC_PI * .5_Kr) / mu * .5_Kr * CTheta2 * (Kappa - CTheta) * U(i)%X
               Uelem(2) = T(iStep) * sqrt(R / PETSC_PI * .5_Kr) / mu * .5_Kr * STheta2 * (Kappa - CTheta) * U(i)%Y
               Uelem(3) = 0.
               Call SectionRealUpdate(USec, MeshTopology%Num_Elems + MeshTopology%Node_Set(iloc)%Node_ID(j)-1, Uelem, INSERT_VALUES, iErr); CHKERRQ(iErr)
               Call SectionRealRestore (CoordSec, MeshTopology%Num_Elems + MeshTopology%Node_Set(iloc)%Node_ID(j)-1, Coordelem, iErr); CHKERRQ(iErr)
            End Do   
         Case(3)
            Do j = 1, MeshTopology%Node_Set(iloc)%Num_Nodes
               Call SectionRealRestrict(CoordSec, MeshTopology%Num_Elems + MeshTopology%Node_Set(iloc)%Node_ID(j)-1, Coordelem, iErr); CHKERRQ(iErr)
               R        = sqrt((CoordElem(1)-T(iStep))**2 + CoordElem(2)**2)
               CTheta   = (CoordElem(1)-T(iStep))/R
               Ctheta2 = sqrt((1.0_Kr + CTheta) * .5_Kr)
               If (CoordElem(2) > 0.) Then
                  Stheta2 = sqrt((1.0_Kr - CTheta) * .5_Kr)
               Else
                  Stheta2 = -sqrt((1.0_Kr - CTheta) * .5_Kr)
               End If
               Uelem(1) = sqrt(R / PETSC_PI * .5_Kr) / mu * .5_Kr * CTheta2 * (Kappa - CTheta) * U(i)%X
               Uelem(2) = sqrt(R / PETSC_PI * .5_Kr) / mu * .5_Kr * STheta2 * (Kappa - CTheta) * U(i)%Y
               Uelem(3) = 0.
               Call SectionRealUpdate(USec, MeshTopology%Num_Elems + MeshTopology%Node_Set(iloc)%Node_ID(j)-1, Uelem, INSERT_VALUES, iErr); CHKERRQ(iErr)
               Call SectionRealRestore (CoordSec, MeshTopology%Num_Elems + MeshTopology%Node_Set(iloc)%Node_ID(j)-1, Coordelem, iErr); CHKERRQ(iErr)
            End Do   
         Case default
            !!! Default is MIL
            Uelem(1) = T(iStep) * U(i)%X
            Uelem(2) = T(iStep) * U(i)%Y
            Uelem(3) = T(iStep) * U(i)%Z
            Do j = 1, MeshTopology%Node_Set(iloc)%Num_Nodes
               Call SectionRealUpdate(USec, MeshTopology%Num_Elems + MeshTopology%Node_Set(iloc)%Node_ID(j)-1, Uelem, INSERT_VALUES, iErr); CHKERRQ(iErr)
            End Do
         End Select
      End Do
      Call Write_EXO_Result_Vertex(MyEXO, MeshTopology, MyEXO%VertVariable(VarFrac_VertVar_DisplacementX)%Offset, iStep, USec)
   End Do
   DeAllocate(Uelem)
   Call SectionRealDestroy(USec, iErr); CHKERRQ(iErr)
   DeAllocate(U)
   !!!
   !!! Compute the value of the fracture field at the vertices
   !!! Here is the place to request additional parameters if needed
   !!!
   Call MeshGetVertexSectionReal(MeshTopology%mesh, 'V', 1, VSec, iErr); CHKERRQ(iErr)
   Allocate(Velem(1))
   Do iStep = 1, NumSteps
      Call SectionRealSet(VSec, 1.0_Kr, iErr); CHKERRQ(iErr)
      Do iloc = 1, MeshTopology%Num_Node_Sets         
         i = MeshTopology%Node_Set(iloc)%ID
         Select Case(iCase)
         !!! Write special cases here 
         Case default
            Velem(1) = V(i)
            Do j = 1, MeshTopology%Node_Set(iloc)%Num_Nodes
               Call SectionRealUpdate(VSec, MeshTopology%Num_Elems + MeshTopology%Node_Set(iloc)%Node_ID(j)-1, Velem, INSERT_VALUES, iErr); CHKERRQ(iErr)
            End Do
         End Select
      End Do
      Call Write_EXO_Result_Vertex(MyEXO, MeshTopology, MyEXO%VertVariable(VarFrac_VertVar_Fracture)%Offset, iStep, VSec) 
   End Do
   DeAllocate(Velem)
   Call SectionRealDestroy(VSec, iErr); CHKERRQ(iErr)
   DeAllocate(V)

   Call PetscRandomDestroy(RandomCtx, iErr); CHKERRQ(iErr)
   Close(BatchUnit)
   DeAllocate(TestCase)
   DeAllocate(T)
   Call SectionRealDestroy(CoordSec, iErr); CHKERRQ(iErr)
   Call MeshTopologyDestroy(MeshTopology)
   Call MEF90_Finalize()

 100 Format('*** Element Block ', T24, I3, '\n')
! 101 Format('*** Side Set      ', T24, I3, '\n')
 102 Format('*** Node Set      ', T24, I3, '\n')
 300 Format('EB', I4.4, ': ', A)
 310 Format('Grains: ', A)
! 301 Format('SS', I4.4, ': ', A)
 302 Format('NS', I4.4, ': ', A)

End Program PrepVarFrac
