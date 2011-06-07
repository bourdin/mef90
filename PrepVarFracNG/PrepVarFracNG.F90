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
   Type(DM)                                     :: Tmp_Mesh
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
   PetscReal, Dimension(:), Pointer             :: V, Theta
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
   PetscReal                                    :: Eeff,nueff,Kappa,Mu
   PetscReal                                    :: CTheta,CTheta2,STheta,STheta2
   PetscReal, Dimension(:), Pointer             :: GlobVars
   PetscInt                                     :: vers
   PetscBool                                    :: IsBatch, HasBatchFile
   PetscInt                                     :: BatchUnit=99
   Character(len=MEF90_MXSTRLEN)                :: BatchFileName
   PetscBool                                    :: EraseBatch
   
   PetscRandom                                  :: RandomCtx
   PetscInt                                     :: Seed
   PetscLogDouble                               :: Time
   PetscBool                                    :: Has_Seed, Has_n
   PetscBool                                    :: saveElemVar, PlaneStrain
   
   PetscReal                                    :: R

   Call MEF90_Initialize()

   NumTestCase = 4
   Allocate(TestCase(NumTestCase))
   Do i = 1, NumTestCase
      TestCase(i)%Index = i
   End Do
   TestCase(1)%Description = "MIL"
   TestCase(2)%Description = "Mode-I scaling experiment"
   TestCase(3)%Description = "Mode-I surfing experiment"
   TestCase(4)%Description = "Single well, constant flux"

   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-verbose', verbose, HasPrefix, iErr)    
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p', prefix, HasPrefix, iErr); CHKERRQ(iErr)
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No mesh prefix given\n", iErr)
      Call MEF90_Finalize()
      STOP
   End If
   
   saveElemVar=PETSC_TRUE
   Call PetscOptionsGetBool(PETSC_NULL_CHARACTER, '-saveelemvar', saveElemVar, HasPrefix, iErr);CHKERRQ(iErr)
   PlaneStrain=PETSC_FALSE
   Call PetscOptionsGetBool(PETSC_NULL_CHARACTER, '-planestrain', PlaneStrain, HasPrefix, iErr);CHKERRQ(iErr)
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-i', BatchFileName, IsBatch, iErr); CHKERRQ(iErr)

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
   
   !!! Initialize random number generator
   Call PetscRandomCreate(PETSC_COMM_WORLD, RandomCtx, iErr); CHKERRQ(iErr)
   Call PetscRandomSetFromOptions(RandomCtx, iErr); CHKERRQ(iErr)
   Call PetscOptionsGetReal(PETSC_NULL_CHARACTER, '-seed', Seed, Has_Seed, iErr); CHKERRQ(iErr)
   If (.NOT. Has_Seed) Then
      Call PetscGetTime(Time, iErr); CHKERRQ(iErr)
      Seed =  Time * (Time - Int(Time))
   Else 
      Write(IOBuffer, *) 'Seeding random number generator with seed ', Seed, '\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If

   Call PetscRandomSetSeed(RandomCtx, Seed, iErr); CHKERRQ(iErr)
   Call PetscRandomSeed(RandomCtx, iErr); CHKERRQ(iErr)
   Call PetscRandomGetSeed(RandomCtx, Seed, iErr); CHKERRQ(iErr)

   
   !!!
   !!! Read mesh, partition
   !!!!
   Call Write_EXO_Case(prefix, '%0.4d', MEF90_NumProcs)
   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'
   
   Call EXO_Check_Numbering(EXO, iErr)
   If (iErr /= 0) Then
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, 'Unsupported numbering of the element blocks, side sets or node sets\n', iErr)
   End If

   !!! Reading and distributing sequential mesh
   If (MEF90_NumProcs == 1) Then
      Call DMMeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, MeshTopology%mesh, ierr); CHKERRQ(iErr)
   Else
      Call DMMeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
      Call DMMeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, MeshTopology%mesh, ierr); CHKERRQ(iErr)
      Call DMDestroy(Tmp_mesh, ierr); CHKERRQ(iErr)
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
   Call DMMeshGetSectionReal(MeshTopology%mesh, 'coordinates', CoordSec, iErr); CHKERRQ(ierr)
   
   !!!
   !!! EB, SS, NS Properties
   !!!
   Call VarFracEXOProperty_Init(MyEXO, MeshTopology)   
   If (verbose > 0) Then
      Write(IOBuffer, *) "Done with VarFracEXOProperty_Init\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   Write(IOBuffer, *) '\nElement Block and Node Set Properties\n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)      

   Call EXOEBProperty_AskWithBatch(MyEXO, MeshTopology, BatchUnit, IsBatch)
   If (verbose > 0) Then
      Write(IOBuffer, *) "Done with EXOEBProperty_AskWithBatch\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   Call EXOSSProperty_AskWithBatch(MyEXO, MeshTopology, BatchUnit, IsBatch)
   If (verbose > 0) Then
      Write(IOBuffer, *) "Done with EXOSSProperty_AskWithBatch\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   Call EXONSProperty_AskWithBatch(MyEXO, MeshTopology, BatchUnit, IsBatch)
   If (verbose > 0) Then
      Write(IOBuffer, *) "Done with EXONSProperty_AskWithBatch\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   Do i = 1, MeshTopology%num_elem_blks
      MeshTopology%elem_blk(i)%Elem_Type = MyEXO%EBProperty( VarFrac_EBProp_Elem_Type )%Value( MeshTopology%elem_blk(i)%ID )
      Call Init_Elem_Blk_Type(MeshTopology%Elem_Blk(i), MeshTopology%num_dim)
   End Do

   Call MeshTopologyWriteGlobal(MeshTopology, MyEXO, PETSC_COMM_WORLD)
   If (verbose > 0) Then
      Write(IOBuffer, *) "Done with MeshTopologyWriteGlobal\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If

   Call EXOProperty_Write(MyEXO)
   If (verbose > 0) Then
      Write(IOBuffer, '(A)') 'Done with EXOProperty_Write\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If

   !!!
   !!! Transient variables
   !!!
   Call VarFracEXOVariable_Init(MyEXO,saveElemVar)
   Call EXOVariable_Write(MyEXO)
   If (verbose > 0) Then
      Write(IOBuffer, '(A)') 'Done with EXOVariable_Write\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If

   !!!
   !!! Now we are done with the geometry, the variable definition and properties, so we can 
   !!! pre-compute the loads, displacement, time steps etc
   !!!
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
      Do i = 1, MeshTopology%Num_Elem_Blks_Global
         Write(IOBuffer, 100) i
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         
         Write(IOBuffer, 300) i, 'Toughness'
         Call AskReal(Toughness, IOBuffer, BatchUnit, IsBatch)
         Write(IOBuffer, 300) i, 'Young Modulus'
         Call AskReal(E, IOBuffer, BatchUnit, IsBatch)
         Write(IOBuffer, 300) i, 'Poisson Ratio'
         Call AskReal(nu, IOBuffer, BatchUnit, IsBatch)
         Write(IOBuffer, 300) i, 'Therm Exp'
         Call AskReal(Therm_ExpScal, IOBuffer, BatchUnit, IsBatch)

         Select Case(MeshTopology%Num_Dim)
         Case(2)
            MatProp2D(i)%Toughness = Toughness
                If (PlaneStrain) Then
                    Call GenHL_Iso2D_EnuPlaneStrain(E, nu, MatProp2D(i)%Hookes_Law)
                Else 
                    Call GenHL_Iso2D_EnuPlaneStress(E, nu, MatProp2D(i)%Hookes_Law)
                End If       
            MatProp2D(i)%Therm_Exp    = 0.0_Kr
            MatProp2D(i)%Therm_Exp%XX = Therm_ExpScal
            MatProp2D(i)%Therm_Exp%YY = Therm_ExpScal
         Case(3)
            MatProp3D(i)%Toughness = Toughness
            Call GenHL_Iso3D_Enu(E, nu, MatProp3D(i)%Hookes_Law)
            MatProp3D(i)%Therm_Exp    = 0.0_Kr
            MatProp3D(i)%Therm_Exp%XX = Therm_ExpScal
            MatProp3D(i)%Therm_Exp%YY = Therm_ExpScal
            MatProp3D(i)%Therm_Exp%ZZ = Therm_ExpScal
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
   Select Case(iCase)
   !!! Write special cases here
   Case Default
      Do i = 1, MeshTopology%Num_Elem_Blks_Global
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
         If (.NOT. IsBatch) Then
            Write(BatchUnit, *)
         End If
      End Do
   End Select
   
   !!!
   !!! Compute the value of the force field at the vertices
   !!! Here is the place to request additional parameters if needed
   !!!
   Call DMMeshGetVertexSectionReal(MeshTopology%mesh, 'F', 3, FSec, iErr); CHKERRQ(iErr)
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
   !!!
   Call DMMeshGetVertexSectionReal(MeshTopology%mesh, 'Theta', 1, ThetaSec, iErr); CHKERRQ(iErr)
   Do_Step_Theta: Do iStep = 1, NumSteps
      Call SectionRealSet(ThetaSec, 0.0_Kr, iErr); CHKERRQ(iErr)
      Do_Blk_Theta: Do iloc = 1, MeshTopology%Num_Elem_Blks
         i = MeshTopology%Elem_Blk(iLoc)%ID
         !!! Initialize the Section
         Num_DoF = MeshTopology%Elem_Blk(iloc)%Num_DoF
         Allocate(Thetaelem(Num_DoF))
         Select Case(iCase)
         !!! Write special cases here
         Case(5) !!! Single well, cst flux
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,'Not implemented yet\n',iErr)
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
   Select Case(iCase)
   !!! Write special cases here
   Case(2,3) !!! Scaling, surfing
      Call AskReal(Eeff,  'E effective (for displacement field)',  BatchUnit, IsBatch)
      Call AskReal(nueff, 'nu effective (for displacement field)', BatchUnit, IsBatch)
      Kappa = (3.0-nu)/(1.0+nu)
      Mu = E / (1.0_Kr + nu) * .5_Kr
   Case Default
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
   End Select
   Write(IOBuffer, *) '   \n\n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   !!!
   !!! Compute the value of the Displacement at the vertices
   !!!
   Call DMMeshGetVertexSectionReal(MeshTopology%mesh, 'U', 3, USec, iErr); CHKERRQ(iErr)
   Allocate(Uelem(3))
   Do iStep = 1, NumSteps
      Call SectionRealSet(USec, 0.0_Kr, iErr); CHKERRQ(iErr)
      Do iloc = 1, MeshTopology%Num_Node_Sets         
         i = MeshTopology%Node_Set(iloc)%ID
         Select Case(iCase)
         !!! Write special cases here 
         Case(2)
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
         Case Default
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
   Call DMMeshGetVertexSectionReal(MeshTopology%mesh, 'V', 1, VSec, iErr); CHKERRQ(iErr)
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
! 310 Format('Grains: ', A)
! 301 Format('SS', I4.4, ': ', A)
 302 Format('NS', I4.4, ': ', A)

End Program PrepVarFrac
