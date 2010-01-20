Program PrepVarFrac

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscmeshdef.h"

   Use m_MEF90
   Use m_VarFrac_Struct

   Use petsc
   Use petscvec
   Use petscmat
   Use petscksp
   Use petscmesh

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
   PetscTruth                                   :: HasPrefix
   PetscInt                                     :: verbose = 0
   PetscInt                                     :: iErr, iloc, i, j, k, iCase, iStep

   PetscInt                                     :: NumTestCase
   Type(TestCase_Type), Dimension(:) , Pointer  :: TestCase
   PetscInt, Parameter                          :: QuadOrder=2
   Type(SectionReal)                            :: USec, FSec, VSec, ThetaSec, CoordSec
   Type(Vect3D), Dimension(:), Pointer          :: U, F
   PetscReal, Dimension(:), Pointer             :: V
   PetscReal, Dimension(:), Pointer             :: Theta
   PetscReal, Dimension(:), Pointer             :: Uelem, Felem, Velem, Thetaelem, Coordelem
   PetscReal                                    :: Tmin, Tmax
   PetscReal, Dimension(:), Pointer             :: T
   PetscReal                                    :: Beta, eta, tau, Y
   PetscReal                                    :: FixedPoint
   PetscInt                                     :: NumSteps
   PetscInt                                     :: Num_DoF
   PetscReal                                    :: RealBuffer
   Type(MatProp2D_Type), Dimension(:), Pointer  :: MatProp2D
   Type(MatProp3D_Type), Dimension(:), Pointer  :: MatProp3D
   PetscReal                                    :: E, nu, Toughness, Therm_ExpScal
   PetscReal, Dimension(:), Pointer             :: GlobVars
   PetscReal                                    :: rDummy
   Character                                    :: cDummy
   PetscInt                                     :: vers
   PetscTruth                                   :: IsBatch, HasBatchFile
   PetscInt                                     :: BatchUnit=99
   Character(len=MEF90_MXSTRLEN)                :: BatchFileName
   PetscTruth                                   :: EraseBatch

   Call MEF90_Initialize()
   
   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-verbose', verbose, HasPrefix, iErr)    
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p', prefix, HasPrefix, iErr); CHKERRQ(iErr)
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No mesh prefix given\n", iErr)
      Call MEF90_Finalize()
      STOP
   End If
   Call PetscOptionsGetTruth(PETSC_NULL_CHARACTER, '-force', EraseBatch, HasPrefix, iErr)    
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-i', BatchFileName, IsBatch, iErr); CHKERRQ(iErr)
   If (IsBatch) Then
      Write(IOBuffer, *) "Processing batch file ", Trim(BatchFileName), "\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Open(Unit=BatchUnit, File=BatchFileName, Status='Old')
      Rewind(BatchUnit)
   Else
      BatchFileName = Trim(prefix)//'.args'
      Inquire(File=BatchFileName, EXIST=HasBatchFile)
      Write(IOBuffer, *) "Running interactively and generating batch file ", trim(BatchFileName), "\n"
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      If (HasBatchFile .AND. (.NOT. EraseBatch)) Then
         Write(IOBuffer, *) "Batch file ", trim(BatchFileName), " already exists. Erase it and restart\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Call MEF90_Finalize()
         STOP
      End If
      Open(Unit=BatchUnit, File=BatchFileName, Status='Unknown')
      Rewind(BatchUnit)
   End If
   
   NumTestCase = 6
   Allocate(TestCase(NumTestCase))
   Do i = 1, NumTestCase
      TestCase(i)%Index = i
   End Do
   TestCase(1)%Description = "MIL 2D/3D elasticity                       "
   TestCase(2)%Description = "Geothermal thermal cracks: proof of concept"
   TestCase(3)%Description = "Cooling: infinite domain, thermal conduction only"
   TestCase(4)%Description = "Cooling: infinite domain, convection and conduction (Newtonian cooling)"
   TestCase(5)%Description = "MIL affine loading"
   TestCase(6)%Description = "Loads given by a polar angle (2D)"
   

   Call Write_EXO_Case(prefix, '%0.4d', MEF90_NumProcs)
   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'
   
   Call EXO_Check_Numbering(EXO, iErr)
   If (iErr /= 0) Then
      SETERRQ(PETSC_ERR_SUP, 'Unsupported numbering of the element blocks, side sets or node sets\n', iErr)
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
   Call EXOProperty_AskWithBatch(MyEXO, MeshTopology, BatchUnit, IsBatch)
   
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
      SETERRQ(PETSC_ERR_SUP, 'Only 2 and 3 dimensional elements are supported', iErr)
   End Select
   If (verbose > 0) Then
      Write(IOBuffer, '(A)') 'Done with ElementInit\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
!!! Initialize Sections   
   Call MeshGetVertexSectionReal(MeshTopology%mesh, 'U', 3, USec, iErr); CHKERRQ(iErr)
   Call MeshGetVertexSectionReal(MeshTopology%mesh, 'F', 3, FSec, iErr); CHKERRQ(iErr)
   Call MeshGetVertexSectionReal(MeshTopology%mesh, 'V', 1, VSec, iErr); CHKERRQ(iErr)
   Call MeshGetVertexSectionReal(MeshTopology%mesh, 'Theta', 1, ThetaSec, iErr); CHKERRQ(iErr)
   
   If (verbose > 0) Then
      Write(IOBuffer, '(A)') 'Done with Initializing Sections\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   
   Do i = 1, NumTestCase
      Write(IOBuffer, "('[',I2.2,'] ',A)"), TestCase(i)%Index, Trim(TestCase(i)%Description)//'\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End Do
   
   Call AskInt(iCase, 'Test Case', BatchUnit, IsBatch)

   Select Case(iCase)
   Case (1,2,3,4,5,6)! MIL, geothermal PoC

      !!! Time Steps
      Write(IOBuffer, *) '\nGlobal Variables\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)      
      Call AskReal(TMin, 'TMin', BatchUnit, IsBatch)
      Call AskReal(TMax, 'TMax', BatchUnit, IsBatch)
      Call AskInt(NumSteps, 'NumSteps', BatchUnit, IsBatch)
      
      If (iCase == 4) Then
         Call AskReal(Beta, 'Beta', BatchUnit, IsBatch)
      End If
      
      If (iCase == 5) Then
         Call AskReal(FixedPoint, 'FixedPt', BatchUnit, IsBatch)
      End If

      If (iCase == 4) Then
         Call AskReal(Beta, 'polar angle', BatchUnit, IsBatch)
      End If

      Allocate(GlobVars(VarFrac_Num_GlobVar))
      GlobVars = 0.0_Kr
      Allocate(T(NumSteps))
      Do i = 1, NumSteps-1
         If ( (iCase == 3) .OR. (iCase == 4) ) Then
            T(i) = Tmin + (Real(i-1) * (sqrt(Tmax-TMin))/Real(NumSteps-1))**2
         Else
            T(i) = Tmin + Real(i-1) * (Tmax-TMin)/Real(NumSteps-1)
         End If
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

     !!! Elem Blocks BC and Variables
   
      Select Case(MeshTopology%Num_Dim)
      Case(2)
         Allocate(MatProp2D(MeshTopology%Num_Elem_Blks_Global))
      Case(3)
         Allocate(MatProp3D(MeshTopology%Num_Elem_Blks_Global))
      End Select
      
      !!! Material Properties
      Write(IOBuffer, *) '\nMaterial Properties\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)      
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
         Call AskReal(THerm_ExpScal, IOBuffer, BatchUnit, IsBatch)

         Select Case(MeshTopology%Num_Dim)
         Case(2)
            MatProp2D(i)%Toughness = Toughness
            Call GenHL_Iso2D_EnuPlaneStress(E, nu, MatProp2D(i)%Hookes_Law)
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
      If (MEF90_MyRank == 0) Then
         Select Case(MeshTopology%Num_Dim)
         Case (2)
            Call MatProp_Write(MeshTopology, MatProp2D, Trim(prefix)//'.CST')
         Case(3)
            Call MatProp_Write(MeshTopology, MatProp3D, Trim(prefix)//'.CST')
         End Select
      End If



      Write(IOBuffer, *) '\nFields and Loads\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)      
      !!! Variable initialized on EB: F and Theta
      Allocate(F(MeshTopology%Num_Elem_Blks_Global))
      Allocate(Theta(MeshTopology%Num_Elem_Blks_Global))
      F%X = 0.0_Kr
      F%Y = 0.0_Kr
      F%Z = 0.0_Kr
      Theta = 0.0_Kr
      
      Do i = 1, MeshTopology%Num_Elem_Blks_Global
         Write(IOBuffer, 100) i
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

         !!! Force
         If (MyEXO%EBProperty(VarFrac_EBProp_HasBForce)%Value(i) /= 0 ) Then
            Write(IOBuffer, 300) i, 'Fx'
            Call AskReal(F(i)%X, IOBuffer, BatchUnit, IsBatch)

            Write(IOBuffer, 300) i, 'Fy'
            Call AskReal(F(i)%Y, IOBuffer, BatchUnit, IsBatch)

            Write(IOBuffer, 300) i, 'FZ'
            Call AskReal(F(i)%Z, IOBuffer, BatchUnit, IsBatch)
         End If
         
         !!! Temperature
         Write(IOBuffer, 300) i, 'Theta'
         Call AskReal(Theta(i), IOBuffer, BatchUnit, IsBatch)
      End Do   

                  
      Do iStep = 1, NumSteps
         Call SectionRealSet(FSec, 0.0_Kr, iErr); CHKERRQ(iErr)
         Call SectionRealSet(ThetaSec, 0.0_Kr, iErr); CHKERRQ(iErr)
         Do iloc = 1, MeshTopology%Num_Elem_Blks
            i = MeshTopology%Elem_Blk(iLoc)%ID
            !!! Initialize the Section
            Num_DoF = MeshTopology%Elem_Blk(iloc)%Num_DoF
            Allocate(Felem(3*Num_DoF))
            Allocate(Thetaelem(Num_DoF))
            Allocate(Coordelem(Num_DoF * MeshTopology%Num_Dim))
            
            !!! Update F
            If ( MyEXO%EBProperty(VarFrac_EBProp_HasBForce)%Value(i) /= 0 ) Then
               Select Case (iCase)
               Case(1,2,3,4,5)
                  Do k = 0, Num_DoF-1
                     Felem(3*k+1) = T(iStep) * F(i)%X
                     Felem(3*k+2) = T(iStep) * F(i)%Y
                     Felem(3*k+3) = T(iStep) * F(i)%Z
                  End Do
                  Do j = 1, MeshTopology%Elem_Blk(iloc)%Num_Elems
                     Call SectionRealUpdateClosure(FSec, MeshTopology%Mesh, MeshTopology%Elem_Blk(iloc)%Elem_ID(j)-1, Felem, INSERT_VALUES, iErr); CHKERRQ(iErr)            
                  End Do
               Case(6)
                  Do k = 0, Num_DoF-1
                     Felem(3*k+1) = cos(T(iStep)*PETSC_PI/180.0_Kr) * F(i)%X
                     Felem(3*k+2) = sin(T(iStep)*PETSC_PI/180.0_Kr) * F(i)%Y
                     Felem(3*k+3) = 0.
                  End Do
                  Do j = 1, MeshTopology%Elem_Blk(iloc)%Num_Elems
                     Call SectionRealUpdateClosure(FSec, MeshTopology%Mesh, MeshTopology%Elem_Blk(iloc)%Elem_ID(j)-1, Felem, INSERT_VALUES, iErr); CHKERRQ(iErr)            
                  End Do
               End Select
            End If
   
            !!! Update Theta
            Select Case (iCase)
            Case(1)
               Thetaelem = T(iSTep) * Theta(i)
               Do j = 1, MeshTopology%Elem_Blk(iloc)%Num_Elems
                  Call SectionRealUpdateClosure(ThetaSec, MeshTopology%Mesh, MeshTopology%Elem_Blk(iloc)%Elem_ID(j)-1, Thetaelem, INSERT_VALUES, iErr); CHKERRQ(iErr) 
               End Do
            Case(2)
               Do j = 1, MeshTopology%Elem_Blk(iloc)%Num_Elems
                  Call SectionRealRestrictClosure(CoordSec, MeshTopology%mesh, MeshTopology%Elem_Blk(iloc)%Elem_ID(j)-1, Num_DoF * MeshTopology%Num_Dim, CoordElem, iErr); CHKERRQ(iErr)
                  Do k = 1, Num_DoF
                     ThetaElem(k) = erf( CoordElem((k-1) * MeshTopology%Num_Dim + 2)**2 / T(iStep) )
                  End Do
                  ThetaElem = Theta(i) * (1.0-ThetaElem)
                  Call SectionRealUpdateClosure(ThetaSec, MeshTopology%Mesh, MeshTopology%Elem_Blk(iloc)%Elem_ID(j)-1, Thetaelem, INSERT_VALUES, iErr); CHKERRQ(iErr) 
               End Do
            Case(3)
               Do j = 1, MeshTopology%Elem_Blk(iloc)%Num_Elems
                  Call SectionRealRestrictClosure(CoordSec, MeshTopology%mesh, MeshTopology%Elem_Blk(iloc)%Elem_ID(j)-1, Num_DoF * MeshTopology%Num_Dim, CoordElem, iErr); CHKERRQ(iErr)
                  Do k = 1, Num_DoF
                     ThetaElem(k) = erf( CoordElem((k-1) * MeshTopology%Num_Dim + 2) / Sqrt(T(iStep)) * 0.5_Kr )
                  End Do
                  ThetaElem = Theta(i) * (1.0-ThetaElem)
                  Call SectionRealUpdateClosure(ThetaSec, MeshTopology%Mesh, MeshTopology%Elem_Blk(iloc)%Elem_ID(j)-1, Thetaelem, INSERT_VALUES, iErr); CHKERRQ(iErr) 
               End Do
            Case(4)
               Do j = 1, MeshTopology%Elem_Blk(iloc)%Num_Elems
                  Call SectionRealRestrictClosure(CoordSec, MeshTopology%mesh, MeshTopology%Elem_Blk(iloc)%Elem_ID(j)-1, Num_DoF * MeshTopology%Num_Dim, CoordElem, iErr); CHKERRQ(iErr)
                  Do k = 1, Num_DoF
                     eta = CoordElem((k-1) * MeshTopology%Num_Dim + 2)
                     tau  = sqrt(T(iStep))
                     ThetaElem(k) = erfc(eta / tau * 0.5_Kr ) - exp(beta * eta + beta**2 * tau**2) * erfc(eta / tau * 0.5_Kr + beta * tau)
                  End Do
                  ThetaElem = Theta(i) * ThetaElem
                  Call SectionRealUpdateClosure(ThetaSec, MeshTopology%Mesh, MeshTopology%Elem_Blk(iloc)%Elem_ID(j)-1, Thetaelem, INSERT_VALUES, iErr); CHKERRQ(iErr) 
               End Do
            End Select
            DeAllocate(Felem)
            DeAllocate(Thetaelem)         
            DeAllocate(CoordElem)
         End Do
         
         Call Write_EXO_Result_Vertex(MyEXO, MeshTopology, MyEXO%VertVariable(VarFrac_VertVar_ForceX)%Offset, iStep, FSec) 
         Call Write_EXO_Result_Vertex(MyEXO, MeshTopology, MyEXO%VertVariable(VarFrac_VertVar_Temperature)%Offset, iStep, ThetaSec) 
      End Do
      DeAllocate(F)
      DeAllocate(Theta)
      

     !!! Variables Initialized at NS
      Allocate(U(MeshTopology%Num_Node_Sets_Global))
      Allocate(V(MeshTopology%Num_Node_Sets_Global))
      U%X = 0.0_Kr
      U%Y = 0.0_Kr
      U%Z = 0.0_Kr
      V   = 1.0_Kr
      Do i = 1, MeshTopology%Num_Node_Sets_Global
         Write(IOBuffer, 102) i
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

         !!! Displacement
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
      End Do
      
      Write(IOBuffer, *) '   \n\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      !!! Initialize the Section
      Allocate(Uelem(3))
      Allocate(Velem(1))
      
      Do iStep = 1, NumSteps
         Call SectionRealSet(USec, 0.0_Kr, iErr); CHKERRQ(iErr)
         Call SectionRealSet(VSec, 1.0_Kr, iErr); CHKERRQ(iErr)
         Do iloc = 1, MeshTopology%Num_Node_Sets         
            i = MeshTopology%Node_Set(iloc)%ID
            
            Select Case (iCase)
            Case(1,2,3,4)
               Uelem(1) = T(iStep) * U(i)%X
               Uelem(2) = T(iStep) * U(i)%Y
               Uelem(3) = T(iStep) * U(i)%Z
               Velem    = V(i)
               Do j = 1, MeshTopology%Node_Set(iloc)%Num_Nodes
                  Call SectionRealUpdateClosure(USec, MeshTopology%Mesh, MeshTopology%Num_Elems + MeshTopology%Node_Set(iloc)%Node_ID(j)-1, Uelem, INSERT_VALUES, iErr); CHKERRQ(iErr)            
                  Call SectionRealUpdateClosure(VSec, MeshTopology%Mesh, MeshTopology%Num_Elems + MeshTopology%Node_Set(iloc)%Node_ID(j)-1, Velem, INSERT_VALUES, iErr); CHKERRQ(iErr)            
               End Do
            Case(5)
               Velem    = V(i)
               Allocate(Coordelem(MeshTopology%Num_Dim))
               Do j = 1, MeshTopology%Node_Set(iloc)%Num_Nodes
                  Call SectionRealRestrictClosure(CoordSec, MeshTopology%mesh, MeshTopology%Num_Elems + MeshTopology%Node_Set(iloc)%Node_ID(j)-1, MeshTopology%Num_Dim, CoordElem, iErr); CHKERRQ(iErr)
                  Uelem(1) = T(iStep) * U(i)%X
                  Uelem(2) = T(iStep) * U(i)%Y * (FixedPoint-CoordElem(1))
                  Uelem(3) = T(iStep) * U(i)%Z
                  Call SectionRealUpdateClosure(USec, MeshTopology%Mesh, MeshTopology%Num_Elems + MeshTopology%Node_Set(iloc)%Node_ID(j)-1, Uelem, INSERT_VALUES, iErr); CHKERRQ(iErr)            
                  Call SectionRealUpdateClosure(VSec, MeshTopology%Mesh, MeshTopology%Num_Elems + MeshTopology%Node_Set(iloc)%Node_ID(j)-1, Velem, INSERT_VALUES, iErr); CHKERRQ(iErr)            
               End Do
               DeAllocate(Coordelem)
            Case(6)
               Uelem(1) = cos(T(iStep)*PETSC_PI/180_Kr) * U(i)%X
               Uelem(2) = sin(T(iStep)*PETSC_PI/180_Kr) * U(i)%Y
               Uelem(3) = 0
               Velem    = V(i)
               Do j = 1, MeshTopology%Node_Set(iloc)%Num_Nodes
                  Call SectionRealUpdateClosure(USec, MeshTopology%Mesh, MeshTopology%Num_Elems + MeshTopology%Node_Set(iloc)%Node_ID(j)-1, Uelem, INSERT_VALUES, iErr); CHKERRQ(iErr)            
                  Call SectionRealUpdateClosure(VSec, MeshTopology%Mesh, MeshTopology%Num_Elems + MeshTopology%Node_Set(iloc)%Node_ID(j)-1, Velem, INSERT_VALUES, iErr); CHKERRQ(iErr)            
               End Do
            End Select
         End Do
         Call Write_EXO_Result_Vertex(MyEXO, MeshTopology, MyEXO%VertVariable(VarFrac_VertVar_Fracture)%Offset, iStep, VSec) 
         Call Write_EXO_Result_Vertex(MyEXO, MeshTopology, MyEXO%VertVariable(VarFrac_VertVar_DisplacementX)%Offset, iStep, USec) 
      End Do
      
      DeAllocate(Uelem)
      DeAllocate(Velem)
      DeAllocate(U)
      DeAllocate(V)
      
   Case Default
      SETERRQ(PETSC_ERR_SUP, 'Unknown test case\n', iErr)      
   End Select

   Close(BatchUnit)
   Call MEF90_Finalize()

 100 Format('*** Element Block ', T24, I3, '\n')
 101 Format('*** Side Set      ', T24, I3, '\n')
 102 Format('*** Node Set      ', T24, I3, '\n')
 !200 Format(A,t60, ': ')
 300 Format('EB', I4.4, ': ', A)
 301 Format('SS', I4.4, ': ', A)
 302 Format('NS', I4.4, ': ', A)

Contains
   Subroutine AskInt(val, msg, ArgUnit, IsBatch)
      PetscInt                                  :: Val
      Character(len=*)                          :: msg 
      PetscInt                                  :: argunit
      PetscTruth                                :: IsBatch

      Character(len=MEF90_MXSTRLEN)             :: prefix, IOBuffer      
      If (IsBatch) Then
         If (MEF90_MyRank == 0) Then
            Read(ArgUnit,*) Val
         End If
         Call MPI_BCast(Val, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, iErr)
      Else
         Write(IOBuffer, "(A, t60,':  ')") Trim(msg)
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         If (MEF90_MyRank == 0) Then
            Read(*,*) Val
            Write(ArgUnit, "(I4, t60, A)") val, Trim(msg)
         End If
         Call MPI_BCast(Val, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, iErr)
      End If
   End Subroutine AskInt   
   
   Subroutine AskReal(val, msg, ArgUnit, IsBatch)
      PetscReal                                 :: Val
      Character(len=*)                          :: msg 
      PetscInt                                  :: argunit
      PetscTruth                                :: IsBatch

      Character(len=MEF90_MXSTRLEN)             :: prefix, IOBuffer      
      If (IsBatch) Then
         If (MEF90_MyRank == 0) Then
            Read(ArgUnit,*) Val
         End If
         Call MPI_BCast(Val, 1, MPIU_SCALAR, 0, PETSC_COMM_WORLD, iErr)
      Else
         Write(IOBuffer, "(A, t60,':  ')") Trim(msg)
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         If (MEF90_MyRank == 0) Then
            Read(*,*) Val
            Write(ArgUnit, "(ES12.5, t60, A)") val, Trim(msg)
         End If
         Call MPI_BCast(Val, 1, MPIU_SCALAR, 0, PETSC_COMM_WORLD, iErr)
      End If
   End Subroutine AskReal

   Subroutine EXOProperty_AskWithBatch(dEXO, dMeshTopology, BatchUnit, IsBatch)
      Type(EXO_Type)                                 :: dEXO
      Type(MeshTopology_Type)                        :: dMeshTopology
      PetscInt                                       :: BatchUnit
      PetscTruth                                     :: IsBatch
      
      PetscInt                                       :: iErr
      PetscInt                                       :: i, j, IntBuffer

      PetscInt                                       :: NumEB, NumSS, NumNS
      PetscInt                                       :: EXO_MyRank
      Character(len=MEF90_MXSTRLEN)                  :: IOBuffer

      Do i = 1, dMeshTopology%Num_Elem_Blks_Global
         Write(IOBuffer, 100) i
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Do j = 1, dEXO%Num_EBProperties
            Write(IOBuffer, 200) i, Trim(dEXO%EBProperty(j)%Name)
            Call AskInt(dEXO%EBProperty(j)%Value(i), IOBuffer, BatchUnit, IsBatch)
         End Do
      End Do
      
      Do i = 1, dMeshTopology%Num_Side_Sets_Global
         Write(IOBuffer, 101) i
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Do j = 1, dEXO%Num_SSProperties
            Write(IOBuffer, 201) i, Trim(dEXO%SSProperty(j)%Name)
            Call AskInt(dEXO%SSProperty(j)%Value(i), IOBuffer, BatchUnit, IsBatch)
         End Do
      End Do

      Do i = 1, dMeshTopology%Num_Node_Sets_Global
         Write(IOBuffer, 102) i
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Do j = 1, dEXO%Num_NSProperties
            Write(IOBuffer, 202) i, Trim(dEXO%NSProperty(j)%Name)
            Call AskInt(dEXO%NSProperty(j)%Value(i), IOBuffer, BatchUnit, IsBatch)
         End Do
      End Do
 100 Format('*** Element Block ', T24, I3, '\n')
 101 Format('*** Side Set      ', T24, I3, '\n')
 102 Format('*** Node Set      ', T24, I3, '\n')
 110 Format(T24, A, T60, ': ')
 200 Format('EB', I4.4, ': ', A)
 201 Format('SS', I4.4, ': ', A)
 202 Format('NS', I4.4, ': ', A)
   End Subroutine EXOProperty_AskWithBatch
      

End Program PrepVarFrac
