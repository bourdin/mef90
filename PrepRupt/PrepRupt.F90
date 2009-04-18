Program PrepRupt

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscmeshdef.h"

   Use m_MEF90
   Use m_RuptStruct

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
   PetscTruth                                   :: HasPrefix, verbose
   PetscInt                                     :: iErr, iloc, i, j, k, iCase, step

   PetscInt, Parameter                          :: NumTestCase=2
   Type(TestCase_Type), Dimension(NumTestCase)  :: TestCase
   PetscInt, Parameter                          :: QuadOrder=2
   Type(SectionReal)                            :: USec, FSec, ThetaSec
   Type(Vect3D), Dimension(:), Pointer          :: U, F
   PetscReal, Dimension(:), Pointer             :: Theta
   PetscReal, Dimension(:), Pointer             :: Uelem, Felem, Thetaelem
   PetscReal                                    :: Tmin, Tmax
   PetscReal, Dimension(:), Pointer             :: T
   PetscInt                                     :: NumSteps
   PetscInt                                     :: Num_DoF
   PetscReal                                    :: RealBuffer
   Type(MatProp2D_Type), Dimension(:), Pointer  :: MatProp2D
   Type(MatProp3D_Type), Dimension(:), Pointer  :: MatProp3D
   PetscReal                                    :: E, nu, Toughness, Therm_ExpScal

   Call MEF90_Initialize()
   
   Do i = 1, NumTestCase
      TestCase(i)%Index = i
   End Do
   TestCase(1)%Description = "MIL 2D/3D elasticity    "
   TestCase(2)%Description = "MIL antiplane elasticity"
   
   Call PetscOptionsHasName(PETSC_NULL_CHARACTER, '-verbose', verbose, iErr)    
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p',       prefix, HasPrefix, iErr); CHKERRQ(iErr)
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No input file prefix given\n", iErr)
      Call MEF90_Finalize()
      STOP
   End If
   
   Call Write_EXO_Case(prefix, '%0.4d', MEF90_NumProcs)
   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'
   
   Call EXO_Check_Numbering(EXO, iErr)
   If (iErr /= 0) Then
      SETERRQ(PETSC_ERR_SUP, 'Unsupported numbering of the element blocks, side sets or node sets\n'c, iErr)
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
   If (verbose) Then
      Write(IOBuffer, *) "Done reading and partitioning the mesh\n"c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If

   MyEXO%comm = PETSC_COMM_SELF
   MyEXO%exoid = EXO%exoid
   Write(MyEXO%filename, 99) trim(prefix), MEF90_MyRank
 99  Format(A, '-', I4.4, '.gen')

   
   Call RuptEXOProperty_Init(MyEXO, MeshTopology)   
   If (verbose) Then
      Write(IOBuffer, *) "Done with RuptEXOProperty_Init\n"c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   Call EXOProperty_Ask(MyEXO, MeshTopology)
   
   Do i = 1, MeshTopology%num_elem_blks
      MeshTopology%elem_blk(i)%Elem_Type = MyEXO%EBProperty( Rupt_EBProp_Elem_Type )%Value( MeshTopology%elem_blk(i)%ID )
      Call Init_Elem_Blk_Type(MeshTopology%Elem_Blk(i), MeshTopology%num_dim)
   End Do



   Call Write_MeshTopologyGlobal(MeshTopology, MyEXO, PETSC_COMM_WORLD)
   If (verbose) Then
      Write(IOBuffer, *) "Done with Write_MeshTopologyGlobal\n"c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If

   Call EXOProperty_Write(MyEXO)
   If (verbose) Then
      Write(IOBuffer, '(A)') 'Done with EXOProperty_Write\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If


   Call RuptEXOVariable_Init(MyEXO)
   Call EXOVariable_Write(MyEXO)
   If (verbose) Then
      Write(IOBuffer, '(A)') 'Done with EXOVariable_Write\n'c
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
   If (verbose) Then
      Write(IOBuffer, '(A)') 'Done with ElementInit\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
!!! Initialize Sections   
   Call MeshGetVertexSectionReal(MeshTopology%mesh, 3, USec, iErr); CHKERRQ(iErr)
   Call MeshGetVertexSectionReal(MeshTopology%mesh, 3, FSec, iErr); CHKERRQ(iErr)
   Call MeshGetVertexSectionReal(MeshTopology%mesh, 1, ThetaSec, iErr); CHKERRQ(iErr)
   If (verbose) Then
      Write(IOBuffer, '(A)') 'Done with Initializing Sections\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   
   Do i = 1, NumTestCase
      Write(IOBuffer, "('[',I2.2,'] ',A)"), TestCase(i)%Index, Trim(TestCase(i)%Description)//'\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End Do
   
   Write(IOBuffer, 200) 'Test Case'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   If (MEF90_MyRank == 0) Then
      Read(*,*) iCase
   End If
   Call MPI_BCast(iCase, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, iErr)

   Select Case(iCase)
   Case (1)! MIL
      !!! Time Steps
      Write(IOBuffer, 200) 'TMin'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      If (MEF90_MyRank == 0) Then
         Read(*,*) TMin
      End If
      Call MPI_BCast(TMin, 1, MPIU_SCALAR, 0, EXO%Comm, iErr)
      Write(IOBuffer, 200) 'TMax'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      If (MEF90_MyRank == 0) Then
         Read(*,*) TMax
      End If
      Call MPI_BCast(TMax, 1, MPIU_SCALAR, 0, EXO%Comm, iErr)
      Write(IOBuffer, 200) 'NumSteps'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      If (MEF90_MyRank == 0) Then
         Read(*,*) NumSteps
      End If
      Call MPI_BCast(NumSteps, 1, MPI_INTEGER, 0, EXO%Comm, iErr)
      
      Allocate(T(NumSteps))
      Do k = 1, NumSteps-1
         T(i) = Tmin + Real(i-1) * (Tmax-TMin)/Real(NumSteps-1)
      End Do
      T(NumSteps) = TMax
         
     !!! Elem Blocks BC and Variables
      Allocate(U(MeshTopology%Num_Elem_Blks_Global))
      Allocate(F(MeshTopology%Num_Elem_Blks_Global))
      Allocate(Theta(MeshTopology%Num_Elem_Blks_Global))
      Select Case (MeshTopology%Num_Dim)
      Case(2)
         Allocate(MatProp2D(MeshTopology%Num_Elem_Blks_Global))
      Case(3)
         Allocate(MatProp3D(MeshTopology%Num_Elem_Blks_Global))
      End Select
      
      U%X = 0.0_Kr
      U%Y = 0.0_Kr
      U%Z = 0.0_Kr
      F%X = 0.0_Kr
      F%Y = 0.0_Kr
      F%Z = 0.0_Kr
      Theta = 0.0_Kr
      Do i = 1, MeshTopology%Num_Elem_Blks_Global
         Write(IOBuffer, 100) i
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

         !!! Displacement
         If (MyEXO%EBProperty(Rupt_EBProp_BCTypeX)%Value(i) /= 0 ) Then
            Write(IOBuffer, 200) 'Ux'
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
            If (MEF90_MyRank == 0) Then
               Read(*,*) U(i)%X
            End If
            Call MPI_BCast(U(i)%X, 1, MPIU_SCALAR, 0, EXO%Comm, iErr)
         End If
         If (MyEXO%EBProperty(Rupt_EBProp_BCTypeY)%Value(i) /= 0 ) Then
            Write(IOBuffer, 200) 'Uy'
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
            If (MEF90_MyRank == 0) Then
               Read(*,*) U(i)%Y
            End If
            Call MPI_BCast(U(i)%Y, 1, MPIU_SCALAR, 0, EXO%Comm, iErr)
         End If
         If (MyEXO%EBProperty(Rupt_EBProp_BCTypeZ)%Value(i) /= 0 ) Then
            Write(IOBuffer, 200) 'Uz'
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
            If (MEF90_MyRank == 0) Then
               Read(*,*) U(i)%Z
            End If
            Call MPI_BCast(U(i)%Z, 1, MPIU_SCALAR, 0, EXO%Comm, iErr)
         End If

         !!! Force
         If (MyEXO%EBProperty(Rupt_EBProp_HasBForce)%Value(i) /= 0 ) Then
            Write(IOBuffer, 200) 'Fx'
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
            If (MEF90_MyRank == 0) Then
               Read(*,*) F(i)%X
            End If
            Call MPI_BCast(F(i)%X, 1, MPIU_SCALAR, 0, EXO%Comm, iErr)

            Write(IOBuffer, 200) 'Fy'
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
            If (MEF90_MyRank == 0) Then
               Read(*,*) F(i)%Y
            End If
            Call MPI_BCast(F(i)%Y, 1, MPIU_SCALAR, 0, EXO%Comm, iErr)

            Write(IOBuffer, 200) 'Fz'
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
            If (MEF90_MyRank == 0) Then
               Read(*,*) F(i)%Z
            End If
            Call MPI_BCast(F(i)%Z, 1, MPIU_SCALAR, 0, EXO%Comm, iErr)
         End If
         
         !!! Temperature
         Write(IOBuffer, 200) 'Theta'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         If (MEF90_MyRank == 0) Then
            Read(*,*) Theta(i)
         End If
         Call MPI_BCast(Theta(i), 1, MPIU_SCALAR, 0, EXO%Comm, iErr)

         !!! Material Properties
         If (MEF90_MyRank == 0) Then
            Write(IOBuffer, 200) 'Toughness'
            Call PetscPrintf(PETSC_COMM_SELF, IOBuffer, iErr); CHKERRQ(iErr)
            Read(*,*) Toughness
            Write(IOBuffer, 200) 'Young modulus'
            Call PetscPrintf(PETSC_COMM_SELF, IOBuffer, iErr); CHKERRQ(iErr)
            Read(*,*) E
            Write(IOBuffer, 200) 'Poisson ratio'
            Call PetscPrintf(PETSC_COMM_SELF, IOBuffer, iErr); CHKERRQ(iErr)
            Read(*,*) nu
            Write(IOBuffer, 200) 'Thermal expansion coefficient'
            Call PetscPrintf(PETSC_COMM_SELF, IOBuffer, iErr); CHKERRQ(iErr)
            Read(*,*) Therm_ExpScal
            
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
         End If        

      End Do   
      
      Do iloc = 1, MeshTopology%Num_Elem_Blks
         i = MeshTopology%Elem_Blk(iLoc)%ID
         !!! Initialize the Section
         Num_DoF = MeshTopology%Elem_Blk(iloc)%Num_DoF
         Allocate(Uelem(3*Num_DoF))
         Allocate(Felem(3*Num_DoF))
         Allocate(Thetaelem(Num_DoF))
         
         !!! Update U 
         If ( (MyEXO%EBProperty(Rupt_EBProp_BCTypeX)%Value(i) /= 0 ) .OR. (MyEXO%EBProperty(Rupt_EBProp_BCTypeY)%Value(i) /= 0 ) .OR. (MyEXO%EBProperty(Rupt_EBProp_BCTypeZ)%Value(i) /= 0 ) ) Then
            Do k = 0, Num_DoF-1
               Uelem(3*k+1) = U(i)%X
               Uelem(3*k+2) = U(i)%Y
               Uelem(3*k+3) = U(i)%Z
            End Do
            Do j = 1, MeshTopology%Elem_Blk(iloc)%Num_Elems
               Call MeshUpdateClosure(MeshTopology%Mesh, USec, MeshTopology%Elem_Blk(iloc)%Elem_ID(j)-1, Uelem, iErr); CHKERRQ(iErr)            
            End Do
         End If
         
         !!! Update F
         If ( MyEXO%EBProperty(Rupt_EBProp_HasBForce)%Value(i) /= 0 ) Then
            Do k = 0, Num_DoF-1
               Felem(3*k+1) = F(i)%X
               Felem(3*k+2) = F(i)%Y
               Felem(3*k+3) = F(i)%Z
            End Do
            Do j = 1, MeshTopology%Elem_Blk(iloc)%Num_Elems
               Call MeshUpdateClosure(MeshTopology%Mesh, FSec, MeshTopology%Elem_Blk(iloc)%Elem_ID(j)-1, Felem, iErr); CHKERRQ(iErr)            
            End Do
         End If

         !!! Update Theta
         Thetaelem = Theta(i)
         Do j = 1, MeshTopology%Elem_Blk(iloc)%Num_Elems
            Call MeshUpdateClosure(MeshTopology%Mesh, ThetaSec, MeshTopology%Elem_Blk(iloc)%Elem_ID(j)-1, Thetaelem, iErr); CHKERRQ(iErr) 
         End Do
         DeAllocate(Uelem)
         DeAllocate(Felem)
         DeAllocate(Thetaelem)         
      End Do
      DeAllocate(U)
      DeAllocate(F)
      DeAllocate(Theta)
      
!!!     !!! NS BC and Variables
!!!      Do i = 1, MeshTopology%Num_Node_Sets_Global
!!!         Write(IOBuffer, 102) i
!!!         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
!!!
!!!         !!! Displacement
!!!         U     = 0.0_Kr
!!!         If (MyEXO%NSProperty(Rupt_NSProp_BCTypeX)%Value(i) /= 0 ) Then
!!!            Write(IOBuffer, 200) 'Ux'
!!!            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
!!!            If (MEF90_MyRank == 0) Then
!!!               Read(*,*) U%X
!!!            End If
!!!            Call MPI_BCast(U%X, 1, MPIU_SCALAR, 0, EXO%Comm, iErr)
!!!         End If
!!!         If (MyEXO%NSProperty(Rupt_NSProp_BCTypeY)%Value(i) /= 0 ) Then
!!!            Write(IOBuffer, 200) 'Uy'
!!!            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
!!!            If (MEF90_MyRank == 0) Then
!!!               Read(*,*) U%Y
!!!            End If
!!!            Call MPI_BCast(U%Y, 1, MPIU_SCALAR, 0, EXO%Comm, iErr)
!!!         End If
!!!         If (MyEXO%NSProperty(Rupt_NSProp_BCTypeZ)%Value(i) /= 0 ) Then
!!!            Write(IOBuffer, 200) 'Uz'
!!!            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
!!!            If (MEF90_MyRank == 0) Then
!!!               Read(*,*) U%Z
!!!            End If
!!!            Call MPI_BCast(U%Z, 1, MPIU_SCALAR, 0, EXO%Comm, iErr)
!!!         End If
!!!
!!!         !!! Force
!!!         F     = 0.0_Kr
!!!         If (MyEXO%NSProperty(Rupt_NSProp_HasPForce)%Value(i) /= 0 ) Then
!!!            Write(IOBuffer, 200) 'Fx'
!!!            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
!!!            If (MEF90_MyRank == 0) Then
!!!               Read(*,*) F%X
!!!            End If
!!!            Call MPI_BCast(F%X, 1, MPIU_SCALAR, 0, EXO%Comm, iErr)
!!!
!!!            Write(IOBuffer, 200) 'Fy'
!!!            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
!!!            If (MEF90_MyRank == 0) Then
!!!               Read(*,*) F%Y
!!!            End If
!!!            Call MPI_BCast(F%Y, 1, MPIU_SCALAR, 0, EXO%Comm, iErr)
!!!
!!!            Write(IOBuffer, 200) 'Fz'
!!!            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
!!!            If (MEF90_MyRank == 0) Then
!!!               Read(*,*) F%Z
!!!            End If
!!!            Call MPI_BCast(F%Z, 1, MPIU_SCALAR, 0, EXO%Comm, iErr)
!!!         End If
!!!         
!!!         !!! Initialize the Section
!!!         Allocate(Uelem(3))
!!!         Allocate(Felem(3))
!!!         
!!!         !!! Update 
!!!         Uelem(1) = U%X
!!!         Uelem(2) = U%Y
!!!         Uelem(3) = U%Z
!!!         Felem(1) = F%X
!!!         Felem(2) = F%Y
!!!         Felem(3) = F%Z
!!!
!!!         Do j = 1, MeshTopology%Node_Set(i)%Num_Nodes
!!!!            Call MeshUpdateClosure(MeshTopology%Mesh, USec, MeshTopology%Num_Elems + MeshTopology%Node_Set(i)%Node_ID(j)-1, Uelem, iErr); CHKERRQ(iErr)            
!!!!            Call MeshUpdateClosure(MeshTopology%Mesh, FSec, MeshTopology%Num_Elems + MeshTopology%Node_Set(i)%Node_ID(j)-1, Felem, iErr); CHKERRQ(iErr)            
!!!            !!! We will need to do a restrict then update here!
!!!         End Do
!!!         DeAllocate(Uelem)
!!!         DeAllocate(Felem)
!!!      End Do

!      Call SectionRealComplete(USec, iErr); CHKERRQ(iErr)
!      Call SectionRealComplete(FSec, iErr); CHKERRQ(iErr)
!      Call SectionRealComplete(ThetaSec, iErr); CHKERRQ(iErr)
      
      ! Save fields ( Sec -> Vec then VecScale then save )
      !That'll be fine for now
      step = 1
      Call Write_EXO_Result_Vertex(MyEXO, MeshTopology, MyEXO%VertVariable(Rupt_VertVar_DisplacementX)%Offset, step, USec) 
      Call Write_EXO_Result_Vertex(MyEXO, MeshTopology, MyEXO%VertVariable(Rupt_VertVar_ForceX)%Offset, step, FSec) 
      Call Write_EXO_Result_Vertex(MyEXO, MeshTopology, MyEXO%VertVariable(Rupt_VertVar_Temperature)%Offset, step, ThetaSec) 
      Call Write_EXO_Result_Global(MyEXO, MyEXO%GlobVariable(Rupt_GlobVar_Load)%Offset, step, T(step))
      
      Select Case(MeshTopology%Num_Dim)
      Case (2)
         Call MatProp_Write(MeshTopology, MatProp2D, Trim(prefix)//'.CST')
      Case(3)
         Call MatProp_Write(MeshTopology, MatProp3D, Trim(prefix)//'.CST')
      End Select

   Case (2)

   Case Default
      SETERRQ(PETSC_ERR_SUP, 'Unknown test case\n'c, iErr)      
   End Select


   Call MEF90_Finalize()

 100 Format('*** Element Block ', T24, I3, '\n'c)
 101 Format('*** Side Set      ', T24, I3, '\n'c)
 102 Format('*** Node Set      ', T24, I3, '\n'c)
 200 Format(A,t60, ': ')
   
End Program PrepRupt