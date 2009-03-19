Module m_RuptStruct

#include "finclude/petscdef.h"

   Use m_MEF90
   Use petsc

   Implicit NONE
   Private

   
   Public :: GenHL_Iso3D_Enu
   Public :: GenHL_Iso2D_EnuPlaneStress
   Public :: GenHL_Iso2D_EnuPlaneStrain
   Public :: GenHL_Ortho2D_LambdaMu
   
   Public :: RuptSchemeParam_View
   Public :: RuptSchemeParam_Load
   Public :: RuptSchemeParam_GetFromOptions
   
   Public :: RuptEXOProperty_Init
   Public :: RuptEXOVariable_Init
   
   Public :: MatProp2D_Type, MatProp3D_Type
   
   Public :: EXOProperty_InitBCFlag2DA

   Interface MatProp_Write
      Module Procedure MatProp2D_Write, MatProp3D_Write
   End Interface
  
   Interface MatProp_Read
      Module Procedure MatProp2D_Read, MatProp3D_Read
   End Interface
   
   Interface GenHL_Iso_LambdaMu
      Module Procedure GenHL_Iso2D_LambdaMu, GenHL_Iso3D_LambdaMu
   End Interface

   PetscInt, Parameter, Public                     :: BC_Type_NONE = 0
   PetscInt, Parameter, Public                     :: BC_Type_DIRI = 1

   PetscInt, Parameter, Public                     :: Init_V_PREV = 0
   PetscInt, Parameter, Public                     :: Init_V_ONE  = 1
   PetscInt, Parameter, Public                     :: Init_V_RND  = 2
   PetscInt, Parameter, Public                     :: Init_V_SPH  = 3
   
   PetscInt, Parameter, Public                     :: Init_U_PREV = 0
   PetscInt, Parameter, Public                     :: Init_U_ZERO = 1
   
   PetscInt, Parameter, Public                     :: Irrev_NONE = 0
   PetscInt, Parameter, Public                     :: irrev_eq   = 1
   PetscInt, Parameter, Public                     :: Irrev_Ineq = 2
   
   PetscInt, Parameter, Public                     :: Rupt_Num_VertVar           = 8
   PetscInt, Parameter, Public                     :: Rupt_VertVar_Fracture      = 1
   PetscInt, Parameter, Public                     :: Rupt_VertVar_DisplacementX = 2   
   PetscInt, Parameter, Public                     :: Rupt_VertVar_DisplacementY = 3
   PetscInt, Parameter, Public                     :: Rupt_VertVar_DisplacementZ = 4   
   PetscInt, Parameter, Public                     :: Rupt_VertVar_ForceX        = 5   
   PetscInt, Parameter, Public                     :: Rupt_VertVar_ForceY        = 6
   PetscInt, Parameter, Public                     :: Rupt_VertVar_ForceZ        = 7   
   PetscInt, Parameter, Public                     :: Rupt_VertVar_Temperature   = 8
   
   PetscInt, Parameter, Public                     :: Rupt_Num_CellVar      = 12
   PetscInt, Parameter, Public                     :: Rupt_CellVar_StrainXX = 1
   PetscInt, Parameter, Public                     :: Rupt_CellVar_StrainYY = 2 
   PetscInt, Parameter, Public                     :: Rupt_CellVar_StrainZZ = 3
   PetscInt, Parameter, Public                     :: Rupt_CellVar_StrainXY = 4
   PetscInt, Parameter, Public                     :: Rupt_CellVar_StrainYZ = 5
   PetscInt, Parameter, Public                     :: Rupt_CellVar_StrainXZ = 6
   PetscInt, Parameter, Public                     :: Rupt_CellVar_StressXX = 7
   PetscInt, Parameter, Public                     :: Rupt_CellVar_StressYY = 8
   PetscInt, Parameter, Public                     :: Rupt_CellVar_StressZZ = 9
   PetscInt, Parameter, Public                     :: Rupt_CellVar_StressXY = 10
   PetscInt, Parameter, Public                     :: Rupt_CellVar_StressYZ = 11
   PetscInt, Parameter, Public                     :: Rupt_CellVar_StressZX = 12

   PetscInt, Parameter, Public                     :: Rupt_Num_GlobVar           = 5
   PetscInt, Parameter, Public                     :: Rupt_GlobVar_BulkEnergy    = 1
   PetscInt, Parameter, Public                     :: Rupt_GlobVar_SurfaceEnergy = 2 
   PetscInt, Parameter, Public                     :: Rupt_GlobVar_KineticEnergy = 3 
   PetscInt, Parameter, Public                     :: Rupt_GlobVar_TotalEnergy   = 4
   PetscInt, Parameter, Public                     :: Rupt_GlobVar_Load          = 5
   
   PetscInt, Parameter, Public                     :: Rupt_Num_EBProperties  = 6
   PetscInt, Parameter, Public                     :: Rupt_EBProp_IsBrittle  = 1
   PetscInt, Parameter, Public                     :: Rupt_EBProp_IsDomain   = 2
   PetscInt, Parameter, Public                     :: Rupt_EBProp_HasBForce  = 3
   PetscInt, Parameter, Public                     :: Rupt_EBProp_BCTypeX    = 4
   PetscInt, Parameter, Public                     :: Rupt_EBProp_BCTypeY    = 5
   PetscInt, Parameter, Public                     :: Rupt_EBProp_BCTypeZ    = 6
   
   PetscInt, Parameter, Public                     :: Rupt_Num_SSProperties = 4
   PetscInt, Parameter, Public                     :: Rupt_SSProp_BCTypeX   = 1
   PetscInt, Parameter, Public                     :: Rupt_SSProp_BCTypeY   = 2
   PetscInt, Parameter, Public                     :: Rupt_SSProp_BCTypeZ   = 3
   PetscInt, Parameter, Public                     :: Rupt_SSProp_HasSForce  = 4

   PetscInt, Parameter, Public                     :: Rupt_Num_NSProperties  = 4
   PetscInt, Parameter, Public                     :: Rupt_NSProp_BCTypeX    = 1
   PetscInt, Parameter, Public                     :: Rupt_NSProp_BCTypeY    = 2
   PetscInt, Parameter, Public                     :: Rupt_NSProp_BCTypeZ    = 3
   PetscInt, Parameter, Public                     :: Rupt_NSProp_HasPForce  = 4
   
   Type MatProp2D_Type
      PetscReal, Dimension(:), Pointer             :: Toughness
      Type(Tens4OS2D), Dimension(:), Pointer       :: Hookes_Law
      PetscReal, Dimension(:), Pointer             :: Therm_Exp      
   End Type MatProp2D_Type
   
   Type MatProp3D_Type
      PetscReal, Dimension(:), Pointer             :: Toughness
      Type(Tens4OS3D), Dimension(:), Pointer       :: Hookes_Law
      PetscReal, Dimension(:), Pointer             :: Therm_Exp      
   End Type MatProp3D_Type
   
   Type RuptSchemeParam_Type
      PetscInt                                     :: DoIrrev
      PetscReal                                    :: IrrevTol
      
      PetscTruth                                   :: DoBT
      PetscReal                                    :: BTTol
      PetscInt                                     :: BTInt
      
      PetscInt                                     :: InitU
      
      PetscInt                                     :: InitV
      PetscInt                                     :: nbCracks
      PetscReal                                    :: MaxCrackLength     
      
      PetscInt                                     :: RelaxMaxIter
      PetscReal                                    :: RelaxTol

      PetscReal                                    :: KSPUrtol
      PetscReal                                    :: KSPVrtol
      
      PetscInt                                     :: SaveInt
      
      PetscReal                                    :: Epsilon
      PetscReal                                    :: KEpsilon

      PetscInt                                     :: ATNum
   End Type RuptSchemeParam_Type
   
 Contains
   Subroutine EXOProperty_InitBCFlag2DA(dEXO, dMeshTopology, dBCFlag)
      Type(EXO_Type)                               :: dEXO
      Type(MeshTopology_Type)                      :: dMeshTopology
      Type(SectionInt)                             :: dBCFlag 
      
      PetscInt                                     :: iErr, NumDoF, i, j
      PetscInt, Dimension(:), Pointer              :: Flag
      
      Call SectionIntZero(dBCFlag, iErr); CHKERRQ(iErr)
      !!! Element Blocks
      Do i = 1, dMeshTopology%Num_Elem_Blks
         NumDoF = dMeshTopology%Elem_Blk(i)%Num_DoF
         Allocate(Flag(NumDoF))
         Flag = dEXO%EBProperty( Rupt_EBProp_BCTypeZ )%Value( dMeshTopology%Elem_Blk(i)%ID )
         Do j = 1, dMeshTopology%Elem_Blk(i)%Num_Elems
            Call MeshUpdateAddClosureInt(dMeshTopology%Mesh, dBCFlag, dMeshTopology%Elem_Blk(i)%Elem_ID(j)-1, Flag, iErr); CHKERRQ(iErr)
         End Do
         DeAllocate(Flag)
      End Do
      
      !!! Side Sets
      !!! To be implemented
      
      !!! Node Sets
      Do i = 1, dMeshTopology%Num_Node_Sets
         Allocate(Flag(1))
         Flag = dEXO%NSProperty( Rupt_NSProp_BCTypeZ )%Value( dMeshTopology%Node_Set(i)%ID )
         Do j = 1, dMeshTopology%Node_Set(i)%Num_Nodes
            Print*, MEF90_MyRank, i, dMeshTopology%Node_Set(i)%ID, j, dMeshTopology%Node_Set(i)%Node_ID(j), Flag
            Call MeshUpdateAddClosureInt(dMeshTopology%Mesh, dBCFlag, dMeshTopology%Node_Set(i)%Node_ID(j) + dMeshTopology%Num_Elems-1, Flag, iErr); CHKERRQ(iErr)
         End Do
         DeAllocate(Flag)
      End Do
      Call SectionIntComplete(dBCFlag, iErr); CHKERRQ(iErr)
   End Subroutine EXOProperty_InitBCFlag2DA
   
   Subroutine MatProp2D_Write(MeshTopology, MatProp, filename)
      Type(MeshTopology_Type)                      :: MeshTopology
      Type(MatProp2D_Type), Dimension(:), Pointer  :: MatProp
      Character(len=*)                             :: filename
      PetscMPIInt                                  :: rank
      PetscInt                                     :: iBlk, Blk_ID
      
      Open(File = filename, Unit = F_OUT, Status = 'Unknown')
      Rewind(F_OUT)
      Write(F_OUT, 110) MeshTopology%Num_Elem_Blks
      Do iBlk = 1, Size(MatProp)
         Blk_ID = MeshTopology%Elem_Blk(iBlk)%ID
         Write(F_OUT,120) Blk_ID, MatProp(iBlk)%Toughness, MatProp(iBlk)%Hookes_Law, MatProp(iBlk)%Therm_Exp
      End Do
      Close(F_OUT)
      
110   Format(I6,'     Toughness    A1111        A1112        A1122        A1212        A1222        A2222        Alpha')
120   Format(I6, 8(ES12.5,' '))   
   End Subroutine MatProp2D_Write
 
   Subroutine MatProp3D_Write(MeshTopology, MatProp, filename)
      Type(MeshTopology_Type)                      :: MeshTopology
      Type(MatProp3D_Type), Dimension(:), Pointer  :: MatProp
      Character(len=*)                             :: filename

      PetscMPIInt                                  :: rank
      PetscInt                                     :: iBlk, Blk_ID
      
      Open(File = filename, Unit = F_OUT, Status = 'Unknown')
      Rewind(F_OUT)
      Write(F_OUT, 110) MeshTopology%Num_Elem_Blks
      Do iBlk = 1, Size(MatProp)
         Blk_ID = MeshTopology%Elem_Blk(iBlk)%ID
         Write(F_OUT,120) Blk_ID, MatProp(iBlk)%Toughness, MatProp(iBlk)%Hookes_Law, MatProp(iBlk)%Therm_Exp
      End Do
      Close(F_OUT)
      
110   Format(I6,' Toughness    A_1111       A_1112       A_1113       A_1122       A_1123       A_1133       A_1212       A_1213       A_1222       A_1223       A_12133      A_1313       A_1322       A_1323       A_1333       A_2222       A_2223       A_2233       A_2323       A_2333       A_3333       Alpha')
120   Format(I6, 23(ES12.5,' '))
   End Subroutine MatProp3D_Write
 
 
   Subroutine MatProp2D_Read(MeshTopology, MatProp, filename)
      Type(MeshTopology_Type)                      :: MeshTopology
      Type(MatProp2D_Type), Dimension(:), Pointer  :: MatProp
      Character(len=*)                             :: filename

      PetscInt                                     :: iBlk, iErr
      
      PetscInt                                     :: NumBlks, IdxMin, IdxMax, Idx
      Type(Tens4OS2D), Dimension(:), Pointer       :: Hookes_Law
      PetscReal                                    :: Toughness, Therm_Exp
   
      Open(File = filename, Unit = F_IN, Status = 'Unknown')
      Rewind(F_IN)
      Read(F_IN, *) NumBlks
      If (NumBlks /= MeshTopology%Num_Elem_Blks) Then
         SETERRQ(PETSC_ERR_ARG_SIZ, 'MatProp2DRead: non matching blocks numbers', iErr)
      End If
      !!! Reading the file once first to get the right number of blocks
      IdxMin = 0
      IdxMax = 0
      Do iBlk = 1, NumBlks
         Read(F_IN, *) Idx
         !!! Check that this will work!
         IdxMin = Min(IdxMin, Idx)
         IdxMax = Max(IdxMax, Idx)
      End Do
      Allocate(MatProp(IdxMin:IdxMax))
      Rewind(F_IN)
      Do iBlk = 1, NumBlks
         Read(F_IN, 120) Idx, Toughness, Hookes_Law, Therm_exp
         MatProp(Idx)%Toughness  = Toughness
         MatProp(Idx)%Hookes_Law = Hookes_Law
         MatProp(Idx)%Therm_Exp  = Therm_Exp
      End Do
      Close(F_IN)
      
120   Format(I6, 8(ES12.5,' '))   
   End Subroutine MatProp2D_Read
   
   Subroutine MatProp3D_Read(MeshTopology, MatProp, filename)
      Type(MeshTopology_Type)                      :: MeshTopology
      Type(MatProp3D_Type), Dimension(:), Pointer  :: MatProp
      Character(len=*)                             :: filename

      PetscInt                                     :: iBlk, Blk_Id, iErr
      
      PetscInt                                     :: NumBlks, IdxMin, IdxMax, Idx
      Type(Tens4OS3D), Dimension(:), Pointer       :: Hookes_Law
      PetscReal                                    :: Toughness, Therm_Exp
   
      Open(File = filename, Unit = F_IN, Status = 'Unknown')
      Rewind(F_IN)
      Read(F_IN, *) NumBlks
      If (NumBlks /= MeshTopology%Num_Elem_Blks) Then
         SETERRQ(PETSC_ERR_ARG_SIZ, 'MatProp3DRead: non matching blocks numbers', iErr)
      End If
      !!! Reading the file once first to get the right number of blocks
      IdxMin = 0
      IdxMax = 0
      Do iBlk = 1, NumBlks
         Read(F_IN, *) Idx
         !!! Check that this will work!
         IdxMin = Min(IdxMin, Idx)
         IdxMax = Max(IdxMax, Idx)
      End Do
      Allocate(MatProp(IdxMin:IdxMax))
      Rewind(F_IN)
      Do iBlk = 1, NumBlks
         Read(F_IN, 120) Idx, Toughness, Hookes_Law, Therm_exp
         MatProp(Idx)%Toughness  = Toughness
         MatProp(Idx)%Hookes_Law = Hookes_Law
         MatProp(Idx)%Therm_Exp  = Therm_Exp
      End Do
      Close(F_IN)
      
120   Format(I6, 23(ES12.5,' '))
   End Subroutine MatProp3D_Read

   Subroutine RuptSchemeParam_View(dSchemeParam, viewer)
      Type(RuptSchemeParam_Type)                   :: dSchemeParam
      Type(PetscViewer)                            :: viewer
      PetscInt                                     :: iErr
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
      
      Write(IOBuffer, "(I1,T32, 'DoIrrev')")             dSchemeParam%DoIrrev 
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(ES12.5,T32, 'IrrevTol')")        dSchemeParam%IrrevTol
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(L1,T32, 'DoBT')")                dSchemeParam%DoBT 
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(ES12.5,T32, 'BTTol')")           dSchemeParam%BTTol
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(I5,T32, 'BTInt')")               dSchemeParam%BTInt 
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(I1,T32, 'InitU')")               dSchemeParam%InitU 
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(I1,T32, 'InitV')")               dSchemeParam%InitV 
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(I5,T32, 'NbCracks')")            dSchemeParam%NbCracks
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(ES12.5,T32, 'MaxCrackLength')")  dSchemeParam%MaxCrackLength
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(I5,T32, 'RelaxMaxIter')")        dSchemeParam%RelaxMaxIter
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(ES12.5,T32, 'RelaxTol')")        dSchemeParam%RelaxTol
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(ES12.5,T32, 'KSPUrTol')")        dSchemeParam%KSPUrTol
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(ES12.5,T32, 'KSPVrTol')")        dSchemeParam%KSPVrTol
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(I5,T32, 'SaveInt')")             dSchemeParam%SaveInt
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(ES12.5,T32, 'Epsilon')")         dSchemeParam%Epsilon
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(ES12.5,T32, 'KEpsilon')")        dSchemeParam%KEpsilon
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(I1,T32, 'ATNum')")               dSchemeParam%ATNum
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
   End Subroutine RuptSchemeParam_View

   Subroutine RuptSchemeParam_Load(dSchemeParam, filename)
      Type(RuptSchemeParam_Type)                   :: dSchemeParam
      Character(len=*)                             :: filename
      
      Open(File = filename, status='old', Unit = F_IN)
      Rewind(F_IN)
      Read(F_IN, *) dSchemeParam%DoIrrev 
      Read(F_IN, *) dSchemeParam%IrrevTol
      Read(F_IN, *) dSchemeParam%DoBT 
      Read(F_IN, *) dSchemeParam%BTTol
      Read(F_IN, *) dSchemeParam%BTInt 
      Read(F_IN, *) dSchemeParam%InitU 
      Read(F_IN, *) dSchemeParam%InitV 
      Read(F_IN, *) dSchemeParam%NbCracks
      Read(F_IN, *) dSchemeParam%MaxCrackLength
      Read(F_IN, *) dSchemeParam%RelaxMaxIter
      Read(F_IN, *) dSchemeParam%RelaxTol
      Read(F_IN, *) dSchemeParam%KSPUrTol
      Read(F_IN, *) dSchemeParam%KSPVrTol
      Read(F_IN, *) dSchemeParam%SaveInt
      Read(F_IN, *) dSchemeParam%Epsilon
      Read(F_IN, *) dSchemeParam%KEpsilon
      Read(F_IN, *) dSchemeParam%ATNum
      Close(F_IN)
   End Subroutine RuptSchemeParam_Load
   
   Subroutine RuptSchemeParam_GetFromOptions(dSchemeParam)
      Type(RuptSchemeParam_Type)                   :: dSchemeParam
      PetscInt                                     :: iErr

      dSchemeParam%DoIrrev        = Irrev_Eq
      dSchemeParam%IrrevTol       = 1.0D-2
      dSchemeParam%DoBT           = PETSC_FALSE
      dSchemeParam%BTTol          = 0.0D0
      dSchemeParam%BTInt          = 0
      dSchemeParam%InitU          = Init_U_PREV
      dSchemeParam%InitV          = Init_V_PREV
      dSchemeParam%nbCracks       = 0
      dSchemeParam%MaxCrackLength = 0.0D0  
      dSchemeParam%RelaxMaxIter   = 1000
      dSchemeParam%RelaxTol       = 1.0D-4
      dSchemeParam%KSPUrtol       = 1.0D-6
      dSchemeParam%KSPVrtol       = 1.0D-6
      dSchemeParam%SaveInt        = 25
      dSchemeParam%Epsilon        = .1
      dSchemeParam%KEpsilon       = 1.0E-6
      dSchemeParam%ATNum          = 2

      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-doirrev',        dSchemeParam%DoIrrev, iErr); CHKERRQ(iErr) 
      Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,  '-irrevtol',       dSchemeParam%IrrevTol, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetTruth(PETSC_NULL_CHARACTER, '-dobt',           dSchemeParam%DoBT, iErr); CHKERRQ(iErr) 
      Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,  '-bttol',          dSchemeParam%BTTol, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-btint',          dSchemeParam%BTInt, iErr); CHKERRQ(iErr) 
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-initu',          dSchemeParam%InitU, iErr); CHKERRQ(iErr) 
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-initv',          dSchemeParam%InitV, iErr); CHKERRQ(iErr) 
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-nbcracks',       dSchemeParam%NbCracks, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,  '-maxcracklength', dSchemeParam%MaxCrackLength, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-relaxmaxiter',   dSchemeParam%RelaxMaxIter, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,  '-relaxtol',       dSchemeParam%RelaxTol, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,  '-kspurol',        dSchemeParam%KSPUrTol, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,  '-kspvrol',        dSchemeParam%KSPVrTol, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-saveint',        dSchemeParam%SaveInt, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,  '-epsilon',        dSchemeParam%Epsilon, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,  '-kepsilon',       dSchemeParam%KEpsilon, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-atum',           dSchemeParam%ATNum, iErr); CHKERRQ(iErr)
   End Subroutine RuptSchemeParam_GetFromOptions
   
   Subroutine RuptEXOProperty_Init(dEXO)
      Type(EXO_Type)                      :: dEXO
      PetscInt                            :: i, vers, iErr
      PetscInt                            :: NumEB, NumSS, NumNS

      Integer                             :: EXO_MyRank
      PetscReal                           :: rDummy
      Character                           :: cDummy
          

      Call MPI_COMM_RANK(dEXO%Comm, EXO_MyRank, iErr)

      If (EXO_MyRank == 0) Then
         dEXO%exoid = EXOPEN(dEXO%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers, ierr)
         
         !!! This is ridiculous... 
         !!! This won't work if the mesh has not be written to disk yet...

         Call EXINQ(dEXO%exoid, EXELBL, NumEB, rDummy, cDummy, iErr)
         Call EXINQ(dEXO%exoid, EXSIDS, NumSS, rDummy, cDummy, iErr)
         Call EXINQ(dEXO%exoid, EXNODS, NumNS, rDummy, cDummy, iErr)
         
         If ( (NumEB == 0) .AND. (NumSS == 0) .AND. (NumSS ==0) ) Then
            Call PetscPrintf(PETSC_COMM_SELF, '[WARNING]: The EXO file contains no EB, SS or NS is this right?\n'c, iErr); CHKERRQ(iErr)
            Call PetscPrintf(PETSC_COMM_SELF, '           Was Write_MeshTopologyGlobal called before RuptEXOProperty_Init?\n'c, iErr); CHKERRQ(iErr)
         End If
         Call EXCLOS(dEXO%exoid, iErr)
         dEXO%exoid = 0
      End If

      Call MPI_BCast(NumEB, 1, MPI_INTEGER, 0, dEXO%Comm, iErr)
      Call MPI_BCast(NumSS, 1, MPI_INTEGER, 0, dEXO%Comm, iErr)
      Call MPI_BCast(NumNS, 1, MPI_INTEGER, 0, dEXO%Comm, iErr)
      
      dEXO%Num_EBProperties = Rupt_Num_EBProperties
      Allocate(dEXO%EBProperty(dEXO%Num_EBProperties))
      dEXO%EBProperty(Rupt_EBProp_IsBrittle)%Name = 'Is_Brittle'
      dEXO%EBProperty(Rupt_EBProp_IsDomain)%Name  = 'Is_Domain'
      dEXO%EBProperty(Rupt_EBProp_HasBForce)%Name = 'Has_BForce'
      dEXO%EBProperty(Rupt_EBProp_BCTypeX)%Name   = 'BC_Type_X'
      dEXO%EBProperty(Rupt_EBProp_BCTypeY)%Name   = 'BC_Type_Y'
      dEXO%EBProperty(Rupt_EBProp_BCTypeZ)%Name   = 'BC_Type_Z'
      Do i = 1, dEXO%Num_EBProperties
         Allocate(dEXO%EBProperty(i)%Value(NumEB))
         dEXO%EBProperty(i)%Value = 0
      End Do
      
      dEXO%Num_SSProperties = Rupt_Num_SSProperties
      Allocate(dEXO%SSProperty(dEXO%Num_SSProperties))
      dEXO%SSProperty(Rupt_SSProp_BCTypeX)%Name   = 'BC_Type_X'
      dEXO%SSProperty(Rupt_SSProp_BCTypeY)%Name   = 'BC_Type_Y'
      dEXO%SSProperty(Rupt_SSProp_BCTypeZ)%Name   = 'BC_Type_Z'
      dEXO%SSProperty(Rupt_SSProp_HasSForce)%Name = 'Has_SForce'
      Do i = 1, dEXO%Num_SSProperties
         Allocate(dEXO%SSProperty(i)%Value(NumSS))
         dEXO%SSProperty(i)%Value = 0
      End Do
      
      dEXO%Num_NSProperties = Rupt_Num_NSProperties
      Allocate(dEXO%NSProperty(dEXO%Num_NSProperties))
      dEXO%NSProperty(Rupt_NSProp_BCTypeX)%Name   = 'BC_Type_X'
      dEXO%NSProperty(Rupt_NSProp_BCTypeY)%Name   = 'BC_Type_Y'
      dEXO%NSProperty(Rupt_NSProp_BCTypeZ)%Name   = 'BC_Type_Z'
      dEXO%NSProperty(Rupt_NSProp_HasPForce)%Name = 'Has_PForce'
      Do i = 1, dEXO%Num_NSProperties
         Allocate(dEXO%NSProperty(i)%Value(NumNS))
         dEXO%NSProperty(i)%Value = 0
      End Do
   End Subroutine RuptEXOProperty_Init   

   Subroutine RuptEXOVariable_Init(dEXO)
      Type(EXO_Type)                      :: dEXO
      PetscInt                            :: i
      
      dEXO%Num_GlobVariables = Rupt_Num_GlobVar
      Allocate(dEXO%GlobVariable(dEXO%Num_GlobVariables))
      dEXO%GlobVariable(Rupt_GlobVar_BulkEnergy)%Name    = 'Bulk energy'
      dEXO%GlobVariable(Rupt_GlobVar_SurfaceEnergy)%Name = 'Surface energy'
      dEXO%GlobVariable(Rupt_GlobVar_KineticEnergy)%Name = 'Kinetic energy'
      dEXO%GlobVariable(Rupt_GlobVar_TotalEnergy)%Name   = 'Total energy'
      dEXO%GlobVariable(Rupt_GlobVar_Load)%Name          = 'Load'
      dEXO%GlobVariable(:)%Offset = (/ (i, i=1,dEXO%Num_GlobVariables) /)
      
      dEXO%Num_CellVariables = Rupt_Num_CellVar
      Allocate(dEXO%CellVariable(dEXO%Num_CellVariables))
      dEXO%CellVariable(Rupt_CellVar_StrainXX)%Name = 'Strain XX'
      dEXO%CellVariable(Rupt_CellVar_StrainYY)%Name = 'Strain YY' 
      dEXO%CellVariable(Rupt_CellVar_StrainZZ)%Name = 'Strain ZZ'
      dEXO%CellVariable(Rupt_CellVar_StrainXY)%Name = 'Strain XY'
      dEXO%CellVariable(Rupt_CellVar_StrainYZ)%Name = 'Strain YZ'
      dEXO%CellVariable(Rupt_CellVar_StrainXZ)%Name = 'Strain XZ'
      dEXO%CellVariable(Rupt_CellVar_StressXX)%Name = 'Stress XX'
      dEXO%CellVariable(Rupt_CellVar_StressYY)%Name = 'Stress YY'
      dEXO%CellVariable(Rupt_CellVar_StressZZ)%Name = 'Stress ZZ'
      dEXO%CellVariable(Rupt_CellVar_StressXY)%Name = 'Stress XY'
      dEXO%CellVariable(Rupt_CellVar_StressYZ)%Name = 'Stress YZ'
      dEXO%CellVariable(Rupt_CellVar_StressZX)%Name = 'Stress ZX'
      dEXO%CellVariable(:)%Offset = (/ (i, i=1,dEXO%Num_CellVariables) /)

      dEXO%Num_VertVariables = Rupt_Num_VertVar
      Allocate(dEXO%VertVariable(dEXO%Num_VertVariables))
      dEXO%VertVariable(Rupt_VertVar_Fracture)%Name      = 'Fracture'
      dEXO%VertVariable(Rupt_VertVar_DisplacementX)%Name = 'Displacement X'   
      dEXO%VertVariable(Rupt_VertVar_DisplacementY)%Name = 'Displacement Y'
      dEXO%VertVariable(Rupt_VertVar_DisplacementZ)%Name = 'Displacement Z'   
      dEXO%VertVariable(Rupt_VertVar_ForceX)%Name        = 'Force X'   
      dEXO%VertVariable(Rupt_VertVar_ForceY)%Name        = 'Force Y'
      dEXO%VertVariable(Rupt_VertVar_ForceZ)%Name        = 'Force Z'   
      dEXO%VertVariable(Rupt_VertVar_Temperature)%Name   = 'Temperature'
      dEXO%VertVariable(:)%Offset = (/ (i, i=1,dEXO%Num_VertVariables) /)
   End Subroutine RuptEXOVariable_Init
  
   Subroutine GenHL_Iso2D_LambdaMu(lambda, mu, A) 
      PetscReal, Intent(IN)               :: lambda, mu
      Type(Tens4OS2D), Intent(OUT)        :: A
      A = 0.0_Kr
      
      A%XXXX = lambda + 2.0_Kr * mu
      A%XXYY = lambda
      A%XYXY = mu
      A%YYYY = lambda + 2.0_Kr * mu
   End Subroutine GenHL_Iso2D_LambdaMu         

   Subroutine GenHL_Iso2D_EnuPlaneStress(E, nu, A) 
      PetscReal, Intent(IN)               :: E, nu
      Type(Tens4OS2D), Intent(OUT)        :: A
      
      PetscReal                           :: Lambda, mu
      
      lambda = E * nu / (1.0_Kr - nu**2) 
      mu     = E / (1.0_Kr + nu) * .5_Kr
      A = 0.0_Kr
      A%XXXX = lambda + 2.0_Kr * mu
      A%XXYY = lambda
      A%XYXY = mu
      A%YYYY = lambda + 2.0_Kr * mu
   End Subroutine GenHL_Iso2D_EnuPlaneStress         
   
   Subroutine GenHL_Iso2D_EnuPlaneStrain(E, nu, A) 
      PetscReal, Intent(IN)               :: E, nu
      Type(Tens4OS2D), Intent(OUT)        :: A
      
      PetscReal                           :: Lambda, mu
      
      lambda = E * nu / (1.0_Kr + nu) / (1.0_Kr - 2.0_Kr * nu)
      mu     = E / (1.0_Kr + nu) * .5_Kr      
      A = 0.0_Kr
      A%XXXX = lambda + 2.0_Kr * mu
      A%XXYY = lambda
      A%XYXY = mu
      A%YYYY = lambda + 2.0_Kr * mu
   End Subroutine GenHL_Iso2D_EnuPlaneStrain         

   Subroutine GenHL_Iso3D_LambdaMu(lambda, mu, A)
      PetscReal, Intent(IN)               :: Lambda, Mu
      Type(Tens4OS3D), Intent(OUT)        :: A
   
      A = 0.0_Kr
      A%XXXX = lambda + mu * 2.0_Kr
      A%XXYY = lambda
      A%XXZZ = lambda
      
      A%XYXY = mu
      
      A%XZXZ = mu
      
      A%YYYY = lambda + mu * 2.0_Kr
      A%YYZZ = lambda
      
      A%YZYZ = mu
      
      A%ZZZZ = lambda + mu * 2.0_Kr
   End Subroutine GenHL_Iso3D_LambdaMu

   Subroutine GenHL_Iso3D_Enu(E, nu, A)
      PetscReal, Intent(IN)               :: E, nu
      Type(Tens4OS3D), Intent(OUT)        :: A
      
      Real(Kind = Kr)                     :: Lambda, mu
   
      lambda = E * nu / (1.0_Kr + nu) / (1 - 2.0_Kr * nu)
      mu     = E / (1.0_Kr + nu) * .5_Kr      
   
      A = 0.0_Kr
      A%XXXX = lambda + mu * 2.0_Kr
      A%XXYY = lambda
      A%XXZZ = lambda
      
      A%XYXY = mu
      
      A%XZXZ = mu
      
      A%YYYY = lambda + mu * 2.0_Kr
      A%YYZZ = lambda
      
      A%YZYZ = mu
      
      A%ZZZZ = lambda + mu * 2.0_Kr
   End Subroutine GenHL_Iso3D_Enu

   Subroutine GenHL_Ortho2D_LambdaMu(lambda, mu1, mu2, theta, A)
      PetscReal, Intent(IN)               :: Lambda, mu1, mu2, theta
      Type(Tens4OS2D), Intent(OUT)        :: A
      
      A = 0.0_Kr
      A%XXXX = lambda + mu1 * (1.0_Kr + (cos(theta))**2) +  mu2 * (sin(theta))**2
      A%XXXY = (mu1-mu2) * cos(theta) * sin(theta)
      A%XXYY = lambda + (mu1-mu2) * (sin(theta))**2

      A%XYXY = mu1 * (sin(theta))**2 + mu2 * (cos(theta))**2
      A%XYYY = -A%XXXY

      A%YYYY =  A%XXXX
   End Subroutine GenHL_Ortho2D_LambdaMu
End Module m_RuptStruct
