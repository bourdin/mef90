Module m_VarFracFilm_Struct

#include "finclude/petscdef.h"

   Use m_MEF90
   Use petsc

   Implicit NONE
   Private

   
   Public :: GenHL_Iso2D_EnuPlaneStress
   Public :: GenHL_Iso2D_EnuPlaneStrain
   Public :: GenHL_Ortho2D_LambdaMu
   
   Public :: VarFracFilmSchemeParam_View
   Public :: VarFracFilmSchemeParam_Load
   Public :: VarFracFilmSchemeParam_GetFromOptions
   
   Public :: VarFracFilmEXOProperty_Init
   Public :: VarFracFilmEXOVariable_Init
   
   Public :: MatProp2D_Type
   Public :: MatProp_Write, MatProp_Read
   
   Public :: EXOProperty_InitBCVFlag
   Public :: EXOProperty_InitBCUFlag2D
   
   Public :: VarFracFilmSchemeParam_Type

   Interface MatProp_Write
      Module Procedure MatProp2D_Write
   End Interface
  
   Interface MatProp_Read
      Module Procedure MatProp2D_Read
   End Interface
   
   Interface GenHL_Iso_LambdaMu
      Module Procedure GenHL_Iso2D_LambdaMu
   End Interface

   PetscInt, Parameter, Public                     :: BC_Type_NONE = 0
   PetscInt, Parameter, Public                     :: BC_Type_DIRI = 1

   PetscInt, Parameter, Public                     :: Init_V_PREV = 0
   PetscInt, Parameter, Public                     :: Init_V_ONE  = 1
   PetscInt, Parameter, Public                     :: Init_V_RND  = 2
   PetscInt, Parameter, Public                     :: Init_V_SPH  = 3
   
   PetscInt, Parameter, Public                     :: Init_PHI_PREV = 0
   PetscInt, Parameter, Public                     :: Init_PHI_ONE  = 1
   PetscInt, Parameter, Public                     :: Init_PHI_RND  = 2
   PetscInt, Parameter, Public                     :: Init_PHI_SPH  = 3
   
   PetscInt, Parameter, Public                     :: Init_U_PREV = 0
   PetscInt, Parameter, Public                     :: Init_U_ZERO = 1
   
   PetscInt, Parameter, Public                     :: Irrev_NONE = 0
   PetscInt, Parameter, Public                     :: Irrev_Eq   = 1
   PetscInt, Parameter, Public                     :: Irrev_Ineq = 2
   
   PetscInt, Parameter, Public                     :: VarFracFilm_Num_VertVar           = 6
   PetscInt, Parameter, Public                     :: VarFracFilm_VertVar_Fracture      = 1
   PetscInt, Parameter, Public                     :: VarFracFilm_VertVar_DisplacementX = 2   
   PetscInt, Parameter, Public                     :: VarFracFilm_VertVar_DisplacementY = 3
   PetscInt, Parameter, Public                     :: VarFracFilm_VertVar_U0X           = 4
   PetscInt, Parameter, Public                     :: VarFracFilm_VertVar_U0Y           = 5
   PetscInt, Parameter, Public                     :: VarFracFilm_VertVar_Temperature   = 6
   
   PetscInt, Parameter, Public                     :: VarFracFilm_Num_CellVar          = 7
   PetscInt, Parameter, Public                     :: VarFracFilm_CellVar_Delamination = 1
   PetscInt, Parameter, Public                     :: VarFracFilm_CellVar_StrainXX     = 2
   PetscInt, Parameter, Public                     :: VarFracFilm_CellVar_StrainYY     = 3 
   PetscInt, Parameter, Public                     :: VarFracFilm_CellVar_StrainXY     = 4
   PetscInt, Parameter, Public                     :: VarFracFilm_CellVar_StressXX     = 5
   PetscInt, Parameter, Public                     :: VarFracFilm_CellVar_StressYY     = 6
   PetscInt, Parameter, Public                     :: VarFracFilm_CellVar_StressXY     = 7

   PetscInt, Parameter, Public                     :: VarFracFilm_Num_GlobVar                    = 6       
   PetscInt, Parameter, Public                     :: VarFracFilm_GlobVar_ElasticBulkEnergy      = 1   
   PetscInt, Parameter, Public                     :: VarFracFilm_GlobVar_ElasticInterEnergy     = 2
   PetscInt, Parameter, Public                     :: VarFracFilm_GlobVar_SurfaceEnergyT         = 3      
   PetscInt, Parameter, Public                     :: VarFracFilm_GlobVar_SurfaceEnergyD         = 4      
   PetscInt, Parameter, Public                     :: VarFracFilm_GlobVar_TotalEnergy            = 5       
   PetscInt, Parameter, Public                     :: VarFracFilm_GlobVar_Load                   = 6       
   
   PetscInt, Parameter, Public                     :: VarFracFilm_Num_EBProperties    = 3
   PetscInt, Parameter, Public                     :: VarFracFilm_EBProp_IsBrittle    = 1
   PetscInt, Parameter, Public                     :: VarFracFilm_EBProp_HasSubstrate = 2
   PetscInt, Parameter, Public                     :: VarFracFilm_EBProp_Elem_Type    = 3
   
   PetscInt, Parameter, Public                     :: VarFracFilm_Num_SSProperties  = 4
   PetscInt, Parameter, Public                     :: VarFracFilm_SSProp_BCUTypeX   = 1
   PetscInt, Parameter, Public                     :: VarFracFilm_SSProp_BCUTypeY   = 2
   PetscInt, Parameter, Public                     :: VarFracFilm_SSProp_BCVType    = 3
   PetscInt, Parameter, Public                     :: VarFracFilm_SSProp_Elem_Type  = 4

   PetscInt, Parameter, Public                     :: VarFracFilm_Num_NSProperties  = 3
   PetscInt, Parameter, Public                     :: VarFracFilm_NSProp_BCUTypeX   = 1
   PetscInt, Parameter, Public                     :: VarFracFilm_NSProp_BCUTypeY   = 2
   PetscInt, Parameter, Public                     :: VarFracFilm_NSProp_BCVType    = 3
   
   Type MatProp2D_Type
      PetscReal                                    :: ToughnessT     ! 1 value
      PetscReal                                    :: ToughnessD     ! 1 value
      Type(Tens4OS2D)                              :: Hookes_Law     ! 6 values
      Type(MatS2D)                                 :: Therm_Exp      ! 3 values
      Type(MatS2D)                                 :: K_interface    ! 3 values
   End Type MatProp2D_Type                                           ! Total: 14 PetscReals
      
   Type VarFracFilmSchemeParam_Type
      PetscInt                                     :: IrrevType
      PetscReal                                    :: IrrevTol
      
      PetscTruth                                   :: DoBT
      PetscReal                                    :: BTTol
      PetscInt                                     :: BTInt
      
      PetscInt                                     :: InitV
      PetscInt                                     :: InitPHI
      PetscInt                                     :: nbCracks
      PetscReal                                    :: MaxCrackLength     
      
      PetscInt                                     :: AltMinMaxIter
      PetscReal                                    :: AltMinTol
      PetscInt                                     :: AltMinSaveInt

      PetscReal                                    :: KSPUrtol
      PetscReal                                    :: KSPVrtol
      
      
      PetscReal                                    :: Epsilon
      PetscReal                                    :: KEpsilonV
      PetscReal                                    :: KEpsilonPhi
      
      PetscInt                                     :: ATNum
      PetscInt                                     :: IntegOrder
      
      PetscTruth                                   :: SaveStress
      PetscTruth                                   :: SaveStrain
   End Type VarFracFilmSchemeParam_Type
   
 Contains
    Subroutine EXOProperty_InitBCVFlag(dEXO, dMeshTopology, dBCFlag)
      Type(EXO_Type)                               :: dEXO
      Type(MeshTopology_Type)                      :: dMeshTopology
      Type(SectionInt)                             :: dBCFlag 
      
      PetscInt                                     :: iErr, NumDoF, i, j
      PetscInt, Dimension(:), Pointer              :: Flag
      
      Call SectionIntZero(dBCFlag, iErr); CHKERRQ(iErr)
      !!! Side Sets
      !!! To be implemented
      
      !!! Node Sets
      Do i = 1, dMeshTopology%Num_Node_Sets
         Allocate(Flag(1))
         Flag = dEXO%NSProperty( VarFracFilm_NSProp_BCVType )%Value( dMeshTopology%Node_Set(i)%ID )
         Do j = 1, dMeshTopology%Node_Set(i)%Num_Nodes
            Call MeshUpdateAddClosureInt(dMeshTopology%Mesh, dBCFlag, dMeshTopology%Node_Set(i)%Node_ID(j) + dMeshTopology%Num_Elems-1, Flag, iErr); CHKERRQ(iErr)
         End Do
         DeAllocate(Flag)
      End Do
      Call SectionIntComplete(dBCFlag, iErr); CHKERRQ(iErr)
   End Subroutine EXOProperty_InitBCVFlag
   
   Subroutine EXOProperty_InitBCUFlag2D(dEXO, dMeshTopology, dBCFlag)
      Type(EXO_Type)                               :: dEXO
      Type(MeshTopology_Type)                      :: dMeshTopology
      Type(SectionInt)                             :: dBCFlag 
      
      PetscInt                                     :: iErr, NumDoF, i, j
      PetscInt, Dimension(:), Pointer              :: Flag
      
      Call SectionIntZero(dBCFlag, iErr); CHKERRQ(iErr)
      !!! Side Sets
      !!! To be implemented
      
      !!! Node Sets
      Do i = 1, dMeshTopology%Num_Node_Sets
         Allocate(Flag(2))
         Flag(1) = dEXO%NSProperty( VarFracFilm_NSProp_BCUTypeX )%Value( dMeshTopology%Node_Set(i)%ID )
         Flag(2) = dEXO%NSProperty( VarFracFilm_NSProp_BCUTypeY )%Value( dMeshTopology%Node_Set(i)%ID )
         Do j = 1, dMeshTopology%Node_Set(i)%Num_Nodes
            Call MeshUpdateAddClosureInt(dMeshTopology%Mesh, dBCFlag, dMeshTopology%Node_Set(i)%Node_ID(j) + dMeshTopology%Num_Elems-1, Flag, iErr); CHKERRQ(iErr)
         End Do
         DeAllocate(Flag)
      End Do
      Call SectionIntComplete(dBCFlag, iErr); CHKERRQ(iErr)
   End Subroutine EXOProperty_InitBCUFlag2D


   Subroutine MatProp2D_Write(MeshTopology, MatProp, filename)
      Type(MeshTopology_Type)                      :: MeshTopology
      Type(MatProp2D_Type), Dimension(:), Pointer  :: MatProp
      Character(len=*)                             :: filename
      PetscMPIInt                                  :: rank
      PetscInt                                     :: iBlk, Blk_ID
      
      Open(File = filename, Unit = F_OUT, Status = 'Unknown')
      Rewind(F_OUT)
      Write(F_OUT, *) MeshTopology%Num_Elem_Blks_Global
      Do iBlk = 1, Size(MatProp)
         Blk_ID = MeshTopology%Elem_Blk(iBlk)%ID
         Write(F_OUT,120) Blk_ID, MatProp(iBlk)%ToughnessT, MatProp(iBlk)%ToughnessD, MatProp(iBlk)%Hookes_Law, MatProp(iBlk)%Therm_Exp, MatProp(iBlk)%K_Interface
      End Do
      Close(F_OUT)
      
120   Format(I6, '      ', 14(ES12.5,' '))   
   End Subroutine MatProp2D_Write
  
 
   Subroutine MatProp2D_Read(MeshTopology, MatProp, filename)
      Type(MeshTopology_Type)                      :: MeshTopology
      Type(MatProp2D_Type), Dimension(:), Pointer  :: MatProp
      Character(len=*)                             :: filename

      PetscInt                                     :: iBlk, iErr
      
      PetscInt                                     :: NumBlks, IdxMin, IdxMax, Idx
      Type(Tens4OS2D)                              :: Hookes_Law
      PetscReal                                    :: ToughnessT
      PetscReal                                    :: ToughnessD
      Type(MatS2D)                                 :: Therm_Exp
      Type(MatS2D)                                 :: K_interface
   
      Open(File = filename, Unit = F_IN, Status = 'Unknown', Action = 'Read')
      Rewind(F_IN)
      Read(F_IN, *) NumBlks
      If (NumBlks /= MeshTopology%Num_Elem_Blks_Global) Then
         SETERRQ(PETSC_ERR_ARG_SIZ, 'MatProp2DRead: non matching blocks numbers', iErr)
      End If
      !!! Reading the file once first to get the right number of blocks
      IdxMin =  100000
      IdxMax = -100000
      Do iBlk = 1, NumBlks
         Read(F_IN, *) Idx
         !!! Check that this will work!
         IdxMin = Min(IdxMin, Idx)
         IdxMax = Max(IdxMax, Idx)
      End Do
      Allocate(MatProp(IdxMin:IdxMax))
      Rewind(F_IN)
      Read(F_IN, *) Idx
      Do iBlk = 1, NumBlks
         Read(F_IN, *) Idx, ToughnessT, ToughnessD, Hookes_Law, Therm_exp, K_interface
         MatProp(Idx)%ToughnessT   = ToughnessT
         MatProp(Idx)%ToughnessD   = ToughnessD
         MatProp(Idx)%Hookes_Law   = Hookes_Law
         MatProp(Idx)%Therm_Exp    = Therm_Exp
         MatProp(Idx)%K_interface  = K_interface
      End Do
      Close(F_IN)
      Return
   End Subroutine MatProp2D_Read
   


   Subroutine VarFracFilmSchemeParam_View(dSchemeParam, viewer)
      Type(VarFracFilmSchemeParam_Type)            :: dSchemeParam
      Type(PetscViewer)                            :: viewer
      PetscInt                                     :: iErr
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
      
      Write(IOBuffer, "(I1,T32, 'IrrevType')")           dSchemeParam%IrrevType 
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(ES12.5,T32, 'IrrevTol')")        dSchemeParam%IrrevTol
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(L1,T32, 'DoBT')")                dSchemeParam%DoBT 
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(ES12.5,T32, 'BTTol')")           dSchemeParam%BTTol
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(I5,T32, 'BTInt')")               dSchemeParam%BTInt 
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(I1,T32, 'InitV')")               dSchemeParam%InitV 
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(I1,T32, 'InitPhi')")               dSchemeParam%InitPhi 
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(I5,T32, 'NbCracks')")            dSchemeParam%NbCracks
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(ES12.5,T32, 'MaxCrackLength')")  dSchemeParam%MaxCrackLength
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(I5,T32, 'AltMinMaxIter')")       dSchemeParam%AltMinMaxIter
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(ES12.5,T32, 'AltMinTol')")        dSchemeParam%AltMinTol
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(ES12.5,T32, 'KSPUrTol')")        dSchemeParam%KSPUrTol
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(I5,T32, 'AltMinSaveInt')")       dSchemeParam%AltMinSaveInt
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(ES12.5,T32, 'KSPVrTol')")        dSchemeParam%KSPVrTol
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(ES12.5,T32, 'Epsilon')")         dSchemeParam%Epsilon
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(ES12.5,T32, 'KEpsilonV')")        dSchemeParam%KEpsilonV
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(ES12.5,T32, 'KEpsilonPhi')")        dSchemeParam%KEpsilonPhi
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(I1,T32, 'ATNum')")               dSchemeParam%ATNum
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(I1,T32, 'IntegOrder')")          dSchemeParam%IntegOrder
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(L1,T32, 'SaveStress')")          dSchemeParam%SaveStress
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "(L1,T32, 'SaveStrain')")          dSchemeParam%SaveStrain
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
   End Subroutine VarFracFilmSchemeParam_View

   Subroutine VarFracFilmSchemeParam_Load(dSchemeParam, filename)
      Type(VarFracFilmSchemeParam_Type)                   :: dSchemeParam
      Character(len=*)                             :: filename
      
      Open(File = filename, status='old', Unit = F_IN)
      Rewind(F_IN)
      Read(F_IN, *) dSchemeParam%IrrevType
      Read(F_IN, *) dSchemeParam%IrrevTol
      Read(F_IN, *) dSchemeParam%DoBT 
      Read(F_IN, *) dSchemeParam%BTTol
      Read(F_IN, *) dSchemeParam%BTInt 
      Read(F_IN, *) dSchemeParam%InitV 
      Read(F_IN, *) dSchemeParam%InitPhi 
      Read(F_IN, *) dSchemeParam%NbCracks
      Read(F_IN, *) dSchemeParam%MaxCrackLength
      Read(F_IN, *) dSchemeParam%AltMinMaxIter
      Read(F_IN, *) dSchemeParam%AltMinTol
      Read(F_IN, *) dSchemeParam%AltMinSaveInt
      Read(F_IN, *) dSchemeParam%KSPUrTol
      Read(F_IN, *) dSchemeParam%KSPVrTol
      Read(F_IN, *) dSchemeParam%Epsilon
      Read(F_IN, *) dSchemeParam%KEpsilonV
      Read(F_IN, *) dSchemeParam%KEpsilonPhi
      Read(F_IN, *) dSchemeParam%ATNum
      Read(F_IN, *) dSchemeParam%IntegOrder
      Read(F_IN, *) dSchemeParam%SaveStress
      Read(F_IN, *) dSchemeParam%SaveStrain
      Close(F_IN)
   End Subroutine VarFracFilmSchemeParam_Load
   
   Subroutine VarFracFilmSchemeParam_GetFromOptions(dSchemeParam)
      Type(VarFracFilmSchemeParam_Type)            :: dSchemeParam
      PetscInt                                     :: iErr
      PetscTruth                                   :: flag

      dSchemeParam%IrrevType        = Irrev_Eq
      dSchemeParam%IrrevTol         = 1.0D-2
      dSchemeParam%DoBT             = PETSC_FALSE
      dSchemeParam%BTTol            = 0.0D0
      dSchemeParam%BTInt            = 0
      dSchemeParam%InitV            = Init_V_PREV
      dSchemeParam%InitPhi          = Init_Phi_PREV
      dSchemeParam%nbCracks         = 0
      dSchemeParam%MaxCrackLength   = 0.0D0  
      dSchemeParam%AltMinMaxIter    = 1000
      dSchemeParam%AltMinTol        = 1.0D-4
      dSchemeParam%AltMinSaveInt    = 25
      dSchemeParam%KSPUrtol         = 1.0D-6
      dSchemeParam%KSPVrtol         = 1.0D-6
      dSchemeParam%Epsilon          = .1
      dSchemeParam%KEpsilonV        = 1.0E-6
      dSchemeParam%KEpsilonPhi      = 1.0E-6
      dSchemeParam%ATNum            = 2
      dSchemeParam%IntegOrder       = 3
      dSchemeParam%SaveStress       = PETSC_FALSE
      dSchemeParam%SaveStrain       = PETSC_FALSE

      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-irrevtype',      dSchemeParam%IrrevType, flag, iErr); CHKERRQ(iErr) 
      Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,  '-irrevtol',       dSchemeParam%IrrevTol, flag, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetTruth(PETSC_NULL_CHARACTER, '-dobt',           dSchemeParam%DoBT, flag, iErr); CHKERRQ(iErr) 
      Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,  '-bttol',          dSchemeParam%BTTol, flag, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-btint',          dSchemeParam%BTInt, flag, iErr); CHKERRQ(iErr) 
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-initv',          dSchemeParam%InitV, flag, iErr); CHKERRQ(iErr) 
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-initphi',        dSchemeParam%InitPhi, flag, iErr); CHKERRQ(iErr) 
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-nbcracks',       dSchemeParam%NbCracks, flag, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,  '-maxcracklength', dSchemeParam%MaxCrackLength, flag, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-altminmaxiter',  dSchemeParam%AltMinMaxIter, flag, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,  '-altmintol',      dSchemeParam%AltMinTol, flag, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-altminsaveint',  dSchemeParam%AltMinSaveInt, flag, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,  '-kspurtol',       dSchemeParam%KSPUrTol, flag, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,  '-kspvrtol',       dSchemeParam%KSPVrTol, flag, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,  '-epsilon',        dSchemeParam%Epsilon, flag, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,  '-kepsilonv',      dSchemeParam%KEpsilonV, flag, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,  '-kepsilonphi',    dSchemeParam%KEpsilonPhi, flag, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-atnum',          dSchemeParam%ATNum, flag, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-integorder',     dSchemeParam%IntegOrder, flag, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetTruth(PETSC_NULL_CHARACTER, '-savestress',     dSchemeParam%SaveStress, flag, iErr); CHKERRQ(iErr) 
      Call PetscOptionsGetTruth(PETSC_NULL_CHARACTER, '-savestrain',     dSchemeParam%SaveStrain, flag, iErr); CHKERRQ(iErr) 
   End Subroutine VarFracFilmSchemeParam_GetFromOptions
   
   Subroutine VarFracFilmEXOProperty_Init(dEXO, dMeshTopology)
      Type(EXO_Type)                      :: dEXO
      Type(MeshTopology_Type)             :: dMeshTopology
      PetscInt                            :: i, iErr
      PetscInt                            :: NumEB, NumSS, NumNS

      Integer                             :: EXO_MyRank
      PetscReal                           :: rDummy
      Character                           :: cDummy
          

      Call MPI_COMM_RANK(dEXO%Comm, EXO_MyRank, iErr)

      NumEB = dMeshTopology%Num_Elem_Blks_Global
      NumSS = dMeshTopology%Num_Side_Sets_Global
      NumNS = dMeshTopology%Num_Node_Sets_Global
      
      If ( (NumEB == 0) .AND. (NumSS == 0) .AND. (NumSS ==0) ) Then
         Call PetscPrintf(PETSC_COMM_WORLD, '[WARNING]: The EXO file contains no EB, SS or NS is this right?\n', iErr); CHKERRQ(iErr)
         Call PetscPrintf(PETSC_COMM_WORLD, '           Was Write_MeshTopologyGlobal called before VarFracFilmEXOProperty_Init?\n', iErr); CHKERRQ(iErr)
      End If

      dEXO%Num_EBProperties = VarFracFilm_Num_EBProperties
      Allocate(dEXO%EBProperty(dEXO%Num_EBProperties))
      dEXO%EBProperty(VarFracFilm_EBProp_IsBrittle)%Name    = 'Is_Brittle'
      dEXO%EBProperty(VarFracFilm_EBProp_HasSubstrate)%Name = 'HasSubstrate'
      dEXO%EBProperty(VarFracFilm_EBProp_Elem_Type)%Name    = 'Elem_Type'
      Do i = 1, dEXO%Num_EBProperties
         Allocate(dEXO%EBProperty(i)%Value(NumEB))
         dEXO%EBProperty(i)%Value = 0
      End Do
      
      dEXO%Num_SSProperties = VarFracFilm_Num_SSProperties
      Allocate(dEXO%SSProperty(dEXO%Num_SSProperties))
      dEXO%SSProperty(VarFracFilm_SSProp_BCUTypeX)%Name  = 'BCU_Type_X'
      dEXO%SSProperty(VarFracFilm_SSProp_BCUTypeY)%Name  = 'BCU_Type_Y'
      dEXO%SSProperty(VarFracFilm_SSProp_BCVType)%Name   = 'BCV_Type'
      dEXO%SSProperty(VarFracFilm_SSProp_Elem_Type)%Name = 'Elem_Type'
      Do i = 1, dEXO%Num_SSProperties
         Allocate(dEXO%SSProperty(i)%Value(NumSS))
         dEXO%SSProperty(i)%Value = 0
      End Do
      
      dEXO%Num_NSProperties = VarFracFilm_Num_NSProperties
      Allocate(dEXO%NSProperty(dEXO%Num_NSProperties))
      dEXO%NSProperty(VarFracFilm_NSProp_BCUTypeX)%Name  = 'BCU_Type_X'
      dEXO%NSProperty(VarFracFilm_NSProp_BCUTypeY)%Name  = 'BCU_Type_Y'
      dEXO%NSProperty(VarFracFilm_NSProp_BCVType)%Name   = 'BCV_Type'
      Do i = 1, dEXO%Num_NSProperties
         Allocate(dEXO%NSProperty(i)%Value(NumNS))
         dEXO%NSProperty(i)%Value = 0
      End Do
   End Subroutine VarFracFilmEXOProperty_Init   

   Subroutine VarFracFilmEXOVariable_Init(dEXO)
      Type(EXO_Type)                      :: dEXO
      PetscInt                            :: i
      
      dEXO%Num_GlobVariables = VarFracFilm_Num_GlobVar
      Allocate(dEXO%GlobVariable(dEXO%Num_GlobVariables))
      dEXO%GlobVariable(VarFracFilm_GlobVar_SurfaceEnergyT)%Name      = 'Surface energy T'
      dEXO%GlobVariable(VarFracFilm_GlobVar_SurfaceEnergyD)%Name      = 'Surface energy D'
      dEXO%GlobVariable(VarFracFilm_GlobVar_ElasticBulkEnergy)%Name   = 'Bulk elastic energy'
      dEXO%GlobVariable(VarFracFilm_GlobVar_ElasticInterEnergy)%Name  = 'Interface elastic energy '
      dEXO%GlobVariable(VarFracFilm_GlobVar_TotalEnergy)%Name         = 'Total energy'
      dEXO%GlobVariable(VarFracFilm_GlobVar_Load)%Name                = 'Load'
      dEXO%GlobVariable(:)%Offset = (/ (i, i=1,dEXO%Num_GlobVariables) /)
      
      dEXO%Num_CellVariables = VarFracFilm_Num_CellVar
      Allocate(dEXO%CellVariable(dEXO%Num_CellVariables))
      dEXO%CellVariable(VarFracFilm_CellVar_Delamination)%Name  = 'Delamination'
      dEXO%CellVariable(VarFracFilm_CellVar_StrainXX)%Name      = 'Strain XX'
      dEXO%CellVariable(VarFracFilm_CellVar_StrainYY)%Name      = 'Strain YY' 
      dEXO%CellVariable(VarFracFilm_CellVar_StrainXY)%Name      = 'Strain XY'
      dEXO%CellVariable(VarFracFilm_CellVar_StressXX)%Name      = 'Stress XX'
      dEXO%CellVariable(VarFracFilm_CellVar_StressYY)%Name      = 'Stress YY'
      dEXO%CellVariable(VarFracFilm_CellVar_StressXY)%Name      = 'Stress XY'
      dEXO%CellVariable(:)%Offset = (/ (i, i=1,dEXO%Num_CellVariables) /)

      dEXO%Num_VertVariables = VarFracFilm_Num_VertVar
      Allocate(dEXO%VertVariable(dEXO%Num_VertVariables))
      dEXO%VertVariable(VarFracFilm_VertVar_Fracture)%Name      = 'Fracture'
      dEXO%VertVariable(VarFracFilm_VertVar_DisplacementX)%Name = 'Displacement X'   
      dEXO%VertVariable(VarFracFilm_VertVar_DisplacementY)%Name = 'Displacement Y'
      dEXO%VertVariable(VarFracFilm_VertVar_U0X)%Name           = 'U0X'
      dEXO%VertVariable(VarFracFilm_VertVar_U0Y)%Name           = 'U0Y'
      dEXO%VertVariable(VarFracFilm_VertVar_Temperature)%Name   = 'Temperature'
      dEXO%VertVariable(:)%Offset = (/ (i, i=1,dEXO%Num_VertVariables) /)
   End Subroutine VarFracFilmEXOVariable_Init
  
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
End Module m_VarFracFilm_Struct
