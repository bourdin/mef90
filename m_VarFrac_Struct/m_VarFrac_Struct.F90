Module m_VarFrac_Struct

#include "finclude/petscdef.h"

   Use m_MEF90
   Use petsc

   Implicit NONE
   Private

   
   Public :: GenHL_Iso3D_Enu
   Public :: GenHL_Iso2D_EnuPlaneStress
   Public :: GenHL_Iso2D_EnuPlaneStrain
   Public :: GenHL_Ortho2D_LambdaMu
   Public :: GenHL_Laminate_LambdaMu
   
   Public :: VarFracSchemeParam_View
   Public :: VarFracSchemeParam_GetFromOptions
   
   Public :: VarFracEXOProperty_Init
   Public :: VarFracEXOVariable_Init
   
   Public :: MatProp2D_Type, MatProp3D_Type
   Public :: MatProp_Write, MatProp_Read
   
   Public :: EXOProperty_InitBCVFlag
   Public :: EXOProperty_InitBCUFlag2DA
   Public :: EXOProperty_InitBCUFlag
   
   Public :: VarFracSchemeParam_Type
   
   Interface MatProp_Write
      Module Procedure MatProp2D_Write, MatProp3D_Write
   End Interface
  
   Interface MatProp_Read
      Module Procedure MatProp2D_Read, MatProp3D_Read
   End Interface
   
   Interface GenHL_Iso_LambdaMu
      Module Procedure GenHL_Iso2D_LambdaMu, GenHL_Iso3D_LambdaMu
   End Interface
   
   Interface GenHL_Laminate_LambdaMu
      Module Procedure GenHL_Laminate2D_LambdaMu, GenHL_Laminate3D_LambdaMu
   End Interface

   PetscInt, Parameter, Public                     :: VarFrac_BC_Type_NONE = 0
   PetscInt, Parameter, Public                     :: VarFrac_BC_Type_DIRI = 1

   PetscInt, Parameter, Public                     :: VarFrac_Init_V_PREV    = 0
   PetscInt, Parameter, Public                     :: VarFrac_Init_V_RND     = 1
   PetscInt, Parameter, Public                     :: VarFrac_Init_V_SPH     = 2
   PetscInt, Parameter, Public                     :: VarFrac_Init_V_CRACKS  = 3
   
   PetscInt, Parameter, Public                     :: VarFrac_Irrev_NONE = 0
   PetscInt, Parameter, Public                     :: VarFrac_Irrev_Eq   = 1
   PetscInt, Parameter, Public                     :: VarFrac_Irrev_Ineq = 2

   PetscInt, Parameter, Public                     :: VarFrac_Unilateral_NONE  = 0
   PetscInt, Parameter, Public                     :: VarFrac_Unilateral_Full  = 1
   PetscInt, Parameter, Public                     :: VarFrac_Unilateral_Shear = 2

   PetscInt, Parameter, Public                     :: VarFrac_Num_VertVar           = 8
   PetscInt, Parameter, Public                     :: VarFrac_VertVar_Fracture      = 1
   PetscInt, Parameter, Public                     :: VarFrac_VertVar_DisplacementX = 2   
   PetscInt, Parameter, Public                     :: VarFrac_VertVar_DisplacementY = 3
   PetscInt, Parameter, Public                     :: VarFrac_VertVar_DisplacementZ = 4   
   PetscInt, Parameter, Public                     :: VarFrac_VertVar_ForceX        = 5   
   PetscInt, Parameter, Public                     :: VarFrac_VertVar_ForceY        = 6
   PetscInt, Parameter, Public                     :: VarFrac_VertVar_ForceZ        = 7   
   PetscInt, Parameter, Public                     :: VarFrac_VertVar_Temperature   = 8
   
   PetscInt, Parameter, Public                     :: VarFrac_Num_CellVar      = 12
   PetscInt, Parameter, Public                     :: VarFrac_CellVar_StrainXX = 1
   PetscInt, Parameter, Public                     :: VarFrac_CellVar_StrainYY = 2 
   PetscInt, Parameter, Public                     :: VarFrac_CellVar_StrainXY = 3
   PetscInt, Parameter, Public                     :: VarFrac_CellVar_StrainZZ = 4
   PetscInt, Parameter, Public                     :: VarFrac_CellVar_StrainYZ = 5
   PetscInt, Parameter, Public                     :: VarFrac_CellVar_StrainXZ = 6
   PetscInt, Parameter, Public                     :: VarFrac_CellVar_StressXX = 7
   PetscInt, Parameter, Public                     :: VarFrac_CellVar_StressYY = 8
   PetscInt, Parameter, Public                     :: VarFrac_CellVar_StressXY = 9
   PetscInt, Parameter, Public                     :: VarFrac_CellVar_StressZZ = 10
   PetscInt, Parameter, Public                     :: VarFrac_CellVar_StressYZ = 11
   PetscInt, Parameter, Public                     :: VarFrac_CellVar_StressZX = 12

   PetscInt, Parameter, Public                     :: VarFrac_Num_GlobVar           = 6
   PetscInt, Parameter, Public                     :: VarFrac_GlobVar_ElasticEnergy = 1
   PetscInt, Parameter, Public                     :: VarFrac_GlobVar_ExtForcesWork = 2
   PetscInt, Parameter, Public                     :: VarFrac_GlobVar_KineticEnergy = 3 
   PetscInt, Parameter, Public                     :: VarFrac_GlobVar_SurfaceEnergy = 4 
   PetscInt, Parameter, Public                     :: VarFrac_GlobVar_TotalEnergy   = 5
   PetscInt, Parameter, Public                     :: VarFrac_GlobVar_Load          = 6
   
   PetscInt, Parameter, Public                     :: VarFrac_Num_EBProperties  = 3
   PetscInt, Parameter, Public                     :: VarFrac_EBProp_IsBrittle  = 1
   PetscInt, Parameter, Public                     :: VarFrac_EBProp_HasBForce  = 2
   PetscInt, Parameter, Public                     :: VarFrac_EBProp_Elem_Type  = 3
   
   PetscInt, Parameter, Public                     :: VarFrac_Num_SSProperties  = 6
   PetscInt, Parameter, Public                     :: VarFrac_SSProp_BCUTypeX   = 1
   PetscInt, Parameter, Public                     :: VarFrac_SSProp_BCUTypeY   = 2
   PetscInt, Parameter, Public                     :: VarFrac_SSProp_BCUTypeZ   = 3
   PetscInt, Parameter, Public                     :: VarFrac_SSProp_BCVType    = 4
   PetscInt, Parameter, Public                     :: VarFrac_SSProp_HasSForce  = 5
   PetscInt, Parameter, Public                     :: VarFrac_SSProp_Elem_Type  = 6

   PetscInt, Parameter, Public                     :: VarFrac_Num_NSProperties  = 5
   PetscInt, Parameter, Public                     :: VarFrac_NSProp_BCUTypeX   = 1
   PetscInt, Parameter, Public                     :: VarFrac_NSProp_BCUTypeY   = 2
   PetscInt, Parameter, Public                     :: VarFrac_NSProp_BCUTypeZ   = 3
   PetscInt, Parameter, Public                     :: VarFrac_NSProp_BCVType    = 4
   PetscInt, Parameter, Public                     :: VarFrac_NSProp_HasPForce  = 5
   
   Type MatProp2D_Type
      PetscReal                                    :: Toughness
      Type(Tens4OS2D)                              :: Hookes_Law
      Type(MatS2D)                                 :: Therm_Exp      
   End Type MatProp2D_Type
   
   Type MatProp3D_Type
      PetscReal                                    :: Toughness
      Type(Tens4OS3D)                              :: Hookes_Law
      Type(MatS3D)                                 :: Therm_Exp      
   End Type MatProp3D_Type
   
   Type VarFracSchemeParam_Type
      PetscInt                                     :: IrrevType
      PetscReal                                    :: IrrevTol
      
      PetscTruth                                   :: DoBT
      PetscReal                                    :: BTTol
      PetscInt                                     :: BTInt
      PetscInt                                     :: BTScope
      
      PetscInt                                     :: Unilateral
      
      PetscInt                                     :: InitV
      PetscInt                                     :: nbCracks
      PetscReal                                    :: MaxCrackLength     
      
      PetscInt                                     :: AltMinMaxIter
      PetscReal                                    :: AltMinTol
      PetscInt                                     :: AltMinSaveInt

      PetscReal                                    :: Epsilon
      PetscReal                                    :: KEpsilon

      PetscInt                                     :: ATNum
      PetscReal                                    :: ATCv ! The c_v constant in Braides-1998 p.48

      PetscInt                                     :: IntegOrder
      
      PetscTruth                                   :: SaveStress
      PetscTruth                                   :: SaveStrain
      
      PetscTruth                                   :: U_UseTao
      PetscTruth                                   :: V_UseTao
   End Type VarFracSchemeParam_Type
   
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
         Flag = dEXO%NSProperty( VarFrac_NSProp_BCVType )%Value( dMeshTopology%Node_Set(i)%ID )
         Do j = 1, dMeshTopology%Node_Set(i)%Num_Nodes
            Call SectionIntUpdate(dBCFlag, dMeshTopology%Node_Set(i)%Node_ID(j) + dMeshTopology%Num_Elems-1, Flag, ADD_VALUES, iErr); CHKERRQ(iErr)
         End Do
         DeAllocate(Flag)
      End Do
      Call SectionIntComplete(dBCFlag, iErr); CHKERRQ(iErr)
   End Subroutine EXOProperty_InitBCVFlag
   
   Subroutine EXOProperty_InitBCUFlag2DA(dBCFlag, dEXO, dMeshTopology)
      Type(EXO_Type)                               :: dEXO
      Type(MeshTopology_Type)                      :: dMeshTopology
      Type(Flag)                                   :: dBCFlag 
      
      PetscInt                                     :: iErr, NumDoF, i, j
      PetscInt, Dimension(:), Pointer              :: Flag
      
      Call SectionIntZero(dBCFlag, iErr); CHKERRQ(iErr)
      !!! Side Sets
      !!! To be implemented
      
      !!! Node Sets
      Do i = 1, dMeshTopology%Num_Node_Sets
         Allocate(Flag(1))
         Flag(1) = dEXO%NSProperty( VarFrac_NSProp_BCUTypeZ )%Value( dMeshTopology%Node_Set(i)%ID )
         Do j = 1, dMeshTopology%Node_Set(i)%Num_Nodes
            Call SectionIntUpdate(dBCFlag%Sec, dMeshTopology%Node_Set(i)%Node_ID(j) + dMeshTopology%Num_Elems-1, Flag, ADD_VALUES, iErr); CHKERRQ(iErr)
         End Do
         DeAllocate(Flag)
      End Do
   End Subroutine EXOProperty_InitBCUFlag2DA
   
   !!! Rewrite as FlagInitNSProperty(Flag, Property, MeshTopology)
   Subroutine EXOProperty_InitBCUFlag(dBCFlag, dEXO, dMeshTopology)
      Type(EXO_Type)                               :: dEXO
      Type(MeshTopology_Type)                      :: dMeshTopology
      Type(Flag)                                   :: dBCFlag 
      
      PetscInt                                     :: iErr, NumDoF, i, j, k
      PetscInt, Dimension(:), Pointer              :: FlagX, FlagY, FlagZ
      
      Call SectionIntZero(dBCFlag, iErr); CHKERRQ(iErr)
      !!! Side Sets
      !!! To be implemented
      
      !!! Node Sets
      Do i = 1, dMeshTopology%Num_Node_Sets
         Allocate(FlagX(1))
         FlagX = dEXO%NSProperty( VarFrac_NSProp_BCUTypeX )%Value( dMeshTopology%Node_Set(i)%ID )
         Allocate(FlagY(1))
         FlagY = dEXO%NSProperty( VarFrac_NSProp_BCUTypeY )%Value( dMeshTopology%Node_Set(i)%ID )
         If (dMeshTopology%num_dim == 3) Then
            Allocate(FlagZ(1))
            FlagZ = dEXO%NSProperty( VarFrac_NSProp_BCUTypeZ )%Value( dMeshTopology%Node_Set(i)%ID )
         End If
         Do j = 1, dMeshTopology%Node_Set(i)%Num_Nodes
            Call SectionIntUpdate(dBCFlag%Component_Sec(1), dMeshTopology%Node_Set(i)%Node_ID(j) + dMeshTopology%Num_Elems-1, FlagX, ADD_VALUES, iErr); CHKERRQ(iErr)
            Call SectionIntUpdate(dBCFlag%Component_Sec(2), dMeshTopology%Node_Set(i)%Node_ID(j) + dMeshTopology%Num_Elems-1, FlagY, ADD_VALUES, iErr); CHKERRQ(iErr)
            If (dMeshTopology%num_dim == 3) Then
               Call SectionIntUpdate(dBCFlag%Component_Sec(3), dMeshTopology%Node_Set(i)%Node_ID(j) + dMeshTopology%Num_Elems-1, FlagZ, ADD_VALUES, iErr); CHKERRQ(iErr)
            End If
         End Do
         DeAllocate(FlagX)
         DeAllocate(FlagX)
         If (dMeshTopology%num_dim == 3) Then
            DeAllocate(FlagZ)
         End If
      End Do
   End Subroutine EXOProperty_InitBCUFlag
   
   Subroutine MatProp2D_Write(MeshTopology, MatProp, filename)
      Type(MeshTopology_Type)                      :: MeshTopology
      Type(MatProp2D_Type), Dimension(:), Pointer  :: MatProp
      Character(len=*)                             :: filename
      PetscMPIInt                                  :: rank
      PetscInt                                     :: iBlk
      
      Open(File = filename, Unit = F_OUT, Status = 'Unknown')
      Rewind(F_OUT)
      Write(F_OUT, *) MeshTopology%Num_Elem_Blks_Global
      Do iBlk = 1, Size(MatProp)
         Write(F_OUT,120) iBlk, MatProp(iBlk)%Toughness, MatProp(iBlk)%Hookes_Law, MatProp(iBlk)%Therm_Exp
      End Do
      Close(F_OUT)
      
110   Format(I6,'      Toughness    A1111        A1112        A1122        A1212        A1222        A2222        Alpha')
120   Format(I6, '      ', 10(ES12.5,' '))   
   End Subroutine MatProp2D_Write
 
   Subroutine MatProp3D_Write(MeshTopology, MatProp, filename)
      Type(MeshTopology_Type)                      :: MeshTopology
      Type(MatProp3D_Type), Dimension(:), Pointer  :: MatProp
      Character(len=*)                             :: filename

      PetscMPIInt                                  :: rank
      PetscInt                                     :: iBlk
      
      Open(File = filename, Unit = F_OUT, Status = 'Unknown')
      Rewind(F_OUT)
      Write(F_OUT, *) MeshTopology%Num_Elem_Blks_Global
      Do iBlk = 1, Size(MatProp)
         Write(F_OUT,120) iBlk, MatProp(iBlk)%Toughness, MatProp(iBlk)%Hookes_Law, MatProp(iBlk)%Therm_Exp
      End Do
      Close(F_OUT)
      
110   Format(I6,'      Toughness    A_1111       A_1112       A_1113       A_1122       A_1123       A_1133       A_1212       A_1213       A_1222       A_1223       A_12133      A_1313       A_1322       A_1323       A_1333       A_2222       A_2223       A_2233       A_2323       A_2333       A_3333       Alpha')
120   Format(I6, '      ', 28(ES12.5,' '))
   End Subroutine MatProp3D_Write
 
 
   Subroutine MatProp2D_Read(MeshTopology, MatProp, filename)
      Type(MeshTopology_Type)                      :: MeshTopology
      Type(MatProp2D_Type), Dimension(:), Pointer  :: MatProp
      Character(len=*)                             :: filename

      PetscInt                                     :: iBlk, iErr
      
      PetscInt                                     :: NumBlks, IdxMin, IdxMax, Idx
      Type(Tens4OS2D)                              :: Hookes_Law
      PetscReal                                    :: Toughness
      Type(MatS2D)                                 :: Therm_Exp
   
      Open(File = filename, Unit = F_IN, Status = 'Unknown', Action = 'Read')
      Rewind(F_IN)
      Read(F_IN, *) NumBlks
      If (NumBlks /= MeshTopology%Num_Elem_Blks_Global) Then
         SETERRQ(PETSC_ERR_ARG_SIZ, 'MatProp2DRead: non matching blocks numbers', iErr)
      End If
      !!! Reading the file once first to get the right number of blocks
      IdxMin =  100000000
      IdxMax = -100000000
      Do iBlk = 1, NumBlks
         Read(F_IN, *) Idx
         IdxMin = Min(IdxMin, Idx)
         IdxMax = Max(IdxMax, Idx)
      End Do
      Allocate(MatProp(IdxMin:IdxMax))
      Rewind(F_IN)
      Read(F_IN, *) Idx
      Do iBlk = 1, NumBlks
         Read(F_IN, *) Idx, Toughness, Hookes_Law, Therm_exp

         MatProp(Idx)%Toughness  = Toughness
         MatProp(Idx)%Hookes_Law = Hookes_Law
         MatProp(Idx)%Therm_Exp  = Therm_Exp
      End Do
      Close(F_IN)
      Return
!120   Format(I6, '      ', 10(ES12.5,' '))   
!120   Format(*)
   End Subroutine MatProp2D_Read
   
   Subroutine MatProp3D_Read(MeshTopology, MatProp, filename)
      Type(MeshTopology_Type)                      :: MeshTopology
      Type(MatProp3D_Type), Dimension(:), Pointer  :: MatProp
      Character(len=*)                             :: filename

      PetscInt                                     :: iBlk, Blk_Id, iErr
      
      PetscInt                                     :: NumBlks, IdxMin, IdxMax, Idx
      Type(Tens4OS3D)                              :: Hookes_Law
      PetscReal                                    :: Toughness
      Type(MatS3D)                                 :: Therm_Exp
   
      Open(File = filename, Unit = F_IN, Status = 'Old', Action = 'Read')
      Rewind(F_IN)
      Read(F_IN, *) NumBlks
      If (NumBlks /= MeshTopology%Num_Elem_Blks_Global) Then
         SETERRQ(PETSC_ERR_ARG_SIZ, 'MatProp3DRead: non matching blocks numbers', iErr)
      End If
      !!! Reading the file once first to get the right number of blocks
      IdxMin =  100000000
      IdxMax = -100000000
      Do iBlk = 1, NumBlks
         Read(F_IN, *) Idx
         IdxMin = Min(IdxMin, Idx)
         IdxMax = Max(IdxMax, Idx)
      End Do
      Allocate(MatProp(IdxMin:IdxMax))
      Rewind(F_IN)
      Read(F_IN, *) Idx
      Do iBlk = 1, NumBlks
         Read(F_IN, *) Idx, Toughness, Hookes_Law, Therm_exp
         MatProp(Idx)%Toughness  = Toughness
         MatProp(Idx)%Hookes_Law = Hookes_Law
         MatProp(Idx)%Therm_Exp  = Therm_Exp
      End Do
      Close(F_IN)
      
220   Format(I6, 28(ES12.5,' '))
   End Subroutine MatProp3D_Read

   Subroutine VarFracSchemeParam_View(dSchemeParam, viewer)
      Type(VarFracSchemeParam_Type)                :: dSchemeParam
      Type(PetscViewer)                            :: viewer
      PetscInt                                     :: iErr
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
      
      Write(IOBuffer, "('-irrev ', I1, A)")               dSchemeParam%IrrevType, '\n'
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "('-irrevtol ', ES12.5, A)")        dSchemeParam%IrrevTol, '\n'
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "('-bt ', L1, A)")                  dSchemeParam%DoBT, '\n' 
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "('-bttol ', ES12.5, A)")           dSchemeParam%BTTol, '\n'
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "('-btint ', I5, A)")               dSchemeParam%BTInt, '\n' 
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "('-btscope ', I5, A)")             dSchemeParam%BTScope, '\n' 
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "('-unilateral ', I5, A)")          dSchemeParam%Unilateral, '\n' 
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "('-initv ', I1, A)")               dSchemeParam%InitV, '\n' 
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "('-nbcracks ', I5, A)")            dSchemeParam%NbCracks, '\n'
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "('-maxcracklength ', ES12.5, A)")  dSchemeParam%MaxCrackLength, '\n'
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "('-altminmaxiter ', I5, A)")       dSchemeParam%AltMinMaxIter, '\n'
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "('-altmintol', ES12.5, A)")        dSchemeParam%AltMinTol, '\n'
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "('-altminsaveint ', I5, A)")       dSchemeParam%AltMinSaveInt, '\n'
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "('-epsilon ', ES12.5, A)")         dSchemeParam%Epsilon, '\n'
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "('-kepsilon ', ES12.5, A)")        dSchemeParam%KEpsilon, '\n'
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "('-atnum ', I1, A)")               dSchemeParam%ATNum, '\n'
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "('-integorder', I1, A)")          dSchemeParam%IntegOrder, '\n'
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "('-savestress ', L1, A)")          dSchemeParam%SaveStress, '\n'
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "('-savestrain ', L1, A)")          dSchemeParam%SaveStrain, '\n'
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "('-u_tao ', L1, A)")               dSchemeParam%U_UseTao, '\n'
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, "('-v_tao ', L1, A)")               dSchemeParam%V_UseTao, '\n'
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)
   End Subroutine VarFracSchemeParam_View

   Subroutine VarFracSchemeParam_GetFromOptions(dSchemeParam)
      Type(VarFracSchemeParam_Type)                :: dSchemeParam
      PetscInt                                     :: iErr
      PetscTruth                                   :: flag

      dSchemeParam%IrrevType        = VarFrac_Irrev_Ineq
      dSchemeParam%IrrevTol         = 1.0D-2
      dSchemeParam%DoBT             = PETSC_FALSE
      dSchemeParam%BTTol            = 1.0D-2
      dSchemeParam%BTInt            = 10
      dSchemeParam%BTScope          = 10000
      dSchemeParam%Unilateral       = 0
      dSchemeParam%InitV            = VarFrac_Init_V_PREV
      dSchemeParam%nbCracks         = 0
      dSchemeParam%MaxCrackLength   = 0.0D0  
      dSchemeParam%AltMinMaxIter    = 1000
      dSchemeParam%AltMinTol        = 1.0D-4
      dSchemeParam%AltMinSaveInt    = 25

      dSchemeParam%Epsilon          = .1
      dSchemeParam%KEpsilon         = 1.0E-6
      dSchemeParam%ATNum            = 1
      dSchemeParam%IntegOrder       = 2
      dSchemeParam%SaveStress       = PETSC_FALSE
      dSchemeParam%SaveStrain       = PETSC_FALSE
      dSchemeParam%U_UseTao         = PETSC_FALSE
      dSchemeParam%V_UseTao         = PETSC_TRUE

      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-irrev',          dSchemeParam%IrrevType, flag, iErr); CHKERRQ(iErr) 
      Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,  '-irrevtol',       dSchemeParam%IrrevTol, flag, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetTruth(PETSC_NULL_CHARACTER, '-bt',             dSchemeParam%DoBT, flag, iErr); CHKERRQ(iErr) 
      Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,  '-bttol',          dSchemeParam%BTTol, flag, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-btint',          dSchemeParam%BTInt, flag, iErr); CHKERRQ(iErr) 
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-btscope',        dSchemeParam%BTScope, flag, iErr); CHKERRQ(iErr) 
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-unilateral',     dSchemeParam%Unilateral, flag, iErr); CHKERRQ(iErr) 
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-initv',          dSchemeParam%InitV, flag, iErr); CHKERRQ(iErr) 
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-nbcracks',       dSchemeParam%NbCracks, flag, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,  '-maxcracklength', dSchemeParam%MaxCrackLength, flag, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-altminmaxiter',  dSchemeParam%AltMinMaxIter, flag, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,  '-altmintol',      dSchemeParam%AltMinTol, flag, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-altminsaveint',  dSchemeParam%AltMinSaveInt, flag, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,  '-epsilon',        dSchemeParam%Epsilon, flag, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,  '-kepsilon',       dSchemeParam%KEpsilon, flag, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-atnum',          dSchemeParam%ATNum, flag, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-integorder',     dSchemeParam%IntegOrder, flag, iErr); CHKERRQ(iErr)
      Call PetscOptionsGetTruth(PETSC_NULL_CHARACTER, '-savestress',     dSchemeParam%SaveStress, flag, iErr); CHKERRQ(iErr) 
      Call PetscOptionsGetTruth(PETSC_NULL_CHARACTER, '-savestrain',     dSchemeParam%SaveStrain, flag, iErr); CHKERRQ(iErr) 
      Call PetscOptionsGetTruth(PETSC_NULL_CHARACTER, '-u_tao',          dSchemeParam%U_UseTao, flag, iErr); CHKERRQ(iErr) 
      Call PetscOptionsGetTruth(PETSC_NULL_CHARACTER, '-v_tao',          dSchemeParam%V_UseTao, flag, iErr); CHKERRQ(iErr) 
      
      Select Case(dSchemeParam%ATNum)
      !!! Braides 2008 p. 48
      !!! \int V(v)/\eps + \eps |\nabla v| Gamma converges to 4c_V H^{N-1){S_u)
      !!! with Cv = \int_0^1 \sqrt(V(s))ds

      Case(1)
         dSchemeParam%ATCv = 2.0_Kr / 3.0_Kr
      Case(2)
         dSchemeParam%ATCV = 0.5_Kr
      Case Default
         SETERRQ(PETSC_ERR_SUP, 'Only AT1 and AT2 are implemented\n', iErr)
      End Select
   End Subroutine VarFracSchemeParam_GetFromOptions
   
   Subroutine VarFracEXOProperty_Init(dEXO, dMeshTopology)
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
         Call PetscPrintf(PETSC_COMM_WORLD, '           Was Write_MeshTopologyGlobal called before VarFracEXOProperty_Init?\n', iErr); CHKERRQ(iErr)
      End If

      dEXO%Num_EBProperties = VarFrac_Num_EBProperties
      Allocate(dEXO%EBProperty(dEXO%Num_EBProperties))
      dEXO%EBProperty(VarFrac_EBProp_IsBrittle)%Name = 'Is_Brittle'
      dEXO%EBProperty(VarFrac_EBProp_HasBForce)%Name = 'Has_BForce'
      dEXO%EBProperty(VarFrac_EBProp_Elem_Type)%Name = 'Elem_Type'
      Do i = 1, dEXO%Num_EBProperties
         Allocate(dEXO%EBProperty(i)%Value(NumEB))
         dEXO%EBProperty(i)%Value = 0
      End Do
      
      dEXO%Num_SSProperties = VarFrac_Num_SSProperties
      Allocate(dEXO%SSProperty(dEXO%Num_SSProperties))
      dEXO%SSProperty(VarFrac_SSProp_BCUTypeX)%Name  = 'BCU_Type_X'
      dEXO%SSProperty(VarFrac_SSProp_BCUTypeY)%Name  = 'BCU_Type_Y'
      dEXO%SSProperty(VarFrac_SSProp_BCUTypeZ)%Name  = 'BCU_Type_Z'
      dEXO%SSProperty(VarFrac_SSProp_BCVType)%Name   = 'BCV_Type'
      dEXO%SSProperty(VarFrac_SSProp_HasSForce)%Name = 'Has_SForce'
      dEXO%SSProperty(VarFrac_SSProp_Elem_Type)%Name = 'Elem_Type'
      Do i = 1, dEXO%Num_SSProperties
         Allocate(dEXO%SSProperty(i)%Value(NumSS))
         dEXO%SSProperty(i)%Value = 0
      End Do
      
      dEXO%Num_NSProperties = VarFrac_Num_NSProperties
      Allocate(dEXO%NSProperty(dEXO%Num_NSProperties))
      dEXO%NSProperty(VarFrac_NSProp_BCUTypeX)%Name  = 'BCU_Type_X'
      dEXO%NSProperty(VarFrac_NSProp_BCUTypeY)%Name  = 'BCU_Type_Y'
      dEXO%NSProperty(VarFrac_NSProp_BCUTypeZ)%Name  = 'BCU_Type_Z'
      dEXO%NSProperty(VarFrac_NSProp_BCVType)%Name   = 'BCV_Type'
      dEXO%NSProperty(VarFrac_NSProp_HasPForce)%Name = 'Has_PForce'
      Do i = 1, dEXO%Num_NSProperties
         Allocate(dEXO%NSProperty(i)%Value(NumNS))
         dEXO%NSProperty(i)%Value = 0
      End Do
   End Subroutine VarFracEXOProperty_Init   

   Subroutine VarFracEXOVariable_Init(dEXO, dSkipElementVariables)
      Type(EXO_Type)                      :: dEXO
      PetscInt                            :: i
      PetscTruth, optional                :: dSkipElementVariables
      
      dEXO%Num_GlobVariables = VarFrac_Num_GlobVar
      Allocate(dEXO%GlobVariable(dEXO%Num_GlobVariables))
      dEXO%GlobVariable(VarFrac_GlobVar_SurfaceEnergy)%Name = 'Surface energy'
      dEXO%GlobVariable(VarFrac_GlobVar_ElasticEnergy)%Name = 'Elastic energy'
      dEXO%GlobVariable(VarFrac_GlobVar_KineticEnergy)%Name = 'Kinetic energy'
      dEXO%GlobVariable(VarFrac_GlobVar_ExtForcesWork)%Name = 'External Forces Work'
      dEXO%GlobVariable(VarFrac_GlobVar_TotalEnergy)%Name   = 'Total energy'
      dEXO%GlobVariable(VarFrac_GlobVar_Load)%Name          = 'Load'
      dEXO%GlobVariable(:)%Offset = (/ (i, i=1,dEXO%Num_GlobVariables) /)
      
      If ( Present(dSkipElementVariables) .AND. (.NOT. dSkipElementVariables)) Then
         dEXO%Num_CellVariables = VarFrac_Num_CellVar
         Allocate(dEXO%CellVariable(dEXO%Num_CellVariables))
         dEXO%CellVariable(VarFrac_CellVar_StrainXX)%Name = 'Strain XX'
         dEXO%CellVariable(VarFrac_CellVar_StrainYY)%Name = 'Strain YY' 
         dEXO%CellVariable(VarFrac_CellVar_StrainZZ)%Name = 'Strain ZZ'
         dEXO%CellVariable(VarFrac_CellVar_StrainXY)%Name = 'Strain XY'
         dEXO%CellVariable(VarFrac_CellVar_StrainYZ)%Name = 'Strain YZ'
         dEXO%CellVariable(VarFrac_CellVar_StrainXZ)%Name = 'Strain XZ'
         dEXO%CellVariable(VarFrac_CellVar_StressXX)%Name = 'Stress XX'
         dEXO%CellVariable(VarFrac_CellVar_StressYY)%Name = 'Stress YY'
         dEXO%CellVariable(VarFrac_CellVar_StressZZ)%Name = 'Stress ZZ'
         dEXO%CellVariable(VarFrac_CellVar_StressXY)%Name = 'Stress XY'
         dEXO%CellVariable(VarFrac_CellVar_StressYZ)%Name = 'Stress YZ'
         dEXO%CellVariable(VarFrac_CellVar_StressZX)%Name = 'Stress ZX'
         dEXO%CellVariable(:)%Offset = (/ (i, i=1,dEXO%Num_CellVariables) /)
      Else
         dEXO%Num_CellVariables = 0
      End If
         
      dEXO%Num_VertVariables = VarFrac_Num_VertVar
      Allocate(dEXO%VertVariable(dEXO%Num_VertVariables))
      dEXO%VertVariable(VarFrac_VertVar_Fracture)%Name      = 'Fracture'
      dEXO%VertVariable(VarFrac_VertVar_DisplacementX)%Name = 'Displacement X'   
      dEXO%VertVariable(VarFrac_VertVar_DisplacementY)%Name = 'Displacement Y'
      dEXO%VertVariable(VarFrac_VertVar_DisplacementZ)%Name = 'Displacement Z'   
      dEXO%VertVariable(VarFrac_VertVar_ForceX)%Name        = 'Force X'   
      dEXO%VertVariable(VarFrac_VertVar_ForceY)%Name        = 'Force Y'
      dEXO%VertVariable(VarFrac_VertVar_ForceZ)%Name        = 'Force Z'   
      dEXO%VertVariable(VarFrac_VertVar_Temperature)%Name   = 'Temperature'
      dEXO%VertVariable(:)%Offset = (/ (i, i=1,dEXO%Num_VertVariables) /)
   End Subroutine VarFracEXOVariable_Init
  
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
   
   PetscReal Function StrainProjectionComponent2D_LambdaMu(xi, k, lambda, mu)
      !!! Compute the projection operator on strain fields as in Allaire 1997 Corollary 1.14.13
      !!! f_B(k) \xi \dot \xi = \frac{1}{\mu} ( |\xi k|^2 - (\xi k \dot k)^2) + \frac{1}{\lamba + 2\mu}(\xi k \dot k)^2
      Type(MatS2D), Intent(IN)            :: xi
      Type(Vect2D), Intent(IN)            :: k
      PetscReal, Intent(IN)               :: lambda, mu
      
      Type(Tens4OS2D)                     :: TmpProj
      Type(Vect2D)                        :: xik, knorm
      PetscReal                           :: xikk, xik2
      
      knorm = k / Norm(k)
      
      xik  = xi * knorm
      xikk = xik .DotP. knorm
      xik2 = xik%X**2 + xik%Y**2
      StrainProjectionComponent2D_LambdaMu = (xik2 - xikk**2) / mu + xikk**2 / (lambda + 2.0_Kr * mu)
   End Function StrainProjectionComponent2D_LambdaMu
   
   PetscReal Function StrainProjectionComponent3D_LambdaMu(xi, k, lambda, mu)
      !!! Compute the projection operator on strain fields as in Allaire 1997 Corollary 1.14.13
      !!! f_B(k) \xi \dot \xi = \frac{1}{\mu} ( |\xi k|^2 - (\xi k \dot k)^2) + \frac{1}{\lamba + 2\mu}(\xi k \dot k)^2
      Type(MatS3D), Intent(IN)            :: xi
      Type(Vect3D), Intent(IN)            :: k
      PetscReal, Intent(IN)               :: lambda, mu
      
      Type(Tens4OS3D)                     :: TmpProj
      Type(Vect3D)                        :: xik, knorm
      PetscReal                           :: xikk, xik2
      
      knorm = k / Norm(k)
      
      xik  = xi * knorm
      xikk = xik .DotP. knorm
      xik2 = xik%X**2 + xik%Y**2 * xik%Z**2
      StrainProjectionComponent3D_LambdaMu = (xik2 - xikk**2) / mu + xikk**2 / (lambda + 2.0_Kr * mu)
   End Function StrainProjectionComponent3D_LambdaMu

   Type(Tens4OS2D) Function StrainProjection2D_LambdaMu(k, lambda, mu)
      Type(Vect2D), Intent(IN)            :: k
      PetscReal, Intent(IN)               :: lambda, mu
      
      Type(MatS2D)                        :: xi
      Type(Tens4OS2D)                     :: fB
      
      fB = 0.0_Kr
      
      xi = 0.0_Kr
      xi%XX = 1.0_Kr
      fB%XXXX = StrainProjectionComponent2D_LambdaMu(xi, k, lambda, mu)
      
      xi = 0.0_Kr
      xi%YY = 1.0_Kr
      fB%YYYY = StrainProjectionComponent2D_LambdaMu(xi, k, lambda, mu)

      xi = 0.0_Kr
      xi%XY = 0.5_Kr
      fB%XYXY = StrainProjectionComponent2D_LambdaMu(xi, k, lambda, mu)

      xi = 0.0_Kr
      xi%XX = 1.0_Kr; xi%YY = 1.0_Kr
      fB%XXYY = (StrainProjectionComponent2D_LambdaMu(xi, k, lambda, mu) - fB%XXXX - fB%YYYY) * 0.5_Kr

      xi = 0.0_Kr
      xi%XX = 1.0_Kr; xi%XY = 0.5_Kr
      fB%XXXY = (StrainProjectionComponent2D_LambdaMu(xi, k, lambda, mu) - fB%XXXX - fB%XYXY) * 0.5_Kr

      xi = 0.0_Kr
      xi%YY = 1.0_Kr; xi%XY = 0.5_Kr
      fB%XYYY = (StrainProjectionComponent2D_LambdaMu(xi, k, lambda, mu) - fB%YYYY - fB%XYXY) * 0.5_Kr
      
      StrainProjection2D_LambdaMu = fB
   End Function StrainProjection2D_LambdaMu
   
   Type(Tens4OS3D) Function StrainProjection3D_LambdaMu(k, lambda, mu)
      Type(Vect3D), Intent(IN)            :: k
      PetscReal, Intent(IN)               :: lambda, mu
      
      Type(MatS3D)                        :: xi
      Type(Tens4OS3D)                     :: fB
      
      
      Write(*,*) 'StrainProjection3D_LambdaMu Not implemented yet'
      STOP
      StrainProjection3D_LambdaMu = 0.0_Kr
   End Function StrainProjection3D_LambdaMu


   Subroutine GenHL_Laminate2D_LambdaMu(A, k, lambda, mu, theta, Astar)
      Type(Tens4OS2D), Intent(IN)         :: A
      Type(Vect2D), Intent(IN)            :: k
      PetscReal, Intent(IN)               :: lambda, mu, theta
      Type(Tens4OS2D), Intent(OUT)        :: Astar
      
      Type(Tens4OS2D)                     :: B
      
      Call GenHL_Iso_LambdaMu(lambda, mu, B)
      
      Astar = Invert(Invert(A-B)  + (1.0_Kr-theta) * StrainProjection2D_LambdaMu(k, lambda, mu)) * theta + B
   End Subroutine GenHL_Laminate2D_LambdaMu

   Subroutine GenHL_Laminate3D_LambdaMu(A, k, lambda, mu, theta, Astar)
      Type(Tens4OS3D), Intent(IN)         :: A
      Type(Vect3D), Intent(IN)            :: k
      PetscReal, Intent(IN)               :: lambda, mu, theta
      Type(Tens4OS3D), Intent(OUT)        :: Astar
   
      Type(Tens4OS3D)                     :: B
   
      Call GenHL_Iso_LambdaMu(lambda, mu, B)
   
      Astar = Invert(Invert(A-B)  + (1.0_Kr-theta) * StrainProjection3D_LambdaMu(k, lambda, mu)) * theta + B
   End Subroutine GenHL_Laminate3D_LambdaMu
      
End Module m_VarFrac_Struct
