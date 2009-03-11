Module m_Rupt_Struct

#include "finclude/petscdef.h"

   Use m_MEF90
   Use petsc

   Implicit NONE
   Private

   
   Public :: GenHL_Iso3D_Enu
   Public :: GenHL_Iso2D_EnuPlaneStress
   Public :: GenHL_Iso2D_EnuPlaneStrain
   Public :: GenHL_Ortho2D_LambdaMu
   
   Public :: SchemeParamView
   Public :: SchemeParamLoad
   Public :: SchemeParamGetFromOptions

   Public :: MatProp2D_Type, MatProp3D_Type
   Public :: EXO_RuptProperties_Type

   
   Interface MatPropWrite
      Module Procedure MatProp2DWrite, MatProp3DWrite
   End Interface
  
   Interface MatPropRead
      Module Procedure MatProp2DRead, MatProp3DRead
   End Interface
   
   Interface GenHL_Iso_LambdaMu
      Module Procedure GenHL_Iso2D_LambdaMu, GenHL_Iso3D_LambdaMu
   End Interface

   Integer, Parameter, Public                      :: BC_Type_NONE = 0
   Integer, Parameter, Public                      :: BC_Type_DIRI = 1
   

!!! Add Side Set properties
!!! Add BC to EB and SS
   Integer, Parameter                              :: Num_Prop_EB  = 3
   Character(len=MXSTLN), Dimension(6), Parameter  :: Prop_Name_EB =          &
       &     (/ 'Is_Brittle    ',                                             &
       &        'Is_Domain     ',                                             &
       &        'Has_BodyForce ',                                             &
       &        'BC_Type_X ',                                                 &
       &        'BC_Type_Y ',                                                 &
       &        'BC_Type_Z ' /) 
   
   Integer, Parameter                              :: Num_Prop_SS  = 3
   Character(len=MXSTLN), Dimension(4), Parameter  :: Prop_Name_SS =          &
       &     (/ 'Has_SurfForce ',                                             &
       &        'BC_Type_X ',                                                 &
       &        'BC_Type_Y ',                                                 &
       &        'BC_Type_Z ' /) 

   Integer, Parameter                              :: Num_Prop_NS  = 3
   Character(len=MXSTLN), Dimension(3), Parameter  :: Prop_Name_NS =          &
       &     (/ 'BC_Type_X ',                                                 &
       &        'BC_Type_Y ',                                                 &
       &        'BC_Type_Z ' /)
   
   Integer, Parameter                              :: Num_Res_G = 5
   Character(len=MXSTLN), Dimension(5), Parameter  :: Res_Name_G =            &
       &     (/ 'Bulk energy   ',                                             &
       &        'Surface energy',                                             &
       &        'Total energy  ',                                             &
       &        'Load          ',                                             &
       &        'Analysis time '/)
   
   Integer, Parameter                              :: Num_Res_N  = 8
   Character(len=MXSTLN), Dimension(8), Parameter  :: Res_Name_N =            &
      &      (/ 'Fracture      ',                                             &
      &         'Displacement X',                                             &
      &         'Displacement Y',                                             &
      &         'Displacement Z',                                             &
      &         'Force X       ',                                             &
      &         'Force Y       ',                                             &
      &         'Force Z       ',                                             &
      &         'Temperature   ' /)
   
   Integer, Parameter                              :: Num_Res_E  = 12
   Character(len=MXSTLN), Dimension(12), Parameter :: Res_Name_E =            &
      &      (/ 'Strain XX',                                                  &
      &         'Strain YY',                                                  &
      &         'Strain ZZ',                                                  &
      &         'Strain XY',                                                  &
      &         'Strain YZ',                                                  &
      &         'Strain XZ',                                                  &
      &         'Stress XX',                                                  &
      &         'Stress YY',                                                  &
      &         'Stress ZZ',                                                  &
      &         'Stress XY',                                                  &
      &         'Stress YZ',                                                  &
      &         'Stress ZX' /)
      
   Type RuptAppCtx
   
   End Type RuptAppCtx

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
   
   Type SchemeParam_Type
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
   End Type SchemeParam_Type
   
   PetscInt, Parameter, Public                     :: Init_V_PREV = 0
   PetscInt, Parameter, Public                     :: Init_V_ONE  = 1
   PetscInt, Parameter, Public                     :: Init_V_RND  = 2
   PetscInt, Parameter, Public                     :: Init_V_SPH  = 3
   
   PetscInt, Parameter, Public                     :: Init_U_PREV = 0
   PetscInt, Parameter, Public                     :: Init_U_ZERO = 1
   
   PetscInt, Parameter, Public                     :: Irrev_NONE = 0
   PetscInt, Parameter, Public                     :: irrev_eq   = 1
   PetscInt, Parameter, Public                     :: Irrev_Ineq = 2
   
   Type EXO_RuptProperties_Type
      !!! Properties stored in the exodus file
      
      ! ELEMENT BLOCK PROPERTIES
      Logical, Dimension(:), Pointer               :: Is_Brittle
      Logical, Dimension(:), Pointer               :: Is_Domain
      Logical, Dimension(:), Pointer               :: Has_BodyForce
      PetscInt, Dimension(:), Pointer              :: EB_BC_Type_X
      PetscInt, Dimension(:), Pointer              :: EB_BC_Type_Y
      PetscInt, Dimension(:), Pointer              :: EB_BC_Type_Z
      
      ! SIDE SETS PROPERTIES
      Logical, Dimension(:), Pointer               :: Has_SurfForce
      PetscInt, Dimension(:), Pointer              :: SS_BC_Type_X
      PetscInt, Dimension(:), Pointer              :: SS_BC_Type_Y
      PetscInt, Dimension(:), Pointer              :: SS_BC_Type_Z
      
      
      ! NODE SETS PROPERTIES 
      PetscInt, Dimension(:), Pointer              :: NS_BC_Type_X
      PetscInt, Dimension(:), Pointer              :: NS_BC_Type_Y
      PetscInt, Dimension(:), Pointer              :: NS_BC_Type_Z
   End Type EXO_RuptProperties_Type

 Contains
   Subroutine MatProp2DWrite(MeshTopology, MatProp, filename)
      Type(MeshTopology_Info)                      :: MeshTopology
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
   End Subroutine MatProp2DWrite
 
   Subroutine MatProp3DWrite(MeshTopology, MatProp, filename)
      Type(MeshTopology_Info)                      :: MeshTopology
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
   End Subroutine MatProp3DWrite
 
 
   Subroutine MatProp2DRead(MeshTopology, MatProp, filename)
      Type(MeshTopology_Info)                      :: MeshTopology
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
   End Subroutine MatProp2DRead
   
   Subroutine MatProp3DRead(MeshTopology, MatProp, filename)
      Type(MeshTopology_Info)                      :: MeshTopology
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
   End Subroutine MatProp3DRead

   Subroutine EXO_RuptFormat(dEXO)
      Type(EXO_Info)                                :: dEXO
      PetscInt                                      :: iErr
      Integer                                       :: exo_ver
      
      dEXO%exoid = EXOPEN(dEXO%filename, EXREAD, exo_cpu_ws, exo_io_ws, exo_ver, iErr)
      
      !!! Write Property names
      Call EXPPN(dEXO%exoid, EXEBLK, Num_Prop_EB, Prop_Name_EB, iErr)
      Call EXPPN(dEXO%exoid, EXSSET, Num_Prop_SS, Prop_Name_SS, iErr)
      Call EXPPN(dEXO%exoid, EXNSET, Num_Prop_NS, Prop_Name_NS, iErr)

      !!! Write Variable Properties (i.e. create variables in the file)
      Call EXPVP (dEXO%exoid, 'g', Num_Res_G, iErr)
      Call EXPVP (dEXO%exoid, 'e', Num_Res_E, iErr)
      Call EXPVP (dEXO%exoid, 'n', Num_Res_N, iErr)

      !!! Write Variable Names
      Call EXPVAN (dEXO%exoid, 'g', Num_Res_G, Res_Name_G, iErr)
      Call EXPVAN (dEXO%exoid, 'e', Num_Res_E, Res_Name_E, iErr)
      Call EXPVAN (dEXO%exoid, 'n', Num_Res_N, Res_Name_N, iErr)

      Call EXCLOS(dEXO%exoid, iErr)
      dEXO%exoid = 0
   End Subroutine EXO_RuptFormat
   
   Subroutine EXO_RuptRead(dEXO, dMeshTopology, dEXO_RuptProperties)
      Type(EXO_Info)                                :: dEXO
      Type(MeshTopology_Info)                       :: dMeshTopology
      Type(EXO_RuptProperties_Type)                 :: dEXO_RuptProperties
      PetscInt                                      :: iErr
      PetscInt                                      :: i, j
      PetscInt                                      :: iRec, NumRec
      PetscInt, Dimension(:), Pointer               :: Tmp_Prop, Ids, GlobalId
      Logical                                       :: Do_IO=.FALSE.
      Integer                                       :: exo_ver
      
      
      If ( ((dEXO%comm == PETSC_COMM_WORLD) .AND. (MEF90_MyRank == 0)) .OR. (dEXO%comm == PETSC_COMM_SELF) ) Then
         Do_IO = .TRUE.
         dEXO%exoid = EXOPEN(dEXO%filename, EXREAD, exo_cpu_ws, exo_io_ws, exo_ver, iErr)
      End If
      
      !!! Element block properties
      ! Get the number of records on cpu 0 of the communicator, and bcast them
      NumRec = dMeshTopology%Num_Elem_Blks
      Call MPI_Bcast(NumRec, 1, MPI_INTEGER, 0, dEXO%comm, iErr)
      Allocate(Tmp_Prop(NumRec))
      Allocate(IDs(NumRec))
      Allocate(GlobalId(dMeshTopology%Num_Elem_Blks))

      If (MEF90_MyRank == 0) Then
         Do iRec = 1, NumRec
            IDs(iRec) = dMeshTopology%Elem_Blk(iRec)%ID
         End Do
      End If
      Call MPI_Bcast(IDs, NumRec, MPI_INTEGER, 0, dEXO%comm, iErr)
      Do i = 1, dMeshTopology%Num_Elem_Blks
         Do j = 1, NumRec
               If (dMeshTopology%Elem_blk(i)%ID == IDs(j)) Then
                  GlobalID(i) = j
                  EXIT
               End If
         End Do
      End Do
      
      !!! Is_Brittle
      Do iRec = 1, dMeshTopology%Num_Elem_Blks
         If (Do_IO) Then
            Call EXGP(dEXO%exoid, EXEBLK, dMeshTopology%Elem_Blk(iRec)%ID, Prop_Name_EB(1), Tmp_Prop(iRec), iErr)
         End If
      End Do
      Call MPI_Bcast(Tmp_Prop, NumRec, MPI_INTEGER, 0, dEXO%comm, iErr)
      Allocate(dEXO_RuptProperties%Is_Brittle(dMeshTopology%Num_Elem_Blks))
      Do iRec = 1, dMeshTopology%Num_Elem_Blks
         dEXO_RuptProperties%Is_Brittle(i) = Tmp_Prop(GlobalId(i))
      End Do
     
      !!! Is_Domain
      Do iRec = 1, dMeshTopology%Num_Elem_Blks
         If (Do_IO) Then
            Call EXGP(dEXO%exoid, EXEBLK, dMeshTopology%Elem_Blk(iRec)%ID, Prop_Name_EB(2), Tmp_Prop(iRec), iErr)
         End If
      End Do
      Call MPI_Bcast(Tmp_Prop, NumRec, MPI_INTEGER, 0, dEXO%comm, iErr)
      Allocate(dEXO_RuptProperties%Is_Domain(dMeshTopology%Num_Elem_Blks))
      Do iRec = 1, dMeshTopology%Num_Elem_Blks
         dEXO_RuptProperties%Is_Domain(i) = Tmp_Prop(GlobalId(i))
      End Do

      !!! Has_BodyForce
      Do iRec = 1, dMeshTopology%Num_Elem_Blks
         If (Do_IO) Then
            Call EXGP(dEXO%exoid, EXEBLK, dMeshTopology%Elem_Blk(iRec)%ID, Prop_Name_EB(3), Tmp_Prop(iRec), iErr)
         End If
      End Do
      Call MPI_Bcast(Tmp_Prop, NumRec, MPI_INTEGER, 0, dEXO%comm, iErr)
      Allocate(dEXO_RuptProperties%Has_BodyForce(dMeshTopology%Num_Elem_Blks))
      Do iRec = 1, dMeshTopology%Num_Elem_Blks
         dEXO_RuptProperties%Has_BodyForce(i) = Tmp_Prop(GlobalId(i))
      End Do

      !!! EB_BC_Type_X
      Do iRec = 1, dMeshTopology%Num_Elem_Blks
         If (Do_IO) Then
            Call EXGP(dEXO%exoid, EXEBLK, dMeshTopology%Elem_Blk(iRec)%ID, Prop_Name_EB(4), Tmp_Prop(iRec), iErr)
         End If
      End Do
      Call MPI_Bcast(Tmp_Prop, NumRec, MPI_INTEGER, 0, dEXO%comm, iErr)
      Allocate(dEXO_RuptProperties%EB_BC_Type_X(dMeshTopology%Num_Elem_Blks))
      Do iRec = 1, dMeshTopology%Num_Elem_Blks
         dEXO_RuptProperties%EB_BC_Type_X(i) = Tmp_Prop(GlobalId(i))
      End Do
      !!! EB_BC_Type_Y
      Do iRec = 1, dMeshTopology%Num_Elem_Blks
         If (Do_IO) Then
            Call EXGP(dEXO%exoid, EXEBLK, dMeshTopology%Elem_Blk(iRec)%ID, Prop_Name_EB(5), Tmp_Prop(iRec), iErr)
         End If
      End Do
      Call MPI_Bcast(Tmp_Prop, NumRec, MPI_INTEGER, 0, dEXO%comm, iErr)
      Allocate(dEXO_RuptProperties%EB_BC_Type_Y(dMeshTopology%Num_Elem_Blks))
      Do iRec = 1, dMeshTopology%Num_Elem_Blks
         dEXO_RuptProperties%EB_BC_Type_Y(i) = Tmp_Prop(GlobalId(i))
      End Do
      !!! EB_BC_Type_Z
      Do iRec = 1, dMeshTopology%Num_Elem_Blks
         If (Do_IO) Then
            Call EXGP(dEXO%exoid, EXEBLK, dMeshTopology%Elem_Blk(iRec)%ID, Prop_Name_EB(6), Tmp_Prop(iRec), iErr)
         End If
      End Do
      Call MPI_Bcast(Tmp_Prop, NumRec, MPI_INTEGER, 0, dEXO%comm, iErr)
      Allocate(dEXO_RuptProperties%EB_BC_Type_Z(dMeshTopology%Num_Elem_Blks))
      Do iRec = 1, dMeshTopology%Num_Elem_Blks
         dEXO_RuptProperties%EB_BC_Type_Z(i) = Tmp_Prop(GlobalId(i))
      End Do
      DeAllocate(Tmp_Prop)
      DeAllocate(IDs)
      DeAllocate(GlobalId)

      !!! Side Sets
      ! To Do
      
      !!! Node Sets
      ! Get the number of records on cpu 0 of the communicator, and bcast them
      NumRec = dMeshTopology%Num_Node_Sets
      Call MPI_Bcast(NumRec, 1, MPI_INTEGER, 0, dEXO%comm, iErr)
      Allocate(Tmp_Prop(NumRec))
      Allocate(IDs(NumRec))
      Allocate(GlobalId(dMeshTopology%Num_Node_Sets))

      If (MEF90_MyRank == 0) Then
         Do iRec = 1, NumRec
            IDs(iRec) = dMeshTopology%Node_Set(iRec)%ID
         End Do
      End If
      Call MPI_Bcast(IDs, NumRec, MPI_INTEGER, 0, dEXO%comm, iErr)
      Do i = 1, dMeshTopology%Num_Node_Sets
         Do j = 1, NumRec
               If (dMeshTopology%Node_Set(i)%ID == IDs(j)) Then
                  GlobalID(i) = j
                  EXIT
               End If
         End Do
      End Do
      
      !!! NS_BC_Type_X
      Do iRec = 1, dMeshTopology%Num_Node_Sets
         If (Do_IO) Then
            Call EXGP(dEXO%exoid, EXNSET, dMeshTopology%Elem_Blk(iRec)%ID, Prop_Name_NS(1), Tmp_Prop(iRec), iErr)
         End If
      End Do
      Call MPI_Bcast(Tmp_Prop, NumRec, MPI_INTEGER, 0, dEXO%comm, iErr)
      Allocate(dEXO_RuptProperties%NS_BC_Type_X(dMeshTopology%Num_Node_Sets))
      Do iRec = 1, dMeshTopology%Num_Node_Sets
         dEXO_RuptProperties%NS_BC_Type_X(i) = Tmp_Prop(GlobalId(i))
      End Do
      !!! NS_BC_Type_Y
      Do iRec = 1, dMeshTopology%Num_Node_Sets
         If (Do_IO) Then
            Call EXGP(dEXO%exoid, EXNSET, dMeshTopology%Elem_Blk(iRec)%ID, Prop_Name_NS(2), Tmp_Prop(iRec), iErr)
         End If
      End Do
      Call MPI_Bcast(Tmp_Prop, NumRec, MPI_INTEGER, 0, dEXO%comm, iErr)
      Allocate(dEXO_RuptProperties%NS_BC_Type_Y(dMeshTopology%Num_Node_Sets))
      Do iRec = 1, dMeshTopology%Num_Node_Sets
         dEXO_RuptProperties%NS_BC_Type_Z(i) = Tmp_Prop(GlobalId(i))
      End Do
      !!! NS_BC_Type_Z
      Do iRec = 1, dMeshTopology%Num_Node_Sets
         If (Do_IO) Then
            Call EXGP(dEXO%exoid, EXNSET, dMeshTopology%Elem_Blk(iRec)%ID, Prop_Name_NS(3), Tmp_Prop(iRec), iErr)
         End If
      End Do
      Call MPI_Bcast(Tmp_Prop, NumRec, MPI_INTEGER, 0, dEXO%comm, iErr)
      Allocate(dEXO_RuptProperties%NS_BC_Type_Z(dMeshTopology%Num_Node_Sets))
      Do iRec = 1, dMeshTopology%Num_Node_Sets
         dEXO_RuptProperties%NS_BC_Type_Z(i) = Tmp_Prop(GlobalId(i))
      End Do
      DeAllocate(Tmp_Prop)
      DeAllocate(IDs)
      DeAllocate(GlobalId)
      
      
      If (Do_IO) Then
         Call EXCLOS(dEXO%exoid, iErr)
         dEXO%exoid = 0
      End If
   End Subroutine EXO_RuptRead

   Subroutine EXO_RuptWrite(dEXO, dMeshTopology, dEXO_RuptProperties)
      Type(EXO_Info)                                :: dEXO
      Type(MeshTopology_Info)                       :: dMeshTopology
      Type(EXO_RuptProperties_Type)                 :: dEXO_RuptProperties
      PetscInt                                      :: iErr
      PetscInt, Dimension(:), Pointer               :: Tmp_Prop

      If (dEXO%comm == PETSC_COMM_WORLD) Then
         SETERRQ(PETSC_ERR_SUP, 'EXO_RuptWrite: Writing properties on EXO files defined on PETSC_COMM_WORLD not implemented', iErr)
      End If
      
      dEXO%exoid = EXOPEN(dEXO%filename, EXWRIT, exo_cpu_ws, exo_io_ws, exo_ver, ierr)
      
      !!! Element block properties
      Allocate(Tmp_Prop(dMeshTopology%Num_Elem_Blks))
      Tmp_Prop = dEXO_RuptProperties%Is_Brittle
      Call EXPPA(dEXO%exoid, EXEBLK, Prop_Name_EB(1), Tmp_Prop, iErr)
      Tmp_Prop = dEXO_RuptProperties%Is_Domain
      Call EXPPA(dEXO%exoid, EXEBLK, Prop_Name_EB(2), Tmp_Prop, iErr)
      Tmp_Prop = dEXO_RuptProperties%Has_BodyForce
      Call EXPPA(dEXO%exoid, EXEBLK, Prop_Name_EB(3), Tmp_Prop, iErr)
      Call EXPPA(dEXO%exoid, EXEBLK, Prop_Name_EB(4), dEXO_RuptProperties%EB_BC_Type_X, iErr)
      Call EXPPA(dEXO%exoid, EXEBLK, Prop_Name_EB(5), dEXO_RuptProperties%EB_BC_Type_Y, iErr)
      Call EXPPA(dEXO%exoid, EXEBLK, Prop_Name_EB(6), dEXO_RuptProperties%EB_BC_Type_Z, iErr)
      DeAllocate(Tmp_Prop)

      !!! Side Set properties
      ! To Do
      
      !! Node Sets properties
      Call EXPPA(dEXO%exoid, EXNSET, Prop_Name_NS(1), dEXO_RuptProperties%NS_BC_Type_X, iErr)
      Call EXPPA(dEXO%exoid, EXNSET, Prop_Name_NS(2), dEXO_RuptProperties%NS_BC_Type_Y, iErr)
      Call EXPPA(dEXO%exoid, EXNSET, Prop_Name_NS(3), dEXO_RuptProperties%NS_BC_Type_Z, iErr)

      Call EXCLOS(dEXO%exoid, iErr)
      dEXO%exoid = 0
   End Subroutine EXO_RuptWrite
   
   Subroutine SchemeParamView(dSchemeParam, viewer)
      Type(SchemeParam_Type)                       :: dSchemeParam
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
   End Subroutine SchemeParamView

   Subroutine SchemeParamLoad(dSchemeParam, filename)
      Type(SchemeParam_Type)                       :: dSchemeParam
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
   End Subroutine SchemeParamLoad
   
   Subroutine SchemeParamGetFromOptions(dSchemeParam)
      Type(SchemeParam_Type)                       :: dSchemeParam
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
   End Subroutine SchemeParamGetFromOptions

  
   Subroutine GenHL_Iso2D_LambdaMu(lambda, mu, A) 
      Real(Kind = Kr), Intent(IN)         :: lambda, mu
      Type(Tens4OS2D), Intent(OUT)        :: A
      A = 0.0_Kr
      
      A%XXXX = lambda + 2.0_Kr * mu
      A%XXYY = lambda
      A%XYXY = mu
      A%YYYY = lambda + 2.0_Kr * mu
   End Subroutine GenHL_Iso2D_LambdaMu         

   Subroutine GenHL_Iso2D_EnuPlaneStress(E, nu, A) 
      Real(Kind = Kr), Intent(IN)         :: E, nu
      Type(Tens4OS2D), Intent(OUT)        :: A
      
      Real(Kind = Kr)                     :: Lambda, mu
      
      lambda = E * nu / (1.0_Kr - nu**2) 
      mu     = E / (1.0_Kr + nu) * .5_Kr
      A = 0.0_Kr
      A%XXXX = lambda + 2.0_Kr * mu
      A%XXYY = lambda
      A%XYXY = mu
      A%YYYY = lambda + 2.0_Kr * mu
   End Subroutine GenHL_Iso2D_EnuPlaneStress         
   
   Subroutine GenHL_Iso2D_EnuPlaneStrain(E, nu, A) 
      Real(Kind = Kr), Intent(IN)         :: E, nu
      Type(Tens4OS2D), Intent(OUT)        :: A
      
      Real(Kind = Kr)                     :: Lambda, mu
      
      lambda = E * nu / (1.0_Kr + nu) / (1.0_Kr - 2.0_Kr * nu)
      mu     = E / (1.0_Kr + nu) * .5_Kr      
      A = 0.0_Kr
      A%XXXX = lambda + 2.0_Kr * mu
      A%XXYY = lambda
      A%XYXY = mu
      A%YYYY = lambda + 2.0_Kr * mu
   End Subroutine GenHL_Iso2D_EnuPlaneStrain         

   Subroutine GenHL_Iso3D_LambdaMu(lambda, mu, A)
      Real(Kind = Kr), Intent(IN)         :: Lambda, Mu
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
      Real(Kind = Kr), Intent(IN)         :: E, nu
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
      Real(Kind = Kr), Intent(IN)         :: Lambda, mu1, mu2, theta
      Type(Tens4OS2D), Intent(OUT)        :: A
      
      A = 0.0_Kr
      A%XXXX = lambda + mu1 * (1.0_Kr + (cos(theta))**2) +  mu2 * (sin(theta))**2
      A%XXXY = (mu1-mu2) * cos(theta) * sin(theta)
      A%XXYY = lambda + (mu1-mu2) * (sin(theta))**2

      A%XYXY = mu1 * (sin(theta))**2 + mu2 * (cos(theta))**2
      A%XYYY = -A%XXXY

      A%YYYY =  A%XXXX
   End Subroutine GenHL_Ortho2D_LambdaMu
End Module m_Rupt_Struct
