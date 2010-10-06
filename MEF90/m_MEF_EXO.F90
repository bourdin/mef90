Module m_MEF_EXO
#include "finclude/petscdef.h"

   Use m_MEF_LinAlg
   Use m_MEF_Parameters
   Use m_MEF_Types
   Use m_MEF_MPI
   Use m_MEF_Elements
   Use m_MEF_Utils
   Use petsc
   
   IMPLICIT NONE

   Integer, Parameter, Public                        :: exo_cpu_ws = 8
   Integer, Parameter, Public                        :: exo_io_ws = 8
   PetscInt, Public                                  :: exo_ver

   
   Public :: Write_EXO_Case
   Public :: EXO_Check_Numbering
   Public :: Write_MeshTopology
   Public :: Write_MeshTopologyGlobal

   Public :: Read_EXO_Result_Global
   Public :: Write_EXO_Result_Global   
   Public :: Write_EXO_AllResult_Global   
   Public :: Read_EXO_Result_Vertex
   Public :: Write_EXO_Result_Vertex
   Public :: Read_EXO_Result_Cell
   Public :: Write_EXO_Result_Cell

   Public :: EXOProperty_Copy
   Public :: EXOProperty_Write
   Public :: EXOProperty_Ask
   Public :: EXOProperty_Read
   
   Public :: EXOVariable_Copy
   Public :: EXOVariable_Write
   Public :: EXOVariable_Read
   
   Interface Read_EXO_Result_Vertex
      Module Procedure Read_EXO_Result_VertexPtrInterlaced, Read_EXO_Result_VertexField, Read_EXO_Result_VertexSection, Read_EXO_Result_VertexVec, Read_EXO_Result_VertexVect2D, Read_EXO_Result_VertexVect3D, Read_EXO_Result_VertexMat2D, Read_EXO_Result_VertexMatS2D, Read_EXO_Result_VertexMat3D, Read_EXO_Result_VertexMatS3D
   End Interface Read_EXO_Result_Vertex
   !!! These are BROKEN! the number of values may not match the number of vertices (non P1)

   Interface Write_EXO_Result_Vertex
      Module Procedure Write_EXO_Result_VertexPtrInterlaced, Write_EXO_Result_VertexField, Write_EXO_Result_VertexSection, Write_EXO_Result_VertexVec, Write_EXO_Result_VertexVect2D, Write_EXO_Result_VertexVect3D, Write_EXO_Result_VertexMat2D, Write_EXO_Result_VertexMatS2D, Write_EXO_Result_VertexMat3D, Write_EXO_Result_VertexMatS3D
   End Interface Write_EXO_Result_Vertex
   
   Interface Read_EXO_Result_Cell
      Module Procedure Read_EXO_Result_CellPtrInterlaced, Read_EXO_Result_CellSection, Read_EXO_Result_CellVec, Read_EXO_Result_CellVect2D, Read_EXO_Result_CellVect3D, Read_EXO_Result_CellMat2D, Read_EXO_Result_CellMatS2D, Read_EXO_Result_CellMat3D, Read_EXO_Result_CellMatS3D
   End Interface Read_EXO_Result_Cell

   Interface Write_EXO_Result_Cell
      Module Procedure Write_EXO_Result_CellPtrInterlaced, Write_EXO_Result_CellSection, Write_EXO_Result_CellVec, Write_EXO_Result_CellVect2D, Write_EXO_Result_CellVect3D, Write_EXO_Result_CellMat2D, Write_EXO_Result_CellMatS2D, Write_EXO_Result_CellMat3D, Write_EXO_Result_CellMatS3D
   End Interface Write_EXO_Result_Cell

 Contains
    Subroutine Write_EXO_Case(prefix, formatstring, numprocs)
      Character(len=*)                               :: prefix, formatstring
      PetscInt                                       :: numprocs
      
      Character(len=MEF90_MXSTRLEN)                  :: casefile
      
      casefile = Trim(prefix)//'.e2c'
      If (MEF90_MyRank == 0) Then
         Open(unit = F_OUT, File = casefile, status = 'Unknown')
         Rewind(F_OUT)
         Write(F_OUT, 100) '#!EXODUS_CASE 1.0'
         Write(F_OUT, 101) numprocs
         Write(F_OUT, 102) Trim(prefix), Trim(formatstring)
         Close(F_OUT)
      End If
100 Format(A)
101 Format('FILES_PER_TIMESET', I4)
102 Format('TIMESET_TEMPLATE "', A, '-', A, '.gen"')
   End Subroutine Write_EXO_Case
   
   Subroutine EXO_Check_Numbering(dEXO, ErrorCode)
      Type(EXO_Type)                                 :: dEXO
      PetscInt, Intent(OUT)                          :: ErrorCode
      
      PetscInt                                       :: iErr
      Character(len=MEF90_MXSTRLEN)                  :: IOBuffer
      PetscInt                                       :: i
      PetscInt                                       :: NumEB, NumSS, NumNS
      PetscInt, Dimension(:), Pointer                :: IDs
      PetscInt                                       :: EXO_MyRank
      PetscReal                                      :: rDummy
      Character                                      :: cDummy
      PetscInt                                       :: vers
      
      ErrorCode = 0
      
      Call MPI_COMM_RANK(dEXO%Comm, EXO_MyRank, iErr)
      If (EXO_MyRank == 0) Then
         dEXO%exoid = EXOPEN(dEXO%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers, ierr)
         Call EXINQ(dEXO%exoid, EXELBL, NumEB, rDummy, cDummy, iErr)
         Call EXINQ(dEXO%exoid, EXSIDS, NumSS, rDummy, cDummy, iErr)
         Call EXINQ(dEXO%exoid, EXNODS, NumNS, rDummy, cDummy, iErr)
         
         If (NumEB > 0) Then
            Allocate(IDs(NumEB))
            Call EXGEBI(dEXO%exoid, IDs, iErr)   
            Do i = 1, NumEB
               if (IDs(i) /= i) Then            
                  ErrorCode = ErrorCode + 1
                  Write(IOBuffer, 100)  i, IDs(i)
                  Call PetscPrintf(PETSC_COMM_SELF, IOBuffer, iErr); CHKERRQ(iErr)
               End If
            End Do
            DeAllocate(IDs)
         End If
         
         If (NumSS > 0) Then
            Allocate(IDs(NumSS))
            Call EXGSSI(dEXO%exoid, IDs, iErr)   
            Do i = 1, NumSS
               if (IDs(i) /= i) Then            
                  ErrorCode = ErrorCode + 1
                  Write(IOBuffer, 101)  i, IDs(i)
                  Call PetscPrintf(PETSC_COMM_SELF, IOBuffer, iErr); CHKERRQ(iErr)
               End If
            End Do
            DeAllocate(IDs)
         End If
         
         If (NumNS > 0) Then         
            Allocate(IDs(NumNS))
            Call EXGNSI(dEXO%exoid, IDs, iErr)   
            Do i = 1, NumSS
               if (IDs(i) /= i) Then            
                  ErrorCode = ErrorCode + 1
                  Write(IOBuffer, 102)  i, IDs(i)
                  Call PetscPrintf(PETSC_COMM_SELF, IOBuffer, iErr); CHKERRQ(iErr)
               End If
            End Do
            DeAllocate(IDs)
         End If
         
         Call EXCLOS(dEXO%exoid, iErr)
         dEXO%exoid = 0
      End If
      Call MPI_BCast(ErrorCode, 1, MPI_INTEGER, 0, dEXO%Comm, iErr)
 100 Format('[ERROR] in MeshTopology_Check Numbering. Element block ', I3, ' index is ', I3, '\n')   
 101 Format('[ERROR] in MeshTopology_Check Numbering. Side Set ', I3, ' index is ', I3, '\n')   
 102 Format('[ERROR] in MeshTopology_Check Numbering. Node Set ', I3, ' index is ', I3, '\n')   
   End Subroutine EXO_Check_Numbering


   
   Subroutine EXOProperty_Copy(dEXO_in, dEXO_out)
      Type(EXO_Type)                                 :: dEXO_in, dEXO_out
      PetscInt                                       :: i
      
      dEXO_out%Num_EBProperties = dEXO_in%Num_EBProperties
      dEXO_out%Num_SSProperties = dEXO_in%Num_SSProperties
      dEXO_out%Num_NSProperties = dEXO_in%Num_NSProperties
      
      Allocate(dEXO_out%EBProperty(dEXO_out%Num_EBProperties))
      Do i = 1, dEXO_out%Num_EBProperties
         dEXO_out%EBProperty(i)%Name = dEXO_in%EBProperty(i)%Name
         Allocate(dEXO_out%EBProperty(i)%Value(Size(dEXO_in%EBProperty(i)%Value)))
         dEXO_out%EBProperty(i)%Value = dEXO_in%EBProperty(i)%Value
      End Do
      
      Allocate(dEXO_out%SSProperty(dEXO_out%Num_SSProperties))
      Do i = 1, dEXO_out%Num_SSProperties
         dEXO_out%SSProperty(i)%Name = dEXO_in%SSProperty(i)%Name
         Allocate(dEXO_out%SSProperty(i)%Value(Size(dEXO_in%SSProperty(i)%Value)))
         dEXO_out%SSProperty(i)%Value = dEXO_in%SSProperty(i)%Value
      End Do

      Allocate(dEXO_out%NSProperty(dEXO_out%Num_NSProperties))
      Do i = 1, dEXO_out%Num_NSProperties
         dEXO_out%NSProperty(i)%Name = dEXO_in%NSProperty(i)%Name
         Allocate(dEXO_out%NSProperty(i)%Value(Size(dEXO_in%NSProperty(i)%Value)))
         dEXO_out%NSProperty(i)%Value = dEXO_in%NSProperty(i)%Value
      End Do

   End Subroutine EXOProperty_Copy
   
   Subroutine EXOProperty_Write(dEXO)
      Type(EXO_Type)                                 :: dEXO
      PetscInt                                       :: vers
      PetscInt                                       :: iErr
      PetscInt                                       :: i

      PetscInt                                       :: NumEB, NumSS, NumNS
      PetscInt                                       :: EXO_MyRank
      PetscReal                                      :: rDummy
      Character                                      :: cDummy

      Call MPI_COMM_RANK(dEXO%Comm, EXO_MyRank, iErr)
      If (EXO_MyRank == 0) Then
         dEXO%exoid = EXOPEN(dEXO%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers, ierr)
         Call EXINQ(dEXO%exoid, EXELBL, NumEB, rDummy, cDummy, iErr)
         Call EXINQ(dEXO%exoid, EXSIDS, NumSS, rDummy, cDummy, iErr)
         Call EXINQ(dEXO%exoid, EXNODS, NumNS, rDummy, cDummy, iErr)

         !!! EB Properties
         If ((dEXO%Num_EBProperties > 0) .AND. (NumEB > 0)) Then
            Call EXPPN(dEXO%exoid, EXEBLK, dEXO%Num_EBProperties, dEXO%EBProperty(:)%Name, iErr)
            Do i = 1, dEXO%Num_EBProperties
               Call EXPPA(dEXO%exoid, EXEBLK, dEXO%EBProperty(i)%Name, dEXO%EBProperty(i)%Value, iErr)
            End Do
         End If
         !!! SS Properties
         If ((dEXO%Num_SSProperties > 0) .AND. (NumSS > 0)) Then
            Call EXPPN(dEXO%exoid, EXSSET, dEXO%Num_SSProperties, dEXO%SSProperty(:)%Name, iErr)
            Do i = 1, dEXO%Num_SSProperties
               Call EXPPA(dEXO%exoid, EXSSET, dEXO%SSProperty(i)%Name, dEXO%SSProperty(i)%Value, iErr)
            End Do
         End If
         !!! NS Properties
         If ((dEXO%Num_NSProperties > 0) .AND. (NumNS > 0)) Then
            Call EXPPN(dEXO%exoid, EXNSET, dEXO%Num_NSProperties, dEXO%NSProperty(:)%Name, iErr)         
            Do i = 1, dEXO%Num_NSProperties
               Call EXPPA(dEXO%exoid, EXNSET, dEXO%NSProperty(i)%Name, dEXO%NSProperty(i)%Value, iErr)
            End Do
         End If
         Call EXCLOS(dEXO%exoid, iErr)
         dEXO%exoid = 0
      End If
   End Subroutine EXOProperty_Write
   
   
   Subroutine EXOProperty_Ask(dEXO, dMeshTopology)
      Type(EXO_Type)                                 :: dEXO
      Type(MeshTopology_Type)                        :: dMeshTopology
      PetscInt                                       :: iErr
      PetscInt                                       :: i, j, IntBuffer

      PetscInt                                       :: NumEB, NumSS, NumNS
      PetscInt                                       :: EXO_MyRank
      Character(len=MEF90_MXSTRLEN)                  :: IOBuffer

      Do i = 1, dMeshTopology%Num_Elem_Blks_Global
         Write(IOBuffer, 100) i
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Do j = 1, dEXO%Num_EBProperties
            Write(IOBuffer, 110) dEXO%EBProperty(j)%Name
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
            If (MEF90_MyRank == 0) Then
               Read(*,*) dEXO%EBProperty(j)%Value(i)
            End If
            Call MPI_BCast(dEXO%EBProperty(j)%Value(i), 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, iErr)
         End Do
      End Do
      
      Do i = 1, dMeshTopology%Num_Side_Sets_Global
         Write(IOBuffer, 101) i
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Do j = 1, dEXO%Num_SSProperties
            Write(IOBuffer, 110) dEXO%SSProperty(j)%Name
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
            If (MEF90_MyRank == 0) Then
               Read(*,*) dEXO%SSProperty(j)%Value(i)
            End If
            Call MPI_BCast(dEXO%SSProperty(j)%Value(i), 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, iErr)
         End Do
      End Do

      Do i = 1, dMeshTopology%Num_Node_Sets_Global
         Write(IOBuffer, 102) i
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Do j = 1, dEXO%Num_NSProperties
            Write(IOBuffer, 110) dEXO%NSProperty(j)%Name
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
            If (MEF90_MyRank == 0) Then
               Read(*,*) dEXO%NSProperty(j)%Value(i)
            End If
            Call MPI_BCast(dEXO%NSProperty(j)%Value(i), 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, iErr)
         End Do
      End Do
 100 Format('*** Element Block ', T24, I3, '\n')
 101 Format('*** Side Set      ', T24, I3, '\n')
 102 Format('*** Node Set      ', T24, I3, '\n')
 110 Format(T24, A, T60, ': ')
   End Subroutine EXOProperty_Ask
      
   Subroutine EXOProperty_Read(dEXO)
      Type(EXO_Type)                                 :: dEXO
      PetscInt                                       :: vers
      PetscInt                                       :: iErr
      PetscInt                                       :: i
      PetscInt                                       :: NumProp, NumVal
      PetscReal                                      :: rDummy
      Character                                      :: cDummy
      
      Character(len=MXSTLN), Dimension(:), Pointer   :: PropName
      Integer                                        :: EXO_MyRank
      
      Call MPI_COMM_RANK(dEXO%Comm, EXO_MyRank, iErr)

      If (EXO_MyRank == 0) Then
         dEXO%exoid = EXOPEN(dEXO%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers, ierr)
         !!! EB Properties
         Call EXINQ(dEXO%exoid, EXNEBP, NumProp, rDummy, cDummy, iErr)
         Allocate(PropName(NumProp))
         Call EXGPN(dEXO%exoid, EXEBLK, PropName, iErr)
         !!!! id is always property 1, but we don.t want to count it
         dEXO%Num_EBProperties = max(0,NumProp-1)
         Allocate(dEXO%EBProperty(dEXO%Num_EBProperties))
         Call EXINQ(dEXO%exoid, EXELBL, NumVal, rDummy, cDummy, iErr)
         Do i = 1, dEXO%Num_EBProperties
            Allocate(dEXO%EBProperty(i)%Value(NumVal))
            dEXO%EBProperty(i)%Value(NumVal) = -1
            Call EXGPA(dEXO%exoid, EXEBLK, PropName(i+1), dEXO%EBProperty(i)%Value, iErr)
            dEXO%EBProperty(i)%Name  = PropName(i+1)
         End Do
         DeAllocate(PropName)
         
         !!! SS Properties
         Call EXINQ(dEXO%exoid, EXNSSP, NumProp, rDummy, cDummy, iErr)
         Allocate(PropName (NumProp))
         Call EXGPN(dEXO%exoid, EXSSET, PropName, iErr)
         dEXO%Num_SSProperties = max(0,NumProp-1)
         Call EXINQ(dEXO%exoid, EXSIDS, NumVal, rDummy, cDummy, iErr)
         Allocate(dEXO%SSProperty(dEXO%Num_SSProperties))
         Do i = 1, dEXO%Num_SSProperties
            Allocate(dEXO%SSProperty(i)%Value(NumVal))
            dEXO%SSProperty(i)%Value(NumVal) = -1
            Call EXGPA(dEXO%exoid, EXSSET, PropName(i+1), dEXO%SSProperty(i)%Value, iErr)
            dEXO%SSProperty(i)%Name  = PropName(i+1)
         End Do
         DeAllocate(PropName)

         !!! NS Properties
         Call EXINQ(dEXO%exoid, EXNNSP, NumProp, rDummy, cDummy, iErr)
         Allocate(PropName (NumProp))
         Call EXGPN(dEXO%exoid, EXNSET, PropName, iErr)
         dEXO%Num_NSProperties = max(0,NumProp-1)
         Call EXINQ(dEXO%exoid, EXNODS, NumVal, rDummy, cDummy, iErr)
         Allocate(dEXO%NSProperty(dEXO%Num_NSProperties))
         Do i = 1, dEXO%Num_NSProperties
            Allocate(dEXO%NSProperty(i)%Value(NumVal))
            Call EXGPA(dEXO%exoid, EXNSET, PropName(i+1), dEXO%NSProperty(i)%Value, iErr)
            dEXO%NSProperty(i)%Name  = PropName(i+1)
         End Do
         DeAllocate(PropName)
         
         Call EXCLOS(dEXO%exoid, iErr)
         dEXO%exoid = 0
      End If
      !!! Broadcast everything now
      Call MPI_BCast(dEXO%Num_EBProperties, 1, MPI_INTEGER, 0, dEXO%Comm, iErr)
      Call MPI_BCast(dEXO%Num_SSProperties, 1, MPI_INTEGER, 0, dEXO%Comm, iErr)
      Call MPI_BCast(dEXO%Num_NSProperties, 1, MPI_INTEGER, 0, dEXO%Comm, iErr)

      If (EXO_MyRank /= 0) Then
         Allocate(dEXO%EBProperty(dEXO%Num_EBProperties))
         Allocate(dEXO%SSProperty(dEXO%Num_SSProperties))
         Allocate(dEXO%NSProperty(dEXO%Num_NSProperties))
      EndIf
      !!! Element Blocks
      Do i = 1, dEXO%Num_EBProperties
         Call MPI_BCast(dEXO%EBProperty(i)%Name, MXSTLN, MPI_CHARACTER, 0, dEXO%Comm, iErr)
         NumProp = Size(dEXO%EBProperty(i)%Value)
         Call MPI_BCast(NumProp, 1, MPI_INTEGER, 0, dEXO%Comm, iErr)
         If (EXO_MyRank /= 0) Then
            Allocate(dEXO%EBProperty(i)%Value(NumProp))
         End If
         Call MPI_BCast(dEXO%EBProperty(i)%Value, NumProp, MPI_INTEGER, 0, dEXO%Comm, iErr)
      End Do
      !!! Side Sets
      Do i = 1, dEXO%Num_SSProperties
         Call MPI_BCast(dEXO%SSProperty(i)%Name, MXSTLN, MPI_CHARACTER, 0, dEXO%Comm, iErr)
         NumProp = Size(dEXO%SSProperty(i)%Value)
         Call MPI_BCast(NumProp, 1, MPI_INTEGER, 0, dEXO%Comm, iErr)
         If (EXO_MyRank /= 0) Then
            Allocate(dEXO%SSProperty(i)%Value(NumProp))
         End If
         Call MPI_BCast(dEXO%SSProperty(i)%Value, NumProp, MPI_INTEGER, 0, dEXO%Comm, iErr)
      End Do
      !!! Node Sets
      Do i = 1, dEXO%Num_NSProperties
         Call MPI_BCast(dEXO%NSProperty(i)%Name, MXSTLN, MPI_CHARACTER, 0, dEXO%Comm, iErr)
         NumProp = Size(dEXO%NSProperty(i)%Value)
         Call MPI_BCast(NumProp, 1, MPI_INTEGER, 0, dEXO%Comm, iErr)
         If (EXO_MyRank /= 0) Then
            Allocate(dEXO%NSProperty(i)%Value(NumProp))
         End If
         Call MPI_BCast(dEXO%NSProperty(i)%Value, NumProp, MPI_INTEGER, 0, dEXO%Comm, iErr)
      End Do
   End Subroutine EXOProperty_Read

   Subroutine EXOVariable_Copy(dEXO_in, dEXO_out)
      Type(EXO_Type)                                 :: dEXO_in, dEXO_out

      dEXO_out%Num_GlobVariables = dEXO_in%Num_GlobVariables
      dEXO_out%Num_CellVariables = dEXO_in%Num_CellVariables
      dEXO_out%Num_VertVariables = dEXO_in%Num_VertVariables
      
      If (dEXO_out%Num_GlobVariables > 0) Then
         Allocate(dEXO_out%GlobVariable(dEXO_out%Num_GlobVariables))
         dEXO_out%GlobVariable(:)%Name   = dEXO_in%GlobVariable(:)%Name
         dEXO_out%GlobVariable(:)%Offset = dEXO_in%GlobVariable(:)%Offset
      End If

      If (dEXO_out%Num_CellVariables > 0) Then
         Allocate(dEXO_out%CellVariable(dEXO_out%Num_CellVariables))
         dEXO_out%CellVariable(:)%Name   = dEXO_in%CellVariable(:)%Name
         dEXO_out%CellVariable(:)%Offset = dEXO_in%CellVariable(:)%Offset
      End If
      
      If (dEXO_out%Num_VertVariables > 0) Then
         Allocate(dEXO_out%VertVariable(dEXO_out%Num_VertVariables))
         dEXO_out%VertVariable(:)%Name   = dEXO_in%VertVariable(:)%Name
         dEXO_out%VertVariable(:)%Offset = dEXO_in%VertVariable(:)%Offset
      End If
   End Subroutine EXOVariable_Copy


   Subroutine EXOVariable_Write(dEXO)
      Type(EXO_Type)                                 :: dEXO
 
      PetscInt                                       :: i, iErr
      Integer                                        :: EXO_MyRank
      PetscInt                                       :: vers
      
      Call MPI_COMM_RANK(dEXO%Comm, EXO_MyRank, iErr)

      If (EXO_MyRank == 0) Then
         dEXO%exoid = EXOPEN(dEXO%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers, ierr)

         If (dEXO%Num_GlobVariables > 0) Then
            Call EXPVP (dEXO%exoid, 'g', dEXO%Num_GlobVariables, iErr)
            Call EXPVAN(dEXO%exoid, 'g', dEXO%Num_GlobVariables, dEXO%GlobVariable(:)%Name, iErr)
         End If
         If (dEXO%Num_CellVariables > 0) Then
            Call EXPVP (dEXO%exoid, 'e', dEXO%Num_CellVariables, iErr)
            Call EXPVAN(dEXO%exoid, 'e', dEXO%Num_CellVariables, dEXO%CellVariable(:)%Name, iErr)
         End If
         If (dEXO%Num_VertVariables > 0) Then
            Call EXPVP (dEXO%exoid, 'n', dEXO%Num_VertVariables, iErr)
            Call EXPVAN(dEXO%exoid, 'n', dEXO%Num_VertVariables, dEXO%VertVariable(:)%Name, iErr)
         End If
         Call EXCLOS(dEXO%exoid, iErr)
         dEXO%exoid = 0
      End If
   End Subroutine EXOVariable_Write 
   
   Subroutine EXOVariable_Read(dEXO)
      Type(EXO_Type)                                 :: dEXO
      PetscInt                                       :: vers
      PetscInt                                       :: iErr
      PetscInt                                       :: i
      PetscInt                                       :: NumProp, NumVal
      PetscReal                                      :: rDummy
      Character                                      :: cDummy
      
      Integer                                        :: EXO_MyRank
      
      Call MPI_COMM_RANK(dEXO%Comm, EXO_MyRank, iErr)

      If (EXO_MyRank == 0) Then
         dEXO%exoid = EXOPEN(dEXO%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers, ierr)

         Call EXGVP (dEXO%exoid, 'g', dEXO%Num_GlobVariables, iErr)
         If (dEXO%Num_GlobVariables > 0) Then
            Allocate(dEXO%GlobVariable(dEXO%Num_GlobVariables))      
            Call EXGVAN(dEXO%exoid, 'g', dEXO%Num_GlobVariables, dEXO%GlobVariable(:)%Name, iErr)
         End If
         
         Call EXGVP (dEXO%exoid, 'e', dEXO%Num_CellVariables, iErr)
         If (dEXO%Num_CellVariables > 0) Then
            Allocate(dEXO%CellVariable(dEXO%Num_CellVariables))      
            Call EXGVAN(dEXO%exoid, 'e', dEXO%Num_CellVariables, dEXO%CellVariable(:)%Name, iErr)
         End If

         Call EXGVP (dEXO%exoid, 'n', dEXO%Num_VertVariables, iErr)
         If (dEXO%Num_VertVariables > 0) Then
            Allocate(dEXO%VertVariable(dEXO%Num_VertVariables))      
            Call EXGVAN(dEXO%exoid, 'n', dEXO%Num_VertVariables, dEXO%VertVariable(:)%Name, iErr)
         End If
         
         Call EXCLOS(dEXO%exoid, iErr)
         dEXO%exoid = 0
      End If

      Call MPI_BCast(dEXO%Num_GlobVariables, 1, MPI_INTEGER, 0, dEXO%Comm, iErr)
      Call MPI_BCast(dEXO%Num_CellVariables, 1, MPI_INTEGER, 0, dEXO%Comm, iErr)
      Call MPI_BCast(dEXO%Num_VertVariables, 1, MPI_INTEGER, 0, dEXO%Comm, iErr)

      If (EXO_MyRank /= 0) Then
         If (dEXO%Num_GlobVariables > 0) Then
            Allocate(dEXO%GlobVariable(dEXO%Num_GlobVariables))      
         End If
         If (dEXO%Num_CellVariables > 0) Then
            Allocate(dEXO%CellVariable(dEXO%Num_CellVariables))      
         End If
         If (dEXO%Num_VertVariables > 0) Then
            Allocate(dEXO%VertVariable(dEXO%Num_VertVariables))      
         End If
      End If      

      Call MPI_BCast(dEXO%GlobVariable(:)%Offset, dEXO%Num_GlobVariables, MPI_INTEGER, 0, dEXO%Comm, iErr)
      Do i = 1, dEXO%Num_GlobVariables
         Call MPI_BCast(dEXO%GlobVariable(i)%Name, MXSTLN, MPI_CHARACTER, 0, dEXO%Comm, iErr)
      End Do
      Call MPI_BCast(dEXO%CellVariable(:)%Offset, dEXO%Num_CellVariables, MPI_INTEGER, 0, dEXO%Comm, iErr)
      Do i = 1, dEXO%Num_CellVariables
         Call MPI_BCast(dEXO%CellVariable(i)%Name, MXSTLN, MPI_CHARACTER, 0, dEXO%Comm, iErr)
      End Do
      Call MPI_BCast(dEXO%VertVariable(:)%Offset, dEXO%Num_VertVariables, MPI_INTEGER, 0, dEXO%Comm, iErr)
      Do i = 1, dEXO%Num_VertVariables
         Call MPI_BCast(dEXO%VertVariable(i)%Name, MXSTLN, MPI_CHARACTER, 0, dEXO%Comm, iErr)
      End Do
   End Subroutine EXOVariable_Read
   
   Subroutine Write_MeshTopology(dMeshTopology, dEXO)
      Type(MeshTopology_Type)                        :: dMeshTopology
      Type(EXO_Type)                                 :: dEXO
      PetscInt                                       :: vers
      PetscInt                                       :: iErr
      PetscInt                                       :: iDummy
      PetscReal                                      :: rDummy
      Character                                      :: cDummy
      
      PetscInt                                       :: iBlk, iSet, i
      PetscInt                                       :: iE, iELoc
      PetscInt                                       :: Num_Attr = 0         
      PetscInt                                       :: Num_Dist_Factor = 0
      
      Character(len=MXSTLN), Dimension(3)            :: Coord_Names, Elem_Type
      PetscReal, Dimension(:,:), Pointer             :: Coordinates
      PetscInt, Dimension(:,:), Pointer              :: ConnectMesh
      PetscInt, Dimension(:), Pointer                :: ConnectBlk
      PetscInt                                       :: offset
      
      Coord_Names(1) = 'X'
      Coord_Names(2) = 'Y'
      Coord_Names(3) = 'Z'
      
      Is_IO: If (dEXO%comm == PETSC_COMM_SELF) Then
         ! Open File
         dEXO%exoid = EXCRE (dEXO%filename, EXCLOB, exo_cpu_ws, exo_io_ws, iErr)
         dEXO%exoid = EXOPEN(dEXO%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers, ierr)
         
         ! Write Mesh Topology Info
         Call EXPINI(dEXO%exoid, dEXO%Title, dMeshTopology%Num_Dim, dMeshTopology%Num_Verts, dMeshTopology%Num_Elems, dMeshTopology%Num_Elem_Blks, dMeshTopology%Num_Node_Sets, dMeshTopology%Num_Side_Sets, iErr)
            
         ! Write Coordinate Names
         Call EXPCON(dEXO%exoid, Coord_Names, iErr)
         
         ! Write Elem blocks informations
         Do iBlk = 1, dMeshTopology%Num_Elem_Blks
            Select Case (dMeshTopology%elem_blk(iBlk)%Elem_Type)
            Case (MEF90_P1_Lagrange)
               Select Case (dMeshTopology%Num_Dim)
               Case (2)
                  Elem_Type = 'TRI3'
               Case (3)
                  Elem_Type = 'TETRA4'
               Case Default
                  SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, 'Only 2 and 3 dimensional elements are supported', iErr)
               End Select
            Case (MEF90_P2_Lagrange)
               Select Case (dMeshTopology%Num_Dim)
               Case (2)
                  Elem_Type = 'TRI6'
               Case (3)
                  Elem_Type = 'TETRA10'
               Case Default
                  SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, 'Only 2 and 3 dimensional elements are supported', iErr)
               End Select
            Case Default
                  SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, 'Only MEF90_P1_Lagrange and MEF90_P2_Lagrange elements are supported', iErr)
            End Select
            Call EXPELB(dEXO%exoid, dMeshTopology%elem_blk(iBlk)%ID, Elem_Type, dMeshTopology%elem_blk(iBlk)%Num_Elems, dMeshTopology%elem_blk(iBlk)%Num_DoF, Num_Attr, iErr)
         End Do
   
         ! Write Side sets informations
         If (dMeshTopology%Num_Side_Sets > 0) Then
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, 'Side sets not supported yet', iErr)
         End If
         
         ! Write Node sets informations
         Do iSet = 1, dMeshTopology%Num_Node_Sets
            Call EXPNP(dEXO%exoid, dMeshTopology%Node_Set(iSet)%ID, dMeshTopology%Node_Set(iSet)%Num_Nodes, Num_Dist_Factor, iErr)
            If (dMeshTopology%Node_Set(iSet)%Num_Nodes>0) Then
               Call EXPNS(dEXO%exoid, dMeshTopology%Node_Set(iSet)%ID, dMeshTopology%Node_Set(iSet)%Node_ID(:), iErr)
            End If
         End Do
   
         ! Write vertex coordinates
         Call MeshGetCoordinatesF90(dMeshTopology%mesh, Coordinates, iErr)
         Call EXPCOR(dEXO%exoid, Coordinates(:,1), Coordinates(:,2), Coordinates(:,3), iErr)
         Call MeshRestoreCoordinatesF90(dMeshTopology%mesh, Coordinates, iErr)
         
          ! Write Connectivity tables
         Call MeshGetElementsF90(dMeshTopology%mesh, ConnectMesh, iErr)
         Do iBlk = 1, dMeshTopology%Num_Elem_Blks
            If (dMeshTopology%Elem_Blk(iBlk)%Num_Elems > 0) Then
               Allocate (ConnectBlk(dMeshTopology%Elem_Blk(iBlk)%Num_Elems * dMeshTopology%Elem_Blk(iBlk)%Num_DoF))
               
               Do iELoc = 1, dMeshTopology%Elem_Blk(iBlk)%Num_Elems
                  iE = dMeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
                   offset = (iEloc-1) * dMeshTopology%Elem_Blk(iBlk)%Num_DoF
                   ConnectBlk(offset+1: offset+dMeshTopology%Elem_Blk(iBlk)%Num_DoF) = ConnectMesh(iE, :)
                End Do
                Call EXPELC (dEXO%exoid, dMeshTopology%Elem_Blk(iBlk)%ID, ConnectBlk, iErr)
                DeAllocate(ConnectBlk)
             End If
          End Do
          Call MeshRestoreElementsF90(dMeshTopology%mesh, ConnectMesh, iErr)

         Call EXCLOS(dEXO%exoid, iErr)
         dEXO%exoid = 0
      Else
         SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, 'Synchronized I/O on PETSC_COMM_WORLD not supported yet', ierr)
      End If Is_IO
   End Subroutine Write_MeshTopology
   
   
   Subroutine Write_MeshTopologyGlobal(dMeshTopology, dEXO, dGlobalComm)
      !!! Reconstruct the topology informatins for of all meshes in dGlobalComm then write it
      Type(MeshTopology_Type)                        :: dMeshTopology
      Type(EXO_Type)                                 :: dEXO
      MPI_Comm                                       :: dGlobalComm
      
      PetscInt                                       :: iBlk, iSet, i, iErr
      PetscInt                                       :: iE, iELoc
      PetscInt                                       :: Num_Attr = 0         
      PetscInt                                       :: Num_Dist_Factor = 0
      
      Type(MeshTopology_Type)                        :: GlobalMeshTopology
      PetscInt, Dimension(:), Pointer                :: Tmp_GlobalID, Tmp_ID


      ! Gather Global Sizes
      GlobalMeshTopology%num_dim   = dMeshTopology%num_dim
      GlobalMeshTopology%num_verts = dMeshTopology%num_verts
      GlobalMeshTopology%num_elems = dMeshTopology%num_elems
      
      Allocate(Tmp_ID(dMeshTopology%num_elem_blks))
      Tmp_ID = dMeshTopology%elem_blk(:)%ID
      Call Uniq(dGlobalComm, Tmp_ID, Tmp_GlobalID)            
      GlobalMeshTopology%Num_elem_blks = Size(Tmp_GlobalID)
      Allocate(GlobalMeshTopology%elem_blk(GlobalMeshTopology%num_elem_blks))
      GlobalMeshTopology%elem_blk(:)%ID = Tmp_GlobalID
      DeAllocate(Tmp_ID)

      GlobalMeshTopology%num_side_sets = 0

      Allocate(Tmp_ID(dMeshTopology%num_node_sets))
      Tmp_ID = dMeshTopology%node_set(:)%ID
      Call Uniq(dGlobalComm, Tmp_ID, Tmp_GlobalID)            
      GlobalMeshTopology%Num_node_sets = Size(Tmp_GlobalID)
      Allocate(GlobalMeshTopology%node_set(GlobalMeshTopology%num_node_sets))
      GlobalMeshTopology%node_set(:)%ID = Tmp_GlobalID
      DeAllocate(Tmp_ID)
      
      ! Element Blocks
      Allocate(Tmp_ID(GlobalMeshTopology%num_elem_blks))
      Tmp_ID = 0
      Do i = 1, GlobalMeshTopology%num_elem_blks
         GlobalMeshTopology%elem_blk(i)%num_Elems = 0
         GlobalMeshTopology%elem_blk(i)%DoF_Location(:) = 0
         GlobalMeshTopology%elem_blk(i)%num_DoF = 0
         Do iBlk = 1, dMeshTopology%num_elem_blks
            If (GlobalMeshTopology%elem_blk(i)%ID == dMeshTopology%elem_blk(iBlk)%ID) Then
               Tmp_ID(i) = dMeshTopology%elem_blk(iBlk)%Elem_Type
               GlobalMeshTopology%elem_blk(i)%num_elems = dMeshTopology%elem_blk(iBlk)%num_elems
               Allocate(GlobalMeshTopology%elem_blk(i)%elem_id(GlobalMeshTopology%elem_blk(i)%num_elems))
               GlobalMeshTopology%elem_blk(i)%elem_id = dMeshTopology%elem_blk(iBlk)%elem_id
               GlobalMeshTopology%elem_blk(i)%DoF_Location = dMeshTopology%elem_blk(iBlk)%DoF_Location
               GlobalMeshTopology%elem_blk(i)%num_DoF = dMeshTopology%elem_blk(iBlk)%num_DoF
            End If
         End Do
      End Do
      Call MPI_AllReduce(MPI_IN_PLACE, Tmp_ID, GlobalMeshTopology%num_elem_blks, MPI_INTEGER, MPI_MAX, dGlobalComm, iErr)
      GlobalMeshTopology%elem_blk(:)%Elem_Type = Tmp_ID
      DeAllocate(Tmp_ID)
      Do iBlk = 1, GlobalMeshTopology%num_elem_blks
         Call MPI_AllReduce(MPI_IN_PLACE, GlobalMeshTopology%elem_blk(iBlk)%DoF_Location, 4, MPI_INTEGER, MPI_MAX, dGlobalComm, iErr)
         Call MPI_AllReduce(MPI_IN_PLACE, GlobalMeshTopology%elem_blk(iBlk)%Num_DoF, 1, MPI_INTEGER, MPI_MAX, dGlobalComm, iErr)
      End Do
      
      ! Node Sets
      GlobalMeshTopology%node_set(:)%num_nodes = 0
      Do i = 1, GlobalMeshTopology%num_node_sets
         Do iSet = 1, dMeshTopology%num_node_sets
            If (GlobalMeshTopology%node_set(i)%ID == dMeshTopology%node_set(iSet)%ID) Then
               GlobalMeshTopology%node_set(i)%num_nodes = dMeshTopology%node_set(iSet)%num_nodes
               Allocate(GlobalMeshTopology%node_set(i)%node_ID(GlobalMeshTopology%node_set(i)%num_nodes))
               GlobalMeshTopology%node_set(i)%node_ID = dMeshTopology%node_set(iSet)%node_ID
            End If
         End Do
      End Do
      Do iSet = 1, GlobalMeshTopology%num_node_sets
      End Do

      GlobalMeshTopology%mesh = dMeshTopology%mesh
      Call Write_MeshTopology(GlobalMeshTopology, dEXO)
      Call MeshTopologyDestroy(GlobalMeshTopology)
   End Subroutine Write_MeshTopologyGlobal
   
!!!
!!! RESULT FILES
!!!

!!! READ
!!! In the sequel, we always assume that we are doing distributed I/O when EXO%comm == PETSC_COMM_SELF (i.e. 1 file per CPU) and sequential (i.e. IO operation on a single file on CPU 0) when EXO%comm == PETSC_COM_WORLD
   
   Subroutine Read_EXO_Result_Global(dEXO, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      PetscReal                                      :: dRes
      
      PetscReal                                      :: MyRes
      PetscInt                                       :: iErr
      PetscReal                                      :: Vers
      PetscReal, Dimension(:), Pointer               :: Tmp_Res
      PetscInt                                       :: Num_Vars, Num_TS
      PetscReal                                      :: fDum
      Character                                      :: cDum
      
      dRes = 0.0_Kr
      If ( ((dEXO%comm == PETSC_COMM_WORLD) .AND. (MEF90_MyRank == 0)) .OR. (dEXO%comm == PETSC_COMM_SELF) ) Then
         dEXO%exoid = EXOPEN(dEXO%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers, ierr)
         ! Get the number of global variables stored in the database    
         Call EXGVP(dEXO%exoid, 'G', Num_Vars, iErr)
         Allocate(Tmp_Res(Num_Vars))
         
         ! Read All the global variables at time step TS and copy the right one
         ! into Res
         Call EXGGV(dEXO%exoid, dTS, Num_Vars, Tmp_Res, iErr)
         MyRes = Tmp_Res(dIdx)
         DeAllocate (Tmp_Res)
         Call EXCLOS(dEXO%exoid, iErr)
         dEXO%exoid = 0
      End If
            
      !!! Broacast if dEXO%comm == PETSC_COMM_WORLD
      If (dEXO%comm == PETSC_COMM_WORLD) Then
         Call MPI_BCast(MyRes, 1, MPIU_SCALAR, 0, dEXO%comm, iErr)
      End If
      dRes = MyRes
   End Subroutine Read_EXO_Result_Global

!!!
!!! READ VERTEX BASED VARIABLES (NODAL VARIABLES)
!!!
   Subroutine Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      PetscReal, Dimension(:), Pointer               :: dRes
      
      PetscInt                                       :: Num_Rec, iRec
      PetscInt                                       :: iErr
      PetscReal                                      :: Vers
      
      dEXO%exoid = EXOPEN(dEXO%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers, ierr)
      
      If (Mod(Size(dRes), dMeshTopology%Num_Verts) /= 0) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Read_EXO_Result_VertexPtrInterlaced: The argument does not match the number of vertices in the mesh', iErr)
      End If
      Num_Rec = Size(dRes) / dMeshTopology%Num_Verts
      
      Do iRec = 1, Num_Rec
         Call EXGNV(dEXO%exoid, dTS, dIdx + iRec-1, dMeshTopology%Num_Verts, dRes(iRec:dMeshTopology%Num_Verts*Num_Rec:Num_Rec), iErr); CHKERRQ(iErr)
      End Do
      Call EXCLOS(dEXO%exoid, iErr)
      dEXO%exoid = 0
   End Subroutine Read_EXO_Result_VertexPtrInterlaced

   Subroutine Read_EXO_Result_VertexField(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Field)                                   :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr
      
      !!! We Assume that the section is initialized and has the proper size
      Call VecGetArrayF90(dRes%LocalVec, Res_Ptr, iErr); CHKERRQ(iErr)
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx, dTS, Res_Ptr)
      Call VecRestoreArrayF90(dRes%LocalVec, Res_Ptr, iErr); CHKERRQ(iErr)
      Call SectionRealToVec(dRes%Sec, dRes%Scatter, SCATTER_FORWARD, dRes%Vec, iErr); CHKERRQ(iErr)
   End Subroutine Read_EXO_Result_VertexField

   Subroutine Read_EXO_Result_VertexSection(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (SectionReal)                             :: dRes
      
      Type (Vec)                                     :: Res_Vec
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr
      
      !!! We Assume that the section is initialized and has the proper size
      Call SectionRealCreateLocalVector(dRes, Res_Vec, iErr); CHKERRQ(iErr)
      Call VecGetArrayF90(Res_Vec, Res_Ptr, iErr); CHKERRQ(iErr)
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx, dTS, Res_Ptr)
      Call VecRestoreArrayF90(Res_Vec, Res_Ptr, iErr); CHKERRQ(iErr)
      Call VecDestroy(Res_Vec, iErr); CHKERRQ(iErr)
   End Subroutine Read_EXO_Result_VertexSection

   Subroutine Read_EXO_Result_VertexVec(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Vec)                                     :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr
      
      !!! We Assume that the vector is initialized and has the proper size
      Call VecGetArrayF90(dRes, Res_Ptr, iErr); CHKERRQ(iErr)
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx, dTS, Res_Ptr)
      Call VecRestoreArrayF90(dRes, Res_Ptr, iErr); CHKERRQ(iErr)
   End Subroutine Read_EXO_Result_VertexVec

   Subroutine Read_EXO_Result_VertexVect2D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Vect2D), Dimension(:), Pointer           :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr
      
      If ( Size(dRes) /= dMeshTopology%Num_Verts) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Read_EXO_Result_VertexVect2D: The argument does not match the number of vertices in the mesh', iErr)
      End If

      Allocate(Res_Ptr(dMeshTopology%Num_Verts))
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      dRes(:)%X = Res_Ptr
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+1, dTS, Res_Ptr)
      dRes(:)%Y = Res_Ptr
      DeAllocate(Res_Ptr)
   End Subroutine Read_EXO_Result_VertexVect2D

   Subroutine Read_EXO_Result_VertexVect3D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Vect3D), Dimension(:), Pointer           :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr
      
      If ( Size(dRes) /= dMeshTopology%Num_Verts) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Read_EXO_Result_VertexVect3D: The argument does not match the number of vertices in the mesh', iErr)
      End If

      Allocate(Res_Ptr(dMeshTopology%Num_Verts))
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      dRes(:)%X = Res_Ptr
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+1, dTS, Res_Ptr)
      dRes(:)%Y = Res_Ptr
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+2, dTS, Res_Ptr)
      dRes(:)%Z = Res_Ptr
      DeAllocate(Res_Ptr)
   End Subroutine Read_EXO_Result_VertexVect3D

   Subroutine Read_EXO_Result_VertexMat2D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Mat2D), Dimension(:), Pointer            :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr
      
      If ( Size(dRes) /= dMeshTopology%Num_Verts) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Read_EXO_Result_VertexMat2D: The argument does not match the number of vertices in the mesh', iErr)
      End If

      Allocate(Res_Ptr(dMeshTopology%Num_Verts))
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      dRes(:)%XX = Res_Ptr
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+1, dTS, Res_Ptr)
      dRes(:)%XY = Res_Ptr
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+2, dTS, Res_Ptr)
      dRes(:)%YX = Res_Ptr
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+3, dTS, Res_Ptr)
      dRes(:)%YY = Res_Ptr
      DeAllocate(Res_Ptr)
   End Subroutine Read_EXO_Result_VertexMat2D

   Subroutine Read_EXO_Result_VertexMatS2D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (MatS2D), Dimension(:), Pointer           :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr
      
      If ( Size(dRes) /= dMeshTopology%Num_Verts) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Read_EXO_Result_VertexMatS2D: The argument does not match the number of vertices in the mesh', iErr)
      End If

      Allocate(Res_Ptr(dMeshTopology%Num_Verts))
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      dRes(:)%XX = Res_Ptr
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+1, dTS, Res_Ptr)
      dRes(:)%YY = Res_Ptr
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+2, dTS, Res_Ptr)
      dRes(:)%XY = Res_Ptr
      DeAllocate(Res_Ptr)
   End Subroutine Read_EXO_Result_VertexMatS2D

   Subroutine Read_EXO_Result_VertexMat3D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Mat3D), Dimension(:), Pointer            :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr
      
      If ( Size(dRes) /= dMeshTopology%Num_Verts) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Read_EXO_Result_VertexMat3D: The argument does not match the number of vertices in the mesh', iErr)
      End If

      Allocate(Res_Ptr(dMeshTopology%Num_Verts))
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      dRes(:)%XX = Res_Ptr
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+1, dTS, Res_Ptr)
      dRes(:)%XY = Res_Ptr
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+2, dTS, Res_Ptr)
      dRes(:)%XZ = Res_Ptr
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+3, dTS, Res_Ptr)
      dRes(:)%YX = Res_Ptr
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+4, dTS, Res_Ptr)
      dRes(:)%YY = Res_Ptr
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+5, dTS, Res_Ptr)
      dRes(:)%YZ = Res_Ptr
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+6, dTS, Res_Ptr)
      dRes(:)%ZX = Res_Ptr
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+7, dTS, Res_Ptr)
      dRes(:)%ZY = Res_Ptr
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+8, dTS, Res_Ptr)
      dRes(:)%ZZ = Res_Ptr
      DeAllocate(Res_Ptr)
   End Subroutine Read_EXO_Result_VertexMat3D

   Subroutine Read_EXO_Result_VertexMatS3D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (MatS3D), Dimension(:), Pointer           :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr
      
      If ( Size(dRes) /= dMeshTopology%Num_Verts) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Read_EXO_Result_VertexMatS3D: The argument does not match the number of vertices in the mesh', iErr)
      End If

      Allocate(Res_Ptr(dMeshTopology%Num_Verts))
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      dRes(:)%XX = Res_Ptr
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+1, dTS, Res_Ptr)
      dRes(:)%YY = Res_Ptr
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+2, dTS, Res_Ptr)
      dRes(:)%ZZ = Res_Ptr
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+3,   dTS, Res_Ptr)
      dRes(:)%YZ = Res_Ptr
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+4, dTS, Res_Ptr)
      dRes(:)%XZ = Res_Ptr
      Call Read_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+5, dTS, Res_Ptr)
      dRes(:)%XY = Res_Ptr
      DeAllocate(Res_Ptr)
   End Subroutine Read_EXO_Result_VertexMatS3D

!!!
!!! READ CELL BASED VARIABLE (ELEMENTAL VARIABLES)
!!!

   Subroutine Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      PetscReal, Dimension(:), Pointer               :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Tmp
      PetscInt                                       :: Num_Rec, iRec, iBlk, iE
      PetscInt                                       :: iErr
      PetscReal                                      :: Vers
      
      dEXO%exoid = EXOPEN(dEXO%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers, ierr)
      
      If (Mod(Size(dRes), dMeshTopology%Num_Elems) /= 0) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Read_EXO_Result_CellPtrInterlaced: The argument does not match the number of cells in the mesh', iErr)
      End If
      Num_Rec = Size(dRes) / dMeshTopology%Num_Elems

      Do_iBlk: Do iBlk = 1, dMeshTopology%Num_Elem_Blks      
         If (dMeshTopology%Elem_Blk(iBlk)%Num_elems > 0) Then
            Allocate(Res_Tmp(dMeshTopology%Elem_Blk(iBlk)%Num_elems))
            Do_iRec: Do iRec = 1, Num_Rec
               Call EXGEV(dEXO%exoid, dTS, dIdx + iRec-1, dMeshTopology%Elem_Blk(iBlk)%ID, dMeshTopology%Elem_Blk(iBlk)%Num_Elems, Res_Tmp, iErr); CHKERRQ(iErr)
               Do_iE: Do iE = 1, dMeshTopology%Elem_Blk(iBlk)%Num_Elems
                  dRes(iRec +  Num_Rec * (dMeshTopology%Elem_Blk(iBlk)%Elem_ID(iE)-1)) = Res_Tmp(iE)
               End Do Do_iE
            End Do Do_iRec
            DeAllocate(Res_Tmp)
         End If
      End Do Do_iBlk
      Call EXCLOS(dEXO%exoid, iErr)
      dEXO%exoid = 0
   End Subroutine Read_EXO_Result_CellPtrInterlaced

   Subroutine Read_EXO_Result_CellSection(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (SectionReal)                             :: dRes
      
      Type (Vec)                                     :: Res_Vec
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr
      
      !!! We Assume that the section is initialized and has the proper size
      Call SectionRealCreateLocalVector(dRes, Res_Vec, iErr); CHKERRQ(iErr)
      Call VecGetArrayF90(Res_Vec, Res_Ptr, iErr); CHKERRQ(iErr)
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx, dTS, Res_Ptr)
      Call VecRestoreArrayF90(Res_Vec, Res_Ptr, iErr); CHKERRQ(iErr)
      Call VecDestroy(Res_Vec, iErr); CHKERRQ(iErr)
   End Subroutine Read_EXO_Result_CellSection

   Subroutine Read_EXO_Result_CellVec(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Vec)                                     :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr
      
      !!! We Assume that the vector is initialized and has the proper size
      Call VecGetArrayF90(dRes, Res_Ptr, iErr); CHKERRQ(iErr)
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx, dTS, Res_Ptr)
      Call VecRestoreArrayF90(dRes, Res_Ptr, iErr); CHKERRQ(iErr)
   End Subroutine Read_EXO_Result_CellVec

   Subroutine Read_EXO_Result_CellVect2D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Vect2D), Dimension(:), Pointer           :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr
      
      If ( Size(dRes) /= dMeshTopology%Num_Verts) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Read_EXO_Result_CellVect2D: The argument does not match the number of cells in the mesh', iErr)
      End If

      Allocate(Res_Ptr(Size(dRes)))
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      dRes(:)%X = Res_Ptr
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+1, dTS, Res_Ptr)
      dRes(:)%Y = Res_Ptr
      DeAllocate(Res_Ptr)
   End Subroutine Read_EXO_Result_CellVect2D

   Subroutine Read_EXO_Result_CellVect3D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Vect3D), Dimension(:), Pointer           :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr
      
      If ( Size(dRes) /= dMeshTopology%Num_Verts) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Read_EXO_Result_CellVect3D: The argument does not match the number of cells in the mesh', iErr)
      End If

      Allocate(Res_Ptr(Size(dRes)))
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      dRes(:)%X = Res_Ptr
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+1, dTS, Res_Ptr)
      dRes(:)%Y = Res_Ptr
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+2, dTS, Res_Ptr)
      dRes(:)%Z = Res_Ptr
      DeAllocate(Res_Ptr)
   End Subroutine Read_EXO_Result_CellVect3D

   Subroutine Read_EXO_Result_CellMat2D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Mat2D), Dimension(:), Pointer            :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr
      
      If ( Size(dRes) /= dMeshTopology%Num_elems) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Read_EXO_Result_VertexMat2D: The argument does not match the number of cells in the mesh', iErr)
      End If

      Allocate(Res_Ptr(Size(dRes)))
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      dRes(:)%XX = Res_Ptr
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+1, dTS, Res_Ptr)
      dRes(:)%XY = Res_Ptr
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+2, dTS, Res_Ptr)
      dRes(:)%YX = Res_Ptr
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+3, dTS, Res_Ptr)
      dRes(:)%YY = Res_Ptr
      DeAllocate(Res_Ptr)
   End Subroutine Read_EXO_Result_CellMat2D

   Subroutine Read_EXO_Result_CellMatS2D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (MatS2D), Dimension(:), Pointer           :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr
      
      If ( Size(dRes) /= dMeshTopology%Num_elems) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Read_EXO_Result_CellMatS2D: The argument does not match the number of cells in the mesh', iErr)
      End If

      Allocate(Res_Ptr(Size(dRes)))
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      dRes(:)%XX = Res_Ptr
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+1, dTS, Res_Ptr)
      dRes(:)%YY = Res_Ptr
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+2, dTS, Res_Ptr)
      dRes(:)%XY = Res_Ptr
      DeAllocate(Res_Ptr)
   End Subroutine Read_EXO_Result_CellMatS2D

   Subroutine Read_EXO_Result_CellMat3D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Mat3D), Dimension(:), Pointer            :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr
      
      If ( Size(dRes) /= dMeshTopology%Num_Elems) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Read_EXO_Result_CellMat3D: The argument does not match the number of cells in the mesh', iErr)
      End If

      Allocate(Res_Ptr(Size(dRes)))
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      dRes(:)%XX = Res_Ptr
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+1, dTS, Res_Ptr)
      dRes(:)%XY = Res_Ptr
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+2, dTS, Res_Ptr)
      dRes(:)%XZ = Res_Ptr
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+3, dTS, Res_Ptr)
      dRes(:)%YX = Res_Ptr
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+4, dTS, Res_Ptr)
      dRes(:)%YY = Res_Ptr
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+5, dTS, Res_Ptr)
      dRes(:)%YZ = Res_Ptr
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+6, dTS, Res_Ptr)
      dRes(:)%ZX = Res_Ptr
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+7, dTS, Res_Ptr)
      dRes(:)%ZY = Res_Ptr
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+8, dTS, Res_Ptr)
      dRes(:)%ZZ = Res_Ptr
      DeAllocate(Res_Ptr)
   End Subroutine Read_EXO_Result_CellMat3D

   Subroutine Read_EXO_Result_CellMatS3D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (MatS3D), Dimension(:), Pointer           :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr
      
      If ( Size(dRes) /= dMeshTopology%Num_Elems) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Read_EXO_Result_CellMatS3D: The argument does not match the number of cells in the mesh', iErr)
      End If

      Allocate(Res_Ptr(Size(dRes)))
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      dRes(:)%XX = Res_Ptr
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+1, dTS, Res_Ptr)
      dRes(:)%YY = Res_Ptr
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+2, dTS, Res_Ptr)
      dRes(:)%ZZ = Res_Ptr
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+3,   dTS, Res_Ptr)
      dRes(:)%YZ = Res_Ptr
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+4, dTS, Res_Ptr)
      dRes(:)%XZ = Res_Ptr
      Call Read_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+5, dTS, Res_Ptr)
      dRes(:)%XY = Res_Ptr
      DeAllocate(Res_Ptr)
   End Subroutine Read_EXO_Result_CellMatS3D

!!! WRITE
   Subroutine Write_EXO_AllResult_Global(dEXO, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      PetscInt                                       :: dTS
      PetscReal, DImension(:), Pointer               :: dRes
      
      PetscInt                                       :: iErr
      PetscReal                                      :: Vers
      PetscInt                                       :: Num_Vars
      PetscReal                                      :: fDum
      Character                                      :: cDum
      
      If ( ((dEXO%comm == PETSC_COMM_WORLD) .AND. (MEF90_MyRank == 0)) .OR. (dEXO%comm == PETSC_COMM_SELF) ) Then
         dEXO%exoid = EXOPEN(dEXO%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers, ierr)
         Call EXGVP(dEXO%exoid, 'G', Num_Vars, iErr)
         If (Num_Vars /= Size(dRes)) Then
            SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Write_EXO_AllResult_Global: The argument does not match the number of global variables in the mesh', iErr)
         End If
         Call EXPGV(dEXO%exoid, dTS, Num_Vars, dRes, iErr)
         Call EXCLOS(dEXO%exoid, iErr)
         dEXO%exoid = 0
      End If
   End Subroutine Write_EXO_AllResult_Global

   Subroutine Write_EXO_Result_Global(dEXO, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      PetscReal                                      :: dRes
      
      PetscInt                                       :: iErr
      PetscReal                                      :: Vers
      PetscReal, Dimension(:), Pointer               :: Tmp_Res
      PetscInt                                       :: Num_Vars
      PetscReal                                      :: fDum
      Character                                      :: cDum
      
      If ( ((dEXO%comm == PETSC_COMM_WORLD) .AND. (MEF90_MyRank == 0)) .OR. (dEXO%comm == PETSC_COMM_SELF) ) Then
         dEXO%exoid = EXOPEN(dEXO%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers, ierr)
         ! Get the number of global variables stored in the database    
         Call EXGVP(dEXO%exoid, 'G', Num_Vars, iErr)
         Allocate(Tmp_Res(Num_Vars))
         Tmp_Res = 0.0_Kr
         
         ! Read All the global variables at time step TS into Tmp_Res
         ! Modify Tmp_Res(Idx) and write everything back...
         Call EXGGV(dEXO%exoid, dTS, Num_Vars, Tmp_Res, iErr)
         Tmp_Res(dIdx) = dRes
         
         Call EXPGV(dEXO%exoid, dTS, Num_Vars, Tmp_Res, iErr)
         DeAllocate (Tmp_Res)
         Call EXCLOS(dEXO%exoid, iErr)
         dEXO%exoid = 0
      End If
   End Subroutine Write_EXO_Result_Global

   Subroutine Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      PetscReal, Dimension(:), Pointer               :: dRes
      
      PetscInt                                       :: Num_Rec, iRec
      PetscInt                                       :: iErr
      PetscReal                                      :: Vers
      
      dEXO%exoid = EXOPEN(dEXO%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers, ierr)
      
      If (Mod(Size(dRes), dMeshTopology%Num_Verts) /= 0) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Write_EXO_Result_VertexPtrInterlaced: The argument does not match the number of vertices in the mesh', iErr)
      End If
      Num_Rec = Size(dRes) / dMeshTopology%Num_Verts
      Do iRec = 1, Num_Rec
         Call EXPNV(dEXO%exoid, dTS, dIdx + iRec-1, dMeshTopology%Num_Verts, dRes(iRec:dMeshTopology%Num_Verts*Num_Rec:Num_Rec), iErr); CHKERRQ(iErr)
      End Do
      Call EXCLOS(dEXO%exoid, iErr)
      dEXO%exoid = 0
   End Subroutine Write_EXO_Result_VertexPtrInterlaced

   Subroutine Write_EXO_Result_VertexField(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Field)                                   :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr

      Call VecGetArrayF90(dRes%LocalVec, Res_Ptr, iErr); CHKERRQ(iErr)
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx, dTS, Res_Ptr)
      Call VecRestoreArrayF90(dRes%LocalVec, Res_Ptr, iErr); CHKERRQ(iErr)
   End Subroutine Write_EXO_Result_VertexField

   Subroutine Write_EXO_Result_VertexSection(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (SectionReal)                             :: dRes
      
      Type (Vec)                                     :: Res_Vec
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr

      Call SectionRealCreateLocalVector(dRes, Res_Vec, iErr); CHKERRQ(iErr)
      Call VecGetArrayF90(Res_Vec, Res_Ptr, iErr); CHKERRQ(iErr)
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx, dTS, Res_Ptr)
      Call VecRestoreArrayF90(Res_Vec, Res_Ptr, iErr); CHKERRQ(iErr)
      Call VecDestroy(Res_Vec, iErr); CHKERRQ(iErr)
   End Subroutine Write_EXO_Result_VertexSection

   Subroutine Write_EXO_Result_VertexVec(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Vec)                                     :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr

      Call VecGetArrayF90(dRes, Res_Ptr, iErr); CHKERRQ(iErr)
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx, dTS, Res_Ptr)
      Call VecRestoreArrayF90(dRes, Res_Ptr, iErr); CHKERRQ(iErr)
   End Subroutine Write_EXO_Result_VertexVec

   Subroutine Write_EXO_Result_VertexVect2D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Vect2D), Dimension(:), Pointer           :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr

      If ( Size(dRes) /= dMeshTopology%Num_Verts ) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Write_EXO_Result_VertexVect2D: The argument size does not match the number of vertices in the mesh', iErr)
      End If
      Allocate(Res_Ptr(dMeshTopology%Num_Verts))
      Res_Ptr = dRes(:)%X
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%Y
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+1, dTS, Res_Ptr)
      DeAllocate(Res_Ptr)
   End Subroutine Write_EXO_Result_VertexVect2D

   Subroutine Write_EXO_Result_VertexVect3D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Vect3D), Dimension(:), Pointer           :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr

      If ( Size(dRes) /= dMeshTopology%Num_Verts ) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Write_EXO_Result_VertexVect2D: The argument size does not match the number of vertices in the mesh', iErr)
      End If
      Allocate(Res_Ptr(dMeshTopology%Num_Verts))
      Res_Ptr = dRes(:)%X
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%Y
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+1, dTS, Res_Ptr)
      Res_Ptr = dRes(:)%Z
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+2, dTS, Res_Ptr)
      DeAllocate(Res_Ptr)
   End Subroutine Write_EXO_Result_VertexVect3D

   Subroutine Write_EXO_Result_VertexMat2D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Mat2D), Dimension(:), Pointer            :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr

      If ( Size(dRes) /= dMeshTopology%Num_Verts ) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Write_EXO_Result_VertexMat2D: The argument size does not match the number of vertices in the mesh', iErr)
      End If
      Allocate(Res_Ptr(dMeshTopology%Num_Verts))
      Res_Ptr = dRes(:)%XX
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%XY
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+1, dTS, Res_Ptr)
      Res_Ptr = dRes(:)%YX
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+2,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%YY
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+3, dTS, Res_Ptr)
      DeAllocate(Res_Ptr)
   End Subroutine Write_EXO_Result_VertexMat2D

   Subroutine Write_EXO_Result_VertexMatS2D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (MatS2D), Dimension(:), Pointer           :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr

      If ( Size(dRes) /= dMeshTopology%Num_Verts ) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Write_EXO_Result_VertexMatS2D: The argument size does not match the number of vertices in the mesh', iErr)
      End If
      Allocate(Res_Ptr(dMeshTopology%Num_Verts))
      Res_Ptr = dRes(:)%XX
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%YY
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+1,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%XY
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+2, dTS, Res_Ptr)
      DeAllocate(Res_Ptr)
   End Subroutine Write_EXO_Result_VertexMatS2D

   Subroutine Write_EXO_Result_VertexMat3D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Mat3D), Dimension(:), Pointer            :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr

      If ( Size(dRes) /= dMeshTopology%Num_Verts ) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Write_EXO_Result_VertexMat3D: The argument size does not match the number of vertices in the mesh', iErr)
      End If
      Allocate(Res_Ptr(dMeshTopology%Num_Verts))
      Res_Ptr = dRes(:)%XX
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%XY
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+1, dTS, Res_Ptr)
      Res_Ptr = dRes(:)%XZ
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+2,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%YX
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+3,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%YY
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+4, dTS, Res_Ptr)
      Res_Ptr = dRes(:)%YZ
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+5,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%ZX
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+6,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%ZY
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+7, dTS, Res_Ptr)
      Res_Ptr = dRes(:)%ZZ
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+8,   dTS, Res_Ptr)
      DeAllocate(Res_Ptr)
   End Subroutine Write_EXO_Result_VertexMat3D

   Subroutine Write_EXO_Result_VertexMatS3D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (MatS3D), Dimension(:), Pointer           :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr

      If ( Size(dRes) /= dMeshTopology%Num_Verts ) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Write_EXO_Result_VertexMatS3D: The argument size does not match the number of vertices in the mesh', iErr)
      End If
      Allocate(Res_Ptr(dMeshTopology%Num_Verts))
      Res_Ptr = dRes(:)%XX
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%YY
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+1,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%ZZ
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+2, dTS, Res_Ptr)
      Res_Ptr = dRes(:)%YZ
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+3,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%XZ
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+4,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%XY
      Call Write_EXO_Result_VertexPtrInterlaced(dExo, dMeshTopology, dIdx+5, dTS, Res_Ptr)
      DeAllocate(Res_Ptr)
   End Subroutine Write_EXO_Result_VertexMatS3D

!!!
!!! WRITE CELL BASED VARIABLES (ELEMENTAL)
!!!   
   Subroutine Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      PetscReal, Dimension(:), Pointer               :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Tmp
      PetscInt                                       :: Num_Rec, iRec, iBlk, iE
      PetscInt                                       :: iErr
      PetscReal                                      :: Vers
      
      dEXO%exoid = EXOPEN(dEXO%filename, EXWRIT, exo_cpu_ws, exo_io_ws, vers, ierr)
      
      If (Mod(Size(dRes), dMeshTopology%Num_Elems) /= 0) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Read_EXO_Result_CellPtrInterlaced: The argument does not match the number of cells in the mesh', iErr)
      End If
      Num_Rec = Size(dRes) / dMeshTopology%Num_Elems
      
      Do_iBlk: Do iBlk = 1, dMeshTopology%Num_Elem_Blks      
         If (dMeshTopology%Elem_Blk(iBlk)%Num_elems > 0) Then
            Allocate(Res_Tmp(dMeshTopology%Elem_Blk(iBlk)%Num_elems))
            Do_iRec: Do iRec = 1, Num_Rec
               Do_iE: Do iE = 1, dMeshTopology%Elem_Blk(iBlk)%Num_Elems
                  Res_Tmp(iE) = dRes(iRec +  Num_Rec * (dMeshTopology%Elem_Blk(iBlk)%Elem_ID(iE)-1))
               End Do Do_iE
               Call EXPEV(dEXO%exoid, dTS, dIdx + iRec-1, dMeshTopology%Elem_Blk(iBlk)%ID, dMeshTopology%Elem_Blk(iBlk)%Num_Elems, Res_Tmp, iErr); CHKERRQ(iErr)
            End Do Do_iRec
            DeAllocate(Res_Tmp)
         End If
      End Do Do_iBlk
      Call EXCLOS(dEXO%exoid, iErr)
      dEXO%exoid = 0
   End Subroutine Write_EXO_Result_CellPtrInterlaced

   Subroutine Write_EXO_Result_CellSection(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (SectionReal)                             :: dRes
      
      Type (Vec)                                     :: Res_Vec
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr

      Call SectionRealCreateLocalVector(dRes, Res_Vec, iErr); CHKERRQ(iErr)
      Call VecGetArrayF90(Res_Vec, Res_Ptr, iErr); CHKERRQ(iErr)
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx, dTS, Res_Ptr)
      Call VecRestoreArrayF90(Res_Vec, Res_Ptr, iErr); CHKERRQ(iErr)
      Call VecDestroy(Res_Vec, iErr); CHKERRQ(iErr)
   End Subroutine Write_EXO_Result_CellSection

   Subroutine Write_EXO_Result_CellVec(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Vec)                                     :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr

      Call VecGetArrayF90(dRes, Res_Ptr, iErr); CHKERRQ(iErr)
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx, dTS, Res_Ptr)
      Call VecRestoreArrayF90(dRes, Res_Ptr, iErr); CHKERRQ(iErr)
   End Subroutine Write_EXO_Result_CellVec

   Subroutine Write_EXO_Result_CellVect2D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Vect2D), Dimension(:), Pointer           :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr

      If ( Size(dRes) /= dMeshTopology%Num_Elems ) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Write_EXO_Result_CellVect2D: The argument size does not match the number of cells in the mesh', iErr)
      End If
      Allocate(Res_Ptr(Size(dRes)))
      Res_Ptr = dRes(:)%X
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%Y
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+1, dTS, Res_Ptr)
      DeAllocate(Res_Ptr)
   End Subroutine Write_EXO_Result_CellVect2D

   Subroutine Write_EXO_Result_CellVect3D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Vect3D), Dimension(:), Pointer           :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr

      If ( Size(dRes) /= dMeshTopology%Num_Elems ) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Write_EXO_Result_CellVect2D: The argument size does not match the number of cells in the mesh', iErr)
      End If
      Allocate(Res_Ptr(Size(dRes)))
      Res_Ptr = dRes(:)%X
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%Y
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+1, dTS, Res_Ptr)
      Res_Ptr = dRes(:)%Z
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+2, dTS, Res_Ptr)
      DeAllocate(Res_Ptr)
   End Subroutine Write_EXO_Result_CellVect3D

   Subroutine Write_EXO_Result_CellMat2D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Mat2D), Dimension(:), Pointer            :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr

      If ( Size(dRes) /= dMeshTopology%Num_Elems ) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Write_EXO_Result_CellMat2D: The argument size does not match the number of cells in the mesh', iErr)
      End If
      Allocate(Res_Ptr(Size(dRes)))
      Res_Ptr = dRes(:)%XX
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%XY
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+1, dTS, Res_Ptr)
      Res_Ptr = dRes(:)%YX
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+2,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%YY
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+3, dTS, Res_Ptr)
      DeAllocate(Res_Ptr)
   End Subroutine Write_EXO_Result_CellMat2D

   Subroutine Write_EXO_Result_CellMatS2D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (MatS2D), Dimension(:), Pointer           :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr

      If ( Size(dRes) /= dMeshTopology%Num_Elems ) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Write_EXO_Result_CellMatS2D: The argument size does not match the number of cells in the mesh', iErr)
      End If
      Allocate(Res_Ptr(Size(dRes)))
      Res_Ptr = dRes(:)%XX
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%YY
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+1,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%XY
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+2, dTS, Res_Ptr)
      DeAllocate(Res_Ptr)
   End Subroutine Write_EXO_Result_CellMatS2D

   Subroutine Write_EXO_Result_CellMat3D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (Mat3D), Dimension(:), Pointer            :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr

      If ( Size(dRes) /= dMeshTopology%Num_Elems ) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Write_EXO_Result_CellMat3D: The argument size does not match the number of cells in the mesh', iErr)
      End If
      Allocate(Res_Ptr(Size(dRes)))
      Res_Ptr = dRes(:)%XX
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%XY
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+1, dTS, Res_Ptr)
      Res_Ptr = dRes(:)%XZ
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+2,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%YX
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+3,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%YY
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+4, dTS, Res_Ptr)
      Res_Ptr = dRes(:)%YZ
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+5,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%ZX
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+6,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%ZY
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+7, dTS, Res_Ptr)
      Res_Ptr = dRes(:)%ZZ
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+8,   dTS, Res_Ptr)
      DeAllocate(Res_Ptr)
   End Subroutine Write_EXO_Result_CellMat3D

   Subroutine Write_EXO_Result_CellMatS3D(dExo, dMeshTopology, dIdx, dTS, dRes)
      Type (EXO_Type), Intent(INOUT)                 :: dEXO
      Type (MeshTopology_Type)                       :: dMeshTopology
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      Type (MatS3D), Dimension(:), Pointer           :: dRes
      
      PetscReal, Dimension(:), Pointer               :: Res_Ptr
      PetscInt                                       :: iErr

      If ( Size(dRes) /= dMeshTopology%Num_Elems ) Then
         SETERRQ(dEXO%Comm, PETSC_ERR_ARG_SIZ, 'Write_EXO_Result_CellMatS3D: The argument size does not match the number of cells in the mesh', iErr)
      End If
      Allocate(Res_Ptr(Size(dRes)))
      Res_Ptr = dRes(:)%XX
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%YY
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+1,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%ZZ
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+2, dTS, Res_Ptr)
      Res_Ptr = dRes(:)%YZ
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+3,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%XZ
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+4,   dTS, Res_Ptr)
      Res_Ptr = dRes(:)%XY
      Call Write_EXO_Result_CellPtrInterlaced(dExo, dMeshTopology, dIdx+5, dTS, Res_Ptr)
      DeAllocate(Res_Ptr)
   End Subroutine Write_EXO_Result_CellMatS3D
End Module m_MEF_EXO

