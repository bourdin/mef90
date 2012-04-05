Module m_MEF_EXO
#include "finclude/petscdef.h"
   Use m_MEF_Parameters
   Use m_MEF_Utils
   Use petsc
   IMPLICIT NONE

   !Integer,Parameter,Public                        :: exo_cpu_ws = 8
   !Integer,Parameter,Public                        :: exo_io_ws = 8
   PetscInt,Public                                 :: exo_ver

   
   Public :: Write_EXO_Case
   Public :: EXO_Check_Numbering

   Public :: Read_EXO_Result_Global
   Public :: Write_EXO_Result_Global   
   Public :: Write_EXO_AllResult_Global   

   Public :: EXOProperty_Copy
   Public :: EXOProperty_Write
   Public :: EXOCellSetProperty_Ask,EXOVertexSetProperty_Ask
   Public :: EXOProperty_Read
   
   Public :: EXOVariable_Copy
   Public :: EXOVariable_Write
   Public :: EXOVariable_Read
   
   Public :: EXO_Type,EXO_Property_Type,EXO_Variable_Type
   
   Public :: EXOView
   Public :: EXOSetElementType_Load      

!!! To do:
!!! 1. Remove EXO type (filename and exoid will be in a exo viewer one day
!!! 2. rename EXO_Property_Type MEF90_RealProperty_Type
!!! 3. rename EXO_Variable_Type, MEF90_Variable_Type
   Type EXO_Type
      MPI_Comm                                     :: comm      
      Integer                                      :: exoid
      Character(len=MXLNLN)                        :: filename  
      Character(len=MXLNLN)                        :: title     
      ! QA DATAS
      !Integer                                      :: num_QA    
      !Character(len=MXSTLN),Dimension(:,:),Pointer :: QA_rec
      ! Properties
      PetscInt                                     :: Num_EBProperties
      Type(EXO_Property_Type),Dimension(:),Pointer :: EBProperty
      PetscInt                                     :: Num_NSProperties
      Type(EXO_Property_Type),Dimension(:),Pointer :: NSProperty
      ! Since I always use EXOProperties as boolean flags, these
      ! Could be replaced by ISes
      ! OR I could get rid of them and use PETSc options instead
      
      ! Variables
      PetscInt                                     :: Num_GlobVariables    
      Type(EXO_Variable_Type),Dimension(:),Pointer :: GlobVariable
      PetscInt                                     :: Num_CellVariables    
      Type(EXO_Variable_Type),Dimension(:),Pointer :: CellVariable
      PetscInt                                     :: Num_VertVariables    
      Type(EXO_Variable_Type),Dimension(:),Pointer :: VertVariable
   End Type EXO_Type
   
   Type EXO_Property_Type
      !!! Derived type used to store the values of a property at ALL NS,EB
      !!! in an EXO file
      Character(MXSTLN)                            :: Name
      PetscInt,Dimension(:),Pointer                :: Value
   End Type EXO_Property_Type
   
   Type EXO_Variable_Type
      !!! Links variable name and offset in an exodus file
      Character(MXSTLN)                            :: Name
      PetscInt                                     :: Offset
      !!! the position of the variable in the exo file
   End Type EXO_Variable_Type 
   
     
Contains
#undef __FUNCT__
#define __FUNCT__ "EXOSetElementType_Load"
   Subroutine EXOSetElementType_Load(exoID,setID,EXOElemType)
      Integer,Intent(IN)                             :: exoID
      PetscInt,Intent(IN)                            :: setID
      Character(len=*),Intent(OUT)                   :: EXOElemType
      
      Integer                                        :: junk1,junk2,junk3,ierr
      
      Call EXGELB(exoID,setID,EXOElemType,junk1,junk2,junk3,ierr)
   End Subroutine EXOSetElementType_Load
      
#undef __FUNCT__
#define __FUNCT__ "Write_EXO_Case"
   Subroutine Write_EXO_Case(prefix,formatstring,numprocs)
      Character(len=*)                               :: prefix,formatstring
      PetscInt                                       :: numprocs
      
      Character(len=MEF90_MXSTRLEN)                  :: casefile
      
      casefile = Trim(prefix)//'.e2c'
      If (MEF90_MyRank == 0) Then
         Open(unit = F_OUT,File = casefile,status = 'Unknown')
         Rewind(F_OUT)
         Write(F_OUT,100) '#!EXODUS_CASE 1.0'
         Write(F_OUT,101) numprocs
         Write(F_OUT,102) Trim(prefix),Trim(formatstring)
         Close(F_OUT)
      End If
100 Format(A)
101 Format('FILES_PER_TIMESET',I4)
102 Format('TIMESET_TEMPLATE "',A,'-',A,'.gen"')
   End Subroutine Write_EXO_Case
   
#undef __FUNCT__
#define __FUNCT__ "EXO_Check_Numbering"
!!! I don't think that this is required anymore
!!! It is only needed if we don;t want to seek for material properties

   Subroutine EXO_Check_Numbering(dEXO,ErrorCode)
      Type(EXO_Type)                                 :: dEXO
      PetscInt,Intent(OUT)                           :: ErrorCode
      
      PetscInt                                       :: ierr
      Character(len=MEF90_MXSTRLEN)                  :: IOBuffer
      PetscInt                                       :: i
      PetscInt                                       :: NumEB,NumNS
      PetscInt,Dimension(:),Pointer                  :: IDs
      PetscInt                                       :: EXO_MyRank
      PetscReal                                      :: rDummy
      Character                                      :: cDummy
      
      ErrorCode = 0
      
      Call MPI_COMM_RANK(dEXO%Comm,EXO_MyRank,ierr)
      If (EXO_MyRank == 0) Then
         If (dEXO%exoid == 0) Then
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_NULL,"[ERROR]: Exodus file not open before IO operations\n",ierr)
         End If
         Call EXINQ(dEXO%exoid,EXELBL,NumEB,rDummy,cDummy,ierr)
         Call EXINQ(dEXO%exoid,EXNODS,NumNS,rDummy,cDummy,ierr)
         
         If (NumEB > 0) Then
            Allocate(IDs(NumEB))
            Call EXGEBI(dEXO%exoid,IDs,ierr)   
            Do i = 1,NumEB
               if (IDs(i) /= i) Then            
                  ErrorCode = ErrorCode + 1
                  Write(IOBuffer,100)  i,IDs(i)
                  Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr)
               End If
            End Do
            DeAllocate(IDs)
         End If
         
         If (NumNS > 0) Then         
            Allocate(IDs(NumNS))
            Call EXGNSI(dEXO%exoid,IDs,ierr)   
            Do i = 1,NumNS
               if (IDs(i) /= i) Then            
                  ErrorCode = ErrorCode + 1
                  Write(IOBuffer,102)  i,IDs(i)
                  Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr)
               End If
            End Do
            DeAllocate(IDs)
         End If
      End If
      Call MPI_BCast(ErrorCode,1,MPIU_INTEGER,0,dEXO%Comm,ierr)
 100 Format('[ERROR] in EXO_Check_Numbering. Element block ',I3,' index is ',I3,'\n')   
 102 Format('[ERROR] in EXO_Check_Numbering. Node Set ',I3,' index is ',I3,'\n')   
   End Subroutine EXO_Check_Numbering


#undef __FUNCT__
#define __FUNCT__ "EXOProperty_Copy"
   Subroutine EXOProperty_Copy(dEXO_in,dEXO_out)
      Type(EXO_Type)                                 :: dEXO_in,dEXO_out

      PetscInt                                       :: i
      
      dEXO_out%Num_EBProperties = dEXO_in%Num_EBProperties
      dEXO_out%Num_NSProperties = dEXO_in%Num_NSProperties
      
      Allocate(dEXO_out%EBProperty(dEXO_out%Num_EBProperties))
      Do i = 1,dEXO_out%Num_EBProperties
         dEXO_out%EBProperty(i)%Name = dEXO_in%EBProperty(i)%Name
         Allocate(dEXO_out%EBProperty(i)%Value(Size(dEXO_in%EBProperty(i)%Value)))
         dEXO_out%EBProperty(i)%Value = dEXO_in%EBProperty(i)%Value
      End Do
      
      Allocate(dEXO_out%NSProperty(dEXO_out%Num_NSProperties))
      Do i = 1,dEXO_out%Num_NSProperties
         dEXO_out%NSProperty(i)%Name = dEXO_in%NSProperty(i)%Name
         Allocate(dEXO_out%NSProperty(i)%Value(Size(dEXO_in%NSProperty(i)%Value)))
         dEXO_out%NSProperty(i)%Value = dEXO_in%NSProperty(i)%Value
      End Do
   End Subroutine EXOProperty_Copy
   
#undef __FUNCT__
#define __FUNCT__ "EXOProperty_Write"
   Subroutine EXOProperty_Write(dEXO)
      Type(EXO_Type)                                 :: dEXO

      PetscInt                                       :: ierr
      PetscInt                                       :: i
      PetscInt                                       :: NumEB,NumNS
      PetscInt                                       :: EXO_MyRank
      PetscReal                                      :: rDummy
      Character                                      :: cDummy

      Call MPI_COMM_RANK(dEXO%Comm,EXO_MyRank,ierr)
      If (EXO_MyRank == 0) Then
         If (dEXO%exoid == 0) Then
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_NULL,"[ERROR]: Exodus file not open before IO operations\n",ierr)
         End If
         Call EXINQ(dEXO%exoid,EXELBL,NumEB,rDummy,cDummy,ierr)
         Call EXINQ(dEXO%exoid,EXNODS,NumNS,rDummy,cDummy,ierr)

         !!! EB Properties
         If ((dEXO%Num_EBProperties > 0) .AND. (NumEB > 0)) Then
            Call EXPPN(dEXO%exoid,EXEBLK,dEXO%Num_EBProperties,dEXO%EBProperty(:)%Name,ierr)
            Do i = 1,dEXO%Num_EBProperties
               Call EXPPA(dEXO%exoid,EXEBLK,dEXO%EBProperty(i)%Name,dEXO%EBProperty(i)%Value,ierr)
            End Do
         End If
         !!! NS Properties
         If ((dEXO%Num_NSProperties > 0) .AND. (NumNS > 0)) Then
            Call EXPPN(dEXO%exoid,EXNSET,dEXO%Num_NSProperties,dEXO%NSProperty(:)%Name,ierr)         
            Do i = 1,dEXO%Num_NSProperties
               Call EXPPA(dEXO%exoid,EXNSET,dEXO%NSProperty(i)%Name,dEXO%NSProperty(i)%Value,ierr)
            End Do
         End If
      End If
   End Subroutine EXOProperty_Write
   
#undef __FUNCT__
#define __FUNCT__ "EXOCellSetProperty_Ask"   
   Subroutine EXOCellSetProperty_Ask(dEXO,cellSetGlobalIS)
      Type(EXO_Type)                                 :: dEXO
      Type(IS)                                       :: cellSetGlobalIS

      PetscInt                                       :: ierr
      PetscInt                                       :: i,j
      PetscInt                                       :: numCellSetGlobal
      Character(len=MEF90_MXSTRLEN)                  :: IOBuffer

      Call ISGetLocalSize(cellSetGlobalIS,numCellSetGlobal,ierr);CHKERRQ(ierr)
      Do i = 1,numCellSetGlobal 
         Write(IOBuffer,100) i
         Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
         Do j = 1,dEXO%Num_EBProperties
            Write(IOBuffer,110) dEXO%EBProperty(j)%Name
            Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
            If (MEF90_MyRank == 0) Then
               Read(*,*) dEXO%EBProperty(j)%Value(i)
            End If
            Call MPI_BCast(dEXO%EBProperty(j)%Value(i),1,MPIU_INTEGER,0,PETSC_COMM_WORLD,ierr)
         End Do
      End Do
 100 Format('*** Cell Set ',T24,I3,'\n')
 110 Format(T24,A,T60,': ')
   End Subroutine EXOCellSetProperty_Ask    
   
#undef __FUNCT__
#define __FUNCT__ "EXOCellSetProperty_AskWithBatch"
   Subroutine EXOCellSetProperty_AskWithBatch(dEXO,cellSetGlobalIS,BatchUnit,IsBatch)
      Type(EXO_Type)                                 :: dEXO
      Type(IS)                                       :: cellSetGlobalIS
      PetscInt                                       :: BatchUnit
      PetscBool                                      :: IsBatch

      PetscInt                                       :: ierr
      PetscInt                                       :: i,j
      PetscInt                                       :: numCellSetGlobal

      Character(len=MEF90_MXSTRLEN)                  :: IOBuffer

      Call ISGetLocalSize(cellSetGlobalIS,numCellSetGlobal,ierr);CHKERRQ(ierr)
      Do i = 1,numCellSetGlobal
         Write(IOBuffer,100) i
         Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
         Do j = 1,dEXO%Num_EBProperties
            Write(IOBuffer,200) i,Trim(dEXO%EBProperty(j)%Name)
            Call MEF90_AskInt(dEXO%EBProperty(j)%Value(i),IOBuffer,BatchUnit,IsBatch)
         End Do
         If ((.NOT. IsBatch) .AND. (MEF90_MyRank == 0)) Then
            Write(BatchUnit,*)
         End If
      End Do
100 Format('    Element Block ',T24,I3,'\n')
200 Format('EB',I4.4,': ',A)
   End Subroutine EXOCellSetProperty_AskWithBatch

#undef __FUNCT__
#define __FUNCT__ "EXOVertexSetProperty_Ask"
   Subroutine EXOVertexSetProperty_Ask(dEXO,vertexSetGlobalIS)
      Type(EXO_Type)                                 :: dEXO
      Type(IS)                                       :: vertexSetGlobalIS

      PetscInt                                       :: ierr
      PetscInt                                       :: i,j
      PetscInt                                       :: numvertexSetGlobal
      Character(len=MEF90_MXSTRLEN)                  :: IOBuffer

      Call ISGetLocalSize(vertexSetGlobalIS,numVertexSetGlobal,ierr);CHKERRQ(ierr)
      Do i = 1,numVertexSetGlobal
         Write(IOBuffer,102) i
         Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
         Do j = 1,dEXO%Num_NSProperties
            Write(IOBuffer,110) dEXO%NSProperty(j)%Name
            Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
            If (MEF90_MyRank == 0) Then
               Read(*,*) dEXO%NSProperty(j)%Value(i)
            End If
            Call MPI_BCast(dEXO%NSProperty(j)%Value(i),1,MPIU_INTEGER,0,PETSC_COMM_WORLD,ierr)
         End Do
      End Do
 102 Format('*** Vertex Set      ',T24,I3,'\n')
 110 Format(T24,A,T60,': ')
   End Subroutine EXOVertexSetProperty_Ask
      
#undef __FUNCT__
#define __FUNCT__ "EXOVertexProperty_AskWithBatch"
   Subroutine EXOVertexProperty_AskWithBatch(dEXO,vertexSetGlobalIS,BatchUnit,IsBatch)
      Type(EXO_Type)                                 :: dEXO
      Type(IS)                                       :: vertexSetGlobalIS
      PetscInt                                       :: BatchUnit
      PetscBool                                      :: IsBatch

      PetscInt                                       :: ierr
      PetscInt                                       :: i,j
      PetscInt                                       :: numVertexSetGlobal

      Character(len=MEF90_MXSTRLEN)                  :: IOBuffer
   
      Call ISGetLocalSize(vertexSetGlobalIS,numVertexSetGlobal,ierr);CHKERRQ(ierr)
      Do i = 1,numVertexSetGlobal
         Write(IOBuffer,102) i
         Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
         Do j = 1,dEXO%Num_NSProperties
            Write(IOBuffer,202) i,Trim(dEXO%NSProperty(j)%Name)
            Call MEF90_AskInt(dEXO%NSProperty(j)%Value(i),IOBuffer,BatchUnit,IsBatch)
         End Do
         If ((.NOT. IsBatch) .AND. (MEF90_MyRank == 0)) Then
            Write(BatchUnit,*)
         End If
      End Do
102 Format('    Node Set      ',T24,I3,'\n')
202 Format('NS',I4.4,': ',A)
   End Subroutine EXOVertexProperty_AskWithBatch

#undef __FUNCT__
#define __FUNCT__ "EXOProperty_Read"
   Subroutine EXOProperty_Read(dEXO)
      Type(EXO_Type)                                 :: dEXO

      PetscInt                                       :: ierr
      PetscInt                                       :: i
      PetscInt                                       :: NumProp,NumVal
      PetscReal                                      :: rDummy
      Character                                      :: cDummy
      PetscInt                                       :: EXO_MyRank
      Character(len=MXSTLN),Dimension(:),Pointer     :: PropName
      
      Call MPI_COMM_RANK(dEXO%Comm,EXO_MyRank,ierr)

      If (EXO_MyRank == 0) Then
         If (dEXO%exoid == 0) Then
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_NULL,"[ERROR]: Exodus file not open before IO operations\n",ierr)
         End If
         !!! EB Properties
         Call EXINQ(dEXO%exoid,EXNEBP,NumProp,rDummy,cDummy,ierr)
         Allocate(PropName(NumProp))
         Call EXGPN(dEXO%exoid,EXEBLK,PropName,ierr)
         !!!! id is always property 1,but we don.t want to count it
         dEXO%Num_EBProperties = max(0,NumProp-1)
         Allocate(dEXO%EBProperty(dEXO%Num_EBProperties))
         Call EXINQ(dEXO%exoid,EXELBL,NumVal,rDummy,cDummy,ierr)
         Do i = 1,dEXO%Num_EBProperties
            Allocate(dEXO%EBProperty(i)%Value(NumVal))
            dEXO%EBProperty(i)%Value(NumVal) = -1
            Call EXGPA(dEXO%exoid,EXEBLK,PropName(i+1),dEXO%EBProperty(i)%Value,ierr)
            dEXO%EBProperty(i)%Name  = PropName(i+1)
         End Do
         DeAllocate(PropName)
         
         !!! NS Properties
         Call EXINQ(dEXO%exoid,EXNNSP,NumProp,rDummy,cDummy,ierr)
         Allocate(PropName (NumProp))
         Call EXGPN(dEXO%exoid,EXNSET,PropName,ierr)
         dEXO%Num_NSProperties = max(0,NumProp-1)
         Call EXINQ(dEXO%exoid,EXNODS,NumVal,rDummy,cDummy,ierr)
         Allocate(dEXO%NSProperty(dEXO%Num_NSProperties))
         Do i = 1,dEXO%Num_NSProperties
            Allocate(dEXO%NSProperty(i)%Value(NumVal))
            Call EXGPA(dEXO%exoid,EXNSET,PropName(i+1),dEXO%NSProperty(i)%Value,ierr)
            dEXO%NSProperty(i)%Name  = PropName(i+1)
         End Do
         DeAllocate(PropName)
      End If
      !!! Broadcast everything now
      Call MPI_BCast(dEXO%Num_EBProperties,1,MPIU_INTEGER,0,dEXO%Comm,ierr)
      Call MPI_BCast(dEXO%Num_NSProperties,1,MPIU_INTEGER,0,dEXO%Comm,ierr)

      If (EXO_MyRank /= 0) Then
         Allocate(dEXO%EBProperty(dEXO%Num_EBProperties))
         Allocate(dEXO%NSProperty(dEXO%Num_NSProperties))
      EndIf
      !!! Element Blocks
      Do i = 1,dEXO%Num_EBProperties
         Call MPI_BCast(dEXO%EBProperty(i)%Name,MXSTLN,MPI_CHARACTER,0,dEXO%Comm,ierr)
         NumProp = Size(dEXO%EBProperty(i)%Value)
         Call MPI_BCast(NumProp,1,MPIU_INTEGER,0,dEXO%Comm,ierr)
         If (EXO_MyRank /= 0) Then
            Allocate(dEXO%EBProperty(i)%Value(NumProp))
         End If
         Call MPI_BCast(dEXO%EBProperty(i)%Value,NumProp,MPIU_INTEGER,0,dEXO%Comm,ierr)
      End Do
      !!! Node Sets
      Do i = 1,dEXO%Num_NSProperties
         Call MPI_BCast(dEXO%NSProperty(i)%Name,MXSTLN,MPI_CHARACTER,0,dEXO%Comm,ierr)
         NumProp = Size(dEXO%NSProperty(i)%Value)
         Call MPI_BCast(NumProp,1,MPIU_INTEGER,0,dEXO%Comm,ierr)
         If (EXO_MyRank /= 0) Then
            Allocate(dEXO%NSProperty(i)%Value(NumProp))
         End If
         Call MPI_BCast(dEXO%NSProperty(i)%Value,NumProp,MPIU_INTEGER,0,dEXO%Comm,ierr)
      End Do
   End Subroutine EXOProperty_Read

#undef __FUNCT__
#define __FUNCT__ "EXOVariable_Copy"
   Subroutine EXOVariable_Copy(dEXO_in,dEXO_out)
      Type(EXO_Type)                                 :: dEXO_in,dEXO_out

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


#undef __FUNCT__
#define __FUNCT__ "EXOVariable_Write"
   Subroutine EXOVariable_Write(dEXO)
      Type(EXO_Type)                                 :: dEXO
 
      PetscInt                                       :: ierr
      Integer                                        :: EXO_MyRank
      
      Call MPI_COMM_RANK(dEXO%Comm,EXO_MyRank,ierr)

      If (EXO_MyRank == 0) Then
         If (dEXO%exoid == 0) Then
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_NULL,"[ERROR]: Exodus file not open before IO operations\n",ierr)
         End If

         If (dEXO%Num_GlobVariables > 0) Then
            Call EXPVP (dEXO%exoid,'g',dEXO%Num_GlobVariables,ierr)
            Call EXPVAN(dEXO%exoid,'g',dEXO%Num_GlobVariables,dEXO%GlobVariable(:)%Name,ierr)
         End If
         If (dEXO%Num_CellVariables > 0) Then
            Call EXPVP (dEXO%exoid,'e',dEXO%Num_CellVariables,ierr)
            Call EXPVAN(dEXO%exoid,'e',dEXO%Num_CellVariables,dEXO%CellVariable(:)%Name,ierr)
         End If
         If (dEXO%Num_VertVariables > 0) Then
            Call EXPVP (dEXO%exoid,'n',dEXO%Num_VertVariables,ierr)
            Call EXPVAN(dEXO%exoid,'n',dEXO%Num_VertVariables,dEXO%VertVariable(:)%Name,ierr)
         End If
      End If
   End Subroutine EXOVariable_Write 
   
#undef __FUNCT__
#define __FUNCT__ "EXOVariable_Read"
   Subroutine EXOVariable_Read(dEXO)
      Type(EXO_Type)                                 :: dEXO

      PetscInt                                       :: ierr
      PetscInt                                       :: i
      Integer                                        :: EXO_MyRank
      
      Call MPI_COMM_RANK(dEXO%Comm,EXO_MyRank,ierr)

      If (EXO_MyRank == 0) Then
         If (dEXO%exoid == 0) Then
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_NULL,"[ERROR]: Exodus file not open before IO operations\n",ierr)
         End If

         Call EXGVP (dEXO%exoid,'g',dEXO%Num_GlobVariables,ierr)
         If (dEXO%Num_GlobVariables > 0) Then
            Allocate(dEXO%GlobVariable(dEXO%Num_GlobVariables))      
            Call EXGVAN(dEXO%exoid,'g',dEXO%Num_GlobVariables,dEXO%GlobVariable(:)%Name,ierr)
         End If
         
         Call EXGVP (dEXO%exoid,'e',dEXO%Num_CellVariables,ierr)
         If (dEXO%Num_CellVariables > 0) Then
            Allocate(dEXO%CellVariable(dEXO%Num_CellVariables))      
            Call EXGVAN(dEXO%exoid,'e',dEXO%Num_CellVariables,dEXO%CellVariable(:)%Name,ierr)
         End If

         Call EXGVP (dEXO%exoid,'n',dEXO%Num_VertVariables,ierr)
         If (dEXO%Num_VertVariables > 0) Then
            Allocate(dEXO%VertVariable(dEXO%Num_VertVariables))      
            Call EXGVAN(dEXO%exoid,'n',dEXO%Num_VertVariables,dEXO%VertVariable(:)%Name,ierr)
         End If
      End If

      Call MPI_BCast(dEXO%Num_GlobVariables,1,MPIU_INTEGER,0,dEXO%Comm,ierr)
      Call MPI_BCast(dEXO%Num_CellVariables,1,MPIU_INTEGER,0,dEXO%Comm,ierr)
      Call MPI_BCast(dEXO%Num_VertVariables,1,MPIU_INTEGER,0,dEXO%Comm,ierr)

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

      If (dEXO%Num_GlobVariables > 0) Then
         Call MPI_BCast(dEXO%GlobVariable(:)%Offset,dEXO%Num_GlobVariables,MPIU_INTEGER,0,dEXO%Comm,ierr)
      End If
      Do i = 1,dEXO%Num_GlobVariables
         Call MPI_BCast(dEXO%GlobVariable(i)%Name,MXSTLN,MPI_CHARACTER,0,dEXO%Comm,ierr)
      End Do

      If (dEXO%Num_CellVariables > 0) Then
         Call MPI_BCast(dEXO%CellVariable(:)%Offset,dEXO%Num_CellVariables,MPIU_INTEGER,0,dEXO%Comm,ierr)
      End If
      Do i = 1,dEXO%Num_CellVariables
         Call MPI_BCast(dEXO%CellVariable(i)%Name,MXSTLN,MPI_CHARACTER,0,dEXO%Comm,ierr)
      End Do

      If (dEXO%Num_VertVariables > 0) Then
         Call MPI_BCast(dEXO%VertVariable(:)%Offset,dEXO%Num_VertVariables,MPIU_INTEGER,0,dEXO%Comm,ierr)
      End If
      Do i = 1,dEXO%Num_VertVariables
         Call MPI_BCast(dEXO%VertVariable(i)%Name,MXSTLN,MPI_CHARACTER,0,dEXO%Comm,ierr)
      End Do
   End Subroutine EXOVariable_Read
   

#undef __FUNCT__
#define __FUNCT__ "Read_EXO_Result_Global"
   Subroutine Read_EXO_Result_Global(dEXO,dIdx,dTS,dRes)
      Type(EXO_Type),Intent(INOUT)                   :: dEXO
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      PetscReal                                      :: dRes
      
      PetscReal                                      :: MyRes
      PetscInt                                       :: ierr
      PetscReal,Dimension(:),Pointer                 :: Tmp_Res
      PetscInt                                       :: Num_Vars
      
      dRes = 0.0_Kr
      If ( ((dEXO%comm == PETSC_COMM_WORLD) .AND. (MEF90_MyRank == 0)) .OR. (dEXO%comm == PETSC_COMM_SELF) ) Then
         If (dEXO%exoid == 0) Then
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_NULL,"[ERROR]: Exodus file not open before IO operations\n",ierr)
         End If
         ! Get the number of global variables stored in the database    
         Call EXGVP(dEXO%exoid,'G',Num_Vars,ierr)
         Allocate(Tmp_Res(Num_Vars))
         
         ! Read All the global variables at time step TS and copy the right one
         ! into Res
         Call EXGGV(dEXO%exoid,dTS,Num_Vars,Tmp_Res,ierr)
         MyRes = Tmp_Res(dIdx)
         DeAllocate (Tmp_Res)
      End If
            
      !!! Broacast if dEXO%comm == PETSC_COMM_WORLD
      If (dEXO%comm == PETSC_COMM_WORLD) Then
         Call MPI_BCast(MyRes,1,MPIU_SCALAR,0,dEXO%comm,ierr)
      End If
      dRes = MyRes
   End Subroutine Read_EXO_Result_Global

!!! WRITE
#undef __FUNCT__
#define __FUNCT__ "Write_EXO_AllResult_Global"
   Subroutine Write_EXO_AllResult_Global(dEXO,dTS,dRes)
      Type(EXO_Type),Intent(INOUT)                   :: dEXO
      PetscInt                                       :: dTS
      PetscReal,Dimension(:),Pointer                 :: dRes
      
      PetscInt                                       :: ierr
      PetscInt                                       :: Num_Vars
      
      If ( ((dEXO%comm == PETSC_COMM_WORLD) .AND. (MEF90_MyRank == 0)) .OR. (dEXO%comm == PETSC_COMM_SELF) ) Then
         If (dEXO%exoid == 0) Then
            SETERRQ(dEXO%Comm,PETSC_ERR_ARG_NULL,"[ERROR]: Exodus file not open before IO operations\n",ierr)
         End If
         Call EXGVP(dEXO%exoid,'G',Num_Vars,ierr)
         If (Num_Vars /= Size(dRes)) Then
            SETERRQ(dEXO%Comm,PETSC_ERR_ARG_SIZ,'Write_EXO_AllResult_Global: The argument does not match the number of global variables in the mesh',ierr)
         End If
         Call EXPGV(dEXO%exoid,dTS,Num_Vars,dRes,ierr)
      End If
   End Subroutine Write_EXO_AllResult_Global

#undef __FUNCT__
#define __FUNCT__ "Write_EXO_Result_Global"
   Subroutine Write_EXO_Result_Global(dEXO,dIdx,dTS,dRes)
      Type(EXO_Type),Intent(INOUT)                   :: dEXO
      PetscInt                                       :: dIdx
      PetscInt                                       :: dTS
      PetscReal                                      :: dRes
      
      PetscInt                                       :: ierr
      PetscReal,Dimension(:),Pointer                 :: Tmp_Res
      PetscInt                                       :: Num_Vars
      
      If ( ((dEXO%comm == PETSC_COMM_WORLD) .AND. (MEF90_MyRank == 0)) .OR. (dEXO%comm == PETSC_COMM_SELF) ) Then
         If (dEXO%exoid == 0) Then
            SETERRQ(dEXO%Comm,PETSC_ERR_ARG_NULL,"[ERROR]: Exodus file not open before IO operations\n",ierr)
         End If
         ! Get the number of global variables stored in the database    
         Call EXGVP(dEXO%exoid,'G',Num_Vars,ierr)
         Allocate(Tmp_Res(Num_Vars))
         Tmp_Res = 0.0_Kr
         
         ! Read All the global variables at time step TS into Tmp_Res
         ! Modify Tmp_Res(Idx) and write everything back...
         Call EXGGV(dEXO%exoid,dTS,Num_Vars,Tmp_Res,ierr)
         Tmp_Res(dIdx) = dRes
         
         Call EXPGV(dEXO%exoid,dTS,Num_Vars,Tmp_Res,ierr)
         DeAllocate (Tmp_Res)
      End If
   End Subroutine Write_EXO_Result_Global

#undef __FUNCT__
#define __FUNCT__ "EXOView"
   Subroutine EXOView(dEXO,viewer)
      Type(EXO_Type)                               :: dEXO
      Type(PetscViewer)                            :: viewer
      
      PetscInt                                     :: i,j,iErr
      Character(len=MEF90_MXSTRLEN)                :: CharBuffer
   
      If (dEXO%comm == PETSC_COMM_WORLD) Then
         Write(CharBuffer,100) 'PETSC_COMM_WORLD'
      ElseIf (dEXO%comm == PETSC_COMM_SELF) Then
         Write(CharBuffer,100) 'PETSC_COMM_SELF'
      Else  
         Write(CharBuffer,105) 'unknown',dEXO%comm
      End If
      Call PetscViewerASCIIPrintf(viewer,CharBuffer,iErr); CHKERRQ(iErr)
      
      Write(CharBuffer,101) dEXO%exoid
      Call PetscViewerASCIIPrintf(viewer,CharBuffer,iErr); CHKERRQ(iErr)
      Write(CharBuffer,102) dEXO%filename
      Call PetscViewerASCIIPrintf(viewer,CharBuffer,iErr); CHKERRQ(iErr)
!      Write(CharBuffer,103) dEXO%title//' '
!      Call PetscViewerASCIIPrintf(viewer,CharBuffer,iErr); CHKERRQ(iErr)
!      Write(CharBuffer,'(A)') '\n'

      Write(CharBuffer,106) dEXO%Num_EBProperties
      Call PetscViewerASCIIPrintf(viewer,CharBuffer,iErr); CHKERRQ(iErr)
      Do i = 1,dEXO%Num_EBProperties
         Write(CharBuffer,109) dEXO%EBProperty(i)%Name,Size(dEXO%EBProperty(i)%Value)
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,iErr); CHKERRQ(iErr)
         Do j = 1,Size(dEXO%EBProperty(i)%Value)
            Write(CharBuffer,201) j,dEXO%EBProperty(i)%Value(j)
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,iErr); CHKERRQ(iErr)
         End Do
      End Do

      Write(CharBuffer,108) dEXO%Num_NSProperties
      Call PetscViewerASCIIPrintf(viewer,CharBuffer,iErr); CHKERRQ(iErr)
      Do i = 1,dEXO%Num_NSProperties
         Write(CharBuffer,109) dEXO%NSProperty(i)%Name,Size(dEXO%NSProperty(i)%Value)
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,iErr); CHKERRQ(iErr)
         Do j = 1,Size(dEXO%NSProperty(i)%Value)
            Write(CharBuffer,203) j,dEXO%NSProperty(i)%Value(j)
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,iErr); CHKERRQ(iErr)
         End Do
      End Do
!      Do i = 1,dEXO%Num_QA
!         Write(CharBuffer,104) i,dEXO%QA_rec(i,:)
!         Call PetscViewerASCIIPrintf(viewer,CharBuffer,iErr); CHKERRQ(iErr)
!      End Do
      Write(CharBuffer,110) dEXO%Num_GlobVariables
      Do i = 1,dEXO%Num_GlobVariables
         Write(CharBuffer,210) dEXO%GlobVariable(i)%Name,dEXO%GlobVariable(i)%offset
      End Do

      Write(CharBuffer,111) dEXO%Num_CellVariables
      Do i = 1,dEXO%Num_CellVariables
         Write(CharBuffer,210) dEXO%CellVariable(i)%Name,dEXO%CellVariable(i)%offset
      End Do

      Write(CharBuffer,112) dEXO%Num_VertVariables
      Do i = 1,dEXO%Num_VertVariables
         Write(CharBuffer,210) dEXO%VertVariable(i)%Name,dEXO%VertVariable(i)%offset
      End Do
    
      
 100 Format('Communicator:       ',A,      '\n')
 101 Format('exo ID:             ',I3,     '\n')
 102 Format('filename:           ',A,      '\n')
 105 Format('Communicator:       ',A,I3,  '\n')
 106 Format('Number of EB Properties: ',I3,'\n')
 108 Format('Number of NS Properties: ',I3,'\n')
 109 Format('   ',A,'num values:',I3,    '\n')
 110 Format('Global Variables: ',I3,       '\n')
 111 Format('Cell Variables:   ',I3,       '\n')
 112 Format('Vertex Variables: ',I3,       '\n')
 201 Format('   Element Block ',I3,' value ',I3,'\n')
 203 Format('   Node Set      ',I3,' value ',I3,'\n')
 210 Format(A,I3,'\n')
   End Subroutine EXOView

End Module m_MEF_EXO

