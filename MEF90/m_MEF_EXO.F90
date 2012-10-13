Module m_MEF_EXO
#include "finclude/petscdef.h"
   Use m_MEF_Parameters
   Use m_MEF_Ctx
   Use m_MEF_Utils
   Use m_MEF_Elements
   Use m_MEF_Sieve
   Use petsc
   IMPLICIT NONE

   !Integer,Parameter,Public                        :: exo_cpu_ws = 8
   !Integer,Parameter,Public                        :: exo_io_ws = 8
   PetscInt,Public                                 :: exo_ver

   
   Public :: MEF90Ctx_GetDMMeshEXO
   Public :: EXOGetCellSetElementType_Scal      
   Public :: EXOGetCellSetElementType_Vect      
   Public :: EXOGetCellSetElementType_Elast      
   Public :: Write_EXO_Case

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90Ctx_GetDMMeshEXO"
!!!
!!!  
!!!  MEF90Ctx_GetDMMeshEXO:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine MEF90Ctx_GetDMMeshEXO(MEF90Ctx,Mesh,ierr)
   Type(MEF90Ctx_Type),Intent(IN)                  :: MEF90Ctx
   Type(DM),Intent(OUT)                            :: Mesh
   PetscErrorCode,Intent(OUT)                      :: ierr

   Character(len=MEF90_MXSTRLEN)                   :: IOBuffer
   Character(len=80) ::filename
   Integer                                         :: cpu_ws,io_ws
   Real                                            :: exoVersion
   Integer                                         :: exoErr,exoUnit
   Type(DM)                                        :: tmpMesh
   
      Type(MEF90GlobalOptions_Type),pointer           :: GlobalOptions      

      Call PetscBagGetDataMEF90GlobalOptions(MEF90Ctx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)
   !!! Open input file
   If (MEF90_MyRank == 0) Then
      cpu_ws = 8
      io_ws = 8
      filename = Trim(MEF90Ctx%prefix)//'.gen'
      exoUnit = EXOPEN(filename,EXREAD,cpu_ws,io_ws,exoVersion,exoErr)
      If (exoerr < 0) Then
         Write(IOBuffer,*) '\n\nError opening EXO file ',trim(filename),'\n\n'
         Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr);
         SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,IOBuffer,ierr);
      EndIf
   End If
   If (MEF90_NumProcs == 1) Then
      Call DMMeshCreateExodusNG(PETSC_COMM_WORLD,exoUnit,Mesh,ierr);CHKERRQ(ierr)
   Else
      Call DMMeshCreateExodusNG(PETSC_COMM_WORLD,exoUnit,tmpMesh,ierr);CHKERRQ(ierr)   
      Call DMMeshDistribute(tmpMesh,PETSC_NULL_CHARACTER,mesh,ierr);CHKERRQ(ierr)
      Call DMDestroy(tmpMesh,ierr);CHKERRQ(ierr)
   End If
   Call EXCLOS(exoUnit,ierr)
End Subroutine MEF90Ctx_GetDMMeshEXO

#undef __FUNCT__
#define __FUNCT__ "MEF90Ctx_OpenEXO"
!!!
!!!  
!!!  MEF90Ctx_OpenEXO:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine MEF90Ctx_OpenEXO(MEF90Ctx,Mesh,ierr)
   Type(MEF90Ctx_Type),Intent(INOUT)               :: MEF90Ctx
   Type(DM), Intent(IN)                            :: Mesh
   PetscErrorCode,Intent(OUT)                      :: ierr

   Integer                                         :: exoUnitIN
   MPI_Comm                                        :: IOComm
   Integer                                         :: IORank
   Character(len=MEF90_MXSTRLEN)                   :: IOBuffer,filename
   Type(MEF90GlobalOptions_Type),pointer           :: GlobalOptions      
   Integer                                         :: cpu_ws,io_ws
   Real                                            :: exo_version
   Integer                                         :: exoerr
   

   Call PetscBagGetDataMEF90GlobalOptions(MEF90Ctx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)

   !!! Get name of output file
   Select Case (GlobalOptions%FileFormat)
   Case (MEF90FileFormat_EXOSplit)
      IOComm = PETSC_COMM_SELF
      Write(filename,100) trim(MEF90Ctx%prefix),MEF90Ctx%rank
   Case (MEF90FileFormat_EXOSingle)   
      IOComm = PETSC_COMM_WORLD
      Write(filename,101) trim(MEF90Ctx%prefix)
   End Select
100 Format(A,'-',I4.4,'.gen')
101 Format(A,'_out.gen')
   Call MPI_Comm_Rank(IOComm,IORank,ierr)

   !!! Open output file or create it and format it depending on loading type
   If (IORank == 0) Then
      cpu_ws = 8
      io_ws = 8
      MEF90Ctx%fileExoUnit = EXOPEN(filename,EXWRIT,cpu_ws,io_ws,exo_version,exoerr)
      If (exoerr < 0) Then
            If (GlobalOptions%verbose > 0) Then    
               Write(IOBuffer,*) 'EXO file ',trim(filename),' does not seem to exist. Creating it.\n'
               Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr);
            EndIf
            MEF90Ctx%fileExoUnit = EXCRE(trim(filename),EXCLOB,cpu_ws,io_ws,ierr)
            Select Case (GlobalOptions%FileFormat)
            Case (MEF90FileFormat_EXOSplit)
               Call DMmeshViewExodusSplit(mesh,MEF90Ctx%fileExoUnit,ierr)
            Case (MEF90FileFormat_EXOSingle)
               Write(filename,102) trim(MEF90Ctx%prefix)
               exoUnitIN = EXOPEN(filename,EXREAD,cpu_ws,io_ws,exo_version,exoerr)
               Call EXCOPY(exoUnitIN,MEF90Ctx%fileExoUnit,ierr)
               Call EXCLOS(exoUnitIN,ierr)
            End Select
      EndIf
   End If
102 Format(A,'.gen')
End Subroutine MEF90Ctx_OpenEXO

#undef __FUNCT__
#define __FUNCT__ "EXOGetCellSetElementType_Scal"
   Subroutine EXOGetCellSetElementType_Scal(comm,exoID,dim,ElemType,ierr)
      MPI_Comm,Intent(IN)                             :: Comm
      Integer,Intent(IN)                              :: exoID
      PetscInt,Intent(IN)                             :: dim
      Type(MEF90Element_Type),Dimension(:),Pointer    :: ElemType
      PetscErrorCode,Intent(OUT)                      :: ierr
      
      PetscInt                                        :: rank
      Character(len=MXSTLN)                           :: EXOElemType
      Integer                                         :: junk1,junk2,junk3,exoerr
      Real                                            :: dummyR
      Character(len=1)                                :: dummyS
      PetscInt                                        :: numSet,set
      PetscInt,Dimension(:),Pointer                   :: setID
      
      Call MPI_Comm_rank(comm,rank,ierr)
      If (rank == 0) Then
         Call EXINQ(exoID,EXELBL,numSet,dummyR,dummyS,exoerr)
      End If
      Call MPI_Bcast(numSet,1,MPIU_INTEGER,0,comm,ierr)
      Allocate(ElemType(numSet))
      Allocate(setID(numSet))
      If (rank == 0) Then       
         Call EXGEBI(exoID,setID,exoerr)
      End If
      Call MPI_Bcast(setID,numSet,MPIU_INTEGER,0,comm,ierr)
      Do set = 1, numSet
         If (rank == 0) Then
            Call EXGELB(exoID,setID(set),EXOElemType,junk1,junk2,junk3,exoerr)
         EndIf
         Call MPI_Bcast(EXOElemType,MXSTLN,MPI_CHAR,0,comm,ierr)
         Call EXO2MEF90ElementType_Scal(EXOElemType,dim,ElemType(set),ierr)
      End Do
      DeAllocate(setID)
   End Subroutine EXOGetCellSetElementType_Scal
      
#undef __FUNCT__
#define __FUNCT__ "EXOGetCellSetElementType_Vect"
   Subroutine EXOGetCellSetElementType_Vect(comm,exoID,dim,ElemType,ierr)
      MPI_Comm,Intent(IN)                             :: Comm
      Integer,Intent(IN)                              :: exoID
      PetscInt,Intent(IN)                             :: dim
      Type(MEF90Element_Type),Dimension(:),Pointer    :: ElemType
      PetscErrorCode,Intent(OUT)                      :: ierr
      
      PetscInt                                        :: rank
      Character(len=MXSTLN)                           :: EXOElemType
      Integer                                         :: junk1,junk2,junk3,exoerr
      Real                                            :: dummyR
      Character(len=1)                                :: dummyS
      PetscInt                                        :: numSet,set
      PetscInt,Dimension(:),Pointer                   :: setID
      
      Call MPI_Comm_rank(comm,rank,ierr)
      If (rank == 0) Then
         Call EXINQ(exoID,EXELBL,numSet,dummyR,dummyS,exoerr)
      End If
      Call MPI_Bcast(numSet,1,MPIU_INTEGER,0,comm,ierr)
      Allocate(ElemType(numSet))
      Allocate(setID(numSet))
      If (rank == 0) Then       
         Call EXGEBI(exoID,setID,exoerr)
      End If
      Call MPI_Bcast(setID,numSet,MPIU_INTEGER,0,comm,ierr)
      Do set = 1, numSet
         If (rank == 0) Then
            Call EXGELB(exoID,setID(set),EXOElemType,junk1,junk2,junk3,exoerr)
         EndIf
         Call MPI_Bcast(EXOElemType,MXSTLN,MPI_CHAR,0,comm,ierr)
         Call EXO2MEF90ElementType_Vect(EXOElemType,dim,ElemType(set),ierr)
      End Do
      DeAllocate(setID)
   End Subroutine EXOGetCellSetElementType_Vect
      
#undef __FUNCT__
#define __FUNCT__ "EXOGetCellSetElementType_Elast"
   Subroutine EXOGetCellSetElementType_Elast(comm,exoID,dim,ElemType,ierr)
      MPI_Comm,Intent(IN)                             :: Comm
      Integer,Intent(IN)                              :: exoID
      PetscInt,Intent(IN)                             :: dim
      Type(MEF90Element_Type),Dimension(:),Pointer    :: ElemType
      PetscErrorCode,Intent(OUT)                      :: ierr
      
      PetscInt                                        :: rank
      Character(len=MXSTLN)                           :: EXOElemType
      Integer                                         :: junk1,junk2,junk3,exoerr
      Real                                            :: dummyR
      Character(len=1)                                :: dummyS
      PetscInt                                        :: numSet,set
      PetscInt,Dimension(:),Pointer                   :: setID
      
      Call MPI_Comm_rank(comm,rank,ierr)
      If (rank == 0) Then
         Call EXINQ(exoID,EXELBL,numSet,dummyR,dummyS,exoerr)
      End If
      Call MPI_Bcast(numSet,1,MPIU_INTEGER,0,comm,ierr)
      Allocate(ElemType(numSet))
      Allocate(setID(numSet))
      If (rank == 0) Then       
         Call EXGEBI(exoID,setID,exoerr)
      End If
      Call MPI_Bcast(setID,numSet,MPIU_INTEGER,0,comm,ierr)
      Do set = 1, numSet
         If (rank == 0) Then
            Call EXGELB(exoID,setID(set),EXOElemType,junk1,junk2,junk3,exoerr)
         EndIf
         Call MPI_Bcast(EXOElemType,MXSTLN,MPI_CHAR,0,comm,ierr)
         Call EXO2MEF90ElementType_Elast(EXOElemType,dim,ElemType(set),ierr)
      End Do
      DeAllocate(setID)
   End Subroutine EXOGetCellSetElementType_Elast
      
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
   

End Module m_MEF_EXO

