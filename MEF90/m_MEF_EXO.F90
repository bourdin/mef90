Module m_MEF_EXO
#include "finclude/petscdef.h"
   Use m_MEF_Parameters
   Use m_MEF_Utils
   Use m_MEF_Elements
   Use petsc
   IMPLICIT NONE

   !Integer,Parameter,Public                        :: exo_cpu_ws = 8
   !Integer,Parameter,Public                        :: exo_io_ws = 8
   PetscInt,Public                                 :: exo_ver

   
   Public :: EXOGetCellSetElementType_Scal      
   Public :: EXOGetCellSetElementType_Vect      
   Public :: EXOGetCellSetElementType_Elast      
   Public :: Write_EXO_Case

Contains
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

