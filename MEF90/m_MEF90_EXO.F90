Module m_MEF90_EXO
#include "petsc/finclude/petsc.h"
   Use m_MEF90_Parameters
   Use m_MEF90_Ctx
   Use m_MEF90_Utils
   Use m_MEF90_Elements
   Use petsc
   IMPLICIT NONE
#include "../mef90version.h"

   Private 
   PetscInt,Public                                 :: exo_ver

   
!   Public :: MEF90CtxOpenEXO
   Public :: MEF90CtxCloseEXO
   Public :: MEF90EXOFormat
   Public :: EXOGetCellSetElementType_Scal      
   Public :: EXOGetCellSetElementType_Vect      
   Public :: EXOGetCellSetElementType_Elast      
   Public :: EXOWriteCase

Contains
<<<<<<< HEAD
#undef __FUNCT__
<<<<<<< HEAD
#define __FUNCT__ "MEF90CtxGetDMMeshEXO"
!!!
!!!  
!!!  MEF90CtxGetDMMeshEXO:
!!!  
!!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90CtxGetDMMeshEXO(MEF90Ctx,Mesh,ierr)
      Type(MEF90Ctx_Type),Intent(IN)                  :: MEF90Ctx
      Type(DM),Intent(OUT)                            :: Mesh
      PetscErrorCode,Intent(OUT)                      :: ierr
   
      Character(len=MEF90_MXSTRLEN)                   :: IOBuffer
      Character(len=MEF90_MXSTRLEN)                   :: filename
      Integer                                         :: cpu_ws,io_ws
      Real                                            :: exoVersion
      Integer                                         :: exoErr,exoUnit
      Type(DM)                                        :: tmpMesh
      Character(len=MXSTLN),Dimension(:),Pointer      :: cellSetName
      
      Type(MEF90CtxGlobalOptions_Type),pointer        :: GlobalOptions      
      PetscSizeT                                      :: sizeofPetscReal
   

   
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)
      !!! Open input file
      If (MEF90_MyRank == 0) Then
         Call PetscDataTypeGetSize(PETSC_REAL,sizeofPetscReal,ierr)
         cpu_ws = int(sizeofPetscReal,kind(cpu_ws))
         io_ws = 0
         filename = Trim(MEF90Ctx%geometryFile)
         exoUnit = EXOPEN(filename,EXREAD,cpu_ws,io_ws,exoVersion,exoErr)
         Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr);

         If (exoerr < 0) Then
            Write(IOBuffer,*) '\n\nError opening EXO file ',trim(filename),'\n\n'
            Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr);
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,IOBuffer,ierr);
         EndIf
      End If
      If (MEF90_NumProcs == 1) Then
         Call DMMeshCreateExodus(PETSC_COMM_WORLD,exoUnit,Mesh,ierr);CHKERRQ(ierr)
      Else
         Call DMMeshCreateExodus(PETSC_COMM_WORLD,exoUnit,tmpMesh,ierr);CHKERRQ(ierr)   
         Call DMMeshDistribute(tmpMesh,PETSC_NULL_CHARACTER,mesh,ierr);CHKERRQ(ierr)
         Call DMDestroy(tmpMesh,ierr);CHKERRQ(ierr)
      End If
      If (MEF90_MyRank == 0) Then
         Call EXCLOS(exoUnit,exoErr)
      End If
   End Subroutine MEF90CtxGetDMMeshEXO


#undef __FUNCT__
#define __FUNCT__ "MEF90CtxOpenEXO"
=======
#define __FUNCT__ "MEF90EXOCtxOpenEXO"
>>>>>>> f87946a (Updated m_MEF90_Materials, started work on elements)
=======
!#undef __FUNCT__
!#define __FUNCT__ "MEF90CtxOpenEXO"
>>>>>>> 7aedf22 (Updated geometry initialization)
!!!
!!!  
!!!  MEF90CtxOpenEXO:
!!!  
!!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!
<<<<<<< HEAD
   Subroutine MEF90CtxOpenEXO(MEF90Ctx,Mesh,ierr)
      Type(MEF90Ctx_Type),Intent(INOUT)               :: MEF90Ctx
      Type(tDM), Intent(IN)                           :: Mesh
      PetscErrorCode,Intent(OUT)                      :: ierr
   
      Integer                                         :: exoUnitIN
      MPI_Comm                                        :: IOComm
      Integer                                         :: IORank
      Character(len=MEF90_MXSTRLEN)                   :: IOBuffer,filename
      Integer,parameter                               :: num_QA_rec=1
      Character(len=MXSTLN)                           :: QA_rec(4)
      Character(len=MXSTLN)                           :: date
      Character(len=MXSTLN)                           :: time

      Type(MEF90CtxGlobalOptions_Type),pointer        :: GlobalOptions      
      Integer                                         :: cpu_ws,io_ws
      Real                                            :: exo_version
      Integer                                         :: exoerr
      Logical                                         :: exoExists
      PetscSizeT                                      :: sizeofPetscReal
   
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)
      !!! Get name of output file
      Select Case (GlobalOptions%FileFormat)
      Case (MEF90FileFormat_EXOSplit)
         IOComm = PETSC_COMM_SELF
         Write(filename,100) trim(MEF90FilePrefix(MEF90Ctx%resultfile)),MEF90Ctx%rank,trim(MEF90FileExtension(MEF90Ctx%resultfile))
      Case (MEF90FileFormat_EXOSingle)   
         IOComm = MEF90Ctx%comm
         filename = MEF90Ctx%resultFile
      End Select
   100 Format(A,'-',I4.4,'.',A)
      Call MPI_Comm_Rank(IOComm,IORank,ierr)
   
      !!! Open output file or create it and format it depending on loading type
      If (IORank == 0) Then
         Call PetscDataTypeGetSize(PETSC_REAL,sizeofPetscReal,ierr)
         cpu_ws = int(sizeofPetscReal,kind(cpu_ws))
         Inquire(file=filename,exist=exoExists)
         If (.NOT. exoExists) Then
            If (GlobalOptions%verbose > 0) Then    
               Write(IOBuffer,*) 'EXO file ',trim(filename),' does not seem to exist. Creating it.\n'
               Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr);
            EndIf
            io_ws = int(sizeofPetscReal,kind(cpu_ws))
            MEF90Ctx%fileExoUnit = EXCRE(trim(filename),EXCLOB,cpu_ws,io_ws,ierr)
            Select Case (GlobalOptions%FileFormat)
            Case (MEF90FileFormat_EXOSplit)
               Call DMmeshViewExodusSplit(mesh,MEF90Ctx%fileExoUnit,ierr)
            Case (MEF90FileFormat_EXOSingle)
               exoUnitIN = EXOPEN(MEF90Ctx%geometryfile,EXREAD,cpu_ws,io_ws,exo_version,exoerr)
               Call EXCOPY(exoUnitIN,MEF90Ctx%fileExoUnit,exoErr)
               Call EXCLOS(exoUnitIN,exoErr)
            End Select
         Else
            io_ws  = 0
            MEF90Ctx%fileExoUnit = EXOPEN(filename,EXWRIT,cpu_ws,io_ws,exo_version,exoerr)
            QA_rec(1) = "mef90"
            QA_rec(2) = MEF90_GITVER
            Call date_and_time(DATE=date,TIME=time)
            write(QA_rec(3),"(a)") date
            write(QA_rec(4),"(a,':',a,':',a)") time(1:2),time(3:4),time(5:6)
            call expqa (MEF90Ctx%fileExoUnit , num_QA_rec, QA_rec, ierr)
         EndIf
      End If
   End Subroutine MEF90CtxOpenEXO
=======
!   Subroutine MEF90CtxOpenEXO(MEF90Ctx,Mesh,ierr)
!      Type(MEF90Ctx_Type),Intent(INOUT)               :: MEF90Ctx
!      Type(tDM), Intent(IN)                           :: Mesh
!      PetscErrorCode,Intent(OUT)                      :: ierr
!   
!      Integer                                         :: exoUnitIN
!      MPI_Comm                                        :: IOComm
!      Integer                                         :: IORank
!      Character(len=MEF90_MXSTRLEN)                   :: IOBuffer,filename
!      Type(MEF90CtxGlobalOptions_Type),pointer        :: GlobalOptions      
!      Integer                                         :: cpu_ws,io_ws
!      Real                                            :: exo_version
!      Integer                                         :: exoerr
!      Logical                                         :: exoExists
!      
!   
!      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)
!   
!      !!! Get name of output file
!      Select Case (GlobalOptions%FileFormat)
!      Case (MEF90FileFormat_EXOSplit)
!         IOComm = PETSC_COMM_SELF
!         Write(filename,100) trim(MEF90Ctx%prefix),MEF90Ctx%rank
!      Case (MEF90FileFormat_EXOSingle)   
!         IOComm = PETSC_COMM_WORLD
!         Write(filename,101) trim(MEF90Ctx%prefix)
!      End Select
!   100 Format(A,'-',I4.4,'.gen')
!   101 Format(A,'_out.gen')
!      Call MPI_Comm_Rank(IOComm,IORank,ierr)
!   
!      !!! Open output file or create it and format it depending on loading type
!      If (IORank == 0) Then
!         Inquire(file=filename,exist=exoExists)
!         cpu_ws = 8
!         io_ws = 8
!         If (.NOT. exoExists) Then
!            If (GlobalOptions%verbose > 0) Then    
!               Write(IOBuffer,*) 'EXO file ',trim(filename),' does not seem to exist. Creating it.\n'
!               Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr);
!            EndIf
!            MEF90Ctx%fileExoUnit = EXCRE(trim(filename),EXCLOB,cpu_ws,io_ws,ierr)
!            Select Case (GlobalOptions%FileFormat)
!            Case (MEF90FileFormat_EXOSplit)
!               Call DMmeshViewExodusSplit(mesh,MEF90Ctx%fileExoUnit,ierr)
!            Case (MEF90FileFormat_EXOSingle)
!               Write(filename,102) trim(MEF90Ctx%prefix)
!               exoUnitIN = EXOPEN(filename,EXREAD,cpu_ws,io_ws,exo_version,exoerr)
!               Call EXCOPY(exoUnitIN,MEF90Ctx%fileExoUnit,exoErr)
!               Call EXCLOS(exoUnitIN,exoErr)
!            End Select
!         Else
!            MEF90Ctx%fileExoUnit = EXOPEN(filename,EXWRIT,cpu_ws,io_ws,exo_version,exoerr)
!         EndIf
!      End If
!   102 Format(A,'.gen')
!   End Subroutine MEF90CtxOpenEXO
>>>>>>> 7aedf22 (Updated geometry initialization)


#undef __FUNCT__
#define __FUNCT__ "MEF90CtxCloseEXO"
!!!
!!!  
!!!  MEF90CtxCloseEXO:
!!!  
!!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90CtxCloseEXO(MEF90Ctx,ierr)
      Type(MEF90Ctx_Type),Intent(INOUT)               :: MEF90Ctx
      PetscErrorCode,Intent(OUT)                      :: ierr
   
      Integer                                         :: exoErr
      MPI_Comm                                        :: IOComm
      Integer                                         :: IORank
      Character(len=MEF90_MXSTRLEN)                   :: filename
      Type(MEF90CtxGlobalOptions_Type),pointer        :: GlobalOptions      
<<<<<<< HEAD
      Real                                            :: exo_version
=======
>>>>>>> f87946a (Updated m_MEF90_Materials, started work on elements)
      
   
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)
      Select Case (GlobalOptions%FileFormat)
      Case (MEF90FileFormat_EXOSplit)
         IOComm = PETSC_COMM_SELF
      Case (MEF90FileFormat_EXOSingle)   
         IOComm = PETSC_COMM_WORLD
      End Select
      Call MPI_Comm_Rank(IOComm,IORank,ierr)   

      If (IORank == 0) Then
         Call EXCLOS(MEF90Ctx%fileExoUnit,exoErr)
         MEF90Ctx%fileExoUnit = 0
      End If
   End Subroutine MEF90CtxCloseEXO

<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 7aedf22 (Updated geometry initialization)
#undef __FUNCT__
#define __FUNCT__ "EXOGetCellSetElementType_Scal"
   Subroutine EXOGetCellSetElementType_Scal(MEF90Ctx,elemType,ierr)
      Type(MEF90Ctx_Type),Intent(IN)                  :: MEF90Ctx
      Type(MEF90Element_Type),Dimension(:),Pointer    :: elemType
      PetscErrorCode,Intent(OUT)                      :: ierr
      
      Character(len=MEF90_MXSTRLEN)                   :: filename,IOBuffer
      Character(len=MXSTLN)                           :: EXOelemType
      Integer                                         :: cpu_ws,io_ws,exoID
      Real                                            :: exoVersion
<<<<<<< HEAD
      Integer                                         :: junk1,junk2,junk3,junk4,exoerr
      Real                                            :: dummyR
=======
      Integer                                         :: junk1,junk2,junk3,exoerr
>>>>>>> 7aedf22 (Updated geometry initialization)
      Character(len=MXLNLN)                           :: dummyS
      Integer                                         :: numSet,set,numDim
      PetscInt,Dimension(:),Pointer                   :: setID
      Type(MEF90CtxGlobalOptions_Type),pointer        :: GlobalOptions      
<<<<<<< HEAD
      PetscSizeT                                      :: sizeofPetscReal

   
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)
      If (MEF90Ctx%rank == 0) Then
         Call PetscDataTypeGetSize(PETSC_REAL,sizeofPetscReal,ierr)
         cpu_ws = int(sizeofPetscReal,kind(cpu_ws))      
         io_ws = 0
         filename = MEF90Ctx%geometryFile
=======
   
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)
      If (MEF90Ctx%rank == 0) Then
         cpu_ws = 8
         io_ws = 8
         filename = Trim(MEF90Ctx%prefix)//'.gen'
>>>>>>> 7aedf22 (Updated geometry initialization)
         exoID = EXOPEN(filename,EXREAD,cpu_ws,io_ws,exoVersion,exoErr)
         If (exoerr < 0) Then
            Write(IOBuffer,*) '\n\nError opening EXO file ',trim(filename),'\n\n'
            Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr);
            STOP
<<<<<<< HEAD
            !SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,IOBuffer,ierr);
=======
            !SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,IOBuffer);
>>>>>>> 7aedf22 (Updated geometry initialization)
         EndIf
         Call EXGINI(exoid,dummyS,numDim,junk1,junk2,numSet,junk3,junk3,exoerr)
      End If
      Call MPI_Bcast(numSet,1,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)
      Call MPI_Bcast(numDim,1,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)

      Allocate(elemType(numSet))
      Allocate(setID(numSet))
      If (MEF90Ctx%rank == 0) Then       
         Call EXGEBI(exoID,setID,exoerr)
      End If
      Call MPI_Bcast(setID,numSet,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)
      Do set = 1, numSet
         If (MEF90Ctx%rank == 0) Then
            Call EXGELB(exoID,setID(set),EXOelemType,junk1,junk2,junk3,exoerr)
         EndIf
         Call MPI_Bcast(EXOelemType,MXSTLN,MPI_CHAR,0,MEF90Ctx%comm,ierr)
         Call EXO2MEF90ElementType_Scal(EXOelemType,numDim,elemType(set),ierr)
      End Do
      DeAllocate(setID)
      If (MEF90Ctx%rank == 0) Then
         Call EXCLOS(exoid,exoErr)
      End If
   End Subroutine EXOGetCellSetElementType_Scal
      
#undef __FUNCT__
#define __FUNCT__ "EXOGetCellSetElementType_Vect"
   Subroutine EXOGetCellSetElementType_Vect(MEF90Ctx,elemType,ierr)
      Type(MEF90Ctx_Type),Intent(IN)                  :: MEF90Ctx
      Type(MEF90Element_Type),Dimension(:),Pointer    :: elemType
      PetscErrorCode,Intent(OUT)                      :: ierr
      
      Character(len=MEF90_MXSTRLEN)                   :: filename,IOBuffer
      Character(len=MXSTLN)                           :: EXOelemType
      Integer                                         :: cpu_ws,io_ws,exoID
      Real                                            :: exoVersion
<<<<<<< HEAD
      Integer                                         :: junk1,junk2,junk3,junk4,exoerr
      Real                                            :: dummyR
      Character(len=MXLNLN)                           :: dummyS
      Integer                                         :: numSet,set,numDim
      PetscInt,Dimension(:),Pointer                   :: setID
      Type(MEF90CtxGlobalOptions_Type),pointer        :: GlobalOptions
      PetscSizeT                                      :: sizeofPetscReal      
   
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)
      If (MEF90Ctx%rank == 0) Then
         Call PetscDataTypeGetSize(PETSC_REAL,sizeofPetscReal,ierr)
         cpu_ws = int(sizeofPetscReal,kind(cpu_ws))      
         io_ws = 0
         filename = MEF90Ctx%geometryFile
=======
      Integer                                         :: junk1,junk2,junk3,exoerr
      Character(len=MXLNLN)                           :: dummyS
      Integer                                         :: numSet,set,numDim
      PetscInt,Dimension(:),Pointer                   :: setID
      Type(MEF90CtxGlobalOptions_Type),pointer        :: GlobalOptions      
   
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)
      If (MEF90Ctx%rank == 0) Then
         cpu_ws = 8
         io_ws = 8
         filename = Trim(MEF90Ctx%prefix)//'.gen'
>>>>>>> 7aedf22 (Updated geometry initialization)
         exoID = EXOPEN(filename,EXREAD,cpu_ws,io_ws,exoVersion,exoErr)
         If (exoerr < 0) Then
            Write(IOBuffer,*) '\n\nError opening EXO file ',trim(filename),'\n\n'
            Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr);
            STOP
<<<<<<< HEAD
            !SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,IOBuffer,ierr);
=======
            !SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,IOBuffer);
>>>>>>> 7aedf22 (Updated geometry initialization)
         EndIf
         Call EXGINI(exoid,dummyS,numDim,junk1,junk2,numSet,junk3,junk3,exoerr)
      End If
      Call MPI_Bcast(numSet,1,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)
      Call MPI_Bcast(numDim,1,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)

      Allocate(elemType(numSet))
      Allocate(setID(numSet))
      If (MEF90Ctx%rank == 0) Then       
         Call EXGEBI(exoID,setID,exoerr)
      End If
      Call MPI_Bcast(setID,numSet,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)
      Do set = 1, numSet
         If (MEF90Ctx%rank == 0) Then
            Call EXGELB(exoID,setID(set),EXOelemType,junk1,junk2,junk3,exoerr)
         EndIf
         Call MPI_Bcast(EXOelemType,MXSTLN,MPI_CHAR,0,MEF90Ctx%comm,ierr)
         Call EXO2MEF90ElementType_Vect(EXOelemType,numDim,elemType(set),ierr)
      End Do
      DeAllocate(setID)
      If (MEF90Ctx%rank == 0) Then
         Call EXCLOS(exoid,exoErr)
      End If
   End Subroutine EXOGetCellSetElementType_Vect
      
#undef __FUNCT__
#define __FUNCT__ "EXOGetCellSetElementType_Elast"
   Subroutine EXOGetCellSetElementType_Elast(MEF90Ctx,elemType,ierr)
      Type(MEF90Ctx_Type),Intent(IN)                  :: MEF90Ctx
      Type(MEF90Element_Type),Dimension(:),Pointer    :: elemType
      PetscErrorCode,Intent(OUT)                      :: ierr
      
      Character(len=MEF90_MXSTRLEN)                   :: filename,IOBuffer
      Character(len=MXSTLN)                           :: EXOelemType
      Integer                                         :: cpu_ws,io_ws,exoID
      Real                                            :: exoVersion
<<<<<<< HEAD
      Integer                                         :: junk1,junk2,junk3,junk4,exoerr
      Real                                            :: dummyR
=======
      Integer                                         :: junk1,junk2,junk3,exoerr
>>>>>>> 7aedf22 (Updated geometry initialization)
      Character(len=MXLNLN)                           :: dummyS
      Integer                                         :: numSet,set,numDim
      PetscInt,Dimension(:),Pointer                   :: setID
      Type(MEF90CtxGlobalOptions_Type),pointer        :: GlobalOptions      
<<<<<<< HEAD
      PetscSizeT                                      :: sizeofPetscReal

      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)
      If (MEF90Ctx%rank == 0) Then
         Call PetscDataTypeGetSize(PETSC_REAL,sizeofPetscReal,ierr)
         cpu_ws = int(sizeofPetscReal,kind(cpu_ws))
         io_ws = 0
         filename = MEF90Ctx%geometryFile
=======
   
      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)
      If (MEF90Ctx%rank == 0) Then
         cpu_ws = 8
         io_ws = 8
         filename = Trim(MEF90Ctx%prefix)//'.gen'
>>>>>>> 7aedf22 (Updated geometry initialization)
         exoID = EXOPEN(filename,EXREAD,cpu_ws,io_ws,exoVersion,exoErr)
         If (exoerr < 0) Then
            Write(IOBuffer,*) '\n\nError opening EXO file ',trim(filename),'\n\n'
            Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr);
            STOP
<<<<<<< HEAD
            !SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,IOBuffer,ierr);
=======
            !SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,IOBuffer);
>>>>>>> 7aedf22 (Updated geometry initialization)
         EndIf
         Call EXGINI(exoid,dummyS,numDim,junk1,junk2,numSet,junk3,junk3,exoerr)
      End If
      Call MPI_Bcast(numSet,1,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)
      Call MPI_Bcast(numDim,1,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)

      Allocate(elemType(numSet))
      Allocate(setID(numSet))
      If (MEF90Ctx%rank == 0) Then       
         Call EXGEBI(exoID,setID,exoerr)
      End If
      Call MPI_Bcast(setID,numSet,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)
      Do set = 1, numSet
         If (MEF90Ctx%rank == 0) Then
            Call EXGELB(exoID,setID(set),EXOelemType,junk1,junk2,junk3,exoerr)
         EndIf
         Call MPI_Bcast(EXOelemType,MXSTLN,MPI_CHAR,0,MEF90Ctx%comm,ierr)
         Call EXO2MEF90ElementType_Elast(EXOelemType,numDim,elemType(set),ierr)
      End Do
      DeAllocate(setID)
      If (MEF90Ctx%rank == 0) Then
         Call EXCLOS(exoid,exoErr)
      End If
   End Subroutine EXOGetCellSetElementType_Elast
<<<<<<< HEAD
=======
!#undef __FUNCT__
!#define __FUNCT__ "EXOGetCellSetElementType_Scal"
!   Subroutine EXOGetCellSetElementType_Scal(MEF90Ctx,elemType,ierr)
!      Type(MEF90Ctx_Type),Intent(IN)                  :: MEF90Ctx
!      Type(MEF90Element_Type),Dimension(:),Pointer    :: elemType
!      PetscErrorCode,Intent(OUT)                      :: ierr
!      
!      Character(len=MEF90_MXSTRLEN)                   :: filename,IOBuffer
!      Character(len=MXSTLN)                           :: EXOelemType
!      Integer                                         :: cpu_ws,io_ws,exoID
!      Real                                            :: exoVersion
!      Integer                                         :: junk1,junk2,junk3,junk4,exoerr
!      Real                                            :: dummyR
!      Character(len=MXLNLN)                           :: dummyS
!      Integer                                         :: numSet,set,numDim
!      PetscInt,Dimension(:),Pointer                   :: setID
!      Type(MEF90CtxGlobalOptions_Type),pointer        :: GlobalOptions      
!   
!      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)
!      If (MEF90Ctx%rank == 0) Then
!         cpu_ws = 8
!         io_ws = 8
!         filename = Trim(MEF90Ctx%prefix)//'.gen'
!         exoID = EXOPEN(filename,EXREAD,cpu_ws,io_ws,exoVersion,exoErr)
!         If (exoerr < 0) Then
!            Write(IOBuffer,*) '\n\nError opening EXO file ',trim(filename),'\n\n'
!            Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr);
!            STOP
!            !SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,IOBuffer);
!         EndIf
!         Call EXGINI(exoid,dummyS,numDim,junk1,junk2,numSet,junk3,junk3,exoerr)
!      End If
!      Call MPI_Bcast(numSet,1,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)
!      Call MPI_Bcast(numDim,1,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)
!
!      Allocate(elemType(numSet))
!      Allocate(setID(numSet))
!      If (MEF90Ctx%rank == 0) Then       
!         Call EXGEBI(exoID,setID,exoerr)
!      End If
!      Call MPI_Bcast(setID,numSet,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)
!      Do set = 1, numSet
!         If (MEF90Ctx%rank == 0) Then
!            Call EXGELB(exoID,setID(set),EXOelemType,junk1,junk2,junk3,exoerr)
!         EndIf
!         Call MPI_Bcast(EXOelemType,MXSTLN,MPI_CHAR,0,MEF90Ctx%comm,ierr)
!         Call EXO2MEF90ElementType_Scal(EXOelemType,numDim,elemType(set),ierr)
!      End Do
!      DeAllocate(setID)
!      If (MEF90Ctx%rank == 0) Then
!         Call EXCLOS(exoid,exoErr)
!      End If
!   End Subroutine EXOGetCellSetElementType_Scal
!      
!#undef __FUNCT__
!#define __FUNCT__ "EXOGetCellSetElementType_Vect"
!   Subroutine EXOGetCellSetElementType_Vect(MEF90Ctx,elemType,ierr)
!      Type(MEF90Ctx_Type),Intent(IN)                  :: MEF90Ctx
!      Type(MEF90Element_Type),Dimension(:),Pointer    :: elemType
!      PetscErrorCode,Intent(OUT)                      :: ierr
!      
!      Character(len=MEF90_MXSTRLEN)                   :: filename,IOBuffer
!      Character(len=MXSTLN)                           :: EXOelemType
!      Integer                                         :: cpu_ws,io_ws,exoID
!      Real                                            :: exoVersion
!      Integer                                         :: junk1,junk2,junk3,junk4,exoerr
!      Real                                            :: dummyR
!      Character(len=MXLNLN)                           :: dummyS
!      Integer                                         :: numSet,set,numDim
!      PetscInt,Dimension(:),Pointer                   :: setID
!      Type(MEF90CtxGlobalOptions_Type),pointer        :: GlobalOptions      
!   
!      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)
!      If (MEF90Ctx%rank == 0) Then
!         cpu_ws = 8
!         io_ws = 8
!         filename = Trim(MEF90Ctx%prefix)//'.gen'
!         exoID = EXOPEN(filename,EXREAD,cpu_ws,io_ws,exoVersion,exoErr)
!         If (exoerr < 0) Then
!            Write(IOBuffer,*) '\n\nError opening EXO file ',trim(filename),'\n\n'
!            Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr);
!            STOP
!            !SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,IOBuffer);
!         EndIf
!         Call EXGINI(exoid,dummyS,numDim,junk1,junk2,numSet,junk3,junk3,exoerr)
!      End If
!      Call MPI_Bcast(numSet,1,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)
!      Call MPI_Bcast(numDim,1,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)
!
!      Allocate(elemType(numSet))
!      Allocate(setID(numSet))
!      If (MEF90Ctx%rank == 0) Then       
!         Call EXGEBI(exoID,setID,exoerr)
!      End If
!      Call MPI_Bcast(setID,numSet,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)
!      Do set = 1, numSet
!         If (MEF90Ctx%rank == 0) Then
!            Call EXGELB(exoID,setID(set),EXOelemType,junk1,junk2,junk3,exoerr)
!         EndIf
!         Call MPI_Bcast(EXOelemType,MXSTLN,MPI_CHAR,0,MEF90Ctx%comm,ierr)
!         Call EXO2MEF90ElementType_Vect(EXOelemType,numDim,elemType(set),ierr)
!      End Do
!      DeAllocate(setID)
!      If (MEF90Ctx%rank == 0) Then
!         Call EXCLOS(exoid,exoErr)
!      End If
!   End Subroutine EXOGetCellSetElementType_Vect
!      
!#undef __FUNCT__
!#define __FUNCT__ "EXOGetCellSetElementType_Elast"
!   Subroutine EXOGetCellSetElementType_Elast(MEF90Ctx,elemType,ierr)
!      Type(MEF90Ctx_Type),Intent(IN)                  :: MEF90Ctx
!      Type(MEF90Element_Type),Dimension(:),Pointer    :: elemType
!      PetscErrorCode,Intent(OUT)                      :: ierr
!      
!      Character(len=MEF90_MXSTRLEN)                   :: filename,IOBuffer
!      Character(len=MXSTLN)                           :: EXOelemType
!      Integer                                         :: cpu_ws,io_ws,exoID
!      Real                                            :: exoVersion
!      Integer                                         :: junk1,junk2,junk3,junk4,exoerr
!      Real                                            :: dummyR
!      Character(len=MXLNLN)                           :: dummyS
!      Integer                                         :: numSet,set,numDim
!      PetscInt,Dimension(:),Pointer                   :: setID
!      Type(MEF90CtxGlobalOptions_Type),pointer        :: GlobalOptions      
!   
!      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)
!      If (MEF90Ctx%rank == 0) Then
!         cpu_ws = 8
!         io_ws = 8
!         filename = Trim(MEF90Ctx%prefix)//'.gen'
!         exoID = EXOPEN(filename,EXREAD,cpu_ws,io_ws,exoVersion,exoErr)
!         If (exoerr < 0) Then
!            Write(IOBuffer,*) '\n\nError opening EXO file ',trim(filename),'\n\n'
!            Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr);
!            STOP
!            !SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,IOBuffer);
!         EndIf
!         Call EXGINI(exoid,dummyS,numDim,junk1,junk2,numSet,junk3,junk3,exoerr)
!      End If
!      Call MPI_Bcast(numSet,1,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)
!      Call MPI_Bcast(numDim,1,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)
!
!      Allocate(elemType(numSet))
!      Allocate(setID(numSet))
!      If (MEF90Ctx%rank == 0) Then       
!         Call EXGEBI(exoID,setID,exoerr)
!      End If
!      Call MPI_Bcast(setID,numSet,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)
!      Do set = 1, numSet
!         If (MEF90Ctx%rank == 0) Then
!            Call EXGELB(exoID,setID(set),EXOelemType,junk1,junk2,junk3,exoerr)
!         EndIf
!         Call MPI_Bcast(EXOelemType,MXSTLN,MPI_CHAR,0,MEF90Ctx%comm,ierr)
!         Call EXO2MEF90ElementType_Elast(EXOelemType,numDim,elemType(set),ierr)
!      End Do
!      DeAllocate(setID)
!      If (MEF90Ctx%rank == 0) Then
!         Call EXCLOS(exoid,exoErr)
!      End If
!   End Subroutine EXOGetCellSetElementType_Elast
>>>>>>> f87946a (Updated m_MEF90_Materials, started work on elements)
=======
>>>>>>> 7aedf22 (Updated geometry initialization)
      
#undef __FUNCT__
#define __FUNCT__ "EXOWriteCase"
   Subroutine EXOWriteCase(prefix,formatstring,numprocs)
      Character(len=*)                               :: prefix,formatstring
      PetscInt                                       :: numprocs
      Integer                                        :: F_out = 999
      
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
   End Subroutine EXOWriteCase
   
#undef __FUNCT__
#define __FUNCT__ "MEF90EXOFormat"
!!!
!!!  
!!!  MEF90EXOFormat:
!!!  
!!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine MEF90EXOFormat(exoid,NameG,NameC,NameV,ierr)
   Integer,Intent(INOUT)                              :: exoid
   Character(len=*),Dimension(:),Pointer,Intent(IN)   :: nameG,nameC,nameV
   PetscErrorCode,Intent(OUT)                         :: ierr
   
   Integer                                            :: numCS
   Real                                               :: rJunk
   Character(len=MXSTLN)                              :: sJunk
   Integer,Dimension(:,:),Pointer                     :: varTruthTable
   
   If (exoID > 0) Then
      !!! Write variable names
      If (size(nameG) > 0) Then
         Call EXPVP (exoid,'g',size(nameG),ierr)
         Call EXPVAN(exoid,'g',size(nameG),NameG,ierr)
      End If
      If (size(nameC) > 0) Then
         Call EXPVP (exoid,'e',size(nameC),ierr)
         Call EXPVAN(exoid,'e',size(nameC),NameC,ierr)
      End If
      If (size(nameV) > 0) Then
         Call EXPVP (exoid,'n',size(nameV),ierr)
         Call EXPVAN(exoid,'n',size(nameV),nameV,ierr)
      End If

      !!! Write truth tables
      Call EXINQ(exoid,EXELBL,numCS,rjunk,sjunk,ierr)
      Allocate(varTruthTable(size(nameC),numCS))
      varTruthTable = 1
      Call EXPVTT(exoid,numCS,size(nameC),varTruthTable,ierr)
      DeAllocate(varTruthTable)
   End If
End Subroutine MEF90EXOFormat
End Module m_MEF90_EXO

