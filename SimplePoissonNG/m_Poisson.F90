#include "SimplePoisson.inc"
Module m_Poisson
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_PoissonGlobalProperties
   Use M_POISSONCELLSETPROPERTIES
   Use m_PoissonVertexSetProperties
   Implicit NONE

   PetscSizeT,protected       :: sizeofPoissonCtx
   
Contains
#undef __FUNCT__
#define __FUNCT__ "m_Poisson_Initialize"
   Subroutine m_Poisson_Initialize(ierr)
      PetscErrorCode,intent(OUT)          :: ierr

      Type(MEF90Ctx_Type),Target          :: PoissonCtx
      character(len=1),pointer            :: dummychar(:)
      PetscSizeT                          :: sizeofchar
   
      Call PoissonGlobalPropertiesInitialize(ierr)
      Call PoissonCellSetPropertiesInitialize(ierr)
      Call PoissonVertexSetPropertiesInitialize(ierr)
      
      Call PetscDataTypeGetSize(PETSC_CHAR,sizeofchar,ierr)
      sizeofPoissonCtx = size(transfer(PoissonCtx,dummychar))*sizeofchar
   End Subroutine m_Poisson_Initialize
   
   
#undef __FUNCT__
#define __FUNCT__ "PoissonCtxDestroy"
!!!
!!!  
!!!  PoissonCtxDestroy:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine PoissonCtxDestroy(PoissonCtx,snesTemp,ierr)
   Type(MEF90Ctx_Type)                                      :: PoissonCtx
   Type(SNES),Intent(IN)                                    :: snesTemp
   PetscErrorCode,Intent(OUT)                               :: ierr

   Type(IS)                                                 :: setIS   
   Type(DM)                                                 :: mesh
   PetscInt                                                 :: e,set,nset
   
   Call SNESGetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)

   Call DMmeshGetLabelIdIS(mesh,'Cell Sets',SetIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr) 
   Call ISGetLocalSize(setIS,nset,ierr);CHKERRQ(ierr)
   
   Do set = 1, nset
      Call PetscBagDestroy(PoissonCtx%CellSetPropertiesBag(set),ierr);CHKERRQ(ierr)
      Call PetscBagDestroy(PoissonCtx%MaterialPropertiesBag(set),ierr);CHKERRQ(ierr)
   End Do
   DeAllocate(PoissonCtx%CellSetPropertiesBag,stat=ierr)
   DeAllocate(PoissonCtx%MaterialPropertiesBag,stat=ierr)
   Call ISDestroy(setIS,ierr);CHKERRQ(ierr)

   Call DMmeshGetLabelIdIS(mesh,'Vertex Sets',SetIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr) 
   Call ISGetLocalSize(setIS,nset,ierr);CHKERRQ(ierr)
   Do set = 1, nset
      Call PetscBagDestroy(PoissonCtx%VertexSetPropertiesBag(set),ierr);CHKERRQ(ierr)
   End Do
   DeAllocate(PoissonCtx%VertexSetPropertiesBag)
   Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   
   Call PetscBagDestroy(PoissonCtx%GlobalPropertiesBag,ierr);CHKERRQ(ierr)
End Subroutine PoissonCtxDestroy 


#undef __FUNCT__
#define __FUNCT__ "SimplePoissonEXOOpenInputFile"
!!!
!!!  
!!!  SimplePoissonEXOOpenInputFile:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine SimplePoissonEXOOpenInputFile(prefix,exoIN,ierr)
   Character(len=*),intent(IN)                     :: prefix
   Integer,intent(OUT)                             :: exoIN
   PetscErrorCode,Intent(OUT)                      :: ierr

   Character(len=MEF90_MXSTRLEN)                   :: IOBuffer,filename
   Integer                                         :: cpu_ws,io_ws
   Real                                            :: exo_version
   Integer                                         :: exoerr
   
   !!! Open input file
   If (MEF90_MyRank == 0) Then
      cpu_ws = 8
      io_ws = 8
      filename = Trim(prefix)//'.gen'
      exoIN = EXOPEN(filename,EXREAD,cpu_ws,io_ws,exo_version,exoerr)
      If (exoerr < 0) Then
         Write(IOBuffer,*) '\n\nError opening EXO file ',trim(filename),'\n\n'
         Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr);
         SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,IOBuffer,ierr);
      EndIf
   End If
End Subroutine SimplePoissonEXOOpenInputFile

#undef __FUNCT__
#define __FUNCT__ "SimplePoissonPrepareOutputEXO"
!!!
!!!  
!!!  SimplePoissonPrepareOutputEXO:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine SimplePoissonPrepareOutputEXO(mesh,prefix,exoIN,exoOUT,Poissonctx,ierr)
   Type(DM),intent(IN)                             :: mesh
   Character(len=*),intent(IN)                     :: prefix
   Integer,intent(IN)                              :: exoIN
   Integer,intent(OUT)                             :: exoOUT
   Type(MEF90Ctx_Type),intent(IN)                  :: PoissonCtx
   PetscErrorCode,Intent(OUT)                      :: ierr

   MPI_Comm                                        :: IOComm
   Integer                                         :: IORank
   Integer                                         :: cpu_ws,io_ws
   Real                                            :: exo_version
   Integer                                         :: exoerr
   Character(len=MEF90_MXSTRLEN)                   :: IOBuffer,filename
   Type(PoissonGlobalProperties_Type),pointer      :: GlobalProperties 
   
   Call PetscBagGetDataPoissonGlobalProperties(PoissonCtx%GlobalPropertiesBag,GlobalProperties,ierr);CHKERRQ(ierr)      
   
   !!! Get name of output file
   If (GlobalProperties%FileFormat == Poisson_EXOSplit) Then
      IOComm = PETSC_COMM_SELF
      Write(filename,100) trim(prefix),MEF90_MyRank
   Else
      IOComm = PETSC_COMM_WORLD
      Write(filename,101) trim(prefix)
   End If
100 Format(A,'-',I4.4,'.gen')
101 Format(A,'_out.gen')
   Call MPI_Comm_Rank(IOComm,IORank,ierr)

   !!! Open output file or create it and format it depending on loading type
   If (IORank == 0) Then
      If (GlobalProperties%LoadingType == Poisson_FILE) Then
         !!! We assume that the output file already exists
         cpu_ws = 8
         io_ws = 8
         exoOUT = EXOPEN(filename,EXWRIT,cpu_ws,io_ws,exo_version,exoerr)
         If (exoerr < 0) Then
            Write(IOBuffer,*) '\n\nError opening EXO file ',trim(filename),'\n\n'
            Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr);
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,IOBuffer,ierr);
         EndIf
      Else
         !!! Create output file, add topology informations, add EXO format
         cpu_ws = 8
         io_ws = 8
         exoOUT = EXCRE(trim(filename),EXCLOB,cpu_ws,io_ws,ierr)
         If (GlobalProperties%FileFormat == Poisson_EXOSplit) Then
            Call DMmeshViewExodusSplit(mesh,exoOUT,ierr)
         Else
            Call EXCOPY(exoIN,exoOUT,ierr)
         End If            
         Call EXPVP (exoOUT,'g',3,ierr)
         Call EXPVAN(exoOUT,'g',3,(/'Energy      ','Fluxes work ','Total Energy'/),ierr)
         Call EXPVP (exoOUT,'n',1,ierr)
         Call EXPVAN(exoOUT,'n',1,(/'U'/),ierr)
         !Call EXPVAN(exoOUT,'e',1,(/'F'/),ierr)
      End If
   End If
End Subroutine SimplePoissonPrepareOutputEXO

#undef __FUNCT__
#define __FUNCT__ "SimplePoissonSaveEXO"
!!!
!!!  
!!!  SimplePoissonSaveEXO:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine SimplePoissonSaveEXO(mesh,EXOout,sol,TimeStepNum,t,energy,work,PoissonCtx,ierr)
   Type(DM),Intent(IN)                             :: mesh
   Integer,Intent(IN)                              :: EXOout
   Type(Vec),Intent(IN)                            :: sol
   PetscInt,Intent(IN)                             :: TimeStepNum
   PetscReal,Intent(IN)                            :: t
   PetscReal,Intent(IN)                            :: energy,work
   Type(MEF90Ctx_Type),intent(IN)                  :: PoissonCtx
   PetscErrorCode,Intent(OUT)                      :: ierr
   
   MPI_Comm                                        :: IOComm
   Integer                                         :: IORank
   Type(Vec)                                       :: locTemp
   Type(SectionReal)                               :: secTemp
   Type(VecScatter)                                :: ScatterSecToVec
   Type(PoissonGlobalProperties_Type),pointer      :: GlobalProperties 
   Character(len=MEF90_MXSTRLEN)                   :: IOBuffer
   Integer                                         :: tempOffset
   
   Call PetscBagGetDataPoissonGlobalProperties(PoissonCtx%GlobalPropertiesBag,GlobalProperties,ierr);CHKERRQ(ierr)      
   If (GlobalProperties%FileFormat == Poisson_EXOSplit) Then
      IOComm = PETSC_COMM_SELF
   Else
      IOComm = PETSC_COMM_WORLD
   End If
   Call MPI_Comm_Rank(IOComm,IORank,ierr)
   
   !!! Just as for global vectors, there are two ways to deal with local vectors
   !!! Obtaining them from the DM, and using DMLocalToGlobalBegin/End
   !Call DMGetLocalVector(mesh,locTemp,ierr);CHKERRQ(ierr)
   !Call DMGlobalToLocalBegin(mesh,solTemp,INSERT_VALUES,locTemp,ierr);CHKERRQ(ierr)
   !Call DMGlobalToLocalEnd(mesh,solTemp,INSERT_VALUES,locTemp,ierr);CHKERRQ(ierr)
   !!! Or by from a SectionReal, and using SectionRealToVec to copy from Vec to Section
   !!! This would require passing the section as argument, or pulling it from the DMMesh
   Call DMMeshGetSectionReal(mesh,'default',secTemp,ierr);CHKERRQ(ierr)
   Call DMMeshCreateGlobalScatter(mesh,secTemp,ScatterSecToVec,ierr);CHKERRQ(ierr)
   Call SectionRealCreateLocalVector(secTemp,locTemp,ierr);CHKERRQ(ierr)   
   Call SectionRealToVec(secTemp,ScatterSecToVec,SCATTER_REVERSE,sol,ierr);CHKERRQ(ierr)
   Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)


   If (GlobalProperties%LoadingType == Poisson_FILE) Then
      tempOffset = GlobalProperties%tempOffset
   Else
      tempOffset = 1
   EndIf
   
   Call VecViewExodusVertex(mesh,locTemp,IOComm,exoOUT,TimeStepNum,tempoffset,ierr);CHKERRQ(ierr)
   If ((GlobalProperties%LoadingType /= Poisson_FILE) .AND. (IORank == 0)) Then
      Call EXPGV(exoOUT,timeStepNum,3,(/ energy,work,energy-work /),ierr)   
      If (GlobalProperties%FileMode == Poisson_Replace) Then
         Call EXPTIM(exoOUT,timeStepNum,t,ierr)
      End If
   End If
   !!!
   !!! Restore local vector if obtained from the DM, 
   !Call DMrestoreLocalVector(mesh,locTemp,ierr);CHKERRQ(ierr)   
   Call SectionRealDestroy(secTemp,ierr);CHKERRQ(ierr)
   Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
   If (GlobalProperties%verbose > 0) Then
      Write(IOBuffer,100) tempOffset,TimeStepNum
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
   End If
100 Format('Wrote temperature at offset ',I2,' for time step ',I4,'\n');CHKERRQ(ierr)
End Subroutine SimplePoissonSaveEXO

#undef __FUNCT__
#define __FUNCT__ "SimplePoissonLoadEXO"
!!!
!!!  
!!!  SimplePoissonLoadEXO:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine SimplePoissonLoadEXO(mesh,EXOIN,bctemp,flux,reftemp,TimeStepNum,PoissonCtx,ierr)
   Type(DM),Intent(IN)                             :: mesh
   Integer,Intent(IN)                              :: EXOIN
   Type(Vec),Intent(IN),optional                   :: bctemp,flux,refTemp
   PetscInt,Intent(IN)                             :: TimeStepNum
   Type(MEF90Ctx_Type),intent(IN)                  :: PoissonCtx
   PetscErrorCode,Intent(OUT)                      :: ierr

   MPI_Comm                                        :: IOComm
   Integer                                         :: IORank
   Type(Vec)                                       :: locVec
   Type(SectionReal)                               :: IOsec
   Type(VecScatter)                                :: ScatterSecToVec
   Type(PoissonGlobalProperties_Type),pointer      :: GlobalProperties 
   Character(len=MEF90_MXSTRLEN)                   :: IOBuffer
   
   Call PetscBagGetDataPoissonGlobalProperties(PoissonCtx%GlobalPropertiesBag,GlobalProperties,ierr);CHKERRQ(ierr)      
   If (GlobalProperties%FileFormat == Poisson_EXOSplit) Then
      IOComm = PETSC_COMM_SELF
   Else
      IOComm = PETSC_COMM_WORLD
   End If
   Call MPI_Comm_Rank(IOComm,IORank,ierr)

   Call DMMeshGetSectionReal(mesh,'default',IOsec,ierr);CHKERRQ(ierr)
   Call DMMeshCreateGlobalScatter(mesh,IOsec,ScatterSecToVec,ierr);CHKERRQ(ierr)
   Call SectionRealCreateLocalVector(IOsec,locVec,ierr);CHKERRQ(ierr)   

   Call VecLoadExodusVertex(mesh,locVec,IOComm,exoIN,TimeStepNum,GlobalProperties%bctempOffset,ierr);CHKERRQ(ierr)
   Call SectionRealToVec(IOsec,ScatterSecToVec,SCATTER_FORWARD,bctemp,ierr);CHKERRQ(ierr)

   Call VecLoadExodusVertex(mesh,locVec,IOComm,exoIN,TimeStepNum,GlobalProperties%fluxOffset,ierr);CHKERRQ(ierr)
   Call SectionRealToVec(IOsec,ScatterSecToVec,SCATTER_FORWARD,flux,ierr);CHKERRQ(ierr)
   
   Call VecLoadExodusVertex(mesh,locVec,IOComm,exoIN,TimeStepNum,GlobalProperties%refTempOffset,ierr);CHKERRQ(ierr)
   Call SectionRealToVec(IOsec,ScatterSecToVec,SCATTER_FORWARD,refTemp,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(IOsec,ierr);CHKERRQ(ierr)
   Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
End Subroutine SimplePoissonLoadEXO

#undef __FUNCT__
#define __FUNCT__ "SimplePoissonFormInitialGuess_Cst"
!!!
!!!  
!!!  SimplePoissonFormInitialGuess_Cst:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine SimplePoissonFormInitialGuess_Cst(snesTemp,x,t,PoissonCtx,ierr)
   Type(SNES),Intent(IN)                           :: snesTemp
   Type(Vec),Intent(IN)                            :: x
   PetscReal,Intent(IN)                            :: t
   Type(MEF90Ctx_Type),intent(IN)                  :: PoissonCtx
   PetscErrorCode,Intent(OUT)                      :: ierr

   Type(PoissonVertexSetProperties_Type),pointer   :: vertexSetProperties
   Type(IS)                                        :: VertexSetGlobalIS,setIS,setISdof
   PetscInt,dimension(:),Pointer                   :: setID
   PetscInt,Dimension(:),Pointer                   :: setIdx
   PetscInt                                        :: set
   PetscReal,Dimension(:),Pointer                  :: BC
   Type(DM)                                        :: mesh
   Type(SectionReal)                               :: xSec
   
   Call SNESGetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)
   Call DMMeshGetSectionReal(mesh,'default',xSec,ierr);CHKERRQ(ierr)

   Call DMmeshGetLabelIdIS(mesh,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
   Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      Call PetscBagGetDataPoissonVertexSetProperties(PoissonCtx%VertexSetPropertiesBag(set),vertexSetProperties,ierr);CHKERRQ(ierr)
      If (vertexSetProperties%Has_BC) Then
         Call DMMeshGetStratumIS(mesh,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call DMMeshISCreateISglobaldof(mesh,xSec,setIS,0,setISdof,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
         Allocate(BC(size(setIdx)),stat=ierr)
         BC = vertexSetProperties%BC * t
         Call VecSetValues(x,size(setIdx),setIdx,BC,INSERT_VALUES,ierr);CHKERRQ(ierr)
         DeAllocate(BC)
         Call ISRestoreIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
      End If
   End Do
   Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)

   Call VecAssemblyBegin(x,ierr);CHKERRQ(ierr)
   Call VecAssemblyEnd(x,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
End Subroutine SimplePoissonFormInitialGuess_Cst

#undef __FUNCT__
#define __FUNCT__ "SimplePoissonFormInitialGuess"
!!!
!!!  
!!!  SimplePoissonFormInitialGuess:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine SimplePoissonFormInitialGuess(snesTemp,x,xbc,PoissonCtx,ierr)
   Type(SNES),Intent(IN)                           :: snesTemp
   Type(Vec),Intent(IN)                            :: x
   Type(Vec),Intent(IN)                            :: xbc
   Type(MEF90Ctx_Type),intent(IN)                  :: PoissonCtx
   PetscErrorCode,Intent(OUT)                      :: ierr

   Type(PoissonVertexSetProperties_Type),pointer   :: vertexSetProperties
   Type(IS)                                        :: VertexSetGlobalIS,setIS,setISdof
   PetscInt,dimension(:),Pointer                   :: setID
   PetscInt,Dimension(:),Pointer                   :: setIdx
   PetscInt                                        :: set
   PetscReal,Dimension(:),Pointer                  :: BC
   Type(DM)                                        :: mesh
   Type(SectionReal)                               :: xSec
   
   Call SNESGetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)
   Call DMMeshGetSectionReal(mesh,'default',xSec,ierr);CHKERRQ(ierr)

   Call DMmeshGetLabelIdIS(mesh,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
   Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      Call PetscBagGetDataPoissonVertexSetProperties(PoissonCtx%VertexSetPropertiesBag(set),vertexSetProperties,ierr);CHKERRQ(ierr)
      If (vertexSetProperties%Has_BC) Then
         Call DMMeshGetStratumIS(mesh,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call DMMeshISCreateISglobaldof(mesh,xSec,setIS,0,setISdof,ierr);CHKERRQ(ierr)
         Call ISGetIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
         Allocate(BC(size(setIdx)),stat=ierr)
         Call VecGetValues(xbc,size(setIdx),setIdx,BC,ierr);CHKERRQ(ierr)
         Call VecSetValues(x,size(setIdx),setIdx,BC,INSERT_VALUES,ierr);CHKERRQ(ierr)
         DeAllocate(BC)
         Call ISRestoreIndicesF90(setISdof,setIdx,ierr);CHKERRQ(ierr)
      End If
   End Do
   Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)

   Call VecAssemblyBegin(x,ierr);CHKERRQ(ierr)
   Call VecAssemblyEnd(x,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
End Subroutine SimplePoissonFormInitialGuess

#undef __FUNCT__
#define __FUNCT__ "SimplePoissonGetTime"
!!!
!!!  
!!!  SimplePoissonGetTime:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine SimplePoissonGetTime(time,exoIN,PoissonCtx,ierr)
   PetscReal,Dimension(:),Pointer                  :: time
   Integer,Intent(IN)                              :: exoIN
   Type(MEF90Ctx_Type),intent(IN)                  :: PoissonCtx
   PetscErrorCode,Intent(OUT)                      :: ierr

   Type(PoissonGlobalProperties_Type),pointer      :: GlobalProperties
   MPI_Comm                                        :: IOComm
   Integer                                         :: IORank
   PetscInt                                        :: NumTimeStep,i
   Real                                            :: rdumm
   Character(len=1)                                :: cdumm
   Integer                                         :: exoerr

   Call PetscBagGetDataPoissonGlobalProperties(PoissonCtx%GlobalPropertiesBag,GlobalProperties,ierr);CHKERRQ(ierr)      
   Select Case (GlobalProperties%LoadingType)
      Case (Poisson_CST)
         Allocate(Time(1),stat=ierr)
         time = 1.0_Kr
      Case (Poisson_MIL)
         Allocate(time(GlobalProperties%numTimeStep))
         Do i = 1,GlobalProperties%numTimeStep-1
            time(i) = GlobalProperties%timeMin + (GlobalProperties%timeMax - GlobalProperties%timeMin) / (GlobalProperties%numTimeStep - 1.0_Kr) * (i-1.0_Kr)
         End Do
         time(GlobalProperties%numTimeStep) = GlobalProperties%timeMax
      Case (Poisson_FILE)
         If (GlobalProperties%FileFormat == Poisson_EXOSplit) Then
            Call EXINQ(exoIN,EXTIMS,numTimeStep,rdumm,cdumm,exoerr)
            Allocate(time(numTimeStep))
            Call EXGATM(exoIN,time,exoerr)
         Else
            If (MEF90_MyRank == 0) Then
               Call EXINQ(exoIN,EXTIMS,numTimeStep,rdumm,cdumm,exoerr)
               Call MPI_Bcast(numTimeStep,1,MPIU_INTEGER,0,PETSC_COMM_WORLD,ierr)
               Allocate(time(numTimeStep))
               Call EXGATM(exoIN,time,exoerr)
               Call MPI_Bcast(time,numTimeStep,MPIU_INTEGER,0,PETSC_COMM_WORLD,ierr)
            Else
               Call MPI_Bcast(numTimeStep,1,MPIU_INTEGER,0,PETSC_COMM_WORLD,ierr)
               Allocate(time(numTimeStep))
               Call MPI_Bcast(time,numTimeStep,MPIU_INTEGER,0,PETSC_COMM_WORLD,ierr)
            End If         
         End If
   End Select
End Subroutine SimplePoissonGetTime

End Module m_Poisson
