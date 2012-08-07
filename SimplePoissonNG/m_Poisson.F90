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
#define __FUNCT__ "PoissonEXOOpenInputFile"
!!!
!!!  
!!!  PoissonEXOOpenInputFile:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine PoissonEXOOpenInputFile(prefix,exoIN,ierr)
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
End Subroutine PoissonEXOOpenInputFile

#undef __FUNCT__
#define __FUNCT__ "PoissonPrepareOutput"
!!!
!!!  
!!!  PoissonPrepareOutput:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine PoissonPrepareOutput(mesh,prefix,exoIN,exoOUT,IOComm,Poissonctx,ierr)
   Type(DM),intent(IN)                             :: mesh
   Character(len=*),intent(IN)                     :: prefix
   Integer,intent(IN)                              :: exoIN
   Integer,intent(OUT)                             :: exoOUT
   Type(MEF90Ctx_Type),intent(IN)                  :: PoissonCtx
   MPI_Comm, intent(OUT)                           :: IOComm
   PetscErrorCode,Intent(OUT)                      :: ierr

   Integer                                         :: cpu_ws,io_ws
   Real                                            :: exo_version
   Integer                                         :: exoerr
   Character(len=MEF90_MXSTRLEN)                   :: IOBuffer,filename
   Type(PoissonGlobalProperties_Type),pointer      :: GlobalProperties 
   PetscBool                                       :: IONode != PETSC_FALSE
   
   IONode = PETSC_FALSE
   Call PetscBagGetDataPoissonGlobalProperties(PoissonCtx%GlobalPropertiesBag,GlobalProperties,ierr);CHKERRQ(ierr)      
   
   !!! Get name of output file
   If (GlobalProperties%FileFormat == Poisson_EXOSplit) Then
      IOComm = PETSC_COMM_SELF
      IONode = PETSC_TRUE
      Write(filename,100) trim(prefix),MEF90_MyRank
   Else
      IOComm = PETSC_COMM_WORLD
      If (MEF90_MyRank == 0) Then
         IONode = PETSC_TRUE
      End If
      Write(filename,101) trim(prefix)
   End If
100 Format(A,'-',I4.4,'.gen')
101 Format(A,'_out.gen')

   !!! Open output file or create it and format it depending on loading type
   If (IONode) Then
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
End Subroutine PoissonPrepareOutput

#undef __FUNCT__
#define __FUNCT__ "SimplePoissonFormInitialGuess"
!!!
!!!  
!!!  SimplePoissonFormInitialGuess:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine SimplePoissonFormInitialGuess(snesTemp,x,PoissonCtx,ierr)
   Type(SNES),Intent(IN)                           :: snesTemp
   Type(Vec),Intent(IN)                            :: x
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
         Allocate(BC(size(setIdx)))
         BC = vertexSetProperties%BC
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


End Module m_Poisson
