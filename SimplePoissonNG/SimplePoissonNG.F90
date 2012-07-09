Program  SimplePoissonNG
#include "SimplePoisson.inc"
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use M_POISSON
   Implicit NONE   


   Character(len=MEF90_MXSTRLEN)                :: prefix,filename
   Type(PoissonCtx_Type)                        :: AppCtx
   PetscErrorCode                               :: iErr
   PetscInt                                     :: exo_step,exo_field
   Character(len=MEF90_MXSTRLEN)                :: IOBuffer
   PetscInt                                     :: numDim
   PetscInt                                     :: Verbose
   PetscBool                                    :: Restart
   PetscBool                                    :: splitIO
   PetscBool                                    :: flg
   Type(SNES)                                   :: snesTemp
   Type(KSP)                                    :: kspTemp
   Type(PC)                                     :: pcTemp
   Type(Mat)                                    :: matTemp
   Type(Vec)                                    :: solTemp,resTemp,locTemp
   Integer                                      :: cpu_ws,io_ws,exoIN=0,exoOUT=0
   Real                                         :: exo_version
   Integer                                      :: exoerr
   MPI_Comm                                     :: IOComm

   Call MEF90_Initialize()
   Call m_Poisson_Initialize(ierr);CHKERRQ(ierr)

   Call PetscOptionsGetString(PETSC_NULL_CHARACTER,'--prefix',prefix,flg,ierr);CHKERRQ(ierr)
   If (.NOT. flg) Then
      Call PetscPrintf(PETSC_COMM_WORLD,"No mesh prefix given\n",ierr)
      Call MEF90_Finalize()
      STOP
   End If

   verbose = 0
   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'--verbose',verbose,flg,ierr);CHKERRQ(ierr)
   restart = PETSC_FALSE
   Call PetscOptionsGetBool(PETSC_NULL_CHARACTER,'--restart',restart,flg,ierr);CHKERRQ(ierr)
   
   splitIO = PETSC_FALSE
   Call PetscOptionsGetBool(PETSC_NULL_CHARACTER,'--splitIO',splitIO,flg,ierr);CHKERRQ(ierr)
   

   If (MEF90_MyRank == 0) Then
      cpu_ws = 8
      io_ws = 8
      filename = Trim(prefix)//'.gen'
      exoIN = EXOPEN(filename,EXREAD,cpu_ws,io_ws,exo_version,exoerr)
      If (verbose >0) Then
         Write(IOBuffer,99) exoERR, exoIN
         Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr)
      End If
99 Format('EXOPEN status: ',I4,' exoIN: ',I4,'\n')         
      If (exoerr < 0) Then
         Write(IOBuffer,*) '\n\nError opening EXO file ', trim(filename),'\n\n'
         Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr);
         SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,IOBuffer,ierr);
      EndIf
   End If
   If (splitIO) Then
      IOComm = PETSC_COMM_SELF
      Write(filename,105) trim(prefix),MEF90_MyRank
   Else
      IOComm = PETSC_COMM_WORLD
      Write(filename,106) trim(prefix)
   End If
105 Format(A,'-',I4.4,'.gen')
106 Format(A,'_out.gen')

   !!! 
   !!! Create PoissonCtx
   !!!
   Call PoissonCtxCreate(AppCtx,exoIN,ierr);CHKERRQ(ierr)
   
   !!!
   !!! format the output file
   !!!
   If (splitIO) Then
      cpu_ws = 8
      io_ws = 8
      exoOUT = EXCRE(trim(filename),EXCLOB,cpu_ws,io_ws,ierr)
      Call DMmeshViewExodusSplit(AppCtx%mesh,exoOUT,ierr)
      Call EXPVP (exoOUT,'g',3,ierr)
      Call EXPVAN(exoOUT,'g',3,(/'Elastic Energy ','Ext Forces work','Total Energy   '/),ierr)
      Call EXPVP (exoOUT,'n',1,ierr)
      Call EXPVAN(exoOUT,'n',1,(/'U'/),ierr)
      Call EXPVAN(exoOUT,'e',1,(/'F'/),ierr)
   Else
      If (MEF90_MyRank == 0) Then
         cpu_ws = 8
         io_ws = 8
         exoOUT = EXCRE(filename,EXCLOB,cpu_ws,io_ws,ierr)
         Call EXCOPY(exoIN,exoOUT,ierr)
         Call EXPVP (exoOUT,'g',3,ierr)
         Call EXPVAN(exoOUT,'g',3,(/'Elastic Energy ','Ext Forces work','Total Energy   '/),ierr)
         Call EXPVP (exoOUT,'n',1,ierr)
         Call EXPVAN(exoOUT,'n',1,(/'U'/),ierr)
         Call EXPVAN(exoOUT,'e',1,(/'F'/),ierr)
      End If
   End If

   !!!
   !!! Get Matrix for the Jacobian / SNES and unknown vector
   !!!
   Call DMMeshSetMaxDof(AppCtx%mesh,1,iErr); CHKERRQ(iErr) 
   Call DMMeshCreateVector(AppCtx%mesh,AppCtx%Section,solTemp,ierr);CHKERRQ(ierr)
   Call SectionRealCreateLocalVector(AppCtx%Section,locTemp,ierr);CHKERRQ(ierr)
   Call VecDuplicate(solTemp,resTemp,ierr);CHKERRQ(ierr)

   Call DMMeshCreateMatrix(AppCtx%mesh,AppCtx%Section,MATMPIAIJ,matTemp,iErr);CHKERRQ(iErr)
   !!! Not sure if this is still needed when using MatZeroRowsColumnsIS for BC handling
   Call MatSetOption(matTemp,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr);CHKERRQ(ierr)
   Call MatSetFromOptions(matTemp,ierr);CHKERRQ(ierr)
   !!!
   !!! Create SNES, 
   !!!
   Call SNESCreate(PETSC_COMM_WORLD,snesTemp,ierr);CHKERRQ(ierr)
   Call SNESSetDM(snesTemp,AppCtx%mesh,ierr);CHKERRQ(ierr)
   Call SNESSetOptionsPrefix(snesTemp,'temp_',ierr);CHKERRQ(ierr)
   Call SNESSetFromOptions(snesTemp,ierr);CHKERRQ(ierr)
   Call SNESSetFunction(snesTemp,resTemp,SimplePoissonNGOperatorAssembly,AppCtx,ierr);CHKERRQ(ierr)
   Call SNESSetJacobian(snesTemp,matTemp,matTemp,SimplePoissonNGBilinearFormAssembly,AppCtx,ierr);CHKERRQ(ierr)

   Call SNESGetKSP(snesTemp,kspTemp,ierr);CHKERRQ(ierr)
   Call KSPSetType(kspTemp,KSPCG,ierr);CHKERRQ(ierr)
   Call KSPSetInitialGuessNonzero(kspTemp,PETSC_TRUE,ierr);CHKERRQ(ierr)
   Call KSPSetFromOptions(kspTemp,ierr);CHKERRQ(ierr)

   Call KSPGetPC(kspTemp,pcTemp,ierr);CHKERRQ(ierr)
   Call PCSetFromOptions(pcTemp,ierr);CHKERRQ(ierr)
   !!! Setup GAMG here (coordinates, in particular)
   If (verbose > 0) Then
      Call SNESView(snesTemp,PETSC_VIEWER_STDOUT_WORLD,ierr)
   End If
   
   
   Call SimplePoissonFormInitialGuess(solTemp,AppCtx,ierr);CHKERRQ(ierr)
   Call SNESSolve(snesTemp,PETSC_NULL_OBJECT,solTemp,ierr);CHKERRQ(ierr)
   
   Call SectionRealToVec(AppCtx%Section,AppCtx%ScatterSecToVec,SCATTER_REVERSE,solTemp,ierr);CHKERRQ(ierr)
   !!! solTemp values have copied to the Section, localTemp values should magically be up to date in the local
   !!! Vector since it shares the same storage space.
   Call VecViewExodusVertex(AppCtx%mesh,locTemp,IOComm,exoOUT,1,1,ierr)
   

   !!!
   !!! Cleanup
   !!!
   Call EXCLOS(exoIN,ierr)
   Call EXCLOS(exoOUT,ierr)
   Call PoissonCtxDestroy(AppCtx,ierr);CHKERRQ(ierr)
   Call SNESDestroy(snesTemp,ierr);CHKERRQ(ierr)
   Call VecDestroy(solTemp,ierr);CHKERRQ(ierr)   
   Call VecDestroy(resTemp,ierr);CHKERRQ(ierr)   
   Call VecDestroy(locTemp,ierr);CHKERRQ(ierr)   
   Call MEF90_Finalize()


End Program  SimplePoissonNG
