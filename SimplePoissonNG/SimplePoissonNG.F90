!!!
!!! Try: mpiexec -n 2 ./Darwin-intel11.1-mef90-g/SimplePoissonNG2D -vs0001_tempBC 1 -vs0001_temp 0 -vs0002_tempBC 1 -vs0002_temp 0 -cs0001_force 0. -cs0002_force 1. -temp_pc_type bjacobi -temp_ksp_monitor -temp_snes_type ksponly --prefix SquareNG2
!!! 
Program  SimplePoissonNG
#include "SimplePoisson.inc"
#include <finclude/petscdef.h>
#include <finclude/petscbagdef.h>
   Use m_MEF90
   Use M_POISSON
   Use petsc
   Implicit NONE   

   Character(len=MEF90_MXSTRLEN)                :: prefix,filename
   Type(MEF90Ctx_Type)                          :: MEF90Ctx
   Type(DM)                                     :: mesh,tmp_mesh
   PetscErrorCode                               :: iErr
   PetscInt                                     :: exo_step,exo_field
   Character(len=MEF90_MXSTRLEN)                :: IOBuffer
   PetscInt                                     :: numDim
   PetscInt                                     :: Verbose
   PetscBool                                    :: splitIO
   PetscBool                                    :: flg
   Type(SNES)                                   :: snesTemp
   Type(KSP)                                    :: kspTemp
   Type(PC)                                     :: pcTemp
   Type(Mat)                                    :: matTemp
   Type(Vec)                                    :: solTemp,resTemp,locTemp
   !Type(Vec)                                    :: ubTemp,lbTemp
   Integer                                      :: cpu_ws,io_ws,exoIN=0,exoOUT=0
   Real                                         :: exo_version
   Integer                                      :: exoerr
   MPI_Comm                                     :: IOComm
   PetscReal,Dimension(:),Pointer               :: energy,work
   PetscInt,dimension(:),Pointer                :: setID
   PetscInt                                     :: set
   Type(MatNullSpace)                           :: nspTemp
   SNESConvergedReason                          :: reasonTemp
   PetscInt                                     :: itsTemp
   Type(SectionReal)                            :: secTemp
   Type(VecScatter)                             :: ScatterSecToVec
   
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
   
   splitIO = PETSC_FALSE
   Call PetscOptionsGetBool(PETSC_NULL_CHARACTER,'--splitIO',splitIO,flg,ierr);CHKERRQ(ierr)

   If (MEF90_MyRank == 0) Then
      cpu_ws = 8
      io_ws = 8
      filename = Trim(prefix)//'.gen'
      exoIN = EXOPEN(filename,EXREAD,cpu_ws,io_ws,exo_version,exoerr)
      If (verbose > 0) Then
         Write(IOBuffer,99) exoERR,exoIN
         Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr)
      End If
99 Format('EXOPEN status: ',I4,' exoIN: ',I4,'\n')         
      If (exoerr < 0) Then
         Write(IOBuffer,*) '\n\nError opening EXO file ',trim(filename),'\n\n'
         Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr);
         SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,IOBuffer,ierr);
      EndIf
   End If
   If (splitIO) Then
      IOComm = PETSC_COMM_SELF
      Write(filename,100) trim(prefix),MEF90_MyRank
   Else
      IOComm = PETSC_COMM_WORLD
      Write(filename,101) trim(prefix)
   End If
100 Format(A,'-',I4.4,'.gen')
101 Format(A,'_out.gen')

   !!!
   !!! Read DMMesh from exoIN
   !!!
   If (MEF90_NumProcs == 1) Then
      Call DMmeshCreateExodusNG(PETSC_COMM_WORLD,exoIN,mesh,ierr);CHKERRQ(ierr)
   Else
      Call DMmeshCreateExodusNG(PETSC_COMM_WORLD,exoIN,Tmp_mesh,ierr);CHKERRQ(ierr)
   
      Call DMmeshDistribute(Tmp_mesh,PETSC_NULL_CHARACTER,mesh,ierr);CHKERRQ(ierr)
      Call DMDestroy(Tmp_mesh,ierr);CHKERRQ(ierr)
   End If

   !!! 
   !!! Create SNES and associate mesh
   !!!
   Call SNESCreate(PETSC_COMM_WORLD,snesTemp,ierr);CHKERRQ(ierr)
   Call SNESSetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)
   Call SNESSetOptionsPrefix(snesTemp,'temp_',ierr);CHKERRQ(ierr)

   !!! 
   !!! Create PoissonCtx
   !!!
   Call PoissonCtxCreate(MEF90Ctx,snesTemp,exoIN,ierr);CHKERRQ(ierr)
   Call DMMeshGetSectionReal(mesh,'default',secTemp,ierr);CHKERRQ(ierr)
   !!!
   !!! Get Matrix for the Jacobian / SNES and unknown vector
   !!!
   Call DMMeshSetMaxDof(mesh,1,iErr); CHKERRQ(iErr) 
   Call DMMeshCreateMatrix(mesh,secTemp,MATAIJ,matTemp,iErr);CHKERRQ(iErr)

   !!! There are two ways to get global and local vectors:
   !!! Through the generic DM interface, which will only pull the layout from the 'default' section
   !Call DMGetGlobalVector(mesh,solTemp,ierr);CHKERRQ(ierr)
   !!! Through DMMeshCreateVector which lets use any layout
   Call DMMeshCreateVector(mesh,secTemp,solTemp,ierr);CHKERRQ(ierr)

   Call VecDuplicate(solTemp,resTemp,ierr);CHKERRQ(ierr)
      
   !!! Adding a null space when some boundary conditions are prescribes breaks everything...
   !!! Need to add a flag and make adding the null space optional
   !Call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,PETSC_NULL_OBJECT,nspTemp,ierr);CHKERRQ(ierr)
   !Call MatSetNullSpace(matTemp,nspTemp,ierr);CHKERRQ(ierr)

   !!! Not sure if this is still needed when using MatZeroRowsColumnsIS for BC handling
   Call MatSetOption(matTemp,MAT_SPD,PETSC_TRUE,ierr);CHKERRQ(ierr)
   Call MatSetOption(matTemp,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr);CHKERRQ(ierr)
   Call MatSetOption(matTemp,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr);CHKERRQ(ierr)
   Call MatSetFromOptions(matTemp,ierr);CHKERRQ(ierr)
   !!!
   !!! Create SNES,
   !!!
   !Call SNESCreate(PETSC_COMM_WORLD,snesTemp,ierr);CHKERRQ(ierr)
   !Call SNESSetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)
   !Call SNESSetOptionsPrefix(snesTemp,'temp_',ierr);CHKERRQ(ierr)
   Call SNESSetFunction(snesTemp,resTemp,SimplePoissonOperator,MEF90Ctx,ierr);CHKERRQ(ierr)
   Call SNESSetJacobian(snesTemp,matTemp,matTemp,SimplePoissonBilinearForm,MEF90Ctx,ierr);CHKERRQ(ierr)

   !!!
   !!! Testing SNESVI: does not look ready for prime time...
   !!!
   !Call VecDuplicate(solTemp,lbTemp,ierr);CHKERRQ(ierr)
   !Call VecSet(lbTemp,-1.0_Kr,ierr);CHKERRQ(ierr)
   !Call VecDuplicate(solTemp,ubTemp,ierr);CHKERRQ(ierr)
   !Call VecSet(ubTemp,0.2_Kr,ierr);CHKERRQ(ierr)
   !Call SNESVISetVariableBounds(snesTemp,lbTemp,ubTemp,ierr);CHKERRQ(ierr)

   Call SNESSetFromOptions(snesTemp,ierr);CHKERRQ(ierr)

   Call SNESGetKSP(snesTemp,kspTemp,ierr);CHKERRQ(ierr)
   Call KSPSetType(kspTemp,KSPCG,ierr);CHKERRQ(ierr)
   Call KSPSetInitialGuessNonzero(kspTemp,PETSC_TRUE,ierr);CHKERRQ(ierr)
!!!   Call KSPSetNullSpace(kspTemp,nspTemp,ierr);CHKERRQ(ierr)
   Call KSPSetFromOptions(kspTemp,ierr);CHKERRQ(ierr)

   Call KSPGetPC(kspTemp,pcTemp,ierr);CHKERRQ(ierr)
   Call PCSetFromOptions(pcTemp,ierr);CHKERRQ(ierr)
   !!! Setup GAMG here (coordinates,in particular)
   
   !!! Solve Poisson Equation
   Call SimplePoissonFormInitialGuess(snesTemp,solTemp,MEF90Ctx,ierr);CHKERRQ(ierr)
   Call SNESSolve(snesTemp,PETSC_NULL_OBJECT,solTemp,ierr);CHKERRQ(ierr)

   !!! Check SNES / KSP convergence
   Call SNESGetConvergedReason(snesTemp,reasonTemp,ierr);CHKERRQ(ierr)
   Call SNESGetIterationNumber(snesTemp,itsTemp,ierr);CHKERRQ(ierr)
   Write(IOBuffer,110) itsTemp,reasonTemp
110 Format('SNESTemp converged in in ',I4,' iterations. SNESConvergedReason is ', I4,'\n')
 
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
   
   !!! Compute energy and work
   Call ISGetIndicesF90(MEF90Ctx%CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Allocate(energy(size(setID)))
   Allocate(work(size(setID)))

   Call SimplePoissonEnergies(snesTemp,solTemp,MEF90Ctx,energy,work,ierr)
   Do set = 1,size(setID)
      Write(IOBuffer,102) setID(set),energy(set),work(set)
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
   End Do
   Call ISRestoreIndicesF90(MEF90Ctx%CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
102 Format('Cell set ',I4.4,' energy: ',ES12.5,' work: ',ES12.5,'\n')


   !!! Save results
   
   !!!
   !!! format the output file
   !!!
   If (splitIO) Then
      cpu_ws = 8
      io_ws = 8
      exoOUT = EXCRE(trim(filename),EXCLOB,cpu_ws,io_ws,ierr)
      Call DMmeshViewExodusSplit(mesh,exoOUT,ierr)
      Call EXPVP (exoOUT,'g',3,ierr)
      Call EXPVAN(exoOUT,'g',3,(/'Energy         ','Ext Forces work','Total Energy   '/),ierr)
      Call EXPVP (exoOUT,'n',1,ierr)
      Call EXPVAN(exoOUT,'n',1,(/'U'/),ierr)
      !Call EXPVAN(exoOUT,'e',1,(/'F'/),ierr)
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
   !!! Write in exo file
   !!!
   !!! Just as for global vectors, there are two ways to deal with local vectors
   !!! Obtaining them from the DM, and using DMLocalToGlobalBegin/End
   !Call DMGetLocalVector(mesh,locTemp,ierr);CHKERRQ(ierr)
   !Call DMGlobalToLocalBegin(mesh,solTemp,INSERT_VALUES,locTemp,ierr);CHKERRQ(ierr)
   !Call DMGlobalToLocalEnd(mesh,solTemp,INSERT_VALUES,locTemp,ierr);CHKERRQ(ierr)
   !!! Or by from a SectionReal, and using SectionRealToVec to copy from Vec to Section
   Call DMMeshCreateGlobalScatter(mesh,secTemp,ScatterSecToVec,ierr);CHKERRQ(ierr)
   Call SectionRealCreateLocalVector(secTemp,locTemp,ierr);CHKERRQ(ierr)   
   Call SectionRealToVec(secTemp,ScatterSecToVec,SCATTER_REVERSE,solTemp,ierr);CHKERRQ(ierr)
   Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)

   Call VecViewExodusVertex(mesh,locTemp,IOComm,exoOUT,1,1,ierr)
   Call EXPGV(exoOUT,1,3,(/ sum(energy),sum(work),sum(energy)-sum(work) /),ierr)   
   Write(IOBuffer,103) sum(energy),sum(work)
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
103 Format('Total:        energy: ',ES12.5,' work: ',ES12.5,'\n')
   Call DMrestoreLocalVector(mesh,locTemp,ierr);CHKERRQ(ierr)

   !!!
   !!! Cleanup
   !!!
   Call EXCLOS(exoIN,ierr)
   Call EXCLOS(exoOUT,ierr)
   Call PoissonCtxDestroy(MEF90Ctx,ierr);CHKERRQ(ierr)
   Call SNESDestroy(snesTemp,ierr);CHKERRQ(ierr)
   Call VecDestroy(solTemp,ierr);CHKERRQ(ierr)
   Call VecDestroy(resTemp,ierr);CHKERRQ(ierr)
   Call VecDestroy(locTemp,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(secTemp,ierr);CHKERRQ(ierr)
   Call DMDestroy(mesh,ierr);CHKERRQ(ierr);

   DeAllocate(work)
   DeAllocate(energy)
   
   If (verbose > 0) Then
      Call PetscOptionsView(PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
   EndIf
   Call MEF90_Finalize()


End Program  SimplePoissonNG
