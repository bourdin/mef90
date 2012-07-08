Program  SimplePoissonNG
#include "SimplePoisson.inc"
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use M_POISSON
   Implicit NONE   


   Character(len=MEF90_MXSTRLEN)                :: prefix
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
   Type(Vec)                                    :: solTemp,resTemp
   
      MatStructure                     :: matflg


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
   
   !!! 
   !!! Create PoissonCtx
   !!!
   Call PoissonCtxCreate(AppCtx,prefix,ierr);CHKERRQ(ierr)
   
   !!!
   !!! Get Matrix for the Jacobian / SNES and unknown vector
   !!!
   Call DMMeshSetMaxDof(AppCtx%mesh,1,iErr); CHKERRQ(iErr) 
   Call DMMeshCreateVector(AppCtx%mesh,AppCtx%Section,solTemp,ierr);CHKERRQ(ierr)
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
   
   
   !!!Call FormInitialSolution(solTemp,AppCtx,ierr);CHKERRQ(ierr)
   Call SNESSolve(snesTemp,PETSC_NULL_OBJECT,solTemp,ierr);CHKERRQ(ierr)
   !!!
   !!! Cleanup
   !!!
   Call PoissonCtxDestroy(AppCtx,ierr);CHKERRQ(ierr)
   Call SNESDestroy(snesTemp,ierr);CHKERRQ(ierr)
   Call VecDestroy(solTemp,ierr);CHKERRQ(ierr)   
   Call VecDestroy(resTemp,ierr);CHKERRQ(ierr)   
   Call MEF90_Finalize()


End Program  SimplePoissonNG
