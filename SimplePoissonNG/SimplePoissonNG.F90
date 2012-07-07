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
   Type(Vec)                                    :: solTemp
   
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
   Call DMMeshCreateMatrix(AppCtx%mesh,AppCtx%Section,MATMPIAIJ,matTemp,iErr);CHKERRQ(iErr)
   Call MatSetFromOptions(matTemp,ierr);CHKERRQ(ierr)
   !!!
   !!! Create SNES, 
   !!!
   Call SNESCreate(PETSC_COMM_WORLD,snesTemp,ierr);CHKERRQ(ierr)
   Call SNESSetDM(snesTemp,AppCtx%mesh,ierr);CHKERRQ(ierr)
   Call SNESSetOptionsPrefix(snesTemp,'temp_',ierr);CHKERRQ(ierr)
   Call SNESSetFromOptions(snesTemp,ierr);CHKERRQ(ierr)
   If (verbose > 0) Then
      Call SNESView(snesTemp,PETSC_VIEWER_STDOUT_WORLD,ierr)
   End If

   Call SimplePoissonNGBilinearFormAssembly(snesTemp,solTemp,matTemp,matTemp,matflg,AppCtx,ierr)   
   Call MatView(matTemp,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)   
   
      
   Call PoissonCtxDestroy(AppCtx,ierr);CHKERRQ(ierr)
   Call MEF90_Finalize()


End Program  SimplePoissonNG
