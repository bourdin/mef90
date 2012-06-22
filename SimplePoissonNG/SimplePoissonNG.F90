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
   
   

   Call MEF90_Initialize()
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
   
   Call SimplePoissonInitialize(AppCtx,prefix,ierr);CHKERRQ(ierr)
   
   
   Call SimplePoissonFinalize(AppCtx,ierr);CHKERRQ(ierr)
   Call MEF90_Finalize()


Contains
#undef __FUNCT__
#define __FUNCT__ "SimplePoissonInitialize"
   Subroutine SimplePoissonInitialize(PoissonCtx,prefix,ierr)
      Type(PoissonCtx_Type),intent(INOUT)       :: PoissonCtx
      Character(len=*),intent(IN)               :: prefix
      PetscErrorCode,Intent(OUT)                :: ierr

      Call PetscLogBegin(ierr);CHKERRQ(ierr)
      Call m_Poisson_Initialize(ierr);CHKERRQ(ierr)
      Call PoissonCtxInit(PoissonCtx,prefix,ierr);CHKERRQ(ierr)
   End Subroutine SimplePoissonInitialize   

#undef __FUNCT__
#define __FUNCT__ "SimplePoissonFinalize"
   Subroutine SimplePoissonFinalize(PoissonCtx,ierr)
      Type(PoissonCtx_Type),intent(IN)          :: PoissonCtx
      PetscErrorCode,Intent(OUT)                :: ierr
      
      ierr = 0
   End Subroutine SimplePoissonFinalize   
End Program  SimplePoissonNG
