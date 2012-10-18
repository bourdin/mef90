Module m_MEF90
#include "finclude/petscdef.h"
   Use petsc
   Use m_MEF_Ctx
   Use m_MEF_Diffusion
   Use m_MEF_Elements 
   Use m_MEF_EXO  
   Use m_MEF_LinAlg
   Use m_MEF_Parameters
   Use m_MEF_MassMatrixInterface
   Use m_MEF_Materials
   Use m_MEF_MPI
   Use m_MEF_Norm
   Use m_MEF_Sieve
   Use m_MEF_Utils

   Implicit NONE
   Public :: MEF90_Initialize
   Public :: MEF90_Finalize
   
   
Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_Initialize"
   Subroutine MEF90_Initialize(ierr)
      PetscInt,Intent(OUT)                               :: ierr
       
      Call PetscLogBegin(ierr);CHKERRQ(ierr)

      
      !!! Individual modules runtime initialization should be called here
      Call MEF90MPI_InitializePrivate(ierr);CHKERRQ(ierr)
      Call MEF90Materials_InitializePrivate(ierr);CHKERRQ(ierr)
      Call MEF90Ctx_InitializePrivate(ierr);CHKERRQ(ierr)
      
   End Subroutine MEF90_Initialize
   
   
#undef __FUNCT__
#define __FUNCT__ "MEF90_Finalize"
   Subroutine MEF90_Finalize(ierr)
      PetscInt,Intent(OUT)                   :: ierr
      
      Call MEF90MPI_FinalizePrivate(ierr);CHKERRQ(ierr)
      !Call PetscFinalize(ierr)
   End Subroutine MEF90_Finalize
End Module m_MEF90

