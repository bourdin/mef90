Module m_MEF90
#include "finclude/petscdef.h"
   Use petsc
   Use m_MEF_Ctx
   Use m_MEF_Diffusion
   Use m_MEF_Elements 
   Use m_MEF_EXO  
   Use m_MEF_LinAlg
   Use m_MEF_Parameters
   Use m_MEF_MassMatrix
   Use m_MEF_Materials
   Use m_MEF_MPI
   Use m_MEF_Norm
   Use m_MEF_Sieve
   Use m_MEF_Utils

   Public :: MEF90_Initialize
   Public :: MEF90_Finalize
   
Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_Initialize"
   Subroutine MEF90_Initialize()
      PetscInt                      :: ierr
       
      Call PetscInitialize(PETSC_NULL_CHARACTER,ierr);CHKERRQ(ierr)
      Call MEF90_MPIInitializePrivate(ierr);CHKERRQ(ierr)

      Call PetscLogBegin(ierr);CHKERRQ(ierr)
      
      !!! Individual modules runtime initialization should be called here
      Call MEF90_MaterialsInitializePrivate(ierr);CHKERRQ(ierr)
      Call MEF90Ctx_InitializePrivate(ierr);CHKERRQ(ierr)
   End Subroutine MEF90_Initialize
   
#undef __FUNCT__
#define __FUNCT__ "MEF90_Finalize"
   Subroutine MEF90_Finalize()
      PetscInt                         :: ierr
      
      Call MEF90_MPIFinalizePrivate(ierr);CHKERRQ(ierr)
      Call PetscFinalize(ierr)
   End Subroutine MEF90_Finalize
End Module m_MEF90

