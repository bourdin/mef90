Module m_MEF90
#include "finclude/petscdef.h"
   Use petsc
   Use m_MEF_Ctx
   Use m_MEF_DiffusionInterface
   Use m_MEF_ElasticityInterface
   Use m_MEF_Elements 
   Use m_MEF_EXO  
   Use m_MEF_LinAlg
   Use m_MEF_Parameters
   Use m_MEF_MassMatrixInterface
   Use m_MEF_Materials
   Use m_MEF_MPI
   Use m_MEF_Norm
   Use m_MEF_Utils

   Implicit NONE
   Public :: MEF90Initialize
   Public :: MEF90Finalize
   
   
Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90Initialize"
!!!
!!!  
!!!  MEF90Initialize:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90Initialize(ierr)
      PetscInt,Intent(OUT)                               :: ierr
       
      Call PetscLogBegin(ierr);CHKERRQ(ierr)

      
      !!! Individual modules runtime initialization should be called here
      Call MEF90MPIInitialize_Private(ierr);CHKERRQ(ierr)
      Call MEF90MaterialsInitialize_Private(ierr);CHKERRQ(ierr)
      Call MEF90CtxInitialize_Private(ierr);CHKERRQ(ierr)
   End Subroutine MEF90Initialize
   
   
#undef __FUNCT__
#define __FUNCT__ "MEF90Finalize"
!!!
!!!  
!!!  MEF90Finalize:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90Finalize(ierr)
      PetscInt,Intent(OUT)                   :: ierr
      
      Call MEF90MPIFinalize_Private(ierr);CHKERRQ(ierr)
   End Subroutine MEF90Finalize
End Module m_MEF90

