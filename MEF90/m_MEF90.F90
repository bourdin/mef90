Module m_MEF90
#include "finclude/petscdef.h"
   Use petsc
   Use m_MEF_Parameters
   Use m_MEF_LinAlg
   Use m_MEF_Utils
   Use m_MEF_MPI
   Use m_MEF_Elements 
   Use m_MEF_Norm
   Use m_MEF_EXO  
   Use m_MEF_Sieve
   Use m_MEF_Materials

   Public :: MEF90_Initialize
   Public :: MEF90_Finalize
   
Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_Initialize"
   Subroutine MEF90_Initialize()
      PetscInt                      :: ierr
       
      Call PetscInitialize(PETSC_NULL_CHARACTER,ierr);CHKERRQ(ierr)
      Call MPIType_Initialize()
      Call MPI_COMM_RANK(MPI_COMM_WORLD,MEF90_MyRank,ierr)
      Call MPI_COMM_SIZE(MPI_COMM_WORLD,MEF90_NumProcs,ierr)
      
      Call MEF90_MaterialsInitialize(ierr);CHKERRQ(ierr)
   End Subroutine MEF90_Initialize
   
#undef __FUNCT__
#define __FUNCT__ "MEF90_Finalize"
   Subroutine MEF90_Finalize()
      PetscInt                         :: ierr
      
      Call MPIType_Finalize()
      Call PetscFinalize(ierr)
   End Subroutine MEF90_Finalize
End Module m_MEF90

