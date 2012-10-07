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

   Implicit NONE
   Public :: MEF90_Initialize
   Public :: MEF90_Finalize
   Public :: MEF90CtxBag
   
   PetscBag                      :: MEF90CtxBag 
   
Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_Initialize"
   Subroutine MEF90_Initialize(MEF90Ctx,default,ierr)
      Type(MEF90Ctx_type),Intent(OUT),pointer   :: MEF90Ctx
      Type(MEF90Ctx_Type),Intent(IN)            :: default
      PetscInt,Intent(OUT)                      :: ierr
       
      !Call PetscInitialize(PETSC_NULL_CHARACTER,ierr);CHKERRQ(ierr)
      Call PetscLogBegin(ierr);CHKERRQ(ierr)

      
      !!! Individual modules runtime initialization should be called here
      Call MEF90MPI_InitializePrivate(ierr);CHKERRQ(ierr)
      Call MEF90Materials_InitializePrivate(ierr);CHKERRQ(ierr)
      Call MEF90Ctx_InitializePrivate(ierr);CHKERRQ(ierr)
      
      Call PetscBagCreate(PETSC_COMM_WORLD,sizeofMEF90Ctx,MEF90CtxBag,ierr)
      Call PetscBagRegisterMEF90Ctx(MEF90CtxBag,"MEF90 Global Context",PETSC_NULL_CHARACTER,default,ierr)
      Call PetscBagGetDataMEF90Ctx(MEF90CtxBag,MEF90Ctx,ierr);CHKERRQ(ierr)  
   End Subroutine MEF90_Initialize
   
   
#undef __FUNCT__
#define __FUNCT__ "MEF90_Finalize"
   Subroutine MEF90_Finalize()
      PetscInt                   :: ierr
      
      Call PetscBagDestroy(MEF90CtxBag,ierr)
      Call MEF90MPI_FinalizePrivate(ierr);CHKERRQ(ierr)
      !Call PetscFinalize(ierr)
   End Subroutine MEF90_Finalize
End Module m_MEF90

