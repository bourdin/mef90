#include "SimplePoisson.inc"
Module m_Poisson
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_PoissonGlobalProperties
   Use M_POISSONCELLSETPROPERTIES
   Use m_PoissonVertexSetProperties
   Implicit NONE
   Private
  
   !Public :: MEF90Ctx_Type
    
   Public :: m_Poisson_Initialize
   !Public :: PoissonCtxCreate
   Public :: PoissonCtxDestroy
   
   PetscSizeT,protected,public    :: sizeofPoissonCtx
   
Contains
#undef __FUNCT__
#define __FUNCT__ "m_Poisson_Initialize"
   Subroutine m_Poisson_Initialize(ierr)
      PetscErrorCode,intent(OUT)          :: ierr

      Type(MEF90Ctx_Type),Target          :: PoissonCtx
      character(len=1),pointer            :: dummychar(:)
      PetscSizeT                          :: sizeofchar
   
      Call PoissonGlobalPropertiesInitialize(ierr)
      Call PoissonCellSetPropertiesInitialize(ierr)
      Call PoissonVertexSetPropertiesInitialize(ierr)
      
      Call PetscDataTypeGetSize(PETSC_CHAR,sizeofchar,ierr)
      sizeofPoissonCtx = size(transfer(PoissonCtx,dummychar))*sizeofchar
   End Subroutine m_Poisson_Initialize
   
   
#undef __FUNCT__
#define __FUNCT__ "PoissonCtxDestroy"
!!!
!!!  
!!!  PoissonCtxDestroy:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
Subroutine PoissonCtxDestroy(PoissonCtx,snesTemp,ierr)
   Type(MEF90Ctx_Type)                                      :: PoissonCtx
   Type(SNES),Intent(IN)                                    :: snesTemp
   PetscErrorCode,Intent(OUT)                               :: ierr

   Type(IS)                                                 :: setIS   
   Type(DM)                                                 :: mesh
   PetscInt                                                 :: e,set,nset
   
   Call SNESGetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)

   Call DMmeshGetLabelIdIS(mesh,'Cell Sets',SetIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr) 
   Call ISGetLocalSize(setIS,nset,ierr);CHKERRQ(ierr)
   Do set = 1, nset
      Call PetscBagDestroy(PoissonCtx%CellSetPropertiesBag(set),ierr);CHKERRQ(ierr)
      Call PetscBagDestroy(PoissonCtx%MaterialPropertiesBag(set),ierr);CHKERRQ(ierr)
   End Do
   DeAllocate(PoissonCtx%CellSetPropertiesBag,stat=ierr)
   DeAllocate(PoissonCtx%MaterialPropertiesBag,stat=ierr)
   Call ISDestroy(setIS,ierr);CHKERRQ(ierr)

   Call DMmeshGetLabelIdIS(mesh,'Vertex Sets',SetIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr) 
   Call ISGetLocalSize(setIS,nset,ierr);CHKERRQ(ierr)
   Do set = 1, nset
      Call PetscBagDestroy(PoissonCtx%VertexSetPropertiesBag(set),ierr);CHKERRQ(ierr)
   End Do
   DeAllocate(PoissonCtx%VertexSetPropertiesBag)
   Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
   
   Call PetscBagDestroy(PoissonCtx%GlobalPropertiesBag,ierr);CHKERRQ(ierr)
End Subroutine PoissonCtxDestroy 

End Module m_Poisson
