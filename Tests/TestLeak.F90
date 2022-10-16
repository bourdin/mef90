Program TestNSP
#include "petsc/finclude/petsc.h"
   Use m_MEF90
   Use petsc
   Implicit NONE   

   PetscErrorCode                      :: ierr
   Type(DM),target                     :: Mesh,MeshClone
   Character(len=MEF90MXSTRLEN)       :: IOBuffer
   Type(SectionReal)                   :: defaultSection1,defaultSection2,copySection
   PetscBool                           :: flg
   PetscInt                            :: i,dim
   PetscReal                           :: val
   Type(MEF90Ctx_Type)                 :: MEF90Ctx
   PetscReal,Dimension(:,:),Pointer    :: Coord      
   Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions
   Type(MEF90CtxGlobalOptions_Type),Parameter         :: MEF90DefaultGlobalOptions = MEF90CtxGlobalOptions_Type( &
                                                         1,                             & ! verbose
                                                         PETSC_FALSE,                   & ! validate
                                                         MEF90TimeInterpolation_linear, & ! timeInterpolation
                                                         0.0_Kr,                        & ! timeMin
                                                         1.0_Kr,                        & ! timeMax
                                                         11,                            & ! timeNumStep
                                                         MEF90FileFormat_EXOSingle,     & ! fileFormat
                                                         1.0_Kr)                         ! frequency



   !!! Initialize MEF90
   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   Call MEF90Initialize(ierr)

   !!! Get all MEF90-wide options
   Call MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90DefaultGlobalOptions,ierr);CHKERRQ(ierr)

   !!! Get DM from mesh
   Call MEF90CtxGetDMMeshEXO(MEF90Ctx,Mesh,ierr);CHKERRQ(ierr)
   Call DMMeshGetDimension(Mesh,dim,ierr);CHKERRQ(ierr)
   Call DMSetBlockSize(Mesh,dim,ierr);CHKERRQ(ierr)
   Call DMMeshSetMaxDof(Mesh,dim,ierr);CHKERRQ(ierr) 
   Call DMMeshGetVertexSectionReal(Mesh,"default",dim,defaultSection1,ierr);CHKERRQ(ierr)
   Call DMMeshSetSectionReal(Mesh,"default",defaultSection1,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(defaultSection1,ierr);CHKERRQ(ierr)

   !Do i = 1, 1000 
   !      if (mod(i,100) == 0) Then
   !         write(*,*) i
   !      end if
   !   Call DMMeshGetSectionReal(Mesh,'default',defaultSection1,ierr);CHKERRQ(ierr)
   !   Call SectionRealDuplicate(defaultSection1,copySection,ierr);CHKERRQ(ierr)
   !   val = i + 1.234
   !   Call SectionRealSet(defaultSection1,val,ierr);CHKERRQ(ierr)
   !   Call SectionRealDestroy(defaultSection1,ierr);CHKERRQ(ierr)
   !   Call SectionRealDestroy(copySection,ierr);CHKERRQ(ierr)
   !End Do

   !Do i = 1, 1000
   !      if (mod(i,100) == 0) Then
   !         write(*,*) i
   !      end if
   !      Call DMMeshGetCoordinatesF90(mesh,Coord,ierr);CHKERRQ(ierr)
   !      val = sum(Coord) * (i-1.)
   !      Call DMMeshRestoreCoordinatesF90(mesh,Coord,ierr);CHKERRQ(ierr)
   !End Do


   Do i = 1, 10000000
      Call DMMeshGetSectionReal(Mesh,'default',defaultSection1,ierr);CHKERRQ(ierr)
      Call DMMeshGetSectionReal(Mesh,'default',defaultSection2,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(defaultSection1,copySection,ierr);CHKERRQ(ierr)
      !Call SectionRealDestroy(defaultSection1,ierr);CHKERRQ(ierr)
      print*, defaultSection1%V, defaultSection2%V, copySection%V
      Call SectionRealDestroy(defaultSection2,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(copySection,ierr);CHKERRQ(ierr)
   End Do

   Call DMDestroy(Mesh,ierr);CHKERRQ(ierr)
   Call MEF90CtxDestroy(MEF90Ctx,ierr);CHKERRQ(ierr)   
   Call MEF90Finalize(ierr)
   Call PetscFinalize()
End Program TestNSP
