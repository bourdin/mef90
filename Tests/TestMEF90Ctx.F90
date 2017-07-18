Program  TestMEF90Ctx
#include <petsc/finclude/petsc.h>
   Use m_MEF90
   Use petsc
   Implicit NONE   

   PetscErrorCode                      :: ierr
   PetscReal,Dimension(:),Pointer      :: time
   Type(MEF90Ctx_Type),target          :: MEF90Ctx
   Type(MEF90CtxGlobalOptions_Type)    :: MEF90GlobalOptions_default
   Type(tDM),target                    :: dm,dmDist
   PetscBool                           :: flg
   PetscBool                           :: interpolate = PETSC_FALSE
   Character(len=MEF90_MXSTRLEN)       :: IOBuffer

   MEF90GlobalOptions_default%verbose           = 1
   MEF90GlobalOptions_default%dryrun            = PETSC_FALSE
   MEF90GlobalOptions_default%timeInterpolation = MEF90TimeInterpolation_linear
   MEF90GlobalOptions_default%timeMin           = 0.0_Kr
   MEF90GlobalOptions_default%timeMax           = 1.0_Kr
   MEF90GlobalOptions_default%timeFrequency     = 0.0_Kr
   MEF90GlobalOptions_default%timeNumStep       = 11
   MEF90GlobalOptions_default%fileFormat        = MEF90FileFormat_EXOSingle

   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   Call MEF90Initialize(ierr)
   Call MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90GlobalOptions_default,ierr)

   Call DMPlexCreateFromFile(PETSC_COMM_WORLD,MEF90Ctx%inputmesh,interpolate,dm,ierr);CHKERRQ(ierr);
   Call DMPlexDistribute(dm,0,PETSC_NULL_SF,dmDist,ierr);CHKERRQ(ierr)

   !!!
   !!! I am not sure everybody would approve of this...
   !!!
   if (dmDist%v /= -1) then
      call DMDestroy(dm,ierr)
      dm%v = dmDist%v
   end if

   Call DMView(dm,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)

   !Call MEF90CtxOpenEXO(MEF90Ctx,Mesh,ierr)
   !Call MEF90CtxGetTime(MEF90Ctx,time,ierr)

   Call MEF90CtxDestroy(MEF90Ctx,ierr)
   
   Call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-verbose",flg,ierr)
   If (flg) Then
      Write(IOBuffer,*) "verbose is set\n"
   Else
      Write(IOBuffer,*) "verbose is NOT set\n"
   End If
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr)
   Call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-time_min",flg,ierr)
   If (flg) Then
      Write(IOBuffer,*) "time_min is set\n"
   Else
      Write(IOBuffer,*) "time_min is NOT set\n"
   End If
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr)
   Call MEF90Finalize(ierr)
   Call PetscFinalize()
End Program  TestMEF90Ctx
