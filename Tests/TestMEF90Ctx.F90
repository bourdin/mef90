Program  TestMEF90Ctx
#include <finclude/petscdef.h>
#include <finclude/petscbagdef.h>
   Use m_MEF90
   Use petsc
   Implicit NONE   

   PetscErrorCode                      :: ierr
   PetscReal,Dimension(:),Pointer      :: time
   Type(MEF90Ctx_Type),target          :: MEF90Ctx
   Type(MEF90CtxGlobalOptions_Type)    :: MEF90GlobalOptions_default
   Type(DM),target                     :: Mesh
   PetscBool                           :: flg
   Character(len=MEF90_MXSTRLEN)       :: IOBuffer

   MEF90GlobalOptions_default%verbose           = 1
   !MEF90GlobalOptions_default%prefix            = '../TestMeshes/SquareNG-tri3'
   MEF90GlobalOptions_default%timeInterpolation = MEF90TimeInterpolation_linear
   MEF90GlobalOptions_default%timeMin           = 0.0_Kr
   MEF90GlobalOptions_default%timeMax           = 1.0_Kr
   MEF90GlobalOptions_default%timeNumStep       = 11
   MEF90GlobalOptions_default%fileFormat        = MEF90FileFormat_EXOSingle

   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   Call MEF90_Initialize(ierr)
   Call MEF90Ctx_Create(PETSC_COMM_WORLD,MEF90Ctx,MEF90GlobalOptions_default,ierr)

   Call MEF90Ctx_GetDMMeshEXO(MEF90Ctx,Mesh,ierr)
   Call DMView(Mesh,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)

   Call MEF90Ctx_OpenEXO(MEF90Ctx,Mesh,ierr)
   Call MEF90Ctx_GetTime(MEF90Ctx,time,ierr)

   Call MEF90Ctx_Destroy(MEF90Ctx,ierr);CHKERRQ(ierr)
   
   Call PetscOptionsHasName(PETSC_NULL_CHARACTER,"-verbose",flg,ierr);CHKERRQ(ierr)
   If (flg) Then
      Write(IOBuffer,*) "verbose is set\n"
   Else
      Write(IOBuffer,*) "verbose is NOT set\n"
   End If
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
   Call PetscOptionsHasName(PETSC_NULL_CHARACTER,"-time_min",flg,ierr);CHKERRQ(ierr)
   If (flg) Then
      Write(IOBuffer,*) "time_min is set\n"
   Else
      Write(IOBuffer,*) "time_min is NOT set\n"
   End If
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
   Call MEF90_Finalize(ierr)
   Call PetscFinalize()
End Program  TestMEF90Ctx
