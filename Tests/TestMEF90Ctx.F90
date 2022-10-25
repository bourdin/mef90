Program  TestMEF90Ctx
#include <petsc/finclude/petsc.h>
   Use m_MEF90
   Use petsc
   Implicit NONE   

   PetscErrorCode                      :: ierr
   !PetscReal,Dimension(:),Pointer      :: time
   Type(MEF90Ctx_Type),target          :: MEF90Ctx
   Type(MEF90CtxGlobalOptions_Type)    :: MEF90GlobalOptions_default
   Type(tDM),target                    :: dm
   PetscBool                           :: flg
   PetscBool                           :: interpolate = PETSC_FALSE
   Character(len=MEF90MXSTRLEN)       :: IOBuffer

   MEF90GlobalOptions_default%verbose           = 1
   MEF90GlobalOptions_default%dryrun            = PETSC_FALSE
   MEF90GlobalOptions_default%timeInterpolation = MEF90TimeInterpolation_linear
   MEF90GlobalOptions_default%timeMin           = 0.0_Kr
   MEF90GlobalOptions_default%timeMax           = 1.0_Kr
   MEF90GlobalOptions_default%timeNumStep       = 11
   MEF90GlobalOptions_default%timeSkip          = 0
   MEF90GlobalOptions_default%timeNumCycle      = 1
   MEF90GlobalOptions_default%elementFamily     = MEF90ElementFamilyLagrange
   MEF90GlobalOptions_default%elementOrder      = 1

   PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
   PetscCallA(MEF90Initialize(ierr))
   PetscCallA(MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90GlobalOptions_default,ierr))

   PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm,MEF90Ctx%geometryfile,PETSC_NULL_CHARACTER,interpolate,dm,ierr))
   PetscCallA(DMPlexDistributeSetDefault(dm,PETSC_FALSE,ierr))
   PetscCallA(DMSetFromOptions(dm,ierr))
   PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-dm_view",ierr))

   distribute: Block 
      Type(tDM),target                    :: dmDist
      If (MEF90Ctx%NumProcs > 1) Then
         PetscCallA(DMPlexDistribute(dm,0,PETSC_NULL_SF,dmDist,ierr))
         PetscCallA(DMDestroy(dm,ierr))
         dm = dmDist
      End If
   End Block distribute
   PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-mef90dm_view",ierr))

   ! PetscCallA(MEF90CtxOpenEXO(MEF90Ctx,dm,ierr))
   ! PetscCallA(MEF90CtxGetTime(MEF90Ctx,time,ierr))

   PetscCallA(MEF90CtxDestroy(MEF90Ctx,ierr))
   
   PetscCallA(PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-verbose",flg,ierr))
   If (flg) Then
      Write(IOBuffer,*) "verbose is set\n"
   Else
      Write(IOBuffer,*) "verbose is NOT set\n"
   End If
   PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
   PetscCallA(PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-time_min",flg,ierr))
   If (flg) Then
      Write(IOBuffer,*) "time_min is set\n"
   Else
      Write(IOBuffer,*) "time_min is NOT set\n"
   End If
   PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))

   PetscCallA(MEF90Finalize(ierr))
   PetscCallA(PetscFinalize(ierr))
End Program  TestMEF90Ctx
