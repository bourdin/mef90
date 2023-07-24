Program  gmsh2exo
#include <petsc/finclude/petsc.h>
    Use m_MEF90
    Use petsc
    Implicit NONE   

    PetscErrorCode                      :: ierr
    Type(MEF90Ctx_Type),target          :: MEF90Ctx
    Type(MEF90CtxGlobalOptions_Type)    :: MEF90GlobalOptions_default
    Type(tDM),target                    :: dm
    PetscBool                           :: interpolate = PETSC_FALSE

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
    PetscCallA(DMSetFromOptions(dm,ierr))
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-dm_view",ierr))

    PetscCallA(PetscViewerDestroy(MEF90Ctx%resultViewer,ierr))
    PetscCallA(MEF90CtxOpenEXO(MEF90Ctx,MEF90Ctx%resultViewer,FILE_MODE_WRITE,ierr))
    PetscCallA(MEF90EXODMView(dm,MEF90Ctx%resultViewer,MEF90GlobalOptions_default%elementOrder,ierr))
    PetscCallA(PetscViewerDestroy(MEF90Ctx%resultViewer,ierr))
    PetscCallA(MEF90CtxDestroy(MEF90Ctx,ierr))
    PetscCallA(MEF90Finalize(ierr))
    PetscCallA(PetscFinalize(ierr))
End Program gmsh2exo 