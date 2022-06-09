Program  TestLabels
#include <petsc/finclude/petsc.h>
Use m_MEF90
Use petsc
Implicit NONE   
    
    PetscErrorCode                      :: ierr
    Type(MEF90Ctx_Type),target          :: MEF90Ctx
    Type(MEF90CtxGlobalOptions_Type)    :: MEF90GlobalOptions_default
    Type(tDM),target                    :: dm,dmDist
    PetscBool                           :: interpolate = PETSC_TRUE
    Character(len=MEF90_MXSTRLEN)       :: setType

    PetscInt                            :: set
    type(tIS)                           :: setIS,pointIS
    PetscInt,Dimension(:),pointer       :: setID,pointID
    PetscInt                            :: i

    MEF90GlobalOptions_default%verbose           = 1
    MEF90GlobalOptions_default%dryrun            = PETSC_FALSE
    MEF90GlobalOptions_default%timeInterpolation = MEF90TimeInterpolation_linear
    MEF90GlobalOptions_default%timeMin           = 0.0_Kr
    MEF90GlobalOptions_default%timeMax           = 1.0_Kr
    MEF90GlobalOptions_default%timeNumStep       = 11
    MEF90GlobalOptions_default%timeSkip          = 0
    MEF90GlobalOptions_default%timeNumCycle      = 1
    MEF90GlobalOptions_default%fileFormat        = MEF90FileFormat_EXOSingle

    PetscCall(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
    
    Call MEF90Initialize(ierr)
    Call MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90GlobalOptions_default,ierr)
    
    PetscCall(PetscPrintf(MEF90Ctx%Comm,MEF90Ctx%geometryfile//'\n',ierr))
    PetscCall(DMPlexCreateFromFile(MEF90Ctx%Comm,MEF90Ctx%geometryfile,PETSC_NULL_CHARACTER,interpolate,dm,ierr))
    PetscCall(DMPlexDistributeSetDefault(dm,PETSC_FALSE,ierr))
    PetscCall(DMSetFromOptions(dm,ierr))
    PetscCall(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-dm_view",ierr))
    
    PetscCall(DMPlexDistribute(dm,0,PETSC_NULL_SF,dmDist,ierr))
    if (MEF90Ctx%NumProcs > 1) then
        PetscCall(DMDestroy(dm,ierr))
        dm = dmDist
    end if
    PetscCall(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-dm_view",ierr))

    Do i = 1, size(MEF90_DMPlexSetTypes)
        PetscCall(PetscPrintf(PETSC_COMM_WORLD,'=== '//trim(MEF90_DMPlexSetTypes(i))//' ===\n',ierr))
        PetscCall(DMGetLabelIdIS(dm,MEF90_DMPlexSetTypes(i), SetIS, ierr))
        PetscCall(ISGetIndicesF90(SetIS,setID,ierr))
        Write(*,*) trim(MEF90_DMPlexSetTypes(i))//' ID: ', setID
        Do set = 1, size(setID)
            PetscCall(DMGetStratumIS(dm,MEF90_DMPlexSetTypes(i),setID(set),pointIS,ierr))
            PetscCall(ISGetIndicesF90(pointIS,pointID,ierr))
            Write(*,*) '   points', pointID
            PetscCall(ISRestoreIndicesF90(pointIS,pointID,ierr))
            PetscCall(ISDestroy(pointIS,ierr))
        End Do
        PetscCall(ISRestoreIndicesF90(SetIS,setID,ierr))
        PetscCall(ISDestroy(SetIS,ierr))
    End Do
    PetscCall(DMDestroy(dm,ierr))
    PetscCall(PetscFinalize(ierr))
    
    Call MEF90CtxDestroy(MEF90Ctx,ierr)   
    Call MEF90Finalize(ierr)
    Call PetscFinalize(ierr)
End Program  TestLabels
