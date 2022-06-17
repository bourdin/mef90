Program  TestLabels
#include <petsc/finclude/petsc.h>
Use m_MEF90
Use petsc
Implicit NONE   
    
    PetscErrorCode                      :: ierr
    Type(MEF90Ctx_Type),target          :: MEF90Ctx
    Type(MEF90CtxGlobalOptions_Type)    :: MEF90GlobalOptions_default
    Type(tDM),target                    :: dm
    PetscBool                           :: interpolate = PETSC_TRUE

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

    PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
    
    Call MEF90Initialize(ierr)
    Call MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90GlobalOptions_default,ierr)
    
    PetscCallA(PetscPrintf(MEF90Ctx%Comm,MEF90Ctx%geometryfile//'\n',ierr))
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
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-dm_view",ierr))

    Do i = 1, size(MEF90_DMPlexSetLabelName)
        PetscCallA(PetscPrintf(PETSC_COMM_WORLD,'=== '//trim(MEF90_DMPlexSetLabelName(i))//' ===\n',ierr))
        PetscCallA(DMGetLabelIdIS(dm,MEF90_DMPlexSetLabelName(i), SetIS, ierr))
        PetscCallA(ISGetIndicesF90(SetIS,setID,ierr))
        Write(*,*) trim(MEF90_DMPlexSetLabelName(i))//' ID: ', setID
        Do set = 1, size(setID)
            PetscCallA(DMGetStratumIS(dm,MEF90_DMPlexSetLabelName(i),setID(set),pointIS,ierr))
            PetscCallA(ISGetIndicesF90(pointIS,pointID,ierr))
            Write(*,*) '   points', pointID
            PetscCallA(ISRestoreIndicesF90(pointIS,pointID,ierr))
            PetscCallA(ISDestroy(pointIS,ierr))
        End Do
        PetscCallA(ISRestoreIndicesF90(SetIS,setID,ierr))
        PetscCallA(ISDestroy(SetIS,ierr))
    End Do
    PetscCallA(DMDestroy(dm,ierr))
    PetscCallA(PetscFinalize(ierr))
    
    Call MEF90CtxDestroy(MEF90Ctx,ierr)   
    Call MEF90Finalize(ierr)
    Call PetscFinalize(ierr)
End Program  TestLabels
