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
    type(tIS)                           :: setIS,setPointIS
    PetscInt,Dimension(:),pointer       :: setID,pointID
    PetscInt                            :: labelSize,i
    Character(len=MEF90MXSTRLEN)        :: IOBuffer

    MEF90GlobalOptions_default%verbose           = 1
    MEF90GlobalOptions_default%dryrun            = PETSC_FALSE
    MEF90GlobalOptions_default%timeInterpolation = MEF90TimeInterpolation_linear
    MEF90GlobalOptions_default%timeMin           = 0.0_Kr
    MEF90GlobalOptions_default%timeMax           = 1.0_Kr
    MEF90GlobalOptions_default%timeNumStep       = 11
    MEF90GlobalOptions_default%timeSkip          = 0
    MEF90GlobalOptions_default%timeNumCycle      = 1

    PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
    
    Call MEF90Initialize(ierr)
    Call MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90GlobalOptions_default,ierr)
    
    PetscCallA(PetscPrintf(MEF90Ctx%Comm,MEF90Ctx%geometryfile//'\n',ierr))
    PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm,MEF90Ctx%geometryfile,PETSC_NULL_CHARACTER,interpolate,dm,ierr))
    PetscCallA(DMPlexDistributeSetDefault(dm,PETSC_FALSE,ierr))
    PetscCallA(DMSetFromOptions(dm,ierr))
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-mef90dm_view",ierr))
    
    distribute: Block 
        Type(tDM),target                    :: dmDist
        PetscInt                            :: ovlp = 0_Ki
        If (MEF90Ctx%NumProcs > 1) Then
            PetscCallA(DMPlexDistribute(dm,ovlp,PETSC_NULL_SF,dmDist,ierr))
            PetscCallA(DMDestroy(dm,ierr))
            dm = dmDist
        End If
    End Block distribute
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-dm_view",ierr))

    Do i = 1, size(MEF90SetLabelName)
        PetscCallA(PetscPrintf(PETSC_COMM_WORLD,'=== '//trim(MEF90SetLabelName(i))//' ===\n',ierr))
        PetscCallA(DMGetLabelSize(dm,MEF90SetLabelName(i),labelSize,ierr))
        Write(IOBuffer,'("label size",I5,"\n")') labelSize
        ! Note that label size refers only to the local dm, not the whole mesh
        PetscCallA(PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
        PetscCallA(PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT,ierr))
        PetscCallA(DMGetLabelIdIS(dm,MEF90SetLabelName(i), setIS, ierr))
        If (setIS /= PETSC_NULL_IS) Then
            PetscCallA(ISGetIndicesF90(setIS,setID,ierr))
            Write(IOBuffer,*) trim(MEF90SetLabelName(i))//' ID: ', setID,'\n'
            PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
            Do set = 1, size(setID)
                PetscCallA(DMGetStratumIS(dm,MEF90SetLabelName(i),setID(set),setPointIS,ierr))
                If (setPointIS /= PETSC_NULL_IS) Then
                    PetscCallA(ISGetIndicesF90(setPointIS,pointID,ierr))
                    Write(*,*) '   points'
                    Write(*,*) pointID
                    PetscCallA(ISRestoreIndicesF90(setPointIS,pointID,ierr))
                End If ! setPointIS
                PetscCallA(ISDestroy(setPointIS,ierr))
            End Do
            PetscCallA(ISRestoreIndicesF90(setIS,setID,ierr))
        End If ! setIS
        PetscCallA(ISDestroy(setIS,ierr))
    End Do
    PetscCallA(DMDestroy(dm,ierr))    
    
    Call MEF90CtxDestroy(MEF90Ctx,ierr)   
    Call MEF90Finalize(ierr)
    Call PetscFinalize(ierr)
End Program  TestLabels
