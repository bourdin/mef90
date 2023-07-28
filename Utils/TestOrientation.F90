Program  TestOrientation
#include <petsc/finclude/petsc.h>
Use m_MEF90
Use petsc
Implicit NONE   
    
    PetscErrorCode                            :: ierr
    Type(MEF90Ctx_Type),target                :: MEF90Ctx
    Type(MEF90CtxGlobalOptions_Type)          :: MEF90GlobalOptions_default
    Type(tDM)                                 :: dm
    PetscBool                                 :: interpolate = PETSC_TRUE
    Character(len=MEF90MXSTRLEN)              :: IOBuffer

    PetscInt                                  :: dim,pStart,pEnd
    PetscInt                                  :: cStart,cEnd,c
    PetscReal,Dimension(:),Pointer            :: v0,BB,BBinv
    PetscReal                                 :: detBBinv

    PetscBool                                 :: reversedCells = PETSC_FALSE
    Type(MEF90CtxGlobalOptions_Type),pointer  :: MEF90GlobalOptions


    PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
    PetscCallA(MEF90Initialize(ierr))

    MEF90GlobalOptions_default%verbose           = 0
    MEF90GlobalOptions_default%dryrun            = PETSC_FALSE
    MEF90GlobalOptions_default%timeMin           = 0.0_Kr
    MEF90GlobalOptions_default%timeMax           = 1.0_Kr
    MEF90GlobalOptions_default%timeNumStep       = 11
    MEF90GlobalOptions_default%timeInterpolation = MEF90TimeInterpolation_linear
    MEF90GlobalOptions_default%timeSkip          = 0
    MEF90GlobalOptions_default%timeNumCycle      = 1
    MEF90GlobalOptions_default%elementFamily     = MEF90ElementFamilyLagrange
    MEF90GlobalOptions_default%elementOrder      = 1

    PetscCallA(MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90GlobalOptions_default,ierr))
    PetscCallA(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr))
    
    PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm,MEF90Ctx%geometryfile,PETSC_NULL_CHARACTER,interpolate,dm,ierr))
    PetscCallA(DMPlexDistributeSetDefault(dm,PETSC_FALSE,ierr))
    PetscCallA(DMSetFromOptions(dm,ierr))
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-dm_view",ierr))
    
    distribute: Block 
        Type(tDM),target                    :: dmDist
        PetscInt                            :: ovlp = 0_Ki
        If (MEF90Ctx%NumProcs > 1) Then
            PetscCallA(DMPlexDistribute(dm,ovlp,PETSC_NULL_SF,dmDist,ierr))
            PetscCallA(DMDestroy(dm,ierr))
            dm = dmDist
        End If
    End Block distribute
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-mef90dm_view",ierr))

    PetscCallA(DMGetDimension(dm,dim,ierr))
    PetscCallA(DMPlexGetChart(dm,pStart,pEnd,ierr))
    PetscCallA(DMPlexGetHeightStratum(dm,0_Ki,cStart,cEnd,ierr))
    Allocate(v0(2))
    Allocate(BB(4))
    Allocate(BBinv(4))
    Do c = cStart,cEnd-1
        PetscCall(DMPlexComputeCellGeometryAffineFEM(dm,c,v0,BB,BBinv,detBBinv,ierr))
        If (detBBinv <= 0.0_Kr) Then
            reversedCells = PETSC_TRUE
            If (MEF90GlobalOptions%verbose > 0) Then
                Write(IOBuffer,'("Reversed cell: ",I0," Jacobian : ",ES12.5,"\n")') c,detBBinv
                PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
            End If
        End If
    End Do
    If (reversedCells) Then
        Write(IOBuffer,*) "WARNING: this mesh has cells with negative Jacobian. \nRun with -verbose 1 to see a list of flipped cells\n"
        PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
    Else
        Write(IOBuffer,*) "No reversed cells\n"
    End If
    PetscCallA(DMDestroy(dm,ierr))
    DeAllocate(v0)
    DeAllocate(BB)
    DeAllocate(BBinv)

    Call MEF90CtxDestroy(MEF90Ctx,ierr)   
    Call MEF90Finalize(ierr)
    Call PetscFinalize(ierr)
End Program  TestOrientation
       
