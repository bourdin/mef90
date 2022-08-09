Program  TestVec
#include <petsc/finclude/petsc.h>
Use petsc
Use m_MEF90
Use m_MEF90_HeatXfer
Implicit NONE   
    
    PetscErrorCode                                     :: ierr
    Type(MEF90Ctx_Type),target                         :: MEF90Ctx
    Type(MEF90CtxGlobalOptions_Type)                   :: MEF90GlobalOptions_default
    Type(MEF90CtxGlobalOptions_Type),Pointer           :: MEF90GlobalOptions
    Type(tDM)                                          :: dm
    PetscBool                                          :: interpolate = PETSC_TRUE
    Character(len=MEF90MXSTRLEN)                       :: name

    Type(tVec)                                         :: v
    PetscBool,Dimension(:,:),Pointer                   :: cellSetBC,faceSetBC,vertexSetBC
    Type(MEF90HeatXferCellSetOptions_Type),pointer     :: cellSetOptions,facesetOptions
    Type(MEF90HeatXferVertexSetOptions_Type),pointer   :: vertexSetOptions
    PetscInt                                           :: set
    Type(tIS)                                          :: setIS
    PetscInt,Dimension(:),Pointer                      :: setID
    PetscEnum                                          :: setType

    Type(MEF90HeatXferCtx_Type)                        :: MEF90HeatXferCtx
    Type(MEF90HeatXferGlobalOptions_Type),Pointer      :: MEF90HeatXferGlobalOptions
    Type(MEF90HeatXferGlobalOptions_Type),Parameter    :: MEF90HeatXferDefaultGlobalOptions = MEF90HeatXferGlobalOptions_Type( &
                                                          MEF90HeatXFer_timeSteppingTypeSteadyState, & ! timeSteppingType
                                                          PETSC_FALSE,         & ! addNullSpace
                                                          0.,                  & ! initialTemperature
                                                          MEF90Scaling_Linear, & ! boundaryTempScaling
                                                          MEF90Scaling_Linear, & ! externalTempScaling
                                                          MEF90Scaling_Linear)   ! fluxScaling
    Type(MEF90HeatXferCellSetOptions_Type),Parameter   :: MEF90HeatXferDefaultCellSetOptions = MEF90HeatXferCellSetOptions_Type( &
                                                          0.0_Kr,        & ! flux
                                                          0.0_Kr,        & ! surfaceThermalConductivity
                                                          0.0_Kr,        & ! externalTemp
                                                          PETSC_FALSE,   & ! Has BC
                                                          0.0_Kr,        & ! boundaryTemp
                                                          [0.0_Kr,0.0_Kr,0.0_Kr]) ! AdvectionVelocity
                                                          
    Type(MEF90HeatXferVertexSetOptions_Type),Parameter :: MEF90HeatXferDefaultVertexSetOptions = MEF90HeatXferVertexSetOptions_Type( &
                                                          PETSC_FALSE,   & ! Has BC
                                                          0.0_Kr)          ! boundaryTemp
 
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
    
    Call MEF90Initialize(ierr)
    Call MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90GlobalOptions_default,ierr)
    PetscCallA(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr))
    
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
    
    !!! Create HeatXfer context, get all HeatXfer options
    Call MEF90HeatXferCtxCreate(MEF90HeatXferCtx,dm,MEF90Ctx,ierr);CHKERRQ(ierr)
    Call MEF90HeatXferCtxSetFromOptions(MEF90HeatXferCtx,PETSC_NULL_CHARACTER,MEF90HeatXferDefaultGlobalOptions, &
                                        MEF90HeatXferDefaultCellSetOptions,MEF90HeatXferDefaultVertexSetOptions,ierr)
  
    setType = MEF90CellSetType
    PetscCall(DMGetLabelIdIS(dm,MEF90SetLabelName(setType),setIS,ierr))
    PetscCall(MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr))
    PetscCall(ISGetIndicesF90(setIS,setID,ierr))
Write(*,*) 'cell sets ', setID
    Allocate(cellSetBC(size(setID),1))
    Do set = 1,size(setID)
        Call PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
        cellSetBC(set,:) = cellSetOptions%Has_BC
    End Do ! set
    PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
    PetscCall(ISDestroy(setIS,ierr))

    setType = MEF90FaceSetType
    PetscCall(DMGetLabelIdIS(dm,MEF90SetLabelName(setType),setIS,ierr))
    PetscCall(MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr))
    PetscCall(ISGetIndicesF90(setIS,setID,ierr))
Write(*,*) 'face sets ', setID
    Allocate(faceSetBC(size(setID),1))
    Do set = 1,size(setID)
        Call PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%FaceSetOptionsBag(set),faceSetOptions,ierr);CHKERRQ(ierr)
        faceSetBC(set,:) = faceSetOptions%Has_BC
    End Do ! set
    PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
    PetscCall(ISDestroy(setIS,ierr))

    setType = MEF90VertexSetType
    PetscCall(DMGetLabelIdIS(dm,MEF90SetLabelName(setType),setIS,ierr))
    PetscCall(MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr))
    PetscCall(ISGetIndicesF90(setIS,setID,ierr))
Write(*,*) 'vertex sets ', setID
    Allocate(vertexSetBC(size(setID),1))
    Do set = 1,size(setID)
        Call PetscBagGetDataMEF90HeatXferCtxVertexSetOptions(MEF90HeatXferCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
        vertexSetBC(set,:) = vertexSetOptions%Has_BC
    End Do ! set
    PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
    PetscCall(ISDestroy(setIS,ierr))

    write(*,*) 'cellSetBC:   ',cellSetBC
    write(*,*) 'faceSetBC:   ',faceSetBC
    write(*,*) 'vertexSetBC: ',vertexSetBC

    name = "Temperature"
    PetscCallA(MEF90VecCreate(dm,MEF90GlobalOptions%elementFamily,MEF90GlobalOptions%elementOrder,1,CellSetBC,FaceSetBC,VertexSetBC,name,V,ierr))
    !PetscCallA(VecViewFromOptions(V,PETSC_NULL_OPTIONS,"-vec_view",ierr))

    DeAllocate(cellSetBC)
    DeAllocate(faceSetBC)
    DeAllocate(vertexSetBC)
    PetscCallA(DMDestroy(dm,ierr))
    Call MEF90CtxDestroy(MEF90Ctx,ierr)   
    Call MEF90Finalize(ierr)
    Call PetscFinalize(ierr)
End Program  TestVec