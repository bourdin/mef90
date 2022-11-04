Program  TestMassMatrix
#include <petsc/finclude/petsc.h>
    Use m_MEF90
    Use petsc
    Implicit NONE   
        
    PetscErrorCode                                 :: ierr
    Type(MEF90Ctx_Type),target                     :: MEF90Ctx
    Type(MEF90CtxGlobalOptions_Type)               :: MEF90GlobalOptions_default
    Type(tDM)                                      :: dm
    PetscBool                                      :: interpolate = PETSC_TRUE
    Character(len=MEF90MXSTRLEN)                   :: IOBuffer,name
    PetscEnum                                      :: setType
    PetscInt                                       :: set
    type(tIS)                                      :: setIS
    PetscInt,Dimension(:),pointer                  :: setID
    PetscInt                                       :: dim
    Type(tVec)                                     :: U,U0
    Type(tDM)                                      :: dmU,dmU0
    Type(tMat)                                     :: M,M0
    Type(MEF90CtxGlobalOptions_Type),pointer       :: MEF90GlobalOptions
    Type(MEF90Element2DScal),dimension(:),Pointer  :: elem2D
    Type(MEF90Element3DScal),dimension(:),Pointer  :: elem3D
    PetscInt                                       :: quadratureOrder
    Type(tIS)                                      :: setPointIS
    PetscInt,Dimension(:),Pointer                  :: setPointID
    Type(MEF90ElementType)                         :: elementType
    DMPolytopeType                                 :: cellType,faceType
    PetscReal                                      :: myL2NormSet,L2NormSet,L2Norm

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
    distribute: Block 
        Type(tDM),target                    :: dmDist
        PetscInt                            :: ovlp = 0
        If (MEF90Ctx%NumProcs > 1) Then
            PetscCallA(DMPlexDistribute(dm,ovlp,PETSC_NULL_SF,dmDist,ierr))
            PetscCallA(DMDestroy(dm,ierr))
            dm = dmDist
        End If
    End Block distribute
    PetscCallA(DMGetDimension(dm,dim,ierr))
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-mef90dm_view",ierr))
    
    BoundingBlock: Block
        PetscReal,dimension(3)         :: bbmin,bbmax

        PetscCallA(DMGetBoundingBox(dm,bbmin,bbmax,ierr))
        Write(IOBUffer,'("DM Bonding box ",(*(ES12.3,"  ")))') bbmin(1:dim),bbmax(1:dim)
        PetscCallA(PetscPrintf(PETSC_COMM_WORLD,trim(IOBuffer)//'\n',ierr))
    End Block BoundingBlock
        
    name = "U"
    PetscCallA(MEF90CreateLocalVector(dm,MEF90GlobalOptions%elementFamily,MEF90GlobalOptions%elementOrder,1_Ki,name,U,ierr))
    PetscCallA(VecSet(U,1.0_Kr,ierr))
    name = "U0"
    PetscCallA(MEF90CreateBoundaryLocalVector(dm,MEF90GlobalOptions%elementFamily,MEF90GlobalOptions%elementOrder,1_Ki,name,U0,ierr))
    PetscCallA(VecSet(U0,1.0_Kr,ierr))

    PetscCallA(VecGetDM(U,dmU,ierr))
    PetscCallA(VecGetDM(U0,dmU0,ierr))

    PetscCallA(DMCreateMatrix(dmU,M,iErr))
    PetscCallA(DMCreateMatrix(dmU0,M0,iErr))

    setType = MEF90CellSetType
    PetscCallA(DMGetLabelIdIS(dmU,MEF90SetLabelName(setType),setIS,ierr))
    PetscCall(MEF90ISAllGatherMerge(MEF90Ctx%Comm,setIS,ierr))
    L2Norm    = 0.0_Kr
    If (setIS /= PETSC_NULL_IS) Then
        PetscCallA(ISGetIndicesF90(setIS,setID,ierr))
        Do set = 1,size(setID)
            myL2NormSet = 0.0_Kr
            PetscCallA(DMGetStratumIS(dmU,MEF90SetLabelName(setType),setID(set),setPointIS,ierr))
            PetscCallA(ISGetIndicesF90(setPointIS,setPointID,ierr))
            PetscCallA(DMPlexGetCellType(dmU,setPointID(1),cellType,ierr))
            PetscCallA(MEF90ElementGetType(MEF90GlobalOptions_default%elementFamily,MEF90GlobalOptions_default%elementOrder,cellType,elementType,ierr))
            quadratureOrder = elementType%order * 2
            If (dim == 2) Then
                PetscCallA(MEF90ElementCreate(dmU,setPointIS,elem2D,QuadratureOrder,elementType,ierr))
                PetscCallA(MEF90_MassMatrixAssembleSet(M,dmU,setType,setID(set),elem2D,elementType,ierr))
                !PetscCallA(MEF90L2DotProductSet(myL2NormSet,U,U,setType,setID(set),elem2D,elementType,ierr))
                PetscCallA(MEF90L2NormSet(myL2NormSet,U,setType,setID(set),elem2D,elementType,ierr))
                PetscCallA(MEF90ElementDestroy(elem2D,ierr))
            Else
                PetscCallA(MEF90ElementCreate(dmU,setPointIS,elem3D,QuadratureOrder,elementType,ierr))
                PetscCallA(MEF90_MassMatrixAssembleSet(M,dmU,setType,setID(set),elem3D,elementType,ierr))
                !PetscCallA(MEF90L2DotProductSet(myL2NormSet,U,U,setType,setID(set),elem3D,elementType,ierr))
                PetscCallA(MEF90L2NormSet(myL2NormSet,U,setType,setID(set),elem3D,elementType,ierr))
                PetscCallA(MEF90ElementDestroy(elem3D,ierr))
            End If
            PetscCallA(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
            PetscCallA(ISDestroy(setPointIS,ierr))
            Call MPI_AllReduce(myL2NormSet,L2NormSet,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD,ierr)
            L2Norm = L2Norm + L2NormSet
            Write(IOBuffer,'("set ",I4," my L2 norm ", ES12.5,"\n")') setID(set), myL2NormSet
            PetscCallA(PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
            PetscCallA(PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT,ierr))
            Write(IOBuffer,'("set ",I4," L2 norm ", ES12.5,"\n")') setID(set), L2NormSet
            PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
        End Do
        PetscCallA(ISRestoreIndicesF90(setIS,setID,ierr))
    End If ! setIS
    PetscCallA(ISDestroy(setIS,ierr))
    PetscCallA(MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY,ierr))
    PetscCallA(MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY,ierr))
    Write(IOBuffer,'(" L2 norm ", ES12.5,"\n")') L2NormSet
    PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
    PetscCallA(MatViewFromOptions(M,PETSC_NULL_OPTIONS,"-mef90matM_view",ierr))

    setType = MEF90FaceSetType
    PetscCallA(DMGetLabelIdIS(dmU0,MEF90SetLabelName(setType),setIS,ierr))
    If (setIS /= PETSC_NULL_IS) Then
        PetscCallA(ISGetIndicesF90(setIS,setID,ierr))
        QuadratureOrder = elementType%order * 2
        Do set = 1,size(setID)
            PetscCallA(DMGetStratumIS(dmU0,MEF90SetLabelName(setType),setID(set),setPointIS,ierr))
            PetscCallA(ISGetIndicesF90(setPointIS,setPointID,ierr))
            PetscCall(DMPlexGetCellType(dmU0,setPointID(1),faceType,ierr))
            PetscCall(MEF90ElementGetTypeBoundary(MEF90GlobalOptions_default%elementFamily,MEF90GlobalOptions_default%elementOrder,faceType,elementType,ierr))
                If (dim == 2) Then
                PetscCallA(MEF90ElementCreate(dmU0,setPointIS,elem2D,QuadratureOrder,elementType,ierr))
                PetscCallA(MEF90_MassMatrixAssembleSet(M0,dmU0,setType,setID(set),elem2D,elementType,ierr))
                PetscCallA(MEF90ElementDestroy(elem2D,ierr))
            Else
                PetscCallA(MEF90ElementCreate(dmU0,setPointIS,elem3D,QuadratureOrder,elementType,ierr))
                PetscCallA(MEF90_MassMatrixAssembleSet(M,dmU0,setType,setID(set),elem3D,elementType,ierr))
                PetscCallA(MEF90ElementDestroy(elem3D,ierr))
            End If
            PetscCallA(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
            PetscCallA(ISDestroy(setPointIS,ierr))
        End Do
        PetscCallA(ISRestoreIndicesF90(setIS,setID,ierr))
    End If ! setIS
    PetscCallA(ISDestroy(setIS,ierr))
    PetscCallA(MatAssemblyBegin(M0,MAT_FINAL_ASSEMBLY,ierr))
    PetscCallA(MatAssemblyEnd(M0,MAT_FINAL_ASSEMBLY,ierr))
    PetscCallA(MatViewFromOptions(M0,PETSC_NULL_OPTIONS,"-mef90matM0_view",ierr))


    TestMass: Block
        Type(tVec)                 :: Uglob,MUglob
        PetscReal                  :: Mass

        PetscCallA(dmCreateGlobalVector(dmU,Uglob,ierr))
        PetscCallA(VecDuplicate(Uglob,MUglob,ierr))
        PetscCallA(VecSet(UGlob,1.0_Kr,ierr))
        PetscCallA(MatMult(M,Uglob,MUglob,ierr))
        PetscCallA(VecDot(UGlob,MUGlob,mass,ierr))
        Write(IOBuffer,'("Total mass", ES12.5," \n")') mass
        PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer, ierr))
        PetscCallA(VecViewFromOptions(Uglob,PETSC_NULL_OPTIONS,"-mef90VecU_view",ierr))
        PetscCallA(VecViewFromOptions(MUglob,PETSC_NULL_OPTIONS,"-mef90VecMU_view",ierr))
        PetscCallA(VecDestroy(Uglob,ierr))
        PetscCallA(VecDestroy(MUglob,ierr))

        PetscCallA(dmCreateGlobalVector(dmU0,Uglob,ierr))
        PetscCallA(VecDuplicate(Uglob,MUglob,ierr))
        PetscCallA(VecSet(UGlob,1.0_Kr,ierr))
        PetscCallA(MatMult(M0,Uglob,MUglob,ierr))
        PetscCallA(VecDot(UGlob,MUGlob,mass,ierr))
        Write(IOBuffer,'("Boundary mass", ES12.5," \n")') mass
        PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer, ierr))
        PetscCallA(VecDestroy(Uglob,ierr))
        PetscCallA(VecDestroy(MUglob,ierr))
    End Block TestMass

    PetscCallA(MatDestroy(M,ierr))
    PetscCallA(MatDestroy(M0,ierr))
    PetscCallA(VecDestroy(U,ierr))
    PetscCallA(VecDestroy(U0,ierr))
    PetscCallA(DMDestroy(dm,ierr))

    Call MEF90CtxDestroy(MEF90Ctx,ierr)   
    Call MEF90Finalize(ierr)
    PetscCallA(PetscFinalize(ierr))
End Program TestMassMatrix
           
    