Program  TestMassMatrix
#include <petsc/finclude/petsc.h>
    Use m_MEF90
    Use petsc
    Implicit NONE   
        
    PetscErrorCode                      :: ierr
    Type(MEF90Ctx_Type),target          :: MEF90Ctx
    Type(MEF90CtxGlobalOptions_Type)    :: MEF90GlobalOptions_default
    Type(tDM)                           :: dm
    PetscBool                           :: interpolate = PETSC_TRUE, hasConstraints = PETSC_FALSE
    Character(len=MEF90MXSTRLEN)       :: IOBuffer
    PetscEnum                           :: setType

    PetscInt                            :: numComponents
    PetscInt                            :: set
    Type(MEF90ElementType)             :: cellSetElementType,faceSetElementType
    type(tIS)                           :: setIS
    PetscInt,Dimension(:),pointer       :: setID
    PetscInt                            :: dim,pStart,pEnd,order = 1
    PetscBool                           :: flg
    Type(tPetscSection)                 :: section
    Logical,Dimension(:,:),Pointer      :: ConstraintTruthTableU
    Logical,Dimension(:),Pointer        :: Constraints
    Type(tMat)                          :: M

    PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
    Call MEF90Initialize(ierr)

    MEF90GlobalOptions_default%verbose           = 0
    MEF90GlobalOptions_default%dryrun            = PETSC_FALSE
    MEF90GlobalOptions_default%timeMin           = 0.0_Kr
    MEF90GlobalOptions_default%timeMax           = 1.0_Kr
    MEF90GlobalOptions_default%timeNumStep       = 11


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
        end If
    end Block distribute
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-mef90dm_view",ierr))

    PetscCallA(DMGetDimension(dm,dim,ierr))
    PetscCallA(DMPlexGetChart(dm,pStart,pEnd,ierr))

    PetscCallA(PetscSectionCreate(MEF90Ctx%Comm,section,ierr))
    PetscCallA(PetscObjectSetName(section,"Section for U",ierr))
    PetscCallA(PetscSectionSetChart(section,pStart,pEnd,ierr))


    numComponents = dim
    !!! Allocate DoF at cell and face sets
    !!! Note that if the face sets corresponds to faces in elements in cell set 
    !!! (which will always be the case in an exodusII mesh), the second call does nothing
    PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,'','-order',order,flg,ierr))
    If (dim == 2) Then
        Select case(order)
        case(1)
            cellSetElementType = MEF90P1Lagrange2D
            faceSetElementType = MEF90_P1_Lagrange_2DBoundary
        case(2)
            cellSetElementType = MEF90_P2_Lagrange_2D
            faceSetElementType = MEF90_P2_Lagrange_2DBoundary
        Case default
            Write(IOBuffer,*) 'ERROR: unimplemented order ', order, '\n'
            SETERRA(MEF90Ctx%Comm,PETSC_ERR_USER,IOBuffer)
        End Select
    Else If (dim == 3) Then
        Select case(order)
            case(1)
                cellSetElementType = MEF90_P1_Lagrange_3D
                faceSetElementType = MEF90_P1_Lagrange_3DBoundary
            case(2)
                cellSetElementType = MEF90_P2_Lagrange_3D
                faceSetElementType = MEF90_P2_Lagrange_3DBoundary
            Case default
            Write(IOBuffer,*) 'ERROR: unimplemented order ', order, '\n'
            SETERRA(MEF90Ctx%Comm,PETSC_ERR_USER,IOBuffer)
        End Select
    End If
    PetscCallA(MEF90_SectionAllocateDof(dm,MEF90_DMPlexcellSetType,cellSetElementType,numComponents,section,ierr))
    PetscCallA(MEF90_SectionAllocateDof(dm,MEF90_DMPlexfaceSetType,faceSetElementType,numComponents,section,ierr))


    PetscCallA(PetscOptionsGetBool(PETSC_NULL_OPTIONS,'','-constraint',hasConstraints,flg,ierr))
    If (hasConstraints) Then

        !!! Allocate constraints.
        !!! This is definitely not optimized. We could create a table of dimension
        !!! # points with dof x # components instead of #points x # components
        !!! We can address this later if needed
        !!! The whole constraint setup takes 2 passes: 
        !!!   1. Fill the constraint truth table (typically from data in CS/FS/ES/VS bag)
        !!!      This is done in MEF90_SetupConstraintTableSet
        !!!   2. Allocate space in the section, call PetscSectionSetup, and set the constraint indices
        !!!      for each constrained dof. This is done in MEF90_SectionAllocateConstraint
        PetscCallA(DMPlexGetChart(dm,pStart,pEnd,ierr))
        Allocate(ConstraintTruthTableU(pEnd,numComponents))
        ConstraintTruthTableU  = .FALSE.

        Allocate(constraints(numComponents))
        constraints = .FALSE.

        setType = MEF90_DMPlexfaceSetType
        PetscCallA(DMGetLabelIdIS(dm,MEF90_DMPlexSetLabelName(setType),setIS,ierr))
        PetscCallA(ISGetIndicesF90(setIS,setID,ierr))
        Do set = 1,size(setID)
            !!! setting the constrained components to an arbitrary value
            !!! In real life, we would get constraint from the CS/FS/ES/VS bag
            constraints = .FALSE.
            constraints(mod(setID(set),numComponents)+1) = .TRUE.
            PetscCallA(MEF90_SetupConstraintTableSet(dm,section,setType,setID(set),constraints,ConstraintTruthTableU,ierr))
        End Do
        PetscCallA(ISRestoreIndicesF90(setIS,setID,ierr))
        PetscCallA(ISDestroy(setIS,ierr))

        setType = MEF90_DMPlexVertexSetType
        PetscCallA(DMGetLabelIdIS(dm,MEF90_DMPlexSetLabelName(setType),setIS,ierr))
        PetscCallA(ISGetIndicesF90(setIS,setID,ierr))
        Do set = 1,size(setID)
            !!! setting the constrained components to an arbitrary value
            !!! In real life, we would get constraint from the CS/FS/ES/VS bag
            constraints = .FALSE.
            constraints(mod(setID(set),numComponents)+1) = .TRUE.
            PetscCallA(MEF90_SetupConstraintTableSet(dm,section,setType,setID(set),constraints,ConstraintTruthTableU,ierr))
        End Do
        PetscCallA(ISRestoreIndicesF90(setIS,setID,ierr))
        PetscCallA(ISDestroy(setIS,ierr))
        DeAllocate(constraints)

        PetscCallA(MEF90_SectionAllocateConstraint(dm,ConstraintTruthTableU,section,ierr))

        DeAllocate(ConstraintTruthTableU)
    Else
        PetscCallA(PetscSectionSetup(section,ierr))
    End If
    PetscCallA(PetscSectionViewFromOptions(section,PETSC_NULL_OPTIONS,"-mef90section_view",ierr))

    PetscCallA(DMSetLocalSection(dm,section,ierr))
    PetscCallA(DMCreateMatrix(dm,M,ierr))

    testClosure: Block
        PetscScalar,dimension(:),Pointer :: velem
        Type(tVec)                       :: v
        PetscInt                         :: point = 0

        PetscCallA(DMCreateGlobalVector(dm,v,ierr))
        PetscCallA(DMPlexVecGetClosure(dm,section,v,point,velem,ierr))
        write(*,*)'---', velem
        PetscCallA(DMPlexVecrestoreClosure(dm,section,v,point,velem,ierr))
        PetscCallA(VecDestroy(v,ierr))

    end block testClosure


    elementCreate: Block
        Type(MEF90Element2DVect),dimension(:),Pointer  :: elem2D
        Type(MEF90Element3DVect),dimension(:),Pointer  :: elem3D
        PetscInt                                       :: quadratureOrder = 2
        Type(tIS)                                      :: setPointIS

        setType = MEF90_DMPlexcellSetType
        PetscCallA(DMGetLabelIdIS(dm,MEF90_DMPlexSetLabelName(setType),setIS,ierr))
        If (setIS /= PETSC_NULL_IS) Then
            PetscCallA(ISGetIndicesF90(setIS,setID,ierr))
            Do set = 1,size(setID)
                PetscCallA(DMGetStratumIS(dm,MEF90_DMPlexSetLabelName(setType),setID(set),setPointIS,ierr))
                If (dim == 2) Then
                    PetscCallA(MEF90ElementCreate(dm,setPointIS,elem2D,QuadratureOrder,cellSetElementType,ierr))
                    PetscCallA(MEF90_MassMatrixAssembleSet(M,dm,setType,setID(set),elem2D,cellSetElementType,ierr))
                Else
                    PetscCallA(MEF90ElementCreate(dm,setPointIS,elem3D,QuadratureOrder,cellSetElementType,ierr))
                    PetscCallA(MEF90_MassMatrixAssembleSet(M,dm,setType,setID(set),elem3D,cellSetElementType,ierr))
                End If
                PetscCallA(ISDestroy(setPointIS,ierr))
            End Do
            PetscCallA(ISRestoreIndicesF90(setIS,setID,ierr))
        End If ! setIS
        PetscCallA(ISDestroy(setIS,ierr))
    End Block elementCreate
    PetscCallA(MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY,ierr))
    PetscCallA(MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY,ierr))

    PetscCallA(MatViewFromOptions(M,PETSC_NULL_OPTIONS,"-mef90mat_view",ierr))

    PetscCallA(MatDestroy(M,ierr))
    !PetscCallA(PetscSectionDestroy(section,ierr))
    PetscCallA(DMDestroy(dm,ierr))

    Call MEF90CtxDestroy(MEF90Ctx,ierr)   
    Call MEF90Finalize(ierr)
    Call PetscFinalize(ierr)
End Program TestMassMatrix
           
    