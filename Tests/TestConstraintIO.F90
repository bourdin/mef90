Program  TestConstraintIO
#include <petsc/finclude/petsc.h>
Use m_MEF90
Use petsc
Implicit NONE   
    
    PetscErrorCode                      :: ierr
    Type(MEF90Ctx_Type),target          :: MEF90Ctx
    Type(MEF90CtxGlobalOptions_Type)    :: MEF90GlobalOptions_default
    Type(tDM),target                    :: dm,dmDist,dmU,dmU0
    PetscBool                           :: interpolate = PETSC_TRUE
    Character(len=MEF90_MXSTRLEN)       :: IOBuffer,setType

    PetscInt                            :: numComponents,numFConstraints = 1,numVConstraints = 3
    PetscInt                            :: set, order = 1
    PetscBool                           :: flg
    type(tIS)                           :: CSIS,FSIS,VSIS,setIS
    PetscInt,Dimension(:),pointer       :: CSID,FSID,VSID,setID
    PetscInt                            :: dim,pStart,pEnd,depth,numdof,i,globalSize,localSize,size
    Type(tPetscSection)                 :: sectionU, sectionU0, gSectionU0
    Type(tPetscSF)                      :: naturalPointSF, natSFU, natSFU0, lioSFU, iolSFU, lioSFU0, iolSFU0, lcSF, clSF
    Logical,Dimension(:,:),Pointer      :: ConstraintTruthTableU, ConstraintTruthTableU0
    Logical,Dimension(:),Pointer        :: constraints
    Type(tVec)                          :: ioU, locCoord, locU, gU, natU, ioU0, locU0, gU0, natU0
    PetscInt                            :: numstep = 1, step
    PetscInt                            :: numNodalVar, numZonalVar = 0, numCS

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
    PetscCallA(DMGetDimension(dm, dim,ierr))
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-dm_view",ierr))
    
    PetscCallA(DMSetUseNatural(dm, PETSC_TRUE, ierr));
    PetscCallA(DMPlexDistribute(dm,0,naturalPointSF,dmDist,ierr))
    if (MEF90Ctx%NumProcs > 1) then
        PetscCallA(DMDestroy(dm,ierr))
        dm = dmDist
        PetscCallA(DMPlexSetMigrationSF(dm,naturalPointSF,ierr))
        PetscCallA(PetscSFDestroy(naturalPointSF,ierr))
    end if
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-dm_view",ierr))

    PetscCallA(PetscSectionCreate(MEF90Ctx%Comm,sectionU,ierr))
    PetscCallA(DMPlexGetChart(dm,pStart,pEnd,ierr))
    PetscCallA(PetscSectionSetChart(sectionU,pStart,pEnd,ierr))
    PetscCallA(PetscSectionCreate(MEF90Ctx%Comm,sectionU0,ierr))
    PetscCall(PetscObjectSetName(SectionU0,"Section for boundary values of U",ierr))
    PetscCallA(PetscSectionSetChart(sectionU0,pStart,pEnd,ierr))
    PetscCallA(DMGetDimension(dm,dim,ierr))

    numComponents = dim
    Call PetscOptionsGetInt(PETSC_NULL_OPTIONS,'','-order',order,flg,ierr);CHKERRQ(ierr);
    If (dim == 2) Then
        Select case(order)
        case(1)
            PetscCallA(MEF90_SectionAllocateDof(dm,"Cell Sets",MEF90_P1_Lagrange_2D,numComponents,sectionU,ierr))
            PetscCallA(MEF90_SectionAllocateDof(dm,"Face Sets",MEF90_P1_Lagrange_2DBoundary,numComponents,sectionU,ierr))
            PetscCallA(MEF90_SectionAllocateDof(dm,"Face Sets",MEF90_P1_Lagrange_2DBoundary,numComponents,sectionU0,ierr))
        case(2)
            PetscCallA(MEF90_SectionAllocateDof(dm,"Cell Sets",MEF90_P2_Lagrange_2D,numComponents,sectionU,ierr))
            PetscCallA(MEF90_SectionAllocateDof(dm,"Face Sets",MEF90_P2_Lagrange_2DBoundary,numComponents,sectionU,ierr))
            PetscCallA(MEF90_SectionAllocateDof(dm,"Face Sets",MEF90_P2_Lagrange_2DBoundary,numComponents,sectionU0,ierr))
        Case default
            Write(IOBuffer,*) 'ERROR: unimplemented order ', order, '\n'
            SETERRA(MEF90Ctx%Comm,PETSC_ERR_USER,IOBuffer)
        End Select
    Else If (dim == 3) Then
        Select case(order)
            case(1)
            PetscCallA(MEF90_SectionAllocateDof(dm,"Cell Sets",MEF90_P1_Lagrange_3D,numComponents,sectionU,ierr))
            PetscCallA(MEF90_SectionAllocateDof(dm,"Face Sets",MEF90_P1_Lagrange_3DBoundary,numComponents,sectionU,ierr))
            PetscCallA(MEF90_SectionAllocateDof(dm,"Face Sets",MEF90_P1_Lagrange_3DBoundary,numComponents,sectionU0,ierr))
        case(2)
            PetscCallA(MEF90_SectionAllocateDof(dm,"Cell Sets",MEF90_P2_Lagrange_3D,numComponents,sectionU,ierr))
            PetscCallA(MEF90_SectionAllocateDof(dm,"Face Sets",MEF90_P2_Lagrange_3DBoundary,numComponents,sectionU,ierr))
            PetscCallA(MEF90_SectionAllocateDof(dm,"Face Sets",MEF90_P2_Lagrange_3DBoundary,numComponents,sectionU0,ierr))
        Case default
            Write(IOBuffer,*) 'ERROR: unimplemented order ', order, '\n'
            SETERRA(MEF90Ctx%Comm,PETSC_ERR_USER,IOBuffer)
        End Select
    End If


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
    Allocate(ConstraintTruthTableU0(pEnd,numComponents))
    ConstraintTruthTableU  = .FALSE.
    ConstraintTruthTableU0 = .FALSE.

    Allocate(constraints(numComponents))
    constraints = .FALSE.

    setType = MEF90_DMPlexFaceSetType
    PetscCallA(DMGetLabelIdIS(dm,setType,SetIS,ierr))
    PetscCallA(ISGetIndicesF90(SetIS,setID,ierr))
    Do set = 1,size(setID)
        !!! setting the constrained components to an arbitrary value
        !!! In real life, we would get constraint from the CS/FS/ES/VS bag
        constraints = .FALSE.
        constraints(mod(setID(set),numComponents)+1) = .TRUE.
        PetscCallA(MEF90_SetupConstraintTableSet(dm,sectionU,setType,setID(set),constraints,ConstraintTruthTableU,ierr))
        PetscCallA(MEF90_SetupConstraintTableSet(dm,sectionU,setType,setID(set),constraints,ConstraintTruthTableU0,ierr))
    End Do
    PetscCallA(ISRestoreIndicesF90(SetIS,setID,ierr))
    PetscCallA(ISDestroy(SetIS,ierr))

    setType = MEF90_DMPlexVertexSetType
    PetscCallA(DMGetLabelIdIS(dm,setType,SetIS,ierr))
    PetscCallA(ISGetIndicesF90(SetIS,setID,ierr))
    Do set = 1,size(setID)
        !!! setting the constrained components to an arbitrary value
        !!! In real life, we would get constraint from the CS/FS/ES/VS bag
        constraints = .FALSE.
        constraints(mod(setID(set),numComponents)+1) = .TRUE.
        PetscCallA(MEF90_SetupConstraintTableSet(dm,sectionU,setType,setID(set),constraints,ConstraintTruthTableU,ierr))
        PetscCallA(MEF90_SetupConstraintTableSet(dm,sectionU0,setType,setID(set),constraints,ConstraintTruthTableU0,ierr))
    End Do
    PetscCallA(ISRestoreIndicesF90(SetIS,setID,ierr))
    PetscCallA(ISDestroy(SetIS,ierr))
    DeAllocate(constraints)

    PetscCallA(MEF90_SectionAllocateConstraint(dm,ConstraintTruthTableU,sectionU,ierr))
    PetscCallA(MEF90_SectionAllocateConstraint(dm,ConstraintTruthTableU0,sectionU0,ierr))

    DeAllocate(ConstraintTruthTableU)
    DeAllocate(ConstraintTruthTableU0)
    PetscCallA(PetscSectionViewFromOptions(SectionU,PETSC_NULL_OPTIONS,"-sectionU_view",ierr))
    PetscCallA(PetscSectionViewFromOptions(SectionU0,PETSC_NULL_OPTIONS,"-sectionU0_view",ierr))

    PetscCall(DMClone(dm,dmU,ierr))
    PetscCallA(DMSetUseNatural(dmU, PETSC_TRUE, ierr));
    PetscCallA(DMSetLocalSection(dmU,sectionU,ierr))

    PetscCall(DMClone(dm,dmU0,ierr))
    PetscCallA(DMSetUseNatural(dmU0, PETSC_TRUE, ierr));
    PetscCallA(DMSetLocalSection(dmU0,sectionU0,ierr))

    ! Creating the GlobalToNatural SF
    if (MEF90Ctx%NumProcs > 1) then
        PetscCallA(DMPlexGetMigrationSF(dm, naturalPointSF, ierr))
        PetscCallA(DMPlexSetMigrationSF(dmU, naturalPointSF, ierr))
        PetscCallA(DMPlexSetMigrationSF(dmU0, naturalPointSF, ierr))
        PetscCallA(DMPlexCreateGlobalToNaturalSF(dmU, PETSC_NULL_SECTION, naturalPointSF, natSFU, ierr))
        PetscCallA(DMSetNaturalSF(dmU, natSFU, ierr))
        PetscCallA(PetscObjectDereference(natSFU, ierr))
        PetscCallA(DMPlexCreateGlobalToNaturalSF(dmU0, PETSC_NULL_SECTION, naturalPointSF, natSFU0, ierr))
        PetscCallA(DMSetNaturalSF(dmU0, natSFU0, ierr))
        PetscCallA(PetscObjectDereference(natSFU0, ierr))
    end if

    PetscCallA(MEF90_CreateLocalToConstraintSF(MEF90Ctx,dmU,dmU0,lcSF,clSF,ierr))

    PetscCallA(PetscSFDestroy(lcSF,ierr))
    PetscCallA(PetscSFDestroy(clSF,ierr))

    PetscCallA(PetscSectionDestroy(sectionU,ierr))
    PetscCallA(PetscSectionDestroy(sectionU0,ierr))

    PetscCallA(DMDestroy(dmU,ierr))
    PetscCallA(DMDestroy(dmU0,ierr))
    PetscCallA(DMDestroy(dm,ierr))
    Call MEF90CtxDestroy(MEF90Ctx,ierr)   
    Call MEF90Finalize(ierr)
    Call PetscFinalize(ierr)
End Program  TestConstraintIO