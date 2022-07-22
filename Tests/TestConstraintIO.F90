Program  TestConstraintIO
#include <petsc/finclude/petsc.h>
Use m_MEF90
Use petsc
Implicit NONE   
    
    PetscErrorCode                          :: ierr
    Type(MEF90Ctx_Type),target              :: MEF90Ctx
    Type(MEF90CtxGlobalOptions_Type)        :: MEF90GlobalOptions_default
    Type(tDM),target                        :: dm,dmU,dmU0
    PetscBool                               :: interpolate = PETSC_TRUE
    Character(len=MEF90_MXSTRLEN)           :: IOBuffer
    PetscEnum                               :: setType

    Type(MEF90Element_Type)                 :: cellSetElementType,faceSetElementType
    PetscInt                                :: numComponents
    PetscInt                                :: set, order = 1
    PetscBool                               :: flg
    PetscMPIInt                             :: rank = 0
    type(tIS)                               :: setIS
    PetscInt,Dimension(:),pointer           :: setID,remoteOffsets
    PetscInt                                :: dim,pStart,pEnd,size,n,p
    Type(tPetscSection)                     :: sectionU, sectionU0
    Type(tPetscSF)                          :: naturalPointSF, natSFU, natSFU0, lcSF, clSF, idSF, testSF
    Type(PetscSFNode),dimension(:),Pointer  :: remote
    Logical,Dimension(:,:),Pointer          :: ConstraintTruthTableU, ConstraintTruthTableU0
    Logical,Dimension(:),Pointer            :: constraints

    PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
    Call MEF90Initialize(ierr)

    MEF90GlobalOptions_default%verbose           = 0
    MEF90GlobalOptions_default%dryrun            = PETSC_FALSE
    MEF90GlobalOptions_default%timeMin           = 0.0_Kr
    MEF90GlobalOptions_default%timeMax           = 1.0_Kr
    MEF90GlobalOptions_default%timeNumStep       = 11


    Call MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90GlobalOptions_default,ierr)
    
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

    PetscCallA(DMGetDimension(dm,dim,ierr))
    PetscCallA(DMPlexGetChart(dm,pStart,pEnd,ierr))

    PetscCallA(PetscSectionCreate(MEF90Ctx%Comm,sectionU,ierr))
    PetscCallA(PetscObjectSetName(SectionU,"Section for U",ierr))
    PetscCallA(PetscSectionSetChart(sectionU,pStart,pEnd,ierr))
    PetscCallA(PetscSectionCreate(MEF90Ctx%Comm,sectionU0,ierr))
    PetscCallA(PetscObjectSetName(SectionU,"Section for boundary values of U",ierr))
    PetscCallA(PetscSectionSetChart(sectionU0,pStart,pEnd,ierr))


    numComponents = dim
    !!! Allocate DoF at cell and face sets
    !!! Note that if the face sets corresponds to faces in elements in cell set 
    !!! (which will always be the case in an exodusII mesh), the second call does nothing
    PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,'','-order',order,flg,ierr))
    If (dim == 2) Then
        Select case(order)
        case(1)
            cellSetElementType = MEF90_P1_Lagrange_2D
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
    PetscCallA(MEF90_SectionAllocateDof(dm,MEF90_DMPlexcellSetType,cellSetElementType,numComponents,sectionU,ierr))
    PetscCallA(MEF90_SectionAllocateDof(dm,MEF90_DMPlexfaceSetType,faceSetElementType,numComponents,sectionU,ierr))
    PetscCallA(MEF90_SectionAllocateDof(dm,MEF90_DMPlexfaceSetType,faceSetElementType,numComponents,sectionU0,ierr))


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
    Allocate(ConstraintTruthTableU(pEnd,numComponents),source=.FALSE.)
    Allocate(ConstraintTruthTableU0(pEnd,numComponents),source=.FALSE.)

    Allocate(constraints(numComponents))

    setType = MEF90_DMPlexFaceSetType
    PetscCallA(DMGetLabelIdIS(dm,MEF90_DMPlexSetLabelName(setType),SetIS,ierr))
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
    PetscCallA(DMGetLabelIdIS(dm,MEF90_DMPlexSetLabelName(setType),SetIS,ierr))
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
    PetscCallA(PetscSectionViewFromOptions(SectionU,PETSC_NULL_OPTIONS,"-section_view",ierr))
    PetscCallA(PetscSectionViewFromOptions(SectionU0,PETSC_NULL_OPTIONS,"-section_view",ierr))

    PetscCallA(DMClone(dm,dmU,ierr))
    PetscCallA(DMSetLocalSection(dmU,sectionU,ierr))
    PetscCallA(DMSetUseNatural(dmU,PETSC_TRUE,ierr))

    PetscCallA(DMClone(dm,dmU0,ierr))
    PetscCallA(DMSetLocalSection(dmU0,sectionU0,ierr))
    PetscCallA(DMSetUseNatural(dmU0,PETSC_TRUE,ierr))

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

    PetscCallA(PetscSectionGetChart(sectionU, pStart, pEnd, ierr))
    n = pEnd-pStart
    Allocate(remote(n))
    Do p = 1,n
        remote(p)%rank = rank
        remote(p)%index = p-1
    End Do
    PetscCallA(PetscSFCreate(MEF90Ctx%Comm, idSF, ierr))
    PetscCallA(PetscSFSetFromOptions(idSF, ierr))
    PetscCallA(PetscSFSetGraph(idSF, n, n, PETSC_NULL_INTEGER, PETSC_COPY_VALUES, remote, PETSC_COPY_VALUES, ierr))
    PetscCallA(PetscSFSetUp(idSF, ierr))
    DeAllocate(remote)

    PetscCallA(PetscSFCreateRemoteOffsetsF90(idSF,sectionU,sectionU0,remoteOffsets,ierr))
    PetscCallA(PetscSFCreateSectionSFF90(idSF,sectionU,remoteOffsets,SectionU0,testSF,ierr))
!    PetscCallA(PetscSFDestroyRemoteOffsetsF90(remoteOffsets,ierr))
PetscCallA(PetscIntArray1dDestroyF90(remoteOffsets,ierr))
    PetscCallA(PetscSFView(testSF,PETSC_NULL_VIEWER,ierr))
    PetscCallA(PetscSFDestroy(testSF,ierr))
    PetscCallA(PetscSFDestroy(idSF,ierr))

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