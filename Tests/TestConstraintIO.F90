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
    Character(len=MEF90MXSTRLEN)            :: IOBuffer
    PetscEnum                               :: setType

    Type(MEF90ElementType)                  :: cellSetElementType,faceSetElementType
    PetscInt                                :: numComponents, numNodalVar = 2, numCellVar = 3, numGVar = 0, numStep = 3
    Character(len=MEF90MXSTRLEN),allocatable:: nodalVarName(:), cellVarName(:), gVarName(:)
    PetscInt                                :: set, order = 1, nroots, nleaves
    PetscBool                               :: flg
    PetscMPIInt                             :: rank = 0
    type(tIS)                               :: setIS
    PetscInt,Dimension(:),pointer           :: setID,remoteOffsets,ilocal
    PetscInt                                :: dim,pStart,pEnd,size,n,p
    Type(tPetscSection)                     :: sectionU, sectionU0
    Type(tPetscSF)                          :: naturalPointSF, natSFU, natSFU0, lcSF, clSF, idSF, testSF, lioSF, iolSF, lioBSF, iolBSF
    Type(tVec)                              :: locCoord, locVec, gVec, locVecB, ioVec
    Type(PetscSFNode),dimension(:),Pointer  :: iremote
    Type(tPetscViewer)                      :: viewer
    Logical,Dimension(:,:),Pointer          :: ConstraintTruthTableU, ConstraintTruthTableU0
    Logical,Dimension(:),Pointer            :: constraints

    PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
    Call MEF90Initialize(ierr)

    MEF90GlobalOptions_default%verbose           = 0
    MEF90GlobalOptions_default%dryrun            = PETSC_FALSE
    MEF90GlobalOptions_default%timeMin           = 0.0_Kr
    MEF90GlobalOptions_default%timeMax           = 1.0_Kr
    MEF90GlobalOptions_default%timeNumStep       = 11
    MEF90GlobalOptions_default%elementFamily     = MEF90ElementFamilyLagrange
    MEF90GlobalOptions_default%elementOrder      = 1
 
    Call MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90GlobalOptions_default,ierr)

    MEF90Ctx%resultFile = "results/TestConstraintIO_out.exo"
    Allocate(nodalVarName(numNodalVar))
    Allocate(cellVarName(numCellVar))
    Allocate(gVarName(numGVar))
    nodalVarName(1:2) = ["U_x ","U_y "]
    cellVarName(1:3) = ["Sigma_11","Sigma_22","Sigma_12"]
    
    PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm,MEF90Ctx%geometryfile,PETSC_NULL_CHARACTER,interpolate,dm,ierr))
    PetscCallA(DMPlexDistributeSetDefault(dm,PETSC_FALSE,ierr))
    PetscCallA(DMSetFromOptions(dm,ierr))
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-dm_view",ierr))

    call MEF90CtxOpenEXO(MEF90Ctx,viewer,ierr)
    call MEF90EXODMView(dm,viewer,ierr)
    call MEF90EXOFormat(viewer,gVarName,cellVarName,nodalVarName,numStep,ierr)

    DeAllocate(nodalVarName)
    DeAllocate(cellVarName)
    DeAllocate(gVarName)

    distribute: Block 
        Type(tDM),target                    :: dmDist
        If (MEF90Ctx%NumProcs > 1) Then
            PetscCallA(DMPlexDistribute(dm,0,naturalPointSF,dmDist,ierr))
            PetscCallA(DMPlexSetMigrationSF(dmDist,naturalPointSF, ierr))
            PetscCallA(PetscSFDestroy(naturalPointSF,ierr))
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
            cellSetElementType = MEF90P1Lagrange2D
            faceSetElementType = MEF90P1Lagrange2DBoundary
        case(2)
            cellSetElementType = MEF90P2Lagrange2D
            faceSetElementType = MEF90P2Lagrange2DBoundary
        Case default
            Write(IOBuffer,*) 'ERROR: unimplemented order ', order, '\n'
            SETERRA(MEF90Ctx%Comm,PETSC_ERR_USER,IOBuffer)
        End Select
    Else If (dim == 3) Then
        Select case(order)
            case(1)
                cellSetElementType = MEF90P1Lagrange3D
                faceSetElementType = MEF90P1Lagrange3DBoundary
            case(2)
                cellSetElementType = MEF90P2Lagrange3D
                faceSetElementType = MEF90P2Lagrange3DBoundary
            Case default
            Write(IOBuffer,*) 'ERROR: unimplemented order ', order, '\n'
            SETERRA(MEF90Ctx%Comm,PETSC_ERR_USER,IOBuffer)
        End Select
    End If
    PetscCallA(MEF90SectionAllocateDof(dm,MEF90CellSetType,cellSetElementType,numComponents,sectionU,ierr))
    PetscCallA(MEF90SectionAllocateDof(dm,MEF90FaceSetType,faceSetElementType,numComponents,sectionU,ierr))
    PetscCallA(MEF90SectionAllocateDof(dm,MEF90CellSetType,cellSetElementType,numComponents,sectionU0,ierr))
    PetscCallA(MEF90SectionAllocateDof(dm,MEF90FaceSetType,faceSetElementType,numComponents,sectionU0,ierr))


    !!! Allocate constraints.
    !!! This is definitely not optimized. We could create a table of dimension
    !!! # points with dof x # components instead of #points x # components
    !!! We can address this later if needed
    !!! The whole constraint setup takes 2 passes: 
    !!!   1. Fill the constraint truth table (typically from data in CS/FS/ES/VS bag)
    !!!      This is done in MEF90SetupConstraintTableSet
    !!!   2. Allocate space in the section, call PetscSectionSetup, and set the constraint indices
    !!!      for each constrained dof. This is done in MEF90SectionAllocateConstraint
    PetscCallA(DMPlexGetChart(dm,pStart,pEnd,ierr))
    Allocate(ConstraintTruthTableU(pEnd,numComponents),source=.FALSE.)
    Allocate(ConstraintTruthTableU0(pEnd,numComponents),source=.FALSE.)

    Allocate(constraints(numComponents))

    setType = MEF90FaceSetType
    PetscCallA(DMGetLabelIdIS(dm,MEF90SetLabelName(setType),SetIS,ierr))
    PetscCallA(ISGetIndicesF90(SetIS,setID,ierr))
    Do set = 1,size(setID)
        !!! setting the constrained components to an arbitrary value
        !!! In real life, we would get constraint from the CS/FS/ES/VS bag
        constraints = .FALSE.
        constraints(mod(setID(set),numComponents)+1) = .TRUE.
        PetscCallA(MEF90SetupConstraintTableSet(dm,sectionU,setType,setID(set),constraints,ConstraintTruthTableU,ierr))
        PetscCallA(MEF90SetupConstraintTableSet(dm,sectionU,setType,setID(set),constraints,ConstraintTruthTableU0,ierr))
    End Do
    PetscCallA(ISRestoreIndicesF90(SetIS,setID,ierr))
    PetscCallA(ISDestroy(SetIS,ierr))

    setType = MEF90VertexSetType
    PetscCallA(DMGetLabelIdIS(dm,MEF90SetLabelName(setType),SetIS,ierr))
    PetscCallA(ISGetIndicesF90(SetIS,setID,ierr))
    Do set = 1,size(setID)
        !!! setting the constrained components to an arbitrary value
        !!! In real life, we would get constraint from the CS/FS/ES/VS bag
        constraints = .FALSE.
        constraints(mod(setID(set),numComponents)+1) = .TRUE.
        PetscCallA(MEF90SetupConstraintTableSet(dm,sectionU,setType,setID(set),constraints,ConstraintTruthTableU,ierr))
        PetscCallA(MEF90SetupConstraintTableSet(dm,sectionU0,setType,setID(set),constraints,ConstraintTruthTableU0,ierr))
    End Do
    PetscCallA(ISRestoreIndicesF90(SetIS,setID,ierr))
    PetscCallA(ISDestroy(SetIS,ierr))
    DeAllocate(constraints)

    PetscCallA(MEF90SectionAllocateConstraint(dm,ConstraintTruthTableU,sectionU,ierr))
    PetscCallA(MEF90SectionAllocateConstraint(dm,ConstraintTruthTableU0,sectionU0,ierr))

    DeAllocate(ConstraintTruthTableU)
    DeAllocate(ConstraintTruthTableU0)
    PetscCallA(PetscSectionViewFromOptions(SectionU,PETSC_NULL_OPTIONS,"-section_view",ierr))
    PetscCallA(PetscSectionViewFromOptions(SectionU0,PETSC_NULL_OPTIONS,"-section_view",ierr))

    PetscCallA(DMSetUseNatural(dm,PETSC_TRUE,ierr))

    PetscCallA(DMClone(dm,dmU,ierr))
    PetscCallA(DMSetLocalSection(dmU,sectionU,ierr))
    PetscCallA(DMSetUseNatural(dmU,PETSC_TRUE,ierr))

    PetscCallA(DMClone(dm,dmU0,ierr))
    PetscCallA(DMSetLocalSection(dmU0,sectionU0,ierr))
    PetscCallA(DMSetUseNatural(dmU0,PETSC_TRUE,ierr))

    ! Creating the GlobalToNatural SF
    ! if (MEF90Ctx%NumProcs > 1) then
    !     PetscCallA(DMPlexGetMigrationSF(dm, naturalPointSF, ierr))
    !     PetscCallA(DMPlexSetMigrationSF(dmU, naturalPointSF, ierr))
    !     PetscCallA(DMPlexSetMigrationSF(dmU0, naturalPointSF, ierr))
    !     PetscCallA(DMPlexCreateGlobalToNaturalSF(dmU, PETSC_NULL_SECTION, naturalPointSF, natSFU, ierr))
    !     PetscCallA(DMSetNaturalSF(dmU, natSFU, ierr))
    !     PetscCallA(PetscObjectDereference(natSFU, ierr))
    !     PetscCallA(DMPlexCreateGlobalToNaturalSF(dmU0, PETSC_NULL_SECTION, naturalPointSF, natSFU0, ierr))
    !     PetscCallA(DMSetNaturalSF(dmU0, natSFU0, ierr))
    !     PetscCallA(PetscObjectDereference(natSFU0, ierr))
    ! end if

    PetscCallA(MEF90CreateLocalToConstraintSF(MEF90Ctx,dmU,dmU0,lcSF,clSF,ierr))
    ! PetscCallA(MEF90CreateLocalToIOSF(MEF90Ctx,dmU,lioSF,ierr))
    ! PetscCallA(MEF90CreateIOToLocalSF(MEF90Ctx,dmU,iolSF,ierr))
    ! PetscCallA(MEF90CreateLocalToIOSF(MEF90Ctx,dmU0,lioBSF,ierr))
    ! PetscCallA(MEF90CreateIOToLocalSF(MEF90Ctx,dmU0,iolBSF,ierr))

    ! PetscCallA(DMGetCoordinatesLocal(dmU,locCoord,ierr))
    ! PetscCallA(DMGetLocalVector(dmU,locVec,ierr))
    ! PetscCallA(DMGetLocalVector(dmU0,locVecB,ierr))
    PetscCallA(DMGetGlobalVector(dmU,gVec,ierr))
    PetscCallA(VecSet(gVec,-111.0_kr,ierr))
    ! PetscCallA(VecSet(locVecB,999.0_kr,ierr))
    ! PetscCallA(VecCopy(locCoord,locVec,ierr))

    ! PetscCallA(DMGlobalToLocalBegin(dmU,gVec,INSERT_VALUES,locVec,ierr))
    ! PetscCallA(DMGlobalToLocalEnd(dmU,gVec,INSERT_VALUES,locVec,ierr))

    ! PetscCallA(VecView(locVec,PETSC_VIEWER_STDOUT_SELF,ierr))

    ! PetscCallA(MEF90VecReorderingSF(locVecB,locVec,clSF,ierr))

    ! PetscCallMPIA(MPI_Barrier(PETSC_COMM_WORLD,ierr))

    ! PetscCallA(VecView(locVec,PETSC_VIEWER_STDOUT_SELF,ierr))

    ! PetscCallA(PetscSFGetGraph(lioSF,nroots,nleaves,ilocal,iremote,ierr))
    ! PetscCallA(VecCreateMPI(MEF90Ctx%Comm,nleaves,PETSC_DETERMINE,iovec,ierr))
    ! PetscCallA(PetscObjectSetName(iovec,"U",ierr))

    ! PetscCallA(MEF90VecReorderingSF(locVec,ioVec,lioSF,ierr))

    ! PetscCallA(VecView(ioVec,PETSC_VIEWER_STDOUT_WORLD,ierr))

    ! PetscCallMPIA(MPI_Barrier(PETSC_COMM_WORLD,ierr))

    PetscCallA(MEF90EXOVecView(gVec,viewer,ierr))

    ! PetscCallA(DMRestoreLocalVector(dmU,locVec,ierr))
    ! PetscCallA(DMRestoreLocalVector(dmU0,locVecB,ierr))
    PetscCallA(DMRestoreGlobalVector(dmU,gVec,ierr))
    ! PetscCallA(VecDestroy(ioVec,ierr))

    PetscCallA(PetscSFDestroy(lcSF,ierr))
    PetscCallA(PetscSFDestroy(clSF,ierr))
    ! PetscCallA(PetscSFDestroy(lioSF,ierr))
    ! PetscCallA(PetscSFDestroy(iolSF,ierr))
    ! PetscCallA(PetscSFDestroy(lioBSF,ierr))
    ! PetscCallA(PetscSFDestroy(iolBSF,ierr))

    PetscCallA(PetscSectionDestroy(sectionU,ierr))
    PetscCallA(PetscSectionDestroy(sectionU0,ierr))

    PetscCallA(DMDestroy(dmU,ierr))
    PetscCallA(DMDestroy(dmU0,ierr))
    PetscCallA(DMDestroy(dm,ierr))
    call MEF90CtxCloseEXO(viewer,ierr)
    Call MEF90CtxDestroy(MEF90Ctx,ierr)   
    Call MEF90Finalize(ierr)
    Call PetscFinalize(ierr)
End Program  TestConstraintIO