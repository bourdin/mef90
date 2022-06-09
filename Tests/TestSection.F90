Program  TestSection
#include <petsc/finclude/petsc.h>
Use m_MEF90
Use petsc
Implicit NONE   
    
    PetscErrorCode                      :: ierr
    Type(MEF90Ctx_Type),target          :: MEF90Ctx
    Type(MEF90CtxGlobalOptions_Type)    :: MEF90GlobalOptions_default
    Type(tDM),target                    :: dm,dmDist
    PetscBool                           :: interpolate = PETSC_TRUE
    Character(len=MEF90_MXSTRLEN)       :: IOBuffer

    PetscInt                            :: numComponents = 3, numFConstraints = 1, numVConstraints = 3
    PetscInt                            :: set
    type(tIS)                           :: CSIS,VSIS
    PetscInt,Dimension(:),pointer       :: setID
    PetscInt                            :: dim,pStart,pEnd,depth,numdof,i
    Type(tPetscSection)                 :: section
    Logical,Dimension(:,:),Allocatable  :: ConstraintTruthTable
    PetscInt,Dimension(:),Pointer       :: ConstraintIndices
    Type(tVec)                          :: v

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
    Else
        dmDist = dm
    end if
    PetscCall(DMViewFromOptions(dmDist,PETSC_NULL_OPTIONS,"-dm_view",ierr))

    PetscCall(PetscSectionCreate(MEF90Ctx%Comm,section,ierr))
    PetscCall(DMPlexGetChart(dmDist, pStart, pEnd,ierr))
    PetscCall(PetscSectionSetChart(section, pStart, pEnd,ierr))

    !!! Allocate DoF at cell sets
    PetscCall(DMGetLabelIdIS(dmDist,MEF90_DMPlexCellSetType, CSIS, ierr))
    PetscCall(ISGetIndicesF90(CSIS,setID,ierr))
    ! Write(IOBuffer,*) 'Sets for MEF90_SectionAllocateDofSet         ', setID, '\n'
    ! PetscCall(PetscPrintf(Mef90Ctx%comm,IOBuffer, ierr))
    Do set = 1, size(setID)
        PetscCall(MEF90_SectionAllocateDofSet(dmDist,MEF90_DMPlexCellSetType,setID(set),MEF90_P1_Lagrange_2D,numComponents,Section, ierr))
    End Do
    PetscCall(ISRestoreIndicesF90(CSIS,setID,ierr))
    PetscCall(ISDestroy(CSIS,ierr))

    !PetscCall(PetscSectionSetup(Section,ierr))

    ! !!! Allocate constraints at face sets
    ! PetscCall(DMGetLabelIdIS(dmDist,MEF90_DMPlexFaceSetType, CSIS, ierr))
    ! PetscCall(ISGetIndicesF90(CSIS,setID,ierr))
    ! ! Write(IOBuffer,*) 'Sets for MEF90_SectionAllocateConstraintSet  ', setID, '\n'
    ! ! PetscCall(PetscPrintf(Mef90Ctx%comm,IOBuffer, ierr))
    ! Do set = 1, size(setID)
    !     PetscCall(MEF90_SectionAllocateConstraintSet(dmDist,MEF90_DMPlexFaceSetType,setID(set),numFConstraints,section,ierr))
    ! End Do
    ! PetscCall(ISRestoreIndicesF90(CSIS,setID,ierr))
    ! PetscCall(ISDestroy(CSIS,ierr))

    !! Allocate constraints at vertex sets
    PetscCall(DMGetLabelIdIS(dmDist,MEF90_DMPlexVertexSetType, VSIS, ierr))
    PetscCall(ISGetIndicesF90(VSIS,setID,ierr))
    Write(IOBuffer,*) 'Sets for MEF90_SectionAllocateConstraintSet  ', setID, '\n'
    PetscCall(PetscPrintf(Mef90Ctx%comm,IOBuffer, ierr))
    Do set = 1, size(setID)
        PetscCall(MEF90_SectionAllocateConstraintSet(dmDist,MEF90_DMPlexVertexSetType,setID(set),numVConstraints,section,ierr))
    End Do
    PetscCall(ISRestoreIndicesF90(VSIS,setID,ierr))
    PetscCall(ISDestroy(VSIS,ierr))

    PetscCall(PetscSectionSetup(Section,ierr))
!     !!! Set constraints at Face sets
!     Allocate(constraintIndices(numFConstraints))
!     constraintIndices = [ (i-1, i = 1,numFConstraints)]

!     PetscCall(DMGetLabelIdIS(dmDist,MEF90_DMPlexFaceSetType, CSIS, ierr))
!     PetscCall(ISGetIndicesF90(CSIS,setID,ierr))
!     Write(IOBuffer,*) 'Sets for MEF90_SectionSetConstraintIndicesSet', setID, '\n'
!     PetscCall(PetscPrintf(Mef90Ctx%comm,IOBuffer, ierr))
!     Do set = 1, size(setID)
!         PetscCall(MEF90_SectionSetConstraintIndicesSet(dmDist,MEF90_DMPlexFaceSetType,setID(set),constraintIndices,section,ierr))
!     End Do
!     PetscCall(ISRestoreIndicesF90(CSIS,setID,ierr))
!     PetscCall(ISDestroy(CSIS,ierr))
!     DeAllocate(constraintIndices)

    PetscCall(PetscSectionGetStorageSize(section,numdof,ierr))
    Write(*,*) '---', numdof / numComponents
    !!! Set constraints at Vertex sets
    numVConstraints = numVConstraints
    Allocate(constraintIndices(numVConstraints))
    constraintIndices = [ (i-1, i = 1,numVConstraints)]

    PetscCall(DMGetLabelIdIS(dmDist,MEF90_DMPlexVertexSetType, CSIS, ierr))
    PetscCall(ISGetIndicesF90(CSIS,setID,ierr))
    Write(IOBuffer,*) 'Sets for MEF90_SectionSetConstraintIndicesSet', setID, '\n'
    PetscCall(PetscPrintf(Mef90Ctx%comm,IOBuffer, ierr))
    Do set = 1, size(setID)
        PetscCall(MEF90_SectionSetConstraintIndicesSet(dmDist,MEF90_DMPlexVertexSetType,setID(set),constraintIndices,section,ierr))
    End Do
    PetscCall(ISRestoreIndicesF90(CSIS,setID,ierr))
    PetscCall(ISDestroy(CSIS,ierr))
    DeAllocate(constraintIndices)

    PetscCall(PetscSectionSetup(Section,ierr))

    PetscCall(PetscSectionViewFromOptions(Section,PETSC_NULL_OPTIONS,"-section_view",ierr))

    PetscCall(DMSetLocalSection(dmDist, section,ierr))
    PetscCall(DMGetGlobalVector(dmDist,v,ierr))
    PetscCall(VecSet(v,0.0_Kr,ierr))
    PetscCall(VecViewFromOptions(v,PETSC_NULL_OPTIONS,"-vec_view",ierr))


    PetscCall(DMPlexGetChart(dmDist,pStart,pEnd,ierr))
    write(*,*)'=== ', pstart,pend
    Allocate(ConstraintTruthTable(pEnd,numComponents))

    PetscCall(DMGetDimension(dmDist,dim,ierr))
    Do depth = 1, dim+1
        PetscCall(DMPlexGetDepthStratum(dmDist,depth-1,pStart,pEnd,ierr))
        write(*,*) MEF90Ctx%rank, '+++', depth, pStart,pEnd
    End Do



    PetscCall(PetscSectionDestroy(section,ierr))
    PetscCall(DMDestroy(dmDist,ierr))

    Call MEF90CtxDestroy(MEF90Ctx,ierr)   
    Call MEF90Finalize(ierr)
    Call PetscFinalize(ierr)
End Program  TestSection
 
       