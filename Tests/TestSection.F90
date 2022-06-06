Program  TestDMPlexComputeGeometry
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

    PetscInt                            :: numComponents = 2
    PetscInt                            :: set!,cell
    type(tIS)                           :: CSIS!,FSIS,CellIS
    PetscInt,Dimension(:),pointer       :: setID!,cellID
    PetscInt                            :: dim,pStart,pEnd
    Type(tPetscSection)                 :: section
    PetscInt,Dimension(:),Pointer       :: ConstraintIndices


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
    if (MEF90Ctx%NumProcs == 1) then
        dmDist = dm
    end if
    PetscCall(DMViewFromOptions(dmDist,PETSC_NULL_OPTIONS,"-dm_view",ierr))

     PetscCall(PetscSectionCreate(MEF90Ctx%Comm,section,ierr))
     PetscCall(DMPlexGetChart(dm, pStart, pEnd,ierr))
     PetscCall(PetscSectionSetChart(section, pStart, pEnd,ierr))

    PetscCall(DMGetLabelIdIS(dm,MEF90_DMPlexCellSetType, CSIS, ierr))
    PetscCall(ISGetIndicesF90(CSIS,setID,ierr))
    Write(IOBuffer,*) 'Sets for MEF90_SectionAllocateDofSet       ', setID, '\n'
    PetscCall(PetscPrintf(Mef90Ctx%comm,IOBuffer, ierr))
    Do set = 1, size(setID)
        PetscCall(MEF90_SectionAllocateDofSet(dm,MEF90_DMPlexCellSetType,setID(set),MEF90_P1_Lagrange_2D,numComponents,Section, ierr))
    End Do
    PetscCall(ISRestoreIndicesF90(CSIS,setID,ierr))
    PetscCall(ISDestroy(CSIS,ierr))

    PetscCall(PetscSectionSetup(Section,ierr))


    Allocate(constraintIndices(numComponents))
    constraintIndices(1) = 0
    constraintIndices(2) = 1


    PetscCall(DMGetLabelIdIS(dm,MEF90_DMPlexFaceSetType, CSIS, ierr))
    PetscCall(ISGetIndicesF90(CSIS,setID,ierr))
    Write(IOBuffer,*) 'Sets for MEF90_SectionAllocateConstraintSet', setID, '\n'
    PetscCall(PetscPrintf(Mef90Ctx%comm,IOBuffer, ierr))
    Do set = 1, size(setID)
        PetscCall(MEF90_SectionAllocateConstraintSet(dm,MEF90_DMPlexFaceSetType,setID(set),numComponents,section,ierr))
    End Do
    PetscCall(ISRestoreIndicesF90(CSIS,setID,ierr))
    PetscCall(ISDestroy(CSIS,ierr))
    PetscCall(PetscSectionSetup(Section,ierr))

    PetscCall(DMGetLabelIdIS(dm,MEF90_DMPlexFaceSetType, CSIS, ierr))
    PetscCall(ISGetIndicesF90(CSIS,setID,ierr))
    Write(IOBuffer,*) 'Sets for MEF90_SectionSetConstraintIndicesSet', setID, '\n'
    PetscCall(PetscPrintf(Mef90Ctx%comm,IOBuffer, ierr))
    Do set = 1, size(setID)
        PetscCall(MEF90_SectionSetConstraintIndicesSet(dm,MEF90_DMPlexFaceSetType,setID(set),constraintIndices,section,ierr))
    End Do
    PetscCall(ISRestoreIndicesF90(CSIS,setID,ierr))
    PetscCall(ISDestroy(CSIS,ierr))
    PetscCall(PetscSectionSetup(Section,ierr))
    DeAllocate(constraintIndices)


    PetscCall(PetscSectionViewFromOptions(Section,PETSC_NULL_OPTIONS,"-section_view",ierr))
    PetscCall(PetscSectionDestroy(section,ierr))
    PetscCall(DMDestroy(dm,ierr))

    Call MEF90CtxDestroy(MEF90Ctx,ierr)   
    Call MEF90Finalize(ierr)
    Call PetscFinalize(ierr)
End Program  TestDMPlexComputeGeometry
 
       