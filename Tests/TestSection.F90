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
    Character(len=MEF90_MXSTRLEN)       :: IOBuffer,setType

    PetscInt                            :: numComponents,numFConstraints = 1,numVConstraints = 3
    PetscInt                            :: set
    type(tIS)                           :: CSIS,FSIS,VSIS,setIS
    PetscInt,Dimension(:),pointer       :: CSID,FSID,VSID,setID
    PetscInt                            :: dim,pStart,pEnd,depth,numdof,i
    Type(tPetscSection)                 :: section
    Logical,Dimension(:,:),Pointer      :: ConstraintTruthTable
    Logical,Dimension(:),Pointer        :: constraints
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

    PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
    
    Call MEF90Initialize(ierr)
    Call MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90GlobalOptions_default,ierr)
    
    PetscCallA(PetscPrintf(MEF90Ctx%Comm,MEF90Ctx%geometryfile//'\n',ierr))
    PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm,MEF90Ctx%geometryfile,PETSC_NULL_CHARACTER,interpolate,dm,ierr))
    PetscCallA(DMPlexDistributeSetDefault(dm,PETSC_FALSE,ierr))
    PetscCallA(DMSetFromOptions(dm,ierr))
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-dm_view",ierr))
    
    PetscCallA(DMPlexDistribute(dm,0,PETSC_NULL_SF,dmDist,ierr))
    if (MEF90Ctx%NumProcs > 1) then
        PetscCallA(DMDestroy(dm,ierr))
        dm = dmDist
    end if
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-dm_view",ierr))

    PetscCallA(PetscSectionCreate(MEF90Ctx%Comm,section,ierr))
    PetscCallA(DMPlexGetChart(dm,pStart,pEnd,ierr))
    PetscCallA(PetscSectionSetChart(section,pStart,pEnd,ierr))
    PetscCallA(DMGetDimension(dm,dim,ierr))

    numComponents = dim
    !!! Allocate DoF at cell and face sets
    !!! Note that if the face sets corresponds to faces in elements in cell set 
    !!! (which will always be the case in an exodusII mesh), the second call does nothing
    If (dim == 2) Then
        PetscCallA(MEF90_SectionAllocateDof(dm,"Cell Sets",MEF90_P1_Lagrange_2D,numComponents,section,ierr))
        PetscCallA(MEF90_SectionAllocateDof(dm,"Face Sets",MEF90_P1_Lagrange_2DBoundary,numComponents,section,ierr))
    Else If (dim == 3) Then
        PetscCallA(MEF90_SectionAllocateDof(dm,"Cell Sets",MEF90_P1_Lagrange_3D,numComponents,section,ierr))
        PetscCallA(MEF90_SectionAllocateDof(dm,"Face Sets",MEF90_P1_Lagrange_3DBoundary,numComponents,section,ierr))
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
    Allocate(ConstraintTruthTable(pEnd,numComponents))
    ConstraintTruthTable = .FALSE.

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
        PetscCallA(MEF90_SetupConstraintTableSet(dm,section,setType,setID(set),constraints,ConstraintTruthTable,ierr))
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
        PetscCallA(MEF90_SetupConstraintTableSet(dm,section,setType,setID(set),constraints,ConstraintTruthTable,ierr))
    End Do
    PetscCallA(ISRestoreIndicesF90(SetIS,setID,ierr))
    PetscCallA(ISDestroy(SetIS,ierr))
    DeAllocate(constraints)

    PetscCallA(MEF90_SectionAllocateConstraint(dm,ConstraintTruthTable,section,ierr))

    DeAllocate(ConstraintTruthTable)
    PetscCallA(PetscSectionViewFromOptions(Section,PETSC_NULL_OPTIONS,"-section_view",ierr))

    PetscCallA(DMSetLocalSection(dm,section,ierr))
    PetscCallA(DMCreateGlobalVector(dm,v,ierr))
    PetscCallA(VecSet(v,0.0_Kr,ierr))
    PetscCallA(VecViewFromOptions(v,PETSC_NULL_OPTIONS,"-vec_view",ierr))

    PetscCallA(VecDestroy(v,ierr))
    PetscCallA(PetscSectionDestroy(section,ierr))
    PetscCallA(DMDestroy(dm,ierr))
    Call MEF90CtxDestroy(MEF90Ctx,ierr)   
    Call MEF90Finalize(ierr)
    Call PetscFinalize(ierr)
End Program  TestSection
 
       
