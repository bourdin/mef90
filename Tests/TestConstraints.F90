module localFunctions
#include <petsc/finclude/petsc.h>
use petsc
use m_MEF90
implicit none

    abstract interface
        pure function f_interface(x,y)
            PetscReal,intent(in) :: x,y
            PetscReal            :: f_interface
        end function f_interface
    end interface
        
contains
#undef __FUNCT__
#define __FUNCT__ "f1"
    pure function f1(x,y)
        PetscReal,intent(in)   :: x,y
        PetscReal              :: f1

        f1 = x**2 + y**2
    end function f1

#undef __FUNCT__
#define __FUNCT__ "f2"
    pure function f2(x,y)
        PetscReal,intent(in)   :: x,y
        PetscReal              :: f2

        f2 = 100.0_Kr * x**2 + y**2
    end function f2

#undef __FUNCT__
#define __FUNCT__
    subroutine project(v,s,f,ierr)
        Type(tVec),intent(IN)              :: v
        Type(tPetscSection),intent(IN)     :: s
        procedure(f_interface)             :: f
        PetscErrorCode,intent(INOUT)       :: ierr

        PetscInt                           :: pStart,pEnd,p,numDof,i,pOffset
        Type(tDM)                          :: dm
        Type(tPetscSection)                :: coordSection
        Type(tVec)                         :: coordVec,vLoc
        PetscScalar,dimension(:),Pointer   :: coordArray,vArray
        PetscScalar,dimension(3)           :: xyz
        PetscInt                           :: dim

        PetscCall(PetscSectionGetChart(s,pStart,pEnd,ierr))
        PetscCall(VecGetDM(v,dm,ierr))
        PetscCall(DMGetLocalVector(dm,vLoc,ierr))
        PetscCall(DMGetCoordinateSection(dm,coordSection,ierr))
        PetscCall(DMGetCoordinatesLocal(dm,coordVec,ierr))
        PetscCall(DmGetDimension(dm,dim,ierr))

        PetscCall(VecGetArrayF90(vLoc,vArray,ierr))
        Do p = pStart,pEnd-1
            PetscCall(PetscSectionGetDof(s,p,numDof,ierr))
            If (numDof > 0) Then
                !!! trick: the coordinate of a point is the average of the coordinates of the points in its closure
                PetscCall(PetscSectionGetOffset(s,p,pOffset,ierr))
                PetscCall(DMPlexVecGetClosure(dm,coordSection,coordVec,p,coordArray,ierr))
                xyz = 0.0_Kr
                Do i = 1,dim
                    xyz(i) = sum(coordArray(i:size(coordArray):dim)) 
                    vArray(pOffset+i) = f(xyz(1),xyz(2)) + 10**(i-1)
                End Do
                xyz = xyz * dim / size(coordArray)
                PetscCall(DMPlexVecrestoreClosure(dm,coordSection,coordVec,p,coordArray,ierr))
            End If
        End Do
        PetscCall(VecRestoreArrayF90(vLoc,vArray,ierr))
        PetscCall(DMLocalToGlobal(dm,vLoc,INSERT_VALUES,v,ierr))
        PetscCall(DMRestoreLocalVector(dm,vLoc,ierr))
    End subroutine project    
end module localFunctions

Program  TestConstraints
#include <petsc/finclude/petsc.h>
Use m_MEF90
Use petsc
Use localFunctions
Implicit NONE   
    
    PetscErrorCode                      :: ierr
    Type(MEF90Ctx_Type),target          :: MEF90Ctx
    Type(MEF90CtxGlobalOptions_Type)    :: MEF90GlobalOptions_default
    Type(tDM),target                    :: dm,dmDist,dmU,dmU0
    PetscBool                           :: interpolate = PETSC_TRUE
    Character(len=MEF90_MXSTRLEN)       :: IOBuffer,setType

    PetscInt                            :: numComponents,numFConstraints = 1,numVConstraints = 3
    PetscInt                            :: set
    type(tIS)                           :: CSIS,FSIS,VSIS,setIS
    PetscInt,Dimension(:),pointer       :: CSID,FSID,VSID,setID
    PetscInt                            :: dim,pStart,pEnd,depth,numdof,i,p
    Type(tPetscSection)                 :: sectionU,sectionU0
    Logical,Dimension(:,:),Pointer      :: ConstraintTruthTableU,ConstraintTruthTableU0
    Logical,Dimension(:),Pointer        :: Constraints
    Type(tVec)                          :: U,U0

    PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
    Call MEF90Initialize(ierr)

    MEF90GlobalOptions_default%verbose           = 1
    MEF90GlobalOptions_default%dryrun            = PETSC_FALSE
    MEF90GlobalOptions_default%timeInterpolation = MEF90TimeInterpolation_linear
    MEF90GlobalOptions_default%timeMin           = 0.0_Kr
    MEF90GlobalOptions_default%timeMax           = 1.0_Kr
    MEF90GlobalOptions_default%timeNumStep       = 11
    MEF90GlobalOptions_default%timeSkip          = 0
    MEF90GlobalOptions_default%timeNumCycle      = 1
    MEF90GlobalOptions_default%fileFormat        = MEF90FileFormat_EXOSingle

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

    PetscCallA(DMGetDimension(dm,dim,ierr))
    PetscCallA(DMPlexGetChart(dm,pStart,pEnd,ierr))

    PetscCallA(PetscSectionCreate(MEF90Ctx%Comm,sectionU,ierr))
    PetscCall(PetscObjectSetName(SectionU,"Section for U",ierr))
    PetscCallA(PetscSectionSetChart(sectionU,pStart,pEnd,ierr))
    PetscCallA(PetscSectionCreate(MEF90Ctx%Comm,sectionU0,ierr))
    PetscCall(PetscObjectSetName(SectionU,"Section for boundary values of U",ierr))
    PetscCallA(PetscSectionSetChart(sectionU0,pStart,pEnd,ierr))


    numComponents = dim
    !!! Allocate DoF at cell and face sets
    !!! Note that if the face sets corresponds to faces in elements in cell set 
    !!! (which will always be the case in an exodusII mesh), the second call does nothing
    If (dim == 2) Then
        PetscCallA(MEF90_SectionAllocateDof(dm,"Cell Sets",MEF90_P1_Lagrange_2D,numComponents,sectionU,ierr))
        PetscCallA(MEF90_SectionAllocateDof(dm,"Face Sets",MEF90_P1_Lagrange_2DBoundary,numComponents,sectionU,ierr))
        PetscCallA(MEF90_SectionAllocateDof(dm,"Face Sets",MEF90_P1_Lagrange_2DBoundary,numComponents,sectionU0,ierr))
    Else If (dim == 3) Then
        PetscCallA(MEF90_SectionAllocateDof(dm,"Cell Sets",MEF90_P1_Lagrange_3D,numComponents,sectionU,ierr))
        PetscCallA(MEF90_SectionAllocateDof(dm,"Face Sets",MEF90_P1_Lagrange_3DBoundary,numComponents,sectionU,ierr))
        PetscCallA(MEF90_SectionAllocateDof(dm,"Face Sets",MEF90_P1_Lagrange_3DBoundary,numComponents,sectionU0,ierr))
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
    PetscCallA(PetscSectionViewFromOptions(SectionU,PETSC_NULL_OPTIONS,"-section_view",ierr))
    PetscCallA(PetscSectionViewFromOptions(SectionU0,PETSC_NULL_OPTIONS,"-section_view",ierr))

    PetscCall(DMClone(dm,dmU,ierr))
    PetscCallA(DMSetLocalSection(dmU,sectionU,ierr))
    !PetscCallA(DMCreateLocalVector(dmU,Uloc,ierr))
    PetscCallA(DMCreateGlobalVector(dmU,U,ierr))
    PetscCall(PetscObjectSetName(U,"U vector",ierr))
    PetscCallA(VecSet(U,0.0_Kr,ierr))

    PetscCall(DMClone(dm,dmU0,ierr))
    PetscCallA(DMSetLocalSection(dmU0,sectionU0,ierr))
    PetscCallA(DMCreateGlobalVector(dmU0,U0,ierr))
    PetscCall(PetscObjectSetName(U0,"U boundary vector U0",ierr))
    PetscCallA(VecSet(U0,0.0_Kr,ierr))

    PetscCall(project(U,sectionU,f1,ierr))
    PetscCall(project(U0,sectionU0,f1,ierr))
    PetscCallA(VecViewFromOptions(U,PETSC_NULL_OPTIONS,"-vec_view",ierr))
    PetscCallA(VecViewFromOptions(U0,PETSC_NULL_OPTIONS,"-vec_view",ierr))

    PetscCallA(VecDestroy(U,ierr))
    PetscCallA(VecDestroy(U0,ierr))
    PetscCallA(PetscSectionDestroy(sectionU,ierr))
    PetscCallA(PetscSectionDestroy(sectionU0,ierr))
    PetscCallA(DMDestroy(dm,ierr))
    PetscCallA(DMDestroy(dmU,ierr))
    PetscCallA(DMDestroy(dmU0,ierr))

    Call MEF90CtxDestroy(MEF90Ctx,ierr)   
    Call MEF90Finalize(ierr)
    Call PetscFinalize(ierr)
End Program  TestConstraints
 
       
