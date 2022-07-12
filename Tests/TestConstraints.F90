Module localFunctions
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
#define __FUNCT__ "project"
    subroutine project(v,s,f,ierr)
        Type(tVec),intent(IN)              :: v
        Type(tPetscSection),intent(IN)     :: s
        procedure(f_interface)             :: f
        PetscErrorCode,intent(INOUT)       :: ierr

        PetscInt                           :: pStart,pEnd,p,numDof,i
        Type(tDM)                          :: dm
        Type(tPetscSection)                :: coordSection
        Type(tVec)                         :: coordVec,vLoc
        PetscScalar,dimension(:),Pointer   :: coordArray,vArray
        PetscScalar,dimension(3)           :: xyz
        PetscInt                           :: dim

        PetscCallA(PetscSectionGetChart(s,pStart,pEnd,ierr))
        PetscCallA(VecGetDM(v,dm,ierr))
        PetscCallA(DMGetLocalVector(dm,vLoc,ierr))
        PetscCallA(DMGetCoordinateSection(dm,coordSection,ierr))
        PetscCallA(DMGetCoordinatesLocal(dm,coordVec,ierr))
        PetscCallA(DMGetDimension(dm,dim,ierr))

        Do p = pStart,pEnd-1
            PetscCallA(PetscSectionGetDof(s,p,numDof,ierr))
            If (numDof > 0) Then
                !!! trick: the coordinate of a point is the average of the coordinates of the points in its closure
                PetscCallA(DMPlexVecGetClosure(dm,coordSection,coordVec,p,coordArray,ierr))
                Do i = 1,dim
                    xyz(i) = sum(coordArray(i:size(coordArray):dim)) * dim / size(coordArray)
                End Do
                PetscCallA(DMPlexVecRestoreClosure(dm,coordSection,coordVec,p,coordArray,ierr))

                PetscCallA(DMPlexVecGetClosure(dm,s,vLoc,p,vArray,ierr))
                Do i = 1,numDof
                    vArray(i) = f(xyz(1),xyz(2)) * 10**i
                End Do
                !!! This is dangerous as I could potentially overwrite the value of 
                !!! v at other points
                !!! In the cirrent state of DMPlexVecSetClosure / DMPlexVecRestoreClosure
                !!! in fortran, it also forces 1 alloc and 1 free 
                PetscCallA(DMPlexVecSetClosure(dm,s,vLoc,p,vArray,INSERT_VALUES,ierr))
                PetscCallA(DMPlexVecRestoreClosure(dm,s,vLoc,p,vArray,ierr))
            End If
        End Do
        PetscCallA(DMLocalToGlobal(dm,vLoc,INSERT_VALUES,v,ierr))
        PetscCallA(DMRestoreLocalVector(dm,vLoc,ierr))
    End subroutine project    

#undef __FUNCT__
#define __FUNCT__ "project2"
    subroutine project2(v,s,f,ierr)
        Type(tVec),intent(IN)              :: v
        Type(tPetscSection),intent(IN)     :: s
        procedure(f_interface)             :: f
        PetscErrorCode,intent(INOUT)       :: ierr

        PetscInt                           :: pStart,pEnd,p,numDof,i
        Type(tDM)                          :: dm
        Type(tPetscSection)                :: coordSection
        Type(tVec)                         :: coordVec,vLoc
        PetscScalar,dimension(:),Pointer   :: coordArray,vArray
        PetscScalar,dimension(3)           :: xyz
        PetscInt                           :: dim

        PetscCallA(PetscSectionGetChart(s,pStart,pEnd,ierr))
        PetscCallA(VecGetDM(v,dm,ierr))
        PetscCallA(DMGetLocalVector(dm,vLoc,ierr))
        PetscCallA(DMGetCoordinateSection(dm,coordSection,ierr))
        PetscCallA(DMGetCoordinatesLocal(dm,coordVec,ierr))
        PetscCallA(DMGetDimension(dm,dim,ierr))

        Do p = pStart,pEnd-1
            PetscCallA(PetscSectionGetDof(s,p,numDof,ierr))
            If (numDof > 0) Then
                !!! trick: the coordinate of a point is the average of the coordinates of the points in its closure
                PetscCallA(DMPlexVecGetClosure(dm,coordSection,coordVec,p,coordArray,ierr))
                Do i = 1,dim
                    xyz(i) = sum(coordArray(i:size(coordArray):dim)) * dim / size(coordArray)
                End Do
                PetscCallA(DMPlexVecRestoreClosure(dm,coordSection,coordVec,p,coordArray,ierr))

                !!! If numdof is constant (where > 0) I could move the allocate
                !!! outside of this loop 
                Allocate(vArray(numDof))
                Do i = 1,numDof
                    vArray(i) = f(xyz(1),xyz(2)) * 10**i
                End Do
                PetscCallA(VecSetValuesSectionF90(vLoc,s,p,vArray,INSERT_VALUES,ierr))
                DeAllocate(vArray)
            End If
        End Do
        PetscCallA(DMLocalToGlobal(dm,vLoc,INSERT_VALUES,v,ierr))
        PetscCallA(DMRestoreLocalVector(dm,vLoc,ierr))
    End subroutine project2    
    
#undef __FUNCT__
#define __FUNCT__ "project3"
    subroutine project3(v,s,f,ierr)
        Type(tVec),intent(IN)              :: v
        Type(tPetscSection),intent(IN)     :: s
        procedure(f_interface)             :: f
        PetscErrorCode,intent(INOUT)       :: ierr

        PetscInt                           :: pStart,pEnd,p,numDof,i
        Type(tDM)                          :: dm
        Type(tPetscSection)                :: coordSection
        Type(tVec)                         :: coordVec,vLoc
        PetscScalar,dimension(:),Pointer   :: coordArray,vArray
        PetscScalar,dimension(3)           :: xyz
        PetscInt                           :: dim,pOffset

        PetscCallA(PetscSectionGetChart(s,pStart,pEnd,ierr))
        PetscCallA(VecGetDM(v,dm,ierr))
        PetscCallA(DMGetLocalVector(dm,vLoc,ierr))
        PetscCallA(DMGetCoordinateSection(dm,coordSection,ierr))
        PetscCallA(DMGetCoordinatesLocal(dm,coordVec,ierr))
        PetscCallA(DMGetDimension(dm,dim,ierr))
        PetscCallA(VecGetArrayF90(vLoc,vArray,ierr))

        !nullify(vArray)
        Do p = pStart,pEnd-1
            PetscCallA(PetscSectionGetDof(s,p,numDof,ierr))
            If (numDof > 0) Then
                !!! trick: the coordinate of a point is the average of the coordinates of the points in its closure
                PetscCallA(DMPlexVecGetClosure(dm,coordSection,coordVec,p,coordArray,ierr))
                Do i = 1,dim
                    xyz(i) = sum(coordArray(i:size(coordArray):dim)) * dim / size(coordArray)
                End Do
                PetscCallA(DMPlexVecRestoreClosure(dm,coordSection,coordVec,p,coordArray,ierr))

                PetscCallA(PetscSectionGetOffset(s,p,pOffset,ierr))
                Do i = 1,numDof
                    vArray(pOffset+i) = f(xyz(1),xyz(2)) * 10**i
                End Do
            End If
        End Do
        PetscCallA(VecRestoreArrayF90(vLoc,vArray,ierr))
        PetscCallA(DMLocalToGlobal(dm,vLoc,INSERT_VALUES,v,ierr))
        PetscCallA(DMRestoreLocalVector(dm,vLoc,ierr))
    End subroutine project3    
End Module localFunctions

Program  TestConstraints
#include <petsc/finclude/petsc.h>
Use m_MEF90
Use petsc
Use localFunctions
Implicit NONE   
    
    PetscErrorCode                      :: ierr
    Type(MEF90Ctx_Type),target          :: MEF90Ctx
    Type(MEF90CtxGlobalOptions_Type)    :: MEF90GlobalOptions_default
    Type(tDM)                           :: dm,dmU,dmU0
    PetscBool                           :: interpolate = PETSC_TRUE
    Character(len=MEF90_MXSTRLEN)       :: IOBuffer
    PetscEnum                           :: setType

    Type(MEF90Element_Type)             :: cellSetElementType,faceSetElementType
    PetscInt                            :: numComponents
    PetscInt                            :: set
    type(tIS)                           :: setIS
    PetscInt,Dimension(:),pointer       :: setID
    PetscInt                            :: dim,pStart,pEnd,order = 1
    PetscBool                           :: flg
    Type(tPetscSection)                 :: sectionU,sectionU0
    Logical,Dimension(:,:),Pointer      :: ConstraintTruthTableU,ConstraintTruthTableU0
    Logical,Dimension(:),Pointer        :: Constraints
    Type(tVec)                          :: U,U0

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
    PetscCallA(DMCreateGlobalVector(dmU,U,ierr))

    PetscCallA(DMClone(dm,dmU0,ierr))
    PetscCallA(DMSetLocalSection(dmU0,sectionU0,ierr))
    PetscCallA(DMCreateGlobalVector(dmU0,U0,ierr))

    PetscCallA(VecSet(U,-1.0_Kr,ierr))
    PetscCallA(VecSet(U0,-1.0_Kr,ierr))
    PetscCallA(project(U,sectionU,f1,ierr))
    PetscCallA(project(U0,sectionU0,f1,ierr))
    PetscCallA(VecViewFromOptions(U,PETSC_NULL_OPTIONS,"-vec_view",ierr))
    PetscCallA(VecViewFromOptions(U0,PETSC_NULL_OPTIONS,"-vec_view",ierr))

    PetscCallA(VecSet(U,-1.0_Kr,ierr))
    PetscCallA(VecSet(U0,-1.0_Kr,ierr))
    PetscCallA(project2(U,sectionU,f1,ierr))
    PetscCallA(project2(U0,sectionU0,f1,ierr))
    PetscCallA(VecViewFromOptions(U,PETSC_NULL_OPTIONS,"-vec_view",ierr))
    PetscCallA(VecViewFromOptions(U0,PETSC_NULL_OPTIONS,"-vec_view",ierr))

    PetscCallA(VecSet(U,-1.0_Kr,ierr))
    PetscCallA(VecSet(U0,-1.0_Kr,ierr))
    PetscCallA(project3(U,sectionU,f1,ierr))
    PetscCallA(project3(U0,sectionU0,f1,ierr))
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
 
       