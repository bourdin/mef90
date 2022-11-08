module mef90Projection
#include <petsc/finclude/petsc.h>
use petsc
use m_MEF90
implicit none

    abstract interface
        function f_interface(x,y,z)
            PetscReal,intent(in) :: x,y,z
            PetscReal            :: f_interface
        end function f_interface
    end interface
        
Contains
#undef __FUNCT__
#define __FUNCT__ "f1"

    function f1(x,y,z)
        PetscReal,intent(in)   :: x,y,z
        PetscReal              :: f1

        PetscInt               :: i,j,k
        PetscBool              :: flg
        PetscErrorCode         :: ierr

        i = 0_Ki
        j = 0_Ki
        k = 0_Ki
        PetscCall(PetscOptionsGetInt(PETSC_NULL_OPTIONS,'','-i',i,flg,ierr))
        PetscCall(PetscOptionsGetInt(PETSC_NULL_OPTIONS,'','-j',j,flg,ierr))
        PetscCall(PetscOptionsGetInt(PETSC_NULL_OPTIONS,'','-k',k,flg,ierr))
        f1 = x**i * y**j * z**k
        !f1 = 1-x-y-z
    end function f1
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
        Type(tVec)                         :: coordVec
        PetscScalar,dimension(:),Pointer   :: coordArray,vArray
        PetscScalar,dimension(3)           :: xyz
        PetscInt                           :: dim

        PetscCallA(PetscSectionGetChart(s,pStart,pEnd,ierr))
        PetscCallA(VecGetDM(v,dm,ierr))
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
                    vArray(i) = f(xyz(1),xyz(2),xyz(3)) * 10**(i-1)
                End Do
                PetscCallA(VecSetValuesSectionF90(v,s,p,vArray,INSERT_ALL_VALUES,ierr))
                !!! This call drops the constrained values if the mode is INSERT_VALUES
                DeAllocate(vArray)
            End If
        End Do
    End subroutine project    
end module mef90Projection
Program  TestMassMatrix
#include <petsc/finclude/petsc.h>
    Use m_MEF90
    Use petsc
    Use mef90Projection
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
    Type(tVec)                                     :: U,U0,One,One0
    Type(tDM)                                      :: dmU,dmU0
    Type(MEF90CtxGlobalOptions_Type),pointer       :: MEF90GlobalOptions
    Type(MEF90Element2DScal),dimension(:),Pointer  :: elem2D
    Type(MEF90Element3DScal),dimension(:),Pointer  :: elem3D
    PetscInt                                       :: quadratureOrder
    Type(tIS)                                      :: setPointIS
    PetscInt,Dimension(:),Pointer                  :: setPointID
    Type(MEF90ElementType)                         :: elementType
    DMPolytopeType                                 :: cellType,faceType
    PetscReal                                      :: myL1NormSet,L1NormSet,L1Norm
    PetscReal                                      :: myL2NormSet,L2NormSet,L2Norm
    PetscReal                                      :: myH1NormSet,H1NormSet,H1Norm
    Type(tPetscSection)                            :: sectionU

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
    
    name = "U"
    PetscCallA(MEF90CreateLocalVector(dm,MEF90GlobalOptions%elementFamily,MEF90GlobalOptions%elementOrder,1_Ki,name,U,ierr))
    PetscCallA(VecGetDM(U,dmU,ierr))
    PetscCallA(DMGetLocalSection(dmU,sectionU,ierr))
    PetscCallA(project(U,sectionU,f1,ierr))
    PetscCallA(VecViewFromOptions(U,PETSC_NULL_OPTIONS,"-mef90VecU_view",ierr))
    PetscCallA(VecDuplicate(U,One,ierr))
    PetscCallA(VecSet(One,1.0_Kr,ierr))

    name = "U0"
    PetscCallA(MEF90CreateBoundaryLocalVector(dm,MEF90GlobalOptions%elementFamily,MEF90GlobalOptions%elementOrder,1_Ki,name,U0,ierr))
    PetscCallA(VecGetDM(U0,dmU0,ierr))
    PetscCallA(VecSet(U0,1.0_Kr,ierr))
    PetscCallA(VecDuplicate(U0,One0,ierr))
    PetscCallA(VecSet(One0,1.0_Kr,ierr))

    setType = MEF90CellSetType
    PetscCallA(DMGetLabelIdIS(dmU,MEF90SetLabelName(setType),setIS,ierr))
    PetscCall(MEF90ISAllGatherMerge(MEF90Ctx%Comm,setIS,ierr))
    L1Norm    = 0.0_Kr
    L2Norm    = 0.0_Kr
    H1Norm    = 0.0_Kr
    If (setIS /= PETSC_NULL_IS) Then
        PetscCallA(ISGetIndicesF90(setIS,setID,ierr))
        Do set = 1,size(setID)
            myL1NormSet = 0.0_Kr
            myL2NormSet = 0.0_Kr
            myH1NormSet = 0.0_Kr
            PetscCallA(DMGetStratumIS(dmU,MEF90SetLabelName(setType),setID(set),setPointIS,ierr))
            PetscCallA(ISGetIndicesF90(setPointIS,setPointID,ierr))
            PetscCallA(DMPlexGetCellType(dmU,setPointID(1),cellType,ierr))
            PetscCallA(MEF90ElementGetType(MEF90GlobalOptions%elementFamily,MEF90GlobalOptions%elementOrder,cellType,elementType,ierr))
            quadratureOrder = elementType%order * 2
            If (dim == 2) Then
                PetscCallA(MEF90ElementCreate(dmU,setPointIS,elem2D,QuadratureOrder,elementType,ierr))
                PetscCallA(MEF90L2DotProductSet(myL1NormSet,One,U,setType,setID(set),elem2D,elementType,ierr))
                PetscCallA(MEF90L2NormSet(myL2NormSet,U,setType,setID(set),elem2D,elementType,ierr))
                PetscCallA(MEF90H1DotProductSet(myH1NormSet,U,U,setType,setID(set),elem2D,elementType,ierr))
                PetscCallA(MEF90ElementDestroy(elem2D,ierr))
            Else
                PetscCallA(MEF90ElementCreate(dmU,setPointIS,elem3D,QuadratureOrder,elementType,ierr))
                PetscCallA(MEF90L2DotProductSet(myL1NormSet,One,U,setType,setID(set),elem3D,elementType,ierr))
                PetscCallA(MEF90L2NormSet(myL2NormSet,U,setType,setID(set),elem3D,elementType,ierr))
                PetscCallA(MEF90H1DotProductSet(myH1NormSet,U,U,setType,setID(set),elem3D,elementType,ierr))
                PetscCallA(MEF90ElementDestroy(elem3D,ierr))
            End If
            PetscCallA(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
            PetscCallA(ISDestroy(setPointIS,ierr))
            Call MPI_AllReduce(myL1NormSet,L1NormSet,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD,ierr)
            Call MPI_AllReduce(myL2NormSet,L2NormSet,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD,ierr)
            Call MPI_AllReduce(myH1NormSet,H1NormSet,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD,ierr)
            L1Norm = L1Norm + L1NormSet
            L2Norm = L2Norm + L2NormSet
            H1Norm = H1Norm + H1NormSet
            ! Write(IOBuffer,'("set ",I4," L1 norm ", ES12.5,"\n")') setID(set), L1NormSet
            ! PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
            ! Write(IOBuffer,'("set ",I4," L2 norm ", ES12.5,"\n")') setID(set), L2NormSet
            ! PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
            ! Write(IOBuffer,'("set ",I4," H1 norm ", ES12.5,"\n")') setID(set), H1NormSet
            ! PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
        End Do
        PetscCallA(ISRestoreIndicesF90(setIS,setID,ierr))
    End If ! setIS
    PetscCallA(ISDestroy(setIS,ierr))
    Write(IOBuffer,'("L1 norm          ", ES12.5,"\n")') L1Norm
    PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
    Write(IOBuffer,'("L2 norm          ", ES12.5,"\n")') L2Norm
    PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
    Write(IOBuffer,'("H1 norm          ", ES12.5,"\n")') H1Norm
    PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))

    ! setType = MEF90FaceSetType
    ! PetscCallA(DMGetLabelIdIS(dmU0,MEF90SetLabelName(setType),setIS,ierr))
    ! L2Norm    = 0.0_Kr
    ! If (setIS /= PETSC_NULL_IS) Then
    !     PetscCallA(ISGetIndicesF90(setIS,setID,ierr))
    !     QuadratureOrder = elementType%order * 2
    !     Do set = 1,size(setID)
    !         myL2NormSet = 0.0_Kr
    !         PetscCallA(DMGetStratumIS(dmU0,MEF90SetLabelName(setType),setID(set),setPointIS,ierr))
    !         PetscCallA(ISGetIndicesF90(setPointIS,setPointID,ierr))
    !         PetscCall(DMPlexGetCellType(dmU0,setPointID(1),faceType,ierr))
    !         PetscCall(MEF90ElementGetTypeBoundary(MEF90GlobalOptions%elementFamily,MEF90GlobalOptions%elementOrder,faceType,elementType,ierr))
    !             If (dim == 2) Then
    !             PetscCallA(MEF90ElementCreate(dmU0,setPointIS,elem2D,QuadratureOrder,elementType,ierr))
    !             PetscCallA(MEF90_MassMatrixAssembleSet(M0,dmU0,setType,setID(set),elem2D,elementType,ierr))
    !             PetscCallA(MEF90L2NormSet(myL2NormSet,U,setType,setID(set),elem2D,elementType,ierr))
    !             PetscCallA(MEF90ElementDestroy(elem2D,ierr))
    !         Else
    !             PetscCallA(MEF90ElementCreate(dmU0,setPointIS,elem3D,QuadratureOrder,elementType,ierr))
    !             PetscCallA(MEF90_MassMatrixAssembleSet(M0,dmU0,setType,setID(set),elem3D,elementType,ierr))
    !             PetscCallA(MEF90L2NormSet(myL2NormSet,U,setType,setID(set),elem3D,elementType,ierr))
    !             PetscCallA(MEF90ElementDestroy(elem3D,ierr))
    !         End If
    !         PetscCallA(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
    !         PetscCallA(ISDestroy(setPointIS,ierr))
    !         Call MPI_AllReduce(myL2NormSet,L2NormSet,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD,ierr)
    !         L2Norm = L2Norm + L2NormSet
    !         Write(IOBuffer,'("set ",I4," L2 norm ", ES12.5,"\n")') setID(set), L2NormSet
    !         PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
    !     End Do
    !     PetscCallA(ISRestoreIndicesF90(setIS,setID,ierr))
    ! End If ! setIS
    ! PetscCallA(ISDestroy(setIS,ierr))
    ! PetscCallA(MatAssemblyBegin(M0,MAT_FINAL_ASSEMBLY,ierr))
    ! PetscCallA(MatAssemblyEnd(M0,MAT_FINAL_ASSEMBLY,ierr))
    ! Write(IOBuffer,'("L2 norm          ", ES12.5,"\n")') L2Norm
    ! PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
    ! PetscCallA(MatViewFromOptions(M0,PETSC_NULL_OPTIONS,"-mef90matM0_view",ierr))


    PetscCallA(VecDestroy(U,ierr))
    PetscCallA(VecDestroy(U0,ierr))
    PetscCallA(VecDestroy(One,ierr))
    PetscCallA(VecDestroy(One0,ierr))
    PetscCallA(DMDestroy(dm,ierr))

    Call MEF90CtxDestroy(MEF90Ctx,ierr)   
    Call MEF90Finalize(ierr)
    PetscCallA(PetscFinalize(ierr))
End Program TestMassMatrix
           
    