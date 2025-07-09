Program  TestVecSetBCFromOptionsExpr
#include <petsc/finclude/petsc.h>
Use petsc
Use m_MEF90
#ifdef MEF90_HAVE_SYMENGINEF90
use symengine
#endif

Implicit NONE   
    
    PetscErrorCode                                     :: ierr
    Type(MEF90Ctx_Type),target                         :: MEF90Ctx
    Type(MEF90CtxGlobalOptions_Type)                   :: MEF90GlobalOptions_default
    Type(MEF90CtxGlobalOptions_Type),Pointer           :: MEF90GlobalOptions
    Type(tDM)                                          :: dm
    PetscBool                                          :: flg,interpolate = PETSC_TRUE
    Character(len=MEF90MXSTRLEN)                       :: name
    PetscInt                                           :: sDim = 1

    Type(tVec)                                         :: v
    PetscReal                                          :: scalingFactor

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
    
    Call MEF90Initialize(PETSC_COMM_WORLD,ierr)
    Call MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90GlobalOptions_default,ierr)
    PetscCallA(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr))
    
    PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm,MEF90Ctx%geometryfile,PETSC_NULL_CHARACTER,interpolate,dm,ierr))
    PetscCallA(DMPlexDistributeSetDefault(dm,PETSC_FALSE,ierr))
    PetscCallA(DMSetFromOptions(dm,ierr))
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OBJECT,"-mef90dm_view",ierr))
    
    distribute: Block 
        Type(tDM),target                    :: dmDist
        PetscInt                            :: ovlp = 0
        If (MEF90Ctx%NumProcs > 1) Then
            PetscCallA(DMPlexDistribute(dm,ovlp,PETSC_NULL_SF,dmDist,ierr))
            PetscCallA(DMDestroy(dm,ierr))
            dm = dmDist
        End If
    End Block distribute

    name = "Temperature"
    PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-sdim',sdim,flg,ierr))
    PetscCallA(MEF90CreateLocalVector(dm,MEF90GlobalOptions%elementFamily,MEF90GlobalOptions%elementOrder,sdim,name,V,ierr))
    PetscCallA(VecSet(v,-1.234_Kr,ierr))
    scalingFactor = 1.0_kI
    PetscCallA(MEF90VecSetBCValuesFromOptions(V,scalingFactor,ierr))
    PetscCallA(VecViewFromOptions(V,PETSC_NULL_OBJECT,"-temperature_vec_view",ierr))

    PetscCallA(VecSet(v,5.678_Kr,ierr))
    scalingFactor = 1.0_kI
    PetscCallA(MEF90VecSetValuesFromOptions(V,scalingFactor,ierr))
    PetscCallA(VecViewFromOptions(V,PETSC_NULL_OBJECT,"-temperature_vec_view",ierr))

    ViewSec: block
        Type(tPetscSection)     :: sectionV
        Type(tDM)               :: dmV

        PetscCallA(VecGetDM(V,dmV,ierr))
        PetscCallA(DMGetLocalSection(dmV,sectionV,ierr))
        PetscCallA(PetscSectionViewFromOptions(sectionV,PETSC_NULL_OBJECT,"-temperature_section_view",ierr))
    end block ViewSec


    PetscCallA(VecDestroy(v,ierr))
    PetscCallA(DMDestroy(dm,ierr))
    Call MEF90CtxDestroy(MEF90Ctx,ierr)   
    Call MEF90Finalize(ierr)
    Call PetscFinalize(ierr)

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90VecSetValuesFromOptionsExpr"
!!!
!!!  
!!!  MEF90VecSetValuesFromOptionsExpr: Fill values of a Vec using command line options
!!!  
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

    Subroutine MEF90VecSetValuesFromOptionsExpr(v,t,ierr)
        Type(tVec),Intent(INOUT)                 :: v
        PetscReal,Intent(IN)                     :: t
        PetscErrorCode,Intent(INOUT)             :: ierr

#ifdef MEF90_HAVE_SYMENGINEF90
        Type(tDM)                                :: dm
        PetscEnum                                :: setType
        PetscInt                                 :: set,point
        Type(tIS)                                :: setIS,pointIS
        PetscInt,Dimension(:),Pointer            :: setID,pointID
        Character(len=MEF90MXSTRLEN)             :: ValueKey,name,ExprStr
        Character(len=MEF90MXSTRLEN),allocatable :: ExprStrComp
        PetscBool                                :: flg
        PetscInt                                 :: dim,numOpt,bs,numDofClosure,i,c
        PetscReal,Dimension(:),pointer           :: Val,vArray
        Type(tPetscSection)                      :: section
        Type(tPetscSection)                      :: coordSection
        Type(tVec)                               :: coordVec
        PetscReal,dimension(:),Pointer           :: coordArray
        PetscReal,dimension(3)                   :: xyz
        type(Basic),dimension(:),allocatable     :: exprs
        type(Basic)                              :: tmpExpr
        type(symbol),dimension(3)                :: vars = [Symbol("x"), Symbol("y"), Symbol("z")]

        PetscCall(VecGetDM(v,dm,ierr))
        PetscCall(PetscObjectGetName(v,name,ierr))
        PetscCall(DMGetLocalSection(dm,section,ierr))

        PetscCall(DMGetDimension(dm,dim,ierr))
        PetscCall(VecGetBlockSize(v,bs,ierr))
        Allocate(Val(bs))

        PetscCallA(DMGetCoordinateSection(dm,coordSection,ierr))
        PetscCallA(DMGetCoordinatesLocal(dm,coordVec,ierr))
        allocate(exprs(dim))

        vals[1] = RealDouble(t)
        Do setType = 1,size(MEF90SetType)
            PetscCall(DMGetLabelIdIS(dm,MEF90SetLabelName(setType),setIS,ierr))
            ! PetscCall(MEF90ISAllGatherMerge(comm,setIS,ierr))
            If (.NOT. PetscObjectIsNull(setIS)) Then
                PetscCall(ISGetIndices(setIS,setID,ierr))
                Do set = 1,size(setID)
                    write(ValueKey,'("-",a2,I4.4,"_",a)') MEF90SetPrefix(setType),setID(set),trim(name)
                    numOpt = bs
                    PetscCall(PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,trim(ValueKey),ExprStr,flg,ierr))
                    numOpts = MEF90StrCount(ExprStr,',')
                    call MEF90StrTokenize(ExprStr,',',ExprStrComp)
                    if (shape(ExprStrComp) < dim) then
                        Write(IOBuffer,"(A,' was expecting ',I2,' expressions but got only ',I2,' for key ',A,'\n')") __FUNCT__, dim, shape(ExprStrComp), ValueKey
                        SETERRQ(MEF90Ctx%Comm,PETSC_ERR_ARG_SIZ,IOBuffer)
                    end if
                    !!! create parsers
                    allocate(exprs(shape(ExprStrComp)))
                    do c = 1, shape(ExprStrComp)
                        exprs(c) = parse(ExprStrComp(c))
                    end do
                    
                    If (numOpt > 0) Then
                        PetscCall(DMGetStratumIS(dm,MEF90SetLabelName(setType),setID(set),pointIS,ierr))
                        !!! Set the values on the closure of the current point
                        If (.NOT. PetscObjectIsNull(pointIS)) Then
                            PetscCall(ISGetIndices(pointIS,pointID,ierr))
                            Do point = 1, size(pointID)
                                PetscCall(MEF90VecGetClosureSize(v,pointID(point),numDofClosure,ierr))
                                    If (numDofClosure > 0) Then
                                    !!! Get the coordinates of the dof associated with the point
                                    !!! trick: the coordinate of a point is the average of the coordinates of the points in its closure
                                    PetscCallA(DMPlexVecGetClosure(dm,coordSection,coordVec,p,PETSC_NULL_INTEGER,coordArray,ierr))
                                    Do i = 1,dim
                                        xyz(i) = sum(coordArray(i:size(coordArray):dim)) * dim / size(coordArray)
                                    End Do
                                    PetscCallA(DMPlexVecRestoreClosure(dm,coordSection,coordVec,p,PETSC_NULL_INTEGER,coordArray,ierr))
                                    PetscCall(DMPlexVecGetClosure(dm,section,v,pointID(point),PETSC_NULL_INTEGER,vArray,ierr))

                                    !!! Parse expressions
                                    do c = 1, numDofClosure
                                        tmpExpr = Exprs(c)
                                        tmpExpr = tmpExpr%subs(Symbol("t"),RealDouble(t))
                                        do i = 1, shape(ExprStrComp)
                                            tmpExpr = tmpExpr%subs(vars(i),RealDouble(xyz(i)))
                                        end do
                                        tmpExpr = tmpExpr%evalf()
                                        vArray(c) = tmpExpr%dbl()
                                    end do
                                    PetscCall(DMPlexVecSetClosure(dm,section,v,pointID(point),vArray,INSERT_ALL_VALUES,ierr))
                                    PetscCall(DMPlexVecRestoreClosure(dm,section,v,pointID(point),PETSC_NULL_INTEGER,vArray,ierr))
                                End If ! numDofClosure
                            End Do ! point
                            PetscCall(ISRestoreIndices(pointIS,pointID,ierr))
                        End If ! pointIS
                        PetscCall(ISDestroy(pointIS,ierr))
                    End If ! numOpt
                    deAllocate(exprs)
                    DeAllocate(ExprStrComp)
                End Do ! set
                PetscCall(ISRestoreIndices(setIS,setID,ierr))
            End If ! setIS
            PetscCall(ISDestroy(setIS,ierr))
        End Do ! setType
        DeAllocate(exprs)
        DeAllocate(Val)
#else
    Write(*,*) "ERROR: ",__FUNC__, " requires symengine-f90 support"
    STOP
#endif
    End Subroutine MEF90VecSetValuesFromOptionsExpr


End Program  TestVecSetBCFromOptionsExpr

!!! TEST
!!! ./TestVecSetBCFromOptions -geometry ../TestMeshes/SquareFaceSet.gen -temperature_section_view    \
!!!     -temperature_vec_view -sdim 2 -fs0021_temperatureBC yes,no -fs0021_boundaryTemperature 10,20 \
!!!     -vs0010_temperatureBC yes,yes -vs0010_boundaryTemperature 100,101 -cs0001_temperature 1000   \
!!!     -vs0010_temperature -2 -fs0020_temperature 20 -fs0021_temperature 21 -vs0010_temperature 100 
