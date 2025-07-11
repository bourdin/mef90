module TestExpr
#include <petsc/finclude/petsc.h>
Use petsc
Use m_MEF90
Use m_MEF90_DMPlex
Use symengine

Implicit NONE   
Contains

#undef __FUNCT__
#define __FUNCT__ "myMEF90VecSetValuesFromOptionsExpr"
!!!
!!!  
!!!  myMEF90VecSetValuesFromOptionsExpr: Fill values of a Vec using expressions passed as petsc options
!!!  
!!!  (c) 2025      Blaise Bourdin bourdin@mcmaster.ca
!!!

    Subroutine myMEF90VecSetValuesFromOptionsExpr(v,t,ierr)
        Type(tVec),Intent(INOUT)                 :: v
        PetscReal,Intent(IN)                     :: t
        PetscErrorCode,Intent(INOUT)             :: ierr

#ifdef MEF90_HAVE_SYMENGINEF90
        Type(tDM)                                :: dm
        PetscEnum                                :: setType
        PetscInt                                 :: set,point,p
        Type(tIS)                                :: setIS,pointIS
        PetscInt,Dimension(:),Pointer            :: setID,pointID
        Character(len=MEF90MXSTRLEN)             :: ValueKey,name,ExprStr,IOBuffer
        PetscBool                                :: flg
        PetscInt                                 :: dim,numOpt,bs,numDofClosure,numDof,i,c,dof
        PetscInt,Dimension(:),Pointer            :: closure
        PetscReal,Dimension(:),pointer           :: vArray
        Type(tPetscSection)                      :: section
        Type(tPetscSection)                      :: coordSection
        Type(tVec)                               :: coordVec
        PetscReal,dimension(:),Pointer           :: coordArray
        PetscReal,dimension(3)                   :: xyz
            type(Basic),dimension(:),allocatable     :: exprs
        type(Basic)                              :: tmpExpr
        type(symbol),dimension(3)                :: vars
        type(realdouble),dimension(3)            :: vals
        character                                :: delim = ';'
        Character(len=MEF90MXSTRLEN),dimension(:),allocatable :: ExprStrComp

        PetscCall(VecGetDM(v,dm,ierr))
        PetscCall(PetscObjectGetName(v,name,ierr))
        PetscCall(DMGetLocalSection(dm,section,ierr))

        PetscCall(DMGetDimension(dm,dim,ierr))
        PetscCall(VecGetBlockSize(v,bs,ierr))

        vars = [Symbol("x"), Symbol("y"), Symbol("z")]
        allocate(exprs(bs))

        PetscCall(DMGetCoordinateSection(dm,coordSection,ierr))
        PetscCall(DMGetCoordinatesLocal(dm,coordVec,ierr))

        vals(1) = RealDouble(t)
        Do setType = 1,size(MEF90SetType)
            PetscCall(DMGetLabelIdIS(dm,MEF90SetLabelName(setType),setIS,ierr))
            ! PetscCall(MEF90ISAllGatherMerge(comm,setIS,ierr))
            If (.NOT. PetscObjectIsNull(setIS)) Then
                PetscCall(ISGetIndices(setIS,setID,ierr))
                Do set = 1,size(setID)
                    write(ValueKey,'("-",a2,I4.4,"_",a)') MEF90SetPrefix(setType),setID(set),trim(name)
                    PetscCall(PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,trim(ValueKey),ExprStr,flg,ierr))
                    if (flg) then
                        call MEF90StrTokenize(ExprStr,delim,ExprStrComp)
                        numOpt = size(ExprStrComp)
                        if (numOpt /= bs) then
                            Write(IOBuffer,"(A,' was expecting ',I2,' expressions but got ',I2,' for key ',A)") __FUNCT__, bs, size(ExprStrComp), trim(ValueKey)
                            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,IOBuffer)
                        end if
                        !!! create parsers
                        do c = 1, bs
                            exprs(c) = parse(ExprStrComp(c))
                            exprs(c) = exprs(c)%subs(Symbol("t"),RealDouble(t))
                        end do
                        PetscCall(DMGetStratumIS(dm,MEF90SetLabelName(setType),setID(set),pointIS,ierr))
                        !!! Set the values on the closure of the current point
                        If (.NOT. PetscObjectIsNull(pointIS)) Then
                            PetscCall(ISGetIndices(pointIS,pointID,ierr))
                            Do point = 1, size(pointID)
                                PetscCall(DMPlexGetTransitiveClosure(dm,pointID(point),PETSC_TRUE,numDofClosure,closure,ierr))
                                If (numDofClosure > 0) Then
                                    PetscCall(DMPlexVecGetClosure(dm,section,v,pointID(point),PETSC_NULL_INTEGER,vArray,ierr))
                                    dof = 0
                                    Do p = 1, size(closure), 2
                                        PetscCall(PetscSectionGetDof(section,closure(p),numDof,ierr))
                                        if (numDof > 0) Then
                                            dof = dof + 1
                                            !!! Get the coordinates of the dof associated with the point
                                            !!! trick: the coordinate of a point is the average of the coordinates of the points in its closure
                                            PetscCall(DMPlexVecGetClosure(dm,coordSection,coordVec,closure(p),PETSC_NULL_INTEGER,coordArray,ierr))
                                            Do i = 1,dim
                                                xyz(i) = sum(coordArray(i:size(coordArray):dim)) * dim / size(coordArray)
                                            End Do
                                            PetscCall(DMPlexVecRestoreClosure(dm,coordSection,coordVec,closure(p),PETSC_NULL_INTEGER,coordArray,ierr))
                                                do c = 1, bs
                                                tmpExpr = Exprs(c)
                                                do i = 1, dim
                                                    tmpExpr = tmpExpr%subs(vars(i),RealDouble(xyz(i)))
                                                end do ! i
                                                tmpExpr = tmpExpr%evalf()
                                                vArray((dof - 1) * bs + c) = tmpExpr%dbl()
                                            end do ! c
                                        end if !bnumDof
                                    End Do !p
                                    PetscCall(DMPlexVecSetClosure(dm,section,v,pointID(point),vArray,INSERT_ALL_VALUES,ierr))
                                    PetscCall(DMPlexVecRestoreClosure(dm,section,v,pointID(point),PETSC_NULL_INTEGER,vArray,ierr))
                                End If ! numDofClosure
                                PetscCall(DMPlexRestoreTransitiveClosure(dm,pointID(point),PETSC_TRUE,numDofClosure,closure,ierr))
                            End Do ! point
                            PetscCall(ISRestoreIndices(pointIS,pointID,ierr))
                        End If ! pointIS
                        PetscCall(ISDestroy(pointIS,ierr))
                        deAllocate(ExprStrComp)
                    end if ! flg
                End Do ! set
                PetscCall(ISRestoreIndices(setIS,setID,ierr))
            End If ! setIS
            PetscCall(ISDestroy(setIS,ierr))
        End Do ! setType
        DeAllocate(exprs)
#else
    Write(*,*) "ERROR: ",__FUNCT__, " requires symengine-f90 support"
    STOP
#endif
    End Subroutine myMEF90VecSetValuesFromOptionsExpr

#undef __FUNCT__
#define __FUNCT__ "myMEF90VecSetBCValuesFromOptionsExpr"
!!!
!!!  
!!!  myMEF90VecSetBCValuesFromOptionsExpr: Fill boundary values of a Vec using expressions passed as petsc options
!!!  
!!!  (c) 2025      Blaise Bourdin bourdin@mcmaster.ca
!!!

    Subroutine myMEF90VecSetBCValuesFromOptionsExpr(v,t,ierr)
        Type(tVec),Intent(INOUT)                 :: v
        PetscReal,Intent(IN)                     :: t
        PetscErrorCode,Intent(INOUT)             :: ierr

#ifdef MEF90_HAVE_SYMENGINEF90
        Type(tDM)                                :: dm
        PetscEnum                                :: setType
        PetscInt                                 :: set,point,c,p,i
        Type(tIS)                                :: setIS,pointIS
        PetscInt,Dimension(:),Pointer            :: setID,pointID
        Character(len=MEF90MXSTRLEN)             :: BCOptionKey,BCValueKey,name,ExprStr,IOBuffer
        PetscBool,Dimension(:),Pointer           :: setBC
        PetscBool                                :: flg
        PetscInt                                 :: dim,numOpt,numBC,bs,numDofClosure,numDof,dof
        PetscInt,Dimension(:),Pointer            :: closure
        PetscReal,Dimension(:),pointer           :: vArray
        Type(tPetscSection)                      :: section
        Type(tPetscSection)                      :: coordSection
        Type(tVec)                               :: coordVec
        PetscReal,dimension(:),Pointer           :: coordArray
        PetscReal,dimension(3)                   :: xyz
        type(Basic),dimension(:),allocatable     :: exprs
        type(Basic)                              :: tmpExpr
        type(symbol),dimension(3)                :: vars
        type(realdouble),dimension(3)            :: vals
        character                                :: delim = ';'
        Character(len=MEF90MXSTRLEN),dimension(:),allocatable :: ExprStrComp

        PetscCall(VecGetDM(v,dm,ierr))
        PetscCall(PetscObjectGetName(v,name,ierr))
        PetscCall(DMGetLocalSection(dm,section,ierr))

        PetscCall(DMGetDimension(dm,dim,ierr))
        PetscCall(VecGetBlockSize(v,bs,ierr))
        Allocate(setBC(bs))

        vars = [Symbol("x"), Symbol("y"), Symbol("z")]
        allocate(exprs(bs))

        PetscCall(DMGetCoordinateSection(dm,coordSection,ierr))
        PetscCall(DMGetCoordinatesLocal(dm,coordVec,ierr))      

        vals(1) = RealDouble(t)
        Do setType = 1,size(MEF90SetType)
            PetscCall(DMGetLabelIdIS(dm,MEF90SetLabelName(setType),setIS,ierr))
            ! PetscCall(MEF90ISAllGatherMerge(comm,setIS,ierr))
            If (.NOT. PetscObjectIsNull(setIS)) Then
                PetscCall(ISGetIndices(setIS,setID,ierr))
                Do set = 1,size(setID)
                    setBC = .FALSE.
                    write(BCOptionKey,'("-",a2,I4.4,"_",a,"BC")') MEF90SetPrefix(setType),setID(set),trim(name)
                    numBC = bs
                    PetscCall(PetscOptionsGetBoolArray(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,trim(BCOptionKey),setBC,numBC,flg,ierr))
                    If (any(setBC)) Then
                        !!! At least 1 dof has a boundary condition
                        !!! Get the unit BC value on the set
                        write(BCValueKey,'("-",a2,I4.4,"_Boundary",a)') MEF90SetPrefix(setType),setID(set),trim(name)
                        PetscCall(PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,trim(BCValueKey),ExprStr,flg,ierr))
                        numOpt = MEF90StrCount(ExprStr,delim)
                        call MEF90StrTokenize(ExprStr,delim,ExprStrComp)
                        if (size(ExprStrComp) /= bs) then
                            Write(IOBuffer,"(A,' was expecting ',I2,' expressions but got ',I2,' for key ',A)") __FUNCT__, bs, size(ExprStrComp), trim(BCValueKey)
                            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,IOBuffer)
                        end if
                        !!! create parsers
                        do c = 1, bs
                            exprs(c) = parse(ExprStrComp(c))
                        end do

                        PetscCall(DMGetStratumIS(dm,MEF90SetLabelName(setType),setID(set),pointIS,ierr))
                        !!! Set the boundary values on the closure of the current point
                        If (.NOT. PetscObjectIsNull(pointIS)) Then
                            PetscCall(ISGetIndices(pointIS,pointID,ierr))
                            Do point = 1, size(pointID)
                                PetscCall(DMPlexGetTransitiveClosure(dm,pointID(point),PETSC_TRUE,numDofClosure,closure,ierr))
                                If (numDofClosure > 0) Then
                                    PetscCall(DMPlexVecGetClosure(dm,section,v,pointID(point),PETSC_NULL_INTEGER,vArray,ierr))
                                    dof = 0
                                    Do p = 1, size(closure), 2
                                        PetscCall(PetscSectionGetDof(section,closure(p),numDof,ierr))
                                        if (numDof > 0) Then
                                            dof = dof + 1
                                            !!! Get the coordinates of the dof associated with the point
                                            !!! trick: the coordinate of a point is the average of the coordinates of the points in its closure
                                            PetscCall(DMPlexVecGetClosure(dm,coordSection,coordVec,closure(p),PETSC_NULL_INTEGER,coordArray,ierr))
                                            Do i = 1,dim
                                                xyz(i) = sum(coordArray(i:size(coordArray):dim)) * dim / size(coordArray)
                                            End Do
                                            PetscCall(DMPlexVecRestoreClosure(dm,coordSection,coordVec,closure(p),PETSC_NULL_INTEGER,coordArray,ierr))
                                            do c = 1, bs
                                                if (setBC(c)) Then
                                                    tmpExpr = Exprs(c)
                                                    do i = 1, dim
                                                        tmpExpr = tmpExpr%subs(vars(i),RealDouble(xyz(i)))
                                                    end do ! i
                                                    tmpExpr = tmpExpr%evalf()
                                                    vArray((dof - 1) * bs + c) = tmpExpr%dbl()
                                                end if
                                            end do ! c
                                        end if !bnumDof
                                    End Do !p
                                    PetscCall(DMPlexVecSetClosure(dm,section,v,pointID(point),vArray,INSERT_ALL_VALUES,ierr))
                                    PetscCall(DMPlexVecRestoreClosure(dm,section,v,pointID(point),PETSC_NULL_INTEGER,vArray,ierr))
                                End If ! numDofClosure
                                PetscCall(DMPlexRestoreTransitiveClosure(dm,pointID(point),PETSC_TRUE,numDofClosure,closure,ierr))
                            End Do ! point
                            PetscCall(ISRestoreIndices(pointIS,pointID,ierr))
                        End If ! pointIS
                        PetscCall(ISDestroy(pointIS,ierr))
                        deAllocate(ExprStrComp)
                    End If ! setBC
                End Do ! set
                PetscCall(ISRestoreIndices(setIS,setID,ierr))
            End If ! setIS
            PetscCall(ISDestroy(setIS,ierr))
        End Do ! setType
        DeAllocate(setBC)
        DeAllocate(exprs)
#else
    Write(*,*) "ERROR: ",__FUNCT__, " requires symengine-f90 support"
    STOP
#endif
    End Subroutine myMEF90VecSetBCValuesFromOptionsExpr

end module TestExpr

Program  TestVecSetBCFromOptionsExpr
#include <petsc/finclude/petsc.h>
Use petsc
Use m_MEF90
Use m_MEF90_DMPlex
Use TestExpr

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
    
    PetscCallA(VecSet(v,5.678_Kr,ierr))
    scalingFactor = 2.0_Kr
    PetscCallA(MEF90VecSetValuesFromOptionsExpr(V,scalingFactor,ierr))
    PetscCallA(VecViewFromOptions(V,PETSC_NULL_OBJECT,"-temperature_vec_view",ierr))

    PetscCallA(VecSet(v,-1.234_Kr,ierr))
    scalingFactor = 1.0_Kr
    PetscCallA(MEF90VecSetBCValuesFromOptionsExpr(V,scalingFactor,ierr))
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
End Program  TestVecSetBCFromOptionsExpr

!!! TEST
!!! ./TestVecSetBCFromOptions -geometry ../TestMeshes/SquareFaceSet.gen -temperature_section_view    \
!!!     -temperature_vec_view -sdim 2 -fs0021_temperatureBC yes,no -fs0021_boundaryTemperature 10,20 \
!!!     -vs0010_temperatureBC yes,yes -vs0010_boundaryTemperature 100,101 -cs0001_temperature 1000   \
!!!     -vs0010_temperature -2 -fs0020_temperature 20 -fs0021_temperature 21 -vs0010_temperature 100 
