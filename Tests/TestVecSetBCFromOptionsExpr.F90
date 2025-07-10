Program  TestVecSetBCFromOptionsExpr
#include <petsc/finclude/petsc.h>
Use petsc
Use m_MEF90
Use m_MEF90_DMPlex

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
    scalingFactor = 1.0_Kr
    PetscCallA(MEF90VecSetBCValuesFromOptionsExpr(V,scalingFactor,ierr))
    PetscCallA(VecViewFromOptions(V,PETSC_NULL_OBJECT,"-temperature_vec_view",ierr))

    PetscCallA(VecSet(v,5.678_Kr,ierr))
    scalingFactor = 2.0_Kr
    PetscCallA(MEF90VecSetValuesFromOptionsExpr(V,scalingFactor,ierr))
    PetscCallA(VecViewFromOptions(V,PETSC_NULL_OBJECT,"-temperature_vec_view",ierr))

    ViewSec: block
        Type(tPetscSection)     :: sectionV
        Type(tDM)               :: dmV

        PetscCallA(VecGetDM(V,dmV,ierr))
        PetscCallA(DMGetLocalSection(dmV,sectionV,ierr))
        PetscCallA(PetscSectionViewFromOptions(sectionV,PETSC_NULL_OBJECT,"-temperature_section_view",ierr))
    end block ViewSec

    ! test: block
    !     character(len=32) :: str
    !     character         :: delim
    !     character(len=2)  :: delim2
    !     character(len=MEF90MXSTRLEN),dimension(:),allocatable :: tok
    !     integer :: i

    !     str = 'x^2 - 2x + 1, x+y,z'
    !     delim = ','
    !     delim2 = ', '
    !     write(*,*) MEF90StrCount(str,delim)
    !     write(*,*) MEF90StrCount(str,delim2)
    !     write(*,*) MEF90StrCount(trim(str),' ')

    !     write(*,*) "MEF90StrTokenize(str,delim,tok)"
    !     call MEF90StrTokenize(str,delim,tok)
    !     do i = 1, MEF90StrCount(str,delim)+1
    !         write(*,*) i, trim(tok(i))
    !     end do
    !     deallocate(tok)

    !     write(*,*) "MEF90StrTokenize(str,delim2,tok)"
    !     call MEF90StrTokenize(str,delim2,tok)
    !     do i = 1, MEF90StrCount(str,delim2)+1
    !         write(*,*) i, trim(tok(i))
    !     end do
    !     deallocate(tok)

    !     write(*,*) "MEF90StrTokenize(str,',',tok)"
    !     call MEF90StrTokenize(str,',',tok)
    !     do i = 1, MEF90StrCount(str,delim2)+1
    !         write(*,*) i, trim(tok(i))
    !     end do
    ! end block test

    ! test2: block
    !     integer :: i
    !     write(*,*) MEF90CellSetType
    !     write(*,*) 'size(MEF90SetType) ', size(MEF90SetType)
    !     do i = 1, size(MEF90SetType)
    !         write(*,*) 'MEF90SetLabelName(i) ', MEF90SetLabelName(i), MEF90SetPrefix(i)
    !     end do
    ! end block test2

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
