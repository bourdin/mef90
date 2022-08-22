Program  TestVecSetBCFromOptions
#include <petsc/finclude/petsc.h>
Use petsc
Use m_MEF90
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
    
    Call MEF90Initialize(ierr)
    Call MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90GlobalOptions_default,ierr)
    PetscCallA(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr))
    
    PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm,MEF90Ctx%geometryfile,PETSC_NULL_CHARACTER,interpolate,dm,ierr))
    PetscCallA(DMPlexDistributeSetDefault(dm,PETSC_FALSE,ierr))
    PetscCallA(DMSetFromOptions(dm,ierr))
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-mef90dm_view",ierr))
    
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
    PetscCallA(MEF90VecSetBCFromOptions(V,scalingFactor,ierr))

    ViewSec: block
        Type(tPetscSection)     :: sectionV
        Type(tDM)               :: dmV

        PetscCallA(VecGetDM(V,dmV,ierr))
        PetscCallA(DMGetLocalSection(dmV,sectionV,ierr))
        PetscCallA(PetscSectionViewFromOptions(sectionV,PETSC_NULL_OPTIONS,"-mef90section_view",ierr))
        PetscCallA(VecViewFromOptions(V,PETSC_NULL_OPTIONS,"-mef90vec_view",ierr))
    end block ViewSec


    PetscCallA(VecDestroy(v,ierr))
    PetscCallA(DMDestroy(dm,ierr))
    Call MEF90CtxDestroy(MEF90Ctx,ierr)   
    Call MEF90Finalize(ierr)
    Call PetscFinalize(ierr)
End Program  TestVecSetBCFromOptions
