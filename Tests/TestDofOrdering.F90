Program  TestDofOrdering
#include <petsc/finclude/petsc.h>
Use m_MEF90
Use petsc
Implicit NONE   
    
    PetscErrorCode                                 :: ierr
    Type(MEF90Ctx_Type),target                     :: MEF90Ctx
    Type(MEF90CtxGlobalOptions_Type)               :: MEF90GlobalOptions_default
    Type(tDM)                                      :: dm,dmU
    Type(tPetscSection)                            :: sectionU
    PetscBool                                      :: interpolate = PETSC_TRUE
    Character(len=MEF90MXSTRLEN)                   :: name
    Type(MEF90CtxGlobalOptions_Type),pointer       :: MEF90GlobalOptions

    PetscInt                                       :: nVal,i,dim
    Type(tVec)                                     :: U
    PetscReal,Dimension(:),Pointer                 :: UArray

    PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
    PetscCallA(MEF90Initialize(ierr))

    MEF90GlobalOptions_default%verbose           = 1
    MEF90GlobalOptions_default%dryrun            = PETSC_FALSE
    MEF90GlobalOptions_default%timeMin           = 0.0_Kr
    MEF90GlobalOptions_default%timeMax           = 1.0_Kr
    MEF90GlobalOptions_default%timeNumStep       = 11
    MEF90GlobalOptions_default%timeInterpolation = MEF90TimeInterpolation_linear
    MEF90GlobalOptions_default%timeSkip          = 0
    MEF90GlobalOptions_default%timeNumCycle      = 1
    MEF90GlobalOptions_default%elementFamily     = MEF90ElementFamilyLagrange
    MEF90GlobalOptions_default%elementOrder      = 1

    PetscCallA(MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90GlobalOptions_default,ierr))
    PetscCallA(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr))    
    
    PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm,MEF90Ctx%geometryfile,PETSC_NULL_CHARACTER,interpolate,dm,ierr))
    PetscCallA(DMPlexDistributeSetDefault(dm,PETSC_FALSE,ierr))
    PetscCallA(DMSetFromOptions(dm,ierr))
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-dm_view",ierr))
    PetscCallA(DMGetDimension(dm,dim,ierr))
    
    distribute: Block 
        Type(tDM),target                    :: dmDist
        PetscInt                            :: ovlp = 0
        If (MEF90Ctx%NumProcs > 1) Then
            PetscCallA(DMPlexDistribute(dm,ovlp,PETSC_NULL_SF,dmDist,ierr))
            PetscCallA(DMDestroy(dm,ierr))
            dm = dmDist
            PetscCallA(DMSetUseNatural(dm,PETSC_TRUE,ierr))
        End If
    End Block distribute
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-mef90dm_view",ierr))

    ! Create nodal local Vec holding coordinates
    name = "U"
    PetscCallA(MEF90CreateLocalVector(dm,MEF90GlobalOptions%elementFamily,MEF90GlobalOptions%elementOrder,1_Ki,name,U,ierr))
    PetscCallA(VecGetDM(U,dmU,ierr))
    PetscCallA(DMGetLocalSection(dmU,sectionU,ierr))
    PetscCallA(PetscSectionViewFromOptions(sectionU,PETSC_NULL_OPTIONS,"-mef90sectionU_view",ierr))
    PetscCallA(VecSet(U,0.0_Kr,ierr))

    PetscCallA(VecGetLocalSize(U,nVal,ierr))
    Do i = 1, nVal
        PetscCallA(VecGetArrayF90(U,UArray,ierr))
        UArray    = 0.0_Kr
        UArray(i) = 1.0_Kr
        PetscCallA(VecRestoreArrayF90(U,UArray,ierr))
        PetscCallA(DMPlexVecGetClosure(dmU,sectionU,U,0_Ki,UArray,ierr))
        Write(*,*) i, UArray
        PetscCallA(DMPlexVecRestoreClosure(dmU,sectionU,U,0_Ki,UArray,ierr))
    End Do

    PetscCallA(VecDestroy(U,ierr))
    PetscCallA(DMDestroy(dm,ierr))

    Call MEF90CtxDestroy(MEF90Ctx,ierr)   
    Call MEF90Finalize(ierr)
    Call PetscFinalize(ierr)
End Program  TestDofOrdering
 
       
