Program  TestSection
#include <petsc/finclude/petsc.h>
Use m_MEF90
Use petsc
Implicit NONE   
    
    PetscErrorCode                      :: ierr
    Type(MEF90Ctx_Type),target          :: MEF90Ctx
    Type(MEF90CtxGlobalOptions_Type)    :: MEF90GlobalOptions_default
    Type(tDM)                           :: dm,dmU
    PetscBool                           :: interpolate = PETSC_TRUE
    Character(len=MEF90MXSTRLEN)        :: name

    PetscInt                            :: dim
    Type(tPetscSection)                 :: section
    Type(tVec)                          :: U,U0,V,V0


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
    
    PetscCallA(PetscPrintf(MEF90Ctx%Comm,MEF90Ctx%geometryfile//'\n',ierr))
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
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-mef90dm_view",ierr))

    PetscCallA(DMGetDimension(dm,dim,ierr))

    name = "U"
    PetscCallA(MEF90CreateLocalVector(dm,MEF90GlobalOptions_default%elementFamily,MEF90GlobalOptions_default%elementOrder,dim,name,U,ierr))
    name = "U0"
    PetscCallA(MEF90CreateBoundaryLocalVector(dm,MEF90GlobalOptions_default%elementFamily,MEF90GlobalOptions_default%elementOrder,dim,name,U0,ierr))

    name = "V"
    PetscCallA(MEF90CreateCellVector(dm,dim,name,V,ierr))
    name = "V0"
    PetscCallA(MEF90CreateBoundaryCellVector(dm,dim,name,V0,ierr))

    PetscCallA(VecGetDM(U,dmU,ierr))
    PetscCallA(DMGetLocalSection(dmU,section,ierr))
    PetscCallA(PetscSectionViewFromOptions(section,PETSC_NULL_OPTIONS,"-sectionU_view",ierr))
    PetscCallA(VecViewFromOptions(U,PETSC_NULL_OPTIONS,"-U_view",ierr))
    PetscCallA(PetscSectionDestroy(section,ierr))
    PetscCallA(DMDestroy(dmU,ierr))

    PetscCallA(VecGetDM(U0,dmU,ierr))
    PetscCallA(DMGetLocalSection(dmU,section,ierr))
    PetscCallA(PetscSectionViewFromOptions(section,PETSC_NULL_OPTIONS,"-sectionU0_view",ierr))
    PetscCallA(VecViewFromOptions(U0,PETSC_NULL_OPTIONS,"-U0_view",ierr))
    PetscCallA(PetscSectionDestroy(section,ierr))
    PetscCallA(DMDestroy(dmU,ierr))

    PetscCallA(VecGetDM(V,dmU,ierr))
    PetscCallA(DMGetLocalSection(dmU,section,ierr))
    PetscCallA(PetscSectionViewFromOptions(section,PETSC_NULL_OPTIONS,"-sectionV_view",ierr))
    PetscCallA(VecViewFromOptions(V,PETSC_NULL_OPTIONS,"-V_view",ierr))
    PetscCallA(PetscSectionDestroy(section,ierr))
    PetscCallA(DMDestroy(dmU,ierr))

    PetscCallA(VecGetDM(V0,dmU,ierr))
    PetscCallA(DMGetLocalSection(dmU,section,ierr))
    PetscCallA(PetscSectionViewFromOptions(section,PETSC_NULL_OPTIONS,"-sectionV0_view",ierr))
    PetscCallA(VecViewFromOptions(V0,PETSC_NULL_OPTIONS,"-V0_view",ierr))
    PetscCallA(PetscSectionDestroy(section,ierr))
    PetscCallA(DMDestroy(dmU,ierr))

    PetscCallA(VecDestroy(U,ierr))
    PetscCallA(VecDestroy(U0,ierr))
    PetscCallA(VecDestroy(V,ierr))
    PetscCallA(VecDestroy(V0,ierr))
    PetscCallA(DMDestroy(dm,ierr))

    Call MEF90CtxDestroy(MEF90Ctx,ierr)   
    Call MEF90Finalize(ierr)
    Call PetscFinalize(ierr)
End Program  TestSection
 
       
