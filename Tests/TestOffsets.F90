Program  TestOffsets
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
    Type(tVec)                          :: U


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
 
    PetscCallA(PetscInitialize(ierr))
    
    Call MEF90Initialize(PETSC_COMM_WORLD,ierr)
    Call MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90GlobalOptions_default,ierr)
    
    PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm,MEF90Ctx%geometryfile,PETSC_NULL_CHARACTER,interpolate,dm,ierr))
    PetscCallA(DMPlexDistributeSetDefault(dm,PETSC_FALSE,ierr))
    PetscCallA(DMSetFromOptions(dm,ierr))
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OBJECT,"-mef90dm_view",ierr))
    distribute: Block 
        Type(tDM),target                    :: dmDist
        PetscInt                            :: ovlp = 2
        If (MEF90Ctx%NumProcs > 1) Then
            PetscCallA(DMPlexDistribute(dm,ovlp,PETSC_NULL_SF,dmDist,ierr))
            PetscCallA(DMDestroy(dm,ierr))
            dm = dmDist
        End If
    End Block distribute
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OBJECT,"-mef90dm_view",ierr))
    PetscCallA(DMGetDimension(dm,dim,ierr))



    name = "U"
    PetscCallA(MEF90CreateLocalVector(dm,MEF90GlobalOptions_default%elementFamily,MEF90GlobalOptions_default%elementOrder,dim,name,U,ierr))

    PetscCallA(VecGetDM(U,dmU,ierr))
    PetscCallA(DMGetLocalSection(dmU,section,ierr))
    PetscCallA(PetscSectionViewFromOptions(section,PETSC_NULL_OBJECT,"-sectionU_view",ierr))
    PetscCallA(VecViewFromOptions(U,PETSC_NULL_OBJECT,"-U_view",ierr))

    PetscCallA(VecDestroy(U,ierr))

    testRemote: block
        Type(tPetscSection)                     :: locsection, gsection
        Type(tPetscSF)                          :: overlapSF, idSF, sf
        PetscInt                                :: pStart,pEnd,n,p
        Type(sPetscSFNode),dimension(:),Pointer :: remote
        PetscInt,dimension(:),Pointer           :: remoteOffsets


        PetscCallA(DMGetLocalSection(dm,locSection,ierr))
        PetscCall(DMGetPointSF(dm,overlapSF,ierr))
        PetscCallA(PetscSectionCreateGlobalSection(locSection,overlapSF,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE,gSection,ierr))

        PetscCall(PetscSectionGetChart(locSection,pStart,pEnd,ierr))
        n = pEnd-pStart
        Allocate(remote(n))
        Do p = 1,n
            remote(p)%rank  = MEF90Ctx%rank
            remote(p)%index = p-1
        End Do

        PetscCall(PetscSFCreate(MEF90Ctx%Comm,idSF,ierr))
        PetscCall(PetscSFSetFromOptions(idSF,ierr))
        PetscCall(PetscSFSetGraph(idSF,n,n,PETSC_NULL_INTEGER_ARRAY,PETSC_COPY_VALUES,remote,PETSC_COPY_VALUES,ierr))
        PetscCall(PetscSFSetUp(idSF,ierr))
        PetscCall(PetscSFCreateRemoteOffsets(idSF,locSection,gSection,remoteOffsets,ierr))
        PetscCallA(PetscSectionViewFromOptions(locsection,PETSC_NULL_OBJECT,"-locsection_view",ierr))
        PetscCallA(PetscSectionViewFromOptions(gsection,PETSC_NULL_OBJECT,"-gsection_view",ierr))
        if (associated(remoteOffsets)) then
            PetscCall(PetscSFCreateSectionSF(idSF,locSection,remoteOffsets,gSection,sf,ierr))    
            PetscCall(PetscSFDestroyRemoteOffsets(remoteOffsets,ierr))
        end if
    end block testRemote


    PetscCallA(DMDestroy(dm,ierr))

    Call MEF90CtxDestroy(MEF90Ctx,ierr)   
    Call MEF90Finalize(ierr)
    Call PetscFinalize(ierr)
End Program  TestOffsets
 
       
! mpirun -np 3 ./TestOffsets -geometry ../TestMeshes/SquareFaceSetCubit2CS.gen -result test.exo -mef90dm_view -log_view