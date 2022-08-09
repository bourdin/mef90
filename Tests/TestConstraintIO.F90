Program  TestConstraintIO
#include <petsc/finclude/petsc.h>
Use m_MEF90
Use petsc
Implicit NONE   
    
    PetscErrorCode                                      :: ierr
    Type(MEF90Ctx_Type),target                          :: MEF90Ctx
    Type(MEF90CtxGlobalOptions_Type)                    :: MEF90GlobalOptions_default
    Type(MEF90CtxGlobalOptions_Type),Pointer            :: MEF90GlobalOptions
    Type(tDM),target                                    :: dm,dmU,dmU0,dmS
    PetscBool                                           :: interpolate = PETSC_TRUE

    PetscInt                                            :: numCellComponents, numNodalVar = 2, numCellVar = 3, numGVar = 0, numStep = 3
    Character(len=MEF90MXSTRLEN),Dimension(:),Pointer   :: nodalVarName, cellVarName, gVarName
    Character(len=MEF90MXSTRLEN)                        :: name
    PetscInt                                            :: nroots, nleaves
    type(tIS)                                           :: cellIS,csIS
    PetscInt,Dimension(:),pointer                       :: ilocal,cellID,csID
    PetscInt                                            :: dim,pStart,pEnd,size,bs,cell,c,set,numCells,numCS
    Type(tPetscSection)                                 :: sectionU, sectionU0, sectionS, coordSection
    Type(tPetscSF)                                      :: naturalPointSF, natSFU, natSFU0, natSFS, lcSF, clSF, lioSF, iolSF, lioBSF, iolBSF, lioSSF, iolSSF
    Type(tVec)                                          :: locCoord, locVec, locVecB, ioVec, S, locS, ioS, ioVecRead, ioSRead
    Type(PetscSFNode),dimension(:),Pointer              :: iremote
    Type(tPetscViewer)                                  :: viewer
    PetscReal                                           :: time = 1.234_Kr,norm
    PetscScalar,dimension(:),pointer                    :: cval,xyz

    PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
    Call MEF90Initialize(ierr)

    MEF90GlobalOptions_default%verbose           = 0
    MEF90GlobalOptions_default%dryrun            = PETSC_FALSE
    MEF90GlobalOptions_default%timeMin           = 0.0_Kr
    MEF90GlobalOptions_default%timeMax           = 1.0_Kr
    MEF90GlobalOptions_default%timeNumStep       = 11
    MEF90GlobalOptions_default%elementFamily     = MEF90ElementFamilyLagrange
    MEF90GlobalOptions_default%elementOrder      = 1
 
    Call MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90GlobalOptions_default,ierr)
    PetscCallA(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr))

    Allocate(nodalVarName(numNodalVar))
    Allocate(cellVarName(numCellVar))
    Allocate(gVarName(numGVar))
    nodalVarName(1:2) = ["U_X","U_Y"]
    cellVarName(1:3) = ["Sigma_11","Sigma_22","Sigma_33"]
    
    ! Create DM from file
    PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm,MEF90Ctx%geometryfile,PETSC_NULL_CHARACTER,interpolate,dm,ierr))
    PetscCallA(DMPlexDistributeSetDefault(dm,PETSC_FALSE,ierr))
    PetscCallA(DMSetFromOptions(dm,ierr))
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-dm_view",ierr))

    ! Open exodus file + write geometry + format the file
    call MEF90CtxOpenEXO(MEF90Ctx,viewer,ierr)
    call MEF90EXODMView(dm,viewer,ierr)
    call MEF90EXOFormat(viewer,gVarName,cellVarName,nodalVarName,numStep,ierr)

    DeAllocate(nodalVarName)
    DeAllocate(cellVarName)
    DeAllocate(gVarName)

    ! Distribute DM
    distribute: Block 
        Type(tDM),target                    :: dmDist
        If (MEF90Ctx%NumProcs > 1) Then
            PetscCallA(DMPlexDistribute(dm,0,naturalPointSF,dmDist,ierr))
            PetscCallA(DMPlexSetMigrationSF(dmDist,naturalPointSF, ierr))
            PetscCallA(PetscSFDestroy(naturalPointSF,ierr))
            PetscCallA(DMDestroy(dm,ierr))
            dm = dmDist
        End If
    End Block distribute
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-dm_view",ierr))

    PetscCallA(DMGetDimension(dm,dim,ierr))
    PetscCallA(DMPlexGetChart(dm,pStart,pEnd,ierr))

    PetscCallA(DMSetUseNatural(dm,PETSC_TRUE,ierr))

    ! Create local Vec holding coordinates
    name = "U"
    PetscCallA(MEF90VecCreate(dm,MEF90GlobalOptions%elementFamily,MEF90GlobalOptions%elementOrder,dim,name,locVec,ierr))
    PetscCallA(VecGetDM(locVec,dmU,ierr))
    PetscCallA(DMGetLocalSection(dmU,sectionU,ierr))
    PetscCallA(DMSetUseNatural(dmU,PETSC_TRUE,ierr))

    ! Create local Vec holding constraints
    name = "U0"
    PetscCallA(MEF90VecCreate(dm,MEF90GlobalOptions%elementFamily,MEF90GlobalOptions%elementOrder,dim,name,locVecB,ierr))
    PetscCallA(VecGetDM(locVecB,dmU0,ierr))
    PetscCallA(DMGetLocalSection(dmU0,sectionU0,ierr))
    PetscCallA(DMSetUseNatural(dmU0,PETSC_TRUE,ierr))

    ! TODO: create cell Vec through MEF90VecCreate
    PetscCallA(DMClone(dm,dmS,ierr))
    PetscCallA(PetscSectionCreate(MEF90Ctx%Comm,sectionS,ierr))
    PetscCallA(PetscObjectSetName(SectionS,"Section for Sigma",ierr))
    PetscCallA(PetscSectionSetChart(sectionS,pStart,pEnd,ierr))
    numCellComponents = 3
    PetscCallA(MEF90CellSectionCreate(dm,numCellComponents,sectionS,ierr))
    PetscCallA(DMSetLocalSection(dmS,sectionS,ierr))
    PetscCallA(DMSetUseNatural(dmS,PETSC_TRUE,ierr))

    ! View all sections
    PetscCallA(PetscSectionViewFromOptions(SectionU,PETSC_NULL_OPTIONS,"-sectionU_view",ierr))
    PetscCallA(PetscSectionViewFromOptions(SectionU0,PETSC_NULL_OPTIONS,"-sectionU0_view",ierr))
    PetscCallA(PetscSectionViewFromOptions(SectionS,PETSC_NULL_OPTIONS,"-sectionS_view",ierr))

    ! Create GlobalToNatural SF
    if (MEF90Ctx%NumProcs > 1) then
        PetscCallA(DMPlexGetMigrationSF(dm, naturalPointSF, ierr))
        PetscCallA(DMPlexSetMigrationSF(dmU, naturalPointSF, ierr))
        PetscCallA(DMPlexSetMigrationSF(dmU0, naturalPointSF, ierr))
        PetscCallA(DMPlexSetMigrationSF(dmS, naturalPointSF, ierr))
        PetscCallA(DMPlexCreateGlobalToNaturalSF(dmU, PETSC_NULL_SECTION, naturalPointSF, natSFU, ierr))
        PetscCallA(DMSetNaturalSF(dmU, natSFU, ierr))
        PetscCallA(PetscSFDestroy(natSFU, ierr))
        PetscCallA(DMPlexCreateGlobalToNaturalSF(dmU0, PETSC_NULL_SECTION, naturalPointSF, natSFU0, ierr))
        PetscCallA(DMSetNaturalSF(dmU0, natSFU0, ierr))
        PetscCallA(PetscSFDestroy(natSFU0, ierr))
        PetscCallA(DMPlexCreateGlobalToNaturalSF(dmS, PETSC_NULL_SECTION, naturalPointSF, natSFS, ierr))
        PetscCallA(DMSetNaturalSF(dmS, natSFS, ierr))
        PetscCallA(PetscSFDestroy(natSFS, ierr))
    end if

    ! Create all SFs for Vec reordering and copying constraint values
    PetscCallA(MEF90LocalToConstraintSFCreate(MEF90Ctx,dmU,dmU0,lcSF,clSF,ierr))
    PetscCallA(MEF90LocalToIOSFCreate(MEF90Ctx,dmU,lioSF,ierr))
    PetscCallA(MEF90IOToLocalSFCreate(MEF90Ctx,dmU,iolSF,ierr))
    PetscCallA(MEF90LocalToIOSFCreate(MEF90Ctx,dmU0,lioBSF,ierr))
    PetscCallA(MEF90IOToLocalSFCreate(MEF90Ctx,dmU0,iolBSF,ierr))
    PetscCallA(MEF90LocalToIOSFCreate(MEF90Ctx,dmS,lioSSF,ierr))
    PetscCallA(MEF90IOToLocalSFCreate(MEF90Ctx,dmS,iolSSF,ierr))

    ! Fill locVec holding coordinates and copy constraint values = 10 from locVecB
    PetscCallA(DMGetCoordinatesLocal(dmU,locCoord,ierr))
    PetscCallA(VecSet(locVecB,10.0_kr,ierr))
    PetscCallA(VecCopy(locCoord,locVec,ierr))
    PetscCallA(MEF90VecReorder(locVecB,locVec,clSF,ierr))

    ! Create IO Vec for writing (ioVec) and reading (ioVecRead)
    PetscCallA(PetscSFGetGraph(lioSF,nroots,nleaves,ilocal,iremote,ierr))
    PetscCallA(VecCreateMPI(MEF90Ctx%Comm,nleaves,PETSC_DETERMINE,iovec,ierr))
    PetscCallA(VecGetBlockSize(locVec,bs,ierr))
    PetscCallA(VecSetBlockSize(iovec,bs,ierr))
    PetscCallA(VecSetDM(iovec,dmU,ierr))
    PetscCallA(PetscObjectSetName(iovec,"U",ierr))
    PetscCallA(VecCreateMPI(MEF90Ctx%Comm,nleaves,PETSC_DETERMINE,ioVecRead,ierr))
    PetscCallA(VecSetBlockSize(ioVecRead,bs,ierr))
    PetscCallA(VecSetDM(ioVecRead,dmU,ierr))
    PetscCallA(PetscObjectSetName(ioVecRead,"U",ierr))

    ! Reorder locVec into ioVec and write ioVec
    PetscCallA(MEF90VecReorder(locVec,ioVec,lioSF,ierr))
    PetscCallA(DMSetOutputSequenceNumber(dmU,0,time,ierr))
    PetscCallA(VecViewFromOptions(ioVec,PETSC_NULL_OPTIONS,"-iovec_view",ierr))
    PetscCallA(MEF90EXOVecView(ioVec,viewer,ierr))

    ! Create Vec with 3 components on each cell: (i) rank, (ii) x geometric center,
    ! (iii) y geometric center
    PetscCallA(DMGetGlobalVector(dmS,S,ierr))
    PetscCallA(DMGetLocalVector(dmS,locS,ierr))
    PetscCallA(DMGetCoordinateSection(dmS, coordSection,ierr))
    PetscCallA(DMGetLabelIdIS(dmS, "Cell Sets", csIS,ierr))
    PetscCallA(DMGetLabelSize(dmS, "Cell Sets",numCS,ierr))
    PetscCallA(ISGetIndicesF90(csIS, csID,ierr))
    do set = 1, numCS
        PetscCallA(DMGetStratumIS(dmS, "Cell Sets", csID(set), cellIS,ierr))
        PetscCallA(ISGetIndicesF90(cellIS, cellID,ierr))
        PetscCallA(ISGetSize(cellIS, numCells,ierr))
        do cell = 1,numCells
            PetscCallA(DMPlexVecGetClosure(dmS, PETSC_NULL_SECTION, S, cellID(cell), cval,ierr))
            PetscCallA(DMPlexVecGetClosure(dmS, coordSection, locCoord, cellID(cell), xyz,ierr))
            cval(1) = MEF90Ctx%rank
            cval(2) = 0.0_Kr
            do c = 1, size(xyz),dim-1
                cval(2) = cval(2) + xyz(c)
            end do
            cval(2) = cval(2) * dim / size(xyz)
            cval(3) = 0.0_Kr
            do c = 1, size(xyz),dim
                cval(3) = cval(3) + xyz(c)
            end do
            cval(3) = cval(3) * dim / size(xyz)
            PetscCallA(DMPlexVecSetClosure(dmS, PETSC_NULL_SECTION, S, cellID(cell), cval, INSERT_ALL_VALUES,ierr))
            PetscCallA(DMPlexVecRestoreClosure(dmS, PETSC_NULL_SECTION, S, cellID(cell), cval,ierr))
            PetscCallA(DMPlexVecRestoreClosure(dmS, coordSection, locCoord, cellID(cell), xyz,ierr))
        end do
        PetscCallA(ISRestoreIndicesF90(cellIS, cellID,ierr))
        PetscCallA(ISDestroy(cellIS,ierr))
    end do
    PetscCallA(ISRestoreIndicesF90(csIS, csID,ierr))
    PetscCallA(ISDestroy(csIS,ierr))
    PetscCallA(DMGlobalToLocalBegin(dmS, S, INSERT_VALUES, locS,ierr))
    PetscCallA(DMGlobalToLocalEnd(dmS, S, INSERT_VALUES, locS,ierr))
    PetscCallA(DMRestoreGlobalVector(dmS,S,ierr))
    PetscCallA(PetscSFGetGraph(lioSSF,nroots,nleaves,ilocal,iremote,ierr))

    ! Create IO Vec for writing (ioS) and reading (ioSRead)
    PetscCallA(VecCreateMPI(MEF90Ctx%Comm,nleaves,PETSC_DETERMINE,ioS,ierr))
    PetscCallA(VecGetBlockSize(locS,bs,ierr))
    PetscCallA(VecSetBlockSize(ioS,bs,ierr))
    PetscCallA(VecSetDM(ioS,dmS,ierr))
    PetscCallA(PetscObjectSetName(ioS,"Sigma",ierr))
    PetscCallA(VecCreateMPI(MEF90Ctx%Comm,nleaves,PETSC_DETERMINE,ioSRead,ierr))
    PetscCallA(VecSetBlockSize(ioSRead,bs,ierr))
    PetscCallA(VecSetDM(ioSRead,dmS,ierr))
    PetscCallA(PetscObjectSetName(ioSRead,"Sigma",ierr))

    ! Reorder and write ioS
    PetscCallA(DMSetOutputSequenceNumber(dmS,0,time,ierr))
    PetscCallA(MEF90VecReorder(locS,ioS,lioSSF,ierr))
    PetscCallA(VecViewFromOptions(ioS,PETSC_NULL_OPTIONS,"-ios_view",ierr))
    PetscCallA(MEF90EXOVecView(ioS,viewer,ierr))

    ! Test read ioVecRead and ioSRead
    PetscCallA(MEF90EXOVecLoad(ioVecRead,viewer,ierr))
    PetscCallA(VecViewFromOptions(ioVecRead,PETSC_NULL_OPTIONS,"-iovec_view",ierr))
    PetscCallA(MEF90EXOVecLoad(ioSRead,viewer,ierr))
    PetscCallA(VecViewFromOptions(ioSRead,PETSC_NULL_OPTIONS,"-ios_view",ierr))

    ! Cleanup Vec
    PetscCallA(VecDestroy(locVec,ierr))
    PetscCallA(VecDestroy(locVecB,ierr))
    PetscCallA(DMRestoreLocalVector(dmS,locS,ierr))
    PetscCallA(VecDestroy(ioVec,ierr))
    PetscCallA(VecDestroy(ioVecRead,ierr))
    PetscCallA(VecDestroy(ioS,ierr))
    PetscCallA(VecDestroy(ioSRead,ierr))

    ! Cleanup SF
    PetscCallA(PetscSFDestroy(lcSF,ierr))
    PetscCallA(PetscSFDestroy(clSF,ierr))
    PetscCallA(PetscSFDestroy(lioSF,ierr))
    PetscCallA(PetscSFDestroy(iolSF,ierr))
    PetscCallA(PetscSFDestroy(lioBSF,ierr))
    PetscCallA(PetscSFDestroy(iolBSF,ierr))
    PetscCallA(PetscSFDestroy(lioSSF,ierr))
    PetscCallA(PetscSFDestroy(iolSSF,ierr))

    ! Cleanup Sections
    PetscCallA(PetscSectionDestroy(sectionU,ierr))
    PetscCallA(PetscSectionDestroy(sectionU0,ierr))
    PetscCallA(PetscSectionDestroy(sectionS,ierr))

    ! Cleanup DMs
    PetscCallA(DMDestroy(dmU,ierr))
    PetscCallA(DMDestroy(dmU0,ierr))
    PetscCallA(DMDestroy(dmS,ierr))
    PetscCallA(DMDestroy(dm,ierr))

    ! Exit nicely
    call MEF90CtxCloseEXO(viewer,ierr)
    Call MEF90CtxDestroy(MEF90Ctx,ierr)   
    Call MEF90Finalize(ierr)
    Call PetscFinalize(ierr)
End Program  TestConstraintIO