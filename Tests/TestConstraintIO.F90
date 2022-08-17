Program  TestConstraintIO
#include <petsc/finclude/petsc.h>
Use m_MEF90
Use petsc
Implicit NONE   
    
    PetscErrorCode                                      :: ierr
    Type(MEF90Ctx_Type),target                          :: MEF90Ctx
    Type(MEF90CtxGlobalOptions_Type)                    :: MEF90GlobalOptions_default
    Type(MEF90CtxGlobalOptions_Type),Pointer            :: MEF90GlobalOptions
    Type(tDM),target                                    :: dm,dmU,dmU0,dmSigma
    PetscBool                                           :: interpolate = PETSC_TRUE

    PetscInt                                            :: numNodalVar = 2, numCellVar = 3, numGVar = 0, numStep = 3
    Character(len=MEF90MXSTRLEN),Dimension(:),Pointer   :: nodalVarName, cellVarName, gVarName
    Character(len=MEF90MXSTRLEN)                        :: name
    type(tIS)                                           :: cellIS,csIS
    PetscInt,Dimension(:),pointer                       :: cellID,csID
    PetscInt                                            :: dim,pStart,pEnd,size,bs,cell,c,set,numCells,numCS
    Type(tPetscSection)                                 :: sectionU, sectionU0, sectionSigma, coordSection
    Type(tPetscSF)                                      :: naturalPointSF, natSFU, natSFU0, natSFS, lcSF, clSF, lioSF, iolSF, lioBSF, iolBSF, lioSSF, iolSSF
    Type(tVec)                                          :: locCoord, locVecU, locVecU0, ioVec, globVecSigma, locVecSigma, ioS, ioVecRead, ioSRead
    Type(tPetscViewer)                                  :: viewer
    PetscScalar,dimension(:),pointer                    :: cval,xyz
    PetscInt                                            :: orderSigma = 0,numCompSigma = 3
    PetscInt                                            :: step = 0

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
    nodalVarName = ["U_X","U_Y"]
    cellVarName  = ["Sigma_11","Sigma_22","Sigma_33"]
    
    ! Create DM from file
    PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm,MEF90Ctx%geometryfile,PETSC_NULL_CHARACTER,interpolate,dm,ierr))
    PetscCallA(DMPlexDistributeSetDefault(dm,PETSC_FALSE,ierr))
    PetscCallA(DMSetUseNatural(dm,PETSC_TRUE,ierr))
    PetscCallA(DMSetFromOptions(dm,ierr))
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-mef90dm_view",ierr))

    ! Open exodus file + write geometry + format the file
    call MEF90CtxOpenEXO(MEF90Ctx,viewer,ierr)
    call MEF90EXODMView(dm,viewer,MEF90GlobalOptions%elementOrder,ierr)
    call MEF90EXOFormat(viewer,gVarName,cellVarName,nodalVarName,numStep,ierr)

    DeAllocate(nodalVarName)
    DeAllocate(cellVarName)
    DeAllocate(gVarName)

    ! Distribute DM
    distribute: Block 
        Type(tDM),target                    :: dmDist
        PetscInt                            :: ovlp = 0
        If (MEF90Ctx%NumProcs > 1) Then
            PetscCallA(DMPlexDistribute(dm,ovlp,naturalPointSF,dmDist,ierr))
            PetscCallA(DMPlexSetMigrationSF(dmDist,naturalPointSF, ierr))
            PetscCallA(PetscSFDestroy(naturalPointSF,ierr))
            PetscCallA(DMDestroy(dm,ierr))
            dm = dmDist
        End If
    End Block distribute
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-mef90dm_view",ierr))

    PetscCallA(DMGetDimension(dm,dim,ierr))
    PetscCallA(DMPlexGetChart(dm,pStart,pEnd,ierr))

    PetscCallA(DMSetUseNatural(dm,PETSC_TRUE,ierr))

    ! Create nodal local Vec holding coordinates
    name = "U"
    PetscCallA(MEF90VecCreate(dm,MEF90GlobalOptions%elementFamily,MEF90GlobalOptions%elementOrder,dim,name,locVecU,ierr))
    PetscCallA(VecGetDM(locVecU,dmU,ierr))
    PetscCallA(DMGetLocalSection(dmU,sectionU,ierr))
    PetscCallA(DMSetUseNatural(dmU,PETSC_TRUE,ierr))

    ! Create nodal local Vec holding constraints
    name = "U0"
    PetscCallA(MEF90VecCreate(dm,MEF90GlobalOptions%elementFamily,MEF90GlobalOptions%elementOrder,dim,name,locVecU0,ierr))
    PetscCallA(VecGetDM(locVecU0,dmU0,ierr))
    PetscCallA(DMGetLocalSection(dmU0,sectionU0,ierr))
    PetscCallA(DMSetUseNatural(dmU0,PETSC_TRUE,ierr))

    ! create cell Vec holding sigma
    name = "Sigma"
    PetscCallA(MEF90VecCreate(dm,MEF90GlobalOptions%elementFamily,orderSigma,numCompSigma,name,locVecSigma,ierr))
    PetscCallA(VecGetDM(locVecSigma,dmSigma,ierr))
    PetscCallA(DMGetLocalSection(dmSigma,sectionSigma,ierr))

    PetscCallA(DMGetGlobalVector(dmSigma,globVecSigma,ierr))
    PetscCallA(DMGetLocalVector(dmSigma,locVecSigma,ierr))
    PetscCallA(DMGetCoordinatesLocal(dmSigma,locCoord,ierr))
    PetscCallA(DMGetCoordinateSection(dmSigma, coordSection,ierr))
    PetscCallA(DMGetLabelIdIS(dmSigma, "Cell Sets", csIS,ierr))
    PetscCallA(DMGetLabelSize(dmSigma, "Cell Sets",numCS,ierr))
    PetscCallA(ISGetIndicesF90(csIS, csID,ierr))
    do set = 1, numCS
        PetscCallA(DMGetStratumIS(dmSigma, "Cell Sets", csID(set), cellIS,ierr))
        PetscCallA(ISGetIndicesF90(cellIS, cellID,ierr))
        PetscCallA(ISGetSize(cellIS, numCells,ierr))
        do cell = 1,numCells
            PetscCallA(DMPlexVecGetClosure(dmSigma, PETSC_NULL_SECTION, globVecSigma, cellID(cell), cval,ierr))
            PetscCallA(DMPlexVecGetClosure(dmSigma, coordSection, locCoord, cellID(cell), xyz,ierr))
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
            PetscCallA(DMPlexVecSetClosure(dmSigma, PETSC_NULL_SECTION, globVecSigma, cellID(cell), cval, INSERT_ALL_VALUES,ierr))
            PetscCallA(DMPlexVecRestoreClosure(dmSigma, PETSC_NULL_SECTION, globVecSigma, cellID(cell), cval,ierr))
            PetscCallA(DMPlexVecRestoreClosure(dmSigma, coordSection, locCoord, cellID(cell), xyz,ierr))
        end do
        PetscCallA(ISRestoreIndicesF90(cellIS, cellID,ierr))
        PetscCallA(ISDestroy(cellIS,ierr))
    end do
    PetscCallA(ISRestoreIndicesF90(csIS, csID,ierr))
    PetscCallA(ISDestroy(csIS,ierr))
    PetscCallA(DMGlobalToLocalBegin(dmSigma, globVecSigma, INSERT_VALUES, locVecSigma,ierr))
    PetscCallA(DMGlobalToLocalEnd(dmSigma, globVecSigma, INSERT_VALUES, locVecSigma,ierr))
    PetscCallA(DMRestoreGlobalVector(dmSigma,globVecSigma,ierr))

    ! View all sections
    PetscCallA(PetscSectionViewFromOptions(SectionU,PETSC_NULL_OPTIONS,"-sectionU_view",ierr))
    PetscCallA(PetscSectionViewFromOptions(SectionU0,PETSC_NULL_OPTIONS,"-sectionU0_view",ierr))
    PetscCallA(PetscSectionViewFromOptions(SectionSigma,PETSC_NULL_OPTIONS,"-sectionSigma_view",ierr))

    ! Create GlobalToNatural SFs
    if (MEF90Ctx%NumProcs > 1) then
        PetscCallA(DMPlexGetMigrationSF(dm, naturalPointSF, ierr))
        PetscCallA(DMPlexSetMigrationSF(dmU, naturalPointSF, ierr))
        PetscCallA(DMPlexSetMigrationSF(dmU0, naturalPointSF, ierr))
        PetscCallA(DMPlexSetMigrationSF(dmSigma, naturalPointSF, ierr))
        PetscCallA(DMPlexCreateGlobalToNaturalSF(dmU, PETSC_NULL_SECTION, naturalPointSF, natSFU, ierr))
        PetscCallA(DMSetNaturalSF(dmU, natSFU, ierr))
        PetscCallA(PetscSFDestroy(natSFU, ierr))
        PetscCallA(DMPlexCreateGlobalToNaturalSF(dmU0, PETSC_NULL_SECTION, naturalPointSF, natSFU0, ierr))
        PetscCallA(DMSetNaturalSF(dmU0, natSFU0, ierr))
        PetscCallA(PetscSFDestroy(natSFU0, ierr))
        PetscCallA(DMPlexCreateGlobalToNaturalSF(dmSigma, PETSC_NULL_SECTION, naturalPointSF, natSFS, ierr))
        PetscCallA(DMSetNaturalSF(dmSigma, natSFS, ierr))
        PetscCallA(PetscSFDestroy(natSFS, ierr))
    end if

    ! Create SFs for copying from/into IO coordinates Vec
    PetscCallA(MEF90IOSFCreate(MEF90Ctx,locVecU,lioSF,iolSF,ierr))
    ! Create SFs for copying from/into IO Constraint Vec
    PetscCallA(MEF90IOSFCreate(MEF90Ctx,locVecU0,lioBSF,iolBSF,ierr))
    ! Create SFs for copying constrained dofs from/into Constraint Vec
    PetscCallA(MEF90ConstraintSFCreate(MEF90Ctx,locVecU,locVecU0,lcSF,clSF,ierr))
    ! Create SFs for copying from/into IO cell Vec
    PetscCallA(MEF90IOSFCreate(MEF90Ctx,locVecSigma,lioSSF,iolSSF,ierr))

    ! Fill locVecU holding coordinates and copy constraint values = 10 from locVecU0
    PetscCallA(VecSet(locVecU0,10.0_kr,ierr))
    ! locCoord is obtained from dmSigma but all DMs have the same coordinates
    PetscCallA(VecCopy(locCoord,locVecU,ierr))
    PetscCallA(MEF90VecCopySF(locVecU0,locVecU,clSF,ierr))

    ! Create IO Vec for writing (ioVec) and reading (ioVecRead)
    PetscCallA(VecGetBlockSize(locVecU,bs,ierr))
    PetscCallA(MEF90VecCreateIO(MEF90Ctx,ioVec,bs,lioSF,ierr))
    PetscCallA(PetscObjectSetName(ioVec,"U",ierr))
    PetscCallA(MEF90VecCreateIO(MEF90Ctx,ioVecRead,bs,lioSF,ierr))
    PetscCallA(PetscObjectSetName(ioVecRead,"U",ierr))

    ! Reorder locVecU into ioVec and write ioVec
    PetscCallA(MEF90VecCopySF(locVecU,ioVec,lioSF,ierr))
    PetscCallA(VecViewFromOptions(ioVec,PETSC_NULL_OPTIONS,"-iovec_view",ierr))
    PetscCallA(MEF90EXOVecView(ioVec,viewer,0,ierr))

    ! Create Vec with 3 components on each cell: (i) rank, (ii) x geometric center,
    ! (iii) y geometric center
    PetscCallA(DMGetGlobalVector(dmSigma,globVecSigma,ierr))
    PetscCallA(DMCreateLocalVector(dmSigma,locVecSigma,ierr))
    PetscCallA(DMGetCoordinateSection(dmSigma, coordSection,ierr))
    PetscCallA(DMGetLabelIdIS(dmSigma, "Cell Sets", csIS,ierr))
    PetscCallA(DMGetLabelSize(dmSigma, "Cell Sets",numCS,ierr))
    PetscCallA(ISGetIndicesF90(csIS, csID,ierr))
    Do set = 1, numCS
        PetscCallA(DMGetStratumIS(dmSigma, "Cell Sets", csID(set), cellIS,ierr))
        If (cellIS /= PETSC_NULL_IS) Then
            PetscCallA(ISGetIndicesF90(cellIS, cellID,ierr))
            PetscCallA(ISGetSize(cellIS, numCells,ierr))
            Do cell = 1,numCells
                PetscCallA(DMPlexVecGetClosure(dmSigma, PETSC_NULL_SECTION, globVecSigma, cellID(cell), cval,ierr))
                PetscCallA(DMPlexVecGetClosure(dmSigma, coordSection, locCoord, cellID(cell), xyz,ierr))
                cval(1) = MEF90Ctx%rank
                cval(2) = 0.0_Kr
                Do c = 1, size(xyz),dim-1
                    cval(2) = cval(2) + xyz(c)
                End Do
                cval(2) = cval(2) * dim / size(xyz)
                cval(3) = 0.0_Kr
                Do c = 1, size(xyz),dim
                    cval(3) = cval(3) + xyz(c)
                End Do
                cval(3) = cval(3) * dim / size(xyz)
                PetscCallA(DMPlexVecSetClosure(dmSigma, PETSC_NULL_SECTION, globVecSigma, cellID(cell), cval, INSERT_ALL_VALUES,ierr))
                PetscCallA(DMPlexVecRestoreClosure(dmSigma, PETSC_NULL_SECTION, globVecSigma, cellID(cell), cval,ierr))
                PetscCallA(DMPlexVecRestoreClosure(dmSigma, coordSection, locCoord, cellID(cell), xyz,ierr))
            End Do
            PetscCallA(ISRestoreIndicesF90(cellIS, cellID,ierr))
        End If ! cellIS
        PetscCallA(ISDestroy(cellIS,ierr))
    End Do
    PetscCallA(ISRestoreIndicesF90(csIS, csID,ierr))
    PetscCallA(ISDestroy(csIS,ierr))
    PetscCallA(DMGlobalToLocalBegin(dmSigma, globVecSigma, INSERT_VALUES, locVecSigma,ierr))
    PetscCallA(DMGlobalToLocalEnd(dmSigma, globVecSigma, INSERT_VALUES, locVecSigma,ierr))
    PetscCallA(DMRestoreGlobalVector(dmSigma,globVecSigma,ierr))
    PetscCallA(PetscSFGetGraph(lioSSF,nroots,nleaves,ilocal,iremote,ierr))

    ! Create IO Vec for writing (ioS) and reading (ioSRead)
    PetscCallA(VecGetBlockSize(locVecSigma,bs,ierr))
    PetscCallA(MEF90VecCreateIO(MEF90Ctx,ioS,bs,lioSSF,ierr))
    PetscCallA(PetscObjectSetName(ioS,"Sigma",ierr))
    PetscCallA(MEF90VecCreateIO(MEF90Ctx,ioSRead,bs,lioSSF,ierr))
    PetscCallA(PetscObjectSetName(ioSRead,"Sigma",ierr))

    ! Reorder and write ioS
    PetscCallA(MEF90VecCopySF(locVecSigma,ioS,lioSSF,ierr))
    PetscCallA(VecViewFromOptions(ioS,PETSC_NULL_OPTIONS,"-ios_view",ierr))
    PetscCallA(MEF90EXOVecView(ioS,viewer,0,ierr))

    ! Test read ioVecRead and ioSRead
    PetscCallA(MEF90EXOVecLoad(ioVecRead,viewer,0,ierr))
    PetscCallA(VecViewFromOptions(ioVecRead,PETSC_NULL_OPTIONS,"-iovec_view",ierr))
    PetscCallA(MEF90EXOVecLoad(ioSRead,viewer,0,ierr))
    PetscCallA(VecViewFromOptions(ioSRead,PETSC_NULL_OPTIONS,"-ios_view",ierr))

    ! Cleanup Vec
    PetscCallA(VecDestroy(locVecU,ierr))
    PetscCallA(VecDestroy(locVecU0,ierr))
    PetscCallA(DMRestoreLocalVector(dmSigma,locVecSigma,ierr))
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
    PetscCallA(PetscSectionDestroy(sectionSigma,ierr))

    ! Cleanup DMs
    PetscCallA(DMDestroy(dmU,ierr))
    PetscCallA(DMDestroy(dmU0,ierr))
    PetscCallA(DMDestroy(dmSigma,ierr))
    PetscCallA(DMDestroy(dm,ierr))

    ! Exit nicely
    Call MEF90CtxCloseEXO(viewer,ierr)
    Call MEF90CtxDestroy(MEF90Ctx,ierr)   
    Call MEF90Finalize(ierr)
    Call PetscFinalize(ierr)
End Program  TestConstraintIO