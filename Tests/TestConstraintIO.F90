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

    PetscInt                                            :: numNodalVar = 2, numCellVar = 3, numGVar = 0
    Character(len=MEF90MXSTRLEN),Dimension(:),Pointer   :: nodalVarName, cellVarName, gVarName
    Character(len=MEF90MXSTRLEN)                        :: name
    type(tIS)                                           :: cellIS,csIS
    PetscInt,Dimension(:),Pointer                       :: cellID,csID
    PetscInt                                            :: dim,pStart,pEnd,size,bs,cell,c,set
    Type(tPetscSection)                                 :: sectionU, sectionU0, sectionSigma, coordSection
    Type(tPetscSF)                                      :: naturalPointSF,lcSF, clSF, lioSF, iolSF, lioBSF, iolBSF, lioSSF, iolSSF
    Type(tVec)                                          :: locCoord, locVecU, locVecU0, ioVec, locVecSigma, ioS, ioVecRead, ioSRead
    PetscReal,Dimension(:),Pointer                      :: cval,xyz
    PetscReal,Dimension(:),Pointer                      :: time

    PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
    PetscCallA(MEF90Initialize(ierr))

    MEF90GlobalOptions_default%verbose           = 1
    MEF90GlobalOptions_default%dryrun            = PETSC_FALSE
    MEF90GlobalOptions_default%timeMin           = 0.0_Kr
    MEF90GlobalOptions_default%timeMax           = 1.0_Kr
    MEF90GlobalOptions_default%timeNumStep       = 11
    MEF90GlobalOptions_default%elementFamily     = MEF90ElementFamilyLagrange
    MEF90GlobalOptions_default%elementOrder      = 1
 
    PetscCallA(MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90GlobalOptions_default,ierr))
    PetscCallA(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr))
    PetscCallA(MEF90CtxGetTime(MEF90Ctx,time,ierr))

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
    PetscCallA(MEF90CtxOpenEXO(MEF90Ctx,MEF90Ctx%resultViewer,ierr))
    PetscCallA(MEF90EXODMView(dm,MEF90Ctx%resultViewer,MEF90GlobalOptions%elementOrder,ierr))
    PetscCallA(MEF90EXOFormat(MEF90Ctx%resultViewer,gVarName,cellVarName,nodalVarName,time,ierr))

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
    PetscCallA(MEF90CreateLocalVector(dm,MEF90GlobalOptions%elementFamily,MEF90GlobalOptions%elementOrder,dim,name,locVecU,ierr))
    PetscCallA(VecGetDM(locVecU,dmU,ierr))
    PetscCallA(DMGetLocalSection(dmU,sectionU,ierr))

    ! Create nodal local Vec holding constraints
    name = "U0"
    PetscCallA(MEF90CreateBoundaryLocalVector(dm,MEF90GlobalOptions%elementFamily,MEF90GlobalOptions%elementOrder,dim,name,locVecU0,ierr))
    PetscCallA(VecGetDM(locVecU0,dmU0,ierr))
    PetscCallA(DMGetLocalSection(dmU0,sectionU0,ierr))

    ! create cell Vec holding sigma
    name = "Sigma"
    PetscCall(MEF90CreateCellVector(dm,3_Ki,name,locVecSigma,ierr))
    PetscCallA(VecGetDM(locVecSigma,dmSigma,ierr))
    PetscCallA(DMGetLocalSection(dmSigma,sectionSigma,ierr))

    ! View all sections
    PetscCallA(PetscSectionViewFromOptions(SectionU,PETSC_NULL_OPTIONS,"-sectionU_view",ierr))
    PetscCallA(PetscSectionViewFromOptions(SectionU0,PETSC_NULL_OPTIONS,"-sectionU0_view",ierr))
    PetscCallA(PetscSectionViewFromOptions(SectionSigma,PETSC_NULL_OPTIONS,"-sectionSigma_view",ierr))

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
    PetscCallA(DMGetCoordinatesLocal(dmSigma,locCoord,ierr))
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
    PetscCallA(MEF90EXOVecView(ioVec,MEF90Ctx%resultViewer,1_Ki,ierr))

    ! Create Vec with 3 components on each cell: (i) rank, (ii) x geometric center,
    ! (iii) y geometric center
    PetscCallA(DMGetCoordinateSection(dmSigma, coordSection,ierr))
    PetscCallA(DMGetLabelIdIS(dmSigma, "Cell Sets", csIS,ierr))
    PetscCallA(ISGetIndicesF90(csIS, csID,ierr))
    Do set = 1, size(csID)
        PetscCallA(DMGetStratumIS(dmSigma, "Cell Sets", csID(set), cellIS,ierr))
        If (cellIS /= PETSC_NULL_IS) Then
            PetscCallA(ISGetIndicesF90(cellIS, cellID,ierr))
            Do cell = 1,size(cellID)
                PetscCallA(DMPlexVecGetClosure(dmSigma, PETSC_NULL_SECTION, locVecSigma, cellID(cell), cval,ierr))
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
                PetscCallA(DMPlexVecRestoreClosure(dmSigma, coordSection, locCoord, cellID(cell), xyz,ierr))
                PetscCallA(DMPlexVecSetClosure(dmSigma, PETSC_NULL_SECTION, locVecSigma, cellID(cell), cval, INSERT_VALUES, ierr))
                PetscCallA(DMPlexVecRestoreClosure(dmSigma, PETSC_NULL_SECTION, locVecSigma, cellID(cell), cval,ierr))
            End Do
            PetscCallA(ISRestoreIndicesF90(cellIS, cellID,ierr))
        End If ! cellIS
        PetscCallA(ISDestroy(cellIS,ierr))
    End Do
    PetscCallA(ISRestoreIndicesF90(csIS, csID,ierr))
    PetscCallA(ISDestroy(csIS,ierr))

    ! Create IO Vec for writing (ioS) and reading (ioSRead)
    PetscCallA(VecGetBlockSize(locVecSigma,bs,ierr))
    PetscCallA(MEF90VecCreateIO(MEF90Ctx,ioS,bs,lioSSF,ierr))
    PetscCallA(PetscObjectSetName(ioS,"Sigma",ierr))
    PetscCallA(MEF90VecCreateIO(MEF90Ctx,ioSRead,bs,lioSSF,ierr))
    PetscCallA(PetscObjectSetName(ioSRead,"Sigma",ierr))

    ! Reorder and write ioS
    PetscCallA(MEF90VecCopySF(locVecSigma,ioS,lioSSF,ierr))
    PetscCallA(VecViewFromOptions(ioS,PETSC_NULL_OPTIONS,"-ios_view",ierr))
    PetscCallA(MEF90EXOVecView(ioS,MEF90Ctx%resultViewer,1_Ki,ierr))

    ! Test read ioVecRead and ioSRead
    PetscCallA(MEF90EXOVecLoad(ioVecRead,MEF90Ctx%resultViewer,1_Ki,ierr))
    PetscCallA(VecViewFromOptions(ioVecRead,PETSC_NULL_OPTIONS,"-iovec_view",ierr))
    PetscCallA(MEF90EXOVecLoad(ioSRead,MEF90Ctx%resultViewer,1_Ki,ierr))
    PetscCallA(VecViewFromOptions(ioSRead,PETSC_NULL_OPTIONS,"-ios_view",ierr))

    ! Cleanup Vec
    PetscCallA(VecDestroy(locVecU,ierr))
    PetscCallA(VecDestroy(locVecU0,ierr))
    PetscCallA(VecDestroy(locVecSigma,ierr))
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

    ! Cleanup DMs
    PetscCallA(DMDestroy(dm,ierr))
    ! Note that I would need to manually destroy these DM no matter what
    PetscCallA(DMDestroy(dmU,ierr))
    PetscCallA(DMDestroy(dmU0,ierr))
    PetscCallA(DMDestroy(dmSigma,ierr))
    PetscCallA(DMDestroy(dm,ierr))

    ! Exit nicely
    PetscCallA(MEF90CtxDestroy(MEF90Ctx,ierr))
    PetscCallA(MEF90Finalize(ierr))
    PetscCallA(PetscFinalize(ierr))
End Program  TestConstraintIO