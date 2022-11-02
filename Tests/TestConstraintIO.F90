Module localFunctions
#include <petsc/finclude/petsc.h>
use petsc
use m_MEF90
implicit none
    
contains

#undef __FUNCT__
#define __FUNCT__ "project"

    subroutine project(v,s,ierr)
        Type(tVec),intent(IN)              :: v
        Type(tPetscSection),intent(IN)     :: s
        PetscErrorCode,intent(INOUT)       :: ierr

        PetscInt                           :: pStart,pEnd,p,numDof,i
        Type(tDM)                          :: dm
        Type(tPetscSection)                :: coordSection
        Type(tVec)                         :: coordVec
        PetscScalar,dimension(:),Pointer   :: coordArray,vArray
        PetscScalar,dimension(3)           :: xyz
        PetscInt                           :: dim,pOffset

        PetscCallA(PetscSectionGetChart(s,pStart,pEnd,ierr))
        PetscCallA(VecGetDM(v,dm,ierr))
        PetscCallA(DMGetCoordinateSection(dm,coordSection,ierr))
        PetscCallA(DMGetCoordinatesLocal(dm,coordVec,ierr))
        PetscCallA(DMGetDimension(dm,dim,ierr))
        PetscCallA(VecGetArrayF90(v,vArray,ierr))

        Do p = pStart,pEnd-1
            PetscCallA(PetscSectionGetDof(s,p,numDof,ierr))
            If (numDof > 0) Then
                !!! trick: the coordinate of a point is the average of the coordinates of the points in its closure
                PetscCallA(DMPlexVecGetClosure(dm,coordSection,coordVec,p,coordArray,ierr))
                Do i = 1,dim
                    xyz(i) = sum(coordArray(i:size(coordArray):dim)) * dim / size(coordArray)
                End Do
                PetscCallA(DMPlexVecRestoreClosure(dm,coordSection,coordVec,p,coordArray,ierr))

                PetscCallA(PetscSectionGetOffset(s,p,pOffset,ierr))
                Do i = 1,numDof
                    vArray(pOffset+i) = xyz(i)
                End Do
            End If
        End Do
        PetscCallA(VecRestoreArrayF90(v,vArray,ierr))
        !!! Of course, this does not use informations from the section, so it does over-write constrained values
    End subroutine project

End Module localFunctions

Program  TestConstraintIO
#include <petsc/finclude/petsc.h>
Use m_MEF90
Use petsc
Use localFunctions
Implicit NONE   
    
    PetscErrorCode                                      :: ierr
    Type(MEF90Ctx_Type),target                          :: MEF90Ctx
    Type(MEF90CtxGlobalOptions_Type)                    :: MEF90GlobalOptions_default
    Type(MEF90CtxGlobalOptions_Type),Pointer            :: MEF90GlobalOptions
    Type(tDM),target                                    :: dm,dmU,dmU0,dmSigma,dmSigma0
    PetscBool                                           :: interpolate = PETSC_TRUE

    PetscInt                                            :: numNodalVar = 2, numCellVar = 3, numGVar = 0, numSideVar = 3
    Character(len=MEF90MXSTRLEN),Dimension(:),Pointer   :: nodalVarName, cellVarName, gVarName, sideVarName
    Character(len=MEF90MXSTRLEN)                        :: name
    type(tIS)                                           :: cellIS,csIS,faceIS,ssIS
    PetscInt,Dimension(:),Pointer                       :: cellID,csID,faceID,ssID
    PetscInt                                            :: dim,pStart,pEnd,size,cell,c,set,face
    Type(tPetscSection)                                 :: sectionU, sectionU0, sectionSigma, sectionSigma0, coordSection,coordSection0
    Type(tPetscSF)                                      :: naturalPointSF,lcSF, clSF, lcSSF, clSSF, lioSF, iolSF, lioBSF, iolBSF, lioSSF, iolSSF, lioBSSF, iolBSSF
    Type(tVec)                                          :: locCoord, locVecU, locVecU0, locVecSigma, locVecSigma0
    PetscReal,Dimension(:),Pointer                      :: cval,xyz
    PetscReal,Dimension(:),Pointer                      :: time
    PetscInt                                            :: step = 1_Ki

    PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
    PetscCallA(MEF90Initialize(ierr))

    MEF90GlobalOptions_default%verbose           = 1
    MEF90GlobalOptions_default%dryrun            = PETSC_FALSE
    MEF90GlobalOptions_default%timeMin           = 0.0_Kr
    MEF90GlobalOptions_default%timeMax           = 1.0_Kr
    MEF90GlobalOptions_default%timeNumStep       = 1
    MEF90GlobalOptions_default%timeSkip          = 0
    MEF90GlobalOptions_default%timeNumCycle      = 1
    MEF90GlobalOptions_default%timeInterpolation = MEF90TimeInterpolation_linear
    MEF90GlobalOptions_default%elementFamily     = MEF90ElementFamilyLagrange
    MEF90GlobalOptions_default%elementOrder      = 1
 
    PetscCallA(MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90GlobalOptions_default,ierr))
    PetscCallA(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr))
    PetscCallA(MEF90CtxGetTime(MEF90Ctx,time,ierr))

    Allocate(nodalVarName(numNodalVar))
    Allocate(cellVarName(numCellVar))
    Allocate(gVarName(numGVar))
    Allocate(sideVarName(numSideVar))
    nodalVarName = ["U_X","U_Y"]
    cellVarName  = ["Sigma_11","Sigma_22","Sigma_33"]
    sideVarName = ["Sigma0_11","Sigma0_22","Sigma0_33"]
    
    ! Create DM from file
    PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm,MEF90Ctx%geometryfile,PETSC_NULL_CHARACTER,interpolate,dm,ierr))
    PetscCallA(DMPlexDistributeSetDefault(dm,PETSC_FALSE,ierr))
    PetscCallA(DMSetUseNatural(dm,PETSC_TRUE,ierr))
    PetscCallA(DMSetFromOptions(dm,ierr))
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-mef90dm_view",ierr))
    
    ! Open exodus file + write geometry + format the file
    PetscCallA(MEF90CtxOpenEXO(MEF90Ctx,MEF90Ctx%resultViewer,FILE_MODE_WRITE,ierr))
    PetscCallA(MEF90EXODMView(dm,MEF90Ctx%resultViewer,MEF90GlobalOptions%elementOrder,ierr))
    PetscCallA(MEF90EXOFormat(MEF90Ctx%resultViewer,gVarName,cellVarName,nodalVarName,sideVarName,time,ierr))

    DeAllocate(nodalVarName)
    DeAllocate(cellVarName)
    DeAllocate(gVarName)
    DeAllocate(sideVarName)

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

    ! create cell Vec holding sigma
    name = "Sigma0"
    PetscCall(MEF90CreateBoundaryCellVector(dm,3_Ki,name,locVecSigma0,ierr))
    PetscCallA(VecGetDM(locVecSigma0,dmSigma0,ierr))
    PetscCallA(DMGetLocalSection(dmSigma0,sectionSigma0,ierr))

    ! View all sections
    PetscCallA(PetscSectionViewFromOptions(SectionU,PETSC_NULL_OPTIONS,"-sectionU_view",ierr))
    PetscCallA(PetscSectionViewFromOptions(SectionU0,PETSC_NULL_OPTIONS,"-sectionU0_view",ierr))
    PetscCallA(PetscSectionViewFromOptions(SectionSigma,PETSC_NULL_OPTIONS,"-sectionSigma_view",ierr))
    PetscCallA(PetscSectionViewFromOptions(SectionSigma0,PETSC_NULL_OPTIONS,"-sectionSigma0_view",ierr))

    ! Create SFs for copying from/into IO coordinates Vec
    PetscCallA(MEF90IOSFCreate(MEF90Ctx,locVecU,lioSF,iolSF,ierr))
    ! Create SFs for copying from/into IO Constraint Vec
    PetscCallA(MEF90IOSFCreate(MEF90Ctx,locVecU0,lioBSF,iolBSF,ierr))
    ! Create SFs for copying constrained dofs from/into Constraint Vec
    PetscCallA(MEF90ConstraintSFCreate(MEF90Ctx,locVecU,locVecU0,lcSF,clSF,ierr))
    ! Create SFs for copying from/into IO cell Vec
    PetscCallA(MEF90IOSFCreate(MEF90Ctx,locVecSigma,lioSSF,iolSSF,ierr))
    ! Create SFs for copying from/into IO cell Vec
    PetscCallA(MEF90BoundaryIOSFCreate(MEF90Ctx,locVecSigma0,lioBSSF,iolBSSF,ierr))
    ! Create SFs for copying constrained dofs from/into Constraint Vec
    PetscCallA(MEF90ConstraintSFCreate(MEF90Ctx,locVecSigma,locVecSigma0,lcSSF,clSSF,ierr))

    ! Fill locVecU holding coordinates and copy constraint values = 10 from locVecU0
    PetscCallA(VecSet(locVecU0,10.0_kr,ierr))
    ! locCoord is obtained from dmSigma but all DMs have the same coordinates
    PetscCallA(DMGetCoordinatesLocal(dmSigma,locCoord,ierr))
    ! Fill locVecU with coordinates
    PetscCall(project(locVecU,sectionU,ierr))
    PetscCallA(MEF90VecCopySF(locVecU0,locVecU,clSF,ierr))

    ! Reorder locVecU into ioVec and write ioVec
    PetscCallA(VecViewFromOptions(locVecU,PETSC_NULL_OPTIONS,"-iovec_view",ierr))
    PetscCallA(MEF90EXOVecView(locVecU,lioSF,iolSF,MEF90Ctx%resultViewer,step,ierr))

    ! Create Vec with 3 components on each cell: (i) rank, (ii) x geometric center,
    ! (iii) y geometric center
    PetscCallA(DMGetCoordinateSection(dmSigma0, coordSection0,ierr))
    PetscCallA(DMGetLabelIdIS(dmSigma0, "Face Sets", ssIS,ierr))
    PetscCallA(ISGetIndicesF90(ssIS, ssID,ierr))
    Do set = 1, size(ssID)
        PetscCallA(DMGetStratumIS(dmSigma0, "Face Sets", ssID(set), faceIS,ierr))
        If (faceIS /= PETSC_NULL_IS) Then
            PetscCallA(ISGetIndicesF90(faceIS, faceID,ierr))
            Do face = 1,size(faceID)
                PetscCallA(DMPlexVecGetClosure(dmSigma0, PETSC_NULL_SECTION, locVecSigma0, faceID(face), cval,ierr))
                PetscCallA(DMPlexVecGetClosure(dmSigma0, coordSection0, locCoord, faceID(face), xyz,ierr))
                cval(1) = MEF90Ctx%rank
                cval(2) = 0.0_Kr
                Do c = 1, size(xyz),dim
                    cval(2) = cval(2) + xyz(c)
                End Do
                cval(2) = cval(2) * dim / size(xyz)
                cval(3) = 0.0_Kr
                Do c = 0, size(xyz),dim
                    cval(3) = cval(3) + xyz(c)
                End Do
                cval(3) = cval(3) * dim / size(xyz)
                !write(*,*) 'RANK = ',MEF90Ctx%rank,' set = ',set,' set ID = ',ssID(set),' face = ',face,' faceID = ',faceID(face),' x = ',cval(2),' y = ',cval(3)
                PetscCallA(DMPlexVecRestoreClosure(dmSigma0, coordSection0, locCoord, faceID(face), xyz,ierr))
                PetscCallA(DMPlexVecSetClosure(dmSigma0, PETSC_NULL_SECTION, locVecSigma0, faceID(face), cval, INSERT_VALUES, ierr))
                PetscCallA(DMPlexVecRestoreClosure(dmSigma0, PETSC_NULL_SECTION, locVecSigma0, faceID(face), cval,ierr))
            End Do
            PetscCallA(ISRestoreIndicesF90(faceIS, faceID,ierr))
        End If ! cellIS
        PetscCallA(ISDestroy(faceIS,ierr))
    End Do
    PetscCallA(ISRestoreIndicesF90(ssIS, ssID,ierr))
    PetscCallA(ISDestroy(ssIS,ierr))
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
                PetscCallA(DMPlexVecSetClosure(dmSigma, PETSC_NULL_SECTION, locVecSigma, cellID(cell), cval, INSERT_ALL_VALUES, ierr))
                PetscCallA(DMPlexVecRestoreClosure(dmSigma, PETSC_NULL_SECTION, locVecSigma, cellID(cell), cval,ierr))
            End Do
            PetscCallA(ISRestoreIndicesF90(cellIS, cellID,ierr))
        End If ! cellIS
        PetscCallA(ISDestroy(cellIS,ierr))
    End Do
    PetscCallA(ISRestoreIndicesF90(csIS, csID,ierr))
    PetscCallA(ISDestroy(csIS,ierr))

    ! Reorder and write ioS
    PetscCallA(VecViewFromOptions(locVecSigma,PETSC_NULL_OPTIONS,"-ios_view",ierr))
    PetscCallA(MEF90EXOVecView(locVecSigma,lioSSF,iolSSF,MEF90Ctx%resultViewer,step,ierr))

    ! Reorder and write ioS0
    PetscCallA(VecViewFromOptions(locVecSigma0,PETSC_NULL_OPTIONS,"-ios0_view",ierr))
    PetscCallA(MEF90EXOVecView(locVecSigma0,lioBSSF,iolBSSF,MEF90Ctx%resultViewer,step,ierr))

    ! Test read ioVecRead and ioSRead
    PetscCallA(VecSet(locVecU,1000.0_kr,ierr))
    PetscCallA(MEF90EXOVecLoad(locVecU,lioSF,iolSF,MEF90Ctx%resultViewer,step,ierr))
    PetscCallA(VecViewFromOptions(locVecU,PETSC_NULL_OPTIONS,"-iovec_view",ierr))
    PetscCallA(VecSet(locVecSigma,1000.0_kr,ierr))
    PetscCallA(MEF90EXOVecLoad(locVecSigma,lioSSF,iolSSF,MEF90Ctx%resultViewer,step,ierr))
    PetscCallA(VecViewFromOptions(locVecSigma,PETSC_NULL_OPTIONS,"-ios_view",ierr))
    PetscCallA(VecSet(locVecSigma0,1000.0_kr,ierr))
    PetscCallA(MEF90EXOVecLoad(locVecSigma0,lioBSSF,iolBSSF,MEF90Ctx%resultViewer,step,ierr))
    PetscCallA(VecViewFromOptions(locVecSigma0,PETSC_NULL_OPTIONS,"-ios0_view",ierr))

    ! Cleanup Vec
    PetscCallA(VecDestroy(locVecU,ierr))
    PetscCallA(VecDestroy(locVecU0,ierr))
    PetscCallA(VecDestroy(locVecSigma,ierr))
    PetscCallA(VecDestroy(locVecSigma0,ierr))

    ! Cleanup SF
    PetscCallA(PetscSFDestroy(lcSF,ierr))
    PetscCallA(PetscSFDestroy(clSF,ierr))
    PetscCallA(PetscSFDestroy(lcSSF,ierr))
    PetscCallA(PetscSFDestroy(clSSF,ierr))
    PetscCallA(PetscSFDestroy(lioSF,ierr))
    PetscCallA(PetscSFDestroy(iolSF,ierr))
    PetscCallA(PetscSFDestroy(lioBSF,ierr))
    PetscCallA(PetscSFDestroy(iolBSF,ierr))
    PetscCallA(PetscSFDestroy(lioSSF,ierr))
    PetscCallA(PetscSFDestroy(iolSSF,ierr))
    PetscCallA(PetscSFDestroy(lioBSSF,ierr))
    PetscCallA(PetscSFDestroy(iolBSSF,ierr))

    ! Cleanup DMs
    DeAllocate(time)
    ! Note that I would need to manually destroy these DM no matter what
    PetscCallA(DMDestroy(dm,ierr))

    ! Exit nicely
    PetscCallA(MEF90CtxDestroy(MEF90Ctx,ierr))
    PetscCallA(MEF90Finalize(ierr))
    PetscCallA(PetscFinalize(ierr))
End Program  TestConstraintIO