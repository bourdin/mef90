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

Program  TestConstraintIO2
#include <petsc/finclude/petsc.h>
Use m_MEF90
Use petsc
Use localFunctions
Implicit NONE   
    
    PetscErrorCode                                      :: ierr
    Type(MEF90Ctx_Type),target                          :: MEF90Ctx
    Type(MEF90CtxGlobalOptions_Type)                    :: MEF90GlobalOptions_default
    Type(MEF90CtxGlobalOptions_Type),Pointer            :: MEF90GlobalOptions
    Type(tDM),target                                    :: dm,dmU,dmU0
    PetscBool                                           :: interpolate = PETSC_TRUE

    PetscInt                                            :: numNodalVar = 3, numCellVar = 0, numGVar = 0, numSideVar = 0
    Character(len=MEF90MXSTRLEN),Dimension(:),Pointer   :: nodalVarName, cellVarName, gVarName, sideVarName
    Character(len=MEF90MXSTRLEN)                        :: name,IOBuffer
    PetscInt                                            :: dim,pStart,pEnd
    Type(tPetscSection)                                 :: sectionU, sectionU0
    Type(tPetscSF)                                      :: naturalPointSF,lcSF, clSF, lcSSF, clSSF, lioSF, iolSF, lioBSF, iolBSF, lioSSF, iolSSF, lioBSSF, iolBSSF
    Type(tVec)                                          :: locCoord, locVecU, locVecU0, U, U0, locVecV, V
    PetscReal,Dimension(:),Pointer                      :: time
    PetscInt                                            :: step = 1_Ki
    PetscReal                                           :: err

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
    nodalVarName = ["U_X","U_Y","U_Z"]
    
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
            PetscCallA(DMSetUseNatural(dm,PETSC_TRUE,ierr))
        End If
    End Block distribute
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-mef90dm_view",ierr))

    PetscCallA(DMGetDimension(dm,dim,ierr))
    PetscCallA(DMPlexGetChart(dm,pStart,pEnd,ierr))

    ! Create nodal local Vec holding coordinates
    name = "U"
    PetscCallA(MEF90CreateLocalVector(dm,MEF90GlobalOptions%elementFamily,MEF90GlobalOptions%elementOrder,dim,name,locVecU,ierr))
    PetscCallA(VecGetDM(locVecU,dmU,ierr))
    PetscCallA(DMGetLocalSection(dmU,sectionU,ierr))
    PetscCallA(DMCreateGlobalVector(dmU,U,ierr))

    ! Create nodal local Vec holding constraints
    name = "U0"
    PetscCallA(MEF90CreateBoundaryLocalVector(dm,MEF90GlobalOptions%elementFamily,MEF90GlobalOptions%elementOrder,dim,name,locVecU0,ierr))
    PetscCallA(VecGetDM(locVecU0,dmU0,ierr))
    PetscCallA(DMGetLocalSection(dmU0,sectionU0,ierr))
    PetscCallA(DMCreateGlobalVector(dmU0,U0,ierr))

    ! Create SFs for copying from/into IO coordinates Vec
    PetscCallA(MEF90IOSFCreate(MEF90Ctx,locVecU,lioSF,iolSF,ierr))
    ! Create SFs for copying from/into IO Constraint Vec
    PetscCallA(MEF90IOSFCreate(MEF90Ctx,locVecU0,lioBSF,iolBSF,ierr))
    ! Create SFs for copying constrained dofs from/into Constraint Vec
    PetscCallA(MEF90ConstraintSFCreate(MEF90Ctx,locVecU,locVecU0,lcSF,clSF,ierr))

    ! Fill locVecU holding coordinates and copy constraint values = 10 from locVecU0
    PetscCallA(VecSet(locVecU0,10.0_kr,ierr))
    PetscCallA(DMGetCoordinatesLocal(dmU0,locCoord,ierr))

    ! Fill locVecU with coordinates
    PetscCallA(project(locVecU,sectionU,ierr))
    PetscCallA(DMLocalToGlobal(dmU,locVecU,INSERT_VALUES,U,ierr))
    PetscCallA(VecViewFromOptions(locVecU,PETSC_NULL_OPTIONS,"-Uloc_view",ierr))
    PetscCallA(VecViewFromOptions(U,PETSC_NULL_OPTIONS,"-U_view",ierr))

    ! Save locVecU in exo file
    PetscCallA(MEF90EXOVecView(locVecU,lioSF,iolSF,MEF90Ctx%resultViewer,step,dim,ierr))

    ! Create new vectors and read them from the file
    PetscCallA(VecDuplicate(locVecU,locVecV,ierr))
    PetscCallA(VecDuplicate(U,V,ierr))
    PetscCallA(VecSet(locVecV,-1000.0_kr,ierr))
    PetscCallA(PetscObjectSetName(locVecV,"U",ierr))
    PetscCallA(MEF90EXOVecLoad(locVecV,lioSF,iolSF,MEF90Ctx%resultViewer,step,dim,ierr))
    PetscCallA(DMLocalToGlobal(dmU,locVecV,INSERT_VALUES,V,ierr))

    PetscCallA(VecViewFromOptions(locVecV,PETSC_NULL_OPTIONS,"-Vloc_view",ierr))
    PetscCallA(VecViewFromOptions(V,PETSC_NULL_OPTIONS,"-V_view",ierr))

    ! Save it again 
    PetscCallA(MEF90EXOVecView(locVecU,lioSF,iolSF,MEF90Ctx%resultViewer,step+1,dim,ierr))



    ! Compute the difference between the LOCAL vector we wrote and the one we read
    PetscCallA(VecAXPY(locVecV,-1.0_Kr,locVecU,ierr))
    PetscCallA(VecNorm(locVecV,NORM_INFINITY,err,ierr))
    Write(IOBuffer,'("Local vector L^infty error: ",ES12.5,"\n")') err
    PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))

    ! Compute the difference between the LOCAL vector we wrote and the one we read
    PetscCallA(VecAXPY(V,-1.0_Kr,U,ierr))
    PetscCallA(VecNorm(V,NORM_INFINITY,err,ierr))
    Write(IOBuffer,'("Global vector L^infty error: ",ES12.5,"\n")') err
    PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))


    
    ! Cleanup Vec
    PetscCallA(VecDestroy(locVecU,ierr))
    PetscCallA(VecDestroy(locVecU0,ierr))
    PetscCallA(VecDestroy(locVecV,ierr))

    PetscCallA(VecDestroy(U,ierr))
    PetscCallA(VecDestroy(U0,ierr))
    PetscCallA(VecDestroy(V,ierr))

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
End Program  TestConstraintIO2
