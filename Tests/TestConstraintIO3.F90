Module localFunctions
#include <petsc/finclude/petsc.h>
use m_MEF90
use m_MEF90_DefMech
use m_vDefDefault
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

Program  TestConstraintIO3
#include <petsc/finclude/petsc.h>
Use m_MEF90
Use petsc
Use localFunctions
Implicit NONE   
    
    PetscErrorCode                                      :: ierr
    Type(MEF90Ctx_Type),target                          :: MEF90Ctx
    Type(MEF90CtxGlobalOptions_Type)                    :: MEF90GlobalOptions_default
    Type(MEF90CtxGlobalOptions_Type),Pointer            :: MEF90GlobalOptions
    Type(tDM),target                                    :: dm,dmU
    PetscBool                                           :: interpolate = PETSC_TRUE

   !!! Defect mechanics contexts
    Type(MEF90DefMechCtx_Type)                          :: MEF90DefMechCtx
    Type(MEF90DefMechGlobalOptions_Type),pointer        :: MEF90DefMechGlobalOptions

    PetscInt                                            :: numNodalVar = 3, numCellVar = 0, numGVar = 0, numSideVar = 0
    Character(len=MEF90MXSTRLEN),Dimension(:),Pointer   :: nodalVarName, cellVarName, gVarName, sideVarName
    Character(len=MEF90MXSTRLEN)                        :: filename
    PetscInt                                            :: dim
    Type(tPetscSection)                                 :: sectionU
    Type(tVec)                                          :: U,V,locVecV
    PetscReal,Dimension(:),Pointer                      :: time
    PetscInt                                            :: i,step = 1_Ki
    PetscReal,Dimension(:),Pointer                      :: DisplacementArray, locVecVArray
    ! PetscReal                                           :: myerr,err

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
    
    ! Create DM from file
    PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm,MEF90Ctx%geometryfile,PETSC_NULL_CHARACTER,interpolate,dm,ierr))
    PetscCallA(DMPlexDistributeSetDefault(dm,PETSC_FALSE,ierr))
    PetscCallA(DMSetUseNatural(dm,PETSC_TRUE,ierr))
    PetscCallA(DMSetFromOptions(dm,ierr))
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-mef90dm_view",ierr))

    PetscCallA(DMGetDimension(dm,dim,ierr))

    ! Open exodus file + write geometry + format the file
    numNodalVar = dim
    numCellVar  = 0
    numGVar     = 0
    numSideVar  = 0

    Allocate(nodalVarName(numNodalVar))
    Allocate(cellVarName(numCellVar))
    Allocate(gVarName(numGVar))
    Allocate(sideVarName(numSideVar))
    If (dim == 2) Then
        nodalVarName = ["Displacement_X","Displacement_Y"]
    Else
        nodalVarName = ["Displacement_X","Displacement_Y","Displacement_Z"]
    End If

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
        Type(tPetscSF)                      :: naturalPointSF
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

    PetscCallA(MEF90DefMechCtxCreate(MEF90DefMechCtx,dm,MEF90Ctx,ierr))
    PetscCallA(MEF90DefMechCtxSetFromOptions(MEF90DefMechCtx,PETSC_NULL_CHARACTER,DefMechDefaultGlobalOptions,DefMechDefaultCellSetOptions,DefMechDefaultFaceSetOptions,DefMechDefaultVertexSetOptions,ierr))
    PetscCallA(PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr))
    PetscCallA(DMDestroy(dm,ierr))

    PetscCallA(VecGetDM(MEF90DefMechCtx%displacementLocal,dmU,ierr))
    PetscCallA(DMGetLocalSection(dmU,sectionU,ierr))
    PetscCallA(DMGetGlobalVector(dmU,U,ierr))
    PetscCallA(DMGetGlobalVector(dmU,V,ierr))

    ! ! Initialize boundary values of MEF90DefMechCtx%displacementLocal with values from the command line
    ! PetscCallA(VecSet(MEF90DefMechCtx%displacementLocal,-1.23456789_Kr,ierr))
    ! PetscCallA(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%displacementLocal,1.0_Kr,ierr))

    ! project a field onto a local vector
    PetscCallA(project(MEF90DefMechCtx%displacementLocal,sectionU,ierr))


    ! Save MEF90DefMechCtx%displacementLocal in exo file
    PetscCallA(MEF90EXOVecView(MEF90DefMechCtx%displacementLocal,MEF90DefMechCtx%displacementToIOSF,MEF90DefMechCtx%IOTodisplacementSF,MEF90Ctx%resultViewer,step,dim,ierr))

    ! ! Create new vectors and read them from the file
    PetscCallA(VecDuplicate(MEF90DefMechCtx%displacementLocal,locVecV,ierr))
    PetscCallA(VecSet(locVecV,-1000.0_kr,ierr))
    PetscCallA(PetscObjectSetName(locVecV,"Displacement",ierr))
    PetscCallA(MEF90EXOVecLoad(locVecV,MEF90DefMechCtx%displacementToIOSF,MEF90DefMechCtx%IOTodisplacementSF,MEF90Ctx%resultViewer,step,dim,ierr))

    write(filename,'("out-",I4.4,".txt")') MEF90Ctx%rank
    Open(file=filename,unit=99)
    Write(*,*) 'Opening ', filename

    PetscCallA(VecGetArrayReadF90(MEF90DefMechCtx%displacementLocal,DisplacementArray,ierr))
    PetscCallA(VecGetArrayReadF90(locVecV,locVecVArray,ierr))

    write(99,*) 'Max diff', maxval(DisplacementArray - locVecVArray), maxloc(DisplacementArray - locVecVArray)
    Do i = 1, size(locVecVArray)
        if (abs(DisplacementArray(i) - locVecVArray(i)) > 1.e-5) then
            write(99,*) i, DisplacementArray(i), locVecVArray(i)
        end if
    End Do
    PetscCallA(VecRestoreArrayReadF90(MEF90DefMechCtx%displacementLocal,DisplacementArray,ierr))
    PetscCallA(VecRestoreArrayReadF90(locVecV,locVecVArray,ierr))

    PetscCallA(VecSet(U,-1.0_Kr,ierr))
    PetscCallA(VecSet(V,1.0_Kr,ierr))
    PetscCallA(DMLocalToGlobal(dmU,MEF90DefMechCtx%displacementLocal,INSERT_VALUES,U,ierr))
    PetscCallA(DMLocalToGlobal(dmU,locVecV,INSERT_VALUES,V,ierr))
    PetscCallA(VecGetArrayReadF90(U,DisplacementArray,ierr))
    PetscCallA(VecGetArrayReadF90(V,locVecVArray,ierr))

    write(99,*) 'Max diff', maxval(DisplacementArray - locVecVArray), maxloc(DisplacementArray - locVecVArray)
    Do i = 1, size(locVecVArray)
        if (abs(DisplacementArray(i) - locVecVArray(i)) > 1.e-5) then
            write(99,*) i, DisplacementArray(i), locVecVArray(i)
        end if
    End Do
    PetscCallA(VecRestoreArrayReadF90(U,DisplacementArray,ierr))
    PetscCallA(VecRestoreArrayReadF90(V,locVecVArray,ierr))
    Close(99)

    ! Cleanup Vec
    PetscCallA(DMRestoreGlobalVector(dmU,U,ierr))
    PetscCallA(DMRestoreGlobalVector(dmU,V,ierr))

    ! Cleanup DMs
    DeAllocate(time)

    ! Exit nicely
    PetscCallA(MEF90CtxDestroy(MEF90Ctx,ierr))
    PetscCallA(MEF90Finalize(ierr))
    PetscCallA(PetscFinalize(ierr))
End Program  TestConstraintIO3
