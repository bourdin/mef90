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
    ! Character(len=MEF90MXSTRLEN)                        :: IOBuffer
    PetscInt                                            :: dim
    Type(tPetscSection)                                 :: sectionU
    Type(tVec)                                          :: U,V,locVecV
    PetscReal,Dimension(:),Pointer                      :: time
    PetscInt                                            :: step = 1_Ki
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

    PetscCallA(VecGetDM(MEF90DefMechCtx%displacementLocal,dmU,ierr))
    PetscCallA(DMGetLocalSection(dmU,sectionU,ierr))
    PetscCallA(DMCreateGlobalVector(dmU,U,ierr))
    PetscCallA(DMCreateGlobalVector(dmU,V,ierr))

    ! Initialize boundary values of MEF90DefMechCtx%displacementLocal with values from the command line
    PetscCallA(VecSet(MEF90DefMechCtx%displacementLocal,-1.23456789_Kr,ierr))
    PetscCallA(MEF90VecSetBCValuesFromOptions(MEF90DefMechCtx%displacementLocal,1.0_Kr,ierr))

    ! PetscCallA(VecView(MEF90DefMechCtx%displacementLocal,PETSC_VIEWER_STDOUT_WORLD,ierr))

    ! Save MEF90DefMechCtx%displacementLocal in exo file
    PetscCallA(MEF90EXOVecView(MEF90DefMechCtx%displacementLocal,MEF90DefMechCtx%displacementToIOSF,MEF90DefMechCtx%IOTodisplacementSF,MEF90Ctx%resultViewer,step,dim,ierr))

    ! ! Create new vectors and read them from the file
    PetscCallA(VecDuplicate(MEF90DefMechCtx%displacementLocal,locVecV,ierr))
    PetscCallA(VecSet(locVecV,-1000.0_kr,ierr))
    PetscCallA(PetscObjectSetName(locVecV,"Displacement",ierr))
    PetscCallA(MEF90EXOVecLoad(locVecV,MEF90DefMechCtx%displacementToIOSF,MEF90DefMechCtx%IOTodisplacementSF,MEF90Ctx%resultViewer,step,dim,ierr))
    PetscCallA(DMLocalToGlobal(dmU,locVecV,INSERT_VALUES,U,ierr))
    PetscCallA(DMGlobalToLocal(dmU,U,INSERT_VALUES,locVecV,ierr))
debug1: block
    PetscReal,Dimension(:),Pointer  :: DisplacementArray, locVecVArray
    Integer                         :: i
    Character(len=MEF90MXSTRLEN)    :: filename
    Type(tPetscSF)                  :: lcgSF,cglSF,llSF
    Type(tDM)                       :: odm
    Type(tVec)                      :: displacementOutput

    ! PetscCallA(DMLocalToGlobal(dmU,locVecV,INSERT_VALUES,V,ierr))
    ! PetscCallA(DMLocalToGlobal(dmU,MEF90DefMechCtx%displacementLocal,INSERT_VALUES,U,ierr))

    write(filename,'("out-",I4.4,".txt")') MEF90Ctx%rank
    Open(file=filename,unit=99)
    Write(*,*) 'Opening ', filename

    ! PetscCallA(VecGetArrayF90(MEF90DefMechCtx%displacementLocal,DisplacementArray,ierr))
    ! PetscCallA(VecGetArrayF90(locVecV,locVecVArray,ierr))

    ! write(99,*) 'Max diff', maxval(DisplacementArray - locVecVArray), maxloc(DisplacementArray - locVecVArray)
    ! Do i = 1, size(locVecVArray)
    !     write(99,*) i, DisplacementArray(i), locVecVArray(i)
    ! End Do
    ! PetscCallA(VecRestoreArrayF90(MEF90DefMechCtx%displacementLocal,DisplacementArray,ierr))
    ! PetscCallA(VecRestoreArrayF90(locVecV,locVecVArray,ierr))

    ! PetscCallA(DMLocalToGlobal(dmU,locVecV,INSERT_VALUES,V,ierr))
    ! PetscCallA(DMLocalToGlobal(dmU,MEF90DefMechCtx%displacementLocal,INSERT_VALUES,U,ierr))

    ! PetscCallA(VecGetArrayF90(U,DisplacementArray,ierr))
    ! PetscCallA(VecGetArrayF90(V,locVecVArray,ierr))

    ! write(99,*) 'Max diff', maxval(DisplacementArray - locVecVArray), maxloc(DisplacementArray - locVecVArray)
    ! Do i = 1, size(locVecVArray)
    !     write(99,*) i, DisplacementArray(i), locVecVArray(i)
    ! End Do

    ! PetscCallA(VecRestoreArrayF90(U,DisplacementArray,ierr))
    ! PetscCallA(VecRestoreArrayF90(V,locVecVArray,ierr))
 
    ! PetscCallA(CreateLocalToCGlobalSF_Private(MEF90Ctx,dmU,lcgSF,ierr))
    ! PetscCallA(CreateCGlobalToLocalSF_Private(MEF90Ctx,dmU,cglSF,ierr))
    ! PetscCallA(PetscSFCompose(lcgSF,cglSF,llSF,ierr))
    ! PetscCallA(PetscSFView(lcgSF,PETSC_VIEWER_STDOUT_WORLD,ierr))
    ! PetscCallA(PetscSFView(cglSF,PETSC_VIEWER_STDOUT_WORLD,ierr))
    ! PetscCallA(MEF90VecCopySF(locVecV,locVecV,llSF,ierr))

    ! PetscCallA(VecGetArrayF90(MEF90DefMechCtx%displacementLocal,DisplacementArray,ierr))
    ! PetscCallA(VecGetArrayF90(locVecV,locVecVArray,ierr))

    ! write(99,*) 'Max diff', maxval(DisplacementArray - locVecVArray), maxloc(DisplacementArray - locVecVArray)
    ! Do i = 1, size(locVecVArray)
    !     write(99,*) i, DisplacementArray(i), locVecVArray(i)
    ! End Do
    ! PetscCallA(VecRestoreArrayF90(MEF90DefMechCtx%displacementLocal,DisplacementArray,ierr))
    ! PetscCallA(VecRestoreArrayF90(locVecV,locVecVArray,ierr))

    PetscCallA(DMGetOutputDM(dmU,oDM,ierr))
    PetscCallA(DMGetGlobalVector(oDM,displacementOutput,ierr))
    PetscCallA(DMLocalToGlobal(oDM,locVecV,INSERT_VALUES,displacementOutput,ierr))
    PetscCallA(DMLocalToGlobal(oDM,MEF90DefMechCtx%displacementLocal,INSERT_VALUES,locVecV,ierr))
    PetscCallA(DMRestoreGlobalVector(oDM,displacementOutput,ierr))

    PetscCallA(VecGetArrayF90(MEF90DefMechCtx%displacementLocal,DisplacementArray,ierr))
    PetscCallA(VecGetArrayF90(locVecV,locVecVArray,ierr))

    write(99,*) 'Max diff', maxval(DisplacementArray - locVecVArray), maxloc(DisplacementArray - locVecVArray)
    Do i = 1, size(locVecVArray)
        write(99,*) i, DisplacementArray(i), locVecVArray(i)
    End Do
    PetscCallA(VecRestoreArrayF90(MEF90DefMechCtx%displacementLocal,DisplacementArray,ierr))
    PetscCallA(VecRestoreArrayF90(locVecV,locVecVArray,ierr))

    Close(99)
end block debug1



! ! Compute the difference between the LOCAL vector we wrote and the one we read
    ! PetscCallA(VecAXPY(locVecV,-1.0_Kr,locVecU,ierr))
    ! PetscCallA(VecNorm(locVecV,NORM_INFINITY,err,ierr))
    ! Write(IOBuffer,'("Local vector L^infty error:  ",ES12.5,"\n")') err
    ! PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
    ! PetscCallA(VecViewFromOptions(locVecV,PETSC_NULL_OPTIONS,"-diffloc_view",ierr))

    ! ! Compute the difference between the LOCAL vector we wrote and the one we read
    ! PetscCallA(VecAXPY(V,-1.0_Kr,U,ierr))
    ! PetscCallA(VecNorm(V,NORM_INFINITY,myerr,ierr))
    ! PetscCallMPIA(MPI_AllReduce(myerr,err,1,MPIU_SCALAR,MPI_MAX,PETSC_COMM_WORLD,ierr))
    ! Write(IOBuffer,'("Global vector L^infty error: ",ES12.5,"\n")') err
    ! PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
    ! PetscCallA(VecViewFromOptions(V,PETSC_NULL_OPTIONS,"-diff_view",ierr))


    ! Cleanup Vec
    PetscCallA(VecDestroy(U,ierr))
    PetscCallA(VecDestroy(V,ierr))

    ! Cleanup DMs
    DeAllocate(time)
    ! Note that I would need to manually destroy these DM no matter what
    PetscCallA(DMDestroy(dm,ierr))

    ! Exit nicely
    PetscCallA(MEF90CtxDestroy(MEF90Ctx,ierr))
    PetscCallA(MEF90Finalize(ierr))
    PetscCallA(PetscFinalize(ierr))
End Program  TestConstraintIO3
