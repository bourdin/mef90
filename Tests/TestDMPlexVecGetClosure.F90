Module localFunctions
#include <petsc/finclude/petsc.h>
use petsc
use m_MEF90
implicit none

contains
    Subroutine MyVecView(v,ierr)
        Type(tVec),Intent(IN)               :: v
        PetscErrorCode,Intent(INOUT)        :: ierr

        Type(tDM)                           :: dm
        PetscInt                            :: p,pStart,pEnd,numDofClosure
        Character(len=MEF90MXSTRLEN)        :: IOBuffer
        PetscScalar,Dimension(:),Pointer    :: vArray
        PetscInt                            :: height = 0


        PetscCall(VecGetDM(v,dm,ierr))

        PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Cell closure\n",ierr))
        PetscCall(DMPlexGetHeightStratum(dm,height,pStart,pEnd,ierr))
        Write(*,*) 'cells: ',pStart,pEnd
        Do p = pStart,pEnd-1
            PetscCall(MEF90VecGetClosureSize(v,p,numDofClosure,ierr))
            If (numDofClosure > 0) Then
                PetscCall(DMPlexVecGetClosure(dm,PETSC_NULL_SECTION,v,p,vArray,ierr))
                Write(IOBuffer,*) p, vArray,"\n"
                PetscCall(PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr))
                PetscCall(DMPlexVecRestoreClosure(dm,PETSC_NULL_SECTION,v,p,vArray,ierr))
            End If
        End Do
    
        PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Point Values\n",ierr))
        PetscCall(DMPlexGetChart(dm,pStart,pEnd,ierr))
        Do p = pStart,pEnd-1
            PetscCall(MEF90VecGetClosureSize(v,p,numDofClosure,ierr))
            If (numDofClosure > 0) Then
                PetscCall(DMPlexVecGetClosure(dm,PETSC_NULL_SECTION,v,p,vArray,ierr))
                Write(IOBuffer,*) p, vArray,"\n"
                PetscCall(PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr))
                PetscCall(DMPlexVecRestoreClosure(dm,PETSC_NULL_SECTION,v,p,vArray,ierr))
            End If
        End Do
    End Subroutine MyVecView
End Module localFunctions

Program  TestDMPlexVecGetClosure
#include <petsc/finclude/petsc.h>
Use m_MEF90
Use petsc
Use localFunctions
Implicit NONE   
    
    PetscErrorCode                      :: ierr
    Type(MEF90Ctx_Type),target          :: MEF90Ctx
    Type(MEF90CtxGlobalOptions_Type)    :: MEF90GlobalOptions_default

    Type(tDM)                           :: dm
    PetscBool                           :: interpolate = PETSC_TRUE
    Type(tVec)                          :: v
    Type(tPetscSection)                 :: section
    PetscInt                            :: pStart,pEnd,p,numDofClosure
    PetscReal,Dimension(:),Pointer      :: time

    PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
    PetscCallA(MEF90Initialize(ierr))

    MEF90GlobalOptions_default%verbose           = 1_Ki
    MEF90GlobalOptions_default%dryrun            = PETSC_FALSE
    MEF90GlobalOptions_default%timeMin           = 0.0_Kr
    MEF90GlobalOptions_default%timeMax           = 1.0_Kr
    MEF90GlobalOptions_default%timeNumStep       = 11_Ki
    MEF90GlobalOptions_default%timeInterpolation = MEF90TimeInterpolation_linear
    MEF90GlobalOptions_default%timeSkip          = 0_Ki
    MEF90GlobalOptions_default%timeNumCycle      = 1_Ki
    MEF90GlobalOptions_default%elementFamily     = MEF90ElementFamilyLagrange
    MEF90GlobalOptions_default%elementOrder      = 1_Ki
 
    PetscCallA(MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90GlobalOptions_default,ierr))
    PetscCallA(MEF90CtxGetTime(MEF90Ctx,time,ierr))

    PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm,MEF90Ctx%geometryfile,PETSC_NULL_CHARACTER,interpolate,dm,ierr))
    PetscCallA(DMPlexDistributeSetDefault(dm,PETSC_FALSE,ierr))
    PetscCallA(DMSetFromOptions(dm,ierr))
    PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-dm_view",ierr))
    
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

    PetscCallA(PetscSectionCreate(PETSC_COMM_WORLD,section,ierr))
    PetscCallA(DMPlexGetChart(dm,pStart,pEnd,ierr))
    PetscCallA(PetscSectionSetChart(section,pStart,pEnd,ierr))

    PetscCallA(PetscSectionSetDof(section,pStart,1_kI,ierr))

    PetscCallA(PetscSectionSetUp(section,ierr))

    PetscCallA(DMSetLocalSection(dm,section,ierr))
    PetscCallA(PetscObjectViewFromOptions(section,PETSC_NULL_OPTIONS,"-dm_section_view",ierr))

    PetscCallA(DMGetLocalVector(dm,v,ierr))

    PetscCallA(VecSet(v,-1.0_kR,ierr))
    PetscCallA(VecViewFromOptions(v,PETSC_NULL_OPTIONS,"-dm_vec_view",ierr))

    Do p = pStart,pEnd-1
        PetscCall(MEF90VecGetClosureSize(v,p,numDofClosure,ierr))
    End Do
    PetscCallA(MyVecView(v,ierr))

    PetscCallA(PetscSectionDestroy(section,ierr))
    PetscCallA(DMRestoreLocalVector(dm,v,ierr))
    PetscCallA(DMDestroy(dm,ierr))

    Call MEF90CtxDestroy(MEF90Ctx,ierr)   
    Call MEF90Finalize(ierr)
    Call PetscFinalize(ierr)
End Program  TestDMPlexVecGetClosure
 
       
