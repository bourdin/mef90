Module m_MEF90_DMPlex
#include "petsc/finclude/petsc.h"
    Use m_MEF90_Elements
    Use m_MEF90_Ctx
    Use petsc
    Use,intrinsic :: iso_c_binding
    IMPLICIT NONE
    
    Enum,bind(c)
        enumerator  :: MEF90_DMPlexCellSetType = 1, &
                       MEF90_DMPlexFaceSetType,     &
                       MEF90_DMPlexEdgeSetType,     &
                       MEF90_DMPlexVertexSetType
    End Enum
    Character(len=MEF90_MXSTRLEN),Parameter :: MEF90_DMPlexCellSetLabelName   = 'Cell Sets  '
    Character(len=MEF90_MXSTRLEN),Parameter :: MEF90_DMPlexFaceSetLabelName   = 'Face Sets  '
    Character(len=MEF90_MXSTRLEN),Parameter :: MEF90_DMPlexEdgeSetLabelName   = 'Edge Sets  '
    Character(len=MEF90_MXSTRLEN),Parameter :: MEF90_DMPlexVertexSetLabelName = 'Vertex Sets'

    Character(len=MEF90_MXSTRLEN),dimension(4),Parameter  :: MEF90_DMPlexSetLabelName = &
              [MEF90_DMPlexCellSetLabelName,MEF90_DMPlexFaceSetLabelName,MEF90_DMPlexEdgeSetLabelName,MEF90_DMPlexVertexSetLabelName]

    private
    public :: MEF90_DMPlexCellSetLabelName,MEF90_DMPlexFaceSetLabelName,             &
              MEF90_DMPlexEdgeSetLabelName,MEF90_DMPlexVertexSetLabelName,           &
              MEF90_DMPlexSetLabelName,                                              &
              MEF90_DMPlexCellSetType,MEF90_DMPlexFaceSetType,                       &
              MEF90_DMPlexEdgeSetType,MEF90_DMPlexVertexSetType,                     &
              MEF90_SectionAllocateDof,MEF90_SectionAllocateDofSet,                  &
              MEF90_SetupConstraintTableSet,MEF90_SectionAllocateConstraint,         &
              MEF90_VecReorderingSF,MEF90_CreateLocalToIOSF,MEF90_CreateIOToLocalSF, &
              MEF90_CreateLocalToConstraintSF,MEF90_VecGlobalToLocalConstraint
Contains

#undef __FUNCT__
#define __FUNCT__ "MEF90_SectionAllocateDof"
!!!
!!!  
!!!  MEF90_SectionAllocateDof: Associates dof to a section
!!!  
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

    Subroutine MEF90_SectionAllocateDof(dm,setType,elemType,numComponents,section,ierr)
        Type(tDM),Intent(IN)               :: dm
        PetscEnum,intent(IN)               :: setType
        Type(MEF90Element_Type),Intent(IN) :: elemType
        PetscInt,Intent(IN)                :: numComponents
        Type(tPetscSection),Intent(INOUT)  :: section
        PetscErrorCode,Intent(INOUT)       :: ierr

        Type(tIS)                          :: setIS
        PetscInt,dimension(:),Pointer      :: setID
        PetscInt                           :: set

        PetscCall(DMGetLabelIdIS(dm,MEF90_DMPlexSetLabelName(setType),setIS,ierr))
        PetscCall(ISGetIndicesF90(setIS,setID,ierr))
        Do set = 1,size(setID)
            PetscCall(MEF90_SectionAllocateDofSet(dm,setType,setID(set),elemType,numComponents,Section,ierr))
        End Do
        PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
        PetscCall(ISDestroy(setIS,ierr))
    End Subroutine MEF90_SectionAllocateDof

#undef __FUNCT__
#define __FUNCT__ "MEF90_SectionAllocateDofSet"
!!!
!!!  
!!!  MEF90_SectionAllocateDofSet: Associates dof to a section in a set
!!!  
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

    Subroutine MEF90_SectionAllocateDofSet(dm,setType,setID,elemType,numComponents,section,ierr)
        Type(tDM),Intent(IN)               :: dm
        PetscEnum,intent(IN)               :: setType
        PetscInt,Intent(IN)                :: setID
        Type(MEF90Element_Type),Intent(IN) :: elemType
        PetscInt,Intent(IN)                :: numComponents
        Type(tPetscSection),Intent(INOUT)  :: section
        PetscErrorCode,Intent(INOUT)       :: ierr

        Type(tIS)                          :: setPointIS
        PetscInt,dimension(:),Pointer      :: setPointID
        PetscInt                           :: point
        PetscInt,dimension(:),pointer      :: closure
        PetscInt                           :: p,depth


        PetscCall(DMGetStratumIS(dm,MEF90_DMPlexSetLabelName(setType),setID,setPointIS,ierr))
        PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
        If (size(setPointID) > 0) Then
            !!! This can probably be optimized by allocating closure outside of the loop
            !!! But I can't figure out how it is done at the moment.
            Nullify(closure)
            Do point = 1,size(setPointID)
                PetscCall(DMPlexGetTransitiveClosure(dm,setPointID(point),PETSC_TRUE,closure,ierr))
                Do p = 1,size(closure),2
                    PetscCall(DMPlexGetPointDepth(dm,closure(p),depth,ierr))
                    If (elemType%numDofs(depth+1) > 0) Then
                        PetscCall(PetscSectionSetDof(section,closure(p),elemType%numDofs(depth+1)*numComponents,ierr))
                        PetscCall(PetscSectionSetFieldDof(section,closure(p),[0],elemType%numDofs(depth+1)*numComponents,ierr))
                    End If
                End Do! p
                PetscCall(DMPlexRestoreTransitiveClosure(dm,setPointID(point),PETSC_TRUE,closure,ierr))
            End Do! cell
        End If
        PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
        PetscCall(ISDestroy(setPointIS,ierr))
    End Subroutine MEF90_SectionAllocateDofSet

#undef __FUNCT__
#define __FUNCT__ "MEF90_SetupConstraintTableSet"
!!!
!!!  
!!!  MEF90_SetupConstraintTableSet: Build the contribution of a set to the constraint table
!!!  
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!
    
    Subroutine MEF90_SetupConstraintTableSet(dm,section,setType,setID,constraints,table,ierr)
        Type(tDM),Intent(IN)               :: dm
        Type(tPetscSection),Intent(INOUT)  :: section
        PetscEnum,intent(IN)               :: setType
        PetscInt,Intent(IN)                :: setID
        Logical,Dimension(:),Pointer       :: constraints
        Logical,Dimension(:,:),Pointer     :: table
        PetscErrorCode,Intent(INOUT)       :: ierr

        Type(tIS)                          :: setPointIS
        PetscInt,dimension(:),Pointer      :: setPointID
        PetscInt                           :: point,numDof
        PetscInt,dimension(:),pointer      :: closure
        PetscInt                           :: p

        PetscCall(DMGetStratumIS(dm,MEF90_DMPlexSetLabelName(setType),setID,setPointIS,ierr))
        PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
        If (size(setPointID) > 0) Then
            !!! This can probably be optimized by allocating closure outside of the loop
            !!! But I can't figure out how it is done at the moment.
            Nullify(closure)
            Do point = 1,size(setPointID)
                PetscCall(DMPlexGetTransitiveClosure(dm,setPointID(point),PETSC_TRUE,closure,ierr))
                Do p = 1,size(closure),2
                    PetscCall(PetscSectionGetDoF(section,closure(p),numDof,ierr))
                    If (numDof > 0) Then
                        table(closure(p)+1,:) = table(closure(p)+1,:) .OR. constraints
                    End If ! numDof
                End Do! p
                PetscCall(DMPlexRestoreTransitiveClosure(dm,setPointID(point),PETSC_TRUE,closure,ierr))
            End Do! cell
        End If
        PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
        PetscCall(ISDestroy(setPointIS,ierr))
    End Subroutine MEF90_SetupConstraintTableSet

#undef __FUNCT__
#define __FUNCT__ "MEF90_SectionAllocateConstraint"
!!!
!!!  
!!!  MEF90_SectionAllocateConstraint: Associates dof to a section
!!!  
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

    Subroutine MEF90_SectionAllocateConstraint(dm,table,section,ierr)
        Type(tDM),Intent(IN)               :: dm
        Logical,Dimension(:,:),Pointer     :: table
        Type(tPetscSection),Intent(INOUT)  :: section
        PetscErrorCode,Intent(INOUT)       :: ierr

        PetscInt                           :: p,i,pStart,pEnd,numConstraints,numComponents
        PetscInt,Dimension(:),Pointer      :: constraints
        PetscInt                           :: field = 0

        PetscCall(DMPlexGetChart(dm,pStart,pEnd,ierr))
        Do p = 1, pEnd
            numConstraints = count(table(p,:))
            If (numConstraints > 0) Then
                PetscCall(PetscSectionSetConstraintDof(section,p-1,numConstraints,ierr))
                PetscCall(PetscSectionSetFieldConstraintDof(section,p-1,[0],numConstraints,ierr))
            End If
        End Do

        PetscCall(PetscSectionSetup(section,ierr))

        numComponents = size(table,2)
        Do p = 1, pEnd
            numConstraints = count(table(p,:))
            If (numConstraints > 0) Then
                Allocate(constraints(numConstraints))
                constraints = pack([ (i-1, i = 1,numComponents) ],table(p,:))
                PetscCall(PetscSectionSetConstraintIndicesF90(section,p-1,constraints,ierr))
                PetscCall(PetscSectionSetFieldConstraintIndicesF90(section,p-1,field,constraints,ierr))
                DeAllocate(constraints)
            End If
        End Do
    End Subroutine MEF90_SectionAllocateConstraint

#undef __FUNCT__
#define __FUNCT__ "MEF90_VecReorderingSF"
!!!
!!!  
!!!  MEF90_VecReorderingSF: rearrange a Vec according to the given SF
!!!  
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!

    Subroutine MEF90_VecReorderingSF(vin,vout,sf,ierr)
        Type(tVec),intent(IN)              :: vin
        Type(tVec),intent(INOUT)           :: vout
        Type(tPetscSF),intent(IN)          :: sf
        PetscErrorCode,intent(INOUT)       :: ierr

        PetscScalar,dimension(:),Pointer   :: arrayin, arrayout

        PetscCallA(VecGetArrayReadF90(vin, arrayin, ierr))
        PetscCallA(VecGetArrayF90(vout, arrayout, ierr))
        PetscCallA(PetscSFBcastBegin(sf, MPIU_SCALAR, arrayin, arrayout, MPI_REPLACE, ierr))
        PetscCallA(PetscSFBcastEnd(sf, MPIU_SCALAR, arrayin, arrayout, MPI_REPLACE, ierr))
        PetscCallA(VecRestoreArrayReadF90(vin, arrayin, ierr))
        PetscCallA(VecRestoreArrayF90(vout, arrayout, ierr))
    End subroutine MEF90_VecReorderingSF

#undef __FUNCT__
#define __FUNCT__ "MEF90_CreateLocalToIOSF"
!!!
!!!  
!!!  MEF90_CreateLocalToIOSF: sf mapping between local to IO ordering and distribution
!!!  
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!

    Subroutine MEF90_CreateLocalToIOSF(MEF90Ctx,dm,sf,ierr)
        Type(tDM),intent(IN)               :: dm
        Type(tPetscSF),intent(OUT)         :: sf
        Type(MEF90Ctx_type),Intent(IN)     :: MEF90Ctx
        PetscErrorCode,intent(INOUT)       :: ierr

        Type(tPetscSF)                     :: naturalSF, ioSF, lcgSF, tempSF

        PetscCallA(CreateLocalToCGlobalSF_internal(MEF90Ctx,dm,lcgSF,ierr))
        If (MEF90Ctx%NumProcs > 1) Then
            PetscCallA(CreateNaturalToIOSF_internal(MEF90Ctx,dm,ioSF,ierr))
            PetscCallA(DMGetNaturalSF(dm,naturalSF,ierr))
            PetscCallA(PetscSFCompose(lcgSF,naturalSF,tempSF,ierr))
            PetscCallA(PetscSFCompose(tempSF,ioSF,sf,ierr))
            PetscCallA(PetscSFDestroy(lcgSF,ierr))
            PetscCallA(PetscSFDestroy(ioSF,ierr))
            PetscCallA(PetscSFDestroy(tempSF,ierr))
            PetscCallA(PetscSFSetUp(sf,ierr))
        Else
            sf = lcgSF
        End If
    End subroutine MEF90_CreateLocalToIOSF

#undef __FUNCT__
#define __FUNCT__ "MEF90_CreateIOToLocalSF"
!!!
!!!  
!!!  MEF90_CreateIOToLocalSF: sf mapping between IO to local ordering and distribution
!!!  
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!

    Subroutine MEF90_CreateIOToLocalSF(MEF90Ctx,dm,sf,ierr)
        Type(tDM),intent(IN)               :: dm
        Type(tPetscSF),intent(OUT)         :: sf
        Type(MEF90Ctx_type),Intent(IN)     :: MEF90Ctx
        PetscErrorCode,intent(INOUT)       :: ierr

        Type(tPetscSF)                     :: naturalSF, ioSF, cglSF, tempSF, invTempSF

        PetscCallA(CreateCGlobalToLocalSF_internal(MEF90Ctx,dm,cglSF,ierr))
        If (MEF90Ctx%NumProcs > 1) Then
            PetscCallA(CreateNaturalToIOSF_internal(MEF90Ctx,dm,ioSF,ierr))
            PetscCallA(DMGetNaturalSF(dm,naturalSF,ierr))
            PetscCallA(PetscSFCompose(naturalSF,ioSF,tempSF,ierr))
            PetscCallA(PetscSFCreateInverseSF(tempSF,invTempSF,ierr))
            PetscCallA(PetscSFCompose(invTempSF,cglSF,sf,ierr))
            PetscCallA(PetscSFDestroy(ioSF,ierr))
            PetscCallA(PetscSFDestroy(cglSF,ierr))
            PetscCallA(PetscSFDestroy(tempSF,ierr))
            PetscCallA(PetscSFDestroy(invTempSF,ierr))
            PetscCallA(PetscSFSetUp(sf,ierr))
        Else
            sf = cglSF
        End If
    End subroutine MEF90_CreateIOToLocalSF

#undef __FUNCT__
#define __FUNCT__ "MEF90_CreateLocalToConstraintSF"
!!!
!!!  
!!!  MEF90_CreateLocalToConstraintSF: sf mapping between local and constraint Vec ordering and distribution
!!!  
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!                Blaise Bourdin  bourdin@mcmaster.ca
!!!

    Subroutine MEF90_CreateLocalToConstraintSF(MEF90Ctx,dm,dmB,sf,invSF,ierr)
        Type(tDM),intent(IN)                    :: dm, dmB
        Type(tPetscSF),intent(OUT)              :: sf, invSF
        Type(MEF90Ctx_type),Intent(IN)          :: MEF90Ctx
        PetscErrorCode,intent(INOUT)            :: ierr

        Type(tPetscSection)                     :: locSection, locBSection
        Type(PetscSFNode),dimension(:),Pointer  :: remote
        Type(PetscInt),dimension(:),Pointer     :: local, cindices
        PetscInt                                :: pStart, pEnd, p, d, nleaves = 0, ldof, loff, cdof, coff, nsize = 0, nroots
        PetscMPIInt                             :: rank

        nleaves = 0
        nsize   = 0
        PetscCallMPI(MPI_Comm_rank(MEF90Ctx%Comm, rank, ierr))
        PetscCall(DMGetLocalSection(dm,locSection,ierr))
        PetscCall(PetscSectionGetStorageSize(locSection,nroots,ierr))
        PetscCall(DMGetLocalSection(dmB,locBSection,ierr))
        PetscCall(PetscSectionGetChart(locBSection, pStart, pEnd, ierr))
        Do p = pStart,pEnd-1
            PetscCall(PetscSectionGetConstraintDof(locBSection,p,cdof,ierr))
            nleaves = nleaves + cdof
        End Do
        Allocate(remote(nleaves))
        Allocate(local(nleaves))
        Do p = pStart,pEnd-1
            PetscCall(PetscSectionGetDoF(locSection,p,ldof,ierr))
            PetscCall(PetscSectionGetOffset(locSection,p,loff,ierr))
            PetscCall(PetscSectionGetConstraintDof(locBSection,p,cdof,ierr))
            If (cdof > 0) then
                PetscCall(PetscSectionGetConstraintIndicesF90(locBSection,p,cindices,ierr))
                PetscCall(PetscSectionGetOffset(locBSection,p,coff,ierr))
                If (coff >= 0) Then
                    Do d=1,cdof
                        local(nsize+1) = coff+cindices(d)
                        remote(nsize+1)%rank = rank
                        remote(nsize+1)%index = loff+cindices(d)
                        nsize = nsize + 1
                    End Do
                End If
                PetscCall(PetscSectionRestoreConstraintIndicesF90(locBSection,p,cindices,ierr))
            End If ! cdof
        End Do
        PetscCall(PetscSFCreate(MEF90Ctx%Comm, sf, ierr))
        PetscCall(PetscSFSetFromOptions(sf, ierr))
        PetscCall(PetscSFSetGraph(sf, nroots, nleaves, local, PETSC_COPY_VALUES, remote, PETSC_COPY_VALUES, ierr))
        DeAllocate(remote)
        DeAllocate(local)
        PetscCall(PetscSFSetUp(sf, ierr))
        PetscCall(PetscSFCreateInverseSF(sf,invSF,ierr))
        PetscCall(PetscSFSetUp(invSF, ierr))
    End Subroutine MEF90_CreateLocalToConstraintSF

#undef __FUNCT__
#define __FUNCT__ "MEF90_VecGlobalToLocalConstraint"
!!!
!!!  
!!!  MEF90_VecGlobalToLocalConstraint: do a VecGlobaltoLocal then copy constrained values
!!!  
!!!  (c) 2022      Blaise Bourdin  bourdin@mcmaster.ca
!!!

    Subroutine MEF90_VecGlobalToLocalConstraint(g,c,l,ierr)
        Type(tVec),Intent(IN)                   :: g,c
        Type(tVec),Intent(INOUT)                :: l
        PetscErrorCode,Intent(INOUT)            :: ierr

        Type(tDM)                               :: dm
        Type(tPetscSection)                     :: s
        PetscInt                                :: numConstraint
        PetscInt                                :: p,pStart,pEnd
        PetscReal,Dimension(:),Pointer          :: vArray

        PetscCall(VecGetDM(g,dm,ierr))
        PetscCall(DMGlobalToLocal(dm,g,INSERT_VALUES,l,ierr))
        If (c /= PETSC_NULL_VEC) Then
            PetscCall(DMGetLocalSection(dm,s,ierr))
            PetscCallA(PetscSectionGetChart(s,pStart,pEnd,ierr))
            Do p = pStart,pEnd-1
                PetscCall(PetscSectionGetConstraintDof(s,p,numConstraint,ierr))
                If (numConstraint > 0) Then
                    !PetscCallA(VecGetValuesSectionF90(c,s,p,vArray,ierr))
                    PetscCallA(VecSetValuesSectionF90(l,s,p,vArray,INSERT_ALL_VALUES,ierr))
                    !PetscCallA(VecRestoreValuesSectionF90(c,s,p,vArray,ierr))
                End If ! numConstraint
            End Do ! p
        End If
    End Subroutine MEF90_VecGlobalToLocalConstraint

#undef __FUNCT__
#define __FUNCT__ "CreateNaturalToIOSF_internal"
!!!
!!!  
!!!  CreateNaturalToIOSF_internal: 
!!!  
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!
    
    subroutine CreateNaturalToIOSF_internal(MEF90Ctx,dm,sf,ierr)
        Type(tDM),intent(IN)                    :: dm
        Type(MEF90Ctx_type),Intent(IN)          :: MEF90Ctx
        Type(tPetscSF),intent(OUT)              :: sf
        PetscErrorCode,intent(INOUT)            :: ierr
    
        Type(tVec)                              :: vnat, vio
        Type(PetscLayout)                       :: natMap, ioMap
        Type(PetscSFNode),dimension(:),Pointer  :: remote
        PetscInt,dimension(:),Pointer           :: ioRange, remoteRange
        PetscMPIInt                             :: rank, remoteRank
        PetscInt                                :: nroots, nleaves, globalIndex, i, globalSize
    
        if (MEF90Ctx%NumProcs > 1) then
!            PetscCallA(DMPlexCreateNaturalVector(dm,vnat,ierr))
            Continue
        else
            PetscCallA(DMCreateGlobalVector(dm,vnat,ierr))
        end if
        PetscCallA(VecGetSize(vnat, globalSize, ierr))
        PetscCallA(VecCreateMPI(MEF90Ctx%Comm, PETSC_DECIDE, globalSize, vio, ierr))
        PetscCallA(VecGetLayout(vnat, natMap, ierr))
        PetscCallA(VecGetLayout(vio, ioMap, ierr))
        PetscCallMPI(MPI_Comm_rank(MEF90Ctx%Comm, rank, ierr))
        PetscCallA(PetscLayoutGetLocalSize(natMap, nroots, ierr))
        PetscCallA(PetscLayoutGetLocalSize(ioMap, nleaves, ierr))
        PetscCallA(PetscLayoutGetRangesF90(ioMap, ioRange, ierr))
        PetscCallA(PetscLayoutGetRangesF90(natMap, remoteRange, ierr))
        Allocate(remote(nleaves))
        Do i = 0,nleaves-1          
            globalIndex = ioRange(rank+1) + i
            PetscCallA(PetscLayoutFindOwner(natMap, globalIndex, remoteRank, ierr))
            remote(i+1)%rank = remoteRank
            remote(i+1)%index = globalIndex - remoteRange(remoteRank+1)
        End Do
        PetscCallA(PetscSFCreate(MEF90Ctx%Comm, sf, ierr))
        PetscCallA(PetscObjectSetName(sf, "Natural-To-IO SF", ierr))
        PetscCallA(PetscSFSetFromOptions(sf, ierr))
        PetscCallA(PetscSFSetGraph(sf, nroots, nleaves, PETSC_NULL_INTEGER, PETSC_COPY_VALUES, remote, PETSC_COPY_VALUES, ierr))
        PetscCallA(PetscSFSetUp(sf, ierr))
        PetscCallA(PetscSFViewFromOptions(sf,PETSC_NULL_OPTIONS,"-naturaltoio_sf_view",ierr))
        PetscCallA(VecDestroy(vio,ierr))
        PetscCallA(VecDestroy(vnat,ierr))
        DeAllocate(remote)
    End subroutine CreateNaturalToIOSF_internal
    
#undef __FUNCT__
#define __FUNCT__ "CreateLocalToCGlobalSF_internal"
!!!
!!!  
!!!  CreateLocalToCGlobalSF_internal: 
!!!  
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!
    
    subroutine CreateLocalToCGlobalSF_internal(MEF90Ctx,dm,sf,ierr)
        Type(tDM),intent(IN)                    :: dm
        Type(MEF90Ctx_type),Intent(IN)          :: MEF90Ctx
        Type(tPetscSF),intent(OUT)              :: sf
        PetscErrorCode,intent(INOUT)            :: ierr
    
        Type(tPetscSection)                     :: locSection, gSection
        Type(tPetscSF)                          :: overlapSF, idSF
        Type(PetscSFNode),dimension(:),Pointer  :: remote
        PetscInt,dimension(:),Pointer           :: remoteOffsets
        PetscInt                                :: pStart, pEnd, p, n = 1
        PetscMPIInt                             :: rank
    
        PetscCallMPI(MPI_Comm_rank(MEF90Ctx%Comm, rank, ierr))
        PetscCallA(DMGetLocalSection(dm,locSection,ierr))
        PetscCallA(DMGetPointSF(dm,overlapSF,ierr))
        PetscCallA(PetscSectionCreateGlobalSection(locSection, overlapSF, PETSC_TRUE, PETSC_TRUE, gSection,ierr))
        PetscCallA(PetscSectionGetChart(locSection, pStart, pEnd, ierr))
        n = pEnd-pStart
        Allocate(remote(n))
        Do p = 1,n
            remote(p)%rank = rank
            remote(p)%index = p-1
        End Do
        PetscCallA(PetscSFCreate(MEF90Ctx%Comm, idSF, ierr))
        PetscCallA(PetscSFSetFromOptions(idSF, ierr))
        PetscCallA(PetscSFSetGraph(idSF, n, n, PETSC_NULL_INTEGER, PETSC_COPY_VALUES, remote, PETSC_COPY_VALUES, ierr))
        PetscCallA(PetscSFSetUp(idSF, ierr))
        PetscCallA(PetscSFCreateRemoteOffsetsF90(idSF,locSection,gSection,remoteOffsets,ierr))
        PetscCallA(PetscSFCreateSectionSFF90(idSF,locSection,remoteOffsets,gSection,sf,ierr))
        PetscCallA(PetscIntArray1dDestroyF90(remoteOffsets,ierr))
        PetscCallA(PetscSFSetUp(sf,ierr))
        PetscCallA(PetscObjectSetName(sf, "Local-To-CGlobal SF", ierr))
        PetscCallA(PetscSFSetFromOptions(sf, ierr))
        PetscCallA(PetscSFViewFromOptions(sf,PETSC_NULL_OPTIONS,"-localtocglobal_sf_view",ierr))
        PetscCallA(PetscSectionDestroy(gSection,ierr))
        PetscCallA(PetscSFDestroy(idSF,ierr))
        DeAllocate(remote)
    
    End subroutine CreateLocalToCGlobalSF_internal
    
#undef __FUNCT__
#define __FUNCT__ "CreateCGlobalToLocalSF_internal"
!!!
!!!  
!!!  CreateCGlobalToLocalSF_internal: 
!!!  
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!
    
    subroutine CreateCGlobalToLocalSF_internal(MEF90Ctx,dm,sf,ierr)
        Type(tDM),intent(IN)                    :: dm
        Type(MEF90Ctx_type),Intent(IN)          :: MEF90Ctx
        Type(tPetscSF),intent(OUT)              :: sf
        PetscErrorCode,intent(INOUT)            :: ierr
    
        Type(tPetscSection)                     :: locSection, gSection
        Type(tPetscSF)                          :: overlapSF, idSF, tempSF, ttempSF
        Type(PetscSFNode),dimension(:),Pointer  :: remote, tempRemote, lgRemote, glRemote
        PetscInt,dimension(:),Pointer           :: tempLocal, lgLocal, glLocal, remoteOffsets
        PetscInt                                :: pStart, pEnd, p, n = 1, lgNRoots, lgNLeaves, tempNRoots, tempNLeaves, glNRoots, glNLeaves
        PetscMPIInt                             :: rank
    
        PetscCallMPI(MPI_Comm_rank(MEF90Ctx%Comm, rank, ierr))
        PetscCallA(DMGetLocalSection(dm,locSection,ierr))
        PetscCallA(DMGetPointSF(dm,overlapSF,ierr))
        PetscCallA(PetscSectionCreateGlobalSection(locSection, overlapSF, PETSC_TRUE, PETSC_TRUE, gSection,ierr))
        PetscCallA(PetscSectionGetChart(locSection, pStart, pEnd, ierr))
        n = pEnd-pStart
        Allocate(remote(n))
        Do p = 1,n
            remote(p)%rank = rank
            remote(p)%index = p-1
        End Do
        PetscCallA(PetscSFCreate(MEF90Ctx%Comm, idSF, ierr))
        PetscCallA(PetscSFSetFromOptions(idSF, ierr))
        PetscCallA(PetscSFSetGraph(idSF, n, n, PETSC_NULL_INTEGER, PETSC_COPY_VALUES, remote, PETSC_COPY_VALUES, ierr))
        PetscCallA(PetscSFSetUp(idSF, ierr))
        PetscCallA(PetscSFCreateRemoteOffsetsF90(idSF,locSection,gSection,remoteOffsets,ierr))
        PetscCallA(PetscSFCreateSectionSFF90(idSF,locSection,remoteOffsets,gSection,tempSF,ierr))
        PetscCallA(PetscSFCreateInverseSF(tempSF,sf,ierr))
        PetscCallA(PetscIntArray1dDestroyF90(remoteOffsets,ierr))
        PetscCallA(PetscSFDestroy(tempSF,ierr))
        If (MEF90Ctx%NumProcs > 1) Then
            PetscCallA(PetscSFGetGraph(sf,lgNRoots,lgNLeaves,lgLocal,lgRemote,ierr))
            PetscCallA(PetscSFCreateRemoteOffsetsF90(overlapSF,locSection,gSection,remoteOffsets,ierr))
            PetscCallA(PetscSFCreateSectionSFF90(overlapSF,locSection,remoteOffsets,gSection,ttempSF,ierr))
            PetscCallA(PetscSFCreateInverseSF(ttempSF,tempSF,ierr))
            PetscCallA(PetscIntArray1dDestroyF90(remoteOffsets,ierr))
            PetscCallA(PetscSFDestroy(ttempSF,ierr))
            PetscCallA(PetscSFGetGraph(tempSF,tempNRoots,tempNLeaves,tempLocal,tempRemote,ierr))
            glNRoots = lgNRoots
            glNLeaves = lgNLeaves + tempNLeaves
            Allocate(glLocal(glNLeaves))
            Allocate(glRemote(glNLeaves))
            If (loc(lgLocal) .ne. loc(PETSC_NULL_INTEGER)) Then
                Do p = 1,lgNLeaves
                    glLocal(p) = lgLocal(p)
                    glRemote(p)%rank = lgRemote(p)%rank
                    glRemote(p)%index = lgRemote(p)%index
                End Do
            Else 
                Do p = 1,lgNLeaves
                    glLocal(p) = p-1
                    glRemote(p)%rank = lgRemote(p)%rank
                    glRemote(p)%index = lgRemote(p)%index
                End Do
            End If
            If (loc(tempLocal) .ne. loc(PETSC_NULL_INTEGER)) Then
                Do p = 1,tempNLeaves
                    glLocal(p+lgNLeaves) = tempLocal(p)
                    glRemote(p+lgNLeaves)%rank = tempRemote(p)%rank
                    glRemote(p+lgNLeaves)%index = tempRemote(p)%index
                End Do
            Else
                Do p = 1,tempNLeaves
                    glLocal(p+lgNLeaves) = p+lgNLeaves-1
                    glRemote(p+lgNLeaves)%rank = tempRemote(p)%rank
                    glRemote(p+lgNLeaves)%index = tempRemote(p)%index
                End Do
            End If
            PetscCallA(PetscSFDestroy(tempSF,ierr))
            PetscCallA(PetscSFDestroy(sf,ierr))
            PetscCallA(PetscSFCreate(MEF90Ctx%Comm, sf, ierr))
            PetscCallA(PetscSFSetFromOptions(sf, ierr))
            PetscCallA(PetscSFSetGraph(sf,glNRoots,glNLeaves,glLocal,PETSC_COPY_VALUES,glRemote,PETSC_COPY_VALUES,ierr))
            PetscCallA(PetscSFSetUp(sf,ierr))
            DeAllocate(glLocal)
            DeAllocate(glRemote)
        End If
        PetscCallA(PetscSectionDestroy(gSection,ierr))
        PetscCallA(PetscSFDestroy(idSF,ierr))
        PetscCallA(PetscObjectSetName(sf, "CGlobal-To-Local SF", ierr))
        PetscCallA(PetscSFViewFromOptions(sf,PETSC_NULL_OPTIONS,"-cglobaltolocal_sf_view",ierr))
        DeAllocate(remote)
    End subroutine CreateCGlobalToLocalSF_internal                    
End Module m_MEF90_DMPlex
