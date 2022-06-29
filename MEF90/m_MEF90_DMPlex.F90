Module m_MEF90_DMPlex
#include "petsc/finclude/petsc.h"
    Use m_MEF90_Elements
    Use m_MEF90_Ctx
    Use petsc
    IMPLICIT NONE
    

    Character(len=MEF90_MXSTRLEN),dimension(4),Parameter  :: MEF90_DMPlexSetTypes = &
              ['Cell Sets  ','Face Sets  ','Edge Sets  ','Vertex Sets']
    Character(len=MEF90_MXSTRLEN),Parameter :: MEF90_DMPlexCellSetType   = 'Cell Sets  '
    Character(len=MEF90_MXSTRLEN),Parameter :: MEF90_DMPlexFaceSetType   = 'Face Sets  '
    Character(len=MEF90_MXSTRLEN),Parameter :: MEF90_DMPlexEdgeSetType   = 'Edge Sets  '
    Character(len=MEF90_MXSTRLEN),Parameter :: MEF90_DMPlexVertexSetType = 'Vertex Sets'

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
    Character(len=*),Intent(IN)        :: setType
    Type(MEF90Element_Type),Intent(IN) :: elemType
    PetscInt,Intent(IN)                :: numComponents
    Type(tPetscSection),Intent(INOUT)  :: section
    PetscErrorCode,Intent(INOUT)       :: ierr

    Type(tIS)                          :: setIS
    PetscInt,dimension(:),Pointer      :: setID
    Type(MEF90Element_Type)            :: CSelemType
    PetscInt                           :: set

    PetscCall(DMGetLabelIdIS(dm,setType,setIS,ierr))
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
        Character(len=*),Intent(IN)        :: setType
        PetscInt,Intent(IN)                :: setID
        Type(MEF90Element_Type),Intent(IN) :: elemType
        PetscInt,Intent(IN)                :: numComponents
        Type(tPetscSection),Intent(INOUT)  :: section
        PetscErrorCode,Intent(INOUT)       :: ierr

        Type(tIS)                          :: setPointIS
        PetscInt,dimension(:),Pointer      :: setPointID
        PetscInt                           :: point
        PetscInt,dimension(:),pointer      :: closure
        PetscInt                           :: p,depth,closureSize

        PetscCall(DMGetStratumIS(dm,setType,setID,setPointIS,ierr))
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
        Character(len=*),Intent(IN)        :: setType
        PetscInt,Intent(IN)                :: setID
        Logical,Dimension(:),Pointer       :: constraints
        Logical,Dimension(:,:),Pointer     :: table
        PetscErrorCode,Intent(INOUT)       :: ierr

        Type(tIS)                          :: setPointIS
        PetscInt,dimension(:),Pointer      :: setPointID
        PetscInt                           :: point,numDof
        PetscInt,dimension(:),pointer      :: closure
        PetscInt                           :: p,c

        PetscCall(DMGetStratumIS(dm,setType,setID,setPointIS,ierr))
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

    PetscInt                           :: p,c,i,pStart,pEnd,numConstraints,numComponents
    PetscInt,Dimension(:),Pointer      :: constraints

    PetscCall(DMPlexGetChart(dm,pStart,pEnd,ierr))
    Do p = 1, pEnd
        numConstraints = count(table(p,:))
        If (numConstraints > 0) Then
            PetscCall(PetscSectionSetConstraintDof(section,p-1,numConstraints,ierr))
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
            DeAllocate(constraints)
        End If
    End Do
End Subroutine MEF90_SectionAllocateConstraint

#undef __FUNCT__
#define __FUNCT__ "MEF90_CreateLocalToConstraintSF"
!!!
!!!  
!!!  MEF90_CreateLocalToConstraintSF: sf mapping between local and constraint Vec odering and distribution
!!!  
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!

subroutine MEF90_CreateLocalToConstraintSF(MEF90Ctx,dm,dmB,sf,invSF,ierr)
    Type(tDM),intent(IN)                    :: dm, dmB
    Type(tPetscSF),intent(OUT)              :: sf, invSF
    Type(MEF90Ctx_type),Intent(IN)          :: MEF90Ctx
    PetscErrorCode,intent(INOUT)            :: ierr

    Type(tPetscSection)                     :: locSection, locBSection, gBSection
    Type(PetscSFNode),dimension(:),Pointer  :: remote
    Type(PetscInt),dimension(:),Pointer     :: local, cindices
    PetscInt                                :: pStart, pEnd, p, d, nleaves = 0, ldof, loff, cdof, coff, nsize = 0, nroots
    PetscMPIInt                             :: rank

    PetscCallMPI(MPI_Comm_rank(MEF90Ctx%Comm, rank, ierr))
    PetscCallA(DMGetLocalSection(dm,locSection,ierr))
    PetscCallA(PetscSectionGetStorageSize(locSection,nroots,ierr))
    PetscCallA(DMGetLocalSection(dmB,locBSection,ierr))
    PetscCallA(PetscSectionGetChart(locBSection, pStart, pEnd, ierr))
    Do p = pStart,pEnd-1
        PetscCallA(PetscSectionGetConstraintDof(locBSection,p,cdof,ierr))
        nleaves = nleaves + cdof
    End Do
    Allocate(remote(nleaves))
    Allocate(local(nleaves))
    Do p = pStart,pEnd-1
        PetscCallA(PetscSectionGetDoF(locSection,p,ldof,ierr))
        PetscCallA(PetscSectionGetOffset(locSection,p,loff,ierr))
        PetscCallA(PetscSectionGetConstraintDof(locBSection,p,cdof,ierr))
        PetscCallA(PetscSectionGetConstraintIndicesF90(locBSection,p,cindices,ierr))
        PetscCallA(PetscSectionGetOffset(locBSection,p,coff,ierr))
        If (coff >= 0) Then
            Do d=1,cdof
                local(nsize+1) = coff+cindices(d)
                remote(nsize+1)%rank = rank
                remote(nsize+1)%index = loff+cindices(d)
                nsize = nsize + 1
            End Do
        End If
        PetscCallA(PetscSectionRestoreConstraintIndicesF90(locBSection,p,cindices,ierr))
    End Do
    PetscCallA(PetscSFCreate(MEF90Ctx%Comm, sf, ierr))
    PetscCallA(PetscSFSetFromOptions(sf, ierr))
    PetscCallA(PetscSFSetGraph(sf, nroots, nleaves, local, PETSC_COPY_VALUES, remote, PETSC_COPY_VALUES, ierr))
    DeAllocate(remote)
    DeAllocate(local)
    PetscCallA(PetscSFSetUp(sf, ierr))
    PetscCallA(PetscObjectSetName(sf, "Local-To-Constraint SF", ierr))
    PetscCallA(PetscSFViewFromOptions(sf,PETSC_NULL_OPTIONS,"-localtoconstraint_sf_view",ierr))
    PetscCallA(PetscSFCreateInverseSF(sf,invSF,ierr))
    PetscCallA(PetscSFSetUp(invSF, ierr))
    PetscCallA(PetscObjectSetName(invSF, "Constraint-To-Local SF", ierr))
    PetscCallA(PetscSFViewFromOptions(invSF,PETSC_NULL_OPTIONS,"-constrainttolocal_sf_view",ierr))
    
End subroutine MEF90_CreateLocalToConstraintSF
End Module m_MEF90_DMPlex