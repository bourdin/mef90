Module m_MEF90_DMPlex
#include "petsc/finclude/petsc.h"
    Use m_MEF90_Elements
    Use m_MEF90_Ctx
    Use petsc
    Use,intrinsic :: iso_c_binding
    IMPLICIT NONE
    

    Character(len=MEF90_MXSTRLEN),dimension(4),Parameter  :: MEF90_DMPlexSetLabelName = &
              ['Cell Sets  ','Face Sets  ','Edge Sets  ','Vertex Sets']
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
!!!  MEF90_CreateLocalToConstraintSF: sf mapping between local and constraint Vec ordering and distribution
!!!  
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!                Blaise Bourdin  bourdin@mcmaster.ca
!!!

    subroutine MEF90_CreateLocalToConstraintSF(MEF90Ctx,dm,dmB,sf,invSF,ierr)
        Type(tDM),intent(IN)                    :: dm, dmB
        Type(tPetscSF),intent(OUT)              :: sf, invSF
        Type(MEF90Ctx_type),Intent(IN)          :: MEF90Ctx
        PetscErrorCode,intent(INOUT)            :: ierr

        Type(tPetscSection)                     :: locSection, locBSection
        Type(PetscSFNode),dimension(:),Pointer  :: remote
        Type(PetscInt),dimension(:),Pointer     :: local, cindices
        PetscInt                                :: pStart, pEnd, p, d, nleaves, ldof, loff, cdof, coff, nsize, nroots
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
        PetscCall(PetscObjectSetName(sf, "Local-To-Constraint SF", ierr))
        !PetscCall(PetscSFViewFromOptions(sf,PETSC_NULL_OPTIONS,"-localtoconstraint_sf_view",ierr))
        PetscCall(PetscSFCreateInverseSF(sf,invSF,ierr))
        PetscCall(PetscSFSetUp(invSF, ierr))
        PetscCall(PetscObjectSetName(invSF, "Constraint-To-Local SF", ierr))
        !PetscCall(PetscSFViewFromOptions(invSF,PETSC_NULL_OPTIONS,"-constrainttolocal_sf_view",ierr))        
    End subroutine MEF90_CreateLocalToConstraintSF

#undef __FUNCT__
#define __FUNCT__ "MEF90_CreateLocalMEF90_VecGlobalToLocalConstraintToConstraintSF"
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
                    PetscCallA(VecGetValuesSectionF90(c,s,p,vArray,ierr))
                    PetscCallA(VecSetValuesSectionF90(l,s,p,vArray,INSERT_ALL_VALUES,ierr))
                    PetscCallA(VecRestoreValuesSectionF90(c,s,p,vArray,ierr))
                End If ! numConstraint
            End Do ! p
        End If
    End Subroutine MEF90_VecGlobalToLocalConstraint
End Module m_MEF90_DMPlex