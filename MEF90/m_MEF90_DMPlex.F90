Module m_MEF90_DMPlex
#include "petsc/finclude/petsc.h"
    Use m_MEF90_Elements
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
#define __FUNCT__ "MEF90_SectionAllocateDofAllSets"
!!!
!!!  
!!!  MEF90_SectionAllocateDofAllSets: Associates dof to a section
!!!  
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

Subroutine MEF90_SectionAllocateDofAllSets(dm,setType,elemType,numComponents,section,ierr)
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
End Subroutine MEF90_SectionAllocateDofAllSets

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
        PetscInt                           :: p,depth,sDim,closureSize
        PetscInt,dimension(:),pointer      :: pStartDepth,pEndDepth


        PetscCall(DMGetDimension(dm,sdim,ierr))
        Allocate(pStartDepth(sdim+1))
        Allocate(pEndDepth(sdim+1))
        Do depth = 1,sdim+1
            PetscCall(DMPlexGetDepthStratum(dm,depth-1,pStartDepth(depth),pEndDepth(depth),ierr))
        End Do
    
        PetscCall(DMGetStratumIS(dm,setType,setID,setPointIS,ierr))
        PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
        If (size(setPointID) > 0) Then
            !!! This can probably be optimize by allocating closure outside of the loop
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
#define __FUNCT__ "MEF90_SectionAllocateConstraintSet"
!!!
!!!  
!!!  MEF90_SectionAllocateConstraintSet: Allocate space for constraints in a section
!!!  
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!
    
    Subroutine MEF90_SectionAllocateConstraintSet(dm,setType,setID,numComponents,section,ierr)
        Type(tDM),Intent(IN)               :: dm
        Character(len=*),Intent(IN)        :: setType
        PetscInt,Intent(IN)                :: setID
        PetscInt,Intent(IN)                :: numComponents
        Type(tPetscSection),Intent(INOUT)  :: section
        PetscErrorCode,Intent(INOUT)       :: ierr

        Type(tIS)                          :: setPointIS
        PetscInt,dimension(:),Pointer      :: setPointID
        PetscInt                           :: point,numDof
        PetscInt,dimension(:),pointer      :: closure
        PetscInt                           :: p,depth,sDim,closureSize
        PetscInt,dimension(:),pointer      :: pStartDepth,pEndDepth


        PetscCall(DMGetDimension(dm,sdim,ierr))
        Allocate(pStartDepth(sdim+1))
        Allocate(pEndDepth(sdim+1))
        Do depth = 1,sdim+1
            PetscCall(DMPlexGetDepthStratum(dm,depth-1,pStartDepth(depth),pEndDepth(depth),ierr))
        End Do
    
        PetscCall(DMGetStratumIS(dm,setType,setID,setPointIS,ierr))
        PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
        If (size(setPointID) > 0) Then
            !!! This can probably be optimize by allocating closure outside of the loop
            !!! But I can't figure out how it is done at the moment.
            Nullify(closure)
            Do point = 1,size(setPointID)
                PetscCall(DMPlexGetTransitiveClosure(dm,setPointID(point),PETSC_TRUE,closure,ierr))
                Do p = 1,size(closure),2
                    PetscCall(PetscSectionGetDoF(section,closure(p),numDof,ierr))
                    If (numDof > 0) Then
                        PetscCall(PetscSectionSetConstraintDof(section,closure(p),numComponents,ierr))
                    End If
                End Do! p
                PetscCall(DMPlexRestoreTransitiveClosure(dm,setPointID(point),PETSC_TRUE,closure,ierr))
            End Do! cell
        End If
        PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
        PetscCall(ISDestroy(setPointIS,ierr))
    End Subroutine MEF90_SectionAllocateConstraintSet

#undef __FUNCT__
#define __FUNCT__ "MEF90_SectionSetConstraintIndicesSet"
!!!
!!!  
!!!  MEF90_SectionSetConstraintIndicesSet: Set constrained points in a section
!!!  
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!
        
    Subroutine MEF90_SectionSetConstraintIndicesSet(dm,setType,setID,ConstraintIndices,section,ierr)
        Type(tDM),Intent(IN)               :: dm
        Character(len=*),Intent(IN)        :: setType
        PetscInt,Intent(IN)                :: setID
        PetscInt,Dimension(:),Pointer      :: ConstraintIndices
        Type(tPetscSection),Intent(INOUT)  :: section
        PetscErrorCode,Intent(INOUT)       :: ierr

        Type(tIS)                          :: setPointIS
        PetscInt,dimension(:),Pointer      :: setPointID
        PetscInt                           :: point,numDof
        PetscInt,dimension(:),pointer      :: closure
        PetscInt                           :: p,depth,sDim,closureSize
        PetscInt,dimension(:),pointer      :: pStartDepth,pEndDepth


        PetscCall(DMGetDimension(dm,sdim,ierr))
        Allocate(pStartDepth(sdim+1))
        Allocate(pEndDepth(sdim+1))
        Do depth = 1,sdim+1
            PetscCall(DMPlexGetDepthStratum(dm,depth-1,pStartDepth(depth),pEndDepth(depth),ierr))
        End Do
    
        PetscCall(DMGetStratumIS(dm,setType,setID,setPointIS,ierr))
        PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
        If (size(setPointID) > 0) Then
            !!! This can probably be optimize by allocating closure outside of the loop
            !!! But I can't figure out how it is done at the moment.
            Nullify(closure)
            Do point = 1,size(setPointID)
                PetscCall(DMPlexGetTransitiveClosure(dm,setPointID(point),PETSC_TRUE,closure,ierr))
                Do p = 1,size(closure),2
                    PetscCall(PetscSectionGetDoF(section,closure(p),numDof,ierr))
                    If (numDof > 0) Then
                        PetscCall(PetscSectionSetConstraintIndicesF90(section,closure(p),ConstraintIndices,ierr))
                    End If
                End Do! p
                PetscCall(DMPlexRestoreTransitiveClosure(dm,setPointID(point),PETSC_TRUE,closure,ierr))
            End Do! cell
        End If
        PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
        PetscCall(ISDestroy(setPointIS,ierr))
    End Subroutine MEF90_SectionSetConstraintIndicesSet    
End Module m_MEF90_DMPlex