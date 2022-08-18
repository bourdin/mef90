Module m_MEF90_DMPlex
#include "petsc/finclude/petsc.h"
    Use m_MEF90_Elements
    Use m_MEF90_Ctx
    Use petsc
    Use,intrinsic :: iso_c_binding
    IMPLICIT NONE
    
    Enum,bind(c)
        enumerator  :: MEF90CellSetType = 1, &
                       MEF90FaceSetType,     &
                       !MEF90EdgeSetType,     &
                       MEF90VertexSetType
    End Enum
    PetscEnum,Dimension(3),Parameter :: MEF90SetType = [ &
        MEF90CellSetType,                                &
        MEF90FaceSetType,                                &
        !MEF90EdgeSetType,                                &
        MEF90VertexSetType ]

    Character(len=MEF90MXSTRLEN),Parameter :: MEF90CellSetLabelName   = 'Cell Sets  '
    Character(len=MEF90MXSTRLEN),Parameter :: MEF90FaceSetLabelName   = 'Face Sets  '
    !Character(len=MEF90MXSTRLEN),Parameter :: MEF90EdgeSetLabelName   = 'Edge Sets  '
    Character(len=MEF90MXSTRLEN),Parameter :: MEF90VertexSetLabelName = 'Vertex Sets'

    Character(len=MEF90MXSTRLEN),dimension(3),Parameter  :: MEF90SetLabelName = &
              [MEF90CellSetLabelName, &
               MEF90FaceSetLabelName, &
               !MEF90EdgeSetLabelName, &
               MEF90VertexSetLabelName]

    Character(len=MEF90MXSTRLEN),Parameter :: MEF90CellSetprefix   = 'cs'
    Character(len=MEF90MXSTRLEN),Parameter :: MEF90FaceSetprefix   = 'fs'
    !Character(len=MEF90MXSTRLEN),Parameter :: MEF90EdgeSetprefix   = 'es'
    Character(len=MEF90MXSTRLEN),Parameter :: MEF90VertexSetprefix = 'vs'

    Character(len=MEF90MXSTRLEN),dimension(3),Parameter  :: MEF90SetPrefix = &
                [MEF90CellSetPrefix, &
                MEF90FaceSetPrefix, &
                !MEF90EdgeSetPrefix, &
                MEF90VertexSetPrefix]

    private
    public :: MEF90CellSetLabelName,                                                 &
              MEF90FaceSetLabelName,                                                 &
              !MEF90EdgeSetLabelName,                                                 &
              MEF90VertexSetLabelName,                                               &
              MEF90SetLabelName,                                                     &
              MEF90CellSetPrefix,                                                    &
              MEF90FaceSetPrefix,                                                    &
              !MEF90EdgeSetPrefix,                                                    &
              MEF90VertexSetPrefix,                                                  &
              MEF90SetPrefix,                                                        &
              MEF90CellSetType,                                                      &
              MEF90FaceSetType,                                                      &
              !MEF90EdgeSetType,                                                      &
              MEF90VertexSetType,                                                    &
              MEF90SectionAllocateDof,MEF90SectionAllocateDofSet,                    &
              MEF90SetupConstraintTableSet,MEF90SectionAllocateConstraint,           &
              MEF90CellSectionCreate,                                                &
              MEF90VecCopySF,MEF90IOSFCreate,                                        &
              MEF90ConstraintSFCreate,MEF90VecGlobalToLocalConstraint,               &
              MEF90VecCreateIO,                                                      &
              MEF90VecCreate
Contains

#undef __FUNCT__
#define __FUNCT__ "MEF90VecCreate"
!!!
!!!  
!!!  MEF90VecCreate: create a Vec associated with a FE space and constraints
!!!      cell   set BC are obtained from the command line option -cs<set ID>_<name>BC [bool], [bool], ...
!!!      face   set BC are obtained from the command line option -fs<set ID>_<name>BC [bool], [bool], ...
!!!      vertex set BC are obtained from the command line option -vs<set ID>_<name>BC [bool], [bool], ...
!!!  
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!
    Subroutine MEF90VecCreate(dm,elemFamily,elemOrder,sDim,name,V,ierr)
        Type(tDM),Intent(IN)                    :: dm
        PetscEnum,Intent(IN)                    :: elemFamily
        PetscInt,Intent(IN)                     :: elemOrder,sDim
        Character(len=MEF90MXSTRLEN),Intent(IN) :: name
        Type(tVec),Intent(OUT)                  :: V
        PetscErrorCode,Intent(INOUT)            :: ierr

        Type(tPetscSection)                     :: sectionV
        Type(tDM)                               :: dmV
        PetscInt,Dimension(1)                   :: fieldV = 0
        PetscInt                                :: pStart,pEnd
        MPI_Comm                                :: comm
        PetscInt                                :: set
        PetscEnum                               :: setType
        PetscInt,Dimension(:),pointer           :: setID,pointID
        Type(tIS)                               :: setIS,pointIS
        Type(MEF90ElementType)                  :: elemType
        DMPolytopeType                          :: cellType
        PetscBool,Dimension(:,:),Pointer        :: constraintTruthTable
        PetscBool,Dimension(:),Pointer          :: setConstraints
        PetscInt                                :: numBC
        PetscBool                               :: flg
        Character(len=MEF90MXSTRLEN)            :: BCOptionName

        PetscCall(DMClone(dm,dmV,ierr))
        PetscCall(PetscObjectSetName(dmv,name,ierr))
        PetscCall(DMGetUseNatural(dm,flg,ierr))
        PetscCall(DMSetUseNatural(dmV,flg,ierr))

        PetscCall(PetscObjectGetComm(dmV,comm,ierr))

        PetscCall(PetscSectionCreate(comm,sectionV,ierr))
        PetscCall(PetscObjectSetName(sectionV,name,ierr))
        PetscCall(PetscSectionSetNumFields(sectionV,sDim,ierr))
        PetscCall(PetscSectionSetFieldName(sectionV,fieldV,trim(name),ierr))
        PetscCall(PetscSectionSetFieldComponents(sectionV,fieldV,sdim,ierr))
        PetscCall(DMPlexGetChart(dmV,pStart,pEnd,ierr))
        PetscCall(PetscSectionSetChart(sectionV,pStart,pEnd,ierr))
    
        PetscCall(DMGetLabelIdIS(dmV,MEF90CellSetLabelName,setIS,ierr))
        !!! Get a GLOBAL cell set IS
        ! PetscCall(MEF90ISAllGatherMerge(comm,setIS,ierr))
        If (setIS /= PETSC_NULL_IS) Then
            PetscCall(ISGetIndicesF90(setIS,setID,ierr))
            Do set = 1,size(setID)
                !!! Get cell type in order to pick the proper element type.
                !!! We assume that all cells in a set have the same type, so all we need it to query the first cell in the set
                PetscCall(DMGetStratumIS(dmV,MEF90CellSetLabelName,setID(set),pointIS,ierr))
                If (pointIS /= PETSC_NULL_IS) Then
                    PetscCall(ISGetIndicesF90(pointIS,pointID,ierr))
                    PetscCall(DMPlexGetCellType(dmV,pointID(1),cellType,ierr))
                    PetscCall(MEF90ElementGetType(elemFamily,elemOrder,cellType,elemType,ierr))
                    PetscCall(MEF90SectionAllocateDofSet(dmV,MEF90CellSetType,setID(set),elemType,sdim,sectionV,ierr))
                    PetscCall(ISRestoreIndicesF90(pointIS,pointID,ierr))
                End If ! pointIS
                PetscCall(ISDestroy(pointIS,ierr))
            End Do ! set
            PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
        End If ! setIS
        PetscCall(ISDestroy(setIS,ierr))

!!! removing DOF allocation at face sets for now.
!!! Instead, I will create Boundary Vecs
        
        ! PetscCall(DMGetLabelIdIS(dmV,MEF90FaceSetLabelName,setIS,ierr))
        ! !!! Get a GLOBAL face set IS
        ! ! PetscCall(MEF90ISAllGatherMerge(comm,setIS,ierr))
        ! If (setIS /= PETSC_NULL_IS) Then
        !     PetscCall(ISGetIndicesF90(setIS,setID,ierr))
        !     Do set = 1,size(setID)
        !         !!! Get cell type in order to pick the proper element type.
        !         !!! We assume that all cells in a set have the same type, so all we need it to query the first cell in the set
        !         PetscCall(DMGetStratumIS(dmV,MEF90FaceSetLabelName,setID(set),pointIS,ierr))
        !         If (pointIS /= PETSC_NULL_IS) Then
        !             PetscCall(ISGetIndicesF90(pointIS,pointID,ierr))
        !             PetscCall(DMPlexGetCellType(dmV,pointID(1),cellType,ierr))
        !             PetscCall(MEF90ElementGetTypeBoundary(elemFamily,elemOrder,cellType,elemType,ierr))
        !             PetscCall(MEF90SectionAllocateDofSet(dmV,MEF90FaceSetType,setID(set),elemType,sdim,sectionV,ierr))
        !             PetscCall(ISRestoreIndicesF90(pointIS,pointID,ierr))
        !         End If ! pointIS
        !         PetscCall(ISDestroy(pointIS,ierr))
        !     End Do ! set
        !     PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
        ! End If ! setIS
        ! PetscCall(ISDestroy(setIS,ierr))

        Allocate(ConstraintTruthTable(pEnd,sDim))
        ConstraintTruthTable = .FALSE.
        Allocate(setConstraints(sDim))

        Do setType = 1,size(MEF90SetType)
            PetscCall(DMGetLabelIdIS(dm,MEF90SetLabelName(setType),setIS,ierr))
            ! PetscCall(MEF90ISAllGatherMerge(comm,setIS,ierr))
            If (setIS /= PETSC_NULL_IS) Then
                PetscCall(ISGetIndicesF90(setIS,setID,ierr))
                Do set = 1,size(setID)
                    setConstraints = .FALSE.
                    write(BCOptionName,'("-",a2,I4.4,"_",a,"BC")') MEF90SetPrefix(setType),setID(set),trim(name)
                    numBC = sDim
                    PetscCall(PetscOptionsGetBoolArray(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,trim(BCOptionName),setConstraints,numBC,flg,ierr))
                    PetscCall(MEF90SetupConstraintTableSet(dmV,sectionV,MEF90SetType(setType),setID(set),setConstraints,ConstraintTruthTable,ierr))
                End Do
                PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
            End If ! setIS
            PetscCall(ISDestroy(setIS,ierr))
         End Do ! setType

        PetscCall(MEF90SectionAllocateConstraint(dmV,ConstraintTruthTable,sectionV,ierr))
        DeAllocate(ConstraintTruthTable)
        DeAllocate(setConstraints)
    
        PetscCall(DMSetLocalSection(dmV,sectionV,ierr))
        PetscCall(DMCreateLocalVector(dmV,V,ierr))
        PetscCall(PetscObjectSetName(V,name,ierr))
    End Subroutine MEF90VecCreate

#undef __FUNCT__
#define __FUNCT__ "MEF90SectionAllocateDof"
!!!
!!!  
!!!  MEF90SectionAllocateDof: Associates dof to a section
!!!  
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

    Subroutine MEF90SectionAllocateDof(dm,setType,elemType,numComponents,section,ierr)
        Type(tDM),Intent(IN)               :: dm
        PetscEnum,intent(IN)               :: setType
        Type(MEF90ElementType),Intent(IN)  :: elemType
        PetscInt,Intent(IN)                :: numComponents
        Type(tPetscSection),Intent(INOUT)  :: section
        PetscErrorCode,Intent(INOUT)       :: ierr

        Type(tIS)                          :: setIS
        PetscInt,dimension(:),Pointer      :: setID
        PetscInt                           :: set

        PetscCall(DMGetLabelIdIS(dm,MEF90SetLabelName(setType),setIS,ierr))
        If (setIS /= PETSC_NULL_IS) Then
            PetscCall(ISGetIndicesF90(setIS,setID,ierr))
            Do set = 1,size(setID)
                PetscCall(MEF90SectionAllocateDofSet(dm,setType,setID(set),elemType,numComponents,Section,ierr))
            End Do
            PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
        End If ! setIS
        PetscCall(ISDestroy(setIS,ierr))
    End Subroutine MEF90SectionAllocateDof

#undef __FUNCT__
#define __FUNCT__ "MEF90_SectionAllocateDofSet"
!!!
!!!  
!!!  MEF90SectionAllocateDofSet: Associates dof to a section in a set
!!!  
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

    Subroutine MEF90SectionAllocateDofSet(dm,setType,setID,elemType,numComponents,section,ierr)
        Type(tDM),Intent(IN)               :: dm
        PetscEnum,intent(IN)               :: setType
        PetscInt,Intent(IN)                :: setID
        Type(MEF90ElementType),Intent(IN)  :: elemType
        PetscInt,Intent(IN)                :: numComponents
        Type(tPetscSection),Intent(INOUT)  :: section
        PetscErrorCode,Intent(INOUT)       :: ierr

        Type(tIS)                          :: setPointIS
        PetscInt,dimension(:),Pointer      :: setPointID
        PetscInt                           :: point
        PetscInt,dimension(:),pointer      :: closure
        PetscInt                           :: p,depth
        PetscInt                           :: field = 0


        PetscCall(DMGetStratumIS(dm,MEF90SetLabelName(setType),setID,setPointIS,ierr))
        If (setPointIS /= PETSC_NULL_IS) Then
            PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
            !!! This can probably be optimized by allocating closure outside of the loop
            !!! But I can't figure out how it is done at the moment.
            Nullify(closure)
            Do point = 1,size(setPointID)
                PetscCall(DMPlexGetTransitiveClosure(dm,setPointID(point),PETSC_TRUE,closure,ierr))
                Do p = 1,size(closure),2
                    PetscCall(DMPlexGetPointDepth(dm,closure(p),depth,ierr))
                    If (elemType%numDofs(depth+1) > 0) Then
                        PetscCall(PetscSectionSetDof(section,closure(p),elemType%numDofs(depth+1)*numComponents,ierr))
                        PetscCall(PetscSectionSetFieldDof(section,closure(p),field,elemType%numDofs(depth+1)*numComponents,ierr))
                    End If
                End Do! p
                PetscCall(DMPlexRestoreTransitiveClosure(dm,setPointID(point),PETSC_TRUE,closure,ierr))
            End Do! cell
            PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
        End If ! setPointIS
    PetscCall(ISDestroy(setPointIS,ierr))
    End Subroutine MEF90SectionAllocateDofSet

#undef __FUNCT__
#define __FUNCT__ "MEF90CellSectionCreate"
!!!
!!!  
!!!  MEF90CellSectionCreate: create a section with numComponent dof at cells
!!!  
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!

    Subroutine MEF90CellSectionCreate(dm,numComponents,section,ierr)
        Type(tDM),Intent(IN)               :: dm
        PetscInt,Intent(IN)                :: numComponents
        Type(tPetscSection),Intent(INOUT)  :: section
        PetscErrorCode,Intent(INOUT)       :: ierr

        Type(tIS)                          :: setIS,setPointIS
        PetscInt,dimension(:),Pointer      :: setID,setPointID
        PetscInt                           :: set,point
        PetscInt,dimension(:),pointer      :: closure
        PetscInt                           :: p,depth,dim

        PetscCallA(DMGetDimension(dm,dim,ierr))
        PetscCall(DMGetLabelIdIS(dm,'Cell Sets',setIS,ierr))
        If (setIS /= PETSC_NULL_IS) Then
            PetscCall(ISGetIndicesF90(setIS,setID,ierr))
            Do set = 1,size(setID)
                PetscCall(DMGetStratumIS(dm,'Cell Sets',setID(set),setPointIS,ierr))
                If (setPointIS /= PETSC_NULL_IS) Then
                    PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
                    !!! This can probably be optimized by allocating closure outside of the loop
                    !!! But I can't figure out how it is done at the moment.
                    Nullify(closure)
                    Do point = 1,size(setPointID)
                        PetscCall(DMPlexGetTransitiveClosure(dm,setPointID(point),PETSC_TRUE,closure,ierr))
                        Do p = 1,size(closure),2
                            PetscCall(DMPlexGetPointDepth(dm,closure(p),depth,ierr))
                            If (depth == dim) Then
                                PetscCall(PetscSectionSetDof(section,closure(p),numComponents,ierr))
                            End If
                        End Do! p
                        PetscCall(DMPlexRestoreTransitiveClosure(dm,setPointID(point),PETSC_TRUE,closure,ierr))
                    End Do! cell
                    PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
                End If ! setPointIS
            PetscCall(ISDestroy(setPointIS,ierr))
            End Do ! set
            PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
        End If ! setIS
        PetscCall(ISDestroy(setIS,ierr))
        PetscCall(PetscSectionSetup(section,ierr))
    End Subroutine MEF90CellSectionCreate

#undef __FUNCT__
#define __FUNCT__ "MEF90SetupConstraintTableSet"
!!!
!!!  
!!!  MEF90SetupConstraintTableSet: Build the contribution of a set to the constraint table
!!!  
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!
    
    Subroutine MEF90SetupConstraintTableSet(dm,section,setType,setID,constraints,table,ierr)
        Type(tDM),Intent(IN)               :: dm
        Type(tPetscSection),Intent(INOUT)  :: section
        PetscEnum,intent(IN)               :: setType
        PetscInt,Intent(IN)                :: setID
        PetscBool,Dimension(:),Pointer     :: constraints
        PetscBool,Dimension(:,:),Pointer   :: table
        PetscErrorCode,Intent(INOUT)       :: ierr

        Type(tIS)                          :: setPointIS
        PetscInt,dimension(:),Pointer      :: setPointID
        PetscInt                           :: point,numDof
        PetscInt,dimension(:),pointer      :: closure
        PetscInt                           :: p

        PetscCall(DMGetStratumIS(dm,MEF90SetLabelName(setType),setID,setPointIS,ierr))
        If (setPointIS /= PETSC_NULL_IS) Then
            PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
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
            PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
        End If ! setPointIS
        PetscCall(ISDestroy(setPointIS,ierr))
    End Subroutine MEF90SetupConstraintTableSet

#undef __FUNCT__
#define __FUNCT__ "MEF90SectionAllocateConstraint"
!!!
!!!  
!!!  MEF90SectionAllocateConstraint: Associates dof to a section
!!!  
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

    Subroutine MEF90SectionAllocateConstraint(dm,table,section,ierr)
        Type(tDM),Intent(IN)               :: dm
        Logical,Dimension(:,:),Pointer     :: table
        Type(tPetscSection),Intent(INOUT)  :: section
        PetscErrorCode,Intent(INOUT)       :: ierr

        PetscInt                           :: p,i,pStart,pEnd,numConstraints,numComponents
        PetscInt,Dimension(:),Pointer      :: constraints
        PetscInt                           :: field = 0

        PetscCall(DMPlexGetChart(dm,pStart,pEnd,ierr))
        Do p = 1,pEnd
            numConstraints = count(table(p,:))
            If (numConstraints > 0) Then
                PetscCall(PetscSectionSetConstraintDof(section,p-1,numConstraints,ierr))
                PetscCall(PetscSectionSetFieldConstraintDof(section,p-1,field,numConstraints,ierr))
            End If
        End Do

        PetscCall(PetscSectionSetup(section,ierr))

        numComponents = size(table,2)
        Do p = 1,pEnd
            numConstraints = count(table(p,:))
            If (numConstraints > 0) Then
                Allocate(constraints(numConstraints))
                constraints = pack([ (i-1,i = 1,numComponents) ],table(p,:))
                PetscCall(PetscSectionSetConstraintIndicesF90(section,p-1,constraints,ierr))
                PetscCall(PetscSectionSetFieldConstraintIndicesF90(section,p-1,field,constraints,ierr))
                DeAllocate(constraints)
            End If
        End Do
    End Subroutine MEF90SectionAllocateConstraint

#undef __FUNCT__
#define __FUNCT__ "MEF90VecCopySF"
!!!
!!!  
!!!  MEF90VecCopySF: rearrange a Vec according to the given SF
!!!  
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!

    Subroutine MEF90VecCopySF(vin,vout,sf,ierr)
        Type(tVec),intent(IN)              :: vin
        Type(tVec),intent(INOUT)           :: vout
        Type(tPetscSF),intent(IN)          :: sf
        PetscErrorCode,intent(INOUT)       :: ierr

        PetscScalar,dimension(:),Pointer   :: arrayin,arrayout

        PetscCall(VecGetArrayReadF90(vin,arrayin,ierr))
        PetscCall(VecGetArrayF90(vout,arrayout,ierr))
        PetscCall(PetscSFBcastBegin(sf,MPIU_SCALAR,arrayin,arrayout,MPI_REPLACE,ierr))
        PetscCall(PetscSFBcastEnd(sf,MPIU_SCALAR,arrayin,arrayout,MPI_REPLACE,ierr))
        PetscCall(VecRestoreArrayReadF90(vin,arrayin,ierr))
        PetscCall(VecRestoreArrayF90(vout,arrayout,ierr))
    End subroutine MEF90VecCopySF

#undef __FUNCT__
#define __FUNCT__ "MEF90IOSFCreate"
!!!
!!!  
!!!  MEF90IOSFCreate: sf mapping between local and IO ordering and distribution
!!!  
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!

    Subroutine MEF90IOSFCreate(MEF90Ctx,v,liosf,iolsf,ierr)
        Type(tVec),intent(IN)              :: v
        Type(tPetscSF),intent(OUT)         :: liosf,iolsf
        Type(MEF90Ctx_type),Intent(IN)     :: MEF90Ctx
        PetscErrorCode,intent(INOUT)       :: ierr

        Type(tPetscSF)                     :: naturalSF,ioSF,lcgSF,cglSF,tempSF,invTempSF
        Type(tDM)                          :: dm

        PetscCallA(VecGetDM(v,dm,ierr))
        PetscCall(CreateLocalToCGlobalSF_Private(MEF90Ctx,dm,lcgSF,ierr))
        PetscCall(CreateCGlobalToLocalSF_Private(MEF90Ctx,dm,cglSF,ierr))
        If (MEF90Ctx%NumProcs > 1) Then
            PetscCall(CreateNaturalToIOSF_Private(MEF90Ctx,dm,ioSF,ierr))
            PetscCall(DMGetNaturalSF(dm,naturalSF,ierr))
            PetscCall(PetscSFCompose(lcgSF,naturalSF,tempSF,ierr))
            PetscCall(PetscSFCompose(tempSF,ioSF,liosf,ierr))
            PetscCall(PetscSFSetUp(liosf,ierr))
            PetscCall(PetscSFDestroy(tempSF,ierr))
            PetscCall(PetscSFCompose(naturalSF,ioSF,tempSF,ierr))
            PetscCall(PetscSFCreateInverseSF(tempSF,invTempSF,ierr))
            PetscCall(PetscSFCompose(invTempSF,cglSF,iolsf,ierr))
            PetscCall(PetscSFSetUp(iolsf,ierr))
            PetscCall(PetscSFDestroy(ioSF,ierr))
            PetscCall(PetscSFDestroy(cglSF,ierr))
            PetscCall(PetscSFDestroy(lcgSF,ierr))
            PetscCall(PetscSFDestroy(tempSF,ierr))
            PetscCall(PetscSFDestroy(invTempSF,ierr))
        Else
            liosf = lcgSF
            iolsf = cglSF
        End If
    End subroutine MEF90IOSFCreate

#undef __FUNCT__
#define __FUNCT__ "MEF90ConstraintSFCreate"
!!!
!!!  
!!!  MEF90ConstraintSFCreate: sf mapping between local and constraint Vec ordering and distribution
!!!  
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!                Blaise Bourdin  bourdin@mcmaster.ca
!!!

    Subroutine MEF90ConstraintSFCreate(MEF90Ctx,v,vB,sf,invSF,ierr)
        Type(tVec),intent(IN)                   :: v,vB
        Type(tPetscSF),intent(OUT)              :: sf,invSF
        Type(MEF90Ctx_type),Intent(IN)          :: MEF90Ctx
        PetscErrorCode,intent(INOUT)            :: ierr

        Type(tDM)                               :: dm,dmB
        Type(tPetscSection)                     :: locSection,locBSection
        Type(PetscSFNode),dimension(:),Pointer  :: remote
        Type(PetscInt),dimension(:),Pointer     :: local,cindices
        PetscInt                                :: pStart,pEnd,p,d,nleaves = 0,ldof,loff,cdof,coff,nsize = 0,nroots

        nleaves = 0
        nsize   = 0
        PetscCall(VecGetDM(v,dm,ierr))
        PetscCall(VecGetDM(vB,dmB,ierr))
        PetscCall(DMGetLocalSection(dm,locSection,ierr))
        PetscCall(PetscSectionGetStorageSize(locSection,nroots,ierr))
        PetscCall(DMGetLocalSection(dmB,locBSection,ierr))
        PetscCall(PetscSectionGetChart(locBSection,pStart,pEnd,ierr))
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
            If (cdof > 0) Then
                PetscCall(PetscSectionGetConstraintIndicesF90(locBSection,p,cindices,ierr))
                PetscCall(PetscSectionGetOffset(locBSection,p,coff,ierr))
                If (coff >= 0) Then
                    Do d=1,cdof
                        local(nsize+1)        = coff+cindices(d)
                        remote(nsize+1)%rank  = MEF90Ctx%rank
                        remote(nsize+1)%index = loff+cindices(d)
                        nsize = nsize + 1
                    End Do ! d
                End If ! coff
                PetscCall(PetscSectionRestoreConstraintIndicesF90(locBSection,p,cindices,ierr))
            End If ! cdof
        End Do
        PetscCall(PetscSFCreate(MEF90Ctx%Comm,sf,ierr))
        PetscCall(PetscSFSetFromOptions(sf,ierr))
        PetscCall(PetscSFSetGraph(sf,nroots,nleaves,local,PETSC_COPY_VALUES,remote,PETSC_COPY_VALUES,ierr))
        DeAllocate(remote)
        DeAllocate(local)
        PetscCall(PetscSFSetUp(sf,ierr))
        PetscCall(PetscSFCreateInverseSF(sf,invSF,ierr))
        PetscCall(PetscSFSetUp(invSF, ierr))
    End Subroutine MEF90ConstraintSFCreate

#undef __FUNCT__
#define __FUNCT__ "MEF90VecGlobalToLocalConstraint"
!!!
!!!  
!!!  MEF90VecGlobalToLocalConstraint: do a VecGlobaltoLocal then copy constrained values
!!!  
!!!  (c) 2022      Blaise Bourdin  bourdin@mcmaster.ca
!!!

    Subroutine MEF90VecGlobalToLocalConstraint(g,c,l,ierr)
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
            PetscCall(PetscSectionGetChart(s,pStart,pEnd,ierr))
            Do p = pStart,pEnd-1
                PetscCall(PetscSectionGetConstraintDof(s,p,numConstraint,ierr))
                If (numConstraint > 0) Then
                    PetscCall(VecGetValuesSectionF90(c,s,p,vArray,ierr))
                    PetscCall(VecSetValuesSectionF90(l,s,p,vArray,INSERT_ALL_VALUES,ierr))
                    PetscCall(VecRestoreValuesSectionF90(c,s,p,vArray,ierr))
                End If ! numConstraint
            End Do ! p
        End If
    End Subroutine MEF90VecGlobalToLocalConstraint

#undef __FUNCT__
#define __FUNCT__ "MEF90VecCreateIO"
!!!
!!!  
!!!  MEF90VecCreateIO: create IO Vec
!!!  
!!!  (c) 2022      Alexis Marboeuf  marboeua@mcmaster.ca
!!!

    Subroutine MEF90VecCreateIO(MEF90Ctx,v,bs,sf,ierr)
        Type(tPetscSF),Intent(IN)               :: sf
        Type(tVec),Intent(INOUT)                :: v
        PetscInt,Intent(IN)                     :: bs
        Type(MEF90Ctx_type),Intent(IN)          :: MEF90Ctx
        PetscErrorCode,Intent(INOUT)            :: ierr

        PetscInt                                :: nroots, nleaves
        Type(PetscSFNode),Dimension(:),Pointer  :: iremote
        PetscInt,Dimension(:),Pointer           :: ilocal

        PetscCallA(PetscSFGetGraph(sf,nroots,nleaves,ilocal,iremote,ierr))
        PetscCallA(VecCreateMPI(MEF90Ctx%Comm,nleaves,PETSC_DETERMINE,v,ierr))
        PetscCallA(VecSetBlockSize(v,bs,ierr))

    End Subroutine MEF90VecCreateIO

#undef __FUNCT__
#define __FUNCT__ "CreateNaturalToIOSF_Private"
!!!
!!!  
!!!  CreateNaturalToIOSF_Private: 
!!!  
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!
    
    subroutine CreateNaturalToIOSF_Private(MEF90Ctx,dm,sf,ierr)
        Type(tDM),intent(IN)                    :: dm
        Type(MEF90Ctx_type),Intent(IN)          :: MEF90Ctx
        Type(tPetscSF),intent(OUT)              :: sf
        PetscErrorCode,intent(INOUT)            :: ierr
    
        Type(tVec)                              :: vnat,vio
        Type(PetscLayout)                       :: natMap,ioMap
        Type(PetscSFNode),dimension(:),Pointer  :: remote
        PetscInt,dimension(:),Pointer           :: ioRange, remoteRange
        PetscMPIInt                             :: remoteRank
        PetscInt                                :: nroots, nleaves, globalIndex, i, globalSize, bs
    
        PetscCallA(DMPlexCreateNaturalVector(dm,vnat,ierr))
        PetscCallA(VecGetSize(vnat, globalSize, ierr))
        PetscCallA(VecGetBlockSize(vnat,bs,ierr))
        PetscCallA(VecCreate(MEF90Ctx%Comm,vio,ierr))
        PetscCallA(VecSetBlockSize(vio,bs,ierr))
        PetscCallA(VecSetSizes(vio,PETSC_DETERMINE,globalSize,ierr))
        PetscCallA(VecSetFromOptions(vio,ierr))
        PetscCallA(VecGetLayout(vnat, natMap, ierr))
        PetscCallA(VecGetLayout(vio, ioMap, ierr))
        PetscCall(PetscLayoutGetLocalSize(natMap, nroots, ierr))
        PetscCall(PetscLayoutGetLocalSize(ioMap, nleaves, ierr))
        PetscCall(PetscLayoutGetRangesF90(ioMap, ioRange, ierr))
        PetscCall(PetscLayoutGetRangesF90(natMap, remoteRange, ierr))
        Allocate(remote(nleaves))
        Do i = 0,nleaves-1          
            globalIndex = ioRange(MEF90Ctx%rank+1) + i
            PetscCall(PetscLayoutFindOwner(natMap,globalIndex,remoteRank,ierr))
            remote(i+1)%rank  = remoteRank
            remote(i+1)%index = globalIndex - remoteRange(remoteRank+1)
        End Do
        PetscCall(PetscSFCreate(MEF90Ctx%Comm,sf,ierr))
        PetscCall(PetscObjectSetName(sf,"Natural-To-IO SF",ierr))
        PetscCall(PetscSFSetFromOptions(sf,ierr))
        PetscCall(PetscSFSetGraph(sf,nroots,nleaves,PETSC_NULL_INTEGER,PETSC_COPY_VALUES,remote,PETSC_COPY_VALUES,ierr))
        PetscCall(PetscSFSetUp(sf,ierr))
        PetscCall(PetscSFViewFromOptions(sf,PETSC_NULL_OPTIONS,"-naturaltoio_sf_view",ierr))
        PetscCall(VecDestroy(vio,ierr))
        PetscCall(VecDestroy(vnat,ierr))
        DeAllocate(remote)
    End subroutine CreateNaturalToIOSF_Private
    
#undef __FUNCT__
#define __FUNCT__ "CreateLocalToCGlobalSF_Private"
!!!
!!!  
!!!  CreateLocalToCGlobalSF_Private: 
!!!  
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!
    
    subroutine CreateLocalToCGlobalSF_Private(MEF90Ctx,dm,sf,ierr)
        Type(tDM),intent(IN)                    :: dm
        Type(MEF90Ctx_type),Intent(IN)          :: MEF90Ctx
        Type(tPetscSF),intent(OUT)              :: sf
        PetscErrorCode,intent(INOUT)            :: ierr
    
        Type(tPetscSection)                     :: locSection,gSection
        Type(tPetscSF)                          :: overlapSF,idSF
        Type(PetscSFNode),dimension(:),Pointer  :: remote
        PetscInt,dimension(:),Pointer           :: remoteOffsets
        PetscInt                                :: pStart,pEnd,p,n = 1
    
        PetscCall(DMGetLocalSection(dm,locSection,ierr))
        PetscCall(DMGetPointSF(dm,overlapSF,ierr))
        PetscCall(PetscSectionCreateGlobalSection(locSection,overlapSF,PETSC_TRUE,PETSC_TRUE,gSection,ierr))
        PetscCall(PetscSectionGetChart(locSection,pStart,pEnd,ierr))
        n = pEnd-pStart
        Allocate(remote(n))
        Do p = 1,n
            remote(p)%rank  = MEF90Ctx%rank
            remote(p)%index = p-1
        End Do
        PetscCall(PetscSFCreate(MEF90Ctx%Comm,idSF,ierr))
        PetscCall(PetscSFSetFromOptions(idSF,ierr))
        PetscCall(PetscSFSetGraph(idSF,n,n,PETSC_NULL_INTEGER,PETSC_COPY_VALUES,remote,PETSC_COPY_VALUES,ierr))
        PetscCall(PetscSFSetUp(idSF,ierr))
        PetscCall(PetscSFCreateRemoteOffsetsF90(idSF,locSection,gSection,remoteOffsets,ierr))
        PetscCall(PetscSFCreateSectionSFF90(idSF,locSection,remoteOffsets,gSection,sf,ierr))
        PetscCall(PetscIntArray1dDestroyF90(remoteOffsets,ierr))
        PetscCall(PetscSFSetUp(sf,ierr))
        PetscCall(PetscObjectSetName(sf,"Local-To-CGlobal SF",ierr))
        PetscCall(PetscSFSetFromOptions(sf,ierr))
        PetscCall(PetscSFViewFromOptions(sf,PETSC_NULL_OPTIONS,"-localtocglobal_sf_view",ierr))
        PetscCall(PetscSectionDestroy(gSection,ierr))
        PetscCall(PetscSFDestroy(idSF,ierr))
        DeAllocate(remote)
    
    End subroutine CreateLocalToCGlobalSF_Private
    
#undef __FUNCT__
#define __FUNCT__ "CreateCGlobalToLocalSF_Private"
!!!
!!!  
!!!  CreateCGlobalToLocalSF_Private: 
!!!  
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!
    
    subroutine CreateCGlobalToLocalSF_Private(MEF90Ctx,dm,sf,ierr)
        Type(tDM),intent(IN)                    :: dm
        Type(MEF90Ctx_type),Intent(IN)          :: MEF90Ctx
        Type(tPetscSF),intent(OUT)              :: sf
        PetscErrorCode,intent(INOUT)            :: ierr
    
        Type(tPetscSection)                     :: locSection,gSection
        Type(tPetscSF)                          :: overlapSF,idSF,tempSF,ttempSF
        Type(PetscSFNode),dimension(:),Pointer  :: remote,tempRemote,lgRemote,glRemote
        PetscInt,dimension(:),Pointer           :: tempLocal,lgLocal,glLocal,remoteOffsets
        PetscInt                                :: pStart,pEnd,p,n = 1,lgNRoots,lgNLeaves,tempNRoots,tempNLeaves,glNRoots,glNLeaves
    
        PetscCall(DMGetLocalSection(dm,locSection,ierr))
        PetscCall(DMGetPointSF(dm,overlapSF,ierr))
        PetscCall(PetscSectionCreateGlobalSection(locSection,overlapSF,PETSC_TRUE,PETSC_TRUE,gSection,ierr))
        PetscCall(PetscSectionGetChart(locSection,pStart,pEnd,ierr))
        n = pEnd-pStart
        Allocate(remote(n))
        Do p = 1,n
            remote(p)%rank  = MEF90Ctx%rank
            remote(p)%index = p-1
        End Do
        PetscCall(PetscSFCreate(MEF90Ctx%Comm,idSF,ierr))
        PetscCall(PetscSFSetFromOptions(idSF,ierr))
        PetscCall(PetscSFSetGraph(idSF,n,n,PETSC_NULL_INTEGER,PETSC_COPY_VALUES,remote,PETSC_COPY_VALUES,ierr))
        PetscCall(PetscSFSetUp(idSF,ierr))
        PetscCall(PetscSFCreateRemoteOffsetsF90(idSF,locSection,gSection,remoteOffsets,ierr))
        PetscCall(PetscSFCreateSectionSFF90(idSF,locSection,remoteOffsets,gSection,tempSF,ierr))
        PetscCall(PetscSFCreateInverseSF(tempSF,sf,ierr))
        PetscCall(PetscIntArray1dDestroyF90(remoteOffsets,ierr))
        PetscCall(PetscSFDestroy(tempSF,ierr))
        If (MEF90Ctx%NumProcs > 1) Then
            PetscCall(PetscSFGetGraph(sf,lgNRoots,lgNLeaves,lgLocal,lgRemote,ierr))
            PetscCall(PetscSFCreateRemoteOffsetsF90(overlapSF,locSection,gSection,remoteOffsets,ierr))
            PetscCall(PetscSFCreateSectionSFF90(overlapSF,locSection,remoteOffsets,gSection,ttempSF,ierr))
            PetscCall(PetscSFCreateInverseSF(ttempSF,tempSF,ierr))
            PetscCall(PetscIntArray1dDestroyF90(remoteOffsets,ierr))
            PetscCall(PetscSFDestroy(ttempSF,ierr))
            PetscCall(PetscSFGetGraph(tempSF,tempNRoots,tempNLeaves,tempLocal,tempRemote,ierr))
            glNRoots = lgNRoots
            glNLeaves = lgNLeaves + tempNLeaves
            Allocate(glLocal(glNLeaves))
            Allocate(glRemote(glNLeaves))
            If (loc(lgLocal) .ne. loc(PETSC_NULL_INTEGER)) Then
                Do p = 1,lgNLeaves
                    glLocal(p) = lgLocal(p)
                    glRemote(p)%rank  = lgRemote(p)%rank
                    glRemote(p)%index = lgRemote(p)%index
                End Do
            Else 
                Do p = 1,lgNLeaves
                    glLocal(p) = p-1
                    glRemote(p)%rank  = lgRemote(p)%rank
                    glRemote(p)%index = lgRemote(p)%index
                End Do
            End If
            If (loc(tempLocal) .ne. loc(PETSC_NULL_INTEGER)) Then
                Do p = 1,tempNLeaves
                    glLocal(p+lgNLeaves)        = tempLocal(p)
                    glRemote(p+lgNLeaves)%rank  = tempRemote(p)%rank
                    glRemote(p+lgNLeaves)%index = tempRemote(p)%index
                End Do
            Else
                Do p = 1,tempNLeaves
                    glLocal(p+lgNLeaves)        = p+lgNLeaves-1
                    glRemote(p+lgNLeaves)%rank  = tempRemote(p)%rank
                    glRemote(p+lgNLeaves)%index = tempRemote(p)%index
                End Do
            End If
            PetscCall(PetscSFDestroy(tempSF,ierr))
            PetscCall(PetscSFDestroy(sf,ierr))
            PetscCall(PetscSFCreate(MEF90Ctx%Comm,sf,ierr))
            PetscCall(PetscSFSetFromOptions(sf,ierr))
            PetscCall(PetscSFSetGraph(sf,glNRoots,glNLeaves,glLocal,PETSC_COPY_VALUES,glRemote,PETSC_COPY_VALUES,ierr))
            PetscCall(PetscSFSetUp(sf,ierr))
            DeAllocate(glLocal)
            DeAllocate(glRemote)
        End If
        PetscCall(PetscSectionDestroy(gSection,ierr))
        PetscCall(PetscSFDestroy(idSF,ierr))
        PetscCall(PetscObjectSetName(sf,"CGlobal-To-Local SF",ierr))
        PetscCall(PetscSFViewFromOptions(sf,PETSC_NULL_OPTIONS,"-cglobaltolocal_sf_view",ierr))
        DeAllocate(remote)
    End subroutine CreateCGlobalToLocalSF_Private                    
End Module m_MEF90_DMPlex
