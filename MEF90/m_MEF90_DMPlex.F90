module m_MEF90_DMPlex
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscsf.h"
   use m_MEF90_Parameters
   use m_MEF90_Utils
   use m_MEF90_LinAlg
   use m_MEF90_Elements
   use m_MEF90_Ctx
   use, intrinsic :: iso_c_binding
#ifdef MEF90_HAVE_SYMENGINEF90
   use symengine
#endif
   implicit none(type)

   enum, bind(c)
      enumerator  :: MEF90CellSetType = 1, &
         MEF90FaceSetType, &
         !MEF90EdgeSetType,     &
         MEF90VertexSetType
   end enum

    !!! Not sure why this has to be made explicitly public
   PetscEnum, dimension(3), parameter, public :: MEF90SetType = [ &
                                                                MEF90CellSetType, &
                                                                MEF90FaceSetType, &
                                                                !MEF90EdgeSetType,                                       &
                                                                MEF90VertexSetType]

   character(len=MEF90MXSTRLEN), parameter :: MEF90CellSetLabelName = 'Cell Sets  '
   character(len=MEF90MXSTRLEN), parameter :: MEF90FaceSetLabelName = 'Face Sets  '
   !Character(len=MEF90MXSTRLEN),Parameter :: MEF90EdgeSetLabelName   = 'Edge Sets  '
   character(len=MEF90MXSTRLEN), parameter :: MEF90VertexSetLabelName = 'Vertex Sets'

   character(len=MEF90MXSTRLEN), dimension(3), parameter  :: MEF90SetLabelName = &
                                                             [MEF90CellSetLabelName, &
                                                              MEF90FaceSetLabelName, &
                                                              !MEF90EdgeSetLabelName, &
                                                              MEF90VertexSetLabelName]

   character(len=MEF90MXSTRLEN), parameter :: MEF90CellSetprefix = 'cs'
   character(len=MEF90MXSTRLEN), parameter :: MEF90FaceSetprefix = 'fs'
   !Character(len=MEF90MXSTRLEN),Parameter :: MEF90EdgeSetprefix   = 'es'
   character(len=MEF90MXSTRLEN), parameter :: MEF90VertexSetprefix = 'vs'

   character(len=MEF90MXSTRLEN), dimension(3), parameter  :: MEF90SetPrefix = &
                                                             [MEF90CellSetPrefix, &
                                                              MEF90FaceSetPrefix, &
                                                              !MEF90EdgeSetPrefix, &
                                                              MEF90VertexSetPrefix]

   private
   public :: MEF90CellSetLabelName, &
             MEF90FaceSetLabelName, &
             !MEF90EdgeSetLabelName,                                                 &
             MEF90VertexSetLabelName, &
             MEF90SetLabelName, &
             MEF90CellSetPrefix, &
             MEF90FaceSetPrefix, &
             !MEF90EdgeSetPrefix,                                                    &
             MEF90VertexSetPrefix, &
             MEF90SetPrefix, &
             MEF90CellSetType, &
             MEF90FaceSetType, &
             !MEF90EdgeSetType,                                                      &
             MEF90VertexSetType, &
             MEF90SectionAllocateDof, MEF90SectionAllocateDofSet, &
             MEF90SetupConstraintTableSet, MEF90SectionAllocateConstraint, &
             MEF90CellSectionCreate, &
             MEF90VecCopySF, MEF90IOSFCreate, MEF90FaceSetIOSFCreate, &
             MEF90ConstraintSFCreate, MEF90VecGlobalToLocalConstraint, &
             MEF90VecCreateIO, &
             MEF90CreateLocalVector, &
             MEF90CreateBoundaryLocalVector, &
             MEF90CreateCellVector, &
             MEF90CreateBoundaryCellVector, &
             MEF90VecSetBCValuesFromOptions, &
             MEF90VecSetValuesFromOptions, &
             MEF90VecSetBCValuesFromOptionsExpr, &
             MEF90VecSetValuesFromOptionsExpr, &
             MEF90VecGetClosureSize
contains

#undef __FUNCT__
#define __FUNCT__ "MEF90CreateLocalVector"
!!!
!!!
!!!  MEF90CreateLocalVector: create a Vec associated with a FE space and constraints
!!!      cell   set BC are obtained from the command line option -cs<set ID>_<name>BC [bool], [bool], ...
!!!      face   set BC are obtained from the command line option -fs<set ID>_<name>BC [bool], [bool], ...
!!!      vertex set BC are obtained from the command line option -vs<set ID>_<name>BC [bool], [bool], ...
!!!
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine MEF90CreateLocalVector(dm, elemFamily, elemOrder, sDim, name, V, ierr)
      type(tDM), intent(IN)                    :: dm
      PetscEnum, intent(IN)                    :: elemFamily
      PetscInt, intent(IN)                     :: elemOrder, sDim
      character(len=MEF90MXSTRLEN), intent(IN) :: name
      type(tVec), intent(OUT)                  :: V
      PetscErrorCode, intent(INOUT)            :: ierr

      type(tPetscSection)                     :: sectionV
      type(tDM)                               :: dmV
      PetscInt                                :: fieldV = 0
      PetscInt                                :: pStart, pEnd
      MPI_Comm                                :: comm
      PetscInt                                :: set
      PetscEnum                               :: setType
      PetscInt, dimension(:), pointer           :: setID, pointID
      type(tIS)                               :: setIS, pointIS
      type(MEF90ElementType)                  :: elemType
      DMPolytopeType                          :: cellType
      PetscBool, dimension(:, :), pointer        :: constraintTruthTable
      PetscBool, dimension(:), pointer          :: setConstraints
      PetscInt                                :: numBC
      PetscBool                               :: flg
      character(len=MEF90MXSTRLEN)            :: BCOptionName
      type(tPetscSF)                          :: naturalPointSF, naturalSF

      PetscCall(DMClone(dm, dmV, ierr))
      PetscCall(PetscObjectSetName(dmv, name, ierr))
      PetscCall(DMGetUseNatural(dm, flg, ierr))
      PetscCall(DMSetUseNatural(dmV, flg, ierr))

      PetscCall(PetscObjectGetComm(dmV, comm, ierr))

      PetscCall(PetscSectionCreate(comm, sectionV, ierr))
      PetscCall(PetscObjectSetName(sectionV, name, ierr))
      PetscCall(PetscSectionSetNumFields(sectionV, 1_ki, ierr))
      PetscCall(PetscSectionSetFieldName(sectionV, fieldV, trim(name), ierr))
      PetscCall(PetscSectionSetFieldComponents(sectionV, fieldV, sdim, ierr))
      PetscCall(DMPlexGetChart(dmV, pStart, pEnd, ierr))
      PetscCall(PetscSectionSetChart(sectionV, pStart, pEnd, ierr))

      PetscCall(DMGetLabelIdIS(dmV, MEF90CellSetLabelName, setIS, ierr))
        !!! Get a GLOBAL cell set IS
      ! PetscCall(MEF90ISAllGatherMerge(comm,setIS,ierr))
      if (.not. PetscObjectIsNull(setIS)) then
         PetscCall(ISGetIndices(setIS, setID, ierr))
         do set = 1, size(setID)
                !!! Get cell type in order to pick the proper element type.
                !!! We assume that all cells in a set have the same type, so all we need it to query the first cell in the set
            PetscCall(DMGetStratumIS(dmV, MEF90CellSetLabelName, setID(set), pointIS, ierr))
            if (.not. PetscObjectIsNull(pointIS)) then
               PetscCall(ISGetIndices(pointIS, pointID, ierr))
               PetscCall(DMPlexGetCellType(dmV, pointID(1), cellType, ierr))
               PetscCall(MEF90ElementGetType(elemFamily, elemOrder, cellType, elemType, ierr))
               PetscCall(MEF90SectionAllocateDofSet(dmV, MEF90CellSetType, setID(set), elemType, sdim, sectionV, ierr))
               PetscCall(ISRestoreIndices(pointIS, pointID, ierr))
            end if ! pointIS
            PetscCall(ISDestroy(pointIS, ierr))
         end do ! set
         PetscCall(ISRestoreIndices(setIS, setID, ierr))
      end if ! setIS
      PetscCall(ISDestroy(setIS, ierr))

      allocate (ConstraintTruthTable(pEnd, sDim))
      ConstraintTruthTable = .false.
      allocate (setConstraints(sDim))

      do setType = 1, size(MEF90SetType)
         PetscCall(DMGetLabelIdIS(dm, MEF90SetLabelName(setType), setIS, ierr))
         ! PetscCall(MEF90ISAllGatherMerge(comm,setIS,ierr))
         if (.not. PetscObjectIsNull(setIS)) then
            PetscCall(ISGetIndices(setIS, setID, ierr))
            do set = 1, size(setID)
               setConstraints = .false.
               write (BCOptionName, '("-",a2,I4.4,"_",a,"BC")') MEF90SetPrefix(setType), setID(set), trim(name)
               numBC = sDim
               PetscCall(PetscOptionsGetBoolArray(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, trim(BCOptionName), setConstraints, numBC, flg, ierr))
               PetscCall(MEF90SetupConstraintTableSet(dmV, sectionV, MEF90SetType(setType), setID(set), setConstraints, ConstraintTruthTable, ierr))
            end do
            PetscCall(ISRestoreIndices(setIS, setID, ierr))
         end if ! setIS
         PetscCall(ISDestroy(setIS, ierr))
      end do ! setType

      PetscCall(MEF90SectionAllocateConstraint(dmV, ConstraintTruthTable, sectionV, ierr))
      deallocate (ConstraintTruthTable)
      deallocate (setConstraints)

      PetscCall(DMSetLocalSection(dmV, sectionV, ierr))

      PetscCall(DMPlexGetMigrationSF(dm, naturalPointSF, ierr))
      if (.not. PetscObjectIsNull(naturalPointSF)) then
         PetscCall(DMPlexSetMigrationSF(dmV, naturalPointSF, ierr))
         PetscCall(DMPlexCreateGlobalToNaturalSF(dmV, PETSC_NULL_SECTION, naturalPointSF, naturalSF, ierr))
         PetscCall(DMSetNaturalSF(dmV, naturalSF, ierr))
      end if

#ifdef PETSC_USE_DEBUG
      write (BCOptionName, '("-",a,"_section_view")') trim(name)
      PetscCall(PetscSectionViewFromOptions(sectionV, PETSC_NULL_OBJECT, BCOptionName, ierr))
      if (.not. PetscObjectIsNull(naturalSF)) then
         write (BCOptionName, '("-",a,"_naturalSF_view")') trim(name)
         PetscCall(PetscSFViewFromOptions(naturalSF, PETSC_NULL_OBJECT, BCOptionName, ierr))
      end if
#endif
      PetscCall(DMCreateLocalVector(dmV, V, ierr))
      PetscCall(PetscObjectSetName(V, name, ierr))
      PetscCall(PetscSectionDestroy(sectionV, ierr))
      if (.not. PetscObjectIsNull(naturalSF)) then
         PetscCall(PetscSFDestroy(naturalSF, ierr))
      end if
      PetscCall(DMDestroy(dmV, ierr))
   end subroutine MEF90CreateLocalVector

#undef __FUNCT__
#define __FUNCT__ "MEF90CreateCellVector"
!!!
!!!
!!!  MEF90CreateCellVector: create a Vec for a cell-based vector
!!!
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine MEF90CreateCellVector(dm, sDim, name, V, ierr)
      type(tDM), intent(IN)                    :: dm
      PetscInt, intent(IN)                     :: sDim
      character(len=MEF90MXSTRLEN), intent(IN) :: name
      type(tVec), intent(OUT)                  :: V
      PetscErrorCode, intent(INOUT)            :: ierr

      PetscInt                                :: dim
      type(tPetscSection)                     :: sectionV
      type(tDM)                               :: dmV
      PetscInt                                :: fieldV = 0
      PetscInt                                :: pStart, pEnd
      MPI_Comm                                :: comm
      PetscInt                                :: set
      PetscInt, dimension(:), pointer           :: setID
      type(tIS)                               :: setIS
      type(MEF90ElementType)                  :: elemType
      PetscBool                               :: flg
      type(tPetscSF)                          :: naturalPointSF, naturalSF

      PetscCall(DMClone(dm, dmV, ierr))
      PetscCall(PetscObjectSetName(dmv, name, ierr))
      PetscCall(DMGetUseNatural(dm, flg, ierr))
      PetscCall(DMSetUseNatural(dmV, flg, ierr))

      PetscCall(PetscObjectGetComm(dmV, comm, ierr))

      PetscCall(PetscSectionCreate(comm, sectionV, ierr))
      PetscCall(PetscObjectSetName(sectionV, name, ierr))
      PetscCall(PetscSectionSetNumFields(sectionV, 1_ki, ierr))
      PetscCall(PetscSectionSetFieldName(sectionV, fieldV, trim(name), ierr))
      PetscCall(PetscSectionSetFieldComponents(sectionV, fieldV, sdim, ierr))
      PetscCall(DMPlexGetChart(dmV, pStart, pEnd, ierr))
      PetscCall(PetscSectionSetChart(sectionV, pStart, pEnd, ierr))

      PetscCall(DMGetDimension(dm, dim, ierr))
      if (dim == 2) then
         elemType = MEF90P0Lagrange2D
      else
         elemType = MEF90P0Lagrange3D
      end if

      PetscCall(DMGetLabelIdIS(dmV, MEF90CellSetLabelName, setIS, ierr))
        !!! Get a GLOBAL cell set IS
      ! PetscCall(MEF90ISAllGatherMerge(comm,setIS,ierr))
      if (.not. PetscObjectIsNull(setIS)) then
         PetscCall(ISGetIndices(setIS, setID, ierr))
         do set = 1, size(setID)
            PetscCall(MEF90SectionAllocateDofSet(dmV, MEF90CellSetType, setID(set), elemType, sdim, sectionV, ierr))
         end do ! set
         PetscCall(ISRestoreIndices(setIS, setID, ierr))
      end if ! setIS
      PetscCall(ISDestroy(setIS, ierr))
      PetscCall(PetscSectionSetup(sectionV, ierr))

      PetscCall(DMSetLocalSection(dmV, sectionV, ierr))

      PetscCall(DMPlexGetMigrationSF(dm, naturalPointSF, ierr))
      if (.not. PetscObjectIsNull(naturalPointSF)) then
         PetscCall(DMPlexSetMigrationSF(dmV, naturalPointSF, ierr))
         PetscCall(DMPlexCreateGlobalToNaturalSF(dmV, PETSC_NULL_SECTION, naturalPointSF, naturalSF, ierr))
         PetscCall(DMSetNaturalSF(dmV, naturalSF, ierr))
      end if

#ifdef PETSC_USE_DEBUG
      debugBlock: block
         character(len=MEF90MXSTRLEN)            :: BCoptionName
         write (BCOptionName, '("-",a,"_section_view")') trim(name)
         PetscCall(PetscSectionViewFromOptions(sectionV, PETSC_NULL_OBJECT, BCOptionName, ierr))
         if (.not. PetscObjectIsNull(naturalSF)) then
            write (BCOptionName, '("-",a,"_naturalSF_view")') trim(name)
            PetscCall(PetscSFViewFromOptions(naturalSF, PETSC_NULL_OBJECT, BCOptionName, ierr))
         end if
      end block debugBlock
#endif
      PetscCall(DMCreateLocalVector(dmV, V, ierr))
      PetscCall(PetscObjectSetName(V, name, ierr))
      PetscCall(DMDestroy(dmV, ierr))
      PetscCall(PetscSectionDestroy(sectionV, ierr))
      if (.not. PetscObjectIsNull(naturalSF)) then
         PetscCall(PetscSFDestroy(naturalSF, ierr))
      end if
   end subroutine MEF90CreateCellVector

#undef __FUNCT__
#define __FUNCT__ "MEF90CreateBoundaryLocalVector"
!!!
!!!
!!!  MEF90CreateBoundaryLocalVector: create a Vec associated with a FE space over the boundary of a domain (face sets)
!!!                                  and constraints
!!!      face   set BC are obtained from the command line option -fs<set ID>_<name>BC [bool], [bool], ...
!!!      vertex set BC are obtained from the command line option -vs<set ID>_<name>BC [bool], [bool], ...
!!!
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine MEF90CreateBoundaryLocalVector(dm, elemFamily, elemOrder, sDim, name, V, ierr)
      type(tDM), intent(IN)                    :: dm
      PetscEnum, intent(IN)                    :: elemFamily
      PetscInt, intent(IN)                     :: elemOrder, sDim
      character(len=MEF90MXSTRLEN), intent(IN) :: name
      type(tVec), intent(OUT)                  :: V
      PetscErrorCode, intent(INOUT)            :: ierr

      type(tPetscSection)                     :: sectionV
      type(tDM)                               :: dmV
      PetscInt                                :: fieldV = 0
      PetscInt                                :: pStart, pEnd
      MPI_Comm                                :: comm
      PetscInt                                :: set
      PetscEnum                               :: setType
      PetscInt, dimension(:), pointer           :: setID, pointID
      type(tIS)                               :: setIS, pointIS
      type(MEF90ElementType)                  :: elemType
      DMPolytopeType                          :: cellType
      PetscBool, dimension(:, :), pointer        :: constraintTruthTable
      PetscBool, dimension(:), pointer          :: setConstraints
      PetscInt                                :: numBC
      PetscBool                               :: flg
      character(len=MEF90MXSTRLEN)            :: BCOptionName

      PetscCall(DMClone(dm, dmV, ierr))
      PetscCall(PetscObjectSetName(dmv, name, ierr))
      PetscCall(DMGetUseNatural(dm, flg, ierr))
      PetscCall(DMSetUseNatural(dmV, flg, ierr))

      PetscCall(PetscObjectGetComm(dmV, comm, ierr))

      PetscCall(PetscSectionCreate(comm, sectionV, ierr))
      PetscCall(PetscObjectSetName(sectionV, name, ierr))
      PetscCall(PetscSectionSetNumFields(sectionV, 1_ki, ierr))
      PetscCall(PetscSectionSetFieldName(sectionV, fieldV, trim(name), ierr))
      PetscCall(PetscSectionSetFieldComponents(sectionV, fieldV, sdim, ierr))
      PetscCall(DMPlexGetChart(dmV, pStart, pEnd, ierr))
      PetscCall(PetscSectionSetChart(sectionV, pStart, pEnd, ierr))

      PetscCall(DMGetLabelIdIS(dmV, MEF90FaceSetLabelName, setIS, ierr))
        !!! Get a GLOBAL cell set IS
      ! PetscCall(MEF90ISAllGatherMerge(comm,setIS,ierr))
      if (.not. PetscObjectIsNull(setIS)) then
         PetscCall(ISGetIndices(setIS, setID, ierr))
         do set = 1, size(setID)
                !!! Get cell type in order to pick the proper element type.
                !!! We assume that all cells in a set have the same type, so all we need it to query the first cell in the set
            PetscCall(DMGetStratumIS(dmV, MEF90FaceSetLabelName, setID(set), pointIS, ierr))
            if (.not. PetscObjectIsNull(pointIS)) then
               PetscCall(ISGetIndices(pointIS, pointID, ierr))
               PetscCall(DMPlexGetCellType(dmV, pointID(1), cellType, ierr))
               PetscCall(MEF90ElementGetTypeBoundary(elemFamily, elemOrder, cellType, elemType, ierr))
               PetscCall(MEF90SectionAllocateDofSet(dmV, MEF90FaceSetType, setID(set), elemType, sdim, sectionV, ierr))
               PetscCall(ISRestoreIndices(pointIS, pointID, ierr))
            end if ! pointIS
            PetscCall(ISDestroy(pointIS, ierr))
         end do ! set
         PetscCall(ISRestoreIndices(setIS, setID, ierr))
      end if ! setIS
      PetscCall(ISDestroy(setIS, ierr))

      allocate (ConstraintTruthTable(pEnd, sDim))
      ConstraintTruthTable = .false.
      allocate (setConstraints(sDim))

      do setType = 1, size(MEF90SetType)
         PetscCall(DMGetLabelIdIS(dm, MEF90SetLabelName(setType), setIS, ierr))
         ! PetscCall(MEF90ISAllGatherMerge(comm,setIS,ierr))
         if (.not. PetscObjectIsNull(setIS)) then
            PetscCall(ISGetIndices(setIS, setID, ierr))
            do set = 1, size(setID)
               setConstraints = .false.
               write (BCOptionName, '("-",a2,I4.4,"_",a,"BC")') MEF90SetPrefix(setType), setID(set), trim(name)
               numBC = sDim
               PetscCall(PetscOptionsGetBoolArray(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, trim(BCOptionName), setConstraints, numBC, flg, ierr))
               PetscCall(MEF90SetupConstraintTableSet(dmV, sectionV, MEF90SetType(setType), setID(set), setConstraints, ConstraintTruthTable, ierr))
            end do
            PetscCall(ISRestoreIndices(setIS, setID, ierr))
         end if ! setIS
         PetscCall(ISDestroy(setIS, ierr))
      end do ! setType

      PetscCall(MEF90SectionAllocateConstraint(dmV, ConstraintTruthTable, sectionV, ierr))
      deallocate (ConstraintTruthTable)
      deallocate (setConstraints)

      PetscCall(DMSetLocalSection(dmV, sectionV, ierr))
#ifdef PETSC_USE_DEBUG
      debugBlock: block
         character(len=MEF90MXSTRLEN)            :: BCoptionName
         write (BCOptionName, '("-",a,"_section_view")') trim(name)
         PetscCall(PetscSectionViewFromOptions(sectionV, PETSC_NULL_OBJECT, BCOptionName, ierr))
      end block debugBlock
#endif
      PetscCall(DMCreateLocalVector(dmV, V, ierr))
      PetscCall(PetscObjectSetName(V, name, ierr))
      PetscCall(PetscSectionDestroy(sectionV, ierr))
      PetscCall(DMDestroy(dmV, ierr))
   end subroutine MEF90CreateBoundaryLocalVector

#undef __FUNCT__
#define __FUNCT__ "MEF90CreateBoundaryCellVector"
!!!
!!!
!!!  MEF90CreateBoundaryCellVector: create a Vec for a cell-based vector
!!!
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine MEF90CreateBoundaryCellVector(dm, sDim, name, V, ierr)
      type(tDM), intent(IN)                    :: dm
      PetscInt, intent(IN)                     :: sDim
      character(len=MEF90MXSTRLEN), intent(IN) :: name
      type(tVec), intent(OUT)                  :: V
      PetscErrorCode, intent(INOUT)            :: ierr

      PetscInt                                :: dim
      type(tPetscSection)                     :: sectionV
      type(tDM)                               :: dmV
      PetscInt                                :: fieldV = 0
      PetscInt                                :: pStart, pEnd
      MPI_Comm                                :: comm
      PetscInt                                :: set
      PetscInt, dimension(:), pointer           :: setID
      type(tIS)                               :: setIS
      type(MEF90ElementType)                  :: elemType
      PetscBool                               :: flg

      PetscCall(DMClone(dm, dmV, ierr))
      PetscCall(PetscObjectSetName(dmv, name, ierr))
      PetscCall(DMGetUseNatural(dm, flg, ierr))
      PetscCall(DMSetUseNatural(dmV, flg, ierr))

      PetscCall(PetscObjectGetComm(dmV, comm, ierr))

      PetscCall(PetscSectionCreate(comm, sectionV, ierr))
      PetscCall(PetscObjectSetName(sectionV, name, ierr))
      PetscCall(PetscSectionSetNumFields(sectionV, 1_ki, ierr))
      PetscCall(PetscSectionSetFieldName(sectionV, fieldV, trim(name), ierr))
      PetscCall(PetscSectionSetFieldComponents(sectionV, fieldV, sdim, ierr))
      PetscCall(DMPlexGetChart(dmV, pStart, pEnd, ierr))
      PetscCall(PetscSectionSetChart(sectionV, pStart, pEnd, ierr))

      PetscCall(DMGetDimension(dm, dim, ierr))
      if (dim == 2) then
         elemType = MEF90P0Lagrange2DBoundary
      else
         elemType = MEF90P0Lagrange3DBoundary
      end if

      PetscCall(DMGetLabelIdIS(dmV, MEF90FaceSetLabelName, setIS, ierr))
        !!! Get a GLOBAL cell set IS
      ! PetscCall(MEF90ISAllGatherMerge(comm,setIS,ierr))
      if (.not. PetscObjectIsNull(setIS)) then
         PetscCall(ISGetIndices(setIS, setID, ierr))
         do set = 1, size(setID)
            PetscCall(MEF90SectionAllocateDofSet(dmV, MEF90FaceSetType, setID(set), elemType, sdim, sectionV, ierr))
         end do ! set
         PetscCall(ISRestoreIndices(setIS, setID, ierr))
      end if ! setIS
      PetscCall(ISDestroy(setIS, ierr))
      PetscCall(PetscSectionSetup(sectionV, ierr))

      PetscCall(DMSetLocalSection(dmV, sectionV, ierr))
#ifdef PETSC_USE_DEBUG
      debugBlock: block
         character(len=MEF90MXSTRLEN)            :: BCoptionName
         write (BCOptionName, '("-",a,"_section_view")') trim(name)
         PetscCall(PetscSectionViewFromOptions(sectionV, PETSC_NULL_OBJECT, BCOptionName, ierr))
      end block debugBlock
#endif
      PetscCall(DMCreateLocalVector(dmV, V, ierr))
      PetscCall(PetscObjectSetName(V, name, ierr))
      PetscCall(PetscSectionDestroy(sectionV, ierr))
      PetscCall(DMDestroy(dmV, ierr))
   end subroutine MEF90CreateBoundaryCellVector

#undef __FUNCT__
#define __FUNCT__ "MEF90VecGetClosureSize"
!!!
!!!
!!!  MEF90VecGetClosureSize: Associates dof to a section
!!!
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90VecGetClosureSize(v, p, clSize, ierr)
      type(tVec), intent(IN)              :: v
      PetscInt, intent(IN)                :: p
      PetscInt, intent(OUT)               :: clSize
      PetscErrorCode, intent(INOUT)       :: ierr

      type(tDM)                          :: dm
      type(tPetscSection)                :: section
      PetscInt, dimension(:), pointer      :: closure
      PetscInt                           :: point, numDof, numClosure

      clSize = 0
      PetscCall(VecGetDM(v, dm, ierr))
      PetscCall(DMGetLocalSection(dm, section, ierr))

      PetscCall(DMPlexGetTransitiveClosure(dm, p, PETSC_TRUE, numClosure, closure, ierr))
      if (size(closure) > 0) then
         do point = 1, size(closure), 2
            PetscCall(PetscSectionGetDof(section, closure(point), numDof, ierr))
            clSize = clSize + numDof
         end do
      end if
      PetscCall(DMPlexRestoreTransitiveClosure(dm, p, PETSC_TRUE, numClosure, closure, ierr))
   end subroutine MEF90VecGetClosureSize

#undef __FUNCT__
#define __FUNCT__ "MEF90SectionAllocateDof"
!!!
!!!
!!!  MEF90SectionAllocateDof: Associates dof to a section
!!!
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90SectionAllocateDof(dm, setType, elemType, numComponents, section, ierr)
      type(tDM), intent(IN)               :: dm
      PetscEnum, intent(IN)               :: setType
      type(MEF90ElementType), intent(IN)  :: elemType
      PetscInt, intent(IN)                :: numComponents
      type(tPetscSection), intent(INOUT)  :: section
      PetscErrorCode, intent(INOUT)       :: ierr

      type(tIS)                          :: setIS
      PetscInt, dimension(:), pointer      :: setID
      PetscInt                           :: set

      PetscCall(DMGetLabelIdIS(dm, MEF90SetLabelName(setType), setIS, ierr))
      if (.not. PetscObjectIsNull(setIS)) then
         PetscCall(ISGetIndices(setIS, setID, ierr))
         do set = 1, size(setID)
            PetscCall(MEF90SectionAllocateDofSet(dm, setType, setID(set), elemType, numComponents, Section, ierr))
         end do
         PetscCall(ISRestoreIndices(setIS, setID, ierr))
      end if ! setIS
      PetscCall(ISDestroy(setIS, ierr))
   end subroutine MEF90SectionAllocateDof

#undef __FUNCT__
#define __FUNCT__ "MEF90_SectionAllocateDofSet"
!!!
!!!
!!!  MEF90SectionAllocateDofSet: Associates dof to a section in a set
!!!
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90SectionAllocateDofSet(dm, setType, setID, elemType, numComponents, section, ierr)
      type(tDM), intent(IN)               :: dm
      PetscEnum, intent(IN)               :: setType
      PetscInt, intent(IN)                :: setID
      type(MEF90ElementType), intent(IN)  :: elemType
      PetscInt, intent(IN)                :: numComponents
      type(tPetscSection), intent(INOUT)  :: section
      PetscErrorCode, intent(INOUT)       :: ierr

      type(tIS)                          :: setPointIS
      PetscInt, dimension(:), pointer      :: setPointID
      PetscInt                           :: point
      PetscInt, dimension(:), pointer      :: closure
      PetscInt                           :: p, depth
      PetscInt                           :: field = 0

      PetscCall(DMGetStratumIS(dm, MEF90SetLabelName(setType), setID, setPointIS, ierr))
      if (.not. PetscObjectIsNull(setPointIS)) then
         PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
            !!! This can probably be optimized by allocating closure outside of the loop
            !!! But I can't figure out how it is done at the moment.
         nullify (closure)
         do point = 1, size(setPointID)
            PetscCall(DMPlexGetTransitiveClosure(dm, setPointID(point), PETSC_TRUE, PETSC_NULL_INTEGER, closure, ierr))
            do p = 1, size(closure), 2
               PetscCall(DMPlexGetPointDepth(dm, closure(p), depth, ierr))
               if (elemType%numDofs(depth + 1) > 0) then
                  PetscCall(PetscSectionSetDof(section, closure(p), elemType%numDofs(depth + 1) * numComponents, ierr))
                  PetscCall(PetscSectionSetFieldDof(section, closure(p), field, elemType%numDofs(depth + 1) * numComponents, ierr))
               end if
            end do! p
            PetscCall(DMPlexRestoreTransitiveClosure(dm, setPointID(point), PETSC_TRUE, PETSC_NULL_INTEGER, closure, ierr))
         end do! cell
         PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
      end if ! setPointIS
      PetscCall(ISDestroy(setPointIS, ierr))
   end subroutine MEF90SectionAllocateDofSet

#undef __FUNCT__
#define __FUNCT__ "MEF90CellSectionCreate"
!!!
!!!
!!!  MEF90CellSectionCreate: create a section with numComponent dof at cells
!!!
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!

   subroutine MEF90CellSectionCreate(dm, numComponents, section, ierr)
      type(tDM), intent(IN)               :: dm
      PetscInt, intent(IN)                :: numComponents
      type(tPetscSection), intent(INOUT)  :: section
      PetscErrorCode, intent(INOUT)       :: ierr

      type(tIS)                          :: setIS, setPointIS
      PetscInt, dimension(:), pointer      :: setID, setPointID
      PetscInt                           :: set, point
      PetscInt, dimension(:), pointer      :: closure
      PetscInt                           :: p, depth, dim

      PetscCall(DMGetDimension(dm, dim, ierr))
      PetscCall(DMGetLabelIdIS(dm, 'Cell Sets', setIS, ierr))
      if (.not. PetscObjectIsNull(setIS)) then
         PetscCall(ISGetIndices(setIS, setID, ierr))
         do set = 1, size(setID)
            PetscCall(DMGetStratumIS(dm, 'Cell Sets', setID(set), setPointIS, ierr))
            if (.not. PetscObjectIsNull(setPointIS)) then
               PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
                    !!! This can probably be optimized by allocating closure outside of the loop
                    !!! But I can't figure out how it is done at the moment.
               nullify (closure)
               do point = 1, size(setPointID)
                  PetscCall(DMPlexGetTransitiveClosure(dm, setPointID(point), PETSC_TRUE, PETSC_NULL_INTEGER, closure, ierr))
                  do p = 1, size(closure), 2
                     PetscCall(DMPlexGetPointDepth(dm, closure(p), depth, ierr))
                     if (depth == dim) then
                        PetscCall(PetscSectionSetDof(section, closure(p), numComponents, ierr))
                     end if
                  end do! p
                  PetscCall(DMPlexRestoreTransitiveClosure(dm, setPointID(point), PETSC_TRUE, PETSC_NULL_INTEGER, closure, ierr))
               end do! cell
               PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
            end if ! setPointIS
            PetscCall(ISDestroy(setPointIS, ierr))
         end do ! set
         PetscCall(ISRestoreIndices(setIS, setID, ierr))
      end if ! setIS
      PetscCall(ISDestroy(setIS, ierr))
      PetscCall(PetscSectionSetup(section, ierr))
   end subroutine MEF90CellSectionCreate

#undef __FUNCT__
#define __FUNCT__ "MEF90SetupConstraintTableSet"
!!!
!!!
!!!  MEF90SetupConstraintTableSet: Build the contribution of a set to the constraint table
!!!
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90SetupConstraintTableSet(dm, section, setType, setID, constraints, table, ierr)
      type(tDM), intent(IN)               :: dm
      type(tPetscSection), intent(INOUT)  :: section
      PetscEnum, intent(IN)               :: setType
      PetscInt, intent(IN)                :: setID
      PetscBool, dimension(:), pointer     :: constraints
      PetscBool, dimension(:, :), pointer   :: table
      PetscErrorCode, intent(INOUT)       :: ierr

      type(tIS)                          :: setPointIS
      PetscInt, dimension(:), pointer      :: setPointID
      PetscInt                           :: point, numDof
      PetscInt, dimension(:), pointer      :: closure
      PetscInt                           :: p

      PetscCall(DMGetStratumIS(dm, MEF90SetLabelName(setType), setID, setPointIS, ierr))
      if (.not. PetscObjectIsNull(setPOintIS)) then
         PetscCall(ISGetIndices(setPointIS, setPointID, ierr))
            !!! This can probably be optimized by allocating closure outside of the loop
            !!! But I can't figure out how it is done at the moment.
         nullify (closure)
         do point = 1, size(setPointID)
            PetscCall(DMPlexGetTransitiveClosure(dm, setPointID(point), PETSC_TRUE, PETSC_NULL_INTEGER, closure, ierr))
            do p = 1, size(closure), 2
               PetscCall(PetscSectionGetDoF(section, closure(p), numDof, ierr))
               if (numDof > 0) then
                  table(closure(p) + 1, :) = table(closure(p) + 1, :) .or. constraints
               end if ! numDof
            end do! p
            PetscCall(DMPlexRestoreTransitiveClosure(dm, setPointID(point), PETSC_TRUE, PETSC_NULL_INTEGER, closure, ierr))
         end do! cell
         PetscCall(ISRestoreIndices(setPointIS, setPointID, ierr))
      end if ! setPointIS
      PetscCall(ISDestroy(setPointIS, ierr))
   end subroutine MEF90SetupConstraintTableSet

#undef __FUNCT__
#define __FUNCT__ "MEF90SectionAllocateConstraint"
!!!
!!!
!!!  MEF90SectionAllocateConstraint: Associates dof to a section
!!!
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90SectionAllocateConstraint(dm, table, section, ierr)
      type(tDM), intent(IN)               :: dm
      PetscBool, dimension(:, :), pointer :: table
      type(tPetscSection), intent(INOUT)  :: section
      PetscErrorCode, intent(INOUT)       :: ierr

      PetscInt                            :: p, i, pStart, pEnd, numConstraints, numComponents
      PetscInt, dimension(:), pointer     :: constraints
      PetscInt                            :: field = 0

      PetscCall(DMPlexGetChart(dm, pStart, pEnd, ierr))
      do p = 1, pEnd
         numConstraints = count(table(p, :))
         if (numConstraints > 0) then
            PetscCall(PetscSectionSetConstraintDof(section, p - 1, numConstraints, ierr))
            PetscCall(PetscSectionSetFieldConstraintDof(section, p - 1, field, numConstraints, ierr))
         end if
      end do

      PetscCall(PetscSectionSetup(section, ierr))

      numComponents = size(table, 2)
      do p = 1, pEnd
         numConstraints = count(table(p, :))
         if (numConstraints > 0) then
            allocate (constraints(numConstraints))
            constraints = pack([(i - 1, i=1, numComponents)], table(p, :))
            PetscCall(PetscSectionSetConstraintIndices(section, p - 1, constraints, ierr))
            PetscCall(PetscSectionSetFieldConstraintIndices(section, p - 1, field, constraints, ierr))
            deallocate (constraints)
         end if
      end do
   end subroutine MEF90SectionAllocateConstraint

#undef __FUNCT__
#define __FUNCT__ "MEF90VecCopySF"
!!!
!!!
!!!  MEF90VecCopySF: rearrange a Vec according to the given SF
!!!
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!

   subroutine MEF90VecCopySF(vin, vout, sf, ierr)
      type(tVec), intent(IN)              :: vin
      type(tVec), intent(INOUT)           :: vout
      type(tPetscSF), intent(IN)          :: sf
      PetscErrorCode, intent(INOUT)       :: ierr

      PetscScalar, dimension(:), pointer   :: arrayin, arrayout

      PetscCall(VecGetArrayRead(vin, arrayin, ierr))
      PetscCall(VecGetArray(vout, arrayout, ierr))
      PetscCall(PetscSFBcastBegin(sf, MPIU_SCALAR, arrayin, arrayout, MPI_REPLACE, ierr))
      PetscCall(PetscSFBcastEnd(sf, MPIU_SCALAR, arrayin, arrayout, MPI_REPLACE, ierr))
      PetscCall(VecRestoreArrayRead(vin, arrayin, ierr))
      PetscCall(VecRestoreArray(vout, arrayout, ierr))
   end subroutine MEF90VecCopySF

#undef __FUNCT__
#define __FUNCT__ "MEF90IOSFCreate"
!!!
!!!
!!!  MEF90IOSFCreate: sf mapping between local and IO ordering and distribution
!!!
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!

   subroutine MEF90IOSFCreate(MEF90Ctx, v, liosf, iolsf, ierr)
      type(tVec), intent(IN)              :: v
      type(tPetscSF), intent(OUT)         :: liosf, iolsf
      type(MEF90Ctx_type), intent(IN)     :: MEF90Ctx
      PetscErrorCode, intent(INOUT)       :: ierr

      type(tPetscSF)                     :: naturalSF, ioSF, lcgSF, cglSF, tempSF, invTempSF
      type(tDM)                          :: dm
      character(len=PETSC_MAX_PATH_LEN)  :: vecname

      PetscCall(VecGetDM(v, dm, ierr))
      PetscCall(CreateLocalToCGlobalSF_Private(MEF90Ctx, dm, lcgSF, ierr))
      PetscCall(CreateCGlobalToLocalSF_Private(MEF90Ctx, dm, cglSF, ierr))
      PetscCall(PetscObjectGetName(v, vecname, ierr))
      if (MEF90Ctx%NumProcs > 1) then
         PetscCall(CreateNaturalToIOSF_Private(MEF90Ctx, dm, ioSF, ierr))
         PetscCall(DMGetNaturalSF(dm, naturalSF, ierr))
         PetscCall(PetscSFCompose(lcgSF, naturalSF, tempSF, ierr))
         PetscCall(PetscSFCompose(tempSF, ioSF, liosf, ierr))
         PetscCall(PetscSFSetUp(liosf, ierr))
         PetscCall(PetscSFDestroy(tempSF, ierr))
         PetscCall(PetscSFCompose(naturalSF, ioSF, tempSF, ierr))
         PetscCall(PetscSFCreateInverseSF(tempSF, invTempSF, ierr))
         PetscCall(PetscSFCompose(invTempSF, cglSF, iolsf, ierr))
         PetscCall(PetscSFSetUp(iolsf, ierr))
         PetscCall(PetscSFDestroy(ioSF, ierr))
         PetscCall(PetscSFDestroy(cglSF, ierr))
         PetscCall(PetscSFDestroy(lcgSF, ierr))
         PetscCall(PetscSFDestroy(tempSF, ierr))
         PetscCall(PetscSFDestroy(invTempSF, ierr))
      else
         liosf = lcgSF
         iolsf = cglSF
      end if
   end subroutine MEF90IOSFCreate

#undef __FUNCT__
#define __FUNCT__ "MEF90FaceSetIOSFCreate"
!!!
!!!
!!!  MEF90FaceSetIOSFCreate: sf mapping between local and IO ordering and distribution for
!!!                           Vec defined on Face Sets
!!!
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!

   subroutine MEF90FaceSetIOSFCreate(MEF90Ctx, v, liosf, iolsf, ierr)
      type(tVec), intent(IN)              :: v
      type(tPetscSF), intent(OUT)         :: liosf, iolsf
      type(MEF90Ctx_type), intent(IN)     :: MEF90Ctx
      PetscErrorCode, intent(INOUT)       :: ierr

      type(tPetscSF)                     :: naturalSF, ioSF, lcgSF, cglSF, tempSF, temp2SF, invTempSF, iosideSF, sideioSF
      type(tDM)                          :: dm
      character(len=PETSC_MAX_PATH_LEN)  :: vecname

      PetscCall(VecGetDM(v, dm, ierr))
      PetscCall(CreateLocalToCGlobalSF_Private(MEF90Ctx, dm, lcgSF, ierr))
      PetscCall(CreateCGlobalToLocalSF_Private(MEF90Ctx, dm, cglSF, ierr))
      PetscCall(CreateSideSF_Private(MEF90Ctx, dm, iosideSF, sideioSF, ierr))
      PetscCall(PetscObjectGetName(v, vecname, ierr))
      if (MEF90Ctx%NumProcs > 1) then
         PetscCall(CreateNaturalToIOSF_Private(MEF90Ctx, dm, ioSF, ierr))
         PetscCall(DMGetNaturalSF(dm, naturalSF, ierr))
         PetscCall(PetscSFCompose(lcgSF, naturalSF, tempSF, ierr))
         PetscCall(PetscSFCompose(tempSF, iosideSF, temp2SF, ierr))
         PetscCall(PetscSFCompose(temp2SF, ioSF, liosf, ierr))
         PetscCall(PetscSFDestroy(temp2SF, ierr))
         PetscCall(PetscSFSetUp(liosf, ierr))
         PetscCall(PetscSFDestroy(tempSF, ierr))
         PetscCall(PetscSFCompose(naturalSF, iosideSF, temp2SF, ierr))
         PetscCall(PetscSFCompose(temp2SF, ioSF, tempSF, ierr))
         PetscCall(PetscSFDestroy(temp2SF, ierr))
         PetscCall(PetscSFCreateInverseSF(tempSF, invTempSF, ierr))
         PetscCall(PetscSFCompose(invTempSF, cglSF, iolsf, ierr))
         PetscCall(PetscSFSetUp(iolsf, ierr))
         PetscCall(PetscSFDestroy(ioSF, ierr))
         PetscCall(PetscSFDestroy(cglSF, ierr))
         PetscCall(PetscSFDestroy(lcgSF, ierr))
         PetscCall(PetscSFDestroy(tempSF, ierr))
         PetscCall(PetscSFDestroy(invTempSF, ierr))
      else
         PetscCall(PetscSFCompose(lcgSF, iosideSF, liosf, ierr))
         PetscCall(PetscSFCompose(sideioSF, cglSF, iolsf, ierr))
         PetscCall(PetscSFDestroy(lcgSF, ierr))
         PetscCall(PetscSFDestroy(cglSF, ierr))
      end if
      PetscCall(PetscSFDestroy(iosideSF, ierr))
      PetscCall(PetscSFDestroy(sideioSF, ierr))
   end subroutine MEF90FaceSetIOSFCreate

#undef __FUNCT__
#define __FUNCT__ "MEF90ConstraintSFCreate"
!!!
!!!
!!!  MEF90ConstraintSFCreate: sf mapping between local and constraint Vec ordering and distribution
!!!
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!                Blaise Bourdin  bourdin@mcmaster.ca
!!!

   subroutine MEF90ConstraintSFCreate(MEF90Ctx, v, vB, sf, invSF, ierr)
      type(tVec), intent(IN)                   :: v, vB
      type(tPetscSF), intent(OUT)              :: sf, invSF
      type(MEF90Ctx_type), intent(IN)          :: MEF90Ctx
      PetscErrorCode, intent(INOUT)            :: ierr

      type(tDM)                               :: dm, dmB
      type(tPetscSection)                     :: locSection, locBSection
      type(sPetscSFNode), dimension(:), pointer :: remote
      PetscInt, dimension(:), pointer           :: local, cindices
      PetscInt                                :: pStart, pEnd, p, d, nleaves, ldof, loff, cdof, coff, nsize, nroots

      nleaves = 0
      nsize = 0
      PetscCall(VecGetDM(v, dm, ierr))
      PetscCall(VecGetDM(vB, dmB, ierr))
      PetscCall(DMGetLocalSection(dm, locSection, ierr))
      PetscCall(PetscSectionGetStorageSize(locSection, nroots, ierr))
      PetscCall(DMGetLocalSection(dmB, locBSection, ierr))
      PetscCall(PetscSectionGetChart(locBSection, pStart, pEnd, ierr))
      do p = pStart, pEnd - 1
         PetscCall(PetscSectionGetConstraintDof(locBSection, p, cdof, ierr))
         nleaves = nleaves + cdof
      end do
      allocate (remote(nleaves))
      allocate (local(nleaves))
      do p = pStart, pEnd - 1
         PetscCall(PetscSectionGetDoF(locSection, p, ldof, ierr))
         PetscCall(PetscSectionGetOffset(locSection, p, loff, ierr))
         PetscCall(PetscSectionGetConstraintDof(locBSection, p, cdof, ierr))
         if (cdof > 0) then
            PetscCall(PetscSectionGetConstraintIndices(locBSection, p, cindices, ierr))
            PetscCall(PetscSectionGetOffset(locBSection, p, coff, ierr))
            if (coff >= 0) then
               do d = 1, cdof
                  local(nsize + 1) = coff + cindices(d)
                  remote(nsize + 1)%rank = MEF90Ctx%rank
                  remote(nsize + 1)%index = loff + cindices(d)
                  nsize = nsize + 1
               end do ! d
            end if ! coff
            PetscCall(PetscSectionRestoreConstraintIndices(locBSection, p, cindices, ierr))
         end if ! cdof
      end do
      PetscCall(PetscSFCreate(MEF90Ctx%Comm, sf, ierr))
      PetscCall(PetscSFSetFromOptions(sf, ierr))
      PetscCall(PetscSFSetGraph(sf, nroots, nleaves, local, PETSC_COPY_VALUES, remote, PETSC_COPY_VALUES, ierr))
      deallocate (remote)
      deallocate (local)
      PetscCall(PetscSFSetUp(sf, ierr))
      PetscCall(PetscSFCreateInverseSF(sf, invSF, ierr))
      PetscCall(PetscSFSetUp(invSF, ierr))
   end subroutine MEF90ConstraintSFCreate

#undef __FUNCT__
#define __FUNCT__ "MEF90VecGlobalToLocalConstraint"
!!!
!!!
!!!  MEF90VecGlobalToLocalConstraint: do a VecGlobaltoLocal then copy constrained values
!!!
!!!  (c) 2022      Blaise Bourdin  bourdin@mcmaster.ca
!!!

   subroutine MEF90VecGlobalToLocalConstraint(g, c, l, ierr)
      type(tVec), intent(IN)                   :: g, c
      type(tVec), intent(INOUT)                :: l
      PetscErrorCode, intent(INOUT)            :: ierr

      type(tDM)                               :: dm
      type(tPetscSection)                     :: s
      PetscInt                                :: numConstraint
      PetscInt                                :: p, pStart, pEnd
      PetscReal, dimension(:), pointer          :: vArray

      PetscCall(VecGetDM(g, dm, ierr))
      PetscCall(DMGlobalToLocal(dm, g, INSERT_VALUES, l, ierr))
      if (.not. PetscObjectIsNull(c)) then
         PetscCall(DMGetLocalSection(dm, s, ierr))
         PetscCall(PetscSectionGetChart(s, pStart, pEnd, ierr))
         do p = pStart, pEnd - 1
            PetscCall(PetscSectionGetConstraintDof(s, p, numConstraint, ierr))
            if (numConstraint > 0) then
               PetscCall(VecGetValuesSection(c, s, p, vArray, ierr))
               PetscCall(VecSetValuesSection(l, s, p, vArray, INSERT_ALL_VALUES, ierr))
               PetscCall(VecRestoreValuesSection(c, s, p, vArray, ierr))
            end if ! numConstraint
         end do ! p
      end if
   end subroutine MEF90VecGlobalToLocalConstraint

#undef __FUNCT__
#define __FUNCT__ "MEF90VecCreateIO"
!!!
!!!
!!!  MEF90VecCreateIO: create IO Vec
!!!
!!!  (c) 2022      Alexis Marboeuf  marboeua@mcmaster.ca
!!!

   subroutine MEF90VecCreateIO(v, bs, sf, ierr)
      type(tPetscSF), intent(IN)               :: sf
      type(tVec), intent(INOUT)                :: v
      PetscInt, intent(IN)                     :: bs
      PetscErrorCode, intent(INOUT)            :: ierr

      MPI_Comm                                :: comm
      PetscInt                                :: nroots, nleaves
      type(sPetscSFNode), dimension(:), pointer :: iremote
      PetscInt, dimension(:), pointer           :: ilocal

      PetscCall(PetscSFGetGraph(sf, nroots, nleaves, ilocal, iremote, ierr))
      PetscCall(PetscObjectGetComm(sf, comm, ierr))
      PetscCall(VecCreateMPI(comm, nleaves, PETSC_DETERMINE, v, ierr))
      PetscCall(VecSetBlockSize(v, bs, ierr))
      PetscCall(PetscSFrestoreGraph(sf, nroots, nleaves, ilocal, iremote, ierr))
      if (associated(iLocal)) then
         deallocate (ilocal)
      end if
   end subroutine MEF90VecCreateIO

#undef __FUNCT__
#define __FUNCT__ "MEF90VecSetValuesFromOptions"
!!!
!!!
!!!  MEF90VecSetValuesFromOptions: Fill values of a Vec using command line options
!!!
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90VecSetValuesFromOptions(v, scalingFactor, ierr)
      type(tVec), intent(INOUT)                :: v
      PetscReal, intent(IN)                    :: scalingFactor
      PetscErrorCode, intent(INOUT)            :: ierr

      type(tDM)                               :: dm
      PetscEnum                               :: setType
      PetscInt                                :: set, point
      type(tIS)                               :: setIS, pointIS
      PetscInt, dimension(:), pointer           :: setID, pointID
      character(len=MEF90MXSTRLEN)            :: ValueKey, name
      PetscBool                               :: flg
      PetscInt                                :: dim, numOpt, bs, numDofClosure, i
      PetscReal, dimension(:), pointer          :: Val, vArray
      type(tPetscSection)                     :: section

      PetscCall(VecGetDM(v, dm, ierr))
      PetscCall(PetscObjectGetName(v, name, ierr))
      PetscCall(DMGetLocalSection(dm, section, ierr))

      PetscCall(DMGetDimension(dm, dim, ierr))
      PetscCall(VecGetBlockSize(v, bs, ierr))
      allocate (Val(bs))

      do setType = 1, size(MEF90SetType)
         PetscCall(DMGetLabelIdIS(dm, MEF90SetLabelName(setType), setIS, ierr))
         ! PetscCall(MEF90ISAllGatherMerge(comm,setIS,ierr))
         if (.not. PetscObjectIsNull(setIS)) then
            PetscCall(ISGetIndices(setIS, setID, ierr))
            do set = 1, size(setID)
               write (ValueKey, '("-",a2,I4.4,"_",a)') MEF90SetPrefix(setType), setID(set), trim(name)
               numOpt = bs
               PetscCall(PetscOptionsGetRealArray(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, trim(ValueKey), Val, numOpt, flg, ierr))
               if (numOpt > 0) then
                  PetscCall(DMGetStratumIS(dm, MEF90SetLabelName(setType), setID(set), pointIS, ierr))
                        !!! Set the values on the closure of the current point
                  if (.not. PetscObjectIsNull(pointIS)) then
                     PetscCall(ISGetIndices(pointIS, pointID, ierr))
                     do point = 1, size(pointID)
                        PetscCall(MEF90VecGetClosureSize(v, pointID(point), numDofClosure, ierr))
                        if (numDofClosure > 0) then
                           PetscCall(DMPlexVecGetClosure(dm, section, v, pointID(point), PETSC_NULL_INTEGER, vArray, ierr))
                           do i = 1, numDofClosure / bs
                              vArray((i - 1) * bs + 1:i * bs) = scalingFactor * Val
                           end do
                           PetscCall(DMPlexVecSetClosure(dm, section, v, pointID(point), vArray, INSERT_ALL_VALUES, ierr))
                           PetscCall(DMPlexVecRestoreClosure(dm, section, v, pointID(point), PETSC_NULL_INTEGER, vArray, ierr))
                        end if ! numDofClosure
                     end do ! point
                     PetscCall(ISRestoreIndices(pointIS, pointID, ierr))
                  end if ! pointIS
                  PetscCall(ISDestroy(pointIS, ierr))
               end if ! numOpt
            end do ! set
            PetscCall(ISRestoreIndices(setIS, setID, ierr))
         end if ! setIS
         PetscCall(ISDestroy(setIS, ierr))
      end do ! setType
      deallocate (Val)
   end subroutine MEF90VecSetValuesFromOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90VecSetBCValuesFromOptions"
!!!
!!!
!!!  MEF90VecSetBCValuesFromOptions: Fill boundary values of a Vec using command line options
!!!
!!!  (c) 2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90VecSetBCValuesFromOptions(v, scalingFactor, ierr)
      type(tVec), intent(INOUT)                :: v
      PetscReal, intent(IN)                    :: scalingFactor
      PetscErrorCode, intent(INOUT)            :: ierr

      type(tDM)                               :: dm
      PetscEnum                               :: setType
      PetscInt                                :: set, point, c
      type(tIS)                               :: setIS, pointIS
      PetscInt, dimension(:), pointer           :: setID, pointID
      character(len=MEF90MXSTRLEN)            :: BCOptionKey, BCValueKey, name
      PetscBool, dimension(:), pointer          :: setBC
      PetscBool                               :: flg
      PetscInt                                :: dim, numBC, bs, numDofClosure
      PetscReal, dimension(:), pointer          :: BCVal, vArray
      type(tPetscSection)                     :: section

      PetscCall(VecGetDM(v, dm, ierr))
      PetscCall(PetscObjectGetName(v, name, ierr))
      PetscCall(DMGetLocalSection(dm, section, ierr))

      PetscCall(DMGetDimension(dm, dim, ierr))
      PetscCall(VecGetBlockSize(v, bs, ierr))
      allocate (setBC(bs))
      allocate (BCVal(bs))

      do setType = 1, size(MEF90SetType)
         PetscCall(DMGetLabelIdIS(dm, MEF90SetLabelName(setType), setIS, ierr))
         ! PetscCall(MEF90ISAllGatherMerge(comm,setIS,ierr))
         if (.not. PetscObjectIsNull(setIS)) then
            PetscCall(ISGetIndices(setIS, setID, ierr))
            do set = 1, size(setID)
               setBC = .false.
               write (BCOptionKey, '("-",a2,I4.4,"_",a,"BC")') MEF90SetPrefix(setType), setID(set), trim(name)
               numBC = bs
               PetscCall(PetscOptionsGetBoolArray(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, trim(BCOptionKey), setBC, numBC, flg, ierr))
               if (any(setBC)) then
                        !!! At least 1 dof has a boundary condition
                        !!! Get the unit BC value on the set
                  write (BCValueKey, '("-",a2,I4.4,"_Boundary",a)') MEF90SetPrefix(setType), setID(set), trim(name)
                  numBC = bs
                  PetscCall(PetscOptionsGetRealArray(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, trim(BCValueKey), BCVal, numBC, flg, ierr))
                  PetscCall(DMGetStratumIS(dm, MEF90SetLabelName(setType), setID(set), pointIS, ierr))
                        !!! Set the boundary values on the closure of the current point
                  if (.not. PetscObjectIsNull(pointIS)) then
                     PetscCall(ISGetIndices(pointIS, pointID, ierr))
                     do point = 1, size(pointID)
                        PetscCall(MEF90VecGetClosureSize(v, pointID(point), numDofClosure, ierr))
                        if (numDofClosure > 0) then
                           PetscCall(DMPlexVecGetClosure(dm, section, v, pointID(point), PETSC_NULL_INTEGER, vArray, ierr))
                           do c = 1, bs
                              if (setBC(c)) then
                                 vArray(c::bs) = scalingFactor * BCVal(c)
                              end if ! setBC
                           end do ! c
                           PetscCall(DMPlexVecSetClosure(dm, section, v, pointID(point), vArray, INSERT_ALL_VALUES, ierr))
                           PetscCall(DMPlexVecRestoreClosure(dm, section, v, pointID(point), PETSC_NULL_INTEGER, vArray, ierr))
                        end if
                     end do ! point
                     PetscCall(ISRestoreIndices(pointIS, pointID, ierr))
                  end if ! pointIS
                  PetscCall(ISDestroy(pointIS, ierr))
               end if ! setBC
            end do ! set
            PetscCall(ISRestoreIndices(setIS, setID, ierr))
         end if ! setIS
         PetscCall(ISDestroy(setIS, ierr))
      end do ! setType
      deallocate (setBC)
      deallocate (BCVal)
   end subroutine MEF90VecSetBCValuesFromOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90VecSetValuesFromOptionsExpr"
!!!
!!!
!!!  MEF90VecSetValuesFromOptionsExpr: Fill values of a Vec using expressions passed as petsc options
!!!
!!!  (c) 2025      Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90VecSetValuesFromOptionsExpr(v, t, ierr)
      type(tVec), intent(INOUT)                 :: v
      PetscReal, intent(IN)                     :: t
      PetscErrorCode, intent(INOUT)             :: ierr

#ifdef MEF90_HAVE_SYMENGINEF90
      type(tDM)                                :: dm
      PetscEnum                                :: setType
      PetscInt                                 :: set, point, p
      type(tIS)                                :: setIS, pointIS
      PetscInt, dimension(:), pointer            :: setID, pointID
      character(len=MEF90MXSTRLEN)             :: ValueKey, name, ExprStr, IOBuffer
      PetscBool                                :: flg
      PetscInt                                 :: dim, numOpt, bs, numDofClosure, numDof, i, c, dof
      PetscInt, dimension(:), pointer            :: closure
      PetscReal, dimension(:), pointer           :: vArray
      type(tPetscSection)                      :: section
      type(tPetscSection)                      :: coordSection
      type(tVec)                               :: coordVec
      PetscReal, dimension(:), pointer           :: coordArray
      PetscReal, dimension(3)                   :: xyz
      type(Basic), dimension(:), allocatable     :: exprs
      type(Basic)                              :: tmpExpr
      type(symbol), dimension(3)                :: vars
      character                                :: delim = ';'
      character(len=MEF90MXSTRLEN), dimension(:), allocatable :: ExprStrComp

      PetscCall(VecGetDM(v, dm, ierr))
      PetscCall(PetscObjectGetName(v, name, ierr))
      PetscCall(DMGetLocalSection(dm, section, ierr))

      PetscCall(DMGetDimension(dm, dim, ierr))
      PetscCall(VecGetBlockSize(v, bs, ierr))

      vars = [Symbol("x"), Symbol("y"), Symbol("z")]
      allocate (exprs(bs))

      PetscCall(DMGetCoordinateSection(dm, coordSection, ierr))
      PetscCall(DMGetCoordinatesLocal(dm, coordVec, ierr))

      do setType = 1, size(MEF90SetType)
         PetscCall(DMGetLabelIdIS(dm, MEF90SetLabelName(setType), setIS, ierr))
         ! PetscCall(MEF90ISAllGatherMerge(comm,setIS,ierr))
         if (.not. PetscObjectIsNull(setIS)) then
            PetscCall(ISGetIndices(setIS, setID, ierr))
            do set = 1, size(setID)
               write (ValueKey, '("-",a2,I4.4,"_",a,"Expr")') MEF90SetPrefix(setType), setID(set), trim(name)
               PetscCall(PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, trim(ValueKey), ExprStr, flg, ierr))
               if (flg) then
                  call MEF90StrTokenize(ExprStr, delim, ExprStrComp)
                  numOpt = size(ExprStrComp)
                  if (numOpt < bs) then
                     write (IOBuffer, "(A,' was expecting ',I2,' expressions but got ',I2,' for key ',A)") __FUNCT__, bs, size(ExprStrComp), trim(ValueKey)
                     SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, IOBuffer)
                  end if
                        !!! create parsers
                  do c = 1, bs
                     exprs(c) = parse(ExprStrComp(c))
                     exprs(c) = exprs(c)%subs(Symbol("t"), RealDouble(t))
                  end do
                  PetscCall(DMGetStratumIS(dm, MEF90SetLabelName(setType), setID(set), pointIS, ierr))
                        !!! Set the values on the closure of the current point
                  if (.not. PetscObjectIsNull(pointIS)) then
                     PetscCall(ISGetIndices(pointIS, pointID, ierr))
                     do point = 1, size(pointID)
                        PetscCall(DMPlexGetTransitiveClosure(dm, pointID(point), PETSC_TRUE, numDofClosure, closure, ierr))
                        if (numDofClosure > 0) then
                           PetscCall(DMPlexVecGetClosure(dm, section, v, pointID(point), PETSC_NULL_INTEGER, vArray, ierr))
                           dof = 0
                           do p = 1, size(closure), 2
                              PetscCall(PetscSectionGetDof(section, closure(p), numDof, ierr))
                              if (numDof > 0) then
                                 dof = dof + 1
                                            !!! Get the coordinates of the dof associated with the point
                                            !!! trick: the coordinate of a point is the average of the coordinates of the points in its closure
                                 PetscCall(DMPlexVecGetClosure(dm, coordSection, coordVec, closure(p), PETSC_NULL_INTEGER, coordArray, ierr))
                                 do i = 1, dim
                                    xyz(i) = sum(coordArray(i:size(coordArray):dim)) * dim / size(coordArray)
                                 end do
                                 PetscCall(DMPlexVecRestoreClosure(dm, coordSection, coordVec, closure(p), PETSC_NULL_INTEGER, coordArray, ierr))
                                 do c = 1, bs
                                    tmpExpr = Exprs(c)
                                    do i = 1, dim
                                       tmpExpr = tmpExpr%subs(vars(i), RealDouble(xyz(i)))
                                    end do ! i
                                    tmpExpr = tmpExpr%evalf()
                                    vArray((dof - 1) * bs + c) = tmpExpr%dbl()
                                 end do ! c
                              end if !bnumDof
                           end do !p
                           PetscCall(DMPlexVecSetClosure(dm, section, v, pointID(point), vArray, INSERT_ALL_VALUES, ierr))
                           PetscCall(DMPlexVecRestoreClosure(dm, section, v, pointID(point), PETSC_NULL_INTEGER, vArray, ierr))
                        end if ! numDofClosure
                        PetscCall(DMPlexRestoreTransitiveClosure(dm, pointID(point), PETSC_TRUE, numDofClosure, closure, ierr))
                     end do ! point
                     PetscCall(ISRestoreIndices(pointIS, pointID, ierr))
                  end if ! pointIS
                  PetscCall(ISDestroy(pointIS, ierr))
                  deallocate (ExprStrComp)
               end if ! flg
            end do ! set
            PetscCall(ISRestoreIndices(setIS, setID, ierr))
         end if ! setIS
         PetscCall(ISDestroy(setIS, ierr))
      end do ! setType
      deallocate (exprs)
#else
      write (*, *) "ERROR: ", __FUNCT__, " requires symengine-f90 support"
      stop
#endif
   end subroutine MEF90VecSetValuesFromOptionsExpr

#undef __FUNCT__
#define __FUNCT__ "MEF90VecSetBCValuesFromOptionsExpr"
!!!
!!!
!!!  MEF90VecSetBCValuesFromOptionsExpr: Fill boundary values of a Vec using expressions passed as petsc options
!!!
!!!  (c) 2025      Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90VecSetBCValuesFromOptionsExpr(v, t, ierr)
      type(tVec), intent(INOUT)                 :: v
      PetscReal, intent(IN)                     :: t
      PetscErrorCode, intent(INOUT)             :: ierr

#ifdef MEF90_HAVE_SYMENGINEF90
      type(tDM)                                :: dm
      PetscEnum                                :: setType
      PetscInt                                 :: set, point, c, p, i
      type(tIS)                                :: setIS, pointIS
      PetscInt, dimension(:), pointer            :: setID, pointID
      character(len=MEF90MXSTRLEN)             :: BCOptionKey, BCValueKey, name, ExprStr, IOBuffer
      PetscBool, dimension(:), pointer           :: setBC
      PetscBool                                :: flg
      PetscInt                                 :: dim, numOpt, numBC, bs, numDofClosure, numDof, dof
      PetscInt, dimension(:), pointer            :: closure
      PetscReal, dimension(:), pointer           :: vArray
      type(tPetscSection)                      :: section
      type(tPetscSection)                      :: coordSection
      type(tVec)                               :: coordVec
      PetscReal, dimension(:), pointer           :: coordArray
      PetscReal, dimension(3)                   :: xyz
      type(Basic), dimension(:), allocatable     :: exprs
      type(Basic)                              :: tmpExpr
      type(symbol), dimension(3)                :: vars
      character                                :: delim = ';'
      character(len=MEF90MXSTRLEN), dimension(:), allocatable :: ExprStrComp

      PetscCall(VecGetDM(v, dm, ierr))
      PetscCall(PetscObjectGetName(v, name, ierr))
      PetscCall(DMGetLocalSection(dm, section, ierr))

      PetscCall(DMGetDimension(dm, dim, ierr))
      PetscCall(VecGetBlockSize(v, bs, ierr))
      allocate (setBC(bs))

      vars = [Symbol("x"), Symbol("y"), Symbol("z")]
      allocate (exprs(bs))

      PetscCall(DMGetCoordinateSection(dm, coordSection, ierr))
      PetscCall(DMGetCoordinatesLocal(dm, coordVec, ierr))

      do setType = 1, size(MEF90SetType)
         PetscCall(DMGetLabelIdIS(dm, MEF90SetLabelName(setType), setIS, ierr))
         ! PetscCall(MEF90ISAllGatherMerge(comm,setIS,ierr))
         if (.not. PetscObjectIsNull(setIS)) then
            PetscCall(ISGetIndices(setIS, setID, ierr))
            do set = 1, size(setID)
               setBC = .false.
               write (BCOptionKey, '("-",a2,I4.4,"_",a,"BC")') MEF90SetPrefix(setType), setID(set), trim(name)
               numBC = bs
               PetscCall(PetscOptionsGetBoolArray(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, trim(BCOptionKey), setBC, numBC, flg, ierr))
               if (any(setBC)) then
                        !!! At least 1 dof has a boundary condition
                        !!! Get the unit BC value on the set
                  write (BCValueKey, '("-",a2,I4.4,"_Boundary",a,"Expr")') MEF90SetPrefix(setType), setID(set), trim(name)
                  PetscCall(PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, trim(BCValueKey), ExprStr, flg, ierr))
                  numOpt = MEF90StrCount(ExprStr, delim)
                  call MEF90StrTokenize(ExprStr, delim, ExprStrComp)
                  if (size(ExprStrComp) < bs) then
                     write (IOBuffer, "(A,' was expecting ',I2,' expressions but got ',I2,' for key ',A)") __FUNCT__, bs, size(ExprStrComp), trim(BCValueKey)
                     SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, IOBuffer)
                  end if
                        !!! create parsers
                  do c = 1, bs
                     exprs(c) = parse(ExprStrComp(c))
                     exprs(c) = exprs(c)%subs(Symbol("t"), RealDouble(t))
                  end do

                  PetscCall(DMGetStratumIS(dm, MEF90SetLabelName(setType), setID(set), pointIS, ierr))
                        !!! Set the boundary values on the closure of the current point
                  if (.not. PetscObjectIsNull(pointIS)) then
                     PetscCall(ISGetIndices(pointIS, pointID, ierr))
                     do point = 1, size(pointID)
                        PetscCall(DMPlexGetTransitiveClosure(dm, pointID(point), PETSC_TRUE, numDofClosure, closure, ierr))
                        if (numDofClosure > 0) then
                           PetscCall(DMPlexVecGetClosure(dm, section, v, pointID(point), PETSC_NULL_INTEGER, vArray, ierr))
                           dof = 0
                           do p = 1, size(closure), 2
                              PetscCall(PetscSectionGetDof(section, closure(p), numDof, ierr))
                              if (numDof > 0) then
                                 dof = dof + 1
                                            !!! Get the coordinates of the dof associated with the point
                                            !!! trick: the coordinate of a point is the average of the coordinates of the points in its closure
                                 PetscCall(DMPlexVecGetClosure(dm, coordSection, coordVec, closure(p), PETSC_NULL_INTEGER, coordArray, ierr))
                                 do i = 1, dim
                                    xyz(i) = sum(coordArray(i:size(coordArray):dim)) * dim / size(coordArray)
                                 end do
                                 PetscCall(DMPlexVecRestoreClosure(dm, coordSection, coordVec, closure(p), PETSC_NULL_INTEGER, coordArray, ierr))
                                 do c = 1, bs
                                    if (setBC(c)) then
                                       tmpExpr = Exprs(c)
                                       do i = 1, dim
                                          tmpExpr = tmpExpr%subs(vars(i), RealDouble(xyz(i)))
                                       end do ! i
                                       tmpExpr = tmpExpr%evalf()
                                       vArray((dof - 1) * bs + c) = tmpExpr%dbl()
                                    end if
                                 end do ! c
                              end if !bnumDof
                           end do !p
                           PetscCall(DMPlexVecSetClosure(dm, section, v, pointID(point), vArray, INSERT_ALL_VALUES, ierr))
                           PetscCall(DMPlexVecRestoreClosure(dm, section, v, pointID(point), PETSC_NULL_INTEGER, vArray, ierr))
                        end if ! numDofClosure
                        PetscCall(DMPlexRestoreTransitiveClosure(dm, pointID(point), PETSC_TRUE, numDofClosure, closure, ierr))
                     end do ! point
                     PetscCall(ISRestoreIndices(pointIS, pointID, ierr))
                  end if ! pointIS
                  PetscCall(ISDestroy(pointIS, ierr))
                  deallocate (ExprStrComp)
               end if ! setBC
            end do ! set
            PetscCall(ISRestoreIndices(setIS, setID, ierr))
         end if ! setIS
         PetscCall(ISDestroy(setIS, ierr))
      end do ! setType
      deallocate (setBC)
      deallocate (exprs)
#else
      write (*, *) "ERROR: ", __FUNCT__, " requires symengine-f90 support"
      stop
#endif
   end subroutine MEF90VecSetBCValuesFromOptionsExpr

!!! Private functions below
#undef __FUNCT__
#define __FUNCT__ "CreateNaturalToIOSF_Private"
!!!
!!!
!!!  CreateNaturalToIOSF_Private:
!!!
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!

   subroutine CreateNaturalToIOSF_Private(MEF90Ctx, dm, sf, ierr)
      type(tDM), intent(IN)                    :: dm
      type(MEF90Ctx_type), intent(IN)          :: MEF90Ctx
      type(tPetscSF), intent(OUT)              :: sf
      PetscErrorCode, intent(INOUT)            :: ierr

      type(tVec)                              :: vnat, vio
      type(tPetscLayout)                      :: ioMap, natMap
      type(sPetscSFNode), dimension(:), pointer :: remote
      PetscInt, dimension(:), pointer           :: ioRange
      PetscInt                                :: nroots, nleaves, globalIndex, i, globalSize, bs

      PetscCall(DMPlexCreateNaturalVector(dm, vnat, ierr))
      PetscCall(VecGetSize(vnat, globalSize, ierr))
      PetscCall(VecGetBlockSize(vnat, bs, ierr))
      PetscCall(VecCreate(MEF90Ctx%Comm, vio, ierr))
      PetscCall(VecSetBlockSize(vio, bs, ierr))
      PetscCall(VecSetSizes(vio, PETSC_DETERMINE, globalSize, ierr))
      PetscCall(VecSetFromOptions(vio, ierr))
      PetscCall(VecGetLayout(vio, ioMap, ierr))
      PetscCall(VecGetLayout(vnat, natMap, ierr))
      PetscCall(PetscLayoutGetLocalSize(natMap, nroots, ierr))
      PetscCall(PetscLayoutGetLocalSize(ioMap, nleaves, ierr))
      PetscCall(PetscLayoutGetRanges(ioMap, ioRange, ierr))
      allocate (remote(nleaves))
      do i = 0, nleaves - 1
         globalIndex = ioRange(MEF90Ctx%rank + 1) + i
         remote(i + 1)%rank = 0
         remote(i + 1)%index = globalIndex
      end do
      PetscCall(PetscSFCreate(MEF90Ctx%Comm, sf, ierr))
      ! PetscCall(PetscObjectSetName(sf,"Natural-To-IO SF",ierr))
      PetscCall(PetscSFSetFromOptions(sf, ierr))
      PetscCall(PetscSFSetGraph(sf, nroots, nleaves, PETSC_NULL_INTEGER_ARRAY, PETSC_COPY_VALUES, remote, PETSC_COPY_VALUES, ierr))
      PetscCall(PetscSFSetUp(sf, ierr))
      ! PetscCall(PetscSFViewFromOptions(sf,PETSC_NULL_OBJECT,"-naturaltoio_sf_view",ierr))
      PetscCall(VecDestroy(vio, ierr))
      PetscCall(VecDestroy(vnat, ierr))
      deallocate (remote)
   end subroutine CreateNaturalToIOSF_Private

#undef __FUNCT__
#define __FUNCT__ "CreateLocalToCGlobalSF_Private"
!!!
!!!
!!!  CreateLocalToCGlobalSF_Private:
!!!
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!

   subroutine CreateLocalToCGlobalSF_Private(MEF90Ctx, dm, sf, ierr)
      type(tDM), intent(IN)                    :: dm
      type(MEF90Ctx_type), intent(IN)          :: MEF90Ctx
      type(tPetscSF), intent(OUT)              :: sf
      PetscErrorCode, intent(INOUT)            :: ierr

      type(tPetscSection)                     :: locSection, gSection
      type(tPetscSF)                          :: overlapSF, idSF
      type(sPetscSFNode), dimension(:), pointer :: remote
      PetscInt, dimension(:), pointer           :: remoteOffsets
      PetscInt                                :: pStart, pEnd, p, n

      PetscCall(DMGetLocalSection(dm, locSection, ierr))
      PetscCall(DMGetPointSF(dm, overlapSF, ierr))
      PetscCall(PetscSectionCreateGlobalSection(locSection, overlapSF, PETSC_TRUE, PETSC_TRUE, PETSC_TRUE, gSection, ierr))
      PetscCall(PetscSectionGetChart(locSection, pStart, pEnd, ierr))
      n = pEnd - pStart
      allocate (remote(n))
      do p = 1, n
         remote(p)%rank = MEF90Ctx%rank
         remote(p)%index = p - 1
      end do
      PetscCall(PetscSFCreate(MEF90Ctx%Comm, idSF, ierr))
      PetscCall(PetscSFSetFromOptions(idSF, ierr))
      PetscCall(PetscSFSetGraph(idSF, n, n, PETSC_NULL_INTEGER_ARRAY, PETSC_COPY_VALUES, remote, PETSC_COPY_VALUES, ierr))
      PetscCall(PetscSFSetUp(idSF, ierr))
      PetscCall(PetscSFCreateRemoteOffsets(idSF, locSection, gSection, remoteOffsets, ierr))
      PetscCall(PetscSFCreateSectionSF(idSF, locSection, remoteOffsets, gSection, sf, ierr))
      if (associated(remoteOffsets)) then
         PetscCall(PetscSFDestroyRemoteOffsets(remoteOffsets, ierr))
      end if
      PetscCall(PetscSFSetUp(sf, ierr))
      ! PetscCall(PetscObjectSetName(sf,"Local-To-CGlobal SF",ierr))
      PetscCall(PetscSFSetFromOptions(sf, ierr))
      ! PetscCall(PetscSFViewFromOptions(sf,PETSC_NULL_OBJECT,"-localtocglobal_sf_view",ierr))
      PetscCall(PetscSectionDestroy(gSection, ierr))
      PetscCall(PetscSFDestroy(idSF, ierr))
      deallocate (remote)
   end subroutine CreateLocalToCGlobalSF_Private

#undef __FUNCT__
#define __FUNCT__ "CreateCGlobalToLocalSF_Private"
!!!
!!!
!!!  CreateCGlobalToLocalSF_Private:
!!!
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!

   subroutine CreateCGlobalToLocalSF_Private(MEF90Ctx, dm, sf, ierr)
      type(tDM), intent(IN)                    :: dm
      type(MEF90Ctx_type), intent(IN)          :: MEF90Ctx
      type(tPetscSF), intent(OUT)              :: sf
      PetscErrorCode, intent(INOUT)            :: ierr

      type(tPetscSection)                     :: locSection, gSection
      type(tPetscSF)                          :: overlapSF, idSF, tempSF, ttempSF
      type(sPetscSFNode), dimension(:), pointer :: remote, tempRemote, lgRemote, glRemote
      PetscInt, dimension(:), pointer           :: tempLocal, lgLocal, glLocal, remoteOffsets
      PetscInt                                :: pStart, pEnd, p, n, lgNRoots, lgNLeaves, tempNRoots, tempNLeaves, glNRoots, glNLeaves

      PetscCall(DMGetLocalSection(dm, locSection, ierr))
      PetscCall(DMGetPointSF(dm, overlapSF, ierr))
      PetscCall(PetscSectionCreateGlobalSection(locSection, overlapSF, PETSC_TRUE, PETSC_TRUE, PETSC_TRUE, gSection, ierr))
      PetscCall(PetscSectionGetChart(locSection, pStart, pEnd, ierr))
      n = pEnd - pStart
      allocate (remote(n))
      do p = 1, n
         remote(p)%rank = MEF90Ctx%rank
         remote(p)%index = p - 1
      end do
      PetscCall(PetscSFCreate(MEF90Ctx%Comm, idSF, ierr))
      PetscCall(PetscSFSetFromOptions(idSF, ierr))
      PetscCall(PetscSFSetGraph(idSF, n, n, PETSC_NULL_INTEGER_ARRAY, PETSC_COPY_VALUES, remote, PETSC_COPY_VALUES, ierr))
      PetscCall(PetscSFSetUp(idSF, ierr))
      PetscCall(PetscSFCreateRemoteOffsets(idSF, locSection, gSection, remoteOffsets, ierr))
      PetscCall(PetscSFCreateSectionSF(idSF, locSection, remoteOffsets, gSection, tempSF, ierr))
      PetscCall(PetscSFCreateInverseSF(tempSF, sf, ierr))
      PetscCall(PetscSFDestroyRemoteOffsets(remoteOffsets, ierr))
      PetscCall(PetscSFDestroy(tempSF, ierr))
      if (MEF90Ctx%NumProcs > 1) then
         PetscCall(PetscSFGetGraph(sf, lgNRoots, lgNLeaves, lgLocal, lgRemote, ierr))
         PetscCall(PetscSFCreateRemoteOffsets(overlapSF, locSection, gSection, remoteOffsets, ierr))
         PetscCall(PetscSFCreateSectionSF(overlapSF, locSection, remoteOffsets, gSection, ttempSF, ierr))
         PetscCall(PetscSFCreateInverseSF(ttempSF, tempSF, ierr))
         PetscCall(PetscSFDestroyRemoteOffsets(remoteOffsets, ierr))
         PetscCall(PetscSFDestroy(ttempSF, ierr))
         PetscCall(PetscSFGetGraph(tempSF, tempNRoots, tempNLeaves, tempLocal, tempRemote, ierr))
         glNRoots = lgNRoots
         glNLeaves = lgNLeaves + tempNLeaves
         allocate (glLocal(glNLeaves))
         allocate (glRemote(glNLeaves))
         if (loc(lgLocal) /= loc(PETSC_NULL_INTEGER)) then
            do p = 1, lgNLeaves
               glLocal(p) = lgLocal(p)
               glRemote(p)%rank = lgRemote(p)%rank
               glRemote(p)%index = lgRemote(p)%index
            end do
         else
            do p = 1, lgNLeaves
               glLocal(p) = p - 1
               glRemote(p)%rank = lgRemote(p)%rank
               glRemote(p)%index = lgRemote(p)%index
            end do
         end if
         if (loc(tempLocal) /= loc(PETSC_NULL_INTEGER)) then
            do p = 1, tempNLeaves
               glLocal(p + lgNLeaves) = tempLocal(p)
               glRemote(p + lgNLeaves)%rank = tempRemote(p)%rank
               glRemote(p + lgNLeaves)%index = tempRemote(p)%index
            end do
         else
            do p = 1, tempNLeaves
               glLocal(p + lgNLeaves) = p + lgNLeaves - 1
               glRemote(p + lgNLeaves)%rank = tempRemote(p)%rank
               glRemote(p + lgNLeaves)%index = tempRemote(p)%index
            end do
         end if
         PetscCall(PetscSFDestroy(sf, ierr))
         PetscCall(PetscSFCreate(MEF90Ctx%Comm, sf, ierr))
         PetscCall(PetscSFSetFromOptions(sf, ierr))
         PetscCall(PetscSFSetGraph(sf, glNRoots, glNLeaves, glLocal, PETSC_COPY_VALUES, glRemote, PETSC_COPY_VALUES, ierr))
         PetscCall(PetscSFSetUp(sf, ierr))
         PetscCall(PetscSFRestoreGraph(tempSF, tempNRoots, tempNLeaves, tempLocal, tempRemote, ierr))
         deallocate (glLocal)
         deallocate (glRemote)
         PetscCall(PetscSFRestoreGraph(tempSF, tempNRoots, tempNLeaves, tempLocal, tempRemote, ierr))
         PetscCall(PetscSFRestoreGraph(sf, lgNRoots, lgNLeaves, lgLocal, lgRemote, ierr))
         PetscCall(PetscSFDestroy(tempSF, ierr))
      end if
      PetscCall(PetscSectionDestroy(gSection, ierr))
      PetscCall(PetscSFDestroy(idSF, ierr))
      ! PetscCall(PetscObjectSetName(sf,"CGlobal-To-Local SF",ierr))
      ! PetscCall(PetscSFViewFromOptions(sf,PETSC_NULL_OBJECT,"-cglobaltolocal_sf_view",ierr))
      deallocate (remote)
   end subroutine CreateCGlobalToLocalSF_Private

#undef __FUNCT__
#define __FUNCT__ "CreateSideSF_Private"
!!!
!!!
!!!  CreateSideSF_Private:
!!!
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!

   subroutine CreateSideSF_Private(MEF90Ctx, dm, sf, invSF, ierr)
      type(tDM), intent(IN)                    :: dm
      type(MEF90Ctx_type), intent(IN)          :: MEF90Ctx
      type(tPetscSF), intent(OUT)              :: sf, invSF
      PetscErrorCode, intent(INOUT)            :: ierr

      type(tIS)                               :: ssIS, gssIS, faceIS, facesIS
      PetscInt, dimension(:), pointer           :: ssID, faceID
      type(tPetscSF)                          :: migrationSF, tempSF
      type(tVec)                              :: localVec
      PetscInt                                :: set, face, nroots, nleaves, i, totalleaves, numComponent, uNumComponent, j, numSS, key, numFaces
      type(sPetscSFNode), dimension(:), pointer :: iremote
      type(tIS), dimension(:), pointer          :: locfacesIS
      PetscInt, dimension(:), pointer           :: ilocal, permIndices, emptyInd, facesID, procSSID

      nleaves = 0_ki
      totalleaves = 0_ki
      PetscCall(DMGetLabelIdIS(dm, "Face Sets", gssIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90Ctx%comm, gssIS, ierr))
      PetscCall(ISGetSize(gssIS, numSS, ierr))
      PetscCall(DMGetLabelIdIS(dm, "Face Sets", ssIS, ierr))
      PetscCall(DMPlexGetMigrationSF(dm, migrationSF, ierr))
      PetscCall(ISGetIndices(ssIS, ssID, ierr))
      allocate (locfacesIS(numSS))
      allocate (procSSID(numSS))
      do set = 1, numSS
         procSSID(set) = -1
      end do
      do set = 1, size(ssID)
         PetscCall(ISLocate(gssIS, ssID(set), key, ierr))
         if (key >= 0) then
            procSSID(key + 1) = ssID(set)
         end if
      end do
      do set = 1, numSS
         PetscCall(DMGetStratumIS(dm, "Face Sets", procSSID(set), faceIS, ierr))
         if (.not. PetscObjectIsNull(faceIS)) then
            PetscCall(ISGetIndices(faceIS, faceID, ierr))
         else
            allocate (faceID(0))
         end if
         numFaces = size(faceID)
         allocate (ilocal(numFaces))
         allocate (iremote(numFaces))
         do face = 1, numFaces
            iremote(face)%rank = MEF90Ctx%Rank
            iremote(face)%index = face - 1
            ilocal(face) = faceID(face)
         end do
         PetscCall(PetscSFCreate(MEF90Ctx%Comm, tempSF, ierr))
         PetscCall(PetscSFSetFromOptions(tempSF, ierr))

         PetscCall(PetscSFSetGraph(tempSF, numFaces, numFaces, ilocal, PETSC_COPY_VALUES, iremote, PETSC_COPY_VALUES, ierr))
         PetscCall(PetscSFSetUp(tempSF, ierr))
         if (MEF90Ctx%NumProcs > 1) then
            PetscCall(PetscSFComposeInverse(migrationSF, tempSF, sf, ierr))
            PetscCall(PetscSFDestroy(tempSF, ierr))
         else
            PetscCall(PetscSFCreateInverseSF(tempSF, sf, ierr))
            PetscCall(PetscSFDestroy(tempSF, ierr))
         end if
         PetscCall(PetscSFSetUp(sf, ierr))
         deallocate (ilocal)
         deallocate (iremote)
         if (.not. PetscObjectIsNull(faceIS)) then
            PetscCall(ISRestoreIndices(faceIS, faceID, ierr))
         else
            deallocate (faceID)
         end if
         PetscCall(PetscSFGetGraph(sf, nroots, nleaves, emptyInd, iremote, ierr))
         allocate (ilocal(nleaves))
         do i = 1, nleaves
            ilocal(i) = iremote(i)%index
         end do
         PetscCall(ISCreateGeneral(MEF90Ctx%comm, nleaves, ilocal, PETSC_COPY_VALUES, locfacesIS(set), ierr))
         PetscCall(ISDestroy(faceIS, ierr))
         PetscCall(PetscSFDestroy(sf, ierr))
         deallocate (ilocal)
         PetscCall(PetscSFRestoreGraph(sf, nroots, nleaves, emptyInd, iremote, ierr))
      end do
      do set = 1, numSS
         PetscCall(MEF90ISAllGatherMerge(MEF90Ctx%comm, locfacesIS(set), ierr))
      end do
      PetscCall(ISConcatenate(MEF90Ctx%comm, numSS, locfacesIS, facesIS, ierr))
      PetscCall(ISGetIndices(facesIS, faceID, ierr))
      if (MEF90Ctx%rank == 0) then
         totalleaves = size(faceID)
      end if
      allocate (facesID(totalleaves))
      allocate (permIndices(totalleaves))
      do i = 1, totalleaves
         permIndices(i) = i - 1
         facesID(i) = faceID(i)
      end do
      PetscCall(ISRestoreIndices(facesIS, faceID, ierr))
      PetscCall(ISRestoreIndices(ssIS, ssID, ierr))
      PetscCall(PetscSortIntWithPermutation(totalleaves, facesID, permIndices, ierr))
      PetscCall(DMGetLocalVector(dm, localVec, ierr))
      PetscCall(VecGetBlockSize(localVec, numComponent, ierr))
      PetscCall(MPI_Allreduce(numComponent, uNumComponent, 1, MPIU_INTEGER, MPI_MAX, MEF90Ctx%comm, ierr))
      if ((numComponent == 1) .and. (uNumComponent > 1)) then
         numComponent = uNumComponent
      end if
      PetscCall(DMRestoreLocalVector(dm, localVec, ierr))
      allocate (iremote(numComponent * totalleaves))
      allocate (ilocal(numComponent * totalleaves))
      do i = 1, totalleaves
         do j = 1, numComponent
            ilocal((i - 1) * numComponent + j) = (i - 1) * numComponent + j - 1
            iremote((i - 1) * numComponent + j)%rank = MEF90Ctx%rank
            iremote((i - 1) * numComponent + j)%index = numComponent * permIndices(i) + j - 1
         end do
      end do
      PetscCall(PetscSFCreate(MEF90Ctx%Comm, invSF, ierr))
      PetscCall(PetscSFSetFromOptions(invSF, ierr))
      PetscCall(PetscSFSetGraph(invSF, numComponent * totalleaves, numComponent * totalleaves, ilocal, PETSC_COPY_VALUES, iremote, PETSC_COPY_VALUES, ierr))
      PetscCall(PetscSFSetUp(invSF, ierr))
      PetscCall(PetscSFCreateInverseSF(invSF, sf, ierr))
      PetscCall(PetscSFSetUp(sf, ierr))
      do set = 1, numSS
         PetscCall(ISDestroy(locfacesIS(set), ierr))
      end do
      deallocate (locfacesIS)
      deallocate (permIndices)
      deallocate (procSSID)
      deallocate (facesID)
      deallocate (ilocal)
      deallocate (iremote)
      PetscCall(ISDestroy(gssIS, ierr))
      PetscCall(ISDestroy(ssIS, ierr))
      PetscCall(ISDestroy(facesIS, ierr))
   end subroutine CreateSideSF_Private
end module m_MEF90_DMPlex
