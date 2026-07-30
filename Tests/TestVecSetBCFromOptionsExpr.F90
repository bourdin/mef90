module TestExpr
#include <petsc/finclude/petsc.h>
   use petsc
   use m_MEF90
   use m_MEF90_DMPlex
   use symengine

   implicit none(type, external)
contains

#undef __FUNCT__
#define __FUNCT__ "myMEF90VecSetValuesFromOptionsExpr"
!!!
!!!
!!!  myMEF90VecSetValuesFromOptionsExpr: Fill values of a Vec using expressions passed as petsc options
!!!
!!!  (c) 2025      Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine myMEF90VecSetValuesFromOptionsExpr(v, t, ierr)
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
      type(realdouble), dimension(3)            :: vals
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

      vals(1) = RealDouble(t)
      do setType = 1, size(MEF90SetType)
         PetscCall(DMGetLabelIdIS(dm, MEF90SetLabelName(setType), setIS, ierr))
         ! PetscCall(MEF90ISAllGatherMerge(comm,setIS,ierr))
         if (.not. PetscObjectIsNull(setIS)) then
            PetscCall(ISGetIndices(setIS, setID, ierr))
            do set = 1, size(setID)
               write (ValueKey, '("-",a2,I4.4,"_",a)') MEF90SetPrefix(setType), setID(set), trim(name)
               PetscCall(PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, trim(ValueKey), ExprStr, flg, ierr))
               if (flg) then
                  call MEF90StrTokenize(ExprStr, delim, ExprStrComp)
                  numOpt = size(ExprStrComp)
                  if (numOpt /= bs) then
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
   end subroutine myMEF90VecSetValuesFromOptionsExpr

#undef __FUNCT__
#define __FUNCT__ "myMEF90VecSetBCValuesFromOptionsExpr"
!!!
!!!
!!!  myMEF90VecSetBCValuesFromOptionsExpr: Fill boundary values of a Vec using expressions passed as petsc options
!!!
!!!  (c) 2025      Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine myMEF90VecSetBCValuesFromOptionsExpr(v, t, ierr)
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
      type(realdouble), dimension(3)            :: vals
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

      vals(1) = RealDouble(t)
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
                  PetscCall(PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, trim(BCValueKey), ExprStr, flg, ierr))
                  numOpt = MEF90StrCount(ExprStr, delim)
                  call MEF90StrTokenize(ExprStr, delim, ExprStrComp)
                  if (size(ExprStrComp) /= bs) then
                     write (IOBuffer, "(A,' was expecting ',I2,' expressions but got ',I2,' for key ',A)") __FUNCT__, bs, size(ExprStrComp), trim(BCValueKey)
                     SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, IOBuffer)
                  end if
                        !!! create parsers
                  do c = 1, bs
                     exprs(c) = parse(ExprStrComp(c))
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
   end subroutine myMEF90VecSetBCValuesFromOptionsExpr

end module TestExpr

program TestVecSetBCFromOptionsExpr
#include <petsc/finclude/petsc.h>
   use petsc
   use m_MEF90
   use m_MEF90_DMPlex
   use TestExpr

   implicit none(type, external)

   PetscErrorCode                                     :: ierr
   type(MEF90Ctx_Type), target                         :: MEF90Ctx
   type(MEF90CtxGlobalOptions_Type)                   :: MEF90GlobalOptions_default
   type(MEF90CtxGlobalOptions_Type)                    :: MEF90GlobalOptions
   type(tDM)                                          :: dm
   PetscBool                                          :: flg, interpolate = PETSC_TRUE
   character(len=MEF90MXSTRLEN)                       :: name
   PetscInt                                           :: sDim = 1

   type(tVec)                                         :: v
   PetscReal                                          :: scalingFactor

   MEF90GlobalOptions_default%verbose = 1
   MEF90GlobalOptions_default%dryrun = PETSC_FALSE
   MEF90GlobalOptions_default%timeInterpolation = MEF90TimeInterpolation_linear
   MEF90GlobalOptions_default%timeMin = 0.0_kr
   MEF90GlobalOptions_default%timeMax = 1.0_kr
   MEF90GlobalOptions_default%timeNumStep = 11
   MEF90GlobalOptions_default%timeSkip = 0
   MEF90GlobalOptions_default%timeNumCycle = 1
   MEF90GlobalOptions_default%elementFamily = MEF90ElementFamilyLagrange
   MEF90GlobalOptions_default%elementOrder = 1

   PetscCallA(PetscInitialize(ierr))

   call MEF90Initialize(PETSC_COMM_WORLD, ierr)
   call MEF90CtxCreate(PETSC_COMM_WORLD, MEF90Ctx, "", ierr)
   MEF90GlobalOptions = MEF90GlobalOptions_default
   call MEF90Ctx%setFromOptions(ierr); CHKERRQ(ierr)
   PetscCallA(MEF90CtxGlobalOptionsSetFromOptions(MEF90Ctx%comm, trim(MEF90Ctx%prefix), MEF90GlobalOptions, ierr))
   PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm, MEF90Ctx%geometryfile, PETSC_NULL_CHARACTER, interpolate, dm, ierr))
   PetscCallA(DMPlexDistributeSetDefault(dm, PETSC_FALSE, ierr))
   PetscCallA(DMSetFromOptions(dm, ierr))
   PetscCallA(DMViewFromOptions(dm, PETSC_NULL_OBJECT, "-mef90dm_view", ierr))

   distribute: block
      type(tDM), target                    :: dmDist
      PetscInt                            :: ovlp = 0
      if (MEF90Ctx%NumProcs > 1) then
         PetscCallA(DMPlexDistribute(dm, ovlp, PETSC_NULL_SF, dmDist, ierr))
         PetscCallA(DMDestroy(dm, ierr))
         dm = dmDist
      end if
   end block distribute

   name = "Temperature"
   PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-sdim', sdim, flg, ierr))
   PetscCallA(MEF90CreateLocalVector(dm, MEF90GlobalOptions%elementFamily, MEF90GlobalOptions%elementOrder, sdim, name, V, ierr))

   PetscCallA(VecSet(v, 5.678_kr, ierr))
   scalingFactor = 2.0_kr
   PetscCallA(MEF90VecSetValuesFromOptionsExpr(V, scalingFactor, ierr))
   PetscCallA(VecViewFromOptions(V, PETSC_NULL_OBJECT, "-temperature_vec_view", ierr))

   PetscCallA(VecSet(v, -1.234_kr, ierr))
   scalingFactor = 1.0_kr
   PetscCallA(MEF90VecSetBCValuesFromOptionsExpr(V, scalingFactor, ierr))
   PetscCallA(VecViewFromOptions(V, PETSC_NULL_OBJECT, "-temperature_vec_view", ierr))

   ViewSec: block
      type(tPetscSection)     :: sectionV
      type(tDM)               :: dmV

      PetscCallA(VecGetDM(V, dmV, ierr))
      PetscCallA(DMGetLocalSection(dmV, sectionV, ierr))
      PetscCallA(PetscSectionViewFromOptions(sectionV, PETSC_NULL_OBJECT, "-temperature_section_view", ierr))
   end block ViewSec

   PetscCallA(VecDestroy(v, ierr))
   PetscCallA(DMDestroy(dm, ierr))
   call MEF90CtxDestroy(MEF90Ctx, ierr)
   call MEF90Finalize(ierr)
   call PetscFinalize(ierr)
end program TestVecSetBCFromOptionsExpr

!!! TEST
!!! ./TestVecSetBCFromOptions -geometry ../TestMeshes/SquareFaceSet.gen -temperature_section_view    \
!!!     -temperature_vec_view -sdim 2 -fs0021_temperatureBC yes,no -fs0021_boundaryTemperature 10,20 \
!!!     -vs0010_temperatureBC yes,yes -vs0010_boundaryTemperature 100,101 -cs0001_temperature 1000   \
!!!     -vs0010_temperature -2 -fs0020_temperature 20 -fs0021_temperature 21 -vs0010_temperature 100
