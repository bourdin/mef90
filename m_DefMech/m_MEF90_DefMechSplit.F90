#include "../MEF90/mef90.inc"
module m_MEF90_DefMechSplit
#include "petsc/finclude/petsc.h"
use m_MEF90_DefMechSplit_class
use m_MEF90_DefMechSplitNone
! use m_MEF90_DefMechSplitHD
! use m_MEF90_DefMechSplitDeviatoric
! use m_MEF90_DefMechSplitHydrostatic
use m_MEF90_DefMechCtx
use petscsys

implicit none(type)
! private
public :: MEF90DefMechGetSplit
! public :: MEF90DefMechSplitDeviatoric
! public :: MEF90DefMechSplitHD
! public :: MEF90DefMechSplitHydrostatic
public :: MEF90DefMechSplitNone

enum, bind(c)
   enumerator :: MEF90DefMechSplitEnumNone = 0, &
      MEF90DefMechSplitEnumHydrostaticDeviatoric, &
      MEF90DefMechSplitEnumHydrostatic, &
      MEF90DefMechSplitEnumDeviatoric
end enum
character(len=MEF90MXSTRLEN), dimension(4), protected :: MEF90DefMechSplitEnumList = [ &
         'None                      ', &
         ! 'HydrostaticDeviatoric     ', &
         ! 'Hydrostatic               ', &
         ! 'Deviatoric                ', &
         'MEF90DefMechSplitEnumList ', &
         '_MEF90DefMechSplitEnumList', &
         '                          '  &
      ]


contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechGetSplit"
!!!
!!!
!!!  MEF90DefMechGetSplit: Return the split object from the cell set options
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!      2026 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine MEF90DefMechGetSplit(comm, prefix, split, ierr)
      MPIU_Comm, intent(in)                              :: comm
      character(len = MEF90MXSTRLEN), intent(in)         :: prefix
      class(MEF90DefMechSplit), allocatable, intent(out) :: split
      PetscErrorCode, intent(inout)                      :: ierr

      PetscEnum :: splitType

      splitType = MEF90DefMechSplitEnumNone
      PetscCall(PetscOptionsBegin(comm, trim(prefix), "Options for MEF90DefMechSplit", "MEF90", ierr))
         PetscCall(PetscOptionsEnum("-unilateralContact_type", "unilateral contact type", "MEF90", MEF90DefMechSplitEnumList, MEF90DefMechSplitEnumNone, splitType, PETSC_NULL_BOOL, ierr))
      PetscCall(PetscOptionsEnd(ierr))

      select case (splitType)
      case (MEF90DefMechSplitEnumNone)
         Split = MEF90DefMechSplitNone(comm = comm, prefix = prefix)
      ! case (MEF90DefMech_unilateralContactTypeHydrostaticDeviatoric)
      !    Split = MEF90DefMechSplitHD(cellSetOptions%unilateralContactHydrostaticDeviatoricGamma)
      ! case (MEF90DefMech_unilateralContactTypeDeviatoric)
      !    Split = MEF90DefMechSplitDeviatoric()
      ! case (MEF90DefMech_unilateralContactTypeHydrostatic)
      !    Split = MEF90DefMechSplitHydrostatic(cellSetOptions%unilateralContactHydrostaticDeviatoricGamma)
      case default
         print *, __FUNCT__, ': Unimplemented split Type', MEF90DefMechSplitEnumNone
         stop
      end select
   end subroutine MEF90DefMechGetSplit
end module m_MEF90_DefMechSplit
