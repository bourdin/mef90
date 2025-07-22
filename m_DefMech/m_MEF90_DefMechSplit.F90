#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
module MEF90_APPEND(m_MEF90_DefMechSplit,MEF90_DIM)D
#include "petsc/finclude/petsc.h"
use MEF90_APPEND(m_MEF90_DefMechSplit_class,MEF90_DIM)D
use MEF90_APPEND(m_MEF90_DefMechSplitNone,MEF90_DIM)D
!use MEF90_APPEND(m_MEF90_DefMechSplitMasonry,MEF90_DIM)D
use MEF90_APPEND(m_MEF90_DefMechSplitHD,MEF90_DIM)D
use MEF90_APPEND(m_MEF90_DefMechSplitDeviatoric,MEF90_DIM)D
use MEF90_APPEND(m_MEF90_DefMechSplitHydrostatic,MEF90_DIM)D
use m_MEF90_DefMechCtx

implicit none(type, external)
! private
public :: MEF90DefMechGetSplit
public :: MEF90_DEFMECHSPLIT
public :: MEF90_DEFMECHSPLITDEVIATORIC
public :: MEF90_DEFMECHSPLITHD
public :: MEF90_DEFMECHSPLITHYDROSTATIC
! public :: MEF90_DEFMECHSPLITMASONRY
public :: MEF90_DEFMECHSPLITNONE

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechGetSplit"
!!!
!!!
!!!  MEF90DefMechGetSplit: Return the split object from the cell set options
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
subroutine MEF90DefMechGetSplit(cellSetOptions, Split, ierr)
   type(MEF90DefMechCellSetOptions_Type), pointer      :: cellSetOptions
   class(MEF90_DEFMECHSPLIT), allocatable, intent(OUT)  :: Split
   PetscErrorCode, intent(INOUT)                       :: ierr

   select case (cellSetOptions%unilateralContactType)
   case (MEF90DefMech_unilateralContactTypeNone)
      Split = MEF90_DEFMECHSPLITNONE()
      !Case(MEF90DefMech_unilateralContactTypeMasonry)
      !   Split = MEF90_DEFMECHSPLITMASONRY()
   case (MEF90DefMech_unilateralContactTypeHydrostaticDeviatoric)
      Split = MEF90_DEFMECHSPLITHD(cellSetOptions%unilateralContactHydrostaticDeviatoricGamma)
   case (MEF90DefMech_unilateralContactTypeDeviatoric)
      Split = MEF90_DEFMECHSPLITDEVIATORIC()
   case (MEF90DefMech_unilateralContactTypeHydrostatic)
      Split = MEF90_DEFMECHSPLITHYDROSTATIC(cellSetOptions%unilateralContactHydrostaticDeviatoricGamma)
   case default
      print *, __FUNCT__, ': Unimplemented split Type', cellSetOptions%unilateralContactType
      stop
   end select
end subroutine MEF90DefMechGetSplit
end module MEF90_APPEND(m_MEF90_DefMechSplit,MEF90_DIM)D
