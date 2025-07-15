#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
module m_MEF90_DefMechAT
#include "petsc/finclude/petsc.h"
   use m_MEF90_DefMechAT_class
   use m_MEF90_DefMechAT1
   use m_MEF90_DefMechAT1exp
   use m_MEF90_DefMechAT2
   use m_MEF90_DefMechATKKL
   use m_MEF90_DefMechATLinSoft
   use m_MEF90_DefMechCtx

   implicit none(type, external)
   ! private
   public :: MEF90DefMechGetATModel
   public :: MEF90DefMechAT_Type
   public :: MEF90DefMechAT1_Type
   public :: MEF90DefMechAT1exp_Type
   public :: MEF90DefMechAT2_Type
   public :: MEF90DefMechATKKL_Type
   public :: MEF90DefMechATLinSoft_Type

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechGetATModel"
!!!
!!!
!!!  MEF90DefMechGetATModel: Return the AT model object from the cell set options
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   subroutine MEF90DefMechGetATModel(cellSetOptions, ATModel, isElastic, ierr)
      type(MEF90DefMechCellSetOptions_Type), pointer      :: cellSetOptions
      class(MEF90DefMechAT_Type), allocatable, intent(OUT) :: ATModel
      PetscBool, intent(OUT)                              :: isElastic
      PetscErrorCode, intent(INOUT)                       :: ierr

      isElastic = .false.
      select case (cellSetOptions % damageType)
      case (MEF90DefMech_damageTypeAT1, MEF90DefMech_damageTypeAT1Elastic)
         ATModel = MEF90DefMechAT1_Type()
      case (MEF90DefMech_damageTypeAT1exp, MEF90DefMech_damageTypeAT1expElastic)
         ATModel = MEF90DefMechAT1exp_Type(cellSetOptions % DamageAT1expb)
      case (MEF90DefMech_damageTypeAT2, MEF90DefMech_damageTypeAT2Elastic)
         ATModel = MEF90DefMechAT2_Type()
      case (MEF90DefMech_damageTypeKKL, MEF90DefMech_damageTypeKKLElastic)
         ATModel = MEF90DefMechATKKL_Type()
      case (MEF90DefMech_damageTypeLinSoft, MEF90DefMech_damageTypeLinSoftElastic)
         ATModel = MEF90DefMechATLinSoft_Type(cellSetOptions % DamageATLinSoftk)
      case default
         print *, __FUNCT__, ': Unimplemented damage Type', cellSetOptions % damageType
         stop
      end select
      !!!  Check if the block is elastic
      select case (cellSetOptions % damageType)
      case (MEF90DefMech_damageTypeAT1Elastic, MEF90DefMech_damageTypeAT1expElastic, MEF90DefMech_damageTypeAT2Elastic, MEF90DefMech_damageTypeKKLElastic, MEF90DefMech_damageTypeLinSoftElastic)
         isElastic = .true.
      end select
   end subroutine MEF90DefMechGetATModel

end module m_MEF90_DefMechAT
