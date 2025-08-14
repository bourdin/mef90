#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
module m_MEF90_DefMechAT
#include "petsc/finclude/petsc.h"
   use m_MEF90_DefMechAT_class
   use m_MEF90_DefMechAT1
   ! use m_MEF90_DefMechAT1exp
   ! use m_MEF90_DefMechAT2
   ! use m_MEF90_DefMechATKKL
   ! use m_MEF90_DefMechATLinSoft
   use m_MEF90_DefMechCtx
   use petscsys
   implicit none(type)

   private
   public :: MEF90DefMechGetATModel
   public :: MEF90DefMechAT_Type
   public :: MEF90DefMechAT1_Type
   ! public :: MEF90DefMechAT1exp_Type
   ! public :: MEF90DefMechAT2_Type
   ! public :: MEF90DefMechATKKL_Type
   ! public :: MEF90DefMechATLinSoft_Type

   public :: MEF90DefMech_damageTypeAT1

   enum, bind(c)
      enumerator :: MEF90DefMech_damageTypeAT1 = 0, &
         MEF90DefMech_damageTypeAT1exp, &
         MEF90DefMech_damageTypeAT2, &
         MEF90DefMech_damageTypeLinSoft, &
         MEF90DefMech_damageTypeKKL, &
         MEf90DefMech_damageTypeAT1Elastic, &
         MEf90DefMech_damageTypeAT1expElastic, &
         MEf90DefMech_damageTypeAT2Elastic, &
         MEF90DefMech_damageTypeLinSoftElastic, &
         MEF90DefMech_damageTypeKKLElastic
   end enum

   character(len=MEF90MXSTRLEN), dimension(13), protected   :: MEF90DefMech_damageTypeList = [ &
      'AT1                     ', &
      'AT1exp                  ', &
      'AT2                     ', &
      'LinSoft                 ', &
      'KKL                     ', &
      'AT1Elastic              ', &
      'AT1expElastic           ', &
      'AT2Elastic              ', &
      'LinSoftElastic          ', &
      'KKLElastic              ', &
      'MEF90DefMech_damageType ', &
      '_MEF90DefMech_damageType', &
      '                        '  &
   ]

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechGetATModel"
!!!
!!!
!!!  MEF90DefMechGetATModel: Return the AT model object from the cell set options
!!!
!!!  (c) 2025 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine MEF90DefMechGetATModel(comm, prefix, ATModel, ierr)
      MPI_Comm, intent(in) :: comm
      character(len = MEF90MXSTRLEN), intent(in) :: prefix
      class(MEF90DefMechAT_Type), allocatable, intent(out) :: ATModel
      PetscErrorCode, intent(inout) :: ierr

      PetscEnum :: damageType
      PetscCall(PetscOptionsBegin(comm, trim(prefix), "Options for MEF90DefMechAT_type", "mef90DefMech", ierr))
         PetscCall(PetscOptionsEnum("-damage_type", "AT model type", "mef90DefMech", MEF90DefMech_damageTypeList, MEF90DefMech_damageTypeAT1, damageType, PETSC_NULL_BOOL, ierr))
      PetscCall(PetscOptionsEnd(ierr))


      select case (damageType)
         case (MEF90DefMech_damageTypeAT1)
            ATModel = MEF90DefMechAT1_Type(comm, prefix, isElastic = PETSC_FALSE)
         case (MEF90DefMech_damageTypeAT1Elastic)
            ATModel = MEF90DefMechAT1_Type(comm, prefix, isElastic = PETSC_TRUE)
         ! case (MEF90DefMech_damageTypeAT1exp, MEF90DefMech_damageTypeAT1expElastic)
         !    ATModel = MEF90DefMechAT1exp_Type(cellSetOptions%DamageAT1expb)
         ! case (MEF90DefMech_damageTypeAT2, MEF90DefMech_damageTypeAT2Elastic)
         !    ATModel = MEF90DefMechAT2_Type()
         ! case (MEF90DefMech_damageTypeKKL, MEF90DefMech_damageTypeKKLElastic)
         !    ATModel = MEF90DefMechATKKL_Type()
         ! case (MEF90DefMech_damageTypeLinSoft, MEF90DefMech_damageTypeLinSoftElastic)
         !    ATModel = MEF90DefMechATLinSoft_Type(cellSetOptions%DamageATLinSoftk)
      case default
         print *, __FUNCT__, ': Unimplemented damage Type', damageType
         stop
      end select
   end subroutine MEF90DefMechGetATModel
end module m_MEF90_DefMechAT
