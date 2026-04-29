#include "../MEF90/mef90.inc"
module m_MEF90_DefMechAT
#include "petsc/finclude/petsc.h"
   use m_MEF90_DefMechAT_class
   use m_MEF90_DefMechAT1
   use m_MEF90_DefMechAT1exp
   use m_MEF90_DefMechAT2
   use m_MEF90_DefMechATKKL
   ! use m_MEF90_DefMechATLinSoft
   use m_MEF90_DefMechCtx
   use petscsys
   implicit none(type)

   private
   public :: MEF90DefMechGetATModel
   public :: MEF90DefMechAT_Type
   public :: MEF90DefMechAT1_Type
   ! public :: MEF90DefMechAT1exp_Type
   public :: MEF90DefMechAT2_Type
   ! public :: MEF90DefMechATKKL_Type
   ! public :: MEF90DefMechATLinSoft_Type

   public :: MEF90DefMech_damageTypeAT1

   enum, bind(c)
      enumerator :: MEF90DefMech_damageTypeAT1 = 0, &
         MEF90DefMech_damageTypeAT1exp, &
         MEF90DefMech_damageTypeAT2, &
         MEF90DefMech_damageTypeLinSoft, &
         MEF90DefMech_damageTypeKKL
   end enum

   character(len=MEF90MXSTRLEN), dimension(8), protected   :: MEF90DefMech_damageTypeList = [ &
      'AT1                     ', &
      'AT1exp                  ', &
      'AT2                     ', &
      'LinSoft                 ', &
      'KKL                     ', &
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
   subroutine MEF90DefMechGetATModel(comm, prefix, dim, ATModel, ierr)
      MPIU_Comm, intent(in) :: comm
      character(len = MEF90MXSTRLEN), intent(in) :: prefix
      PetscInt, intent(IN) :: dim
      class(MEF90DefMechAT_Type), allocatable, intent(out) :: ATModel
      PetscErrorCode, intent(inout) :: ierr

      PetscEnum :: damageType
      PetscCall(PetscOptionsBegin(comm, trim(prefix), "Options for MEF90DefMechAT_type", "mef90DefMech", ierr))
         PetscCall(PetscOptionsEnum("-damage_type", "AT model type", "mef90DefMech", MEF90DefMech_damageTypeList, MEF90DefMech_damageTypeAT1, damageType, PETSC_NULL_BOOL, ierr))
      PetscCall(PetscOptionsEnd(ierr))


      select case (damageType)
         case (MEF90DefMech_damageTypeAT1)
            ATModel = MEF90DefMechAT1_Type(comm, prefix)
            case (MEF90DefMech_damageTypeAT1exp)
               ATModel = MEF90DefMechAT1exp_Type(comm, prefix)
            case (MEF90DefMech_damageTypeAT2)
               ATModel = MEF90DefMechAT2_Type(comm, prefix)
            case (MEF90DefMech_damageTypeKKL)
               ATModel = MEF90DefMechATKKL_Type(comm, prefix)
            ! case (MEF90DefMech_damageTypeLinSoft)
            !    ATModel = MEF90DefMechATLinSoft_Type(comm, prefix)
         case default
            print *, __FUNCT__, ': Unimplemented damage Type', damageType
            stop
      end select
      select case(dim)
         case(2)
            ATModel%toughnessAnisotropyMatrix = MEF90MatS2DIdentity
         case(3)
            ATModel%toughnessAnisotropyMatrix = MEF90MatS3DIdentity
         case default
            write(*,"('In ', A, ' dim = ', I2, 'and should be 2 or 3')") __FUNCT__, dim
            STOP
      end select
   end subroutine MEF90DefMechGetATModel
end module m_MEF90_DefMechAT
