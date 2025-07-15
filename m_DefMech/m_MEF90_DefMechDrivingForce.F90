#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
module m_MEF90_DefMechDrivingForce
#include "petsc/finclude/petsc.h"
   use m_MEF90_DefMechDrivingForceDruckerPrager
   implicite none(type, external)

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechGetDrivingForce"
!!!
!!!
!!!  MEF90DefMechGetATModel: Return the AT model object from the cell set options
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   subroutine MEF90DefMechGetDrivingForce(cellSetOptions, DrivingForce)
      type(MEF90DefMechCellSetOptions_Type), pointer, intent(IN)      :: cellSetOptions
      class(MEF90_DefMechDrivingForce_Type), allocatable, intent(OUT) :: DrivingForce

      select case (cellSetOptions % drivingForceType)
      case (MEF90DefMech_drivingForceTypeDruckerPrager)
         DrivingForce = MEF90_DefMechDrivingForceDruckerPrager_Type()
      case (MEF90DefMech_drivingForceTypeDruckerPrager2)
         DrivingForce = MEF90_DefMechDrivingForceDruckerPrager2_Type()
      case default
         print *, __FUNCT__, ': Unimplemented damage Type, only DP and DP2 implemented', cellSetOptions % drivingForceType
         stop
      end select
   end subroutine MEF90DefMechGetDrivingForce
end module m_MEF90_DefMechDrivingForce
