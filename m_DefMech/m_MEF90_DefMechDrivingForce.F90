#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
module m_MEF90_DefMechDrivingForce
#include "finclude/petscdef.h"
   use m_MEF90_DefMechDrivingForceDruckerPrager
   use m_MEF90_DefMechCtx

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechGetDrivingForce"
!!!
!!!  
!!!  MEF90DefMechGetATModel: Return the AT model object from the cell set options
!!!  
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechGetDrivingForce(cellSetOptions,DrivingForce)
      Type(MEF90DefMechCellSetOptions_Type),Pointer                 :: cellSetOptions
      Class(MEF90_DefMechDrivingForce_Type),Allocatable,Intent(OUT) :: DrivingForce

      Select Case (cellSetOptions%drivingForceType)
      Case (MEF90DefMech_drivingForceTypeDruckerPrager)
         DrivingForce = MEF90_DefMechDrivingForceDruckerPrager_Type()
      Case (MEF90DefMech_drivingForceTypeDruckerPrager2)
         DrivingForce = MEF90_DefMechDrivingForceDruckerPrager2_Type()
      Case default
         Print*,__FUNCT__,': Unimplemented damage Type, only DP and DP2 implemented',cellSetOptions%drivingForceType
         STOP  
      End Select
   End Subroutine MEF90DefMechGetDrivingForce
End module m_MEF90_DefMechDrivingForce
