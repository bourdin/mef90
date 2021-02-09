#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
module m_MEF90_DefMechAT
#include "finclude/petscdef.h"
   use m_MEF90_DefMechAT1
   use m_MEF90_DefMechAT1exp
   use m_MEF90_DefMechAT2
   use m_MEF90_DefMechATKKL
   use m_MEF90_DefMechATLinSoft
   use m_MEF90_DefMechCtx

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechGetATModel"
!!!
!!!  
!!!  MEF90DefMechGetATModel: Return the AT model object from the cell set options
!!!  
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechGetATModel(cellSetOptions,ATModel,isElastic)
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Class(MEF90_DefMechAT_Type),Allocatable,Intent(OUT):: ATModel
      PetscBool,Intent(OUT)                              :: isElastic

      isElastic = .FALSE.
      Select Case (cellSetOptions%damageType)
      Case (MEF90DefMech_damageTypeAT1,MEF90DefMech_damageTypeAT1Elastic)
         ATModel = MEF90_DefMechAT1_Type()
      Case (MEF90DefMech_damageTypeAT1exp,MEF90DefMech_damageTypeAT1expElastic)
         ATModel = MEF90_DefMechAT1exp_Type(cellSetOptions%DamageAT1expb)
      Case (MEF90DefMech_damageTypeAT2,MEF90DefMech_damageTypeAT2Elastic)
         ATModel = MEF90_DefMechAT2_Type()
      Case (MEF90DefMech_damageTypeKKL,MEF90DefMech_damageTypeKKLElastic)
         ATModel = MEF90_DefMechATKKL_Type()
      Case (MEF90DefMech_damageTypeLinSoft,MEF90DefMech_damageTypeLinSoftElastic)
         ATModel = MEF90_DefMechATLinSoft_Type(cellSetOptions%DamageATLinSoftk)
      Case default
         Print*,__FUNCT__,': Unimplemented damage Type',cellSetOptions%damageType
         STOP  
      End Select
      !!!  Check if the block is elastic 
      Select Case (cellSetOptions%damageType)
      Case (MEF90DefMech_damageTypeAT1Elastic,MEF90DefMech_damageTypeAT1expElastic,MEF90DefMech_damageTypeAT2Elastic,MEF90DefMech_damageTypeKKLElastic,MEF90DefMech_damageTypeLinSoftElastic)
         isElastic = .TRUE.
      End Select
   End Subroutine MEF90DefMechGetATModel

End module m_MEF90_DefMechAT
