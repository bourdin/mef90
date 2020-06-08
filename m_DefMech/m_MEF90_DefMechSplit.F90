#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
module MEF90_APPEND(m_MEF90_DefMechSplit,MEF90_DIM)D
   use MEF90_APPEND(m_MEF90_DefMechNone,MEF90_DIM)D
   use MEF90_APPEND(m_MEF90_DefMechMasonry,MEF90_DIM)D
   Use m_MEF90_DefMechCtx

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechGetSplit"
!!!
!!!  
!!!  MEF90DefMechGetSplit: Return the split object from the cell set options
!!!  
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechGetSplit(cellSetOptions,Split)
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Class(MEF90_DEFMECHSPLIT),Allocatable,intent(OUT)  :: Split

      Select Case(cellSetOptions%unilateralContactType)
      Case(MEF90DefMech_unilateralContactTypeNone)
         Split = MEF90_DEFMECHNONESPLIT()
      Case(MEF90DefMech_unilateralContactTypeMasonry)
         Split = MEF90_DEFMECHMASONRYSPLIT()
      Case default
         Print*,__FUNCT__,': Unimplemented split Type, only NONE and Masonry implemented',cellSetOptions%unilateralContactType
         STOP  
      End Select
   End Subroutine MEF90DefMechGetSplit
End module MEF90_APPEND(m_MEF90_DefMechSplit,MEF90_DIM)D
