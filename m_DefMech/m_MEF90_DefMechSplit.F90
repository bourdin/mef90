#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
module MEF90_APPEND(m_MEF90_DefMechSplit,MEF90_DIM)D
   use MEF90_APPEND(m_MEF90_DefMechSplitNone,MEF90_DIM)D
   !use MEF90_APPEND(m_MEF90_DefMechSplitMasonry,MEF90_DIM)D
   use MEF90_APPEND(m_MEF90_DefMechSplitHD,MEF90_DIM)D
   use MEF90_APPEND(m_MEF90_DefMechSplitDeviatoric,MEF90_DIM)D
   use MEF90_APPEND(m_MEF90_DefMechSplitHydrostatic,MEF90_DIM)D
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
         Split = MEF90_DEFMECHSPLITNONE()
      !Case(MEF90DefMech_unilateralContactTypeMasonry)
      !   Split = MEF90_DEFMECHSPLITMASONRY()
      Case(MEF90DefMech_unilateralContactTypeHydrostaticDeviatoric)
         Split = MEF90_DEFMECHSPLITHD(cellSetOptions%unilateralContactHydrostaticDeviatoricGamma)
      Case(MEF90DefMech_unilateralContactTypeDeviatoric)
         Split = MEF90_DEFMECHSPLITDEVIATORIC()
      Case(MEF90DefMech_unilateralContactTypeHydrostatic)
         Split = MEF90_DEFMECHSPLITHYDROSTATIC(cellSetOptions%unilateralContactHydrostaticDeviatoricGamma)
      Case default
         Print*,__FUNCT__,': Unimplemented split Type',cellSetOptions%unilateralContactType
         STOP  
      End Select
   End Subroutine MEF90DefMechGetSplit
End module MEF90_APPEND(m_MEF90_DefMechSplit,MEF90_DIM)D
