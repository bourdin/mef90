#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
module MEF90_APPEND(m_MEF90_DefMechSplit,MEF90_DIM)D
   use MEF90_APPEND(m_MEF90_DefMechNone,MEF90_DIM)D
   use MEF90_APPEND(m_MEF90_DefMechMasonry,MEF90_DIM)D
End module MEF90_APPEND(m_MEF90_DefMechSplit,MEF90_DIM)D
