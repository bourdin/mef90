#include "../MEF90/mef90.inc"
Module MEF90_APPEND(m_MEF90_GradDamageImplementation_,MEF90_DIM)D
#include "finclude/petscdef.h"
   Use m_MEF90_LinAlg
   Use m_MEF90_Parameters
   Use m_MEF90_Elements
   Use m_MEF90_Utils
   Use m_MEF90_Materials
   IMPLICIT NONE
   Private

   Public :: MEF90GradDamageCellAverage

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageCellAverage"
   !!!
   !!!  
   !!!  MEF90GradDamageCellAverage: Compute a cell--average section of the average damage in each cell
   !!!  
   !!!  (c) 2015 Erwan TANNE erwan.tanne@gmail.com
   !!!
   
   Subroutine MEF90GradDamageCellAverage(damageCellAvg,damageLoc,elemDamage,elemDamageType,ierr)
      PetscReal,Intent(OUT)                              :: damageCellAvg
      PetscReal,Dimension(:),pointer                     :: damageLoc
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      Type(MEF90Element_Type),Intent(IN)                 :: elemDamageType
      PetscErrorCode,Intent(OUT)                         :: ierr

      PetscReal                                          :: cellSize
      PetscInt                                           :: iDoF1,iGauss
      PetscLogDouble                                     :: flops
     
      cellSize = 0.0_Kr
      Do iGauss = 1,size(elemDamage%Gauss_C)
         damageCellAvg = 0.0_Kr
         cellSize = 0.0_Kr
         Do iDoF1 = 1,elemDamageType%numDof

            damageCellAvg = damageCellAvg + damageLoc(iDof1) * elemDamage%BF(iDof1,iGauss) * elemDamage%Gauss_C(iGauss)
            cellSize = cellSize +  elemDamage%BF(iDof1,iGauss) * elemDamage%Gauss_C(iGauss)
         End Do
      End Do
      damageCellAvg = damageCellAvg / cellSize

      flops = 5 * elemDamageType%numDof * size(elemDamage%Gauss_C) + 1 
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageCellAverage

End Module MEF90_APPEND(m_MEF90_GradDamageImplementation_,MEF90_DIM)D


