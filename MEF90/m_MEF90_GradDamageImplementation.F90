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

   Public :: MEF90GradDamageElasticEnergySet   
   Public :: MEF90GradDamageSurfaceEnergySetAT1
   Public :: MEF90GradDamageSurfaceEnergySetAT2
   Public :: MEF90GradDamageSurfaceEnergySetLinSoft
   Public :: MEF90GradDamageCellAverage

!!! The Euler-Lagrange equation for the alpha problem for AT1 is
!!!
!!! \int_\Omega 2W(u,T,p) \alpha \beta + \frac{G_c\varepsilon}{2c_V}\nabla \alpha \cdot \nabla \beta \, dx = 
!!! \int_\Omega \left(2W(u,T,p)-\frac{G_c}{4c_V\varepsilon}\right) \beta \, dx = 0 for all \beta,
!!! where u denotes the displacement field, T the temperature field, p the plastic strain, and W is the elastic energy density.
!!!
!!! For AT2, the E-L equation is
!!!
!!! \int_\Omega 2\left(W(u,T,p) + \frac{G_c}{4c_V\varepsilon}\right) \alpha \beta + \frac{G_c\varepsilon}{2c_V}\nabla \alpha \cdot \nabla \beta \, dx = 
!!! \int_\Omega 2W(u,T,p) \beta \, dx = 0 for all \beta,
!!!

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageElasticEnergySet"
   !!!
   !!!  
   !!!  MEF90GradDamageElasticEnergySet:  Contribution of a cell set to elastic energy. 
   !!!                        It is assumed that the temperature is interpolated on the FE space while the plastic strain 
   !!!                        is cell-based
   !!!  
   !!!  (c) 2012-2019 Blaise Bourdin bourdin@lsu.edu
   !!!
   Subroutine MEF90GradDamageElasticEnergySet(energy,x,alpha,plasticStrain,temperature,mesh,meshScal,cellIS,eta,HookesLaw,CoefficientLinSoft,ThermalExpansion,elemDisplacement,elemDisplacementType,elemScal,elemScalType,ierr)
      PetscReal,Intent(OUT)                              :: energy
      Type(SectionReal),Intent(IN)                       :: x,alpha,plasticStrain,temperature
      Type(DM),Intent(IN)                                :: mesh,meshScal
      Type(IS),Intent(IN)                                :: cellIS
      PetscReal                                          :: eta
      Type(MEF90_HOOKESLAW),Intent(IN)                   :: HookesLaw
      PetscReal,Intent(IN)                               :: CoefficientLinSoft
      Type(MEF90_MATS),Intent(IN)                        :: ThermalExpansion
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elemScal
      Type(MEF90Element_Type),Intent(IN)                 :: elemDisplacementType,elemScalType
      PetscErrorCode,Intent(OUT)                         :: ierr

      PetscReal,Dimension(:),Pointer                     :: xloc,alphaLoc,temperatureLoc,plasticStrainLoc
      PetscReal                                          :: alphaElem,SalphaElem,temperatureElem,Stress_ZZ_planeStrain
      Type(MEF90_MATS)                                   :: StrainElem,StressElem,plasticStrainElem
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscInt                                           :: iDoF1,iGauss
      PetscLogDouble                                     :: flops
     
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(xloc(elemDisplacementType%numDof))
         Allocate(temperatureloc(elemScalType%numDof))
         Allocate(alphaLoc(elemScalType%numDof))
         Do cell = 1,size(cellID)   
            Call SectionRealRestrictClosure(x,mesh,cellID(cell),elemDisplacementType%numDof,xloc,ierr);CHKERRQ(ierr)
            If (temperature%v /= 0) Then
               Call SectionRealRestrictClosure(temperature,meshScal,cellID(cell),elemScalType%numDof,temperatureLoc,ierr);CHKERRQ(ierr)
            End If
            If (plasticStrain%v /= 0) Then
               Call SectionRealRestrict(plasticStrain,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
            End If
            If (alpha%v /= 0) Then
               Call SectionRealRestrictClosure(alpha,meshScal,cellID(cell),elemScalType%numDof,alphaLoc,ierr);CHKERRQ(ierr)
            End If
            Do iGauss = 1,size(elemDisplacement(cell)%Gauss_C)
               alphaElem = 0.0_Kr
               If (alpha%v /= 0) Then
                  Do iDoF1 = 1, elemScalType%numDof
                     alphaElem = alphaElem + alphaLoc(iDoF1) * elemScal(cell)%BF(iDoF1,iGauss)
                  End Do ! iDoF1
               End If

               If (CoefficientLinSoft == 0) Then
                  SalphaElem = eta + (1.0_Kr - alphaElem)**2
               Else 
                  SalphaElem = eta + (1.0_Kr - alphaElem)**2 / ( 1.0_Kr + ( CoefficientLinSoft - 1.0_Kr )*(1.0_Kr - (1.0_Kr - alphaElem)**2 ) )
               End IF


               strainElem = 0.0_Kr
               Do iDoF1 = 1,elemDisplacementType%numDof
                  strainElem = strainElem + xloc(iDof1) * elemDisplacement(cell)%GradS_BF(iDof1,iGauss)
               End Do

               temperatureElem   = 0.0_Kr
               If (temperature%v /= 0) Then
                  Do iDoF1 = 1,elemScalType%numDof
                     temperatureElem = temperatureElem + temperatureLoc(iDof1) * elemScal(cell)%BF(iDof1,iGauss)
                  End Do
                  strainElem = strainElem - (temperatureElem * thermalExpansion)
               End If

               plasticStrainElem = 0.0_Kr
               If (plasticStrain%v /= 0) Then
                  plasticStrainElem = plasticStrainLoc
               End If
               stressElem = HookesLaw * (strainElem - plasticStrainElem)
               energy = energy + SalphaElem * ( (strainElem - plasticStrainElem) .dotP. stressElem) * elemDisplacement(cell)%Gauss_C(iGauss) * 0.5_Kr

#if MEF90_DIM == 2
               if ( HookesLaw%isPlaneStress .eqv. .FALSE. ) then
                  Stress_ZZ_planeStrain = ( HookesLaw%YoungsModulus - 2.0_Kr*HookesLaw%PoissonRatio*HookesLaw%mu )*trace(plasticStrainElem) + HookesLaw%lambda*trace(strainElem)
                  energy = energy + 0.5_Kr * elemDisplacement(cell)%Gauss_C(iGauss)* SalphaElem * ( Stress_ZZ_planeStrain + HookesLaw%lambda*trace(strainElem - plasticStrainElem) ) * trace(plasticStrainElem)
               endif
#endif

            End Do ! Gauss
            If (plasticStrain%v /= 0) Then
               Call SectionRealRestore(plasticStrain,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
            End If
         End Do ! cell
         flops = 7 * size(elemDisplacement(1)%Gauss_C) * size(cellID) 
         If (alpha%v /= 0) Then
            flops = flops + 2 * elemScalType%numDof * size(elemScal(1)%Gauss_C) * size(cellID) 
         End If
         If (temperature%v /= 0) Then
            flops = flops + 2 * elemScalType%numDof * size(elemScal(1)%Gauss_C) * size(cellID) 
         End If
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(xloc)
         DeAllocate(temperatureloc)
         DeAllocate(alphaloc)
      End If   
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageElasticEnergySet


#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageSurfaceEnergySetAT1"
!!!
!!!  
!!!  MEF90GradDamageSurfaceEnergySetAT1:
!!!  
!!!  (c) 2014-2019 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90GradDamageSurfaceEnergySetAT1(energy,alpha,meshScal,cellIS,internalLength,fractureToughness,elemScal,elemScalType,ierr)
      PetscReal,Intent(OUT)                              :: energy
      Type(SectionReal),Intent(IN)                       :: alpha
      Type(DM),Intent(IN)                                :: meshScal
      Type(IS),Intent(IN)                                :: cellIS
      PetscReal                                          :: internalLength,fractureToughness
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elemScal
      Type(MEF90Element_Type),Intent(IN)                 :: elemScalType
      PetscErrorCode,Intent(OUT)                         :: ierr

      PetscReal,Dimension(:),Pointer                     :: alphaLoc
      PetscReal                                          :: alphaElem
      Type(MEF90_VECT)                                   :: gradAlphaElem
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscInt                                           :: iDoF1,iGauss
      PetscReal                                          :: C1, C2
      PetscLogDouble                                     :: flops

      
      C1 = fractureToughness / internalLength * .375_Kr
      C2 = fractureToughness * internalLength * .375_Kr
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(alphaloc(elemScalType%numDof))
         Do cell = 1,size(cellID)   
            Call SectionRealRestrictClosure(alpha,meshScal,cellID(cell),elemScalType%numDof,alphaLoc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elemScal(cell)%Gauss_C)
               alphaElem     = 0.0_Kr
               gradAlphaElem = 0.0_Kr
               Do iDoF1 = 1,elemScalType%numDof
                  alphaElem     = alphaElem     + alphaLoc(iDof1) * elemScal(cell)%BF(iDoF1,iGauss)
                  gradAlphaElem = gradAlphaElem + alphaLoc(iDof1) * elemScal(cell)%Grad_BF(iDoF1,iGauss)
               End Do
               energy = energy + elemScal(cell)%Gauss_C(iGauss) * (alphaElem * C1 + (gradAlphaElem .dotP. gradAlphaElem) * C2)
            End Do ! Gauss
         End Do ! cell
         flops = (2 * elemScalType%numDof + 5 )* size(elemScal(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(alphaloc)
      End If   
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageSurfaceEnergySetAT1

#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageSurfaceEnergySetAT2"
!!!
!!!  
!!!  MEF90GradDamageSurfaceEnergySetAT2:
!!!  
!!!  (c) 2014-2019 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90GradDamageSurfaceEnergySetAT2(energy,alpha,meshScal,cellIS,internalLength,fractureToughness,elemScal,elemScalType,ierr)
      PetscReal,Intent(OUT)                              :: energy
      Type(SectionReal),Intent(IN)                       :: alpha
      Type(DM),Intent(IN)                                :: meshScal
      Type(IS),Intent(IN)                                :: cellIS
      PetscReal                                          :: internalLength,fractureToughness
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elemScal
      Type(MEF90Element_Type),Intent(IN)                 :: elemScalType
      PetscErrorCode,Intent(OUT)                         :: ierr

      PetscReal,Dimension(:),Pointer                     :: alphaLoc
      PetscReal                                          :: alphaElem
      Type(MEF90_VECT)                                   :: gradAlphaElem
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscInt                                           :: iDoF1,iGauss
      PetscReal                                          :: C1, C2
      PetscLogDouble                                     :: flops

      
      C1 = fractureToughness / internalLength * .5_Kr
      C2 = fractureToughness * internalLength * .5_Kr
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(alphaloc(elemScalType%numDof))
         Do cell = 1,size(cellID)   
            Call SectionRealRestrictClosure(alpha,meshScal,cellID(cell),elemScalType%numDof,alphaLoc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elemScal(cell)%Gauss_C)
               alphaElem     = 0.0_Kr
               gradAlphaElem = 0.0_Kr
               Do iDoF1 = 1,elemScalType%numDof
                  alphaElem     = alphaElem     + alphaLoc(iDof1) * elemScal(cell)%BF(iDoF1,iGauss)
                  gradAlphaElem = gradAlphaElem + alphaLoc(iDof1) * elemScal(cell)%Grad_BF(iDoF1,iGauss)
               End Do
               energy = energy + elemScal(cell)%Gauss_C(iGauss) * (alphaElem**2 * C1 + (gradAlphaElem .dotP. gradAlphaElem) * C2)
            End Do ! Gauss
         End Do ! cell
         flops = (2 * elemScalType%numDof + 6)* size(elemScal(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(alphaloc)
      End If   
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageSurfaceEnergySetAT2


#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageSurfaceEnergySetLinSoft"
!!!
!!!  
!!!  MEF90GradDamageSurfaceEnergySetLinSoft:
!!!  
!!!  (c) 2014-2019 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90GradDamageSurfaceEnergySetLinSoft(energy,alpha,meshScal,cellIS,internalLength,fractureToughness,elemScal,elemScalType,ierr)
      PetscReal,Intent(OUT)                              :: energy
      Type(SectionReal),Intent(IN)                       :: alpha
      Type(DM),Intent(IN)                                :: meshScal
      Type(IS),Intent(IN)                                :: cellIS
      PetscReal                                          :: internalLength,fractureToughness
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elemScal
      Type(MEF90Element_Type),Intent(IN)                 :: elemScalType
      PetscErrorCode,Intent(OUT)                         :: ierr

      PetscReal,Dimension(:),Pointer                     :: alphaLoc
      PetscReal                                          :: alphaElem
      Type(MEF90_VECT)                                   :: gradAlphaElem
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscInt                                           :: iDoF1,iGauss
      PetscReal                                          :: C1, C2
      PetscLogDouble                                     :: flops

      
      C1 = fractureToughness / internalLength * .3183098862_Kr
      C2 = fractureToughness * internalLength * .3183098862_Kr
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(alphaloc(elemScalType%numDof))
         Do cell = 1,size(cellID)   
            Call SectionRealRestrictClosure(alpha,meshScal,cellID(cell),elemScalType%numDof,alphaLoc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elemScal(cell)%Gauss_C)
               alphaElem     = 0.0_Kr
               gradAlphaElem = 0.0_Kr
               Do iDoF1 = 1,elemScalType%numDof
                  alphaElem     = alphaElem     + alphaLoc(iDof1) * elemScal(cell)%BF(iDoF1,iGauss)
                  gradAlphaElem = gradAlphaElem + alphaLoc(iDof1) * elemScal(cell)%Grad_BF(iDoF1,iGauss)
               End Do
               energy = energy + elemScal(cell)%Gauss_C(iGauss) * ((1.0_Kr - (1.0_Kr- alphaElem)**2 )* C1 + (gradAlphaElem .dotP. gradAlphaElem) * C2)
            End Do ! Gauss
         End Do ! cell
         flops = (2 * elemScalType%numDof + 5 )* size(elemScal(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(alphaloc)
      End If   
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageSurfaceEnergySetLinSoft

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


