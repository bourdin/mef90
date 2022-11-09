#include "../MEF90/mef90.inc"
#include "../m_DefMech/mef90DefMech.inc"
Module MEF90_APPEND(m_MEF90_DefMechAssembly,MEF90_DIM)D
#include "petsc/finclude/petsc.h"
#define MEF90_HDRegularization 0.01_Kr

   Use m_MEF90
   Use m_MEF90_DefMechCtx
   Use MEF90_APPEND(m_MEF90_DefMechSplit,MEF90_DIM)D
   Use m_MEF90_DefMechAT
   
   Implicit none
   Private
   Public MEF90DefMechOperatorDisplacement,     &
          MEF90DefMechBilinearFormDisplacement, &
          MEF90DefMechWork,                     &
          MEF90DefMechCohesiveEnergy,           &
          MEF90DefMechPlasticDissipation,       &
          MEF90DefMechElasticEnergy,            &
          MEF90DefMechOperatorDamage,           &
          MEF90DefMechBilinearFormDamage,       &
          MEF90DefMechSurfaceEnergy,            &
          MEF90DefMechCrackVolume,              &
          MEF90DefMechStress

   PetscReal,Parameter,Private                         :: cwLinSoft = 0.7853981634_Kr
Contains


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDisplacement"
!!!
!!!  
!!!  MEF90DefMechOperatorDisplacement: Build the operator. When called in SNES, the solution time should always match the target time, 
!!!                                    so there is no need for interpolation of the forcees, external, and boundary values
!!!  
!!!  (c) 2012-20 Blaise Bourdin bourdin@lsu.edu, Erwan Tanne erwan.tanne@gmail.com
!!!         2022 Alexis Marboeuf, marboeua@mcmaster.ca
!!!

   Subroutine MEF90DefMechOperatorDisplacement(snesDisplacement,displacement,residual,MEF90DefMechCtx,ierr)
      Type(tSNES),Intent(IN)                             :: snesDisplacement
      Type(tVec),Intent(IN)                              :: displacement
      Type(tVec),Intent(INOUT)                           :: residual
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(INOUT)                       :: ierr

   End Subroutine MEF90DefMechOperatorDisplacement

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDisplacement"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDisplacement:
!!!  
!!!  (c) 2012-20 Blaise Bourdin bourdin@lsu.edu,Erwan Tanne erwan.tanne@gmail.com
!!!

   Subroutine MEF90DefMechBilinearFormDisplacement(snesDisplacement,displacement,A,M,flg,MEF90DefMechCtx,ierr)
      Type(tSNES),Intent(IN)                             :: snesDisplacement
      Type(tVec),Intent(IN)                              :: displacement
      Type(tMat),Intent(INOUT)                           :: A,M
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr  
   

   End Subroutine MEF90DefMechBilinearFormDisplacement

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechWork"
!!!
!!!  
!!!  MEF90DefMechWork:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechWork(xVec,MEF90DefMechCtx,work,ierr)
      Type(tVec),Intent(IN)                              :: xVec
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,Dimension(:),Pointer                     :: work
      PetscErrorCode,Intent(OUT)                         :: ierr
   
   End Subroutine MEF90DefMechWork

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCohesiveEnergy"
!!!
!!!  
!!!  MEF90DefMechCohesiveEnergy:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechCohesiveEnergy(xVec,MEF90DefMechCtx,cohesiveEnergy,ierr)
      Type(tVec),Intent(IN)                              :: xVec
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,Dimension(:),Pointer                     :: cohesiveEnergy
      PetscErrorCode,Intent(OUT)                         :: ierr
   
   End Subroutine MEF90DefMechCohesiveEnergy


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechPlasticDissipation"
!!!
!!!  
!!!  MEF90DefMechPlasticDissipation: 
!!!  
!!!  (c) 2015 Erwan TANNE erwan.tanne@gmail.com
!!!

   Subroutine MEF90DefMechPlasticDissipation(x,MEF90DefMechCtx,plasticStrainOld,energy,ierr)
      Type(tVec),Intent(IN)                              :: x
      Type(tVec),Intent(IN)                              :: plasticStrainOld
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,dimension(:),Pointer                     :: energy
      PetscErrorCode,Intent(OUT)                         :: ierr
   
   End Subroutine MEF90DefMechPlasticDissipation



#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechElasticEnergy"
!!!
!!!  
!!!  MEF90DefMechElasticEnergy: 
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechElasticEnergy(x,MEF90DefMechCtx,energy,ierr)
      Type(tVec),Intent(IN)                              :: x
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,dimension(:),Pointer                     :: energy
      PetscErrorCode,Intent(OUT)                         :: ierr
   
   End Subroutine MEF90DefMechElasticEnergy


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechStress"
!!!
!!!  
!!!  MEF90DefMechStress: 
!!!  
!!!  (c) 2012-20 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechStress(displacement,MEF90DefMechCtx,stress,ierr)
      Type(tVec),Intent(IN)                              :: displacement
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      Type(tVec),Intent(IN)                              :: stress
         
   End Subroutine MEF90DefMechStress


! #undef __FUNCT__
! #define __FUNCT__ "MEF90DefMechBilinearFormDamageLinSoftLoc"
! !!!
! !!!  
! !!!  MEF90DefMechBilinearFormDamageLinSoftLoc: Assembles the bilinear form for LinSoft, that is:
! !!!    <K(alpha) beta , delta> = \int a''(\alpha) W(e(u)) beta \delta 
! !!!                            + Gc/pi/\ell w''(alpha) beta \delta
! !!!                            + 2 Gc/pi*\ell \nabla beta \nabla \delta
! !!! with a(\alpha) = (1-w(\alpha)) / ( 1 + (k-1) w(\alpha)) and w(\alpha) = 1 - (1-\alpha)^2
! !!! i.e. a''(\alpha) = 2 k * ( k + 3(k-1) v^2) / [k + (1-k) v^2]^3, w''(\alpha) = -2
! !!!  
! !!!  (c) 2014-2015 Erwan Tanne erwan.tanne@gmail.com, Blaise Bourdin bourdin@lsu.edu
! !!!

!    Subroutine MEF90DefMechBilinearFormDamageLinSoftLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
!       PetscReal,Dimension(:,:),Pointer                   :: Aloc
!       PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
!       Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
!       PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
!       Type(MEF90_MATPROP),Intent(IN)                     :: matprop
!       Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
!       Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

!       PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
!       PetscReal                                          :: elasticEnergyDensityGauss,temperatureGauss,damageGauss
!       Type(MEF90_MATS)                                   :: inelasticStrainGauss
!       PetscReal                                          :: C2,C1,N,k,C0
!       PetscLogDouble                                     :: flops
!       PetscErrorCode                                     :: ierr

!       numDofDisplacement = size(elemDisplacement%BF,1)
!       numDofDamage = size(elemDamage%BF,1)
!       numGauss = size(elemDamage%BF,2)
      
!       k = matprop%CoefficientLinSoft
!       C1 = 2.0_Kr * matprop%fractureToughness / matprop%internalLength / PETSC_PI
!       C2 = 2.0_Kr * matprop%fractureToughness * matprop%internalLength / PETSC_PI
!       N  = matprop%DuctileCouplingPower

!       Do iGauss = 1,numGauss
!          temperatureGauss = 0.0_Kr
!          If (Associated(temperatureDof)) Then
!             Do iDoF1 = 1,numDofDamage
!                temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
!             End Do
!          End If

!          inelasticStrainGauss = 0.0_Kr
!          Do iDoF1 = 1,numDofDisplacement
!             inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
!          End Do
!          inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion - plasticStrainCell
!          elasticEnergyDensityGauss = (matprop%HookesLaw * inelasticStrainGauss) .DotP. inelasticStrainGauss

!          damageGauss = 0.0_Kr
!          Do iDof1 = 1, numDofDamage
!             damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
!          End Do
!          C0 = k * (k + 3.0_Kr * (k - 1.0_Kr) * (1.0_Kr - damageGauss)**2) / &
!               (k + (1.0_Kr - k) * (1.0_Kr - damageGauss)**2)**3
!          !!! This is really twice the elastic energy density

!          Do iDoF1 = 1,numDofDamage
!             Do iDoF2 = 1,numDofDamage
!                Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) * ( & 
!                                   (elasticEnergyDensityGauss * C0 - C1 + &
!                                   (N * (N - 1.0) *( (1.0_Kr-damageGauss)**(N - 2.0) ) ) * cumulatedDissipatedPlasticEnergyCell ) * elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss) &
!                                   + C2 * ((matprop%toughnessAnisotropyMatrix * elemDamage%Grad_BF(iDoF1,iGauss)) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
!             End Do
!          End Do
!       End Do
!       !flops = 2 * numGauss * numDofDisplacement**2
!       !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
!    End Subroutine MEF90DefMechBilinearFormDamageLinSoftLoc

! #undef __FUNCT__
! #define __FUNCT__ "MEF90DefMechOperatorDamageLinSoftLoc"
! !!!
! !!!  
! !!!  MEF90DefMechOperatorDamageLinSoftLoc: Assembles the operator for LinSoft:
! !!!    F(\alpha)beta = \int [k(\alpha-1)W(e(u))] / [k-(k-1)(1-\alpha)^2]^2 beta 
! !!!                   + 2Gc/pi (1-\alpha) beta / \ell 
! !!!                   + 2Gc/pi \nabla \alpha \nabla beta * \ell
! !!!  
! !!!  (c) 2014-2015 Erwan Tanne erwan.tanne@gmail.com, Blaise Bourdin bourdin@lsu.edu
! !!!

!    Subroutine MEF90DefMechOperatorDamageLinSoftLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
!       PetscReal,Dimension(:),Pointer                     :: residualLoc
!       PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
!       Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
!       PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
!       Type(MEF90_MATPROP),Intent(IN)                     :: matprop
!       Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
!       Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
!       PetscReal,Intent(IN)                               :: CrackPressureCell

!       PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
!       PetscReal                                          :: elasticEnergyDensityGauss,temperatureGauss,damageGauss
!       Type(MEF90_MATS)                                   :: inelasticStrainGauss
!       Type(MEF90_VECT)                                   :: gradientDamageGauss
!       PetscReal                                          :: C1,C2,k,C0Operator,N
!       PetscLogDouble                                     :: flops
!       PetscErrorCode                                     :: ierr


!       numDofDisplacement = size(elemDisplacement%BF,1)
!       numDofDamage = size(elemDamage%BF,1)
!       numGauss = size(elemDamage%BF,2)
      
!       C1 = 2.0_Kr * matprop%fractureToughness / matprop%internalLength / PETSC_PI
!       C2 = 2.0_Kr * matprop%fractureToughness * matprop%internalLength / PETSC_PI
!       k = matprop%CoefficientLinSoft
!       N = matprop%DuctileCouplingPower

!       Do iGauss = 1,numGauss
!          temperatureGauss = 0.0_Kr
!          If (Associated(temperatureDof)) Then
!             Do iDoF1 = 1,numDofDamage
!                temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
!             End Do
!          End If

!          inelasticStrainGauss = 0.0_Kr
!          Do iDoF1 = 1,numDofDisplacement
!             inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
!          End Do
!          inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion - plasticStrainCell
!          elasticEnergyDensityGauss = (matprop%HookesLaw * inelasticStrainGauss) .DotP. inelasticStrainGauss
!          !!! This is really twice the elastic energy density

!          damageGauss = 0.0_Kr
!          gradientDamageGauss = 0.0_Kr
!          Do iDoF1 = 1,numDofDamage
!             damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
!             gradientDamageGauss = gradientDamageGauss + elemDamage%Grad_BF(iDoF1,iGauss) * xDof(iDoF1)
!          End Do
!          C0Operator = k * (damageGauss - 1.0_Kr) * elasticEnergyDensityGauss / &
!                      (k + (1.0_Kr - k) * (1.0_kr - damageGauss)**2)**2

!          Do iDoF2 = 1,numDofDamage
!             residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * ( &
!                                  (C0Operator + C1 * (1.0_Kr - damageGauss) &
!                                  - ( N * cumulatedDissipatedPlasticEnergyCell * (1.0_Kr - damageGauss)**(DBLE(N - 1.0))  ) ) * elemDamage%BF(iDoF2,iGauss) & 
!                                   + C2 * ((matprop%toughnessAnisotropyMatrix * gradientDamageGauss) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
!          End Do
!       End Do
!       !flops = 2 * numGauss * numDofDisplacement**2
!       !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
!    End Subroutine MEF90DefMechOperatorDamageLinSoftLoc



#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamage"
!!!
!!!  
!!!  MEF90DefMechOperatorDamage: Build the operator. When called in SNES, the solution time should always match the target time, 
!!!                                    so there is no need for interpolation of the forcees, external, and boundary values
!!!  
!!!  (c) 2012-20 Blaise Bourdin bourdin@lsu.edu, Erwan Tanne erwan.tanne@gmail.com
!!!

   Subroutine  MEF90DefMechOperatorDamage(snesDamage,damage,residual,MEF90DefMechCtx,ierr)
      Type(tSNES),Intent(IN)                             :: snesDamage
      Type(tVec),Intent(IN)                              :: damage
      Type(tVec),Intent(INOUT)                           :: residual
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
   
   End Subroutine MEF90DefMechOperatorDamage

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamage"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamage:
!!!  
!!!  (c) 2012-19 Blaise Bourdin bourdin@lsu.edu, Erwan Tanne erwan.tanne@gmail.com
!!!

   Subroutine MEF90DefMechBilinearFormDamage(snesDamage,damage,A,M,flg,MEF90DefMechCtx,ierr)
      Type(tSNES),Intent(IN)                             :: snesDamage
      Type(tVec),Intent(IN)                              :: damage
      Type(tMat),Intent(INOUT)                           :: A,M
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr  
      
   End Subroutine MEF90DefMechBilinearFormDamage


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSurfaceEnergy"
!!!
!!!  
!!!  MEF90DefMechSurfaceEnergy:
!!!  
!!!  (c) 2014-2020 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechSurfaceEnergy(Damage,MEF90DefMechCtx,energy,ierr)
      Type(tVec),Intent(IN)                              :: Damage
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,dimension(:),Pointer                     :: energy
      PetscErrorCode,Intent(OUT)                         :: ierr
   
   End Subroutine MEF90DefMechSurfaceEnergy


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCrackVolume"
!!!
!!!  
!!!  MEF90DefMechCrackVolume:
!!!  
!!!  (c) 2016-2021 Erwan erwan.tanne@gmail.com, Blaise Bourdin bourdin@lsu.edu  
!!!

   Subroutine MEF90DefMechCrackVolume(xVec,MEF90DefMechCtx,CrackVolume,ierr)
      Type(tVec),Intent(IN)                              :: xVec
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,Dimension(:),Pointer                     :: CrackVolume
      PetscErrorCode,Intent(OUT)                         :: ierr
   
   End Subroutine MEF90DefMechCrackVolume
End Module MEF90_APPEND(m_MEF90_DefMechAssembly,MEF90_DIM)D
