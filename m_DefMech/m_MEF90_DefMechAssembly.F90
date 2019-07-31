#include "../MEF90/mef90.inc"
Module MEF90_APPEND(m_MEF90_DefMechAssembly,MEF90_DIM)D
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
#define MEF90_HDRegularization 0.01_Kr

   Use m_MEF90
   Use m_MEF90_DefMechCtx
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
          MEF90DefMechCrackPressureRescaling,   &
          MEF90DefMechStress

   Abstract Interface
      Subroutine MEF90DefMechBilinearFormLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
         Use m_MEF90_DefMechCtx
         PetscReal,Dimension(:,:),Pointer                   :: Aloc
         PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
         Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
         PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
         Type(MEF90_MATPROP),Intent(IN)                     :: matprop
         Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
         Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      End Subroutine MEF90DefMechBilinearFormLoc
   End Interface

   Abstract Interface
      Subroutine MEF90DefMechOperatorLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
         Use m_MEF90_DefMechCtx
         PetscReal,Dimension(:),Pointer                     :: residualLoc
         PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
         Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
         PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
         Type(MEF90_MATPROP),Intent(IN)                     :: matprop
         Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
         Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
         PetscReal,Intent(IN)                               :: CrackPressureCell
      End Subroutine MEF90DefMechOperatorLoc
   End Interface

   Abstract Interface
      Subroutine MEF90DefMechRHSLoc(residualLoc,xDof,forceCell,pressureCell,matprop,elemDisplacement,damageDof,elemDamage,CrackPressureCell)
         Use m_MEF90_DefMechCtx
         PetscReal,Dimension(:),Pointer                     :: residualLoc
         PetscReal,Dimension(:),Pointer                     :: xDof,damageDof
         Type(MEF90_VECT),Intent(IN)                        :: forceCell
         PetscReal,Intent(IN)                               :: pressureCell
         Type(MEF90_MATPROP),Intent(IN)                     :: matprop
         Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
         Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
         PetscReal,Intent(IN)                               :: CrackPressureCell
      End Subroutine MEF90DefMechRHSLoc
   End Interface
   
   Abstract Interface
      Subroutine MEF90DefMechEnergyLoc(energyLoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
         Use m_MEF90_DefMechCtx
         PetscReal,Dimension(:),Pointer                     :: residualLoc
         PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
         Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
         PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
         Type(MEF90_MATPROP),Intent(IN)                     :: matprop
         Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
         Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      End Subroutine MEF90DefMechEnergyLoc
   End Interface

   PetscReal,Parameter,Private                         :: cwKKL = 0.7165753016381484_Kr
Contains
#undef __FUNCT__
#define __FUNCT__ "gKKL"
!!!
!!!  
!!!  MEF90gKKL: g(\alpha) = 4(1-\alpha)^3 - 3 (1-\alpha)^4
!!!  The KKL (v2) model can be rewritten as a regular GD model with 
!!!  a(\alpha) = g(\alpha) and w(\alpha) = 1-g(\alpha)
!!!  
!!!  (c) 2019 Blaise Bourdin bourdin@lsu.edu
!!!

   Function MEF90gKKL(s)
      PetscReal,Intent(IN)                             :: s
      PetscReal                                        :: MEF90gKKL
      PetscLogDouble                                   :: flops
      PetscErrorCode                                   :: ierr

      MEF90gKKL = 4.0_Kr * (1.0_Kr - s)**3 - 3.0_Kr * (1.0_Kr - s)**4
      flops = 7.0
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Function MEF90gKKL

#undef __FUNCT__
#define __FUNCT__ "DgKKL"
!!!
!!!  
!!!  MEF90DgKKL: g'(\alpha) = -12(1-\alpha)^2 + 12 (1-\alpha)^3
!!!  The KKL (v2) model can be rewritten as a regular GD model with 
!!!  a(\alpha) = g(\alpha) and w(\alpha) = 1-g(\alpha)
!!!  
!!!  (c) 2019 Blaise Bourdin bourdin@lsu.edu
!!!

   Function MEF90DgKKL(s)
      PetscReal,Intent(IN)                             :: s
      PetscReal                                        :: MEF90DgKKL
      PetscLogDouble                                   :: flops
      PetscErrorCode                                   :: ierr

      MEF90DgKKL = -12.0_Kr * ((1.0_Kr - s)**2 - (1.0_Kr - s)**3)
      flops = 6.0
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Function MEF90DgKKL

#undef __FUNCT__
#define __FUNCT__ "D2gKKL"
!!!
!!!  
!!!  MEF90D2gKKL: g''(\alpha) = 24(1-\alpha)^2 - 36 (1-\alpha)^3
!!!  The KKL (v2) model can be rewritten as a regular GD model with 
!!!  a(\alpha) = g(\alpha) and w(\alpha) = 1-g(\alpha)
!!!  
!!!  (c) 2019 Blaise Bourdin bourdin@lsu.edu
!!!

   Function MEF90D2gKKL(s)
      PetscReal,Intent(IN)                             :: s
      PetscReal                                        :: MEF90D2gKKL
      PetscLogDouble                                   :: flops
      PetscErrorCode                                   :: ierr

      MEF90D2gKKL = 24.0_Kr * (1.0_Kr - s) - 36.0_Kr * (1.0_Kr - s)**2
      flops = 6.0
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Function MEF90D2gKKL

#undef __FUNCT__
#define __FUNCT__ "MEF90HydrostaticDeviatoricPenaltyFunction"
!!!
!!!  
!!!  MEF90HydrostaticDeviatoricPenaltyFunction: a low order polynomial C^2 regularization of (max(x,0))**2, defined by
!!!       0                         if x \le 0
!!!       x^3/3/\gamma             if 0 < x \le \gamma
!!!       (x-\gamma/2)^2+gamma^2/12 otherwise
!!!    
!!!  (c) 2018 Blaise Bourdin bourdin@lsu.edu
!!!

   Pure Function MEF90HydrostaticDeviatoricPenaltyFunction(x)
      PetscReal,Intent(IN) :: x
      PetscReal            :: MEF90HydrostaticDeviatoricPenaltyFunction

      If (x <= 0.0_Kr) Then
         MEF90HydrostaticDeviatoricPenaltyFunction = 0.0_Kr
      Else if (x <= MEF90_HDRegularization) Then
         MEF90HydrostaticDeviatoricPenaltyFunction = x**3/3.0_Kr / MEF90_HDRegularization
      Else
         MEF90HydrostaticDeviatoricPenaltyFunction = (x-MEF90_HDRegularization/2.0_Kr)**2 + MEF90_HDRegularization**2/12.0_Kr
      End If
   End Function MEF90HydrostaticDeviatoricPenaltyFunction

#undef __FUNCT__
#define __FUNCT__ "MEF90HydrostaticDeviatoricDPenaltyFunction"
!!!
!!!  
!!!  MEF90HydrostaticDeviatoricDPenaltyFunction: the first derivative of MEF90HydrostaticDeviatoricPenaltyFunction
!!!    
!!!  (c) 2018 Blaise Bourdin bourdin@lsu.edu
!!!

   Pure Function MEF90HydrostaticDeviatoricDPenaltyFunction(x)
      PetscReal,Intent(IN) :: x
      PetscReal            :: MEF90HydrostaticDeviatoricDPenaltyFunction

      If (x <= 0.0_Kr) Then
         MEF90HydrostaticDeviatoricDPenaltyFunction = 0.0_Kr
      Else if (x <= MEF90_HDRegularization) Then
         MEF90HydrostaticDeviatoricDPenaltyFunction = x**2 / MEF90_HDRegularization
      Else
         MEF90HydrostaticDeviatoricDPenaltyFunction = 2.0_Kr * x - MEF90_HDRegularization
      End If
   End Function MEF90HydrostaticDeviatoricDPenaltyFunction

#undef __FUNCT__
#define __FUNCT__ "MEF90HydrostaticDeviatoricD2PenaltyFunction"
!!!
!!!  
!!!  MEF90HydrostaticDeviatoricD2PenaltyFunction: the second derivative of MEF90HydrostaticDeviatoricPenaltyFunction
!!!    
!!!  (c) 2018 Blaise Bourdin bourdin@lsu.edu
!!!

   Pure Function MEF90HydrostaticDeviatoricD2PenaltyFunction(x)
      PetscReal,Intent(IN) :: x
      PetscReal            :: MEF90HydrostaticDeviatoricD2PenaltyFunction

      If (x <= 0.0_Kr) Then
         MEF90HydrostaticDeviatoricD2PenaltyFunction = 0.0_Kr
      Else if (x <= MEF90_HDRegularization) Then
         MEF90HydrostaticDeviatoricD2PenaltyFunction = 2.0_Kr * x / MEF90_HDRegularization
      Else
         MEF90HydrostaticDeviatoricD2PenaltyFunction = 2.0_Kr 
      End If
   End Function MEF90HydrostaticDeviatoricD2PenaltyFunction



#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorNull"
!!!
!!!  
!!!  MEF90DefMechOperatorNull:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorNull(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
      Use m_MEF90_DefMechCtx
      PetscReal,Dimension(:),Pointer                        :: residualLoc
      PetscReal,Dimension(:),Pointer                        :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
         Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
         PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
         Type(MEF90_MATPROP),Intent(IN)                     :: matprop
         Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
         Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
         PetscReal,Intent(IN)                               :: CrackPressureCell

      Print*,__FUNCT__,': Unimplemented Operator local assembly function'
      STOP
   End Subroutine MEF90DefMechOperatorNull


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormNull"
!!!
!!!  
!!!  MEF90DefMechBilinearFormNull:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormNull(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      
      Print*,__FUNCT__,': Unimplemented Bilinear Form local assembly function'
      STOP 
   End Subroutine MEF90DefMechBilinearFormNull
   
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDisplacementElasticLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDisplacementElasticLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDisplacementElasticLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numGauss
      Type(MEF90_MATS)                                   :: sigma
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numGauss = size(elemDisplacement%BF,2)
      !Aloc = 0.0_Kr
      Do iGauss = 1,numGauss
         Do iDoF1 = 1,numDofDisplacement
            sigma = matProp%HookesLaw * elemDisplacement%GradS_BF(iDoF1,iGauss) 
            Do iDoF2 = 1,numDofDisplacement
               ALoc(iDoF2,iDoF1) = ALoc(iDoF2,iDoF1) + elemDisplacement%Gauss_C(iGauss) * (sigma .DotP. elemDisplacement%GradS_BF(iDoF2,iGauss))
               If (matprop%cohesiveStiffness /= 0.0_Kr) Then
                  ALoc(iDoF2,iDoF1) = ALoc(iDoF2,iDoF1) + elemDisplacement%Gauss_C(iGauss) * matprop%cohesivestiffness * (elemDisplacement%BF(iDoF1,iGauss) .DotP. elemDisplacement%BF(iDoF2,iGauss))
               End If            
            End Do
         End Do
      End Do
      flops = 2 * numGauss * numDofDisplacement**2
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDisplacementElasticLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDisplacementATLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDisplacementATLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDisplacementATLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: ALoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: stiffness
      Type(MEF90_MATS)                                   :: sigma
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDisplacement%BF,2)
      !ALoc = 0.0_Kr
      Do iGauss = 1,numGauss
         stiffness = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            stiffness = stiffness + damageDof(iDoF1) * elemDamage%BF(iDoF1,iGauss)
         End Do
         stiffness = (1.0_Kr - stiffness)**2 + matProp%residualStiffness
         Do iDoF1 = 1,numDofDisplacement
            sigma = stiffness * (matProp%HookesLaw * elemDisplacement%GradS_BF(iDoF1,iGauss))
            Do iDoF2 = 1,numDofDisplacement
               ALoc(iDoF2,iDoF1) = ALoc(iDoF2,iDoF1) + elemDisplacement%Gauss_C(iGauss) * (sigma .DotP. elemDisplacement%GradS_BF(iDoF2,iGauss))
               If (matprop%cohesiveStiffness /= 0.0_Kr) Then
                  ALoc(iDoF2,iDoF1) = ALoc(iDoF2,iDoF1) + elemDisplacement%Gauss_C(iGauss) * matprop%cohesivestiffness * (elemDisplacement%BF(iDoF1,iGauss) .DotP. elemDisplacement%BF(iDoF2,iGauss))
               End If            
            End Do
         End Do
      End Do
      flops = numGauss * ( 4. * numDofDisplacement**2 + 2. * numDofDamage + 3.)
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDisplacementATLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDisplacementKKLLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDisplacementKKLLoc:
!!!  
!!!  (c) 2019 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDisplacementKKLLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: ALoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: stiffness,damageGauss
      Type(MEF90_MATS)                                   :: sigma
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDisplacement%BF,2)
      !ALoc = 0.0_Kr
      Do iGauss = 1,numGauss
         damageGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            damageGauss = damageGauss + damageDof(iDoF1) * damageDof(iDof1) * elemDamage%BF(iDoF1,iGauss)
         End Do
         stiffness = MEF90gKKL(damageGauss) + matProp%residualStiffness
         Do iDoF1 = 1,numDofDisplacement
            sigma = stiffness * (matProp%HookesLaw * elemDisplacement%GradS_BF(iDoF1,iGauss))
            Do iDoF2 = 1,numDofDisplacement
               ALoc(iDoF2,iDoF1) = ALoc(iDoF2,iDoF1) + elemDisplacement%Gauss_C(iGauss) * (sigma .DotP. elemDisplacement%GradS_BF(iDoF2,iGauss))
               If (matprop%cohesiveStiffness /= 0.0_Kr) Then
                  ALoc(iDoF2,iDoF1) = ALoc(iDoF2,iDoF1) + elemDisplacement%Gauss_C(iGauss) * matprop%cohesivestiffness * (elemDisplacement%BF(iDoF1,iGauss) .DotP. elemDisplacement%BF(iDoF2,iGauss))
               End If            
            End Do
         End Do
      End Do
      flops = numGauss * ( 4. * numDofDisplacement**2 + 2. * numDofDamage + 3.)
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDisplacementKKLLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDisplacementLinSoftLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDisplacementElasticLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDisplacementLinSoftLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: ALoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: stiffness
      Type(MEF90_MATS)                                   :: sigma
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr
      PetscReal                                          :: k

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDisplacement%BF,2)
      !ALoc = 0.0_Kr
      k = matprop%CoefficientLinSoft
      Do iGauss = 1,numGauss
         stiffness = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            stiffness = stiffness + damageDof(iDoF1) * elemDamage%BF(iDoF1,iGauss)
         End Do
         stiffness = (1.0_Kr - stiffness)**2 /(1.0_Kr + ( k - 1.0_Kr ) * ( 1.0_Kr - (1.0_Kr - stiffness)**2 ) ) + matProp%residualStiffness
         Do iDoF1 = 1,numDofDisplacement
            sigma = stiffness * (matProp%HookesLaw * elemDisplacement%GradS_BF(iDoF1,iGauss))
            Do iDoF2 = 1,numDofDisplacement
               ALoc(iDoF2,iDoF1) = ALoc(iDoF2,iDoF1) + elemDisplacement%Gauss_C(iGauss) * (sigma .DotP. elemDisplacement%GradS_BF(iDoF2,iGauss))
               If (matprop%cohesiveStiffness /= 0.0_Kr) Then
                  ALoc(iDoF2,iDoF1) = ALoc(iDoF2,iDoF1) + elemDisplacement%Gauss_C(iGauss) * matprop%cohesivestiffness * (elemDisplacement%BF(iDoF1,iGauss) .DotP. elemDisplacement%BF(iDoF2,iGauss))
               End If            
            End Do
         End Do
      End Do
      !flops = numGauss * ( 2. * numDofDisplacement**2 + 2. * numDofDamage + 3.)
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDisplacementLinSoftLoc
   
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDisplacementATUnilateralHDLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDisplacementATUnilateralHDLoc:
!!!  
!!!  (c) 2014 - 2018 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDisplacementATUnilateralHDLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: ALoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: stiffnessD
      PetscReal                                          :: inelasticStrainTrace
      Type(MEF90_MATS)                                   :: sigmaD
      PetscReal                                          :: sigmaH
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr
      PetscReal                                          :: kappa

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDisplacement%BF,2)
      kappa = ((matProp%HookesLaw * MEF90_MATS_IDENTITY) .dotP. MEF90_MATS_IDENTITY) / MEF90_DIM**2
      !!! This assumes that the Hooke's law is such that it preserves hydrostatic / deviatoric orthogonality

      !ALoc = 0.0_Kr
      Do iGauss = 1,numGauss
         inelasticStrainTrace = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               inelasticStrainTrace = inelasticStrainTrace - temperatureDof(iDoF1) * elemDamage%BF(iDoF1,iGauss)
            End Do         
            inelasticStrainTrace = inelasticStrainTrace * trace(matProp%LinearThermalExpansion)
         End If
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrainTrace = inelasticStrainTrace + xDof(iDoF1) * trace(elemDisplacement%GradS_BF(iDoF1,iGauss))
         End Do
         inelasticStrainTrace = inelasticStrainTrace - Trace(plasticStrainCell)

         stiffnessD = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            stiffnessD = stiffnessD + damageDof(iDoF1) * elemDamage%BF(iDoF1,iGauss)
         End Do
         stiffnessD = (1.0_Kr - stiffnessD)**2 + matProp%residualStiffness

         Do iDoF1 = 1,numDofDisplacement
            sigmaD = stiffnessD * (matProp%HookesLaw * DeviatoricPart(elemDisplacement%GradS_BF(iDoF1,iGauss)))
            sigmaH = kappa * (stiffnessD * MEF90HydrostaticDeviatoricD2PenaltyFunction(inelasticStrainTrace)   &
                                         + MEF90HydrostaticDeviatoricD2PenaltyFunction(-inelasticStrainTrace)) &
                           * Trace(elemDisplacement%GradS_BF(iDoF1,iGauss)) * 0.5_Kr
            Do iDoF2 = 1,numDofDisplacement
               ALoc(iDoF2,iDoF1) = ALoc(iDoF2,iDoF1) + elemDisplacement%Gauss_C(iGauss) * ((sigmaD .DotP. elemDisplacement%GradS_BF(iDoF2,iGauss)) + &
                                                                                           sigmaH * Trace(elemDisplacement%GradS_BF(iDoF2,iGauss)))
            End Do
         End Do
      End Do
      !flops = numGauss * ( 2. * numDofDisplacement**2 + 3. * numDofDamage + 2.)
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDisplacementATUnilateralHDLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDisplacementATUnilateralDeviatoricLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDisplacementATUnilateralDeviatoricLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDisplacementATUnilateralDeviatoricLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: ALoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: stiffnessD
      PetscReal                                          :: inelasticStrainTrace
      Type(MEF90_MATS)                                   :: plasticStrain,sigma
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDisplacement%BF,2)

      !ALoc = 0.0_Kr
      Do iGauss = 1,numGauss
         inelasticStrainTrace = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               inelasticStrainTrace = inelasticStrainTrace - temperatureDof(iDoF1) * elemDamage%BF(iDoF1,iGauss)
            End Do         
            inelasticStrainTrace = inelasticStrainTrace * trace(matProp%LinearThermalExpansion)
         End If
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrainTrace = inelasticStrainTrace + xDof(iDoF1) * trace(elemDisplacement%GradS_BF(iDoF1,iGauss))
         End Do
         stiffnessD = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            stiffnessD = stiffnessD + damageDof(iDoF1) * elemDamage%BF(iDoF1,iGauss)
         End Do
         stiffnessD = (1.0_Kr - stiffnessD)**2 + matProp%residualStiffness

         Do iDoF1 = 1,numDofDisplacement
            sigma = hydrostaticPart(elemDisplacement%GradS_BF(iDoF1,iGauss)) + stiffnessD * deviatoricPart(elemDisplacement%GradS_BF(iDoF1,iGauss))
            sigma = matProp%HookesLaw * sigma
            Do iDoF2 = 1,numDofDisplacement
               ALoc(iDoF2,iDoF1) = ALoc(iDoF2,iDoF1) + elemDisplacement%Gauss_C(iGauss) * (sigma .DotP. elemDisplacement%GradS_BF(iDoF2,iGauss))
            End Do
         End Do
      End Do
      !flops = numGauss * ( 2. * numDofDisplacement**2 + 3. * numDofDamage + 2.)
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDisplacementATUnilateralDeviatoricLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDisplacementATUnilateralPHDLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDisplacementATUnilateralPHDLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDisplacementATUnilateralPHDLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: ALoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: i,iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: stiffness
      PetscReal                                          :: inelasticStrainTrace
      Type(MEF90_MATS)                                   :: plasticStrain,sigma
      PetscReal, Dimension(MEF90_DIM)                    :: PpalStrain
      Type(MEF90_MATS),Dimension(MEF90_DIM)              :: PpalDirection
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDisplacement%BF,2)

      !ALoc = 0.0_Kr
      Do iGauss = 1,numGauss
         inelasticStrainTrace = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               inelasticStrainTrace = inelasticStrainTrace - temperatureDof(iDoF1) * elemDamage%BF(iDoF1,iGauss)
            End Do         
            inelasticStrainTrace = inelasticStrainTrace * trace(matProp%LinearThermalExpansion)
         End If
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrainTrace = inelasticStrainTrace + xDof(iDoF1) * trace(elemDisplacement%GradS_BF(iDoF1,iGauss))
         End Do
         If (inelasticStrainTrace < 0.0_Kr) Then
            stiffness = 1.0_Kr
         Else
            stiffness = 0.0_Kr
            Do iDoF1 = 1,numDofDamage
               stiffness = stiffness + damageDof(iDoF1) * elemDamage%BF(iDoF1,iGauss)
            End Do
            stiffness = (1.0_Kr - stiffness)**2 + matProp%residualStiffness
         End If
         
         Do iDoF1 = 1,numDofDisplacement
            sigma = stiffness * elemDisplacement%GradS_BF(iDoF1,iGauss)
            sigma = matProp%HookesLaw * sigma
            Do iDoF2 = 1,numDofDisplacement
               ALoc(iDoF2,iDoF1) = ALoc(iDoF2,iDoF1) + elemDisplacement%Gauss_C(iGauss) * (sigma .DotP. elemDisplacement%GradS_BF(iDoF2,iGauss))
            End Do
         End Do
      End Do
      !flops = numGauss * ( 2. * numDofDisplacement**2 + 3. * numDofDamage + 2.)
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDisplacementATUnilateralPHDLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDisplacementATUnilateralPSLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDisplacementATUnilateralPSLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDisplacementATUnilateralPSLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: ALoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: i,iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: stiffness
      Type(MEF90_MATS)                                   :: sigma,inelasticStrain
      PetscReal, Dimension(MEF90_DIM)                    :: PpalStrain
      Type(MEF90_MATS),Dimension(MEF90_DIM)              :: PpalDirection
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDisplacement%BF,2)

      !ALoc = 0.0_Kr
      Do iGauss = 1,numGauss
         stiffness = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            stiffness = stiffness + damageDof(iDoF1) * elemDamage%BF(iDoF1,iGauss)
         End Do
         stiffness = (1.0_Kr - stiffness + matProp%residualStiffness)**2
   
         inelasticStrain = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               inelasticStrain = inelasticStrain - temperatureDof(iDoF1) * elemDamage%BF(iDoF1,iGauss) * matProp%LinearThermalExpansion
            End Do         
         End If
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrain = inelasticStrain + xDof(iDoF1) * elemDisplacement%GradS_BF(iDoF1,iGauss)
         End Do
         inelasticStrain = inelasticStrain - plasticStrainCell

         Call SpectralDecomposition(inelasticStrain,PpalStrain,PpalDirection)
         Do iDoF1 = 1,numDofDisplacement
            sigma = 0.0_Kr
            Do i = 1,MEF90_DIM 
               If (PpalStrain(i) < 0.0_Kr) Then
                  sigma = sigma + (PpalStrain(i) * PpalDirection(i))
               Else
                  sigma = sigma + Stiffness * (PpalStrain(i) * PpalDirection(i))
               End If
            End Do
            sigma = matProp%HookesLaw * sigma
            Do iDoF2 = 1,numDofDisplacement
               ALoc(iDoF2,iDoF1) = ALoc(iDoF2,iDoF1) + elemDisplacement%Gauss_C(iGauss) * (sigma .DotP. elemDisplacement%GradS_BF(iDoF2,iGauss))
            End Do
         End Do
      End Do
      !flops = numGauss * ( 2. * numDofDisplacement**2 + 3. * numDofDamage + 2.)
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDisplacementATUnilateralPSLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDisplacementElasticLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDisplacementElasticLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDisplacementElasticLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      PetscReal,Intent(IN)                               :: CrackPressureCell

      Type(MEF90_MATS)                                   :: sigma
      Type(MEF90_VECT)                                   :: UU0
      PetscReal                                          :: temperature
      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage       = size(elemDamage%BF,1)
      numGauss = size(elemDisplacement%BF,2)
      residualLoc = 0.0_Kr
      Do iGauss = 1,numGauss
         sigma = 0.0_Kr
         Do iDoF1 = 1,numDofDisplacement
            sigma = sigma + xDof(iDoF1) * elemDisplacement%GradS_BF(iDoF1,iGauss)
         End Do

         temperature = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperature = temperature + temperatureDoF(iDoF1) * elemDamage%BF(iDoF1,iGauss)
            End Do
         End If

         sigma = matProp%HookesLaw * (sigma - temperature * matProp%LinearThermalExpansion - plasticStrainCell) 
         Do iDoF2 = 1,numDofDisplacement
            residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement%Gauss_C(iGauss) * (sigma .DotP. elemDisplacement%GradS_BF(iDoF2,iGauss))
         End Do
   
         If (Associated(boundaryDisplacementDof) .AND. (matprop%cohesiveStiffness /= 0.0_Kr)) Then
            UU0 = 0.0_Kr
            Do iDoF1 = 1,numDofDisplacement
               UU0 = UU0 + (xDof(iDoF1) - boundaryDisplacementDof(iDoF1)) * elemDisplacement%BF(iDoF1,iGauss)
            End Do
            UU0 = UU0 * matprop%cohesiveStiffness
            Do iDoF2 = 1,numDofDisplacement
               residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement%Gauss_C(iGauss) * (UU0 .DotP. elemDisplacement%BF(iDoF2,iGauss))
            End Do
         End If
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDisplacementElasticLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDisplacementATLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDisplacementATLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDisplacementATLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      PetscReal,Intent(IN)                               :: CrackPressureCell

      Type(MEF90_MATS)                                   :: sigma
      Type(MEF90_VECT)                                   :: UU0
      PetscReal                                          :: stiffness,temperature
      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDisplacement%BF,2)

      !residualLoc = 0.0_Kr
      Do iGauss = 1,numGauss
         stiffness = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            stiffness = stiffness + damageDof(iDoF1) * elemDamage%BF(iDoF1,iGauss)
         End Do
         stiffness = (1.0_Kr - stiffness)**2 + matProp%residualStiffness

         sigma = 0.0_Kr
         Do iDoF1 = 1,numDofDisplacement
            sigma = sigma + xDof(iDoF1) * elemDisplacement%GradS_BF(iDoF1,iGauss)
         End Do

         temperature = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperature = temperature + temperatureDoF(iDoF1) * elemDamage%BF(iDoF1,iGauss)
            End Do
         End If


         sigma = stiffness * (matProp%HookesLaw * (sigma - temperature * matProp%LinearThermalExpansion - plasticStrainCell) )

#if MEF90_DIM == 2
         !!! Adding terms in planestrain for plasticity with tr(p) = 0
         If (.NOT. matProp%HookesLaw%isPlaneStress) Then
            sigma = sigma +  stiffness * ( matProp%HookesLaw%lambda*trace(plasticStrainCell)*MEF90MatS2DIdentity )
         End If
#endif

         Do iDoF2 = 1,numDofDisplacement
            residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement%Gauss_C(iGauss) * (sigma .DotP. elemDisplacement%GradS_BF(iDoF2,iGauss))
         End Do
         If (Associated(boundaryDisplacementDof) .AND. (matprop%cohesiveStiffness /= 0.0_Kr)) Then
            UU0 = 0.0_Kr
            Do iDoF1 = 1,numDofDisplacement
               UU0 = UU0 + (xDof(iDoF1) - boundaryDisplacementDof(iDoF1)) * elemDisplacement%BF(iDoF1,iGauss)
            End Do
            UU0 = UU0 * matprop%cohesiveStiffness
            Do iDoF2 = 1,numDofDisplacement
               residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement%Gauss_C(iGauss) * (UU0 .DotP. elemDisplacement%BF(iDoF2,iGauss))
            End Do
         End If
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDisplacementATLoc


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDisplacementKKLLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDisplacementKKLLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDisplacementKKLLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      PetscReal,Intent(IN)                               :: CrackPressureCell

      Type(MEF90_MATS)                                   :: sigma
      Type(MEF90_VECT)                                   :: UU0
      PetscReal                                          :: stiffness,temperature,damageGauss
      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDisplacement%BF,2)

      !residualLoc = 0.0_Kr
      Do iGauss = 1,numGauss
         damageGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            damageGauss = damageGauss + damageDof(iDoF1) * damageDof(iDof1) * elemDamage%BF(iDoF1,iGauss)
         End Do
         stiffness = MEF90gKKL(damageGauss) + matProp%residualStiffness

         sigma = 0.0_Kr
         Do iDoF1 = 1,numDofDisplacement
            sigma = sigma + xDof(iDoF1) * elemDisplacement%GradS_BF(iDoF1,iGauss)
         End Do

         temperature = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperature = temperature + temperatureDoF(iDoF1) * elemDamage%BF(iDoF1,iGauss)
            End Do
         End If


         sigma = stiffness * (matProp%HookesLaw * (sigma - temperature * matProp%LinearThermalExpansion - plasticStrainCell) )

#if MEF90_DIM == 2
         !!! Adding terms in planestrain for plasticity with tr(p) = 0
         If (.NOT. matProp%HookesLaw%isPlaneStress) Then
            sigma = sigma +  stiffness * ( matProp%HookesLaw%lambda*trace(plasticStrainCell)*MEF90MatS2DIdentity )
         End If
#endif

         Do iDoF2 = 1,numDofDisplacement
            residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement%Gauss_C(iGauss) * (sigma .DotP. elemDisplacement%GradS_BF(iDoF2,iGauss))
         End Do
         If (Associated(boundaryDisplacementDof) .AND. (matprop%cohesiveStiffness /= 0.0_Kr)) Then
            UU0 = 0.0_Kr
            Do iDoF1 = 1,numDofDisplacement
               UU0 = UU0 + (xDof(iDoF1) - boundaryDisplacementDof(iDoF1)) * elemDisplacement%BF(iDoF1,iGauss)
            End Do
            UU0 = UU0 * matprop%cohesiveStiffness
            Do iDoF2 = 1,numDofDisplacement
               residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement%Gauss_C(iGauss) * (UU0 .DotP. elemDisplacement%BF(iDoF2,iGauss))
            End Do
         End If
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDisplacementKKLLoc


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDisplacementLinSoftLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDisplacementLinSoftLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDisplacementLinSoftLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      PetscReal,Intent(IN)                               :: CrackPressureCell

      Type(MEF90_MATS)                                   :: sigma
      Type(MEF90_VECT)                                   :: UU0
      PetscReal                                          :: stiffness,temperature
      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr
      PetscReal                                          :: k

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDisplacement%BF,2)
      k = matprop%CoefficientLinSoft
      residualLoc = 0.0_Kr
      Do iGauss = 1,numGauss
         stiffness = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            stiffness = stiffness + damageDof(iDoF1) * elemDamage%BF(iDoF1,iGauss)
         End Do
         stiffness = (1.0_Kr - stiffness)**2 /( 1.0_Kr + ( k - 1.0_Kr )*(1.0_Kr - (1.0_Kr - stiffness)**2 ) )+ matProp%residualStiffness

         sigma = 0.0_Kr
         Do iDoF1 = 1,numDofDisplacement
            sigma = sigma + xDof(iDoF1) * elemDisplacement%GradS_BF(iDoF1,iGauss)
         End Do

         temperature = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperature = temperature + temperatureDoF(iDoF1) * elemDamage%BF(iDoF1,iGauss)
            End Do
         End If
         sigma = stiffness * (matProp%HookesLaw * (sigma - temperature * matProp%LinearThermalExpansion - plasticStrainCell))

         Do iDoF2 = 1,numDofDisplacement
            residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement%Gauss_C(iGauss) * (sigma .DotP. elemDisplacement%GradS_BF(iDoF2,iGauss))
         End Do

         UU0 = 0.0_Kr
         If (Associated(boundaryDisplacementDof) .AND. (matprop%cohesiveStiffness /= 0.0_Kr)) Then
            Do iDoF1 = 1,numDofDisplacement
               UU0 = UU0 + (xDof(iDoF1) - boundaryDisplacementDof(iDoF1)) * elemDisplacement%BF(iDoF1,iGauss)
            End Do
            UU0 = UU0 * matprop%cohesiveStiffness

            Do iDoF2 = 1,numDofDisplacement
               residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement%Gauss_C(iGauss) * (UU0 .DotP. elemDisplacement%BF(iDoF2,iGauss))
            End Do
         End If
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDisplacementLinSoftLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDisplacementATUnilateralHDLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDisplacementATUnilateralHDLoc:
!!!  
!!!  (c) 2016-2018 Blaise Bourdin bourdin@lsu.edu, Erwan TANNE erwan.tanne@gmail.com 
!!!

   Subroutine MEF90DefMechOperatorDisplacementATUnilateralHDLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      PetscReal,Intent(IN)                               :: CrackPressureCell

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: stiffnessD
      Type(MEF90_MATS)                                   :: inelasticStrain,sigmaD
      PetscReal                                          :: sigmaH
      Type(MEF90_VECT)                                   :: UU0
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr
      PetscScalar                                        :: kappa

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDisplacement%BF,2)
      kappa = ((matProp%HookesLaw * MEF90_MATS_IDENTITY) .dotP. MEF90_MATS_IDENTITY) / MEF90_DIM**2
      !!! This assumes that the Hooke's law is such that it preserves hydrostatic / deviatoric orthogonality

      residualLoc = 0.0_Kr
      Do iGauss = 1,numGauss
         inelasticStrain = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               inelasticStrain = inelasticStrain - matProp%LinearThermalExpansion * temperatureDof(iDoF1) * elemDamage%BF(iDoF1,iGauss)
            End Do         
         End If
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrain = inelasticStrain + xDof(iDoF1) * elemDisplacement%GradS_BF(iDoF1,iGauss)
         End Do
         inelasticStrain = inelasticStrain - plasticStrainCell

         stiffnessD = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            stiffnessD = stiffnessD + damageDof(iDoF1) * elemDamage%BF(iDoF1,iGauss)
         End Do
         stiffnessD = (1.0_Kr - stiffnessD)**2 + matProp%residualStiffness

         sigmaD = stiffnessD * (matProp%HookesLaw * deviatoricPart(inelasticStrain))
         sigmaH = kappa * (stiffnessD * MEF90HydrostaticDeviatoricDPenaltyFunction(Trace(inelasticStrain)) - & 
                                        MEF90HydrostaticDeviatoricDPenaltyFunction(-Trace(inelasticStrain))) * 0.5_Kr
         Do iDoF2 = 1,numDofDisplacement
            residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement%Gauss_C(iGauss) * ((sigmaD .DotP. elemDisplacement%GradS_BF(iDoF2,iGauss)) + &
                                                                                            sigmaH * Trace(elemDisplacement%GradS_BF(iDoF2,iGauss)))
         End Do ! iDof2

         !!! Contribution of the cohesive energy
         If (Associated(boundaryDisplacementDof) .AND. (matprop%cohesiveStiffness /= 0.0_Kr)) Then
            UU0 = 0.0_Kr
            Do iDoF1 = 1,numDofDisplacement
               UU0 = UU0 + (xDof(iDoF1) - boundaryDisplacementDof(iDoF1)) * elemDisplacement%BF(iDoF1,iGauss)
            End Do
            Do iDoF2 = 1,numDofDisplacement
               residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement%Gauss_C(iGauss) * matprop%cohesiveStiffness * (UU0 .DotP. elemDisplacement%BF(iDoF2,iGauss))
            End Do ! iDof2
         End If            
      End Do ! iGauss
      !flops = numGauss * ( 4. * numDofDisplacement**2 + 2. * numDofDamage + 3.)
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDisplacementATUnilateralHDLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDisplacementATUnilateralPHDLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDisplacementATUnilateralPHDLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDisplacementATUnilateralPHDLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      PetscReal,Intent(IN)                               :: CrackPressureCell

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: temperature,stiffness
      Type(MEF90_MATS)                                   :: plasticStrain,inelasticStrain,sigma
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDisplacement%BF,2)

      residualLoc = 0.0_Kr
      Do iGauss = 1,numGauss
         temperature = 0.0_Kr
         inelasticStrain = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperature = temperature - temperatureDof(iDoF1) * elemDamage%BF(iDoF1,iGauss)
            End Do         
            inelasticStrain = matProp%LinearThermalExpansion * temperature - plasticStrainCell
         End If
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrain = inelasticStrain + xDof(iDoF1) * elemDisplacement%GradS_BF(iDoF1,iGauss)
         End Do

         If (trace(inelasticStrain) < 0.0_Kr) Then
            stiffness = 1.0_Kr
         Else
            stiffness = 0.0_Kr
            Do iDoF1 = 1,numDofDamage
               stiffness = stiffness + damageDof(iDoF1) * elemDamage%BF(iDoF1,iGauss)
            End Do
            stiffness = (1.0_Kr - stiffness)**2 + matProp%residualStiffness
         End If
         
         sigma = stiffness * inelasticStrain
         sigma = matProp%HookesLaw * sigma
         Do iDoF2 = 1,numDofDisplacement
            residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement%Gauss_C(iGauss) * (sigma .DotP. elemDisplacement%GradS_BF(iDoF2,iGauss))
         End Do
      End Do
      !flops = numGauss * ( 2. * numDofDisplacement**2 + 3. * numDofDamage + 2.)
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDisplacementATUnilateralPHDLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDisplacementATUnilateralDeviatoricLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDisplacementATUnilateralDeviatoricLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDisplacementATUnilateralDeviatoricLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      PetscReal,Intent(IN)                               :: CrackPressureCell

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: temperature,stiffnessD
      Type(MEF90_MATS)                                   :: plasticStrain,inelasticStrain,sigma
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDisplacement%BF,2)

      residualLoc = 0.0_Kr
      Do iGauss = 1,numGauss
         temperature = 0.0_Kr
         inelasticStrain = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperature = temperature - temperatureDof(iDoF1) * elemDamage%BF(iDoF1,iGauss)
            End Do         
            inelasticStrain = matProp%LinearThermalExpansion * temperature
         End If
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrain = inelasticStrain + xDof(iDoF1) * elemDisplacement%GradS_BF(iDoF1,iGauss) - plasticStrainCell
         End Do

         stiffnessD = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            stiffnessD = stiffnessD + damageDof(iDoF1) * elemDamage%BF(iDoF1,iGauss)
         End Do
         stiffnessD = (1.0_Kr - stiffnessD)**2 + matProp%residualStiffness

         sigma = hydrostaticPart(inelasticStrain) + stiffnessD * deviatoricPart(inelasticStrain)
         sigma = matProp%HookesLaw * sigma
         Do iDoF2 = 1,numDofDisplacement
            residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement%Gauss_C(iGauss) * (sigma .DotP. elemDisplacement%GradS_BF(iDoF2,iGauss))
         End Do
      End Do
      !flops = numGauss * ( 2. * numDofDisplacement**2 + 3. * numDofDamage + 2.)
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDisplacementATUnilateralDeviatoricLoc


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDisplacementATUnilateralPSLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDisplacementATUnilateralPSLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDisplacementATUnilateralPSLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      PetscReal,Intent(IN)                               :: CrackPressureCell

      PetscInt                                           :: i,iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: temperature,stiffness
      Type(MEF90_MATS)                                   :: inelasticStrain,sigma
      PetscReal, Dimension(MEF90_DIM)                    :: PpalStrain
      Type(MEF90_MATS),Dimension(MEF90_DIM)              :: PpalDirection
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDisplacement%BF,2)

      residualLoc = 0.0_Kr
      Do iGauss = 1,numGauss
         temperature = 0.0_Kr
         inelasticStrain = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperature = temperature - temperatureDof(iDoF1) * elemDamage%BF(iDoF1,iGauss)
            End Do         
            inelasticStrain = matProp%LinearThermalExpansion * temperature
         End If
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrain = inelasticStrain + xDof(iDoF1) * elemDisplacement%GradS_BF(iDoF1,iGauss)
         End Do
         inelasticStrain = inelasticStrain - plasticStrainCell
         
         stiffness = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            stiffness = stiffness + damageDof(iDoF1) * elemDamage%BF(iDoF1,iGauss)
         End Do
         stiffness = (1.0_Kr - stiffness + matProp%residualStiffness)**2
   
         Call SpectralDecomposition(inelasticStrain,PpalStrain,PpalDirection)
         sigma = 0.0_Kr
         Do i = 1,MEF90_DIM 
            If (PpalStrain(i) < 0.0_Kr) Then
               sigma = sigma + (PpalStrain(i) * PpalDirection(i))
            Else
               sigma = sigma + Stiffness * (PpalStrain(i) * PpalDirection(i))
            End If
         End Do
         sigma = matProp%HookesLaw * sigma
         Do iDoF2 = 1,numDofDisplacement
            residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement%Gauss_C(iGauss) * (sigma .DotP. elemDisplacement%GradS_BF(iDoF2,iGauss))
         End Do
      End Do
      !flops = numGauss * ( 2. * numDofDisplacement**2 + 3. * numDofDamage + 2.)
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDisplacementATUnilateralPSLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDisplacementATUnilateralMasonryLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDisplacementATUnilateralMasonryLoc:
!!!  
!!!  (c) 2016 Blaise Bourdin bourdin@lsu.edu, Erwan TANNE erwan.tanne@gmail.com 
!!!

   Subroutine MEF90DefMechOperatorDisplacementATUnilateralMasonryLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      PetscReal,Intent(IN)                               :: CrackPressureCell

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: damageElem
      Type(MEF90_MATS)                                   :: inelasticStrain,positiveInelasticStrain,negativeInelasticStrain,stress
      Type(MEF90_MAT)                                    :: Pinv
      Type(MEF90_VECT)                                   :: UU0
      PetscReal                                          :: C1
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      If (matprop%HookesLaw%type /= MEF90HookesLawTypeIsotropic) Then
         SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unilateral Masonry not implemented for non isotropic Hooke laws: "//__FUNCT__,ierr)
      End If

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDisplacement%BF,2)
      C1 = matprop%HookesLaw%lambda / (matprop%HookesLaw%lambda + 2.0_Kr * matprop%HookesLaw%mu)

      residualLoc = 0.0_Kr
      Do iGauss = 1,numGauss
         damageElem = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            damageElem = damageElem + damageDof(iDoF1) * elemDamage%BF(iDoF1,iGauss)
         End Do

         inelasticStrain = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               inelasticStrain = inelasticStrain - temperatureDof(iDoF1)  * elemDamage%BF(iDoF1,iGauss) * matProp%LinearThermalExpansion
            End Do         
         End If
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrain = inelasticStrain + xDof(iDoF1) * elemDisplacement%GradS_BF(iDoF1,iGauss)
         End Do
         inelasticStrain = inelasticStrain - plasticStrainCell
         Call MasonryProjection(inelasticStrain,matprop%HookesLaw,positiveInelasticStrain,negativeInelasticStrain)
         inelasticStrain = (( (1.0_Kr - damageElem)**2 + matprop%residualStiffness) * positiveInelasticStrain) + &
                         & ((1.0_Kr + matprop%residualStiffness) * negativeInelasticStrain)
         stress = matprop%HookesLaw * inelasticStrain

         Do iDoF2 = 1,numDofDisplacement
            residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement%Gauss_C(iGauss) * (stress .DotP. elemDisplacement%GradS_BF(iDoF2,iGauss))
         End Do ! iDof2

         !!! Contribution of the cohesive energy
         If (Associated(boundaryDisplacementDof) .AND. (matprop%cohesiveStiffness /= 0.0_Kr)) Then
            UU0 = 0.0_Kr
            Do iDoF1 = 1,numDofDisplacement
               UU0 = UU0 + (xDof(iDoF1) - boundaryDisplacementDof(iDoF1)) * elemDisplacement%BF(iDoF1,iGauss)
            End Do
            Do iDoF2 = 1,numDofDisplacement
               residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement%Gauss_C(iGauss) * matprop%cohesiveStiffness * (UU0 .DotP. elemDisplacement%BF(iDoF2,iGauss))
            End Do ! iDof2
         End If            
      End Do ! iGauss
      !flops = numGauss * ( 4. * numDofDisplacement**2 + 2. * numDofDamage + 3.)
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDisplacementATUnilateralMasonryLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechRHSDisplacementLoc"
!!!
!!!  
!!!  MEF90DefMechRHSDisplacementLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechRHSDisplacementLoc(residualLoc,xDof,forceCell,pressureForceCell,matprop,elemDisplacement,damageDof,elemDamage,CrackPressureCell)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,damageDof
      Type(MEF90_VECT),Intent(IN)                        :: forceCell
      PetscReal,Intent(IN)                               :: pressureForceCell
      PetscReal,Intent(IN)                               :: CrackPressureCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      
      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numGauss,numDofDamage
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr
      Type(MEF90_VECT)                                   :: totalForce,GradientDamageElem


      totalForce = forceCell + pressureForceCell * elemDisplacement%outerNormal
      If (Norm(totalForce) > 0.0_Kr .or. CrackPressureCell/=0 ) Then
         numDofDisplacement = size(elemDisplacement%BF,1)
         numGauss = size(elemDisplacement%BF,2)
         numDofDamage = size(elemDamage%BF,1)
         Do iGauss = 1,numGauss
            GradientDamageElem=0.0_Kr
            Do iDoF1 = 1,numDofDamage
               GradientDamageElem = GradientDamageElem + damageDof(iDoF1) * elemDamage%Grad_BF(iDoF1,iGauss) 
            End Do
            Do iDoF2 = 1,numDofDisplacement
               residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement%Gauss_C(iGauss) * &
                                  ( (totalForce + GradientDamageElem*CrackPressureCell ) .DotP. elemDisplacement%BF(iDoF2,iGauss))
            End Do
         End Do
         !flops = 2 * numGauss * numDofDisplacement**2
         !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      End If
   End Subroutine MEF90DefMechRHSDisplacementLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDisplacementDampingLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDisplacementDampingLoc: Contribution of the damping term \int delta (u-u_{n-1}) 
!!!                                              in an element
!!!  
!!!  (c) 2018 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDisplacementDampingLoc(residualLoc,xDof,xpreviousStepDof,delta,elemDisplacement)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,xpreviousStepDof
      PetscReal                                          :: delta
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      
      Type(MEF90_VECT)                                   :: contrElem
      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numGauss
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numGauss = size(elemDisplacement%BF,2)
      Do iGauss = 1,numGauss
         contrElem = 0.0_Kr
         Do iDoF1 = 1,numDofDisplacement
            contrElem = contrElem + (xDof(iDof1) - xpreviousStepDof(iDoF1)) * elemDisplacement%BF(iDoF1,iGauss) 
         End Do
         Do iDoF2 = 1,numDofDisplacement
            residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement%Gauss_C(iGauss) * &
                                 delta * ( contrElem .DotP. elemDisplacement%BF(iDoF2,iGauss) )
         End Do
      End Do
      !All flops logged
   End Subroutine MEF90DefMechOperatorDisplacementDampingLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDisplacement"
!!!
!!!  
!!!  MEF90DefMechOperatorDisplacement: Build the operator. When called in SNES, the solution time should always match the target time, 
!!!                                    so there is no need for interpolation of the forcees, external, and boundary values
!!!  
!!!  (c) 2012-18 Blaise Bourdin bourdin@lsu.edu, Erwan Tanne erwan.tanne@gmail.com
!!!

   Subroutine MEF90DefMechOperatorDisplacement(snesDisplacement,displacement,residual,MEF90DefMechCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesDisplacement
      Type(Vec),Intent(IN)                               :: displacement
      Type(Vec),Intent(INOUT)                            :: residual
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(DM)                                           :: mesh
      Type(SectionReal)                                  :: displacementSec,displacementPreviousStepSec,residualSec
      Type(SectionReal)                                  :: plasticStrainSec,temperatureSec,damageSec
      Type(SectionReal)                                  :: cumulatedDissipatedPlasticEnergySec
      Type(SectionReal)                                  :: boundaryDisplacementSec,forceSec,pressureForceSec,CrackPressureSec
      PetscReal,Dimension(:),Pointer                     :: displacementPtr,residualPtr
      PetscReal,Dimension(:),Pointer                     :: boundaryDisplacementPtr
      PetscReal,Dimension(:),Pointer                     :: boundaryDisplacementDof,displacementDof,displacementPreviousStepDof,damageDof,temperatureDof,residualLoc
      PetscReal,Dimension(:),Pointer                     :: pressureForceLoc,CrackPressureLoc
      PetscReal                                          :: pressureForceCell,CrackPressureCell
      PetscReal,Dimension(:),Pointer                     :: forceLoc
      Type(MEF90_VECT)                                   :: forceCell
      PetscReal,Dimension(:),Pointer                     :: plasticStrainLoc
      Type(MEF90_MATS)                                   :: plasticStrainCell
      PetscReal,Dimension(:),Pointer                     :: cumulatedDissipatedPlasticEnergyLoc
      PetscReal                                          :: cumulatedDissipatedPlasticEnergyCell
      Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof,bcIS
      PetscInt,dimension(:),Pointer                      :: setID,cellID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder,cell
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90DefMechVertexSetOptions_Type),Pointer    :: vertexSetOptions
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemDamage
      Type(MEF90Element_Type)                            :: elemDisplacementType,elemDamageType
      PetscInt                                           :: p,dof,c,dim,numDof
      PetscInt                                           :: nVal
      PetscReal,Dimension(:),Pointer                     :: nullPtr
      
      Procedure(MEF90DefMechOperatorLoc),pointer         :: localOperatorFunction
      Procedure(MEF90DefMechRHSLoc),pointer              :: localRHSFunction
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: globalOptions
      Type(Vec)                                          :: damageOld
      PetscReal                                          :: damageMin,damageMax,damageMaxChange
      Character(len=MEF90_MXSTRLEN)                      :: IOBuffer
    
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)
      If ( (GlobalOptions%solverType == MEF90DefMech_SolverTypeQuasiNewton1) .OR. (GlobalOptions%solverType == MEF90DefMech_SolverTypeQuasiNewton2)) Then
         Call VecDuplicate(MEF90DefMechCtx%damage,damageOld,ierr);CHKERRQ(ierr)
         Call VecCopy(MEF90DefMechCtx%damage,damageOld,ierr);CHKERRQ(ierr)
         Call SNESSolve(MEF90DefMechCtx%SNESdamage,PETSC_NULL_OBJECT,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)
         Call VecMin(MEF90DefMechCtx%damage,PETSC_NULL_INTEGER,damageMin,ierr);CHKERRQ(ierr)
         Call VecMax(MEF90DefMechCtx%damage,PETSC_NULL_INTEGER,damageMax,ierr);CHKERRQ(ierr)
         Call VecAxPy(damageOld,-1.0_Kr,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)
         Call VecNorm(damageOld,NORM_INFINITY,damageMaxChange,ierr);CHKERRQ(ierr)
         Call VecDestroy(damageOld,ierr);CHKERRQ(ierr)
       End If

      localOperatorFunction =>MEF90DefMechOperatorNull      
      Call SNESGetDM(snesDisplacement,mesh,ierr);CHKERRQ(ierr)

      !!! Get Section and scatter associated with each field
      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,displacementSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(displacementSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,displacement,ierr);CHKERRQ(ierr) 
      If (globalOptions%dampingCoefficientDisplacement * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,displacementpreviousStepSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(displacementpreviousStepSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,MEF90DefMechCtx%displacementpreviousStep,ierr);CHKERRQ(ierr) 
      EndIf

      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,boundaryDisplacementSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(boundaryDisplacementSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,MEF90DefMechCtx%boundaryDisplacement,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,residualSec,ierr);CHKERRQ(ierr)


      Call SectionRealSet(residualSec,0.0_Kr,ierr);CHKERRQ(ierr)
      Call VecSet(residual,0.0_Kr,ierr);CHKERRQ(ierr)

      Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,damageSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(damageSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)   

      Nullify(TemperatureDof)
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,temperatureSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(temperatureSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr)   
      Else
         temperatureSec%v = 0
      End If

      Nullify(ForceLoc)
      If (Associated(MEF90DefMechCtx%force)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%cellDMVectSec,forceSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(forceSec,MEF90DefMechCtx%cellDMVectScatter,SCATTER_REVERSE,MEF90DefMechCtx%force,ierr);CHKERRQ(ierr)
      Else
         ForceSec%v = 0
      End If
      
      Nullify(pressureForceLoc)
      If (Associated(MEF90DefMechCtx%pressureForce)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%cellDMScalSec,pressureForceSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(pressureForceSec,MEF90DefMechCtx%CellDMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%pressureForce,ierr);CHKERRQ(ierr)
      Else
         pressureForceSec%v = 0
      End If

      Nullify(CrackPressureLoc)
      If (Associated(MEF90DefMechCtx%CrackPressure)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%cellDMScalSec,CrackPressureSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(CrackPressureSec,MEF90DefMechCtx%CellDMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%CrackPressure,ierr);CHKERRQ(ierr)
      Else
         CrackPressureSec%v = 0
      End If


      Nullify(plasticStrainLoc)
      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%CellDMMatSSec,plasticStrainSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticStrainSec,MEF90DefMechCtx%CellDMMatSScatter,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)
         
         Call SectionRealDuplicate(MEF90DefMechCtx%CellDMScalSec,cumulatedDissipatedPlasticEnergySec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(cumulatedDissipatedPlasticEnergySec,MEF90DefMechCtx%CellDMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%cumulatedDissipatedPlasticEnergy,ierr);CHKERRQ(ierr)
      Else
         PlasticStrainSec%v = 0
         cumulatedDissipatedPlasticEnergySec%v = 0
      End If

      !!! get IS for cell sets
      Call DMMeshGetDimension(mesh,dim,ierr);CHKERRQ(ierr)
      Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

      !!!
      !!! We loop over all element twice. The first time in order to assembly all non BC cell sets
      !!! In the second pass, we only update the BC where necessary
      !!! vertex set BC are updated last, so that they override cell set BC
      !!!
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
         elemDamageType       = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)

         Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
         If (Size(cellID) > 0) Then
            Select Case (cellSetOptions%damageType)
            Case (MEF90DefMech_damageTypeAT1Elastic,MEF90DefMech_damageTypeAT2Elastic,MEF90DefMech_damageTypeKKLElastic)
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  QuadratureOrder = elemDisplacementType%order - 1 + max(elemDisplacementType%order - 1, elemDamageType%order)
               Else
                  QuadratureOrder = 2 * (elemDisplacementType%order - 1)
               End If
               localOperatorFunction => MEF90DefMechOperatorDisplacementElasticLoc
               localRHSFunction => MEF90DefMechRHSDisplacementLoc

            Case (MEF90DefMech_damageTypeLinSoft)
               Select Case(cellSetOptions%unilateralContactType)
               Case (MEF90DefMech_unilateralContactTypeNone)
                  QuadratureOrder = 5 !* (elemDisplacementType%order - 1) + 5 * elemDamageType%order
                  localOperatorFunction => MEF90DefMechOperatorDisplacementLinSoftLoc
               Case default
                  localOperatorFunction => MEF90DefMechOperatorNull
               End Select
               localRHSFunction => MEF90DefMechRHSDisplacementLoc

            Case (MEF90DefMech_damageTypeAT1,MEF90DefMech_damageTypeAT2)
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  QuadratureOrder = elemDisplacementType%order - 1 + max(elemDisplacementType%order - 1, elemDamageType%order) + 2 * elemDamageType%order
               Else
                  QuadratureOrder = 2 * (elemDisplacementType%order - 1) + 2 * elemDamageType%order
               End If
               Select Case(cellSetOptions%unilateralContactType)
               Case (MEF90DefMech_unilateralContactTypeNone,MEF90DefMech_unilateralContactTypeBrittleDuctile,MEF90DefMech_unilateralContactTypeHybridHydrostaticDeviatoric)
                  localOperatorFunction => MEF90DefMechOperatorDisplacementATLoc
               Case (MEF90DefMech_unilateralContactTypeHydrostaticDeviatoric)
                  localOperatorFunction => MEF90DefMechOperatorDisplacementATUnilateralHDLoc
               Case (MEF90DefMech_unilateralContactTypePositiveHydrostatic)
                  localOperatorFunction => MEF90DefMechOperatorDisplacementATUnilateralPHDLoc
               Case (MEF90DefMech_unilateralContactTypeDeviatoric)
                  localOperatorFunction => MEF90DefMechOperatorDisplacementATUnilateralDeviatoricLoc
               Case (MEF90DefMech_unilateralContactTypePrincipalStrains)
                  localOperatorFunction => MEF90DefMechOperatorDisplacementATUnilateralPSLoc
               Case (MEF90DefMech_unilateralContactTypeMasonry)
                  localOperatorFunction => MEF90DefMechOperatorDisplacementATUnilateralMasonryLoc
               End Select
               localRHSFunction => MEF90DefMechRHSDisplacementLoc

            Case (MEF90DefMech_damageTypeKKL)
               QuadratureOrder = 2 * (elemDisplacementType%order - 1) + 3
               Select Case(cellSetOptions%unilateralContactType)
               Case (MEF90DefMech_unilateralContactTypeNone)
                  localOperatorFunction => MEF90DefMechOperatorDisplacementKKLLoc
               End Select
               localRHSFunction => MEF90DefMechRHSDisplacementLoc
            End Select

            !!! Allocate storage for fields at dof and Gauss points
            Allocate(residualLoc(ElemDisplacementType%numDof))
            Allocate(displacementDof(ElemDisplacementType%numDof))
            If (globalOptions%dampingCoefficientDisplacement * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
               Allocate(displacementPreviousStepDof(ElemDisplacementType%numDof))
            Else
               Nullify(displacementpreviousStepDof)
            End If
            If (matpropSet%cohesiveStiffness /= 0.0_Kr) Then
               Allocate(boundaryDisplacementDof(ElemDisplacementType%numDof))
            Else
               Nullify(boundaryDisplacementDof)
            End If
            Allocate(damageDof(ElemDamageType%numDof))
            If (Associated(MEF90DefMechCtx%temperature)) Then
               Allocate(temperatureDof(ElemDamageType%numDof))
            Else
               Nullify(temperatureDof)
            End If
            !!! Allocate elements 
            Call MEF90Element_Create(mesh,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
            Call MEF90Element_Create(mesh,setIS,elemDamage,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
            Do cell = 1,size(cellID)
               !!! Get value of each field at each Dof of the local element
               Call SectionRealRestrictClosure(displacementSec,MEF90DefMechCtx%DMVect,cellID(cell),elemDisplacementType%numDof,displacementDof,ierr);CHKERRQ(ierr)
               Call SectionRealRestrictClosure(damageSec,MEF90DefMechCtx%DMScal,cellID(cell),elemDamageType%numDof,damageDof,ierr);CHKERRQ(ierr)
               If (matpropSet%cohesiveStiffness /= 0.0_Kr) Then
                  Call SectionRealRestrictClosure(boundaryDisplacementSec,MEF90DefMechCtx%DMVect,cellID(cell),elemDisplacementType%numDof,boundaryDisplacementDof,ierr);CHKERRQ(ierr)
               End If
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  Call SectionRealRestrictClosure(temperatureSec,MEF90DefMechCtx%DMScal,cellID(cell),elemDamageType%numDof,temperatureDof,ierr);CHKERRQ(ierr)
               End If
               If (Associated(MEF90DefMechCtx%force)) Then
                  Call SectionRealRestrict(forceSec,cellID(cell),forceLoc,ierr);CHKERRQ(ierr)
                  forceCell = -forceLoc
               Else
                  forceCell = 0.0_Kr
               End If
               If (Associated(MEF90DefMechCtx%pressureForce)) Then
                  Call SectionRealRestrict(pressureForceSec,cellID(cell),pressureForceLoc,ierr);CHKERRQ(ierr)
                  pressureForceCell = -pressureForceLoc(1)
               Else
                  pressureForceCell = 0.0_Kr
               End If

               If (Associated(MEF90DefMechCtx%CrackPressure)) Then
                  Call SectionRealRestrict(CrackPressureSec,cellID(cell),CrackPressureLoc,ierr);CHKERRQ(ierr)
                  CrackPressureCell = CrackPressureLoc(1)
               Else
                  CrackPressureCell = 0.0_Kr
               End If

               If (Associated(MEF90DefMechCtx%plasticStrain)) Then
                  Call SectionRealRestrict(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
                  Call SectionRealRestrict(cumulatedDissipatedPlasticEnergySec,cellID(cell),cumulatedDissipatedPlasticEnergyLoc,ierr);CHKERRQ(ierr)
                  plasticStrainCell = plasticStrainLoc
                  cumulatedDissipatedPlasticEnergyCell = cumulatedDissipatedPlasticEnergyLoc(1)
               Else
                  plasticStrainCell = 0.0_Kr
                  cumulatedDissipatedPlasticEnergyCell = 0.0_Kr
               End If
                  
               residualLoc = 0.0_Kr
               If (globalOptions%dampingCoefficientDisplacement * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
                  Call SectionRealRestrictClosure(displacementpreviousStepSec,MEF90DefMechCtx%DMVect,cellID(cell),elemDisplacementType%numDof,displacementpreviousStepDof,ierr);CHKERRQ(ierr)
                  call MEF90DefMechOperatorDisplacementDampingLoc(residualLoc,displacementDof,displacementpreviousStepDof,globalOptions%dampingCoefficientDisplacement / MEF90DefMechCtx%timeStep,elemDisplacement(cell))
               End If
               Call localOperatorFunction   (residualLoc,displacementDof,nullPtr,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matpropSet,elemDisplacement(cell),elemDamage(cell),CrackPressureCell)
               Call localRHSFunction        (residualLoc,displacementDof,forceCell,pressureForceCell,matpropSet,elemDisplacement(cell),damageDof,elemDamage(cell),CrackPressureCell)
               Call SectionRealUpdateClosure(residualSec,MEF90DefMechCtx%DMVect,cellID(cell),residualLoc,ADD_VALUES,ierr);CHKERRQ(iErr)

               If (Associated(MEF90DefMechCtx%force)) Then
                  Call SectionRealRestore(forceSec,cellID(cell),forceLoc,ierr);CHKERRQ(ierr)
               End If
               If (Associated(MEF90DefMechCtx%pressureForce)) Then
                  Call SectionRealRestore(pressureForceSec,cellID(cell),pressureForceLoc,ierr);CHKERRQ(ierr)
               End If
               If (Associated(MEF90DefMechCtx%CrackPressure)) Then
                  Call SectionRealRestore(CrackPressureSec,cellID(cell),CrackPressureLoc,ierr);CHKERRQ(ierr)
               End If
               If (Associated(MEF90DefMechCtx%plasticStrain)) Then
                  Call SectionRealRestore(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
                  Call SectionRealRestore(cumulatedDissipatedPlasticEnergySec,cellID(cell),cumulatedDissipatedPlasticEnergyLoc,ierr);CHKERRQ(ierr)
               End If
            End Do

            DeAllocate(displacementDof)
            If (globalOptions%dampingCoefficientDisplacement * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
               DeAllocate(displacementPreviousStepDof)
            End If
            If (matpropSet%cohesiveStiffness /= 0.0_Kr) Then
               DeAllocate(boundaryDisplacementDof)
            End If
            DeAllocate(damageDof)
            If (Associated(MEF90DefMechCtx%temperature)) Then
               DeAllocate(temperatureDof)
            End If
            DeAllocate(residualLoc)

            Call MEF90Element_Destroy(elemDisplacement,ierr)
            Call MEF90Element_Destroy(elemDamage,ierr)
         End If 
         Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do ! set

      !!! "Ghost update" for the residual Section
      Call SectionRealComplete(residualSec,ierr);CHKERRQ(ierr)
      !!! Scatter back from SectionReal to Vec
      Call SectionRealToVec(residualSec,MEF90DefMechCtx%DMVectScatter,SCATTER_FORWARD,residual,ierr);CHKERRQ(ierr)

      !!!
      !!! Cell set BC
      !!!
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Do c = 1, dim
            If (cellSetOptions%Has_DisplacementBC(c)) Then
               Call DMMeshGetStratumIS(MEF90DefMechCtx%DMVect,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
               Call MEF90ISCreateCelltoVertex(MEF90DefMechCtx%DMVect,PETSC_COMM_WORLD,setIS,bcIS,ierr)
               Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
               Allocate(displacementPtr(nval),stat=ierr)
               Allocate(residualPtr(nval),stat=ierr)
               Allocate(boundaryDisplacementPtr(nval),stat=ierr)
               Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMVect,displacement,displacementPtr,bcIS,c,ierr)
               Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMVect,MEF90DefMechCtx%boundaryDisplacement,boundaryDisplacementPtr,bcIS,c,ierr)
               residualPtr = displacementPtr - boundaryDisplacementPtr
               Call MEF90VecSetValuesISdof(MEF90DefMechCtx%DMVect,residual,residualPtr,bcIS,c,INSERT_VALUES,ierr)
               DeAllocate(boundaryDisplacementPtr)
               DeAllocate(residualPtr)
               DeAllocate(displacementPtr)
               Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
               Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
            End If ! cellSetOptions%Has_BC
         End Do
      End Do ! set
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)

      !!!
      !!! Vertex set BC
      !!!
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxVertexSetOptions(MEF90DefMechCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
         Do c = 1, dim
            If (vertexSetOptions%Has_DisplacementBC(c)) Then
               Call DMMeshGetStratumIS(MEF90DefMechCtx%DMVect,'Vertex Sets',setID(set),bcIS,ierr);CHKERRQ(iErr)
               Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
               Allocate(displacementPtr(nval),stat=ierr)
               Allocate(residualPtr(nval),stat=ierr)
               Allocate(boundaryDisplacementPtr(nval),stat=ierr)
               Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMVect,displacement,displacementPtr,bcIS,c,ierr)
               Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMVect,MEF90DefMechCtx%boundaryDisplacement,boundaryDisplacementPtr,bcIS,c,ierr)
               residualPtr = displacementPtr - boundaryDisplacementPtr
               Call MEF90VecSetValuesISdof(MEF90DefMechCtx%DMVect,residual,residualPtr,bcIS,c,INSERT_VALUES,ierr)
               DeAllocate(boundaryDisplacementPtr)
               DeAllocate(residualPtr)
               DeAllocate(displacementPtr)
               Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
            End If ! vertexSetOptions%Has_BC
         End Do
      End Do ! set
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
   
      Call VecAssemblyBegin(residual,ierr)
      Call VecAssemblyEnd(residual,ierr)

      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
         Call SectionRealDestroy(cumulatedDissipatedPlasticEnergySec,ierr);CHKERRQ(ierr)
      End If
      
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
      End If

      If (Associated(MEF90DefMechCtx%pressureForce)) Then
         Call SectionRealDestroy(pressureForceSec,ierr);CHKERRQ(ierr)
      End If

      If (Associated(MEF90DefMechCtx%CrackPressure)) Then
         Call SectionRealDestroy(CrackPressureSec,ierr);CHKERRQ(ierr)
      End If

      If (Associated(MEF90DefMechCtx%force)) Then
         Call SectionRealDestroy(forceSec,ierr);CHKERRQ(ierr)
      End If

      Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(boundaryDisplacementSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(residualSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(displacementSec,ierr);CHKERRQ(ierr)
      If (globalOptions%dampingCoefficientDisplacement * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
         Call SectionRealDestroy(displacementPreviousStepSec,ierr);CHKERRQ(ierr)
      End If
   End Subroutine MEF90DefMechOperatorDisplacement

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDisplacement"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDisplacement:
!!!  
!!!  (c) 2012-16 Blaise Bourdin bourdin@lsu.edu,Erwan Tanne erwan.tanne@gmail.com
!!!

   Subroutine MEF90DefMechBilinearFormDisplacement(snesDisplacement,displacement,A,M,flg,MEF90DefMechCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesDisplacement
      Type(Vec),Intent(IN)                               :: displacement
      Type(Mat),Intent(INOUT)                            :: A,M
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr  
         
      Type(SectionReal)                                  :: defaultSec
      Type(SectionReal)                                  :: displacementSec,damageSec,temperatureSec,plasticStrainSec
      Type(SectionReal)                                  :: cumulatedDissipatedPlasticEnergySec
      PetscReal,Dimension(:),Pointer                     :: displacementDOF,damageDof,temperatureDof,nullPtr
      PetscReal,Dimension(:),Pointer                     :: plasticStrainLoc
      Type(MEF90_MATS)                                   :: plasticStrainCell
      PetscReal,Dimension(:),Pointer                     :: cumulatedDissipatedPlasticEnergyLoc
      PetscReal                                          :: cumulatedDissipatedPlasticEnergyCell
      Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90DefMechVertexSetOptions_Type),Pointer    :: vertexSetOptions
      Type(DM)                                           :: mesh
      Type(MEF90Element_Type)                            :: ElemDisplacementType,ElemDamageType
      PetscInt,Dimension(:),Pointer                      :: setIdx,bcIdx,Cone
      Type(IS)                                           :: bcIS
      PetscInt                                           :: cell,v,numBC,numDoF,numCell,c,dim
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscInt                                           :: iGauss,iDoF1,iDoF2
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemDamage
      
      Procedure(MEF90DefMechBilinearFormLoc),pointer     :: localAssemblyFunction
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: globalOptions
      Type(Vec)                                          :: damageOld
      PetscReal                                          :: damageMin,damageMax,damageMaxChange,cohesiveStiffness
      Character(len=MEF90_MXSTRLEN)                      :: IOBuffer
    
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)
      If (GlobalOptions%solverType == MEF90DefMech_SolverTypeQuasiNewton2) Then
         Call VecDuplicate(MEF90DefMechCtx%damage,damageOld,ierr);CHKERRQ(ierr)
         Call VecCopy(MEF90DefMechCtx%damage,damageOld,ierr);CHKERRQ(ierr)
         Call SNESSolve(MEF90DefMechCtx%SNESdamage,PETSC_NULL_OBJECT,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)
         Call VecMin(MEF90DefMechCtx%damage,PETSC_NULL_INTEGER,damageMin,ierr);CHKERRQ(ierr)
         Call VecMax(MEF90DefMechCtx%damage,PETSC_NULL_INTEGER,damageMax,ierr);CHKERRQ(ierr)
         Call VecAxPy(damageOld,-1.0_Kr,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)
         Call VecNorm(damageOld,NORM_INFINITY,damageMaxChange,ierr);CHKERRQ(ierr)
         Call VecDestroy(damageOld,ierr);CHKERRQ(ierr)
       End If
      
      localAssemblyFunction =>MEF90DefMechBilinearFormNull      
      Call MatZeroEntries(A,ierr);CHKERRQ(ierr)
      Call SNESGetDM(snesDisplacement,mesh,ierr);CHKERRQ(ierr)
      
      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,displacementSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(displacementSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,displacement,ierr);CHKERRQ(ierr) 

      Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,damageSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(damageSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)   

      Nullify(TemperatureDof)
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,temperatureSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(temperatureSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr)   
      Else
         temperatureSec%v = 0
      End If

      Nullify(plasticStrainLoc)
      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%CellDMMatSSec,plasticStrainSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticStrainSec,MEF90DefMechCtx%CellDMMatSScatter,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)
         
         Call SectionRealDuplicate(MEF90DefMechCtx%CellDMScalSec,cumulatedDissipatedPlasticEnergySec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(cumulatedDissipatedPlasticEnergySec,MEF90DefMechCtx%CellDMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%cumulatedDissipatedPlasticEnergy,ierr);CHKERRQ(ierr)
      Else
         PlasticStrainSec%v = 0
         cumulatedDissipatedPlasticEnergySec%v = 0
      End If

      !!! get IS for cell sets
      Call DMMeshGetDimension(mesh,dim,ierr);CHKERRQ(ierr)
      Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID) 
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         ElemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
         ElemDamageType       = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)
         !!! Go through cells and call local assembly functions
         Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
         If ((Size(cellID) > 0) .AND. (ElemDisplacementType%coDim == 0)) Then
         !!! Call proper local assembly depending on the type of damage law

            !!! select the proper local assembly routine, compute proper integration order
            Select Case (cellSetOptions%damageType)
            Case (MEF90DefMech_damageTypeAT1Elastic,MEF90DefMech_damageTypeAT2Elastic,MEF90DefMech_damageTypeKKLElastic)
               QuadratureOrder = 2 * (elemDisplacementType%order - 1)
               localAssemblyFunction => MEF90DefMechBilinearFormDisplacementElasticLoc
            Case (MEF90DefMech_damageTypeLinSoft)
               QuadratureOrder = 5 !* (elemDisplacementType%order - 1) + 5 * ElemDamageType%order
               localAssemblyFunction => MEF90DefMechBilinearFormDisplacementLinSoftLoc
            Case (MEF90DefMech_damageTypeKKL)
               QuadratureOrder = 5 !* (elemDisplacementType%order - 1) + 5 * ElemDamageType%order
               Select Case(cellSetOptions%unilateralContactType)
               Case (MEF90DefMech_unilateralContactTypeNone)
                  localAssemblyFunction => MEF90DefMechBilinearFormDisplacementKKLLoc
               End Select
            Case (MEF90DefMech_damageTypeAT1,MEF90DefMech_damageTypeAT2)
               QuadratureOrder = 2 * (elemDisplacementType%order - 1) + 2 * ElemDamageType%order
               Select Case(cellSetOptions%unilateralContactType)
               Case (MEF90DefMech_unilateralContactTypeNone,MEF90DefMech_unilateralContactTypeBrittleDuctile,MEF90DefMech_unilateralContactTypeHybridHydrostaticDeviatoric)
                  localAssemblyFunction => MEF90DefMechBilinearFormDisplacementATLoc
               Case (MEF90DefMech_unilateralContactTypeHydrostaticDeviatoric)
                  localAssemblyFunction => MEF90DefMechBilinearFormDisplacementATUnilateralHDLoc
               Case (MEF90DefMech_unilateralContactTypePositiveHydrostatic)
                  localAssemblyFunction => MEF90DefMechBilinearFormDisplacementATUnilateralPHDLoc
               Case (MEF90DefMech_unilateralContactTypeDeviatoric)
                  localAssemblyFunction => MEF90DefMechBilinearFormDisplacementATUnilateralDeviatoricLoc
               Case (MEF90DefMech_unilateralContactTypeMasonry)
                  ! due to the lag technique used for the implementation of Masonry, the Bilinear form is that of an elasticity problem
                  localAssemblyFunction => MEF90DefMechBilinearFormDisplacementElasticLoc
               Case (MEF90DefMech_unilateralContactTypePrincipalStrains)
                  localAssemblyFunction => MEF90DefMechBilinearFormDisplacementATUnilateralPSLoc
               End Select
            End Select
            !!! Allocate storage for fields at dof and Gauss points
            !!! Leaving plasticstrain aside until I can remember how it is supposed to be dealt with
            Allocate(Aloc(ElemDisplacementType%numDof,ElemDisplacementType%numDof))
            Allocate(displacementDof(ElemDisplacementType%numDof))
            Allocate(damageDof(ElemDamageType%numDof))
            If (Associated(MEF90DefMechCtx%temperature)) Then
               Allocate(temperatureDof(ElemDamageType%numDof))
            End If
            
            !!! Allocate elements 
            Call MEF90Element_Create(mesh,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
            Call MEF90Element_Create(mesh,setIS,elemDamage,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
            Do cell = 1,size(cellID)
               !!! Get value of each field at each Dof of the local element
               Call SectionRealRestrictClosure(displacementSec,MEF90DefMechCtx%DMVect,cellID(cell),elemDisplacementType%numDof,displacementDof,ierr);CHKERRQ(ierr)
               Call SectionRealRestrictClosure(damageSec,MEF90DefMechCtx%DMScal,cellID(cell),elemDamageType%numDof,damageDof,ierr);CHKERRQ(ierr)
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  Call SectionRealRestrictClosure(temperatureSec,MEF90DefMechCtx%DMScal,cellID(cell),elemDamageType%numDof,temperatureDof,ierr);CHKERRQ(ierr)
               End If
               If (Associated(MEF90DefMechCtx%plasticStrain)) Then
                  Call SectionRealRestrict(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
                  Call SectionRealRestrict(cumulatedDissipatedPlasticEnergySec,cellID(cell),cumulatedDissipatedPlasticEnergyLoc,ierr);CHKERRQ(ierr)
                  plasticStrainCell = plasticStrainLoc
                  cumulatedDissipatedPlasticEnergyCell = cumulatedDissipatedPlasticEnergyLoc(1)
               Else
                  plasticStrainCell = 0.0_Kr
                  cumulatedDissipatedPlasticEnergyCell = 0.0_Kr
               End If
                  
               Aloc = 0.0_Kr
               If (globalOptions%dampingCoefficientDisplacement * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
                  Call MEF90_MassMatrixAssembleLoc(Aloc,globalOptions%dampingCoefficientDisplacement / MEF90DefMechCtx%timeStep,elemDisplacement(cell))
               End If
               Call localAssemblyFunction(Aloc,displacementDof,nullPtr,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matpropSet,elemDisplacement(cell),elemDamage(cell))
               Call DMmeshAssembleMatrix(A,mesh,displacementSec,cellID(cell),ALoc,ADD_VALUES,ierr);CHKERRQ(ierr)

               If (Associated(MEF90DefMechCtx%plasticStrain)) Then
                  Call SectionRealRestore(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
                  Call SectionRealRestore(cumulatedDissipatedPlasticEnergySec,cellID(cell),cumulatedDissipatedPlasticEnergyLoc,ierr);CHKERRQ(ierr)
               End If
            End Do

            DeAllocate(displacementDof)
            DeAllocate(damageDof)
            If (Associated(MEF90DefMechCtx%temperature)) Then
               DeAllocate(temperatureDof)
            End If
            Call MEF90Element_Destroy(elemDisplacement,ierr)
            Call MEF90Element_Destroy(elemDamage,ierr)
            DeAllocate(Aloc)    
         End If 

         Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do ! set
      Call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
      Call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
      !!!
      !!! Boundary conditions at cell sets
      !!!
      Call DMmeshGetStratumSize(mesh,"height",0,numCell,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call MEF90ISCreateCelltoVertex(mesh,PETSC_COMM_WORLD,setIS,bcIS,ierr)
         Do c = 1, dim
            If (cellSetOptions%Has_displacementBC(c)) Then
               Call DMMeshISCreateISglobaldof(mesh,bcIS,c-1,setISdof,ierr);CHKERRQ(ierr)
               Call MatZeroRowsColumnsIS(A,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
               Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
            End If ! cellSetOptions%Has_displacementBC
         End Do ! c
         Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do     
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
   
      !!!
      !!! Boundary conditions at vertex sets
      !!!
      Call DMmeshGetLabelIdIS(mesh,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxVertexSetOptions(MEF90DefMechCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
         Do c = 1, dim
            If (vertexSetOptions%Has_displacementBC(c)) Then
               Call DMMeshGetStratumIS(mesh,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
               Call DMMeshISCreateISglobaldof(mesh,setIS,c-1,setISdof,ierr);CHKERRQ(ierr)
               Call MatZeroRowsColumnsIS(A,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
               Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
               Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
            End If
         End Do
      End Do
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      
      Call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
      Call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)

      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
         Call SectionRealDestroy(cumulatedDissipatedPlasticEnergySec,ierr);CHKERRQ(ierr)
      End If
      
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
      End If

      Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(displacementSec,ierr);CHKERRQ(ierr)
      flg = SAME_NONZERO_PATTERN
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
      Type(Vec),Intent(IN)                               :: xVec
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,Dimension(:),Pointer                     :: work
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(SectionReal)                                  :: xSec
      Type(SectionReal)                                  :: forceSec,pressureForceSec
      Type(IS)                                           :: CellSetGlobalIS,setIS
      PetscInt,Dimension(:),Pointer                      :: setID
      PetscInt                                           :: set,QuadratureOrder,numCell
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90Element_Type)                            :: elemDisplacementType
      Type(MEF90DefMechGlobalOptions_Type),Pointer       :: globalOptions
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      PetscReal                                          :: myWork,myWorkSet
      
      
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%globalOptionsBag,globalOptions,ierr);CHKERRQ(ierr)

      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,xSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(xSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,xVec,ierr);CHKERRQ(ierr)

      Call SectionRealDuplicate(MEF90DefMechCtx%cellDMVectSec,forceSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(forceSec,MEF90DefMechCtx%cellDMVectScatter,SCATTER_REVERSE,MEF90DefMechCtx%force,ierr);CHKERRQ(ierr)

      Call SectionRealDuplicate(MEF90DefMechCtx%cellDMScalSec,pressureForceSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(pressureForceSec,MEF90DefMechCtx%CellDMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%pressureForce,ierr);CHKERRQ(ierr)

      Work   = 0.0_Kr
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         mywork = 0.0_Kr
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)

         !!! Call proper local assembly depending on the type of damage law
         QuadratureOrder = elemDisplacementType%order
         Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)

         !!! Force work
         If (globalOptions%forceScaling /= MEF90Scaling_Null) Then
            Call MEF90ElasticityWorkSetCell(mywork,xSec,MEF90DefMechCtx%CellDMVect,ForceSec,setIS,elemDisplacement,elemDisplacementType,ierr)
         End If

         !!! pressure force work
         If ((globalOptions%pressureForceScaling /= MEF90Scaling_Null) .AND. (elemDisplacementType%coDim > 0)) Then
            Call MEF90ElasticityPressureWorkSetCell(mywork,xSec,MEF90DefMechCtx%CellDMVect,pressureForceSec,setIS,elemDisplacement,elemDisplacementType,ierr)
         End If
         
         Call MEF90Element_Destroy(elemDisplacement,ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)

         Call MPI_AllReduce(myWork,myWorkSet,1,MPIU_SCALAR,MPI_SUM,MEF90DefMechCtx%MEF90Ctx%comm,ierr);CHKERRQ(ierr)
         work(set) = work(set) + myWorkSet
      End Do ! set
      Call ISrestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr) 

      Call SectionRealDestroy(pressureForceSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(forceSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
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
      Type(Vec),Intent(IN)                               :: xVec
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,Dimension(:),Pointer                     :: cohesiveEnergy
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(SectionReal)                                  :: xSec
      Type(SectionReal)                                  :: boundaryDisplacementSec
      Type(IS)                                           :: CellSetGlobalIS,setIS
      PetscInt,Dimension(:),Pointer                      :: setID
      PetscInt                                           :: set,QuadratureOrder,numCell
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90Element_Type)                            :: elemDisplacementType
      Type(MEF90DefMechGlobalOptions_Type),Pointer       :: globalOptions
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      PetscReal                                          :: mycohesiveEnergy,mycohesiveEnergySet
      
      
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%globalOptionsBag,globalOptions,ierr);CHKERRQ(ierr)
      !!! Create dof-based sections
      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,xSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(xSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,MEF90DefMechCtx%boundaryDisplacement,ierr);CHKERRQ(ierr)

      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,boundaryDisplacementSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(boundaryDisplacementSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,xVec,ierr);CHKERRQ(ierr)

      cohesiveEnergy   = 0.0_Kr
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         mycohesiveEnergy = 0.0_Kr
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)

         !!! Call proper local assembly depending on the type of damage law
         QuadratureOrder = elemDisplacementType%order
         Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)

         Call MEF90ElasticityCohesiveEnergySet(mycohesiveEnergy,xSec,MEF90DefMechCtx%CellDMVect,boundaryDisplacementSec,setIS,elemDisplacement,elemDisplacementType,ierr)

         Call MEF90Element_Destroy(elemDisplacement,ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)

         Call MPI_AllReduce(mycohesiveEnergy,mycohesiveEnergySet,1,MPIU_SCALAR,MPI_SUM,MEF90DefMechCtx%MEF90Ctx%comm,ierr);CHKERRQ(ierr)
         cohesiveEnergy(set) = cohesiveEnergy(set) + mycohesiveEnergySet * matpropSet%cohesiveStiffness
      End Do ! set
      Call ISrestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr) 

      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(boundaryDisplacementSec,ierr);CHKERRQ(ierr)
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
      Type(Vec),Intent(IN)                               :: x
      Type(Vec),Intent(IN)                               :: plasticStrainOld
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,dimension(:),Pointer                     :: energy
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(SectionReal)                                  :: xSec
      Type(SectionReal)                                  :: damageSec,plasticStrainSec,plasticStrainOldSec,temperatureSec,cumulatedDissipatedPlasticEnergySec
      Type(IS)                                           :: CellSetGlobalIS,setIS,setISdof,bcIS
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemScal
      Type(MEF90Element_Type)                            :: elemDisplacementType,elemScalType
      PetscReal                                          :: myenergy,myenergySet

      !!! Get Section and scatter associated with each field
      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,xSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(xSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr) 

      Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,damageSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(damageSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)   

      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,temperatureSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(temperatureSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr)   
      Else
         temperatureSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%CellDMMatSSec,plasticStrainSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticStrainSec,MEF90DefMechCtx%CellDMMatSScatter,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)
         
         Call SectionRealDuplicate(MEF90DefMechCtx%CellDMScalSec,cumulatedDissipatedPlasticEnergySec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(cumulatedDissipatedPlasticEnergySec,MEF90DefMechCtx%CellDMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%cumulatedDissipatedPlasticEnergy,ierr);CHKERRQ(ierr)
      Else
         PlasticStrainSec%v = 0
         cumulatedDissipatedPlasticEnergySec%v = 0
      End If

      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)      
      Do set = 1,size(setID)
         myenergy = 0.0_Kr
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DMVect,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
         elemScalType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)


         !!! Call proper local assembly depending on the type of damage law
         Select Case (cellSetOptions%damageType)
         Case (MEF90DefMech_damageTypeAT1Elastic,MEF90DefMech_damageTypeAT2Elastic,MEF90DefMech_damageTypeAT1,MEF90DefMech_damageTypeLinSoft)
            If (elemDisplacementType%coDim == 0) Then
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  QuadratureOrder = 2 * max(elemDisplacementType%order - 1, elemScalType%order)
               Else
                  QuadratureOrder = 2 * (elemDisplacementType%order - 1)
               End If
               Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
               Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
               Call MEF90PlasticityEnergySet(myenergy,xSec,damageSec,cumulatedDissipatedPlasticEnergySec,MEF90DefMechCtx%DMVect,MEF90DefMechCtx%DMScal,setIS, &
                                             elemDisplacement,elemDisplacementType,elemScal,elemScalType,matpropSet,ierr)
               Call MEF90Element_Destroy(elemDisplacement,ierr)
               Call MEF90Element_Destroy(elemScal,ierr)
            End If
         Case default
            Print*,__FUNCT__,': Unimplemented damage Type, only AT1Elastic and AT2Elastic implement',cellSetOptions%damageType
            STOP  
         End Select


         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         Call MPI_AllReduce(myenergy,myenergySet,1,MPIU_SCALAR,MPI_SUM,MEF90DefMechCtx%MEF90Ctx%comm,ierr);CHKERRQ(ierr)
         energy(set) = energy(set) + myenergySet
      End Do ! set
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
      End If
      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
         Call SectionRealDestroy(cumulatedDissipatedPlasticEnergySec,ierr);CHKERRQ(ierr)
      End If
      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechPlasticDissipation



#undef __FUNCT__
#define __FUNCT__ "MEF90DefmechElasticEnergyKKLSet"
   !!!
   !!!  
   !!!  MEF90DefmechElasticEnergyKKLSet:  Contribution of a cell set to elastic energy. 
   !!!                        It is assumed that the temperature is interpolated on the FE space while the plastic strain 
   !!!                        is cell-based
   !!!  
   !!!  (c) 2012-2019 Blaise Bourdin bourdin@lsu.edu
   !!!
   Subroutine MEF90DefmechElasticEnergyKKLSet(energy,x,damage,plasticStrain,temperature,mesh,meshScal,cellIS,matprop,elemDisplacement,elemDisplacementType,elemDamage,elemDamageType,ierr)
      PetscReal,Intent(OUT)                              :: energy
      Type(SectionReal),Intent(IN)                       :: x,damage,plasticStrain,temperature
      Type(DM),Intent(IN)                                :: mesh,meshScal
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elemDamage
      Type(MEF90Element_Type),Intent(IN)                 :: elemDisplacementType,elemDamageType
      PetscErrorCode,Intent(OUT)                         :: ierr

      PetscReal,Dimension(:),Pointer                     :: xloc,damageLoc,temperatureLoc,plasticStrainLoc
      Type(MEF90_MATS)                                   :: plasticStrainCell
      PetscReal                                          :: damageGauss,temperatureGauss,Stress_ZZ_planeStrain
      PetscReal                                          :: elasticEnergyDensityGauss,stiffness
      Type(MEF90_MATS)                                   :: inelasticStrainGauss,plasticStrainGauss
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscInt                                           :: iDoF1,iGauss
      PetscLogDouble                                     :: flops
     
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(xloc(elemDisplacementType%numDof))
         Allocate(temperatureloc(elemDamageType%numDof))
         Allocate(damageLoc(elemDamageType%numDof))
         Do cell = 1,size(cellID)   
            Call SectionRealRestrictClosure(x,mesh,cellID(cell),elemDisplacementType%numDof,xloc,ierr);CHKERRQ(ierr)
            If (temperature%v /= 0) Then
               Call SectionRealRestrictClosure(temperature,meshScal,cellID(cell),elemDamageType%numDof,temperatureLoc,ierr);CHKERRQ(ierr)
            End If
            plasticStrainCell = 0.0_Kr
            If (plasticStrain%v /= 0) Then
               Call SectionRealRestrict(plasticStrain,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
               plasticStrainCell = plasticStrainLoc
            End If
            If (damage%v /= 0) Then
               Call SectionRealRestrictClosure(damage,meshScal,cellID(cell),elemDamageType%numDof,damageLoc,ierr);CHKERRQ(ierr)
            End If
            Do iGauss = 1,size(elemDisplacement(cell)%Gauss_C)
               stiffness = 1.0_Kr
               If (damage%v /= 0) Then
                  damageGauss = 0.0_Kr
                  Do iDoF1 = 1, elemDamageType%numDof
                     damageGauss = damageGauss + damageLoc(iDoF1) * elemDamage(cell)%BF(iDoF1,iGauss)
                  End Do ! iDoF1
                  stiffness = MEF90gKKL(damageGauss) + matprop%residualStiffness
               End If

               inelasticStrainGauss = 0.0_Kr
               Do iDoF1 = 1,elemDisplacementType%numDof
                  inelasticStrainGauss = inelasticStrainGauss + xloc(iDof1) * elemDisplacement(cell)%GradS_BF(iDof1,iGauss)
               End Do

               temperatureGauss   = 0.0_Kr
               If (temperature%v /= 0) Then
                  Do iDoF1 = 1,elemDamageType%numDof
                     temperatureGauss = temperatureGauss + temperatureLoc(iDof1) * elemDamage(cell)%BF(iDof1,iGauss)
                  End Do
                  inelasticStrainGauss = inelasticStrainGauss - (temperatureGauss * matprop%LinearThermalExpansion)
               End If
               elasticEnergyDensityGauss = 0.5_Kr * (matprop%HookesLaw * (inelasticStrainGauss - plasticStrainCell) .dotP. (inelasticStrainGauss - plasticStrainCell))

#if MEF90_DIM == 2
               if (.NOT. matprop%HookesLaw%isPlaneStress) then
                  Stress_ZZ_planeStrain = ( matprop%HookesLaw%YoungsModulus - 2.0_Kr*matprop%HookesLaw%PoissonRatio*matprop%HookesLaw%mu )*trace(plasticStrainCell) + matprop%HookesLaw%lambda*trace(inelasticStrainGauss)
                  elasticEnergyDensityGauss = elasticEnergyDensityGauss + ( Stress_ZZ_planeStrain + matprop%HookesLaw%lambda*trace(inelasticStrainGauss - plasticStrainCell) ) * trace(inelasticStrainGauss)
               endif
#endif
            energy = energy + elemDamage(cell)%Gauss_C(iGauss) * elasticEnergyDensityGauss
            End Do ! Gauss
            If (plasticStrain%v /= 0) Then
               Call SectionRealRestore(plasticStrain,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
            End If
         End Do ! cell
         flops = 7 * size(elemDisplacement(1)%Gauss_C) * size(cellID) 
         If (damage%v /= 0) Then
            flops = flops + 2 * elemDamageType%numDof * size(elemDamage(1)%Gauss_C) * size(cellID) 
         End If
         If (temperature%v /= 0) Then
            flops = flops + 2 * elemDamageType%numDof * size(elemDamage(1)%Gauss_C) * size(cellID) 
         End If
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(xloc)
         DeAllocate(temperatureloc)
         DeAllocate(damageloc)
      End If   
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefmechElasticEnergyKKLSet


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechElasticEnergy"
!!!
!!!  
!!!  MEF90DefMechElasticEnergy: 
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechElasticEnergy(x,MEF90DefMechCtx,energy,ierr)
      Type(Vec),Intent(IN)                               :: x
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,dimension(:),Pointer                     :: energy
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(SectionReal)                                  :: xSec
      Type(SectionReal)                                  :: damageSec,plasticStrainSec,temperatureSec
      Type(IS)                                           :: CellSetGlobalIS,setIS,setISdof,bcIS
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemScal
      Type(MEF90Element_Type)                            :: elemDisplacementType,elemScalType
      PetscScalar                                        :: myenergy,myenergySet

      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,xSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(xSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr) 

      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%CellDMMatSSec,plasticStrainSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticStrainSec,MEF90DefMechCtx%CellDMMatSScatter,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr) 
      Else
         PlasticStrainSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,TemperatureSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(TemperatureSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr) 
      Else
         temperatureSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%damage)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,damageSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(damageSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr) 
      Else
         damageSec%v = 0
      End If

      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)      
      Do set = 1,size(setID)
         myenergy = 0.0_Kr
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DMVect,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
         elemScalType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)

         !!! Call proper local assembly depending on the type of damage law
         Select Case (cellSetOptions%damageType)
         Case (MEF90DefMech_damageTypeAT1Elastic,MEF90DefMech_damageTypeAT2Elastic,MEF90DefMech_damageTypeKKLElastic)
            If (elemDisplacementType%coDim == 0) Then
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  QuadratureOrder = 2 * max(elemDisplacementType%order - 1, elemScalType%order)
               Else
                  QuadratureOrder = 2 * (elemDisplacementType%order - 1)
               End If
               Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
               Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
               Call MEF90ElasticityEnergySet(myenergy,xSec,plasticStrainSec,temperatureSec,MEF90DefMechCtx%DMVect,MEF90DefMechCtx%DMScal,setIS, &
                                             matpropSet%HookesLaw,matpropSet%LinearThermalExpansion,elemDisplacement,elemDisplacementType,elemScal,elemScalType,ierr)
               Call MEF90Element_Destroy(elemDisplacement,ierr)
               Call MEF90Element_Destroy(elemScal,ierr)
            End If
         Case (MEF90DefMech_damageTypeLinSoft)
            If (elemDisplacementType%coDim == 0) Then
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  QuadratureOrder = 5 * max(elemDisplacementType%order - 1, elemScalType%order)
               Else
                  QuadratureOrder = 5 ! * (elemDisplacementType%order - 1)
               End If
               Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
               Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
               Call MEF90GradDamageElasticEnergySet(myenergy,xSec,damageSec,plasticStrainSec,temperatureSec,MEF90DefMechCtx%DMVect,MEF90DefMechCtx%DMScal,setIS, &
                                             matPropSet%residualStiffness,matpropSet%HookesLaw,matpropSet%CoefficientLinSoft,matpropSet%LinearThermalExpansion,elemDisplacement,elemDisplacementType,elemScal,elemScalType,ierr)
               Call MEF90Element_Destroy(elemDisplacement,ierr)
               Call MEF90Element_Destroy(elemScal,ierr)
            End If
         Case (MEF90DefMech_damageTypeAT1,MEF90DefMech_damageTypeAT2)
            If (elemDisplacementType%coDim == 0) Then
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  QuadratureOrder = 2 * max(elemDisplacementType%order - 1, elemScalType%order) + 2 * elemScalType%order
               Else
                  QuadratureOrder = 2 * (elemDisplacementType%order - 1) + 2 * elemScalType%order
               End If
               Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
               Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
               
               matpropSet%CoefficientLinSoft=0
               Call MEF90GradDamageElasticEnergySet(myenergy,xSec,damageSec,plasticStrainSec,temperatureSec,MEF90DefMechCtx%DMVect,MEF90DefMechCtx%DMScal,setIS, &
                                             matPropSet%residualStiffness,matpropSet%HookesLaw,matpropSet%CoefficientLinSoft,matpropSet%LinearThermalExpansion,elemDisplacement,elemDisplacementType,elemScal,elemScalType,ierr)
               Call MEF90Element_Destroy(elemDisplacement,ierr)
               Call MEF90Element_Destroy(elemScal,ierr)
            End If
         Case (MEF90DefMech_damageTypeKKL)
            If (elemDisplacementType%coDim == 0) Then
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  QuadratureOrder = 2 * max(elemDisplacementType%order - 1, elemScalType%order) + 2 * elemScalType%order
               Else
                  QuadratureOrder = 2 * (elemDisplacementType%order - 1) + 2 * elemScalType%order
               End If
               Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
               Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
               
               matpropSet%CoefficientLinSoft=0
               Call MEF90DefmechElasticEnergyKKLSet(myenergy,xSec,damageSec,plasticStrainSec,temperatureSec,MEF90DefMechCtx%DMVect,MEF90DefMechCtx%DMScal,setIS, &
                                             matPropSet,elemDisplacement,elemDisplacementType,elemScal,elemScalType,ierr)
               Call MEF90Element_Destroy(elemDisplacement,ierr)
               Call MEF90Element_Destroy(elemScal,ierr)
            End If
         !Case default
         !   Print*,__FUNCT__,': Unimplemented damage Type',cellSetOptions%damageType
         !   STOP  
         End Select
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         Call MPI_AllReduce(myenergy,myenergySet,1,MPIU_SCALAR,MPI_SUM,MEF90DefMechCtx%MEF90Ctx%comm,ierr);CHKERRQ(ierr)
         energy(set) = energy(set) + myenergySet
      End Do ! set
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      If (Associated(MEF90DefMechCtx%damage)) Then
         Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)
      End If
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
      End If
      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
      End If
      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechElasticEnergy



#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechStress"
!!!
!!!  
!!!  MEF90DefMechStress: 
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechStress(x,MEF90DefMechCtx,stress,ierr)
      Type(Vec),Intent(IN)                               :: x
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      Type(Vec),Intent(IN)                               :: stress
         
      Type(SectionReal)                                  :: xSec,stressSec
      Type(SectionReal)                                  :: plasticStrainSec,temperatureSec,damageSec
      Type(IS)                                           :: CellSetGlobalIS,setIS,setISdof,bcIS
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemScal
      Type(MEF90Element_Type)                            :: elemDisplacementType,elemScalType

      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,xSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(xSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr) 

      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,temperatureSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(temperatureSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr)   
      Else
         temperatureSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%CellDMMatSSec,plasticStrainSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticStrainSec,MEF90DefMechCtx%CellDMMatSScatter,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)
      Else
         PlasticStrainSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%damage)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,damageSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(damageSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)   
      Else
         damageSec%v = 0
      End If

      Call SectionRealDuplicate(MEF90DefMechCtx%CellDMMatSSec,stressSec,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)      
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DMVect,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
         elemScalType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)

         !!! Call proper local assembly depending on the type of damage law
         Select Case (cellSetOptions%damageType)

         Case (MEF90DefMech_damageTypeLinSoft)
            If (elemDisplacementType%coDim == 0) Then
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  QuadratureOrder = max(elemDisplacementType%order - 1, elemScalType%order)
               Else
                  QuadratureOrder = 5 !*elemDisplacementType%order - 1
               End If
               Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
               Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
               
               Call MEF90ElasticityStressSet(stressSec,xSec,plasticStrainSec,temperatureSec,damageSec,matpropSet%CoefficientLinSoft,MEF90DefMechCtx%DMVect,MEF90DefMechCtx%DMScal,setIS, &
                                             matpropSet%HookesLaw,matpropSet%LinearThermalExpansion,elemDisplacement,elemDisplacementType,elemScal,elemScalType,ierr)

               Call MEF90Element_Destroy(elemDisplacement,ierr)
               Call MEF90Element_Destroy(elemScal,ierr)
            End If

         Case (MEF90DefMech_damageTypeAT1Elastic,MEF90DefMech_damageTypeAT2Elastic, &
               MEF90DefMech_damageTypeAT1,MEF90DefMech_damageTypeAT2)
            If (elemDisplacementType%coDim == 0) Then
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  QuadratureOrder = max(elemDisplacementType%order - 1, elemScalType%order)
               Else
                  QuadratureOrder = elemDisplacementType%order - 1
               End If
               Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
               Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)

               matpropSet%CoefficientLinSoft = 0
               Call MEF90ElasticityStressSet(stressSec,xSec,plasticStrainSec,temperatureSec,damageSec,matpropSet%CoefficientLinSoft,MEF90DefMechCtx%DMVect,MEF90DefMechCtx%DMScal,setIS, &
                                             matpropSet%HookesLaw,matpropSet%LinearThermalExpansion,elemDisplacement,elemDisplacementType,elemScal,elemScalType,ierr)

               Call MEF90Element_Destroy(elemDisplacement,ierr)
               Call MEF90Element_Destroy(elemScal,ierr)
            End If
         Case (MEF90DefMech_damageTypeKKLElastic,MEF90DefMech_damageTypeKKL)
            Print*,__FUNCT__,': Unimplemented damage Type',cellSetOptions%damageType
         Case default
            Print*,__FUNCT__,': Unimplemented damage Type',cellSetOptions%damageType
            STOP  
         End Select
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do ! set

      !!! No need to complete since stressSec is cell-based
      Call SectionRealToVec(stressSec,MEF90DefMechCtx%CellDMMatSScatter,SCATTER_FORWARD,stress,ierr);CHKERRQ(ierr)          
      
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
      End If

      If (Associated(MEF90DefMechCtx%damage)) Then
         Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)
      End If

      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
      End If
      Call SectionRealDestroy(stressSec,ierr);CHKERRQ(ierr)

   End Subroutine MEF90DefMechStress

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamageAT1Loc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamageAT1Loc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDamageAT1Loc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: elasticEnergyDensityGauss,temperatureGauss,damageGauss,Stress_ZZ_planeStrain
      Type(MEF90_MATS)                                   :: inelasticStrainGauss
      PetscReal                                          :: C2,N
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)
      
      N = matprop%DuctileCouplingPower
      C2 = matprop%fractureToughness * matprop%internalLength * .75
      Aloc = 0.0_Kr
      Do iGauss = 1,numGauss
         temperatureGauss = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
            End Do
         End If

         inelasticStrainGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
         End Do
         inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion 
         elasticEnergyDensityGauss = ( matprop%HookesLaw * (inelasticStrainGauss - plasticStrainCell) ) .DotP. (inelasticStrainGauss - plasticStrainCell)

#if MEF90_DIM == 2
         !!! Adding terms in planestrain for plasticity with tr(p) = 0
         If (.NOT. matProp%HookesLaw%isPlaneStress) Then
            Stress_ZZ_planeStrain = ( matprop%HookesLaw%YoungsModulus - 2.0_Kr*matprop%HookesLaw%PoissonRatio*matprop%HookesLaw%mu )*trace(plasticStrainCell) + matprop%HookesLaw%lambda*trace(inelasticStrainGauss)
            elasticEnergyDensityGauss = elasticEnergyDensityGauss + ( Stress_ZZ_planeStrain + matprop%HookesLaw%lambda*trace(inelasticStrainGauss - plasticStrainCell) ) * trace(plasticStrainCell)
         End If
#endif
         !!! This is really twice the elastic energy density

         If ((N == 1.0_Kr) .OR. (cumulatedDissipatedPlasticEnergyCell == 0.0_Kr)) Then
            Do iDoF1 = 1,numDofDamage
               Do iDoF2 = 1,numDofDamage
                  Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) * ( &
                                         (elasticEnergyDensityGauss  ) * elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss) + &
                                         C2 * ((matprop%toughnessAnisotropyMatrix * elemDamage%Grad_BF(iDoF1,iGauss)) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
               End Do
            End Do
         Else
            damageGauss = 0.0_Kr
            Do iDoF1 = 1,numDofDamage
               damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
            End Do
            Do iDoF1 = 1,numDofDamage
               Do iDoF2 = 1,numDofDamage
                  Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) * ( &
                                         (elasticEnergyDensityGauss + &
                                         ( N * (N - 1.0_Kr) *( (1.0_Kr-damageGauss)**(N - 2.0_Kr) ) ) * cumulatedDissipatedPlasticEnergyCell ) * elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss) + &
                                         C2 * ((matprop%toughnessAnisotropyMatrix * elemDamage%Grad_BF(iDoF1,iGauss)) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
               End Do
            End Do
         EndIf

      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDamageAT1Loc


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamageLinSoftLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamageLinSoftLoc: Assembles the bilinear form for LinSoft, that is:
!!!    <K(alpha) beta , delta> = \int a''(\alpha) W(e(u)) \beta \delta 
!!!                            + Gc/pi/\ell w''(alpha) \beta \delta
!!!                            + 2 Gc/pi*\ell \nabla \beta \nabla \delta
!!! with a(\alpha) = (1-w(\alpha)) / ( 1 + (k-1) w(\alpha)) and w(\alpha) = 1 - (1-\alpha)^2
!!! i.e. a''(\alpha) = 2 k * ( k + 3(k-1) v^2) / [k + (1-k) v^2]^3, w''(\alpha) = -2
!!!  
!!!  (c) 2014-2015 Erwan Tanne erwan.tanne@gmail.com, Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDamageLinSoftLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: elasticEnergyDensityGauss,temperatureGauss,damageGauss
      Type(MEF90_MATS)                                   :: inelasticStrainGauss
      PetscReal                                          :: C2,C1,N,k,C0
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)
      
      k = matprop%CoefficientLinSoft
      C1 = 2.0_Kr * matprop%fractureToughness / matprop%internalLength / PETSC_PI
      C2 = 2.0_Kr * matprop%fractureToughness * matprop%internalLength / PETSC_PI
      N  = matprop%DuctileCouplingPower

      Aloc = 0.0_Kr
      Do iGauss = 1,numGauss
         temperatureGauss = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
            End Do
         End If

         inelasticStrainGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
         End Do
         inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion - plasticStrainCell
         elasticEnergyDensityGauss = (matprop%HookesLaw * inelasticStrainGauss) .DotP. inelasticStrainGauss

         damageGauss = 0.0_Kr
         Do iDof1 = 1, numDofDamage
            damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
         End Do
         C0 = k * (k + 3.0_Kr * (k - 1.0_Kr) * (1.0_Kr - damageGauss)**2) / &
              (k + (1.0_Kr - k) * (1.0_Kr - damageGauss)**2)**3
         !!! This is really twice the elastic energy density

         Do iDoF1 = 1,numDofDamage
            Do iDoF2 = 1,numDofDamage
               Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) * ( & 
                                  (elasticEnergyDensityGauss * C0 - C1 + &
                                  (N * (N - 1.0) *( (1.0_Kr-damageGauss)**(N - 2.0) ) ) * cumulatedDissipatedPlasticEnergyCell ) * elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss) &
                                  + C2 * (elemDamage%Grad_BF(iDoF1,iGauss) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
            End Do
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDamageLinSoftLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamageAT1ElasticLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamageAT1ElasticLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDamageAT1ElasticLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDamage,numGauss
      PetscReal                                          :: C2 
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDisplacement%BF,2)

      C2 = matprop%fractureToughness * matprop%internalLength * .75_Kr
      Aloc = 0.0_Kr
      Do iGauss = 1,numGauss
         Do iDoF1 = 1,numDofDamage
            Do iDoF2 = 1,numDofDamage
               Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) *  &
                                      C2 * (elemDamage%Grad_BF(iDoF1,iGauss) .dotP. elemDamage%Grad_BF(iDoF2,iGauss))
            End Do
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDamageAT1ElasticLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamageKKLLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamageKKLLoc:
!!!  
!!!  (c) 2019 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDamageKKLLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: elasticEnergyDensityGauss,temperatureGauss,damageGauss,Stress_ZZ_planeStrain
      Type(MEF90_MATS)                                   :: inelasticStrainGauss
      PetscReal                                          :: C1,C2,N,D2gDamageGauss
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)
      
      N = matprop%DuctileCouplingPower
      C1 = matprop%fractureToughness / matprop%internalLength / 4.0_Kr / cwKKL
      C2 = matprop%fractureToughness * matprop%internalLength / 2.0_Kr / cwKKL
      Aloc = 0.0_Kr
      Do iGauss = 1,numGauss
         temperatureGauss = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
            End Do
         End If

         inelasticStrainGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
         End Do
         inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion 
         elasticEnergyDensityGauss = ( matprop%HookesLaw * (inelasticStrainGauss - plasticStrainCell) ) .DotP. (inelasticStrainGauss - plasticStrainCell)

#if MEF90_DIM == 2
         !!! Adding terms in planestrain for plasticity with tr(p) = 0
         If (.NOT. matProp%HookesLaw%isPlaneStress) Then
            Stress_ZZ_planeStrain = ( matprop%HookesLaw%YoungsModulus - 2.0_Kr*matprop%HookesLaw%PoissonRatio*matprop%HookesLaw%mu )*trace(plasticStrainCell) + matprop%HookesLaw%lambda*trace(inelasticStrainGauss)
            elasticEnergyDensityGauss = elasticEnergyDensityGauss + ( Stress_ZZ_planeStrain + matprop%HookesLaw%lambda*trace(inelasticStrainGauss - plasticStrainCell) ) * trace(plasticStrainCell)
         End If
#endif
         elasticEnergyDensityGauss = elasticEnergyDensityGauss * 0.5_Kr

         damageGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
         End Do
         D2gDamageGauss = MEF90D2gKKL(damageGauss)

         If ((N == 1.0_Kr) .OR. (cumulatedDissipatedPlasticEnergyCell == 0.0_Kr)) Then
            Do iDoF1 = 1,numDofDamage
               Do iDoF2 = 1,numDofDamage
                  Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) * ( &
                                         (elasticEnergyDensityGauss - C1) * D2gDamageGauss * elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss) + &
                                         C2 * ((matprop%toughnessAnisotropyMatrix * elemDamage%Grad_BF(iDoF1,iGauss)) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
               End Do
            End Do
         Else
            Do iDoF1 = 1,numDofDamage
               Do iDoF2 = 1,numDofDamage
                  Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) * ( &
                                         ((elasticEnergyDensityGauss - C1) * D2gDamageGauss +  &
                                         ( N * (N - 1.0_Kr) *( (1.0_Kr-damageGauss)**(N - 2.0_Kr) ) ) * cumulatedDissipatedPlasticEnergyCell ) * elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss) + &
                                         C2 * ((matprop%toughnessAnisotropyMatrix * elemDamage%Grad_BF(iDoF1,iGauss)) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
               End Do
            End Do
         EndIf

      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDamageKKLLoc


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamageKKLElasticLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamageKKLElasticLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDamageKKLElasticLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDamage,numGauss
      PetscReal                                          :: C1,C2,damageGauss,D2gDamageGauss
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDisplacement%BF,2)

      C1 = matprop%fractureToughness / matprop%internalLength / 4.0_Kr / cwKKL
      C2 = matprop%fractureToughness * matprop%internalLength / 2.0_Kr / cwKKL
      Aloc = 0.0_Kr
      Do iGauss = 1,numGauss
         damageGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
         End Do
         D2gDamageGauss = MEF90D2gKKL(damageGauss)
         Do iDoF1 = 1,numDofDamage
            Do iDoF2 = 1,numDofDamage
               Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) *  ( &
                                      -C1 * D2gDamageGauss * elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss) + &
                                      C2 * (elemDamage%Grad_BF(iDoF1,iGauss) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
            End Do
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDamageKKLElasticLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamageAT1UnilateralHDLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamageAT1UnilateralHDLoc:
!!!  
!!!  (c) 2016 - 2018 Blaise Bourdin bourdin@lsu.edu, Erwan Tanne erwan.tanne@gmail.com
!!!

   Subroutine MEF90DefMechBilinearFormDamageAT1UnilateralHDLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: elasticEnergyDensityGauss,damageGauss,temperatureGauss,Stress_ZZ_planeStrain
      Type(MEF90_MATS)                                   :: inelasticStrainGauss
      PetscReal                                          :: C2,N,kappa
      PetscLogDouble                                     :: flops

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)
      kappa = ((matProp%HookesLaw * MEF90_MATS_IDENTITY) .dotP. MEF90_MATS_IDENTITY) / MEF90_DIM**2
      
      N = matprop%DuctileCouplingPower
      C2 = matprop%fractureToughness * matprop%internalLength * .75_Kr
      Aloc = 0.0_Kr
      Do iGauss = 1,numGauss
         temperatureGauss = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
            End Do
         End If

         inelasticStrainGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
         End Do
         inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion 
         elasticEnergyDensityGauss = ((matprop%HookesLaw * DeviatoricPart(inelasticStrainGauss - plasticStrainCell)) .DotP. (inelasticStrainGauss - plasticStrainCell)) + &
                                      kappa * MEF90HydrostaticDeviatoricPenaltyFunction(Trace(inelasticStrainGauss - plasticStrainCell))
         !!! This is really twice the elastic energy density

#if MEF90_DIM == 2
         !!! Adding terms in planestrain for plasticity with tr(p) = 0
         !!! Not completely sure about this in the context of HD decomposition
         !If (.NOT. matProp%HookesLaw%isPlaneStress) Then
         !   Stress_ZZ_planeStrain = ( matprop%HookesLaw%YoungsModulus - 2.0_Kr*matprop%HookesLaw%PoissonRatio*matprop%HookesLaw%mu )*trace(plasticStrainCell) + matprop%HookesLaw%lambda*trace(inelasticStrainGauss)
         !   elasticEnergyDensityGauss = elasticEnergyDensityGauss + ( Stress_ZZ_planeStrain + matprop%HookesLaw%lambda*trace(inelasticStrainGauss - plasticStrainCell) ) * trace(plasticStrainCell)
         !End If
#endif

         If ((N == 1.0_Kr) .OR. (cumulatedDissipatedPlasticEnergyCell == 0.0_Kr)) Then
            Do iDoF1 = 1,numDofDamage
               Do iDoF2 = 1,numDofDamage
                  Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) * ( &
                                         elasticEnergyDensityGauss * elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss) + &
                                         C2 * (elemDamage%Grad_BF(iDoF1,iGauss) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
               End Do
            End Do
         Else
            damageGauss = 0.0_Kr
            Do iDoF1 = 1,numDofDamage
               damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
            End Do
            Do iDoF1 = 1,numDofDamage
               Do iDoF2 = 1,numDofDamage
                  Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) * ( &
                                         (elasticEnergyDensityGauss + &
                                         ( N * (N - 1.0_Kr) *( (1.0_Kr-damageGauss)**(N - 2.0_Kr) ) ) * cumulatedDissipatedPlasticEnergyCell ) * elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss) + &
                                         C2 * (elemDamage%Grad_BF(iDoF1,iGauss) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
               End Do
            End Do
         EndIf

      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDamageAT1UnilateralHDLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamageAT1UnilateralDeviatoricLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamageAT1UnilateralDeviatoricLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDamageAT1UnilateralDeviatoricLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: elasticEnergyDensityGauss,temperatureGauss
      Type(MEF90_MATS)                                   :: inelasticStrainGauss
      PetscReal                                          :: C2
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)
      
      C2 = matprop%fractureToughness * matprop%internalLength * .75
      Aloc = 0.0_Kr
      Do iGauss = 1,numGauss
         temperatureGauss = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
            End Do
         End If

         inelasticStrainGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
         End Do
         inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion - plasticStrainCell

         elasticEnergyDensityGauss = (matprop%HookesLaw * inelasticStrainGauss) .DotP. DeviatoricPart(inelasticStrainGauss)
         !!! This is really twice the elastic energy density
         
         Do iDoF1 = 1,numDofDamage
            Do iDoF2 = 1,numDofDamage
               Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) * ( &
                                      elasticEnergyDensityGauss * elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss) + &
                                      C2 * (elemDamage%Grad_BF(iDoF1,iGauss) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
            End Do
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDamageAT1UnilateralDeviatoricLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamageAT1UnilateralPHDLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamageAT1UnilateralPHDLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDamageAT1UnilateralPHDLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: elasticEnergyDensityGauss,temperatureGauss
      Type(MEF90_MATS)                                   :: inelasticStrainGauss
      PetscReal                                          :: C2
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)
      
      C2 = matprop%fractureToughness * matprop%internalLength * .75
      Aloc = 0.0_Kr
      Do iGauss = 1,numGauss
         temperatureGauss = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
            End Do
         End If

         inelasticStrainGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
         End Do
         inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion - plasticStrainCell

         If (Trace(inelasticStrainGauss) < 0.0_Kr) Then
            elasticEnergyDensityGauss = 0.0_Kr
         Else
            elasticEnergyDensityGauss = (matprop%HookesLaw * inelasticStrainGauss) .DotP. inelasticStrainGauss
         End If
         !!! This is really twice the elastic energy density
         
         Do iDoF1 = 1,numDofDamage
            Do iDoF2 = 1,numDofDamage
               Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) * ( &
                                      elasticEnergyDensityGauss * elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss) + &
                                      C2 * (elemDamage%Grad_BF(iDoF1,iGauss) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
            End Do
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDamageAT1UnilateralPHDLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamageAT1UnilateralMasonryLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamageAT1UnilateralMasonryLoc:
!!!  
!!!  (c) 2016 Blaise Bourdin bourdin@lsu.edu, Erwan Tanne erwan.tanne@gmail.com
!!!

   Subroutine MEF90DefMechBilinearFormDamageAT1UnilateralMasonryLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: positiveElasticEnergyDensityGauss
      Type(MEF90_MATS)                                   :: inelasticStrain,positiveInelasticStrain,negativeInelasticStrain
      PetscReal                                          :: C2
      PetscLogDouble                                     :: flops

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)
      
      C2 = matprop%fractureToughness * matprop%internalLength * .75
      Aloc = 0.0_Kr
      Do iGauss = 1,numGauss
         inelasticStrain = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               inelasticStrain = inelasticStrain - temperatureDof(iDoF1) * (matProp%LinearThermalExpansion * elemDamage%BF(iDoF1,iGauss))
            End Do         
         End If
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrain = inelasticStrain + displacementDof(iDoF1) * elemDisplacement%GradS_BF(iDoF1,iGauss)
         End Do
         inelasticStrain = inelasticStrain - plasticStrainCell
         Call MasonryProjection(inelasticStrain,matprop%HookesLaw,positiveInelasticStrain,negativeInelasticStrain)
         positiveElasticEnergyDensityGauss = (matprop%HookesLaw * positiveInelasticStrain) .dotP. positiveInelasticStrain
         !!! This is really twice the elastic energy density
         Do iDoF1 = 1,numDofDamage
            Do iDoF2 = 1,numDofDamage
               Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) * ( &
                                      positiveElasticEnergyDensityGauss * elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss) + &
                                      C2 * (elemDamage%Grad_BF(iDoF1,iGauss) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
            End Do
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDamageAT1UnilateralMasonryLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamageAT2Loc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamageAT2Loc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDamageAT2Loc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: elasticEnergyDensityGauss,temperatureGauss
      Type(MEF90_MATS)                                   :: inelasticStrainGauss
      PetscReal                                          :: C1,C2
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDisplacement%BF,2)
      
      C1 = matprop%fractureToughness / matprop%internalLength 
      C2 = matprop%fractureToughness * matprop%internalLength
      Aloc = 0.0_Kr
      Do iGauss = 1,numGauss
         temperatureGauss = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
            End Do
         End If

         inelasticStrainGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
         End Do
         inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion - plasticStrainCell
         elasticEnergyDensityGauss = (matprop%HookesLaw * inelasticStrainGauss) .DotP. inelasticStrainGauss
         !!! Again, this is really twice the elastic energy density

         Do iDoF1 = 1,numDofDamage
            Do iDoF2 = 1,numDofDamage
               Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) * ( &
                                      (elasticEnergyDensityGauss + C1) * elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss) + &
                                       C2 * ((matprop%toughnessAnisotropyMatrix * elemDamage%Grad_BF(iDoF1,iGauss)) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
            End Do
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDamageAT2Loc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamageAT2ElasticLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamageAT2ElasticLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDamageAT2ElasticLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDamage,numGauss
      PetscReal                                          :: C1,C2
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDisplacement%BF,2)
      
      C1 = matprop%fractureToughness / matprop%internalLength 
      C2 = matprop%fractureToughness * matprop%internalLength
      Aloc = 0.0_Kr
      Do iGauss = 1,numGauss
         Do iDoF1 = 1,numDofDamage
            Do iDoF2 = 1,numDofDamage
               Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) * ( &
                                       C1 * elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss) + & 
                                       C2 * ((matprop%toughnessAnisotropyMatrix * elemDamage%Grad_BF(iDoF1,iGauss)) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
            End Do
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDamageAT2ElasticLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamageAT2UnilateralHDLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamageAT2UnilateralHDLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDamageAT2UnilateralHDLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: elasticEnergyDensityGauss,temperatureGauss,damageGauss,Stress_ZZ_planeStrain
      Type(MEF90_MATS)                                   :: inelasticStrainGauss
      PetscReal                                          :: C1,C2,N,kappa
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDisplacement%BF,2)
      kappa = ((matProp%HookesLaw * MEF90_MATS_IDENTITY) .dotP. MEF90_MATS_IDENTITY) / MEF90_DIM**2
      
      N = matprop%DuctileCouplingPower
      C1 = matprop%fractureToughness / matprop%internalLength 
      C2 = matprop%fractureToughness * matprop%internalLength
      Aloc = 0.0_Kr
      Do iGauss = 1,numGauss
         temperatureGauss = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
            End Do
         End If

         inelasticStrainGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
         End Do
         inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion - plasticStrainCell

         inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion 
         elasticEnergyDensityGauss = ((matprop%HookesLaw * DeviatoricPart(inelasticStrainGauss - plasticStrainCell)) .DotP. (inelasticStrainGauss - plasticStrainCell)) + &
                                      kappa * MEF90HydrostaticDeviatoricPenaltyFunction(Trace(inelasticStrainGauss - plasticStrainCell))
         !!! This is really twice the elastic energy density

#if MEF90_DIM == 2
         !!! Adding terms in planestrain for plasticity with tr(p) = 0
         !!! Not completely sure about this in the context of HD decomposition
         !If (.NOT. matProp%HookesLaw%isPlaneStress) Then
         !   Stress_ZZ_planeStrain = ( matprop%HookesLaw%YoungsModulus - 2.0_Kr*matprop%HookesLaw%PoissonRatio*matprop%HookesLaw%mu )*trace(plasticStrainCell) + matprop%HookesLaw%lambda*trace(inelasticStrainGauss)
         !   elasticEnergyDensityGauss = elasticEnergyDensityGauss + ( Stress_ZZ_planeStrain + matprop%HookesLaw%lambda*trace(inelasticStrainGauss - plasticStrainCell) ) * trace(plasticStrainCell)
         !End If
#endif

         If ((N == 1.0_Kr) .OR. (cumulatedDissipatedPlasticEnergyCell == 0.0_Kr)) Then
            Do iDoF1 = 1,numDofDamage
               Do iDoF2 = 1,numDofDamage
                  Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) * ( &
                                         (elasticEnergyDensityGauss + C1) * elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss) + &
                                         C2 * (elemDamage%Grad_BF(iDoF1,iGauss) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
               End Do
            End Do
         Else
            damageGauss = 0.0_Kr
            Do iDoF1 = 1,numDofDamage
               damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
            End Do
            Do iDoF1 = 1,numDofDamage
               Do iDoF2 = 1,numDofDamage
                  Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) * ( &
                                         (elasticEnergyDensityGauss + C1 + &
                                         ( N * (N - 1.0_Kr) *( (1.0_Kr-damageGauss)**(N - 2.0_Kr) ) ) * cumulatedDissipatedPlasticEnergyCell ) * elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss) + &
                                         C2 * (elemDamage%Grad_BF(iDoF1,iGauss) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
               End Do
            End Do
         EndIf
      End Do

      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDamageAT2UnilateralHDLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamageAT2UnilateralDeviatoricLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamageAT2UnilateralDeviatoricLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDamageAT2UnilateralDeviatoricLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: elasticEnergyDensityGauss,temperatureGauss
      Type(MEF90_MATS)                                   :: inelasticStrainGauss
      PetscReal                                          :: C1,C2
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDisplacement%BF,2)
      
      C1 = matprop%fractureToughness / matprop%internalLength 
      C2 = matprop%fractureToughness * matprop%internalLength
      Aloc = 0.0_Kr
      Do iGauss = 1,numGauss
         temperatureGauss = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
            End Do
         End If

         inelasticStrainGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
         End Do
         inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion - plasticStrainCell

         elasticEnergyDensityGauss = (matprop%HookesLaw * inelasticStrainGauss) .DotP. DeviatoricPart(inelasticStrainGauss)
         !!! Again, this is really twice the elastic energy density

         Do iDoF1 = 1,numDofDamage
            Do iDoF2 = 1,numDofDamage
               Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) * ( &
                                      (elasticEnergyDensityGauss + C1) * elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss) + &
                                       C2 * (elemDamage%Grad_BF(iDoF1,iGauss) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
            End Do
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDamageAT2UnilateralDeviatoricLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamageAT2UnilateralPHDLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamageAT2UnilateralPHDLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDamageAT2UnilateralPHDLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: elasticEnergyDensityGauss,temperatureGauss
      Type(MEF90_MATS)                                   :: inelasticStrainGauss
      PetscReal                                          :: C1,C2
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDisplacement%BF,2)
      
      C1 = matprop%fractureToughness / matprop%internalLength 
      C2 = matprop%fractureToughness * matprop%internalLength
      Aloc = 0.0_Kr
      Do iGauss = 1,numGauss
         temperatureGauss = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
            End Do
         End If

         inelasticStrainGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
         End Do
         inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion - plasticStrainCell

         If (Trace(inelasticStrainGauss) < 0.0_Kr) Then
            elasticEnergyDensityGauss = 0.0_Kr
         Else
            elasticEnergyDensityGauss = (matprop%HookesLaw * inelasticStrainGauss) .DotP. inelasticStrainGauss
         End If
         !!! Again, this is really twice the elastic energy density

         Do iDoF1 = 1,numDofDamage
            Do iDoF2 = 1,numDofDamage
               Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) * ( &
                                      (elasticEnergyDensityGauss + C1) * elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss) + &
                                       C2 * (elemDamage%Grad_BF(iDoF1,iGauss) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
            End Do
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDamageAT2UnilateralPHDLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamageAT1Loc"
!!!
!!!  
!!!  MEF90DefMechOperatorDamageAT1Loc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDamageAT1Loc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      PetscReal,Intent(IN)                               :: CrackPressureCell

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: elasticEnergyDensityGauss,temperatureGauss,damageGauss,Stress_ZZ_planeStrain
      Type(MEF90_MATS)                                   :: inelasticStrainGauss
      Type(MEF90_VECT)                                   :: gradientDamageGauss,Displacement
      PetscReal                                          :: C1,C2,N
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)
      
      N  = matprop%DuctileCouplingPower
      C1 = matprop%fractureToughness / matprop%internalLength * .375
      C2 = matprop%fractureToughness * matprop%internalLength * .75

      residualLoc = 0.0_Kr
      Do iGauss = 1,numGauss
         temperatureGauss = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
            End Do
         End If

         Displacement = 0.0_Kr
         inelasticStrainGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDisplacement
            Displacement = Displacement + elemDisplacement%BF(iDoF1,iGauss) * displacementDof(iDoF1)
            inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
         End Do
         inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion 
         elasticEnergyDensityGauss = (matprop%HookesLaw * (inelasticStrainGauss - plasticStrainCell) ) .DotP. (inelasticStrainGauss - plasticStrainCell)

#if MEF90_DIM == 2
         !!! Adding terms in planestrain for plasticity with tr(p) = 0
         If (.NOT. matProp%HookesLaw%isPlaneStress) Then
           Stress_ZZ_planeStrain = ( matprop%HookesLaw%YoungsModulus - 2.0_Kr*matprop%HookesLaw%PoissonRatio*matprop%HookesLaw%mu )*trace(plasticStrainCell) + matprop%HookesLaw%lambda*trace(inelasticStrainGauss)
           elasticEnergyDensityGauss = elasticEnergyDensityGauss + ( Stress_ZZ_planeStrain + matprop%HookesLaw%lambda*trace(inelasticStrainGauss - plasticStrainCell) ) * trace(plasticStrainCell)
         End If
#endif
         
         !!! This is really twice the elastic energy density
         damageGauss = 0.0_Kr
         gradientDamageGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
            gradientDamageGauss = gradientDamageGauss + elemDamage%Grad_BF(iDoF1,iGauss) * xDof(iDoF1)
         End Do

         If (N==1.0_Kr) Then
            Do iDoF2 = 1,numDofDamage
               residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * ( &
                                    ( elasticEnergyDensityGauss * (damageGauss - 1.0_Kr) &
                                    - ( cumulatedDissipatedPlasticEnergyCell ) + C1 ) * elemDamage%BF(iDoF2,iGauss) + &
                                     (( C2 * matprop%toughnessAnisotropyMatrix * gradientDamageGauss +  Displacement*CrackPressureCell ) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
            End Do
         Else
            Do iDoF2 = 1,numDofDamage
               residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * ( &
                                    ( elasticEnergyDensityGauss * (damageGauss - 1.0_Kr) &
                                    - ( N * cumulatedDissipatedPlasticEnergyCell * (1.0_Kr - damageGauss)**(N - 1.0_Kr)  ) + C1 ) * elemDamage%BF(iDoF2,iGauss) + &
                                     (( C2 * matprop%toughnessAnisotropyMatrix * gradientDamageGauss +  Displacement*CrackPressureCell ) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
            End Do
         EndIf
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDamageAT1Loc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamageLinSoftLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDamageLinSoftLoc: Assembles the operator for LinSoft:
!!!    F(\alpha)\beta = \int [k(\alpha-1)W(e(u))] / [k-(k-1)(1-\alpha)^2]^2 \beta 
!!!                   + 2Gc/pi (1-\alpha) \beta / \ell 
!!!                   + 2Gc/pi \nabla \alpha \nabla \beta * \ell
!!!  
!!!  (c) 2014-2015 Erwan Tanne erwan.tanne@gmail.com, Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDamageLinSoftLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      PetscReal,Intent(IN)                               :: CrackPressureCell

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: elasticEnergyDensityGauss,temperatureGauss,damageGauss
      Type(MEF90_MATS)                                   :: inelasticStrainGauss
      Type(MEF90_VECT)                                   :: gradientDamageGauss
      PetscReal                                          :: C1,C2,k,C0Operator,N
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr


      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)
      
      C1 = 2.0_Kr * matprop%fractureToughness / matprop%internalLength / PETSC_PI
      C2 = 2.0_Kr * matprop%fractureToughness * matprop%internalLength / PETSC_PI
      k = matprop%CoefficientLinSoft
      N = matprop%DuctileCouplingPower

      residualLoc = 0.0_Kr
      Do iGauss = 1,numGauss
         temperatureGauss = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
            End Do
         End If

         inelasticStrainGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
         End Do
         inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion - plasticStrainCell
         elasticEnergyDensityGauss = (matprop%HookesLaw * inelasticStrainGauss) .DotP. inelasticStrainGauss
         !!! This is really twice the elastic energy density

         damageGauss = 0.0_Kr
         gradientDamageGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
            gradientDamageGauss = gradientDamageGauss + elemDamage%Grad_BF(iDoF1,iGauss) * xDof(iDoF1)
         End Do
         C0Operator = k * (damageGauss - 1.0_Kr) * elasticEnergyDensityGauss / &
                     (k + (1.0_Kr - k) * (1.0_kr - damageGauss)**2)**2

         Do iDoF2 = 1,numDofDamage
            residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * ( &
                                 (C0Operator + C1 * (1.0_Kr - damageGauss) &
                                 - ( N * cumulatedDissipatedPlasticEnergyCell * (1.0_Kr - damageGauss)**(DBLE(N - 1.0))  ) ) * elemDamage%BF(iDoF2,iGauss) & 
                                  + C2 * (gradientDamageGauss .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDamageLinSoftLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamageAT1ElasticLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDamageAT1ElasticLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDamageAT1ElasticLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      PetscReal,Intent(IN)                               :: CrackPressureCell

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDamage,numGauss
      Type(MEF90_VECT)                                   :: gradientDamageGauss
      PetscReal                                          :: C1,C2
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)
      
      C1 = matprop%fractureToughness / matprop%internalLength * .375
      C2 = matprop%fractureToughness * matprop%internalLength * .75
      residualLoc = 0.0_Kr
      Do iGauss = 1,numGauss
         gradientDamageGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            gradientDamageGauss = gradientDamageGauss + elemDamage%Grad_BF(iDoF1,iGauss) * xDof(iDoF1)
         End Do

         Do iDoF2 = 1,numDofDamage
            residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * ( &
                                  C1 * elemDamage%BF(iDoF2,iGauss) + &
                                  C2 * ((matprop%toughnessAnisotropyMatrix * gradientDamageGauss) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDamageAT1ElasticLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamageKKLLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDamageKKLLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDamageKKLLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      PetscReal,Intent(IN)                               :: CrackPressureCell

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: elasticEnergyDensityGauss,temperatureGauss,damageGauss,Stress_ZZ_planeStrain
      Type(MEF90_MATS)                                   :: inelasticStrainGauss
      Type(MEF90_VECT)                                   :: gradientDamageGauss,Displacement
      PetscReal                                          :: C1,C2,DgDamageGauss,N
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)

      N  = matprop%DuctileCouplingPower
      C1 = matprop%fractureToughness / matprop%internalLength / 4.0 / cwKKL
      C2 = matprop%fractureToughness * matprop%internalLength * 2.0 / cwKKL

      residualLoc = 0.0_Kr
      Do iGauss = 1,numGauss
         temperatureGauss = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
            End Do
         End If

         Displacement = 0.0_Kr
         inelasticStrainGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDisplacement
            Displacement = Displacement + elemDisplacement%BF(iDoF1,iGauss) * displacementDof(iDoF1)
            inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
         End Do
         inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion 
         elasticEnergyDensityGauss = (matprop%HookesLaw * (inelasticStrainGauss - plasticStrainCell) ) .DotP. (inelasticStrainGauss - plasticStrainCell)

#if MEF90_DIM == 2
         !!! Adding terms in planestrain for plasticity with tr(p) = 0
         If (.NOT. matProp%HookesLaw%isPlaneStress) Then
           Stress_ZZ_planeStrain = ( matprop%HookesLaw%YoungsModulus - 2.0_Kr*matprop%HookesLaw%PoissonRatio*matprop%HookesLaw%mu )*trace(plasticStrainCell) + matprop%HookesLaw%lambda*trace(inelasticStrainGauss)
           elasticEnergyDensityGauss = elasticEnergyDensityGauss + ( Stress_ZZ_planeStrain + matprop%HookesLaw%lambda*trace(inelasticStrainGauss - plasticStrainCell) ) * trace(plasticStrainCell)
         End If
#endif
         elasticEnergyDensityGauss = elasticEnergyDensityGauss * 0.5_Kr

         damageGauss = 0.0_Kr
         gradientDamageGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
            gradientDamageGauss = gradientDamageGauss + elemDamage%Grad_BF(iDoF1,iGauss) * xDof(iDoF1)
         End Do
         DgDamageGauss = MEF90DgKKL(damageGauss)

         If (N==1.0_Kr) Then
            Do iDoF2 = 1,numDofDamage
               residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * ( &
                                    ( (elasticEnergyDensityGauss - C1) * DgDamageGauss &
                                    - cumulatedDissipatedPlasticEnergyCell ) * elemDamage%BF(iDoF2,iGauss) + &
                                     (( C2 * matprop%toughnessAnisotropyMatrix * gradientDamageGauss +  Displacement*CrackPressureCell ) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
            End Do
         Else
            Do iDoF2 = 1,numDofDamage
               residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * ( &
                                    ( (elasticEnergyDensityGauss - C1) * DgDamageGauss &
                                    - N * cumulatedDissipatedPlasticEnergyCell * (1.0_Kr - damageGauss)**(N - 1.0_Kr) ) * elemDamage%BF(iDoF2,iGauss) + &
                                     (( C2 * matprop%toughnessAnisotropyMatrix * gradientDamageGauss +  Displacement*CrackPressureCell ) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
            End Do
         EndIf
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDamageKKLLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamageKKLElasticLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDamageKKLElasticLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDamageKKLElasticLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      PetscReal,Intent(IN)                               :: CrackPressureCell

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDamage,numGauss
      PetscReal                                          :: damageGauss
      Type(MEF90_VECT)                                   :: gradientDamageGauss
      PetscReal                                          :: C1,C2,DgDamageGauss
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)

      C1 = matprop%fractureToughness / matprop%internalLength / 4.0 / cwKKL
      C2 = matprop%fractureToughness * matprop%internalLength * 2.0 / cwKKL

      residualLoc = 0.0_Kr
      Do iGauss = 1,numGauss
         damageGauss = 0.0_Kr
         gradientDamageGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
            gradientDamageGauss = gradientDamageGauss + elemDamage%Grad_BF(iDoF1,iGauss) * xDof(iDoF1)
         End Do
         DgDamageGauss = MEF90DgKKL(damageGauss)

         Do iDoF2 = 1,numDofDamage
            residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * ( &
                                  - C1 * DgDamageGauss * elemDamage%BF(iDoF2,iGauss) + &
                                  (C2 * matprop%toughnessAnisotropyMatrix * gradientDamageGauss .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDamageKKLElasticLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamageAT1UnilateralHDLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDamageAT1UnilateralHDLoc:
!!!  
!!!  (c) 2014 - 2018 Blaise Bourdin bourdin@lsu.edu, Erwan Tanne erwan.tanne@gmail.com
!!!

   Subroutine MEF90DefMechOperatorDamageAT1UnilateralHDLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      PetscReal,Intent(IN)                               :: CrackPressureCell

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: elasticEnergyDensityGauss,temperatureGauss,damageGauss,Stress_ZZ_planeStrain
      Type(MEF90_MATS)                                   :: inelasticStrainGauss
      Type(MEF90_VECT)                                   :: gradientDamageGauss,DisplacementGauss
      PetscReal                                          :: C1,C2,N,kappa
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr


      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)
      kappa = ((matProp%HookesLaw * MEF90_MATS_IDENTITY) .dotP. MEF90_MATS_IDENTITY) / MEF90_DIM**2
      
      N  = matprop%DuctileCouplingPower
      C1 = matprop%fractureToughness / matprop%internalLength * .375
      C2 = matprop%fractureToughness * matprop%internalLength * .75

      residualLoc = 0.0_Kr
      Do iGauss = 1,numGauss
         temperatureGauss = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
            End Do
         End If

         DisplacementGauss = 0.0_Kr
         inelasticStrainGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDisplacement
            DisplacementGauss = DisplacementGauss + elemDisplacement%BF(iDoF1,iGauss) * displacementDof(iDoF1)
            inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
         End Do
         inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion 
         elasticEnergyDensityGauss = ((matprop%HookesLaw * DeviatoricPart(inelasticStrainGauss - plasticStrainCell)) .DotP. (inelasticStrainGauss - plasticStrainCell)) + &
                                      kappa * MEF90HydrostaticDeviatoricPenaltyFunction(Trace(inelasticStrainGauss - plasticStrainCell))
         !elasticEnergyDensityGauss = (matprop%HookesLaw * (inelasticStrainGauss - plasticStrainCell)) .DotP. (inelasticStrainGauss - plasticStrainCell)
         !!! This is really twice the elastic energy density

#if MEF90_DIM == 2
         !!! Adding terms in planestrain for plasticity with tr(p) = 0
         !!! Not completely sure about this in the context of HD decomposition
         !If (.NOT. matProp%HookesLaw%isPlaneStress) Then
         !  Stress_ZZ_planeStrain = ( matprop%HookesLaw%YoungsModulus - 2.0_Kr*matprop%HookesLaw%PoissonRatio*matprop%HookesLaw%mu )*trace(plasticStrainCell) + matprop%HookesLaw%lambda*trace(inelasticStrainGauss)
         !  elasticEnergyDensityGauss = elasticEnergyDensityGauss + ( Stress_ZZ_planeStrain + matprop%HookesLaw%lambda*trace(inelasticStrainGauss - plasticStrainCell) ) * trace(plasticStrainCell)
         !End If
#endif
         
         damageGauss = 0.0_Kr
         gradientDamageGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
            gradientDamageGauss = gradientDamageGauss + elemDamage%Grad_BF(iDoF1,iGauss) * xDof(iDoF1)
         End Do

         If (N==1.0_Kr) Then
            Do iDoF2 = 1,numDofDamage
               residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * ( &
                                    ( elasticEnergyDensityGauss * (damageGauss - 1.0_Kr) &
                                    - ( cumulatedDissipatedPlasticEnergyCell ) + C1 ) * elemDamage%BF(iDoF2,iGauss) + &
                                     ( ( C2 * gradientDamageGauss +  DisplacementGauss * CrackPressureCell ) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
            End Do
         Else
            Do iDoF2 = 1,numDofDamage
               residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * ( &
                                    ( elasticEnergyDensityGauss * (damageGauss - 1.0_Kr) &
                                    - ( N * cumulatedDissipatedPlasticEnergyCell * (1.0_Kr - damageGauss)**(N - 1.0_Kr)  ) + C1 ) * elemDamage%BF(iDoF2,iGauss) + &
                                     ( ( C2 * gradientDamageGauss +  DisplacementGauss * CrackPressureCell ) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
            End Do
         EndIf
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDamageAT1UnilateralHDLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamageAT1UnilateralDeviatoricLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDamageAT1UnilateralDeviatoricLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDamageAT1UnilateralDeviatoricLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      PetscReal,Intent(IN)                               :: CrackPressureCell

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: elasticEnergyDensityGauss,temperatureGauss,damageGauss
      Type(MEF90_MATS)                                   :: inelasticStrainGauss
      Type(MEF90_VECT)                                   :: gradientDamageGauss
      PetscReal                                          :: C1,C2
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)
      
      C1 = matprop%fractureToughness / matprop%internalLength * .375
      C2 = matprop%fractureToughness * matprop%internalLength * .75
      residualLoc = 0.0_Kr
      Do iGauss = 1,numGauss
         temperatureGauss = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
            End Do
         End If

         inelasticStrainGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
         End Do
         inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion - plasticStrainCell

         elasticEnergyDensityGauss = (matprop%HookesLaw * inelasticStrainGauss) .DotP. DeviatoricPart(inelasticStrainGauss)
         !!! This is really twice the elastic energy density

         damageGauss = 0.0_Kr
         gradientDamageGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
            gradientDamageGauss = gradientDamageGauss + elemDamage%Grad_BF(iDoF1,iGauss) * xDof(iDoF1)
         End Do

         Do iDoF2 = 1,numDofDamage
            residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * ( &
                                  (elasticEnergyDensityGauss * (damageGauss - 1.0_Kr) + C1) * elemDamage%BF(iDoF2,iGauss) + &
                                  C2 * (gradientDamageGauss .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDamageAT1UnilateralDeviatoricLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamageAT1UnilateralPHDLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDamageAT1UnilateralPHDLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDamageAT1UnilateralPHDLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      PetscReal,Intent(IN)                               :: CrackPressureCell

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: elasticEnergyDensityGauss,temperatureGauss,damageGauss
      Type(MEF90_MATS)                                   :: inelasticStrainGauss
      Type(MEF90_VECT)                                   :: gradientDamageGauss
      PetscReal                                          :: C1,C2
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)
      
      C1 = matprop%fractureToughness / matprop%internalLength * .375
      C2 = matprop%fractureToughness * matprop%internalLength * .75
      residualLoc = 0.0_Kr
      Do iGauss = 1,numGauss
         temperatureGauss = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
            End Do
         End If

         inelasticStrainGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
         End Do
         inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion - plasticStrainCell

         If (Trace(inelasticStrainGauss) < 0.0_Kr) Then
            elasticEnergyDensityGauss = 0.0_Kr
         Else
            elasticEnergyDensityGauss = (matprop%HookesLaw * inelasticStrainGauss) .DotP. inelasticStrainGauss
         End If
         !!! This is really twice the elastic energy density

         damageGauss = 0.0_Kr
         gradientDamageGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
            gradientDamageGauss = gradientDamageGauss + elemDamage%Grad_BF(iDoF1,iGauss) * xDof(iDoF1)
         End Do

         Do iDoF2 = 1,numDofDamage
            residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * ( &
                                  (elasticEnergyDensityGauss * (damageGauss - 1.0_Kr) + C1) * elemDamage%BF(iDoF2,iGauss) + &
                                  C2 * (gradientDamageGauss .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDamageAT1UnilateralPHDLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamageAT1UnilateralMasonryLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDamageAT1UnilateralMasonryLoc:
!!!  
!!!  (c) 2016 Erwan TANNE erwan.tanne@gmail.com, Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDamageAT1UnilateralMasonryLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      PetscReal,Intent(IN)                               :: CrackPressureCell

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: positiveElasticEnergyDensityGauss
      Type(MEF90_MATS)                                   :: inelasticStrain,positiveInelasticStrain,negativeInelasticStrain
      Type(MEF90_VECT)                                   :: gradientDamageGauss
      PetscReal                                          :: damageGauss
      PetscReal                                          :: C1,C2,N
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)
      
      N  = matprop%DuctileCouplingPower
      C1 = matprop%fractureToughness / matprop%internalLength * .375
      C2 = matprop%fractureToughness * matprop%internalLength * .75

      residualLoc = 0.0_Kr
      Do iGauss = 1,numGauss
         inelasticStrain = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               inelasticStrain = inelasticStrain - temperatureDof(iDoF1) * (matProp%LinearThermalExpansion * elemDamage%BF(iDoF1,iGauss))
            End Do         
         End If
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrain = inelasticStrain + displacementDof(iDoF1) * elemDisplacement%GradS_BF(iDoF1,iGauss)
         End Do
         inelasticStrain = inelasticStrain - plasticStrainCell
         Call MasonryProjection(inelasticStrain,matprop%HookesLaw,positiveInelasticStrain,negativeInelasticStrain)
         positiveElasticEnergyDensityGauss = (matprop%HookesLaw * positiveInelasticStrain) .dotP. positiveInelasticStrain
         !!! This is really twice the elastic energy density

         damageGauss = 0.0_Kr
         gradientDamageGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
            gradientDamageGauss = gradientDamageGauss + elemDamage%Grad_BF(iDoF1,iGauss) * xDof(iDoF1)
         End Do

         Do iDoF2 = 1,numDofDamage
            residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * ( &
                                 ( positiveElasticEnergyDensityGauss * (damageGauss - 1.0_Kr) &
                                 - ( N * cumulatedDissipatedPlasticEnergyCell * (1.0_Kr - damageGauss)**(N - 1)  ) + C1 ) * elemDamage%BF(iDoF2,iGauss) + &
                                 C2 * (gradientDamageGauss .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDamageAT1UnilateralMasonryLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamageAT2Loc"
!!!
!!!  
!!!  MEF90DefMechOperatorDamageAT2Loc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDamageAT2Loc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      PetscReal,Intent(IN)                               :: CrackPressureCell

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: elasticEnergyDensityGauss,temperatureGauss,damageGauss
      Type(MEF90_MATS)                                   :: inelasticStrainGauss
      Type(MEF90_VECT)                                   :: gradientDamageGauss
      PetscReal                                          :: C1,C2
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)
      
      C1 = matprop%fractureToughness / matprop%internalLength
      C2 = matprop%fractureToughness * matprop%internalLength
      residualLoc = 0.0_Kr
      Do iGauss = 1,numGauss
         temperatureGauss = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
            End Do
         End If

         inelasticStrainGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
         End Do

         inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion - plasticStrainCell
         elasticEnergyDensityGauss = (matprop%HookesLaw * inelasticStrainGauss) .DotP. inelasticStrainGauss
         !!! This is really twice the elastic energy density

         damageGauss = 0.0_Kr
         gradientDamageGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
            gradientDamageGauss = gradientDamageGauss + elemDamage%Grad_BF(iDoF1,iGauss) * xDof(iDoF1)
         End Do

         Do iDoF2 = 1,numDofDamage
            residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * ( &
                                  (elasticEnergyDensityGauss * (damageGauss - 1.0_Kr)+ C1 * damageGauss) * elemDamage%BF(iDoF2,iGauss) + &
                                  C2 * ((matprop%toughnessAnisotropyMatrix * gradientDamageGauss) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDamageAT2Loc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamageAT2ElasticLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDamageAT2ElasticLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDamageAT2ElasticLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      PetscReal,Intent(IN)                               :: CrackPressureCell

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDamage,numGauss
      PetscReal                                          :: damageGauss
      Type(MEF90_VECT)                                   :: gradientDamageGauss
      PetscReal                                          :: C1,C2
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)
      
      C1 = matprop%fractureToughness / matprop%internalLength
      C2 = matprop%fractureToughness * matprop%internalLength
      residualLoc = 0.0_Kr
      Do iGauss = 1,numGauss
         damageGauss = 0.0_Kr
         gradientDamageGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) *  xDof(iDoF1)
            gradientDamageGauss = gradientDamageGauss + elemDamage%Grad_BF(iDoF1,iGauss) * xDof(iDoF1)
         End Do

         Do iDoF2 = 1,numDofDamage
            residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * ( &
                                  C1 * damageGauss * elemDamage%BF(iDoF2,iGauss) + &
                                  C2 * ((matprop%toughnessAnisotropyMatrix * gradientDamageGauss) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDamageAT2ElasticLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamageAT2UnilateralHDLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDamageAT2UnilateralHDLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDamageAT2UnilateralHDLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      PetscReal,Intent(IN)                               :: CrackPressureCell

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: elasticEnergyDensityGauss,temperatureGauss,damageGauss
      Type(MEF90_MATS)                                   :: inelasticStrainGauss
      Type(MEF90_VECT)                                   :: gradientDamageGauss,displacementGauss
      PetscReal                                          :: C1,C2,N,kappa
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)
      kappa = ((matProp%HookesLaw * MEF90_MATS_IDENTITY) .dotP. MEF90_MATS_IDENTITY) / MEF90_DIM**2
 
      N  = matprop%DuctileCouplingPower
      C1 = matprop%fractureToughness / matprop%internalLength
      C2 = matprop%fractureToughness * matprop%internalLength
      residualLoc = 0.0_Kr
      Do iGauss = 1,numGauss
         temperatureGauss = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
            End Do
         End If

         DisplacementGauss = 0.0_Kr
         inelasticStrainGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDisplacement
            DisplacementGauss = DisplacementGauss + elemDisplacement%BF(iDoF1,iGauss) * displacementDof(iDoF1)
            inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
         End Do
         inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion - plasticStrainCell
         elasticEnergyDensityGauss = ((matprop%HookesLaw * DeviatoricPart(inelasticStrainGauss - plasticStrainCell)) .DotP. (inelasticStrainGauss - plasticStrainCell)) + &
                                      kappa * MEF90HydrostaticDeviatoricPenaltyFunction(Trace(inelasticStrainGauss - plasticStrainCell))
         !!! This is really twice the elastic energy density

         damageGauss = 0.0_Kr
         gradientDamageGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
            gradientDamageGauss = gradientDamageGauss + elemDamage%Grad_BF(iDoF1,iGauss) * xDof(iDoF1)
         End Do

         If ((N == 1.0_Kr) .OR. (cumulatedDissipatedPlasticEnergyCell == 0.0_Kr)) Then
            Do iDoF2 = 1,numDofDamage
               residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * (elemDamage%BF(iDoF2,iGauss) *  &
                                                (elasticEnergyDensityGauss * (damageGauss - 1.0_Kr) - cumulatedDissipatedPlasticEnergyCell + C1 * damageGauss) &
                                             + ((C2 * gradientDamageGauss +  DisplacementGauss * CrackPressureCell) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
            End Do
         Else
            Do iDoF2 = 1,numDofDamage
               residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * ( &
                                    ( elasticEnergyDensityGauss * (damageGauss - 1.0_Kr) &
                                    - ( N * cumulatedDissipatedPlasticEnergyCell * (1.0_Kr - damageGauss)**(N - 1.0_Kr)  ) + C1 * damageGauss) * elemDamage%BF(iDoF2,iGauss) + &
                                     ( ( C2 * gradientDamageGauss +  DisplacementGauss * CrackPressureCell ) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
            End Do
         EndIf
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDamageAT2UnilateralHDLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamageAT2UnilateralDeviatoricLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDamageAT2UnilateralDeviatoricLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDamageAT2UnilateralDeviatoricLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      PetscReal,Intent(IN)                               :: CrackPressureCell

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: elasticEnergyDensityGauss,temperatureGauss,damageGauss
      Type(MEF90_MATS)                                   :: inelasticStrainGauss
      Type(MEF90_VECT)                                   :: gradientDamageGauss
      PetscReal                                          :: C1,C2
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)
      
      C1 = matprop%fractureToughness / matprop%internalLength
      C2 = matprop%fractureToughness * matprop%internalLength
      residualLoc = 0.0_Kr
      Do iGauss = 1,numGauss
         temperatureGauss = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
            End Do
         End If

         inelasticStrainGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
         End Do
         inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion - plasticStrainCell

         elasticEnergyDensityGauss = (matprop%HookesLaw * inelasticStrainGauss) .DotP. DeviatoricPart(inelasticStrainGauss)
         !!! This is really twice the elastic energy density

         damageGauss = 0.0_Kr
         gradientDamageGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
            gradientDamageGauss = gradientDamageGauss + elemDamage%Grad_BF(iDoF1,iGauss) * xDof(iDoF1)
         End Do

         Do iDoF2 = 1,numDofDamage
            residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * ( &
                                  (elasticEnergyDensityGauss * (damageGauss - 1.0_Kr)+ C1 * damageGauss) * elemDamage%BF(iDoF2,iGauss) + &
                                  C2 * (gradientDamageGauss .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDamageAT2UnilateralDeviatoricLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamageAT2UnilateralPHDLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDamageAT2UnilateralPHDLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDamageAT2UnilateralPHDLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      PetscReal,Intent(IN)                               :: CrackPressureCell

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: elasticEnergyDensityGauss,temperatureGauss,damageGauss
      Type(MEF90_MATS)                                   :: inelasticStrainGauss
      Type(MEF90_VECT)                                   :: gradientDamageGauss
      PetscReal                                          :: C1,C2
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)
      
      C1 = matprop%fractureToughness / matprop%internalLength
      C2 = matprop%fractureToughness * matprop%internalLength
      residualLoc = 0.0_Kr
      Do iGauss = 1,numGauss
         temperatureGauss = 0.0_Kr
         If (Associated(temperatureDof)) Then
            Do iDoF1 = 1,numDofDamage
               temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
            End Do
         End If

         inelasticStrainGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDisplacement
            inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
         End Do
         inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion - plasticStrainCell

         If (Trace(inelasticStrainGauss) < 0.0_Kr) Then
            elasticEnergyDensityGauss = 0.0_Kr
         Else
            elasticEnergyDensityGauss = (matprop%HookesLaw * inelasticStrainGauss) .DotP. inelasticStrainGauss
         End If
         !!! This is really twice the elastic energy density

         damageGauss = 0.0_Kr
         gradientDamageGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
            gradientDamageGauss = gradientDamageGauss + elemDamage%Grad_BF(iDoF1,iGauss) * xDof(iDoF1)
         End Do

         Do iDoF2 = 1,numDofDamage
            residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * ( &
                                  (elasticEnergyDensityGauss * (damageGauss - 1.0_Kr)+ C1 * damageGauss) * elemDamage%BF(iDoF2,iGauss) + &
                                  C2 * (gradientDamageGauss .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDamageAT2UnilateralPHDLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamageDampingLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDamageDampingLoc: Contribution of the damping term \int delta (u-u_{n-1}) 
!!!                                              in an element
!!!  
!!!  (c) 2018 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDamageDampingLoc(residualLoc,xDof,xpreviousStepDof,delta,elemDamage)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,xpreviousStepDof
      PetscReal                                          :: delta
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      
      PetscReal                                          :: contrElem
      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDamage,numGauss
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)
      Do iGauss = 1,numGauss
         contrElem = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            contrElem = contrElem + (xDof(iDof1) - xpreviousStepDof(iDoF1)) * elemDamage%BF(iDoF1,iGauss) 
         End Do
         Do iDoF2 = 1,numDofDamage
            residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * &
                                 delta * ( contrElem * elemDamage%BF(iDoF2,iGauss) )
         End Do
      End Do
      flops = 7 * numDofDamage * numGauss
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDamageDampingLoc


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamage"
!!!
!!!  
!!!  MEF90DefMechOperatorDamage: Build the operator. When called in SNES, the solution time should always match the target time, 
!!!                                    so there is no need for interpolation of the forcees, external, and boundary values
!!!  
!!!  (c) 2012-16 Blaise Bourdin bourdin@lsu.edu, Erwan Tanne erwan.tanne@gmail.com
!!!

   Subroutine MEF90DefMechOperatorDamage(snesDamage,damage,residual,MEF90DefMechCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesDamage
      Type(Vec),Intent(IN)                               :: damage
      Type(Vec),Intent(INOUT)                            :: residual
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(DM)                                           :: mesh
      Type(SectionReal)                                  :: displacementSec,residualSec
      Type(SectionReal)                                  :: plasticStrainSec,temperatureSec,damageSec,damagePreviousStepSec
      Type(SectionReal)                                  :: boundaryDamageSec,forceSec,pressureForceSec
      Type(SectionReal)                                  :: CrackPressureSec,cumulatedDissipatedPlasticEnergySec
      PetscReal,Dimension(:),Pointer                     :: damagePtr,residualPtr
      PetscReal,Dimension(:),Pointer                     :: boundaryDamagePtr
      PetscReal,Dimension(:),Pointer                     :: displacementDof,damageDof,damagePreviousStepDof,temperatureDof,residualLoc
      PetscReal,Dimension(:),Pointer                     :: plasticStrainLoc
      PetscReal,Dimension(:),Pointer                     :: CrackPressureLoc,cumulatedDissipatedPlasticEnergyLoc
      Type(MEF90_MATS)                                   :: plasticStrainCell
      Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof,bcIS
      PetscInt,Dimension(:),Pointer                      :: setID,cellID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder,cell
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90DefMechVertexSetOptions_Type),Pointer    :: vertexSetOptions
      PetscReal                                          :: CrackPressureCell,cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemDamage
      Type(MEF90Element_Type)                            :: elemDisplacementType,elemDamageType
      PetscInt                                           :: p,dof,c,dim,numDof
      PetscInt                                           :: nVal
      PetscReal,Dimension(:),Pointer                     :: nullPtr
      
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: globalOptions
      Procedure(MEF90DefMechOperatorLoc),pointer         :: localOperatorFunction
      
      localOperatorFunction =>MEF90DefMechOperatorNull
      Call SNESGetDM(snesDamage,mesh,ierr);CHKERRQ(ierr)

      !!! Get Section and scatter associated with each field


      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,displacementSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(displacementSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,MEF90DefMechCtx%displacement,ierr);CHKERRQ(ierr) 

      If (globalOptions%dampingCoefficientDamage * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,damagepreviousStepSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(damagepreviousStepSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%damagepreviousStep,ierr);CHKERRQ(ierr) 
      EndIf
      
      Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,damageSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(damageSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,damage,ierr);CHKERRQ(ierr) 

      Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,boundaryDamageSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(boundaryDamageSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%boundaryDamage,ierr);CHKERRQ(ierr) 


      Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,residualSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(residualSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,residual,ierr);CHKERRQ(ierr) 
      Call SectionRealSet(residualSec,0.0_Kr,ierr);CHKERRQ(ierr)
      Call VecSet(residual,0.0_Kr,ierr);CHKERRQ(ierr)

      Nullify(TemperatureDof)
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,temperatureSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(temperatureSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr) 
      Else
         temperatureSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%CrackPressure)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%cellDMScalSec,CrackPressureSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(CrackPressureSec,MEF90DefMechCtx%CellDMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%CrackPressure,ierr);CHKERRQ(ierr)
      Else
         CrackPressureSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%cellDMMatSSec,plasticstrainSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticstrainSec,MEF90DefMechCtx%CellDMMatSScatter,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr) 
      
         Call SectionRealDuplicate(MEF90DefMechCtx%cellDMScalSec,cumulatedDissipatedPlasticEnergySec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(cumulatedDissipatedPlasticEnergySec,MEF90DefMechCtx%CellDMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%cumulatedDissipatedPlasticEnergy,ierr);CHKERRQ(ierr) 
      Else
         PlasticStrainSec%v = 0
         cumulatedDissipatedPlasticEnergySec%v = 0
      End If

      !!! get IS for cell sets
      Call DMMeshGetDimension(mesh,dim,ierr);CHKERRQ(ierr)
      Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      !!!
      !!! We loop over all element twice. The first time in order to assembly all non BC cell sets
      !!! In the second pass, we only update the BC where necessary
      !!! vertex set BC are updated last, so that they override cell set BC
      !!!
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
         elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
         elemDamageType       = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)
         If ((Size(cellID) > 0) .AND. (elemDamageType%coDim == 0)) Then
            Select Case (cellSetOptions%damageType)
            Case (MEF90DefMech_damageTypeAT1Elastic)
               QuadratureOrder = 2 * (elemDamageType%order - 1)
               localOperatorFunction => MEF90DefMechOperatorDamageAT1ElasticLoc

            Case (MEF90DefMech_damageTypeAT2Elastic)
               QuadratureOrder = 2 * elemDamageType%order
               localOperatorFunction => MEF90DefMechOperatorDamageAT2ElasticLoc

            Case (MEF90DefMech_damageTypeKKLElastic)
               QuadratureOrder = max(3, 2 * (elemDamageType%order-1))
               localOperatorFunction => MEF90DefMechOperatorDamageKKLElasticLoc

            Case (MEF90DefMech_damageTypeAT1)
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  QuadratureOrder = 2 * elemDamageType%order + 2 * max (elemDisplacementType%order - 1, elemDamageType%order)
               Else
                  QuadratureOrder = 2 * elemDamageType%order + 2 * (elemDisplacementType%order - 1)
               End If
               Select Case(cellSetOptions%unilateralContactType)
               Case (MEF90DefMech_unilateralContactTypeNone,MEF90DefMech_unilateralContactTypeBrittleDuctile)
                  localOperatorFunction => MEF90DefMechOperatorDamageAT1Loc
               Case (MEF90DefMech_unilateralContactTypeHydrostaticDeviatoric,MEF90DefMech_unilateralContactTypeHybridHydrostaticDeviatoric)
                  localOperatorFunction => MEF90DefMechOperatorDamageAT1UnilateralHDLoc
               Case (MEF90DefMech_unilateralContactTypePositiveHydrostatic)
                  localOperatorFunction => MEF90DefMechOperatorDamageAT1UnilateralPHDLoc
               Case (MEF90DefMech_unilateralContactTypeDeviatoric)
                  localOperatorFunction => MEF90DefMechOperatorDamageAT1UnilateralDeviatoricLoc
               Case (MEF90DefMech_unilateralContactTypePrincipalStrains)
                  localOperatorFunction => MEF90DefMechOperatorNull
               Case (MEF90DefMech_unilateralContactTypeMasonry)
                  localOperatorFunction => MEF90DefMechOperatorDamageAT1UnilateralMasonryLoc
               End Select

            Case (MEF90DefMech_damageTypeLinSoft)
               Select Case(cellSetOptions%unilateralContactType)
               Case (MEF90DefMech_unilateralContactTypeNone)
                  QuadratureOrder = 9! 5 * elemDamageType%order + 5 * (elemDisplacementType%order - 1)
                  localOperatorFunction => MEF90DefMechOperatorDamageLinSoftLoc
               End Select

            Case (MEF90DefMech_damageTypeAT2)
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  QuadratureOrder = 2 * elemDamageType%order + 2 * max (elemDisplacementType%order - 1, elemDamageType%order)
               Else
                  QuadratureOrder = 2 * elemDamageType%order + 2 * (elemDisplacementType%order - 1)
               End If
               Select Case(cellSetOptions%unilateralContactType)
               Case (MEF90DefMech_unilateralContactTypeNone)
                  localOperatorFunction => MEF90DefMechOperatorDamageAT2Loc
               Case (MEF90DefMech_unilateralContactTypeHydrostaticDeviatoric,MEF90DefMech_unilateralContactTypeHybridHydrostaticDeviatoric)
                  localOperatorFunction => MEF90DefMechOperatorDamageAT2UnilateralHDLoc
               Case (MEF90DefMech_unilateralContactTypePositiveHydrostatic)
                  localOperatorFunction => MEF90DefMechOperatorDamageAT2UnilateralPHDLoc
               Case (MEF90DefMech_unilateralContactTypeDeviatoric)
                  localOperatorFunction => MEF90DefMechOperatorDamageAT2UnilateralDeviatoricLoc
               Case (MEF90DefMech_unilateralContactTypePrincipalStrains,MEF90DefMech_unilateralContactTypeMasonry)
                  localOperatorFunction => MEF90DefMechOperatorNull
               End Select

            Case (MEF90DefMech_damageTypeKKL)
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  QuadratureOrder = 3 + 2 * max (elemDisplacementType%order - 1, elemDamageType%order)
               Else
                  QuadratureOrder = 3 + 2 * (elemDisplacementType%order - 1)
               End If
               Select Case(cellSetOptions%unilateralContactType)
               Case (MEF90DefMech_unilateralContactTypeNone)
                  localOperatorFunction => MEF90DefMechOperatorDamageKKLLoc
               End Select
            End Select

            !!! Allocate storage for fields at dof and Gauss points
            !!! Leaving plasticstrain aside until I can remember how it is supposed to be dealt with
            Allocate(residualLoc(ElemDamageType%numDof))
            Allocate(displacementDof(ElemDisplacementType%numDof))
            Allocate(damageDof(ElemDamageType%numDof))
            If (globalOptions%dampingCoefficientDamage * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
               Allocate(damagePreviousStepDof(ElemDamageType%numDof))
            Else
               Nullify(damagepreviousStepDof)
            End If
            If (Associated(MEF90DefMechCtx%temperature)) Then
               Allocate(temperatureDof(ElemDamageType%numDof))
            End If

            !!! Allocate elements 
            Call MEF90Element_Create(mesh,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
            Call MEF90Element_Create(mesh,setIS,elemDamage,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
            Do cell = 1,size(cellID)
               !!! Get value of each field at each Dof of the local element
               Call SectionRealRestrictClosure(displacementSec,MEF90DefMechCtx%DMVect,cellID(cell),elemDisplacementType%numDof,displacementDof,ierr);CHKERRQ(ierr)
               Call SectionRealRestrictClosure(damageSec,MEF90DefMechCtx%DMScal,cellID(cell),elemDamageType%numDof,damageDof,ierr);CHKERRQ(ierr)
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  Call SectionRealRestrictClosure(temperatureSec,MEF90DefMechCtx%DMScal,cellID(cell),elemDamageType%numDof,temperatureDof,ierr);CHKERRQ(ierr)
               End If
               If (Associated(MEF90DefMechCtx%CrackPressure)) Then
                  Call SectionRealRestrict(CrackPressureSec,cellID(cell),CrackPressureLoc,ierr);CHKERRQ(ierr)
                  CrackPressureCell = CrackPressureLoc(1)
               Else
                  CrackPressureCell = 0.0_Kr
               End If
               If (Associated(MEF90DefMechCtx%plasticStrain)) Then
                  Call SectionRealRestrict(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
                  Call SectionRealRestrict(cumulatedDissipatedPlasticEnergySec,cellID(cell),cumulatedDissipatedPlasticEnergyLoc,ierr);CHKERRQ(ierr)
                  plasticStrainCell = plasticStrainLoc
                  cumulatedDissipatedPlasticEnergyCell = cumulatedDissipatedPlasticEnergyLoc(1)
               Else
                  plasticStrainCell = 0.0_Kr
                  cumulatedDissipatedPlasticEnergyCell = 0.0_Kr
               End If
                  
               residualLoc = 0.0_Kr
               Call localOperatorFunction   (residualLoc,damageDof,displacementDof,nullPtr,nullPtr,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matpropSet,elemDisplacement(cell),elemDamage(cell),CrackPressureCell)
               If (globalOptions%dampingCoefficientDamage * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
                  Call SectionRealRestrictClosure(damagepreviousStepSec,MEF90DefMechCtx%DMScal,cellID(cell),elemDamageType%numDof,damagepreviousStepDof,ierr);CHKERRQ(ierr)
                  call MEF90DefMechOperatorDamageDampingLoc(residualLoc,damageDof,damagepreviousStepDof,globalOptions%dampingCoefficientDamage / MEF90DefMechCtx%timeStep,elemDamage(cell))
               End If
               Call SectionRealUpdateClosure(residualSec,MEF90DefMechCtx%DMScal,cellID(cell),residualLoc,ADD_VALUES,ierr);CHKERRQ(iErr)

               If (Associated(MEF90DefMechCtx%plasticStrain)) Then
                  Call SectionRealRestore(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
                  Call SectionRealRestore(cumulatedDissipatedPlasticEnergySec,cellID(cell),cumulatedDissipatedPlasticEnergyLoc,ierr);CHKERRQ(ierr)
               End If

               If (Associated(MEF90DefMechCtx%CrackPressure)) Then
                  Call SectionRealRestore(CrackPressureSec,cellID(cell),CrackPressureLoc,ierr);CHKERRQ(ierr)
               EndIf


            End Do
            DeAllocate(displacementDof)
            DeAllocate(damageDof)
            If (globalOptions%dampingCoefficientDamage * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
               DeAllocate(damagePreviousStepDof)
            End If
            If (Associated(MEF90DefMechCtx%temperature)) Then
               DeAllocate(temperatureDof)
            End If
            DeAllocate(residualLoc)

            Call MEF90Element_Destroy(elemDisplacement,ierr)
            Call MEF90Element_Destroy(elemDamage,ierr)
         End If 
         Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do ! set

      !!! "Ghost update" for the residual Section
      Call SectionRealComplete(residualSec,ierr);CHKERRQ(ierr)
      !!! Scatter back from SectionReal to Vec
      Call SectionRealToVec(residualSec,MEF90DefMechCtx%DMScalScatter,SCATTER_FORWARD,residual,ierr);CHKERRQ(ierr) 

      !!!
      !!! Cell set BC
      !!!
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         If (cellSetOptions%Has_DamageBC) Then
            Call DMMeshGetStratumIS(MEF90DefMechCtx%DMScal,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
            Call MEF90ISCreateCelltoVertex(MEF90DefMechCtx%DMScal,PETSC_COMM_WORLD,setIS,bcIS,ierr)
            Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
            Allocate(damagePtr(nval),stat=ierr)
            Allocate(residualPtr(nval),stat=ierr)
            Allocate(boundaryDamagePtr(nval),stat=ierr)
            Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMScal,damage,damagePtr,bcIS,1,ierr)
            Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMScal,MEF90DefMechCtx%boundaryDamage,boundaryDamagePtr,bcIS,1,ierr)
            residualPtr = damagePtr - boundaryDamagePtr
            Call MEF90VecSetValuesISdof(MEF90DefMechCtx%DMScal,residual,residualPtr,bcIS,1,INSERT_VALUES,ierr)
            DeAllocate(boundaryDamagePtr)
            DeAllocate(residualPtr)
            DeAllocate(damagePtr)
            Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If ! cellSetOptions%Has_BC
      End Do ! set
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
      
      !!!
      !!! Vertex set BC
      !!!
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxVertexSetOptions(MEF90DefMechCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
         If (vertexSetOptions%Has_DamageBC) Then
            Call DMMeshGetStratumIS(MEF90DefMechCtx%DMScal,'Vertex Sets',setID(set),bcIS,ierr);CHKERRQ(iErr)
            Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
            Allocate(damagePtr(nval),stat=ierr)
            Allocate(residualPtr(nval),stat=ierr)
            Allocate(boundaryDamagePtr(nval),stat=ierr)
            Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMScal,damage,damagePtr,bcIS,1,ierr)
            Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMScal,MEF90DefMechCtx%boundaryDamage,boundaryDamagePtr,bcIS,1,ierr)
            residualPtr = damagePtr - boundaryDamagePtr
            Call MEF90VecSetValuesISdof(MEF90DefMechCtx%DMScal,residual,residualPtr,bcIS,1,INSERT_VALUES,ierr)
            DeAllocate(boundaryDamagePtr)
            DeAllocate(residualPtr)
            DeAllocate(damagePtr)
            Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
         End If ! vertexSetOptions%Has_BC
      End Do ! set
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
   
      Call VecAssemblyBegin(residual,ierr)
      Call VecAssemblyEnd(residual,ierr)

      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
         Call SectionRealDestroy(cumulatedDissipatedPlasticEnergySec,ierr);CHKERRQ(ierr)
      End If

      If (Associated(MEF90DefMechCtx%CrackPressure)) Then
         Call SectionRealDestroy(CrackPressureSec,ierr);CHKERRQ(ierr)
      End If
      
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
      End If

      Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)
      If (globalOptions%dampingCoefficientDamage * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
         Call SectionRealDestroy(damagePreviousStepSec,ierr);CHKERRQ(ierr)
      End If
      Call SectionRealDestroy(boundaryDamageSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(residualSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(displacementSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDamage

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamage"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamage:
!!!  
!!!  (c) 2012-16 Blaise Bourdin bourdin@lsu.edu, Erwan Tanne erwan.tanne@gmail.com
!!!

   Subroutine MEF90DefMechBilinearFormDamage(snesDamage,damage,A,M,flg,MEF90DefMechCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesDamage
      Type(Vec),Intent(IN)                               :: damage
      Type(Mat),Intent(INOUT)                            :: A,M
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr  
      
      Type(SectionReal)                                  :: displacementSec,damageSec,temperatureSec,plasticStrainSec
      Type(SectionReal)                                  :: cumulatedDissipatedPlasticEnergySec
      PetscReal,Dimension(:),Pointer                     :: displacementDOF,damageDof,temperatureDof,nullPtr
      PetscReal,Dimension(:),Pointer                     :: plasticStrainLoc
      Type(MEF90_MATS)                                   :: plasticStrainCell
      PetscReal,Dimension(:),Pointer                     :: cumulatedDissipatedPlasticEnergyLoc
      PetscReal                                          :: cumulatedDissipatedPlasticEnergyCell
      Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90DefMechVertexSetOptions_Type),Pointer    :: vertexSetOptions
      Type(DM)                                           :: mesh
      Type(MEF90Element_Type)                            :: ElemDisplacementType,ElemDamageType
      PetscInt,Dimension(:),Pointer                      :: setIdx,bcIdx,Cone
      Type(IS)                                           :: bcIS
      PetscInt                                           :: cell,v,numBC,numDoF,numCell,c,dim
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscInt                                           :: iGauss,iDoF1,iDoF2
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemDamage
      
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: globalOptions
      Procedure(MEF90DefMechBilinearFormLoc),pointer     :: localAssemblyFunction
      
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)
      localAssemblyFunction =>MEF90DefMechBilinearFormNull
      Call MatZeroEntries(A,ierr);CHKERRQ(ierr)
      Call SNESGetDM(snesDamage,mesh,ierr);CHKERRQ(ierr)
      
      !!! Get Section and scatter associated with each field
      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,displacementSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(displacementSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,MEF90DefMechCtx%displacement,ierr);CHKERRQ(ierr)   

      Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,damageSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(damageSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,damage,ierr);CHKERRQ(ierr)   

      Nullify(temperatureDof)
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,temperatureSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(temperatureSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr)   
      Else
         temperatureSec%v = 0
      End If
      
      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%CellDMMatSSec,plasticStrainSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticStrainSec,MEF90DefMechCtx%CellDMMatSScatter,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)

         Call SectionRealDuplicate(MEF90DefMechCtx%CellDMScalSec,cumulatedDissipatedPlasticEnergySec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(cumulatedDissipatedPlasticEnergySec,MEF90DefMechCtx%CellDMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%cumulatedDissipatedPlasticEnergy,ierr);CHKERRQ(ierr)
      Else
         PlasticStrainSec%v = 0
         cumulatedDissipatedPlasticEnergySec%v = 0
      End If

      !!! get IS for cell sets
      Call DMMeshGetDimension(mesh,dim,ierr);CHKERRQ(ierr)
      Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID) 
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         ElemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
         ElemDamageType       = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)
         !!! Go through cells and call local assembly functions
         Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
         If ((Size(cellID) > 0) .AND. (ElemDisplacementType%coDim == 0)) Then
         !!! Call proper local assembly depending on the type of damage law
            !!! select the proper local assembly routine, compute proper integration order
            Select Case (cellSetOptions%damageType)
            Case (MEF90DefMech_damageTypeAT1Elastic)
               QuadratureOrder = max(elemDamageType%order, 2 * (elemDamageType%order - 1))
               localAssemblyFunction => MEF90DefMechBilinearFormDamageAT1ElasticLoc

            Case (MEF90DefMech_damageTypeAT2Elastic)
               QuadratureOrder = 2 * elemDamageType%order
               localAssemblyFunction => MEF90DefMechBilinearFormDamageAT2ElasticLoc

            Case (MEF90DefMech_damageTypeKKLElastic)
               QuadratureOrder = max(2, 2 * (elemDamageType%order-1))
               localAssemblyFunction => MEF90DefMechBilinearFormDamageKKLElasticLoc

            Case (MEF90DefMech_damageTypeLinSoft)
               QuadratureOrder = 9 !* (elemDisplacementType%order - 1) + 5 * elemDamageType%order
               localAssemblyFunction => MEF90DefMechBilinearFormDamageLinSoftLoc
               Select Case(cellSetOptions%unilateralContactType)
               Case (MEF90DefMech_unilateralContactTypeNone)
                  QuadratureOrder = 9 !* (elemDisplacementType%order - 1) + 5 * elemDamageType%order
                  localAssemblyFunction => MEF90DefMechBilinearFormDamageLinSoftLoc
               End Select

            Case (MEF90DefMech_damageTypeAT1)
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  QuadratureOrder = 2 * max(elemDisplacementType%order - 1, elemDamageType%order) + 2 * elemDamageType%order
               Else
                  QuadratureOrder = 2 * (elemDisplacementType%order - 1) + 2 * elemDamageType%order
               End If
               Select Case(cellSetOptions%unilateralContactType)
               Case (MEF90DefMech_unilateralContactTypeNone,MEF90DefMech_unilateralContactTypeBrittleDuctile)
                  localAssemblyFunction => MEF90DefMechBilinearFormDamageAT1Loc
               Case (MEF90DefMech_unilateralContactTypeHydrostaticDeviatoric,MEF90DefMech_unilateralContactTypeHybridHydrostaticDeviatoric)
                  localAssemblyFunction => MEF90DefMechBilinearFormDamageAT1UnilateralHDLoc
               Case (MEF90DefMech_unilateralContactTypePositiveHydrostatic)
                  localAssemblyFunction => MEF90DefMechBilinearFormDamageAT1UnilateralPHDLoc
               Case (MEF90DefMech_unilateralContactTypeDeviatoric)
                  localAssemblyFunction => MEF90DefMechBilinearFormDamageAT1UnilateralDeviatoricLoc
               Case (MEF90DefMech_unilateralContactTypePrincipalStrains)
                  localAssemblyFunction => MEF90DefMechBilinearFormNull
               Case (MEF90DefMech_unilateralContactTypeMasonry)
                  localAssemblyFunction => MEF90DefMechBilinearFormDamageAT1UnilateralMasonryLoc
               End Select

            Case (MEF90DefMech_damageTypeAT2)
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  QuadratureOrder = 2 * max(elemDisplacementType%order - 1, elemDamageType%order) + 2 * elemDamageType%order
               Else
                  QuadratureOrder = 2 * (elemDisplacementType%order - 1) + 2 * elemDamageType%order
               End If
               Select Case(cellSetOptions%unilateralContactType)
               Case (MEF90DefMech_unilateralContactTypeNone)
                  localAssemblyFunction => MEF90DefMechBilinearFormDamageAT2Loc
               Case (MEF90DefMech_unilateralContactTypeHydrostaticDeviatoric,MEF90DefMech_unilateralContactTypeHybridHydrostaticDeviatoric)
                  localAssemblyFunction => MEF90DefMechBilinearFormDamageAT2UnilateralHDLoc
               Case (MEF90DefMech_unilateralContactTypePositiveHydrostatic)
                  localAssemblyFunction => MEF90DefMechBilinearFormDamageAT2UnilateralPHDLoc
               Case (MEF90DefMech_unilateralContactTypeDeviatoric)
                  localAssemblyFunction => MEF90DefMechBilinearFormDamageAT2UnilateralDeviatoricLoc
               Case (MEF90DefMech_unilateralContactTypePrincipalStrains,MEF90DefMech_unilateralContactTypeMasonry)
                  localAssemblyFunction => MEF90DefMechBilinearFormNull
               End Select

            Case (MEF90DefMech_damageTypeKKL)
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  QuadratureOrder = 2 * max(elemDisplacementType%order - 1, elemDamageType%order) + 2 
               Else
                  QuadratureOrder = 2 * (elemDisplacementType%order - 1) + 2 
               End If
               Select Case(cellSetOptions%unilateralContactType)
               Case (MEF90DefMech_unilateralContactTypeNone)
                  localAssemblyFunction => MEF90DefMechBilinearFormDamageKKLLoc
               End Select
            End Select

            !!! Allocate storage for fields at dof and Gauss points
            Allocate(Aloc(ElemDamageType%numDof,ElemDamageType%numDof))
            Allocate(displacementDof(ElemDisplacementType%numDof))
            Allocate(damageDof(ElemDamageType%numDof))
            If (Associated(MEF90DefMechCtx%temperature)) Then
               Allocate(temperatureDof(ElemDamageType%numDof))
            End If
            
            !!! Allocate elements 
            Call MEF90Element_Create(mesh,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
            Call MEF90Element_Create(mesh,setIS,elemDamage,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
            Do cell = 1,size(cellID)
               !!! Get value of each field at each Dof of the local element
               Call SectionRealRestrictClosure(displacementSec,MEF90DefMechCtx%DMVect,cellID(cell),elemDisplacementType%numDof,displacementDof,ierr);CHKERRQ(ierr)
               Call SectionRealRestrictClosure(damageSec,MEF90DefMechCtx%DMScal,cellID(cell),elemDamageType%numDof,damageDof,ierr);CHKERRQ(ierr)
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  Call SectionRealRestrictClosure(temperatureSec,MEF90DefMechCtx%DMScal,cellID(cell),elemDamageType%numDof,temperatureDof,ierr);CHKERRQ(ierr)
               End If
               If (Associated(MEF90DefMechCtx%plasticStrain)) Then
                  Call SectionRealRestrict(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
                  Call SectionRealRestrict(cumulatedDissipatedPlasticEnergySec,cellID(cell),cumulatedDissipatedPlasticEnergyLoc,ierr);CHKERRQ(ierr)
                  plasticStrainCell = plasticStrainLoc
                  cumulatedDissipatedPlasticEnergyCell = cumulatedDissipatedPlasticEnergyLoc(1)
               Else
                 plasticStrainCell = 0.0_Kr
                 cumulatedDissipatedPlasticEnergyCell = 0.0_Kr
               End If
                  
               Aloc = 0.0_Kr
               Call localAssemblyFunction(Aloc,damageDof,displacementDof,nullPtr,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matpropSet,elemDisplacement(cell),elemDamage(cell))
               If (globalOptions%dampingCoefficientDamage * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
                  Call MEF90_MassMatrixAssembleLoc(Aloc,globalOptions%dampingCoefficientDamage / MEF90DefMechCtx%timeStep,elemDamage(cell))
               End If
               Call DMmeshAssembleMatrix(A,mesh,damageSec,cellID(cell),ALoc,ADD_VALUES,ierr);CHKERRQ(ierr)

               If (Associated(MEF90DefMechCtx%plasticStrain)) Then
                  Call SectionRealRestore(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
                  Call SectionRealRestore(cumulatedDissipatedPlasticEnergySec,cellID(cell),cumulatedDissipatedPlasticEnergyLoc,ierr);CHKERRQ(ierr)
               End If
            End Do
            DeAllocate(displacementDof)
            DeAllocate(damageDof)
            If (Associated(MEF90DefMechCtx%temperature)) Then
               DeAllocate(temperatureDof)
            End If
            Call MEF90Element_Destroy(elemDisplacement,ierr)
            Call MEF90Element_Destroy(elemDamage,ierr)
            DeAllocate(Aloc)    
         End If 
         Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do ! set
      Call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
      Call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
      !!!
      !!! Boundary conditions at cell sets
      !!!
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         If (cellSetOptions%Has_damageBC) Then
            Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
            Call MEF90ISCreateCelltoVertex(mesh,PETSC_COMM_WORLD,setIS,bcIS,ierr)
            Call DMMeshISCreateISglobaldof(mesh,bcIS,0,setISdof,ierr);CHKERRQ(ierr)
            Call MatZeroRowsColumnsIS(A,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
            Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
            Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If ! cellSetOptions%Has_damageBC
      End Do     
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
   
      !!!
      !!! Boundary conditions at vertex sets
      !!!
      Call DMmeshGetLabelIdIS(mesh,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxVertexSetOptions(MEF90DefMechCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
         If (vertexSetOptions%Has_damageBC) Then
            Call DMMeshGetStratumIS(mesh,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
            Call DMMeshISCreateISglobaldof(mesh,setIS,0,setISdof,ierr);CHKERRQ(ierr)
            Call MatZeroRowsColumnsIS(A,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
            Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If
      End Do
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      
      Call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
      Call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)


      Call SectionRealDestroy(displacementSec,ierr);CHKERRQ(ierr)   
      Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)

      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
      End If
      
      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
         Call SectionRealDestroy(cumulatedDissipatedPlasticEnergySec,ierr);CHKERRQ(ierr)
      End If
      flg = SAME_NONZERO_PATTERN
   End Subroutine MEF90DefMechBilinearFormDamage

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSurfaceEnergy"
!!!
!!!  
!!!  MEF90DefMechSurfaceEnergy:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechSurfaceEnergy(alpha,MEF90DefMechCtx,energy,ierr)
      Type(Vec),Intent(IN)                               :: alpha
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,dimension(:),Pointer                     :: energy
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(IS)                                           :: CellSetGlobalIS,setIS,setISdof,bcIS
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(SectionReal)                                  :: alphaSec
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemScal
      Type(MEF90Element_Type)                            :: elemScalType
      PetscReal                                          :: myenergy,myenergySet

      !!! Scatter data from Vec to Sec, or initialize
      Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,alphaSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(alphaSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,alpha,ierr);CHKERRQ(ierr) 

      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         myenergy = 0.0_Kr
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         elemScalType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)

         !!! Call proper local assembly depending on the type of damage law
         Select Case (cellSetOptions%damageType)
         Case (MEF90DefMech_damageTypeAT1,MEF90DefMech_damageTypeAT1Elastic)
            If (elemScalType%coDim == 0) Then
               QuadratureOrder = 2 * (elemScalType%order - 1)
               Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
               Call MEF90GradDamageSurfaceEnergySetAT1(myenergy,alphaSec,MEF90DefMechCtx%DMScal,setIS,matpropSet%internalLength,matpropSet%fractureToughness,elemScal,elemScalType,ierr)
               Call MEF90Element_Destroy(elemScal,ierr)
            End If

         Case (MEF90DefMech_damageTypeAT2,MEF90DefMech_damageTypeAT2Elastic)
            If (elemScalType%coDim == 0) Then
               QuadratureOrder = 2*elemScalType%order
               Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
               Call MEF90GradDamageSurfaceEnergySetAT2(myenergy,alphaSec,MEF90DefMechCtx%DMScal,setIS,matpropSet%internalLength,matpropSet%fractureToughness,elemScal,elemScalType,ierr)
               Call MEF90Element_Destroy(elemScal,ierr)
            End If

         Case (MEF90DefMech_damageTypeKKL,MEF90DefMech_damageTypeKKLElastic)
            If (elemScalType%coDim == 0) Then
               QuadratureOrder = 2*elemScalType%order
               Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
               Call MEF90DefMechSurfaceEnergySetKKL(myenergy,alphaSec,MEF90DefMechCtx%DMScal,setIS,matpropSet%internalLength,matpropSet%fractureToughness,elemScal,elemScalType,ierr)
               Call MEF90Element_Destroy(elemScal,ierr)
            End If

         Case (MEF90DefMech_damageTypeLinSoft)
            If (elemScalType%coDim == 0) Then
               ! Which QuadratureOdre chose for LinSoft
               QuadratureOrder = 5 !* (elemScalType%order - 1)
               Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
               Call MEF90GradDamageSurfaceEnergySetLinSoft(myenergy,alphaSec,MEF90DefMechCtx%DMScal,setIS,matpropSet%internalLength,matpropSet%fractureToughness,elemScal,elemScalType,ierr)
               Call MEF90Element_Destroy(elemScal,ierr)
            End If

         End Select
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         Call MPI_AllReduce(myenergy,myenergySet,1,MPIU_SCALAR,MPI_SUM,MEF90DefMechCtx%MEF90Ctx%comm,ierr);CHKERRQ(ierr)
         energy(set) = energy(set) + myenergySet
      End Do ! set
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr) 

      Call SectionRealDestroy(alphaSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechSurfaceEnergy

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSurfaceEnergySetKKL"
!!!
!!!  
!!!  MEF90DefMechSurfaceEnergySetKKL:
!!!  
!!!  (c) 2014-2019 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechSurfaceEnergySetKKL(energy,alpha,meshScal,cellIS,internalLength,fractureToughness,elemScal,elemScalType,ierr)
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

      C1 = fractureToughness / internalLength / 4.0_Kr / cwKKL
      C2 = fractureToughness * internalLength / 4.0_Kr / cwKKL
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
               energy = energy + elemScal(cell)%Gauss_C(iGauss) * (MEF90gKKL(alphaElem) * C1 + (gradAlphaElem .dotP. gradAlphaElem) * C2)
            End Do ! Gauss
         End Do ! cell
         flops = (2 * elemScalType%numDof + 5 )* size(elemScal(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(alphaloc)
      End If   
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechSurfaceEnergySetKKL

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCrackVolume"
!!!
!!!  
!!!  MEF90DefMechCrackVolume:
!!!  
!!!  (c) 2016 Erwan erwan.tanne@gmail.com  
!!!

   Subroutine MEF90DefMechCrackVolume(xVec,MEF90DefMechCtx,CrackVolume,ierr)
      Type(Vec),Intent(IN)                               :: xVec
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,Dimension(:),Pointer                     :: CrackVolume
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(SectionReal)                                  :: xSec,alphaSec,CrackPressureSec
      Type(IS)                                           :: CellSetGlobalIS,setIS
      PetscInt,Dimension(:),Pointer                      :: setID
      PetscInt                                           :: set,QuadratureOrder,numCell
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elemScal
      Type(MEF90Element_Type)                            :: elemDisplacementType,elemScalType
      Type(MEF90DefMechGlobalOptions_Type),Pointer       :: globalOptions
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      PetscReal                                          :: myCrackVolume,myCrackVolumeSet

      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscReal,Dimension(:),Pointer                     :: xloc,alphaloc,CrackPressureLoc
      PetscInt                                           :: cell
      PetscReal                                          :: CrackPressureElem,alphaElem
      PetscInt                                           :: iDoF1,iGauss
      Type(MEF90_VECT)                                   :: gradAlphaElem,DisplacementElem
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elemDisp

      
      
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%globalOptionsBag,globalOptions,ierr);CHKERRQ(ierr)

      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,xSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(xSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,xVec,ierr);CHKERRQ(ierr)
      
      Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,alphaSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(alphaSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%Damage,ierr);CHKERRQ(ierr) 
      
      CrackVolume   = 0.0_Kr
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         myCrackVolume = 0.0_Kr
         
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
         elemScalType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)

         QuadratureOrder = elemDisplacementType%order
         Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
         QuadratureOrder = 2 * (elemScalType%order - 1)
         Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)

         If (elemScalType%coDim == 0) Then
            
            Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
            If (Size(cellID) > 0) Then
               Allocate(xloc(elemDisplacementType%numDof))
               Allocate(alphaloc(elemScalType%numDof))
               Do cell = 1,size(cellID)

                  Call SectionRealRestrictClosure(xSec,MEF90DefMechCtx%DMVect,cellID(cell),elemDisplacementType%numDof,xloc,ierr);CHKERRQ(ierr)
                  Call SectionRealRestrictClosure(alphaSec,MEF90DefMechCtx%DMScal,cellID(cell),elemScalType%numDof,alphaLoc,ierr);CHKERRQ(ierr)

                  Do iGauss = 1,size(elemScal(cell)%Gauss_C)

                     DisplacementElem = 0.0_Kr
                     Do iDoF1 = 1,elemDisplacementType%numDof
                        DisplacementElem = DisplacementElem + xloc(iDof1) * elemDisplacement(cell)%BF(iDof1,iGauss)
                     End Do
                     gradAlphaElem = 0.0_Kr
                     alphaElem = 0.0_Kr
                     Do iDoF1 = 1,elemScalType%numDof
                        alphaElem     = alphaElem     + alphaLoc(iDof1) * elemScal(cell)%BF(iDoF1,iGauss)
                        gradAlphaElem = gradAlphaElem + alphaLoc(iDof1) * elemScal(cell)%Grad_BF(iDoF1,iGauss)
                     End Do
                     myCrackVolume = myCrackVolume - (elemDisplacement(cell)%Gauss_C(iGauss) *( ( gradAlphaElem) .dotP. DisplacementElem  ) )
                  End Do ! Gauss
               End Do ! cell

               DeAllocate(xloc)
               DeAllocate(alphaloc)
               End If
            Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
         End If

         Call MEF90Element_Destroy(elemScal,ierr)
         Call MEF90Element_Destroy(elemDisplacement,ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         Call MPI_AllReduce(myCrackVolume,myCrackVolumeSet,1,MPIU_SCALAR,MPI_SUM,MEF90DefMechCtx%MEF90Ctx%comm,ierr);CHKERRQ(ierr)
         CrackVolume(set) = CrackVolume(set) + myCrackVolumeSet
      End Do ! set
      Call ISrestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr) 

      Call SectionRealDestroy(alphaSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechCrackVolume



#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCrackPressureRescaling"
!!!
!!!  
!!!  MEF90DefMechCrackPressureRescaling:
!!!  
!!!  (c) 2016 Erwan erwan.tanne@gmail.com  
!!!

   Subroutine MEF90DefMechCrackPressureRescaling(xVec,CrackPressureVec,MEF90DefMechCtx,CrackVolumeSet,ActivatedCrackPressureBlocksList,timeStep,ierr)
      Type(Vec),Intent(IN)                               :: xVec
      Type(Vec),Intent(INOUT)                            :: CrackPressureVec
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      PetscReal,Dimension(:),Pointer,Intent(IN)          :: CrackVolumeSet
      PetscReal,Intent(IN)                               :: timeStep
      PetscBool,Dimension(:),Pointer,Intent(IN)          :: ActivatedCrackPressureBlocksList


      Type(MEF90DefMechGlobalOptions_Type),Pointer       :: globalOptions
      Type(SectionReal)                                  :: xSec,CrackPressureSec,alphaSec
      Type(IS)                                           :: CellSetGlobalIS,setIS
      PetscInt,Dimension(:),Pointer                      :: setID,cellID
      PetscInt                                           :: set,QuadratureOrder,numCell
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elemScal
      Type(MEF90Element_Type)                            :: elemDisplacementType,elemScalType
      PetscReal                                          :: cellSize,CrackPressureElem,alphaElem,PressureVolumeElem,CrackVolumeControlled
      Type(MEF90_MATS)                                   :: StrainElem
      Type(MEF90_VECT)                                   :: gradAlphaElem,DisplacementElem
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      PetscInt                                           :: iDoF1,iGauss,cell
      PetscReal,Dimension(:),Pointer                     :: xloc,alphaloc,CrackPressureLoc
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elemDisp
      Type(MEF90_MATPROP),pointer                        :: matpropSet


      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%globalOptionsBag,globalOptions,ierr);CHKERRQ(ierr)

      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,xSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(xSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,xVec,ierr);CHKERRQ(ierr)

      Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,alphaSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(alphaSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%Damage,ierr);CHKERRQ(ierr) 

      Call SectionRealDuplicate(MEF90DefMechCtx%cellDMScalSec,CrackPressureSec,ierr)
      Call SectionRealToVec(CrackPressureSec,MEF90DefMechCtx%CellDMScalScatter,SCATTER_REVERSE,CrackPressureVec,ierr);CHKERRQ(ierr) 

      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)

         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
         elemScalType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)
         !Allocate(ActivatedCrackPressureBlocksList(size(MEF90DefMechCtx%CellSetOptionsBag)))

         QuadratureOrder = elemDisplacementType%order
         Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
         QuadratureOrder = 2 * (elemScalType%order - 1)
         Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)

         If (elemScalType%coDim == 0) Then

            Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
            If (Size(cellID) > 0) Then
               Allocate(xloc(elemDisplacementType%numDof))
               Allocate(alphaloc(elemScalType%numDof))

               PressureVolumeElem = 0.0_Kr
               alphaElem = 0.0_Kr
               Do cell = 1,size(cellID)
                  cellSize = 0.0_Kr
                  PressureVolumeElem = 0.0_Kr

                  Call SectionRealRestrictClosure(xSec,MEF90DefMechCtx%DMVect,cellID(cell),elemDisplacementType%numDof,xloc,ierr);CHKERRQ(ierr)
                  Call SectionRealRestrictClosure(alphaSec,MEF90DefMechCtx%DMScal,cellID(cell),elemScalType%numDof,alphaLoc,ierr);CHKERRQ(ierr)
                  Call SectionRealRestrict(CrackPressureSec,cellID(cell),CrackPressureLoc,ierr);CHKERRQ(ierr)
                  Do iGauss = 1,size(elemScal(cell)%Gauss_C)
                     StrainElem = 0.0_Kr
                     Do iDoF1 = 1,elemDisplacementType%numDof
                        StrainElem = StrainElem + xloc(iDof1) *  elemDisplacement(cell)%GradS_BF(iDof1,iGauss)
                     End Do

                     DisplacementElem = 0.0_Kr
                     Do iDoF1 = 1,elemDisplacementType%numDof
                        DisplacementElem = DisplacementElem + xloc(iDof1) * elemDisplacement(cell)%BF(iDof1,iGauss)
                     End Do

                     gradAlphaElem = 0.0_Kr
                     alphaElem = 0.0_Kr
                     Do iDoF1 = 1,elemScalType%numDof
                        alphaElem     = alphaElem     + alphaLoc(iDof1) * elemScal(cell)%BF(iDoF1,iGauss)
                        gradAlphaElem = gradAlphaElem + alphaLoc(iDof1) * elemScal(cell)%Grad_BF(iDoF1,iGauss)
                     End Do
                     cellSize = cellSize + elemDisplacement(cell)%Gauss_C(iGauss)
                     PressureVolumeElem = PressureVolumeElem - (elemDisplacement(cell)%Gauss_C(iGauss) *( gradAlphaElem .dotP. DisplacementElem  ) )
                  End Do ! Gauss

                  !! Case where the total volume is controlled
                  If (ActivatedCrackPressureBlocksList(set)) then
                     !! rescal the total 
                     CrackVolumeControlled = sum(CrackVolumeSet,MASK=ActivatedCrackPressureBlocksList)
                     CrackPressureLoc = ( timeStep/CrackVolumeControlled )*CrackPressureLoc
                  Else 
                     If (trace(StrainElem)<0 .and. (alphaElem > 0.0_Kr)) Then
                        CrackPressureLoc = (1.0_Kr-alphaElem)**2*((matpropSet%HookesLaw*StrainElem)*gradAlphaElem)*gradAlphaElem / (gradAlphaElem .dotP. gradAlphaElem)
                     Else
                        CrackPressureLoc = 0.0_Kr
                     EndIf
                  EndIf
                  Call SectionRealRestore(CrackPressureSec,cellID(cell),CrackPressureLoc,ierr);CHKERRQ(ierr)
               End Do ! cell

               DeAllocate(xloc)
               DeAllocate(alphaloc)
               End If
            Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
         End If

         Call MEF90Element_Destroy(elemScal,ierr)
         Call MEF90Element_Destroy(elemDisplacement,ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do ! set
      Call ISrestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call SectionRealToVec(CrackPressureSec,MEF90DefMechCtx%CellDMScalScatter,SCATTER_FORWARD,CrackPressureVec,ierr);CHKERRQ(ierr) 
      
      Call SectionRealDestroy(alphaSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(CrackpressureSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechCrackPressureRescaling

End Module MEF90_APPEND(m_MEF90_DefMechAssembly,MEF90_DIM)D
