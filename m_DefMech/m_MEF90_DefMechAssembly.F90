#include "../MEF90/mef90.inc"
Module MEF90_APPEND(m_MEF90_DefMechAssembly,MEF90_DIM)D
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
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
          MEF90DefMechStress

   Abstract Interface
      Subroutine MEF90DefMechBilinearFormLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
         Use m_MEF90_DefMechCtx
         PetscReal,Dimension(:,:),Pointer                   :: Aloc
         PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
         Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
         Type(MEF90_MATPROP),Intent(IN)                     :: matprop
         Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
         Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      End Subroutine MEF90DefMechBilinearFormLoc
   End Interface

   Abstract Interface
      Subroutine MEF90DefMechOperatorLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
         Use m_MEF90_DefMechCtx
         PetscReal,Dimension(:),Pointer                     :: residualLoc
         PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
         Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
         Type(MEF90_MATPROP),Intent(IN)                     :: matprop
         Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
         Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      End Subroutine MEF90DefMechOperatorLoc
   End Interface

   Abstract Interface
      Subroutine MEF90DefMechRHSLoc(residualLoc,xDof,forceCell,pressureCell,matprop,elemDisplacement)
         Use m_MEF90_DefMechCtx
         PetscReal,Dimension(:),Pointer                     :: residualLoc
         PetscReal,Dimension(:),Pointer                     :: xDof
         Type(MEF90_VECT),Intent(IN)                        :: forceCell
         PetscReal,Intent(IN)                               :: pressureCell
         Type(MEF90_MATPROP),Intent(IN)                     :: matprop
         Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      End Subroutine MEF90DefMechRHSLoc
   End Interface
   
   Abstract Interface
      Subroutine MEF90DefMechEnergyLoc(energyLoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
         Use m_MEF90_DefMechCtx
         PetscReal,Dimension(:),Pointer                     :: residualLoc
         PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
         Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
         Type(MEF90_MATPROP),Intent(IN)                     :: matprop
         Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
         Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
      End Subroutine MEF90DefMechEnergyLoc
   End Interface

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorNull"
   Subroutine MEF90DefMechOperatorNull(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      Use m_MEF90_DefMechCtx
      PetscReal,Dimension(:),Pointer                        :: residualLoc
      PetscReal,Dimension(:),Pointer                        :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
         Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
         Type(MEF90_MATPROP),Intent(IN)                     :: matprop
         Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
         Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      Print*,__FUNCT__,': Unimplemented Operator local assembly function'
      STOP
   End Subroutine MEF90DefMechOperatorNull

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormNull"
   Subroutine MEF90DefMechBilinearFormNull(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
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

   Subroutine MEF90DefMechBilinearFormDisplacementElasticLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numGauss
      Type(MEF90_MATS)                                   :: sigma
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numGauss = size(elemDisplacement%BF,2)
      Aloc = 0.0_Kr
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
!!!  MEF90DefMechBilinearFormDisplacementElasticLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDisplacementATLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: ALoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
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
      ALoc = 0.0_Kr
      Do iGauss = 1,numGauss
         stiffness = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            stiffness = stiffness + damageDof(iDoF1) * elemDamage%BF(iDoF1,iGauss)
         End Do
         stiffness = (1.0_Kr - stiffness)**2 + matProp%residualStiffness
         Do iDoF1 = 1,numDofDisplacement
            sigma = stiffness * (matProp%HookesLaw * elemDisplacement%GradS_BF(iDoF1,iGauss))
            Do iDoF2 = 1,numDofDisplacement
               ALoc(iDoF2,iDoF1) = ALoc(iDoF2,iDoF1) + elemDisplacement%Gauss_C(iGauss) * ( (sigma .DotP. elemDisplacement%GradS_BF(iDoF2,iGauss)) + matprop%cohesiveStiffness * elemDisplacement%BF(iDoF2,iGauss) * elemDisplacement%BF(iDoF1,iGauss))
               If (matprop%cohesiveStiffness /= 0.0_Kr) Then
                  ALoc(iDoF2,iDoF1) = ALoc(iDoF2,iDoF1) + elemDisplacement%Gauss_C(iGauss) * matprop%cohesivestiffness * (elemDisplacement%BF(iDoF1,iGauss) .DotP. elemDisplacement%BF(iDoF2,iGauss))
               End If            
            End Do
         End Do
      End Do
      flops = numGauss * ( 4. * numDofDisplacement**2 + 2. * numDofDamage + 3.)
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDisplacementATLoc
   
!!erwan-->!!
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDisplacementATkLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDisplacementElasticLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDisplacementATkLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: ALoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
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
      ALoc = 0.0_Kr
      k = matprop%k_for_ATk
      Do iGauss = 1,numGauss
         stiffness = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            stiffness = stiffness + damageDof(iDoF1) * elemDamage%BF(iDoF1,iGauss)
         End Do
         stiffness = (1.0_Kr - stiffness)**2 /(1.0_Kr + ( k - 1.0_Kr ) * ( 1.0_Kr - (1.0_Kr - stiffness)**2 ) ) + matProp%residualStiffness
         Do iDoF1 = 1,numDofDisplacement
            sigma = stiffness * (matProp%HookesLaw * elemDisplacement%GradS_BF(iDoF1,iGauss))
            Do iDoF2 = 1,numDofDisplacement
               ALoc(iDoF2,iDoF1) = ALoc(iDoF2,iDoF1) + elemDisplacement%Gauss_C(iGauss) * ( (sigma .DotP. elemDisplacement%GradS_BF(iDoF2,iGauss)) + matprop%cohesiveStiffness * elemDisplacement%BF(iDoF2,iGauss) * elemDisplacement%BF(iDoF1,iGauss))
               If (matprop%cohesiveStiffness /= 0.0_Kr) Then
                  ALoc(iDoF2,iDoF1) = ALoc(iDoF2,iDoF1) + elemDisplacement%Gauss_C(iGauss) * matprop%cohesivestiffness * (elemDisplacement%BF(iDoF1,iGauss) .DotP. elemDisplacement%BF(iDoF2,iGauss))
               End If            
            End Do
         End Do
      End Do
      !flops = numGauss * ( 2. * numDofDisplacement**2 + 2. * numDofDamage + 3.)
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDisplacementATkLoc
   
!!<--erwan!!




#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDisplacementATUnilateralHDLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDisplacementATUnilateralHDLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDisplacementATUnilateralHDLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: ALoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: stiffnessH,stiffnessD
      PetscReal                                          :: inelasticStrainTrace
      Type(MEF90_MATS)                                   :: plasticStrain,sigma
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDisplacement%BF,2)

      ALoc = 0.0_Kr
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

         If (inelasticStrainTrace < 0.0_Kr) Then
            stiffnessH = 1.0_Kr
         Else
            stiffnessH = stiffnessD
         End If
         
         Do iDoF1 = 1,numDofDisplacement
            sigma = stiffnessH * hydrostaticPart(elemDisplacement%GradS_BF(iDoF1,iGauss)) + stiffnessD * deviatoricPart(elemDisplacement%GradS_BF(iDoF1,iGauss))
            sigma = matProp%HookesLaw * sigma
            Do iDoF2 = 1,numDofDisplacement
               ALoc(iDoF2,iDoF1) = ALoc(iDoF2,iDoF1) + elemDisplacement%Gauss_C(iGauss) * (sigma .DotP. elemDisplacement%GradS_BF(iDoF2,iGauss))
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

   Subroutine MEF90DefMechBilinearFormDisplacementATUnilateralDeviatoricLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: ALoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
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

      ALoc = 0.0_Kr
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

   Subroutine MEF90DefMechBilinearFormDisplacementATUnilateralPHDLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: ALoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
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

      ALoc = 0.0_Kr
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

   Subroutine MEF90DefMechBilinearFormDisplacementATUnilateralPSLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: ALoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
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

      ALoc = 0.0_Kr
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

   Subroutine MEF90DefMechOperatorDisplacementElasticLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

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
   
         If (Associated(boundaryDisplacementDof)) Then
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

   Subroutine MEF90DefMechOperatorDisplacementATLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      Type(MEF90_MATS)                                   :: sigma
      Type(MEF90_VECT)                                   :: UU0
      PetscReal                                          :: stiffness,temperature
      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDisplacement%BF,2)

      residualLoc = 0.0_Kr
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
         sigma = stiffness * (matProp%HookesLaw * (sigma - temperature * matProp%LinearThermalExpansion - plasticStrainCell))

         Do iDoF2 = 1,numDofDisplacement
            residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement%Gauss_C(iGauss) * (sigma .DotP. elemDisplacement%GradS_BF(iDoF2,iGauss))
         End Do
         If (Associated(boundaryDisplacementDof)) Then
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

!!erwan-->!!
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDisplacementATkLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDisplacementATkLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDisplacementATkLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

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
      k = matprop%k_for_ATk
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
         If (Associated(boundaryDisplacementDof)) Then
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
   End Subroutine MEF90DefMechOperatorDisplacementATkLoc
!!<--erwan!!



#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDisplacementATUnilateralHDLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDisplacementATUnilateralHDLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDisplacementATUnilateralHDLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: temperature,stiffnessH,stiffnessD
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

         stiffnessD = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            stiffnessD = stiffnessD + damageDof(iDoF1) * elemDamage%BF(iDoF1,iGauss)
         End Do
         stiffnessD = (1.0_Kr - stiffnessD)**2 + matProp%residualStiffness

         If (Trace(inelasticStrain) < 0.0_Kr) Then
            stiffnessH = 1.0_Kr
         Else
            stiffnessH = stiffnessD
         End If
         
         sigma = stiffnessH * hydrostaticPart(inelasticStrain) + stiffnessD * deviatoricPart(inelasticStrain)
         sigma = matProp%HookesLaw * sigma
         Do iDoF2 = 1,numDofDisplacement
            residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement%Gauss_C(iGauss) * (sigma .DotP. elemDisplacement%GradS_BF(iDoF2,iGauss))
         End Do
      End Do
      !flops = numGauss * ( 2. * numDofDisplacement**2 + 3. * numDofDamage + 2.)
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

   Subroutine MEF90DefMechOperatorDisplacementATUnilateralPHDLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

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

   Subroutine MEF90DefMechOperatorDisplacementATUnilateralDeviatoricLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

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

   Subroutine MEF90DefMechOperatorDisplacementATUnilateralPSLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

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
#define __FUNCT__ "MEF90DefMechRHSDisplacementLoc"
!!!
!!!  
!!!  MEF90DefMechRHSDisplacementLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechRHSDisplacementLoc(residualLoc,xDof,forceCell,pressureForceCell,matprop,elemDisplacement)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof
      Type(MEF90_VECT),Intent(IN)                        :: forceCell
      PetscReal,Intent(IN)                               :: pressureForceCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      
      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numGauss
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr
      Type(MEF90_VECT)                                   :: totalForce

      totalForce = forceCell + pressureForceCell * elemDisplacement%outerNormal
      If (Norm(totalForce) > 0.0_Kr) Then
         numDofDisplacement = size(elemDisplacement%BF,1)
         numGauss = size(elemDisplacement%BF,2)
         residualLoc = 0.0_Kr
         Do iGauss = 1,numGauss
            Do iDoF2 = 1,numDofDisplacement
               residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement%Gauss_C(iGauss) * &
                                  (totalForce .DotP. elemDisplacement%BF(iDoF2,iGauss))
            End Do
         End Do
         !flops = 2 * numGauss * numDofDisplacement**2
         !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      End If
   End Subroutine MEF90DefMechRHSDisplacementLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDisplacement"
!!!
!!!  
!!!  MEF90DefMechOperatorDisplacement: Build the operator. When called in SNES, the solution time should always match the target time, 
!!!                                    so there is no need for interpolation of the forcees, external, and boundary values
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDisplacement(snesDisplacement,displacement,residual,MEF90DefMechCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesDisplacement
      Type(Vec),Intent(IN)                               :: displacement
      Type(Vec),Intent(INOUT)                            :: residual
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(DM)                                           :: mesh
      Type(SectionReal)                                  :: displacementSec,residualSec
      Type(SectionReal)                                  :: plasticStrainSec,temperatureSec,damageSec
      Type(SectionReal)                                  :: boundaryDisplacementSec,forceSec,pressureForceSec
      PetscReal,Dimension(:),Pointer                     :: displacementPtr,residualPtr
      PetscReal,Dimension(:),Pointer                     :: boundaryDisplacementPtr
      PetscReal,Dimension(:),Pointer                     :: boundaryDisplacementDof,displacementDof,damageDof,temperatureDof,residualLoc
      PetscReal,Dimension(:),Pointer                     :: pressureForceLoc
      PetscReal                                          :: pressureForceCell
      PetscReal,Dimension(:),Pointer                     :: forceLoc
      Type(MEF90_VECT)                                   :: forceCell
      PetscReal,Dimension(:),Pointer                     :: plasticStrainLoc
      Type(MEF90_MATS)                                   :: plasticStrainCell
      Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof,bcIS
      PetscInt,dimension(:),Pointer                      :: setID,cellID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder,cell
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90DefMechVertexSetOptions_Type),Pointer    :: vertexSetOptions
      Type(VecScatter)                                   :: ScatterSecToVec,ScatterSecToVecScal
      Type(VecScatter)                                   :: ScatterSecToVecCell,ScatterSecToVecCellScal,ScatterSecToVecCellMatS
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemDamage
      Type(MEF90Element_Type)                            :: elemDisplacementType,elemDamageType
      PetscInt                                           :: p,dof,c,dim,numDof
      PetscInt                                           :: nVal
      PetscReal,Dimension(:),Pointer                     :: nullPtr
      
      Procedure(MEF90DefMechOperatorLoc),pointer         :: localOperatorFunction
      Procedure(MEF90DefMechRHSLoc),pointer              :: localRHSFunction
      
      localOperatorFunction =>MEF90DefMechOperatorNull      
      Call SNESGetDM(snesDisplacement,mesh,ierr);CHKERRQ(ierr)

      !!! Get Section and scatter associated with each field
      Call DMMeshGetSectionReal(MEF90DefMechCtx%DMVect,'default',displacementSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMVect,displacementSec,ScatterSecToVec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(displacementSec,residualSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(displacementSec,boundaryDisplacementSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(displacementSec,ScatterSecToVec,SCATTER_REVERSE,displacement,ierr);CHKERRQ(ierr) 
      Call SectionRealToVec(boundaryDisplacementSec,ScatterSecToVec,SCATTER_REVERSE,MEF90DefMechCtx%boundaryDisplacement,ierr);CHKERRQ(ierr)


      Call SectionRealSet(residualSec,0.0_Kr,ierr);CHKERRQ(ierr)
      Call VecSet(residual,0.0_Kr,ierr);CHKERRQ(ierr)

      Call DMMeshGetSectionReal(MEF90DefMechCtx%DMScal,'default',damageSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMScal,damageSec,ScatterSecToVecScal,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(damageSec,ScatterSecToVecScal,SCATTER_REVERSE,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)   

      Nullify(TemperatureDof)
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDuplicate(damageSec,temperatureSec,ierr);CHKERRQ(ierr)
         Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMScal,temperatureSec,ScatterSecToVecScal,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(temperatureSec,ScatterSecToVecScal,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr)   
      Else
         temperatureSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%force)) Then
         Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMVect,'default',forceSec,ierr);CHKERRQ(ierr)
         Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMVect,forceSec,ScatterSecToVecCell,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(forceSec,ScatterSecToVecCell,SCATTER_REVERSE,MEF90DefMechCtx%force,ierr);CHKERRQ(ierr)
      Else
         ForceSec%v = 0
      End If
      
      If (Associated(MEF90DefMechCtx%pressureForce)) Then
         Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMScal,'default',pressureForceSec,ierr);CHKERRQ(ierr)
         Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMScal,pressureForceSec,ScatterSecToVecCellScal,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(pressureForceSec,ScatterSecToVecCellScal,SCATTER_REVERSE,MEF90DefMechCtx%pressureForce,ierr);CHKERRQ(ierr)
      Else
         pressureForceSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMMatS,'default',plasticStrainSec,ierr);CHKERRQ(ierr)
         Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMMatS,plasticStrainSec,ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticStrainSec,ScatterSecToVecCellMatS,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)          
      Else
         PlasticStrainSec%v = 0
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
            Case (MEF90DefMech_damageTypeAT1Elastic,MEF90DefMech_damageTypeAT2Elastic)
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  QuadratureOrder = elemDisplacementType%order - 1 + max(elemDisplacementType%order - 1, elemDamageType%order)
               Else
                  QuadratureOrder = 2 * (elemDisplacementType%order - 1)
               End If
               localOperatorFunction => MEF90DefMechOperatorDisplacementElasticLoc
               localRHSFunction => MEF90DefMechRHSDisplacementLoc
!!erwan-->!!
            Case (MEF90DefMech_damageTypeATk)
               QuadratureOrder = 5 !* (elemDisplacementType%order - 1) + 5 * elemDamageType%order
               localOperatorFunction => MEF90DefMechOperatorDisplacementATkLoc
               localRHSFunction => MEF90DefMechRHSDisplacementLoc
!!<--erwan!!
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
                  localOperatorFunction => MEF90DefMechOperatorNull
               End Select
               localRHSFunction => MEF90DefMechRHSDisplacementLoc
            End Select

            !!! Allocate storage for fields at dof and Gauss points
            Allocate(residualLoc(ElemDisplacementType%numDof))
            Allocate(displacementDof(ElemDisplacementType%numDof))
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
               If (Associated(MEF90DefMechCtx%plasticStrain)) Then
                  Call SectionRealRestrict(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
                  plasticStrainCell = plasticStrainLoc
               Else
                  plasticStrainCell = 0.0_Kr
               End If
                  
               residualLoc = 0.0_Kr
               Call localOperatorFunction   (residualLoc,displacementDof,nullPtr,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,matpropSet,elemDisplacement(cell),elemDamage(cell))
               Call localRHSFunction        (residualLoc,displacementDof,forceCell,pressureForceCell,matpropSet,elemDisplacement(cell))
               Call SectionRealUpdateClosure(residualSec,MEF90DefMechCtx%DMVect,cellID(cell),residualLoc,ADD_VALUES,ierr);CHKERRQ(iErr)

               If (Associated(MEF90DefMechCtx%force)) Then
                  Call SectionRealRestore(forceSec,cellID(cell),forceLoc,ierr);CHKERRQ(ierr)
               End If
               If (Associated(MEF90DefMechCtx%pressureForce)) Then
                  Call SectionRealRestore(pressureForceSec,cellID(cell),pressureForceLoc,ierr);CHKERRQ(ierr)
               End If
               If (Associated(MEF90DefMechCtx%plasticStrain)) Then
                  Call SectionRealRestore(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
               End If
            End Do

            DeAllocate(displacementDof)
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
            Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If 
      End Do ! set

      !!! "Ghost update" for the residual Section
      Call SectionRealComplete(residualSec,ierr);CHKERRQ(ierr)
      !!! Scatter back from SectionReal to Vec
      Call SectionRealToVec(residualSec,ScatterSecToVec,SCATTER_FORWARD,residual,ierr);CHKERRQ(ierr)

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
         Call VecScatterDestroy(ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)      
      End If
      
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
      End If

      If (Associated(MEF90DefMechCtx%pressureForce)) Then
         Call SectionRealDestroy(pressureForceSec,ierr);CHKERRQ(ierr)
         Call VecScatterDestroy(ScatterSecToVecCellScal,ierr);CHKERRQ(ierr)      
      End If
      
      If (Associated(MEF90DefMechCtx%damage)) Then
         Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)
      End If
      
      Call SectionRealDestroy(forceSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(pressureForceSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(boundaryDisplacementSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(residualSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(displacementSec,ierr);CHKERRQ(ierr)

      Call VecScatterDestroy(ScatterSecToVecScal,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVecCell,ierr);CHKERRQ(ierr)      
      Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDisplacement

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDisplacement"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDisplacement:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDisplacement(snesDispl,displacement,A,M,flg,MEF90DefMechCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesDispl
      Type(Vec),Intent(IN)                               :: displacement
      Type(Mat),Intent(INOUT)                            :: A,M
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr  
         
      Type(SectionReal)                                  :: displacementSec,damageSec,temperatureSec,plasticStrainSec
      PetscReal,Dimension(:),Pointer                     :: displacementDOF,damageDof,temperatureDof,nullPtr
      Type(VecScatter)                                   :: ScatterSecToVec,ScatterSecToVecScal,ScatterSecToVecCellMatS
      PetscReal,Dimension(:),Pointer                     :: plasticStrainLoc
      Type(MEF90_MATS)                                   :: plasticStrainCell
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
      
      localAssemblyFunction =>MEF90DefMechBilinearFormNull      
      Call MatZeroEntries(A,ierr);CHKERRQ(ierr)
      Call SNESGetDM(snesDispl,mesh,ierr);CHKERRQ(ierr)
      
      !!! Get Section and scatter associated with each field
      Call DMMeshGetSectionReal(MEF90DefMechCtx%DMVect,'default',displacementSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMVect,displacementSec,ScatterSecToVec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(displacementSec,ScatterSecToVec,SCATTER_REVERSE,displacement,ierr);CHKERRQ(ierr)   

      Call DMMeshGetSectionReal(MEF90DefMechCtx%DMScal,'default',damageSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMScal,damageSec,ScatterSecToVecScal,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(damageSec,ScatterSecToVecScal,SCATTER_REVERSE,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)   

      Nullify(temperatureDof)
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDuplicate(damageSec,temperatureSec,ierr);CHKERRQ(ierr)
         !Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMScal,temperatureSec,ScatterSecToVecScal,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(temperatureSec,ScatterSecToVecScal,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr)   
      Else
         temperatureSec%v = 0
      End If
      
      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMMatS,'default',plasticStrainSec,ierr);CHKERRQ(ierr)
         Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMMatS,plasticStrainSec,ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticStrainSec,ScatterSecToVecCellMatS,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)          
      Else
         PlasticStrainSec%v = 0
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
            Case (MEF90DefMech_damageTypeAT1Elastic,MEF90DefMech_damageTypeAT2Elastic)
               QuadratureOrder = 2 * (elemDisplacementType%order - 1)
               localAssemblyFunction => MEF90DefMechBilinearFormDisplacementElasticLoc
!!erwan-->!!
            Case (MEF90DefMech_damageTypeATk)
               QuadratureOrder = 5 !* (elemDisplacementType%order - 1) + 5 * ElemDamageType%order
               localAssemblyFunction => MEF90DefMechBilinearFormDisplacementATkLoc
!!<--erwan!!
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
               Case (MEF90DefMech_unilateralContactTypePrincipalStrains)
                  localAssemblyFunction => MEF90DefMechBilinearFormNull
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
                  plasticStrainCell = plasticStrainLoc
               Else
                  plasticStrainCell = 0.0_Kr
               End If
                  
               Aloc = 0.0_Kr
               Call localAssemblyFunction(Aloc,displacementDof,nullPtr,damageDof,temperatureDof,plasticStrainCell,matpropSet,elemDisplacement(cell),elemDamage(cell))
               Call DMmeshAssembleMatrix(A,mesh,displacementSec,cellID(cell),ALoc,ADD_VALUES,ierr);CHKERRQ(ierr)

               If (Associated(MEF90DefMechCtx%plasticStrain)) Then
                  Call SectionRealRestore(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
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
            Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If 
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

      Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVecScal,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(displacementSec,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)
      End If      
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
      Type(VecScatter)                                   :: ScatterSecToVec
      Type(VecScatter)                                   :: ScatterSecToVecCell,ScatterSecToVecCellScal
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90Element_Type)                            :: elemDisplacementType
      Type(MEF90DefMechGlobalOptions_Type),Pointer       :: globalOptions
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      PetscReal                                          :: myWork,myWorkSet
      
      
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%globalOptionsBag,globalOptions,ierr);CHKERRQ(ierr)
      !!! Create dof-based sections
      Call DMMeshGetSectionReal(MEF90DefMechCtx%DMVect,'default',xSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMVect,xSec,ScatterSecToVec,ierr);CHKERRQ(ierr)

      !!! Create cell based sections, and allocate required Pointers
      !!! I could check if any of these vectors are null 
      Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMVect,'default',ForceSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMVect,forceSec,ScatterSecToVecCell,ierr);CHKERRQ(ierr)
      Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMScal,'default',pressureForceSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMScal,pressureForceSec,ScatterSecToVecCellScal,ierr);CHKERRQ(ierr)
   
      !!! Scatter data from Vec to Sec, or initialize
      Call SectionRealToVec(xSec,ScatterSecToVec,SCATTER_REVERSE,xVec,ierr);CHKERRQ(ierr)         
      Call SectionRealToVec(forceSec,ScatterSecToVecCell,SCATTER_REVERSE,MEF90DefMechCtx%force,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(pressureForceSec,ScatterSecToVecCellScal,SCATTER_REVERSE,MEF90DefMechCtx%pressureForce,ierr);CHKERRQ(ierr)

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

      Call VecScatterDestroy(ScatterSecToVecCellScal,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(pressureForceSec,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVecCell,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(forceSec,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
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
      Type(VecScatter)                                   :: ScatterSecToVec
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90Element_Type)                            :: elemDisplacementType
      Type(MEF90DefMechGlobalOptions_Type),Pointer       :: globalOptions
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      PetscReal                                          :: mycohesiveEnergy,mycohesiveEnergySet
      
      
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%globalOptionsBag,globalOptions,ierr);CHKERRQ(ierr)
      !!! Create dof-based sections
      Call DMMeshGetSectionReal(MEF90DefMechCtx%DMVect,'default',xSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMVect,xSec,ScatterSecToVec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(xSec,boundaryDisplacementSec,ierr);CHKERRQ(ierr)
      
      !!! Scatter data from Vec to Sec, or initialize
      Call SectionRealToVec(xSec,ScatterSecToVec,SCATTER_REVERSE,xVec,ierr);CHKERRQ(ierr)         
      Call SectionRealToVec(boundaryDisplacementSec,ScatterSecToVec,SCATTER_REVERSE,MEF90DefMechCtx%boundaryDisplacement,ierr);CHKERRQ(ierr)

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

      Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
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
      Type(SectionReal)                                  :: damageSec,plasticStrainSec,plasticStrainOldSec,temperatureSec
      Type(IS)                                           :: CellSetGlobalIS,setIS,setISdof,bcIS
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(VecScatter)                                   :: ScatterSecToVec,ScatterSecToVecScal,ScatterSecToVecCellMatS
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemScal
      Type(MEF90Element_Type)                            :: elemDisplacementType,elemScalType
      PetscReal                                          :: myenergy,myenergySet

      !!! Create dof-based sections
      Call DMMeshGetSectionReal(MEF90DefMechCtx%DMVect,'default',xSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMVect,xSec,ScatterSecToVec,ierr);CHKERRQ(ierr)

      !!! Scatter data from Vec to Sec, or initialize
      Call SectionRealToVec(xSec,ScatterSecToVec,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr) 

      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMMatS,'default',plasticStrainSec,ierr);CHKERRQ(ierr)
         Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMMatS,plasticStrainSec,ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticStrainSec,ScatterSecToVecCellMatS,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)       
         
         Call SectionRealDuplicate(plasticStrainSec,plasticStrainOldSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticStrainOldSec,ScatterSecToVecCellMatS,SCATTER_REVERSE,plasticStrainOld,ierr);CHKERRQ(ierr)          
      Else
         PlasticStrainSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call DMMeshGetSectionReal(MEF90DefMechCtx%DMScal,'default',temperatureSec,ierr);CHKERRQ(ierr)
         Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMScal,temperatureSec,ScatterSecToVecScal,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(temperatureSec,ScatterSecToVecScal,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr)          
      Else
         temperatureSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%damage)) Then
         If (Associated(MEF90DefMechCtx%temperature)) Then
            Call SectionRealDuplicate(temperatureSec,damageSec,ierr);CHKERRQ(ierr)
         Else
            Call DMMeshGetSectionReal(MEF90DefMechCtx%DMScal,'default',damageSec,ierr);CHKERRQ(ierr)
            Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMScal,damageSec,ScatterSecToVecScal,ierr);CHKERRQ(ierr)
         End If
         Call SectionRealToVec(damageSec,ScatterSecToVecScal,SCATTER_REVERSE,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)          
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
         Case (MEF90DefMech_damageTypeAT1Elastic,MEF90DefMech_damageTypeAT2Elastic,MEF90DefMech_damageTypeAT1,MEF90DefMech_damageTypeATk)
            If (elemDisplacementType%coDim == 0) Then
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  QuadratureOrder = 2 * max(elemDisplacementType%order - 1, elemScalType%order)
               Else
                  QuadratureOrder = 2 * (elemDisplacementType%order - 1)
               End If
               Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
               Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
               Call MEF90PlasticityEnergySet(myenergy,xSec,damageSec,plasticStrainSec,plasticStrainOldSec,temperatureSec,MEF90DefMechCtx%DMVect,MEF90DefMechCtx%DMScal,setIS, &
                                             matpropSet%HookesLaw,matpropSet%LinearThermalExpansion,elemDisplacement,elemDisplacementType,elemScal,elemScalType,ierr)
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
      If (Associated(MEF90DefMechCtx%damage)) Then
         Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)
      End If
      Call VecScatterDestroy(ScatterSecToVecScal,ierr);CHKERRQ(ierr)
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
      End If
      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
         Call SectionRealDestroy(plasticStrainOldSec,ierr);CHKERRQ(ierr)
         Call VecScatterDestroy(ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)                   
      End If
      Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
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
      Type(VecScatter)                                   :: ScatterSecToVec,ScatterSecToVecScal,ScatterSecToVecCellMatS
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemScal
      Type(MEF90Element_Type)                            :: elemDisplacementType,elemScalType
      PetscScalar                                        :: myenergy,myenergySet

      !!! Create dof-based sections
      Call DMMeshGetSectionReal(MEF90DefMechCtx%DMVect,'default',xSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMVect,xSec,ScatterSecToVec,ierr);CHKERRQ(ierr)

      !!! Scatter data from Vec to Sec, or initialize
      Call SectionRealToVec(xSec,ScatterSecToVec,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr) 

      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMMatS,'default',plasticStrainSec,ierr);CHKERRQ(ierr)
         Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMMatS,plasticStrainSec,ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticStrainSec,ScatterSecToVecCellMatS,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)          
      Else
         PlasticStrainSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call DMMeshGetSectionReal(MEF90DefMechCtx%DMScal,'default',temperatureSec,ierr);CHKERRQ(ierr)
         Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMScal,temperatureSec,ScatterSecToVecScal,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(temperatureSec,ScatterSecToVecScal,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr)          
      Else
         temperatureSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%damage)) Then
         If (Associated(MEF90DefMechCtx%temperature)) Then
            Call SectionRealDuplicate(temperatureSec,damageSec,ierr);CHKERRQ(ierr)
         Else
            Call DMMeshGetSectionReal(MEF90DefMechCtx%DMScal,'default',damageSec,ierr);CHKERRQ(ierr)
            Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMScal,damageSec,ScatterSecToVecScal,ierr);CHKERRQ(ierr)
         End If
         Call SectionRealToVec(damageSec,ScatterSecToVecScal,SCATTER_REVERSE,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)          
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
         Case (MEF90DefMech_damageTypeAT1Elastic,MEF90DefMech_damageTypeAT2Elastic)
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
!!erwan->!!
         Case (MEF90DefMech_damageTypeATk)
            If (elemDisplacementType%coDim == 0) Then
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  QuadratureOrder = 5 * max(elemDisplacementType%order - 1, elemScalType%order)
               Else
                  QuadratureOrder = 5 ! * (elemDisplacementType%order - 1)
               End If
               Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
               Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
               Call MEF90ElasticityEnergySet(myenergy,xSec,plasticStrainSec,temperatureSec,MEF90DefMechCtx%DMVect,MEF90DefMechCtx%DMScal,setIS, &
                                             matpropSet%HookesLaw,matpropSet%LinearThermalExpansion,elemDisplacement,elemDisplacementType,elemScal,elemScalType,ierr)
               Call MEF90Element_Destroy(elemDisplacement,ierr)
               Call MEF90Element_Destroy(elemScal,ierr)
            End If
!!<--erwan!!
         Case (MEF90DefMech_damageTypeAT1,MEF90DefMech_damageTypeAT2)
            If (elemDisplacementType%coDim == 0) Then
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  QuadratureOrder = 2 * max(elemDisplacementType%order - 1, elemScalType%order) + 2 * elemScalType%order
               Else
                  QuadratureOrder = 2 * (elemDisplacementType%order - 1) + 2 * elemScalType%order
               End If
               Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
               Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
               Call MEF90GradDamageElasticEnergySet(myenergy,xSec,damageSec,plasticStrainSec,temperatureSec,MEF90DefMechCtx%DMVect,MEF90DefMechCtx%DMScal,setIS, &
                                             matPropSet%residualStiffness,matpropSet%HookesLaw,matpropSet%LinearThermalExpansion,elemDisplacement,elemDisplacementType,elemScal,elemScalType,ierr)
               Call MEF90Element_Destroy(elemDisplacement,ierr)
               Call MEF90Element_Destroy(elemScal,ierr)
            End If
         Case default
            Print*,__FUNCT__,': Unimplemented damage Type',cellSetOptions%damageType
            STOP  
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
      Call VecScatterDestroy(ScatterSecToVecScal,ierr);CHKERRQ(ierr)
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
      End If
      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
         Call VecScatterDestroy(ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)                   
      End If
      Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
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
      Type(SectionReal)                                  :: plasticStrainSec,temperatureSec
      Type(IS)                                           :: CellSetGlobalIS,setIS,setISdof,bcIS
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(VecScatter)                                   :: ScatterSecToVec,ScatterSecToVecScal,ScatterSecToVecCellMatS
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemScal
      Type(MEF90Element_Type)                            :: elemDisplacementType,elemScalType

      Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMMatS,'default',stressSec,ierr);CHKERRQ(ierr)
      Call SectionRealSet(stressSec,0.0_Kr,ierr);CHKERRQ(ierr)
      
      Call DMMeshGetSectionReal(MEF90DefMechCtx%DMVect,'default',xSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMVect,xSec,ScatterSecToVec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(xSec,ScatterSecToVec,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr) 
      
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMMatS,stressSec,ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)

      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMMatS,'default',plasticStrainSec,ierr);CHKERRQ(ierr)
         Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMMatS,plasticStrainSec,ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticStrainSec,ScatterSecToVecCellMatS,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)          
      Else
         PlasticStrainSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call DMMeshGetSectionReal(MEF90DefMechCtx%DMScal,'default',temperatureSec,ierr);CHKERRQ(ierr)
         Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMScal,temperatureSec,ScatterSecToVecScal,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(temperatureSec,ScatterSecToVecScal,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr)          
      Else
         temperatureSec%v = 0
      End If

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
!!erwan-->!!
         Case (MEF90DefMech_damageTypeATk)
            If (elemDisplacementType%coDim == 0) Then
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  QuadratureOrder = max(elemDisplacementType%order - 1, elemScalType%order)
               Else
                  QuadratureOrder = 5 !*elemDisplacementType%order - 1
               End If
               Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
               Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
               Call MEF90ElasticityStressSet(stressSec,xSec,plasticStrainSec,temperatureSec,MEF90DefMechCtx%DMVect,MEF90DefMechCtx%DMScal,setIS, &
                                             matpropSet%HookesLaw,matpropSet%LinearThermalExpansion,elemDisplacement,elemDisplacementType,elemScal,elemScalType,ierr)
               Call MEF90Element_Destroy(elemDisplacement,ierr)
               Call MEF90Element_Destroy(elemScal,ierr)
            End If
!!<--erwan!!
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
               Call MEF90ElasticityStressSet(stressSec,xSec,plasticStrainSec,temperatureSec,MEF90DefMechCtx%DMVect,MEF90DefMechCtx%DMScal,setIS, &
                                             matpropSet%HookesLaw,matpropSet%LinearThermalExpansion,elemDisplacement,elemDisplacementType,elemScal,elemScalType,ierr)
               Call MEF90Element_Destroy(elemDisplacement,ierr)
               Call MEF90Element_Destroy(elemScal,ierr)
            End If
         Case default
            Print*,__FUNCT__,': Unimplemented damage Type',cellSetOptions%damageType
            STOP  
         End Select
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do ! set
      !!! No need to complete since stressSec is cell-based
      Call SectionRealToVec(stressSec,ScatterSecToVecCellMatS,SCATTER_FORWARD,stress,ierr);CHKERRQ(ierr)          
      
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
         Call VecScatterDestroy(ScatterSecToVecScal,ierr);CHKERRQ(ierr)
      End If
      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
         Call VecScatterDestroy(ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)                   
      End If
      Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
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

   Subroutine MEF90DefMechBilinearFormDamageAT1Loc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      TYPE(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: elasticEnergyDensityGauss,PlasticDissipationDensityGauss,temperatureGauss
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
         elasticEnergyDensityGauss = (matprop%HookesLaw * inelasticStrainGauss) .DotP. inelasticStrainGauss

         !!! This is really twice the plastic energy density
         !!! This is not the PlasticDissipationDensityGauss properly, we should write 2*(matprop%HookesLaw * inelasticStrainGauss .DotP. (plasticStrainCell - plasticStrainoldCell ) )
         PlasticDissipationDensityGauss = 2.0_Kr*(matprop%HookesLaw * inelasticStrainGauss .DotP. plasticStrainCell)
         !!! This is really twice the elastic energy density
         Do iDoF1 = 1,numDofDamage
            Do iDoF2 = 1,numDofDamage
               Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) * ( &
                                      (elasticEnergyDensityGauss + PlasticDissipationDensityGauss) * elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss) + &
                                      C2 * (elemDamage%Grad_BF(iDoF1,iGauss) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
            End Do
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDamageAT1Loc

!!erwan_bilinear-->!!

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamageATkLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamageATkLoc: Assembles the bilinear form for ATk, that is:
!!!    <K(alpha) beta , delta> = \int a''(\alpha) W(e(u)) \beta \delta 
!!!                            + Gc/pi/\ell w''(alpha) \beta \delta
!!!                            + 2 Gc/pi*\ell \nabla \beta \nabla \delta
!!! with a(\alpha) = (1-w(\alpha)) / ( 1 + (k-1) w(\alpha)) and w(\alpha) = 1 - (1-\alpha)^2
!!! i.e. a''(\alpha) = 2 k * ( k + 3(k-1) v^2) / [k + (1-k) v^2]^3, w''(\alpha) = -2
!!!  
!!!  (c) 2014-2015 Erwan Tanne erwan.tanne@gmail.com, Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDamageATkLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      TYPE(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: elasticEnergyDensityGauss,PlasticDissipationDensityGauss,temperatureGauss,damageGauss
      Type(MEF90_MATS)                                   :: inelasticStrainGauss
      PetscReal                                          :: C2,C1
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr
      PetscReal                                          :: k,C0

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)
      
      k = matprop%k_for_ATk
      C1 = 2.0_Kr * matprop%fractureToughness / matprop%internalLength / PETSC_PI
      C2 = 2.0_Kr * matprop%fractureToughness * matprop%internalLength / PETSC_PI
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

         !!! This is really twice the plastic energy density
         !!! This is not the PlasticDissipationDensityGauss properly, we should write 2*(matprop%HookesLaw * inelasticStrainGauss .DotP. (plasticStrainCell - plasticStrainoldCell ) )
         PlasticDissipationDensityGauss = 2.0_Kr*(matprop%HookesLaw * inelasticStrainGauss .DotP. plasticStrainCell)

         Do iDoF1 = 1,numDofDamage
            Do iDoF2 = 1,numDofDamage
               Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) * ( & 
                                  (elasticEnergyDensityGauss * C0 - C1 + PlasticDissipationDensityGauss) * elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss) &
                              + C2 * (elemDamage%Grad_BF(iDoF1,iGauss) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
            End Do
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDamageATkLoc
!!<--erwan!!

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamageAT1ElasticLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamageAT1ElasticLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDamageAT1ElasticLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      TYPE(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
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
#define __FUNCT__ "MEF90DefMechBilinearFormDamageAT1UnilateralHDLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamageAT1UnilateralHDLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDamageAT1UnilateralHDLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      TYPE(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
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
            elasticEnergyDensityGauss = (matprop%HookesLaw * inelasticStrainGauss) .DotP. DeviatoricPart(inelasticStrainGauss)
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
   End Subroutine MEF90DefMechBilinearFormDamageAT1UnilateralHDLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamageAT1UnilateralDeviatoricLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamageAT1UnilateralDeviatoricLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDamageAT1UnilateralDeviatoricLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      TYPE(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
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

   Subroutine MEF90DefMechBilinearFormDamageAT1UnilateralPHDLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      TYPE(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
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
#define __FUNCT__ "MEF90DefMechBilinearFormDamageAT2Loc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamageAT2Loc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDamageAT2Loc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      TYPE(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
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
                                       C2 * (elemDamage%Grad_BF(iDoF1,iGauss) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
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

   Subroutine MEF90DefMechBilinearFormDamageAT2ElasticLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      TYPE(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
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
                                       C2 * (elemDamage%Grad_BF(iDoF1,iGauss) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
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

   Subroutine MEF90DefMechBilinearFormDamageAT2UnilateralHDLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      TYPE(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
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
            elasticEnergyDensityGauss = (matprop%HookesLaw * inelasticStrainGauss) .DotP. DeviatoricPart(inelasticStrainGauss)
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
   End Subroutine MEF90DefMechBilinearFormDamageAT2UnilateralHDLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamageAT2UnilateralDeviatoricLoc"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamageAT2UnilateralDeviatoricLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDamageAT2UnilateralDeviatoricLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      TYPE(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
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

   Subroutine MEF90DefMechBilinearFormDamageAT2UnilateralPHDLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      TYPE(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
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

   Subroutine MEF90DefMechOperatorDamageAT1Loc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      TYPE(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: elasticEnergyDensityGauss,PlasticDissipationDensityGauss,temperatureGauss,damageGauss
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
         elasticEnergyDensityGauss = (matprop%HookesLaw * inelasticStrainGauss) .DotP. inelasticStrainGauss
         
         !!! Be careful It is not true in all cases, only for monotone increasing plastic strain
         !!! PlasticDissipationDensityGauss = 2*(matprop%HookesLaw * inelasticStrainGauss) .DotP. plasticStrainCell
         PlasticDissipationDensityGauss = 2.0_Kr* ((matprop%HookesLaw * inelasticStrainGauss) .DotP. plasticStrainCell)
         !!! This is really twice the elastic energy density

         damageGauss = 0.0_Kr
         gradientDamageGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
            gradientDamageGauss = gradientDamageGauss + elemDamage%Grad_BF(iDoF1,iGauss) * xDof(iDoF1)
         End Do

         Do iDoF2 = 1,numDofDamage
            residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * ( &
                                  ( ( elasticEnergyDensityGauss + PlasticDissipationDensityGauss ) * (damageGauss - 1.0_Kr) + C1) * elemDamage%BF(iDoF2,iGauss) + &
                                  C2 * (gradientDamageGauss .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
         End Do

      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDamageAT1Loc


!!erwan_operator-->!!
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamageATkLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDamageATkLoc: Assembles the operator for ATk:
!!!    F(\alpha)\beta = \int [k(\alpha-1)W(e(u))] / [k-(k-1)(1-\alpha)^2]^2 \beta 
!!!                   + 2Gc/pi (1-\alpha) \beta / \ell 
!!!                   + 2Gc/pi \nabla \alpha \nabla \beta * \ell
!!!  
!!!  (c) 2014-2015 Erwan Tanne erwan.tanne@gmail.com, Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDamageATkLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      TYPE(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      PetscReal                                          :: elasticEnergyDensityGauss,PlasticDissipationDensityGauss,temperatureGauss,damageGauss
      Type(MEF90_MATS)                                   :: inelasticStrainGauss
      Type(MEF90_VECT)                                   :: gradientDamageGauss
      PetscReal                                          :: C1,C2,k,C0Operator
      PetscLogDouble                                     :: flops
      PetscErrorCode                                     :: ierr


      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)
      
      C1 = 2.0_Kr * matprop%fractureToughness / matprop%internalLength / PETSC_PI
      C2 = 2.0_Kr * matprop%fractureToughness * matprop%internalLength / PETSC_PI
      k = matprop%k_for_ATk
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


         !!! PlasticDissipationDensityGauss = 2*(matprop%HookesLaw * inelasticStrainGauss) .DotP. plasticStrainCell
         PlasticDissipationDensityGauss = 2.0_Kr* ((matprop%HookesLaw * inelasticStrainGauss) .DotP. plasticStrainCell)
         !!! This is really twice the plastic energy density

         damageGauss = 0.0_Kr
         gradientDamageGauss = 0.0_Kr
         Do iDoF1 = 1,numDofDamage
            damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
            gradientDamageGauss = gradientDamageGauss + elemDamage%Grad_BF(iDoF1,iGauss) * xDof(iDoF1)
         End Do
         C0Operator = k * (damageGauss - 1.0_Kr) * elasticEnergyDensityGauss / &
                     (k + (1.0_Kr - k) * (1.0_kr - damageGauss)**2)**2

         Do iDoF2 = 1,numDofDamage
            residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * (                       &
                                  (C0Operator + C1 * (1.0_Kr - damageGauss) + PlasticDissipationDensityGauss * (damageGauss - 1.0_Kr) ) * elemDamage%BF(iDoF2,iGauss) & 
                                  + C2 * (gradientDamageGauss .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDamageATkLoc
!!<--erwan!!

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamageAT1ElasticLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDamageAT1ElasticLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDamageAT1ElasticLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      TYPE(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

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
                                  C2 * (gradientDamageGauss .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDamageAT1ElasticLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamageAT1UnilateralHDLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDamageAT1UnilateralHDLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDamageAT1UnilateralHDLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      TYPE(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

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
            elasticEnergyDensityGauss = (matprop%HookesLaw * inelasticStrainGauss) .DotP. DeviatoricPart(inelasticStrainGauss)
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
   End Subroutine MEF90DefMechOperatorDamageAT1UnilateralHDLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamageAT1UnilateralDeviatoricLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDamageAT1UnilateralDeviatoricLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDamageAT1UnilateralDeviatoricLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      TYPE(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

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

   Subroutine MEF90DefMechOperatorDamageAT1UnilateralPHDLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      TYPE(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

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
#define __FUNCT__ "MEF90DefMechOperatorDamageAT2Loc"
!!!
!!!  
!!!  MEF90DefMechOperatorDamageAT2Loc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDamageAT2Loc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      TYPE(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

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
                                  C2 * (gradientDamageGauss .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
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
   Subroutine MEF90DefMechOperatorDamageAT2ElasticLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      TYPE(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

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
                                  C2 * (gradientDamageGauss .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
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
   Subroutine MEF90DefMechOperatorDamageAT2UnilateralHDLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      TYPE(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

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
            elasticEnergyDensityGauss = (matprop%HookesLaw * inelasticStrainGauss) .DotP. DeviatoricPart(inelasticStrainGauss)
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
   End Subroutine MEF90DefMechOperatorDamageAT2UnilateralHDLoc

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamageAT2UnilateralDeviatoricLoc"
!!!
!!!  
!!!  MEF90DefMechOperatorDamageAT2UnilateralDeviatoricLoc:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90DefMechOperatorDamageAT2UnilateralDeviatoricLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      TYPE(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

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

   Subroutine MEF90DefMechOperatorDamageAT2UnilateralPHDLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:),Pointer                     :: residualLoc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
      TYPE(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

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
#define __FUNCT__ "MEF90DefMechOperatorDamage"
!!!
!!!  
!!!  MEF90DefMechOperatorDamage: Build the operator. When called in SNES, the solution time should always match the target time, 
!!!                                    so there is no need for interpolation of the forcees, external, and boundary values
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechOperatorDamage(snesDamage,damage,residual,MEF90DefMechCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesDamage
      Type(Vec),Intent(IN)                               :: damage
      Type(Vec),Intent(INOUT)                            :: residual
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(DM)                                           :: mesh
      Type(SectionReal)                                  :: displacementSec,residualSec
      Type(SectionReal)                                  :: plasticStrainSec,temperatureSec,damageSec
      Type(SectionReal)                                  :: boundaryDamageSec,forceSec,pressureForceSec
      PetscReal,Dimension(:),Pointer                     :: damagePtr,residualPtr
      PetscReal,Dimension(:),Pointer                     :: boundaryDamagePtr
      PetscReal,Dimension(:),Pointer                     :: displacementDof,damageDof,temperatureDof,residualLoc
      PetscReal,Dimension(:),Pointer                     :: plasticStrainLoc
      Type(MEF90_MATS)                                   :: plasticStrainCell
      Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof,bcIS
      PetscInt,Dimension(:),Pointer                      :: setID,cellID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder,cell
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90DefMechVertexSetOptions_Type),Pointer    :: vertexSetOptions
      Type(VecScatter)                                   :: ScatterSecToVec,ScatterSecToVecScal
      Type(VecScatter)                                   :: ScatterSecToVecCellMatS
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemDamage
      Type(MEF90Element_Type)                            :: elemDisplacementType,elemDamageType
      PetscInt                                           :: p,dof,c,dim,numDof
      PetscInt                                           :: nVal
      PetscReal,Dimension(:),Pointer                     :: nullPtr
      
      Procedure(MEF90DefMechOperatorLoc),pointer         :: localOperatorFunction
      
      localOperatorFunction =>MEF90DefMechOperatorNull
      Call SNESGetDM(snesDamage,mesh,ierr);CHKERRQ(ierr)

      !!! Get Section and scatter associated with each field
      Call DMMeshGetSectionReal(MEF90DefMechCtx%DMVect,'default',displacementSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMVect,displacementSec,ScatterSecToVec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(displacementSec,ScatterSecToVec,SCATTER_REVERSE,MEF90DefMechCtx%displacement,ierr);CHKERRQ(ierr) 

      Call DMMeshGetSectionReal(MEF90DefMechCtx%DMScal,'default',damageSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMScal,damageSec,ScatterSecToVecScal,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(damageSec,ScatterSecToVecScal,SCATTER_REVERSE,damage,ierr);CHKERRQ(ierr)   
      Call SectionRealDuplicate(damageSec,boundaryDamageSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(boundaryDamageSec,ScatterSecToVecScal,SCATTER_REVERSE,MEF90DefMechCtx%boundaryDamage,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(damageSec,residualSec,ierr);CHKERRQ(ierr)
      Call SectionRealSet(residualSec,0.0_Kr,ierr);CHKERRQ(ierr)
      Call VecSet(residual,0.0_Kr,ierr);CHKERRQ(ierr)

      Nullify(TemperatureDof)
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDuplicate(damageSec,temperatureSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(temperatureSec,ScatterSecToVecScal,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr)   
      Else
         temperatureSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMMatS,'default',plasticStrainSec,ierr);CHKERRQ(ierr)
         Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMMatS,plasticStrainSec,ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticStrainSec,ScatterSecToVecCellMatS,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)          
      Else
         PlasticStrainSec%v = 0
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
         If ((Size(cellID) > 0) .AND. (elemDamageType%coDim == 0)) Then
            Select Case (cellSetOptions%damageType)
            Case (MEF90DefMech_damageTypeAT1Elastic)
               QuadratureOrder = 2 * (elemDamageType%order - 1)
               localOperatorFunction => MEF90DefMechOperatorDamageAT1ElasticLoc
            Case (MEF90DefMech_damageTypeAT2Elastic)
               QuadratureOrder = 2 * elemDamageType%order
               localOperatorFunction => MEF90DefMechOperatorDamageAT2ElasticLoc
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
               End Select
!!erwan-->!!
            Case (MEF90DefMech_damageTypeATk)
               QuadratureOrder = 9! 5 * elemDamageType%order + 5 * (elemDisplacementType%order - 1)
               localOperatorFunction => MEF90DefMechOperatorDamageATkLoc
!!<--erwan!!
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
               Case (MEF90DefMech_unilateralContactTypePrincipalStrains)
                  localOperatorFunction => MEF90DefMechOperatorNull
               End Select
            End Select

            !!! Allocate storage for fields at dof and Gauss points
            !!! Leaving plasticstrain aside until I can remember how it is supposed to be dealt with
            Allocate(residualLoc(ElemDamageType%numDof))
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
                  plasticStrainCell = plasticStrainLoc
               Else
                  plasticStrainCell = 0.0_Kr
               End If
                  
               residualLoc = 0.0_Kr
               Call localOperatorFunction   (residualLoc,damageDof,displacementDof,nullPtr,nullPtr,temperatureDof,plasticStrainCell,matpropSet,elemDisplacement(cell),elemDamage(cell))
               Call SectionRealUpdateClosure(residualSec,MEF90DefMechCtx%DMScal,cellID(cell),residualLoc,ADD_VALUES,ierr);CHKERRQ(iErr)

               If (Associated(MEF90DefMechCtx%plasticStrain)) Then
                  Call SectionRealRestore(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
               End If
            End Do

            DeAllocate(displacementDof)
            DeAllocate(damageDof)
            If (Associated(MEF90DefMechCtx%temperature)) Then
               DeAllocate(temperatureDof)
            End If
            DeAllocate(residualLoc)

            Call MEF90Element_Destroy(elemDisplacement,ierr)
            Call MEF90Element_Destroy(elemDamage,ierr)
            Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If 
      End Do ! set

      !!! "Ghost update" for the residual Section
      Call SectionRealComplete(residualSec,ierr);CHKERRQ(ierr)
      !!! Scatter back from SectionReal to Vec
      Call SectionRealToVec(residualSec,ScatterSecToVecScal,SCATTER_FORWARD,residual,ierr);CHKERRQ(ierr)
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
         Call VecScatterDestroy(ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)      
      End If
      
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
      End If

      Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(boundaryDamageSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(residualSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(displacementSec,ierr);CHKERRQ(ierr)

      Call VecScatterDestroy(ScatterSecToVecScal,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDamage

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamage"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamage:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDamage(snesDamage,damage,A,M,flg,MEF90DefMechCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesDamage
      Type(Vec),Intent(IN)                               :: damage
      Type(Mat),Intent(INOUT)                            :: A,M
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr  
         
      Type(SectionReal)                                  :: displacementSec,damageSec,temperatureSec,plasticStrainSec
      PetscReal,Dimension(:),Pointer                     :: displacementDOF,damageDof,temperatureDof,nullPtr
      Type(VecScatter)                                   :: ScatterSecToVec,ScatterSecToVecScal,ScatterSecToVecCellMatS
      PetscReal,Dimension(:),Pointer                     :: plasticStrainLoc
      Type(MEF90_MATS)                                   :: plasticStrainCell
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
      
      localAssemblyFunction =>MEF90DefMechBilinearFormNull      
      Call MatZeroEntries(A,ierr);CHKERRQ(ierr)
      Call SNESGetDM(snesDamage,mesh,ierr);CHKERRQ(ierr)
      
      !!! Get Section and scatter associated with each field
      Call DMMeshGetSectionReal(MEF90DefMechCtx%DMVect,'default',displacementSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMVect,displacementSec,ScatterSecToVec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(displacementSec,ScatterSecToVec,SCATTER_REVERSE,MEF90DefMechCtx%displacement,ierr);CHKERRQ(ierr)   

      Call DMMeshGetSectionReal(MEF90DefMechCtx%DMScal,'default',damageSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMScal,damageSec,ScatterSecToVecScal,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(damageSec,ScatterSecToVecScal,SCATTER_REVERSE,damage,ierr);CHKERRQ(ierr)   

      Nullify(temperatureDof)
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDuplicate(damageSec,temperatureSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(temperatureSec,ScatterSecToVecScal,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr)   
      Else
         temperatureSec%v = 0
      End If
      
      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMMatS,'default',plasticStrainSec,ierr);CHKERRQ(ierr)
         Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMMatS,plasticStrainSec,ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticStrainSec,ScatterSecToVecCellMatS,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)          
      Else
         PlasticStrainSec%v = 0
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
!!erwan-->!!
            Case (MEF90DefMech_damageTypeATk)
               QuadratureOrder = 9 !* (elemDisplacementType%order - 1) + 5 * elemDamageType%order
               localAssemblyFunction => MEF90DefMechBilinearFormDamageATkLoc
!!<--erwan!!
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
               Case (MEF90DefMech_unilateralContactTypePrincipalStrains)
                  localAssemblyFunction => MEF90DefMechBilinearFormNull
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
                  plasticStrainCell = plasticStrainLoc
               Else
                 plasticStrainCell = 0.0_Kr
               End If
                  
               Aloc = 0.0_Kr
               Call localAssemblyFunction(Aloc,damageDof,displacementDof,nullPtr,temperatureDof,plasticStrainCell,matpropSet,elemDisplacement(cell),elemDamage(cell))
               Call DMmeshAssembleMatrix(A,mesh,damageSec,cellID(cell),ALoc,ADD_VALUES,ierr);CHKERRQ(ierr)

               If (Associated(MEF90DefMechCtx%plasticStrain)) Then
                  Call SectionRealRestore(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
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
            Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If 
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

      Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVecScal,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(displacementSec,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)
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
   
      Type(SectionReal)                                  :: alphaSec
      Type(IS)                                           :: CellSetGlobalIS,setIS,setISdof,bcIS
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(VecScatter)                                   :: ScatterSecToVecScal
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemScal
      Type(MEF90Element_Type)                            :: elemScalType
      PetscReal                                          :: myenergy,myenergySet

      !!! Create dof-based sections
      Call DMMeshGetSectionReal(MEF90DefMechCtx%DMScal,'default',alphaSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMScal,alphaSec,ScatterSecToVecScal,ierr);CHKERRQ(ierr)

      !!! Scatter data from Vec to Sec, or initialize
      Call SectionRealToVec(alphaSec,ScatterSecToVecScal,SCATTER_REVERSE,alpha,ierr);CHKERRQ(ierr) 

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
         Case (MEF90DefMech_damageTypeAT1)
            If (elemScalType%coDim == 0) Then
               QuadratureOrder = 2 * (elemScalType%order - 1)
               Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
               Call MEF90GradDamageSurfaceEnergySetAT1(myenergy,alphaSec,MEF90DefMechCtx%DMScal,setIS,matpropSet%internalLength,matpropSet%fractureToughness,elemScal,elemScalType,ierr)
               Call MEF90Element_Destroy(elemScal,ierr)
            End If
         Case (MEF90DefMech_damageTypeAT2)
            If (elemScalType%coDim == 0) Then
               QuadratureOrder = 2*elemScalType%order
               Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
               Call MEF90GradDamageSurfaceEnergySetAT2(myenergy,alphaSec,MEF90DefMechCtx%DMScal,setIS,matpropSet%internalLength,matpropSet%fractureToughness,elemScal,elemScalType,ierr)
               Call MEF90Element_Destroy(elemScal,ierr)
            End If
!!erwan-->!!
         Case (MEF90DefMech_damageTypeATk)
            If (elemScalType%coDim == 0) Then
               ! Which QuadratureOdre chose for ATk
               QuadratureOrder = 5 !* (elemScalType%order - 1)
               Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
               Call MEF90GradDamageSurfaceEnergySetATk(myenergy,alphaSec,MEF90DefMechCtx%DMScal,setIS,matpropSet%internalLength,matpropSet%fractureToughness,elemScal,elemScalType,ierr)
               Call MEF90Element_Destroy(elemScal,ierr)
            End If
!!<--erwan!!
         End Select
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         Call MPI_AllReduce(myenergy,myenergySet,1,MPIU_SCALAR,MPI_SUM,MEF90DefMechCtx%MEF90Ctx%comm,ierr);CHKERRQ(ierr)
         energy(set) = energy(set) + myenergySet
      End Do ! set
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr) 

      Call VecScatterDestroy(ScatterSecToVecScal,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(alphaSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechSurfaceEnergy

End Module MEF90_APPEND(m_MEF90_DefMechAssembly,MEF90_DIM)D
