#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
module m_MEF90_DefMechDrivingForceDruckerPrager
#include "finclude/petscdef.h"
   Use m_MEF90
   Use m_MEF90_DefMechDrivingForceDruckerPrager_class
   implicit none

   Type, extends(MEF90_DefMechDrivingForce_Type)       :: MEF90_DefMechDrivingForceDruckerPrager_Type
   Contains
      Procedure, pass(self)                            :: df  => dfDruckerPrager
      Procedure, pass(self)                            :: Ddf => dfDruckerPrager
   end Type

   interface MEF90_DefMechDrivingForceDruckerPrager_Type
      module procedure MEF90_DefMechDrivingForceDruckerPrager_Constructor
   end interface

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_DefMechDrivingForceDruckerPrager_Type"
!!!
!!!  
!!!  MEF90_DefMechDrivingForceDruckerPrager_Type: the default constructor for a MEF90_DefMechAT2_Type
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   Type(MEF90_DefMechAT2_Type) Function MEF90_DefMechDrivingForceDruckerPrager_Type()
      MEF90_DefMechDrivingForceDruckerPrager_Type%type              = 'MEF90_DefMechDrivingForceTypeDruckerPrager'
   End Function MEF90_DefMechDrivingForceDruckerPrager_Type

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamageDrivingForceAT1DruckerPrager"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamageDrivingForceAT1DruckerPrager: 
!!!               The nucleation driving force derived in Kumar, Bourdin, Francfort, Lopez-Pamies
!!!               Revisiting Nucleation in the Phase-Field Approach to Brittle Fracture
!!!               Initial version without corrections:
!!!               de = 1/(1+\beta_3 I_1^2) (\beta_0 + \beta_1 I_1 + \beta_2 \sqrt(J_2))
!!!  
!!!  (c) 2019 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDamageDrivingForceAT1DruckerPrager(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDamage,numDofDisplacement,numGauss
      Type(MEF90_MATS)                                   :: AeGauss
      Type(MatS3D)                                       :: Ae3DGauss
      PetscReal                                          :: I1Gauss,sqrtJ2Gauss,DalphaI1Gauss,DalphasqrtJ2Gauss
      PetscReal                                          :: Dbeta1,Dbeta2,D
      PetscReal                                          :: beta0,beta1,beta2,beta3,DalphadeGauss
      PetscReal                                          :: E,Gc,delta,gamma,ell,sigma_cs,sigma_ts
      PetscReal                                          :: damageGauss,stiffnessMultGauss,DalphaStiffnessMultGauss
      PetscLogDouble                                     :: flops
      Character(len=MEF90MXSTRLEN)                      :: IOBuffer
      PetscErrorCode                                     :: ierr

      If (matprop%HookesLaw%Type /= MEF90HookesLawTypeIsotropic) Then
         Write(IOBuffer,'(''ERROR: unimplemented Hookes law type in'' (A),''\n'')') __FUNCT__
         SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,IOBuffer,ierr)
      End If

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)

      E        = matprop%HookesLaw%YoungsModulus
      Gc       = matprop%FractureToughness
      ell      = matprop%internalLength
      delta    = matprop%drivingForceDelta
      gamma    = matprop%drivingForceGamma
      sigma_cs = matprop%drivingForceCompressiveStrength
      sigma_ts = matprop%drivingForceTensileStrength

      beta3  = gamma * ell  
      beta0  = 3.0_Kr * Gc / 8.0_Kr / ell * delta

      D      = 2.0_Kr * sigma_cs * sigma_ts
      Dbeta1 = -Gc * (1.+delta) / 8.0_Kr * (sigma_cs - sigma_ts)
      Dbeta2 = -Gc * (1.+delta) / 8.0_Kr * (sigma_cs + sigma_ts)
      beta1  = Dbeta1 / D
      beta2  = Dbeta2 / D

      Do iGauss = 1,numGauss
         !!! Compute the value of the damage field at the Gauss points
         damageGauss = 0.0_Kr
         Do iDof1 = 1,numDofDamage
            damageGauss = damageGauss + elemDamage%BF(iDof1,iGauss) * xDof(iDof1)
         End Do
         stiffnessMultGauss       = MEF90aAT(damageGauss) + matprop%residualStiffness
         DalphaStiffnessMultGauss = MEF90DaAT(damageGauss)
         AeGauss = 0.0_Kr
         Do iDof1 = 1,numDofDisplacement
            AeGauss = AeGauss + elemDisplacement%GradS_BF(iDof1,iGauss) * displacementDof(iDof1)
         End Do
         AeGauss           = matprop%HookesLaw * AeGauss
#if MEF90_DIM == 3
         I1Gauss           = Moment(1,AeGauss) * stiffnessMultGauss
         DalphaI1Gauss     = Moment(1,AeGauss) * DalphaStiffnessMultGauss
         sqrtJ2Gauss       = sqrt(Moment(2,DeviatoricPart(AeGauss))) * stiffnessMultGauss
         DalphaSqrtJ2Gauss = sqrt(Moment(2,DeviatoricPart(AeGauss))) * DalphaStiffnessMultGauss
#else
         If (matProp%HookesLaw%isPlaneStress) Then
            Ae3DGauss    = 0.0_Kr
            Ae3DGauss%XX = AeGauss%XX
            Ae3DGauss%YY = AeGauss%YY
            Ae3DGauss%XY = AeGauss%XY
         Else 
            Ae3DGauss    = 0.0_Kr
            Ae3DGauss%XX = AeGauss%XX
            Ae3DGauss%YY = AeGauss%YY
            Ae3DGauss%XY = AeGauss%XY
            Ae3DGauss%ZZ = matprop%HookesLaw%lambda / (matprop%HookesLaw%lambda + matprop%HookesLaw%mu) * 0.5_Kr * &
                           (AeGauss%XX + AeGauss%YY)
         End If
         I1Gauss           = Moment(1,Ae3DGauss) * stiffnessMultGauss
         DalphaI1Gauss     = Moment(1,Ae3DGauss) * DalphaStiffnessMultGauss
         sqrtJ2Gauss       = sqrt(Moment(2,DeviatoricPart(Ae3DGauss))) * stiffnessMultGauss
         DalphaSqrtJ2Gauss = sqrt(Moment(2,DeviatoricPart(Ae3DGauss))) * DalphaStiffnessMultGauss
#endif
         DalphadeGauss = -(beta0 + beta1 * I1Gauss + beta2 * sqrtJ2Gauss) * 2.0_Kr * beta3 * I1Gauss * DalphaI1Gauss / (1.0_Kr + beta3 * I1Gauss**2)**2 & 
                         +(beta1 * DalphaI1Gauss + beta2 * DalphaSqrtJ2Gauss)  / (1.0_Kr + beta3 * I1Gauss**2)
         Do iDoF1 = 1,numDofDamage
            Do iDoF2 = 1,numDofDamage
               Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) * DalphadeGauss * &
                                        elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss)
            End Do
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDamageDrivingForceAT1DruckerPrager

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamageDrivingForceAT1DruckerPrager2"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamageDrivingForceAT1DruckerPrager2: 
!!!               The nucleation driving force derived in Kumar, Bourdin, Francfort, Lopez-Pamies
!!!               Revisiting Nucleation in the Phase-Field Approach to Brittle Fracture
!!!               Version with corrections:
!!!               de = 1/(1+\beta_3 I_1^2) (\beta_0 + \beta_1 I_1 + \beta_2 \sqrt(J_2))
!!!  
!!!  (c) 2019 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechBilinearFormDamageDrivingForceAT1DruckerPrager2(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
      Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
      PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
      Type(MEF90_MATPROP),Intent(IN)                     :: matprop
      Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDamage,numDofDisplacement,numGauss
      Type(MEF90_MATS)                                   :: AeGauss
      Type(MatS3D)                                       :: Ae3DGauss
      PetscReal                                          :: I1Gauss,sqrtJ2Gauss,DalphaI1Gauss,DalphasqrtJ2Gauss
      PetscReal                                          :: beta0,beta1,beta2,beta3,DalphadeGauss
      PetscReal                                          :: E,Gc,delta,gamma,ell,sigma_cs,sigma_ts
      PetscReal                                          :: Dbeta1,Dbeta2,D
      PetscReal                                          :: damageGauss,stiffnessMultGauss,DalphaStiffnessMultGauss
      PetscLogDouble                                     :: flops
      Character(len=MEF90MXSTRLEN)                      :: IOBuffer
      PetscErrorCode                                     :: ierr

      If (matprop%HookesLaw%Type /= MEF90HookesLawTypeIsotropic) Then
         Write(IOBuffer,'(''ERROR: unimplemented Hookes law type in'' (A),''\n'')') __FUNCT__
         SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,IOBuffer,ierr)
      End If

      numDofDisplacement = size(elemDisplacement%BF,1)
      numDofDamage = size(elemDamage%BF,1)
      numGauss = size(elemDamage%BF,2)

      E        = matprop%HookesLaw%YoungsModulus
      Gc       = matprop%FractureToughness
      ell      = matprop%internalLength
      delta    = matprop%drivingForceDelta
      gamma    = matprop%drivingForceGamma
      sigma_cs = matprop%drivingForceCompressiveStrength
      sigma_ts = matprop%drivingForceTensileStrength

      beta3  = gamma * ell  
      beta0  = 3.0_Kr * Gc / 8.0_Kr / ell * delta

      D      = 2.0_Kr * sigma_ts * sigma_cs / sqrt(3.0_Kr)
      Dbeta1 = sigma_cs / sqrt(3.0_Kr) * ((1.0_Kr + gamma * ell  * sigma_ts**2) * (sigma_ts**2/E - 3.0_Kr * Gc / 8.0_Kr / ell) - 3.0_Kr * Gc * delta / 8.0_Kr / ell) & 
             - sigma_ts / sqrt(3.0_Kr) * ((1.0_Kr + gamma * ell  * sigma_cs**2) * (sigma_cs**2/E - 3.0_Kr * Gc / 8.0_Kr / ell) - 3.0_Kr * Gc * delta / 8.0_Kr / ell)
      Dbeta2 = sigma_ts * ((1.0_Kr+gamma * ell  * sigma_cs**2) * (sigma_cs**2/E - 3.0_Kr * Gc / 8.0_Kr / ell) - 3.0_Kr * Gc * delta / 8.0_Kr / ell) &
             + sigma_cs * ((1.0_Kr+gamma * ell  * sigma_ts**2) * (sigma_ts**2/E - 3.0_Kr * Gc / 8.0_Kr / ell) - 3.0_Kr * Gc * delta / 8.0_Kr / ell)
      beta1 = Dbeta1 / D
      beta2 = Dbeta2 / D

      Do iGauss = 1,numGauss
         !!! Compute the value of the damage field at the Gauss points
         damageGauss = 0.0_Kr
         Do iDof1 = 1,numDofDamage
            damageGauss = damageGauss + elemDamage%BF(iDof1,iGauss) * xDof(iDof1)
         End Do
         stiffnessMultGauss       = MEF90aAT(damageGauss) + matprop%residualStiffness
         DalphaStiffnessMultGauss = MEF90DaAT(damageGauss)
         AeGauss = 0.0_Kr
         Do iDof1 = 1,numDofDisplacement
            AeGauss = AeGauss + elemDisplacement%GradS_BF(iDof1,iGauss) * displacementDof(iDof1)
         End Do
         AeGauss           = matprop%HookesLaw * AeGauss
#if MEF90_DIM == 3
         I1Gauss           = Moment(1,AeGauss) * stiffnessMultGauss
         DalphaI1Gauss     = Moment(1,AeGauss) * DalphaStiffnessMultGauss
         sqrtJ2Gauss       = sqrt(Moment(2,DeviatoricPart(AeGauss))) * stiffnessMultGauss
         DalphaSqrtJ2Gauss = sqrt(Moment(2,DeviatoricPart(AeGauss))) * DalphaStiffnessMultGauss
#else
         If (matProp%HookesLaw%isPlaneStress) Then
            Ae3DGauss    = 0.0_Kr
            Ae3DGauss%XX = AeGauss%XX
            Ae3DGauss%YY = AeGauss%YY
            Ae3DGauss%XY = AeGauss%XY
         Else 
            Ae3DGauss    = 0.0_Kr
            Ae3DGauss%XX = AeGauss%XX
            Ae3DGauss%YY = AeGauss%YY
            Ae3DGauss%XY = AeGauss%XY
            Ae3DGauss%ZZ = matprop%HookesLaw%lambda / (matprop%HookesLaw%lambda + matprop%HookesLaw%mu) * 0.5_Kr * &
                           (AeGauss%XX + AeGauss%YY)
         End If
         I1Gauss           = Moment(1,Ae3DGauss) * stiffnessMultGauss
         DalphaI1Gauss     = Moment(1,Ae3DGauss) * DalphaStiffnessMultGauss
         sqrtJ2Gauss       = sqrt(Moment(2,DeviatoricPart(Ae3DGauss))) * stiffnessMultGauss
         DalphaSqrtJ2Gauss = sqrt(Moment(2,DeviatoricPart(Ae3DGauss))) * DalphaStiffnessMultGauss
#endif
         DalphadeGauss = -(beta0 + beta1 * I1Gauss + beta2 * sqrtJ2Gauss) * 2.0_Kr * beta3 * I1Gauss * DalphaI1Gauss / (1.0_Kr + beta3 * I1Gauss**2)**2 & 
                         +(beta1 * DalphaI1Gauss + beta2 * DalphaSqrtJ2Gauss)  / (1.0_Kr + beta3 * I1Gauss**2)
         Do iDoF1 = 1,numDofDamage
            Do iDoF2 = 1,numDofDamage
               Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) * DalphadeGauss * &
                                        elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss)
            End Do
         End Do
      End Do
      !flops = 2 * numGauss * numDofDisplacement**2
      !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechBilinearFormDamageDrivingForceAT1DruckerPrager2




End module m_MEF90_DefMechDrivingForceDruckerPrager
