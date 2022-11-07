Module m_MEF90_Materials_Types
#include "petsc/finclude/petsc.h"
   Use m_MEF90_LinAlg
   Use m_MEF90_Ctx
   Use m_MEF90_DMPlex
   Use petsc
   IMPLICIT NONE

   Type MEF90HookesLaw2D
      Type(Tens4OS2D)    :: fullTensor
      Type(Tens4OS3D)    :: fullTensorLocal,fullTensor3D
      PetscReal          :: lambda,mu,YoungsModulus,PoissonRatio,BulkModulus
      PetscEnum          :: type
      PetscBool          :: isPlaneStress
   End Type MEF90HookesLaw2D

   Type MEF90HookesLaw3D
      Type(Tens4OS3D)    :: fullTensor,fullTensorLocal
      PetscReal          :: lambda,mu,YoungsModulus,PoissonRatio,BulkModulus
      PetscEnum          :: type
#if (PETSC_SIZEOF_INT == 4)
      ! With 4-byte integers, this declared type is 4-bytes shy of being
      ! 8-byte aligned, which can be problematic for arrays of this type.
      ! The MEF90HookesLaw2D has an extra PetscBool that keeps it 8-byte
      ! aligned, so we'll mimic that:
      PetscBool          :: padding = PETSC_FALSE
#endif
   End Type MEF90HookesLaw3D

   Type MEF90RotationMatrix3D
      ! The Bunge (passive) convention is used. The rotation matrix transforms a vector from the global frame to the local frame.
      ! - rotation of a vector:              [V_local]_i    = [R]_ij . [V_global]_j
      ! - rotation of a 2nd order tensor:    [M_local]_ij   = [R]_ik . [M_global]_kl . [R^T]_lj
      ! - rotation of a fourth order tensor: [C_local]_ijkl = [R]_ip . [R]_jq . [R]_kr . [R]_ls . [C_global]_pqrs 
      Type(MAT3D)        :: fullTensor
      PetscReal          :: phi1,Phi,phi2
      Type(Vect3D)       :: V1,V2,V3
      PetscBool          :: fromEuler
   End Type MEF90RotationMatrix3D

   Type MEF90MatProp2D_Type
      PetscReal                     :: Density                                          ! rho
      PetscReal                     :: FractureToughness                                ! Gc
      Type(MatS2D)                  :: toughnessAnisotropyMatrix                        ! 
      PetscReal                     :: SpecificHeat                                     ! Cp
      Type(MatS2D)                  :: ThermalConductivity                              ! K
      Type(MatS2D)                  :: LinearThermalExpansion                           ! alpha
      Type(MEF90HookesLaw2D)        :: HookesLaw                                        ! A
      PetscReal                     :: internalLength                                   ! l
      PetscReal                     :: residualStiffness                                ! eta
      PetscReal                     :: yieldStress                                      ! yield stress
      PetscReal                     :: residualYieldStress                              ! residual yield stress
      PetscReal                     :: DuctileCouplingPower                             ! Coupling power between plasticity and damage
      PetscReal                     :: CoefficientCapModel0                             ! C0 in CapModel
      PetscReal                     :: CoefficientCapModel1                             ! C1 in CapModel
      PetscReal                     :: CoefficientCapModel2                             ! C2 in CapModel
      PetscReal                     :: CoefficientCapModelD                             ! CD in CapModel
      PetscReal                     :: CoefficientDruckerPrager                         ! k  in DruckerPrager
      PetscReal                     :: CoeffF                                           ! coefficient F Hill matrix
      PetscReal                     :: CoeffG                                           ! coefficient G Hill matrix
      PetscReal                     :: CoeffH                                           ! coefficient H Hill matrix
      PetscReal                     :: CoeffM                                           ! coefficient M Hill matrix
      PetscReal                     :: CoeffN                                           ! coefficient N Hill matrix
      PetscReal                     :: CoeffL                                           ! coefficient L Hill matrix
      PetscReal                     :: YieldTau0                                        ! Hill yield stress
      PetscReal                     :: residualYieldTau0                                ! Hill residual yield stress
      PetscReal                     :: phi1                                             ! Bunge Euler angle phi1
      PetscReal                     :: phi2                                             ! Bunge Euler angle phi2
      PetscReal                     :: Phi                                              ! Bunge Euler angle Phi
      PetscReal                     :: delta                                            ! residual Gurson and Green
      PetscReal                     :: cohesiveStiffness
      PetscReal                     :: drivingForceTensileStrength                      ! tensile strength in Drucker-Prager driving force
      PetscReal                     :: drivingForceCompressiveStrength                  ! compressive strength in Drucker-Prager driving force
      PetscReal                     :: drivingForceDelta                                ! delta parameter in Drucker-Prager driving force
      PetscReal                     :: drivingForceGamma                                ! gamma parameter in Drucker-Prager driving force
      PetscBool                     :: isLinearIsotropicHardening
      PetscBool                     :: isNoPlCoupling
      Type(MEF90RotationMatrix3D)   :: RotationMatrix                                   ! rotation matrix from the global frame to the material frame: X_local = R . X_global
      PetscBool                     :: isViscousPlasticity                              ! boolean telling if crystal plasticity is viscous or rate-independent
      PetscReal                     :: ViscosityGamma0                                  ! viscosity reference slip rate
      PetscReal                     :: ViscosityN                                       ! viscosity exponent
      PetscReal                     :: Viscositydt                                      ! time step size
      PetscReal                     :: m                                                ! equivalent stress exponent for rate-independent crystal plasticity
      Character(len=MEF90MXSTRLEN)  :: Name
   End Type MEF90MatProp2D_Type

   Type MEF90MatProp3D_Type
      PetscReal                     :: Density                                          ! rho
      PetscReal                     :: FractureToughness                                ! Gc
      Type(MatS3D)                  :: toughnessAnisotropyMatrix                        ! 
      PetscReal                     :: SpecificHeat                                     ! Cp
      Type(MatS3D)                  :: ThermalConductivity                              ! K
      Type(MatS3D)                  :: LinearThermalExpansion                           ! alpha
      Type(MEF90HookesLaw3D)        :: HookesLaw                                        ! A
      PetscReal                     :: internalLength                                   ! l
      PetscReal                     :: residualStiffness                                ! eta
      PetscReal                     :: yieldStress                                      ! yield stress
      PetscReal                     :: residualYieldStress                              ! residual yield stress
      PetscReal                     :: DuctileCouplingPower                             ! Coupling power between plasticity and damage
      PetscReal                     :: CoefficientCapModel0                             ! C0 in CapModel
      PetscReal                     :: CoefficientCapModel1                             ! C1 in CapModel
      PetscReal                     :: CoefficientCapModel2                             ! C2 in CapModel
      PetscReal                     :: CoefficientCapModelD                             ! CD in CapModel
      PetscReal                     :: CoefficientDruckerPrager                         ! k  in DruckerPrager
      PetscReal                     :: CoeffF                                           ! coefficient F Hill matrix
      PetscReal                     :: CoeffG                                           ! coefficient G Hill matrix
      PetscReal                     :: CoeffH                                           ! coefficient H Hill matrix
      PetscReal                     :: CoeffM                                           ! coefficient M Hill matrix
      PetscReal                     :: CoeffN                                           ! coefficient N Hill matrix
      PetscReal                     :: CoeffL                                           ! coefficient L Hill matrix
      PetscReal                     :: YieldTau0                                        ! Hill yield stress
      PetscReal                     :: residualYieldTau0                                ! Hill residual yield stress
      PetscReal                     :: phi1                                             ! Bunge Euler angle phi1
      PetscReal                     :: phi2                                             ! Bunge Euler angle phi2
      PetscReal                     :: Phi                                              ! Bunge Euler angle Phi
      PetscReal                     :: delta                                            ! residual Gurson and Green
      PetscReal                     :: cohesiveStiffness
      PetscReal                     :: drivingForceTensileStrength                      ! tensile strength in Drucker-Prager driving force
      PetscReal                     :: drivingForceCompressiveStrength                  ! compressive strength in Drucker-Prager driving force
      PetscReal                     :: drivingForceDelta                                ! delta parameter in Drucker-Prager driving force
      PetscReal                     :: drivingForceGamma                                ! gamma parameter in Drucker-Prager driving force
      PetscBool                     :: isLinearIsotropicHardening
      PetscBool                     :: isNoPlCoupling
      Type(MEF90RotationMatrix3D)   :: RotationMatrix                                   ! rotation matrix from the global frame to the material frame: X_local = R . X_global
      PetscBool                     :: isViscousPlasticity                              ! boolean telling if crystal plasticity is viscous or rate-independent
      PetscReal                     :: ViscosityGamma0                                  ! viscosity reference slip rate
      PetscReal                     :: ViscosityN                                       ! viscosity exponent
      PetscReal                     :: Viscositydt                                      ! time step size
      PetscReal                     :: m                                                ! equivalent stress exponent for rate-independent crystal plasticity
      Character(len=MEF90MXSTRLEN)  :: Name
   End Type MEF90MatProp3D_Type

   Enum,bind(c)
      enumerator :: MEF90HookesLawTypeFull = 0, &
                    MEF90HookesLawTypeIsotropic
   End Enum

   !!! The Mathium is a bogus isotropic material whose material properties are all 1
   !!! except for a Poisson Ratio of 0.3
   Type(MEF90MatProp2D_Type),Parameter     :: MEF90Mathium2D = MEF90MatProp2D_Type(    &
      1.0_Kr,                                                                          & ! Density
      1.0_Kr,                                                                          & ! FractureToughness
      MEF90MatS2DIdentity,                                                             & ! toughnessAnisotropyMatrix
      1.0_Kr,                                                                          & ! SpecificHeat
      MEF90MatS2DIdentity,                                                             & ! ThermalConductivity
      MEF90MatS2DIdentity,                                                             & ! LinearThermalExpansion
      MEF90HookesLaw2D(                                                                &
         Tens4OS2D(1.09890_Kr,0.32967_Kr,0.00000_Kr,                                   & ! HookesLaw XXXX,XXYY,XXXY
                              1.09890_Kr,0.00000_Kr,                                   & !                YYYY,YYXY
                                         0.38462_Kr),                                  & !                     XYXY
         Tens4OS3D(1.34615_Kr,0.57692_Kr,0.57692_Kr,0.00000_Kr,0.00000_Kr,0.00000_Kr,  & ! HookesLaw Local XXXX,XXYY,XXZZ,XXYZ,XXXZ,XXXY
                              1.34615_Kr,0.57692_Kr,0.00000_Kr,0.00000_Kr,0.00000_Kr,  & !                      YYYY,YYZZ,YYYZ,YYXZ,YYXY
                                         1.34615_Kr,0.00000_Kr,0.00000_Kr,0.00000_Kr,  & !                           ZZZZ ZZYZ,ZZXZ,ZZXY
                                                    0.38462_Kr,0.00000_Kr,0.00000_Kr,  & !                                YXYX,YZXZ,YZXY
                                                               0.38462_Kr,0.00000_Kr,  & !                                     XZXZ,XZXY
                                                                          0.38462_Kr), & !        
         Tens4OS3D(1.34615_Kr,0.57692_Kr,0.57692_Kr,0.00000_Kr,0.00000_Kr,0.00000_Kr,  & ! HookesLaw 3D    XXXX,XXYY,XXZZ,XXYZ,XXXZ,XXXY
                              1.34615_Kr,0.57692_Kr,0.00000_Kr,0.00000_Kr,0.00000_Kr,  & !                      YYYY,YYZZ,YYYZ,YYXZ,YYXY
                                         1.34615_Kr,0.00000_Kr,0.00000_Kr,0.00000_Kr,  & !                           ZZZZ ZZYZ,ZZXZ,ZZXY
                                                    0.38462_Kr,0.00000_Kr,0.00000_Kr,  & !                                YXYX,YZXZ,YZXY
                                                               0.38462_Kr,0.00000_Kr,  & !                                     XZXZ,XZXY
                                                                          0.38462_Kr), & ! 
         0.0_Kr,0.0_Kr,1.0_Kr,.3_Kr,0.0_Kr,                                            & ! lambda, mu, E, nu, kappa (lambda, mu, kappa will be recomputed)
         MEF90HookesLawTypeIsotropic,                                                  & ! type
         .FALSE.),                                                                     & ! isPlaneStress
      1.0_Kr,                                                                          & ! Internal Length
      1.0e-9,                                                                          & ! Residual Stiffness
      1.0_Kr,                                                                          & ! Yield Stress
      0.0_Kr,                                                                          & ! Residual Yield Stress
      2.0_Kr,                                                                          & ! Coupling power between damage and plasticity
      -0.3_Kr,                                                                         & ! C0 in CapModel
      0.4_Kr,                                                                          & ! C1 in CapModel
      1.0_Kr,                                                                          & ! C2 in CapModel
      1.0_Kr,                                                                          & ! CD in CapModel
      -0.5_Kr,                                                                         & ! k  in DruckerPrager
      0.3_Kr,                                                                          & ! coefficient F Hill matrix
      0.3_Kr,                                                                          & ! coefficient G Hill matrix
      0.3_Kr,                                                                          & ! coefficient H Hill matrix
      1.0_Kr,                                                                          & ! coefficient M Hill matrix
      1.0_Kr,                                                                          & ! coefficient N Hill matrix
      1.0_Kr,                                                                          & ! coefficient L Hill matrix
      1.0_Kr,                                                                          & ! Hill yield stress
      0.0_Kr,                                                                          & ! Hill residual yield stress
      0.0_Kr,                                                                          & ! Bunge Euler angle phi1
      0.0_Kr,                                                                          & ! Bunge Euler angle phi2
      0.0_Kr,                                                                          & ! Bunge Euler angle Phi
      0.0001_Kr,                                                                       & ! Residual Gurson and Green
      0.0_Kr,                                                                          & ! cohesive stiffness
      0.0_Kr,                                                                          & ! drivingForceTensileStrength
      0.0_Kr,                                                                          & ! drivingForceCompressiveStrength
      0.0_Kr,                                                                          & ! drivingForceDelta
      0.0_Kr,                                                                          & ! drivingForceGamma
      .FALSE.,                                                                         & ! isLinearIsotropicHardening
      .FALSE.,                                                                         & ! isNoPlCoupling
      MEF90RotationMatrix3D(MEF90Mat3DIdentity,0.0_Kr,0.0_Kr,0.0_Kr,                   & ! RotationMatrix, phi1, Phi, phi2,
      Vect3D(1.0_Kr,0.0_Kr,0.0_Kr),Vect3D(0.0_Kr,1.0_Kr,0.0_Kr),                       & ! V1, V2,
      Vect3D(0.0_Kr,0.0_Kr,1.0_Kr),.False.),                                           & ! V3, fromEuler
      .FALSE.,                                                                         & ! isViscousPlasticity
      1.0_Kr,                                                                          & ! ViscosityGamma0
      1.0_Kr,                                                                          & ! ViscosityN
      1.0_Kr,                                                                          & ! Viscositydt
      1.0_Kr,                                                                          & ! m
      "MEF90Mathium2D")

   Type(MEF90MatProp3D_Type),Parameter     :: MEF90Mathium3D = MEF90MatProp3D_Type(    &
      1.0_Kr,                                                                          & ! Density
      1.0_Kr,                                                                          & ! FractureToughness
      MEF90MatS3DIdentity,                                                             & ! toughnessAnisotropyMatrix
      1.0_Kr,                                                                          & ! SpecificHeat
      MEF90MatS3DIdentity,                                                             & ! ThermalConductivity
      MEF90MatS3DIdentity,                                                             & ! LinearThermalExpansion
      MEF90HookesLaw3D(                                                                &
         Tens4OS3D(1.34615_Kr,0.57692_Kr,0.57692_Kr,0.00000_Kr,0.00000_Kr,0.00000_Kr,  & ! HookesLaw XXXX,XXYY,XXZZ,XXYZ,XXXZ,XXXY
                              1.34615_Kr,0.57692_Kr,0.00000_Kr,0.00000_Kr,0.00000_Kr,  & !                YYYY,YYZZ,YYYZ,YYXZ,YYXY
                                         1.34615_Kr,0.00000_Kr,0.00000_Kr,0.00000_Kr,  & !                     ZZZZ ZZYZ,ZZXZ,ZZXY
                                                    0.38462_Kr,0.00000_Kr,0.00000_Kr,  & !                          YXYX,YZXZ,YZXY
                                                               0.38462_Kr,0.00000_Kr,  & !                               XZXZ,XZXY
                                                                          0.38462_Kr), & !                                    XYXY
         Tens4OS3D(1.34615_Kr,0.57692_Kr,0.57692_Kr,0.00000_Kr,0.00000_Kr,0.00000_Kr,  & ! HookesLaw Local XXXX,XXYY,XXZZ,XXYZ,XXXZ,XXXY
                              1.34615_Kr,0.57692_Kr,0.00000_Kr,0.00000_Kr,0.00000_Kr,  & !                      YYYY,YYZZ,YYYZ,YYXZ,YYXY
                                         1.34615_Kr,0.00000_Kr,0.00000_Kr,0.00000_Kr,  & !                           ZZZZ ZZYZ,ZZXZ,ZZXY
                                                    0.38462_Kr,0.00000_Kr,0.00000_Kr,  & !                                YXYX,YZXZ,YZXY
                                                               0.38462_Kr,0.00000_Kr,  & !                                     XZXZ,XZXY
                                                                          0.38462_Kr), & ! 
         0.0_Kr,0.0_Kr,1.0_Kr,.3_Kr,0.0_Kr,                                            & ! lambda, mu, E, nu, kappa (lambda, mu, kappa will be recomputed)
         MEF90HookesLawTypeIsotropic,                                                  & ! type
#if (PETSC_SIZEOF_INT == 4)
         PETSC_FALSE                                                                   & ! padding
#endif
      ),                                                                               & ! 
      1.0_Kr,                                                                          & ! Internal Length
      1.0e-9,                                                                          & ! Residual Stiffness
      1.0_Kr,                                                                          & ! Yield Stress
      0.0_Kr,                                                                          & ! Residual Yield Stress
      2.0_Kr,                                                                          & ! Coupling power between damage and plasticity
      -0.3_Kr,                                                                         & ! C0 in CapModel
      0.4_Kr,                                                                          & ! C1 in CapModel
      1.0_Kr,                                                                          & ! C2 in CapModel
      1.0_Kr,                                                                          & ! CD in CapModel
      -0.5_Kr,                                                                         & ! k  in DruckerPrager
      0.3_Kr,                                                                          & ! coefficient F Hill matrix
      0.3_Kr,                                                                          & ! coefficient G Hill matrix
      0.3_Kr,                                                                          & ! coefficient H Hill matrix
      1.0_Kr,                                                                          & ! coefficient M Hill matrix
      1.0_Kr,                                                                          & ! coefficient N Hill matrix
      1.0_Kr,                                                                          & ! coefficient L Hill matrix
      1.0_Kr,                                                                          & ! Hill yield stress
      0.0_Kr,                                                                          & ! Hill residual yield stress
      0.0_Kr,                                                                          & ! Bunge Euler angle phi1
      0.0_Kr,                                                                          & ! Bunge Euler angle phi2
      0.0_Kr,                                                                          & ! Bunge Euler angle Phi
      0.0001_Kr,                                                                       & ! Residual Gurson and Green
      0.0_Kr,                                                                          & ! cohesive stiffness
      0.0_Kr,                                                                          & ! drivingForceTensileStrength
      0.0_Kr,                                                                          & ! drivingForceCompressiveStrength
      1.0_Kr,                                                                          & ! drivingForceBeta
      2,                                                                               & ! drivingForcep
      .FALSE.,                                                                         & ! isLinearIsotropicHardening
      .FALSE.,                                                                         & ! isNoPlCoupling
      MEF90RotationMatrix3D(MEF90Mat3DIdentity,0.0_Kr,0.0_Kr,0.0_Kr,                   & ! RotationMatrix, phi1, Phi, phi2,
      Vect3D(1.0_Kr,0.0_Kr,0.0_Kr),Vect3D(0.0_Kr,1.0_Kr,0.0_Kr),                       & ! V1, V2,
      Vect3D(0.0_Kr,0.0_Kr,1.0_Kr),.False.),                                           & ! V3, fromEuler
      .FALSE.,                                                                         & ! isViscousPlasticity
      1.0_Kr,                                                                          & ! ViscosityGamma0
      1.0_Kr,                                                                          & ! ViscosityN
      1.0_Kr,                                                                          & ! Viscositydt
      1.0_Kr,                                                                          & ! m
      "MEF90Mathium3D")
End Module m_MEF90_Materials_Types

Module m_MEF90_Materials_Interface2D
#include "petsc/finclude/petsc.h"
   Use petsc
   Use m_MEF90_Materials_Types
   Implicit NONE
   Private
   Public :: PetscBagGetDataMEF90MatProp2D

   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF90_Materials_Types
         PetscBag                             :: bag
         type(MEF90MatProp2D_Type),pointer    :: data
         PetscErrorCode                       :: ierr
      End subroutine PetscBagGetData
   End interface
Contains
#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90MatProp2D"
   Subroutine PetscBagGetDataMEF90MatProp2D(bag,data,ierr)
      PetscBag                                :: bag
      type(MEF90MatProp2D_Type),pointer       :: data
      PetscErrorCode                          :: ierr
      PetscReal                               :: normV1,normV2,normV3
      Type(Tens4OS3D)                         :: HookesLaw3D

      PetscCall(PetscBagGetData(bag,data,ierr))
      Select case(data%HookesLaw%type)
         Case(MEF90HookesLawTypeFull)
            Continue
         Case(MEF90HookesLawTypeIsotropic)
            If (data%HookesLaw%isPlaneStress) Then
               data%HookesLaw%lambda      = data%HookesLaw%YoungsModulus * data%HookesLaw%PoissonRatio / (1.0_Kr - data%HookesLaw%PoissonRatio**2)
               data%HookesLaw%mu          = data%HookesLaw%YoungsModulus / (1.0_Kr + data%HookesLaw%PoissonRatio) * .5_Kr
               data%HookesLaw%BulkModulus = data%HookesLaw%lambda + 2.0_Kr * data%HookesLaw%mu / 3.0_Kr
            Else
               data%HookesLaw%lambda      = data%HookesLaw%YoungsModulus * data%HookesLaw%PoissonRatio / (1.0_Kr + data%HookesLaw%PoissonRatio) / (1.0_Kr - 2.0_Kr * data%HookesLaw%PoissonRatio)
               data%HookesLaw%mu          = data%HookesLaw%YoungsModulus / (1.0_Kr + data%HookesLaw%PoissonRatio) * .5_Kr
               data%HookesLaw%BulkModulus = data%HookesLaw%lambda + data%HookesLaw%mu
            End If
      End Select
      If (data%RotationMatrix%fromEuler) Then
         data%RotationMatrix%fullTensor%XX =  cos(data%RotationMatrix%phi1)*cos(data%RotationMatrix%phi2)-sin(data%RotationMatrix%phi1)*sin(data%RotationMatrix%phi2)*cos(data%RotationMatrix%Phi)
         data%RotationMatrix%fullTensor%XY =  sin(data%RotationMatrix%phi1)*cos(data%RotationMatrix%phi2)+cos(data%RotationMatrix%phi1)*sin(data%RotationMatrix%phi2)*cos(data%RotationMatrix%Phi)
         data%RotationMatrix%fullTensor%YX = -cos(data%RotationMatrix%phi1)*sin(data%RotationMatrix%phi2)-sin(data%RotationMatrix%phi1)*cos(data%RotationMatrix%phi2)*cos(data%RotationMatrix%Phi)
         data%RotationMatrix%fullTensor%YY = -sin(data%RotationMatrix%phi1)*sin(data%RotationMatrix%phi2)+cos(data%RotationMatrix%phi1)*cos(data%RotationMatrix%phi2)*cos(data%RotationMatrix%Phi)
         data%RotationMatrix%fullTensor%XZ =  sin(data%RotationMatrix%phi2)*sin(data%RotationMatrix%Phi)
         data%RotationMatrix%fullTensor%YZ =  cos(data%RotationMatrix%phi2)*sin(data%RotationMatrix%Phi)
         data%RotationMatrix%fullTensor%ZX =  sin(data%RotationMatrix%phi1)*sin(data%RotationMatrix%Phi)
         data%RotationMatrix%fullTensor%ZY = -cos(data%RotationMatrix%phi1)*sin(data%RotationMatrix%Phi)
         data%RotationMatrix%fullTensor%ZZ =  cos(data%RotationMatrix%Phi)
      Else
         normV1 = SQRT(data%RotationMatrix%V1%X*data%RotationMatrix%V1%X + data%RotationMatrix%V1%Y*data%RotationMatrix%V1%Y + data%RotationMatrix%V1%Z*data%RotationMatrix%V1%Z)
         normV2 = SQRT(data%RotationMatrix%V2%X*data%RotationMatrix%V2%X + data%RotationMatrix%V2%Y*data%RotationMatrix%V2%Y + data%RotationMatrix%V2%Z*data%RotationMatrix%V2%Z)
         normV3 = SQRT(data%RotationMatrix%V3%X*data%RotationMatrix%V3%X + data%RotationMatrix%V3%Y*data%RotationMatrix%V3%Y + data%RotationMatrix%V3%Z*data%RotationMatrix%V3%Z)
         data%RotationMatrix%fullTensor%XX = data%RotationMatrix%V1%X / normV1
         data%RotationMatrix%fullTensor%YX = data%RotationMatrix%V1%Y / normV1
         data%RotationMatrix%fullTensor%ZX = data%RotationMatrix%V1%Z / normV1
         data%RotationMatrix%fullTensor%XY = data%RotationMatrix%V2%X / normV2
         data%RotationMatrix%fullTensor%YY = data%RotationMatrix%V2%Y / normV2
         data%RotationMatrix%fullTensor%ZY = data%RotationMatrix%V2%Z / normV2
         data%RotationMatrix%fullTensor%XZ = data%RotationMatrix%V3%X / normV3
         data%RotationMatrix%fullTensor%YZ = data%RotationMatrix%V3%Y / normV3
         data%RotationMatrix%fullTensor%ZZ = data%RotationMatrix%V3%Z / normV3
      End If
      HookesLaw3D = Tens4OSTransform(data%HookesLaw%fullTensorLocal,transpose(data%RotationMatrix%fullTensor))
      data%HookesLaw%fullTensor3D = HookesLaw3D
      data%HookesLaw%fullTensor%XXXX = HookesLaw3D%XXXX
      data%HookesLaw%fullTensor%XXYY = HookesLaw3D%XXYY
      data%HookesLaw%fullTensor%XXXY = HookesLaw3D%XXXY
      data%HookesLaw%fullTensor%YYYY = HookesLaw3D%YYYY
      data%HookesLaw%fullTensor%YYXY = HookesLaw3D%YYXY
      data%HookesLaw%fullTensor%XYXY = HookesLaw3D%XYXY
   End Subroutine PetscBagGetDataMEF90MatProp2D
End Module m_MEF90_Materials_Interface2D

Module m_MEF90_Materials_Interface3D
#include "petsc/finclude/petsc.h"
   Use petsc
   Use m_MEF90_Materials_Types
   Implicit NONE
   Private
   Public :: PetscBagGetDataMEF90MatProp3D

   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF90_Materials_Types
         PetscBag                             :: bag
         type(MEF90MatProp3D_Type),pointer    :: data
         PetscErrorCode                       :: ierr
      End subroutine PetscBagGetData
   End interface
Contains
#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataMEF90MatProp3D"
   Subroutine PetscBagGetDataMEF90MatProp3D(bag,data,ierr)
      PetscBag                                :: bag
      type(MEF90MatProp3D_Type),pointer       :: data
      PetscErrorCode                          :: ierr
      PetscReal                               :: normV1,normV2,normV3

      PetscCall(PetscBagGetData(bag,data,ierr))
      Select case(data%HookesLaw%type)
         Case(MEF90HookesLawTypeFull)
            Continue
         Case(MEF90HookesLawTypeIsotropic)
            data%HookesLaw%lambda       = data%HookesLaw%YoungsModulus * data%HookesLaw%PoissonRatio / (1.0_Kr + data%HookesLaw%PoissonRatio) / (1.0_Kr - 2.0_Kr * data%HookesLaw%PoissonRatio)
            data%HookesLaw%mu           = data%HookesLaw%YoungsModulus / (1.0_Kr + data%HookesLaw%PoissonRatio) * .5_Kr
            data%HookesLaw%BulkModulus  = data%HookesLaw%lambda + 2.0_Kr * data%HookesLaw%mu / 3.0_Kr
      End Select
      If (data%RotationMatrix%fromEuler) Then
         data%RotationMatrix%fullTensor%XX =  cos(data%RotationMatrix%phi1)*cos(data%RotationMatrix%phi2)-sin(data%RotationMatrix%phi1)*sin(data%RotationMatrix%phi2)*cos(data%RotationMatrix%Phi)
         data%RotationMatrix%fullTensor%XY =  sin(data%RotationMatrix%phi1)*cos(data%RotationMatrix%phi2)+cos(data%RotationMatrix%phi1)*sin(data%RotationMatrix%phi2)*cos(data%RotationMatrix%Phi)
         data%RotationMatrix%fullTensor%YX = -cos(data%RotationMatrix%phi1)*sin(data%RotationMatrix%phi2)-sin(data%RotationMatrix%phi1)*cos(data%RotationMatrix%phi2)*cos(data%RotationMatrix%Phi)
         data%RotationMatrix%fullTensor%YY = -sin(data%RotationMatrix%phi1)*sin(data%RotationMatrix%phi2)+cos(data%RotationMatrix%phi1)*cos(data%RotationMatrix%phi2)*cos(data%RotationMatrix%Phi)
         data%RotationMatrix%fullTensor%XZ =  sin(data%RotationMatrix%phi2)*sin(data%RotationMatrix%Phi)
         data%RotationMatrix%fullTensor%YZ =  cos(data%RotationMatrix%phi2)*sin(data%RotationMatrix%Phi)
         data%RotationMatrix%fullTensor%ZX =  sin(data%RotationMatrix%phi1)*sin(data%RotationMatrix%Phi)
         data%RotationMatrix%fullTensor%ZY = -cos(data%RotationMatrix%phi1)*sin(data%RotationMatrix%Phi)
         data%RotationMatrix%fullTensor%ZZ =  cos(data%RotationMatrix%Phi)
      Else
         normV1 = SQRT(data%RotationMatrix%V1%X*data%RotationMatrix%V1%X + data%RotationMatrix%V1%Y*data%RotationMatrix%V1%Y + data%RotationMatrix%V1%Z*data%RotationMatrix%V1%Z)
         normV2 = SQRT(data%RotationMatrix%V2%X*data%RotationMatrix%V2%X + data%RotationMatrix%V2%Y*data%RotationMatrix%V2%Y + data%RotationMatrix%V2%Z*data%RotationMatrix%V2%Z)
         normV3 = SQRT(data%RotationMatrix%V3%X*data%RotationMatrix%V3%X + data%RotationMatrix%V3%Y*data%RotationMatrix%V3%Y + data%RotationMatrix%V3%Z*data%RotationMatrix%V3%Z)
         data%RotationMatrix%fullTensor%XX = data%RotationMatrix%V1%X / normV1
         data%RotationMatrix%fullTensor%YX = data%RotationMatrix%V1%Y / normV1
         data%RotationMatrix%fullTensor%ZX = data%RotationMatrix%V1%Z / normV1
         data%RotationMatrix%fullTensor%XY = data%RotationMatrix%V2%X / normV2
         data%RotationMatrix%fullTensor%YY = data%RotationMatrix%V2%Y / normV2
         data%RotationMatrix%fullTensor%ZY = data%RotationMatrix%V2%Z / normV2
         data%RotationMatrix%fullTensor%XZ = data%RotationMatrix%V3%X / normV3
         data%RotationMatrix%fullTensor%YZ = data%RotationMatrix%V3%Y / normV3
         data%RotationMatrix%fullTensor%ZZ = data%RotationMatrix%V3%Z / normV3
      End If
      data%HookesLaw%fullTensor = Tens4OSTransform(data%HookesLaw%fullTensorLocal,transpose(data%RotationMatrix%fullTensor))
   End Subroutine PetscBagGetDataMEF90MatProp3D
End Module m_MEF90_Materials_Interface3D

Module m_MEF90_Materials
#include "petsc/finclude/petsc.h"
   Use petsc
   Use m_MEF90_Materials_Types
   Use m_MEF90_Materials_Interface2D
   Use m_MEF90_Materials_Interface3D
   Implicit NONE

   Interface PetscBagGetDataMEF90MatProp
      Module Procedure PetscBagGetDataMEF90MatProp2D,PetscBagGetDataMEF90MatProp3D
   End Interface

   Interface PetscBagRegisterMEF90MatProp
      Module Procedure PetscBagRegisterMEF90MatProp2D,PetscBagRegisterMEF90MatProp3D
   End Interface

   Interface MEF90MatPropBagSetFromOptions
      Module Procedure MEF90MatPropBagSetFromOptions2D,MEF90MatPropBagSetFromOptions3D
   End Interface

   Interface  Operator (+)
      Module Procedure MEF90HookesLaw2DSum,MEF90HookesLaw3DSum
   End Interface

   Interface  Operator (-)
      Module Procedure MEF90HookesLaw2DDiff,MEF90HookesLaw3DDiff
   End Interface

   Interface  Operator (*)
      Module Procedure MEF90HookesLaw2DXMatS2D,MEF90HookesLaw3DXMatS3D,MEF90HookesLaw2DXMat2D,MEF90HookesLaw3DXMat3D, ScalarXMEF90HookesLaw2D, ScalarXMEF90HookesLaw3D
   End Interface

   PetscSizeT,protected   :: sizeofMEF90MatProp2D
   PetscSizeT,protected   :: sizeofMEF90MatProp3D

   PetscSizeT,protected   :: sizeofMEF90HookesLaw2D
   PetscSizeT,protected   :: sizeofMEF90HookesLaw3D

   Character(len =  MEF90MXSTRLEN) ,Dimension(5),protected   :: MEF90HookesLawTypeList

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90MaterialsInitialize_Private"
!!!
!!!
!!!  MEF90MaterialsInitialize_Private:
!!!
!!!  (c) 2013-2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90MaterialsInitialize_Private(ierr)
      PetscErrorCode,intent(OUT)          :: ierr

      Type(MEF90MatProp2D_Type),Target    :: matProp2D
      Type(MEF90MatProp3D_Type),Target    :: matProp3D
      Type(MEF90HookesLaw2D),Target       :: HookesLaw2D
      Type(MEF90HookesLaw3D),Target       :: HookesLaw3D
      character(len=1),pointer            :: dummychar(:)
      PetscSizeT                          :: sizeofchar

      PetscCall(PetscDataTypeGetSize(PETSC_CHAR,sizeofchar,ierr))
      sizeofMEF90MatProp2D = size(transfer(matProp2D,dummychar))*sizeofchar
      sizeofMEF90MatProp3D = size(transfer(matProp3D,dummychar))*sizeofchar
      sizeofMEF90HookesLaw2D = size(transfer(HookesLaw2D,dummychar))*sizeofchar
      sizeofMEF90HookesLaw3D = size(transfer(HookesLaw3D,dummychar))*sizeofchar

      MEF90HookesLawTypeList(1) = 'Full'
      MEF90HookesLawTypeList(2) = 'Isotropic'
      MEF90HookesLawTypeList(3) = 'MEF90HookesLawTypeList'
      MEF90HookesLawTypeList(4) = '_MEF90HookesLawTypeList'
      MEF90HookesLawTypeList(5) = ''
   End Subroutine MEF90MaterialsInitialize_Private

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90MatProp2D"
!!!
!!!
!!!  PetscBagRegisterMEF90MatProp2D:
!!!
!!!  (c) 2013-2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine PetscBagRegisterMEF90MatProp2D(bag,name,prefix,default,ierr)
      PetscBag                               :: bag
      Character(len=*),intent(IN)            :: prefix,name
      type(MEF90MatProp2D_Type),intent(IN)   :: default
      PetscErrorCode,intent(OUT)             :: ierr

      Type(MEF90MatProp2D_Type),pointer      :: matprop

      PetscCall(PetscBagGetDataMEF90MatProp2D(bag,matprop,ierr))
      PetscCall(PetscBagSetName(bag,trim(name),"MatProp2D object: material properties",ierr))
      PetscCall(PetscBagSetOptionsPrefix(bag,trim(prefix), ierr))

      PetscCall(PetscBagRegisterString(bag,matprop%name,trim(default%name),'Name','',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%density,default%density,'Density','[kg.m^(-2)] (rho) Density',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%FractureToughness,default%FractureToughness,'FractureToughness','[N.m^(-1)] (G_c) Fracture toughness',ierr))
      matprop%toughnessAnisotropyMatrix = default%toughnessAnisotropyMatrix
      PetscCall(PetscBagRegisterRealArray(bag,matprop%toughnessAnisotropyMatrix,3_Ki,'toughnessAnisotropyMatrix','[] toughness Anisotropy Matrix',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%SpecificHeat,default%SpecificHeat,'SpecificHeat','[J.kg^(-1).K^(-1)] (Cp) Specific heat',ierr))
      matprop%ThermalConductivity = default%ThermalConductivity
      PetscCall(PetscBagRegisterRealArray(bag,matprop%ThermalConductivity,3_Ki,'ThermalConductivity','[J.m^(-1).s^(-1).K^(-1)] (K) Thermal conductivity',ierr))
      matprop%LinearThermalExpansion = default%LinearThermalExpansion
      PetscCall(PetscBagRegisterRealArray(bag,matprop%LinearThermalExpansion,3_Ki,'LinearThermalExpansion','[K^(-1)] (alpha) Linear thermal expansion matrix',ierr))

      PetscCall(PetscBagRegisterEnum(bag,matprop%HookesLaw%type,MEF90HookesLawTypeList,default%HookesLaw%type,'hookeslaw_type','Type of Hooke''s law',ierr))
      Select case(matprop%HookesLaw%type)
         Case (MEF90HookesLawTypeFull)
            matprop%HookesLaw%fullTensorLocal = default%HookesLaw%fullTensorLocal
            PetscCall(PetscBagRegisterRealArray(bag,matprop%HookesLaw%fullTensorLocal,21_Ki,'HookesLaw_tensor','[N.m^(-2)] (A) Hooke''s law in the local frame',ierr))
         Case(MEF90HookesLawTypeIsotropic)
            PetscCall(PetscBagRegisterReal(bag,matprop%HookesLaw%YoungsModulus,default%HookesLaw%YoungsModulus,'hookeslaw_YoungsModulus','[N.m^(-2)] (E) Young''s Modulus',ierr))
            PetscCall(PetscBagRegisterReal(bag,matprop%HookesLaw%PoissonRatio,default%HookesLaw%PoissonRatio,'hookeslaw_PoissonRatio','[] (nu) Poisson Modulus',ierr))
            PetscCall(PetscBagRegisterBool(bag,matprop%HookesLaw%isPlaneStress,default%HookesLaw%isPlaneStress,'hookeslaw_planeStress','Use plane stress elasticity',ierr))
            matprop%HookesLaw%fulltensor = -1.D+30
            matprop%HookesLaw%fulltensorLocal = -1.D+30
      End Select

      PetscCall(PetscBagRegisterReal(bag,matprop%internalLength,default%internalLength,'internalLength','[m] (l) Internal Length',ierr))

      PetscCall(PetscBagRegisterReal(bag,matprop%yieldStress,default%yieldStress,'yieldStress','[N.m^(-2)] (sigma_y) stress threshold for plasticity',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%residualYieldStress,default%residualYieldStress,'residualyieldStress','[unit-less] (eta) residual yield stress',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%DuctileCouplingPower,default%DuctileCouplingPower,'DuctileCouplingPower','[] power of the coupling between the damage and the plasticity',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%CoefficientCapModel0,default%CoefficientCapModel0,'CoefficientCapModel0','C0 in the Yield function: CD || dev(stress) || - C2 tr(stress)^2 - C1 tr(stress) - C0 <= 0',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%CoefficientCapModel1,default%CoefficientCapModel1,'CoefficientCapModel1','C1 in the Yield function: CD || dev(stress) || - C2 tr(stress)^2 - C1 tr(stress) - C0 <= 0',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%CoefficientCapModel2,default%CoefficientCapModel2,'CoefficientCapModel2','C2 in the Yield function: CD || dev(stress) || - C2 tr(stress)^2 - C1 tr(stress) - C0 <= 0',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%CoefficientCapModelD,default%CoefficientCapModelD,'CoefficientCapModelD','CD in the Yield function: CD || dev(stress) || - C2 tr(stress)^2 - C1 tr(stress) - C0 <= 0',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%CoefficientDruckerPrager,default%CoefficientDruckerPrager,'CoefficientDruckerPrager','k in the Yield function: || dev(stress) || - k tr(stress) - yieldStress <= 0',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%CoeffF,default%CoeffF,'CoeffF','[unit-less] (F) coefficient F in the Hill yield criterion',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%CoeffG,default%CoeffG,'CoeffG','[unit-less] (G) coefficient G in the Hill yield criterion',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%CoeffH,default%CoeffH,'CoeffH','[unit-less] (H) coefficient H in the Hill yield criterion',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%CoeffM,default%CoeffM,'CoeffM','[unit-less] (M) coefficient M in the Hill yield criterion',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%CoeffN,default%CoeffN,'CoeffN','[unit-less] (N) coefficient N in the Hill yield criterion',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%CoeffL,default%CoeffL,'CoeffL','[unit-less] (L) coefficient L in the Hill yield criterion',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%YieldTau0,default%YieldTau0,'YieldTau0','[N.m^(-2)] (tau_0) stress threshold in the Hill yield criterion',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%residualYieldTau0,default%residualYieldTau0,'residualYieldTau0','[unit-less] residual stress threshold in the Hill yield criterion',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%phi1,default%phi1,'phi1','[radians] Bunge-Euler angle in the Hill yield criterion',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%phi2,default%phi2,'phi2','[radians] Bunge-Euler angle in the Hill yield criterion',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%Phi,default%Phi,'Phi','[radians] Bunge-Euler angle in the Hill yield criterion',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%delta,default%delta,'delta','[unit-less] residual in the definition of the porosity, Gurson and Green criteria',ierr))

      PetscCall(PetscBagRegisterReal(bag,matprop%cohesiveStiffness,default%cohesiveStiffness,'cohesiveStiffness','[N.m^(-4)] (k) cohesive stiffness in Winkler-type models',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%residualStiffness,default%residualStiffness,'residualStiffness','[unit-less] (eta) residual stiffness',ierr))

      PetscCall(PetscBagRegisterReal(bag,matprop%drivingForceTensileStrength,default%drivingForceTensileStrength,'drivingForce_tensileStrength','[N.m^(-2)] (\sigma_{ts}) tensile strength in Drucker-Prager driving Force',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%drivingForceCompressiveStrength,default%drivingForceCompressiveStrength,'drivingForce_CompressiveStrength','[N.m^(-2)] (\sigma_{cs}) compressive strength in Drucker-Prager driving Force',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%drivingForceDelta,default%drivingForceDelta,'drivingForce_Delta','[unit-less] (\delta) delta parameter in Drucker-Prager driving Force',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%drivingForceGamma,default%drivingForceGamma,'drivingForce_Gamma','[m^3.N^(-2)] (\gamma) gamma parameter in Drucker-Prager driving Force',ierr))

      PetscCall(PetscBagRegisterBool(bag,matprop%isLinearIsotropicHardening,default%isLinearIsotropicHardening,'isLinearIsotropicHardening','[bool] Plasticity with Linear Isotropic Hardening',ierr))
      PetscCall(PetscBagRegisterBool(bag,matprop%isNoPlCoupling,default%isNoPlCoupling,'isNoPlCoupling','[bool] Coupling between damage and plastic dissipation',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%RotationMatrix%phi1,default%RotationMatrix%phi1,'RotationMatrix_phi1','[radians] (phi1) First Bunge-Euler angle',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%RotationMatrix%Phi,default%RotationMatrix%Phi,'RotationMatrix_Phi','[radians] (Phi) Second Bunge-Euler angle',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%RotationMatrix%phi2,default%RotationMatrix%phi2,'RotationMatrix_phi2','[radians] (phi2) Third Bunge-Euler angle',ierr))
      matprop%RotationMatrix%V1 = default%RotationMatrix%V1
      PetscCall(PetscBagRegisterRealArray(bag,matprop%RotationMatrix%V1,3_Ki,'RotationMatrix_V1','[] (V1) First column of the rotation matrix',ierr))
      matprop%RotationMatrix%V2 = default%RotationMatrix%V2
      PetscCall(PetscBagRegisterRealArray(bag,matprop%RotationMatrix%V2,3_Ki,'RotationMatrix_V2','[] (V2) Second column of the rotation matrix',ierr))
      matprop%RotationMatrix%V3 = default%RotationMatrix%V3
      PetscCall(PetscBagRegisterRealArray(bag,matprop%RotationMatrix%V3,3_Ki,'RotationMatrix_V3','[] (V3) Third column of the rotation matrix',ierr))    
      PetscCall(PetscBagRegisterBool(bag,matprop%RotationMatrix%fromEuler,default%RotationMatrix%fromEuler,'RotationMatrix_fromEuler','Define rotation matrix from Bunge-Euler angles',ierr))
      PetscCall(PetscBagRegisterBool(bag,matprop%isViscousPlasticity,default%isViscousPlasticity,'isViscousPlasticity','[bool] Viscous plastic potential',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%ViscosityGamma0,default%ViscosityGamma0,'ViscosityGamma0','[s^(-1)] Reference plastic deformation rate',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%ViscosityN,default%ViscosityN,'ViscosityN','[unit-less] Viscosity exponent',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%Viscositydt,default%Viscositydt,'Viscositydt','[s] Viscosity time step size',ierr))

      PetscCall(PetscBagRegisterReal(bag,matprop%m,default%m,'m','[unit-less] Equivalent stress exponent for rate-independent crystal plasticity',ierr))
   End Subroutine PetscBagRegisterMEF90MatProp2D

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90MatProp3D"
!!!
!!!
!!!  PetscBagRegisterMEF90MatProp3D:
!!!
!!!  (c) 2013-2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine PetscBagRegisterMEF90MatProp3D(bag,name,prefix,default,ierr)
      PetscBag                               :: bag
      Character(len=*),intent(IN)            :: prefix,name
      type(MEF90MatProp3D_Type),intent(IN)   :: default
      PetscErrorCode,intent(OUT)             :: ierr

      Type(MEF90MatProp3D_Type),pointer      :: matprop
      PetscCall(PetscBagGetDataMEF90MatProp3D(bag,matprop,ierr))
      PetscCall(PetscBagSetName(bag,trim(name),"MatProp3D object: material properties",ierr))
      PetscCall(PetscBagSetOptionsPrefix(bag,trim(prefix),ierr))

      PetscCall(PetscBagRegisterString(bag,matprop%name,trim(default%name),'Name','',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%density,default%density,'Density','[kg.m^(-3)] (rho) Density',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%FractureToughness,default%FractureToughness,'FractureToughness','[N.m^(-1)] (G_c) Fracture toughness',ierr))
      matprop%toughnessAnisotropyMatrix = default%toughnessAnisotropyMatrix
      PetscCall(PetscBagRegisterRealArray(bag,matprop%toughnessAnisotropyMatrix,6_Ki,'toughnessAnisotropyMatrix','[] toughness Anisotropy Matrix',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%SpecificHeat,default%SpecificHeat,'SpecificHeat','[J.kg^(-1).K^(-1)] (Cp) Specific heat',ierr))
      matprop%ThermalConductivity = default%ThermalConductivity
      PetscCall(PetscBagRegisterRealArray(bag,matprop%ThermalConductivity,6_Ki,'ThermalConductivity','[J.m^(-1).s^(-1).K^(-1)] (K) Thermal conductivity',ierr))
      matprop%LinearThermalExpansion = default%LinearThermalExpansion
      PetscCall(PetscBagRegisterRealArray(bag,matprop%LinearThermalExpansion,6_Ki,'LinearThermalExpansion','[K^(-1)] (alpha) Linear thermal expansion matrix',ierr))

      PetscCall(PetscBagRegisterEnum(bag,matprop%HookesLaw%type,MEF90HookesLawTypeList,default%HookesLaw%type,'hookeslaw_type','Type of Hooke''s law',ierr))
      Select case(matprop%HookesLaw%type)
         Case (MEF90HookesLawTypeFull)
            matprop%HookesLaw%fullTensorLocal = default%HookesLaw%fullTensorLocal
            PetscCall(PetscBagRegisterRealArray(bag,matprop%HookesLaw%fullTensorLocal,21_Ki,'HookesLaw_tensor','[N.m^(-2)] (A) Hooke''s law in the local frame',ierr))
         Case(MEF90HookesLawTypeIsotropic)
            PetscCall(PetscBagRegisterReal(bag,matprop%HookesLaw%YoungsModulus,default%HookesLaw%YoungsModulus,'hookeslaw_YoungsModulus','[N.m^(-2)] (E) Young''s Modulus',ierr))
            PetscCall(PetscBagRegisterReal(bag,matprop%HookesLaw%PoissonRatio,default%HookesLaw%PoissonRatio,'hookeslaw_PoissonRatio','[] (nu) Poisson Modulus',ierr))
            matprop%HookesLaw%fulltensor = -1.D+30
      End Select
      PetscCall(PetscBagRegisterReal(bag,matprop%internalLength,default%internalLength,'internalLength','[m] (l) Internal Length',ierr))

      PetscCall(PetscBagRegisterReal(bag,matprop%yieldStress,default%yieldStress,'yieldStress','[N.m^(-2)] (sigma_y) stress threshold for plasticity',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%residualYieldStress,default%residualYieldStress,'residualyieldStress','[unit-less] percentage of the yield stress',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%DuctileCouplingPower,default%DuctileCouplingPower,'DuctileCouplingPower','[] power of the coupling between the damage and the plasticity',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%CoefficientCapModel0,default%CoefficientCapModel0,'CoefficientCapModel0','C0 in the Yield function: CD || dev(stress) || + C2 tr(stress)^2 + C1 tr(stress) - C0 <= 0',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%CoefficientCapModel1,default%CoefficientCapModel1,'CoefficientCapModel1','C1 in the Yield function: CD || dev(stress) || + C2 tr(stress)^2 + C1 tr(stress) - C0 <= 0',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%CoefficientCapModel2,default%CoefficientCapModel2,'CoefficientCapModel2','C2 in the Yield function: CD || dev(stress) || + C2 tr(stress)^2 + C1 tr(stress) - C0 <= 0',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%CoefficientCapModelD,default%CoefficientCapModelD,'CoefficientCapModelD','CD in the Yield function: CD || dev(stress) || + C2 tr(stress)^2 + C1 tr(stress) - C0 <= 0',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%CoefficientDruckerPrager,default%CoefficientDruckerPrager,'CoefficientDruckerPrager','k in the Yield function: || dev(stress) || - k tr(stress) - yieldStress <= 0',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%CoeffF,default%CoeffF,'CoeffF','[unit-less] (F) coefficient F in the Hill yield criterion',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%CoeffG,default%CoeffG,'CoeffG','[unit-less] (G) coefficient G in the Hill yield criterion',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%CoeffH,default%CoeffH,'CoeffH','[unit-less] (H) coefficient H in the Hill yield criterion',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%CoeffM,default%CoeffM,'CoeffM','[unit-less] (M) coefficient M in the Hill yield criterion',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%CoeffN,default%CoeffN,'CoeffN','[unit-less] (N) coefficient N in the Hill yield criterion',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%CoeffL,default%CoeffL,'CoeffL','[unit-less] (L) coefficient L in the Hill yield criterion',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%YieldTau0,default%YieldTau0,'YieldTau0','[N.m^(-2)] (tau_0) stress threshold in the Hill yield criterion',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%residualYieldTau0,default%residualYieldTau0,'residualYieldTau0','[unit-less] residual stress threshold in the Hill yield criterion',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%phi1,default%phi1,'phi1','[radians] Bunge-Euler angle in the Hill yield criterion',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%phi2,default%phi2,'phi2','[radians] Bunge-Euler angle in the Hill yield criterion',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%Phi,default%Phi,'Phi','[radians] Bunge-Euler angle in the Hill yield criterion',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%delta,default%delta,'delta','[unit-less] residual in the definition of the porosity, Gurson and Green criteria',ierr))

      PetscCall(PetscBagRegisterReal(bag,matprop%cohesiveStiffness,default%cohesiveStiffness,'cohesiveStiffness','[N.m^(-4)] (k) cohesive stiffness in Winkler-type models',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%residualStiffness,default%residualStiffness,'residualStiffness','[unit-less] (eta) residual stiffness',ierr))

      PetscCall(PetscBagRegisterReal(bag,matprop%drivingForceTensileStrength,default%drivingForceTensileStrength,'drivingForce_tensileStrength','[N.m^(-2)] (\sigma_{ts}) tensile strength in Drucker-Prager driving Force',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%drivingForceCompressiveStrength,default%drivingForceCompressiveStrength,'drivingForce_CompressiveStrength','[N.m^(-2)] (\sigma_{cs}) compressive strength in Drucker-Prager driving Force',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%drivingForceDelta,default%drivingForceDelta,'drivingForce_Delta','[unit-less] (\delta) delta parameter in Drucker-Prager driving Force',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%drivingForceGamma,default%drivingForceGamma,'drivingForce_Gamma','[m^3.N^(-2)] (\gamma) gamma parameter in Drucker-Prager driving Force',ierr))

      PetscCall(PetscBagRegisterBool(bag,matprop%isLinearIsotropicHardening,default%isLinearIsotropicHardening,'isLinearIsotropicHardening','[bool] Plasticity with Linear Isotropic Hardening',ierr))
      PetscCall(PetscBagRegisterBool(bag,matprop%isNoPlCoupling,default%isNoPlCoupling,'isNoPlCoupling','[bool] Coupling between damage and plastic dissipation',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%RotationMatrix%phi1,default%RotationMatrix%phi1,'RotationMatrix_phi1','[radians] (phi1) First Bunge-Euler angle',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%RotationMatrix%Phi,default%RotationMatrix%Phi,'RotationMatrix_Phi','[radians] (Phi) Second Bunge-Euler angle',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%RotationMatrix%phi2,default%RotationMatrix%phi2,'RotationMatrix_phi2','[radians] (phi2) Third Bunge-Euler angle',ierr))
      matprop%RotationMatrix%V1 = default%RotationMatrix%V1
      PetscCall(PetscBagRegisterRealArray(bag,matprop%RotationMatrix%V1,3_Ki,'RotationMatrix_V1','[] (V1) First column of the rotation matrix',ierr))
      matprop%RotationMatrix%V2 = default%RotationMatrix%V2
      PetscCall(PetscBagRegisterRealArray(bag,matprop%RotationMatrix%V2,3_Ki,'RotationMatrix_V2','[] (V2) Second column of the rotation matrix',ierr))
      matprop%RotationMatrix%V3 = default%RotationMatrix%V3
      PetscCall(PetscBagRegisterRealArray(bag,matprop%RotationMatrix%V3,3_Ki,'RotationMatrix_V3','[] (V3) Third column of the rotation matrix',ierr))    
      PetscCall(PetscBagRegisterBool(bag,matprop%RotationMatrix%fromEuler,default%RotationMatrix%fromEuler,'RotationMatrix_fromEuler','Define rotation matrix from Bunge-Euler angles',ierr))
      PetscCall(PetscBagRegisterBool(bag,matprop%isViscousPlasticity,default%isViscousPlasticity,'isViscousPlasticity','[bool] Viscous plastic potential',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%ViscosityGamma0,default%ViscosityGamma0,'ViscosityGamma0','[s^(-1)] Reference plastic deformation rate',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%ViscosityN,default%ViscosityN,'ViscosityN','[unit-less] Viscosity exponent',ierr))
      PetscCall(PetscBagRegisterReal(bag,matprop%Viscositydt,default%Viscositydt,'Viscositydt','[s] Viscosity time step size',ierr))

      PetscCall(PetscBagRegisterReal(bag,matprop%m,default%m,'m','[unit-less] Equivalent stress exponent for rate-independent crystal plasticity',ierr))
   End Subroutine PetscBagRegisterMEF90MatProp3D


#undef __FUNCT__
#define __FUNCT__ "MEF90MatPropBagSetFromOptions2D"
!!!
!!!
!!!  MEF90MatPropBagSetFromOptionsierr:
!!!  MEF90MatPropBagSetFromOptions2D
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90MatPropBagSetFromOptions2D(MEF90MatPropBag,dm,defaultMaterial,MEF90Ctx,ierr)
      PetscBag,Dimension(:),Pointer                   :: MEF90MatPropBag
      Type(tDM),Intent(IN)                            :: dm
      Type(MEF90MatProp2D_Type),intent(IN)            :: defaultMaterial
      PetscErrorCode,Intent(OUT)                      :: ierr
      Type(MEF90Ctx_Type),Intent(IN)                  :: MEF90Ctx

      Type(tIS)                                       :: setIS
      PetscInt,Dimension(:),Pointer                   :: setID
      PetscInt                                        :: numSet,set
      Character(len=MEF90MXSTRLEN)                    :: setName,setprefix,IOBuffer
      Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90GlobalOptions

      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr))
      PetscCall(DMGetLabelIdIS(dm,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr))
      PetscCall(ISGetLocalSize(setIS,numSet,ierr))
      PetscCall(ISGetIndicesF90(setIS,setID,ierr))
      Allocate(MEF90MatPropBag(numSet))
      Do set = 1,numSet
         Write(setName,100) setID(set)
         Write(setprefix,101) setID(set)
         If (MEF90GlobalOptions%verbose > 0) Then
            Write(IOBuffer,102) setID(set),trim(setprefix)
            PetscCall(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
         End If
         PetscCall(PetscBagCreate(PETSC_COMM_WORLD,sizeofMEF90MatProp2D,MEF90MatPropBag(set),ierr))
         PetscCall(PetscBagRegisterMEF90MatProp(MEF90MatPropBag(set),setName,setprefix,defaultMaterial,ierr))
         If (MEF90GlobalOptions%verbose > 0) Then
            PetscCall(PetscBagView(MEF90MatPropBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr))
            PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n",ierr))
         End If
      End Do
      PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      PetscCall(ISDestroy(setIS,ierr))
100 Format('Cell set ',I4)
101 Format('cs',I4.4,'_')
102 Format('Registering materials properties for cell set ', I4,': ',A,'\n')
   End Subroutine MEF90MatPropBagSetFromOptions2D

#undef __FUNCT__
#define __FUNCT__ "MEF90MatPropBagSetFromOptions3D"
!!!
!!!
!!!  MEF90MatPropBagSetFromOptions3D:
!!!
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90MatPropBagSetFromOptions3D(MEF90MatPropBag,dm,defaultMaterial,MEF90Ctx,ierr)
      PetscBag,Dimension(:),Pointer                   :: MEF90MatPropBag
      Type(tDM),Intent(IN)                            :: dm
      Type(MEF90MatProp3D_Type),intent(IN)            :: defaultMaterial
      PetscErrorCode,Intent(OUT)                      :: ierr
      Type(MEF90Ctx_Type),Intent(IN)                  :: MEF90Ctx

      Type(tIS)                                       :: setIS
      PetscInt,Dimension(:),Pointer                   :: setID
      PetscInt                                        :: numSet,set
      Character(len=MEF90MXSTRLEN)                    :: setName,setprefix,IOBuffer
      Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90GlobalOptions

      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr))
      PetscCall(DMGetLabelIdIS(dm,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr))
      PetscCall(ISGetLocalSize(setIS,numSet,ierr))
      PetscCall(ISGetIndicesF90(setIS,setID,ierr))
      Allocate(MEF90MatPropBag(numSet))
      Do set = 1,numSet
         Write(setName,100) setID(set)
         Write(setprefix,101) setID(set)
         If (MEF90GlobalOptions%verbose > 0) Then
            Write(IOBuffer,102) setID(set),trim(setprefix)
            PetscCall(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
         End If
         PetscCall(PetscBagCreate(PETSC_COMM_WORLD,sizeofMEF90MatProp3D,MEF90MatPropBag(set),ierr))
         PetscCall(PetscBagRegisterMEF90MatProp(MEF90MatPropBag(set),setName,setprefix,defaultMaterial,ierr))
         If (MEF90GlobalOptions%verbose > 0) Then
            PetscCall(PetscBagView(MEF90MatPropBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr))
            PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n",ierr))
         End If
      End Do
      PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      PetscCall(ISDestroy(setIS,ierr))
100 Format('Cell set ',I4)
101 Format('cs',I4.4,'_')
102 Format('Registering materials properties for cell set ', I4,': ',A,'\n')
   End Subroutine MEF90MatPropBagSetFromOptions3D


!!! Subroutine generating various types of Hooke's laws
#undef __FUNCT__
#define __FUNCT__ "MEF90HookeLawIsoLambdaMu2D"
   Subroutine MEF90HookeLawIsoLambdaMu2D(A,lambda,mu)
      Type(Tens4OS2D),Intent(OUT)         :: A
      PetscReal,Intent(IN)                :: lambda,mu
      A = 0.0_Kr

      A%XXXX = lambda + 2.0_Kr * mu
      A%XXYY = lambda
      A%XYXY = mu
      A%YYYY = lambda + 2.0_Kr * mu
   End Subroutine MEF90HookeLawIsoLambdaMu2D

#undef __FUNCT__
#define __FUNCT__ "MEF90HookeLawIsoEnu2DPlaneStress"
   Subroutine MEF90HookeLawIsoEnu2DPlaneStress(A,E,nu)
      PetscReal,Intent(IN)                :: E,nu
      Type(Tens4OS2D),Intent(OUT)         :: A

      PetscReal                           :: Lambda,mu

      lambda = E * nu / (1.0_Kr - nu**2)
      mu     = E / (1.0_Kr + nu) * .5_Kr
      A = 0.0_Kr
      A%XXXX = lambda + 2.0_Kr * mu
      A%XXYY = lambda
      A%XYXY = mu
      A%YYYY = lambda + 2.0_Kr * mu
   End Subroutine MEF90HookeLawIsoEnu2DPlaneStress

#undef __FUNCT__
#define __FUNCT__ "MEF90HookeLawIsoEnu2DPlaneStrain"
   Subroutine MEF90HookeLawIsoEnu2DPlaneStrain(A,E,nu)
      PetscReal,Intent(IN)                :: E,nu
      Type(Tens4OS2D),Intent(OUT)         :: A

      PetscReal                           :: Lambda,mu

      lambda = E * nu / (1.0_Kr + nu) / (1.0_Kr - 2.0_Kr * nu)
      mu     = E / (1.0_Kr + nu) * .5_Kr
      A = 0.0_Kr
      A%XXXX = lambda + 2.0_Kr * mu
      A%XXYY = lambda
      A%XYXY = mu
      A%YYYY = lambda + 2.0_Kr * mu
   End Subroutine MEF90HookeLawIsoEnu2DPlaneStrain

#undef __FUNCT__
#define __FUNCT__ "MEF90HookeLawIsoLambdaMu3D"
   Subroutine MEF90HookeLawIsoLambdaMu3D(A,lambda,mu)
      PetscReal,Intent(IN)                :: lambda,mu
      Type(Tens4OS3D),Intent(OUT)         :: A

      A = 0.0_Kr
      A%XXXX = lambda + mu * 2.0_Kr
      A%XXYY = lambda
      A%XXZZ = lambda

      A%XYXY = mu

      A%XZXZ = mu

      A%YYYY = lambda + mu * 2.0_Kr
      A%YYZZ = lambda

      A%YZYZ = mu

      A%ZZZZ = lambda + mu * 2.0_Kr
   End Subroutine MEF90HookeLawIsoLambdaMu3D

#undef __FUNCT__
#define __FUNCT__ "MEF90HookeLawIsoENu3D"
   Subroutine MEF90HookeLawIsoENu3D(A,E,nu)
      PetscReal,Intent(IN)                :: E,nu
      Type(Tens4OS3D),Intent(OUT)         :: A

      Real(Kind = Kr)                     :: Lambda,mu

      lambda = E * nu / (1.0_Kr + nu) / (1 - 2.0_Kr * nu)
      mu     = E / (1.0_Kr + nu) * .5_Kr

      A = 0.0_Kr
      A%XXXX = lambda + mu * 2.0_Kr
      A%XXYY = lambda
      A%XXZZ = lambda

      A%XYXY = mu

      A%XZXZ = mu

      A%YYYY = lambda + mu * 2.0_Kr
      A%YYZZ = lambda

      A%YZYZ = mu

      A%ZZZZ = lambda + mu * 2.0_Kr
   End Subroutine MEF90HookeLawIsoENu3D

!!! Overloading linear algebra functions with Hookes Laws.
#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLaw2DSum"
!!!
!!!
!!!  MEF90HookesLaw2DSum:
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!

   Function MEF90HookesLaw2DSum(A,B)
      Type(MEF90HookesLaw2D), Intent(IN)           :: A,B
      Type(MEF90HookesLaw2D)                       :: MEF90HookesLaw2DSum

      Character(len=MEF90MXSTRLEN)                 :: IOBuffer
      PetscErrorCode                               :: ierr

      If ((A%type == MEF90HookesLawTypeIsotropic) .AND. (B%type == MEF90HookesLawTypeIsotropic)) Then
         MEF90HookesLaw2DSum%type             = MEF90HookesLawTypeIsotropic
         If (A%isPlaneStress .EQV. B%isPlaneStress) Then
            MEF90HookesLaw2DSum%lambda        = A%lambda + B%lambda
            MEF90HookesLaw2DSum%mu            = A%mu + B%mu
            MEF90HookesLaw2DSum%isPlaneStress = A%isPlaneStress
         Else
            Write(IOBuffer,*) "Incompatible planar Hooke law type in "//__FUNCT__//'\n'
            PetscCall(PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr))
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,IOBuffer)
         End If
         If (A%isPlaneStress) Then
            MEF90HookesLaw2DSum%PoissonRatio  = MEF90HookesLaw2DSum%lambda / (MEF90HookesLaw2DSum%lambda + MEF90HookesLaw2DSum%mu) * 0.5_Kr
            MEF90HookesLaw2DSum%YoungsModulus = 2.0_Kr * MEF90HookesLaw2DSum%mu * (1.0_Kr + MEF90HookesLaw2DSum%PoissonRatio)
            MEF90HookesLaw2DSum%BulkModulus   = MEF90HookesLaw2DSum%lambda + MEF90HookesLaw2DSum%mu
            PetscCall(PetscLogFlops(9._pflop,ierr))
         Else
            MEF90HookesLaw2DSum%PoissonRatio  = MEF90HookesLaw2DSum%lambda / (MEF90HookesLaw2DSum%lambda + 2.0_Kr * MEF90HookesLaw2DSum%mu) * 0.5_Kr
            MEF90HookesLaw2DSum%YoungsModulus = 2.0_Kr * MEF90HookesLaw2DSum%mu * (1.0_Kr + MEF90HookesLaw2DSum%PoissonRatio)
            MEF90HookesLaw2DSum%BulkModulus   = MEF90HookesLaw2DSum%lambda + MEF90HookesLaw2DSum%mu
            PetscCall(PetscLogFlops(10._pflop,ierr))
         End If
      Else If ((A%type == MEF90HookesLawTypeFull) .AND. (B%type == MEF90HookesLawTypeFull)) Then
         MEF90HookesLaw2DSum%type       = MEF90HookesLawTypeFull
         MEF90HookesLaw2DSum%fullTensor = A%fullTensor + B%fullTensor
      Else
            Write(IOBuffer,*) "Incompatible planar Hooke law type in "//__FUNCT__//'\n'
            PetscCall(PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr))
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,IOBuffer)
      End If
   End Function MEF90HookesLaw2DSum

#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLaw3DSum"
!!!
!!!
!!!  MEF90HookesLaw3DSum:
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!

   Function MEF90HookesLaw3DSum(A,B)
      Type(MEF90HookesLaw3D), Intent(IN)           :: A,B
      Type(MEF90HookesLaw3D)                       :: MEF90HookesLaw3DSum

      Character(len=MEF90MXSTRLEN)                 :: IOBuffer
      PetscErrorCode                               :: ierr

      If ((A%type == MEF90HookesLawTypeIsotropic) .AND. (B%type == MEF90HookesLawTypeIsotropic)) Then
         MEF90HookesLaw3DSum%type          = MEF90HookesLawTypeIsotropic
         MEF90HookesLaw3DSum%lambda        = A%lambda + B%lambda
         MEF90HookesLaw3DSum%mu            = A%mu + B%mu
         MEF90HookesLaw3DSum%PoissonRatio  = MEF90HookesLaw3DSum%lambda / (MEF90HookesLaw3DSum%lambda + MEF90HookesLaw3DSum%mu) * 0.5_Kr
         MEF90HookesLaw3DSum%YoungsModulus = MEF90HookesLaw3DSum%mu * (3.0_Kr * MEF90HookesLaw3DSum%lambda + 2.0_Kr * MEF90HookesLaw3DSum%mu) / (MEF90HookesLaw3DSum%lambda + MEF90HookesLaw3DSum%mu)
         MEF90HookesLaw3DSum%BulkModulus   = MEF90HookesLaw3DSum%lambda + MEF90HookesLaw3DSum%mu * 2.0_Kr / 3.0_Kr
         PetscCall(PetscLogFlops(14._pflop,ierr))
      Else If ((A%type == MEF90HookesLawTypeFull) .AND. (B%type == MEF90HookesLawTypeFull)) Then
         MEF90HookesLaw3DSum%type       = MEF90HookesLawTypeFull
         MEF90HookesLaw3DSum%fullTensor = A%fullTensor + B%fullTensor
      Else
         Write(IOBuffer,*) "Incompatible planar Hooke law type in "//__FUNCT__//'\n'
         PetscCall(PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr))
         SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,IOBuffer)
      End If
   End Function MEF90HookesLaw3DSum

#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLaw2DDiff"
!!!
!!!
!!!  MEF90HookesLaw2DDiff:
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!

   Function MEF90HookesLaw2DDiff(A,B)
      Type(MEF90HookesLaw2D), Intent(IN)           :: A,B
      Type(MEF90HookesLaw2D)                       :: MEF90HookesLaw2DDiff

      Character(len=MEF90MXSTRLEN)                 :: IOBuffer
      PetscErrorCode                               :: ierr

      If ((A%type == MEF90HookesLawTypeIsotropic) .AND. (B%type == MEF90HookesLawTypeIsotropic)) Then
         MEF90HookesLaw2DDiff%type             = MEF90HookesLawTypeIsotropic
         If (A%isPlaneStress .EQV. B%isPlaneStress) Then
            MEF90HookesLaw2DDiff%lambda        = A%lambda - B%lambda
            MEF90HookesLaw2DDiff%mu            = A%mu - B%mu
            MEF90HookesLaw2DDiff%isPlaneStress = A%isPlaneStress
         Else
            Write(IOBuffer,*) "Incompatible planar Hooke law type in "//__FUNCT__//'\n'
            PetscCall(PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr))
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,IOBuffer)
         End If
         If (A%isPlaneStress) Then
            MEF90HookesLaw2DDiff%PoissonRatio  = MEF90HookesLaw2DDiff%lambda / (MEF90HookesLaw2DDiff%lambda + MEF90HookesLaw2DDiff%mu) * 0.5_Kr
            MEF90HookesLaw2DDiff%YoungsModulus = 2.0_Kr * MEF90HookesLaw2DDiff%mu * (1.0_Kr + MEF90HookesLaw2DDiff%PoissonRatio)
            MEF90HookesLaw2DDiff%BulkModulus   = MEF90HookesLaw2DDiff%lambda + MEF90HookesLaw2DDiff%mu
            PetscCall(PetscLogFlops(9._pflop,ierr))
         Else
            MEF90HookesLaw2DDiff%PoissonRatio  = MEF90HookesLaw2DDiff%lambda / (MEF90HookesLaw2DDiff%lambda + 2.0_Kr * MEF90HookesLaw2DDiff%mu) * 0.5_Kr
            MEF90HookesLaw2DDiff%YoungsModulus = 2.0_Kr * MEF90HookesLaw2DDiff%mu * (1.0_Kr + MEF90HookesLaw2DDiff%PoissonRatio)
            MEF90HookesLaw2DDiff%BulkModulus   = MEF90HookesLaw2DDiff%lambda + MEF90HookesLaw2DDiff%mu
            PetscCall(PetscLogFlops(10._pflop,ierr))
         End If
      Else If ((A%type == MEF90HookesLawTypeFull) .AND. (B%type == MEF90HookesLawTypeFull)) Then
         MEF90HookesLaw2DDiff%type       = MEF90HookesLawTypeFull
         MEF90HookesLaw2DDiff%fullTensor = A%fullTensor - B%fullTensor
      Else
            Write(IOBuffer,*) "Incompatible planar Hooke law type in "//__FUNCT__//'\n'
            PetscCall(PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr))
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,IOBuffer)
      End If
   End Function MEF90HookesLaw2DDiff

#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLaw3DDiff"
!!!
!!!
!!!  MEF90HookesLaw3DDiff:
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!

   Function MEF90HookesLaw3DDiff(A,B)
      Type(MEF90HookesLaw3D), Intent(IN)           :: A,B
      Type(MEF90HookesLaw3D)                       :: MEF90HookesLaw3DDiff

      Character(len=MEF90MXSTRLEN)                 :: IOBuffer
      PetscErrorCode                               :: ierr

      If ((A%type == MEF90HookesLawTypeIsotropic) .AND. (B%type == MEF90HookesLawTypeIsotropic)) Then
         MEF90HookesLaw3DDiff%type          = MEF90HookesLawTypeIsotropic
         MEF90HookesLaw3DDiff%lambda        = A%lambda - B%lambda
         MEF90HookesLaw3DDiff%mu            = A%mu - B%mu
         MEF90HookesLaw3DDiff%PoissonRatio  = MEF90HookesLaw3DDiff%lambda / (MEF90HookesLaw3DDiff%lambda + MEF90HookesLaw3DDiff%mu) * 0.5_Kr
         MEF90HookesLaw3DDiff%YoungsModulus = MEF90HookesLaw3DDiff%mu * (3.0_Kr * MEF90HookesLaw3DDiff%lambda + 2.0_Kr * MEF90HookesLaw3DDiff%mu) / (MEF90HookesLaw3DDiff%lambda + MEF90HookesLaw3DDiff%mu)
         MEF90HookesLaw3DDiff%BulkModulus   = MEF90HookesLaw3DDiff%lambda + MEF90HookesLaw3DDiff%mu * 2.0_Kr / 3.0_Kr
         PetscCall(PetscLogFlops(14._pflop,ierr))
      Else If ((A%type == MEF90HookesLawTypeFull) .AND. (B%type == MEF90HookesLawTypeFull)) Then
         MEF90HookesLaw3DDiff%type       = MEF90HookesLawTypeFull
         MEF90HookesLaw3DDiff%fullTensor = A%fullTensor - B%fullTensor
      Else
         Write(IOBuffer,*) "Incompatible planar Hooke law type in "//__FUNCT__//'\n'
         PetscCall(PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr))
         SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,IOBuffer)
      End If
   End Function MEF90HookesLaw3DDiff

#undef __FUNCT__
#define __FUNCT__ "ScalarXMEF90HookesLaw2D"
!!!
!!!
!!!  ScalarXMEF90HookesLaw2D:
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!

   Function ScalarXMEF90HookesLaw2D(t,A)
      PetscReal, Intent(IN)                        :: t
      Type(MEF90HookesLaw2D), Intent(IN)           :: A
      Type(MEF90HookesLaw2D)                       :: ScalarXMEF90HookesLaw2D

      PetscErrorCode                               :: ierr

      If (A%type == MEF90HookesLawTypeIsotropic) Then
         ScalarXMEF90HookesLaw2D%type          = MEF90HookesLawTypeIsotropic
         ScalarXMEF90HookesLaw2D%lambda        = t * A%lambda
         ScalarXMEF90HookesLaw2D%mu            = t * A%mu
         ScalarXMEF90HookesLaw2D%isPlaneStress = A%isPlaneStress
         If (A%isPlaneStress) Then
            ScalarXMEF90HookesLaw2D%PoissonRatio  = ScalarXMEF90HookesLaw2D%lambda / (ScalarXMEF90HookesLaw2D%lambda + ScalarXMEF90HookesLaw2D%mu) * 0.5_Kr
            ScalarXMEF90HookesLaw2D%YoungsModulus = 2.0_Kr * ScalarXMEF90HookesLaw2D%mu * (1.0_Kr + ScalarXMEF90HookesLaw2D%PoissonRatio)
            ScalarXMEF90HookesLaw2D%BulkModulus   = ScalarXMEF90HookesLaw2D%lambda + ScalarXMEF90HookesLaw2D%mu
            PetscCall(PetscLogFlops(9._pflop,ierr))
         Else
            ScalarXMEF90HookesLaw2D%PoissonRatio  = ScalarXMEF90HookesLaw2D%lambda / (ScalarXMEF90HookesLaw2D%lambda + 2.0_Kr * ScalarXMEF90HookesLaw2D%mu) * 0.5_Kr
            ScalarXMEF90HookesLaw2D%YoungsModulus = 2.0_Kr * ScalarXMEF90HookesLaw2D%mu * (1.0_Kr + ScalarXMEF90HookesLaw2D%PoissonRatio)
            ScalarXMEF90HookesLaw2D%BulkModulus   = ScalarXMEF90HookesLaw2D%lambda + ScalarXMEF90HookesLaw2D%mu
            PetscCall(PetscLogFlops(10._pflop,ierr))
         End If
      Else 
         ScalarXMEF90HookesLaw2D%type       = MEF90HookesLawTypeFull
         ScalarXMEF90HookesLaw2D%fullTensor = t * A%fullTensor
      End If
   End Function ScalarXMEF90HookesLaw2D


#undef __FUNCT__
#define __FUNCT__ "ScalarXMEF90HookesLaw3D"
!!!
!!!
!!!  ScalarXMEF90HookesLaw3D:
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!

   Function ScalarXMEF90HookesLaw3D(t,A)
      PetscReal, Intent(IN)                        :: t
      Type(MEF90HookesLaw3D), Intent(IN)           :: A
      Type(MEF90HookesLaw3D)                       :: ScalarXMEF90HookesLaw3D

      PetscErrorCode                               :: ierr

      If (A%type == MEF90HookesLawTypeIsotropic) Then
         ScalarXMEF90HookesLaw3D%type          = MEF90HookesLawTypeIsotropic
         ScalarXMEF90HookesLaw3D%lambda        = t * A%lambda
         ScalarXMEF90HookesLaw3D%mu            = t * A%mu
         ScalarXMEF90HookesLaw3D%PoissonRatio  = ScalarXMEF90HookesLaw3D%lambda / (ScalarXMEF90HookesLaw3D%lambda + ScalarXMEF90HookesLaw3D%mu) * 0.5_Kr
         ScalarXMEF90HookesLaw3D%YoungsModulus = ScalarXMEF90HookesLaw3D%mu * (3.0_Kr * ScalarXMEF90HookesLaw3D%lambda + 2.0_Kr * ScalarXMEF90HookesLaw3D%mu) / (ScalarXMEF90HookesLaw3D%lambda + ScalarXMEF90HookesLaw3D%mu)
         ScalarXMEF90HookesLaw3D%BulkModulus   = ScalarXMEF90HookesLaw3D%lambda + ScalarXMEF90HookesLaw3D%mu * 2.0_Kr / 3.0_Kr
         PetscCall(PetscLogFlops(14._pflop,ierr))
      Else 
         ScalarXMEF90HookesLaw3D%type       = MEF90HookesLawTypeFull
         ScalarXMEF90HookesLaw3D%fullTensor = t * A%fullTensor
      End If
   End Function ScalarXMEF90HookesLaw3D

#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLaw2DXMatS2D"
!!!
!!!
!!!  MEF90HookesLaw2DXMatS2D:
!!!
!!!  (c) 2016 Blaise Bourdin bourdin@lsu.edu
!!!

   Function MEF90HookesLaw2DXMatS2D(A,X)
      Type(MEF90HookesLaw2D), Intent(IN)           :: A
      Type(MatS2D), Intent(IN)                     :: X
      Type(MatS2D)                                 :: MEF90HookesLaw2DXMatS2D

      PetscErrorCode                               :: ierr
      Real(Kind = Kr)                              :: C1, C2

      Select case(A%type)
         Case(MEF90HookesLawTypeIsotropic)
            C1 = A%lambda + 2.0_Kr * A%mu
            C2 = 2.0_Kr * A%mu
            MEF90HookesLaw2DXMatS2D%XX = C1 * X%XX       + A%lambda * X%YY
            MEF90HookesLaw2DXMatS2D%YY = A%lambda * X%XX + C1 * X%YY
            MEF90HookesLaw2DXMatS2D%XY = C2 * X%XY
            PetscCall(PetscLogFlops(10._pflop,ierr))
         Case(MEF90HookesLawTypeFull)
            MEF90HookesLaw2DXMatS2D = A%fullTensor * X
            ! flops are counted in m_MEF90_LinAlg
      End Select
   End Function MEF90HookesLaw2DXMatS2D

#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLaw2DXMatS3D"
!!!
!!!
!!!  MEF90HookesLaw2DXMatS3D:
!!!
!!!  (c) 2016 Blaise Bourdin bourdin@lsu.edu
!!!

   Function MEF90HookesLaw3DXMatS3D(A,X)
      Type(MEF90HookesLaw3D), Intent(IN)           :: A
      Type(MatS3D), Intent(IN)                     :: X
      Type(MatS3D)                                 :: MEF90HookesLaw3DXMatS3D

      PetscErrorCode                               :: ierr
      Real(Kind = Kr)                              :: C1, C2

      Select case(A%type)
         Case(MEF90HookesLawTypeIsotropic)
            C1 = A%lambda + 2.0_Kr * A%mu
            C2 = 2.0_Kr * A%mu
            MEF90HookesLaw3DXMatS3D%XX = C1 * X%XX       + A%lambda * X%YY + A%lambda * X%ZZ
            MEF90HookesLaw3DXMatS3D%YY = A%lambda * X%XX + C1 * X%YY       + A%lambda * X%ZZ
            MEF90HookesLaw3DXMatS3D%ZZ = A%lambda * X%XX + A%lambda * X%YY + C1 * X%ZZ
            MEF90HookesLaw3DXMatS3D%YZ = C2 * X%YZ
            MEF90HookesLaw3DXMatS3D%XZ = C2 * X%XZ
            MEF90HookesLaw3DXMatS3D%XY = C2 * X%XY

            PetscCall(PetscLogFlops(21._pflop,ierr))
         Case(MEF90HookesLawTypeFull)
            MEF90HookesLaw3DXMatS3D = A%fullTensor * X
            ! flops are counted in m_MEF90_LinAlg
      End Select
   End Function MEF90HookesLaw3DXMatS3D

#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLaw2DXMat2D"
!!!
!!!
!!!  MEF90HookesLaw2DXMat2D:
!!!
!!!  (c) 2016 Blaise Bourdin bourdin@lsu.edu
!!!
   Function MEF90HookesLaw2DXMat2D(A,X)
      Type(MEF90HookesLaw2D), Intent(IN)           :: A
      Type(Mat2D), Intent(IN)                      :: X
      Type(Mat2D)                                  :: MEF90HookesLaw2DXMat2D

      PetscErrorCode                               :: ierr
      Real(Kind = Kr)                              :: C1, C2

      Select case(A%type)
         Case(MEF90HookesLawTypeIsotropic)
            C1 = A%lambda + 2.0_Kr * A%mu
            C2 = 2.0_Kr * A%mu
            MEF90HookesLaw2DXMat2D%XX = C1 * X%XX       + A%lambda * X%YY
            MEF90HookesLaw2DXMat2D%XY = C2 * X%XY
            MEF90HookesLaw2DXMat2D%YY = A%lambda * X%XX + C1 * X%YY
            MEF90HookesLaw2DXMat2D%YX = C2 * X%YX
            PetscCall(PetscLogFlops(11._pflop,ierr))
         Case(MEF90HookesLawTypeFull)
            MEF90HookesLaw2DXMat2D = A%fullTensor * X
            ! flops are counted in m_MEF90_LinAlg
      End Select
   End Function MEF90HookesLaw2DXMat2D

#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLaw2DXMat3D"
!!!
!!!
!!!  MEF90HookesLaw2DXMat3D:
!!!
!!!  (c) 2016 Blaise Bourdin bourdin@lsu.edu
!!!
   Function MEF90HookesLaw3DXMat3D(A,X)
      Type(MEF90HookesLaw3D), Intent(IN)           :: A
      Type(Mat3D), Intent(IN)                      :: X
      Type(Mat3D)                                  :: MEF90HookesLaw3DXMat3D

      PetscErrorCode                               :: ierr
      Real(Kind = Kr)                              :: C1, C2

      Select case(A%type)
         Case(MEF90HookesLawTypeIsotropic)
            C1 = A%lambda + 2.0_Kr * A%mu
            C2 = 2.0_Kr * A%mu
            MEF90HookesLaw3DXMat3D%XX = C1 * X%XX       + A%lambda * X%YY + A%lambda * X%ZZ
            MEF90HookesLaw3DXMat3D%XY = C2 * X%XY
            MEF90HookesLaw3DXMat3D%XZ = C2 * X%XZ

            MEF90HookesLaw3DXMat3D%YX = C2 * X%YX
            MEF90HookesLaw3DXMat3D%YY = A%lambda * X%XX + C1 * X%YY       + A%lambda * X%ZZ
            MEF90HookesLaw3DXMat3D%YZ = C2 * X%YZ

            MEF90HookesLaw3DXMat3D%XZ = C2 * X%XZ
            MEF90HookesLaw3DXMat3D%YZ = C2 * X%YZ
            MEF90HookesLaw3DXMat3D%ZZ = A%lambda * X%XX + A%lambda * X%YY + C1 * X%ZZ
            PetscCall(PetscLogFlops(24._pflop,ierr))
         Case(MEF90HookesLawTypeFull)
            MEF90HookesLaw3DXMat3D = A%fullTensor * X
            ! flops are counted in m_MEF90_LinAlg
      End Select
   End Function MEF90HookesLaw3DXMat3D

End Module m_MEF90_Materials
