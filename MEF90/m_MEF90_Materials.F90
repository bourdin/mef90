Module m_MEF90_Materials_Types
#include "finclude/petscdef.h"
   Use m_MEF90_LinAlg
   Use m_MEF90_Ctx
   Use petsc
   IMPLICIT NONE

   Type MEF90HookesLaw2D
      Sequence
      Type(Tens4OS2D)    :: fullTensor
      PetscReal          :: lambda,mu,YoungsModulus,PoissonRatio,BulkModulus
      PetscEnum          :: type
      PetscBool          :: isPlaneStress
   End Type MEF90HookesLaw2D

   Type MEF90HookesLaw3D
      Sequence
      Type(Tens4OS3D)    :: fullTensor
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

   Type MEF90MatProp2D_Type
      Sequence
      PetscReal                     :: Density                                          ! rho
      PetscReal                     :: FractureToughness                                ! Gc
      Type(MatS2D)                  :: toughnessAnisotropyMatrix                        ! 
      PetscReal                     :: SpecificHeat                                     ! Cp
      Type(MatS2D)                  :: ThermalConductivity                              ! K
      Type(MatS2D)                  :: LinearThermalExpansion                           ! alpha
      Type(MEF90HookesLaw2D)        :: HookesLaw                                        ! A
      PetscReal                     :: internalLength                                   ! l
      PetscReal                     :: CoefficientLinSoft                               ! k
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
      PetscReal                     :: drivingForceBeta                                 ! beta parameter in Drucker-Prager driving force
      PetscInt                      :: drivingForceP                                    ! p parameter in Drucker-Prager driving force
      PetscBool                     :: isLinearIsotropicHardening
      PetscBool                     :: isNoPlCoupling
      Character(len=MEF90_MXSTRLEN) :: Name
   End Type MEF90MatProp2D_Type

   Type MEF90MatProp3D_Type
      Sequence
      PetscReal                     :: Density                                          ! rho
      PetscReal                     :: FractureToughness                                ! Gc
      Type(MatS3D)                  :: toughnessAnisotropyMatrix                        ! 
      PetscReal                     :: SpecificHeat                                     ! Cp
      Type(MatS3D)                  :: ThermalConductivity                              ! K
      Type(MatS3D)                  :: LinearThermalExpansion                           ! alpha
      Type(MEF90HookesLaw3D)        :: HookesLaw                                        ! A
      PetscReal                     :: internalLength                                   ! l
      PetscReal                     :: CoefficientLinSoft                               ! k
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
      PetscReal                     :: drivingForceBeta                                 ! beta parameter in Drucker-Prager driving force
      PetscInt                      :: drivingForceP                                    ! p parameter in Drucker-Prager driving force
      PetscBool                     :: isLinearIsotropicHardening
      PetscBool                     :: isNoPlCoupling
      Character(len=MEF90_MXSTRLEN) :: Name
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
         0.0_Kr,0.0_Kr,1.0_Kr,.3_Kr,0.0_Kr,                                            & ! lambda, mu, E, nu, kappa (lambda, mu, kappa will be recomputed)
         MEF90HookesLawTypeIsotropic,                                                  & ! type
         .FALSE.),                                                                     & ! isPlaneStress
      1.0_Kr,                                                                          & ! Internal Length
      0.0_Kr,                                                                          & ! CoefficientLinSoft
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
      "MEF90Mathium2D")

   Type(MEF90MatProp3D_Type),Parameter     :: MEF90Mathium3D = MEF90MatProp3D_Type(    &
      1.0_Kr,                                                                          & ! Density
      1.0_Kr,                                                                          & ! FractureToughness
      MEF90MatS3DIdentity,                                                             & ! toughnessAnisotropyMatrix
      1.0_Kr,                                                                          & ! SpecificHeat
      MEF90MatS3DIdentity,                                                             & ! ThermalConductivity
      MEF90MatS3DIdentity,                                                             & ! LinearThermalExpansion
      MEF90HookesLaw3D(                                                                &
         Tens4OS3D(1.34615_Kr,0.57692_Kr,0.57692_Kr,0.00000_Kr,0.00000_Kr,0.00000_Kr,  & ! XXXX,XXYY,XXZZ,XXYZ,XXXZ,XXXY
                              1.34615_Kr,0.57692_Kr,0.00000_Kr,0.00000_Kr,0.00000_Kr,  & !      YYYY,YYZZ,YYYZ,YYXZ,YYXY
                                         1.34615_Kr,0.00000_Kr,0.00000_Kr,0.00000_Kr,  & !           ZZZZ ZZYZ,ZZXZ,ZZXY
                                                    0.38462_Kr,0.00000_Kr,0.00000_Kr,  & !                YXYX,YZXZ,YZXY
                                                               0.38462_Kr,0.00000_Kr,  & !                     XZXZ,XZXY
                                                                          0.38462_Kr), & !                          XYXY
         0.0_Kr,0.0_Kr,1.0_Kr,.3_Kr,0.0_Kr,                                            & ! lambda, mu, E, nu, kappa (lambda, mu, kappa will be recomputed)
         MEF90HookesLawTypeIsotropic,                                                  & ! type
#if (PETSC_SIZEOF_INT == 4)
         PETSC_FALSE                                                                   & ! padding
#endif
      ),                                                                               & ! 
      1.0_Kr,                                                                          & ! Internal Length
      2.0_Kr,                                                                          & ! CoefficientLinSoft
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
      "MEF90Mathium3D")
End Module m_MEF90_Materials_Types

Module m_MEF90_Materials_Interface2D
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
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

      Call PetscBagGetData(bag,data,ierr)
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
   End Subroutine PetscBagGetDataMEF90MatProp2D
End Module m_MEF90_Materials_Interface2D

Module m_MEF90_Materials_Interface3D
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
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

      Call PetscBagGetData(bag,data,ierr)
      Select case(data%HookesLaw%type)
         Case(MEF90HookesLawTypeFull)
            Continue
         Case(MEF90HookesLawTypeIsotropic)
            data%HookesLaw%lambda       = data%HookesLaw%YoungsModulus * data%HookesLaw%PoissonRatio / (1.0_Kr + data%HookesLaw%PoissonRatio) / (1.0_Kr - 2.0_Kr * data%HookesLaw%PoissonRatio)
            data%HookesLaw%mu           = data%HookesLaw%YoungsModulus / (1.0_Kr + data%HookesLaw%PoissonRatio) * .5_Kr
            data%HookesLaw%BulkModulus  = data%HookesLaw%lambda + 2.0_Kr * data%HookesLaw%mu / 3.0_Kr
      End Select
   End Subroutine PetscBagGetDataMEF90MatProp3D
End Module m_MEF90_Materials_Interface3D

Module m_MEF90_Materials
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
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

   Interface  Operator (*)
      Module Procedure MEF90HookesLaw2DXMatS2D,MEF90HookesLaw3DXMatS3D,MEF90HookesLaw2DXMat2D,MEF90HookesLaw3DXMat3D
   End Interface

   Interface MasonryProjection
      Module Procedure MasonryProjection2D,MasonryProjection3D
   End Interface

   Interface HDProjection
      Module Procedure HDProjection2D,HDProjection3D
   End Interface

   PetscSizeT,protected   :: sizeofMEF90MatProp2D
   PetscSizeT,protected   :: sizeofMEF90MatProp3D

   PetscSizeT,protected   :: sizeofMEF90HookesLaw2D
   PetscSizeT,protected   :: sizeofMEF90HookesLaw3D

   Character(len =  MEF90_MXSTRLEN),Dimension(5),protected   :: MEF90HookesLawTypeList

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

      Call PetscDataTypeGetSize(PETSC_CHAR,sizeofchar,ierr)
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

      Call PetscBagGetDataMEF90MatProp2D(bag,matprop,ierr)
      Call PetscBagSetName(bag,trim(name),"MatProp2D object: material properties",ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix), ierr)

      Call PetscBagRegisterString(bag,matprop%name,trim(default%name),'Name','',ierr)
      Call PetscBagRegisterReal(bag,matprop%density,default%density,'Density','[kg.m^(-2)] (rho) Density',ierr)
      Call PetscBagRegisterReal(bag,matprop%FractureToughness,default%FractureToughness,'FractureToughness','[N.m^(-1)] (G_c) Fracture toughness',ierr)
      matprop%toughnessAnisotropyMatrix = default%toughnessAnisotropyMatrix
      Call PetscBagRegisterRealArray(bag,matprop%toughnessAnisotropyMatrix,3,'toughnessAnisotropyMatrix','[] toughness Anisotropy Matrix',ierr)
      Call PetscBagRegisterReal(bag,matprop%SpecificHeat,default%SpecificHeat,'SpecificHeat','[J.kg^(-1).K^(-1)] (Cp) Specific heat',ierr)
      matprop%ThermalConductivity = default%ThermalConductivity
      Call PetscBagRegisterRealArray(bag,matprop%ThermalConductivity,3,'ThermalConductivity','[J.m^(-1).s^(-1).K^(-1)] (K) Thermal conductivity',ierr)
      matprop%LinearThermalExpansion = default%LinearThermalExpansion
      Call PetscBagRegisterRealArray(bag,matprop%LinearThermalExpansion,3,'LinearThermalExpansion','[K^(-1)] (alpha) Linear thermal expansion matrix',ierr)

      Call PetscBagRegisterEnum(bag,matprop%HookesLaw%type,MEF90HookesLawTypeList,default%HookesLaw%type,'hookeslaw_type','Type of Hooke''s law',ierr);CHKERRQ(ierr)
      Select case(matprop%HookesLaw%type)
         Case (MEF90HookesLawTypeFull)
            matprop%HookesLaw%fullTensor = default%HookesLaw%fullTensor
            Call PetscBagRegisterRealArray(bag,matprop%HookesLaw%fullTensor,6,'HookesLaw_tensor','[N.m^(-2)] (A) Hooke''s law',ierr)
         Case(MEF90HookesLawTypeIsotropic)
            Call PetscBagRegisterReal(bag,matprop%HookesLaw%YoungsModulus,default%HookesLaw%YoungsModulus,'hookeslaw_YoungsModulus','[N.m^(-2)] (E) Young''s Modulus',ierr)
            Call PetscBagRegisterReal(bag,matprop%HookesLaw%PoissonRatio,default%HookesLaw%PoissonRatio,'hookeslaw_PoissonRatio','[] (nu) Poisson Modulus',ierr)
            Call PetscBagRegisterBool(bag,matprop%HookesLaw%isPlaneStress,default%HookesLaw%isPlaneStress,'hookeslaw_planeStress','Use plane stress elasticity',ierr);CHKERRQ(ierr)
            matprop%HookesLaw%fulltensor = -1.D+30
      End Select

      Call PetscBagRegisterReal(bag,matprop%internalLength,default%internalLength,'internalLength','[m] (l) Internal Length',ierr)
      Call PetscBagRegisterReal(bag,matprop%CoefficientLinSoft,default%CoefficientLinSoft,'CoefficientLinSoft','[] (k) Linear softening coefficient for LinSoft',ierr)

      Call PetscBagRegisterReal(bag,matprop%yieldStress,default%yieldStress,'yieldStress','[N.m^(-2)] (sigma_y) stress threshold for plasticity',ierr)
      Call PetscBagRegisterReal(bag,matprop%residualYieldStress,default%residualYieldStress,'residualyieldStress','[unit-less] (eta) residual yield stress',ierr)
      Call PetscBagRegisterReal(bag,matprop%DuctileCouplingPower,default%DuctileCouplingPower,'DuctileCouplingPower','[] power of the coupling between the damage and the plasticity',ierr)
      Call PetscBagRegisterReal(bag,matprop%CoefficientCapModel0,default%CoefficientCapModel0,'CoefficientCapModel0','C0 in the Yield function: CD || dev(stress) || - C2 tr(stress)^2 - C1 tr(stress) - C0 <= 0',ierr)
      Call PetscBagRegisterReal(bag,matprop%CoefficientCapModel1,default%CoefficientCapModel1,'CoefficientCapModel1','C1 in the Yield function: CD || dev(stress) || - C2 tr(stress)^2 - C1 tr(stress) - C0 <= 0',ierr)
      Call PetscBagRegisterReal(bag,matprop%CoefficientCapModel2,default%CoefficientCapModel2,'CoefficientCapModel2','C2 in the Yield function: CD || dev(stress) || - C2 tr(stress)^2 - C1 tr(stress) - C0 <= 0',ierr)
      Call PetscBagRegisterReal(bag,matprop%CoefficientCapModelD,default%CoefficientCapModelD,'CoefficientCapModelD','CD in the Yield function: CD || dev(stress) || - C2 tr(stress)^2 - C1 tr(stress) - C0 <= 0',ierr)
      Call PetscBagRegisterReal(bag,matprop%CoefficientDruckerPrager,default%CoefficientDruckerPrager,'CoefficientDruckerPrager','k in the Yield function: || dev(stress) || - k tr(stress) - yieldStress <= 0',ierr)
      Call PetscBagRegisterReal(bag,matprop%CoeffF,default%CoeffF,'CoeffF','[unit-less] (F) coefficient F in the Hill yield criterion',ierr)
      Call PetscBagRegisterReal(bag,matprop%CoeffG,default%CoeffG,'CoeffG','[unit-less] (G) coefficient G in the Hill yield criterion',ierr)
      Call PetscBagRegisterReal(bag,matprop%CoeffH,default%CoeffH,'CoeffH','[unit-less] (H) coefficient H in the Hill yield criterion',ierr)
      Call PetscBagRegisterReal(bag,matprop%CoeffM,default%CoeffM,'CoeffM','[unit-less] (M) coefficient M in the Hill yield criterion',ierr)
      Call PetscBagRegisterReal(bag,matprop%CoeffN,default%CoeffN,'CoeffN','[unit-less] (N) coefficient N in the Hill yield criterion',ierr)
      Call PetscBagRegisterReal(bag,matprop%CoeffL,default%CoeffL,'CoeffL','[unit-less] (L) coefficient L in the Hill yield criterion',ierr)
      Call PetscBagRegisterReal(bag,matprop%YieldTau0,default%YieldTau0,'YieldTau0','[N.m^(-2)] (tau_0) stress threshold in the Hill yield criterion',ierr)
      Call PetscBagRegisterReal(bag,matprop%residualYieldTau0,default%residualYieldTau0,'residualYieldTau0','[unit-less] residual stress threshold in the Hill yield criterion',ierr)
      Call PetscBagRegisterReal(bag,matprop%phi1,default%phi1,'phi1','[radians] Bunge-Euler angle in the Hill yield criterion',ierr)
      Call PetscBagRegisterReal(bag,matprop%phi2,default%phi2,'phi2','[radians] Bunge-Euler angle in the Hill yield criterion',ierr)
      Call PetscBagRegisterReal(bag,matprop%Phi,default%Phi,'Phi','[radians] Bunge-Euler angle in the Hill yield criterion',ierr)
      Call PetscBagRegisterReal(bag,matprop%delta,default%delta,'delta','[unit-less] residual in the definition of the porosity, Gurson and Green criteria',ierr)

      Call PetscBagRegisterReal(bag,matprop%cohesiveStiffness,default%cohesiveStiffness,'cohesiveStiffness','[N.m^(-4)] (k) cohesive stiffness in Winkler-type models',ierr)
      Call PetscBagRegisterReal(bag,matprop%residualStiffness,default%residualStiffness,'residualStiffness','[unit-less] (eta) residual stiffness',ierr)

      Call PetscBagRegisterReal(bag,matprop%drivingForceTensileStrength,default%drivingForceTensileStrength,'drivingForce_tensileStrength','[N.m^(-2)] (\sigma_{ts}) tensile strength in Drucker-Prager driving Force',ierr)
      Call PetscBagRegisterReal(bag,matprop%drivingForceCompressiveStrength,default%drivingForceCompressiveStrength,'drivingForce_CompressiveStrength','[N.m^(-2)] (\sigma_{cs}) compressive strength in Drucker-Prager driving Force',ierr)
      Call PetscBagRegisterReal(bag,matprop%drivingForceBeta,default%drivingForceBeta,'drivingForce_Beta','[unit-less] (\beta) alpha parameter in Drucker-Prager driving Force',ierr)
      Call PetscBagRegisterInt(bag,matprop%drivingForcep,default%drivingForcep,'drivingForce_p','[unit-less] (p) p parameter in Drucker-Prager driving Force',ierr)

      Call PetscBagRegisterBool(bag,matprop%isLinearIsotropicHardening,default%isLinearIsotropicHardening,'isLinearIsotropicHardening','[bool] Plasticity with Linear Isotropic Hardening',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterBool(bag,matprop%isNoPlCoupling,default%isNoPlCoupling,'isNoPlCoupling','[bool] Coupling between damage and plastic dissipation',ierr);CHKERRQ(ierr)
      !Call PetscBagSetFromOptions(bag,ierr)
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
      Call PetscBagGetDataMEF90MatProp3D(bag,matprop,ierr)
      Call PetscBagSetName(bag,trim(name),"MatProp3D object: material properties",ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix), ierr)

      Call PetscBagRegisterString(bag,matprop%name,trim(default%name),'Name','',ierr)
      Call PetscBagRegisterReal(bag,matprop%density,default%density,'Density','[kg.m^(-3)] (rho) Density',ierr)
      Call PetscBagRegisterReal(bag,matprop%FractureToughness,default%FractureToughness,'FractureToughness','[N.m^(-1)] (G_c) Fracture toughness',ierr)
      matprop%toughnessAnisotropyMatrix = default%toughnessAnisotropyMatrix
      Call PetscBagRegisterRealArray(bag,matprop%toughnessAnisotropyMatrix,6,'toughnessAnisotropyMatrix','[] toughness Anisotropy Matrix',ierr)
      Call PetscBagRegisterReal(bag,matprop%SpecificHeat,default%SpecificHeat,'SpecificHeat','[J.kg^(-1).K^(-1)] (Cp) Specific heat',ierr)
      matprop%ThermalConductivity = default%ThermalConductivity
      Call PetscBagRegisterRealArray(bag,matprop%ThermalConductivity,6,'ThermalConductivity','[J.m^(-1).s^(-1).K^(-1)] (K) Thermal conductivity',ierr)
      matprop%LinearThermalExpansion = default%LinearThermalExpansion
      Call PetscBagRegisterRealArray(bag,matprop%LinearThermalExpansion,6,'LinearThermalExpansion','[K^(-1)] (alpha) Linear thermal expansion matrix',ierr)

      Call PetscBagRegisterEnum(bag,matprop%HookesLaw%type,MEF90HookesLawTypeList,default%HookesLaw%type,'hookeslaw_type','Type of Hooke''s law',ierr);CHKERRQ(ierr)
      Select case(matprop%HookesLaw%type)
         Case (MEF90HookesLawTypeFull)
            matprop%HookesLaw%fullTensor = default%HookesLaw%fullTensor
            Call PetscBagRegisterRealArray(bag,matprop%HookesLaw%fullTensor,21,'HookesLaw','[N.m^(-2)] (A) Hooke''s law',ierr)
         Case(MEF90HookesLawTypeIsotropic)
            Call PetscBagRegisterReal(bag,matprop%HookesLaw%YoungsModulus,default%HookesLaw%YoungsModulus,'hookeslaw_YoungsModulus','[N.m^(-2)] (E) Young''s Modulus',ierr)
            Call PetscBagRegisterReal(bag,matprop%HookesLaw%PoissonRatio,default%HookesLaw%PoissonRatio,'hookeslaw_PoissonRatio','[] (nu) Poisson Modulus',ierr)
            matprop%HookesLaw%fulltensor = -1.D+30
      End Select
      Call PetscBagRegisterReal(bag,matprop%internalLength,default%internalLength,'internalLength','[m] (l) Internal Length',ierr)
      Call PetscBagRegisterReal(bag,matprop%CoefficientLinSoft,default%CoefficientLinSoft,'CoefficientLinSoft','[] (k) Linear softening coefficient for LinSoft',ierr)

      Call PetscBagRegisterReal(bag,matprop%yieldStress,default%yieldStress,'yieldStress','[N.m^(-2)] (sigma_y) stress threshold for plasticity',ierr)
      Call PetscBagRegisterReal(bag,matprop%residualYieldStress,default%residualYieldStress,'residualyieldStress','[unit-less] percentage of the yield stress',ierr)
      Call PetscBagRegisterReal(bag,matprop%DuctileCouplingPower,default%DuctileCouplingPower,'DuctileCouplingPower','[] power of the coupling between the damage and the plasticity',ierr)
      Call PetscBagRegisterReal(bag,matprop%CoefficientCapModel0,default%CoefficientCapModel0,'CoefficientCapModel0','C0 in the Yield function: CD || dev(stress) || + C2 tr(stress)^2 + C1 tr(stress) - C0 <= 0',ierr)
      Call PetscBagRegisterReal(bag,matprop%CoefficientCapModel1,default%CoefficientCapModel1,'CoefficientCapModel1','C1 in the Yield function: CD || dev(stress) || + C2 tr(stress)^2 + C1 tr(stress) - C0 <= 0',ierr)
      Call PetscBagRegisterReal(bag,matprop%CoefficientCapModel2,default%CoefficientCapModel2,'CoefficientCapModel2','C2 in the Yield function: CD || dev(stress) || + C2 tr(stress)^2 + C1 tr(stress) - C0 <= 0',ierr)
      Call PetscBagRegisterReal(bag,matprop%CoefficientCapModelD,default%CoefficientCapModelD,'CoefficientCapModelD','CD in the Yield function: CD || dev(stress) || + C2 tr(stress)^2 + C1 tr(stress) - C0 <= 0',ierr)
      Call PetscBagRegisterReal(bag,matprop%CoefficientDruckerPrager,default%CoefficientDruckerPrager,'CoefficientDruckerPrager','k in the Yield function: || dev(stress) || - k tr(stress) - yieldStress <= 0',ierr)
      Call PetscBagRegisterReal(bag,matprop%CoeffF,default%CoeffF,'CoeffF','[unit-less] (F) coefficient F in the Hill yield criterion',ierr)
      Call PetscBagRegisterReal(bag,matprop%CoeffG,default%CoeffG,'CoeffG','[unit-less] (G) coefficient G in the Hill yield criterion',ierr)
      Call PetscBagRegisterReal(bag,matprop%CoeffH,default%CoeffH,'CoeffH','[unit-less] (H) coefficient H in the Hill yield criterion',ierr)
      Call PetscBagRegisterReal(bag,matprop%CoeffM,default%CoeffM,'CoeffM','[unit-less] (M) coefficient M in the Hill yield criterion',ierr)
      Call PetscBagRegisterReal(bag,matprop%CoeffN,default%CoeffN,'CoeffN','[unit-less] (N) coefficient N in the Hill yield criterion',ierr)
      Call PetscBagRegisterReal(bag,matprop%CoeffL,default%CoeffL,'CoeffL','[unit-less] (L) coefficient L in the Hill yield criterion',ierr)
      Call PetscBagRegisterReal(bag,matprop%YieldTau0,default%YieldTau0,'YieldTau0','[N.m^(-2)] (tau_0) stress threshold in the Hill yield criterion',ierr)
      Call PetscBagRegisterReal(bag,matprop%residualYieldTau0,default%residualYieldTau0,'residualYieldTau0','[unit-less] residual stress threshold in the Hill yield criterion',ierr)
      Call PetscBagRegisterReal(bag,matprop%phi1,default%phi1,'phi1','[radians] Bunge-Euler angle in the Hill yield criterion',ierr)
      Call PetscBagRegisterReal(bag,matprop%phi2,default%phi2,'phi2','[radians] Bunge-Euler angle in the Hill yield criterion',ierr)
      Call PetscBagRegisterReal(bag,matprop%Phi,default%Phi,'Phi','[radians] Bunge-Euler angle in the Hill yield criterion',ierr)
      Call PetscBagRegisterReal(bag,matprop%delta,default%delta,'delta','[unit-less] residual in the definition of the porosity, Gurson and Green criteria',ierr)

      Call PetscBagRegisterReal(bag,matprop%cohesiveStiffness,default%cohesiveStiffness,'cohesiveStiffness','[N.m^(-4)] (k) cohesive stiffness in Winkler-type models',ierr)
      Call PetscBagRegisterReal(bag,matprop%residualStiffness,default%residualStiffness,'residualStiffness','[unit-less] (eta) residual stiffness',ierr)

      Call PetscBagRegisterReal(bag,matprop%drivingForceTensileStrength,default%drivingForceTensileStrength,'drivingForce_tensileStrength','[N.m^(-2)] (\sigma_{ts}) tensile strength in Drucker-Prager driving Force',ierr)
      Call PetscBagRegisterReal(bag,matprop%drivingForceCompressiveStrength,default%drivingForceCompressiveStrength,'drivingForce_CompressiveStrength','[N.m^(-2)] (\sigma_{cs}) compressive strength in Drucker-Prager driving Force',ierr)
      Call PetscBagRegisterReal(bag,matprop%drivingForceBeta,default%drivingForceBeta,'drivingForce_Beta','[unit-less] (\beta) alpha parameter in Drucker-Prager driving Force',ierr)
      Call PetscBagRegisterInt(bag,matprop%drivingForcep,default%drivingForcep,'drivingForce_p','[unit-less] (p) p parameter in Drucker-Prager driving Force',ierr)

      Call PetscBagRegisterBool(bag,matprop%isLinearIsotropicHardening,default%isLinearIsotropicHardening,'isLinearIsotropicHardening','[bool] Plasticity with Linear Isotropic Hardening',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterBool(bag,matprop%isNoPlCoupling,default%isNoPlCoupling,'isNoPlCoupling','[bool] Coupling between damage and plastic dissipation',ierr);CHKERRQ(ierr)
      !Call PetscBagSetFromOptions(bag,ierr)
   End Subroutine PetscBagRegisterMEF90MatProp3D


#undef __FUNCT__
#define __FUNCT__ "MEF90MatPropBagSetFromOptions2D"
!!!
!!!
!!!  MEF90MatPropBagSetFromOptionsierr:
!!!  MEF90MatPropBagSetFromOptions2D
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90MatPropBagSetFromOptions2D(MEF90MatPropBag,Mesh,defaultMaterial,MEF90Ctx,ierr)
      PetscBag,Dimension(:),Pointer                   :: MEF90MatPropBag
      Type(DM),Intent(IN)                             :: Mesh
      Type(MEF90MatProp2D_Type),intent(IN)            :: defaultMaterial
      PetscErrorCode,Intent(OUT)                      :: ierr
      Type(MEF90Ctx_Type),Intent(IN)                  :: MEF90Ctx

      Type(IS)                                        :: setIS
      PetscInt,Dimension(:),Pointer                   :: setID
      PetscInt                                        :: numSet,set
      Character(len=MEF90_MXSTRLEN)                   :: setName,setprefix,IOBuffer
      Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90GlobalOptions

      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)
      Call DMmeshGetLabelIdIS(Mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr)
      Call ISGetLocalSize(setIS,numSet,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Allocate(MEF90MatPropBag(numSet))
      Do set = 1,numSet
         Write(setName,100) setID(set)
         Write(setprefix,101) setID(set)
         If (MEF90GlobalOptions%verbose > 0) Then
            Write(IOBuffer,102) setID(set),trim(setprefix)
            Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
         End If
         Call PetscBagCreate(PETSC_COMM_WORLD,sizeofMEF90MatProp2D,MEF90MatPropBag(set),ierr)
         Call PetscBagRegisterMEF90MatProp(MEF90MatPropBag(set),setName,setprefix,defaultMaterial,ierr);CHKERRQ(ierr)
         If (MEF90GlobalOptions%verbose > 0) Then
            Call PetscBagView(MEF90MatPropBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
            Call PetscPrintf(PETSC_COMM_WORLD,"\n",ierr);CHKERRQ(ierr)
         End If
      End Do
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
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
   Subroutine MEF90MatPropBagSetFromOptions3D(MEF90MatPropBag,Mesh,defaultMaterial,MEF90Ctx,ierr)
      PetscBag,Dimension(:),Pointer                   :: MEF90MatPropBag
      Type(DM),Intent(IN)                             :: Mesh
      Type(MEF90MatProp3D_Type),intent(IN)            :: defaultMaterial
      PetscErrorCode,Intent(OUT)                      :: ierr
      Type(MEF90Ctx_Type),Intent(IN)                  :: MEF90Ctx

      Type(IS)                                        :: setIS
      PetscInt,Dimension(:),Pointer                   :: setID
      PetscInt                                        :: numSet,set
      Character(len=MEF90_MXSTRLEN)                   :: setName,setprefix,IOBuffer
      Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90GlobalOptions

      Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)
      Call DMmeshGetLabelIdIS(Mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr)
      Call ISGetLocalSize(setIS,numSet,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Allocate(MEF90MatPropBag(numSet))
      Do set = 1,numSet
         Write(setName,100) setID(set)
         Write(setprefix,101) setID(set)
         If (MEF90GlobalOptions%verbose > 0) Then
            Write(IOBuffer,102) setID(set),trim(setprefix)
            Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
         End If
         Call PetscBagCreate(PETSC_COMM_WORLD,sizeofMEF90MatProp3D,MEF90MatPropBag(set),ierr)
         Call PetscBagRegisterMEF90MatProp(MEF90MatPropBag(set),setName,setprefix,defaultMaterial,ierr);CHKERRQ(ierr)
         If (MEF90GlobalOptions%verbose > 0) Then
            Call PetscBagView(MEF90MatPropBag(set),PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
            Call PetscPrintf(PETSC_COMM_WORLD,"\n",ierr);CHKERRQ(ierr)
         End If
      End Do
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
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
            Call PetscLogFlops(10._pflop,ierr);CHKERRQ(ierr)
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

            Call PetscLogFlops(21._pflop,ierr);CHKERRQ(ierr)
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
            Call PetscLogFlops(11._pflop,ierr);CHKERRQ(ierr)
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
            Call PetscLogFlops(24._pflop,ierr);CHKERRQ(ierr)
         Case(MEF90HookesLawTypeFull)
            MEF90HookesLaw3DXMat3D = A%fullTensor * X
            ! flops are counted in m_MEF90_LinAlg
      End Select
   End Function MEF90HookesLaw3DXMat3D

#undef __FUNCT__
#define __FUNCT__ "MasonryProjection2D"
!!!
!!!
!!!  MasonryProjection2D:
!!!
!!!  (c) 2016 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MasonryProjection2D(Epsilon,A,PositivePart,NegativePart)
      Type(MatS2D),Intent(IN)                     :: Epsilon
      Type(MEF90HookesLaw2D),Intent(IN)           :: A
      Type(MatS2D),Intent(OUT)                    :: PositivePart,NegativePart

      Type(MatS2D)                                :: D
      Type(Mat2D)                                 :: Pinv
      PetscErrorCode                              :: ierr

      If (A%type /= MEF90HookesLawTypeIsotropic) Then
         SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Masonry projection not implemented for non isotropic Hooke laws: "//__FUNCT__,ierr)
      End If

if ((epsilon .dotP. epsilon) < 1.e-8) Then
   PositivePart = 0.0_Kr
   NegativePart = 0.0_Kr
Else
      Call Diagonalize(Epsilon,Pinv,D)
      If (D%XX >= 0.0_Kr) Then
         PositivePart = D
      Else If (A%lambda * D%XX + (A%lambda + 2.0_Kr * A%mu) * D%YY >= 0.0_Kr ) Then
         PositivePart = 0.0_Kr
         PositivePart%YY = A%lambda / (A%lambda + 2.0_Kr * A%mu) * D%XX + D%YY
      Else
         PositivePart = 0.0_Kr
      End If
      PositivePart = MatRaRt(PositivePart,Pinv)
      NegativePart = Epsilon - PositivePart
End If
   End Subroutine MasonryProjection2D

#undef __FUNCT__
#define __FUNCT__ "MasonryProjection3D"
!!!
!!!
!!!  MasonryProjection3D:
!!!
!!!  (c) 2016 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MasonryProjection3D(Epsilon,A,PositivePart,NegativePart)
      Type(MatS3D),Intent(IN)                     :: Epsilon
      Type(MEF90HookesLaw3D),Intent(IN)           :: A
      Type(MatS3D),Intent(OUT)                    :: PositivePart,NegativePart

      Type(MatS3D)                                :: D
      Type(Mat3D)                                 :: Pinv
      PetscReal                                   :: nu
      PetscErrorCode                              :: ierr

      If (A%type /= MEF90HookesLawTypeIsotropic) Then
         SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Masonry projection not implemented for non isotropic Hooke laws: "//__FUNCT__,ierr)
      End If
      Call Diagonalize(Epsilon,Pinv,D)
      nu = A%lambda / (A%lambda + A%mu) * 0.5_Kr
      If (D%XX >= 0.0_Kr) Then
         PositivePart = D
      Else If (nu * D%XX + D%YY >= 0.0_Kr ) Then
         PositivePart = 0.0_Kr
         PositivePart%YY = nu * D%XX + D%YY
         PositivePart%ZZ = nu * D%XX + D%ZZ
      Else If (nu * (D%XX + D%YY) + (1.0_Kr - nu) * D%ZZ >= 0.0_Kr ) Then
         PositivePart = 0.0_Kr
         PositivePart%ZZ = nu / (1.0_Kr - nu) * (D%XX + D%YY) + D%ZZ
      Else
         PositivePart = 0.0_Kr
      End If
      PositivePart = MatRaRt(PositivePart,Pinv)
      NegativePart = Epsilon - PositivePart
   End Subroutine MasonryProjection3D

#undef __FUNCT__
#define __FUNCT__ "HDProjection2D"
!!!
!!!
!!!  HDProjection2D:
!!!
!!!  (c) 2016 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine HDProjection2D(Epsilon,A,PositivePart,NegativePart)
      Type(MatS2D),Intent(IN)                     :: Epsilon
      Type(MEF90HookesLaw2D),Intent(IN)           :: A
      Type(MatS2D),Intent(OUT)                    :: PositivePart,NegativePart

      PetscErrorCode                              :: ierr

      If (A%type /= MEF90HookesLawTypeIsotropic) Then
         SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Hydrostatic-Deviatoric projection not implemented for non isotropic Hooke laws: "//__FUNCT__,ierr)
      End If

      If (trace(Epsilon) >= 0.0_Kr) Then
         PositivePart = Epsilon
         NegativePart = 0.0_Kr
      Else
         PositivePart = deviatoricPart(Epsilon)
         NegativePart = hydrostaticPart(Epsilon)
      End If
   End Subroutine HDProjection2D

#undef __FUNCT__
#define __FUNCT__ "HDProjection3D"
!!!
!!!
!!!  HDProjection3D:
!!!
!!!  (c) 2016 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine HDProjection3D(Epsilon,A,PositivePart,NegativePart)
      Type(MatS3D),Intent(IN)                     :: Epsilon
      Type(MEF90HookesLaw3D),Intent(IN)           :: A
      Type(MatS3D),Intent(OUT)                    :: PositivePart,NegativePart

      PetscReal                                   :: tr
      PetscErrorCode                              :: ierr

      If (A%type /= MEF90HookesLawTypeIsotropic) Then
         SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Hydrostatic-Deviatoric projection not implemented for non isotropic Hooke laws: "//__FUNCT__,ierr)
      End If

      If (trace(Epsilon) >= 0.0_Kr) Then
         PositivePart = Epsilon
         NegativePart = 0.0_Kr
      Else
         PositivePart = deviatoricPart(Epsilon)
         NegativePart = hydrostaticPart(Epsilon)
      End If
   End Subroutine HDProjection3D

End Module m_MEF90_Materials
