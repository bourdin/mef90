Module m_MEF_Materials_Types
#include "finclude/petscdef.h"
   Use m_MEF_LinAlg
   Use m_MEF_Ctx
   Use petsc
   IMPLICIT NONE

   Type MEF90MatProp2D_Type
      Sequence
      PetscReal                     :: Density                    ! rho
      PetscReal                     :: FractureToughness          ! Gc
      PetscReal                     :: SpecificHeat               ! Cp
      Type(MatS2D)                  :: ThermalConductivity        ! K
      Type(MatS2D)                  :: LinearThermalExpansion     ! alpha
      Type(Tens4OS2D)               :: HookesLaw                  ! A
      Character(len=MEF90_MXSTRLEN) :: Name
   End Type MEF90MatProp2D_Type

   Type MEF90MatProp3D_Type
      Sequence
      PetscReal                     :: Density                    ! rho
      PetscReal                     :: FractureToughness          ! Gc
      PetscReal                     :: SpecificHeat               ! Cp
      Type(MatS3D)                  :: ThermalConductivity        ! K
      Type(MatS3D)                  :: LinearThermalExpansion     ! alpha
      Type(Tens4OS3D)               :: HookesLaw                  ! A
      Character(len=MEF90_MXSTRLEN) :: Name
   End Type MEF90MatProp3D_Type
   
   !!! The Mathium is a bogus isotropic material whose material properties are all 1 
   !!! except for a Poisson Ratio of 0.3
   Type(MEF90MatProp2D_Type),Parameter     :: MEF90_Mathium2D = MEF90MatProp2D_Type ( &
      1.0_Kr,                                          & ! Density
      1.0_Kr,                                          & ! FractureToughness
      1.0_Kr,                                          & ! SpecificHeat
      MEF90_MatS2D_Identity,                           & ! ThermalConductivity
      MEF90_MatS2D_Identity,                           & ! LinearThermalExpansion
      Tens4OS2D( 1.0989010989_Kr,     & ! A%XXXX
                   0.0_Kr,            & ! A%XXXY
                   0.32967032967_Kr,  & ! A%XXYY
                   0.384615384615_Kr, & ! A%XYXY
                   0.0_Kr,            & ! A%XYYY
                   1.0989010989_Kr),  & ! A%YYYY        & ! HookesLaw
      !MEF90_Tens4OS2D_Identity,                          & ! HookesLaw
      "MEF90_Mathium2D")  

   Type(MEF90MatProp3D_Type),Parameter     :: MEF90_Mathium3D = MEF90MatProp3D_Type ( &
      1.0_Kr,                                          & ! Density
      1.0_Kr,                                          & ! FractureToughness
      1.0_Kr,                                          & ! SpecificHeat
      MEF90_MatS3D_Identity,                           & ! ThermalConductivity
      MEF90_MatS3D_Identity,                           & ! LinearThermalExpansion
      MEF90_Tens4OS3D_Identity,                        & ! HookesLaw
      "MEF90_Mathium3D")  
End Module m_MEF_Materials_Types

Module m_MEF_Materials_Interface2D
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use petsc
   Use m_MEF_Materials_Types
   Implicit NONE
   Private
   Public :: PetscBagGetDataMEF90MatProp2D

   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF_Materials_Types
         PetscBag                             :: bag
         type(MEF90MatProp2D_Type),pointer    :: data
         PetscErrorCode                       :: ierr
      End subroutine PetscBagGetData
   End interface
Contains
   Subroutine PetscBagGetDataMEF90MatProp2D(bag,data,ierr)
      PetscBag                              :: bag
      type(MEF90MatProp2D_Type),pointer     :: data
      PetscErrorCode                        :: ierr
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataMEF90MatProp2D
End Module m_MEF_Materials_Interface2D

Module m_MEF_Materials_Interface3D
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use petsc
   Use m_MEF_Materials_Types
   Implicit NONE
   Private
   Public :: PetscBagGetDataMEF90MatProp3D
   
   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF_Materials_Types
         PetscBag                             :: bag
         type(MEF90MatProp3D_Type),pointer    :: data
         PetscErrorCode                       :: ierr
      End subroutine PetscBagGetData
   End interface
Contains
   Subroutine PetscBagGetDataMEF90MatProp3D(bag,data,ierr)
      PetscBag                              :: bag
      type(MEF90MatProp3D_Type),pointer     :: data
      PetscErrorCode                        :: ierr
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataMEF90MatProp3D
End Module m_MEF_Materials_Interface3D

Module m_MEF_Materials
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use petsc
   Use m_MEF_Materials_Types
   Use m_MEF_Materials_Interface2D
   Use m_MEF_Materials_Interface3D
   Implicit NONE
   
   Interface PetscBagGetDataMEF90MatProp
      Module Procedure PetscBagGetDataMEF90MatProp2D,PetscBagGetDataMEF90MatProp3D
   End Interface 
   
   Interface PetscBagRegisterMEF90MatProp
      Module Procedure PetscBagRegisterMEF90MatProp2D,PetscBagRegisterMEF90MatProp3D
   End Interface
   
   Interface MEF90MatPropBag_SetFromOptions
      Module Procedure MEF90MatPropBag_SetFromOptions2D,MEF90MatPropBag_SetFromOptions3D
   End Interface
   
   PetscSizeT,protected   :: sizeofMEF90MatProp2D
   PetscSizeT,protected   :: sizeofMEF90MatProp3D


Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90Materials_InitializePrivate"
   Subroutine MEF90Materials_InitializePrivate(ierr)
      PetscErrorCode,intent(OUT)          :: ierr

      Type(MEF90MatProp2D_Type),Target    :: matProp2D
      Type(MEF90MatProp3D_Type),Target    :: matProp3D
      character(len=1),pointer            :: dummychar(:)
      PetscSizeT                          :: sizeofchar
   
      Call PetscDataTypeGetSize(PETSC_CHAR,sizeofchar,ierr)
      sizeofMEF90MatProp2D = size(transfer(matProp2D,dummychar))*sizeofchar
      sizeofMEF90MatProp3D = size(transfer(matProp3D,dummychar))*sizeofchar
   End Subroutine MEF90Materials_InitializePrivate

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90MatProp2D"
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
      Call PetscBagRegisterReal(bag,matprop%density,default%density,'Density','[kg.m^(-3)] (rho) Density',ierr)
      Call PetscBagRegisterReal(bag,matprop%FractureToughness,default%FractureToughness,'FractureToughness','[N.m^(-1)] (G_c) Fracture toughness',ierr)
      Call PetscBagRegisterReal(bag,matprop%SpecificHeat,default%SpecificHeat,'SpecificHeat','[J.kg^(-1).K^(-1)] (Cp) Specific heat',ierr)
      matprop%ThermalConductivity = default%ThermalConductivity
      Call PetscBagRegisterRealArray(bag,matprop%ThermalConductivity,3,'ThermalConductivity','[J.m^(-1).s^(-1).K^(-1)] (K) Thermal conductivity',ierr)
      matprop%LinearThermalExpansion = default%LinearThermalExpansion
      Call PetscBagRegisterRealArray(bag,matprop%LinearThermalExpansion,3,'LinearThermalExpansion','[K^(-1)] (alpha) Linear thermal expansion matrix',ierr)
      matprop%HookesLaw = default%HookesLaw
      Call PetscBagRegisterRealArray(bag,matprop%HookesLaw,6,'HookesLaw','[N.m^(-2)] (A) Hooke''s law',ierr)
      !Call PetscBagSetFromOptions(bag,ierr)
   End Subroutine PetscBagRegisterMEF90MatProp2D

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90MatProp3D"
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
      Call PetscBagRegisterReal(bag,matprop%SpecificHeat,default%SpecificHeat,'SpecificHeat','[J.kg^(-1).K^(-1)] (Cp) Specific heat',ierr)
      matprop%ThermalConductivity = default%ThermalConductivity
      Call PetscBagRegisterRealArray(bag,matprop%ThermalConductivity,6,'ThermalConductivity','[J.m^(-1).s^(-1).K^(-1)] (K) Thermal conductivity',ierr)
      matprop%LinearThermalExpansion = default%LinearThermalExpansion
      Call PetscBagRegisterRealArray(bag,matprop%LinearThermalExpansion,6,'LinearThermalExpansion','[K^(-1)] (alpha) Linear thermal expansion matrix',ierr)
      matprop%HookesLaw = default%HookesLaw
      Call PetscBagRegisterRealArray(bag,matprop%HookesLaw,21,'HookesLaw','[N.m^(-2)] (A) Hooke''s law',ierr)
      !Call PetscBagSetFromOptions(bag,ierr)
   End Subroutine PetscBagRegisterMEF90MatProp3D


#undef __FUNCT__
#define __FUNCT__ "MEF90MatPropBag_SetFromOptions2D"
!!!
!!!  
!!!  MEF90MatPropBag_SetFromOptions2D:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90MatPropBag_SetFromOptions2D(MEF90MatPropBag,Mesh,defaultMaterial,MEF90Ctx,ierr)
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
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr) 
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
   End Subroutine MEF90MatPropBag_SetFromOptions2D

#undef __FUNCT__
#define __FUNCT__ "MEF90MatPropBag_SetFromOptions3D"
!!!
!!!  
!!!  MEF90MatPropBag_SetFromOptions3D:
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90MatPropBag_SetFromOptions3D(MEF90MatPropBag,Mesh,defaultMaterial,MEF90Ctx,ierr)
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
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS,ierr);CHKERRQ(ierr) 
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
   End Subroutine MEF90MatPropBag_SetFromOptions3D


!!! Subroutine generating various types of Hooke's laws 
#undef __FUNCT__
#define __FUNCT__ "MEF90_HookeLawIsoLambdaMu_2D"
   Subroutine MEF90_HookeLawIsoLambdaMu_2D(A,lambda,mu) 
      Type(Tens4OS2D),Intent(OUT)         :: A
      PetscReal,Intent(IN)                :: lambda,mu
      A = 0.0_Kr
      
      A%XXXX = lambda + 2.0_Kr * mu
      A%XXYY = lambda
      A%XYXY = mu
      A%YYYY = lambda + 2.0_Kr * mu
   End Subroutine MEF90_HookeLawIsoLambdaMu_2D         

#undef __FUNCT__
#define __FUNCT__ "MEF90_HookeLawIsoEnu_2DPlaneStress"
   Subroutine MEF90_HookeLawIsoEnu_2DPlaneStress(A,E,nu) 
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
   End Subroutine MEF90_HookeLawIsoEnu_2DPlaneStress         
   
#undef __FUNCT__
#define __FUNCT__ "MEF90_HookeLawIsoEnu_2DPlaneStrain"
   Subroutine MEF90_HookeLawIsoEnu_2DPlaneStrain(A,E,nu) 
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
   End Subroutine MEF90_HookeLawIsoEnu_2DPlaneStrain         

#undef __FUNCT__
#define __FUNCT__ "MEF90_HookeLawIsoLambdaMu_3D"
   Subroutine MEF90_HookeLawIsoLambdaMu_3D(A,lambda,mu)
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
   End Subroutine MEF90_HookeLawIsoLambdaMu_3D

#undef __FUNCT__
#define __FUNCT__ "MEF90_HookeLawIsoENu_3D"
   Subroutine MEF90_HookeLawIsoENu_3D(A,E,nu)
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
   End Subroutine MEF90_HookeLawIsoENu_3D
End Module m_MEF_Materials