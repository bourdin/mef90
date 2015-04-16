Module m_MEF90_Materials_Types
#include "finclude/petscdef.h"
   Use m_MEF90_LinAlg
   Use m_MEF90_Ctx
   Use petsc
   IMPLICIT NONE

   Type MEF90HookesLaw2D
      Sequence
      Type(Tens4OS2D)    :: fullTensor
      PetscReal          :: lambda,mu
      PetscEnum          :: type
      PetscBool          :: isPlaneStress
   End Type MEF90HookesLaw2D
   
   Type MEF90HookesLaw3D
      Sequence
      Type(Tens4OS3D)    :: fullTensor
      PetscReal          :: lambda,mu
      PetscEnum          :: type
   End Type MEF90HookesLaw3D

   Type MEF90MatProp2D_Type
      Sequence
      PetscReal                     :: Density                    ! rho
      PetscReal                     :: FractureToughness          ! Gc
      PetscReal                     :: SpecificHeat               ! Cp
      Type(MatS2D)                  :: ThermalConductivity        ! K
      Type(MatS2D)                  :: LinearThermalExpansion     ! alpha
      Type(MEF90HookesLaw2D)        :: HookesLaw                  ! A
      PetscReal                     :: internalLength             ! l
      PetscReal                     :: residualStiffness          ! eta
      Character(len=MEF90_MXSTRLEN) :: Name
   End Type MEF90MatProp2D_Type

   Type MEF90MatProp3D_Type
      Sequence
      PetscReal                     :: Density                    ! rho
      PetscReal                     :: FractureToughness          ! Gc
      PetscReal                     :: SpecificHeat               ! Cp
      Type(MatS3D)                  :: ThermalConductivity        ! K
      Type(MatS3D)                  :: LinearThermalExpansion     ! alpha
      Type(MEF90HookesLaw3D)        :: HookesLaw                  ! A
      PetscReal                     :: internalLength             ! l
      PetscReal                     :: residualStiffness          ! eta
      Character(len=MEF90_MXSTRLEN) :: Name
   End Type MEF90MatProp3D_Type

   Enum,bind(c)
      enumerator :: MEF90HookesLawTypeFull = 0, &
                    MEF90HookesLawTypeIsotropic
   End Enum

   !!! The Mathium is a bogus isotropic material whose material properties are all 1 
   !!! except for a Poisson Ratio of 0.3
   Type(MEF90MatProp2D_Type),Parameter     :: MEF90Mathium2D = MEF90MatProp2D_Type( &
      1.0_Kr,                                                                       & ! Density
      1.0_Kr,                                                                       & ! FractureToughness
      1.0_Kr,                                                                       & ! SpecificHeat
      MEF90MatS2DIdentity,                                                          & ! ThermalConductivity
      MEF90MatS2DIdentity,                                                          & ! LinearThermalExpansion
      MEF90HookesLaw2D(                                                             & 
         Tens4OS2D(1.09890_Kr,0.32967_Kr,0.00000_Kr,                                & ! HookesLaw XXXX,XXYY,XXXY
                              1.09890_Kr,0.00000_Kr,                                & !                YYYY,YYXY
                                         0.38462_Kr),                               & !                     XYXY        
         !!! EVIL HACK: we store the lame coefficients but query E,nu from the command line
         !!!            the proper conversion is done in PetscBagGetDataMEF90MatProp2D
         1.0_Kr, .3_Kr,                                                             & ! lambda, mu
         MEF90HookesLawTypeIsotropic,                                               & ! type
         .FALSE.),                                                                  & ! isPlaneStress
      1.0_Kr,                                                                       & ! Internal Length
      1.0D-9,                                                                       & ! Residual Stiffness
      "MEF90Mathium2D")  

   Type(MEF90MatProp3D_Type),Parameter     :: MEF90Mathium3D = MEF90MatProp3D_Type(    &  
      1.0_Kr,                                                                          & ! Density
      1.0_Kr,                                                                          & ! FractureToughness
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
         !!! EVIL HACK: we store the lame coefficients but query E,nu from the command line
         !!!            the proper conversion is done in PetscBagGetDataMEF90MatProp3D
         1.0_Kr, .3_Kr,                                                                & ! lambda, mu
         MEF90HookesLawTypeIsotropic),                                                 & ! type
      1.0_Kr,                                                                          & ! Internal Length
      1.0D-9,                                                                          & ! Residual Stiffness
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
   Subroutine PetscBagGetDataMEF90MatProp2D(bag,data,ierr)
      PetscBag                                :: bag
      type(MEF90MatProp2D_Type),pointer       :: data
      PetscErrorCode                          :: ierr
      PetscReal                               :: E,nu
      
      Call PetscBagGetData(bag,data,ierr)
      Select case(data%HookesLaw%type)
         Case(MEF90HookesLawTypeFull)
            Continue
         Case(MEF90HookesLawTypeIsotropic)
            E  = data%HookesLaw%lambda
            nu = data%HookesLaw%mu
            If (data%HookesLaw%isPlaneStress) Then
               data%HookesLaw%lambda = E * nu / (1.0_Kr - nu**2) 
               data%HookesLaw%mu     = E / (1.0_Kr + nu) * .5_Kr
            Else
               data%HookesLaw%lambda = E * nu / (1.0_Kr + nu) / (1.0_Kr - 2.0_Kr * nu)
               data%HookesLaw%mu     = E / (1.0_Kr + nu) * .5_Kr      
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
   Subroutine PetscBagGetDataMEF90MatProp3D(bag,data,ierr)
      PetscBag                                :: bag
      type(MEF90MatProp3D_Type),pointer       :: data
      PetscErrorCode                          :: ierr
      
      PetscReal                               :: E, nu
      Call PetscBagGetData(bag,data,ierr)
      Select case(data%HookesLaw%type)
         Case(MEF90HookesLawTypeFull)
            Continue
         Case(MEF90HookesLawTypeIsotropic)
            E  = data%HookesLaw%lambda
            nu = data%HookesLaw%mu
            data%HookesLaw%lambda = E * nu / (1.0_Kr + nu) / (1.0_Kr - 2.0_Kr * nu)
            data%HookesLaw%mu     = E / (1.0_Kr + nu) * .5_Kr      
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
      Module Procedure MEF90HookesLaw2DXMatS2D,MEF90HookesLaw3DXMatS3D
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
      Call PetscBagRegisterReal(bag,matprop%SpecificHeat,default%SpecificHeat,'SpecificHeat','[J.kg^(-1).K^(-1)] (Cp) Specific heat',ierr)
      matprop%ThermalConductivity = default%ThermalConductivity
      Call PetscBagRegisterRealArray(bag,matprop%ThermalConductivity,3,'ThermalConductivity','[J.m^(-1).s^(-1).K^(-1)] (K) Thermal conductivity',ierr)
      matprop%LinearThermalExpansion = default%LinearThermalExpansion
      Call PetscBagRegisterRealArray(bag,matprop%LinearThermalExpansion,3,'LinearThermalExpansion','[K^(-1)] (alpha) Linear thermal expansion matrix',ierr)
      
      Call PetscBagRegisterEnum(bag,matprop%HookesLaw%type,MEF90HookesLawTypeList,default%HookesLaw%type,'hookeslaw_type','Type of Hooke''s law',ierr);CHKERRQ(ierr)
      Select case(matprop%HookesLaw%type)
         Case (MEF90HookesLawTypeFull)
            matprop%HookesLaw = default%HookesLaw
            Call PetscBagRegisterRealArray(bag,matprop%HookesLaw,6,'HookesLaw','[N.m^(-2)] (A) Hooke''s law',ierr)
         Case(MEF90HookesLawTypeIsotropic)
            !!! EVIL HACK: we store the lame coefficients but query E,nu from the command line
            !!!            the proper conversion is done in PetscBagGetDataMEF90MatProp2D
            Call PetscBagRegisterReal(bag,matprop%HookesLaw%lambda,default%HookesLaw%lambda,'YoungsModulus','[N.m^(-2)] (E) Young''s Modulus',ierr)
            Call PetscBagRegisterReal(bag,matprop%HookesLaw%mu,default%HookesLaw%mu,'PoissonRatio','[] (nu) Poisson Modulus',ierr)
            Call PetscBagRegisterBool(bag,matprop%HookesLaw%isPlaneStress,default%HookesLaw%isPlaneStress,'planeStress','Use plane stress elasticity',ierr);CHKERRQ(ierr)
            matprop%HookesLaw%fulltensor = -1.D+30
      End Select
      Call PetscBagRegisterReal(bag,matprop%internalLength,default%internalLength,'internalLength','[m] (l) Internal Length',ierr)
      Call PetscBagRegisterReal(bag,matprop%residualStiffness,default%residualStiffness,'residualStiffness','[unit-less] (eta) residual stiffness',ierr)
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
      Call PetscBagRegisterReal(bag,matprop%SpecificHeat,default%SpecificHeat,'SpecificHeat','[J.kg^(-1).K^(-1)] (Cp) Specific heat',ierr)
      matprop%ThermalConductivity = default%ThermalConductivity
      Call PetscBagRegisterRealArray(bag,matprop%ThermalConductivity,6,'ThermalConductivity','[J.m^(-1).s^(-1).K^(-1)] (K) Thermal conductivity',ierr)
      matprop%LinearThermalExpansion = default%LinearThermalExpansion
      Call PetscBagRegisterRealArray(bag,matprop%LinearThermalExpansion,6,'LinearThermalExpansion','[K^(-1)] (alpha) Linear thermal expansion matrix',ierr)
      Select case(matprop%HookesLaw%type)
         Case (MEF90HookesLawTypeFull)
            matprop%HookesLaw = default%HookesLaw
            Call PetscBagRegisterRealArray(bag,matprop%HookesLaw,21,'HookesLaw','[N.m^(-2)] (A) Hooke''s law',ierr)
         Case(MEF90HookesLawTypeIsotropic)
            !!! EVIL HACK: we store the lame coefficients but query E,nu from the command line
            !!!            the proper conversion is done in PetscBagGetDataMEF90MatProp3D
            Call PetscBagRegisterReal(bag,matprop%HookesLaw%lambda,default%HookesLaw%lambda,'YoungsModulus','[N.m^(-2)] (E) Young''s Modulus',ierr)
            Call PetscBagRegisterReal(bag,matprop%HookesLaw%mu,default%HookesLaw%mu,'PoissonRatio','[] (nu) Poisson Modulus',ierr)
            matprop%HookesLaw%fulltensor = -1.D+30
      End Select
      Call PetscBagRegisterReal(bag,matprop%internalLength,default%internalLength,'internalLength','[m] (l) Internal Length',ierr)
      Call PetscBagRegisterReal(bag,matprop%residualStiffness,default%residualStiffness,'residualStiffness','[unit-less] (eta) residual stiffness',ierr)
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
!!!  (c) 2015 Blaise Bourdin bourdin@lsu.edu
!!!
   Function MEF90HookesLaw2DXMatS2D(A,X)
      Type(MEF90HookesLaw2D), Intent(IN)           :: A
      Type(MatS2D), Intent(IN)                     :: X
      Type(MatS2D)                                 :: MEF90HookesLaw2DXMatS2D
      
      
      PetscErrorCode                               :: ierr
      Real(Kind = Kr)                              :: C1, C2
      PetscLogDouble                               :: flops

      Select case(A%type)
         Case(MEF90HookesLawTypeIsotropic)
            C1 = A%lambda + 2.0_Kr * A%mu
            C2 = 2.0_Kr * A%mu
            MEF90HookesLaw2DXMatS2D%XX = C1 * X%XX + A%lambda * X%YY
            MEF90HookesLaw2DXMatS2D%YY = A%lambda * X%XX + C1 * X%YY
            MEF90HookesLaw2DXMatS2D%XY = C2 * X%XY
            flops = 10.0
            Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
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
!!!  (c) 2015 Blaise Bourdin bourdin@lsu.edu
!!!
   Function MEF90HookesLaw3DXMatS3D(A,X)
      Type(MEF90HookesLaw3D), Intent(IN)           :: A
      Type(MatS3D), Intent(IN)                     :: X
      Type(MatS3D)                                 :: MEF90HookesLaw3DXMatS3D
      
      PetscErrorCode                               :: ierr
      Real(Kind = Kr)                              :: C1, C2
      PetscLogDouble                               :: flops

      Select case(A%type)
         Case(MEF90HookesLawTypeIsotropic)
            C1 = A%lambda + 2.0_Kr * A%mu
            C2 = 2.0_Kr * A%mu
            MEF90HookesLaw3DXMatS3D%XX = C1 * X%XX + A%lambda * X%YY + A%lambda * X%ZZ
            MEF90HookesLaw3DXMatS3D%YY = A%lambda * X%XX + C1 * X%YY * A%lambda * X%ZZ
            MEF90HookesLaw3DXMatS3D%YY = A%lambda * X%XX + A%lambda * X%YY * C1 * X%ZZ
            MEF90HookesLaw3DXMatS3D%YZ = C2 * X%YZ
            MEF90HookesLaw3DXMatS3D%XZ = C2 * X%XZ
            MEF90HookesLaw3DXMatS3D%XY = C2 * X%XY
            flops = 21.0
            Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         Case(MEF90HookesLawTypeFull)
            MEF90HookesLaw3DXMatS3D = A%fullTensor * X
            ! flops are counted in m_MEF90_LinAlg
      End Select
   End Function MEF90HookesLaw3DXMatS3D

End Module m_MEF90_Materials