Module m_MEF_Materials_Types
#include "finclude/petscdef.h"
   Use m_MEF_LinAlg
   Use petsc
   IMPLICIT NONE

   Type MEF90_MatProp2D_Type
      Sequence
      PetscReal                     :: Density                 ! rho
      PetscReal                     :: FractureToughness       ! Gc
      Type(MatS2D)                  :: ThermalDiffusivity      ! 
      Type(MatS2D)                  :: LinearThermalExpansion  ! 
      Type(Tens4OS2D)               :: HookesLaw         
      Character(len=MEF90_MXSTRLEN) :: Name
   End Type MEF90_MatProp2D_Type

   Type MEF90_MatProp3D_Type
      Sequence
      PetscReal                     :: Density                 ! rho
      PetscReal                     :: FractureToughness       ! Gc
      Type(MatS3D)                  :: ThermalDiffusivity      ! 
      Type(MatS3D)                  :: LinearThermalExpansion  ! 
      Type(Tens4OS3D)               :: HookesLaw         
      Character(len=MEF90_MXSTRLEN) :: Name
   End Type MEF90_MatProp3D_Type
   
   !!! The Mathium is a bogus isotropic material whose material properties are all 1 
   !!! except for a Poisson Ratio of 0
   Type(MEF90_MatProp2D_Type),Parameter     :: MEF90_Mathium2D = MEF90_MatProp2D_Type ( &
      1.0_Kr,                                          & ! Density
      1.0_Kr,                                          & ! FractureToughness
      MEF90_MatS2D_Identity,                           & ! ThermalDiffusivity
      MEF90_MatS2D_Identity,                           & ! LinearThermalExpansion
      MEF90_Tens4OS2D_Identity,                        & ! HookesLaw
      "MEF90_Mathium2D")  

   Type(MEF90_MatProp3D_Type),Parameter     :: MEF90_Mathium3D = MEF90_MatProp3D_Type ( &
      1.0_Kr,                                          & ! Density
      1.0_Kr,                                          & ! FractureToughness
      MEF90_MatS3D_Identity,                           & ! ThermalDiffusivity
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
   Public :: PetscBagGetDataMEF90_MatProp2D

   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF_Materials_Types
         PetscBag                             :: bag
         type(MEF90_MatProp2D_Type),pointer   :: data
         PetscErrorCode                       :: ierr
      End subroutine PetscBagGetData
   End interface
Contains
   Subroutine PetscBagGetDataMEF90_MatProp2D(bag,data,ierr)
      PetscBag                              :: bag
      type(MEF90_MatProp2D_Type),pointer    :: data
      PetscErrorCode                        :: ierr
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataMEF90_MatProp2D
End Module m_MEF_Materials_Interface2D

Module m_MEF_Materials_Interface3D
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use petsc
   Use m_MEF_Materials_Types
   Implicit NONE
   Private
   Public :: PetscBagGetDataMEF90_MatProp3D
   
   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF_Materials_Types
         PetscBag                             :: bag
         type(MEF90_MatProp3D_Type),pointer   :: data
         PetscErrorCode                       :: ierr
      End subroutine PetscBagGetData
   End interface
Contains
   Subroutine PetscBagGetDataMEF90_MatProp3D(bag,data,ierr)
      PetscBag                              :: bag
      type(MEF90_MatProp3D_Type),pointer    :: data
      PetscErrorCode                        :: ierr
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataMEF90_MatProp3D
End Module m_MEF_Materials_Interface3D

Module m_MEF_Materials
   Use m_MEF_Materials_Types
   Use m_MEF_Materials_Interface2D
   Use m_MEF_Materials_Interface3D
   Implicit NONE
   
   Interface PetscBagGetDataMEF90_MatProp
      Module Procedure PetscBagGetDataMEF90_MatProp2D,PetscBagGetDataMEF90_MatProp3D
   End Interface 
   
   Interface PetscBagRegisterMEF90_MatProp
      Module Procedure PetscBagRegisterMEF90_MatProp2D,PetscBagRegisterMEF90_MatProp3D
   End Interface
   
   Interface MEF90_MatPropView
     Module procedure MEF90_MatProp2DView,MEF90_MatProp3DView
   End Interface
   
   PetscSizeT,protected   :: sizeofMEF90_MatProp2D,sizeofMEF90_MatProp3D


Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_MaterialsInitialize"
   Subroutine MEF90_MaterialsInitialize(ierr)
      PetscErrorCode,intent(OUT)          :: ierr

      Type(MEF90_MatProp2D_Type),Target   :: matProp2D
      Type(MEF90_MatProp3D_Type),Target   :: matProp3D
      character(len=1),pointer            :: dummychar(:)
      PetscSizeT                          :: sizeofchar
   
      Call PetscDataTypeGetSize(PETSC_CHAR,sizeofchar,ierr)
      sizeofMEF90_MatProp2D = size(transfer(matProp2D,dummychar))*sizeofchar
      sizeofMEF90_MatProp3D = size(transfer(matProp3D,dummychar))*sizeofchar
   End Subroutine MEF90_MaterialsInitialize

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90_MatProp2D"
   Subroutine PetscBagRegisterMEF90_MatProp2D(bag,name,prefix,default,ierr)
      PetscBag                               :: bag
      Character(len=*),intent(IN)            :: prefix,name
      type(MEF90_MatProp2D_Type),intent(IN)  :: default
      PetscErrorCode,intent(OUT)             :: ierr

      Type(MEF90_MatProp2D_Type),pointer     :: matprop
      
      Call PetscBagGetDataMEF90_MatProp2D(bag,matprop,ierr)
      Call PetscBagSetName(bag,trim(name),"MatProp2D object: material properties",ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix), ierr)
      Call PetscBagRegisterReal(bag,matprop%density,default%density,'Density','density (rho)',ierr)
      Call PetscBagRegisterReal(bag,matprop%FractureToughness,default%FractureToughness,'FractureToughness','Fracture toughness (G_c)',ierr)
      matprop%ThermalDiffusivity = default%ThermalDiffusivity
      Call PetscBagRegisterRealArray(bag,matprop%ThermalDiffusivity,3,'ThermalDiffusivity','ThermalDiffusivity (kappa)',ierr)
      matprop%LinearThermalExpansion = default%LinearThermalExpansion
      Call PetscBagRegisterRealArray(bag,matprop%LinearThermalExpansion,3,'LinearThermalExpansion','Linear Thermal Expansion (alpha)',ierr)
      matprop%HookesLaw = default%HookesLaw
      Call PetscBagRegisterRealArray(bag,matprop%HookesLaw,6,'HookesLaw','Hooke''s law (A)',ierr)
      Call PetscBagRegisterString(bag,matprop%name,trim(default%name),'Name','Material name',ierr)
      !Call PetscBagSetFromOptions(bag,ierr)
   End Subroutine PetscBagRegisterMEF90_MatProp2D

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterMEF90_MatProp3D"
   Subroutine PetscBagRegisterMEF90_MatProp3D(bag,name,prefix,default,ierr)
      PetscBag                               :: bag
      Character(len=*),intent(IN)            :: prefix,name
      type(MEF90_MatProp3D_Type),intent(IN)  :: default
      PetscErrorCode,intent(OUT)             :: ierr

      Type(MEF90_MatProp3D_Type),pointer     :: matprop
      Call PetscBagGetDataMEF90_MatProp3D(bag,matprop,ierr)
      Call PetscBagSetName(bag,trim(name),"MatProp3D object: material properties",ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix), ierr)
      Call PetscBagRegisterReal(bag,matprop%density,default%density,'Density','density (rho)',ierr)
      Call PetscBagRegisterReal(bag,matprop%FractureToughness,default%FractureToughness,'FractureToughness','Fracture toughness (G_c)',ierr)
      matprop%ThermalDiffusivity = default%ThermalDiffusivity
      Call PetscBagRegisterRealArray(bag,matprop%ThermalDiffusivity,6,'ThermalDiffusivity','ThermalDiffusivity (kappa)',ierr)
      matprop%LinearThermalExpansion = default%LinearThermalExpansion
      Call PetscBagRegisterRealArray(bag,matprop%LinearThermalExpansion,6,'LinearThermalExpansion','Linear Thermal Expansion (alpha)',ierr)
      matprop%HookesLaw = default%HookesLaw
      Call PetscBagRegisterRealArray(bag,matprop%HookesLaw,21,'HookesLaw','Hooke''s law (A)',ierr)
      Call PetscBagRegisterString(bag,matprop%name,trim(default%name),'Name','Material name',ierr)
      !Call PetscBagSetFromOptions(bag,ierr)
   End Subroutine PetscBagRegisterMEF90_MatProp3D



#undef __FUNCT__
#define __FUNCT__ "MEF90_MatProp2DView"
   Subroutine MEF90_MatProp2DView(MatProp,viewer,ierr)
      Type(MEF90_MatProp2D_Type),Intent(IN) :: MatProp
      Type(PetscViewer),Intent(IN)          :: viewer
      PetscErrorCode,Intent(OUT)            :: ierr
      
      Character(len=MEF90_MXSTRLEN+100)         :: IOBuffer

      Write(IOBuffer,100) trim(MatProp%Name)
      Call PetscViewerASCIIPrintf(viewer,IOBuffer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,*) "\n"
      Call PetscViewerASCIIPrintf(viewer,IOBuffer,ierr);CHKERRQ(ierr)

      Write(IOBuffer,101) MatProp%Density
      Call PetscViewerASCIIPrintf(IOBuffer,viewer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,102) MatProp%FractureToughness
      Call PetscViewerASCIIPrintf(IOBuffer,viewer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,103) MatProp%ThermalDiffusivity
      Call PetscViewerASCIIPrintf(IOBuffer,viewer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,104) MatProp%LinearThermalExpansion
      Call PetscViewerASCIIPrintf(IOBuffer,viewer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,105) MatProp%HookesLaw
      Call PetscViewerASCIIPrintf(IOBuffer,viewer,ierr);CHKERRQ(ierr)
      
100 Format("Name:...................",A)
101 Format("Density:................",ES12.5)
102 Format("FractureToughness:......",ES12.5)
103 Format("ThermalDiffusivity:.....",3(ES12.5,'  '))
104 Format("LinearThermalExpansion:.",3(ES12.5,'  '))
105 Format("HookesLaw:..............",6(ES12.5,'  '))
   End Subroutine MEF90_MatProp2DView  
   
#undef __FUNCT__
#define __FUNCT_ "MEF90_MatProp3DView"
   Subroutine MEF90_MatProp3DView(MatProp,viewer,ierr)
      Type(MEF90_MatProp3D_Type),Intent(IN) :: MatProp
      Type(PetscViewer),Intent(IN)          :: viewer
      PetscErrorCode,Intent(OUT)            :: ierr
      
      Character(len=MEF90_MXSTRLEN)         :: IOBuffer
      Write(IOBuffer,100) Trim(MatProp%Name)
      Call PetscViewerASCIIPrintf(IOBuffer,viewer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,101) MatProp%Density
      Call PetscViewerASCIIPrintf(IOBuffer,viewer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,102) MatProp%FractureToughness
      Call PetscViewerASCIIPrintf(IOBuffer,viewer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,103) MatProp%ThermalDiffusivity
      Call PetscViewerASCIIPrintf(IOBuffer,viewer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,104) MatProp%LinearThermalExpansion
      Call PetscViewerASCIIPrintf(IOBuffer,viewer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,105) MatProp%HookesLaw
      Call PetscViewerASCIIPrintf(IOBuffer,viewer,ierr);CHKERRQ(ierr)
      
100 Format("Name:...................",A)
101 Format("Density:................",ES12.5)
102 Format("FractureToughness:......",ES12.5)
103 Format("ThermalDiffusivity:.....",6(ES12.5,'  '))
104 Format("LinearThermalExpansion:.",6(ES12.5,'  '))
105 Format("HookesLaw:..............",21(ES12.5,'  '))
   End Subroutine MEF90_MatProp3DView  
   

!!! Subroutine generating various types of Hooke's laws 
#undef __FUNCT__
#define __FUNCT__ "MEF90_HookeLawIsoLambdaMu_2D"
   Subroutine MEF90_HookeLawIsoLambdaMu_2D(A,lambda,mu) 
      Type(Tens4OS2D),Intent(OUT)        :: A
      PetscReal,Intent(IN)               :: lambda,mu
      A = 0.0_Kr
      
      A%XXXX = lambda + 2.0_Kr * mu
      A%XXYY = lambda
      A%XYXY = mu
      A%YYYY = lambda + 2.0_Kr * mu
   End Subroutine MEF90_HookeLawIsoLambdaMu_2D         

#undef __FUNCT__
#define __FUNCT__ "MEF90_HookeLawIsoEnu_2DPlaneStress"
   Subroutine MEF90_HookeLawIsoEnu_2DPlaneStress(A,E,nu) 
      PetscReal,Intent(IN)               :: E,nu
      Type(Tens4OS2D),Intent(OUT)        :: A
      
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
      PetscReal,Intent(IN)               :: E,nu
      Type(Tens4OS2D),Intent(OUT)        :: A
      
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
      PetscReal,Intent(IN)               :: lambda,mu
      Type(Tens4OS3D),Intent(OUT)        :: A
   
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
      PetscReal,Intent(IN)               :: E,nu
      Type(Tens4OS3D),Intent(OUT)        :: A
      
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