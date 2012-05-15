Module m_MEF_Materials_Types
#include "finclude/petscdef.h"
   Use m_MEF_LinAlg
   Use petsc
   IMPLICIT NONE

   Type MEF90_MatProp2D_Type
      Sequence
      PetscReal                     :: Density                 ! rho
      PetscReal                     :: FractureToughness       ! Gc
      PetscReal                     :: LatentHeat              ! Cp
      Type(MatS2D)                  :: HeatConductivity        ! Kappa
      Type(MatS2D)                  :: LinearThermalExpansion  ! alpha
      Type(Tens4OS2D)               :: HookesLaw         
      Character(len=MEF90_MXSTRLEN) :: Name
   End Type MEF90_MatProp2D_Type

   Type MEF90_MatProp3D_Type
      Sequence
      PetscReal                     :: Density                 ! rho
      PetscReal                     :: FractureToughness       ! Gc
      PetscReal                     :: LatentHeat              ! Cp
      Type(MatS3D)                  :: HeatConductivity        ! Kappa
      Type(MatS3D)                  :: LinearThermalExpansion  ! alpha
      Type(Tens4OS3D)               :: HookesLaw         
      Character(len=MEF90_MXSTRLEN) :: Name
   End Type MEF90_MatProp3D_Type
   
   !!! The Mathium is a bogus isotropic material whose material properties are all 1 
   !!! except for a Poisson Ratio of 0
   Type(MEF90_MatProp2D_Type),Parameter     :: MEF90_Mathium2D = MEF90_MatProp2D_Type ( &
      1.0_Kr,                                          & ! Density
      1.0_Kr,                                          & ! FractureToughness
      1.0_Kr,                                          & ! LatentHeat
      MEF90_MatS2D_Identity,                           & ! HeatConductivity
      MEF90_MatS2D_Identity,                           & ! LinearThermalExpansion
      MEF90_Tens4OS2D_Identity,                        & ! HookesLaw
      "MEF90_Mathium2D")  

   Type(MEF90_MatProp3D_Type),Parameter     :: MEF90_Mathium3D = MEF90_MatProp3D_Type ( &
      1.0_Kr,                                          & ! Density
      1.0_Kr,                                          & ! FractureToughness
      1.0_Kr,                                          & ! LatentHeat
      MEF90_MatS3D_Identity,                           & ! HeatConductivity
      MEF90_MatS3D_Identity,                           & ! LinearThermalExpansion
      MEF90_Tens4OS3D_Identity,                        & ! HookesLaw
      "MEF90_Mathium3D")  
End Module m_MEF_Materials_Types

Module m_MEF_Materials_Interface2D
#include "finclude/petscdef.h"
#include <finclude/petscbagdef.h>
   Use petsc
   Use m_MEF_Materials_Types
   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF_Materials_Types
         PetscBag                             :: bag
         type(MEF90_MatProp2D_Type),pointer   :: data
         PetscErrorCode                       :: ierr
      End subroutine PetscBagGetData
   End interface
Contains
   Subroutine PetscBagGetData_MEF90_MatProp2D_Type(bag,data,ierr)
      PetscBag                              :: bag
      type(MEF90_MatProp2D_Type),pointer    :: data
      PetscErrorCode                        :: ierr
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetData_MEF90_MatProp2D_Type
End Module m_MEF_Materials_Interface2D

Module m_MEF_Materials_Interface3D
#include "finclude/petscdef.h"
#include <finclude/petscbagdef.h>
   Use petsc
   Use m_MEF_Materials_Types
   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_MEF_Materials_Types
         PetscBag                             :: bag
         type(MEF90_MatProp3D_Type),pointer   :: data
         PetscErrorCode                       :: ierr
      End subroutine PetscBagGetData
   End interface
Contains
   Subroutine PetscBagGetData_MEF90_MatProp3D_Type(bag,data,ierr)
      PetscBag                              :: bag
      type(MEF90_MatProp3D_Type),pointer    :: data
      PetscErrorCode                        :: ierr
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetData_MEF90_MatProp3D_Type
End Module m_MEF_Materials_Interface3D

Module m_MEF_Materials
   Use m_MEF_Materials_Types
   Use m_MEF_Materials_Interface2D
   Use m_MEF_Materials_Interface3D
   
   Interface PetscBagGetData_MEF90_MatProp
      Module Procedure PetscBagGetData_MEF90_MatProp2D_Type,PetscBagGetData_MEF90_MatProp3D_Type
   End Interface 
   
   Interface MEF90_MatPropBagRegister
      Module Procedure MEF90_MatProp2DBagRegister,MEF90_MatProp3DBagRegister
   End Interface
   
   Interface MEF90_MatPropView
     Module procedure MEF90_MatPropView_2D,MEF90_MatPropView_3D 
   End Interface
Contains

#undef __FUNCT__
#define __FUNCT__ "MEF90_MatProp2DBagRegister"
   Subroutine MEF90_MatProp2DBagRegister(bag,name,prefix,default,ierr)
      PetscBag                               :: bag
      Character(len=*),intent(IN)            :: prefix,name
      type(MEF90_MatProp2D_Type),intent(IN)  :: default
      PetscErrorCode,intent(OUT)             :: ierr

      Type(MEF90_MatProp2D_Type),pointer     :: matprop
      Call PetscBagGetData_MEF90_MatProp(bag,matprop,ierr)
      Call PetscBagSetName(bag,trim(name),'Material name',ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix), ierr)
      Call PetscBagRegisterInt(bag,matprop%density,default%density,'density','density (rho)',ierr)
      Call PetscBagRegisterReal(bag,matprop%FractureToughness,default%FractureToughness,'Fracture toughness','Fracture toughness (G_c)',ierr)
      Call PetscBagRegisterReal(bag,matprop%LatentHeat,default%LatentHeat,'Latent Heat','Latent Heat (C_p)',ierr)
      Call PetscBagRegisterRealArray(bag,matprop%HeatConductivity,2,'Heat Conductivityn','HeatConductivity (kappa)',ierr)
      Call PetscBagRegisterRealArray(bag,matprop%LinearThermalExpansion,2,'Linear Thermal Expansion','Linear Thermal Expansion (alpha)',ierr)
      Call PetscBagRegisterRealArray(bag,matprop%HookesLaw,6,'Hooke''s law','Hooke''s law (A)',ierr)
      Call PetscBagSetFromOptions(bag,ierr)
   End Subroutine MEF90_MatProp2DBagRegister

#undef __FUNCT__
#define __FUNCT__ "MEF90_MatProp3DBagRegister"
   Subroutine MEF90_MatProp3DBagRegister(bag,name,prefix,default,ierr)
      PetscBag                               :: bag
      Character(len=*),intent(IN)            :: prefix,name
      type(MEF90_MatProp3D_Type),intent(IN)  :: default
      PetscErrorCode,intent(OUT)             :: ierr

      Type(MEF90_MatProp3D_Type),pointer     :: matprop
      Call PetscBagGetData_MEF90_MatProp(bag,matprop,ierr)
      Call PetscBagSetName(bag,trim(name),'Material name',ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix), ierr)
      Call PetscBagRegisterInt(bag,matprop%density,default%density,'density','density (rho)',ierr)
      Call PetscBagRegisterReal(bag,matprop%FractureToughness,default%FractureToughness,'Fracture toughness','Fracture toughness (G_c)',ierr)
      Call PetscBagRegisterReal(bag,matprop%LatentHeat,default%LatentHeat,'Latent Heat','Latent Heat (C_p)',ierr)
      Call PetscBagRegisterRealArray(bag,matprop%HeatConductivity,6,'Heat Conductivityn','HeatConductivity (kappa)',ierr)
      Call PetscBagRegisterRealArray(bag,matprop%LinearThermalExpansion,6,'Linear Thermal Expansion','Linear Thermal Expansion (alpha)',ierr)
      Call PetscBagRegisterRealArray(bag,matprop%HookesLaw,21,'Hooke''s law','Hooke''s law (A)',ierr)
      Call PetscBagSetFromOptions(bag,ierr)
   End Subroutine MEF90_MatProp3DBagRegister



#undef __FUNCT__
#define __FUNCT__ "MEF90_MatPropView_2D"
   Subroutine MEF90_MatPropView_2D(MatProp,viewer,ierr)
      Type(MEF90_MatProp2D_Type),Intent(IN) :: MatProp
      Type(PetscViewer),Intent(IN)          :: viewer
      PetscErrorCode,Intent(OUT)            :: ierr
      
      Character(len=MEF90_MXSTRLEN)         :: IOBuffer
      Write(IOBuffer,100) Trim(MatProp%Name)
      Call PetscViewerASCIIPrintf(IOBuffer,viewer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,101) MatProp%Density
      Call PetscViewerASCIIPrintf(IOBuffer,viewer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,101) MatProp%FractureToughness
      Call PetscViewerASCIIPrintf(IOBuffer,viewer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,101) MatProp%LatentHeat
      Call PetscViewerASCIIPrintf(IOBuffer,viewer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,101) MatProp%HeatConductivity
      Call PetscViewerASCIIPrintf(IOBuffer,viewer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,101) MatProp%LinearThermalExpansion
      Call PetscViewerASCIIPrintf(IOBuffer,viewer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,101) MatProp%HookesLaw
      Call PetscViewerASCIIPrintf(IOBuffer,viewer,ierr);CHKERRQ(ierr)
      
100 Format("Name:...................",A)
101 Format("Density:................",ES12.5)
102 Format("FractureToughness:......",ES12.5)
103 Format("LatentHeat:.............",ES12.5)
104 Format("HeatConductivity:.......",3(ES12.5,'  '))
105 Format("LinearThermatExpansion:.",3(ES12.5,'  '))
106 Format("HookesLaw:..............",6(ES12.5,'  '))
   End Subroutine MEF90_MatPropView_2D  
   
#undef __FUNCT__
#define __FUNCT__ "MEF90_MatPropView_3D"
   Subroutine MEF90_MatPropView_3D(MatProp,viewer,ierr)
      Type(MEF90_MatProp3D_Type),Intent(IN) :: MatProp
      Type(PetscViewer),Intent(IN)          :: viewer
      PetscErrorCode,Intent(OUT)            :: ierr
      
      Character(len=MEF90_MXSTRLEN)         :: IOBuffer
      Write(IOBuffer,100) Trim(MatProp%Name)
      Call PetscViewerASCIIPrintf(IOBuffer,viewer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,101) MatProp%Density
      Call PetscViewerASCIIPrintf(IOBuffer,viewer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,101) MatProp%FractureToughness
      Call PetscViewerASCIIPrintf(IOBuffer,viewer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,101) MatProp%LatentHeat
      Call PetscViewerASCIIPrintf(IOBuffer,viewer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,101) MatProp%HeatConductivity
      Call PetscViewerASCIIPrintf(IOBuffer,viewer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,101) MatProp%LinearThermalExpansion
      Call PetscViewerASCIIPrintf(IOBuffer,viewer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,101) MatProp%HookesLaw
      Call PetscViewerASCIIPrintf(IOBuffer,viewer,ierr);CHKERRQ(ierr)
      
100 Format("Name:...................",A)
101 Format("Density:................",ES12.5)
102 Format("FractureToughness:......",ES12.5)
103 Format("LatentHeat:.............",ES12.5)
104 Format("HeatConductivity:.......",6(ES12.5,'  '))
105 Format("LinearThermatExpansion:.",6(ES12.5,'  '))
106 Format("HookesLaw:..............",21(ES12.5,'  '))
   End Subroutine MEF90_MatPropView_3D  
   
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