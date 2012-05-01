Module m_MEF_Materials
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
   !!! except for a Poisson Ratio of .3
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

Contains
#undef __FUNCT__
#define __FUNCT_ "MEF90_MatPropView_2D"
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
#define __FUNCT__ "MEF90_PetscOptionsGetMatProp_2D"
   Subroutine MEF90_PetscOptionsGetMatProp_2D(prefix,MatProp,prefix,ierr)
      Type(MEF90_MatProp2D_Type),Intent(OUT)         :: MatProp
      Character(len=*),Intent(IN)                    :: prefix 
      PetscErrorCode,Intent(OUT)                     :: ierr

      PetscTruth                                     :: flg

      Call PetscOptionsGetReal(prefix,'-density',MatProp%Density,flag,ierr);CHKERRQ(ierr)
   End Subroutine MEF90_PetscOptionsGetMatProp_2D

#undef __FUNCT__
#define __FUNCT__ "MEF90_MatPropHelp"  
   Subroutine MEF90_MatPropHelp(viewer,ierr)
      Type(PetscViewer),Intent(IN)          :: viewer 
      PetscErrorCode,Intent(OUT)            :: ierr
      
   End Subroutine MEF90_MatPropHelp
   
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