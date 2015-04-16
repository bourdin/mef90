#include "../MEF90/mef90.inc"
Module MEF90_APPEND(m_MEF90_GradDamageImplementation_,MEF90_DIM)D
#include "finclude/petscdef.h"
   Use m_MEF90_LinAlg
   Use m_MEF90_Parameters
   Use m_MEF90_Elements
   Use m_MEF90_Utils
   Use m_MEF90_Materials
   IMPLICIT NONE
   Private

   Public :: MEF90GradDamageDispBilinearFormSet
   Public :: MEF90GradDamageDispOperatorSet
   Public :: MEF90GradDamageDispInelasticStrainRHSSetVertex2
   Public :: MEF90GradDamageDispInelasticStrainRHSSetCell
   Public :: MEF90GradDamageElasticEnergySet   

   Public :: MEF90GradDamageDamageBilinearFormSetAT1Elastic
   Public :: MEF90GradDamageDamageOperatorSetAT1Elastic
   Public :: MEF90GradDamageDamageRHSSetAT1Elastic
   Public :: MEF90GradDamageDamageBilinearFormSetAT1
   Public :: MEF90GradDamageDamageOperatorSetAT1
   Public :: MEF90GradDamageDamageRHSSetAT1
   Public :: MEF90GradDamageSurfaceEnergySetAT1
   Public :: MEF90GradDamageDamageBilinearFormSetAT2Elastic
   Public :: MEF90GradDamageDamageOperatorSetAT2Elastic
   Public :: MEF90GradDamageDamageBilinearFormSetAT2
   Public :: MEF90GradDamageDamageOperatorSetAT2
   Public :: MEF90GradDamageDamageRHSSetAT2
   Public :: MEF90GradDamageSurfaceEnergySetAT2

!!! The Euler-Lagrange equation for the alpha problem for AT1 is
!!!
!!! \int_\Omega 2W(u,T,p) \alpha \beta + \frac{G_c\varepsilon}{2c_V}\nabla \alpha \cdot \nabla \beta \, dx = 
!!! \int_\Omega \left(2W(u,T,p)-\frac{G_c}{4c_V\varepsilon}\right) \beta \, dx = 0 for all \beta,
!!! where u denotes the displacement field, T the temperature field, p the plastic strain, and W is the elastic energy density.
!!!
!!! For AT2, the E-L equation is
!!!
!!! \int_\Omega 2\left(W(u,T,p) + \frac{G_c}{4c_V\varepsilon}\right) \alpha \beta + \frac{G_c\varepsilon}{2c_V}\nabla \alpha \cdot \nabla \beta \, dx = 
!!! \int_\Omega 2W(u,T,p) \beta \, dx = 0 for all \beta,
!!!

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageDispBilinearFormSet"
!!!
!!!  
!!!  MEF90GradDamageDispBilinearFormSet:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90GradDamageDispBilinearFormSet(K,mesh,meshScal,alpha,eta,cellIS,A,elem,elemType,elemScal,elemScalType,ierr)
      Type(Mat),Intent(IN)                               :: K 
      Type(DM),Intent(IN)                                :: mesh,meshScal
      Type(SectionReal),Intent(IN)                       :: alpha
      PetscReal,Intent(IN)                               :: eta
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_HOOKESLAW),Intent(IN)                   :: A
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elemScal
      Type(MEF90Element_Type),Intent(IN)                 :: elemScalType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscReal,Dimension(:,:),Pointer                   :: MatElem
      PetscReal,Dimension(:),Pointer                     :: alphaLoc
      PetscReal                                          :: alphaElem,SalphaElem
      PetscInt                                           :: iDoF1,iDoF2,iGauss
      Type(SectionReal)                                  :: defaultSection
      PetscLogDouble                                     :: flops
     
      Call DMMeshGetSectionReal(mesh,'default',defaultSection,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(alphaLoc(elemScalType%numDof))
         Allocate(MatElem(elemType%numDof,elemType%numDof))
         Do cell = 1,size(cellID)   
            MatElem = 0.0_Kr
            Call SectionRealRestrictClosure(alpha,meshScal,cellID(cell),elemScalType%numDof,alphaLoc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               alphaElem = 0.0_Kr
               Do iDoF1 = 1, elemScalType%numDof
                  alphaElem = alphaElem + alphaLoc(iDoF1) * elemScal(cell)%BF(iDoF1,iGauss)
               End Do ! iDoF1
               SalphaElem = eta + (1.0_Kr - alphaElem)**2
               Do iDoF1 = 1,elemType%numDof
                  Do iDoF2 = 1,elemType%numDof
                     MatElem(iDoF2,iDoF1) = MatElem(iDoF2,iDoF1) + elem(cell)%Gauss_C(iGauss) * SalphaElem * &
                                    ((A * elem(cell)%GradS_BF(iDoF1,iGauss)) .DotP. elem(cell)%GradS_BF(iDoF2,iGauss))
                  End Do
               End Do
            End Do
            Call DMmeshAssembleMatrix(K,mesh,defaultSection,cellID(cell),MatElem,ADD_VALUES,ierr);CHKERRQ(ierr)
         End Do
         flops = (3 * elemType%numDof**2 + 2 * elemType%numDof + 3) * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(alphaLoc)
         DeAllocate(MatElem)
      End If   
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(defaultSection,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageDispBilinearFormSet
   
#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageDispOperatorSet"   
!!!
!!!  
!!!  MEF90GradDamageDispOperatorSet:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90GradDamageDispOperatorSet(G,mesh,meshScal,U,alpha,eta,cellIS,A,elem,elemType,elemScal,elemScalType,ierr)
      Type(SectionReal),Intent(IN)                       :: G
      Type(DM),Intent(IN)                                :: mesh,meshScal
      Type(SectionReal),Intent(IN)                       :: U,alpha
      PetscReal,Intent(IN)                               :: eta
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_HOOKESLAW),Intent(IN)                     :: A
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elem
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elemScal
      Type(MEF90Element_Type),Intent(IN)                 :: elemType,elemScalType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscReal,Dimension(:),Pointer                     :: Gloc
      PetscReal,Dimension(:),Pointer                     :: Uloc
      PetscReal,Dimension(:),Pointer                     :: alphaloc
      Type(MEF90_MATS)                                   :: GradSUelem
      PetscReal                                          :: alphaElem,SalphaElem
      PetscInt                                           :: iDoF1,iGauss
      PetscLogDouble                                     :: flops
           
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(Uloc(elemType%numDof))
         Allocate(Gloc(elemType%numDof))
         Allocate(alphaLoc(elemScalType%numDof))
         Do cell = 1,size(cellID)      
            Gloc = 0.0_Kr
            Call SectionRealRestrictClosure(U,mesh,cellID(cell),elemType%numDof,Uloc,ierr);CHKERRQ(ierr)
            Call SectionRealRestrictClosure(alpha,meshScal,cellID(cell),elemScalType%numDof,alphaLoc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               alphaElem = 0.0_Kr
               Do iDoF1 = 1, elemScalType%numDof
                  alphaElem = alphaElem + alphaLoc(iDoF1) * elemScal(cell)%BF(iDoF1,iGauss)
               End Do ! iDoF1
               SalphaElem = eta + (1.0_Kr - alphaElem)**2
               GradSUelem = 0.0_Kr
               Do iDoF1 = 1,elemType%numDof
                  GradSUelem = GradSUelem + Uloc(iDof1) * elem(cell)%GradS_BF(iDoF1,iGauss)
               End Do
               Do iDoF1 = 1,elemType%numDof
                  Gloc(iDoF1) = Gloc(iDoF1) + elem(cell)%Gauss_C(iGauss) * SalphaElem * &
                                 (( (A * elem(cell)%GradS_BF(iDoF1,iGauss)) .DotP. GradSUelem)) 
               End Do
            End Do
            Call SectionRealUpdateClosure(G,mesh,cellID(cell),Gloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do      
         DeAllocate(Gloc)
         DeAllocate(Uloc)
         DeAllocate(alphaLoc)
         flops = (5 * elemType%numDof + 3) * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      End If
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageDispOperatorSet

#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageDispInelasticStrainRHSSetVertex2"
!!!
!!!  
!!!  MEF90GradDamageDispInelasticStrainRHSSetVertex2: Contribution of a vertex-based inelastic strain e0 K0
!!!                                          to a RHS in a cell set
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90GradDamageDispInelasticStrainRHSSetVertex2(RHS,mesh,meshScal,e0,K0,alpha,eta,cellIS,elem,elemType,elemScal,elemScalType,ierr)
      Type(SectionReal),Intent(IN)                       :: RHS
      Type(DM),Intent(IN)                                :: mesh,meshScal
      Type(SectionReal),Intent(IN)                       :: e0
      Type(MEF90_MATS),Intent(IN)                        :: K0 ! The inelastic strain is e0.K0
      Type(SectionReal),Intent(IN)                       :: alpha
      PetscReal,Intent(IN)                               :: eta
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elem
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elemScal
      Type(MEF90Element_Type),Intent(IN)                 :: elemType,elemScalType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscReal,Dimension(:),Pointer                     :: RHSloc,e0loc,alphaloc
      PetscReal                                          :: e0elem
      PetscReal                                          :: alphaElem,SalphaElem
      PetscInt                                           :: iDoF1,iGauss
      PetscLogDouble                                     :: flops
           
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(e0loc(elemScalType%numDof))
         Allocate(RHSloc(elemType%numDof))
         Allocate(alphaloc(elemScalType%numDof))
         Do cell = 1,size(cellID)      
            RHSloc = 0.0_Kr
            Call SectionRealRestrictClosure(e0,mesh,cellID(cell),elemScalType%numDof,e0loc,ierr);CHKERRQ(ierr)
            Call SectionRealRestrictClosure(alpha,meshScal,cellID(cell),elemScalType%numDof,alphaLoc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               alphaElem = 0.0_Kr
               Do iDoF1 = 1, elemScalType%numDof
                  alphaElem = alphaElem + alphaLoc(iDoF1) * elemScal(cell)%BF(iDoF1,iGauss)
               End Do ! iDoF1
               SalphaElem = eta + (1.0_Kr - alphaElem)**2
               e0elem = 0.0_Kr
               Do iDoF1 = 1,elemScalType%numDof
                  e0elem = e0elem + e0loc(iDof1) * elemScal(cell)%BF(iDoF1,iGauss)
               End Do
               Do iDoF1 = 1,elemType%numDof
                  RHSloc(iDoF1) = RHSloc(iDoF1) - elem(cell)%Gauss_C(iGauss) * SalphaElem * &
                                 e0elem * (elem(cell)%GradS_BF(iDoF1,iGauss) .DotP. K0)
               End Do
            End Do
            Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do
      
         DeAllocate(RHSloc)
         DeAllocate(e0loc)
         DeAllocate(alphaloc)
         flops = (7 * elemType%numDof + 3) * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      End If
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageDispInelasticStrainRHSSetVertex2

#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageDispInelasticStrainRHSSetCell"
!!!
!!!  
!!!  MEF90GradDamageDispInelasticStrainRHSSetCell: Contribution of a cell-based inelastic strain e0
!!!                                       to a RHS in a cell set
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90GradDamageDispInelasticStrainRHSSetCell(RHS,mesh,meshMatS,meshScal,e0,HookesLaw,alpha,eta,cellIS,elem,elemType,elemScal,elemScalType,ierr)
      Type(SectionReal),Intent(IN)                       :: RHS
      Type(DM),Intent(IN)                                :: mesh,meshMatS,meshScal
      Type(SectionReal),Intent(IN)                       :: e0
      Type(MEF90_HOOKESLAW),Intent(IN)                   :: HookesLaw
      Type(SectionReal),Intent(IN)                       :: alpha
      PetscReal,Intent(IN)                               :: eta
      Type(IS),Intent(IN)                                :: cellIS
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elem
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elemScal
      Type(MEF90Element_Type),Intent(IN)                 :: elemType,elemScalType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      Type(MEF90_MATS)                                   :: e0_elem
      PetscReal,Dimension(:),Pointer                     :: e0_loc
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscReal,Dimension(:),Pointer                     :: RHSloc,alphaloc
      PetscReal                                          :: alphaElem,SalphaElem
      PetscInt                                           :: iDoF1,iGauss,i
      PetscLogDouble                                     :: flops
           
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(RHSloc(elemType%numDof))
         Allocate(alphaloc(elemScalType%numDof))
         Do cell = 1,size(cellID)      
            RHSloc = 0.0_Kr
            Call SectionRealRestrict(e0,cellID(cell),e0_loc,ierr);CHKERRQ(ierr)
            Call SectionRealRestrictClosure(alpha,meshScal,cellID(cell),elemScalType%numDof,alphaLoc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               alphaElem = 0.0_Kr
               Do iDoF1 = 1, elemScalType%numDof
                  alphaElem = alphaElem + alphaLoc(iDoF1) * elemScal(cell)%BF(iDoF1,iGauss)
               End Do ! iDoF1
               SalphaElem = eta + (1.0_Kr - alphaElem)**2
               e0_elem = e0_loc
               Do iDoF1 = 1,elemType%numDof
                  RHSloc(iDoF1) = RHSloc(iDoF1) - elem(cell)%Gauss_C(iGauss) * SalphaElem * &
                                 (elem(cell)%GradS_BF(iDoF1,iGauss) .DotP. (HookesLaw * e0_elem))
               End Do
            End Do
            Call SectionRealUpdateClosure(RHS,mesh,cellID(cell),RHSloc,ADD_VALUES,ierr);CHKERRQ(iErr)
            Call SectionRealRestore(e0,cellID(cell),e0_loc,ierr);CHKERRQ(ierr)
         End Do
         DeAllocate(alphaloc)
         DeAllocate(RHSloc)
         flops = (5 * elemType%numDof + 3) * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      End If
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageDispInelasticStrainRHSSetCell

#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageElasticEnergySet"
   !!!
   !!!  
   !!!  MEF90GradDamageElasticEnergySet:  Contribution of a cell set to elastic energy. 
   !!!                        It is assumed that the temperature is interpolated on the FE space while the plastic strain 
   !!!                        is cell-based
   !!!  
   !!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
   !!!
   Subroutine MEF90GradDamageElasticEnergySet(energy,x,alpha,plasticStrain,temperature,mesh,meshScal,cellIS,eta,HookesLaw,ThermalExpansion,elemDisplacement,elemDisplacementType,elemScal,elemScalType,ierr)
      PetscReal,Intent(OUT)                              :: energy
      Type(SectionReal),Intent(IN)                       :: x,alpha,plasticStrain,temperature
      Type(DM),Intent(IN)                                :: mesh,meshScal
      Type(IS),Intent(IN)                                :: cellIS
      PetscReal                                          :: eta
      Type(MEF90_HOOKESLAW),Intent(IN)                   :: HookesLaw
      Type(MEF90_MATS),Intent(IN)                        :: ThermalExpansion
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elemScal
      Type(MEF90Element_Type),Intent(IN)                 :: elemDisplacementType,elemScalType
      PetscErrorCode,Intent(OUT)                         :: ierr

      PetscReal,Dimension(:),Pointer                     :: xloc,alphaLoc,temperatureLoc,plasticStrainLoc
      PetscReal                                          :: alphaElem,SalphaElem,temperatureElem
      Type(MEF90_MATS)                                   :: StrainElem,StressElem,plasticStrainElem
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscInt                                           :: iDoF1,iGauss
      PetscLogDouble                                     :: flops
     
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(xloc(elemDisplacementType%numDof))
         Allocate(temperatureloc(elemScalType%numDof))
         Allocate(alphaLoc(elemScalType%numDof))
         Do cell = 1,size(cellID)   
            Call SectionRealRestrictClosure(x,mesh,cellID(cell),elemDisplacementType%numDof,xloc,ierr);CHKERRQ(ierr)
            If (temperature%v /= 0) Then
               Call SectionRealRestrictClosure(temperature,meshScal,cellID(cell),elemScalType%numDof,temperatureLoc,ierr);CHKERRQ(ierr)
            End If
            If (plasticStrain%v /= 0) Then
               Call SectionRealRestrict(plasticStrain,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
            End If
            If (alpha%v /= 0) Then
               Call SectionRealRestrictClosure(alpha,meshScal,cellID(cell),elemScalType%numDof,alphaLoc,ierr);CHKERRQ(ierr)
            End If
            Do iGauss = 1,size(elemDisplacement(cell)%Gauss_C)
               alphaElem = 0.0_Kr
               If (alpha%v /= 0) Then
                  Do iDoF1 = 1, elemScalType%numDof
                     alphaElem = alphaElem + alphaLoc(iDoF1) * elemScal(cell)%BF(iDoF1,iGauss)
                  End Do ! iDoF1
               End If
               SalphaElem = eta + (1.0_Kr - alphaElem)**2

               strainElem = 0.0_Kr
               Do iDoF1 = 1,elemDisplacementType%numDof
                  strainElem = strainElem + xloc(iDof1) * elemDisplacement(cell)%GradS_BF(iDof1,iGauss)
               End Do

               temperatureElem   = 0.0_Kr
               If (temperature%v /= 0) Then
                  Do iDoF1 = 1,elemScalType%numDof
                     temperatureElem = temperatureElem + temperatureLoc(iDof1) * elemScal(cell)%BF(iDof1,iGauss)
                  End Do
                  strainElem = strainElem - (temperatureElem * thermalExpansion)
               End If

               plasticStrainElem = 0.0_Kr
               If (plasticStrain%v /= 0) Then
                  plasticStrainElem = plasticStrainLoc
                  strainElem = strainElem - plasticStrainElem
               End If
               stressElem = HookesLaw * strainElem
               energy = energy + SalphaElem * (strainElem .dotP. stressElem) * elemDisplacement(cell)%Gauss_C(iGauss) * 0.5_Kr
            End Do ! Gauss
            If (plasticStrain%v /= 0) Then
               Call SectionRealRestore(plasticStrain,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
            End If
         End Do ! cell
         flops = 7 * size(elemDisplacement(1)%Gauss_C) * size(cellID) 
         If (alpha%v /= 0) Then
            flops = flops + 2 * elemScalType%numDof * size(elemScal(1)%Gauss_C) * size(cellID) 
         End If
         If (temperature%v /= 0) Then
            flops = flops + 2 * elemScalType%numDof * size(elemScal(1)%Gauss_C) * size(cellID) 
         End If
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(xloc)
         DeAllocate(temperatureloc)
         DeAllocate(alphaloc)
      End If   
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageElasticEnergySet

!!! ============== Functions for the alpha problem ==============
#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageDamageBilinearFormSetAT1Elastic"
!!!
!!!  
!!!  MEF90GradDamageDamageBilinearFormSetAT1Elastic: the bilinear form for the surface energy term only
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90GradDamageDamageBilinearFormSetAT1Elastic(K,mesh,cellIS,internalLength,FractureToughness,elem,elemType,ierr)
      Type(Mat),Intent(IN)                               :: K 
      Type(DM),Intent(IN)                                :: mesh
      Type(IS),Intent(IN)                                :: cellIS
      PetscReal,Intent(IN)                               :: internalLength,FractureToughness
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      Type(SectionReal)                                  :: defaultSection
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscInt                                           :: iDoF1,iDoF2,iGauss
      PetscReal,Dimension(:,:),Pointer                   :: MatElem
      PetscReal                                          :: C2
      PetscLogDouble                                     :: flops

      C2 = FractureToughness * internalLength * .75_Kr

      Call DMMeshGetSectionReal(mesh,'default',defaultSection,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(MatElem(elemType%numDof,elemType%numDof))
         Do cell = 1,size(cellID) 
            MatElem = 0.0_Kr  
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               Do iDoF1 = 1,elemType%numDof
                  Do iDoF2 = 1,elemType%numDof
                     MatElem(iDoF2,iDoF1) = MatElem(iDoF2,iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                    C2 * (elem(cell)%Grad_BF(iDof1,iGauss) .dotP. elem(cell)%Grad_BF(iDof2,iGauss))
                  End Do
               End Do
            End Do
            Call DMmeshAssembleMatrix(K,mesh,defaultSection,cellID(cell),MatElem,ADD_VALUES,ierr);CHKERRQ(ierr)
         End Do ! cell
         flops = 3 * elemType%numDof**2 * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(MatElem)
      End If   
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)      
      Call SectionRealDestroy(defaultSection,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageDamageBilinearFormSetAT1Elastic

#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageDamageOperatorSetAT1Elastic"
!!!
!!!  
!!!  MEF90GradDamageDamageOperatorSetAT1Elastic:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90GradDamageDamageOperatorSetAT1Elastic(G,mesh,alpha,cellIS,internalLength,FractureToughness,elem,elemType,ierr)
      Type(SectionReal),Intent(INOUT)                    :: G
      Type(DM),Intent(IN)                                :: mesh
      Type(SectionReal),Intent(IN)                       :: alpha
      Type(IS),Intent(IN)                                :: cellIS
      PetscReal,Intent(IN)                               :: internalLength,FractureToughness
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      Type(SectionReal)                                  :: defaultSection
      PetscReal,Dimension(:),Pointer                     :: Gloc
      PetscReal,Dimension(:),Pointer                     :: alphaloc
      Type(MEF90_VECT)                                   :: gradAlphaElem
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscInt                                           :: iDoF1,iGauss
      PetscReal                                          :: C2
      PetscLogDouble                                     :: flops
           
      !C1 = FractureToughness / internalLength * .75_Kr 
      C2 = FractureToughness * internalLength * .75_Kr

      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(alphaLoc(elemType%numDof))
         Allocate(Gloc(elemType%numDof))
         Do cell = 1,size(cellID)     
            Gloc = 0.0_Kr
            Call SectionRealRestrictClosure(alpha,mesh,cellID(cell),elemType%numDof,alphaLoc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               gradAlphaElem = 0.0_Kr
               Do iDoF1 = 1,elemType%numDof
                  gradAlphaElem = gradAlphaElem + alphaLoc(iDof1) * elem(cell)%Grad_BF(iDoF1,iGauss)
               End Do
               Do iDoF1 = 1,elemType%numDof
                  Gloc(iDoF1) = Gloc(iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                  C2 * (elem(cell)%Grad_BF(iDoF1,iGauss) .DotP. gradAlphaElem)
               End Do
            End Do ! iGauss
            Call SectionRealUpdateClosure(G,mesh,cellID(cell),Gloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do ! cell
      
         flops = 3 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(Gloc)
         DeAllocate(alphaLoc)
      End If   
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageDamageOperatorSetAT1Elastic
   
#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageDamageRHSSetAT1Elastic"
!!!
!!!  
!!!  MEF90GradDamageDamageRHSSetAT1Elastic:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90GradDamageDamageRHSSetAT1Elastic(G,scaling,mesh,cellIS,internalLength,FractureToughness,elem,elemType,ierr)
      Type(SectionReal),Intent(INOUT)                    :: G
      PetscReal,Intent(IN)                               :: scaling
      Type(DM),Intent(IN)                                :: mesh
      Type(IS),Intent(IN)                                :: cellIS
      PetscReal,Intent(IN)                               :: internalLength
      PetscReal,Intent(IN)                               :: FractureToughness
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      Type(SectionReal)                                  :: defaultSection
      PetscReal,Dimension(:),Pointer                     :: Gloc
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscInt                                           :: iDoF1,iGauss
      PetscReal                                          :: C2
      PetscLogDouble                                     :: flops
           
      C2 = FractureToughness / internalLength * .375_Kr

      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(Gloc(elemType%numDof))
         Do cell = 1,size(cellID)      
            Gloc = 0.0_Kr
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               Do iDoF1 = 1,elemType%numDof
                  Gloc(iDoF1) = Gloc(iDoF1) - elem(cell)%Gauss_C(iGauss) * C2 * elem(cell)%BF(iDoF1,iGauss)
               End Do
            End Do ! iGauss
            Gloc = Gloc * scaling
            Call SectionRealUpdateClosure(G,mesh,cellID(cell),Gloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do ! cell
      
         flops = (3 * elemType%numDof + 1) * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(Gloc)
      End If   
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageDamageRHSSetAT1Elastic

#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageDamageBilinearFormSetAT1"
!!!
!!!  
!!!  MEF90GradDamageDamageBilinearFormSetAT1:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90GradDamageDamageBilinearFormSetAT1(K,mesh,meshDisp,cellIS,displacement,temperature,plasticStrain,internalLength,HookesLaw,LinearThermalExpansion,FractureToughness,elem,elemType,elemDisp,elemDispType,ierr)
      Type(Mat),Intent(IN)                               :: K 
      Type(DM),Intent(IN)                                :: mesh,meshDisp
      Type(IS),Intent(IN)                                :: cellIS
      Type(SectionReal),Intent(IN)                       :: displacement,temperature,plasticStrain
      PetscReal,Intent(IN)                               :: internalLength
      Type(MEF90_HOOKESLAW),Intent(IN)                   :: HookesLaw
      Type(MEF90_MATS),Intent(IN)                        :: LinearThermalExpansion
      PetscReal,Intent(IN)                               :: FractureToughness
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elem
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elemDisp
      Type(MEF90Element_Type),Intent(IN)                 :: elemType,elemDispType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      Type(SectionReal)                                  :: defaultSection
      PetscReal,Dimension(:),Pointer                     :: displacementloc,temperatureLoc,plasticStrainLoc
      PetscReal                                          :: temperatureElem
      Type(MEF90_MATS)                                   :: StrainElem,StressElem,plasticStrainElem
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscInt                                           :: iDoF1,iDoF2,iGauss
      PetscReal,Dimension(:,:),Pointer                   :: MatElem
      PetscReal                                          :: C2
      PetscLogDouble                                     :: flops

      C2 = FractureToughness * internalLength * .75_Kr

      Call DMMeshGetSectionReal(mesh,'default',defaultSection,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(displacementloc(elemDispType%numDof))
         Allocate(temperatureloc(elemType%numDof))
         Allocate(MatElem(elemType%numDof,elemType%numDof))
         Do cell = 1,size(cellID) 
            MatElem = 0.0_Kr  
            Call SectionRealRestrictClosure(displacement,meshDisp,cellID(cell),elemDispType%numDof,displacementloc,ierr);CHKERRQ(ierr)
            If (temperature%v /= 0) Then
               Call SectionRealRestrictClosure(temperature,mesh,cellID(cell),elemType%numDof,temperatureLoc,ierr);CHKERRQ(ierr)
            End If
            If (plasticStrain%v /= 0) Then
               Call SectionRealRestrict(plasticStrain,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
            End If
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               temperatureElem   = 0.0_Kr
               plasticStrainElem = 0.0_Kr
               strainElem        = 0.0_Kr
               Do iDoF1 = 1,elemDispType%numDof
                  strainElem = strainElem + displacementloc(iDof1) * elemDisp(cell)%GradS_BF(iDof1,iGauss)
               End Do
               If (temperature%v /= 0) Then
                  Do iDoF1 = 1,elemType%numDof
                     temperatureElem = temperatureElem + temperatureLoc(iDof1) * elem(cell)%BF(iDof1,iGauss)
                  End Do
                  strainElem = strainElem - (temperatureElem * LinearThermalExpansion)
               End If
               If (plasticStrain%v /= 0) Then
                  plasticStrainElem = plasticStrainLoc
                  strainElem = strainElem - plasticStrainElem
               End If
               stressElem = HookesLaw * strainElem
               Do iDoF1 = 1,elemType%numDof
                  Do iDoF2 = 1,elemType%numDof
                     MatElem(iDoF2,iDoF1) = MatElem(iDoF2,iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                    ( (stressElem .DotP. StrainElem) * elem(cell)%BF(iDof1,iGauss) * elem(cell)%BF(iDof2,iGauss) + & 
                                    C2 * (elem(cell)%Grad_BF(iDof1,iGauss) .dotP. elem(cell)%Grad_BF(iDof2,iGauss)))
                  End Do
               End Do
            End Do
            Call DMmeshAssembleMatrix(K,mesh,defaultSection,cellID(cell),MatElem,ADD_VALUES,ierr);CHKERRQ(ierr)
            If (plasticStrain%v /= 0) Then
               Call SectionRealRestore(plasticStrain,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
            End If
         End Do ! cell
         flops = 6 * elemType%numDof**2 * size(elem(1)%Gauss_C) * size(cellID) 
         If (temperature%v /= 0) Then
            flops = flops + 2 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID)
         End If
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(displacementloc)
         DeAllocate(temperatureloc)
         DeAllocate(MatElem)
      End If   
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)      
      Call SectionRealDestroy(defaultSection,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageDamageBilinearFormSetAT1

#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageDamageOperatorSetAT1"
!!!
!!!  
!!!  MEF90GradDamageDamageOperatorSetAT1:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90GradDamageDamageOperatorSetAT1(G,mesh,meshDisp,alpha,cellIS,displacement,temperature,plasticStrain,internalLength,HookesLaw,LinearThermalExpansion,FractureToughness,elem,elemType,elemDisp,elemDispType,ierr)
      Type(SectionReal),Intent(INOUT)                    :: G
      Type(DM),Intent(IN)                                :: mesh,meshDisp
      Type(SectionReal),Intent(IN)                       :: alpha
      Type(IS),Intent(IN)                                :: cellIS
      Type(SectionReal),Intent(IN)                       :: displacement,temperature,plasticStrain
      PetscReal,Intent(IN)                               :: internalLength
      Type(MEF90_HOOKESLAW),Intent(IN)                   :: HookesLaw
      Type(MEF90_MATS),Intent(IN)                        :: LinearThermalExpansion
      PetscReal,Intent(IN)                               :: FractureToughness
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elem
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elemDisp
      Type(MEF90Element_Type),Intent(IN)                 :: elemType,elemDispType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      Type(SectionReal)                                  :: defaultSection
      PetscReal,Dimension(:),Pointer                     :: Gloc
      PetscReal,Dimension(:),Pointer                     :: alphaloc
      PetscReal                                          :: alphaElem
      Type(MEF90_VECT)                                   :: gradAlphaElem
      PetscReal,Dimension(:),Pointer                     :: displacementloc,temperatureLoc,plasticStrainLoc
      PetscReal                                          :: temperatureElem
      Type(MEF90_MATS)                                   :: StrainElem,StressElem,plasticStrainElem
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscInt                                           :: iDoF1,iGauss
      PetscReal                                          :: C2
      PetscLogDouble                                     :: flops
           
      !C1 = FractureToughness / internalLength * .75_Kr 
      C2 = FractureToughness * internalLength * .75_Kr

      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(alphaLoc(elemType%numDof))
         Allocate(Gloc(elemType%numDof))
         Allocate(displacementLoc(elemDispType%numDof))
         Allocate(temperatureLoc(elemType%numDof))
         Do cell = 1,size(cellID)     
            Gloc = 0.0_Kr
            Call SectionRealRestrictClosure(alpha,mesh,cellID(cell),elemType%numDof,alphaLoc,ierr);CHKERRQ(ierr)
            Call SectionRealRestrictClosure(displacement,meshDisp,cellID(cell),elemDispType%numDof,displacementLoc,ierr);CHKERRQ(ierr)
            If (temperature%v /= 0) Then
               Call SectionRealRestrictClosure(temperature,mesh,cellID(cell),elemType%numDof,temperatureLoc,ierr);CHKERRQ(ierr)
            End If
            If (plasticStrain%v /= 0) Then
               Call SectionRealRestrict(plasticStrain,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
            End If
            
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               temperatureElem   = 0.0_Kr
               plasticStrainElem = 0.0_Kr
               strainElem        = 0.0_Kr
               Do iDoF1 = 1,elemDispType%numDof
                  strainElem = strainElem + displacementloc(iDof1) * elemDisp(cell)%GradS_BF(iDof1,iGauss)
               End Do
               If (temperature%v /= 0) Then
                  Do iDoF1 = 1,elemType%numDof
                     temperatureElem = temperatureElem + temperatureLoc(iDof1) * elem(cell)%BF(iDof1,iGauss)
                  End Do
                  strainElem = strainElem - (temperatureElem * LinearThermalExpansion)
               End If
               If (plasticStrain%v /= 0) Then
                  plasticStrainElem = plasticStrainLoc
                  strainElem = strainElem - plasticStrainElem
               End If
               stressElem = HookesLaw * strainElem

               alphaElem     = 0.0_Kr
               gradAlphaElem = 0.0_Kr
               Do iDoF1 = 1,elemType%numDof
                  alphaElem     = alphaElem     + alphaLoc(iDof1) * elem(cell)%BF(iDoF1,iGauss)
                  gradAlphaElem = gradAlphaElem + alphaLoc(iDof1) * elem(cell)%Grad_BF(iDoF1,iGauss)
               End Do
               Do iDoF1 = 1,elemType%numDof
                  Gloc(iDoF1) = Gloc(iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                 ((strainElem .dotP. stressElem) * elem(cell)%BF(iDoF1,iGauss) * alphaElem + &
                                  ( C2 * (elem(cell)%Grad_BF(iDoF1,iGauss) .DotP. gradAlphaElem)))
               End Do
            End Do ! iGauss
            Call SectionRealUpdateClosure(G,mesh,cellID(cell),Gloc,ADD_VALUES,ierr);CHKERRQ(iErr)
            If (plasticStrain%v /= 0) Then
               Call SectionRealRestore(plasticStrain,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
            End If
         End Do ! cell
      
         flops = 8 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         If (temperature%v /= 0) Then
            flops = flops + 2 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID)
         End If
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(temperatureLoc)
         DeAllocate(displacementLoc)
         DeAllocate(Gloc)
         DeAllocate(alphaLoc)
      End If   
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageDamageOperatorSetAT1
   
#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageDamageRHSSetAT1"
!!!
!!!  
!!!  MEF90GradDamageDamageRHSSetAT1:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90GradDamageDamageRHSSetAT1(G,scaling,mesh,meshDisp,cellIS,displacement,temperature,plasticStrain,internalLength,HookesLaw,LinearThermalExpansion,FractureToughness,elem,elemType,elemDisp,elemDispType,ierr)
      Type(SectionReal),Intent(INOUT)                    :: G
      PetscReal,Intent(IN)                               :: scaling
      Type(DM),Intent(IN)                                :: mesh,meshDisp
      Type(IS),Intent(IN)                                :: cellIS
      Type(SectionReal),Intent(IN)                       :: displacement,temperature,plasticStrain
      PetscReal,Intent(IN)                               :: internalLength
      Type(MEF90_HOOKESLAW),Intent(IN)                   :: HookesLaw
      Type(MEF90_MATS),Intent(IN)                        :: LinearThermalExpansion
      PetscReal,Intent(IN)                               :: FractureToughness
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elem
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elemDisp
      Type(MEF90Element_Type),Intent(IN)                 :: elemType,elemDispType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      Type(SectionReal)                                  :: defaultSection
      PetscReal,Dimension(:),Pointer                     :: Gloc
      PetscReal,Dimension(:),Pointer                     :: displacementloc,temperatureLoc,plasticStrainLoc
      PetscReal                                          :: temperatureElem
      Type(MEF90_MATS)                                   :: StrainElem,StressElem,plasticStrainElem
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscInt                                           :: iDoF1,iGauss
      PetscReal                                          :: C2
      PetscLogDouble                                     :: flops
           
      C2 = FractureToughness / internalLength * .375_Kr

      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(Gloc(elemType%numDof))
         Allocate(displacementLoc(elemDispType%numDof))
         Allocate(temperatureLoc(elemType%numDof))
         Do cell = 1,size(cellID)      
            Gloc = 0.0_Kr
            If (temperature%v /= 0) Then
               Call SectionRealRestrictClosure(temperature,mesh,cellID(cell),elemType%numDof,temperatureLoc,ierr);CHKERRQ(ierr)
            End If
            If (plasticStrain%v /= 0) Then
               Call SectionRealRestrict(plasticStrain,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
            End If
            Call SectionRealRestrictClosure(displacement,meshDisp,cellID(cell),elemDispType%numDof,displacementLoc,ierr);CHKERRQ(ierr)
            
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               temperatureElem   = 0.0_Kr
               plasticStrainElem = 0.0_Kr
               strainElem        = 0.0_Kr
               Do iDoF1 = 1,elemDispType%numDof
                  strainElem = strainElem + displacementloc(iDof1) * elemDisp(cell)%GradS_BF(iDof1,iGauss)
               End Do
               If (temperature%v /= 0) Then
                  Do iDoF1 = 1,elemType%numDof
                     temperatureElem = temperatureElem + temperatureLoc(iDof1) * elem(cell)%BF(iDof1,iGauss)
                  End Do
                  strainElem = strainElem - (temperatureElem * LinearThermalExpansion)
               End If
               If (plasticStrain%v /= 0) Then
                  plasticStrainElem = plasticStrainLoc
                  strainElem = strainElem - plasticStrainElem
               End If
               stressElem = HookesLaw * strainElem

               Do iDoF1 = 1,elemType%numDof
                  Gloc(iDoF1) = Gloc(iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                 ((strainElem .dotP. stressElem) - C2) * elem(cell)%BF(iDoF1,iGauss)
               End Do
            End Do ! iGauss
            Gloc = Gloc * scaling
            Call SectionRealUpdateClosure(G,mesh,cellID(cell),Gloc,ADD_VALUES,ierr);CHKERRQ(iErr)
            If (plasticStrain%v /= 0) Then
               Call SectionRealRestore(plasticStrain,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
            End If
         End Do ! cell
      
         flops = (4 * elemType%numDof * size(elem(1)%Gauss_C) + 1) * size(cellID) 
         If (temperature%v /= 0) Then
            flops = flops + 2 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         End If
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(temperatureLoc)
         DeAllocate(displacementLoc)
         DeAllocate(Gloc)
      End If   
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageDamageRHSSetAT1

#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageSurfaceEnergySetAT1"
!!!
!!!  
!!!  MEF90GradDamageSurfaceEnergySetAT1:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90GradDamageSurfaceEnergySetAT1(energy,alpha,meshScal,cellIS,internalLength,fractureToughness,elemScal,elemScalType,ierr)
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

      
      C1 = fractureToughness / internalLength * .375_Kr
      C2 = fractureToughness * internalLength * .375_Kr
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
               energy = energy + elemScal(cell)%Gauss_C(iGauss) * (alphaElem * C1 + (gradAlphaElem .dotP. gradAlphaElem) * C2)
            End Do ! Gauss
         End Do ! cell
         flops = (2 * elemScalType%numDof + 5 )* size(elemScal(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(alphaloc)
      End If   
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageSurfaceEnergySetAT1

#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageDamageBilinearFormSetAT2Elastic"
!!!
!!!  
!!!  MEF90GradDamageDamageBilinearFormSetAT2Elastic:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90GradDamageDamageBilinearFormSetAT2Elastic(K,mesh,cellIS,internalLength,FractureToughness,elem,elemType,ierr)
      Type(Mat),Intent(IN)                               :: K 
      Type(DM),Intent(IN)                                :: mesh
      Type(IS),Intent(IN)                                :: cellIS
      PetscReal,Intent(IN)                               :: internalLength
      PetscReal,Intent(IN)                               :: FractureToughness
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      Type(SectionReal)                                  :: defaultSection
      PetscReal,Dimension(:),Pointer                     :: Gloc
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscInt                                           :: iDoF1,iDoF2,iGauss
      PetscReal,Dimension(:,:),Pointer                   :: MatElem
      PetscReal                                          :: C1,C2
      PetscLogDouble                                     :: flops

      C1 = FractureToughness / internalLength
      C2 = FractureToughness * internalLength

      Call DMMeshGetSectionReal(mesh,'default',defaultSection,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(MatElem(elemType%numDof,elemType%numDof))
         Do cell = 1,size(cellID)   
            MatElem = 0.0_Kr
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               Do iDoF1 = 1,elemType%numDof
                  Do iDoF2 = 1,elemType%numDof
                     MatElem(iDoF2,iDoF1) = MatElem(iDoF2,iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                    (C1 * elem(cell)%BF(iDof1,iGauss) * elem(cell)%BF(iDof2,iGauss) + & 
                                     C2 * (elem(cell)%Grad_BF(iDof1,iGauss) .dotP. elem(cell)%Grad_BF(iDof2,iGauss)))
                  End Do
               End Do
            End Do
            Call DMmeshAssembleMatrix(K,mesh,defaultSection,cellID(cell),MatElem,ADD_VALUES,ierr);CHKERRQ(ierr)
         End Do ! cell
         flops = 6 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(MatElem)
      End If   
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)      
      Call SectionRealDestroy(defaultSection,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageDamageBilinearFormSetAT2Elastic

#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageDamageOperatorSetAT2Elastic"
!!!
!!!  
!!!  MEF90GradDamageDamageOperatorSetAT2Elastic:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90GradDamageDamageOperatorSetAT2Elastic(G,mesh,alpha,cellIS,internalLength,FractureToughness,elem,elemType,ierr)
      Type(SectionReal),Intent(INOUT)                    :: G
      Type(DM),Intent(IN)                                :: mesh
      Type(SectionReal),Intent(IN)                       :: alpha
      Type(IS),Intent(IN)                                :: cellIS
      PetscReal,Intent(IN)                               :: internalLength,FractureToughness
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elem
      Type(MEF90Element_Type),Intent(IN)                 :: elemType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      Type(SectionReal)                                  :: defaultSection
      PetscReal,Dimension(:),Pointer                     :: Gloc
      PetscReal,Dimension(:),Pointer                     :: alphaloc
      PetscReal                                          :: alphaElem
      Type(MEF90_VECT)                                   :: gradAlphaElem
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscInt                                           :: iDoF1,iGauss
      PetscReal                                          :: C1,C2
      PetscLogDouble                                     :: flops
           
      C1 = FractureToughness / internalLength
      C2 = FractureToughness * internalLength

      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(alphaLoc(elemType%numDof))
         Allocate(Gloc(elemType%numDof))
         Do cell = 1,size(cellID)      
            Gloc = 0.0_Kr
            Call SectionRealRestrictClosure(alpha,mesh,cellID(cell),elemType%numDof,alphaLoc,ierr);CHKERRQ(ierr)
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               alphaElem     = 0.0_Kr
               gradAlphaElem = 0.0_Kr
               Do iDoF1 = 1,elemType%numDof
                  alphaElem     = alphaElem     + alphaLoc(iDof1) * elem(cell)%BF(iDoF1,iGauss)
                  gradAlphaElem = gradAlphaElem + alphaLoc(iDof1) * elem(cell)%Grad_BF(iDoF1,iGauss)
               End Do
               Do iDoF1 = 1,elemType%numDof
                  Gloc(iDoF1) = Gloc(iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                 ( C1 * elem(cell)%BF(iDoF1,iGauss) * alphaElem + &
                                   C2 * (elem(cell)%Grad_BF(iDoF1,iGauss) .DotP. gradAlphaElem))
               End Do
            End Do ! iGauss
            Call SectionRealUpdateClosure(G,mesh,cellID(cell),Gloc,ADD_VALUES,ierr);CHKERRQ(iErr)
         End Do ! cell
      
         flops = 8 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(Gloc)
         DeAllocate(alphaLoc)
      End If   
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageDamageOperatorSetAT2Elastic

#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageDamageBilinearFormSetAT2"
!!!
!!!  
!!!  MEF90GradDamageDamageBilinearFormSetAT2:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90GradDamageDamageBilinearFormSetAT2(K,mesh,meshDisp,cellIS,displacement,temperature,plasticStrain,internalLength,HookesLaw,LinearThermalExpansion,FractureToughness,elem,elemType,elemDisp,elemDispType,ierr)
      Type(Mat),Intent(IN)                               :: K 
      Type(DM),Intent(IN)                                :: mesh,meshDisp
      Type(IS),Intent(IN)                                :: cellIS
      Type(SectionReal),Intent(IN)                       :: displacement,temperature,plasticStrain
      PetscReal,Intent(IN)                               :: internalLength
      Type(MEF90_HOOKESLAW),Intent(IN)                   :: HookesLaw
      Type(MEF90_MATS),Intent(IN)                        :: LinearThermalExpansion
      PetscReal,Intent(IN)                               :: FractureToughness
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elem
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elemDisp
      Type(MEF90Element_Type),Intent(IN)                 :: elemType,elemDispType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      Type(SectionReal)                                  :: defaultSection
      PetscReal,Dimension(:),Pointer                     :: displacementloc,temperatureLoc,plasticStrainLoc
      PetscReal                                          :: temperatureElem
      Type(MEF90_MATS)                                   :: StrainElem,StressElem,plasticStrainElem
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscInt                                           :: iDoF1,iDoF2,iGauss
      PetscReal,Dimension(:,:),Pointer                   :: MatElem
      PetscReal                                          :: C1,C2
      PetscLogDouble                                     :: flops

      C1 = FractureToughness / internalLength
      C2 = FractureToughness * internalLength

      Call DMMeshGetSectionReal(mesh,'default',defaultSection,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(displacementloc(elemDispType%numDof))
         Allocate(temperatureloc(elemType%numDof))
         Allocate(MatElem(elemType%numDof,elemType%numDof))
         Do cell = 1,size(cellID)   
            MatElem = 0.0_Kr
            Call SectionRealRestrictClosure(displacement,meshDisp,cellID(cell),elemDispType%numDof,displacementloc,ierr);CHKERRQ(ierr)
            If (temperature%v /= 0) Then
               Call SectionRealRestrictClosure(temperature,mesh,cellID(cell),elemType%numDof,temperatureLoc,ierr);CHKERRQ(ierr)
            End If
            If (plasticStrain%v /= 0) Then
               Call SectionRealRestrict(plasticStrain,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
            End If
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               temperatureElem   = 0.0_Kr
               plasticStrainElem = 0.0_Kr
               strainElem        = 0.0_Kr
               Do iDoF1 = 1,elemDispType%numDof
                  strainElem = strainElem + displacementloc(iDof1) * elemDisp(cell)%GradS_BF(iDof1,iGauss)
               End Do
               If (temperature%v /= 0) Then
                  Do iDoF1 = 1,elemType%numDof
                     temperatureElem = temperatureElem + temperatureLoc(iDof1) * elem(cell)%BF(iDof1,iGauss)
                  End Do
                  strainElem = strainElem - (temperatureElem * LinearThermalExpansion)
               End If
               If (plasticStrain%v /= 0) Then
                  plasticStrainElem = plasticStrainLoc
                  strainElem = strainElem - plasticStrainElem
               End If
               stressElem = HookesLaw * strainElem
               Do iDoF1 = 1,elemType%numDof
                  Do iDoF2 = 1,elemType%numDof
                     MatElem(iDoF2,iDoF1) = MatElem(iDoF2,iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                    (((stressElem .DotP. StrainElem) + C1) * elem(cell)%BF(iDof1,iGauss) * elem(cell)%BF(iDof2,iGauss) + & 
                                    C2 * (elem(cell)%Grad_BF(iDof1,iGauss) .dotP. elem(cell)%Grad_BF(iDof2,iGauss)))
                  End Do
               End Do
            End Do
            Call DMmeshAssembleMatrix(K,mesh,defaultSection,cellID(cell),MatElem,ADD_VALUES,ierr);CHKERRQ(ierr)
            If (plasticStrain%v /= 0) Then
               Call SectionRealRestore(plasticStrain,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
            End If
         End Do ! cell
         flops = 7 * elemType%numDof**2 * size(elem(1)%Gauss_C) * size(cellID) 
         If (temperature%v /= 0) Then
            flops = flops + 2 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         End If
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(displacementloc)
         DeAllocate(temperatureloc)
         DeAllocate(MatElem)
      End If   
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)      
      Call SectionRealDestroy(defaultSection,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageDamageBilinearFormSetAT2

#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageDamageOperatorSetAT2"
!!!
!!!  
!!!  MEF90GradDamageDamageOperatorSetAT2:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90GradDamageDamageOperatorSetAT2(G,mesh,meshDisp,alpha,cellIS,displacement,temperature,plasticStrain,internalLength,HookesLaw,LinearThermalExpansion,FractureToughness,elem,elemType,elemDisp,elemDispType,ierr)
      Type(SectionReal),Intent(INOUT)                    :: G
      Type(DM),Intent(IN)                                :: mesh,meshDisp
      Type(SectionReal),Intent(IN)                       :: alpha
      Type(IS),Intent(IN)                                :: cellIS
      Type(SectionReal),Intent(IN)                       :: displacement,temperature,plasticStrain
      PetscReal,Intent(IN)                               :: internalLength
      Type(MEF90_HOOKESLAW),Intent(IN)                   :: HookesLaw
      Type(MEF90_MATS),Intent(IN)                        :: LinearThermalExpansion
      PetscReal,Intent(IN)                               :: FractureToughness
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elem
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elemDisp
      Type(MEF90Element_Type),Intent(IN)                 :: elemType,elemDispType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      Type(SectionReal)                                  :: defaultSection
      PetscReal,Dimension(:),Pointer                     :: Gloc
      PetscReal,Dimension(:),Pointer                     :: alphaloc
      PetscReal                                          :: alphaElem
      Type(MEF90_VECT)                                   :: gradAlphaElem
      PetscReal,Dimension(:),Pointer                     :: displacementloc,temperatureLoc,plasticStrainLoc
      PetscReal                                          :: temperatureElem
      Type(MEF90_MATS)                                   :: StrainElem,StressElem,plasticStrainElem
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscInt                                           :: iDoF1,iGauss
      PetscReal                                          :: C1,C2
      PetscLogDouble                                     :: flops
           
      C1 = FractureToughness / internalLength
      C2 = FractureToughness * internalLength

      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(alphaLoc(elemType%numDof))
         Allocate(Gloc(elemType%numDof))
         Allocate(displacementLoc(elemDispType%numDof))
         Allocate(temperatureLoc(elemType%numDof))
         Do cell = 1,size(cellID)      
            Gloc = 0.0_Kr
            Call SectionRealRestrictClosure(alpha,mesh,cellID(cell),elemType%numDof,alphaLoc,ierr);CHKERRQ(ierr)
            Call SectionRealRestrictClosure(displacement,meshDisp,cellID(cell),elemDispType%numDof,displacementLoc,ierr);CHKERRQ(ierr)
            If (temperature%v /= 0) Then
               Call SectionRealRestrictClosure(temperature,mesh,cellID(cell),elemType%numDof,temperatureLoc,ierr);CHKERRQ(ierr)
            End If
            If (plasticStrain%v /= 0) Then
               Call SectionRealRestrict(plasticStrain,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
            End If
            
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               temperatureElem   = 0.0_Kr
               plasticStrainElem = 0.0_Kr
               strainElem        = 0.0_Kr
               Do iDoF1 = 1,elemDispType%numDof
                  strainElem = strainElem + displacementloc(iDof1) * elemDisp(cell)%GradS_BF(iDof1,iGauss)
               End Do
               If (temperature%v /= 0) Then
                  Do iDoF1 = 1,elemType%numDof
                     temperatureElem = temperatureElem + temperatureLoc(iDof1) * elem(cell)%BF(iDof1,iGauss)
                  End Do
                  strainElem = strainElem - (temperatureElem * LinearThermalExpansion)
               End If
               If (plasticStrain%v /= 0) Then
                  plasticStrainElem = plasticStrainLoc
                  strainElem = strainElem - plasticStrainElem
               End If
               stressElem = HookesLaw * strainElem

               alphaElem     = 0.0_Kr
               gradAlphaElem = 0.0_Kr
               Do iDoF1 = 1,elemType%numDof
                  alphaElem     = alphaElem     + alphaLoc(iDof1) * elem(cell)%BF(iDoF1,iGauss)
                  gradAlphaElem = gradAlphaElem + alphaLoc(iDof1) * elem(cell)%Grad_BF(iDoF1,iGauss)
               End Do
               Do iDoF1 = 1,elemType%numDof
                  Gloc(iDoF1) = Gloc(iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                 (((strainElem .dotP. stressElem) + C1)* elem(cell)%BF(iDoF1,iGauss) * alphaElem + &
                                  ( C2 * (elem(cell)%Grad_BF(iDoF1,iGauss) .DotP. gradAlphaElem)) )
               End Do
            End Do ! iGauss
            Call SectionRealUpdateClosure(G,mesh,cellID(cell),Gloc,ADD_VALUES,ierr);CHKERRQ(iErr)
            If (plasticStrain%v /= 0) Then
               Call SectionRealRestore(plasticStrain,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
            End If
         End Do ! cell
      
         flops = 9 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         If (temperature%v /= 0) Then
            flops = flops + 2 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         End If
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(temperatureLoc)
         DeAllocate(displacementLoc)
         DeAllocate(Gloc)
         DeAllocate(alphaLoc)
      End If   
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageDamageOperatorSetAT2

#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageDamageRHSSetAT2"
!!!
!!!  
!!!  MEF90GradDamageDamageRHSSetAT2:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90GradDamageDamageRHSSetAT2(G,scaling,mesh,meshDisp,cellIS,displacement,temperature,plasticStrain,internalLength,HookesLaw,LinearThermalExpansion,FractureToughness,elem,elemType,elemDisp,elemDispType,ierr)
      Type(SectionReal),Intent(INOUT)                    :: G
      PetscReal,Intent(IN)                               :: scaling
      Type(DM),Intent(IN)                                :: mesh,meshDisp
      Type(IS),Intent(IN)                                :: cellIS
      Type(SectionReal),Intent(IN)                       :: displacement,temperature,plasticStrain
      PetscReal,Intent(IN)                               :: internalLength
      Type(MEF90_HOOKESLAW),Intent(IN)                   :: HookesLaw
      Type(MEF90_MATS),Intent(IN)                        :: LinearThermalExpansion
      PetscReal,Intent(IN)                               :: FractureToughness
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elem
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elemDisp
      Type(MEF90Element_Type),Intent(IN)                 :: elemType,elemDispType
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      Type(SectionReal)                                  :: defaultSection
      PetscReal,Dimension(:),Pointer                     :: Gloc
      Type(MEF90_VECT)                                   :: gradAlphaElem
      PetscReal,Dimension(:),Pointer                     :: displacementloc,temperatureLoc,plasticStrainLoc
      PetscReal                                          :: temperatureElem
      Type(MEF90_MATS)                                   :: StrainElem,StressElem,plasticStrainElem
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscInt                                           :: iDoF1,iGauss
      PetscLogDouble                                     :: flops
           
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      If (Size(cellID) > 0) Then
         Allocate(Gloc(elemType%numDof))
         Allocate(displacementLoc(elemDispType%numDof))
         Allocate(temperatureLoc(elemType%numDof))
         Do cell = 1,size(cellID)      
            Gloc = 0.0_Kr
            If (temperature%v /= 0) Then
               Call SectionRealRestrictClosure(temperature,mesh,cellID(cell),elemType%numDof,temperatureLoc,ierr);CHKERRQ(ierr)
            End If
            If (plasticStrain%v /= 0) Then
               Call SectionRealRestrict(plasticStrain,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
            End If
            Call SectionRealRestrictClosure(displacement,meshDisp,cellID(cell),elemDIspType%numDof,displacementLoc,ierr);CHKERRQ(ierr)
            
            Do iGauss = 1,size(elem(cell)%Gauss_C)
               temperatureElem   = 0.0_Kr
               plasticStrainElem = 0.0_Kr
               strainElem        = 0.0_Kr
               Do iDoF1 = 1,elemDispType%numDof
                  strainElem = strainElem + displacementloc(iDof1) * elemDisp(cell)%GradS_BF(iDof1,iGauss)
               End Do
               If (temperature%v /= 0) Then
                  Do iDoF1 = 1,elemType%numDof
                     temperatureElem = temperatureElem + temperatureLoc(iDof1) * elem(cell)%BF(iDof1,iGauss)
                  End Do
                  strainElem = strainElem - (temperatureElem * LinearThermalExpansion)
               End If
               If (plasticStrain%v /= 0) Then
                  plasticStrainElem = plasticStrainLoc
                  strainElem = strainElem - plasticStrainElem
               End If
               stressElem = HookesLaw * strainElem

               Do iDoF1 = 1,elemType%numDof
                  Gloc(iDoF1) = Gloc(iDoF1) + elem(cell)%Gauss_C(iGauss) * &
                                (strainElem .dotP. stressElem) * elem(cell)%BF(iDoF1,iGauss) 
               End Do
            End Do ! iGauss
            Gloc = Gloc * scaling
            Call SectionRealUpdateClosure(G,mesh,cellID(cell),Gloc,ADD_VALUES,ierr);CHKERRQ(iErr)
            If (plasticStrain%v /= 0) Then
               Call SectionRealRestore(plasticStrain,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
            End If
         End Do ! cell
      
         flops = (3 * elemType%numDof * size(elem(1)%Gauss_C) + 1) * size(cellID) 
         If (temperature%v /= 0) Then
            flops = flops + 2 * elemType%numDof * size(elem(1)%Gauss_C) * size(cellID) 
         End If
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(temperatureLoc)
         DeAllocate(displacementLoc)
         DeAllocate(Gloc)
      End If   
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageDamageRHSSetAT2
   
#undef __FUNCT__
#define __FUNCT__ "MEF90GradDamageSurfaceEnergySetAT2"
!!!
!!!  
!!!  MEF90GradDamageSurfaceEnergySetAT2:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90GradDamageSurfaceEnergySetAT2(energy,alpha,meshScal,cellIS,internalLength,fractureToughness,elemScal,elemScalType,ierr)
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

      
      C1 = fractureToughness / internalLength * .5_Kr
      C2 = fractureToughness * internalLength * .5_Kr
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
               energy = energy + elemScal(cell)%Gauss_C(iGauss) * (alphaElem**2 * C1 + (gradAlphaElem .dotP. gradAlphaElem) * C2)
            End Do ! Gauss
         End Do ! cell
         flops = (2 * elemScalType%numDof + 6)* size(elemScal(1)%Gauss_C) * size(cellID) 
         Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
         DeAllocate(alphaloc)
      End If   
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
   End Subroutine MEF90GradDamageSurfaceEnergySetAT2
End Module MEF90_APPEND(m_MEF90_GradDamageImplementation_,MEF90_DIM)D
