#include "../MEF90/mef90.inc"
Module MEF90_APPEND(m_MEF90_DefMechPlasticity,MEF90_DIM)D
#include "finclude/petscdef.h"
   use m_MEF90
   use m_MEF90_DefMechCtx
   implicit NONE
   private
   public MEF90DefMechPlasticStrainUpdate

   !!! note that this type is NOT C interoperable, which is not an issue, since we only
   !!! need SNLP to carry its address

   type :: MEF90DefMechPlasticityCtx
      Type(MEF90_HOOKESLAW)       :: HookesLaw
      real(Kind = Kr)             :: YieldStress
      real(Kind = Kr)             :: DuctileCouplingPower
      Type(MEF90_MATS)            :: InelasticStrain
      Type(MEF90_MATS)            :: PlasticStrainOld
      real(Kind = Kr)             :: Damage
      real(Kind = Kr)             :: residualStiffness
      real(Kind = Kr)             :: CoefficientLinSoft
      real(Kind = Kr)             :: CoefficientDruckerPrager
      real(Kind = Kr)             :: CoefficientDruckerPragerCapModel1
      real(Kind = Kr)             :: CoefficientDruckerPragerCapModel2
      real(Kind = Kr)             :: CoefficientDruckerPragerCapModel3
      logical(Kind = Kr)          :: isPlaneStress
   end type MEF90DefMechPlasticityCtx

contains
!!! Since these functions are C interoperable, fortran cannot rename them in a use statement
!!! We play with the pre-processor in order to avoid duplicate symbols.

#define FHG_NONE MEF90_APPEND(fhg_None,MEF90_DIM)D

#define FHG_VONMISES MEF90_APPEND(fhg_VonMises,MEF90_DIM)D

#define FHG_VONMISESPLANESTRESS MEF90_APPEND(fhg_VonMisesPlaneStress,MEF90_DIM)D

#define FHG_VONMISES1D MEF90_APPEND(fhg_VonMises1D,MEF90_DIM)D

#define FHG_DRUCKERPRAGER MEF90_APPEND(fhg_DruckerPrager,MEF90_DIM)D

#define FHG_DRUCKERPRAGERCAPMODEL MEF90_APPEND(fhg_DruckerPragerCapModel,MEF90_DIM)D

#define FHG_TRESCA MEF90_APPEND(fhg_Tresca,MEF90_DIM)D

#undef __FUNCT__
#define __FUNCT__ "FHG_NONE"
!!!
!!!  
!!!  fhg: VonMises
!!!  
!!!  (c) 2015 Erwan Tanne : erwan.tanne@gmail.com
!!!
!!!

   subroutine FHG_NONE(x,f,h,g,myctx) bind(c)
      use,intrinsic :: iso_c_binding
      use m_MEF90

      real(kind=c_double)                       :: x(*)
      real(kind=c_double)                       :: f(*)
      real(kind=c_double)                       :: h(*)
      real(kind=c_double)                       :: g(*)
      type(c_ptr),intent(in),value              :: myctx

      CONTINUE
   end subroutine FHG_NONE

#undef __FUNCT__
#define __FUNCT__ "FHG_VONMISES"
!!!
!!!  
!!!  fhg: VonMises
!!!  
!!!  (c) 2015 Erwan Tanne : erwan.tanne@gmail.com
!!!
!!!

   subroutine FHG_VONMISES(x,f,h,g,myctx) bind(c)
      use,intrinsic :: iso_c_binding
      use m_MEF90

      real(kind=c_double)                       :: x(*)
      real(kind=c_double)                       :: f(*)
      real(kind=c_double)                       :: h(*)
      real(kind=c_double)                       :: g(*)
      real(Kind = Kr)                           :: StiffnessA         !! Stiffness = a(alpha)/b(alpha)
      real(Kind = Kr)                           :: StiffnessB
      type(c_ptr),intent(in),value              :: myctx

      type(MEF90DefMechPlasticityCtx),pointer   :: myctx_ptr
      type(MEF90_MATS)                          :: xMatS
      type(MEF90_MATS)                          :: Stress


      xMatS = x(1:SIZEOFMEF90_MATS)
      !!! This is the fortran equivalent of casting ctx into a c_ptr
      call c_f_pointer(myctx,myctx_ptr)


      !!! Select which softening young model
      if (myctx_ptr%CoefficientLinSoft==0) then
         StiffnessA = (1.0_Kr - myctx_ptr%Damage)**2 + myctx_ptr%residualStiffness
         StiffnessB = (1.0_Kr - myctx_ptr%Damage)**myctx_ptr%DuctileCouplingPower + myctx_ptr%residualStiffness
      else 
         StiffnessA = ( (1.0_Kr - myctx_ptr%Damage)**2 /( 1.0_Kr + ( myctx_ptr%CoefficientLinSoft - 1.0_Kr )*(1.0_Kr - (1.0_Kr - myctx_ptr%Damage)**2 ) ) ) + myctx_ptr%residualStiffness
         StiffnessB = (1.0_Kr - myctx_ptr%Damage)**myctx_ptr%DuctileCouplingPower + myctx_ptr%residualStiffness
      endif

      Stress=myctx_ptr%HookesLaw*(myctx_ptr%InelasticStrain-xMatS)
      f(1) = ( (myctx_ptr%HookesLaw *(xMatS-myctx_ptr%PlasticStrainOld)) .DotP. (xMatS-myctx_ptr%PlasticStrainOld) ) * StiffnessA / 2.0
      g(1) = StiffnessA * sqrt( MEF90_DIM / (MEF90_DIM - 1.0_kr)  * ( deviatoricPart(Stress)  .DotP.  deviatoricPart(Stress) ))  - myctx_ptr%YieldStress*StiffnessB
      h(1) = Trace(xMatS)
   end subroutine FHG_VONMISES



#undef __FUNCT__
#define __FUNCT__ "FHG_VONMISESPLANESTRESS"
!!!
!!!  
!!!  fhg: VonMises
!!!  
!!!  (c) 2015 Erwan Tanne : erwan.tanne@gmail.com
!!!
!!!

   subroutine FHG_VONMISESPLANESTRESS(x,f,h,g,myctx) bind(c)
      use,intrinsic :: iso_c_binding
      use m_MEF90

      real(kind=c_double)                       :: x(*)
      real(kind=c_double)                       :: f(*)
      real(kind=c_double)                       :: h(*)
      real(kind=c_double)                       :: g(*)
      real(Kind = Kr)                           :: StiffnessA         !! Stiffness = a(alpha)/b(alpha)
      real(Kind = Kr)                           :: StiffnessB
      type(c_ptr),intent(in),value              :: myctx

      type(MEF90DefMechPlasticityCtx),pointer   :: myctx_ptr
      type(MEF90_MATS)                          :: xMatS
      real(Kind = Kr)                           :: lambda,mu,E,nu
      type(MatS3D)                              :: Strain,PlasticStrainFlow,Stress

      xMatS = x(1:SIZEOFMEF90_MATS)
      !!! This is the fortran equivalent of casting ctx into a c_ptr
      call c_f_pointer(myctx,myctx_ptr)

      !!! Select which softening young model
      if (myctx_ptr%CoefficientLinSoft==0) then
         StiffnessA = (1.0_Kr - myctx_ptr%Damage)**2
         StiffnessB = (1.0_Kr - myctx_ptr%Damage)**myctx_ptr%DuctileCouplingPower
      else 
         StiffnessA = ( (1.0_Kr - myctx_ptr%Damage)**2 /( 1.0_Kr + ( myctx_ptr%CoefficientLinSoft - 1.0_Kr )*(1.0_Kr - (1.0_Kr - myctx_ptr%Damage)**2 ) ) )
         StiffnessB = (1.0_Kr - myctx_ptr%Damage)**myctx_ptr%DuctileCouplingPower
      endif

      E     = myctx_ptr%HookesLaw%YoungsModulus
      nu    = myctx_ptr%HookesLaw%PoissonRatio
      mu    = E / (1.0_Kr + nu) * .5_Kr
      lambda  = E * nu / (1.0_Kr + nu) / (1 - 2.0_Kr * nu)

      Strain      = 0.0_Kr
      Strain%XX   = myctx_ptr%InelasticStrain%XX
      Strain%YY   = myctx_ptr%InelasticStrain%YY
      Strain%XY   = myctx_ptr%InelasticStrain%XY
      Strain%ZZ   = (-lambda*Trace(myctx_ptr%InelasticStrain) - 2*mu*Trace(myctx_ptr%PlasticStrainOld))/(lambda + 2*mu)

      PlasticStrainFlow    = 0.0_Kr
      PlasticStrainFlow%XX = xMatS%XX-myctx_ptr%PlasticStrainOld%XX
      PlasticStrainFlow%YY = xMatS%YY-myctx_ptr%PlasticStrainOld%YY
      PlasticStrainFlow%XY = xMatS%XY-myctx_ptr%PlasticStrainOld%XY
      PlasticStrainFlow%ZZ = -( PlasticStrainFlow%XX + PlasticStrainFlow%YY )

      Stress    = 0.0_Kr
      Stress%XX = lambda*(Trace(Strain)) + 2*mu*(Strain%XX-xMatS%XX)
      Stress%YY = lambda*(Trace(Strain)) + 2*mu*(Strain%YY-xMatS%YY)
      Stress%XY = 2*mu*(Strain%XY-xMatS%XY)
      Stress%ZZ = lambda*(Trace(Strain)) + 2*mu*(Strain%ZZ + Trace(xMatS) )

      f(1) = mu*(  (PlasticStrainFlow .DotP. PlasticStrainFlow) ) * StiffnessA 
      g(1) = StiffnessA * sqrt( (3.0/2.0)*( deviatoricPart(Stress) .dotP. deviatoricPart(Stress) ) ) - myctx_ptr%YieldStress*StiffnessB

   end subroutine FHG_VONMISESPLANESTRESS


#undef __FUNCT__
#define __FUNCT__ "FHG_VONMISES1D"
!!!
!!!  
!!!  fhg: VonMises
!!!  
!!!  (c) 2015 Erwan Tanne : erwan.tanne@gmail.com
!!!
!!!

   subroutine FHG_VONMISES1D(x,f,h,g,myctx) bind(c)
      use,intrinsic :: iso_c_binding
      use m_MEF90

      real(kind=c_double)                       :: x(*)
      real(kind=c_double)                       :: f(*)
      real(kind=c_double)                       :: h(*)
      real(kind=c_double)                       :: g(*)
      real(Kind = Kr)                           :: StiffnessA         !! Stiffness = a(alpha)/b(alpha)
      real(Kind = Kr)                           :: StiffnessB
      real(Kind = Kr)                           :: Stiffness
      type(c_ptr),intent(in),value              :: myctx

      type(MEF90DefMechPlasticityCtx),pointer   :: myctx_ptr
      type(MEF90_MATS)                          :: xMatS,Stress

      xMatS = x(1:SIZEOFMEF90_MATS)
      !!! This is the fortran equivalent of casting ctx into a c_ptr
      call c_f_pointer(myctx,myctx_ptr)

       !!! Select which softening young model
      if (myctx_ptr%CoefficientLinSoft==0) then
         StiffnessA = (1.0_Kr - myctx_ptr%Damage)**2 + myctx_ptr%residualStiffness
         StiffnessB = (1.0_Kr - myctx_ptr%Damage)**myctx_ptr%DuctileCouplingPower + myctx_ptr%residualStiffness
      else 
         StiffnessA = ( (1.0_Kr - myctx_ptr%Damage)**2 /( 1.0_Kr + ( myctx_ptr%CoefficientLinSoft - 1.0_Kr )*(1.0_Kr - (1.0_Kr - myctx_ptr%Damage)**2 ) ) ) + myctx_ptr%residualStiffness
         StiffnessB = (1.0_Kr - myctx_ptr%Damage)**myctx_ptr%DuctileCouplingPower + myctx_ptr%residualStiffness
      endif

      Stress = myctx_ptr%HookesLaw*(myctx_ptr%InelasticStrain-xMatS)

      f(1) = ( (myctx_ptr%HookesLaw *(xMatS-myctx_ptr%PlasticStrainOld)) .DotP. (xMatS-myctx_ptr%PlasticStrainOld) ) * StiffnessA / 2.0
      g(1) = StiffnessA * sqrt( MEF90_DIM / (MEF90_DIM - 1.0_kr)  * ( deviatoricPart(Stress)  .DotP.  deviatoricPart(Stress) ))  - myctx_ptr%YieldStress*StiffnessB
      h(1) = xMatS%YY
      h(2) = xMatS%XY
   end subroutine FHG_VONMISES1D




#undef __FUNCT__
#define __FUNCT__ "FHG_DRUCKERPRAGER"
!!!
!!!  
!!!  fhg: VonMises
!!!  
!!!  (c) 2015 Erwan Tanne : erwan.tanne@gmail.com
!!!
!!!

   subroutine FHG_DRUCKERPRAGER(x,f,h,g,myctx) bind(c)
      use,intrinsic :: iso_c_binding
      use m_MEF90

      real(kind=c_double)                       :: x(*)
      real(kind=c_double)                       :: f(*)
      real(kind=c_double)                       :: h(*)
      real(kind=c_double)                       :: g(*)
      type(c_ptr),intent(in),value              :: myctx

      type(MEF90DefMechPlasticityCtx),pointer   :: myctx_ptr
      type(MEF90_MATS)                          :: xMatS

      xMatS = x(1:SIZEOFMEF90_MATS)
      !!! This is the fortran equivalent of casting ctx into a c_ptr
      call c_f_pointer(myctx,myctx_ptr)


      f(1) = ( (myctx_ptr%HookesLaw * (xMatS-myctx_ptr%PlasticStrainOld)) .DotP. (xMatS-myctx_ptr%PlasticStrainOld) ) /2.0 
      g(1) =  sqrt( MEF90_DIM / (MEF90_DIM - 1.0_kr) * ( deviatoricPart(myctx_ptr%HookesLaw*(myctx_ptr%InelasticStrain-xMatS))  .DotP.  deviatoricPart(myctx_ptr%HookesLaw*(myctx_ptr%InelasticStrain-xMatS)) ))  - myctx_ptr%YieldStress -  myctx_ptr%CoefficientDruckerPrager*Trace(myctx_ptr%HookesLaw*(myctx_ptr%InelasticStrain-xMatS))

   end subroutine FHG_DRUCKERPRAGER


#undef __FUNCT__
#define __FUNCT__ "FHG_DRUCKERPRAGERCAPMODEL"
!!!
!!!  
!!!  fhg: VonMises
!!!  
!!!  (c) 2015 Erwan Tanne : erwan.tanne@gmail.com
!!!
!!!

   subroutine FHG_DRUCKERPRAGERCAPMODEL(x,f,h,g,myctx) bind(c)
      use,intrinsic :: iso_c_binding
      use m_MEF90

      real(kind=c_double)                       :: x(*)
      real(kind=c_double)                       :: f(*)
      real(kind=c_double)                       :: h(*)
      real(kind=c_double)                       :: g(*)
      type(c_ptr),intent(in),value              :: myctx

      type(MEF90DefMechPlasticityCtx),pointer   :: myctx_ptr
      type(MEF90_MATS)                          :: xMatS

      xMatS = x(1:SIZEOFMEF90_MATS)
      !!! This is the fortran equivalent of casting ctx into a c_ptr
      call c_f_pointer(myctx,myctx_ptr)


      f(1) = ( (myctx_ptr%HookesLaw * (xMatS-myctx_ptr%PlasticStrainOld)) .DotP. (xMatS-myctx_ptr%PlasticStrainOld) ) /2.0 
      g(1) =  sqrt( MEF90_DIM / (MEF90_DIM - 1.0_kr) * ( deviatoricPart(myctx_ptr%HookesLaw*(myctx_ptr%InelasticStrain-xMatS))  .DotP.  deviatoricPart(myctx_ptr%HookesLaw*(myctx_ptr%InelasticStrain-xMatS)) ))  - myctx_ptr%YieldStress -  myctx_ptr%CoefficientDruckerPrager*Trace(myctx_ptr%HookesLaw*(myctx_ptr%InelasticStrain-xMatS))
      g(2) =  sqrt( MEF90_DIM / (MEF90_DIM - 1.0_kr) * ( deviatoricPart(myctx_ptr%HookesLaw*(myctx_ptr%InelasticStrain-xMatS))  .DotP.  deviatoricPart(myctx_ptr%HookesLaw*(myctx_ptr%InelasticStrain-xMatS)) ))  - myctx_ptr%CoefficientDruckerPragerCapModel3 -  myctx_ptr%CoefficientDruckerPragerCapModel2*Trace(myctx_ptr%HookesLaw*(myctx_ptr%InelasticStrain-xMatS)) - myctx_ptr%CoefficientDruckerPragerCapModel1*Trace(myctx_ptr%HookesLaw*(myctx_ptr%InelasticStrain-xMatS))**2.0

   end subroutine FHG_DRUCKERPRAGERCAPMODEL





#undef __FUNCT__
#define __FUNCT__ "FHG_TRESCA"

!!!
!!!  
!!!  fhg: Tresca
!!!  
!!!  (c) 2015 Erwan Tanne : erwan.tanne@gmail.com
!!!
!!!

   subroutine FHG_TRESCA(x,f,h,g,myctx) bind(c)
      use,intrinsic :: iso_c_binding
      use m_MEF90

      real(kind=c_double)                       :: x(*)
      real(kind=c_double)                       :: f(*)
      real(kind=c_double)                       :: h(*)
      real(kind=c_double)                       :: g(*)
      type(c_ptr),intent(in),value              :: myctx

      type(MEF90DefMechPlasticityCtx),pointer   :: myctx_ptr
      type(MEF90_MATS)                          :: xMatS
      type(MEF90_MAT)                           :: MatProj
      type(MEF90_MATS)                          :: MatDiag
      type(MEF90_MAT)                           :: MatPrincipal
      type(MEF90_MATS)                          :: gMatS

      xMatS = x(1:SIZEOFMEF90_MATS)
      !!! This is the fortran equivalent of casting ctx into a c_ptr
      call c_f_pointer(myctx,myctx_ptr)


      !write(*,*) 'A.e(u):         ', myctx_ptr%HookesLaw*myctx_ptr%Strain
      ! D=P^(-1).A.P 
      call EigenVectorValues(deviatoricPart(myctx_ptr%HookesLaw*myctx_ptr%InelasticStrain),MatProj,MatDiag)
      MatPrincipal = Transpose(MatProj)*MatSymToMat(deviatoricPart(myctx_ptr%HookesLaw*myctx_ptr%InelasticStrain) - xMatS)*MatProj

      f(1) = ( (myctx_ptr%HookesLaw * (xMatS-myctx_ptr%PlasticStrainOld)) .DotP. (xMatS-myctx_ptr%PlasticStrainOld) ) /2.
      h(1) = Trace(xMatS)

#if MEF90_DIM == 2
      g(1) = +(MatPrincipal%XX-MatPrincipal%YY) - myctx_ptr%YieldStress
      g(2) = -(MatPrincipal%XX-MatPrincipal%YY) - myctx_ptr%YieldStress
      g(3) = +(MatPrincipal%YY)                 - myctx_ptr%YieldStress
      g(4) = -(MatPrincipal%YY)                 - myctx_ptr%YieldStress
      g(5) = +(MatPrincipal%XX)                 - myctx_ptr%YieldStress
      g(6) = -(MatPrincipal%XX)                 - myctx_ptr%YieldStress
#else
      Write(*,*) 'Tresca3D is NOT implemented'
#endif
   end subroutine FHG_TRESCA


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechPlasticStrainUpdate"
!!!
!!!  
!!!  MEF90DefMechPlasticStrainUpdate:
!!!  
!!!  (c) 2015 Erwan Tanne : erwan.tanne@gmail.com
!!!

   Subroutine MEF90DefMechPlasticStrainUpdate(MEF90DefMechCtx,plasticStrain,x,PlasticStrainOld,plasticStrainPrevious,cumulatedDissipatedPlasticEnergyVariation,ierr)
      use,intrinsic :: iso_c_binding
#ifdef MEF90_HAVE_SNLP
      use SNLPF90
#endif

      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      Type(Vec),Intent(INOUT)                            :: plasticStrain
      Type(Vec),Intent(IN)                               :: x,PlasticStrainOld,plasticStrainPrevious,cumulatedDissipatedPlasticEnergyVariation
      PetscErrorCode,Intent(OUT)                         :: ierr

#ifdef MEF90_HAVE_SNLP
      Type(DM)                                           :: Mesh
      Type(SectionReal)                                  :: plasticStrainSec,plasticStrainOldSec,inelasticStrainSec,plasticStrainPreviousSec,cumulatedDissipatedPlasticEnergyVariationSec
      PetscReal,Dimension(:),Pointer                     :: plasticStrainLoc,plasticStrainOldLoc,inelasticStrainLoc,damageLoc,plasticStrainPreviousLoc
      PetscReal,Dimension(:),Pointer                     :: cumulatedDissipatedPlasticEnergyVariationLoc
      type(c_funptr)                                     :: snlp_fhg,snlp_Dfhg
      integer(kind=c_int)                                :: snlp_n,snlp_m,snlp_p
      type(SNLP),pointer                                 :: s
      integer                                            :: i,j
      integer(kind=c_int)                                :: exit_code

      type(MEF90DefMechPlasticityCtx),target             :: PlasticityCtx
      type(c_ptr)                                        :: snlp_ctx


      PetscInt                                           :: dim,set,cell,QuadratureOrder
      Type(IS)                                           :: cellSetGlobalIS,setIS
      PetscInt,dimension(:),Pointer                      :: setID,cellID
      Type(MEF90DefMechCellSetOptions_Type),pointer      :: cellSetOptions
      Type(MEF90Element_Type)                            :: elemDisplacementType
      Type(MEF90_MATPROP),pointer                        :: matPropSet
      Type(SectionReal)                                  :: xSec,temperatureSec
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemScal
      Type(MEF90Element_Type)                            :: elemScalType

      Type(SectionReal)                                  :: damageSec
      PetscReal                                          :: damageCellAvg,Stiffness

      Type(MEF90_MATS)                                   :: PlasticStrainMatS


      Call SectionRealDuplicate(MEF90DefMechCtx%cellDMMatSSec,plasticStrainSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(plasticStrainSec,MEF90DefMechCtx%cellDMMatSScatter,SCATTER_REVERSE,plasticStrain,ierr);CHKERRQ(ierr) 

      Call SectionRealDuplicate(MEF90DefMechCtx%cellDMMatSSec,plasticStrainOldSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(plasticStrainOldSec,MEF90DefMechCtx%cellDMMatSScatter,SCATTER_REVERSE,plasticStrainOld,ierr);CHKERRQ(ierr) 

      Call SectionRealDuplicate(MEF90DefMechCtx%cellDMMatSSec,inelasticStrainSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(plasticStrainOldSec,MEF90DefMechCtx%cellDMMatSScatter,SCATTER_REVERSE,plasticStrainOld,ierr);CHKERRQ(ierr) 

      Call SectionRealDuplicate(MEF90DefMechCtx%cellDMMatSSec,plasticStrainPreviousSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(plasticStrainOldSec,MEF90DefMechCtx%cellDMMatSScatter,SCATTER_REVERSE,plasticStrainOld,ierr);CHKERRQ(ierr) 

      Call SectionRealDuplicate(MEF90DefMechCtx%cellDMScalSec,cumulatedDissipatedPlasticEnergyVariationSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(cumulatedDissipatedPlasticEnergyVariationSec,MEF90DefMechCtx%cellDMScalScatter,SCATTER_REVERSE,cumulatedDissipatedPlasticEnergyVariation,ierr);CHKERRQ(ierr) 

      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,xSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(xSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr)

      If (Associated(MEF90DefMechCtx%damage)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,damageSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(damageSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)   
      Else
         damageSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,temperatureSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(temperatureSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr)   
      Else
         temperatureSec%v = 0
      End If

      Call DMMeshGetDimension(MEF90DefMechCtx%DM,dim,ierr);CHKERRQ(ierr)
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         
         !!GET DAMAGE TYPE
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         If (cellSetOptions%plasticityType /= MEF90DefMech_plasticityTypeNONE) Then
            elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)      
            elemScalType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)

            Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(ierr)
            Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
            If ((Size(cellID) > 0) .AND. (elemDisplacementType%coDim == 0)) Then
               !!! Call proper local assembly depending on the type of damage law
               Select Case (cellSetOptions%plasticityType)
                  case(MEF90DefMech_plasticityTypeVonMises)
                     snlp_Dfhg = c_null_funptr
                     snlp_fhg  = c_funloc(FHG_VONMISES)
                     snlp_n    = SIZEOFMEF90_MATS
                     snlp_m    = 1
                     snlp_p    = 1
                     snlp_ctx  = c_loc(PlasticityCtx)

                     case(MEF90DefMech_plasticityTypeVonMisesPlaneStress)
                     snlp_Dfhg = c_null_funptr
                     snlp_fhg  = c_funloc(FHG_VONMISESPLANESTRESS)
                     snlp_n    = SIZEOFMEF90_MATS
                     snlp_m    = 0
                     snlp_p    = 1
                     snlp_ctx  = c_loc(PlasticityCtx)

                     case(MEF90DefMech_plasticityTypeVonMises1D)
                     snlp_Dfhg = c_null_funptr
                     snlp_fhg  = c_funloc(FHG_VONMISES1D)
                     snlp_n    = SIZEOFMEF90_MATS
                     snlp_m    = 2
                     snlp_p    = 1
                     snlp_ctx  = c_loc(PlasticityCtx)

                  case(MEF90DefMech_plasticityTypeTresca)
                     snlp_Dfhg = c_null_funptr
                     snlp_fhg  = c_funloc(FHG_TRESCA)
                     snlp_n    = SIZEOFMEF90_MATS
                     snlp_m    = 1
                     snlp_p    = 2*SIZEOFMEF90_MATS
                     snlp_ctx  = c_loc(PlasticityCtx)

                  case(MEF90DefMech_plasticityTypeDruckerPrager)
                     snlp_Dfhg = c_null_funptr
                     snlp_fhg  = c_funloc(FHG_DRUCKERPRAGER)
                     snlp_n    = SIZEOFMEF90_MATS
                     snlp_m    = 0
                     snlp_p    = 1
                     snlp_ctx  = c_loc(PlasticityCtx)

                  case(MEF90DefMech_plasticityTypeDruckerPragerCapModel)
                     snlp_Dfhg = c_null_funptr
                     snlp_fhg  = c_funloc(FHG_DRUCKERPRAGERCAPMODEL)
                     snlp_n    = SIZEOFMEF90_MATS
                     snlp_m    = 0
                     snlp_p    = 2
                     snlp_ctx  = c_loc(PlasticityCtx)

                  case(MEF90DefMech_plasticityTypeNONE)
                     snlp_Dfhg = c_null_funptr
                     snlp_fhg  = c_funloc(FHG_NONE)
                     snlp_n    = 1
                     snlp_m    = 0
                     snlp_p    = 0
                     snlp_ctx  = c_loc(PlasticityCtx)
                  case default
                     Print*,__FUNCT__,': Unimplemented plasticity Type',cellSetOptions%PlasticityType
                     STOP 
               End select

               If (cellSetOptions%damageType == MEF90DefMech_damageTypeLinSoft) Then
                  PlasticityCtx%CoefficientLinSoft = matpropSet%CoefficientLinSoft
               Else
                  PlasticityCtx%CoefficientLinSoft = 0.0_Kr
               End If
               PlasticityCtx%HookesLaw = matpropSet%HookesLaw
               PlasticityCtx%residualStiffness = matpropSet%residualStiffness
               PlasticityCtx%YieldStress = matpropSet%YieldStress
               PlasticityCtx%DuctileCouplingPower = matpropSet%DuctileCouplingPower
               PlasticityCtx%CoefficientDruckerPrager = matpropSet%CoefficientDruckerPrager
               PlasticityCtx%CoefficientDruckerPragerCapModel1 = matpropSet%CoefficientDruckerPragerCapModel1
               PlasticityCtx%CoefficientDruckerPragerCapModel2 = matpropSet%CoefficientDruckerPragerCapModel2
               PlasticityCtx%CoefficientDruckerPragerCapModel3 = matpropSet%CoefficientDruckerPragerCapModel3


#if MEF90_DIM == 2
   PlasticityCtx%isPlaneStress = matPropSet%HookesLaw%isPlaneStress
#else 
   PlasticityCtx%isPlaneStress = .FALSE.
#endif


               Call SNLPNew(s,snlp_n,snlp_m,snlp_p,snlp_fhg,snlp_Dfhg,snlp_ctx)
               QuadratureOrder = 2 * (elemDisplacementType%order - 1)
               Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
               Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
               Call MEF90InelasticStrainSet(inelasticStrainSec,xSec,temperatureSec,MEF90DefMechCtx%DMVect,MEF90DefMechCtx%DMScal,setIS,matpropSet%LinearThermalExpansion, &
                                            elemDisplacement,elemDisplacementType,elemScal,elemScalType,ierr)
               Allocate(damageloc(elemScalType%numDof))
               Do cell = 1,size(cellID)
                  !! actualiser le ctx (  HookesLaw ,InelasticStrainSec, plasticStrainStrainSec, plasticStrainOldSec  )
                  Call SectionRealRestrict(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
                  Call SectionRealRestrict(plasticStrainOldSec,cellID(cell),plasticStrainOldLoc,ierr);CHKERRQ(ierr)
                  Call SectionRealRestrict(inelasticStrainSec,cellID(cell),inelasticStrainLoc,ierr);CHKERRQ(ierr)
                  Call SectionRealRestrict(cumulatedDissipatedPlasticEnergyVariationSec,cellID(cell),cumulatedDissipatedPlasticEnergyVariationLoc,ierr);CHKERRQ(ierr)

                  If (Associated(MEF90DefMechCtx%damage)) Then
                     Call SectionRealRestrictClosure(damageSec,MEF90DefMechCtx%DMScal,cellID(cell),elemScalType%numDof,damageLoc,ierr);CHKERRQ(ierr)
                     Call MEF90GradDamageCellAverage(damageCellAvg,damageLoc,elemScal(cell),elemScalType,ierr)
                  Else
                     damageCellAvg = 0.0_Kr
                  End If

                  PlasticityCtx%Damage = damageCellAvg 
                  PlasticityCtx%PlasticStrainOld = plasticStrainOldLoc
                  PlasticityCtx%InelasticStrain = InelasticStrainLoc
                  s%show_progress = 0
         
                  !!! This is a bit dangerous:
                  !!! If PetscReal is not the same as c_double, this call will fail
                  !!! Brittle in traction, Ductile in compression
                  Select Case(cellSetOptions%unilateralContactType)
                  Case (MEF90DefMech_unilateralContactTypeBrittleDuctile)
                        If (Trace(PlasticityCtx%InelasticStrain) > 0.0_Kr ) Then
                           plasticStrainLoc = plasticStrainOldLoc
                           Call SectionRealRestrict(plasticStrainPreviousSec,cellID(cell),plasticStrainPreviousLoc,ierr);CHKERRQ(ierr)
                           plasticStrainPreviousLoc = plasticStrainLoc
                           Call SectionRealRestore(plasticStrainPreviousSec,cellID(cell),plasticStrainPreviousLoc,ierr);CHKERRQ(ierr)
                        Else 
                           exit_code = SNLPL1SQP(s,plasticStrainLoc)
                        End if
                  Case default
                     if (cellSetOptions%plasticityType /= MEF90DefMech_plasticityTypeNONE) then
                        exit_code = SNLPL1SQP(s,plasticStrainLoc)
                     end if
                  End Select

                  !!! cumulatedDissipatedPlasticEnergy
                  PlasticStrainMatS = plasticStrainLoc
                  stiffness = 0.0_Kr
                  Select Case (cellSetOptions%damageType)
                     Case (MEF90DefMech_damageTypeAT1Elastic,MEF90DefMech_damageTypeAT2Elastic)
                        Stiffness = 1.0_Kr
                     Case (MEF90DefMech_damageTypeAT1,MEF90DefMech_damageTypeAT2)
                        Stiffness = (1.0_Kr - PlasticityCtx%Damage)**(2.0_Kr-PlasticityCtx%DuctileCouplingPower)  + PlasticityCtx%residualStiffness**(2.0_Kr-PlasticityCtx%DuctileCouplingPower)
                     Case (MEF90DefMech_damageTypeLinSoft)
                        Stiffness = ( (1.0_Kr - PlasticityCtx%Damage)**(2.0_Kr - PlasticityCtx%DuctileCouplingPower) / ( 1.0_Kr + ( PlasticityCtx%CoefficientLinSoft - 1.0_Kr )*(1.0_Kr - (1.0_Kr - PlasticityCtx%Damage)**2.0_kr ) ) ) + PlasticityCtx%residualStiffness**(2.0_Kr-PlasticityCtx%DuctileCouplingPower)
                     Case default
                        Print*,__FUNCT__,': Unimplemented damage Type, only AT1Elastic and AT2Elastic implement',cellSetOptions%damageType
                        STOP  
                  End Select
                  cumulatedDissipatedPlasticEnergyVariationLoc(1) = Stiffness * ( PlasticityCtx%HookesLaw * ( PlasticityCtx%InelasticStrain - PlasticStrainMatS ) ) .dotP. ( PlasticStrainMatS - PlasticityCtx%plasticStrainOld )
                  Call SectionRealRestore(cumulatedDissipatedPlasticEnergyVariationSec,cellID(cell),cumulatedDissipatedPlasticEnergyVariationLoc,ierr);CHKERRQ(ierr)

                  Call SectionRealRestore(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
                  Call SectionRealRestore(plasticStrainOldSec,cellID(cell),plasticStrainOldLoc,ierr);CHKERRQ(ierr)
                  Call SectionRealRestore(inelasticStrainSec,cellID(cell),inelasticStrainLoc,ierr);CHKERRQ(ierr)
               End Do !cell
               Call MEF90Element_Destroy(elemDisplacement,ierr)
               Call MEF90Element_Destroy(elemScal,ierr)
               Call SNLPDelete(s)
               DeAllocate(damageLoc)
            End If ! set 
            Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If !damageType /= NONE
      End Do !! set
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)

      !!! forward data plasticStrain & cumulatedDissipatedPlasticEnergy
      Call SectionRealToVec(plasticStrainSec,MEF90DefMechCtx%cellDMMatSScatter,SCATTER_FORWARD,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(cumulatedDissipatedPlasticEnergyVariationSec,MEF90DefMechCtx%cellDMScalScatter,SCATTER_FORWARD,cumulatedDissipatedPlasticEnergyVariation,ierr);CHKERRQ(ierr)

      Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(plasticStrainOldSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(inelasticStrainSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(plasticStrainPreviousSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(cumulatedDissipatedPlasticEnergyVariationSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
      If (Associated(MEF90DefMechCtx%damage)) Then
         Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)
      End If

      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
      End If
#else
      write(*,*) 'This example needs SNLP'
#endif
      End Subroutine MEF90DefMechPlasticStrainUpdate
End Module MEF90_APPEND(m_MEF90_DefMechPlasticity,MEF90_DIM)D
