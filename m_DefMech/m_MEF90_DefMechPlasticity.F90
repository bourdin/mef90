#include "../MEF90/mef90.inc"
#include "mef90DefMech.inc"
Module MEF90_APPEND(m_MEF90_DefMechPlasticity,MEF90_DIM)D
#include "petsc/finclude/petsc.h"
   use m_MEF90
   use m_MEF90_DefMechCtx
   use MEF90_APPEND(m_MEF90_DefMechPlasticityCtx,MEF90_DIM)D
   use MEF90_APPEND(m_MEF90_DefMechPlasticityNone,MEF90_DIM)D
   use MEF90_APPEND(m_MEF90_DefMechPlasticityCap,MEF90_DIM)D
   use MEF90_APPEND(m_MEF90_DefMechPlasticityCrystalBCC,MEF90_DIM)D
   use MEF90_APPEND(m_MEF90_DefMechPlasticityCrystalSingleSlip,MEF90_DIM)D
   use MEF90_APPEND(m_MEF90_DefMechPlasticityDruckerPragerCap,MEF90_DIM)D
   use MEF90_APPEND(m_MEF90_DefMechPlasticityGreen,MEF90_DIM)D
   use MEF90_APPEND(m_MEF90_DefMechPlasticityGurson,MEF90_DIM)D
   use MEF90_APPEND(m_MEF90_DefMechPlasticityHill,MEF90_DIM)D
   use MEF90_APPEND(m_MEF90_DefMechPlasticityTresca,MEF90_DIM)D
   use MEF90_APPEND(m_MEF90_DefMechPlasticityVonMises,MEF90_DIM)D

   implicit NONE
   private
   public MEF90DefMechPlasticStrainUpdate

contains

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechPlasticStrainUpdate"
!!!
!!!
!!!  MEF90DefMechPlasticStrainUpdate:
!!!
!!!  (c) 2015 Erwan Tanne : erwan.tanne@gmail.com
!!!  (c) 2017 Stella Brach : brach.ste@gmail.com
!!!  (c) 2022 Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90DefMechPlasticStrainUpdate(MEF90DefMechCtx,plasticStrain,x,PlasticStrainOld,plasticStrainPrevious,cumulatedDissipatedPlasticEnergyVariation,cumulatedDissipatedPlasticEnergyOld,ierr)
      use,intrinsic :: iso_c_binding
#ifdef MEF90_HAVE_SNLP
      use SNLPF90
#endif

      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      Type(tVec),Intent(INOUT)                           :: plasticStrain
      Type(tVec),Intent(IN)                              :: x,PlasticStrainOld,plasticStrainPrevious,cumulatedDissipatedPlasticEnergyVariation,cumulatedDissipatedPlasticEnergyOld
      PetscErrorCode,Intent(INOUT)                       :: ierr

! #ifdef MEF90_HAVE_SNLP
!       Type(DM)                                           :: Mesh
!       Type(SectionReal)                                  :: plasticStrainSec,plasticStrainOldSec,totalStrainSec,plasticStrainPreviousSec,cumulatedDissipatedPlasticEnergyVariationSec,cumulatedDissipatedPlasticEnergyOldSec
!       PetscReal,Dimension(:),Pointer                     :: plasticStrainLoc,plasticStrainOldLoc,totalStrainLoc,damageLoc,plasticStrainPreviousLoc
!       PetscReal,Dimension(:),Pointer                     :: cumulatedDissipatedPlasticEnergyVariationLoc,cumulatedDissipatedPlasticEnergyOldLoc
!       type(c_funptr)                                     :: snlp_fhg,snlp_Dfhg
!       integer(kind=c_int)                                :: snlp_n,snlp_m,snlp_p
!       type(SNLP),pointer                                 :: s
!       integer                                            :: i,j
!       integer(kind=c_int)                                :: exit_code

!       type(MEF90DefMechPlasticityCtx),target             :: PlasticityCtx
!       type(c_ptr)                                        :: snlp_ctx

!       PetscInt                                           :: dim,set,cell,QuadratureOrder
!       Type(tIS)                                           :: cellSetGlobalIS,setIS
!       PetscInt,dimension(:),Pointer                      :: setID,cellID
!       Type(MEF90DefMechCellSetOptions_Type),pointer      :: cellSetOptions
!       Type(MEF90Element_Type)                            :: elemDisplacementType
!       Type(MEF90_MATPROP),pointer                        :: matPropSet
!       Type(SectionReal)                                  :: xSec,temperatureSec
!       Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
!       Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemScal
!       Type(MEF90Element_Type)                            :: elemScalType

!       Type(SectionReal)                                  :: damageSec
!       PetscReal                                          :: damageCellAvg,Stiffness

!       Type(MEF90_MATS)                                   :: PlasticStrainMatS
!       PetscReal                                          :: Sigma_33_PlaneStrain

!       Call SectionRealDuplicate(MEF90DefMechCtx%cellDMMatSSec,plasticStrainSec,ierr);CHKERRQ(ierr)
!       Call SectionRealToVec(plasticStrainSec,MEF90DefMechCtx%cellDMMatSScatter,SCATTER_REVERSE,plasticStrain,ierr);CHKERRQ(ierr)

!       Call SectionRealDuplicate(MEF90DefMechCtx%cellDMMatSSec,plasticStrainOldSec,ierr);CHKERRQ(ierr)
!       Call SectionRealToVec(plasticStrainOldSec,MEF90DefMechCtx%cellDMMatSScatter,SCATTER_REVERSE,plasticStrainOld,ierr);CHKERRQ(ierr)

!       Call SectionRealDuplicate(MEF90DefMechCtx%cellDMMatSSec,totalStrainSec,ierr);CHKERRQ(ierr)
!       Call SectionRealToVec(plasticStrainOldSec,MEF90DefMechCtx%cellDMMatSScatter,SCATTER_REVERSE,plasticStrainOld,ierr);CHKERRQ(ierr)

!       Call SectionRealDuplicate(MEF90DefMechCtx%cellDMMatSSec,plasticStrainPreviousSec,ierr);CHKERRQ(ierr)
!       Call SectionRealToVec(plasticStrainOldSec,MEF90DefMechCtx%cellDMMatSScatter,SCATTER_REVERSE,plasticStrainOld,ierr);CHKERRQ(ierr)

!       Call SectionRealDuplicate(MEF90DefMechCtx%cellDMScalSec,cumulatedDissipatedPlasticEnergyVariationSec,ierr);CHKERRQ(ierr)
!       Call SectionRealToVec(cumulatedDissipatedPlasticEnergyVariationSec,MEF90DefMechCtx%cellDMScalScatter,SCATTER_REVERSE,cumulatedDissipatedPlasticEnergyVariation,ierr);CHKERRQ(ierr)

!       Call SectionRealDuplicate(MEF90DefMechCtx%cellDMScalSec,cumulatedDissipatedPlasticEnergyOldSec,ierr);CHKERRQ(ierr)
!       Call SectionRealToVec(cumulatedDissipatedPlasticEnergyOldSec,MEF90DefMechCtx%cellDMScalScatter,SCATTER_REVERSE,cumulatedDissipatedPlasticEnergyOld,ierr);CHKERRQ(ierr)

!       Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,xSec,ierr);CHKERRQ(ierr)
!       Call SectionRealToVec(xSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr)

!       If (Associated(MEF90DefMechCtx%damage)) Then
!          Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,damageSec,ierr);CHKERRQ(ierr)
!          Call SectionRealToVec(damageSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)
!       Else
!          damageSec%v = 0
!       End If

!       If (Associated(MEF90DefMechCtx%temperature)) Then
!          Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,temperatureSec,ierr);CHKERRQ(ierr)
!          Call SectionRealToVec(temperatureSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr)
!       Else
!          temperatureSec%v = 0
!       End If

!       Call DMMeshGetDimension(MEF90DefMechCtx%DM,dim,ierr);CHKERRQ(ierr)
!       Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
!       Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)
!       Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

!       Do set = 1,size(setID)
!          Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)

!          !!GET DAMAGE TYPE
!          Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
!          If (cellSetOptions%plasticityType /= MEF90DefMech_plasticityTypeNONE) Then
!             elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
!             elemScalType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)

!             Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(ierr)
!             Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
!             If ((Size(cellID) > 0) .AND. (elemDisplacementType%coDim == 0)) Then
!                !!! Call proper local assembly depending on the type of damage law
!                Select Case (cellSetOptions%plasticityType)
!                   case(MEF90DefMech_plasticityTypeVonMises)
!                      snlp_Dfhg = c_null_funptr
!                      snlp_fhg  = c_funloc(FHG_VONMISES)
!                      snlp_n    = SIZEOFMEF90_MATS
!                      snlp_m    = 1
!                      snlp_p    = 2
!                      snlp_ctx  = c_loc(PlasticityCtx)

!                      case(MEF90DefMech_plasticityTypeVonMisesPlaneTheory)
!                      snlp_Dfhg = c_null_funptr
!                      snlp_fhg  = c_funloc(FHG_VONMISESPLANETHEORY)
!                      snlp_n    = SIZEOFMEF90_MATS
!                      snlp_m    = 0
!                      snlp_p    = 1
!                      snlp_ctx  = c_loc(PlasticityCtx)

!                      case(MEF90DefMech_plasticityTypeVonMises1D)
!                      snlp_Dfhg = c_null_funptr
!                      snlp_fhg  = c_funloc(FHG_VONMISES1D)
!                      snlp_n    = SIZEOFMEF90_MATS
!                      snlp_m    = 2
!                      snlp_p    = 2
!                      snlp_ctx  = c_loc(PlasticityCtx)

!                   case(MEF90DefMech_plasticityTypeTresca)
!                      snlp_Dfhg = c_null_funptr
!                      snlp_fhg  = c_funloc(FHG_TRESCA)
!                      snlp_n    = SIZEOFMEF90_MATS
!                      snlp_m    = 1
!                      snlp_p    = 2*SIZEOFMEF90_MATS
!                      snlp_ctx  = c_loc(PlasticityCtx)

!                   case(MEF90DefMech_plasticityTypeCapModel)
!                      snlp_Dfhg = c_null_funptr
!                      snlp_fhg  = c_funloc(FHG_CAPMODEL)
!                      snlp_n    = SIZEOFMEF90_MATS
!                      snlp_m    = 0
!                      snlp_p    = 1
!                      snlp_ctx  = c_loc(PlasticityCtx)

!                   case(MEF90DefMech_plasticityTypeDruckerPragerCapModel)
!                      snlp_Dfhg = c_null_funptr
!                      snlp_fhg  = c_funloc(FHG_DRUCKERPRAGERCAPMODEL)
!                      snlp_n    = SIZEOFMEF90_MATS
!                      snlp_m    = 0
!                      snlp_p    = 2
!                      snlp_ctx  = c_loc(PlasticityCtx)

!                   case(MEF90DefMech_plasticityTypeHillPlaneTheory)
!                      snlp_Dfhg = c_null_funptr
!                      snlp_fhg  = c_funloc(FHG_HILLPLANETHEORY)
!                      snlp_n    = SIZEOFMEF90_MATS
!                      snlp_m    = 0
!                      snlp_p    = 1
!                      snlp_ctx  = c_loc(PlasticityCtx)

!                   case(MEF90DefMech_plasticityTypeGreen)
!                      snlp_Dfhg = c_null_funptr
!                      snlp_fhg  = c_funloc(FHG_GREEN)
!                      snlp_n    = SIZEOFMEF90_MATS
!                      snlp_m    = 1
!                      snlp_p    = 1
!                      snlp_ctx  = c_loc(PlasticityCtx)

!                   case(MEF90DefMech_plasticityTypeGurson)
!                      snlp_Dfhg = c_null_funptr
!                      snlp_fhg  = c_funloc(FHG_GURSON)
!                      snlp_n    = SIZEOFMEF90_MATS
!                      snlp_m    = 0
!                      snlp_p    = 1
!                      snlp_ctx  = c_loc(PlasticityCtx)

!                   case(MEF90DefMech_plasticityTypeNONE)
!                      snlp_Dfhg = c_null_funptr
!                      snlp_fhg  = c_funloc(FHG_NONE)
!                      snlp_n    = 1
!                      snlp_m    = 0
!                      snlp_p    = 0
!                      snlp_ctx  = c_loc(PlasticityCtx)

!                   case(MEF90DefMech_plasticityTypeCrystalSingleSlip)
!                      snlp_Dfhg = c_null_funptr
!                      snlp_fhg  = c_funloc(FHG_CRYSTALSINGLESLIP)
!                      snlp_n    = SIZEOFMEF90_MATS
! #if MEF90_DIM==2
!                      if (matPropSet%HookesLaw%isPlaneStress) then
!                         write(*,*) "Plane stress CrystalSingleSlip is not implemented"
!                      end if
!                      snlp_m    = 3
! #elif MEF90_DIM==3
!                      snlp_m    = 1
! #endif
!                      if (matpropSet%isViscousPlasticity) then
!                         snlp_p = 0
!                      else
!                         snlp_p = 0
!                      end if
!                      snlp_ctx  = c_loc(PlasticityCtx)

!                   case(MEF90DefMech_plasticityTypeCrystalBCC)
!                      snlp_Dfhg = c_null_funptr
!                      snlp_fhg  = c_funloc(FHG_CRYSTALBCC)
!                      snlp_n    = SIZEOFMEF90_MATS
! #if MEF90_DIM==2
!                      if (matPropSet%HookesLaw%isPlaneStress) then
!                         write(*,*) "Plane stress CrystalBCC is not implemented"
!                      end if
!                      if (matpropSet%isViscousPlasticity) then
!                         snlp_m    = 3
!                      else
!                         snlp_m    = 1
!                      end if
! #elif MEF90_DIM==3
!                      snlp_m    = 1
! #endif
!                      if (matpropSet%isViscousPlasticity) then
!                         snlp_p = 0
!                      else   
!                         snlp_p = 1
!                      end if
!                      snlp_ctx  = c_loc(PlasticityCtx)
!                   case default
!                      Print*,__FUNCT__,': Unimplemented plasticity Type',cellSetOptions%PlasticityType
!                      STOP
!                End select

!                If (cellSetOptions%damageType == MEF90DefMech_damageTypeLinSoft) Then
!                   PlasticityCtx%CoefficientLinSoft = cellSetOptions%DamageATLinSoftk
!                Else
!                   PlasticityCtx%CoefficientLinSoft = 0.0_Kr
!                End If
!                PlasticityCtx%HookesLaw = matpropSet%HookesLaw
!                PlasticityCtx%residualStiffness = matpropSet%residualStiffness
!                PlasticityCtx%YieldStress = matpropSet%YieldStress
!                PlasticityCtx%DuctileCouplingPower = matpropSet%DuctileCouplingPower
!                PlasticityCtx%isNoPlCoupling = matpropSet%isNoPlCoupling
!                PlasticityCtx%CoefficientDruckerPrager = matpropSet%CoefficientDruckerPrager
!                PlasticityCtx%CoefficientCapModel0 = matpropSet%CoefficientCapModel0
!                PlasticityCtx%CoefficientCapModel1 = matpropSet%CoefficientCapModel1
!                PlasticityCtx%CoefficientCapModel2 = matpropSet%CoefficientCapModel2
!                PlasticityCtx%CoefficientCapModelD = matpropSet%CoefficientCapModelD
!                PlasticityCtx%residualYieldStress = matpropSet%residualYieldStress
!                PlasticityCtx%isLinearIsotropicHardening = matpropSet%isLinearIsotropicHardening
!                PlasticityCtx%CoeffF = matpropSet%CoeffF
!                PlasticityCtx%CoeffG = matpropSet%CoeffG
!                PlasticityCtx%CoeffH = matpropSet%CoeffH
!                PlasticityCtx%CoeffM = matpropSet%CoeffM
!                PlasticityCtx%CoeffN = matpropSet%CoeffN
!                PlasticityCtx%CoeffL = matpropSet%CoeffL
!                PlasticityCtx%YieldTau0 = matpropSet%YieldTau0
!                PlasticityCtx%residualYieldTau0 = matpropSet%residualYieldTau0
!                PlasticityCtx%phi1 = matpropSet%phi1
!                PlasticityCtx%phi2 = matpropSet%phi2
!                PlasticityCtx%Phi = matpropSet%Phi
!                PlasticityCtx%RotationMatrix3D = matpropSet%RotationMatrix
!                PlasticityCtx%isViscousPlasticity = matpropSet%isViscousPlasticity
!                PlasticityCtx%ViscosityGamma0 = matpropSet%ViscosityGamma0
!                PlasticityCtx%ViscosityN = matpropSet%ViscosityN
!                PlasticityCtx%Viscositydt = matpropSet%Viscositydt
!                PlasticityCtx%m = matpropSet%m


! #if MEF90_DIM == 2
!    PlasticityCtx%isPlaneStress = matPropSet%HookesLaw%isPlaneStress
! #else
!    PlasticityCtx%isPlaneStress = .FALSE.
! #endif


!                Call SNLPNew(s,snlp_n,snlp_m,snlp_p,snlp_fhg,snlp_Dfhg,snlp_ctx)
!                QuadratureOrder = 2 * (elemDisplacementType%order - 1)
!                Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
!                Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
!                Call MEF90totalStrainSet(totalStrainSec,xSec,temperatureSec,MEF90DefMechCtx%DMVect,MEF90DefMechCtx%DMScal,setIS,matpropSet%LinearThermalExpansion, &
!                                             elemDisplacement,elemDisplacementType,elemScal,elemScalType,ierr)
!                Allocate(damageloc(elemScalType%numDof))
!                Do cell = 1,size(cellID)
!                   !! actualiser le ctx (  HookesLaw ,totalStrainSec, plasticStrainStrainSec, plasticStrainOldSec  )
!                   Call SectionRealRestrict(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
!                   Call SectionRealRestrict(plasticStrainOldSec,cellID(cell),plasticStrainOldLoc,ierr);CHKERRQ(ierr)
!                   Call SectionRealRestrict(totalStrainSec,cellID(cell),totalStrainLoc,ierr);CHKERRQ(ierr)
!                   Call SectionRealRestrict(cumulatedDissipatedPlasticEnergyVariationSec,cellID(cell),cumulatedDissipatedPlasticEnergyVariationLoc,ierr);CHKERRQ(ierr)
!                   Call SectionRealRestrict(cumulatedDissipatedPlasticEnergyOldSec,cellID(cell),cumulatedDissipatedPlasticEnergyOldLoc,ierr);CHKERRQ(ierr)


!                   If (Associated(MEF90DefMechCtx%damage)) Then
!                      Call SectionRealRestrictClosure(damageSec,MEF90DefMechCtx%DMScal,cellID(cell),elemScalType%numDof,damageLoc,ierr);CHKERRQ(ierr)
!                      Call MEF90GradDamageCellAverage(damageCellAvg,damageLoc,elemScal(cell),elemScalType,ierr)
!                   Else
!                      damageCellAvg = 0.0_Kr
!                   End If

!                   PlasticityCtx%Damage = damageCellAvg
!                   PlasticityCtx%PlasticStrainOld = plasticStrainOldLoc
!                   PlasticityCtx%totalStrain = totalStrainLoc
!                   PlasticityCtx%plasticStrainPrevious = plasticStrainLoc
!                   PlasticityCtx%cumulatedDissipatedPlasticEnergy = cumulatedDissipatedPlasticEnergyOldLoc(1)

!                   s%show_progress = 0

!                   !!! This is a bit dangerous:
!                   !!! If PetscReal is not the same as c_double, this call will fail
!                   !!! Brittle in traction, Ductile in compression
!                   Select Case(cellSetOptions%unilateralContactType)
!                   ! Case (MEF90DefMech_unilateralContactTypeBrittleDuctile)
!                   !       If (Trace(PlasticityCtx%totalStrain) > 0.0_Kr ) Then
!                   !          plasticStrainLoc = plasticStrainOldLoc
!                   !          Call SectionRealRestrict(plasticStrainPreviousSec,cellID(cell),plasticStrainPreviousLoc,ierr);CHKERRQ(ierr)
!                   !          plasticStrainPreviousLoc = plasticStrainLoc
!                   !          Call SectionRealRestore(plasticStrainPreviousSec,cellID(cell),plasticStrainPreviousLoc,ierr);CHKERRQ(ierr)
!                   !       Else
!                   !          exit_code = SNLPL1SQP(s,plasticStrainLoc)
!                   !       End if
!                   Case default
!                      if (cellSetOptions%plasticityType /= MEF90DefMech_plasticityTypeNONE) then
!                         exit_code = SNLPL1SQP(s,plasticStrainLoc)
!                      end if
!                   End Select

!                   !!! cumulatedDissipatedPlasticEnergy
!                   PlasticStrainMatS = plasticStrainLoc
!                   stiffness = 0.0_Kr
!                   Select Case (cellSetOptions%damageType)
!                      Case (MEF90DefMech_damageTypeAT1Elastic,MEF90DefMech_damageTypeAT2Elastic)
!                         Stiffness = 1.0_Kr
!                      Case (MEF90DefMech_damageTypeAT1,MEF90DefMech_damageTypeAT2)
!                         Stiffness = (1.0_Kr - PlasticityCtx%Damage)**(2.0_Kr-PlasticityCtx%DuctileCouplingPower) + PlasticityCtx%residualStiffness
!                      Case (MEF90DefMech_damageTypeLinSoft)
!                         Stiffness = ( (1.0_Kr - PlasticityCtx%Damage)**(2.0_Kr - PlasticityCtx%DuctileCouplingPower) / ( 1.0_Kr + ( PlasticityCtx%CoefficientLinSoft - 1.0_Kr )*(1.0_Kr - (1.0_Kr - PlasticityCtx%Damage)**2.0_kr ) ) ) + PlasticityCtx%residualStiffness
!                      Case default
!                         Print*,__FUNCT__,': Unimplemented damage Type, only AT1Elastic and AT2Elastic implement',cellSetOptions%damageType
!                         STOP
!                   End Select
!                   if (.NOT. PlasticityCtx%isViscousPlasticity) then
!                      cumulatedDissipatedPlasticEnergyVariationLoc(1) = Stiffness * ( PlasticityCtx%HookesLaw * ( PlasticityCtx%totalStrain - PlasticStrainMatS ) ) .dotP. ( PlasticStrainMatS - PlasticityCtx%plasticStrainOld )
!                   else
!                      cumulatedDissipatedPlasticEnergyVariationLoc(1) = Stiffness * PlasticityCtx%viscousCumulatedDissipatedPlasticEnergyVariation
!                   endif

! #if MEF90_DIM == 2
!                   if (PlasticityCtx%isPlaneStress .eqv. .FALSE.) then
!                      Sigma_33_PlaneStrain = (PlasticityCtx%HookesLaw%YoungsModulus - 2.0_Kr*PlasticityCtx%HookesLaw%PoissonRatio*PlasticityCtx%HookesLaw%mu)*trace(PlasticStrainMatS) + PlasticityCtx%HookesLaw%lambda*trace(PlasticityCtx%totalStrain)
!                      cumulatedDissipatedPlasticEnergyVariationLoc(1) = cumulatedDissipatedPlasticEnergyVariationLoc(1) + Stiffness * ( PlasticityCtx%HookesLaw%lambda *trace(PlasticStrainMatS) - Sigma_33_PlaneStrain ) * trace( PlasticStrainMatS - PlasticityCtx%plasticStrainOld )
!                   endif
! #endif

!                   Call SectionRealRestore(cumulatedDissipatedPlasticEnergyVariationSec,cellID(cell),cumulatedDissipatedPlasticEnergyVariationLoc,ierr);CHKERRQ(ierr)
!                   Call SectionRealRestore(cumulatedDissipatedPlasticEnergyOldSec,cellID(cell),cumulatedDissipatedPlasticEnergyOldLoc,ierr);CHKERRQ(ierr)

!                   Call SectionRealRestore(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
!                   Call SectionRealRestore(plasticStrainOldSec,cellID(cell),plasticStrainOldLoc,ierr);CHKERRQ(ierr)
!                   Call SectionRealRestore(totalStrainSec,cellID(cell),totalStrainLoc,ierr);CHKERRQ(ierr)
!                End Do !cell
!                Call MEF90Element_Destroy(elemDisplacement,ierr)
!                Call MEF90Element_Destroy(elemScal,ierr)
!                Call SNLPDelete(s)
!                DeAllocate(damageLoc)
!             End If ! set
!             Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
!             Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
!          End If !damageType /= NONE
!       End Do !! set
!       Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
!       Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)

!       !!! forward data plasticStrain & cumulatedDissipatedPlasticEnergy
!       Call SectionRealToVec(plasticStrainSec,MEF90DefMechCtx%cellDMMatSScatter,SCATTER_FORWARD,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)
!       Call SectionRealToVec(cumulatedDissipatedPlasticEnergyVariationSec,MEF90DefMechCtx%cellDMScalScatter,SCATTER_FORWARD,cumulatedDissipatedPlasticEnergyVariation,ierr);CHKERRQ(ierr)

!       Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
!       Call SectionRealDestroy(plasticStrainOldSec,ierr);CHKERRQ(ierr)
!       Call SectionRealDestroy(totalStrainSec,ierr);CHKERRQ(ierr)
!       Call SectionRealDestroy(plasticStrainPreviousSec,ierr);CHKERRQ(ierr)
!       Call SectionRealDestroy(cumulatedDissipatedPlasticEnergyVariationSec,ierr);CHKERRQ(ierr)
!       Call SectionRealDestroy(cumulatedDissipatedPlasticEnergyOldSec,ierr);CHKERRQ(ierr)
!       Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
!       If (Associated(MEF90DefMechCtx%damage)) Then
!          Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)
!       End If

!       If (Associated(MEF90DefMechCtx%temperature)) Then
!          Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
!       End If
! #else
!       write(*,*) 'This example needs SNLP'

!       ! Indicate a generic error:
!       ierr = 1
! #endif
   End Subroutine MEF90DefMechPlasticStrainUpdate
End Module MEF90_APPEND(m_MEF90_DefMechPlasticity,MEF90_DIM)D
