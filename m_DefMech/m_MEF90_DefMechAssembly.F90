#include "../MEF90/mef90.inc"
#include "../m_DefMech/mef90DefMech.inc"
Module MEF90_APPEND(m_MEF90_DefMechAssembly,MEF90_DIM)D
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
#define MEF90_HDRegularization 0.01_Kr

   Use m_MEF90
   Use m_MEF90_DefMechCtx
   Use MEF90_APPEND(m_MEF90_DefMechSplit,MEF90_DIM)D
   Use m_MEF90_DefMechAT
   
   Implicit none
   Private
   Public MEF90DefMechOperatorDisplacement,     &
          MEF90DefMechBilinearFormDisplacement, &
          MEF90DefMechWork,                     &
          MEF90DefMechCohesiveEnergy,           &
          MEF90DefMechPlasticDissipation,       &
          MEF90DefMechElasticEnergy,            &
          MEF90DefMechOperatorDamage,           &
          MEF90DefMechBilinearFormDamage,       &
          MEF90DefMechSurfaceEnergy,            &
          MEF90DefMechCrackVolume,              &
          MEF90DefMechCrackPressureRescaling,   &
          MEF90DefMechStress

   PetscReal,Parameter,Private                         :: cwLinSoft = 0.7853981634_Kr
Contains

! #undef __FUNCT__
! #define __FUNCT__ "MEF90DefMechOperatorDisplacementLinSoftLoc"
! !!!
! !!!  
! !!!  MEF90DefMechOperatorDisplacementLinSoftLoc:
! !!!  
! !!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
! !!!

!    Subroutine MEF90DefMechOperatorDisplacementLinSoftLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
!       PetscReal,Dimension(:),Pointer                     :: residualLoc
!       PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
!       Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
!       PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
!       Type(MEF90_MATPROP),Intent(IN)                     :: matprop
!       Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
!       Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
!       PetscReal,Intent(IN)                               :: CrackPressureCell

!       Type(MEF90_MATS)                                   :: sigma
!       Type(MEF90_VECT)                                   :: UU0
!       PetscReal                                          :: stiffness,temperature
!       PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
!       PetscLogDouble                                     :: flops
!       PetscErrorCode                                     :: ierr
!       PetscReal                                          :: k

!       numDofDisplacement = size(elemDisplacement%BF,1)
!       numDofDamage = size(elemDamage%BF,1)
!       numGauss = size(elemDisplacement%BF,2)
!       k = matprop%CoefficientLinSoft

!       Do iGauss = 1,numGauss
!          stiffness = 0.0_Kr
!          Do iDoF1 = 1,numDofDamage
!             stiffness = stiffness + damageDof(iDoF1) * elemDamage%BF(iDoF1,iGauss)
!          End Do
!          stiffness = (1.0_Kr - stiffness)**2 /( 1.0_Kr + ( k - 1.0_Kr )*(1.0_Kr - (1.0_Kr - stiffness)**2 ) )+ matProp%residualStiffness

!          sigma = 0.0_Kr
!          Do iDoF1 = 1,numDofDisplacement
!             sigma = sigma + xDof(iDoF1) * elemDisplacement%GradS_BF(iDoF1,iGauss)
!          End Do

!          temperature = 0.0_Kr
!          If (Associated(temperatureDof)) Then
!             Do iDoF1 = 1,numDofDamage
!                temperature = temperature + temperatureDoF(iDoF1) * elemDamage%BF(iDoF1,iGauss)
!             End Do
!          End If
!          sigma = stiffness * (matProp%HookesLaw * (sigma - temperature * matProp%LinearThermalExpansion - plasticStrainCell))

!          Do iDoF2 = 1,numDofDisplacement
!             residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement%Gauss_C(iGauss) * (sigma .DotP. elemDisplacement%GradS_BF(iDoF2,iGauss))
!          End Do

!          UU0 = 0.0_Kr
!          If (Associated(boundaryDisplacementDof) .AND. (matprop%cohesiveStiffness /= 0.0_Kr)) Then
!             Do iDoF1 = 1,numDofDisplacement
!                UU0 = UU0 + (xDof(iDoF1) - boundaryDisplacementDof(iDoF1)) * elemDisplacement%BF(iDoF1,iGauss)
!             End Do
!             UU0 = UU0 * matprop%cohesiveStiffness

!             Do iDoF2 = 1,numDofDisplacement
!                residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement%Gauss_C(iGauss) * (UU0 .DotP. elemDisplacement%BF(iDoF2,iGauss))
!             End Do
!          End If
!       End Do
!       !flops = 2 * numGauss * numDofDisplacement**2
!       !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
!    End Subroutine MEF90DefMechOperatorDisplacementLinSoftLoc



#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDisplacement"
!!!
!!!  
!!!  MEF90DefMechOperatorDisplacement: Build the operator. When called in SNES, the solution time should always match the target time, 
!!!                                    so there is no need for interpolation of the forcees, external, and boundary values
!!!  
!!!  (c) 2012-20 Blaise Bourdin bourdin@lsu.edu, Erwan Tanne erwan.tanne@gmail.com
!!!

   Subroutine MEF90DefMechOperatorDisplacement(snesDisplacement,displacement,residual,MEF90DefMechCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesDisplacement
      Type(Vec),Intent(IN)                               :: displacement
      Type(Vec),Intent(INOUT)                            :: residual
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(SectionReal)                                  :: displacementSec,displacementPreviousStepSec,residualSec
      Type(SectionReal)                                  :: plasticStrainSec,temperatureSec,damageSec
      Type(SectionReal)                                  :: boundaryDisplacementSec,forceSec,pressureForceSec,CrackPressureSec
      PetscReal,Dimension(:),Pointer                     :: displacementPtr,residualPtr
      PetscReal,Dimension(:),Pointer                     :: boundaryDisplacementPtr
      PetscReal,Dimension(:),Pointer                     :: boundaryDisplacementDof,displacementDof,displacementPreviousStepDof,damageDof,temperatureDof,residualLoc
      PetscReal,Dimension(:),Pointer                     :: pressureForceLoc,CrackPressureLoc
      PetscReal                                          :: pressureForceCell,CrackPressureCell
      PetscReal,Dimension(:),Pointer                     :: forceLoc
      Type(MEF90_VECT)                                   :: forceCell
      PetscReal,Dimension(:),Pointer                     :: plasticStrainLoc
      Type(MEF90_MATS)                                   :: plasticStrainCell
      Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof,bcIS
      PetscInt,dimension(:),Pointer                      :: setID,cellID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder,cell
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90DefMechVertexSetOptions_Type),Pointer    :: vertexSetOptions
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemDamage
      Type(MEF90Element_Type)                            :: elemDisplacementType,elemDamageType
      PetscInt                                           :: p,dof,c,dim,numDof
      PetscInt                                           :: nVal
      PetscReal,Dimension(:),Pointer                     :: nullPtr
      
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: globalOptions
      Type(Vec)                                          :: damageOld
      PetscReal                                          :: damageMin,damageMax,damageMaxChange
      Character(len=MEF90_MXSTRLEN)                      :: IOBuffer
      Class(MEF90_DefMechAT_Type),Allocatable            :: ATModel
      Class(MEF90_DEFMECHSPLIT),Allocatable              :: Split
      PetscBool                                          :: cellIsElastic,cellHasForce
      Type(MEF90_MATS)                                   :: inelasticStrainGauss,stressGaussPlus,stressGaussMinus,stressGauss
      Type(MEF90_VECT)                                   :: U0Gauss,gradientDamageGauss,displacementDampingGauss
      PetscReal                                          :: damageGauss,temperatureGauss
      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
    
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)
      If ( (GlobalOptions%solverType == MEF90DefMech_SolverTypeQuasiNewton1) .OR. (GlobalOptions%solverType == MEF90DefMech_SolverTypeQuasiNewton2)) Then
         Call VecDuplicate(MEF90DefMechCtx%damage,damageOld,ierr);CHKERRQ(ierr)
         Call VecCopy(MEF90DefMechCtx%damage,damageOld,ierr);CHKERRQ(ierr)
         Call SNESSolve(MEF90DefMechCtx%SNESdamage,PETSC_NULL_OBJECT,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)
         Call VecMin(MEF90DefMechCtx%damage,PETSC_NULL_INTEGER,damageMin,ierr);CHKERRQ(ierr)
         Call VecMax(MEF90DefMechCtx%damage,PETSC_NULL_INTEGER,damageMax,ierr);CHKERRQ(ierr)
         Call VecAxPy(damageOld,-1.0_Kr,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)
         Call VecNorm(damageOld,NORM_INFINITY,damageMaxChange,ierr);CHKERRQ(ierr)
         Call VecDestroy(damageOld,ierr);CHKERRQ(ierr)
       End If

      !!! Get Section and scatter associated with each field
      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,displacementSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(displacementSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,displacement,ierr);CHKERRQ(ierr) 
      If (globalOptions%dampingCoefficientDisplacement * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,displacementpreviousStepSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(displacementpreviousStepSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,MEF90DefMechCtx%displacementpreviousStep,ierr);CHKERRQ(ierr) 
      EndIf

      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,boundaryDisplacementSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(boundaryDisplacementSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,MEF90DefMechCtx%boundaryDisplacement,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,residualSec,ierr);CHKERRQ(ierr)


      Call SectionRealSet(residualSec,0.0_Kr,ierr);CHKERRQ(ierr)
      Call VecSet(residual,0.0_Kr,ierr);CHKERRQ(ierr)

      Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,damageSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(damageSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)   

      Nullify(TemperatureDof)
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,temperatureSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(temperatureSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr)   
      Else
         temperatureSec%v = 0
      End If

      Nullify(ForceLoc)
      If (Associated(MEF90DefMechCtx%force)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%cellDMVectSec,forceSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(forceSec,MEF90DefMechCtx%cellDMVectScatter,SCATTER_REVERSE,MEF90DefMechCtx%force,ierr);CHKERRQ(ierr)
      Else
         ForceSec%v = 0
      End If
      
      Nullify(pressureForceLoc)
      If (Associated(MEF90DefMechCtx%pressureForce)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%cellDMScalSec,pressureForceSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(pressureForceSec,MEF90DefMechCtx%CellDMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%pressureForce,ierr);CHKERRQ(ierr)
      Else
         pressureForceSec%v = 0
      End If

      Nullify(CrackPressureLoc)
      If (Associated(MEF90DefMechCtx%CrackPressure)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%cellDMScalSec,CrackPressureSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(CrackPressureSec,MEF90DefMechCtx%CellDMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%CrackPressure,ierr);CHKERRQ(ierr)
      Else
         CrackPressureSec%v = 0
      End If

      Nullify(plasticStrainLoc)
      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%CellDMMatSSec,plasticStrainSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticStrainSec,MEF90DefMechCtx%CellDMMatSScatter,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)
      Else
         PlasticStrainSec%v = 0
      End If

      !!! get IS for cell sets
      Call DMMeshGetDimension(MEF90DefMechCtx%DM,dim,ierr);CHKERRQ(ierr)
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

      !!!
      !!! We loop over all element twice. The first time in order to assemble all non BC cell sets
      !!! In the second pass, we only update the BC where necessary
      !!! vertex set BC are updated last, so that they override cell set BC
      !!!
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
         elemDamageType       = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)

         Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
         If (Size(cellID) > 0) Then
            !!! get the ATModel and split objects
            Call MEF90DefMechGetATModel(cellSetOptions,ATModel,cellIsElastic)
            Call MEF90DefMechGetSplit(cellSetOptions,Split)

            !!! Allocate elements 
            QuadratureOrder = max(2*elemDisplacementType%order,split%damageOrder + split%strainOrder)
            Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
            Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemDamage,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)

            !!! Allocate storage for fields at dof and Gauss points
            Allocate(residualLoc(ElemDisplacementType%numDof))
            Allocate(displacementDof(ElemDisplacementType%numDof))
            Allocate(damageDof(ElemDamageType%numDof))
            If (globalOptions%dampingCoefficientDisplacement * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
               Allocate(displacementPreviousStepDof(ElemDisplacementType%numDof))
            Else
               Nullify(displacementpreviousStepDof)
            End If
            If (matpropSet%cohesiveStiffness /= 0.0_Kr) Then
               Allocate(boundaryDisplacementDof(ElemDisplacementType%numDof))
            Else
               Nullify(boundaryDisplacementDof)
            End If
            If (Associated(MEF90DefMechCtx%temperature)) Then
               Allocate(temperatureDof(ElemDamageType%numDof))
            Else
               Nullify(temperatureDof)
            End If
            Do cell = 1,size(cellID)
               !!! Get value of each field at each Dof of the local element
               Call SectionRealRestrictClosure(displacementSec,MEF90DefMechCtx%DMVect,cellID(cell),elemDisplacementType%numDof,displacementDof,ierr);CHKERRQ(ierr)
               Call SectionRealRestrictClosure(damageSec,MEF90DefMechCtx%DMScal,cellID(cell),elemDamageType%numDof,damageDof,ierr);CHKERRQ(ierr)
               If (matpropSet%cohesiveStiffness /= 0.0_Kr) Then
                  Call SectionRealRestrictClosure(boundaryDisplacementSec,MEF90DefMechCtx%DMVect,cellID(cell),elemDisplacementType%numDof,boundaryDisplacementDof,ierr);CHKERRQ(ierr)
               End If
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  Call SectionRealRestrictClosure(temperatureSec,MEF90DefMechCtx%DMScal,cellID(cell),elemDamageType%numDof,temperatureDof,ierr);CHKERRQ(ierr)
               End If

               forceCell = 0.0_Kr
               If (Associated(MEF90DefMechCtx%force)) Then
                  Call SectionRealRestrict(forceSec,cellID(cell),forceLoc,ierr);CHKERRQ(ierr)
                  forceCell = forceLoc
               End If
               If ((Associated(MEF90DefMechCtx%pressureForce)) .AND. (ElemDisplacementType%coDim == 0)) Then
                  Call SectionRealRestrict(pressureForceSec,cellID(cell),pressureForceLoc,ierr);CHKERRQ(ierr)
                  ForceCell = forceCell + pressureForceLoc(1) * elemDisplacement(cell)%outerNormal
               End If
               cellHasForce = .FALSE.
               If (norm(forceCell) > 0.0_Kr) Then
                  cellHasForce = .TRUE.
               End If

               If (Associated(MEF90DefMechCtx%CrackPressure)) Then
                  Call SectionRealRestrict(CrackPressureSec,cellID(cell),CrackPressureLoc,ierr);CHKERRQ(ierr)
                  CrackPressureCell = CrackPressureLoc(1)
               Else
                  CrackPressureCell = 0.0_Kr
               End If

               If (Associated(MEF90DefMechCtx%plasticStrain)) Then
                  Call SectionRealRestrict(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
                  plasticStrainCell = plasticStrainLoc
               Else
                  plasticStrainCell = 0.0_Kr
               End If
                  
               If (globalOptions%dampingCoefficientDisplacement * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
                  Call SectionRealRestrictClosure(displacementpreviousStepSec,MEF90DefMechCtx%DMVect,cellID(cell),elemDisplacementType%numDof,displacementPreviousStepDof,ierr);CHKERRQ(ierr)
               End If

               !!! Loop over Gauss point, evaluate fields and assemble the cell contribution to the residual
               numDofDisplacement = size(elemDisplacement(cell)%BF,1)
               numDofDamage = size(elemDamage(cell)%BF,1)
               numGauss = size(elemDisplacement(cell)%BF,2)

               residualLoc = 0.0_Kr
               Do iGauss = 1,numGauss

                  !!! Main term: [a(\alpha) sigma^+(u) + sigma^-(u)] . e(v) 
                  If (.NOT. cellIsElastic) Then
                     damageGauss = 0.0_Kr
                     Do iDof1 = 1,numDofDamage
                        damageGauss = damageGauss + damageDof(iDoF1) * elemDamage(cell)%BF(iDoF1,iGauss)
                     End Do ! iDof1
                  End If

                  inelasticStrainGauss = 0.0_Kr
                  Do iDof1 = 1,numDofDisplacement
                     inelasticStrainGauss = inelasticStrainGauss + displacementDof(iDoF1) * elemDisplacement(cell)%GradS_BF(iDoF1,iGauss)
                  End Do

                  temperatureGauss = 0.0_Kr
                  If (Associated(temperatureDof)) Then
                     Do iDoF1 = 1,numDofDamage
                        temperatureGauss = temperatureGauss + temperatureDoF(iDoF1) * elemDamage(cell)%BF(iDoF1,iGauss)
                     End Do ! iDof1
                     inelasticStrainGauss = inelasticStrainGauss - (temperatureGauss * matpropSet%LinearThermalExpansion)
                  End If ! Associated(temperatureDof)
#if MEF90_DIM == 2
!!! We need something along these lines
                  !!! Adding terms in planestrain for plasticity with tr(p) = 0
                  ! If (.NOT. matProp%HookesLaw%isPlaneStress) Then
                  !    stressGauss = stressGauss +  stiffness * ( matProp%HookesLaw%lambda*trace(plasticStrainCell)*MEF90MatS2DIdentity )
                  ! End If
#endif

                  Call Split%DEED(inelasticStrainGauss - plasticStrainCell,matpropSet%HookesLaw,stressGaussPlus,stressGaussMinus)
                  If (cellIsElastic) Then
                     stressGauss = stressGaussPlus + stressGaussMinus
                  Else
                     stressGauss = ATModel%a(damageGauss) * stressGaussPlus + stressGaussMinus
                  End If
                  Do iDoF2 = 1,numDofDisplacement
                     residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement(cell)%Gauss_C(iGauss) * (stressGauss .DotP. elemDisplacement(cell)%GradS_BF(iDoF2,iGauss))
                  End Do ! iDof2

                  If (ElemDisplacementType%coDim == 0) Then
                     !!! Damping
                     If (globalOptions%dampingCoefficientDisplacement * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
                        displacementDampingGauss = 0.0_Kr
                        Do iDof1 = 1,numDofDisplacement
                           displacementDampingGauss = displacementDampingGauss + ((displacementDof(iDoF1) - displacementPreviousStepDof(iDoF1)) * elemDisplacement(cell)%BF(iDoF1,iGauss))
                        End Do
                        displacementDampingGauss = displacementDampingGauss * globalOptions%dampingCoefficientDisplacement / MEF90DefMechCtx%timeStep
                        Do iDoF2 = 1,numDofDisplacement
                           residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement(cell)%Gauss_C(iGauss) * &
                                                ( displacementDampingGauss .DotP. elemDisplacement(cell)%BF(iDoF2,iGauss) )
                        End Do
                     End If ! damping

                     !!! Cohesive force
                     If (Associated(boundaryDisplacementDof) .AND. (matpropSet%cohesiveStiffness /= 0.0_Kr)) Then
                        U0Gauss = 0.0_Kr
                        Do iDof1 = 1,numDofDisplacement
                           U0Gauss = U0Gauss + (displacementDof(iDoF1) - boundaryDisplacementDof(iDoF1)) * elemDisplacement(cell)%BF(iDoF1,iGauss)
                        End Do ! iDof1
                        U0Gauss = U0Gauss * matpropSet%cohesiveStiffness
                        Do iDoF2 = 1,numDofDisplacement
                           residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement(cell)%Gauss_C(iGauss) * (U0Gauss .DotP. elemDisplacement(cell)%BF(iDoF2,iGauss))
                        End Do ! iDof2
                     End If ! cohesiveStiffness

                     !!! crack pressure
                     If (Associated(MEF90DefMechCtx%CrackPressure) .AND. (CrackPressureCell /= 0.0_Kr)) Then
                        gradientDamageGauss = 0.0_Kr
                        Do iDoF1 = 1,numDofDamage
                           GradientDamageGauss = GradientDamageGauss + damageDof(iDoF1) * elemDamage(cell)%Grad_BF(iDoF1,iGauss) 
                        End Do ! iDof1
                        Do iDoF2 = 1,numDofDisplacement
                           residualLoc(iDoF2) = residualLoc(iDoF2) + elemDisplacement(cell)%Gauss_C(iGauss) * &
                                                CrackPressureCell * (GradientDamageGauss  .DotP. elemDisplacement(cell)%BF(iDoF2,iGauss))
                        End Do ! iDof2
                     End If ! crack Pressure
                  End If ! coDim == 0

                  !!! force
                  If (cellHasForce) Then
                     Do iDoF2 = 1,numDofDisplacement
                        residualLoc(iDoF2) = residualLoc(iDoF2) - elemDisplacement(cell)%Gauss_C(iGauss) * &
                                           ( forceCell .DotP. elemDisplacement(cell)%BF(iDoF2,iGauss))
                     End Do ! iDof2 
                  End If ! cellHasForce

               End Do ! iGauss

               Call SectionRealUpdateClosure(residualSec,MEF90DefMechCtx%DMVect,cellID(cell),residualLoc,ADD_VALUES,ierr);CHKERRQ(iErr)

               If (Associated(MEF90DefMechCtx%force)) Then
                  Call SectionRealRestore(forceSec,cellID(cell),forceLoc,ierr);CHKERRQ(ierr)
               End If
               If ((Associated(MEF90DefMechCtx%pressureForce)) .AND. (ElemDisplacementType%coDim == 0)) Then
                  Call SectionRealRestore(pressureForceSec,cellID(cell),pressureForceLoc,ierr);CHKERRQ(ierr)
               End If
               If (Associated(MEF90DefMechCtx%CrackPressure)) Then
                  Call SectionRealRestore(CrackPressureSec,cellID(cell),CrackPressureLoc,ierr);CHKERRQ(ierr)
               End If
               If (Associated(MEF90DefMechCtx%plasticStrain)) Then
                  Call SectionRealRestore(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
               End If
            End Do ! cell

            DeAllocate(displacementDof)
            If (globalOptions%dampingCoefficientDisplacement * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
               DeAllocate(displacementPreviousStepDof)
            End If
            If (matpropSet%cohesiveStiffness /= 0.0_Kr) Then
               DeAllocate(boundaryDisplacementDof)
            End If
            DeAllocate(damageDof)
            If (Associated(MEF90DefMechCtx%temperature)) Then
               DeAllocate(temperatureDof)
            End If
            DeAllocate(residualLoc)

            Call MEF90Element_Destroy(elemDisplacement,ierr)
            Call MEF90Element_Destroy(elemDamage,ierr)
         End If 
         Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do ! set

      !!! "Ghost update" for the residual Section
      Call SectionRealComplete(residualSec,ierr);CHKERRQ(ierr)
      !!! Scatter back from SectionReal to Vec
      Call SectionRealToVec(residualSec,MEF90DefMechCtx%DMVectScatter,SCATTER_FORWARD,residual,ierr);CHKERRQ(ierr)

      !!!
      !!! Cell set BC
      !!!
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Do c = 1, dim
            If (cellSetOptions%Has_DisplacementBC(c)) Then
               Call DMMeshGetStratumIS(MEF90DefMechCtx%DMVect,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
               Call MEF90ISCreateCelltoVertex(MEF90DefMechCtx%DMVect,PETSC_COMM_WORLD,setIS,bcIS,ierr)
               Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
               Allocate(displacementPtr(nval),stat=ierr)
               Allocate(residualPtr(nval),stat=ierr)
               Allocate(boundaryDisplacementPtr(nval),stat=ierr)
               Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMVect,displacement,displacementPtr,bcIS,c,ierr)
               Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMVect,MEF90DefMechCtx%boundaryDisplacement,boundaryDisplacementPtr,bcIS,c,ierr)
               residualPtr = displacementPtr - boundaryDisplacementPtr
               Call MEF90VecSetValuesISdof(MEF90DefMechCtx%DMVect,residual,residualPtr,bcIS,c,INSERT_VALUES,ierr)
               DeAllocate(boundaryDisplacementPtr)
               DeAllocate(residualPtr)
               DeAllocate(displacementPtr)
               Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
               Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
            End If ! cellSetOptions%Has_BC
         End Do
      End Do ! set
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)

      !!!
      !!! Vertex set BC
      !!!
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxVertexSetOptions(MEF90DefMechCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
         Do c = 1, dim
            If (vertexSetOptions%Has_DisplacementBC(c)) Then
               Call DMMeshGetStratumIS(MEF90DefMechCtx%DMVect,'Vertex Sets',setID(set),bcIS,ierr);CHKERRQ(iErr)
               Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
               Allocate(displacementPtr(nval),stat=ierr)
               Allocate(residualPtr(nval),stat=ierr)
               Allocate(boundaryDisplacementPtr(nval),stat=ierr)
               Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMVect,displacement,displacementPtr,bcIS,c,ierr)
               Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMVect,MEF90DefMechCtx%boundaryDisplacement,boundaryDisplacementPtr,bcIS,c,ierr)
               residualPtr = displacementPtr - boundaryDisplacementPtr
               Call MEF90VecSetValuesISdof(MEF90DefMechCtx%DMVect,residual,residualPtr,bcIS,c,INSERT_VALUES,ierr)
               DeAllocate(boundaryDisplacementPtr)
               DeAllocate(residualPtr)
               DeAllocate(displacementPtr)
               Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
            End If ! vertexSetOptions%Has_BC
         End Do
      End Do ! set
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
   
      Call VecAssemblyBegin(residual,ierr)
      Call VecAssemblyEnd(residual,ierr)

      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
      End If
      
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
      End If

      If (Associated(MEF90DefMechCtx%pressureForce)) Then
         Call SectionRealDestroy(pressureForceSec,ierr);CHKERRQ(ierr)
      End If

      If (Associated(MEF90DefMechCtx%CrackPressure)) Then
         Call SectionRealDestroy(CrackPressureSec,ierr);CHKERRQ(ierr)
      End If

      If (Associated(MEF90DefMechCtx%force)) Then
         Call SectionRealDestroy(forceSec,ierr);CHKERRQ(ierr)
      End If

      Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(boundaryDisplacementSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(residualSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(displacementSec,ierr);CHKERRQ(ierr)
      If (globalOptions%dampingCoefficientDisplacement * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
         Call SectionRealDestroy(displacementPreviousStepSec,ierr);CHKERRQ(ierr)
      End If
   End Subroutine MEF90DefMechOperatorDisplacement

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDisplacement"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDisplacement:
!!!  
!!!  (c) 2012-20 Blaise Bourdin bourdin@lsu.edu,Erwan Tanne erwan.tanne@gmail.com
!!!

   Subroutine MEF90DefMechBilinearFormDisplacement(snesDisplacement,displacement,A,M,flg,MEF90DefMechCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesDisplacement
      Type(Vec),Intent(IN)                               :: displacement
      Type(Mat),Intent(INOUT)                            :: A,M
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr  
         
      Type(SectionReal)                                  :: defaultSec
      Type(SectionReal)                                  :: displacementSec,damageSec,temperatureSec,plasticStrainSec
      PetscReal,Dimension(:),Pointer                     :: displacementDOF,damageDof,temperatureDof,nullPtr
      PetscReal,Dimension(:),Pointer                     :: plasticStrainLoc
      Type(MEF90_MATS)                                   :: plasticStrainCell
      Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90DefMechVertexSetOptions_Type),Pointer    :: vertexSetOptions
      Type(MEF90Element_Type)                            :: ElemDisplacementType,ElemDamageType
      PetscInt,Dimension(:),Pointer                      :: setIdx,bcIdx,Cone
      Type(IS)                                           :: bcIS
      PetscInt                                           :: cell,v,numBC,numDoF,numCell,c,dim
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemDamage
      
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: globalOptions
      Type(Vec)                                          :: damageOld
      PetscReal                                          :: damageMin,damageMax,damageMaxChange,cohesiveStiffness
      Character(len=MEF90_MXSTRLEN)                      :: IOBuffer
      Class(MEF90_DefMechAT_Type),Allocatable            :: ATModel
      Class(MEF90_DEFMECHSPLIT),Allocatable              :: Split
      PetscBool                                          :: cellIsElastic
      PetscReal                                          :: cellDampingCoefficient
      Type(MEF90_HOOKESLAW)                              :: AGaussPlus,AGaussMinus,AGauss
      Type(MEF90_MATS)                                   :: inelasticStrainGauss,AGradS_BF1
      PetscReal                                          :: damageGauss,temperatureGauss
      Type(MEF90_VECT)                                   :: BF1
    
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)
      If (GlobalOptions%solverType == MEF90DefMech_SolverTypeQuasiNewton2) Then
         Call VecDuplicate(MEF90DefMechCtx%damage,damageOld,ierr);CHKERRQ(ierr)
         Call VecCopy(MEF90DefMechCtx%damage,damageOld,ierr);CHKERRQ(ierr)
         Call SNESSolve(MEF90DefMechCtx%SNESdamage,PETSC_NULL_OBJECT,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)
         Call VecMin(MEF90DefMechCtx%damage,PETSC_NULL_INTEGER,damageMin,ierr);CHKERRQ(ierr)
         Call VecMax(MEF90DefMechCtx%damage,PETSC_NULL_INTEGER,damageMax,ierr);CHKERRQ(ierr)
         Call VecAxPy(damageOld,-1.0_Kr,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)
         Call VecNorm(damageOld,NORM_INFINITY,damageMaxChange,ierr);CHKERRQ(ierr)
         Call VecDestroy(damageOld,ierr);CHKERRQ(ierr)
       End If
      
      Call MatZeroEntries(A,ierr);CHKERRQ(ierr)
      Call SNESGetDM(snesDisplacement,MEF90DefMechCtx%DM,ierr);CHKERRQ(ierr)
      
      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,displacementSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(displacementSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,displacement,ierr);CHKERRQ(ierr) 

      Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,damageSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(damageSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)   

      Nullify(TemperatureDof)
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,temperatureSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(temperatureSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr)   
      Else
         temperatureSec%v = 0
      End If

      Nullify(plasticStrainLoc)
      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%CellDMMatSSec,plasticStrainSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticStrainSec,MEF90DefMechCtx%CellDMMatSScatter,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)
      Else
         PlasticStrainSec%v = 0
      End If

      !!! get IS for cell sets
      Call DMMeshGetDimension(MEF90DefMechCtx%DM,dim,ierr);CHKERRQ(ierr)
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID) 
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         ElemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
         ElemDamageType       = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)
         !!! Go through cells and call local assembly functions
         Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
         If ((Size(cellID) > 0) .AND. (ElemDisplacementType%coDim == 0)) Then
            !!! get the ATModel and split objects
            Call MEF90DefMechGetATModel(cellSetOptions,ATModel,cellIsElastic)
            Call MEF90DefMechGetSplit(cellSetOptions,Split)
            QuadratureOrder = split%damageOrder + split%strainOrder

            !!! Allocate elements 
            Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
            Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemDamage,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)

            !!! Allocate storage for fields at dof and Gauss points
            Allocate(Aloc(ElemDisplacementType%numDof,ElemDisplacementType%numDof))
            Allocate(displacementDof(ElemDisplacementType%numDof))
            Allocate(damageDof(ElemDamageType%numDof))
            If (Associated(MEF90DefMechCtx%temperature)) Then
               Allocate(temperatureDof(ElemDamageType%numDof))
            End If

            cellDampingCoefficient = globalOptions%dampingCoefficientDisplacement * MEF90DefMechCtx%timeStep + matpropSet%cohesiveStiffness
            
            Do cell = 1,size(cellID)
               !!! Get value of each field at each Dof of the local element
               Call SectionRealRestrictClosure(displacementSec,MEF90DefMechCtx%DMVect,cellID(cell),elemDisplacementType%numDof,displacementDof,ierr);CHKERRQ(ierr)
               Call SectionRealRestrictClosure(damageSec,MEF90DefMechCtx%DMScal,cellID(cell),elemDamageType%numDof,damageDof,ierr);CHKERRQ(ierr)
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  Call SectionRealRestrictClosure(temperatureSec,MEF90DefMechCtx%DMScal,cellID(cell),elemDamageType%numDof,temperatureDof,ierr);CHKERRQ(ierr)
               End If
               If (Associated(MEF90DefMechCtx%plasticStrain)) Then
                  Call SectionRealRestrict(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
                  plasticStrainCell = plasticStrainLoc
               Else
                  plasticStrainCell = 0.0_Kr
               End If
                  
               !!! Loop over Gauss point, evaluate fields and assemble the cell contribution to the bilinear form
               numDofDisplacement = size(elemDisplacement(cell)%BF,1)
               numDofDamage = size(elemDamage(cell)%BF,1)
               numGauss = size(elemDisplacement(cell)%BF,2)

               Aloc = 0.0_Kr
               Do iGauss = 1,numGauss
                  If (.NOT. cellIsElastic) Then
                     damageGauss = 0.0_Kr
                     Do iDof1 = 1,numDofDamage
                        damageGauss = damageGauss + damageDof(iDoF1) * elemDamage(cell)%BF(iDoF1,iGauss)
                     End Do ! iDof1
                  End If ! cellIsElastic

                  inelasticStrainGauss = 0.0_Kr
                  Do iDof1 = 1,numDofDisplacement
                     inelasticStrainGauss = inelasticStrainGauss + displacementDof(iDoF1) * elemDisplacement(cell)%GradS_BF(iDoF1,iGauss)
                  End Do ! iDoF1

                  temperatureGauss = 0.0_Kr
                  If (Associated(temperatureDof)) Then
                     Do iDoF1 = 1,numDofDamage
                        temperatureGauss = temperatureGauss + temperatureDoF(iDoF1) * elemDamage(cell)%BF(iDoF1,iGauss)
                     End Do ! iDof1
                     inelasticStrainGauss = inelasticStrainGauss - (temperatureGauss * matpropSet%LinearThermalExpansion)
                  End If ! Associated(temperatureDof)
#if MEF90_DIM == 2
!!! We need something along these lines
               !!! Adding terms in planestrain for plasticity with tr(p) = 0
               ! If (.NOT. matProp%HookesLaw%isPlaneStress) Then
               !    stressGauss = stressGauss +  stiffness * ( matProp%HookesLaw%lambda*trace(plasticStrainCell)*MEF90MatS2DIdentity )
               ! End If
#endif
                  If (cellIsElastic) Then
                     AGauss = matpropSet%HookesLaw
                  Else
                     Call Split%D2EED(inelasticStrainGauss - plasticStrainCell,matpropSet%HookesLaw,AGaussPlus,AGaussMinus)
                     AGauss = ATModel%a(damageGauss) * AGaussPlus + AGaussMinus
                  End If ! cellIsElastic
                  Do iDoF1 = 1,numDofDisplacement
                     AGradS_BF1 = AGauss * elemDisplacement(cell)%GradS_BF(iDoF1,iGauss)
                     Do iDoF2 = 1,numDofDisplacement
                        ALoc(iDoF2,iDoF1) = ALoc(iDoF2,iDoF1) + elemDisplacement(cell)%Gauss_C(iGauss) * (AGradS_BF1 .DotP. elemDisplacement(cell)%GradS_BF(iDoF2,iGauss))
                     End Do ! iDoF2
                  End Do ! iDoF1
                  If (cellDampingCoefficient /= 0.0_Kr) Then
                     Do iDoF1 = 1,numDofDisplacement
                        BF1 = cellDampingCoefficient * elemDisplacement(cell)%BF(iDoF1,iGauss)
                        Do iDoF2 = 1,numDofDisplacement
                              ALoc(iDoF2,iDoF1) = ALoc(iDoF2,iDoF1) + elemDisplacement(cell)%Gauss_C(iGauss) * (BF1 .DotP. elemDisplacement(cell)%BF(iDoF2,iGauss))
                        End Do ! iDoF1
                     End Do ! iDoF2
                  End If ! cellDampingCoefficient
               End Do ! iGauss
               Call DMmeshAssembleMatrix(A,MEF90DefMechCtx%DM,displacementSec,cellID(cell),ALoc,ADD_VALUES,ierr);CHKERRQ(ierr)
               If (Associated(MEF90DefMechCtx%plasticStrain)) Then
                  Call SectionRealRestore(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
               End If
            End Do ! cell

            DeAllocate(displacementDof)
            DeAllocate(damageDof)
            If (Associated(MEF90DefMechCtx%temperature)) Then
               DeAllocate(temperatureDof)
            End If
            Call MEF90Element_Destroy(elemDisplacement,ierr)
            Call MEF90Element_Destroy(elemDamage,ierr)
            DeAllocate(Aloc)    
         End If ! coDim

         Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do ! set
      Call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
      Call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
      !!!
      !!! Boundary conditions at cell sets
      !!!
      Call DMmeshGetStratumSize(MEF90DefMechCtx%DM,"height",0,numCell,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call MEF90ISCreateCelltoVertex(MEF90DefMechCtx%DM,PETSC_COMM_WORLD,setIS,bcIS,ierr)
         Do c = 1, dim
            If (cellSetOptions%Has_displacementBC(c)) Then
               Call DMMeshISCreateISglobaldof(MEF90DefMechCtx%DM,bcIS,c-1,setISdof,ierr);CHKERRQ(ierr)
               Call MatZeroRowsColumnsIS(A,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
               Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
            End If ! cellSetOptions%Has_displacementBC
         End Do ! c
         Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do     
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
   
      !!!
      !!! Boundary conditions at vertex sets
      !!!
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DM,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxVertexSetOptions(MEF90DefMechCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
         Do c = 1, dim
            If (vertexSetOptions%Has_displacementBC(c)) Then
               Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
               Call DMMeshISCreateISglobaldof(MEF90DefMechCtx%DM,setIS,c-1,setISdof,ierr);CHKERRQ(ierr)
               Call MatZeroRowsColumnsIS(A,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
               Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
               Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
            End If
         End Do
      End Do
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      
      Call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
      Call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)

      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
      End If
      
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
      End If

      Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(displacementSec,ierr);CHKERRQ(ierr)
      flg = SAME_NONZERO_PATTERN
   End Subroutine MEF90DefMechBilinearFormDisplacement

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechWork"
!!!
!!!  
!!!  MEF90DefMechWork:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechWork(xVec,MEF90DefMechCtx,work,ierr)
      Type(Vec),Intent(IN)                               :: xVec
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,Dimension(:),Pointer                     :: work
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(SectionReal)                                  :: xSec
      Type(SectionReal)                                  :: forceSec,pressureForceSec
      Type(IS)                                           :: CellSetGlobalIS,setIS
      PetscInt,Dimension(:),Pointer                      :: setID
      PetscInt                                           :: set,QuadratureOrder,numCell
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90Element_Type)                            :: elemDisplacementType
      Type(MEF90DefMechGlobalOptions_Type),Pointer       :: globalOptions
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      PetscReal                                          :: myWork,myWorkSet
      
      
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%globalOptionsBag,globalOptions,ierr);CHKERRQ(ierr)

      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,xSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(xSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,xVec,ierr);CHKERRQ(ierr)

      Call SectionRealDuplicate(MEF90DefMechCtx%cellDMVectSec,forceSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(forceSec,MEF90DefMechCtx%cellDMVectScatter,SCATTER_REVERSE,MEF90DefMechCtx%force,ierr);CHKERRQ(ierr)

      Call SectionRealDuplicate(MEF90DefMechCtx%cellDMScalSec,pressureForceSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(pressureForceSec,MEF90DefMechCtx%CellDMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%pressureForce,ierr);CHKERRQ(ierr)

      Work   = 0.0_Kr
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         mywork = 0.0_Kr
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)

         !!! Call proper local assembly depending on the type of damage law
         QuadratureOrder = elemDisplacementType%order
         Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)

         !!! Force work
         If (globalOptions%forceScaling /= MEF90Scaling_Null) Then
            Call MEF90ElasticityWorkSetCell(mywork,xSec,MEF90DefMechCtx%CellDMVect,ForceSec,setIS,elemDisplacement,elemDisplacementType,ierr)
         End If

         !!! pressure force work
         If ((globalOptions%pressureForceScaling /= MEF90Scaling_Null) .AND. (elemDisplacementType%coDim > 0)) Then
            Call MEF90ElasticityPressureWorkSetCell(mywork,xSec,MEF90DefMechCtx%CellDMVect,pressureForceSec,setIS,elemDisplacement,elemDisplacementType,ierr)
         End If
         
         Call MEF90Element_Destroy(elemDisplacement,ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)

         Call MPI_AllReduce(myWork,myWorkSet,1,MPIU_SCALAR,MPI_SUM,MEF90DefMechCtx%MEF90Ctx%comm,ierr);CHKERRQ(ierr)
         work(set) = work(set) + myWorkSet
      End Do ! set
      Call ISrestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr) 

      Call SectionRealDestroy(pressureForceSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(forceSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechWork

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCohesiveEnergy"
!!!
!!!  
!!!  MEF90DefMechCohesiveEnergy:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechCohesiveEnergy(xVec,MEF90DefMechCtx,cohesiveEnergy,ierr)
      Type(Vec),Intent(IN)                               :: xVec
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,Dimension(:),Pointer                     :: cohesiveEnergy
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(SectionReal)                                  :: xSec
      Type(SectionReal)                                  :: boundaryDisplacementSec
      Type(IS)                                           :: CellSetGlobalIS,setIS
      PetscInt,Dimension(:),Pointer                      :: setID
      PetscInt                                           :: set,QuadratureOrder,numCell
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90Element_Type)                            :: elemDisplacementType
      Type(MEF90DefMechGlobalOptions_Type),Pointer       :: globalOptions
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      PetscReal                                          :: mycohesiveEnergy,mycohesiveEnergySet
      
      
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%globalOptionsBag,globalOptions,ierr);CHKERRQ(ierr)
      !!! Create dof-based sections
      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,xSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(xSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,MEF90DefMechCtx%boundaryDisplacement,ierr);CHKERRQ(ierr)

      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,boundaryDisplacementSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(boundaryDisplacementSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,xVec,ierr);CHKERRQ(ierr)

      cohesiveEnergy   = 0.0_Kr
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         mycohesiveEnergy = 0.0_Kr
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)

         !!! Call proper local assembly depending on the type of damage law
         QuadratureOrder = elemDisplacementType%order
         Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)

         Call MEF90ElasticityCohesiveEnergySet(mycohesiveEnergy,xSec,MEF90DefMechCtx%CellDMVect,boundaryDisplacementSec,setIS,elemDisplacement,elemDisplacementType,ierr)

         Call MEF90Element_Destroy(elemDisplacement,ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)

         Call MPI_AllReduce(mycohesiveEnergy,mycohesiveEnergySet,1,MPIU_SCALAR,MPI_SUM,MEF90DefMechCtx%MEF90Ctx%comm,ierr);CHKERRQ(ierr)
         cohesiveEnergy(set) = cohesiveEnergy(set) + mycohesiveEnergySet * matpropSet%cohesiveStiffness
      End Do ! set
      Call ISrestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr) 

      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(boundaryDisplacementSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechCohesiveEnergy


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechPlasticDissipation"
!!!
!!!  
!!!  MEF90DefMechPlasticDissipation: 
!!!  
!!!  (c) 2015 Erwan TANNE erwan.tanne@gmail.com
!!!

   Subroutine MEF90DefMechPlasticDissipation(x,MEF90DefMechCtx,plasticStrainOld,energy,ierr)
      Type(Vec),Intent(IN)                               :: x
      Type(Vec),Intent(IN)                               :: plasticStrainOld
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,dimension(:),Pointer                     :: energy
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(SectionReal)                                  :: xSec
      Type(SectionReal)                                  :: damageSec,plasticStrainSec,plasticStrainOldSec,temperatureSec,cumulatedDissipatedPlasticEnergySec
      Type(IS)                                           :: CellSetGlobalIS,setIS,setISdof,bcIS
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemScal
      Type(MEF90Element_Type)                            :: elemDisplacementType,elemScalType
      PetscReal                                          :: myenergy,myenergySet

      !!! Get Section and scatter associated with each field
      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,xSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(xSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr) 

      Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,damageSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(damageSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)   

      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,temperatureSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(temperatureSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr)   
      Else
         temperatureSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%CellDMMatSSec,plasticStrainSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticStrainSec,MEF90DefMechCtx%CellDMMatSScatter,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)
         
         Call SectionRealDuplicate(MEF90DefMechCtx%CellDMScalSec,cumulatedDissipatedPlasticEnergySec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(cumulatedDissipatedPlasticEnergySec,MEF90DefMechCtx%CellDMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%cumulatedDissipatedPlasticEnergy,ierr);CHKERRQ(ierr)
      Else
         PlasticStrainSec%v = 0
         cumulatedDissipatedPlasticEnergySec%v = 0
      End If

      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)      
      Do set = 1,size(setID)
         myenergy = 0.0_Kr
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DMVect,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
         elemScalType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)


         !!! Call proper local assembly depending on the type of damage law
         Select Case (cellSetOptions%damageType)
         Case (MEF90DefMech_damageTypeAT1Elastic,MEF90DefMech_damageTypeAT2Elastic,MEF90DefMech_damageTypeAT1,MEF90DefMech_damageTypeLinSoft)
            If (elemDisplacementType%coDim == 0) Then
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  QuadratureOrder = 2 * max(elemDisplacementType%order - 1, elemScalType%order)
               Else
                  QuadratureOrder = 2 * (elemDisplacementType%order - 1)
               End If
               Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
               Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
               Call MEF90PlasticityEnergySet(myenergy,xSec,damageSec,cumulatedDissipatedPlasticEnergySec,MEF90DefMechCtx%DMVect,MEF90DefMechCtx%DMScal,setIS, &
                                             elemDisplacement,elemDisplacementType,elemScal,elemScalType,matpropSet,ierr)
               Call MEF90Element_Destroy(elemDisplacement,ierr)
               Call MEF90Element_Destroy(elemScal,ierr)
            End If
         Case default
            Print*,__FUNCT__,': Unimplemented damage Type, only AT1Elastic and AT2Elastic implement',cellSetOptions%damageType
            STOP  
         End Select


         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         Call MPI_AllReduce(myenergy,myenergySet,1,MPIU_SCALAR,MPI_SUM,MEF90DefMechCtx%MEF90Ctx%comm,ierr);CHKERRQ(ierr)
         energy(set) = energy(set) + myenergySet
      End Do ! set
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
      End If
      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
         Call SectionRealDestroy(cumulatedDissipatedPlasticEnergySec,ierr);CHKERRQ(ierr)
      End If
      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechPlasticDissipation



#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechElasticEnergy"
!!!
!!!  
!!!  MEF90DefMechElasticEnergy: 
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechElasticEnergy(x,MEF90DefMechCtx,energy,ierr)
      Type(Vec),Intent(IN)                               :: x
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,dimension(:),Pointer                     :: energy
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Class(MEF90_DefMechAT_Type),Allocatable            :: ATModel
      Class(MEF90_DEFMECHSPLIT),Allocatable              :: Split
      PetscBool                                          :: isElastic
      Type(SectionReal)                                  :: xSec
      Type(SectionReal)                                  :: damageSec,plasticStrainSec,temperatureSec
      Type(IS)                                           :: CellSetGlobalIS,cellIS
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemScal
      Type(MEF90Element_Type)                            :: elemDisplacementType,elemScalType
      PetscScalar                                        :: myenergy,myenergySet
      PetscReal,Dimension(:),Pointer                     :: xloc,damageLoc,temperatureLoc,plasticStrainLoc
      PetscReal                                          :: EEDPlus,EEDMinus
      Type(MEF90_MATS)                                   :: plasticStrainCell
      PetscReal                                          :: damageGauss,temperatureGauss,Stress_ZZ_planeStrain
      PetscReal                                          :: elasticEnergyDensityGauss
      Type(MEF90_MATS)                                   :: inelasticStrainGauss,plasticStrainGauss
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscInt                                           :: iDoF1,iGauss

      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,xSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(xSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr) 

      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%CellDMMatSSec,plasticStrainSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticStrainSec,MEF90DefMechCtx%CellDMMatSScatter,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr) 
      Else
         PlasticStrainSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,TemperatureSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(TemperatureSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr) 
      Else
         temperatureSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%damage)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,damageSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(damageSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr) 
      Else
         damageSec%v = 0
      End If

      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)      
      Do set = 1,size(setID)
         myenergy = 0.0_Kr
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
         elemScalType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)
         If (elemDisplacementType%coDim == 0) Then
            Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),cellIS,ierr);CHKERRQ(iErr)

            Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
            !!! compute the contribution of each element
            If (Size(cellID) > 0) Then
               !!! get the ATModel and split objects
               Call MEF90DefMechGetATModel(cellSetOptions,ATModel,isElastic)
               Call MEF90DefMechGetSplit(cellSetOptions,Split)
               quadratureOrder = max(2 * elemScalType%order,ATModel%aOrder + Split%strainOrder + Split%damageOrder)
               Call MEF90Element_Create(MEF90DefMechCtx%DMVect,cellIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
               Call MEF90Element_Create(MEF90DefMechCtx%DMScal,cellIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)

               !!! Allocate storage for fields at dof and Gauss points
               Allocate(xloc(elemDisplacementType%numDof))
               Allocate(temperatureloc(elemScalType%numDof))
               Allocate(damageLoc(elemScalType%numDof))
               Do cell = 1,size(cellID)   
                  !!! Get value of each field at each Dof of the local element
                  Call SectionRealRestrictClosure(xSec,MEF90DefMechCtx%DMVect,cellID(cell),elemDisplacementType%numDof,xloc,ierr);CHKERRQ(ierr)
                  If (temperatureSec%v /= 0) Then
                     Call SectionRealRestrictClosure(temperatureSec,MEF90DefMechCtx%DMScal,cellID(cell),elemScalType%numDof,temperatureLoc,ierr);CHKERRQ(ierr)
                  End If
                  plasticStrainCell = 0.0_Kr
                  If (plasticStrainSec%v /= 0) Then
                     Call SectionRealRestrict(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
                     plasticStrainCell = plasticStrainLoc
                  End If
                  If (damageSec%v /= 0) Then
                     Call SectionRealRestrictClosure(damageSec,MEF90DefMechCtx%DMScal,cellID(cell),elemScalType%numDof,damageLoc,ierr);CHKERRQ(ierr)
                  End If
                  Do iGauss = 1,size(elemDisplacement(cell)%Gauss_C)
                     damageGauss = 0.0_Kr
                     If (damageSec%v /= 0) Then
                        Do iDoF1 = 1, elemScalType%numDof
                           damageGauss = damageGauss + damageLoc(iDoF1) * elemScal(cell)%BF(iDoF1,iGauss)
                        End Do ! iDoF1
                     End If

                     inelasticStrainGauss = 0.0_Kr
                     Do iDoF1 = 1,elemDisplacementType%numDof
                        inelasticStrainGauss = inelasticStrainGauss + xloc(iDof1) * elemDisplacement(cell)%GradS_BF(iDof1,iGauss)
                     End Do

                     temperatureGauss   = 0.0_Kr
                     If (temperatureSec%v /= 0) Then
                        Do iDoF1 = 1,elemScalType%numDof
                           temperatureGauss = temperatureGauss + temperatureLoc(iDof1) * elemScal(cell)%BF(iDof1,iGauss)
                        End Do
                        inelasticStrainGauss = inelasticStrainGauss - (temperatureGauss * matpropSet%LinearThermalExpansion)
                     End If
                     Call Split%EED(inelasticStrainGauss - plasticStrainCell,matpropSet%HookesLaw,EEDPlus,EEDMinus)
                     If (isElastic) Then
                        elasticEnergyDensityGauss = EEDPlus + EEDMinus
                     Else
                        elasticEnergyDensityGauss = (ATmodel%a(damageGauss) + matpropSet%residualStiffness) * EEDPlus + EEDMinus
                     End If
      ! #if MEF90_DIM == 2
      !                if (.NOT. matpropSet%HookesLaw%isPlaneStress) then
      !                   Stress_ZZ_planeStrain = ( matpropSet%HookesLaw%YoungsModulus - 2.0_Kr*matpropSet%HookesLaw%PoissonRatio*matpropSet%HookesLaw%mu )*trace(plasticStrainCell) + matpropSet%HookesLaw%lambda*trace(inelasticStrainGauss)
      !                   elasticEnergyDensityGauss = elasticEnergyDensityGauss + ( Stress_ZZ_planeStrain + matpropSet%HookesLaw%lambda*trace(inelasticStrainGauss - plasticStrainCell) ) * trace(inelasticStrainGauss)
      !                endif
      ! #endif
                  myenergy = myenergy + elemScal(cell)%Gauss_C(iGauss) * elasticEnergyDensityGauss
                  End Do ! Gauss
                  If (plasticStrainSec%v /= 0) Then
                     Call SectionRealRestore(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
                  End If
               End Do ! cell
               ! flops = 7 * size(elemDisplacement(1)%Gauss_C) * size(cellID) 
               ! If (damage%v /= 0) Then
               !    flops = flops + 2 * elemDamageType%numDof * size(elemDamage(1)%Gauss_C) * size(cellID) 
               ! End If
               ! If (temperature%v /= 0) Then
               !    flops = flops + 2 * elemDamageType%numDof * size(elemDamage(1)%Gauss_C) * size(cellID) 
               ! End If
               ! Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
               DeAllocate(xloc)
               DeAllocate(temperatureloc)
               DeAllocate(damageloc)
               Call MEF90Element_Destroy(elemDisplacement,ierr)
               Call MEF90Element_Destroy(elemScal,ierr)
            End If   
            Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
            Call ISDestroy(cellIS,ierr);CHKERRQ(ierr)
            Call MPI_AllReduce(myenergy,myenergySet,1,MPIU_SCALAR,MPI_SUM,MEF90DefMechCtx%MEF90Ctx%comm,ierr);CHKERRQ(ierr)
            energy(set) = energy(set) + myenergySet
         End If
      End Do ! set
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      If (Associated(MEF90DefMechCtx%damage)) Then
         Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)
      End If
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
      End If
      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
      End If
      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechElasticEnergy


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechStress"
!!!
!!!  
!!!  MEF90DefMechStress: 
!!!  
!!!  (c) 2012-20 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechStress(displacement,MEF90DefMechCtx,stress,ierr)
      Type(Vec),Intent(IN)                               :: displacement
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      Type(Vec),Intent(IN)                               :: stress
         
      Type(SectionReal)                                  :: displacementSec,stressSec
      Type(SectionReal)                                  :: plasticStrainSec,temperatureSec,damageSec
      Type(IS)                                           :: CellSetGlobalIS,setIS,setISdof,bcIS
      PetscInt,dimension(:),Pointer                      :: setID,cellID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder,cell
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemDamage
      Type(MEF90Element_Type)                            :: elemDisplacementType,elemDamageType
      PetscReal,Dimension(:),Pointer                     :: displacementDof,damageDof,temperatureDof
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: globalOptions
      Character(len=MEF90_MXSTRLEN)                      :: IOBuffer
      Class(MEF90_DefMechAT_Type),Allocatable            :: ATModel
      Class(MEF90_DEFMECHSPLIT),Allocatable              :: Split
      PetscBool                                          :: cellIsElastic,cellHasForce
      Type(MEF90_MATS)                                   :: inelasticStrainGauss,stressGaussPlus,stressGaussMinus,stressGauss,stressCell
      PetscReal                                          :: cellSize
      PetscReal,Dimension(:),Pointer                     :: stressCellPtr
      PetscReal                                          :: damageGauss,temperatureGauss
      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss


      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)

      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,displacementSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(displacementSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,displacement,ierr);CHKERRQ(ierr) 

      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,temperatureSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(temperatureSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr)   
      Else
         temperatureSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%CellDMMatSSec,plasticStrainSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticStrainSec,MEF90DefMechCtx%CellDMMatSScatter,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)
      Else
         PlasticStrainSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%damage)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,damageSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(damageSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)   
      Else
         damageSec%v = 0
      End If

      Call SectionRealDuplicate(MEF90DefMechCtx%CellDMMatSSec,stressSec,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)      
      Allocate(stressCellPtr(SIZEOFMEF90_MATS))
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DMVect,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
         elemDamageType       = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)

         Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
         If (Size(cellID) > 0) Then
            !!! get the ATModel and split objects
            Call MEF90DefMechGetATModel(cellSetOptions,ATModel,cellIsElastic)
            Call MEF90DefMechGetSplit(cellSetOptions,Split)

            !!! Allocate elements 
            QuadratureOrder = max(2*elemDisplacementType%order,split%damageOrder + split%strainOrder)
            Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
            Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemDamage,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)

            !!! Allocate storage for fields at dof and Gauss points
            Allocate(displacementDof(ElemDisplacementType%numDof))
            Allocate(damageDof(ElemDamageType%numDof))
            If (Associated(MEF90DefMechCtx%temperature)) Then
               Allocate(temperatureDof(ElemDamageType%numDof))
            Else
               Nullify(temperatureDof)
            End If

            Do cell = 1,size(cellID)
               stressCell = 0.0_Kr
               cellSize   = 0.0_Kr
               !!! Get value of each field at each Dof of the local element
               Call SectionRealRestrictClosure(displacementSec,MEF90DefMechCtx%DMVect,cellID(cell),elemDisplacementType%numDof,displacementDof,ierr);CHKERRQ(ierr)
               Call SectionRealRestrictClosure(damageSec,MEF90DefMechCtx%DMScal,cellID(cell),elemDamageType%numDof,damageDof,ierr);CHKERRQ(ierr)

               !!! Loop over Gauss point, evaluate fields and assemble the cell contribution to the residual
               numDofDisplacement = size(elemDisplacement(cell)%BF,1)
               numDofDamage       = size(elemDamage(cell)%BF,1)
               numGauss           = size(elemDisplacement(cell)%BF,2)

               Do iGauss = 1,numGauss
                  If (.NOT. cellIsElastic) Then
                     damageGauss = 0.0_Kr
                     Do iDof1 = 1,numDofDamage
                        damageGauss = damageGauss + damageDof(iDoF1) * elemDamage(cell)%BF(iDoF1,iGauss)
                     End Do ! iDof1
                  End If

                  inelasticStrainGauss = 0.0_Kr
                  Do iDof1 = 1,numDofDisplacement
                     inelasticStrainGauss = inelasticStrainGauss + displacementDof(iDoF1) * elemDisplacement(cell)%GradS_BF(iDoF1,iGauss)
                  End Do

                  temperatureGauss = 0.0_Kr
                  If (Associated(temperatureDof)) Then
                     Do iDoF1 = 1,numDofDamage
                        temperatureGauss = temperatureGauss + temperatureDoF(iDoF1) * elemDamage(cell)%BF(iDoF1,iGauss)
                     End Do ! iDof1
                     inelasticStrainGauss = inelasticStrainGauss - (temperatureGauss * matpropSet%LinearThermalExpansion)
                  End If ! Associated(temperatureDof)

                  Call Split%DEED(inelasticStrainGauss,matpropSet%HookesLaw,stressGaussPlus,stressGaussMinus)
                  If (cellIsElastic) Then
                     stressCell = stressCell + stressGaussPlus + stressGaussMinus * elemDisplacement(cell)%Gauss_C(iGauss)
                  Else
                     stressCell = stressCell + (ATModel%a(damageGauss) * stressGaussPlus + stressGaussMinus) * elemDisplacement(cell)%Gauss_C(iGauss)
                  End If
                  cellSize   = cellSize + elemDisplacement(cell)%Gauss_C(iGauss)
               End Do ! iGauss

               stressCell    = stressCell / cellSize
               stressCellPtr = stressCell
               Call SectionRealUpdate(stressSec,cellID(cell),stressCellPtr,INSERT_VALUES,ierr);CHKERRQ(ierr)
            End Do ! cell
         End If ! cell coDim
      End Do ! set
      deAllocate(stressCellPtr)
      !!! No need to complete since stressSec is cell-based
      Call SectionRealToVec(stressSec,MEF90DefMechCtx%CellDMMatSScatter,SCATTER_FORWARD,stress,ierr);CHKERRQ(ierr)          
      
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call SectionRealDestroy(stressSec,ierr);CHKERRQ(ierr)
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
      End If

      If (Associated(MEF90DefMechCtx%damage)) Then
         Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)
      End If

      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
      End If
      Call SectionRealDestroy(stressSec,ierr);CHKERRQ(ierr)

   End Subroutine MEF90DefMechStress


! #undef __FUNCT__
! #define __FUNCT__ "MEF90DefMechBilinearFormDamageLinSoftLoc"
! !!!
! !!!  
! !!!  MEF90DefMechBilinearFormDamageLinSoftLoc: Assembles the bilinear form for LinSoft, that is:
! !!!    <K(alpha) beta , delta> = \int a''(\alpha) W(e(u)) beta \delta 
! !!!                            + Gc/pi/\ell w''(alpha) beta \delta
! !!!                            + 2 Gc/pi*\ell \nabla beta \nabla \delta
! !!! with a(\alpha) = (1-w(\alpha)) / ( 1 + (k-1) w(\alpha)) and w(\alpha) = 1 - (1-\alpha)^2
! !!! i.e. a''(\alpha) = 2 k * ( k + 3(k-1) v^2) / [k + (1-k) v^2]^3, w''(\alpha) = -2
! !!!  
! !!!  (c) 2014-2015 Erwan Tanne erwan.tanne@gmail.com, Blaise Bourdin bourdin@lsu.edu
! !!!

!    Subroutine MEF90DefMechBilinearFormDamageLinSoftLoc(ALoc,xDof,displacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage)
!       PetscReal,Dimension(:,:),Pointer                   :: Aloc
!       PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,damageDof,temperatureDof
!       Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
!       PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
!       Type(MEF90_MATPROP),Intent(IN)                     :: matprop
!       Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
!       Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage

!       PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
!       PetscReal                                          :: elasticEnergyDensityGauss,temperatureGauss,damageGauss
!       Type(MEF90_MATS)                                   :: inelasticStrainGauss
!       PetscReal                                          :: C2,C1,N,k,C0
!       PetscLogDouble                                     :: flops
!       PetscErrorCode                                     :: ierr

!       numDofDisplacement = size(elemDisplacement%BF,1)
!       numDofDamage = size(elemDamage%BF,1)
!       numGauss = size(elemDamage%BF,2)
      
!       k = matprop%CoefficientLinSoft
!       C1 = 2.0_Kr * matprop%fractureToughness / matprop%internalLength / PETSC_PI
!       C2 = 2.0_Kr * matprop%fractureToughness * matprop%internalLength / PETSC_PI
!       N  = matprop%DuctileCouplingPower

!       Do iGauss = 1,numGauss
!          temperatureGauss = 0.0_Kr
!          If (Associated(temperatureDof)) Then
!             Do iDoF1 = 1,numDofDamage
!                temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
!             End Do
!          End If

!          inelasticStrainGauss = 0.0_Kr
!          Do iDoF1 = 1,numDofDisplacement
!             inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
!          End Do
!          inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion - plasticStrainCell
!          elasticEnergyDensityGauss = (matprop%HookesLaw * inelasticStrainGauss) .DotP. inelasticStrainGauss

!          damageGauss = 0.0_Kr
!          Do iDof1 = 1, numDofDamage
!             damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
!          End Do
!          C0 = k * (k + 3.0_Kr * (k - 1.0_Kr) * (1.0_Kr - damageGauss)**2) / &
!               (k + (1.0_Kr - k) * (1.0_Kr - damageGauss)**2)**3
!          !!! This is really twice the elastic energy density

!          Do iDoF1 = 1,numDofDamage
!             Do iDoF2 = 1,numDofDamage
!                Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) * ( & 
!                                   (elasticEnergyDensityGauss * C0 - C1 + &
!                                   (N * (N - 1.0) *( (1.0_Kr-damageGauss)**(N - 2.0) ) ) * cumulatedDissipatedPlasticEnergyCell ) * elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss) &
!                                   + C2 * ((matprop%toughnessAnisotropyMatrix * elemDamage%Grad_BF(iDoF1,iGauss)) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
!             End Do
!          End Do
!       End Do
!       !flops = 2 * numGauss * numDofDisplacement**2
!       !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
!    End Subroutine MEF90DefMechBilinearFormDamageLinSoftLoc

! #undef __FUNCT__
! #define __FUNCT__ "MEF90DefMechOperatorDamageLinSoftLoc"
! !!!
! !!!  
! !!!  MEF90DefMechOperatorDamageLinSoftLoc: Assembles the operator for LinSoft:
! !!!    F(\alpha)beta = \int [k(\alpha-1)W(e(u))] / [k-(k-1)(1-\alpha)^2]^2 beta 
! !!!                   + 2Gc/pi (1-\alpha) beta / \ell 
! !!!                   + 2Gc/pi \nabla \alpha \nabla beta * \ell
! !!!  
! !!!  (c) 2014-2015 Erwan Tanne erwan.tanne@gmail.com, Blaise Bourdin bourdin@lsu.edu
! !!!

!    Subroutine MEF90DefMechOperatorDamageLinSoftLoc(residualLoc,xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matprop,elemDisplacement,elemDamage,CrackPressureCell)
!       PetscReal,Dimension(:),Pointer                     :: residualLoc
!       PetscReal,Dimension(:),Pointer                     :: xDof,displacementDof,boundaryDisplacementDof,damageDof,temperatureDof
!       Type(MEF90_MATS),Intent(IN)                        :: plasticStrainCell
!       PetscReal,Intent(IN)                               :: cumulatedDissipatedPlasticEnergyCell
!       Type(MEF90_MATPROP),Intent(IN)                     :: matprop
!       Type(MEF90_ELEMENT_ELAST),Intent(IN)               :: elemDisplacement
!       Type(MEF90_ELEMENT_SCAL),Intent(IN)                :: elemDamage
!       PetscReal,Intent(IN)                               :: CrackPressureCell

!       PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
!       PetscReal                                          :: elasticEnergyDensityGauss,temperatureGauss,damageGauss
!       Type(MEF90_MATS)                                   :: inelasticStrainGauss
!       Type(MEF90_VECT)                                   :: gradientDamageGauss
!       PetscReal                                          :: C1,C2,k,C0Operator,N
!       PetscLogDouble                                     :: flops
!       PetscErrorCode                                     :: ierr


!       numDofDisplacement = size(elemDisplacement%BF,1)
!       numDofDamage = size(elemDamage%BF,1)
!       numGauss = size(elemDamage%BF,2)
      
!       C1 = 2.0_Kr * matprop%fractureToughness / matprop%internalLength / PETSC_PI
!       C2 = 2.0_Kr * matprop%fractureToughness * matprop%internalLength / PETSC_PI
!       k = matprop%CoefficientLinSoft
!       N = matprop%DuctileCouplingPower

!       Do iGauss = 1,numGauss
!          temperatureGauss = 0.0_Kr
!          If (Associated(temperatureDof)) Then
!             Do iDoF1 = 1,numDofDamage
!                temperatureGauss = temperatureGauss + elemDamage%BF(iDoF1,iGauss) * temperatureDof(iDoF1)
!             End Do
!          End If

!          inelasticStrainGauss = 0.0_Kr
!          Do iDoF1 = 1,numDofDisplacement
!             inelasticStrainGauss = inelasticStrainGauss + elemDisplacement%GradS_BF(iDoF1,iGauss) * displacementDof(iDoF1)
!          End Do
!          inelasticStrainGauss = inelasticStrainGauss - temperatureGauss * matprop%linearThermalExpansion - plasticStrainCell
!          elasticEnergyDensityGauss = (matprop%HookesLaw * inelasticStrainGauss) .DotP. inelasticStrainGauss
!          !!! This is really twice the elastic energy density

!          damageGauss = 0.0_Kr
!          gradientDamageGauss = 0.0_Kr
!          Do iDoF1 = 1,numDofDamage
!             damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
!             gradientDamageGauss = gradientDamageGauss + elemDamage%Grad_BF(iDoF1,iGauss) * xDof(iDoF1)
!          End Do
!          C0Operator = k * (damageGauss - 1.0_Kr) * elasticEnergyDensityGauss / &
!                      (k + (1.0_Kr - k) * (1.0_kr - damageGauss)**2)**2

!          Do iDoF2 = 1,numDofDamage
!             residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * ( &
!                                  (C0Operator + C1 * (1.0_Kr - damageGauss) &
!                                  - ( N * cumulatedDissipatedPlasticEnergyCell * (1.0_Kr - damageGauss)**(DBLE(N - 1.0))  ) ) * elemDamage%BF(iDoF2,iGauss) & 
!                                   + C2 * ((matprop%toughnessAnisotropyMatrix * gradientDamageGauss) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
!          End Do
!       End Do
!       !flops = 2 * numGauss * numDofDisplacement**2
!       !Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
!    End Subroutine MEF90DefMechOperatorDamageLinSoftLoc



#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechOperatorDamage"
!!!
!!!  
!!!  MEF90DefMechOperatorDamage: Build the operator. When called in SNES, the solution time should always match the target time, 
!!!                                    so there is no need for interpolation of the forcees, external, and boundary values
!!!  
!!!  (c) 2012-20 Blaise Bourdin bourdin@lsu.edu, Erwan Tanne erwan.tanne@gmail.com
!!!

   Subroutine  MEF90DefMechOperatorDamage(snesDamage,damage,residual,MEF90DefMechCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesDamage
      Type(Vec),Intent(IN)                               :: damage
      Type(Vec),Intent(INOUT)                            :: residual
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(SectionReal)                                  :: displacementSec,residualSec
      Type(SectionReal)                                  :: plasticStrainSec,temperatureSec,damageSec,damagePreviousStepSec
      Type(SectionReal)                                  :: boundaryDamageSec,forceSec,pressureForceSec
      Type(SectionReal)                                  :: CrackPressureSec
      PetscReal,Dimension(:),Pointer                     :: damagePtr,residualPtr
      PetscReal,Dimension(:),Pointer                     :: boundaryDamagePtr
      PetscReal,Dimension(:),Pointer                     :: displacementDof,damageDof,damagePreviousStepDof,temperatureDof,residualLoc
      PetscReal,Dimension(:),Pointer                     :: plasticStrainLoc
      PetscReal,Dimension(:),Pointer                     :: CrackPressureLoc
      Type(MEF90_MATS)                                   :: plasticStrainCell
      Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof,bcIS
      PetscInt,Dimension(:),Pointer                      :: setID,cellID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder,cell
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90DefMechVertexSetOptions_Type),Pointer    :: vertexSetOptions
      PetscReal                                          :: CrackPressureCell
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemDamage
      Type(MEF90Element_Type)                            :: elemDisplacementType,elemDamageType
      PetscInt                                           :: p,dof,c,dim,numDof
      PetscInt                                           :: nVal
      PetscReal,Dimension(:),Pointer                     :: nullPtr
      
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: globalOptions
      Class(MEF90_DefMechAT_Type),Allocatable            :: ATModel
      Class(MEF90_DEFMECHSPLIT),Allocatable              :: Split
      !Class(MEF90_DefMechDrivingForce_Type), Allocatable :: drivingForce
      PetscBool                                          :: cellIsElastic
      PetscReal                                          :: damageGauss,damageDampingGauss,temperatureGauss,EEDGaussPlus,EEDGaussMinus,C1,C3,cellDampingCoefficient
      Type(MEF90_VECT)                                   :: gradDamageGauss,displacementGauss,C4
      Type(MEF90_MATS)                                   :: inelasticStrainGauss,C2
      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss


      Call SNESGetDM(snesDamage,MEF90DefMechCtx%DM,ierr);CHKERRQ(ierr)
      !!! Get Section and scatter associated with each field

      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,displacementSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(displacementSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,MEF90DefMechCtx%displacement,ierr);CHKERRQ(ierr) 

      If (globalOptions%dampingCoefficientDamage * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,damagepreviousStepSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(damagepreviousStepSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%damagepreviousStep,ierr);CHKERRQ(ierr) 
      EndIf
      
      Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,damageSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(damageSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,damage,ierr);CHKERRQ(ierr) 

      Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,boundaryDamageSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(boundaryDamageSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%boundaryDamage,ierr);CHKERRQ(ierr) 


      Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,residualSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(residualSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,residual,ierr);CHKERRQ(ierr) 
      Call SectionRealSet(residualSec,0.0_Kr,ierr);CHKERRQ(ierr)
      Call VecSet(residual,0.0_Kr,ierr);CHKERRQ(ierr)

      Nullify(TemperatureDof)
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,temperatureSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(temperatureSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr) 
      Else
         temperatureSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%CrackPressure)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%cellDMScalSec,CrackPressureSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(CrackPressureSec,MEF90DefMechCtx%CellDMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%CrackPressure,ierr);CHKERRQ(ierr)
      Else
         CrackPressureSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%cellDMMatSSec,plasticstrainSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticstrainSec,MEF90DefMechCtx%CellDMMatSScatter,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr) 
      Else
         PlasticStrainSec%v = 0
      End If

      !!! get IS for cell sets
      Call DMMeshGetDimension(MEF90DefMechCtx%DM,dim,ierr);CHKERRQ(ierr)
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      !!!
      !!! We loop over all element twice. The first time in order to assembly all non BC cell sets
      !!! In the second pass, we only update the BC where necessary
      !!! vertex set BC are updated last, so that they override cell set BC
      !!!
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
         elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
         elemDamageType       = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)
         If ((Size(cellID) > 0) .AND. (elemDamageType%coDim == 0)) Then
            !!! get the ATModel and split objects
            Call MEF90DefMechGetATModel(cellSetOptions,ATModel,cellIsElastic)
            Call MEF90DefMechGetSplit(cellSetOptions,Split)
            !Call MEF90DefMechGetDrivingForce(cellSetOption,drivingForce)
            C1 = matpropSet%fractureToughness / ATModel%cw * 0.25_Kr / matpropSet%internalLength
            C2 = matpropSet%fractureToughness / ATModel%cw * 0.5_Kr * matpropSet%internalLength * matpropSet%toughnessAnisotropyMatrix
            cellDampingCoefficient = globalOptions%dampingCoefficientDamage * MEF90DefMechCtx%timeStep

            !!! Allocate storage for fields at dof and Gauss points
            !!! Leaving plasticstrain aside until I can remember how it is supposed to be dealt with
            Allocate(residualLoc(ElemDamageType%numDof))
            Allocate(displacementDof(ElemDisplacementType%numDof))
            Allocate(damageDof(ElemDamageType%numDof))
            If (globalOptions%dampingCoefficientDamage * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
               Allocate(damagePreviousStepDof(ElemDamageType%numDof))
            Else
               Nullify(damagepreviousStepDof)
            End If
            If (Associated(MEF90DefMechCtx%temperature)) Then
               Allocate(temperatureDof(ElemDamageType%numDof))
            End If

            !!! Allocate elements 
            quadratureOrder = max(2 * elemDamageType%order,ATModel%aOrder + Split%strainOrder + Split%damageOrder)
            Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
            Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemDamage,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
            Do cell = 1,size(cellID)
               !!! Get value of each field at each Dof of the local element
               Call SectionRealRestrictClosure(displacementSec,MEF90DefMechCtx%DMVect,cellID(cell),elemDisplacementType%numDof,displacementDof,ierr);CHKERRQ(ierr)
               Call SectionRealRestrictClosure(damageSec,MEF90DefMechCtx%DMScal,cellID(cell),elemDamageType%numDof,damageDof,ierr);CHKERRQ(ierr)
               If (globalOptions%dampingCoefficientDamage * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
                  Call SectionRealRestrictClosure(damagePreviousStepSec,MEF90DefMechCtx%DMScal,cellID(cell),elemDamageType%numDof,damagePreviousStepDof,ierr);CHKERRQ(ierr)
               End If

               If (Associated(MEF90DefMechCtx%temperature)) Then
                  Call SectionRealRestrictClosure(temperatureSec,MEF90DefMechCtx%DMScal,cellID(cell),elemDamageType%numDof,temperatureDof,ierr);CHKERRQ(ierr)
               End If

               If (Associated(MEF90DefMechCtx%CrackPressure)) Then
                  Call SectionRealRestrict(CrackPressureSec,cellID(cell),CrackPressureLoc,ierr);CHKERRQ(ierr)
                  CrackPressureCell = CrackPressureLoc(1)
               Else
                  CrackPressureCell = 0.0_Kr
               End If

               If (Associated(MEF90DefMechCtx%plasticStrain)) Then
                  Call SectionRealRestrict(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
                  plasticStrainCell = plasticStrainLoc
               Else
                  plasticStrainCell = 0.0_Kr
               End If
                  
               !!! Loop over Gauss point, evaluate fields and assemble the cell contribution to the residual
               numDofDisplacement = size(elemDisplacement(cell)%BF,1)
               numDofDamage = size(elemDamage(cell)%BF,1)
               numGauss = size(elemDisplacement(cell)%BF,2)

               residualLoc = 0.0_Kr
               Do iGauss = 1,numGauss
                  !!! Main term: r(\alpha)\beta = a'(\alpha) W(e(u)) \beta + Gc/4/cw [w'(\alpha) / ell \beta + 2 \ell \nabla \alpha \nabla \beta]
                  damageGauss = 0.0_Kr
                  gradDamageGauss = 0.0_Kr
                  Do iDof1 = 1,numDofDamage
                     damageGauss     = damageGauss     + damageDof(iDoF1) * elemDamage(cell)%BF(iDoF1,iGauss)
                     gradDamageGauss = gradDamageGauss + damageDof(iDoF1) * elemDamage(cell)%Grad_BF(iDoF1,iGauss)
                  End Do ! iDof1

                  inelasticStrainGauss = 0.0_Kr
                  displacementGauss    = 0.0_Kr
                  Do iDof1 = 1,numDofDisplacement
                     inelasticStrainGauss = inelasticStrainGauss + displacementDof(iDoF1) * elemDisplacement(cell)%GradS_BF(iDoF1,iGauss)
                     displacementGauss    = displacementGauss    + displacementDof(iDoF1) * elemDisplacement(cell)%BF(iDoF1,iGauss)
                  End Do

                  If (Associated(temperatureDof)) Then
                     temperatureGauss = 0.0_Kr
                     Do iDoF1 = 1,numDofDamage
                        temperatureGauss = temperatureGauss + temperatureDoF(iDoF1) * elemDamage(cell)%BF(iDoF1,iGauss)
                     End Do ! iDof1
                     inelasticStrainGauss = inelasticStrainGauss - (temperatureGauss * matpropSet%LinearThermalExpansion)
                  End If ! Associated(temperatureDof)
#if MEF90_DIM == 2
!!! We need something along these lines
                  !!! Adding terms in planestrain for plasticity with tr(p) = 0
                  ! If (.NOT. matProp%HookesLaw%isPlaneStress) Then
                  !    stressGauss = stressGauss +  stiffness * ( matProp%HookesLaw%lambda*trace(plasticStrainCell)*MEF90MatS2DIdentity )
                  ! End If
#endif

                  If (.NOT. cellIsElastic) Then
                     Call Split%EED(inelasticStrainGauss - plasticStrainCell,matpropSet%HookesLaw,EEDGaussPlus,EEDGaussMinus)
                  Else
                     EEDGaussPlus = 0.0_Kr
                  End If

                  C3 = ATModel%Da(damageGauss) * EEDGaussPlus + C1 * ATModel%Dw(damageGauss)
                  C4 = C2 * gradDamageGauss + CrackPressureCell * displacementGauss
                  Do iDoF2 = 1,numDofDamage
                     residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage(cell)%Gauss_C(iGauss) * ( &
                                          C3 * elemDamage(cell)%BF(iDoF2,iGauss) + ( C4 .dotP. elemDamage(cell)%Grad_BF(iDoF2,iGauss)) )
                  End Do ! iDof2

                  !!! Damping
                  If (cellDampingCoefficient /= 0.0_Kr) Then
                     damageDampingGauss = 0.0_Kr
                     Do iDof1 = 1,numDofDamage
                        damageDampingGauss = damageDampingGauss + ((damageDof(iDoF1) - damagePreviousStepDof(iDoF1)) * elemDamage(cell)%BF(iDoF1,iGauss))
                     End Do
                     damageDampingGauss = damageDampingGauss * cellDampingCoefficient
                     Do iDoF2 = 1,numDofDamage
                        residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage(cell)%Gauss_C(iGauss) * &
                                             damageDampingGauss * elemDamage(cell)%BF(iDoF2,iGauss)
                     End Do
                  End If ! damping

!!! Add driving force here
               ! If (hasLocalOperatorFunction2) Then
               !    Call localOperatorFunction2(residualLoc,damageDof,displacementDof,nullPtr,nullPtr,temperatureDof,plasticStrainCell,cumulatedDissipatedPlasticEnergyCell,matpropSet,elemDisplacement(cell),elemDamage(cell),CrackPressureCell)
               ! End If
!!! Add plastic Dissipation here
! #if MEF90_DIM == 2
!          !!! Adding terms in planestrain for plasticity with tr(p) = 0
!          If (.NOT. matProp%HookesLaw%isPlaneStress) Then
!            Stress_ZZ_planeStrain = ( matprop%HookesLaw%YoungsModulus - 2.0_Kr*matprop%HookesLaw%PoissonRatio*matprop%HookesLaw%mu )*trace(plasticStrainCell) + matprop%HookesLaw%lambda*trace(inelasticStrainGauss)
!            elasticEnergyDensityGauss = elasticEnergyDensityGauss + ( Stress_ZZ_planeStrain + matprop%HookesLaw%lambda*trace(inelasticStrainGauss - plasticStrainCell) ) * trace(plasticStrainCell)
!          End If
! #endif
!          If (N==1.0_Kr) Then
!             Do iDoF2 = 1,numDofDamage
!                residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * ( &
!                                     ( elasticEnergyDensityGauss * (damageGauss - 1.0_Kr) &
!                                     - cumulatedDissipatedPlasticEnergyCell + C1 ) * elemDamage%BF(iDoF2,iGauss) + &
!                                      (( C2 * matprop%toughnessAnisotropyMatrix * gradientDamageGauss +  Displacement*CrackPressureCell ) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
!             End Do
!          Else
!             Do iDoF2 = 1,numDofDamage
!                residualLoc(iDoF2) = residualLoc(iDoF2) + elemDamage%Gauss_C(iGauss) * ( &
!                                     ( elasticEnergyDensityGauss * (damageGauss - 1.0_Kr) &
!                                     - ( N * cumulatedDissipatedPlasticEnergyCell * (1.0_Kr - damageGauss)**(N - 1.0_Kr)  ) + C1 ) * elemDamage%BF(iDoF2,iGauss) + &
!                                      (( C2 * matprop%toughnessAnisotropyMatrix * gradientDamageGauss +  Displacement*CrackPressureCell ) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)) )
!             End Do
!          EndIf

               End Do ! iGauss
               Call SectionRealUpdateClosure(residualSec,MEF90DefMechCtx%DMScal,cellID(cell),residualLoc,ADD_VALUES,ierr);CHKERRQ(iErr)

               If (Associated(MEF90DefMechCtx%plasticStrain)) Then
                  Call SectionRealRestore(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
               End If

               If (Associated(MEF90DefMechCtx%CrackPressure)) Then
                  Call SectionRealRestore(CrackPressureSec,cellID(cell),CrackPressureLoc,ierr);CHKERRQ(ierr)
               EndIf
            End Do ! cell
            DeAllocate(displacementDof)
            DeAllocate(damageDof)
            If (globalOptions%dampingCoefficientDamage * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
               DeAllocate(damagePreviousStepDof)
            End If
            If (Associated(MEF90DefMechCtx%temperature)) Then
               DeAllocate(temperatureDof)
            End If
            DeAllocate(residualLoc)

            Call MEF90Element_Destroy(elemDisplacement,ierr)
            Call MEF90Element_Destroy(elemDamage,ierr)
         End If ! size(cellID) >0
         Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do ! set

      !!! "Ghost update" for the residual Section
      Call SectionRealComplete(residualSec,ierr);CHKERRQ(ierr)
      !!! Scatter back from SectionReal to Vec
      Call SectionRealToVec(residualSec,MEF90DefMechCtx%DMScalScatter,SCATTER_FORWARD,residual,ierr);CHKERRQ(ierr) 

      !!!
      !!! Cell set BC
      !!!
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         If (cellSetOptions%Has_DamageBC) Then
            Call DMMeshGetStratumIS(MEF90DefMechCtx%DMScal,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
            Call MEF90ISCreateCelltoVertex(MEF90DefMechCtx%DMScal,PETSC_COMM_WORLD,setIS,bcIS,ierr)
            Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
            Allocate(damagePtr(nval),stat=ierr)
            Allocate(residualPtr(nval),stat=ierr)
            Allocate(boundaryDamagePtr(nval),stat=ierr)
            Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMScal,damage,damagePtr,bcIS,1,ierr)
            Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMScal,MEF90DefMechCtx%boundaryDamage,boundaryDamagePtr,bcIS,1,ierr)
            residualPtr = damagePtr - boundaryDamagePtr
            Call MEF90VecSetValuesISdof(MEF90DefMechCtx%DMScal,residual,residualPtr,bcIS,1,INSERT_VALUES,ierr)
            DeAllocate(boundaryDamagePtr)
            DeAllocate(residualPtr)
            DeAllocate(damagePtr)
            Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If ! cellSetOptions%Has_BC
      End Do ! set
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
      
      !!!
      !!! Vertex set BC
      !!!
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxVertexSetOptions(MEF90DefMechCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
         If (vertexSetOptions%Has_DamageBC) Then
            Call DMMeshGetStratumIS(MEF90DefMechCtx%DMScal,'Vertex Sets',setID(set),bcIS,ierr);CHKERRQ(iErr)
            Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
            Allocate(damagePtr(nval),stat=ierr)
            Allocate(residualPtr(nval),stat=ierr)
            Allocate(boundaryDamagePtr(nval),stat=ierr)
            Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMScal,damage,damagePtr,bcIS,1,ierr)
            Call MEF90VecGetValuesISdof(MEF90DefMechCtx%DMScal,MEF90DefMechCtx%boundaryDamage,boundaryDamagePtr,bcIS,1,ierr)
            residualPtr = damagePtr - boundaryDamagePtr
            Call MEF90VecSetValuesISdof(MEF90DefMechCtx%DMScal,residual,residualPtr,bcIS,1,INSERT_VALUES,ierr)
            DeAllocate(boundaryDamagePtr)
            DeAllocate(residualPtr)
            DeAllocate(damagePtr)
            Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
         End If ! vertexSetOptions%Has_BC
      End Do ! set
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
   
      Call VecAssemblyBegin(residual,ierr)
      Call VecAssemblyEnd(residual,ierr)

      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
      End If

      If (Associated(MEF90DefMechCtx%CrackPressure)) Then
         Call SectionRealDestroy(CrackPressureSec,ierr);CHKERRQ(ierr)
      End If
      
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
      End If

      Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)
      If (globalOptions%dampingCoefficientDamage * MEF90DefMechCtx%timeStep /= 0.0_Kr) Then
         Call SectionRealDestroy(damagePreviousStepSec,ierr);CHKERRQ(ierr)
      End If
      Call SectionRealDestroy(boundaryDamageSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(residualSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(displacementSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechOperatorDamage

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechBilinearFormDamage"
!!!
!!!  
!!!  MEF90DefMechBilinearFormDamage:
!!!  
!!!  (c) 2012-19 Blaise Bourdin bourdin@lsu.edu, Erwan Tanne erwan.tanne@gmail.com
!!!

   Subroutine MEF90DefMechBilinearFormDamage(snesDamage,damage,A,M,flg,MEF90DefMechCtx,ierr)
      Type(SNES),Intent(IN)                              :: snesDamage
      Type(Vec),Intent(IN)                               :: damage
      Type(Mat),Intent(INOUT)                            :: A,M
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr  
      
      Type(SectionReal)                                  :: displacementSec,damageSec,temperatureSec,plasticStrainSec
      PetscReal,Dimension(:),Pointer                     :: displacementDOF,damageDof,temperatureDof,nullPtr
      PetscReal,Dimension(:),Pointer                     :: plasticStrainLoc
      Type(MEF90_MATS)                                   :: plasticStrainCell
      Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(MEF90DefMechVertexSetOptions_Type),Pointer    :: vertexSetOptions
      Type(MEF90Element_Type)                            :: ElemDisplacementType,ElemDamageType
      PetscInt,Dimension(:),Pointer                      :: setIdx,bcIdx,Cone
      Type(IS)                                           :: bcIS
      PetscInt                                           :: cell,v,numBC,numDoF,numCell,c,dim
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscReal,Dimension(:,:),Pointer                   :: Aloc
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemDamage
      
      Type(MEF90DefMechGlobalOptions_Type),pointer       :: globalOptions
      Class(MEF90_DefMechAT_Type),Allocatable            :: ATModel
      Class(MEF90_DEFMECHSPLIT),Allocatable              :: Split
      !Class(MEF90_DefMechDrivingForce_Type), Allocatable :: drivingForce
      PetscBool                                          :: cellIsElastic
      PetscReal                                          :: damageGauss,damageDampingGauss,temperatureGauss,EEDGaussPlus,EEDGaussMinus,C1,C3,cellDampingCoefficient,BF1
      Type(MEF90_VECT)                                   :: gradDamageGauss,displacementGauss,grad_BF1
      Type(MEF90_MATS)                                   :: inelasticStrainGauss,C2
      PetscInt                                           :: iDoF1,iDoF2,iGauss,numDofDisplacement,numDofDamage,numGauss
      
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,GlobalOptions,ierr);CHKERRQ(ierr)
      Call MatZeroEntries(A,ierr);CHKERRQ(ierr)
      Call SNESGetDM(snesDamage,MEF90DefMechCtx%DM,ierr);CHKERRQ(ierr)
      
      !!! Get Section and scatter associated with each field
      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,displacementSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(displacementSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,MEF90DefMechCtx%displacement,ierr);CHKERRQ(ierr)   

      Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,damageSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(damageSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,damage,ierr);CHKERRQ(ierr)   

      Nullify(temperatureDof)
      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,temperatureSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(temperatureSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr)   
      Else
         temperatureSec%v = 0
      End If
      
      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDuplicate(MEF90DefMechCtx%CellDMMatSSec,plasticStrainSec,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(plasticStrainSec,MEF90DefMechCtx%CellDMMatSScatter,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)
      Else
         PlasticStrainSec%v = 0
      End If

      !!! get IS for cell sets
      Call DMMeshGetDimension(MEF90DefMechCtx%DM,dim,ierr);CHKERRQ(ierr)
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID) 
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         ElemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
         ElemDamageType       = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)
         !!! Go through cells and call local assembly functions
         Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
         If ((Size(cellID) > 0) .AND. (ElemDisplacementType%coDim == 0)) Then
            !!! get the ATModel and split objects
            Call MEF90DefMechGetATModel(cellSetOptions,ATModel,cellIsElastic)
            Call MEF90DefMechGetSplit(cellSetOptions,Split)
            !Call MEF90DefMechGetDrivingForce(cellSetOption,drivingForce)
            C1 = matpropSet%fractureToughness / ATModel%cw * 0.25_Kr / matpropSet%internalLength
            C2 = matpropSet%fractureToughness / ATModel%cw * 0.5_Kr * matpropSet%internalLength * matpropSet%toughnessAnisotropyMatrix
            cellDampingCoefficient = globalOptions%dampingCoefficientDamage * MEF90DefMechCtx%timeStep

            !!! Allocate storage for fields at dof and Gauss points
            Allocate(Aloc(ElemDamageType%numDof,ElemDamageType%numDof))
            Allocate(displacementDof(ElemDisplacementType%numDof))
            Allocate(damageDof(ElemDamageType%numDof))
            If (Associated(MEF90DefMechCtx%temperature)) Then
               Allocate(temperatureDof(ElemDamageType%numDof))
            End If
            
            !!! Allocate elements 
            quadratureOrder = max(2 * elemDamageType%order,ATModel%aOrder + Split%strainOrder + Split%damageOrder)
            Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
            Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemDamage,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
            Do cell = 1,size(cellID)
               !!! Get value of each field at each Dof of the local element
               Call SectionRealRestrictClosure(displacementSec,MEF90DefMechCtx%DMVect,cellID(cell),elemDisplacementType%numDof,displacementDof,ierr);CHKERRQ(ierr)
               Call SectionRealRestrictClosure(damageSec,MEF90DefMechCtx%DMScal,cellID(cell),elemDamageType%numDof,damageDof,ierr);CHKERRQ(ierr)
               If (Associated(MEF90DefMechCtx%temperature)) Then
                  Call SectionRealRestrictClosure(temperatureSec,MEF90DefMechCtx%DMScal,cellID(cell),elemDamageType%numDof,temperatureDof,ierr);CHKERRQ(ierr)
               End If
               If (Associated(MEF90DefMechCtx%plasticStrain)) Then
                  Call SectionRealRestrict(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
                  plasticStrainCell = plasticStrainLoc
               Else
                 plasticStrainCell = 0.0_Kr
               End If
                  
               !!! Loop over Gauss point, evaluate fields and assemble the cell contribution to the bilinear form
               numDofDisplacement = size(elemDisplacement(cell)%BF,1)
               numDofDamage = size(elemDamage(cell)%BF,1)
               numGauss = size(elemDisplacement(cell)%BF,2)
               Aloc = 0.0_Kr
               Do iGauss = 1,numGauss
                  !!! Main term: A(\alpha)\beta,\gamma = a''(\alpha) W(e(u)) \beta \gamma + Gc/4/cw [w''(\alpha)/ ell \beta \gamma + 2 \ell \nabla \beta \cdot \nabla \gamma]
                  damageGauss = 0.0_Kr
                  gradDamageGauss = 0.0_Kr
                  Do iDof1 = 1,numDofDamage
                     damageGauss     = damageGauss     + damageDof(iDoF1) * elemDamage(cell)%BF(iDoF1,iGauss)
                     gradDamageGauss = gradDamageGauss + damageDof(iDoF1) * elemDamage(cell)%Grad_BF(iDoF1,iGauss)
                  End Do ! iDof1

                  inelasticStrainGauss = 0.0_Kr
                  displacementGauss    = 0.0_Kr
                  Do iDof1 = 1,numDofDisplacement
                     inelasticStrainGauss = inelasticStrainGauss + displacementDof(iDoF1) * elemDisplacement(cell)%GradS_BF(iDoF1,iGauss)
                     displacementGauss    = displacementGauss    + displacementDof(iDoF1) * elemDisplacement(cell)%BF(iDoF1,iGauss)
                  End Do

                  If (Associated(temperatureDof)) Then
                     temperatureGauss = 0.0_Kr
                     Do iDoF1 = 1,numDofDamage
                        temperatureGauss = temperatureGauss + temperatureDoF(iDoF1) * elemDamage(cell)%BF(iDoF1,iGauss)
                     End Do ! iDof1
                     inelasticStrainGauss = inelasticStrainGauss - (temperatureGauss * matpropSet%LinearThermalExpansion)
                  End If ! Associated(temperatureDof)
#if MEF90_DIM == 2
!!! We need something along these lines
                  !!! Adding terms in planestrain for plasticity with tr(p) = 0
                  ! If (.NOT. matProp%HookesLaw%isPlaneStress) Then
                  !    stressGauss = stressGauss +  stiffness * ( matProp%HookesLaw%lambda*trace(plasticStrainCell)*MEF90MatS2DIdentity )
                  ! End If
#endif
                  If (.NOT. cellIsElastic) Then
                     Call Split%EED(inelasticStrainGauss - plasticStrainCell,matpropSet%HookesLaw,EEDGaussPlus,EEDGaussMinus)
                  Else
                     EEDGaussPlus = 0.0_Kr
                  End If

                  C3 = ATModel%D2a(damageGauss) * EEDGaussPlus + C1 * ATModel%D2w(damageGauss)
                  Do iDoF1 = 1,numDofDamage
                     BF1      = C3 * elemDamage(cell)%BF(iDoF1,iGauss)
                     grad_BF1 = C2 * elemDamage(cell)%Grad_BF(iDoF1,iGauss)
                     Do iDoF2 = 1,numDofDamage
                        Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage(cell)%Gauss_C(iGauss) * ( &
                                           BF1 * elemDamage(cell)%BF(iDoF2,iGauss) + &
                                           (grad_BF1 .dotP. elemDamage(cell)%Grad_BF(iDoF2,iGauss)) )
                     End Do
                  End Do
                  If (cellDampingCoefficient /= 0.0_Kr) Then
                     Do iDoF1 = 1,numDofDamage
                        BF1 = cellDampingCoefficient * elemDamage(cell)%BF(iDoF1,iGauss)
                        Do iDoF2 = 1,numDofDamage
                              ALoc(iDoF2,iDoF1) = ALoc(iDoF2,iDoF1) + elemDamage(cell)%Gauss_C(iGauss) * BF1 * elemDamage(cell)%BF(iDoF2,iGauss)
                        End Do ! iDoF1
                     End Do ! iDoF2
                  End If ! cellDampingCoefficient
!!! Need to add plastic Dissipation
! #if MEF90_DIM == 2
!          !!! Adding terms in planestrain for plasticity with tr(p) = 0
!          If (.NOT. matProp%HookesLaw%isPlaneStress) Then
!             Stress_ZZ_planeStrain = ( matprop%HookesLaw%YoungsModulus - 2.0_Kr*matprop%HookesLaw%PoissonRatio*matprop%HookesLaw%mu )*trace(plasticStrainCell) + matprop%HookesLaw%lambda*trace(inelasticStrainGauss)
!             elasticEnergyDensityGauss = elasticEnergyDensityGauss + ( Stress_ZZ_planeStrain + matprop%HookesLaw%lambda*trace(inelasticStrainGauss - plasticStrainCell) ) * trace(plasticStrainCell)
!          End If
! #endif
!          !!! This is really twice the elastic energy density

!          If ((N == 1.0_Kr) .OR. (cumulatedDissipatedPlasticEnergyCell == 0.0_Kr)) Then
!             Do iDoF1 = 1,numDofDamage
!                Do iDoF2 = 1,numDofDamage
!                   Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) * ( &
!                                          elasticEnergyDensityGauss * elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss) + &
!                                          C2 * ( (matprop%toughnessAnisotropyMatrix * elemDamage%Grad_BF(iDoF1,iGauss)) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
!                End Do
!             End Do
!          Else
!             damageGauss = 0.0_Kr
!             Do iDoF1 = 1,numDofDamage
!                damageGauss = damageGauss + elemDamage%BF(iDoF1,iGauss) * xDof(iDoF1)
!             End Do
!             Do iDoF1 = 1,numDofDamage
!                Do iDoF2 = 1,numDofDamage
!                   Aloc(iDoF2,iDoF1) = Aloc(iDoF2,iDoF1) + elemDamage%Gauss_C(iGauss) * ( &
!                                          (elasticEnergyDensityGauss + &
!                                          ( N * (N - 1.0_Kr) *( (1.0_Kr-damageGauss)**(N - 2.0_Kr) ) ) * cumulatedDissipatedPlasticEnergyCell ) * elemDamage%BF(iDoF1,iGauss) * elemDamage%BF(iDoF2,iGauss) + &
!                                          C2 * ((matprop%toughnessAnisotropyMatrix * elemDamage%Grad_BF(iDoF1,iGauss)) .dotP. elemDamage%Grad_BF(iDoF2,iGauss)))
!                End Do
!             End Do
!          EndIf

               End Do ! iGauss
               Call DMmeshAssembleMatrix(A,MEF90DefMechCtx%DM,damageSec,cellID(cell),ALoc,ADD_VALUES,ierr);CHKERRQ(ierr)

               If (Associated(MEF90DefMechCtx%plasticStrain)) Then
                  Call SectionRealRestore(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
               End If
            End Do ! cell
            DeAllocate(displacementDof)
            DeAllocate(damageDof)
            If (Associated(MEF90DefMechCtx%temperature)) Then
               DeAllocate(temperatureDof)
            End If
            Call MEF90Element_Destroy(elemDisplacement,ierr)
            Call MEF90Element_Destroy(elemDamage,ierr)
            DeAllocate(Aloc)    
         End If 
         Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do ! set
      Call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
      Call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)

      !!!
      !!! Boundary conditions at cell sets
      !!!
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         If (cellSetOptions%Has_damageBC) Then
            Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
            Call MEF90ISCreateCelltoVertex(MEF90DefMechCtx%DM,PETSC_COMM_WORLD,setIS,bcIS,ierr)
            Call DMMeshISCreateISglobaldof(MEF90DefMechCtx%DM,bcIS,0,setISdof,ierr);CHKERRQ(ierr)
            Call MatZeroRowsColumnsIS(A,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
            Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
            Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If ! cellSetOptions%Has_damageBC
      End Do     
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
   
      !!!
      !!! Boundary conditions at vertex sets
      !!!
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DM,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxVertexSetOptions(MEF90DefMechCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
         If (vertexSetOptions%Has_damageBC) Then
            Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
            Call DMMeshISCreateISglobaldof(MEF90DefMechCtx%DM,setIS,0,setISdof,ierr);CHKERRQ(ierr)
            Call MatZeroRowsColumnsIS(A,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
            Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
            Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         End If
      End Do
      Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      
      Call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
      Call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)

      Call SectionRealDestroy(displacementSec,ierr);CHKERRQ(ierr)   
      Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)

      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
      End If
      
      If (Associated(MEF90DefMechCtx%plasticStrain)) Then
         Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
      End If
      flg = SAME_NONZERO_PATTERN
   End Subroutine MEF90DefMechBilinearFormDamage


! #undef __FUNCT__
! #define __FUNCT__ "MEF90DefMechSurfaceEnergySetLinSoft"
! !!!
! !!!  
! !!!  MEF90DefMechSurfaceEnergySetLinSoft:
! !!!  
! !!!  (c) 2014-2020 Blaise Bourdin bourdin@lsu.edu
! !!!
!    Subroutine MEF90DefMechSurfaceEnergySetLinSoft(energy,Damage,meshScal,cellIS,matprop,elemScal,elemScalType,ierr)
!       PetscReal,Intent(INOUT)                            :: energy
!       Type(SectionReal),Intent(IN)                       :: Damage
!       Type(DM),Intent(IN)                                :: meshScal
!       Type(IS),Intent(IN)                                :: cellIS
!       Type(MEF90_MATPROP),Intent(IN)                     :: matprop
!       Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elemScal
!       Type(MEF90Element_Type),Intent(IN)                 :: elemScalType
!       PetscErrorCode,Intent(OUT)                         :: ierr

!       PetscReal,Dimension(:),Pointer                     :: DamageLoc
!       PetscReal                                          :: DamageElem
!       Type(MEF90_VECT)                                   :: gradDamageElem
!       PetscInt,Dimension(:),Pointer                      :: cellID
!       PetscInt                                           :: cell
!       PetscInt                                           :: iDoF1,iGauss
!       PetscReal                                          :: C1, C2
!       PetscLogDouble                                     :: flops

      
!       C1 = matprop%fractureToughness / matprop%internalLength / 4.0_Kr / cwLinSoft
!       C2 = matprop%fractureToughness * matprop%internalLength / 4.0_Kr / cwLinSoft
!       Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
!       If (Size(cellID) > 0) Then
!          Allocate(Damageloc(elemScalType%numDof))
!          Do cell = 1,size(cellID)   
!             Call SectionRealRestrictClosure(Damage,meshScal,cellID(cell),elemScalType%numDof,DamageLoc,ierr);CHKERRQ(ierr)
!             Do iGauss = 1,size(elemScal(cell)%Gauss_C)
!                DamageElem     = 0.0_Kr
!                gradDamageElem = 0.0_Kr
!                Do iDoF1 = 1,elemScalType%numDof
!                   DamageElem     = DamageElem     + DamageLoc(iDof1) * elemScal(cell)%BF(iDoF1,iGauss)
!                   gradDamageElem = gradDamageElem + DamageLoc(iDof1) * elemScal(cell)%Grad_BF(iDoF1,iGauss)
!                End Do
!                energy = energy + elemScal(cell)%Gauss_C(iGauss) * ((1.0_Kr - (1.0_Kr- DamageElem)**2 )* C1 + ((matprop%toughnessAnisotropyMatrix * gradDamageElem) .dotP. gradDamageElem) * C2)
!             End Do ! Gauss
!          End Do ! cell
!          flops = (2 * elemScalType%numDof + 5 )* size(elemScal(1)%Gauss_C) * size(cellID) 
!          Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
!          DeAllocate(Damageloc)
!       End If   
!       Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
!    End Subroutine MEF90DefMechSurfaceEnergySetLinSoft


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSurfaceEnergy"
!!!
!!!  
!!!  MEF90DefMechSurfaceEnergy:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90DefMechSurfaceEnergy(Damage,MEF90DefMechCtx,energy,ierr)
      Type(Vec),Intent(IN)                               :: Damage
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,dimension(:),Pointer                     :: energy
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Class(MEF90_DefMechAT_Type),Allocatable            :: ATModel
      PetscBool                                          :: isElastic
      Type(IS)                                           :: CellSetGlobalIS,cellIS
      PetscInt,dimension(:),Pointer                      :: setID
      PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      PetscInt                                           :: set,QuadratureOrder
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(SectionReal)                                  :: DamageSec
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemScal
      Type(MEF90Element_Type)                            :: elemScalType
      PetscReal                                          :: myenergy,myenergySet
      PetscReal,Dimension(:),Pointer                     :: DamageLoc
      PetscReal                                          :: DamageGauss
      Type(MEF90_VECT)                                   :: gradDamageGauss
      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscInt                                           :: cell
      PetscInt                                           :: iDoF1,iGauss
      PetscReal                                          :: C1
      Type(MEF90_MATS)                                   :: C2

      !!! Scatter data from Vec to Sec, or initialize
      Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,DamageSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(DamageSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,Damage,ierr);CHKERRQ(ierr) 

      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      energy = 0.0_Kr
      Do set = 1,size(setID)
         myenergy = 0.0_Kr
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         elemScalType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)

         !!! Call proper local assembly depending on the type of damage law
         If (elemScalType%coDim == 0) Then
            Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),cellIS,ierr);CHKERRQ(iErr)
            Call MEF90DefMechGetATModel(cellSetOptions,ATModel,isElastic)
            QuadratureOrder = max(ATmodel%wOrder,2*(elemScalType%order-1))
            Call MEF90Element_Create(MEF90DefMechCtx%DMScal,cellIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
            !!! Do the cell based computations
            C1 =  matpropSet%fractureToughness / matpropSet%internalLength * 0.25_Kr / ATModel%cw
            C2 = (matpropSet%fractureToughness * matpropSet%internalLength * 0.25_Kr / ATModel%cw) * matPropSet%toughnessAnisotropyMatrix
            Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)

            If (Size(cellID) > 0) Then
               Allocate(Damageloc(elemScalType%numDof))
               Do cell = 1,size(cellID)   
                  Call SectionRealRestrictClosure(DamageSec,MEF90DefMechCtx%DMScal,cellID(cell),elemScalType%numDof,DamageLoc,ierr);CHKERRQ(ierr)
                  Do iGauss = 1,size(elemScal(cell)%Gauss_C)
                     DamageGauss     = 0.0_Kr
                     gradDamageGauss = 0.0_Kr
                     Do iDoF1 = 1,elemScalType%numDof
                        DamageGauss     = DamageGauss     + DamageLoc(iDof1) * elemScal(cell)%BF(iDoF1,iGauss)
                        gradDamageGauss = gradDamageGauss + DamageLoc(iDof1) * elemScal(cell)%Grad_BF(iDoF1,iGauss)
                     End Do
                     myenergy = myenergy + elemScal(cell)%Gauss_C(iGauss) * ( C1 * ATModel%w(DamageGauss) + ((C2 * gradDamageGauss) .dotP. gradDamageGauss) )
                  End Do ! Gauss
               End Do ! cell
               DeAllocate(Damageloc)
            End If   

            Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
            Call MEF90Element_Destroy(elemScal,ierr)
            Call ISDestroy(cellIS,ierr);CHKERRQ(ierr)
            Call MPI_AllReduce(myenergy,myenergySet,1,MPIU_SCALAR,MPI_SUM,MEF90DefMechCtx%MEF90Ctx%comm,ierr);CHKERRQ(ierr)
            energy(set) = myenergySet
         End If ! coDim
      End Do ! set
      Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr) 

      Call SectionRealDestroy(DamageSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechSurfaceEnergy


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCrackVolume"
!!!
!!!  
!!!  MEF90DefMechCrackVolume:
!!!  
!!!  (c) 2016 Erwan erwan.tanne@gmail.com  
!!!

   Subroutine MEF90DefMechCrackVolume(xVec,MEF90DefMechCtx,CrackVolume,ierr)
      Type(Vec),Intent(IN)                               :: xVec
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscReal,Dimension(:),Pointer                     :: CrackVolume
      PetscErrorCode,Intent(OUT)                         :: ierr
   
      Type(SectionReal)                                  :: xSec,alphaSec,CrackPressureSec
      Type(IS)                                           :: CellSetGlobalIS,setIS
      PetscInt,Dimension(:),Pointer                      :: setID
      PetscInt                                           :: set,QuadratureOrder,numCell
      Type(MEF90_MATPROP),Pointer                        :: matpropSet
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elemScal
      Type(MEF90Element_Type)                            :: elemDisplacementType,elemScalType
      Type(MEF90DefMechGlobalOptions_Type),Pointer       :: globalOptions
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      PetscReal                                          :: myCrackVolume,myCrackVolumeSet

      PetscInt,Dimension(:),Pointer                      :: cellID
      PetscReal,Dimension(:),Pointer                     :: xloc,alphaloc,CrackPressureLoc
      PetscInt                                           :: cell
      PetscReal                                          :: CrackPressureElem,alphaElem
      PetscInt                                           :: iDoF1,iGauss
      Type(MEF90_VECT)                                   :: gradAlphaElem,DisplacementElem
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elemDisp

      
      
      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%globalOptionsBag,globalOptions,ierr);CHKERRQ(ierr)

      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,xSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(xSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,xVec,ierr);CHKERRQ(ierr)
      
      Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,alphaSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(alphaSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%Damage,ierr);CHKERRQ(ierr) 
      
      CrackVolume   = 0.0_Kr
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         myCrackVolume = 0.0_Kr
         
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
         elemScalType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)

         QuadratureOrder = 2 * (elemScalType%order - 1)
         Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
         Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)

         If (elemScalType%coDim == 0) Then
            
            Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
            If (Size(cellID) > 0) Then
               Allocate(xloc(elemDisplacementType%numDof))
               Allocate(alphaloc(elemScalType%numDof))
               Do cell = 1,size(cellID)

                  Call SectionRealRestrictClosure(xSec,MEF90DefMechCtx%DMVect,cellID(cell),elemDisplacementType%numDof,xloc,ierr);CHKERRQ(ierr)
                  Call SectionRealRestrictClosure(alphaSec,MEF90DefMechCtx%DMScal,cellID(cell),elemScalType%numDof,alphaLoc,ierr);CHKERRQ(ierr)

                  Do iGauss = 1,size(elemScal(cell)%Gauss_C)

                     DisplacementElem = 0.0_Kr
                     Do iDoF1 = 1,elemDisplacementType%numDof
                        DisplacementElem = DisplacementElem + xloc(iDof1) * elemDisplacement(cell)%BF(iDof1,iGauss)
                     End Do
                     gradAlphaElem = 0.0_Kr
                     alphaElem = 0.0_Kr
                     Do iDoF1 = 1,elemScalType%numDof
                        alphaElem     = alphaElem     + alphaLoc(iDof1) * elemScal(cell)%BF(iDoF1,iGauss)
                        gradAlphaElem = gradAlphaElem + alphaLoc(iDof1) * elemScal(cell)%Grad_BF(iDoF1,iGauss)
                     End Do
                     myCrackVolume = myCrackVolume - (elemDisplacement(cell)%Gauss_C(iGauss) *( ( gradAlphaElem) .dotP. DisplacementElem  ) )
                  End Do ! Gauss
               End Do ! cell

               DeAllocate(xloc)
               DeAllocate(alphaloc)
               End If
            Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
         End If

         Call MEF90Element_Destroy(elemScal,ierr)
         Call MEF90Element_Destroy(elemDisplacement,ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
         Call MPI_AllReduce(myCrackVolume,myCrackVolumeSet,1,MPIU_SCALAR,MPI_SUM,MEF90DefMechCtx%MEF90Ctx%comm,ierr);CHKERRQ(ierr)
         CrackVolume(set) = CrackVolume(set) + myCrackVolumeSet
      End Do ! set
      Call ISrestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr) 

      Call SectionRealDestroy(alphaSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechCrackVolume



#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCrackPressureRescaling"
!!!
!!!  
!!!  MEF90DefMechCrackPressureRescaling:
!!!  
!!!  (c) 2016 Erwan erwan.tanne@gmail.com  
!!!

   Subroutine MEF90DefMechCrackPressureRescaling(xVec,CrackPressureVec,MEF90DefMechCtx,CrackVolumeSet,ActivatedCrackPressureBlocksList,timeStep,ierr)
      Type(Vec),Intent(IN)                               :: xVec
      Type(Vec),Intent(INOUT)                            :: CrackPressureVec
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      PetscReal,Dimension(:),Pointer,Intent(IN)          :: CrackVolumeSet
      PetscReal,Intent(IN)                               :: timeStep
      PetscBool,Dimension(:),Pointer,Intent(IN)          :: ActivatedCrackPressureBlocksList


      Type(MEF90DefMechGlobalOptions_Type),Pointer       :: globalOptions
      Type(SectionReal)                                  :: xSec,CrackPressureSec,alphaSec
      Type(IS)                                           :: CellSetGlobalIS,setIS
      PetscInt,Dimension(:),Pointer                      :: setID,cellID
      PetscInt                                           :: set,QuadratureOrder,numCell
      Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
      Type(MEF90_ELEMENT_SCAL), Dimension(:), Pointer    :: elemScal
      Type(MEF90Element_Type)                            :: elemDisplacementType,elemScalType
      PetscReal                                          :: cellSize,CrackPressureElem,alphaElem,PressureVolumeElem,CrackVolumeControlled
      Type(MEF90_MATS)                                   :: StrainElem
      Type(MEF90_VECT)                                   :: gradAlphaElem,DisplacementElem
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      PetscInt                                           :: iDoF1,iGauss,cell
      PetscReal,Dimension(:),Pointer                     :: xloc,alphaloc,CrackPressureLoc
      Type(MEF90_ELEMENT_ELAST), Dimension(:), Pointer   :: elemDisp
      Type(MEF90_MATPROP),pointer                        :: matpropSet


      Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%globalOptionsBag,globalOptions,ierr);CHKERRQ(ierr)

      Call SectionRealDuplicate(MEF90DefMechCtx%DMVectSec,xSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(xSec,MEF90DefMechCtx%DMVectScatter,SCATTER_REVERSE,xVec,ierr);CHKERRQ(ierr)

      Call SectionRealDuplicate(MEF90DefMechCtx%DMScalSec,alphaSec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(alphaSec,MEF90DefMechCtx%DMScalScatter,SCATTER_REVERSE,MEF90DefMechCtx%Damage,ierr);CHKERRQ(ierr) 

      Call SectionRealDuplicate(MEF90DefMechCtx%cellDMScalSec,CrackPressureSec,ierr)
      Call SectionRealToVec(CrackPressureSec,MEF90DefMechCtx%CellDMScalScatter,SCATTER_REVERSE,CrackPressureVec,ierr);CHKERRQ(ierr) 

      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)

         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
         elemScalType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)
         !Allocate(ActivatedCrackPressureBlocksList(size(MEF90DefMechCtx%CellSetOptionsBag)))

         QuadratureOrder = 2 * (elemScalType%order - 1)
         Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
         Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)

         If (elemScalType%coDim == 0) Then

            Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
            If (Size(cellID) > 0) Then
               Allocate(xloc(elemDisplacementType%numDof))
               Allocate(alphaloc(elemScalType%numDof))

               PressureVolumeElem = 0.0_Kr
               alphaElem = 0.0_Kr
               Do cell = 1,size(cellID)
                  cellSize = 0.0_Kr
                  PressureVolumeElem = 0.0_Kr

                  Call SectionRealRestrictClosure(xSec,MEF90DefMechCtx%DMVect,cellID(cell),elemDisplacementType%numDof,xloc,ierr);CHKERRQ(ierr)
                  Call SectionRealRestrictClosure(alphaSec,MEF90DefMechCtx%DMScal,cellID(cell),elemScalType%numDof,alphaLoc,ierr);CHKERRQ(ierr)
                  Call SectionRealRestrict(CrackPressureSec,cellID(cell),CrackPressureLoc,ierr);CHKERRQ(ierr)
                  Do iGauss = 1,size(elemScal(cell)%Gauss_C)
                     StrainElem = 0.0_Kr
                     Do iDoF1 = 1,elemDisplacementType%numDof
                        StrainElem = StrainElem + xloc(iDof1) *  elemDisplacement(cell)%GradS_BF(iDof1,iGauss)
                     End Do

                     DisplacementElem = 0.0_Kr
                     Do iDoF1 = 1,elemDisplacementType%numDof
                        DisplacementElem = DisplacementElem + xloc(iDof1) * elemDisplacement(cell)%BF(iDof1,iGauss)
                     End Do

                     gradAlphaElem = 0.0_Kr
                     alphaElem = 0.0_Kr
                     Do iDoF1 = 1,elemScalType%numDof
                        alphaElem     = alphaElem     + alphaLoc(iDof1) * elemScal(cell)%BF(iDoF1,iGauss)
                        gradAlphaElem = gradAlphaElem + alphaLoc(iDof1) * elemScal(cell)%Grad_BF(iDoF1,iGauss)
                     End Do
                     cellSize = cellSize + elemDisplacement(cell)%Gauss_C(iGauss)
                     PressureVolumeElem = PressureVolumeElem - (elemDisplacement(cell)%Gauss_C(iGauss) *( gradAlphaElem .dotP. DisplacementElem  ) )
                  End Do ! Gauss

                  !! Case where the total volume is controlled
                  If (ActivatedCrackPressureBlocksList(set)) then
                     !! rescal the total 
                     CrackVolumeControlled = sum(CrackVolumeSet,MASK=ActivatedCrackPressureBlocksList)
                     CrackPressureLoc = ( timeStep/CrackVolumeControlled )*CrackPressureLoc
                  Else 
                     If (trace(StrainElem)<0 .and. (alphaElem > 0.0_Kr)) Then
                        CrackPressureLoc = (1.0_Kr-alphaElem)**2*((matpropSet%HookesLaw*StrainElem)*gradAlphaElem)*gradAlphaElem / (gradAlphaElem .dotP. gradAlphaElem)
                     Else
                        CrackPressureLoc = 0.0_Kr
                     EndIf
                  EndIf
                  Call SectionRealRestore(CrackPressureSec,cellID(cell),CrackPressureLoc,ierr);CHKERRQ(ierr)
               End Do ! cell

               DeAllocate(xloc)
               DeAllocate(alphaloc)
               End If
            Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
         End If

         Call MEF90Element_Destroy(elemScal,ierr)
         Call MEF90Element_Destroy(elemDisplacement,ierr)
         Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      End Do ! set
      Call ISrestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call SectionRealToVec(CrackPressureSec,MEF90DefMechCtx%CellDMScalScatter,SCATTER_FORWARD,CrackPressureVec,ierr);CHKERRQ(ierr) 
      
      Call SectionRealDestroy(alphaSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(CrackpressureSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90DefMechCrackPressureRescaling

End Module MEF90_APPEND(m_MEF90_DefMechAssembly,MEF90_DIM)D
