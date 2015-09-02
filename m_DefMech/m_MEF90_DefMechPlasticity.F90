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
      Type(MEF90_MATS)            :: InelasticStrain
      Type(MEF90_MATS)            :: PlasticStrainOld
      real(Kind = Kr)             :: Damage
      real(Kind = Kr)             :: ValueOfk
   end type MEF90DefMechPlasticityCtx

contains
!!! Since these functions are C interoperable, fortran cannot rename them in a use statement
!!! We play with the pre-processor in order to avoid duplicate symbols.

#define FHG_VONMISES MEF90_APPEND(fhg_VonMises,MEF90_DIM)D

#define FHG_TRESCA MEF90_APPEND(fhg_Tresca,MEF90_DIM)D

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
      type(c_ptr),intent(in),value              :: myctx

      type(MEF90DefMechPlasticityCtx),pointer   :: myctx_ptr
      type(MEF90_MATS)                          :: xMatS

      xMatS = x(1:SIZEOFMEF90_MATS)
      !!! This is the fortran equivalent of casting ctx into a c_ptr
      call c_f_pointer(myctx,myctx_ptr)

      if (myctx_ptr%ValueOfk==0) then
         f(1) = ( (myctx_ptr%HookesLaw * (xMatS-myctx_ptr%PlasticStrainOld)) .DotP. (xMatS-myctx_ptr%PlasticStrainOld) ) /2.0 * (1-myctx_ptr%Damage)**2
         g(1) =  sqrt( MEF90_DIM * ( deviatoricPart(myctx_ptr%HookesLaw*(myctx_ptr%InelasticStrain-xMatS))  .DotP.  deviatoricPart(myctx_ptr%HookesLaw*(myctx_ptr%InelasticStrain-xMatS)) ))  - myctx_ptr%YieldStress
      else 
         f(1) = ( (myctx_ptr%HookesLaw * (xMatS-myctx_ptr%PlasticStrainOld)) .DotP. (xMatS-myctx_ptr%PlasticStrainOld) ) /2.0 * ((1.0_Kr - myctx_ptr%Damage)**2 /( 1.0_Kr + ( myctx_ptr%ValueOfk - 1.0_Kr )*(1.0_Kr - (1.0_Kr - myctx_ptr%Damage)**2 ) ))
         g(1) =  sqrt( MEF90_DIM * ( deviatoricPart(myctx_ptr%HookesLaw*(myctx_ptr%InelasticStrain-xMatS))  .DotP.  deviatoricPart(myctx_ptr%HookesLaw*(myctx_ptr%InelasticStrain-xMatS)) ))*((1.0_Kr - myctx_ptr%Damage)**2 /( 1.0_Kr + ( myctx_ptr%ValueOfk - 1.0_Kr )*(1.0_Kr - (1.0_Kr - myctx_ptr%Damage)**2 ) ))  - myctx_ptr%YieldStress* (1-myctx_ptr%Damage)**2
      end if
      h(1) = Trace(xMatS)
   end subroutine FHG_VONMISES





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
      call EigenVectorValues(deviatoricPart(myctx_ptr%HookesLaw*myctx_ptr%InelasticStrain),MatProj,MatDiag)
      ! D=P^(-1).A.P 
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

   Subroutine MEF90DefMechPlasticStrainUpdate(MEF90DefMechCtx,plasticStrain,x,PlasticStrainOld,ierr)
      use,intrinsic :: iso_c_binding
#ifdef MEF90_HAVE_SNLP
      use SNLPF90
#endif

      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      Type(Vec),Intent(INOUT)                            :: plasticStrain
      Type(Vec),Intent(IN)                               :: x,PlasticStrainOld
      PetscErrorCode,Intent(OUT)                         :: ierr

#ifdef MEF90_HAVE_SNLP
      Type(DM)                                           :: Mesh
      Type(SectionReal)                                  :: plasticStrainSec,plasticStrainOldSec,inelasticStrainSec
      PetscReal,Dimension(:),Pointer                     :: plasticStrainLoc,plasticStrainOldLoc,inelasticStrainLoc,damageLoc
      type(c_funptr)                                     :: snlp_fhg,snlp_Dfhg
      integer(kind=c_int)                                :: snlp_n,snlp_m,snlp_p
      type(SNLP),pointer                                 :: s
      integer                                            :: i,j
      integer(kind=c_int)                                :: exit_code

      type(MEF90DefMechPlasticityCtx),target             :: PlasticityCtx
      type(c_ptr)                                        :: snlp_ctx


      type(VecScatter)                                   :: ScatterSecToVecCellMatS,ScatterSecToVec,ScatterSecToVecScal
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
      PetscReal                                          :: damageElem

      !PetscReal,Dimension(:),Pointer                     :: damageLoc
      !Type(MEF90Element_Type)                            :: elemDamageType

      Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMMatS,'default',plasticStrainSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(plasticStrainSec,plasticStrainOldSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(plasticStrainSec,inelasticStrainSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMMatS,plasticStrainSec,ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(plasticStrainSec,ScatterSecToVecCellMatS,SCATTER_REVERSE,plasticStrain,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(plasticStrainOldSec,ScatterSecToVecCellMatS,SCATTER_REVERSE,plasticStrainOld,ierr);CHKERRQ(ierr)
   
      Call DMMeshGetSectionReal(MEF90DefMechCtx%DMVect,'default',xSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMVect,xSec,ScatterSecToVec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(xSec,ScatterSecToVec,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr)         

      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call DMMeshGetSectionReal(MEF90DefMechCtx%DMScal,'default',temperatureSec,ierr);CHKERRQ(ierr)
         Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMScal,temperatureSec,ScatterSecToVecScal,ierr);CHKERRQ(ierr)
         Call SectionRealToVec(temperatureSec,ScatterSecToVecScal,SCATTER_REVERSE,MEF90DefMechCtx%temperature,ierr);CHKERRQ(ierr)          
      Else
         temperatureSec%v = 0
      End If

      If (Associated(MEF90DefMechCtx%damage)) Then
         If (Associated(MEF90DefMechCtx%temperature)) Then
            Call SectionRealDuplicate(temperatureSec,damageSec,ierr);CHKERRQ(ierr)
         Else
            Call DMMeshGetSectionReal(MEF90DefMechCtx%DMScal,'default',damageSec,ierr);CHKERRQ(ierr)
            Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMScal,damageSec,ScatterSecToVecScal,ierr);CHKERRQ(ierr)
         End If
         Call SectionRealToVec(damageSec,ScatterSecToVecScal,SCATTER_REVERSE,MEF90DefMechCtx%damage,ierr);CHKERRQ(ierr)          
      Else
         damageSec%v = 0
      End If
!write(*,*)'Damage',damageSec
!   cellSetOptions%damageType

!write(*,*)'Damage',cellSetOptions%damageType


      Call DMMeshGetDimension(MEF90DefMechCtx%DM,dim,ierr);CHKERRQ(ierr)
      Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
         
         !!GET DAMAGE TYPE
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Select Case (cellSetOptions%damageType)
            Case (MEF90DefMech_damageTypeAT1,MEF90DefMech_damageTypeAT2)
            PlasticityCtx%ValueOfk=0
            Case (MEF90DefMech_damageTypeATk)
            PlasticityCtx%ValueOfk=matpropSet%k_for_ATk
         End Select


         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(ierr)
         elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
      
         elemScalType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)

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

               case(MEF90DefMech_plasticityTypeTresca)
                  snlp_Dfhg = c_null_funptr
                  snlp_fhg  = c_funloc(FHG_TRESCA)
                  snlp_n    = SIZEOFMEF90_MATS
                  snlp_m    = 1
                  snlp_p    = 2*SIZEOFMEF90_MATS
                  snlp_ctx  = c_loc(PlasticityCtx)

               case(MEF90DefMech_plasticityTypeNONE)
                  return

               case default
                  Print*,__FUNCT__,': Unimplemented plasticity Type',cellSetOptions%PlasticityType
                  STOP 

            End select
   
            !! Remplissage du Ctx

            PlasticityCtx%YieldStress = matpropSet%YieldStress
            PlasticityCtx%HookesLaw = matpropSet%HookesLaw

            Call SNLPNew(s,snlp_n,snlp_m,snlp_p,snlp_fhg,snlp_Dfhg,snlp_ctx)

            QuadratureOrder = 2 * (elemDisplacementType%order - 1)
            Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
            Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
            Call MEF90InelasticStrainSet(inelasticStrainSec,xSec,temperatureSec,MEF90DefMechCtx%DMVect,MEF90DefMechCtx%DMScal,setIS,matpropSet%LinearThermalExpansion, &
                                         elemDisplacement,elemDisplacementType,elemScal,elemScalType,ierr)
            
            !Call MEF90DamageSet(damageSec,xSec,MEF90DefMechCtx%DMVect,MEF90DefMechCtx%DMScal,setIS, &
            !                    elemDisplacement,elemDisplacementType,elemScal,elemScalType,ierr)

            Call MEF90Element_Destroy(elemDisplacement,ierr)
            Call MEF90Element_Destroy(elemScal,ierr)





            Do cell = 1,size(cellID)
               !! actualiser le ctx (  HookesLaw ,InelasticStrainSec, plasticStrainStrainSec, plasticStrainOldSec  )
               Call SectionRealRestrict(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
               Call SectionRealRestrict(plasticStrainOldSec,cellID(cell),plasticStrainOldLoc,ierr);CHKERRQ(ierr)
               Call SectionRealRestrict(inelasticStrainSec,cellID(cell),inelasticStrainLoc,ierr);CHKERRQ(ierr)
               !Call SectionRealRestrict(damageSec,cellID(cell),damageLoc,ierr);CHKERRQ(ierr)
               

               !!! Damage on one element is not the good way to do that, HAVE TO CHANGE THAT
               Call MEF90Element_Create(MEF90DefMechCtx%DMVect,setIS,elemDisplacement,QuadratureOrder,CellSetOptions%elemTypeShortIDDisplacement,ierr);CHKERRQ(ierr)
               Call MEF90Element_Create(MEF90DefMechCtx%DMScal,setIS,elemScal,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
               !Call MEF90InelasticStrainSet(inelasticStrainSec,xSec,temperatureSec,MEF90DefMechCtx%DMVect,MEF90DefMechCtx%DMScal,setIS,matpropSet%LinearThermalExpansion, &
               !                             elemDisplacement,elemDisplacementType,elemScal,elemScalType,ierr)
            
               Call MEF90DamageSet(damageSec,damageElem,xSec,MEF90DefMechCtx%DMVect,MEF90DefMechCtx%DMScal,setIS,cell, &
                                elemDisplacement,elemDisplacementType,elemScal,elemScalType,ierr)

               Call MEF90Element_Destroy(elemDisplacement,ierr)
               Call MEF90Element_Destroy(elemScal,ierr)


               PlasticityCtx%Damage = damageElem 
!write(*,*)'DamageCtx',PlasticityCtx%Damage


               PlasticityCtx%PlasticStrainOld = plasticStrainOldLoc
               PlasticityCtx%InelasticStrain = InelasticStrainLoc
            
               !!! This is just for testing
               s%show_progress = 0
      
!write(*,*) 'Plastic Strain Old step: ', PlasticStrainOldLoc
!write(*,*) 'Plastic Strain before:   ', PlasticStrainLoc
!write(*,*) 'Inelastic Strain:        ', InelasticStrainLoc
               !!! This is a bit dangerous:
               !!! If PetscReal is not the same as c_double, this call will fail




               !!! Brittle in traction, Ductile in compression
               Select Case(cellSetOptions%unilateralContactType)
               Case (MEF90DefMech_unilateralContactTypeBrittleDuctile)
                     If (Trace(PlasticityCtx%InelasticStrain) > 0.0_Kr ) Then
                        plasticStrainLoc = plasticStrainOldLoc
                     Else 
                        exit_code = SNLPL1SQP(s,plasticStrainLoc)
                     End if
               Case default
                  exit_code = SNLPL1SQP(s,plasticStrainLoc)
               End Select



!write(*,*) 'Plastic Strain after:  ', PlasticStrainLoc
               Call SectionRealRestore(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
               Call SectionRealRestore(plasticStrainOldSec,cellID(cell),plasticStrainOldLoc,ierr);CHKERRQ(ierr)
               Call SectionRealRestore(inelasticStrainSec,cellID(cell),inelasticStrainLoc,ierr);CHKERRQ(ierr)
               Call SectionRealRestore(damageSec,cellID(cell),damageLoc,ierr);CHKERRQ(ierr)
      
            End Do !cell
            call SNLPDelete(s)
         End If ! set 
      End Do !! set

      Call SectionRealToVec(plasticStrainSec,ScatterSecToVecCellMatS,SCATTER_FORWARD,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)


      Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(plasticStrainOldSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(inelasticStrainSec,ierr);CHKERRQ(ierr)

      Call VecScatterDestroy(ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)

      If (Associated(MEF90DefMechCtx%damage)) Then
         Call SectionRealDestroy(damageSec,ierr);CHKERRQ(ierr)
         Call VecScatterDestroy(ScatterSecToVecScal,ierr);CHKERRQ(ierr)
      End If

      If (Associated(MEF90DefMechCtx%temperature)) Then
         Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
         Call VecScatterDestroy(ScatterSecToVecScal,ierr);CHKERRQ(ierr)
      End If
#else
      write(*,*) 'This example needs SNLP'
#endif
      End Subroutine MEF90DefMechPlasticStrainUpdate
End Module MEF90_APPEND(m_MEF90_DefMechPlasticity,MEF90_DIM)D
