#include "../MEF90/mef90.inc"
Module m_MEF90_DefMechPlasticity
#include "finclude/petscdef.h"
   use m_MEF90
   use m_MEF90_DefMechCtx
   implicit NONE
   private
   public MEF90DefMechPlasticStrainUpdate

   !!! note that this type is NOT C interoperable, which is not an issue, since we only
   !!! need SNLP to carry its address

   type :: ctx 
      Type(MEF90HookesLaw2D)  :: HookesLaw
      real(Kind = Kr)         :: YieldStress
      Type(MatS2D)            :: InelasticStrain
      Type(MatS2D)            :: PlasticStrainOld
      !Type(MatS2D)            :: PlasticStrain
   end type ctx


contains
#undef __FUNCT__
#define __FUNCT__ "fhg_VonMises2D"

!!!
!!!  
!!!  fhg:
!!!  
!!!  (c) 2015 Erwan Tanne : erwan.tanne@gmail.com
!!!
!!!

   subroutine fhg_VonMises2D(x,f,h,g,myctx) bind(c)
      use,intrinsic :: iso_c_binding
      use m_MEF90

      real(kind=c_double)           :: x(*)
      real(kind=c_double)           :: f(*)
      real(kind=c_double)           :: h(*)
      real(kind=c_double)           :: g(*)
      type(c_ptr),intent(in),value  :: myctx

      type(ctx),pointer             :: myctx_ptr
      type(MatS2D)                  :: x2D

      x2D = x(1:3)
      !!! This is the fortran equivalent of casting ctx into a c_ptr
      call c_f_pointer(myctx,myctx_ptr)
      write(*,*) 'FHG: x:              ',x(1:3)
      !write(*,*) 'FHG: HookesLaw:      ',myctx_ptr%HookesLaw
      !write(*,*) 'FHG: x2D:            ',x2D
      !write(*,*) 'FHG: PlasticStrainOld', myctx_ptr%PlasticStrainOld
      f(1) = ( (myctx_ptr%HookesLaw * (x2D-myctx_ptr%PlasticStrainOld)) .DotP. (x2D-myctx_ptr%PlasticStrainOld) ) /2.
      h(1) = Trace(x2D)
      g(1) = sqrt( 2.0*trace(  deviatoricPart(myctx_ptr%HookesLaw*(myctx_ptr%InelasticStrain-x2D))  *  deviatoricPart(myctx_ptr%HookesLaw*(myctx_ptr%InelasticStrain-x2D)) )) - myctx_ptr%YieldStress
      write(*,*) 'FHG: f,h,g           ', f(1),h(1),g(1)
      !write(*,*)
   end subroutine fhg_VonMises2D


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechPlasticStrainUpdate"
!!!
!!!  
!!!  MEF90DefMechPlasticStrainUpdate:
!!!  
!!!  (c) 2015 Erwan Tanne : erwan.tanne@gmail.com
!!!

   Subroutine MEF90DefMechPlasticStrainUpdate(MEF90DefMechCtx,PlasticStrain,PlasticStrainOld,ierr)
   use,intrinsic :: iso_c_binding
#ifdef MEF90_HAVE_SNLP
   use SNLPF90

   PetscErrorCode,Intent(OUT)                         :: ierr
   Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
   Type(Vec)                                          :: PlasticStrain,PlasticStrainOld

   Type(DM)                                           :: Mesh
   Type(SectionReal)                                  :: plasticStrainSec,plasticStrainOldSec
   PetscReal,Dimension(:),Pointer                     :: plasticStrainLoc,plasticStrainOldLoc
   type(c_funptr)                                     :: snlp_fhg,snlp_Dfhg
   integer(kind=c_int)                                :: snlp_n,snlp_m,snlp_p
   type(c_ptr)                                        :: snlp_ctx
   type(SNLP),pointer                                 :: s
   integer                                            :: i,j
   integer(kind=c_int)                                :: exit_code
   type(ctx),target                                   :: ctx_ptr
   type(VecScatter)                                   :: ScatterSecToVecCellMatS
   PetscInt                                           :: dim,set,cell
   Type(IS)                                           :: cellSetGlobalIS,setIS
   PetscInt,dimension(:),Pointer                      :: setID,cellID
   Type(MEF90DefMechCellSetOptions_Type),pointer      :: cellSetOptions
   !Type(MEF90_MATPROP),Pointer                        :: matpropSet
   Type(MEF90Element_Type)                            :: elemDisplacementType
   Type(MEF90MatProp2D_Type),pointer                  :: matProp2D
   !Adding type needs for inelastic strain
   PetscReal,Dimension(:),Pointer                     :: InelasticStrainLoc
   Type(SectionReal)                                  :: xSec,temperatureSec
   Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
   Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemScal
   Type(MEF90Element_Type)                            :: elemScalType

   Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMMatS,'default',plasticStrainSec,ierr);CHKERRQ(ierr)
   Call SectionRealDuplicate(plasticStrainSec,plasticStrainOldSec,ierr);CHKERRQ(ierr)
   Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMMatS,plasticStrainSec,ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)
   Call SectionRealToVec(plasticStrainSec,ScatterSecToVecCellMatS,SCATTER_REVERSE,plasticStrain,ierr);CHKERRQ(ierr)
   Call SectionRealToVec(plasticStrainOldSec,ScatterSecToVecCellMatS,SCATTER_REVERSE,plasticStrainOld,ierr);CHKERRQ(ierr)
   
   Call DMMeshGetDimension(MEF90DefMechCtx%DM,dim,ierr);CHKERRQ(ierr)
   Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
   Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

   Do set = 1,size(setID)
      Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matprop2D,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
      elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
      
      Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
      If ((Size(cellID) > 0) .AND. (elemDisplacementType%coDim == 0)) Then
         !!! Call proper local assembly depending on the type of damage law
         Select Case (cellSetOptions%plasticityType)
            case(MEF90DefMech_plasticityTypeVonMises)
               ! Faire bloucle pour dim2 et dim 3
               ctx_ptr%HookesLaw = MatProp2D%HookesLaw
               snlp_Dfhg = c_null_funptr
               snlp_fhg  = c_funloc(fhg_VonMises2D)
               snlp_n    = 3
               snlp_m    = 1!1
               snlp_p    = 1!1
               snlp_ctx  = c_loc(ctx_ptr)
            case default
               Print*,__FUNCT__,': Unimplemented plasticity Type',cellSetOptions%PlasticityType
               STOP  

         End select
   
         !! Remplissage du Ctx

         ctx_ptr%YieldStress = 1.0_Kr
         ctx_ptr%HookesLaw = matprop2D%HookesLaw
      
         
         !if (dim == 2) then
            Call SNLPNew(s,snlp_n,snlp_m,snlp_p,snlp_fhg,snlp_Dfhg,snlp_ctx)
         !else
         !   Call SNLPNew(s,snlp_n,snlp_m,snlp_p,snlp_fhg,snlp_Dfhg,snlp_ctx)
         !End If

         Do cell = 1,size(cellID)
            !! actualiser le ctx (  HookesLaw ,InelasticStrainSec, plasticStrainStrainSec, plasticStrainOldSec  )
            Call SectionRealRestrict(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
            Call SectionRealRestrict(plasticStrainOldSec,cellID(cell),plasticStrainOldLoc,ierr);CHKERRQ(ierr)
            !! need inelastic strain e(u)= \nabla_s u - e_therm
            !!! You don't need the Section, only the local version.
            !!! Modify your function so that it computes the inelastic strain in a given cell, given local (stress, strain, temp etc) fields
            !Call MEF90InelasticStrainCell(InelasticStrainSec,xSec,temperatureSec,MEF90DefMechCtx%DMVect,MEF90DefMechCtx%DMScal,setIS,cell &
            !                                 ,matpropSet%LinearThermalExpansion,elemDisplacement,elemDisplacementType,elemScal,elemScalType,ierr)
            ctx_ptr%InelasticStrain    =   [0.0_Kr,0.0_Kr,0.0_Kr]

            ctx_ptr%PlasticStrainOld   =   plasticStrainOldLoc
      
            
            !!! This is just for testing
            s%show_progress = 1
            plasticStrainLoc = [0.1_Kr,-1._Kr,0.0_Kr]*100.0_Kr
      

            write(*,*) 'Plastic Strain before: ', PlasticStrainLoc
            !write(*,*) 'Plastic Strain: Old    ', ctx_ptr%PlasticStrainOld
            !write(*,*) 'Inelastic Strain:      ', ctx_ptr%InelasticStrain
            !Write(*,*) 'HookesLaw               ', ctx_ptr%HookesLaw
         
            !!! This is a bit dangerous:
            !!! If PetscReal is not the same as c_double, this call will fail
            exit_code = SNLPL1SQP(s,plasticStrainLoc)
            write(*,*) 'Plastic Strain after:  ', PlasticStrainLoc
            write(*,*) 
            Call SectionRealRestore(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
            Call SectionRealRestore(plasticStrainOldSec,cellID(cell),plasticStrainOldLoc,ierr);CHKERRQ(ierr)
      
         End Do !cell
         call SNLPDelete(s)
      End If ! set 
   End Do !! set

   Call SectionRealToVec(plasticStrainSec,ScatterSecToVecCellMatS,SCATTER_FORWARD,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)

   Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(plasticStrainOldSec,ierr);CHKERRQ(ierr)
   Call VecScatterDestroy(ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)
#else
   write(*,*) 'This example needs SNLP'
#endif
   End Subroutine MEF90DefMechPlasticStrainUpdate

End module m_MEF90_DefMechPlasticity