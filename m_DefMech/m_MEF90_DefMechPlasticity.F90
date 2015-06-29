#include "../MEF90/mef90.inc"
Module m_MEF90_DefMechPlasticity_Type
#include "finclude/petscdef.h"
   use m_MEF90
   Use m_MEF90_DefMechCtx
   implicit NONE
   !private
   !public MEF90DefMechPlasticStrainUpdate

   !!! note that this type is NOT C interoperable, which is not an issue, since we only
   !!! need SNLP to carry its address

   type :: ctx 
      Type(MEF90HookesLaw2D)  :: HookesLaw
      real(Kind = Kr)         :: YieldStress
      Type(MatS2D)            :: Strain
      Type(MatS2D)            :: PlasticStrainOld
      Type(MatS2D)            :: PlasticStrain
   end type ctx


!#ifdef REMOVETHIS

contains


#undef __FUNCT__
#define __FUNCT__ "fhg"

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
      type(MatS2D)                  :: x3D
      x3D%XX = x(1)
      x3D%YY = x(2)
      x3D%XY = x(3)
      !!! This is the fortran equivalent of casting ctx into a c_ptr
      call c_f_pointer(myctx,myctx_ptr)      
      f(1) = ( (myctx_ptr%HookesLaw * (x3D-myctx_ptr%PlasticStrainOld)) .DotP. (x3D-myctx_ptr%PlasticStrainOld) ) /2.
      h(1) = Trace(x3D)
      g(1) = sqrt( 2.0*trace(  deviatoricPart(myctx_ptr%HookesLaw*(myctx_ptr%Strain-x3D))  *  deviatoricPart(myctx_ptr%HookesLaw*(myctx_ptr%Strain-x3D)) )) - myctx_ptr%YieldStress
   end subroutine fhg_VonMises2D


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechPlasticStrainUpdate"
!!!
!!!  
!!!  MEF90DefMechPlasticStrainUpdate:
!!!  
!!!  (c) 2015 Erwan Tanne : erwan.tanne@gmail.com
!!!

   Subroutine MEF90DefMechPlasticStrainUpdate(MEF90DefMechCtx,Mesh,ierr)
   use,intrinsic :: iso_c_binding
#ifdef MEF90_HAVE_SNLP
   use SNLPF90

   PetscErrorCode,Intent(OUT)                         :: ierr
   Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
   Type(DM),target,Intent(IN)                               :: Mesh

   Type(Vec)                                          :: PlasticStrainOld !!! p_{i-1}
   Type(SectionReal)                                  :: plasticStrainSec,plasticStrainOldSec
   PetscReal,Dimension(:),Pointer                     :: plasticStrainLoc,plasticStrainOldLoc
   type(c_funptr)                                     :: snlp_fhg,snlp_Dfhg
   integer(kind=c_int)                                :: snlp_n,snlp_m,snlp_p
   type(c_ptr)                                        :: snlp_ctx
   type(SNLP),pointer                                 :: s
   integer                                            :: i,j
   integer(kind=c_int)                                :: exit_code
   real(kind=c_double),dimension(:),pointer           :: x
   type(ctx),target                                   :: ctx_ptr
   type(VecScatter)                                   :: ScatterSecToVecCellMatS
   PetscInt                                           :: dim,set,cell
   Type(IS)                                           :: cellSetGlobalIS,setIS
   PetscInt,dimension(:),Pointer                      :: setID,cellID
   Type(MEF90DefMechCellSetOptions_Type),pointer      :: cellSetOptions
   Type(MEF90_MATPROP),Pointer                        :: matpropSet
   Type(MEF90Element_Type)                            :: elemDisplacementType
   Type(MEF90MatProp2D_Type),pointer                  :: matProp2D

   Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMMatS,'default',plasticStrainSec,ierr);CHKERRQ(ierr)
   Call SectionRealDuplicate(plasticStrainSec,plasticStrainOldSec,ierr);CHKERRQ(ierr)
   Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMMatS,plasticStrainSec,ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)
   Call SectionRealToVec(plasticStrainSec,ScatterSecToVecCellMatS,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)
   Call SectionRealToVec(plasticStrainOldSec,ScatterSecToVecCellMatS,SCATTER_REVERSE,plasticStrainOld,ierr);CHKERRQ(ierr)


   Call DMMeshGetDimension(mesh,dim,ierr);CHKERRQ(ierr)
   Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
   Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

   Do set = 1,size(setID)
      Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
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
         snlp_m    = 1
         snlp_p    = 1
         snlp_ctx  = c_loc(ctx_ptr)
      End select

      ctx_ptr%YieldStress = 1.0_Kr
      ctx_ptr%PlasticStrain = plasticStrainLoc
      ctx_ptr%Strain = 0.0_Kr
   
   
      ctx_ptr%PlasticStrainOld = 0.0_Kr
      ctx_ptr%PlasticStrain = 0.0_Kr
   
      
      allocate(x(snlp_n))
      x = ctx_ptr%PlasticStrain
      

   
      call SNLPNew(s,snlp_n,snlp_m,snlp_p,snlp_fhg,snlp_Dfhg,snlp_ctx)
      s%show_progress = 1
      
      !exit_code = SNLPL1SQP(s,x)
      write(*,*) 'exit_code: ',exit_code
      write(*,*) 'x:         ',x
      deallocate(x)

      if (dim == 2) then
         Call SNLPNew(s,snlp_n,snlp_m,snlp_p,snlp_fhg,snlp_Dfhg,snlp_ctx)
      else
         Call SNLPNew(s,snlp_n,snlp_m,snlp_p,snlp_fhg,snlp_Dfhg,snlp_ctx)
      End If
      End If

      Do cell = 1,size(cellID)
         Call SectionRealRestrict(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
         Call SectionRealRestrict(plasticStrainOldSec,cellID(cell),plasticStrainOldLoc,ierr);CHKERRQ(ierr)
         ctx_ptr%PlasticStrainOld = plasticStrainOldLoc
         !call SNLPSolve(plasticStrainLoc,snlp_fhg,snlp_Dfhg,snlp_ctx)
         !call SNLPL1SQP(s,x)
         Call SectionRealRestore(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
         Call SectionRealRestore(plasticStrainOldSec,cellID(cell),plasticStrainOldLoc,ierr);CHKERRQ(ierr)
      End Do !cell
      call SNLPDelete(s)
   End Do !! set

   Call SectionRealToVec(plasticStrainSec,ScatterSecToVecCellMatS,SCATTER_FORWARD,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(plasticStrainOldSec,ierr);CHKERRQ(ierr)
   Call VecScatterDestroy(ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)


#else
   write(*,*) 'This example needs SNLP'
#endif
   End Subroutine MEF90DefMechPlasticStrainUpdate
!#endif

End module m_MEF90_DefMechPlasticity_Type