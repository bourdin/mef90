#include "../MEF90/mef90.inc"
Module m_MEF90_DefMechPlasticity_Type_mod
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
      Type(MatS2D)            :: InelasticStrain
      Type(MatS2D)            :: PlasticStrainOld
      Type(MatS2D)            :: PlasticStrain
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
      type(MatS2D)                  :: x3D
      x3D%XX = x(1)
      x3D%YY = x(2)
      x3D%XY = x(3)
      !!! This is the fortran equivalent of casting ctx into a c_ptr
      call c_f_pointer(myctx,myctx_ptr)
      f(1) = ( (myctx_ptr%HookesLaw * (x3D-myctx_ptr%PlasticStrainOld)) .DotP. (x3D-myctx_ptr%PlasticStrainOld) ) /2.
      h(1) = Trace(x3D)
      g(1) = sqrt( 2.0*trace(  deviatoricPart(myctx_ptr%HookesLaw*(myctx_ptr%InelasticStrain-x3D))  *  deviatoricPart(myctx_ptr%HookesLaw*(myctx_ptr%InelasticStrain-x3D)) )) - myctx_ptr%YieldStress
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
   Type(DM),target,Intent(IN)                         :: Mesh

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
   !Adding type needs for inelastic strain
   PetscReal,Dimension(:),Pointer                     :: InelasticStrainLoc
   Type(SectionReal)                                  :: InelasticStrainSec
   Type(SectionReal)                                  :: xSec,temperatureSec
   Type(MEF90_ELEMENT_ELAST),Dimension(:),Pointer     :: elemDisplacement
   Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elemScal
   Type(MEF90Element_Type)                            :: elemScalType

   Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMMatS,'default',plasticStrainSec,ierr);CHKERRQ(ierr)
   Call SectionRealDuplicate(plasticStrainSec,plasticStrainOldSec,ierr);CHKERRQ(ierr)
   Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMMatS,plasticStrainSec,ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)
   Call SectionRealToVec(plasticStrainSec,ScatterSecToVecCellMatS,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)
!   Call SectionRealToVec(plasticStrainOldSec,ScatterSecToVecCellMatS,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrainOld,ierr);CHKERRQ(ierr)

   Call DMMeshGetDimension(mesh,dim,ierr);CHKERRQ(ierr)
   Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
   Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

   Do set = 1,size(setID)
!      Call PetscBagGetDataMEF90MatProp(MEF90DefMechCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
      Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
      Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
      elemDisplacementType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDisplacement)
      Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)
      If ((Size(cellID) > 0) .AND. (elemDisplacementType%coDim == 0)) Then

!!       Call proper local assembly depending on the type of damage law
         Select Case (cellSetOptions%plasticityType)
            case(MEF90DefMech_plasticityTypeVonMises)
!!          ! Faire bloucle pour dim2 et dim 3
            ctx_ptr%HookesLaw = matprop2D%HookesLaw
            snlp_Dfhg = c_null_funptr
            snlp_fhg  = c_funloc(fhg_VonMises2D)
            snlp_n    = 3
            snlp_m    = 1
            snlp_p    = 1
            snlp_ctx  = c_loc(ctx_ptr)
         End select

         !! Remplissage du Ctx

         ctx_ptr%YieldStress      = 1.0_Kr
         ctx_ptr%PlasticStrain    = 0.0_Kr
         ctx_ptr%InelasticStrain  = 0.0_Kr
         ctx_ptr%PlasticStrainOld = 0.0_Kr

         allocate(x(snlp_n))
         x = ctx_ptr%PlasticStrain

         if (dim == 2) then
            Call SNLPNew(s,snlp_n,snlp_m,snlp_p,snlp_fhg,snlp_Dfhg,snlp_ctx)
         else
            Call SNLPNew(s,snlp_n,snlp_m,snlp_p,snlp_fhg,snlp_Dfhg,snlp_ctx)
         End If
      End If

      Do cell = 1,size(cellID)

         !! actualiser le ctx (  HookesLaw ,InelasticStrainSec, plasticStrainStrainSec, plasticStrainOldSec  )

         Call SectionRealRestrict(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
         Call SectionRealRestrict(plasticStrainOldSec,cellID(cell),plasticStrainOldLoc,ierr);CHKERRQ(ierr)
         !! need inelastic strain e(u)= \nabla_s u - e_therm
         !Call SectionRealRestrict(InelasticStrainSec,cellID(cell),InelasticStrainLoc,ierr);CHKERRQ(ierr)
         !Call MEF90InelasticStrainCell(InelasticStrainSec,xSec,temperatureSec,MEF90DefMechCtx%DMVect,MEF90DefMechCtx%DMScal,setIS,cell &
         !                                 ,matpropSet%LinearThermalExpansion,elemDisplacement,elemDisplacementType,elemScal,elemScalType,ierr)

         !problem hookes law
         ctx_ptr%PlasticStrainOld   =   plasticStrainOldLoc
         ctx_ptr%PlasticStrain      =   plasticStrainLoc
         ctx_ptr%InelasticStrain%XX =   1.0_Kr
         write(*,*) 'HookeLaw: ',ctx_ptr%HookesLaw

         s%show_progress = 1
         exit_code = SNLPL1SQP(s,x)

         write(*,*) 'x: ',x

         Call SectionRealRestore(InelasticStrainSec,cellID(cell),InelasticStrainLoc,ierr);CHKERRQ(ierr)
         Call SectionRealRestore(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
!         Call SectionRealRestore(plasticStrainOldSec,cellID(cell),plasticStrainOldLoc,ierr);CHKERRQ(ierr)
      End Do !cell

      call SNLPDelete(s)
      DeAllocate(x)
   End Do !! set

   Call SectionRealToVec(plasticStrainSec,ScatterSecToVecCellMatS,SCATTER_FORWARD,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
!   Call SectionRealDestroy(plasticStrainOldSec,ierr);CHKERRQ(ierr)
   Call VecScatterDestroy(ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)

#else
   write(*,*) 'This example needs SNLP'
#endif
   End Subroutine MEF90DefMechPlasticStrainUpdate
End module m_MEF90_DefMechPlasticity_Type_mod










Program TestInelasticStrainCell
#include "../MEF90/mef90.inc"
#include <finclude/petscdef.h>
   Use petsc
   Use m_MEF90
   !Use m_vDef
   Use SNLPF90
   Use m_MEF90_DefMechCtx
   Use m_MEF90_DefMech
   Use m_MEF90_DefMechPlasticity_Type_mod
   Implicit NONE


   PetscErrorCode                      :: ierr
   Type(DM),target                     :: Mesh
   Character(len=MEF90_MXSTRLEN)       :: IOBuffer
   PetscInt                            :: dim
   Type(Vec)                           :: xVec
   PetscReal,Dimension(:,:),Pointer    :: coordPtr
   PetscReal,Dimension(:),Pointer      :: xPtr
   PetscBool                           :: flg
   Integer                             :: i,j,k,v

   PetscReal                           :: scal,sol
   Type(Vect2D)                        :: v2D
   Type(Vect3D)                        :: v3D
   Type(MatS2D)                        :: M2D
   Type(MatS3D)                        :: M3D

   Integer                             :: numVertex
   Integer                             :: QuadOrderMax,QuadOrder
   PetscInt                            :: localSize,globalSize
   Type(SectionReal)                   :: defaultSection


   Type(MEF90Ctx_Type),target                         :: MEF90Ctx
   Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions
   Type(MEF90CtxGlobalOptions_Type),Parameter         :: MEF90DefaultGlobalOptions = MEF90CtxGlobalOptions_Type( &
                                                         1,                             & ! verbose
                                                         PETSC_FALSE,                   & ! validate
                                                         MEF90TimeInterpolation_linear, & ! timeInterpolation
                                                         0.0_Kr,                        & ! timeMin
                                                         1.0_Kr,                        & ! timeMax
                                                         11,                            & ! timeNumStep
                                                         MEF90FileFormat_EXOSingle)       ! fileFormat

   !!! Defect mechanics contexts
   Type(MEF90DefMechCtx_Type)                         :: MEF90DefMechCtx
   Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
   Type(MEF90DefMechGlobalOptions_Type),Parameter     :: MEF90DefMechDefaultGlobalOptions2D = MEF90DefMechGlobalOptions_Type( &
                                                         MEF90DefMech_ModeQuasiStatic, & ! mode
                                                         PETSC_TRUE,              & ! disp_addNullSpace
                                                         3,                       & ! DisplacementOffset
                                                         2,                       & ! DamageOffset
                                                         3,                       & ! boundaryDisplacementOffset
                                                         0,                       & ! boundaryDamageOffset
                                                         1,                       & ! temperatureOffset
                                                         4,                       & ! ForceOffset
                                                         3,                       & ! pressureForceOffset
                                                         0,                       & ! plasticStrainOffset
                                                         6,                       & ! StressOffset
                                                         MEF90Scaling_Linear,     & ! boundaryDisplacementScaling
                                                         MEF90Scaling_CST,        & ! boundaryDamageScaling
                                                         MEF90Scaling_Linear,     & ! ForceScaling
                                                         MEF90Scaling_Linear,     & ! pressureForceScaling
                                                         1e-4,                    & ! damage_atol
                                                         1000,                    & ! maxit
                                                         0.,                      & ! irrevThres 
                                                         MEF90DefMech_BTTypeNULL, & ! BTType
                                                         -1,                      & ! BTInt
                                                         -1,                      & ! BTScope
                                                         1.0e-2)                    ! BTTol
   Type(MEF90DefMechGlobalOptions_Type),Parameter     :: MEF90DefMechDefaultGlobalOptions3D = MEF90DefMechGlobalOptions_Type( &
                                                         MEF90DefMech_ModeQuasiStatic, & ! mode
                                                         PETSC_TRUE,              & ! disp_addNullSpace
                                                         3,                       & ! DisplacementOffset
                                                         2,                       & ! DamageOffset
                                                         3,                       & ! boundaryDisplacementOffset
                                                         0,                       & ! boundaryDamageOffset
                                                         1,                       & ! temperatureOffset
                                                         4,                       & ! ForceOffset
                                                         3,                       & ! pressureForceOffset
                                                         0,                       & ! plasticStrainOffset
                                                         7,                       & ! StressOffset
                                                         MEF90Scaling_Linear,     & ! boundaryDisplacementScaling
                                                         MEF90Scaling_CST,        & ! boundaryDamageScaling
                                                         MEF90Scaling_Linear,     & ! ForceScaling
                                                         MEF90Scaling_Linear,     & ! pressureForceScaling
                                                         1e-4,                    & ! damage_atol
                                                         1000,                    & ! maxit
                                                         0.,                      & ! irrevThres 
                                                         MEF90DefMech_BTTypeNULL, & ! BTType
                                                         -1,                      & ! BTInt
                                                         -1,                      & ! BTScope
                                                         1.0e-2)                    ! BTTol

   Type(MEF90DefMechCellSetOptions_Type),Parameter    :: MEF90DefMechDefaultCellSetOptions = MEF90DefMechCellSetOptions_Type( &
                                                         -1,                                      & ! elemTypeShortIDDispl will be overriden
                                                         -1,                                      & ! elemTypeShortIDDamage will be overriden
                                                         [0.0_Kr,0.0_Kr,0.0_Kr],                  & ! force
                                                         0.0_Kr,                                  & ! pressureForce
                                                         MEF90DefMech_damageTypeAT1,              & ! damageType
                                                         MEF90DefMech_plasticityTypeVonMises,     & ! plasticityType
                                                         MEF90DefMech_unilateralContactTypeNone,  & ! unilateralContactType
                                                         [PETSC_FALSE,PETSC_FALSE,PETSC_FALSE],   & ! Has Displacement BC
                                                         [0.0_Kr,0.0_Kr,0.0_Kr],                  & ! boundary Displacement
                                                         PETSC_FALSE,                             & ! Has Damage BC
                                                         0._Kr)                                     ! Boundary Damage
   Type(MEF90DefMechVertexSetOptions_Type),Parameter  :: MEF90DefMechDefaultVertexSetOptions = MEF90DefMechVertexSetOptions_Type( &
                                                         [PETSC_FALSE,PETSC_FALSE,PETSC_FALSE],   & ! Has Displacement BC
                                                         [0.0_Kr,0.0_Kr,0.0_Kr],                  & ! boundary Displacement
                                                         PETSC_FALSE,                             & ! Has Damage BC
                                                         0.0_Kr)                                    ! boundary Damage




   !!! Initialize MEF90
   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   Call MEF90Initialize(ierr)
   

   !!! Get all MEF90-wide options
   Call MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90DefaultGlobalOptions,ierr);CHKERRQ(ierr)
   Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)

   !!! Get DM from mesh
   Call MEF90CtxGetDMMeshEXO(MEF90Ctx,Mesh,ierr);CHKERRQ(ierr)
   Call DMMeshGetDimension(Mesh,dim,ierr);CHKERRQ(ierr)
   Call DMMeshSetMaxDof(Mesh,dim,ierr);CHKERRQ(ierr) 
   Call DMSetBlockSize(Mesh,dim,ierr);CHKERRQ(ierr)

   !!! Create DefMech context, get all DefMech options
   Call MEF90DefMechCtxCreate(MEF90DefMechCtx,Mesh,MEF90Ctx,ierr);CHKERRQ(ierr)
   If (dim == 2) Then
      Call MEF90DefMechCtxSetFromOptions(MEF90DefMechCtx,PETSC_NULL_CHARACTER,MEF90DefMechDefaultGlobalOptions2D, &
                                         MEF90DefMechDefaultCellSetOptions,MEF90DefMechDefaultVertexSetOptions,ierr)
   Else
      Call MEF90DefMechCtxSetFromOptions(MEF90DefMechCtx,PETSC_NULL_CHARACTER,MEF90DefMechDefaultGlobalOptions3D, &
                                         MEF90DefMechDefaultCellSetOptions,MEF90DefMechDefaultVertexSetOptions,ierr)
   End If
   Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)

   !!! Create sections, vectors, and solvers for DefMech Context
   Call MEF90DefMechCtxSetSections(MEF90DefMechCtx,ierr)
   Call MEF90DefMechCtxCreateVectors(MEF90DefMechCtx,ierr)

   !!! test de la function DefMechPlasticStrainUpdate

   Call MEF90DefMechPlasticStrainUpdate(MEF90DefMechCtx,Mesh,ierr)





End Program TestInelasticStrainCell