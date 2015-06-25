#include "../MEF90/mef90.inc"
Module m_MEF90_DefMechPlasticity_Type
#include "finclude/petscdef.h"
   Use m_MEF90
   Implicit none
   !private
   !public MEF90DefMechPlasticStrainUpdate

   Type MEF90_PlasticityProjectionCtx2D
      Type(MEF90HookesLaw2D)     :: HookesLaw
      Type(MatS3D)               :: PlasticStrainPrevious
   End Type

#ifdef REMOVETHIS
Contains
#undef __FUNCT__
#define __FUNCT__ "fhg_tresca2D"
!!!
!!!  
!!!  fhg_tresca2D:
!!!  
!!!  (c) 2015 Erwan Tanne 
!!!
   Subroutine fhg_tresca2D(,ierr)
      PetscErrorCode,Intent(OUT)                         :: ierr

      
   End Subroutine fhg_tresca2D


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechPlasticStrainUpdate"
!!!
!!!  
!!!  MEF90DefMechPlasticStrainUpdate:
!!!  
!!!  (c) 2015 Erwan Tanne
!!!
   Subroutine MEF90DefMechPlasticStrainUpdate(MEF90DefMechCtx,PlasticStrainOld,ierr)
      use,intrinsic :: iso_c_binding
#ifdef MEF90_HAVE_SNLP   
      use snlpF90
#endif
      PetscErrorCode,Intent(OUT)                         :: ierr
      Type(Vec)                                          :: PlasticStrainOld !!! p_{i-1}
      !!! p_{i}^k is MEF90DefMechCtx%PlasticStrain

      type(c_func_ptr)                                   :: fhg,Dfhg
      Type(SectionReal)                                  :: plasticStrainSec,plasticStrainOldSec      
      PetscReal,Dimension(:),Pointer                     :: plasticStrainLoc
      Type(SNLP),pointer                                 :: SNLPPlasticStrain
      Type(MEF90_PlasticityProjectionCtx2D)              :: ctx2D
      Type(MEF90_PlasticityProjectionCtx3D)              :: ctx3D
      Type(MEF90MatProp2D_Type),pointer                  :: matProp2D
      Type(MEF90MatProp3D_Type),pointer                  :: matProp3D

      fhg  = c_null_functionpointer
      Dfhg = c_null_functionpointer

#ifdef MEF90_HAVE_SNLP   
      Call DMMeshGetSectionReal(MEF90DefMechCtx%CellDMMatS,'default',plasticStrainSec,ierr);CHKERRQ(ierr)
      Call SectionRealDuplicate(plasticStrainSec,plasticStrainOldSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%CellDMMatS,plasticStrainSec,ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(plasticStrainSec,ScatterSecToVecCellMatS,SCATTER_REVERSE,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)          
      Call SectionRealToVec(plasticStrainOldSec,ScatterSecToVecCellMatS,SCATTER_REVERSE,plasticStrainOld,ierr);CHKERRQ(ierr)          

      !!! get IS for cell sets
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
            case(tresca)
               if (dim == 2) then
                  snlp_fhg  => c_fun_loc(fhg_tresca2D)
                  snlp_Dfhg =  c_null_functionpointer
                  snlp_n = 
                  snlp_m = 
                  snlp_p = 
                  ctx2D%HookesLaw = MatProp2D%HookesLaw
                  :
                  :         
               end if
         end select
         !Setup ctx from matpropSet
         if (dim == 2) then
            Call SNLPNew(s,fhg,Dfhg,ctx2D)
         else
            Call SNLPNew(s,fhg,Dfhg,ctx3D)
         end if
         
         Do cell = 1,size(cellID)
            Call SectionRealRestrict(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
            Call SectionRealRestrict(plasticStrainOldSec,cellID(cedll),plasticStrainOldLoc,ierr);CHKERRQ(ierr)
            ctx%PlasticStrainOld = plasticStrainOldLoc
            call SNLPSolve(plasticStrainLoc,fhg,Dfhg,ctx)
            Call SectionRealRestore(plasticStrainSec,cellID(cell),plasticStrainLoc,ierr);CHKERRQ(ierr)
            Call SectionRealRestore(plasticStrainOldSec,cellID(cell),plasticStrainOldLoc,ierr);CHKERRQ(ierr)
         End Do !cell
         call SNLPDelete(s)
      End Do !set
      Call SectionRealToVec(plasticStrainSec,ScatterSecToVecCellMatS,SCATTER_FORWARD,MEF90DefMechCtx%plasticStrain,ierr);CHKERRQ(ierr)          
      Call SectionRealDestroy(plasticStrainSec,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(plasticStrainOldSec,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVecCellMatS,ierr);CHKERRQ(ierr)    
#else
      Write(*,*) ' MEF90DefMechPlasticStrainUpdate requires SNLP'
#endif
   End Subroutine MEF90DefMechPlasticStrainUpdate
#endif

End Module m_MEF90_DefMechPlasticity_Type