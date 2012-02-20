Module m_VarFilmQS_Post
#include "finclude/petscdef.h"

   Use m_VarFilmQS_Types
   Use m_MEF90
   Use m_Film_Struct

   Implicit NONE   
   
Contains
#undef __FUNC__ 
#define __FUNC__ "FilmEnergy_Assembly"
Subroutine FilmEnergy_Assembly(FilmEnergy, FilmEnergyBlock, AppCtx)     
   PetscReal, Intent(OUT)                       :: FilmEnergy
   PetscReal, Dimension(:), Pointer             :: FilmEnergyBlock
   Type(AppCtx_Type)                            :: AppCtx
   
   PetscInt                                     :: iBlk, iBlkId, iErr
   PetscReal                                    :: MyFilmEnergy
   PetscReal, Dimension(:), Pointer             :: MyFilmEnergyBlock
  

   Call PetscLogStagePush(AppCtx%LogInfo%PostProc_Stage, iErr); CHKERRQ(iErr)
   Call PetscLogEventBegin(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
   
   MyFilmEnergy = 0.0_Kr
   Allocate(MyFilmEnergyBlock(AppCtx%MeshTopology%Num_Elem_Blks_Global))
   MyFilmEnergyBlock = 0.0_Kr
   Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
      If (AppCtx%MyEXO%EBProperty(VarFrac_EBProp_IsBrittle)%Value(iBlkID) /= 0) Then
         Select Case (AppCtx%VarFracSchemeParam%Unilateral)
         Case (VarFrac_Unilateral_NONE)
            Call FilmEnergy_AssemblyBlk_Brittle(MyFilmEnergyBlock(iBlkID), iBlk, AppCtx%U%Sec, AppCtx%Theta%Sec, AppCtx%V%Sec, AppCtx)
         End Select
      Else
         Call FilmEnergy_AssemblyBlk_NonBrittle(MyFilmEnergyBlock(iBlkId), iBlk, AppCtx%U%Sec, AppCtx%Theta%Sec, AppCtx)
      End If
      MyFilmEnergy = MyFilmEnergy + MyFilmEnergyBlock(iBlkId)
   End Do

   Call MPI_AllReduce(MyFilmEnergy, FilmEnergy, 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
   Call MPI_AllReduce(MyFilmEnergyBlock, FilmEnergyBlock, AppCtx%MeshTopology%Num_Elem_Blks_Global, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
   DeAllocate(MyFilmEnergyBlock)
   Call PetscLogEventEnd(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
   Call PetscLogStagePop(iErr); CHKERRQ(iErr)
End Subroutine FilmEnergy_Assembly

#undef __FUNC__ 
#define __FUNC__ "ExtForcesWork_Assembly"
Subroutine ExtForcesWork_Assembly(ExtForcesWork, ExtForcesWorkBlock, AppCtx)     
   PetscReal, Intent(OUT)                       :: ExtForcesWork
   PetscReal, Dimension(:), Pointer             :: ExtForcesWorkBlock
   Type(AppCtx_Type)                            :: AppCtx
   
   PetscInt                                     :: iBlk, iBlkId, iErr
   PetscReal                                    :: MyExtForcesWork
   PetscReal, Dimension(:), Pointer             :: MyExtForcesWorkBlock
  

   Call PetscLogStagePush(AppCtx%LogInfo%PostProc_Stage, iErr); CHKERRQ(iErr)
   Call PetscLogEventBegin(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
   
   MyExtForcesWork = 0.0_Kr
   Allocate(MyExtForcesWorkBlock(AppCtx%MeshTopology%Num_Elem_Blks_Global))
   MyExtForcesWorkBlock = 0.0_Kr
   Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
   End Do

   Call MPI_AllReduce(MyExtForcesWork, ExtForcesWork, 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
   Call MPI_AllReduce(MyExtForcesWorkBlock, ExtForcesWorkBlock, AppCtx%MeshTopology%Num_Elem_Blks_Global, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
   DeAllocate(MyExtForcesWorkBlock)
   Call PetscLogEventEnd(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
   Call PetscLogStagePop(iErr); CHKERRQ(iErr)
End Subroutine ExtForcesWork_Assembly

#undef __FUNC__ 
#define __FUNC__ "FractureEnergy_Assembly"
Subroutine FractureEnergy_Assembly(FractureEnergy, FractureEnergyBlock, AppCtx)     
   PetscReal, Intent(OUT)                       :: FractureEnergy
   PetscReal, Dimension(:), Pointer             :: FractureEnergyBlock
   Type(AppCtx_Type)                            :: AppCtx
   
   PetscInt                                     :: iBlk, iBlkId, iErr
   PetscReal                                    :: MyFractureEnergy
   PetscReal, Dimension(:), Pointer             :: MyFractureEnergyBlock
  

   Call PetscLogStagePush(AppCtx%LogInfo%PostProc_Stage, iErr); CHKERRQ(iErr)
   Call PetscLogEventBegin(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
   
   MyFractureEnergy = 0.0_Kr
   Allocate(MyFractureEnergyBlock(AppCtx%MeshTopology%Num_Elem_Blks_Global))
   Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
      Select Case (AppCtx%VarFracSchemeParam%AtNum)
      Case(1)
         Call FractureEnergy_AssemblyBlk_AT1(MyFractureEnergyBlock(iBlkID), iBlk, AppCtx%V%Sec, AppCtx)
      Case(2)
         Call FractureEnergy_AssemblyBlk_AT2(MyFractureEnergyBlock(iBlkID), iBlk, AppCtx%V%Sec, AppCtx)
      Case Default
       SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP, 'Only AT1 and AT2 are implemented\n', iErr)
      End Select
      MyFractureEnergy = MyFractureEnergy + MyFractureEnergyBlock(iBlkID)
   End Do

   Call MPI_AllReduce(MyFractureEnergy, FractureEnergy, 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
   Call MPI_AllReduce(MyFractureEnergyBlock, FractureEnergyBlock, AppCtx%MeshTopology%Num_Elem_Blks_Global, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
   DeAllocate(MyFractureEnergyBlock)
   Call PetscLogEventEnd(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
   Call PetscLogStagePop(iErr); CHKERRQ(iErr)
End Subroutine FractureEnergy_Assembly

#undef __FUNC__ 
#define __FUNC__ "DelaminationEnergy_Assembly"
Subroutine DelaminationEnergy_Assembly(DelaminationEnergy, DelaminationEnergyBlock, AppCtx)     
	PetscReal, Intent(OUT)                       :: DelaminationEnergy
	PetscReal, Dimension(:), Pointer             :: DelaminationEnergyBlock
	Type(AppCtx_Type)                            :: AppCtx
	
	PetscInt                                     :: iBlk, iBlkId, iErr
	PetscReal                                    :: MyDelaminationEnergy
	PetscReal, Dimension(:), Pointer             :: MyDelaminationEnergyBlock
	
	
	Call PetscLogStagePush(AppCtx%LogInfo%PostProc_Stage, iErr); CHKERRQ(iErr)
	Call PetscLogEventBegin(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
	
	MyDelaminationEnergy = 0.0_Kr
	Allocate(MyDelaminationEnergyBlock(AppCtx%MeshTopology%Num_Elem_Blks_Global))
	Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
		iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
		Call DelaminationEnergy_AssemblyBlk(MyDelaminationEnergyBlock(iBlkID), iBlk, AppCtx%W%Sec, AppCtx)
		MyDelaminationEnergy = MyDelaminationEnergy + MyDelaminationEnergyBlock(iBlkID)
	End Do
	
	Call MPI_AllReduce(MyDelaminationEnergy, DelaminationEnergy, 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
	Call MPI_AllReduce(MyDelaminationEnergyBlock, DelaminationEnergyBlock, AppCtx%MeshTopology%Num_Elem_Blks_Global, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
	DeAllocate(MyDelaminationEnergyBlock)
	Call PetscLogEventEnd(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
	Call PetscLogStagePop(iErr); CHKERRQ(iErr)
End Subroutine DelaminationEnergy_Assembly

#undef __FUNC__ 
#define __FUNC__ "BondingLayerEnergy_Assembly"
Subroutine BondingLayerEnergy_Assembly(BondingLayerEnergy, BondingLayerEnergyBlock, AppCtx)     
   PetscReal, Intent(OUT)                       :: BondingLayerEnergy
   PetscReal, Dimension(:), Pointer             :: BondingLayerEnergyBlock
   Type(AppCtx_Type)                            :: AppCtx
   
   PetscInt                                     :: iBlk, iBlkId, iErr
   PetscReal                                    :: MyBondingLayerEnergy
   PetscReal, Dimension(:), Pointer             :: MyBondingLayerEnergyBlock
   
   
   Call PetscLogStagePush(AppCtx%LogInfo%PostProc_Stage, iErr); CHKERRQ(iErr)
   Call PetscLogEventBegin(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
   
   MyBondingLayerEnergy = 0.0_Kr
   Allocate(MyBondingLayerEnergyBlock(AppCtx%MeshTopology%Num_Elem_Blks_Global))
   Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
      Call BondingLayerEnergy_AssemblyBlk(MyBondingLayerEnergyBlock(iBlkID), iBlk, AppCtx%U%Sec, AppCtx%U0%Sec, AppCtx%W%Sec, AppCtx)
      MyBondingLayerEnergy = MyBondingLayerEnergy + MyBondingLayerEnergyBlock(iBlkID)
   End Do

   BondingLayerEnergy = 0.0_Kr
   
   Call MPI_AllReduce(MyBondingLayerEnergy, BondingLayerEnergy, 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
   Call MPI_AllReduce(MyBondingLayerEnergyBlock, BondingLayerEnergyBlock, AppCtx%MeshTopology%Num_Elem_Blks_Global, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
   DeAllocate(MyBondingLayerEnergyBlock)
   Call PetscLogEventEnd(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
   Call PetscLogStagePop(iErr); CHKERRQ(iErr)
End Subroutine BondingLayerEnergy_Assembly


!!!
!!! Block Assembly Routines
!!!
#undef __FUNC__ 
#define __FUNC__ "FilmEnergy_AssemblyBlk_Brittle"
Subroutine FilmEnergy_AssemblyBlk_Brittle(FilmEnergyBlock, iBlk, U_Sec, Theta_Sec, V_Sec, AppCtx)
   PetscReal, Intent(OUT)                       :: FilmEnergyBlock
   PetscInt                                     :: iBlk
   Type(SectionReal)                            :: U_Sec, Theta_Sec, V_Sec
   Type(AppCtx_Type)                            :: AppCtx

   !!!   _Loc are restriction of fields to local patch (the element)
   !!!   _Elem are local contribution over the element (u_Elem = \sum_i U_Loc(i) BF(i))
   PetscReal, Dimension(:), Pointer             :: U_Loc, V_Loc, Theta_Loc
   Type(Vect2D)                                 :: U_Elem
   Type(Mats2D)                                 :: Strain_Elem, EffectiveStrain_Elem
   PetscReal                                    :: V_Elem, Theta_Elem
   PetscInt                                     :: iE, iEloc, iBlkId, iErr
   PetscInt                                     :: NumDoFScal, NumDoFVect
   PetscInt                                     :: iDoF1, iDoF2, iGauss
   PetscLogDouble                               :: flops       

   flops = 0.0
   FilmEnergyBlock = 0.0_Kr

   NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
   NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF

   Allocate(U_Loc(NumDoFVect))
   Allocate(Theta_Loc(NumDoFScal))
   Allocate(V_Loc(NumDoFScal))

   iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID

   Do_iEloc: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
      iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
      Call SectionRealRestrictClosure(U_Sec,     AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U_Loc,     iErr); CHKERRQ(ierr)
      Call SectionRealRestrictClosure(Theta_Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, Theta_Loc, iErr); CHKERRQ(ierr)
      Call SectionRealRestrictClosure(V_Sec,     AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, V_Loc,     iErr); CHKERRQ(ierr)
      Do_iGauss: Do iGauss = 1, size(AppCtx%ElemVect(iE)%Gauss_C)
         Strain_Elem = 0.0_Kr
         Do iDoF2 = 1, NumDoFVect
            Strain_Elem = Strain_Elem + U_Loc(iDoF2) * AppCtx%ElemVect(iE)%GradS_BF(iDoF2, iGauss)
         End Do
         Theta_Elem = 0.0_Kr
         V_Elem = 0.0_Kr
         Do iDoF2 = 1, NumDoFScal
            Theta_Elem = Theta_Elem + AppCtx%ElemScal(iE)%BF(iDoF2, iGauss) * Theta_Loc(iDoF2)
            V_Elem     = V_Elem     + AppCtx%ElemScal(iE)%BF(iDoF2, iGauss) * V_Loc(iDoF2)
            flops = flops + 2.0
         End Do
         EffectiveStrain_Elem = Strain_Elem - AppCtx%MatProp(iBlkId)%Therm_Exp * Theta_Elem
         FilmEnergyBlock = FilmEnergyblock + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * (V_Elem**2 + AppCtx%VarFracSchemeParam%KEpsilon) * ((AppCtx%MatProp(iBlkId)%Hookes_Law * EffectiveStrain_Elem) .DotP. EffectiveStrain_Elem ) * 0.5_Kr
         flops = flops + 6.0 
      End Do Do_iGauss
   End Do Do_iEloc

   DeAllocate(U_Loc)
   DeAllocate(Theta_Loc)
   DeAllocate(V_Loc)
   Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
End Subroutine FilmEnergy_AssemblyBlk_Brittle

#undef __FUNC__ 
#define __FUNC__ "FilmEnergy_AssemblyBlk_NonBrittle"
Subroutine FilmEnergy_AssemblyBlk_NonBrittle(ElasticEnergyBlock, iBlk, U_Sec, Theta_Sec, AppCtx)
      PetscReal, Intent(OUT)                       :: ElasticEnergyBlock
      PetscInt                                     :: iBlk
      Type(SectionReal)                            :: U_Sec, Theta_Sec
      Type(AppCtx_Type)                            :: AppCtx

      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_ELem = \sum_i U_Loc(i) BF(i))
      PetscReal, Dimension(:), Pointer             :: U_Loc, Theta_Loc
      Type(Mats2D)                                 :: Strain_Elem, EffectiveStrain_Elem
      PetscReal                                    :: Theta_Elem
      PetscInt                                     :: iE, iEloc, iBlkId, iErr
      PetscInt                                     :: NumDoFScal, NumDoFVect
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      PetscLogDouble                               :: flops       
      
      flops = 0.0
      ElasticEnergyBlock = 0.0_Kr

      NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
      NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF

      Allocate(U_Loc(NumDoFVect))
      Allocate(Theta_Loc(NumDoFScal))

      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID

      Do_iEloc: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         Call SectionRealRestrictClosure(U_Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U_Loc, iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(Theta_Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, Theta_Loc, iErr); CHKERRQ(ierr)
         Do_iGauss: Do iGauss = 1, size(AppCtx%ElemVect(iE)%Gauss_C)
            Strain_Elem = 0.0_Kr
            Do iDoF2 = 1, NumDoFVect
               Strain_Elem = Strain_Elem + U_Loc(iDoF2) * AppCtx%ElemVect(iE)%GradS_BF(iDoF2, iGauss)
            End Do
            Theta_Elem = 0.0_Kr
            Do iDoF2 = 1, NumDoFScal
               Theta_Elem = Theta_Elem + AppCtx%ElemScal(iE)%BF(iDoF2, iGauss) * Theta_Loc(iDoF2)
               flops = flops + 2.0
            End Do
            EffectiveStrain_Elem = Strain_Elem - AppCtx%MatProp(iBlkId)%Therm_Exp * Theta_Elem
            
            ElasticEnergyBLock = ElasticEnergyblock + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * ((AppCtx%MatProp(iBlkId)%Hookes_Law * EffectiveStrain_Elem) .DotP. EffectiveStrain_Elem ) * 0.5_Kr 
         End Do Do_iGauss
      End Do Do_iEloc

      DeAllocate(U_Loc)
      DeAllocate(Theta_Loc)
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
   End Subroutine FilmEnergy_AssemblyBlk_NonBrittle

Subroutine DelaminationEnergy_AssemblyBlk(DelaminationEnergyBlock, iBlk, W_Sec, AppCtx)
   PetscReal, Intent(OUT)                       :: DelaminationEnergyBlock
   PetscInt                                     :: iBlk
   Type(SectionReal)                            :: W_Sec
   Type(AppCtx_Type)                            :: AppCtx

   !!!   _Loc are restriction of fields to local patch (the element)
   !!!   _Elem are local contribution over the element (W_Elem = \sum_i W_Loc(i) BF(i))
   PetscReal, Dimension(:), Pointer             :: W_Loc
   PetscReal                                    :: W_Elem, Delamination_Elem
   PetscInt                                     :: iE, iEloc, iBlkId, iErr
   PetscInt                                     :: NumDoFScal
   PetscInt                                     :: iDoF, iGauss
   PetscLogDouble                               :: flops       
   
   flops = 0.0
   DelaminationEnergyBlock = 0.0_Kr

   NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
   
   Allocate(W_Loc(NumDoFScal))
   
   iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
   
   Do_iEloc: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
      iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
      Call SectionRealRestrictClosure(W_Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, W_Loc, iErr); CHKERRQ(ierr)
      
      Do_iGauss: Do iGauss = 1, size(AppCtx%ElemVect(iE)%Gauss_C)
         Delamination_Elem = 0.0_Kr
         Do iDoF = 1, NumDoFScal
            Delamination_Elem = Delamination_Elem + AppCtx%ElemScal(iE)%BF(iDoF, iGauss) * (1.0_Kr - W_Loc(iDoF))
            flops = flops + 1.0
         End Do
         DelaminationEnergyBlock = DelaminationEnergyBlock + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * AppCtx%MatProp(iBlkId)%DelamToughness * Delamination_Elem 
      End Do Do_iGauss
   End Do Do_iEloc

   DeAllocate(W_Loc)
   Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
End Subroutine DelaminationEnergy_AssemblyBlk

Subroutine BondingLayerEnergy_AssemblyBlk(BondingLayerEnergyBlock, iBlk, U_Sec, U0_Sec, W_Sec, AppCtx)
   PetscReal, Intent(OUT)                       :: BondingLayerEnergyBlock
   PetscInt                                     :: iBlk
   Type(SectionReal)                            :: U_Sec, U0_Sec, W_Sec
   Type(AppCtx_Type)                            :: AppCtx

   !!!   _Loc are restriction of fields to local patch (the element)
   !!!   _Elem are local contribution over the element (u_ELem = \sum_i U_Loc(i) BF(i))
   PetscReal, Dimension(:), Pointer             :: U_Loc, U0_Loc, W_Loc
   PetscReal                                    :: Cohesive_Elem, Delamination_Elem
   PetscInt                                     :: iE, iEloc, iBlkId, iErr
   PetscInt                                     :: NumDoFVect, NumDoFScal
   PetscInt                                     :: iDoF, iGauss
   PetscLogDouble                               :: flops       
   Type(Vect2D)               :: U_eff_elem
   
   flops = 0.0
   BondingLayerEnergyBlock = 0.0_Kr
   
   NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
   NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF 
   
   Allocate(U_Loc(NumDoFVect))
   Allocate(W_Loc(NumDoFScal))
   Allocate(U0_Loc(NumDoFVect))
   
   iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
   BondingLayerEnergyBlock = 0.0_Kr
   
   Do_iEloc: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
      iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)

!  REM:  W=1 ==> bonded
!  REM:  V=1 ==> sound

      Call SectionRealRestrictClosure(U_Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U_Loc, iErr); CHKERRQ(ierr)
      Call SectionRealRestrictClosure(W_Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, W_Loc, iErr); CHKERRQ(ierr)
      Call SectionRealRestrictClosure(U0_Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U0_Loc, iErr); CHKERRQ(ierr)

      Do_iGauss: Do iGauss = 1, size(AppCtx%ElemVect(iE)%Gauss_C)
         Delamination_Elem = 0.0_Kr
         U_eff_elem        = 0.0_Kr

         Do iDoF = 1, NumDoFVect
            U_eff_elem=U_eff_elem+AppCtx%ElemVect(iE)%BF(iDoF, iGauss) * (U_loc(iDoF)-U0_loc(iDoF))
         End Do
         
         Do iDoF = 1, NumDoFScal
            Delamination_Elem = Delamination_Elem + AppCtx%ElemScal(iE)%BF(iDoF, iGauss) * (W_Loc(iDoF))
            flops = flops + 1.0
         End Do
         
         BondingLayerEnergyBlock = BondingLayerEnergyBlock + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * AppCtx%MatProp(iBlkId)%Ksubst * (U_eff_elem .DotP. U_eff_elem) * Delamination_Elem * 0.5_Kr 

        End Do Do_iGauss
   End Do Do_iEloc

   DeAllocate(U_Loc)
   DeAllocate(U0_Loc)
   DeAllocate(W_Loc)

   Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
End Subroutine BondingLayerEnergy_AssemblyBlk

Subroutine ExtForcesWork_AssemblyBlk(ExtForcesWorkBlock, iBlk, U_Sec, F_Sec, AppCtx)
      PetscReal, Intent(OUT)                       :: ExtForcesWorkBlock
      PetscInt                                     :: iBlk
      Type(SectionReal)                            :: U_Sec, F_Sec
      Type(AppCtx_Type)                            :: AppCtx

      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_Elem = \sum_i U_Loc(i) BF(i))
      PetscReal, Dimension(:), Pointer             :: U_Loc, F_Loc
      Type(Vect2D)                                 :: U_Elem, F_Elem
      PetscInt                                     :: iE, iEloc, iErr
      PetscInt                                     :: NumDoFVect
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      PetscLogDouble                               :: flops       
      
      flops = 0.0
      ExtForcesWorkBlock = 0.0_Kr

      NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim

      Allocate(U_Loc(NumDoFVect))
      Allocate(F_Loc(NumDoFVect))

!      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID

      Do_iEloc: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         Call SectionRealRestrictClosure(U_Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U_Loc, iErr); CHKERRQ(ierr)
         Call SectionRealRestrictClosure(F_Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, F_Loc, iErr); CHKERRQ(ierr)
         Do_iGauss: Do iGauss = 1, size(AppCtx%ElemVect(iE)%Gauss_C)
            U_Elem = 0.0_Kr
            F_Elem = 0.0_Kr
            Do iDoF2 = 1, NumDoFVect
               U_Elem = U_Elem + U_Loc(iDoF2) * AppCtx%ElemVect(iE)%BF(iDoF2, iGauss)
               F_Elem = F_Elem + F_Loc(iDoF2) * AppCtx%ElemVect(iE)%BF(iDoF2, iGauss)
            End Do
            
            ExtForcesWorkBlock = ExtForcesWorkBlock + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * (U_Elem .DotP. F_Elem) 
         End Do Do_iGauss
      End Do Do_iEloc

      DeAllocate(U_Loc)
      DeAllocate(F_Loc)
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
End Subroutine ExtForcesWork_AssemblyBlk

Subroutine FractureEnergy_AssemblyBlk_AT2(FractureEnergyBlock, iBlk, V_Sec, AppCtx)
   PetscReal, Intent(OUT)                       :: FractureEnergyBlock 
!  PetscReal, Intent(OUT)                       :: DelaminationEnergyBlock 
      PetscInt                                     :: iBlk
      Type(SectionReal)                            :: V_Sec
      Type(AppCtx_Type)                            :: AppCtx

      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_Elem = \sum_i U_Loc(i) BF(i))
      PetscReal, Dimension(:), Pointer             :: V_Loc
      PetscReal                                    :: V_Elem
      Type(Vect2D)                                 :: GradV_Elem
      PetscInt                                     :: iE, iEloc, iBlkId, iErr
      PetscInt                                     :: NumDoFScal
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      PetscLogDouble                               :: flops       
      
      flops = 0.0
      FractureEnergyBlock = 0.0_Kr
!     DelaminationEnergyBlock = 0.0_Kr

      NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF

      Allocate(V_Loc(NumDoFScal))

      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID

      Do_iEloc: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         Call SectionRealRestrictClosure(V_Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, V_Loc, iErr); CHKERRQ(ierr)
         Do_iGauss: Do iGauss = 1, size(AppCtx%ElemVect(iE)%Gauss_C)
            V_Elem = 0.0_Kr
            GradV_Elem = 0.0_Kr
            Do iDoF2 = 1, NumDoFScal
               V_Elem     = V_Elem     + V_Loc(iDoF2) * AppCtx%ElemScal(iE)%BF(iDoF2, iGauss)
               GradV_Elem = GradV_Elem + V_Loc(iDoF2) * AppCtx%ElemScal(iE)%Grad_BF(iDoF2, iGauss)
            End Do
            
            FractureEnergyBlock = FractureEnergyBlock + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * ( (1.0_Kr-V_Elem)**2 / AppCtx%VarFracSchemeParam%Epsilon + AppCtx%VarFracSchemeParam%Epsilon * (GradV_Elem .DotP. GradV_Elem)  )
         End Do Do_iGauss
      End Do Do_iEloc
      FractureEnergyBlock = FractureEnergyBlock * AppCtx%MatProp(iBlkID)%FracToughness / AppCtx%VarFracSchemeParam%ATCv * 0.25_Kr

      DeAllocate(V_Loc)
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
End Subroutine FractureEnergy_AssemblyBlk_AT2

Subroutine FractureEnergy_AssemblyBlk_AT1(FractureEnergyBlock, iBlk, V_Sec, AppCtx)
      PetscReal, Intent(OUT)                       :: FractureEnergyBlock
!     PetscReal, Intent(OUT)                       :: DelaminationEnergyBlock
      PetscInt                                     :: iBlk
      Type(SectionReal)                            :: V_Sec
      Type(AppCtx_Type)                            :: AppCtx

      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_Elem = \sum_i U_Loc(i) BF(i))
      PetscReal, Dimension(:), Pointer             :: V_Loc
      PetscReal                                    :: V_Elem
      Type(Vect2D)                                 :: GradV_Elem
      PetscInt                                     :: iE, iEloc, iBlkId, iErr
      PetscInt                                     :: NumDoFScal
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      PetscLogDouble                               :: flops       
      
      flops = 0.0
      FractureEnergyBlock = 0.0_Kr
!     DelaminationEnergyBlock = 0.0_Kr

      NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF

      Allocate(V_Loc(NumDoFScal))

      iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID

      Do_iEloc: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         Call SectionRealRestrictClosure(V_Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, V_Loc, iErr); CHKERRQ(ierr)
         Do_iGauss: Do iGauss = 1, size(AppCtx%ElemVect(iE)%Gauss_C)
            V_Elem = 0.0_Kr
            GradV_Elem = 0.0_Kr
            Do iDoF2 = 1, NumDoFScal
               V_Elem     = V_Elem     + V_Loc(iDoF2) * AppCtx%ElemScal(iE)%BF(iDoF2, iGauss)
               GradV_Elem = GradV_Elem + V_Loc(iDoF2) * AppCtx%ElemScal(iE)%Grad_BF(iDoF2, iGauss)
            End Do
            
            FractureEnergyBlock = FractureEnergyBlock + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * (1.0_Kr-V_Elem) / AppCtx%VarFracSchemeParam%Epsilon
            FractureEnergyBlock = FractureEnergyBlock + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * (GradV_Elem .DotP. GradV_Elem) * AppCtx%VarFracSchemeParam%Epsilon 
         End Do Do_iGauss
      End Do Do_iEloc
      FractureEnergyBlock = FractureEnergyBlock * AppCtx%MatProp(iBlkID)%FracToughness / AppCtx%VarFracSchemeParam%ATCv * 0.25_Kr 

      DeAllocate(V_Loc)
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
End Subroutine FractureEnergy_AssemblyBlk_AT1

!----------------------------------------------------------------------------------------!      
! ComputeStrainStress (CM)  
!----------------------------------------------------------------------------------------!      
Subroutine ComputeStrainStress(AppCtx)
     Type(AppCtx_Type)                             :: AppCtx
      
      PetscInt                                     :: iErr
     Type(MatS2D)                                  :: Strain_Elem, Stress_Elem 
     Type(Vect2D)                                  :: F_Elem, U_Elem  
      PetscReal                                    :: Theta_Elem
      PetscReal                                    :: Vol
      PetscInt                                     :: NumDoFVect, NumDofScal
      PetscReal, Dimension(:), Pointer             :: U, Theta
      PetscInt                                     :: iBlk, iBlkID, iELoc, iE, iDoF, iGauss
      PetscReal, Dimension(:), Pointer             :: Stress_Ptr, Strain_Ptr
      PetscLogDouble                               :: flops       
        
      Call PetscLogStagePush (AppCtx%LogInfo%PostProc_Stage, iErr); CHKERRQ(iErr)
      Call PetscLogEventBegin(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)

      flops = 0.0
      Allocate(Stress_Ptr( AppCtx%MeshTopology%Num_Dim * ( AppCtx%MeshTopology%Num_Dim-1 ) / 2))
      Allocate(Strain_Ptr( AppCtx%MeshTopology%Num_Dim * ( AppCtx%MeshTopology%Num_Dim-1 ) / 2))
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         NumDoFVect = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF * AppCtx%MeshTopology%Num_Dim
         NumDoFScal = AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_DoF
         iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
         
         Allocate(U(NumDoFVect))
         Allocate(Theta(NumDoFScal))

         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)

            Call SectionRealRestrictClosure(AppCtx%U%Sec,     AppCtx%MeshTopology%mesh, iE-1, NumDoFVect, U,     iErr); CHKERRQ(ierr)
            Call SectionRealRestrictClosure(AppCtx%Theta%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoFScal, Theta, iErr); CHKERRQ(ierr)
            Strain_Elem       = 0.0_Kr
            Stress_Elem       = 0.0_Kr
            Theta_Elem        = 0.0_Kr
            Vol               = 0.0_Kr
            Do iGauss = 1, Size(AppCtx%ElemVect(iE)%Gauss_C)
               Do iDof = 1, NumDoFScal
                     Theta_Elem  = Theta_Elem + AppCtx%ElemScal(iE)%Gauss_C(iGauss) * AppCtx%ElemScal(iE)%BF(iDoF, iGauss) * Theta(iDoF)
                     Vol         = Vol        + AppCtx%ElemScal(iE)%Gauss_C(iGauss) * AppCtx%ElemScal(iE)%BF(iDoF, iGauss)
                     flops = flops + 5.0
               End Do
               Do iDoF = 1, NumDoFVect
                  Strain_Elem = Strain_Elem + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * AppCtx%ElemVect(iE)%GradS_BF(iDoF, iGauss) * U(iDoF)
               End Do
            End Do
            Strain_Elem = Strain_Elem / Vol
            Theta_Elem  = Theta_Elem / Vol
            Stress_Elem = AppCtx%MatProp(iBlkID)%Hookes_Law * ( Strain_Elem - Theta_Elem * (AppCtx%MatProp(iBlkId)%Therm_Exp) )
            flops = flops + 2.0
            Stress_Ptr = (/ Stress_Elem%XX, Stress_Elem%YY, Stress_Elem%XY /)
            Strain_Ptr = (/ Strain_Elem%XX, Strain_Elem%YY, Strain_Elem%XY /)
            ! Update the Sections with the local values
            Call SectionRealUpdateClosure(AppCtx%StressU, AppCtx%MeshTopology%Mesh, iE-1, Stress_Ptr, INSERT_VALUES, iErr)
            Call SectionRealUpdateClosure(AppCtx%StrainU, AppCtx%MeshTopology%Mesh, iE-1, Strain_Ptr, INSERT_VALUES, iErr)
         End Do Do_Elem_iE
         DeAllocate(U)
         DeAllocate(Theta)
      End Do Do_Elem_iBlk
      DeAllocate(Stress_Ptr)
      DeAllocate(Strain_Ptr)
      Call SectionRealComplete(AppCtx%StressU, iErr); CHKERRQ(iErr)
      Call SectionRealComplete(AppCtx%StrainU, iErr); CHKERRQ(iErr)

      Call PetscLogFlops(flops, iErr); CHKERRQ(iErr)
      Call PetscLogEventEnd  (AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
End Subroutine ComputeStrainStress

End Module m_VarFilmQS_Post