#if defined PB_2D
Module m_VarFracQS_Post2D
#elif defined PB_3D
Module m_VarFracQS_Post3D
#endif

#include "finclude/petscdef.h"
#include "finclude/petscmeshdef.h"
#include "finclude/petscviewerdef.h"

#if defined PB_2D
   Use m_VarFracQS_Types2D
#elif defined PB_3D
   Use m_VarFracQS_Types3D   
#endif   
   Use m_MEF90
   Use m_VarFrac_Struct
   Use petsc
   Use petscmesh

   Implicit NONE   
   
Contains
   Subroutine ElasticEnergy_Assembly(ElasticEnergy, ElasticEnergyBlock, AppCtx)     
      PetscReal, Intent(OUT)                       :: ElasticEnergy
      PetscReal, Dimension(:), Pointer             :: ElasticEnergyBlock
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iBlk, iBlkId, iErr
      PetscReal                                    :: MyElasticEnergy
      PetscReal, Dimension(:), Pointer             :: MyElasticEnergyBlock
     

      Call PetscLogStagePush(AppCtx%LogInfo%PostProc_Stage, iErr); CHKERRQ(iErr)
      Call PetscLogEventBegin(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      
      MyElasticEnergy = 0.0_Kr
      Allocate(MyELasticEnergyBlock(AppCtx%MeshTopology%Num_Elem_Blks_Global))
      MyElasticEnergyBlock = 0.0_Kr
      Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
         If (AppCtx%MyEXO%EBProperty(VarFrac_EBProp_IsBrittle)%Value(iBlkID) /= 0) Then
            Select Case (AppCtx%VarFracSchemeParam%Unilateral)
            Case (VarFrac_Unilateral_NONE)
               Call ElasticEnergy_AssemblyBlk_Brittle(MyElasticEnergyBlock(iBlkID), iBlk, AppCtx%U%Sec, AppCtx%Theta%Sec, AppCtx%V%Sec, AppCtx)
            Case (VarFrac_Unilateral_Full)
               Call ElasticEnergy_AssemblyBlk_BrittleUnilateralFull(MyElasticEnergyBlock(iBlkId), iBlk, AppCtx%U%Sec, AppCtx%Theta%Sec, AppCtx%V%Sec, AppCtx)
            Case (VarFrac_Unilateral_Shear)
               Call ElasticEnergy_AssemblyBlk_BrittleUnilateralShear(MyElasticEnergyBlock(iBlkId), iBlk, AppCtx%U%Sec, AppCtx%Theta%Sec, AppCtx%V%Sec, AppCtx)
            End Select
         Else
            Call ElasticEnergy_AssemblyBlk_NonBrittle(MyElasticEnergyBlock(iBlkId), iBlk, AppCtx%U%Sec, AppCtx%Theta%Sec, AppCtx)
         End If
         MyElasticEnergy = MyElasticEnergy + MyElasticEnergyBlock(iBlkId)
      End Do

      Call MPI_AllReduce(MyElasticEnergy, ElasticEnergy, 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      Call MPI_AllReduce(MyElasticEnergyBlock, ElasticEnergyBlock, AppCtx%MeshTopology%Num_Elem_Blks_Global, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      !Do iBlkID = 1, AppCtx%MeshTopology%Num_Elem_Blks_Global
      !   Call PetscGlobalSum(MyElasticEnergyBlock(iBlkID), ElasticEnergyBlock(iBlkID), PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      !End Do
      DeAllocate(MyElasticEnergyBlock)
      Call PetscLogEventEnd(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine ElasticEnergy_Assembly
   
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
         If (AppCtx%MyEXO%EBProperty(VarFrac_EBProp_HasBForce)%Value(iBlkID) /= 0) Then
            Call ExtForcesWork_AssemblyBlk(MyExtForcesWorkBlock(iBlkID), iBlk, AppCtx%U%Sec, AppCtx%F%Sec, AppCtx)
            MyExtForcesWork = MyExtForcesWork + MyExtForcesWorkBlock(iBlkID)
         End If
      End Do

      !Call PetscGlobalSum(MyExtForcesWork, ExtForcesWork, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      !Do iBlkID = 1, AppCtx%MeshTopology%Num_Elem_Blks_Global
      !   Call PetscGlobalSum(MyExtForcesWorkBlock(iBlkID), ExtForcesWorkBlock(iBlkID), PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      !End Do
      Call MPI_AllReduce(MyExtForcesWork, ExtForcesWork, 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      Call MPI_AllReduce(MyExtForcesWorkBlock, ExtForcesWorkBlock, AppCtx%MeshTopology%Num_Elem_Blks_Global, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      DeAllocate(MyExtForcesWorkBlock)
      Call PetscLogEventEnd(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine ExtForcesWork_Assembly
   
   Subroutine SurfaceEnergy_Assembly(SurfaceEnergy, SurfaceEnergyBlock, AppCtx)     
      PetscReal, Intent(OUT)                       :: SurfaceEnergy
      PetscReal, Dimension(:), Pointer             :: SurfaceEnergyBlock
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iBlk, iBlkId, iErr
      PetscReal                                    :: MySurfaceEnergy
      PetscReal, Dimension(:), Pointer             :: MySurfaceEnergyBlock
     

      Call PetscLogStagePush(AppCtx%LogInfo%PostProc_Stage, iErr); CHKERRQ(iErr)
      Call PetscLogEventBegin(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      
      MySurfaceEnergy = 0.0_Kr
      Allocate(MySurfaceEnergyBlock(AppCtx%MeshTopology%Num_Elem_Blks_Global))
      Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         iBlkID = AppCtx%MeshTopology%Elem_Blk(iBlk)%ID
         Select Case (AppCtx%VarFracSchemeParam%AtNum)
         Case(1)
            Call SurfaceEnergy_AssemblyBlk_AT1(MySurfaceEnergyBlock(iBlkID), iBlk, AppCtx%V%Sec, AppCtx)
         Case(2)
            Call SurfaceEnergy_AssemblyBlk_AT2(MySurfaceEnergyBlock(iBlkID), iBlk, AppCtx%V%Sec, AppCtx)
         Case Default
          SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP, 'Only AT1 and AT2 are implemented\n', iErr)
         End Select
         MySurfaceEnergy = MySurfaceEnergy + MySurfaceEnergyBlock(iBlkID)
      End Do

!      Call PetscGlobalSum(MySurfaceEnergy, SurfaceEnergy, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
!      Do iBlkID = 1, AppCtx%MeshTopology%Num_Elem_Blks_Global
!         Call PetscGlobalSum(MySurfaceEnergyBlock(iBlkID), SurfaceEnergyBlock(iBlkID), PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
!      End Do
      Call MPI_AllReduce(MySurfaceEnergy, SurfaceEnergy, 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      Call MPI_AllReduce(MySurfaceEnergyBlock, SurfaceEnergyBlock, AppCtx%MeshTopology%Num_Elem_Blks_Global, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      DeAllocate(MySurfaceEnergyBlock)
      Call PetscLogEventEnd(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine SurfaceEnergy_Assembly
   
!!!
!!! Block Assembly Routines
!!!
   Subroutine ElasticEnergy_AssemblyBlk_Brittle(ElasticEnergyBlock, iBlk, U_Sec, Theta_Sec, V_Sec, AppCtx)
      PetscReal, Intent(OUT)                       :: ElasticEnergyBlock
      PetscInt                                     :: iBlk
      Type(SectionReal)                            :: U_Sec, Theta_Sec, V_Sec
      Type(AppCtx_Type)                            :: AppCtx

      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_Elem = \sum_i U_Loc(i) BF(i))
      PetscReal, Dimension(:), Pointer             :: U_Loc, V_Loc, Theta_Loc
#if defined PB_2D
      Type(Vect2D)                                 :: U_Elem
      Type(Mats2D)                                 :: Strain_Elem, EffectiveStrain_Elem
#elif defined PB_3D  
      Type(Vect3D)                                 :: U_Elem    
      Type(Mats3D)                                 :: Strain_Elem, EffectiveStrain_Elem
#endif
      PetscReal                                    :: V_Elem, Theta_Elem
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
            ElasticEnergyBlock = ElasticEnergyblock + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * (V_Elem**2 + AppCtx%VarFracSchemeParam%KEpsilon) * ((AppCtx%MatProp(iBlkId)%Hookes_Law * EffectiveStrain_Elem) .DotP. EffectiveStrain_Elem ) * 0.5_Kr
            flops = flops + 6.0 
         End Do Do_iGauss
      End Do Do_iEloc

      DeAllocate(U_Loc)
      DeAllocate(Theta_Loc)
      DeAllocate(V_Loc)
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
   End Subroutine ElasticEnergy_AssemblyBlk_Brittle

   Subroutine ElasticEnergy_AssemblyBlk_BrittleUnilateralFull(ElasticEnergyBlock, iBlk, U_Sec, Theta_Sec, V_Sec, AppCtx)
      PetscReal, Intent(OUT)                       :: ElasticEnergyBlock
      PetscInt                                     :: iBlk
      Type(SectionReal)                            :: U_Sec, Theta_Sec, V_Sec
      Type(AppCtx_Type)                            :: AppCtx

      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_Elem = \sum_i U_Loc(i) BF(i))
      PetscReal, Dimension(:), Pointer             :: U_Loc, V_Loc, Theta_Loc
#if defined PB_2D
      Type(Vect2D)                                 :: U_Elem
      Type(Mats2D)                                 :: Strain_Elem, EffectiveStrain_Elem
      Type(MatS2D)                                 :: EffectiveStrain_Elem_S, EffectiveStrain_Elem_D
#elif defined PB_3D  
      Type(Vect3D)                                 :: U_Elem    
      Type(Mats3D)                                 :: Strain_Elem, EffectiveStrain_Elem
      Type(MatS3D)                                 :: EffectiveStrain_Elem_S, EffectiveStrain_Elem_D
#endif
      PetscReal                                    :: Strain_Trace
      PetscReal                                    :: V_Elem, Theta_Elem
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
            Strain_Trace = Trace(Strain_Elem)
            If (Strain_Trace >= 0.0_Kr) Then
               ElasticEnergyBlock = ElasticEnergyblock + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * (V_Elem**2 + AppCtx%VarFracSchemeParam%KEpsilon) * ((AppCtx%MatProp(iBlkId)%Hookes_Law * EffectiveStrain_Elem) .DotP. EffectiveStrain_Elem ) * 0.5_Kr
               flops = flops + 6.0 
            Else
               EffectiveStrain_Elem_S = SphericalPart(EffectiveStrain_Elem)
               EffectiveStrain_Elem_D = DeviatoricPart(EffectiveStrain_Elem)

               ElasticEnergyBlock = ElasticEnergyblock + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * ((V_Elem**2 + AppCtx%VarFracSchemeParam%KEpsilon) * ((AppCtx%MatProp(iBlkId)%Hookes_Law * EffectiveStrain_Elem_D) .DotP. EffectiveStrain_Elem_D )  + ((AppCtx%MatProp(iBlkId)%Hookes_Law * EffectiveStrain_Elem_S) .DotP. EffectiveStrain_Elem_S ) ) * 0.5_Kr
               flops = flops + 6.0 
            End If
         End Do Do_iGauss
      End Do Do_iEloc

      DeAllocate(U_Loc)
      DeAllocate(Theta_Loc)
      DeAllocate(V_Loc)
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
   End Subroutine ElasticEnergy_AssemblyBlk_BrittleUnilateralFull

   Subroutine ElasticEnergy_AssemblyBlk_BrittleUnilateralShear(ElasticEnergyBlock, iBlk, U_Sec, Theta_Sec, V_Sec, AppCtx)
      PetscReal, Intent(OUT)                       :: ElasticEnergyBlock
      PetscInt                                     :: iBlk
      Type(SectionReal)                            :: U_Sec, Theta_Sec, V_Sec
      Type(AppCtx_Type)                            :: AppCtx

      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_Elem = \sum_i U_Loc(i) BF(i))
      PetscReal, Dimension(:), Pointer             :: U_Loc, V_Loc, Theta_Loc
#if defined PB_2D
      Type(Vect2D)                                 :: U_Elem
      Type(Mats2D)                                 :: Strain_Elem, EffectiveStrain_Elem
      Type(MatS2D)                                 :: EffectiveStrain_Elem_S, EffectiveStrain_Elem_D
#elif defined PB_3D  
      Type(Vect3D)                                 :: U_Elem    
      Type(Mats3D)                                 :: Strain_Elem, EffectiveStrain_Elem
      Type(MatS3D)                                 :: EffectiveStrain_Elem_S, EffectiveStrain_Elem_D
#endif
      PetscReal                                    :: Strain_Trace
      PetscReal                                    :: V_Elem, Theta_Elem
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
            EffectiveStrain_Elem   = Strain_Elem - AppCtx%MatProp(iBlkId)%Therm_Exp * Theta_Elem            
            EffectiveStrain_Elem_S = SphericalPart(EffectiveStrain_Elem)
            EffectiveStrain_Elem_D = DeviatoricPart(EffectiveStrain_Elem)

            ElasticEnergyBlock = ElasticEnergyblock + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * ((V_Elem**2 + AppCtx%VarFracSchemeParam%KEpsilon) * ((AppCtx%MatProp(iBlkId)%Hookes_Law * EffectiveStrain_Elem_D) .DotP. EffectiveStrain_Elem_D )  + ((AppCtx%MatProp(iBlkId)%Hookes_Law * EffectiveStrain_Elem_S) .DotP. EffectiveStrain_Elem_S ) ) * 0.5_Kr
            flops = flops + 6.0 
         End Do Do_iGauss
      End Do Do_iEloc

      DeAllocate(U_Loc)
      DeAllocate(Theta_Loc)
      DeAllocate(V_Loc)
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
   End Subroutine ElasticEnergy_AssemblyBlk_BrittleUnilateralShear

   Subroutine ElasticEnergy_AssemblyBlk_NonBrittle(ElasticEnergyBlock, iBlk, U_Sec, Theta_Sec, AppCtx)
      PetscReal, Intent(OUT)                       :: ElasticEnergyBlock
      PetscInt                                     :: iBlk
      Type(SectionReal)                            :: U_Sec, Theta_Sec
      Type(AppCtx_Type)                            :: AppCtx

      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_ELem = \sum_i U_Loc(i) BF(i))
      PetscReal, Dimension(:), Pointer             :: U_Loc, Theta_Loc
#if defined PB_2D
      Type(Mats2D)                                 :: Strain_Elem, EffectiveStrain_Elem
#elif defined PB_3D  
      Type(Mats3D)                                 :: Strain_Elem, EffectiveStrain_Elem
#endif
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
   End Subroutine ElasticEnergy_AssemblyBlk_NonBrittle

   Subroutine ExtForcesWork_AssemblyBlk(ExtForcesWorkBlock, iBlk, U_Sec, F_Sec, AppCtx)
      PetscReal, Intent(OUT)                       :: ExtForcesWorkBlock
      PetscInt                                     :: iBlk
      Type(SectionReal)                            :: U_Sec, F_Sec
      Type(AppCtx_Type)                            :: AppCtx

      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_Elem = \sum_i U_Loc(i) BF(i))
      PetscReal, Dimension(:), Pointer             :: U_Loc, F_Loc
#if defined PB_2D
      Type(Vect2D)                                 :: U_Elem, F_Elem
#elif defined PB_3D  
      Type(Vect3D)                                 :: U_Elem, F_Elem
#endif
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

   Subroutine SurfaceEnergy_AssemblyBlk_AT2(SurfaceEnergyBLock, iBlk, V_Sec, AppCtx)
      PetscReal, Intent(OUT)                       :: SurfaceEnergyBlock
      PetscInt                                     :: iBlk
      Type(SectionReal)                            :: V_Sec
      Type(AppCtx_Type)                            :: AppCtx

      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_Elem = \sum_i U_Loc(i) BF(i))
      PetscReal, Dimension(:), Pointer             :: V_Loc
      PetscReal                                    :: V_Elem
#if defined PB_2D
      Type(Vect2D)                                 :: GradV_Elem
#elif defined PB_3D  
      Type(Vect3D)                                 :: GradV_Elem
#endif
      PetscInt                                     :: iE, iEloc, iBlkId, iErr
      PetscInt                                     :: NumDoFScal
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      PetscLogDouble                               :: flops       
      
      flops = 0.0
      SurfaceEnergyBlock = 0.0_Kr

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
            
            SurfaceEnergyBlock = SurfaceEnergyBlock + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * ( (1.0_Kr-V_Elem)**2 / AppCtx%VarFracSchemeParam%Epsilon + AppCtx%VarFracSchemeParam%Epsilon * (GradV_Elem .DotP. GradV_Elem)  )
         End Do Do_iGauss
      End Do Do_iEloc
      SurfaceEnergyBlock = SurfaceEnergyBlock * AppCtx%MatProp(iBlkID)%Toughness / AppCtx%VarFracSchemeParam%ATCv * 0.25_Kr

      DeAllocate(V_Loc)
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
   End Subroutine SurfaceEnergy_AssemblyBlk_AT2

   Subroutine SurfaceEnergy_AssemblyBlk_AT1(SurfaceEnergyBlock, iBlk, V_Sec, AppCtx)
      PetscReal, Intent(OUT)                       :: SurfaceEnergyBlock
      PetscInt                                     :: iBlk
      Type(SectionReal)                            :: V_Sec
      Type(AppCtx_Type)                            :: AppCtx

      !!!   _Loc are restriction of fields to local patch (the element)
      !!!   _Elem are local contribution over the element (u_Elem = \sum_i U_Loc(i) BF(i))
      PetscReal, Dimension(:), Pointer             :: V_Loc
      PetscReal                                    :: V_Elem
#if defined PB_2D
      Type(Vect2D)                                 :: GradV_Elem
#elif defined PB_3D  
      Type(Vect3D)                                 :: GradV_Elem
#endif
      PetscInt                                     :: iE, iEloc, iBlkId, iErr
      PetscInt                                     :: NumDoFScal
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      PetscLogDouble                               :: flops       
      
      flops = 0.0
      SurfaceEnergyBlock = 0.0_Kr

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
            
            SurfaceEnergyBlock = SurfaceEnergyBlock + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * (1.0_Kr-V_Elem) / AppCtx%VarFracSchemeParam%Epsilon
            SurfaceEnergyBlock = SurfaceEnergyBlock + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * (GradV_Elem .DotP. GradV_Elem) * AppCtx%VarFracSchemeParam%Epsilon 
         End Do Do_iGauss
      End Do Do_iEloc
      SurfaceEnergyBlock = SurfaceEnergyBlock * AppCtx%MatProp(iBlkID)%Toughness / AppCtx%VarFracSchemeParam%ATCv * 0.25_Kr 

      DeAllocate(V_Loc)
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
   End Subroutine SurfaceEnergy_AssemblyBlk_AT1

!----------------------------------------------------------------------------------------!      
! ComputeStrainStress (CM)  
!----------------------------------------------------------------------------------------!      
   Subroutine ComputeStrainStress(AppCtx)
     Type(AppCtx_Type)                             :: AppCtx
      
      PetscInt                                     :: iErr
#if defined PB_2D
     Type(MatS2D)                                  :: Strain_Elem, Stress_Elem 
     Type(Vect2D)                                  :: F_Elem, U_Elem  
#elif defined PB_3D 
     Type(MatS3D)                                  :: Strain_Elem, Stress_Elem 
     Type(Vect3D)                                  :: F_Elem, U_Elem  
#endif
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
#if defined PB_2D
            Stress_Ptr = (/ Stress_Elem%XX, Stress_Elem%YY, Stress_Elem%XY /)
            Strain_Ptr = (/ Strain_Elem%XX, Strain_Elem%YY, Strain_Elem%XY /)
#elif defined PB_3D
            Stress_Ptr = (/ Stress_Elem%XX, Stress_Elem%YY, Stress_Elem%ZZ, Stress_Elem%YZ, Stress_Elem%XZ, Stress_Elem%XY /)
            Strain_Ptr = (/ Strain_Elem%XX, Strain_Elem%YY, Strain_Elem%ZZ, Strain_Elem%YZ, Strain_Elem%XZ, Strain_Elem%XY  /)
#endif
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

#if defined PB_2D
End Module m_VarFracQS_Post2D
#elif defined PB_3D
End Module m_VarFracQS_Post3D
#endif
