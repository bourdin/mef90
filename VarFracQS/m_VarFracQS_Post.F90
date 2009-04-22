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
   Use m_VarFracQS_U2D
   Use m_VarFracQS_V2D
#elif defined PB_3D
   Use m_VarFracQS_Types3D   
   Use m_VarFracQS_U3D
   Use m_VarFracQS_V3D
#endif   
   Use m_MEF90
   Use m_VarFrac_Struct
   Use petsc
   Use petscmesh

   Implicit NONE   
   
Contains
!----------------------------------------------------------------------------------------!      
! ComputeEnergy (CM)  
!----------------------------------------------------------------------------------------!      
   Subroutine ComputeEnergy(AppCtx)
      !!! CM
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoFVect, NumDoFScal, NumGauss
      PetscReal, Dimension(:), Pointer             :: F, U, V, Theta
      PetscInt                                     :: iBlk, iELoc, iE
      PetscInt                                     :: iDoF, iGauss
      PetscReal                                    :: Theta_Elem, V_Elem
#if defined PB_2D
      Type(MatS2D)                                 :: Strain_Elem, Stress_Elem 
      Type(MatS2D)                                 :: Effective_Strain_Elem
      Type(Vect2D)                                 :: F_Elem, U_Elem, GradV_Elem 

#elif defined PB_3D
      Type(MatS3D)                                 :: Strain_Elem, Stress_Elem 
      Type(MatS3D)                                 :: Effective_Strain_Elem
      Type(Vect3D)                                 :: F_Elem, U_Elem,  GradV_Elem  
#endif
      PetscReal                                    :: MyExtForcesWork, MyElasticEnergy, MySurfaceEnergy
     

      Call PetscLogStagePush(AppCtx%LogInfo%PostProc_Stage, iErr); CHKERRQ(iErr)
      Call PetscLogEventBegin(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      
      MyElasticEnergy = 0.0_Kr
      MyExtForcesWork = 0.0_Kr
      MySUrfaceEnergy = 0.0_Kr
      !---------------------------------------------------------------------
      ! Surface Energy      : Gc (1/ (4 eps) (1-v)^2 + eps GradV * GradV)
      ! Elastic energy      : 1/2 v^2 (A*EffectiveStrain) * EffectiveStrain 
      ! Work of body forces : - F*U   
      !---------------------------------------------------------------------
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            ! Define the indices
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            NumDoFVect = Size(AppCtx%ElemVect(iE)%BF,1)
            NumDoFScal = Size(AppCtx%ElemScal(iE)%BF,1)
            NumGauss   = Size(AppCtx%ElemVect(iE)%BF,2)
            ! Allocate the variables
            Allocate(F(NumDoFVect))
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%F, iE-1, NumDoFVect, F, iErr); CHKERRQ(ierr)
            Allocate(U(NumDoFVect))
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%U, iE-1, NumDoFVect, U, iErr); CHKERRQ(ierr)
            Allocate(Theta(NumDoFScal))
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%Theta, iE-1, NumDoFScal, Theta, iErr); CHKERRQ(ierr)
            Allocate(V(NumDoFScal))
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%V, iE-1, NumDoFScal, V, iErr); CHKERRQ(ierr)
            Do iGauss = 1, NumGauss
            ! Inizialize the variables at the gauss points
               Strain_Elem            = 0.0_Kr
               Stress_Elem            = 0.0_Kr
               Theta_Elem             = 0.0_Kr
               V_Elem                 = 0.0_Kr
               GradV_Elem             = 0.0_Kr 
               Effective_Strain_Elem  = 0.0_Kr
               F_Elem                 = 0.0_Kr
               U_Elem                 = 0.0_Kr
            ! Calculate the variables at the gauss points
               Do iDof = 1, NumDoFScal
                  Theta_Elem          = Theta_Elem + AppCtx%ElemScal(iE)%BF(iDoF, iGauss) * Theta(iDoF)
                  V_Elem              = V_Elem     + AppCtx%ElemScal(iE)%BF(iDoF, iGauss) * V(iDoF)
                  GradV_Elem          = V_Elem     + AppCtx%ElemScal(iE)%BF(iDoF, iGauss) * V(iDoF)
               End Do
               Do iDoF = 1, NumDoFVect
                  Strain_Elem         = Strain_Elem + AppCtx%ElemVect(iE)%GradS_BF(iDoF, iGauss) * U(iDoF)
                  F_Elem              = F_Elem      + AppCtx%ElemVect(iE)%BF(iDoF, iGauss)       * F(iDoF) 
                  U_Elem              = U_Elem      + AppCtx%ElemVect(iE)%BF(iDoF, iGauss)       * U(iDoF)
               End Do            
               Effective_Strain_Elem  = Strain_Elem - Theta_Elem * AppCtx%MatProp(iBlk)%Therm_Exp   ;
               Stress_Elem            = AppCtx%MatProp( AppCtx%MeshTopology%Elem_Blk(iBlk)%ID )%Hookes_Law * Effective_Strain_Elem
            ! Calculate the elastic energy
               MyElasticEnergy  = MyElasticEnergy  + AppCtx%ElemVect(iE)%Gauss_C(iGauss) *  V_Elem**2 * (Stress_Elem .DotP. Effective_Strain_Elem) * 0.5_Kr
            ! Calculate the work of body forces
               MyExtForcesWork = MyExtForcesWork + AppCtx%ElemVect(iE)%Gauss_C(iGauss) *  (F_Elem .DotP. U_Elem)
            ! Calculate the suface energy
               MySurfaceEnergy  = MySurfaceEnergy  + AppCtx%MatProp(iBlk)%Toughness * AppCtx%ElemVect(iE)%Gauss_C(iGauss) *  0.25_Kr / AppCtx%VarFracSchemeParam%Epsilon *  ( 1.0_Kr - V_Elem)**2 +  AppCtx%VarFracSchemeParam%Epsilon * (GradV_Elem .DotP. GradV_Elem)
            Call PetscLogFlops(AppCtx%MeshTopology%Num_Dim+4, iErr)
            End Do
            ! DeAllocate the variables
            DeAllocate(F)
            DeAllocate(U)
            DeAllocate(Theta)
            DeAllocate(V)
         End Do Do_Elem_iE
      End Do Do_Elem_iBlk
      ! Global sum among  
      Call PetscGlobalSum(MyElasticEnergy, AppCtx%ElasticEnergy, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      Call PetscGlobalSum(MyExtForcesWork, AppCtx%ExtForcesWork, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)           
      Call PetscGlobalSum(MySurfaceEnergy, AppCtx%SurfaceEnergy, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)    
      AppCtx%TotalEnergy = AppCtx%ElasticEnergy - AppCtx%ExtForcesWork + AppCtx%SurfaceEnergy
      Call PetscLogEventEnd  (AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine ComputeEnergy
   
!----------------------------------------------------------------------------------------!      
! ComputeStrainStress (CM)  
!----------------------------------------------------------------------------------------!      
   Subroutine ComputeStrainStress(AppCtx)
     Type(AppCtx_Type)                             :: AppCtx
      
      PetscInt                                     :: iErr
#if defined PB_2D
     Type(MatS2D)                                 :: Strain_Elem, Stress_Elem 
     Type(Vect2D)                                 :: F_Elem, U_Elem  
#elif defined PB_3D 
     Type(MatS3D)                                 :: Strain_Elem, Stress_Elem 
     Type(Vect3D)                          :: F_Elem, U_Elem  
#endif
      PetscReal                                    :: Theta_Elem, V_Elem
      PetscReal                                    :: Vol
      PetscInt                                     :: NumDoFVect, NumDofScal, NumGauss
      PetscReal, Dimension(:), Pointer             :: U, Theta, V
      PetscInt                                     :: iBlk, iELoc, iE
      PetscInt                                     :: iDoF, iGauss
      PetscReal, Dimension(:), Pointer             :: Stress_Ptr, Strain_Ptr
       
        
!      Call PetscLogEventBegin(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
!      Call PetscLogStagePush (AppCtx%LogInfo%PostProc_Stage, iErr); CHKERRQ(iErr)
      Allocate(Stress_Ptr( AppCtx%MeshTopology%Num_Dim * ( AppCtx%MeshTopology%Num_Dim-1 ) / 2))
      Allocate(Strain_Ptr( AppCtx%MeshTopology%Num_Dim * ( AppCtx%MeshTopology%Num_Dim-1 ) / 2))
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            NumDoFScal   = Size(AppCtx%ElemScal(iE)%BF,1)
            NumDoFVect   = Size(AppCtx%ElemVect(iE)%BF,1)
            NumGauss = Size(AppCtx%ElemVect(iE)%BF,2)
            Allocate(U(NumDoFVect))
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%U, iE-1, NumDoFVect, U, iErr); CHKERRQ(ierr)
            Allocate(Theta(NumDoFScal))
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%Theta, iE-1, NumDoFScal, Theta, iErr); CHKERRQ(ierr)
            Allocate(V(NumDoFScal))
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%V, iE-1, NumDoFScal, V, iErr); CHKERRQ(ierr)
            V_Elem                = 0.0_Kr
            Strain_Elem           = 0.0_Kr
            Stress_Elem           = 0.0_Kr
            Theta_Elem            = 0.0_Kr
            Vol  = 0.0_Kr
            Do iDof = 1, NumDoFScal
                  V_Elem          = V_Elem     + AppCtx%ElemScal(iE)%BF(iDoF, iGauss) * V(iDoF)
                  Theta_Elem      = Theta_Elem + AppCtx%ElemScal(iE)%BF(iDoF, iGauss) * Theta(iDoF)
                  Vol             = Vol + AppCtx%ElemScal(iE)%Gauss_C(iGauss) * AppCtx%ElemScal(iE)%BF(iDoF, iGauss)
            End Do
            Do iGauss = 1, NumGauss
               Do iDoF = 1, NumDoFVect
                  Strain_Elem            = Strain_Elem + AppCtx%ElemVect(iE)%GradS_BF(iDoF, iGauss) * U(iDoF)
                  Stress_Elem            = Stress_Elem + (V_Elem**2) * AppCtx%MatProp( AppCtx%MeshTopology%Elem_Blk(iBlk)%ID )%Hookes_Law * ( Strain_Elem - Theta_Elem * (AppCtx%MatProp(iBlk)%Therm_Exp) )
!             Call PetscLogFlops(3*AppCtx%MeshTopology%Num_Dim+2, iErr)
               End Do
            End Do
            Strain_Elem = Strain_Elem / Vol
            Stress_Elem = Stress_Elem / Vol
            Call PetscLogFlops(AppCtx%MeshTopology%Num_Dim, iErr)
#if defined PB_2D
            Stress_Ptr = (/ Stress_Elem%XX, Stress_Elem%YY, Stress_Elem%XY /)
            Strain_Ptr = (/ Strain_Elem%XX, Strain_Elem%YY, Strain_Elem%XY /)
#elif defined PB_3D
            Stress_Ptr = (/ Stress_Elem%XX, Stress_Elem%YY, Stress_Elem%ZZ, Stress_Elem%YZ, Stress_Elem%XZ, Stress_Elem%XY /)
            Strain_Ptr = (/ Strain_Elem%XX, Strain_Elem%YY, Strain_Elem%ZZ, Strain_Elem%YZ, Strain_Elem%XZ, Strain_Elem%XY  /)
#endif
            ! Update the Sections with the local values
            Call MeshUpdateClosure(AppCtx%MeshTopology%Mesh, AppCtx%StressU, iE-1, Stress_Ptr, iErr)
            Call MeshUpdateClosure(AppCtx%MeshTopology%Mesh, AppCtx%StrainU, iE-1, Strain_Ptr, iErr)
            DeAllocate(U)
            DeAllocate(Theta)
         End Do Do_Elem_iE
      End Do Do_Elem_iBlk
      DeAllocate(Stress_Ptr)
      DeAllocate(Strain_Ptr)
      Call SectionRealComplete(AppCtx%StressU, iErr); CHKERRQ(iErr)
      Call SectionRealComplete(AppCtx%StrainU, iErr); CHKERRQ(iErr)

!      Call PetscLogEventEnd  (AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
!      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine ComputeStrainStress
!----------------------------------------------------------------------------------------!

#if defined PB_2D
End Module m_VarFracQS_Post2D
#elif defined PB_3D
End Module m_VarFracQS_Post3D
#endif
