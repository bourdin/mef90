Module m_VarFracFilm_Post

#include "finclude/petscdef.h"
#include "finclude/petscmeshdef.h"
#include "finclude/petscviewerdef.h"

   Use m_VarFracFilm_Types
   Use m_VarFracFilm_U
   Use m_VarFracFilm_V
   Use m_VarFracFilm_Phi
   Use m_MEF90
   Use m_VarFracFilm_Struct
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
      PetscReal, Dimension(:), Pointer             :: U, U0, V, Theta, Phi
      PetscInt                                     :: iBlk, iELoc, iE
      PetscInt                                     :: iDoF, iGauss
      PetscReal                                    :: Theta_Elem, V_Elem
      Type(MatS2D)                                 :: Strain_Elem, Stress_Elem 
      Type(MatS2D)                                 :: Effective_Strain_Elem
      Type(Vect2D)                                 :: DU_Elem, GradV_Elem 
      PetscReal                                    :: MyElasticBulkEnergy, MyElasticInterEnergy, MySurfaceEnergyT, MySurfaceEnergyD
     

      Call PetscLogStagePush(AppCtx%LogInfo%PostProc_Stage, iErr); CHKERRQ(iErr)
      Call PetscLogEventBegin(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      
      MyElasticBulkEnergy  = 0.0_Kr
      MyElasticInterEnergy = 0.0_Kr
      MySurfaceEnergyT     = 0.0_Kr
      MySurfaceEnergyD     = 0.0_Kr
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
            Allocate(U0(NumDoFVect))
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%U0, iE-1, NumDoFVect, U0, iErr); CHKERRQ(ierr)
            Allocate(U(NumDoFVect))
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%U, iE-1, NumDoFVect, U, iErr); CHKERRQ(ierr)
            Allocate(Theta(NumDoFScal))
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%Theta, iE-1, NumDoFScal, Theta, iErr); CHKERRQ(ierr)
            Allocate(V(NumDoFScal))
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%V, iE-1, NumDoFScal, V, iErr); CHKERRQ(ierr)
            Allocate(Phi(1))
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%Phi, iE-1, NumDoFScal, Phi, iErr); CHKERRQ(ierr)
            Do iGauss = 1, NumGauss
            ! Inizialize the variables at the gauss points
               Strain_Elem            = 0.0_Kr
               Stress_Elem            = 0.0_Kr
               Theta_Elem             = 0.0_Kr
               V_Elem                 = 0.0_Kr
               GradV_Elem             = 0.0_Kr 
               Effective_Strain_Elem  = 0.0_Kr
               DU_Elem                = 0.0_Kr
            ! Calculate the variables at the gauss points
               Do iDof = 1, NumDoFScal
                  Theta_Elem          = Theta_Elem + AppCtx%ElemScal(iE)%BF(iDoF, iGauss) * Theta(iDoF)
                  V_Elem              = V_Elem     + AppCtx%ElemScal(iE)%BF(iDoF, iGauss) * V(iDoF)
                  GradV_Elem          = GradV_Elem + AppCtx%ElemScal(iE)%Grad_BF(iDoF, iGauss) * V(iDoF)
               End Do
               Do iDoF = 1, NumDoFVect
                  Strain_Elem         = Strain_Elem + AppCtx%ElemVect(iE)%GradS_BF(iDoF, iGauss) * U(iDoF)
                  DU_Elem             = DU_Elem     + AppCtx%ElemVect(iE)%BF(iDoF, iGauss)       * (U(iDoF) - U0(iDof))
               End Do            
               Effective_Strain_Elem  = Strain_Elem - Theta_Elem * AppCtx%MatProp(iBlk)%Therm_Exp   ;
               Stress_Elem            = AppCtx%MatProp( AppCtx%MeshTopology%Elem_Blk(iBlk)%ID )%Hookes_Law * Effective_Strain_Elem
            ! Should we include the energy due to residual stiffness ?????? Add a further contribution ?
            ! Calculate the bulk elastic energy
               MyElasticBulkEnergy  = MyElasticBulkEnergy  + AppCtx%ElemVect(iE)%Gauss_C(iGauss) *  V_Elem**2 * (Stress_Elem .DotP. Effective_Strain_Elem) * 0.5_Kr
            ! Calculate the interface elastic energy
               MyElasticInterEnergy = MyElasticInterEnergy +  0.5_Kr * AppCtx%ElemVect(iE)%Gauss_C(iGauss) *  (( AppCtx%MatProp( AppCtx%MeshTopology%Elem_Blk(iBlk)%ID )%K_interface * DU_Elem ) .DotP. DU_Elem)
            ! Calculate the suface energy of transverse crack
               MySurfaceEnergyT  = MySurfaceEnergyT  + Phi(1) * AppCtx%MatProp(iBlk)%ToughnessT * AppCtx%ElemVect(iE)%Gauss_C(iGauss) *  ((0.25_Kr / AppCtx%VarFracFilmSchemeParam%Epsilon) *  ( 1.0_Kr - V_Elem)**2 + ( AppCtx%VarFracFilmSchemeParam%Epsilon * (GradV_Elem .DotP. GradV_Elem)) )
            ! Calculate the suface energy due to delamination
            ! see the comment on the volume for calculate the stresses
               MySurfaceEnergyD  = MySurfaceEnergyD + AppCtx%MatProp(iBlk)%ToughnessD * AppCtx%ElemVect(iE)%Gauss_C(iGauss) * (1-Phi(1))
            ! For Blaise: what's the line below????                 
            Call PetscLogFlops(AppCtx%MeshTopology%Num_Dim+4, iErr)
            
            End Do
               
            ! DeAllocate 
            DeAllocate(U0)
            DeAllocate(U)
            DeAllocate(Theta)
            DeAllocate(V)
            DeAllocate(Phi)
         End Do Do_Elem_iE
      End Do Do_Elem_iBlk
      
      ! Global sum of the energies
      Call PetscGlobalSum(MyElasticBulkEnergy, AppCtx%ElasticBulkEnergy, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      Call PetscGlobalSum(MyElasticInterEnergy, AppCtx%ElasticInterEnergy, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)           
      Call PetscGlobalSum(MySurfaceEnergyT, AppCtx%SurfaceEnergyT, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)    
      Call PetscGlobalSum(MySurfaceEnergyD, AppCtx%SurfaceEnergyD, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)    
      AppCtx%TotalEnergy = AppCtx%ElasticBulkEnergy+ AppCtx%ElasticInterEnergy + AppCtx%SurfaceEnergyT + AppCtx%SurfaceEnergyD
      Call PetscLogEventEnd  (AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine ComputeEnergy
   
!----------------------------------------------------------------------------------------!      
! ComputeStrainStress (CM)  
!----------------------------------------------------------------------------------------!      
   Subroutine ComputeStrainStress(AppCtx)
     Type(AppCtx_Type)                             :: AppCtx
      
      PetscInt                                     :: iErr
     Type(MatS2D)                                 :: Strain_Elem, Stress_Elem 
     Type(Vect2D)                                 :: F_Elem, U_Elem  
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
    !! I think the volume here is wrong. I moved it in the gauss points cycle
    !             Vol             = Vol + AppCtx%ElemScal(iE)%Gauss_C(iGauss) * AppCtx%ElemScal(iE)%BF(iDoF, iGauss)
            End Do
            Do iGauss = 1, NumGauss
               Do iDoF = 1, NumDoFVect
                  Strain_Elem            = Strain_Elem + AppCtx%ElemVect(iE)%GradS_BF(iDoF, iGauss) * U(iDoF)
                  Stress_Elem            = Stress_Elem + (V_Elem**2) * AppCtx%MatProp( AppCtx%MeshTopology%Elem_Blk(iBlk)%ID )%Hookes_Law * ( Strain_Elem - Theta_Elem * (AppCtx%MatProp(iBlk)%Therm_Exp) )
!             Call PetscLogFlops(3*AppCtx%MeshTopology%Num_Dim+2, iErr)
                  Vol                    = Vol + AppCtx%ElemScal(iE)%Gauss_C(iGauss) 
               End Do
            End Do
            Strain_Elem = Strain_Elem / Vol
            Stress_Elem = Stress_Elem / Vol
            Call PetscLogFlops(AppCtx%MeshTopology%Num_Dim, iErr)
            Stress_Ptr = (/ Stress_Elem%XX, Stress_Elem%YY, Stress_Elem%XY /)
            Strain_Ptr = (/ Strain_Elem%XX, Strain_Elem%YY, Strain_Elem%XY /)
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

End Module m_VarFracFilm_Post
