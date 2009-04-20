#if defined PB_2D
module m_AmbrosioTortorelli2D_Post
#elif defined PB_3D
Module m_AmbrosioTortorelli3D_Post
#endif

#include "finclude/petscdef.h"
#include "finclude/petscmeshdef.h"
#include "finclude/petscviewerdef.h"

#if defined PB_2D
   Use m_AmbrosioTortorelli_Types2D
   Use m_AmbrosioTortorelli_U2D
   Use m_AmbrosioTortorelli_V2D
#elif defined PB_3D
   Use m_AmbrosioTortorelli_Types3D   
   Use m_AmbrosioTortorelli_U3D
   Use m_AmbrosioTortorelli_V3D
#endif   
   Use m_MEF90
   Use m_RuptStruct
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
      PetscReal, Dimension(:), Pointer             :: F, U, Theta
      PetscInt                                     :: iBlk, iELoc, iE
      PetscInt                                     :: iDoF, iGauss
      PetscReal                                    :: Theta_Elem
#if defined PB_2D
     Type(MatS2D)                                  :: Strain_Elem, Stress_Elem 
     Type(MatS2D)                                  :: Effective_Strain_Elem
     Type(Vect2D)                                  :: F_Elem, U_Elem  

#elif defined PB_3D
     Type(MatS3D)                                  :: Strain_Elem, Stress_Elem 
     Type(MatS3D)                                  :: Effective_Strain_Elem
     Type(Vect3D)                                  :: F_Elem, U_Elem  
#endif
      PetscReal                                    :: MyBulkEnergy
     

      Call PetscLogStagePush(AppCtx%LogInfo%PostProc_Stage, iErr); CHKERRQ(iErr)
      Call PetscLogEventBegin(AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      
      MyBulkEnergy = 0.0_Kr
         
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            NumDoFVect = Size(AppCtx%ElemVect(iE)%BF,1)
            NumDoFScal = Size(AppCtx%ElemScal(iE)%BF,1)
            NumGauss   = Size(AppCtx%ElemVect(iE)%BF,2)
            Allocate(F(NumDoFVect))
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%F, iE-1, NumDoFVect, F, iErr); CHKERRQ(ierr)
            Allocate(U(NumDoFVect))
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%U, iE-1, NumDoFVect, U, iErr); CHKERRQ(ierr)
            Allocate(Theta(NumDoFScal))
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%Theta, iE-1, NumDoFScal, Theta, iErr); CHKERRQ(ierr)
            Do iGauss = 1, NumGauss
               Strain_Elem           = 0.0_Kr
               Stress_Elem           = 0.0_Kr
               Theta_Elem            = 0.0_Kr
               Effective_Strain_Elem = 0.0_Kr
               F_Elem                = 0.0_Kr
               U_Elem                = 0.0_Kr
               Do iDof = 1, NumDoFScal
                  Theta_Elem = Theta_Elem + AppCtx%ElemScal(iE)%BF(iDoF, iGauss) * Theta(iDoF)
               End Do
               Do iDoF = 1, NumDoFVect
                  Strain_Elem = Strain_Elem + AppCtx%ElemVect(iE)%GradS_BF(iDoF, iGauss) * U(iDoF)
                  F_Elem     = F_Elem + F(iDoF) * AppCtx%ElemVect(iE)%BF(iDoF, iGauss) 
                  U_Elem     = U_Elem + U(iDoF) * AppCtx%ElemVect(iE)%BF(iDoF, iGauss) 
               End Do            
               Effective_Strain_Elem  = Strain_Elem - Theta_Elem * AppCtx%MatProp(iBlk)%Therm_Exp   ;
               Stress_Elem            = AppCtx%MatProp( AppCtx%MeshTopology%Elem_Blk(iBlk)%ID )%Hookes_Law * Effective_Strain_Elem
               MyBulkEnergy           = MyBulkEnergy + AppCtx%ElemVect(iE)%Gauss_C(iGauss) * ( (Stress_Elem .DotP. Effective_Strain_Elem) * 0.5_Kr - (F_Elem .DotP. U_Elem))
            Call PetscLogFlops(AppCtx%MeshTopology%Num_Dim+4, iErr)
            End Do
            DeAllocate(F)
            DeAllocate(U)
            DeAllocate(Theta)
         End Do Do_Elem_iE
      End Do Do_Elem_iBlk
      Call PetscGlobalSum(MyBulkEnergy, AppCtx%BulkEnergy, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
               
      Call PetscLogEventEnd  (AppCtx%LogInfo%PostProc_Event, iErr); CHKERRQ(iErr)
      Call PetscLogStagePop(iErr); CHKERRQ(iErr)
   End Subroutine ComputeEnergy
   
!----------------------------------------------------------------------------------------!      
! ComputeStrainStress (CM)  
!----------------------------------------------------------------------------------------!      
   Subroutine ComputeStrainStress(AppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
#if defined PB_2D
      Type(MatS2D)                                 :: Strain_Elem, Stress_Elem 
      Type(MatS2D)                                 :: Effective_Strain_Elem
      Type(Vect2D)                                 :: F_Elem, U_Elem  
#elif defined PB_3D 
      Type(MatS3D)                                 :: Strain_Elem, Stress_Elem 
      Type(MatS3D)                                 :: Effective_Strain_Elem
      Type(Vect3D)                                 :: F_Elem, U_Elem  
#endif
      PetscReal                                    :: Theta_Elem
      PetscReal                                    :: Vol
      PetscInt                                     :: NumDoFVect, NumDofScal, NumGauss
      PetscReal, Dimension(:), Pointer             :: U, Theta
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
            NumDoFScal = Size(AppCtx%ElemScal(iE)%BF,1)
            NumDoFVect = Size(AppCtx%ElemVect(iE)%BF,1)
            NumGauss   = Size(AppCtx%ElemVect(iE)%BF,2)
            Allocate(U(NumDoFVect))
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%U, iE-1, NumDoFVect, U, iErr); CHKERRQ(ierr)
            Allocate(Theta(NumDoFScal))
            Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, AppCtx%Theta, iE-1, NumDoFScal, Theta, iErr); CHKERRQ(ierr)

            Strain_Elem           = 0.0_Kr
            Stress_Elem           = 0.0_Kr
            Theta_Elem            = 0.0_Kr
            Effective_Strain_Elem = 0.0_Kr
            Vol  = 0.0_Kr
            Do iDof = 1, NumDoFScal
                  Theta_Elem = Theta_Elem + AppCtx%ElemScal(iE)%BF(iDoF, iGauss) * Theta(iDoF)
                  Vol = Vol + AppCtx%ElemScal(iE)%Gauss_C(iGauss) * AppCtx%ElemScal(iE)%BF(iDoF, iGauss)
               End Do
            Do iGauss = 1, NumGauss
               Do iDoF = 1, NumDoFVect
                  Strain_Elem = Strain_Elem + AppCtx%ElemVect(iE)%GradS_BF(iDoF, iGauss) * U(iDoF)
!             Call PetscLogFlops(3*AppCtx%MeshTopology%Num_Dim+2, iErr)
               End Do
            End Do
            Strain_Elem = Strain_Elem / Vol
            Effective_Strain_Elem  = Strain_Elem - Theta_Elem * (AppCtx%MatProp(iBlk)%Therm_Exp)   ;
            Stress_Elem = AppCtx%MatProp( AppCtx%MeshTopology%Elem_Blk(iBlk)%ID )%Hookes_Law * Effective_Strain_Elem
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
End Module m_AmbrosioTortorelli2D_Post
#elif defined PB_3D
End Module m_AmbrosioTortorelli3D_Post
#endif