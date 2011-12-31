#if defined PB_2D
Module m_Poisson2D
#elif defined PB_3D
Module m_Poisson3D
#endif

#include "finclude/petscdef.h"

   Use m_MEF90
   Use m_Heat_Struct 
   
   Implicit NONE   

   
   
Contains

   
   
 

#undef __FUNCT__
#define __FUNCT__ "ComputeEnergy"
   Subroutine ComputeEnergy(AppCtx)
      Type(Heat_AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoF, NumGauss
      PetscReal, Dimension(:), Pointer             :: F, U
      PetscInt                                     :: iBlk, iELoc, iE
      PetscInt                                     :: iDoF, iGauss
#if defined PB_2D
      Type(Vect2D)                                 :: Strain_Elem, Stress_Elem      
#elif defined PB_3D
      Type(Vect3D)                                 :: Strain_Elem, Stress_Elem      
#endif
      PetscReal                                    :: F_Elem, U_Elem
      PetscReal                                    :: MyElasticEnergy, MyExtForcesWork
      PetscLogDouble                               :: flops

      MyElasticEnergy = 0.0_Kr
      MyExtForcesWork = 0.0_Kr
      flops = 0.0_Kr
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            NumDoF   = Size(AppCtx%Elem(iE)%BF,1)
            NumGauss = Size(AppCtx%Elem(iE)%BF,2)
            Allocate(F(NumDoF))
            Call SectionRealRestrictClosure(AppCtx%F%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoF, F, iErr); CHKERRQ(ierr)
            Allocate(U(NumDoF))
            Call SectionRealRestrictClosure(AppCtx%U%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoF, U, iErr); CHKERRQ(ierr)
            Do iGauss = 1, NumGauss
               Strain_Elem = 0.0_Kr
               Stress_Elem = 0.0_Kr
               F_Elem      = 0.0_Kr
               U_Elem      = 0.0_Kr
               Do iDoF = 1, NumDoF
                  Stress_Elem = Stress_Elem + AppCtx%Elem(iE)%Grad_BF(iDoF, iGauss) * U(iDoF)
                  Strain_Elem = Strain_Elem + AppCtx%Elem(iE)%Grad_BF(iDoF, iGauss) * U(iDoF)
                  F_Elem = F_Elem + AppCtx%Elem(iE)%BF(iDoF, iGauss) * F(iDoF)
                  U_Elem = U_Elem + AppCtx%Elem(iE)%BF(iDoF, iGauss) * U(iDoF)
                  flops = flops + 4
               End Do
               MyElasticEnergy = MyElasticEnergy + AppCtx%Elem(iE)%Gauss_C(iGauss) * ( (Strain_Elem .DotP. Strain_Elem) * 0.5_Kr)
               MyExtForcesWork = MyExtForcesWork - AppCtx%Elem(iE)%Gauss_C(iGauss) * F_Elem * U_Elem
               flops = flops + 5
            End Do
            DeAllocate(U)
            DeAllocate(F)
         End Do Do_Elem_iE
      End Do Do_Elem_iBlk
      Call MPI_AllReduce(MyElasticEnergy, AppCtx%ElasticEnergy, 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      Call MPI_AllReduce(MyExtForcesWork, AppCtx%ExtForcesWork, 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      AppCtx%TotalEnergy = AppCtx%ElasticEnergy + AppCtx%ExtForcesWork
   End Subroutine ComputeEnergy
 
#if defined PB_2D
End Module m_Poisson2D
#elif defined PB_3D
End Module m_Poisson3D
#endif
