Module m_MEF_Integrate

#include "finclude/petscdef.h"

   Use m_MEF90

IMPLICIT NONE

Interface IntegrateScalLp
! returns an array of lp norm of the V_section by block  
   Module Procedure IntegrateScalLp_2D, IntegrateScalLp_3D
!TODO add IntegrateScalLp_3D, IntegrateVectLp_2D, IntegrateVectLp_3D
End Interface IntegrateScalLp

Contains 
   
#undef __FUNCT__
#define __FUNCT__ "IntegrateScalLp_2D"
   Subroutine IntegrateScalLp_2D(MeshTopology,Elem,V_Sec,IntLp,p)
      Type(MeshTopology_Type)                      :: MeshTopology
      Type(Element2D_Scal), Dimension(:), Pointer  :: Elem
      Type(SectionReal)                            :: V_Sec
      PetscReal, Dimension(:), Pointer             :: IntLp
      PetscInt                                     :: p
      
      PetscInt                                     :: iBlk
      PetscInt                                     :: iErr, NumDoFScal
      PetscInt                                     :: iE, iEloc, IGauss, IDoF1
      PetscReal                                    :: IntRes
      PetscReal, Dimension(:), Pointer             :: MyIntLp
      PetscReal, Dimension(:), Pointer             :: V_Loc

      Allocate(MyIntLp(MeshTopology%Num_Elem_Blks_Global))
      Do iBlk = 1, MeshTopology%Num_Elem_Blks
         NumDoFScal = MeshTopology%Elem_Blk(iBlk)%Num_DoF
         Allocate(V_Loc(NumDoFScal))
         MyIntLp(iBlk) = 0
         Do_iELoc: Do iELoc = 1, MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            Call SectionRealRestrictClosure(V_Sec, MeshTopology%mesh, iE-1, NumDoFScal, V_Loc, iErr); CHKERRQ(ierr)
            Do iGauss = 1, size(Elem(iE)%Gauss_C)
            IntRes = 0
               Do iDoF1 = 1, NumDoFScal
                     IntRes = IntRes   + V_Loc(iDoF1) * Elem(iE)%BF(iDoF1, iGauss)
               End Do
               MyIntLp(iBlk) = MyIntLp(iBlk) + Elem(iE)%Gauss_C(iGauss)*(Intres**p)
            End DO
         End Do Do_IEloc 
         DeAllocate(V_Loc)
      End DO
      Call MPI_AllReduce(MyIntLp, IntLp, MeshTopology%Num_Elem_Blks_Global, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      Do iBlk = 1, MeshTopology%Num_Elem_Blks
         IntLp(iBlk) =  IntLp(Iblk)**(1/p)
      End Do
      DeAllocate(MyIntLp)
   End Subroutine IntegrateScalLp_2D

#undef __FUNCT__
#define __FUNCT__ "IntegrateScalLp_3D"
   Subroutine IntegrateScalLp_3D(MeshTopology,Elem,V_Sec,IntLp,p)
      Type(MeshTopology_Type)                      :: MeshTopology
      Type(Element3D_Scal), Dimension(:), Pointer  :: Elem
      Type(SectionReal)                            :: V_Sec
      PetscReal, Dimension(:), Pointer             :: IntLp
      PetscInt                                     :: p
      
      PetscInt                                     :: iBlk
      PetscInt                                     :: iErr, NumDoFScal
      PetscInt                                     :: iE, iEloc, IGauss, IDoF1
      PetscReal                                    :: IntRes
      PetscReal, Dimension(:), Pointer             :: MyIntLp
      PetscReal, Dimension(:), Pointer             :: V_Loc

      Allocate(MyIntLp(MeshTopology%Num_Elem_Blks_Global))
      Do iBlk = 1, MeshTopology%Num_Elem_Blks
         NumDoFScal = MeshTopology%Elem_Blk(iBlk)%Num_DoF
         Allocate(V_Loc(NumDoFScal))
         MyIntLp(iBlk) = 0
         Do_iELoc: Do iELoc = 1, MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            Call SectionRealRestrictClosure(V_Sec, MeshTopology%mesh, iE-1, NumDoFScal, V_Loc, iErr); CHKERRQ(ierr)
            Do iGauss = 1, size(Elem(iE)%Gauss_C)
            IntRes = 0
               Do iDoF1 = 1, NumDoFScal
                     IntRes = IntRes   + V_Loc(iDoF1) * Elem(iE)%BF(iDoF1, iGauss)
               End Do
               MyIntLp(iBlk) = MyIntLp(iBlk) + Elem(iE)%Gauss_C(iGauss)*(Intres**p)
            End DO
         End Do Do_IEloc 
         DeAllocate(V_Loc)
      End DO
      Call MPI_AllReduce(MyIntLp, IntLp, MeshTopology%Num_Elem_Blks_Global, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      Do iBlk = 1, MeshTopology%Num_Elem_Blks
         IntLp(iBlk) =  IntLp(Iblk)**(1/p)
      End Do
      DeAllocate(MyIntLp)
   End Subroutine IntegrateScalLp_3D

#undef __FUNCT__
#define __FUNCT__ "IntegrateScalSemiH1"
   Subroutine IntegrateScalSemiH1(MeshTopology,iBlk,Elem,V_Sec,IntSemiH1)
      Type(SectionReal)                            :: V_Sec
      Type(MeshTopology_Type)                      :: MeshTopology
      PetscInt                                     :: iBlk
      PetscInt                                     :: iErr, NumDoFScal
      PetscInt                                     :: iE, iEloc, IGauss, IDoF1
      PetscReal                                    :: IntSemiH1
      PetscReal, Dimension(:), Pointer             :: V_Loc
#if defined PB_2D
      Type(Element2D_Scal), Dimension(:), Pointer  :: Elem
       Type(Vect2D)                                 :: GradV_Elem
#elif defined PB_3D 
      Type(Element3D_Scal), Dimension(:), Pointer  :: Elem
      Type(Vect3D)                                 :: GradV_Elem
#endif   

      NumDoFScal = MeshTopology%Elem_Blk(iBlk)%Num_DoF
      Allocate(V_Loc(NumDoFScal))
      GradV_Elem = 0.0_Kr
      IntSemiH1 = 0.0_Kr
      Do_iELoc: Do iELoc = 1, MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         Call SectionRealRestrictClosure(V_Sec, MeshTopology%mesh, iE-1, NumDoFScal, V_Loc, iErr); CHKERRQ(ierr)
         Do iGauss = 1, size(Elem(iE)%Gauss_C)
            Do iDoF1 = 1, NumDoFScal
                  GradV_Elem = GradV_Elem  + (V_Loc(iDoF1) * Elem(iE)%Grad_BF(iDoF1, iGauss))
            End Do
         IntSemiH1 = IntSemiH1 + (GradV_Elem .DotP. GradV_Elem)
         End DO
      End Do Do_IEloc 
      IntSemiH1 = sqrt(IntSemiH1)
   DeAllocate(V_Loc)

   End Subroutine IntegrateScalSemiH1

End Module m_MEF_Integrate
