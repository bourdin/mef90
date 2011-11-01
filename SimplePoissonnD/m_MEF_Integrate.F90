Module m_MEF_Integrate

#include "finclude/petscdef.h"

   Use m_MEF90

IMPLICIT NONE

Contains 
   
#undef __FUNCT__
#define __FUNCT__ "IntegrateScalL2"
   Subroutine IntegrateScalL2(MeshTopology,iBlk,Elem,V_Sec,IntL1)
      Type(SectionReal)                            :: V_Sec
      Type(MeshTopology_Type)                      :: MeshTopology
      PetscInt                                     :: iBlk
      PetscInt                                     :: iErr, NumDoFScal
      PetscInt                                     :: iE, iEloc, IGauss, IDoF1
      PetscReal                                    :: IntRes, IntL1
      PetscReal, Dimension(:), Pointer             :: V_Loc
#if defined PB_2D
      Type(Element2D_Scal), Dimension(:), Pointer  :: Elem
#elif defined PB_3D 
      Type(Element3D_Scal), Dimension(:), Pointer  :: Elem
#endif   

      NumDoFScal = MeshTopology%Elem_Blk(iBlk)%Num_DoF
      Allocate(V_Loc(NumDoFScal))
      IntL1 = 0
      Do_iELoc: Do iELoc = 1, MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         Call SectionRealRestrictClosure(V_Sec, MeshTopology%mesh, iE-1, NumDoFScal, V_Loc, iErr); CHKERRQ(ierr)
         Do iGauss = 1, size(Elem(iE)%Gauss_C)
         IntRes = 0
            Do iDoF1 = 1, NumDoFScal
                  IntRes = IntRes     + V_Loc(iDoF1) *   Elem(iE)%BF(iDoF1, iGauss)
            End Do
            IntL1 = IntL1 + Elem(iE)%Gauss_C(iGauss)*(Intres**2) 
         End DO
      End Do Do_IEloc 
      IntL1 = sqrt(IntL1)
   DeAllocate(V_Loc)

   End Subroutine IntegrateScalL2
   
#undef __FUNCT__
#define __FUNCT__ "IntegrateScalL1"
   Subroutine IntegrateScalL1(MeshTopology,iBlk,Elem,V_Sec,IntL1)
      Type(SectionReal)                            :: V_Sec
      Type(MeshTopology_Type)                      :: MeshTopology
      PetscInt                                     :: iBlk
      PetscInt                                     :: iErr, NumDoFScal
      PetscInt                                     :: iE, iEloc, IGauss, IDoF1
      PetscReal                                    :: IntRes, IntL1
      PetscReal, Dimension(:), Pointer             :: V_Loc
#if defined PB_2D
      Type(Element2D_Scal), Dimension(:), Pointer  :: Elem
#elif defined PB_3D 
      Type(Element3D_Scal), Dimension(:), Pointer  :: Elem
#endif   

      NumDoFScal = MeshTopology%Elem_Blk(iBlk)%Num_DoF
      Allocate(V_Loc(NumDoFScal))
      IntL1 = 0
      Do_iELoc: Do iELoc = 1, MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         Call SectionRealRestrictClosure(V_Sec, MeshTopology%mesh, iE-1, NumDoFScal, V_Loc, iErr); CHKERRQ(ierr)
         Do iGauss = 1, size(Elem(iE)%Gauss_C)
         IntRes = 0
            Do iDoF1 = 1, NumDoFScal
                  IntRes = IntRes     + V_Loc(iDoF1) * Elem(iE)%BF(iDoF1, iGauss)
            End Do
            IntL1 = IntL1 + Elem(iE)%Gauss_C(iGauss)*Intres 
         End DO
      End Do Do_IEloc 
      IntL1 = sqrt(IntL1)
   DeAllocate(V_Loc)

   End Subroutine IntegrateScalL1

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
