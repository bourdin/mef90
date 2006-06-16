#ifdef PB_2D
Module m_OvSch2D_Vars
#else
Module m_OvSch3D_Vars
#endif

  Use m_MEF90
  Use m_Poisson_Struct

  Implicit NONE
  PRIVATE

#include "include/finclude/petsc.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petscao.h"

  Character(len=80), Public                             :: Method, PC_Name
  Integer, Public                                       :: Num_Method
  Mat, Public                                           :: MR
  Mat, Public                                           :: MyMR
  Mat, Public                                           :: MyMR_Neu
  Mat, Public                                           :: MyMR_Re
  Vec, Public                                           :: RHS, SOL, CoSOL, CoSOL_Master
  Integer, Public                                       :: MyRank
  Integer, Public                                       :: MySize_Elem
  Integer, Public                                       :: MyIdxMin_Node, MyIdxMax_Node
  Integer, Public                                       :: MySize_Node
  Integer, Public                                       :: NumProcs

  KSP, Public                                           :: KSP_MR
  PC, Public                                            :: PC_MR
  KSPConvergedReason, Public                            :: KSP_Reason
#ifdef PB_2D
  Type(Element2D_Scal), Dimension(:), Pointer, Public   :: Elem_db
  Type(Node2D), Dimension(:), Pointer, Public           :: Node_db
#else
  Type(Element3D_Scal), Dimension(:), Pointer, Public   :: Elem_db
  Type(Node3D), Dimension(:), Pointer, Public           :: Node_db
#endif
  Real(Kind = Kr), Dimension(:), Pointer, Public        :: Disp_Sol
  Real(Kind = Kr), Dimension(:), Pointer, Public        :: Load

  Integer, Public                                       :: GaussOrder=3

  Type(EXO_Geom_Info), Public                           :: Geom
  Type(Poisson_Params), Public                          :: Params
  Integer, Public                                       :: ierr
  Integer, Public                                       :: TimeStep

  VecScatter, Public                                    :: SOL_DistToMaster
  Character(len=128), Public                            :: CharBuffer
  Integer,Dimension(:), Pointer,Public                  :: All_list
  Integer, Public                                       :: L_c
  Integer, Dimension(:), Pointer, Public                :: MyElem
  Integer, Dimension(:), Pointer, Public                :: MyNode
  Integer, Dimension(:), Pointer, Public                :: Elem_Owner
  Integer, Dimension(:), Pointer, Public                :: Node_Owner

  IS, Public                                            :: Mangled_IS
  IS, Public                                            :: Layout_IS
  AO, Public                                            :: DeMangle_AO

  VecScatter, Public                                    :: MangledToMaster
  Vec, Public                                           :: Sol_Master
!  Vec, Public                                           :: Sol_Mangled
  Integer, Dimension(:), Pointer, Public                :: Renum
!  Integer, Dimension(:), Pointer, Public                :: Mangle


  PetscLogDouble, Public                                :: AssembTS, AssembTF
  PetscLogDouble, Public                                :: SolveTS, SolveTF
  PetscLogDouble, Public                                :: TotalTS, TotalTF

!  VecScatter, Public                                    :: Sol_ManToL2 
!  Vec, Public                                           :: Sol_Comp
  Real(Kind = Kr), Dimension(:), Pointer, Public        :: Disp_L2
  Real(Kind = Kr), Dimension(:), Pointer, Public        :: D_Vec
  Real(Kind = Kr), Public                               :: MyL2
  Real(Kind = Kr), Public                               :: MyLIf
  Real(Kind = Kr), Public                               :: FinalL2
  Real(Kind = Kr), Public                               :: MySEMIH2
  Real(Kind = Kr), Public                               :: FinalSEMIH2
  Integer, Public                                       :: Size_Node_Diri_Bd
  Integer, Dimension(:), Pointer, Public                :: Node_Diri_Bd
  Integer, Dimension(:), Pointer, Public                :: MyElem_Ovlp
  Integer, Dimension(:), Pointer, Public                :: MyNode_Ovlp
  Integer, Dimension(:), Pointer, Public                :: MyElem_Ghost
  Integer, Dimension(:), Pointer, Public                :: MyNode_Ghost
  Integer, Dimension(:), Pointer, Public                :: MyNode_Bd
  Integer, Dimension(:), Pointer, Public                :: MyNode_Old
  Real(Kind = Kr), Dimension(:), Pointer , Public       :: MyPU, MyPU2
  Integer, Dimension(:), Pointer, Public                :: MyCount_function
  Integer, Dimension(:), Pointer, Public                :: MyCount_ID_Elem
  Integer, Dimension(:), Pointer, Public                :: MyCount_ID_Node
  Integer, Dimension(:), Pointer, Public                :: MyCount_ID_Elem_Ghost
  Integer, Dimension(:), Pointer, Public                :: MyCount_ID_Node_Ghost
  Integer, Public                                       :: MySize_Elem_Ovlp
  Integer, Public                                       :: MySize_Node_Ovlp
  Integer, Public                                       :: MySize_Elem_Ghost
  Integer, Public                                       :: MySize_Node_Ghost
  Integer, Public                                       :: MySize_Node_METIS
  Integer, Public                                       :: MySize_Node_Bd
  Integer, Public                                       :: MySize_Node_Old
  
  Integer, Dimension(:), Pointer, Public                :: MyCount_ID_Node_METIS
  Integer, Dimension(:), Pointer, Public                :: MyCount_ID_Node_Old
  Integer, Dimension(:), Pointer, Public                :: MyCount_ID_Node_Bd
  Integer, Dimension(:), Pointer, Public                :: MyCount_ID_Node_Ghost_METIS
  Integer, Public                                       :: MySize_Node_Ovlp_METIS
  Integer, Public                                       :: MySize_Node_Ghost_METIS
  Integer, Public                                       :: ovlp
  Integer,Dimension(:), Pointer,Public                  :: MyNode_Ovlp_list
  Integer,Dimension(:), Pointer,Public                  :: MyNode_Bd_Local_list
   Vec, Public                                          :: RHS_Ov, SOL_Ov
   Vec, Public                                          :: RES_Ov, IMP_Ov
   Vec, Public                                          :: MyRHS_Ov, MySOL_Ov
   Vec, Public                                          :: MyRES_Ov, MyIMP_Ov
   Vec, Public                                          :: SOL_Ov_Master
   Vec, Public                                          :: Update_Ov, MyUpdate_Ov
  Integer, Public                                       :: MyIdxMin_Node_Ov, MyIdxMax_Node_Ov
  Integer, Dimension(:), Pointer, Public                :: Renum2
  PetscScalar, Dimension(:), Pointer, Public            :: RHS_Ptr, RES_Ptr
  PetscScalar, Dimension(:), Pointer, Public            :: MySOL_Ov_Ptr
   PetscScalar, Dimension(:), Pointer, Public           :: MyRHS_Ov_Ptr
  PetscScalar, Dimension(:), Pointer, Public            :: MyIMP_Ov_Ptr
  PetscScalar, Dimension(:), Pointer, Public            :: MyRES_Ov_Ptr
  KSP, Public                                           :: KSP_MyMR
  PC, Public                                            :: PC_MyMR
  Real(Kind = Kr), Dimension(:), Pointer, Public        :: Disp_Ov_L2
  Real(Kind = Kr), Dimension(:), Pointer, Public        :: Disp_Sol_Ov
  PetscScalar, Public                                   :: initial


   PetscLogDouble, Public                                :: ASMTS, ASMTF, ASMT, ASMSETUPT
   Integer, Public                                       :: ASMiter

   Integer, Public                                       :: MySize_Neighbor
   Integer, Public                                       :: MySize_Group 
   Integer,Dimension(:), Pointer,Public                  :: MyNeighbor_list 
   Integer,Dimension(:), Pointer,Public                  :: MyGroup_list
   Integer, Dimension(:), Pointer, Public                :: MyNeighbor
   Integer, Dimension(:), Pointer, Public                :: MyGroup
   Integer, Dimension(:), Pointer, Public                :: MyCount_Neighbor

   Integer, Public                                       :: MySize_NeighborTmp
   Integer, Public                                       :: MySize_GroupTmp 
   Integer,Dimension(:), Pointer,Public                  :: MyNeighborTmp_list 
   Integer,Dimension(:), Pointer,Public                  :: MyGroupTmp_list
   Integer, Dimension(:), Pointer, Public                :: MyNeighborTmp
   Integer, Dimension(:), Pointer, Public                :: MyGroupTmp
   Integer, Dimension(:), Pointer, Public                :: MyCount_NeighborTmp
   Integer, Dimension(:), Pointer,Public                 :: TmpNeighbor

   Real(Kind = Kr), Dimension(:,:), Pointer , Public     :: GroupPU
   Mat, Public                                           :: CoMR
   Mat, Public                                           :: MyAdMR
   Integer,Dimension(:), Pointer,Public                  :: PU_node_list
   Integer,Dimension(:), Pointer,Public                  :: PU_node_size
   Integer,Dimension(:), Pointer,Public                  :: PU_node_init

   Vec,Dimension(:), Pointer,Public                      :: PUvec
   Vec,Dimension(:), Pointer,Public                      :: MRPUvec
   Vec,Dimension(:), Pointer,Public                      :: PUvec_Master
   VecScatter,Dimension(:), Pointer,Public               :: PUMaterToMangled
   Integer, Dimension(:), Pointer, Public                :: PUnum

   Mat, Public                                           :: PUCoMR
   Mat, Public                                           :: CoarMR
   Mat, Public                                           :: MyMRT
   Mat, Public                                           :: MRC, MRCO, MRCC

   PetscLogDouble, Public                                :: CoTS, CoTF
   KSP, Public                                           :: KSP_CO
   PC, Public                                            :: PC_CO

   Integer, Dimension(:), Pointer, Public                :: Renum3
   PetscLogDouble, Public                                :: CoSoTS, CoSoTF, CoSoTT
   Vec, Public                                           :: MyPUSec
   Vec,Dimension(:), Pointer,Public                      :: PUVec_Local

   Integer, Dimension(:), Pointer, Public                :: MRCOnum
   Mat, Public                                           :: MRCO_S
   KSP, Public                                           :: KSP_CO_Seq
   PC, Public                                            :: PC_CO_Seq
   VecScatter, Public                                    :: CoVecScatter

   PetscScalar, Dimension(:), Pointer,Public             :: SOL_CoOut_Ptr
   Vec, Public                                           :: CoOut
   Vec, Public                                           :: CoIn_Seq
   Vec, Public                                           :: SOL_CoOut
   Vec, Public                                           :: SOL_CoIn_Seq

   PetscScalar, Public                                   :: FinalRe
   KSPConvergedReason, Public                            :: Reason
   Vec, Public                                           :: CoPU_Ov, MyCoPU_Ov
   PetscScalar, Dimension(:), Pointer,Public             :: MyCoPU_Ov_Ptr
   KSP, Public                                           :: KSP_MyMR_Neu
   PC, Public                                            :: PC_MyMR_Neu
   Vec, Public                                           :: Orth
   MatNullSpace, Public                                  :: NeuNullSpace

   PetscLogDouble, Public                                :: T1S, T1F, TT1
   PetscLogDouble, Public                                :: T2S, T2F, TT2
   PetscLogDouble, Public                                :: T3S, T3F, TT3
   PetscLogDouble, Public                                :: T4S, T4F, TT4
   PetscLogDouble, Public                                :: T5S, T5F, TT5
   PetscLogDouble, Public                                :: T6S, T6F, TT6
   PetscLogDouble, Public                                :: T7S, T7F, TT7
   PetscLogDouble, Public                                :: T8S, T8F, TT8
   PetscLogDouble, Public                                :: T9S, T9F, TT9
   PetscLogDouble, Public                                :: T10S, T10F, TT10
   PetscLogDouble, Public                                :: T11S, T11F, TT11
   PetscLogDouble, Public                                :: T12S, T12F, TT12
   PetscLogDouble, Public                                :: T13S, T13F, TT13
   PetscLogDouble, Public                                :: T14S, T14F, TT14
   PetscLogDouble, Public                                :: T15S, T15F, TT15
   PetscLogDouble, Public                                :: T16S, T16F, TT16
   PetscLogDouble, Public                                :: T17S, T17F, TT17
   PetscLogDouble, Public                                :: T18S, T18F, TT18
   PetscLogDouble, Public                                :: T19S, T19F, TT19
   PetscLogDouble, Public                                :: T20S, T20F, TT20
   PetscLogDouble, Public                                :: T21S, T21F, TT21
   PetscLogDouble, Public                                :: TP1S, TP1F, TP1
   PetscLogDouble, Public                                :: TP2S, TP2F, TP2
   PetscLogDouble, Public                                :: TP3S, TP3F, TP3
   PetscLogDouble, Public                                :: TTTS, TTTF, TTTT
   
!   PetscScalar, Public                              :: zz, vv
   PetscScalar, Public                              :: none
   PetscScalar, Public                              :: pone

 Integer, Dimension(:), Pointer, Public              :: stages

#ifdef PB_2D
End Module m_OvSch2D_Vars
#else
End Module m_OvSch3D_Vars
#endif

