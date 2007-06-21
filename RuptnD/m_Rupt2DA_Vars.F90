#if defined PB_2D
Module m_Rupt2D_Vars
#elif defined PB_3D
Module m_Rupt3D_Vars
#else 
Module m_Rupt2DA_Vars
#endif

  Use m_MEF90
  Use m_Rupt_Struct

  Implicit NONE
  PRIVATE


#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscsys.h"


  Integer, Public                                       :: SaveInt = 25
  Integer, Public                                       :: BTInt = 5
  
  Mat, Public                                           :: MR_U
  Mat, Public                                           :: MR_V
 

  Vec, Public                                           :: RHS_U
  Vec, Public                                           :: U_Dist, U_Loc, U_Master
  Integer, Public                                       :: UMin, UMax

  Vec, Public                                           :: RHS_V
  Vec, Public                                           :: V_Dist, V_Loc, V_Master
  Vec, Public                                           :: V_Old, V_Change
  Real(Kind = Kr), Public                               :: VMin, VMax, ErrV

  Vec, Public                                           :: BCU_Dist, BCU_Loc, BCU_Master

  Vec, Public                                           :: Temp_Dist, Temp_Loc, Temp_Master

  Vec, Public                                           :: F_Dist, F_Loc, F_Master

  KSP, Public                                           :: KSP_U, KSP_V
  PC, Public                                            :: PC_U, PC_V
  KSPConvergedReason, Public                            :: KSP_Reason
  Integer, Public                                       :: NbIterKSP

#if defined PB_2D
  Type(Node2D), Dimension(:), Pointer, Public           :: Node_db_U
  Type(Element2D_Elast), Dimension(:), Pointer, Public  :: Elem_db_U
  Type(Node2D), Dimension(:), Pointer, Public           :: Node_db_V
  Type(Element2D_Scal), Dimension(:), Pointer, Public   :: Elem_db_V
  Type(Rupt_Params2D), Public                           :: Params
#elif defined PB_3D
  Type(Node3D), Dimension(:), Pointer, Public           :: Node_db_U
  Type(Element3D_Elast), Dimension(:), Pointer, Public  :: Elem_db_U
  Type(Node3D), Dimension(:), Pointer, Public           :: Node_db_V
  Type(Element3D_Scal), Dimension(:), Pointer, Public   :: Elem_db_V
  Type(Rupt_Params3D), Public                           :: Params
#else
  ! 2DA
  Type(Node2D), Dimension(:), Pointer, Public           :: Node_db_U
  Type(Element2D_Scal), Dimension(:), Pointer, Public   :: Elem_db_U
  Type(Node2D), Dimension(:), Pointer, Public           :: Node_db_V
  Type(Element2D_Scal), Dimension(:), Pointer, Public   :: Elem_db_V
  Type(Rupt_Params2D), Public                           :: Params
#endif

  Type(EXO_Geom_Info), Public                           :: Geom
  Type(SD_Info), Public                                 :: MySD_U
  Type(SD_Info), Public                                 :: MySD_V

  Integer, Public                                       :: iErr
  Integer, Public                                       :: iIter
  Integer, Public                                       :: iTS

  Character(len=128), Public                            :: CharBuffer

  PetscLogDouble, Public                                :: TotalTS, TotalTF
  PetscLogDouble, Public                                :: TotalT
  PetscLogDouble, Public                                :: InitTS, InitTF
  PetscLogDouble, Public                                :: AssembTS, AssembTF
  PetscLogDouble, Public                                :: RHSTS, RHSTF
  PetscLogDouble, Public                                :: SolveTS, SolveTF

  Real(Kind = Kr),Dimension(:), Pointer, Public         :: Bulk_Ener
  Real(Kind = Kr),Dimension(:), Pointer, Public         :: Surf_Ener
  Real(Kind = Kr),Dimension(:), Pointer, Public         :: Tot_Ener

  Character(len=128), Public                            :: Ener_Str
  Integer, Parameter, Public                            :: Ener_Unit = 99
  
  Character(len=128), Public                            :: Log_Str
  Integer, Public                                       :: Log_Unit = 98

  KSPConvergedReason, Public                            :: KSP_TestCVG
  Logical, Public                                       :: Is_Interactive
  PetscTruth, Public                                    :: Is_Restarting
  Integer, Public                                       :: TimeStep
  
  Integer, Public                                       :: LogStage_Init
  Integer, Public                                       :: LogStage_IO
  Integer, Public                                       :: LogStage_Assembly
  Integer, Public                                       :: LogStage_Solve
  
#if defined PB_2D
End Module m_Rupt2D_Vars
#elif defined PB_3D
End Module m_Rupt3D_Vars
#else 
End Module m_Rupt2DA_Vars
#endif

  
 
