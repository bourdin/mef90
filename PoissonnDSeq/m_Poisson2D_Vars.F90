#if defined PB_2D
Module m_Poisson2D_Vars
#elif defined PB_3D
Module m_Poisson3D_Vars
#endif

  Use m_MEF90
  Use m_Poisson_Struct

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


  Integer, Public                                       :: SaveInt = 5
  Mat, Public                                           :: MR

  Vec, Public                                           :: RHS_U
  Vec, Public                                           :: U_Dist, U_Loc, U_Master
  Real(Kind = Kr), Dimension(:), Pointer, Public        :: U_Ptr
  Integer, Public                                       :: UMin, UMax

  Vec, Public                                           :: Load_Dist, Load_Loc, Load_Master
  Real(Kind = Kr), Dimension(:), Pointer, Public        :: Load_Ptr

  KSP, Public                                           :: KSP_U
  PC, Public                                            :: PC_U
  KSPConvergedReason, Public                            :: KSP_Reason

#if defined PB_2D
  Type(Node2D), Dimension(:), Pointer, Public           :: Node_db
  Type(Element2D_Scal), Dimension(:), Pointer, Public   :: Elem_db
#elif defined PB_3D
  Type(Node3D), Dimension(:), Pointer, Public           :: Node_db
  Type(Element3D_Scal), Dimension(:), Pointer, Public   :: Elem_db
#endif

  Integer, Public                                       :: GaussOrder = 2

  Type(EXO_Geom_Info), Public                           :: Geom
  Type(Poisson_Params), Public                          :: Params
  Type(SD_Info), Public                                 :: MySD

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

  Real(Kind = Kr), Public                               :: Ener

  PetscReal, Parameter, Public                          :: VLV = 1.0e+20

  Integer, Public                                       :: NbIter
#if defined PB_2D
End Module m_Poisson2D_Vars
#elif defined PB_3D
End Module m_Poisson3D_Vars
#endif

  
