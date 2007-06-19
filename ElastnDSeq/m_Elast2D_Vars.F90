#ifdef PB_2D
Module m_Elast2D_Vars
#else
Module m_Elast3D_Vars
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


  Mat, Public                                           :: MR
  Vec, Public                                           :: RHS

  Vec, Public                                           :: SOL                                                        
  Vec, Public                                           :: BC
  Vec, Public                                           :: F
  Vec, Public                                           :: Temp

  KSP, Public                                           :: KSP_MR
  PC, Public                                            :: PC_MR
  KSPConvergedReason, Public                            :: KSP_Reason

                                                        
#ifdef PB_2D
  Type(Element2D_Elast), Dimension(:), Pointer, Public  :: Elem_db
  Type(Node2D), Dimension(:), Pointer, Public           :: Node_db
  Type(Element2D_Scal), Dimension(:), Pointer, Public   :: Elem_Scal
  Type(Node2D), Dimension(:), Pointer, Public           :: Node_Scal
#else
  Type(Element3D_Elast), Dimension(:), Pointer, Public  :: Elem_db
  Type(Node3D), Dimension(:), Pointer, Public           :: Node_db
  Type(Element3D_Scal), Dimension(:), Pointer, Public   :: Elem_Scal
  Type(Node3D), Dimension(:), Pointer, Public           :: Node_Scal
#endif


  Integer, Public, Parameter                            :: GaussOrder=2

  Type(EXO_Geom_Info), Public                           :: Geom
  Type(Rupt_Params), Public                             :: Params

  Integer, Public                                       :: iErr
  Integer, Public                                       :: TimeStep

  Character(len=128), Public                            :: CharBuffer

  PetscLogDouble, Public                                :: TotalTS, TotalTF
  PetscLogDouble, Public                                :: TotalT

  Real(Kind = Kr), Public                               :: Ener_Elast
  PetscReal, Public, Parameter                          :: VLV = 1.0e+20

#ifdef PB_2D
End Module m_Elast2D_Vars
#else
End Module m_Elast3D_Vars
#endif
