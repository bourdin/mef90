Module m_Elast2DA_Vars
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
   Vec, Public                                           :: RHS, RHS_Loc
   Vec, Public                                           :: IO_N
   
   Vec, Public                                           :: SOL
                                                        
   Vec, Public                                           :: BC, BC_Local
   
   Vec, Public                                           :: F, F_Local
   
   KSP, Public                                           :: KSP_MR
   PC, Public                                            :: PC_MR
   KSPConvergedReason, Public                            :: KSP_Reason
   
                                                        
   Type(Element2D_Scal), Dimension(:), Pointer, Public   :: Elem_db
   Type(Node2D), Dimension(:), Pointer, Public           :: Node_db
   
   
   Integer, Public                                       :: GaussOrder=2
   
   Type(EXO_Geom_Info), Public                           :: Geom, MyGeom
   Type(Rupt_Params2D), Public                           :: Params
   Type (Layout_Info), Public                            :: Layout
   
   
   Integer, Public                                       :: iErr
   Integer, Public                                       :: TimeStep
   
   Character(len=128), Public                            :: CharBuffer
   
   PetscLogDouble, Public                                :: TotalTS, TotalTF
   PetscLogDouble, Public                                :: TotalT
   PetscLogDouble, Public                                :: InitTS, InitTF
   PetscLogDouble, Public                                :: AssembTS, AssembTF
   PetscLogDouble, Public                                :: RHSTS, RHSTF
   PetscLogDouble, Public                                :: SolveTS, SolveTF
   PetscLogDouble, Public                                :: ExportTS, ExportTF
   
   Real(Kind = Kr), Public                               :: Ener_Elast
   PetscReal, Public, Parameter                          :: VLV = 1.0e+30
End Module m_Elast2DA_Vars
