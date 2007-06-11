#ifdef PB_2D
Module m_Elast2D_Debug
#else
Module m_Elast3D_Debug
#endif
  Use m_MEF90
  Use m_Rupt_Struct
  
#ifdef PB_2D
  Use m_Elast2D_Vars
#else
  Use m_Elast3D_Vars
#endif  
  Implicit NONE
  PRIVATE
  
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
  
  Public :: CompInt

  
Contains
  Subroutine CompInt(Geom, Elems, Nodes)
    Type(EXO_Geom_Info)                                 :: Geom
#ifdef PB_2D
    Type(Element2D_Elast), Dimension(:), Pointer        :: Elems
    Type(Node2D), DImension(:), Pointer                 :: Nodes
#else
    Type(Element3D_Elast), Dimension(:), Pointer        :: Elems
    Type(Node3D), DImension(:), Pointer                 :: Nodes
#endif
    Integer                                             :: iBlk, iE, iEBlk
    Integer                                             :: iSL, iSG
    Integer                                             :: iG
    Real(Kind = Kr)                                     :: Vol_Blk

    DoiBlk: Do iBlk = 1, Geom%Num_Elem_Blks
       Vol_Blk = 0.0_Kr
       DoiEBlk: Do iEBlk = 1, Geom%Elem_Blk(iBlk)%Num_Elems
          iE = Geom%Elem_blk(iBlk)%Elem_ID(iEBlk)
          DoiG: Do iG = 1, Elems(iE)%Nb_Gauss
             Vol_Blk = Vol_Blk + Elems(iE)%Gauss_C(iG)
          End Do DoiG
       End Do DoiEBlk
       Write(*,100) iBlk, Vol_Blk
    End Do DoiBlk

100 Format ('Block: ', I4, 'Volume/Area: ', ES12.5)
  End Subroutine CompInt

!!! This is not updated
!!$
!!$  Subroutine Test_EigenValues(Geom, Elems, Nodes)
!!$    Type(EXO_Geom_Info)                                 :: Geom
!!$    Type(Element3D_Elast), Dimension(:), Pointer        :: Elems
!!$    Type(Node3D), Dimension(:), Pointer                 :: Nodes
!!$    
!!$    Integer                                             :: Nb_Nodes, iS
!!$
!!$    Nb_Nodes = Size(Nodes)/3
!!$
!!$
!!$    TimeStep = 1
!!$    Write(*,*) 'Step 1: Uniform compression'
!!$    Do iS = 1, Nb_Nodes
!!$       BC(iS) = Nodes(iS)%Coord%X
!!$       BC(iS + Nb_Nodes) = Nodes(iS)%Coord%Y
!!$       BC(iS + 2 * Nb_Nodes) = Nodes(iS)%Coord%Z
!!$    End Do
!!$    Call Assemb_RHS_Elast(RHS, Geom, Params, Elems, Nodes, BC)
!!$    Call KSPSolve(KSP_MR, iErr)
!!$    Call Export()
!!$
!!$    TimeStep = TimeStep + 1
!!$    Write(*,*) 'Step 2: Shear 1'
!!$    Do iS = 1, Nb_Nodes
!!$       BC(iS) = Nodes(iS)%Coord%Y
!!$       BC(iS + Nb_Nodes) = Nodes(iS)%Coord%X
!!$       BC(iS + 2 * Nb_Nodes) = 0.0_Kr
!!$    End Do
!!$    Call Assemb_RHS_Elast(RHS, Geom, Params, Elems, Nodes, BC)
!!$    Call KSPSolve(KSP_MR, iErr)
!!$    Call Export()
!!$
!!$    TimeStep = TimeStep + 1
!!$    Write(*,*) 'Step 3: Shear 2'
!!$    Do iS = 1, Nb_Nodes
!!$       BC(iS) = Nodes(iS)%Coord%Z
!!$       BC(iS + Nb_Nodes) = 0.0_Kr
!!$       BC(iS + 2 * Nb_Nodes) = Nodes(iS)%Coord%X
!!$    End Do
!!$    Call Assemb_RHS_Elast(RHS, Geom, Params, Elems, Nodes, BC)
!!$    Call KSPSolve(KSP_MR, iErr)
!!$    Call Export()
!!$
!!$    TimeStep = TimeStep + 1
!!$    Write(*,*) 'Step 4: Shear 3'
!!$    Do iS = 1, Nb_Nodes
!!$       BC(iS) = 0.0_Kr
!!$       BC(iS + Nb_Nodes) = Nodes(iS)%Coord%Z
!!$       BC(iS + 2 * Nb_Nodes) = Nodes(iS)%Coord%Y
!!$    End Do
!!$    Call Assemb_RHS_Elast(RHS, Geom, Params, Elems, Nodes, BC)
!!$    Call KSPSolve(KSP_MR, iErr)
!!$    Call Export()
!!$
!!$    TimeStep = TimeStep + 1
!!$    Write(*,*) 'Step 5: Shear 4'
!!$    Do iS = 1, Nb_Nodes
!!$       BC(iS) = Nodes(iS)%Coord%X
!!$       BC(iS + Nb_Nodes) = 0.0_Kr
!!$       BC(iS + 2 * Nb_Nodes) = -Nodes(iS)%Coord%Z
!!$    End Do
!!$    Call Assemb_RHS_Elast(RHS, Geom, Params, Elems, Nodes, BC)
!!$    Call KSPSolve(KSP_MR, iErr)
!!$    Call Export()
!!$
!!$    TimeStep = TimeStep + 1
!!$    Write(*,*) 'Step 6: Shear 5'
!!$    Do iS = 1, Nb_Nodes
!!$       BC(iS) = Nodes(iS)%Coord%X
!!$       BC(iS + Nb_Nodes) = -Nodes(iS)%Coord%Y
!!$       BC(iS + 2 * Nb_Nodes) = 0.0_Kr
!!$    End Do
!!$    Call Assemb_RHS_Elast(RHS, Geom, Params, Elems, Nodes, BC)
!!$    Call KSPSolve(KSP_MR, iErr)
!!$    Call Export()
!!$  End Subroutine Test_EigenValues
!!$
!!$  Subroutine Test_Hookes_Law()
!!$    Integer      :: iS
!!$    PetscScalar  :: value
!!$    
!!$  TimeStep = 1
!!$  Write(*,*)
!!$  Write(*,*) 'TIMESTEP ', TimeStep
!!$  Call Vecset(0.0_Kr, Sol, iErr)
!!$  Do iS = 1, Size(Node_db)/3
!!$     value = Node_db(is)%Coord%X
!!$     Call VecSetValue(Sol, iS-1, value, INSERT_VALUES, iErr)
!!$  End Do
!!$  Call VecAssemblyBegin(Sol, iErr)
!!$  Call VecAssemblyEnd(Sol, iErr)
!!$  Call Export()
!!$
!!$  TimeStep = 2
!!$  Write(*,*)
!!$  Write(*,*) 'TIMESTEP ', TimeStep
!!$  Call Vecset(0.0_Kr, Sol, iErr)
!!$  Do iS = 1, Size(Node_db)/3
!!$     value = Node_db(is)%Coord%Y
!!$     Call VecSetValue(Sol, iS-1, value, INSERT_VALUES, iErr)
!!$  End Do
!!$  Call VecAssemblyBegin(Sol, iErr)
!!$  Call VecAssemblyEnd(Sol, iErr)
!!$  Call Export()
!!$
!!$  TimeStep = 3
!!$  Write(*,*)
!!$  Write(*,*) 'TIMESTEP ', TimeStep
!!$  Call Vecset(0.0_Kr, Sol, iErr)
!!$  Do iS = 1, Size(Node_db)/3
!!$     value = Node_db(is)%Coord%Z
!!$     Call VecSetValue(Sol, iS-1, value, INSERT_VALUES, iErr)
!!$  End Do
!!$  Call VecAssemblyBegin(Sol, iErr)
!!$  Call VecAssemblyEnd(Sol, iErr)
!!$  Call Export()
!!$
!!$  TimeStep = 4
!!$  Write(*,*)
!!$  Write(*,*) 'TIMESTEP ', TimeStep
!!$  Call Vecset(0.0_Kr, Sol, iErr)
!!$  Do iS = Size(Node_db)/3 + 1, 2 * Size(Node_db)/3
!!$     value = Node_db(is)%Coord%X
!!$     Call VecSetValue(Sol, iS-1, value, INSERT_VALUES, iErr)
!!$  End Do
!!$  Call VecAssemblyBegin(Sol, iErr)
!!$  Call VecAssemblyEnd(Sol, iErr)
!!$  Call Export()
!!$
!!$  TimeStep = 5
!!$  Write(*,*)
!!$  Write(*,*) 'TIMESTEP ', TimeStep
!!$  Call Vecset(0.0_Kr, Sol, iErr)
!!$  Do iS = Size(Node_db)/3 + 1, 2 * Size(Node_db)/3
!!$     value = Node_db(is)%Coord%Y
!!$     Call VecSetValue(Sol, iS-1, value, INSERT_VALUES, iErr)
!!$  End Do
!!$  Call VecAssemblyBegin(Sol, iErr)
!!$  Call VecAssemblyEnd(Sol, iErr)
!!$  Call Export()
!!$
!!$  TimeStep = 6
!!$  Write(*,*)
!!$  Write(*,*) 'TIMESTEP ', TimeStep
!!$  Call Vecset(0.0_Kr, Sol, iErr)
!!$  Do iS = Size(Node_db)/3 + 1, 2 * Size(Node_db)/3
!!$     value = Node_db(is)%Coord%Z
!!$     Call VecSetValue(Sol, iS-1, value, INSERT_VALUES, iErr)
!!$  End Do
!!$  Call VecAssemblyBegin(Sol, iErr)
!!$  Call VecAssemblyEnd(Sol, iErr)
!!$  Call Export()
!!$
!!$  TimeStep = 7
!!$  Write(*,*)
!!$  Write(*,*) 'TIMESTEP ', TimeStep
!!$  Call Vecset(0.0_Kr, Sol, iErr)
!!$  Do iS = 2 * Size(Node_db)/3 + 1, Size(Node_db)
!!$     value = Node_db(is)%Coord%X
!!$     Call VecSetValue(Sol, iS-1, value, INSERT_VALUES, iErr)
!!$  End Do
!!$  Call VecAssemblyBegin(Sol, iErr)
!!$  Call VecAssemblyEnd(Sol, iErr)
!!$  Call Export()
!!$
!!$  TimeStep = 8
!!$  Write(*,*)
!!$  Write(*,*) 'TIMESTEP ', TimeStep
!!$  Call Vecset(0.0_Kr, Sol, iErr)
!!$  Do iS = 2 * Size(Node_db)/3 + 1, Size(Node_db)
!!$     value = Node_db(is)%Coord%Y
!!$     Call VecSetValue(Sol, iS-1, value, INSERT_VALUES, iErr)
!!$  End Do
!!$  Call VecAssemblyBegin(Sol, iErr)
!!$  Call VecAssemblyEnd(Sol, iErr)
!!$  Call Export()
!!$
!!$  TimeStep = 9
!!$  Write(*,*)
!!$  Write(*,*) 'TIMESTEP ', TimeStep
!!$  Call Vecset(0.0_Kr, Sol, iErr)
!!$  Do iS = 2 * Size(Node_db)/3 + 1, Size(Node_db)
!!$     value = Node_db(is)%Coord%Z
!!$     Call VecSetValue(Sol, iS-1, value, INSERT_VALUES, iErr)
!!$  End Do
!!$  Call VecAssemblyBegin(Sol, iErr)
!!$  Call VecAssemblyEnd(Sol, iErr)
!!$  Call Export()
!!$
!!$End Subroutine Test_Hookes_Law


#ifdef PB_2D
End Module m_Elast2D_Debug
#else
End Module m_Elast3D_Debug
#endif
