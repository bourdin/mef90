Program TestRuptStruct

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscmeshdef.h"

   Use m_MEF90
   Use m_RuptStruct

   Use petsc
   Use petscvec
   Use petscmat
   Use petscksp
   Use petscmesh

   Implicit NONE   
   
   Character(len=MEF90_MXSTRLEN)                :: prefix
   Type(MeshTopology_Type)                      :: MeshTopology
   Type(EXO_Type)                               :: EXO, MyEXO
   Type(Mesh)                                   :: Tmp_Mesh
   PetscTruth                                   :: HasPrefix
   PetscInt                                     :: iErr, i
   Type(EXO_RuptProperties_Type)                :: EXO_RuptProperties
   
   Call MEF90_Initialize()
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p',       prefix, HasPrefix, iErr); CHKERRQ(iErr)

   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'

   Call MeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
   Call MeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, MeshTopology%mesh, ierr); CHKERRQ(iErr)
   Call MeshDestroy(Tmp_mesh, iErr); CHKERRQ(iErr)

   Call MeshTopologyReadEXO(MeshTopology, EXO)
   Do i = 1, MeshTopology%Num_Elem_Blks
      MeshTopology%Elem_Blk(i)%Elem_Type = MEF90_P1_Lagrange
      Call Init_Elem_Blk_Type(MeshTopology%Elem_Blk(i), MeshTopology%num_dim)
   End Do

      
   MyEXO%comm = PETSC_COMM_SELF
   MyEXO%exoid = EXO%exoid
   
   Write(MyEXO%filename, 200) trim(prefix), MEF90_MyRank
   200 Format(A, '-', I4.4, '.gen')

   Call Write_MeshTopologyGlobal(MeshTopology, MyEXO, PETSC_COMM_WORLD)
   
!!!   MyEXO%exoid = EXOPEN(MyEXO%filename, EXWRIT, exo_cpu_ws, exo_io_ws, exo_ver, ierr)
!!!   Do i = 1, MeshTopology%Num_Elem_Blks
!!!      Call EXPP(MyEXO%exoid, EXEBLK, MeshTopology%Elem_Blk(i)%ID, 'MeshTopology%Elem_Blk(i)%ID', MeshTopology%Elem_Blk(i)%ID, iErr)
!!!      Call EXPP(MyEXO%exoid, EXEBLK, MeshTopology%Elem_Blk(i)%ID, 'i', i, iErr)
!!!   End Do
!!!   Do i = 1, MeshTopology%Num_Node_Sets
!!!      Call EXPP(MyEXO%exoid, EXNSET, MeshTopology%Node_Set(i)%ID, 'MeshTopology%Node_Set(i)%ID', MeshTopology%Node_Set(i)%ID, iErr)
!!!      Call EXPP(MyEXO%exoid, EXNSET, MeshTopology%Node_Set(i)%ID, 'i', i, iErr)
!!!   End Do
!!!   
!!!   Call EXCLOS(MyEXO%exoid, iErr)
!!!   MyEXO%exoid = 0
!   Call EXO_RuptFormat(MyEXO)

   Allocate(EXO_RuptProperties%Is_Brittle(MeshTopology%Num_Elem_Blks))
   Allocate(EXO_RuptProperties%Is_Domain(MeshTopology%Num_Elem_Blks))
   Allocate(EXO_RuptProperties%Has_BodyForce(MeshTopology%Num_Elem_Blks))
   Allocate(EXO_RuptProperties%EB_BC_Type_X(MeshTopology%Num_Elem_Blks))
   Allocate(EXO_RuptProperties%EB_BC_Type_Y(MeshTopology%Num_Elem_Blks))
   Allocate(EXO_RuptProperties%EB_BC_Type_Z(MeshTopology%Num_Elem_Blks))
   Allocate(EXO_RuptProperties%SS_BC_Type_X(MeshTopology%Num_Side_Sets))
   Allocate(EXO_RuptProperties%SS_BC_Type_Y(MeshTopology%Num_Side_Sets))
   Allocate(EXO_RuptProperties%SS_BC_Type_Z(MeshTopology%Num_Side_Sets))
   Allocate(EXO_RuptProperties%NS_BC_Type_X(MeshTopology%Num_Node_Sets))
   Allocate(EXO_RuptProperties%NS_BC_Type_Y(MeshTopology%Num_Node_Sets))
   Allocate(EXO_RuptProperties%NS_BC_Type_Z(MeshTopology%Num_Node_Sets))

   EXO_RuptProperties%Is_Brittle     = PETSC_TRUE
   EXO_RuptProperties%Is_Domain      = PETSC_TRUE
   EXO_RuptProperties%Has_BodyForce  = PETSC_FALSE
   EXO_RuptProperties%EB_BC_Type_X   = -27
   EXO_RuptProperties%EB_BC_Type_Y   = MeshTopology%Elem_Blk(:)%ID
   EXO_RuptProperties%EB_BC_Type_Z   = BC_TYPE_DIRI
   EXO_RuptProperties%SS_BC_Type_X   = BC_TYPE_DIRI
   EXO_RuptProperties%SS_BC_Type_Y   = BC_TYPE_DIRI
   EXO_RuptProperties%SS_BC_Type_Z   = BC_TYPE_DIRI
   EXO_RuptProperties%NS_BC_Type_X   = MeshTopology%Node_Set(:)%ID
   EXO_RuptProperties%NS_BC_Type_Y   = BC_TYPE_DIRI
   EXO_RuptProperties%NS_BC_Type_Z   = BC_TYPE_DIRI
   
   
   Call  EXO_RuptPropertiesWrite(MyEXO, MeshTopology, EXO_RuptProperties)
!   Call MeshDestroy(MeshTopology%Mesh, iErr); CHKERRQ(ierr)
   
   Call MEF90_Finalize()

End Program  TestRuptStruct
