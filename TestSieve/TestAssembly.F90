Program TestAssembly

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscmeshdef.h"

   Use m_MEF90
   Use petsc
   Use petscvec
   Use petscmat
   Use petscmesh

   Implicit NONE   

   Type (MeshTopology_Info)                     :: MeshTopology
   Type (EXO_Info)                              :: EXO, MyEXO
   Type(Element2D_Scal), Dimension(:), Pointer  :: Elem2DA
   Type(Vect3D), Dimension(:), Pointer          :: Coords
   
   PetscTruth                                   :: HasPrefix
   PetscTruth                                   :: verbose
   PetscErrorCode                               :: iErr
   PetscInt                                     :: iBlk
   Character(len=256)                           :: CharBuffer, IOBuffer, filename
   Character(len=256)                           :: prefix
   Type(PetscViewer)                            :: viewer, myviewer
   Type(SectionReal)                            :: U_Sec
   Type(Vec)                                    :: U_VecG, U_VecL
   PetscReal, Dimension(:), Pointer             :: val
   PetscReal, DImension(:,:), Pointer           :: Kelem
   Type(Mat)                                    :: K
   Type(VecScatter)                             :: scatter
   PetscInt                                     :: exo_ver
   PetscInt                                     :: i
   PetscReal                                    :: T
   Type(Vect2D), Dimension(:), Pointer          :: V2D
   Type(Vect3D), Dimension(:), Pointer          :: V3D
   PetscReal, Dimension(:), Pointer             :: V_Ptr
     
   Call MEF90_Initialize()
   Call PetscOptionsHasName(PETSC_NULL_CHARACTER, '-verbose', verbose, iErr)    
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p', prefix, HasPrefix, iErr)    
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No input file prefix given\n", iErr)
      Call MEF90_Finalize()
      STOP
   End If

   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'


   Call MeshTopologyReadEXO(MeshTopology, Coords, Elem2DA, EXO)
   
   MeshTopology%Elem_Blk%Elem_Type    = MEF90_P1_Lagrange
   Do iBlk = 1, MeshTopology%Num_Elem_Blks
      Call Init_Elem_Blk_Info(MeshTopology%Elem_Blk(iBlk), MeshTopology%num_dim)
   End Do
   
   MyEXO%comm = PETSC_COMM_SELF
   MyEXO%exoid = EXO%exoid
   Write(MyEXO%filename, 200) trim(prefix), MEF90_MyRank
200 Format(A, '-', I4.4, '.gen')
 
   MyEXO%title = trim(EXO%title)
   MyEXO%Num_QA = EXO%Num_QA
   
   Call Write_MeshTopologyGlobal(MeshTopology, MyEXO, PETSC_COMM_WORLD)

   Call MeshGetVertexSectionReal(MeshTopology%mesh, 1, U_Sec, iErr); CHKERRQ(iErr)
   Call MeshCreateGlobalScatter(MeshTopology%mesh, U_Sec, scatter, iErr); CHKERRQ(iErr)
   Call MeshCreateVector(MeshTopology%mesh, U_Sec, U_VecG, iErr); CHKERRQ(iErr)
   Call MeshCreateMatrix(MeshTopology%mesh, U_Sec, MATMPIAIJ, K, iErr); CHKERRQ(iErr)
   Call SectionRealCreateLocalVector(U_Sec, U_VecL, iErr); CHKERRQ(iErr)

   Call VecSet(U_VecG, 1.0_Kr, iErr); CHKERRQ(iErr)
   Call VecView(U_VecG, PetscViewer(PETSC_VIEWER_STDOUT_WORLD), iErr); CHKERRQ(iErr)
   

   Allocate(val(3))
   Allocate(Kelem(3,3))
   val = 1.0_Kr
   Kelem = 1.0_Kr
   Do i = 1, MeshTopology%Num_Elems
      Call assemblevector(U_VecG, i, val, ADD_VALUES, iErr); CHKERRQ(iErr)
!      Call assembleMatrix(K, i-1, Kelem, ADD_VALUES, iErr); CHKERRQ(iErr)
   End Do
   DeAllocate(val)
   DeAllocate(Kelem)
   Call VecView(U_VecG, PetscViewer(PETSC_VIEWER_STDOUT_WORLD), iErr); CHKERRQ(iErr)

   Call VecDestroy(U_VecG, iErr); CHKERRQ(iErr)
   Call VecDestroy(U_VecL, iErr); CHKERRQ(iErr)

   Call MeshTopologyDestroy(MeshTopology)
   Call SectionRealDestroy(U_Sec, iErr); CHKERRQ(iErr)
   DeAllocate(V2D)
   
   
   If (verbose) Then
      Call PetscViewerFlush(myviewer, iErr); CHKERRQ(iErr)
      Call PetscViewerDestroy(myviewer, iErr); CHKERRQ(iErr)
      Call PetscViewerFlush(viewer, iErr); CHKERRQ(iErr)
      Call PetscViewerDestroy(viewer, iErr); CHKERRQ(iErr)
      Call VecScatterDestroy(scatter, iErr); CHKERRQ(iErr)
   End If
   
300 Format(A)   
   Call MEF90_Finalize()
End Program TestAssembly
