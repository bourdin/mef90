Program TestSectionInt

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
   
   PetscTruth                                   :: HasPrefix
   PetscTruth                                   :: verbose
   PetscErrorCode                               :: iErr, iBlk, iE
   Character(len=256)                           :: CharBuffer, IOBuffer, filename
   Character(len=256)                           :: prefix
   Type(PetscViewer)                            :: viewer, myviewer
   Type(SectionReal)                            :: U_Sec
   PetscReal, Dimension(:), Pointer             :: Values
   PetscInt, Dimension(:), Pointer              :: IntValues
   Type(SectionInt)                             :: Flag_Sec
     
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


   Call MeshTopologyReadEXO(MeshTopology, EXO)
   
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

   Allocate(Values(1))
   Call MeshGetVertexSectionInt(MeshTopology%mesh, 1, Flag_Sec, iErr); CHKERRQ(iErr)
   Do iE = 1, MeshTopology%Num_Verts
      Values = iE
!      Call MeshUpdateClosureInt(MeshTopology%Mesh, Flag_Sec, iE-1+MeshTopology%Num_Elems, Values, iErr)
   End Do  
!!!   Call MeshGetCellSectionInt(MeshTopology%mesh, 1, Flag_Sec, iErr); CHKERRQ(iErr)
!!!   Do iE = 1, MeshTopology%Num_Elems
!!!      IntValues = iE
!!!      Call MeshUpdateClosureInt(MeshTopology%Mesh, Flag_Sec, iE-1, IntValues, iErr)
!!!   End Do  
   DeAllocate(Values)
   Call SectionIntView(Flag_Sec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)

   Allocate(Values(1))
   Call MeshGetVertexSectionReal(MeshTopology%mesh, 1, U_Sec, iErr); CHKERRQ(iErr)
   Do iE = 1, MeshTopology%Num_Verts
      Values = iE
      Call MeshUpdateClosure(MeshTopology%Mesh, U_Sec, iE-1+MeshTopology%Num_Elems, Values, iErr)
   End Do  
!!!   Call MeshGetCellSectionReal(MeshTopology%mesh, 1, U_Sec, iErr); CHKERRQ(iErr)
!!!   Do iE = 1, MeshTopology%Num_Elems
!!!      Values = iE
!!!      Call MeshUpdateClosure(MeshTopology%Mesh, U_Sec, iE-1, Values, iErr)
!!!   End Do  
   DeAllocate(Values)
!   Call SectionRealView(U_Sec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)

   Allocate(Values(3))
   Do iE = 1, MeshTopology%Num_Elems
      Call MeshRestrictClosure(MeshTopology%Mesh, U_Sec, iE-1, 3, Values, iErr)
   End Do
   DeAllocate(Values)

   
   Call MeshTopologyDestroy(MeshTopology)
!   Call SectionIntDestroy(Flag_Sec, iErr); CHKERRQ(iErr)
   Call SectionRealDestroy(U_Sec, iErr); CHKERRQ(iErr)
   
   
   If (verbose) Then
      Call PetscViewerFlush(myviewer, iErr); CHKERRQ(iErr)
      Call PetscViewerDestroy(myviewer, iErr); CHKERRQ(iErr)
      Call PetscViewerFlush(viewer, iErr); CHKERRQ(iErr)
      Call PetscViewerDestroy(viewer, iErr); CHKERRQ(iErr)
   End If
   
300 Format(A)   
   Call MEF90_Finalize()
End Program TestSectionInt
