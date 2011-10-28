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

   Type (MeshTopology_Type)                     :: MeshTopology
   Type (EXO_Type)                              :: EXO, MyEXO
   Type(Element2D_Scal), Dimension(:), Pointer  :: Elem2DA
   
   PetscBool                                   :: HasPrefix
   PetscBool                                   :: verbose
   PetscErrorCode                               :: iErr, iBlk, iE, i
   Character(len=256)                           :: CharBuffer, IOBuffer, filename
   Character(len=256)                           :: prefix
   Type(PetscViewer)                            :: viewer, myviewer
   Type(SectionReal)                            :: U_Sec
   Type(Vec)                                    :: U_Vec
   PetscReal, Dimension(:), Pointer             :: Values
   PetscInt, Dimension(:), Pointer              :: IntValues
   Type(SectionInt)                             :: Flag_Sec
   Type(SectionInt)                             :: Sec1, Sec2
   Type(Mesh)                                   :: Tmp_Mesh
   PetscInt                                     :: Junk, Num_Verts, num_Elems
     
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


   Call MeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
   Call MeshExodusGetInfo(Tmp_mesh, Junk, Num_Verts, Num_Elems, Junk, Junk, iErr); CHKERRQ(iErr)


   !!! Create the Section on cpu0 
   Call MeshGetSectionInt(Tmp_Mesh, prefix, Sec1, iErr); CHKERRQ(iErr)
   Do i = 1, num_verts
      Call SectionIntSetFiberDimension(Sec1, i+Num_Elems-1, 1, iErr); CHKERRQ(iErr)
   End Do 
   Call SectionIntAllocate(Sec1, iErr); CHKERRQ(iErr)

   !!! Setup and initialize an internal SectionInt
   Allocate(IntValues(1))
   Do i = 1, num_verts
      IntValues = i
      Call SectionIntUpdateClosure(Sec1, Tmp_mesh, i+Num_Elems-1, IntValues, INSERT_VALUES, iErr); CHKERRQ(iErr)
   End Do 
   DeAllocate(IntValues)
   
   Call PetscPrintf(PETSC_COMM_WORLD, "Sec1 before MeshDistribute\n"c, iErr); CHKERRQ(iErr)
   Call SectionIntView(Sec1, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   
  
   Call MeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, MeshTopology%mesh, ierr); CHKERRQ(iErr)
   Call MeshDestroy(Tmp_mesh, ierr); CHKERRQ(iErr)
   

   !!! Create Sec2
   Call MeshGetSectionInt(MeshTopology%mesh, prefix, Sec2, iErr); CHKERRQ(iErr)
!!!   Do i = 1, MeshTopology%Num_Verts
!!!      Call SectionIntSetFiberDimension(Sec2, i+MeshTopology%Num_Elems-1, 1, iErr); CHKERRQ(iErr)
!!!   End Do 
!!!   Call SectionIntAllocate(Sec2, iErr); CHKERRQ(iErr)
!!!
!!!   Call SectionIntDistribute(Sec1, MeshTopology%mesh, Sec2, iErr); CHKERRQ(iErr)

   Call PetscPrintf(PETSC_COMM_WORLD, "Sec2 after SectionIntDistribute\n"c, iErr); CHKERRQ(iErr)
   Call SectionIntView(Sec2, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)


   Call MeshTopologyGetInfo(MeshTopology, PETSC_COMM_WORLD)
   
   MeshTopology%Elem_Blk%Elem_Type    = MEF90_P1_Lagrange
   Do iBlk = 1, MeshTopology%Num_Elem_Blks
      Call Init_Elem_Blk_Type(MeshTopology%Elem_Blk(iBlk), MeshTopology%num_dim)
   End Do
   
   MyEXO%comm = PETSC_COMM_SELF
   MyEXO%exoid = EXO%exoid
   Write(MyEXO%filename, 200) trim(prefix), MEF90_MyRank
200 Format(A, '-', I4.4, '.gen')
 
   MyEXO%title = trim(EXO%title)
   MyEXO%Num_QA = EXO%Num_QA
   
   Call Write_MeshTopologyGlobal(MeshTopology, MyEXO, PETSC_COMM_WORLD)

   Allocate(IntValues(3))

!   Call MeshGetVertexSectionInt(MeshTopology%mesh, 1, Flag_Sec, iErr); CHKERRQ(iErr)
   iE = 1
   IntValues = 1
!   Call SectionIntUpdateClosure(Flag_Sec, MeshTopology%Mesh, iE-1, IntValues, ADD_VALUES, iErr)
   IntValues = 34

   Call MEF90_Finalize()
   STOP

   Call SectionIntUpdateClosure(Flag_Sec, MeshTopology%Mesh, iE-1, IntValues, ADD_VALUES, iErr)
   iE = 2
   IntValues = 10
!      Call MeshUpdateClosureInt(MeshTopology%Mesh, Flag_Sec, iE-1, IntValues, iErr)
!   End Do  
!!!   Call MeshGetCellSectionInt(MeshTopology%mesh, 1, Flag_Sec, iErr); CHKERRQ(iErr)
!!!   Do iE = 1, MeshTopology%Num_Elems
!!!      IntValues = iE
!!!      Call MeshUpdateClosureInt(MeshTopology%Mesh, Flag_Sec, iE-1, IntValues, iErr)
!!!   End Do  
!   DeAllocate(IntValues)
   Call SectionIntView(Flag_Sec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   
   
   Call MEF90_Finalize()
   STOP
   
   
   Call MeshGetVertexSectionReal(MeshTopology%mesh, 1, U_Sec, iErr); CHKERRQ(iErr)

!   Allocate(Values(1))
!   Do iE = 1, MeshTopology%Num_Verts
!      Values = 1
!      Call MeshUpdateClosure(MeshTopology%Mesh, U_Sec, iE-1+MeshTopology%Num_Elems, Values, iErr)
!      Call SectionRealComplete(U_Sec, iErr)
!   End Do  
!!!   Call MeshGetCellSectionReal(MeshTopology%mesh, 1, U_Sec, iErr); CHKERRQ(iErr)
!!!   Do iE = 1, MeshTopology%Num_Elems
!!!      Values = iE
!!!      Call MeshUpdateClosure(MeshTopology%Mesh, U_Sec, iE-1, Values, iErr)
!!!   End Do  
!   DeAllocate(Values)

   Allocate(Values(3))
   iE = 5
   Values = 1.0   
   Call SectionRealUpdateClosure(U_Sec, MeshTopology%Mesh, iE-1, Values, INSERT_VALUES, iErr)
   Values = 10.0
   Call SectionRealUpdateClosure(U_Sec, MeshTopology%Mesh, iE-1, Values, ADD_VALUES, iErr)
   DeAllocate(Values)
   Call SectionRealView(U_Sec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   
!   Allocate(Values(3))
!   Do iE = 1, MeshTopology%Num_Elems
!      Call MeshRestrictClosure(MeshTopology%Mesh, U_Sec, iE-1, 3, Values, iErr)
!   End Do
!   DeAllocate(Values)

   
   Call MeshTopologyDestroy(MeshTopology)
!   Call SectionIntDestroy(Flag_Sec, iErr); CHKERRQ(iErr)
   Call SectionRealDestroy(U_Sec, iErr); CHKERRQ(iErr)
   
   
   If (verbose) Then
      Call PetscViewerFlush(myviewer, iErr); CHKERRQ(iErr)
      Call PetscViewerDestroy(myviewer, iErr); CHKERRQ(iErr)
      Call PetscViewerFlush(viewer, iErr); CHKERRQ(iErr)
      Call PetscViewerDestroy(viewer, iErr); CHKERRQ(iErr)
   End If

   Write(*,*) 'CALLING FINALIZE'   
300 Format(A)   
   Call MEF90_Finalize()
End Program TestSectionInt
