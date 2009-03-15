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
   PetscErrorCode                               :: iErr, iBlk, iE, i
   Character(len=256)                           :: CharBuffer, IOBuffer, filename
   Character(len=256)                           :: prefix
   Type(PetscViewer)                            :: viewer, myviewer
   Type(SectionReal)                            :: U_Sec
   Type(Vec)                                    :: U_Vec
   PetscReal, Dimension(:), Pointer             :: Values
   PetscInt, Dimension(:), Pointer              :: IntValues
   Type(SectionInt)                             :: Flag_Sec
   Type(SectionInt)                             :: Sec1, Sec2, Distributed_Sec
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


   !!! Create the section using a shortcut. How do I set its name?
!   Call MeshGetVertexSectionInt(Tmp_mesh, 1, Sec1, iErr); CHKERRQ(iErr)
!   Call MeshSetSectionInt(Tmp_mesh, Sec1, iErr); CHKERRQ(iErr)
   !!! The fortran bindings don't exist, but it doesn't really matter as there is currently no way to set the section name from fortran!...
!   Call PetscObjectSetName(Sec1, 'Sec1')
   !!! This would be nice, except that it won't work with the fortran types...

   !!! Do it the most generic way
   Call MeshGetSectionInt(Tmp_Mesh, prefix, Sec2, iErr); CHKERRQ(iErr)
   Do i = 1, num_verts
      Call SectionIntSetFiberDimension(Sec2, i+Num_Elems-1, 1, iErr); CHKERRQ(iErr)
   End Do 
   Call SectionIntAllocate(Sec2, iErr); CHKERRQ(iErr)

   !!! Setup and initialize an internal SectionInt
   Allocate(IntValues(1))
   Do i = 1, num_verts
      IntValues = i
!      Call MeshUpdateClosureInt(Tmp_mesh, Sec1, i+Num_Elems-1, IntValues, iErr); CHKERRQ(iErr)
      Call MeshUpdateClosureInt(Tmp_mesh, Sec2, i+Num_Elems-1, IntValues, iErr); CHKERRQ(iErr)
   End Do 
   DeAllocate(IntValues)
   
!   Call PetscPrintf(PETSC_COMM_WORLD, "Sec1 before MeshDistribute\n"c, iErr); CHKERRQ(iErr)
!   Call SectionIntView(Sec1, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   Call PetscPrintf(PETSC_COMM_WORLD, "Sec2 before MeshDistribute\n"c, iErr); CHKERRQ(iErr)
   Call SectionIntView(Sec2, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   
   
   Call MeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, MeshTopology%mesh, ierr); CHKERRQ(iErr)
   Call MeshDestroy(Tmp_mesh, ierr); CHKERRQ(iErr)
   

!   Call PetscPrintf(PETSC_COMM_WORLD, "Sec1 after MeshDistribute\n"c, iErr); CHKERRQ(iErr)
!   Call SectionIntView(Sec1, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   Call PetscPrintf(PETSC_COMM_WORLD, "Sec2 after MeshDistribute\n"c, iErr); CHKERRQ(iErr)
   Call SectionIntView(Sec2, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   Call MEF90_Finalize()
   STOP


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

   Allocate(IntValues(3))
   Call MeshGetVertexSectionInt(MeshTopology%mesh, 1, Flag_Sec, iErr); CHKERRQ(iErr)
!   Do iE = 1, MeshTopology%Num_Verts
      iE = 5
      IntValues = 1
      Call MeshUpdateClosureInt(MeshTopology%Mesh, Flag_Sec, iE-1, IntValues, iErr)
      iE = 9
      IntValues = 10
      Call MeshUpdateClosureInt(MeshTopology%Mesh, Flag_Sec, iE-1, IntValues, iErr)
!   End Do  
!!!   Call MeshGetCellSectionInt(MeshTopology%mesh, 1, Flag_Sec, iErr); CHKERRQ(iErr)
!!!   Do iE = 1, MeshTopology%Num_Elems
!!!      IntValues = iE
!!!      Call MeshUpdateClosureInt(MeshTopology%Mesh, Flag_Sec, iE-1, IntValues, iErr)
!!!   End Do  
   DeAllocate(IntValues)
   Call SectionIntView(Flag_Sec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)

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
   Call MeshUpdateClosure(MeshTopology%Mesh, U_Sec, iE-1, Values, iErr)
   Values = 10.0
   Call MeshUpdateAddClosure(MeshTopology%Mesh, U_Sec, iE-1, Values, iErr)
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
   
300 Format(A)   
   Call MEF90_Finalize()
End Program TestSectionInt
