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
   Type(Mesh)                                   :: Tmp_Mesh
   Type (EXO_Type)                              :: EXO, MyEXO
   Type (Element2D_Scal), Dimension(:), Pointer :: Elem2DA
   
   PetscTruth                                   :: HasPrefix
   PetscTruth                                   :: verbose
   PetscErrorCode                               :: iErr, iBlk, iE, i
   Character(len=256)                           :: CharBuffer, IOBuffer, filename
   Character(len=256)                           :: prefix
   Type(SectionReal)                            :: rSec
   Type(SectionInt)                             :: iSec

   PetscInt                                     :: Point
   PetscInt                                     :: SecSize

   PetscReal                                    :: rVal
   PetscReal, Dimension(:), Pointer             :: rVals
   PetscInt                                     :: iVal
   PetscInt, Dimension(:), Pointer              :: iVals
   
   PetscInt                                     :: dof = 1

     
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


   If (MEF90_NumProcs == 1) Then
      Call MeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, MeshTopology%mesh, ierr); CHKERRQ(iErr)
   Else
      Call MeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
      Call MeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, MeshTopology%mesh, ierr); CHKERRQ(iErr)
      Call MeshDestroy(Tmp_mesh, ierr); CHKERRQ(iErr)
   End If
   
   Call MeshTopologyReadEXO(MeshTopology, EXO)
   
   Do iBlk = 1, MeshTopology%Num_Elem_Blks
      MeshTopology%Elem_Blk(iBlk)%Elem_Type = MEF90_P1_Lagrange
      Call Init_Elem_Blk_Type(MeshTopology%Elem_Blk(iBlk), MeshTopology%num_dim)
   End Do

   Call MeshGetVertexSectionReal(MeshTopology%Mesh, 'rSec', dof, rSec, iErr); CHKERRQ(iErr)
   Call MeshGetVertexSectionInt (MeshTopology%Mesh, 'iSec', dof, iSec, iErr); CHKERRQ(iErr)
   
   Write(*,*) 'Number of vertices: ', MeshTopology%Num_Verts
   Write(*,*) 'Number of elements: ', MeshTopology%Num_elems


!!! Initializing Sections
   Write(IOBuffer, *) '\n\n === Initializing Sections ===\n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   rVal = 1.00_Kr
   Call SectionRealSet(rSec, rVal, iErr); CHKERRQ(iErr)
   Call SectionIntZero (iSec, iErr); CHKERRQ(iErr)

   Write(IOBuffer, *) 'rSec after SectionRealSet: \n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call SectionRealView(rSec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)

   Write(IOBuffer, *) 'iSec after SectionIntSet:  \n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call SectionIntView(iSec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)

!!! Testing SectionRealUpdate / SectionIntUpdate
   Write(IOBuffer, *) '\n\n === Testing SectionRealUpdate / SectionIntUpdate ===\n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   
   Point = MeshTopology%Num_Elems
   Allocate(rVals(dof))
   rVals = -2.00_Kr
   Call SectionRealUpdate(rSec, Point, rVals, INSERT_VALUES, iErr); CHKERRQ(iErr)

   Point = MeshTopology%Num_Elems+1
   rVals = 2.00_Kr
   Call SectionRealUpdate(rSec, Point, rVals, ADD_VALUES, iErr); CHKERRQ(iErr)

   Write(IOBuffer, *) 'rSec after SectionRealUpdate:  \n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call SectionRealView(rSec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   DeAllocate(rVals)

   Point = MeshTopology%Num_Elems
   Allocate(iVals(dof))
   iVals = -7
   Call SectionIntUpdate(iSec, Point, iVals, INSERT_VALUES, iErr); CHKERRQ(iErr)

   Point = MeshTopology%Num_Elems+1
   iVals = 4
   Call SectionIntUpdate(iSec, Point, iVals, ADD_VALUES, iErr); CHKERRQ(iErr)
   Write(IOBuffer, *) 'iSec after SectionIntUpdate:  \n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call SectionIntView(iSec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   DeAllocate(iVals)
   
!!! Testing SectionRealRestrict / SectionIntRestrict
   Write(IOBuffer, *) '\n\n === SectionRealRestrict / SectionIntRestrict ===\n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

   Point = MeshTopology%Num_Elems
   Call SectionRealRestrict(rSec, Point, rVals, iErr); CHKERRQ(iErr)
   Write(IOBuffer, *) 'rSec(0)=', rVals   
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call SectionRealRestore(rSec, Point, rVals, iErr); CHKERRQ(iErr)

   Point = MeshTopology%Num_Elems+1
   Call SectionRealRestrict(rSec, Point, rVals, iErr); CHKERRQ(iErr)
   Write(IOBuffer, *) ' rSec(1)=', rVals, '\n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call SectionRealRestore(rSec, Point, rVals, iErr); CHKERRQ(iErr)
   

   Point = MeshTopology%Num_Elems
   Call SectionIntRestrict(iSec, Point, iVals, iErr); CHKERRQ(iErr)
   Write(IOBuffer, *) 'iSec(0)=', iVals   
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call SectionIntRestore(iSec, Point, iVals, iErr); CHKERRQ(iErr)

   Point = MeshTopology%Num_Elems+1
   Call SectionIntRestrict(iSec, Point, iVals, iErr); CHKERRQ(iErr)
   Write(IOBuffer, *) ' iSec(1)=', iVals, '\n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call SectionIntRestore(iSec, Point, iVals, iErr); CHKERRQ(iErr)


!!! Testing SectionRealUpdateClosure / SectionIntUpdateClosure
   Write(IOBuffer, *) '\n\n === SectionRealUpdateClosure / SectionIntUpdateClosure ===\n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

   Allocate(rVals(3*dof))
   Point = 2
   rVals = 6.5_Kr
   Call SectionRealUpdateClosure(rSec, MeshTopology%Mesh, Point, rVals, INSERT_VALUES, iErr); CHKERRQ(iErr)

   Point = 18
   rVals = 3.0_Kr
   Call SectionRealUpdateClosure(rSec, MeshTopology%Mesh, Point, rVals, ADD_VALUES, iErr); CHKERRQ(iErr)
   DeAllocate(rVals)

   Write(IOBuffer, *) 'rSec after SectionRealUpdateClosure:  \n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call SectionRealView(rSec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)

   Allocate(iVals(3*dof))
   Point = 2
   iVals = -9
   Call SectionIntUpdateClosure(iSec, MeshTopology%Mesh, Point, iVals, INSERT_VALUES, iErr); CHKERRQ(iErr)

   Point = 18
   iVals = 4
   Call SectionIntUpdateClosure(iSec, MeshTopology%Mesh, Point, iVals, ADD_VALUES, iErr); CHKERRQ(iErr)
   Write(IOBuffer, *) 'iSec after SectionIntUpdate:  \n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call SectionIntView(iSec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   DeAllocate(iVals)

!!! Testing SectionRealRestrictClosure / SectionIntRestrictClosure
   Write(IOBuffer, *) '\n\n === SectionRealRestrictClosure / SectionIntRestrictClosure ===\n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

   Allocate(rVals(3*dof))
   Point = 0
   Call SectionRealRestrictClosure(rSec, MeshTopology%Mesh, Point, 3*dof, rVals, iErr); CHKERRQ(iErr)
   Write(IOBuffer, *) 'rSec(0)=', rVals   
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

   Point = 14
   Call SectionRealRestrictClosure(rSec, MeshTopology%Mesh, Point, 3*dof, rVals, iErr); CHKERRQ(iErr)
   Write(IOBuffer, *) ' rSec(1)=', rVals, '\n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   DeAllocate(rVals)
   

   Allocate(iVals(3*dof))
   Point = 0
   Call SectionIntRestrictClosure(iSec, MeshTopology%Mesh, Point, 3*dof, iVals, iErr); CHKERRQ(iErr)
   Write(IOBuffer, *) 'iSec(0)=', iVals   
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

   Point = 14
   Call SectionIntRestrictClosure(iSec, MeshTopology%Mesh, Point, 3*dof, iVals, iErr); CHKERRQ(iErr)
   Write(IOBuffer, *) ' iSec(1)=', iVals, '\n'
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   DeAllocate(iVals)
   
!!! Testing SectionRealGetSize
   Call SectionRealGetSize(rSec, SecSize, iErr); CHKERRQ(iErr)
   Write(IOBuffer, *) 'size(rSec)=', SecSize, '\n'   
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
   

   Call MeshTopologyDestroy(MeshTopology)
   Call SectionRealDestroy(rSec, iErr); CHKERRQ(iErr)
   Call SectionIntDestroy (iSec, iErr); CHKERRQ(iErr)
   
   
   Call MEF90_Finalize()
End Program TestSectionInt
