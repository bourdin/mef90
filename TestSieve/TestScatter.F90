Program TestScatter

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscmeshdef.h"
#include "finclude/petscmatdef.h"

   Use m_MEF90
   Use petsc
   Use petscvec
   Use petscmesh
   Use petscmat

   Implicit NONE   

   Type (MeshTopology_Type)                     :: MeshTopology
   Type (EXO_Type)                              :: EXO, MyEXO
   Type (Element2D_Scal), Dimension(:), Pointer :: Elem2DA
   Type (Vect3D), Dimension(:), Pointer         :: Coords
   
   Type(Mesh)                                   :: Tmp_Mesh
   PetscTruth                                   :: HasPrefix
   PetscTruth                                   :: verbose
   PetscErrorCode                               :: iErr
   PetscInt                                     :: i
   Character(len=256)                           :: CharBuffer, IOBuffer, filename
   Character(len=256)                           :: prefix
   Type(PetscViewer)                            :: viewer, myviewer
   Type(SectionReal)                            :: Sscal, Svect
   Type(Vec)                                    :: Vscal, Vvect
   Type(VecScatter)                             :: ScatterScal, ScatterVect
   Type(Mat)                                    :: KScal, KVect
   PetscInt                                     :: dof = 2
   PetscReal, Dimension(:), Pointer             :: ValScal, ValVect
   PetscInt                                     :: NumDofScal, NumDofVect
     
   Call MEF90_Initialize()
   Call PetscOptionsHasName(PETSC_NULL_CHARACTER, '-verbose', verbose, iErr); CHKERRQ(iErr)
   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-dof', dof,  HasPrefix, iErr); CHKERRQ(iErr) 
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p', prefix, HasPrefix, iErr); CHKERRQ(iErr)    
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No input file prefix given\n", iErr); CHKERRQ(iErr)
      Call MEF90_Finalize()
      STOP
   End If

   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'


   Call MeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
   Call MeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, MeshTopology%mesh, ierr); CHKERRQ(iErr)
   Call MeshDestroy(Tmp_mesh, iErr); CHKERRQ(iErr)
   Call MeshTopologyReadEXO(MeshTopology, EXO)

   Call MeshGetVertexSectionReal(MeshTopology%mesh, 'vect', dof, SVect, ierr); CHKERRQ(iErr)
   Call MeshCreateGlobalScatter(MeshTopology%mesh, SVect, ScatterVect, iErr); CHKERRQ(iErr)
   Call MeshCreateVector(MeshTopology%mesh, SVect, VVect, iErr); CHKERRQ(iErr)

   Call MeshGetVertexSectionReal(MeshTopology%mesh, 'scal', 1, SScal, ierr); CHKERRQ(iErr)
   Call MeshCreateGlobalScatter(MeshTopology%mesh, SScal, ScatterScal, iErr); CHKERRQ(iErr)
   Call MeshCreateVector(MeshTopology%mesh, SScal, VScal, iErr); CHKERRQ(iErr)

   Call MeshSetMaxDof(MeshTopology%Mesh, dof+1, iErr); CHKERRQ(iErr) 

!   MeshTopology%Elem_Blk%Elem_Type    = MEF90_P1_Lagrange
!   Do iBlk = 1, MeshTopology%Num_Elem_Blks
!      Call Init_Elem_Blk_Type(MeshTopology%Elem_Blk(iBlk), MeshTopology%num_dim)
!   End Do

   Call MeshCreateMatrix(MeshTopology%mesh, SVect, MATMPIAIJ, KVect, iErr); CHKERRQ(iErr)
   Allocate(ValVect( (3 * dof)**2))
   ValVect = 1.0_Kr
   Do i = 1, MeshTopology%Num_Elems
      Call assembleMatrix(KVect, MeshTopology%mesh, SVect, i-1, ValVect, ADD_VALUES, iErr);! CHKERRQ(iErr)
   End Do
   Call MatAssemblyBegin(KVect, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
   Call MatAssemblyEnd  (KVect, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
!   Call MatView(KVect, PetscViewer(PETSC_VIEWER_STDOUT_WORLD), iErr); CHKERRQ(iErr)

   Call MeshCreateMatrix(MeshTopology%mesh, SScal, MATMPIAIJ, KScal, iErr); CHKERRQ(iErr)
   Allocate(ValScal(9))
   ValScal = 1.0_Kr
   Do i = 1, MeshTopology%Num_Elems
      Call assembleMatrix(KScal, MeshTopology%mesh, Sscal, i-1, ValScal, ADD_VALUES, iErr); !CHKERRQ(iErr)
   End Do
   Call MatAssemblyBegin(KScal, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
   Call MatAssemblyEnd  (KScal, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
!   Call MatView(KScal, PetscViewer(PETSC_VIEWER_STDOUT_WORLD), iErr); CHKERRQ(iErr)

   
!   Call VecScatterView(ScatterScal, PETSC_NULL_OBJECT, iErr)

   
!   Do i = 1, MeshTopology%Num_Verts
!      Val = i
!      Call MeshUpdateClosure(MeshTopology%mesh, S, MeshTopology%Num_Elems+i-1, Val, iErr); CHKERRQ(iErr)
!   End Do
   Call SectionRealSet(SScal, 2.0_Kr, iErr); CHKERRQ(iErr);
!   Call SectionRealComplete(SScal, iErr); CHKERRQ(iErr)
   
   Write(IOBuffer, *) '\n\nSec Scal: \n'
!   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   Call SectionRealView(SScal, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)

   Call SectionRealCreateLocalVector(SScal, VScal, iErr); CHKERRQ(iErr)
   Call VecView(VScal, PETSC_VIEWER_STDOUT_SELF, iErr)  
   Call VecSet(VScal, -2.0_Kr, iErr); CHKERRQ(iErr)
   Call VecDestroy(VScal, iErr); CHKERRQ(iErr)
   Call SectionRealView(SScal, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
!   Call SectionRealToVec(SScal, ScatterScal, SCATTER_FORWARD, VScal, iErr); CHKERRQ(iErr)

!   Write(IOBuffer, *) '\n\nVec Scal: \n'
!   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
!   Call VecView(VScal, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   
   
   Call MEF90_Finalize()
End Program TestScatter
