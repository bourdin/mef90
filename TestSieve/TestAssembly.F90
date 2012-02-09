Program TestAssembly

#include "finclude/petscdef.h"

   Use m_MEF90
   Use petsc

   Implicit NONE   

   Type(MeshTopology_Type)                      :: MeshTopology
   Type(DM)                                     :: dmFS
   Type(EXO_Type)                               :: EXO, MyEXO
   Type(Element2D_Scal), Dimension(:), Pointer  :: Elem2DA
   
   PetscBool                                    :: HasPrefix
   PetscBool                                    :: verbose
   PetscErrorCode                               :: iErr
   Character(len=256)                           :: filename
   Character(len=256)                           :: prefix
   Character(len = MEF90_MXSTRLEN)              :: IOBuffer
   Type(SectionReal)                            :: U_Sec
   Type(Vec)                                    :: U_Vec
   Type(VecScatter)                             :: U_Scatter
   
   PetscReal, Dimension(:), Pointer             :: val
   PetscReal, Dimension(:,:), Pointer           :: Kelem
   Type(Mat)                                    :: K
   PetscInt                                     :: cellset,cellIdx,cell
   Type(IS)                                     :: my_setsIS, my_cellsIS
   PetscInt, Dimension(:), Pointer              :: my_sets,my_cells
   PetscInt                                     :: conesize,dof
   PetscInt, Dimension(:), Pointer              :: iRows
   PetscReal                                    :: one = 1.0_Kr
     
   Call MEF90_Initialize()
   dof=1
   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-dof', dof, HasPrefix, iErr)    
   Call PetscOptionsHasName(PETSC_NULL_CHARACTER, '-verbose', verbose, iErr)    
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p', prefix, HasPrefix, iErr)    
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No input file prefix given\n", iErr)
      Call MEF90_Finalize()
      STOP
   End If

   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'

   Call DMMeshCreateExodusNG(PETSC_COMM_WORLD, EXO%filename, MeshTopology%mesh, MeshTopology%meshFS,ierr); CHKERRQ(iErr)


   Call DMView(MeshTopology%mesh,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
   Call DMView(MeshTopology%meshFS,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
   
   Call MeshTopologyGetInfo(MeshTopology, PETSC_COMM_WORLD)
   Call DMMeshSetMaxDof(MeshTopology%Mesh, dof, iErr); CHKERRQ(iErr) 

   Call DMMeshGetVertexSectionReal(MeshTopology%mesh, 'U', dof, U_Sec, iErr); CHKERRQ(iErr)
   Call DMMeshCreateGlobalScatter(MeshTopology%mesh,U_Sec,U_Scatter,ierr);CHKERRQ(ierr)
   Call DMMeshCreateVector(MeshTopology%mesh,U_Sec,U_Vec,ierr);CHKERRQ(ierr)

   Call DMMeshCreateMatrix(MeshTopology%mesh, U_Sec, MATMPIAIJ, K, iErr); CHKERRQ(iErr)

   !!! 
   !!! Semi-Synchronized assembly: each cpu goes through all cell sets 
   !!! even if the set does not intersect its subdomain. 
   !!!
   Call PetscPrintf(PETSC_COMM_WORLD,"Doing semi-synchronous assembly\n",ierr);CHKERRQ(ierr)
   Call MatZeroEntries(K,ierr);CHKERRQ(ierr)
   Call SectionRealSet(U_Sec,0.0_Kr,iErr);CHKERRQ(ierr)
   Do cellSet = 1, MeshTopology%num_elem_blks
      If (MeshTopology%elem_blk(cellSet)%Num_Elems > 0) Then  
         Call ISGetIndicesF90(MeshTopology%elem_blk(cellSet)%Cell_IS,my_cells,iErr);CHKERRQ(iErr)
         Write(IOBuffer,*) MEF90_MyRank, " assembling ", size(my_cells), " cells in cell set ", MeshTopology%elem_blk(cellSet)%ID, "\n"
         Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr)
         Do cellIdx = 1, size(my_cells)
            cell = my_cells(cellIdx)
            Call DMMeshGetConeSize(MeshTopology%mesh,cell,conesize,ierr);CHKERRQ(ierr);
            Allocate(val(conesize*dof))
            Allocate(Kelem(3*conesize,3*conesize))
            val = 1.0_Kr
            Kelem = 1.0_Kr
   
            Call DMMeshAssembleMatrix(K, MeshTopology%mesh, U_Sec, cell, Kelem, ADD_VALUES, iErr); CHKERRQ(iErr)
            Call SectionRealUpdateClosure(U_Sec, MeshTopology%mesh, cell, val, ADD_VALUES, iErr); CHKERRQ(iErr)  
            DeAllocate(val)
            DeAllocate(Kelem)
         End Do
         Call ISRestoreIndicesF90(MeshTopology%elem_blk(cellSet)%Cell_IS,my_cells,iErr);CHKERRQ(iErr)
      Else
         Write(IOBuffer,*) MEF90_MyRank, " skipping cell set ", MeshTopology%elem_blk(cellSet)%ID, "\n"
         Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr)
      End If
   End Do
   Call MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
   Call MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
   Call SectionRealComplete(U_Sec,ierr);CHKERRQ(ierr)
   Call SectionRealToVec(U_Sec,U_Scatter,SCATTER_FORWARD,U_Vec,ierr);CHKERRQ(iErr)
   
   Call PetscPrintf(PETSC_COMM_WORLD, "Matrix K: \n", iErr); CHKERRQ(iErr)
   Call MatView(K, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   Call PetscPrintf(PETSC_COMM_WORLD, "\n\nSection U_Sec: \n", iErr); CHKERRQ(iErr)
   Call SectionRealView(U_Sec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   Call PetscPrintf(PETSC_COMM_WORLD, "\n\nVec U_Vec: \n", iErr); CHKERRQ(iErr)
   Call VecView(U_Vec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)

   !!! 
   !!! Fully asynchronous assembly: each cpu goes only the cell sets that 
   !!! intersect its subdomain. 
   !!! Note that this does not uses MeshTopology at all. The number of dof per cell
   !!! is obtained using DMMeshGetConeSize. It could be cached using the element_type of the current block
   !!!
   Call PetscPrintf(PETSC_COMM_WORLD,"Doing fully asynchronous assembly\n",ierr);CHKERRQ(ierr)
   Call MatZeroEntries(K,ierr);CHKERRQ(ierr)
   Call SectionRealSet(U_Sec,0.0_Kr,iErr);CHKERRQ(ierr)
   
   Call DMMeshGetLabelIdIS(MeshTopology%mesh,'Cell Sets',my_setsIS,iErr);CHKERRQ(iErr)
   Call ISGetIndicesF90(my_setsIS,my_sets,iErr);CHKERRQ(ierr)
   Do cellSet = 1, size(my_sets)
      Call DMMeshGetStratumIS(MeshTopology%mesh,'Cell Sets',my_sets(cellSet),my_cellsIS,ierr); CHKERRQ(iErr)
      Call ISGetIndicesF90(my_cellsIS,my_cells,iErr);CHKERRQ(iErr)
      Write(IOBuffer,*) MEF90_MyRank, " assembling ", size(my_cells), " cells in cell set ", my_sets(cellSet), "\n"
      Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr)
      Do cellIdx = 1, size(my_cells)
         cell = my_cells(cellIdx)
         Call DMMeshGetConeSize(MeshTopology%mesh,cell,conesize,ierr);CHKERRQ(ierr);
         Allocate(val(conesize*dof))
         Allocate(Kelem(3*conesize,3*conesize))
         val = 1.0_Kr
         Kelem = 1.0_Kr

         Call DMMeshAssembleMatrix(K, MeshTopology%mesh, U_Sec, cell, Kelem, ADD_VALUES, iErr); CHKERRQ(iErr)
         Call SectionRealUpdateClosure(U_Sec, MeshTopology%mesh, cell, val, ADD_VALUES, iErr); CHKERRQ(iErr)  
         DeAllocate(val)
         DeAllocate(Kelem)
      End Do
      Call ISRestoreIndicesF90(my_cellsIS,my_cells,iErr);CHKERRQ(iErr)
   End Do
   Call ISGetIndicesF90(my_setsIS,my_sets,iErr);CHKERRQ(ierr)
   
   Call MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
   Call MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
   Call SectionRealComplete(U_Sec,ierr);CHKERRQ(ierr)
   Call SectionRealToVec(U_Sec,U_Scatter,SCATTER_FORWARD,U_Vec,ierr);CHKERRQ(iErr)
   
   Call PetscPrintf(PETSC_COMM_WORLD, "Matrix K: \n", iErr); CHKERRQ(iErr)
   Call MatView(K, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   Call PetscPrintf(PETSC_COMM_WORLD, "\n\nSection U_Sec: \n", iErr); CHKERRQ(iErr)
   Call SectionRealView(U_Sec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   Call PetscPrintf(PETSC_COMM_WORLD, "\n\nVec U_Vec: \n", iErr); CHKERRQ(iErr)
   Call VecView(U_Vec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)

   Call MeshTopologyDestroy(MeshTopology)
   Call SectionRealDestroy(U_Sec, iErr); CHKERRQ(iErr)
   Call VecDestroy(U_Vec,iErr);CHKERRQ(ierr)
   Call VecScatterDestroy(U_Scatter,iErr);CHKERRQ(ierr)
   Call MEF90_Finalize()
End Program TestAssembly
