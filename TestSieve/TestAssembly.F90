Program TestAssembly
#include "finclude/petscdef.h"
   Use m_MEF90
   Use petsc
   Implicit NONE   

   Type(DM)                                     :: mesh
   Type(IS)                                     :: cellSetGlobalIS
   Type(EXO_Type)                               :: EXO,MyEXO
   Type(Element2D_Scal),Dimension(:),Pointer    :: Elem2DA
   
   PetscBool                                    :: HasPrefix
   PetscBool                                    :: verbose
   PetscErrorCode                               :: iErr
   Character(len = 256)                         :: filename
   Character(len = 256)                         :: prefix
   Character(len = MEF90_MXSTRLEN)              :: IOBuffer
   Type(SectionReal)                            :: U_Sec
   Type(Vec)                                    :: U_Vec
   Type(VecScatter)                             :: U_Scatter
   
   PetscReal,Dimension(:),Pointer               :: val
   PetscReal,Dimension(:,:),Pointer             :: Kelem
   Type(Mat)                                    :: K
   
   Type(IS)                                     :: setIS,cellIS
   PetscInt,Dimension(:),Pointer                :: setID,cellID
   PetscInt                                     :: set,cell
   
   PetscInt                                     :: conesize,dof
   Type(DM)                                     :: tmpDM
   PetscReal                                    :: one = 1.0_Kr
   Integer                                      :: cpu_ws = 0
   Integer                                      :: io_ws = 0
   PetscReal                                    :: vers
     
   Call MEF90_Initialize()
   dof=1
   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-dof',dof,HasPrefix,iErr)    
   Call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-verbose',verbose,iErr)    
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-p',prefix,HasPrefix,iErr)    
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD,"No input file prefix given\n",iErr)
      Call MEF90_Finalize()
      STOP
   End If

   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'
   If (MEF90_Myrank == 0) Then
      EXO%exoid = EXOPEN(EXO%filename,EXREAD,cpu_ws,io_ws,vers,ierr)
   End If
   If (MEF90_numprocs == 1) Then
      Call DMMeshCreateExodusNG(PETSC_COMM_WORLD,EXO%exoid,mesh,ierr);CHKERRQ(ierr)
   Else
      Call DMMeshCreateExodusNG(PETSC_COMM_WORLD,EXO%exoid,tmpDM,ierr);CHKERRQ(ierr)
      Call DMMeshDistribute(tmpDM,PETSC_NULL_CHARACTER,mesh,ierr);CHKERRQ(iErr)
      Call DMDestroy(tmpDM,ierr);CHKERRQ(iErr)
   End If
   If (MEF90_Myrank == 0) Then
      Call EXCLOS(EXO%exoid,ierr)
   End If
   Call DMMeshSetMaxDof(mesh,dof,iErr); CHKERRQ(iErr) 

   Call DMMeshGetVertexSectionReal(mesh,'U',dof,U_Sec,iErr); CHKERRQ(iErr)
   Call DMMeshCreateGlobalScatter(mesh,U_Sec,U_Scatter,ierr);CHKERRQ(ierr)
   Call DMMeshCreateVector(mesh,U_Sec,U_Vec,ierr);CHKERRQ(ierr)

   Call DMMeshCreateMatrix(mesh,U_Sec,MATMPIAIJ,K,iErr); CHKERRQ(iErr)

   !!! 
   !!! Semi-Synchronized assembly: each cpu goes through all cell sets 
   !!! even if the set does not intersect its subdomain. 
   !!!
   Call PetscPrintf(PETSC_COMM_WORLD,"Doing semi-synchronous assembly\n",ierr);CHKERRQ(ierr)
   Call MatZeroEntries(K,ierr);CHKERRQ(ierr)
   Call SectionRealSet(U_Sec,0.0_Kr,iErr);CHKERRQ(ierr)
   
   Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,setIS)
   Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      !!!
      !!! Get the IS containing the cell ID of the cells in the cell set setID(set)
      !!! and the associated F90 array
      !!!
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),cellIS,ierr); CHKERRQ(iErr)
      Call ISGetIndicesF90(cellIS,cellID,iErr);CHKERRQ(iErr)
      Do cell = 1,size(cellID)
         !!! 
         !!! Get the cone size of the local cell i.e. its number of vertices (+ number of edges in an interpolated mesh)
         !!!
         Call DMMeshGetConeSize(mesh,cellID(cell),conesize,ierr);CHKERRQ(ierr);
            
         Allocate(val(conesize*dof))
         Allocate(Kelem(conesize*dof,conesize*dof))
         val = 1.0_Kr
         Kelem = 1.0_Kr
         Call DMMeshAssembleMatrix(K,mesh,U_Sec,cellID(cell),Kelem,ADD_VALUES,iErr); CHKERRQ(iErr)
         Call SectionRealUpdateClosure(U_Sec,mesh,cellID(cell),val,ADD_VALUES,iErr); CHKERRQ(iErr)  
         DeAllocate(val)
         DeAllocate(Kelem)
      End Do
      !!! The IS must be destroyed and the associated F90 array restored
      !!! There is a bug in ISRestoreIndicesF90 when the F90 array has zero size, 
      !!! so we check beforehand
      If (size(cellID) > 0) Then
         Call ISRestoreIndicesF90(cellIS,cellID,iErr);CHKERRQ(iErr)
      End If
      Call ISDestroy(cellIS,ierr);CHKERRQ(ierr)
   End Do
   Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
   Call ISDestroy(setIS,ierr);CHKERRQ(ierr)

   Call MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY,iErr); CHKERRQ(iErr)
   Call MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY,iErr); CHKERRQ(iErr)
   Call SectionRealComplete(U_Sec,ierr);CHKERRQ(ierr)
   Call SectionRealToVec(U_Sec,U_Scatter,SCATTER_FORWARD,U_Vec,ierr);CHKERRQ(iErr)

   
   Call PetscPrintf(PETSC_COMM_WORLD,"Matrix K: \n",iErr); CHKERRQ(iErr)
   Call MatView(K,PETSC_VIEWER_STDOUT_WORLD,iErr); CHKERRQ(iErr)
   Call PetscPrintf(PETSC_COMM_WORLD,"\n\nSection U_Sec: \n",iErr); CHKERRQ(iErr)
   Call SectionRealView(U_Sec,PETSC_VIEWER_STDOUT_WORLD,iErr); CHKERRQ(iErr)
   Call PetscPrintf(PETSC_COMM_WORLD,"\n\nVec U_Vec: \n",iErr); CHKERRQ(iErr)
   Call VecView(U_Vec,PETSC_VIEWER_STDOUT_WORLD,iErr); CHKERRQ(iErr)

   Call MatZeroEntries(K,ierr);CHKERRQ(ierr)
   Call SectionRealSet(U_Sec,0.0,ierr);CHKERRQ(ierr)

   !!! 
   !!! Fully asynchronous assembly: each cpu goes only the cell sets that 
   !!! intersect its subdomain. 
   !!! The only difference from the previous case is that the outer
   !!! loop is given by cellIS obtained with
   !!!    Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)

   Call DMMeshGetLabelIdIS(mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
   Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      !!!
      !!! Get the IS containing the cell ID of the cells in the cell set setID(set)
      !!! and the associated F90 array
      !!!
      Call DMMeshGetStratumIS(mesh,'Cell Sets',setID(set),cellIS,ierr); CHKERRQ(iErr)
      Call ISGetIndicesF90(cellIS,cellID,iErr);CHKERRQ(iErr)
      Do cell = 1,size(cellID)
         !!! 
         !!! Get the cone size of the local cell i.e. its number of vertices (+ number of edges in an interpolated mesh)
         !!! Since we know that in a cell set all cells are of the same type
         !!! We could have allocated just once outside of the loop
         !!!
         Call DMMeshGetConeSize(mesh,cellID(cell),conesize,ierr);CHKERRQ(ierr);
         Allocate(val(conesize*dof))
         Allocate(Kelem(conesize*dof,conesize*dof))
         val = 1.0_Kr
         Kelem = 1.0_Kr
         Call DMMeshAssembleMatrix(K,mesh,U_Sec,cellID(cell),Kelem,ADD_VALUES,iErr); CHKERRQ(iErr)
         Call SectionRealUpdateClosure(U_Sec,mesh,cellID(cell),val,ADD_VALUES,iErr); CHKERRQ(iErr)  
         DeAllocate(val)
         DeAllocate(Kelem)
      End Do
      !!! The IS must be destroyed and the associated F90 array restored
      !!! There is a bug in ISRestoreIndicesF90 when the F90 array has zero size, 
      !!! so we check beforehand
      If (size(cellID) > 0) Then
         Call ISRestoreIndicesF90(cellIS,cellID,iErr);CHKERRQ(iErr)
      End If
      Call ISDestroy(cellIS,ierr);CHKERRQ(ierr)
   End Do
   If (size(setID) > 0) Then
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
   End If
   Call ISDestroy(cellIS,ierr);CHKERRQ(ierr)

   Call MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY,iErr); CHKERRQ(iErr)
   Call MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY,iErr); CHKERRQ(iErr)
   Call SectionRealComplete(U_Sec,ierr);CHKERRQ(ierr)
   Call SectionRealToVec(U_Sec,U_Scatter,SCATTER_FORWARD,U_Vec,ierr);CHKERRQ(iErr)

   
   Call PetscPrintf(PETSC_COMM_WORLD,"Matrix K: \n",iErr); CHKERRQ(iErr)
   Call MatView(K,PETSC_VIEWER_STDOUT_WORLD,iErr); CHKERRQ(iErr)
   Call PetscPrintf(PETSC_COMM_WORLD,"\n\nSection U_Sec: \n",iErr); CHKERRQ(iErr)
   Call SectionRealView(U_Sec,PETSC_VIEWER_STDOUT_WORLD,iErr); CHKERRQ(iErr)
   Call PetscPrintf(PETSC_COMM_WORLD,"\n\nVec U_Vec: \n",iErr); CHKERRQ(iErr)
   Call VecView(U_Vec,PETSC_VIEWER_STDOUT_WORLD,iErr); CHKERRQ(iErr)

   Call DMMeshGetStratumIS(mesh,'height',0,cellIS,ierr);CHKERRQ(ierr)
   Call ISGetIndicesF90(cellIS,cellID,iErr);CHKERRQ(iErr)
   write(*,*) cellID
   Call ISRestoreIndicesF90(cellIS,cellID,iErr);CHKERRQ(iErr)
   Call ISDestroy(cellIS,ierr);CHKERRQ(ierr)
   
   Call MeshTopologyDestroy(MeshTopology)
   Call SectionRealDestroy(U_Sec,iErr); CHKERRQ(iErr)
   Call VecDestroy(U_Vec,iErr);CHKERRQ(ierr)
   Call VecScatterDestroy(U_Scatter,iErr);CHKERRQ(ierr)
   Call MEF90_Finalize()
End Program TestAssembly
