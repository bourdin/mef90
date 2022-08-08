Program TestUpdateSets
#include "finclude/petscdef.h"
   Use m_MEF90
   Use m_MEF90_Materials
   Use petsc
   Implicit NONE   

   PetscErrorCode                      :: ierr
   Type(DM),target                     :: Mesh
   Character(len=MEF90MXSTRLEN)       :: IOBuffer
   PetscInt                            :: dim
   Type(Vec)                           :: VecIn,VecOut
   Type(SectionReal)                   :: Sec
   Type(IS)                            :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof,closureIS
   PetscInt,Dimension(:),Pointer       :: VertexSetGlobalIdx,CellSetGlobalIdx,setIdx,setIdxdof,closureIdx
   PetscInt                            :: set,p,c,nval
   Type(VecScatter)                    :: ScatterSecToVec
   Integer                             :: rank
   PetscReal                           :: val
   PetscReal,Dimension(:),Pointer      :: Ptr
   PetscBool                           :: flg

   Type(MEF90Ctx_Type),target                         :: MEF90Ctx
   Type(MEF90CtxGlobalOptions_Type),Parameter         :: MEF90DefaultGlobalOptions = MEF90CtxGlobalOptions_Type( &
                                                         0,                             & ! verbose
                                                         PETSC_FALSE,                   & ! DryRun 
                                                         MEF90TimeInterpolation_linear, & ! timeInterpolation
                                                         0.0_Kr,                        & ! timeMin
                                                         1.0_Kr,                        & ! timeMax
                                                         11,                            & ! timeNumStep
                                                         MEF90FileFormat_EXOSingle,     & ! fileFormat
                                                         1)                               ! timeFrequency
   Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions

   !!! Initialize MEF90
   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   Call MEF90Initialize(ierr)
   Call MPI_Comm_Rank(PETSC_COMM_WORLD,rank,ierr)

   !!! Get all MEF90-wide options
   Call MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90DefaultGlobalOptions,ierr);CHKERRQ(ierr)
   Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)

   !!! Get DM from mesh
   Call MEF90CtxGetDMMeshEXO(MEF90Ctx,Mesh,ierr);CHKERRQ(ierr)
   dim = 1
   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-dim',dim,flg,ierr);CHKERRQ(ierr);
   Call DMMeshSetMaxDof(Mesh,dim,ierr);CHKERRQ(ierr) 
   
   Call DMMeshGetVertexSectionReal(Mesh,"default",dim,Sec,ierr);CHKERRQ(ierr)
   Call DMMeshSetSectionReal(Mesh,"default",Sec,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(Sec,ierr);CHKERRQ(ierr)

   Call DMCreateGlobalVector(Mesh,VecIn,ierr);CHKERRQ(ierr)
   Call DMCreateGlobalVector(Mesh,VecOut,ierr);CHKERRQ(ierr)
   Call DMMeshGetSectionReal(Mesh,'default',Sec,ierr);CHKERRQ(ierr)
   Call DMMeshCreateGlobalScatter(Mesh,Sec,ScatterSecToVec,ierr);CHKERRQ(ierr)

   !!! Copy From Vertex sets:
   Call DMmeshGetLabelIdIS(Mesh,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
   Call ISGetIndicesF90(VertexSetGlobalIS,VertexSetGlobalIdx,ierr);CHKERRQ(ierr)   
   Write(IOBuffer,*) "VertexSetGlobalIdx: ",VertexSetGlobalIdx,"\n"
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)

   Call VecSet(VecOut,300.0_Kr,ierr);CHKERRQ(ierr)
   Call VecSet(VecIn,100.0_Kr,ierr);CHKERRQ(ierr)
   Do set = 1,size(VertexSetGlobalIdx)
      Call DMMeshGetStratumIS(Mesh,'Vertex Sets',VertexSetGlobalIdx(set),setIS,ierr);CHKERRQ(iErr)
      Call ISGetSize(setIS,nval,ierr);CHKERRQ(ierr)
      Allocate(Ptr(nval))
      Do c = 1, dim
         Call MEF90VecGetValuesISdof(Mesh,VecIn,Ptr,setIS,c,ierr)
         Ptr = Ptr+set*10.0_Kr+c
         Call MEF90VecSetValuesISdof(Mesh,VecIn,Ptr,setIS,c,INSERT_VALUES,ierr)
      End Do !c
      DeAllocate(Ptr)
      Call ISDestroy(SetIS,ierr);CHKERRQ(ierr)
   End Do! set
   Call ISRestoreIndicesF90(VertexSetGlobalIS,VertexSetGlobalIdx,ierr);CHKERRQ(ierr)   
   Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)

   Call VecAssemblyBegin(VecIn,ierr);CHKERRQ(ierr)
   Call VecAssemblyEnd(VecIn,ierr);CHKERRQ(ierr)
   Call VecView(VecIn,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
   Call SectionRealToVec(Sec,ScatterSecToVec,SCATTER_REVERSE,VecIn,ierr);CHKERRQ(ierr) 
   Call SectionRealView(Sec,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)

   !!! Same thing with cell sets
   Call DMmeshGetLabelIdIS(Mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
   Call ISGetIndicesF90(CellSetGlobalIS,CellSetGlobalIdx,ierr);CHKERRQ(ierr)   
   Write(IOBuffer,*) "CellSetGlobalIdx: ",CellSetGlobalIdx,"\n"
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)

   Call VecSet(VecIn,5000.0_Kr,ierr);CHKERRQ(ierr)
   Do set = 1,size(CellSetGlobalIdx)
      Call DMMeshGetStratumIS(Mesh,'Cell Sets',CellSetGlobalIdx(set),setIS,ierr);CHKERRQ(iErr)
      Call MEF90ISCreateCelltoVertex(mesh,PETSC_COMM_WORLD,setIS,closureIS,ierr)
      Call ISGetIndicesF90(closureIS,closureIdx,ierr);CHKERRQ(ierr)
      Write(IOBuffer,*) set,closureIdx,'\n'
      Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Call PetscSynchronizedFlush(PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)
      
      Call ISGetSize(closureIS,nval,ierr);CHKERRQ(ierr)
      Allocate(Ptr(nval))
      Do c = 1, dim
         Call MEF90VecGetValuesISdof(Mesh,VecIn,Ptr,closureIS,c,ierr)
         Ptr = Ptr+set*10.0_Kr+c
         Call MEF90VecSetValuesISdof(Mesh,VecIn,Ptr,closureIS,c,INSERT_VALUES,ierr)
      End Do !c
      DeAllocate(Ptr)

      Call ISRestoreIndicesF90(closureIS,closureIdx,ierr);CHKERRQ(ierr)      
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)   
      Call ISDestroy(closureIS,ierr);CHKERRQ(ierr)   
   End Do
   Call ISRestoreIndicesF90(CellSetGlobalIS,CellSetGlobalIdx,ierr);CHKERRQ(ierr)   
   Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)

   Call VecAssemblyBegin(VecIn,ierr);CHKERRQ(ierr)
   Call VecAssemblyEnd(VecIn,ierr);CHKERRQ(ierr)
   Call VecView(VecIn,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
   Call SectionRealToVec(Sec,ScatterSecToVec,SCATTER_REVERSE,VecIn,ierr);CHKERRQ(ierr) 
   Call SectionRealView(Sec,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)

   Call PetscPrintf(PETSC_COMM_WORLD,"\n\n===== DMMeshISCreateISglobaldof for vertex sets ====\n",ierr)
   Call DMmeshGetLabelIdIS(Mesh,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
   Call ISGetIndicesF90(VertexSetGlobalIS,VertexSetGlobalIdx,ierr);CHKERRQ(ierr)   
   Do set = 1,size(VertexSetGlobalIdx)
      Call DMMeshGetStratumIS(Mesh,'Vertex Sets',VertexSetGlobalIdx(set),setIS,ierr);CHKERRQ(iErr)
      Write(IOBuffer,*) 'set ', set,'\n'
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Call ISView(setIS,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
      Do c = 1, dim
         Call DMMeshISCreateISglobaldof(Mesh,setIS,c-1,setISdof,ierr);CHKERRQ(ierr)
         Write(IOBuffer,*) '   dof ', c,' setISdof\n'
         Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
         Call ISView(setISdof,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
         Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
      End Do
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)   
   End Do
   Call ISRestoreIndicesF90(VertexSetGlobalIS,VertexSetGlobalIdx,ierr);CHKERRQ(ierr)   
   Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)

   Call PetscPrintf(PETSC_COMM_WORLD,"\n\n===== DMMeshISCreateISglobaldof for cell sets ====\n",ierr)
   Call DMmeshGetLabelIdIS(Mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
   Call ISGetIndicesF90(CellSetGlobalIS,CellSetGlobalIdx,ierr);CHKERRQ(ierr)   
   Write(IOBuffer,*) "CellSetGlobalIdx: ",CellSetGlobalIdx,"\n"
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)

   Do set = 1,size(CellSetGlobalIdx)
      Call DMMeshGetStratumIS(Mesh,'Cell Sets',CellSetGlobalIdx(set),setIS,ierr);CHKERRQ(iErr)
      Write(IOBuffer,*) 'set ', set,'\n'
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Call ISView(setIS,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
      Do c = 1, dim
         Call DMMeshISCreateISglobaldof(Mesh,setIS,c-1,setISdof,ierr);CHKERRQ(ierr)
         Write(IOBuffer,*) '   dof ', c,' setISdof\n'
         Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
         Call ISView(setISdof,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
         Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
      End Do
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)   
   End Do
   Call ISRestoreIndicesF90(CellSetGlobalIS,CellSetGlobalIdx,ierr);CHKERRQ(ierr)   
   Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)


   Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(Sec,ierr);CHKERRQ(ierr)
   Call VecDestroy(VecIn,ierr);CHKERRQ(ierr)
   Call VecDestroy(VecOut,ierr);CHKERRQ(ierr)
   Call DMDestroy(Mesh,ierr);CHKERRQ(ierr)
   Call MEF90CtxDestroy(MEF90Ctx,ierr);CHKERRQ(ierr)   
   Call MEF90Finalize(ierr)
   Call PetscFinalize()

End Program TestUpdateSets
