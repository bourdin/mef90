Program TestNSP
#include "finclude/petscdef.h"
   Use m_MEF90
   Use m_MEF_Materials
   Use petsc
   Implicit NONE   

   PetscErrorCode                      :: ierr
   Type(DM),target                     :: Mesh,MeshClone
   Character(len=MEF90_MXSTRLEN)       :: IOBuffer
   PetscInt                            :: dim,bs
   Type(Mat)                           :: matDisp
   Type(MatNullSpace)                  :: nspDisp
   Type(Vec)                           :: coordVec,tmpVec
   PetscReal,Dimension(:,:),Pointer    :: coordPtr
   PetscReal,Dimension(:),Pointer      :: coordPCPtr
   Type(SectionReal)                   :: defaultSection,coordSec
   Type(VecScatter)                    :: ScatterSecToVec
   PetscBool                           :: flg
   Integer                             :: rank,numproc,i,r,c
   PetscInt                            :: localSize,globalSize


   Type(MEF90Ctx_Type),target                         :: MEF90Ctx
   Type(MEF90CtxGlobalOptions_Type),Parameter         :: MEF90DefaultGlobalOptions = MEF90CtxGlobalOptions_Type( &
                                                         0,                             & ! verbose
                                                         MEF90TimeInterpolation_linear, & ! timeInterpolation
                                                         0.0_Kr,                        & ! timeMin
                                                         1.0_Kr,                        & ! timeMax
                                                         11,                            & ! timeNumStep
                                                         MEF90FileFormat_EXOSingle)       ! fileFormat
   Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions

   !!! Initialize MEF90
   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   Call MEF90_Initialize(ierr)
   Call MPI_Comm_Rank(PETSC_COMM_WORLD,rank,ierr)
   Call MPI_Comm_Size(PETSC_COMM_WORLD,numproc,ierr)

   !!! Get all MEF90-wide options
   Call MEF90Ctx_Create(PETSC_COMM_WORLD,MEF90Ctx,MEF90DefaultGlobalOptions,ierr);CHKERRQ(ierr)
   Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)

   !!! Get DM from mesh
   Call MEF90Ctx_GetDMMeshEXO(MEF90Ctx,Mesh,ierr);CHKERRQ(ierr)
   Call DMMeshGetDimension(Mesh,dim,ierr);CHKERRQ(ierr)
   Call DMSetBlockSize(Mesh,dim,ierr);CHKERRQ(ierr)
   Call DMMeshSetMaxDof(Mesh,dim,ierr);CHKERRQ(ierr) 
   Call DMMeshGetVertexSectionReal(Mesh,"default",dim,defaultSection,ierr);CHKERRQ(ierr)
   Call DMMeshSetSectionReal(Mesh,"default",defaultSection,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(defaultSection,ierr);CHKERRQ(ierr)

   Call DMCreateMatrix(Mesh,MATAIJ,matDisp,iErr);CHKERRQ(ierr)

   Call DMMeshGetCoordinatesF90(mesh,coordPtr,ierr);CHKERRQ(ierr)
   Do r = 0, numproc-1
      If (r == rank) Then
         Write(IOBuffer,100) rank,"coordPtr \n"
         Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr)
         Do i = 1, size(CoordPtr,1)
            If (dim ==2) Then
               Write(IOBuffer,102) rank,i, CoordPtr(i,:)
            Else
               Write(IOBuffer,103) rank,i, CoordPtr(i,:)
            End If
            Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr)
         End Do
      End If
      Call MPI_Barrier(PETSC_COMM_WORLD,ierr)
   End Do
100 Format("[",I3.3,"]: ",A,"\n")
102 Format("[",I3.3,"]: ",I3.3,"  ", 2(F5.2,"  "),"\n")
103 Format("[",I3.3,"]: ",I3.3,"  ", 3(F5.2,"  "),"\n")
   Allocate(coordPCPtr(size(CoordPtr)))
   coordPCPtr = reshape(transpose(coordPtr),[size(CoordPtr)])
   Call DMMeshRestoreCoordinatesF90(mesh,coordPtr,ierr);



   Call DMMeshGetSectionReal(Mesh,'coordinates',coordSec,ierr);CHKERRQ(ierr)
   Call DMMeshCreateGlobalScatter(Mesh,coordSec,ScatterSecToVec,ierr);CHKERRQ(ierr)
   Call DMCreateGlobalVector(Mesh,coordVec,ierr)
   Call SectionRealToVec(coordSec,ScatterSecToVec,SCATTER_FORWARD,coordVec,ierr);CHKERRQ(ierr)
   !Call SectionRealView(coordSec,PETSC_VIEWER_STDOUT_WORLD,ierr)
   !Call VecView(coordVec,PETSC_VIEWER_STDOUT_WORLD,ierr)

   Call MatNullSpaceCreateRigidBody(coordVec,nspDisp,ierr);CHKERRQ(ierr)
   Call MatNullSpaceView(nspDisp,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
   Call MatSetNearNullSpace(matDisp,nspDisp,ierr);CHKERRQ(ierr)
   Call MatNullSpaceDestroy(nspDisp,ierr);CHKERRQ(ierr)
   
   Call DMMeshClone(Mesh,MeshClone,ierr);CHKERRQ(ierr)
   Call DMSetBlockSize(MeshClone,1,ierr);CHKERRQ(ierr)
   Call DMMeshSetMaxDof(MeshClone,1,ierr);CHKERRQ(ierr) 
   Call DMMeshGetVertexSectionReal(MeshClone,"default",dim,defaultSection,ierr);CHKERRQ(ierr)
   Call DMMeshSetSectionReal(MeshClone,"default",defaultSection,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(defaultSection,ierr);CHKERRQ(ierr)

   Call DMCreateGlobalVector(MeshClone,coordVec,ierr)
   Call VecGetBlockSize(coordVec,bs,ierr);CHKERRQ(ierr)
   Write(IOBuffer,*) 'Block size for Vec created from MeshClone: ',bs,'\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
   Call VecDestroy(coordVec,ierr);CHKERRQ(ierr)
   Call DMCreateGlobalVector(Mesh,coordVec,ierr)
   Call VecGetBlockSize(coordVec,bs,ierr);CHKERRQ(ierr)
   Write(IOBuffer,*) 'Block size for Vec created from Mesh: ',bs,'\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
   Call VecDestroy(coordVec,ierr);CHKERRQ(ierr)
   
   
   
   Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
   Call VecDestroy(coordVec,ierr);CHKERRQ(ierr)
   Call DMDestroy(Mesh,ierr);CHKERRQ(ierr)
   Call MEF90Ctx_Destroy(MEF90Ctx,ierr);CHKERRQ(ierr)   
   Call MEF90_Finalize(ierr)
   Call PetscFinalize()
End Program TestNSP
