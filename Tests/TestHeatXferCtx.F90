Program  TestHeatXferCtx
#include "petsc/finclude/petsc.h"
   Use m_MEF90
   Use m_MEF90_HeatXfer
   Use m_MEF90_HeatXferDefault
   Use petsc
   Implicit NONE   
   
   PetscErrorCode                                     :: ierr
   Type(MEF90Ctx_Type),target                         :: MEF90Ctx
   Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions
   Type(MEF90HeatXferCtx_Type)                        :: MEF90HeatXferCtx
   Type(MEF90HeatXferGlobalOptions_Type),Pointer      :: MEF90HeatXferGlobalOptions
   Type(tDM)                                          :: dm

   Type(tIS)                                          :: cellSetIS,faceSetIS
   PetscInt,Dimension(:),Pointer                      :: cellSetID,faceSetID
   PetscInt                                           :: step,dim,set,numCellSet,numFaceSet
   PetscReal,Dimension(:),Pointer                     :: time
   Character(len=MEF90MXSTRLEN)                       :: IOBuffer
   PetscBool                                          :: flg
   PetscReal,dimension(:),Pointer                     :: energy,bodyWork,surfaceWork


   !!! Initialize MEF90
   PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
   PetscCallA(MEF90Initialize(ierr))

   !!! Get all MEF90-wide options
   PetscCallA(MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90DefaultGlobalOptions,ierr))
   PetscCallA(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr))

   PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm,MEF90Ctx%geometryFile,PETSC_NULL_CHARACTER,PETSC_TRUE,dm,ierr))
   PetscCallA(DMPlexDistributeSetDefault(dm,PETSC_FALSE,ierr))
   PetscCallA(DMSetUseNatural(dm,PETSC_TRUE,ierr))
   PetscCallA(DMSetFromOptions(dm,ierr))
   PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-heatXfer_dm_view",ierr))

   PetscCallA(MEF90CtxGetTime(MEF90Ctx,time,ierr))

   Inquire(file=MEF90Ctx%resultFile,exist=flg)
   If (flg) Then
      ! we assume that the output file exists and is formatted
      PetscCallA(MEF90CtxOpenEXO(MEF90Ctx,MEF90Ctx%resultViewer,FILE_MODE_APPEND,ierr))
   Else
      ! we need to create the output file
      EXOFormat: block
         PetscInt                                            :: numNodalVar = 1, numCellVar = 1, numFaceVar = 2, numGVar = 0
         Character(len=MEF90MXSTRLEN),Dimension(:),Pointer   :: faceVarName, nodalVarName, cellVarName, gVarName

         Allocate(nodalVarName(numNodalVar))
         Allocate(cellVarName(numCellVar))
         Allocate(faceVarName(numFaceVar))
         Allocate(gVarName(numGVar))
         nodalVarName = ["Temperature        "]
         cellVarName  = ["Flux               "]
         faceVarName  = ["boundaryFlux       ","externalTemperature"]
         PetscCallA(MEF90CtxOpenEXO(MEF90Ctx,MEF90Ctx%resultViewer,FILE_MODE_WRITE,ierr))
         PetscCallA(MEF90EXODMView(dm,MEF90Ctx%resultViewer,MEF90GlobalOptions%elementOrder,ierr))
         PetscCallA(MEF90EXOFormat(MEF90Ctx%resultViewer,gVarName,cellVarName,nodalVarName,faceVarName,time,ierr))
         DeAllocate(nodalVarName)
         DeAllocate(cellVarName)
         DeAllocate(faceVarName)
         DeAllocate(gVarName)
      End block EXOFormat
   End If
   distribute: Block 
       Type(tDM),target                    :: dmDist
       PetscInt                            :: ovlp = 0
       Type(tPetscSF)                      :: naturalPointSF

       If (MEF90Ctx%NumProcs > 1) Then
           PetscCallA(DMPlexDistribute(dm,ovlp,naturalPointSF,dmDist,ierr))
           PetscCallA(DMPlexSetMigrationSF(dmDist,naturalPointSF, ierr))
           PetscCallA(PetscSFDestroy(naturalPointSF,ierr))
           PetscCallA(DMDestroy(dm,ierr))
           dm = dmDist
       End If
   End Block distribute
   PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-heatXfer_dm_view",ierr))

   !!! Create HeatXfer context, get all HeatXfer options
   PetscCallA(MEF90HeatXferCtxCreate(MEF90HeatXferCtx,dm,MEF90Ctx,ierr))
   PetscCallA(MEF90HeatXferCtxSetFromOptions(MEF90HeatXferCtx,PETSC_NULL_CHARACTER,MEF90HeatXferDefaultGlobalOptions,MEF90HeatXferDefaultCellSetOptions,MEF90HeatXferDefaultFaceSetOptions,MEF90HeatXferDefaultVertexSetOptions,ierr))
   !!! We no longer need the DM. We have the megaDM in MEF90HeatXferCtx
   PetscCallA(DMDestroy(dm,ierr))
   PetscCallA(PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr))

   !!! Get parse all materials data from the command line
   PetscCallA(DMGetDimension(MEF90HeatXferCtx%megaDM,dim,ierr))
   If (dim == 2) Then
      PetscCallA(MEF90MatPropBagSetFromOptions(MEF90HeatXferCtx%MaterialPropertiesBag,MEF90HeatXferCtx%megaDM,MEF90Mathium2D,MEF90Ctx,ierr))
   Else
      PetscCallA(MEF90MatPropBagSetFromOptions(MEF90HeatXferCtx%MaterialPropertiesBag,MEF90HeatXferCtx%megaDM,MEF90Mathium3D,MEF90Ctx,ierr))
   End If   

   PetscCallA(DMGetLabelIdIS(MEF90HeatXferCtx%megaDM,MEF90CellSetLabelName,cellSetIS,ierr))
   PetscCallA(MEF90ISAllGatherMerge(PETSC_COMM_WORLD,cellSetIS,ierr))
   PetscCallA(IsGetSize(cellSetIS,numCellSet,ierr))
   PetscCallA(ISGetIndicesF90(cellSetIS,cellSetID,ierr))

   PetscCallA(DMGetLabelIdIS(MEF90HeatXferCtx%megaDM,MEF90FaceSetLabelName,faceSetIS,ierr))
   PetscCallA(MEF90ISAllGatherMerge(PETSC_COMM_WORLD,faceSetIS,ierr))
   PetscCallA(IsGetSize(faceSetIS,numFaceSet,ierr))
   PetscCallA(ISGetIndicesF90(faceSetIS,faceSetID,ierr))

   Allocate(energy(numCellSet))
   Allocate(bodyWork(numCellSet))
   Allocate(surfaceWork(numFaceSet))

   !!! Analysis loop:
   Do step = 1, size(time)
      Write(IOBuffer,'("Step: ",I4," Analysis time: ",ES12.5,"\n")') step,time(step)
      PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
      PetscCallA(VecSet(MEF90HeatXferCtx%temperatureLocal,time(step),ierr))
      PetscCallA(MEF90HeatXferSetTransients(MEF90HeatXferCtx,step,time(step),ierr))
      PetscCallA(MEF90HeatXFerEnergy(MEF90HeatXferCtx,energy,bodyWork,surfaceWork,ierr))
      Do set = 1, numCellSet
         Write(*,'("   Cell set ",I4," energy:       ",ES12.5," body work: ",ES12.5)') cellSetID(set),energy(set),bodyWork(set)
      End Do
      Do set = 1, numFaceSet
         Write(*,'("   Face set ",I4," surface work: ",ES12.5)') faceSetID(set),surfaceWork(set)
      End Do
      PetscCallA(MEF90HeatXferViewEXO(MEF90HeatXferCtx,step,ierr))
   End Do ! step
   DeAllocate(energy)
   DeAllocate(bodyWork)
   DeAllocate(surfaceWork)
   PetscCallA(ISGetIndicesF90(faceSetIS,faceSetID,ierr))
   PetscCallA(ISDestroy(faceSetIS,ierr))
   PetscCallA(ISGetIndicesF90(cellSetIS,cellSetID,ierr))
   PetscCallA(ISDestroy(cellSetIS,ierr))

   PetscCallA(MEF90HeatXferCtxDestroy(MEF90HeatXferCtx,ierr))
   PetscCallA(MEF90CtxDestroy(MEF90Ctx,ierr))
   PetscCallA(MEF90Finalize(ierr))
   PetscCallA(PetscFinalize(ierr))
End Program  TestHeatXferCtx