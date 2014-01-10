#include "../MEF90/mef90.inc"
Program ThermoElasticity
#include <finclude/petscdef.h>
   Use m_MEF90
   Use m_MEF90_DefMechCtx
   Use m_MEF90_DefMech
   Use petsc
   Implicit NONE   

   PetscErrorCode                                     :: ierr
   Type(MEF90Ctx_Type),target                         :: MEF90Ctx
   Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions
   Type(MEF90CtxGlobalOptions_Type),Parameter         :: MEF90DefaultGlobalOptions = MEF90CtxGlobalOptions_Type( &
                                                         1,                             & ! verbose
                                                         MEF90TimeInterpolation_linear, & ! timeInterpolation
                                                         0.0_Kr,                        & ! timeMin
                                                         1.0_Kr,                        & ! timeMax
                                                         11,                            & ! timeNumStep
                                                         MEF90FileFormat_EXOSingle)       ! fileFormat

   Type(MEF90DefMechCtx_Type)                         :: MEF90DefMechCtx
   Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
   Type(MEF90DefMechGlobalOptions_Type),Parameter     :: MEF90DefMechDefaultGlobalOptions2D = MEF90DefMechGlobalOptions_Type( &
                                                         MEF90DefMech_ModeQuasiStatic, & ! mode
                                                         PETSC_TRUE,          & ! disp_addNullSpace
                                                         3,                   & ! DisplacementOffset
                                                         2,                   & ! DamageOffset
                                                         0,                   & ! boundaryDisplacementOffset
                                                         0,                   & ! boundaryDamageOffset
                                                         0,                   & ! temperatureOffset
                                                         4,                   & ! ForceOffset
                                                         3,                   & ! pressureForceOffset
                                                         0,                   & ! plasticStrainOffset
                                                         0,                   & ! StressOffset
                                                         MEF90Scaling_Linear, & ! boundaryDisplacementScaling
                                                         MEF90Scaling_Linear, & ! ForceScaling
                                                         MEF90Scaling_Linear)   ! pressureForceScaling
   Type(MEF90DefMechGlobalOptions_Type),Parameter     :: MEF90DefMechDefaultGlobalOptions3D = MEF90DefMechGlobalOptions_Type( &
                                                         MEF90DefMech_ModeQuasiStatic, & ! mode
                                                         PETSC_TRUE,          & ! disp_addNullSpace
                                                         3,                   & ! DisplacementOffset
                                                         2,                   & ! DamageOffset
                                                         0,                   & ! boundaryDisplacementOffset
                                                         0,                   & ! boundaryDamageOffset
                                                         0,                   & ! temperatureOffset
                                                         4,                   & ! ForceOffset
                                                         3,                   & ! pressureForceOffset
                                                         0,                   & ! plasticStrainOffset
                                                         0,                   & ! StressOffset
                                                         MEF90Scaling_Linear, & ! boundaryDisplacementScaling
                                                         MEF90Scaling_Linear, & ! ForceScaling
                                                         MEF90Scaling_Linear)   ! pressureForceScaling

   Type(MEF90DefMechCellSetOptions_Type),Parameter    :: MEF90DefMechDefaultCellSetOptions = MEF90DefMechCellSetOptions_Type( &
                                                         -1,                                      & ! elemTypeShortIDDispl will be overriden
                                                         -1,                                      & ! elemTypeShortIDDamage will be overriden
                                                         (/0.0_Kr,0.0_Kr,0.0_Kr/),                & ! force
                                                         0.0_Kr,                                  & ! pressureForce
                                                         MEF90DefMech_defectLawElasticity,        & ! defect law
                                                         (/PETSC_FALSE,PETSC_FALSE,PETSC_FALSE/), & ! Has Displacement BC
                                                         0.0_Kr,                                  & ! boundary Displacement
                                                         PETSC_FALSE,                             & ! Has Damage BC
                                                         0.0_Kr)                                    ! Boundary Damage
   Type(MEF90DefMechVertexSetOptions_Type),Parameter  :: MEF90DefMechDefaultVertexSetOptions = MEF90DefMechVertexSetOptions_Type( &
                                                         (/PETSC_FALSE,PETSC_FALSE,PETSC_FALSE/), & ! Has Displacement BC
                                                         0.0_Kr,                                  & ! boundary Displacement
                                                         PETSC_FALSE,                             & ! Has Damage BC
                                                         0.0_Kr)                                    ! boundary Damage
   PetscBag,dimension(:),pointer                      :: MEF90MatPropBag                                                      

   Type(DM),target                                    :: Mesh
   Type(IS)                                           :: setIS,CellSetGlobalIS
   PetscInt,Dimension(:),Pointer                      :: setID
   PetscInt                                           :: numset,set
   PetscReal,Dimension(:),Pointer                     :: time,energy,work

   Type(SNES)                                         :: snesDisp
   SNESConvergedReason                                :: snesDispConvergedReason
   Type(Vec)                                          :: residualDisp
          
   PetscBool                                          :: flg
   Character(len=MEF90_MXSTRLEN)                      :: IOBuffer
   Integer                                            :: numfield
   
   Integer                                            :: step
   PetscInt                                           :: dim
      
   !!! Initialize MEF90
   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   Call MEF90Initialize(ierr)

   !!! Get all MEF90-wide options
   Call MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90DefaultGlobalOptions,ierr);CHKERRQ(ierr)
   Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)

   !!! Get DM from mesh
   Call MEF90CtxGetDMMeshEXO(MEF90Ctx,Mesh,ierr);CHKERRQ(ierr)
   Call DMMeshGetDimension(Mesh,dim,ierr);CHKERRQ(ierr)
   Call DMMeshSetMaxDof(Mesh,dim,ierr);CHKERRQ(ierr) 
   Call DMSetBlockSize(Mesh,dim,ierr);CHKERRQ(ierr)
   
   !!! Open output file
   Call MEF90CtxOpenEXO(MEF90Ctx,Mesh,ierr)
   
   !!! Create DefMech context, get all DefMech options
   Call MEF90DefMechCtxCreate(MEF90DefMechCtx,Mesh,MEF90Ctx,ierr);CHKERRQ(ierr)
   If (dim == 2) Then
      Call MEF90DefMechCtxSetFromOptions(MEF90DefMechCtx,PETSC_NULL_CHARACTER,MEF90DefMechDefaultGlobalOptions2D, &
                                         MEF90DefMechDefaultCellSetOptions,MEF90DefMechDefaultVertexSetOptions,ierr)
   Else
      Call MEF90DefMechCtxSetFromOptions(MEF90DefMechCtx,PETSC_NULL_CHARACTER,MEF90DefMechDefaultGlobalOptions3D, &
                                         MEF90DefMechDefaultCellSetOptions,MEF90DefMechDefaultVertexSetOptions,ierr)
   End If
   Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)
   
   !!! Get material properties bags
   If (dim == 2) Then
      Call MEF90MatPropBagSetFromOptions(MEF90MatPropBag,MEF90DefMechCtx%DM,MEF90_Mathium2D,MEF90Ctx,ierr)
   Else
      Call MEF90MatPropBagSetFromOptions(MEF90DefMechCtx%MaterialPropertiesBag,MEF90DefMechCtx%DM,MEF90_Mathium3D,MEF90Ctx,ierr)
   End If   
   MEF90DefMechCtx%MaterialPropertiesBag => MEF90MatPropBag
   Call MEF90CtxGetTime(MEF90Ctx,time,ierr)

   !!! Set the data layout
   Call MEF90DefMechCtxSetSections(MEF90DefMechCtx,ierr)

   !!! Create vectors
   Call MEF90DefMechCtxCreateVectors(MEF90DefMechCtx,ierr)

   !!! 
   !!! Create SNES or TS, Mat and set KSP default options
   !!!
   Call VecDuplicate(MEF90DefMechCtx%displacement,residualDisp,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(residualDisp,"residualDisp",ierr);CHKERRQ(ierr)
   Call MEF90DefMechCreateSolvers(MEF90DefMechCtx,snesDisp,residualDisp,ierr)

   !!! 
   !!! Allocate array of works and energies
   !!!
   Allocate(energy(size(MEF90DefMechCtx%CellSetOptionsBag)))
   energy = 0.0_Kr
   Allocate(work(size(MEF90DefMechCtx%CellSetOptionsBag)))
   work = 0.0_Kr

   !!!
   !!! Try to figure out if the file was formatted
   !!!
   If (MEF90Ctx%rank == 0) Then
      Call EXGVP(MEF90Ctx%fileExoUnit,"N",numfield,ierr)
   End If
   Call MPI_Bcast(numfield,1,MPIU_INTEGER,0,MEF90Ctx%comm,ierr)   
   If (numfield == 0) Then
      Call MEF90DefMechFormatEXO(MEF90DefMechCtx,ierr)
   End If
   
   !!!
   !!! Actual computations / time stepping
   !!!
   If (MEF90DefMechGlobalOptions%mode == MEF90DefMech_ModeQuasiStatic) Then
      Do step = 1,MEF90GlobalOptions%timeNumStep
         Write(IOBuffer,100) step,time(step)
         Call PetscPrintf(MEF90Ctx%comm,IOBuffer,ierr);CHKERRQ(ierr)

         !!! Update fields
         Call MEF90DefMechSetTransients(MEF90DefMechCtx,step,time(step),ierr)
         Call MEF90DefMechUpdateboundaryDisplacement(MEF90DefMechCtx%displacement,MEF90DefMechCtx,ierr)

         !!! Solve SNES
         Call SNESSolve(snesDisp,PETSC_NULL_OBJECT,MEF90DefMechCtx%displacement,ierr);CHKERRQ(ierr)
         Call SNESGetConvergedReason(snesDisp,snesDispConvergedReason,ierr);CHKERRQ(ierr)
         Write(IOBuffer,*) "SNESConvergedReason returned ",snesDispConvergedReason,"\n"
         Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
         
         !!! Compute energies
         energy = 0.0_Kr
         work = 0.0_Kr
         Call MEF90DefMechWork(MEF90DefMechCtx%displacement,MEF90DefMechCtx,work,ierr);CHKERRQ(ierr)
         Call MEF90DefMechElasticEnergy(MEF90DefMechCtx%displacement,MEF90DefMechCtx,energy,ierr);CHKERRQ(ierr)
         Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
         Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
         Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
         Do set = 1, size(setID)
            Write(IOBuffer,101) setID(set),energy(set),work(set),energy(set)-work(set)
            Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
         End Do
         Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
         Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
         Write(IOBuffer,102) sum(energy),sum(work),sum(energy)-sum(work)
         Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
     
         !!! Save results and boundary Values
         Call MEF90DefMechViewEXO(MEF90DefMechCtx,step,ierr)
      End Do
   End If
100 Format("Solving steady state step ",I4,", t=",ES12.5,"\n")
101 Format("cell set ",I4," elastic energy: ",ES12.5," work: ",ES12.5," total: ",ES12.5,"\n")
102 Format("======= Total elastic energy: ",ES12.5," work: ",ES12.5," total: ",ES12.5,"\n")

   !!! Clean up and exit nicely
   If (MEF90DefMechGlobalOptions%mode == MEF90DefMech_ModeQuasiStatic) Then
      Call SNESDestroy(snesDisp,ierr);CHKERRQ(ierr)
      Call VecDestroy(residualDisp,ierr);CHKERRQ(ierr)
   End If

   Call MEF90DefMechCtxDestroyVectors(MEF90DefMechCtx,ierr)

   Call DMDestroy(Mesh,ierr);CHKERRQ(ierr)

   DeAllocate(time)
   DeAllocate(energy)
   DeAllocate(work)
   Call MEF90DefMechCtxDestroy(MEF90DefMechCtx,ierr);CHKERRQ(ierr)
   Call MEF90CtxCloseEXO(MEF90Ctx,ierr)

   Call PetscLogView(MEF90Ctx%logViewer,ierr);CHKERRQ(ierr)
   Call MEF90CtxDestroy(MEF90Ctx,ierr)
   Call MEF90Finalize(ierr)
   Call PetscFinalize(ierr)
End Program ThermoElasticity
