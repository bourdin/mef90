#include "../MEF90/mef90.inc"
Program ThermoElasticity
#include <finclude/petscdef.h>
   Use m_MEF90
   Use m_MEF90_DefMechCtx
   Use m_MEF90_DefMech
   Use petsc
   Implicit NONE   
   PetscErrorCode                                     :: ierr
   Type(MEF90DefMechCtx_Type)                         :: MEF90DefMechCtx
   Type(MEF90DefMechGlobalOptions_Type),Parameter     :: MEF90DefMechDefaultGlobalOptions2D = MEF90DefMechGlobalOptions_Type( &
                                                         MEF90DefMech_ModeQuasiStatic, & ! mode
                                                         PETSC_TRUE,          & ! disp_addNullSpace
                                                         1,                   & ! DisplacementOffset
                                                         3,                   & ! DamageOffset
                                                         4,                   & ! boundaryDisplacementOffset
                                                         5,                   & ! boundaryDamageOffset
                                                         6,                   & ! inelasticStrainOffset
                                                         1,                   & ! ForceOffset
                                                         3,                   & ! pressureForceOffset
                                                         4,                   & ! inelasticStrainCellOffset
                                                         8,                   & ! StressOffset
                                                         MEF90Scaling_Linear, & ! boundaryDisplacementScaling
                                                         MEF90Scaling_Linear, & ! ForceScaling
                                                         MEF90Scaling_Linear)   ! pressureForceScaling
   Type(MEF90DefMechGlobalOptions_Type),Parameter     :: MEF90DefMechDefaultGlobalOptions3D = MEF90DefMechGlobalOptions_Type( &
                                                         MEF90DefMech_ModeQuasiStatic, & ! mode
                                                         PETSC_TRUE,          & ! disp_addNullSpace
                                                         1,                   & ! DisplacementOffset
                                                         4,                   & ! DamageOffset
                                                         5,                   & ! boundaryDisplacementOffset
                                                         8,                   & ! boundaryDamageOffset
                                                         9,                   & ! inelasticStrainOffset
                                                         1,                   & ! ForceOffset
                                                         3,                   & ! pressureForceOffset
                                                         4,                   & ! inelasticStrainCellOffset
                                                         11,                  & ! StressOffset
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
   Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
                                                         
   Type(MEF90Ctx_Type),target                         :: MEF90Ctx
   Type(MEF90CtxGlobalOptions_Type),Parameter         :: MEF90DefaultGlobalOptions = MEF90CtxGlobalOptions_Type( &
                                                         1,                             & ! verbose
                                                         MEF90TimeInterpolation_linear, & ! timeInterpolation
                                                         0.0_Kr,                        & ! timeMin
                                                         1.0_Kr,                        & ! timeMax
                                                         11,                            & ! timeNumStep
                                                         MEF90FileFormat_EXOSingle)       ! fileFormat
   Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions
   PetscBag,dimension(:),pointer                      :: MEF90MatPropBag

   Type(DM),target                                    :: Mesh
   Type(IS)                                           :: setIS,cellIS,CellSetGlobalIS
   PetscInt,Dimension(:),Pointer                      :: setID
   PetscInt                                           :: numset,set
   Type(SectionReal)                                  :: defaultSection,coordSec
   Type(Vec)                                          :: residualDisplacement
   Type(Vec)                                          :: coordVec
   PetscReal,Dimension(:),Pointer                     :: time,energy,work
   Type(VecScatter)                                   :: ScatterSecToVec
   PetscReal,Dimension(:,:),Pointer                   :: coordPtr
   PetscReal,Dimension(:),Pointer                     :: coordPCPtr

   Type(SNES)                                         :: snesDisp
   SNESConvergedReason                                :: snesDispConvergedReason
   Type(KSP)                                          :: kspDisp
   Type(PC)                                           :: pcDisp
   Type(Mat)                                          :: matDisp
   Type(MatNullSpace)                                 :: nspDisp
   PetscReal                                          :: rtol,dtol,atol
          
   PetscBool                                          :: flg
   Character(len=MEF90_MXSTRLEN)                      :: IOBuffer
   Character(len=MEF90_MXSTRLEN)                      :: setName,setprefix
   Character(len=MXSTLN),Dimension(:),Pointer         :: nameG,nameV,nameC
   Integer                                            :: numfield
   
   Integer                                            :: step
   Type(Vec)                                          :: localVec
   PetscInt                                           :: dim
      
   !!! Initialize MEF90
   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   Call MEF90_Initialize(ierr)

   !!! Get all MEF90-wide options
   Call MEF90Ctx_Create(PETSC_COMM_WORLD,MEF90Ctx,MEF90DefaultGlobalOptions,ierr);CHKERRQ(ierr)
   Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)

   !!! Get DM from mesh
   Call MEF90Ctx_GetDMMeshEXO(MEF90Ctx,Mesh,ierr);CHKERRQ(ierr)
   Call DMMeshGetDimension(Mesh,dim,ierr);CHKERRQ(ierr)
   Call DMMeshSetMaxDof(Mesh,dim,ierr);CHKERRQ(ierr) 
   Call DMSetBlockSize(Mesh,dim,ierr);CHKERRQ(ierr)
   
   !!! Open output file
   Call MEF90Ctx_OpenEXO(MEF90Ctx,Mesh,ierr)
   
   !!! Create DefMech context, get all DefMech options
   Call MEF90DefMechCtx_Create(MEF90DefMechCtx,Mesh,MEF90Ctx,ierr);CHKERRQ(ierr)
   If (dim == 2) Then
      Call MEF90DefMechCtx_SetFromOptions(MEF90DefMechCtx,PETSC_NULL_CHARACTER,MEF90DefMechDefaultGlobalOptions2D, &
                                        MEF90DefMechDefaultCellSetOptions,MEF90DefMechDefaultVertexSetOptions,ierr)
   Else
      Call MEF90DefMechCtx_SetFromOptions(MEF90DefMechCtx,PETSC_NULL_CHARACTER,MEF90DefMechDefaultGlobalOptions3D, &
                                        MEF90DefMechDefaultCellSetOptions,MEF90DefMechDefaultVertexSetOptions,ierr)
   End If
   Call PetscBagGetDataMEF90DefMechCtxGlobalOptions(MEF90DefMechCtx%GlobalOptionsBag,MEF90DefMechGlobalOptions,ierr);CHKERRQ(ierr)
   
   !!! Get material properties bags
   If (dim == 2) Then
      Call MEF90MatPropBag_SetFromOptions(MEF90MatPropBag,MEF90DefMechCtx%DMVect,MEF90_Mathium2D,MEF90Ctx,ierr)
   Else
      Call MEF90MatPropBag_SetFromOptions(MEF90MatPropBag,MEF90DefMechCtx%DMVect,MEF90_Mathium3D,MEF90Ctx,ierr)
   End If   
   MEF90DefMechCtx%MaterialPropertiesBag => MEF90MatPropBag

   Call MEF90Ctx_GetTime(MEF90Ctx,time,ierr)

   !!! Set the data layout
   Call MEF90DefMechCtx_SetSections(MEF90DefMechCtx,ierr)

   !!! Create vectors
   Call MEF90DefMechCtx_CreateVectors(MEF90DefMechCtx,ierr)

   !!! 
   !!! Create SNES or TS, Mat and set KSP default options
   !!!
   Call DMCreateMatrix(MEF90DefMechCtx%DMVect,MATAIJ,matDisp,iErr);CHKERRQ(iErr)
   Call MatSetOptionsPrefix(matDisp,"Disp_",ierr);CHKERRQ(ierr)
   Call MatSetOption(matDisp,MAT_SPD,PETSC_TRUE,ierr);CHKERRQ(ierr)
   Call MatSetOption(matDisp,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr);CHKERRQ(ierr)
   Call MatSetOption(matDisp,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr);CHKERRQ(ierr)
   If (MEF90DefMechGlobalOptions%addDisplacementNullSpace) Then
      Call DMMeshGetSectionReal(MEF90DefMechCtx%DMVect,'coordinates',coordSec,ierr);CHKERRQ(ierr)
      Call DMMeshCreateGlobalScatter(MEF90DefMechCtx%DMVect,coordSec,ScatterSecToVec,ierr);CHKERRQ(ierr)
      Call DMCreateGlobalVector(MEF90DefMechCtx%DMVect,coordVec,ierr)
      Call SectionRealToVec(coordSec,ScatterSecToVec,SCATTER_FORWARD,coordVec,ierr);CHKERRQ(ierr)
      Call MatNullSpaceCreateRigidBody(coordVec,nspDisp,ierr);CHKERRQ(ierr)
      Call MatSetNearNullSpace(matDisp,nspDisp,ierr);CHKERRQ(ierr)
      !!!Call MatSetNullSpace(matDisp,nspDisp,ierr);CHKERRQ(ierr)
      Call MatNullSpaceDestroy(nspDisp,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(coordSec,ierr);CHKERRQ(ierr)
      Call VecDestroy(coordVec,ierr);CHKERRQ(ierr)
      Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)
   End If


   Call MatSetFromOptions(matDisp,ierr);CHKERRQ(ierr)

   If (MEF90DefMechGlobalOptions%mode == MEF90DefMech_ModeQuasiStatic) Then
      Call SNESCreate(PETSC_COMM_WORLD,snesDisp,ierr);CHKERRQ(ierr)
      Call SNESSetApplicationContext(snesDisp,MEF90DefMechCtx,ierr);CHKERRQ(ierr)
      Call SNESSetDM(snesDisp,MEF90DefMechCtx%DMVect,ierr);CHKERRQ(ierr)
      Call SNESSetOptionsPrefix(snesDisp,'Disp_',ierr);CHKERRQ(ierr)
      Call SNESSetType(snesDisp,SNESKSPONLY,ierr);CHKERRQ(ierr)

      Call SNESSetFunction(snesDisp,residualDisplacement,MEF90DefMechOperator,MEF90DefMechCtx,ierr);CHKERRQ(ierr)
      Call SNESSetJacobian(snesDisp,matDisp,matDisp,MEF90DefMechBilinearForm,MEF90DefMechCtx,ierr);CHKERRQ(ierr)
      !atol = 1.0D-10
      !rtol = 1.0D-10
      !Call SNESSetTolerances(snesDisp,atol,PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_INTEGER,ierr);CHKERRQ(ierr)
      Call SNESSetFromOptions(snesDisp,ierr);CHKERRQ(ierr)
      If (MEF90GlobalOptions%verbose > 0) Then
         Call SNESView(snesDisp,PETSC_VIEWER_STDOUT_WORLD,ierr)
      End If
   End If
   !!! 
   !!! Set some KSP options
   !!!
   Call SNESGetKSP(snesDisp,kspDisp,ierr);CHKERRQ(ierr)
   Call KSPSetType(kspDisp,KSPCG,ierr);CHKERRQ(ierr)
   Call KSPSetInitialGuessNonzero(kspDisp,PETSC_TRUE,ierr);CHKERRQ(ierr)
   rtol = 1.0D-8
   atol = 1.0D-8
   dtol = 1.0D+10
   Call KSPSetTolerances(kspDisp,rtol,atol,dtol,PETSC_DEFAULT_INTEGER,ierr);CHKERRQ(ierr)
   Call KSPSetFromOptions(kspDisp,ierr);CHKERRQ(ierr)


   ! set coordinates in PC for GAMG
   Call KSPGetPC(kspDisp,pcDisp,ierr);CHKERRQ(ierr)
   Call DMMeshGetCoordinatesF90(MEF90DefMechCtx%DMVect,coordPtr,ierr);CHKERRQ(ierr)
   Allocate(coordPCPtr(size(CoordPtr)))
   coordPCPtr = reshape(transpose(coordPtr),[size(CoordPtr)])
   coordPCPtr = reshape((coordPtr),[size(CoordPtr)])
   !Call PCSetCoordinates(pcDisp,dim,size(coordPtr),coordPCPtr,ierr);CHKERRQ(ierr)
   DeAllocate(coordPCPtr)
   Call DMMeshRestoreCoordinatesF90(MEF90DefMechCtx%DMVect,coordPtr,ierr);CHKERRQ(ierr)
   Call PCSetFromOptions(pcDisp,ierr);CHKERRQ(ierr)
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
      Allocate(nameG(2))
      nameG(1) = "Energy"
      nameG(2) = "Work"
   
      numfield = max(MEF90DefMechGlobalOptions%displacementOffset+dim-1, &
                     MEF90DefMechGlobalOptions%damageOffset,&
                     MEF90DefMechGlobalOptions%boundaryDisplacementOffset+dim-1,&
                     MEF90DefMechGlobalOptions%boundaryDamageOffset,&
                     MEF90DefMechGlobalOptions%temperatureOffset)
      Allocate(nameV(numfield))

      nameV = "empty"
      nameV(MEF90DefMechGlobalOptions%displacementOffset+0)            = "Displacement_X"
      nameV(MEF90DefMechGlobalOptions%displacementOffset+1)            = "Displacement_Y"
      nameV(MEF90DefMechGlobalOptions%boundaryDisplacementOffset+0)    = "Boundary_Displacement_X"
      nameV(MEF90DefMechGlobalOptions%boundaryDisplacementOffset+1)    = "Boundary_Displacement_Y"
      If (dim == 3) Then
         nameV(MEF90DefMechGlobalOptions%displacementOffset+2)         = "Displacement_Z"
         nameV(MEF90DefMechGlobalOptions%boundaryDisplacementOffset+2) = "Boundary_Displacement_Z"
      End If
      nameV(MEF90DefMechGlobalOptions%damageOffset)                    = "Damage"
      nameV(MEF90DefMechGlobalOptions%boundaryDamageOffset)            = "Boundary_Damage"
      nameV(MEF90DefMechGlobalOptions%temperatureOffset)               = "Temperature"
                     
      numfield = max(MEF90DefMechGlobalOptions%forceOffset+dim-1,&
                     MEF90DefMechGlobalOptions%pressureForceOffset,&
                     MEF90DefMechGlobalOptions%StressOffset+(dim*(dim+1))/2-1,&
                     MEF90DefMechGlobalOptions%plasticStrainOffset+(dim*(dim+1))/2-1)
      Allocate(nameC(numfield))
      nameC = "empty"
      nameC(MEF90DefMechGlobalOptions%forceOffset+0)                 = "Force_X"
      nameC(MEF90DefMechGlobalOptions%forceOffset+1)                 = "Force_Y"
      If (dim == 3) Then
         nameC(MEF90DefMechGlobalOptions%forceOffset+2)              = "Force_Z"
      End If

      nameC(MEF90DefMechGlobalOptions%pressureForceOffset)           = "Pressure_Force"
      If (dim == 2) Then
         nameC(MEF90DefMechGlobalOptions%stressOffset+0)             = "Stress_XX"
         nameC(MEF90DefMechGlobalOptions%stressOffset+1)             = "Stress_YY"
         nameC(MEF90DefMechGlobalOptions%stressOffset+2)             = "Stress_XY"
         nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+0)      = "plasticStrainOffset_XX"
         nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+1)      = "plasticStrainOffset_YY"
         nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+2)      = "plasticStrainOffset_XY"
      Else
         nameC(MEF90DefMechGlobalOptions%stressOffset+0)             = "Stress_XX"
         nameC(MEF90DefMechGlobalOptions%stressOffset+1)             = "Stress_YY"
         nameC(MEF90DefMechGlobalOptions%stressOffset+2)             = "Stress_ZZ"
         nameC(MEF90DefMechGlobalOptions%stressOffset+3)             = "Stress_YZ"
         nameC(MEF90DefMechGlobalOptions%stressOffset+4)             = "Stress_XZ"
         nameC(MEF90DefMechGlobalOptions%stressOffset+5)             = "Stress_XY"
         nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+0)      = "plasticStrain_XX"
         nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+1)      = "plasticStrain_YY"
         nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+2)      = "plasticStrain_ZZ"
         nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+3)      = "plasticStrain_YZ"
         nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+4)      = "plasticStrain_XZ"
         nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+5)      = "plasticStrain_XY"
      End If
      Call MEF90EXOFormat(MEF90Ctx%fileEXOUNIT,nameG,nameC,nameV,ierr)
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
   End If

   Call MEF90DefMechCtx_DestroyVectors(MEF90DefMechCtx,ierr)
!   Call VecDestroy(residualDisplacement,ierr);CHKERRQ(ierr)   

   Call DMDestroy(Mesh,ierr);CHKERRQ(ierr)
   DeAllocate(time)
   DeAllocate(energy)
   DeAllocate(work)
   Call PetscLogView(MEF90Ctx%logViewer,ierr);CHKERRQ(ierr)
   Call MEF90DefMechCtx_Destroy(MEF90DefMechCtx,ierr);CHKERRQ(ierr)
   Call MEF90Ctx_CloseEXO(MEF90Ctx,ierr)
   Call MEF90_Finalize(ierr)
   Call PetscFinalize(ierr)
End Program ThermoElasticity
