#include "../MEF90/mef90.inc"
Program ThermoElasticity
#include <finclude/petscdef.h>
   Use m_MEF90
   Use m_MEF90_DefMech
   Use petsc
   Implicit NONE   
PetscInt n1,n2
   PetscErrorCode                                     :: ierr
   Type(MEF90DefMechCtx_Type)                         :: MEF90DefMechCtx
   Type(MEF90DefMechGlobalOptions_Type),Parameter     :: MEF90DefMechDefaultGlobalOptions2D = MEF90DefMechGlobalOptions_Type( &
                                                         MEF90DefMech_ModeQuasiStatic, & ! mode
                                                         PETSC_FALSE,         & ! addDisplacementNullSpace
                                                         1,                   & ! DisplacementOffset
                                                         3,                   & ! DamageOffset
                                                         4,                   & ! boundaryDisplacementOffset
                                                         5,                   & ! boundaryDamageOffset
                                                         1,                   & ! ForceOffset
                                                         3,                   & ! pressureForceOffset
                                                         4,                   & ! plasticStrainOffset
                                                         7,                   & ! StressOffset
                                                         MEF90Scaling_Linear, & ! boundaryDisplacementScaling
                                                         MEF90Scaling_Linear, & ! ForceScaling
                                                         MEF90Scaling_Linear)   ! pressureForceScaling
   Type(MEF90DefMechGlobalOptions_Type),Parameter     :: MEF90DefMechDefaultGlobalOptions3D = MEF90DefMechGlobalOptions_Type( &
                                                         MEF90DefMech_ModeQuasiStatic, & ! mode
                                                         PETSC_FALSE,         & ! addDisplacementNullSpace
                                                         1,                   & ! DisplacementOffset
                                                         4,                   & ! DamageOffset
                                                         5,                   & ! boundaryDisplacementOffset
                                                         8,                   & ! boundaryDamageOffset
                                                         1,                   & ! ForceOffset
                                                         4,                   & ! pressureForceOffset
                                                         5,                   & ! plasticStrainOffset
                                                         11,                  & ! StressOffset
                                                         MEF90Scaling_Linear, & ! boundaryDisplacementScaling
                                                         MEF90Scaling_Linear, & ! ForceScaling
                                                         MEF90Scaling_Linear)   ! pressureForceScaling
   Type(MEF90DefMechCellSetOptions_Type),Parameter    :: MEF90DefMechDefaultCellSetOptions = MEF90DefMechCellSetOptions_Type( &
                                                         -1,                                   & ! elemTypeShortID will be overriden
                                                         (/0.0_Kr,0.0_Kr,0.0_Kr/),             & ! force
                                                         0.0_Kr,                               & ! pressureForce
                                                         MEF90DefMech_defectLawElasticity,     & ! defect law
                                                         (/PETSC_TRUE,PETSC_TRUE,PETSC_TRUE/), & ! Has Displacement BC
                                                         0.0_Kr,                               & ! boundary Displacement
                                                         PETSC_FALSE,                          & ! Has Damage BC
                                                         0.0_Kr)                                 ! Boundary Damage
   Type(MEF90DefMechVertexSetOptions_Type),Parameter  :: MEF90DefMechDefaultVertexSetOptions = MEF90DefMechVertexSetOptions_Type( &
                                                         (/PETSC_TRUE,PETSC_TRUE,PETSC_TRUE/),     & ! Has Displacement BC
                                                         0.0_Kr,        & ! boundary Displacement
                                                         PETSC_FALSE,   & ! Has Damage BC
                                                         0.0_Kr)          ! boundary Damage
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
   Type(SectionReal)                                  :: defaultSection
   Type(Vec),target                                   :: displacement
   Type(Vec),target                                   :: boundaryDisplacement
   Type(Vec),target                                   :: force
   Type(Vec),target                                   :: pressureForce
   Type(Vec)                                          :: residualDisplacement
   PetscReal,Dimension(:),Pointer                     :: time,energy,work

   Type(SNES)                                         :: snesDisp
   Type(KSP)                                          :: kspDisp
   Type(PC)                                           :: pcDisp
   Type(Mat)                                          :: matDisp
   Type(MatNullSpace)                                 :: nspDisp
   PetscReal                                          :: rtol
          
   PetscBool                                          :: flg
   Character(len=MEF90_MXSTRLEN)                      :: IOBuffer
   Character(len=MEF90_MXSTRLEN)                      :: setName,setprefix
   Character(len=MXSTLN),Dimension(:),Pointer         :: nameG,nameC,nameV
   Integer                                            :: numfield
   
   Integer                                            :: step
   Type(Vec)                                          :: localVec
   PetscInt                                           :: dim
      
   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   Call MEF90_Initialize(ierr)

   !!! Initialize MEF90
   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   Call MEF90_Initialize(ierr)

   !!! Get all MEF90-wide options
   Call MEF90Ctx_Create(PETSC_COMM_WORLD,MEF90Ctx,MEF90DefaultGlobalOptions,ierr);CHKERRQ(ierr)
   Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr);CHKERRQ(ierr)

   !!! Get DM from mesh
   Call MEF90Ctx_GetDMMeshEXO(MEF90Ctx,Mesh,ierr);CHKERRQ(ierr)
   Call DMMeshGetDimension(Mesh,dim,ierr);CHKERRQ(ierr)
   Call DMMeshSetMaxDof(mesh,dim,ierr);CHKERRQ(ierr) 
   
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

   !!! Create default section matching element type
   Call DMMeshGetVertexSectionReal(MEF90DefMechCtx%DMVect,"default",dim,defaultSection,ierr);CHKERRQ(ierr)
   Call DMMeshSetSectionReal(MEF90DefMechCtx%DMVect,"default",defaultSection,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(defaultSection,ierr);CHKERRQ(ierr)
   
   Call DMMeshGetCellSectionReal(MEF90DefMechCtx%cellDMVect,"default",dim,defaultSection,ierr);CHKERRQ(ierr)
   Call DMMeshSetSectionReal(MEF90DefMechCtx%cellDMVect,"default",defaultSection,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(defaultSection,ierr);CHKERRQ(ierr)
      
   Call DMMeshGetVertexSectionReal(MEF90DefMechCtx%DMScal,"default",1,defaultSection,ierr);CHKERRQ(ierr)
   Call DMMeshSetSectionReal(MEF90DefMechCtx%DMScal,"default",defaultSection,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(defaultSection,ierr);CHKERRQ(ierr)

   Call DMMeshGetCellSectionReal(MEF90DefMechCtx%cellDMScal,"default",1,defaultSection,ierr);CHKERRQ(ierr)
   Call DMMeshSetSectionReal(MEF90DefMechCtx%cellDMScal,"default",defaultSection,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(defaultSection,ierr);CHKERRQ(ierr)

   Call DMMeshGetCellSectionReal(MEF90DefMechCtx%cellDMMatS,"default",(dim*(dim+1))/2,defaultSection,ierr);CHKERRQ(ierr)
   Call DMMeshSetSectionReal(MEF90DefMechCtx%cellDMMatS,"default",defaultSection,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(defaultSection,ierr);CHKERRQ(ierr)

   !!! Create vectors
   Call DMCreateGlobalVector(MEF90DefMechCtx%DMVect,Displacement,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(Displacement,"Displacement",ierr);CHKERRQ(ierr)
   MEF90DefMechCtx%Displacement => Displacement

   Call DMCreateGlobalVector(MEF90DefMechCtx%DMVect,boundaryDisplacement,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(boundaryDisplacement,"boundary Displacement",ierr);CHKERRQ(ierr)
   MEF90DefMechCtx%boundaryDisplacement => boundaryDisplacement
   
   Call DMCreateGlobalVector(MEF90DefMechCtx%cellDMVect,force,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(force,"Force",ierr);CHKERRQ(ierr)
   MEF90DefMechCtx%Force => Force

   Call DMCreateGlobalVector(MEF90DefMechCtx%cellDMScal,pressureForce,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(pressureForce,"Pressure Force",ierr);CHKERRQ(ierr)
   MEF90DefMechCtx%pressureForce => pressureForce

   Call DMCreateGlobalVector(MEF90DefMechCtx%DMVect,residualDisplacement,ierr);CHKERRQ(ierr)
   Call PetscObjectSetName(residualDisplacement,"residualDisplacement",ierr);CHKERRQ(ierr)
   !!! 
   !!! Create SNES or TS, Mat and set KSP default options
   !!!
   Call DMCreateMatrix(MEF90DefMechCtx%DMVect,MATAIJ,matDisp,iErr);CHKERRQ(iErr)
   Call MatSetOptionsPrefix(matDisp,"Disp_",ierr);CHKERRQ(ierr)
   Call MatSetOption(matDisp,MAT_SPD,PETSC_TRUE,ierr);CHKERRQ(ierr)
   Call MatSetOption(matDisp,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr);CHKERRQ(ierr)
   Call MatSetOption(matDisp,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr);CHKERRQ(ierr)
   Call MatSetFromOptions(matDisp,ierr);CHKERRQ(ierr)

   If (MEF90DefMechGlobalOptions%mode == MEF90DefMech_ModeQuasiSTatic) Then
      Call SNESCreate(PETSC_COMM_WORLD,snesDisp,ierr);CHKERRQ(ierr)
      Call SNESSetApplicationContext(snesDisp,MEF90DefMechCtx,ierr);CHKERRQ(ierr)
      Call SNESSetDM(snesDisp,MEF90DefMechCtx%DMVect,ierr);CHKERRQ(ierr)
      Call SNESSetOptionsPrefix(snesDisp,'Disp_',ierr);CHKERRQ(ierr)

      !!!Call SNESSetFunction(snesDisp,residual,MEF90DefMechOperator,MEF90DefMechCtx,ierr);CHKERRQ(ierr)
      !!!Call SNESSetJacobian(snesDisp,matDisp,matDisp,MEF90DefMechBilinearForm,MEF90DefMechCtx,ierr);CHKERRQ(ierr)
      Call SNESSetFromOptions(snesDisp,ierr);CHKERRQ(ierr)
      If (MEF90GlobalOptions%verbose > 0) Then
         Call SNESView(snesDisp,PETSC_VIEWER_STDOUT_WORLD,ierr)
      End If
   End If
   If (MEF90DefMechGlobalOptions%addDisplacementNullSpace) Then
      !!!Call MatNullSpaceCreateRigidBody(Vec coords,nspDisp,ierr);CHKERRQ(ierr))
      !!!Call MatSetNullSpace(matDisp,nspDisp,ierr);CHKERRQ(ierr)
   End If
   !!! 
   !!! Set some KSP options
   !!!
   Call SNESGetKSP(snesDisp,kspDisp,ierr);CHKERRQ(ierr)
   Call KSPSetType(kspDisp,KSPCG,ierr);CHKERRQ(ierr)
   Call KSPSetInitialGuessNonzero(kspDisp,PETSC_TRUE,ierr);CHKERRQ(ierr)
   If (MEF90DefMechGlobalOptions%addDisplacementNullSpace) Then
      !!!Call KSPSetNullSpace(kspDisp,nspDisp,ierr);CHKERRQ(ierr)
   End If
   rtol = 1.0D-8
   Call KSPSetTolerances(kspDisp,rtol,PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_INTEGER,ierr);CHKERRQ(ierr)
   Call KSPSetFromOptions(kspDisp,ierr);CHKERRQ(ierr)
   
   !!! 
   !!! Allocate array of works and energies
   !!!
   Allocate(energy(size(MEF90DefMechCtx%CellSetOptionsBag)))
   Allocate(work(size(MEF90DefMechCtx%CellSetOptionsBag)))

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
   
      numfield = max(MEF90DefMechGlobalOptions%displacementOffset+dim, &
                     MEF90DefMechGlobalOptions%damageOffset,&
                     MEF90DefMechGlobalOptions%boundaryDisplacementOffset+dim,&
                     MEF90DefMechGlobalOptions%boundaryDamageOffset)-1
      Allocate(nameV(numfield))
      nameV = "empty"
      nameV(MEF90DefMechGlobalOptions%displacementOffset+0)         = "Displacement X"
      nameV(MEF90DefMechGlobalOptions%displacementOffset+1)         = "Displacement Y"
      nameV(MEF90DefMechGlobalOptions%boundaryDisplacementOffset+0) = "Boundary Displacement X"
      nameV(MEF90DefMechGlobalOptions%boundaryDisplacementOffset+1) = "Boundary Displacement Y"
      If (dim == 3) Then
         nameV(MEF90DefMechGlobalOptions%displacementOffset+2)         = "Displacement Z"
         nameV(MEF90DefMechGlobalOptions%boundaryDisplacementOffset+2) = "Boundary Displacement Z"
      End If
      
      nameV(MEF90DefMechGlobalOptions%damageOffset)                    = "Damage"
      nameV(MEF90DefMechGlobalOptions%boundaryDamageOffset)            = "Boundary Damage"
                     
      numfield = max(MEF90DefMechGlobalOptions%forceOffset+dim,&
                     MEF90DefMechGlobalOptions%pressureForceOffset,&
                     MEF90DefMechGlobalOptions%StressOffset+dim*(dim+1)/2,&
                     MEF90DefMechGlobalOptions%plasticStrainOffset+dim*(dim+1)/2)-1
      Allocate(nameC(numfield))
      nameC = "empty"
      nameC(MEF90DefMechGlobalOptions%forceOffset+0)                 = "Force X"
      nameC(MEF90DefMechGlobalOptions%forceOffset+1)                 = "Force Y"
      If (dim == 3) Then
         nameC(MEF90DefMechGlobalOptions%forceOffset+2)              = "Force Z"
      End If

      nameC(MEF90DefMechGlobalOptions%pressureForceOffset)           = "Pressure Force"
      If (dim == 2) Then
         nameC(MEF90DefMechGlobalOptions%stressOffset+0)             = "Stress XX"
         nameC(MEF90DefMechGlobalOptions%stressOffset+1)             = "Stress YY"
         nameC(MEF90DefMechGlobalOptions%stressOffset+2)             = "Stress XY"
         nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+0)      = "plastic Strain XX"
         nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+1)      = "plastic Strain YY"
         nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+2)      = "plastic Strain XY"
      Else
         nameC(MEF90DefMechGlobalOptions%stressOffset+0)             = "Stress XX"
         nameC(MEF90DefMechGlobalOptions%stressOffset+1)             = "Stress YY"
         nameC(MEF90DefMechGlobalOptions%stressOffset+2)             = "Stress ZZ"
         nameC(MEF90DefMechGlobalOptions%stressOffset+3)             = "Stress YZ"
         nameC(MEF90DefMechGlobalOptions%stressOffset+4)             = "Stress XZ"
         nameC(MEF90DefMechGlobalOptions%stressOffset+5)             = "Stress XY"
         nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+0)      = "plastic Strain XX"
         nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+1)      = "plastic Strain YY"
         nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+2)      = "plastic Strain ZZ"
         nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+3)      = "plastic Strain YZ"
         nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+4)      = "plastic Strain XZ"
         nameC(MEF90DefMechGlobalOptions%plasticStrainOffset+5)      = "plastic Strain XY"
      End If
      Call MEF90EXOFormat(MEF90Ctx%fileEXOUNIT,nameG,nameC,nameV,ierr)
   End If
   
   !!!
   !!! Actual computations / time stepping
   !!!
   If (MEF90DefMechGlobalOptions%mode == MEF90DefMech_ModeQuasiStatic) Then
      Call MEF90DefMechSetTransients(MEF90DefMechCtx,1,time(1),ierr)
      Do step = 1,MEF90GlobalOptions%timeNumStep
         Write(IOBuffer,100) step,time(step)
         Call PetscPrintf(MEF90Ctx%comm,IOBuffer,ierr);CHKERRQ(ierr)

         !!! Update fields
         Call MEF90DefMechSetTransients(MEF90DefMechCtx,step,time(step),ierr)

         !!! Solve SNES
         !!!Call SNESSolve(snesDisp,PETSC_NULL_OBJECT,Displacement,ierr);CHKERRQ(ierr)
         
         !!! Compute energies
         !Call MEF90DefMechEnergy(Displacement,time(step),MEF90DefMechCtx,energy,work,ierr);CHKERRQ(ierr)
         Call DMmeshGetLabelIdIS(MEF90DefMechCtx%DMVect,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
         Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
         Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
         Do set = 1, size(setID)
            Write(IOBuffer,101) setID(set),energy(set),work(set)
            Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
         End Do
         Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
         Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
         Write(IOBuffer,102) sum(energy),sum(work)
         Call PetscPrintf(MEF90Ctx%Comm,IOBuffer,ierr);CHKERRQ(ierr)
     
         
         !!! Save results and boundary Values
         Call DMGetLocalVector(MEF90DefMechCtx%DMVect,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalBegin(MEF90DefMechCtx%DMVect,MEF90DefMechCtx%Displacement,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalEnd(MEF90DefMechCtx%DMVect,MEF90DefMechCtx%Displacement,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call VecViewExodusVertex(MEF90DefMechCtx%DMVect,localVec,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
                                  MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%displacementOffset,ierr);CHKERRQ(ierr)
         Call DMRestoreLocalVector(MEF90DefMechCtx%DMVect,localVec,ierr);CHKERRQ(ierr)

         Call DMGetLocalVector(MEF90DefMechCtx%DMVect,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalBegin(MEF90DefMechCtx%DMVect,MEF90DefMechCtx%boundaryDisplacement,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalEnd(MEF90DefMechCtx%DMVect,MEF90DefMechCtx%boundaryDisplacement,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call VecViewExodusVertex(MEF90DefMechCtx%DMVect,localVec,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
                                  MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%boundaryDisplacementOffset,ierr);CHKERRQ(ierr)
         Call DMRestoreLocalVector(MEF90DefMechCtx%DMVect,localVec,ierr);CHKERRQ(ierr)

         Call DMGetLocalVector(MEF90DefMechCtx%cellDMVect,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalBegin(MEF90DefMechCtx%cellDMVect,MEF90DefMechCtx%Force,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call DMGlobalToLocalEnd(MEF90DefMechCtx%cellDMVect,MEF90DefMechCtx%Force,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         Call VecViewExodusCell(MEF90DefMechCtx%cellDMVect,localVec,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
                                MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%forceOffset,ierr);CHKERRQ(ierr)
         Call DMRestoreLocalVector(MEF90DefMechCtx%cellDMVect,localVec,ierr);CHKERRQ(ierr)

         Call DMGetLocalVector(MEF90DefMechCtx%cellDMScal,localVec,ierr);CHKERRQ(ierr)
Call VecGetSize(localvec,n1, ierr)
Call VecGetBlockSize(localvec,n2, ierr)
Write(*,*) "localvec size: ", n1,n2
Call VecGetSize(MEF90DefMechCtx%pressureForce,n1, ierr)
Call VecGetBlockSize(MEF90DefMechCtx%pressureForce,n2, ierr)
Write(*,*) "MEF90DefMechCtx%pressureForce size: ", n1,n2
Call DMGetBlockSize(MEF90DefMechCtx%cellDMScal,n1,ierr)
Write(*,*) "MEF90DefMechCtx%cellDMScal size: ", n1

Call DMView(MEF90DefMechCtx%cellDMScal,PETSC_VIEWER_STDOUT_WORLD,ierr)

!         Call DMGlobalToLocalBegin(MEF90DefMechCtx%cellDMScal,MEF90DefMechCtx%pressureForce,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
!         Call DMGlobalToLocalEnd(MEF90DefMechCtx%cellDMScal,MEF90DefMechCtx%pressureForce,INSERT_VALUES,localVec,ierr);CHKERRQ(ierr)
         !Call VecViewExodusCell(MEF90DefMechCtx%cellDMScal,localVec,MEF90DefMechCtx%MEF90Ctx%IOcomm, &
         !                       MEF90DefMechCtx%MEF90Ctx%fileExoUnit,step,MEF90DefMechGlobalOptions%pressureForceOffset,ierr);CHKERRQ(ierr)
         Call DMRestoreLocalVector(MEF90DefMechCtx%cellDMScal,localVec,ierr);CHKERRQ(ierr)
      End Do
   End If
100 Format("Solving steady state step ",I4,", t=",ES12.5,"\n")
101 Format("cell set ",I4," thermal energy: ",ES12.5," fluxes work: ",ES12.5,"\n")
102 Format("======= Total thermal energy: ",ES12.5," fluxes work: ",ES12.5,"\n")
   !!! Clean up and exit nicely
   If (MEF90DefMechGlobalOptions%mode == MEF90DefMech_ModeQuasiStatic) Then
      Call SNESDestroy(snesDisp,ierr);CHKERRQ(ierr)
   End If

   If (Associated(MEF90DefMechCtx%boundaryDisplacement)) Then 
      Call VecDestroy(MEF90DefMechCtx%boundaryDisplacement,ierr);CHKERRQ(ierr)
      Nullify(MEF90DefMechCtx%boundaryDisplacement)
   End If   
   If (Associated(MEF90DefMechCtx%boundaryDamage)) Then 
      Call VecDestroy(MEF90DefMechCtx%boundaryDamage,ierr);CHKERRQ(ierr)
      Nullify(MEF90DefMechCtx%boundaryDamage)
   End If
   
   If (Associated(MEF90DefMechCtx%force)) Then 
      Call VecDestroy(MEF90DefMechCtx%force,ierr);CHKERRQ(ierr)
      Nullify(MEF90DefMechCtx%force)
   End If
   If (Associated(MEF90DefMechCtx%pressureForce)) Then 
      Call VecDestroy(MEF90DefMechCtx%pressureForce,ierr);CHKERRQ(ierr)
      Nullify(MEF90DefMechCtx%pressureForce)
   End If

   Call VecDestroy(Displacement,ierr);CHKERRQ(ierr)
   Call VecDestroy(boundaryDisplacement,ierr);CHKERRQ(ierr)
   Call VecDestroy(residualDisplacement,ierr);CHKERRQ(ierr)   
   Call VecDestroy(force,ierr);CHKERRQ(ierr)   
   Call VecDestroy(pressureForce,ierr);CHKERRQ(ierr)   
   Call DMDestroy(Mesh,ierr);CHKERRQ(ierr)
   DeAllocate(time)
   DeAllocate(energy)
   DeAllocate(work)
   Call PetscLogView(MEF90Ctx%logViewer,ierr);CHKERRQ(ierr)
   Call MEF90DefMechCtx_Destroy(MEF90DefMechCtx,ierr);CHKERRQ(ierr)
   Call MEF90Ctx_CloseEXO(MEF90Ctx,ierr)
   Call MEF90_Finalize(ierr)
   Call PetscFinalize()
End Program ThermoElasticity
