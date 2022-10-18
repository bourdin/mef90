#include "../MEF90/mef90.inc"
Module m_MEF90_HeatXfer
#include "petsc/finclude/petsc.h"
   Use m_MEF90
   Use m_MEF90_HeatXferCtx
   Use m_MEF90_HeatXferAssembly2D, &
      MEF90HeatXferEnergy2D       => MEF90HeatXferEnergy, &
      MEF90HeatXferOperator2D     => MEF90HeatXferOperator, &
      MEF90HeatXferBilinearForm2D => MEF90HeatXferBilinearForm, &
      MEF90HeatXferIFunction2D    => MEF90HeatXferIFunction, &
      MEF90HeatXferIJacobian2D    => MEF90HeatXferIJacobian
   Use m_MEF90_HeatXferAssembly3D, &
      MEF90HeatXferEnergy3D       => MEF90HeatXFerEnergy, &
      MEF90HeatXferOperator3D     => MEF90HeatXferOperator, &
      MEF90HeatXferBilinearForm3D => MEF90HeatXferBilinearForm, &
      MEF90HeatXferIFunction3D    => MEF90HeatXferIFunction, &
      MEF90HeatXferIJacobian3D    => MEF90HeatXferIJacobian

   Implicit none

   !Private
   Public MEF90HeatXferOperator
   Public MEF90HeatXferBilinearForm
   Public MEF90HeatXferEnergy
   Public MEF90HeatXferUpdateTransients
   Public MEF90HeatXferIFunction
   Public MEF90HeatXferIJacobian
   Public MEF90HeatXferViewEXO
   Public MEF90HeatXferCreateSNES
   Public MEF90HeatXferCreateTS
Contains

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferUpdateTransients"
!!!
!!!  
!!!  MEF90HeatXferUpdateTransients: Update all transient data (boundary / external temperature and fluxes)
!!!                              using the proper scaling law
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90HeatXferUpdateTransients(MEF90HeatXferCtx,step,time,ierr)
      Type(MEF90HeatXferCtx_Type),Intent(INOUT)       :: MEF90HeatXferCtx
      PetscInt,Intent(IN)                             :: step
      PetscReal,Intent(IN)                            :: time
      PetscErrorCode,Intent(OUT)                      :: ierr
   
      Type(MEF90HeatXferGlobalOptions_Type),pointer   :: MEF90HeatXferGlobalOptions
      Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90GlobalOptions
      Type(tDM)                                       :: dm
      Type(tVec)                                      :: tmpVec

      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90HeatXferCtx%MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr))
      PetscCall(PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr))

      PetscCall(VecGetDM(MEF90HeatXferCtx%temperatureLocal,dm,ierr))

      Select case (MEF90HeatXferGlobalOptions%boundaryTemperatureScaling)
      Case (MEF90Scaling_File)
         PetscCall(DMGetLocalVector(dm,tmpVec,ierr))
         PetscCall(PetscObjectSetName(tmpVec,"Temperature",ierr))
         PetscCall(MEF90EXOVecLoad(tmpVec,MEF90HeatXferCtx%temperatureToIOSF,MEF90HeatXferCtx%IOToTemperatureSF,MEF90HeatXferCtx%MEF90Ctx%resultViewer,step,ierr))
         PetscCall(MEF90VecCopySF(tmpVec,MEF90HeatXferCtx%temperatureLocal,MEF90HeatXferCtx%boundaryToTemperatureSF,ierr))
         PetscCall(DMRestoreLocalVector(dm,tmpVec,ierr))
      Case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90HeatXferCtx%temperatureLocal,time,ierr))
      Case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetBCValuesFromOptions(MEF90HeatXferCtx%temperatureLocal,1.0_Kr,ierr))
      End Select

      Select case (MEF90HeatXferGlobalOptions%externalTemperatureScaling)
      Case (MEF90Scaling_File)
         PetscCall(MEF90EXOVecLoad(MEF90HeatXferCtx%externalTemperatureLocal,MEF90HeatXferCtx%externalTemperatureToIOSF,MEF90HeatXferCtx%IOToExternalTemperatureSF,MEF90HeatXferCtx%MEF90Ctx%resultViewer,step,ierr))    
      Case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90HeatXferCtx%externalTemperatureLocal,time,ierr))
      Case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90HeatXferCtx%externalTemperatureLocal,1.0_Kr,ierr))
      End Select

      Select case (MEF90HeatXferGlobalOptions%fluxScaling)
      Case (MEF90Scaling_File)
         PetscCall(MEF90EXOVecLoad(MEF90HeatXferCtx%fluxLocal,MEF90HeatXferCtx%fluxToIOSF,MEF90HeatXferCtx%IOToFluxSF,MEF90HeatXferCtx%MEF90Ctx%resultViewer,step,ierr))
      Case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90HeatXferCtx%fluxLocal,time,ierr))
      Case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90HeatXferCtx%fluxLocal,1.0_Kr,ierr))
      End Select

      Select case (MEF90HeatXferGlobalOptions%boundaryFluxScaling)
      Case (MEF90Scaling_File)
         PetscCall(MEF90EXOVecLoad(MEF90HeatXferCtx%boundaryFluxLocal,MEF90HeatXferCtx%boundaryFluxToIOSF,MEF90HeatXferCtx%IOToBoundaryFluxSF,MEF90HeatXferCtx%MEF90Ctx%resultViewer,step,ierr))      
      Case (MEF90Scaling_Linear)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90HeatXferCtx%boundaryFluxLocal,time,ierr))
      Case (MEF90Scaling_CST)
         PetscCall(MEF90VecSetValuesFromOptions(MEF90HeatXferCtx%boundaryFluxLocal,1.0_Kr,ierr))
      End Select
   End Subroutine MEF90HeatXferUpdateTransients
   
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferOperator"
!!!
!!!  
!!!  MEF90HeatXferOperator: wraps calls to MEF90HeatXferOperator from m_MEF90_HeatXferAssembly
!!!                         since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90HeatXferOperator(snesTemp,x,residual,MEF90HeatXferCtx,ierr)
      Type(tSNES),Intent(IN)                             :: snesTemp
      Type(tVec),Intent(IN)                              :: x
      Type(tVec),Intent(INOUT)                           :: residual
      Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      PetscInt                                           :: dim      
      PetscCall(DMGetDimension(MEF90HeatXferCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         PetscCall(MEF90HeatXferOperator2D(snesTemp,x,residual,MEF90HeatXferCtx,ierr))
      Else If (dim == 3) Then
         PetscCall(MEF90HeatXferOperator3D(snesTemp,x,residual,MEF90HeatXferCtx,ierr))
      End If      
   End Subroutine MEF90HeatXferOperator
   
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferBilinearForm"
!!!
!!!  
!!!  MEF90HeatXferBilinearForm: wraps calls to MEF90HeatXferBilinearForm from m_MEF90_HeatXferAssembly
!!!                             since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90HeatXferBilinearForm(snesTemp,x,A,M,MEF90HeatXferCtx,ierr)
      Type(tSNES),Intent(IN)                             :: snesTemp
      Type(tVec),Intent(IN)                              :: x
      Type(tMat),Intent(INOUT)                           :: A,M
      Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
      PetscErrorCode,Intent(OUT)                         :: ierr  

      PetscInt                                           :: dim      

      PetscCall(DMGetDimension(MEF90HeatXferCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         PetscCall(MEF90HeatXferBilinearForm2D(snesTemp,x,A,M,MEF90HeatXferCtx,ierr))
      Else If (dim == 3) Then
         PetscCall(MEF90HeatXferBilinearForm3D(snesTemp,x,A,M,MEF90HeatXferCtx,ierr))
      End If      
   End Subroutine MEF90HeatXferBilinearForm

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXFerEnergy"
!!!
!!!  
!!!  MEF90HeatXFerEnergy: wraps calls to MEF90HeatXferEnergy from m_MEF90_HeatXferAssembly
!!!                       since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90HeatXFerEnergy(MEF90HeatXferCtx,energy,bodyWork,surfaceWork,ierr)
      Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
      PetscReal,Dimension(:),Pointer                     :: energy,bodyWork,surfaceWork
      PetscErrorCode,Intent(INOUT)                       :: ierr
   
      PetscInt                                           :: dim      

      PetscCall(DMGetDimension(MEF90HeatXferCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         PetscCall(MEF90HeatXFerEnergy2D(MEF90HeatXferCtx,energy,bodyWork,surfaceWork,ierr))
      Else If (dim == 3) Then
         PetscCall(MEF90HeatXFerEnergy3D(MEF90HeatXferCtx,energy,bodyWork,surfaceWork,ierr))
      End If      
   End Subroutine MEF90HeatXFerEnergy

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXFerIFunction"
!!!
!!!  
!!!  MEF90HeatXFerIFunction: wraps calls to MEF90HeatXFerIFunction from m_MEF90_HeatXferAssembly
!!!                         since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90HeatXFerIFunction(tempTS,time,x,xdot,F,MEF90HeatXferCtx,ierr)
      Type(tTS),Intent(IN)                               :: tempTS
      PetscReal,Intent(IN)                               :: time
      Type(tVec),Intent(IN)                              :: x,xdot
      Type(tVec),Intent(INOUT)                           :: F
      Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
      PetscErrorCode,Intent(OUT)                         :: ierr
      
      PetscInt                                           :: dim      

      PetscCall(DMGetDimension(MEF90HeatXferCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         !PetscCall(MEF90HeatXFerIFunction2D(tempTS,time,x,xdot,F,MEF90HeatXferCtx,ierr))
      Else If (dim == 3) Then
         !PetscCall(MEF90HeatXFerIFunction3D(tempTS,time,x,xdot,F,MEF90HeatXferCtx,ierr))
      End If      
   End Subroutine MEF90HeatXFerIFunction
   
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferIJacobian"
!!!
!!!  
!!!  MEF90HeatXferIJacobian: wraps calls to MEF90HeatXferIJacobian from m_MEF90_HeatXferAssembly
!!!                         since overloading cannot be used here
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90HeatXferIJacobian(tempTS,t,x,xdot,shift,A,M,flg,MEF90HeatXferCtx,ierr)
      Type(tTS),Intent(IN)                               :: tempTS
      PetscReal,Intent(IN)                               :: t
      Type(tVec),Intent(IN)                              :: x,xdot
      PetscReal,Intent(IN)                               :: shift
      Type(tMat),Intent(INOUT)                           :: A,M
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
      PetscErrorCode,Intent(OUT)                         :: ierr  
      
      PetscInt                                           :: dim      

      PetscCall(DMGetDimension(MEF90HeatXferCtx%megaDM,dim,ierr))
      If (dim == 2) Then
         !PetscCall(MEF90HeatXferIJacobian2D(tempTS,t,x,xdot,shift,A,M,flg,MEF90HeatXferCtx,ierr))
      Else If (dim == 3) Then
         !PetscCall(MEF90HeatXferIJacobian3D(tempTS,t,x,xdot,shift,A,M,flg,MEF90HeatXferCtx,ierr))
      End If      
   End Subroutine MEF90HeatXferIJacobian
   
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferViewEXO"
!!!
!!!  
!!!  MEF90HeatXferViewEXO:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90HeatXferViewEXO(MEF90HeatXferCtx,step,ierr)
      Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
      PetscInt,Intent(IN)                                :: step
      PetscErrorCode,Intent(OUT)                         :: ierr

      PetscCall(MEF90EXOVecView(MEF90HeatXferCtx%temperatureLocal,MEF90HeatXferCtx%temperatureToIOSF,MEF90HeatXferCtx%IOToTemperatureSF,MEF90HeatXferCtx%MEF90Ctx%resultViewer,step,ierr))
   End Subroutine MEF90HeatXferViewEXO

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferCreateSNES"
!!!
!!!  
!!!  MEF90HeatXferCreateSNES:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90HeatXferCreateSNES(MEF90HeatXferCtx,snesTemp,residual,ierr)
      Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
      Type(tSNES),Intent(OUT)                            :: snesTemp
      Type(tVec),Intent(IN)                              :: residual
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(MEF90HeatXferGlobalOptions_Type),pointer      :: MEF90HeatXferGlobalOptions
      Type(tDM)                                          :: dm
      Type(tMat)                                         :: matTemp
      Type(tMatNullSpace)                                :: nspTemp
      Type(tKSP)                                         :: kspTemp
      PetscReal                                          :: rtol,dtol

      PetscCall(PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr))
      PetscCall(VecGetDM(MEF90HeatXferCtx%temperatureLocal,dm,ierr))
      PetscCall(DMCreateMatrix(dm,matTemp,iErr))
      PetscCall(MatSetOptionsPrefix(matTemp,"Temperature_",ierr))
      !!! The matrix is not symmetric if the advection vector is /= 0
      PetscCall(MatSetOption(matTemp,MAT_SPD,PETSC_TRUE,ierr))
      PetscCall(MatSetOption(matTemp,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr))
      PetscCall(MatSetOption(matTemp,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr))
      If (MEF90HeatXferGlobalOptions%addNullSpace) Then
         PetscCall(MatNullSpaceCreate(MEF90HeatXferCtx%MEF90Ctx%Comm,PETSC_TRUE,0_Ki,PETSC_NULL_VEC,nspTemp,ierr))
         PetscCall(MatSetNullSpace(matTemp,nspTemp,ierr))
      End If
      PetscCall(MatSetFromOptions(matTemp,ierr))

      PetscCall(SNESCreate(MEF90HeatXferCtx%MEF90Ctx%Comm,snesTemp,ierr))
      PetscCall(SNESSetApplicationContext(snesTemp,MEF90HeatXferCtx,ierr))
      PetscCall(SNESSetDM(snesTemp,dm,ierr))
      PetscCall(SNESSetType(snesTemp,SNESKSPONLY,ierr))
      PetscCall(SNESSetOptionsPrefix(snesTemp,'Temperature_',ierr))

      PetscCall(SNESSetFunction(snesTemp,residual,MEF90HeatXferOperator,MEF90HeatXferCtx,ierr))
      PetscCall(SNESSetJacobian(snesTemp,matTemp,matTemp,MEF90HeatXferBilinearForm,MEF90HeatXferCtx,ierr))
      PetscCall(SNESSetFromOptions(snesTemp,ierr))
      !!! 
      !!! Set some KSP options
      !!!
      PetscCall(SNESGetKSP(snesTemp,kspTemp,ierr))
      PetscCall(KSPSetType(kspTemp,KSPCG,ierr))
      PetscCall(KSPSetInitialGuessNonzero(kspTemp,PETSC_TRUE,ierr))
      rtol = 1.0D-8
      dtol = 1.0D+10
      PetscCall(KSPSetTolerances(kspTemp,rtol,PETSC_DEFAULT_REAL,dtol,PETSC_DEFAULT_INTEGER,ierr))
      PetscCall(KSPSetFromOptions(kspTemp,ierr))
      PetscCall(MatDestroy(matTemp,ierr))
   End Subroutine MEF90HeatXferCreateSNES

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferCreateTS"
!!!
!!!  
!!!  MEF90HeatXferCreateTS:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90HeatXferCreateTS(MEF90HeatXferCtx,tsTemp,residual,initialTime,initialStep,ierr)

      Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
      Type(tTS),Intent(OUT)                              :: tsTemp
      Type(tVec),Intent(IN)                              :: residual
      PetscReal,Intent(IN)                               :: initialTime,initialStep
      PetscErrorCode,Intent(OUT)                         :: ierr

      Type(MEF90HeatXferGlobalOptions_Type),pointer      :: MEF90HeatXferGlobalOptions
      Type(tDM)                                          :: dm
      Type(tMat)                                         :: matTemp
      Type(tMatNullSpace)                                :: nspTemp
      Type(tSNES)                                        :: snesTemp
      Type(tKSP)                                         :: kspTemp
      PetscReal                                          :: rtol,dtol

      PetscCall(PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr))
      PetscCall(VecGetDM(MEF90HeatXferCtx%temperatureLocal,dm,ierr))
      PetscCall(DMCreateMatrix(dm,matTemp,iErr))
      PetscCall(MatSetOptionsPrefix(matTemp,"Temperature_",ierr))
      PetscCall(MatSetOption(matTemp,MAT_SPD,PETSC_TRUE,ierr))
      PetscCall(MatSetOption(matTemp,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr))
      PetscCall(MatSetOption(matTemp,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr))
      If (MEF90HeatXferGlobalOptions%addNullSpace) Then
         PetscCall(MatNullSpaceCreate(MEF90HeatXferCtx%MEF90Ctx%Comm,PETSC_TRUE,0_kI,PETSC_NULL_VEC,nspTemp,ierr))
         PetscCall(MatSetNullSpace(matTemp,nspTemp,ierr))
      End If
      PetscCall(MatSetFromOptions(matTemp,ierr))

      PetscCall(TSCreate(MEF90HeatXferCtx%MEF90Ctx%Comm,tsTemp,ierr))
      PetscCall(TSSetDM(tsTemp,dm,ierr))
      PetscCall(TSSetOptionsPrefix(tsTemp,'Temperature_',ierr))
      PetscCall(TSGetSNES(tsTemp,snesTemp,ierr))

      PetscCall(TSSetIFunction(tsTemp,residual,MEF90HeatXFerIFunction,MEF90HeatXferCtx,ierr))
      PetscCall(TSSetIJacobian(tsTemp,matTemp,matTemp,MEF90HeatXFerIJacobian,MEF90HeatXferCtx,ierr))

      PetscCall(TSSetType(tsTemp,'rosw',ierr))
      PetscCall(TSRosWSetType(tsTemp,'ra3pw',ierr))
      PetscCall(TSSetProblemType(tsTemp,TS_LINEAR,ierr))
      PetscCall(VecSet(MEF90HeatXferCtx%temperatureLocal,MEF90HeatXferGlobalOptions%initialTemperature,ierr))
      PetscCall(TSSetSolution(tsTemp,MEF90HeatXferCtx%temperatureLocal,ierr))
      PetscCall(TSSetInitialTimeStep(tsTemp,initialTime,initialStep,ierr))
      !PetscCall(TSSetExactFinalTime(tsTemp,PETSC_TRUE,ierr))
      PetscCall(TSSetFromOptions(tsTemp,ierr))
      !!! 
      !!! Set some KSP options
      !!!
      PetscCall(SNESGetKSP(snesTemp,kspTemp,ierr))
      PetscCall(KSPSetType(kspTemp,KSPCG,ierr))
      PetscCall(KSPSetInitialGuessNonzero(kspTemp,PETSC_TRUE,ierr))
      rtol = 1.0D-8
      dtol = 1.0D+10
      PetscCall(KSPSetTolerances(kspTemp,rtol,PETSC_DEFAULT_REAL,dtol,PETSC_DEFAULT_INTEGER,ierr))
      PetscCall(KSPSetFromOptions(kspTemp,ierr))
   End Subroutine MEF90HeatXferCreateTS
End Module m_MEF90_HeatXfer
