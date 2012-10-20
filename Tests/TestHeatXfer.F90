Program TestHeatXfer
#include <finclude/petscdef.h>
   Use m_MEF90
   Use m_MEF90_HeatXfer
   Use petsc
   Implicit NONE   

   PetscErrorCode                                     :: ierr
   Type(MEF90HeatXferCtx_Type)                        :: MEF90HeatXferCtx
   Type(MEF90HeatXferGlobalOptions_Type),Parameter    :: MEF90HeatXferDefaultGlobalOptions = MEF90HeatXferGlobalOptions_Type( &
                                                         MEF90HeatXFer_ModeSteadyState, & ! mode
                                                         PETSC_FALSE,   & ! addNullSpace
                                                         1,             & ! tempOffset
                                                         PETSC_FALSE,   & ! boundaryTempCst
                                                         2,             & ! boundaryTempOffset
                                                         PETSC_TRUE,    & ! externalTempCst
                                                         3,             & ! externalTempOffset
                                                         PETSC_TRUE,    & ! fluxCst
                                                         4)               ! fluxOffset
   Type(MEF90HeatXferCellSetOptions_Type),Parameter   :: MEF90HeatXferDefaultCellSetOptions = MEF90HeatXferCellSetOptions_Type( &
                                                         -1,            & ! elemTypeShortID will be overriden
                                                         0.0_Kr,        & ! flux
                                                         0.0_Kr,        & ! surfaceThermalConductivity
                                                         0.0_Kr)          ! externalTemp
   Type(MEF90HeatXferVertexSetOptions_Type),Parameter :: MEF90HeatXferDefaultVertexSetOptions = MEF90HeatXferVertexSetOptions_Type( &
                                                         PETSC_TRUE,    & ! Has BC
                                                         0.0_Kr)          ! boundaryTemp
                                                         
   Type(MEF90Ctx_Type),target                         :: MEF90Ctx
   Type(MEF90CtxGlobalOptions_Type)                   :: MEF90GlobalOptions_default
   Type(DM),target                                    :: Mesh
   PetscBool                                          :: flg
   Character(len=MEF90_MXSTRLEN)                      :: IOBuffer
   

   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   Call MEF90_Initialize(ierr)

   MEF90GlobalOptions_default%verbose           = 1
   MEF90GlobalOptions_default%timeInterpolation = MEF90TimeInterpolation_linear
   MEF90GlobalOptions_default%timeMin           = 0.0_Kr
   MEF90GlobalOptions_default%timeMax           = 1.0_Kr
   MEF90GlobalOptions_default%timeNumStep       = 11
   MEF90GlobalOptions_default%fileFormat        = MEF90FileFormat_EXOSingle
   
   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   Call MEF90_Initialize(ierr)

   Call MEF90Ctx_Create(PETSC_COMM_WORLD,MEF90Ctx,MEF90GlobalOptions_default,ierr);CHKERRQ(ierr)
   Call MEF90Ctx_GetDMMeshEXO(MEF90Ctx,Mesh,ierr);CHKERRQ(ierr)
   Call DMView(Mesh,PETSC_VIEWER_STDOUT_WORLD,ierr)
   
   Call MEF90Ctx_OpenEXO(MEF90Ctx,Mesh,ierr)

   Call MEF90HeatXferCtx_Create(MEF90HeatXferCtx,Mesh,MEF90Ctx,ierr);CHKERRQ(ierr)
   Call MEF90HeatXferCtx_SetFromOptions(MEF90HeatXferCtx,PETSC_NULL_CHARACTER,Mesh,MEF90Ctx,MEF90HeatXferDefaultGlobalOptions, &
                                        MEF90HeatXferDefaultCellSetOptions,MEF90HeatXferDefaultVertexSetOptions,ierr)
   
   
   Call MEF90HeatXferCtx_Destroy(MEF90HeatXferCtx,ierr);CHKERRQ(ierr)
   
   Call MEF90Ctx_CloseEXO(MEF90Ctx,ierr)
   Call MEF90_Finalize(ierr)
   Call PetscFinalize()
End Program TestHeatXfer