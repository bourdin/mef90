Program TestQuadrature
#include "finclude/petscdef.h"
   Use m_MEF90
   Use m_MEF90_DefMechCtx
   Use m_MEF90_DefMech
   Use petsc
   Implicit NONE   

   PetscErrorCode                      :: ierr
   Type(DM),target                     :: Mesh
   Character(len=MEF90_MXSTRLEN)       :: IOBuffer
   PetscInt                            :: dim
   Type(Vec)                           :: xVec
   PetscReal,Dimension(:,:),Pointer    :: coordPtr
   PetscReal,Dimension(:),Pointer      :: xPtr
   PetscBool                           :: flg
   Integer                             :: i,j,k,v

   PetscReal                           :: scal,sol
   Type(Vect2D)                        :: v2D
   Type(Vect3D)                        :: v3D
   Type(MatS2D)                        :: M2D
   Type(MatS3D)                        :: M3D

   Integer                             :: numVertex
   Integer                             :: QuadOrderMax,QuadOrder
   PetscInt                            :: localSize,globalSize
   Type(SectionReal)                   :: defaultSection


   Type(MEF90Ctx_Type),target                         :: MEF90Ctx
   Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions
   Type(MEF90CtxGlobalOptions_Type),Parameter         :: MEF90DefaultGlobalOptions = MEF90CtxGlobalOptions_Type( &
                                                         1,                             & ! verbose
                                                         PETSC_FALSE,                   & ! validate
                                                         MEF90TimeInterpolation_linear, & ! timeInterpolation
                                                         0.0_Kr,                        & ! timeMin
                                                         1.0_Kr,                        & ! timeMax
                                                         11,                            & ! timeNumStep
                                                         MEF90FileFormat_EXOSingle)       ! fileFormat

   !!! Defect mechanics contexts
   Type(MEF90DefMechCtx_Type)                         :: MEF90DefMechCtx
   Type(MEF90DefMechGlobalOptions_Type),pointer       :: MEF90DefMechGlobalOptions
   Type(MEF90DefMechGlobalOptions_Type),Parameter     :: MEF90DefMechDefaultGlobalOptions2D = MEF90DefMechGlobalOptions_Type( &
                                                         MEF90DefMech_ModeQuasiStatic, & ! mode
                                                         PETSC_TRUE,              & ! disp_addNullSpace
                                                         3,                       & ! DisplacementOffset
                                                         2,                       & ! DamageOffset
                                                         3,                       & ! boundaryDisplacementOffset
                                                         0,                       & ! boundaryDamageOffset
                                                         1,                       & ! temperatureOffset
                                                         4,                       & ! ForceOffset
                                                         3,                       & ! pressureForceOffset
                                                         0,                       & ! plasticStrainOffset
                                                         6,                       & ! StressOffset
                                                         MEF90Scaling_Linear,     & ! boundaryDisplacementScaling
                                                         MEF90Scaling_CST,        & ! boundaryDamageScaling
                                                         MEF90Scaling_Linear,     & ! ForceScaling
                                                         MEF90Scaling_Linear,     & ! pressureForceScaling
                                                         1e-4,                    & ! damage_atol
                                                         1000,                    & ! maxit
                                                         0.,                      & ! irrevThres 
                                                         MEF90DefMech_BTTypeNULL, & ! BTType
                                                         -1,                      & ! BTInt
                                                         -1,                      & ! BTScope
                                                         1.0e-2)                    ! BTTol
   Type(MEF90DefMechGlobalOptions_Type),Parameter     :: MEF90DefMechDefaultGlobalOptions3D = MEF90DefMechGlobalOptions_Type( &
                                                         MEF90DefMech_ModeQuasiStatic, & ! mode
                                                         PETSC_TRUE,              & ! disp_addNullSpace
                                                         3,                       & ! DisplacementOffset
                                                         2,                       & ! DamageOffset
                                                         3,                       & ! boundaryDisplacementOffset
                                                         0,                       & ! boundaryDamageOffset
                                                         1,                       & ! temperatureOffset
                                                         4,                       & ! ForceOffset
                                                         3,                       & ! pressureForceOffset
                                                         0,                       & ! plasticStrainOffset
                                                         7,                       & ! StressOffset
                                                         MEF90Scaling_Linear,     & ! boundaryDisplacementScaling
                                                         MEF90Scaling_CST,        & ! boundaryDamageScaling
                                                         MEF90Scaling_Linear,     & ! ForceScaling
                                                         MEF90Scaling_Linear,     & ! pressureForceScaling
                                                         1e-4,                    & ! damage_atol
                                                         1000,                    & ! maxit
                                                         0.,                      & ! irrevThres 
                                                         MEF90DefMech_BTTypeNULL, & ! BTType
                                                         -1,                      & ! BTInt
                                                         -1,                      & ! BTScope
                                                         1.0e-2)                    ! BTTol

   Type(MEF90DefMechCellSetOptions_Type),Parameter    :: MEF90DefMechDefaultCellSetOptions = MEF90DefMechCellSetOptions_Type( &
                                                         -1,                                      & ! elemTypeShortIDDispl will be overriden
                                                         -1,                                      & ! elemTypeShortIDDamage will be overriden
                                                         [0.0_Kr,0.0_Kr,0.0_Kr],                  & ! force
                                                         0.0_Kr,                                  & ! pressureForce
                                                         MEF90DefMech_damageTypeAT1,              & ! damageType
                                                         MEF90DefMech_plasticityTypeNone,         & ! plasticityType
                                                         MEF90DefMech_unilateralContactTypeNone,  & ! unilateralContactType
                                                         [PETSC_FALSE,PETSC_FALSE,PETSC_FALSE],   & ! Has Displacement BC
                                                         [0.0_Kr,0.0_Kr,0.0_Kr],                  & ! boundary Displacement
                                                         PETSC_FALSE,                             & ! Has Damage BC
                                                         0._Kr)                                     ! Boundary Damage
   Type(MEF90DefMechVertexSetOptions_Type),Parameter  :: MEF90DefMechDefaultVertexSetOptions = MEF90DefMechVertexSetOptions_Type( &
                                                         [PETSC_FALSE,PETSC_FALSE,PETSC_FALSE],   & ! Has Displacement BC
                                                         [0.0_Kr,0.0_Kr,0.0_Kr],                  & ! boundary Displacement
                                                         PETSC_FALSE,                             & ! Has Damage BC
                                                         0.0_Kr)                                    ! boundary Damage

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

   !!! Create sections, vectors, and solvers for DefMech Context
   Call MEF90DefMechCtxSetSections(MEF90DefMechCtx,ierr)
   Call MEF90DefMechCtxCreateVectors(MEF90DefMechCtx,ierr)

   Call DMMeshGetStratumSize(mesh,"depth",0,numVertex,ierr);CHKERRQ(ierr)
   QuadOrderMax = 4
   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-order',QuadOrderMax,flg,ierr);CHKERRQ(ierr);
   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-i',i,flg,ierr);CHKERRQ(ierr);
   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-j',j,flg,ierr);CHKERRQ(ierr);
   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-k',k,flg,ierr);CHKERRQ(ierr);
   Do QuadOrder = QuadOrderMax, QuadOrderMax
      Do k = 0, (dim-2) * QuadOrderMax
         Do j = 0, QuadOrderMax
            Do i = 0, QuadOrderMax
               !!! Initialize a field
               If (i+j+k <= QuadOrder) Then
                  sol = 1.0_Kr / (1.0_kr + i) / (1.0_Kr + j) / (1.0_Kr + k)
                  !!! Integrate
                  If (dim ==2) Then
                     Call Integrate2D_Scal(MEF90DefMechCtx,i,j,QuadOrder,Scal,v2d,ierr)
                  Else
                     Call Integrate3D_Scal(MEF90DefMechCtx,i,j,k,QuadOrder,Scal,v3d,ierr)
                  End If
                  Write(IOBuffer,100) i,j,k,QuadOrder,Scal, sol, (Scal - sol)/Scal
                  Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr) 
               End If
            End Do
         End Do
      End Do
   End Do
!   Do QuadOrder = 0, QuadOrderMax
!      If (dim ==2) Then
!         Call Integrate2D_Scal(MEF90DefMechCtx,i,j,QuadOrder,Scal,v2d,ierr)
!      Else
!         Call Integrate3D_Scal(MEF90DefMechCtx,i,j,k,QuadOrder,Scal,v3d,ierr)
!      End If
!      Write(IOBuffer,200) i,j,k,QuadOrder,Scal
!      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)   
!   End Do
   
   100 Format('Integrating x^',I1,' * Y^',I1,' * Z^', I1,' at order ',I4,' : ',2ES12.5,': relative error ',ES12.5,"\n")
   200 Format('Integrating x^',I1,' * Y^',I1,' * Z^', I1,' at order ',I4,' : ',ES12.5,'\n')

   
   Call MEF90DefMechCtxDestroy(MEF90DefMechCtx,ierr);CHKERRQ(ierr)
   Call MEF90CtxDestroy(MEF90Ctx,ierr);CHKERRQ(ierr)   
   Call MEF90Finalize(ierr)
   Call PetscFinalize()
   
Contains
   Subroutine Integrate3D_Scal(MEF90DefMechCtx,i,j,k,QuadratureOrder,i1,i2,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscInt                                           :: i,j,k,QuadratureOrder
      PetscReal,Intent(OUT)                              :: i1
      Type(Vect3D),Intent(OUT)                           :: i2
      PetscErrorCode,Intent(OUT)                         :: ierr      
                     
      Integer                                            :: iDof,numDof
      Integer                                            :: iGauss,numGauss
      Type(SectionReal)                                  :: coordSec
      PetscReal,Dimension(:),Pointer                     :: coordDof
      PetscReal                                          :: X,Y,Z
      
      Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(IS)                                           :: CellSetGlobalIS,setIS
      PetscInt,Dimension(:),Pointer                      :: setID,cellID
      PetscInt                                           :: set,cell
      Type(MEF90Element_Type)                            :: ElemType
      Type(MEF90Element3D_Scal),Dimension(:),Pointer     :: Elem

      i1 = 0.0_Kr
      i2 = 0.0_Kr
      
      Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   
      Call DMMeshGetSectionReal(mesh,'coordinates',coordSec,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)

         ElemType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)
         Call MEF90Element_Create(mesh,setIS,Elem,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
         !!! Integrate
         numDof   = ElemType%numDof
         numGauss = size(Elem(1)%BF,2)
         Allocate(coordDof(numDof*3))
         Do cell = 1,size(cellID)
            !!! Get value of each field at each Dof of the local element
            Call SectionRealRestrictClosure(coordSec,MEF90DefMechCtx%DMVect,cellID(cell),numDof*3,coordDof,ierr);CHKERRQ(ierr)
            Do iGauss = 1, numGauss
               X = 0.0_Kr
               Y = 0.0_Kr
               Z = 0.0_Kr
               Do iDof = 1, numDof
                  X = X + Elem(cell)%BF(iDoF,iGauss) * coordDof(3*(iDof-1)+1)
                  Y = Y + Elem(cell)%BF(iDoF,iGauss) * coordDof(3*(iDof-1)+2)
                  Z = Z + Elem(cell)%BF(iDoF,iGauss) * coordDof(3*(iDof-1)+3)
               End Do
               i1 = i1 + Elem(cell)%Gauss_C(iGauss) * X**i * Y**j * Z**k
            End Do
         End Do
         Call MEF90Element_Destroy(elem,ierr)
         DeAllocate(coordDof)
         Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)

      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(coordSec,ierr);CHKERRQ(ierr)
   End Subroutine Integrate3D_Scal

   Subroutine Integrate2D_Scal(MEF90DefMechCtx,i,j,QuadratureOrder,i1,i2,ierr)
      Type(MEF90DefMechCtx_Type),Intent(IN)              :: MEF90DefMechCtx
      PetscInt                                           :: i,j,QuadratureOrder
      PetscReal,Intent(OUT)                              :: i1
      Type(Vect2D),Intent(OUT)                           :: i2
      PetscErrorCode,Intent(OUT)                         :: ierr      
                     
      Integer                                            :: iDof,numDof
      Integer                                            :: iGauss,numGauss
      Type(SectionReal)                                  :: coordSec
      PetscReal,Dimension(:),Pointer                     :: coordDof
      PetscReal                                          :: X,Y
      
      Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions
      Type(MEF90DefMechCellSetOptions_Type),Pointer      :: cellSetOptions
      Type(IS)                                           :: CellSetGlobalIS,setIS
      PetscInt,Dimension(:),Pointer                      :: setID,cellID
      PetscInt                                           :: set,cell
      Type(MEF90Element_Type)                            :: ElemType
      Type(MEF90Element2D_Scal),Dimension(:),Pointer     :: Elem

      i1 = 0.0_Kr
      i2 = 0.0_Kr
      
      Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   
      Call DMMeshGetSectionReal(mesh,'coordinates',coordSec,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call PetscBagGetDataMEF90DefMechCtxCellSetOptions(MEF90DefMechCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         Call DMMeshGetStratumIS(MEF90DefMechCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
         Call ISGetIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)

         ElemType = MEF90KnownElements(cellSetOptions%elemTypeShortIDDamage)
         Call MEF90Element_Create(mesh,setIS,Elem,QuadratureOrder,CellSetOptions%elemTypeShortIDDamage,ierr);CHKERRQ(ierr)
         !!! Integrate
         numDof   = ElemType%numDof
         numGauss = size(Elem(1)%BF,2)
         Allocate(coordDof(numDof*2))
         Do cell = 1,size(cellID)
            !!! Get value of each field at each Dof of the local element
            Call SectionRealRestrictClosure(coordSec,MEF90DefMechCtx%DMVect,cellID(cell),numDof*2,coordDof,ierr);CHKERRQ(ierr)
            Do iGauss = 1, numGauss
               X = 0.0_Kr
               Y = 0.0_Kr
               Do iDof = 1, numDof
                  X = X + Elem(cell)%BF(iDoF,iGauss) * coordDof(2*(iDof-1)+1)
                  Y = Y + Elem(cell)%BF(iDoF,iGauss) * coordDof(2*(iDof-1)+2)
               End Do
               i1 = i1 + Elem(cell)%Gauss_C(iGauss) * X**i * Y**j
            End Do
         End Do
         Call MEF90Element_Destroy(elem,ierr)
         DeAllocate(coordDof)
         Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(setIS,cellID,ierr);CHKERRQ(ierr)

      Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call SectionRealDestroy(coordSec,ierr);CHKERRQ(ierr)
   End Subroutine Integrate2D_Scal
End Program TestQuadrature
