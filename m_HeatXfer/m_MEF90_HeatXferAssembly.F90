#include "../MEF90/mef90.inc"
Module MEF90_APPEND(m_MEF90_HeatXferAssembly,MEF90_DIM)D
#include "petsc/finclude/petsc.h"
   Use m_MEF90
   Use m_MEF90_HeatXferCtx
   Implicit none
   
   Private
   Public MEF90HeatXferOperator, &
          MEF90HeatXferBilinearForm, &
          MEF90HeatXFerEnergy, &
          MEF90HeatXFerIFunction, &
          MEF90HeatXferIJacobian

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferOperator"
!!!
!!!  
!!!  MEF90HeatXferOperator: Build the operator. When called in SNES, the solution time should always match the target time, 
!!!                         so there is no need for interpolation of the fluxes, external, and boundary values
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90HeatXferOperator(snesTemp,x,residual,MEF90HeatXferCtx,ierr)
      Type(tSNES),Intent(IN)                             :: snesTemp
      Type(tVec),Intent(IN)                              :: x
      Type(tVec),Intent(INOUT)                           :: residual
      Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
      PetscErrorCode,Intent(INOUT)                       :: ierr
   
      ! !Type(DM)                                           :: mesh
      ! Type(SectionReal)                                  :: xSec,residualSec,modifiedFluxSec
      ! Type(SectionReal)                                  :: fluxSec,externalTemperatureSec,boundaryTemperatureSec
      ! PetscReal,Dimension(:),Pointer                     :: xPtr,residualPtr,modifiedFluxPtr
      ! PetscReal,Dimension(:),Pointer                     :: fluxPtr,externalTemperaturePtr,boundaryTemperaturePtr
      ! Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,bcIS
      ! PetscInt,dimension(:),Pointer                      :: setID
      ! PetscInt,Dimension(:),Pointer                      :: setIdx,setdofIdx
      ! PetscInt                                           :: set,QuadratureOrder
      ! Type(MEF90_MATPROP),pointer                        :: matpropSet
      ! Type(MEF90HeatXferCellSetOptions_Type),pointer     :: cellSetOptions
      ! Type(MEF90HeatXferVertexSetOptions_Type),pointer   :: vertexSetOptions
      ! Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elem
      ! Type(MEF90ElementType)                             :: elemType
      ! PetscInt                                           :: cell,dof,nVal
      ! Type(MEF90_VECT)                                   :: advectionVector

      ! !!! Create dof-based sections
      ! Call SectionRealDuplicate(MEF90HeatXferCtx%DMScalSec,xSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealToVec(xSec,MEF90HeatXferCtx%DMScalScatter,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr) 

      ! Call SectionRealDuplicate(MEF90HeatXferCtx%DMScalSec,residualSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealSet(residualSec,0.0_Kr,ierr);CHKERRQ(ierr)
      ! Call VecSet(residual,0.0_Kr,ierr);CHKERRQ(ierr)

      ! Call SectionRealDuplicate(MEF90HeatXferCtx%DMScalSec,boundaryTemperatureSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealToVec(boundaryTemperatureSec,MEF90HeatXferCtx%DMScalScatter,SCATTER_REVERSE,MEF90HeatXferCtx%boundaryTemperature,ierr);CHKERRQ(ierr)

      ! !!! Create cell based sections, and allocate required pointers
      ! Allocate(modifiedFluxPtr(1))
      ! Call SectionRealDuplicate(MEF90HeatXferCtx%cellDMScalSec,modifiedFluxSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealSet(modifiedFluxSec,0.0_Kr,ierr);CHKERRQ(ierr)

      ! Call SectionRealDuplicate(MEF90HeatXferCtx%cellDMScalSec,fluxSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealToVec(fluxSec,MEF90HeatXferCtx%cellDMScalScatter,SCATTER_REVERSE,MEF90HeatXferCtx%flux,ierr);CHKERRQ(ierr) 

      ! Call SectionRealDuplicate(MEF90HeatXferCtx%cellDMScalSec,externalTemperatureSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealToVec(externalTemperatureSec,MEF90HeatXferCtx%cellDMScalScatter,SCATTER_REVERSE,MEF90HeatXferCtx%externalTemperature,ierr);CHKERRQ(ierr) 

      ! Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      ! Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      ! Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)

      ! !!!
      ! !!! We loop over all element twice. The first time in order to assembly all not BC cell sets
      ! !!! In the second pass, we only update the BC where necessary
      ! !!! vertex set BC are updated last, so that they override cell set BC
      ! !!!
      ! Do set = 1,size(setID)
      !    Call PetscBagGetDataMEF90MatProp(MEF90HeatXferCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
      !    Call PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         
      !    If (.not. cellSetOptions%Has_BC) Then
      !       Call DMMeshGetStratumIS(MEF90HeatXferCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
      !       Call ISGetIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
         
      !       elemType = MEF90knownElements(cellSetOptions%ElemTypeShortID)
      !       QuadratureOrder = elemType%order * 2
      !       Call MEF90ElementCreate(MEF90HeatXferCtx%DM,setIS,elem,QuadratureOrder,CellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)
      !       If (elemType%Codim == 0) Then
      !          Call MEF90DiffusionOperatorSet(residualSec,MEF90HeatXferCtx%DM,xSec,setIS,matpropSet%ThermalConductivity,cellSetOptions%SurfaceThermalConductivity,elem,elemType,ierr);CHKERRQ(ierr)
      !          advectionVector = -matPropSet%density * matPropSet%SpecificHeat * cellSetOptions%advectionVector(1:MEF90_DIM)
      !          If (norm(adVectionVector) /= 0.0_Kr) Then
      !             Call MEF90DiffusionOperatorAdvectionSet(residualsec,MEF90HeatXferCtx%DM,xSec,setIS,advectionVector,elem,elemType,ierr);CHKERRQ(ierr)
      !          End If
      !       End If

      !       !!! Modified flux is flux + surfaceThermalConductivity * refTemp      
      !       !!! I _could_ use a SecAXPY, but this would summ all values at all cells for each block
      !       !!! I _could_ also create Sections restricted to the cell set only
      !       !!! But this is going away with sections anyway...
      !       Do cell = 1,size(setIdx)
      !          Call SectionRealRestrict(fluxSec,setIdx(cell),fluxPtr,ierr);CHKERRQ(ierr)
      !          modifiedFluxPtr = - fluxPtr
      !          If (cellSetOptions%SurfaceThermalConductivity /= 0.0_Kr) Then
      !             Call SectionRealRestrict(externalTemperatureSec,setIdx(cell),externalTemperaturePtr,ierr);CHKERRQ(ierr)
      !             modifiedFluxPtr = modifiedFluxPtr - cellSetOptions%SurfaceThermalConductivity * externalTemperaturePtr
      !             Call SectionRealRestore(externalTemperatureSec,setIdx(cell),externalTemperaturePtr,ierr);CHKERRQ(ierr)
      !          End If
      !          Call SectionRealRestore(fluxSec,setIdx(cell),fluxPtr,ierr);CHKERRQ(ierr)
      !          Call SectionRealUpdate(modifiedFluxSec,setIdx(cell),modifiedFluxPtr,ADD_VALUES,ierr);CHKERRQ(ierr)
      !       End Do
      !       Call MEF90DiffusionRHSSetCell(residualsec,MEF90HeatXferCtx%DM,modifiedFluxSec,setIS,elem,elemType,ierr);CHKERRQ(ierr)
          
      !       Call MEF90ElementDestroy(elem,ierr);CHKERRQ(ierr)
      !       Call ISRestoreIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
      !       Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      !    End If ! cellSetOptions%Has_BC
      ! End Do ! set
      
      ! !!! "Ghost update" for the residual Section
      ! Call SectionRealComplete(residualSec,ierr);CHKERRQ(ierr)
      ! !!! Scatter back from SectionReal to Vec
      ! Call SectionRealToVec(residualSec,MEF90HeatXferCtx%DMScalScatter,SCATTER_FORWARD,residual,ierr);CHKERRQ(ierr)

      ! !!!
      ! !!! Cell set BC
      ! !!!
      ! Do set = 1,size(setID)
      !    Call PetscBagGetDataMEF90MatProp(MEF90HeatXferCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
      !    Call PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
      !    If (cellSetOptions%Has_BC) Then
      !       Call DMMeshGetStratumIS(MEF90HeatXferCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
      !       Call MEF90ISCreateCelltoVertex(MEF90HeatXferCtx%DMScal,PETSC_COMM_WORLD,setIS,bcIS,ierr)
      !       Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
      !       Allocate(xPtr(nval),stat=ierr)
      !       Allocate(residualPtr(nval),stat=ierr)
      !       Allocate(boundaryTemperaturePtr(nval),stat=ierr)
      !       Call MEF90VecGetValuesISdof(MEF90HeatXferCtx%DMScal,x,xPtr,bcIS,1,ierr)
      !       Call MEF90VecGetValuesISdof(MEF90HeatXferCtx%DMScal,MEF90HeatXferCtx%boundaryTemperature,boundaryTemperaturePtr,bcIS,1,ierr)
      !       residualPtr = xPtr - boundaryTemperaturePtr
      !       Call MEF90VecSetValuesISdof(MEF90HeatXferCtx%DMScal,residual,residualPtr,bcIS,1,INSERT_VALUES,ierr)
      !       DeAllocate(boundaryTemperaturePtr)
      !       DeAllocate(residualPtr)
      !       DeAllocate(xPtr)
      !       Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
      !       Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      !    End If ! cellSetOptions%Has_BC
      ! End Do ! set
      ! Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      ! Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
            
      ! !!!
      ! !!! Vertex set BC
      ! !!!
      ! Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      ! Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      ! Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      ! Do set = 1,size(setID)
      !    Call PetscBagGetDataMEF90HeatXferCtxVertexSetOptions(MEF90HeatXferCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
      !    If (vertexSetOptions%Has_BC) Then
      !       Call DMMeshGetStratumIS(MEF90HeatXferCtx%DMScal,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
      !       Call ISGetSize(setIS,nval,ierr);CHKERRQ(ierr)
      !       Allocate(xPtr(nval),stat=ierr)
      !       Allocate(residualPtr(nval),stat=ierr)
      !       Allocate(boundaryTemperaturePtr(nval),stat=ierr)
      !       Call MEF90VecGetValuesISdof(MEF90HeatXferCtx%DMScal,x,xPtr,setIS,1,ierr)
      !       Call MEF90VecGetValuesISdof(MEF90HeatXferCtx%DMScal,MEF90HeatXferCtx%boundaryTemperature,boundaryTemperaturePtr,setIS,1,ierr)
      !       residualPtr = xPtr - boundaryTemperaturePtr
      !       Call MEF90VecSetValuesISdof(MEF90HeatXferCtx%DMScal,residual,residualPtr,setIS,1,INSERT_VALUES,ierr)
      !       DeAllocate(boundaryTemperaturePtr)
      !       DeAllocate(residualPtr)
      !       DeAllocate(xPtr)
      !       Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      !    EndIf
      ! End Do ! set
      ! Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      ! Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
   
      ! Call VecAssemblyBegin(residual,ierr)
      ! Call VecAssemblyEnd(residual,ierr)
      
      ! Call SectionRealDestroy(residualSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealDestroy(modifiedfluxSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealDestroy(fluxSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealDestroy(externalTemperatureSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealDestroy(boundaryTemperatureSec,ierr);CHKERRQ(ierr)
      ! DeAllocate(modifiedFluxPtr)
   End Subroutine MEF90HeatXferOperator
   
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferBilinearForm"
!!!
!!!  
!!!  MEF90HeatXferBilinearForm:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90HeatXferBilinearForm(snesTemp,x,A,M,flg,MEF90HeatXferCtx,ierr)
      Type(tSNES),Intent(IN)                             :: snesTemp
      Type(tVec),Intent(IN)                              :: x
      Type(tMat),Intent(INOUT)                           :: A,M
      MatStructure,Intent(INOUT)                         :: flg
      Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
      PetscErrorCode,Intent(INOUT)                       :: ierr  
         
      ! Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof
      ! PetscInt,dimension(:),Pointer                      :: setID
      ! PetscInt                                           :: set,QuadratureOrder
      ! Type(MEF90_MATPROP),pointer                        :: matpropSet
      ! Type(MEF90HeatXferCellSetOptions_Type),pointer     :: cellSetOptions
      ! Type(MEF90HeatXferVertexSetOptions_Type),pointer   :: vertexSetOptions
      ! Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elem
      ! Type(MEF90ElementType)                            :: elemType
      ! PetscInt,Dimension(:),Pointer                      :: setIdx,bcIdx,Cone
      ! Type(IS)                                           :: bcIS
      ! PetscInt                                           :: cell,v,numBC,numDoF,numCell
      ! Type(MEF90_VECT)                                   :: adVectionVector
      
      ! Call MatZeroEntries(A,ierr);CHKERRQ(ierr)
      ! Call DMMeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      ! Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      ! Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      ! Do set = 1,size(setID)
      !    Call PetscBagGetDataMEF90MatProp(MEF90HeatXferCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
      !    Call PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
      !    If (.not. cellSetOptions%Has_BC) Then
      !       elemType = MEF90knownElements(cellSetOptions%ElemTypeShortID)
      !       If (elemType%codim == 0) Then
      !          Call DMMeshGetStratumIS(MEF90HeatXferCtx%DMScal,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
      !          QuadratureOrder = 2 * elemType%order
      !          Call MEF90ElementCreate(MEF90HeatXferCtx%DMScal,setIS,elem,QuadratureOrder,CellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)
      !          Call MEF90DiffusionBilinearFormSet(A,MEF90HeatXferCtx%DMScal,setIS,matpropSet%ThermalConductivity,cellSetOptions%SurfaceThermalConductivity,elem,elemType,ierr);CHKERRQ(ierr)
      !          advectionVector = -matPropSet%density * matPropSet%SpecificHeat * cellSetOptions%advectionVector(1:MEF90_DIM)
      !          If (norm(adVectionVector) /= 0.0_Kr) Then
      !             Call  MEF90DiffusionBilinearFormAdvectionSet(A,MEF90HeatXferCtx%DMScal,setIS,advectionVector,elem,elemType,ierr);CHKERRQ(ierr)
      !          End If
      !          Call MEF90ElementDestroy(elem,ierr)
      !          Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      !       End If
      !    End If ! cellSetOptions%Has_BC
      ! End Do     
      ! Call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
      ! Call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)

      ! !!!
      ! !!! Boundary conditions at cell sets
      ! !!!
      ! Call DMmeshGetStratumSize(MEF90HeatXferCtx%DM,"height",0,numCell,ierr);CHKERRQ(ierr)
      ! Do set = 1,size(setID)
      !    Call PetscBagGetDataMEF90MatProp(MEF90HeatXferCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
      !    Call PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
      !    If (cellSetOptions%Has_BC) Then
      !       Call DMMeshGetStratumIS(MEF90HeatXferCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
      !       Call MEF90ISCreateCelltoVertex(MEF90HeatXferCtx%DMScal,PETSC_COMM_WORLD,setIS,bcIS,ierr)
      !       Call DMMeshISCreateISglobaldof(MEF90HeatXferCtx%DMScal,bcIS,0,setISdof,ierr);CHKERRQ(ierr)
      !       Call MatZeroRowsColumnsIS(A,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
      !       Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
      !       Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      !       Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
      !    End If ! cellSetOptions%Has_BC
      ! End Do     
      ! Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      ! Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
   
      ! !!!
      ! !!! Boundary conditions at vertex sets
      ! !!!
      ! Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      ! Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      ! Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      ! Do set = 1,size(setID)
      !    Call PetscBagGetDataMEF90HeatXferCtxVertexSetOptions(MEF90HeatXferCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
      !    If (vertexSetOptions%Has_BC) Then
      !       Call DMMeshGetStratumIS(MEF90HeatXferCtx%DMScal,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
      !       Call DMMeshISCreateISglobaldof(MEF90HeatXferCtx%DMScal,setIS,0,setISdof,ierr);CHKERRQ(ierr)
      !       Call MatZeroRowsColumnsIS(A,setISdof,123456.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
      !       Call ISDestroy(setISdof,ierr);CHKERRQ(iErr)
      !       Call ISDestroy(setIS,ierr);CHKERRQ(iErr)
      !    End If
      ! End Do
      ! Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      ! Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      ! flg = SAME_NONZERO_PATTERN
   End Subroutine MEF90HeatXferBilinearForm

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXFerEnergy"
!!!
!!!  
!!!  MEF90HeatXFerEnergy:
!!!  
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90HeatXFerEnergy(temperatureVec,t,MEF90HeatXferCtx,energy,work,ierr)
      Type(tVec),Intent(IN)                           :: temperatureVec
      PetscReal,Intent(IN)                            :: t
      Type(MEF90HeatXferCtx_Type),Intent(IN)          :: MEF90HeatXferCtx
      PetscReal,Dimension(:),Pointer                  :: energy,work
      PetscErrorCode,Intent(INOUT)                    :: ierr
   
      ! Type(SectionReal)                               :: temperatureSec,fluxSec,externalTemperatureSec
      ! Type(IS)                                        :: CellSetGlobalIS,setIS,setISdof
      ! PetscInt,dimension(:),Pointer                   :: setID,SetIdx
      ! PetscInt                                        :: set,QuadratureOrder,numCell
      ! Type(MEF90HeatXferCellSetOptions_Type),pointer  :: cellSetOptions
      ! Type(MEF90_MATPROP),pointer                     :: matpropSet
      ! Type(VecScatter)                                :: ScatterSecToVec,ScatterSecToVecCell
      ! Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer   :: elem
      ! Type(MEF90ElementType)                         :: elemType
      ! Type(MEF90CtxGlobalOptions_Type),pointer        :: MEF90CtxGlobalOptions
      ! Type(MEF90HeatXferGlobalOptions_Type),pointer   :: MEF90HeatXferGlobalOptions
      ! PetscReal                                       :: myEnergy,myWork
      
      
      ! energy = 0.0_Kr
      ! work   = 0.0_Kr
      ! Call PetscBagGetDataMEF90CtxGlobalOptions(MEF90HeatXferCtx%MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr);CHKERRQ(ierr)
      ! Call PetscBagGetDataMEF90HeatXferCtxGlobalOptions(MEF90HeatXferCtx%GlobalOptionsBag,MEF90HeatXferGlobalOptions,ierr);CHKERRQ(ierr)
   
      ! !!! Create dof-based sections
      ! Call SectionRealDuplicate(MEF90HeatXferCtx%DMScalSec,temperatureSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealToVec(temperatureSec,MEF90HeatXferCtx%DMScalScatter,SCATTER_REVERSE,temperatureVec,ierr);CHKERRQ(ierr) 

      ! !!! Create cell based sections, and allocate required pointers
      ! Call SectionRealDuplicate(MEF90HeatXferCtx%cellDMScalSec,fluxSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealToVec(fluxSec,MEF90HeatXferCtx%cellDMScalScatter,SCATTER_REVERSE,MEF90HeatXferCtx%flux,ierr);CHKERRQ(ierr) 

      ! Call SectionRealDuplicate(MEF90HeatXferCtx%cellDMScalSec,externalTemperatureSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealToVec(externalTemperatureSec,MEF90HeatXferCtx%cellDMScalScatter,SCATTER_REVERSE,MEF90HeatXferCtx%externalTemperature,ierr);CHKERRQ(ierr) 
      
      ! Call DMMeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      ! Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      ! Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      ! Do set = 1,size(setID)
      !    myenergy = 0.0_Kr
      !    mywork = 0.0_Kr
      !    Call DMMeshGetStratumIS(MEF90HeatXferCtx%DMScal,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
      !    Call PetscBagGetDataMEF90MatProp(MEF90HeatXferCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
      !    Call PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
         
      !    elemType = MEF90knownElements(cellSetOptions%ElemTypeShortID)
      !    QuadratureOrder = elemType%order * 2
      !    Call MEF90ElementCreate(MEF90HeatXferCtx%DMScal,setIS,elem,QuadratureOrder,CellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)   
         
      !    Call MEF90DiffusionEnergySet(myenergy,temperatureSec,externalTemperatureSec,MEF90HeatXferCtx%DMScal,matpropSet%ThermalConductivity,cellSetOptions%SurfaceThermalConductivity,setIS,elem,elemType,ierr);CHKERRQ(ierr)
      !    Call MPI_AllReduce(myEnergy,energy(set),1,MPIU_SCALAR,MPI_SUM,MEF90HeatXferCtx%MEF90Ctx%comm,ierr);CHKERRQ(ierr)
         
      !    If (MEF90HeatXferGlobalOptions%fluxScaling /= MEF90Scaling_Null) Then
      !       Call MEF90DiffusionWorkSetCell(mywork,temperatureSec,MEF90HeatXferCtx%DMScal,fluxSec,setIS,elem,elemType,ierr);CHKERRQ(ierr)
      !       Call MPI_AllReduce(myWork,work(set),1,MPIU_SCALAR,MPI_SUM,MEF90HeatXferCtx%MEF90Ctx%comm,ierr);CHKERRQ(ierr)
      !    End If
      !    Call MEF90ElementDestroy(elem,ierr)
      !    Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      ! End Do
      ! Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      ! Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
      
      ! Call SectionRealDestroy(externalTemperatureSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealDestroy(fluxSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealDestroy(temperatureSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90HeatXFerEnergy

#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXFerIFunction"
!!!
!!!  
!!!  MEF90HeatXFerIFunction:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022    Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90HeatXFerIFunction(tempTS,time,x,xdot,F,MEF90HeatXferCtx,ierr)
      Type(tTS),Intent(IN)                               :: tempTS
      PetscReal,Intent(IN)                               :: time
      Type(tVec),Intent(IN)                              :: x,xdot
      Type(tVec),Intent(INOUT)                           :: F
      Type(MEF90HeatXferCtx_Type),Intent(IN)             :: MEF90HeatXferCtx
      PetscErrorCode,Intent(INOUT)                       :: ierr
         
      ! Type(SNES)                                         :: tempSNES
      ! Type(SectionReal)                                  :: FSec,xdotSec,xSec,modifiedFluxSec
      ! Type(SectionReal)                                  :: fluxSec,externalTemperatureSec,boundaryTemperatureSec
      ! Type(IS)                                           :: CellSetGlobalIS,setIS,bcIS,VertexSetGlobalIS
      ! PetscInt,Dimension(:),Pointer                      :: setID,setIdx
      ! PetscInt                                           :: set,nval
      ! Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elem
      ! Type(MEF90ElementType)                            :: elemType
      ! PetscInt                                           :: QuadratureOrder,cell
      ! Type(MEF90_MATPROP),pointer                        :: matpropSet
      ! Type(MEF90HeatXferCellSetOptions_Type),pointer     :: cellSetOptions
      ! Type(MEF90HeatXferVertexSetOptions_Type),pointer   :: vertexSetOptions
      ! PetscReal,Dimension(:),Pointer                     :: xPtr,FPtr,modifiedFluxPtr
      ! PetscReal,Dimension(:),Pointer                     :: fluxPtr,externalTemperaturePtr,boundaryTemperaturePtr
      ! Type(MEF90_VECT)                                   :: advectionVector
      
      ! Call TSGetSNES(tempTS,tempSNES,ierr);CHKERRQ(ierr)

      ! Call SectionRealDuplicate(MEF90HeatXferCtx%DMScalSec,FSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealSet(Fsec,0.0_Kr,ierr);CHKERRQ(ierr)

      ! Call SectionRealDuplicate(MEF90HeatXferCtx%DMScalSec,xSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealToVec(xSec,MEF90HeatXferCtx%DMScalScatter,SCATTER_REVERSE,x,ierr);CHKERRQ(ierr) 

      ! Call SectionRealDuplicate(MEF90HeatXferCtx%DMScalSec,xdotSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealToVec(xdotSec,MEF90HeatXferCtx%DMScalScatter,SCATTER_REVERSE,xdot,ierr);CHKERRQ(ierr) 

      ! Call SectionRealDuplicate(MEF90HeatXferCtx%DMScalSec,boundaryTemperatureSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealToVec(boundaryTemperatureSec,MEF90HeatXferCtx%DMScalScatter,SCATTER_REVERSE,MEF90HeatXferCtx%boundaryTemperature,ierr);CHKERRQ(ierr) 

      ! !!! Create cell based sections, and allocate required pointers
      ! Allocate(modifiedFluxPtr(1))
      ! Call SectionRealDuplicate(MEF90HeatXferCtx%cellDMScalSec,modifiedFluxSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealSet(modifiedFluxSec,0.0_Kr,ierr);CHKERRQ(ierr)

      ! Call SectionRealDuplicate(MEF90HeatXferCtx%cellDMScalSec,fluxSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealToVec(fluxSec,MEF90HeatXferCtx%cellDMScalScatter,SCATTER_REVERSE,MEF90HeatXferCtx%flux,ierr);CHKERRQ(ierr) 

      ! Call SectionRealDuplicate(MEF90HeatXferCtx%cellDMScalSec,externalTemperatureSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealToVec(externalTemperatureSec,MEF90HeatXferCtx%cellDMScalScatter,SCATTER_REVERSE,MEF90HeatXferCtx%externalTemperature,ierr);CHKERRQ(ierr) 

      ! Call DMMeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      ! Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      ! Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      ! Do set = 1,size(setID) 
      !    Call PetscBagGetDataMEF90MatProp(MEF90HeatXferCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
      !    Call PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)         

      !    If (.not. cellSetOptions%Has_BC) Then
      !       Call DMMeshGetStratumIS(MEF90HeatXferCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
      !       Call ISGetIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)

      !       elemType = MEF90knownElements(cellSetOptions%ElemTypeShortID)
      !       QuadratureOrder = elemType%order * 2
      !       Call MEF90ElementCreate(MEF90HeatXferCtx%DMScal,setIS,elem,QuadratureOrder,CellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)
      !       If (elemType%codim == 0) Then
      !          Call MEF90DiffusionOperatorSet(FSec,MEF90HeatXferCtx%DMScal,xSec,setIS,matpropSet%ThermalConductivity,cellSetOptions%SurfaceThermalConductivity,elem,elemType,ierr);CHKERRQ(ierr)
      !          Call MEF90DiffusionOperatorAddTransientTermSet(FSec,MEF90HeatXferCtx%DMScal,xdotSec,setIS,matpropSet%density*matpropSet%SpecificHeat,elem,elemType,ierr)
      !          advectionVector = -matPropSet%density * matPropSet%SpecificHeat * cellSetOptions%advectionVector(1:MEF90_DIM)
      !          If (norm(advectionVector) /= 0.0_Kr) Then
      !             Call MEF90DiffusionOperatorAdvectionSet(Fsec,MEF90HeatXferCtx%DM,xSec,setIS,advectionVector,elem,elemType,ierr);CHKERRQ(ierr)
      !          End If
      !       End If
      !       !!! Modified flux is flux + surfaceThermalConductivity * refTemp      
      !       !!! I _could_ use a SecAXPY, but this would summ all values at all cells for each block
      !       !!! I _could_ also create Sections restricted to the cell set only
      !       !!! But this is going away with sections anyway...
      !       Do cell = 1,size(setIdx)
      !          Call SectionRealRestrict(fluxSec,setIdx(cell),fluxPtr,ierr);CHKERRQ(ierr)
      !          modifiedFluxPtr = -fluxPtr
      !          Call SectionRealUpdate(modifiedFluxSec,setIdx(cell),modifiedFluxPtr,ADD_VALUES,ierr);CHKERRQ(ierr)
      !          Call SectionRealRestore(fluxSec,setIdx(cell),fluxPtr,ierr);CHKERRQ(ierr)
      !       End Do ! cell
      !       Do cell = 1,size(setIdx)
      !          Call SectionRealRestrict(externalTemperatureSec,setIdx(cell),externalTemperaturePtr,ierr);CHKERRQ(ierr)
      !          modifiedFluxPtr = -cellSetOptions%SurfaceThermalConductivity * externalTemperaturePtr
      !          Call SectionRealUpdate(modifiedFluxSec,setIdx(cell),modifiedFluxPtr,ADD_VALUES,ierr);CHKERRQ(ierr)
      !          Call SectionRealRestore(externalTemperatureSec,setIdx(cell),externalTemperaturePtr,ierr);CHKERRQ(ierr)
      !       End Do ! cell
      !       Call MEF90DiffusionRHSSetCell(FSec,MEF90HeatXferCtx%DMScal,modifiedFluxSec,setIS,elem,elemType,ierr);CHKERRQ(ierr)
      !       Call MEF90ElementDestroy(elem,ierr)
      !       Call ISRestoreIndicesF90(setIS,setIdx,ierr);CHKERRQ(ierr)
      !       Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      !    End If ! cellSetOptions%Has_BC
      ! End Do ! set
      ! DeAllocate(modifiedFluxPtr)

      ! !!! "Ghost update" for the residual Section
      ! Call SectionRealComplete(FSec,ierr);CHKERRQ(ierr)
      ! !!! Scatter back from SectionReal to Vec
      ! Call SectionRealToVec(FSec,MEF90HeatXferCtx%DMScalScatter,SCATTER_FORWARD,F,ierr);CHKERRQ(ierr)

      ! !!!
      ! !!! Cell set BC
      ! !!!
      ! Do set = 1,size(setID)
      !    Call PetscBagGetDataMEF90MatProp(MEF90HeatXferCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
      !    Call PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
      !    If (cellSetOptions%Has_BC) Then
      !       Call DMMeshGetStratumIS(MEF90HeatXferCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
      !       Call MEF90ISCreateCelltoVertex(MEF90HeatXferCtx%DMScal,PETSC_COMM_WORLD,setIS,bcIS,ierr)
      !       Call ISGetSize(bcIS,nval,ierr);CHKERRQ(ierr)
      !       Allocate(xPtr(nval),stat=ierr)
      !       Allocate(FPtr(nval),stat=ierr)
      !       Allocate(boundaryTemperaturePtr(nval),stat=ierr)
      !       Call MEF90VecGetValuesISdof(MEF90HeatXferCtx%DMScal,x,xPtr,bcIS,1,ierr)
      !       Call MEF90VecGetValuesISdof(MEF90HeatXferCtx%DMScal,MEF90HeatXferCtx%boundaryTemperature,boundaryTemperaturePtr,bcIS,1,ierr)
      !       FPtr = xPtr - boundaryTemperaturePtr
      !       Call MEF90VecSetValuesISdof(MEF90HeatXferCtx%DMScal,F,FPtr,bcIS,1,INSERT_VALUES,ierr)
      !       DeAllocate(boundaryTemperaturePtr)
      !       DeAllocate(FPtr)
      !       DeAllocate(xPtr)
      !       Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
      !       Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      !    End If ! cellSetOptions%Has_BC
      ! End Do ! set
      ! Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      ! Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
            
      ! !!!
      ! !!! Vertex set BC
      ! !!!
      ! Call DMMeshGetLabelIdIS(MEF90HeatXferCtx%DMScal,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      ! Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      ! Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      ! Do set = 1,size(setID)
      !    Call PetscBagGetDataMEF90HeatXferCtxVertexSetOptions(MEF90HeatXferCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
      !    If (vertexSetOptions%Has_BC) Then
      !       Call DMMeshGetStratumIS(MEF90HeatXferCtx%DM,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
      !       Call ISGetSize(setIS,nval,ierr);CHKERRQ(ierr)
      !       Allocate(xPtr(nval),stat=ierr)
      !       Allocate(FPtr(nval),stat=ierr)
      !       Allocate(boundaryTemperaturePtr(nval),stat=ierr)
      !       Call MEF90VecGetValuesISdof(MEF90HeatXferCtx%DMScal,x,xPtr,setIS,1,ierr)
      !       Call MEF90VecGetValuesISdof(MEF90HeatXferCtx%DMScal,MEF90HeatXferCtx%boundaryTemperature,boundaryTemperaturePtr,setIS,1,ierr)
      !       FPtr = xPtr - boundaryTemperaturePtr
      !       Call MEF90VecSetValuesISdof(MEF90HeatXferCtx%DMScal,F,FPtr,setIS,1,INSERT_VALUES,ierr)
      !       DeAllocate(boundaryTemperaturePtr)
      !       DeAllocate(Fptr)
      !       DeAllocate(xPtr)
      !       Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      !    EndIf
      ! End Do ! set
      ! Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      ! Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
           
      ! Call VecAssemblyBegin(F,ierr);CHKERRQ(ierr)
      ! Call VecAssemblyEnd(F,ierr);CHKERRQ(ierr)
      ! Call SectionRealDestroy(FSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealDestroy(xdotSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealDestroy(xSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealDestroy(modifiedFluxSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealDestroy(fluxSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealDestroy(externalTemperatureSec,ierr);CHKERRQ(ierr)
      ! Call SectionRealDestroy(boundaryTemperatureSec,ierr);CHKERRQ(ierr)
   End Subroutine MEF90HeatXFerIFunction
   
#undef __FUNCT__
#define __FUNCT__ "MEF90HeatXferIJacobian"
!!!
!!!  
!!!  MEF90HeatXferIJacobian:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
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
      PetscErrorCode,Intent(INOUT)                       :: ierr  

      ! Type(SNES)                                         :: tempSNES
      ! Type(IS)                                           :: VertexSetGlobalIS,CellSetGlobalIS,setIS,setISdof
      ! PetscInt,dimension(:),Pointer                      :: setID
      ! PetscInt                                           :: set,QuadratureOrder
      ! Type(MEF90_MATPROP),pointer                        :: matpropSet
      ! Type(MEF90HeatXferCellSetOptions_Type),pointer     :: cellSetOptions
      ! Type(MEF90HeatXferVertexSetOptions_Type),pointer   :: vertexSetOptions
      ! !Type(DM)                                           :: mesh
      ! Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elem
      ! Type(MEF90ElementType)                            :: elemType
      ! PetscInt,Dimension(:),Pointer                      :: setIdx,bcIdx,Cone
      ! Type(IS)                                           :: bcIS
      ! PetscInt                                           :: cell,v,numBC,numDoF,numCell
      ! Type(MEF90_VECT)                                   :: adVectionVector

      ! Call MatZeroEntries(A,ierr);CHKERRQ(ierr)
      ! Call TSGetSNES(tempTS,tempSNES,ierr);CHKERRQ(ierr)

      ! Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DMScal,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
      ! Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr) 
      ! Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      ! Do set = 1,size(setID)
      !    Call PetscBagGetDataMEF90MatProp(MEF90HeatXferCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
      !    Call PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
      !    If (.not. cellSetOptions%Has_BC) Then
      !       elemType = MEF90knownElements(cellSetOptions%ElemTypeShortID)
      !       If (elemType%codim == 0) Then
      !          Call DMMeshGetStratumIS(MEF90HeatXferCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
      !          QuadratureOrder = 2 * elemType%order
      !          Call MEF90ElementCreate(MEF90HeatXferCtx%DMScal,setIS,elem,QuadratureOrder,CellSetOptions%ElemTypeShortID,ierr);CHKERRQ(ierr)
      !          Call MEF90DiffusionBilinearFormSet(A,MEF90HeatXferCtx%DMScal,setIS,matpropSet%ThermalConductivity,cellSetOptions%SurfaceThermalConductivity,elem,elemType,ierr);CHKERRQ(ierr)
      !          Call MEF90_MassMatrixAssembleSet(A,MEF90HeatXferCtx%DMScal,setIS,shift*matpropSet%density*matpropSet%SpecificHeat,elem,elemType,ierr)
      !          advectionVector = -matPropSet%density * matPropSet%SpecificHeat * cellSetOptions%advectionVector(1:MEF90_DIM)
      !          If (norm(adVectionVector) /= 0.0_Kr) Then
      !             Call  MEF90DiffusionBilinearFormAdvectionSet(A,MEF90HeatXferCtx%DMScal,setIS,advectionVector,elem,elemType,ierr);CHKERRQ(ierr)
      !          End If
      !          Call MEF90ElementDestroy(elem,ierr)
      !          Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      !       End If
      !    End If ! cellSetOptions%Has_BC
      ! End Do     
      ! Call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
      ! Call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,iErr);CHKERRQ(iErr)
      ! !!!
      ! !!! Boundary conditions at cell sets
      ! !!!
      ! Call DMmeshGetStratumSize(MEF90HeatXferCtx%DMScal,"height",0,numCell,ierr);CHKERRQ(ierr)
      ! Do set = 1,size(setID)
      !    Call PetscBagGetDataMEF90MatProp(MEF90HeatXferCtx%MaterialPropertiesBag(set),matpropSet,ierr);CHKERRQ(ierr)
      !    Call PetscBagGetDataMEF90HeatXferCtxCellSetOptions(MEF90HeatXferCtx%CellSetOptionsBag(set),cellSetOptions,ierr);CHKERRQ(ierr)
      !    If (cellSetOptions%Has_BC) Then
      !       Call DMMeshGetStratumIS(MEF90HeatXferCtx%DM,'Cell Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
      !       Call MEF90ISCreateCelltoVertex(MEF90HeatXferCtx%DMScal,PETSC_COMM_WORLD,setIS,bcIS,ierr)
      !       Call DMMeshISCreateISglobaldof(MEF90HeatXferCtx%DMScal,bcIS,0,setISdof,ierr);CHKERRQ(ierr)
      !       Call MatZeroRowsColumnsIS(A,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
      !       Call ISDestroy(bcIS,ierr);CHKERRQ(ierr)
      !       Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      !       Call ISDestroy(setISdof,ierr);CHKERRQ(ierr)
      !    End If ! cellSetOptions%Has_BC
      ! End Do     
      ! Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      ! Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
   
      ! !!!
      ! !!! Boundary conditions at vertex sets
      ! !!!
      ! Call DMmeshGetLabelIdIS(MEF90HeatXferCtx%DM,'Vertex Sets',VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      ! Call MEF90ISAllGatherMerge(PETSC_COMM_WORLD,VertexSetGlobalIS,ierr);CHKERRQ(ierr) 
      ! Call ISGetIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      ! Do set = 1,size(setID)
      !    Call PetscBagGetDataMEF90HeatXferCtxVertexSetOptions(MEF90HeatXferCtx%VertexSetOptionsBag(set),vertexSetOptions,ierr);CHKERRQ(ierr)
      !    If (vertexSetOptions%Has_BC) Then
      !       Call DMMeshGetStratumIS(MEF90HeatXferCtx%DM,'Vertex Sets',setID(set),setIS,ierr);CHKERRQ(iErr)
      !       Call DMMeshISCreateISglobaldof(MEF90HeatXferCtx%DM,setIS,0,setISdof,ierr);CHKERRQ(ierr)
      !       Call MatZeroRowsColumnsIS(A,setISdof,1.0_Kr,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
      !    End If
      ! End Do
      ! Call ISRestoreIndicesF90(VertexSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      ! Call ISDestroy(VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      
      ! flg = SAME_NONZERO_PATTERN
   End Subroutine MEF90HeatXferIJacobian
End Module MEF90_APPEND(m_MEF90_HeatXferAssembly,MEF90_DIM)D
