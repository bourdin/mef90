#include "../MEF90/mef90.inc"
Program TestQuadrature
#include "petsc/finclude/petsc.h"
   Use m_MEF90
   Use petsc
   Implicit NONE   

   PetscErrorCode                      :: ierr
   Type(tDM),target                    :: dm,dmU
   Type(tPetscSection),target          :: sectionU
   Character(len=MEF90MXSTRLEN)        :: IOBuffer
   PetscInt                            :: dim
   PetscBool                           :: flg
   PetscInt                            :: i,j,k

   PetscReal                           :: scal,sol,xr=1.0_Kr,xl=0.0_Kr,yr=1.0_Kr,yl=0.0_Kr,zr=1.0_Kr,zl=0.0_Kr
   Character(len=MEF90MXSTRLEN)        :: name
   Type(tVec)                          :: locVecU
   Type(Vect2D)                        :: v2D
   Type(Vect3D)                        :: v3D
   PetscInt                            :: QuadOrderMax,QuadOrder


   Type(MEF90Ctx_Type),target                         :: MEF90Ctx
   Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90GlobalOptions
   Type(MEF90CtxGlobalOptions_Type)                   :: MEF90GlobalOptions_default

   !!! Initialize MEF90
   PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
   PetscCallA(MEF90Initialize(ierr))

   MEF90GlobalOptions_default%verbose           = 1
   MEF90GlobalOptions_default%dryrun            = PETSC_FALSE
   MEF90GlobalOptions_default%timeMin           = 0.0_Kr
   MEF90GlobalOptions_default%timeMax           = 1.0_Kr
   MEF90GlobalOptions_default%timeNumStep       = 1
   MEF90GlobalOptions_default%timeSkip          = 0
   MEF90GlobalOptions_default%timeNumCycle      = 1
   MEF90GlobalOptions_default%timeInterpolation = MEF90TimeInterpolation_linear
   MEF90GlobalOptions_default%elementFamily     = MEF90ElementFamilyLagrange
   MEF90GlobalOptions_default%elementOrder      = 1

   !!! Get all MEF90-wide options
   PetscCallA(MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90GlobalOptions_default,ierr))
   PetscCallA(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90GlobalOptions,ierr))

   PetscCallA(DMPlexCreateFromFile(MEF90Ctx%Comm,MEF90Ctx%geometryFile,PETSC_NULL_CHARACTER,PETSC_TRUE,dm,ierr))
   PetscCallA(DMSetFromOptions(dm,ierr))
   PetscCallA(DMViewFromOptions(dm,PETSC_NULL_OPTIONS,"-dm_view",ierr))
   PetscCall(DMGetDimension(dm,dim,ierr))

   QuadOrderMax = 4
   PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-order",QuadOrderMax,flg,ierr))
   PetscCallA(PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-xl",xl,flg,ierr))
   PetscCallA(PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-xr",xr,flg,ierr))
   PetscCallA(PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-yl",yl,flg,ierr))
   PetscCallA(PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-yr",yr,flg,ierr))
   PetscCallA(PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-zl",zl,flg,ierr))
   PetscCallA(PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-zr",zr,flg,ierr))

   ! Create nodal local Vec holding coordinates
   name = "U"
   PetscCallA(MEF90CreateLocalVector(dm,MEF90GlobalOptions%elementFamily,MEF90GlobalOptions%elementOrder,dim,name,locVecU,ierr))
   PetscCallA(VecGetDM(locVecU,dmU,ierr))
   PetscCallA(DMGetLocalSection(dmU,sectionU,ierr))
   PetscCall(project(locVecU,sectionU,ierr))

   Do QuadOrder = QuadOrderMax, QuadOrderMax
      Do k = 0, (dim-2) * QuadOrderMax
         Do j = 0, QuadOrderMax
            Do i = 0, QuadOrderMax
               !!! Initialize a field
               If (i+j+k <= QuadOrder) Then
                  !!! Integrate
                  If (dim ==2) Then
                     sol = ((xr**(1.0_kr + i) - xl**(1.0_kr + i)) / (1.0_kr + i)) * ((yr**(1.0_kr + j) - yl**(1.0_kr + j)) / (1.0_kr + j))
                     PetscCallA(Integrate2D_Scal(MEF90Ctx,locVecU,i,j,QuadOrder,Scal,v2d,ierr))
                  Else
                     sol = ((xr**(1.0_kr + i) - xl**(1.0_kr + i)) / (1.0_kr + i)) * ((yr**(1.0_kr + j) - yl**(1.0_kr + j)) / (1.0_kr + j)) * ((zr**(1.0_kr + k) - zl**(1.0_kr + k)) / (1.0_kr + k))
                     PetscCallA(Integrate3D_Scal(MEF90Ctx,locVecU,i,j,k,QuadOrder,Scal,v3d,ierr))
                  End If
                  Write(IOBuffer,100) i,j,k,QuadOrder,Scal, sol, abs(Scal - sol),  abs(Scal - sol)/sol
                  PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
               End If
            End Do
         End Do
      End Do
   End Do
   
   100 Format('Integrating X^',I1,' * Y^',I1,' * Z^', I1,' at order ',I4,' : ',2ES12.5,': absolute error ',ES12.5,': relative error ',ES12.5,"\n")
   
   PetscCallA(VecDestroy(locVecU,ierr))
   PetscCall(DMDestroy(dm,ierr))
   PetscCall(MEF90CtxDestroy(MEF90Ctx,ierr))
   PetscCall(MEF90Finalize(ierr))
   PetscCall(PetscFinalize(ierr))
   
Contains
   Subroutine project(v,s,ierr)
      Type(tVec),intent(IN)              :: v
      Type(tPetscSection),intent(IN)     :: s
      PetscErrorCode,intent(INOUT)       :: ierr

      PetscInt                           :: pStart,pEnd,p,numDof,i
      Type(tDM)                          :: dm
      Type(tPetscSection)                :: coordSection
      Type(tVec)                         :: coordVec
      PetscScalar,dimension(:),Pointer   :: coordArray,vArray
      PetscScalar,dimension(3)           :: xyz
      PetscInt                           :: dim,pOffset

      PetscCallA(PetscSectionGetChart(s,pStart,pEnd,ierr))
      PetscCallA(VecGetDM(v,dm,ierr))
      PetscCallA(DMGetCoordinateSection(dm,coordSection,ierr))
      PetscCallA(DMGetCoordinatesLocal(dm,coordVec,ierr))
      PetscCallA(DMGetDimension(dm,dim,ierr))
      PetscCallA(VecGetArrayF90(v,vArray,ierr))

      Do p = pStart,pEnd-1
         PetscCallA(PetscSectionGetDof(s,p,numDof,ierr))
         If (numDof > 0) Then
            !!! trick: the coordinate of a point is the average of the coordinates of the points in its closure
            PetscCallA(DMPlexVecGetClosure(dm,coordSection,coordVec,p,coordArray,ierr))
            Do i = 1,dim
                  xyz(i) = sum(coordArray(i:size(coordArray):dim)) * dim / size(coordArray)
            End Do
            PetscCallA(DMPlexVecRestoreClosure(dm,coordSection,coordVec,p,coordArray,ierr))

            PetscCallA(PetscSectionGetOffset(s,p,pOffset,ierr))
            Do i = 1,numDof
                  vArray(pOffset+i) = xyz(i)
            End Do
         End If
      End Do
      PetscCallA(VecRestoreArrayF90(v,vArray,ierr))
      !!! Of course, this does not use informations from the section, so it does over-write constrained values
   End Subroutine project

   Subroutine Integrate3D_Scal(MEF90Ctx,v,i,j,k,QuadratureOrder,i1,i2,ierr)
      Type(MEF90Ctx_Type),Intent(IN)                     :: MEF90Ctx
      Type(tVec),Intent(IN)                              :: v
      PetscInt,Intent(IN)                                :: QuadratureOrder,i,j,k
      PetscReal,Intent(OUT)                              :: i1
      Type(Vect3D),Intent(OUT)                           :: i2
      PetscErrorCode,Intent(OUT)                         :: ierr      
                     
      Type(tDM)                                          :: dm
      Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90CtxGlobalOptions
      Type(tIS)                                          :: setIS,setPointIS
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elem
      Type(MEF90ElementType)                             :: elementType
      DMPolytopeType                                     :: cellType
      PetscInt                                           :: iDof,iGauss,cell,set,dim
      PetscInt,Dimension(:),Pointer                      :: setID,setPointID
      PetscReal                                          :: X,Y,Z
      PetscReal,Dimension(:),Pointer                     :: coordDof

      i1 = 0.0_Kr
      i2 = 0.0_Kr
      
      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr))
      PetscCall(VecGetDM(v,dm,ierr))
      PetscCall(DMGetDimension(dm,dim,ierr))
      PetscCall(DMGetLabelIdIS(dm,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90Ctx%comm,setIS,ierr))
      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         Do set = 1,size(setID)
            PetscCall(DMGetStratumIS(dm,MEF90CellSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then

               PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(DMPlexGetCellType(dm,setPointID(1),cellType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellType,elementType,ierr))
               PetscCall(MEF90ElementCreate(dm,setPointIS,elem,QuadratureOrder,elementType,ierr))

               Do cell = 1, size(setPointID)
                  PetscCall(DMPlexVecGetClosure(dm,PETSC_NULL_SECTION,v,setPointID(cell),coordDof,ierr))
                  Do iGauss = 1,size(elem(cell)%Gauss_C)
                     X = 0.0_Kr
                     Y = 0.0_Kr
                     Z = 0.0_Kr
                     Do iDof = 1, size(elem(cell)%BF(:,1))
                        X = X + Elem(cell)%BF(iDoF,iGauss) * coordDof(3*(iDof-1)+1)
                        Y = Y + Elem(cell)%BF(iDoF,iGauss) * coordDof(3*(iDof-1)+2)
                        Z = Z + Elem(cell)%BF(iDoF,iGauss) * coordDof(3*(iDof-1)+3)
                     End Do
                     i1 = i1 + X**i * Y**j * Z**k * elem(cell)%Gauss_C(iGauss)
                  End Do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dm,PETSC_NULL_SECTION,v,setPointID(cell),coordDof,ierr))
               End Do ! cell

               PetscCall(MEF90ElementDestroy(elem,ierr))
               PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
            End If ! pointIS
            PetscCall(ISDestroy(setPointIS,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))
   End Subroutine Integrate3D_Scal

   Subroutine Integrate2D_Scal(MEF90Ctx,v,i,j,QuadratureOrder,i1,i2,ierr)
      Type(MEF90Ctx_Type),Intent(IN)                     :: MEF90Ctx
      Type(tVec),Intent(IN)                              :: v
      PetscInt,Intent(IN)                                :: QuadratureOrder,i,j
      PetscReal,Intent(OUT)                              :: i1
      Type(Vect2D),Intent(OUT)                           :: i2
      PetscErrorCode,Intent(OUT)                         :: ierr      
                     
      Type(tDM)                                          :: dm
      Type(MEF90CtxGlobalOptions_Type),pointer           :: MEF90CtxGlobalOptions
      Type(tIS)                                          :: setIS,setPointIS
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer      :: elem
      Type(MEF90ElementType)                             :: elementType
      DMPolytopeType                                     :: cellType
      PetscInt                                           :: iDof,iGauss,cell,set,dim
      PetscInt,Dimension(:),Pointer                      :: setID,setPointID
      PetscReal                                          :: X,Y
      PetscReal,Dimension(:),Pointer                     :: coordDof

      i1 = 0.0_Kr
      i2 = 0.0_Kr
      
      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(MEF90Ctx%GlobalOptionsBag,MEF90CtxGlobalOptions,ierr))
      PetscCall(VecGetDM(v,dm,ierr))
      PetscCall(DMGetDimension(dm,dim,ierr))
      PetscCall(DMGetLabelIdIS(dm,MEF90CellSetLabelName,setIS,ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90Ctx%comm,setIS,ierr))
      If (setIS /= PETSC_NULL_IS) Then
         PetscCall(ISGetIndicesF90(setIS,setID,ierr))
         Do set = 1,size(setID)
            PetscCall(DMGetStratumIS(dm,MEF90CellSetLabelName,setID(set),setPointIS,ierr))
            If (setPointIS /= PETSC_NULL_IS) Then

               PetscCall(ISGetIndicesF90(setPointIS,setPointID,ierr))
               PetscCall(DMPlexGetCellType(dm,setPointID(1),cellType,ierr))
               PetscCall(MEF90ElementGetType(MEF90CtxGlobalOptions%elementFamily,MEF90CtxGlobalOptions%elementOrder,cellType,elementType,ierr))
               PetscCall(MEF90ElementCreate(dm,setPointIS,elem,QuadratureOrder,elementType,ierr))

               Do cell = 1, size(setPointID)
                  PetscCall(DMPlexVecGetClosure(dm,PETSC_NULL_SECTION,v,setPointID(cell),coordDof,ierr))
                  Do iGauss = 1,size(elem(cell)%Gauss_C)
                     X = 0.0_Kr
                     Y = 0.0_Kr
                     Do iDof = 1, size(elem(cell)%BF(:,1))
                        X = X + Elem(cell)%BF(iDoF,iGauss) * coordDof(2*(iDof-1)+1)
                        Y = Y + Elem(cell)%BF(iDoF,iGauss) * coordDof(2*(iDof-1)+2)
                     End Do
                     i1 = i1 + X**i * Y**j * elem(cell)%Gauss_C(iGauss)
                  End Do ! iGauss
                  PetscCall(DMPlexVecRestoreClosure(dm,PETSC_NULL_SECTION,v,setPointID(cell),coordDof,ierr))
               End Do ! cell

               PetscCall(MEF90ElementDestroy(elem,ierr))
               PetscCall(ISRestoreIndicesF90(setPointIS,setPointID,ierr))
            End If ! pointIS
            PetscCall(ISDestroy(setPointIS,ierr))
         End Do ! set
         PetscCall(ISRestoreIndicesF90(setIS,setID,ierr))
      End If ! setIS
      PetscCall(ISDestroy(setIS,ierr))
   End Subroutine Integrate2D_Scal
End Program TestQuadrature
