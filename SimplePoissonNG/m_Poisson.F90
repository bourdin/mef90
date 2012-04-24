#if MEF90_DIM == 3
#define M_POISSON m_Poisson3D
#else
#define M_POISSON m_Poisson2D
#endif

Module M_POISSON
#include "../MEF90/mef90.inc"
#include "finclude/petscdef.h"
   Use m_MEF90
   Implicit none
   Private  
   Public :: PoissonCtx_Type
   
   !Public :: PoissonTS_IFunction
   !Public :: PoissonTS_IJacobian
   !Public :: PoissonTS_RHS
   
   !Public :: PoissonSNES_Function
   !Public :: PoissonSNES_Jacobian
   
   !Public :: PoissonTAO_Function
   !Public :: PoissonTAO_Gradient
   !Public :: PoissonTAO_Hessian
   
   !Public :: PoissonKSP_Matrix
   !Public :: PoissonKSP_RHS
   
   Type PoissonCtx_Type
      Type(MEF90_ELEMENT_SCAL),Dimension(:),Pointer   :: Elem
      Type(Element_Type),Dimension(:),Pointer         :: ElementType
      Type(Field)                                     :: BCU
      Type(Field)                                     :: F
      Type(Flag)                                      :: BCFlag
   End Type PoissonCtx_Type
   
   
!!!
!!!   
   
Contains
#undef __FUNCT__
#define __FUNCT__ "PoissonMatAssembly"
!!! Make interface compatible with SNES
   Subroutine PoissonMatAssembly(PoissonSNES,X,K,Kpre,mStruct,AppCtx,ierr)   
      Type(SNES),intent(IN)                        :: PoissonSNES
      Type(Vec),intent(IN)                         :: X
      Type(Mat),intent(IN)                         :: K,Kpre
      MatStructure,intent(IN)                      :: mStruct
      Type(PoissonCtx_Type),Intent(IN)             :: AppCtx
      PetscInt,Intent(OUT)                         :: ierr
      
      PetscInt                                     :: iBlk
      Type(IS)                                     :: setIS
      PetscInt,Dimension(:),Pointer                :: setID
      PetscInt                                     :: set
      
      Call MatInsertVertexBoundaryValues(K,AppCtx%U,AppCtx%BCFlag,AppCtx%mesh)
      Call MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      Call MatAssemblyEnd  (K,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(AppCtx%mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         !!! Test if block corresponds to a volume or surface element
         If (AppCtx%ElementType(set)%codim == 0) Then
            Call PoissonMatAssemblyBlock(K,setID(set),AppCtx)
         End If
      End Do
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      Call MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      Call MatAssemblyEnd  (K,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
   End Subroutine PoissonMatAssembly
      
#undef __FUNCT__
#define __FUNCT__ "PoissonMatAssemblyBlock"
   Subroutine PoissonMatAssemblyBlock(K,setID,AppCtx)
      Type(Mat),Intent(IN)                         :: K
      PetscInt,Intent(IN)                          :: setID
      Type(PoissonCtx_Type),Intent(IN)             :: AppCtx
      
      Type(IS)                                     :: cellIS
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscInt                                     :: numDof
      PetscInt                                     :: ierr,i
      PetscReal,Dimension(:,:),Pointer             :: MatElem
      PetscInt,Dimension(:),Pointer                :: BCFlag
      PetscInt                                     :: iDoF1,iDoF2,iGauss
      PetscLogDouble                               :: flops = 0
      PetscReal                                    :: lDiff
      PetscReal,Dimension(:),Pointer               :: T_Loc
      PetscReal                                    :: T_Elem
     
      Call PetscOptions
      numDof = AppCtx%ElementType(iBlk)%numDof
      Allocate(MatElem(numDof,numDof))
      Allocate(BCFlag(numDof))
      Allocate(T_Loc(numDof))
      Call DMmeshGetStratumIS(AppCtx%mesh,'Cell Sets',setID,cellIS,ierr); CHKERRQ(ierr)
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Do cell = 1,size(cellID)      
         MatElem = 0.0_Kr
         BCFlag = 0
         Call SectionIntRestrictClosure(AppCtx%BCFlag%Sec,AppCtx%mesh,cellID(cell),numDof,BCFlag,ierr);CHKERRQ(ierr)
         Do iGauss = 1,size(AppCtx%Elem(cellID(cell)+1)%Gauss_C)
            Select Case(AppCtx%AppParam%TestCase)
            Case(2)
               lDiff = AppCtx%Diff(iBlk) 
               !!! THIS ASSUMES THAT BLOCKS ARE NUMBERED SEQUENTIALLY
               !!! Use PetscOptions to def lDiff, instead
            Case(3)
               ! Mensi Law D(C) = A exp (B*C) 
               ! 1988 Mensi-Acker-Attolou Mater Struct
               ! C : water content, A and B material parameters
               Call SectionRealRestrictClosure(AppCtx%U%Sec,AppCtx%mesh,cellID(cell),numDoF,T_Loc,ierr);CHKERRQ(ierr)
               T_Elem = 0.0_Kr
               Do iDoF1 = 1,numDof
                  T_Elem = T_Elem + AppCtx%Elem(cellID(cell)+1)%BF(iDoF1,iGauss) * T_Loc(iDoF1)
               End DO
               lDiff = AppCtx%Diff(iBlk)*exp(AppCtx%B_Mensi(iBlk)*T_Elem) 
               !!! SAME AS ABOVE
            Case Default
               ldiff =1 
            End Select
            Do iDoF1 = 1,numDof
               If (BCFlag(iDoF1) == 0) Then
                  Do iDoF2 = 1,numDof
                    ! MatElem(iDoF1,iDoF1) = 1./2.
                     MatElem(iDoF2,iDoF1) = MatElem(iDoF2,iDoF1) + lDiff * AppCtx%Elem(cellID(cell)+1)%Gauss_C(iGauss) * &
                                           (AppCtx%Elem(cellID(cell)+1)%Grad_BF(iDoF1,iGauss) .DotP. AppCtx%Elem(cellID(cell)+1)%Grad_BF(iDoF2,iGauss) )
                  End Do
                  flops = flops + 3 * numDof
               End If
            End Do
         End Do
         Call DMmeshAssembleMatrix(K,AppCtx%mesh,AppCtx%U%Sec,cellID(cell),MatElem,ADD_VALUES,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Call ISDestroy(cellIS,ierr);CHKERRQ(ierr)
   
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      DeAllocate(BCFlag)
      DeAllocate(MatElem)
   End Subroutine PoissonMatAssemblyBlock


#undef __FUNCT__
#define __FUNCT__ "PoissonRHSAssembly"
!!! Change interface.
   Subroutine PoissonRHSAssembly(SNESU,RHSVec,UVec,AppCtx)
      Type(SNES),intent(IN)                        :: SNESU
      Type(Vec),intent(IN)                         :: RHSVec
      !!! Technically, RHS is an INOUT arg, but the value of RHSVec
      !!! should not be modified by this function
      Type(Vec),intent(IN)                         :: UVec
      Type(PoissonCtx_Type)                    :: AppCtx

      PetscInt                                     :: ierr
      PetscInt                                     :: set
      PetscInt,Dimension(:),Pointer                :: setID
      Type(IS)                                     :: setIS

      Call SectionRealZero(AppCtx%RHS%Sec,ierr);CHKERRQ(ierr)
      
      Call DMmeshGetLabelIdIS(AppCtx%mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         If (AppCtx%ElementType(set)%codim > 0) Then
            !!! Only surface forces
            Call RHSAssemblyBlock(AppCtx%RHS,setID(set),AppCtx)
         End If
      End Do ! set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)

      Call SectionRealComplete(AppCtx%RHS%Sec,ierr);CHKERRQ(ierr)

      Call FieldInsertVertexBoundaryValues(AppCtx%RHS,AppCtx%U,AppCtx%BCFlag,AppCtx%mesh)
      Call SectionRealToVec(AppCtx%RHS%Sec,AppCtx%RHS%Scatter,SCATTER_FORWARD,RHSVec,ierr);CHKERRQ(ierr)
      !Call VecCopy(AppCtx%RHS%Vec,RHSVec,ierr);CHKERRQ(ierr)
   End Subroutine PoissonRHSAssembly

#undef __FUNCT__
#define __FUNCT__ "PoissonRHSAssemblyBlock"
   Subroutine PoissonRHSAssemblyBlock(RHS,setID,AppCtx)
      Type(Field),Intent(IN)                       :: RHS
      PetscInt,Intent(IN)                          :: setID
      Type(PoissonCtx_Type),Intent(IN)             :: AppCtx

      PetscInt                                     :: cell
      PetscInt,Dimension(:),Pointer                :: cellID
      Type(IS)                                     :: cellIS
      PetscInt                                     :: ierr
      PetscInt                                     :: numDof,iDoF
      PetscInt                                     :: iGauss
      PetscInt,Dimension(:),Pointer                :: BCFlag_Loc
      PetscReal,Dimension(:),Pointer               :: RHS_Loc,F_Loc
      PetscReal                                    :: F_Elem
      PetscLogDouble                               :: flops 
      
      Call PetscLogEventBegin(AppCtx%LogInfo%RHSAssemblyBlock_Event,ierr);CHKERRQ(ierr)
      flops = 0.0

      numDof = AppCtx%ElementType(setID)%numDof
      Allocate(F_Loc(numDof))
      Allocate(RHS_Loc(numDof))
      Allocate(BCFlag_Loc(numDof))

      Call DMmeshGetStratumIS(AppCtx%mesh,'Cell Sets',setID,cellIS,ierr); CHKERRQ(ierr)
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Do cell = 1,size(cellID)      
         RHS_Loc = 0.0_Kr
         Call SectionRealRestrictClosure(AppCtx%F%Sec,AppCtx%mesh,cellID(cell),numDof,F_Loc,ierr);CHKERRQ(ierr)
         Call SectionIntRestrictClosure(AppCtx%BCFlag%Sec,AppCtx%mesh,cellID(cell),numDof,BCFlag_Loc,ierr);CHKERRQ(ierr)
         Do iGauss = 1,Size(AppCtx%Elem(cellID(cell)+1)%Gauss_C)
            F_Elem = 0.0_Kr
            Do iDoF = 1,numDof
               F_Elem = F_Elem + AppCtx%Elem(cellID(cell)+1)%BF(iDoF,iGauss) * F_Loc(iDoF)
               flops = flops + 2.0
            End Do
            Do iDoF = 1,numDof
               If (BCFlag_Loc(iDoF) == 0) Then
                  RHS_Loc(iDoF) = RHS_Loc(iDoF) + AppCtx%Elem(cellID(cell)+1)%Gauss_C(iGauss) * ( F_Elem * AppCtx%Elem(cellID(cell)+1)%BF(iDoF,iGauss) )
                  flops = flops + 3.0
               End If
            End Do
         End Do
         Call SectionRealUpdateClosure(AppCtx%RHS%Sec,AppCtx%mesh,cellID(cell),RHS_Loc,ADD_VALUES,ierr);CHKERRQ(ierr)
      End Do ! cell
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Call ISDestroy(cellIS,ierr);CHKERRQ(ierr)

      DeAllocate(BCFlag_Loc)
      DeAllocate(RHS_Loc)
      DeAllocate(F_Loc)
      Call PetscLogFlops(flops ,ierr);CHKERRQ(ierr)
      Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyBlock_Event,ierr);CHKERRQ(ierr)
   End Subroutine PoissonRHSAssemblyBlock

#undef __FUNCT__
#define __FUNCT__ "PoissonEnergy"
!!! Change interface to make compatible with TAO
   Subroutine PoissonEnergy(AppCtx)
      Type(PoissonCtx_Type)                    :: AppCtx
      
      PetscInt                                     :: ierr
      PetscInt                                     :: NumDoF
      PetscReal,Dimension(:),Pointer               :: F,U
      PetscInt                                     :: set,cell
      PetscInt,Dimension(:),Pointer                :: setID,cellID
      Type(IS)                                     :: setIS,cellIS
      PetscInt                                     :: iDoF,iGauss
      Type(MEF90_VECT)                             :: Strain_Elem,Stress_Elem      
      PetscReal                                    :: F_Elem,U_Elem
      PetscReal                                    :: MyElasticEnergy,MyExtForcesWork
      PetscLogDouble                               :: flops

      MyElasticEnergy = 0.0_Kr
      MyExtForcesWork = 0.0_Kr
      flops = 0.0_Kr
      Call DMmeshGetLabelIdIS(AppCtx%mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         numDof = AppCtx%ElementType(setID(set))%numDof
         Call DMmeshGetStratumIS(AppCtx%mesh,'Cell Sets',setID(set),cellIS,ierr); CHKERRQ(ierr)
         Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         Do cell = 1,size(cellID)      
            Allocate(F(NumDoF))
            Call SectionRealRestrictClosure(AppCtx%F%Sec,AppCtx%mesh,cellID(cell),NumDoF,F,ierr);CHKERRQ(ierr)
            Allocate(U(NumDoF))
            Call SectionRealRestrictClosure(AppCtx%U%Sec,AppCtx%mesh,cellID(cell),NumDoF,U,ierr);CHKERRQ(ierr)
            Do iGauss = 1,Size(AppCtx%Elem(cellID(cell)+1)%BF,2)
               Strain_Elem = 0.0_Kr
               Stress_Elem = 0.0_Kr
               F_Elem      = 0.0_Kr
               U_Elem      = 0.0_Kr
               Do iDoF = 1,NumDoF
                  If (AppCtx%ElementType(set)%codim == 0) Then
                     Stress_Elem = Stress_Elem + AppCtx%Elem(cellID(cell)+1)%Grad_BF(iDoF,iGauss) * U(iDoF)
                     Strain_Elem = Strain_Elem + AppCtx%Elem(cellID(cell)+1)%Grad_BF(iDoF,iGauss) * U(iDoF)
                  End If
                  F_Elem = F_Elem + AppCtx%Elem(cellID(cell)+1)%BF(iDoF,iGauss) * F(iDoF)
                  U_Elem = U_Elem + AppCtx%Elem(cellID(cell)+1)%BF(iDoF,iGauss) * U(iDoF)
                  flops = flops + 4
               End Do
               MyElasticEnergy = MyElasticEnergy + AppCtx%Elem(cellID(cell)+1)%Gauss_C(iGauss) * ( (Strain_Elem .DotP. Strain_Elem) * 0.5_Kr)
               MyExtForcesWork = MyExtForcesWork - AppCtx%Elem(cellID(cell)+1)%Gauss_C(iGauss) * F_Elem * U_Elem
               flops = flops + 5
            End Do
            DeAllocate(U)
            DeAllocate(F)
         End Do ! cell
         Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         Call ISDestroy(cellIS,ierr);CHKERRQ(ierr)
      End Do ! set
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      Call MPI_AllReduce(MyElasticEnergy,AppCtx%ElasticEnergy,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)
      Call MPI_AllReduce(MyExtForcesWork,AppCtx%ExtForcesWork,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)
      AppCtx%TotalEnergy = AppCtx%ElasticEnergy + AppCtx%ExtForcesWork
      Call PetscLogFlops(flops,ierr)
   End Subroutine PoissonEnergy
 
End Module M_POISSON
