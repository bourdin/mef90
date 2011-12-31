#if defined PB_2D
Module m_Poisson2D
#elif defined PB_3D
Module m_Poisson3D
#endif

#include "finclude/petscdef.h"

   Use m_MEF90
   Use m_Heat_Struct 
   
   Implicit NONE   

   Type AppParam_Type
      PetscBool                                    :: Restart
      PetscInt                                     :: Verbose
      PetscInt                                     :: TestCase
      Character(len=MEF90_MXSTRLEN)                :: prefix
      Type(PetscViewer)                            :: LogViewer, MyLogViewer
   End Type AppParam_Type


   Type Heat_AppCtx_Type
      Type (MeshTopology_Type)                     :: MeshTopology
      Type (EXO_Type)                              :: EXO, MyEXO
#if defined PB_2D
      Type(Element2D_Scal), Dimension(:), Pointer  :: Elem
#elif defined PB_3D
      Type(Element3D_Scal), Dimension(:), Pointer  :: Elem
#endif
      Type(SectionReal)                            :: GradU
      PetscReal                                    :: ElasticEnergy
      Type(Field)                                  :: U
      Type(Field)                                  :: UBC
      Type(Field)                                  :: F
      PetscReal                                    :: ExtForcesWork
      PetscReal                                    :: TotalEnergy
      Type(Flag)                                   :: BCFlag
      Type(Mat)                                    :: K
      Type(Mat)                                    :: M ! Mass matrix 
      Type(Mat)                                    :: Jac ! jacobian
      Type(Field)                                  :: RHS
      Type(KSP)                                    :: KSP
      Type(PC)                                     :: PC
      Type(AppParam_Type)                          :: AppParam
   !For TS
      Type(TS)                                     :: TS
      PetscInt                                     :: NumSteps
      PetscReal                                    :: maxtime
      PetscInt                                     :: VertVar_Temperature 
      Type(MatHeat_Type), Dimension(:), Pointer    :: MatProp
      Type(HeatSchemeParam_Type)                   :: HeatSchemeParam 
   End Type Heat_AppCtx_Type
   
   
Contains

   
#undef __FUNCT__
#define __FUNCT__ "ExoFormat_SimplePoisson"
   Subroutine EXOFormat_SimplePoisson(AppCtx)
      Type(Heat_AppCtx_Type)                            :: AppCtx
      PetscInt                                     :: iErr
   
      Call EXPVP (AppCtx%MyEXO%exoid, 'g', 3, iErr)
      Call EXPVAN(AppCtx%MyEXO%exoid, 'g', 3, (/'Elastic Energy ', 'Ext Forces work', 'Total Energy   '/), iErr)
      Call EXPVP (AppCtx%MyEXO%exoid, 'n', 2, iErr)
      Call EXPVAN(AppCtx%MyEXO%exoid, 'n', 2, (/'U', 'F'/), iErr)
#if defined PB_2D
      Call EXPVP (AppCtx%MyEXO%exoid, 'e', 2, iErr)
      Call EXPVAN(AppCtx%MyEXO%exoid, 'e', 2, (/'Grad U_X', 'Grad U_Y'/), iErr)
#elif defined PB_3D
      Call EXPVP (AppCtx%MyEXO%exoid, 'e', 3, iErr)
      Call EXPVAN(AppCtx%MyEXO%exoid, 'e', 3, (/'Grad U_X', 'Grad U_Y', 'Grad U_Z'/), iErr)
#endif
      Call EXPTIM(AppCtx%MyEXO%exoid, 1, 1.0_Kr, iErr)
   End Subroutine EXOFormat_SimplePoisson
   
 
#undef __FUNCT__
#define __FUNCT__ "HeatMatAssembly"
   Subroutine HeatMatAssembly(AppCtx, MeshTopology, ExtraField)   
      Type(Heat_AppCtx_Type)                             :: AppCtx
      Type (MeshTopology_Type)                           :: MeshTopology
      PetscInt                                           :: iBlk, iErr
      Type(Field)                                        :: ExtraField
      
      Call MatInsertVertexBoundaryValues(AppCtx%K, AppCtx%U, AppCtx%BCFlag, MeshTopology)
      Call MatAssemblyBegin(AppCtx%K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (AppCtx%K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)

      Do_iBlk: Do iBlk = 1, MeshTopology%Num_Elem_Blks
         Call HeatMatAssemblyBlock(iBlk, AppCtx, MeshTopology, ExtraField)
      End Do Do_iBlk
      Call MatAssemblyBegin(AppCtx%K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)
      Call MatAssemblyEnd  (AppCtx%K, MAT_FINAL_ASSEMBLY, iErr); CHKERRQ(iErr)

   End Subroutine HeatMatAssembly
      
   Subroutine HeatMatAssemblyBlock(iBlk, AppCtx, MeshTopology, ExtraField)
      Type(Heat_AppCtx_Type)                       :: AppCtx
      Type (MeshTopology_Type)                     :: MeshTopology
      PetscInt                                     :: iBlk
      Type(Field)                                  :: ExtraField
      
      PetscInt                                     :: iE, iELoc, iErr, i
      PetscReal, Dimension(:,:), Pointer           :: MatElem
      PetscInt, Dimension(:), Pointer              :: BCFlag
      PetscInt                                     :: iDoF1, iDoF2, iGauss
      PetscInt                                     :: NumDoFScal
      PetscLogDouble                               :: flops = 0
      PetscReal                                    :: lDiff
      PetscReal, Dimension(:), Pointer             :: T_Loc
      PetscReal                                    :: T_Elem
     
      NumDoFScal = MeshTopology%Elem_Blk(iBlk)%Num_DoF
      Allocate(MatElem(MeshTopology%Elem_Blk(iBlk)%Num_DoF, MeshTopology%Elem_Blk(iBlk)%Num_DoF))
      Allocate(BCFlag(MeshTopology%Elem_Blk(iBlk)%Num_DoF))
      Allocate(T_Loc(NumDoFScal))

      Do_iELoc: Do iELoc = 1, MeshTopology%Elem_Blk(iBlk)%Num_Elems
         i = MeshTopology%Elem_Blk(iBlk)%ID
         iE = MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         MatElem = 0.0_Kr
         BCFlag = 0
         Call SectionIntRestrictClosure(AppCtx%BCFlag%Sec, MeshTopology%mesh, iE-1, MeshTopology%Elem_Blk(iBlk)%Num_DoF, BCFlag, iErr); CHKERRQ(ierr)
         Do iGauss = 1, size(AppCtx%Elem(iE)%Gauss_C)
         Select Case(AppCtx%MatProp(i)%Type_Law)
            Case(Heat_Constant_ID)
            ! Constant Diffusion
               lDiff = AppCtx%MatProp(i)%Diffusivity(1) 
            Case(Heat_Increasing_ID )
            !Increasing diffusion 
               ! Mensi Law D(C) = A exp (B*C) 
               ! 1988 Mensi-Acker-Attolou Mater Struct
               ! C : water content,  A and B material parameters
               Call SectionRealRestrictClosure(AppCtx%U%Sec, MeshTopology%mesh,  iE-1, NumDoFScal, T_Loc, iErr); CHKERRQ(ierr)
               T_Elem = 0.0_Kr
               Do iDoF1 = 1, NumDoFScal
                  T_Elem = T_Elem + AppCtx%Elem(iE)%BF(iDoF1, iGauss) * T_Loc(iDoF1)
               End DO
               lDiff = AppCtx%MatProp(i)%Diffusivity(1)*exp(AppCtx%MatProp(i)%Diffusivity(2)*T_Elem) 
            Case(Heat_Decreasing_ID)
      !Diffusion is decreasing with the variable
               Call SectionRealRestrictClosure(AppCtx%U%Sec, MeshTopology%mesh,  iE-1, NumDoFScal, T_Loc, iErr); CHKERRQ(ierr)
               T_Elem = 0.0_Kr
               Do iDoF1 = 1, NumDoFScal
                  T_Elem = T_Elem + AppCtx%Elem(iE)%BF(iDoF1, iGauss) * T_Loc(iDoF1)
               End DO
               lDiff =  AppCtx%MatProp(i)%Diffusivity(1)*(AppCtx%MatProp(i)%Diffusivity(2)-T_Elem)**AppCtx%MatProp(i)%Diffusivity(3) 
            Case(Heat_Non_Monotonic_ID)
      !Diffusion is non monotonic with the variable
               Call SectionRealRestrictClosure(AppCtx%U%Sec, MeshTopology%mesh,  iE-1, NumDoFScal, T_Loc, iErr); CHKERRQ(ierr)
               T_Elem = 0.0_Kr
               Do iDoF1 = 1, NumDoFScal
                  T_Elem = T_Elem + AppCtx%Elem(iE)%BF(iDoF1, iGauss) * T_Loc(iDoF1)
               End DO
      !TODO implement 2007_CCR_baroghel-bouny_I and II         
               lDiff = AppCtx%MatProp(i)%Diffusivity(1)*(AppCtx%MatProp(i)%Diffusivity(2)-T_Elem)**AppCtx%MatProp(i)%Diffusivity(3)+AppCtx%MatProp(i)%Diffusivity(4)*exp(AppCtx%MatProp(i)%Diffusivity(5)*T_Elem) 
           Case(Heat_Damage_ID)
      ! Diffusion Depends on the damaging variable
      ! The problem is that we do not have access to that field here    !!!!!!
      ! - The function should take a section as argument. This section could then  be any field -- > Try this 
      ! - merge both AppCtx
      ! - copy to HeatAppCtx the damaging and displacement field 
               Call SectionRealRestrictClosure(ExtraField%Sec, MeshTopology%mesh,  iE-1, NumDoFScal, T_Loc, iErr); CHKERRQ(ierr)
               T_Elem = 0.001_Kr
               Do iDoF1 = 1, NumDoFScal
                  T_Elem = T_Elem + AppCtx%Elem(iE)%BF(iDoF1, iGauss) * T_Loc(iDoF1)
               End DO
               lDiff = min(AppCtx%MatProp(i)%Diffusivity(1)/(T_Elem+0.001), AppCtx%MatProp(i)%Diffusivity(2))
      ! TODO : Test that ExtraField is defined to evoid Segment Fault error
            Case(Heat_Crack_Open_ID)
      ! Diffusion depends on crack opening
            Case Default
               ldiff =1 
            End Select
            Do iDoF1 = 1, MeshTopology%Elem_Blk(iBlk)%Num_DoF
               If (BCFlag(iDoF1) == 0) Then
                  Do iDoF2 = 1, MeshTopology%Elem_Blk(iBlk)%Num_DoF
                    ! MatElem(iDoF1, iDoF1) = 1./2.
                     MatElem(iDoF2, iDoF1) = MatElem(iDoF2, iDoF1) + lDiff * AppCtx%Elem(iE)%Gauss_C(iGauss) * ( AppCtx%Elem(iE)%Grad_BF(iDoF1, iGauss) .DotP. AppCtx%Elem(iE)%Grad_BF(iDoF2, iGauss) )
                     flops = flops + 1
                  End Do
               End If
            End Do
         End Do
         Call DMMeshAssembleMatrix(AppCtx%K, MeshTopology%mesh, AppCtx%U%Sec, iE-1, MatElem, ADD_VALUES, iErr); CHKERRQ(iErr)
         !Call DMMeshAssembleMatrix(AppCtx%K, MeshTopology%mesh, AppCtx%U%Sec, iE-1, MatElem, INSERT_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iELoc
   
      Call PetscLogFlops(flops, iErr);CHKERRQ(iErr)
      DeAllocate(BCFlag)
      DeAllocate(MatElem)
   End Subroutine HeatMatAssemblyBlock


#undef __FUNCT__
#define __FUNCT__ "RHSAssembly"
   Subroutine RHSAssembly(AppCtx, MeshTopology, MyExo)
      Type(Heat_AppCtx_Type)                       :: AppCtx
      Type (MeshTopology_Type)                     :: MeshTopology
      Type (EXO_Type)                              :: MyEXO

      PetscInt                                     :: iErr
      PetscInt                                     :: iBlk

      Call SectionRealZero(AppCtx%RHS%Sec, iErr); CHKERRQ(iErr)
      
      Do_iBlk: Do iBlk = 1, MeshTopology%Num_Elem_Blks
         Call RHSAssemblyBlock(iBlk, AppCtx, MeshTopology)
      End Do Do_iBlk

      Call SectionRealComplete(AppCtx%RHS%Sec, iErr); CHKERRQ(iErr)

      !!! Set Dirichlet Boundary Values
!Suppose that loading is contant (replace second to last by AppCtx%Timestep otherwise)
      Call Read_EXO_Result_Vertex(MyEXO, MeshTopology, AppCtx%VertVar_Temperature, 1, AppCtx%UBC)
      Call FieldInsertVertexBoundaryValues(AppCtx%RHS, AppCtx%UBC, AppCtx%BCFlag, MeshTopology)

      Call SectionRealToVec(AppCtx%RHS%Sec, AppCtx%RHS%Scatter, SCATTER_FORWARD, AppCtx%RHS%Vec, iErr); CHKERRQ(iErr)
      Call FieldInsertVertexBoundaryValues(AppCtx%RHS, AppCtx%UBC, AppCtx%BCFlag, MeshTopology)
      !!! VERY important! This is the equivalent of a ghost update
   End Subroutine RHSAssembly

#undef __FUNCT__
#define __FUNCT__ "RHSAssemblyBlock"
   Subroutine RHSAssemblyBlock(iBlk, AppCtx, MeshTopology)
      PetscInt                                     :: iBlk
      Type(Heat_AppCtx_Type)                       :: AppCtx
      Type (MeshTopology_Type)                     :: MeshTopology

      PetscInt                                     :: iErr
      PetscInt                                     :: iBlkId
      PetscInt                                     :: Num_DoF, iDoF
      PetscInt                                     :: iEloc, iE, iGauss
      PetscInt, Dimension(:), Pointer              :: BCFlag_Loc
      PetscReal, Dimension(:), Pointer             :: RHS_Loc, F_Loc
      PetscReal                                    :: F_Elem
      PetscLogDouble                               :: flops 
      
      flops = 0.0

      Num_DoF = MeshTopology%Elem_Blk(iBlk)%Num_DoF
      Allocate(F_Loc(Num_DoF))
      Allocate(RHS_Loc(Num_DoF))
      Allocate(BCFlag_Loc(Num_DoF))

      iBlkID = MeshTopology%Elem_Blk(iBlk)%ID
      Do_iEloc: Do iELoc = 1, MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         RHS_Loc = 0.0_Kr
         Call SectionRealRestrictClosure(AppCtx%F%Sec, MeshTopology%mesh, iE-1, Num_DoF, F_Loc, iErr); CHKERRQ(ierr)
         Call SectionIntRestrictClosure(AppCtx%BCFlag%Sec, MeshTopology%mesh, iE-1, Num_DoF, BCFlag_Loc, iErr); CHKERRQ(ierr)
         Do iGauss = 1, Size(AppCtx%Elem(iE)%Gauss_C)
            F_Elem = 0.0_Kr
            Do iDoF = 1, Num_DoF
               F_Elem = F_Elem + AppCtx%Elem(iE)%BF(iDoF, iGauss) * F_Loc(iDoF)
               flops = flops + 2.0
            End Do
            Do iDoF = 1, Num_DoF
               If (BCFlag_Loc(iDoF) == 0) Then
                  RHS_Loc(iDoF) = RHS_Loc(iDoF) + AppCtx%Elem(iE)%Gauss_C(iGauss) * ( F_Elem * AppCtx%Elem(iE)%BF(iDoF, iGauss) )
                  flops = flops + 3.0
               End If
            End Do
         End Do
         Call SectionRealUpdateClosure(AppCtx%RHS%Sec, MeshTopology%Mesh, iE-1, RHS_Loc, ADD_VALUES, iErr); CHKERRQ(iErr)
      End Do Do_iEloc

      DeAllocate(BCFlag_Loc)
      DeAllocate(RHS_Loc)
      DeAllocate(F_Loc)
   End Subroutine RHSAssemblyBlock

#undef __FUNCT__
#define __FUNCT__ "ComputeEnergy"
   Subroutine ComputeEnergy(AppCtx)
      Type(Heat_AppCtx_Type)                            :: AppCtx
      
      PetscInt                                     :: iErr
      PetscInt                                     :: NumDoF, NumGauss
      PetscReal, Dimension(:), Pointer             :: F, U
      PetscInt                                     :: iBlk, iELoc, iE
      PetscInt                                     :: iDoF, iGauss
#if defined PB_2D
      Type(Vect2D)                                 :: Strain_Elem, Stress_Elem      
#elif defined PB_3D
      Type(Vect3D)                                 :: Strain_Elem, Stress_Elem      
#endif
      PetscReal                                    :: F_Elem, U_Elem
      PetscReal                                    :: MyElasticEnergy, MyExtForcesWork
      PetscLogDouble                               :: flops

      MyElasticEnergy = 0.0_Kr
      MyExtForcesWork = 0.0_Kr
      flops = 0.0_Kr
      Do_Elem_iBlk: Do iBlk = 1, AppCtx%MeshTopology%Num_Elem_Blks
         Do_Elem_iE: Do iELoc = 1, AppCtx%MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iE = AppCtx%MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
            NumDoF   = Size(AppCtx%Elem(iE)%BF,1)
            NumGauss = Size(AppCtx%Elem(iE)%BF,2)
            Allocate(F(NumDoF))
            Call SectionRealRestrictClosure(AppCtx%F%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoF, F, iErr); CHKERRQ(ierr)
            Allocate(U(NumDoF))
            Call SectionRealRestrictClosure(AppCtx%U%Sec, AppCtx%MeshTopology%mesh, iE-1, NumDoF, U, iErr); CHKERRQ(ierr)
            Do iGauss = 1, NumGauss
               Strain_Elem = 0.0_Kr
               Stress_Elem = 0.0_Kr
               F_Elem      = 0.0_Kr
               U_Elem      = 0.0_Kr
               Do iDoF = 1, NumDoF
                  Stress_Elem = Stress_Elem + AppCtx%Elem(iE)%Grad_BF(iDoF, iGauss) * U(iDoF)
                  Strain_Elem = Strain_Elem + AppCtx%Elem(iE)%Grad_BF(iDoF, iGauss) * U(iDoF)
                  F_Elem = F_Elem + AppCtx%Elem(iE)%BF(iDoF, iGauss) * F(iDoF)
                  U_Elem = U_Elem + AppCtx%Elem(iE)%BF(iDoF, iGauss) * U(iDoF)
                  flops = flops + 4
               End Do
               MyElasticEnergy = MyElasticEnergy + AppCtx%Elem(iE)%Gauss_C(iGauss) * ( (Strain_Elem .DotP. Strain_Elem) * 0.5_Kr)
               MyExtForcesWork = MyExtForcesWork - AppCtx%Elem(iE)%Gauss_C(iGauss) * F_Elem * U_Elem
               flops = flops + 5
            End Do
            DeAllocate(U)
            DeAllocate(F)
         End Do Do_Elem_iE
      End Do Do_Elem_iBlk
      Call MPI_AllReduce(MyElasticEnergy, AppCtx%ElasticEnergy, 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      Call MPI_AllReduce(MyExtForcesWork, AppCtx%ExtForcesWork, 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)
      AppCtx%TotalEnergy = AppCtx%ElasticEnergy + AppCtx%ExtForcesWork
   End Subroutine ComputeEnergy
 
#if defined PB_2D
End Module m_Poisson2D
#elif defined PB_3D
End Module m_Poisson3D
#endif
