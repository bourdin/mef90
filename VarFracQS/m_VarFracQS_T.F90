#if defined PB_2D
Module m_VarFracQS_T2D
#elif defined PB_3D
Module m_VarFracQS_T3D
#endif
#include "finclude/petscdef.h"

#if defined PB_2D
   Use m_TransientHeat2D
   Use m_VarFracQS_Types2D
#elif defined PB_3D
   Use m_TransientHeat3D
   Use m_VarFracQS_Types3D 
#endif   

   Use m_MEF90
   Use m_VarFrac_Struct

   Implicit NONE

Contains

#undef __FUNCT__
#define __FUNCT__ "VarFracHeat_Init"
   Subroutine VarFracHeat_Init(AppCtx, HeatAppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      Type(Heat_AppCtx_Type)                       :: HeatAppCtx
      
      Character(len=MEF90_MXSTRLEN)                :: TCST_FileName
      PetscInt, Dimension(:), Pointer              :: SizeScal
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
      PetscInt                                     :: iErr

      !Read Heat material parameters from *.TSCT file 
      TCST_FileName = trim(AppCtx%AppParam%prefix)//'.TCST'
      Call MatHeat_Read(AppCtx%MeshTopology, HeatAppCtx%MatProp, TCST_FileName)
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) "Done with MatHeat_Read\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If

      Call HeatSchemeParam_GetFromOptions(HeatAppCtx%HeatSchemeParam)
      If (AppCtx%AppParam%verbose > 0) Then
         Call HeatSchemeParam_View(HeatAppCtx%HeatSchemeParam, PetscViewer(PETSC_VIEWER_STDOUT_WORLD))
      End If
   
      HeatAppCtx%VertVar_Temperature = AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Temperature)%Offset
      Call ElementInit(AppCtx%MeshTopology, HeatAppCtx%Elem, 2)
  
      Call  HeatFieldCreate(HeatAppCtx, AppCtx%MeshTopology)

   !Read Initial Temerature Field 
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology,  AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Temperature)%Offset, 1, HeatAppCtx%T)
      Call SectionRealToVec(HeatAppCtx%T%Sec, HeatAppCtx%T%Scatter,  SCATTER_REVERSE, HeatAppCtx%T%Vec, ierr); CHKERRQ(ierr)
   
      !Set BC in Temperature
      Call SectionIntZero(HeatAppCtx%BCFlag%Sec, iErr); CHKERRQ(iErr)
      Call SectionIntAddNSProperty(HeatAppCtx%BCFlag%Sec,   AppCtx%MyEXO%NSProperty(VarFrac_NSProp_HasPForce), AppCtx%MeshTopology)

      Call Poisson_TSSetUp(HeatAppCtx, AppCtx%MeshTopology)
      Call RHSAssembly(HeatAppCtx, AppCtx%MeshTopology, AppCtx%MyExo)
      Call MatMassAssembly(HeatAppCtx, AppCtx%MeshTopology)

   End Subroutine VarFracHeat_Init

#undef __FUNCT__
#define __FUNCT__ "VarFracHeat_Step_Compute"
   Subroutine VarFracHeat_Step_Compute(AppCtx, HeatAppCtx)
      Type(AppCtx_Type)                            :: AppCtx
      Type(Heat_AppCtx_Type)                       :: HeatAppCtx
      
      PetscInt                                     :: iErr
      
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_ForceTemp)%Offset, AppCtx%TimeStep, HeatAppCtx%F)
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Temperature)%Offset,  AppCtx%TimeStep, HeatAppCtx%TBC) 

! Compute the crack opening here
      !!Compute Heat Field 
      Call Read_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_Load)%Offset, AppCtx%TimeStep-1, AppCtx%Load(AppCtx%TimeStep-1))
      Call Read_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_Load)%Offset, AppCtx%TimeStep, AppCtx%Load(AppCtx%TimeStep))
      Call MatZeroEntries(HeatAppCtx%K, iErr); CHKERRQ(iErr)
      Call HeatMatAssembly(HeatAppCtx, AppCtx%MeshTopology, AppCtx%V)
      Call SolveTransientStep(HeatAppCtx, AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%Load(AppCtx%TimeStep-1), AppCtx%Load(AppCtx%TimeStep), AppCtx%TimeStep-1)
      Call Write_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_Load)%Offset, AppCtx%TimeStep, AppCtx%Load(AppCtx%TimeStep))
      Call EXPTIM(AppCtx%MyEXO%exoid, AppCtx%TimeStep, AppCtx%Load(AppCtx%TimeStep), iErr)
!AppCtx%Load is the list of time steps. Quasi-Static ..... 
   End Subroutine VarFracHeat_Step_Compute

#if defined PB_2D
End Module m_VarFracQS_T2D
#elif defined PB_3D
End Module m_VarFracQS_T3D
#endif
