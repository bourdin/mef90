#if defined PB_2D
Module m_VarFracQS_T2D
#elif defined PB_3D
Module m_VarFracQS_T3D
#endif
#include "finclude/petscdef.h"

#if defined PB_2D
   Use m_Poisson2D
   Use m_TransientHeat2D
   Use m_VarFracQS_Types2D
#elif defined PB_3D
   Use m_Poisson3D
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
      !!! Read Mat Properties from the CST file
      TCST_FileName = trim(AppCtx%AppParam%prefix)//'.TCST'
      Call MatHeat_Read(AppCtx%MeshTopology, HeatAppCtx%MatProp, TCST_FileName)
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) "Done with MatHeat_Read\n"
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If


   !Init For HEAT
      Call HeatSchemeParam_GetFromOptions(HeatAppCtx%HeatSchemeParam)
      If (AppCtx%AppParam%verbose > 0) Then
         Call HeatSchemeParam_View(HeatAppCtx%HeatSchemeParam, PetscViewer(PETSC_VIEWER_STDOUT_WORLD))
      End If
   
      HeatAppCtx%VertVar_Temperature = AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Temperature)%Offset
      Call ElementInit(AppCtx%MeshTopology, HeatAppCtx%Elem, 2)
  
      Allocate(SizeScal(1)) 
      SizeScal=1
      Call FieldCreateVertex(HeatAppCtx%U,     'T',   AppCtx%MeshTopology, SizeScal)
      Call FieldCreateVertex(HeatAppCtx%F,     'F',  AppCtx%MeshTopology,  SizeScal)
      Call FieldCreateVertex(HeatAppCtx%RHS,   'RHS', AppCtx%MeshTopology, SizeScal)
      Call FlagCreateVertex(HeatAppCtx%BCFlag, 'BC',   AppCtx%MeshTopology, SizeScal)
      Call FieldCreateVertex(HeatAppCtx%UBC,    'TBC',      AppCtx%MeshTopology, SizeScal)
      DeAllocate(SizeScal)

   !Read Initial Temerature Field 
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology,  AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Temperature)%Offset, 1, HeatAppCtx%U)
      Call SectionRealToVec(HeatAppCtx%U%Sec, HeatAppCtx%U%Scatter,  SCATTER_REVERSE, HeatAppCtx%U%Vec, ierr); CHKERRQ(ierr)
   
   !      Call HeatSetInitial(HeatAppCtx, MeshTopology, ValT_Init,ValT_F)
      !Set BC in Temperature
!      Call HeatSetBC(HeatAppCtx, T_BC, MyExo, MeshTopology,VarFrac_NSProp_HasPForce)
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
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Temperature)%Offset,  AppCtx%TimeStep, HeatAppCtx%UBC) 

! Compute the crack opening here

      !!Compute Heat Field 
      If (AppCtx%TimeStep < 2) Then
!         Call SolveTransientStep(HeatAppCtx, AppCtx%MyEXO, AppCtx%MeshTopology, 0.0_Kr , AppCtx%Load(AppCtx%TimeStep), AppCtx%TimeStep)
      else 
         Call MatZeroEntries(HeatAppCtx%K, iErr); CHKERRQ(iErr)
         Call HeatMatAssembly(HeatAppCtx, AppCtx%MeshTopology, AppCtx%V)
         Call SolveTransientStep(HeatAppCtx, AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%Load(AppCtx%TimeStep-1), AppCtx%Load(AppCtx%TimeStep), AppCtx%TimeStep)
      End If 
!AppCtx%Load is the list of time steps. Quasi-Static ..... 
   End Subroutine VarFracHeat_Step_Compute

   
#if defined PB_2D
End Module m_VarFracQS_T2D
#elif defined PB_3D
End Module m_VarFracQS_T3D
#endif
