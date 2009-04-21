Program  Elast

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscmeshdef.h"

   Use m_MEF90
#if defined PB_2D
   Use m_Elast2D
#elif defined PB_3D 
   Use m_Elast3D
#endif

   Use petsc
   Use petscvec
   Use petscmat
   Use petscksp
   Use petscmesh

   Implicit NONE   


   Type(AppCtx_Type)                            :: AppCtx
   PetscInt                                     :: i, iErr
   Character(len=MEF90_MXSTRLEN)                :: IOBuffer

   Call ElastInit(AppCtx)
   
   If (AppCtx%AppParam%verbose) Then
      Call EXOView(AppCtx%EXO, AppCtx%AppParam%LogViewer)
      Call EXOView(AppCtx%MyEXO, AppCtx%AppParam%MyLogViewer)
      Call MeshTopologyView(AppCtx%MeshTopology, AppCtx%AppParam%MyLogViewer)
   End If   


   If (AppCtx%AppParam%verbose) Then
      Write(IOBuffer, *) 'Assembling the matrix\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If   
   Call MatAssembly(AppCtx)
   
   Do i = 1, AppCtx%NumTimeSteps
      AppCtx%TimeStep = i
      !!! Read U, F, and Temperature
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(Rupt_VertVar_DisplacementX)%Offset, AppCtx%TimeStep, AppCtx%U) 
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(Rupt_VertVar_ForceX)%Offset, AppCtx%TimeStep, AppCtx%F) 
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(Rupt_VertVar_Temperature)%Offset, AppCtx%TimeStep, AppCtx%Theta) 
   
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) 'Assembling the RHS\n'c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      Call RHSAssembly(AppCtx)
   
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) 'Calling KSPSolve\n'c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      Call Solve(AppCtx)
   
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) 'Computing Elastic energy, strains and stresses\n'c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      Call ComputeEnergy(AppCtx)
      Call ComputeStrainStress(AppCtx)
      
      Write(IOBuffer, 108) AppCtx%TimeStep, AppCtx%ElasticEnergy, AppCtx%ExtForcesWork, AppCtx%ElasticEnergy - AppCtx%ExtForcesWork
108 Format('TS ',G, ' Elast ', G, ' Work ', G, ' Total ', G, '\n'c)
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

      Write(IOBuffer, 110) AppCtx%TimeStep, AppCtx%ElasticEnergy, AppCtx%ExtForcesWork, AppCtx%ElasticEnergy - AppCtx%ExtForcesWork
110 Format(4(G), '\n'c)
      Call PetscViewerASCIIPrintf(AppCtx%AppParam%EnergyViewer, IOBuffer, iErr); CHKERRQ(iErr)

      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) 'Saving results\n'c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If      
      Call Save(AppCtx)
   End Do

   Call ElastFinalize(AppCtx)
End Program  Elast
