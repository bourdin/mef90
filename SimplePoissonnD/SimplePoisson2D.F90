#if defined PB_2D
Program  SimplePoisson2D
#elif defined PB_3D
Program SimplePoisson3D
#endif

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscmeshdef.h"

   Use m_MEF90
#if defined PB_2D
   Use m_SimplePoisson2D
#elif defined PB_3D
   Use m_SimplePoisson3D
#endif

   Use petsc
   Use petscvec
   Use petscmat
   Use petscksp
   Use petscmesh

   Implicit NONE   


   Type(AppCtx_Type)                            :: AppCtx
   PetscInt                                     :: iErr
   Type(Vec)                                    :: F
   Character(len=MEF90_MXSTRLEN)                :: IOBuffer
   Type(Vec)                                    :: V

   Call SimplePoissonInit(AppCtx)
   
   If (AppCtx%AppParam%verbose) Then
      Call EXOView(AppCtx%EXO, AppCtx%AppParam%LogViewer)
      Call EXOView(AppCtx%MyEXO, AppCtx%AppParam%MyLogViewer)
!      Call MeshTopologyView(AppCtx%MeshTopology, AppCtx%AppParam%MyLogViewer)
   End If   

   If (AppCtx%AppParam%verbose) Then
      Write(IOBuffer, *) 'Assembling the matrix\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   Call MatAssembly(AppCtx)
   
   If (AppCtx%AppParam%verbose) Then
      Write(IOBuffer, *) 'Assembling the RHS\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   Call SectionRealCreateLocalVector(AppCtx%F, F, iErr); CHKERRQ(iErr)
   Call VecSet(F, 1.0_Kr, iErr); CHKERRQ(iErr);
   Call VecDestroy(F, iErr); CHKERRQ(iErr)
   
   Call RHSAssembly(AppCtx)
   
   If (AppCtx%AppParam%verbose) Then
      Write(IOBuffer, *) 'Calling KSPSolve\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   Call Solve(AppCtx)
   Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
   Call Write_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, 1, 1, AppCtx%U) 
   Call PetscLogStagePop (AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)

   If (AppCtx%AppParam%verbose) Then
      Write(IOBuffer, *) 'Computing energy\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   Call ComputeEnergy(AppCtx)
   Write(IOBuffer, 100) AppCtx%Energy
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

100 Format('Total energy: ', ES12.5, '\n'c)    
   If (AppCtx%AppParam%verbose) Then
      Write(IOBuffer, *) 'Saving U\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   Call PetscLogStagePush(AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)
   Call Write_EXO_Result_Global(AppCtx%MyExo, 1, 1, AppCtx%Energy)
   Call PetscLogStagePop (AppCtx%LogInfo%IO_Stage, iErr); CHKERRQ(iErr)

   Call SimplePoissonFinalize(AppCtx)
#if defined PB_2D
End Program  SimplePoisson2D
#elif defined PB_3D
End Program SimplePoisson3D
#endif
