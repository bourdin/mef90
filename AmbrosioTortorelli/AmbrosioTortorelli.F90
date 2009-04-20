Program  Elast

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscmeshdef.h"

#if defined PB_2D
   Use m_AmbrosioTortorelli2D
#elif defined PB_3D
   Use m_AmbrosioTortorelli2D
#endif   
   Use m_MEF90
   Use m_RuptStruct
   Use petsc
   Use petscvec
   Use petscmat
   Use petscksp
   Use petscmesh
   
   Implicit NONE   


   Type(AppCtx_Type)                            :: AppCtx
   PetscInt                                     :: iErr
   Character(len=MEF90_MXSTRLEN)                :: IOBuffer

   Call AmbrosioTortorelliInit(AppCtx)
   
   If (AppCtx%AppParam%verbose) Then
      Call EXOView(AppCtx%EXO, AppCtx%AppParam%LogViewer)
      Call EXOView(AppCtx%MyEXO, AppCtx%AppParam%MyLogViewer)
      Call MeshTopologyView(AppCtx%MeshTopology, AppCtx%AppParam%MyLogViewer)
   End If   

   If (AppCtx%AppParam%verbose) Then
      Write(IOBuffer, *) 'Assembling the matrix of the U problem\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   Call MatAssembly_U(AppCtx)
   
!   Call MatView(AppCtx%KU, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   
   If (AppCtx%AppParam%verbose) Then
      Write(IOBuffer, *) 'Assembling the RHS of the U problem \n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   Call RHSAssembly_U(AppCtx)
   Call VecView(AppCtx%RHSU, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   
   If (AppCtx%AppParam%verbose) Then
      Write(IOBuffer, *) 'Calling KSPSolve for the U problem\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   AppCtx%TimeStep = AppCtx%TimeStep+1
   Call Solve_U(AppCtx)
   
   If (AppCtx%AppParam%verbose) Then
      Write(IOBuffer, *) 'Computing bulk energy, strains and stresses\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   Call ComputeEnergy(AppCtx)

   Write(IOBuffer, 100) AppCtx%BulkEnergy
100 Format('Total energy: ', ES12.5, '\n'c)    
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

   Call ComputeStrainStress(AppCtx)

   If (AppCtx%AppParam%verbose) Then
      Write(IOBuffer, *) 'Saving results\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   Call Save(AppCtx)

   Call AmbrosioTortorelliFinalize(AppCtx)
End Program  Elast
