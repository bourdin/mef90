Program  VarFracFilm

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscmeshdef.h"

   Use m_VarFracFilm
   Use m_MEF90
   Use m_VarFracFilm_Struct
   Use petsc
   Use petscvec
   Use petscmat
   Use petscksp
   Use petscmesh
   
   Implicit NONE   


   Type(AppCtx_Type)                            :: AppCtx
   PetscInt                                     :: iErr
   Character(len=MEF90_MXSTRLEN)                :: IOBuffer
   PetscInt                                     :: iter

   Call VarFracFilmInit(AppCtx)
   
   If (AppCtx%AppParam%verbose) Then
      Call EXOView(AppCtx%EXO, AppCtx%AppParam%LogViewer)
      Call EXOView(AppCtx%MyEXO, AppCtx%AppParam%MyLogViewer) 
      Call MeshTopologyView(AppCtx%MeshTopology, AppCtx%AppParam%MyLogViewer) 
   End If   
   
   AppCtx%TimeStep = 1
   Do iter=1, AppCtx%VarFracFilmSchemeParam%AltMinMaxIter
      Write(IOBuffer, "('Iteration ', I4,A)") iter, '\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

      !------------------------------------------------------------------- 
      ! Problem for U
      !-------------------------------------------------------------------
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) 'Assembling the Matrix and RHS for  the U-subproblem \n'c 
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      Call RHSU_Assembly(AppCtx)
   !   Call VecView(AppCtx%RHSU, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
      Call MatU_Assembly(AppCtx)
   !   Call MatView(AppCtx%KU, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
      
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) 'Calling KSPSolve for the U-subproblem\n'c 
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
      End If
      Call Solve_U(AppCtx)
      
      !------------------------------------------------------------------- 
      ! Problem for V
      !-------------------------------------------------------------------
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) 'Assembling the Matrix and RHS for the V-subproblem \n'c 
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
      End If
      Call RHSV_Assembly(AppCtx)
   !   Call VecView(AppCtx%RHSV, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
      Call MatV_Assembly(AppCtx)
   !   Call MatView(AppCtx%KV, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
      
      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) 'Calling KSPSolve for the V-subproblem\n'c 
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
      End If
      Call Solve_V(AppCtx)
      
      !------------------------------------------------------------------- 
      ! Problem for phi
      !-------------------------------------------------------------------

      If (AppCtx%AppParam%verbose) Then
         Write(IOBuffer, *) 'Calling the Solver for the PHI-subproblem\n'c 
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
      End If
      Call Solve_PHI(AppCtx)
   
      

      !------------------------------------------------------------------- 
      ! Check the exit condition: tolerance on the error in V 
      !------------------------------------------------------------------- 
      If ((AppCtx%ErrV.LT.AppCtx%VarFracFilmSchemeParam%AltMinTol).AND.(AppCtx%ErrPHI.LT.AppCtx%VarFracFilmSchemeParam%AltMinTol)) then 
         EXIT 
      End If
      If (Mod(iter, AppCtx%VarFracFilmSchemeParam%AltMinSaveInt) == 0) Then
         If (AppCtx%AppParam%verbose) Then
            Write(IOBuffer, *) 'Saving U and V\n'c
            Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
         End If   
         Call Save_U(AppCtx)
         Call Save_V(AppCtx)
         Call Save_PHI(AppCtx)
         Call ComputeEnergy(AppCtx)
      
         Write(IOBuffer, 100) AppCtx%ElasticBulkEnergy
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Write(IOBuffer, 101) AppCtx%ElasticInterEnergy
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Write(IOBuffer, 102) AppCtx%SurfaceEnergyT
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Write(IOBuffer, 102) AppCtx%SurfaceEnergyD
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Write(IOBuffer, 103) AppCtx%TotalEnergy
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

!         Call ComputeStrainStress(AppCtx)
!         Call Save_StrainStress(AppCtx)
      End If
   End Do
   
   !------------------------------------------------------------------- ! Save the results
   Call Save_U(AppCtx)
   Call Save_V(AppCtx)
   Call Save_PHI(AppCtx)
   !-------------------------------------------------------------------
   
   If (AppCtx%AppParam%verbose) Then
      Write(IOBuffer, *) 'Computing bulk energy, strains and stresses\n'c 
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr) 
   End If
   
   Call ComputeEnergy(AppCtx)
   Call Save_Ener(AppCtx)

    Write(IOBuffer, 100) AppCtx%ElasticBulkEnergy
    Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
    Write(IOBuffer, 101) AppCtx%ElasticInterEnergy
    Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
    Write(IOBuffer, 102) AppCtx%SurfaceEnergyT
    Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
    Write(IOBuffer, 102) AppCtx%SurfaceEnergyD
    Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
    Write(IOBuffer, 103) AppCtx%TotalEnergy
    Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

100 Format('Elastic energy:       ', ES12.5, '\n'c)    
101 Format('External Forces Work: ', ES12.5, '\n'c)    
102 Format('Surface energy:       ', ES12.5, '\n'c)    
103 Format('Total energy:         ', ES12.5, '\n'c)    

!   Call ComputeStrainStress(AppCtx)

   Call VarFracFilmFinalize(AppCtx)
End Program  VarFracFilm
