Program  SimplePoisson

#include "finclude/petscdef.h"

   Use m_MEF90
#if defined PB_2D
   Use m_Poisson2D
#elif defined PB_3D 
   Use m_Poisson3D
#endif

   Implicit NONE   


   Type(Poisson_AppCtx_Type)                    :: AppCtx
   PetscErrorCode                               :: iErr
   PetscInt                                     :: exo_step,exo_field
   Character(len=MEF90_MXSTRLEN)                :: IOBuffer
   Type(Vec)                                    :: GradULocalVec
   PetscInt                                     :: numDim
   PetscReal                                    :: Umin,Umax
   PetscInt                                     :: imin,imax
   
   
   Call SimplePoissonInit(AppCtx)
   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Assembling the matrix\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   Call PoissonMatAssembly(AppCtx%SNESU,AppCtx%U%Vec,AppCtx%K,AppCtx%K,SAME_NONZERO_PATTERN,AppCtx,ierr)
   If (AppCtx%AppParam%verbose > 3) Then
      Write(IOBuffer, *) 'Matrix\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call MatView(AppCtx%K, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   End If

   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Assembling the RHS\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   Call RHSAssembly(AppCtx%SNESU,AppCtx%RHS%Vec,AppCtx%U%Vec,AppCtx)
   If (AppCtx%AppParam%verbose > 2) Then
      Write(IOBuffer, *) 'RHS\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call SectionRealView(AppCtx%RHS, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   End If
   
   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Calling Solve\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   Call Solve(AppCtx)
   
   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Computing Energies\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   Call ComputeEnergy(AppCtx)

   Write(IOBuffer, 100) AppCtx%ElasticEnergy, AppCtx%ExtForcesWork, AppCtx%TotalEnergy
100 Format('Elastic energy: ', ES12.5, ' Forces Work: ', ES12.5, ' Total: ', ES12.5, '\n')    
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

   Call ComputeGradU(AppCtx)

   If (AppCtx%AppParam%verbose > 0) Then
      Write(IOBuffer, *) 'Saving results\n'
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
   End If
   
   Call VecViewExodusVertex(AppCtx%mesh,AppCtx%U%LocalVec,AppCtx%EXO%comm,AppCtx%EXO%exoid,1,1,ierr)
   Call VecViewExodusVertex(AppCtx%mesh,AppCtx%F%LocalVec,AppCtx%EXO%comm,AppCtx%EXO%exoid,1,2,ierr)

   Call SectionRealCreateLocalVector(AppCtx%GradU,GradULocalVec,ierr);CHKERRQ(ierr)
   Call DMmeshGetDimension(AppCtx%mesh,numDim,ierr);CHKERRQ(ierr)
   Call VecSetBlockSize(GradULocalVec,numDim,ierr);CHKERRQ(ierr)
   Call VecViewExodusCell(AppCtx%mesh,GradULocalVec,AppCtx%EXO%comm,AppCtx%EXO%exoid,1,1,ierr)
   Call VecDestroy(GradULocalVec,ierr);CHKERRQ(ierr)
   
   Call Write_EXO_Result_Global(AppCtx%Exo, 1, 1, AppCtx%ElasticEnergy)
   Call Write_EXO_Result_Global(AppCtx%Exo, 2, 1, AppCtx%ExtForcesWork)
   Call Write_EXO_Result_Global(AppCtx%Exo, 3, 1, AppCtx%TotalEnergy)
  
   Call VecMin(AppCtx%U%Vec,imin,Umin,ierr);CHKERRQ(ierr)
   Write(*,*) 'U_Min: ', umin, 'at vertex ',imin
   Call VecMax(AppCtx%U%Vec,imax,Umax,ierr);CHKERRQ(ierr)
   Write(*,*) 'U_Max: ', umax, 'at vertex ',imax

   if (AppCtx%EXO%exoid > 0) Then
      Call EXCLOS(AppCtx%EXO%exoid, iErr)
   End If
   AppCtx%EXO%exoid = 0

   Call SimplePoissonFinalize(AppCtx)
End Program  SimplePoisson
