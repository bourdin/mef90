Program  VarFracQSRecomputeEnergy
#include "finclude/petscdef.h"

   Use m_VarFilmQS
   Use m_MEF90
   Use m_Film_Struct
   
   Implicit NONE   

   Type(AppCtx_Type)                            :: AppCtx
   PetscInt                                     :: iErr, iBlk, i
   Character(len=MEF90_MXSTRLEN)                :: IOBuffer
   PetscInt                                     :: AltMinIter
   Character(len=MEF90_MXSTRLEN)                :: filename

   Call VarFracQSInit(AppCtx)
   
   TimeStep: Do i = 1, AppCtx%NumTimeSteps
      AppCtx%TimeStep = i
      Write(IOBuffer, 99) AppCtx%TimeStep
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
99    Format('\n=== Solving time step ', I4, '\n\n')

      !!! Init the fields:
      Call Init_TS_Loads(AppCtx)      
      !!! Update U at fixed nodes
      Call FieldInsertVertexBoundaryValues(AppCtx%U, AppCtx%UBC, AppCtx%BCUFlag, AppCtx%MeshTopology)
      
      !------------------------------------------------------------------- 
      ! Problem for U
      !-------------------------------------------------------------------
      Call Init_TS_Loads(AppCtx)      
      Call FieldInsertVertexBoundaryValues(AppCtx%U, AppCtx%UBC, AppCtx%BCUFlag, AppCtx%MeshTopology)

      Call Step_U(AppCtx)
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Fracture)%Offset, AppCtx%TimeStep, AppCtx%V)
      Call ComputeEnergies(AppCtx)


      Write(IOBuffer, 104) AppCtx%Load(AppCtx%TimeStep)
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, 100) AppCtx%ElasticEnergy(AppCtx%TimeStep)
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, 101) AppCtx%ExtForcesWork(AppCtx%TimeStep)
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, 102) AppCtx%FractureEnergy(AppCtx%TimeStep)
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Write(IOBuffer, 103) AppCtx%TotalEnergy(AppCtx%TimeStep)
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Call Save_Ener(AppCtx)
   End Do TimeStep

100   Format('Elastic energy:             ', ES12.5, '\n')    
101   Format('External Forces Work:       ', ES12.5, '\n')    
102   Format('Fracture energy:       ', ES12.5, '\n')    
103   Format('Total energy:         ', ES12.5, '\n')    
104   Format('Load:                 ', ES12.5, '\n')    

   Call VarFracQSFinalize(AppCtx)
End Program  VarFracQSRecomputeEnergy
