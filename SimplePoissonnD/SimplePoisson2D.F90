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
   Character(len=MXSTLN)                        :: IOBuffer
   KSPConvergedReason                           :: reason

   Call SimplePoissonInit(AppCtx)

   Write(IOBuffer, *) 'Assembling the matrix\n'c
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

   Call MatAssembly(AppCtx)
   
   Write(IOBuffer, *) 'Assembling the RHS\n'c
   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

   Call SectionRealCreateLocalVector(AppCtx%F, F, iErr); CHKERRQ(iErr)
   Call VecSet(F, 1.0_Kr, iErr); CHKERRQ(iErr);
   Call VecDestroy(F, iErr); CHKERRQ(iErr)
   
   Call RHSAssembly(AppCtx)
   
   Call KSPSolve(AppCtx%KSP, AppCtx%RHS, AppCtx%RHS, iErr); CHKERRQ(iErr)
   Call VecView(AppCtx%RHS, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)   
   
!   Call Write_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, 1, 1, AppCtx%RHS) 
   !!! Why is this crashing?
!   Call KSPGetConvergedReason(AppCtx%KSP, reason, iErr); CHKERRQ(iErr)
!   Write(IOBuffer, *) 'KSPGetConvergedReason returned ', reason, '\n'c
!   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

   Call SimplePoissonFinalize(AppCtx)
#if defined PB_2D
End Program  SimplePoisson2D
#elif defined PB_3D
End Program SimplePoisson3D
#endif
