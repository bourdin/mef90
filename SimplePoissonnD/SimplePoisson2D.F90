#if defined PB_2D
Program  SimplePoisson2D
#elif defined PB_3D
Program SimplePoisson3D
#endif

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
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
   Use petscmesh

   Implicit NONE   


   Type(AppCtx_Type)                            :: AppCtx
   PetscInt                                     :: iErr
   
   Call SimplePoissonInit(AppCtx)
   
   
   Call SimplePoissonFinalize(AppCtx)
#if defined PB_2D
End Program  SimplePoisson2D
#elif defined PB_3D
End Program SimplePoisson3D
#endif
