#if defined PB_2D
Module m_TSPoisson2D
#elif defined PB_3D
Module m_TSPoisson3D
#endif

#include "finclude/petscdef.h"

   Use m_MEF90

#if defined PB_2D
   Use m_Poisson2D
#elif defined PB_3D 
   Use m_Poisson3D
#endif     

   Implicit NONE   


   
   
Contains

  

#if defined PB_2D
End Module m_TSPoisson2D
#elif defined PB_3D
End Module m_TSPoisson3D
#endif
