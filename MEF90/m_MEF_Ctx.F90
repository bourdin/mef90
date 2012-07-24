Module m_MEF_Ctx
#include "finclude/petscdef.h"
   Implicit none
   Private  
   Public :: MEF90Ctx_Type
   
   Type MEF90Ctx_Type
      PetscBag                                        :: GlobalPropertiesBag
      PetscBag,Dimension(:),Pointer                   :: CellSetPropertiesBag
      PetscBag,Dimension(:),Pointer                   :: VertexSetPropertiesBag
      PetscBag,Dimension(:),Pointer                   :: MaterialPropertiesBag
   End Type MEF90Ctx_Type
End Module m_MEF_Ctx

